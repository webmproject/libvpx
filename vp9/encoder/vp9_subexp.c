/*
 *  Copyright (c) 2013 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include "vp9/common/vp9_common.h"
#include "vp9/common/vp9_entropy.h"

#include "vp9/encoder/vp9_boolhuff.h"
#include "vp9/encoder/vp9_treewriter.h"

#define vp9_cost_upd  ((int)(vp9_cost_one(upd) - vp9_cost_zero(upd)) >> 8)
#define vp9_cost_upd256  ((int)(vp9_cost_one(upd) - vp9_cost_zero(upd)))

static int update_bits[255];

static int count_uniform(int v, int n) {
  int l = get_unsigned_bits(n);
  int m;
  if (l == 0) return 0;
  m = (1 << l) - n;
  if (v < m)
    return l - 1;
  else
    return l;
}

static int split_index(int i, int n, int modulus) {
  int max1 = (n - 1 - modulus / 2) / modulus + 1;
  if (i % modulus == modulus / 2)
    i = i / modulus;
  else
    i = max1 + i - (i + modulus - modulus / 2) / modulus;
  return i;
}

static int recenter_nonneg(int v, int m) {
  if (v > (m << 1))
    return v;
  else if (v >= m)
    return ((v - m) << 1);
  else
    return ((m - v) << 1) - 1;
}

static int remap_prob(int v, int m) {
  const int n = 255;
  const int modulus = MODULUS_PARAM;
  int i;
  v--;
  m--;
  if ((m << 1) <= n)
    i = recenter_nonneg(v, m) - 1;
  else
    i = recenter_nonneg(n - 1 - v, n - 1 - m) - 1;

  i = split_index(i, n - 1, modulus);
  return i;
}

static int count_term_subexp(int word, int k, int num_syms) {
  int count = 0;
  int i = 0;
  int mk = 0;
  while (1) {
    int b = (i ? k + i - 1 : k);
    int a = (1 << b);
    if (num_syms <= mk + 3 * a) {
      count += count_uniform(word - mk, num_syms - mk);
      break;
    } else {
      int t = (word >= mk + a);
      count++;
      if (t) {
        i = i + 1;
        mk += a;
      } else {
        count += b;
        break;
      }
    }
  }
  return count;
}

static int prob_diff_update_cost(vp9_prob newp, vp9_prob oldp) {
  int delp = remap_prob(newp, oldp);
  return update_bits[delp] * 256;
}

static void encode_uniform(vp9_writer *w, int v, int n) {
  int l = get_unsigned_bits(n);
  int m;
  if (l == 0)
    return;
  m = (1 << l) - n;
  if (v < m) {
    vp9_write_literal(w, v, l - 1);
  } else {
    vp9_write_literal(w, m + ((v - m) >> 1), l - 1);
    vp9_write_literal(w, (v - m) & 1, 1);
  }
}

static void encode_term_subexp(vp9_writer *w, int word, int k, int num_syms) {
  int i = 0;
  int mk = 0;
  while (1) {
    int b = (i ? k + i - 1 : k);
    int a = (1 << b);
    if (num_syms <= mk + 3 * a) {
      encode_uniform(w, word - mk, num_syms - mk);
      break;
    } else {
      int t = (word >= mk + a);
      vp9_write_literal(w, t, 1);
      if (t) {
        i = i + 1;
        mk += a;
      } else {
        vp9_write_literal(w, word - mk, b);
        break;
      }
    }
  }
}

void vp9_write_prob_diff_update(vp9_writer *w, vp9_prob newp, vp9_prob oldp) {
  const int delp = remap_prob(newp, oldp);
  encode_term_subexp(w, delp, SUBEXP_PARAM, 255);
}

void vp9_compute_update_table() {
  int i;
  for (i = 0; i < 254; i++)
    update_bits[i] = count_term_subexp(i, SUBEXP_PARAM, 255);
}

int vp9_prob_diff_update_savings_search(const unsigned int *ct,
                                        vp9_prob oldp, vp9_prob *bestp,
                                        vp9_prob upd) {
  const int old_b = cost_branch256(ct, oldp);
  int bestsavings = 0;
  vp9_prob newp, bestnewp = oldp;
  const int step = *bestp > oldp ? -1 : 1;

  for (newp = *bestp; newp != oldp; newp += step) {
    const int new_b = cost_branch256(ct, newp);
    const int update_b = prob_diff_update_cost(newp, oldp) + vp9_cost_upd256;
    const int savings = old_b - new_b - update_b;
    if (savings > bestsavings) {
      bestsavings = savings;
      bestnewp = newp;
    }
  }
  *bestp = bestnewp;
  return bestsavings;
}

int vp9_prob_diff_update_savings_search_model(const unsigned int *ct,
                                              const vp9_prob *oldp,
                                              vp9_prob *bestp,
                                              vp9_prob upd,
                                              int b, int r) {
  int i, old_b, new_b, update_b, savings, bestsavings, step;
  int newp;
  vp9_prob bestnewp, newplist[ENTROPY_NODES], oldplist[ENTROPY_NODES];
  vp9_model_to_full_probs(oldp, oldplist);
  vpx_memcpy(newplist, oldp, sizeof(vp9_prob) * UNCONSTRAINED_NODES);
  for (i = UNCONSTRAINED_NODES, old_b = 0; i < ENTROPY_NODES; ++i)
    old_b += cost_branch256(ct + 2 * i, oldplist[i]);
  old_b += cost_branch256(ct + 2 * PIVOT_NODE, oldplist[PIVOT_NODE]);

  bestsavings = 0;
  bestnewp = oldp[PIVOT_NODE];

  step = (*bestp > oldp[PIVOT_NODE] ? -1 : 1);

  for (newp = *bestp; newp != oldp[PIVOT_NODE]; newp += step) {
    if (newp < 1 || newp > 255)
      continue;
    newplist[PIVOT_NODE] = newp;
    vp9_model_to_full_probs(newplist, newplist);
    for (i = UNCONSTRAINED_NODES, new_b = 0; i < ENTROPY_NODES; ++i)
      new_b += cost_branch256(ct + 2 * i, newplist[i]);
    new_b += cost_branch256(ct + 2 * PIVOT_NODE, newplist[PIVOT_NODE]);
    update_b = prob_diff_update_cost(newp, oldp[PIVOT_NODE]) +
        vp9_cost_upd256;
    savings = old_b - new_b - update_b;
    if (savings > bestsavings) {
      bestsavings = savings;
      bestnewp = newp;
    }
  }
  *bestp = bestnewp;
  return bestsavings;
}

void vp9_cond_prob_diff_update(vp9_writer *w, vp9_prob *oldp,
                               vp9_prob upd, unsigned int *ct) {
  vp9_prob newp = get_binary_prob(ct[0], ct[1]);
  const int savings = vp9_prob_diff_update_savings_search(ct, *oldp, &newp,
                                                          upd);
  assert(newp >= 1);
  if (savings > 0) {
    vp9_write(w, 1, upd);
    vp9_write_prob_diff_update(w, newp, *oldp);
    *oldp = newp;
  } else {
    vp9_write(w, 0, upd);
  }
}
