/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include <assert.h>
#include <stdio.h>
#include <limits.h>

#include "vpx/vpx_encoder.h"
#include "vpx_mem/vpx_mem.h"

#include "vp9/common/vp9_entropymode.h"
#include "vp9/common/vp9_entropymv.h"
#include "vp9/common/vp9_findnearmv.h"
#include "vp9/common/vp9_tile_common.h"
#include "vp9/common/vp9_seg_common.h"
#include "vp9/common/vp9_pred_common.h"
#include "vp9/common/vp9_entropy.h"
#include "vp9/common/vp9_entropymv.h"
#include "vp9/common/vp9_mvref_common.h"
#include "vp9/common/vp9_treecoder.h"
#include "vp9/common/vp9_systemdependent.h"
#include "vp9/common/vp9_pragmas.h"

#include "vp9/encoder/vp9_mcomp.h"
#include "vp9/encoder/vp9_encodemv.h"
#include "vp9/encoder/vp9_bitstream.h"
#include "vp9/encoder/vp9_segmentation.h"
#include "vp9/encoder/vp9_write_bit_buffer.h"


#if defined(SECTIONBITS_OUTPUT)
unsigned __int64 Sectionbits[500];
#endif

#ifdef ENTROPY_STATS
int intra_mode_stats[VP9_INTRA_MODES]
                    [VP9_INTRA_MODES]
                    [VP9_INTRA_MODES];
vp9_coeff_stats tree_update_hist[TX_SIZE_MAX_SB][BLOCK_TYPES];

extern unsigned int active_section;
#endif

#define vp9_cost_upd  ((int)(vp9_cost_one(upd) - vp9_cost_zero(upd)) >> 8)
#define vp9_cost_upd256  ((int)(vp9_cost_one(upd) - vp9_cost_zero(upd)))

#ifdef MODE_STATS
int64_t tx_count_32x32p_stats[TX_SIZE_CONTEXTS][TX_SIZE_MAX_SB];
int64_t tx_count_16x16p_stats[TX_SIZE_CONTEXTS][TX_SIZE_MAX_SB - 1];
int64_t tx_count_8x8p_stats[TX_SIZE_CONTEXTS][TX_SIZE_MAX_SB - 2];
int64_t switchable_interp_stats[VP9_SWITCHABLE_FILTERS+1]
                               [VP9_SWITCHABLE_FILTERS];

void init_tx_count_stats() {
  vp9_zero(tx_count_32x32p_stats);
  vp9_zero(tx_count_16x16p_stats);
  vp9_zero(tx_count_8x8p_stats);
}

void init_switchable_interp_stats() {
  vp9_zero(switchable_interp_stats);
}

static void update_tx_count_stats(VP9_COMMON *cm) {
  int i, j;
  for (i = 0; i < TX_SIZE_CONTEXTS; i++) {
    for (j = 0; j < TX_SIZE_MAX_SB; j++) {
      tx_count_32x32p_stats[i][j] += cm->fc.tx_count_32x32p[i][j];
    }
  }
  for (i = 0; i < TX_SIZE_CONTEXTS; i++) {
    for (j = 0; j < TX_SIZE_MAX_SB - 1; j++) {
      tx_count_16x16p_stats[i][j] += cm->fc.tx_count_16x16p[i][j];
    }
  }
  for (i = 0; i < TX_SIZE_CONTEXTS; i++) {
    for (j = 0; j < TX_SIZE_MAX_SB - 2; j++) {
      tx_count_8x8p_stats[i][j] += cm->fc.tx_count_8x8p[i][j];
    }
  }
}

static void update_switchable_interp_stats(VP9_COMMON *cm) {
  int i, j;
  for (i = 0; i < VP9_SWITCHABLE_FILTERS+1; ++i)
    for (j = 0; j < VP9_SWITCHABLE_FILTERS; ++j) {
      switchable_interp_stats[i][j] += cm->fc.switchable_interp_count[i][j];
    }
}

void write_tx_count_stats() {
  int i, j;
  FILE *fp = fopen("tx_count.bin", "wb");
  fwrite(tx_count_32x32p_stats, sizeof(tx_count_32x32p_stats), 1, fp);
  fwrite(tx_count_16x16p_stats, sizeof(tx_count_16x16p_stats), 1, fp);
  fwrite(tx_count_8x8p_stats, sizeof(tx_count_8x8p_stats), 1, fp);
  fclose(fp);

  printf(
      "vp9_default_tx_count_32x32p[TX_SIZE_CONTEXTS][TX_SIZE_MAX_SB] = {\n");
  for (i = 0; i < TX_SIZE_CONTEXTS; i++) {
    printf("  { ");
    for (j = 0; j < TX_SIZE_MAX_SB; j++) {
      printf("%"PRId64", ", tx_count_32x32p_stats[i][j]);
    }
    printf("},\n");
  }
  printf("};\n");
  printf(
      "vp9_default_tx_count_16x16p[TX_SIZE_CONTEXTS][TX_SIZE_MAX_SB-1] = {\n");
  for (i = 0; i < TX_SIZE_CONTEXTS; i++) {
    printf("  { ");
    for (j = 0; j < TX_SIZE_MAX_SB - 1; j++) {
      printf("%"PRId64", ", tx_count_16x16p_stats[i][j]);
    }
    printf("},\n");
  }
  printf("};\n");
  printf(
      "vp9_default_tx_count_8x8p[TX_SIZE_CONTEXTS][TX_SIZE_MAX_SB-2] = {\n");
  for (i = 0; i < TX_SIZE_CONTEXTS; i++) {
    printf("  { ");
    for (j = 0; j < TX_SIZE_MAX_SB - 2; j++) {
      printf("%"PRId64", ", tx_count_8x8p_stats[i][j]);
    }
    printf("},\n");
  }
  printf("};\n");
}

void write_switchable_interp_stats() {
  int i, j;
  FILE *fp = fopen("switchable_interp.bin", "wb");
  fwrite(switchable_interp_stats, sizeof(switchable_interp_stats), 1, fp);
  fclose(fp);

  printf(
      "vp9_default_switchable_filter_count[VP9_SWITCHABLE_FILTERS+1]"
      "[VP9_SWITCHABLE_FILTERS] = {\n");
  for (i = 0; i < VP9_SWITCHABLE_FILTERS+1; i++) {
    printf("  { ");
    for (j = 0; j < VP9_SWITCHABLE_FILTERS; j++) {
      printf("%"PRId64", ", switchable_interp_stats[i][j]);
    }
    printf("},\n");
  }
  printf("};\n");
}
#endif

static int update_bits[255];

static INLINE void write_be32(uint8_t *p, int value) {
  p[0] = value >> 24;
  p[1] = value >> 16;
  p[2] = value >> 8;
  p[3] = value;
}



int recenter_nonneg(int v, int m) {
  if (v > (m << 1))
    return v;
  else if (v >= m)
    return ((v - m) << 1);
  else
    return ((m - v) << 1) - 1;
}

static int get_unsigned_bits(unsigned num_values) {
  int cat = 0;
  if ((num_values--) <= 1) return 0;
  while (num_values > 0) {
    cat++;
    num_values >>= 1;
  }
  return cat;
}

void vp9_encode_unsigned_max(struct vp9_write_bit_buffer *wb,
                             int data, int max) {
  vp9_wb_write_literal(wb, data, get_unsigned_bits(max));
}

void encode_uniform(vp9_writer *w, int v, int n) {
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

int count_uniform(int v, int n) {
  int l = get_unsigned_bits(n);
  int m;
  if (l == 0) return 0;
  m = (1 << l) - n;
  if (v < m)
    return l - 1;
  else
    return l;
}

void encode_term_subexp(vp9_writer *w, int word, int k, int num_syms) {
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

int count_term_subexp(int word, int k, int num_syms) {
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

static void compute_update_table() {
  int i;
  for (i = 0; i < 254; i++)
    update_bits[i] = count_term_subexp(i, SUBEXP_PARAM, 255);
}

static int split_index(int i, int n, int modulus) {
  int max1 = (n - 1 - modulus / 2) / modulus + 1;
  if (i % modulus == modulus / 2) i = i / modulus;
  else i = max1 + i - (i + modulus - modulus / 2) / modulus;
  return i;
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

static void write_prob_diff_update(vp9_writer *w,
                                   vp9_prob newp, vp9_prob oldp) {
  int delp = remap_prob(newp, oldp);
  encode_term_subexp(w, delp, SUBEXP_PARAM, 255);
}

static int prob_diff_update_cost(vp9_prob newp, vp9_prob oldp) {
  int delp = remap_prob(newp, oldp);
  return update_bits[delp] * 256;
}

static int prob_update_savings(const unsigned int *ct,
                               const vp9_prob oldp, const vp9_prob newp,
                               const vp9_prob upd) {
  const int old_b = cost_branch256(ct, oldp);
  const int new_b = cost_branch256(ct, newp);
  const int update_b = 2048 + vp9_cost_upd256;
  return old_b - new_b - update_b;
}

static int prob_diff_update_savings_search(const unsigned int *ct,
                                           const vp9_prob oldp, vp9_prob *bestp,
                                           const vp9_prob upd) {
  const int old_b = cost_branch256(ct, oldp);
  int new_b, update_b, savings, bestsavings, step;
  vp9_prob newp, bestnewp;

  bestsavings = 0;
  bestnewp = oldp;

  step = (*bestp > oldp ? -1 : 1);
  for (newp = *bestp; newp != oldp; newp += step) {
    new_b = cost_branch256(ct, newp);
    update_b = prob_diff_update_cost(newp, oldp) + vp9_cost_upd256;
    savings = old_b - new_b - update_b;
    if (savings > bestsavings) {
      bestsavings = savings;
      bestnewp = newp;
    }
  }
  *bestp = bestnewp;
  return bestsavings;
}

static int prob_diff_update_savings_search_model(const unsigned int *ct,
                                                 const vp9_prob *oldp,
                                                 vp9_prob *bestp,
                                                 const vp9_prob upd,
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
  newp = *bestp;
  for (; newp != oldp[PIVOT_NODE]; newp += step) {
    if (newp < 1 || newp > 255) continue;
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

static void vp9_cond_prob_update(vp9_writer *bc, vp9_prob *oldp, vp9_prob upd,
                                 unsigned int *ct) {
  vp9_prob newp;
  int savings;
  newp = get_binary_prob(ct[0], ct[1]);
  assert(newp >= 1);
  savings = prob_update_savings(ct, *oldp, newp, upd);
  if (savings > 0) {
    vp9_write(bc, 1, upd);
    vp9_write_prob(bc, newp);
    *oldp = newp;
  } else {
    vp9_write(bc, 0, upd);
  }
}

static void vp9_cond_prob_diff_update(vp9_writer *bc, vp9_prob *oldp,
                                      vp9_prob upd,
                                      unsigned int *ct) {
  vp9_prob newp;
  int savings;
  newp = get_binary_prob(ct[0], ct[1]);
  assert(newp >= 1);
  savings = prob_diff_update_savings_search(ct, *oldp, &newp, upd);
  if (savings > 0) {
    vp9_write(bc, 1, upd);
    write_prob_diff_update(bc, newp, *oldp);
    *oldp = newp;
  } else {
    vp9_write(bc, 0, upd);
  }
}

static void update_mode(
  vp9_writer *w,
  int n,
  const struct vp9_token tok[/* n */],
  vp9_tree tree,
  vp9_prob Pnew[/* n-1 */],
  vp9_prob Pcur[/* n-1 */],
  unsigned int bct[/* n-1 */] [2],
  const unsigned int num_events[/* n */]
) {
  int i = 0;

  vp9_tree_probs_from_distribution(tree, Pnew, bct, num_events, 0);
  n--;

  for (i = 0; i < n; ++i) {
    vp9_cond_prob_diff_update(w, &Pcur[i], VP9_MODE_UPDATE_PROB, bct[i]);
  }
}

static void update_mbintra_mode_probs(VP9_COMP* const cpi,
                                      vp9_writer* const bc) {
  VP9_COMMON *const cm = &cpi->common;
  int j;
  vp9_prob pnew[VP9_INTRA_MODES - 1];
  unsigned int bct[VP9_INTRA_MODES - 1][2];

  for (j = 0; j < BLOCK_SIZE_GROUPS; j++)
    update_mode(bc, VP9_INTRA_MODES, vp9_intra_mode_encodings,
                vp9_intra_mode_tree, pnew,
                cm->fc.y_mode_prob[j], bct,
                (unsigned int *)cpi->y_mode_count[j]);
}

void vp9_update_skip_probs(VP9_COMP *cpi, vp9_writer *bc) {
  VP9_COMMON *const pc = &cpi->common;
  int k;

  for (k = 0; k < MBSKIP_CONTEXTS; ++k) {
    vp9_cond_prob_diff_update(bc, &pc->fc.mbskip_probs[k],
                              VP9_MODE_UPDATE_PROB, pc->fc.mbskip_count[k]);
  }
}

static void write_intra_mode(vp9_writer *bc, int m, const vp9_prob *p) {
  write_token(bc, vp9_intra_mode_tree, p, vp9_intra_mode_encodings + m);
}

static void update_switchable_interp_probs(VP9_COMP *const cpi,
                                           vp9_writer* const bc) {
  VP9_COMMON *const pc = &cpi->common;
  unsigned int branch_ct[VP9_SWITCHABLE_FILTERS + 1]
                        [VP9_SWITCHABLE_FILTERS - 1][2];
  vp9_prob new_prob[VP9_SWITCHABLE_FILTERS + 1][VP9_SWITCHABLE_FILTERS - 1];
  int i, j;
  for (j = 0; j <= VP9_SWITCHABLE_FILTERS; ++j) {
    vp9_tree_probs_from_distribution(
        vp9_switchable_interp_tree,
        new_prob[j], branch_ct[j],
        pc->fc.switchable_interp_count[j], 0);
  }
  for (j = 0; j <= VP9_SWITCHABLE_FILTERS; ++j) {
    for (i = 0; i < VP9_SWITCHABLE_FILTERS - 1; ++i) {
      vp9_cond_prob_diff_update(bc, &pc->fc.switchable_interp_prob[j][i],
                                VP9_MODE_UPDATE_PROB, branch_ct[j][i]);
    }
  }
#ifdef MODE_STATS
  if (!cpi->dummy_packing)
    update_switchable_interp_stats(pc);
#endif
}

static void update_inter_mode_probs(VP9_COMMON *pc, vp9_writer* const bc) {
  int i, j;

  for (i = 0; i < INTER_MODE_CONTEXTS; i++) {
    for (j = 0; j < VP9_INTER_MODES - 1; j++) {
      vp9_cond_prob_diff_update(bc, &pc->fc.inter_mode_probs[i][j],
                                VP9_MODE_UPDATE_PROB,
                                pc->fc.inter_mode_counts[i][j]);
    }
  }
}

static void pack_mb_tokens(vp9_writer* const bc,
                           TOKENEXTRA **tp,
                           const TOKENEXTRA *const stop) {
  TOKENEXTRA *p = *tp;

  while (p < stop) {
    const int t = p->token;
    const struct vp9_token *const a = vp9_coef_encodings + t;
    const vp9_extra_bit *const b = vp9_extra_bits + t;
    int i = 0;
    const vp9_prob *pp;
    int v = a->value;
    int n = a->len;
    vp9_prob probs[ENTROPY_NODES];

    if (t == EOSB_TOKEN) {
      ++p;
      break;
    }
    if (t >= TWO_TOKEN) {
      vp9_model_to_full_probs(p->context_tree, probs);
      pp = probs;
    } else {
      pp = p->context_tree;
    }
    assert(pp != 0);

    /* skip one or two nodes */
#if !CONFIG_BALANCED_COEFTREE
    if (p->skip_eob_node) {
      n -= p->skip_eob_node;
      i = 2 * p->skip_eob_node;
    }
#endif

    do {
      const int bb = (v >> --n) & 1;
#if CONFIG_BALANCED_COEFTREE
      if (i == 2 && p->skip_eob_node) {
        i += 2;
        assert(bb == 1);
        continue;
      }
#endif
      vp9_write(bc, bb, pp[i >> 1]);
      i = vp9_coef_tree[i + bb];
    } while (n);

    if (b->base_val) {
      const int e = p->extra, l = b->len;

      if (l) {
        const unsigned char *pb = b->prob;
        int v = e >> 1;
        int n = l;              /* number of bits in v, assumed nonzero */
        int i = 0;

        do {
          const int bb = (v >> --n) & 1;
          vp9_write(bc, bb, pb[i >> 1]);
          i = b->tree[i + bb];
        } while (n);
      }

      vp9_write_bit(bc, e & 1);
    }
    ++p;
  }

  *tp = p;
}

static void write_sb_mv_ref(vp9_writer *bc, MB_PREDICTION_MODE m,
                            const vp9_prob *p) {
#if CONFIG_DEBUG
  assert(NEARESTMV <= m && m <= NEWMV);
#endif
  write_token(bc, vp9_sb_mv_ref_tree, p,
              vp9_sb_mv_ref_encoding_array - NEARESTMV + m);
}

// This function writes the current macro block's segnment id to the bitstream
// It should only be called if a segment map update is indicated.
static void write_mb_segid(vp9_writer *bc,
                           const MB_MODE_INFO *mi, const MACROBLOCKD *xd) {
  if (xd->segmentation_enabled && xd->update_mb_segmentation_map)
    treed_write(bc, vp9_segment_tree, xd->mb_segment_tree_probs,
                mi->segment_id, 3);
}

// This function encodes the reference frame
static void encode_ref_frame(VP9_COMP *cpi, vp9_writer *bc) {
  VP9_COMMON *const pc = &cpi->common;
  MACROBLOCK *const x = &cpi->mb;
  MACROBLOCKD *const xd = &x->e_mbd;
  MB_MODE_INFO *mi = &xd->mode_info_context->mbmi;
  const int segment_id = mi->segment_id;
  int seg_ref_active = vp9_segfeature_active(xd, segment_id,
                                             SEG_LVL_REF_FRAME);
  // If segment level coding of this signal is disabled...
  // or the segment allows multiple reference frame options
  if (!seg_ref_active) {
    // does the feature use compound prediction or not
    // (if not specified at the frame/segment level)
    if (pc->comp_pred_mode == HYBRID_PREDICTION) {
      vp9_write(bc, mi->ref_frame[1] > INTRA_FRAME,
                vp9_get_pred_prob(pc, xd, PRED_COMP_INTER_INTER));
    } else {
      assert((mi->ref_frame[1] <= INTRA_FRAME) ==
                 (pc->comp_pred_mode == SINGLE_PREDICTION_ONLY));
    }

    if (mi->ref_frame[1] > INTRA_FRAME) {
      vp9_write(bc, mi->ref_frame[0] == GOLDEN_FRAME,
                vp9_get_pred_prob(pc, xd, PRED_COMP_REF_P));
    } else {
      vp9_write(bc, mi->ref_frame[0] != LAST_FRAME,
                vp9_get_pred_prob(pc, xd, PRED_SINGLE_REF_P1));
      if (mi->ref_frame[0] != LAST_FRAME)
        vp9_write(bc, mi->ref_frame[0] != GOLDEN_FRAME,
                  vp9_get_pred_prob(pc, xd, PRED_SINGLE_REF_P2));
    }
  } else {
    assert(mi->ref_frame[1] <= INTRA_FRAME);
    assert(vp9_get_segdata(xd, segment_id, SEG_LVL_REF_FRAME) ==
           mi->ref_frame[0]);
  }

  // if using the prediction mdoel we have nothing further to do because
  // the reference frame is fully coded by the segment
}

static void pack_inter_mode_mvs(VP9_COMP *cpi, MODE_INFO *m,
                                vp9_writer *bc, int mi_row, int mi_col) {
  VP9_COMMON *const pc = &cpi->common;
  const nmv_context *nmvc = &pc->fc.nmvc;
  MACROBLOCK *const x = &cpi->mb;
  MACROBLOCKD *const xd = &x->e_mbd;
  MB_MODE_INFO *const mi = &m->mbmi;
  const MV_REFERENCE_FRAME rf = mi->ref_frame[0];
  const MB_PREDICTION_MODE mode = mi->mode;
  const int segment_id = mi->segment_id;
  int skip_coeff;

  xd->prev_mode_info_context = pc->prev_mi + (m - pc->mi);
  x->partition_info = x->pi + (m - pc->mi);

#ifdef ENTROPY_STATS
  active_section = 9;
#endif

  if (cpi->mb.e_mbd.update_mb_segmentation_map) {
    // Is temporal coding of the segment map enabled
    if (pc->temporal_update) {
      unsigned char prediction_flag = vp9_get_pred_flag(xd, PRED_SEG_ID);
      vp9_prob pred_prob = vp9_get_pred_prob(pc, xd, PRED_SEG_ID);

      // Code the segment id prediction flag for this mb
      vp9_write(bc, prediction_flag, pred_prob);

      // If the mb segment id wasn't predicted code explicitly
      if (!prediction_flag)
        write_mb_segid(bc, mi, &cpi->mb.e_mbd);
    } else {
      // Normal unpredicted coding
      write_mb_segid(bc, mi, &cpi->mb.e_mbd);
    }
  }

  if (vp9_segfeature_active(xd, segment_id, SEG_LVL_SKIP)) {
    skip_coeff = 1;
  } else {
    skip_coeff = m->mbmi.mb_skip_coeff;
    vp9_write(bc, skip_coeff,
              vp9_get_pred_prob(pc, xd, PRED_MBSKIP));
  }

  if (!vp9_segfeature_active(xd, segment_id, SEG_LVL_REF_FRAME))
    vp9_write(bc, rf != INTRA_FRAME,
              vp9_get_pred_prob(pc, xd, PRED_INTRA_INTER));

  if (mi->sb_type >= BLOCK_SIZE_SB8X8 && pc->txfm_mode == TX_MODE_SELECT &&
      !(rf != INTRA_FRAME &&
        (skip_coeff || vp9_segfeature_active(xd, segment_id, SEG_LVL_SKIP)))) {
    TX_SIZE sz = mi->txfm_size;
    const vp9_prob *tx_probs = vp9_get_pred_probs(pc, xd, PRED_TX_SIZE);
    vp9_write(bc, sz != TX_4X4, tx_probs[0]);
    if (mi->sb_type >= BLOCK_SIZE_MB16X16 && sz != TX_4X4) {
      vp9_write(bc, sz != TX_8X8, tx_probs[1]);
      if (mi->sb_type >= BLOCK_SIZE_SB32X32 && sz != TX_8X8)
        vp9_write(bc, sz != TX_16X16, tx_probs[2]);
    }
  }

  if (rf == INTRA_FRAME) {
#ifdef ENTROPY_STATS
    active_section = 6;
#endif

    if (m->mbmi.sb_type >= BLOCK_SIZE_SB8X8) {
      const BLOCK_SIZE_TYPE bsize = xd->mode_info_context->mbmi.sb_type;
      const int bwl = b_width_log2(bsize), bhl = b_height_log2(bsize);
      const int bsl = MIN(bwl, bhl);
      write_intra_mode(bc, mode, pc->fc.y_mode_prob[MIN(3, bsl)]);
    } else {
      int idx, idy;
      int bw = 1 << b_width_log2(mi->sb_type);
      int bh = 1 << b_height_log2(mi->sb_type);
      for (idy = 0; idy < 2; idy += bh)
        for (idx = 0; idx < 2; idx += bw) {
          MB_PREDICTION_MODE bm = m->bmi[idy * 2 + idx].as_mode.first;
          write_intra_mode(bc, bm, pc->fc.y_mode_prob[0]);
        }
    }
    write_intra_mode(bc, mi->uv_mode,
                     pc->fc.uv_mode_prob[mode]);
  } else {
    vp9_prob *mv_ref_p;

    encode_ref_frame(cpi, bc);

    mv_ref_p = cpi->common.fc.inter_mode_probs[mi->mb_mode_context[rf]];

#ifdef ENTROPY_STATS
    active_section = 3;
#endif

    // If segment skip is not enabled code the mode.
    if (!vp9_segfeature_active(xd, segment_id, SEG_LVL_SKIP)) {
      if (mi->sb_type >= BLOCK_SIZE_SB8X8) {
        write_sb_mv_ref(bc, mode, mv_ref_p);
        vp9_accum_mv_refs(&cpi->common, mode, mi->mb_mode_context[rf]);
      }
    }

    if (cpi->common.mcomp_filter_type == SWITCHABLE) {
      write_token(bc, vp9_switchable_interp_tree,
                  vp9_get_pred_probs(&cpi->common, xd,
                                     PRED_SWITCHABLE_INTERP),
                  vp9_switchable_interp_encodings +
                  vp9_switchable_interp_map[mi->interp_filter]);
    } else {
      assert(mi->interp_filter == cpi->common.mcomp_filter_type);
    }

    if (xd->mode_info_context->mbmi.sb_type < BLOCK_SIZE_SB8X8) {
      int j;
      MB_PREDICTION_MODE blockmode;
      int_mv blockmv;
      int bwl = b_width_log2(mi->sb_type), bw = 1 << bwl;
      int bhl = b_height_log2(mi->sb_type), bh = 1 << bhl;
      int idx, idy;
      for (idy = 0; idy < 2; idy += bh) {
        for (idx = 0; idx < 2; idx += bw) {
          j = idy * 2 + idx;
          blockmode = cpi->mb.partition_info->bmi[j].mode;
          blockmv = cpi->mb.partition_info->bmi[j].mv;
          write_sb_mv_ref(bc, blockmode, mv_ref_p);
          vp9_accum_mv_refs(&cpi->common, blockmode, mi->mb_mode_context[rf]);
          if (blockmode == NEWMV) {
#ifdef ENTROPY_STATS
            active_section = 11;
#endif
            vp9_encode_mv(bc, &blockmv.as_mv, &mi->best_mv.as_mv,
                          nmvc, xd->allow_high_precision_mv);

            if (mi->ref_frame[1] > INTRA_FRAME)
              vp9_encode_mv(bc,
                            &cpi->mb.partition_info->bmi[j].second_mv.as_mv,
                            &mi->best_second_mv.as_mv,
                            nmvc, xd->allow_high_precision_mv);
          }
        }
      }
    } else if (mode == NEWMV) {
#ifdef ENTROPY_STATS
      active_section = 5;
#endif
      vp9_encode_mv(bc,
                    &mi->mv[0].as_mv, &mi->best_mv.as_mv,
                    nmvc, xd->allow_high_precision_mv);

      if (mi->ref_frame[1] > INTRA_FRAME)
        vp9_encode_mv(bc,
                      &mi->mv[1].as_mv, &mi->best_second_mv.as_mv,
                      nmvc, xd->allow_high_precision_mv);
    }
  }
}

static void write_mb_modes_kf(const VP9_COMP *cpi,
                              MODE_INFO *m,
                              vp9_writer *bc, int mi_row, int mi_col) {
  const VP9_COMMON *const c = &cpi->common;
  const MACROBLOCKD *const xd = &cpi->mb.e_mbd;
  const int ym = m->mbmi.mode;
  const int mis = c->mode_info_stride;
  const int segment_id = m->mbmi.segment_id;
  int skip_coeff;

  if (xd->update_mb_segmentation_map)
    write_mb_segid(bc, &m->mbmi, xd);

  if (vp9_segfeature_active(xd, segment_id, SEG_LVL_SKIP)) {
    skip_coeff = 1;
  } else {
    skip_coeff = m->mbmi.mb_skip_coeff;
    vp9_write(bc, skip_coeff, vp9_get_pred_prob(c, xd, PRED_MBSKIP));
  }

  if (m->mbmi.sb_type >= BLOCK_SIZE_SB8X8 && c->txfm_mode == TX_MODE_SELECT) {
    TX_SIZE sz = m->mbmi.txfm_size;
    const vp9_prob *tx_probs = vp9_get_pred_probs(c, xd, PRED_TX_SIZE);
    vp9_write(bc, sz != TX_4X4, tx_probs[0]);
    if (m->mbmi.sb_type >= BLOCK_SIZE_MB16X16 && sz != TX_4X4) {
      vp9_write(bc, sz != TX_8X8, tx_probs[1]);
      if (m->mbmi.sb_type >= BLOCK_SIZE_SB32X32 && sz != TX_8X8)
        vp9_write(bc, sz != TX_16X16, tx_probs[2]);
    }
  }

  if (m->mbmi.sb_type >= BLOCK_SIZE_SB8X8) {
    const MB_PREDICTION_MODE A = above_block_mode(m, 0, mis);
    const MB_PREDICTION_MODE L = xd->left_available ?
                                 left_block_mode(m, 0) : DC_PRED;
    write_intra_mode(bc, ym, c->kf_y_mode_prob[A][L]);
  } else {
    int idx, idy;
    int bw = 1 << b_width_log2(m->mbmi.sb_type);
    int bh = 1 << b_height_log2(m->mbmi.sb_type);
    for (idy = 0; idy < 2; idy += bh) {
      for (idx = 0; idx < 2; idx += bw) {
        int i = idy * 2 + idx;
        const MB_PREDICTION_MODE A = above_block_mode(m, i, mis);
        const MB_PREDICTION_MODE L = (xd->left_available || idx) ?
                                     left_block_mode(m, i) : DC_PRED;
        const int bm = m->bmi[i].as_mode.first;
#ifdef ENTROPY_STATS
        ++intra_mode_stats[A][L][bm];
#endif
        write_intra_mode(bc, bm, c->kf_y_mode_prob[A][L]);
      }
    }
  }

  write_intra_mode(bc, m->mbmi.uv_mode, c->kf_uv_mode_prob[ym]);
}

static void write_modes_b(VP9_COMP *cpi, MODE_INFO *m, vp9_writer *bc,
                          TOKENEXTRA **tok, TOKENEXTRA *tok_end,
                          int mi_row, int mi_col) {
  VP9_COMMON *const cm = &cpi->common;
  MACROBLOCKD *const xd = &cpi->mb.e_mbd;

  if (m->mbmi.sb_type < BLOCK_SIZE_SB8X8)
    if (xd->ab_index > 0)
      return;
  xd->mode_info_context = m;
  set_mi_row_col(&cpi->common, xd, mi_row,
                 1 << mi_height_log2(m->mbmi.sb_type),
                 mi_col, 1 << mi_width_log2(m->mbmi.sb_type));
  if ((cm->frame_type == KEY_FRAME) || cm->intra_only) {
    write_mb_modes_kf(cpi, m, bc, mi_row, mi_col);
#ifdef ENTROPY_STATS
    active_section = 8;
#endif
  } else {
    pack_inter_mode_mvs(cpi, m, bc, mi_row, mi_col);
#ifdef ENTROPY_STATS
    active_section = 1;
#endif
  }

  assert(*tok < tok_end);
  pack_mb_tokens(bc, tok, tok_end);
}

static void write_modes_sb(VP9_COMP *cpi, MODE_INFO *m, vp9_writer *bc,
                           TOKENEXTRA **tok, TOKENEXTRA *tok_end,
                           int mi_row, int mi_col,
                           BLOCK_SIZE_TYPE bsize) {
  VP9_COMMON *const cm = &cpi->common;
  MACROBLOCKD *xd = &cpi->mb.e_mbd;
  const int mis = cm->mode_info_stride;
  int bwl, bhl;
  int bsl = b_width_log2(bsize);
  int bs = (1 << bsl) / 4;  // mode_info step for subsize
  int n;
  PARTITION_TYPE partition;
  BLOCK_SIZE_TYPE subsize;

  if (mi_row >= cm->mi_rows || mi_col >= cm->mi_cols)
    return;

  bwl = b_width_log2(m->mbmi.sb_type);
  bhl = b_height_log2(m->mbmi.sb_type);

  // parse the partition type
  if ((bwl == bsl) && (bhl == bsl))
    partition = PARTITION_NONE;
  else if ((bwl == bsl) && (bhl < bsl))
    partition = PARTITION_HORZ;
  else if ((bwl < bsl) && (bhl == bsl))
    partition = PARTITION_VERT;
  else if ((bwl < bsl) && (bhl < bsl))
    partition = PARTITION_SPLIT;
  else
    assert(0);

  if (bsize < BLOCK_SIZE_SB8X8)
    if (xd->ab_index > 0)
      return;

  if (bsize >= BLOCK_SIZE_SB8X8) {
    int pl;
    int idx = check_bsize_coverage(cm, xd, mi_row, mi_col, bsize);
    xd->left_seg_context = cm->left_seg_context + (mi_row & MI_MASK);
    xd->above_seg_context = cm->above_seg_context + mi_col;
    pl = partition_plane_context(xd, bsize);
    // encode the partition information
    if (idx == 0)
      write_token(bc, vp9_partition_tree,
                  cm->fc.partition_prob[cm->frame_type][pl],
                  vp9_partition_encodings + partition);
    else if (idx > 0)
      vp9_write(bc, partition == PARTITION_SPLIT,
                cm->fc.partition_prob[cm->frame_type][pl][idx]);
  }

  subsize = get_subsize(bsize, partition);
  *(get_sb_index(xd, subsize)) = 0;

  switch (partition) {
    case PARTITION_NONE:
      write_modes_b(cpi, m, bc, tok, tok_end, mi_row, mi_col);
      break;
    case PARTITION_HORZ:
      write_modes_b(cpi, m, bc, tok, tok_end, mi_row, mi_col);
      *(get_sb_index(xd, subsize)) = 1;
      if ((mi_row + bs) < cm->mi_rows)
        write_modes_b(cpi, m + bs * mis, bc, tok, tok_end, mi_row + bs, mi_col);
      break;
    case PARTITION_VERT:
      write_modes_b(cpi, m, bc, tok, tok_end, mi_row, mi_col);
      *(get_sb_index(xd, subsize)) = 1;
      if ((mi_col + bs) < cm->mi_cols)
        write_modes_b(cpi, m + bs, bc, tok, tok_end, mi_row, mi_col + bs);
      break;
    case PARTITION_SPLIT:
      for (n = 0; n < 4; n++) {
        int j = n >> 1, i = n & 0x01;
        *(get_sb_index(xd, subsize)) = n;
        write_modes_sb(cpi, m + j * bs * mis + i * bs, bc, tok, tok_end,
                       mi_row + j * bs, mi_col + i * bs, subsize);
      }
      break;
    default:
      assert(0);
  }

  // update partition context
  if (bsize >= BLOCK_SIZE_SB8X8 &&
      (bsize == BLOCK_SIZE_SB8X8 || partition != PARTITION_SPLIT)) {
    set_partition_seg_context(cm, xd, mi_row, mi_col);
    update_partition_context(xd, subsize, bsize);
  }
}

static void write_modes(VP9_COMP *cpi, vp9_writer* const bc,
                        TOKENEXTRA **tok, TOKENEXTRA *tok_end) {
  VP9_COMMON *const c = &cpi->common;
  const int mis = c->mode_info_stride;
  MODE_INFO *m, *m_ptr = c->mi;
  int mi_row, mi_col;

  m_ptr += c->cur_tile_mi_col_start + c->cur_tile_mi_row_start * mis;

  for (mi_row = c->cur_tile_mi_row_start;
       mi_row < c->cur_tile_mi_row_end;
       mi_row += 8, m_ptr += 8 * mis) {
    m = m_ptr;
    vpx_memset(c->left_seg_context, 0, sizeof(c->left_seg_context));
    for (mi_col = c->cur_tile_mi_col_start;
         mi_col < c->cur_tile_mi_col_end;
         mi_col += 64 / MI_SIZE, m += 64 / MI_SIZE)
      write_modes_sb(cpi, m, bc, tok, tok_end, mi_row, mi_col,
                     BLOCK_SIZE_SB64X64);
  }
}

/* This function is used for debugging probability trees. */
static void print_prob_tree(vp9_coeff_probs *coef_probs, int block_types) {
  /* print coef probability tree */
  int i, j, k, l, m;
  FILE *f = fopen("enc_tree_probs.txt", "a");
  fprintf(f, "{\n");
  for (i = 0; i < block_types; i++) {
    fprintf(f, "  {\n");
    for (j = 0; j < REF_TYPES; ++j) {
      fprintf(f, "  {\n");
      for (k = 0; k < COEF_BANDS; k++) {
        fprintf(f, "    {\n");
        for (l = 0; l < PREV_COEF_CONTEXTS; l++) {
          fprintf(f, "      {");
          for (m = 0; m < ENTROPY_NODES; m++) {
            fprintf(f, "%3u, ",
                    (unsigned int)(coef_probs[i][j][k][l][m]));
          }
        }
        fprintf(f, " }\n");
      }
      fprintf(f, "    }\n");
    }
    fprintf(f, "  }\n");
  }
  fprintf(f, "}\n");
  fclose(f);
}

static void build_tree_distribution(VP9_COMP *cpi, TX_SIZE txfm_size) {
  vp9_coeff_probs_model *coef_probs = cpi->frame_coef_probs[txfm_size];
  vp9_coeff_count *coef_counts = cpi->coef_counts[txfm_size];
  unsigned int (*eob_branch_ct)[REF_TYPES][COEF_BANDS][PREV_COEF_CONTEXTS] =
      cpi->common.fc.eob_branch_counts[txfm_size];
  vp9_coeff_stats *coef_branch_ct = cpi->frame_branch_ct[txfm_size];
  vp9_prob full_probs[ENTROPY_NODES];
  int i, j, k, l;

  for (i = 0; i < BLOCK_TYPES; ++i) {
    for (j = 0; j < REF_TYPES; ++j) {
      for (k = 0; k < COEF_BANDS; ++k) {
        for (l = 0; l < PREV_COEF_CONTEXTS; ++l) {
          if (l >= 3 && k == 0)
            continue;
          vp9_tree_probs_from_distribution(vp9_coef_tree,
                                           full_probs,
                                           coef_branch_ct[i][j][k][l],
                                           coef_counts[i][j][k][l], 0);
          vpx_memcpy(coef_probs[i][j][k][l], full_probs,
                     sizeof(vp9_prob) * UNCONSTRAINED_NODES);
#if CONFIG_BALANCED_COEFTREE
          coef_branch_ct[i][j][k][l][1][1] = eob_branch_ct[i][j][k][l] -
                                             coef_branch_ct[i][j][k][l][1][0];
          coef_probs[i][j][k][l][1] =
              get_binary_prob(coef_branch_ct[i][j][k][l][1][0],
                              coef_branch_ct[i][j][k][l][1][1]);
#else
          coef_branch_ct[i][j][k][l][0][1] = eob_branch_ct[i][j][k][l] -
                                             coef_branch_ct[i][j][k][l][0][0];
          coef_probs[i][j][k][l][0] =
              get_binary_prob(coef_branch_ct[i][j][k][l][0][0],
                              coef_branch_ct[i][j][k][l][0][1]);
#endif
#ifdef ENTROPY_STATS
          if (!cpi->dummy_packing) {
            int t;
            for (t = 0; t < MAX_ENTROPY_TOKENS; ++t)
              context_counters[txfm_size][i][j][k][l][t] +=
                  coef_counts[i][j][k][l][t];
            context_counters[txfm_size][i][j][k][l][MAX_ENTROPY_TOKENS] +=
                eob_branch_ct[i][j][k][l];
          }
#endif
        }
      }
    }
  }
}

static void build_coeff_contexts(VP9_COMP *cpi) {
  TX_SIZE t;
  for (t = TX_4X4; t <= TX_32X32; t++)
    build_tree_distribution(cpi, t);
}

static void update_coef_probs_common(vp9_writer* const bc, VP9_COMP *cpi,
                                     TX_SIZE tx_size) {
  vp9_coeff_probs_model *new_frame_coef_probs = cpi->frame_coef_probs[tx_size];
  vp9_coeff_probs_model *old_frame_coef_probs =
      cpi->common.fc.coef_probs[tx_size];
  vp9_coeff_stats *frame_branch_ct = cpi->frame_branch_ct[tx_size];
  int i, j, k, l, t;
  int update[2] = {0, 0};
  int savings;

  const int entropy_nodes_update = UNCONSTRAINED_NODES;

  const int tstart = 0;
  /* dry run to see if there is any udpate at all needed */
  savings = 0;
  for (i = 0; i < BLOCK_TYPES; ++i) {
    for (j = 0; j < REF_TYPES; ++j) {
      for (k = 0; k < COEF_BANDS; ++k) {
        // int prev_coef_savings[ENTROPY_NODES] = {0};
        for (l = 0; l < PREV_COEF_CONTEXTS; ++l) {
          for (t = tstart; t < entropy_nodes_update; ++t) {
            vp9_prob newp = new_frame_coef_probs[i][j][k][l][t];
            const vp9_prob oldp = old_frame_coef_probs[i][j][k][l][t];
            const vp9_prob upd = VP9_COEF_UPDATE_PROB;
            int s;
            int u = 0;

            if (l >= 3 && k == 0)
              continue;
            if (t == PIVOT_NODE)
              s = prob_diff_update_savings_search_model(
                  frame_branch_ct[i][j][k][l][0],
                  old_frame_coef_probs[i][j][k][l], &newp, upd, i, j);
            else
              s = prob_diff_update_savings_search(
                  frame_branch_ct[i][j][k][l][t], oldp, &newp, upd);
            if (s > 0 && newp != oldp)
              u = 1;
            if (u)
              savings += s - (int)(vp9_cost_zero(upd));
            else
              savings -= (int)(vp9_cost_zero(upd));
            update[u]++;
          }
        }
      }
    }
  }

  // printf("Update %d %d, savings %d\n", update[0], update[1], savings);
  /* Is coef updated at all */
  if (update[1] == 0 || savings < 0) {
    vp9_write_bit(bc, 0);
    return;
  }
  vp9_write_bit(bc, 1);
  for (i = 0; i < BLOCK_TYPES; ++i) {
    for (j = 0; j < REF_TYPES; ++j) {
      for (k = 0; k < COEF_BANDS; ++k) {
        // int prev_coef_savings[ENTROPY_NODES] = {0};
        for (l = 0; l < PREV_COEF_CONTEXTS; ++l) {
          // calc probs and branch cts for this frame only
          for (t = tstart; t < entropy_nodes_update; ++t) {
            vp9_prob newp = new_frame_coef_probs[i][j][k][l][t];
            vp9_prob *oldp = old_frame_coef_probs[i][j][k][l] + t;
            const vp9_prob upd = VP9_COEF_UPDATE_PROB;
            int s;
            int u = 0;
            if (l >= 3 && k == 0)
              continue;
            if (t == PIVOT_NODE)
              s = prob_diff_update_savings_search_model(
                  frame_branch_ct[i][j][k][l][0],
                  old_frame_coef_probs[i][j][k][l], &newp, upd, i, j);
            else
              s = prob_diff_update_savings_search(
                  frame_branch_ct[i][j][k][l][t],
                  *oldp, &newp, upd);
            if (s > 0 && newp != *oldp)
              u = 1;
            vp9_write(bc, u, upd);
#ifdef ENTROPY_STATS
            if (!cpi->dummy_packing)
              ++tree_update_hist[tx_size][i][j][k][l][t][u];
#endif
            if (u) {
              /* send/use new probability */
              write_prob_diff_update(bc, newp, *oldp);
              *oldp = newp;
            }
          }
        }
      }
    }
  }
}

static void update_coef_probs(VP9_COMP* const cpi, vp9_writer* const bc) {
  const TXFM_MODE txfm_mode = cpi->common.txfm_mode;

  vp9_clear_system_state();

  // Build the cofficient contexts based on counts collected in encode loop
  build_coeff_contexts(cpi);

  update_coef_probs_common(bc, cpi, TX_4X4);

  // do not do this if not even allowed
  if (txfm_mode > ONLY_4X4)
    update_coef_probs_common(bc, cpi, TX_8X8);

  if (txfm_mode > ALLOW_8X8)
    update_coef_probs_common(bc, cpi, TX_16X16);

  if (txfm_mode > ALLOW_16X16)
    update_coef_probs_common(bc, cpi, TX_32X32);
}

static void encode_loopfilter(VP9_COMMON *pc, MACROBLOCKD *xd,
                              struct vp9_write_bit_buffer *wb) {
  int i;

  // Encode the loop filter level and type
  vp9_wb_write_literal(wb, pc->filter_level, 6);
  vp9_wb_write_literal(wb, pc->sharpness_level, 3);

  // Write out loop filter deltas applied at the MB level based on mode or
  // ref frame (if they are enabled).
  vp9_wb_write_bit(wb, xd->mode_ref_lf_delta_enabled);

  if (xd->mode_ref_lf_delta_enabled) {
    // Do the deltas need to be updated
    vp9_wb_write_bit(wb, xd->mode_ref_lf_delta_update);
    if (xd->mode_ref_lf_delta_update) {
      // Send update
      for (i = 0; i < MAX_REF_LF_DELTAS; i++) {
        const int delta = xd->ref_lf_deltas[i];

        // Frame level data
        if (delta != xd->last_ref_lf_deltas[i]) {
          xd->last_ref_lf_deltas[i] = delta;
          vp9_wb_write_bit(wb, 1);

          assert(delta != 0);
          vp9_wb_write_literal(wb, abs(delta) & 0x3F, 6);
          vp9_wb_write_bit(wb, delta < 0);
        } else {
          vp9_wb_write_bit(wb, 0);
        }
      }

      // Send update
      for (i = 0; i < MAX_MODE_LF_DELTAS; i++) {
        const int delta = xd->mode_lf_deltas[i];
        if (delta != xd->last_mode_lf_deltas[i]) {
          xd->last_mode_lf_deltas[i] = delta;
          vp9_wb_write_bit(wb, 1);

          assert(delta != 0);
          vp9_wb_write_literal(wb, abs(delta) & 0x3F, 6);
          vp9_wb_write_bit(wb, delta < 0);
        } else {
          vp9_wb_write_bit(wb, 0);
        }
      }
    }
  }
}

static void write_delta_q(struct vp9_write_bit_buffer *wb, int delta_q) {
  if (delta_q != 0) {
    vp9_wb_write_bit(wb, 1);
    vp9_wb_write_literal(wb, abs(delta_q), 4);
    vp9_wb_write_bit(wb, delta_q < 0);
  } else {
    vp9_wb_write_bit(wb, 0);
  }
}

static void encode_quantization(VP9_COMMON *cm,
                                struct vp9_write_bit_buffer *wb) {
  vp9_wb_write_literal(wb, cm->base_qindex, QINDEX_BITS);
  write_delta_q(wb, cm->y_dc_delta_q);
  write_delta_q(wb, cm->uv_dc_delta_q);
  write_delta_q(wb, cm->uv_ac_delta_q);
}


static void encode_segmentation(VP9_COMP *cpi,
                               struct vp9_write_bit_buffer *wb) {
  int i, j;
  VP9_COMMON *const cm = &cpi->common;
  MACROBLOCKD *const xd = &cpi->mb.e_mbd;

  vp9_wb_write_bit(wb, xd->segmentation_enabled);
  if (!xd->segmentation_enabled)
    return;

  // Segmentation map
  vp9_wb_write_bit(wb, xd->update_mb_segmentation_map);
  if (xd->update_mb_segmentation_map) {
    // Select the coding strategy (temporal or spatial)
    vp9_choose_segmap_coding_method(cpi);
    // Write out probabilities used to decode unpredicted  macro-block segments
    for (i = 0; i < MB_SEG_TREE_PROBS; i++) {
      const int prob = xd->mb_segment_tree_probs[i];
      const int update = prob != MAX_PROB;
      vp9_wb_write_bit(wb, update);
      if (update)
        vp9_wb_write_literal(wb, prob, 8);
    }

    // Write out the chosen coding method.
    vp9_wb_write_bit(wb, cm->temporal_update);
    if (cm->temporal_update) {
      for (i = 0; i < PREDICTION_PROBS; i++) {
        const int prob = cm->segment_pred_probs[i];
        const int update = prob != MAX_PROB;
        vp9_wb_write_bit(wb, update);
        if (update)
          vp9_wb_write_literal(wb, prob, 8);
      }
    }
  }

  // Segmentation data
  vp9_wb_write_bit(wb, xd->update_mb_segmentation_data);
  if (xd->update_mb_segmentation_data) {
    vp9_wb_write_bit(wb, xd->mb_segment_abs_delta);

    for (i = 0; i < MAX_MB_SEGMENTS; i++) {
      for (j = 0; j < SEG_LVL_MAX; j++) {
        const int active = vp9_segfeature_active(xd, i, j);
        vp9_wb_write_bit(wb, active);
        if (active) {
          const int data = vp9_get_segdata(xd, i, j);
          const int data_max = vp9_seg_feature_data_max(j);

          if (vp9_is_segfeature_signed(j)) {
            vp9_encode_unsigned_max(wb, abs(data), data_max);
            vp9_wb_write_bit(wb, data < 0);
          } else {
            vp9_encode_unsigned_max(wb, data, data_max);
          }
        }
      }
    }
  }
}


static void encode_txfm_probs(VP9_COMP *cpi, vp9_writer *w) {
  VP9_COMMON *const cm = &cpi->common;

  // Mode
  vp9_write_literal(w, MIN(cm->txfm_mode, ALLOW_32X32), 2);
  if (cm->txfm_mode >= ALLOW_32X32)
    vp9_write_bit(w, cm->txfm_mode == TX_MODE_SELECT);

  // Probabilities
  if (cm->txfm_mode == TX_MODE_SELECT) {
    int i, j;
    unsigned int ct_8x8p[TX_SIZE_MAX_SB - 3][2];
    unsigned int ct_16x16p[TX_SIZE_MAX_SB - 2][2];
    unsigned int ct_32x32p[TX_SIZE_MAX_SB - 1][2];


    for (i = 0; i < TX_SIZE_CONTEXTS; i++) {
      tx_counts_to_branch_counts_8x8(cm->fc.tx_count_8x8p[i],
                                     ct_8x8p);
      for (j = 0; j < TX_SIZE_MAX_SB - 3; j++) {
        vp9_cond_prob_diff_update(w, &cm->fc.tx_probs_8x8p[i][j],
                                  VP9_MODE_UPDATE_PROB, ct_8x8p[j]);
      }
    }
    for (i = 0; i < TX_SIZE_CONTEXTS; i++) {
      tx_counts_to_branch_counts_16x16(cm->fc.tx_count_16x16p[i],
                                       ct_16x16p);
      for (j = 0; j < TX_SIZE_MAX_SB - 2; j++) {
        vp9_cond_prob_diff_update(w, &cm->fc.tx_probs_16x16p[i][j],
                                  VP9_MODE_UPDATE_PROB, ct_16x16p[j]);
      }
    }
    for (i = 0; i < TX_SIZE_CONTEXTS; i++) {
      tx_counts_to_branch_counts_32x32(cm->fc.tx_count_32x32p[i],
                                       ct_32x32p);
      for (j = 0; j < TX_SIZE_MAX_SB - 1; j++) {
        vp9_cond_prob_diff_update(w, &cm->fc.tx_probs_32x32p[i][j],
                                  VP9_MODE_UPDATE_PROB, ct_32x32p[j]);
      }
    }
#ifdef MODE_STATS
    if (!cpi->dummy_packing)
      update_tx_count_stats(cm);
#endif
  }
}

static void write_interp_filter_type(INTERPOLATIONFILTERTYPE type,
                                     struct vp9_write_bit_buffer *wb) {
  vp9_wb_write_bit(wb, type == SWITCHABLE);
  if (type != SWITCHABLE)
    vp9_wb_write_literal(wb, type, 2);
}

static void fix_mcomp_filter_type(VP9_COMP *cpi) {
  VP9_COMMON *const cm = &cpi->common;

  if (cm->mcomp_filter_type == SWITCHABLE) {
    // Check to see if only one of the filters is actually used
    int count[VP9_SWITCHABLE_FILTERS];
    int i, j, c = 0;
    for (i = 0; i < VP9_SWITCHABLE_FILTERS; ++i) {
      count[i] = 0;
      for (j = 0; j <= VP9_SWITCHABLE_FILTERS; ++j)
        count[i] += cm->fc.switchable_interp_count[j][i];
      c += (count[i] > 0);
    }
    if (c == 1) {
      // Only one filter is used. So set the filter at frame level
      for (i = 0; i < VP9_SWITCHABLE_FILTERS; ++i) {
        if (count[i]) {
          cm->mcomp_filter_type = vp9_switchable_interp[i];
          break;
        }
      }
    }
  }
}

static void write_tile_info(VP9_COMMON *cm, struct vp9_write_bit_buffer *wb) {
  int min_log2_tiles, delta_log2_tiles, n_tile_bits, n;
  vp9_get_tile_n_bits(cm, &min_log2_tiles, &delta_log2_tiles);
  n_tile_bits = cm->log2_tile_columns - min_log2_tiles;
  for (n = 0; n < delta_log2_tiles; n++) {
    if (n_tile_bits--) {
      vp9_wb_write_bit(wb, 1);
    } else {
      vp9_wb_write_bit(wb, 0);
      break;
    }
  }

  vp9_wb_write_bit(wb, cm->log2_tile_rows != 0);
  if (cm->log2_tile_rows != 0)
    vp9_wb_write_bit(wb, cm->log2_tile_rows != 1);
}

static int get_refresh_mask(VP9_COMP *cpi) {
    // Should the GF or ARF be updated using the transmitted frame or buffer
#if CONFIG_MULTIPLE_ARF
    if (!cpi->multi_arf_enabled && cpi->refresh_golden_frame &&
        !cpi->refresh_alt_ref_frame) {
#else
    if (cpi->refresh_golden_frame && !cpi->refresh_alt_ref_frame) {
#endif
      // Preserve the previously existing golden frame and update the frame in
      // the alt ref slot instead. This is highly specific to the use of
      // alt-ref as a forward reference, and this needs to be generalized as
      // other uses are implemented (like RTC/temporal scaling)
      //
      // gld_fb_idx and alt_fb_idx need to be swapped for future frames, but
      // that happens in vp9_onyx_if.c:update_reference_frames() so that it can
      // be done outside of the recode loop.
      return (cpi->refresh_last_frame << cpi->lst_fb_idx) |
             (cpi->refresh_golden_frame << cpi->alt_fb_idx);
    } else {
      int arf_idx = cpi->alt_fb_idx;
#if CONFIG_MULTIPLE_ARF
      // Determine which ARF buffer to use to encode this ARF frame.
      if (cpi->multi_arf_enabled) {
        int sn = cpi->sequence_number;
        arf_idx = (cpi->frame_coding_order[sn] < 0) ?
            cpi->arf_buffer_idx[sn + 1] :
            cpi->arf_buffer_idx[sn];
      }
#endif
      return (cpi->refresh_last_frame << cpi->lst_fb_idx) |
             (cpi->refresh_golden_frame << cpi->gld_fb_idx) |
             (cpi->refresh_alt_ref_frame << arf_idx);
    }
}

static void write_display_size(VP9_COMP *cpi, struct vp9_write_bit_buffer *wb) {
  VP9_COMMON *const cm = &cpi->common;

  const int scaling_active = cm->width != cm->display_width ||
                             cm->height != cm->display_height;
  vp9_wb_write_bit(wb, scaling_active);
  if (scaling_active) {
    vp9_wb_write_literal(wb, cm->display_width - 1, 16);
    vp9_wb_write_literal(wb, cm->display_height - 1, 16);
  }
}

static void write_frame_size(VP9_COMP *cpi,
                             struct vp9_write_bit_buffer *wb) {
  VP9_COMMON *const cm = &cpi->common;
  vp9_wb_write_literal(wb, cm->width - 1, 16);
  vp9_wb_write_literal(wb, cm->height - 1, 16);

  write_display_size(cpi, wb);
}

static void write_frame_size_with_refs(VP9_COMP *cpi,
                                       struct vp9_write_bit_buffer *wb) {
  VP9_COMMON *const cm = &cpi->common;
  int refs[ALLOWED_REFS_PER_FRAME] = {cpi->lst_fb_idx, cpi->gld_fb_idx,
                                      cpi->alt_fb_idx};
  int i, found = 0;

  for (i = 0; i < ALLOWED_REFS_PER_FRAME; ++i) {
    YV12_BUFFER_CONFIG *cfg = &cm->yv12_fb[cm->ref_frame_map[refs[i]]];
    found = cm->width == cfg->y_crop_width &&
            cm->height == cfg->y_crop_height;
    vp9_wb_write_bit(wb, found);
    if (found)
      break;
  }

  if (!found) {
    vp9_wb_write_literal(wb, cm->width - 1, 16);
    vp9_wb_write_literal(wb, cm->height - 1, 16);
  }

  write_display_size(cpi, wb);
}

static void write_sync_code(struct vp9_write_bit_buffer *wb) {
  vp9_wb_write_literal(wb, SYNC_CODE_0, 8);
  vp9_wb_write_literal(wb, SYNC_CODE_1, 8);
  vp9_wb_write_literal(wb, SYNC_CODE_2, 8);
}

static void write_uncompressed_header(VP9_COMP *cpi,
                                      struct vp9_write_bit_buffer *wb) {
  VP9_COMMON *const cm = &cpi->common;
  MACROBLOCKD *const xd = &cpi->mb.e_mbd;

  // frame marker bits
  vp9_wb_write_literal(wb, 0x2, 2);

  // bitstream version.
  // 00 - profile 0. 4:2:0 only
  // 10 - profile 1. adds 4:4:4, 4:2:2, alpha
  vp9_wb_write_bit(wb, cm->version);
  vp9_wb_write_bit(wb, 0);

  vp9_wb_write_bit(wb, 0);
  vp9_wb_write_bit(wb, cm->frame_type);
  vp9_wb_write_bit(wb, cm->show_frame);
  vp9_wb_write_bit(wb, cm->error_resilient_mode);

  if (cm->frame_type == KEY_FRAME) {
    write_sync_code(wb);
    // colorspaces
    // 000 - Unknown
    // 001 - BT.601
    // 010 - BT.709
    // 011 - SMPTE-170
    // 100 - SMPTE-240
    // 101 - Reserved
    // 110 - Reserved
    // 111 - sRGB (RGB)
    vp9_wb_write_literal(wb, 0, 3);
    if (1 /* colorspace != sRGB */) {
      vp9_wb_write_bit(wb, 0);  // 0: [16, 235] (i.e. xvYCC), 1: [0, 255]
      if (cm->version == 1) {
        vp9_wb_write_bit(wb, cm->subsampling_x);
        vp9_wb_write_bit(wb, cm->subsampling_y);
        vp9_wb_write_bit(wb, 0);  // has extra plane
      }
    } else {
      assert(cm->version == 1);
      vp9_wb_write_bit(wb, 0);  // has extra plane
    }

    write_frame_size(cpi, wb);
  } else {
    const int refs[ALLOWED_REFS_PER_FRAME] = {cpi->lst_fb_idx, cpi->gld_fb_idx,
                                              cpi->alt_fb_idx};
    if (!cm->show_frame)
      vp9_wb_write_bit(wb, cm->intra_only);

    if (!cm->error_resilient_mode)
      vp9_wb_write_literal(wb, cm->reset_frame_context, 2);

    if (cm->intra_only) {
      write_sync_code(wb);

      vp9_wb_write_literal(wb, get_refresh_mask(cpi), NUM_REF_FRAMES);
      write_frame_size(cpi, wb);
    } else {
      int i;
      vp9_wb_write_literal(wb, get_refresh_mask(cpi), NUM_REF_FRAMES);
      for (i = 0; i < ALLOWED_REFS_PER_FRAME; ++i) {
        vp9_wb_write_literal(wb, refs[i], NUM_REF_FRAMES_LG2);
        vp9_wb_write_bit(wb, cm->ref_frame_sign_bias[LAST_FRAME + i]);
      }

      write_frame_size_with_refs(cpi, wb);

      vp9_wb_write_bit(wb, xd->allow_high_precision_mv);

      fix_mcomp_filter_type(cpi);
      write_interp_filter_type(cm->mcomp_filter_type, wb);
    }
  }

  if (!cm->error_resilient_mode) {
    vp9_wb_write_bit(wb, cm->refresh_frame_context);
    vp9_wb_write_bit(wb, cm->frame_parallel_decoding_mode);
  }

  vp9_wb_write_literal(wb, cm->frame_context_idx, NUM_FRAME_CONTEXTS_LG2);

  encode_loopfilter(cm, xd, wb);
  encode_quantization(cm, wb);
  encode_segmentation(cpi, wb);

  write_tile_info(cm, wb);
}

void vp9_pack_bitstream(VP9_COMP *cpi, uint8_t *dest, unsigned long *size) {
  int i, bytes_packed;
  VP9_COMMON *const pc = &cpi->common;
  vp9_writer header_bc, residual_bc;
  MACROBLOCKD *const xd = &cpi->mb.e_mbd;

  uint8_t *cx_data = dest;
  struct vp9_write_bit_buffer wb = {dest, 0};
  struct vp9_write_bit_buffer first_partition_size_wb;

  write_uncompressed_header(cpi, &wb);
  first_partition_size_wb = wb;
  vp9_wb_write_literal(&wb, 0, 16);  // don't know in advance first part. size

  bytes_packed = vp9_rb_bytes_written(&wb);
  cx_data += bytes_packed;

  compute_update_table();

  vp9_start_encode(&header_bc, cx_data);

#ifdef ENTROPY_STATS
  if (pc->frame_type == INTER_FRAME)
    active_section = 0;
  else
    active_section = 7;
#endif

  vp9_clear_system_state();  // __asm emms;

  vp9_copy(pc->fc.pre_coef_probs, pc->fc.coef_probs);
  vp9_copy(pc->fc.pre_y_mode_prob, pc->fc.y_mode_prob);
  vp9_copy(pc->fc.pre_uv_mode_prob, pc->fc.uv_mode_prob);
  vp9_copy(pc->fc.pre_partition_prob, pc->fc.partition_prob[INTER_FRAME]);
  pc->fc.pre_nmvc = pc->fc.nmvc;
  vp9_copy(pc->fc.pre_switchable_interp_prob, pc->fc.switchable_interp_prob);
  vp9_copy(pc->fc.pre_inter_mode_probs, pc->fc.inter_mode_probs);
  vp9_copy(pc->fc.pre_intra_inter_prob, pc->fc.intra_inter_prob);
  vp9_copy(pc->fc.pre_comp_inter_prob, pc->fc.comp_inter_prob);
  vp9_copy(pc->fc.pre_comp_ref_prob, pc->fc.comp_ref_prob);
  vp9_copy(pc->fc.pre_single_ref_prob, pc->fc.single_ref_prob);
  vp9_copy(pc->fc.pre_tx_probs_8x8p, pc->fc.tx_probs_8x8p);
  vp9_copy(pc->fc.pre_tx_probs_16x16p, pc->fc.tx_probs_16x16p);
  vp9_copy(pc->fc.pre_tx_probs_32x32p, pc->fc.tx_probs_32x32p);
  vp9_copy(pc->fc.pre_mbskip_probs, pc->fc.mbskip_probs);

  if (xd->lossless) {
    pc->txfm_mode = ONLY_4X4;
  } else {
    encode_txfm_probs(cpi, &header_bc);
  }

  update_coef_probs(cpi, &header_bc);

#ifdef ENTROPY_STATS
  active_section = 2;
#endif

  vp9_update_skip_probs(cpi, &header_bc);

  if (pc->frame_type != KEY_FRAME) {
#ifdef ENTROPY_STATS
    active_section = 1;
#endif

    update_inter_mode_probs(pc, &header_bc);
    vp9_zero(cpi->common.fc.inter_mode_counts);

    if (pc->mcomp_filter_type == SWITCHABLE)
      update_switchable_interp_probs(cpi, &header_bc);

    for (i = 0; i < INTRA_INTER_CONTEXTS; i++)
      vp9_cond_prob_diff_update(&header_bc, &pc->fc.intra_inter_prob[i],
                                VP9_MODE_UPDATE_PROB,
                                cpi->intra_inter_count[i]);

    if (pc->allow_comp_inter_inter) {
      const int comp_pred_mode = cpi->common.comp_pred_mode;
      const int use_compound_pred = (comp_pred_mode != SINGLE_PREDICTION_ONLY);
      const int use_hybrid_pred = (comp_pred_mode == HYBRID_PREDICTION);

      vp9_write_bit(&header_bc, use_compound_pred);
      if (use_compound_pred) {
        vp9_write_bit(&header_bc, use_hybrid_pred);
        if (use_hybrid_pred) {
          for (i = 0; i < COMP_INTER_CONTEXTS; i++)
            vp9_cond_prob_diff_update(&header_bc, &pc->fc.comp_inter_prob[i],
                                      VP9_MODE_UPDATE_PROB,
                                      cpi->comp_inter_count[i]);
        }
      }
    }

    if (pc->comp_pred_mode != COMP_PREDICTION_ONLY) {
      for (i = 0; i < REF_CONTEXTS; i++) {
        vp9_cond_prob_diff_update(&header_bc, &pc->fc.single_ref_prob[i][0],
                                  VP9_MODE_UPDATE_PROB,
                                  cpi->single_ref_count[i][0]);
        vp9_cond_prob_diff_update(&header_bc, &pc->fc.single_ref_prob[i][1],
                                  VP9_MODE_UPDATE_PROB,
                                  cpi->single_ref_count[i][1]);
      }
    }

    if (pc->comp_pred_mode != SINGLE_PREDICTION_ONLY) {
      for (i = 0; i < REF_CONTEXTS; i++)
        vp9_cond_prob_diff_update(&header_bc, &pc->fc.comp_ref_prob[i],
                                  VP9_MODE_UPDATE_PROB,
                                  cpi->comp_ref_count[i]);
    }

    update_mbintra_mode_probs(cpi, &header_bc);

    for (i = 0; i < NUM_PARTITION_CONTEXTS; ++i) {
      vp9_prob Pnew[PARTITION_TYPES - 1];
      unsigned int bct[PARTITION_TYPES - 1][2];
      update_mode(&header_bc, PARTITION_TYPES, vp9_partition_encodings,
                  vp9_partition_tree, Pnew,
                  pc->fc.partition_prob[pc->frame_type][i], bct,
                  (unsigned int *)cpi->partition_count[i]);
    }

    vp9_write_nmv_probs(cpi, xd->allow_high_precision_mv, &header_bc);
  }


  vp9_stop_encode(&header_bc);


  // first partition size
  assert(header_bc.pos <= 0xffff);
  vp9_wb_write_literal(&first_partition_size_wb, header_bc.pos, 16);
  *size = bytes_packed + header_bc.pos;

  {
    int tile_row, tile_col, total_size = 0;
    unsigned char *data_ptr = cx_data + header_bc.pos;
    TOKENEXTRA *tok[4][1 << 6], *tok_end;

    vpx_memset(cpi->common.above_seg_context, 0, sizeof(PARTITION_CONTEXT) *
               mi_cols_aligned_to_sb(&cpi->common));
    tok[0][0] = cpi->tok;
    for (tile_row = 0; tile_row < pc->tile_rows; tile_row++) {
      if (tile_row) {
        tok[tile_row][0] = tok[tile_row - 1][pc->tile_columns - 1] +
                           cpi->tok_count[tile_row - 1][pc->tile_columns - 1];
      }
      for (tile_col = 1; tile_col < pc->tile_columns; tile_col++) {
        tok[tile_row][tile_col] = tok[tile_row][tile_col - 1] +
                                  cpi->tok_count[tile_row][tile_col - 1];
      }
    }

    for (tile_row = 0; tile_row < pc->tile_rows; tile_row++) {
      vp9_get_tile_row_offsets(pc, tile_row);
      for (tile_col = 0; tile_col < pc->tile_columns; tile_col++) {
        vp9_get_tile_col_offsets(pc, tile_col);
        tok_end = tok[tile_row][tile_col] + cpi->tok_count[tile_row][tile_col];

        if (tile_col < pc->tile_columns - 1 || tile_row < pc->tile_rows - 1)
          vp9_start_encode(&residual_bc, data_ptr + total_size + 4);
        else
          vp9_start_encode(&residual_bc, data_ptr + total_size);
        write_modes(cpi, &residual_bc, &tok[tile_row][tile_col], tok_end);
        assert(tok[tile_row][tile_col] == tok_end);
        vp9_stop_encode(&residual_bc);
        if (tile_col < pc->tile_columns - 1 || tile_row < pc->tile_rows - 1) {
          // size of this tile
          write_be32(data_ptr + total_size, residual_bc.pos);
          total_size += 4;
        }

        total_size += residual_bc.pos;
      }
    }

    *size += total_size;
  }
}

#ifdef ENTROPY_STATS
static void print_tree_update_for_type(FILE *f,
                                       vp9_coeff_stats *tree_update_hist,
                                       int block_types, const char *header) {
  int i, j, k, l, m;

  fprintf(f, "const vp9_coeff_prob %s = {\n", header);
  for (i = 0; i < block_types; i++) {
    fprintf(f, "  { \n");
    for (j = 0; j < REF_TYPES; j++) {
      fprintf(f, "  { \n");
      for (k = 0; k < COEF_BANDS; k++) {
        fprintf(f, "    {\n");
        for (l = 0; l < PREV_COEF_CONTEXTS; l++) {
          fprintf(f, "      {");
          for (m = 0; m < ENTROPY_NODES; m++) {
            fprintf(f, "%3d, ",
                    get_binary_prob(tree_update_hist[i][j][k][l][m][0],
                                    tree_update_hist[i][j][k][l][m][1]));
          }
          fprintf(f, "},\n");
        }
        fprintf(f, "},\n");
      }
      fprintf(f, "    },\n");
    }
    fprintf(f, "  },\n");
  }
  fprintf(f, "};\n");
}

void print_tree_update_probs() {
  FILE *f = fopen("coefupdprob.h", "w");
  fprintf(f, "\n/* Update probabilities for token entropy tree. */\n\n");

  print_tree_update_for_type(f, tree_update_hist[TX_4X4],   BLOCK_TYPES,
                             "vp9_coef_update_probs_4x4[BLOCK_TYPES]");
  print_tree_update_for_type(f, tree_update_hist[TX_8X8],   BLOCK_TYPES,
                             "vp9_coef_update_probs_8x8[BLOCK_TYPES]");
  print_tree_update_for_type(f, tree_update_hist[TX_16X16], BLOCK_TYPES,
                             "vp9_coef_update_probs_16x16[BLOCK_TYPES]");
  print_tree_update_for_type(f, tree_update_hist[TX_32X32], BLOCK_TYPES,
                             "vp9_coef_update_probs_32x32[BLOCK_TYPES]");

  fclose(f);
  f = fopen("treeupdate.bin", "wb");
  fwrite(tree_update_hist, sizeof(tree_update_hist), 1, f);
  fclose(f);
}
#endif
