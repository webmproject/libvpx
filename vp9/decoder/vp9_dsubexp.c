/*
  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include "vp9/common/vp9_entropy.h"

#include "vp9/decoder/vp9_dsubexp.h"

static int inv_recenter_nonneg(int v, int m) {
  if (v > 2 * m)
    return v;

  return v % 2 ? m - (v + 1) / 2 : m + v / 2;
}

static int decode_uniform(vp9_reader *r, int n) {
  int v;
  const int l = get_unsigned_bits(n);
  const int m = (1 << l) - n;
  if (!l)
    return 0;

  v = vp9_read_literal(r, l - 1);
  return v < m ?  v : (v << 1) - m + vp9_read_bit(r);
}


static int merge_index(int v, int n, int modulus) {
  int max1 = (n - 1 - modulus / 2) / modulus + 1;
  if (v < max1) {
    v = v * modulus + modulus / 2;
  } else {
    int w;
    v -= max1;
    w = v;
    v += (v + modulus - modulus / 2) / modulus;
    while (v % modulus == modulus / 2 ||
           w != v - (v + modulus - modulus / 2) / modulus) v++;
  }
  return v;
}

static int inv_remap_prob(int v, int m) {
  const int n = 255;

  v = merge_index(v, n - 1, MODULUS_PARAM);
  m--;
  if ((m << 1) <= n) {
    return 1 + inv_recenter_nonneg(v + 1, m);
  } else {
    return n - inv_recenter_nonneg(v + 1, n - 1 - m);
  }
}

static int decode_term_subexp(vp9_reader *r, int k, int num_syms) {
  int i = 0, mk = 0, word;
  while (1) {
    const int b = i ? k + i - 1 : k;
    const int a = 1 << b;
    if (num_syms <= mk + 3 * a) {
      word = decode_uniform(r, num_syms - mk) + mk;
      break;
    } else {
      if (vp9_read_bit(r)) {
        i++;
        mk += a;
      } else {
        word = vp9_read_literal(r, b) + mk;
        break;
      }
    }
  }
  return word;
}

void vp9_diff_update_prob(vp9_reader *r, vp9_prob* p) {
  int delp = decode_term_subexp(r, SUBEXP_PARAM, 255);
  *p = (vp9_prob)inv_remap_prob(delp, *p);
}
