/*
 * Copyright (c) 2016, Alliance for Open Media. All rights reserved
 *
 * This source code is subject to the terms of the BSD 2 Clause License and
 * the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
 * was not distributed with this source code in the LICENSE file, you can
 * obtain it at www.aomedia.org/license/software. If the Alliance for Open
 * Media Patent License 1.0 was not distributed with this source code in the
 * PATENTS file, you can obtain it at www.aomedia.org/license/patent.
 */

#include <assert.h>
#include "./aom_config.h"
#include "aom/aom_integer.h"
#include "aom_dsp/ans.h"
#include "aom_dsp/prob.h"

static int find_largest(const aom_cdf_prob *const pdf_tab, int num_syms) {
  int largest_idx = -1;
  int largest_p = -1;
  int i;
  for (i = 0; i < num_syms; ++i) {
    int p = pdf_tab[i];
    if (p > largest_p) {
      largest_p = p;
      largest_idx = i;
    }
  }
  return largest_idx;
}

void aom_rans_merge_prob8_pdf(aom_cdf_prob *const out_pdf,
                              const AnsP8 node_prob,
                              const aom_cdf_prob *const src_pdf, int in_syms) {
  int i;
  int adjustment = RANS_PRECISION;
  const int round_fact = ANS_P8_PRECISION >> 1;
  const AnsP8 p1 = ANS_P8_PRECISION - node_prob;
  const int out_syms = in_syms + 1;
  assert(src_pdf != out_pdf);

  out_pdf[0] = node_prob << (RANS_PROB_BITS - ANS_P8_SHIFT);
  adjustment -= out_pdf[0];
  for (i = 0; i < in_syms; ++i) {
    int p = (p1 * src_pdf[i] + round_fact) >> ANS_P8_SHIFT;
    p = AOMMIN(p, (int)RANS_PRECISION - in_syms);
    p = AOMMAX(p, 1);
    out_pdf[i + 1] = p;
    adjustment -= p;
  }

  // Adjust probabilities so they sum to the total probability
  if (adjustment > 0) {
    i = find_largest(out_pdf, out_syms);
    out_pdf[i] += adjustment;
  } else {
    while (adjustment < 0) {
      i = find_largest(out_pdf, out_syms);
      --out_pdf[i];
      assert(out_pdf[i] > 0);
      adjustment++;
    }
  }
}
