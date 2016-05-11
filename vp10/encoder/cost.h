/*
 *  Copyright (c) 2014 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#ifndef VP10_ENCODER_COST_H_
#define VP10_ENCODER_COST_H_

#include "vpx_dsp/prob.h"
#include "vpx/vpx_integer.h"
#if CONFIG_ANS
#include "vp10/common/ans.h"
#endif  // CONFIG_ANS

#ifdef __cplusplus
extern "C" {
#endif

extern const uint16_t vp10_prob_cost[256];

// The factor to scale from cost in bits to cost in vp10_prob_cost units.
#define VP9_PROB_COST_SHIFT 9

#define vp10_cost_zero(prob) (vp10_prob_cost[prob])

#define vp10_cost_one(prob) vp10_cost_zero(256 - (prob))

#define vp10_cost_bit(prob, bit) vp10_cost_zero((bit) ? 256 - (prob) \
                                                    : (prob))

// Cost of coding an n bit literal, using 128 (i.e. 50%) probability
// for each bit.
#define vp10_cost_literal(n) ((n) * (1 << VP9_PROB_COST_SHIFT))

static INLINE unsigned int cost_branch256(const unsigned int ct[2],
                                          vpx_prob p) {
  return ct[0] * vp10_cost_zero(p) + ct[1] * vp10_cost_one(p);
}

static INLINE int treed_cost(vpx_tree tree, const vpx_prob *probs,
                             int bits, int len) {
  int cost = 0;
  vpx_tree_index i = 0;

  do {
    const int bit = (bits >> --len) & 1;
    cost += vp10_cost_bit(probs[i >> 1], bit);
    i = tree[i + bit];
  } while (len);

  return cost;
}

void vp10_cost_tokens(int *costs, const vpx_prob *probs, vpx_tree tree);
void vp10_cost_tokens_skip(int *costs, const vpx_prob *probs, vpx_tree tree);

#if CONFIG_ANS
void vp10_cost_tokens_ans(int *costs, const vpx_prob *tree_probs,
                          const rans_dec_lut token_cdf, int skip_eob);
#endif

#ifdef __cplusplus
}  // extern "C"
#endif

#endif  // VP10_ENCODER_COST_H_
