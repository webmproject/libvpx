/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#ifndef VP10_ENCODER_TREEWRITER_H_
#define VP10_ENCODER_TREEWRITER_H_

#include "vpx_dsp/bitwriter.h"
#include "vp10/common/ans.h"

#ifdef __cplusplus
extern "C" {
#endif

void vp10_tree_probs_from_distribution(vpx_tree tree,
                                      unsigned int branch_ct[ /* n - 1 */ ][2],
                                      const unsigned int num_events[ /* n */ ]);

struct vp10_token {
  int value;
  int len;
};

void vp10_tokens_from_tree(struct vp10_token*, const vpx_tree_index *);

// TODO: CHECK MAX REVERSIBLE TREE SIZE, security concerns
#define VP10_TOKEN_SCRATCH_LEN 32
static INLINE void vp10_write_tree_r(struct AnsCoder *const ans,
                                     const vpx_tree_index *const tree,
                                     const vpx_prob *const probs,
                                     int bits, int len,
                                     vpx_tree_index tidx) {
  int i;
  struct { uint8_t bit; vpx_prob prob; } scratch[VP10_TOKEN_SCRATCH_LEN];
  // assert(len < VP10_TOKEN_SCRATCH_LEN);

  for (i = len - 1; i >= 0; --i) {
    const int bit = (bits >> i) & 1;
    scratch[i].bit = bit;
    scratch[i].prob = probs[tidx >> 1];
    tidx = tree[tidx + bit];
  }
  for (i = 0; i < len; ++i) {
    rabs_write(ans, scratch[i].bit, scratch[i].prob);
  }
}
#undef VP10_TOKEN_SCRATCH_LEN

static INLINE void vp10_write_tree(vpx_writer *w, const vpx_tree_index *tree,
                                  const vpx_prob *probs, int bits, int len,
                                  vpx_tree_index i) {
  do {
    const int bit = (bits >> --len) & 1;
    vpx_write(w, bit, probs[i >> 1]);
    i = tree[i + bit];
  } while (len);
}

static INLINE void vp10_write_token(vpx_writer *w, const vpx_tree_index *tree,
                                   const vpx_prob *probs,
                                   const struct vp10_token *token) {
  vp10_write_tree(w, tree, probs, token->value, token->len, 0);
}

#ifdef __cplusplus
}  // extern "C"
#endif

#endif  // VP10_ENCODER_TREEWRITER_H_
