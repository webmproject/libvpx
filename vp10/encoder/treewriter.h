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

#ifdef VP10_FORCE_VPXBOOL_TREEWRITER
#include "vpx_dsp/bitwriter.h"
#define tree_writer vpx_writer
#define tree_bit_write vpx_write
#else
#include "vp10/encoder/bitwriter.h"
#define tree_writer vp10_writer
#define tree_bit_write vp10_write
#endif

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

static INLINE void vp10_write_tree(tree_writer *w, const vpx_tree_index *tree,
                                  const vpx_prob *probs, int bits, int len,
                                  vpx_tree_index i) {
  do {
    const int bit = (bits >> --len) & 1;
    tree_bit_write(w, bit, probs[i >> 1]);
    i = tree[i + bit];
  } while (len);
}

static INLINE void vp10_write_token(tree_writer *w, const vpx_tree_index *tree,
                                   const vpx_prob *probs,
                                   const struct vp10_token *token) {
  vp10_write_tree(w, tree, probs, token->value, token->len, 0);
}

#undef tree_writer
#undef tree_bit_write
#ifdef __cplusplus
}  // extern "C"
#endif

#endif  // VP10_ENCODER_TREEWRITER_H_
