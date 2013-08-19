/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#ifndef VP9_COMMON_VP9_ENTROPYMODE_H_
#define VP9_COMMON_VP9_ENTROPYMODE_H_

#include "vp9/common/vp9_blockd.h"
#include "vp9/common/vp9_treecoder.h"

#define SUBMVREF_COUNT 5
#define TX_SIZE_CONTEXTS 2
#define VP9_MODE_UPDATE_PROB  252
#define VP9_SWITCHABLE_FILTERS 3   // number of switchable filters

#if CONFIG_INTERINTRA
#define VP9_UPD_INTERINTRA_PROB 248
#define SEPARATE_INTERINTRA_UV  0
#if CONFIG_MASKED_COMPOUND
#define VP9_UPD_MASKED_INTERINTRA_PROB 248
#endif
#endif

#if CONFIG_MASKED_COMPOUND
#define VP9_UPD_MASKED_COMPOUND_PROB 248
#endif

// #define MODE_STATS

struct VP9Common;

struct tx_probs {
  vp9_prob p32x32[TX_SIZE_CONTEXTS][TX_SIZES - 1];
  vp9_prob p16x16[TX_SIZE_CONTEXTS][TX_SIZES - 2];
  vp9_prob p8x8[TX_SIZE_CONTEXTS][TX_SIZES - 3];
};

struct tx_counts {
  unsigned int p32x32[TX_SIZE_CONTEXTS][TX_SIZES];
  unsigned int p16x16[TX_SIZE_CONTEXTS][TX_SIZES - 1];
  unsigned int p8x8[TX_SIZE_CONTEXTS][TX_SIZES - 2];
};

extern const vp9_prob vp9_kf_uv_mode_prob[VP9_INTRA_MODES][VP9_INTRA_MODES - 1];
extern const vp9_prob vp9_kf_y_mode_prob[VP9_INTRA_MODES][VP9_INTRA_MODES]
                                        [VP9_INTRA_MODES - 1];

extern const vp9_tree_index vp9_intra_mode_tree[];
extern const vp9_tree_index vp9_inter_mode_tree[];

extern struct vp9_token vp9_intra_mode_encodings[VP9_INTRA_MODES];
extern struct vp9_token vp9_inter_mode_encodings[VP9_INTER_MODES];

// probability models for partition information
extern const vp9_tree_index vp9_partition_tree[];
extern struct vp9_token vp9_partition_encodings[PARTITION_TYPES];

extern const vp9_tree_index vp9_switchable_interp_tree
                 [2 * (VP9_SWITCHABLE_FILTERS - 1)];

extern struct vp9_token vp9_switchable_interp_encodings[VP9_SWITCHABLE_FILTERS];

void vp9_entropy_mode_init();

void vp9_setup_past_independence(struct VP9Common *cm, MACROBLOCKD *xd);

void vp9_init_mbmode_probs(struct VP9Common *x);

void vp9_adapt_mode_probs(struct VP9Common *);

void tx_counts_to_branch_counts_32x32(unsigned int *tx_count_32x32p,
                                      unsigned int (*ct_32x32p)[2]);
void tx_counts_to_branch_counts_16x16(unsigned int *tx_count_16x16p,
                                      unsigned int (*ct_16x16p)[2]);
void tx_counts_to_branch_counts_8x8(unsigned int *tx_count_8x8p,
                                    unsigned int (*ct_8x8p)[2]);

#endif  // VP9_COMMON_VP9_ENTROPYMODE_H_
