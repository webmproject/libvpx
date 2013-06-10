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

// #define MODE_STATS

extern int vp9_mv_cont(const int_mv *l, const int_mv *a);


extern const vp9_prob vp9_kf_default_bmode_probs[VP9_INTRA_MODES]
                                                [VP9_INTRA_MODES]
                                                [VP9_INTRA_MODES - 1];

extern const vp9_tree_index vp9_intra_mode_tree[];
extern const vp9_tree_index  vp9_sb_mv_ref_tree[];

extern struct vp9_token vp9_intra_mode_encodings[VP9_INTRA_MODES];

/* Inter mode values do not start at zero */

extern struct vp9_token vp9_sb_mv_ref_encoding_array[VP9_INTER_MODES];

// probability models for partition information
extern const vp9_tree_index  vp9_partition_tree[];
extern struct vp9_token vp9_partition_encodings[PARTITION_TYPES];
extern const vp9_prob vp9_partition_probs[NUM_FRAME_TYPES]
                                         [NUM_PARTITION_CONTEXTS]
                                         [PARTITION_TYPES - 1];

void vp9_entropy_mode_init(void);

struct VP9Common;

/* sets up common features to forget past dependence */
void vp9_setup_past_independence(struct VP9Common *cm, MACROBLOCKD *xd);

void vp9_init_mbmode_probs(struct VP9Common *x);

extern void vp9_init_mode_contexts(struct VP9Common *pc);

extern void vp9_adapt_mode_context(struct VP9Common *pc);

extern void vp9_accum_mv_refs(struct VP9Common *pc,
                              MB_PREDICTION_MODE m,
                              const int context);

void vp9_adapt_mode_probs(struct VP9Common *);

#define VP9_SWITCHABLE_FILTERS 3 /* number of switchable filters */

extern const  INTERPOLATIONFILTERTYPE vp9_switchable_interp
                  [VP9_SWITCHABLE_FILTERS];

extern const  int vp9_switchable_interp_map[SWITCHABLE + 1];

extern const  int vp9_is_interpolating_filter[SWITCHABLE + 1];

extern const  vp9_tree_index vp9_switchable_interp_tree
                  [2 * (VP9_SWITCHABLE_FILTERS - 1)];

extern struct vp9_token vp9_switchable_interp_encodings[VP9_SWITCHABLE_FILTERS];

extern const  vp9_prob vp9_switchable_interp_prob[VP9_SWITCHABLE_FILTERS + 1]
                                                 [VP9_SWITCHABLE_FILTERS - 1];

extern const vp9_prob vp9_default_tx_probs_32x32p[TX_SIZE_CONTEXTS]
                                                 [TX_SIZE_MAX_SB - 1];
extern const vp9_prob vp9_default_tx_probs_16x16p[TX_SIZE_CONTEXTS]
                                                 [TX_SIZE_MAX_SB - 2];
extern const vp9_prob vp9_default_tx_probs_8x8p[TX_SIZE_CONTEXTS]
                                               [TX_SIZE_MAX_SB - 3];

extern void tx_counts_to_branch_counts_32x32(unsigned int *tx_count_32x32p,
                                             unsigned int (*ct_32x32p)[2]);
extern void tx_counts_to_branch_counts_16x16(unsigned int *tx_count_16x16p,
                                             unsigned int (*ct_16x16p)[2]);
extern void tx_counts_to_branch_counts_8x8(unsigned int *tx_count_8x8p,
                                           unsigned int (*ct_8x8p)[2]);
#endif  // VP9_COMMON_VP9_ENTROPYMODE_H_
