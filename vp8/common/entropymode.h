/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */


#ifndef __INC_ENTROPYMODE_H
#define __INC_ENTROPYMODE_H

#include "blockd.h"
#include "treecoder.h"

#define SUBMVREF_COUNT 5
#define VP8_NUMMBSPLITS 4

typedef const int vp9_mbsplit[16];

extern vp9_mbsplit vp9_mbsplits [VP8_NUMMBSPLITS];

extern const int vp9_mbsplit_count [VP8_NUMMBSPLITS];    /* # of subsets */

extern const vp8_prob vp9_mbsplit_probs [VP8_NUMMBSPLITS - 1];

extern int vp9_mv_cont(const int_mv *l, const int_mv *a);

extern const vp8_prob vp9_sub_mv_ref_prob [VP8_SUBMVREFS - 1];
extern const vp8_prob vp9_sub_mv_ref_prob2 [SUBMVREF_COUNT][VP8_SUBMVREFS - 1];


extern const unsigned int vp9_kf_default_bmode_counts[VP8_BINTRAMODES][VP8_BINTRAMODES][VP8_BINTRAMODES];


extern const vp8_tree_index vp9_bmode_tree[];

extern const vp8_tree_index  vp9_ymode_tree[];
extern const vp8_tree_index  vp9_kf_ymode_tree[];
extern const vp8_tree_index  vp9_uv_mode_tree[];
#define vp8_sb_ymode_tree vp9_uv_mode_tree
extern const vp8_tree_index  vp9_i8x8_mode_tree[];
extern const vp8_tree_index  vp9_mbsplit_tree[];
extern const vp8_tree_index  vp9_mv_ref_tree[];
extern const vp8_tree_index  vp9_sb_mv_ref_tree[];
extern const vp8_tree_index  vp9_sub_mv_ref_tree[];

extern struct vp8_token_struct vp9_bmode_encodings   [VP8_BINTRAMODES];
extern struct vp8_token_struct vp9_ymode_encodings   [VP8_YMODES];
extern struct vp8_token_struct vp9_sb_kf_ymode_encodings [VP8_I32X32_MODES];
extern struct vp8_token_struct vp9_kf_ymode_encodings [VP8_YMODES];
extern struct vp8_token_struct vp9_i8x8_mode_encodings  [VP8_I8X8_MODES];
extern struct vp8_token_struct vp9_uv_mode_encodings  [VP8_UV_MODES];
extern struct vp8_token_struct vp9_mbsplit_encodings  [VP8_NUMMBSPLITS];

/* Inter mode values do not start at zero */

extern struct vp8_token_struct vp9_mv_ref_encoding_array    [VP8_MVREFS];
extern struct vp8_token_struct vp9_sb_mv_ref_encoding_array    [VP8_MVREFS];
extern struct vp8_token_struct vp9_sub_mv_ref_encoding_array [VP8_SUBMVREFS];

void vp9_entropy_mode_init(void);

struct VP8Common;
void vp9_init_mbmode_probs(struct VP8Common *x);
extern void vp9_init_mode_contexts(struct VP8Common *pc);
extern void vp9_update_mode_context(struct VP8Common *pc);;
extern void vp9_accum_mv_refs(struct VP8Common *pc,
                              MB_PREDICTION_MODE m,
                              const int ct[4]);

void vp9_default_bmode_probs(vp8_prob dest [VP8_BINTRAMODES - 1]);
void vp9_kf_default_bmode_probs(vp8_prob dest [VP8_BINTRAMODES] [VP8_BINTRAMODES] [VP8_BINTRAMODES - 1]);

void vp9_adapt_mode_probs(struct VP8Common *);

#define VP8_SWITCHABLE_FILTERS 2 /* number of switchable filters */
extern const  INTERPOLATIONFILTERTYPE vp9_switchable_interp
                  [VP8_SWITCHABLE_FILTERS];
extern const  int vp9_switchable_interp_map[SWITCHABLE + 1];
extern const  vp8_tree_index vp9_switchable_interp_tree
                  [2*(VP8_SWITCHABLE_FILTERS - 1)];
extern struct vp8_token_struct vp9_switchable_interp_encodings
                  [VP8_SWITCHABLE_FILTERS];
extern const  vp8_prob vp9_switchable_interp_prob
                  [VP8_SWITCHABLE_FILTERS + 1][VP8_SWITCHABLE_FILTERS - 1];
#endif
