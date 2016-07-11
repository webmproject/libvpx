/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#ifndef VP10_COMMON_ENTROPYMODE_H_
#define VP10_COMMON_ENTROPYMODE_H_

#include "vp10/common/entropy.h"
#include "vp10/common/entropymv.h"
#include "vp10/common/filter.h"
#include "vp10/common/seg_common.h"
#include "vpx_dsp/vpx_filter.h"

#ifdef __cplusplus
extern "C" {
#endif

#define BLOCK_SIZE_GROUPS 4

#define TX_SIZE_CONTEXTS 2

#define INTER_OFFSET(mode) ((mode) - NEARESTMV)
#if CONFIG_EXT_INTER
#define INTER_COMPOUND_OFFSET(mode) ((mode) - NEAREST_NEARESTMV)
#endif  // CONFIG_EXT_INTER

#define PALETTE_COLOR_CONTEXTS 16
#define PALETTE_MAX_SIZE 8
#define PALETTE_BLOCK_SIZES (BLOCK_LARGEST - BLOCK_8X8 + 1)
#define PALETTE_Y_MODE_CONTEXTS 3
#define PALETTE_MAX_BLOCK_SIZE (64 * 64)

struct VP10Common;

struct seg_counts {
  unsigned int tree_total[MAX_SEGMENTS];
  unsigned int tree_mispred[MAX_SEGMENTS];
  unsigned int pred[PREDICTION_PROBS][2];
};

typedef struct frame_contexts {
  vpx_prob y_mode_prob[BLOCK_SIZE_GROUPS][INTRA_MODES - 1];
  vpx_prob uv_mode_prob[INTRA_MODES][INTRA_MODES - 1];
#if CONFIG_EXT_PARTITION_TYPES
  vpx_prob partition_prob[PARTITION_CONTEXTS][EXT_PARTITION_TYPES - 1];
#else
  vpx_prob partition_prob[PARTITION_CONTEXTS][PARTITION_TYPES - 1];
#endif
  vp10_coeff_probs_model coef_probs[TX_SIZES][PLANE_TYPES];
#if CONFIG_ANS
  coeff_cdf_model coef_cdfs[TX_SIZES][PLANE_TYPES];
#endif
  vpx_prob switchable_interp_prob[SWITCHABLE_FILTER_CONTEXTS]
                                 [SWITCHABLE_FILTERS - 1];

#if CONFIG_REF_MV
  vpx_prob newmv_prob[NEWMV_MODE_CONTEXTS];
  vpx_prob zeromv_prob[ZEROMV_MODE_CONTEXTS];
  vpx_prob refmv_prob[REFMV_MODE_CONTEXTS];
  vpx_prob drl_prob[DRL_MODE_CONTEXTS];

#if CONFIG_EXT_INTER
  vpx_prob new2mv_prob;
#endif  // CONFIG_EXT_INTER
#endif  // CONFIG_REF_MV

  vpx_prob inter_mode_probs[INTER_MODE_CONTEXTS][INTER_MODES - 1];
#if CONFIG_EXT_INTER
  vpx_prob inter_compound_mode_probs[INTER_MODE_CONTEXTS]
                                    [INTER_COMPOUND_MODES - 1];
  vpx_prob interintra_prob[BLOCK_SIZE_GROUPS];
  vpx_prob interintra_mode_prob[BLOCK_SIZE_GROUPS][INTERINTRA_MODES - 1];
  vpx_prob wedge_interintra_prob[BLOCK_SIZES];
  vpx_prob wedge_interinter_prob[BLOCK_SIZES];
#endif  // CONFIG_EXT_INTER
#if CONFIG_OBMC || CONFIG_WARPED_MOTION
  vpx_prob motvar_prob[BLOCK_SIZES][MOTION_VARIATIONS - 1];
#endif  // CONFIG_OBMC || CONFIG_WARPED_MOTION
  vpx_prob intra_inter_prob[INTRA_INTER_CONTEXTS];
  vpx_prob comp_inter_prob[COMP_INTER_CONTEXTS];
  vpx_prob single_ref_prob[REF_CONTEXTS][SINGLE_REFS-1];
#if CONFIG_EXT_REFS
  vpx_prob comp_ref_prob[REF_CONTEXTS][FWD_REFS-1];
  vpx_prob comp_bwdref_prob[REF_CONTEXTS][BWD_REFS-1];
#else
  vpx_prob comp_ref_prob[REF_CONTEXTS][COMP_REFS-1];
#endif  // CONFIG_EXT_REFS
  vpx_prob tx_size_probs[TX_SIZES - 1][TX_SIZE_CONTEXTS][TX_SIZES - 1];
#if CONFIG_VAR_TX
  vpx_prob txfm_partition_prob[TXFM_PARTITION_CONTEXTS];
#endif
  vpx_prob skip_probs[SKIP_CONTEXTS];
#if CONFIG_REF_MV
  nmv_context nmvc[NMV_CONTEXTS];
#else
  nmv_context nmvc;
#endif
  int initialized;
#if CONFIG_EXT_TX
  vpx_prob inter_ext_tx_prob[EXT_TX_SETS_INTER][EXT_TX_SIZES][TX_TYPES - 1];
  vpx_prob intra_ext_tx_prob[EXT_TX_SETS_INTRA][EXT_TX_SIZES][INTRA_MODES]
                            [TX_TYPES - 1];
#else
  vpx_prob intra_ext_tx_prob[EXT_TX_SIZES][TX_TYPES][TX_TYPES - 1];
  vpx_prob inter_ext_tx_prob[EXT_TX_SIZES][TX_TYPES - 1];
#endif  // CONFIG_EXT_TX
#if CONFIG_SUPERTX
  vpx_prob supertx_prob[PARTITION_SUPERTX_CONTEXTS][TX_SIZES];
#endif  // CONFIG_SUPERTX
  struct segmentation_probs seg;
#if CONFIG_EXT_INTRA
  vpx_prob ext_intra_probs[PLANE_TYPES];
  vpx_prob intra_filter_probs[INTRA_FILTERS + 1][INTRA_FILTERS - 1];
#endif  // CONFIG_EXT_INTRA
#if CONFIG_GLOBAL_MOTION
  vpx_prob global_motion_types_prob[GLOBAL_MOTION_TYPES - 1];
#endif  // CONFIG_GLOBAL_MOTION
} FRAME_CONTEXT;

typedef struct FRAME_COUNTS {
  // Note: This structure should only contain 'unsigned int' fields, or
  // aggregates built solely from 'unsigned int' fields/elements
  unsigned int kf_y_mode[INTRA_MODES][INTRA_MODES][INTRA_MODES];
  unsigned int y_mode[BLOCK_SIZE_GROUPS][INTRA_MODES];
  unsigned int uv_mode[INTRA_MODES][INTRA_MODES];
#if CONFIG_EXT_PARTITION_TYPES
  unsigned int partition[PARTITION_CONTEXTS][EXT_PARTITION_TYPES];
#else
  unsigned int partition[PARTITION_CONTEXTS][PARTITION_TYPES];
#endif
  vp10_coeff_count_model coef[TX_SIZES][PLANE_TYPES];
  unsigned int eob_branch[TX_SIZES][PLANE_TYPES][REF_TYPES]
                         [COEF_BANDS][COEFF_CONTEXTS];
  unsigned int switchable_interp[SWITCHABLE_FILTER_CONTEXTS]
                                [SWITCHABLE_FILTERS];
#if CONFIG_REF_MV
  unsigned int newmv_mode[NEWMV_MODE_CONTEXTS][2];
  unsigned int zeromv_mode[ZEROMV_MODE_CONTEXTS][2];
  unsigned int refmv_mode[REFMV_MODE_CONTEXTS][2];
  unsigned int drl_mode[DRL_MODE_CONTEXTS][2];
#if CONFIG_EXT_INTER
  unsigned int new2mv_mode[2];
#endif  // CONFIG_EXT_INTER
#endif

  unsigned int inter_mode[INTER_MODE_CONTEXTS][INTER_MODES];
#if CONFIG_EXT_INTER
  unsigned int inter_compound_mode[INTER_MODE_CONTEXTS][INTER_COMPOUND_MODES];
  unsigned int interintra[BLOCK_SIZE_GROUPS][2];
  unsigned int interintra_mode[BLOCK_SIZE_GROUPS][INTERINTRA_MODES];
  unsigned int wedge_interintra[BLOCK_SIZES][2];
  unsigned int wedge_interinter[BLOCK_SIZES][2];
#endif  // CONFIG_EXT_INTER
#if CONFIG_OBMC || CONFIG_WARPED_MOTION
  unsigned int motvar[BLOCK_SIZES][MOTION_VARIATIONS];
#endif  // CONFIG_OBMC || CONFIG_WARPED_MOTION
  unsigned int intra_inter[INTRA_INTER_CONTEXTS][2];
  unsigned int comp_inter[COMP_INTER_CONTEXTS][2];
  unsigned int single_ref[REF_CONTEXTS][SINGLE_REFS-1][2];
#if CONFIG_EXT_REFS
  unsigned int comp_ref[REF_CONTEXTS][FWD_REFS-1][2];
  unsigned int comp_bwdref[REF_CONTEXTS][BWD_REFS-1][2];
#else
  unsigned int comp_ref[REF_CONTEXTS][COMP_REFS-1][2];
#endif  // CONFIG_EXT_REFS
  unsigned int tx_size_totals[TX_SIZES];
  unsigned int tx_size[TX_SIZES - 1][TX_SIZE_CONTEXTS][TX_SIZES];
#if CONFIG_VAR_TX
  unsigned int txfm_partition[TXFM_PARTITION_CONTEXTS][2];
#endif
  unsigned int skip[SKIP_CONTEXTS][2];
#if CONFIG_REF_MV
  nmv_context_counts mv[NMV_CONTEXTS];
#else
  nmv_context_counts mv;
#endif
#if CONFIG_EXT_TX
  unsigned int inter_ext_tx[EXT_TX_SETS_INTER][EXT_TX_SIZES][TX_TYPES];
  unsigned int intra_ext_tx[EXT_TX_SETS_INTRA][EXT_TX_SIZES][INTRA_MODES]
                           [TX_TYPES];
#else
  unsigned int intra_ext_tx[EXT_TX_SIZES][TX_TYPES][TX_TYPES];
  unsigned int inter_ext_tx[EXT_TX_SIZES][TX_TYPES];
#endif  // CONFIG_EXT_TX
#if CONFIG_SUPERTX
  unsigned int supertx[PARTITION_SUPERTX_CONTEXTS][TX_SIZES][2];
  unsigned int supertx_size[TX_SIZES];
#endif  // CONFIG_SUPERTX
  struct seg_counts seg;
#if CONFIG_EXT_INTRA
  unsigned int ext_intra[PLANE_TYPES][2];
  unsigned int intra_filter[INTRA_FILTERS + 1][INTRA_FILTERS];
#endif  // CONFIG_EXT_INTRA
} FRAME_COUNTS;

extern const vpx_prob vp10_kf_y_mode_prob[INTRA_MODES][INTRA_MODES]
                                        [INTRA_MODES - 1];
extern const vpx_prob
vp10_default_palette_y_mode_prob[PALETTE_BLOCK_SIZES][PALETTE_Y_MODE_CONTEXTS];
extern const vpx_prob vp10_default_palette_uv_mode_prob[2];
extern const vpx_prob
vp10_default_palette_y_size_prob[PALETTE_BLOCK_SIZES][PALETTE_SIZES - 1];
extern const vpx_prob
vp10_default_palette_uv_size_prob[PALETTE_BLOCK_SIZES][PALETTE_SIZES - 1];
extern const vpx_prob vp10_default_palette_y_color_prob
[PALETTE_MAX_SIZE - 1][PALETTE_COLOR_CONTEXTS][PALETTE_COLORS - 1];
extern const vpx_prob vp10_default_palette_uv_color_prob
[PALETTE_MAX_SIZE - 1][PALETTE_COLOR_CONTEXTS][PALETTE_COLORS - 1];

extern const vpx_tree_index vp10_intra_mode_tree[TREE_SIZE(INTRA_MODES)];
extern const vpx_tree_index vp10_inter_mode_tree[TREE_SIZE(INTER_MODES)];
#if CONFIG_EXT_INTER
extern const vpx_tree_index vp10_interintra_mode_tree
                            [TREE_SIZE(INTERINTRA_MODES)];
extern const vpx_tree_index vp10_inter_compound_mode_tree
                            [TREE_SIZE(INTER_COMPOUND_MODES)];
#endif  // CONFIG_EXT_INTER
extern const vpx_tree_index vp10_partition_tree[TREE_SIZE(PARTITION_TYPES)];
#if CONFIG_EXT_PARTITION_TYPES
extern const vpx_tree_index vp10_ext_partition_tree
                                [TREE_SIZE(EXT_PARTITION_TYPES)];
#endif
extern const vpx_tree_index vp10_switchable_interp_tree
                                [TREE_SIZE(SWITCHABLE_FILTERS)];
extern const vpx_tree_index vp10_palette_size_tree[TREE_SIZE(PALETTE_SIZES)];
extern const vpx_tree_index
vp10_palette_color_tree[PALETTE_MAX_SIZE - 1][TREE_SIZE(PALETTE_COLORS)];
extern const vpx_tree_index
vp10_tx_size_tree[TX_SIZES - 1][TREE_SIZE(TX_SIZES)];
#if CONFIG_EXT_INTRA
extern const vpx_tree_index vp10_intra_filter_tree[TREE_SIZE(INTRA_FILTERS)];
#endif  // CONFIG_EXT_INTRA
#if CONFIG_EXT_TX
extern const vpx_tree_index
    vp10_ext_tx_inter_tree[EXT_TX_SETS_INTER][TREE_SIZE(TX_TYPES)];
extern const vpx_tree_index
    vp10_ext_tx_intra_tree[EXT_TX_SETS_INTRA][TREE_SIZE(TX_TYPES)];
#else
extern const vpx_tree_index
    vp10_ext_tx_tree[TREE_SIZE(TX_TYPES)];
#endif  // CONFIG_EXT_TX
#if CONFIG_OBMC || CONFIG_WARPED_MOTION
extern const vpx_tree_index vp10_motvar_tree[TREE_SIZE(MOTION_VARIATIONS)];
#endif  // CONFIG_OBMC || CONFIG_WARPED_MOTION

void vp10_setup_past_independence(struct VP10Common *cm);

void vp10_adapt_intra_frame_probs(struct VP10Common *cm);
void vp10_adapt_inter_frame_probs(struct VP10Common *cm);

static INLINE int vp10_ceil_log2(int n) {
  int i = 1, p = 2;
  while (p < n) {
    i++;
    p = p << 1;
  }
  return i;
}

int vp10_get_palette_color_context(const uint8_t *color_map, int cols,
                                   int r, int c, int n, int *color_order);

#ifdef __cplusplus
}  // extern "C"
#endif

#endif  // VP10_COMMON_ENTROPYMODE_H_
