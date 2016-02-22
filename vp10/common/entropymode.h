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
#define PALETTE_BLOCK_SIZES (BLOCK_64X64 - BLOCK_8X8 + 1)
#define PALETTE_Y_MODE_CONTEXTS 3

struct VP10Common;

struct tx_probs {
  vpx_prob p32x32[TX_SIZE_CONTEXTS][TX_SIZES - 1];
  vpx_prob p16x16[TX_SIZE_CONTEXTS][TX_SIZES - 2];
  vpx_prob p8x8[TX_SIZE_CONTEXTS][TX_SIZES - 3];
};

struct tx_counts {
  unsigned int p32x32[TX_SIZE_CONTEXTS][TX_SIZES];
  unsigned int p16x16[TX_SIZE_CONTEXTS][TX_SIZES - 1];
  unsigned int p8x8[TX_SIZE_CONTEXTS][TX_SIZES - 2];
  unsigned int tx_totals[TX_SIZES];
};

struct seg_counts {
  unsigned int tree_total[MAX_SEGMENTS];
  unsigned int tree_mispred[MAX_SEGMENTS];
  unsigned int pred[PREDICTION_PROBS][2];
};

typedef struct frame_contexts {
  vpx_prob y_mode_prob[BLOCK_SIZE_GROUPS][INTRA_MODES - 1];
  vpx_prob uv_mode_prob[INTRA_MODES][INTRA_MODES - 1];
  vpx_prob partition_prob[PARTITION_CONTEXTS][PARTITION_TYPES - 1];
  vp10_coeff_probs_model coef_probs[TX_SIZES][PLANE_TYPES];
  vpx_prob switchable_interp_prob[SWITCHABLE_FILTER_CONTEXTS]
                                 [SWITCHABLE_FILTERS - 1];

#if CONFIG_REF_MV
  vpx_prob newmv_prob[NEWMV_MODE_CONTEXTS];
  vpx_prob zeromv_prob[ZEROMV_MODE_CONTEXTS];
  vpx_prob refmv_prob[REFMV_MODE_CONTEXTS];
  vpx_prob drl_prob0[DRL_MODE_CONTEXTS];
  vpx_prob drl_prob1[DRL_MODE_CONTEXTS];

#if CONFIG_EXT_INTER
  vpx_prob new2mv_prob;
#endif  // CONFIG_EXT_INTER
#endif

  vpx_prob inter_mode_probs[INTER_MODE_CONTEXTS][INTER_MODES - 1];
#if CONFIG_EXT_INTER
  vpx_prob inter_compound_mode_probs[INTER_MODE_CONTEXTS]
                                    [INTER_COMPOUND_MODES - 1];
#endif  // CONFIG_EXT_INTER
#if CONFIG_OBMC
  vpx_prob obmc_prob[BLOCK_SIZES];
#endif  // CONFIG_OBMC
  vpx_prob intra_inter_prob[INTRA_INTER_CONTEXTS];
  vpx_prob comp_inter_prob[COMP_INTER_CONTEXTS];
  vpx_prob single_ref_prob[REF_CONTEXTS][SINGLE_REFS-1];
  vpx_prob comp_ref_prob[REF_CONTEXTS][COMP_REFS-1];
  struct tx_probs tx_probs;
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
} FRAME_CONTEXT;

typedef struct FRAME_COUNTS {
  unsigned int kf_y_mode[INTRA_MODES][INTRA_MODES][INTRA_MODES];
  unsigned int y_mode[BLOCK_SIZE_GROUPS][INTRA_MODES];
  unsigned int uv_mode[INTRA_MODES][INTRA_MODES];
  unsigned int partition[PARTITION_CONTEXTS][PARTITION_TYPES];
  vp10_coeff_count_model coef[TX_SIZES][PLANE_TYPES];
  unsigned int eob_branch[TX_SIZES][PLANE_TYPES][REF_TYPES]
                         [COEF_BANDS][COEFF_CONTEXTS];
  unsigned int switchable_interp[SWITCHABLE_FILTER_CONTEXTS]
                                [SWITCHABLE_FILTERS];
#if CONFIG_REF_MV
  unsigned int newmv_mode[NEWMV_MODE_CONTEXTS][2];
  unsigned int zeromv_mode[ZEROMV_MODE_CONTEXTS][2];
  unsigned int refmv_mode[REFMV_MODE_CONTEXTS][2];
  unsigned int drl_mode0[DRL_MODE_CONTEXTS][2];
  unsigned int drl_mode1[DRL_MODE_CONTEXTS][2];
#if CONFIG_EXT_INTER
  unsigned int new2mv_mode[2];
#endif  // CONFIG_EXT_INTER
#endif

  unsigned int inter_mode[INTER_MODE_CONTEXTS][INTER_MODES];
#if CONFIG_EXT_INTER
  unsigned int inter_compound_mode[INTER_MODE_CONTEXTS][INTER_COMPOUND_MODES];
#endif  // CONFIG_EXT_INTER
#if CONFIG_OBMC
  unsigned int obmc[BLOCK_SIZES][2];
#endif  // CONFIG_OBMC
  unsigned int intra_inter[INTRA_INTER_CONTEXTS][2];
  unsigned int comp_inter[COMP_INTER_CONTEXTS][2];
  unsigned int single_ref[REF_CONTEXTS][SINGLE_REFS-1][2];
  unsigned int comp_ref[REF_CONTEXTS][COMP_REFS-1][2];
  struct tx_counts tx;
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
extern const vpx_tree_index vp10_inter_compound_mode_tree
                            [TREE_SIZE(INTER_COMPOUND_MODES)];
#endif  // CONFIG_EXT_INTER
extern const vpx_tree_index vp10_partition_tree[TREE_SIZE(PARTITION_TYPES)];
extern const vpx_tree_index vp10_switchable_interp_tree
                                [TREE_SIZE(SWITCHABLE_FILTERS)];
extern const vpx_tree_index vp10_palette_size_tree[TREE_SIZE(PALETTE_SIZES)];
extern const vpx_tree_index
vp10_palette_color_tree[PALETTE_MAX_SIZE - 1][TREE_SIZE(PALETTE_COLORS)];
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

void vp10_setup_past_independence(struct VP10Common *cm);

void vp10_adapt_intra_frame_probs(struct VP10Common *cm);
void vp10_adapt_inter_frame_probs(struct VP10Common *cm);

void vp10_tx_counts_to_branch_counts_32x32(const unsigned int *tx_count_32x32p,
                                      unsigned int (*ct_32x32p)[2]);
void vp10_tx_counts_to_branch_counts_16x16(const unsigned int *tx_count_16x16p,
                                      unsigned int (*ct_16x16p)[2]);
void vp10_tx_counts_to_branch_counts_8x8(const unsigned int *tx_count_8x8p,
                                    unsigned int (*ct_8x8p)[2]);

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
