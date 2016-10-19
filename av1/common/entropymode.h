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

#ifndef AV1_COMMON_ENTROPYMODE_H_
#define AV1_COMMON_ENTROPYMODE_H_

#include "av1/common/entropy.h"
#include "av1/common/entropymv.h"
#include "av1/common/filter.h"
#include "av1/common/seg_common.h"
#include "aom_dsp/aom_filter.h"

#ifdef __cplusplus
extern "C" {
#endif

#define BLOCK_SIZE_GROUPS 4

#define TX_SIZE_CONTEXTS 2

#define INTER_OFFSET(mode) ((mode)-NEARESTMV)
#if CONFIG_EXT_INTER
#define INTER_COMPOUND_OFFSET(mode) ((mode)-NEAREST_NEARESTMV)
#endif  // CONFIG_EXT_INTER

#if CONFIG_PALETTE
#define PALETTE_COLOR_CONTEXTS 16
#define PALETTE_MAX_SIZE 8
#define PALETTE_BLOCK_SIZES (BLOCK_LARGEST - BLOCK_8X8 + 1)
#define PALETTE_Y_MODE_CONTEXTS 3
#define PALETTE_MAX_BLOCK_SIZE (64 * 64)
#endif  // CONFIG_PALETTE

struct AV1Common;

typedef struct {
  const int16_t *scan;
  const int16_t *iscan;
  const int16_t *neighbors;
} SCAN_ORDER;

struct seg_counts {
  unsigned int tree_total[MAX_SEGMENTS];
  unsigned int tree_mispred[MAX_SEGMENTS];
  unsigned int pred[PREDICTION_PROBS][2];
};

typedef struct frame_contexts {
  aom_prob y_mode_prob[BLOCK_SIZE_GROUPS][INTRA_MODES - 1];
  aom_prob uv_mode_prob[INTRA_MODES][INTRA_MODES - 1];
#if CONFIG_EXT_PARTITION_TYPES
  aom_prob partition_prob[PARTITION_CONTEXTS][EXT_PARTITION_TYPES - 1];
#else
  aom_prob partition_prob[PARTITION_CONTEXTS][PARTITION_TYPES - 1];
#endif
  av1_coeff_probs_model coef_probs[TX_SIZES][PLANE_TYPES];
#if CONFIG_EC_MULTISYMBOL
  coeff_cdf_model coef_cdfs[TX_SIZES][PLANE_TYPES];
#endif  // CONFIG_EC_MULTISYMBOL
  aom_prob switchable_interp_prob[SWITCHABLE_FILTER_CONTEXTS]
                                 [SWITCHABLE_FILTERS - 1];

#if CONFIG_ADAPT_SCAN
  // TODO(angiebird): try aom_prob
  uint32_t non_zero_prob_4X4[TX_TYPES][16];
  uint32_t non_zero_prob_8X8[TX_TYPES][64];
  uint32_t non_zero_prob_16X16[TX_TYPES][256];
  uint32_t non_zero_prob_32X32[TX_TYPES][1024];

  DECLARE_ALIGNED(16, int16_t, scan_4X4[TX_TYPES][16]);
  DECLARE_ALIGNED(16, int16_t, scan_8X8[TX_TYPES][64]);
  DECLARE_ALIGNED(16, int16_t, scan_16X16[TX_TYPES][256]);
  DECLARE_ALIGNED(16, int16_t, scan_32X32[TX_TYPES][1024]);

  DECLARE_ALIGNED(16, int16_t, iscan_4X4[TX_TYPES][16]);
  DECLARE_ALIGNED(16, int16_t, iscan_8X8[TX_TYPES][64]);
  DECLARE_ALIGNED(16, int16_t, iscan_16X16[TX_TYPES][256]);
  DECLARE_ALIGNED(16, int16_t, iscan_32X32[TX_TYPES][1024]);

  int16_t nb_4X4[TX_TYPES][(16 + 1) * 2];
  int16_t nb_8X8[TX_TYPES][(64 + 1) * 2];
  int16_t nb_16X16[TX_TYPES][(256 + 1) * 2];
  int16_t nb_32X32[TX_TYPES][(1024 + 1) * 2];

  SCAN_ORDER sc[TX_SIZES][TX_TYPES];
#endif  // CONFIG_ADAPT_SCAN

#if CONFIG_REF_MV
  aom_prob newmv_prob[NEWMV_MODE_CONTEXTS];
  aom_prob zeromv_prob[ZEROMV_MODE_CONTEXTS];
  aom_prob refmv_prob[REFMV_MODE_CONTEXTS];
  aom_prob drl_prob[DRL_MODE_CONTEXTS];

#if CONFIG_EXT_INTER
  aom_prob new2mv_prob;
#endif  // CONFIG_EXT_INTER
#endif  // CONFIG_REF_MV

  aom_prob inter_mode_probs[INTER_MODE_CONTEXTS][INTER_MODES - 1];
#if CONFIG_EXT_INTER
  aom_prob inter_compound_mode_probs[INTER_MODE_CONTEXTS]
                                    [INTER_COMPOUND_MODES - 1];
  aom_prob interintra_prob[BLOCK_SIZE_GROUPS];
  aom_prob interintra_mode_prob[BLOCK_SIZE_GROUPS][INTERINTRA_MODES - 1];
  aom_prob wedge_interintra_prob[BLOCK_SIZES];
  aom_prob wedge_interinter_prob[BLOCK_SIZES];
#endif  // CONFIG_EXT_INTER
#if CONFIG_MOTION_VAR || CONFIG_WARPED_MOTION
  aom_prob motion_mode_prob[BLOCK_SIZES][MOTION_MODES - 1];
#endif  // CONFIG_MOTION_VAR || CONFIG_WARPED_MOTION
  aom_prob intra_inter_prob[INTRA_INTER_CONTEXTS];
  aom_prob comp_inter_prob[COMP_INTER_CONTEXTS];
  aom_prob single_ref_prob[REF_CONTEXTS][SINGLE_REFS - 1];
#if CONFIG_EXT_REFS
  aom_prob comp_ref_prob[REF_CONTEXTS][FWD_REFS - 1];
  aom_prob comp_bwdref_prob[REF_CONTEXTS][BWD_REFS - 1];
#else
  aom_prob comp_ref_prob[REF_CONTEXTS][COMP_REFS - 1];
#endif  // CONFIG_EXT_REFS
  aom_prob tx_size_probs[MAX_TX_DEPTH][TX_SIZE_CONTEXTS][MAX_TX_DEPTH];
#if CONFIG_VAR_TX
  aom_prob txfm_partition_prob[TXFM_PARTITION_CONTEXTS];
#if CONFIG_EXT_TX && CONFIG_RECT_TX
  // TODO(yuec) make this flag harmonize with the original syntax
  aom_prob rect_tx_prob[TX_SIZES - 1];
#endif  // CONFIG_EXT_TX && CONFIG_RECT_TX
#endif
  aom_prob skip_probs[SKIP_CONTEXTS];
#if CONFIG_REF_MV
  nmv_context nmvc[NMV_CONTEXTS];
#else
  nmv_context nmvc;
#endif
  int initialized;
#if CONFIG_EXT_TX
  aom_prob inter_ext_tx_prob[EXT_TX_SETS_INTER][EXT_TX_SIZES][TX_TYPES - 1];
  aom_prob intra_ext_tx_prob[EXT_TX_SETS_INTRA][EXT_TX_SIZES][INTRA_MODES]
                            [TX_TYPES - 1];
#else
  aom_prob intra_ext_tx_prob[EXT_TX_SIZES][TX_TYPES][TX_TYPES - 1];
  aom_prob inter_ext_tx_prob[EXT_TX_SIZES][TX_TYPES - 1];
#endif  // CONFIG_EXT_TX
#if CONFIG_SUPERTX
  aom_prob supertx_prob[PARTITION_SUPERTX_CONTEXTS][TX_SIZES];
#endif  // CONFIG_SUPERTX
  struct segmentation_probs seg;
#if CONFIG_EXT_INTRA
  aom_prob intra_filter_probs[INTRA_FILTERS + 1][INTRA_FILTERS - 1];
#endif  // CONFIG_EXT_INTRA
#if CONFIG_FILTER_INTRA
  aom_prob filter_intra_probs[PLANE_TYPES];
#endif  // CONFIG_FILTER_INTRA
#if CONFIG_GLOBAL_MOTION
  aom_prob global_motion_types_prob[GLOBAL_MOTION_TYPES - 1];
#endif  // CONFIG_GLOBAL_MOTION
#if CONFIG_LOOP_RESTORATION
  aom_prob switchable_restore_prob[RESTORE_SWITCHABLE_TYPES - 1];
#endif  // CONFIG_LOOP_RESTORATION
#if CONFIG_DAALA_EC
  aom_cdf_prob y_mode_cdf[BLOCK_SIZE_GROUPS][INTRA_MODES];
  aom_cdf_prob uv_mode_cdf[INTRA_MODES][INTRA_MODES];
  aom_cdf_prob partition_cdf[PARTITION_CONTEXTS][PARTITION_TYPES];
  aom_cdf_prob switchable_interp_cdf[SWITCHABLE_FILTER_CONTEXTS]
                                    [SWITCHABLE_FILTERS];
  aom_cdf_prob inter_mode_cdf[INTER_MODE_CONTEXTS][INTER_MODES];
  aom_cdf_prob intra_ext_tx_cdf[EXT_TX_SIZES][TX_TYPES][TX_TYPES];
  aom_cdf_prob inter_ext_tx_cdf[EXT_TX_SIZES][TX_TYPES];
#endif
#if CONFIG_DELTA_Q
  aom_prob delta_q_prob[DELTA_Q_CONTEXTS];
#endif
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
  av1_coeff_count_model coef[TX_SIZES][PLANE_TYPES];
  unsigned int eob_branch[TX_SIZES][PLANE_TYPES][REF_TYPES][COEF_BANDS]
                         [COEFF_CONTEXTS];
  unsigned int switchable_interp[SWITCHABLE_FILTER_CONTEXTS]
                                [SWITCHABLE_FILTERS];
#if CONFIG_ADAPT_SCAN
  unsigned int non_zero_count_4X4[TX_TYPES][16];
  unsigned int non_zero_count_8X8[TX_TYPES][64];
  unsigned int non_zero_count_16X16[TX_TYPES][256];
  unsigned int non_zero_count_32X32[TX_TYPES][1024];
  unsigned int txb_count[TX_SIZES][TX_TYPES];
#endif

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
#if CONFIG_MOTION_VAR || CONFIG_WARPED_MOTION
  unsigned int motion_mode[BLOCK_SIZES][MOTION_MODES];
#endif  // CONFIG_MOTION_VAR || CONFIG_WARPED_MOTION
  unsigned int intra_inter[INTRA_INTER_CONTEXTS][2];
  unsigned int comp_inter[COMP_INTER_CONTEXTS][2];
  unsigned int single_ref[REF_CONTEXTS][SINGLE_REFS - 1][2];
#if CONFIG_EXT_REFS
  unsigned int comp_ref[REF_CONTEXTS][FWD_REFS - 1][2];
  unsigned int comp_bwdref[REF_CONTEXTS][BWD_REFS - 1][2];
#else
  unsigned int comp_ref[REF_CONTEXTS][COMP_REFS - 1][2];
#endif  // CONFIG_EXT_REFS
  // TODO(any): tx_size_totals is only used by the encoder to decide whether
  // to use forward updates for the coeff probs, and as such it does not really
  // belong into this structure.
  unsigned int tx_size_totals[TX_SIZES];
  unsigned int tx_size[MAX_TX_DEPTH][TX_SIZE_CONTEXTS][TX_SIZES];
#if CONFIG_VAR_TX
  unsigned int txfm_partition[TXFM_PARTITION_CONTEXTS][2];
#if CONFIG_EXT_TX && CONFIG_RECT_TX
  unsigned int rect_tx[TX_SIZES - 1][2];
#endif  // CONFIG_EXT_TX && CONFIG_RECT_TX
#endif
  unsigned int skip[SKIP_CONTEXTS][2];
#if CONFIG_REF_MV
  nmv_context_counts mv[NMV_CONTEXTS];
#else
  nmv_context_counts mv;
#endif
#if CONFIG_DELTA_Q
  unsigned int delta_q[DELTA_Q_CONTEXTS][2];
#endif
#if CONFIG_EXT_TX
#if CONFIG_RECT_TX
  unsigned int tx_size_implied[TX_SIZES][TX_SIZES];
#endif  // CONFIG_RECT_TX
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
  unsigned int intra_filter[INTRA_FILTERS + 1][INTRA_FILTERS];
#endif  // CONFIG_EXT_INTRA
#if CONFIG_FILTER_INTRA
  unsigned int filter_intra[PLANE_TYPES][2];
#endif  // CONFIG_FILTER_INTRA
} FRAME_COUNTS;

extern const aom_prob av1_kf_y_mode_prob[INTRA_MODES][INTRA_MODES]
                                        [INTRA_MODES - 1];
#if CONFIG_DAALA_EC
extern aom_cdf_prob av1_kf_y_mode_cdf[INTRA_MODES][INTRA_MODES][INTRA_MODES];
#endif
#if CONFIG_PALETTE
extern const aom_prob av1_default_palette_y_mode_prob[PALETTE_BLOCK_SIZES]
                                                     [PALETTE_Y_MODE_CONTEXTS];
extern const aom_prob av1_default_palette_uv_mode_prob[2];
extern const aom_prob av1_default_palette_y_size_prob[PALETTE_BLOCK_SIZES]
                                                     [PALETTE_SIZES - 1];
extern const aom_prob av1_default_palette_uv_size_prob[PALETTE_BLOCK_SIZES]
                                                      [PALETTE_SIZES - 1];
extern const aom_prob av1_default_palette_y_color_prob[PALETTE_MAX_SIZE - 1]
                                                      [PALETTE_COLOR_CONTEXTS]
                                                      [PALETTE_COLORS - 1];
extern const aom_prob av1_default_palette_uv_color_prob[PALETTE_MAX_SIZE - 1]
                                                       [PALETTE_COLOR_CONTEXTS]
                                                       [PALETTE_COLORS - 1];
#endif  // CONFIG_PALETTE

extern const aom_tree_index av1_intra_mode_tree[TREE_SIZE(INTRA_MODES)];
extern const aom_tree_index av1_inter_mode_tree[TREE_SIZE(INTER_MODES)];
#if CONFIG_DAALA_EC
extern int av1_intra_mode_ind[INTRA_MODES];
extern int av1_intra_mode_inv[INTRA_MODES];
extern int av1_inter_mode_ind[INTER_MODES];
extern int av1_inter_mode_inv[INTER_MODES];
#endif
#if CONFIG_EXT_INTER
extern const aom_tree_index
    av1_interintra_mode_tree[TREE_SIZE(INTERINTRA_MODES)];
extern const aom_tree_index
    av1_inter_compound_mode_tree[TREE_SIZE(INTER_COMPOUND_MODES)];
#endif  // CONFIG_EXT_INTER
extern const aom_tree_index av1_partition_tree[TREE_SIZE(PARTITION_TYPES)];
#if CONFIG_EXT_PARTITION_TYPES
extern const aom_tree_index
    av1_ext_partition_tree[TREE_SIZE(EXT_PARTITION_TYPES)];
#endif
extern const aom_tree_index
    av1_switchable_interp_tree[TREE_SIZE(SWITCHABLE_FILTERS)];
#if CONFIG_PALETTE
extern const aom_tree_index av1_palette_size_tree[TREE_SIZE(PALETTE_SIZES)];
extern const aom_tree_index av1_palette_color_tree[PALETTE_MAX_SIZE - 1]
                                                  [TREE_SIZE(PALETTE_COLORS)];
#endif  // CONFIG_PALETTE
extern const aom_tree_index av1_tx_size_tree[MAX_TX_DEPTH][TREE_SIZE(TX_SIZES)];
#if CONFIG_EXT_INTRA
extern const aom_tree_index av1_intra_filter_tree[TREE_SIZE(INTRA_FILTERS)];
#endif  // CONFIG_EXT_INTRA
#if CONFIG_EXT_TX
extern const aom_tree_index av1_ext_tx_inter_tree[EXT_TX_SETS_INTER]
                                                 [TREE_SIZE(TX_TYPES)];
extern const aom_tree_index av1_ext_tx_intra_tree[EXT_TX_SETS_INTRA]
                                                 [TREE_SIZE(TX_TYPES)];
#else
extern const aom_tree_index av1_ext_tx_tree[TREE_SIZE(TX_TYPES)];
#endif  // CONFIG_EXT_TX
#if CONFIG_MOTION_VAR || CONFIG_WARPED_MOTION
extern const aom_tree_index av1_motion_mode_tree[TREE_SIZE(MOTION_MODES)];
#endif  // CONFIG_MOTION_VAR || CONFIG_WARPED_MOTION

#if CONFIG_LOOP_RESTORATION
#define RESTORE_NONE_BILATERAL_PROB 16
#define RESTORE_NONE_WIENER_PROB 64
extern const aom_tree_index
    av1_switchable_restore_tree[TREE_SIZE(RESTORE_SWITCHABLE_TYPES)];
#endif  // CONFIG_LOOP_RESTORATION
#if CONFIG_DAALA_EC
extern int av1_switchable_interp_ind[SWITCHABLE_FILTERS];
extern int av1_switchable_interp_inv[SWITCHABLE_FILTERS];

void av1_set_mode_cdfs(struct AV1Common *cm);
#endif

void av1_setup_past_independence(struct AV1Common *cm);

void av1_adapt_intra_frame_probs(struct AV1Common *cm);
void av1_adapt_inter_frame_probs(struct AV1Common *cm);
#if CONFIG_DAALA_EC
extern int av1_ext_tx_ind[TX_TYPES];
extern int av1_ext_tx_inv[TX_TYPES];
#endif

static INLINE int av1_ceil_log2(int n) {
  int i = 1, p = 2;
  while (p < n) {
    i++;
    p = p << 1;
  }
  return i;
}

#if CONFIG_PALETTE
int av1_get_palette_color_context(const uint8_t *color_map, int cols, int r,
                                  int c, int n, uint8_t *color_order,
                                  int *color_idx);
#endif  // CONFIG_PALETTE

#ifdef __cplusplus
}  // extern "C"
#endif

#endif  // AV1_COMMON_ENTROPYMODE_H_
