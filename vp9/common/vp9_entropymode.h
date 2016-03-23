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
#include "vp9/common/vp9_entropy.h"
#include "vp9/common/vp9_entropymv.h"

#ifdef __cplusplus
extern "C" {
#endif

#if CONFIG_INTRABC
#define INTRABC_PROB 192
#endif  // CONFIG_INTRABC

#define TX_SIZE_CONTEXTS 2

struct VP9Common;

struct tx_probs {
#if CONFIG_TX64X64
  vp9_prob p64x64[TX_SIZE_CONTEXTS][4];
#endif
  vp9_prob p32x32[TX_SIZE_CONTEXTS][3];
  vp9_prob p16x16[TX_SIZE_CONTEXTS][2];
  vp9_prob p8x8[TX_SIZE_CONTEXTS][1];
};

struct tx_counts {
#if CONFIG_TX64X64
  // counter for entropy coding
  // (In some case, tx_size may be not written to the bitstream)
  unsigned int p64x64[TX_SIZE_CONTEXTS][5];
#endif
  unsigned int p32x32[TX_SIZE_CONTEXTS][4];
  unsigned int p16x16[TX_SIZE_CONTEXTS][3];
  unsigned int p8x8[TX_SIZE_CONTEXTS][2];
#if CONFIG_SR_MODE
  // counter for tx_size actual usage (to determine ALLOW_16X16, etc.)
#if CONFIG_TX64X64
  unsigned int real_p64x64[TX_SIZE_CONTEXTS][5];
#endif  // CONFIG_TX64X64
  unsigned int real_p32x32[TX_SIZE_CONTEXTS][4];
  unsigned int real_p16x16[TX_SIZE_CONTEXTS][3];
  unsigned int real_p8x8[TX_SIZE_CONTEXTS][2];
#endif  // CONFIG_SR_MODE
};

typedef struct frame_contexts {
  vp9_prob y_mode_prob[BLOCK_SIZE_GROUPS][INTRA_MODES - 1];
  vp9_prob uv_mode_prob[INTRA_MODES][INTRA_MODES - 1];
#if CONFIG_EXT_PARTITION
  vp9_prob partition_prob[PARTITION_CONTEXTS][EXT_PARTITION_TYPES - 1];
#else
  vp9_prob partition_prob[PARTITION_CONTEXTS][PARTITION_TYPES - 1];
#endif
  vp9_coeff_probs_model coef_probs[TX_SIZES][PLANE_TYPES];
  vp9_prob switchable_interp_prob[SWITCHABLE_FILTER_CONTEXTS]
                                 [SWITCHABLE_FILTERS - 1];
  vp9_prob inter_mode_probs[INTER_MODE_CONTEXTS][INTER_MODES - 1];
#if CONFIG_NEW_INTER
  vp9_prob inter_compound_mode_probs[INTER_MODE_CONTEXTS]
                                    [INTER_COMPOUND_MODES - 1];
#endif  // CONFIG_NEW_INTER
  vp9_prob intra_inter_prob[INTRA_INTER_CONTEXTS];
  vp9_prob comp_inter_prob[COMP_INTER_CONTEXTS];
  vp9_prob single_ref_probs[REF_CONTEXTS][SINGLE_REFS - 1];
  vp9_prob comp_ref_probs[REF_CONTEXTS][COMP_REFS - 1];
  struct tx_probs tx_probs;
  vp9_prob skip_probs[SKIP_CONTEXTS];
  nmv_context nmvc;
#if CONFIG_SR_MODE
  vp9_prob sr_probs[SR_CONTEXTS];
#if SR_USE_MULTI_F
  vp9_prob sr_usfilter_probs[SR_USFILTER_CONTEXTS][SR_USFILTER_NUM - 1];
#endif  // SR_USE_MULTI_F
#endif  // CONFIG_SR_MODE
#if CONFIG_INTRABC
  nmv_context ndvc;
#endif  // CONFIG_INTRABC
#if CONFIG_FILTERINTRA
  vp9_prob filterintra_prob[TX_SIZES][INTRA_MODES];
#endif  // CONFIG_FILTERINTRA
#if CONFIG_EXT_TX
#if CONFIG_WAVELETS
  vp9_prob ext_tx_prob[TX_SIZES][EXT_TX_TYPES - 1];
#else
  vp9_prob ext_tx_prob[3][EXT_TX_TYPES - 1];
#endif  // CONFIG_WAVELETS
#endif  // CONFIG_EXT_TX
#if CONFIG_PALETTE
  vp9_prob palette_enabled_prob[10][3];
  vp9_prob palette_uv_enabled_prob[2];
  vp9_prob palette_size_prob[10][PALETTE_SIZES - 1];
  vp9_prob palette_uv_size_prob[10][PALETTE_SIZES - 1];
  vp9_prob palette_color_prob[PALETTE_MAX_SIZE - 1][PALETTE_COLOR_CONTEXTS]
                                                    [PALETTE_COLORS - 1];
  vp9_prob palette_uv_color_prob[PALETTE_MAX_SIZE - 1][PALETTE_COLOR_CONTEXTS]
                                                       [PALETTE_COLORS - 1];
#endif  // CONFIG_PALETTE
#if CONFIG_SUPERTX
  vp9_prob supertx_prob[PARTITION_SUPERTX_CONTEXTS][TX_SIZES];
#endif  // CONFIG_SUPERTX
#if CONFIG_TX_SKIP
  vp9_prob y_tx_skip_prob[2];
  vp9_prob uv_tx_skip_prob[2];
  vp9_coeff_probs_pxd coef_probs_pxd[TX_SIZES][PLANE_TYPES];
#endif  // CONFIG_TX_SKIP
#if CONFIG_COPY_MODE
  vp9_prob copy_noref_prob[COPY_MODE_CONTEXTS][BLOCK_SIZES];
  vp9_prob copy_mode_probs_l2[COPY_MODE_CONTEXTS][1];
  vp9_prob copy_mode_probs[COPY_MODE_CONTEXTS][COPY_MODE_COUNT - 2];
#endif  // CONFIG_COPY_MODE
#if CONFIG_INTERINTRA
  vp9_prob interintra_prob[BLOCK_SIZES];
#if CONFIG_WEDGE_PARTITION
  vp9_prob wedge_interintra_prob[BLOCK_SIZES];
#endif  // CONFIG_WEDGE_PARTITION
#endif  // CONFIG_INTERINTRA
#if CONFIG_WEDGE_PARTITION
  vp9_prob wedge_interinter_prob[BLOCK_SIZES];
#endif  // CONFIG_WEDGE_PARTITION
#if CONFIG_GLOBAL_MOTION
  vp9_prob global_motion_types_prob[GLOBAL_MOTION_TYPES - 1];
#endif  // CONFIG_GLOBAL_MOTION
#if CONFIG_NEW_QUANT && QUANT_PROFILES > 1 && !Q_CTX_BASED_PROFILES
  vp9_prob dq_profile_prob[QUANT_PROFILES - 1];
#endif  // CONFIG_NEW_QUANT && QUANT_PROFILES > 1 && !Q_CTX_BASED_PROFILES
} FRAME_CONTEXT;

typedef struct {
  unsigned int y_mode[BLOCK_SIZE_GROUPS][INTRA_MODES];
  unsigned int uv_mode[INTRA_MODES][INTRA_MODES];
#if CONFIG_EXT_PARTITION
  unsigned int partition[PARTITION_CONTEXTS][EXT_PARTITION_TYPES];
#else
  unsigned int partition[PARTITION_CONTEXTS][PARTITION_TYPES];
#endif
  vp9_coeff_count_model coef[TX_SIZES][PLANE_TYPES];
  unsigned int eob_branch[TX_SIZES][PLANE_TYPES][REF_TYPES]
                         [COEF_BANDS][COEFF_CONTEXTS];
  unsigned int switchable_interp[SWITCHABLE_FILTER_CONTEXTS]
                                [SWITCHABLE_FILTERS];
  unsigned int inter_mode[INTER_MODE_CONTEXTS][INTER_MODES];
#if CONFIG_NEW_INTER
  unsigned int inter_compound_mode[INTER_MODE_CONTEXTS][INTER_COMPOUND_MODES];
#endif  // CONFIG_NEW_INTER
  unsigned int intra_inter[INTRA_INTER_CONTEXTS][2];
  unsigned int comp_inter[COMP_INTER_CONTEXTS][2];
  unsigned int single_ref[REF_CONTEXTS][SINGLE_REFS-1][2];
  unsigned int comp_ref[REF_CONTEXTS][COMP_REFS-1][2];
  struct tx_counts tx;
  unsigned int skip[SKIP_CONTEXTS][2];
  nmv_context_counts mv;
#if CONFIG_SR_MODE
  unsigned int sr[SR_CONTEXTS][2];
#if SR_USE_MULTI_F
  unsigned int sr_usfilters[SR_USFILTER_CONTEXTS][SR_USFILTER_NUM];
#endif  // SR_USE_MULTI_F
#endif  // CONFIG_SR_MODE
#if CONFIG_INTRABC
  nmv_context_counts dv;
#endif  // CONFIG_INTRABC
#if CONFIG_FILTERINTRA
  unsigned int filterintra[TX_SIZES][INTRA_MODES][2];
#endif  // CONFIG_FILTERINTRA
#if CONFIG_EXT_TX
#if CONFIG_WAVELETS
  unsigned int ext_tx[TX_SIZES][EXT_TX_TYPES];
#else
  unsigned int ext_tx[3][EXT_TX_TYPES];
#endif  // CONFIG_WAVELETS
#endif  // CONFIG_EXT_TX
#if CONFIG_SUPERTX
  unsigned int supertx[PARTITION_SUPERTX_CONTEXTS][TX_SIZES][2];
  unsigned int supertx_size[BLOCK_SIZES];
#endif  // CONFIG_SUPERTX
#if CONFIG_TX_SKIP
  unsigned int y_tx_skip[2][2];
  unsigned int uv_tx_skip[2][2];
  vp9_coeff_counts_pxd coef_pxd[TX_SIZES][PLANE_TYPES];
  unsigned int eob_branch_pxd[TX_SIZES][PLANE_TYPES][REF_TYPES][COEFF_CONTEXTS];
#endif  // CONFIG_TX_SKIP
#if CONFIG_COPY_MODE
  unsigned int copy_noref[COPY_MODE_CONTEXTS][BLOCK_SIZES][2];
  unsigned int copy_mode_l2[COPY_MODE_CONTEXTS][2];
  unsigned int copy_mode[COPY_MODE_CONTEXTS][COPY_MODE_COUNT - 1];
#endif  // CONFIG_COPY_MODE
#if CONFIG_INTERINTRA
  unsigned int interintra[BLOCK_SIZES][2];
#if CONFIG_WEDGE_PARTITION
  unsigned int wedge_interintra[BLOCK_SIZES][2];
#endif  // CONFIG_WEDGE_PARTITION
#endif  // CONFIG_INTERINTRA
#if CONFIG_WEDGE_PARTITION
  unsigned int wedge_interinter[BLOCK_SIZES][2];
#endif  // CONFIG_WEDGE_PARTITION
#if CONFIG_PALETTE
  unsigned int y_palette_enabled[10][3][2];
  unsigned int uv_palette_enabled[2][2];
  unsigned int y_palette_size[10][PALETTE_SIZES];
  unsigned int uv_palette_size[10][PALETTE_SIZES];
#endif  // CONFIG_PALETTE
#if CONFIG_GLOBAL_MOTION
  unsigned int global_motion_types[GLOBAL_MOTION_TYPES];
#endif  // CONFIG_GLOBAL_MOTION
#if CONFIG_NEW_QUANT && QUANT_PROFILES > 1 && !Q_CTX_BASED_PROFILES
  unsigned int dq_profile[QUANT_PROFILES];
#endif  // CONFIG_NEW_QUANT && QUANT_PROFILES > 1 && !Q_CTX_BASED_PROFILES
} FRAME_COUNTS;

extern const vp9_prob vp9_kf_uv_mode_prob[INTRA_MODES][INTRA_MODES - 1];
extern const vp9_prob vp9_kf_y_mode_prob[INTRA_MODES][INTRA_MODES]
                                        [INTRA_MODES - 1];
#if CONFIG_EXT_PARTITION
extern const vp9_prob vp9_kf_partition_probs[PARTITION_CONTEXTS]
                                            [EXT_PARTITION_TYPES - 1];
#else
extern const vp9_prob vp9_kf_partition_probs[PARTITION_CONTEXTS]
                                            [PARTITION_TYPES - 1];
#endif
extern const vp9_tree_index vp9_intra_mode_tree[TREE_SIZE(INTRA_MODES)];
extern const vp9_tree_index vp9_inter_mode_tree[TREE_SIZE(INTER_MODES)];
extern const vp9_tree_index vp9_partition_tree[TREE_SIZE(PARTITION_TYPES)];
#if CONFIG_SR_MODE && SR_USE_MULTI_F
extern const vp9_tree_index vp9_sr_usfilter_tree[TREE_SIZE(SR_USFILTER_NUM)];
#endif  // CONFIG_SR_MODE && SR_USE_MULTI_F
#if CONFIG_EXT_PARTITION
extern const vp9_tree_index vp9_ext_partition_tree
                                [TREE_SIZE(EXT_PARTITION_TYPES)];
#endif
extern const vp9_tree_index vp9_switchable_interp_tree
                                [TREE_SIZE(SWITCHABLE_FILTERS)];
#if CONFIG_EXT_TX
extern const vp9_tree_index vp9_ext_tx_tree[TREE_SIZE(EXT_TX_TYPES)];
#if CONFIG_WAVELETS
extern const
    vp9_tree_index vp9_ext_tx_large_tree[TREE_SIZE(EXT_TX_TYPES_LARGE)];
#endif  // CONFIG_WAVELETS
#endif  // CONFIG_EXT_TX
#if CONFIG_PALETTE
extern const vp9_tree_index vp9_palette_size_tree[TREE_SIZE(PALETTE_SIZES)];
extern const vp9_tree_index vp9_palette_color_tree[TREE_SIZE(PALETTE_COLORS)];
#endif  // CONFIG_PALETTE
#if CONFIG_COPY_MODE
extern const vp9_tree_index vp9_copy_mode_tree_l2[TREE_SIZE(2)];
extern const vp9_tree_index vp9_copy_mode_tree[TREE_SIZE(COPY_MODE_COUNT - 1)];
#endif  // CONFIG_COPY_MODE

#if CONFIG_NEW_INTER
extern const vp9_tree_index vp9_inter_compound_mode_tree
                                [TREE_SIZE(INTER_COMPOUND_MODES)];
#endif  // CONFIG_NEW_INTER

#if CONFIG_NEW_QUANT && QUANT_PROFILES > 1 && !Q_CTX_BASED_PROFILES
extern const vp9_tree_index vp9_dq_profile_tree[TREE_SIZE(QUANT_PROFILES)];
#endif  // CONFIG_NEW_QUANT && QUANT_PROFILES > 1 && !Q_CTX_BASED_PROFILES

void vp9_setup_past_independence(struct VP9Common *cm);
#if CONFIG_ROW_TILE
void vp9_dec_setup_past_independence(struct VP9Common *cm,
                                     int dec_tile_row, int dec_tile_col);
#endif

void vp9_init_mode_probs(FRAME_CONTEXT *fc);

void vp9_adapt_mode_probs(struct VP9Common *cm);

#if CONFIG_TX64X64
void tx_counts_to_branch_counts_64x64(const unsigned int *tx_count_64x64p,
                                      unsigned int (*ct_64x64p)[2]);
#endif
void tx_counts_to_branch_counts_32x32(const unsigned int *tx_count_32x32p,
                                      unsigned int (*ct_32x32p)[2]);
void tx_counts_to_branch_counts_16x16(const unsigned int *tx_count_16x16p,
                                      unsigned int (*ct_16x16p)[2]);
void tx_counts_to_branch_counts_8x8(const unsigned int *tx_count_8x8p,
                                    unsigned int (*ct_8x8p)[2]);

static INLINE const vp9_prob *get_y_mode_probs(const MODE_INFO *mi,
                                               const MODE_INFO *above_mi,
                                               const MODE_INFO *left_mi,
                                               int block) {
  PREDICTION_MODE above = vp9_above_block_mode(mi, above_mi, block);
  PREDICTION_MODE left = vp9_left_block_mode(mi, left_mi, block);
#if CONFIG_INTRABC
  if (is_intrabc_mode(above)) above = DC_PRED;
  if (is_intrabc_mode(left)) left = DC_PRED;
#endif
  return vp9_kf_y_mode_prob[above][left];
}

#ifdef __cplusplus
}  // extern "C"
#endif

#endif  // VP9_COMMON_VP9_ENTROPYMODE_H_
