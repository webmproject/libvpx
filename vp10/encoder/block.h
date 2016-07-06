/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#ifndef VP10_ENCODER_BLOCK_H_
#define VP10_ENCODER_BLOCK_H_

#include "vp10/common/entropymv.h"
#include "vp10/common/entropy.h"
#if CONFIG_REF_MV
#include "vp10/common/mvref_common.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
  unsigned int sse;
  int sum;
  unsigned int var;
} diff;

typedef struct macroblock_plane {
  DECLARE_ALIGNED(16, int16_t, src_diff[MAX_SB_SQUARE]);
  tran_low_t *qcoeff;
  tran_low_t *coeff;
  uint16_t *eobs;
  struct buf_2d src;

  // Quantizer setings
  const int16_t *quant_fp;
  const int16_t *round_fp;
  const int16_t *quant;
  const int16_t *quant_shift;
  const int16_t *zbin;
  const int16_t *round;
#if CONFIG_NEW_QUANT
  cuml_bins_type_nuq *cuml_bins_nuq[QUANT_PROFILES];
#endif  // CONFIG_NEW_QUANT

  int64_t quant_thred[2];
} MACROBLOCK_PLANE;

/* The [2] dimension is for whether we skip the EOB node (i.e. if previous
 * coefficient in this block was zero) or not. */
typedef unsigned int vp10_coeff_cost[PLANE_TYPES][REF_TYPES][COEF_BANDS][2]
                                   [COEFF_CONTEXTS][ENTROPY_TOKENS];

typedef struct {
  int_mv ref_mvs[MODE_CTX_REF_FRAMES][MAX_MV_REF_CANDIDATES];
  int16_t mode_context[MODE_CTX_REF_FRAMES];
#if CONFIG_REF_MV
  uint8_t ref_mv_count[MODE_CTX_REF_FRAMES];
  CANDIDATE_MV ref_mv_stack[MODE_CTX_REF_FRAMES][MAX_REF_MV_STACK_SIZE];
#if CONFIG_EXT_INTER
  int16_t compound_mode_context[MODE_CTX_REF_FRAMES];
#endif  // CONFIG_EXT_INTER
#endif
} MB_MODE_INFO_EXT;

typedef struct {
  uint8_t best_palette_color_map[MAX_SB_SQUARE];
  float kmeans_data_buf[2 * MAX_SB_SQUARE];
  uint8_t kmeans_indices_buf[MAX_SB_SQUARE];
  uint8_t kmeans_pre_indices_buf[MAX_SB_SQUARE];
} PALETTE_BUFFER;

typedef struct macroblock MACROBLOCK;
struct macroblock {
  struct macroblock_plane plane[MAX_MB_PLANE];

  MACROBLOCKD e_mbd;
  MB_MODE_INFO_EXT *mbmi_ext;
  int skip_block;
  int select_tx_size;
  int skip_optimize;
  int q_index;

  // The equivalent error at the current rdmult of one whole bit (not one
  // bitcost unit).
  int errorperbit;
  // The equivalend SAD error of one (whole) bit at the current quantizer
  // for large blocks.
  int sadperbit16;
  // The equivalend SAD error of one (whole) bit at the current quantizer
  // for sub-8x8 blocks.
  int sadperbit4;
  int rddiv;
  int rdmult;
  int mb_energy;
  int * m_search_count_ptr;
  int * ex_search_count_ptr;

  // These are set to their default values at the beginning, and then adjusted
  // further in the encoding process.
  BLOCK_SIZE min_partition_size;
  BLOCK_SIZE max_partition_size;

  int mv_best_ref_index[MAX_REF_FRAMES];
  unsigned int max_mv_context[MAX_REF_FRAMES];
  unsigned int source_variance;
  unsigned int recon_variance;
  unsigned int pred_sse[MAX_REF_FRAMES];
  int pred_mv_sad[MAX_REF_FRAMES];

#if CONFIG_REF_MV
  int *nmvjointcost;
  int nmv_vec_cost[NMV_CONTEXTS][MV_JOINTS];
  int *nmvcost[NMV_CONTEXTS][2];
  int *nmvcost_hp[NMV_CONTEXTS][2];
  int **mv_cost_stack[NMV_CONTEXTS];
  int *nmvjointsadcost;
  int zero_rmv_cost[NMV_CONTEXTS][2];
  int comp_rmv_cost[2];
#else
  int nmvjointcost[MV_JOINTS];
  int *nmvcost[2];
  int *nmvcost_hp[2];
  int nmvjointsadcost[MV_JOINTS];
#endif

  int **mvcost;
  int *nmvsadcost[2];
  int *nmvsadcost_hp[2];
  int **mvsadcost;

  PALETTE_BUFFER *palette_buffer;

  // These define limits to motion vector components to prevent them
  // from extending outside the UMV borders
  int mv_col_min;
  int mv_col_max;
  int mv_row_min;
  int mv_row_max;

  // Notes transform blocks where no coefficents are coded.
  // Set during mode selection. Read during block encoding.
  uint8_t zcoeff_blk[TX_SIZES][MAX_MIB_SIZE * MAX_MIB_SIZE * 4];
#if CONFIG_VAR_TX
  uint8_t blk_skip[MAX_MB_PLANE][MAX_MIB_SIZE * MAX_MIB_SIZE * 4];
#if CONFIG_REF_MV
  uint8_t blk_skip_drl[MAX_MB_PLANE][MAX_MIB_SIZE * MAX_MIB_SIZE * 4];
#endif
#endif

  int skip;

  int encode_breakout;

  // note that token_costs is the cost when eob node is skipped
  vp10_coeff_cost token_costs[TX_SIZES];

  int optimize;

  // indicate if it is in the rd search loop or encoding process
  int use_lp32x32fdct;

  // Used to store sub partition's choices.
  MV pred_mv[MAX_REF_FRAMES];

  // Store the best motion vector during motion search
  int_mv best_mv;

  // Strong color activity detection. Used in RTC coding mode to enhance
  // the visual quality at the boundary of moving color objects.
  uint8_t color_sensitivity[2];

  // use default transform and skip transform type search for intra modes
  int use_default_intra_tx_type;
  // use default transform and skip transform type search for inter modes
  int use_default_inter_tx_type;
};

#ifdef __cplusplus
}  // extern "C"
#endif

#endif  // VP10_ENCODER_BLOCK_H_
