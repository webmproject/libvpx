/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#ifndef VP9_ENCODER_VP9_BLOCK_H_
#define VP9_ENCODER_VP9_BLOCK_H_

#include "vp9/common/vp9_entropymv.h"
#include "vp9/common/vp9_entropy.h"
#include "vpx_ports/mem.h"
#include "vp9/common/vp9_onyxc_int.h"

#ifdef __cplusplus
extern "C" {
#endif

// motion search site
typedef struct {
  MV mv;
  int offset;
} search_site;

// Structure to hold snapshot of coding context during the mode picking process
typedef struct {
  MODE_INFO mic;
  uint8_t *zcoeff_blk;
  int16_t *coeff[MAX_MB_PLANE][3];
  int16_t *qcoeff[MAX_MB_PLANE][3];
  int16_t *dqcoeff[MAX_MB_PLANE][3];
  uint16_t *eobs[MAX_MB_PLANE][3];

  // dual buffer pointers, 0: in use, 1: best in store
  int16_t *coeff_pbuf[MAX_MB_PLANE][3];
  int16_t *qcoeff_pbuf[MAX_MB_PLANE][3];
  int16_t *dqcoeff_pbuf[MAX_MB_PLANE][3];
  uint16_t *eobs_pbuf[MAX_MB_PLANE][3];

  int is_coded;
  int num_4x4_blk;
  int skip;
  int_mv best_ref_mv[2];
  int_mv ref_mvs[MAX_REF_FRAMES][MAX_MV_REF_CANDIDATES];
  int rate;
  int distortion;
  int best_mode_index;
  int rddiv;
  int rdmult;
  int hybrid_pred_diff;
  int comp_pred_diff;
  int single_pred_diff;
  int64_t tx_rd_diff[TX_MODES];
  int64_t best_filter_diff[SWITCHABLE_FILTER_CONTEXTS];

  // motion vector cache for adaptive motion search control in partition
  // search loop
  int_mv pred_mv[MAX_REF_FRAMES];
  INTERP_FILTER pred_interp_filter;
} PICK_MODE_CONTEXT;

struct macroblock_plane {
  DECLARE_ALIGNED(16, int16_t, src_diff[64 * 64]);
  int16_t *qcoeff;
  int16_t *coeff;
  uint16_t *eobs;
  struct buf_2d src;

  // Quantizer setings
  int16_t *quant;
  int16_t *quant_shift;
  int16_t *zbin;
  int16_t *round;

  // Zbin Over Quant value
  int16_t zbin_extra;
};
typedef struct PC_TREE {
  int index;
  PARTITION_TYPE partitioning;
  BLOCK_SIZE block_size;
  PICK_MODE_CONTEXT none;
  PICK_MODE_CONTEXT horizontal[2];
  PICK_MODE_CONTEXT vertical[2];
  union {
    struct PC_TREE *split[4];
    PICK_MODE_CONTEXT *leaf_split[4];
  };
} PC_TREE;

/* The [2] dimension is for whether we skip the EOB node (i.e. if previous
 * coefficient in this block was zero) or not. */
typedef unsigned int vp9_coeff_cost[PLANE_TYPES][REF_TYPES][COEF_BANDS][2]
                                   [COEFF_CONTEXTS][ENTROPY_TOKENS];

typedef struct macroblock MACROBLOCK;
struct macroblock {
  struct macroblock_plane plane[MAX_MB_PLANE];

  MACROBLOCKD e_mbd;
  int skip_block;
  int select_txfm_size;
  int skip_recode;
  int skip_optimize;
  int q_index;

  search_site *ss;
  int ss_count;
  int searches_per_step;

  int errorperbit;
  int sadperbit16;
  int sadperbit4;
  int rddiv;
  int rdmult;
  unsigned int mb_energy;

  int mv_best_ref_index[MAX_REF_FRAMES];
  unsigned int max_mv_context[MAX_REF_FRAMES];
  unsigned int source_variance;
  unsigned int pred_sse[MAX_REF_FRAMES];
  int pred_mv_sad[MAX_REF_FRAMES];

  int nmvjointcost[MV_JOINTS];
  int nmvcosts[2][MV_VALS];
  int *nmvcost[2];
  int nmvcosts_hp[2][MV_VALS];
  int *nmvcost_hp[2];
  int **mvcost;

  int nmvjointsadcost[MV_JOINTS];
  int nmvsadcosts[2][MV_VALS];
  int *nmvsadcost[2];
  int nmvsadcosts_hp[2][MV_VALS];
  int *nmvsadcost_hp[2];
  int **mvsadcost;

  int mbmode_cost[INTRA_MODES];
  unsigned inter_mode_cost[INTER_MODE_CONTEXTS][INTER_MODES];
  int intra_uv_mode_cost[FRAME_TYPES][INTRA_MODES];
  int y_mode_costs[INTRA_MODES][INTRA_MODES][INTRA_MODES];
  int switchable_interp_costs[SWITCHABLE_FILTER_CONTEXTS][SWITCHABLE_FILTERS];

  // These define limits to motion vector components to prevent them
  // from extending outside the UMV borders
  int mv_col_min;
  int mv_col_max;
  int mv_row_min;
  int mv_row_max;

  uint8_t zcoeff_blk[TX_SIZES][256];
  int skip;

  int encode_breakout;

  int in_active_map;

  // note that token_costs is the cost when eob node is skipped
  vp9_coeff_cost token_costs[TX_SIZES];

  int in_static_area;

  int optimize;

  // indicate if it is in the rd search loop or encoding process
  int use_lp32x32fdct;
  int skip_encode;

  // Used to store sub partition's choices.
  int_mv pred_mv[MAX_REF_FRAMES];

  PICK_MODE_CONTEXT *leaf_tree;
  PC_TREE *pc_tree;
  PC_TREE *pc_root;
  int partition_cost[PARTITION_CONTEXTS][PARTITION_TYPES];

  void (*fwd_txm4x4)(const int16_t *input, int16_t *output, int stride);
};


#ifdef __cplusplus
}  // extern "C"
#endif

#endif  // VP9_ENCODER_VP9_BLOCK_H_
