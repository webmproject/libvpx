/*
 *  Copyright (c) 2014 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#ifndef VP10_ENCODER_CONTEXT_TREE_H_
#define VP10_ENCODER_CONTEXT_TREE_H_

#include "vp10/common/blockd.h"
#include "vp10/encoder/block.h"

#ifdef __cplusplus
extern "C" {
#endif

struct VP10_COMP;
struct VP10Common;
struct ThreadData;

// Structure to hold snapshot of coding context during the mode picking process
typedef struct {
  MODE_INFO mic;
  MB_MODE_INFO_EXT mbmi_ext;
  uint8_t *color_index_map[2];
#if CONFIG_VAR_TX
  uint8_t *blk_skip[MAX_MB_PLANE];
#endif
  tran_low_t *coeff[MAX_MB_PLANE][3];
  tran_low_t *qcoeff[MAX_MB_PLANE][3];
  tran_low_t *dqcoeff[MAX_MB_PLANE][3];
  uint16_t *eobs[MAX_MB_PLANE][3];

  // dual buffer pointers, 0: in use, 1: best in store
  tran_low_t *coeff_pbuf[MAX_MB_PLANE][3];
  tran_low_t *qcoeff_pbuf[MAX_MB_PLANE][3];
  tran_low_t *dqcoeff_pbuf[MAX_MB_PLANE][3];
  uint16_t *eobs_pbuf[MAX_MB_PLANE][3];

  int is_coded;
  int num_4x4_blk;
  int skip;
  int pred_pixel_ready;
  // For current partition, only if all Y, U, and V transform blocks'
  // coefficients are quantized to 0, skippable is set to 0.
  int skippable;
  int best_mode_index;
  int hybrid_pred_diff;
  int comp_pred_diff;
  int single_pred_diff;

  // TODO(jingning) Use RD_COST struct here instead. This involves a boarder
  // scope of refactoring.
  int rate;
  int64_t dist;

  // motion vector cache for adaptive motion search control in partition
  // search loop
  MV pred_mv[TOTAL_REFS_PER_FRAME];
  INTERP_FILTER pred_interp_filter;
#if CONFIG_EXT_PARTITION_TYPES
  PARTITION_TYPE partition;
#endif
} PICK_MODE_CONTEXT;

typedef struct PC_TREE {
  int index;
  PARTITION_TYPE partitioning;
  BLOCK_SIZE block_size;
  PICK_MODE_CONTEXT none;
  PICK_MODE_CONTEXT horizontal[2];
  PICK_MODE_CONTEXT vertical[2];
#if CONFIG_EXT_PARTITION_TYPES
  PICK_MODE_CONTEXT horizontala[3];
  PICK_MODE_CONTEXT horizontalb[3];
  PICK_MODE_CONTEXT verticala[3];
  PICK_MODE_CONTEXT verticalb[3];
#endif
  union {
    struct PC_TREE *split[4];
    PICK_MODE_CONTEXT *leaf_split[4];
  };
#ifdef CONFIG_SUPERTX
  PICK_MODE_CONTEXT horizontal_supertx;
  PICK_MODE_CONTEXT vertical_supertx;
  PICK_MODE_CONTEXT split_supertx;
#if CONFIG_EXT_PARTITION_TYPES
  PICK_MODE_CONTEXT horizontala_supertx;
  PICK_MODE_CONTEXT horizontalb_supertx;
  PICK_MODE_CONTEXT verticala_supertx;
  PICK_MODE_CONTEXT verticalb_supertx;
#endif
#endif
} PC_TREE;

void vp10_setup_pc_tree(struct VP10Common *cm, struct ThreadData *td);
void vp10_free_pc_tree(struct ThreadData *td);

#ifdef __cplusplus
}  // extern "C"
#endif

#endif /* VP10_ENCODER_CONTEXT_TREE_H_ */
