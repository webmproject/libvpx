/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#ifndef VP9_ENCODER_VP9_RDOPT_H_
#define VP9_ENCODER_VP9_RDOPT_H_

#include "vp9/common/vp9_blockd.h"

#include "vp9/encoder/vp9_block.h"
#include "vp9/encoder/vp9_context_tree.h"

#ifdef __cplusplus
extern "C" {
#endif

struct TileInfo;
struct VP9_COMP;
struct macroblock;
struct RD_COST;

void vp9_rd_pick_intra_mode_sb(struct VP9_COMP *cpi, struct macroblock *x,
#if CONFIG_INTRABC
                               int mi_row, int mi_col,
#endif  // CONFIG_INTRABC
                               struct RD_COST *rd_cost, BLOCK_SIZE bsize,
                               PICK_MODE_CONTEXT *ctx, int64_t best_rd);

void vp9_rd_pick_inter_mode_sb(struct VP9_COMP *cpi, struct macroblock *x,
                               const struct TileInfo *const tile,
                               int mi_row, int mi_col,
                               struct RD_COST *rd_cost,
#if CONFIG_SUPERTX
                               int *returnrate_nocoef,
#endif
                               BLOCK_SIZE bsize, PICK_MODE_CONTEXT *ctx,
                               int64_t best_rd_so_far);

void vp9_rd_pick_inter_mode_sb_seg_skip(struct VP9_COMP *cpi,
                                        struct macroblock *x,
                                        struct RD_COST *rd_cost,
                                        BLOCK_SIZE bsize,
                                        PICK_MODE_CONTEXT *ctx,
                                        int64_t best_rd_so_far);

void vp9_rd_pick_inter_mode_sub8x8(struct VP9_COMP *cpi,
                                   struct macroblock *x,
                                   const struct TileInfo *const tile,
                                   int mi_row, int mi_col,
                                   struct RD_COST *rd_cost,
#if CONFIG_SUPERTX
                                   int *returnrate_nocoef,
#endif
                                   BLOCK_SIZE bsize, PICK_MODE_CONTEXT *ctx,
                                   int64_t best_rd_so_far);

#if CONFIG_SUPERTX
void txfm_rd_in_plane_supertx(MACROBLOCK *x,
                              int *rate, int64_t *distortion,
                              int *skippable, int64_t *sse,
                              int64_t ref_best_rd, int plane,
                              BLOCK_SIZE bsize, TX_SIZE tx_size,
                              int use_fast_coef_casting);
void txfm_rd_in_plane(MACROBLOCK *x,
                      int *rate, int64_t *distortion,
                      int *skippable, int64_t *sse,
                      int64_t ref_best_rd, int plane,
                      BLOCK_SIZE bsize, TX_SIZE tx_size,
                      int use_fast_coef_casting);
#endif
#ifdef __cplusplus
}  // extern "C"
#endif

#endif  // VP9_ENCODER_VP9_RDOPT_H_
