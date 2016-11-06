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

#ifndef AV1_ENCODER_RDOPT_H_
#define AV1_ENCODER_RDOPT_H_

#include "av1/common/blockd.h"

#include "av1/encoder/block.h"
#include "av1/encoder/context_tree.h"

#ifdef __cplusplus
extern "C" {
#endif

struct TileInfo;
struct AV1_COMP;
struct macroblock;
struct RD_COST;

#if CONFIG_VAR_TX
static INLINE void av1_init_rd_stats(RD_STATS *rd_stats) {
#if CONFIG_RD_DEBUG
  int plane;
#endif
  rd_stats->rate = 0;
  rd_stats->dist = 0;
  rd_stats->sse = 0;
  rd_stats->skip = 1;
#if CONFIG_RD_DEBUG
  for (plane = 0; plane < MAX_MB_PLANE; ++plane) {
    int r, c;
    rd_stats->txb_coeff_cost[plane] = 0;
    for (r = 0; r < TXB_COEFF_COST_MAP_SIZE; ++r)
      for (c = 0; c < TXB_COEFF_COST_MAP_SIZE; ++c)
        rd_stats->txb_coeff_cost_map[plane][r][c] = 0;
  }
#endif
}

static INLINE void av1_invalid_rd_stats(RD_STATS *rd_stats) {
#if CONFIG_RD_DEBUG
  int plane;
#endif
  rd_stats->rate = INT_MAX;
  rd_stats->dist = INT64_MAX;
  rd_stats->sse = INT64_MAX;
  rd_stats->skip = 0;
#if CONFIG_RD_DEBUG
  for (plane = 0; plane < MAX_MB_PLANE; ++plane) {
    int r, c;
    rd_stats->txb_coeff_cost[plane] = INT_MAX;
    for (r = 0; r < TXB_COEFF_COST_MAP_SIZE; ++r)
      for (c = 0; c < TXB_COEFF_COST_MAP_SIZE; ++c)
        rd_stats->txb_coeff_cost_map[plane][r][c] = INT_MAX;
  }
#endif
}

static INLINE void av1_merge_rd_stats(RD_STATS *rd_stats_dst,
                                      const RD_STATS *rd_stats_src) {
#if CONFIG_RD_DEBUG
  int plane;
#endif
  rd_stats_dst->rate += rd_stats_src->rate;
  rd_stats_dst->dist += rd_stats_src->dist;
  rd_stats_dst->sse += rd_stats_src->sse;
  rd_stats_dst->skip &= rd_stats_src->skip;
#if CONFIG_RD_DEBUG
  for (plane = 0; plane < MAX_MB_PLANE; ++plane) {
    int r, c;
    int ref_txb_coeff_cost = 0;
    rd_stats_dst->txb_coeff_cost[plane] += rd_stats_src->txb_coeff_cost[plane];
    // TODO(angiebird): optimize this part
    for (r = 0; r < TXB_COEFF_COST_MAP_SIZE; ++r)
      for (c = 0; c < TXB_COEFF_COST_MAP_SIZE; ++c) {
        rd_stats_dst->txb_coeff_cost_map[plane][r][c] +=
            rd_stats_src->txb_coeff_cost_map[plane][r][c];
        ref_txb_coeff_cost += rd_stats_dst->txb_coeff_cost_map[plane][r][c];
      }
    assert(ref_txb_coeff_cost == rd_stats_dst->txb_coeff_cost[plane]);
  }
#endif
}
#endif

int av1_cost_coeffs(const AV1_COMMON *const cm, MACROBLOCK *x, int plane,
                    int block, int coeff_ctx, TX_SIZE tx_size,
                    const int16_t *scan, const int16_t *nb,
                    int use_fast_coef_costing);
void av1_rd_pick_intra_mode_sb(const struct AV1_COMP *cpi, struct macroblock *x,
                               struct RD_COST *rd_cost, BLOCK_SIZE bsize,
                               PICK_MODE_CONTEXT *ctx, int64_t best_rd);

unsigned int av1_get_sby_perpixel_variance(const AV1_COMP *cpi,
                                           const struct buf_2d *ref,
                                           BLOCK_SIZE bs);
#if CONFIG_AOM_HIGHBITDEPTH
unsigned int av1_high_get_sby_perpixel_variance(const AV1_COMP *cpi,
                                                const struct buf_2d *ref,
                                                BLOCK_SIZE bs, int bd);
#endif

void av1_rd_pick_inter_mode_sb(const struct AV1_COMP *cpi,
                               struct TileDataEnc *tile_data,
                               struct macroblock *x, int mi_row, int mi_col,
                               struct RD_COST *rd_cost,
#if CONFIG_SUPERTX
                               int *returnrate_nocoef,
#endif  // CONFIG_SUPERTX
                               BLOCK_SIZE bsize, PICK_MODE_CONTEXT *ctx,
                               int64_t best_rd_so_far);

void av1_rd_pick_inter_mode_sb_seg_skip(
    const struct AV1_COMP *cpi, struct TileDataEnc *tile_data,
    struct macroblock *x, struct RD_COST *rd_cost, BLOCK_SIZE bsize,
    PICK_MODE_CONTEXT *ctx, int64_t best_rd_so_far);

int av1_internal_image_edge(const struct AV1_COMP *cpi);
int av1_active_h_edge(const struct AV1_COMP *cpi, int mi_row, int mi_step);
int av1_active_v_edge(const struct AV1_COMP *cpi, int mi_col, int mi_step);
int av1_active_edge_sb(const struct AV1_COMP *cpi, int mi_row, int mi_col);

void av1_rd_pick_inter_mode_sub8x8(const struct AV1_COMP *cpi,
                                   struct TileDataEnc *tile_data,
                                   struct macroblock *x, int mi_row, int mi_col,
                                   struct RD_COST *rd_cost,
#if CONFIG_SUPERTX
                                   int *returnrate_nocoef,
#endif  // CONFIG_SUPERTX
                                   BLOCK_SIZE bsize, PICK_MODE_CONTEXT *ctx,
                                   int64_t best_rd_so_far);

#if CONFIG_SUPERTX
#if CONFIG_VAR_TX
void av1_tx_block_rd_b(const AV1_COMP *cpi, MACROBLOCK *x, TX_SIZE tx_size,
                       int blk_row, int blk_col, int plane, int block,
                       int plane_bsize, int coeff_ctx, RD_STATS *rd_stats);
#endif

void av1_txfm_rd_in_plane_supertx(MACROBLOCK *x, const AV1_COMP *cpi, int *rate,
                                  int64_t *distortion, int *skippable,
                                  int64_t *sse, int64_t ref_best_rd, int plane,
                                  BLOCK_SIZE bsize, TX_SIZE tx_size,
                                  int use_fast_coef_casting);
#endif  // CONFIG_SUPERTX

#ifdef __cplusplus
}  // extern "C"
#endif

#endif  // AV1_ENCODER_RDOPT_H_
