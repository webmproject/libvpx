/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include "vp10/common/tile_common.h"
#include "vp10/common/onyxc_int.h"
#include "vpx_dsp/vpx_dsp_common.h"

void vp10_tile_set_row(TileInfo *tile, const VP10_COMMON *cm, int row) {
  tile->mi_row_start = row * cm->tile_height;
  tile->mi_row_end   = VPXMIN(tile->mi_row_start + cm->tile_height,
                              cm->mi_rows);
}

void vp10_tile_set_col(TileInfo *tile, const VP10_COMMON *cm, int col) {
  tile->mi_col_start = col * cm->tile_width;
  tile->mi_col_end   = VPXMIN(tile->mi_col_start + cm->tile_width,
                              cm->mi_cols);
}

void vp10_tile_init(TileInfo *tile, const VP10_COMMON *cm, int row, int col) {
  vp10_tile_set_row(tile, cm, row);
  vp10_tile_set_col(tile, cm, col);
}

#if !CONFIG_EXT_TILE

# if CONFIG_EXT_PARTITION
#   define MIN_TILE_WIDTH_MAX_SB 2
#   define MAX_TILE_WIDTH_MAX_SB 32
# else
#   define MIN_TILE_WIDTH_MAX_SB 4
#   define MAX_TILE_WIDTH_MAX_SB 64
# endif  // CONFIG_EXT_PARTITION

static int get_min_log2_tile_cols(const int max_sb_cols) {
  int min_log2 = 0;
  while ((MAX_TILE_WIDTH_MAX_SB << min_log2) < max_sb_cols)
    ++min_log2;
  return min_log2;
}

static int get_max_log2_tile_cols(const int max_sb_cols) {
  int max_log2 = 1;
  while ((max_sb_cols >> max_log2) >= MIN_TILE_WIDTH_MAX_SB)
    ++max_log2;
  return max_log2 - 1;
}

void vp10_get_tile_n_bits(const int mi_cols,
                          int *min_log2_tile_cols, int *max_log2_tile_cols) {
  const int max_sb_cols =
      ALIGN_POWER_OF_TWO(mi_cols, MAX_MIB_SIZE_LOG2) >> MAX_MIB_SIZE_LOG2;
  *min_log2_tile_cols = get_min_log2_tile_cols(max_sb_cols);
  *max_log2_tile_cols = get_max_log2_tile_cols(max_sb_cols);
  assert(*min_log2_tile_cols <= *max_log2_tile_cols);
}
#endif  // !CONFIG_EXT_TILE
