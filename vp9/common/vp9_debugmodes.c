/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include <stdio.h>
#include "vp9/common/vp9_onyxc_int.h"
#include "vp9/common/vp9_blockd.h"
#include "vp9/common/vp9_tile_common.h"
typedef struct {
  char *debug_array;
  int w;
  int h;
} DEBUG_MODE_STRUCT;

static void draw_rect(int r, int c, int w, int h, DEBUG_MODE_STRUCT *da) {
  int i;
  da->debug_array[r / 2 * da->w + c] = '+';
  for (i = r / 2 + 1; i < r / 2 + h / 2; i++) {
    da->debug_array[i * da->w + c] = '|';
  }
  for (i = c + 1; i < c + w; i++) {
    da->debug_array[r / 2 * da->w + i] = '-';
  }
}
static void debug_partitioning(VP9_COMMON * cm, MODE_INFO *m, int mi_row,
                               int mi_col, BLOCK_SIZE_TYPE bsize,
                               DEBUG_MODE_STRUCT *da) {
  const int mis = cm->mode_info_stride;
  int bwl, bhl;
  int bw, bh;
  int bsl = mi_width_log2(bsize), bs = (1 << bsl) / 2;
  int n;
  PARTITION_TYPE partition;
  BLOCK_SIZE_TYPE subsize;

  if (mi_row >= cm->mi_rows || mi_col >= cm->mi_cols)
    return;

  bwl = mi_width_log2(m->mbmi.sb_type);
  bhl = mi_height_log2(m->mbmi.sb_type);
  bw = 1 << bwl;
  bh = 1 << bhl;

  // parse the partition type
  if ((bwl == bsl) && (bhl == bsl))
    partition = PARTITION_NONE;
  else if ((bwl == bsl) && (bhl < bsl))
    partition = PARTITION_HORZ;
  else if ((bwl < bsl) && (bhl == bsl))
    partition = PARTITION_VERT;
  else if ((bwl < bsl) && (bhl < bsl))
    partition = PARTITION_SPLIT;
  else
    assert(0);

#if CONFIG_AB4X4
  if (bsize == BLOCK_SIZE_SB8X8 && m->mbmi.sb_type < BLOCK_SIZE_SB8X8)
  partition = PARTITION_SPLIT;
  if (bsize < BLOCK_SIZE_SB8X8)
  return;
#endif

#if CONFIG_AB4X4
  if (bsize >= BLOCK_SIZE_SB8X8) {
#else
  if (bsize > BLOCK_SIZE_SB8X8) {
#endif
  }

  subsize = get_subsize(bsize, partition);
  switch (partition) {
    case PARTITION_NONE:
      draw_rect(mi_row * 8, mi_col * 8, bw * 8, bh * 8, da);
      break;
    case PARTITION_HORZ:
      draw_rect(mi_row * 8, mi_col * 8, bw * 8, bh * 8, da);
      if ((mi_row + bh) < cm->mi_rows)
        draw_rect(8 * bs + mi_row * 8, mi_col * 8, bw * 8, bh * 8, da);
      break;
    case PARTITION_VERT:
      draw_rect(mi_row * 8, mi_col * 8, bw * 8, bh * 8, da);
      if ((mi_col + bw) < cm->mi_cols)
        draw_rect(mi_row * 8, 8 * bs + mi_col * 8, bw * 8, bh * 8, da);
      break;
    case PARTITION_SPLIT:
      for (n = 0; n < 4; n++) {
        int j = n >> 1, i = n & 0x01;
        debug_partitioning(cm, m + j * bs * mis + i * bs, mi_row + j * bs,
                           mi_col + i * bs, subsize, da);
      }
      break;
    default:
      assert(0);
  }
}
static void debug_partitionings(VP9_COMMON *c, DEBUG_MODE_STRUCT *da) {
  const int mis = c->mode_info_stride;
  MODE_INFO *m, *m_ptr = c->mi;
  int mi_row, mi_col;

  m_ptr += c->cur_tile_mi_col_start + c->cur_tile_mi_row_start * mis;

  for (mi_row = c->cur_tile_mi_row_start; mi_row < c->cur_tile_mi_row_end;
      mi_row += 8, m_ptr += 8 * mis) {
    m = m_ptr;
    for (mi_col = c->cur_tile_mi_col_start; mi_col < c->cur_tile_mi_col_end;
        mi_col += 8, m += 8) {
      debug_partitioning(c, m, mi_row, mi_col, BLOCK_SIZE_SB64X64, da);
    }
  }
}
void vp9_debug_tile_partitionings(VP9_COMMON *pc) {
  int tile_row, tile_col;
  DEBUG_MODE_STRUCT da;

  da.w = pc->width;
  da.h = pc->height / 2;
  da.debug_array = vpx_malloc(da.h * da.w);
  vpx_memset(da.debug_array, ' ', da.h * da.w);
  for (tile_row = 0; tile_row < pc->tile_rows; tile_row++) {
    vp9_get_tile_row_offsets(pc, tile_row);
    for (tile_col = 0; tile_col < pc->tile_columns; tile_col++) {
      vp9_get_tile_col_offsets(pc, tile_col);

      debug_partitionings(pc, &da);
    }
  }
  {
    FILE *f = fopen("partitionings.txt", "a");
    int i, j;
    fprintf(f, "\n\n\nFrame: %d \n", pc->current_video_frame);
    for (i = 0; i < da.h; i++) {
      for (j = 0; j < da.w; j++) {
        fprintf(f, "%c", da.debug_array[i * da.w + j]);
      }
      fprintf(f, "\n");
    }
    fclose(f);
  }
  vpx_free(da.debug_array);
}
