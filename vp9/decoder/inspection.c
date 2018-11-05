/*
 * Copyright (c) 2017, Alliance for Open Media. All rights reserved
 *
 * This source code is subject to the terms of the BSD 2 Clause License and
 * the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
 * was not distributed with this source code in the LICENSE file, you can
 * obtain it at www.aomedia.org/license/software. If the Alliance for Open
 * Media Patent License 1.0 was not distributed with this source code in the
 * PATENTS file, you can obtain it at www.aomedia.org/license/patent.
 */
#include "vp9/decoder/vp9_decoder.h"
#include "vp9/decoder/inspection.h"
#include "vp9/common/vp9_enums.h"

static void ifd_init_mi_rc(insp_frame_data *fd, int mi_cols, int mi_rows) {
  fd->mi_cols = mi_cols;
  fd->mi_rows = mi_rows;
  fd->mi_grid = (insp_mi_data *)vpx_malloc(sizeof(insp_mi_data) * fd->mi_rows *
                                           fd->mi_cols);
}

void ifd_init(insp_frame_data *fd, int frame_width, int frame_height) {
  int mi_cols = ALIGN_POWER_OF_TWO(frame_width, 3) >> MI_SIZE_LOG2;
  int mi_rows = ALIGN_POWER_OF_TWO(frame_height, 3) >> MI_SIZE_LOG2;
  ifd_init_mi_rc(fd, mi_cols, mi_rows);
}

void ifd_clear(insp_frame_data *fd) {
  vpx_free(fd->mi_grid);
  fd->mi_grid = NULL;
}

/* TODO(negge): This function may be called by more than one thread when using
               a multi-threaded decoder and this may cause a data race. */
int ifd_inspect(insp_frame_data *fd, void *decoder) {
  struct VP9Decoder *pbi = (struct VP9Decoder *)decoder;
  VP9_COMMON *const cm = &pbi->common;
  if (fd->mi_rows != cm->mi_rows || fd->mi_cols != cm->mi_cols) {
    ifd_clear(fd);
    ifd_init_mi_rc(fd, cm->mi_rows, cm->mi_cols);
  }
  fd->show_frame = cm->show_frame;
  fd->frame_type = cm->frame_type;
  fd->base_qindex = cm->base_qindex;
  // TODO(jimbankoski): copy tile data
  fd->show_existing_frame = cm->show_existing_frame;

  // fd->tile_mi_cols = cm->tile_width;
  // fd->tile_mi_rows = cm->tile_height;
#if CONFIG_ACCOUNTING
  fd->accounting = &pbi->accounting;
#endif
  int i, j;
  for (i = 0; i < MAX_SEGMENTS; i++) {
    for (j = 0; j < 2; j++) {
      fd->y_dequant[i][j] = cm->y_dequant[i][j];
      fd->uv_dequant[i][j] = cm->uv_dequant[i][j];
    }
  }
  for (j = 0; j < cm->mi_rows; j++) {
    for (i = 0; i < cm->mi_cols; i++) {
      const MODE_INFO *bmi = cm->mi_grid_visible[j * cm->mi_stride + i];
      insp_mi_data *mi = &fd->mi_grid[j * cm->mi_cols + i];
      // Segment
      mi->segment_id = bmi->segment_id;
      // Motion Vectors
      mi->mv[0].row = bmi->mv[0].as_mv.row;
      mi->mv[0].col = bmi->mv[0].as_mv.col;

      if (bmi->ref_frame[1] == -1) {
        mi->mv[1].row = 0;
        mi->mv[1].col = 0;
      } else {
        mi->mv[1].row = bmi->mv[1].as_mv.row;
        mi->mv[1].col = bmi->mv[1].as_mv.col;
      }
      // Reference Frames
      mi->ref_frame[0] = bmi->ref_frame[0];
      mi->ref_frame[1] = bmi->ref_frame[1];
      // Prediction Mode
      mi->mode = bmi->mode;
      // Prediction Mode for Chromatic planes
      // TODO(jbb): mbx handle UV_MODE_INVALID case
      // if (mi->mode < INTRA_MODES) {
      mi->uv_mode = bmi->uv_mode;
      // } else {
      //   mi->uv_mode = UV_MODE_INVALID;
      // }
      // Block Size
      mi->sb_type = bmi->sb_type;
      // Skip Flag
      mi->skip = bmi->skip;
      // Filters
      mi->filter = bmi->interp_filter;
      // Transform
      // TODO(jbb): mi->tx_type = bmi->tx_type;
      mi->tx_size = bmi->tx_size;
      // delta_q
      // TODO(jbb): mi->current_qindex = bmi->current_q_index;
    }
  }
  return 1;
}
