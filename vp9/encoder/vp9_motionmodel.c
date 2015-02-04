/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include <limits.h>

#include "vpx_mem/vpx_mem.h"
#include "vp9/encoder/vp9_segmentation.h"
#include "vp9/encoder/vp9_mcomp.h"
#include "vp9/common/vp9_blockd.h"
#include "vp9/common/vp9_reconinter.h"
#include "vp9/common/vp9_reconintra.h"
#include "vp9/common/vp9_systemdependent.h"
#include "vp9/encoder/vp9_global_motion.h"


static unsigned int do_motion_iteration(VP9_COMP *cpi,
                                        const MV *ref_mv,
                                        MV *dst_mv,
                                        int bsize,
                                        int mb_row,
                                        int mb_col,
                                        unsigned int *sse) {
  MACROBLOCK *const x = &cpi->mb;
  MACROBLOCKD *const xd = &x->e_mbd;

  PREDICTION_MODE tmp_mode = xd->mi[0].src_mi->mbmi.mode;
  MV tmp_mv = xd->mi[0].src_mi->mbmi.mv[0].as_mv;
  int tmp_frame = xd->mi[0].src_mi->mbmi.ref_frame[1];
  struct macroblockd_plane *const tmp_pd = &xd->plane[0];
  struct macroblockd_plane otherpd;

  const MV_SPEED_FEATURES *const mv_sf = &cpi->sf.mv;

  const int tmp_col_min = x->mv_col_min;
  const int tmp_col_max = x->mv_col_max;
  const int tmp_row_min = x->mv_row_min;
  const int tmp_row_max = x->mv_row_max;
  MV ref_full;
  int cost_list[5];
  int sad = INT32_MAX;
  uint8_t tmpbuf[4096];
  BLOCK_SIZE block = bsize == 16 ? BLOCK_16X16 : BLOCK_8X8;
  const vp9_variance_fn_ptr_t v_fn_ptr = cpi->fn_ptr[block];

  // Further step/diamond searches as necessary
  int step_param = mv_sf->reduce_first_step_size;
  step_param = MIN(step_param, MAX_MVSEARCH_STEPS - 2);

  otherpd.dst.buf = tmpbuf;
  xd->plane[0] = otherpd;

  vp9_set_mv_search_range(x, ref_mv);

  ref_full.col = ref_mv->col >> 3;
  ref_full.row = ref_mv->row >> 3;

  /*cpi->sf.search_method == HEX*/
  vp9_hex_search(x, &ref_full, step_param, x->errorperbit, 0,
                 cond_cost_list(cpi, cost_list),
                 &v_fn_ptr, 0, ref_mv, dst_mv);


  // Try sub-pixel MC
  // if (bestsme > error_thresh && bestsme < INT_MAX)
  {
    int distortion;
    unsigned int sse;
    cpi->find_fractional_mv_step(
        x, dst_mv, ref_mv, cpi->common.allow_high_precision_mv, x->errorperbit,
        &v_fn_ptr, 0, mv_sf->subpel_iters_per_step,
        cond_cost_list(cpi, cost_list),
        NULL, NULL,
        &distortion, &sse, NULL, 0, 0);
  }

#if CONFIG_COMPOUND_MODES
  if (has_second_ref(&xd->mi[0].src_mi->mbmi)) {
    xd->mi[0].src_mi->mbmi.mode = NEW_NEWMV;
  } else {
#endif
  xd->mi[0].src_mi->mbmi.mode = NEWMV;
#if CONFIG_COMPOUND_MODES
  }
#endif
  xd->mi[0].src_mi->mbmi.mv[0].as_mv = *dst_mv;
#if CONFIG_INTERINTRA
  xd->mi[0].src_mi->mbmi.ref_frame[1] = NONE;
#endif

  vp9_build_inter_predictors_sby(xd, mb_row, mb_col, block);

  /* restore UMV window */
  x->mv_col_min = tmp_col_min;
  x->mv_col_max = tmp_col_max;
  x->mv_row_min = tmp_row_min;
  x->mv_row_max = tmp_row_max;

  if (bsize == 16) {
    sad = vp9_sad16x16(x->plane[0].src.buf, x->plane[0].src.stride,
                       xd->plane[0].dst.buf, xd->plane[0].dst.stride);
    vp9_variance16x16(x->plane[0].src.buf, x->plane[0].src.stride,
                      xd->plane[0].dst.buf, xd->plane[0].dst.stride, sse);
  } else if (bsize == 8) {
    sad = vp9_sad8x8(x->plane[0].src.buf, x->plane[0].src.stride,
                     xd->plane[0].dst.buf, xd->plane[0].dst.stride);
    vp9_variance8x8(x->plane[0].src.buf, x->plane[0].src.stride,
                    xd->plane[0].dst.buf, xd->plane[0].dst.stride, sse);
  }
  xd->mi[0].src_mi->mbmi.mode = tmp_mode;
  xd->mi[0].src_mi->mbmi.mv[0].as_mv = tmp_mv;
  xd->mi[0].src_mi->mbmi.ref_frame[1] = tmp_frame;

  xd->plane[0] = *tmp_pd;

  return sad;
}

static int do_motion_search(VP9_COMP *cpi, const MV *ref_mv, int bsize,
                            int_mv *dst_mv, int mb_row, int mb_col,
                            unsigned int *sse) {
  MACROBLOCK *const x = &cpi->mb;
  MACROBLOCKD *const xd = &x->e_mbd;
  unsigned int err, tmp_err;
  MV tmp_mv;

  // Try zero MV first
  // FIXME should really use something like near/nearest MV and/or MV prediction
  if (bsize == 16) {
    err = vp9_sad16x16(x->plane[0].src.buf, x->plane[0].src.stride,
                       xd->plane[0].pre[0].buf, xd->plane[0].pre[0].stride);
  } else {
    err = vp9_sad8x8(x->plane[0].src.buf, x->plane[0].src.stride,
                     xd->plane[0].pre[0].buf, xd->plane[0].pre[0].stride);
  }
  dst_mv->as_int = 0;

  // Test last reference frame using the previous best mv as the
  // starting point (best reference) for the search
  tmp_err = do_motion_iteration(cpi, ref_mv, &tmp_mv,
                                      bsize, mb_row, mb_col, sse);
  if (tmp_err < err) {
    err = tmp_err;
    dst_mv->as_mv = tmp_mv;
  }

  // If the current best reference mv is not centered on 0,0 then do a 0,0
  // based search as well.
  if (ref_mv->row != 0 || ref_mv->col != 0) {
    unsigned int tmp_err;
    MV zero_ref_mv = {0, 0}, tmp_mv;

    tmp_err = do_motion_iteration(cpi, &zero_ref_mv, &tmp_mv, bsize,
                                        mb_row, mb_col, sse);
    if (tmp_err < err) {
      dst_mv->as_mv = tmp_mv;
      err = tmp_err;
    }
  }

  return err;
}

static void get_mb_motionfield(VP9_COMP *cpi,
                               YV12_BUFFER_CONFIG *buf,
                               int mb_y_offset,
                               YV12_BUFFER_CONFIG *ref,
                               const MV *prev_ref_mv,
                               int bsize,
                               int mb_row,
                               int mb_col,
                               MV *mv,
                               double *confidence) {
  MACROBLOCK *const x = &cpi->mb;
  MACROBLOCKD *const xd = &x->e_mbd;
  VP9_COMMON *cm = &cpi->common;
  uint8_t *tmp_buf = x->plane[0].src.buf;
  int tmp_stride = x->plane[0].src.stride;
  uint8_t *tmp_dst_buf = xd->plane[0].dst.buf;
  int tmp_dst_stride = xd->plane[0].dst.stride;

  // FIXME in practice we're completely ignoring chroma here
  x->plane[0].src.buf = buf->y_buffer + mb_y_offset;
  x->plane[0].src.stride = buf->y_stride;

  xd->plane[0].dst.buf = get_frame_new_buffer(cm)->y_buffer + mb_y_offset;
  xd->plane[0].dst.stride = get_frame_new_buffer(cm)->y_stride;

  // Golden frame MV search, if it exists and is different than last frame
  if (ref) {
    int_mv intmv;
    unsigned int sse, sad;
    xd->plane[0].pre[0].buf = ref->y_buffer + mb_y_offset;
    xd->plane[0].pre[0].stride = ref->y_stride;
    sad = do_motion_search(cpi,
                           prev_ref_mv,
                           bsize,
                           &intmv,
                           mb_row, mb_col, &sse);
    *confidence = (sse)/(sad+1);
    *mv = intmv.as_mv;
  }

  x->plane[0].src.buf = tmp_buf;
  x->plane[0].src.stride = tmp_stride;

  xd->plane[0].dst.buf = tmp_dst_buf;
  xd->plane[0].dst.stride = tmp_dst_stride;
}

static void get_frame_motionfield(VP9_COMP *cpi,
                                           YV12_BUFFER_CONFIG *buf,
                                           YV12_BUFFER_CONFIG *ref,
                                           int blocksize,
                                           MV *motionfield,
                                           double *confidence) {
  MACROBLOCK *const x = &cpi->mb;
  MACROBLOCKD *const xd = &x->e_mbd;
  VP9_COMMON *const cm = &cpi->common;
  int mb_col, mb_row, offset = 0;
  int mb_y_offset = 0, ref_y_offset = 0;
  int tmp_mv_row_min = x->mv_row_min, tmp_mv_row_max = x->mv_row_max;
  int tmp_up_available = xd->up_available;
  int tmp_left_available = xd->left_available;
  int tmp_y_dst_stride = xd->plane[0].dst.stride;
  int tmp_y_pre_stride = xd->plane[0].pre[0].stride;
  int tmp_uv_dst_stride = xd->plane[1].dst.stride;

  int bsize = blocksize;
  int border = BORDER_MV_PIXELS_B16;

  MV ref_top_mv = {0, 0};
  MODE_INFO mi_local;
  MODE_INFO *tmp_mi = xd->mi[0].src_mi;
  vp9_zero(mi_local);
  // Set up limit values for motion vectors to prevent them extending outside
//  // the UMV borders.
  x->mv_row_min     = -border;
  x->mv_row_max     = (cm->mb_rows - 1) * (bsize/2) + border;
  xd->up_available  = 0;
  xd->plane[0].dst.stride  = buf->y_stride;
  xd->plane[0].pre[0].stride  = buf->y_stride;
  xd->plane[1].dst.stride = buf->uv_stride;
  xd->mi[0].src_mi = &mi_local;
  mi_local.mbmi.sb_type = bsize == 16 ? BLOCK_16X16 : BLOCK_8X8;
  mi_local.mbmi.ref_frame[0] = LAST_FRAME;
  mi_local.mbmi.ref_frame[1] = NONE;
  for (mb_row = 0; mb_row < cm->mb_rows; mb_row++) {
    MV ref_left_mv = ref_top_mv;
    int mb_y_in_offset  = mb_y_offset;
    int ref_y_in_offset = ref_y_offset;

    // Set up limit values for motion vectors to prevent them extending outside
    // the UMV borders.
    x->mv_col_min      = -border;
    x->mv_col_max      = (cm->mb_cols - 1) * (bsize/2) + border;
    xd->left_available = 0;

    for (mb_col = 0; mb_col < cm->mb_cols; mb_col++) {
      MV mv;
      get_mb_motionfield(cpi, buf, mb_y_in_offset,
                                  ref, &ref_left_mv,
                                  blocksize,
                                  mb_row, mb_col, &mv,
                                  &confidence[mb_row*cm->mb_cols + mb_col]);
      motionfield[mb_row*cm->mb_cols + mb_col] = mv;
      if (mb_col == 0) {
        ref_top_mv = ref_left_mv;
      }
      xd->left_available = 1;
      mb_y_in_offset    += bsize;
      ref_y_in_offset   += bsize;
      x->mv_col_min     -= bsize;
      x->mv_col_max     -= bsize;
    }
    xd->up_available = 1;
    mb_y_offset     += buf->y_stride * bsize;
    ref_y_offset    += ref->y_stride * bsize;
    x->mv_row_min   -= bsize;
    x->mv_row_max   -= bsize;
    offset          += cm->mb_cols;
  }
  xd->mi[0].src_mi           = tmp_mi;
  x->mv_row_min              = tmp_mv_row_min;
  x->mv_row_max              = tmp_mv_row_max;
  xd->up_available           = tmp_up_available;
  xd->left_available         = tmp_left_available;
  xd->plane[0].dst.stride    = tmp_y_dst_stride;
  xd->plane[0].pre[0].stride = tmp_y_pre_stride;
  xd->plane[1].dst.stride    = tmp_uv_dst_stride;
}

void vp9_get_motionfield(VP9_COMP *cpi, int ref, int blocksize,
                         MV *motionfield, double *confidence) {
  YV12_BUFFER_CONFIG *ref_buf = get_ref_frame_buffer(cpi, ref);
  struct lookahead_entry *q_cur = vp9_lookahead_peek(cpi->lookahead, 0);

  if (q_cur) {
    get_frame_motionfield(cpi, &q_cur->img, ref_buf,
                                   blocksize, motionfield, confidence);
  }
}
