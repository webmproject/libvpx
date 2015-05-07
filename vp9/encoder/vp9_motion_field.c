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

unsigned int get_sby_perpixel_variance(VP9_COMP *cpi,
                                       const struct buf_2d *ref,
                                       BLOCK_SIZE bs);
unsigned int get_sby_perpixel_ssd(VP9_COMP *cpi,
                                  const struct buf_2d *ref,
                                  BLOCK_SIZE bs);
#if CONFIG_VP9_HIGHBITDEPTH
unsigned int high_get_sby_perpixel_ssd(
    VP9_COMP *cpi, const struct buf_2d *ref, BLOCK_SIZE bs, int bd);
unsigned int high_get_sby_perpixel_variance(
    VP9_COMP *cpi, const struct buf_2d *ref, BLOCK_SIZE bs, int bd);
#endif

static int do_motion_iteration(MACROBLOCK *const x,
                               const vp9_variance_fn_ptr_t *v_fn_ptr,
                               const MV *ref_mv,
                               MV *dst_mv) {
  MV ref_full;
  unsigned int sse;
  int besterr, distortion;

  // Further step/diamond searches as necessary
  // int step_param = mv_sf->reduce_first_step_size;
  // step_param = MIN(step_param, MAX_MVSEARCH_STEPS - 2);
  int step_param = 0;
  int subpel_iters_per_step = 2;
  int allow_high_precision_mv = 1;

  vp9_set_mv_search_range(x, ref_mv);

  ref_full.col = ref_mv->col >> 3;
  ref_full.row = ref_mv->row >> 3;

  /*cpi->sf.search_method == HEX*/
  vp9_hex_search(x, &ref_full, step_param, x->errorperbit, 0, NULL,
                 v_fn_ptr, 0, ref_mv, dst_mv);

  besterr = vp9_find_best_sub_pixel_tree(
      x, dst_mv, ref_mv, allow_high_precision_mv, x->errorperbit,
      v_fn_ptr, 0, subpel_iters_per_step,
      NULL, NULL, NULL, &distortion, &sse, NULL, 0, 0);

  return besterr;
}

static int do_motion_search(MACROBLOCK *const x,
                            const vp9_variance_fn_ptr_t *v_fn_ptr,
                            const MV *ref_mv,
                            int_mv *dst_mv) {
  int err = do_motion_iteration(x, v_fn_ptr,
                                ref_mv, &dst_mv->as_mv);

  // If the current best reference mv is not centered on 0,0 then do a 0,0
  // based search as well.
  if (ref_mv->row != 0 || ref_mv->col != 0) {
    int tmp_err;
    MV zero_ref_mv = {0, 0}, tmp_mv;

    tmp_err = do_motion_iteration(x, v_fn_ptr,
                                  &zero_ref_mv, &tmp_mv);
    if (tmp_err < err) {
      dst_mv->as_mv = tmp_mv;
      err = tmp_err;
    }
  }
  return err;
}

void vp9_get_frame_motionfield(struct VP9_COMP *cpi,
                               YV12_BUFFER_CONFIG *buf,
                               YV12_BUFFER_CONFIG *ref,
                               BLOCK_SIZE bsize,
                               MV *motionfield,
                               double *confidence) {
  MACROBLOCK mx;
  MACROBLOCK *const x = &mx;
  MACROBLOCKD *const xd = &x->e_mbd;
  VP9_COMMON *const cm = &cpi->common;
  int mb_col, mb_row;
  int mb_y_offset = 0, ref_y_offset = 0;

  int border = BORDER_MV_PIXELS_B16;
  int bwidth = num_4x4_blocks_wide_lookup[bsize] << 2;
  int bheight = num_4x4_blocks_high_lookup[bsize] << 2;

  MV ref_top_mv = {0, 0};

  x->errorperbit =
      vp9_compute_rd_mult(cpi, cm->base_qindex + cm->y_dc_delta_q) /
      64;
  x->errorperbit += (x->errorperbit == 0);

  // the UMV borders.
  x->mv_row_min     = -border;
  x->mv_row_max     = cm->mi_rows * 8 + border;
  xd->up_available  = 0;
  xd->plane[0].dst.stride  = buf->y_stride;
  xd->plane[1].dst.stride = buf->uv_stride;
  xd->plane[0].pre[0].stride  = buf->y_stride;
  for (mb_row = 0; mb_row < cm->mb_rows; mb_row++) {
    MV ref_left_mv = ref_top_mv;
    int mb_y_in_offset  = mb_y_offset;
    int ref_y_in_offset = ref_y_offset;

    // Set up limit values for motion vectors to prevent them extending outside
    // the UMV borders.
    x->mv_col_min      = -border;
    x->mv_col_max      = cm->mi_cols * 8 + border;
    xd->left_available = 0;

    for (mb_col = 0; mb_col < cm->mb_cols; mb_col++) {
      int_mv intmv;
      unsigned int ssd, err;

      x->plane[0].src.buf = buf->y_buffer + mb_y_in_offset;
      x->plane[0].src.stride = buf->y_stride;
      xd->plane[0].pre[0].buf = ref->y_buffer + ref_y_in_offset;
      xd->plane[0].pre[0].stride = ref->y_stride;
      ssd = get_sby_perpixel_ssd(cpi, &x->plane[0].src, bsize);
      err =
          do_motion_search(x,
                           &cpi->fn_ptr[bsize],
                           &ref_left_mv,
                           &intmv);
      confidence[mb_row * cm->mb_cols + mb_col] = (double)ssd / (err + 1);
      motionfield[mb_row * cm->mb_cols + mb_col] = intmv.as_mv;
      if (mb_col == 0) {
        ref_top_mv = ref_left_mv;
      }
      xd->left_available = 1;
      mb_y_in_offset    += bwidth;
      ref_y_in_offset   += bwidth;
      x->mv_col_min     -= bwidth;
      x->mv_col_max     -= bwidth;
    }
    xd->up_available = 1;
    mb_y_offset     += buf->y_stride * bheight;
    ref_y_offset    += ref->y_stride * bheight;
    x->mv_row_min   -= bheight;
    x->mv_row_max   -= bheight;
  }
}

void vp9_get_ref_motionfield(VP9_COMP *cpi, int ref, BLOCK_SIZE bsize,
                             MV *motionfield, double *confidence) {
  YV12_BUFFER_CONFIG *ref_buf = get_ref_frame_buffer(cpi, ref);

  if (cpi->Source && ref_buf) {
    vp9_get_frame_motionfield(cpi, cpi->Source, ref_buf,
                              bsize, motionfield, confidence);
  }
}
