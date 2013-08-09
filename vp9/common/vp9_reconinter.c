/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include <assert.h>

#include "./vpx_config.h"
#include "vpx/vpx_integer.h"
#include "vp9/common/vp9_blockd.h"
#include "vp9/common/vp9_filter.h"
#include "vp9/common/vp9_reconinter.h"
#include "vp9/common/vp9_reconintra.h"
#include "./vpx_scale_rtcd.h"

static int scale_value_x_with_scaling(int val,
                                      const struct scale_factors *scale) {
  return (val * scale->x_scale_fp >> VP9_REF_SCALE_SHIFT);
}

static int scale_value_y_with_scaling(int val,
                                      const struct scale_factors *scale) {
  return (val * scale->y_scale_fp >> VP9_REF_SCALE_SHIFT);
}

static int unscaled_value(int val, const struct scale_factors *scale) {
  (void) scale;
  return val;
}

static MV32 mv_with_scaling(const MV *mv,
                               const struct scale_factors *scale) {
  const MV32 res = {
    (mv->row * scale->y_scale_fp >> VP9_REF_SCALE_SHIFT) + scale->y_offset_q4,
    (mv->col * scale->x_scale_fp >> VP9_REF_SCALE_SHIFT) + scale->x_offset_q4
  };
  return res;
}

static MV32 mv_without_scaling(const MV *mv,
                               const struct scale_factors *scale) {
  const MV32 res = {
    mv->row,
    mv->col
  };
  return res;
}

static void set_offsets_with_scaling(struct scale_factors *scale,
                                     int row, int col) {
  const int x_q4 = 16 * col;
  const int y_q4 = 16 * row;

  scale->x_offset_q4 = (x_q4 * scale->x_scale_fp >> VP9_REF_SCALE_SHIFT) & 0xf;
  scale->y_offset_q4 = (y_q4 * scale->y_scale_fp >> VP9_REF_SCALE_SHIFT) & 0xf;
}

static void set_offsets_without_scaling(struct scale_factors *scale,
                                        int row, int col) {
  scale->x_offset_q4 = 0;
  scale->y_offset_q4 = 0;
}

static int get_fixed_point_scale_factor(int other_size, int this_size) {
  // Calculate scaling factor once for each reference frame
  // and use fixed point scaling factors in decoding and encoding routines.
  // Hardware implementations can calculate scale factor in device driver
  // and use multiplication and shifting on hardware instead of division.
  return (other_size << VP9_REF_SCALE_SHIFT) / this_size;
}

void vp9_setup_scale_factors_for_frame(struct scale_factors *scale,
                                       int other_w, int other_h,
                                       int this_w, int this_h) {
  scale->x_scale_fp = get_fixed_point_scale_factor(other_w, this_w);
  scale->x_offset_q4 = 0;  // calculated per-mb
  scale->x_step_q4 = (16 * scale->x_scale_fp >> VP9_REF_SCALE_SHIFT);

  scale->y_scale_fp = get_fixed_point_scale_factor(other_h, this_h);
  scale->y_offset_q4 = 0;  // calculated per-mb
  scale->y_step_q4 = (16 * scale->y_scale_fp >> VP9_REF_SCALE_SHIFT);

  if ((other_w == this_w) && (other_h == this_h)) {
    scale->scale_value_x = unscaled_value;
    scale->scale_value_y = unscaled_value;
    scale->set_scaled_offsets = set_offsets_without_scaling;
    scale->scale_mv = mv_without_scaling;
  } else {
    scale->scale_value_x = scale_value_x_with_scaling;
    scale->scale_value_y = scale_value_y_with_scaling;
    scale->set_scaled_offsets = set_offsets_with_scaling;
    scale->scale_mv = mv_with_scaling;
  }

  // TODO(agrange): Investigate the best choice of functions to use here
  // for EIGHTTAP_SMOOTH. Since it is not interpolating, need to choose what
  // to do at full-pel offsets. The current selection, where the filter is
  // applied in one direction only, and not at all for 0,0, seems to give the
  // best quality, but it may be worth trying an additional mode that does
  // do the filtering on full-pel.
  if (scale->x_step_q4 == 16) {
    if (scale->y_step_q4 == 16) {
      // No scaling in either direction.
      scale->predict[0][0][0] = vp9_convolve_copy;
      scale->predict[0][0][1] = vp9_convolve_avg;
      scale->predict[0][1][0] = vp9_convolve8_vert;
      scale->predict[0][1][1] = vp9_convolve8_avg_vert;
      scale->predict[1][0][0] = vp9_convolve8_horiz;
      scale->predict[1][0][1] = vp9_convolve8_avg_horiz;
    } else {
      // No scaling in x direction. Must always scale in the y direction.
      scale->predict[0][0][0] = vp9_convolve8_vert;
      scale->predict[0][0][1] = vp9_convolve8_avg_vert;
      scale->predict[0][1][0] = vp9_convolve8_vert;
      scale->predict[0][1][1] = vp9_convolve8_avg_vert;
      scale->predict[1][0][0] = vp9_convolve8;
      scale->predict[1][0][1] = vp9_convolve8_avg;
    }
  } else {
    if (scale->y_step_q4 == 16) {
      // No scaling in the y direction. Must always scale in the x direction.
      scale->predict[0][0][0] = vp9_convolve8_horiz;
      scale->predict[0][0][1] = vp9_convolve8_avg_horiz;
      scale->predict[0][1][0] = vp9_convolve8;
      scale->predict[0][1][1] = vp9_convolve8_avg;
      scale->predict[1][0][0] = vp9_convolve8_horiz;
      scale->predict[1][0][1] = vp9_convolve8_avg_horiz;
    } else {
      // Must always scale in both directions.
      scale->predict[0][0][0] = vp9_convolve8;
      scale->predict[0][0][1] = vp9_convolve8_avg;
      scale->predict[0][1][0] = vp9_convolve8;
      scale->predict[0][1][1] = vp9_convolve8_avg;
      scale->predict[1][0][0] = vp9_convolve8;
      scale->predict[1][0][1] = vp9_convolve8_avg;
    }
  }
  // 2D subpel motion always gets filtered in both directions
  scale->predict[1][1][0] = vp9_convolve8;
  scale->predict[1][1][1] = vp9_convolve8_avg;
}

void vp9_setup_interp_filters(MACROBLOCKD *xd,
                              INTERPOLATIONFILTERTYPE mcomp_filter_type,
                              VP9_COMMON *cm) {
  if (xd->mode_info_context) {
    MB_MODE_INFO *mbmi = &xd->mode_info_context->mbmi;

    set_scale_factors(xd, mbmi->ref_frame[0] - 1, mbmi->ref_frame[1] - 1,
                      cm->active_ref_scale);
  }

  switch (mcomp_filter_type) {
    case EIGHTTAP:
    case SWITCHABLE:
      xd->subpix.filter_x = xd->subpix.filter_y = vp9_sub_pel_filters_8;
      break;
    case EIGHTTAP_SMOOTH:
      xd->subpix.filter_x = xd->subpix.filter_y = vp9_sub_pel_filters_8lp;
      break;
    case EIGHTTAP_SHARP:
      xd->subpix.filter_x = xd->subpix.filter_y = vp9_sub_pel_filters_8s;
      break;
    case BILINEAR:
      xd->subpix.filter_x = xd->subpix.filter_y = vp9_bilinear_filters;
      break;
  }
  assert(((intptr_t)xd->subpix.filter_x & 0xff) == 0);
}

void vp9_build_inter_predictor(const uint8_t *src, int src_stride,
                               uint8_t *dst, int dst_stride,
                               const MV *src_mv,
                               const struct scale_factors *scale,
                               int w, int h, int weight,
                               const struct subpix_fn_table *subpix,
                               enum mv_precision precision) {
  const int is_q4 = precision == MV_PRECISION_Q4;
  const MV mv_q4 = { is_q4 ? src_mv->row : src_mv->row << 1,
                     is_q4 ? src_mv->col : src_mv->col << 1 };
  const MV32 mv = scale->scale_mv(&mv_q4, scale);
  const int subpel_x = mv.col & SUBPEL_MASK;
  const int subpel_y = mv.row & SUBPEL_MASK;

  src += (mv.row >> SUBPEL_BITS) * src_stride + (mv.col >> SUBPEL_BITS);
  scale->predict[subpel_x != 0][subpel_y != 0][weight](
      src, src_stride, dst, dst_stride,
      subpix->filter_x[subpel_x], scale->x_step_q4,
      subpix->filter_y[subpel_y], scale->y_step_q4,
      w, h);
}

static INLINE int round_mv_comp_q4(int value) {
  return (value < 0 ? value - 2 : value + 2) / 4;
}

static MV mi_mv_pred_q4(const MODE_INFO *mi, int idx) {
  MV res = { round_mv_comp_q4(mi->bmi[0].as_mv[idx].as_mv.row +
                              mi->bmi[1].as_mv[idx].as_mv.row +
                              mi->bmi[2].as_mv[idx].as_mv.row +
                              mi->bmi[3].as_mv[idx].as_mv.row),
             round_mv_comp_q4(mi->bmi[0].as_mv[idx].as_mv.col +
                              mi->bmi[1].as_mv[idx].as_mv.col +
                              mi->bmi[2].as_mv[idx].as_mv.col +
                              mi->bmi[3].as_mv[idx].as_mv.col) };
  return res;
}

// TODO(jkoleszar): yet another mv clamping function :-(
MV clamp_mv_to_umv_border_sb(const MACROBLOCKD *xd, const MV *src_mv,
                             int bw, int bh, int ss_x, int ss_y) {
  // If the MV points so far into the UMV border that no visible pixels
  // are used for reconstruction, the subpel part of the MV can be
  // discarded and the MV limited to 16 pixels with equivalent results.
  const int spel_left = (VP9_INTERP_EXTEND + bw) << SUBPEL_BITS;
  const int spel_right = spel_left - SUBPEL_SHIFTS;
  const int spel_top = (VP9_INTERP_EXTEND + bh) << SUBPEL_BITS;
  const int spel_bottom = spel_top - SUBPEL_SHIFTS;
  MV clamped_mv = {
    src_mv->row << (1 - ss_y),
    src_mv->col << (1 - ss_x)
  };
  assert(ss_x <= 1);
  assert(ss_y <= 1);

  clamp_mv(&clamped_mv, (xd->mb_to_left_edge << (1 - ss_x)) - spel_left,
                        (xd->mb_to_right_edge << (1 - ss_x)) + spel_right,
                        (xd->mb_to_top_edge << (1 - ss_y)) - spel_top,
                        (xd->mb_to_bottom_edge << (1 - ss_y)) + spel_bottom);

  return clamped_mv;
}

struct build_inter_predictors_args {
  MACROBLOCKD *xd;
  int x;
  int y;
  struct buf_2d *dst[MAX_MB_PLANE];
  struct buf_2d *pre[2][MAX_MB_PLANE];
};
static void build_inter_predictors(int plane, int block,
                                   BLOCK_SIZE_TYPE bsize,
                                   int pred_w, int pred_h,
                                   void *argv) {
  const struct build_inter_predictors_args* const arg = argv;
  MACROBLOCKD *const xd = arg->xd;
  struct macroblockd_plane *const pd = &xd->plane[plane];
  const int bwl = b_width_log2(bsize) - pd->subsampling_x;
  const int bw = 4 << bwl;
  const int bh = plane_block_height(bsize, pd);
  const int x = 4 * (block & ((1 << bwl) - 1));
  const int y = 4 * (block >> bwl);
  const MODE_INFO *const mi = xd->mode_info_context;
  const int use_second_ref = mi->mbmi.ref_frame[1] > 0;
  int which_mv;

  assert(x < bw);
  assert(y < bh);
  assert(mi->mbmi.sb_type < BLOCK_8X8 || 4 << pred_w == bw);
  assert(mi->mbmi.sb_type < BLOCK_8X8 || 4 << pred_h == bh);

  for (which_mv = 0; which_mv < 1 + use_second_ref; ++which_mv) {
    struct scale_factors *const scale = &xd->scale_factor[which_mv];
    struct buf_2d *const pre_buf = arg->pre[which_mv][plane];
    struct buf_2d *const dst_buf = arg->dst[plane];

    const uint8_t *const pre = pre_buf->buf + scaled_buffer_offset(x, y,
                               pre_buf->stride, scale);

    uint8_t *const dst = dst_buf->buf + dst_buf->stride * y + x;

    // TODO(jkoleszar): All chroma MVs in SPLITMV mode are taken as the
    // same MV (the average of the 4 luma MVs) but we could do something
    // smarter for non-4:2:0. Just punt for now, pending the changes to get
    // rid of SPLITMV mode entirely.
    const MV mv = mi->mbmi.sb_type < BLOCK_8X8
               ? (plane == 0 ? mi->bmi[block].as_mv[which_mv].as_mv
                             : mi_mv_pred_q4(mi, which_mv))
               : mi->mbmi.mv[which_mv].as_mv;

    // TODO(jkoleszar): This clamping is done in the incorrect place for the
    // scaling case. It needs to be done on the scaled MV, not the pre-scaling
    // MV. Note however that it performs the subsampling aware scaling so
    // that the result is always q4.
    const MV res_mv = clamp_mv_to_umv_border_sb(xd, &mv, bw, bh,
                                                pd->subsampling_x,
                                                pd->subsampling_y);

    scale->set_scaled_offsets(scale, arg->y + y, arg->x + x);
    vp9_build_inter_predictor(pre, pre_buf->stride, dst, dst_buf->stride,
                              &res_mv, scale,
                              4 << pred_w, 4 << pred_h, which_mv,
                              &xd->subpix, MV_PRECISION_Q4);
  }
}
void vp9_build_inter_predictors_sby(MACROBLOCKD *xd, int mi_row, int mi_col,
                                    BLOCK_SIZE_TYPE bsize) {
  struct build_inter_predictors_args args = {
    xd, mi_col * MI_SIZE, mi_row * MI_SIZE,
    {&xd->plane[0].dst, NULL, NULL},
    {{&xd->plane[0].pre[0], NULL, NULL},
     {&xd->plane[0].pre[1], NULL, NULL}},
  };

  foreach_predicted_block_in_plane(xd, bsize, 0, build_inter_predictors, &args);
}
void vp9_build_inter_predictors_sbuv(MACROBLOCKD *xd, int mi_row, int mi_col,
                                     BLOCK_SIZE_TYPE bsize) {
  struct build_inter_predictors_args args = {
    xd, mi_col * MI_SIZE, mi_row * MI_SIZE,
#if CONFIG_ALPHA
    {NULL, &xd->plane[1].dst, &xd->plane[2].dst, &xd->plane[3].dst},
    {{NULL, &xd->plane[1].pre[0], &xd->plane[2].pre[0], &xd->plane[3].pre[0]},
     {NULL, &xd->plane[1].pre[1], &xd->plane[2].pre[1], &xd->plane[3].pre[1]}},
#else
    {NULL, &xd->plane[1].dst, &xd->plane[2].dst},
    {{NULL, &xd->plane[1].pre[0], &xd->plane[2].pre[0]},
     {NULL, &xd->plane[1].pre[1], &xd->plane[2].pre[1]}},
#endif
  };
  foreach_predicted_block_uv(xd, bsize, build_inter_predictors, &args);
}
void vp9_build_inter_predictors_sb(MACROBLOCKD *xd,
                                   int mi_row, int mi_col,
                                   BLOCK_SIZE_TYPE bsize) {

  vp9_build_inter_predictors_sby(xd, mi_row, mi_col, bsize);
  vp9_build_inter_predictors_sbuv(xd, mi_row, mi_col, bsize);
}

// TODO(dkovalev: find better place for this function)
void vp9_setup_scale_factors(VP9_COMMON *cm, int i) {
  const int ref = cm->active_ref_idx[i];
  struct scale_factors *const sf = &cm->active_ref_scale[i];
  if (ref >= NUM_YV12_BUFFERS) {
    vp9_zero(*sf);
  } else {
    YV12_BUFFER_CONFIG *const fb = &cm->yv12_fb[ref];
    vp9_setup_scale_factors_for_frame(sf,
                                      fb->y_crop_width, fb->y_crop_height,
                                      cm->width, cm->height);

    if (sf->x_scale_fp != VP9_REF_NO_SCALE ||
        sf->y_scale_fp != VP9_REF_NO_SCALE)
      vp9_extend_frame_borders(fb, cm->subsampling_x, cm->subsampling_y);
  }
}

