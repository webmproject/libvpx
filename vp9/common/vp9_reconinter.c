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

static MV32 mv_q3_to_q4_with_scaling(const MV *mv,
                                     const struct scale_factors *scale) {
  const MV32 res = {
    ((mv->row << 1) * scale->y_scale_fp >> VP9_REF_SCALE_SHIFT)
        + scale->y_offset_q4,
    ((mv->col << 1) * scale->x_scale_fp >> VP9_REF_SCALE_SHIFT)
        + scale->x_offset_q4
  };
  return res;
}

static MV32 mv_q3_to_q4_without_scaling(const MV *mv,
                                        const struct scale_factors *scale) {
  const MV32 res = {
     mv->row << 1,
     mv->col << 1
  };
  return res;
}

static MV32 mv_q4_with_scaling(const MV *mv,
                               const struct scale_factors *scale) {
  const MV32 res = {
    (mv->row * scale->y_scale_fp >> VP9_REF_SCALE_SHIFT) + scale->y_offset_q4,
    (mv->col * scale->x_scale_fp >> VP9_REF_SCALE_SHIFT) + scale->x_offset_q4
  };
  return res;
}

static MV32 mv_q4_without_scaling(const MV *mv,
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
    scale->scale_mv_q3_to_q4 = mv_q3_to_q4_without_scaling;
    scale->scale_mv_q4 = mv_q4_without_scaling;
  } else {
    scale->scale_value_x = scale_value_x_with_scaling;
    scale->scale_value_y = scale_value_y_with_scaling;
    scale->set_scaled_offsets = set_offsets_with_scaling;
    scale->scale_mv_q3_to_q4 = mv_q3_to_q4_with_scaling;
    scale->scale_mv_q4 = mv_q4_with_scaling;
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
  const MV32 mv = precision == MV_PRECISION_Q4
                     ? scale->scale_mv_q4(src_mv, scale)
                     : scale->scale_mv_q3_to_q4(src_mv, scale);
  const int subpel_x = mv.col & 15;
  const int subpel_y = mv.row & 15;

  src += (mv.row >> 4) * src_stride + (mv.col >> 4);
  scale->predict[!!subpel_x][!!subpel_y][weight](
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
MV clamp_mv_to_umv_border_sb(const MV *src_mv,
    int bwl, int bhl, int ss_x, int ss_y,
    int mb_to_left_edge, int mb_to_top_edge,
    int mb_to_right_edge, int mb_to_bottom_edge) {
  // If the MV points so far into the UMV border that no visible pixels
  // are used for reconstruction, the subpel part of the MV can be
  // discarded and the MV limited to 16 pixels with equivalent results.
  const int spel_left = (VP9_INTERP_EXTEND + (4 << bwl)) << 4;
  const int spel_right = spel_left - (1 << 4);
  const int spel_top = (VP9_INTERP_EXTEND + (4 << bhl)) << 4;
  const int spel_bottom = spel_top - (1 << 4);
  MV clamped_mv = {
    src_mv->row << (1 - ss_y),
    src_mv->col << (1 - ss_x)
  };
  assert(ss_x <= 1);
  assert(ss_y <= 1);

  clamp_mv(&clamped_mv, (mb_to_left_edge << (1 - ss_x)) - spel_left,
                        (mb_to_right_edge << (1 - ss_x)) + spel_right,
                        (mb_to_top_edge << (1 - ss_y)) - spel_top,
                        (mb_to_bottom_edge << (1 - ss_y)) + spel_bottom);

  return clamped_mv;
}

#if CONFIG_MASKED_COMPOUND
#define MASK_WEIGHT_BITS 6

static int get_masked_weight(int m) {
  #define SMOOTHER_LEN  32
  static const uint8_t smoothfn[2 * SMOOTHER_LEN + 1] = {
      0,  0,  0,  0,  0,  0,  0,  0,
      0,  0,  0,  0,  0,  1,  1,  1,
      1,  1,  2,  2,  3,  4,  5,  6,
      8,  9, 12, 14, 17, 21, 24, 28,
      32,
      36, 40, 43, 47, 50, 52, 55, 56,
      58, 59, 60, 61, 62, 62, 63, 63,
      63, 63, 63, 64, 64, 64, 64, 64,
      64, 64, 64, 64, 64, 64, 64, 64,
  };
  if (m < -SMOOTHER_LEN)
    return 0;
  else if (m > SMOOTHER_LEN)
    return (1 << MASK_WEIGHT_BITS);
  else
    return smoothfn[m + SMOOTHER_LEN];
}

static int get_hard_mask(int m) {
  return m > 0;
}

// Equation of line: f(x, y) = a[0]*(x - a[2]*w/4) + a[1]*(y - a[3]*h/4) = 0
// The soft mask is obtained by computing f(x, y) and then calling
// get_masked_weight(f(x, y)).
static const int mask_params_sml[1 << MASK_BITS_SML][4] = {
  {-1,  2, 2, 2},
  { 1, -2, 2, 2},
  {-2,  1, 2, 2},
  { 2, -1, 2, 2},
  { 2,  1, 2, 2},
  {-2, -1, 2, 2},
  { 1,  2, 2, 2},
  {-1, -2, 2, 2},
};

static const int mask_params_med_hgtw[1 << MASK_BITS_MED][4] = {
  {-1,  2, 2, 2},
  { 1, -2, 2, 2},
  {-2,  1, 2, 2},
  { 2, -1, 2, 2},
  { 2,  1, 2, 2},
  {-2, -1, 2, 2},
  { 1,  2, 2, 2},
  {-1, -2, 2, 2},

  {-1,  2, 2, 1},
  { 1, -2, 2, 1},
  {-1,  2, 2, 3},
  { 1, -2, 2, 3},
  { 1,  2, 2, 1},
  {-1, -2, 2, 1},
  { 1,  2, 2, 3},
  {-1, -2, 2, 3},
};

static const int mask_params_med_hltw[1 << MASK_BITS_MED][4] = {
  {-1,  2, 2, 2},
  { 1, -2, 2, 2},
  {-2,  1, 2, 2},
  { 2, -1, 2, 2},
  { 2,  1, 2, 2},
  {-2, -1, 2, 2},
  { 1,  2, 2, 2},
  {-1, -2, 2, 2},

  {-2,  1, 1, 2},
  { 2, -1, 1, 2},
  {-2,  1, 3, 2},
  { 2, -1, 3, 2},
  { 2,  1, 1, 2},
  {-2, -1, 1, 2},
  { 2,  1, 3, 2},
  {-2, -1, 3, 2},
};

static const int mask_params_med_heqw[1 << MASK_BITS_MED][4] = {
  {-1,  2, 2, 2},
  { 1, -2, 2, 2},
  {-2,  1, 2, 2},
  { 2, -1, 2, 2},
  { 2,  1, 2, 2},
  {-2, -1, 2, 2},
  { 1,  2, 2, 2},
  {-1, -2, 2, 2},

  { 0,  2, 0, 1},
  { 0, -2, 0, 1},
  { 0,  2, 0, 3},
  { 0, -2, 0, 3},
  { 2,  0, 1, 0},
  {-2,  0, 1, 0},
  { 2,  0, 3, 0},
  {-2,  0, 3, 0},
};

static const int mask_params_big_hgtw[1 << MASK_BITS_BIG][4] = {
  {-1,  2, 2, 2},
  { 1, -2, 2, 2},
  {-2,  1, 2, 2},
  { 2, -1, 2, 2},
  { 2,  1, 2, 2},
  {-2, -1, 2, 2},
  { 1,  2, 2, 2},
  {-1, -2, 2, 2},

  {-1,  2, 2, 1},
  { 1, -2, 2, 1},
  {-1,  2, 2, 3},
  { 1, -2, 2, 3},
  { 1,  2, 2, 1},
  {-1, -2, 2, 1},
  { 1,  2, 2, 3},
  {-1, -2, 2, 3},

  {-2,  1, 1, 2},
  { 2, -1, 1, 2},
  {-2,  1, 3, 2},
  { 2, -1, 3, 2},
  { 2,  1, 1, 2},
  {-2, -1, 1, 2},
  { 2,  1, 3, 2},
  {-2, -1, 3, 2},

  { 0,  2, 0, 1},
  { 0, -2, 0, 1},
  { 0,  2, 0, 2},
  { 0, -2, 0, 2},
  { 0,  2, 0, 3},
  { 0, -2, 0, 3},
  { 2,  0, 2, 0},
  {-2,  0, 2, 0},
};

static const int mask_params_big_hltw[1 << MASK_BITS_BIG][4] = {
  {-1,  2, 2, 2},
  { 1, -2, 2, 2},
  {-2,  1, 2, 2},
  { 2, -1, 2, 2},
  { 2,  1, 2, 2},
  {-2, -1, 2, 2},
  { 1,  2, 2, 2},
  {-1, -2, 2, 2},

  {-1,  2, 2, 1},
  { 1, -2, 2, 1},
  {-1,  2, 2, 3},
  { 1, -2, 2, 3},
  { 1,  2, 2, 1},
  {-1, -2, 2, 1},
  { 1,  2, 2, 3},
  {-1, -2, 2, 3},

  {-2,  1, 1, 2},
  { 2, -1, 1, 2},
  {-2,  1, 3, 2},
  { 2, -1, 3, 2},
  { 2,  1, 1, 2},
  {-2, -1, 1, 2},
  { 2,  1, 3, 2},
  {-2, -1, 3, 2},

  { 0,  2, 0, 2},
  { 0, -2, 0, 2},
  { 2,  0, 1, 0},
  {-2,  0, 1, 0},
  { 2,  0, 2, 0},
  {-2,  0, 2, 0},
  { 2,  0, 3, 0},
  {-2,  0, 3, 0},
};

static const int mask_params_big_heqw[1 << MASK_BITS_BIG][4] = {
  {-1,  2, 2, 2},
  { 1, -2, 2, 2},
  {-2,  1, 2, 2},
  { 2, -1, 2, 2},
  { 2,  1, 2, 2},
  {-2, -1, 2, 2},
  { 1,  2, 2, 2},
  {-1, -2, 2, 2},

  {-1,  2, 2, 1},
  { 1, -2, 2, 1},
  {-1,  2, 2, 3},
  { 1, -2, 2, 3},
  { 1,  2, 2, 1},
  {-1, -2, 2, 1},
  { 1,  2, 2, 3},
  {-1, -2, 2, 3},

  {-2,  1, 1, 2},
  { 2, -1, 1, 2},
  {-2,  1, 3, 2},
  { 2, -1, 3, 2},
  { 2,  1, 1, 2},
  {-2, -1, 1, 2},
  { 2,  1, 3, 2},
  {-2, -1, 3, 2},

  { 0,  2, 0, 1},
  { 0, -2, 0, 1},
  { 0,  2, 0, 3},
  { 0, -2, 0, 3},
  { 2,  0, 1, 0},
  {-2,  0, 1, 0},
  { 2,  0, 3, 0},
  {-2,  0, 3, 0},
};

static const int *get_mask_params(int mask_index,
                                  BLOCK_SIZE_TYPE sb_type,
                                  int h, int w) {
  const int *a;
  const int mask_bits = get_mask_bits(sb_type);

  if (mask_index == MASK_NONE)
    return NULL;

  if (mask_bits == MASK_BITS_SML) {
    a = mask_params_sml[mask_index];
  } else if (mask_bits == MASK_BITS_MED) {
    if (h > w)
      a = mask_params_med_hgtw[mask_index];
    else if (h < w)
      a = mask_params_med_hltw[mask_index];
    else
      a = mask_params_med_heqw[mask_index];
  } else if (mask_bits == MASK_BITS_BIG) {
    if (h > w)
      a = mask_params_big_hgtw[mask_index];
    else if (h < w)
      a = mask_params_big_hltw[mask_index];
    else
      a = mask_params_big_heqw[mask_index];
  } else {
    assert(0);
  }
  return a;
}

void vp9_generate_masked_weight(int mask_index,
                                BLOCK_SIZE_TYPE sb_type,
                                int h, int w,
                                uint8_t *mask, int stride) {
  int i, j;
  const int *a = get_mask_params(mask_index, sb_type, h, w);
  if (!a) return;
  for (i = 0; i < h; ++i)
    for (j = 0; j < w; ++j) {
      int x = (j - (a[2] * w) / 4);
      int y = (i - (a[3] * h) / 4);
      int m = a[0] * x + a[1] * y;
      mask[i * stride + j] = get_masked_weight(m);
    }
}

void vp9_generate_hard_mask(int mask_index, BLOCK_SIZE_TYPE sb_type,
                            int h, int w, uint8_t *mask, int stride) {
  int i, j;
  const int *a = get_mask_params(mask_index, sb_type, h, w);
  if (!a) return;
  for (i = 0; i < h; ++i)
    for (j = 0; j < w; ++j) {
      int x = (j - (a[2] * w) / 4);
      int y = (i - (a[3] * h) / 4);
      int m = a[0] * x + a[1] * y;
      mask[i * stride + j] = get_hard_mask(m);
    }
}

static void build_masked_compound(uint8_t *dst, int dst_stride,
                                  uint8_t *dst2, int dst2_stride,
                                  int mask_index, BLOCK_SIZE_TYPE sb_type,
                                  int h, int w) {
  int i, j;
  uint8_t mask[4096];
  vp9_generate_masked_weight(mask_index, sb_type, h, w, mask, 64);
  for (i = 0; i < h; ++i)
    for (j = 0; j < w; ++j) {
      int m = mask[i * 64 + j];
      dst[i * dst_stride + j] =  (dst[i * dst_stride + j] * m +
                                  dst2[i * dst2_stride + j] *
                                  ((1 << MASK_WEIGHT_BITS) - m) +
                                  (1 << (MASK_WEIGHT_BITS - 1))) >>
                                 MASK_WEIGHT_BITS;
    }
}
#endif

struct build_inter_predictors_args {
  MACROBLOCKD *xd;
  int x;
  int y;
  uint8_t* dst[MAX_MB_PLANE];
  int dst_stride[MAX_MB_PLANE];
  uint8_t* pre[2][MAX_MB_PLANE];
  int pre_stride[2][MAX_MB_PLANE];
};
static void build_inter_predictors(int plane, int block,
                                   BLOCK_SIZE_TYPE bsize,
                                   int pred_w, int pred_h,
                                   void *argv) {
  const struct build_inter_predictors_args* const arg = argv;
  MACROBLOCKD * const xd = arg->xd;
  const int bwl = b_width_log2(bsize) - xd->plane[plane].subsampling_x;
  const int bhl = b_height_log2(bsize) - xd->plane[plane].subsampling_y;
  const int x = 4 * (block & ((1 << bwl) - 1)), y = 4 * (block >> bwl);
  const MODE_INFO *const mi = xd->mode_info_context;
  const int use_second_ref = mi->mbmi.ref_frame[1] > 0;
  int which_mv;

  assert(x < (4 << bwl));
  assert(y < (4 << bhl));
  assert(mi->mbmi.sb_type < BLOCK_8X8 || 4 << pred_w == (4 << bwl));
  assert(mi->mbmi.sb_type < BLOCK_8X8 || 4 << pred_h == (4 << bhl));

  for (which_mv = 0; which_mv < 1 + use_second_ref; ++which_mv) {
    // source
    const uint8_t * const base_pre = arg->pre[which_mv][plane];
    const int pre_stride = arg->pre_stride[which_mv][plane];
    const uint8_t *const pre = base_pre +
        scaled_buffer_offset(x, y, pre_stride, &xd->scale_factor[which_mv]);
    struct scale_factors * const scale = &xd->scale_factor[which_mv];

    // dest
    uint8_t *const dst = arg->dst[plane] + arg->dst_stride[plane] * y + x;

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
    const MV res_mv = clamp_mv_to_umv_border_sb(&mv, bwl, bhl,
                                                xd->plane[plane].subsampling_x,
                                                xd->plane[plane].subsampling_y,
                                                xd->mb_to_left_edge,
                                                xd->mb_to_top_edge,
                                                xd->mb_to_right_edge,
                                                xd->mb_to_bottom_edge);
    scale->set_scaled_offsets(scale, arg->y + y, arg->x + x);

#if CONFIG_MASKED_COMPOUND
    if (which_mv && xd->mode_info_context->mbmi.use_masked_compound) {
      uint8_t tmp_dst[4096];
      vp9_build_inter_predictor(pre, pre_stride,
                                tmp_dst, 64,
                                &res_mv, &xd->scale_factor[which_mv],
                                4 << pred_w, 4 << pred_h, 0,
                                &xd->subpix, MV_PRECISION_Q4);
      build_masked_compound(dst, arg->dst_stride[plane],
                            tmp_dst, 64,
                            xd->mode_info_context->mbmi.mask_index,
                            xd->mode_info_context->mbmi.sb_type,
                            (4 << pred_h), (4 << pred_w));

    } else {
#endif
    vp9_build_inter_predictor(pre, pre_stride,
                              dst, arg->dst_stride[plane],
                              &res_mv, &xd->scale_factor[which_mv],
                              4 << pred_w, 4 << pred_h, which_mv,
                              &xd->subpix, MV_PRECISION_Q4);
#if CONFIG_MASKED_COMPOUND
    }
#endif
  }
}
void vp9_build_inter_predictors_sby(MACROBLOCKD *xd,
                                    int mi_row,
                                    int mi_col,
                                    BLOCK_SIZE_TYPE bsize) {
  struct build_inter_predictors_args args = {
    xd, mi_col * MI_SIZE, mi_row * MI_SIZE,
    {xd->plane[0].dst.buf, NULL, NULL}, {xd->plane[0].dst.stride, 0, 0},
    {{xd->plane[0].pre[0].buf, NULL, NULL},
     {xd->plane[0].pre[1].buf, NULL, NULL}},
    {{xd->plane[0].pre[0].stride, 0, 0}, {xd->plane[0].pre[1].stride, 0, 0}},
  };

  foreach_predicted_block_in_plane(xd, bsize, 0, build_inter_predictors, &args);
}
void vp9_build_inter_predictors_sbuv(MACROBLOCKD *xd,
                                     int mi_row,
                                     int mi_col,
                                     BLOCK_SIZE_TYPE bsize) {
  struct build_inter_predictors_args args = {
    xd, mi_col * MI_SIZE, mi_row * MI_SIZE,
#if CONFIG_ALPHA
    {NULL, xd->plane[1].dst.buf, xd->plane[2].dst.buf,
     xd->plane[3].dst.buf},
    {0, xd->plane[1].dst.stride, xd->plane[1].dst.stride,
     xd->plane[3].dst.stride},
    {{NULL, xd->plane[1].pre[0].buf, xd->plane[2].pre[0].buf,
      xd->plane[3].pre[0].buf},
     {NULL, xd->plane[1].pre[1].buf, xd->plane[2].pre[1].buf,
      xd->plane[3].pre[1].buf}},
    {{0, xd->plane[1].pre[0].stride, xd->plane[1].pre[0].stride,
      xd->plane[3].pre[0].stride},
     {0, xd->plane[1].pre[1].stride, xd->plane[1].pre[1].stride,
      xd->plane[3].pre[1].stride}},
#else
    {NULL, xd->plane[1].dst.buf, xd->plane[2].dst.buf},
    {0, xd->plane[1].dst.stride, xd->plane[1].dst.stride},
    {{NULL, xd->plane[1].pre[0].buf, xd->plane[2].pre[0].buf},
     {NULL, xd->plane[1].pre[1].buf, xd->plane[2].pre[1].buf}},
    {{0, xd->plane[1].pre[0].stride, xd->plane[1].pre[0].stride},
     {0, xd->plane[1].pre[1].stride, xd->plane[1].pre[1].stride}},
#endif
  };
  foreach_predicted_block_uv(xd, bsize, build_inter_predictors, &args);
}
void vp9_build_inter_predictors_sb(MACROBLOCKD *xd,
                                   int mi_row, int mi_col,
                                   BLOCK_SIZE_TYPE bsize) {

#if CONFIG_INTERINTRA
  uint8_t *const y = xd->plane[0].dst.buf;
  uint8_t *const u = xd->plane[1].dst.buf;
  uint8_t *const v = xd->plane[2].dst.buf;
  const int y_stride = xd->plane[0].dst.stride;
  const int uv_stride = xd->plane[1].dst.stride;
#endif
  vp9_build_inter_predictors_sby(xd, mi_row, mi_col, bsize);
  vp9_build_inter_predictors_sbuv(xd, mi_row, mi_col, bsize);
#if CONFIG_INTERINTRA
  if (xd->mode_info_context->mbmi.ref_frame[1] == INTRA_FRAME
      && is_interintra_allowed(xd->mode_info_context->mbmi.sb_type)) {
    xd->right_available = 0;
    vp9_build_interintra_predictors(xd, y, u, v,
                                    y_stride, uv_stride, bsize);
  }
#endif
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

