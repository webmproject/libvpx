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

#include "./vpx_scale_rtcd.h"
#include "./vpx_config.h"

#include "vpx/vpx_integer.h"

#include "vp9/common/vp9_blockd.h"
#include "vp9/common/vp9_filter.h"
#include "vp9/common/vp9_reconinter.h"
#include "vp9/common/vp9_reconintra.h"
#if CONFIG_GLOBAL_MOTION
#include "vp9/common/vp9_motion_model.h"
#endif  // CONFIG_GLOBAL_MOTION

static void build_mc_border(const uint8_t *src, int src_stride,
                            uint8_t *dst, int dst_stride,
                            int x, int y, int b_w, int b_h, int w, int h) {
  // Get a pointer to the start of the real data for this row.
  const uint8_t *ref_row = src - x - y * src_stride;

  if (y >= h)
    ref_row += (h - 1) * src_stride;
  else if (y > 0)
    ref_row += y * src_stride;

  do {
    int right = 0, copy;
    int left = x < 0 ? -x : 0;

    if (left > b_w)
      left = b_w;

    if (x + b_w > w)
      right = x + b_w - w;

    if (right > b_w)
      right = b_w;

    copy = b_w - left - right;

    if (left)
      memset(dst, ref_row[0], left);

    if (copy)
      memcpy(dst + left, ref_row + x + left, copy);

    if (right)
      memset(dst + left + copy, ref_row[w - 1], right);

    dst += dst_stride;
    ++y;

    if (y > 0 && y < h)
      ref_row += src_stride;
  } while (--b_h);
}

#if CONFIG_VP9_HIGHBITDEPTH
static void high_build_mc_border(const uint8_t *src8, int src_stride,
                                 uint16_t *dst, int dst_stride,
                                 int x, int y, int b_w, int b_h,
                                 int w, int h) {
  // Get a pointer to the start of the real data for this row.
  const uint16_t *src = CONVERT_TO_SHORTPTR(src8);
  const uint16_t *ref_row = src - x - y * src_stride;

  if (y >= h)
    ref_row += (h - 1) * src_stride;
  else if (y > 0)
    ref_row += y * src_stride;

  do {
    int right = 0, copy;
    int left = x < 0 ? -x : 0;

    if (left > b_w)
      left = b_w;

    if (x + b_w > w)
      right = x + b_w - w;

    if (right > b_w)
      right = b_w;

    copy = b_w - left - right;

    if (left)
      vpx_memset16(dst, ref_row[0], left);

    if (copy)
      memcpy(dst + left, ref_row + x + left, copy * sizeof(uint16_t));

    if (right)
      vpx_memset16(dst + left + copy, ref_row[w - 1], right);

    dst += dst_stride;
    ++y;

    if (y > 0 && y < h)
      ref_row += src_stride;
  } while (--b_h);
}
#endif  // CONFIG_VP9_HIGHBITDEPTH

static void inter_predictor(const uint8_t *src, int src_stride,
                            uint8_t *dst, int dst_stride,
                            const int subpel_x,
                            const int subpel_y,
                            const struct scale_factors *sf,
                            int w, int h, int ref,
                            const InterpKernel *kernel,
                            int xs, int ys) {
#if CONFIG_EXT_CODING_UNIT_SIZE
  int i, j;
  for (j = 0; j < h; j += 64) {
    int y = subpel_y + j * ys;
    int frac_y = y & SUBPEL_MASK;
    int floor_y = y >> SUBPEL_BITS;
    for (i = 0; i < w; i += 64) {
      int x = subpel_x + i * xs;
      int frac_x = x & SUBPEL_MASK;
      int floor_x = x >> SUBPEL_BITS;
      sf->predict[frac_x != 0][frac_y != 0][ref](
          src + floor_y * src_stride + floor_x, src_stride,
          dst + j * dst_stride + i, dst_stride,
          kernel[frac_x], xs, kernel[frac_y], ys,
          w < 64 ? w : 64, h < 64 ? h : 64);
    }
  }
#else
  sf->predict[subpel_x != 0][subpel_y != 0][ref](
      src, src_stride, dst, dst_stride,
      kernel[subpel_x], xs, kernel[subpel_y], ys, w, h);
#endif
}

void vp9_build_inter_predictor(const uint8_t *src, int src_stride,
                               uint8_t *dst, int dst_stride,
                               const MV *src_mv,
                               const struct scale_factors *sf,
                               int w, int h, int ref,
                               const InterpKernel *kernel,
                               enum mv_precision precision,
                               int x, int y) {
  const int is_q4 = precision == MV_PRECISION_Q4;
  const MV mv_q4 = { is_q4 ? src_mv->row : src_mv->row * 2,
                     is_q4 ? src_mv->col : src_mv->col * 2 };
  MV32 mv = vp9_scale_mv(&mv_q4, x, y, sf);
  const int subpel_x = mv.col & SUBPEL_MASK;
  const int subpel_y = mv.row & SUBPEL_MASK;

  src += (mv.row >> SUBPEL_BITS) * src_stride + (mv.col >> SUBPEL_BITS);

  inter_predictor(src, src_stride, dst, dst_stride, subpel_x, subpel_y,
                  sf, w, h, ref, kernel, sf->x_step_q4, sf->y_step_q4);
}

#if CONFIG_VP9_HIGHBITDEPTH
static void highbd_inter_predictor(const uint8_t *src, int src_stride,
                                   uint8_t *dst, int dst_stride,
                                   const int subpel_x,
                                   const int subpel_y,
                                   const struct scale_factors *sf,
                                   int w, int h, int ref,
                                   const InterpKernel *kernel,
                                   int xs, int ys, int bd) {
#if CONFIG_EXT_CODING_UNIT_SIZE
  int i, j;
  for (j = 0; j < h; j += 64) {
    int y = subpel_y + j * ys;
    int frac_y = y & SUBPEL_MASK;
    int floor_y = y >> SUBPEL_BITS;
    for (i = 0; i < w; i += 64) {
      int x = subpel_x + i * xs;
      int frac_x = x & SUBPEL_MASK;
      int floor_x = x >> SUBPEL_BITS;
      sf->highbd_predict[frac_x != 0][frac_y != 0][ref](
          src + floor_y * src_stride + floor_x, src_stride,
          dst + j * dst_stride + i, dst_stride,
          kernel[frac_x], xs, kernel[frac_y], ys,
          w < 64 ? w : 64, h < 64 ? h : 64, bd);
    }
  }
#else
  sf->highbd_predict[subpel_x != 0][subpel_y != 0][ref](
      src, src_stride, dst, dst_stride,
      kernel[subpel_x], xs, kernel[subpel_y], ys, w, h, bd);
#endif
}

void vp9_highbd_build_inter_predictor(const uint8_t *src, int src_stride,
                                      uint8_t *dst, int dst_stride,
                                      const MV *src_mv,
                                      const struct scale_factors *sf,
                                      int w, int h, int ref,
                                      const InterpKernel *kernel,
                                      enum mv_precision precision,
                                      int x, int y, int bd) {
  const int is_q4 = precision == MV_PRECISION_Q4;
  const MV mv_q4 = { is_q4 ? src_mv->row : src_mv->row * 2,
                     is_q4 ? src_mv->col : src_mv->col * 2 };
  MV32 mv = vp9_scale_mv(&mv_q4, x, y, sf);
  const int subpel_x = mv.col & SUBPEL_MASK;
  const int subpel_y = mv.row & SUBPEL_MASK;

  src += (mv.row >> SUBPEL_BITS) * src_stride + (mv.col >> SUBPEL_BITS);

  highbd_inter_predictor(src, src_stride, dst, dst_stride, subpel_x, subpel_y,
                         sf, w, h, ref, kernel, sf->x_step_q4, sf->y_step_q4,
                         bd);
}
#endif  // CONFIG_VP9_HIGHBITDEPTH

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

static INLINE int round_mv_comp_q2(int value) {
  return (value < 0 ? value - 1 : value + 1) / 2;
}

static MV mi_mv_pred_q2(const MODE_INFO *mi, int idx, int block0, int block1) {
  MV res = { round_mv_comp_q2(mi->bmi[block0].as_mv[idx].as_mv.row +
                              mi->bmi[block1].as_mv[idx].as_mv.row),
             round_mv_comp_q2(mi->bmi[block0].as_mv[idx].as_mv.col +
                              mi->bmi[block1].as_mv[idx].as_mv.col) };
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
    src_mv->row * (1 << (1 - ss_y)),
    src_mv->col * (1 << (1 - ss_x))
  };
  assert(ss_x <= 1);
  assert(ss_y <= 1);

  clamp_mv(&clamped_mv,
           xd->mb_to_left_edge * (1 << (1 - ss_x)) - spel_left,
           xd->mb_to_right_edge * (1 << (1 - ss_x)) + spel_right,
           xd->mb_to_top_edge * (1 << (1 - ss_y)) - spel_top,
           xd->mb_to_bottom_edge * (1 << (1 - ss_y)) + spel_bottom);

  return clamped_mv;
}

static MV average_split_mvs(const struct macroblockd_plane *pd,
                            const MODE_INFO *mi, int ref, int block) {
  const int ss_idx = ((pd->subsampling_x > 0) << 1) | (pd->subsampling_y > 0);
  MV res = {0, 0};
  switch (ss_idx) {
    case 0:
      res = mi->bmi[block].as_mv[ref].as_mv;
      break;
    case 1:
      res = mi_mv_pred_q2(mi, ref, block, block + 2);
      break;
    case 2:
      res = mi_mv_pred_q2(mi, ref, block, block + 1);
      break;
    case 3:
      res = mi_mv_pred_q4(mi, ref);
      break;
    default:
      assert(ss_idx <= 3 || ss_idx >= 0);
  }
  return res;
}

#if CONFIG_WEDGE_PARTITION
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
    return (1 << WEDGE_WEIGHT_BITS);
  else
    return smoothfn[m + SMOOTHER_LEN];
}

// [negative][transpose][reverse]
DECLARE_ALIGNED(16, static uint8_t,
                wedge_mask_obl[2][2][2][MASK_MASTER_SIZE * MASK_MASTER_SIZE]);
// [negative][transpose]
DECLARE_ALIGNED(16, static uint8_t,
                wedge_mask_str[2][2][MASK_MASTER_SIZE * MASK_MASTER_SIZE]);

void vp9_init_wedge_masks() {
  int i, j;
  const int w = MASK_MASTER_SIZE;
  const int h = MASK_MASTER_SIZE;
  const int stride = MASK_MASTER_STRIDE;
  const int a[4] = {2, 1, 2, 2};
  for (i = 0; i < h; ++i)
    for (j = 0; j < w; ++j) {
      int x = (2 * j + 1 - (a[2] * w) / 2);
      int y = (2 * i + 1 - (a[3] * h) / 2);
      int m = (a[0] * x + a[1] * y) / 2;
      wedge_mask_obl[0][0][0][i * stride + j] =
          wedge_mask_obl[0][1][0][j * stride + i] =
          wedge_mask_obl[0][0][1][i * stride + w - 1 - j] =
          wedge_mask_obl[0][1][1][(w - 1 - j) * stride + i] =
          get_masked_weight(m);
      wedge_mask_obl[1][0][0][i * stride + j] =
          wedge_mask_obl[1][1][0][j * stride + i] =
          wedge_mask_obl[1][0][1][i * stride + w - 1 - j] =
          wedge_mask_obl[1][1][1][(w - 1 - j) * stride + i] =
          (1 << WEDGE_WEIGHT_BITS) - get_masked_weight(m);
      wedge_mask_str[0][0][i * stride + j] =
          wedge_mask_str[0][1][j * stride + i] =
          get_masked_weight(x);
      wedge_mask_str[1][0][i * stride + j] =
          wedge_mask_str[1][1][j * stride + i] =
          (1 << WEDGE_WEIGHT_BITS) - get_masked_weight(x);
    }
}

static void get_wedge_mask_from_array(const int *a,
                                      int h, int w,
                                      uint8_t *mask, int stride) {
  const int woff = (a[2] * w) >> 2;
  const int hoff = (a[3] * h) >> 2;
  const int oblique = (abs(a[0]) + abs(a[1]) == 3);
  const uint8_t *master;
  int transpose, reverse, negative;
  int i;

  if (oblique) {
    negative = (a[0] < 0);
    transpose = (abs(a[0]) == 1);
    reverse = (a[0] < 0) ^ (a[1] < 0);
  } else {
    negative = (a[0] < 0 || a[1] < 0);
    transpose = (a[0] == 0);
    reverse = 0;
  }
  master = (oblique ?
            wedge_mask_obl[negative][transpose][reverse] :
            wedge_mask_str[negative][transpose]) +
      MASK_MASTER_STRIDE * (MASK_MASTER_SIZE / 2 - hoff) +
      MASK_MASTER_SIZE / 2 - woff;
  for (i = 0; i < h; ++i)
    vpx_memcpy(mask + i * stride, master + i * MASK_MASTER_STRIDE, w);
}

static const uint8_t *get_wedge_mask_inplace(const int *a,
                                             int h, int w) {
  const int woff = (a[2] * w) >> 2;
  const int hoff = (a[3] * h) >> 2;
  const int oblique = (abs(a[0]) + abs(a[1]) == 3);
  const uint8_t *master;
  int transpose, reverse, negative;
  if (oblique) {
    negative = (a[0] < 0);
    transpose = (abs(a[0]) == 1);
    reverse = (a[0] < 0) ^ (a[1] < 0);
  } else {
    negative = (a[0] < 0 || a[1] < 0);
    transpose = (a[0] == 0);
    reverse = 0;
  }
  master = (oblique ?
            wedge_mask_obl[negative][transpose][reverse] :
            wedge_mask_str[negative][transpose]) +
      MASK_MASTER_STRIDE * (MASK_MASTER_SIZE / 2 - hoff) +
      MASK_MASTER_SIZE / 2 - woff;
  return master;
}

// Equation of line: f(x, y) = a[0]*(x - a[2]*w/4) + a[1]*(y - a[3]*h/4) = 0
// The soft mask is obtained by computing f(x, y) and then calling
// get_masked_weight(f(x, y)).
static const int wedge_params_sml[1 << WEDGE_BITS_SML][4] = {
  {-1,  2, 2, 2},
  { 1, -2, 2, 2},
  {-2,  1, 2, 2},
  { 2, -1, 2, 2},
  {-2, -1, 2, 2},
  { 2,  1, 2, 2},
  {-1, -2, 2, 2},
  { 1,  2, 2, 2},
};

static const int wedge_params_med_hgtw[1 << WEDGE_BITS_MED][4] = {
  {-1,  2, 2, 2},
  { 1, -2, 2, 2},
  {-2,  1, 2, 2},
  { 2, -1, 2, 2},
  {-2, -1, 2, 2},
  { 2,  1, 2, 2},
  {-1, -2, 2, 2},
  { 1,  2, 2, 2},

  {-1,  2, 2, 1},
  { 1, -2, 2, 1},
  {-1,  2, 2, 3},
  { 1, -2, 2, 3},
  {-1, -2, 2, 1},
  { 1,  2, 2, 1},
  {-1, -2, 2, 3},
  { 1,  2, 2, 3},
};

static const int wedge_params_med_hltw[1 << WEDGE_BITS_MED][4] = {
  {-1,  2, 2, 2},
  { 1, -2, 2, 2},
  {-2,  1, 2, 2},
  { 2, -1, 2, 2},
  {-2, -1, 2, 2},
  { 2,  1, 2, 2},
  {-1, -2, 2, 2},
  { 1,  2, 2, 2},

  {-2,  1, 1, 2},
  { 2, -1, 1, 2},
  {-2,  1, 3, 2},
  { 2, -1, 3, 2},
  {-2, -1, 1, 2},
  { 2,  1, 1, 2},
  {-2, -1, 3, 2},
  { 2,  1, 3, 2},
};

static const int wedge_params_med_heqw[1 << WEDGE_BITS_MED][4] = {
  {-1,  2, 2, 2},
  { 1, -2, 2, 2},
  {-2,  1, 2, 2},
  { 2, -1, 2, 2},
  {-2, -1, 2, 2},
  { 2,  1, 2, 2},
  {-1, -2, 2, 2},
  { 1,  2, 2, 2},

  { 0, -2, 0, 1},
  { 0,  2, 0, 1},
  { 0, -2, 0, 3},
  { 0,  2, 0, 3},
  {-2,  0, 1, 0},
  { 2,  0, 1, 0},
  {-2,  0, 3, 0},
  { 2,  0, 3, 0},
};

static const int wedge_params_big_hgtw[1 << WEDGE_BITS_BIG][4] = {
  {-1,  2, 2, 2},
  { 1, -2, 2, 2},
  {-2,  1, 2, 2},
  { 2, -1, 2, 2},
  {-2, -1, 2, 2},
  { 2,  1, 2, 2},
  {-1, -2, 2, 2},
  { 1,  2, 2, 2},

  {-1,  2, 2, 1},
  { 1, -2, 2, 1},
  {-1,  2, 2, 3},
  { 1, -2, 2, 3},
  {-1, -2, 2, 1},
  { 1,  2, 2, 1},
  {-1, -2, 2, 3},
  { 1,  2, 2, 3},

  {-2,  1, 1, 2},
  { 2, -1, 1, 2},
  {-2,  1, 3, 2},
  { 2, -1, 3, 2},
  {-2, -1, 1, 2},
  { 2,  1, 1, 2},
  {-2, -1, 3, 2},
  { 2,  1, 3, 2},

  { 0, -2, 0, 1},
  { 0,  2, 0, 1},
  { 0, -2, 0, 2},
  { 0,  2, 0, 2},
  { 0, -2, 0, 3},
  { 0,  2, 0, 3},
  {-2,  0, 2, 0},
  { 2,  0, 2, 0},
};

static const int wedge_params_big_hltw[1 << WEDGE_BITS_BIG][4] = {
  {-1,  2, 2, 2},
  { 1, -2, 2, 2},
  {-2,  1, 2, 2},
  { 2, -1, 2, 2},
  {-2, -1, 2, 2},
  { 2,  1, 2, 2},
  {-1, -2, 2, 2},
  { 1,  2, 2, 2},

  {-1,  2, 2, 1},
  { 1, -2, 2, 1},
  {-1,  2, 2, 3},
  { 1, -2, 2, 3},
  {-1, -2, 2, 1},
  { 1,  2, 2, 1},
  {-1, -2, 2, 3},
  { 1,  2, 2, 3},

  {-2,  1, 1, 2},
  { 2, -1, 1, 2},
  {-2,  1, 3, 2},
  { 2, -1, 3, 2},
  {-2, -1, 1, 2},
  { 2,  1, 1, 2},
  {-2, -1, 3, 2},
  { 2,  1, 3, 2},

  { 0, -2, 0, 2},
  { 0,  2, 0, 2},
  {-2,  0, 1, 0},
  { 2,  0, 1, 0},
  {-2,  0, 2, 0},
  { 2,  0, 2, 0},
  {-2,  0, 3, 0},
  { 2,  0, 3, 0},
};

static const int wedge_params_big_heqw[1 << WEDGE_BITS_BIG][4] = {
  {-1,  2, 2, 2},
  { 1, -2, 2, 2},
  {-2,  1, 2, 2},
  { 2, -1, 2, 2},
  {-2, -1, 2, 2},
  { 2,  1, 2, 2},
  {-1, -2, 2, 2},
  { 1,  2, 2, 2},

  {-1,  2, 2, 1},
  { 1, -2, 2, 1},
  {-1,  2, 2, 3},
  { 1, -2, 2, 3},
  {-1, -2, 2, 1},
  { 1,  2, 2, 1},
  {-1, -2, 2, 3},
  { 1,  2, 2, 3},

  {-2,  1, 1, 2},
  { 2, -1, 1, 2},
  {-2,  1, 3, 2},
  { 2, -1, 3, 2},
  {-2, -1, 1, 2},
  { 2,  1, 1, 2},
  {-2, -1, 3, 2},
  { 2,  1, 3, 2},

  { 0, -2, 0, 1},
  { 0,  2, 0, 1},
  { 0, -2, 0, 3},
  { 0,  2, 0, 3},
  {-2,  0, 1, 0},
  { 2,  0, 1, 0},
  {-2,  0, 3, 0},
  { 2,  0, 3, 0},
};

static const int *get_wedge_params(int wedge_index,
                                   BLOCK_SIZE sb_type,
                                   int h, int w) {
  const int *a = NULL;
  const int wedge_bits = get_wedge_bits(sb_type);

  if (wedge_index == WEDGE_NONE)
    return NULL;

  if (wedge_bits == WEDGE_BITS_SML) {
    a = wedge_params_sml[wedge_index];
  } else if (wedge_bits == WEDGE_BITS_MED) {
    if (h > w)
      a = wedge_params_med_hgtw[wedge_index];
    else if (h < w)
      a = wedge_params_med_hltw[wedge_index];
    else
      a = wedge_params_med_heqw[wedge_index];
  } else if (wedge_bits == WEDGE_BITS_BIG) {
    if (h > w)
      a = wedge_params_big_hgtw[wedge_index];
    else if (h < w)
      a = wedge_params_big_hltw[wedge_index];
    else
      a = wedge_params_big_heqw[wedge_index];
  } else {
    assert(0);
  }
  return a;
}

const uint8_t *vp9_get_soft_mask(int wedge_index,
                                 BLOCK_SIZE sb_type,
                                 int h, int w) {
  const int *a = get_wedge_params(wedge_index, sb_type, h, w);
  if (a) {
    return get_wedge_mask_inplace(a, h, w);
  } else {
    return NULL;
  }
}

// To be deprecated
void vp9_generate_soft_mask(int wedge_index,
                            BLOCK_SIZE sb_type,
                            int h, int w,
                            uint8_t *mask, int stride) {
  const int *a = get_wedge_params(wedge_index, sb_type, h, w);
  if (a) {
    get_wedge_mask_from_array(a, h, w, mask, stride);
  }
}

// To be deprecated
void vp9_generate_hard_mask(int wedge_index, BLOCK_SIZE sb_type,
                            int h, int w, uint8_t *mask, int stride) {
  int i, j;
  const int *a = get_wedge_params(wedge_index, sb_type, h, w);
  if (a) {
    get_wedge_mask_from_array(a, h, w, mask, stride);
    for (i = 0; i < h; ++i)
      for (j = 0; j < w; ++j) {
        mask[i * stride + j] = mask[i * stride + j] > 0;
      }
  }
}

static void build_masked_compound(uint8_t *dst, int dst_stride,
                                  uint8_t *dst2, int dst2_stride,
                                  int wedge_index, BLOCK_SIZE sb_type,
                                  int h, int w) {
  int i, j;
  const uint8_t *mask = vp9_get_soft_mask(wedge_index, sb_type, h, w);
  for (i = 0; i < h; ++i)
    for (j = 0; j < w; ++j) {
      int m = mask[i * MASK_MASTER_STRIDE + j];
      dst[i * dst_stride + j] = (dst[i * dst_stride + j] * m +
                                 dst2[i * dst2_stride + j] *
                                 ((1 << WEDGE_WEIGHT_BITS) - m) +
                                 (1 << (WEDGE_WEIGHT_BITS - 1))) >>
                                 WEDGE_WEIGHT_BITS;
    }
}

#if CONFIG_VP9_HIGHBITDEPTH
static void build_masked_compound_highbd(uint8_t *dst_8, int dst_stride,
                                         uint8_t *dst2_8, int dst2_stride,
                                         int wedge_index, BLOCK_SIZE sb_type,
                                         int h, int w) {
  int i, j;
  const uint8_t *mask = vp9_get_soft_mask(wedge_index, sb_type, h, w);
  uint16_t *dst = CONVERT_TO_SHORTPTR(dst_8);
  uint16_t *dst2 = CONVERT_TO_SHORTPTR(dst2_8);
  for (i = 0; i < h; ++i)
    for (j = 0; j < w; ++j) {
      int m = mask[i * MASK_MASTER_STRIDE + j];
      dst[i * dst_stride + j] = (dst[i * dst_stride + j] * m +
                                 dst2[i * dst2_stride + j] *
                                 ((1 << WEDGE_WEIGHT_BITS) - m) +
                                 (1 << (WEDGE_WEIGHT_BITS - 1))) >>
                                 WEDGE_WEIGHT_BITS;
    }
}
#endif  // CONFIG_VP9_HIGHBITDEPTH

#if CONFIG_SUPERTX
const uint8_t *get_soft_mask_extend(int wedge_index, int plane,
                                    BLOCK_SIZE sb_type,
                                    int h, int w,
                                    int wedge_offset_y,
                                    int wedge_offset_x) {
  int subh = (plane ? 2 : 4) << b_height_log2_lookup[sb_type];
  int subw = (plane ? 2 : 4) << b_width_log2_lookup[sb_type];
  const int *a = get_wedge_params(wedge_index, sb_type, subh, subw);
  (void) h;
  (void) w;
  if (a) {
    const uint8_t *mask = get_wedge_mask_inplace(a, subh, subw);
    mask -= (wedge_offset_x + wedge_offset_y * MASK_MASTER_STRIDE);
    return mask;
  } else {
    return NULL;
  }
}

// To be deprecated
static void generate_soft_mask_extend(int wedge_index, int plane,
                                      BLOCK_SIZE sb_type, int h, int w,
                                      int wedge_offset_y,
                                      int wedge_offset_x,
                                      uint8_t *mask, int stride) {
  int i, j;
  int subh = (plane ? 2 : 4) << b_height_log2_lookup[sb_type];
  int subw = (plane ? 2 : 4) << b_width_log2_lookup[sb_type];
  const int *a = get_wedge_params(wedge_index, sb_type, subh, subw);
  if (!a) return;
  for (i = 0; i < h; ++i)
    for (j = 0; j < w; ++j) {
      int x = (2 * j + 1 - (a[2] * subw) / 2 - 2 * wedge_offset_x);
      int y = (2 * i + 1 - (a[3] * subh) / 2 - 2 * wedge_offset_y);
      int m = (a[0] * x + a[1] * y) / 2;
      mask[i * stride + j] = get_masked_weight(m);
    }
}

static void build_masked_compound_extend(uint8_t *dst, int dst_stride,
                                         uint8_t *dst2, int dst2_stride,
                                         int plane,
                                         int wedge_index, BLOCK_SIZE sb_type,
                                         int wedge_offset_y, int wedge_offset_x,
                                         int h, int w) {
  int i, j;
  const uint8_t *mask = get_soft_mask_extend(
     wedge_index, plane, sb_type, h, w, wedge_offset_y, wedge_offset_x);
  for (i = 0; i < h; ++i)
    for (j = 0; j < w; ++j) {
      int m = mask[i * MASK_MASTER_STRIDE + j];
      dst[i * dst_stride + j] = (dst[i * dst_stride + j] * m +
                                 dst2[i * dst2_stride + j] *
                                 ((1 << WEDGE_WEIGHT_BITS) - m) +
                                 (1 << (WEDGE_WEIGHT_BITS - 1))) >>
                                 WEDGE_WEIGHT_BITS;
    }
}

#if CONFIG_VP9_HIGHBITDEPTH
static void build_masked_compound_extend_highbd(
    uint8_t *dst_8, int dst_stride,
    uint8_t *dst2_8, int dst2_stride, int plane,
    int wedge_index, BLOCK_SIZE sb_type,
    int wedge_offset_y, int wedge_offset_x,
    int h, int w) {
  int i, j;
  const uint8_t *mask = get_soft_mask_extend(
      wedge_index, plane, sb_type, h, w, wedge_offset_y, wedge_offset_x);
  uint16_t *dst = CONVERT_TO_SHORTPTR(dst_8);
  uint16_t *dst2 = CONVERT_TO_SHORTPTR(dst2_8);
  for (i = 0; i < h; ++i)
    for (j = 0; j < w; ++j) {
      int m = mask[i * MASK_MASTER_STRIDE + j];
      dst[i * dst_stride + j] = (dst[i * dst_stride + j] * m +
                                 dst2[i * dst2_stride + j] *
                                 ((1 << WEDGE_WEIGHT_BITS) - m) +
                                 (1 << (WEDGE_WEIGHT_BITS - 1))) >>
                                 WEDGE_WEIGHT_BITS;
    }
}
#endif  // CONFIG_VP9_HIGHBITDEPTH
#endif  // CONFIG_SUPERTX
#endif  // CONFIG_WEDGE_PARTITION

static void build_inter_predictors(MACROBLOCKD *xd, int plane, int block,
                                   int bw, int bh,
                                   int x, int y, int w, int h,
#if CONFIG_SUPERTX && CONFIG_WEDGE_PARTITION
                                   int wedge_offset_x, int wedge_offset_y,
#endif  // CONFIG_SUPERTX && CONFIG_WEDGE_PARTITION
                                   int mi_x, int mi_y) {
  struct macroblockd_plane *const pd = &xd->plane[plane];
  const MODE_INFO *mi = xd->mi[0].src_mi;
  const int is_compound = has_second_ref(&mi->mbmi);
#if CONFIG_INTRABC
  const int is_intrabc = is_intrabc_mode(mi->mbmi.mode);
#endif  // CONFIG_INTRABC
  const InterpKernel *kernel = vp9_get_interp_kernel(mi->mbmi.interp_filter);
  int ref;
#if CONFIG_GLOBAL_MOTION
  Global_Motion_Params *gm[2];
  int is_global;
  gm[0] = &xd->global_motion[mi->mbmi.ref_frame[0]][0];
  if (is_compound)
    gm[1] = &xd->global_motion[mi->mbmi.ref_frame[1]][0];
#endif  // CONFIG_GLOBAL_MOTION
#if CONFIG_INTRABC
  assert(!is_intrabc || mi->mbmi.interp_filter == BILINEAR);
#endif  // CONFIG_INTRABC

  for (ref = 0; ref < 1 + is_compound; ++ref) {
    const struct scale_factors *const sf = &xd->block_refs[ref]->sf;
    struct buf_2d *const dst_buf = &pd->dst;
    struct buf_2d *const pre_buf =
#if CONFIG_INTRABC
        is_intrabc ? dst_buf :
#endif  // CONFIG_INTRABC
        &pd->pre[ref];
    uint8_t *const dst = dst_buf->buf + dst_buf->stride * y + x;
    const MV mv = mi->mbmi.sb_type < BLOCK_8X8
               ? average_split_mvs(pd, mi, ref, block)
               : mi->mbmi.mv[ref].as_mv;

    // TODO(jkoleszar): This clamping is done in the incorrect place for the
    // scaling case. It needs to be done on the scaled MV, not the pre-scaling
    // MV. Note however that it performs the subsampling aware scaling so
    // that the result is always q4.
    // mv_precision precision is MV_PRECISION_Q4.
    const MV mv_q4 = clamp_mv_to_umv_border_sb(xd, &mv, bw, bh,
                                               pd->subsampling_x,
                                               pd->subsampling_y);

    uint8_t *pre;
    MV32 scaled_mv;
    int xs, ys, subpel_x, subpel_y;
    const int is_scaled = vp9_is_scaled(sf);

#if CONFIG_GLOBAL_MOTION
    is_global = (get_y_mode(mi, block) == ZEROMV &&
#if CONFIG_INTRABC
                 !is_intrabc &&
#endif
                 get_gmtype(gm[ref]) == GLOBAL_ROTZOOM);
#endif  // CONFIG_GLOBAL_MOTION

    if (is_scaled) {
#if CONFIG_INTRABC
      assert(!is_intrabc);
#endif  // CONFIG_INTRABC
      pre = pre_buf->buf + scaled_buffer_offset(x, y, pre_buf->stride, sf);
      scaled_mv = vp9_scale_mv(&mv_q4, mi_x + x, mi_y + y, sf);
      xs = sf->x_step_q4;
      ys = sf->y_step_q4;
    } else {
      pre = pre_buf->buf + (y * pre_buf->stride + x);
      scaled_mv.row = mv_q4.row;
      scaled_mv.col = mv_q4.col;
      xs = ys = 16;
    }
    subpel_x = scaled_mv.col & SUBPEL_MASK;
    subpel_y = scaled_mv.row & SUBPEL_MASK;
    pre += (scaled_mv.row >> SUBPEL_BITS) * pre_buf->stride
        + (scaled_mv.col >> SUBPEL_BITS);

#if CONFIG_WEDGE_PARTITION
    if (ref && get_wedge_bits(mi->mbmi.sb_type)
        && mi->mbmi.use_wedge_interinter) {
#if CONFIG_VP9_HIGHBITDEPTH
      uint8_t tmp_dst_[2 * CODING_UNIT_SIZE * CODING_UNIT_SIZE];
      uint8_t *tmp_dst =
          (xd->cur_buf->flags & YV12_FLAG_HIGHBITDEPTH) ?
          CONVERT_TO_BYTEPTR(tmp_dst_) : tmp_dst_;
#else
      uint8_t tmp_dst[CODING_UNIT_SIZE * CODING_UNIT_SIZE];
#endif
#if CONFIG_GLOBAL_MOTION
      if (is_global) {
        vp9_warp_plane(gm[ref], pre_buf->buf0,
                       pre_buf->width, pre_buf->height, pre_buf->stride,
                       tmp_dst, (mi_x >> pd->subsampling_x) + x,
                       (mi_y >> pd->subsampling_y) + y, w, h, CODING_UNIT_SIZE,
                       pd->subsampling_x, pd->subsampling_y, xs, ys);
      } else {
#endif  // CONFIG_GLOBAL_MOTION
#if CONFIG_VP9_HIGHBITDEPTH
        if (xd->cur_buf->flags & YV12_FLAG_HIGHBITDEPTH) {
          highbd_inter_predictor(pre, pre_buf->stride, tmp_dst,
                                 CODING_UNIT_SIZE, subpel_x, subpel_y, sf, w, h,
                                 0, kernel, xs, ys, xd->bd);
        } else {
          inter_predictor(pre, pre_buf->stride, tmp_dst, CODING_UNIT_SIZE,
                          subpel_x, subpel_y, sf, w, h, 0, kernel, xs, ys);
        }
#else
        inter_predictor(pre, pre_buf->stride, tmp_dst, CODING_UNIT_SIZE,
                        subpel_x, subpel_y, sf, w, h, 0, kernel, xs, ys);
#endif  // CONFIG_VP9_HIGHBITDEPTH
#if CONFIG_GLOBAL_MOTION
      }
#endif  // CONFIG_GLOBAL_MOTION
#if CONFIG_SUPERTX
#if CONFIG_VP9_HIGHBITDEPTH
      if (xd->cur_buf->flags & YV12_FLAG_HIGHBITDEPTH) {
        build_masked_compound_extend_highbd(
            dst, dst_buf->stride, tmp_dst, CODING_UNIT_SIZE, plane,
            mi->mbmi.interinter_wedge_index,
            mi->mbmi.sb_type,
            wedge_offset_y, wedge_offset_x, h, w);
      } else {
        build_masked_compound_extend(
            dst, dst_buf->stride, tmp_dst, CODING_UNIT_SIZE, plane,
            mi->mbmi.interinter_wedge_index,
            mi->mbmi.sb_type,
            wedge_offset_y, wedge_offset_x, h, w);
      }
#else
      build_masked_compound_extend(dst, dst_buf->stride, tmp_dst,
                                   CODING_UNIT_SIZE, plane,
                                   mi->mbmi.interinter_wedge_index,
                                   mi->mbmi.sb_type,
                                   wedge_offset_y, wedge_offset_x, h, w);
#endif  // CONFIG_VP9_HIGHBITDEPTH
#else   // CONFIG_SUPERTX
#if CONFIG_VP9_HIGHBITDEPTH
      if (xd->cur_buf->flags & YV12_FLAG_HIGHBITDEPTH)
        build_masked_compound_highbd(dst, dst_buf->stride, tmp_dst,
                                     CODING_UNIT_SIZE,
                                     mi->mbmi.interinter_wedge_index,
                                     mi->mbmi.sb_type, h, w);
      else
#endif  // CONFIG_VP9_HIGHBITDEPTH
        build_masked_compound(dst, dst_buf->stride, tmp_dst, CODING_UNIT_SIZE,
                              mi->mbmi.interinter_wedge_index, mi->mbmi.sb_type,
                              h, w);
#endif  // CONFIG_SUPERTX
    } else {
#if CONFIG_GLOBAL_MOTION
      if (is_global) {
        vp9_warp_plane(gm[ref], pre_buf->buf0,
                       pre_buf->width, pre_buf->height, pre_buf->stride, dst,
                       (mi_x >> pd->subsampling_x) + x,
                       (mi_y >> pd->subsampling_y) + y, w, h, dst_buf->stride,
                       pd->subsampling_x, pd->subsampling_y, xs, ys);
      } else {
#endif  // CONFIG_GLOBAL_MOTION
#if CONFIG_VP9_HIGHBITDEPTH
        if (xd->cur_buf->flags & YV12_FLAG_HIGHBITDEPTH)
          highbd_inter_predictor(pre, pre_buf->stride, dst, dst_buf->stride,
                                 subpel_x, subpel_y, sf, w, h, ref, kernel,
                                 xs, ys, xd->bd);
        else
#endif  // CONFIG_VP9_HIGHBITDEPTH
          inter_predictor(pre, pre_buf->stride, dst, dst_buf->stride,
                          subpel_x, subpel_y, sf, w, h, ref, kernel, xs, ys);
#if CONFIG_GLOBAL_MOTION
      }
#endif  // CONFIG_GLOBAL_MOTION
    }

#else  // CONFIG_WEDGE_PARTITION

#if CONFIG_GLOBAL_MOTION
    if (is_global) {
      vp9_warp_plane(gm[ref], pre_buf->buf0,
                     pre_buf->width, pre_buf->height, pre_buf->stride, dst,
                     (mi_x >> pd->subsampling_x) + x,
                     (mi_y >> pd->subsampling_y) + y, w, h, dst_buf->stride,
                     pd->subsampling_x, pd->subsampling_y, xs, ys);
    } else {
#endif  // CONFIG_GLOBAL_MOTION
#if CONFIG_VP9_HIGHBITDEPTH
      if (xd->cur_buf->flags & YV12_FLAG_HIGHBITDEPTH) {
        highbd_inter_predictor(pre, pre_buf->stride, dst, dst_buf->stride,
                               subpel_x, subpel_y, sf, w, h, ref, kernel,
                               xs, ys, xd->bd);
      } else {
        inter_predictor(pre, pre_buf->stride, dst, dst_buf->stride,
                        subpel_x, subpel_y, sf, w, h, ref, kernel, xs, ys);
      }
#else
      inter_predictor(pre, pre_buf->stride, dst, dst_buf->stride,
                      subpel_x, subpel_y, sf, w, h, ref, kernel, xs, ys);
#endif  // CONFIG_VP9_HIGHBITDEPTH
#if CONFIG_GLOBAL_MOTION
    }
#endif  // CONFIG_GLOBAL_MOTION
#endif  // CONFIG_WEDGE_PARTITION
  }
}

static void build_inter_predictors_for_planes(MACROBLOCKD *xd, BLOCK_SIZE bsize,
                                              int mi_row, int mi_col,
                                              int plane_from, int plane_to) {
  int plane;
  const int mi_x = mi_col * MI_SIZE;
  const int mi_y = mi_row * MI_SIZE;
  for (plane = plane_from; plane <= plane_to; ++plane) {
    const BLOCK_SIZE plane_bsize = get_plane_block_size(bsize,
                                                        &xd->plane[plane]);
    const int num_4x4_w = num_4x4_blocks_wide_lookup[plane_bsize];
    const int num_4x4_h = num_4x4_blocks_high_lookup[plane_bsize];
    const int bw = 4 * num_4x4_w;
    const int bh = 4 * num_4x4_h;

    if (xd->mi[0].src_mi->mbmi.sb_type < BLOCK_8X8) {
      int i = 0, x, y;
      assert(bsize == BLOCK_8X8);
      for (y = 0; y < num_4x4_h; ++y)
        for (x = 0; x < num_4x4_w; ++x)
           build_inter_predictors(xd, plane, i++, bw, bh,
                                  4 * x, 4 * y, 4, 4,
#if CONFIG_SUPERTX && CONFIG_WEDGE_PARTITION
                                  0, 0,
#endif
                                  mi_x, mi_y);
    } else {
      build_inter_predictors(xd, plane, 0, bw, bh,
                             0, 0, bw, bh,
#if CONFIG_SUPERTX && CONFIG_WEDGE_PARTITION
                             0, 0,
#endif
                             mi_x, mi_y);
    }
  }
}

void vp9_build_inter_predictors_sby(MACROBLOCKD *xd, int mi_row, int mi_col,
                                    BLOCK_SIZE bsize) {
  build_inter_predictors_for_planes(xd, bsize, mi_row, mi_col, 0, 0);
#if CONFIG_INTERINTRA
  if (xd->mi[0].src_mi->mbmi.ref_frame[1] == INTRA_FRAME &&
#if CONFIG_INTRABC
      xd->mi[0].src_mi->mbmi.ref_frame[0] != INTRA_FRAME &&
#endif  // CONFIG_INTRABC
      is_interintra_allowed(xd->mi[0].src_mi->mbmi.sb_type))
    vp9_build_interintra_predictors_sby(xd, xd->plane[0].dst.buf,
                                        xd->plane[0].dst.stride, bsize);
#endif  // CONFIG_INTERINTRA
}

void vp9_build_inter_predictors_sbuv(MACROBLOCKD *xd, int mi_row, int mi_col,
                                     BLOCK_SIZE bsize) {
  build_inter_predictors_for_planes(xd, bsize, mi_row, mi_col, 1,
                                    MAX_MB_PLANE - 1);
#if CONFIG_INTERINTRA
  if (xd->mi[0].src_mi->mbmi.ref_frame[1] == INTRA_FRAME &&
#if CONFIG_INTRABC
      xd->mi[0].src_mi->mbmi.ref_frame[0] != INTRA_FRAME &&
#endif  // CONFIG_INTRABC
      is_interintra_allowed(xd->mi[0].src_mi->mbmi.sb_type))
    vp9_build_interintra_predictors_sbuv(xd, xd->plane[1].dst.buf,
                                         xd->plane[2].dst.buf,
                                         xd->plane[1].dst.stride,
                                         xd->plane[2].dst.stride, bsize);
#endif  // CONFIG_INTERINTRA
}

void vp9_build_inter_predictors_sb(MACROBLOCKD *xd, int mi_row, int mi_col,
                                   BLOCK_SIZE bsize) {
  build_inter_predictors_for_planes(xd, bsize, mi_row, mi_col, 0,
                                    MAX_MB_PLANE - 1);
#if CONFIG_INTERINTRA
  if (xd->mi[0].src_mi->mbmi.ref_frame[1] == INTRA_FRAME &&
#if CONFIG_INTRABC
      xd->mi[0].src_mi->mbmi.ref_frame[0] != INTRA_FRAME &&
#endif  // CONFIG_INTRABC
      is_interintra_allowed(xd->mi[0].src_mi->mbmi.sb_type))
    vp9_build_interintra_predictors(xd, xd->plane[0].dst.buf,
                                    xd->plane[1].dst.buf, xd->plane[2].dst.buf,
                                    xd->plane[0].dst.stride,
                                    xd->plane[1].dst.stride,
                                    xd->plane[2].dst.stride, bsize);
#endif  // CONFIG_INTERINTRA
}

#if CONFIG_SUPERTX

static const uint8_t mask_8[8] = {
  64, 64, 62, 52, 12,  2,  0,  0
};

static const uint8_t mask_16[16] = {
  63, 62, 60, 58, 55, 50, 43, 36, 28, 21, 14, 9, 6, 4, 2, 1
};

static const uint8_t mask_32[32] = {
  64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 63, 61, 57, 52, 45, 36,
  28, 19, 12,  7,  3,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0
};

#if CONFIG_TX64X64
static const uint8_t mask_64[64] = {
  64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64,
  64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 63, 61, 57, 52, 45, 36,
  28, 19, 12,  7,  3,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  0,   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0
};
#endif

static const uint8_t mask_8_uv[8] = {
  64, 64, 62, 52,  12,  2,  0,  0
};

static const uint8_t mask_16_uv[16] = {
  64, 64, 64, 64, 61, 53, 45, 36, 28, 19, 11, 3, 0,  0,  0,  0
};

static const uint8_t mask_32_uv[32] = {
  64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 60, 54, 46, 36,
  28, 18, 10,  4,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0
};

#if CONFIG_TX64X64
static const uint8_t mask_64_uv[64] = {
  64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64,
  64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 60, 54, 46, 36,
  28, 18, 10,  4,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  0,   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0
};
#endif

static void generate_1dmask(int length, uint8_t *mask, int plane) {
  switch (length) {
    case 8:
      vpx_memcpy(mask, plane ? mask_8_uv : mask_8, length);
      break;
    case 16:
      vpx_memcpy(mask, plane ? mask_16_uv : mask_16, length);
      break;
    case 32:
      vpx_memcpy(mask, plane ? mask_32_uv : mask_32, length);
      break;
#if CONFIG_TX64X64
    case 64:
      vpx_memcpy(mask, plane ? mask_64_uv : mask_64, length);
      break;
#endif
    default:
      assert(0);
  }
}


void vp9_build_masked_inter_predictor_complex(
    MACROBLOCKD *xd,
    uint8_t *dst, int dst_stride, uint8_t *dst2, int dst2_stride,
    const struct macroblockd_plane *pd, int mi_row, int mi_col,
    int mi_row_ori, int mi_col_ori, BLOCK_SIZE bsize, BLOCK_SIZE top_bsize,
    PARTITION_TYPE partition, int plane) {
  int i, j;
  uint8_t mask[MAXTXLEN];
  int top_w = 4 << b_width_log2_lookup[top_bsize],
      top_h = 4 << b_height_log2_lookup[top_bsize];
  int w = 4 << b_width_log2_lookup[bsize], h = 4 << b_height_log2_lookup[bsize];
  int w_offset = (mi_col - mi_col_ori) << 3,
      h_offset = (mi_row - mi_row_ori) << 3;

#if CONFIG_VP9_HIGHBITDEPTH
  uint16_t *dst16= CONVERT_TO_SHORTPTR(dst);
  uint16_t *dst216 = CONVERT_TO_SHORTPTR(dst2);
  int b_hdb = (xd->cur_buf->flags & YV12_FLAG_HIGHBITDEPTH) ? 1 : 0;
#endif  // CONFIG_VP9_HIGHBITDEPTH

  top_w >>= pd->subsampling_x;
  top_h >>= pd->subsampling_y;
  w >>= pd->subsampling_x;
  h >>= pd->subsampling_y;
  w_offset >>= pd->subsampling_x;
  h_offset >>= pd->subsampling_y;

  switch (partition) {
    case PARTITION_HORZ:
    {
#if CONFIG_VP9_HIGHBITDEPTH
      if (b_hdb) {
        uint16_t *dst_tmp = dst16 + h_offset * dst_stride;
        uint16_t *dst2_tmp = dst216 + h_offset * dst2_stride;
        generate_1dmask(h, mask + h_offset,
                        plane && xd->plane[plane].subsampling_y);

        for (i = h_offset; i < h_offset + h; i++) {
          for (j = 0; j < top_w; j++) {
            const int m = mask[i];  assert(m >= 0 && m <= 64);
            if (m == 64)
              continue;

            if (m == 0)
              dst_tmp[j] = dst2_tmp[j];
            else
              dst_tmp[j] = (dst_tmp[j] * m + dst2_tmp[j] * (64 - m) + 32) >> 6;
          }
          dst_tmp += dst_stride;
          dst2_tmp += dst2_stride;
        }

        for (; i < top_h; i ++) {
          vpx_memcpy(dst_tmp, dst2_tmp, top_w * sizeof(uint16_t));
          dst_tmp += dst_stride;
          dst2_tmp += dst2_stride;
        }
      } else {
#endif  // CONFIG_VP9_HIGHBITDEPTH
        uint8_t *dst_tmp = dst + h_offset * dst_stride;
        uint8_t *dst2_tmp = dst2 + h_offset * dst2_stride;
        generate_1dmask(h, mask + h_offset,
                        plane && xd->plane[plane].subsampling_y);

        for (i = h_offset; i < h_offset + h; i++) {
          for (j = 0; j < top_w; j++) {
            const int m = mask[i];  assert(m >= 0 && m <= 64);
            if (m == 64)
              continue;

            if (m == 0)
              dst_tmp[j] = dst2_tmp[j];
            else
              dst_tmp[j] = (dst_tmp[j] * m + dst2_tmp[j] * (64 - m) + 32) >> 6;
          }
          dst_tmp += dst_stride;
          dst2_tmp += dst2_stride;
        }

        for (; i < top_h; i ++) {
          vpx_memcpy(dst_tmp, dst2_tmp, top_w * sizeof(uint8_t));
          dst_tmp += dst_stride;
          dst2_tmp += dst2_stride;
        }
#if CONFIG_VP9_HIGHBITDEPTH
      }
#endif  // CONFIG_VP9_HIGHBITDEPTH
    }

      break;
    case PARTITION_VERT:
    {
#if CONFIG_VP9_HIGHBITDEPTH
      if (b_hdb) {
        uint16_t *dst_tmp = dst16;
        uint16_t *dst2_tmp = dst216;
        generate_1dmask(w, mask + w_offset,
                        plane && xd->plane[plane].subsampling_x);

        for (i = 0; i < top_h; i++) {
          for (j = w_offset; j < w_offset + w; j++) {
            const int m = mask[j];   assert(m >= 0 && m <= 64);
            if (m == 64)
              continue;

            if (m == 0)
              dst_tmp[j] = dst2_tmp[j];
            else
              dst_tmp[j] = (dst_tmp[j] * m + dst2_tmp[j] * (64 - m) + 32) >> 6;
          }
          vpx_memcpy(dst_tmp + j, dst2_tmp + j,
                     (top_w - w_offset - w) * sizeof(uint16_t));
          dst_tmp += dst_stride;
          dst2_tmp += dst2_stride;
        }
      } else {
#endif  // CONFIG_VP9_HIGHBITDEPTH
        uint8_t *dst_tmp = dst;
        uint8_t *dst2_tmp = dst2;
        generate_1dmask(w, mask + w_offset,
                        plane && xd->plane[plane].subsampling_x);

        for (i = 0; i < top_h; i++) {
          for (j = w_offset; j < w_offset + w; j++) {
            const int m = mask[j];   assert(m >= 0 && m <= 64);
            if (m == 64)
              continue;

            if (m == 0)
              dst_tmp[j] = dst2_tmp[j];
            else
              dst_tmp[j] = (dst_tmp[j] * m + dst2_tmp[j] * (64 - m) + 32) >> 6;
          }
            vpx_memcpy(dst_tmp + j, dst2_tmp + j,
                       (top_w - w_offset - w) * sizeof(uint8_t));
          dst_tmp += dst_stride;
          dst2_tmp += dst2_stride;
        }
#if CONFIG_VP9_HIGHBITDEPTH
      }
#endif  // CONFIG_VP9_HIGHBITDEPTH
    }
      break;
    default:
      assert(0);
  }
  (void) xd;
}

void vp9_build_inter_predictors_sb_sub8x8(MACROBLOCKD *xd,
                                          int mi_row, int mi_col,
                                          BLOCK_SIZE bsize, int block) {
  // Prediction function used in supertx:
  // Use the mv at current block (which is less than 8x8)
  // to get prediction of a block located at (mi_row, mi_col) at size of bsize
  // bsize can be larger than 8x8.
  // block (0-3): the sub8x8 location of current block
  int plane;
  const int mi_x = mi_col * MI_SIZE;
  const int mi_y = mi_row * MI_SIZE;

  // For sub8x8 uv:
  // Skip uv prediction in supertx except the first block (block = 0)
  int max_plane = block ? 1 : MAX_MB_PLANE;

  for (plane = 0; plane < max_plane; plane++) {
    const BLOCK_SIZE plane_bsize = get_plane_block_size(bsize,
                                                        &xd->plane[plane]);
    const int num_4x4_w = num_4x4_blocks_wide_lookup[plane_bsize];
    const int num_4x4_h = num_4x4_blocks_high_lookup[plane_bsize];
    const int bw = 4 * num_4x4_w;
    const int bh = 4 * num_4x4_h;

    build_inter_predictors(xd, plane, block, bw, bh,
                           0, 0, bw, bh,
#if CONFIG_SUPERTX && CONFIG_WEDGE_PARTITION
                            0, 0,
#endif
                           mi_x, mi_y);
  }
#if CONFIG_INTERINTRA
  if (xd->mi[0].src_mi->mbmi.ref_frame[1] == INTRA_FRAME &&
#if CONFIG_INTRABC
      xd->mi[0].src_mi->mbmi.ref_frame[0] != INTRA_FRAME &&
#endif  // CONFIG_INTRABC
      is_interintra_allowed(xd->mi[0].src_mi->mbmi.sb_type))
    vp9_build_interintra_predictors(xd, xd->plane[0].dst.buf,
                                    xd->plane[1].dst.buf, xd->plane[2].dst.buf,
                                    xd->plane[0].dst.stride,
                                    xd->plane[1].dst.stride,
                                    xd->plane[2].dst.stride, bsize);
#endif  // CONFIG_INTERINTRA
}

#if CONFIG_WEDGE_PARTITION
void vp9_build_inter_predictors_sb_extend(MACROBLOCKD *xd,
                                          int mi_row, int mi_col,
                                          int mi_row_ori, int mi_col_ori,
                                          BLOCK_SIZE bsize) {
  int plane;
  const int mi_x = mi_col_ori * MI_SIZE;
  const int mi_y = mi_row_ori * MI_SIZE;
  const int wedge_offset_x = (mi_col - mi_col_ori) * MI_SIZE;
  const int wedge_offset_y = (mi_row - mi_row_ori) * MI_SIZE;
  for (plane = 0; plane < MAX_MB_PLANE; ++plane) {
    const BLOCK_SIZE plane_bsize = get_plane_block_size(bsize,
                                                        &xd->plane[plane]);
    const int num_4x4_w = num_4x4_blocks_wide_lookup[plane_bsize];
    const int num_4x4_h = num_4x4_blocks_high_lookup[plane_bsize];
    const int bw = 4 * num_4x4_w;
    const int bh = 4 * num_4x4_h;

    if (xd->mi[0].src_mi->mbmi.sb_type < BLOCK_8X8) {
      int i = 0, x, y;
      assert(bsize == BLOCK_8X8);
      for (y = 0; y < num_4x4_h; ++y)
        for (x = 0; x < num_4x4_w; ++x)
           build_inter_predictors(xd, plane, i++, bw, bh, 4 * x, 4 * y, 4, 4,
                            wedge_offset_x >> (xd->plane[plane].subsampling_x),
                            wedge_offset_y >> (xd->plane[plane].subsampling_y),
                            mi_x, mi_y);
    } else {
      build_inter_predictors(xd, plane, 0, bw, bh, 0, 0, bw, bh,
                            wedge_offset_x >> (xd->plane[plane].subsampling_x),
                            wedge_offset_y >> (xd->plane[plane].subsampling_y),
                            mi_x, mi_y);
    }
  }
}
void vp9_build_inter_predictors_sb_sub8x8_extend(
    MACROBLOCKD *xd,
    int mi_row, int mi_col,
    int mi_row_ori, int mi_col_ori,
    BLOCK_SIZE bsize, int block) {
  // Sub8x8 prediction for wedge partition in supertx
  int plane;
  const int mi_x = mi_col_ori * MI_SIZE;
  const int mi_y = mi_row_ori * MI_SIZE;
  const int wedge_offset_x = (mi_col - mi_col_ori) * MI_SIZE;
  const int wedge_offset_y = (mi_row - mi_row_ori) * MI_SIZE;

  // For sub8x8 uv:
  // Skip uv prediction in supertx except the first block (block = 0)
  int max_plane = block ? 1 : MAX_MB_PLANE;

  for (plane = 0; plane < max_plane; ++plane) {
    const BLOCK_SIZE plane_bsize = get_plane_block_size(bsize,
                                                        &xd->plane[plane]);
    const int num_4x4_w = num_4x4_blocks_wide_lookup[plane_bsize];
    const int num_4x4_h = num_4x4_blocks_high_lookup[plane_bsize];
    const int bw = 4 * num_4x4_w;
    const int bh = 4 * num_4x4_h;

    build_inter_predictors(xd, plane, block, bw, bh, 0, 0, bw, bh,
                           wedge_offset_x >> (xd->plane[plane].subsampling_x),
                           wedge_offset_y >> (xd->plane[plane].subsampling_y),
                           mi_x, mi_y);
  }
}
#endif  // CONFIG_WEDGE_PARTITION
#endif  // CONFIG_SUPERTX

// TODO(jingning): This function serves as a placeholder for decoder prediction
// using on demand border extension. It should be moved to /decoder/ directory.
static void dec_build_inter_predictors(MACROBLOCKD *xd, int plane, int block,
                                       int bw, int bh,
                                       int x, int y, int w, int h,
#if CONFIG_SUPERTX && CONFIG_WEDGE_PARTITION
                                       int wedge_offset_x, int wedge_offset_y,
#endif
                                       int mi_x, int mi_y) {
  struct macroblockd_plane *const pd = &xd->plane[plane];
  const MODE_INFO *mi = xd->mi[0].src_mi;
  const int is_compound = has_second_ref(&mi->mbmi);
  const InterpKernel *kernel = vp9_get_interp_kernel(mi->mbmi.interp_filter);
  int ref;
#if CONFIG_GLOBAL_MOTION
  Global_Motion_Params *gm[2];
  int is_global;
#endif  // CONFIG_GLOBAL_MOTION
#if CONFIG_INTRABC
  const int is_intrabc = is_intrabc_mode(mi->mbmi.mode);
  struct scale_factors sf1;
#if CONFIG_VP9_HIGHBITDEPTH
  int use_highbitdepth = (xd->cur_buf->flags & YV12_FLAG_HIGHBITDEPTH);
  vp9_setup_scale_factors_for_frame(&sf1, 64, 64, 64, 64, use_highbitdepth);
#else
  vp9_setup_scale_factors_for_frame(&sf1, 64, 64, 64, 64);
#endif  // CONFIG_VP9_HIGHBITDEPTH

  assert(!is_intrabc || !is_compound);
#endif  // CONFIG_INTRABC
#if CONFIG_GLOBAL_MOTION
  gm[0] = &xd->global_motion[mi->mbmi.ref_frame[0]][0];
  if (is_compound)
    gm[1] = &xd->global_motion[mi->mbmi.ref_frame[1]][0];
#endif  // CONFIG_GLOBAL_MOTION

  for (ref = 0; ref < 1 + is_compound; ++ref) {
    struct buf_2d *const dst_buf = &pd->dst;
#if CONFIG_INTRABC
    const struct scale_factors *const sf =
        is_intrabc ? &sf1 : &xd->block_refs[ref]->sf;
    struct buf_2d *const pre_buf =
        is_intrabc ? dst_buf : &pd->pre[ref];
    const YV12_BUFFER_CONFIG *ref_buf =
        is_intrabc ? xd->cur_buf : xd->block_refs[ref]->buf;
#else
    const struct scale_factors *const sf = &xd->block_refs[ref]->sf;
    struct buf_2d *const pre_buf = &pd->pre[ref];
    const YV12_BUFFER_CONFIG *ref_buf = xd->block_refs[ref]->buf;
#endif  // CONFIG_INTRABC
    uint8_t *const dst = dst_buf->buf + dst_buf->stride * y + x;
    const MV mv = mi->mbmi.sb_type < BLOCK_8X8
        ? average_split_mvs(pd, mi, ref, block)
        : mi->mbmi.mv[ref].as_mv;

    const MV mv_q4 = clamp_mv_to_umv_border_sb(xd, &mv, bw, bh,
                                               pd->subsampling_x,
                                               pd->subsampling_y);

    MV32 scaled_mv;
    int xs, ys, x0, y0, x0_16, y0_16, frame_width, frame_height, buf_stride,
        subpel_x, subpel_y;
    uint8_t *ref_frame, *buf_ptr;
    const int is_scaled = vp9_is_scaled(sf);

#if CONFIG_GLOBAL_MOTION
    is_global = (get_y_mode(mi, block) == ZEROMV &&
#if CONFIG_INTRABC
                 !is_intrabc &&
#endif
                 get_gmtype(gm[ref]) == GLOBAL_ROTZOOM);
#endif  // CONFIG_GLOBAL_MOTION

    // Get reference frame pointer, width and height.
    if (plane == 0) {
      frame_width = ref_buf->y_crop_width;
      frame_height = ref_buf->y_crop_height;
      ref_frame = ref_buf->y_buffer;
    } else {
      frame_width = ref_buf->uv_crop_width;
      frame_height = ref_buf->uv_crop_height;
      ref_frame = plane == 1 ? ref_buf->u_buffer : ref_buf->v_buffer;
    }

    if (is_scaled) {
      // Co-ordinate of containing block to pixel precision.
      int x_start = (-xd->mb_to_left_edge >> (3 + pd->subsampling_x));
      int y_start = (-xd->mb_to_top_edge >> (3 + pd->subsampling_y));

#if CONFIG_INTRABC
      assert(!is_intrabc);
#endif  // CONFIG_INTRABC

      // Co-ordinate of the block to 1/16th pixel precision.
      x0_16 = (x_start + x) << SUBPEL_BITS;
      y0_16 = (y_start + y) << SUBPEL_BITS;

      // Co-ordinate of current block in reference frame
      // to 1/16th pixel precision.
      x0_16 = sf->scale_value_x(x0_16, sf);
      y0_16 = sf->scale_value_y(y0_16, sf);

      // Map the top left corner of the block into the reference frame.
      x0 = sf->scale_value_x(x_start + x, sf);
      y0 = sf->scale_value_y(y_start + y, sf);

      // Scale the MV and incorporate the sub-pixel offset of the block
      // in the reference frame.
      scaled_mv = vp9_scale_mv(&mv_q4, mi_x + x, mi_y + y, sf);
      xs = sf->x_step_q4;
      ys = sf->y_step_q4;
    } else {
      // Co-ordinate of containing block to pixel precision.
      x0 = (-xd->mb_to_left_edge >> (3 + pd->subsampling_x)) + x;
      y0 = (-xd->mb_to_top_edge >> (3 + pd->subsampling_y)) + y;

      // Co-ordinate of the block to 1/16th pixel precision.
      x0_16 = x0 << SUBPEL_BITS;
      y0_16 = y0 << SUBPEL_BITS;

      scaled_mv.row = mv_q4.row;
      scaled_mv.col = mv_q4.col;
      xs = ys = 16;
    }
    subpel_x = scaled_mv.col & SUBPEL_MASK;
    subpel_y = scaled_mv.row & SUBPEL_MASK;

    // Calculate the top left corner of the best matching block in the
    // reference frame.
    x0 += scaled_mv.col >> SUBPEL_BITS;
    y0 += scaled_mv.row >> SUBPEL_BITS;
    x0_16 += scaled_mv.col;
    y0_16 += scaled_mv.row;

    // Get reference block pointer.
    buf_ptr = ref_frame + y0 * pre_buf->stride + x0;
    buf_stride = pre_buf->stride;

    // Do border extension if there is motion or the
    // width/height is not a multiple of 8 pixels.
    if (is_scaled || scaled_mv.col || scaled_mv.row ||
         (frame_width & 0x7) || (frame_height & 0x7)) {
      // Get reference block bottom right coordinate.
      int x1 = ((x0_16 + (w - 1) * xs) >> SUBPEL_BITS) + 1;
      int y1 = ((y0_16 + (h - 1) * ys) >> SUBPEL_BITS) + 1;
      int x_pad = 0, y_pad = 0;

      if (subpel_x || (sf && sf->x_step_q4 != SUBPEL_SHIFTS)) {
        x0 -= VP9_INTERP_EXTEND - 1;
        x1 += VP9_INTERP_EXTEND;
        x_pad = 1;
      }

      if (subpel_y || (sf && sf->y_step_q4 != SUBPEL_SHIFTS)) {
        y0 -= VP9_INTERP_EXTEND - 1;
        y1 += VP9_INTERP_EXTEND;
        y_pad = 1;
      }

      // Skip border extension if block is inside the frame.
      if (x0 < 0 || x0 > frame_width - 1 || x1 < 0 || x1 > frame_width - 1 ||
          y0 < 0 || y0 > frame_height - 1 || y1 < 0 || y1 > frame_height - 1) {
        uint8_t *buf_ptr1 = ref_frame + y0 * pre_buf->stride + x0;
        // Extend the border.
#if CONFIG_VP9_HIGHBITDEPTH
        if (xd->cur_buf->flags & YV12_FLAG_HIGHBITDEPTH) {
          high_build_mc_border(buf_ptr1,
                               pre_buf->stride,
                               xd->mc_buf_high,
                               x1 - x0 + 1,
                               x0,
                               y0,
                               x1 - x0 + 1,
                               y1 - y0 + 1,
                               frame_width,
                               frame_height);
          buf_stride = x1 - x0 + 1;
          buf_ptr = CONVERT_TO_BYTEPTR(xd->mc_buf_high) +
              y_pad * 3 * buf_stride + x_pad * 3;
        } else {
          build_mc_border(buf_ptr1,
                          pre_buf->stride,
                          xd->mc_buf,
                          x1 - x0 + 1,
                          x0,
                          y0,
                          x1 - x0 + 1,
                          y1 - y0 + 1,
                          frame_width,
                          frame_height);
          buf_stride = x1 - x0 + 1;
          buf_ptr = xd->mc_buf + y_pad * 3 * buf_stride + x_pad * 3;
        }
#else
        build_mc_border(buf_ptr1,
                        pre_buf->stride,
                        xd->mc_buf,
                        x1 - x0 + 1,
                        x0,
                        y0,
                        x1 - x0 + 1,
                        y1 - y0 + 1,
                        frame_width,
                        frame_height);
        buf_stride = x1 - x0 + 1;
        buf_ptr = xd->mc_buf + y_pad * 3 * buf_stride + x_pad * 3;
#endif  // CONFIG_VP9_HIGHBITDEPTH
      }
    }

#if CONFIG_WEDGE_PARTITION
    if (ref && get_wedge_bits(mi->mbmi.sb_type)
        && mi->mbmi.use_wedge_interinter) {
#if CONFIG_VP9_HIGHBITDEPTH
      uint8_t tmp_dst_[2 * CODING_UNIT_SIZE * CODING_UNIT_SIZE];
      uint8_t *tmp_dst =
          (xd->cur_buf->flags & YV12_FLAG_HIGHBITDEPTH) ?
          CONVERT_TO_BYTEPTR(tmp_dst_) : tmp_dst_;
#else
      uint8_t tmp_dst[CODING_UNIT_SIZE * CODING_UNIT_SIZE];
#endif
#if CONFIG_GLOBAL_MOTION
      if (is_global) {
        vp9_warp_plane(gm[ref], pre_buf->buf0,
                       pre_buf->width, pre_buf->height, pre_buf->stride,
                       tmp_dst, (mi_x >> pd->subsampling_x) + x,
                       (mi_y >> pd->subsampling_y) + y, w, h, CODING_UNIT_SIZE,
                       pd->subsampling_x, pd->subsampling_y, xs, ys);
      } else {
#endif  // CONFIG_GLOBAL_MOTION
#if CONFIG_VP9_HIGHBITDEPTH
        if (xd->cur_buf->flags & YV12_FLAG_HIGHBITDEPTH) {
          highbd_inter_predictor(buf_ptr, buf_stride, tmp_dst, CODING_UNIT_SIZE,
                                 subpel_x, subpel_y, sf, w, h, 0, kernel,
                                 xs, ys, xd->bd);
        } else {
          inter_predictor(buf_ptr, buf_stride, tmp_dst, CODING_UNIT_SIZE,
                          subpel_x, subpel_y, sf, w, h, 0, kernel, xs, ys);
        }
#else
        inter_predictor(buf_ptr, buf_stride, tmp_dst, CODING_UNIT_SIZE,
                        subpel_x, subpel_y, sf, w, h, 0, kernel, xs, ys);
#endif  // CONFIG_VP9_HIGHBITDEPTH
#if CONFIG_GLOBAL_MOTION
      }
#endif  // CONFIG_GLOBAL_MOTION
#if CONFIG_SUPERTX
#if CONFIG_VP9_HIGHBITDEPTH
      if (xd->cur_buf->flags & YV12_FLAG_HIGHBITDEPTH) {
        build_masked_compound_extend_highbd(
            dst, dst_buf->stride, tmp_dst, CODING_UNIT_SIZE, plane,
            mi->mbmi.interinter_wedge_index,
            mi->mbmi.sb_type,
            wedge_offset_y, wedge_offset_x, h, w);
      } else {
        build_masked_compound_extend(dst, dst_buf->stride, tmp_dst,
                                     CODING_UNIT_SIZE, plane,
                                     mi->mbmi.interinter_wedge_index,
                                     mi->mbmi.sb_type,
                                     wedge_offset_y, wedge_offset_x, h, w);
      }
#else
      build_masked_compound_extend(dst, dst_buf->stride, tmp_dst,
                                   CODING_UNIT_SIZE, plane,
                                   mi->mbmi.interinter_wedge_index,
                                   mi->mbmi.sb_type,
                                   wedge_offset_y, wedge_offset_x, h, w);
#endif  // CONFIG_VP9_HIGHBITDEPTH
#else   // CONFIG_SUPERTX
#if CONFIG_VP9_HIGHBITDEPTH
      if (xd->cur_buf->flags & YV12_FLAG_HIGHBITDEPTH) {
        build_masked_compound_highbd(dst, dst_buf->stride, tmp_dst,
                                     CODING_UNIT_SIZE,
                                     mi->mbmi.interinter_wedge_index,
                                     mi->mbmi.sb_type, h, w);
      } else {
        build_masked_compound(dst, dst_buf->stride, tmp_dst, CODING_UNIT_SIZE,
                              mi->mbmi.interinter_wedge_index, mi->mbmi.sb_type,
                              h, w);
      }
#else
      build_masked_compound(dst, dst_buf->stride, tmp_dst, CODING_UNIT_SIZE,
                            mi->mbmi.interinter_wedge_index, mi->mbmi.sb_type,
                            h, w);
#endif  // CONFIG_VP9_HIGHBITDEPTH
#endif  // CONFIG_SUPERTX
    } else {
#if CONFIG_GLOBAL_MOTION
      if (is_global) {
        vp9_warp_plane(gm[ref], pre_buf->buf0,
                       pre_buf->width, pre_buf->height, pre_buf->stride, dst,
                       (mi_x >> pd->subsampling_x) + x,
                       (mi_y >> pd->subsampling_y) + y, w, h, dst_buf->stride,
                       pd->subsampling_x, pd->subsampling_y, xs, ys);
      } else {
#endif  // CONFIG_GLOBAL_MOTION
#if CONFIG_VP9_HIGHBITDEPTH
        if (xd->cur_buf->flags & YV12_FLAG_HIGHBITDEPTH)
          highbd_inter_predictor(buf_ptr, buf_stride, dst, dst_buf->stride,
                                 subpel_x, subpel_y, sf, w, h, ref, kernel,
                                 xs, ys, xd->bd);
        else
#endif  // CONFIG_VP9_HIGHBITDEPTH
          inter_predictor(buf_ptr, buf_stride, dst, dst_buf->stride, subpel_x,
                          subpel_y, sf, w, h, ref, kernel, xs, ys);
#if CONFIG_GLOBAL_MOTION
      }
#endif  // CONFIG_GLOBAL_MOTION
    }

#else  // CONFIG_WEDGE_PARTITION

#if CONFIG_GLOBAL_MOTION
    if (is_global) {
      vp9_warp_plane(gm[ref], pre_buf->buf0,
                     pre_buf->width, pre_buf->height, pre_buf->stride, dst,
                     (mi_x >> pd->subsampling_x) + x,
                     (mi_y >> pd->subsampling_y) + y, w, h, dst_buf->stride,
                     pd->subsampling_x, pd->subsampling_y, xs, ys);
    } else {
#endif  // CONFIG_GLOBAL_MOTION
#if CONFIG_VP9_HIGHBITDEPTH
      if (xd->cur_buf->flags & YV12_FLAG_HIGHBITDEPTH) {
        highbd_inter_predictor(buf_ptr, buf_stride, dst, dst_buf->stride,
                               subpel_x, subpel_y, sf, w, h, ref, kernel,
                               xs, ys, xd->bd);
      } else {
        inter_predictor(buf_ptr, buf_stride, dst, dst_buf->stride, subpel_x,
                        subpel_y, sf, w, h, ref, kernel, xs, ys);
      }
#else
      inter_predictor(buf_ptr, buf_stride, dst, dst_buf->stride, subpel_x,
                      subpel_y, sf, w, h, ref, kernel, xs, ys);
#endif  // CONFIG_VP9_HIGHBITDEPTH
#if CONFIG_GLOBAL_MOTION
    }
#endif  // CONFIG_GLOBAL_MOTION
#endif  // CONFIG_WEDGE_PARTITION
  }
}

void vp9_dec_build_inter_predictors_sb(MACROBLOCKD *xd, int mi_row, int mi_col,
                                       BLOCK_SIZE bsize) {
  int plane;
  const int mi_x = mi_col * MI_SIZE;
  const int mi_y = mi_row * MI_SIZE;
  for (plane = 0; plane < MAX_MB_PLANE; ++plane) {
    const BLOCK_SIZE plane_bsize = get_plane_block_size(bsize,
                                                        &xd->plane[plane]);
    const int num_4x4_w = num_4x4_blocks_wide_lookup[plane_bsize];
    const int num_4x4_h = num_4x4_blocks_high_lookup[plane_bsize];
    const int bw = 4 * num_4x4_w;
    const int bh = 4 * num_4x4_h;

    if (xd->mi[0].src_mi->mbmi.sb_type < BLOCK_8X8) {
      int i = 0, x, y;
      assert(bsize == BLOCK_8X8);
      for (y = 0; y < num_4x4_h; ++y)
        for (x = 0; x < num_4x4_w; ++x)
          dec_build_inter_predictors(xd, plane, i++, bw, bh,
                                     4 * x, 4 * y, 4, 4,
#if CONFIG_SUPERTX && CONFIG_WEDGE_PARTITION
                                     0, 0,
#endif
                                     mi_x, mi_y);
    } else {
      dec_build_inter_predictors(xd, plane, 0, bw, bh,
                                 0, 0, bw, bh,
#if CONFIG_SUPERTX && CONFIG_WEDGE_PARTITION
                                 0, 0,
#endif
                                 mi_x, mi_y);
    }
  }
#if CONFIG_INTERINTRA
  if (xd->mi[0].src_mi->mbmi.ref_frame[1] == INTRA_FRAME &&
      is_interintra_allowed(xd->mi[0].src_mi->mbmi.sb_type))
    vp9_build_interintra_predictors(xd, xd->plane[0].dst.buf,
                                    xd->plane[1].dst.buf, xd->plane[2].dst.buf,
                                    xd->plane[0].dst.stride,
                                    xd->plane[1].dst.stride,
                                    xd->plane[2].dst.stride, bsize);
#endif  // CONFIG_INTERINTRA
}

#if CONFIG_SUPERTX
void vp9_dec_build_inter_predictors_sb_sub8x8(MACROBLOCKD *xd,
                                       int mi_row, int mi_col,
                                       BLOCK_SIZE bsize, int block) {
  // Prediction function used in supertx:
  // Use the mv at current block (which is less than 8x8)
  // to get prediction of a block located at (mi_row, mi_col) at size of bsize
  // bsize can be larger than 8x8.
  // block (0-3): the sub8x8 location of current block
  int plane;
  const int mi_x = mi_col * MI_SIZE;
  const int mi_y = mi_row * MI_SIZE;

  // For sub8x8 uv:
  // Skip uv prediction in supertx except the first block (block = 0)
  int max_plane = (block) ? 1 : MAX_MB_PLANE;

  for (plane = 0; plane < max_plane; plane++) {
    const BLOCK_SIZE plane_bsize = get_plane_block_size(bsize,
                                                        &xd->plane[plane]);
    const int num_4x4_w = num_4x4_blocks_wide_lookup[plane_bsize];
    const int num_4x4_h = num_4x4_blocks_high_lookup[plane_bsize];
    const int bw = 4 * num_4x4_w;
    const int bh = 4 * num_4x4_h;

    dec_build_inter_predictors(xd, plane, block, bw, bh,
                                     0, 0, bw, bh,
#if CONFIG_SUPERTX && CONFIG_WEDGE_PARTITION
                                     0, 0,
#endif
                                     mi_x, mi_y);
  }
#if CONFIG_INTERINTRA  // not sure
  if (xd->mi[0].src_mi->mbmi.ref_frame[1] == INTRA_FRAME &&
      is_interintra_allowed(xd->mi[0].src_mi->mbmi.sb_type))
    vp9_build_interintra_predictors(xd, xd->plane[0].dst.buf,
                                    xd->plane[1].dst.buf, xd->plane[2].dst.buf,
                                    xd->plane[0].dst.stride,
                                    xd->plane[1].dst.stride,
                                    xd->plane[2].dst.stride, bsize);
#endif  // CONFIG_INTERINTRA
}

#if CONFIG_WEDGE_PARTITION
void vp9_dec_build_inter_predictors_sb_extend(MACROBLOCKD *xd,
                                              int mi_row, int mi_col,
                                              int mi_row_ori, int mi_col_ori,
                                              BLOCK_SIZE bsize) {
  int plane;
  const int mi_x = mi_col_ori * MI_SIZE;
  const int mi_y = mi_row_ori * MI_SIZE;
  const int wedge_offset_x = (mi_col - mi_col_ori) * MI_SIZE;
  const int wedge_offset_y = (mi_row - mi_row_ori) * MI_SIZE;
  for (plane = 0; plane < MAX_MB_PLANE; ++plane) {
    const BLOCK_SIZE plane_bsize = get_plane_block_size(bsize,
                                                        &xd->plane[plane]);
    const int num_4x4_w = num_4x4_blocks_wide_lookup[plane_bsize];
    const int num_4x4_h = num_4x4_blocks_high_lookup[plane_bsize];
    const int bw = 4 * num_4x4_w;
    const int bh = 4 * num_4x4_h;

    if (xd->mi[0].src_mi->mbmi.sb_type < BLOCK_8X8) {
      int i = 0, x, y;
      assert(bsize == BLOCK_8X8);
      for (y = 0; y < num_4x4_h; ++y)
        for (x = 0; x < num_4x4_w; ++x)
          dec_build_inter_predictors(
              xd, plane, i++, bw, bh, 4 * x, 4 * y, 4, 4,
              wedge_offset_x >> (xd->plane[plane].subsampling_x),
              wedge_offset_y >> (xd->plane[plane].subsampling_y),
              mi_x, mi_y);
    } else {
      dec_build_inter_predictors(
          xd, plane, 0, bw, bh, 0, 0, bw, bh,
          wedge_offset_x >> (xd->plane[plane].subsampling_x),
          wedge_offset_y >> (xd->plane[plane].subsampling_y),
          mi_x, mi_y);
    }
  }
}

void vp9_dec_build_inter_predictors_sb_sub8x8_extend(
    MACROBLOCKD *xd,
    int mi_row, int mi_col,
    int mi_row_ori, int mi_col_ori,
    BLOCK_SIZE bsize, int block) {
  // Sub8x8 prediction for wedge partition in supertx
  int plane;
  const int mi_x = mi_col_ori * MI_SIZE;
  const int mi_y = mi_row_ori * MI_SIZE;
  const int wedge_offset_x = (mi_col - mi_col_ori) * MI_SIZE;
  const int wedge_offset_y = (mi_row - mi_row_ori) * MI_SIZE;

  // For sub8x8 uv:
  // Skip uv prediction in supertx except the first block (block = 0)
  int max_plane = block ? 1 : MAX_MB_PLANE;

  for (plane = 0; plane < max_plane; ++plane) {
    const BLOCK_SIZE plane_bsize = get_plane_block_size(bsize,
                                                        &xd->plane[plane]);
    const int num_4x4_w = num_4x4_blocks_wide_lookup[plane_bsize];
    const int num_4x4_h = num_4x4_blocks_high_lookup[plane_bsize];
    const int bw = 4 * num_4x4_w;
    const int bh = 4 * num_4x4_h;

    dec_build_inter_predictors(
        xd, plane, block, bw, bh, 0, 0, bw, bh,
        wedge_offset_x >> (xd->plane[plane].subsampling_x),
        wedge_offset_y >> (xd->plane[plane].subsampling_y),
        mi_x, mi_y);
  }
}

#endif  // CONFIG_WEDGE_PARTITION

#endif  // CONFIG_SUPERTX

void vp9_setup_dst_planes(struct macroblockd_plane planes[MAX_MB_PLANE],
                          const YV12_BUFFER_CONFIG *src,
                          int mi_row, int mi_col) {
  uint8_t *const buffers[4] = {src->y_buffer, src->u_buffer, src->v_buffer,
                               src->alpha_buffer};
  const int strides[4] = {src->y_stride, src->uv_stride, src->uv_stride,
                          src->alpha_stride};
  const int widths[4] = {src->y_crop_width, src->uv_crop_width,
                         src->uv_crop_width, src->alpha_width};
  const int heights[4] = {src->y_crop_height, src->uv_crop_height,
                          src->uv_crop_height, src->alpha_height};
  int i;

  for (i = 0; i < MAX_MB_PLANE; ++i) {
    struct macroblockd_plane *const pd = &planes[i];
    setup_pred_plane(&pd->dst, widths[i], heights[i],
                     buffers[i], strides[i], mi_row, mi_col, NULL,
                     pd->subsampling_x, pd->subsampling_y);
  }
}

void vp9_setup_pre_planes(MACROBLOCKD *xd, int idx,
                          const YV12_BUFFER_CONFIG *src,
                          int mi_row, int mi_col,
                          const struct scale_factors *sf) {
  if (src != NULL) {
    int i;
    uint8_t *const buffers[4] = {src->y_buffer, src->u_buffer, src->v_buffer,
                                 src->alpha_buffer};
    const int strides[4] = {src->y_stride, src->uv_stride, src->uv_stride,
                            src->alpha_stride};
    const int widths[4] = {src->y_crop_width, src->uv_crop_width,
                           src->uv_crop_width, src->alpha_width};
    const int heights[4] = {src->y_crop_height, src->uv_crop_height,
                            src->uv_crop_height, src->alpha_height};

    for (i = 0; i < MAX_MB_PLANE; ++i) {
      struct macroblockd_plane *const pd = &xd->plane[i];
      setup_pred_plane(&pd->pre[idx], widths[i], heights[i],
                       buffers[i], strides[i], mi_row, mi_col,
                       sf, pd->subsampling_x, pd->subsampling_y);
    }
  }
}

#if CONFIG_WEDGE_PARTITION
// Builds the inter-predictor for the single ref case
// for use in the encoder to search the wedges efficiently.
static void build_inter_predictors_single_buf(MACROBLOCKD *xd,
                                              int plane, int block,
                                              int bw, int bh,
                                              int x, int y, int w, int h,
                                              int mi_x, int mi_y, int ref,
                                              uint8_t *const ext_dst,
                                              int ext_dst_stride) {
  struct macroblockd_plane *const pd = &xd->plane[plane];
  const MODE_INFO *mi = xd->mi[0].src_mi;
#if CONFIG_INTRABC
  const int is_intrabc = is_intrabc_mode(mi->mbmi.mode);
#endif  // CONFIG_INTRABC
  const InterpKernel *kernel = vp9_get_interp_kernel(mi->mbmi.interp_filter);
#if CONFIG_GLOBAL_MOTION
  Global_Motion_Params *gm = &xd->global_motion[mi->mbmi.ref_frame[ref]][0];
  int is_global;
#endif  // CONFIG_GLOBAL_MOTION

  const struct scale_factors *const sf = &xd->block_refs[ref]->sf;
  struct buf_2d *const dst_buf = &pd->dst;
  struct buf_2d *const pre_buf =
#if CONFIG_INTRABC
      is_intrabc ? dst_buf :
#endif  // CONFIG_INTRABC
      &pd->pre[ref];
#if CONFIG_VP9_HIGHBITDEPTH
  uint8_t *const dst =
      (xd->cur_buf->flags & YV12_FLAG_HIGHBITDEPTH ?
      CONVERT_TO_BYTEPTR(ext_dst) : ext_dst) + ext_dst_stride * y + x;
#else
  uint8_t *const dst = ext_dst + ext_dst_stride * y + x;
#endif
  const MV mv = mi->mbmi.sb_type < BLOCK_8X8
      ? average_split_mvs(pd, mi, ref, block)
      : mi->mbmi.mv[ref].as_mv;

  // TODO(jkoleszar): This clamping is done in the incorrect place for the
  // scaling case. It needs to be done on the scaled MV, not the pre-scaling
  // MV. Note however that it performs the subsampling aware scaling so
  // that the result is always q4.
  // mv_precision precision is MV_PRECISION_Q4.
  const MV mv_q4 = clamp_mv_to_umv_border_sb(xd, &mv, bw, bh,
                                             pd->subsampling_x,
                                             pd->subsampling_y);

  uint8_t *pre;
  MV32 scaled_mv;
  int xs, ys, subpel_x, subpel_y;
  const int is_scaled = vp9_is_scaled(sf);
  (void) dst_buf;
#if CONFIG_INTRABC
  assert(!is_intrabc || mi->mbmi.interp_filter == BILINEAR);
#endif  // CONFIG_INTRABC

#if CONFIG_GLOBAL_MOTION
  is_global = (get_y_mode(mi, block) == ZEROMV &&
#if CONFIG_INTRABC
               !is_intrabc &&
#endif
               get_gmtype(gm) == GLOBAL_ROTZOOM);
#endif  // CONFIG_GLOBAL_MOTION

  if (is_scaled) {
#if CONFIG_INTRABC
    assert(!is_intrabc);
#endif  // CONFIG_INTRABC
    pre = pre_buf->buf + scaled_buffer_offset(x, y, pre_buf->stride, sf);
    scaled_mv = vp9_scale_mv(&mv_q4, mi_x + x, mi_y + y, sf);
    xs = sf->x_step_q4;
    ys = sf->y_step_q4;
  } else {
    pre = pre_buf->buf + (y * pre_buf->stride + x);
    scaled_mv.row = mv_q4.row;
    scaled_mv.col = mv_q4.col;
    xs = ys = 16;
  }
  subpel_x = scaled_mv.col & SUBPEL_MASK;
  subpel_y = scaled_mv.row & SUBPEL_MASK;
  pre += (scaled_mv.row >> SUBPEL_BITS) * pre_buf->stride
      + (scaled_mv.col >> SUBPEL_BITS);

#if CONFIG_GLOBAL_MOTION
  if (is_global) {
    vp9_warp_plane(gm, pre_buf->buf0,
                   pre_buf->width, pre_buf->height, pre_buf->stride, dst,
                   (mi_x >> pd->subsampling_x) + x,
                   (mi_y >> pd->subsampling_y) + y, w, h, ext_dst_stride,
                   pd->subsampling_x, pd->subsampling_y, xs, ys);
  } else {
#endif  // CONFIG_GLOBAL_MOTION
#if CONFIG_VP9_HIGHBITDEPTH
    if (xd->cur_buf->flags & YV12_FLAG_HIGHBITDEPTH) {
      highbd_inter_predictor(pre, pre_buf->stride, dst, ext_dst_stride,
                             subpel_x, subpel_y, sf, w, h, 0, kernel,
                             xs, ys, xd->bd);
    } else {
      inter_predictor(pre, pre_buf->stride, dst, ext_dst_stride,
                      subpel_x, subpel_y, sf, w, h, 0, kernel, xs, ys);
    }
#else
    inter_predictor(pre, pre_buf->stride, dst, ext_dst_stride,
                    subpel_x, subpel_y, sf, w, h, 0, kernel, xs, ys);
#endif  // CONFIG_VP9_HIGHBITDEPTH
#if CONFIG_GLOBAL_MOTION
  }
#endif  // CONFIG_GLOBAL_MOTION
}

void vp9_build_inter_predictors_for_planes_single_buf(
    MACROBLOCKD *xd, BLOCK_SIZE bsize,
    int mi_row, int mi_col, int ref,
    uint8_t *ext_dst[3], int ext_dst_stride[3]) {
  const int plane_from = 0;
  const int plane_to = 2;
  int plane;
  const int mi_x = mi_col * MI_SIZE;
  const int mi_y = mi_row * MI_SIZE;
  for (plane = plane_from; plane <= plane_to; ++plane) {
    const BLOCK_SIZE plane_bsize = get_plane_block_size(bsize,
                                                        &xd->plane[plane]);
    const int num_4x4_w = num_4x4_blocks_wide_lookup[plane_bsize];
    const int num_4x4_h = num_4x4_blocks_high_lookup[plane_bsize];
    const int bw = 4 * num_4x4_w;
    const int bh = 4 * num_4x4_h;

    if (xd->mi[0].src_mi->mbmi.sb_type < BLOCK_8X8) {
      int i = 0, x, y;
      assert(bsize == BLOCK_8X8);
      for (y = 0; y < num_4x4_h; ++y)
        for (x = 0; x < num_4x4_w; ++x)
          build_inter_predictors_single_buf(xd, plane, i++, bw, bh,
                                            4 * x, 4 * y, 4, 4,
                                            mi_x, mi_y, ref,
                                            ext_dst[plane],
                                            ext_dst_stride[plane]);
    } else {
      build_inter_predictors_single_buf(xd, plane, 0, bw, bh,
                                        0, 0, bw, bh,
                                        mi_x, mi_y, ref,
                                        ext_dst[plane],
                                        ext_dst_stride[plane]);
    }
  }
}

static void build_wedge_inter_predictor_from_buf(MACROBLOCKD *xd, int plane,
                                                 int block, int bw, int bh,
                                                 int x, int y, int w, int h,
#if CONFIG_SUPERTX
                                                 int wedge_offset_x,
                                                 int wedge_offset_y,
#endif  // CONFIG_SUPERTX
                                                 int mi_x, int mi_y,
                                                 uint8_t *ext_dst0,
                                                 int ext_dst_stride0,
                                                 uint8_t *ext_dst1,
                                                 int ext_dst_stride1) {
  struct macroblockd_plane *const pd = &xd->plane[plane];
  const MODE_INFO *mi = xd->mi[0].src_mi;
  const int is_compound = has_second_ref(&mi->mbmi);
  int ref;
#if CONFIG_INTRABC
  const int is_intrabc = is_intrabc_mode(mi->mbmi.mode);
#endif  // CONFIG_INTRABC
#if CONFIG_GLOBAL_MOTION
  Global_Motion_Params *gm[2];
  gm[0] = &xd->global_motion[mi->mbmi.ref_frame[0]][0];
  if (is_compound)
    gm[1] = &xd->global_motion[mi->mbmi.ref_frame[1]][0];
#endif  // CONFIG_GLOBAL_MOTION
#if CONFIG_INTRABC
  (void) is_intrabc;
  assert(!is_intrabc || mi->mbmi.interp_filter == BILINEAR);
#endif  // CONFIG_INTRABC
  (void) block;
  (void) bw;
  (void) bh;
  (void) mi_x;
  (void) mi_y;

  for (ref = 0; ref < 1 + is_compound; ++ref) {
    struct buf_2d *const dst_buf = &pd->dst;
    uint8_t *const dst = dst_buf->buf + dst_buf->stride * y + x;
#if CONFIG_GLOBAL_MOTION
    const struct scale_factors *const sf = &xd->block_refs[ref]->sf;
    const int is_scaled = vp9_is_scaled(sf);
    int xs, ys;
    struct buf_2d *const pre_buf =
#if CONFIG_INTRABC
        is_intrabc ? dst_buf :
#endif  // CONFIG_INTRABC
        &pd->pre[ref];

    int is_global = (get_y_mode(mi, block) == ZEROMV &&
#if CONFIG_INTRABC
                     !is_intrabc &&
#endif
                     get_gmtype(gm[ref]) == GLOBAL_ROTZOOM);

    if (is_scaled) {
#if CONFIG_INTRABC
      assert(!is_intrabc);
#endif  // CONFIG_INTRABC
      xs = sf->x_step_q4;
      ys = sf->y_step_q4;
    } else {
      xs = ys = 16;
    }
#endif  // CONFIG_GLOBAL_MOTION

    if (ref && get_wedge_bits(mi->mbmi.sb_type)
        && mi->mbmi.use_wedge_interinter) {
#if CONFIG_VP9_HIGHBITDEPTH
      uint8_t tmp_dst_[2 * CODING_UNIT_SIZE * CODING_UNIT_SIZE];
      uint8_t *tmp_dst =
          (xd->cur_buf->flags & YV12_FLAG_HIGHBITDEPTH) ?
          CONVERT_TO_BYTEPTR(tmp_dst_) : tmp_dst_;
#else
      uint8_t tmp_dst[CODING_UNIT_SIZE * CODING_UNIT_SIZE];
#endif  // CONFIG_VP9_HIGHBITDEPTH
#if CONFIG_GLOBAL_MOTION
      if (is_global) {
        vp9_warp_plane(gm[ref], pre_buf->buf0,
                       pre_buf->width, pre_buf->height, pre_buf->stride,
                       tmp_dst, (mi_x >> pd->subsampling_x) + x,
                       (mi_y >> pd->subsampling_y) + y, w, h, CODING_UNIT_SIZE,
                       pd->subsampling_x, pd->subsampling_y, xs, ys);
      } else {
#endif  // CONFIG_GLOBAL_MOTION
#if CONFIG_VP9_HIGHBITDEPTH
        if (xd->cur_buf->flags & YV12_FLAG_HIGHBITDEPTH) {
          int k;
          for (k = 0; k < h; ++k)
            vpx_memcpy(tmp_dst_ + 2 * CODING_UNIT_SIZE * k, ext_dst1 +
                       ext_dst_stride1 * 2 * k, w * 2);
        } else {
          int k;
          for (k = 0; k < h; ++k)
            vpx_memcpy(tmp_dst_ + CODING_UNIT_SIZE * k, ext_dst1 +
                       ext_dst_stride1 * k, w);
        }
#else
        {
          int k;
          for (k = 0; k < h; ++k)
            vpx_memcpy(tmp_dst + CODING_UNIT_SIZE * k, ext_dst1 +
                       ext_dst_stride1 * k, w);
        }
#endif  // CONFIG_VP9_HIGHBITDEPTH
#if CONFIG_GLOBAL_MOTION
      }
#endif  // CONFIG_GLOBAL_MOTION
#if CONFIG_SUPERTX
#if CONFIG_VP9_HIGHBITDEPTH
      if (xd->cur_buf->flags & YV12_FLAG_HIGHBITDEPTH) {
        build_masked_compound_extend_highbd(
            dst, dst_buf->stride, tmp_dst, CODING_UNIT_SIZE, plane,
            mi->mbmi.interinter_wedge_index,
            mi->mbmi.sb_type,
            wedge_offset_y, wedge_offset_x, h, w);
      } else {
        build_masked_compound_extend(
            dst, dst_buf->stride, tmp_dst, CODING_UNIT_SIZE, plane,
            mi->mbmi.interinter_wedge_index,
            mi->mbmi.sb_type,
            wedge_offset_y, wedge_offset_x, h, w);
      }
#else
      build_masked_compound_extend(dst, dst_buf->stride, tmp_dst,
                                   CODING_UNIT_SIZE, plane,
                                   mi->mbmi.interinter_wedge_index,
                                   mi->mbmi.sb_type,
                                   wedge_offset_y, wedge_offset_x, h, w);
#endif  // CONFIG_VP9_HIGHBITDEPTH
#else   // CONFIG_SUPERTX
#if CONFIG_VP9_HIGHBITDEPTH
      if (xd->cur_buf->flags & YV12_FLAG_HIGHBITDEPTH)
        build_masked_compound_highbd(dst, dst_buf->stride, tmp_dst,
                                     CODING_UNIT_SIZE,
                                     mi->mbmi.interinter_wedge_index,
                                     mi->mbmi.sb_type, h, w);
      else
#endif  // CONFIG_VP9_HIGHBITDEPTH
        build_masked_compound(dst, dst_buf->stride, tmp_dst, CODING_UNIT_SIZE,
                              mi->mbmi.interinter_wedge_index,
                              mi->mbmi.sb_type, h, w);
#endif  // CONFIG_SUPERTX
    } else {
#if CONFIG_GLOBAL_MOTION
      if (is_global) {
        vp9_warp_plane(gm[ref], pre_buf->buf0,
                       pre_buf->width, pre_buf->height, pre_buf->stride, dst,
                       (mi_x >> pd->subsampling_x) + x,
                       (mi_y >> pd->subsampling_y) + y, w, h, dst_buf->stride,
                       pd->subsampling_x, pd->subsampling_y, xs, ys);
      } else {
#endif  // CONFIG_GLOBAL_MOTION
#if CONFIG_VP9_HIGHBITDEPTH
        if (xd->cur_buf->flags & YV12_FLAG_HIGHBITDEPTH) {
          int k;
          for (k = 0; k < h; ++k)
            vpx_memcpy(CONVERT_TO_SHORTPTR(dst + dst_buf->stride * k),
                       ext_dst0 + ext_dst_stride0 * 2 * k, w * 2);
        } else {
          int k;
          for (k = 0; k < h; ++k)
            vpx_memcpy(dst + dst_buf->stride * k,
                       ext_dst0 + ext_dst_stride0 * k, w);
        }
#else
        {
          int k;
          for (k = 0; k < h; ++k)
            vpx_memcpy(dst + dst_buf->stride * k,
                       ext_dst0 + ext_dst_stride0 * k, w);
        }
#endif  // CONFIG_VP9_HIGHBITDEPTH
#if CONFIG_GLOBAL_MOTION
      }
#endif  // CONFIG_GLOBAL_MOTION
    }
  }
}

void vp9_build_wedge_inter_predictor_from_buf(
    MACROBLOCKD *xd, BLOCK_SIZE bsize,
    int mi_row, int mi_col,
    uint8_t *ext_dst0[3], int ext_dst_stride0[3],
    uint8_t *ext_dst1[3], int ext_dst_stride1[3]) {
  const int plane_from = 0;
  const int plane_to = 2;
  int plane;
  const int mi_x = mi_col * MI_SIZE;
  const int mi_y = mi_row * MI_SIZE;
  for (plane = plane_from; plane <= plane_to; ++plane) {
    const BLOCK_SIZE plane_bsize = get_plane_block_size(bsize,
                                                        &xd->plane[plane]);
    const int num_4x4_w = num_4x4_blocks_wide_lookup[plane_bsize];
    const int num_4x4_h = num_4x4_blocks_high_lookup[plane_bsize];
    const int bw = 4 * num_4x4_w;
    const int bh = 4 * num_4x4_h;

    if (xd->mi[0].src_mi->mbmi.sb_type < BLOCK_8X8) {
      int i = 0, x, y;
      assert(bsize == BLOCK_8X8);
      for (y = 0; y < num_4x4_h; ++y)
        for (x = 0; x < num_4x4_w; ++x)
          build_wedge_inter_predictor_from_buf(xd, plane, i++, bw, bh,
                                               4 * x, 4 * y, 4, 4,
#if CONFIG_SUPERTX
                                               0, 0,
#endif
                                               mi_x, mi_y,
                                               ext_dst0[plane],
                                               ext_dst_stride0[plane],
                                               ext_dst1[plane],
                                               ext_dst_stride1[plane]);
    } else {
      build_wedge_inter_predictor_from_buf(xd, plane, 0, bw, bh,
                                           0, 0, bw, bh,
#if CONFIG_SUPERTX
                                           0, 0,
#endif
                                           mi_x, mi_y,
                                           ext_dst0[plane],
                                           ext_dst_stride0[plane],
                                           ext_dst1[plane],
                                           ext_dst_stride1[plane]);
    }
  }
}
#endif  // CONFIG_WEDGE_PARTITION
