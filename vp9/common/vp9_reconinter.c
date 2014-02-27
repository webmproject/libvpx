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

static void build_mc_border(const uint8_t *src, uint8_t *dst, int stride,
                             int x, int y, int b_w, int b_h, int w, int h) {
  // Get a pointer to the start of the real data for this row.
  const uint8_t *ref_row = src - x - y * stride;

  if (y >= h)
    ref_row += (h - 1) * stride;
  else if (y > 0)
    ref_row += y * stride;

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
      memmove(dst + left, ref_row + x + left, copy);

    if (right)
      memset(dst + left + copy, ref_row[w - 1], right);

    dst += stride;
    ++y;

    if (y > 0 && y < h)
      ref_row += stride;
  } while (--b_h);
}

static void inter_predictor(const uint8_t *src, int src_stride,
                            uint8_t *dst, int dst_stride,
                            const int subpel_x,
                            const int subpel_y,
                            const struct scale_factors *scale,
                            int w, int h, int ref,
                            const struct subpix_fn_table *subpix,
                            int xs, int ys) {
  scale->sfc->predict[subpel_x != 0][subpel_y != 0][ref](
      src, src_stride, dst, dst_stride,
      subpix->filter_x[subpel_x], xs,
      subpix->filter_y[subpel_y], ys,
      w, h);
}

void vp9_build_inter_predictor(const uint8_t *src, int src_stride,
                               uint8_t *dst, int dst_stride,
                               const MV *src_mv,
                               const struct scale_factors *scale,
                               int w, int h, int ref,
                               const struct subpix_fn_table *subpix,
                               enum mv_precision precision) {
  const int is_q4 = precision == MV_PRECISION_Q4;
  const MV mv_q4 = { is_q4 ? src_mv->row : src_mv->row * 2,
                     is_q4 ? src_mv->col : src_mv->col * 2 };
  const struct scale_factors_common *sfc = scale->sfc;
  const MV32 mv = sfc->scale_mv(&mv_q4, scale);
  const int subpel_x = mv.col & SUBPEL_MASK;
  const int subpel_y = mv.row & SUBPEL_MASK;
  src += (mv.row >> SUBPEL_BITS) * src_stride + (mv.col >> SUBPEL_BITS);

  inter_predictor(src, src_stride, dst, dst_stride, subpel_x, subpel_y,
                  scale, w, h, ref, subpix, sfc->x_step_q4, sfc->y_step_q4);
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

#if CONFIG_MASKED_INTERINTER
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
                                  BLOCK_SIZE sb_type,
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
                                BLOCK_SIZE sb_type,
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

void vp9_generate_hard_mask(int mask_index, BLOCK_SIZE sb_type,
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
                                  int mask_index, BLOCK_SIZE sb_type,
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

// TODO(jkoleszar): In principle, pred_w, pred_h are unnecessary, as we could
// calculate the subsampled BLOCK_SIZE, but that type isn't defined for
// sizes smaller than 16x16 yet.
static void build_inter_predictors(MACROBLOCKD *xd, int plane, int block,
                                   int bw, int bh,
                                   int x, int y, int w, int h,
                                   int mi_x, int mi_y) {
  struct macroblockd_plane *const pd = &xd->plane[plane];
  const MODE_INFO *mi = xd->mi_8x8[0];
  const int is_compound = has_second_ref(&mi->mbmi);
  int ref;

  for (ref = 0; ref < 1 + is_compound; ++ref) {
    struct scale_factors *const scale = &xd->scale_factor[ref];
    struct buf_2d *const pre_buf = &pd->pre[ref];
    struct buf_2d *const dst_buf = &pd->dst;
    uint8_t *const dst = dst_buf->buf + dst_buf->stride * y + x;

    // TODO(jkoleszar): All chroma MVs in SPLITMV mode are taken as the
    // same MV (the average of the 4 luma MVs) but we could do something
    // smarter for non-4:2:0. Just punt for now, pending the changes to get
    // rid of SPLITMV mode entirely.
    const MV mv = mi->mbmi.sb_type < BLOCK_8X8
               ? (plane == 0 ? mi->bmi[block].as_mv[ref].as_mv
                             : mi_mv_pred_q4(mi, ref))
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

    if (vp9_is_scaled(scale->sfc)) {
      pre = pre_buf->buf + scaled_buffer_offset(x, y, pre_buf->stride, scale);
      scale->sfc->set_scaled_offsets(scale, mi_y + y, mi_x + x);
      scaled_mv = scale->sfc->scale_mv(&mv_q4, scale);
      xs = scale->sfc->x_step_q4;
      ys = scale->sfc->y_step_q4;
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

#if CONFIG_MASKED_INTERINTER
    if (ref && get_mask_bits(mi->mbmi.sb_type)
        && mi->mbmi.use_masked_interinter) {
      uint8_t tmp_dst[4096];
      inter_predictor(pre, pre_buf->stride, tmp_dst, 64,
                     subpel_x, subpel_y, scale, w, h, 0, &xd->subpix, xs, ys);
      build_masked_compound(dst, dst_buf->stride, tmp_dst, 64,
                            mi->mbmi.mask_index, mi->mbmi.sb_type, h, w);
    } else {
#endif
    inter_predictor(pre, pre_buf->stride, dst, dst_buf->stride,
                    subpel_x, subpel_y, scale, w, h, ref, &xd->subpix, xs, ys);
#if CONFIG_MASKED_INTERINTER
    }
#endif
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

    if (xd->mi_8x8[0]->mbmi.sb_type < BLOCK_8X8) {
      int i = 0, x, y;
      assert(bsize == BLOCK_8X8);
      for (y = 0; y < num_4x4_h; ++y)
        for (x = 0; x < num_4x4_w; ++x)
           build_inter_predictors(xd, plane, i++, bw, bh,
                                  4 * x, 4 * y, 4, 4, mi_x, mi_y);
    } else {
      build_inter_predictors(xd, plane, 0, bw, bh,
                             0, 0, bw, bh, mi_x, mi_y);
    }
  }
}

void vp9_build_inter_predictors_sby(MACROBLOCKD *xd, int mi_row, int mi_col,
                                    BLOCK_SIZE bsize) {
  build_inter_predictors_for_planes(xd, bsize, mi_row, mi_col, 0, 0);
}
void vp9_build_inter_predictors_sbuv(MACROBLOCKD *xd, int mi_row, int mi_col,
                                     BLOCK_SIZE bsize) {
  build_inter_predictors_for_planes(xd, bsize, mi_row, mi_col, 1,
                                    MAX_MB_PLANE - 1);
}
void vp9_build_inter_predictors_sb(MACROBLOCKD *xd, int mi_row, int mi_col,
                                   BLOCK_SIZE bsize) {
#if CONFIG_INTERINTRA
  uint8_t *const y = xd->plane[0].dst.buf;
  uint8_t *const u = xd->plane[1].dst.buf;
  uint8_t *const v = xd->plane[2].dst.buf;
  const int y_stride = xd->plane[0].dst.stride;
  const int uv_stride = xd->plane[1].dst.stride;
#endif
  build_inter_predictors_for_planes(xd, bsize, mi_row, mi_col, 0,
                                    MAX_MB_PLANE - 1);
#if CONFIG_INTERINTRA
  if (xd->mi_8x8[0]->mbmi.ref_frame[1] == INTRA_FRAME &&
      is_interintra_allowed(xd->mi_8x8[0]->mbmi.sb_type)) {
    vp9_build_interintra_predictors(xd, y, u, v,
                                    y_stride, uv_stride, bsize);
  }
#endif
}

// TODO(jingning): This function serves as a placeholder for decoder prediction
// using on demand border extension. It should be moved to /decoder/ directory.
static void dec_build_inter_predictors(MACROBLOCKD *xd, int plane, int block,
                                       int bw, int bh,
                                       int x, int y, int w, int h,
                                       int mi_x, int mi_y) {
  struct macroblockd_plane *const pd = &xd->plane[plane];
  const MODE_INFO *mi = xd->mi_8x8[0];
  const int is_compound = has_second_ref(&mi->mbmi);
  int ref;

  for (ref = 0; ref < 1 + is_compound; ++ref) {
    struct scale_factors *const scale = &xd->scale_factor[ref];
    struct buf_2d *const pre_buf = &pd->pre[ref];
    struct buf_2d *const dst_buf = &pd->dst;
    uint8_t *const dst = dst_buf->buf + dst_buf->stride * y + x;

    // TODO(jkoleszar): All chroma MVs in SPLITMV mode are taken as the
    // same MV (the average of the 4 luma MVs) but we could do something
    // smarter for non-4:2:0. Just punt for now, pending the changes to get
    // rid of SPLITMV mode entirely.
    const MV mv = mi->mbmi.sb_type < BLOCK_8X8
               ? (plane == 0 ? mi->bmi[block].as_mv[ref].as_mv
                             : mi_mv_pred_q4(mi, ref))
               : mi->mbmi.mv[ref].as_mv;

    // TODO(jkoleszar): This clamping is done in the incorrect place for the
    // scaling case. It needs to be done on the scaled MV, not the pre-scaling
    // MV. Note however that it performs the subsampling aware scaling so
    // that the result is always q4.
    // mv_precision precision is MV_PRECISION_Q4.
    const MV mv_q4 = clamp_mv_to_umv_border_sb(xd, &mv, bw, bh,
                                               pd->subsampling_x,
                                               pd->subsampling_y);

    MV32 scaled_mv;
    int xs, ys, x0, y0, x0_16, y0_16, x1, y1, frame_width,
        frame_height, subpel_x, subpel_y;
    uint8_t *ref_frame, *buf_ptr;
    const YV12_BUFFER_CONFIG *ref_buf = xd->ref_buf[ref];

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

    // Get block position in current frame.
    x0 = (-xd->mb_to_left_edge >> (3 + pd->subsampling_x)) + x;
    y0 = (-xd->mb_to_top_edge >> (3 + pd->subsampling_y)) + y;

    // Precision of x0_16 and y0_16 is 1/16th pixel.
    x0_16 = x0 << SUBPEL_BITS;
    y0_16 = y0 << SUBPEL_BITS;

    if (vp9_is_scaled(scale->sfc)) {
      scale->sfc->set_scaled_offsets(scale, mi_y + y, mi_x + x);
      scaled_mv = scale->sfc->scale_mv(&mv_q4, scale);
      xs = scale->sfc->x_step_q4;
      ys = scale->sfc->y_step_q4;
      // Get block position in the scaled reference frame.
      x0 = scale->sfc->scale_value_x(x0, scale->sfc);
      y0 = scale->sfc->scale_value_y(y0, scale->sfc);
      x0_16 = scale->sfc->scale_value_x(x0_16, scale->sfc);
      y0_16 = scale->sfc->scale_value_y(y0_16, scale->sfc);
    } else {
      scaled_mv.row = mv_q4.row;
      scaled_mv.col = mv_q4.col;
      xs = ys = 16;
    }
    subpel_x = scaled_mv.col & SUBPEL_MASK;
    subpel_y = scaled_mv.row & SUBPEL_MASK;

    // Get reference block top left coordinate.
    x0 += scaled_mv.col >> SUBPEL_BITS;
    y0 += scaled_mv.row >> SUBPEL_BITS;
    x0_16 += scaled_mv.col;
    y0_16 += scaled_mv.row;

    // Get reference block bottom right coordinate.
    x1 = ((x0_16 + (w - 1) * xs) >> SUBPEL_BITS) + 1;
    y1 = ((y0_16 + (h - 1) * xs) >> SUBPEL_BITS) + 1;

    // Get reference block pointer.
    buf_ptr = ref_frame + y0 * pre_buf->stride + x0;

    // Do border extension if there is motion or
    // width/height is not a multiple of 8 pixels.
    if (scaled_mv.col || scaled_mv.row ||
        (frame_width & 0x7) || (frame_height & 0x7)) {

      if (subpel_x) {
        x0 -= VP9_INTERP_EXTEND - 1;
        x1 += VP9_INTERP_EXTEND;
      }

      if (subpel_y) {
        y0 -= VP9_INTERP_EXTEND - 1;
        y1 += VP9_INTERP_EXTEND;
      }

      // Skip border extension if block is inside the frame.
      if (x0 < 0 || x0 > frame_width - 1 || x1 < 0 || x1 > frame_width ||
          y0 < 0 || y0 > frame_height - 1 || y1 < 0 || y1 > frame_height - 1) {
        uint8_t *buf_ptr1 = ref_frame + y0 * pre_buf->stride + x0;
        // Extend the border.
        build_mc_border(buf_ptr1, buf_ptr1, pre_buf->stride, x0, y0, x1 - x0,
                        y1 - y0, frame_width, frame_height);
      }
    }

#if CONFIG_MASKED_INTERINTER
    if (ref && get_mask_bits(mi->mbmi.sb_type)
        && mi->mbmi.use_masked_interinter) {
      uint8_t tmp_dst[4096];
      inter_predictor(buf_ptr, pre_buf->stride, tmp_dst, 64,
                     subpel_x, subpel_y, scale, w, h, 0, &xd->subpix, xs, ys);
      build_masked_compound(dst, dst_buf->stride, tmp_dst, 64,
                            mi->mbmi.mask_index, mi->mbmi.sb_type, h, w);
    } else {
#endif
    inter_predictor(buf_ptr, pre_buf->stride, dst, dst_buf->stride, subpel_x,
                    subpel_y, scale, w, h, ref, &xd->subpix, xs, ys);
#if CONFIG_MASKED_INTERINTER
    }
#endif
  }
}

void vp9_dec_build_inter_predictors_sb(MACROBLOCKD *xd, int mi_row, int mi_col,
                                       BLOCK_SIZE bsize) {
  int plane;
  const int mi_x = mi_col * MI_SIZE;
  const int mi_y = mi_row * MI_SIZE;
#if CONFIG_INTERINTRA
  uint8_t *const y = xd->plane[0].dst.buf;
  uint8_t *const u = xd->plane[1].dst.buf;
  uint8_t *const v = xd->plane[2].dst.buf;
  const int y_stride = xd->plane[0].dst.stride;
  const int uv_stride = xd->plane[1].dst.stride;
#endif
  for (plane = 0; plane < MAX_MB_PLANE; ++plane) {
    const BLOCK_SIZE plane_bsize = get_plane_block_size(bsize,
                                                        &xd->plane[plane]);
    const int num_4x4_w = num_4x4_blocks_wide_lookup[plane_bsize];
    const int num_4x4_h = num_4x4_blocks_high_lookup[plane_bsize];
    const int bw = 4 * num_4x4_w;
    const int bh = 4 * num_4x4_h;

    if (xd->mi_8x8[0]->mbmi.sb_type < BLOCK_8X8) {
      int i = 0, x, y;
      assert(bsize == BLOCK_8X8);
      for (y = 0; y < num_4x4_h; ++y)
        for (x = 0; x < num_4x4_w; ++x)
          dec_build_inter_predictors(xd, plane, i++, bw, bh,
                                     4 * x, 4 * y, 4, 4, mi_x, mi_y);
    } else {
      dec_build_inter_predictors(xd, plane, 0, bw, bh,
                                 0, 0, bw, bh, mi_x, mi_y);
    }
  }
#if CONFIG_INTERINTRA
  if (xd->mi_8x8[0]->mbmi.ref_frame[1] == INTRA_FRAME &&
      is_interintra_allowed(xd->mi_8x8[0]->mbmi.sb_type)) {
    vp9_build_interintra_predictors(xd, y, u, v,
                                    y_stride, uv_stride, bsize);
  }
#endif
}

// TODO(dkovalev: find better place for this function)
void vp9_setup_scale_factors(VP9_COMMON *cm, int i) {
  const int ref = cm->active_ref_idx[i];
  struct scale_factors *const sf = &cm->active_ref_scale[i];
  struct scale_factors_common *const sfc = &cm->active_ref_scale_comm[i];
  if (ref >= cm->fb_count) {
    vp9_zero(*sf);
    vp9_zero(*sfc);
  } else {
    YV12_BUFFER_CONFIG *const fb = &cm->yv12_fb[ref];
    vp9_setup_scale_factors_for_frame(sf, sfc,
                                      fb->y_crop_width, fb->y_crop_height,
                                      cm->width, cm->height);
  }
}

