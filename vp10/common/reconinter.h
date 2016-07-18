/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#ifndef VP10_COMMON_RECONINTER_H_
#define VP10_COMMON_RECONINTER_H_

#include "vp10/common/filter.h"
#include "vp10/common/onyxc_int.h"
#include "vp10/common/vp10_convolve.h"
#include "vpx/vpx_integer.h"

#ifdef __cplusplus
extern "C" {
#endif

static INLINE void inter_predictor(const uint8_t *src, int src_stride,
                                   uint8_t *dst, int dst_stride,
                                   const int subpel_x,
                                   const int subpel_y,
                                   const struct scale_factors *sf,
                                   int w, int h, int ref_idx,
#if CONFIG_DUAL_FILTER
                                   const INTERP_FILTER *interp_filter,
#else
                                   const INTERP_FILTER interp_filter,
#endif
                                   int xs, int ys) {
#if CONFIG_DUAL_FILTER
  InterpFilterParams interp_filter_params_x =
      vp10_get_interp_filter_params(interp_filter[1 + 2 * ref_idx]);
  InterpFilterParams interp_filter_params_y =
      vp10_get_interp_filter_params(interp_filter[0 + 2 * ref_idx]);
#else
  InterpFilterParams interp_filter_params =
      vp10_get_interp_filter_params(interp_filter);
#endif

#if CONFIG_DUAL_FILTER
  if (interp_filter_params_x.taps == SUBPEL_TAPS &&
      interp_filter_params_y.taps == SUBPEL_TAPS &&
      w > 2 && h > 2) {
    const int16_t *kernel_x =
        vp10_get_interp_filter_subpel_kernel(interp_filter_params_x, subpel_x);
    const int16_t *kernel_y =
        vp10_get_interp_filter_subpel_kernel(interp_filter_params_y, subpel_y);
#else
  if (interp_filter_params.taps == SUBPEL_TAPS) {
    const int16_t *kernel_x =
        vp10_get_interp_filter_subpel_kernel(interp_filter_params, subpel_x);
    const int16_t *kernel_y =
        vp10_get_interp_filter_subpel_kernel(interp_filter_params, subpel_y);
#endif
#if CONFIG_EXT_INTERP && SUPPORT_NONINTERPOLATING_FILTERS
    if (IsInterpolatingFilter(interp_filter)) {
      // Interpolating filter
      sf->predict[subpel_x != 0][subpel_y != 0][ref](
          src, src_stride, dst, dst_stride,
          kernel_x, xs, kernel_y, ys, w, h);
    } else {
      sf->predict_ni[subpel_x != 0][subpel_y != 0][ref](
          src, src_stride, dst, dst_stride,
          kernel_x, xs, kernel_y, ys, w, h);
    }
#else
    sf->predict[subpel_x != 0][subpel_y != 0][ref_idx](
        src, src_stride, dst, dst_stride,
        kernel_x, xs, kernel_y, ys, w, h);
#endif  // CONFIG_EXT_INTERP && SUPPORT_NONINTERPOLATING_FILTERS
  } else {
    // ref_idx > 0 means this is the second reference frame
    // first reference frame's prediction result is already in dst
    // therefore we need to average the first and second results
    vp10_convolve(src, src_stride, dst, dst_stride, w, h, interp_filter,
                  subpel_x, xs, subpel_y, ys, ref_idx);
  }
}

#if CONFIG_VP9_HIGHBITDEPTH
static INLINE void highbd_inter_predictor(const uint8_t *src, int src_stride,
                                          uint8_t *dst, int dst_stride,
                                          const int subpel_x,
                                          const int subpel_y,
                                          const struct scale_factors *sf,
                                          int w, int h, int ref,
#if CONFIG_DUAL_FILTER
                                          const INTERP_FILTER *interp_filter,
#else
                                          const INTERP_FILTER interp_filter,
#endif
                                          int xs, int ys, int bd) {
#if CONFIG_DUAL_FILTER
  InterpFilterParams interp_filter_params_x =
      vp10_get_interp_filter_params(interp_filter[1 + 2 * ref]);
  InterpFilterParams interp_filter_params_y =
      vp10_get_interp_filter_params(interp_filter[0 + 2 * ref]);
#else
  InterpFilterParams interp_filter_params =
      vp10_get_interp_filter_params(interp_filter);
#endif

#if CONFIG_DUAL_FILTER
  if (interp_filter_params_x.taps == SUBPEL_TAPS &&
      interp_filter_params_y.taps == SUBPEL_TAPS &&
      w > 2 && h > 2) {
    const int16_t *kernel_x =
        vp10_get_interp_filter_subpel_kernel(interp_filter_params_x, subpel_x);
    const int16_t *kernel_y =
        vp10_get_interp_filter_subpel_kernel(interp_filter_params_y, subpel_y);
#else
  if (interp_filter_params.taps == SUBPEL_TAPS) {
    const int16_t *kernel_x =
        vp10_get_interp_filter_subpel_kernel(interp_filter_params, subpel_x);
    const int16_t *kernel_y =
        vp10_get_interp_filter_subpel_kernel(interp_filter_params, subpel_y);
#endif  // CONFIG_DUAL_FILTER
#if CONFIG_EXT_INTERP && SUPPORT_NONINTERPOLATING_FILTERS
    if (IsInterpolatingFilter(interp_filter)) {
      // Interpolating filter
      sf->highbd_predict[subpel_x != 0][subpel_y != 0][ref](
          src, src_stride, dst, dst_stride,
          kernel_x, xs, kernel_y, ys, w, h, bd);
    } else {
      sf->highbd_predict_ni[subpel_x != 0][subpel_y != 0][ref](
          src, src_stride, dst, dst_stride,
          kernel_x, xs, kernel_y, ys, w, h, bd);
    }
#else
    sf->highbd_predict[subpel_x != 0][subpel_y != 0][ref](
        src, src_stride, dst, dst_stride,
        kernel_x, xs, kernel_y, ys, w, h, bd);
#endif  // CONFIG_EXT_INTERP && SUPPORT_NONINTERPOLATING_FILTERS
  } else {
    // ref > 0 means this is the second reference frame
    // first reference frame's prediction result is already in dst
    // therefore we need to average the first and second results
    int avg = ref > 0;
    vp10_highbd_convolve(src, src_stride, dst, dst_stride, w, h,
                         interp_filter, subpel_x, xs, subpel_y, ys, avg,
                         bd);
  }
}
#endif  // CONFIG_VP9_HIGHBITDEPTH

#if CONFIG_EXT_INTER
// Set to one to use larger codebooks
#define USE_LARGE_WEDGE_CODEBOOK  0

#if USE_LARGE_WEDGE_CODEBOOK
#define MAX_WEDGE_TYPES   (1 << 5)
#else
#define MAX_WEDGE_TYPES   (1 << 4)
#endif

#define MAX_WEDGE_SIZE_LOG2   5   // 32x32
#define MAX_WEDGE_SIZE        (1 << MAX_WEDGE_SIZE_LOG2)
#define MAX_WEDGE_SQUARE      (MAX_WEDGE_SIZE * MAX_WEDGE_SIZE)

#define WEDGE_WEIGHT_BITS 6

#define WEDGE_NONE       -1

// Angles are with respect to horizontal anti-clockwise
typedef enum {
  WEDGE_HORIZONTAL = 0,
  WEDGE_VERTICAL = 1,
  WEDGE_OBLIQUE27 = 2,
  WEDGE_OBLIQUE63 = 3,
  WEDGE_OBLIQUE117 = 4,
  WEDGE_OBLIQUE153 = 5,
  WEDGE_DIRECTIONS
} WedgeDirectionType;

// 3-tuple: {direction, x_offset, y_offset}
typedef struct {
  WedgeDirectionType direction;
  int x_offset;
  int y_offset;
} wedge_code_type;

typedef uint8_t *wedge_masks_type[MAX_WEDGE_TYPES];

typedef struct {
  int bits;
  const wedge_code_type *codebook;
  uint8_t *signflip;
  int smoother;
  wedge_masks_type *masks;
} wedge_params_type;

extern const wedge_params_type wedge_params_lookup[BLOCK_SIZES];

static INLINE int get_wedge_bits_lookup(BLOCK_SIZE sb_type) {
  return wedge_params_lookup[sb_type].bits;
}

static INLINE int is_interinter_wedge_used(BLOCK_SIZE sb_type) {
  (void) sb_type;
  return wedge_params_lookup[sb_type].bits > 0;
}

static INLINE int get_interinter_wedge_bits(BLOCK_SIZE sb_type) {
  const int wbits = wedge_params_lookup[sb_type].bits;
  return (wbits > 0) ? wbits + 1 : 0;
}

static INLINE int is_interintra_wedge_used(BLOCK_SIZE sb_type) {
  (void) sb_type;
  return wedge_params_lookup[sb_type].bits > 0;
}

static INLINE int get_interintra_wedge_bits(BLOCK_SIZE sb_type) {
  return wedge_params_lookup[sb_type].bits;
}
#endif  // CONFIG_EXT_INTER

void build_inter_predictors(MACROBLOCKD *xd, int plane,
#if CONFIG_OBMC
                            int mi_col_offset, int mi_row_offset,
#endif  // CONFIG_OBMC
                            int block,
                            int bw, int bh,
                            int x, int y, int w, int h,
#if CONFIG_SUPERTX && CONFIG_EXT_INTER
                            int wedge_offset_x, int wedge_offset_y,
#endif  // CONFIG_SUPERTX && CONFIG_EXT_INTER
                            int mi_x, int mi_y);

static INLINE void vp10_make_inter_predictor(
    const uint8_t *src,
    int src_stride,
    uint8_t *dst,
    int dst_stride,
    const int subpel_x,
    const int subpel_y,
    const struct scale_factors *sf,
    int w, int h, int ref,
#if CONFIG_DUAL_FILTER
    const INTERP_FILTER *interp_filter,
#else
    const INTERP_FILTER interp_filter,
#endif
    int xs, int ys,
    const MACROBLOCKD *xd) {
  (void) xd;
#if CONFIG_VP9_HIGHBITDEPTH
  if (xd->cur_buf->flags & YV12_FLAG_HIGHBITDEPTH)
    highbd_inter_predictor(src, src_stride, dst, dst_stride,
                           subpel_x, subpel_y, sf, w, h, ref,
                           interp_filter, xs, ys, xd->bd);
  else
#endif  // CONFIG_VP9_HIGHBITDEPTH
    inter_predictor(src, src_stride, dst, dst_stride,
                    subpel_x, subpel_y, sf, w, h, ref,
                    interp_filter, xs, ys);
}

#if CONFIG_EXT_INTER
void vp10_make_masked_inter_predictor(
    const uint8_t *pre,
    int pre_stride,
    uint8_t *dst,
    int dst_stride,
    const int subpel_x,
    const int subpel_y,
    const struct scale_factors *sf,
    int w, int h,
#if CONFIG_DUAL_FILTER
    const INTERP_FILTER *interp_filter,
#else
    const INTERP_FILTER interp_filter,
#endif
    int xs, int ys,
#if CONFIG_SUPERTX
    int wedge_offset_x, int wedge_offset_y,
#endif  // CONFIG_SUPERTX
    const MACROBLOCKD *xd);
#endif  // CONFIG_EXT_INTER

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
static INLINE MV clamp_mv_to_umv_border_sb(const MACROBLOCKD *xd,
                                           const MV *src_mv,
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

static INLINE MV average_split_mvs(const struct macroblockd_plane *pd,
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
      assert(ss_idx <= 3 && ss_idx >= 0);
  }
  return res;
}

void vp10_build_inter_predictor_sub8x8(MACROBLOCKD *xd, int plane,
                                       int i, int ir, int ic,
                                       int mi_row, int mi_col);

void vp10_build_inter_predictors_sby(MACROBLOCKD *xd, int mi_row, int mi_col,
                                     BLOCK_SIZE bsize);

void vp10_build_inter_predictors_sbp(MACROBLOCKD *xd, int mi_row, int mi_col,
                                     BLOCK_SIZE bsize, int plane);

void vp10_build_inter_predictors_sbuv(MACROBLOCKD *xd, int mi_row, int mi_col,
                                      BLOCK_SIZE bsize);

void vp10_build_inter_predictors_sb(MACROBLOCKD *xd, int mi_row, int mi_col,
                                    BLOCK_SIZE bsize);

#if CONFIG_SUPERTX
void vp10_build_inter_predictors_sb_sub8x8_extend(
    MACROBLOCKD *xd,
#if CONFIG_EXT_INTER
    int mi_row_ori, int mi_col_ori,
#endif  // CONFIG_EXT_INTER
    int mi_row, int mi_col,
    BLOCK_SIZE bsize, int block);

void vp10_build_inter_predictors_sb_extend(
    MACROBLOCKD *xd,
#if CONFIG_EXT_INTER
    int mi_row_ori, int mi_col_ori,
#endif  // CONFIG_EXT_INTER
    int mi_row, int mi_col,
    BLOCK_SIZE bsize);
struct macroblockd_plane;
void vp10_build_masked_inter_predictor_complex(
    MACROBLOCKD *xd,
    uint8_t *dst, int dst_stride,
    const uint8_t *pre, int pre_stride,
    int mi_row, int mi_col, int mi_row_ori, int mi_col_ori,
    BLOCK_SIZE bsize, BLOCK_SIZE top_bsize,
    PARTITION_TYPE partition, int plane);
#endif  // CONFIG_SUPERTX

void vp10_build_inter_predictor(const uint8_t *src, int src_stride,
                               uint8_t *dst, int dst_stride,
                               const MV *mv_q3,
                               const struct scale_factors *sf,
                               int w, int h, int do_avg,
#if CONFIG_DUAL_FILTER
                               const INTERP_FILTER *interp_filter,
#else
                               const INTERP_FILTER interp_filter,
#endif
                               enum mv_precision precision,
                               int x, int y);

#if CONFIG_VP9_HIGHBITDEPTH
void vp10_highbd_build_inter_predictor(const uint8_t *src, int src_stride,
                                      uint8_t *dst, int dst_stride,
                                      const MV *mv_q3,
                                      const struct scale_factors *sf,
                                      int w, int h, int do_avg,
#if CONFIG_DUAL_FILTER
                                      const INTERP_FILTER *interp_filter,
#else
                                      const INTERP_FILTER interp_filter,
#endif
                                      enum mv_precision precision,
                                      int x, int y, int bd);
#endif

static INLINE int scaled_buffer_offset(int x_offset, int y_offset, int stride,
                                       const struct scale_factors *sf) {
  const int x = sf ? sf->scale_value_x(x_offset, sf) : x_offset;
  const int y = sf ? sf->scale_value_y(y_offset, sf) : y_offset;
  return y * stride + x;
}

static INLINE void setup_pred_plane(struct buf_2d *dst,
                                    uint8_t *src, int width,
                                    int height, int stride,
                                    int mi_row, int mi_col,
                                    const struct scale_factors *scale,
                                    int subsampling_x, int subsampling_y) {
  const int x = (MI_SIZE * mi_col) >> subsampling_x;
  const int y = (MI_SIZE * mi_row) >> subsampling_y;
  dst->buf = src + scaled_buffer_offset(x, y, stride, scale);
  dst->buf0 = src;
  dst->width = width;
  dst->height = height;
  dst->stride = stride;
}

void vp10_setup_dst_planes(struct macroblockd_plane planes[MAX_MB_PLANE],
                          const YV12_BUFFER_CONFIG *src,
                          int mi_row, int mi_col);

void vp10_setup_pre_planes(MACROBLOCKD *xd, int idx,
                          const YV12_BUFFER_CONFIG *src, int mi_row, int mi_col,
                          const struct scale_factors *sf);

#if CONFIG_DUAL_FILTER
// Detect if the block have sub-pixel level motion vectors
// per component.
static INLINE int has_subpel_mv_component(const MODE_INFO *const mi,
                                          const MACROBLOCKD *const xd,
                                          int dir) {
  const MB_MODE_INFO *const mbmi = &mi->mbmi;
  const BLOCK_SIZE bsize = mbmi->sb_type;
  int plane;
  int ref = (dir >> 1);

  if (bsize >= BLOCK_8X8) {
    if (dir & 0x01) {
      if (mbmi->mv[ref].as_mv.col & SUBPEL_MASK)
        return 1;
    } else {
      if (mbmi->mv[ref].as_mv.row & SUBPEL_MASK)
        return 1;
    }
  } else {
    for (plane = 0; plane < MAX_MB_PLANE; ++plane) {
      const PARTITION_TYPE bp = BLOCK_8X8 - bsize;
      const struct macroblockd_plane *const pd = &xd->plane[plane];
      const int have_vsplit = bp != PARTITION_HORZ;
      const int have_hsplit = bp != PARTITION_VERT;
      const int num_4x4_w = 2 >> ((!have_vsplit) | pd->subsampling_x);
      const int num_4x4_h = 2 >> ((!have_hsplit) | pd->subsampling_y);

      int x, y;
      for (y = 0; y < num_4x4_h; ++y) {
        for (x = 0; x < num_4x4_w; ++x) {
          const MV mv = average_split_mvs(pd, mi, ref, y * 2 + x);
          if (dir & 0x01) {
            if (mv.col & SUBPEL_MASK)
              return 1;
          } else {
            if (mv.row & SUBPEL_MASK)
              return 1;
          }
        }
      }
    }
  }

  return 0;
}
#endif

#if CONFIG_EXT_INTERP
static INLINE int vp10_is_interp_needed(const MACROBLOCKD *const xd) {
  MODE_INFO *const mi = xd->mi[0];
  MB_MODE_INFO *const mbmi = &mi->mbmi;
  const BLOCK_SIZE bsize = mbmi->sb_type;
  const int is_compound = has_second_ref(mbmi);
  int intpel_mv = 1;
  int plane;

#if SUPPORT_NONINTERPOLATING_FILTERS
  // TODO(debargha): This is is currently only for experimentation
  // with non-interpolating filters. Remove later.
  // If any of the filters are non-interpolating, then indicate the
  // interpolation filter always.
  int i;
  for (i = 0; i < SWITCHABLE_FILTERS; ++i) {
    if (!IsInterpolatingFilter(i)) return 1;
  }
#endif

  // For scaled references, interpolation filter is indicated all the time.
  if (vp10_is_scaled(&xd->block_refs[0]->sf))
    return 1;
  if (is_compound && vp10_is_scaled(&xd->block_refs[1]->sf))
    return 1;

  if (bsize < BLOCK_8X8) {
    for (plane = 0; plane < MAX_MB_PLANE; ++plane) {
      const PARTITION_TYPE bp = BLOCK_8X8 - bsize;
      const struct macroblockd_plane *const pd = &xd->plane[plane];
      const int have_vsplit = bp != PARTITION_HORZ;
      const int have_hsplit = bp != PARTITION_VERT;
      const int num_4x4_w = 2 >> ((!have_vsplit) | pd->subsampling_x);
      const int num_4x4_h = 2 >> ((!have_hsplit) | pd->subsampling_y);
      int ref;
      for (ref = 0; ref < 1 + is_compound; ++ref) {
        int x, y;
        for (y = 0; y < num_4x4_h; ++y)
          for (x = 0; x < num_4x4_w; ++x) {
            const MV mv = average_split_mvs(pd, mi, ref, y * 2 + x);
            if (mv_has_subpel(&mv))
              return 1;
          }
      }
    }
    return 0;
  } else {
    intpel_mv = !mv_has_subpel(&mbmi->mv[0].as_mv);
    if (is_compound && intpel_mv) {
      intpel_mv &= !mv_has_subpel(&mbmi->mv[1].as_mv);
    }
  }
  return !intpel_mv;
}
#endif  // CONFIG_EXT_INTERP

#if CONFIG_OBMC
const uint8_t* vp10_get_obmc_mask(int length);
void vp10_build_obmc_inter_prediction(VP10_COMMON *cm,
                                      MACROBLOCKD *xd, int mi_row, int mi_col,
                                      uint8_t *above[MAX_MB_PLANE],
                                      int above_stride[MAX_MB_PLANE],
                                      uint8_t *left[MAX_MB_PLANE],
                                      int left_stride[MAX_MB_PLANE]);
void vp10_build_prediction_by_above_preds(VP10_COMMON *cm,
                                          MACROBLOCKD *xd,
                                          int mi_row, int mi_col,
                                          uint8_t *tmp_buf[MAX_MB_PLANE],
                                          int tmp_width[MAX_MB_PLANE],
                                          int tmp_height[MAX_MB_PLANE],
                                          int tmp_stride[MAX_MB_PLANE]);
void vp10_build_prediction_by_left_preds(VP10_COMMON *cm,
                                         MACROBLOCKD *xd,
                                         int mi_row, int mi_col,
                                         uint8_t *tmp_buf[MAX_MB_PLANE],
                                         int tmp_width[MAX_MB_PLANE],
                                         int tmp_height[MAX_MB_PLANE],
                                         int tmp_stride[MAX_MB_PLANE]);
#endif  // CONFIG_OBMC

#if CONFIG_EXT_INTER
#define MASK_MASTER_SIZE   (2 * MAX_SB_SIZE)
#define MASK_MASTER_STRIDE (2 * MAX_SB_SIZE)

void vp10_init_wedge_masks();

static INLINE const uint8_t *vp10_get_contiguous_soft_mask(int wedge_index,
                                                           int wedge_sign,
                                                           BLOCK_SIZE sb_type) {
  return wedge_params_lookup[sb_type].masks[wedge_sign][wedge_index];
}

const uint8_t *vp10_get_soft_mask(int wedge_index,
                                  int wedge_sign,
                                  BLOCK_SIZE sb_type,
                                  int wedge_offset_x,
                                  int wedge_offset_y);

void vp10_build_interintra_predictors(MACROBLOCKD *xd,
                                      uint8_t *ypred,
                                      uint8_t *upred,
                                      uint8_t *vpred,
                                      int ystride,
                                      int ustride,
                                      int vstride,
                                      BLOCK_SIZE bsize);
void vp10_build_interintra_predictors_sby(MACROBLOCKD *xd,
                                          uint8_t *ypred,
                                          int ystride,
                                          BLOCK_SIZE bsize);
void vp10_build_interintra_predictors_sbc(MACROBLOCKD *xd,
                                          uint8_t *upred,
                                          int ustride,
                                          int plane,
                                          BLOCK_SIZE bsize);
void vp10_build_interintra_predictors_sbuv(MACROBLOCKD *xd,
                                           uint8_t *upred,
                                           uint8_t *vpred,
                                           int ustride, int vstride,
                                           BLOCK_SIZE bsize);

void vp10_build_intra_predictors_for_interintra(
    MACROBLOCKD *xd,
    BLOCK_SIZE bsize, int plane,
    uint8_t *intra_pred, int intra_stride);
void vp10_combine_interintra(
    MACROBLOCKD *xd,
    BLOCK_SIZE bsize, int plane,
    const uint8_t *inter_pred, int inter_stride,
    const uint8_t *intra_pred, int intra_stride);
void vp10_build_interintra_predictors_sbuv(MACROBLOCKD *xd,
                                           uint8_t *upred,
                                           uint8_t *vpred,
                                           int ustride, int vstride,
                                           BLOCK_SIZE bsize);
void vp10_build_interintra_predictors_sby(MACROBLOCKD *xd,
                                          uint8_t *ypred,
                                          int ystride,
                                          BLOCK_SIZE bsize);

// Encoder only
void vp10_build_inter_predictors_for_planes_single_buf(
    MACROBLOCKD *xd, BLOCK_SIZE bsize,
    int plane_from, int plane_to,
    int mi_row, int mi_col, int ref,
                uint8_t *ext_dst[3], int ext_dst_stride[3]);
void vp10_build_wedge_inter_predictor_from_buf(
    MACROBLOCKD *xd, BLOCK_SIZE bsize,
    int plane_from, int plane_to,
    uint8_t *ext_dst0[3], int ext_dst_stride0[3],
    uint8_t *ext_dst1[3], int ext_dst_stride1[3]);
#endif  // CONFIG_EXT_INTER

#ifdef __cplusplus
}  // extern "C"
#endif

#endif  // VP10_COMMON_RECONINTER_H_
