/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#ifndef VP9_COMMON_VP9_RECONINTER_H_
#define VP9_COMMON_VP9_RECONINTER_H_

#include "vpx/vpx_integer.h"
#include "vp9/common/vp9_onyxc_int.h"

struct subpix_fn_table;
void vp9_build_inter_predictors_sby(MACROBLOCKD *xd,
                                    int mb_row,
                                    int mb_col,
                                    BLOCK_SIZE_TYPE bsize);

void vp9_build_inter_predictors_sbuv(MACROBLOCKD *xd,
                                     int mb_row,
                                     int mb_col,
                                     BLOCK_SIZE_TYPE bsize);

void vp9_build_inter_predictors_sb(MACROBLOCKD *mb,
                                   int mb_row, int mb_col,
                                   BLOCK_SIZE_TYPE bsize);

void vp9_setup_interp_filters(MACROBLOCKD *xd,
                              INTERPOLATIONFILTERTYPE filter,
                              VP9_COMMON *cm);

void vp9_setup_scale_factors_for_frame(struct scale_factors *scale,
                                       int other_w, int other_h,
                                       int this_w, int this_h);

void vp9_build_inter_predictor(const uint8_t *src, int src_stride,
                               uint8_t *dst, int dst_stride,
                               const int_mv *mv_q3,
                               const struct scale_factors *scale,
                               int w, int h, int do_avg,
                               const struct subpix_fn_table *subpix);

void vp9_build_inter_predictor_q4(const uint8_t *src, int src_stride,
                                  uint8_t *dst, int dst_stride,
                                  const int_mv *mv_q4,
                                  const struct scale_factors *scale,
                                  int w, int h, int do_avg,
                                  const struct subpix_fn_table *subpix);

static int scaled_buffer_offset(int x_offset, int y_offset, int stride,
                                const struct scale_factors *scale) {
  const int x = scale ? scale->scale_value_x(x_offset, scale) : x_offset;
  const int y = scale ? scale->scale_value_y(y_offset, scale) : y_offset;
  return y * stride + x;
}

static void setup_pred_plane(struct buf_2d *dst,
                             uint8_t *src, int stride,
                             int mi_row, int mi_col,
                             const struct scale_factors *scale,
                             int subsampling_x, int subsampling_y) {
  const int x = (MI_SIZE * mi_col) >> subsampling_x;
  const int y = (MI_SIZE * mi_row) >> subsampling_y;
  dst->buf = src + scaled_buffer_offset(x, y, stride, scale);
  dst->stride = stride;
}

// TODO(jkoleszar): audit all uses of this that don't set mb_row, mb_col
static void setup_dst_planes(MACROBLOCKD *xd,
                             const YV12_BUFFER_CONFIG *src,
                             int mi_row, int mi_col) {
  uint8_t *buffers[4] = {src->y_buffer, src->u_buffer, src->v_buffer,
                         src->alpha_buffer};
  int strides[4] = {src->y_stride, src->uv_stride, src->uv_stride,
                    src->alpha_stride};
  int i;

  for (i = 0; i < MAX_MB_PLANE; ++i) {
    struct macroblockd_plane *pd = &xd->plane[i];
    setup_pred_plane(&pd->dst, buffers[i], strides[i], mi_row, mi_col, NULL,
                     pd->subsampling_x, pd->subsampling_y);
  }
}

static void setup_pre_planes(MACROBLOCKD *xd,
                             const YV12_BUFFER_CONFIG *src0,
                             const YV12_BUFFER_CONFIG *src1,
                             int mi_row, int mi_col,
                             const struct scale_factors *scale,
                             const struct scale_factors *scale_uv) {
  const YV12_BUFFER_CONFIG *srcs[2] = {src0, src1};
  int i, j;

  for (i = 0; i < 2; ++i) {
    const YV12_BUFFER_CONFIG *src = srcs[i];
    if (src) {
      uint8_t* buffers[4] = {src->y_buffer, src->u_buffer, src->v_buffer,
                             src->alpha_buffer};
      int strides[4] = {src->y_stride, src->uv_stride, src->uv_stride,
                        src->alpha_stride};

      for (j = 0; j < MAX_MB_PLANE; ++j) {
        struct macroblockd_plane *pd = &xd->plane[j];
        const struct scale_factors *sf = j ? scale_uv : scale;
        setup_pred_plane(&pd->pre[i],
                         buffers[j], strides[j],
                         mi_row, mi_col, sf ? &sf[i] : NULL,
                         pd->subsampling_x, pd->subsampling_y);
      }
    }
  }
}

static void set_scale_factors(MACROBLOCKD *xd,
    int ref0, int ref1,
    struct scale_factors scale_factor[MAX_REF_FRAMES]) {

  xd->scale_factor[0] = scale_factor[ref0 >= 0 ? ref0 : 0];
  xd->scale_factor[1] = scale_factor[ref1 >= 0 ? ref1 : 0];
  xd->scale_factor_uv[0] = xd->scale_factor[0];
  xd->scale_factor_uv[1] = xd->scale_factor[1];
}

void vp9_setup_scale_factors(VP9_COMMON *cm, int i);

#endif  // VP9_COMMON_VP9_RECONINTER_H_
