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

void vp9_setup_interp_filters(MACROBLOCKD *xd,
                              INTERPOLATIONFILTERTYPE mcomp_filter_type,
                              VP9_COMMON *cm) {
  // TODO(agrange): Investigate the best choice of functions to use here
  // for EIGHTTAP_SMOOTH. Since it is not interpolating, need to choose what
  // to do at full-pel offsets. The current selection, where the filter is
  // applied in one direction only, and not at all for 0,0, seems to give the
  // best quality, but it may be worth trying an additional mode that does
  // do the filtering on full-pel.
  xd->subpix.predict[0][0][0] = vp9_convolve_copy;
  xd->subpix.predict[0][0][1] = vp9_convolve_avg;
  xd->subpix.predict[0][1][0] = vp9_convolve8_vert;
  xd->subpix.predict[0][1][1] = vp9_convolve8_avg_vert;
  xd->subpix.predict[1][0][0] = vp9_convolve8_horiz;
  xd->subpix.predict[1][0][1] = vp9_convolve8_avg_horiz;
  xd->subpix.predict[1][1][0] = vp9_convolve8;
  xd->subpix.predict[1][1][1] = vp9_convolve8_avg;

  xd->subpix.x_step_q4 = 16;
  xd->subpix.y_step_q4 = 16;
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
#if CONFIG_ENABLE_6TAP
    case SIXTAP:
      xd->subpix.filter_x = xd->subpix.filter_y = vp9_sub_pel_filters_6;
      break;
#endif
  }
}

void vp9_copy_mem16x16_c(const uint8_t *src,
                         int src_stride,
                         uint8_t *dst,
                         int dst_stride) {
  int r;

  for (r = 0; r < 16; r++) {
#if !(CONFIG_FAST_UNALIGNED)
    dst[0] = src[0];
    dst[1] = src[1];
    dst[2] = src[2];
    dst[3] = src[3];
    dst[4] = src[4];
    dst[5] = src[5];
    dst[6] = src[6];
    dst[7] = src[7];
    dst[8] = src[8];
    dst[9] = src[9];
    dst[10] = src[10];
    dst[11] = src[11];
    dst[12] = src[12];
    dst[13] = src[13];
    dst[14] = src[14];
    dst[15] = src[15];

#else
    ((uint32_t *)dst)[0] = ((const uint32_t *)src)[0];
    ((uint32_t *)dst)[1] = ((const uint32_t *)src)[1];
    ((uint32_t *)dst)[2] = ((const uint32_t *)src)[2];
    ((uint32_t *)dst)[3] = ((const uint32_t *)src)[3];

#endif
    src += src_stride;
    dst += dst_stride;
  }
}

void vp9_copy_mem8x8_c(const uint8_t *src,
                       int src_stride,
                       uint8_t *dst,
                       int dst_stride) {
  int r;

  for (r = 0; r < 8; r++) {
#if !(CONFIG_FAST_UNALIGNED)
    dst[0] = src[0];
    dst[1] = src[1];
    dst[2] = src[2];
    dst[3] = src[3];
    dst[4] = src[4];
    dst[5] = src[5];
    dst[6] = src[6];
    dst[7] = src[7];
#else
    ((uint32_t *)dst)[0] = ((const uint32_t *)src)[0];
    ((uint32_t *)dst)[1] = ((const uint32_t *)src)[1];
#endif
    src += src_stride;
    dst += dst_stride;
  }
}

void vp9_copy_mem8x4_c(const uint8_t *src,
                       int src_stride,
                       uint8_t *dst,
                       int dst_stride) {
  int r;

  for (r = 0; r < 4; r++) {
#if !(CONFIG_FAST_UNALIGNED)
    dst[0] = src[0];
    dst[1] = src[1];
    dst[2] = src[2];
    dst[3] = src[3];
    dst[4] = src[4];
    dst[5] = src[5];
    dst[6] = src[6];
    dst[7] = src[7];
#else
    ((uint32_t *)dst)[0] = ((const uint32_t *)src)[0];
    ((uint32_t *)dst)[1] = ((const uint32_t *)src)[1];
#endif
    src += src_stride;
    dst += dst_stride;
  }
}

void vp9_build_inter_predictors_b(BLOCKD *d, int pitch,
                                  struct subpix_fn_table *subpix) {
  uint8_t *ptr_base;
  uint8_t *ptr;
  uint8_t *pred_ptr = d->predictor;
  int_mv mv;

  ptr_base = *(d->base_pre);
  mv.as_int = d->bmi.as_mv.first.as_int;
  ptr = ptr_base + d->pre + (mv.as_mv.row >> 3) * d->pre_stride +
        (mv.as_mv.col >> 3);

  subpix->predict[!!(mv.as_mv.col & 7)][!!(mv.as_mv.row & 7)][0](
      ptr, d->pre_stride, pred_ptr, pitch,
      subpix->filter_x[(mv.as_mv.col & 7) << 1], subpix->x_step_q4,
      subpix->filter_y[(mv.as_mv.row & 7) << 1], subpix->y_step_q4,
      4, 4);
}

/*
 * Similar to vp9_build_inter_predictors_b(), but instead of storing the
 * results in d->predictor, we average the contents of d->predictor (which
 * come from an earlier call to vp9_build_inter_predictors_b()) with the
 * predictor of the second reference frame / motion vector.
 */
void vp9_build_2nd_inter_predictors_b(BLOCKD *d, int pitch,
                                      struct subpix_fn_table *subpix) {
  uint8_t *ptr_base;
  uint8_t *ptr;
  uint8_t *pred_ptr = d->predictor;
  int_mv mv;

  ptr_base = *(d->base_second_pre);
  mv.as_int = d->bmi.as_mv.second.as_int;
  ptr = ptr_base + d->pre + (mv.as_mv.row >> 3) * d->pre_stride +
        (mv.as_mv.col >> 3);

  subpix->predict[!!(mv.as_mv.col & 7)][!!(mv.as_mv.row & 7)][1](
      ptr, d->pre_stride, pred_ptr, pitch,
      subpix->filter_x[(mv.as_mv.col & 7) << 1], subpix->x_step_q4,
      subpix->filter_y[(mv.as_mv.row & 7) << 1], subpix->y_step_q4,
      4, 4);
}

void vp9_build_inter_predictors4b(MACROBLOCKD *xd, BLOCKD *d, int pitch) {
  uint8_t *ptr_base;
  uint8_t *ptr;
  uint8_t *pred_ptr = d->predictor;
  int_mv mv;

  ptr_base = *(d->base_pre);
  mv.as_int = d->bmi.as_mv.first.as_int;
  ptr = ptr_base + d->pre + (mv.as_mv.row >> 3) * d->pre_stride +
        (mv.as_mv.col >> 3);

  xd->subpix.predict[!!(mv.as_mv.col & 7)][!!(mv.as_mv.row & 7)][0](
      ptr, d->pre_stride, pred_ptr, pitch,
      xd->subpix.filter_x[(mv.as_mv.col & 7) << 1], xd->subpix.x_step_q4,
      xd->subpix.filter_y[(mv.as_mv.row & 7) << 1], xd->subpix.y_step_q4,
      8, 8);
}

/*
 * Similar to build_inter_predictors_4b(), but instead of storing the
 * results in d->predictor, we average the contents of d->predictor (which
 * come from an earlier call to build_inter_predictors_4b()) with the
 * predictor of the second reference frame / motion vector.
 */
void vp9_build_2nd_inter_predictors4b(MACROBLOCKD *xd,
                                      BLOCKD *d, int pitch) {
  uint8_t *ptr_base;
  uint8_t *ptr;
  uint8_t *pred_ptr = d->predictor;
  int_mv mv;

  ptr_base = *(d->base_second_pre);
  mv.as_int = d->bmi.as_mv.second.as_int;
  ptr = ptr_base + d->pre + (mv.as_mv.row >> 3) * d->pre_stride +
        (mv.as_mv.col >> 3);

  xd->subpix.predict[!!(mv.as_mv.col & 7)][!!(mv.as_mv.row & 7)][1](
      ptr, d->pre_stride, pred_ptr, pitch,
      xd->subpix.filter_x[(mv.as_mv.col & 7) << 1], xd->subpix.x_step_q4,
      xd->subpix.filter_y[(mv.as_mv.row & 7) << 1], xd->subpix.y_step_q4,
      8, 8);
}

static void build_inter_predictors2b(MACROBLOCKD *xd, BLOCKD *d, int pitch) {
  uint8_t *ptr_base;
  uint8_t *ptr;
  uint8_t *pred_ptr = d->predictor;
  int_mv mv;

  ptr_base = *(d->base_pre);
  mv.as_int = d->bmi.as_mv.first.as_int;
  ptr = ptr_base + d->pre + (mv.as_mv.row >> 3) * d->pre_stride +
        (mv.as_mv.col >> 3);

  xd->subpix.predict[!!(mv.as_mv.col & 7)][!!(mv.as_mv.row & 7)][0](
      ptr, d->pre_stride, pred_ptr, pitch,
      xd->subpix.filter_x[(mv.as_mv.col & 7) << 1], xd->subpix.x_step_q4,
      xd->subpix.filter_y[(mv.as_mv.row & 7) << 1], xd->subpix.y_step_q4,
      8, 4);
}

/*encoder only*/
void vp9_build_inter4x4_predictors_mbuv(MACROBLOCKD *xd) {
  int i, j;
  BLOCKD *blockd = xd->block;

  /* build uv mvs */
  for (i = 0; i < 2; i++) {
    for (j = 0; j < 2; j++) {
      int yoffset = i * 8 + j * 2;
      int uoffset = 16 + i * 2 + j;
      int voffset = 20 + i * 2 + j;
      int temp;

      temp = blockd[yoffset  ].bmi.as_mv.first.as_mv.row
             + blockd[yoffset + 1].bmi.as_mv.first.as_mv.row
             + blockd[yoffset + 4].bmi.as_mv.first.as_mv.row
             + blockd[yoffset + 5].bmi.as_mv.first.as_mv.row;

      if (temp < 0) temp -= 4;
      else temp += 4;

      xd->block[uoffset].bmi.as_mv.first.as_mv.row = (temp / 8) &
        xd->fullpixel_mask;

      temp = blockd[yoffset  ].bmi.as_mv.first.as_mv.col
             + blockd[yoffset + 1].bmi.as_mv.first.as_mv.col
             + blockd[yoffset + 4].bmi.as_mv.first.as_mv.col
             + blockd[yoffset + 5].bmi.as_mv.first.as_mv.col;

      if (temp < 0) temp -= 4;
      else temp += 4;

      blockd[uoffset].bmi.as_mv.first.as_mv.col = (temp / 8) &
        xd->fullpixel_mask;

      blockd[voffset].bmi.as_mv.first.as_mv.row =
        blockd[uoffset].bmi.as_mv.first.as_mv.row;
      blockd[voffset].bmi.as_mv.first.as_mv.col =
        blockd[uoffset].bmi.as_mv.first.as_mv.col;

      if (xd->mode_info_context->mbmi.second_ref_frame > 0) {
        temp = blockd[yoffset  ].bmi.as_mv.second.as_mv.row
               + blockd[yoffset + 1].bmi.as_mv.second.as_mv.row
               + blockd[yoffset + 4].bmi.as_mv.second.as_mv.row
               + blockd[yoffset + 5].bmi.as_mv.second.as_mv.row;

        if (temp < 0) {
          temp -= 4;
        } else {
          temp += 4;
        }

        blockd[uoffset].bmi.as_mv.second.as_mv.row = (temp / 8) &
          xd->fullpixel_mask;

        temp = blockd[yoffset  ].bmi.as_mv.second.as_mv.col
               + blockd[yoffset + 1].bmi.as_mv.second.as_mv.col
               + blockd[yoffset + 4].bmi.as_mv.second.as_mv.col
               + blockd[yoffset + 5].bmi.as_mv.second.as_mv.col;

        if (temp < 0) {
          temp -= 4;
        } else {
          temp += 4;
        }

        blockd[uoffset].bmi.as_mv.second.as_mv.col = (temp / 8) &
          xd->fullpixel_mask;

        blockd[voffset].bmi.as_mv.second.as_mv.row =
          blockd[uoffset].bmi.as_mv.second.as_mv.row;
        blockd[voffset].bmi.as_mv.second.as_mv.col =
          blockd[uoffset].bmi.as_mv.second.as_mv.col;
      }
    }
  }

  for (i = 16; i < 24; i += 2) {
    BLOCKD *d0 = &blockd[i];
    BLOCKD *d1 = &blockd[i + 1];

    if (d0->bmi.as_mv.first.as_int == d1->bmi.as_mv.first.as_int)
      build_inter_predictors2b(xd, d0, 8);
    else {
      vp9_build_inter_predictors_b(d0, 8, &xd->subpix);
      vp9_build_inter_predictors_b(d1, 8, &xd->subpix);
    }

    if (xd->mode_info_context->mbmi.second_ref_frame > 0) {
      vp9_build_2nd_inter_predictors_b(d0, 8, &xd->subpix);
      vp9_build_2nd_inter_predictors_b(d1, 8, &xd->subpix);
    }
  }
}

static void clamp_mv_to_umv_border(MV *mv, const MACROBLOCKD *xd) {
  /* If the MV points so far into the UMV border that no visible pixels
   * are used for reconstruction, the subpel part of the MV can be
   * discarded and the MV limited to 16 pixels with equivalent results.
   *
   * This limit kicks in at 19 pixels for the top and left edges, for
   * the 16 pixels plus 3 taps right of the central pixel when subpel
   * filtering. The bottom and right edges use 16 pixels plus 2 pixels
   * left of the central pixel when filtering.
   */
  if (mv->col < (xd->mb_to_left_edge - ((16 + VP9_INTERP_EXTEND) << 3)))
    mv->col = xd->mb_to_left_edge - (16 << 3);
  else if (mv->col > xd->mb_to_right_edge + ((15 + VP9_INTERP_EXTEND) << 3))
    mv->col = xd->mb_to_right_edge + (16 << 3);

  if (mv->row < (xd->mb_to_top_edge - ((16 + VP9_INTERP_EXTEND) << 3)))
    mv->row = xd->mb_to_top_edge - (16 << 3);
  else if (mv->row > xd->mb_to_bottom_edge + ((15 + VP9_INTERP_EXTEND) << 3))
    mv->row = xd->mb_to_bottom_edge + (16 << 3);
}

/* A version of the above function for chroma block MVs.*/
static void clamp_uvmv_to_umv_border(MV *mv, const MACROBLOCKD *xd) {
  const int extend = VP9_INTERP_EXTEND;

  mv->col = (2 * mv->col < (xd->mb_to_left_edge - ((16 + extend) << 3))) ?
            (xd->mb_to_left_edge - (16 << 3)) >> 1 : mv->col;
  mv->col = (2 * mv->col > xd->mb_to_right_edge + ((15 + extend) << 3)) ?
            (xd->mb_to_right_edge + (16 << 3)) >> 1 : mv->col;

  mv->row = (2 * mv->row < (xd->mb_to_top_edge - ((16 + extend) << 3))) ?
            (xd->mb_to_top_edge - (16 << 3)) >> 1 : mv->row;
  mv->row = (2 * mv->row > xd->mb_to_bottom_edge + ((15 + extend) << 3)) ?
            (xd->mb_to_bottom_edge + (16 << 3)) >> 1 : mv->row;
}

/*encoder only*/
void vp9_build_1st_inter16x16_predictors_mby(MACROBLOCKD *xd,
                                             uint8_t *dst_y,
                                             int dst_ystride,
                                             int clamp_mvs) {
  uint8_t *ptr_base = xd->pre.y_buffer;
  uint8_t *ptr;
  int pre_stride = xd->block[0].pre_stride;
  int_mv ymv;

  ymv.as_int = xd->mode_info_context->mbmi.mv[0].as_int;

  if (clamp_mvs)
    clamp_mv_to_umv_border(&ymv.as_mv, xd);

  ptr = ptr_base + (ymv.as_mv.row >> 3) * pre_stride + (ymv.as_mv.col >> 3);

  xd->subpix.predict[!!(ymv.as_mv.col & 7)][!!(ymv.as_mv.row & 7)][0](
      ptr, pre_stride, dst_y, dst_ystride,
      xd->subpix.filter_x[(ymv.as_mv.col & 7) << 1], xd->subpix.x_step_q4,
      xd->subpix.filter_y[(ymv.as_mv.row & 7) << 1], xd->subpix.y_step_q4,
      16, 16);
}

void vp9_build_1st_inter16x16_predictors_mbuv(MACROBLOCKD *xd,
                                              uint8_t *dst_u,
                                              uint8_t *dst_v,
                                              int dst_uvstride) {
  int offset;
  uint8_t *uptr, *vptr;
  int pre_stride = xd->block[0].pre_stride;
  int_mv _o16x16mv;
  int_mv _16x16mv;

  _16x16mv.as_int = xd->mode_info_context->mbmi.mv[0].as_int;

  if (xd->mode_info_context->mbmi.need_to_clamp_mvs)
    clamp_mv_to_umv_border(&_16x16mv.as_mv, xd);

  _o16x16mv = _16x16mv;
  /* calc uv motion vectors */
  if (_16x16mv.as_mv.row < 0)
    _16x16mv.as_mv.row -= 1;
  else
    _16x16mv.as_mv.row += 1;

  if (_16x16mv.as_mv.col < 0)
    _16x16mv.as_mv.col -= 1;
  else
    _16x16mv.as_mv.col += 1;

  _16x16mv.as_mv.row /= 2;
  _16x16mv.as_mv.col /= 2;

  _16x16mv.as_mv.row &= xd->fullpixel_mask;
  _16x16mv.as_mv.col &= xd->fullpixel_mask;

  pre_stride >>= 1;
  offset = (_16x16mv.as_mv.row >> 3) * pre_stride + (_16x16mv.as_mv.col >> 3);
  uptr = xd->pre.u_buffer + offset;
  vptr = xd->pre.v_buffer + offset;

  xd->subpix.predict[!!(_o16x16mv.as_mv.col & 15)]
                    [!!(_o16x16mv.as_mv.row & 15)][0](
      uptr, pre_stride, dst_u, dst_uvstride,
      xd->subpix.filter_x[_o16x16mv.as_mv.col & 15], xd->subpix.x_step_q4,
      xd->subpix.filter_y[_o16x16mv.as_mv.row & 15], xd->subpix.y_step_q4,
      8, 8);

  xd->subpix.predict[!!(_o16x16mv.as_mv.col & 15)]
                    [!!(_o16x16mv.as_mv.row & 15)][0](
      vptr, pre_stride, dst_v, dst_uvstride,
      xd->subpix.filter_x[_o16x16mv.as_mv.col & 15], xd->subpix.x_step_q4,
      xd->subpix.filter_y[_o16x16mv.as_mv.row & 15], xd->subpix.y_step_q4,
      8, 8);
}


void vp9_build_1st_inter16x16_predictors_mb(MACROBLOCKD *xd,
                                            uint8_t *dst_y,
                                            uint8_t *dst_u,
                                            uint8_t *dst_v,
                                            int dst_ystride, int dst_uvstride) {
  vp9_build_1st_inter16x16_predictors_mby(xd, dst_y, dst_ystride,
      xd->mode_info_context->mbmi.need_to_clamp_mvs);
  vp9_build_1st_inter16x16_predictors_mbuv(xd, dst_u, dst_v, dst_uvstride);
}

void vp9_build_inter32x32_predictors_sb(MACROBLOCKD *x,
                                        uint8_t *dst_y,
                                        uint8_t *dst_u,
                                        uint8_t *dst_v,
                                        int dst_ystride,
                                        int dst_uvstride) {
  uint8_t *y1 = x->pre.y_buffer, *u1 = x->pre.u_buffer, *v1 = x->pre.v_buffer;
  uint8_t *y2 = x->second_pre.y_buffer, *u2 = x->second_pre.u_buffer,
          *v2 = x->second_pre.v_buffer;
  int edge[4], n;

  edge[0] = x->mb_to_top_edge;
  edge[1] = x->mb_to_bottom_edge;
  edge[2] = x->mb_to_left_edge;
  edge[3] = x->mb_to_right_edge;

  for (n = 0; n < 4; n++) {
    const int x_idx = n & 1, y_idx = n >> 1;

    x->mb_to_top_edge    = edge[0] -      ((y_idx  * 16) << 3);
    x->mb_to_bottom_edge = edge[1] + (((1 - y_idx) * 16) << 3);
    x->mb_to_left_edge   = edge[2] -      ((x_idx  * 16) << 3);
    x->mb_to_right_edge  = edge[3] + (((1 - x_idx) * 16) << 3);

    x->pre.y_buffer = y1 + y_idx * 16 * x->pre.y_stride  + x_idx * 16;
    x->pre.u_buffer = u1 + y_idx *  8 * x->pre.uv_stride + x_idx *  8;
    x->pre.v_buffer = v1 + y_idx *  8 * x->pre.uv_stride + x_idx *  8;

    vp9_build_1st_inter16x16_predictors_mb(x,
      dst_y + y_idx * 16 * dst_ystride  + x_idx * 16,
      dst_u + y_idx *  8 * dst_uvstride + x_idx *  8,
      dst_v + y_idx *  8 * dst_uvstride + x_idx *  8,
      dst_ystride, dst_uvstride);
    if (x->mode_info_context->mbmi.second_ref_frame > 0) {
      x->second_pre.y_buffer = y2 + y_idx * 16 * x->pre.y_stride  + x_idx * 16;
      x->second_pre.u_buffer = u2 + y_idx *  8 * x->pre.uv_stride + x_idx *  8;
      x->second_pre.v_buffer = v2 + y_idx *  8 * x->pre.uv_stride + x_idx *  8;

      vp9_build_2nd_inter16x16_predictors_mb(x,
        dst_y + y_idx * 16 * dst_ystride  + x_idx * 16,
        dst_u + y_idx *  8 * dst_uvstride + x_idx *  8,
        dst_v + y_idx *  8 * dst_uvstride + x_idx *  8,
        dst_ystride, dst_uvstride);
    }
  }

  x->mb_to_top_edge    = edge[0];
  x->mb_to_bottom_edge = edge[1];
  x->mb_to_left_edge   = edge[2];
  x->mb_to_right_edge  = edge[3];

  x->pre.y_buffer = y1;
  x->pre.u_buffer = u1;
  x->pre.v_buffer = v1;

  if (x->mode_info_context->mbmi.second_ref_frame > 0) {
    x->second_pre.y_buffer = y2;
    x->second_pre.u_buffer = u2;
    x->second_pre.v_buffer = v2;
  }

#if CONFIG_COMP_INTERINTRA_PRED
  if (x->mode_info_context->mbmi.second_ref_frame == INTRA_FRAME) {
    vp9_build_interintra_32x32_predictors_sb(
        x, dst_y, dst_u, dst_v, dst_ystride, dst_uvstride);
  }
#endif
}

void vp9_build_inter64x64_predictors_sb(MACROBLOCKD *x,
                                        uint8_t *dst_y,
                                        uint8_t *dst_u,
                                        uint8_t *dst_v,
                                        int dst_ystride,
                                        int dst_uvstride) {
  uint8_t *y1 = x->pre.y_buffer, *u1 = x->pre.u_buffer, *v1 = x->pre.v_buffer;
  uint8_t *y2 = x->second_pre.y_buffer, *u2 = x->second_pre.u_buffer,
          *v2 = x->second_pre.v_buffer;
  int edge[4], n;

  edge[0] = x->mb_to_top_edge;
  edge[1] = x->mb_to_bottom_edge;
  edge[2] = x->mb_to_left_edge;
  edge[3] = x->mb_to_right_edge;

  for (n = 0; n < 4; n++) {
    const int x_idx = n & 1, y_idx = n >> 1;

    x->mb_to_top_edge    = edge[0] -      ((y_idx  * 32) << 3);
    x->mb_to_bottom_edge = edge[1] + (((1 - y_idx) * 32) << 3);
    x->mb_to_left_edge   = edge[2] -      ((x_idx  * 32) << 3);
    x->mb_to_right_edge  = edge[3] + (((1 - x_idx) * 32) << 3);

    x->pre.y_buffer = y1 + y_idx * 32 * x->pre.y_stride  + x_idx * 32;
    x->pre.u_buffer = u1 + y_idx * 16 * x->pre.uv_stride + x_idx * 16;
    x->pre.v_buffer = v1 + y_idx * 16 * x->pre.uv_stride + x_idx * 16;

    if (x->mode_info_context->mbmi.second_ref_frame > 0) {
      x->second_pre.y_buffer = y2 + y_idx * 32 * x->pre.y_stride  + x_idx * 32;
      x->second_pre.u_buffer = u2 + y_idx * 16 * x->pre.uv_stride + x_idx * 16;
      x->second_pre.v_buffer = v2 + y_idx * 16 * x->pre.uv_stride + x_idx * 16;
    }

    vp9_build_inter32x32_predictors_sb(x,
        dst_y + y_idx * 32 * dst_ystride  + x_idx * 32,
        dst_u + y_idx * 16 * dst_uvstride + x_idx * 16,
        dst_v + y_idx * 16 * dst_uvstride + x_idx * 16,
        dst_ystride, dst_uvstride);
  }

  x->mb_to_top_edge    = edge[0];
  x->mb_to_bottom_edge = edge[1];
  x->mb_to_left_edge   = edge[2];
  x->mb_to_right_edge  = edge[3];

  x->pre.y_buffer = y1;
  x->pre.u_buffer = u1;
  x->pre.v_buffer = v1;

  if (x->mode_info_context->mbmi.second_ref_frame > 0) {
    x->second_pre.y_buffer = y2;
    x->second_pre.u_buffer = u2;
    x->second_pre.v_buffer = v2;
  }

#if CONFIG_COMP_INTERINTRA_PRED
  if (x->mode_info_context->mbmi.second_ref_frame == INTRA_FRAME) {
    vp9_build_interintra_64x64_predictors_sb(x, dst_y, dst_u, dst_v,
                                             dst_ystride, dst_uvstride);
  }
#endif
}

/*
 * The following functions should be called after an initial
 * call to vp9_build_1st_inter16x16_predictors_mb() or _mby()/_mbuv().
 * It will run a second filter on a (different) ref
 * frame and average the result with the output of the
 * first filter. The second reference frame is stored
 * in x->second_pre (the reference frame index is in
 * x->mode_info_context->mbmi.second_ref_frame). The second
 * motion vector is x->mode_info_context->mbmi.second_mv.
 *
 * This allows blending prediction from two reference frames
 * which sometimes leads to better prediction than from a
 * single reference framer.
 */
void vp9_build_2nd_inter16x16_predictors_mby(MACROBLOCKD *xd,
                                             uint8_t *dst_y,
                                             int dst_ystride) {
  uint8_t *ptr;

  int_mv _16x16mv;
  int mv_row;
  int mv_col;

  uint8_t *ptr_base = xd->second_pre.y_buffer;
  int pre_stride = xd->block[0].pre_stride;

  _16x16mv.as_int = xd->mode_info_context->mbmi.mv[1].as_int;

  if (xd->mode_info_context->mbmi.need_to_clamp_secondmv)
    clamp_mv_to_umv_border(&_16x16mv.as_mv, xd);

  mv_row = _16x16mv.as_mv.row;
  mv_col = _16x16mv.as_mv.col;

  ptr = ptr_base + (mv_row >> 3) * pre_stride + (mv_col >> 3);

  xd->subpix.predict[!!(mv_col & 7)][!!(mv_row & 7)][1](
      ptr, pre_stride, dst_y, dst_ystride,
      xd->subpix.filter_x[(mv_col & 7) << 1], xd->subpix.x_step_q4,
      xd->subpix.filter_y[(mv_row & 7) << 1], xd->subpix.y_step_q4,
      16, 16);
}

void vp9_build_2nd_inter16x16_predictors_mbuv(MACROBLOCKD *xd,
                                              uint8_t *dst_u,
                                              uint8_t *dst_v,
                                              int dst_uvstride) {
  int offset;
  uint8_t *uptr, *vptr;

  int_mv _16x16mv;
  int mv_row;
  int mv_col;
  int omv_row, omv_col;

  int pre_stride = xd->block[0].pre_stride;

  _16x16mv.as_int = xd->mode_info_context->mbmi.mv[1].as_int;

  if (xd->mode_info_context->mbmi.need_to_clamp_secondmv)
    clamp_mv_to_umv_border(&_16x16mv.as_mv, xd);

  mv_row = _16x16mv.as_mv.row;
  mv_col = _16x16mv.as_mv.col;

  /* calc uv motion vectors */
  omv_row = mv_row;
  omv_col = mv_col;
  mv_row = (mv_row + (mv_row > 0)) >> 1;
  mv_col = (mv_col + (mv_col > 0)) >> 1;

  mv_row &= xd->fullpixel_mask;
  mv_col &= xd->fullpixel_mask;

  pre_stride >>= 1;
  offset = (mv_row >> 3) * pre_stride + (mv_col >> 3);
  uptr = xd->second_pre.u_buffer + offset;
  vptr = xd->second_pre.v_buffer + offset;

  xd->subpix.predict[!!(omv_col & 15)][!!(omv_row & 15)][1](
      uptr, pre_stride, dst_u, dst_uvstride,
      xd->subpix.filter_x[omv_col & 15], xd->subpix.x_step_q4,
      xd->subpix.filter_y[omv_row & 15], xd->subpix.y_step_q4,
      8, 8);

  xd->subpix.predict[!!(omv_col & 15)][!!(omv_row & 15)][1](
      vptr, pre_stride, dst_v, dst_uvstride,
      xd->subpix.filter_x[omv_col & 15], xd->subpix.x_step_q4,
      xd->subpix.filter_y[omv_row & 15], xd->subpix.y_step_q4,
      8, 8);
}

void vp9_build_2nd_inter16x16_predictors_mb(MACROBLOCKD *xd,
                                            uint8_t *dst_y,
                                            uint8_t *dst_u,
                                            uint8_t *dst_v,
                                            int dst_ystride,
                                            int dst_uvstride) {
  vp9_build_2nd_inter16x16_predictors_mby(xd, dst_y, dst_ystride);
  vp9_build_2nd_inter16x16_predictors_mbuv(xd, dst_u, dst_v, dst_uvstride);
}

static void build_inter4x4_predictors_mb(MACROBLOCKD *xd) {
  int i;
  MB_MODE_INFO * mbmi = &xd->mode_info_context->mbmi;
  BLOCKD *blockd = xd->block;

  if (xd->mode_info_context->mbmi.partitioning != PARTITIONING_4X4) {
    blockd[ 0].bmi = xd->mode_info_context->bmi[ 0];
    blockd[ 2].bmi = xd->mode_info_context->bmi[ 2];
    blockd[ 8].bmi = xd->mode_info_context->bmi[ 8];
    blockd[10].bmi = xd->mode_info_context->bmi[10];

    if (mbmi->need_to_clamp_mvs) {
      clamp_mv_to_umv_border(&blockd[ 0].bmi.as_mv.first.as_mv, xd);
      clamp_mv_to_umv_border(&blockd[ 2].bmi.as_mv.first.as_mv, xd);
      clamp_mv_to_umv_border(&blockd[ 8].bmi.as_mv.first.as_mv, xd);
      clamp_mv_to_umv_border(&blockd[10].bmi.as_mv.first.as_mv, xd);
      if (mbmi->second_ref_frame > 0) {
        clamp_mv_to_umv_border(&blockd[ 0].bmi.as_mv.second.as_mv, xd);
        clamp_mv_to_umv_border(&blockd[ 2].bmi.as_mv.second.as_mv, xd);
        clamp_mv_to_umv_border(&blockd[ 8].bmi.as_mv.second.as_mv, xd);
        clamp_mv_to_umv_border(&blockd[10].bmi.as_mv.second.as_mv, xd);
      }
    }


    vp9_build_inter_predictors4b(xd, &blockd[ 0], 16);
    vp9_build_inter_predictors4b(xd, &blockd[ 2], 16);
    vp9_build_inter_predictors4b(xd, &blockd[ 8], 16);
    vp9_build_inter_predictors4b(xd, &blockd[10], 16);

    if (mbmi->second_ref_frame > 0) {
      vp9_build_2nd_inter_predictors4b(xd, &blockd[ 0], 16);
      vp9_build_2nd_inter_predictors4b(xd, &blockd[ 2], 16);
      vp9_build_2nd_inter_predictors4b(xd, &blockd[ 8], 16);
      vp9_build_2nd_inter_predictors4b(xd, &blockd[10], 16);
    }
  } else {
    for (i = 0; i < 16; i += 2) {
      BLOCKD *d0 = &blockd[i];
      BLOCKD *d1 = &blockd[i + 1];

      blockd[i + 0].bmi = xd->mode_info_context->bmi[i + 0];
      blockd[i + 1].bmi = xd->mode_info_context->bmi[i + 1];

      if (mbmi->need_to_clamp_mvs) {
        clamp_mv_to_umv_border(&blockd[i + 0].bmi.as_mv.first.as_mv, xd);
        clamp_mv_to_umv_border(&blockd[i + 1].bmi.as_mv.first.as_mv, xd);
        if (mbmi->second_ref_frame > 0) {
          clamp_mv_to_umv_border(&blockd[i + 0].bmi.as_mv.second.as_mv, xd);
          clamp_mv_to_umv_border(&blockd[i + 1].bmi.as_mv.second.as_mv, xd);
        }
      }

      if (d0->bmi.as_mv.first.as_int == d1->bmi.as_mv.first.as_int)
        build_inter_predictors2b(xd, d0, 16);
      else {
        vp9_build_inter_predictors_b(d0, 16, &xd->subpix);
        vp9_build_inter_predictors_b(d1, 16, &xd->subpix);
      }

      if (mbmi->second_ref_frame > 0) {
        vp9_build_2nd_inter_predictors_b(d0, 16, &xd->subpix);
        vp9_build_2nd_inter_predictors_b(d1, 16, &xd->subpix);
      }
    }
  }

  for (i = 16; i < 24; i += 2) {
    BLOCKD *d0 = &blockd[i];
    BLOCKD *d1 = &blockd[i + 1];

    if (d0->bmi.as_mv.first.as_int == d1->bmi.as_mv.first.as_int)
      build_inter_predictors2b(xd, d0, 8);
    else {
      vp9_build_inter_predictors_b(d0, 8, &xd->subpix);
      vp9_build_inter_predictors_b(d1, 8, &xd->subpix);
    }

    if (mbmi->second_ref_frame > 0) {
      vp9_build_2nd_inter_predictors_b(d0, 8, &xd->subpix);
      vp9_build_2nd_inter_predictors_b(d1, 8, &xd->subpix);
    }
  }
}

static
void build_4x4uvmvs(MACROBLOCKD *xd) {
  int i, j;
  BLOCKD *blockd = xd->block;

  for (i = 0; i < 2; i++) {
    for (j = 0; j < 2; j++) {
      int yoffset = i * 8 + j * 2;
      int uoffset = 16 + i * 2 + j;
      int voffset = 20 + i * 2 + j;

      int temp;

      temp = xd->mode_info_context->bmi[yoffset + 0].as_mv.first.as_mv.row
             + xd->mode_info_context->bmi[yoffset + 1].as_mv.first.as_mv.row
             + xd->mode_info_context->bmi[yoffset + 4].as_mv.first.as_mv.row
             + xd->mode_info_context->bmi[yoffset + 5].as_mv.first.as_mv.row;

      if (temp < 0) temp -= 4;
      else temp += 4;

      blockd[uoffset].bmi.as_mv.first.as_mv.row = (temp / 8) &
                                                  xd->fullpixel_mask;

      temp = xd->mode_info_context->bmi[yoffset + 0].as_mv.first.as_mv.col
             + xd->mode_info_context->bmi[yoffset + 1].as_mv.first.as_mv.col
             + xd->mode_info_context->bmi[yoffset + 4].as_mv.first.as_mv.col
             + xd->mode_info_context->bmi[yoffset + 5].as_mv.first.as_mv.col;

      if (temp < 0) temp -= 4;
      else temp += 4;

      blockd[uoffset].bmi.as_mv.first.as_mv.col = (temp / 8) &
        xd->fullpixel_mask;

      // if (x->mode_info_context->mbmi.need_to_clamp_mvs)
      clamp_uvmv_to_umv_border(&blockd[uoffset].bmi.as_mv.first.as_mv, xd);

      // if (x->mode_info_context->mbmi.need_to_clamp_mvs)
      clamp_uvmv_to_umv_border(&blockd[uoffset].bmi.as_mv.first.as_mv, xd);

      blockd[voffset].bmi.as_mv.first.as_mv.row =
        blockd[uoffset].bmi.as_mv.first.as_mv.row;
      blockd[voffset].bmi.as_mv.first.as_mv.col =
        blockd[uoffset].bmi.as_mv.first.as_mv.col;

      if (xd->mode_info_context->mbmi.second_ref_frame > 0) {
        temp = xd->mode_info_context->bmi[yoffset + 0].as_mv.second.as_mv.row
               + xd->mode_info_context->bmi[yoffset + 1].as_mv.second.as_mv.row
               + xd->mode_info_context->bmi[yoffset + 4].as_mv.second.as_mv.row
               + xd->mode_info_context->bmi[yoffset + 5].as_mv.second.as_mv.row;

        if (temp < 0) {
          temp -= 4;
        } else {
          temp += 4;
        }

       blockd[uoffset].bmi.as_mv.second.as_mv.row = (temp / 8) &
                                                    xd->fullpixel_mask;

        temp = xd->mode_info_context->bmi[yoffset + 0].as_mv.second.as_mv.col
               + xd->mode_info_context->bmi[yoffset + 1].as_mv.second.as_mv.col
               + xd->mode_info_context->bmi[yoffset + 4].as_mv.second.as_mv.col
               + xd->mode_info_context->bmi[yoffset + 5].as_mv.second.as_mv.col;

        if (temp < 0) {
          temp -= 4;
        } else {
          temp += 4;
        }

        blockd[uoffset].bmi.as_mv.second.as_mv.col = (temp / 8) &
                                                        xd->fullpixel_mask;

        // if (mbmi->need_to_clamp_mvs)
        clamp_uvmv_to_umv_border(
          &blockd[uoffset].bmi.as_mv.second.as_mv, xd);

        // if (mbmi->need_to_clamp_mvs)
        clamp_uvmv_to_umv_border(
          &blockd[uoffset].bmi.as_mv.second.as_mv, xd);

        blockd[voffset].bmi.as_mv.second.as_mv.row =
          blockd[uoffset].bmi.as_mv.second.as_mv.row;
        blockd[voffset].bmi.as_mv.second.as_mv.col =
          blockd[uoffset].bmi.as_mv.second.as_mv.col;
      }
    }
  }
}

void vp9_build_inter_predictors_mb(MACROBLOCKD *xd) {
  if (xd->mode_info_context->mbmi.mode != SPLITMV) {
    vp9_build_1st_inter16x16_predictors_mb(xd, xd->predictor,
                                           &xd->predictor[256],
                                           &xd->predictor[320], 16, 8);

    if (xd->mode_info_context->mbmi.second_ref_frame > 0) {
      /* 256 = offset of U plane in Y+U+V buffer;
       * 320 = offset of V plane in Y+U+V buffer.
       * (256=16x16, 320=16x16+8x8). */
      vp9_build_2nd_inter16x16_predictors_mb(xd, xd->predictor,
                                             &xd->predictor[256],
                                             &xd->predictor[320], 16, 8);
    }
#if CONFIG_COMP_INTERINTRA_PRED
    else if (xd->mode_info_context->mbmi.second_ref_frame == INTRA_FRAME) {
      vp9_build_interintra_16x16_predictors_mb(xd, xd->predictor,
                                               &xd->predictor[256],
                                               &xd->predictor[320], 16, 8);
    }
#endif
  } else {
    build_4x4uvmvs(xd);
    build_inter4x4_predictors_mb(xd);
  }
}
