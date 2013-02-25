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

void vp9_setup_scale_factors_for_frame(struct scale_factors *scale,
                                       YV12_BUFFER_CONFIG *other,
                                       int this_w, int this_h) {
  int other_w, other_h;

  other_h = other->y_height;
  other_w = other->y_width;
  scale->x_num = other_w;
  scale->x_den = this_w;
  scale->x_offset_q4 = 0;  // calculated per-mb
  scale->x_step_q4 = 16 * other_w / this_w;
  scale->y_num = other_h;
  scale->y_den = this_h;
  scale->y_offset_q4 = 0;  // calculated per-mb
  scale->y_step_q4 = 16 * other_h / this_h;

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
  int i;

  /* Calculate scaling factors for each of the 3 available references */
  for (i = 0; i < 3; ++i) {
    if (cm->active_ref_idx[i] >= NUM_YV12_BUFFERS) {
      memset(&cm->active_ref_scale[i], 0, sizeof(cm->active_ref_scale[i]));
      continue;
    }

    vp9_setup_scale_factors_for_frame(&cm->active_ref_scale[i],
                                      &cm->yv12_fb[cm->active_ref_idx[i]],
                                      cm->mb_cols * 16, cm->mb_rows * 16);
  }

  if (xd->mode_info_context) {
    MB_MODE_INFO *mbmi = &xd->mode_info_context->mbmi;

    set_scale_factors(xd,
                      mbmi->ref_frame - 1,
                      mbmi->second_ref_frame - 1,
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

static void set_scaled_offsets(struct scale_factors *scale,
                               int row, int col) {
  const int x_q4 = 16 * col;
  const int y_q4 = 16 * row;

  scale->x_offset_q4 = (x_q4 * scale->x_num / scale->x_den) & 0xf;
  scale->y_offset_q4 = (y_q4 * scale->y_num / scale->y_den) & 0xf;
}

static int32_t scale_motion_vector_component_q3(int mv_q3,
                                                int num,
                                                int den,
                                                int offset_q4) {
  // returns the scaled and offset value of the mv component.
  const int32_t mv_q4 = mv_q3 << 1;

  /* TODO(jkoleszar): make fixed point, or as a second multiply? */
  return mv_q4 * num / den + offset_q4;
}

static int32_t scale_motion_vector_component_q4(int mv_q4,
                                                int num,
                                                int den,
                                                int offset_q4) {
  // returns the scaled and offset value of the mv component.

  /* TODO(jkoleszar): make fixed point, or as a second multiply? */
  return mv_q4 * num / den + offset_q4;
}

static int_mv32 scale_motion_vector_q3_to_q4(
    const int_mv *src_mv,
    const struct scale_factors *scale) {
  // returns mv * scale + offset
  int_mv32 result;

  result.as_mv.row = scale_motion_vector_component_q3(src_mv->as_mv.row,
                                                      scale->y_num,
                                                      scale->y_den,
                                                      scale->y_offset_q4);
  result.as_mv.col = scale_motion_vector_component_q3(src_mv->as_mv.col,
                                                      scale->x_num,
                                                      scale->x_den,
                                                      scale->x_offset_q4);
  return result;
}

void vp9_build_inter_predictor(const uint8_t *src, int src_stride,
                               uint8_t *dst, int dst_stride,
                               const int_mv *mv_q3,
                               const struct scale_factors *scale,
                               int w, int h, int do_avg,
                               const struct subpix_fn_table *subpix) {
  int_mv32 mv;

  mv = scale_motion_vector_q3_to_q4(mv_q3, scale);
  src = src + (mv.as_mv.row >> 4) * src_stride + (mv.as_mv.col >> 4);

  scale->predict[!!(mv.as_mv.col & 15)][!!(mv.as_mv.row & 15)][do_avg](
      src, src_stride, dst, dst_stride,
      subpix->filter_x[mv.as_mv.col & 15], scale->x_step_q4,
      subpix->filter_y[mv.as_mv.row & 15], scale->y_step_q4,
      w, h);
}

/* Like vp9_build_inter_predictor, but takes the full-pel part of the
 * mv separately, and the fractional part as a q4.
 */
void vp9_build_inter_predictor_q4(const uint8_t *src, int src_stride,
                                  uint8_t *dst, int dst_stride,
                                  const int_mv *fullpel_mv_q3,
                                  const int_mv *frac_mv_q4,
                                  const struct scale_factors *scale,
                                  int w, int h, int do_avg,
                                  const struct subpix_fn_table *subpix) {
  const int mv_row_q4 = ((fullpel_mv_q3->as_mv.row >> 3) << 4)
                        + (frac_mv_q4->as_mv.row & 0xf);
  const int mv_col_q4 = ((fullpel_mv_q3->as_mv.col >> 3) << 4)
                        + (frac_mv_q4->as_mv.col & 0xf);
  const int scaled_mv_row_q4 =
      scale_motion_vector_component_q4(mv_row_q4, scale->y_num, scale->y_den,
                                       scale->y_offset_q4);
  const int scaled_mv_col_q4 =
      scale_motion_vector_component_q4(mv_col_q4, scale->x_num, scale->x_den,
                                       scale->x_offset_q4);
  const int subpel_x = scaled_mv_col_q4 & 15;
  const int subpel_y = scaled_mv_row_q4 & 15;

  src = src + (scaled_mv_row_q4 >> 4) * src_stride + (scaled_mv_col_q4 >> 4);
  scale->predict[!!subpel_x][!!subpel_y][do_avg](
      src, src_stride, dst, dst_stride,
      subpix->filter_x[subpel_x], scale->x_step_q4,
      subpix->filter_y[subpel_y], scale->y_step_q4,
      w, h);
}

static void build_2x1_inter_predictor(const BLOCKD *d0, const BLOCKD *d1,
                                      const struct scale_factors *scale,
                                      int block_size, int stride, int which_mv,
                                      const struct subpix_fn_table *subpix) {
  assert(d1->predictor - d0->predictor == block_size);
  assert(d1->pre == d0->pre + block_size);

  if (d0->bmi.as_mv[which_mv].as_int == d1->bmi.as_mv[which_mv].as_int) {
    uint8_t **base_pre = which_mv ? d0->base_second_pre : d0->base_pre;

    vp9_build_inter_predictor(*base_pre + d0->pre,
                              d0->pre_stride,
                              d0->predictor, stride,
                              &d0->bmi.as_mv[which_mv],
                              &scale[which_mv],
                              2 * block_size, block_size, which_mv,
                              subpix);

  } else {
    uint8_t **base_pre0 = which_mv ? d0->base_second_pre : d0->base_pre;
    uint8_t **base_pre1 = which_mv ? d1->base_second_pre : d1->base_pre;

    vp9_build_inter_predictor(*base_pre0 + d0->pre,
                              d0->pre_stride,
                              d0->predictor, stride,
                              &d0->bmi.as_mv[which_mv],
                              &scale[which_mv],
                              block_size, block_size, which_mv,
                              subpix);
    vp9_build_inter_predictor(*base_pre1 + d1->pre,
                              d1->pre_stride,
                              d1->predictor, stride,
                              &d1->bmi.as_mv[which_mv],
                              &scale[which_mv],
                              block_size, block_size, which_mv,
                              subpix);
  }
}

/*encoder only*/
void vp9_build_inter4x4_predictors_mbuv(MACROBLOCKD *xd,
                                        int mb_row,
                                        int mb_col) {
  int i, j;
  BLOCKD *blockd = xd->block;

  /* build uv mvs */
  for (i = 0; i < 2; i++) {
    for (j = 0; j < 2; j++) {
      int yoffset = i * 8 + j * 2;
      int uoffset = 16 + i * 2 + j;
      int voffset = 20 + i * 2 + j;
      int temp;

      temp = blockd[yoffset  ].bmi.as_mv[0].as_mv.row
             + blockd[yoffset + 1].bmi.as_mv[0].as_mv.row
             + blockd[yoffset + 4].bmi.as_mv[0].as_mv.row
             + blockd[yoffset + 5].bmi.as_mv[0].as_mv.row;

      if (temp < 0) temp -= 4;
      else temp += 4;

      xd->block[uoffset].bmi.as_mv[0].as_mv.row = (temp / 8) &
        xd->fullpixel_mask;

      temp = blockd[yoffset  ].bmi.as_mv[0].as_mv.col
             + blockd[yoffset + 1].bmi.as_mv[0].as_mv.col
             + blockd[yoffset + 4].bmi.as_mv[0].as_mv.col
             + blockd[yoffset + 5].bmi.as_mv[0].as_mv.col;

      if (temp < 0) temp -= 4;
      else temp += 4;

      blockd[uoffset].bmi.as_mv[0].as_mv.col = (temp / 8) &
        xd->fullpixel_mask;

      blockd[voffset].bmi.as_mv[0].as_mv.row =
        blockd[uoffset].bmi.as_mv[0].as_mv.row;
      blockd[voffset].bmi.as_mv[0].as_mv.col =
        blockd[uoffset].bmi.as_mv[0].as_mv.col;

      if (xd->mode_info_context->mbmi.second_ref_frame > 0) {
        temp = blockd[yoffset  ].bmi.as_mv[1].as_mv.row
               + blockd[yoffset + 1].bmi.as_mv[1].as_mv.row
               + blockd[yoffset + 4].bmi.as_mv[1].as_mv.row
               + blockd[yoffset + 5].bmi.as_mv[1].as_mv.row;

        if (temp < 0) {
          temp -= 4;
        } else {
          temp += 4;
        }

        blockd[uoffset].bmi.as_mv[1].as_mv.row = (temp / 8) &
          xd->fullpixel_mask;

        temp = blockd[yoffset  ].bmi.as_mv[1].as_mv.col
               + blockd[yoffset + 1].bmi.as_mv[1].as_mv.col
               + blockd[yoffset + 4].bmi.as_mv[1].as_mv.col
               + blockd[yoffset + 5].bmi.as_mv[1].as_mv.col;

        if (temp < 0) {
          temp -= 4;
        } else {
          temp += 4;
        }

        blockd[uoffset].bmi.as_mv[1].as_mv.col = (temp / 8) &
          xd->fullpixel_mask;

        blockd[voffset].bmi.as_mv[1].as_mv.row =
          blockd[uoffset].bmi.as_mv[1].as_mv.row;
        blockd[voffset].bmi.as_mv[1].as_mv.col =
          blockd[uoffset].bmi.as_mv[1].as_mv.col;
      }
    }
  }

  for (i = 16; i < 24; i += 2) {
    const int use_second_ref = xd->mode_info_context->mbmi.second_ref_frame > 0;
    const int x = 4 * (i & 1);
    const int y = ((i - 16) >> 1) * 4;

    int which_mv;
    BLOCKD *d0 = &blockd[i];
    BLOCKD *d1 = &blockd[i + 1];

    for (which_mv = 0; which_mv < 1 + use_second_ref; ++which_mv) {
      set_scaled_offsets(&xd->scale_factor_uv[which_mv],
                         mb_row * 8 + y, mb_col * 8 + x);

      build_2x1_inter_predictor(d0, d1, xd->scale_factor_uv, 4, 8, which_mv,
                                &xd->subpix);
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
void vp9_build_inter16x16_predictors_mby(MACROBLOCKD *xd,
                                         uint8_t *dst_y,
                                         int dst_ystride,
                                         int mb_row,
                                         int mb_col) {
  const int use_second_ref = xd->mode_info_context->mbmi.second_ref_frame > 0;
  int which_mv;

  for (which_mv = 0; which_mv < 1 + use_second_ref; ++which_mv) {
    const int clamp_mvs =
        which_mv ? xd->mode_info_context->mbmi.need_to_clamp_secondmv
                 : xd->mode_info_context->mbmi.need_to_clamp_mvs;
    uint8_t *base_pre;
    int_mv ymv;
    int pre_stride;

    ymv.as_int = xd->mode_info_context->mbmi.mv[which_mv].as_int;
    base_pre = which_mv ? xd->second_pre.y_buffer
                        : xd->pre.y_buffer;
    pre_stride = which_mv ? xd->second_pre.y_stride
                          : xd->pre.y_stride;
    if (clamp_mvs)
      clamp_mv_to_umv_border(&ymv.as_mv, xd);

    set_scaled_offsets(&xd->scale_factor[which_mv], mb_row * 16, mb_col * 16);

    vp9_build_inter_predictor(base_pre, pre_stride,
                              dst_y, dst_ystride,
                              &ymv, &xd->scale_factor[which_mv],
                              16, 16, which_mv, &xd->subpix);
  }
}

void vp9_build_inter16x16_predictors_mbuv(MACROBLOCKD *xd,
                                          uint8_t *dst_u,
                                          uint8_t *dst_v,
                                          int dst_uvstride,
                                          int mb_row,
                                          int mb_col) {
  const int use_second_ref = xd->mode_info_context->mbmi.second_ref_frame > 0;
  int which_mv;

  for (which_mv = 0; which_mv < 1 + use_second_ref; ++which_mv) {
    const int clamp_mvs =
        which_mv ? xd->mode_info_context->mbmi.need_to_clamp_secondmv
                 : xd->mode_info_context->mbmi.need_to_clamp_mvs;
    uint8_t *uptr, *vptr;
    int pre_stride = which_mv ? xd->second_pre.y_stride
                              : xd->pre.y_stride;
    int_mv _o16x16mv;
    int_mv _16x16mv;

    _16x16mv.as_int = xd->mode_info_context->mbmi.mv[which_mv].as_int;

    if (clamp_mvs)
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
    uptr = (which_mv ? xd->second_pre.u_buffer : xd->pre.u_buffer);
    vptr = (which_mv ? xd->second_pre.v_buffer : xd->pre.v_buffer);

    set_scaled_offsets(&xd->scale_factor_uv[which_mv],
                       mb_row * 16, mb_col * 16);

    vp9_build_inter_predictor_q4(uptr, pre_stride,
                                 dst_u, dst_uvstride,
                                 &_16x16mv, &_o16x16mv,
                                 &xd->scale_factor_uv[which_mv],
                                 8, 8, which_mv, &xd->subpix);

    vp9_build_inter_predictor_q4(vptr, pre_stride,
                                 dst_v, dst_uvstride,
                                 &_16x16mv, &_o16x16mv,
                                 &xd->scale_factor_uv[which_mv],
                                 8, 8, which_mv, &xd->subpix);
  }
}

void vp9_build_inter32x32_predictors_sb(MACROBLOCKD *x,
                                        uint8_t *dst_y,
                                        uint8_t *dst_u,
                                        uint8_t *dst_v,
                                        int dst_ystride,
                                        int dst_uvstride,
                                        int mb_row,
                                        int mb_col) {
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
    int scaled_uv_offset;

    x->mb_to_top_edge    = edge[0] -      ((y_idx  * 16) << 3);
    x->mb_to_bottom_edge = edge[1] + (((1 - y_idx) * 16) << 3);
    x->mb_to_left_edge   = edge[2] -      ((x_idx  * 16) << 3);
    x->mb_to_right_edge  = edge[3] + (((1 - x_idx) * 16) << 3);

    x->pre.y_buffer = y1 + scaled_buffer_offset(x_idx * 16,
                                                y_idx * 16,
                                                x->pre.y_stride,
                                                &x->scale_factor[0]);
    scaled_uv_offset = scaled_buffer_offset(x_idx * 8,
                                            y_idx * 8,
                                            x->pre.uv_stride,
                                            &x->scale_factor_uv[0]);
    x->pre.u_buffer = u1 + scaled_uv_offset;
    x->pre.v_buffer = v1 + scaled_uv_offset;

    if (x->mode_info_context->mbmi.second_ref_frame > 0) {
      x->second_pre.y_buffer = y2 +
          scaled_buffer_offset(x_idx * 16,
                               y_idx * 16,
                               x->second_pre.y_stride,
                               &x->scale_factor[1]);
      scaled_uv_offset = scaled_buffer_offset(x_idx * 8,
                                              y_idx * 8,
                                              x->second_pre.uv_stride,
                                              &x->scale_factor_uv[1]);
      x->second_pre.u_buffer = u2 + scaled_uv_offset;
      x->second_pre.v_buffer = v2 + scaled_uv_offset;
    }

    vp9_build_inter16x16_predictors_mb(x,
        dst_y + y_idx * 16 * dst_ystride  + x_idx * 16,
        dst_u + y_idx *  8 * dst_uvstride + x_idx *  8,
        dst_v + y_idx *  8 * dst_uvstride + x_idx *  8,
        dst_ystride, dst_uvstride, mb_row + y_idx, mb_col + x_idx);
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
                                        int dst_uvstride,
                                        int mb_row,
                                        int mb_col) {
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
    int scaled_uv_offset;

    x->mb_to_top_edge    = edge[0] -      ((y_idx  * 32) << 3);
    x->mb_to_bottom_edge = edge[1] + (((1 - y_idx) * 32) << 3);
    x->mb_to_left_edge   = edge[2] -      ((x_idx  * 32) << 3);
    x->mb_to_right_edge  = edge[3] + (((1 - x_idx) * 32) << 3);

    x->pre.y_buffer = y1 + scaled_buffer_offset(x_idx * 32,
                                                y_idx * 32,
                                                x->pre.y_stride,
                                                &x->scale_factor[0]);
    scaled_uv_offset = scaled_buffer_offset(x_idx * 16,
                                            y_idx * 16,
                                            x->pre.uv_stride,
                                            &x->scale_factor_uv[0]);
    x->pre.u_buffer = u1 + scaled_uv_offset;
    x->pre.v_buffer = v1 + scaled_uv_offset;

    if (x->mode_info_context->mbmi.second_ref_frame > 0) {
      x->second_pre.y_buffer = y2 +
          scaled_buffer_offset(x_idx * 32,
                               y_idx * 32,
                               x->second_pre.y_stride,
                               &x->scale_factor[1]);
      scaled_uv_offset = scaled_buffer_offset(x_idx * 16,
                                              y_idx * 16,
                                              x->second_pre.uv_stride,
                                              &x->scale_factor_uv[1]);
      x->second_pre.u_buffer = u2 + scaled_uv_offset;
      x->second_pre.v_buffer = v2 + scaled_uv_offset;
    }

    vp9_build_inter32x32_predictors_sb(x,
        dst_y + y_idx * 32 * dst_ystride  + x_idx * 32,
        dst_u + y_idx * 16 * dst_uvstride + x_idx * 16,
        dst_v + y_idx * 16 * dst_uvstride + x_idx * 16,
        dst_ystride, dst_uvstride, mb_row + y_idx * 2, mb_col + x_idx * 2);
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

static void build_inter4x4_predictors_mb(MACROBLOCKD *xd) {
  int i;
  MB_MODE_INFO * mbmi = &xd->mode_info_context->mbmi;
  BLOCKD *blockd = xd->block;
  int which_mv = 0;
  const int use_second_ref = mbmi->second_ref_frame > 0;

  if (xd->mode_info_context->mbmi.partitioning != PARTITIONING_4X4) {
    for (i = 0; i < 16; i += 8) {
      BLOCKD *d0 = &blockd[i];
      BLOCKD *d1 = &blockd[i + 2];

      blockd[i + 0].bmi = xd->mode_info_context->bmi[i + 0];
      blockd[i + 2].bmi = xd->mode_info_context->bmi[i + 2];

      for (which_mv = 0; which_mv < 1 + use_second_ref; ++which_mv) {
        if (mbmi->need_to_clamp_mvs) {
          clamp_mv_to_umv_border(&blockd[i + 0].bmi.as_mv[which_mv].as_mv, xd);
          clamp_mv_to_umv_border(&blockd[i + 2].bmi.as_mv[which_mv].as_mv, xd);
        }

        /* TODO(jkoleszar): Enabling this for EIGHTTAP_SMOOTH changes the
         * result slightly, for reasons that are not immediately obvious to me.
         * It probably makes sense to enable this for all filter types to be
         * consistent with the way we do 8x4 below. Leaving disabled for now.
         */
        if (mbmi->interp_filter != EIGHTTAP_SMOOTH) {
          build_2x1_inter_predictor(d0, d1, xd->scale_factor, 8, 16,
                                    which_mv, &xd->subpix);
        } else {
          uint8_t **base_pre0 = which_mv ? d0->base_second_pre : d0->base_pre;
          uint8_t **base_pre1 = which_mv ? d1->base_second_pre : d1->base_pre;

          vp9_build_inter_predictor(*base_pre0 + d0->pre,
                                    d0->pre_stride,
                                    d0->predictor, 16,
                                    &d0->bmi.as_mv[which_mv],
                                    &xd->scale_factor[which_mv],
                                    8, 8, which_mv, &xd->subpix);
          vp9_build_inter_predictor(*base_pre1 + d1->pre,
                                    d1->pre_stride,
                                    d1->predictor, 16,
                                    &d1->bmi.as_mv[which_mv],
                                    &xd->scale_factor[which_mv],
                                    8, 8, which_mv, &xd->subpix);
        }
      }
    }
  } else {
    for (i = 0; i < 16; i += 2) {
      BLOCKD *d0 = &blockd[i];
      BLOCKD *d1 = &blockd[i + 1];

      blockd[i + 0].bmi = xd->mode_info_context->bmi[i + 0];
      blockd[i + 1].bmi = xd->mode_info_context->bmi[i + 1];

      for (which_mv = 0; which_mv < 1 + use_second_ref; ++which_mv) {
        build_2x1_inter_predictor(d0, d1, xd->scale_factor, 4, 16,
                                  which_mv, &xd->subpix);
      }
    }
  }

  for (i = 16; i < 24; i += 2) {
    BLOCKD *d0 = &blockd[i];
    BLOCKD *d1 = &blockd[i + 1];

    for (which_mv = 0; which_mv < 1 + use_second_ref; ++which_mv) {
      build_2x1_inter_predictor(d0, d1, xd->scale_factor_uv, 4, 8,
                                which_mv, &xd->subpix);
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

      temp = xd->mode_info_context->bmi[yoffset + 0].as_mv[0].as_mv.row
             + xd->mode_info_context->bmi[yoffset + 1].as_mv[0].as_mv.row
             + xd->mode_info_context->bmi[yoffset + 4].as_mv[0].as_mv.row
             + xd->mode_info_context->bmi[yoffset + 5].as_mv[0].as_mv.row;

      if (temp < 0) temp -= 4;
      else temp += 4;

      blockd[uoffset].bmi.as_mv[0].as_mv.row = (temp / 8) &
                                                  xd->fullpixel_mask;

      temp = xd->mode_info_context->bmi[yoffset + 0].as_mv[0].as_mv.col
             + xd->mode_info_context->bmi[yoffset + 1].as_mv[0].as_mv.col
             + xd->mode_info_context->bmi[yoffset + 4].as_mv[0].as_mv.col
             + xd->mode_info_context->bmi[yoffset + 5].as_mv[0].as_mv.col;

      if (temp < 0) temp -= 4;
      else temp += 4;

      blockd[uoffset].bmi.as_mv[0].as_mv.col = (temp / 8) &
        xd->fullpixel_mask;

      // if (x->mode_info_context->mbmi.need_to_clamp_mvs)
      clamp_uvmv_to_umv_border(&blockd[uoffset].bmi.as_mv[0].as_mv, xd);

      // if (x->mode_info_context->mbmi.need_to_clamp_mvs)
      clamp_uvmv_to_umv_border(&blockd[uoffset].bmi.as_mv[0].as_mv, xd);

      blockd[voffset].bmi.as_mv[0].as_mv.row =
        blockd[uoffset].bmi.as_mv[0].as_mv.row;
      blockd[voffset].bmi.as_mv[0].as_mv.col =
        blockd[uoffset].bmi.as_mv[0].as_mv.col;

      if (xd->mode_info_context->mbmi.second_ref_frame > 0) {
        temp = xd->mode_info_context->bmi[yoffset + 0].as_mv[1].as_mv.row
               + xd->mode_info_context->bmi[yoffset + 1].as_mv[1].as_mv.row
               + xd->mode_info_context->bmi[yoffset + 4].as_mv[1].as_mv.row
               + xd->mode_info_context->bmi[yoffset + 5].as_mv[1].as_mv.row;

        if (temp < 0) {
          temp -= 4;
        } else {
          temp += 4;
        }

       blockd[uoffset].bmi.as_mv[1].as_mv.row = (temp / 8) &
                                                    xd->fullpixel_mask;

        temp = xd->mode_info_context->bmi[yoffset + 0].as_mv[1].as_mv.col
               + xd->mode_info_context->bmi[yoffset + 1].as_mv[1].as_mv.col
               + xd->mode_info_context->bmi[yoffset + 4].as_mv[1].as_mv.col
               + xd->mode_info_context->bmi[yoffset + 5].as_mv[1].as_mv.col;

        if (temp < 0) {
          temp -= 4;
        } else {
          temp += 4;
        }

        blockd[uoffset].bmi.as_mv[1].as_mv.col = (temp / 8) &
                                                        xd->fullpixel_mask;

        // if (mbmi->need_to_clamp_mvs)
        clamp_uvmv_to_umv_border(
          &blockd[uoffset].bmi.as_mv[1].as_mv, xd);

        // if (mbmi->need_to_clamp_mvs)
        clamp_uvmv_to_umv_border(
          &blockd[uoffset].bmi.as_mv[1].as_mv, xd);

        blockd[voffset].bmi.as_mv[1].as_mv.row =
          blockd[uoffset].bmi.as_mv[1].as_mv.row;
        blockd[voffset].bmi.as_mv[1].as_mv.col =
          blockd[uoffset].bmi.as_mv[1].as_mv.col;
      }
    }
  }
}

void vp9_build_inter16x16_predictors_mb(MACROBLOCKD *xd,
                                        uint8_t *dst_y,
                                        uint8_t *dst_u,
                                        uint8_t *dst_v,
                                        int dst_ystride,
                                        int dst_uvstride,
                                        int mb_row,
                                        int mb_col) {
  vp9_build_inter16x16_predictors_mby(xd, dst_y, dst_ystride, mb_row, mb_col);
  vp9_build_inter16x16_predictors_mbuv(xd, dst_u, dst_v, dst_uvstride,
                                       mb_row, mb_col);
}


void vp9_build_inter_predictors_mb(MACROBLOCKD *xd,
                                   int mb_row,
                                   int mb_col) {
  if (xd->mode_info_context->mbmi.mode != SPLITMV) {
    vp9_build_inter16x16_predictors_mb(xd, xd->predictor,
                                       &xd->predictor[256],
                                       &xd->predictor[320], 16, 8,
                                       mb_row, mb_col);

#if CONFIG_COMP_INTERINTRA_PRED
    if (xd->mode_info_context->mbmi.second_ref_frame == INTRA_FRAME) {
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
