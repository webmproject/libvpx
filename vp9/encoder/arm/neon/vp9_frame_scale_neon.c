/*
 *  Copyright (c) 2017 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include "./vp9_rtcd.h"
#include "./vpx_dsp_rtcd.h"
#include "./vpx_scale_rtcd.h"
#include "vp9/common/vp9_blockd.h"
#include "vpx_dsp/arm/transpose_neon.h"
#include "vpx_dsp/arm/vpx_convolve8_neon.h"
#include "vpx_dsp/vpx_filter.h"
#include "vpx_scale/yv12config.h"

// Note: The scaling functions could write extra rows and columns in dst, which
// exceed the right and bottom boundaries of the destination frame. We rely on
// the following frame extension function to fix these rows and columns.

static INLINE void scale_plane_2_to_1_phase_0_neon(const uint8_t *src,
                                                   const int src_stride,
                                                   uint8_t *dst,
                                                   const int dst_stride,
                                                   const int w, const int h) {
  const int max_width = (w + 15) & ~15;
  int y = h;

  assert(w && h);

  do {
    int x = max_width;
    do {
      const uint8x16x2_t s = vld2q_u8(src);
      vst1q_u8(dst, s.val[0]);
      src += 32;
      dst += 16;
      x -= 16;
    } while (x);
    src += 2 * (src_stride - max_width);
    dst += dst_stride - max_width;
  } while (--y);
}

static INLINE void scale_plane_2_to_1_bilinear_phase_non_0_neon(
    const uint8_t *const src, const int src_stride, uint8_t *dst,
    const int dst_stride, const int w, const int h, const int16_t c0,
    const int16_t c1) {
  const int max_width = (w + 7) & ~7;
  const uint8_t *src0 = src;
  const uint8_t *src1 = src + src_stride;
  const uint8x8_t coef0 = vdup_n_u8(c0);
  const uint8x8_t coef1 = vdup_n_u8(c1);
  int y = h;

  assert(w && h);

  do {
    int x = max_width;
    do {
      // 00 02 04 06 08 0A 0C 0E  01 03 05 07 09 0B 0D 0F
      const uint8x8x2_t s0 = vld2_u8(src0);
      // 10 12 14 16 18 1A 1C 1E  11 13 15 17 19 1B 1D 1F
      const uint8x8x2_t s1 = vld2_u8(src1);
      const uint16x8_t h0 = vmull_u8(s0.val[0], coef0);
      const uint16x8_t h1 = vmull_u8(s1.val[0], coef0);
      const uint16x8_t h2 = vmlal_u8(h0, s0.val[1], coef1);
      const uint16x8_t h3 = vmlal_u8(h1, s1.val[1], coef1);

      const uint8x8_t hor0 = vrshrn_n_u16(h2, 7);  // 00 01 02 03 04 05 06 07
      const uint8x8_t hor1 = vrshrn_n_u16(h3, 7);  // 10 11 12 13 14 15 16 17
      const uint16x8_t v0 = vmull_u8(hor0, coef0);
      const uint16x8_t v1 = vmlal_u8(v0, hor1, coef1);
      const uint8x8_t d = vrshrn_n_u16(v1, 7);  // 00 01 02 03 04 05 06 07
      vst1_u8(dst, d);
      src0 += 16;
      src1 += 16;
      dst += 8;
      x -= 8;
    } while (x);
    src0 += 2 * (src_stride - max_width);
    src1 += 2 * (src_stride - max_width);
    dst += dst_stride - max_width;
  } while (--y);
}

static INLINE uint8x8_t scale_filter(const uint8x8_t *const s,
                                     const int16x8_t filters) {
  const int16x8_t filter3 = vdupq_lane_s16(vget_low_s16(filters), 3);
  const int16x8_t filter4 = vdupq_lane_s16(vget_high_s16(filters), 0);
  int16x8_t ss[8];

  ss[0] = vreinterpretq_s16_u16(vmovl_u8(s[0]));
  ss[1] = vreinterpretq_s16_u16(vmovl_u8(s[1]));
  ss[2] = vreinterpretq_s16_u16(vmovl_u8(s[2]));
  ss[3] = vreinterpretq_s16_u16(vmovl_u8(s[3]));
  ss[4] = vreinterpretq_s16_u16(vmovl_u8(s[4]));
  ss[5] = vreinterpretq_s16_u16(vmovl_u8(s[5]));
  ss[6] = vreinterpretq_s16_u16(vmovl_u8(s[6]));
  ss[7] = vreinterpretq_s16_u16(vmovl_u8(s[7]));

  return convolve8_8(ss[0], ss[1], ss[2], ss[3], ss[4], ss[5], ss[6], ss[7],
                     filters, filter3, filter4);
}

static void scale_plane_2_to_1_general_neon(const uint8_t *src,
                                            const int src_stride, uint8_t *dst,
                                            const int dst_stride, const int w,
                                            const int h,
                                            const int16_t *const coef,
                                            uint8_t *const temp_buffer) {
  const int max_width_hor = (w + 3) & ~3;
  const int max_width_ver = (w + 7) & ~7;
  const int max_height_hor = (2 * h + SUBPEL_TAPS - 2 + 7) & ~7;
  const int max_height_ver = (h + 3) & ~3;
  const int16x8_t filters = vld1q_s16(coef);
  int x, y = max_height_hor;
  uint8_t *t = temp_buffer;
  uint8x8_t s[14], d[4];

  assert(w && h);

  src -= (SUBPEL_TAPS / 2 - 1) * src_stride + SUBPEL_TAPS / 2 + 1;

  // horizontal
  // Note: processing 8x8 is about 20% faster than processing row by row using
  // vld4_u8().
  do {
    load_u8_8x8(src + 2, src_stride, &s[0], &s[1], &s[2], &s[3], &s[4], &s[5],
                &s[6], &s[7]);
    transpose_u8_8x8(&s[0], &s[1], &s[2], &s[3], &s[4], &s[5], &s[6], &s[7]);
    x = max_width_hor;

    do {
      src += 8;
      load_u8_8x8(src, src_stride, &s[6], &s[7], &s[8], &s[9], &s[10], &s[11],
                  &s[12], &s[13]);
      transpose_u8_8x8(&s[6], &s[7], &s[8], &s[9], &s[10], &s[11], &s[12],
                       &s[13]);

      d[0] = scale_filter(&s[0], filters);
      d[1] = scale_filter(&s[2], filters);
      d[2] = scale_filter(&s[4], filters);
      d[3] = scale_filter(&s[6], filters);
      transpose_u8_8x4(&d[0], &d[1], &d[2], &d[3]);
      vst1_lane_u32((uint32_t *)(t + 0 * max_width_hor),
                    vreinterpret_u32_u8(d[0]), 0);
      vst1_lane_u32((uint32_t *)(t + 1 * max_width_hor),
                    vreinterpret_u32_u8(d[1]), 0);
      vst1_lane_u32((uint32_t *)(t + 2 * max_width_hor),
                    vreinterpret_u32_u8(d[2]), 0);
      vst1_lane_u32((uint32_t *)(t + 3 * max_width_hor),
                    vreinterpret_u32_u8(d[3]), 0);
      vst1_lane_u32((uint32_t *)(t + 4 * max_width_hor),
                    vreinterpret_u32_u8(d[0]), 1);
      vst1_lane_u32((uint32_t *)(t + 5 * max_width_hor),
                    vreinterpret_u32_u8(d[1]), 1);
      vst1_lane_u32((uint32_t *)(t + 6 * max_width_hor),
                    vreinterpret_u32_u8(d[2]), 1);
      vst1_lane_u32((uint32_t *)(t + 7 * max_width_hor),
                    vreinterpret_u32_u8(d[3]), 1);

      s[0] = s[8];
      s[1] = s[9];
      s[2] = s[10];
      s[3] = s[11];
      s[4] = s[12];
      s[5] = s[13];

      t += 4;
      x -= 4;
    } while (x);
    src += 8 * src_stride - 2 * max_width_hor;
    t += 7 * max_width_hor;
    y -= 8;
  } while (y);

  // vertical
  x = max_width_ver;
  t = temp_buffer;
  do {
    load_u8_8x8(t, max_width_hor, &s[0], &s[1], &s[2], &s[3], &s[4], &s[5],
                &s[6], &s[7]);
    t += 6 * max_width_hor;
    y = max_height_ver;

    do {
      load_u8_8x8(t, max_width_hor, &s[6], &s[7], &s[8], &s[9], &s[10], &s[11],
                  &s[12], &s[13]);
      t += 8 * max_width_hor;

      d[0] = scale_filter(&s[0], filters);
      d[1] = scale_filter(&s[2], filters);
      d[2] = scale_filter(&s[4], filters);
      d[3] = scale_filter(&s[6], filters);
      vst1_u8(dst + 0 * dst_stride, d[0]);
      vst1_u8(dst + 1 * dst_stride, d[1]);
      vst1_u8(dst + 2 * dst_stride, d[2]);
      vst1_u8(dst + 3 * dst_stride, d[3]);

      s[0] = s[8];
      s[1] = s[9];
      s[2] = s[10];
      s[3] = s[11];
      s[4] = s[12];
      s[5] = s[13];

      dst += 4 * dst_stride;
      y -= 4;
    } while (y);
    t -= max_width_hor * (2 * max_height_ver + 6);
    t += 8;
    dst -= max_height_ver * dst_stride;
    dst += 8;
    x -= 8;
  } while (x);
}

void vp9_scale_and_extend_frame_neon(const YV12_BUFFER_CONFIG *src,
                                     YV12_BUFFER_CONFIG *dst,
                                     INTERP_FILTER filter_type,
                                     int phase_scaler) {
  const int src_w = src->y_crop_width;
  const int src_h = src->y_crop_height;
  const int dst_w = dst->y_crop_width;
  const int dst_h = dst->y_crop_height;
  const int dst_uv_w = dst_w / 2;
  const int dst_uv_h = dst_h / 2;
  int scaled = 0;

  if (2 * dst_w == src_w && 2 * dst_h == src_h) {
    scaled = 1;
    if (phase_scaler == 0) {
      scale_plane_2_to_1_phase_0_neon(src->y_buffer, src->y_stride,
                                      dst->y_buffer, dst->y_stride, dst_w,
                                      dst_h);
      scale_plane_2_to_1_phase_0_neon(src->u_buffer, src->uv_stride,
                                      dst->u_buffer, dst->uv_stride, dst_uv_w,
                                      dst_uv_h);
      scale_plane_2_to_1_phase_0_neon(src->v_buffer, src->uv_stride,
                                      dst->v_buffer, dst->uv_stride, dst_uv_w,
                                      dst_uv_h);
    } else if (filter_type == BILINEAR) {
      const int16_t c0 = vp9_filter_kernels[BILINEAR][phase_scaler][3];
      const int16_t c1 = vp9_filter_kernels[BILINEAR][phase_scaler][4];
      scale_plane_2_to_1_bilinear_phase_non_0_neon(src->y_buffer, src->y_stride,
                                                   dst->y_buffer, dst->y_stride,
                                                   dst_w, dst_h, c0, c1);
      scale_plane_2_to_1_bilinear_phase_non_0_neon(
          src->u_buffer, src->uv_stride, dst->u_buffer, dst->uv_stride,
          dst_uv_w, dst_uv_h, c0, c1);
      scale_plane_2_to_1_bilinear_phase_non_0_neon(
          src->v_buffer, src->uv_stride, dst->v_buffer, dst->uv_stride,
          dst_uv_w, dst_uv_h, c0, c1);
    } else {
      const int buffer_max_width = (dst_w + 3) & ~3;
      const int buffer_max_height = (2 * dst_h + SUBPEL_TAPS - 2 + 7) & ~7;
      uint8_t *const temp_buffer =
          (uint8_t *)malloc(buffer_max_width * buffer_max_height);
      if (temp_buffer) {
        scale_plane_2_to_1_general_neon(
            src->y_buffer, src->y_stride, dst->y_buffer, dst->y_stride, dst_w,
            dst_h, vp9_filter_kernels[filter_type][phase_scaler], temp_buffer);
        scale_plane_2_to_1_general_neon(
            src->u_buffer, src->uv_stride, dst->u_buffer, dst->uv_stride,
            dst_uv_w, dst_uv_h, vp9_filter_kernels[filter_type][phase_scaler],
            temp_buffer);
        scale_plane_2_to_1_general_neon(
            src->v_buffer, src->uv_stride, dst->v_buffer, dst->uv_stride,
            dst_uv_w, dst_uv_h, vp9_filter_kernels[filter_type][phase_scaler],
            temp_buffer);
        free(temp_buffer);
      } else {
        scaled = 0;
      }
    }
  }

  if (scaled) {
    vpx_extend_frame_borders(dst);
  } else {
    vp9_scale_and_extend_frame_c(src, dst, filter_type, phase_scaler);
  }
}
