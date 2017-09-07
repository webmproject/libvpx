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

static INLINE void scale_plane_4_to_1_phase_0_neon(const uint8_t *src,
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
      const uint8x16x4_t s = vld4q_u8(src);
      vst1q_u8(dst, s.val[0]);
      src += 64;
      dst += 16;
      x -= 16;
    } while (x);
    src += 4 * (src_stride - max_width);
    dst += dst_stride - max_width;
  } while (--y);
}

static INLINE void scale_plane_bilinear_phase_non_0_kernel(
    const uint8x16_t in0, const uint8x16_t in1, const uint8x16_t in2,
    const uint8x16_t in3, const uint8x8_t coef0, const uint8x8_t coef1,
    uint8_t *const dst) {
  const uint16x8_t h0 = vmull_u8(vget_low_u8(in0), coef0);
  const uint16x8_t h1 = vmull_u8(vget_high_u8(in0), coef0);
  const uint16x8_t h2 = vmull_u8(vget_low_u8(in2), coef0);
  const uint16x8_t h3 = vmull_u8(vget_high_u8(in2), coef0);
  const uint16x8_t h4 = vmlal_u8(h0, vget_low_u8(in1), coef1);
  const uint16x8_t h5 = vmlal_u8(h1, vget_high_u8(in1), coef1);
  const uint16x8_t h6 = vmlal_u8(h2, vget_low_u8(in3), coef1);
  const uint16x8_t h7 = vmlal_u8(h3, vget_high_u8(in3), coef1);

  const uint8x8_t hor0 = vrshrn_n_u16(h4, 7);  // temp: 00 01 02 03 04 05 06 07
  const uint8x8_t hor1 = vrshrn_n_u16(h5, 7);  // temp: 08 09 0A 0B 0C 0D 0E 0F
  const uint8x8_t hor2 = vrshrn_n_u16(h6, 7);  // temp: 10 11 12 13 14 15 16 17
  const uint8x8_t hor3 = vrshrn_n_u16(h7, 7);  // temp: 18 19 1A 1B 1C 1D 1E 1F
  const uint16x8_t v0 = vmull_u8(hor0, coef0);
  const uint16x8_t v1 = vmull_u8(hor1, coef0);
  const uint16x8_t v2 = vmlal_u8(v0, hor2, coef1);
  const uint16x8_t v3 = vmlal_u8(v1, hor3, coef1);
  // dst: 0 1 2 3 4 5 6 7  8 9 A B C D E F
  const uint8x16_t d = vcombine_u8(vrshrn_n_u16(v2, 7), vrshrn_n_u16(v3, 7));
  vst1q_u8(dst, d);
}

static INLINE void scale_plane_2_to_1_bilinear_phase_non_0_neon(
    const uint8_t *const src, const int src_stride, uint8_t *dst,
    const int dst_stride, const int w, const int h, const int16_t c0,
    const int16_t c1) {
  const int max_width = (w + 15) & ~15;
  const uint8_t *src0 = src;
  const uint8_t *src1 = src + src_stride;
  const uint8x8_t coef0 = vdup_n_u8(c0);
  const uint8x8_t coef1 = vdup_n_u8(c1);
  int y = h;

  assert(w && h);

  do {
    int x = max_width;
    do {
      // 000 002 004 006 008 00A 00C 00E  010 012 014 016 018 01A 01C 01E
      // 001 003 005 007 009 00B 00D 00F  011 013 015 017 019 01B 01D 01F
      const uint8x16x2_t s0 = vld2q_u8(src0);
      // 100 102 104 106 108 10A 10C 10E  110 112 114 116 118 11A 11C 11E
      // 101 103 105 107 109 10B 10D 10F  111 113 115 117 119 11B 11D 11F
      const uint8x16x2_t s1 = vld2q_u8(src1);
      scale_plane_bilinear_phase_non_0_kernel(s0.val[0], s0.val[1], s1.val[0],
                                              s1.val[1], coef0, coef1, dst);
      src0 += 32;
      src1 += 32;
      dst += 16;
      x -= 16;
    } while (x);
    src0 += 2 * (src_stride - max_width);
    src1 += 2 * (src_stride - max_width);
    dst += dst_stride - max_width;
  } while (--y);
}

static INLINE void scale_plane_4_to_1_bilinear_phase_non_0_neon(
    const uint8_t *const src, const int src_stride, uint8_t *dst,
    const int dst_stride, const int w, const int h, const int16_t c0,
    const int16_t c1) {
  const int max_width = (w + 15) & ~15;
  const uint8_t *src0 = src;
  const uint8_t *src1 = src + src_stride;
  const uint8x8_t coef0 = vdup_n_u8(c0);
  const uint8x8_t coef1 = vdup_n_u8(c1);
  int y = h;

  assert(w && h);

  do {
    int x = max_width;
    do {
      // (*) -- useless
      // 000 004 008 00C 010 014 018 01C  020 024 028 02C 030 034 038 03C
      // 001 005 009 00D 011 015 019 01D  021 025 029 02D 031 035 039 03D
      // 002 006 00A 00E 012 016 01A 01E  022 026 02A 02E 032 036 03A 03E (*)
      // 003 007 00B 00F 013 017 01B 01F  023 027 02B 02F 033 037 03B 03F (*)
      const uint8x16x4_t s0 = vld4q_u8(src0);
      // 100 104 108 10C 110 114 118 11C  120 124 128 12C 130 134 138 13C
      // 101 105 109 10D 111 115 119 11D  121 125 129 12D 131 135 139 13D
      // 102 106 10A 10E 112 116 11A 11E  122 126 12A 12E 132 136 13A 13E (*)
      // 103 107 10B 10F 113 117 11B 11F  123 127 12B 12F 133 137 13B 13F (*)
      const uint8x16x4_t s1 = vld4q_u8(src1);
      scale_plane_bilinear_phase_non_0_kernel(s0.val[0], s0.val[1], s1.val[0],
                                              s1.val[1], coef0, coef1, dst);
      src0 += 64;
      src1 += 64;
      dst += 16;
      x -= 16;
    } while (x);
    src0 += 4 * (src_stride - max_width);
    src1 += 4 * (src_stride - max_width);
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
  // Note: processing 4x8 is about 20% faster than processing row by row using
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

      d[0] = scale_filter(&s[0], filters);  // 00 10 20 30 40 50 60 70
      d[1] = scale_filter(&s[2], filters);  // 01 11 21 31 41 51 61 71
      d[2] = scale_filter(&s[4], filters);  // 02 12 22 32 42 52 62 72
      d[3] = scale_filter(&s[6], filters);  // 03 13 23 33 43 53 63 73
      // 00 01 02 03 40 41 42 43
      // 10 11 12 13 50 51 52 53
      // 20 21 22 23 60 61 62 63
      // 30 31 32 33 70 71 72 73
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

static void scale_plane_4_to_1_general_neon(const uint8_t *src,
                                            const int src_stride, uint8_t *dst,
                                            const int dst_stride, const int w,
                                            const int h,
                                            const int16_t *const coef,
                                            uint8_t *const temp_buffer) {
  const int max_width_hor = (w + 1) & ~1;
  const int max_width_ver = (w + 7) & ~7;
  const int max_height_hor = (4 * h + SUBPEL_TAPS - 2 + 7) & ~7;
  const int max_height_ver = (h + 1) & ~1;
  const int16x8_t filters = vld1q_s16(coef);
  int x, y = max_height_hor;
  uint8_t *t = temp_buffer;
  uint8x8_t s[12], d[2];

  assert(w && h);

  src -= (SUBPEL_TAPS / 2 - 1) * src_stride + SUBPEL_TAPS / 2 + 3;

  // horizontal
  // Note: processing 2x8 is about 20% faster than processing row by row using
  // vld4_u8().
  do {
    load_u8_8x8(src + 4, src_stride, &s[0], &s[1], &s[2], &s[3], &s[4], &s[5],
                &s[6], &s[7]);
    transpose_u8_4x8(&s[0], &s[1], &s[2], &s[3], s[4], s[5], s[6], s[7]);
    x = max_width_hor;

    do {
      uint8x8x2_t dd;
      src += 8;
      load_u8_8x8(src, src_stride, &s[4], &s[5], &s[6], &s[7], &s[8], &s[9],
                  &s[10], &s[11]);
      transpose_u8_8x8(&s[4], &s[5], &s[6], &s[7], &s[8], &s[9], &s[10],
                       &s[11]);

      d[0] = scale_filter(&s[0], filters);  // 00 10 20 30 40 50 60 70
      d[1] = scale_filter(&s[4], filters);  // 01 11 21 31 41 51 61 71
      // dd.val[0]: 00 01 20 21 40 41 60 61
      // dd.val[1]: 10 11 30 31 50 51 70 71
      dd = vtrn_u8(d[0], d[1]);
      vst1_lane_u16((uint16_t *)(t + 0 * max_width_hor),
                    vreinterpret_u16_u8(dd.val[0]), 0);
      vst1_lane_u16((uint16_t *)(t + 1 * max_width_hor),
                    vreinterpret_u16_u8(dd.val[1]), 0);
      vst1_lane_u16((uint16_t *)(t + 2 * max_width_hor),
                    vreinterpret_u16_u8(dd.val[0]), 1);
      vst1_lane_u16((uint16_t *)(t + 3 * max_width_hor),
                    vreinterpret_u16_u8(dd.val[1]), 1);
      vst1_lane_u16((uint16_t *)(t + 4 * max_width_hor),
                    vreinterpret_u16_u8(dd.val[0]), 2);
      vst1_lane_u16((uint16_t *)(t + 5 * max_width_hor),
                    vreinterpret_u16_u8(dd.val[1]), 2);
      vst1_lane_u16((uint16_t *)(t + 6 * max_width_hor),
                    vreinterpret_u16_u8(dd.val[0]), 3);
      vst1_lane_u16((uint16_t *)(t + 7 * max_width_hor),
                    vreinterpret_u16_u8(dd.val[1]), 3);

      s[0] = s[8];
      s[1] = s[9];
      s[2] = s[10];
      s[3] = s[11];

      t += 2;
      x -= 2;
    } while (x);
    src += 8 * src_stride - 4 * max_width_hor;
    t += 7 * max_width_hor;
    y -= 8;
  } while (y);

  // vertical
  x = max_width_ver;
  t = temp_buffer;
  do {
    load_u8_8x4(t, max_width_hor, &s[0], &s[1], &s[2], &s[3]);
    t += 4 * max_width_hor;
    y = max_height_ver;

    do {
      load_u8_8x8(t, max_width_hor, &s[4], &s[5], &s[6], &s[7], &s[8], &s[9],
                  &s[10], &s[11]);
      t += 8 * max_width_hor;

      d[0] = scale_filter(&s[0], filters);
      d[1] = scale_filter(&s[4], filters);
      vst1_u8(dst + 0 * dst_stride, d[0]);
      vst1_u8(dst + 1 * dst_stride, d[1]);

      s[0] = s[8];
      s[1] = s[9];
      s[2] = s[10];
      s[3] = s[11];

      dst += 2 * dst_stride;
      y -= 2;
    } while (y);
    t -= max_width_hor * (4 * max_height_ver + 4);
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
  } else if (4 * dst_w == src_w && 4 * dst_h == src_h) {
    scaled = 1;
    if (phase_scaler == 0) {
      scale_plane_4_to_1_phase_0_neon(src->y_buffer, src->y_stride,
                                      dst->y_buffer, dst->y_stride, dst_w,
                                      dst_h);
      scale_plane_4_to_1_phase_0_neon(src->u_buffer, src->uv_stride,
                                      dst->u_buffer, dst->uv_stride, dst_uv_w,
                                      dst_uv_h);
      scale_plane_4_to_1_phase_0_neon(src->v_buffer, src->uv_stride,
                                      dst->v_buffer, dst->uv_stride, dst_uv_w,
                                      dst_uv_h);
    } else if (filter_type == BILINEAR) {
      const int16_t c0 = vp9_filter_kernels[BILINEAR][phase_scaler][3];
      const int16_t c1 = vp9_filter_kernels[BILINEAR][phase_scaler][4];
      scale_plane_4_to_1_bilinear_phase_non_0_neon(src->y_buffer, src->y_stride,
                                                   dst->y_buffer, dst->y_stride,
                                                   dst_w, dst_h, c0, c1);
      scale_plane_4_to_1_bilinear_phase_non_0_neon(
          src->u_buffer, src->uv_stride, dst->u_buffer, dst->uv_stride,
          dst_uv_w, dst_uv_h, c0, c1);
      scale_plane_4_to_1_bilinear_phase_non_0_neon(
          src->v_buffer, src->uv_stride, dst->v_buffer, dst->uv_stride,
          dst_uv_w, dst_uv_h, c0, c1);
    } else {
      const int buffer_max_width = (dst_w + 1) & ~1;
      const int buffer_max_height = (4 * dst_h + SUBPEL_TAPS - 2 + 7) & ~7;
      uint8_t *const temp_buffer =
          (uint8_t *)malloc(buffer_max_width * buffer_max_height);
      if (temp_buffer) {
        scale_plane_4_to_1_general_neon(
            src->y_buffer, src->y_stride, dst->y_buffer, dst->y_stride, dst_w,
            dst_h, vp9_filter_kernels[filter_type][phase_scaler], temp_buffer);
        scale_plane_4_to_1_general_neon(
            src->u_buffer, src->uv_stride, dst->u_buffer, dst->uv_stride,
            dst_uv_w, dst_uv_h, vp9_filter_kernels[filter_type][phase_scaler],
            temp_buffer);
        scale_plane_4_to_1_general_neon(
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
