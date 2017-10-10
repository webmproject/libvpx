/*
 *  Copyright (c) 2016 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include <tmmintrin.h>  // SSSE3

#include "./vp9_rtcd.h"
#include "./vpx_dsp_rtcd.h"
#include "./vpx_scale_rtcd.h"
#include "vpx_dsp/x86/convolve_ssse3.h"
#include "vpx_dsp/x86/mem_sse2.h"
#include "vpx_dsp/x86/transpose_sse2.h"
#include "vpx_scale/yv12config.h"

static void scale_plane_2_to_1_phase_0(const uint8_t *src,
                                       const ptrdiff_t src_stride, uint8_t *dst,
                                       const ptrdiff_t dst_stride,
                                       const int dst_w, const int dst_h) {
  const int max_width = (dst_w + 15) & ~15;
  const __m128i mask = _mm_set1_epi16(0x00FF);
  int y = dst_h;

  do {
    int x = max_width;
    do {
      const __m128i a = _mm_loadu_si128((const __m128i *)(src + 0));
      const __m128i b = _mm_loadu_si128((const __m128i *)(src + 16));
      const __m128i a_and = _mm_and_si128(a, mask);
      const __m128i b_and = _mm_and_si128(b, mask);
      const __m128i c = _mm_packus_epi16(a_and, b_and);
      _mm_storeu_si128((__m128i *)dst, c);
      src += 32;
      dst += 16;
      x -= 16;
    } while (x);
    src += 2 * (src_stride - max_width);
    dst += dst_stride - max_width;
  } while (--y);
}

static INLINE __m128i scale_plane_bilinear_kernel(const __m128i *const s,
                                                  const __m128i c0c1) {
  const __m128i k_64 = _mm_set1_epi16(1 << 6);
  const __m128i t0 = _mm_maddubs_epi16(s[0], c0c1);
  const __m128i t1 = _mm_maddubs_epi16(s[1], c0c1);
  // round and shift by 7 bit each 16 bit
  const __m128i t2 = _mm_adds_epi16(t0, k_64);
  const __m128i t3 = _mm_adds_epi16(t1, k_64);
  const __m128i t4 = _mm_srai_epi16(t2, 7);
  const __m128i t5 = _mm_srai_epi16(t3, 7);
  return _mm_packus_epi16(t4, t5);
}

static void scale_plane_2_to_1_bilinear(const uint8_t *src,
                                        const ptrdiff_t src_stride,
                                        uint8_t *dst,
                                        const ptrdiff_t dst_stride,
                                        const int dst_w, const int dst_h,
                                        const __m128i c0c1) {
  const int max_width = (dst_w + 15) & ~15;
  int y = dst_h;

  do {
    int x = max_width;
    do {
      __m128i s[2], d[2];

      // Horizontal
      // Even rows
      s[0] = _mm_loadu_si128((const __m128i *)(src + 0));
      s[1] = _mm_loadu_si128((const __m128i *)(src + 16));
      d[0] = scale_plane_bilinear_kernel(s, c0c1);

      // odd rows
      s[0] = _mm_loadu_si128((const __m128i *)(src + src_stride + 0));
      s[1] = _mm_loadu_si128((const __m128i *)(src + src_stride + 16));
      d[1] = scale_plane_bilinear_kernel(s, c0c1);

      // Vertical
      s[0] = _mm_unpacklo_epi8(d[0], d[1]);
      s[1] = _mm_unpackhi_epi8(d[0], d[1]);
      d[0] = scale_plane_bilinear_kernel(s, c0c1);

      _mm_storeu_si128((__m128i *)dst, d[0]);
      src += 32;
      dst += 16;
      x -= 16;
    } while (x);
    src += 2 * (src_stride - max_width);
    dst += dst_stride - max_width;
  } while (--y);
}

static void scale_plane_2_to_1_general(const uint8_t *src, const int src_stride,
                                       uint8_t *dst, const int dst_stride,
                                       const int w, const int h,
                                       const int16_t *const coef,
                                       uint8_t *const temp_buffer) {
  const int width_hor = (w + 3) & ~3;
  const int width_ver = (w + 7) & ~7;
  const int height_hor = (2 * h + SUBPEL_TAPS - 2 + 7) & ~7;
  const int height_ver = (h + 3) & ~3;
  int x, y = height_hor;
  uint8_t *t = temp_buffer;
  __m128i s[11], d[4];
  __m128i f[4];

  assert(w && h);

  shuffle_filter_ssse3(coef, f);
  src -= (SUBPEL_TAPS / 2 - 1) * src_stride + SUBPEL_TAPS / 2 + 1;

  // horizontal 4x8
  do {
    load_8bit_8x8(src + 2, src_stride, s);
    // 00 01 10 11 20 21 30 31  40 41 50 51 60 61 70 71
    // 02 03 12 13 22 23 32 33  42 43 52 53 62 63 72 73
    // 04 05 14 15 24 25 34 35  44 45 54 55 64 65 74 75
    // 06 07 16 17 26 27 36 37  46 47 56 57 66 67 76 77 (overlapped)
    transpose_16bit_4x8(s, s);
    x = width_hor;

    do {
      src += 8;
      load_8bit_8x8(src, src_stride, &s[3]);
      // 06 07 16 17 26 27 36 37  46 47 56 57 66 67 76 77
      // 08 09 18 19 28 29 38 39  48 49 58 59 68 69 78 79
      // 0A 0B 1A 1B 2A 2B 3A 3B  4A 4B 5A 5B 6A 6B 7A 7B
      // 0C 0D 1C 1D 2C 2D 3C 3D  4C 4D 5C 5D 6C 6D 7C 7D
      transpose_16bit_4x8(&s[3], &s[3]);

      d[0] = convolve8_8_ssse3(&s[0], f);  // 00 10 20 30 40 50 60 70
      d[1] = convolve8_8_ssse3(&s[1], f);  // 01 11 21 31 41 51 61 71
      d[2] = convolve8_8_ssse3(&s[2], f);  // 02 12 22 32 42 52 62 72
      d[3] = convolve8_8_ssse3(&s[3], f);  // 03 13 23 33 43 53 63 73

      // 00 10 20 30 40 50 60 70  02 12 22 32 42 52 62 72
      // 01 11 21 31 41 51 61 71  03 13 23 33 43 53 63 73
      d[0] = _mm_packus_epi16(d[0], d[2]);
      d[1] = _mm_packus_epi16(d[1], d[3]);
      // 00 10 01 11 20 30 21 31  40 50 41 51 60 70 61 71
      // 02 12 03 13 22 32 23 33  42 52 43 53 62 72 63 73
      d[2] = _mm_unpacklo_epi16(d[0], d[1]);
      d[3] = _mm_unpackhi_epi16(d[0], d[1]);
      store_8bit_4x4_sse2(d[2], t + 0, 2 * width_hor);
      store_8bit_4x4_sse2(d[3], t + 4, 2 * width_hor);

      s[0] = s[4];
      s[1] = s[5];
      s[2] = s[6];

      t += 8;
      x -= 4;
    } while (x);
    src += 8 * src_stride - 2 * width_hor;
    t += 6 * width_hor;
    y -= 8;
  } while (y);

  // vertical 8x4
  x = width_ver;
  t = temp_buffer;
  do {
    // 00 10 01 11 02 12 03 13  04 14 05 15 06 16 07 17
    // 20 30 21 31 22 32 23 33  24 34 25 35 26 36 27 37
    // 40 50 41 51 42 52 43 53  44 54 45 55 46 56 47 57
    // 60 70 61 71 62 72 63 73  64 74 65 75 66 76 67 77 (overlapped)
    loadu_8bit_16x4(t, 2 * width_hor, s);
    t += 6 * width_hor;
    y = height_ver;

    do {
      // 60 70 61 71 62 72 63 73  64 74 65 75 66 76 67 77
      // 80 90 81 91 82 92 83 93  84 94 85 95 86 96 87 77
      // A0 B0 A1 B1 A2 B2 A3 B3  A4 B4 A5 B5 A6 B6 A7 77
      // C0 D0 C1 D1 C2 D2 C3 D3  C4 D4 C5 D5 C6 D6 C7 77
      loadu_8bit_16x4(t, 2 * width_hor, &s[3]);
      t += 8 * width_hor;

      d[0] = convolve8_8_ssse3(&s[0], f);
      d[1] = convolve8_8_ssse3(&s[1], f);
      d[2] = convolve8_8_ssse3(&s[2], f);
      d[3] = convolve8_8_ssse3(&s[3], f);

      // 00 01 02 03 04 05 06 07  10 11 12 13 14 15 16 17
      // 20 21 22 23 24 25 26 27  30 31 32 33 34 35 36 37
      d[0] = _mm_packus_epi16(d[0], d[1]);
      d[1] = _mm_packus_epi16(d[2], d[3]);
      _mm_storel_epi64((__m128i *)(dst + 0 * dst_stride), d[0]);
      _mm_storeh_epi64((__m128i *)(dst + 1 * dst_stride), d[0]);
      _mm_storel_epi64((__m128i *)(dst + 2 * dst_stride), d[1]);
      _mm_storeh_epi64((__m128i *)(dst + 3 * dst_stride), d[1]);

      s[0] = s[4];
      s[1] = s[5];
      s[2] = s[6];

      dst += 4 * dst_stride;
      y -= 4;
    } while (y);
    t -= width_hor * (2 * height_ver + 6);
    t += 16;
    dst -= height_ver * dst_stride;
    dst += 8;
    x -= 8;
  } while (x);
}

static INLINE __m128i scale_1_to_2_phase_0_kernel(const __m128i *const s,
                                                  const __m128i *const f) {
  __m128i ss[4], temp;

  ss[0] = _mm_unpacklo_epi8(s[0], s[1]);
  ss[1] = _mm_unpacklo_epi8(s[2], s[3]);
  ss[2] = _mm_unpacklo_epi8(s[4], s[5]);
  ss[3] = _mm_unpacklo_epi8(s[6], s[7]);
  temp = convolve8_8_ssse3(ss, f);
  return _mm_packus_epi16(temp, temp);
}

// Only calculate odd columns since even columns are just src pixels' copies.
static void scale_1_to_2_phase_0_row(const uint8_t *src, uint8_t *dst,
                                     const int w, const __m128i *const f) {
  int x = w;

  do {
    __m128i s[8], temp;
    s[0] = _mm_loadl_epi64((const __m128i *)(src + 0));
    s[1] = _mm_loadl_epi64((const __m128i *)(src + 1));
    s[2] = _mm_loadl_epi64((const __m128i *)(src + 2));
    s[3] = _mm_loadl_epi64((const __m128i *)(src + 3));
    s[4] = _mm_loadl_epi64((const __m128i *)(src + 4));
    s[5] = _mm_loadl_epi64((const __m128i *)(src + 5));
    s[6] = _mm_loadl_epi64((const __m128i *)(src + 6));
    s[7] = _mm_loadl_epi64((const __m128i *)(src + 7));
    temp = scale_1_to_2_phase_0_kernel(s, f);
    _mm_storel_epi64((__m128i *)dst, temp);
    src += 8;
    dst += 8;
    x -= 8;
  } while (x);
}

static void scale_plane_1_to_2_phase_0(const uint8_t *src,
                                       const ptrdiff_t src_stride, uint8_t *dst,
                                       const ptrdiff_t dst_stride,
                                       const int src_w, const int src_h,
                                       const int16_t *const coef,
                                       uint8_t *const temp_buffer) {
  int max_width;
  int y;
  uint8_t *tmp[9];
  __m128i f[4];

  max_width = (src_w + 7) & ~7;
  tmp[0] = temp_buffer + 0 * max_width;
  tmp[1] = temp_buffer + 1 * max_width;
  tmp[2] = temp_buffer + 2 * max_width;
  tmp[3] = temp_buffer + 3 * max_width;
  tmp[4] = temp_buffer + 4 * max_width;
  tmp[5] = temp_buffer + 5 * max_width;
  tmp[6] = temp_buffer + 6 * max_width;
  tmp[7] = temp_buffer + 7 * max_width;

  shuffle_filter_ssse3(coef, f);

  scale_1_to_2_phase_0_row(src - 3 * src_stride - 3, tmp[0], max_width, f);
  scale_1_to_2_phase_0_row(src - 2 * src_stride - 3, tmp[1], max_width, f);
  scale_1_to_2_phase_0_row(src - 1 * src_stride - 3, tmp[2], max_width, f);
  scale_1_to_2_phase_0_row(src + 0 * src_stride - 3, tmp[3], max_width, f);
  scale_1_to_2_phase_0_row(src + 1 * src_stride - 3, tmp[4], max_width, f);
  scale_1_to_2_phase_0_row(src + 2 * src_stride - 3, tmp[5], max_width, f);
  scale_1_to_2_phase_0_row(src + 3 * src_stride - 3, tmp[6], max_width, f);

  y = src_h;
  do {
    int x;
    scale_1_to_2_phase_0_row(src + 4 * src_stride - 3, tmp[7], max_width, f);
    for (x = 0; x < max_width; x += 8) {
      __m128i s[8], C, D, CD;

      // Even rows
      const __m128i a = _mm_loadl_epi64((const __m128i *)(src + x));
      const __m128i b = _mm_loadl_epi64((const __m128i *)(tmp[3] + x));
      const __m128i ab = _mm_unpacklo_epi8(a, b);
      _mm_storeu_si128((__m128i *)(dst + 2 * x), ab);

      // Odd rows
      // Even columns
      load_8bit_8x8(src + x - 3 * src_stride, src_stride, s);
      C = scale_1_to_2_phase_0_kernel(s, f);

      // Odd columns
      s[0] = _mm_loadl_epi64((const __m128i *)(tmp[0] + x));
      s[1] = _mm_loadl_epi64((const __m128i *)(tmp[1] + x));
      s[2] = _mm_loadl_epi64((const __m128i *)(tmp[2] + x));
      s[3] = _mm_loadl_epi64((const __m128i *)(tmp[3] + x));
      s[4] = _mm_loadl_epi64((const __m128i *)(tmp[4] + x));
      s[5] = _mm_loadl_epi64((const __m128i *)(tmp[5] + x));
      s[6] = _mm_loadl_epi64((const __m128i *)(tmp[6] + x));
      s[7] = _mm_loadl_epi64((const __m128i *)(tmp[7] + x));
      D = scale_1_to_2_phase_0_kernel(s, f);

      CD = _mm_unpacklo_epi8(C, D);
      _mm_storeu_si128((__m128i *)(dst + dst_stride + 2 * x), CD);
    }

    src += src_stride;
    dst += 2 * dst_stride;
    tmp[8] = tmp[0];
    tmp[0] = tmp[1];
    tmp[1] = tmp[2];
    tmp[2] = tmp[3];
    tmp[3] = tmp[4];
    tmp[4] = tmp[5];
    tmp[5] = tmp[6];
    tmp[6] = tmp[7];
    tmp[7] = tmp[8];
  } while (--y);
}

void vp9_scale_and_extend_frame_ssse3(const YV12_BUFFER_CONFIG *src,
                                      YV12_BUFFER_CONFIG *dst,
                                      uint8_t filter_type, int phase_scaler) {
  const int src_w = src->y_crop_width;
  const int src_h = src->y_crop_height;
  const int dst_w = dst->y_crop_width;
  const int dst_h = dst->y_crop_height;
  int scaled = 0;

  // phase_scaler is usually 0 or 8.
  assert(phase_scaler >= 0 && phase_scaler < 16);

  if (dst_w * 2 == src_w && dst_h * 2 == src_h) {
    // 2 to 1
    const int dst_uv_w = dst_w / 2;
    const int dst_uv_h = dst_h / 2;
    scaled = 1;

    if (phase_scaler == 0) {
      scale_plane_2_to_1_phase_0(src->y_buffer, src->y_stride, dst->y_buffer,
                                 dst->y_stride, dst_w, dst_h);
      scale_plane_2_to_1_phase_0(src->u_buffer, src->uv_stride, dst->u_buffer,
                                 dst->uv_stride, dst_uv_w, dst_uv_h);
      scale_plane_2_to_1_phase_0(src->v_buffer, src->uv_stride, dst->v_buffer,
                                 dst->uv_stride, dst_uv_w, dst_uv_h);
    } else if (filter_type == BILINEAR) {
      const int16_t c0 = vp9_filter_kernels[BILINEAR][phase_scaler][3];
      const int16_t c1 = vp9_filter_kernels[BILINEAR][phase_scaler][4];
      const __m128i c0c1 = _mm_set1_epi16(c0 | (c1 << 8));  // c0 and c1 >= 0
      scale_plane_2_to_1_bilinear(src->y_buffer, src->y_stride, dst->y_buffer,
                                  dst->y_stride, dst_w, dst_h, c0c1);
      scale_plane_2_to_1_bilinear(src->u_buffer, src->uv_stride, dst->u_buffer,
                                  dst->uv_stride, dst_uv_w, dst_uv_h, c0c1);
      scale_plane_2_to_1_bilinear(src->v_buffer, src->uv_stride, dst->v_buffer,
                                  dst->uv_stride, dst_uv_w, dst_uv_h, c0c1);
    } else {
      const int buffer_stride = (dst_w + 3) & ~3;
      const int buffer_height = (2 * dst_h + SUBPEL_TAPS - 2 + 7) & ~7;
      uint8_t *const temp_buffer =
          (uint8_t *)malloc(buffer_stride * buffer_height);
      if (temp_buffer) {
        scale_plane_2_to_1_general(
            src->y_buffer, src->y_stride, dst->y_buffer, dst->y_stride, dst_w,
            dst_h, vp9_filter_kernels[filter_type][phase_scaler], temp_buffer);
        scale_plane_2_to_1_general(
            src->u_buffer, src->uv_stride, dst->u_buffer, dst->uv_stride,
            dst_uv_w, dst_uv_h, vp9_filter_kernels[filter_type][phase_scaler],
            temp_buffer);
        scale_plane_2_to_1_general(
            src->v_buffer, src->uv_stride, dst->v_buffer, dst->uv_stride,
            dst_uv_w, dst_uv_h, vp9_filter_kernels[filter_type][phase_scaler],
            temp_buffer);
        free(temp_buffer);
      } else {
        scaled = 0;
      }
    }
  } else if (dst_w == src_w * 2 && dst_h == src_h * 2 && phase_scaler == 0) {
    // 1 to 2
    uint8_t *const temp_buffer = (uint8_t *)malloc(8 * ((src_w + 7) & ~7));
    if (temp_buffer) {
      scaled = 1;
      scale_plane_1_to_2_phase_0(
          src->y_buffer, src->y_stride, dst->y_buffer, dst->y_stride, src_w,
          src_h, vp9_filter_kernels[filter_type][8], temp_buffer);
      scale_plane_1_to_2_phase_0(src->u_buffer, src->uv_stride, dst->u_buffer,
                                 dst->uv_stride, src_w / 2, src_h / 2,
                                 vp9_filter_kernels[filter_type][8],
                                 temp_buffer);
      scale_plane_1_to_2_phase_0(src->v_buffer, src->uv_stride, dst->v_buffer,
                                 dst->uv_stride, src_w / 2, src_h / 2,
                                 vp9_filter_kernels[filter_type][8],
                                 temp_buffer);
      free(temp_buffer);
    }
  }

  if (scaled) {
    vpx_extend_frame_borders(dst);
  } else {
    // Call c version for all other scaling ratios.
    vp9_scale_and_extend_frame_c(src, dst, filter_type, phase_scaler);
  }
}
