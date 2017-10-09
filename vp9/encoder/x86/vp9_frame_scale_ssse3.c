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
#include "vpx_scale/yv12config.h"

static void scale_plane_2_to_1_phase_0(const uint8_t *src,
                                       const ptrdiff_t src_stride, uint8_t *dst,
                                       const ptrdiff_t dst_stride,
                                       const int dst_w, const int dst_h) {
  const __m128i mask = _mm_set1_epi16(0x00FF);
  const int max_width = (dst_w + 15) & ~15;
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

  if (dst_w * 2 == src_w && dst_h * 2 == src_h && phase_scaler == 0) {
    // 2 to 1
    const int dst_uv_w = dst_w / 2;
    const int dst_uv_h = dst_h / 2;
    scaled = 1;
    scale_plane_2_to_1_phase_0(src->y_buffer, src->y_stride, dst->y_buffer,
                               dst->y_stride, dst_w, dst_h);
    scale_plane_2_to_1_phase_0(src->u_buffer, src->uv_stride, dst->u_buffer,
                               dst->uv_stride, dst_uv_w, dst_uv_h);
    scale_plane_2_to_1_phase_0(src->v_buffer, src->uv_stride, dst->v_buffer,
                               dst->uv_stride, dst_uv_w, dst_uv_h);
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
