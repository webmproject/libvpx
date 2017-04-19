/*
 *  Copyright (c) 2017 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include <immintrin.h>

#include "./vpx_dsp_rtcd.h"
#include "vpx_dsp/x86/convolve.h"

// -----------------------------------------------------------------------------
// Copy and average

void vpx_highbd_convolve_copy_avx2(const uint8_t *src8, ptrdiff_t src_stride,
                                   uint8_t *dst8, ptrdiff_t dst_stride,
                                   const int16_t *filter_x, int filter_x_stride,
                                   const int16_t *filter_y, int filter_y_stride,
                                   int width, int h, int bd) {
  const uint16_t *src = CAST_TO_SHORTPTR(src8);
  uint16_t *dst = CAST_TO_SHORTPTR(dst8);
  (void)filter_x;
  (void)filter_y;
  (void)filter_x_stride;
  (void)filter_y_stride;
  (void)bd;

  assert(width % 4 == 0);
  if (width > 32) {  // width = 64
    do {
      const __m256i p0 = _mm256_loadu_si256((const __m256i *)src);
      const __m256i p1 = _mm256_loadu_si256((const __m256i *)(src + 16));
      const __m256i p2 = _mm256_loadu_si256((const __m256i *)(src + 32));
      const __m256i p3 = _mm256_loadu_si256((const __m256i *)(src + 48));
      src += src_stride;
      _mm256_storeu_si256((__m256i *)dst, p0);
      _mm256_storeu_si256((__m256i *)(dst + 16), p1);
      _mm256_storeu_si256((__m256i *)(dst + 32), p2);
      _mm256_storeu_si256((__m256i *)(dst + 48), p3);
      dst += dst_stride;
      h--;
    } while (h > 0);
  } else if (width > 16) {  // width = 32
    do {
      const __m256i p0 = _mm256_loadu_si256((const __m256i *)src);
      const __m256i p1 = _mm256_loadu_si256((const __m256i *)(src + 16));
      src += src_stride;
      _mm256_storeu_si256((__m256i *)dst, p0);
      _mm256_storeu_si256((__m256i *)(dst + 16), p1);
      dst += dst_stride;
      h--;
    } while (h > 0);
  } else if (width > 8) {  // width = 16
    __m256i p0, p1;
    do {
      p0 = _mm256_loadu_si256((const __m256i *)src);
      src += src_stride;
      p1 = _mm256_loadu_si256((const __m256i *)src);
      src += src_stride;

      _mm256_storeu_si256((__m256i *)dst, p0);
      dst += dst_stride;
      _mm256_storeu_si256((__m256i *)dst, p1);
      dst += dst_stride;
      h -= 2;
    } while (h > 0);
  } else if (width > 4) {  // width = 8
    __m128i p0, p1;
    do {
      p0 = _mm_loadu_si128((const __m128i *)src);
      src += src_stride;
      p1 = _mm_loadu_si128((const __m128i *)src);
      src += src_stride;

      _mm_storeu_si128((__m128i *)dst, p0);
      dst += dst_stride;
      _mm_storeu_si128((__m128i *)dst, p1);
      dst += dst_stride;
      h -= 2;
    } while (h > 0);
  } else {  // width = 4
    __m128i p0, p1;
    do {
      p0 = _mm_loadl_epi64((const __m128i *)src);
      src += src_stride;
      p1 = _mm_loadl_epi64((const __m128i *)src);
      src += src_stride;

      _mm_storel_epi64((__m128i *)dst, p0);
      dst += dst_stride;
      _mm_storel_epi64((__m128i *)dst, p1);
      dst += dst_stride;
      h -= 2;
    } while (h > 0);
  }
}

void vpx_highbd_convolve_avg_avx2(const uint8_t *src8, ptrdiff_t src_stride,
                                  uint8_t *dst8, ptrdiff_t dst_stride,
                                  const int16_t *filter_x, int filter_x_stride,
                                  const int16_t *filter_y, int filter_y_stride,
                                  int width, int h, int bd) {
  uint16_t *src = CAST_TO_SHORTPTR(src8);
  uint16_t *dst = CAST_TO_SHORTPTR(dst8);
  (void)filter_x;
  (void)filter_y;
  (void)filter_x_stride;
  (void)filter_y_stride;
  (void)bd;

  assert(width % 4 == 0);
  if (width > 32) {  // width = 64
    __m256i p0, p1, p2, p3, u0, u1, u2, u3;
    do {
      p0 = _mm256_loadu_si256((const __m256i *)src);
      p1 = _mm256_loadu_si256((const __m256i *)(src + 16));
      p2 = _mm256_loadu_si256((const __m256i *)(src + 32));
      p3 = _mm256_loadu_si256((const __m256i *)(src + 48));
      src += src_stride;
      u0 = _mm256_loadu_si256((const __m256i *)dst);
      u1 = _mm256_loadu_si256((const __m256i *)(dst + 16));
      u2 = _mm256_loadu_si256((const __m256i *)(dst + 32));
      u3 = _mm256_loadu_si256((const __m256i *)(dst + 48));
      _mm256_storeu_si256((__m256i *)dst, _mm256_avg_epu16(p0, u0));
      _mm256_storeu_si256((__m256i *)(dst + 16), _mm256_avg_epu16(p1, u1));
      _mm256_storeu_si256((__m256i *)(dst + 32), _mm256_avg_epu16(p2, u2));
      _mm256_storeu_si256((__m256i *)(dst + 48), _mm256_avg_epu16(p3, u3));
      dst += dst_stride;
      h--;
    } while (h > 0);
  } else if (width > 16) {  // width = 32
    __m256i p0, p1, u0, u1;
    do {
      p0 = _mm256_loadu_si256((const __m256i *)src);
      p1 = _mm256_loadu_si256((const __m256i *)(src + 16));
      src += src_stride;
      u0 = _mm256_loadu_si256((const __m256i *)dst);
      u1 = _mm256_loadu_si256((const __m256i *)(dst + 16));
      _mm256_storeu_si256((__m256i *)dst, _mm256_avg_epu16(p0, u0));
      _mm256_storeu_si256((__m256i *)(dst + 16), _mm256_avg_epu16(p1, u1));
      dst += dst_stride;
      h--;
    } while (h > 0);
  } else if (width > 8) {  // width = 16
    __m256i p0, p1, u0, u1;
    do {
      p0 = _mm256_loadu_si256((const __m256i *)src);
      p1 = _mm256_loadu_si256((const __m256i *)(src + src_stride));
      src += src_stride << 1;
      u0 = _mm256_loadu_si256((const __m256i *)dst);
      u1 = _mm256_loadu_si256((const __m256i *)(dst + dst_stride));

      _mm256_storeu_si256((__m256i *)dst, _mm256_avg_epu16(p0, u0));
      _mm256_storeu_si256((__m256i *)(dst + dst_stride),
                          _mm256_avg_epu16(p1, u1));
      dst += dst_stride << 1;
      h -= 2;
    } while (h > 0);
  } else if (width > 4) {  // width = 8
    __m128i p0, p1, u0, u1;
    do {
      p0 = _mm_loadu_si128((const __m128i *)src);
      p1 = _mm_loadu_si128((const __m128i *)(src + src_stride));
      src += src_stride << 1;
      u0 = _mm_loadu_si128((const __m128i *)dst);
      u1 = _mm_loadu_si128((const __m128i *)(dst + dst_stride));

      _mm_storeu_si128((__m128i *)dst, _mm_avg_epu16(p0, u0));
      _mm_storeu_si128((__m128i *)(dst + dst_stride), _mm_avg_epu16(p1, u1));
      dst += dst_stride << 1;
      h -= 2;
    } while (h > 0);
  } else {  // width = 4
    __m128i p0, p1, u0, u1;
    do {
      p0 = _mm_loadl_epi64((const __m128i *)src);
      p1 = _mm_loadl_epi64((const __m128i *)(src + src_stride));
      src += src_stride << 1;
      u0 = _mm_loadl_epi64((const __m128i *)dst);
      u1 = _mm_loadl_epi64((const __m128i *)(dst + dst_stride));

      _mm_storel_epi64((__m128i *)dst, _mm_avg_epu16(u0, p0));
      _mm_storel_epi64((__m128i *)(dst + dst_stride), _mm_avg_epu16(u1, p1));
      dst += dst_stride << 1;
      h -= 2;
    } while (h > 0);
  }
}
