/*
 *  Copyright (c) 2012 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include <assert.h>
#include <emmintrin.h>  // SSE2
#include "./vpx_config.h"
#include "vpx/vpx_integer.h"
#include "vp9/common/vp9_common.h"
#include "vp9/common/vp9_idct.h"

#if HAVE_SSE2

void vp9_add_residual_4x4_sse2(const int16_t *diff, const uint8_t *pred,
                               int pitch, uint8_t *dest, int stride) {
  const int width = 4;
  const __m128i zero = _mm_setzero_si128();

  // Diff data
  const __m128i d0 = _mm_loadl_epi64((const __m128i *)(diff + 0 * width));
  const __m128i d1 = _mm_loadl_epi64((const __m128i *)(diff + 1 * width));
  const __m128i d2 = _mm_loadl_epi64((const __m128i *)(diff + 2 * width));
  const __m128i d3 = _mm_loadl_epi64((const __m128i *)(diff + 3 * width));

  // Prediction data.
  __m128i p0 = _mm_cvtsi32_si128(*(const int *)(pred + 0 * pitch));
  __m128i p1 = _mm_cvtsi32_si128(*(const int *)(pred + 1 * pitch));
  __m128i p2 = _mm_cvtsi32_si128(*(const int *)(pred + 2 * pitch));
  __m128i p3 = _mm_cvtsi32_si128(*(const int *)(pred + 3 * pitch));

  p0 = _mm_unpacklo_epi8(p0, zero);
  p1 = _mm_unpacklo_epi8(p1, zero);
  p2 = _mm_unpacklo_epi8(p2, zero);
  p3 = _mm_unpacklo_epi8(p3, zero);

  p0 = _mm_add_epi16(p0, d0);
  p1 = _mm_add_epi16(p1, d1);
  p2 = _mm_add_epi16(p2, d2);
  p3 = _mm_add_epi16(p3, d3);

  p0 = _mm_packus_epi16(p0, p1);
  p2 = _mm_packus_epi16(p2, p3);

  *(int *)dest = _mm_cvtsi128_si32(p0);
  dest += stride;

  p0 = _mm_srli_si128(p0, 8);
  *(int *)dest = _mm_cvtsi128_si32(p0);
  dest += stride;

  *(int *)dest = _mm_cvtsi128_si32(p2);
  dest += stride;

  p2 = _mm_srli_si128(p2, 8);
  *(int *)dest = _mm_cvtsi128_si32(p2);
}

void vp9_add_residual_8x8_sse2(const int16_t *diff, const uint8_t *pred,
                               int pitch, uint8_t *dest, int stride) {
  const int width = 8;
  const __m128i zero = _mm_setzero_si128();

  // Diff data
  const __m128i d0 = _mm_loadu_si128((const __m128i *)(diff + 0 * width));
  const __m128i d1 = _mm_loadu_si128((const __m128i *)(diff + 1 * width));
  const __m128i d2 = _mm_loadu_si128((const __m128i *)(diff + 2 * width));
  const __m128i d3 = _mm_loadu_si128((const __m128i *)(diff + 3 * width));
  const __m128i d4 = _mm_loadu_si128((const __m128i *)(diff + 4 * width));
  const __m128i d5 = _mm_loadu_si128((const __m128i *)(diff + 5 * width));
  const __m128i d6 = _mm_loadu_si128((const __m128i *)(diff + 6 * width));
  const __m128i d7 = _mm_loadu_si128((const __m128i *)(diff + 7 * width));

  // Prediction data.
  __m128i p0 = _mm_loadl_epi64((const __m128i *)(pred + 0 * pitch));
  __m128i p1 = _mm_loadl_epi64((const __m128i *)(pred + 1 * pitch));
  __m128i p2 = _mm_loadl_epi64((const __m128i *)(pred + 2 * pitch));
  __m128i p3 = _mm_loadl_epi64((const __m128i *)(pred + 3 * pitch));
  __m128i p4 = _mm_loadl_epi64((const __m128i *)(pred + 4 * pitch));
  __m128i p5 = _mm_loadl_epi64((const __m128i *)(pred + 5 * pitch));
  __m128i p6 = _mm_loadl_epi64((const __m128i *)(pred + 6 * pitch));
  __m128i p7 = _mm_loadl_epi64((const __m128i *)(pred + 7 * pitch));

  p0 = _mm_unpacklo_epi8(p0, zero);
  p1 = _mm_unpacklo_epi8(p1, zero);
  p2 = _mm_unpacklo_epi8(p2, zero);
  p3 = _mm_unpacklo_epi8(p3, zero);
  p4 = _mm_unpacklo_epi8(p4, zero);
  p5 = _mm_unpacklo_epi8(p5, zero);
  p6 = _mm_unpacklo_epi8(p6, zero);
  p7 = _mm_unpacklo_epi8(p7, zero);

  p0 = _mm_add_epi16(p0, d0);
  p1 = _mm_add_epi16(p1, d1);
  p2 = _mm_add_epi16(p2, d2);
  p3 = _mm_add_epi16(p3, d3);
  p4 = _mm_add_epi16(p4, d4);
  p5 = _mm_add_epi16(p5, d5);
  p6 = _mm_add_epi16(p6, d6);
  p7 = _mm_add_epi16(p7, d7);

  p0 = _mm_packus_epi16(p0, p1);
  p2 = _mm_packus_epi16(p2, p3);
  p4 = _mm_packus_epi16(p4, p5);
  p6 = _mm_packus_epi16(p6, p7);

  _mm_storel_epi64((__m128i *)(dest + 0 * stride), p0);
  p0 = _mm_srli_si128(p0, 8);
  _mm_storel_epi64((__m128i *)(dest + 1 * stride), p0);

  _mm_storel_epi64((__m128i *)(dest + 2 * stride), p2);
  p2 = _mm_srli_si128(p2, 8);
  _mm_storel_epi64((__m128i *)(dest + 3 * stride), p2);

  _mm_storel_epi64((__m128i *)(dest + 4 * stride), p4);
  p4 = _mm_srli_si128(p4, 8);
  _mm_storel_epi64((__m128i *)(dest + 5 * stride), p4);

  _mm_storel_epi64((__m128i *)(dest + 6 * stride), p6);
  p6 = _mm_srli_si128(p6, 8);
  _mm_storel_epi64((__m128i *)(dest + 7 * stride), p6);
}

void vp9_add_residual_16x16_sse2(const int16_t *diff, const uint8_t *pred,
                             int pitch, uint8_t *dest, int stride) {
  const int width = 16;
  int i = 4;
  const __m128i zero = _mm_setzero_si128();

  // Diff data
  __m128i d0, d1, d2, d3, d4, d5, d6, d7;
  __m128i p0, p1, p2, p3, p4, p5, p6, p7;

  do {
    d0 = _mm_loadu_si128((const __m128i *)(diff + 0 * width));
    d1 = _mm_loadu_si128((const __m128i *)(diff + 0 * width + 8));
    d2 = _mm_loadu_si128((const __m128i *)(diff + 1 * width));
    d3 = _mm_loadu_si128((const __m128i *)(diff + 1 * width + 8));
    d4 = _mm_loadu_si128((const __m128i *)(diff + 2 * width));
    d5 = _mm_loadu_si128((const __m128i *)(diff + 2 * width + 8));
    d6 = _mm_loadu_si128((const __m128i *)(diff + 3 * width));
    d7 = _mm_loadu_si128((const __m128i *)(diff + 3 * width + 8));

    // Prediction data.
    p1 = _mm_load_si128((const __m128i *)(pred + 0 * pitch));
    p3 = _mm_load_si128((const __m128i *)(pred + 1 * pitch));
    p5 = _mm_load_si128((const __m128i *)(pred + 2 * pitch));
    p7 = _mm_load_si128((const __m128i *)(pred + 3 * pitch));

    p0 = _mm_unpacklo_epi8(p1, zero);
    p1 = _mm_unpackhi_epi8(p1, zero);
    p2 = _mm_unpacklo_epi8(p3, zero);
    p3 = _mm_unpackhi_epi8(p3, zero);
    p4 = _mm_unpacklo_epi8(p5, zero);
    p5 = _mm_unpackhi_epi8(p5, zero);
    p6 = _mm_unpacklo_epi8(p7, zero);
    p7 = _mm_unpackhi_epi8(p7, zero);

    p0 = _mm_add_epi16(p0, d0);
    p1 = _mm_add_epi16(p1, d1);
    p2 = _mm_add_epi16(p2, d2);
    p3 = _mm_add_epi16(p3, d3);
    p4 = _mm_add_epi16(p4, d4);
    p5 = _mm_add_epi16(p5, d5);
    p6 = _mm_add_epi16(p6, d6);
    p7 = _mm_add_epi16(p7, d7);

    p0 = _mm_packus_epi16(p0, p1);
    p1 = _mm_packus_epi16(p2, p3);
    p2 = _mm_packus_epi16(p4, p5);
    p3 = _mm_packus_epi16(p6, p7);

    _mm_store_si128((__m128i *)(dest + 0 * stride), p0);
    _mm_store_si128((__m128i *)(dest + 1 * stride), p1);
    _mm_store_si128((__m128i *)(dest + 2 * stride), p2);
    _mm_store_si128((__m128i *)(dest + 3 * stride), p3);

    diff += 4 * width;
    pred += 4 * pitch;
    dest += 4 * stride;
  } while (--i);
}

void vp9_add_residual_32x32_sse2(const int16_t *diff, const uint8_t *pred,
                             int pitch, uint8_t *dest, int stride) {
  const int width = 32;
  int i = 16;
  const __m128i zero = _mm_setzero_si128();

  // Diff data
  __m128i d0, d1, d2, d3, d4, d5, d6, d7;
  __m128i p0, p1, p2, p3, p4, p5, p6, p7;

  do {
    d0 = _mm_loadu_si128((const __m128i *)(diff + 0 * width));
    d1 = _mm_loadu_si128((const __m128i *)(diff + 0 * width + 8));
    d2 = _mm_loadu_si128((const __m128i *)(diff + 0 * width + 16));
    d3 = _mm_loadu_si128((const __m128i *)(diff + 0 * width + 24));
    d4 = _mm_loadu_si128((const __m128i *)(diff + 1 * width));
    d5 = _mm_loadu_si128((const __m128i *)(diff + 1 * width + 8));
    d6 = _mm_loadu_si128((const __m128i *)(diff + 1 * width + 16));
    d7 = _mm_loadu_si128((const __m128i *)(diff + 1 * width + 24));

    // Prediction data.
    p1 = _mm_load_si128((const __m128i *)(pred + 0 * pitch));
    p3 = _mm_load_si128((const __m128i *)(pred + 0 * pitch + 16));
    p5 = _mm_load_si128((const __m128i *)(pred + 1 * pitch));
    p7 = _mm_load_si128((const __m128i *)(pred + 1 * pitch + 16));

    p0 = _mm_unpacklo_epi8(p1, zero);
    p1 = _mm_unpackhi_epi8(p1, zero);
    p2 = _mm_unpacklo_epi8(p3, zero);
    p3 = _mm_unpackhi_epi8(p3, zero);
    p4 = _mm_unpacklo_epi8(p5, zero);
    p5 = _mm_unpackhi_epi8(p5, zero);
    p6 = _mm_unpacklo_epi8(p7, zero);
    p7 = _mm_unpackhi_epi8(p7, zero);

    p0 = _mm_add_epi16(p0, d0);
    p1 = _mm_add_epi16(p1, d1);
    p2 = _mm_add_epi16(p2, d2);
    p3 = _mm_add_epi16(p3, d3);
    p4 = _mm_add_epi16(p4, d4);
    p5 = _mm_add_epi16(p5, d5);
    p6 = _mm_add_epi16(p6, d6);
    p7 = _mm_add_epi16(p7, d7);

    p0 = _mm_packus_epi16(p0, p1);
    p1 = _mm_packus_epi16(p2, p3);
    p2 = _mm_packus_epi16(p4, p5);
    p3 = _mm_packus_epi16(p6, p7);

    _mm_store_si128((__m128i *)(dest + 0 * stride), p0);
    _mm_store_si128((__m128i *)(dest + 0 * stride + 16), p1);
    _mm_store_si128((__m128i *)(dest + 1 * stride), p2);
    _mm_store_si128((__m128i *)(dest + 1 * stride + 16), p3);

    diff += 2 * width;
    pred += 2 * pitch;
    dest += 2 * stride;
  } while (--i);
}
#endif
