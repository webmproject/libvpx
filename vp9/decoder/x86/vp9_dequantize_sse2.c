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

void vp9_add_constant_residual_8x8_sse2(const int16_t diff, uint8_t *dest,
                                        int stride) {
  uint8_t abs_diff;
  __m128i d;

  // Prediction data.
  __m128i p0 = _mm_loadl_epi64((const __m128i *)(dest + 0 * stride));
  __m128i p1 = _mm_loadl_epi64((const __m128i *)(dest + 1 * stride));
  __m128i p2 = _mm_loadl_epi64((const __m128i *)(dest + 2 * stride));
  __m128i p3 = _mm_loadl_epi64((const __m128i *)(dest + 3 * stride));
  __m128i p4 = _mm_loadl_epi64((const __m128i *)(dest + 4 * stride));
  __m128i p5 = _mm_loadl_epi64((const __m128i *)(dest + 5 * stride));
  __m128i p6 = _mm_loadl_epi64((const __m128i *)(dest + 6 * stride));
  __m128i p7 = _mm_loadl_epi64((const __m128i *)(dest + 7 * stride));

  p0 = _mm_unpacklo_epi64(p0, p1);
  p2 = _mm_unpacklo_epi64(p2, p3);
  p4 = _mm_unpacklo_epi64(p4, p5);
  p6 = _mm_unpacklo_epi64(p6, p7);

  // Clip diff value to [0, 255] range. Then, do addition or subtraction
  // according to its sign.
  if (diff >= 0) {
    abs_diff = (diff > 255) ? 255 : diff;
    d = _mm_shuffle_epi32(_mm_cvtsi32_si128((int)(abs_diff * 0x01010101u)), 0);

    p0 = _mm_adds_epu8(p0, d);
    p2 = _mm_adds_epu8(p2, d);
    p4 = _mm_adds_epu8(p4, d);
    p6 = _mm_adds_epu8(p6, d);
  } else {
    abs_diff = (diff < -255) ? 255 : -diff;
    d = _mm_shuffle_epi32(_mm_cvtsi32_si128((int)(abs_diff * 0x01010101u)), 0);

    p0 = _mm_subs_epu8(p0, d);
    p2 = _mm_subs_epu8(p2, d);
    p4 = _mm_subs_epu8(p4, d);
    p6 = _mm_subs_epu8(p6, d);
  }

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

void vp9_add_constant_residual_16x16_sse2(const int16_t diff, uint8_t *dest,
                                          int stride) {
  uint8_t abs_diff;
  __m128i d;

  // Prediction data.
  __m128i p0 = _mm_load_si128((const __m128i *)(dest + 0 * stride));
  __m128i p1 = _mm_load_si128((const __m128i *)(dest + 1 * stride));
  __m128i p2 = _mm_load_si128((const __m128i *)(dest + 2 * stride));
  __m128i p3 = _mm_load_si128((const __m128i *)(dest + 3 * stride));
  __m128i p4 = _mm_load_si128((const __m128i *)(dest + 4 * stride));
  __m128i p5 = _mm_load_si128((const __m128i *)(dest + 5 * stride));
  __m128i p6 = _mm_load_si128((const __m128i *)(dest + 6 * stride));
  __m128i p7 = _mm_load_si128((const __m128i *)(dest + 7 * stride));
  __m128i p8 = _mm_load_si128((const __m128i *)(dest + 8 * stride));
  __m128i p9 = _mm_load_si128((const __m128i *)(dest + 9 * stride));
  __m128i p10 = _mm_load_si128((const __m128i *)(dest + 10 * stride));
  __m128i p11 = _mm_load_si128((const __m128i *)(dest + 11 * stride));
  __m128i p12 = _mm_load_si128((const __m128i *)(dest + 12 * stride));
  __m128i p13 = _mm_load_si128((const __m128i *)(dest + 13 * stride));
  __m128i p14 = _mm_load_si128((const __m128i *)(dest + 14 * stride));
  __m128i p15 = _mm_load_si128((const __m128i *)(dest + 15 * stride));

  // Clip diff value to [0, 255] range. Then, do addition or subtraction
  // according to its sign.
  if (diff >= 0) {
    abs_diff = (diff > 255) ? 255 : diff;
    d = _mm_shuffle_epi32(_mm_cvtsi32_si128((int)(abs_diff * 0x01010101u)), 0);

    p0 = _mm_adds_epu8(p0, d);
    p1 = _mm_adds_epu8(p1, d);
    p2 = _mm_adds_epu8(p2, d);
    p3 = _mm_adds_epu8(p3, d);
    p4 = _mm_adds_epu8(p4, d);
    p5 = _mm_adds_epu8(p5, d);
    p6 = _mm_adds_epu8(p6, d);
    p7 = _mm_adds_epu8(p7, d);
    p8 = _mm_adds_epu8(p8, d);
    p9 = _mm_adds_epu8(p9, d);
    p10 = _mm_adds_epu8(p10, d);
    p11 = _mm_adds_epu8(p11, d);
    p12 = _mm_adds_epu8(p12, d);
    p13 = _mm_adds_epu8(p13, d);
    p14 = _mm_adds_epu8(p14, d);
    p15 = _mm_adds_epu8(p15, d);
  } else {
    abs_diff = (diff < -255) ? 255 : -diff;
    d = _mm_shuffle_epi32(_mm_cvtsi32_si128((int)(abs_diff * 0x01010101u)), 0);

    p0 = _mm_subs_epu8(p0, d);
    p1 = _mm_subs_epu8(p1, d);
    p2 = _mm_subs_epu8(p2, d);
    p3 = _mm_subs_epu8(p3, d);
    p4 = _mm_subs_epu8(p4, d);
    p5 = _mm_subs_epu8(p5, d);
    p6 = _mm_subs_epu8(p6, d);
    p7 = _mm_subs_epu8(p7, d);
    p8 = _mm_subs_epu8(p8, d);
    p9 = _mm_subs_epu8(p9, d);
    p10 = _mm_subs_epu8(p10, d);
    p11 = _mm_subs_epu8(p11, d);
    p12 = _mm_subs_epu8(p12, d);
    p13 = _mm_subs_epu8(p13, d);
    p14 = _mm_subs_epu8(p14, d);
    p15 = _mm_subs_epu8(p15, d);
  }

  // Store results
  _mm_store_si128((__m128i *)(dest + 0 * stride), p0);
  _mm_store_si128((__m128i *)(dest + 1 * stride), p1);
  _mm_store_si128((__m128i *)(dest + 2 * stride), p2);
  _mm_store_si128((__m128i *)(dest + 3 * stride), p3);
  _mm_store_si128((__m128i *)(dest + 4 * stride), p4);
  _mm_store_si128((__m128i *)(dest + 5 * stride), p5);
  _mm_store_si128((__m128i *)(dest + 6 * stride), p6);
  _mm_store_si128((__m128i *)(dest + 7 * stride), p7);
  _mm_store_si128((__m128i *)(dest + 8 * stride), p8);
  _mm_store_si128((__m128i *)(dest + 9 * stride), p9);
  _mm_store_si128((__m128i *)(dest + 10 * stride), p10);
  _mm_store_si128((__m128i *)(dest + 11 * stride), p11);
  _mm_store_si128((__m128i *)(dest + 12 * stride), p12);
  _mm_store_si128((__m128i *)(dest + 13 * stride), p13);
  _mm_store_si128((__m128i *)(dest + 14 * stride), p14);
  _mm_store_si128((__m128i *)(dest + 15 * stride), p15);
}

void vp9_add_constant_residual_32x32_sse2(const int16_t diff, uint8_t *dest,
                                          int stride) {
  uint8_t abs_diff;
  __m128i d;
  int i = 8;

  if (diff >= 0) {
    abs_diff = (diff > 255) ? 255 : diff;
    d = _mm_shuffle_epi32(_mm_cvtsi32_si128((int)(abs_diff * 0x01010101u)), 0);
  } else {
    abs_diff = (diff < -255) ? 255 : -diff;
    d = _mm_shuffle_epi32(_mm_cvtsi32_si128((int)(abs_diff * 0x01010101u)), 0);
  }

  do {
    // Prediction data.
    __m128i p0 = _mm_load_si128((const __m128i *)(dest + 0 * stride));
    __m128i p1 = _mm_load_si128((const __m128i *)(dest + 0 * stride + 16));
    __m128i p2 = _mm_load_si128((const __m128i *)(dest + 1 * stride));
    __m128i p3 = _mm_load_si128((const __m128i *)(dest + 1 * stride + 16));
    __m128i p4 = _mm_load_si128((const __m128i *)(dest + 2 * stride));
    __m128i p5 = _mm_load_si128((const __m128i *)(dest + 2 * stride + 16));
    __m128i p6 = _mm_load_si128((const __m128i *)(dest + 3 * stride));
    __m128i p7 = _mm_load_si128((const __m128i *)(dest + 3 * stride + 16));

    // Clip diff value to [0, 255] range. Then, do addition or subtraction
    // according to its sign.
    if (diff >= 0) {
      p0 = _mm_adds_epu8(p0, d);
      p1 = _mm_adds_epu8(p1, d);
      p2 = _mm_adds_epu8(p2, d);
      p3 = _mm_adds_epu8(p3, d);
      p4 = _mm_adds_epu8(p4, d);
      p5 = _mm_adds_epu8(p5, d);
      p6 = _mm_adds_epu8(p6, d);
      p7 = _mm_adds_epu8(p7, d);
    } else {
      p0 = _mm_subs_epu8(p0, d);
      p1 = _mm_subs_epu8(p1, d);
      p2 = _mm_subs_epu8(p2, d);
      p3 = _mm_subs_epu8(p3, d);
      p4 = _mm_subs_epu8(p4, d);
      p5 = _mm_subs_epu8(p5, d);
      p6 = _mm_subs_epu8(p6, d);
      p7 = _mm_subs_epu8(p7, d);
    }

    // Store results
    _mm_store_si128((__m128i *)(dest + 0 * stride), p0);
    _mm_store_si128((__m128i *)(dest + 0 * stride + 16), p1);
    _mm_store_si128((__m128i *)(dest + 1 * stride), p2);
    _mm_store_si128((__m128i *)(dest + 1 * stride + 16), p3);
    _mm_store_si128((__m128i *)(dest + 2 * stride), p4);
    _mm_store_si128((__m128i *)(dest + 2 * stride + 16), p5);
    _mm_store_si128((__m128i *)(dest + 3 * stride), p6);
    _mm_store_si128((__m128i *)(dest + 3 * stride + 16), p7);

    dest += 4 * stride;
  } while (--i);
}
