/*
 *  Copyright (c) 2026 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include <immintrin.h>

#include "./vpx_config.h"
#include "./vpx_dsp_rtcd.h"
#include "vpx/vpx_integer.h"

void vpx_tm_predictor_8x8_avx2(uint8_t *dst, ptrdiff_t stride,
                               const uint8_t *above, const uint8_t *left) {
  __m128i top_left = _mm_set1_epi16(above[-1]);
  __m128i A = _mm_sub_epi16(_mm_cvtepu8_epi16(_mm_loadu_si64(above)), top_left);

  for (int i = 0; i < 4; i++) {
    __m128i L0 = _mm_set1_epi16(left[0]);
    __m128i D = _mm_add_epi16(A, L0);
    __m128i packed = _mm_packus_epi16(D, _mm_undefined_si128());
    _mm_storeu_si64(dst, packed);

    __m128i L1 = _mm_set1_epi16(left[1]);
    D = _mm_add_epi16(A, L1);
    packed = _mm_packus_epi16(D, _mm_undefined_si128());
    _mm_storeu_si64(dst + stride, packed);
    dst += 2 * stride;
    left += 2;
  }
}

void vpx_tm_predictor_16x16_avx2(uint8_t *dst, ptrdiff_t stride,
                                 const uint8_t *above, const uint8_t *left) {
  __m256i top_left = _mm256_set1_epi16(above[-1]);
  __m256i A =
      _mm256_sub_epi16(_mm256_cvtepu8_epi16(*(const __m128i *)above), top_left);

  for (int i = 0; i < 8; i++) {
    __m256i L0 = _mm256_set1_epi16(left[0]);
    __m256i D = _mm256_add_epi16(A, L0);
    __m128i packed = _mm_packus_epi16(_mm256_castsi256_si128(D),
                                      _mm256_extracti128_si256(D, 1));
    _mm_store_si128((__m128i *)dst, packed);

    __m256i L1 = _mm256_set1_epi16(left[1]);
    D = _mm256_add_epi16(A, L1);
    packed = _mm_packus_epi16(_mm256_castsi256_si128(D),
                              _mm256_extracti128_si256(D, 1));
    _mm_store_si128((__m128i *)(dst + stride), packed);
    dst += 2 * stride;
    left += 2;
  }
}

void vpx_tm_predictor_32x32_avx2(uint8_t *dst, ptrdiff_t stride,
                                 const uint8_t *above, const uint8_t *left) {
  __m256i top_left = _mm256_set1_epi16(above[-1]);
  __m256i zero = _mm256_setzero_si256();
  __m256i A = _mm256_load_si256((const __m256i *)above);
  // 0-7 || 16-23
  __m256i A0 = _mm256_sub_epi16(_mm256_unpacklo_epi8(A, zero), top_left);
  // 8-15 || 24-31
  __m256i A1 = _mm256_sub_epi16(_mm256_unpackhi_epi8(A, zero), top_left);

  for (int i = 0; i < 16; i++) {
    __m256i L0 = _mm256_set1_epi16(left[0]);
    __m256i D0 = _mm256_add_epi16(A0, L0);
    __m256i D1 = _mm256_add_epi16(A1, L0);

    // 256-bit packuswb packs within 128-bit lanes
    // bringing back 0-15 || 16-31
    __m256i packed = _mm256_packus_epi16(D0, D1);
    _mm256_store_si256((__m256i *)dst, packed);

    __m256i L1 = _mm256_set1_epi16(left[1]);
    D0 = _mm256_add_epi16(A0, L1);
    D1 = _mm256_add_epi16(A1, L1);

    packed = _mm256_packus_epi16(D0, D1);
    _mm256_store_si256((__m256i *)(dst + stride), packed);

    dst += 2 * stride;
    left += 2;
  }
}
