/*
 *  Copyright (c) 2022 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include <assert.h>
#include <immintrin.h>

#include "./vpx_dsp_rtcd.h"
#include "vpx/vpx_integer.h"

static VPX_FORCE_INLINE void subtract32_avx2(int16_t *diff_ptr,
                                             const uint8_t *src_ptr,
                                             const uint8_t *pred_ptr) {
  const __m256i s = _mm256_lddqu_si256((const __m256i *)src_ptr);
  const __m256i p = _mm256_lddqu_si256((const __m256i *)pred_ptr);
  const __m256i s_0 = _mm256_cvtepu8_epi16(_mm256_castsi256_si128(s));
  const __m256i s_1 = _mm256_cvtepu8_epi16(_mm256_extracti128_si256(s, 1));
  const __m256i p_0 = _mm256_cvtepu8_epi16(_mm256_castsi256_si128(p));
  const __m256i p_1 = _mm256_cvtepu8_epi16(_mm256_extracti128_si256(p, 1));
  const __m256i d_0 = _mm256_sub_epi16(s_0, p_0);
  const __m256i d_1 = _mm256_sub_epi16(s_1, p_1);
  _mm256_storeu_si256((__m256i *)diff_ptr, d_0);
  _mm256_storeu_si256((__m256i *)(diff_ptr + 16), d_1);
}

static VPX_FORCE_INLINE void subtract_block_16xn_avx2(
    int rows, int16_t *diff_ptr, ptrdiff_t diff_stride, const uint8_t *src_ptr,
    ptrdiff_t src_stride, const uint8_t *pred_ptr, ptrdiff_t pred_stride) {
  int j;
  for (j = 0; j < rows; ++j) {
    const __m128i s = _mm_lddqu_si128((const __m128i *)src_ptr);
    const __m128i p = _mm_lddqu_si128((const __m128i *)pred_ptr);
    const __m256i s_0 = _mm256_cvtepu8_epi16(s);
    const __m256i p_0 = _mm256_cvtepu8_epi16(p);
    const __m256i d_0 = _mm256_sub_epi16(s_0, p_0);
    _mm256_storeu_si256((__m256i *)diff_ptr, d_0);
    src_ptr += src_stride;
    pred_ptr += pred_stride;
    diff_ptr += diff_stride;
  }
}

static VPX_FORCE_INLINE void subtract_block_32xn_avx2(
    int rows, int16_t *diff_ptr, ptrdiff_t diff_stride, const uint8_t *src_ptr,
    ptrdiff_t src_stride, const uint8_t *pred_ptr, ptrdiff_t pred_stride) {
  int j;
  for (j = 0; j < rows; ++j) {
    subtract32_avx2(diff_ptr, src_ptr, pred_ptr);
    src_ptr += src_stride;
    pred_ptr += pred_stride;
    diff_ptr += diff_stride;
  }
}

static VPX_FORCE_INLINE void subtract_block_64xn_avx2(
    int rows, int16_t *diff_ptr, ptrdiff_t diff_stride, const uint8_t *src_ptr,
    ptrdiff_t src_stride, const uint8_t *pred_ptr, ptrdiff_t pred_stride) {
  int j;
  for (j = 0; j < rows; ++j) {
    subtract32_avx2(diff_ptr, src_ptr, pred_ptr);
    subtract32_avx2(diff_ptr + 32, src_ptr + 32, pred_ptr + 32);
    src_ptr += src_stride;
    pred_ptr += pred_stride;
    diff_ptr += diff_stride;
  }
}

void vpx_subtract_block_avx2(int rows, int cols, int16_t *diff_ptr,
                             ptrdiff_t diff_stride, const uint8_t *src_ptr,
                             ptrdiff_t src_stride, const uint8_t *pred_ptr,
                             ptrdiff_t pred_stride) {
  switch (cols) {
    case 16:
      subtract_block_16xn_avx2(rows, diff_ptr, diff_stride, src_ptr, src_stride,
                               pred_ptr, pred_stride);
      break;
    case 32:
      subtract_block_32xn_avx2(rows, diff_ptr, diff_stride, src_ptr, src_stride,
                               pred_ptr, pred_stride);
      break;
    case 64:
      subtract_block_64xn_avx2(rows, diff_ptr, diff_stride, src_ptr, src_stride,
                               pred_ptr, pred_stride);
      break;
    default:
      vpx_subtract_block_sse2(rows, cols, diff_ptr, diff_stride, src_ptr,
                              src_stride, pred_ptr, pred_stride);
      break;
  }
}
