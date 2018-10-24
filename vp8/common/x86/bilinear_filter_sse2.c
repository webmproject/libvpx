/*
 *  Copyright (c) 2018 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include <assert.h>
#include <xmmintrin.h>

#include "./vp8_rtcd.h"
#include "vp8/common/filter.h"
#include "vpx_ports/mem.h"

static INLINE void horizontal_8xN(uint8_t *src, const int stride, uint16_t *dst,
                                  const int xoffset, const int height) {
  int h;
  const __m128i zero = _mm_setzero_si128();

  if (xoffset == 0) {
    for (h = 0; h < height; ++h) {
      const __m128i a = _mm_loadl_epi64((__m128i *)src);
      const __m128i a_u16 = _mm_unpacklo_epi8(a, zero);
      _mm_store_si128((__m128i *)dst, a_u16);
      src += stride;
      dst += 8;
    }
    return;
  }

  {
    const __m128i round_factor = _mm_set1_epi16(1 << (VP8_FILTER_SHIFT - 1));
    const __m128i hfilter_0 = _mm_set1_epi16(vp8_bilinear_filters[xoffset][0]);
    const __m128i hfilter_1 = _mm_set1_epi16(vp8_bilinear_filters[xoffset][1]);

    // Filter horizontally. Rather than load the whole array and transpose, load
    // 16 values (overreading) and shift to set up the second value. Do an
    // "extra" 9th line so the vertical pass has the necessary context.
    for (h = 0; h < height; ++h) {
      const __m128i a = _mm_loadu_si128((__m128i *)src);
      const __m128i b = _mm_srli_si128(a, 1);
      const __m128i a_u16 = _mm_unpacklo_epi8(a, zero);
      const __m128i b_u16 = _mm_unpacklo_epi8(b, zero);
      const __m128i a_filtered = _mm_mullo_epi16(a_u16, hfilter_0);
      const __m128i b_filtered = _mm_mullo_epi16(b_u16, hfilter_1);
      const __m128i sum = _mm_add_epi16(a_filtered, b_filtered);
      const __m128i compensated = _mm_add_epi16(sum, round_factor);
      const __m128i shifted = _mm_srai_epi16(compensated, VP8_FILTER_SHIFT);
      _mm_store_si128((__m128i *)dst, shifted);
      src += stride;
      dst += 8;
    }
  }
}

static INLINE void vertical_8xN(uint16_t *src, uint8_t *dst, const int stride,
                                const int yoffset, const int height) {
  int h;

  if (yoffset == 0) {
    for (h = 0; h < height; ++h) {
      const __m128i row = _mm_load_si128((__m128i *)src);
      const __m128i packed = _mm_packus_epi16(row, row);
      _mm_storel_epi64((__m128i *)dst, packed);
      src += 8;
      dst += stride;
    }
    return;
  }

  {
    const __m128i round_factor = _mm_set1_epi16(1 << (VP8_FILTER_SHIFT - 1));
    const __m128i vfilter_0 = _mm_set1_epi16(vp8_bilinear_filters[yoffset][0]);
    const __m128i vfilter_1 = _mm_set1_epi16(vp8_bilinear_filters[yoffset][1]);

    __m128i row_0 = _mm_load_si128((__m128i *)src);
    src += 8;
    for (h = 0; h < height; ++h) {
      const __m128i row_1 = _mm_load_si128((__m128i *)src);
      const __m128i row_0_filtered = _mm_mullo_epi16(row_0, vfilter_0);
      const __m128i row_1_filtered = _mm_mullo_epi16(row_1, vfilter_1);
      const __m128i sum = _mm_add_epi16(row_0_filtered, row_1_filtered);
      const __m128i compensated = _mm_add_epi16(sum, round_factor);
      const __m128i shifted = _mm_srai_epi16(compensated, VP8_FILTER_SHIFT);
      const __m128i packed = _mm_packus_epi16(shifted, shifted);
      _mm_storel_epi64((__m128i *)dst, packed);
      row_0 = row_1;
      src += 8;
      dst += stride;
    }
  }
}

void vp8_bilinear_predict8x8_sse2(uint8_t *src_ptr, int src_pixels_per_line,
                                  int xoffset, int yoffset, uint8_t *dst_ptr,
                                  int dst_pitch) {
  uint16_t FData[8 * 9];

  assert((xoffset | yoffset) != 0);

  horizontal_8xN(src_ptr, src_pixels_per_line, FData, xoffset, 9);

  vertical_8xN(FData, dst_ptr, dst_pitch, yoffset, 8);
}

void vp8_bilinear_predict8x4_sse2(uint8_t *src_ptr, int src_pixels_per_line,
                                  int xoffset, int yoffset, uint8_t *dst_ptr,
                                  int dst_pitch) {
  uint16_t FData[8 * 5];

  assert((xoffset | yoffset) != 0);

  horizontal_8xN(src_ptr, src_pixels_per_line, FData, xoffset, 5);

  vertical_8xN(FData, dst_ptr, dst_pitch, yoffset, 4);
}
