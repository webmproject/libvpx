/*
 *  Copyright (c) 2017 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include <emmintrin.h>  // SSE2

#include "./vpx_config.h"
#include "./vpx_dsp_rtcd.h"
#include "vpx/vpx_integer.h"

// -----------------------------------------------------------------------------

void vpx_highbd_h_predictor_4x4_sse2(uint16_t *dst, ptrdiff_t stride,
                                     const uint16_t *above,
                                     const uint16_t *left, int bd) {
  const __m128i left_u16 = _mm_loadl_epi64((const __m128i *)left);
  const __m128i row0 = _mm_shufflelo_epi16(left_u16, 0x0);
  const __m128i row1 = _mm_shufflelo_epi16(left_u16, 0x55);
  const __m128i row2 = _mm_shufflelo_epi16(left_u16, 0xaa);
  const __m128i row3 = _mm_shufflelo_epi16(left_u16, 0xff);
  (void)above;
  (void)bd;
  _mm_storel_epi64((__m128i *)dst, row0);
  dst += stride;
  _mm_storel_epi64((__m128i *)dst, row1);
  dst += stride;
  _mm_storel_epi64((__m128i *)dst, row2);
  dst += stride;
  _mm_storel_epi64((__m128i *)dst, row3);
}

void vpx_highbd_h_predictor_8x8_sse2(uint16_t *dst, ptrdiff_t stride,
                                     const uint16_t *above,
                                     const uint16_t *left, int bd) {
  const __m128i left_u16 = _mm_load_si128((const __m128i *)left);
  const __m128i row0 = _mm_shufflelo_epi16(left_u16, 0x0);
  const __m128i row1 = _mm_shufflelo_epi16(left_u16, 0x55);
  const __m128i row2 = _mm_shufflelo_epi16(left_u16, 0xaa);
  const __m128i row3 = _mm_shufflelo_epi16(left_u16, 0xff);
  const __m128i row4 = _mm_shufflehi_epi16(left_u16, 0x0);
  const __m128i row5 = _mm_shufflehi_epi16(left_u16, 0x55);
  const __m128i row6 = _mm_shufflehi_epi16(left_u16, 0xaa);
  const __m128i row7 = _mm_shufflehi_epi16(left_u16, 0xff);
  (void)above;
  (void)bd;
  _mm_store_si128((__m128i *)dst, _mm_unpacklo_epi64(row0, row0));
  dst += stride;
  _mm_store_si128((__m128i *)dst, _mm_unpacklo_epi64(row1, row1));
  dst += stride;
  _mm_store_si128((__m128i *)dst, _mm_unpacklo_epi64(row2, row2));
  dst += stride;
  _mm_store_si128((__m128i *)dst, _mm_unpacklo_epi64(row3, row3));
  dst += stride;
  _mm_store_si128((__m128i *)dst, _mm_unpackhi_epi64(row4, row4));
  dst += stride;
  _mm_store_si128((__m128i *)dst, _mm_unpackhi_epi64(row5, row5));
  dst += stride;
  _mm_store_si128((__m128i *)dst, _mm_unpackhi_epi64(row6, row6));
  dst += stride;
  _mm_store_si128((__m128i *)dst, _mm_unpackhi_epi64(row7, row7));
}

static INLINE void h_store_16_unpacklo(uint16_t **dst, const ptrdiff_t stride,
                                       const __m128i *row) {
  const __m128i val = _mm_unpacklo_epi64(*row, *row);
  _mm_store_si128((__m128i *)*dst, val);
  _mm_store_si128((__m128i *)(*dst + 8), val);
  *dst += stride;
}

static INLINE void h_store_16_unpackhi(uint16_t **dst, const ptrdiff_t stride,
                                       const __m128i *row) {
  const __m128i val = _mm_unpackhi_epi64(*row, *row);
  _mm_store_si128((__m128i *)(*dst), val);
  _mm_store_si128((__m128i *)(*dst + 8), val);
  *dst += stride;
}

void vpx_highbd_h_predictor_16x16_sse2(uint16_t *dst, ptrdiff_t stride,
                                       const uint16_t *above,
                                       const uint16_t *left, int bd) {
  int i;
  (void)above;
  (void)bd;

  for (i = 0; i < 2; i++, left += 8) {
    const __m128i left_u16 = _mm_load_si128((const __m128i *)left);
    const __m128i row0 = _mm_shufflelo_epi16(left_u16, 0x0);
    const __m128i row1 = _mm_shufflelo_epi16(left_u16, 0x55);
    const __m128i row2 = _mm_shufflelo_epi16(left_u16, 0xaa);
    const __m128i row3 = _mm_shufflelo_epi16(left_u16, 0xff);
    const __m128i row4 = _mm_shufflehi_epi16(left_u16, 0x0);
    const __m128i row5 = _mm_shufflehi_epi16(left_u16, 0x55);
    const __m128i row6 = _mm_shufflehi_epi16(left_u16, 0xaa);
    const __m128i row7 = _mm_shufflehi_epi16(left_u16, 0xff);
    h_store_16_unpacklo(&dst, stride, &row0);
    h_store_16_unpacklo(&dst, stride, &row1);
    h_store_16_unpacklo(&dst, stride, &row2);
    h_store_16_unpacklo(&dst, stride, &row3);
    h_store_16_unpackhi(&dst, stride, &row4);
    h_store_16_unpackhi(&dst, stride, &row5);
    h_store_16_unpackhi(&dst, stride, &row6);
    h_store_16_unpackhi(&dst, stride, &row7);
  }
}

static INLINE void h_store_32_unpacklo(uint16_t **dst, const ptrdiff_t stride,
                                       const __m128i *row) {
  const __m128i val = _mm_unpacklo_epi64(*row, *row);
  _mm_store_si128((__m128i *)(*dst), val);
  _mm_store_si128((__m128i *)(*dst + 8), val);
  _mm_store_si128((__m128i *)(*dst + 16), val);
  _mm_store_si128((__m128i *)(*dst + 24), val);
  *dst += stride;
}

static INLINE void h_store_32_unpackhi(uint16_t **dst, const ptrdiff_t stride,
                                       const __m128i *row) {
  const __m128i val = _mm_unpackhi_epi64(*row, *row);
  _mm_store_si128((__m128i *)(*dst), val);
  _mm_store_si128((__m128i *)(*dst + 8), val);
  _mm_store_si128((__m128i *)(*dst + 16), val);
  _mm_store_si128((__m128i *)(*dst + 24), val);
  *dst += stride;
}

void vpx_highbd_h_predictor_32x32_sse2(uint16_t *dst, ptrdiff_t stride,
                                       const uint16_t *above,
                                       const uint16_t *left, int bd) {
  int i;
  (void)above;
  (void)bd;

  for (i = 0; i < 4; i++, left += 8) {
    const __m128i left_u16 = _mm_load_si128((const __m128i *)left);
    const __m128i row0 = _mm_shufflelo_epi16(left_u16, 0x0);
    const __m128i row1 = _mm_shufflelo_epi16(left_u16, 0x55);
    const __m128i row2 = _mm_shufflelo_epi16(left_u16, 0xaa);
    const __m128i row3 = _mm_shufflelo_epi16(left_u16, 0xff);
    const __m128i row4 = _mm_shufflehi_epi16(left_u16, 0x0);
    const __m128i row5 = _mm_shufflehi_epi16(left_u16, 0x55);
    const __m128i row6 = _mm_shufflehi_epi16(left_u16, 0xaa);
    const __m128i row7 = _mm_shufflehi_epi16(left_u16, 0xff);
    h_store_32_unpacklo(&dst, stride, &row0);
    h_store_32_unpacklo(&dst, stride, &row1);
    h_store_32_unpacklo(&dst, stride, &row2);
    h_store_32_unpacklo(&dst, stride, &row3);
    h_store_32_unpackhi(&dst, stride, &row4);
    h_store_32_unpackhi(&dst, stride, &row5);
    h_store_32_unpackhi(&dst, stride, &row6);
    h_store_32_unpackhi(&dst, stride, &row7);
  }
}
