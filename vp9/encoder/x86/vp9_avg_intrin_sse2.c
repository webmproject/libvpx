/*
 *  Copyright (c) 2014 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include <emmintrin.h>
#include "vpx_ports/mem.h"


unsigned int vp9_avg_8x8_sse2(const uint8_t *s, int p) {
  __m128i s0, s1, u0;
  unsigned int avg = 0;
  u0  = _mm_setzero_si128();
  s0 = _mm_unpacklo_epi8(_mm_loadl_epi64((const __m128i *)(s)), u0);
  s1 = _mm_unpacklo_epi8(_mm_loadl_epi64((const __m128i *)(s + p)), u0);
  s0 = _mm_adds_epu16(s0, s1);
  s1 = _mm_unpacklo_epi8(_mm_loadl_epi64((const __m128i *)(s + 2 * p)), u0);
  s0 = _mm_adds_epu16(s0, s1);
  s1 = _mm_unpacklo_epi8(_mm_loadl_epi64((const __m128i *)(s + 3 * p)), u0);
  s0 = _mm_adds_epu16(s0, s1);
  s1 = _mm_unpacklo_epi8(_mm_loadl_epi64((const __m128i *)(s + 4 * p)), u0);
  s0 = _mm_adds_epu16(s0, s1);
  s1 = _mm_unpacklo_epi8(_mm_loadl_epi64((const __m128i *)(s + 5 * p)), u0);
  s0 = _mm_adds_epu16(s0, s1);
  s1 = _mm_unpacklo_epi8(_mm_loadl_epi64((const __m128i *)(s + 6 * p)), u0);
  s0 = _mm_adds_epu16(s0, s1);
  s1 = _mm_unpacklo_epi8(_mm_loadl_epi64((const __m128i *)(s + 7 * p)), u0);
  s0 = _mm_adds_epu16(s0, s1);

  s0 = _mm_adds_epu16(s0, _mm_srli_si128(s0, 8));
  s0 = _mm_adds_epu16(s0, _mm_srli_epi64(s0, 32));
  s0 = _mm_adds_epu16(s0, _mm_srli_epi64(s0, 16));
  avg = _mm_extract_epi16(s0, 0);
  return (avg + 32) >> 6;
}

unsigned int vp9_avg_4x4_sse2(const uint8_t *s, int p) {
  __m128i s0, s1, u0;
  unsigned int avg = 0;
  u0  = _mm_setzero_si128();
  s0 = _mm_unpacklo_epi8(_mm_loadl_epi64((const __m128i *)(s)), u0);
  s1 = _mm_unpacklo_epi8(_mm_loadl_epi64((const __m128i *)(s + p)), u0);
  s0 = _mm_adds_epu16(s0, s1);
  s1 = _mm_unpacklo_epi8(_mm_loadl_epi64((const __m128i *)(s + 2 * p)), u0);
  s0 = _mm_adds_epu16(s0, s1);
  s1 = _mm_unpacklo_epi8(_mm_loadl_epi64((const __m128i *)(s + 3 * p)), u0);
  s0 = _mm_adds_epu16(s0, s1);

  s0 = _mm_adds_epu16(s0, _mm_srli_si128(s0, 4));
  s0 = _mm_adds_epu16(s0, _mm_srli_epi64(s0, 16));
  avg = _mm_extract_epi16(s0, 0);
  return (avg + 8) >> 4;
}

void vp9_int_pro_row_sse2(int16_t *hbuf, uint8_t const*ref,
                          const int ref_stride, const int height) {
  int idx;
  __m128i zero = _mm_setzero_si128();
  __m128i src_line = _mm_load_si128((const __m128i *)ref);
  __m128i s0 = _mm_unpacklo_epi8(src_line, zero);
  __m128i s1 = _mm_unpackhi_epi8(src_line, zero);
  __m128i t0, t1;
  int height_1 = height - 1;
  ref += ref_stride;

  for (idx = 1; idx < height_1; idx += 2) {
    src_line = _mm_load_si128((const __m128i *)ref);
    t0 = _mm_unpacklo_epi8(src_line, zero);
    t1 = _mm_unpackhi_epi8(src_line, zero);
    s0 = _mm_adds_epu16(s0, t0);
    s1 = _mm_adds_epu16(s1, t1);
    ref += ref_stride;

    src_line = _mm_load_si128((const __m128i *)ref);
    t0 = _mm_unpacklo_epi8(src_line, zero);
    t1 = _mm_unpackhi_epi8(src_line, zero);
    s0 = _mm_adds_epu16(s0, t0);
    s1 = _mm_adds_epu16(s1, t1);
    ref += ref_stride;
  }

  src_line = _mm_load_si128((const __m128i *)ref);
  t0 = _mm_unpacklo_epi8(src_line, zero);
  t1 = _mm_unpackhi_epi8(src_line, zero);
  s0 = _mm_adds_epu16(s0, t0);
  s1 = _mm_adds_epu16(s1, t1);

  _mm_store_si128((__m128i *)hbuf, s0);
  hbuf += 8;
  _mm_store_si128((__m128i *)hbuf, s1);
}

int16_t vp9_int_pro_col_sse2(uint8_t const *ref, const int width) {
  __m128i zero = _mm_setzero_si128();
  __m128i src_line = _mm_load_si128((const __m128i *)ref);
  __m128i s0 = _mm_sad_epu8(src_line, zero);
  __m128i s1;
  (void) width;  // width = 64

  ref += 16;
  src_line = _mm_load_si128((const __m128i *)ref);
  s1 = _mm_sad_epu8(src_line, zero);
  s0 = _mm_adds_epu16(s0, s1);

  ref += 16;
  src_line = _mm_load_si128((const __m128i *)ref);
  s1 = _mm_sad_epu8(src_line, zero);
  s0 = _mm_adds_epu16(s0, s1);

  ref += 16;
  src_line = _mm_load_si128((const __m128i *)ref);
  s1 = _mm_sad_epu8(src_line, zero);
  s0 = _mm_adds_epu16(s0, s1);

  s1 = _mm_srli_si128(s0, 8);
  s0 = _mm_adds_epu16(s0, s1);

  return _mm_extract_epi16(s0, 0);
}

int vp9_vector_sad_sse2(int16_t const *ref, int16_t const *src,
                        const int width) {
  int idx;
  __m128i zero = _mm_setzero_si128();
  __m128i sum;
  __m128i v0 = _mm_loadu_si128((const __m128i *)ref);
  __m128i v1 = _mm_load_si128((const __m128i *)src);
  __m128i diff = _mm_subs_epi16(v0, v1);
  __m128i sign = _mm_srai_epi16(diff, 15);

  diff = _mm_xor_si128(diff, sign);
  sum = _mm_sub_epi16(diff, sign);

  (void) width;  // width = 64;

  ref += 8;
  src += 8;

  v0 = _mm_unpacklo_epi16(sum, zero);
  v1 = _mm_unpackhi_epi16(sum, zero);
  sum = _mm_add_epi32(v0, v1);

  for (idx = 1; idx < 8; ++idx) {
    v0 = _mm_loadu_si128((const __m128i *)ref);
    v1 = _mm_load_si128((const __m128i *)src);
    diff = _mm_subs_epi16(v0, v1);
    sign = _mm_srai_epi16(diff, 15);
    diff = _mm_xor_si128(diff, sign);
    diff = _mm_sub_epi16(diff, sign);

    v0 = _mm_unpacklo_epi16(diff, zero);
    v1 = _mm_unpackhi_epi16(diff, zero);

    sum = _mm_add_epi32(sum, v0);
    sum = _mm_add_epi32(sum, v1);

    ref += 8;
    src += 8;
  }

  v0 = _mm_srli_si128(sum, 8);
  sum = _mm_add_epi32(sum, v0);
  v0 = _mm_srli_epi64(sum, 32);
  sum = _mm_add_epi32(sum, v0);

  return _mm_cvtsi128_si32(sum);
}
