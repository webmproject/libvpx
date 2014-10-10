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
