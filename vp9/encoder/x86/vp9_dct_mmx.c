/*
 *  Copyright (c) 2014 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include <mmintrin.h>
#include <stdint.h>

#include "./vpx_config.h"

static void INLINE transpose_4x4_mmx(__m64* a, __m64* b, __m64* c, __m64* d) {
  __m64 w, x, y, z;
  w = _mm_unpacklo_pi16(*a, *b);
  x = _mm_unpackhi_pi16(*a, *b);
  y = _mm_unpacklo_pi16(*c, *d);
  z = _mm_unpackhi_pi16(*c, *d);
  *a = _mm_unpacklo_pi32(w, y);
  *b = _mm_unpackhi_pi32(w, y);
  *c = _mm_unpacklo_pi32(x, z);
  *d = _mm_unpackhi_pi32(x, z);
}

static void INLINE fwht_4x4_cols(__m64* out0,
                                 __m64* out1,
                                 __m64* out2,
                                 __m64* out3,
                                 __m64 a1,
                                 __m64 b1,
                                 __m64 c1,
                                 __m64 d1) {
  __m64 e1;

  a1 = _mm_add_pi16(a1, b1);
  d1 = _mm_sub_pi16(d1, c1);
  e1 = _mm_sub_pi16(a1, d1);
  e1 = _mm_srai_pi16(e1, 1);
  b1 = _mm_sub_pi16(e1, b1);
  c1 = _mm_sub_pi16(e1, c1);
  a1 = _mm_sub_pi16(a1, c1);
  d1 = _mm_add_pi16(d1, b1);
  *out0 = a1;
  *out1 = c1;
  *out2 = d1;
  *out3 = b1;
}

void vp9_fwht4x4_mmx(const int16_t* input, int16_t* output, int stride) {
  __m64 a1 = *(const __m64*)input;
  __m64 b1 = *(const __m64*)(input + stride);
  __m64 c1 = *(const __m64*)(input + 2 * stride);
  __m64 d1 = *(const __m64*)(input + 3 * stride);

  fwht_4x4_cols(&a1, &b1, &c1, &d1, a1, b1, c1, d1);
  transpose_4x4_mmx(&a1, &b1, &c1, &d1);
  fwht_4x4_cols(&a1, &b1, &c1, &d1, a1, b1, c1, d1);
  transpose_4x4_mmx(&a1, &b1, &c1, &d1);

  *(__m64*)output = _mm_slli_pi16(a1, 2);
  *(__m64*)(output + 4) = _mm_slli_pi16(b1, 2);
  *(__m64*)(output + 8) = _mm_slli_pi16(c1, 2);
  *(__m64*)(output + 12) = _mm_slli_pi16(d1, 2);
}
