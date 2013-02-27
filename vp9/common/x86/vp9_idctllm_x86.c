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
// In order to improve performance, clip absolute diff values to [0, 255],
// which allows to keep the additions/subtractions in 8 bits.
void vp9_dc_only_idct_add_sse2(int input_dc, uint8_t *pred_ptr,
                               uint8_t *dst_ptr, int pitch, int stride) {
  int a1;
  int16_t out;
  uint8_t abs_diff;
  __m128i p0, p1, p2, p3;
  unsigned int extended_diff;
  __m128i diff;

  out = dct_const_round_shift(input_dc * cospi_16_64);
  out = dct_const_round_shift(out * cospi_16_64);
  a1 = ROUND_POWER_OF_TWO(out, 4);

  // Read prediction data.
  p0 = _mm_cvtsi32_si128 (*(const int *)(pred_ptr + 0 * pitch));
  p1 = _mm_cvtsi32_si128 (*(const int *)(pred_ptr + 1 * pitch));
  p2 = _mm_cvtsi32_si128 (*(const int *)(pred_ptr + 2 * pitch));
  p3 = _mm_cvtsi32_si128 (*(const int *)(pred_ptr + 3 * pitch));

  // Unpack prediction data, and store 4x4 array in 1 XMM register.
  p0 = _mm_unpacklo_epi32(p0, p1);
  p2 = _mm_unpacklo_epi32(p2, p3);
  p0 = _mm_unpacklo_epi64(p0, p2);

  // Clip dc value to [0, 255] range. Then, do addition or subtraction
  // according to its sign.
  if (a1 >= 0) {
    abs_diff = (a1 > 255) ? 255 : a1;
    extended_diff = abs_diff * 0x01010101u;
    diff = _mm_shuffle_epi32(_mm_cvtsi32_si128((int)extended_diff), 0);

    p1 = _mm_adds_epu8(p0, diff);
  } else {
    abs_diff = (a1 < -255) ? 255 : -a1;
    extended_diff = abs_diff * 0x01010101u;
    diff = _mm_shuffle_epi32(_mm_cvtsi32_si128((int)extended_diff), 0);

    p1 = _mm_subs_epu8(p0, diff);
  }

  // Store results to dst.
  *(int *)dst_ptr = _mm_cvtsi128_si32(p1);
  dst_ptr += stride;

  p1 = _mm_srli_si128(p1, 4);
  *(int *)dst_ptr = _mm_cvtsi128_si32(p1);
  dst_ptr += stride;

  p1 = _mm_srli_si128(p1, 4);
  *(int *)dst_ptr = _mm_cvtsi128_si32(p1);
  dst_ptr += stride;

  p1 = _mm_srli_si128(p1, 4);
  *(int *)dst_ptr = _mm_cvtsi128_si32(p1);
}
#endif
