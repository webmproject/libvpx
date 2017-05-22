/*
 *  Copyright (c) 2015 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#ifndef VPX_DSP_X86_HIGHBD_INV_TXFM_SSE2_H_
#define VPX_DSP_X86_HIGHBD_INV_TXFM_SSE2_H_

#include <emmintrin.h>  // SSE2
#include "./vpx_config.h"
#include "vpx/vpx_integer.h"
#include "vpx_dsp/inv_txfm.h"
#include "vpx_dsp/x86/txfm_common_sse2.h"

static INLINE __m128i add_dc_clamp(const __m128i *const min,
                                   const __m128i *const max,
                                   const __m128i *const dc,
                                   const __m128i *const in) {
  __m128i out;
  out = _mm_adds_epi16(*in, *dc);
  out = _mm_max_epi16(out, *min);
  out = _mm_min_epi16(out, *max);
  return out;
}

static INLINE void highbd_idct_1_add_kernel(const tran_low_t *input,
                                            uint16_t *dest, int stride, int bd,
                                            const int size) {
  const __m128i zero = _mm_setzero_si128();
  // Faster than _mm_set1_epi16((1 << bd) - 1).
  const __m128i one = _mm_set1_epi16(1);
  const __m128i max = _mm_sub_epi16(_mm_slli_epi16(one, bd), one);
  int a1, i, j;
  tran_low_t out;
  __m128i dc, d;

  out = HIGHBD_WRAPLOW(dct_const_round_shift(input[0] * cospi_16_64), bd);
  out = HIGHBD_WRAPLOW(dct_const_round_shift(out * cospi_16_64), bd);
  a1 = ROUND_POWER_OF_TWO(out, (size == 8) ? 5 : 6);
  dc = _mm_set1_epi16(a1);

  for (i = 0; i < size; ++i) {
    for (j = 0; j < (size >> 3); ++j) {
      d = _mm_load_si128((const __m128i *)(&dest[j * 8]));
      d = add_dc_clamp(&zero, &max, &dc, &d);
      _mm_store_si128((__m128i *)(&dest[j * 8]), d);
    }
    dest += stride;
  }
}

static INLINE __m128i clamp_high_sse2(__m128i value, int bd) {
  __m128i ubounded, retval;
  const __m128i zero = _mm_set1_epi16(0);
  const __m128i one = _mm_set1_epi16(1);
  const __m128i max = _mm_sub_epi16(_mm_slli_epi16(one, bd), one);
  ubounded = _mm_cmpgt_epi16(value, max);
  retval = _mm_andnot_si128(ubounded, value);
  ubounded = _mm_and_si128(ubounded, max);
  retval = _mm_or_si128(retval, ubounded);
  retval = _mm_and_si128(retval, _mm_cmpgt_epi16(retval, zero));
  return retval;
}

#endif  // VPX_DSP_X86_HIGHBD_INV_TXFM_SSE2_H_
