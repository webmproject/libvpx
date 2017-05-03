/*
 *  Copyright (c) 2015 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include "./vpx_dsp_rtcd.h"
#include "vpx_dsp/x86/inv_txfm_sse2.h"
#include "vpx_dsp/x86/transpose_sse2.h"
#include "vpx_dsp/x86/txfm_common_sse2.h"

void vpx_highbd_idct32x32_1_add_sse2(const tran_low_t *input, uint16_t *dest,
                                     int stride, int bd) {
  __m128i dc_value, d;
  const __m128i zero = _mm_setzero_si128();
  const __m128i one = _mm_set1_epi16(1);
  const __m128i max = _mm_sub_epi16(_mm_slli_epi16(one, bd), one);
  int a, i, j;
  tran_low_t out;

  out = HIGHBD_WRAPLOW(dct_const_round_shift(input[0] * cospi_16_64), bd);
  out = HIGHBD_WRAPLOW(dct_const_round_shift(out * cospi_16_64), bd);
  a = ROUND_POWER_OF_TWO(out, 6);

  d = _mm_set1_epi32(a);
  dc_value = _mm_packs_epi32(d, d);
  for (i = 0; i < 32; ++i) {
    for (j = 0; j < 4; ++j) {
      d = _mm_loadu_si128((const __m128i *)(&dest[j * 8]));
      d = _mm_adds_epi16(d, dc_value);
      d = _mm_max_epi16(d, zero);
      d = _mm_min_epi16(d, max);
      _mm_storeu_si128((__m128i *)(&dest[j * 8]), d);
    }
    dest += stride;
  }
}
