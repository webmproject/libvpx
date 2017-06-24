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

static INLINE __m128i wraplow_16bit(const __m128i in0, const __m128i in1,
                                    const __m128i rounding) {
  __m128i temp[2];
  temp[0] = _mm_add_epi32(in0, rounding);
  temp[1] = _mm_add_epi32(in1, rounding);
  temp[0] = _mm_srai_epi32(temp[0], 4);
  temp[1] = _mm_srai_epi32(temp[1], 4);
  return _mm_packs_epi32(temp[0], temp[1]);
}

static INLINE __m128i dct_const_round_shift_64bit(const __m128i in) {
  const __m128i t = _mm_add_epi64(
      in,
      _mm_setr_epi32(DCT_CONST_ROUNDING << 2, 0, DCT_CONST_ROUNDING << 2, 0));
  return _mm_srli_si128(t, 2);
}

static INLINE __m128i pack_4(const __m128i in0, const __m128i in1) {
  const __m128i t0 = _mm_unpacklo_epi32(in0, in1);  // 0, 2
  const __m128i t1 = _mm_unpackhi_epi32(in0, in1);  // 1, 3
  return _mm_unpacklo_epi32(t0, t1);                // 0, 1, 2, 3
}

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

static INLINE void recon_and_store_4(uint16_t *const dest,
                                     const __m128i *const io, const int stride,
                                     int bd) {
  __m128i d0 = _mm_loadl_epi64((const __m128i *)dest);
  __m128i d2 = _mm_loadl_epi64((const __m128i *)(dest + stride * 2));
  d0 =
      _mm_unpacklo_epi64(d0, _mm_loadl_epi64((const __m128i *)(dest + stride)));
  d2 = _mm_unpacklo_epi64(
      d2, _mm_loadl_epi64((const __m128i *)(dest + stride * 3)));
  d0 = clamp_high_sse2(_mm_adds_epi16(d0, io[0]), bd);
  d2 = clamp_high_sse2(_mm_adds_epi16(d2, io[1]), bd);
  _mm_storel_epi64((__m128i *)dest, d0);
  d0 = _mm_srli_si128(d0, 8);
  _mm_storel_epi64((__m128i *)(dest + stride), d0);
  _mm_storel_epi64((__m128i *)(dest + stride * 2), d2);
  d2 = _mm_srli_si128(d2, 8);
  _mm_storel_epi64((__m128i *)(dest + stride * 3), d2);
}

#endif  // VPX_DSP_X86_HIGHBD_INV_TXFM_SSE2_H_
