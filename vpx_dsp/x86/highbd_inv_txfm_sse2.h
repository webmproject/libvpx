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

static INLINE void extend_64bit(const __m128i in,
                                __m128i *const out /*out[2]*/) {
  out[0] = _mm_unpacklo_epi32(in, in);  // 0, 0, 1, 1
  out[1] = _mm_unpackhi_epi32(in, in);  // 2, 2, 3, 3
}

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

static INLINE __m128i add_clamp(const __m128i in0, const __m128i in1,
                                const int bd) {
  const __m128i zero = _mm_set1_epi16(0);
  // Faster than _mm_set1_epi16((1 << bd) - 1).
  const __m128i one = _mm_set1_epi16(1);
  const __m128i max = _mm_sub_epi16(_mm_slli_epi16(one, bd), one);
  __m128i d;

  d = _mm_adds_epi16(in0, in1);
  d = _mm_max_epi16(d, zero);
  d = _mm_min_epi16(d, max);

  return d;
}

static INLINE void highbd_idct_1_add_kernel(const tran_low_t *input,
                                            uint16_t *dest, int stride, int bd,
                                            const int size) {
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
      d = add_clamp(d, dc, bd);
      _mm_store_si128((__m128i *)(&dest[j * 8]), d);
    }
    dest += stride;
  }
}

static INLINE void recon_and_store_4_dual(const __m128i in,
                                          uint16_t *const dest,
                                          const int stride, const int bd) {
  __m128i d;

  d = _mm_loadl_epi64((const __m128i *)(dest + 0 * stride));
  d = _mm_castps_si128(
      _mm_loadh_pi(_mm_castsi128_ps(d), (const __m64 *)(dest + 1 * stride)));
  d = add_clamp(d, in, bd);
  _mm_storel_epi64((__m128i *)(dest + 0 * stride), d);
  _mm_storeh_pi((__m64 *)(dest + 1 * stride), _mm_castsi128_ps(d));
}

static INLINE void recon_and_store_4(const __m128i *const in, uint16_t *dest,
                                     const int stride, const int bd) {
  recon_and_store_4_dual(in[0], dest, stride, bd);
  dest += 2 * stride;
  recon_and_store_4_dual(in[1], dest, stride, bd);
}

#endif  // VPX_DSP_X86_HIGHBD_INV_TXFM_SSE2_H_
