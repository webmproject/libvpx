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
#include "vpx_dsp/x86/highbd_inv_txfm_sse2.h"
#include "vpx_dsp/x86/inv_txfm_sse2.h"
#include "vpx_dsp/x86/transpose_sse2.h"

static INLINE __m128i dct_const_round_shift_4_sse2(const __m128i in0,
                                                   const __m128i in1) {
  const __m128i t0 = _mm_unpacklo_epi32(in0, in1);  // 0, 1
  const __m128i t1 = _mm_unpackhi_epi32(in0, in1);  // 2, 3
  const __m128i t2 = _mm_unpacklo_epi64(t0, t1);    // 0, 1, 2, 3
  return dct_const_round_shift_sse2(t2);
}

static INLINE void highbd_idct4_small_sse2(__m128i *const io) {
  const __m128i cospi_p16_p16 = _mm_setr_epi32(cospi_16_64, 0, cospi_16_64, 0);
  const __m128i cospi_p08_p08 = _mm_setr_epi32(cospi_8_64, 0, cospi_8_64, 0);
  const __m128i cospi_p24_p24 = _mm_setr_epi32(cospi_24_64, 0, cospi_24_64, 0);
  __m128i temp1[4], temp2[4], step[4];

  transpose_32bit_4x4(io, io);

  // Note: There is no 32-bit signed multiply SIMD instruction in SSE2.
  //       _mm_mul_epu32() is used which can only guarantee the lower 32-bit
  //       (signed) result is meaningful, which is enough in this function.

  // stage 1
  temp1[0] = _mm_add_epi32(io[0], io[2]);             // input[0] + input[2]
  temp2[0] = _mm_sub_epi32(io[0], io[2]);             // input[0] - input[2]
  temp1[1] = _mm_srli_si128(temp1[0], 4);             // 1, 3
  temp2[1] = _mm_srli_si128(temp2[0], 4);             // 1, 3
  temp1[0] = _mm_mul_epu32(temp1[0], cospi_p16_p16);  // ([0] + [2])*cospi_16_64
  temp1[1] = _mm_mul_epu32(temp1[1], cospi_p16_p16);  // ([0] + [2])*cospi_16_64
  temp2[0] = _mm_mul_epu32(temp2[0], cospi_p16_p16);  // ([0] - [2])*cospi_16_64
  temp2[1] = _mm_mul_epu32(temp2[1], cospi_p16_p16);  // ([0] - [2])*cospi_16_64
  step[0] = dct_const_round_shift_4_sse2(temp1[0], temp1[1]);
  step[1] = dct_const_round_shift_4_sse2(temp2[0], temp2[1]);

  temp1[3] = _mm_srli_si128(io[1], 4);
  temp2[3] = _mm_srli_si128(io[3], 4);
  temp1[0] = _mm_mul_epu32(io[1], cospi_p24_p24);     // input[1] * cospi_24_64
  temp1[1] = _mm_mul_epu32(temp1[3], cospi_p24_p24);  // input[1] * cospi_24_64
  temp2[0] = _mm_mul_epu32(io[1], cospi_p08_p08);     // input[1] * cospi_8_64
  temp2[1] = _mm_mul_epu32(temp1[3], cospi_p08_p08);  // input[1] * cospi_8_64
  temp1[2] = _mm_mul_epu32(io[3], cospi_p08_p08);     // input[3] * cospi_8_64
  temp1[3] = _mm_mul_epu32(temp2[3], cospi_p08_p08);  // input[3] * cospi_8_64
  temp2[2] = _mm_mul_epu32(io[3], cospi_p24_p24);     // input[3] * cospi_24_64
  temp2[3] = _mm_mul_epu32(temp2[3], cospi_p24_p24);  // input[3] * cospi_24_64
  temp1[0] = _mm_sub_epi64(temp1[0], temp1[2]);  // [1]*cospi_24 - [3]*cospi_8
  temp1[1] = _mm_sub_epi64(temp1[1], temp1[3]);  // [1]*cospi_24 - [3]*cospi_8
  temp2[0] = _mm_add_epi64(temp2[0], temp2[2]);  // [1]*cospi_8 + [3]*cospi_24
  temp2[1] = _mm_add_epi64(temp2[1], temp2[3]);  // [1]*cospi_8 + [3]*cospi_24
  step[2] = dct_const_round_shift_4_sse2(temp1[0], temp1[1]);
  step[3] = dct_const_round_shift_4_sse2(temp2[0], temp2[1]);

  // stage 2
  io[0] = _mm_add_epi32(step[0], step[3]);  // step[0] + step[3]
  io[1] = _mm_add_epi32(step[1], step[2]);  // step[1] + step[2]
  io[2] = _mm_sub_epi32(step[1], step[2]);  // step[1] - step[2]
  io[3] = _mm_sub_epi32(step[0], step[3]);  // step[0] - step[3]
}

static INLINE void abs_extend_64bit_sse2(const __m128i in,
                                         __m128i *const out /*out[2]*/,
                                         __m128i *const sign /*sign[2]*/) {
  sign[0] = _mm_srai_epi32(in, 31);
  out[0] = _mm_xor_si128(in, sign[0]);
  out[0] = _mm_sub_epi32(out[0], sign[0]);
  sign[1] = _mm_unpackhi_epi32(sign[0], sign[0]);  // 64-bit sign of 2, 3
  sign[0] = _mm_unpacklo_epi32(sign[0], sign[0]);  // 64-bit sign of 0, 1
  out[1] = _mm_unpackhi_epi32(out[0], out[0]);     // 2, 3
  out[0] = _mm_unpacklo_epi32(out[0], out[0]);     // 0, 1
}

static INLINE __m128i multiply_apply_sign_sse2(const __m128i in,
                                               const __m128i sign,
                                               const __m128i cospi) {
  __m128i out = _mm_mul_epu32(in, cospi);
  out = _mm_xor_si128(out, sign);
  return _mm_sub_epi64(out, sign);
}

static INLINE void highbd_idct4_large_sse2(__m128i *const io) {
  const __m128i cospi_p16_p16 =
      _mm_setr_epi32(cospi_16_64 << 2, 0, cospi_16_64 << 2, 0);
  const __m128i cospi_p08_p08 =
      _mm_setr_epi32(cospi_8_64 << 2, 0, cospi_8_64 << 2, 0);
  const __m128i cospi_p24_p24 =
      _mm_setr_epi32(cospi_24_64 << 2, 0, cospi_24_64 << 2, 0);
  __m128i temp1[4], temp2[4], step[4], sign1[4], sign2[4];

  transpose_32bit_4x4(io, io);

  // stage 1
  temp1[0] = _mm_add_epi32(io[0], io[2]);  // input[0] + input[2]
  temp2[0] = _mm_sub_epi32(io[0], io[2]);  // input[0] - input[2]
  abs_extend_64bit_sse2(temp1[0], temp1, sign1);
  abs_extend_64bit_sse2(temp2[0], temp2, sign2);
  temp1[0] = multiply_apply_sign_sse2(temp1[0], sign1[0], cospi_p16_p16);
  temp1[1] = multiply_apply_sign_sse2(temp1[1], sign1[1], cospi_p16_p16);
  temp2[0] = multiply_apply_sign_sse2(temp2[0], sign2[0], cospi_p16_p16);
  temp2[1] = multiply_apply_sign_sse2(temp2[1], sign2[1], cospi_p16_p16);
  temp1[0] = dct_const_round_shift_64bit(temp1[0]);
  temp1[1] = dct_const_round_shift_64bit(temp1[1]);
  temp2[0] = dct_const_round_shift_64bit(temp2[0]);
  temp2[1] = dct_const_round_shift_64bit(temp2[1]);
  step[0] = pack_4(temp1[0], temp1[1]);
  step[1] = pack_4(temp2[0], temp2[1]);

  abs_extend_64bit_sse2(io[1], temp1, sign1);
  abs_extend_64bit_sse2(io[3], temp2, sign2);
  temp1[2] = multiply_apply_sign_sse2(temp1[0], sign1[0], cospi_p08_p08);
  temp1[3] = multiply_apply_sign_sse2(temp1[1], sign1[1], cospi_p08_p08);
  temp1[0] = multiply_apply_sign_sse2(temp1[0], sign1[0], cospi_p24_p24);
  temp1[1] = multiply_apply_sign_sse2(temp1[1], sign1[1], cospi_p24_p24);
  temp2[2] = multiply_apply_sign_sse2(temp2[0], sign2[0], cospi_p24_p24);
  temp2[3] = multiply_apply_sign_sse2(temp2[1], sign2[1], cospi_p24_p24);
  temp2[0] = multiply_apply_sign_sse2(temp2[0], sign2[0], cospi_p08_p08);
  temp2[1] = multiply_apply_sign_sse2(temp2[1], sign2[1], cospi_p08_p08);
  temp1[0] = _mm_sub_epi64(temp1[0], temp2[0]);  // [1]*cospi_24 - [3]*cospi_8
  temp1[1] = _mm_sub_epi64(temp1[1], temp2[1]);  // [1]*cospi_24 - [3]*cospi_8
  temp2[0] = _mm_add_epi64(temp1[2], temp2[2]);  // [1]*cospi_8 + [3]*cospi_24
  temp2[1] = _mm_add_epi64(temp1[3], temp2[3]);  // [1]*cospi_8 + [3]*cospi_24
  temp1[0] = dct_const_round_shift_64bit(temp1[0]);
  temp1[1] = dct_const_round_shift_64bit(temp1[1]);
  temp2[0] = dct_const_round_shift_64bit(temp2[0]);
  temp2[1] = dct_const_round_shift_64bit(temp2[1]);
  step[2] = pack_4(temp1[0], temp1[1]);
  step[3] = pack_4(temp2[0], temp2[1]);

  // stage 2
  io[0] = _mm_add_epi32(step[0], step[3]);  // step[0] + step[3]
  io[1] = _mm_add_epi32(step[1], step[2]);  // step[1] + step[2]
  io[2] = _mm_sub_epi32(step[1], step[2]);  // step[1] - step[2]
  io[3] = _mm_sub_epi32(step[0], step[3]);  // step[0] - step[3]
}

void vpx_highbd_idct4x4_16_add_sse2(const tran_low_t *input, uint16_t *dest,
                                    int stride, int bd) {
  int16_t max = 0, min = 0;
  __m128i io[4], io_short[2];

  io[0] = _mm_load_si128((const __m128i *)(input + 0));
  io[1] = _mm_load_si128((const __m128i *)(input + 4));
  io[2] = _mm_load_si128((const __m128i *)(input + 8));
  io[3] = _mm_load_si128((const __m128i *)(input + 12));

  io_short[0] = _mm_packs_epi32(io[0], io[1]);
  io_short[1] = _mm_packs_epi32(io[2], io[3]);

  if (bd != 8) {
    __m128i max_input, min_input;

    max_input = _mm_max_epi16(io_short[0], io_short[1]);
    min_input = _mm_min_epi16(io_short[0], io_short[1]);
    max_input = _mm_max_epi16(max_input, _mm_srli_si128(max_input, 8));
    min_input = _mm_min_epi16(min_input, _mm_srli_si128(min_input, 8));
    max_input = _mm_max_epi16(max_input, _mm_srli_si128(max_input, 4));
    min_input = _mm_min_epi16(min_input, _mm_srli_si128(min_input, 4));
    max_input = _mm_max_epi16(max_input, _mm_srli_si128(max_input, 2));
    min_input = _mm_min_epi16(min_input, _mm_srli_si128(min_input, 2));
    max = _mm_extract_epi16(max_input, 0);
    min = _mm_extract_epi16(min_input, 0);
  }

  if (bd == 8 || (max < 4096 && min >= -4096)) {
    idct4_sse2(io_short);
    idct4_sse2(io_short);
    io_short[0] = _mm_add_epi16(io_short[0], _mm_set1_epi16(8));
    io_short[1] = _mm_add_epi16(io_short[1], _mm_set1_epi16(8));
    io[0] = _mm_srai_epi16(io_short[0], 4);
    io[1] = _mm_srai_epi16(io_short[1], 4);
  } else {
    if (max < 32767 && min > -32768) {
      highbd_idct4_small_sse2(io);
      highbd_idct4_small_sse2(io);
    } else {
      highbd_idct4_large_sse2(io);
      highbd_idct4_large_sse2(io);
    }
    io[0] = wraplow_16bit(io[0], io[1], _mm_set1_epi32(8));
    io[1] = wraplow_16bit(io[2], io[3], _mm_set1_epi32(8));
  }

  recon_and_store_4(io, dest, stride, bd);
}

void vpx_highbd_idct4x4_1_add_sse2(const tran_low_t *input, uint16_t *dest,
                                   int stride, int bd) {
  int a1, i;
  tran_low_t out;
  __m128i dc, d;

  out = HIGHBD_WRAPLOW(dct_const_round_shift(input[0] * cospi_16_64), bd);
  out = HIGHBD_WRAPLOW(dct_const_round_shift(out * cospi_16_64), bd);
  a1 = ROUND_POWER_OF_TWO(out, 4);
  dc = _mm_set1_epi16(a1);

  for (i = 0; i < 4; ++i) {
    d = _mm_loadl_epi64((const __m128i *)dest);
    d = add_clamp(d, dc, bd);
    _mm_storel_epi64((__m128i *)dest, d);
    dest += stride;
  }
}
