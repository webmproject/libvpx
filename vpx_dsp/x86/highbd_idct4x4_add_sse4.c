/*
 *  Copyright (c) 2017 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include <smmintrin.h>

#include "./vpx_dsp_rtcd.h"
#include "vpx_dsp/x86/highbd_inv_txfm_sse2.h"
#include "vpx_dsp/x86/inv_txfm_sse2.h"
#include "vpx_dsp/x86/transpose_sse2.h"

static INLINE void extend_64bit(const __m128i in,
                                __m128i *const out /*out[2]*/) {
  out[0] = _mm_unpacklo_epi32(in, in);  // 0, 0, 1, 1
  out[1] = _mm_unpackhi_epi32(in, in);  // 2, 2, 3, 3
}

static INLINE void highbd_idct4(__m128i *const io) {
  const __m128i cospi_p16_p16 =
      _mm_setr_epi32(cospi_16_64 << 2, 0, cospi_16_64 << 2, 0);
  const __m128i cospi_p08_p08 =
      _mm_setr_epi32(cospi_8_64 << 2, 0, cospi_8_64 << 2, 0);
  const __m128i cospi_p24_p24 =
      _mm_setr_epi32(cospi_24_64 << 2, 0, cospi_24_64 << 2, 0);
  __m128i temp1[4], temp2[4], step[4];

  transpose_32bit_4x4(&io[0], &io[1], &io[2], &io[3]);

  // stage 1
  temp1[0] = _mm_add_epi32(io[0], io[2]);  // input[0] + input[2]
  temp2[0] = _mm_sub_epi32(io[0], io[2]);  // input[0] - input[2]
  extend_64bit(temp1[0], temp1);
  extend_64bit(temp2[0], temp2);
  temp1[0] = _mm_mul_epi32(temp1[0], cospi_p16_p16);
  temp1[1] = _mm_mul_epi32(temp1[1], cospi_p16_p16);
  temp2[0] = _mm_mul_epi32(temp2[0], cospi_p16_p16);
  temp2[1] = _mm_mul_epi32(temp2[1], cospi_p16_p16);
  temp1[0] = dct_const_round_shift_64bit(temp1[0]);
  temp1[1] = dct_const_round_shift_64bit(temp1[1]);
  temp2[0] = dct_const_round_shift_64bit(temp2[0]);
  temp2[1] = dct_const_round_shift_64bit(temp2[1]);
  step[0] = pack_4(temp1[0], temp1[1]);
  step[1] = pack_4(temp2[0], temp2[1]);

  extend_64bit(io[1], temp1);
  extend_64bit(io[3], temp2);
  temp1[2] = _mm_mul_epi32(temp1[0], cospi_p08_p08);
  temp1[3] = _mm_mul_epi32(temp1[1], cospi_p08_p08);
  temp1[0] = _mm_mul_epi32(temp1[0], cospi_p24_p24);
  temp1[1] = _mm_mul_epi32(temp1[1], cospi_p24_p24);
  temp2[2] = _mm_mul_epi32(temp2[0], cospi_p24_p24);
  temp2[3] = _mm_mul_epi32(temp2[1], cospi_p24_p24);
  temp2[0] = _mm_mul_epi32(temp2[0], cospi_p08_p08);
  temp2[1] = _mm_mul_epi32(temp2[1], cospi_p08_p08);
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

void vpx_highbd_idct4x4_16_add_sse4_1(const tran_low_t *input, uint16_t *dest,
                                      int stride, int bd) {
  __m128i io[4];

  io[0] = _mm_load_si128((const __m128i *)(input + 0));
  io[1] = _mm_load_si128((const __m128i *)(input + 4));
  io[2] = _mm_load_si128((const __m128i *)(input + 8));
  io[3] = _mm_load_si128((const __m128i *)(input + 12));

  if (bd == 8) {
    __m128i io_short[2];

    io_short[0] = _mm_packs_epi32(io[0], io[1]);
    io_short[1] = _mm_packs_epi32(io[2], io[3]);
    idct4_sse2(io_short);
    idct4_sse2(io_short);
    io_short[0] = _mm_add_epi16(io_short[0], _mm_set1_epi16(8));
    io_short[1] = _mm_add_epi16(io_short[1], _mm_set1_epi16(8));
    io[0] = _mm_srai_epi16(io_short[0], 4);
    io[1] = _mm_srai_epi16(io_short[1], 4);
  } else {
    highbd_idct4(io);
    highbd_idct4(io);
    io[0] = wraplow_16bit(io[0], io[1], _mm_set1_epi32(8));
    io[1] = wraplow_16bit(io[2], io[3], _mm_set1_epi32(8));
  }

  recon_and_store_4(dest, io, stride, bd);
}
