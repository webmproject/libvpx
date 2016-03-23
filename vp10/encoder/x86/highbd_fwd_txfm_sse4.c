/*
 *  Copyright (c) 2016 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include <assert.h>
#include <smmintrin.h> /* SSE4.1 */

#include "./vp10_rtcd.h"
#include "./vpx_config.h"
#include "vpx_dsp/txfm_common.h"
#include "vpx_ports/mem.h"

static INLINE void load_buffer_4x4(const int16_t *input, __m128i *in,
                                   int stride, int flipud, int fliplr) {
  const __m128i k__nonzero_bias_a = _mm_setr_epi32(0, 1, 1, 1);
  const __m128i k__nonzero_bias_b = _mm_setr_epi32(1, 0, 0, 0);
  __m128i mask;

  if (!flipud) {
    in[0] = _mm_loadl_epi64((const __m128i *)(input + 0 * stride));
    in[1] = _mm_loadl_epi64((const __m128i *)(input + 1 * stride));
    in[2] = _mm_loadl_epi64((const __m128i *)(input + 2 * stride));
    in[3] = _mm_loadl_epi64((const __m128i *)(input + 3 * stride));
  } else {
    in[0] = _mm_loadl_epi64((const __m128i *)(input + 3 * stride));
    in[1] = _mm_loadl_epi64((const __m128i *)(input + 2 * stride));
    in[2] = _mm_loadl_epi64((const __m128i *)(input + 1 * stride));
    in[3] = _mm_loadl_epi64((const __m128i *)(input + 0 * stride));
  }

  if (fliplr) {
    in[0] = _mm_shufflelo_epi16(in[0], 0x1b);
    in[1] = _mm_shufflelo_epi16(in[1], 0x1b);
    in[2] = _mm_shufflelo_epi16(in[2], 0x1b);
    in[3] = _mm_shufflelo_epi16(in[3], 0x1b);
  }

  in[0] = _mm_cvtepi16_epi32(in[0]);
  in[1] = _mm_cvtepi16_epi32(in[1]);
  in[2] = _mm_cvtepi16_epi32(in[2]);
  in[3] = _mm_cvtepi16_epi32(in[3]);

  in[0] = _mm_slli_epi32(in[0], 4);
  in[1] = _mm_slli_epi32(in[1], 4);
  in[2] = _mm_slli_epi32(in[2], 4);
  in[3] = _mm_slli_epi32(in[3], 4);

  mask = _mm_cmpeq_epi32(in[0], k__nonzero_bias_a);
  in[0] = _mm_add_epi32(in[0], mask);
  in[0] = _mm_add_epi32(in[0], k__nonzero_bias_b);
}

static void fdct4x4_sse4_1(__m128i *in) {
  const __m128i k__cospi_p16_p16 = _mm_set1_epi64x(cospi_16_64);
  const __m128i k__cospi_m16_m16 = _mm_set1_epi64x(-cospi_16_64);
  const __m128i k__cospi_p08_p08 = _mm_set1_epi64x(cospi_8_64);
  const __m128i k__cospi_m08_m08 = _mm_set1_epi64x(-cospi_8_64);
  const __m128i k__cospi_p24_p24 = _mm_set1_epi64x(cospi_24_64);
  const __m128i k__DCT_CONST_ROUNDING = _mm_set1_epi64x(DCT_CONST_ROUNDING);

  __m128i s[8];
  __m128i u0, u1, u2, u3;
  __m128i v0, v1, v2, v3, v4, v5, v6, v7;

  s[0] = _mm_add_epi32(in[0], in[3]);
  s[1] = _mm_add_epi32(in[1], in[2]);
  s[2] = _mm_sub_epi32(in[1], in[2]);
  s[3] = _mm_sub_epi32(in[0], in[3]);

  v0 = _mm_cvtepi32_epi64(s[0]);  // s01 s00
  v1 = _mm_cvtepi32_epi64(_mm_unpackhi_epi64(s[0], s[0]));  // s03 s02
  v2 = _mm_cvtepi32_epi64(s[1]);  // s11 s10
  v3 = _mm_cvtepi32_epi64(_mm_unpackhi_epi64(s[1], s[1]));  // s13 s12
  v4 = _mm_cvtepi32_epi64(s[2]);  // s21 s20
  v5 = _mm_cvtepi32_epi64(_mm_unpackhi_epi64(s[2], s[2]));  // s23 s22
  v6 = _mm_cvtepi32_epi64(s[3]);  // s31 s30
  v7 = _mm_cvtepi32_epi64(_mm_unpackhi_epi64(s[3], s[3]));  // s33 s32

  u0 = _mm_mul_epi32(v0, k__cospi_p16_p16);
  u1 = _mm_mul_epi32(v1, k__cospi_p16_p16);
  u2 = _mm_mul_epi32(v2, k__cospi_p16_p16);
  u3 = _mm_mul_epi32(v3, k__cospi_p16_p16);

  s[0] = _mm_add_epi64(u0, u2);  // y10 y00
  s[1] = _mm_add_epi64(u1, u3);  // y30 y20

  u2 = _mm_mul_epi32(v2, k__cospi_m16_m16);
  u3 = _mm_mul_epi32(v3, k__cospi_m16_m16);

  s[2] = _mm_add_epi64(u0, u2);  // y12 y02
  s[3] = _mm_add_epi64(u1, u3);  // y32 y22

  u0 = _mm_mul_epi32(v6, k__cospi_p08_p08);
  u1 = _mm_mul_epi32(v5, k__cospi_p24_p24);
  u2 = _mm_mul_epi32(v4, k__cospi_p24_p24);
  u3 = _mm_mul_epi32(v7, k__cospi_p08_p08);

  s[4] = _mm_add_epi64(u0, u2);  // y11 y01
  s[5] = _mm_add_epi64(u1, u3);  // y31 y21

  u0 = _mm_mul_epi32(v4, k__cospi_m08_m08);
  u1 = _mm_mul_epi32(v5, k__cospi_m08_m08);
  u2 = _mm_mul_epi32(v6, k__cospi_p24_p24);
  u3 = _mm_mul_epi32(v7, k__cospi_p24_p24);

  s[6] = _mm_add_epi64(u0, u2);  // y13 y03
  s[7] = _mm_add_epi64(u1, u3);  // y33 y23

  s[0] = _mm_add_epi64(s[0], k__DCT_CONST_ROUNDING);
  s[1] = _mm_add_epi64(s[1], k__DCT_CONST_ROUNDING);
  s[2] = _mm_add_epi64(s[2], k__DCT_CONST_ROUNDING);
  s[3] = _mm_add_epi64(s[3], k__DCT_CONST_ROUNDING);
  s[4] = _mm_add_epi64(s[4], k__DCT_CONST_ROUNDING);
  s[5] = _mm_add_epi64(s[5], k__DCT_CONST_ROUNDING);
  s[6] = _mm_add_epi64(s[6], k__DCT_CONST_ROUNDING);
  s[7] = _mm_add_epi64(s[7], k__DCT_CONST_ROUNDING);

  s[0] = _mm_srli_epi64(s[0], DCT_CONST_BITS);
  s[1] = _mm_srli_epi64(s[1], DCT_CONST_BITS);
  s[2] = _mm_srli_epi64(s[2], DCT_CONST_BITS);
  s[3] = _mm_srli_epi64(s[3], DCT_CONST_BITS);
  s[4] = _mm_srli_epi64(s[4], DCT_CONST_BITS);
  s[5] = _mm_srli_epi64(s[5], DCT_CONST_BITS);
  s[6] = _mm_srli_epi64(s[6], DCT_CONST_BITS);
  s[7] = _mm_srli_epi64(s[7], DCT_CONST_BITS);

  s[0] = _mm_shuffle_epi32(s[0], 0x88);
  s[1] = _mm_shuffle_epi32(s[1], 0x88);
  s[2] = _mm_shuffle_epi32(s[2], 0x88);
  s[3] = _mm_shuffle_epi32(s[3], 0x88);
  s[4] = _mm_shuffle_epi32(s[4], 0x88);
  s[5] = _mm_shuffle_epi32(s[5], 0x88);
  s[6] = _mm_shuffle_epi32(s[6], 0x88);
  s[7] = _mm_shuffle_epi32(s[7], 0x88);

  v0 = _mm_unpacklo_epi32(s[0], s[4]);
  v1 = _mm_unpacklo_epi32(s[2], s[6]);
  v2 = _mm_unpacklo_epi32(s[1], s[5]);
  v3 = _mm_unpacklo_epi32(s[3], s[7]);

  in[0] = _mm_unpacklo_epi64(v0, v1);
  in[1] = _mm_unpackhi_epi64(v0, v1);
  in[2] = _mm_unpacklo_epi64(v2, v3);
  in[3] = _mm_unpackhi_epi64(v2, v3);
}

static INLINE void write_buffer_4x4(tran_low_t *output, __m128i *res) {
  const __m128i kOne = _mm_set1_epi32(1);
  res[0] = _mm_add_epi32(res[0], kOne);
  res[1] = _mm_add_epi32(res[1], kOne);
  res[2] = _mm_add_epi32(res[2], kOne);
  res[3] = _mm_add_epi32(res[3], kOne);
  res[0] = _mm_srai_epi32(res[0], 2);
  res[1] = _mm_srai_epi32(res[1], 2);
  res[2] = _mm_srai_epi32(res[2], 2);
  res[3] = _mm_srai_epi32(res[3], 2);
  _mm_store_si128((__m128i *)(output + 0 * 4), res[0]);
  _mm_store_si128((__m128i *)(output + 1 * 4), res[1]);
  _mm_store_si128((__m128i *)(output + 2 * 4), res[2]);
  _mm_store_si128((__m128i *)(output + 3 * 4), res[3]);
}

void vp10_highbd_fht4x4_sse4_1(const int16_t *input, tran_low_t *output,
                               int stride, int tx_type) {
  __m128i in[4];

  switch (tx_type) {
    case DCT_DCT:
      load_buffer_4x4(input, in, stride, 0, 0);
      fdct4x4_sse4_1(in);
      fdct4x4_sse4_1(in);
      write_buffer_4x4(output, in);
      break;
    case ADST_DCT:
    case DCT_ADST:
    case ADST_ADST:
      vp10_highbd_fht4x4_c(input, output, stride, tx_type);
      break;
#if CONFIG_EXT_TX
    case FLIPADST_DCT:
    case DCT_FLIPADST:
    case FLIPADST_FLIPADST:
    case ADST_FLIPADST:
    case FLIPADST_ADST:
      vp10_highbd_fht4x4_c(input, output, stride, tx_type);
      break;
    case V_DCT:
    case H_DCT:
    case V_ADST:
    case H_ADST:
    case V_FLIPADST:
    case H_FLIPADST:
      vp10_highbd_fht4x4_c(input, output, stride, tx_type);
      break;
#endif  // CONFIG_EXT_TX
    default:
      assert(0);
  }
}
