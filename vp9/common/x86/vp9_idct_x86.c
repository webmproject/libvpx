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

void vp9_short_idct4x4_sse2(int16_t *input, int16_t *output, int pitch) {
  const __m128i zero = _mm_setzero_si128();
  const __m128i eight = _mm_set1_epi16(8);
  const __m128i cst = _mm_setr_epi16((int16_t)cospi_16_64, (int16_t)cospi_16_64,
                                    (int16_t)cospi_16_64, (int16_t)-cospi_16_64,
                                    (int16_t)cospi_24_64, (int16_t)-cospi_8_64,
                                    (int16_t)cospi_8_64, (int16_t)cospi_24_64);
  const __m128i rounding = _mm_set1_epi32(DCT_CONST_ROUNDING);
  const int half_pitch = pitch >> 1;
  __m128i input0, input1, input2, input3;

  // Rows
  input0 = _mm_loadl_epi64((__m128i *)input);
  input1 = _mm_loadl_epi64((__m128i *)(input + 4));
  input2 = _mm_loadl_epi64((__m128i *)(input + 8));
  input3 = _mm_loadl_epi64((__m128i *)(input + 12));

  // Construct i3, i1, i3, i1, i2, i0, i2, i0
  input0 = _mm_shufflelo_epi16(input0, 0xd8);
  input1 = _mm_shufflelo_epi16(input1, 0xd8);
  input2 = _mm_shufflelo_epi16(input2, 0xd8);
  input3 = _mm_shufflelo_epi16(input3, 0xd8);

  input0 = _mm_unpacklo_epi32(input0, input0);
  input1 = _mm_unpacklo_epi32(input1, input1);
  input2 = _mm_unpacklo_epi32(input2, input2);
  input3 = _mm_unpacklo_epi32(input3, input3);

  // Stage 1
  input0 = _mm_madd_epi16(input0, cst);
  input1 = _mm_madd_epi16(input1, cst);
  input2 = _mm_madd_epi16(input2, cst);
  input3 = _mm_madd_epi16(input3, cst);

  input0 = _mm_add_epi32(input0, rounding);
  input1 = _mm_add_epi32(input1, rounding);
  input2 = _mm_add_epi32(input2, rounding);
  input3 = _mm_add_epi32(input3, rounding);

  input0 = _mm_srai_epi32(input0, DCT_CONST_BITS);
  input1 = _mm_srai_epi32(input1, DCT_CONST_BITS);
  input2 = _mm_srai_epi32(input2, DCT_CONST_BITS);
  input3 = _mm_srai_epi32(input3, DCT_CONST_BITS);

  // Stage 2
  input0 = _mm_packs_epi32(input0, zero);
  input1 = _mm_packs_epi32(input1, zero);
  input2 = _mm_packs_epi32(input2, zero);
  input3 = _mm_packs_epi32(input3, zero);

  // Transpose
  input1 = _mm_unpacklo_epi16(input0, input1);
  input3 = _mm_unpacklo_epi16(input2, input3);
  input0 = _mm_unpacklo_epi32(input1, input3);
  input1 = _mm_unpackhi_epi32(input1, input3);

  // Switch column2, column 3, and then, we got:
  // input2: column1, column 0;  input3: column2, column 3.
  input1 = _mm_shuffle_epi32(input1, 0x4e);
  input2 = _mm_add_epi16(input0, input1);
  input3 = _mm_sub_epi16(input0, input1);

  // Columns
  // Construct i3, i1, i3, i1, i2, i0, i2, i0
  input0 = _mm_shufflelo_epi16(input2, 0xd8);
  input1 = _mm_shufflehi_epi16(input2, 0xd8);
  input2 = _mm_shufflehi_epi16(input3, 0xd8);
  input3 = _mm_shufflelo_epi16(input3, 0xd8);

  input0 = _mm_unpacklo_epi32(input0, input0);
  input1 = _mm_unpackhi_epi32(input1, input1);
  input2 = _mm_unpackhi_epi32(input2, input2);
  input3 = _mm_unpacklo_epi32(input3, input3);

  // Stage 1
  input0 = _mm_madd_epi16(input0, cst);
  input1 = _mm_madd_epi16(input1, cst);
  input2 = _mm_madd_epi16(input2, cst);
  input3 = _mm_madd_epi16(input3, cst);

  input0 = _mm_add_epi32(input0, rounding);
  input1 = _mm_add_epi32(input1, rounding);
  input2 = _mm_add_epi32(input2, rounding);
  input3 = _mm_add_epi32(input3, rounding);

  input0 = _mm_srai_epi32(input0, DCT_CONST_BITS);
  input1 = _mm_srai_epi32(input1, DCT_CONST_BITS);
  input2 = _mm_srai_epi32(input2, DCT_CONST_BITS);
  input3 = _mm_srai_epi32(input3, DCT_CONST_BITS);

  // Stage 2
  input0 = _mm_packs_epi32(input0, zero);
  input1 = _mm_packs_epi32(input1, zero);
  input2 = _mm_packs_epi32(input2, zero);
  input3 = _mm_packs_epi32(input3, zero);

  // Transpose
  input1 = _mm_unpacklo_epi16(input0, input1);
  input3 = _mm_unpacklo_epi16(input2, input3);
  input0 = _mm_unpacklo_epi32(input1, input3);
  input1 = _mm_unpackhi_epi32(input1, input3);

  // Switch column2, column 3, and then, we got:
  // input2: column1, column 0;  input3: column2, column 3.
  input1 = _mm_shuffle_epi32(input1, 0x4e);
  input2 = _mm_add_epi16(input0, input1);
  input3 = _mm_sub_epi16(input0, input1);

  // Final round and shift
  input2 = _mm_add_epi16(input2, eight);
  input3 = _mm_add_epi16(input3, eight);

  input2 = _mm_srai_epi16(input2, 4);
  input3 = _mm_srai_epi16(input3, 4);

  // Store results
  _mm_storel_epi64((__m128i *)output, input2);
  input2 = _mm_srli_si128(input2, 8);
  _mm_storel_epi64((__m128i *)(output + half_pitch), input2);

  _mm_storel_epi64((__m128i *)(output + 3 * half_pitch), input3);
  input3 = _mm_srli_si128(input3, 8);
  _mm_storel_epi64((__m128i *)(output + 2 * half_pitch), input3);
}

void vp9_idct4_1d_sse2(int16_t *input, int16_t *output) {
  const __m128i zero = _mm_setzero_si128();
  const __m128i c1 = _mm_setr_epi16((int16_t)cospi_16_64, (int16_t)cospi_16_64,
                                    (int16_t)cospi_16_64, (int16_t)-cospi_16_64,
                                    (int16_t)cospi_24_64, (int16_t)-cospi_8_64,
                                    (int16_t)cospi_8_64, (int16_t)cospi_24_64);
  const __m128i c2 = _mm_setr_epi16(1, 1, 1, 1, 1, -1, 1, -1);

  const __m128i rounding = _mm_set1_epi32(DCT_CONST_ROUNDING);
  __m128i in, temp;

  // Load input data.
  in = _mm_loadl_epi64((__m128i *)input);

  // Construct i3, i1, i3, i1, i2, i0, i2, i0
  in = _mm_shufflelo_epi16(in, 0xd8);
  in = _mm_unpacklo_epi32(in, in);

  // Stage 1
  in = _mm_madd_epi16(in, c1);
  in = _mm_add_epi32(in, rounding);
  in = _mm_srai_epi32(in, DCT_CONST_BITS);
  in = _mm_packs_epi32(in, zero);

  // Stage 2
  temp = _mm_shufflelo_epi16(in, 0x9c);
  in = _mm_shufflelo_epi16(in, 0xc9);
  in = _mm_unpacklo_epi64(temp, in);
  in = _mm_madd_epi16(in, c2);
  in = _mm_packs_epi32(in, zero);

  // Store results
  _mm_storel_epi64((__m128i *)output, in);
}


#define IDCT8x8_1D                                             \
  /* Stage1 */                                                 \
  {                                                            \
    const __m128i lo_17 = _mm_unpacklo_epi16(in1, in7);        \
    const __m128i hi_17 = _mm_unpackhi_epi16(in1, in7);        \
    const __m128i lo_35 = _mm_unpacklo_epi16(in3, in5);        \
    const __m128i hi_35 = _mm_unpackhi_epi16(in3, in5);        \
                                                               \
    tmp0 = _mm_madd_epi16(lo_17, stg1_0);                      \
    tmp1 = _mm_madd_epi16(hi_17, stg1_0);                      \
    tmp2 = _mm_madd_epi16(lo_17, stg1_1);                      \
    tmp3 = _mm_madd_epi16(hi_17, stg1_1);                      \
    tmp4 = _mm_madd_epi16(lo_35, stg1_2);                      \
    tmp5 = _mm_madd_epi16(hi_35, stg1_2);                      \
    tmp6 = _mm_madd_epi16(lo_35, stg1_3);                      \
    tmp7 = _mm_madd_epi16(hi_35, stg1_3);                      \
                                                               \
    tmp0 = _mm_add_epi32(tmp0, rounding);                      \
    tmp1 = _mm_add_epi32(tmp1, rounding);                      \
    tmp2 = _mm_add_epi32(tmp2, rounding);                      \
    tmp3 = _mm_add_epi32(tmp3, rounding);                      \
    tmp4 = _mm_add_epi32(tmp4, rounding);                      \
    tmp5 = _mm_add_epi32(tmp5, rounding);                      \
    tmp6 = _mm_add_epi32(tmp6, rounding);                      \
    tmp7 = _mm_add_epi32(tmp7, rounding);                      \
                                                               \
    tmp0 = _mm_srai_epi32(tmp0, DCT_CONST_BITS);               \
    tmp1 = _mm_srai_epi32(tmp1, DCT_CONST_BITS);               \
    tmp2 = _mm_srai_epi32(tmp2, DCT_CONST_BITS);               \
    tmp3 = _mm_srai_epi32(tmp3, DCT_CONST_BITS);               \
    tmp4 = _mm_srai_epi32(tmp4, DCT_CONST_BITS);               \
    tmp5 = _mm_srai_epi32(tmp5, DCT_CONST_BITS);               \
    tmp6 = _mm_srai_epi32(tmp6, DCT_CONST_BITS);               \
    tmp7 = _mm_srai_epi32(tmp7, DCT_CONST_BITS);               \
                                                               \
    stp1_4 = _mm_packs_epi32(tmp0, tmp1);                      \
    stp1_7 = _mm_packs_epi32(tmp2, tmp3);                      \
    stp1_5 = _mm_packs_epi32(tmp4, tmp5);                      \
    stp1_6 = _mm_packs_epi32(tmp6, tmp7);                      \
  }                                                            \
                                                               \
  /* Stage2 */                                                 \
  {                                                            \
    const __m128i lo_04 = _mm_unpacklo_epi16(in0, in4);        \
    const __m128i hi_04 = _mm_unpackhi_epi16(in0, in4);        \
    const __m128i lo_26 = _mm_unpacklo_epi16(in2, in6);        \
    const __m128i hi_26 = _mm_unpackhi_epi16(in2, in6);        \
                                                               \
    tmp0 = _mm_madd_epi16(lo_04, stg2_0);                      \
    tmp1 = _mm_madd_epi16(hi_04, stg2_0);                      \
    tmp2 = _mm_madd_epi16(lo_04, stg2_1);                      \
    tmp3 = _mm_madd_epi16(hi_04, stg2_1);                      \
    tmp4 = _mm_madd_epi16(lo_26, stg2_2);                      \
    tmp5 = _mm_madd_epi16(hi_26, stg2_2);                      \
    tmp6 = _mm_madd_epi16(lo_26, stg2_3);                      \
    tmp7 = _mm_madd_epi16(hi_26, stg2_3);                      \
                                                               \
    tmp0 = _mm_add_epi32(tmp0, rounding);                      \
    tmp1 = _mm_add_epi32(tmp1, rounding);                      \
    tmp2 = _mm_add_epi32(tmp2, rounding);                      \
    tmp3 = _mm_add_epi32(tmp3, rounding);                      \
    tmp4 = _mm_add_epi32(tmp4, rounding);                      \
    tmp5 = _mm_add_epi32(tmp5, rounding);                      \
    tmp6 = _mm_add_epi32(tmp6, rounding);                      \
    tmp7 = _mm_add_epi32(tmp7, rounding);                      \
                                                               \
    tmp0 = _mm_srai_epi32(tmp0, DCT_CONST_BITS);               \
    tmp1 = _mm_srai_epi32(tmp1, DCT_CONST_BITS);               \
    tmp2 = _mm_srai_epi32(tmp2, DCT_CONST_BITS);               \
    tmp3 = _mm_srai_epi32(tmp3, DCT_CONST_BITS);               \
    tmp4 = _mm_srai_epi32(tmp4, DCT_CONST_BITS);               \
    tmp5 = _mm_srai_epi32(tmp5, DCT_CONST_BITS);               \
    tmp6 = _mm_srai_epi32(tmp6, DCT_CONST_BITS);               \
    tmp7 = _mm_srai_epi32(tmp7, DCT_CONST_BITS);               \
                                                               \
    stp2_0 = _mm_packs_epi32(tmp0, tmp1);                      \
    stp2_1 = _mm_packs_epi32(tmp2, tmp3);                      \
    stp2_2 = _mm_packs_epi32(tmp4, tmp5);                      \
    stp2_3 = _mm_packs_epi32(tmp6, tmp7);                      \
                                                               \
    stp2_4 = _mm_adds_epi16(stp1_4, stp1_5);                   \
    stp2_5 = _mm_subs_epi16(stp1_4, stp1_5);                   \
    stp2_6 = _mm_subs_epi16(stp1_7, stp1_6);                   \
    stp2_7 = _mm_adds_epi16(stp1_7, stp1_6);                   \
  }                                                            \
                                                               \
  /* Stage3 */                                                 \
  {                                                            \
    const __m128i lo_56 = _mm_unpacklo_epi16(stp2_6, stp2_5);  \
    const __m128i hi_56 = _mm_unpackhi_epi16(stp2_6, stp2_5);  \
                                                               \
    stp1_0 = _mm_adds_epi16(stp2_0, stp2_3);                   \
    stp1_1 = _mm_adds_epi16(stp2_1, stp2_2);                   \
    stp1_2 = _mm_subs_epi16(stp2_1, stp2_2);                   \
    stp1_3 = _mm_subs_epi16(stp2_0, stp2_3);                   \
                                                               \
    tmp0 = _mm_madd_epi16(lo_56, stg2_1);                      \
    tmp1 = _mm_madd_epi16(hi_56, stg2_1);                      \
    tmp2 = _mm_madd_epi16(lo_56, stg2_0);                      \
    tmp3 = _mm_madd_epi16(hi_56, stg2_0);                      \
                                                               \
    tmp0 = _mm_add_epi32(tmp0, rounding);                      \
    tmp1 = _mm_add_epi32(tmp1, rounding);                      \
    tmp2 = _mm_add_epi32(tmp2, rounding);                      \
    tmp3 = _mm_add_epi32(tmp3, rounding);                      \
                                                               \
    tmp0 = _mm_srai_epi32(tmp0, DCT_CONST_BITS);               \
    tmp1 = _mm_srai_epi32(tmp1, DCT_CONST_BITS);               \
    tmp2 = _mm_srai_epi32(tmp2, DCT_CONST_BITS);               \
    tmp3 = _mm_srai_epi32(tmp3, DCT_CONST_BITS);               \
                                                               \
    stp1_5 = _mm_packs_epi32(tmp0, tmp1);                      \
    stp1_6 = _mm_packs_epi32(tmp2, tmp3);                      \
  }                                                            \
                                                               \
  /* Stage4  */                                                \
  in0 = _mm_adds_epi16(stp1_0, stp2_7);                        \
  in1 = _mm_adds_epi16(stp1_1, stp1_6);                        \
  in2 = _mm_adds_epi16(stp1_2, stp1_5);                        \
  in3 = _mm_adds_epi16(stp1_3, stp2_4);                        \
  in4 = _mm_subs_epi16(stp1_3, stp2_4);                        \
  in5 = _mm_subs_epi16(stp1_2, stp1_5);                        \
  in6 = _mm_subs_epi16(stp1_1, stp1_6);                        \
  in7 = _mm_subs_epi16(stp1_0, stp2_7);

void vp9_short_idct8x8_sse2(int16_t *input, int16_t *output, int pitch) {
  const int half_pitch = pitch >> 1;
  const __m128i rounding = _mm_set1_epi32(DCT_CONST_ROUNDING);
  const __m128i final_rounding = _mm_set1_epi16(1<<4);
  const __m128i stg1_0 = pair_set_epi16(cospi_28_64, -cospi_4_64);
  const __m128i stg1_1 = pair_set_epi16(cospi_4_64, cospi_28_64);
  const __m128i stg1_2 = pair_set_epi16(-cospi_20_64, cospi_12_64);
  const __m128i stg1_3 = pair_set_epi16(cospi_12_64, cospi_20_64);
  const __m128i stg2_0 = pair_set_epi16(cospi_16_64, cospi_16_64);
  const __m128i stg2_1 = pair_set_epi16(cospi_16_64, -cospi_16_64);
  const __m128i stg2_2 = pair_set_epi16(cospi_24_64, -cospi_8_64);
  const __m128i stg2_3 = pair_set_epi16(cospi_8_64, cospi_24_64);

  __m128i in0, in1, in2, in3, in4, in5, in6, in7;
  __m128i stp1_0, stp1_1, stp1_2, stp1_3, stp1_4, stp1_5, stp1_6, stp1_7;
  __m128i stp2_0, stp2_1, stp2_2, stp2_3, stp2_4, stp2_5, stp2_6, stp2_7;
  __m128i tmp0, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7;
  int i;

  // Load input data.
  in0 = _mm_load_si128((__m128i *)input);
  in1 = _mm_load_si128((__m128i *)(input + 8 * 1));
  in2 = _mm_load_si128((__m128i *)(input + 8 * 2));
  in3 = _mm_load_si128((__m128i *)(input + 8 * 3));
  in4 = _mm_load_si128((__m128i *)(input + 8 * 4));
  in5 = _mm_load_si128((__m128i *)(input + 8 * 5));
  in6 = _mm_load_si128((__m128i *)(input + 8 * 6));
  in7 = _mm_load_si128((__m128i *)(input + 8 * 7));

  // 2-D
  for (i = 0; i < 2; i++) {
    // 8x8 Transpose is copied from vp9_short_fdct8x8_sse2()
    {
      const __m128i tr0_0 = _mm_unpacklo_epi16(in0, in1);
      const __m128i tr0_1 = _mm_unpacklo_epi16(in2, in3);
      const __m128i tr0_2 = _mm_unpackhi_epi16(in0, in1);
      const __m128i tr0_3 = _mm_unpackhi_epi16(in2, in3);
      const __m128i tr0_4 = _mm_unpacklo_epi16(in4, in5);
      const __m128i tr0_5 = _mm_unpacklo_epi16(in6, in7);
      const __m128i tr0_6 = _mm_unpackhi_epi16(in4, in5);
      const __m128i tr0_7 = _mm_unpackhi_epi16(in6, in7);

      const __m128i tr1_0 = _mm_unpacklo_epi32(tr0_0, tr0_1);
      const __m128i tr1_1 = _mm_unpacklo_epi32(tr0_2, tr0_3);
      const __m128i tr1_2 = _mm_unpackhi_epi32(tr0_0, tr0_1);
      const __m128i tr1_3 = _mm_unpackhi_epi32(tr0_2, tr0_3);
      const __m128i tr1_4 = _mm_unpacklo_epi32(tr0_4, tr0_5);
      const __m128i tr1_5 = _mm_unpacklo_epi32(tr0_6, tr0_7);
      const __m128i tr1_6 = _mm_unpackhi_epi32(tr0_4, tr0_5);
      const __m128i tr1_7 = _mm_unpackhi_epi32(tr0_6, tr0_7);

      in0 = _mm_unpacklo_epi64(tr1_0, tr1_4);
      in1 = _mm_unpackhi_epi64(tr1_0, tr1_4);
      in2 = _mm_unpacklo_epi64(tr1_2, tr1_6);
      in3 = _mm_unpackhi_epi64(tr1_2, tr1_6);
      in4 = _mm_unpacklo_epi64(tr1_1, tr1_5);
      in5 = _mm_unpackhi_epi64(tr1_1, tr1_5);
      in6 = _mm_unpacklo_epi64(tr1_3, tr1_7);
      in7 = _mm_unpackhi_epi64(tr1_3, tr1_7);
    }

    // 4-stage 1D idct8x8
    IDCT8x8_1D
  }

  // Final rounding and shift
  in0 = _mm_add_epi16(in0, final_rounding);
  in1 = _mm_add_epi16(in1, final_rounding);
  in2 = _mm_add_epi16(in2, final_rounding);
  in3 = _mm_add_epi16(in3, final_rounding);
  in4 = _mm_add_epi16(in4, final_rounding);
  in5 = _mm_add_epi16(in5, final_rounding);
  in6 = _mm_add_epi16(in6, final_rounding);
  in7 = _mm_add_epi16(in7, final_rounding);

  in0 = _mm_srai_epi16(in0, 5);
  in1 = _mm_srai_epi16(in1, 5);
  in2 = _mm_srai_epi16(in2, 5);
  in3 = _mm_srai_epi16(in3, 5);
  in4 = _mm_srai_epi16(in4, 5);
  in5 = _mm_srai_epi16(in5, 5);
  in6 = _mm_srai_epi16(in6, 5);
  in7 = _mm_srai_epi16(in7, 5);

  // Store results
  _mm_store_si128((__m128i *)output, in0);
  _mm_store_si128((__m128i *)(output + half_pitch * 1), in1);
  _mm_store_si128((__m128i *)(output + half_pitch * 2), in2);
  _mm_store_si128((__m128i *)(output + half_pitch * 3), in3);
  _mm_store_si128((__m128i *)(output + half_pitch * 4), in4);
  _mm_store_si128((__m128i *)(output + half_pitch * 5), in5);
  _mm_store_si128((__m128i *)(output + half_pitch * 6), in6);
  _mm_store_si128((__m128i *)(output + half_pitch * 7), in7);
}

void vp9_short_idct10_8x8_sse2(int16_t *input, int16_t *output, int pitch) {
  const int half_pitch = pitch >> 1;
  const __m128i zero = _mm_setzero_si128();
  const __m128i rounding = _mm_set1_epi32(DCT_CONST_ROUNDING);
  const __m128i final_rounding = _mm_set1_epi16(1<<4);
  const __m128i stg1_0 = pair_set_epi16(cospi_28_64, -cospi_4_64);
  const __m128i stg1_1 = pair_set_epi16(cospi_4_64, cospi_28_64);
  const __m128i stg1_2 = pair_set_epi16(-cospi_20_64, cospi_12_64);
  const __m128i stg1_3 = pair_set_epi16(cospi_12_64, cospi_20_64);
  const __m128i stg2_0 = pair_set_epi16(cospi_16_64, cospi_16_64);
  const __m128i stg2_1 = pair_set_epi16(cospi_16_64, -cospi_16_64);
  const __m128i stg2_2 = pair_set_epi16(cospi_24_64, -cospi_8_64);
  const __m128i stg2_3 = pair_set_epi16(cospi_8_64, cospi_24_64);
  const __m128i stg3_0 = pair_set_epi16(-cospi_16_64, cospi_16_64);

  __m128i in0, in1, in2, in3, in4, in5, in6, in7;
  __m128i stp1_0, stp1_1, stp1_2, stp1_3, stp1_4, stp1_5, stp1_6, stp1_7;
  __m128i stp2_0, stp2_1, stp2_2, stp2_3, stp2_4, stp2_5, stp2_6, stp2_7;
  __m128i tmp0, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7;

  // Rows. Load 4-row input data.
  in0 = _mm_load_si128((__m128i *)input);
  in1 = _mm_load_si128((__m128i *)(input + 8 * 1));
  in2 = _mm_load_si128((__m128i *)(input + 8 * 2));
  in3 = _mm_load_si128((__m128i *)(input + 8 * 3));

  // 8x4 Transpose
  {
    const __m128i tr0_0 = _mm_unpacklo_epi16(in0, in1);
    const __m128i tr0_1 = _mm_unpacklo_epi16(in2, in3);
    const __m128i tr0_2 = _mm_unpackhi_epi16(in0, in1);
    const __m128i tr0_3 = _mm_unpackhi_epi16(in2, in3);

    in0 = _mm_unpacklo_epi32(tr0_0, tr0_1);  // i1 i0
    in1 = _mm_unpacklo_epi32(tr0_2, tr0_3);  // i5 i4
    in2 = _mm_unpackhi_epi32(tr0_0, tr0_1);  // i3 i2
    in3 = _mm_unpackhi_epi32(tr0_2, tr0_3);  // i7 i6
  }

  // Stage1
  {
    const __m128i lo_17 = _mm_unpackhi_epi16(in0, in3);
    const __m128i lo_35 = _mm_unpackhi_epi16(in2, in1);

    tmp0 = _mm_madd_epi16(lo_17, stg1_0);
    tmp2 = _mm_madd_epi16(lo_17, stg1_1);
    tmp4 = _mm_madd_epi16(lo_35, stg1_2);
    tmp6 = _mm_madd_epi16(lo_35, stg1_3);

    tmp0 = _mm_add_epi32(tmp0, rounding);
    tmp2 = _mm_add_epi32(tmp2, rounding);
    tmp4 = _mm_add_epi32(tmp4, rounding);
    tmp6 = _mm_add_epi32(tmp6, rounding);
    tmp0 = _mm_srai_epi32(tmp0, DCT_CONST_BITS);
    tmp2 = _mm_srai_epi32(tmp2, DCT_CONST_BITS);
    tmp4 = _mm_srai_epi32(tmp4, DCT_CONST_BITS);
    tmp6 = _mm_srai_epi32(tmp6, DCT_CONST_BITS);

    stp1_4 = _mm_packs_epi32(tmp0, zero);
    stp1_7 = _mm_packs_epi32(tmp2, zero);
    stp1_5 = _mm_packs_epi32(tmp4, zero);
    stp1_6 = _mm_packs_epi32(tmp6, zero);
  }

  // Stage2
  {
    const __m128i lo_04 = _mm_unpacklo_epi16(in0, in1);
    const __m128i lo_26 = _mm_unpacklo_epi16(in2, in3);

    tmp0 = _mm_madd_epi16(lo_04, stg2_0);
    tmp2 = _mm_madd_epi16(lo_04, stg2_1);
    tmp4 = _mm_madd_epi16(lo_26, stg2_2);
    tmp6 = _mm_madd_epi16(lo_26, stg2_3);

    tmp0 = _mm_add_epi32(tmp0, rounding);
    tmp2 = _mm_add_epi32(tmp2, rounding);
    tmp4 = _mm_add_epi32(tmp4, rounding);
    tmp6 = _mm_add_epi32(tmp6, rounding);
    tmp0 = _mm_srai_epi32(tmp0, DCT_CONST_BITS);
    tmp2 = _mm_srai_epi32(tmp2, DCT_CONST_BITS);
    tmp4 = _mm_srai_epi32(tmp4, DCT_CONST_BITS);
    tmp6 = _mm_srai_epi32(tmp6, DCT_CONST_BITS);

    stp2_0 = _mm_packs_epi32(tmp0, zero);
    stp2_1 = _mm_packs_epi32(tmp2, zero);
    stp2_2 = _mm_packs_epi32(tmp4, zero);
    stp2_3 = _mm_packs_epi32(tmp6, zero);

    stp2_4 = _mm_adds_epi16(stp1_4, stp1_5);
    stp2_5 = _mm_subs_epi16(stp1_4, stp1_5);
    stp2_6 = _mm_subs_epi16(stp1_7, stp1_6);
    stp2_7 = _mm_adds_epi16(stp1_7, stp1_6);
  }

  // Stage3
  {
    const __m128i lo_56 = _mm_unpacklo_epi16(stp2_5, stp2_6);
    stp1_0 = _mm_adds_epi16(stp2_0, stp2_3);
    stp1_1 = _mm_adds_epi16(stp2_1, stp2_2);
    stp1_2 = _mm_subs_epi16(stp2_1, stp2_2);
    stp1_3 = _mm_subs_epi16(stp2_0, stp2_3);

    tmp0 = _mm_madd_epi16(lo_56, stg3_0);
    tmp2 = _mm_madd_epi16(lo_56, stg2_0);  // stg3_1 = stg2_0

    tmp0 = _mm_add_epi32(tmp0, rounding);
    tmp2 = _mm_add_epi32(tmp2, rounding);
    tmp0 = _mm_srai_epi32(tmp0, DCT_CONST_BITS);
    tmp2 = _mm_srai_epi32(tmp2, DCT_CONST_BITS);

    stp1_5 = _mm_packs_epi32(tmp0, zero);
    stp1_6 = _mm_packs_epi32(tmp2, zero);
  }

  // Stage4
  in0 = _mm_adds_epi16(stp1_0, stp2_7);
  in1 = _mm_adds_epi16(stp1_1, stp1_6);
  in2 = _mm_adds_epi16(stp1_2, stp1_5);
  in3 = _mm_adds_epi16(stp1_3, stp2_4);
  in4 = _mm_subs_epi16(stp1_3, stp2_4);
  in5 = _mm_subs_epi16(stp1_2, stp1_5);
  in6 = _mm_subs_epi16(stp1_1, stp1_6);
  in7 = _mm_subs_epi16(stp1_0, stp2_7);

  // Columns. 4x8 Transpose
  {
    const __m128i tr0_0 = _mm_unpacklo_epi16(in0, in1);
    const __m128i tr0_1 = _mm_unpacklo_epi16(in2, in3);
    const __m128i tr0_4 = _mm_unpacklo_epi16(in4, in5);
    const __m128i tr0_5 = _mm_unpacklo_epi16(in6, in7);

    const __m128i tr1_0 = _mm_unpacklo_epi32(tr0_0, tr0_1);
    const __m128i tr1_2 = _mm_unpackhi_epi32(tr0_0, tr0_1);
    const __m128i tr1_4 = _mm_unpacklo_epi32(tr0_4, tr0_5);
    const __m128i tr1_6 = _mm_unpackhi_epi32(tr0_4, tr0_5);

    in0 = _mm_unpacklo_epi64(tr1_0, tr1_4);
    in1 = _mm_unpackhi_epi64(tr1_0, tr1_4);
    in2 = _mm_unpacklo_epi64(tr1_2, tr1_6);
    in3 = _mm_unpackhi_epi64(tr1_2, tr1_6);
    in4 = _mm_setzero_si128();
    in5 = _mm_setzero_si128();
    in6 = _mm_setzero_si128();
    in7 = _mm_setzero_si128();
  }

  // 1D idct8x8
  IDCT8x8_1D

  // Final rounding and shift
  in0 = _mm_add_epi16(in0, final_rounding);
  in1 = _mm_add_epi16(in1, final_rounding);
  in2 = _mm_add_epi16(in2, final_rounding);
  in3 = _mm_add_epi16(in3, final_rounding);
  in4 = _mm_add_epi16(in4, final_rounding);
  in5 = _mm_add_epi16(in5, final_rounding);
  in6 = _mm_add_epi16(in6, final_rounding);
  in7 = _mm_add_epi16(in7, final_rounding);

  in0 = _mm_srai_epi16(in0, 5);
  in1 = _mm_srai_epi16(in1, 5);
  in2 = _mm_srai_epi16(in2, 5);
  in3 = _mm_srai_epi16(in3, 5);
  in4 = _mm_srai_epi16(in4, 5);
  in5 = _mm_srai_epi16(in5, 5);
  in6 = _mm_srai_epi16(in6, 5);
  in7 = _mm_srai_epi16(in7, 5);

  // Store results
  _mm_store_si128((__m128i *)output, in0);
  _mm_store_si128((__m128i *)(output + half_pitch * 1), in1);
  _mm_store_si128((__m128i *)(output + half_pitch * 2), in2);
  _mm_store_si128((__m128i *)(output + half_pitch * 3), in3);
  _mm_store_si128((__m128i *)(output + half_pitch * 4), in4);
  _mm_store_si128((__m128i *)(output + half_pitch * 5), in5);
  _mm_store_si128((__m128i *)(output + half_pitch * 6), in6);
  _mm_store_si128((__m128i *)(output + half_pitch * 7), in7);
}
#endif
