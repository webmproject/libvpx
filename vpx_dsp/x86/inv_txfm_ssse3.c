/*
 *  Copyright (c) 2017 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include <tmmintrin.h>

#include "./vpx_dsp_rtcd.h"
#include "vpx_dsp/x86/inv_txfm_sse2.h"
#include "vpx_dsp/x86/txfm_common_sse2.h"

void vpx_idct8x8_64_add_ssse3(const tran_low_t *input, uint8_t *dest,
                              int stride) {
  const __m128i zero = _mm_setzero_si128();
  const __m128i rounding = _mm_set1_epi32(DCT_CONST_ROUNDING);
  const __m128i final_rounding = _mm_set1_epi16(1 << 4);
  const __m128i stg1_0 = pair_set_epi16(cospi_28_64, -cospi_4_64);
  const __m128i stg1_1 = pair_set_epi16(cospi_4_64, cospi_28_64);
  const __m128i stg1_2 = pair_set_epi16(-cospi_20_64, cospi_12_64);
  const __m128i stg1_3 = pair_set_epi16(cospi_12_64, cospi_20_64);
  const __m128i stg2_0 = pair_set_epi16(2 * cospi_16_64, 2 * cospi_16_64);
  const __m128i stg2_2 = pair_set_epi16(cospi_24_64, -cospi_8_64);
  const __m128i stg2_3 = pair_set_epi16(cospi_8_64, cospi_24_64);

  __m128i in0, in1, in2, in3, in4, in5, in6, in7;
  __m128i stp1_0, stp1_1, stp1_2, stp1_3, stp1_4, stp1_5, stp1_6, stp1_7;
  __m128i stp2_0, stp2_1, stp2_2, stp2_3, stp2_4, stp2_5, stp2_6, stp2_7;
  __m128i tmp0, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7;
  int i;

  // Load input data.
  in0 = load_input_data(input);
  in1 = load_input_data(input + 8 * 1);
  in2 = load_input_data(input + 8 * 2);
  in3 = load_input_data(input + 8 * 3);
  in4 = load_input_data(input + 8 * 4);
  in5 = load_input_data(input + 8 * 5);
  in6 = load_input_data(input + 8 * 6);
  in7 = load_input_data(input + 8 * 7);

  // 2-D
  for (i = 0; i < 2; i++) {
    // 8x8 Transpose is copied from vpx_fdct8x8_sse2()
    TRANSPOSE_8X8(in0, in1, in2, in3, in4, in5, in6, in7, in0, in1, in2, in3,
                  in4, in5, in6, in7);

    // 4-stage 1D idct8x8
    {
      /* Stage1 */
      {
        const __m128i lo_17 = _mm_unpacklo_epi16(in1, in7);
        const __m128i hi_17 = _mm_unpackhi_epi16(in1, in7);
        const __m128i lo_35 = _mm_unpacklo_epi16(in3, in5);
        const __m128i hi_35 = _mm_unpackhi_epi16(in3, in5);

        {
          tmp0 = _mm_madd_epi16(lo_17, stg1_0);
          tmp1 = _mm_madd_epi16(hi_17, stg1_0);
          tmp2 = _mm_madd_epi16(lo_17, stg1_1);
          tmp3 = _mm_madd_epi16(hi_17, stg1_1);
          tmp4 = _mm_madd_epi16(lo_35, stg1_2);
          tmp5 = _mm_madd_epi16(hi_35, stg1_2);
          tmp6 = _mm_madd_epi16(lo_35, stg1_3);
          tmp7 = _mm_madd_epi16(hi_35, stg1_3);

          tmp0 = _mm_add_epi32(tmp0, rounding);
          tmp1 = _mm_add_epi32(tmp1, rounding);
          tmp2 = _mm_add_epi32(tmp2, rounding);
          tmp3 = _mm_add_epi32(tmp3, rounding);
          tmp4 = _mm_add_epi32(tmp4, rounding);
          tmp5 = _mm_add_epi32(tmp5, rounding);
          tmp6 = _mm_add_epi32(tmp6, rounding);
          tmp7 = _mm_add_epi32(tmp7, rounding);

          tmp0 = _mm_srai_epi32(tmp0, 14);
          tmp1 = _mm_srai_epi32(tmp1, 14);
          tmp2 = _mm_srai_epi32(tmp2, 14);
          tmp3 = _mm_srai_epi32(tmp3, 14);
          tmp4 = _mm_srai_epi32(tmp4, 14);
          tmp5 = _mm_srai_epi32(tmp5, 14);
          tmp6 = _mm_srai_epi32(tmp6, 14);
          tmp7 = _mm_srai_epi32(tmp7, 14);

          stp1_4 = _mm_packs_epi32(tmp0, tmp1);
          stp1_7 = _mm_packs_epi32(tmp2, tmp3);
          stp1_5 = _mm_packs_epi32(tmp4, tmp5);
          stp1_6 = _mm_packs_epi32(tmp6, tmp7);
        }
      }

      /* Stage2 */
      {
        const __m128i lo_26 = _mm_unpacklo_epi16(in2, in6);
        const __m128i hi_26 = _mm_unpackhi_epi16(in2, in6);

        {
          tmp0 = _mm_add_epi16(in0, in4);
          tmp1 = _mm_sub_epi16(in0, in4);
          stp2_0 = _mm_mulhrs_epi16(tmp0, stg2_0);
          stp2_1 = _mm_mulhrs_epi16(tmp1, stg2_0);

          tmp0 = _mm_madd_epi16(lo_26, stg2_2);
          tmp1 = _mm_madd_epi16(hi_26, stg2_2);
          tmp2 = _mm_madd_epi16(lo_26, stg2_3);
          tmp3 = _mm_madd_epi16(hi_26, stg2_3);

          tmp0 = _mm_add_epi32(tmp0, rounding);
          tmp1 = _mm_add_epi32(tmp1, rounding);
          tmp2 = _mm_add_epi32(tmp2, rounding);
          tmp3 = _mm_add_epi32(tmp3, rounding);

          tmp0 = _mm_srai_epi32(tmp0, 14);
          tmp1 = _mm_srai_epi32(tmp1, 14);
          tmp2 = _mm_srai_epi32(tmp2, 14);
          tmp3 = _mm_srai_epi32(tmp3, 14);

          stp2_2 = _mm_packs_epi32(tmp0, tmp1);
          stp2_3 = _mm_packs_epi32(tmp2, tmp3);
        }

        stp2_4 = _mm_add_epi16(stp1_4, stp1_5);
        stp2_5 = _mm_sub_epi16(stp1_4, stp1_5);
        stp2_6 = _mm_sub_epi16(stp1_7, stp1_6);
        stp2_7 = _mm_add_epi16(stp1_7, stp1_6);
      }

      /* Stage3 */
      {
        stp1_0 = _mm_add_epi16(stp2_0, stp2_3);
        stp1_1 = _mm_add_epi16(stp2_1, stp2_2);
        stp1_2 = _mm_sub_epi16(stp2_1, stp2_2);
        stp1_3 = _mm_sub_epi16(stp2_0, stp2_3);

        tmp0 = _mm_sub_epi16(stp2_6, stp2_5);
        tmp2 = _mm_add_epi16(stp2_6, stp2_5);
        stp1_5 = _mm_mulhrs_epi16(tmp0, stg2_0);
        stp1_6 = _mm_mulhrs_epi16(tmp2, stg2_0);
      }

      /* Stage4  */
      in0 = _mm_add_epi16(stp1_0, stp2_7);
      in1 = _mm_add_epi16(stp1_1, stp1_6);
      in2 = _mm_add_epi16(stp1_2, stp1_5);
      in3 = _mm_add_epi16(stp1_3, stp2_4);
      in4 = _mm_sub_epi16(stp1_3, stp2_4);
      in5 = _mm_sub_epi16(stp1_2, stp1_5);
      in6 = _mm_sub_epi16(stp1_1, stp1_6);
      in7 = _mm_sub_epi16(stp1_0, stp2_7);
    }
  }

  // Final rounding and shift
  in0 = _mm_adds_epi16(in0, final_rounding);
  in1 = _mm_adds_epi16(in1, final_rounding);
  in2 = _mm_adds_epi16(in2, final_rounding);
  in3 = _mm_adds_epi16(in3, final_rounding);
  in4 = _mm_adds_epi16(in4, final_rounding);
  in5 = _mm_adds_epi16(in5, final_rounding);
  in6 = _mm_adds_epi16(in6, final_rounding);
  in7 = _mm_adds_epi16(in7, final_rounding);

  in0 = _mm_srai_epi16(in0, 5);
  in1 = _mm_srai_epi16(in1, 5);
  in2 = _mm_srai_epi16(in2, 5);
  in3 = _mm_srai_epi16(in3, 5);
  in4 = _mm_srai_epi16(in4, 5);
  in5 = _mm_srai_epi16(in5, 5);
  in6 = _mm_srai_epi16(in6, 5);
  in7 = _mm_srai_epi16(in7, 5);

  RECON_AND_STORE(dest + 0 * stride, in0);
  RECON_AND_STORE(dest + 1 * stride, in1);
  RECON_AND_STORE(dest + 2 * stride, in2);
  RECON_AND_STORE(dest + 3 * stride, in3);
  RECON_AND_STORE(dest + 4 * stride, in4);
  RECON_AND_STORE(dest + 5 * stride, in5);
  RECON_AND_STORE(dest + 6 * stride, in6);
  RECON_AND_STORE(dest + 7 * stride, in7);
}

void vpx_idct8x8_12_add_ssse3(const tran_low_t *input, uint8_t *dest,
                              int stride) {
  const __m128i zero = _mm_setzero_si128();
  const __m128i final_rounding = _mm_set1_epi16(1 << 4);
  const __m128i stg1_0 = pair_set_epi16(2 * cospi_28_64, 2 * cospi_28_64);
  const __m128i stg1_1 = pair_set_epi16(2 * cospi_4_64, 2 * cospi_4_64);
  const __m128i stg1_2 = pair_set_epi16(-2 * cospi_20_64, -2 * cospi_20_64);
  const __m128i stg1_3 = pair_set_epi16(2 * cospi_12_64, 2 * cospi_12_64);
  const __m128i stg2_0 = pair_set_epi16(2 * cospi_16_64, 2 * cospi_16_64);
  const __m128i stg2_2 = pair_set_epi16(2 * cospi_24_64, 2 * cospi_24_64);
  const __m128i stg2_3 = pair_set_epi16(2 * cospi_8_64, 2 * cospi_8_64);

  __m128i in0, in1, in2, in3, in4, in5, in6, in7;
  __m128i stp1_0, stp1_1, stp1_2, stp1_3, stp1_4, stp1_5, stp1_6, stp1_7;
  __m128i stp2_0, stp2_1, stp2_2, stp2_3, stp2_4, stp2_5, stp2_6, stp2_7;
  __m128i tmp0, tmp1, tmp2, tmp3;

  // Rows. Load 4-row input data.
  in0 = load_input_data(input);
  in1 = load_input_data(input + 8 * 1);
  in2 = load_input_data(input + 8 * 2);
  in3 = load_input_data(input + 8 * 3);

  // 8x4 Transpose
  TRANSPOSE_8X8_10(in0, in1, in2, in3, in0, in1);

  // Stage1
  tmp0 = _mm_mulhrs_epi16(in0, stg1_0);
  tmp1 = _mm_mulhrs_epi16(in0, stg1_1);
  tmp2 = _mm_mulhrs_epi16(in1, stg1_2);
  tmp3 = _mm_mulhrs_epi16(in1, stg1_3);

  stp1_4 = _mm_unpackhi_epi64(tmp0, tmp1);
  stp1_5 = _mm_unpackhi_epi64(tmp2, tmp3);

  // Stage2
  tmp0 = _mm_mulhrs_epi16(in0, stg2_0);
  stp2_0 = _mm_unpacklo_epi64(tmp0, tmp0);

  tmp1 = _mm_mulhrs_epi16(in1, stg2_2);
  tmp2 = _mm_mulhrs_epi16(in1, stg2_3);
  stp2_2 = _mm_unpacklo_epi64(tmp2, tmp1);

  tmp0 = _mm_add_epi16(stp1_4, stp1_5);
  tmp1 = _mm_sub_epi16(stp1_4, stp1_5);

  stp2_4 = tmp0;
  stp2_5 = _mm_unpacklo_epi64(tmp1, zero);
  stp2_6 = _mm_unpackhi_epi64(tmp1, zero);

  // Stage3
  tmp2 = _mm_add_epi16(stp2_0, stp2_2);
  tmp3 = _mm_sub_epi16(stp2_0, stp2_2);

  stp1_2 = _mm_unpackhi_epi64(tmp3, tmp2);
  stp1_3 = _mm_unpacklo_epi64(tmp3, tmp2);

  tmp0 = _mm_sub_epi16(stp2_6, stp2_5);
  tmp1 = _mm_add_epi16(stp2_6, stp2_5);

  tmp2 = _mm_mulhrs_epi16(tmp0, stg2_0);
  tmp3 = _mm_mulhrs_epi16(tmp1, stg2_0);
  stp1_5 = _mm_unpacklo_epi64(tmp2, tmp3);

  // Stage4
  tmp0 = _mm_add_epi16(stp1_3, stp2_4);
  tmp1 = _mm_add_epi16(stp1_2, stp1_5);
  tmp2 = _mm_sub_epi16(stp1_3, stp2_4);
  tmp3 = _mm_sub_epi16(stp1_2, stp1_5);

  TRANSPOSE_4X8_10(tmp0, tmp1, tmp2, tmp3, in0, in1, in2, in3)

  /* Stage1 */
  stp1_4 = _mm_mulhrs_epi16(in1, stg1_0);
  stp1_7 = _mm_mulhrs_epi16(in1, stg1_1);
  stp1_5 = _mm_mulhrs_epi16(in3, stg1_2);
  stp1_6 = _mm_mulhrs_epi16(in3, stg1_3);

  /* Stage2 */
  stp2_0 = _mm_mulhrs_epi16(in0, stg2_0);
  stp2_1 = _mm_mulhrs_epi16(in0, stg2_0);

  stp2_2 = _mm_mulhrs_epi16(in2, stg2_2);
  stp2_3 = _mm_mulhrs_epi16(in2, stg2_3);

  stp2_4 = _mm_add_epi16(stp1_4, stp1_5);
  stp2_5 = _mm_sub_epi16(stp1_4, stp1_5);
  stp2_6 = _mm_sub_epi16(stp1_7, stp1_6);
  stp2_7 = _mm_add_epi16(stp1_7, stp1_6);

  /* Stage3 */
  stp1_0 = _mm_add_epi16(stp2_0, stp2_3);
  stp1_1 = _mm_add_epi16(stp2_1, stp2_2);
  stp1_2 = _mm_sub_epi16(stp2_1, stp2_2);
  stp1_3 = _mm_sub_epi16(stp2_0, stp2_3);

  tmp0 = _mm_add_epi16(stp2_6, stp2_5);
  tmp1 = _mm_sub_epi16(stp2_6, stp2_5);
  stp1_6 = _mm_mulhrs_epi16(tmp0, stg2_0);
  stp1_5 = _mm_mulhrs_epi16(tmp1, stg2_0);

  /* Stage4  */
  in0 = _mm_add_epi16(stp1_0, stp2_7);
  in1 = _mm_add_epi16(stp1_1, stp1_6);
  in2 = _mm_add_epi16(stp1_2, stp1_5);
  in3 = _mm_add_epi16(stp1_3, stp2_4);
  in4 = _mm_sub_epi16(stp1_3, stp2_4);
  in5 = _mm_sub_epi16(stp1_2, stp1_5);
  in6 = _mm_sub_epi16(stp1_1, stp1_6);
  in7 = _mm_sub_epi16(stp1_0, stp2_7);

  // Final rounding and shift
  in0 = _mm_adds_epi16(in0, final_rounding);
  in1 = _mm_adds_epi16(in1, final_rounding);
  in2 = _mm_adds_epi16(in2, final_rounding);
  in3 = _mm_adds_epi16(in3, final_rounding);
  in4 = _mm_adds_epi16(in4, final_rounding);
  in5 = _mm_adds_epi16(in5, final_rounding);
  in6 = _mm_adds_epi16(in6, final_rounding);
  in7 = _mm_adds_epi16(in7, final_rounding);

  in0 = _mm_srai_epi16(in0, 5);
  in1 = _mm_srai_epi16(in1, 5);
  in2 = _mm_srai_epi16(in2, 5);
  in3 = _mm_srai_epi16(in3, 5);
  in4 = _mm_srai_epi16(in4, 5);
  in5 = _mm_srai_epi16(in5, 5);
  in6 = _mm_srai_epi16(in6, 5);
  in7 = _mm_srai_epi16(in7, 5);

  RECON_AND_STORE(dest + 0 * stride, in0);
  RECON_AND_STORE(dest + 1 * stride, in1);
  RECON_AND_STORE(dest + 2 * stride, in2);
  RECON_AND_STORE(dest + 3 * stride, in3);
  RECON_AND_STORE(dest + 4 * stride, in4);
  RECON_AND_STORE(dest + 5 * stride, in5);
  RECON_AND_STORE(dest + 6 * stride, in6);
  RECON_AND_STORE(dest + 7 * stride, in7);
}

static INLINE void idct32_34(const __m128i *in, __m128i *stp1) {
  const __m128i rounding = _mm_set1_epi32(DCT_CONST_ROUNDING);
  // idct constants for each stage
  const __m128i stk1_0 = pair_set_epi16(2 * cospi_31_64, 2 * cospi_31_64);
  const __m128i stk1_1 = pair_set_epi16(2 * cospi_1_64, 2 * cospi_1_64);
  const __m128i stk1_6 = pair_set_epi16(-2 * cospi_25_64, -2 * cospi_25_64);
  const __m128i stk1_7 = pair_set_epi16(2 * cospi_7_64, 2 * cospi_7_64);
  const __m128i stk1_8 = pair_set_epi16(2 * cospi_27_64, 2 * cospi_27_64);
  const __m128i stk1_9 = pair_set_epi16(2 * cospi_5_64, 2 * cospi_5_64);
  const __m128i stk1_14 = pair_set_epi16(-2 * cospi_29_64, -2 * cospi_29_64);
  const __m128i stk1_15 = pair_set_epi16(2 * cospi_3_64, 2 * cospi_3_64);

  const __m128i stk2_0 = pair_set_epi16(2 * cospi_30_64, 2 * cospi_30_64);
  const __m128i stk2_1 = pair_set_epi16(2 * cospi_2_64, 2 * cospi_2_64);
  const __m128i stk2_6 = pair_set_epi16(-2 * cospi_26_64, -2 * cospi_26_64);
  const __m128i stk2_7 = pair_set_epi16(2 * cospi_6_64, 2 * cospi_6_64);

  const __m128i stk3_0 = pair_set_epi16(2 * cospi_28_64, 2 * cospi_28_64);
  const __m128i stk3_1 = pair_set_epi16(2 * cospi_4_64, 2 * cospi_4_64);
  const __m128i stg3_4 = pair_set_epi16(-cospi_4_64, cospi_28_64);
  const __m128i stg3_5 = pair_set_epi16(cospi_28_64, cospi_4_64);
  const __m128i stg3_6 = pair_set_epi16(-cospi_28_64, -cospi_4_64);
  const __m128i stg3_8 = pair_set_epi16(-cospi_20_64, cospi_12_64);
  const __m128i stg3_9 = pair_set_epi16(cospi_12_64, cospi_20_64);
  const __m128i stg3_10 = pair_set_epi16(-cospi_12_64, -cospi_20_64);

  const __m128i stg4_0 = pair_set_epi16(cospi_16_64, cospi_16_64);
  const __m128i stk4_0 = pair_set_epi16(2 * cospi_16_64, 2 * cospi_16_64);
  const __m128i stg4_1 = pair_set_epi16(cospi_16_64, -cospi_16_64);
  const __m128i stg4_4 = pair_set_epi16(-cospi_8_64, cospi_24_64);
  const __m128i stg4_5 = pair_set_epi16(cospi_24_64, cospi_8_64);
  const __m128i stg4_6 = pair_set_epi16(-cospi_24_64, -cospi_8_64);

  const __m128i stg6_0 = pair_set_epi16(-cospi_16_64, cospi_16_64);
  __m128i stp2_0, stp2_1, stp2_2, stp2_3, stp2_4, stp2_5, stp2_6, stp2_7,
      stp2_8, stp2_9, stp2_10, stp2_11, stp2_12, stp2_13, stp2_14, stp2_15,
      stp2_16, stp2_17, stp2_18, stp2_19, stp2_20, stp2_21, stp2_22, stp2_23,
      stp2_24, stp2_25, stp2_26, stp2_27, stp2_28, stp2_29, stp2_30, stp2_31;
  __m128i tmp0, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7;

  /* Stage1 */

  stp1[16] = _mm_mulhrs_epi16(in[1], stk1_0);
  stp1[31] = _mm_mulhrs_epi16(in[1], stk1_1);

  stp1[19] = _mm_mulhrs_epi16(in[7], stk1_6);
  stp1[28] = _mm_mulhrs_epi16(in[7], stk1_7);

  stp1[20] = _mm_mulhrs_epi16(in[5], stk1_8);
  stp1[27] = _mm_mulhrs_epi16(in[5], stk1_9);

  stp1[23] = _mm_mulhrs_epi16(in[3], stk1_14);
  stp1[24] = _mm_mulhrs_epi16(in[3], stk1_15);

  /* Stage2 */

  stp2_8 = _mm_mulhrs_epi16(in[2], stk2_0);
  stp2_15 = _mm_mulhrs_epi16(in[2], stk2_1);

  stp2_11 = _mm_mulhrs_epi16(in[6], stk2_6);
  stp2_12 = _mm_mulhrs_epi16(in[6], stk2_7);

  /* Stage3 */
  {
    const __m128i lo_17_30 = _mm_unpacklo_epi16(stp1[16], stp1[31]);
    const __m128i hi_17_30 = _mm_unpackhi_epi16(stp1[16], stp1[31]);
    const __m128i lo_18_29 = _mm_unpacklo_epi16(stp1[19], stp1[28]);
    const __m128i hi_18_29 = _mm_unpackhi_epi16(stp1[19], stp1[28]);

    const __m128i lo_21_26 = _mm_unpacklo_epi16(stp1[20], stp1[27]);
    const __m128i hi_21_26 = _mm_unpackhi_epi16(stp1[20], stp1[27]);
    const __m128i lo_22_25 = _mm_unpacklo_epi16(stp1[23], stp1[24]);
    const __m128i hi_22_25 = _mm_unpackhi_epi16(stp1[23], stp1[24]);

    stp1[4] = _mm_mulhrs_epi16(in[4], stk3_0);
    stp1[7] = _mm_mulhrs_epi16(in[4], stk3_1);

    MULTIPLICATION_AND_ADD(lo_17_30, hi_17_30, lo_18_29, hi_18_29, stg3_4,
                           stg3_5, stg3_6, stg3_4, stp1[17], stp1[30], stp1[18],
                           stp1[29])
    MULTIPLICATION_AND_ADD(lo_21_26, hi_21_26, lo_22_25, hi_22_25, stg3_8,
                           stg3_9, stg3_10, stg3_8, stp1[21], stp1[26],
                           stp1[22], stp1[25])
  }

  /* Stage4 */
  {
    const __m128i lo_9_14 = _mm_unpacklo_epi16(stp2_8, stp2_15);
    const __m128i hi_9_14 = _mm_unpackhi_epi16(stp2_8, stp2_15);
    const __m128i lo_10_13 = _mm_unpacklo_epi16(stp2_11, stp2_12);
    const __m128i hi_10_13 = _mm_unpackhi_epi16(stp2_11, stp2_12);

    stp1[0] = _mm_mulhrs_epi16(in[0], stk4_0);
    stp1[1] = _mm_mulhrs_epi16(in[0], stk4_0);  // stk4_1 = stk4_0
    stp1[2] = stp1[0];
    stp1[3] = stp1[1];

    MULTIPLICATION_AND_ADD(lo_9_14, hi_9_14, lo_10_13, hi_10_13, stg4_4, stg4_5,
                           stg4_6, stg4_4, stp2_9, stp2_14, stp2_10, stp2_13)

    stp2_16 = _mm_add_epi16(stp1[16], stp1[19]);
    stp2_17 = _mm_add_epi16(stp1[17], stp1[18]);
    stp2_18 = _mm_sub_epi16(stp1[17], stp1[18]);
    stp2_19 = _mm_sub_epi16(stp1[16], stp1[19]);
    stp2_20 = _mm_sub_epi16(stp1[23], stp1[20]);
    stp2_21 = _mm_sub_epi16(stp1[22], stp1[21]);
    stp2_22 = _mm_add_epi16(stp1[22], stp1[21]);
    stp2_23 = _mm_add_epi16(stp1[23], stp1[20]);

    stp2_24 = _mm_add_epi16(stp1[24], stp1[27]);
    stp2_25 = _mm_add_epi16(stp1[25], stp1[26]);
    stp2_26 = _mm_sub_epi16(stp1[25], stp1[26]);
    stp2_27 = _mm_sub_epi16(stp1[24], stp1[27]);
    stp2_28 = _mm_sub_epi16(stp1[31], stp1[28]);
    stp2_29 = _mm_sub_epi16(stp1[30], stp1[29]);
    stp2_30 = _mm_add_epi16(stp1[29], stp1[30]);
    stp2_31 = _mm_add_epi16(stp1[28], stp1[31]);
  }

  /* Stage5 */
  {
//  Note:
//  #define AVOID_OVERFLOW = 0, code would be faster.  But it can't pass
//  SingleExtreme test.  The MaxSupportedCoeff/MinSupportedCoeff must drop
//  to 23198 and -23197, respectively.
#define AVOID_OVERFLOW (1)

#if AVOID_OVERFLOW
    const __m128i lo_6_5 = _mm_unpacklo_epi16(stp1[7], stp1[4]);
    const __m128i hi_6_5 = _mm_unpackhi_epi16(stp1[7], stp1[4]);
#endif
    const __m128i lo_18_29 = _mm_unpacklo_epi16(stp2_18, stp2_29);
    const __m128i hi_18_29 = _mm_unpackhi_epi16(stp2_18, stp2_29);

    const __m128i lo_19_28 = _mm_unpacklo_epi16(stp2_19, stp2_28);
    const __m128i hi_19_28 = _mm_unpackhi_epi16(stp2_19, stp2_28);
    const __m128i lo_20_27 = _mm_unpacklo_epi16(stp2_20, stp2_27);
    const __m128i hi_20_27 = _mm_unpackhi_epi16(stp2_20, stp2_27);

    const __m128i lo_21_26 = _mm_unpacklo_epi16(stp2_21, stp2_26);
    const __m128i hi_21_26 = _mm_unpackhi_epi16(stp2_21, stp2_26);

#if AVOID_OVERFLOW
    tmp0 = _mm_madd_epi16(lo_6_5, stg4_1);
    tmp1 = _mm_madd_epi16(hi_6_5, stg4_1);
    tmp2 = _mm_madd_epi16(lo_6_5, stg4_0);
    tmp3 = _mm_madd_epi16(hi_6_5, stg4_0);

    tmp0 = _mm_add_epi32(tmp0, rounding);
    tmp1 = _mm_add_epi32(tmp1, rounding);
    tmp2 = _mm_add_epi32(tmp2, rounding);
    tmp3 = _mm_add_epi32(tmp3, rounding);

    tmp0 = _mm_srai_epi32(tmp0, DCT_CONST_BITS);
    tmp1 = _mm_srai_epi32(tmp1, DCT_CONST_BITS);
    tmp2 = _mm_srai_epi32(tmp2, DCT_CONST_BITS);
    tmp3 = _mm_srai_epi32(tmp3, DCT_CONST_BITS);

    stp1[5] = _mm_packs_epi32(tmp0, tmp1);
    stp1[6] = _mm_packs_epi32(tmp2, tmp3);
#else
    tmp0 = _mm_sub_epi16(stp1[7], stp1[4]);
    tmp1 = _mm_adds_epi16(stp1[7], stp1[4]);
    stp1[5] = _mm_mulhrs_epi16(tmp0, stk4_0);
    stp1[6] = _mm_mulhrs_epi16(tmp1, stk4_0);
#endif

    stp1[8] = _mm_add_epi16(stp2_8, stp2_11);
    stp1[9] = _mm_add_epi16(stp2_9, stp2_10);
    stp1[10] = _mm_sub_epi16(stp2_9, stp2_10);
    stp1[11] = _mm_sub_epi16(stp2_8, stp2_11);
    stp1[12] = _mm_sub_epi16(stp2_15, stp2_12);
    stp1[13] = _mm_sub_epi16(stp2_14, stp2_13);
    stp1[14] = _mm_add_epi16(stp2_14, stp2_13);
    stp1[15] = _mm_add_epi16(stp2_15, stp2_12);

    MULTIPLICATION_AND_ADD(lo_18_29, hi_18_29, lo_19_28, hi_19_28, stg4_4,
                           stg4_5, stg4_4, stg4_5, stp1[18], stp1[29], stp1[19],
                           stp1[28])
    MULTIPLICATION_AND_ADD(lo_20_27, hi_20_27, lo_21_26, hi_21_26, stg4_6,
                           stg4_4, stg4_6, stg4_4, stp1[20], stp1[27], stp1[21],
                           stp1[26])

    stp1[16] = stp2_16;
    stp1[17] = stp2_17;
    stp1[22] = stp2_22;
    stp1[23] = stp2_23;
    stp1[24] = stp2_24;
    stp1[25] = stp2_25;
    stp1[30] = stp2_30;
    stp1[31] = stp2_31;
  }

  /* Stage6 */
  {
#if AVOID_OVERFLOW
    const __m128i lo_10_13 = _mm_unpacklo_epi16(stp1[10], stp1[13]);
    const __m128i hi_10_13 = _mm_unpackhi_epi16(stp1[10], stp1[13]);
    const __m128i lo_11_12 = _mm_unpacklo_epi16(stp1[11], stp1[12]);
    const __m128i hi_11_12 = _mm_unpackhi_epi16(stp1[11], stp1[12]);
#endif

    stp2_0 = _mm_add_epi16(stp1[0], stp1[7]);
    stp2_1 = _mm_add_epi16(stp1[1], stp1[6]);
    stp2_2 = _mm_add_epi16(stp1[2], stp1[5]);
    stp2_3 = _mm_add_epi16(stp1[3], stp1[4]);
    stp2_4 = _mm_sub_epi16(stp1[3], stp1[4]);
    stp2_5 = _mm_sub_epi16(stp1[2], stp1[5]);
    stp2_6 = _mm_sub_epi16(stp1[1], stp1[6]);
    stp2_7 = _mm_sub_epi16(stp1[0], stp1[7]);

    stp2_8 = stp1[8];
    stp2_9 = stp1[9];
    stp2_14 = stp1[14];
    stp2_15 = stp1[15];

#if AVOID_OVERFLOW
    MULTIPLICATION_AND_ADD(lo_10_13, hi_10_13, lo_11_12, hi_11_12, stg6_0,
                           stg4_0, stg6_0, stg4_0, stp2_10, stp2_13, stp2_11,
                           stp2_12)
#else
    tmp0 = _mm_add_epi16(stp1[10], stp1[13]);
    tmp1 = _mm_sub_epi16(stp1[13], stp1[10]);
    tmp2 = _mm_add_epi16(stp1[11], stp1[12]);
    tmp3 = _mm_sub_epi16(stp1[12], stp1[11]);

    stp2_10 = _mm_mulhrs_epi16(tmp1, stk4_0);
    stp2_13 = _mm_mulhrs_epi16(tmp0, stk4_0);
    stp2_11 = _mm_mulhrs_epi16(tmp3, stk4_0);
    stp2_12 = _mm_mulhrs_epi16(tmp2, stk4_0);

#endif

    stp2_16 = _mm_add_epi16(stp1[16], stp1[23]);
    stp2_17 = _mm_add_epi16(stp1[17], stp1[22]);
    stp2_18 = _mm_add_epi16(stp1[18], stp1[21]);
    stp2_19 = _mm_add_epi16(stp1[19], stp1[20]);
    stp2_20 = _mm_sub_epi16(stp1[19], stp1[20]);
    stp2_21 = _mm_sub_epi16(stp1[18], stp1[21]);
    stp2_22 = _mm_sub_epi16(stp1[17], stp1[22]);
    stp2_23 = _mm_sub_epi16(stp1[16], stp1[23]);

    stp2_24 = _mm_sub_epi16(stp1[31], stp1[24]);
    stp2_25 = _mm_sub_epi16(stp1[30], stp1[25]);
    stp2_26 = _mm_sub_epi16(stp1[29], stp1[26]);
    stp2_27 = _mm_sub_epi16(stp1[28], stp1[27]);
    stp2_28 = _mm_add_epi16(stp1[27], stp1[28]);
    stp2_29 = _mm_add_epi16(stp1[26], stp1[29]);
    stp2_30 = _mm_add_epi16(stp1[25], stp1[30]);
    stp2_31 = _mm_add_epi16(stp1[24], stp1[31]);
  }

  /* Stage7 */
  {
#if AVOID_OVERFLOW
    const __m128i lo_20_27 = _mm_unpacklo_epi16(stp2_20, stp2_27);
    const __m128i hi_20_27 = _mm_unpackhi_epi16(stp2_20, stp2_27);
    const __m128i lo_21_26 = _mm_unpacklo_epi16(stp2_21, stp2_26);
    const __m128i hi_21_26 = _mm_unpackhi_epi16(stp2_21, stp2_26);

    const __m128i lo_22_25 = _mm_unpacklo_epi16(stp2_22, stp2_25);
    const __m128i hi_22_25 = _mm_unpackhi_epi16(stp2_22, stp2_25);
    const __m128i lo_23_24 = _mm_unpacklo_epi16(stp2_23, stp2_24);
    const __m128i hi_23_24 = _mm_unpackhi_epi16(stp2_23, stp2_24);
#endif
    stp1[0] = _mm_add_epi16(stp2_0, stp2_15);
    stp1[1] = _mm_add_epi16(stp2_1, stp2_14);
    stp1[2] = _mm_add_epi16(stp2_2, stp2_13);
    stp1[3] = _mm_add_epi16(stp2_3, stp2_12);
    stp1[4] = _mm_add_epi16(stp2_4, stp2_11);
    stp1[5] = _mm_add_epi16(stp2_5, stp2_10);
    stp1[6] = _mm_add_epi16(stp2_6, stp2_9);
    stp1[7] = _mm_add_epi16(stp2_7, stp2_8);
    stp1[8] = _mm_sub_epi16(stp2_7, stp2_8);
    stp1[9] = _mm_sub_epi16(stp2_6, stp2_9);
    stp1[10] = _mm_sub_epi16(stp2_5, stp2_10);
    stp1[11] = _mm_sub_epi16(stp2_4, stp2_11);
    stp1[12] = _mm_sub_epi16(stp2_3, stp2_12);
    stp1[13] = _mm_sub_epi16(stp2_2, stp2_13);
    stp1[14] = _mm_sub_epi16(stp2_1, stp2_14);
    stp1[15] = _mm_sub_epi16(stp2_0, stp2_15);

    stp1[16] = stp2_16;
    stp1[17] = stp2_17;
    stp1[18] = stp2_18;
    stp1[19] = stp2_19;

#if AVOID_OVERFLOW
    MULTIPLICATION_AND_ADD(lo_20_27, hi_20_27, lo_21_26, hi_21_26, stg6_0,
                           stg4_0, stg6_0, stg4_0, stp1[20], stp1[27], stp1[21],
                           stp1[26])
    MULTIPLICATION_AND_ADD(lo_22_25, hi_22_25, lo_23_24, hi_23_24, stg6_0,
                           stg4_0, stg6_0, stg4_0, stp1[22], stp1[25], stp1[23],
                           stp1[24])
#else
    tmp0 = _mm_add_epi16(stp2_20, stp2_27);
    tmp1 = _mm_sub_epi16(stp2_27, stp2_20);
    tmp2 = _mm_add_epi16(stp2_21, stp2_26);
    tmp3 = _mm_sub_epi16(stp2_26, stp2_21);

    stp1[20] = _mm_mulhrs_epi16(tmp1, stk4_0);
    stp1[27] = _mm_mulhrs_epi16(tmp0, stk4_0);
    stp1[21] = _mm_mulhrs_epi16(tmp3, stk4_0);
    stp1[26] = _mm_mulhrs_epi16(tmp2, stk4_0);

    tmp0 = _mm_add_epi16(stp2_22, stp2_25);
    tmp1 = _mm_sub_epi16(stp2_25, stp2_22);
    tmp2 = _mm_add_epi16(stp2_23, stp2_24);
    tmp3 = _mm_sub_epi16(stp2_24, stp2_23);

    stp1[22] = _mm_mulhrs_epi16(tmp1, stk4_0);
    stp1[25] = _mm_mulhrs_epi16(tmp0, stk4_0);
    stp1[23] = _mm_mulhrs_epi16(tmp3, stk4_0);
    stp1[24] = _mm_mulhrs_epi16(tmp2, stk4_0);
#endif

    stp1[28] = stp2_28;
    stp1[29] = stp2_29;
    stp1[30] = stp2_30;
    stp1[31] = stp2_31;
  }
#undef AVOID_OVERFLOW
}

// Only upper-left 8x8 has non-zero coeff
void vpx_idct32x32_34_add_ssse3(const tran_low_t *input, uint8_t *dest,
                                int stride) {
  const __m128i zero = _mm_setzero_si128();
  const __m128i final_rounding = _mm_set1_epi16(1 << 5);
  __m128i in[32], col[32];
  __m128i stp1[32];
  int i;

  // Load input data. Only need to load the top left 8x8 block.
  in[0] = load_input_data(input);
  in[1] = load_input_data(input + 32);
  in[2] = load_input_data(input + 64);
  in[3] = load_input_data(input + 96);
  in[4] = load_input_data(input + 128);
  in[5] = load_input_data(input + 160);
  in[6] = load_input_data(input + 192);
  in[7] = load_input_data(input + 224);

  array_transpose_8x8(in, in);
  idct32_34(in, stp1);

  // 1_D: Store 32 intermediate results for each 8x32 block.
  col[0] = _mm_add_epi16(stp1[0], stp1[31]);
  col[1] = _mm_add_epi16(stp1[1], stp1[30]);
  col[2] = _mm_add_epi16(stp1[2], stp1[29]);
  col[3] = _mm_add_epi16(stp1[3], stp1[28]);
  col[4] = _mm_add_epi16(stp1[4], stp1[27]);
  col[5] = _mm_add_epi16(stp1[5], stp1[26]);
  col[6] = _mm_add_epi16(stp1[6], stp1[25]);
  col[7] = _mm_add_epi16(stp1[7], stp1[24]);
  col[8] = _mm_add_epi16(stp1[8], stp1[23]);
  col[9] = _mm_add_epi16(stp1[9], stp1[22]);
  col[10] = _mm_add_epi16(stp1[10], stp1[21]);
  col[11] = _mm_add_epi16(stp1[11], stp1[20]);
  col[12] = _mm_add_epi16(stp1[12], stp1[19]);
  col[13] = _mm_add_epi16(stp1[13], stp1[18]);
  col[14] = _mm_add_epi16(stp1[14], stp1[17]);
  col[15] = _mm_add_epi16(stp1[15], stp1[16]);
  col[16] = _mm_sub_epi16(stp1[15], stp1[16]);
  col[17] = _mm_sub_epi16(stp1[14], stp1[17]);
  col[18] = _mm_sub_epi16(stp1[13], stp1[18]);
  col[19] = _mm_sub_epi16(stp1[12], stp1[19]);
  col[20] = _mm_sub_epi16(stp1[11], stp1[20]);
  col[21] = _mm_sub_epi16(stp1[10], stp1[21]);
  col[22] = _mm_sub_epi16(stp1[9], stp1[22]);
  col[23] = _mm_sub_epi16(stp1[8], stp1[23]);
  col[24] = _mm_sub_epi16(stp1[7], stp1[24]);
  col[25] = _mm_sub_epi16(stp1[6], stp1[25]);
  col[26] = _mm_sub_epi16(stp1[5], stp1[26]);
  col[27] = _mm_sub_epi16(stp1[4], stp1[27]);
  col[28] = _mm_sub_epi16(stp1[3], stp1[28]);
  col[29] = _mm_sub_epi16(stp1[2], stp1[29]);
  col[30] = _mm_sub_epi16(stp1[1], stp1[30]);
  col[31] = _mm_sub_epi16(stp1[0], stp1[31]);
  for (i = 0; i < 4; i++) {
    int j;
    // Transpose 32x8 block to 8x32 block
    array_transpose_8x8(col + i * 8, in);
    idct32_34(in, stp1);

    // 2_D: Calculate the results and store them to destination.
    in[0] = _mm_add_epi16(stp1[0], stp1[31]);
    in[1] = _mm_add_epi16(stp1[1], stp1[30]);
    in[2] = _mm_add_epi16(stp1[2], stp1[29]);
    in[3] = _mm_add_epi16(stp1[3], stp1[28]);
    in[4] = _mm_add_epi16(stp1[4], stp1[27]);
    in[5] = _mm_add_epi16(stp1[5], stp1[26]);
    in[6] = _mm_add_epi16(stp1[6], stp1[25]);
    in[7] = _mm_add_epi16(stp1[7], stp1[24]);
    in[8] = _mm_add_epi16(stp1[8], stp1[23]);
    in[9] = _mm_add_epi16(stp1[9], stp1[22]);
    in[10] = _mm_add_epi16(stp1[10], stp1[21]);
    in[11] = _mm_add_epi16(stp1[11], stp1[20]);
    in[12] = _mm_add_epi16(stp1[12], stp1[19]);
    in[13] = _mm_add_epi16(stp1[13], stp1[18]);
    in[14] = _mm_add_epi16(stp1[14], stp1[17]);
    in[15] = _mm_add_epi16(stp1[15], stp1[16]);
    in[16] = _mm_sub_epi16(stp1[15], stp1[16]);
    in[17] = _mm_sub_epi16(stp1[14], stp1[17]);
    in[18] = _mm_sub_epi16(stp1[13], stp1[18]);
    in[19] = _mm_sub_epi16(stp1[12], stp1[19]);
    in[20] = _mm_sub_epi16(stp1[11], stp1[20]);
    in[21] = _mm_sub_epi16(stp1[10], stp1[21]);
    in[22] = _mm_sub_epi16(stp1[9], stp1[22]);
    in[23] = _mm_sub_epi16(stp1[8], stp1[23]);
    in[24] = _mm_sub_epi16(stp1[7], stp1[24]);
    in[25] = _mm_sub_epi16(stp1[6], stp1[25]);
    in[26] = _mm_sub_epi16(stp1[5], stp1[26]);
    in[27] = _mm_sub_epi16(stp1[4], stp1[27]);
    in[28] = _mm_sub_epi16(stp1[3], stp1[28]);
    in[29] = _mm_sub_epi16(stp1[2], stp1[29]);
    in[30] = _mm_sub_epi16(stp1[1], stp1[30]);
    in[31] = _mm_sub_epi16(stp1[0], stp1[31]);

    for (j = 0; j < 32; ++j) {
      // Final rounding and shift
      in[j] = _mm_adds_epi16(in[j], final_rounding);
      in[j] = _mm_srai_epi16(in[j], 6);
      RECON_AND_STORE(dest + j * stride, in[j]);
    }

    dest += 8;
  }
}

// in0[16] represents the left 8x16 block
// in1[16] represents the right 8x16 block
static void load_buffer_16x16(const tran_low_t *input, __m128i *in0,
                              __m128i *in1) {
  int i;
  for (i = 0; i < 16; i++) {
    in0[i] = load_input_data(input);
    in1[i] = load_input_data(input + 8);
    input += 32;
  }
}

static void array_transpose_16x16_2(__m128i *in0, __m128i *in1, __m128i *out0,
                                    __m128i *out1) {
  array_transpose_8x8(in0, out0);
  array_transpose_8x8(&in0[8], out1);
  array_transpose_8x8(in1, &out0[8]);
  array_transpose_8x8(&in1[8], &out1[8]);
}

// For each 8x16 block __m128i in[16], output __m128i col[32]
static void idct32_8x16_135(__m128i *in) {
  const __m128i rounding = _mm_set1_epi32(DCT_CONST_ROUNDING);
  const __m128i stk1_0 = pair_set_epi16(2 * cospi_31_64, 2 * cospi_31_64);
  const __m128i stk1_1 = pair_set_epi16(2 * cospi_1_64, 2 * cospi_1_64);
  const __m128i stk1_2 = pair_set_epi16(-2 * cospi_17_64, -2 * cospi_17_64);
  const __m128i stk1_3 = pair_set_epi16(2 * cospi_15_64, 2 * cospi_15_64);

  const __m128i stk1_4 = pair_set_epi16(2 * cospi_23_64, 2 * cospi_23_64);
  const __m128i stk1_5 = pair_set_epi16(2 * cospi_9_64, 2 * cospi_9_64);
  const __m128i stk1_6 = pair_set_epi16(-2 * cospi_25_64, -2 * cospi_25_64);
  const __m128i stk1_7 = pair_set_epi16(2 * cospi_7_64, 2 * cospi_7_64);

  const __m128i stk1_8 = pair_set_epi16(2 * cospi_27_64, 2 * cospi_27_64);
  const __m128i stk1_9 = pair_set_epi16(2 * cospi_5_64, 2 * cospi_5_64);
  const __m128i stk1_10 = pair_set_epi16(-2 * cospi_21_64, -2 * cospi_21_64);
  const __m128i stk1_11 = pair_set_epi16(2 * cospi_11_64, 2 * cospi_11_64);

  const __m128i stk1_12 = pair_set_epi16(2 * cospi_19_64, 2 * cospi_19_64);
  const __m128i stk1_13 = pair_set_epi16(2 * cospi_13_64, 2 * cospi_13_64);
  const __m128i stk1_14 = pair_set_epi16(-2 * cospi_29_64, -2 * cospi_29_64);
  const __m128i stk1_15 = pair_set_epi16(2 * cospi_3_64, 2 * cospi_3_64);

  const __m128i stk2_0 = pair_set_epi16(2 * cospi_30_64, 2 * cospi_30_64);
  const __m128i stk2_1 = pair_set_epi16(2 * cospi_2_64, 2 * cospi_2_64);
  const __m128i stk2_2 = pair_set_epi16(-2 * cospi_18_64, -2 * cospi_18_64);
  const __m128i stk2_3 = pair_set_epi16(2 * cospi_14_64, 2 * cospi_14_64);

  const __m128i stk2_4 = pair_set_epi16(2 * cospi_22_64, 2 * cospi_22_64);
  const __m128i stk2_5 = pair_set_epi16(2 * cospi_10_64, 2 * cospi_10_64);
  const __m128i stk2_6 = pair_set_epi16(-2 * cospi_26_64, -2 * cospi_26_64);
  const __m128i stk2_7 = pair_set_epi16(2 * cospi_6_64, 2 * cospi_6_64);

  const __m128i stk3_0 = pair_set_epi16(2 * cospi_28_64, 2 * cospi_28_64);
  const __m128i stk3_1 = pair_set_epi16(2 * cospi_4_64, 2 * cospi_4_64);
  const __m128i stk3_2 = pair_set_epi16(-2 * cospi_20_64, -2 * cospi_20_64);
  const __m128i stk3_3 = pair_set_epi16(2 * cospi_12_64, 2 * cospi_12_64);

  const __m128i stg3_4 = pair_set_epi16(-cospi_4_64, cospi_28_64);
  const __m128i stg3_5 = pair_set_epi16(cospi_28_64, cospi_4_64);
  const __m128i stg3_6 = pair_set_epi16(-cospi_28_64, -cospi_4_64);
  const __m128i stg3_8 = pair_set_epi16(-cospi_20_64, cospi_12_64);
  const __m128i stg3_9 = pair_set_epi16(cospi_12_64, cospi_20_64);
  const __m128i stg3_10 = pair_set_epi16(-cospi_12_64, -cospi_20_64);

  const __m128i stg4_0 = pair_set_epi16(cospi_16_64, cospi_16_64);
  const __m128i stg4_1 = pair_set_epi16(cospi_16_64, -cospi_16_64);
  const __m128i stk4_0 = pair_set_epi16(2 * cospi_16_64, 2 * cospi_16_64);
  const __m128i stk4_2 = pair_set_epi16(2 * cospi_24_64, 2 * cospi_24_64);
  const __m128i stk4_3 = pair_set_epi16(2 * cospi_8_64, 2 * cospi_8_64);

  const __m128i stg4_4 = pair_set_epi16(-cospi_8_64, cospi_24_64);
  const __m128i stg4_5 = pair_set_epi16(cospi_24_64, cospi_8_64);
  const __m128i stg4_6 = pair_set_epi16(-cospi_24_64, -cospi_8_64);

  const __m128i stg6_0 = pair_set_epi16(-cospi_16_64, cospi_16_64);
  __m128i stp1_0, stp1_1, stp1_2, stp1_3, stp1_4, stp1_5, stp1_6, stp1_7,
      stp1_8, stp1_9, stp1_10, stp1_11, stp1_12, stp1_13, stp1_14, stp1_15,
      stp1_16, stp1_17, stp1_18, stp1_19, stp1_20, stp1_21, stp1_22, stp1_23,
      stp1_24, stp1_25, stp1_26, stp1_27, stp1_28, stp1_29, stp1_30, stp1_31;
  __m128i stp2_0, stp2_1, stp2_2, stp2_3, stp2_4, stp2_5, stp2_6, stp2_7,
      stp2_8, stp2_9, stp2_10, stp2_11, stp2_12, stp2_13, stp2_14, stp2_15,
      stp2_16, stp2_17, stp2_18, stp2_19, stp2_20, stp2_21, stp2_22, stp2_23,
      stp2_24, stp2_25, stp2_26, stp2_27, stp2_28, stp2_29, stp2_30, stp2_31;
  __m128i tmp0, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7;

  /* Stage1 */
  stp1_16 = _mm_mulhrs_epi16(in[1], stk1_0);
  stp1_31 = _mm_mulhrs_epi16(in[1], stk1_1);
  stp1_17 = _mm_mulhrs_epi16(in[15], stk1_2);
  stp1_30 = _mm_mulhrs_epi16(in[15], stk1_3);

  stp1_18 = _mm_mulhrs_epi16(in[9], stk1_4);
  stp1_29 = _mm_mulhrs_epi16(in[9], stk1_5);
  stp1_19 = _mm_mulhrs_epi16(in[7], stk1_6);
  stp1_28 = _mm_mulhrs_epi16(in[7], stk1_7);

  stp1_20 = _mm_mulhrs_epi16(in[5], stk1_8);
  stp1_27 = _mm_mulhrs_epi16(in[5], stk1_9);
  stp1_21 = _mm_mulhrs_epi16(in[11], stk1_10);
  stp1_26 = _mm_mulhrs_epi16(in[11], stk1_11);

  stp1_22 = _mm_mulhrs_epi16(in[13], stk1_12);
  stp1_25 = _mm_mulhrs_epi16(in[13], stk1_13);
  stp1_23 = _mm_mulhrs_epi16(in[3], stk1_14);
  stp1_24 = _mm_mulhrs_epi16(in[3], stk1_15);

  /* Stage2 */
  stp2_8 = _mm_mulhrs_epi16(in[2], stk2_0);
  stp2_15 = _mm_mulhrs_epi16(in[2], stk2_1);
  stp2_9 = _mm_mulhrs_epi16(in[14], stk2_2);
  stp2_14 = _mm_mulhrs_epi16(in[14], stk2_3);

  stp2_10 = _mm_mulhrs_epi16(in[10], stk2_4);
  stp2_13 = _mm_mulhrs_epi16(in[10], stk2_5);
  stp2_11 = _mm_mulhrs_epi16(in[6], stk2_6);
  stp2_12 = _mm_mulhrs_epi16(in[6], stk2_7);

  stp2_16 = _mm_add_epi16(stp1_16, stp1_17);
  stp2_17 = _mm_sub_epi16(stp1_16, stp1_17);
  stp2_18 = _mm_sub_epi16(stp1_19, stp1_18);
  stp2_19 = _mm_add_epi16(stp1_19, stp1_18);

  stp2_20 = _mm_add_epi16(stp1_20, stp1_21);
  stp2_21 = _mm_sub_epi16(stp1_20, stp1_21);
  stp2_22 = _mm_sub_epi16(stp1_23, stp1_22);
  stp2_23 = _mm_add_epi16(stp1_23, stp1_22);

  stp2_24 = _mm_add_epi16(stp1_24, stp1_25);
  stp2_25 = _mm_sub_epi16(stp1_24, stp1_25);
  stp2_26 = _mm_sub_epi16(stp1_27, stp1_26);
  stp2_27 = _mm_add_epi16(stp1_27, stp1_26);

  stp2_28 = _mm_add_epi16(stp1_28, stp1_29);
  stp2_29 = _mm_sub_epi16(stp1_28, stp1_29);
  stp2_30 = _mm_sub_epi16(stp1_31, stp1_30);
  stp2_31 = _mm_add_epi16(stp1_31, stp1_30);

  /* Stage3 */
  {
    const __m128i lo_17_30 = _mm_unpacklo_epi16(stp2_17, stp2_30);
    const __m128i hi_17_30 = _mm_unpackhi_epi16(stp2_17, stp2_30);
    const __m128i lo_18_29 = _mm_unpacklo_epi16(stp2_18, stp2_29);
    const __m128i hi_18_29 = _mm_unpackhi_epi16(stp2_18, stp2_29);

    const __m128i lo_21_26 = _mm_unpacklo_epi16(stp2_21, stp2_26);
    const __m128i hi_21_26 = _mm_unpackhi_epi16(stp2_21, stp2_26);
    const __m128i lo_22_25 = _mm_unpacklo_epi16(stp2_22, stp2_25);
    const __m128i hi_22_25 = _mm_unpackhi_epi16(stp2_22, stp2_25);

    stp1_4 = _mm_mulhrs_epi16(in[4], stk3_0);
    stp1_7 = _mm_mulhrs_epi16(in[4], stk3_1);
    stp1_5 = _mm_mulhrs_epi16(in[12], stk3_2);
    stp1_6 = _mm_mulhrs_epi16(in[12], stk3_3);

    stp2_0 = _mm_mulhrs_epi16(in[0], stk4_0);
    stp2_1 = _mm_mulhrs_epi16(in[0], stk4_0);  // stk4_1 = stk4_0
    stp2_2 = _mm_mulhrs_epi16(in[8], stk4_2);
    stp2_3 = _mm_mulhrs_epi16(in[8], stk4_3);

    stp1_8 = _mm_add_epi16(stp2_8, stp2_9);
    stp1_9 = _mm_sub_epi16(stp2_8, stp2_9);
    stp1_10 = _mm_sub_epi16(stp2_11, stp2_10);
    stp1_11 = _mm_add_epi16(stp2_11, stp2_10);
    stp1_12 = _mm_add_epi16(stp2_12, stp2_13);
    stp1_13 = _mm_sub_epi16(stp2_12, stp2_13);
    stp1_14 = _mm_sub_epi16(stp2_15, stp2_14);
    stp1_15 = _mm_add_epi16(stp2_15, stp2_14);

    MULTIPLICATION_AND_ADD(lo_17_30, hi_17_30, lo_18_29, hi_18_29, stg3_4,
                           stg3_5, stg3_6, stg3_4, stp1_17, stp1_30, stp1_18,
                           stp1_29)
    MULTIPLICATION_AND_ADD(lo_21_26, hi_21_26, lo_22_25, hi_22_25, stg3_8,
                           stg3_9, stg3_10, stg3_8, stp1_21, stp1_26, stp1_22,
                           stp1_25)

    stp1_16 = stp2_16;
    stp1_31 = stp2_31;
    stp1_19 = stp2_19;
    stp1_20 = stp2_20;
    stp1_23 = stp2_23;
    stp1_24 = stp2_24;
    stp1_27 = stp2_27;
    stp1_28 = stp2_28;
  }

  /* Stage4 */
  {
    const __m128i lo_9_14 = _mm_unpacklo_epi16(stp1_9, stp1_14);
    const __m128i hi_9_14 = _mm_unpackhi_epi16(stp1_9, stp1_14);
    const __m128i lo_10_13 = _mm_unpacklo_epi16(stp1_10, stp1_13);
    const __m128i hi_10_13 = _mm_unpackhi_epi16(stp1_10, stp1_13);

    stp2_4 = _mm_add_epi16(stp1_4, stp1_5);
    stp2_5 = _mm_sub_epi16(stp1_4, stp1_5);
    stp2_6 = _mm_sub_epi16(stp1_7, stp1_6);
    stp2_7 = _mm_add_epi16(stp1_7, stp1_6);

    MULTIPLICATION_AND_ADD(lo_9_14, hi_9_14, lo_10_13, hi_10_13, stg4_4, stg4_5,
                           stg4_6, stg4_4, stp2_9, stp2_14, stp2_10, stp2_13)

    stp2_8 = stp1_8;
    stp2_15 = stp1_15;
    stp2_11 = stp1_11;
    stp2_12 = stp1_12;

    stp2_16 = _mm_add_epi16(stp1_16, stp1_19);
    stp2_17 = _mm_add_epi16(stp1_17, stp1_18);
    stp2_18 = _mm_sub_epi16(stp1_17, stp1_18);
    stp2_19 = _mm_sub_epi16(stp1_16, stp1_19);
    stp2_20 = _mm_sub_epi16(stp1_23, stp1_20);
    stp2_21 = _mm_sub_epi16(stp1_22, stp1_21);
    stp2_22 = _mm_add_epi16(stp1_22, stp1_21);
    stp2_23 = _mm_add_epi16(stp1_23, stp1_20);

    stp2_24 = _mm_add_epi16(stp1_24, stp1_27);
    stp2_25 = _mm_add_epi16(stp1_25, stp1_26);
    stp2_26 = _mm_sub_epi16(stp1_25, stp1_26);
    stp2_27 = _mm_sub_epi16(stp1_24, stp1_27);
    stp2_28 = _mm_sub_epi16(stp1_31, stp1_28);
    stp2_29 = _mm_sub_epi16(stp1_30, stp1_29);
    stp2_30 = _mm_add_epi16(stp1_29, stp1_30);
    stp2_31 = _mm_add_epi16(stp1_28, stp1_31);
  }

  /* Stage5 */
  {
    const __m128i lo_6_5 = _mm_unpacklo_epi16(stp2_6, stp2_5);
    const __m128i hi_6_5 = _mm_unpackhi_epi16(stp2_6, stp2_5);
    const __m128i lo_18_29 = _mm_unpacklo_epi16(stp2_18, stp2_29);
    const __m128i hi_18_29 = _mm_unpackhi_epi16(stp2_18, stp2_29);

    const __m128i lo_19_28 = _mm_unpacklo_epi16(stp2_19, stp2_28);
    const __m128i hi_19_28 = _mm_unpackhi_epi16(stp2_19, stp2_28);
    const __m128i lo_20_27 = _mm_unpacklo_epi16(stp2_20, stp2_27);
    const __m128i hi_20_27 = _mm_unpackhi_epi16(stp2_20, stp2_27);

    const __m128i lo_21_26 = _mm_unpacklo_epi16(stp2_21, stp2_26);
    const __m128i hi_21_26 = _mm_unpackhi_epi16(stp2_21, stp2_26);

    stp1_0 = _mm_add_epi16(stp2_0, stp2_3);
    stp1_1 = _mm_add_epi16(stp2_1, stp2_2);
    stp1_2 = _mm_sub_epi16(stp2_1, stp2_2);
    stp1_3 = _mm_sub_epi16(stp2_0, stp2_3);

    tmp0 = _mm_madd_epi16(lo_6_5, stg4_1);
    tmp1 = _mm_madd_epi16(hi_6_5, stg4_1);
    tmp2 = _mm_madd_epi16(lo_6_5, stg4_0);
    tmp3 = _mm_madd_epi16(hi_6_5, stg4_0);

    tmp0 = _mm_add_epi32(tmp0, rounding);
    tmp1 = _mm_add_epi32(tmp1, rounding);
    tmp2 = _mm_add_epi32(tmp2, rounding);
    tmp3 = _mm_add_epi32(tmp3, rounding);

    tmp0 = _mm_srai_epi32(tmp0, DCT_CONST_BITS);
    tmp1 = _mm_srai_epi32(tmp1, DCT_CONST_BITS);
    tmp2 = _mm_srai_epi32(tmp2, DCT_CONST_BITS);
    tmp3 = _mm_srai_epi32(tmp3, DCT_CONST_BITS);

    stp1_5 = _mm_packs_epi32(tmp0, tmp1);
    stp1_6 = _mm_packs_epi32(tmp2, tmp3);

    stp1_4 = stp2_4;
    stp1_7 = stp2_7;

    stp1_8 = _mm_add_epi16(stp2_8, stp2_11);
    stp1_9 = _mm_add_epi16(stp2_9, stp2_10);
    stp1_10 = _mm_sub_epi16(stp2_9, stp2_10);
    stp1_11 = _mm_sub_epi16(stp2_8, stp2_11);
    stp1_12 = _mm_sub_epi16(stp2_15, stp2_12);
    stp1_13 = _mm_sub_epi16(stp2_14, stp2_13);
    stp1_14 = _mm_add_epi16(stp2_14, stp2_13);
    stp1_15 = _mm_add_epi16(stp2_15, stp2_12);

    stp1_16 = stp2_16;
    stp1_17 = stp2_17;

    MULTIPLICATION_AND_ADD(lo_18_29, hi_18_29, lo_19_28, hi_19_28, stg4_4,
                           stg4_5, stg4_4, stg4_5, stp1_18, stp1_29, stp1_19,
                           stp1_28)
    MULTIPLICATION_AND_ADD(lo_20_27, hi_20_27, lo_21_26, hi_21_26, stg4_6,
                           stg4_4, stg4_6, stg4_4, stp1_20, stp1_27, stp1_21,
                           stp1_26)

    stp1_22 = stp2_22;
    stp1_23 = stp2_23;
    stp1_24 = stp2_24;
    stp1_25 = stp2_25;
    stp1_30 = stp2_30;
    stp1_31 = stp2_31;
  }

  /* Stage6 */
  {
    const __m128i lo_10_13 = _mm_unpacklo_epi16(stp1_10, stp1_13);
    const __m128i hi_10_13 = _mm_unpackhi_epi16(stp1_10, stp1_13);
    const __m128i lo_11_12 = _mm_unpacklo_epi16(stp1_11, stp1_12);
    const __m128i hi_11_12 = _mm_unpackhi_epi16(stp1_11, stp1_12);

    stp2_0 = _mm_add_epi16(stp1_0, stp1_7);
    stp2_1 = _mm_add_epi16(stp1_1, stp1_6);
    stp2_2 = _mm_add_epi16(stp1_2, stp1_5);
    stp2_3 = _mm_add_epi16(stp1_3, stp1_4);
    stp2_4 = _mm_sub_epi16(stp1_3, stp1_4);
    stp2_5 = _mm_sub_epi16(stp1_2, stp1_5);
    stp2_6 = _mm_sub_epi16(stp1_1, stp1_6);
    stp2_7 = _mm_sub_epi16(stp1_0, stp1_7);

    stp2_8 = stp1_8;
    stp2_9 = stp1_9;
    stp2_14 = stp1_14;
    stp2_15 = stp1_15;

    MULTIPLICATION_AND_ADD(lo_10_13, hi_10_13, lo_11_12, hi_11_12, stg6_0,
                           stg4_0, stg6_0, stg4_0, stp2_10, stp2_13, stp2_11,
                           stp2_12)

    stp2_16 = _mm_add_epi16(stp1_16, stp1_23);
    stp2_17 = _mm_add_epi16(stp1_17, stp1_22);
    stp2_18 = _mm_add_epi16(stp1_18, stp1_21);
    stp2_19 = _mm_add_epi16(stp1_19, stp1_20);
    stp2_20 = _mm_sub_epi16(stp1_19, stp1_20);
    stp2_21 = _mm_sub_epi16(stp1_18, stp1_21);
    stp2_22 = _mm_sub_epi16(stp1_17, stp1_22);
    stp2_23 = _mm_sub_epi16(stp1_16, stp1_23);

    stp2_24 = _mm_sub_epi16(stp1_31, stp1_24);
    stp2_25 = _mm_sub_epi16(stp1_30, stp1_25);
    stp2_26 = _mm_sub_epi16(stp1_29, stp1_26);
    stp2_27 = _mm_sub_epi16(stp1_28, stp1_27);
    stp2_28 = _mm_add_epi16(stp1_27, stp1_28);
    stp2_29 = _mm_add_epi16(stp1_26, stp1_29);
    stp2_30 = _mm_add_epi16(stp1_25, stp1_30);
    stp2_31 = _mm_add_epi16(stp1_24, stp1_31);
  }

  /* Stage7 */
  {
    const __m128i lo_20_27 = _mm_unpacklo_epi16(stp2_20, stp2_27);
    const __m128i hi_20_27 = _mm_unpackhi_epi16(stp2_20, stp2_27);
    const __m128i lo_21_26 = _mm_unpacklo_epi16(stp2_21, stp2_26);
    const __m128i hi_21_26 = _mm_unpackhi_epi16(stp2_21, stp2_26);

    const __m128i lo_22_25 = _mm_unpacklo_epi16(stp2_22, stp2_25);
    const __m128i hi_22_25 = _mm_unpackhi_epi16(stp2_22, stp2_25);
    const __m128i lo_23_24 = _mm_unpacklo_epi16(stp2_23, stp2_24);
    const __m128i hi_23_24 = _mm_unpackhi_epi16(stp2_23, stp2_24);

    stp1_0 = _mm_add_epi16(stp2_0, stp2_15);
    stp1_1 = _mm_add_epi16(stp2_1, stp2_14);
    stp1_2 = _mm_add_epi16(stp2_2, stp2_13);
    stp1_3 = _mm_add_epi16(stp2_3, stp2_12);
    stp1_4 = _mm_add_epi16(stp2_4, stp2_11);
    stp1_5 = _mm_add_epi16(stp2_5, stp2_10);
    stp1_6 = _mm_add_epi16(stp2_6, stp2_9);
    stp1_7 = _mm_add_epi16(stp2_7, stp2_8);
    stp1_8 = _mm_sub_epi16(stp2_7, stp2_8);
    stp1_9 = _mm_sub_epi16(stp2_6, stp2_9);
    stp1_10 = _mm_sub_epi16(stp2_5, stp2_10);
    stp1_11 = _mm_sub_epi16(stp2_4, stp2_11);
    stp1_12 = _mm_sub_epi16(stp2_3, stp2_12);
    stp1_13 = _mm_sub_epi16(stp2_2, stp2_13);
    stp1_14 = _mm_sub_epi16(stp2_1, stp2_14);
    stp1_15 = _mm_sub_epi16(stp2_0, stp2_15);

    stp1_16 = stp2_16;
    stp1_17 = stp2_17;
    stp1_18 = stp2_18;
    stp1_19 = stp2_19;

    MULTIPLICATION_AND_ADD(lo_20_27, hi_20_27, lo_21_26, hi_21_26, stg6_0,
                           stg4_0, stg6_0, stg4_0, stp1_20, stp1_27, stp1_21,
                           stp1_26)
    MULTIPLICATION_AND_ADD(lo_22_25, hi_22_25, lo_23_24, hi_23_24, stg6_0,
                           stg4_0, stg6_0, stg4_0, stp1_22, stp1_25, stp1_23,
                           stp1_24)

    stp1_28 = stp2_28;
    stp1_29 = stp2_29;
    stp1_30 = stp2_30;
    stp1_31 = stp2_31;
  }

  in[0] = _mm_add_epi16(stp1_0, stp1_31);
  in[1] = _mm_add_epi16(stp1_1, stp1_30);
  in[2] = _mm_add_epi16(stp1_2, stp1_29);
  in[3] = _mm_add_epi16(stp1_3, stp1_28);
  in[4] = _mm_add_epi16(stp1_4, stp1_27);
  in[5] = _mm_add_epi16(stp1_5, stp1_26);
  in[6] = _mm_add_epi16(stp1_6, stp1_25);
  in[7] = _mm_add_epi16(stp1_7, stp1_24);
  in[8] = _mm_add_epi16(stp1_8, stp1_23);
  in[9] = _mm_add_epi16(stp1_9, stp1_22);
  in[10] = _mm_add_epi16(stp1_10, stp1_21);
  in[11] = _mm_add_epi16(stp1_11, stp1_20);
  in[12] = _mm_add_epi16(stp1_12, stp1_19);
  in[13] = _mm_add_epi16(stp1_13, stp1_18);
  in[14] = _mm_add_epi16(stp1_14, stp1_17);
  in[15] = _mm_add_epi16(stp1_15, stp1_16);
  in[16] = _mm_sub_epi16(stp1_15, stp1_16);
  in[17] = _mm_sub_epi16(stp1_14, stp1_17);
  in[18] = _mm_sub_epi16(stp1_13, stp1_18);
  in[19] = _mm_sub_epi16(stp1_12, stp1_19);
  in[20] = _mm_sub_epi16(stp1_11, stp1_20);
  in[21] = _mm_sub_epi16(stp1_10, stp1_21);
  in[22] = _mm_sub_epi16(stp1_9, stp1_22);
  in[23] = _mm_sub_epi16(stp1_8, stp1_23);
  in[24] = _mm_sub_epi16(stp1_7, stp1_24);
  in[25] = _mm_sub_epi16(stp1_6, stp1_25);
  in[26] = _mm_sub_epi16(stp1_5, stp1_26);
  in[27] = _mm_sub_epi16(stp1_4, stp1_27);
  in[28] = _mm_sub_epi16(stp1_3, stp1_28);
  in[29] = _mm_sub_epi16(stp1_2, stp1_29);
  in[30] = _mm_sub_epi16(stp1_1, stp1_30);
  in[31] = _mm_sub_epi16(stp1_0, stp1_31);
}

static INLINE void store_buffer_8x32(__m128i *in, uint8_t *dst, int stride) {
  const __m128i final_rounding = _mm_set1_epi16(1 << 5);
  const __m128i zero = _mm_setzero_si128();
  int j = 0;
  while (j < 32) {
    in[j] = _mm_adds_epi16(in[j], final_rounding);
    in[j + 1] = _mm_adds_epi16(in[j + 1], final_rounding);

    in[j] = _mm_srai_epi16(in[j], 6);
    in[j + 1] = _mm_srai_epi16(in[j + 1], 6);

    RECON_AND_STORE(dst, in[j]);
    dst += stride;
    RECON_AND_STORE(dst, in[j + 1]);
    dst += stride;
    j += 2;
  }
}

static INLINE void recon_and_store(__m128i *in0, __m128i *in1, uint8_t *dest,
                                   int stride) {
  store_buffer_8x32(in0, dest, stride);
  store_buffer_8x32(in1, dest + 8, stride);
}

static INLINE void idct32_135(__m128i *col0, __m128i *col1) {
  idct32_8x16_135(col0);
  idct32_8x16_135(col1);
}

// Only upper-left 16x16 has non-zero coeff
void vpx_idct32x32_135_add_ssse3(const tran_low_t *input, uint8_t *dest,
                                 int stride) {
  __m128i col0[32], col1[32], col2[32], col3[32];

  // Load input data. Only need to load the top left 16x16 block.
  load_buffer_16x16(input, col2, col3);

  // columns
  array_transpose_16x16_2(col2, col3, col0, col1);
  idct32_135(col0, col1);

  // rows
  array_transpose_16x16_2(col0, col1, col2, col3);
  idct32_135(col2, col3);
  recon_and_store(col2, col3, dest, stride);

  array_transpose_16x16_2(&col0[16], &col1[16], col2, col3);
  idct32_135(col2, col3);
  recon_and_store(col2, col3, dest + 16, stride);
}
