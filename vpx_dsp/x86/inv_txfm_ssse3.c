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
