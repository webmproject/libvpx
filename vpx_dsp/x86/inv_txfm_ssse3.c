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
#include "vpx_dsp/x86/transpose_sse2.h"
#include "vpx_dsp/x86/txfm_common_sse2.h"

void vpx_idct8x8_12_add_ssse3(const tran_low_t *input, uint8_t *dest,
                              int stride) {
  const __m128i zero = _mm_setzero_si128();
  const __m128i rounding = _mm_set1_epi32(DCT_CONST_ROUNDING);
  const __m128i stg1_0 = pair_set_epi16(2 * cospi_28_64, 2 * cospi_28_64);
  const __m128i stg1_1 = pair_set_epi16(2 * cospi_4_64, 2 * cospi_4_64);
  const __m128i stg1_2 = pair_set_epi16(-2 * cospi_20_64, -2 * cospi_20_64);
  const __m128i stg1_3 = pair_set_epi16(2 * cospi_12_64, 2 * cospi_12_64);
  const __m128i stg2_0 = pair_set_epi16(2 * cospi_16_64, 2 * cospi_16_64);
  const __m128i stk2_0 = pair_set_epi16(cospi_16_64, cospi_16_64);
  const __m128i stk2_1 = pair_set_epi16(cospi_16_64, -cospi_16_64);
  const __m128i stg2_2 = pair_set_epi16(2 * cospi_24_64, 2 * cospi_24_64);
  const __m128i stg2_3 = pair_set_epi16(2 * cospi_8_64, 2 * cospi_8_64);
  const __m128i stg3_0 = pair_set_epi16(-cospi_16_64, cospi_16_64);

  __m128i in[8];
  __m128i stp1_0, stp1_1, stp1_2, stp1_3, stp1_4, stp1_5, stp1_6, stp1_7;
  __m128i stp2_0, stp2_1, stp2_2, stp2_3, stp2_4, stp2_5, stp2_6, stp2_7;
  __m128i tmp[4];

  // Rows. Load 4-row input data.
  in[0] = load_input_data(input);
  in[1] = load_input_data(input + 8 * 1);
  in[2] = load_input_data(input + 8 * 2);
  in[3] = load_input_data(input + 8 * 3);

  // 4x4 Transpose
  transpose_16bit_4x4(in, in);

  // Stage1
  tmp[0] = _mm_mulhrs_epi16(in[0], stg1_0);
  tmp[1] = _mm_mulhrs_epi16(in[0], stg1_1);
  tmp[2] = _mm_mulhrs_epi16(in[1], stg1_2);
  tmp[3] = _mm_mulhrs_epi16(in[1], stg1_3);

  stp1_4 = _mm_unpackhi_epi64(tmp[0], tmp[1]);
  stp1_5 = _mm_unpackhi_epi64(tmp[2], tmp[3]);

  // Stage2
  tmp[0] = _mm_mulhrs_epi16(in[0], stg2_0);
  stp2_0 = _mm_unpacklo_epi64(tmp[0], tmp[0]);

  tmp[1] = _mm_mulhrs_epi16(in[1], stg2_2);
  tmp[2] = _mm_mulhrs_epi16(in[1], stg2_3);
  stp2_2 = _mm_unpacklo_epi64(tmp[2], tmp[1]);

  tmp[0] = _mm_add_epi16(stp1_4, stp1_5);
  tmp[1] = _mm_sub_epi16(stp1_4, stp1_5);

  stp2_4 = tmp[0];
  stp2_5 = _mm_unpacklo_epi64(tmp[1], zero);
  stp2_6 = _mm_unpackhi_epi64(tmp[1], zero);

  tmp[0] = _mm_unpacklo_epi16(stp2_5, stp2_6);
  tmp[1] = _mm_madd_epi16(tmp[0], stg3_0);
  tmp[2] = _mm_madd_epi16(tmp[0], stk2_0);  // stg3_1 = stk2_0

  tmp[1] = _mm_add_epi32(tmp[1], rounding);
  tmp[2] = _mm_add_epi32(tmp[2], rounding);
  tmp[1] = _mm_srai_epi32(tmp[1], DCT_CONST_BITS);
  tmp[2] = _mm_srai_epi32(tmp[2], DCT_CONST_BITS);

  stp1_5 = _mm_packs_epi32(tmp[1], tmp[2]);

  // Stage3
  tmp[2] = _mm_add_epi16(stp2_0, stp2_2);
  tmp[3] = _mm_sub_epi16(stp2_0, stp2_2);

  stp1_2 = _mm_unpackhi_epi64(tmp[3], tmp[2]);
  stp1_3 = _mm_unpacklo_epi64(tmp[3], tmp[2]);

  // Stage4
  tmp[0] = _mm_add_epi16(stp1_3, stp2_4);
  tmp[1] = _mm_add_epi16(stp1_2, stp1_5);
  tmp[2] = _mm_sub_epi16(stp1_3, stp2_4);
  tmp[3] = _mm_sub_epi16(stp1_2, stp1_5);

  idct8x8_12_transpose_16bit_4x8(tmp, in);

  /* Stage1 */
  stp1_4 = _mm_mulhrs_epi16(in[1], stg1_0);
  stp1_7 = _mm_mulhrs_epi16(in[1], stg1_1);
  stp1_5 = _mm_mulhrs_epi16(in[3], stg1_2);
  stp1_6 = _mm_mulhrs_epi16(in[3], stg1_3);

  /* Stage2 */
  stp2_0 = _mm_mulhrs_epi16(in[0], stg2_0);
  stp2_1 = _mm_mulhrs_epi16(in[0], stg2_0);

  stp2_2 = _mm_mulhrs_epi16(in[2], stg2_2);
  stp2_3 = _mm_mulhrs_epi16(in[2], stg2_3);

  stp2_4 = _mm_add_epi16(stp1_4, stp1_5);
  stp2_5 = _mm_sub_epi16(stp1_4, stp1_5);
  stp2_6 = _mm_sub_epi16(stp1_7, stp1_6);
  stp2_7 = _mm_add_epi16(stp1_7, stp1_6);

  /* Stage3 */
  stp1_0 = _mm_add_epi16(stp2_0, stp2_3);
  stp1_1 = _mm_add_epi16(stp2_1, stp2_2);
  stp1_2 = _mm_sub_epi16(stp2_1, stp2_2);
  stp1_3 = _mm_sub_epi16(stp2_0, stp2_3);

  tmp[0] = _mm_unpacklo_epi16(stp2_6, stp2_5);
  tmp[1] = _mm_unpackhi_epi16(stp2_6, stp2_5);

  tmp[2] = _mm_madd_epi16(tmp[0], stk2_0);
  tmp[3] = _mm_madd_epi16(tmp[1], stk2_0);
  tmp[2] = _mm_add_epi32(tmp[2], rounding);
  tmp[3] = _mm_add_epi32(tmp[3], rounding);
  tmp[2] = _mm_srai_epi32(tmp[2], DCT_CONST_BITS);
  tmp[3] = _mm_srai_epi32(tmp[3], DCT_CONST_BITS);
  stp1_6 = _mm_packs_epi32(tmp[2], tmp[3]);

  tmp[2] = _mm_madd_epi16(tmp[0], stk2_1);
  tmp[3] = _mm_madd_epi16(tmp[1], stk2_1);
  tmp[2] = _mm_add_epi32(tmp[2], rounding);
  tmp[3] = _mm_add_epi32(tmp[3], rounding);
  tmp[2] = _mm_srai_epi32(tmp[2], DCT_CONST_BITS);
  tmp[3] = _mm_srai_epi32(tmp[3], DCT_CONST_BITS);
  stp1_5 = _mm_packs_epi32(tmp[2], tmp[3]);

  /* Stage4  */
  in[0] = _mm_add_epi16(stp1_0, stp2_7);
  in[1] = _mm_add_epi16(stp1_1, stp1_6);
  in[2] = _mm_add_epi16(stp1_2, stp1_5);
  in[3] = _mm_add_epi16(stp1_3, stp2_4);
  in[4] = _mm_sub_epi16(stp1_3, stp2_4);
  in[5] = _mm_sub_epi16(stp1_2, stp1_5);
  in[6] = _mm_sub_epi16(stp1_1, stp1_6);
  in[7] = _mm_sub_epi16(stp1_0, stp2_7);

  write_buffer_8x8(in, dest, stride);
}

static void idct32_34_first_half(const __m128i *in, __m128i *stp1) {
  const __m128i stk2_0 = pair_set_epi16(2 * cospi_30_64, 2 * cospi_30_64);
  const __m128i stk2_1 = pair_set_epi16(2 * cospi_2_64, 2 * cospi_2_64);
  const __m128i stk2_6 = pair_set_epi16(-2 * cospi_26_64, -2 * cospi_26_64);
  const __m128i stk2_7 = pair_set_epi16(2 * cospi_6_64, 2 * cospi_6_64);

  const __m128i stk3_0 = pair_set_epi16(2 * cospi_28_64, 2 * cospi_28_64);
  const __m128i stk3_1 = pair_set_epi16(2 * cospi_4_64, 2 * cospi_4_64);

  const __m128i stg4_0 = pair_set_epi16(cospi_16_64, cospi_16_64);
  const __m128i stk4_0 = pair_set_epi16(2 * cospi_16_64, 2 * cospi_16_64);
  const __m128i stg4_1 = pair_set_epi16(cospi_16_64, -cospi_16_64);
  const __m128i stg4_4 = pair_set_epi16(-cospi_8_64, cospi_24_64);
  const __m128i stg4_5 = pair_set_epi16(cospi_24_64, cospi_8_64);
  const __m128i stg4_6 = pair_set_epi16(-cospi_24_64, -cospi_8_64);

  const __m128i stg6_0 = pair_set_epi16(-cospi_16_64, cospi_16_64);
  __m128i u0, u1, u2, u3, u4, u5, u6, u7;
  __m128i x0, x1, x4, x5, x6, x7;
  __m128i v0, v1, v2, v3, v4, v5, v6, v7, v8, v9, v10, v11, v12, v13, v14, v15;

  // phase 1

  // 0, 15
  u2 = _mm_mulhrs_epi16(in[2], stk2_1);  // stp2_15
  u3 = _mm_mulhrs_epi16(in[6], stk2_7);  // stp2_12
  v15 = _mm_add_epi16(u2, u3);
  // in[0], in[4]
  x0 = _mm_mulhrs_epi16(in[0], stk4_0);  // stp1[0]
  x7 = _mm_mulhrs_epi16(in[4], stk3_1);  // stp1[7]
  v0 = _mm_add_epi16(x0, x7);            // stp2_0
  stp1[0] = _mm_add_epi16(v0, v15);
  stp1[15] = _mm_sub_epi16(v0, v15);

  // in[2], in[6]
  u0 = _mm_mulhrs_epi16(in[2], stk2_0);             // stp2_8
  u1 = _mm_mulhrs_epi16(in[6], stk2_6);             // stp2_11
  butterfly(&u0, &u2, &stg4_4, &stg4_5, &u4, &u5);  // stp2_9, stp2_14
  butterfly(&u1, &u3, &stg4_6, &stg4_4, &u6, &u7);  // stp2_10, stp2_13

  v8 = _mm_add_epi16(u0, u1);
  v9 = _mm_add_epi16(u4, u6);
  v10 = _mm_sub_epi16(u4, u6);
  v11 = _mm_sub_epi16(u0, u1);
  v12 = _mm_sub_epi16(u2, u3);
  v13 = _mm_sub_epi16(u5, u7);
  v14 = _mm_add_epi16(u5, u7);

  butterfly_self(&v10, &v13, &stg6_0, &stg4_0);
  butterfly_self(&v11, &v12, &stg6_0, &stg4_0);

  // 1, 14
  x1 = _mm_mulhrs_epi16(in[0], stk4_0);  // stp1[1], stk4_1 = stk4_0
  // stp1[2] = stp1[0], stp1[3] = stp1[1]
  x4 = _mm_mulhrs_epi16(in[4], stk3_0);  // stp1[4]
  butterfly(&x7, &x4, &stg4_1, &stg4_0, &x5, &x6);
  v1 = _mm_add_epi16(x1, x6);  // stp2_1
  v2 = _mm_add_epi16(x0, x5);  // stp2_2
  stp1[1] = _mm_add_epi16(v1, v14);
  stp1[14] = _mm_sub_epi16(v1, v14);

  stp1[2] = _mm_add_epi16(v2, v13);
  stp1[13] = _mm_sub_epi16(v2, v13);

  v3 = _mm_add_epi16(x1, x4);  // stp2_3
  v4 = _mm_sub_epi16(x1, x4);  // stp2_4

  v5 = _mm_sub_epi16(x0, x5);  // stp2_5

  v6 = _mm_sub_epi16(x1, x6);  // stp2_6
  v7 = _mm_sub_epi16(x0, x7);  // stp2_7
  stp1[3] = _mm_add_epi16(v3, v12);
  stp1[12] = _mm_sub_epi16(v3, v12);

  stp1[6] = _mm_add_epi16(v6, v9);
  stp1[9] = _mm_sub_epi16(v6, v9);

  stp1[7] = _mm_add_epi16(v7, v8);
  stp1[8] = _mm_sub_epi16(v7, v8);

  stp1[4] = _mm_add_epi16(v4, v11);
  stp1[11] = _mm_sub_epi16(v4, v11);

  stp1[5] = _mm_add_epi16(v5, v10);
  stp1[10] = _mm_sub_epi16(v5, v10);
}

static void idct32_34_second_half(const __m128i *in, __m128i *stp1) {
  const __m128i stk1_0 = pair_set_epi16(2 * cospi_31_64, 2 * cospi_31_64);
  const __m128i stk1_1 = pair_set_epi16(2 * cospi_1_64, 2 * cospi_1_64);
  const __m128i stk1_6 = pair_set_epi16(-2 * cospi_25_64, -2 * cospi_25_64);
  const __m128i stk1_7 = pair_set_epi16(2 * cospi_7_64, 2 * cospi_7_64);
  const __m128i stk1_8 = pair_set_epi16(2 * cospi_27_64, 2 * cospi_27_64);
  const __m128i stk1_9 = pair_set_epi16(2 * cospi_5_64, 2 * cospi_5_64);
  const __m128i stk1_14 = pair_set_epi16(-2 * cospi_29_64, -2 * cospi_29_64);
  const __m128i stk1_15 = pair_set_epi16(2 * cospi_3_64, 2 * cospi_3_64);
  const __m128i stg3_4 = pair_set_epi16(-cospi_4_64, cospi_28_64);
  const __m128i stg3_5 = pair_set_epi16(cospi_28_64, cospi_4_64);
  const __m128i stg3_6 = pair_set_epi16(-cospi_28_64, -cospi_4_64);
  const __m128i stg3_8 = pair_set_epi16(-cospi_20_64, cospi_12_64);
  const __m128i stg3_9 = pair_set_epi16(cospi_12_64, cospi_20_64);
  const __m128i stg3_10 = pair_set_epi16(-cospi_12_64, -cospi_20_64);

  const __m128i stg4_0 = pair_set_epi16(cospi_16_64, cospi_16_64);
  const __m128i stg4_4 = pair_set_epi16(-cospi_8_64, cospi_24_64);
  const __m128i stg4_5 = pair_set_epi16(cospi_24_64, cospi_8_64);
  const __m128i stg4_6 = pair_set_epi16(-cospi_24_64, -cospi_8_64);

  const __m128i stg6_0 = pair_set_epi16(-cospi_16_64, cospi_16_64);
  __m128i v16, v17, v18, v19, v20, v21, v22, v23;
  __m128i v24, v25, v26, v27, v28, v29, v30, v31;
  __m128i u16, u17, u18, u19, u20, u21, u22, u23;
  __m128i u24, u25, u26, u27, u28, u29, u30, u31;

  v16 = _mm_mulhrs_epi16(in[1], stk1_0);
  v31 = _mm_mulhrs_epi16(in[1], stk1_1);

  v19 = _mm_mulhrs_epi16(in[7], stk1_6);
  v28 = _mm_mulhrs_epi16(in[7], stk1_7);

  v20 = _mm_mulhrs_epi16(in[5], stk1_8);
  v27 = _mm_mulhrs_epi16(in[5], stk1_9);

  v23 = _mm_mulhrs_epi16(in[3], stk1_14);
  v24 = _mm_mulhrs_epi16(in[3], stk1_15);

  butterfly(&v16, &v31, &stg3_4, &stg3_5, &v17, &v30);
  butterfly(&v19, &v28, &stg3_6, &stg3_4, &v18, &v29);
  butterfly(&v20, &v27, &stg3_8, &stg3_9, &v21, &v26);
  butterfly(&v23, &v24, &stg3_10, &stg3_8, &v22, &v25);

  u16 = _mm_add_epi16(v16, v19);
  u17 = _mm_add_epi16(v17, v18);
  u18 = _mm_sub_epi16(v17, v18);
  u19 = _mm_sub_epi16(v16, v19);
  u20 = _mm_sub_epi16(v23, v20);
  u21 = _mm_sub_epi16(v22, v21);
  u22 = _mm_add_epi16(v22, v21);
  u23 = _mm_add_epi16(v23, v20);
  u24 = _mm_add_epi16(v24, v27);
  u27 = _mm_sub_epi16(v24, v27);
  u25 = _mm_add_epi16(v25, v26);
  u26 = _mm_sub_epi16(v25, v26);
  u28 = _mm_sub_epi16(v31, v28);
  u31 = _mm_add_epi16(v28, v31);
  u29 = _mm_sub_epi16(v30, v29);
  u30 = _mm_add_epi16(v29, v30);

  butterfly_self(&u18, &u29, &stg4_4, &stg4_5);
  butterfly_self(&u19, &u28, &stg4_4, &stg4_5);
  butterfly_self(&u20, &u27, &stg4_6, &stg4_4);
  butterfly_self(&u21, &u26, &stg4_6, &stg4_4);

  stp1[16] = _mm_add_epi16(u16, u23);
  stp1[23] = _mm_sub_epi16(u16, u23);

  stp1[17] = _mm_add_epi16(u17, u22);
  stp1[22] = _mm_sub_epi16(u17, u22);

  stp1[18] = _mm_add_epi16(u18, u21);
  stp1[21] = _mm_sub_epi16(u18, u21);

  stp1[19] = _mm_add_epi16(u19, u20);
  stp1[20] = _mm_sub_epi16(u19, u20);

  stp1[24] = _mm_sub_epi16(u31, u24);
  stp1[31] = _mm_add_epi16(u24, u31);

  stp1[25] = _mm_sub_epi16(u30, u25);
  stp1[30] = _mm_add_epi16(u25, u30);

  stp1[26] = _mm_sub_epi16(u29, u26);
  stp1[29] = _mm_add_epi16(u26, u29);

  stp1[27] = _mm_sub_epi16(u28, u27);
  stp1[28] = _mm_add_epi16(u27, u28);

  butterfly_self(&stp1[20], &stp1[27], &stg6_0, &stg4_0);
  butterfly_self(&stp1[21], &stp1[26], &stg6_0, &stg4_0);
  butterfly_self(&stp1[22], &stp1[25], &stg6_0, &stg4_0);
  butterfly_self(&stp1[23], &stp1[24], &stg6_0, &stg4_0);
}

// Only upper-left 8x8 has non-zero coeff
void vpx_idct32x32_34_add_ssse3(const tran_low_t *input, uint8_t *dest,
                                int stride) {
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

  transpose_16bit_8x8(in, in);
  idct32_34_first_half(in, stp1);
  idct32_34_second_half(in, stp1);

  // 1_D: Store 32 intermediate results for each 8x32 block.
  add_sub_butterfly(stp1, col, 32);
  for (i = 0; i < 4; i++) {
    int j;
    // Transpose 32x8 block to 8x32 block
    transpose_16bit_8x8(col + i * 8, in);
    idct32_34_first_half(in, stp1);
    idct32_34_second_half(in, stp1);

    // 2_D: Calculate the results and store them to destination.
    add_sub_butterfly(stp1, in, 32);
    for (j = 0; j < 32; ++j) {
      // Final rounding and shift
      in[j] = _mm_adds_epi16(in[j], final_rounding);
      in[j] = _mm_srai_epi16(in[j], 6);
      recon_and_store(dest + j * stride, in[j]);
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

// Group the coefficient calculation into smaller functions
// to prevent stack spillover:
// quarter_1: 0-7
// quarter_2: 8-15
// quarter_3_4: 16-23, 24-31
static void idct32_8x32_135_quarter_1(const __m128i *in /*in[16]*/,
                                      __m128i *out /*out[8]*/) {
  __m128i u0, u1, u2, u3, u4, u5, u6, u7;
  __m128i v0, v1, v2, v3, v4, v5, v6, v7;

  {
    const __m128i stk4_0 = pair_set_epi16(2 * cospi_16_64, 2 * cospi_16_64);
    const __m128i stk4_2 = pair_set_epi16(2 * cospi_24_64, 2 * cospi_24_64);
    const __m128i stk4_3 = pair_set_epi16(2 * cospi_8_64, 2 * cospi_8_64);
    u0 = _mm_mulhrs_epi16(in[0], stk4_0);
    u2 = _mm_mulhrs_epi16(in[8], stk4_2);
    u3 = _mm_mulhrs_epi16(in[8], stk4_3);
    u1 = u0;
  }

  v0 = _mm_add_epi16(u0, u3);
  v1 = _mm_add_epi16(u1, u2);
  v2 = _mm_sub_epi16(u1, u2);
  v3 = _mm_sub_epi16(u0, u3);

  {
    const __m128i stk3_0 = pair_set_epi16(2 * cospi_28_64, 2 * cospi_28_64);
    const __m128i stk3_1 = pair_set_epi16(2 * cospi_4_64, 2 * cospi_4_64);
    const __m128i stk3_2 = pair_set_epi16(-2 * cospi_20_64, -2 * cospi_20_64);
    const __m128i stk3_3 = pair_set_epi16(2 * cospi_12_64, 2 * cospi_12_64);
    u4 = _mm_mulhrs_epi16(in[4], stk3_0);
    u7 = _mm_mulhrs_epi16(in[4], stk3_1);
    u5 = _mm_mulhrs_epi16(in[12], stk3_2);
    u6 = _mm_mulhrs_epi16(in[12], stk3_3);
  }

  v4 = _mm_add_epi16(u4, u5);
  v5 = _mm_sub_epi16(u4, u5);
  v6 = _mm_sub_epi16(u7, u6);
  v7 = _mm_add_epi16(u7, u6);

  {
    const __m128i stg4_0 = pair_set_epi16(cospi_16_64, cospi_16_64);
    const __m128i stg4_1 = pair_set_epi16(cospi_16_64, -cospi_16_64);
    butterfly(&v6, &v5, &stg4_1, &stg4_0, &v5, &v6);
  }

  out[0] = _mm_add_epi16(v0, v7);
  out[1] = _mm_add_epi16(v1, v6);
  out[2] = _mm_add_epi16(v2, v5);
  out[3] = _mm_add_epi16(v3, v4);
  out[4] = _mm_sub_epi16(v3, v4);
  out[5] = _mm_sub_epi16(v2, v5);
  out[6] = _mm_sub_epi16(v1, v6);
  out[7] = _mm_sub_epi16(v0, v7);
}

static void idct32_8x32_135_quarter_2(const __m128i *in /*in[16]*/,
                                      __m128i *out /*out[8]*/) {
  __m128i u8, u9, u10, u11, u12, u13, u14, u15;
  __m128i v8, v9, v10, v11, v12, v13, v14, v15;

  {
    const __m128i stk2_0 = pair_set_epi16(2 * cospi_30_64, 2 * cospi_30_64);
    const __m128i stk2_1 = pair_set_epi16(2 * cospi_2_64, 2 * cospi_2_64);
    const __m128i stk2_2 = pair_set_epi16(-2 * cospi_18_64, -2 * cospi_18_64);
    const __m128i stk2_3 = pair_set_epi16(2 * cospi_14_64, 2 * cospi_14_64);
    const __m128i stk2_4 = pair_set_epi16(2 * cospi_22_64, 2 * cospi_22_64);
    const __m128i stk2_5 = pair_set_epi16(2 * cospi_10_64, 2 * cospi_10_64);
    const __m128i stk2_6 = pair_set_epi16(-2 * cospi_26_64, -2 * cospi_26_64);
    const __m128i stk2_7 = pair_set_epi16(2 * cospi_6_64, 2 * cospi_6_64);
    u8 = _mm_mulhrs_epi16(in[2], stk2_0);
    u15 = _mm_mulhrs_epi16(in[2], stk2_1);
    u9 = _mm_mulhrs_epi16(in[14], stk2_2);
    u14 = _mm_mulhrs_epi16(in[14], stk2_3);
    u10 = _mm_mulhrs_epi16(in[10], stk2_4);
    u13 = _mm_mulhrs_epi16(in[10], stk2_5);
    u11 = _mm_mulhrs_epi16(in[6], stk2_6);
    u12 = _mm_mulhrs_epi16(in[6], stk2_7);
  }

  v8 = _mm_add_epi16(u8, u9);
  v9 = _mm_sub_epi16(u8, u9);
  v10 = _mm_sub_epi16(u11, u10);
  v11 = _mm_add_epi16(u11, u10);
  v12 = _mm_add_epi16(u12, u13);
  v13 = _mm_sub_epi16(u12, u13);
  v14 = _mm_sub_epi16(u15, u14);
  v15 = _mm_add_epi16(u15, u14);

  {
    const __m128i stg4_4 = pair_set_epi16(-cospi_8_64, cospi_24_64);
    const __m128i stg4_5 = pair_set_epi16(cospi_24_64, cospi_8_64);
    const __m128i stg4_6 = pair_set_epi16(-cospi_24_64, -cospi_8_64);
    butterfly_self(&v9, &v14, &stg4_4, &stg4_5);
    butterfly_self(&v10, &v13, &stg4_6, &stg4_4);
  }

  out[0] = _mm_add_epi16(v8, v11);
  out[1] = _mm_add_epi16(v9, v10);
  out[2] = _mm_sub_epi16(v9, v10);
  out[3] = _mm_sub_epi16(v8, v11);
  out[4] = _mm_sub_epi16(v15, v12);
  out[5] = _mm_sub_epi16(v14, v13);
  out[6] = _mm_add_epi16(v14, v13);
  out[7] = _mm_add_epi16(v15, v12);

  {
    const __m128i stg4_0 = pair_set_epi16(cospi_16_64, cospi_16_64);
    const __m128i stg6_0 = pair_set_epi16(-cospi_16_64, cospi_16_64);
    butterfly_self(&out[2], &out[5], &stg6_0, &stg4_0);
    butterfly_self(&out[3], &out[4], &stg6_0, &stg4_0);
  }
}

// 8x32 block even indexed 8 inputs of in[16],
// output first half 16 to out[32]
static void idct32_8x32_quarter_1_2(const __m128i *in /*in[16]*/,
                                    __m128i *out /*out[32]*/) {
  __m128i temp[16];
  idct32_8x32_135_quarter_1(in, temp);
  idct32_8x32_135_quarter_2(in, &temp[8]);
  add_sub_butterfly(temp, out, 16);
}

// 8x32 block odd indexed 8 inputs of in[16],
// output second half 16 to out[32]
static void idct32_8x32_quarter_3_4(const __m128i *in /*in[16]*/,
                                    __m128i *out /*out[32]*/) {
  __m128i v16, v17, v18, v19, v20, v21, v22, v23;
  __m128i v24, v25, v26, v27, v28, v29, v30, v31;
  __m128i u16, u17, u18, u19, u20, u21, u22, u23;
  __m128i u24, u25, u26, u27, u28, u29, u30, u31;

  {
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
    u16 = _mm_mulhrs_epi16(in[1], stk1_0);
    u31 = _mm_mulhrs_epi16(in[1], stk1_1);
    u17 = _mm_mulhrs_epi16(in[15], stk1_2);
    u30 = _mm_mulhrs_epi16(in[15], stk1_3);

    u18 = _mm_mulhrs_epi16(in[9], stk1_4);
    u29 = _mm_mulhrs_epi16(in[9], stk1_5);
    u19 = _mm_mulhrs_epi16(in[7], stk1_6);
    u28 = _mm_mulhrs_epi16(in[7], stk1_7);

    u20 = _mm_mulhrs_epi16(in[5], stk1_8);
    u27 = _mm_mulhrs_epi16(in[5], stk1_9);
    u21 = _mm_mulhrs_epi16(in[11], stk1_10);
    u26 = _mm_mulhrs_epi16(in[11], stk1_11);

    u22 = _mm_mulhrs_epi16(in[13], stk1_12);
    u25 = _mm_mulhrs_epi16(in[13], stk1_13);
    u23 = _mm_mulhrs_epi16(in[3], stk1_14);
    u24 = _mm_mulhrs_epi16(in[3], stk1_15);
  }

  v16 = _mm_add_epi16(u16, u17);
  v17 = _mm_sub_epi16(u16, u17);
  v18 = _mm_sub_epi16(u19, u18);
  v19 = _mm_add_epi16(u19, u18);

  v20 = _mm_add_epi16(u20, u21);
  v21 = _mm_sub_epi16(u20, u21);
  v22 = _mm_sub_epi16(u23, u22);
  v23 = _mm_add_epi16(u23, u22);

  v24 = _mm_add_epi16(u24, u25);
  v25 = _mm_sub_epi16(u24, u25);
  v26 = _mm_sub_epi16(u27, u26);
  v27 = _mm_add_epi16(u27, u26);

  v28 = _mm_add_epi16(u28, u29);
  v29 = _mm_sub_epi16(u28, u29);
  v30 = _mm_sub_epi16(u31, u30);
  v31 = _mm_add_epi16(u31, u30);

  {
    const __m128i stg3_4 = pair_set_epi16(-cospi_4_64, cospi_28_64);
    const __m128i stg3_5 = pair_set_epi16(cospi_28_64, cospi_4_64);
    const __m128i stg3_6 = pair_set_epi16(-cospi_28_64, -cospi_4_64);
    const __m128i stg3_8 = pair_set_epi16(-cospi_20_64, cospi_12_64);
    const __m128i stg3_9 = pair_set_epi16(cospi_12_64, cospi_20_64);
    const __m128i stg3_10 = pair_set_epi16(-cospi_12_64, -cospi_20_64);

    butterfly_self(&v17, &v30, &stg3_4, &stg3_5);
    butterfly_self(&v18, &v29, &stg3_6, &stg3_4);
    butterfly_self(&v21, &v26, &stg3_8, &stg3_9);
    butterfly_self(&v22, &v25, &stg3_10, &stg3_8);
  }

  u16 = _mm_add_epi16(v16, v19);
  u17 = _mm_add_epi16(v17, v18);
  u18 = _mm_sub_epi16(v17, v18);
  u19 = _mm_sub_epi16(v16, v19);
  u20 = _mm_sub_epi16(v23, v20);
  u21 = _mm_sub_epi16(v22, v21);
  u22 = _mm_add_epi16(v22, v21);
  u23 = _mm_add_epi16(v23, v20);

  u24 = _mm_add_epi16(v24, v27);
  u25 = _mm_add_epi16(v25, v26);
  u26 = _mm_sub_epi16(v25, v26);
  u27 = _mm_sub_epi16(v24, v27);
  u28 = _mm_sub_epi16(v31, v28);
  u29 = _mm_sub_epi16(v30, v29);
  u30 = _mm_add_epi16(v29, v30);
  u31 = _mm_add_epi16(v28, v31);

  {
    const __m128i stg4_4 = pair_set_epi16(-cospi_8_64, cospi_24_64);
    const __m128i stg4_5 = pair_set_epi16(cospi_24_64, cospi_8_64);
    const __m128i stg4_6 = pair_set_epi16(-cospi_24_64, -cospi_8_64);
    butterfly_self(&u18, &u29, &stg4_4, &stg4_5);
    butterfly_self(&u19, &u28, &stg4_4, &stg4_5);
    butterfly_self(&u20, &u27, &stg4_6, &stg4_4);
    butterfly_self(&u21, &u26, &stg4_6, &stg4_4);
  }

  out[0] = _mm_add_epi16(u16, u23);
  out[1] = _mm_add_epi16(u17, u22);
  out[2] = _mm_add_epi16(u18, u21);
  out[3] = _mm_add_epi16(u19, u20);
  v20 = _mm_sub_epi16(u19, u20);
  v21 = _mm_sub_epi16(u18, u21);
  v22 = _mm_sub_epi16(u17, u22);
  v23 = _mm_sub_epi16(u16, u23);

  v24 = _mm_sub_epi16(u31, u24);
  v25 = _mm_sub_epi16(u30, u25);
  v26 = _mm_sub_epi16(u29, u26);
  v27 = _mm_sub_epi16(u28, u27);
  out[12] = _mm_add_epi16(u27, u28);
  out[13] = _mm_add_epi16(u26, u29);
  out[14] = _mm_add_epi16(u25, u30);
  out[15] = _mm_add_epi16(u24, u31);

  {
    const __m128i stg4_0 = pair_set_epi16(cospi_16_64, cospi_16_64);
    const __m128i stg6_0 = pair_set_epi16(-cospi_16_64, cospi_16_64);
    butterfly(&v20, &v27, &stg6_0, &stg4_0, &out[4], &out[11]);
    butterfly(&v21, &v26, &stg6_0, &stg4_0, &out[5], &out[10]);
    butterfly(&v22, &v25, &stg6_0, &stg4_0, &out[6], &out[9]);
    butterfly(&v23, &v24, &stg6_0, &stg4_0, &out[7], &out[8]);
  }
}

// 8x16 block, input __m128i in[16], output __m128i in[32]
static void idct32_8x32_135(__m128i *in /*in[32]*/) {
  __m128i out[32];
  idct32_8x32_quarter_1_2(in, out);
  idct32_8x32_quarter_3_4(in, &out[16]);
  add_sub_butterfly(out, in, 32);
}

static INLINE void recon_and_store_ssse3(__m128i *in0, __m128i *in1,
                                         uint8_t *dest, int stride) {
  store_buffer_8x32(in0, dest, stride);
  store_buffer_8x32(in1, dest + 8, stride);
}

static INLINE void idct32_135(__m128i *col0, __m128i *col1) {
  idct32_8x32_135(col0);
  idct32_8x32_135(col1);
}

typedef enum { left_16, right_16 } ColsIndicator;

static void transpose_and_copy_16x16(__m128i *in0, __m128i *in1, __m128i *store,
                                     ColsIndicator cols) {
  switch (cols) {
    case left_16: {
      int i;
      transpose_16bit_16x16(in0, in1);
      for (i = 0; i < 16; ++i) {
        store[i] = in0[16 + i];
        store[16 + i] = in1[16 + i];
      }
      break;
    }
    case right_16: {
      transpose_16bit_8x8(store, in0);
      transpose_16bit_8x8(&store[8], in1);
      transpose_16bit_8x8(&store[16], &in0[8]);
      transpose_16bit_8x8(&store[24], &in1[8]);
      break;
    }
    default: { assert(0); }
  }
}

// Only upper-left 16x16 has non-zero coeff
void vpx_idct32x32_135_add_ssse3(const tran_low_t *input, uint8_t *dest,
                                 int stride) {
  // Each array represents an 8x32 block
  __m128i col0[32], col1[32];
  // This array represents a 16x16 block
  __m128i temp[32];

  // Load input data. Only need to load the top left 16x16 block.
  load_buffer_16x16(input, col0, col1);

  // columns
  transpose_16bit_16x16(col0, col1);
  idct32_135(col0, col1);

  // rows
  transpose_and_copy_16x16(col0, col1, temp, left_16);
  idct32_135(col0, col1);
  recon_and_store_ssse3(col0, col1, dest, stride);

  transpose_and_copy_16x16(col0, col1, temp, right_16);
  idct32_135(col0, col1);
  recon_and_store_ssse3(col0, col1, dest + 16, stride);
}
