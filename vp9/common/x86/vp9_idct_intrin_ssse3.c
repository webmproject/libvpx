/*
 *  Copyright (c) 2014 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include <tmmintrin.h>  // SSSE3
#include "vp9/common/x86/vp9_idct_intrin_sse2.h"

static void idct16_8col(__m128i *in) {
  const __m128i k__cospi_p30_m02 = pair_set_epi16(cospi_30_64, -cospi_2_64);
  const __m128i k__cospi_p02_p30 = pair_set_epi16(cospi_2_64, cospi_30_64);
  const __m128i k__cospi_p14_m18 = pair_set_epi16(cospi_14_64, -cospi_18_64);
  const __m128i k__cospi_p18_p14 = pair_set_epi16(cospi_18_64, cospi_14_64);
  const __m128i k__cospi_p22_m10 = pair_set_epi16(cospi_22_64, -cospi_10_64);
  const __m128i k__cospi_p10_p22 = pair_set_epi16(cospi_10_64, cospi_22_64);
  const __m128i k__cospi_p06_m26 = pair_set_epi16(cospi_6_64, -cospi_26_64);
  const __m128i k__cospi_p26_p06 = pair_set_epi16(cospi_26_64, cospi_6_64);
  const __m128i k__cospi_p28_m04 = pair_set_epi16(cospi_28_64, -cospi_4_64);
  const __m128i k__cospi_p04_p28 = pair_set_epi16(cospi_4_64, cospi_28_64);
  const __m128i k__cospi_p12_m20 = pair_set_epi16(cospi_12_64, -cospi_20_64);
  const __m128i k__cospi_p20_p12 = pair_set_epi16(cospi_20_64, cospi_12_64);
  const __m128i k__cospi_p24_m08 = pair_set_epi16(cospi_24_64, -cospi_8_64);
  const __m128i k__cospi_p08_p24 = pair_set_epi16(cospi_8_64, cospi_24_64);
  const __m128i k__cospi_m08_p24 = pair_set_epi16(-cospi_8_64, cospi_24_64);
  const __m128i k__cospi_p24_p08 = pair_set_epi16(cospi_24_64, cospi_8_64);
  const __m128i k__cospi_m24_m08 = pair_set_epi16(-cospi_24_64, -cospi_8_64);
  const __m128i k__DCT_CONST_ROUNDING = _mm_set1_epi32(DCT_CONST_ROUNDING);
  const __m128i k__cospi_p16_p16_x2 = pair_set_epi16(23170, 23170);

  __m128i v[16], u[16], s[16], t[16];

  // stage 1
  s[0] = in[0];
  s[1] = in[8];
  s[2] = in[4];
  s[3] = in[12];
  s[4] = in[2];
  s[5] = in[10];
  s[6] = in[6];
  s[7] = in[14];
  s[8] = in[1];
  s[9] = in[9];
  s[10] = in[5];
  s[11] = in[13];
  s[12] = in[3];
  s[13] = in[11];
  s[14] = in[7];
  s[15] = in[15];

  // stage 2
  u[0] = _mm_unpacklo_epi16(s[8], s[15]);
  u[1] = _mm_unpackhi_epi16(s[8], s[15]);
  u[2] = _mm_unpacklo_epi16(s[9], s[14]);
  u[3] = _mm_unpackhi_epi16(s[9], s[14]);
  u[4] = _mm_unpacklo_epi16(s[10], s[13]);
  u[5] = _mm_unpackhi_epi16(s[10], s[13]);
  u[6] = _mm_unpacklo_epi16(s[11], s[12]);
  u[7] = _mm_unpackhi_epi16(s[11], s[12]);

  v[0] = _mm_madd_epi16(u[0], k__cospi_p30_m02);
  v[1] = _mm_madd_epi16(u[1], k__cospi_p30_m02);
  v[2] = _mm_madd_epi16(u[0], k__cospi_p02_p30);
  v[3] = _mm_madd_epi16(u[1], k__cospi_p02_p30);
  v[4] = _mm_madd_epi16(u[2], k__cospi_p14_m18);
  v[5] = _mm_madd_epi16(u[3], k__cospi_p14_m18);
  v[6] = _mm_madd_epi16(u[2], k__cospi_p18_p14);
  v[7] = _mm_madd_epi16(u[3], k__cospi_p18_p14);
  v[8] = _mm_madd_epi16(u[4], k__cospi_p22_m10);
  v[9] = _mm_madd_epi16(u[5], k__cospi_p22_m10);
  v[10] = _mm_madd_epi16(u[4], k__cospi_p10_p22);
  v[11] = _mm_madd_epi16(u[5], k__cospi_p10_p22);
  v[12] = _mm_madd_epi16(u[6], k__cospi_p06_m26);
  v[13] = _mm_madd_epi16(u[7], k__cospi_p06_m26);
  v[14] = _mm_madd_epi16(u[6], k__cospi_p26_p06);
  v[15] = _mm_madd_epi16(u[7], k__cospi_p26_p06);

  u[0] = _mm_add_epi32(v[0], k__DCT_CONST_ROUNDING);
  u[1] = _mm_add_epi32(v[1], k__DCT_CONST_ROUNDING);
  u[2] = _mm_add_epi32(v[2], k__DCT_CONST_ROUNDING);
  u[3] = _mm_add_epi32(v[3], k__DCT_CONST_ROUNDING);
  u[4] = _mm_add_epi32(v[4], k__DCT_CONST_ROUNDING);
  u[5] = _mm_add_epi32(v[5], k__DCT_CONST_ROUNDING);
  u[6] = _mm_add_epi32(v[6], k__DCT_CONST_ROUNDING);
  u[7] = _mm_add_epi32(v[7], k__DCT_CONST_ROUNDING);
  u[8] = _mm_add_epi32(v[8], k__DCT_CONST_ROUNDING);
  u[9] = _mm_add_epi32(v[9], k__DCT_CONST_ROUNDING);
  u[10] = _mm_add_epi32(v[10], k__DCT_CONST_ROUNDING);
  u[11] = _mm_add_epi32(v[11], k__DCT_CONST_ROUNDING);
  u[12] = _mm_add_epi32(v[12], k__DCT_CONST_ROUNDING);
  u[13] = _mm_add_epi32(v[13], k__DCT_CONST_ROUNDING);
  u[14] = _mm_add_epi32(v[14], k__DCT_CONST_ROUNDING);
  u[15] = _mm_add_epi32(v[15], k__DCT_CONST_ROUNDING);

  u[0] = _mm_srai_epi32(u[0], DCT_CONST_BITS);
  u[1] = _mm_srai_epi32(u[1], DCT_CONST_BITS);
  u[2] = _mm_srai_epi32(u[2], DCT_CONST_BITS);
  u[3] = _mm_srai_epi32(u[3], DCT_CONST_BITS);
  u[4] = _mm_srai_epi32(u[4], DCT_CONST_BITS);
  u[5] = _mm_srai_epi32(u[5], DCT_CONST_BITS);
  u[6] = _mm_srai_epi32(u[6], DCT_CONST_BITS);
  u[7] = _mm_srai_epi32(u[7], DCT_CONST_BITS);
  u[8] = _mm_srai_epi32(u[8], DCT_CONST_BITS);
  u[9] = _mm_srai_epi32(u[9], DCT_CONST_BITS);
  u[10] = _mm_srai_epi32(u[10], DCT_CONST_BITS);
  u[11] = _mm_srai_epi32(u[11], DCT_CONST_BITS);
  u[12] = _mm_srai_epi32(u[12], DCT_CONST_BITS);
  u[13] = _mm_srai_epi32(u[13], DCT_CONST_BITS);
  u[14] = _mm_srai_epi32(u[14], DCT_CONST_BITS);
  u[15] = _mm_srai_epi32(u[15], DCT_CONST_BITS);

  s[8]  = _mm_packs_epi32(u[0], u[1]);
  s[15] = _mm_packs_epi32(u[2], u[3]);
  s[9]  = _mm_packs_epi32(u[4], u[5]);
  s[14] = _mm_packs_epi32(u[6], u[7]);
  s[10] = _mm_packs_epi32(u[8], u[9]);
  s[13] = _mm_packs_epi32(u[10], u[11]);
  s[11] = _mm_packs_epi32(u[12], u[13]);
  s[12] = _mm_packs_epi32(u[14], u[15]);

  // stage 3
  t[0] = s[0];
  t[1] = s[1];
  t[2] = s[2];
  t[3] = s[3];
  u[0] = _mm_unpacklo_epi16(s[4], s[7]);
  u[1] = _mm_unpackhi_epi16(s[4], s[7]);
  u[2] = _mm_unpacklo_epi16(s[5], s[6]);
  u[3] = _mm_unpackhi_epi16(s[5], s[6]);

  v[0] = _mm_madd_epi16(u[0], k__cospi_p28_m04);
  v[1] = _mm_madd_epi16(u[1], k__cospi_p28_m04);
  v[2] = _mm_madd_epi16(u[0], k__cospi_p04_p28);
  v[3] = _mm_madd_epi16(u[1], k__cospi_p04_p28);
  v[4] = _mm_madd_epi16(u[2], k__cospi_p12_m20);
  v[5] = _mm_madd_epi16(u[3], k__cospi_p12_m20);
  v[6] = _mm_madd_epi16(u[2], k__cospi_p20_p12);
  v[7] = _mm_madd_epi16(u[3], k__cospi_p20_p12);

  u[0] = _mm_add_epi32(v[0], k__DCT_CONST_ROUNDING);
  u[1] = _mm_add_epi32(v[1], k__DCT_CONST_ROUNDING);
  u[2] = _mm_add_epi32(v[2], k__DCT_CONST_ROUNDING);
  u[3] = _mm_add_epi32(v[3], k__DCT_CONST_ROUNDING);
  u[4] = _mm_add_epi32(v[4], k__DCT_CONST_ROUNDING);
  u[5] = _mm_add_epi32(v[5], k__DCT_CONST_ROUNDING);
  u[6] = _mm_add_epi32(v[6], k__DCT_CONST_ROUNDING);
  u[7] = _mm_add_epi32(v[7], k__DCT_CONST_ROUNDING);

  u[0] = _mm_srai_epi32(u[0], DCT_CONST_BITS);
  u[1] = _mm_srai_epi32(u[1], DCT_CONST_BITS);
  u[2] = _mm_srai_epi32(u[2], DCT_CONST_BITS);
  u[3] = _mm_srai_epi32(u[3], DCT_CONST_BITS);
  u[4] = _mm_srai_epi32(u[4], DCT_CONST_BITS);
  u[5] = _mm_srai_epi32(u[5], DCT_CONST_BITS);
  u[6] = _mm_srai_epi32(u[6], DCT_CONST_BITS);
  u[7] = _mm_srai_epi32(u[7], DCT_CONST_BITS);

  t[4] = _mm_packs_epi32(u[0], u[1]);
  t[7] = _mm_packs_epi32(u[2], u[3]);
  t[5] = _mm_packs_epi32(u[4], u[5]);
  t[6] = _mm_packs_epi32(u[6], u[7]);
  t[8] = _mm_add_epi16(s[8], s[9]);
  t[9] = _mm_sub_epi16(s[8], s[9]);
  t[10] = _mm_sub_epi16(s[11], s[10]);
  t[11] = _mm_add_epi16(s[10], s[11]);
  t[12] = _mm_add_epi16(s[12], s[13]);
  t[13] = _mm_sub_epi16(s[12], s[13]);
  t[14] = _mm_sub_epi16(s[15], s[14]);
  t[15] = _mm_add_epi16(s[14], s[15]);

  // stage 4
  u[0] = _mm_add_epi16(t[0], t[1]);
  u[1] = _mm_sub_epi16(t[0], t[1]);
  u[2] = _mm_unpacklo_epi16(t[2], t[3]);
  u[3] = _mm_unpackhi_epi16(t[2], t[3]);
  u[4] = _mm_unpacklo_epi16(t[9], t[14]);
  u[5] = _mm_unpackhi_epi16(t[9], t[14]);
  u[6] = _mm_unpacklo_epi16(t[10], t[13]);
  u[7] = _mm_unpackhi_epi16(t[10], t[13]);

  s[0] = _mm_mulhrs_epi16(u[0], k__cospi_p16_p16_x2);
  s[1] = _mm_mulhrs_epi16(u[1], k__cospi_p16_p16_x2);
  v[4] = _mm_madd_epi16(u[2], k__cospi_p24_m08);
  v[5] = _mm_madd_epi16(u[3], k__cospi_p24_m08);
  v[6] = _mm_madd_epi16(u[2], k__cospi_p08_p24);
  v[7] = _mm_madd_epi16(u[3], k__cospi_p08_p24);
  v[8] = _mm_madd_epi16(u[4], k__cospi_m08_p24);
  v[9] = _mm_madd_epi16(u[5], k__cospi_m08_p24);
  v[10] = _mm_madd_epi16(u[4], k__cospi_p24_p08);
  v[11] = _mm_madd_epi16(u[5], k__cospi_p24_p08);
  v[12] = _mm_madd_epi16(u[6], k__cospi_m24_m08);
  v[13] = _mm_madd_epi16(u[7], k__cospi_m24_m08);
  v[14] = _mm_madd_epi16(u[6], k__cospi_m08_p24);
  v[15] = _mm_madd_epi16(u[7], k__cospi_m08_p24);

  u[4] = _mm_add_epi32(v[4], k__DCT_CONST_ROUNDING);
  u[5] = _mm_add_epi32(v[5], k__DCT_CONST_ROUNDING);
  u[6] = _mm_add_epi32(v[6], k__DCT_CONST_ROUNDING);
  u[7] = _mm_add_epi32(v[7], k__DCT_CONST_ROUNDING);
  u[8] = _mm_add_epi32(v[8], k__DCT_CONST_ROUNDING);
  u[9] = _mm_add_epi32(v[9], k__DCT_CONST_ROUNDING);
  u[10] = _mm_add_epi32(v[10], k__DCT_CONST_ROUNDING);
  u[11] = _mm_add_epi32(v[11], k__DCT_CONST_ROUNDING);
  u[12] = _mm_add_epi32(v[12], k__DCT_CONST_ROUNDING);
  u[13] = _mm_add_epi32(v[13], k__DCT_CONST_ROUNDING);
  u[14] = _mm_add_epi32(v[14], k__DCT_CONST_ROUNDING);
  u[15] = _mm_add_epi32(v[15], k__DCT_CONST_ROUNDING);

  u[4] = _mm_srai_epi32(u[4], DCT_CONST_BITS);
  u[5] = _mm_srai_epi32(u[5], DCT_CONST_BITS);
  u[6] = _mm_srai_epi32(u[6], DCT_CONST_BITS);
  u[7] = _mm_srai_epi32(u[7], DCT_CONST_BITS);
  u[8] = _mm_srai_epi32(u[8], DCT_CONST_BITS);
  u[9] = _mm_srai_epi32(u[9], DCT_CONST_BITS);
  u[10] = _mm_srai_epi32(u[10], DCT_CONST_BITS);
  u[11] = _mm_srai_epi32(u[11], DCT_CONST_BITS);
  u[12] = _mm_srai_epi32(u[12], DCT_CONST_BITS);
  u[13] = _mm_srai_epi32(u[13], DCT_CONST_BITS);
  u[14] = _mm_srai_epi32(u[14], DCT_CONST_BITS);
  u[15] = _mm_srai_epi32(u[15], DCT_CONST_BITS);

  s[2] = _mm_packs_epi32(u[4], u[5]);
  s[3] = _mm_packs_epi32(u[6], u[7]);
  s[4] = _mm_add_epi16(t[4], t[5]);
  s[5] = _mm_sub_epi16(t[4], t[5]);
  s[6] = _mm_sub_epi16(t[7], t[6]);
  s[7] = _mm_add_epi16(t[6], t[7]);
  s[8] = t[8];
  s[15] = t[15];
  s[9]  = _mm_packs_epi32(u[8], u[9]);
  s[14] = _mm_packs_epi32(u[10], u[11]);
  s[10] = _mm_packs_epi32(u[12], u[13]);
  s[13] = _mm_packs_epi32(u[14], u[15]);
  s[11] = t[11];
  s[12] = t[12];

  // stage 5
  t[0] = _mm_add_epi16(s[0], s[3]);
  t[1] = _mm_add_epi16(s[1], s[2]);
  t[2] = _mm_sub_epi16(s[1], s[2]);
  t[3] = _mm_sub_epi16(s[0], s[3]);
  t[4] = s[4];
  t[7] = s[7];

  u[0] = _mm_sub_epi16(s[6], s[5]);
  u[1] = _mm_add_epi16(s[6], s[5]);
  t[5] = _mm_mulhrs_epi16(u[0], k__cospi_p16_p16_x2);
  t[6] = _mm_mulhrs_epi16(u[1], k__cospi_p16_p16_x2);

  t[8] = _mm_add_epi16(s[8], s[11]);
  t[9] = _mm_add_epi16(s[9], s[10]);
  t[10] = _mm_sub_epi16(s[9], s[10]);
  t[11] = _mm_sub_epi16(s[8], s[11]);
  t[12] = _mm_sub_epi16(s[15], s[12]);
  t[13] = _mm_sub_epi16(s[14], s[13]);
  t[14] = _mm_add_epi16(s[13], s[14]);
  t[15] = _mm_add_epi16(s[12], s[15]);

  // stage 6
  s[0] = _mm_add_epi16(t[0], t[7]);
  s[1] = _mm_add_epi16(t[1], t[6]);
  s[2] = _mm_add_epi16(t[2], t[5]);
  s[3] = _mm_add_epi16(t[3], t[4]);
  s[4] = _mm_sub_epi16(t[3], t[4]);
  s[5] = _mm_sub_epi16(t[2], t[5]);
  s[6] = _mm_sub_epi16(t[1], t[6]);
  s[7] = _mm_sub_epi16(t[0], t[7]);
  s[8] = t[8];
  s[9] = t[9];

  u[0] = _mm_sub_epi16(t[13], t[10]);
  u[1] = _mm_add_epi16(t[13], t[10]);
  u[2] = _mm_sub_epi16(t[12], t[11]);
  u[3] = _mm_add_epi16(t[12], t[11]);

  s[10] = _mm_mulhrs_epi16(u[0], k__cospi_p16_p16_x2);
  s[13] = _mm_mulhrs_epi16(u[1], k__cospi_p16_p16_x2);
  s[11] = _mm_mulhrs_epi16(u[2], k__cospi_p16_p16_x2);
  s[12] = _mm_mulhrs_epi16(u[3], k__cospi_p16_p16_x2);
  s[14] = t[14];
  s[15] = t[15];

  // stage 7
  in[0] = _mm_add_epi16(s[0], s[15]);
  in[1] = _mm_add_epi16(s[1], s[14]);
  in[2] = _mm_add_epi16(s[2], s[13]);
  in[3] = _mm_add_epi16(s[3], s[12]);
  in[4] = _mm_add_epi16(s[4], s[11]);
  in[5] = _mm_add_epi16(s[5], s[10]);
  in[6] = _mm_add_epi16(s[6], s[9]);
  in[7] = _mm_add_epi16(s[7], s[8]);
  in[8] = _mm_sub_epi16(s[7], s[8]);
  in[9] = _mm_sub_epi16(s[6], s[9]);
  in[10] = _mm_sub_epi16(s[5], s[10]);
  in[11] = _mm_sub_epi16(s[4], s[11]);
  in[12] = _mm_sub_epi16(s[3], s[12]);
  in[13] = _mm_sub_epi16(s[2], s[13]);
  in[14] = _mm_sub_epi16(s[1], s[14]);
  in[15] = _mm_sub_epi16(s[0], s[15]);
}

static void idct16_sse2(__m128i *in0, __m128i *in1) {
  array_transpose_16x16(in0, in1);
  idct16_8col(in0);
  idct16_8col(in1);
}

void vp9_idct16x16_256_add_ssse3(const int16_t *input, uint8_t *dest,
                                int stride) {
  __m128i in0[16], in1[16];

  load_buffer_8x16(input, in0);
  input += 8;
  load_buffer_8x16(input, in1);

  idct16_sse2(in0, in1);
  idct16_sse2(in0, in1);

  write_buffer_8x16(dest, in0, stride);
  dest += 8;
  write_buffer_8x16(dest, in1, stride);
}
