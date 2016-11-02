/*
 * Copyright (c) 2016, Alliance for Open Media. All rights reserved
 *
 * This source code is subject to the terms of the BSD 2 Clause License and
 * the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
 * was not distributed with this source code in the LICENSE file, you can
 * obtain it at www.aomedia.org/license/software. If the Alliance for Open
 * Media Patent License 1.0 was not distributed with this source code in the
 * PATENTS file, you can obtain it at www.aomedia.org/license/patent.
 */

#include "./av1_rtcd.h"
#include "aom_dsp/x86/inv_txfm_sse2.h"
#include "aom_dsp/x86/synonyms.h"
#include "aom_dsp/x86/txfm_common_sse2.h"
#include "aom_ports/mem.h"
#include "av1/common/enums.h"

#if CONFIG_EXT_TX
static INLINE void fliplr_4x4(__m128i in[2]) {
  in[0] = _mm_shufflelo_epi16(in[0], 0x1b);
  in[0] = _mm_shufflehi_epi16(in[0], 0x1b);
  in[1] = _mm_shufflelo_epi16(in[1], 0x1b);
  in[1] = _mm_shufflehi_epi16(in[1], 0x1b);
}

static INLINE void fliplr_8x8(__m128i in[8]) {
  in[0] = mm_reverse_epi16(in[0]);
  in[1] = mm_reverse_epi16(in[1]);
  in[2] = mm_reverse_epi16(in[2]);
  in[3] = mm_reverse_epi16(in[3]);

  in[4] = mm_reverse_epi16(in[4]);
  in[5] = mm_reverse_epi16(in[5]);
  in[6] = mm_reverse_epi16(in[6]);
  in[7] = mm_reverse_epi16(in[7]);
}

static INLINE void fliplr_16x8(__m128i in[16]) {
  fliplr_8x8(&in[0]);
  fliplr_8x8(&in[8]);
}

#define FLIPLR_16x16(in0, in1) \
  do {                         \
    __m128i *tmp;              \
    fliplr_16x8(in0);          \
    fliplr_16x8(in1);          \
    tmp = (in0);               \
    (in0) = (in1);             \
    (in1) = tmp;               \
  } while (0)

#define FLIPUD_PTR(dest, stride, size)       \
  do {                                       \
    (dest) = (dest) + ((size)-1) * (stride); \
    (stride) = -(stride);                    \
  } while (0)
#endif

void av1_iht4x4_16_add_sse2(const tran_low_t *input, uint8_t *dest, int stride,
                            int tx_type) {
  __m128i in[2];
  const __m128i zero = _mm_setzero_si128();
  const __m128i eight = _mm_set1_epi16(8);

  in[0] = load_input_data(input);
  in[1] = load_input_data(input + 8);

  switch (tx_type) {
    case DCT_DCT:
      aom_idct4_sse2(in);
      aom_idct4_sse2(in);
      break;
    case ADST_DCT:
      aom_idct4_sse2(in);
      aom_iadst4_sse2(in);
      break;
    case DCT_ADST:
      aom_iadst4_sse2(in);
      aom_idct4_sse2(in);
      break;
    case ADST_ADST:
      aom_iadst4_sse2(in);
      aom_iadst4_sse2(in);
      break;
#if CONFIG_EXT_TX
    case FLIPADST_DCT:
      aom_idct4_sse2(in);
      aom_iadst4_sse2(in);
      FLIPUD_PTR(dest, stride, 4);
      break;
    case DCT_FLIPADST:
      aom_iadst4_sse2(in);
      aom_idct4_sse2(in);
      fliplr_4x4(in);
      break;
    case FLIPADST_FLIPADST:
      aom_iadst4_sse2(in);
      aom_iadst4_sse2(in);
      FLIPUD_PTR(dest, stride, 4);
      fliplr_4x4(in);
      break;
    case ADST_FLIPADST:
      aom_iadst4_sse2(in);
      aom_iadst4_sse2(in);
      fliplr_4x4(in);
      break;
    case FLIPADST_ADST:
      aom_iadst4_sse2(in);
      aom_iadst4_sse2(in);
      FLIPUD_PTR(dest, stride, 4);
      break;
#endif  // CONFIG_EXT_TX
    default: assert(0); break;
  }

  // Final round and shift
  in[0] = _mm_add_epi16(in[0], eight);
  in[1] = _mm_add_epi16(in[1], eight);

  in[0] = _mm_srai_epi16(in[0], 4);
  in[1] = _mm_srai_epi16(in[1], 4);

  // Reconstruction and Store
  {
    __m128i d0 = _mm_cvtsi32_si128(*(const int *)(dest + stride * 0));
    __m128i d1 = _mm_cvtsi32_si128(*(const int *)(dest + stride * 1));
    __m128i d2 = _mm_cvtsi32_si128(*(const int *)(dest + stride * 2));
    __m128i d3 = _mm_cvtsi32_si128(*(const int *)(dest + stride * 3));
    d0 = _mm_unpacklo_epi32(d0, d1);
    d2 = _mm_unpacklo_epi32(d2, d3);
    d0 = _mm_unpacklo_epi8(d0, zero);
    d2 = _mm_unpacklo_epi8(d2, zero);
    d0 = _mm_add_epi16(d0, in[0]);
    d2 = _mm_add_epi16(d2, in[1]);
    d0 = _mm_packus_epi16(d0, d2);
    // store result[0]
    *(int *)dest = _mm_cvtsi128_si32(d0);
    // store result[1]
    d0 = _mm_srli_si128(d0, 4);
    *(int *)(dest + stride) = _mm_cvtsi128_si32(d0);
    // store result[2]
    d0 = _mm_srli_si128(d0, 4);
    *(int *)(dest + stride * 2) = _mm_cvtsi128_si32(d0);
    // store result[3]
    d0 = _mm_srli_si128(d0, 4);
    *(int *)(dest + stride * 3) = _mm_cvtsi128_si32(d0);
  }
}

void av1_iht8x8_64_add_sse2(const tran_low_t *input, uint8_t *dest, int stride,
                            int tx_type) {
  __m128i in[8];
  const __m128i zero = _mm_setzero_si128();
  const __m128i final_rounding = _mm_set1_epi16(1 << 4);

  // load input data
  in[0] = load_input_data(input);
  in[1] = load_input_data(input + 8 * 1);
  in[2] = load_input_data(input + 8 * 2);
  in[3] = load_input_data(input + 8 * 3);
  in[4] = load_input_data(input + 8 * 4);
  in[5] = load_input_data(input + 8 * 5);
  in[6] = load_input_data(input + 8 * 6);
  in[7] = load_input_data(input + 8 * 7);

  switch (tx_type) {
    case DCT_DCT:
      aom_idct8_sse2(in);
      aom_idct8_sse2(in);
      break;
    case ADST_DCT:
      aom_idct8_sse2(in);
      aom_iadst8_sse2(in);
      break;
    case DCT_ADST:
      aom_iadst8_sse2(in);
      aom_idct8_sse2(in);
      break;
    case ADST_ADST:
      aom_iadst8_sse2(in);
      aom_iadst8_sse2(in);
      break;
#if CONFIG_EXT_TX
    case FLIPADST_DCT:
      aom_idct8_sse2(in);
      aom_iadst8_sse2(in);
      FLIPUD_PTR(dest, stride, 8);
      break;
    case DCT_FLIPADST:
      aom_iadst8_sse2(in);
      aom_idct8_sse2(in);
      fliplr_8x8(in);
      break;
    case FLIPADST_FLIPADST:
      aom_iadst8_sse2(in);
      aom_iadst8_sse2(in);
      FLIPUD_PTR(dest, stride, 8);
      fliplr_8x8(in);
      break;
    case ADST_FLIPADST:
      aom_iadst8_sse2(in);
      aom_iadst8_sse2(in);
      fliplr_8x8(in);
      break;
    case FLIPADST_ADST:
      aom_iadst8_sse2(in);
      aom_iadst8_sse2(in);
      FLIPUD_PTR(dest, stride, 8);
      break;
#endif  // CONFIG_EXT_TX
    default: assert(0); break;
  }

  // Final rounding and shift
  in[0] = _mm_adds_epi16(in[0], final_rounding);
  in[1] = _mm_adds_epi16(in[1], final_rounding);
  in[2] = _mm_adds_epi16(in[2], final_rounding);
  in[3] = _mm_adds_epi16(in[3], final_rounding);
  in[4] = _mm_adds_epi16(in[4], final_rounding);
  in[5] = _mm_adds_epi16(in[5], final_rounding);
  in[6] = _mm_adds_epi16(in[6], final_rounding);
  in[7] = _mm_adds_epi16(in[7], final_rounding);

  in[0] = _mm_srai_epi16(in[0], 5);
  in[1] = _mm_srai_epi16(in[1], 5);
  in[2] = _mm_srai_epi16(in[2], 5);
  in[3] = _mm_srai_epi16(in[3], 5);
  in[4] = _mm_srai_epi16(in[4], 5);
  in[5] = _mm_srai_epi16(in[5], 5);
  in[6] = _mm_srai_epi16(in[6], 5);
  in[7] = _mm_srai_epi16(in[7], 5);

  RECON_AND_STORE(dest + 0 * stride, in[0]);
  RECON_AND_STORE(dest + 1 * stride, in[1]);
  RECON_AND_STORE(dest + 2 * stride, in[2]);
  RECON_AND_STORE(dest + 3 * stride, in[3]);
  RECON_AND_STORE(dest + 4 * stride, in[4]);
  RECON_AND_STORE(dest + 5 * stride, in[5]);
  RECON_AND_STORE(dest + 6 * stride, in[6]);
  RECON_AND_STORE(dest + 7 * stride, in[7]);
}

#if CONFIG_EXT_TX
static void iidtx16_8col(__m128i *in) {
  const __m128i k__zero_epi16 = _mm_set1_epi16((int16_t)0);
  const __m128i k__sqrt2_epi16 = _mm_set1_epi16((int16_t)Sqrt2);
  const __m128i k__DCT_CONST_ROUNDING = _mm_set1_epi32(DCT_CONST_ROUNDING);

  __m128i v0, v1, v2, v3, v4, v5, v6, v7;
  __m128i u0, u1, u2, u3, u4, u5, u6, u7;
  __m128i x0, x1, x2, x3, x4, x5, x6, x7;
  __m128i y0, y1, y2, y3, y4, y5, y6, y7;

  in[0] = _mm_slli_epi16(in[0], 1);
  in[1] = _mm_slli_epi16(in[1], 1);
  in[2] = _mm_slli_epi16(in[2], 1);
  in[3] = _mm_slli_epi16(in[3], 1);
  in[4] = _mm_slli_epi16(in[4], 1);
  in[5] = _mm_slli_epi16(in[5], 1);
  in[6] = _mm_slli_epi16(in[6], 1);
  in[7] = _mm_slli_epi16(in[7], 1);
  in[8] = _mm_slli_epi16(in[8], 1);
  in[9] = _mm_slli_epi16(in[9], 1);
  in[10] = _mm_slli_epi16(in[10], 1);
  in[11] = _mm_slli_epi16(in[11], 1);
  in[12] = _mm_slli_epi16(in[12], 1);
  in[13] = _mm_slli_epi16(in[13], 1);
  in[14] = _mm_slli_epi16(in[14], 1);
  in[15] = _mm_slli_epi16(in[15], 1);

  v0 = _mm_unpacklo_epi16(in[0], k__zero_epi16);
  v1 = _mm_unpacklo_epi16(in[1], k__zero_epi16);
  v2 = _mm_unpacklo_epi16(in[2], k__zero_epi16);
  v3 = _mm_unpacklo_epi16(in[3], k__zero_epi16);
  v4 = _mm_unpacklo_epi16(in[4], k__zero_epi16);
  v5 = _mm_unpacklo_epi16(in[5], k__zero_epi16);
  v6 = _mm_unpacklo_epi16(in[6], k__zero_epi16);
  v7 = _mm_unpacklo_epi16(in[7], k__zero_epi16);

  u0 = _mm_unpacklo_epi16(in[8], k__zero_epi16);
  u1 = _mm_unpacklo_epi16(in[9], k__zero_epi16);
  u2 = _mm_unpacklo_epi16(in[10], k__zero_epi16);
  u3 = _mm_unpacklo_epi16(in[11], k__zero_epi16);
  u4 = _mm_unpacklo_epi16(in[12], k__zero_epi16);
  u5 = _mm_unpacklo_epi16(in[13], k__zero_epi16);
  u6 = _mm_unpacklo_epi16(in[14], k__zero_epi16);
  u7 = _mm_unpacklo_epi16(in[15], k__zero_epi16);

  x0 = _mm_unpackhi_epi16(in[0], k__zero_epi16);
  x1 = _mm_unpackhi_epi16(in[1], k__zero_epi16);
  x2 = _mm_unpackhi_epi16(in[2], k__zero_epi16);
  x3 = _mm_unpackhi_epi16(in[3], k__zero_epi16);
  x4 = _mm_unpackhi_epi16(in[4], k__zero_epi16);
  x5 = _mm_unpackhi_epi16(in[5], k__zero_epi16);
  x6 = _mm_unpackhi_epi16(in[6], k__zero_epi16);
  x7 = _mm_unpackhi_epi16(in[7], k__zero_epi16);

  y0 = _mm_unpackhi_epi16(in[8], k__zero_epi16);
  y1 = _mm_unpackhi_epi16(in[9], k__zero_epi16);
  y2 = _mm_unpackhi_epi16(in[10], k__zero_epi16);
  y3 = _mm_unpackhi_epi16(in[11], k__zero_epi16);
  y4 = _mm_unpackhi_epi16(in[12], k__zero_epi16);
  y5 = _mm_unpackhi_epi16(in[13], k__zero_epi16);
  y6 = _mm_unpackhi_epi16(in[14], k__zero_epi16);
  y7 = _mm_unpackhi_epi16(in[15], k__zero_epi16);

  v0 = _mm_madd_epi16(v0, k__sqrt2_epi16);
  v1 = _mm_madd_epi16(v1, k__sqrt2_epi16);
  v2 = _mm_madd_epi16(v2, k__sqrt2_epi16);
  v3 = _mm_madd_epi16(v3, k__sqrt2_epi16);
  v4 = _mm_madd_epi16(v4, k__sqrt2_epi16);
  v5 = _mm_madd_epi16(v5, k__sqrt2_epi16);
  v6 = _mm_madd_epi16(v6, k__sqrt2_epi16);
  v7 = _mm_madd_epi16(v7, k__sqrt2_epi16);

  x0 = _mm_madd_epi16(x0, k__sqrt2_epi16);
  x1 = _mm_madd_epi16(x1, k__sqrt2_epi16);
  x2 = _mm_madd_epi16(x2, k__sqrt2_epi16);
  x3 = _mm_madd_epi16(x3, k__sqrt2_epi16);
  x4 = _mm_madd_epi16(x4, k__sqrt2_epi16);
  x5 = _mm_madd_epi16(x5, k__sqrt2_epi16);
  x6 = _mm_madd_epi16(x6, k__sqrt2_epi16);
  x7 = _mm_madd_epi16(x7, k__sqrt2_epi16);

  u0 = _mm_madd_epi16(u0, k__sqrt2_epi16);
  u1 = _mm_madd_epi16(u1, k__sqrt2_epi16);
  u2 = _mm_madd_epi16(u2, k__sqrt2_epi16);
  u3 = _mm_madd_epi16(u3, k__sqrt2_epi16);
  u4 = _mm_madd_epi16(u4, k__sqrt2_epi16);
  u5 = _mm_madd_epi16(u5, k__sqrt2_epi16);
  u6 = _mm_madd_epi16(u6, k__sqrt2_epi16);
  u7 = _mm_madd_epi16(u7, k__sqrt2_epi16);

  y0 = _mm_madd_epi16(y0, k__sqrt2_epi16);
  y1 = _mm_madd_epi16(y1, k__sqrt2_epi16);
  y2 = _mm_madd_epi16(y2, k__sqrt2_epi16);
  y3 = _mm_madd_epi16(y3, k__sqrt2_epi16);
  y4 = _mm_madd_epi16(y4, k__sqrt2_epi16);
  y5 = _mm_madd_epi16(y5, k__sqrt2_epi16);
  y6 = _mm_madd_epi16(y6, k__sqrt2_epi16);
  y7 = _mm_madd_epi16(y7, k__sqrt2_epi16);

  v0 = _mm_add_epi32(v0, k__DCT_CONST_ROUNDING);
  v1 = _mm_add_epi32(v1, k__DCT_CONST_ROUNDING);
  v2 = _mm_add_epi32(v2, k__DCT_CONST_ROUNDING);
  v3 = _mm_add_epi32(v3, k__DCT_CONST_ROUNDING);
  v4 = _mm_add_epi32(v4, k__DCT_CONST_ROUNDING);
  v5 = _mm_add_epi32(v5, k__DCT_CONST_ROUNDING);
  v6 = _mm_add_epi32(v6, k__DCT_CONST_ROUNDING);
  v7 = _mm_add_epi32(v7, k__DCT_CONST_ROUNDING);

  x0 = _mm_add_epi32(x0, k__DCT_CONST_ROUNDING);
  x1 = _mm_add_epi32(x1, k__DCT_CONST_ROUNDING);
  x2 = _mm_add_epi32(x2, k__DCT_CONST_ROUNDING);
  x3 = _mm_add_epi32(x3, k__DCT_CONST_ROUNDING);
  x4 = _mm_add_epi32(x4, k__DCT_CONST_ROUNDING);
  x5 = _mm_add_epi32(x5, k__DCT_CONST_ROUNDING);
  x6 = _mm_add_epi32(x6, k__DCT_CONST_ROUNDING);
  x7 = _mm_add_epi32(x7, k__DCT_CONST_ROUNDING);

  u0 = _mm_add_epi32(u0, k__DCT_CONST_ROUNDING);
  u1 = _mm_add_epi32(u1, k__DCT_CONST_ROUNDING);
  u2 = _mm_add_epi32(u2, k__DCT_CONST_ROUNDING);
  u3 = _mm_add_epi32(u3, k__DCT_CONST_ROUNDING);
  u4 = _mm_add_epi32(u4, k__DCT_CONST_ROUNDING);
  u5 = _mm_add_epi32(u5, k__DCT_CONST_ROUNDING);
  u6 = _mm_add_epi32(u6, k__DCT_CONST_ROUNDING);
  u7 = _mm_add_epi32(u7, k__DCT_CONST_ROUNDING);

  y0 = _mm_add_epi32(y0, k__DCT_CONST_ROUNDING);
  y1 = _mm_add_epi32(y1, k__DCT_CONST_ROUNDING);
  y2 = _mm_add_epi32(y2, k__DCT_CONST_ROUNDING);
  y3 = _mm_add_epi32(y3, k__DCT_CONST_ROUNDING);
  y4 = _mm_add_epi32(y4, k__DCT_CONST_ROUNDING);
  y5 = _mm_add_epi32(y5, k__DCT_CONST_ROUNDING);
  y6 = _mm_add_epi32(y6, k__DCT_CONST_ROUNDING);
  y7 = _mm_add_epi32(y7, k__DCT_CONST_ROUNDING);

  v0 = _mm_srai_epi32(v0, DCT_CONST_BITS);
  v1 = _mm_srai_epi32(v1, DCT_CONST_BITS);
  v2 = _mm_srai_epi32(v2, DCT_CONST_BITS);
  v3 = _mm_srai_epi32(v3, DCT_CONST_BITS);
  v4 = _mm_srai_epi32(v4, DCT_CONST_BITS);
  v5 = _mm_srai_epi32(v5, DCT_CONST_BITS);
  v6 = _mm_srai_epi32(v6, DCT_CONST_BITS);
  v7 = _mm_srai_epi32(v7, DCT_CONST_BITS);

  x0 = _mm_srai_epi32(x0, DCT_CONST_BITS);
  x1 = _mm_srai_epi32(x1, DCT_CONST_BITS);
  x2 = _mm_srai_epi32(x2, DCT_CONST_BITS);
  x3 = _mm_srai_epi32(x3, DCT_CONST_BITS);
  x4 = _mm_srai_epi32(x4, DCT_CONST_BITS);
  x5 = _mm_srai_epi32(x5, DCT_CONST_BITS);
  x6 = _mm_srai_epi32(x6, DCT_CONST_BITS);
  x7 = _mm_srai_epi32(x7, DCT_CONST_BITS);

  u0 = _mm_srai_epi32(u0, DCT_CONST_BITS);
  u1 = _mm_srai_epi32(u1, DCT_CONST_BITS);
  u2 = _mm_srai_epi32(u2, DCT_CONST_BITS);
  u3 = _mm_srai_epi32(u3, DCT_CONST_BITS);
  u4 = _mm_srai_epi32(u4, DCT_CONST_BITS);
  u5 = _mm_srai_epi32(u5, DCT_CONST_BITS);
  u6 = _mm_srai_epi32(u6, DCT_CONST_BITS);
  u7 = _mm_srai_epi32(u7, DCT_CONST_BITS);

  y0 = _mm_srai_epi32(y0, DCT_CONST_BITS);
  y1 = _mm_srai_epi32(y1, DCT_CONST_BITS);
  y2 = _mm_srai_epi32(y2, DCT_CONST_BITS);
  y3 = _mm_srai_epi32(y3, DCT_CONST_BITS);
  y4 = _mm_srai_epi32(y4, DCT_CONST_BITS);
  y5 = _mm_srai_epi32(y5, DCT_CONST_BITS);
  y6 = _mm_srai_epi32(y6, DCT_CONST_BITS);
  y7 = _mm_srai_epi32(y7, DCT_CONST_BITS);

  in[0] = _mm_packs_epi32(v0, x0);
  in[1] = _mm_packs_epi32(v1, x1);
  in[2] = _mm_packs_epi32(v2, x2);
  in[3] = _mm_packs_epi32(v3, x3);
  in[4] = _mm_packs_epi32(v4, x4);
  in[5] = _mm_packs_epi32(v5, x5);
  in[6] = _mm_packs_epi32(v6, x6);
  in[7] = _mm_packs_epi32(v7, x7);

  in[8] = _mm_packs_epi32(u0, y0);
  in[9] = _mm_packs_epi32(u1, y1);
  in[10] = _mm_packs_epi32(u2, y2);
  in[11] = _mm_packs_epi32(u3, y3);
  in[12] = _mm_packs_epi32(u4, y4);
  in[13] = _mm_packs_epi32(u5, y5);
  in[14] = _mm_packs_epi32(u6, y6);
  in[15] = _mm_packs_epi32(u7, y7);
}

static void iidtx16_sse2(__m128i *in0, __m128i *in1) {
  array_transpose_16x16(in0, in1);
  iidtx16_8col(in0);
  iidtx16_8col(in1);
}
#endif

void av1_iht16x16_256_add_sse2(const tran_low_t *input, uint8_t *dest,
                               int stride, int tx_type) {
  __m128i in[32];
  __m128i *in0 = &in[0];
  __m128i *in1 = &in[16];

  load_buffer_8x16(input, in0);
  input += 8;
  load_buffer_8x16(input, in1);

  switch (tx_type) {
    case DCT_DCT:
      aom_idct16_sse2(in0, in1);
      aom_idct16_sse2(in0, in1);
      break;
    case ADST_DCT:
      aom_idct16_sse2(in0, in1);
      aom_iadst16_sse2(in0, in1);
      break;
    case DCT_ADST:
      aom_iadst16_sse2(in0, in1);
      aom_idct16_sse2(in0, in1);
      break;
    case ADST_ADST:
      aom_iadst16_sse2(in0, in1);
      aom_iadst16_sse2(in0, in1);
      break;
#if CONFIG_EXT_TX
    case FLIPADST_DCT:
      aom_idct16_sse2(in0, in1);
      aom_iadst16_sse2(in0, in1);
      FLIPUD_PTR(dest, stride, 16);
      break;
    case DCT_FLIPADST:
      aom_iadst16_sse2(in0, in1);
      aom_idct16_sse2(in0, in1);
      FLIPLR_16x16(in0, in1);
      break;
    case FLIPADST_FLIPADST:
      aom_iadst16_sse2(in0, in1);
      aom_iadst16_sse2(in0, in1);
      FLIPUD_PTR(dest, stride, 16);
      FLIPLR_16x16(in0, in1);
      break;
    case ADST_FLIPADST:
      aom_iadst16_sse2(in0, in1);
      aom_iadst16_sse2(in0, in1);
      FLIPLR_16x16(in0, in1);
      break;
    case FLIPADST_ADST:
      aom_iadst16_sse2(in0, in1);
      aom_iadst16_sse2(in0, in1);
      FLIPUD_PTR(dest, stride, 16);
      break;
    case V_DCT:
      iidtx16_sse2(in0, in1);
      aom_idct16_sse2(in0, in1);
      break;
    case H_DCT:
      aom_idct16_sse2(in0, in1);
      iidtx16_sse2(in0, in1);
      break;
    case V_ADST:
      iidtx16_sse2(in0, in1);
      aom_iadst16_sse2(in0, in1);
      break;
    case H_ADST:
      aom_iadst16_sse2(in0, in1);
      iidtx16_sse2(in0, in1);
      break;
    case V_FLIPADST:
      iidtx16_sse2(in0, in1);
      aom_iadst16_sse2(in0, in1);
      FLIPUD_PTR(dest, stride, 16);
      break;
    case H_FLIPADST:
      aom_iadst16_sse2(in0, in1);
      iidtx16_sse2(in0, in1);
      FLIPLR_16x16(in0, in1);
      break;
#endif  // CONFIG_EXT_TX
    default: assert(0); break;
  }

  write_buffer_8x16(dest, in0, stride);
  dest += 8;
  write_buffer_8x16(dest, in1, stride);
}

#if CONFIG_EXT_TX
static void iidtx8_sse2(__m128i *in) {
  in[0] = _mm_slli_epi16(in[0], 1);
  in[1] = _mm_slli_epi16(in[1], 1);
  in[2] = _mm_slli_epi16(in[2], 1);
  in[3] = _mm_slli_epi16(in[3], 1);
  in[4] = _mm_slli_epi16(in[4], 1);
  in[5] = _mm_slli_epi16(in[5], 1);
  in[6] = _mm_slli_epi16(in[6], 1);
  in[7] = _mm_slli_epi16(in[7], 1);
}

static INLINE void iidtx4_sse2(__m128i *in) {
  const __m128i v_scale_w = _mm_set1_epi16(Sqrt2);

  const __m128i v_p0l_w = _mm_mullo_epi16(in[0], v_scale_w);
  const __m128i v_p0h_w = _mm_mulhi_epi16(in[0], v_scale_w);
  const __m128i v_p1l_w = _mm_mullo_epi16(in[1], v_scale_w);
  const __m128i v_p1h_w = _mm_mulhi_epi16(in[1], v_scale_w);

  const __m128i v_p0a_d = _mm_unpacklo_epi16(v_p0l_w, v_p0h_w);
  const __m128i v_p0b_d = _mm_unpackhi_epi16(v_p0l_w, v_p0h_w);
  const __m128i v_p1a_d = _mm_unpacklo_epi16(v_p1l_w, v_p1h_w);
  const __m128i v_p1b_d = _mm_unpackhi_epi16(v_p1l_w, v_p1h_w);

  in[0] = _mm_packs_epi32(xx_roundn_epi32_unsigned(v_p0a_d, DCT_CONST_BITS),
                          xx_roundn_epi32_unsigned(v_p0b_d, DCT_CONST_BITS));
  in[1] = _mm_packs_epi32(xx_roundn_epi32_unsigned(v_p1a_d, DCT_CONST_BITS),
                          xx_roundn_epi32_unsigned(v_p1b_d, DCT_CONST_BITS));
}

// load 8x8 array
static INLINE void flip_buffer_lr_8x8(__m128i *in) {
  in[0] = mm_reverse_epi16(in[0]);
  in[1] = mm_reverse_epi16(in[1]);
  in[2] = mm_reverse_epi16(in[2]);
  in[3] = mm_reverse_epi16(in[3]);
  in[4] = mm_reverse_epi16(in[4]);
  in[5] = mm_reverse_epi16(in[5]);
  in[6] = mm_reverse_epi16(in[6]);
  in[7] = mm_reverse_epi16(in[7]);
}
#endif  // CONFIG_EXT_TX

static INLINE void scale_sqrt2_8x4(__m128i *in) {
  // Implements 'ROUND_POWER_OF_TWO(input * Sqrt2, DCT_CONST_BITS)'
  // for each element
  const __m128i v_scale_w = _mm_set1_epi16(Sqrt2);

  const __m128i v_p0l_w = _mm_mullo_epi16(in[0], v_scale_w);
  const __m128i v_p0h_w = _mm_mulhi_epi16(in[0], v_scale_w);
  const __m128i v_p1l_w = _mm_mullo_epi16(in[1], v_scale_w);
  const __m128i v_p1h_w = _mm_mulhi_epi16(in[1], v_scale_w);
  const __m128i v_p2l_w = _mm_mullo_epi16(in[2], v_scale_w);
  const __m128i v_p2h_w = _mm_mulhi_epi16(in[2], v_scale_w);
  const __m128i v_p3l_w = _mm_mullo_epi16(in[3], v_scale_w);
  const __m128i v_p3h_w = _mm_mulhi_epi16(in[3], v_scale_w);

  const __m128i v_p0a_d = _mm_unpacklo_epi16(v_p0l_w, v_p0h_w);
  const __m128i v_p0b_d = _mm_unpackhi_epi16(v_p0l_w, v_p0h_w);
  const __m128i v_p1a_d = _mm_unpacklo_epi16(v_p1l_w, v_p1h_w);
  const __m128i v_p1b_d = _mm_unpackhi_epi16(v_p1l_w, v_p1h_w);
  const __m128i v_p2a_d = _mm_unpacklo_epi16(v_p2l_w, v_p2h_w);
  const __m128i v_p2b_d = _mm_unpackhi_epi16(v_p2l_w, v_p2h_w);
  const __m128i v_p3a_d = _mm_unpacklo_epi16(v_p3l_w, v_p3h_w);
  const __m128i v_p3b_d = _mm_unpackhi_epi16(v_p3l_w, v_p3h_w);

  in[0] = _mm_packs_epi32(xx_roundn_epi32_unsigned(v_p0a_d, DCT_CONST_BITS),
                          xx_roundn_epi32_unsigned(v_p0b_d, DCT_CONST_BITS));
  in[1] = _mm_packs_epi32(xx_roundn_epi32_unsigned(v_p1a_d, DCT_CONST_BITS),
                          xx_roundn_epi32_unsigned(v_p1b_d, DCT_CONST_BITS));
  in[2] = _mm_packs_epi32(xx_roundn_epi32_unsigned(v_p2a_d, DCT_CONST_BITS),
                          xx_roundn_epi32_unsigned(v_p2b_d, DCT_CONST_BITS));
  in[3] = _mm_packs_epi32(xx_roundn_epi32_unsigned(v_p3a_d, DCT_CONST_BITS),
                          xx_roundn_epi32_unsigned(v_p3b_d, DCT_CONST_BITS));
}

static INLINE void scale_sqrt2_8x8(__m128i *in) {
  // Implements 'ROUND_POWER_OF_TWO_SIGNED(input * Sqrt2, DCT_CONST_BITS)'
  // for each element
  const __m128i v_scale_w = _mm_set1_epi16(Sqrt2);

  const __m128i v_p0l_w = _mm_mullo_epi16(in[0], v_scale_w);
  const __m128i v_p0h_w = _mm_mulhi_epi16(in[0], v_scale_w);
  const __m128i v_p1l_w = _mm_mullo_epi16(in[1], v_scale_w);
  const __m128i v_p1h_w = _mm_mulhi_epi16(in[1], v_scale_w);
  const __m128i v_p2l_w = _mm_mullo_epi16(in[2], v_scale_w);
  const __m128i v_p2h_w = _mm_mulhi_epi16(in[2], v_scale_w);
  const __m128i v_p3l_w = _mm_mullo_epi16(in[3], v_scale_w);
  const __m128i v_p3h_w = _mm_mulhi_epi16(in[3], v_scale_w);
  const __m128i v_p4l_w = _mm_mullo_epi16(in[4], v_scale_w);
  const __m128i v_p4h_w = _mm_mulhi_epi16(in[4], v_scale_w);
  const __m128i v_p5l_w = _mm_mullo_epi16(in[5], v_scale_w);
  const __m128i v_p5h_w = _mm_mulhi_epi16(in[5], v_scale_w);
  const __m128i v_p6l_w = _mm_mullo_epi16(in[6], v_scale_w);
  const __m128i v_p6h_w = _mm_mulhi_epi16(in[6], v_scale_w);
  const __m128i v_p7l_w = _mm_mullo_epi16(in[7], v_scale_w);
  const __m128i v_p7h_w = _mm_mulhi_epi16(in[7], v_scale_w);

  const __m128i v_p0a_d = _mm_unpacklo_epi16(v_p0l_w, v_p0h_w);
  const __m128i v_p0b_d = _mm_unpackhi_epi16(v_p0l_w, v_p0h_w);
  const __m128i v_p1a_d = _mm_unpacklo_epi16(v_p1l_w, v_p1h_w);
  const __m128i v_p1b_d = _mm_unpackhi_epi16(v_p1l_w, v_p1h_w);
  const __m128i v_p2a_d = _mm_unpacklo_epi16(v_p2l_w, v_p2h_w);
  const __m128i v_p2b_d = _mm_unpackhi_epi16(v_p2l_w, v_p2h_w);
  const __m128i v_p3a_d = _mm_unpacklo_epi16(v_p3l_w, v_p3h_w);
  const __m128i v_p3b_d = _mm_unpackhi_epi16(v_p3l_w, v_p3h_w);
  const __m128i v_p4a_d = _mm_unpacklo_epi16(v_p4l_w, v_p4h_w);
  const __m128i v_p4b_d = _mm_unpackhi_epi16(v_p4l_w, v_p4h_w);
  const __m128i v_p5a_d = _mm_unpacklo_epi16(v_p5l_w, v_p5h_w);
  const __m128i v_p5b_d = _mm_unpackhi_epi16(v_p5l_w, v_p5h_w);
  const __m128i v_p6a_d = _mm_unpacklo_epi16(v_p6l_w, v_p6h_w);
  const __m128i v_p6b_d = _mm_unpackhi_epi16(v_p6l_w, v_p6h_w);
  const __m128i v_p7a_d = _mm_unpacklo_epi16(v_p7l_w, v_p7h_w);
  const __m128i v_p7b_d = _mm_unpackhi_epi16(v_p7l_w, v_p7h_w);

  in[0] = _mm_packs_epi32(xx_roundn_epi32_unsigned(v_p0a_d, DCT_CONST_BITS),
                          xx_roundn_epi32_unsigned(v_p0b_d, DCT_CONST_BITS));
  in[1] = _mm_packs_epi32(xx_roundn_epi32_unsigned(v_p1a_d, DCT_CONST_BITS),
                          xx_roundn_epi32_unsigned(v_p1b_d, DCT_CONST_BITS));
  in[2] = _mm_packs_epi32(xx_roundn_epi32_unsigned(v_p2a_d, DCT_CONST_BITS),
                          xx_roundn_epi32_unsigned(v_p2b_d, DCT_CONST_BITS));
  in[3] = _mm_packs_epi32(xx_roundn_epi32_unsigned(v_p3a_d, DCT_CONST_BITS),
                          xx_roundn_epi32_unsigned(v_p3b_d, DCT_CONST_BITS));
  in[4] = _mm_packs_epi32(xx_roundn_epi32_unsigned(v_p4a_d, DCT_CONST_BITS),
                          xx_roundn_epi32_unsigned(v_p4b_d, DCT_CONST_BITS));
  in[5] = _mm_packs_epi32(xx_roundn_epi32_unsigned(v_p5a_d, DCT_CONST_BITS),
                          xx_roundn_epi32_unsigned(v_p5b_d, DCT_CONST_BITS));
  in[6] = _mm_packs_epi32(xx_roundn_epi32_unsigned(v_p6a_d, DCT_CONST_BITS),
                          xx_roundn_epi32_unsigned(v_p6b_d, DCT_CONST_BITS));
  in[7] = _mm_packs_epi32(xx_roundn_epi32_unsigned(v_p7a_d, DCT_CONST_BITS),
                          xx_roundn_epi32_unsigned(v_p7b_d, DCT_CONST_BITS));
}

static INLINE void scale_sqrt2_8x16(__m128i *in) {
  scale_sqrt2_8x8(in);
  scale_sqrt2_8x8(in + 8);
}

void av1_iht8x16_128_add_sse2(const tran_low_t *input, uint8_t *dest,
                              int stride, int tx_type) {
  __m128i in[16];

  in[0] = load_input_data(input + 0 * 8);
  in[1] = load_input_data(input + 1 * 8);
  in[2] = load_input_data(input + 2 * 8);
  in[3] = load_input_data(input + 3 * 8);
  in[4] = load_input_data(input + 4 * 8);
  in[5] = load_input_data(input + 5 * 8);
  in[6] = load_input_data(input + 6 * 8);
  in[7] = load_input_data(input + 7 * 8);

  in[8] = load_input_data(input + 8 * 8);
  in[9] = load_input_data(input + 9 * 8);
  in[10] = load_input_data(input + 10 * 8);
  in[11] = load_input_data(input + 11 * 8);
  in[12] = load_input_data(input + 12 * 8);
  in[13] = load_input_data(input + 13 * 8);
  in[14] = load_input_data(input + 14 * 8);
  in[15] = load_input_data(input + 15 * 8);

  // Row transform
  switch (tx_type) {
    case DCT_DCT:
    case ADST_DCT:
#if CONFIG_EXT_TX
    case FLIPADST_DCT:
    case H_DCT:
#endif
      aom_idct8_sse2(in);
      array_transpose_8x8(in, in);
      aom_idct8_sse2(in + 8);
      array_transpose_8x8(in + 8, in + 8);
      break;
    case DCT_ADST:
    case ADST_ADST:
#if CONFIG_EXT_TX
    case DCT_FLIPADST:
    case FLIPADST_FLIPADST:
    case ADST_FLIPADST:
    case FLIPADST_ADST:
    case H_ADST:
    case H_FLIPADST:
#endif
      aom_iadst8_sse2(in);
      array_transpose_8x8(in, in);
      aom_iadst8_sse2(in + 8);
      array_transpose_8x8(in + 8, in + 8);
      break;
#if CONFIG_EXT_TX
    case V_FLIPADST:
    case V_ADST:
    case V_DCT:
    case IDTX:
      iidtx8_sse2(in);
      iidtx8_sse2(in + 8);
      break;
#endif
    default: assert(0); break;
  }
  scale_sqrt2_8x8(in);
  scale_sqrt2_8x8(in + 8);

  // Column transform
  switch (tx_type) {
    case DCT_DCT:
    case DCT_ADST:
#if CONFIG_EXT_TX
    case DCT_FLIPADST:
    case V_DCT:
#endif
      idct16_8col(in);
      break;
    case ADST_DCT:
    case ADST_ADST:
#if CONFIG_EXT_TX
    case FLIPADST_ADST:
    case ADST_FLIPADST:
    case FLIPADST_FLIPADST:
    case FLIPADST_DCT:
    case V_ADST:
    case V_FLIPADST:
#endif
      iadst16_8col(in);
      break;
#if CONFIG_EXT_TX
    case H_DCT:
    case H_ADST:
    case H_FLIPADST:
    case IDTX: iidtx16_8col(in); break;
#endif
    default: assert(0); break;
  }

  switch (tx_type) {
    case DCT_DCT:
    case ADST_DCT:
#if CONFIG_EXT_TX
    case H_DCT:
#endif
    case DCT_ADST:
    case ADST_ADST:
#if CONFIG_EXT_TX
    case H_ADST:
    case V_ADST:
    case V_DCT:
    case IDTX:
#endif
      write_buffer_8x16(dest, in, stride);
      break;
#if CONFIG_EXT_TX
    case FLIPADST_DCT:
    case FLIPADST_ADST:
    case V_FLIPADST: write_buffer_8x16(dest + stride * 15, in, -stride); break;
    case DCT_FLIPADST:
    case ADST_FLIPADST:
    case H_FLIPADST:
      flip_buffer_lr_8x8(in);
      flip_buffer_lr_8x8(in + 8);
      write_buffer_8x16(dest, in, stride);
      break;
    case FLIPADST_FLIPADST:
      flip_buffer_lr_8x8(in);
      flip_buffer_lr_8x8(in + 8);
      write_buffer_8x16(dest + stride * 15, in, -stride);
      break;
#endif
    default: assert(0); break;
  }
}

static INLINE void write_buffer_8x8_round6(uint8_t *dest, __m128i *in,
                                           int stride) {
  const __m128i final_rounding = _mm_set1_epi16(1 << 5);
  const __m128i zero = _mm_setzero_si128();
  // Final rounding and shift
  in[0] = _mm_adds_epi16(in[0], final_rounding);
  in[1] = _mm_adds_epi16(in[1], final_rounding);
  in[2] = _mm_adds_epi16(in[2], final_rounding);
  in[3] = _mm_adds_epi16(in[3], final_rounding);
  in[4] = _mm_adds_epi16(in[4], final_rounding);
  in[5] = _mm_adds_epi16(in[5], final_rounding);
  in[6] = _mm_adds_epi16(in[6], final_rounding);
  in[7] = _mm_adds_epi16(in[7], final_rounding);

  in[0] = _mm_srai_epi16(in[0], 6);
  in[1] = _mm_srai_epi16(in[1], 6);
  in[2] = _mm_srai_epi16(in[2], 6);
  in[3] = _mm_srai_epi16(in[3], 6);
  in[4] = _mm_srai_epi16(in[4], 6);
  in[5] = _mm_srai_epi16(in[5], 6);
  in[6] = _mm_srai_epi16(in[6], 6);
  in[7] = _mm_srai_epi16(in[7], 6);

  RECON_AND_STORE(dest + 0 * stride, in[0]);
  RECON_AND_STORE(dest + 1 * stride, in[1]);
  RECON_AND_STORE(dest + 2 * stride, in[2]);
  RECON_AND_STORE(dest + 3 * stride, in[3]);
  RECON_AND_STORE(dest + 4 * stride, in[4]);
  RECON_AND_STORE(dest + 5 * stride, in[5]);
  RECON_AND_STORE(dest + 6 * stride, in[6]);
  RECON_AND_STORE(dest + 7 * stride, in[7]);
}

void av1_iht16x8_128_add_sse2(const tran_low_t *input, uint8_t *dest,
                              int stride, int tx_type) {
  __m128i in[16];

  // Transpose 16x8 input into in[]
  in[0] = load_input_data(input + 0 * 16);
  in[1] = load_input_data(input + 1 * 16);
  in[2] = load_input_data(input + 2 * 16);
  in[3] = load_input_data(input + 3 * 16);
  in[4] = load_input_data(input + 4 * 16);
  in[5] = load_input_data(input + 5 * 16);
  in[6] = load_input_data(input + 6 * 16);
  in[7] = load_input_data(input + 7 * 16);
  array_transpose_8x8(in, in);

  in[8] = load_input_data(input + 8 + 0 * 16);
  in[9] = load_input_data(input + 8 + 1 * 16);
  in[10] = load_input_data(input + 8 + 2 * 16);
  in[11] = load_input_data(input + 8 + 3 * 16);
  in[12] = load_input_data(input + 8 + 4 * 16);
  in[13] = load_input_data(input + 8 + 5 * 16);
  in[14] = load_input_data(input + 8 + 6 * 16);
  in[15] = load_input_data(input + 8 + 7 * 16);
  array_transpose_8x8(in + 8, in + 8);

  // Row transform
  switch (tx_type) {
    case DCT_DCT:
    case ADST_DCT:
#if CONFIG_EXT_TX
    case FLIPADST_DCT:
    case H_DCT:
#endif
      idct16_8col(in);
      break;
    case DCT_ADST:
    case ADST_ADST:
#if CONFIG_EXT_TX
    case DCT_FLIPADST:
    case FLIPADST_FLIPADST:
    case ADST_FLIPADST:
    case FLIPADST_ADST:
    case H_ADST:
    case H_FLIPADST:
#endif
      iadst16_8col(in);
      break;
#if CONFIG_EXT_TX
    case V_FLIPADST:
    case V_ADST:
    case V_DCT:
    case IDTX: iidtx16_8col(in); break;
#endif
    default: assert(0); break;
  }

  // Scale
  scale_sqrt2_8x8(in);
  scale_sqrt2_8x8(in + 8);

  // Column transform
  switch (tx_type) {
    case DCT_DCT:
    case DCT_ADST:
#if CONFIG_EXT_TX
    case DCT_FLIPADST:
    case V_DCT:
#endif
      aom_idct8_sse2(in);
      aom_idct8_sse2(in + 8);
      break;
    case ADST_DCT:
    case ADST_ADST:
#if CONFIG_EXT_TX
    case FLIPADST_ADST:
    case ADST_FLIPADST:
    case FLIPADST_FLIPADST:
    case FLIPADST_DCT:
    case V_ADST:
    case V_FLIPADST:
#endif
      aom_iadst8_sse2(in);
      aom_iadst8_sse2(in + 8);
      break;
#if CONFIG_EXT_TX
    case H_DCT:
    case H_ADST:
    case H_FLIPADST:
    case IDTX:
      array_transpose_8x8(in, in);
      array_transpose_8x8(in + 8, in + 8);
      iidtx8_sse2(in);
      iidtx8_sse2(in + 8);
      break;
#endif
    default: assert(0); break;
  }

  switch (tx_type) {
    case DCT_DCT:
    case ADST_DCT:
    case DCT_ADST:
    case ADST_ADST:
#if CONFIG_EXT_TX
    case H_DCT:
    case H_ADST:
    case V_ADST:
    case V_DCT:
    case IDTX:
#endif
      write_buffer_8x8_round6(dest, in, stride);
      write_buffer_8x8_round6(dest + 8, in + 8, stride);
      break;
#if CONFIG_EXT_TX
    case FLIPADST_DCT:
    case FLIPADST_ADST:
    case V_FLIPADST:
      write_buffer_8x8_round6(dest + stride * 7, in, -stride);
      write_buffer_8x8_round6(dest + stride * 7 + 8, in + 8, -stride);
      break;
    case DCT_FLIPADST:
    case ADST_FLIPADST:
    case H_FLIPADST:
      flip_buffer_lr_8x8(in);
      flip_buffer_lr_8x8(in + 8);
      write_buffer_8x8_round6(dest, in + 8, stride);
      write_buffer_8x8_round6(dest + 8, in, stride);
      break;
    case FLIPADST_FLIPADST:
      flip_buffer_lr_8x8(in);
      flip_buffer_lr_8x8(in + 8);
      write_buffer_8x8_round6(dest + stride * 7, in + 8, -stride);
      write_buffer_8x8_round6(dest + stride * 7 + 8, in, -stride);
      break;
#endif
    default: assert(0); break;
  }
}

static INLINE void write_buffer_8x4_round5(uint8_t *dest, __m128i *in,
                                           int stride) {
  const __m128i final_rounding = _mm_set1_epi16(1 << 4);
  const __m128i zero = _mm_setzero_si128();
  // Final rounding and shift
  in[0] = _mm_adds_epi16(in[0], final_rounding);
  in[1] = _mm_adds_epi16(in[1], final_rounding);
  in[2] = _mm_adds_epi16(in[2], final_rounding);
  in[3] = _mm_adds_epi16(in[3], final_rounding);

  in[0] = _mm_srai_epi16(in[0], 5);
  in[1] = _mm_srai_epi16(in[1], 5);
  in[2] = _mm_srai_epi16(in[2], 5);
  in[3] = _mm_srai_epi16(in[3], 5);

  RECON_AND_STORE(dest + 0 * stride, in[0]);
  RECON_AND_STORE(dest + 1 * stride, in[1]);
  RECON_AND_STORE(dest + 2 * stride, in[2]);
  RECON_AND_STORE(dest + 3 * stride, in[3]);
}

void av1_iht8x4_32_add_sse2(const tran_low_t *input, uint8_t *dest, int stride,
                            int tx_type) {
  __m128i in[8];

  in[0] = load_input_data(input + 0 * 8);
  in[1] = load_input_data(input + 1 * 8);
  in[2] = load_input_data(input + 2 * 8);
  in[3] = load_input_data(input + 3 * 8);

  // Row transform
  switch (tx_type) {
    case DCT_DCT:
    case ADST_DCT:
#if CONFIG_EXT_TX
    case FLIPADST_DCT:
    case H_DCT:
#endif
      aom_idct8_sse2(in);
      break;
    case DCT_ADST:
    case ADST_ADST:
#if CONFIG_EXT_TX
    case DCT_FLIPADST:
    case FLIPADST_FLIPADST:
    case ADST_FLIPADST:
    case FLIPADST_ADST:
    case H_ADST:
    case H_FLIPADST: aom_iadst8_sse2(in); break;
    case V_FLIPADST:
    case V_ADST:
    case V_DCT:
    case IDTX: iidtx8_sse2(in); array_transpose_8x8(in, in);
#endif
      break;
    default: assert(0); break;
  }

  scale_sqrt2_8x8(in);

  // Repack data. We pack into the bottom half of 'in'
  // so that the next repacking stage can pack into the
  // top half without overwriting anything
  in[7] = _mm_unpacklo_epi64(in[6], in[7]);
  in[6] = _mm_unpacklo_epi64(in[4], in[5]);
  in[5] = _mm_unpacklo_epi64(in[2], in[3]);
  in[4] = _mm_unpacklo_epi64(in[0], in[1]);

  // Column transform
  switch (tx_type) {
    case DCT_DCT:
    case DCT_ADST:
#if CONFIG_EXT_TX
    case DCT_FLIPADST:
    case V_DCT:
#endif
      aom_idct4_sse2(in + 4);
      aom_idct4_sse2(in + 6);
      break;
    case ADST_DCT:
    case ADST_ADST:
#if CONFIG_EXT_TX
    case FLIPADST_ADST:
    case ADST_FLIPADST:
    case FLIPADST_FLIPADST:
    case FLIPADST_DCT:
    case V_ADST:
    case V_FLIPADST:
#endif
      aom_iadst4_sse2(in + 4);
      aom_iadst4_sse2(in + 6);
      break;
#if CONFIG_EXT_TX
    case H_DCT:
    case H_ADST:
    case H_FLIPADST:
    case IDTX:
      iidtx4_sse2(in + 4);
      array_transpose_4x4(in + 4);
      iidtx4_sse2(in + 6);
      array_transpose_4x4(in + 6);
      break;
#endif
    default: assert(0); break;
  }

  // Repack data
  in[0] = _mm_unpacklo_epi64(in[4], in[6]);
  in[1] = _mm_unpackhi_epi64(in[4], in[6]);
  in[2] = _mm_unpacklo_epi64(in[5], in[7]);
  in[3] = _mm_unpackhi_epi64(in[5], in[7]);

  switch (tx_type) {
    case DCT_DCT:
    case ADST_DCT:
    case DCT_ADST:
    case ADST_ADST:
#if CONFIG_EXT_TX
    case H_DCT:
    case H_ADST:
    case V_ADST:
    case V_DCT:
    case IDTX: break;
    case FLIPADST_DCT:
    case FLIPADST_ADST:
    case V_FLIPADST: FLIPUD_PTR(dest, stride, 4); break;
    case DCT_FLIPADST:
    case ADST_FLIPADST:
    case H_FLIPADST:
      in[0] = mm_reverse_epi16(in[0]);
      in[1] = mm_reverse_epi16(in[1]);
      in[2] = mm_reverse_epi16(in[2]);
      in[3] = mm_reverse_epi16(in[3]);
      break;
    case FLIPADST_FLIPADST:
      in[0] = mm_reverse_epi16(in[0]);
      in[1] = mm_reverse_epi16(in[1]);
      in[2] = mm_reverse_epi16(in[2]);
      in[3] = mm_reverse_epi16(in[3]);
      FLIPUD_PTR(dest, stride, 4);
#endif
      break;
    default: assert(0); break;
  }
  write_buffer_8x4_round5(dest, in, stride);
}

static INLINE void write_buffer_4x8_round5(uint8_t *dest, __m128i *in,
                                           int stride) {
  const __m128i final_rounding = _mm_set1_epi16(1 << 4);
  const __m128i zero = _mm_setzero_si128();
  // Final rounding and shift
  in[0] = _mm_adds_epi16(in[0], final_rounding);
  in[1] = _mm_adds_epi16(in[1], final_rounding);
  in[2] = _mm_adds_epi16(in[2], final_rounding);
  in[3] = _mm_adds_epi16(in[3], final_rounding);

  in[0] = _mm_srai_epi16(in[0], 5);
  in[1] = _mm_srai_epi16(in[1], 5);
  in[2] = _mm_srai_epi16(in[2], 5);
  in[3] = _mm_srai_epi16(in[3], 5);

  // Reconstruction and Store
  {
    __m128i d0 = _mm_cvtsi32_si128(*(const int *)(dest + stride * 0));
    __m128i d1 = _mm_cvtsi32_si128(*(const int *)(dest + stride * 1));
    __m128i d2 = _mm_cvtsi32_si128(*(const int *)(dest + stride * 2));
    __m128i d3 = _mm_cvtsi32_si128(*(const int *)(dest + stride * 3));
    __m128i d4 = _mm_cvtsi32_si128(*(const int *)(dest + stride * 4));
    __m128i d5 = _mm_cvtsi32_si128(*(const int *)(dest + stride * 5));
    __m128i d6 = _mm_cvtsi32_si128(*(const int *)(dest + stride * 6));
    __m128i d7 = _mm_cvtsi32_si128(*(const int *)(dest + stride * 7));

    d0 = _mm_unpacklo_epi32(d0, d1);
    d2 = _mm_unpacklo_epi32(d2, d3);
    d4 = _mm_unpacklo_epi32(d4, d5);
    d6 = _mm_unpacklo_epi32(d6, d7);
    d0 = _mm_unpacklo_epi8(d0, zero);
    d2 = _mm_unpacklo_epi8(d2, zero);
    d4 = _mm_unpacklo_epi8(d4, zero);
    d6 = _mm_unpacklo_epi8(d6, zero);
    d0 = _mm_add_epi16(d0, in[0]);
    d2 = _mm_add_epi16(d2, in[1]);
    d4 = _mm_add_epi16(d4, in[2]);
    d6 = _mm_add_epi16(d6, in[3]);

    d0 = _mm_packus_epi16(d0, d2);
    *(int *)dest = _mm_cvtsi128_si32(d0);
    d0 = _mm_srli_si128(d0, 4);
    *(int *)(dest + stride) = _mm_cvtsi128_si32(d0);
    d0 = _mm_srli_si128(d0, 4);
    *(int *)(dest + stride * 2) = _mm_cvtsi128_si32(d0);
    d0 = _mm_srli_si128(d0, 4);
    *(int *)(dest + stride * 3) = _mm_cvtsi128_si32(d0);
    d0 = _mm_packus_epi16(d4, d6);
    *(int *)(dest + stride * 4) = _mm_cvtsi128_si32(d0);
    d0 = _mm_srli_si128(d0, 4);
    *(int *)(dest + stride * 5) = _mm_cvtsi128_si32(d0);
    d0 = _mm_srli_si128(d0, 4);
    *(int *)(dest + stride * 6) = _mm_cvtsi128_si32(d0);
    d0 = _mm_srli_si128(d0, 4);
    *(int *)(dest + stride * 7) = _mm_cvtsi128_si32(d0);
  }
}

void av1_iht4x8_32_add_sse2(const tran_low_t *input, uint8_t *dest, int stride,
                            int tx_type) {
  __m128i in[8];

  // Load rows, packed two per element of 'in'.
  // We pack into the bottom half of 'in' so that the
  // later repacking stage can pack into the
  // top half without overwriting anything
  in[4] = load_input_data(input + 0 * 8);
  in[5] = load_input_data(input + 1 * 8);
  in[6] = load_input_data(input + 2 * 8);
  in[7] = load_input_data(input + 3 * 8);

  scale_sqrt2_8x4(in + 4);

  // Row transform
  switch (tx_type) {
    case DCT_DCT:
    case ADST_DCT:
#if CONFIG_EXT_TX
    case FLIPADST_DCT:
    case H_DCT:
#endif
      aom_idct4_sse2(in + 4);
      aom_idct4_sse2(in + 6);
      break;
    case DCT_ADST:
    case ADST_ADST:
#if CONFIG_EXT_TX
    case DCT_FLIPADST:
    case FLIPADST_FLIPADST:
    case ADST_FLIPADST:
    case FLIPADST_ADST:
    case H_ADST:
    case H_FLIPADST:
#endif
      aom_iadst4_sse2(in + 4);
      aom_iadst4_sse2(in + 6);
      break;
#if CONFIG_EXT_TX
    case V_FLIPADST:
    case V_ADST:
    case V_DCT:
    case IDTX:
      iidtx4_sse2(in + 4);
      array_transpose_4x4(in + 4);
      iidtx4_sse2(in + 6);
      array_transpose_4x4(in + 6);
      break;
#endif
    default: assert(0); break;
  }

  // Repack data
  in[0] = _mm_unpacklo_epi64(in[4], in[6]);
  in[1] = _mm_unpackhi_epi64(in[4], in[6]);
  in[2] = _mm_unpacklo_epi64(in[5], in[7]);
  in[3] = _mm_unpackhi_epi64(in[5], in[7]);

  // Column transform
  switch (tx_type) {
    case DCT_DCT:
    case DCT_ADST:
#if CONFIG_EXT_TX
    case DCT_FLIPADST:
    case V_DCT:
#endif
      aom_idct8_sse2(in);
      break;
    case ADST_DCT:
    case ADST_ADST:
#if CONFIG_EXT_TX
    case FLIPADST_ADST:
    case ADST_FLIPADST:
    case FLIPADST_FLIPADST:
    case FLIPADST_DCT:
    case V_ADST:
    case V_FLIPADST:
#endif
      aom_iadst8_sse2(in);
      break;
#if CONFIG_EXT_TX
    case H_DCT:
    case H_ADST:
    case H_FLIPADST:
    case IDTX:
      iidtx8_sse2(in);
      array_transpose_8x8(in, in);
      break;
#endif
    default: assert(0); break;
  }

  switch (tx_type) {
    case DCT_DCT:
    case ADST_DCT:
    case DCT_ADST:
    case ADST_ADST:
#if CONFIG_EXT_TX
    case H_DCT:
    case H_ADST:
    case V_ADST:
    case V_DCT:
    case IDTX:
#endif
      break;
#if CONFIG_EXT_TX
    case FLIPADST_DCT:
    case FLIPADST_ADST:
    case V_FLIPADST: FLIPUD_PTR(dest, stride, 8); break;
    case DCT_FLIPADST:
    case ADST_FLIPADST:
    case H_FLIPADST:
      in[0] = _mm_shufflelo_epi16(in[0], 0x1b);
      in[1] = _mm_shufflelo_epi16(in[1], 0x1b);
      in[2] = _mm_shufflelo_epi16(in[2], 0x1b);
      in[3] = _mm_shufflelo_epi16(in[3], 0x1b);
      in[4] = _mm_shufflelo_epi16(in[4], 0x1b);
      in[5] = _mm_shufflelo_epi16(in[5], 0x1b);
      in[6] = _mm_shufflelo_epi16(in[6], 0x1b);
      in[7] = _mm_shufflelo_epi16(in[7], 0x1b);
      break;
    case FLIPADST_FLIPADST:
      in[0] = _mm_shufflelo_epi16(in[0], 0x1b);
      in[1] = _mm_shufflelo_epi16(in[1], 0x1b);
      in[2] = _mm_shufflelo_epi16(in[2], 0x1b);
      in[3] = _mm_shufflelo_epi16(in[3], 0x1b);
      in[4] = _mm_shufflelo_epi16(in[4], 0x1b);
      in[5] = _mm_shufflelo_epi16(in[5], 0x1b);
      in[6] = _mm_shufflelo_epi16(in[6], 0x1b);
      in[7] = _mm_shufflelo_epi16(in[7], 0x1b);
      FLIPUD_PTR(dest, stride, 8);
      break;
#endif
    default: assert(0); break;
  }
  in[0] = _mm_unpacklo_epi64(in[0], in[1]);
  in[1] = _mm_unpacklo_epi64(in[2], in[3]);
  in[2] = _mm_unpacklo_epi64(in[4], in[5]);
  in[3] = _mm_unpacklo_epi64(in[6], in[7]);
  write_buffer_4x8_round5(dest, in, stride);
}

// Note: The 16-column 32-element transforms take input in the form of four
// 8x16 blocks (each stored as a __m128i[16]), which are the four quadrants
// of the overall 16x32 input buffer.
static INLINE void idct32_16col(__m128i *tl, __m128i *tr, __m128i *bl,
                                __m128i *br) {
  array_transpose_16x16(tl, tr);
  array_transpose_16x16(bl, br);
  idct32_8col(tl, bl);
  idct32_8col(tr, br);
}

static INLINE void ihalfright32_16col(__m128i *tl, __m128i *tr, __m128i *bl,
                                      __m128i *br) {
  __m128i tmpl[16], tmpr[16];
  int i;

  // Copy the top half of the input to temporary storage
  for (i = 0; i < 16; ++i) {
    tmpl[i] = tl[i];
    tmpr[i] = tr[i];
  }

  // Generate the top half of the output
  for (i = 0; i < 16; ++i) {
    tl[i] = _mm_slli_epi16(bl[i], 2);
    tr[i] = _mm_slli_epi16(br[i], 2);
  }
  array_transpose_16x16(tl, tr);

  // Copy the temporary storage back to the bottom half of the input
  for (i = 0; i < 16; ++i) {
    bl[i] = tmpl[i];
    br[i] = tmpr[i];
  }

  // Generate the bottom half of the output
  scale_sqrt2_8x16(bl);
  scale_sqrt2_8x16(br);
  aom_idct16_sse2(bl, br);  // Includes a transposition
}

#if CONFIG_EXT_TX
static INLINE void iidtx32_16col(__m128i *tl, __m128i *tr, __m128i *bl,
                                 __m128i *br) {
  int i;
  array_transpose_16x16(tl, tr);
  array_transpose_16x16(bl, br);
  for (i = 0; i < 16; ++i) {
    tl[i] = _mm_slli_epi16(tl[i], 2);
    tr[i] = _mm_slli_epi16(tr[i], 2);
    bl[i] = _mm_slli_epi16(bl[i], 2);
    br[i] = _mm_slli_epi16(br[i], 2);
  }
}
#endif  // CONFIG_EXT_TX

static INLINE void write_buffer_16x32_round6(uint8_t *dest, __m128i *intl,
                                             __m128i *intr, __m128i *inbl,
                                             __m128i *inbr, int stride) {
  const __m128i zero = _mm_setzero_si128();
  const __m128i final_rounding = _mm_set1_epi16(1 << 5);
  int i;

  for (i = 0; i < 16; ++i) {
    intl[i] = _mm_adds_epi16(intl[i], final_rounding);
    intr[i] = _mm_adds_epi16(intr[i], final_rounding);
    inbl[i] = _mm_adds_epi16(inbl[i], final_rounding);
    inbr[i] = _mm_adds_epi16(inbr[i], final_rounding);
    intl[i] = _mm_srai_epi16(intl[i], 6);
    intr[i] = _mm_srai_epi16(intr[i], 6);
    inbl[i] = _mm_srai_epi16(inbl[i], 6);
    inbr[i] = _mm_srai_epi16(inbr[i], 6);
    RECON_AND_STORE(dest + i * stride + 0, intl[i]);
    RECON_AND_STORE(dest + i * stride + 8, intr[i]);
    RECON_AND_STORE(dest + (i + 16) * stride + 0, inbl[i]);
    RECON_AND_STORE(dest + (i + 16) * stride + 8, inbr[i]);
  }
}

void av1_iht16x32_512_add_sse2(const tran_low_t *input, uint8_t *dest,
                               int stride, int tx_type) {
  __m128i intl[16], intr[16], inbl[16], inbr[16];

  int i;
  for (i = 0; i < 16; ++i) {
    intl[i] = load_input_data(input + i * 16 + 0);
    intr[i] = load_input_data(input + i * 16 + 8);
    inbl[i] = load_input_data(input + (i + 16) * 16 + 0);
    inbr[i] = load_input_data(input + (i + 16) * 16 + 8);
  }

  // Row transform
  switch (tx_type) {
    case DCT_DCT:
    case ADST_DCT:
#if CONFIG_EXT_TX
    case FLIPADST_DCT:
    case H_DCT:
#endif
      aom_idct16_sse2(intl, intr);
      aom_idct16_sse2(inbl, inbr);
      break;
    case DCT_ADST:
    case ADST_ADST:
#if CONFIG_EXT_TX
    case DCT_FLIPADST:
    case FLIPADST_FLIPADST:
    case ADST_FLIPADST:
    case FLIPADST_ADST:
    case H_ADST:
    case H_FLIPADST:
#endif
      aom_iadst16_sse2(intl, intr);
      aom_iadst16_sse2(inbl, inbr);
      break;
#if CONFIG_EXT_TX
    case V_FLIPADST:
    case V_ADST:
    case V_DCT:
    case IDTX:
      iidtx16_sse2(intl, intr);
      iidtx16_sse2(inbl, inbr);
      break;
#endif
    default: assert(0); break;
  }

  scale_sqrt2_8x16(intl);
  scale_sqrt2_8x16(intr);
  scale_sqrt2_8x16(inbl);
  scale_sqrt2_8x16(inbr);

  // Column transform
  switch (tx_type) {
    case DCT_DCT:
    case DCT_ADST:
#if CONFIG_EXT_TX
    case DCT_FLIPADST:
    case V_DCT:
#endif
      idct32_16col(intl, intr, inbl, inbr);
      break;
    case ADST_DCT:
    case ADST_ADST:
#if CONFIG_EXT_TX
    case FLIPADST_ADST:
    case ADST_FLIPADST:
    case FLIPADST_FLIPADST:
    case FLIPADST_DCT:
    case V_ADST:
    case V_FLIPADST:
#endif
      ihalfright32_16col(intl, intr, inbl, inbr);
      break;
#if CONFIG_EXT_TX
    case H_DCT:
    case H_ADST:
    case H_FLIPADST:
    case IDTX: iidtx32_16col(intl, intr, inbl, inbr); break;
#endif
    default: assert(0); break;
  }

  switch (tx_type) {
    case DCT_DCT:
    case ADST_DCT:
    case DCT_ADST:
    case ADST_ADST:
#if CONFIG_EXT_TX
    case H_DCT:
    case H_ADST:
    case V_ADST:
    case V_DCT:
    case IDTX:
#endif
      break;
#if CONFIG_EXT_TX
    case FLIPADST_DCT:
    case FLIPADST_ADST:
    case V_FLIPADST: FLIPUD_PTR(dest, stride, 32); break;
    case DCT_FLIPADST:
    case ADST_FLIPADST:
    case H_FLIPADST:
      for (i = 0; i < 16; ++i) {
        __m128i tmp = intl[i];
        intl[i] = mm_reverse_epi16(intr[i]);
        intr[i] = mm_reverse_epi16(tmp);
        tmp = inbl[i];
        inbl[i] = mm_reverse_epi16(inbr[i]);
        inbr[i] = mm_reverse_epi16(tmp);
      }
      break;
    case FLIPADST_FLIPADST:
      for (i = 0; i < 16; ++i) {
        __m128i tmp = intl[i];
        intl[i] = mm_reverse_epi16(intr[i]);
        intr[i] = mm_reverse_epi16(tmp);
        tmp = inbl[i];
        inbl[i] = mm_reverse_epi16(inbr[i]);
        inbr[i] = mm_reverse_epi16(tmp);
      }
      FLIPUD_PTR(dest, stride, 32);
      break;
#endif
    default: assert(0); break;
  }
  write_buffer_16x32_round6(dest, intl, intr, inbl, inbr, stride);
}

static INLINE void write_buffer_32x16_round6(uint8_t *dest, __m128i *in0,
                                             __m128i *in1, __m128i *in2,
                                             __m128i *in3, int stride) {
  const __m128i zero = _mm_setzero_si128();
  const __m128i final_rounding = _mm_set1_epi16(1 << 5);
  int i;

  for (i = 0; i < 16; ++i) {
    in0[i] = _mm_adds_epi16(in0[i], final_rounding);
    in1[i] = _mm_adds_epi16(in1[i], final_rounding);
    in2[i] = _mm_adds_epi16(in2[i], final_rounding);
    in3[i] = _mm_adds_epi16(in3[i], final_rounding);
    in0[i] = _mm_srai_epi16(in0[i], 6);
    in1[i] = _mm_srai_epi16(in1[i], 6);
    in2[i] = _mm_srai_epi16(in2[i], 6);
    in3[i] = _mm_srai_epi16(in3[i], 6);
    RECON_AND_STORE(dest + i * stride + 0, in0[i]);
    RECON_AND_STORE(dest + i * stride + 8, in1[i]);
    RECON_AND_STORE(dest + i * stride + 16, in2[i]);
    RECON_AND_STORE(dest + i * stride + 24, in3[i]);
  }
}

void av1_iht32x16_512_add_sse2(const tran_low_t *input, uint8_t *dest,
                               int stride, int tx_type) {
  __m128i in0[16], in1[16], in2[16], in3[16];
  int i;

  for (i = 0; i < 16; ++i) {
    in0[i] = load_input_data(input + i * 32 + 0);
    in1[i] = load_input_data(input + i * 32 + 8);
    in2[i] = load_input_data(input + i * 32 + 16);
    in3[i] = load_input_data(input + i * 32 + 24);
  }

  // Row transform
  switch (tx_type) {
    case DCT_DCT:
    case ADST_DCT:
#if CONFIG_EXT_TX
    case FLIPADST_DCT:
    case H_DCT:
#endif
      idct32_16col(in0, in1, in2, in3);
      break;
    case DCT_ADST:
    case ADST_ADST:
#if CONFIG_EXT_TX
    case DCT_FLIPADST:
    case FLIPADST_FLIPADST:
    case ADST_FLIPADST:
    case FLIPADST_ADST:
    case H_ADST:
    case H_FLIPADST:
#endif
      ihalfright32_16col(in0, in1, in2, in3);
      break;
#if CONFIG_EXT_TX
    case V_FLIPADST:
    case V_ADST:
    case V_DCT:
    case IDTX: iidtx32_16col(in0, in1, in2, in3); break;
#endif
    default: assert(0); break;
  }

  scale_sqrt2_8x16(in0);
  scale_sqrt2_8x16(in1);
  scale_sqrt2_8x16(in2);
  scale_sqrt2_8x16(in3);

  // Column transform
  switch (tx_type) {
    case DCT_DCT:
    case DCT_ADST:
#if CONFIG_EXT_TX
    case DCT_FLIPADST:
    case V_DCT:
#endif
      aom_idct16_sse2(in0, in1);
      aom_idct16_sse2(in2, in3);
      break;
    case ADST_DCT:
    case ADST_ADST:
#if CONFIG_EXT_TX
    case FLIPADST_ADST:
    case ADST_FLIPADST:
    case FLIPADST_FLIPADST:
    case FLIPADST_DCT:
    case V_ADST:
    case V_FLIPADST:
#endif
      aom_iadst16_sse2(in0, in1);
      aom_iadst16_sse2(in2, in3);
      break;
#if CONFIG_EXT_TX
    case H_DCT:
    case H_ADST:
    case H_FLIPADST:
    case IDTX:
      iidtx16_sse2(in0, in1);
      iidtx16_sse2(in2, in3);
      break;
#endif
    default: assert(0); break;
  }

  switch (tx_type) {
    case DCT_DCT:
    case ADST_DCT:
    case DCT_ADST:
    case ADST_ADST:
#if CONFIG_EXT_TX
    case H_DCT:
    case H_ADST:
    case V_ADST:
    case V_DCT:
    case IDTX:
#endif
      break;
#if CONFIG_EXT_TX
    case FLIPADST_DCT:
    case FLIPADST_ADST:
    case V_FLIPADST: FLIPUD_PTR(dest, stride, 16); break;
    case DCT_FLIPADST:
    case ADST_FLIPADST:
    case H_FLIPADST:
      for (i = 0; i < 16; ++i) {
        __m128i tmp1 = in0[i];
        __m128i tmp2 = in1[i];
        in0[i] = mm_reverse_epi16(in3[i]);
        in1[i] = mm_reverse_epi16(in2[i]);
        in2[i] = mm_reverse_epi16(tmp2);
        in3[i] = mm_reverse_epi16(tmp1);
      }
      break;
    case FLIPADST_FLIPADST:
      for (i = 0; i < 16; ++i) {
        __m128i tmp1 = in0[i];
        __m128i tmp2 = in1[i];
        in0[i] = mm_reverse_epi16(in3[i]);
        in1[i] = mm_reverse_epi16(in2[i]);
        in2[i] = mm_reverse_epi16(tmp2);
        in3[i] = mm_reverse_epi16(tmp1);
      }
      FLIPUD_PTR(dest, stride, 16);
      break;
#endif
    default: assert(0); break;
  }
  write_buffer_32x16_round6(dest, in0, in1, in2, in3, stride);
}
