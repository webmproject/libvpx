/*
 *  Copyright (c) 2012 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include "./vp10_rtcd.h"
#include "vpx_dsp/x86/inv_txfm_sse2.h"
#include "vpx_dsp/x86/txfm_common_sse2.h"
#include "vpx_ports/mem.h"
#include "vp10/common/enums.h"

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

#define FLIPLR_16x16(in0, in1) do {             \
  __m128i *tmp;                                 \
  fliplr_16x8(in0);                             \
  fliplr_16x8(in1);                             \
  tmp = (in0);                                  \
  (in0) = (in1);                                \
  (in1) = tmp;                                  \
} while (0)

#define FLIPUD_PTR(dest, stride, size) do {     \
    (dest) = (dest) + ((size) - 1) * (stride);  \
    (stride) = - (stride);                      \
} while (0)
#endif

void vp10_iht4x4_16_add_sse2(const tran_low_t *input, uint8_t *dest, int stride,
                             int tx_type) {
  __m128i in[2];
  const __m128i zero = _mm_setzero_si128();
  const __m128i eight = _mm_set1_epi16(8);

  in[0] = load_input_data(input);
  in[1] = load_input_data(input + 8);

  switch (tx_type) {
    case DCT_DCT:
      idct4_sse2(in);
      idct4_sse2(in);
      break;
    case ADST_DCT:
      idct4_sse2(in);
      iadst4_sse2(in);
      break;
    case DCT_ADST:
      iadst4_sse2(in);
      idct4_sse2(in);
      break;
    case ADST_ADST:
      iadst4_sse2(in);
      iadst4_sse2(in);
      break;
#if CONFIG_EXT_TX
    case FLIPADST_DCT:
      idct4_sse2(in);
      iadst4_sse2(in);
      FLIPUD_PTR(dest, stride, 4);
      break;
    case DCT_FLIPADST:
      iadst4_sse2(in);
      idct4_sse2(in);
      fliplr_4x4(in);
      break;
    case FLIPADST_FLIPADST:
      iadst4_sse2(in);
      iadst4_sse2(in);
      FLIPUD_PTR(dest, stride, 4);
      fliplr_4x4(in);
      break;
    case ADST_FLIPADST:
      iadst4_sse2(in);
      iadst4_sse2(in);
      fliplr_4x4(in);
      break;
    case FLIPADST_ADST:
      iadst4_sse2(in);
      iadst4_sse2(in);
      FLIPUD_PTR(dest, stride, 4);
      break;
#endif  // CONFIG_EXT_TX
    default:
      assert(0);
      break;
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

void vp10_iht8x8_64_add_sse2(const tran_low_t *input, uint8_t *dest, int stride,
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
      idct8_sse2(in);
      idct8_sse2(in);
      break;
    case ADST_DCT:
      idct8_sse2(in);
      iadst8_sse2(in);
      break;
    case DCT_ADST:
      iadst8_sse2(in);
      idct8_sse2(in);
      break;
    case ADST_ADST:
      iadst8_sse2(in);
      iadst8_sse2(in);
      break;
#if CONFIG_EXT_TX
    case FLIPADST_DCT:
      idct8_sse2(in);
      iadst8_sse2(in);
      FLIPUD_PTR(dest, stride, 8);
      break;
    case DCT_FLIPADST:
      iadst8_sse2(in);
      idct8_sse2(in);
      fliplr_8x8(in);
      break;
    case FLIPADST_FLIPADST:
      iadst8_sse2(in);
      iadst8_sse2(in);
      FLIPUD_PTR(dest, stride, 8);
      fliplr_8x8(in);
      break;
    case ADST_FLIPADST:
      iadst8_sse2(in);
      iadst8_sse2(in);
      fliplr_8x8(in);
      break;
    case FLIPADST_ADST:
      iadst8_sse2(in);
      iadst8_sse2(in);
      FLIPUD_PTR(dest, stride, 8);
      break;
#endif  // CONFIG_EXT_TX
    default:
      assert(0);
      break;
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

void vp10_iht16x16_256_add_sse2(const tran_low_t *input, uint8_t *dest,
                                int stride, int tx_type) {
  __m128i in[32];
  __m128i *in0 = &in[0];
  __m128i *in1 = &in[16];

  load_buffer_8x16(input, in0);
  input += 8;
  load_buffer_8x16(input, in1);

  switch (tx_type) {
    case DCT_DCT:
      idct16_sse2(in0, in1);
      idct16_sse2(in0, in1);
      break;
    case ADST_DCT:
      idct16_sse2(in0, in1);
      iadst16_sse2(in0, in1);
      break;
    case DCT_ADST:
      iadst16_sse2(in0, in1);
      idct16_sse2(in0, in1);
      break;
    case ADST_ADST:
      iadst16_sse2(in0, in1);
      iadst16_sse2(in0, in1);
      break;
#if CONFIG_EXT_TX
    case FLIPADST_DCT:
      idct16_sse2(in0, in1);
      iadst16_sse2(in0, in1);
      FLIPUD_PTR(dest, stride, 16);
      break;
    case DCT_FLIPADST:
      iadst16_sse2(in0, in1);
      idct16_sse2(in0, in1);
      FLIPLR_16x16(in0, in1);
      break;
    case FLIPADST_FLIPADST:
      iadst16_sse2(in0, in1);
      iadst16_sse2(in0, in1);
      FLIPUD_PTR(dest, stride, 16);
      FLIPLR_16x16(in0, in1);
      break;
    case ADST_FLIPADST:
      iadst16_sse2(in0, in1);
      iadst16_sse2(in0, in1);
      FLIPLR_16x16(in0, in1);
      break;
    case FLIPADST_ADST:
      iadst16_sse2(in0, in1);
      iadst16_sse2(in0, in1);
      FLIPUD_PTR(dest, stride, 16);
      break;
#endif  // CONFIG_EXT_TX
    default:
      assert(0);
      break;
  }

  write_buffer_8x16(dest, in0, stride);
  dest += 8;
  write_buffer_8x16(dest, in1, stride);
}
