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

void vp9_short_idct4x4llm_sse2(int16_t *input, int16_t *output, int pitch) {
  const __m128i zero = _mm_setzero_si128();
  const __m128i eight = _mm_set1_epi16(8);
  const __m128i cst = _mm_setr_epi16((short)cospi_16_64, (short)cospi_16_64,
                                     (short)cospi_16_64, (short)-cospi_16_64,
                                     (short)cospi_24_64, (short)-cospi_8_64,
                                     (short)cospi_8_64, (short)cospi_24_64);
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
#endif
