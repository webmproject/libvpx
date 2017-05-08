/*
 *  Copyright (c) 2015 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include "./vpx_dsp_rtcd.h"
#include "vpx_dsp/x86/highbd_inv_txfm_sse2.h"
#include "vpx_dsp/x86/inv_txfm_sse2.h"
#include "vpx_dsp/x86/transpose_sse2.h"
#include "vpx_dsp/x86/txfm_common_sse2.h"

void vpx_highbd_idct4x4_16_add_sse2(const tran_low_t *input, uint16_t *dest,
                                    int stride, int bd) {
  tran_low_t out[4 * 4];
  tran_low_t *outptr = out;
  int i, j;
  __m128i inptr[4];
  __m128i sign_bits[2];
  __m128i temp_mm, min_input, max_input;
  int test;
  int optimised_cols = 0;
  const __m128i zero = _mm_set1_epi16(0);
  const __m128i eight = _mm_set1_epi16(8);
  const __m128i max = _mm_set1_epi16(12043);
  const __m128i min = _mm_set1_epi16(-12043);
  // Load input into __m128i
  inptr[0] = _mm_loadu_si128((const __m128i *)input);
  inptr[1] = _mm_loadu_si128((const __m128i *)(input + 4));
  inptr[2] = _mm_loadu_si128((const __m128i *)(input + 8));
  inptr[3] = _mm_loadu_si128((const __m128i *)(input + 12));

  // Pack to 16 bits
  inptr[0] = _mm_packs_epi32(inptr[0], inptr[1]);
  inptr[1] = _mm_packs_epi32(inptr[2], inptr[3]);

  max_input = _mm_max_epi16(inptr[0], inptr[1]);
  min_input = _mm_min_epi16(inptr[0], inptr[1]);
  max_input = _mm_cmpgt_epi16(max_input, max);
  min_input = _mm_cmplt_epi16(min_input, min);
  temp_mm = _mm_or_si128(max_input, min_input);
  test = _mm_movemask_epi8(temp_mm);

  if (!test) {
    // Do the row transform
    idct4_sse2(inptr);

    // Check the min & max values
    max_input = _mm_max_epi16(inptr[0], inptr[1]);
    min_input = _mm_min_epi16(inptr[0], inptr[1]);
    max_input = _mm_cmpgt_epi16(max_input, max);
    min_input = _mm_cmplt_epi16(min_input, min);
    temp_mm = _mm_or_si128(max_input, min_input);
    test = _mm_movemask_epi8(temp_mm);

    if (test) {
      transpose_4x4(inptr);
      sign_bits[0] = _mm_cmplt_epi16(inptr[0], zero);
      sign_bits[1] = _mm_cmplt_epi16(inptr[1], zero);
      inptr[3] = _mm_unpackhi_epi16(inptr[1], sign_bits[1]);
      inptr[2] = _mm_unpacklo_epi16(inptr[1], sign_bits[1]);
      inptr[1] = _mm_unpackhi_epi16(inptr[0], sign_bits[0]);
      inptr[0] = _mm_unpacklo_epi16(inptr[0], sign_bits[0]);
      _mm_storeu_si128((__m128i *)outptr, inptr[0]);
      _mm_storeu_si128((__m128i *)(outptr + 4), inptr[1]);
      _mm_storeu_si128((__m128i *)(outptr + 8), inptr[2]);
      _mm_storeu_si128((__m128i *)(outptr + 12), inptr[3]);
    } else {
      // Set to use the optimised transform for the column
      optimised_cols = 1;
    }
  } else {
    // Run the un-optimised row transform
    for (i = 0; i < 4; ++i) {
      vpx_highbd_idct4_c(input, outptr, bd);
      input += 4;
      outptr += 4;
    }
  }

  if (optimised_cols) {
    idct4_sse2(inptr);

    // Final round and shift
    inptr[0] = _mm_add_epi16(inptr[0], eight);
    inptr[1] = _mm_add_epi16(inptr[1], eight);

    inptr[0] = _mm_srai_epi16(inptr[0], 4);
    inptr[1] = _mm_srai_epi16(inptr[1], 4);

    // Reconstruction and Store
    {
      __m128i d0 = _mm_loadl_epi64((const __m128i *)dest);
      __m128i d2 = _mm_loadl_epi64((const __m128i *)(dest + stride * 2));
      d0 = _mm_unpacklo_epi64(
          d0, _mm_loadl_epi64((const __m128i *)(dest + stride)));
      d2 = _mm_unpacklo_epi64(
          d2, _mm_loadl_epi64((const __m128i *)(dest + stride * 3)));
      d0 = clamp_high_sse2(_mm_adds_epi16(d0, inptr[0]), bd);
      d2 = clamp_high_sse2(_mm_adds_epi16(d2, inptr[1]), bd);
      // store input0
      _mm_storel_epi64((__m128i *)dest, d0);
      // store input1
      d0 = _mm_srli_si128(d0, 8);
      _mm_storel_epi64((__m128i *)(dest + stride), d0);
      // store input2
      _mm_storel_epi64((__m128i *)(dest + stride * 2), d2);
      // store input3
      d2 = _mm_srli_si128(d2, 8);
      _mm_storel_epi64((__m128i *)(dest + stride * 3), d2);
    }
  } else {
    // Run the un-optimised column transform
    tran_low_t temp_in[4], temp_out[4];
    // Columns
    for (i = 0; i < 4; ++i) {
      for (j = 0; j < 4; ++j) temp_in[j] = out[j * 4 + i];
      vpx_highbd_idct4_c(temp_in, temp_out, bd);
      for (j = 0; j < 4; ++j) {
        dest[j * stride + i] = highbd_clip_pixel_add(
            dest[j * stride + i], ROUND_POWER_OF_TWO(temp_out[j], 4), bd);
      }
    }
  }
}
