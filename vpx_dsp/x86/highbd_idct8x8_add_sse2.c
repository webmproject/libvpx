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

void vpx_highbd_idct8x8_64_add_sse2(const tran_low_t *input, uint16_t *dest,
                                    int stride, int bd) {
  tran_low_t out[8 * 8];
  tran_low_t *outptr = out;
  int i, j, test;
  __m128i inptr[8];
  __m128i min_input, max_input, temp1, temp2, sign_bits;
  const __m128i zero = _mm_set1_epi16(0);
  const __m128i sixteen = _mm_set1_epi16(16);
  const __m128i max = _mm_set1_epi16(6201);
  const __m128i min = _mm_set1_epi16(-6201);
  int optimised_cols = 0;

  // Load input into __m128i & pack to 16 bits
  for (i = 0; i < 8; i++) {
    temp1 = _mm_loadu_si128((const __m128i *)(input + 8 * i));
    temp2 = _mm_loadu_si128((const __m128i *)(input + 8 * i + 4));
    inptr[i] = _mm_packs_epi32(temp1, temp2);
  }

  // Find the min & max for the row transform
  max_input = _mm_max_epi16(inptr[0], inptr[1]);
  min_input = _mm_min_epi16(inptr[0], inptr[1]);
  for (i = 2; i < 8; i++) {
    max_input = _mm_max_epi16(max_input, inptr[i]);
    min_input = _mm_min_epi16(min_input, inptr[i]);
  }
  max_input = _mm_cmpgt_epi16(max_input, max);
  min_input = _mm_cmplt_epi16(min_input, min);
  temp1 = _mm_or_si128(max_input, min_input);
  test = _mm_movemask_epi8(temp1);

  if (!test) {
    // Do the row transform
    idct8_sse2(inptr);

    // Find the min & max for the column transform
    max_input = _mm_max_epi16(inptr[0], inptr[1]);
    min_input = _mm_min_epi16(inptr[0], inptr[1]);
    for (i = 2; i < 8; i++) {
      max_input = _mm_max_epi16(max_input, inptr[i]);
      min_input = _mm_min_epi16(min_input, inptr[i]);
    }
    max_input = _mm_cmpgt_epi16(max_input, max);
    min_input = _mm_cmplt_epi16(min_input, min);
    temp1 = _mm_or_si128(max_input, min_input);
    test = _mm_movemask_epi8(temp1);

    if (test) {
      array_transpose_8x8(inptr, inptr);
      for (i = 0; i < 8; i++) {
        sign_bits = _mm_cmplt_epi16(inptr[i], zero);
        temp1 = _mm_unpackhi_epi16(inptr[i], sign_bits);
        temp2 = _mm_unpacklo_epi16(inptr[i], sign_bits);
        _mm_storeu_si128((__m128i *)(outptr + 4 * (2 * i + 1)), temp1);
        _mm_storeu_si128((__m128i *)(outptr + 4 * (2 * i)), temp2);
      }
    } else {
      // Set to use the optimised transform for the column
      optimised_cols = 1;
    }
  } else {
    // Run the un-optimised row transform
    for (i = 0; i < 8; ++i) {
      vpx_highbd_idct8_c(input, outptr, bd);
      input += 8;
      outptr += 8;
    }
  }

  if (optimised_cols) {
    idct8_sse2(inptr);

    // Final round & shift and Reconstruction and Store
    {
      __m128i d[8];
      for (i = 0; i < 8; i++) {
        inptr[i] = _mm_add_epi16(inptr[i], sixteen);
        d[i] = _mm_loadu_si128((const __m128i *)(dest + stride * i));
        inptr[i] = _mm_srai_epi16(inptr[i], 5);
        d[i] = clamp_high_sse2(_mm_adds_epi16(d[i], inptr[i]), bd);
        // Store
        _mm_storeu_si128((__m128i *)(dest + stride * i), d[i]);
      }
    }
  } else {
    // Run the un-optimised column transform
    tran_low_t temp_in[8], temp_out[8];
    for (i = 0; i < 8; ++i) {
      for (j = 0; j < 8; ++j) temp_in[j] = out[j * 8 + i];
      vpx_highbd_idct8_c(temp_in, temp_out, bd);
      for (j = 0; j < 8; ++j) {
        dest[j * stride + i] = highbd_clip_pixel_add(
            dest[j * stride + i], ROUND_POWER_OF_TWO(temp_out[j], 5), bd);
      }
    }
  }
}

void vpx_highbd_idct8x8_12_add_sse2(const tran_low_t *input, uint16_t *dest,
                                    int stride, int bd) {
  tran_low_t out[8 * 8] = { 0 };
  tran_low_t *outptr = out;
  int i, j, test;
  __m128i inptr[8];
  __m128i min_input, max_input, temp1, temp2, sign_bits;
  const __m128i zero = _mm_set1_epi16(0);
  const __m128i sixteen = _mm_set1_epi16(16);
  const __m128i max = _mm_set1_epi16(6201);
  const __m128i min = _mm_set1_epi16(-6201);
  int optimised_cols = 0;

  // Load input into __m128i & pack to 16 bits
  for (i = 0; i < 8; i++) {
    temp1 = _mm_loadu_si128((const __m128i *)(input + 8 * i));
    temp2 = _mm_loadu_si128((const __m128i *)(input + 8 * i + 4));
    inptr[i] = _mm_packs_epi32(temp1, temp2);
  }

  // Find the min & max for the row transform
  // only first 4 row has non-zero coefs
  max_input = _mm_max_epi16(inptr[0], inptr[1]);
  min_input = _mm_min_epi16(inptr[0], inptr[1]);
  for (i = 2; i < 4; i++) {
    max_input = _mm_max_epi16(max_input, inptr[i]);
    min_input = _mm_min_epi16(min_input, inptr[i]);
  }
  max_input = _mm_cmpgt_epi16(max_input, max);
  min_input = _mm_cmplt_epi16(min_input, min);
  temp1 = _mm_or_si128(max_input, min_input);
  test = _mm_movemask_epi8(temp1);

  if (!test) {
    // Do the row transform
    idct8_sse2(inptr);

    // Find the min & max for the column transform
    // N.B. Only first 4 cols contain non-zero coeffs
    max_input = _mm_max_epi16(inptr[0], inptr[1]);
    min_input = _mm_min_epi16(inptr[0], inptr[1]);
    for (i = 2; i < 8; i++) {
      max_input = _mm_max_epi16(max_input, inptr[i]);
      min_input = _mm_min_epi16(min_input, inptr[i]);
    }
    max_input = _mm_cmpgt_epi16(max_input, max);
    min_input = _mm_cmplt_epi16(min_input, min);
    temp1 = _mm_or_si128(max_input, min_input);
    test = _mm_movemask_epi8(temp1);

    if (test) {
      // Use fact only first 4 rows contain non-zero coeffs
      array_transpose_4X8(inptr, inptr);
      for (i = 0; i < 4; i++) {
        sign_bits = _mm_cmplt_epi16(inptr[i], zero);
        temp1 = _mm_unpackhi_epi16(inptr[i], sign_bits);
        temp2 = _mm_unpacklo_epi16(inptr[i], sign_bits);
        _mm_storeu_si128((__m128i *)(outptr + 4 * (2 * i + 1)), temp1);
        _mm_storeu_si128((__m128i *)(outptr + 4 * (2 * i)), temp2);
      }
    } else {
      // Set to use the optimised transform for the column
      optimised_cols = 1;
    }
  } else {
    // Run the un-optimised row transform
    for (i = 0; i < 4; ++i) {
      vpx_highbd_idct8_c(input, outptr, bd);
      input += 8;
      outptr += 8;
    }
  }

  if (optimised_cols) {
    idct8_sse2(inptr);

    // Final round & shift and Reconstruction and Store
    {
      __m128i d[8];
      for (i = 0; i < 8; i++) {
        inptr[i] = _mm_add_epi16(inptr[i], sixteen);
        d[i] = _mm_loadu_si128((const __m128i *)(dest + stride * i));
        inptr[i] = _mm_srai_epi16(inptr[i], 5);
        d[i] = clamp_high_sse2(_mm_adds_epi16(d[i], inptr[i]), bd);
        // Store
        _mm_storeu_si128((__m128i *)(dest + stride * i), d[i]);
      }
    }
  } else {
    // Run the un-optimised column transform
    tran_low_t temp_in[8], temp_out[8];
    for (i = 0; i < 8; ++i) {
      for (j = 0; j < 8; ++j) temp_in[j] = out[j * 8 + i];
      vpx_highbd_idct8_c(temp_in, temp_out, bd);
      for (j = 0; j < 8; ++j) {
        dest[j * stride + i] = highbd_clip_pixel_add(
            dest[j * stride + i], ROUND_POWER_OF_TWO(temp_out[j], 5), bd);
      }
    }
  }
}

void vpx_highbd_idct8x8_1_add_sse2(const tran_low_t *input, uint16_t *dest,
                                   int stride, int bd) {
  highbd_idct_1_add_kernel(input, dest, stride, bd, 8);
}
