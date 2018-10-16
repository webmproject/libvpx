/*
 *  Copyright (c) 2018 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include <emmintrin.h>

#include "./vpx_dsp_rtcd.h"
#include "vpx/vpx_integer.h"
#include "vpx_dsp/x86/convolve.h"
#include "vpx_dsp/x86/convolve_sse2.h"
#include "vpx_ports/mem.h"

void vpx_filter_block1d16_h4_sse2(const uint8_t *src_ptr, ptrdiff_t src_stride,
                                  uint8_t *dst_ptr, ptrdiff_t dst_stride,
                                  uint32_t height, const int16_t *kernel) {
  __m128i kernel_reg;                         // Kernel
  __m128i kernel_reg_23, kernel_reg_45;       // Segments of the kernel used
  const __m128i reg_32 = _mm_set1_epi16(32);  // Used for rounding
  int h;

  __m128i src_reg, src_reg_shift_1, src_reg_shift_2, src_reg_shift_3;
  __m128i dst_first, dst_second;
  __m128i even, odd;

  // Start one pixel before as we need tap/2 - 1 = 1 sample from the past
  src_ptr -= 1;

  // Load Kernel
  kernel_reg = _mm_loadu_si128((const __m128i *)kernel);
  kernel_reg = _mm_srai_epi16(kernel_reg, 1);
  kernel_reg_23 = extract_quarter_2_epi16_sse2(&kernel_reg);
  kernel_reg_45 = extract_quarter_3_epi16_sse2(&kernel_reg);

  for (h = height; h > 0; --h) {
    // We will load multiple shifted versions of the row and shuffle them into
    // 16-bit words of the form
    // ... s[2] s[1] s[0] s[-1]
    // ... s[4] s[3] s[2] s[1]
    // Then we call multiply and add to get partial results
    // s[2]k[3]+s[1]k[2] s[0]k[3]s[-1]k[2]
    // s[4]k[5]+s[3]k[4] s[2]k[5]s[1]k[4]
    // The two results are then added together for the first half of even
    // output.
    // Repeat multiple times to get the whole outoput
    src_reg = _mm_loadu_si128((const __m128i *)src_ptr);
    src_reg_shift_1 = _mm_srli_si128(src_reg, 1);
    src_reg_shift_2 = _mm_srli_si128(src_reg, 2);
    src_reg_shift_3 = _mm_srli_si128(src_reg, 3);

    // Output 6 4 2 0
    even = pad_multiply_add_add_epi8_sse2(&src_reg, &src_reg_shift_2,
                                          &kernel_reg_23, &kernel_reg_45);

    // Output 7 5 3 1
    odd = pad_multiply_add_add_epi8_sse2(&src_reg_shift_1, &src_reg_shift_3,
                                         &kernel_reg_23, &kernel_reg_45);

    // Combine to get the first half of the dst
    dst_first = combine_epi32_sse2(&even, &odd);

    // Do again to get the second half of dst
    src_reg = _mm_loadu_si128((const __m128i *)(src_ptr + 8));
    src_reg_shift_1 = _mm_srli_si128(src_reg, 1);
    src_reg_shift_2 = _mm_srli_si128(src_reg, 2);
    src_reg_shift_3 = _mm_srli_si128(src_reg, 3);

    // Output 14 12 10 8
    even = pad_multiply_add_add_epi8_sse2(&src_reg, &src_reg_shift_2,
                                          &kernel_reg_23, &kernel_reg_45);

    // Output 15 13 11 9
    odd = pad_multiply_add_add_epi8_sse2(&src_reg_shift_1, &src_reg_shift_3,
                                         &kernel_reg_23, &kernel_reg_45);

    // Combine to get the second half of the dst
    dst_second = combine_epi32_sse2(&even, &odd);

    // Round each result
    dst_first = round_epi16_sse2(&dst_first, &reg_32, 6);
    dst_second = round_epi16_sse2(&dst_second, &reg_32, 6);

    // Finally combine to get the final dst
    dst_first = _mm_packus_epi16(dst_first, dst_second);
    _mm_store_si128((__m128i *)dst_ptr, dst_first);

    src_ptr += src_stride;
    dst_ptr += dst_stride;
  }
}

/* The macro used to generate functions shifts the src_ptr up by 3 rows already
 * */

void vpx_filter_block1d16_v4_sse2(const uint8_t *src_ptr, ptrdiff_t src_stride,
                                  uint8_t *dst_ptr, ptrdiff_t dst_stride,
                                  uint32_t height, const int16_t *kernel) {
  // Register for source s[-1:3, :]
  __m128i src_reg_m1, src_reg_0, src_reg_1, src_reg_2, src_reg_3;
  // Interleaved rows of the source. lo is first half, hi second
  __m128i src_reg_m10_lo, src_reg_m10_hi, src_reg_01_lo, src_reg_01_hi;
  __m128i src_reg_12_lo, src_reg_12_hi, src_reg_23_lo, src_reg_23_hi;
  // Half of half of the interleaved rows
  __m128i src_reg_m10_lo_1, src_reg_m10_lo_2, src_reg_m10_hi_1,
      src_reg_m10_hi_2;
  __m128i src_reg_01_lo_1, src_reg_01_lo_2, src_reg_01_hi_1, src_reg_01_hi_2;
  __m128i src_reg_12_lo_1, src_reg_12_lo_2, src_reg_12_hi_1, src_reg_12_hi_2;
  __m128i src_reg_23_lo_1, src_reg_23_lo_2, src_reg_23_hi_1, src_reg_23_hi_2;

  __m128i kernel_reg;                    // Kernel
  __m128i kernel_reg_23, kernel_reg_45;  // Segments of the kernel used

  // Result after multiply and add
  __m128i res_reg_m10_lo, res_reg_01_lo, res_reg_12_lo, res_reg_23_lo;
  __m128i res_reg_m10_hi, res_reg_01_hi, res_reg_12_hi, res_reg_23_hi;
  __m128i res_reg_m1012, res_reg_0123;
  __m128i res_reg_m1012_lo, res_reg_0123_lo, res_reg_m1012_hi, res_reg_0123_hi;

  const __m128i reg_32 = _mm_set1_epi16(32);  // Used for rounding

  // We will compute the result two rows at a time
  const ptrdiff_t src_stride_unrolled = src_stride << 1;
  const ptrdiff_t dst_stride_unrolled = dst_stride << 1;
  int h;

  // We only need to go num_taps/2 - 1 row above the souce, so we move
  // 3 - (num_taps/2 - 1) = 4 - num_taps/2 = 2 back down
  src_ptr += src_stride_unrolled;

  // Load Kernel
  kernel_reg = _mm_loadu_si128((const __m128i *)kernel);
  kernel_reg = _mm_srai_epi16(kernel_reg, 1);
  kernel_reg_23 = extract_quarter_2_epi16_sse2(&kernel_reg);
  kernel_reg_45 = extract_quarter_3_epi16_sse2(&kernel_reg);

  // We will load two rows of pixels as 8-bit words, rearrange them as 16-bit
  // words,
  // shuffle the data into the form
  // ... s[0,1] s[-1,1] s[0,0] s[-1,0]
  // ... s[0,7] s[-1,7] s[0,6] s[-1,6]
  // ... s[0,9] s[-1,9] s[0,8] s[-1,8]
  // ... s[0,13] s[-1,13] s[0,12] s[-1,12]
  // so that we can call multiply and add with the kernel to get 32-bit words of
  // the form
  // ... s[0,1]k[3]+s[-1,1]k[2] s[0,0]k[3]+s[-1,0]k[2]
  // Finally, we can add multiple rows together to get the desired output.

  // First shuffle the data
  src_reg_m1 = _mm_loadu_si128((const __m128i *)src_ptr);
  src_reg_0 = _mm_loadu_si128((const __m128i *)(src_ptr + src_stride));
  src_reg_m10_lo = _mm_unpacklo_epi8(src_reg_m1, src_reg_0);
  src_reg_m10_hi = _mm_unpackhi_epi8(src_reg_m1, src_reg_0);
  src_reg_m10_lo_1 = _mm_unpacklo_epi8(src_reg_m10_lo, _mm_setzero_si128());
  src_reg_m10_lo_2 = _mm_unpackhi_epi8(src_reg_m10_lo, _mm_setzero_si128());
  src_reg_m10_hi_1 = _mm_unpacklo_epi8(src_reg_m10_hi, _mm_setzero_si128());
  src_reg_m10_hi_2 = _mm_unpackhi_epi8(src_reg_m10_hi, _mm_setzero_si128());

  // More shuffling
  src_reg_1 = _mm_loadu_si128((const __m128i *)(src_ptr + src_stride * 2));
  src_reg_01_lo = _mm_unpacklo_epi8(src_reg_0, src_reg_1);
  src_reg_01_hi = _mm_unpackhi_epi8(src_reg_0, src_reg_1);
  src_reg_01_lo_1 = _mm_unpacklo_epi8(src_reg_01_lo, _mm_setzero_si128());
  src_reg_01_lo_2 = _mm_unpackhi_epi8(src_reg_01_lo, _mm_setzero_si128());
  src_reg_01_hi_1 = _mm_unpacklo_epi8(src_reg_01_hi, _mm_setzero_si128());
  src_reg_01_hi_2 = _mm_unpackhi_epi8(src_reg_01_hi, _mm_setzero_si128());

  for (h = height; h > 1; h -= 2) {
    src_reg_2 = _mm_loadu_si128((const __m128i *)(src_ptr + src_stride * 3));

    src_reg_12_lo = _mm_unpacklo_epi8(src_reg_1, src_reg_2);
    src_reg_12_hi = _mm_unpackhi_epi8(src_reg_1, src_reg_2);

    src_reg_3 = _mm_loadu_si128((const __m128i *)(src_ptr + src_stride * 4));

    src_reg_23_lo = _mm_unpacklo_epi8(src_reg_2, src_reg_3);
    src_reg_23_hi = _mm_unpackhi_epi8(src_reg_2, src_reg_3);

    // Partial output from first half
    res_reg_m10_lo = multiply_add_packs_epi16_sse2(
        &src_reg_m10_lo_1, &src_reg_m10_lo_2, &kernel_reg_23);

    res_reg_01_lo = multiply_add_packs_epi16_sse2(
        &src_reg_01_lo_1, &src_reg_01_lo_2, &kernel_reg_23);

    src_reg_12_lo_1 = _mm_unpacklo_epi8(src_reg_12_lo, _mm_setzero_si128());
    src_reg_12_lo_2 = _mm_unpackhi_epi8(src_reg_12_lo, _mm_setzero_si128());
    res_reg_12_lo = multiply_add_packs_epi16_sse2(
        &src_reg_12_lo_1, &src_reg_12_lo_2, &kernel_reg_45);

    src_reg_23_lo_1 = _mm_unpacklo_epi8(src_reg_23_lo, _mm_setzero_si128());
    src_reg_23_lo_2 = _mm_unpackhi_epi8(src_reg_23_lo, _mm_setzero_si128());
    res_reg_23_lo = multiply_add_packs_epi16_sse2(
        &src_reg_23_lo_1, &src_reg_23_lo_2, &kernel_reg_45);

    // Add to get first half of the results
    res_reg_m1012_lo = _mm_adds_epi16(res_reg_m10_lo, res_reg_12_lo);
    res_reg_0123_lo = _mm_adds_epi16(res_reg_01_lo, res_reg_23_lo);

    // Now repeat everything again for the second half
    // Partial output for second half
    res_reg_m10_hi = multiply_add_packs_epi16_sse2(
        &src_reg_m10_hi_1, &src_reg_m10_hi_2, &kernel_reg_23);

    res_reg_01_hi = multiply_add_packs_epi16_sse2(
        &src_reg_01_hi_1, &src_reg_01_hi_2, &kernel_reg_23);

    src_reg_12_hi_1 = _mm_unpacklo_epi8(src_reg_12_hi, _mm_setzero_si128());
    src_reg_12_hi_2 = _mm_unpackhi_epi8(src_reg_12_hi, _mm_setzero_si128());
    res_reg_12_hi = multiply_add_packs_epi16_sse2(
        &src_reg_12_hi_1, &src_reg_12_hi_2, &kernel_reg_45);

    src_reg_23_hi_1 = _mm_unpacklo_epi8(src_reg_23_hi, _mm_setzero_si128());
    src_reg_23_hi_2 = _mm_unpackhi_epi8(src_reg_23_hi, _mm_setzero_si128());
    res_reg_23_hi = multiply_add_packs_epi16_sse2(
        &src_reg_23_hi_1, &src_reg_23_hi_2, &kernel_reg_45);

    // Second half of the results
    res_reg_m1012_hi = _mm_adds_epi16(res_reg_m10_hi, res_reg_12_hi);
    res_reg_0123_hi = _mm_adds_epi16(res_reg_01_hi, res_reg_23_hi);

    // Round the words
    res_reg_m1012_lo = round_epi16_sse2(&res_reg_m1012_lo, &reg_32, 6);
    res_reg_0123_lo = round_epi16_sse2(&res_reg_0123_lo, &reg_32, 6);
    res_reg_m1012_hi = round_epi16_sse2(&res_reg_m1012_hi, &reg_32, 6);
    res_reg_0123_hi = round_epi16_sse2(&res_reg_0123_hi, &reg_32, 6);

    // Combine to get the result
    res_reg_m1012 = _mm_packus_epi16(res_reg_m1012_lo, res_reg_m1012_hi);
    res_reg_0123 = _mm_packus_epi16(res_reg_0123_lo, res_reg_0123_hi);

    _mm_store_si128((__m128i *)dst_ptr, res_reg_m1012);
    _mm_store_si128((__m128i *)(dst_ptr + dst_stride), res_reg_0123);

    // Update the source by two rows
    src_ptr += src_stride_unrolled;
    dst_ptr += dst_stride_unrolled;

    src_reg_m10_lo_1 = src_reg_12_lo_1;
    src_reg_m10_lo_2 = src_reg_12_lo_2;
    src_reg_m10_hi_1 = src_reg_12_hi_1;
    src_reg_m10_hi_2 = src_reg_12_hi_2;
    src_reg_01_lo_1 = src_reg_23_lo_1;
    src_reg_01_lo_2 = src_reg_23_lo_2;
    src_reg_01_hi_1 = src_reg_23_hi_1;
    src_reg_01_hi_2 = src_reg_23_hi_2;
    src_reg_1 = src_reg_3;
  }
}
