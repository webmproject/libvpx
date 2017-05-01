/*
 *  Copyright (c) 2014 The WebM project authors. All Rights Reserved.
 *
 *  Usee of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include <assert.h>
#include <immintrin.h>

#include "./vp9_rtcd.h"
#include "vpx/vpx_integer.h"
#include "vpx_dsp/vpx_dsp_common.h"
#include "vpx_dsp/x86/bitdepth_conversion_avx2.h"

int64_t vp9_block_error_avx2(const tran_low_t *coeff, const tran_low_t *dqcoeff,
                             intptr_t block_size, int64_t *ssz) {
  __m256i sse_reg, ssz_reg;
  __m256i exp_dqcoeff_lo, exp_dqcoeff_hi, exp_coeff_lo, exp_coeff_hi;
  __m256i sse_reg_64hi, ssz_reg_64hi;
  __m128i sse_reg128, ssz_reg128;
  int64_t sse;
  const __m256i zero_reg = _mm256_setzero_si256();

  // If the block size is 16 then the results will fit in 32 bits.
  if (block_size == 16) {
    __m256i coeff_reg, dqcoeff_reg, coeff_reg_hi, dqcoeff_reg_hi;
    // Load 16 elements for coeff and dqcoeff.
    coeff_reg = load_tran_low(coeff);
    dqcoeff_reg = load_tran_low(dqcoeff);
    // dqcoeff - coeff
    dqcoeff_reg = _mm256_sub_epi16(dqcoeff_reg, coeff_reg);
    // madd (dqcoeff - coeff)
    dqcoeff_reg = _mm256_madd_epi16(dqcoeff_reg, dqcoeff_reg);
    // madd coeff
    coeff_reg = _mm256_madd_epi16(coeff_reg, coeff_reg);
    // Save the higher 64 bit of each 128 bit lane.
    dqcoeff_reg_hi = _mm256_srli_si256(dqcoeff_reg, 8);
    coeff_reg_hi = _mm256_srli_si256(coeff_reg, 8);
    // Add the higher 64 bit to the low 64 bit.
    dqcoeff_reg = _mm256_add_epi32(dqcoeff_reg, dqcoeff_reg_hi);
    coeff_reg = _mm256_add_epi32(coeff_reg, coeff_reg_hi);
    // Expand each double word in the lower 64 bits to quad word.
    sse_reg = _mm256_unpacklo_epi32(dqcoeff_reg, zero_reg);
    ssz_reg = _mm256_unpacklo_epi32(coeff_reg, zero_reg);
  } else {
    int i;
    assert(block_size % 32 == 0);
    sse_reg = zero_reg;
    ssz_reg = zero_reg;

    for (i = 0; i < block_size; i += 32) {
      __m256i coeff_reg_0, coeff_reg_1, dqcoeff_reg_0, dqcoeff_reg_1;
      // Load 32 elements for coeff and dqcoeff.
      coeff_reg_0 = load_tran_low(coeff + i);
      dqcoeff_reg_0 = load_tran_low(dqcoeff + i);
      coeff_reg_1 = load_tran_low(coeff + i + 16);
      dqcoeff_reg_1 = load_tran_low(dqcoeff + i + 16);
      // dqcoeff - coeff
      dqcoeff_reg_0 = _mm256_sub_epi16(dqcoeff_reg_0, coeff_reg_0);
      dqcoeff_reg_1 = _mm256_sub_epi16(dqcoeff_reg_1, coeff_reg_1);
      // madd (dqcoeff - coeff)
      dqcoeff_reg_0 = _mm256_madd_epi16(dqcoeff_reg_0, dqcoeff_reg_0);
      dqcoeff_reg_1 = _mm256_madd_epi16(dqcoeff_reg_1, dqcoeff_reg_1);
      // madd coeff
      coeff_reg_0 = _mm256_madd_epi16(coeff_reg_0, coeff_reg_0);
      coeff_reg_1 = _mm256_madd_epi16(coeff_reg_1, coeff_reg_1);
      // Add the first madd (dqcoeff - coeff) with the second.
      dqcoeff_reg_0 = _mm256_add_epi32(dqcoeff_reg_0, dqcoeff_reg_1);
      // Add the first madd (coeff) with the second.
      coeff_reg_0 = _mm256_add_epi32(coeff_reg_0, coeff_reg_1);
      // Expand each double word of madd (dqcoeff - coeff) to quad word.
      exp_dqcoeff_lo = _mm256_unpacklo_epi32(dqcoeff_reg_0, zero_reg);
      exp_dqcoeff_hi = _mm256_unpackhi_epi32(dqcoeff_reg_0, zero_reg);
      // expand each double word of madd (coeff) to quad word
      exp_coeff_lo = _mm256_unpacklo_epi32(coeff_reg_0, zero_reg);
      exp_coeff_hi = _mm256_unpackhi_epi32(coeff_reg_0, zero_reg);
      // Add each quad word of madd (dqcoeff - coeff) and madd (coeff).
      sse_reg = _mm256_add_epi64(sse_reg, exp_dqcoeff_lo);
      ssz_reg = _mm256_add_epi64(ssz_reg, exp_coeff_lo);
      sse_reg = _mm256_add_epi64(sse_reg, exp_dqcoeff_hi);
      ssz_reg = _mm256_add_epi64(ssz_reg, exp_coeff_hi);
    }
  }
  // Save the higher 64 bit of each 128 bit lane.
  sse_reg_64hi = _mm256_srli_si256(sse_reg, 8);
  ssz_reg_64hi = _mm256_srli_si256(ssz_reg, 8);
  // Add the higher 64 bit to the low 64 bit.
  sse_reg = _mm256_add_epi64(sse_reg, sse_reg_64hi);
  ssz_reg = _mm256_add_epi64(ssz_reg, ssz_reg_64hi);

  // Add each 64 bit from each of the 128 bit lane of the 256 bit.
  sse_reg128 = _mm_add_epi64(_mm256_castsi256_si128(sse_reg),
                             _mm256_extractf128_si256(sse_reg, 1));

  ssz_reg128 = _mm_add_epi64(_mm256_castsi256_si128(ssz_reg),
                             _mm256_extractf128_si256(ssz_reg, 1));

  // Store the results.
  _mm_storel_epi64((__m128i *)(&sse), sse_reg128);

  _mm_storel_epi64((__m128i *)(ssz), ssz_reg128);
  return sse;
}
