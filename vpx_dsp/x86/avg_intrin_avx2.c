/*
 *  Copyright (c) 2017 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include <immintrin.h>

#include "./vpx_dsp_rtcd.h"
#include "vpx/vpx_integer.h"
#include "vpx_dsp/x86/bitdepth_conversion_avx2.h"
#include "vpx_ports/mem.h"

static void hadamard_col8x2_avx2(__m256i *in, int iter) {
  __m256i a0 = in[0];
  __m256i a1 = in[1];
  __m256i a2 = in[2];
  __m256i a3 = in[3];
  __m256i a4 = in[4];
  __m256i a5 = in[5];
  __m256i a6 = in[6];
  __m256i a7 = in[7];

  __m256i b0 = _mm256_add_epi16(a0, a1);
  __m256i b1 = _mm256_sub_epi16(a0, a1);
  __m256i b2 = _mm256_add_epi16(a2, a3);
  __m256i b3 = _mm256_sub_epi16(a2, a3);
  __m256i b4 = _mm256_add_epi16(a4, a5);
  __m256i b5 = _mm256_sub_epi16(a4, a5);
  __m256i b6 = _mm256_add_epi16(a6, a7);
  __m256i b7 = _mm256_sub_epi16(a6, a7);

  a0 = _mm256_add_epi16(b0, b2);
  a1 = _mm256_add_epi16(b1, b3);
  a2 = _mm256_sub_epi16(b0, b2);
  a3 = _mm256_sub_epi16(b1, b3);
  a4 = _mm256_add_epi16(b4, b6);
  a5 = _mm256_add_epi16(b5, b7);
  a6 = _mm256_sub_epi16(b4, b6);
  a7 = _mm256_sub_epi16(b5, b7);

  if (iter == 0) {
    b0 = _mm256_add_epi16(a0, a4);
    b7 = _mm256_add_epi16(a1, a5);
    b3 = _mm256_add_epi16(a2, a6);
    b4 = _mm256_add_epi16(a3, a7);
    b2 = _mm256_sub_epi16(a0, a4);
    b6 = _mm256_sub_epi16(a1, a5);
    b1 = _mm256_sub_epi16(a2, a6);
    b5 = _mm256_sub_epi16(a3, a7);

    a0 = _mm256_unpacklo_epi16(b0, b1);
    a1 = _mm256_unpacklo_epi16(b2, b3);
    a2 = _mm256_unpackhi_epi16(b0, b1);
    a3 = _mm256_unpackhi_epi16(b2, b3);
    a4 = _mm256_unpacklo_epi16(b4, b5);
    a5 = _mm256_unpacklo_epi16(b6, b7);
    a6 = _mm256_unpackhi_epi16(b4, b5);
    a7 = _mm256_unpackhi_epi16(b6, b7);

    b0 = _mm256_unpacklo_epi32(a0, a1);
    b1 = _mm256_unpacklo_epi32(a4, a5);
    b2 = _mm256_unpackhi_epi32(a0, a1);
    b3 = _mm256_unpackhi_epi32(a4, a5);
    b4 = _mm256_unpacklo_epi32(a2, a3);
    b5 = _mm256_unpacklo_epi32(a6, a7);
    b6 = _mm256_unpackhi_epi32(a2, a3);
    b7 = _mm256_unpackhi_epi32(a6, a7);

    in[0] = _mm256_unpacklo_epi64(b0, b1);
    in[1] = _mm256_unpackhi_epi64(b0, b1);
    in[2] = _mm256_unpacklo_epi64(b2, b3);
    in[3] = _mm256_unpackhi_epi64(b2, b3);
    in[4] = _mm256_unpacklo_epi64(b4, b5);
    in[5] = _mm256_unpackhi_epi64(b4, b5);
    in[6] = _mm256_unpacklo_epi64(b6, b7);
    in[7] = _mm256_unpackhi_epi64(b6, b7);
  } else {
    in[0] = _mm256_add_epi16(a0, a4);
    in[7] = _mm256_add_epi16(a1, a5);
    in[3] = _mm256_add_epi16(a2, a6);
    in[4] = _mm256_add_epi16(a3, a7);
    in[2] = _mm256_sub_epi16(a0, a4);
    in[6] = _mm256_sub_epi16(a1, a5);
    in[1] = _mm256_sub_epi16(a2, a6);
    in[5] = _mm256_sub_epi16(a3, a7);
  }
}

static void hadamard_8x8x2_avx2(int16_t const *src_diff, int src_stride,
                                tran_low_t *coeff) {
  __m256i src[8];
  src[0] = _mm256_loadu_si256((const __m256i *)src_diff);
  src[1] = _mm256_loadu_si256((const __m256i *)(src_diff += src_stride));
  src[2] = _mm256_loadu_si256((const __m256i *)(src_diff += src_stride));
  src[3] = _mm256_loadu_si256((const __m256i *)(src_diff += src_stride));
  src[4] = _mm256_loadu_si256((const __m256i *)(src_diff += src_stride));
  src[5] = _mm256_loadu_si256((const __m256i *)(src_diff += src_stride));
  src[6] = _mm256_loadu_si256((const __m256i *)(src_diff += src_stride));
  src[7] = _mm256_loadu_si256((const __m256i *)(src_diff += src_stride));

  hadamard_col8x2_avx2(src, 0);
  hadamard_col8x2_avx2(src, 1);

  store_tran_low(_mm256_castsi256_si128(src[0]), coeff);
  coeff += 8;
  store_tran_low(_mm256_castsi256_si128(src[1]), coeff);
  coeff += 8;
  store_tran_low(_mm256_castsi256_si128(src[2]), coeff);
  coeff += 8;
  store_tran_low(_mm256_castsi256_si128(src[3]), coeff);
  coeff += 8;
  store_tran_low(_mm256_castsi256_si128(src[4]), coeff);
  coeff += 8;
  store_tran_low(_mm256_castsi256_si128(src[5]), coeff);
  coeff += 8;
  store_tran_low(_mm256_castsi256_si128(src[6]), coeff);
  coeff += 8;
  store_tran_low(_mm256_castsi256_si128(src[7]), coeff);
  coeff += 8;

  src[0] = _mm256_castsi128_si256(_mm256_extractf128_si256(src[0], 1));
  src[1] = _mm256_castsi128_si256(_mm256_extractf128_si256(src[1], 1));
  src[2] = _mm256_castsi128_si256(_mm256_extractf128_si256(src[2], 1));
  src[3] = _mm256_castsi128_si256(_mm256_extractf128_si256(src[3], 1));
  src[4] = _mm256_castsi128_si256(_mm256_extractf128_si256(src[4], 1));
  src[5] = _mm256_castsi128_si256(_mm256_extractf128_si256(src[5], 1));
  src[6] = _mm256_castsi128_si256(_mm256_extractf128_si256(src[6], 1));
  src[7] = _mm256_castsi128_si256(_mm256_extractf128_si256(src[7], 1));

  store_tran_low(_mm256_castsi256_si128(src[0]), coeff);
  coeff += 8;
  store_tran_low(_mm256_castsi256_si128(src[1]), coeff);
  coeff += 8;
  store_tran_low(_mm256_castsi256_si128(src[2]), coeff);
  coeff += 8;
  store_tran_low(_mm256_castsi256_si128(src[3]), coeff);
  coeff += 8;
  store_tran_low(_mm256_castsi256_si128(src[4]), coeff);
  coeff += 8;
  store_tran_low(_mm256_castsi256_si128(src[5]), coeff);
  coeff += 8;
  store_tran_low(_mm256_castsi256_si128(src[6]), coeff);
  coeff += 8;
  store_tran_low(_mm256_castsi256_si128(src[7]), coeff);
}

void vpx_hadamard_16x16_avx2(int16_t const *src_diff, int src_stride,
                             tran_low_t *coeff) {
  int idx;
  for (idx = 0; idx < 2; ++idx) {
    int16_t const *src_ptr = src_diff + idx * 8 * src_stride;
    hadamard_8x8x2_avx2(src_ptr, src_stride, coeff + (idx * 64 * 2));
  }

  for (idx = 0; idx < 64; idx += 16) {
    const __m256i coeff0 = load_tran_low(coeff);
    const __m256i coeff1 = load_tran_low(coeff + 64);
    const __m256i coeff2 = load_tran_low(coeff + 128);
    const __m256i coeff3 = load_tran_low(coeff + 192);
    __m256i b0 = _mm256_add_epi16(coeff0, coeff1);
    __m256i b1 = _mm256_sub_epi16(coeff0, coeff1);
    __m256i b2 = _mm256_add_epi16(coeff2, coeff3);
    __m256i b3 = _mm256_sub_epi16(coeff2, coeff3);

    b0 = _mm256_srai_epi16(b0, 1);
    b1 = _mm256_srai_epi16(b1, 1);
    b2 = _mm256_srai_epi16(b2, 1);
    b3 = _mm256_srai_epi16(b3, 1);

    store_tran_low_256(_mm256_add_epi16(b0, b2), coeff);
    store_tran_low_256(_mm256_add_epi16(b1, b3), coeff + 64);
    store_tran_low_256(_mm256_sub_epi16(b0, b2), coeff + 128);
    store_tran_low_256(_mm256_sub_epi16(b1, b3), coeff + 192);

    coeff += 16;
  }
}
