/*
 *  Copyright (c) 2016 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#ifndef VPX_DSP_X86_SYNONYS_H_
#define VPX_DSP_X86_SYNONYS_H_

#include <immintrin.h>

#include "./vpx_config.h"
#include "vpx/vpx_integer.h"

/**
 * Various reusable shorthands for x86 SIMD intrinsics.
 *
 * Intrinsics prefixed with xx_ operate on or return 128bit XMM registers.
 * Intrinsics prefixed with yy_ operate on or return 256bit YMM registers.
 */

// Loads and stores to do away with the tedium of casting the address
// to the right type.
static INLINE __m128i xx_loadl_32(const void *a) {
  return _mm_cvtsi32_si128(*(const uint32_t*)a);
}

static INLINE __m128i xx_loadl_64(const void *a) {
  return _mm_loadl_epi64((const __m128i*)a);
}

static INLINE __m128i xx_load_128(const void *a) {
  return _mm_load_si128((const __m128i*)a);
}

static INLINE __m128i xx_loadu_128(const void *a) {
  return _mm_loadu_si128((const __m128i*)a);
}

static INLINE void xx_storel_32(void *const a, const __m128i v) {
  *(uint32_t*)a = _mm_cvtsi128_si32(v);
}

static INLINE void xx_storel_64(void *const a, const __m128i v) {
  _mm_storel_epi64((__m128i*)a, v);
}

static INLINE void xx_store_128(void *const a, const __m128i v) {
  _mm_store_si128((__m128i*)a, v);
}

static INLINE void xx_storeu_128(void *const a, const __m128i v) {
  _mm_storeu_si128((__m128i*)a, v);
}

static INLINE __m128i xx_round_epu16(__m128i v_val_w) {
  return _mm_avg_epu16(v_val_w, _mm_setzero_si128());
}

static INLINE __m128i xx_roundn_epu16(__m128i v_val_w, int bits) {
  const __m128i v_s_w =_mm_srli_epi16(v_val_w, bits-1);
  return _mm_avg_epu16(v_s_w, _mm_setzero_si128());
}

static INLINE __m128i xx_roundn_epu32(__m128i v_val_d, int bits) {
  const __m128i v_bias_d = _mm_set1_epi32(1 << (bits - 1));
  const __m128i v_tmp_d = _mm_add_epi32(v_val_d, v_bias_d);
  return _mm_srli_epi32(v_tmp_d, bits);
}

#ifdef __SSSE3__
static INLINE int32_t xx_hsum_epi32_si32(__m128i v_d) {
  v_d = _mm_hadd_epi32(v_d, v_d);
  v_d = _mm_hadd_epi32(v_d, v_d);
  return _mm_cvtsi128_si32(v_d);
}
#endif  // __SSSE3__

#endif  // VPX_DSP_X86_SYNONYS_H_
