/*
 *  Copyright (c) 2017 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#ifndef VPX_DSP_X86_CONVOLVE_AVX2_H_
#define VPX_DSP_X86_CONVOLVE_AVX2_H_

#include <immintrin.h>  // AVX2

#include "./vpx_config.h"

#if defined(__clang__)
#if (__clang_major__ > 0 && __clang_major__ < 3) ||            \
    (__clang_major__ == 3 && __clang_minor__ <= 3) ||          \
    (defined(__APPLE__) && defined(__apple_build_version__) && \
     ((__clang_major__ == 4 && __clang_minor__ <= 2) ||        \
      (__clang_major__ == 5 && __clang_minor__ == 0)))
#define MM256_BROADCASTSI128_SI256(x) \
  _mm_broadcastsi128_si256((__m128i const *)&(x))
#else  // clang > 3.3, and not 5.0 on macosx.
#define MM256_BROADCASTSI128_SI256(x) _mm256_broadcastsi128_si256(x)
#endif  // clang <= 3.3
#elif defined(__GNUC__)
#if __GNUC__ < 4 || (__GNUC__ == 4 && __GNUC_MINOR__ <= 6)
#define MM256_BROADCASTSI128_SI256(x) \
  _mm_broadcastsi128_si256((__m128i const *)&(x))
#elif __GNUC__ == 4 && __GNUC_MINOR__ == 7
#define MM256_BROADCASTSI128_SI256(x) _mm_broadcastsi128_si256(x)
#else  // gcc > 4.7
#define MM256_BROADCASTSI128_SI256(x) _mm256_broadcastsi128_si256(x)
#endif  // gcc <= 4.6
#else   // !(gcc || clang)
#define MM256_BROADCASTSI128_SI256(x) _mm256_broadcastsi128_si256(x)
#endif  // __clang__

static INLINE void shuffle_filter_avx2(const int16_t *const filter,
                                       __m256i *const f) {
  const __m256i f_values =
      MM256_BROADCASTSI128_SI256(_mm_load_si128((const __m128i *)filter));
  // pack and duplicate the filter values
  f[0] = _mm256_shuffle_epi8(f_values, _mm256_set1_epi16(0x0200u));
  f[1] = _mm256_shuffle_epi8(f_values, _mm256_set1_epi16(0x0604u));
  f[2] = _mm256_shuffle_epi8(f_values, _mm256_set1_epi16(0x0a08u));
  f[3] = _mm256_shuffle_epi8(f_values, _mm256_set1_epi16(0x0e0cu));
}

static INLINE __m256i convolve8_16_avx2(const __m256i *const s,
                                        const __m256i *const f) {
  // multiply 2 adjacent elements with the filter and add the result
  const __m256i k_64 = _mm256_set1_epi16(1 << 6);
  const __m256i x0 = _mm256_maddubs_epi16(s[0], f[0]);
  const __m256i x1 = _mm256_maddubs_epi16(s[1], f[1]);
  const __m256i x2 = _mm256_maddubs_epi16(s[2], f[2]);
  const __m256i x3 = _mm256_maddubs_epi16(s[3], f[3]);
  // add and saturate the results together
  const __m256i min_x2x1 = _mm256_min_epi16(x2, x1);
  const __m256i max_x2x1 = _mm256_max_epi16(x2, x1);
  __m256i temp = _mm256_adds_epi16(x0, x3);
  temp = _mm256_adds_epi16(temp, min_x2x1);
  temp = _mm256_adds_epi16(temp, max_x2x1);
  // round and shift by 7 bit each 16 bit
  temp = _mm256_adds_epi16(temp, k_64);
  temp = _mm256_srai_epi16(temp, 7);
  return temp;
}

static INLINE __m128i convolve8_8_avx2(const __m256i *const s,
                                       const __m256i *const f) {
  // multiply 2 adjacent elements with the filter and add the result
  const __m128i k_64 = _mm_set1_epi16(1 << 6);
  const __m128i x0 = _mm_maddubs_epi16(_mm256_castsi256_si128(s[0]),
                                       _mm256_castsi256_si128(f[0]));
  const __m128i x1 = _mm_maddubs_epi16(_mm256_castsi256_si128(s[1]),
                                       _mm256_castsi256_si128(f[1]));
  const __m128i x2 = _mm_maddubs_epi16(_mm256_castsi256_si128(s[2]),
                                       _mm256_castsi256_si128(f[2]));
  const __m128i x3 = _mm_maddubs_epi16(_mm256_castsi256_si128(s[3]),
                                       _mm256_castsi256_si128(f[3]));
  // add and saturate the results together
  const __m128i min_x2x1 = _mm_min_epi16(x2, x1);
  const __m128i max_x2x1 = _mm_max_epi16(x2, x1);
  __m128i temp = _mm_adds_epi16(x0, x3);
  temp = _mm_adds_epi16(temp, min_x2x1);
  temp = _mm_adds_epi16(temp, max_x2x1);
  // round and shift by 7 bit each 16 bit
  temp = _mm_adds_epi16(temp, k_64);
  temp = _mm_srai_epi16(temp, 7);
  return temp;
}

#undef MM256_BROADCASTSI128_SI256

#endif  // VPX_DSP_X86_CONVOLVE_AVX2_H_
