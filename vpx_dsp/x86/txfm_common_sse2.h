/*
 *  Copyright (c) 2015 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#ifndef VPX_DSP_X86_TXFM_COMMON_SSE2_H_
#define VPX_DSP_X86_TXFM_COMMON_SSE2_H_

#include <emmintrin.h>
#include "vpx/vpx_integer.h"

#define pair_set_epi16(a, b) \
  _mm_set_epi16((int16_t)(b), (int16_t)(a), (int16_t)(b), (int16_t)(a), \
                (int16_t)(b), (int16_t)(a), (int16_t)(b), (int16_t)(a))

#define dual_set_epi16(a, b) \
  _mm_set_epi16((int16_t)(b), (int16_t)(b), (int16_t)(b), (int16_t)(b), \
                (int16_t)(a), (int16_t)(a), (int16_t)(a), (int16_t)(a))

#define octa_set_epi16(a, b, c, d, e, f, g, h) \
  _mm_setr_epi16((int16_t)(a), (int16_t)(b), (int16_t)(c), (int16_t)(d), \
                 (int16_t)(e), (int16_t)(f), (int16_t)(g), (int16_t)(h))

// Reverse the 8 16 bit words in __m128i
static INLINE __m128i mm_reverse_epi16(const __m128i x) {
  const __m128i a = _mm_shufflelo_epi16(x, 0x1b);
  const __m128i b = _mm_shufflehi_epi16(a, 0x1b);
  return _mm_shuffle_epi32(b, 0x4e);
}

#endif  // VPX_DSP_X86_TXFM_COMMON_SSE2_H_
