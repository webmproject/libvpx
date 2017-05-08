/*
 *  Copyright (c) 2015 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#ifndef VPX_DSP_X86_HIGHBD_INV_TXFM_SSE2_H_
#define VPX_DSP_X86_HIGHBD_INV_TXFM_SSE2_H_

#include <emmintrin.h>  // SSE2
#include "./vpx_config.h"
#include "vpx/vpx_integer.h"
#include "vpx_dsp/inv_txfm.h"
#include "vpx_dsp/x86/txfm_common_sse2.h"

static INLINE __m128i clamp_high_sse2(__m128i value, int bd) {
  __m128i ubounded, retval;
  const __m128i zero = _mm_set1_epi16(0);
  const __m128i one = _mm_set1_epi16(1);
  const __m128i max = _mm_sub_epi16(_mm_slli_epi16(one, bd), one);
  ubounded = _mm_cmpgt_epi16(value, max);
  retval = _mm_andnot_si128(ubounded, value);
  ubounded = _mm_and_si128(ubounded, max);
  retval = _mm_or_si128(retval, ubounded);
  retval = _mm_and_si128(retval, _mm_cmpgt_epi16(retval, zero));
  return retval;
}

#endif  // VPX_DSP_X86_HIGHBD_INV_TXFM_SSE2_H_
