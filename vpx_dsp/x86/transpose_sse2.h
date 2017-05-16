/*
 *  Copyright (c) 2015 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#ifndef VPX_DSP_X86_TRANSPOSE_SSE2_H_
#define VPX_DSP_X86_TRANSPOSE_SSE2_H_

#include "./vpx_dsp_rtcd.h"
#include "vpx_dsp/x86/inv_txfm_sse2.h"
#include "vpx_dsp/x86/txfm_common_sse2.h"

static INLINE void transpose_16bit_4x4(__m128i *res) {
  const __m128i tr0_0 = _mm_unpacklo_epi16(res[0], res[1]);
  const __m128i tr0_1 = _mm_unpackhi_epi16(res[0], res[1]);

  res[0] = _mm_unpacklo_epi16(tr0_0, tr0_1);
  res[1] = _mm_unpackhi_epi16(tr0_0, tr0_1);
}

static INLINE void transpose_32bit_4x4(__m128i *const a0, __m128i *const a1,
                                       __m128i *const a2, __m128i *const a3) {
  // Unpack 32 bit elements. Goes from:
  // a0: 00 01 02 03
  // a1: 10 11 12 13
  // a2: 20 21 22 23
  // a3: 30 31 32 33
  // to:
  // b0: 00 10 01 11
  // b1: 20 30 21 31
  // b2: 02 12 03 13
  // b3: 22 32 23 33

  const __m128i b0 = _mm_unpacklo_epi32(*a0, *a1);
  const __m128i b1 = _mm_unpacklo_epi32(*a2, *a3);
  const __m128i b2 = _mm_unpackhi_epi32(*a0, *a1);
  const __m128i b3 = _mm_unpackhi_epi32(*a2, *a3);

  // Unpack 64 bit elements resulting in:
  // a0: 00 10 20 30
  // a1: 01 11 21 31
  // a2: 02 12 22 32
  // a3: 03 13 23 33
  *a0 = _mm_unpacklo_epi64(b0, b1);
  *a1 = _mm_unpackhi_epi64(b0, b1);
  *a2 = _mm_unpacklo_epi64(b2, b3);
  *a3 = _mm_unpackhi_epi64(b2, b3);
}

#endif  // VPX_DSP_X86_TRANSPOSE_SSE2_H_
