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

static INLINE void transpose_16bit_4x4(const __m128i *const in,
                                       __m128i *const out) {
  // Unpack 16 bit elements. Goes from:
  // in[0]: 00 01 02 03  XX XX XX XX
  // in[1]: 10 11 12 13  XX XX XX XX
  // in[2]: 20 21 22 23  XX XX XX XX
  // in[3]: 30 31 32 33  XX XX XX XX
  // to:
  // tr0_0: 00 10 01 11  02 12 03 13
  // tr0_1: 20 30 21 31  22 32 23 33
  const __m128i tr0_0 = _mm_unpacklo_epi16(in[0], in[1]);
  const __m128i tr0_1 = _mm_unpacklo_epi16(in[2], in[3]);

  // Unpack 32 bit elements resulting in:
  // out[0]: 00 10 20 30  01 11 21 31
  // out[1]: 02 12 22 32  03 13 23 33
  out[0] = _mm_unpacklo_epi32(tr0_0, tr0_1);
  out[1] = _mm_unpackhi_epi32(tr0_0, tr0_1);
}

static INLINE void transpose_16bit_8x8(const __m128i *const in,
                                       __m128i *const out) {
  // Unpack 16 bit elements. Goes from:
  // in[0]: 00 01 02 03  04 05 06 07
  // in[1]: 10 11 12 13  14 15 16 17
  // in[2]: 20 21 22 23  24 25 26 27
  // in[3]: 30 31 32 33  34 35 36 37
  // in[4]: 40 41 42 43  44 45 46 47
  // in[5]: 50 51 52 53  54 55 56 57
  // in[6]: 60 61 62 63  64 65 66 67
  // in[7]: 70 71 72 73  74 75 76 77
  // to:
  // tr0_0: 00 10 01 11  02 12 03 13
  // tr0_1: 20 30 21 31  22 32 23 33
  // tr0_2: 40 50 41 51  42 52 43 53
  // tr0_3: 60 70 61 71  62 72 63 73
  // tr0_4: 04 14 05 15  06 16 07 17
  // tr0_5: 24 34 25 35  26 36 27 37
  // tr0_6: 44 54 45 55  46 56 47 57
  // tr0_7: 64 74 65 75  66 76 67 77
  const __m128i tr0_0 = _mm_unpacklo_epi16(in[0], in[1]);
  const __m128i tr0_1 = _mm_unpacklo_epi16(in[2], in[3]);
  const __m128i tr0_2 = _mm_unpacklo_epi16(in[4], in[5]);
  const __m128i tr0_3 = _mm_unpacklo_epi16(in[6], in[7]);
  const __m128i tr0_4 = _mm_unpackhi_epi16(in[0], in[1]);
  const __m128i tr0_5 = _mm_unpackhi_epi16(in[2], in[3]);
  const __m128i tr0_6 = _mm_unpackhi_epi16(in[4], in[5]);
  const __m128i tr0_7 = _mm_unpackhi_epi16(in[6], in[7]);

  // Unpack 32 bit elements resulting in:
  // tr1_0: 00 10 20 30  01 11 21 31
  // tr1_1: 40 50 60 70  41 51 61 71
  // tr1_2: 04 14 24 34  05 15 25 35
  // tr1_3: 44 54 64 74  45 55 65 75
  // tr1_4: 02 12 22 32  03 13 23 33
  // tr1_5: 42 52 62 72  43 53 63 73
  // tr1_6: 06 16 26 36  07 17 27 37
  // tr1_7: 46 56 66 76  47 57 67 77
  const __m128i tr1_0 = _mm_unpacklo_epi32(tr0_0, tr0_1);
  const __m128i tr1_1 = _mm_unpacklo_epi32(tr0_2, tr0_3);
  const __m128i tr1_2 = _mm_unpacklo_epi32(tr0_4, tr0_5);
  const __m128i tr1_3 = _mm_unpacklo_epi32(tr0_6, tr0_7);
  const __m128i tr1_4 = _mm_unpackhi_epi32(tr0_0, tr0_1);
  const __m128i tr1_5 = _mm_unpackhi_epi32(tr0_2, tr0_3);
  const __m128i tr1_6 = _mm_unpackhi_epi32(tr0_4, tr0_5);
  const __m128i tr1_7 = _mm_unpackhi_epi32(tr0_6, tr0_7);

  // Unpack 64 bit elements resulting in:
  // out[0]: 00 10 20 30  40 50 60 70
  // out[1]: 01 11 21 31  41 51 61 71
  // out[2]: 02 12 22 32  42 52 62 72
  // out[3]: 03 13 23 33  43 53 63 73
  // out[4]: 04 14 24 34  44 54 64 74
  // out[5]: 05 15 25 35  45 55 65 75
  // out[6]: 06 16 26 36  46 56 66 76
  // out[7]: 07 17 27 37  47 57 67 77
  out[0] = _mm_unpacklo_epi64(tr1_0, tr1_1);
  out[1] = _mm_unpackhi_epi64(tr1_0, tr1_1);
  out[2] = _mm_unpacklo_epi64(tr1_4, tr1_5);
  out[3] = _mm_unpackhi_epi64(tr1_4, tr1_5);
  out[4] = _mm_unpacklo_epi64(tr1_2, tr1_3);
  out[5] = _mm_unpackhi_epi64(tr1_2, tr1_3);
  out[6] = _mm_unpacklo_epi64(tr1_6, tr1_7);
  out[7] = _mm_unpackhi_epi64(tr1_6, tr1_7);
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
