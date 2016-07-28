/*
 *  Copyright (c) 2016 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#ifndef VPX_DSP_ARM_TRANSPOSE_NEON_H_
#define VPX_DSP_ARM_TRANSPOSE_NEON_H_

#include <arm_neon.h>

#include "./vpx_config.h"

// Transpose 64 bit elements as follows:
// a0: 00 01 02 03 04 05 06 07
// a1: 16 17 18 19 20 21 22 23
//
// b0.val[0]: 00 01 02 03 16 17 18 19
// b0.val[1]: 04 05 06 07 20 21 22 23
static INLINE int16x8x2_t vpx_vtrnq_s64(int32x4_t a0, int32x4_t a1) {
  int16x8x2_t b0;
  b0.val[0] = vcombine_s16(vreinterpret_s16_s32(vget_low_s32(a0)),
                           vreinterpret_s16_s32(vget_low_s32(a1)));
  b0.val[1] = vcombine_s16(vreinterpret_s16_s32(vget_high_s32(a0)),
                           vreinterpret_s16_s32(vget_high_s32(a1)));
  return b0;
}

static INLINE void transpose_s16_8x8(int16x8_t *a0, int16x8_t *a1,
                                     int16x8_t *a2, int16x8_t *a3,
                                     int16x8_t *a4, int16x8_t *a5,
                                     int16x8_t *a6, int16x8_t *a7) {
  // Swap 16 bit elements. Goes from:
  // a0: 00 01 02 03 04 05 06 07
  // a1: 08 09 10 11 12 13 14 15
  // a2: 16 17 18 19 20 21 22 23
  // a3: 24 25 26 27 28 29 30 31
  // a4: 32 33 34 35 36 37 38 39
  // a5: 40 41 42 43 44 45 46 47
  // a6: 48 49 50 51 52 53 54 55
  // a7: 56 57 58 59 60 61 62 63
  // to:
  // b0.val[0]: 00 08 02 10 04 12 06 14
  // b0.val[1]: 01 09 03 11 05 13 07 15
  // b1.val[0]: 16 24 18 26 20 28 22 30
  // b1.val[1]: 17 25 19 27 21 29 23 31
  // b2.val[0]: 32 40 34 42 36 44 38 46
  // b2.val[1]: 33 41 35 43 37 45 39 47
  // b3.val[0]: 48 56 50 58 52 60 54 62
  // b3.val[1]: 49 57 51 59 53 61 55 63

  const int16x8x2_t b0 = vtrnq_s16(*a0, *a1);
  const int16x8x2_t b1 = vtrnq_s16(*a2, *a3);
  const int16x8x2_t b2 = vtrnq_s16(*a4, *a5);
  const int16x8x2_t b3 = vtrnq_s16(*a6, *a7);

  // Swap 32 bit elements resulting in:
  // c0.val[0]: 00 08 16 24 04 12 20 28
  // c0.val[1]: 02 10 18 26 06 14 22 30
  // c1.val[0]: 01 09 17 25 05 13 21 29
  // c1.val[1]: 03 11 19 27 07 15 23 31
  // c2.val[0]: 32 40 48 56 36 44 52 60
  // c2.val[1]: 34 42 50 58 38 46 54 62
  // c3.val[0]: 33 41 49 57 37 45 53 61
  // c3.val[1]: 35 43 51 59 39 47 55 63

  const int32x4x2_t c0 = vtrnq_s32(vreinterpretq_s32_s16(b0.val[0]),
                                   vreinterpretq_s32_s16(b1.val[0]));
  const int32x4x2_t c1 = vtrnq_s32(vreinterpretq_s32_s16(b0.val[1]),
                                   vreinterpretq_s32_s16(b1.val[1]));
  const int32x4x2_t c2 = vtrnq_s32(vreinterpretq_s32_s16(b2.val[0]),
                                   vreinterpretq_s32_s16(b3.val[0]));
  const int32x4x2_t c3 = vtrnq_s32(vreinterpretq_s32_s16(b2.val[1]),
                                   vreinterpretq_s32_s16(b3.val[1]));

  // Swap 64 bit elements resulting in:
  // d0.val[0]: 00 08 16 24 32 40 48 56
  // d0.val[1]: 04 12 20 28 36 44 52 60
  // d1.val[0]: 01 09 17 25 33 41 49 57
  // d1.val[1]: 05 13 21 29 37 45 53 61
  // d2.val[0]: 02 10 18 26 34 42 50 58
  // d2.val[1]: 06 14 22 30 38 46 54 62
  // d3.val[0]: 03 11 19 27 35 43 51 59
  // d3.val[1]: 07 15 23 31 39 47 55 63
  const int16x8x2_t d0 = vpx_vtrnq_s64(c0.val[0], c2.val[0]);
  const int16x8x2_t d1 = vpx_vtrnq_s64(c1.val[0], c3.val[0]);
  const int16x8x2_t d2 = vpx_vtrnq_s64(c0.val[1], c2.val[1]);
  const int16x8x2_t d3 = vpx_vtrnq_s64(c1.val[1], c3.val[1]);

  *a0 = d0.val[0];
  *a1 = d1.val[0];
  *a2 = d2.val[0];
  *a3 = d3.val[0];
  *a4 = d0.val[1];
  *a5 = d1.val[1];
  *a6 = d2.val[1];
  *a7 = d3.val[1];
}

#endif  // VPX_DSP_ARM_TRANSPOSE_NEON_H_
