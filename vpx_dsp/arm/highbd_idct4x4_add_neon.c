/*
 *  Copyright (c) 2016 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include <arm_neon.h>

#include "./vpx_dsp_rtcd.h"
#include "vpx_dsp/arm/idct_neon.h"
#include "vpx_dsp/inv_txfm.h"

void vpx_highbd_idct4x4_1_add_neon(const tran_low_t *input, uint8_t *dest8,
                                   int stride, int bd) {
  int i;
  const int16x8_t max = vdupq_n_s16((1 << bd) - 1);
  const tran_low_t out0 = dct_const_round_shift(input[0] * cospi_16_64);
  const tran_low_t out1 = dct_const_round_shift(out0 * cospi_16_64);
  const int16_t a1 = ROUND_POWER_OF_TWO(out1, 4);
  const int16x8_t dc = vdupq_n_s16(a1);
  uint16_t *dest = CONVERT_TO_SHORTPTR(dest8);
  int16x8_t a;
  uint16x8_t b;
  uint16x4_t d0, d1;

  for (i = 0; i < 2; i++) {
    d0 = vld1_u16(dest);
    d1 = vld1_u16(dest + stride);
    a = vreinterpretq_s16_u16(vcombine_u16(d0, d1));
    a = vaddq_s16(dc, a);
    a = vminq_s16(a, max);
    b = vqshluq_n_s16(a, 0);
    vst1_u16(dest, vget_low_u16(b));
    dest += stride;
    vst1_u16(dest, vget_high_u16(b));
    dest += stride;
  }
}

static INLINE void idct4x4_16_kernel_bd10(const int32x4_t cospis,
                                          int32x4_t *const a0,
                                          int32x4_t *const a1,
                                          int32x4_t *const a2,
                                          int32x4_t *const a3) {
  int32x4_t b0, b1, b2, b3;

  transpose_s32_4x4(a0, a1, a2, a3);
  b0 = vaddq_s32(*a0, *a2);
  b1 = vsubq_s32(*a0, *a2);
  b0 = vmulq_lane_s32(b0, vget_high_s32(cospis), 0);
  b1 = vmulq_lane_s32(b1, vget_high_s32(cospis), 0);
  b2 = vmulq_lane_s32(*a1, vget_high_s32(cospis), 1);
  b3 = vmulq_lane_s32(*a1, vget_low_s32(cospis), 1);
  b2 = vmlsq_lane_s32(b2, *a3, vget_low_s32(cospis), 1);
  b3 = vmlaq_lane_s32(b3, *a3, vget_high_s32(cospis), 1);
  b0 = vrshrq_n_s32(b0, 14);
  b1 = vrshrq_n_s32(b1, 14);
  b2 = vrshrq_n_s32(b2, 14);
  b3 = vrshrq_n_s32(b3, 14);
  *a0 = vaddq_s32(b0, b3);
  *a1 = vaddq_s32(b1, b2);
  *a2 = vsubq_s32(b1, b2);
  *a3 = vsubq_s32(b0, b3);
}

static INLINE void idct4x4_16_kernel_bd12(const int32x4_t cospis,
                                          int32x4_t *const a0,
                                          int32x4_t *const a1,
                                          int32x4_t *const a2,
                                          int32x4_t *const a3) {
  int32x4_t b0, b1, b2, b3;
  int64x2_t c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11;

  transpose_s32_4x4(a0, a1, a2, a3);
  b0 = vaddq_s32(*a0, *a2);
  b1 = vsubq_s32(*a0, *a2);
  c0 = vmull_lane_s32(vget_low_s32(b0), vget_high_s32(cospis), 0);
  c1 = vmull_lane_s32(vget_high_s32(b0), vget_high_s32(cospis), 0);
  c2 = vmull_lane_s32(vget_low_s32(b1), vget_high_s32(cospis), 0);
  c3 = vmull_lane_s32(vget_high_s32(b1), vget_high_s32(cospis), 0);
  c4 = vmull_lane_s32(vget_low_s32(*a1), vget_high_s32(cospis), 1);
  c5 = vmull_lane_s32(vget_high_s32(*a1), vget_high_s32(cospis), 1);
  c6 = vmull_lane_s32(vget_low_s32(*a1), vget_low_s32(cospis), 1);
  c7 = vmull_lane_s32(vget_high_s32(*a1), vget_low_s32(cospis), 1);
  c8 = vmull_lane_s32(vget_low_s32(*a3), vget_low_s32(cospis), 1);
  c9 = vmull_lane_s32(vget_high_s32(*a3), vget_low_s32(cospis), 1);
  c10 = vmull_lane_s32(vget_low_s32(*a3), vget_high_s32(cospis), 1);
  c11 = vmull_lane_s32(vget_high_s32(*a3), vget_high_s32(cospis), 1);
  c4 = vsubq_s64(c4, c8);
  c5 = vsubq_s64(c5, c9);
  c6 = vaddq_s64(c6, c10);
  c7 = vaddq_s64(c7, c11);
  b0 = vcombine_s32(vrshrn_n_s64(c0, 14), vrshrn_n_s64(c1, 14));
  b1 = vcombine_s32(vrshrn_n_s64(c2, 14), vrshrn_n_s64(c3, 14));
  b2 = vcombine_s32(vrshrn_n_s64(c4, 14), vrshrn_n_s64(c5, 14));
  b3 = vcombine_s32(vrshrn_n_s64(c6, 14), vrshrn_n_s64(c7, 14));
  *a0 = vaddq_s32(b0, b3);
  *a1 = vaddq_s32(b1, b2);
  *a2 = vsubq_s32(b1, b2);
  *a3 = vsubq_s32(b0, b3);
}

void vpx_highbd_idct4x4_16_add_neon(const tran_low_t *input, uint8_t *dest8,
                                    int stride, int bd) {
  DECLARE_ALIGNED(16, static const int32_t, kCospi32[4]) = { 0, 15137, 11585,
                                                             6270 };
  const int16x8_t max = vdupq_n_s16((1 << bd) - 1);
  int32x4_t c0 = vld1q_s32(input);
  int32x4_t c1 = vld1q_s32(input + 4);
  int32x4_t c2 = vld1q_s32(input + 8);
  int32x4_t c3 = vld1q_s32(input + 12);
  uint16_t *dest = CONVERT_TO_SHORTPTR(dest8);
  const uint16_t *dst = dest;
  int16x8_t a0, a1, d01, d32;
  int16x4_t d0, d1, d2, d3;
  uint16x8_t d01_u16, d32_u16;

  if (bd == 8) {
    const int16x4_t cospis = vld1_s16(kCospi);

    // Rows
    a0 = vcombine_s16(vmovn_s32(c0), vmovn_s32(c1));
    a1 = vcombine_s16(vmovn_s32(c2), vmovn_s32(c3));
    idct4x4_16_kernel_bd8(cospis, &a0, &a1);

    // Columns
    a1 = vcombine_s16(vget_high_s16(a1), vget_low_s16(a1));
    idct4x4_16_kernel_bd8(cospis, &a0, &a1);
    a0 = vrshrq_n_s16(a0, 4);
    a1 = vrshrq_n_s16(a1, 4);
  } else {
    const int32x4_t cospis = vld1q_s32(kCospi32);

    if (bd == 10) {
      idct4x4_16_kernel_bd10(cospis, &c0, &c1, &c2, &c3);
      idct4x4_16_kernel_bd10(cospis, &c0, &c1, &c2, &c3);
    } else {
      idct4x4_16_kernel_bd12(cospis, &c0, &c1, &c2, &c3);
      idct4x4_16_kernel_bd12(cospis, &c0, &c1, &c2, &c3);
    }
    // Note: In some profile tests, a0 and a1 are quite close to +/-32767.
    // We use saturating narrow shift in case they could be even larger.
    a0 = vcombine_s16(vqrshrn_n_s32(c0, 4), vqrshrn_n_s32(c1, 4));
    a1 = vcombine_s16(vqrshrn_n_s32(c3, 4), vqrshrn_n_s32(c2, 4));
  }

  d0 = vreinterpret_s16_u16(vld1_u16(dst));
  dst += stride;
  d1 = vreinterpret_s16_u16(vld1_u16(dst));
  dst += stride;
  d2 = vreinterpret_s16_u16(vld1_u16(dst));
  dst += stride;
  d3 = vreinterpret_s16_u16(vld1_u16(dst));
  d01 = vcombine_s16(d0, d1);
  d32 = vcombine_s16(d3, d2);

  // Note: In some profile tests, a0 and a1 is quite close to +/-32767.
  // We use saturating addition.
  d01 = vqaddq_s16(a0, d01);
  d32 = vqaddq_s16(a1, d32);
  d01 = vminq_s16(d01, max);
  d32 = vminq_s16(d32, max);
  d01_u16 = vqshluq_n_s16(d01, 0);
  d32_u16 = vqshluq_n_s16(d32, 0);

  vst1_u16(dest, vget_low_u16(d01_u16));
  dest += stride;
  vst1_u16(dest, vget_high_u16(d01_u16));
  dest += stride;
  vst1_u16(dest, vget_high_u16(d32_u16));
  dest += stride;
  vst1_u16(dest, vget_low_u16(d32_u16));
}
