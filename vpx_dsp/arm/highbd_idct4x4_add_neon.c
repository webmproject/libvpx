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

static INLINE void highbd_idct4x4_1_add_kernel1(uint16_t **dest,
                                                const int stride,
                                                const int16x8_t res,
                                                const int16x8_t max) {
  const uint16x4_t a0 = vld1_u16(*dest);
  const uint16x4_t a1 = vld1_u16(*dest + stride);
  const int16x8_t a = vreinterpretq_s16_u16(vcombine_u16(a0, a1));
  // Note: In some profile tests, res is quite close to +/-32767.
  // We use saturating addition.
  const int16x8_t b = vqaddq_s16(res, a);
  const int16x8_t c = vminq_s16(b, max);
  const uint16x8_t d = vqshluq_n_s16(c, 0);
  vst1_u16(*dest, vget_low_u16(d));
  *dest += stride;
  vst1_u16(*dest, vget_high_u16(d));
  *dest += stride;
}

// res is in reverse row order
static INLINE void highbd_idct4x4_1_add_kernel2(uint16_t **dest,
                                                const int stride,
                                                const int16x8_t res,
                                                const int16x8_t max) {
  const uint16x4_t a0 = vld1_u16(*dest);
  const uint16x4_t a1 = vld1_u16(*dest + stride);
  const int16x8_t a = vreinterpretq_s16_u16(vcombine_u16(a1, a0));
  // Note: In some profile tests, res is quite close to +/-32767.
  // We use saturating addition.
  const int16x8_t b = vqaddq_s16(res, a);
  const int16x8_t c = vminq_s16(b, max);
  const uint16x8_t d = vqshluq_n_s16(c, 0);
  vst1_u16(*dest, vget_high_u16(d));
  *dest += stride;
  vst1_u16(*dest, vget_low_u16(d));
  *dest += stride;
}

void vpx_highbd_idct4x4_1_add_neon(const tran_low_t *input, uint16_t *dest,
                                   int stride, int bd) {
  const int16x8_t max = vdupq_n_s16((1 << bd) - 1);
  const tran_low_t out0 = HIGHBD_WRAPLOW(
      dct_const_round_shift(input[0] * (tran_high_t)cospi_16_64), bd);
  const tran_low_t out1 = HIGHBD_WRAPLOW(
      dct_const_round_shift(out0 * (tran_high_t)cospi_16_64), bd);
  const int16_t a1 = ROUND_POWER_OF_TWO(out1, 4);
  const int16x8_t dc = vdupq_n_s16(a1);

  highbd_idct4x4_1_add_kernel1(&dest, stride, dc, max);
  highbd_idct4x4_1_add_kernel1(&dest, stride, dc, max);
}

static INLINE void idct4x4_16_kernel_bd10(const int32x4_t cospis,
                                          int32x4_t *const a) {
  int32x4_t b0, b1, b2, b3;

  transpose_s32_4x4(&a[0], &a[1], &a[2], &a[3]);
  b0 = vaddq_s32(a[0], a[2]);
  b1 = vsubq_s32(a[0], a[2]);
  b0 = vmulq_lane_s32(b0, vget_high_s32(cospis), 0);
  b1 = vmulq_lane_s32(b1, vget_high_s32(cospis), 0);
  b2 = vmulq_lane_s32(a[1], vget_high_s32(cospis), 1);
  b3 = vmulq_lane_s32(a[1], vget_low_s32(cospis), 1);
  b2 = vmlsq_lane_s32(b2, a[3], vget_low_s32(cospis), 1);
  b3 = vmlaq_lane_s32(b3, a[3], vget_high_s32(cospis), 1);
  b0 = vrshrq_n_s32(b0, DCT_CONST_BITS);
  b1 = vrshrq_n_s32(b1, DCT_CONST_BITS);
  b2 = vrshrq_n_s32(b2, DCT_CONST_BITS);
  b3 = vrshrq_n_s32(b3, DCT_CONST_BITS);
  a[0] = vaddq_s32(b0, b3);
  a[1] = vaddq_s32(b1, b2);
  a[2] = vsubq_s32(b1, b2);
  a[3] = vsubq_s32(b0, b3);
}

static INLINE void idct4x4_16_kernel_bd12(const int32x4_t cospis,
                                          int32x4_t *const a) {
  int32x4_t b0, b1, b2, b3;
  int64x2_t c[12];

  transpose_s32_4x4(&a[0], &a[1], &a[2], &a[3]);
  b0 = vaddq_s32(a[0], a[2]);
  b1 = vsubq_s32(a[0], a[2]);
  c[0] = vmull_lane_s32(vget_low_s32(b0), vget_high_s32(cospis), 0);
  c[1] = vmull_lane_s32(vget_high_s32(b0), vget_high_s32(cospis), 0);
  c[2] = vmull_lane_s32(vget_low_s32(b1), vget_high_s32(cospis), 0);
  c[3] = vmull_lane_s32(vget_high_s32(b1), vget_high_s32(cospis), 0);
  c[4] = vmull_lane_s32(vget_low_s32(a[1]), vget_high_s32(cospis), 1);
  c[5] = vmull_lane_s32(vget_high_s32(a[1]), vget_high_s32(cospis), 1);
  c[6] = vmull_lane_s32(vget_low_s32(a[1]), vget_low_s32(cospis), 1);
  c[7] = vmull_lane_s32(vget_high_s32(a[1]), vget_low_s32(cospis), 1);
  c[8] = vmull_lane_s32(vget_low_s32(a[3]), vget_low_s32(cospis), 1);
  c[9] = vmull_lane_s32(vget_high_s32(a[3]), vget_low_s32(cospis), 1);
  c[10] = vmull_lane_s32(vget_low_s32(a[3]), vget_high_s32(cospis), 1);
  c[11] = vmull_lane_s32(vget_high_s32(a[3]), vget_high_s32(cospis), 1);
  c[4] = vsubq_s64(c[4], c[8]);
  c[5] = vsubq_s64(c[5], c[9]);
  c[6] = vaddq_s64(c[6], c[10]);
  c[7] = vaddq_s64(c[7], c[11]);
  b0 = vcombine_s32(vrshrn_n_s64(c[0], DCT_CONST_BITS),
                    vrshrn_n_s64(c[1], DCT_CONST_BITS));
  b1 = vcombine_s32(vrshrn_n_s64(c[2], DCT_CONST_BITS),
                    vrshrn_n_s64(c[3], DCT_CONST_BITS));
  b2 = vcombine_s32(vrshrn_n_s64(c[4], DCT_CONST_BITS),
                    vrshrn_n_s64(c[5], DCT_CONST_BITS));
  b3 = vcombine_s32(vrshrn_n_s64(c[6], DCT_CONST_BITS),
                    vrshrn_n_s64(c[7], DCT_CONST_BITS));
  a[0] = vaddq_s32(b0, b3);
  a[1] = vaddq_s32(b1, b2);
  a[2] = vsubq_s32(b1, b2);
  a[3] = vsubq_s32(b0, b3);
}

void vpx_highbd_idct4x4_16_add_neon(const tran_low_t *input, uint16_t *dest,
                                    int stride, int bd) {
  const int16x8_t max = vdupq_n_s16((1 << bd) - 1);
  int16x8_t a[2];
  int32x4_t c[4];

  c[0] = vld1q_s32(input);
  c[1] = vld1q_s32(input + 4);
  c[2] = vld1q_s32(input + 8);
  c[3] = vld1q_s32(input + 12);

  if (bd == 8) {
    const int16x4_t cospis = vld1_s16(kCospi);

    // Rows
    a[0] = vcombine_s16(vmovn_s32(c[0]), vmovn_s32(c[1]));
    a[1] = vcombine_s16(vmovn_s32(c[2]), vmovn_s32(c[3]));
    idct4x4_16_kernel_bd8(cospis, a);

    // Columns
    a[1] = vcombine_s16(vget_high_s16(a[1]), vget_low_s16(a[1]));
    idct4x4_16_kernel_bd8(cospis, a);
    a[0] = vrshrq_n_s16(a[0], 4);
    a[1] = vrshrq_n_s16(a[1], 4);
  } else {
    const int32x4_t cospis = vld1q_s32(kCospi32);

    if (bd == 10) {
      idct4x4_16_kernel_bd10(cospis, c);
      idct4x4_16_kernel_bd10(cospis, c);
    } else {
      idct4x4_16_kernel_bd12(cospis, c);
      idct4x4_16_kernel_bd12(cospis, c);
    }
    a[0] = vcombine_s16(vqrshrn_n_s32(c[0], 4), vqrshrn_n_s32(c[1], 4));
    a[1] = vcombine_s16(vqrshrn_n_s32(c[3], 4), vqrshrn_n_s32(c[2], 4));
  }

  highbd_idct4x4_1_add_kernel1(&dest, stride, a[0], max);
  highbd_idct4x4_1_add_kernel2(&dest, stride, a[1], max);
}
