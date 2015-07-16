/*
 *  Copyright (c) 2015 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include <arm_neon.h>
#include "./vp9_rtcd.h"
#include "./vpx_config.h"

#include "vpx/vpx_integer.h"

static INLINE unsigned int horizontal_add_u16x8(const uint16x8_t v_16x8) {
  const uint32x4_t a = vpaddlq_u16(v_16x8);
  const uint64x2_t b = vpaddlq_u32(a);
  const uint32x2_t c = vadd_u32(vreinterpret_u32_u64(vget_low_u64(b)),
                                vreinterpret_u32_u64(vget_high_u64(b)));
  return vget_lane_u32(c, 0);
}

unsigned int vp9_avg_8x8_neon(const uint8_t *s, int p) {
  uint8x8_t v_s0 = vld1_u8(s);
  const uint8x8_t v_s1 = vld1_u8(s + p);
  uint16x8_t v_sum = vaddl_u8(v_s0, v_s1);

  v_s0 = vld1_u8(s + 2 * p);
  v_sum = vaddw_u8(v_sum, v_s0);

  v_s0 = vld1_u8(s + 3 * p);
  v_sum = vaddw_u8(v_sum, v_s0);

  v_s0 = vld1_u8(s + 4 * p);
  v_sum = vaddw_u8(v_sum, v_s0);

  v_s0 = vld1_u8(s + 5 * p);
  v_sum = vaddw_u8(v_sum, v_s0);

  v_s0 = vld1_u8(s + 6 * p);
  v_sum = vaddw_u8(v_sum, v_s0);

  v_s0 = vld1_u8(s + 7 * p);
  v_sum = vaddw_u8(v_sum, v_s0);

  return (horizontal_add_u16x8(v_sum) + 32) >> 6;
}

void vp9_int_pro_row_neon(int16_t hbuf[16], uint8_t const *ref,
                          const int ref_stride, const int height) {
  int i;
  uint16x8_t vec_sum_lo = vdupq_n_u16(0);
  uint16x8_t vec_sum_hi = vdupq_n_u16(0);
  const int shift_factor = ((height >> 5) + 3) * -1;
  const int16x8_t vec_shift = vdupq_n_s16(shift_factor);

  for (i = 0; i < height; i += 8) {
    const uint8x16_t vec_row1 = vld1q_u8(ref);
    const uint8x16_t vec_row2 = vld1q_u8(ref + ref_stride);
    const uint8x16_t vec_row3 = vld1q_u8(ref + ref_stride * 2);
    const uint8x16_t vec_row4 = vld1q_u8(ref + ref_stride * 3);
    const uint8x16_t vec_row5 = vld1q_u8(ref + ref_stride * 4);
    const uint8x16_t vec_row6 = vld1q_u8(ref + ref_stride * 5);
    const uint8x16_t vec_row7 = vld1q_u8(ref + ref_stride * 6);
    const uint8x16_t vec_row8 = vld1q_u8(ref + ref_stride * 7);

    vec_sum_lo = vaddw_u8(vec_sum_lo, vget_low_u8(vec_row1));
    vec_sum_hi = vaddw_u8(vec_sum_hi, vget_high_u8(vec_row1));

    vec_sum_lo = vaddw_u8(vec_sum_lo, vget_low_u8(vec_row2));
    vec_sum_hi = vaddw_u8(vec_sum_hi, vget_high_u8(vec_row2));

    vec_sum_lo = vaddw_u8(vec_sum_lo, vget_low_u8(vec_row3));
    vec_sum_hi = vaddw_u8(vec_sum_hi, vget_high_u8(vec_row3));

    vec_sum_lo = vaddw_u8(vec_sum_lo, vget_low_u8(vec_row4));
    vec_sum_hi = vaddw_u8(vec_sum_hi, vget_high_u8(vec_row4));

    vec_sum_lo = vaddw_u8(vec_sum_lo, vget_low_u8(vec_row5));
    vec_sum_hi = vaddw_u8(vec_sum_hi, vget_high_u8(vec_row5));

    vec_sum_lo = vaddw_u8(vec_sum_lo, vget_low_u8(vec_row6));
    vec_sum_hi = vaddw_u8(vec_sum_hi, vget_high_u8(vec_row6));

    vec_sum_lo = vaddw_u8(vec_sum_lo, vget_low_u8(vec_row7));
    vec_sum_hi = vaddw_u8(vec_sum_hi, vget_high_u8(vec_row7));

    vec_sum_lo = vaddw_u8(vec_sum_lo, vget_low_u8(vec_row8));
    vec_sum_hi = vaddw_u8(vec_sum_hi, vget_high_u8(vec_row8));

    ref += ref_stride * 8;
  }

  vec_sum_lo = vshlq_u16(vec_sum_lo, vec_shift);
  vec_sum_hi = vshlq_u16(vec_sum_hi, vec_shift);

  vst1q_s16(hbuf, vreinterpretq_s16_u16(vec_sum_lo));
  hbuf += 8;
  vst1q_s16(hbuf, vreinterpretq_s16_u16(vec_sum_hi));
}

int16_t vp9_int_pro_col_neon(uint8_t const *ref, const int width) {
  int i;
  uint16x8_t vec_sum = vdupq_n_u16(0);

  for (i = 0; i < width; i += 16) {
    const uint8x16_t vec_row = vld1q_u8(ref);
    vec_sum = vaddw_u8(vec_sum, vget_low_u8(vec_row));
    vec_sum = vaddw_u8(vec_sum, vget_high_u8(vec_row));
    ref += 16;
  }

  return horizontal_add_u16x8(vec_sum);
}
