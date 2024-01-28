/*
 *  Copyright (c) 2014 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#ifndef VPX_VPX_DSP_ARM_VPX_CONVOLVE8_NEON_H_
#define VPX_VPX_DSP_ARM_VPX_CONVOLVE8_NEON_H_

#include <arm_neon.h>

#include "./vpx_config.h"
#include "./vpx_dsp_rtcd.h"
#include "vpx_dsp/vpx_filter.h"

static INLINE int16x4_t convolve4_4(const int16x4_t s0, const int16x4_t s1,
                                    const int16x4_t s2, const int16x4_t s3,
                                    const int16x4_t filters) {
  int16x4_t sum = vmul_lane_s16(s0, filters, 0);
  sum = vmla_lane_s16(sum, s1, filters, 1);
  sum = vmla_lane_s16(sum, s2, filters, 2);
  sum = vmla_lane_s16(sum, s3, filters, 3);
  return sum;
}

static INLINE uint8x8_t convolve4_8(const int16x8_t s0, const int16x8_t s1,
                                    const int16x8_t s2, const int16x8_t s3,
                                    const int16x4_t filters) {
  int16x8_t sum = vmulq_lane_s16(s0, filters, 0);
  sum = vmlaq_lane_s16(sum, s1, filters, 1);
  sum = vmlaq_lane_s16(sum, s2, filters, 2);
  sum = vmlaq_lane_s16(sum, s3, filters, 3);
  /* We halved the filter values so -1 from right shift. */
  return vqrshrun_n_s16(sum, FILTER_BITS - 1);
}

static INLINE int16x4_t convolve8_4(const int16x4_t s0, const int16x4_t s1,
                                    const int16x4_t s2, const int16x4_t s3,
                                    const int16x4_t s4, const int16x4_t s5,
                                    const int16x4_t s6, const int16x4_t s7,
                                    const int16x8_t filters) {
  const int16x4_t filters_lo = vget_low_s16(filters);
  const int16x4_t filters_hi = vget_high_s16(filters);
  int16x4_t sum;

  sum = vmul_lane_s16(s0, filters_lo, 0);
  sum = vmla_lane_s16(sum, s1, filters_lo, 1);
  sum = vmla_lane_s16(sum, s2, filters_lo, 2);
  sum = vmla_lane_s16(sum, s5, filters_hi, 1);
  sum = vmla_lane_s16(sum, s6, filters_hi, 2);
  sum = vmla_lane_s16(sum, s7, filters_hi, 3);
  sum = vqadd_s16(sum, vmul_lane_s16(s3, filters_lo, 3));
  sum = vqadd_s16(sum, vmul_lane_s16(s4, filters_hi, 0));
  return sum;
}

static INLINE uint8x8_t convolve8_8(const int16x8_t s0, const int16x8_t s1,
                                    const int16x8_t s2, const int16x8_t s3,
                                    const int16x8_t s4, const int16x8_t s5,
                                    const int16x8_t s6, const int16x8_t s7,
                                    const int16x8_t filters) {
  const int16x4_t filters_lo = vget_low_s16(filters);
  const int16x4_t filters_hi = vget_high_s16(filters);
  int16x8_t sum;

  sum = vmulq_lane_s16(s0, filters_lo, 0);
  sum = vmlaq_lane_s16(sum, s1, filters_lo, 1);
  sum = vmlaq_lane_s16(sum, s2, filters_lo, 2);
  sum = vmlaq_lane_s16(sum, s5, filters_hi, 1);
  sum = vmlaq_lane_s16(sum, s6, filters_hi, 2);
  sum = vmlaq_lane_s16(sum, s7, filters_hi, 3);
  sum = vqaddq_s16(sum, vmulq_lane_s16(s3, filters_lo, 3));
  sum = vqaddq_s16(sum, vmulq_lane_s16(s4, filters_hi, 0));
  return vqrshrun_n_s16(sum, FILTER_BITS);
}

static INLINE uint8x8_t scale_filter_8(const uint8x8_t *const s,
                                       const int16x8_t filters) {
  int16x8_t ss[8];

  ss[0] = vreinterpretq_s16_u16(vmovl_u8(s[0]));
  ss[1] = vreinterpretq_s16_u16(vmovl_u8(s[1]));
  ss[2] = vreinterpretq_s16_u16(vmovl_u8(s[2]));
  ss[3] = vreinterpretq_s16_u16(vmovl_u8(s[3]));
  ss[4] = vreinterpretq_s16_u16(vmovl_u8(s[4]));
  ss[5] = vreinterpretq_s16_u16(vmovl_u8(s[5]));
  ss[6] = vreinterpretq_s16_u16(vmovl_u8(s[6]));
  ss[7] = vreinterpretq_s16_u16(vmovl_u8(s[7]));

  return convolve8_8(ss[0], ss[1], ss[2], ss[3], ss[4], ss[5], ss[6], ss[7],
                     filters);
}

#endif  // VPX_VPX_DSP_ARM_VPX_CONVOLVE8_NEON_H_
