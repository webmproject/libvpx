/*
 *  Copyright (c) 2018 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#ifndef VP9_COMMON_ARM_NEON_VP9_IHT_NEON_H_
#define VP9_COMMON_ARM_NEON_VP9_IHT_NEON_H_

#include <arm_neon.h>

#include "./vp9_rtcd.h"
#include "./vpx_config.h"
#include "vp9/common/vp9_common.h"
#include "vpx_dsp/arm/idct_neon.h"
#include "vpx_dsp/arm/mem_neon.h"
#include "vpx_dsp/txfm_common.h"

static INLINE void iadst4(int16x8_t *const io) {
  const int32x4_t c3 = vdupq_n_s32(sinpi_3_9);
  int16x4_t x[4];
  int32x4_t s[8], output[4];
  const int16x4_t c =
      create_s16x4_neon(sinpi_1_9, sinpi_2_9, sinpi_3_9, sinpi_4_9);

  x[0] = vget_low_s16(io[0]);
  x[1] = vget_low_s16(io[1]);
  x[2] = vget_high_s16(io[0]);
  x[3] = vget_high_s16(io[1]);

  s[0] = vmull_lane_s16(x[0], c, 0);
  s[1] = vmull_lane_s16(x[0], c, 1);
  s[2] = vmull_lane_s16(x[1], c, 2);
  s[3] = vmull_lane_s16(x[2], c, 3);
  s[4] = vmull_lane_s16(x[2], c, 0);
  s[5] = vmull_lane_s16(x[3], c, 1);
  s[6] = vmull_lane_s16(x[3], c, 3);
  s[7] = vaddl_s16(x[0], x[3]);
  s[7] = vsubw_s16(s[7], x[2]);

  s[0] = vaddq_s32(s[0], s[3]);
  s[0] = vaddq_s32(s[0], s[5]);
  s[1] = vsubq_s32(s[1], s[4]);
  s[1] = vsubq_s32(s[1], s[6]);
  s[3] = s[2];
  s[2] = vmulq_s32(c3, s[7]);

  output[0] = vaddq_s32(s[0], s[3]);
  output[1] = vaddq_s32(s[1], s[3]);
  output[2] = s[2];
  output[3] = vaddq_s32(s[0], s[1]);
  output[3] = vsubq_s32(output[3], s[3]);
  dct_const_round_shift_low_8_dual(output, &io[0], &io[1]);
}

#endif  // VP9_COMMON_ARM_NEON_VP9_IHT_NEON_H_
