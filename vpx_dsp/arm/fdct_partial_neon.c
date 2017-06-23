/*
 *  Copyright (c) 2017 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include <arm_neon.h>

#include "./vpx_dsp_rtcd.h"
#include "./vpx_config.h"

void vpx_fdct8x8_1_neon(const int16_t *input, tran_low_t *output, int stride) {
  int r;
  int16x8_t sum = vld1q_s16(&input[0]);
  for (r = 1; r < 8; ++r) {
    const int16x8_t input_00 = vld1q_s16(&input[r * stride]);
    sum = vaddq_s16(sum, input_00);
  }
  {
    const int32x4_t a = vpaddlq_s16(sum);
    const int64x2_t b = vpaddlq_s32(a);
    const int32x2_t c = vadd_s32(vreinterpret_s32_s64(vget_low_s64(b)),
                                 vreinterpret_s32_s64(vget_high_s64(b)));
#if CONFIG_VP9_HIGHBITDEPTH
    output[0] = vget_lane_s32(c, 0);
#else
    output[0] = vget_lane_s16(vreinterpret_s16_s32(c), 0);
#endif
    output[1] = 0;
  }
}
