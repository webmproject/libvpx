/*
 *  Copyright (c) 2014 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include <arm_neon.h>
#include "./vp9_rtcd.h"
#include "./vpx_dsp_rtcd.h"
#include "./vpx_config.h"

#include "vp9/common/vp9_blockd.h"
#include "vp9/common/vp9_idct.h"

void vp9_fdct8x8_1_neon(const int16_t *input, int16_t *output, int stride) {
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
    output[0] = vget_lane_s16(vreinterpret_s16_s32(c), 0);
    output[1] = 0;
  }
}

void vp9_fdct8x8_quant_neon(const int16_t *input, int stride,
                            int16_t* coeff_ptr, intptr_t n_coeffs,
                            int skip_block, const int16_t* zbin_ptr,
                            const int16_t* round_ptr, const int16_t* quant_ptr,
                            const int16_t* quant_shift_ptr,
                            int16_t* qcoeff_ptr, int16_t* dqcoeff_ptr,
                            const int16_t* dequant_ptr, uint16_t* eob_ptr,
                            const int16_t* scan_ptr,
                            const int16_t* iscan_ptr) {
  int16_t temp_buffer[64];
  (void)coeff_ptr;

  vp9_fdct8x8_neon(input, temp_buffer, stride);
  vp9_quantize_fp_neon(temp_buffer, n_coeffs, skip_block, zbin_ptr, round_ptr,
                       quant_ptr, quant_shift_ptr, qcoeff_ptr, dqcoeff_ptr,
                       dequant_ptr, eob_ptr, scan_ptr, iscan_ptr);
}
