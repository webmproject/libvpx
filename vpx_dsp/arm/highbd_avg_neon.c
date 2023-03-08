/*
 *  Copyright (c) 2023 The WebM project authors. All Rights Reserved.
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

#include "vpx_dsp/arm/mem_neon.h"
#include "vpx_dsp/arm/sum_neon.h"

// coeff: 32 bits, dynamic range [-2147483648, 2147483647].
// length: value range {16, 64, 256, 1024}.
// satd: 42 bits, dynamic range [-2147483648 * 1024, 2147483647 * 1024]
int vpx_highbd_satd_neon(const tran_low_t *coeff, int length) {
  int64x2_t sum_s64[2] = { vdupq_n_s64(0), vdupq_n_s64(0) };

  do {
    int32x4_t abs0, abs1;
    const int32x4_t s0 = load_tran_low_to_s32q(coeff);
    const int32x4_t s1 = load_tran_low_to_s32q(coeff + 4);

    abs0 = vabsq_s32(s0);
    sum_s64[0] = vpadalq_s32(sum_s64[0], abs0);
    abs1 = vabsq_s32(s1);
    sum_s64[1] = vpadalq_s32(sum_s64[1], abs1);

    length -= 8;
    coeff += 8;
  } while (length != 0);

  return (int)horizontal_add_int64x2(vaddq_s64(sum_s64[0], sum_s64[1]));
}
