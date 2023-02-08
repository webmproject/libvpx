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

void vpx_highbd_comp_avg_pred_neon(uint16_t *comp_pred, const uint16_t *pred,
                                   int width, int height, const uint16_t *ref,
                                   int ref_stride) {
  int i, j;
  uint32x4_t one_u32 = vdupq_n_u32(1);
  if (width >= 8) {
    for (i = 0; i < height; ++i) {
      for (j = 0; j < width; j += 8) {
        const uint16x8_t pred_u16 = vld1q_u16(&pred[j]);
        const uint16x8_t ref_u16 = vld1q_u16(&ref[j]);
        const uint32x4_t sum1_u32 =
            vaddl_u16(vget_low_u16(pred_u16), vget_low_u16(ref_u16));
        const uint32x4_t sum2_u32 =
            vaddl_u16(vget_high_u16(pred_u16), vget_high_u16(ref_u16));
        const uint16x4_t sum1_u16 =
            vshrn_n_u32(vaddq_u32(sum1_u32, one_u32), 1);
        const uint16x4_t sum2_u16 =
            vshrn_n_u32(vaddq_u32(sum2_u32, one_u32), 1);
        const uint16x8_t vcomp_pred = vcombine_u16(sum1_u16, sum2_u16);
        vst1q_u16(&comp_pred[j], vcomp_pred);
      }
      comp_pred += width;
      pred += width;
      ref += ref_stride;
    }
  } else {
    assert(width >= 4);
    for (i = 0; i < height; ++i) {
      for (j = 0; j < width; j += 4) {
        const uint16x4_t pred_u16 = vld1_u16(&pred[j]);
        const uint16x4_t ref_u16 = vld1_u16(&ref[j]);
        const uint32x4_t sum_u32 = vaddl_u16(pred_u16, ref_u16);
        const uint16x4_t vcomp_pred =
            vshrn_n_u32(vaddq_u32(sum_u32, one_u32), 1);
        vst1_u16(&comp_pred[j], vcomp_pred);
      }
      comp_pred += width;
      pred += width;
      ref += ref_stride;
    }
  }
}
