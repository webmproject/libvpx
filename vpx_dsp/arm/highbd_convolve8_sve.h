/*
 *  Copyright (c) 2024 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#ifndef VPX_VPX_DSP_ARM_HIGHBD_CONVOLVE8_SVE_H_
#define VPX_VPX_DSP_ARM_HIGHBD_CONVOLVE8_SVE_H_

#include <arm_neon.h>

#include "vpx_dsp/arm/vpx_neon_sve_bridge.h"

static INLINE uint16x4_t highbd_convolve8_4(const int16x8_t s[4],
                                            const int16x8_t filter,
                                            const uint16x4_t max) {
  int64x2_t sum[4];

  sum[0] = vpx_dotq_s16(vdupq_n_s64(0), s[0], filter);
  sum[1] = vpx_dotq_s16(vdupq_n_s64(0), s[1], filter);
  sum[2] = vpx_dotq_s16(vdupq_n_s64(0), s[2], filter);
  sum[3] = vpx_dotq_s16(vdupq_n_s64(0), s[3], filter);

  sum[0] = vpaddq_s64(sum[0], sum[1]);
  sum[2] = vpaddq_s64(sum[2], sum[3]);

  int32x4_t res_s32 = vcombine_s32(vmovn_s64(sum[0]), vmovn_s64(sum[2]));

  uint16x4_t res_u16 = vqrshrun_n_s32(res_s32, FILTER_BITS);
  return vmin_u16(res_u16, max);
}

static INLINE uint16x8_t highbd_convolve8_8(const int16x8_t s[8],
                                            const int16x8_t filter,
                                            const uint16x8_t max) {
  int64x2_t sum[8];

  sum[0] = vpx_dotq_s16(vdupq_n_s64(0), s[0], filter);
  sum[1] = vpx_dotq_s16(vdupq_n_s64(0), s[1], filter);
  sum[2] = vpx_dotq_s16(vdupq_n_s64(0), s[2], filter);
  sum[3] = vpx_dotq_s16(vdupq_n_s64(0), s[3], filter);
  sum[4] = vpx_dotq_s16(vdupq_n_s64(0), s[4], filter);
  sum[5] = vpx_dotq_s16(vdupq_n_s64(0), s[5], filter);
  sum[6] = vpx_dotq_s16(vdupq_n_s64(0), s[6], filter);
  sum[7] = vpx_dotq_s16(vdupq_n_s64(0), s[7], filter);

  int64x2_t sum01 = vpaddq_s64(sum[0], sum[1]);
  int64x2_t sum23 = vpaddq_s64(sum[2], sum[3]);
  int64x2_t sum45 = vpaddq_s64(sum[4], sum[5]);
  int64x2_t sum67 = vpaddq_s64(sum[6], sum[7]);

  int32x4_t res0 = vcombine_s32(vmovn_s64(sum01), vmovn_s64(sum23));
  int32x4_t res1 = vcombine_s32(vmovn_s64(sum45), vmovn_s64(sum67));

  uint16x8_t res = vcombine_u16(vqrshrun_n_s32(res0, FILTER_BITS),
                                vqrshrun_n_s32(res1, FILTER_BITS));
  return vminq_u16(res, max);
}

#endif  // VPX_VPX_DSP_ARM_HIGHBD_CONVOLVE8_SVE_H_
