/*
 *  Copyright (c) 2017 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#ifndef VPX_DSP_ARM_MEM_NEON_H_
#define VPX_DSP_ARM_MEM_NEON_H_

#include <arm_neon.h>
#include <assert.h>
#include <string.h>

#include "./vpx_config.h"
#include "vpx/vpx_integer.h"
#include "vpx_dsp/vpx_dsp_common.h"

// Helper functions used to load tran_low_t into int16, narrowing if necessary.
static INLINE int16x8x2_t load_tran_low_to_s16x2q(const tran_low_t *buf) {
#if CONFIG_VP9_HIGHBITDEPTH
  const int32x4x2_t v0 = vld2q_s32(buf);
  const int32x4x2_t v1 = vld2q_s32(buf + 8);
  const int16x4_t s0 = vmovn_s32(v0.val[0]);
  const int16x4_t s1 = vmovn_s32(v0.val[1]);
  const int16x4_t s2 = vmovn_s32(v1.val[0]);
  const int16x4_t s3 = vmovn_s32(v1.val[1]);
  int16x8x2_t res;
  res.val[0] = vcombine_s16(s0, s2);
  res.val[1] = vcombine_s16(s1, s3);
  return res;
#else
  return vld2q_s16(buf);
#endif
}

static INLINE int16x8_t load_tran_low_to_s16q(const tran_low_t *buf) {
#if CONFIG_VP9_HIGHBITDEPTH
  const int32x4_t v0 = vld1q_s32(buf);
  const int32x4_t v1 = vld1q_s32(buf + 4);
  const int16x4_t s0 = vmovn_s32(v0);
  const int16x4_t s1 = vmovn_s32(v1);
  return vcombine_s16(s0, s1);
#else
  return vld1q_s16(buf);
#endif
}

static INLINE int16x4_t load_tran_low_to_s16d(const tran_low_t *buf) {
#if CONFIG_VP9_HIGHBITDEPTH
  const int32x4_t v0 = vld1q_s32(buf);
  return vmovn_s32(v0);
#else
  return vld1_s16(buf);
#endif
}

static INLINE void store_s16q_to_tran_low(tran_low_t *buf, const int16x8_t a) {
#if CONFIG_VP9_HIGHBITDEPTH
  const int32x4_t v0 = vmovl_s16(vget_low_s16(a));
  const int32x4_t v1 = vmovl_s16(vget_high_s16(a));
  vst1q_s32(buf, v0);
  vst1q_s32(buf + 4, v1);
#else
  vst1q_s16(buf, a);
#endif
}
#endif  // VPX_DSP_ARM_MEM_NEON_H_
