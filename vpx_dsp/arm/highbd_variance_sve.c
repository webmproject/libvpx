/*
 * Copyright (c) 2024 The WebM project authors. All Rights Reserved.
 *
 * Use of this source code is governed by a BSD-style license
 * that can be found in the LICENSE file in the root of the source
 * tree. An additional intellectual property rights grant can be found
 * in the file PATENTS.  All contributing project authors may
 * be found in the AUTHORS file in the root of the source tree.
 */

#include <arm_neon.h>

#include "./vpx_dsp_rtcd.h"
#include "./vpx_config.h"

#include "vpx_dsp/arm/dot_neon_sve_bridge.h"
#include "vpx_dsp/arm/sum_neon.h"

static INLINE uint32_t highbd_mse_wxh_sve(const uint16_t *src_ptr,
                                          int src_stride,
                                          const uint16_t *ref_ptr,
                                          int ref_stride, int w, int h) {
  uint64x2_t sse = vdupq_n_u64(0);

  do {
    int j = 0;
    do {
      uint16x8_t s = vld1q_u16(src_ptr + j);
      uint16x8_t r = vld1q_u16(ref_ptr + j);

      uint16x8_t diff = vabdq_u16(s, r);

      sse = vpx_dotq_u16(sse, diff, diff);

      j += 8;
    } while (j < w);

    src_ptr += src_stride;
    ref_ptr += ref_stride;
  } while (--h != 0);

  return (uint32_t)horizontal_add_uint64x2(sse);
}

#define HIGHBD_MSE_WXH_SVE(w, h)                                      \
  uint32_t vpx_highbd_10_mse##w##x##h##_sve(                          \
      const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr, \
      int ref_stride, uint32_t *sse) {                                \
    uint16_t *src = CONVERT_TO_SHORTPTR(src_ptr);                     \
    uint16_t *ref = CONVERT_TO_SHORTPTR(ref_ptr);                     \
    uint32_t sse_tmp =                                                \
        highbd_mse_wxh_sve(src, src_stride, ref, ref_stride, w, h);   \
    sse_tmp = ROUND_POWER_OF_TWO(sse_tmp, 4);                         \
    *sse = sse_tmp;                                                   \
    return sse_tmp;                                                   \
  }                                                                   \
                                                                      \
  uint32_t vpx_highbd_12_mse##w##x##h##_sve(                          \
      const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr, \
      int ref_stride, uint32_t *sse) {                                \
    uint16_t *src = CONVERT_TO_SHORTPTR(src_ptr);                     \
    uint16_t *ref = CONVERT_TO_SHORTPTR(ref_ptr);                     \
    uint32_t sse_tmp =                                                \
        highbd_mse_wxh_sve(src, src_stride, ref, ref_stride, w, h);   \
    sse_tmp = ROUND_POWER_OF_TWO(sse_tmp, 8);                         \
    *sse = sse_tmp;                                                   \
    return sse_tmp;                                                   \
  }

HIGHBD_MSE_WXH_SVE(16, 16)
HIGHBD_MSE_WXH_SVE(16, 8)
HIGHBD_MSE_WXH_SVE(8, 16)
HIGHBD_MSE_WXH_SVE(8, 8)

#undef HIGHBD_MSE_WXH_SVE
