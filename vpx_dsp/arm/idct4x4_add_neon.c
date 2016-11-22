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
#include <assert.h>

#include "./vpx_dsp_rtcd.h"
#include "vpx_dsp/arm/idct_neon.h"
#include "vpx_dsp/arm/transpose_neon.h"
#include "vpx_dsp/txfm_common.h"

static INLINE void idct4x4_16_kernel(const int16x4_t cospis, int16x8_t *a0,
                                     int16x8_t *a1) {
  int16x4_t b0, b1, b2, b3;
  int32x4_t c0, c1, c2, c3;
  int16x8_t d0, d1;

  transpose_s16_4x4q(a0, a1);
  b0 = vget_low_s16(*a0);
  b1 = vget_high_s16(*a0);
  b2 = vget_low_s16(*a1);
  b3 = vget_high_s16(*a1);
  c0 = vmull_lane_s16(b0, cospis, 2);
  c2 = vmull_lane_s16(b1, cospis, 2);
  c1 = vsubq_s32(c0, c2);
  c0 = vaddq_s32(c0, c2);
  c2 = vmull_lane_s16(b2, cospis, 3);
  c3 = vmull_lane_s16(b2, cospis, 1);
  c2 = vmlsl_lane_s16(c2, b3, cospis, 1);
  c3 = vmlal_lane_s16(c3, b3, cospis, 3);
  b0 = vrshrn_n_s32(c0, 14);
  b1 = vrshrn_n_s32(c1, 14);
  b2 = vrshrn_n_s32(c2, 14);
  b3 = vrshrn_n_s32(c3, 14);
  d0 = vcombine_s16(b0, b1);
  d1 = vcombine_s16(b3, b2);
  *a0 = vaddq_s16(d0, d1);
  *a1 = vsubq_s16(d0, d1);
}

void vpx_idct4x4_16_add_neon(const tran_low_t *input, uint8_t *dest,
                             int dest_stride) {
  DECLARE_ALIGNED(16, static const int16_t, cospi[4]) = {
    0, (int16_t)cospi_8_64, (int16_t)cospi_16_64, (int16_t)cospi_24_64
  };
  const uint8_t *dst = dest;
  const int16x4_t cospis = vld1_s16(cospi);
  uint32x2_t dest01_u32 = vdup_n_u32(0);
  uint32x2_t dest32_u32 = vdup_n_u32(0);
  int16x8_t a0, a1;
  uint8x8_t d01, d32;
  uint16x8_t d01_u16, d32_u16;

  assert(!((intptr_t)dest % sizeof(uint32_t)));
  assert(!(dest_stride % sizeof(uint32_t)));

  // Rows
  a0 = load_tran_low_to_s16(input);
  a1 = load_tran_low_to_s16(input + 8);
  idct4x4_16_kernel(cospis, &a0, &a1);

  // Columns
  a1 = vcombine_s16(vget_high_s16(a1), vget_low_s16(a1));
  idct4x4_16_kernel(cospis, &a0, &a1);
  a0 = vrshrq_n_s16(a0, 4);
  a1 = vrshrq_n_s16(a1, 4);

  dest01_u32 = vld1_lane_u32((const uint32_t *)dst, dest01_u32, 0);
  dst += dest_stride;
  dest01_u32 = vld1_lane_u32((const uint32_t *)dst, dest01_u32, 1);
  dst += dest_stride;
  dest32_u32 = vld1_lane_u32((const uint32_t *)dst, dest32_u32, 1);
  dst += dest_stride;
  dest32_u32 = vld1_lane_u32((const uint32_t *)dst, dest32_u32, 0);

  d01_u16 =
      vaddw_u8(vreinterpretq_u16_s16(a0), vreinterpret_u8_u32(dest01_u32));
  d32_u16 =
      vaddw_u8(vreinterpretq_u16_s16(a1), vreinterpret_u8_u32(dest32_u32));
  d01 = vqmovun_s16(vreinterpretq_s16_u16(d01_u16));
  d32 = vqmovun_s16(vreinterpretq_s16_u16(d32_u16));

  vst1_lane_u32((uint32_t *)dest, vreinterpret_u32_u8(d01), 0);
  dest += dest_stride;
  vst1_lane_u32((uint32_t *)dest, vreinterpret_u32_u8(d01), 1);
  dest += dest_stride;
  vst1_lane_u32((uint32_t *)dest, vreinterpret_u32_u8(d32), 1);
  dest += dest_stride;
  vst1_lane_u32((uint32_t *)dest, vreinterpret_u32_u8(d32), 0);
}
