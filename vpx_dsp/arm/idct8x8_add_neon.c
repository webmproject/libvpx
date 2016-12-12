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

#include "./vpx_config.h"
#include "./vpx_dsp_rtcd.h"
#include "vpx_dsp/arm/idct_neon.h"
#include "vpx_dsp/arm/transpose_neon.h"
#include "vpx_dsp/txfm_common.h"

static INLINE void IDCT8x8_1D(const int16x4_t cospis0, const int16x4_t cospis1,
                              int16x8_t *a0, int16x8_t *a1, int16x8_t *a2,
                              int16x8_t *a3, int16x8_t *a4, int16x8_t *a5,
                              int16x8_t *a6, int16x8_t *a7) {
  const int16x4_t a0l = vget_low_s16(*a0);
  const int16x4_t a0h = vget_high_s16(*a0);
  const int16x4_t a1l = vget_low_s16(*a1);
  const int16x4_t a1h = vget_high_s16(*a1);
  const int16x4_t a2l = vget_low_s16(*a2);
  const int16x4_t a2h = vget_high_s16(*a2);
  const int16x4_t a3l = vget_low_s16(*a3);
  const int16x4_t a3h = vget_high_s16(*a3);
  const int16x4_t a4l = vget_low_s16(*a4);
  const int16x4_t a4h = vget_high_s16(*a4);
  const int16x4_t a5l = vget_low_s16(*a5);
  const int16x4_t a5h = vget_high_s16(*a5);
  const int16x4_t a6l = vget_low_s16(*a6);
  const int16x4_t a6h = vget_high_s16(*a6);
  const int16x4_t a7l = vget_low_s16(*a7);
  const int16x4_t a7h = vget_high_s16(*a7);
  int32x4_t b0, b1, b2, b3;
  int16x4_t c0, c1, c2, c3;
  int16x8_t d0, d1, d2, d3, d4, d5, d6, d7, e0, e1, e2, e3;

  b0 = vmull_lane_s16(a1l, cospis1, 3);
  b1 = vmull_lane_s16(a1h, cospis1, 3);
  b2 = vmull_lane_s16(a5l, cospis1, 1);
  b3 = vmull_lane_s16(a5h, cospis1, 1);
  b0 = vmlsl_lane_s16(b0, a7l, cospis1, 0);
  b1 = vmlsl_lane_s16(b1, a7h, cospis1, 0);
  b2 = vmlal_lane_s16(b2, a3l, cospis1, 2);
  b3 = vmlal_lane_s16(b3, a3h, cospis1, 2);
  c0 = vrshrn_n_s32(b0, 14);
  c1 = vrshrn_n_s32(b1, 14);
  c2 = vrshrn_n_s32(b2, 14);
  c3 = vrshrn_n_s32(b3, 14);
  d4 = vcombine_s16(c0, c1);
  d5 = vcombine_s16(c2, c3);

  b0 = vmull_lane_s16(a1l, cospis1, 0);
  b1 = vmull_lane_s16(a1h, cospis1, 0);
  b2 = vmull_lane_s16(a3l, cospis1, 1);
  b3 = vmull_lane_s16(a3h, cospis1, 1);
  b0 = vmlal_lane_s16(b0, a7l, cospis1, 3);
  b1 = vmlal_lane_s16(b1, a7h, cospis1, 3);
  b2 = vmlsl_lane_s16(b2, a5l, cospis1, 2);
  b3 = vmlsl_lane_s16(b3, a5h, cospis1, 2);
  c0 = vrshrn_n_s32(b0, 14);
  c1 = vrshrn_n_s32(b1, 14);
  c2 = vrshrn_n_s32(b2, 14);
  c3 = vrshrn_n_s32(b3, 14);
  d6 = vcombine_s16(c2, c3);
  d7 = vcombine_s16(c0, c1);

  b2 = vmull_lane_s16(a0l, cospis0, 2);
  b3 = vmull_lane_s16(a0h, cospis0, 2);
  b0 = vmlal_lane_s16(b2, a4l, cospis0, 2);
  b1 = vmlal_lane_s16(b3, a4h, cospis0, 2);
  b2 = vmlsl_lane_s16(b2, a4l, cospis0, 2);
  b3 = vmlsl_lane_s16(b3, a4h, cospis0, 2);
  c0 = vrshrn_n_s32(b0, 14);
  c1 = vrshrn_n_s32(b1, 14);
  c2 = vrshrn_n_s32(b2, 14);
  c3 = vrshrn_n_s32(b3, 14);
  e0 = vcombine_s16(c0, c1);
  e1 = vcombine_s16(c2, c3);

  b0 = vmull_lane_s16(a2l, cospis0, 3);
  b1 = vmull_lane_s16(a2h, cospis0, 3);
  b2 = vmull_lane_s16(a2l, cospis0, 1);
  b3 = vmull_lane_s16(a2h, cospis0, 1);
  b0 = vmlsl_lane_s16(b0, a6l, cospis0, 1);
  b1 = vmlsl_lane_s16(b1, a6h, cospis0, 1);
  b2 = vmlal_lane_s16(b2, a6l, cospis0, 3);
  b3 = vmlal_lane_s16(b3, a6h, cospis0, 3);
  c0 = vrshrn_n_s32(b0, 14);
  c1 = vrshrn_n_s32(b1, 14);
  c2 = vrshrn_n_s32(b2, 14);
  c3 = vrshrn_n_s32(b3, 14);
  e2 = vcombine_s16(c0, c1);
  e3 = vcombine_s16(c2, c3);

  d0 = vaddq_s16(e0, e3);
  d1 = vaddq_s16(e1, e2);
  d2 = vsubq_s16(e1, e2);
  d3 = vsubq_s16(e0, e3);

  e0 = vsubq_s16(d4, d5);
  e1 = vsubq_s16(d7, d6);
  d4 = vaddq_s16(d4, d5);
  d7 = vaddq_s16(d7, d6);
  c0 = vget_low_s16(e0);
  c1 = vget_high_s16(e0);
  c2 = vget_low_s16(e1);
  c3 = vget_high_s16(e1);

  b2 = vmull_lane_s16(c2, cospis0, 2);
  b3 = vmull_lane_s16(c3, cospis0, 2);
  b0 = vmlsl_lane_s16(b2, c0, cospis0, 2);
  b1 = vmlsl_lane_s16(b3, c1, cospis0, 2);
  b2 = vmlal_lane_s16(b2, c0, cospis0, 2);
  b3 = vmlal_lane_s16(b3, c1, cospis0, 2);
  c0 = vrshrn_n_s32(b0, 14);
  c1 = vrshrn_n_s32(b1, 14);
  c2 = vrshrn_n_s32(b2, 14);
  c3 = vrshrn_n_s32(b3, 14);
  d5 = vcombine_s16(c0, c1);
  d6 = vcombine_s16(c2, c3);

  *a0 = vaddq_s16(d0, d7);
  *a1 = vaddq_s16(d1, d6);
  *a2 = vaddq_s16(d2, d5);
  *a3 = vaddq_s16(d3, d4);
  *a4 = vsubq_s16(d3, d4);
  *a5 = vsubq_s16(d2, d5);
  *a6 = vsubq_s16(d1, d6);
  *a7 = vsubq_s16(d0, d7);
}

static INLINE void add8x8(int16x8_t a0, int16x8_t a1, int16x8_t a2,
                          int16x8_t a3, int16x8_t a4, int16x8_t a5,
                          int16x8_t a6, int16x8_t a7, uint8_t *dest,
                          const int stride) {
  const uint8_t *dst = dest;
  uint8x8_t d0, d1, d2, d3, d4, d5, d6, d7;
  uint16x8_t d0_u16, d1_u16, d2_u16, d3_u16, d4_u16, d5_u16, d6_u16, d7_u16;

  a0 = vrshrq_n_s16(a0, 5);
  a1 = vrshrq_n_s16(a1, 5);
  a2 = vrshrq_n_s16(a2, 5);
  a3 = vrshrq_n_s16(a3, 5);
  a4 = vrshrq_n_s16(a4, 5);
  a5 = vrshrq_n_s16(a5, 5);
  a6 = vrshrq_n_s16(a6, 5);
  a7 = vrshrq_n_s16(a7, 5);

  d0 = vld1_u8(dst);
  dst += stride;
  d1 = vld1_u8(dst);
  dst += stride;
  d2 = vld1_u8(dst);
  dst += stride;
  d3 = vld1_u8(dst);
  dst += stride;
  d4 = vld1_u8(dst);
  dst += stride;
  d5 = vld1_u8(dst);
  dst += stride;
  d6 = vld1_u8(dst);
  dst += stride;
  d7 = vld1_u8(dst);

  d0_u16 = vaddw_u8(vreinterpretq_u16_s16(a0), d0);
  d1_u16 = vaddw_u8(vreinterpretq_u16_s16(a1), d1);
  d2_u16 = vaddw_u8(vreinterpretq_u16_s16(a2), d2);
  d3_u16 = vaddw_u8(vreinterpretq_u16_s16(a3), d3);
  d4_u16 = vaddw_u8(vreinterpretq_u16_s16(a4), d4);
  d5_u16 = vaddw_u8(vreinterpretq_u16_s16(a5), d5);
  d6_u16 = vaddw_u8(vreinterpretq_u16_s16(a6), d6);
  d7_u16 = vaddw_u8(vreinterpretq_u16_s16(a7), d7);

  d0 = vqmovun_s16(vreinterpretq_s16_u16(d0_u16));
  d1 = vqmovun_s16(vreinterpretq_s16_u16(d1_u16));
  d2 = vqmovun_s16(vreinterpretq_s16_u16(d2_u16));
  d3 = vqmovun_s16(vreinterpretq_s16_u16(d3_u16));
  d4 = vqmovun_s16(vreinterpretq_s16_u16(d4_u16));
  d5 = vqmovun_s16(vreinterpretq_s16_u16(d5_u16));
  d6 = vqmovun_s16(vreinterpretq_s16_u16(d6_u16));
  d7 = vqmovun_s16(vreinterpretq_s16_u16(d7_u16));

  vst1_u8(dest, d0);
  dest += stride;
  vst1_u8(dest, d1);
  dest += stride;
  vst1_u8(dest, d2);
  dest += stride;
  vst1_u8(dest, d3);
  dest += stride;
  vst1_u8(dest, d4);
  dest += stride;
  vst1_u8(dest, d5);
  dest += stride;
  vst1_u8(dest, d6);
  dest += stride;
  vst1_u8(dest, d7);
}

void vpx_idct8x8_64_add_neon(const tran_low_t *input, uint8_t *dest,
                             int stride) {
  const int16x8_t cospis = vld1q_s16(kCospi);
  const int16x4_t cospis0 = vget_low_s16(cospis);   // cospi 0, 8, 16, 24
  const int16x4_t cospis1 = vget_high_s16(cospis);  // cospi 4, 12, 20, 28
  int16x8_t a0, a1, a2, a3, a4, a5, a6, a7;

  a0 = load_tran_low_to_s16q(input);
  a1 = load_tran_low_to_s16q(input + 8);
  a2 = load_tran_low_to_s16q(input + 16);
  a3 = load_tran_low_to_s16q(input + 24);
  a4 = load_tran_low_to_s16q(input + 32);
  a5 = load_tran_low_to_s16q(input + 40);
  a6 = load_tran_low_to_s16q(input + 48);
  a7 = load_tran_low_to_s16q(input + 56);

  transpose_s16_8x8(&a0, &a1, &a2, &a3, &a4, &a5, &a6, &a7);
  IDCT8x8_1D(cospis0, cospis1, &a0, &a1, &a2, &a3, &a4, &a5, &a6, &a7);
  transpose_s16_8x8(&a0, &a1, &a2, &a3, &a4, &a5, &a6, &a7);
  IDCT8x8_1D(cospis0, cospis1, &a0, &a1, &a2, &a3, &a4, &a5, &a6, &a7);
  add8x8(a0, a1, a2, a3, a4, a5, a6, a7, dest, stride);
}

static INLINE void IDCT8x4_1D(const int16x4_t cospis0, const int16x4_t cospisd0,
                              const int16x4_t cospisd1, int16x8_t *a0,
                              int16x8_t *a1, int16x8_t *a2, int16x8_t *a3,
                              int16x8_t *a4, int16x8_t *a5, int16x8_t *a6,
                              int16x8_t *a7) {
  int32x4_t b0, b1, b2, b3;
  int16x4_t c0, c1, c2, c3;
  int16x8_t d0, d1, d2, d3, d4, d5, d6, d7, e0, e1, e2, e3;

  d4 = vqrdmulhq_lane_s16(*a1, cospisd1, 3);
  d5 = vqrdmulhq_lane_s16(*a3, cospisd1, 2);
  d6 = vqrdmulhq_lane_s16(*a3, cospisd1, 1);
  d7 = vqrdmulhq_lane_s16(*a1, cospisd1, 0);
  e0 = vqrdmulhq_lane_s16(*a0, cospisd0, 2);
  e2 = vqrdmulhq_lane_s16(*a2, cospisd0, 3);
  e3 = vqrdmulhq_lane_s16(*a2, cospisd0, 1);

  d0 = vaddq_s16(e0, e3);
  d1 = vaddq_s16(e0, e2);
  d2 = vsubq_s16(e0, e2);
  d3 = vsubq_s16(e0, e3);

  e0 = vsubq_s16(d4, d5);
  e1 = vsubq_s16(d7, d6);
  d4 = vaddq_s16(d4, d5);
  d7 = vaddq_s16(d7, d6);
  c0 = vget_low_s16(e0);
  c1 = vget_high_s16(e0);
  c2 = vget_low_s16(e1);
  c3 = vget_high_s16(e1);

  b2 = vmull_lane_s16(c2, cospis0, 2);
  b3 = vmull_lane_s16(c3, cospis0, 2);
  b0 = vmlsl_lane_s16(b2, c0, cospis0, 2);
  b1 = vmlsl_lane_s16(b3, c1, cospis0, 2);
  b2 = vmlal_lane_s16(b2, c0, cospis0, 2);
  b3 = vmlal_lane_s16(b3, c1, cospis0, 2);
  c0 = vrshrn_n_s32(b0, 14);
  c1 = vrshrn_n_s32(b1, 14);
  c2 = vrshrn_n_s32(b2, 14);
  c3 = vrshrn_n_s32(b3, 14);
  d5 = vcombine_s16(c0, c1);
  d6 = vcombine_s16(c2, c3);

  *a0 = vaddq_s16(d0, d7);
  *a1 = vaddq_s16(d1, d6);
  *a2 = vaddq_s16(d2, d5);
  *a3 = vaddq_s16(d3, d4);
  *a4 = vsubq_s16(d3, d4);
  *a5 = vsubq_s16(d2, d5);
  *a6 = vsubq_s16(d1, d6);
  *a7 = vsubq_s16(d0, d7);
}

void vpx_idct8x8_12_add_neon(const tran_low_t *input, uint8_t *dest,
                             int stride) {
  const int16x8_t cospis = vld1q_s16(kCospi);
  const int16x8_t cospisd = vaddq_s16(cospis, cospis);
  const int16x4_t cospis0 = vget_low_s16(cospis);     // cospi 0, 8, 16, 24
  const int16x4_t cospisd0 = vget_low_s16(cospisd);   // doubled 0, 8, 16, 24
  const int16x4_t cospisd1 = vget_high_s16(cospisd);  // doubled 4, 12, 20, 28
  int16x8_t a0, a1, a2, a3, a4, a5, a6, a7;
  int16x4_t b0, b1, b2, b3, b4, b5, b6, b7, b8, b9, b10, b11;
  int32x4_t c0, c1;

  b8 = load_tran_low_to_s16d(input);
  b9 = load_tran_low_to_s16d(input + 8);
  b10 = load_tran_low_to_s16d(input + 16);
  b11 = load_tran_low_to_s16d(input + 24);

  transpose_s16_4x4d(&b8, &b9, &b10, &b11);

  // First transform rows
  // stage 1
  b4 = vqrdmulh_lane_s16(b9, cospisd1, 3);
  b5 = vqrdmulh_lane_s16(b11, cospisd1, 2);
  b6 = vqrdmulh_lane_s16(b11, cospisd1, 1);
  b7 = vqrdmulh_lane_s16(b9, cospisd1, 0);

  // stage 2 & stage 3 - even half
  b8 = vqrdmulh_lane_s16(b8, cospisd0, 2);
  b11 = vqrdmulh_lane_s16(b10, cospisd0, 3);
  b10 = vqrdmulh_lane_s16(b10, cospisd0, 1);

  // stage 3 -odd half
  b0 = vadd_s16(b8, b10);
  b1 = vadd_s16(b8, b11);
  b2 = vsub_s16(b8, b11);
  b3 = vsub_s16(b8, b10);

  // stage 2 - odd half
  b8 = vsub_s16(b4, b5);
  b4 = vadd_s16(b4, b5);
  b9 = vsub_s16(b7, b6);
  b7 = vadd_s16(b7, b6);

  c1 = vmull_lane_s16(b9, cospis0, 2);
  c0 = vmlsl_lane_s16(c1, b8, cospis0, 2);
  c1 = vmlal_lane_s16(c1, b8, cospis0, 2);
  b5 = vrshrn_n_s32(c0, 14);
  b6 = vrshrn_n_s32(c1, 14);

  // stage 4
  b8 = vadd_s16(b0, b7);
  b9 = vadd_s16(b1, b6);
  b10 = vadd_s16(b2, b5);
  b11 = vadd_s16(b3, b4);
  b4 = vsub_s16(b3, b4);
  b5 = vsub_s16(b2, b5);
  b6 = vsub_s16(b1, b6);
  b7 = vsub_s16(b0, b7);

  transpose_s16_4x8(b8, b9, b10, b11, b4, b5, b6, b7, &a0, &a1, &a2, &a3);
  IDCT8x4_1D(cospis0, cospisd0, cospisd1, &a0, &a1, &a2, &a3, &a4, &a5, &a6,
             &a7);
  add8x8(a0, a1, a2, a3, a4, a5, a6, a7, dest, stride);
}
