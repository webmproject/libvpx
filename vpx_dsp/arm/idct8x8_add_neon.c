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

static INLINE void idct8x8_64_1d_bd8(const int16x4_t cospis0,
                                     const int16x4_t cospis1,
                                     int16x8_t *const io0, int16x8_t *const io1,
                                     int16x8_t *const io2, int16x8_t *const io3,
                                     int16x8_t *const io4, int16x8_t *const io5,
                                     int16x8_t *const io6,
                                     int16x8_t *const io7) {
  int16x4_t input_1l, input_1h, input_3l, input_3h, input_5l, input_5h,
      input_7l, input_7h;
  int16x4_t step1l[4], step1h[4];
  int16x8_t step1[8], step2[8];
  int32x4_t t32[8];
  int16x4_t t16[8];

  transpose_s16_8x8(io0, io1, io2, io3, io4, io5, io6, io7);

  // stage 1
  input_1l = vget_low_s16(*io1);
  input_1h = vget_high_s16(*io1);
  input_3l = vget_low_s16(*io3);
  input_3h = vget_high_s16(*io3);
  input_5l = vget_low_s16(*io5);
  input_5h = vget_high_s16(*io5);
  input_7l = vget_low_s16(*io7);
  input_7h = vget_high_s16(*io7);
  step1l[0] = vget_low_s16(*io0);
  step1h[0] = vget_high_s16(*io0);
  step1l[1] = vget_low_s16(*io2);
  step1h[1] = vget_high_s16(*io2);
  step1l[2] = vget_low_s16(*io4);
  step1h[2] = vget_high_s16(*io4);
  step1l[3] = vget_low_s16(*io6);
  step1h[3] = vget_high_s16(*io6);

  t32[0] = vmull_lane_s16(input_1l, cospis1, 3);
  t32[1] = vmull_lane_s16(input_1h, cospis1, 3);
  t32[2] = vmull_lane_s16(input_3l, cospis1, 2);
  t32[3] = vmull_lane_s16(input_3h, cospis1, 2);
  t32[4] = vmull_lane_s16(input_3l, cospis1, 1);
  t32[5] = vmull_lane_s16(input_3h, cospis1, 1);
  t32[6] = vmull_lane_s16(input_1l, cospis1, 0);
  t32[7] = vmull_lane_s16(input_1h, cospis1, 0);
  t32[0] = vmlsl_lane_s16(t32[0], input_7l, cospis1, 0);
  t32[1] = vmlsl_lane_s16(t32[1], input_7h, cospis1, 0);
  t32[2] = vmlal_lane_s16(t32[2], input_5l, cospis1, 1);
  t32[3] = vmlal_lane_s16(t32[3], input_5h, cospis1, 1);
  t32[4] = vmlsl_lane_s16(t32[4], input_5l, cospis1, 2);
  t32[5] = vmlsl_lane_s16(t32[5], input_5h, cospis1, 2);
  t32[6] = vmlal_lane_s16(t32[6], input_7l, cospis1, 3);
  t32[7] = vmlal_lane_s16(t32[7], input_7h, cospis1, 3);
  t16[0] = vrshrn_n_s32(t32[0], 14);
  t16[1] = vrshrn_n_s32(t32[1], 14);
  t16[2] = vrshrn_n_s32(t32[2], 14);
  t16[3] = vrshrn_n_s32(t32[3], 14);
  t16[4] = vrshrn_n_s32(t32[4], 14);
  t16[5] = vrshrn_n_s32(t32[5], 14);
  t16[6] = vrshrn_n_s32(t32[6], 14);
  t16[7] = vrshrn_n_s32(t32[7], 14);
  step1[4] = vcombine_s16(t16[0], t16[1]);
  step1[5] = vcombine_s16(t16[2], t16[3]);
  step1[6] = vcombine_s16(t16[4], t16[5]);
  step1[7] = vcombine_s16(t16[6], t16[7]);

  // stage 2
  t32[2] = vmull_lane_s16(step1l[0], cospis0, 2);
  t32[3] = vmull_lane_s16(step1h[0], cospis0, 2);
  t32[4] = vmull_lane_s16(step1l[1], cospis0, 3);
  t32[5] = vmull_lane_s16(step1h[1], cospis0, 3);
  t32[6] = vmull_lane_s16(step1l[1], cospis0, 1);
  t32[7] = vmull_lane_s16(step1h[1], cospis0, 1);
  t32[0] = vmlal_lane_s16(t32[2], step1l[2], cospis0, 2);
  t32[1] = vmlal_lane_s16(t32[3], step1h[2], cospis0, 2);
  t32[2] = vmlsl_lane_s16(t32[2], step1l[2], cospis0, 2);
  t32[3] = vmlsl_lane_s16(t32[3], step1h[2], cospis0, 2);
  t32[4] = vmlsl_lane_s16(t32[4], step1l[3], cospis0, 1);
  t32[5] = vmlsl_lane_s16(t32[5], step1h[3], cospis0, 1);
  t32[6] = vmlal_lane_s16(t32[6], step1l[3], cospis0, 3);
  t32[7] = vmlal_lane_s16(t32[7], step1h[3], cospis0, 3);
  t16[0] = vrshrn_n_s32(t32[0], 14);
  t16[1] = vrshrn_n_s32(t32[1], 14);
  t16[2] = vrshrn_n_s32(t32[2], 14);
  t16[3] = vrshrn_n_s32(t32[3], 14);
  t16[4] = vrshrn_n_s32(t32[4], 14);
  t16[5] = vrshrn_n_s32(t32[5], 14);
  t16[6] = vrshrn_n_s32(t32[6], 14);
  t16[7] = vrshrn_n_s32(t32[7], 14);
  step2[0] = vcombine_s16(t16[0], t16[1]);
  step2[1] = vcombine_s16(t16[2], t16[3]);
  step2[2] = vcombine_s16(t16[4], t16[5]);
  step2[3] = vcombine_s16(t16[6], t16[7]);

  step2[4] = vaddq_s16(step1[4], step1[5]);
  step2[5] = vsubq_s16(step1[4], step1[5]);
  step2[6] = vsubq_s16(step1[7], step1[6]);
  step2[7] = vaddq_s16(step1[7], step1[6]);

  // stage 3
  step1[0] = vaddq_s16(step2[0], step2[3]);
  step1[1] = vaddq_s16(step2[1], step2[2]);
  step1[2] = vsubq_s16(step2[1], step2[2]);
  step1[3] = vsubq_s16(step2[0], step2[3]);

  t32[2] = vmull_lane_s16(vget_low_s16(step2[6]), cospis0, 2);
  t32[3] = vmull_lane_s16(vget_high_s16(step2[6]), cospis0, 2);
  t32[0] = vmlsl_lane_s16(t32[2], vget_low_s16(step2[5]), cospis0, 2);
  t32[1] = vmlsl_lane_s16(t32[3], vget_high_s16(step2[5]), cospis0, 2);
  t32[2] = vmlal_lane_s16(t32[2], vget_low_s16(step2[5]), cospis0, 2);
  t32[3] = vmlal_lane_s16(t32[3], vget_high_s16(step2[5]), cospis0, 2);
  t16[0] = vrshrn_n_s32(t32[0], 14);
  t16[1] = vrshrn_n_s32(t32[1], 14);
  t16[2] = vrshrn_n_s32(t32[2], 14);
  t16[3] = vrshrn_n_s32(t32[3], 14);
  step1[5] = vcombine_s16(t16[0], t16[1]);
  step1[6] = vcombine_s16(t16[2], t16[3]);

  // stage 4
  *io0 = vaddq_s16(step1[0], step2[7]);
  *io1 = vaddq_s16(step1[1], step1[6]);
  *io2 = vaddq_s16(step1[2], step1[5]);
  *io3 = vaddq_s16(step1[3], step2[4]);
  *io4 = vsubq_s16(step1[3], step2[4]);
  *io5 = vsubq_s16(step1[2], step1[5]);
  *io6 = vsubq_s16(step1[1], step1[6]);
  *io7 = vsubq_s16(step1[0], step2[7]);
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
  int16x8_t a0 = load_tran_low_to_s16q(input);
  int16x8_t a1 = load_tran_low_to_s16q(input + 8);
  int16x8_t a2 = load_tran_low_to_s16q(input + 16);
  int16x8_t a3 = load_tran_low_to_s16q(input + 24);
  int16x8_t a4 = load_tran_low_to_s16q(input + 32);
  int16x8_t a5 = load_tran_low_to_s16q(input + 40);
  int16x8_t a6 = load_tran_low_to_s16q(input + 48);
  int16x8_t a7 = load_tran_low_to_s16q(input + 56);

  idct8x8_64_1d_bd8(cospis0, cospis1, &a0, &a1, &a2, &a3, &a4, &a5, &a6, &a7);
  idct8x8_64_1d_bd8(cospis0, cospis1, &a0, &a1, &a2, &a3, &a4, &a5, &a6, &a7);
  add8x8(a0, a1, a2, a3, a4, a5, a6, a7, dest, stride);
}

static INLINE void idct8x8_12_pass1_bd8(
    const int16x4_t cospis0, const int16x4_t cospisd0, const int16x4_t cospisd1,
    int16x4_t *const io0, int16x4_t *const io1, int16x4_t *const io2,
    int16x4_t *const io3, int16x4_t *const io4, int16x4_t *const io5,
    int16x4_t *const io6, int16x4_t *const io7) {
  int16x4_t step1[8], step2[8];
  int32x4_t t32[2];

  transpose_s16_4x4d(io0, io1, io2, io3);

  // stage 1
  step1[4] = vqrdmulh_lane_s16(*io1, cospisd1, 3);
  step1[5] = vqrdmulh_lane_s16(*io3, cospisd1, 2);
  step1[6] = vqrdmulh_lane_s16(*io3, cospisd1, 1);
  step1[7] = vqrdmulh_lane_s16(*io1, cospisd1, 0);

  // stage 2
  step2[0] = vqrdmulh_lane_s16(*io0, cospisd0, 2);
  step2[2] = vqrdmulh_lane_s16(*io2, cospisd0, 3);
  step2[3] = vqrdmulh_lane_s16(*io2, cospisd0, 1);

  step2[4] = vadd_s16(step1[4], step1[5]);
  step2[5] = vsub_s16(step1[4], step1[5]);
  step2[6] = vsub_s16(step1[7], step1[6]);
  step2[7] = vadd_s16(step1[7], step1[6]);

  // stage 3
  step1[0] = vadd_s16(step2[0], step2[3]);
  step1[1] = vadd_s16(step2[0], step2[2]);
  step1[2] = vsub_s16(step2[0], step2[2]);
  step1[3] = vsub_s16(step2[0], step2[3]);

  t32[1] = vmull_lane_s16(step2[6], cospis0, 2);
  t32[0] = vmlsl_lane_s16(t32[1], step2[5], cospis0, 2);
  t32[1] = vmlal_lane_s16(t32[1], step2[5], cospis0, 2);
  step1[5] = vrshrn_n_s32(t32[0], 14);
  step1[6] = vrshrn_n_s32(t32[1], 14);

  // stage 4
  *io0 = vadd_s16(step1[0], step2[7]);
  *io1 = vadd_s16(step1[1], step1[6]);
  *io2 = vadd_s16(step1[2], step1[5]);
  *io3 = vadd_s16(step1[3], step2[4]);
  *io4 = vsub_s16(step1[3], step2[4]);
  *io5 = vsub_s16(step1[2], step1[5]);
  *io6 = vsub_s16(step1[1], step1[6]);
  *io7 = vsub_s16(step1[0], step2[7]);
}

static INLINE void idct8x8_12_pass2_bd8(
    const int16x4_t cospis0, const int16x4_t cospisd0, const int16x4_t cospisd1,
    const int16x4_t input0, const int16x4_t input1, const int16x4_t input2,
    const int16x4_t input3, const int16x4_t input4, const int16x4_t input5,
    const int16x4_t input6, const int16x4_t input7, int16x8_t *const output0,
    int16x8_t *const output1, int16x8_t *const output2,
    int16x8_t *const output3, int16x8_t *const output4,
    int16x8_t *const output5, int16x8_t *const output6,
    int16x8_t *const output7) {
  int16x8_t in[4];
  int16x8_t step1[8], step2[8];
  int32x4_t t32[8];
  int16x4_t t16[8];

  transpose_s16_4x8(input0, input1, input2, input3, input4, input5, input6,
                    input7, &in[0], &in[1], &in[2], &in[3]);

  // stage 1
  step1[4] = vqrdmulhq_lane_s16(in[1], cospisd1, 3);
  step1[5] = vqrdmulhq_lane_s16(in[3], cospisd1, 2);
  step1[6] = vqrdmulhq_lane_s16(in[3], cospisd1, 1);
  step1[7] = vqrdmulhq_lane_s16(in[1], cospisd1, 0);

  // stage 2
  step2[0] = vqrdmulhq_lane_s16(in[0], cospisd0, 2);
  step2[2] = vqrdmulhq_lane_s16(in[2], cospisd0, 3);
  step2[3] = vqrdmulhq_lane_s16(in[2], cospisd0, 1);

  step2[4] = vaddq_s16(step1[4], step1[5]);
  step2[5] = vsubq_s16(step1[4], step1[5]);
  step2[6] = vsubq_s16(step1[7], step1[6]);
  step2[7] = vaddq_s16(step1[7], step1[6]);

  // stage 3
  step1[0] = vaddq_s16(step2[0], step2[3]);
  step1[1] = vaddq_s16(step2[0], step2[2]);
  step1[2] = vsubq_s16(step2[0], step2[2]);
  step1[3] = vsubq_s16(step2[0], step2[3]);

  t32[2] = vmull_lane_s16(vget_low_s16(step2[6]), cospis0, 2);
  t32[3] = vmull_lane_s16(vget_high_s16(step2[6]), cospis0, 2);
  t32[0] = vmlsl_lane_s16(t32[2], vget_low_s16(step2[5]), cospis0, 2);
  t32[1] = vmlsl_lane_s16(t32[3], vget_high_s16(step2[5]), cospis0, 2);
  t32[2] = vmlal_lane_s16(t32[2], vget_low_s16(step2[5]), cospis0, 2);
  t32[3] = vmlal_lane_s16(t32[3], vget_high_s16(step2[5]), cospis0, 2);
  t16[0] = vrshrn_n_s32(t32[0], 14);
  t16[1] = vrshrn_n_s32(t32[1], 14);
  t16[2] = vrshrn_n_s32(t32[2], 14);
  t16[3] = vrshrn_n_s32(t32[3], 14);
  step1[5] = vcombine_s16(t16[0], t16[1]);
  step1[6] = vcombine_s16(t16[2], t16[3]);

  // stage 4
  *output0 = vaddq_s16(step1[0], step2[7]);
  *output1 = vaddq_s16(step1[1], step1[6]);
  *output2 = vaddq_s16(step1[2], step1[5]);
  *output3 = vaddq_s16(step1[3], step2[4]);
  *output4 = vsubq_s16(step1[3], step2[4]);
  *output5 = vsubq_s16(step1[2], step1[5]);
  *output6 = vsubq_s16(step1[1], step1[6]);
  *output7 = vsubq_s16(step1[0], step2[7]);
}

void vpx_idct8x8_12_add_neon(const tran_low_t *input, uint8_t *dest,
                             int stride) {
  const int16x8_t cospis = vld1q_s16(kCospi);
  const int16x8_t cospisd = vaddq_s16(cospis, cospis);
  const int16x4_t cospis0 = vget_low_s16(cospis);     // cospi 0, 8, 16, 24
  const int16x4_t cospisd0 = vget_low_s16(cospisd);   // doubled 0, 8, 16, 24
  const int16x4_t cospisd1 = vget_high_s16(cospisd);  // doubled 4, 12, 20, 28
  int16x4_t a0, a1, a2, a3, a4, a5, a6, a7;
  int16x8_t b0, b1, b2, b3, b4, b5, b6, b7;

  a0 = load_tran_low_to_s16d(input);
  a1 = load_tran_low_to_s16d(input + 8);
  a2 = load_tran_low_to_s16d(input + 16);
  a3 = load_tran_low_to_s16d(input + 24);

  idct8x8_12_pass1_bd8(cospis0, cospisd0, cospisd1, &a0, &a1, &a2, &a3, &a4,
                       &a5, &a6, &a7);
  idct8x8_12_pass2_bd8(cospis0, cospisd0, cospisd1, a0, a1, a2, a3, a4, a5, a6,
                       a7, &b0, &b1, &b2, &b3, &b4, &b5, &b6, &b7);
  add8x8(b0, b1, b2, b3, b4, b5, b6, b7, dest, stride);
}
