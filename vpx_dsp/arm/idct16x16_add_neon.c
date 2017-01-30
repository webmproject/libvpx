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

#include "./vpx_dsp_rtcd.h"
#include "vpx_dsp/arm/idct_neon.h"
#include "vpx_dsp/txfm_common.h"

static INLINE void wrap_low_4x2(const int32x4_t *const t32, int16x4_t *const d0,
                                int16x4_t *const d1) {
  *d0 = vrshrn_n_s32(t32[0], 14);
  *d1 = vrshrn_n_s32(t32[1], 14);
}

static INLINE void idct_cospi_2_30(const int16x8_t s0, const int16x8_t s1,
                                   const int16x4_t cospi_2_30_10_22,
                                   int16x8_t *const d0, int16x8_t *const d1) {
  int32x4_t t32[4];

  t32[0] = vmull_lane_s16(vget_low_s16(s0), cospi_2_30_10_22, 1);
  t32[1] = vmull_lane_s16(vget_high_s16(s0), cospi_2_30_10_22, 1);
  t32[2] = vmull_lane_s16(vget_low_s16(s1), cospi_2_30_10_22, 1);
  t32[3] = vmull_lane_s16(vget_high_s16(s1), cospi_2_30_10_22, 1);
  t32[0] = vmlsl_lane_s16(t32[0], vget_low_s16(s1), cospi_2_30_10_22, 0);
  t32[1] = vmlsl_lane_s16(t32[1], vget_high_s16(s1), cospi_2_30_10_22, 0);
  t32[2] = vmlal_lane_s16(t32[2], vget_low_s16(s0), cospi_2_30_10_22, 0);
  t32[3] = vmlal_lane_s16(t32[3], vget_high_s16(s0), cospi_2_30_10_22, 0);
  idct16x16_add_wrap_low_8x2(t32, d0, d1);
}

static INLINE void idct_cospi_4_28(const int16x8_t s0, const int16x8_t s1,
                                   const int16x4_t cospi_4_12_20N_28,
                                   int16x8_t *const d0, int16x8_t *const d1) {
  int32x4_t t32[4];

  t32[0] = vmull_lane_s16(vget_low_s16(s0), cospi_4_12_20N_28, 3);
  t32[1] = vmull_lane_s16(vget_high_s16(s0), cospi_4_12_20N_28, 3);
  t32[2] = vmull_lane_s16(vget_low_s16(s1), cospi_4_12_20N_28, 3);
  t32[3] = vmull_lane_s16(vget_high_s16(s1), cospi_4_12_20N_28, 3);
  t32[0] = vmlsl_lane_s16(t32[0], vget_low_s16(s1), cospi_4_12_20N_28, 0);
  t32[1] = vmlsl_lane_s16(t32[1], vget_high_s16(s1), cospi_4_12_20N_28, 0);
  t32[2] = vmlal_lane_s16(t32[2], vget_low_s16(s0), cospi_4_12_20N_28, 0);
  t32[3] = vmlal_lane_s16(t32[3], vget_high_s16(s0), cospi_4_12_20N_28, 0);
  idct16x16_add_wrap_low_8x2(t32, d0, d1);
}

static INLINE void idct_cospi_6_26(const int16x8_t s0, const int16x8_t s1,
                                   const int16x4_t cospi_6_26_14_18N,
                                   int16x8_t *const d0, int16x8_t *const d1) {
  int32x4_t t32[4];

  t32[0] = vmull_lane_s16(vget_low_s16(s0), cospi_6_26_14_18N, 0);
  t32[1] = vmull_lane_s16(vget_high_s16(s0), cospi_6_26_14_18N, 0);
  t32[2] = vmull_lane_s16(vget_low_s16(s1), cospi_6_26_14_18N, 0);
  t32[3] = vmull_lane_s16(vget_high_s16(s1), cospi_6_26_14_18N, 0);
  t32[0] = vmlal_lane_s16(t32[0], vget_low_s16(s1), cospi_6_26_14_18N, 1);
  t32[1] = vmlal_lane_s16(t32[1], vget_high_s16(s1), cospi_6_26_14_18N, 1);
  t32[2] = vmlsl_lane_s16(t32[2], vget_low_s16(s0), cospi_6_26_14_18N, 1);
  t32[3] = vmlsl_lane_s16(t32[3], vget_high_s16(s0), cospi_6_26_14_18N, 1);
  idct16x16_add_wrap_low_8x2(t32, d0, d1);
}

static INLINE void idct_cospi_8_24_d_kernel(const int16x4_t s0,
                                            const int16x4_t s1,
                                            const int16x4_t cospi_0_8_16_24,
                                            int32x4_t *const t32) {
  t32[0] = vmull_lane_s16(s0, cospi_0_8_16_24, 3);
  t32[1] = vmull_lane_s16(s1, cospi_0_8_16_24, 3);
  t32[0] = vmlsl_lane_s16(t32[0], s1, cospi_0_8_16_24, 1);
  t32[1] = vmlal_lane_s16(t32[1], s0, cospi_0_8_16_24, 1);
}

static INLINE void idct_cospi_8_24_d(const int16x4_t s0, const int16x4_t s1,
                                     const int16x4_t cospi_0_8_16_24,
                                     int16x4_t *const d0, int16x4_t *const d1) {
  int32x4_t t32[2];

  idct_cospi_8_24_d_kernel(s0, s1, cospi_0_8_16_24, t32);
  wrap_low_4x2(t32, d0, d1);
}

static INLINE void idct_cospi_8_24_neg_d(const int16x4_t s0, const int16x4_t s1,
                                         const int16x4_t cospi_0_8_16_24,
                                         int16x4_t *const d0,
                                         int16x4_t *const d1) {
  int32x4_t t32[2];

  idct_cospi_8_24_d_kernel(s0, s1, cospi_0_8_16_24, t32);
  t32[1] = vnegq_s32(t32[1]);
  wrap_low_4x2(t32, d0, d1);
}

static INLINE void idct_cospi_10_22(const int16x8_t s0, const int16x8_t s1,
                                    const int16x4_t cospi_2_30_10_22,
                                    int16x8_t *const d0, int16x8_t *const d1) {
  int32x4_t t32[4];

  t32[0] = vmull_lane_s16(vget_low_s16(s0), cospi_2_30_10_22, 3);
  t32[1] = vmull_lane_s16(vget_high_s16(s0), cospi_2_30_10_22, 3);
  t32[2] = vmull_lane_s16(vget_low_s16(s1), cospi_2_30_10_22, 3);
  t32[3] = vmull_lane_s16(vget_high_s16(s1), cospi_2_30_10_22, 3);
  t32[0] = vmlsl_lane_s16(t32[0], vget_low_s16(s1), cospi_2_30_10_22, 2);
  t32[1] = vmlsl_lane_s16(t32[1], vget_high_s16(s1), cospi_2_30_10_22, 2);
  t32[2] = vmlal_lane_s16(t32[2], vget_low_s16(s0), cospi_2_30_10_22, 2);
  t32[3] = vmlal_lane_s16(t32[3], vget_high_s16(s0), cospi_2_30_10_22, 2);
  idct16x16_add_wrap_low_8x2(t32, d0, d1);
}

static INLINE void idct_cospi_12_20(const int16x8_t s0, const int16x8_t s1,
                                    const int16x4_t cospi_4_12_20N_28,
                                    int16x8_t *const d0, int16x8_t *const d1) {
  int32x4_t t32[4];

  t32[0] = vmull_lane_s16(vget_low_s16(s0), cospi_4_12_20N_28, 1);
  t32[1] = vmull_lane_s16(vget_high_s16(s0), cospi_4_12_20N_28, 1);
  t32[2] = vmull_lane_s16(vget_low_s16(s1), cospi_4_12_20N_28, 1);
  t32[3] = vmull_lane_s16(vget_high_s16(s1), cospi_4_12_20N_28, 1);
  t32[0] = vmlal_lane_s16(t32[0], vget_low_s16(s1), cospi_4_12_20N_28, 2);
  t32[1] = vmlal_lane_s16(t32[1], vget_high_s16(s1), cospi_4_12_20N_28, 2);
  t32[2] = vmlsl_lane_s16(t32[2], vget_low_s16(s0), cospi_4_12_20N_28, 2);
  t32[3] = vmlsl_lane_s16(t32[3], vget_high_s16(s0), cospi_4_12_20N_28, 2);
  idct16x16_add_wrap_low_8x2(t32, d0, d1);
}

static INLINE void idct_cospi_14_18(const int16x8_t s0, const int16x8_t s1,
                                    const int16x4_t cospi_6_26_14_18N,
                                    int16x8_t *const d0, int16x8_t *const d1) {
  int32x4_t t32[4];

  t32[0] = vmull_lane_s16(vget_low_s16(s0), cospi_6_26_14_18N, 2);
  t32[1] = vmull_lane_s16(vget_high_s16(s0), cospi_6_26_14_18N, 2);
  t32[2] = vmull_lane_s16(vget_low_s16(s1), cospi_6_26_14_18N, 2);
  t32[3] = vmull_lane_s16(vget_high_s16(s1), cospi_6_26_14_18N, 2);
  t32[0] = vmlal_lane_s16(t32[0], vget_low_s16(s1), cospi_6_26_14_18N, 3);
  t32[1] = vmlal_lane_s16(t32[1], vget_high_s16(s1), cospi_6_26_14_18N, 3);
  t32[2] = vmlsl_lane_s16(t32[2], vget_low_s16(s0), cospi_6_26_14_18N, 3);
  t32[3] = vmlsl_lane_s16(t32[3], vget_high_s16(s0), cospi_6_26_14_18N, 3);
  idct16x16_add_wrap_low_8x2(t32, d0, d1);
}

static INLINE void idct_cospi_16_16_d(const int16x4_t s0, const int16x4_t s1,
                                      const int16x4_t cospi_0_8_16_24,
                                      int16x4_t *const d0,
                                      int16x4_t *const d1) {
  int32x4_t t32[3];

  t32[2] = vmull_lane_s16(s1, cospi_0_8_16_24, 2);
  t32[0] = vmlsl_lane_s16(t32[2], s0, cospi_0_8_16_24, 2);
  t32[1] = vmlal_lane_s16(t32[2], s0, cospi_0_8_16_24, 2);
  wrap_low_4x2(t32, d0, d1);
}

static void idct16x16_256_add_half1d(const void *const input, int16_t *output,
                                     uint8_t *dest, int stride) {
  const int16x8_t cospis0 = vld1q_s16(kCospi);
  const int16x8_t cospis1 = vld1q_s16(kCospi + 8);
  const int16x4_t cospi_0_8_16_24 = vget_low_s16(cospis0);
  const int16x4_t cospi_4_12_20N_28 = vget_high_s16(cospis0);
  const int16x4_t cospi_2_30_10_22 = vget_low_s16(cospis1);
  const int16x4_t cospi_6_26_14_18N = vget_high_s16(cospis1);
  int16x8_t in[16], step1[16], step2[16], out[16];

  // Load input (16x8)
  if (output) {
    const tran_low_t *inputT = (const tran_low_t *)input;
    in[0] = load_tran_low_to_s16q(inputT);
    inputT += 8;
    in[8] = load_tran_low_to_s16q(inputT);
    inputT += 8;
    in[1] = load_tran_low_to_s16q(inputT);
    inputT += 8;
    in[9] = load_tran_low_to_s16q(inputT);
    inputT += 8;
    in[2] = load_tran_low_to_s16q(inputT);
    inputT += 8;
    in[10] = load_tran_low_to_s16q(inputT);
    inputT += 8;
    in[3] = load_tran_low_to_s16q(inputT);
    inputT += 8;
    in[11] = load_tran_low_to_s16q(inputT);
    inputT += 8;
    in[4] = load_tran_low_to_s16q(inputT);
    inputT += 8;
    in[12] = load_tran_low_to_s16q(inputT);
    inputT += 8;
    in[5] = load_tran_low_to_s16q(inputT);
    inputT += 8;
    in[13] = load_tran_low_to_s16q(inputT);
    inputT += 8;
    in[6] = load_tran_low_to_s16q(inputT);
    inputT += 8;
    in[14] = load_tran_low_to_s16q(inputT);
    inputT += 8;
    in[7] = load_tran_low_to_s16q(inputT);
    inputT += 8;
    in[15] = load_tran_low_to_s16q(inputT);
  } else {
    const int16_t *inputT = (const int16_t *)input;
    in[0] = vld1q_s16(inputT);
    inputT += 8;
    in[8] = vld1q_s16(inputT);
    inputT += 8;
    in[1] = vld1q_s16(inputT);
    inputT += 8;
    in[9] = vld1q_s16(inputT);
    inputT += 8;
    in[2] = vld1q_s16(inputT);
    inputT += 8;
    in[10] = vld1q_s16(inputT);
    inputT += 8;
    in[3] = vld1q_s16(inputT);
    inputT += 8;
    in[11] = vld1q_s16(inputT);
    inputT += 8;
    in[4] = vld1q_s16(inputT);
    inputT += 8;
    in[12] = vld1q_s16(inputT);
    inputT += 8;
    in[5] = vld1q_s16(inputT);
    inputT += 8;
    in[13] = vld1q_s16(inputT);
    inputT += 8;
    in[6] = vld1q_s16(inputT);
    inputT += 8;
    in[14] = vld1q_s16(inputT);
    inputT += 8;
    in[7] = vld1q_s16(inputT);
    inputT += 8;
    in[15] = vld1q_s16(inputT);
  }

  // Transpose
  transpose_s16_8x8(&in[0], &in[1], &in[2], &in[3], &in[4], &in[5], &in[6],
                    &in[7]);
  transpose_s16_8x8(&in[8], &in[9], &in[10], &in[11], &in[12], &in[13], &in[14],
                    &in[15]);

  // stage 1
  step1[0] = in[0 / 2];
  step1[1] = in[16 / 2];
  step1[2] = in[8 / 2];
  step1[3] = in[24 / 2];
  step1[4] = in[4 / 2];
  step1[5] = in[20 / 2];
  step1[6] = in[12 / 2];
  step1[7] = in[28 / 2];
  step1[8] = in[2 / 2];
  step1[9] = in[18 / 2];
  step1[10] = in[10 / 2];
  step1[11] = in[26 / 2];
  step1[12] = in[6 / 2];
  step1[13] = in[22 / 2];
  step1[14] = in[14 / 2];
  step1[15] = in[30 / 2];

  // stage 2
  step2[0] = step1[0];
  step2[1] = step1[1];
  step2[2] = step1[2];
  step2[3] = step1[3];
  step2[4] = step1[4];
  step2[5] = step1[5];
  step2[6] = step1[6];
  step2[7] = step1[7];
  idct_cospi_2_30(step1[8], step1[15], cospi_2_30_10_22, &step2[8], &step2[15]);
  idct_cospi_14_18(step1[9], step1[14], cospi_6_26_14_18N, &step2[9],
                   &step2[14]);
  idct_cospi_10_22(step1[10], step1[13], cospi_2_30_10_22, &step2[10],
                   &step2[13]);
  idct_cospi_6_26(step1[11], step1[12], cospi_6_26_14_18N, &step2[11],
                  &step2[12]);

  // stage 3
  step1[0] = step2[0];
  step1[1] = step2[1];
  step1[2] = step2[2];
  step1[3] = step2[3];
  idct_cospi_4_28(step2[4], step2[7], cospi_4_12_20N_28, &step1[4], &step1[7]);
  idct_cospi_12_20(step2[5], step2[6], cospi_4_12_20N_28, &step1[5], &step1[6]);
  step1[8] = vaddq_s16(step2[8], step2[9]);
  step1[9] = vsubq_s16(step2[8], step2[9]);
  step1[10] = vsubq_s16(step2[11], step2[10]);
  step1[11] = vaddq_s16(step2[11], step2[10]);
  step1[12] = vaddq_s16(step2[12], step2[13]);
  step1[13] = vsubq_s16(step2[12], step2[13]);
  step1[14] = vsubq_s16(step2[15], step2[14]);
  step1[15] = vaddq_s16(step2[15], step2[14]);

  // stage 4
  idct_cospi_16_16_q(step1[1], step1[0], cospi_0_8_16_24, &step2[1], &step2[0]);
  idct_cospi_8_24_q(step1[2], step1[3], cospi_0_8_16_24, &step2[2], &step2[3]);
  step2[4] = vaddq_s16(step1[4], step1[5]);
  step2[5] = vsubq_s16(step1[4], step1[5]);
  step2[6] = vsubq_s16(step1[7], step1[6]);
  step2[7] = vaddq_s16(step1[7], step1[6]);
  step2[8] = step1[8];
  idct_cospi_8_24_q(step1[14], step1[9], cospi_0_8_16_24, &step2[9],
                    &step2[14]);
  idct_cospi_8_24_neg_q(step1[13], step1[10], cospi_0_8_16_24, &step2[13],
                        &step2[10]);
  step2[11] = step1[11];
  step2[12] = step1[12];
  step2[15] = step1[15];

  // stage 5
  step1[0] = vaddq_s16(step2[0], step2[3]);
  step1[1] = vaddq_s16(step2[1], step2[2]);
  step1[2] = vsubq_s16(step2[1], step2[2]);
  step1[3] = vsubq_s16(step2[0], step2[3]);
  step1[4] = step2[4];
  idct_cospi_16_16_q(step2[5], step2[6], cospi_0_8_16_24, &step1[5], &step1[6]);
  step1[7] = step2[7];
  step1[8] = vaddq_s16(step2[8], step2[11]);
  step1[9] = vaddq_s16(step2[9], step2[10]);
  step1[10] = vsubq_s16(step2[9], step2[10]);
  step1[11] = vsubq_s16(step2[8], step2[11]);
  step1[12] = vsubq_s16(step2[15], step2[12]);
  step1[13] = vsubq_s16(step2[14], step2[13]);
  step1[14] = vaddq_s16(step2[14], step2[13]);
  step1[15] = vaddq_s16(step2[15], step2[12]);

  // stage 6
  step2[0] = vaddq_s16(step1[0], step1[7]);
  step2[1] = vaddq_s16(step1[1], step1[6]);
  step2[2] = vaddq_s16(step1[2], step1[5]);
  step2[3] = vaddq_s16(step1[3], step1[4]);
  step2[4] = vsubq_s16(step1[3], step1[4]);
  step2[5] = vsubq_s16(step1[2], step1[5]);
  step2[6] = vsubq_s16(step1[1], step1[6]);
  step2[7] = vsubq_s16(step1[0], step1[7]);
  idct_cospi_16_16_q(step1[10], step1[13], cospi_0_8_16_24, &step2[10],
                     &step2[13]);
  idct_cospi_16_16_q(step1[11], step1[12], cospi_0_8_16_24, &step2[11],
                     &step2[12]);
  step2[8] = step1[8];
  step2[9] = step1[9];
  step2[14] = step1[14];
  step2[15] = step1[15];

  // stage 7
  out[0] = vaddq_s16(step2[0], step2[15]);
  out[1] = vaddq_s16(step2[1], step2[14]);
  out[2] = vaddq_s16(step2[2], step2[13]);
  out[3] = vaddq_s16(step2[3], step2[12]);
  out[4] = vaddq_s16(step2[4], step2[11]);
  out[5] = vaddq_s16(step2[5], step2[10]);
  out[6] = vaddq_s16(step2[6], step2[9]);
  out[7] = vaddq_s16(step2[7], step2[8]);
  out[8] = vsubq_s16(step2[7], step2[8]);
  out[9] = vsubq_s16(step2[6], step2[9]);
  out[10] = vsubq_s16(step2[5], step2[10]);
  out[11] = vsubq_s16(step2[4], step2[11]);
  out[12] = vsubq_s16(step2[3], step2[12]);
  out[13] = vsubq_s16(step2[2], step2[13]);
  out[14] = vsubq_s16(step2[1], step2[14]);
  out[15] = vsubq_s16(step2[0], step2[15]);

  if (output) {
    // pass 1: save the result into output
    vst1q_s16(output, out[0]);
    output += 16;
    vst1q_s16(output, out[1]);
    output += 16;
    vst1q_s16(output, out[2]);
    output += 16;
    vst1q_s16(output, out[3]);
    output += 16;
    vst1q_s16(output, out[4]);
    output += 16;
    vst1q_s16(output, out[5]);
    output += 16;
    vst1q_s16(output, out[6]);
    output += 16;
    vst1q_s16(output, out[7]);
    output += 16;
    vst1q_s16(output, out[8]);
    output += 16;
    vst1q_s16(output, out[9]);
    output += 16;
    vst1q_s16(output, out[10]);
    output += 16;
    vst1q_s16(output, out[11]);
    output += 16;
    vst1q_s16(output, out[12]);
    output += 16;
    vst1q_s16(output, out[13]);
    output += 16;
    vst1q_s16(output, out[14]);
    output += 16;
    vst1q_s16(output, out[15]);
  } else {
    // pass 2: add the result to dest.
    idct16x16_add8x1(out[0], &dest, stride);
    idct16x16_add8x1(out[1], &dest, stride);
    idct16x16_add8x1(out[2], &dest, stride);
    idct16x16_add8x1(out[3], &dest, stride);
    idct16x16_add8x1(out[4], &dest, stride);
    idct16x16_add8x1(out[5], &dest, stride);
    idct16x16_add8x1(out[6], &dest, stride);
    idct16x16_add8x1(out[7], &dest, stride);
    idct16x16_add8x1(out[8], &dest, stride);
    idct16x16_add8x1(out[9], &dest, stride);
    idct16x16_add8x1(out[10], &dest, stride);
    idct16x16_add8x1(out[11], &dest, stride);
    idct16x16_add8x1(out[12], &dest, stride);
    idct16x16_add8x1(out[13], &dest, stride);
    idct16x16_add8x1(out[14], &dest, stride);
    idct16x16_add8x1(out[15], &dest, stride);
  }
}

static void idct16x16_10_add_half1d_pass1(const tran_low_t *input,
                                          int16_t *output) {
  const int16x8_t cospis0 = vld1q_s16(kCospi);
  const int16x8_t cospis1 = vld1q_s16(kCospi + 8);
  const int16x8_t cospisd0 = vaddq_s16(cospis0, cospis0);
  const int16x8_t cospisd1 = vaddq_s16(cospis1, cospis1);
  const int16x4_t cospi_0_8_16_24 = vget_low_s16(cospis0);
  const int16x4_t cospid_0_8_16_24 = vget_low_s16(cospisd0);
  const int16x4_t cospid_4_12_20N_28 = vget_high_s16(cospisd0);
  const int16x4_t cospid_2_30_10_22 = vget_low_s16(cospisd1);
  const int16x4_t cospid_6_26_14_18N = vget_high_s16(cospisd1);
  int16x4_t in[4], step1[16], step2[16], out[16];

  // Load input (4x4)
  in[0] = load_tran_low_to_s16d(input);
  input += 16;
  in[1] = load_tran_low_to_s16d(input);
  input += 16;
  in[2] = load_tran_low_to_s16d(input);
  input += 16;
  in[3] = load_tran_low_to_s16d(input);

  // Transpose
  transpose_s16_4x4d(&in[0], &in[1], &in[2], &in[3]);

  // stage 1
  step1[0] = in[0 / 2];
  step1[4] = in[4 / 2];
  step1[8] = in[2 / 2];
  step1[12] = in[6 / 2];

  // stage 2
  step2[0] = step1[0];
  step2[4] = step1[4];
  step2[8] = vqrdmulh_lane_s16(step1[8], cospid_2_30_10_22, 1);
  step2[11] = vqrdmulh_lane_s16(step1[12], cospid_6_26_14_18N, 1);
  step2[12] = vqrdmulh_lane_s16(step1[12], cospid_6_26_14_18N, 0);
  step2[15] = vqrdmulh_lane_s16(step1[8], cospid_2_30_10_22, 0);

  // stage 3
  step1[0] = step2[0];
  step1[4] = vqrdmulh_lane_s16(step2[4], cospid_4_12_20N_28, 3);
  step1[7] = vqrdmulh_lane_s16(step2[4], cospid_4_12_20N_28, 0);
  step1[8] = step2[8];
  step1[9] = step2[8];
  step1[10] = step2[11];
  step1[11] = step2[11];
  step1[12] = step2[12];
  step1[13] = step2[12];
  step1[14] = step2[15];
  step1[15] = step2[15];

  // stage 4
  step2[0] = step2[1] = vqrdmulh_lane_s16(step1[0], cospid_0_8_16_24, 2);
  step2[4] = step1[4];
  step2[5] = step1[4];
  step2[6] = step1[7];
  step2[7] = step1[7];
  step2[8] = step1[8];
  idct_cospi_8_24_d(step1[14], step1[9], cospi_0_8_16_24, &step2[9],
                    &step2[14]);
  idct_cospi_8_24_neg_d(step1[13], step1[10], cospi_0_8_16_24, &step2[13],
                        &step2[10]);
  step2[11] = step1[11];
  step2[12] = step1[12];
  step2[15] = step1[15];

  // stage 5
  step1[0] = step2[0];
  step1[1] = step2[1];
  step1[2] = step2[1];
  step1[3] = step2[0];
  step1[4] = step2[4];
  idct_cospi_16_16_d(step2[5], step2[6], cospi_0_8_16_24, &step1[5], &step1[6]);
  step1[7] = step2[7];
  step1[8] = vadd_s16(step2[8], step2[11]);
  step1[9] = vadd_s16(step2[9], step2[10]);
  step1[10] = vsub_s16(step2[9], step2[10]);
  step1[11] = vsub_s16(step2[8], step2[11]);
  step1[12] = vsub_s16(step2[15], step2[12]);
  step1[13] = vsub_s16(step2[14], step2[13]);
  step1[14] = vadd_s16(step2[14], step2[13]);
  step1[15] = vadd_s16(step2[15], step2[12]);

  // stage 6
  step2[0] = vadd_s16(step1[0], step1[7]);
  step2[1] = vadd_s16(step1[1], step1[6]);
  step2[2] = vadd_s16(step1[2], step1[5]);
  step2[3] = vadd_s16(step1[3], step1[4]);
  step2[4] = vsub_s16(step1[3], step1[4]);
  step2[5] = vsub_s16(step1[2], step1[5]);
  step2[6] = vsub_s16(step1[1], step1[6]);
  step2[7] = vsub_s16(step1[0], step1[7]);
  idct_cospi_16_16_d(step1[10], step1[13], cospi_0_8_16_24, &step2[10],
                     &step2[13]);
  idct_cospi_16_16_d(step1[11], step1[12], cospi_0_8_16_24, &step2[11],
                     &step2[12]);
  step2[8] = step1[8];
  step2[9] = step1[9];
  step2[14] = step1[14];
  step2[15] = step1[15];

  // stage 7
  out[0] = vadd_s16(step2[0], step2[15]);
  out[1] = vadd_s16(step2[1], step2[14]);
  out[2] = vadd_s16(step2[2], step2[13]);
  out[3] = vadd_s16(step2[3], step2[12]);
  out[4] = vadd_s16(step2[4], step2[11]);
  out[5] = vadd_s16(step2[5], step2[10]);
  out[6] = vadd_s16(step2[6], step2[9]);
  out[7] = vadd_s16(step2[7], step2[8]);
  out[8] = vsub_s16(step2[7], step2[8]);
  out[9] = vsub_s16(step2[6], step2[9]);
  out[10] = vsub_s16(step2[5], step2[10]);
  out[11] = vsub_s16(step2[4], step2[11]);
  out[12] = vsub_s16(step2[3], step2[12]);
  out[13] = vsub_s16(step2[2], step2[13]);
  out[14] = vsub_s16(step2[1], step2[14]);
  out[15] = vsub_s16(step2[0], step2[15]);

  // pass 1: save the result into output
  vst1_s16(output, out[0]);
  output += 4;
  vst1_s16(output, out[1]);
  output += 4;
  vst1_s16(output, out[2]);
  output += 4;
  vst1_s16(output, out[3]);
  output += 4;
  vst1_s16(output, out[4]);
  output += 4;
  vst1_s16(output, out[5]);
  output += 4;
  vst1_s16(output, out[6]);
  output += 4;
  vst1_s16(output, out[7]);
  output += 4;
  vst1_s16(output, out[8]);
  output += 4;
  vst1_s16(output, out[9]);
  output += 4;
  vst1_s16(output, out[10]);
  output += 4;
  vst1_s16(output, out[11]);
  output += 4;
  vst1_s16(output, out[12]);
  output += 4;
  vst1_s16(output, out[13]);
  output += 4;
  vst1_s16(output, out[14]);
  output += 4;
  vst1_s16(output, out[15]);
}

static void idct16x16_10_add_half1d_pass2(const int16_t *input, int16_t *output,
                                          uint8_t *dest, int stride) {
  const int16x8_t cospis0 = vld1q_s16(kCospi);
  const int16x8_t cospis1 = vld1q_s16(kCospi + 8);
  const int16x8_t cospisd0 = vaddq_s16(cospis0, cospis0);
  const int16x8_t cospisd1 = vaddq_s16(cospis1, cospis1);
  const int16x4_t cospi_0_8_16_24 = vget_low_s16(cospis0);
  const int16x4_t cospid_0_8_16_24 = vget_low_s16(cospisd0);
  const int16x4_t cospid_4_12_20N_28 = vget_high_s16(cospisd0);
  const int16x4_t cospid_2_30_10_22 = vget_low_s16(cospisd1);
  const int16x4_t cospid_6_26_14_18N = vget_high_s16(cospisd1);
  int16x4_t ind[8];
  int16x8_t in[4], step1[16], step2[16], out[16];

  // Load input (4x8)
  ind[0] = vld1_s16(input);
  input += 4;
  ind[1] = vld1_s16(input);
  input += 4;
  ind[2] = vld1_s16(input);
  input += 4;
  ind[3] = vld1_s16(input);
  input += 4;
  ind[4] = vld1_s16(input);
  input += 4;
  ind[5] = vld1_s16(input);
  input += 4;
  ind[6] = vld1_s16(input);
  input += 4;
  ind[7] = vld1_s16(input);

  // Transpose
  transpose_s16_4x8(ind[0], ind[1], ind[2], ind[3], ind[4], ind[5], ind[6],
                    ind[7], &in[0], &in[1], &in[2], &in[3]);

  // stage 1
  step1[0] = in[0 / 2];
  step1[4] = in[4 / 2];
  step1[8] = in[2 / 2];
  step1[12] = in[6 / 2];

  // stage 2
  step2[0] = step1[0];
  step2[4] = step1[4];
  step2[8] = vqrdmulhq_lane_s16(step1[8], cospid_2_30_10_22, 1);
  step2[11] = vqrdmulhq_lane_s16(step1[12], cospid_6_26_14_18N, 1);
  step2[12] = vqrdmulhq_lane_s16(step1[12], cospid_6_26_14_18N, 0);
  step2[15] = vqrdmulhq_lane_s16(step1[8], cospid_2_30_10_22, 0);

  // stage 3
  step1[0] = step2[0];
  step1[4] = vqrdmulhq_lane_s16(step2[4], cospid_4_12_20N_28, 3);
  step1[7] = vqrdmulhq_lane_s16(step2[4], cospid_4_12_20N_28, 0);
  step1[8] = step2[8];
  step1[9] = step2[8];
  step1[10] = step2[11];
  step1[11] = step2[11];
  step1[12] = step2[12];
  step1[13] = step2[12];
  step1[14] = step2[15];
  step1[15] = step2[15];

  // stage 4
  step2[0] = step2[1] = vqrdmulhq_lane_s16(step1[0], cospid_0_8_16_24, 2);
  step2[4] = step1[4];
  step2[5] = step1[4];
  step2[6] = step1[7];
  step2[7] = step1[7];
  step2[8] = step1[8];
  idct_cospi_8_24_q(step1[14], step1[9], cospi_0_8_16_24, &step2[9],
                    &step2[14]);
  idct_cospi_8_24_neg_q(step1[13], step1[10], cospi_0_8_16_24, &step2[13],
                        &step2[10]);
  step2[11] = step1[11];
  step2[12] = step1[12];
  step2[15] = step1[15];

  // stage 5
  step1[0] = step2[0];
  step1[1] = step2[1];
  step1[2] = step2[1];
  step1[3] = step2[0];
  step1[4] = step2[4];
  idct_cospi_16_16_q(step2[5], step2[6], cospi_0_8_16_24, &step1[5], &step1[6]);
  step1[7] = step2[7];
  step1[8] = vaddq_s16(step2[8], step2[11]);
  step1[9] = vaddq_s16(step2[9], step2[10]);
  step1[10] = vsubq_s16(step2[9], step2[10]);
  step1[11] = vsubq_s16(step2[8], step2[11]);
  step1[12] = vsubq_s16(step2[15], step2[12]);
  step1[13] = vsubq_s16(step2[14], step2[13]);
  step1[14] = vaddq_s16(step2[14], step2[13]);
  step1[15] = vaddq_s16(step2[15], step2[12]);

  // stage 6
  step2[0] = vaddq_s16(step1[0], step1[7]);
  step2[1] = vaddq_s16(step1[1], step1[6]);
  step2[2] = vaddq_s16(step1[2], step1[5]);
  step2[3] = vaddq_s16(step1[3], step1[4]);
  step2[4] = vsubq_s16(step1[3], step1[4]);
  step2[5] = vsubq_s16(step1[2], step1[5]);
  step2[6] = vsubq_s16(step1[1], step1[6]);
  step2[7] = vsubq_s16(step1[0], step1[7]);
  idct_cospi_16_16_q(step1[10], step1[13], cospi_0_8_16_24, &step2[10],
                     &step2[13]);
  idct_cospi_16_16_q(step1[11], step1[12], cospi_0_8_16_24, &step2[11],
                     &step2[12]);
  step2[8] = step1[8];
  step2[9] = step1[9];
  step2[14] = step1[14];
  step2[15] = step1[15];

  // stage 7
  out[0] = vaddq_s16(step2[0], step2[15]);
  out[1] = vaddq_s16(step2[1], step2[14]);
  out[2] = vaddq_s16(step2[2], step2[13]);
  out[3] = vaddq_s16(step2[3], step2[12]);
  out[4] = vaddq_s16(step2[4], step2[11]);
  out[5] = vaddq_s16(step2[5], step2[10]);
  out[6] = vaddq_s16(step2[6], step2[9]);
  out[7] = vaddq_s16(step2[7], step2[8]);
  out[8] = vsubq_s16(step2[7], step2[8]);
  out[9] = vsubq_s16(step2[6], step2[9]);
  out[10] = vsubq_s16(step2[5], step2[10]);
  out[11] = vsubq_s16(step2[4], step2[11]);
  out[12] = vsubq_s16(step2[3], step2[12]);
  out[13] = vsubq_s16(step2[2], step2[13]);
  out[14] = vsubq_s16(step2[1], step2[14]);
  out[15] = vsubq_s16(step2[0], step2[15]);

  if (output) {
    // pass 1: save the result into output
    vst1q_s16(output, out[0]);
    output += 16;
    vst1q_s16(output, out[1]);
    output += 16;
    vst1q_s16(output, out[2]);
    output += 16;
    vst1q_s16(output, out[3]);
    output += 16;
    vst1q_s16(output, out[4]);
    output += 16;
    vst1q_s16(output, out[5]);
    output += 16;
    vst1q_s16(output, out[6]);
    output += 16;
    vst1q_s16(output, out[7]);
    output += 16;
    vst1q_s16(output, out[8]);
    output += 16;
    vst1q_s16(output, out[9]);
    output += 16;
    vst1q_s16(output, out[10]);
    output += 16;
    vst1q_s16(output, out[11]);
    output += 16;
    vst1q_s16(output, out[12]);
    output += 16;
    vst1q_s16(output, out[13]);
    output += 16;
    vst1q_s16(output, out[14]);
    output += 16;
    vst1q_s16(output, out[15]);
  } else {
    // pass 2: add the result to dest.
    idct16x16_add8x1(out[0], &dest, stride);
    idct16x16_add8x1(out[1], &dest, stride);
    idct16x16_add8x1(out[2], &dest, stride);
    idct16x16_add8x1(out[3], &dest, stride);
    idct16x16_add8x1(out[4], &dest, stride);
    idct16x16_add8x1(out[5], &dest, stride);
    idct16x16_add8x1(out[6], &dest, stride);
    idct16x16_add8x1(out[7], &dest, stride);
    idct16x16_add8x1(out[8], &dest, stride);
    idct16x16_add8x1(out[9], &dest, stride);
    idct16x16_add8x1(out[10], &dest, stride);
    idct16x16_add8x1(out[11], &dest, stride);
    idct16x16_add8x1(out[12], &dest, stride);
    idct16x16_add8x1(out[13], &dest, stride);
    idct16x16_add8x1(out[14], &dest, stride);
    idct16x16_add8x1(out[15], &dest, stride);
  }
}

void vpx_idct16x16_256_add_neon(const tran_low_t *input, uint8_t *dest,
                                int stride) {
  int16_t row_idct_output[16 * 16];

  // pass 1
  // Parallel idct on the upper 8 rows
  idct16x16_256_add_half1d(input, row_idct_output, dest, stride);

  // Parallel idct on the lower 8 rows
  idct16x16_256_add_half1d(input + 8 * 16, row_idct_output + 8, dest, stride);

  // pass 2
  // Parallel idct to get the left 8 columns
  idct16x16_256_add_half1d(row_idct_output, NULL, dest, stride);

  // Parallel idct to get the right 8 columns
  idct16x16_256_add_half1d(row_idct_output + 16 * 8, NULL, dest + 8, stride);
}

void vpx_idct16x16_10_add_neon(const tran_low_t *input, uint8_t *dest,
                               int stride) {
  int16_t row_idct_output[4 * 16];

  // pass 1
  // Parallel idct on the upper 8 rows
  idct16x16_10_add_half1d_pass1(input, row_idct_output);

  // pass 2
  // Parallel idct to get the left 8 columns
  idct16x16_10_add_half1d_pass2(row_idct_output, NULL, dest, stride);

  // Parallel idct to get the right 8 columns
  idct16x16_10_add_half1d_pass2(row_idct_output + 4 * 8, NULL, dest + 8,
                                stride);
}
