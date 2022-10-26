/*
 *  Copyright (c) 2017 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#ifndef VPX_VPX_DSP_ARM_FDCT16X16_NEON_H_
#define VPX_VPX_DSP_ARM_FDCT16X16_NEON_H_

#include <arm_neon.h>

#include "fdct_neon.h"

static INLINE void load(const int16_t *a, int stride, int16x8_t *b /*[16]*/) {
  b[0] = vld1q_s16(a);
  a += stride;
  b[1] = vld1q_s16(a);
  a += stride;
  b[2] = vld1q_s16(a);
  a += stride;
  b[3] = vld1q_s16(a);
  a += stride;
  b[4] = vld1q_s16(a);
  a += stride;
  b[5] = vld1q_s16(a);
  a += stride;
  b[6] = vld1q_s16(a);
  a += stride;
  b[7] = vld1q_s16(a);
  a += stride;
  b[8] = vld1q_s16(a);
  a += stride;
  b[9] = vld1q_s16(a);
  a += stride;
  b[10] = vld1q_s16(a);
  a += stride;
  b[11] = vld1q_s16(a);
  a += stride;
  b[12] = vld1q_s16(a);
  a += stride;
  b[13] = vld1q_s16(a);
  a += stride;
  b[14] = vld1q_s16(a);
  a += stride;
  b[15] = vld1q_s16(a);
}

// Store 8 16x8 values, assuming stride == 16.
static INLINE void store(tran_low_t *a, const int16x8_t *b /*[8]*/) {
  store_s16q_to_tran_low(a, b[0]);
  a += 16;
  store_s16q_to_tran_low(a, b[1]);
  a += 16;
  store_s16q_to_tran_low(a, b[2]);
  a += 16;
  store_s16q_to_tran_low(a, b[3]);
  a += 16;
  store_s16q_to_tran_low(a, b[4]);
  a += 16;
  store_s16q_to_tran_low(a, b[5]);
  a += 16;
  store_s16q_to_tran_low(a, b[6]);
  a += 16;
  store_s16q_to_tran_low(a, b[7]);
}

// Load step of each pass. Add and subtract clear across the input, requiring
// all 16 values to be loaded. For the first pass it also multiplies by 4.

// To maybe reduce register usage this could be combined with the load() step to
// get the first 4 and last 4 values, cross those, then load the middle 8 values
// and cross them.
static INLINE void scale_input(const int16x8_t *a /*[16]*/,
                               int16x8_t *b /*[16]*/) {
  b[0] = vshlq_n_s16(a[0], 2);
  b[1] = vshlq_n_s16(a[1], 2);
  b[2] = vshlq_n_s16(a[2], 2);
  b[3] = vshlq_n_s16(a[3], 2);
  b[4] = vshlq_n_s16(a[4], 2);
  b[5] = vshlq_n_s16(a[5], 2);
  b[6] = vshlq_n_s16(a[6], 2);
  b[7] = vshlq_n_s16(a[7], 2);

  b[8] = vshlq_n_s16(a[8], 2);
  b[9] = vshlq_n_s16(a[9], 2);
  b[10] = vshlq_n_s16(a[10], 2);
  b[11] = vshlq_n_s16(a[11], 2);
  b[12] = vshlq_n_s16(a[12], 2);
  b[13] = vshlq_n_s16(a[13], 2);
  b[14] = vshlq_n_s16(a[14], 2);
  b[15] = vshlq_n_s16(a[15], 2);
}

static INLINE void cross_input(const int16x8_t *a /*[16]*/,
                               int16x8_t *b /*[16]*/) {
  b[0] = vaddq_s16(a[0], a[15]);
  b[1] = vaddq_s16(a[1], a[14]);
  b[2] = vaddq_s16(a[2], a[13]);
  b[3] = vaddq_s16(a[3], a[12]);
  b[4] = vaddq_s16(a[4], a[11]);
  b[5] = vaddq_s16(a[5], a[10]);
  b[6] = vaddq_s16(a[6], a[9]);
  b[7] = vaddq_s16(a[7], a[8]);

  b[8] = vsubq_s16(a[7], a[8]);
  b[9] = vsubq_s16(a[6], a[9]);
  b[10] = vsubq_s16(a[5], a[10]);
  b[11] = vsubq_s16(a[4], a[11]);
  b[12] = vsubq_s16(a[3], a[12]);
  b[13] = vsubq_s16(a[2], a[13]);
  b[14] = vsubq_s16(a[1], a[14]);
  b[15] = vsubq_s16(a[0], a[15]);
}

static INLINE void load_cross(const int16_t *a, int stride,
                              int16x8_t *b /*[16]*/) {
  b[0] = vaddq_s16(vld1q_s16(a + 0 * stride), vld1q_s16(a + 15 * stride));
  b[1] = vaddq_s16(vld1q_s16(a + 1 * stride), vld1q_s16(a + 14 * stride));
  b[2] = vaddq_s16(vld1q_s16(a + 2 * stride), vld1q_s16(a + 13 * stride));
  b[3] = vaddq_s16(vld1q_s16(a + 3 * stride), vld1q_s16(a + 12 * stride));
  b[4] = vaddq_s16(vld1q_s16(a + 4 * stride), vld1q_s16(a + 11 * stride));
  b[5] = vaddq_s16(vld1q_s16(a + 5 * stride), vld1q_s16(a + 10 * stride));
  b[6] = vaddq_s16(vld1q_s16(a + 6 * stride), vld1q_s16(a + 9 * stride));
  b[7] = vaddq_s16(vld1q_s16(a + 7 * stride), vld1q_s16(a + 8 * stride));

  b[8] = vsubq_s16(vld1q_s16(a + 7 * stride), vld1q_s16(a + 8 * stride));
  b[9] = vsubq_s16(vld1q_s16(a + 6 * stride), vld1q_s16(a + 9 * stride));
  b[10] = vsubq_s16(vld1q_s16(a + 5 * stride), vld1q_s16(a + 10 * stride));
  b[11] = vsubq_s16(vld1q_s16(a + 4 * stride), vld1q_s16(a + 11 * stride));
  b[12] = vsubq_s16(vld1q_s16(a + 3 * stride), vld1q_s16(a + 12 * stride));
  b[13] = vsubq_s16(vld1q_s16(a + 2 * stride), vld1q_s16(a + 13 * stride));
  b[14] = vsubq_s16(vld1q_s16(a + 1 * stride), vld1q_s16(a + 14 * stride));
  b[15] = vsubq_s16(vld1q_s16(a + 0 * stride), vld1q_s16(a + 15 * stride));
}

// Quarter round at the beginning of the second pass. Can't use vrshr (rounding)
// because this only adds 1, not 1 << 2.
static INLINE void partial_round_shift(int16x8_t *a /*[16]*/) {
  const int16x8_t one = vdupq_n_s16(1);
  a[0] = vshrq_n_s16(vaddq_s16(a[0], one), 2);
  a[1] = vshrq_n_s16(vaddq_s16(a[1], one), 2);
  a[2] = vshrq_n_s16(vaddq_s16(a[2], one), 2);
  a[3] = vshrq_n_s16(vaddq_s16(a[3], one), 2);
  a[4] = vshrq_n_s16(vaddq_s16(a[4], one), 2);
  a[5] = vshrq_n_s16(vaddq_s16(a[5], one), 2);
  a[6] = vshrq_n_s16(vaddq_s16(a[6], one), 2);
  a[7] = vshrq_n_s16(vaddq_s16(a[7], one), 2);
  a[8] = vshrq_n_s16(vaddq_s16(a[8], one), 2);
  a[9] = vshrq_n_s16(vaddq_s16(a[9], one), 2);
  a[10] = vshrq_n_s16(vaddq_s16(a[10], one), 2);
  a[11] = vshrq_n_s16(vaddq_s16(a[11], one), 2);
  a[12] = vshrq_n_s16(vaddq_s16(a[12], one), 2);
  a[13] = vshrq_n_s16(vaddq_s16(a[13], one), 2);
  a[14] = vshrq_n_s16(vaddq_s16(a[14], one), 2);
  a[15] = vshrq_n_s16(vaddq_s16(a[15], one), 2);
}

// Main body of fdct16x16.
static void vpx_fdct8x16_body(const int16x8_t *in /*[16]*/,
                              int16x8_t *out /*[16]*/) {
  int16x8_t s[8];
  int16x8_t x[4];
  int16x8_t step[8];

  // stage 1
  // From fwd_txfm.c: Work on the first eight values; fdct8(input,
  // even_results);"
  s[0] = vaddq_s16(in[0], in[7]);
  s[1] = vaddq_s16(in[1], in[6]);
  s[2] = vaddq_s16(in[2], in[5]);
  s[3] = vaddq_s16(in[3], in[4]);
  s[4] = vsubq_s16(in[3], in[4]);
  s[5] = vsubq_s16(in[2], in[5]);
  s[6] = vsubq_s16(in[1], in[6]);
  s[7] = vsubq_s16(in[0], in[7]);

  // fdct4(step, step);
  x[0] = vaddq_s16(s[0], s[3]);
  x[1] = vaddq_s16(s[1], s[2]);
  x[2] = vsubq_s16(s[1], s[2]);
  x[3] = vsubq_s16(s[0], s[3]);

  // out[0] = fdct_round_shift((x0 + x1) * cospi_16_64)
  // out[8] = fdct_round_shift((x0 - x1) * cospi_16_64)
  butterfly_one_coeff_s16_s32_fast_narrow(x[0], x[1], cospi_16_64, &out[0],
                                          &out[8]);
  // out[4]  = fdct_round_shift(x3 * cospi_8_64  + x2 * cospi_24_64);
  // out[12] = fdct_round_shift(x3 * cospi_24_64 - x2 * cospi_8_64);
  butterfly_two_coeff(x[3], x[2], cospi_8_64, cospi_24_64, &out[4], &out[12]);

  //  Stage 2
  // Re-using source s5/s6
  // s5 = fdct_round_shift((s6 - s5) * cospi_16_64)
  // s6 = fdct_round_shift((s6 + s5) * cospi_16_64)
  butterfly_one_coeff_s16_fast(s[6], s[5], cospi_16_64, &s[6], &s[5]);

  //  Stage 3
  x[0] = vaddq_s16(s[4], s[5]);
  x[1] = vsubq_s16(s[4], s[5]);
  x[2] = vsubq_s16(s[7], s[6]);
  x[3] = vaddq_s16(s[7], s[6]);

  // Stage 4
  // out[2]  = fdct_round_shift(x3 * cospi_4_64  + x0 * cospi_28_64)
  // out[14] = fdct_round_shift(x3 * cospi_28_64 - x0 * cospi_4_64)
  butterfly_two_coeff(x[3], x[0], cospi_4_64, cospi_28_64, &out[2], &out[14]);
  // out[6]  = fdct_round_shift(x2 * cospi_20_64 + x1 * cospi_12_64)
  // out[10] = fdct_round_shift(x2 * cospi_12_64 - x1 * cospi_20_64)
  butterfly_two_coeff(x[2], x[1], cospi_20_64, cospi_12_64, &out[10], &out[6]);

  // step 2
  // From fwd_txfm.c: Work on the next eight values; step1 -> odd_results"
  // That file distinguished between "in_high" and "step1" but the only
  // difference is that "in_high" is the first 8 values and "step 1" is the
  // second. Here, since they are all in one array, "step1" values are += 8.

  // step2[2] = fdct_round_shift((step1[5] - step1[2]) * cospi_16_64)
  // step2[3] = fdct_round_shift((step1[4] - step1[3]) * cospi_16_64)
  // step2[4] = fdct_round_shift((step1[4] + step1[3]) * cospi_16_64)
  // step2[5] = fdct_round_shift((step1[5] + step1[2]) * cospi_16_64)
  butterfly_one_coeff_s16_fast(in[13], in[10], cospi_16_64, &s[5], &s[2]);
  butterfly_one_coeff_s16_fast(in[12], in[11], cospi_16_64, &s[4], &s[3]);

  // step 3
  s[0] = vaddq_s16(in[8], s[3]);
  s[1] = vaddq_s16(in[9], s[2]);
  x[0] = vsubq_s16(in[9], s[2]);
  x[1] = vsubq_s16(in[8], s[3]);
  x[2] = vsubq_s16(in[15], s[4]);
  x[3] = vsubq_s16(in[14], s[5]);
  s[6] = vaddq_s16(in[14], s[5]);
  s[7] = vaddq_s16(in[15], s[4]);

  // step 4
  // step2[6] = fdct_round_shift(step3[6] * cospi_8_64  + step3[1] *
  // cospi_24_64) step2[1] = fdct_round_shift(step3[6] * cospi_24_64 - step3[1]
  // * cospi_8_64)
  butterfly_two_coeff(s[6], s[1], cospi_8_64, cospi_24_64, &s[6], &s[1]);

  // step2[2] = fdct_round_shift(step3[2] * cospi_24_64 + step3[5] * cospi_8_64)
  // step2[5] = fdct_round_shift(step3[2] * cospi_8_64  - step3[5] *
  // cospi_24_64)
  butterfly_two_coeff(x[0], x[3], cospi_24_64, cospi_8_64, &s[2], &s[5]);

  // step 5
  step[0] = vaddq_s16(s[0], s[1]);
  step[1] = vsubq_s16(s[0], s[1]);
  step[2] = vaddq_s16(x[1], s[2]);
  step[3] = vsubq_s16(x[1], s[2]);
  step[4] = vsubq_s16(x[2], s[5]);
  step[5] = vaddq_s16(x[2], s[5]);
  step[6] = vsubq_s16(s[7], s[6]);
  step[7] = vaddq_s16(s[7], s[6]);

  // step 6
  // out[9] = fdct_round_shift(step1[6] * cospi_18_64 + step1[1] * cospi_14_64)
  // out[7] = fdct_round_shift(step1[6] * cospi_14_64 - step1[1] * cospi_18_64)
  butterfly_two_coeff(step[6], step[1], cospi_18_64, cospi_14_64, &out[9],
                      &out[7]);
  // out[1]  = fdct_round_shift(step1[7] * cospi_2_64  + step1[0] * cospi_30_64)
  // out[15] = fdct_round_shift(step1[7] * cospi_30_64 - step1[0] * cospi_2_64)
  butterfly_two_coeff(step[7], step[0], cospi_2_64, cospi_30_64, &out[1],
                      &out[15]);

  // out[13] = fdct_round_shift(step1[4] * cospi_26_64 + step1[3] * cospi_6_64)
  // out[3]  = fdct_round_shift(step1[4] * cospi_6_64  - step1[3] * cospi_26_64)
  butterfly_two_coeff(step[4], step[3], cospi_26_64, cospi_6_64, &out[13],
                      &out[3]);

  // out[5]  = fdct_round_shift(step1[5] * cospi_10_64 + step1[2] * cospi_22_64)
  // out[11] = fdct_round_shift(step1[5] * cospi_22_64 - step1[2] * cospi_10_64)
  butterfly_two_coeff(step[5], step[2], cospi_10_64, cospi_22_64, &out[5],
                      &out[11]);
}

#if CONFIG_VP9_HIGHBITDEPTH

static INLINE void highbd_scale_input(const int16x8_t *a /*[16]*/,
                                      int32x4_t *left /*[16]*/,
                                      int32x4_t *right /* [16] */) {
  left[0] = vshll_n_s16(vget_low_s16(a[0]), 2);
  left[1] = vshll_n_s16(vget_low_s16(a[1]), 2);
  left[2] = vshll_n_s16(vget_low_s16(a[2]), 2);
  left[3] = vshll_n_s16(vget_low_s16(a[3]), 2);
  left[4] = vshll_n_s16(vget_low_s16(a[4]), 2);
  left[5] = vshll_n_s16(vget_low_s16(a[5]), 2);
  left[6] = vshll_n_s16(vget_low_s16(a[6]), 2);
  left[7] = vshll_n_s16(vget_low_s16(a[7]), 2);
  left[8] = vshll_n_s16(vget_low_s16(a[8]), 2);
  left[9] = vshll_n_s16(vget_low_s16(a[9]), 2);
  left[10] = vshll_n_s16(vget_low_s16(a[10]), 2);
  left[11] = vshll_n_s16(vget_low_s16(a[11]), 2);
  left[12] = vshll_n_s16(vget_low_s16(a[12]), 2);
  left[13] = vshll_n_s16(vget_low_s16(a[13]), 2);
  left[14] = vshll_n_s16(vget_low_s16(a[14]), 2);
  left[15] = vshll_n_s16(vget_low_s16(a[15]), 2);

  right[0] = vshll_n_s16(vget_high_s16(a[0]), 2);
  right[1] = vshll_n_s16(vget_high_s16(a[1]), 2);
  right[2] = vshll_n_s16(vget_high_s16(a[2]), 2);
  right[3] = vshll_n_s16(vget_high_s16(a[3]), 2);
  right[4] = vshll_n_s16(vget_high_s16(a[4]), 2);
  right[5] = vshll_n_s16(vget_high_s16(a[5]), 2);
  right[6] = vshll_n_s16(vget_high_s16(a[6]), 2);
  right[7] = vshll_n_s16(vget_high_s16(a[7]), 2);
  right[8] = vshll_n_s16(vget_high_s16(a[8]), 2);
  right[9] = vshll_n_s16(vget_high_s16(a[9]), 2);
  right[10] = vshll_n_s16(vget_high_s16(a[10]), 2);
  right[11] = vshll_n_s16(vget_high_s16(a[11]), 2);
  right[12] = vshll_n_s16(vget_high_s16(a[12]), 2);
  right[13] = vshll_n_s16(vget_high_s16(a[13]), 2);
  right[14] = vshll_n_s16(vget_high_s16(a[14]), 2);
  right[15] = vshll_n_s16(vget_high_s16(a[15]), 2);
}

static INLINE void highbd_cross_input(const int32x4_t *a_left /*[16]*/,
                                      int32x4_t *a_right /*[16]*/,
                                      int32x4_t *b_left /*[16]*/,
                                      int32x4_t *b_right /*[16]*/) {
  b_left[0] = vaddq_s32(a_left[0], a_left[15]);
  b_left[1] = vaddq_s32(a_left[1], a_left[14]);
  b_left[2] = vaddq_s32(a_left[2], a_left[13]);
  b_left[3] = vaddq_s32(a_left[3], a_left[12]);
  b_left[4] = vaddq_s32(a_left[4], a_left[11]);
  b_left[5] = vaddq_s32(a_left[5], a_left[10]);
  b_left[6] = vaddq_s32(a_left[6], a_left[9]);
  b_left[7] = vaddq_s32(a_left[7], a_left[8]);

  b_right[0] = vaddq_s32(a_right[0], a_right[15]);
  b_right[1] = vaddq_s32(a_right[1], a_right[14]);
  b_right[2] = vaddq_s32(a_right[2], a_right[13]);
  b_right[3] = vaddq_s32(a_right[3], a_right[12]);
  b_right[4] = vaddq_s32(a_right[4], a_right[11]);
  b_right[5] = vaddq_s32(a_right[5], a_right[10]);
  b_right[6] = vaddq_s32(a_right[6], a_right[9]);
  b_right[7] = vaddq_s32(a_right[7], a_right[8]);

  b_left[8] = vsubq_s32(a_left[7], a_left[8]);
  b_left[9] = vsubq_s32(a_left[6], a_left[9]);
  b_left[10] = vsubq_s32(a_left[5], a_left[10]);
  b_left[11] = vsubq_s32(a_left[4], a_left[11]);
  b_left[12] = vsubq_s32(a_left[3], a_left[12]);
  b_left[13] = vsubq_s32(a_left[2], a_left[13]);
  b_left[14] = vsubq_s32(a_left[1], a_left[14]);
  b_left[15] = vsubq_s32(a_left[0], a_left[15]);

  b_right[8] = vsubq_s32(a_right[7], a_right[8]);
  b_right[9] = vsubq_s32(a_right[6], a_right[9]);
  b_right[10] = vsubq_s32(a_right[5], a_right[10]);
  b_right[11] = vsubq_s32(a_right[4], a_right[11]);
  b_right[12] = vsubq_s32(a_right[3], a_right[12]);
  b_right[13] = vsubq_s32(a_right[2], a_right[13]);
  b_right[14] = vsubq_s32(a_right[1], a_right[14]);
  b_right[15] = vsubq_s32(a_right[0], a_right[15]);
}

static INLINE void highbd_partial_round_shift(int32x4_t *left /*[16]*/,
                                              int32x4_t *right /* [16] */) {
  const int32x4_t one = vdupq_n_s32(1);
  left[0] = vshrq_n_s32(vaddq_s32(left[0], one), 2);
  left[1] = vshrq_n_s32(vaddq_s32(left[1], one), 2);
  left[2] = vshrq_n_s32(vaddq_s32(left[2], one), 2);
  left[3] = vshrq_n_s32(vaddq_s32(left[3], one), 2);
  left[4] = vshrq_n_s32(vaddq_s32(left[4], one), 2);
  left[5] = vshrq_n_s32(vaddq_s32(left[5], one), 2);
  left[6] = vshrq_n_s32(vaddq_s32(left[6], one), 2);
  left[7] = vshrq_n_s32(vaddq_s32(left[7], one), 2);
  left[8] = vshrq_n_s32(vaddq_s32(left[8], one), 2);
  left[9] = vshrq_n_s32(vaddq_s32(left[9], one), 2);
  left[10] = vshrq_n_s32(vaddq_s32(left[10], one), 2);
  left[11] = vshrq_n_s32(vaddq_s32(left[11], one), 2);
  left[12] = vshrq_n_s32(vaddq_s32(left[12], one), 2);
  left[13] = vshrq_n_s32(vaddq_s32(left[13], one), 2);
  left[14] = vshrq_n_s32(vaddq_s32(left[14], one), 2);
  left[15] = vshrq_n_s32(vaddq_s32(left[15], one), 2);

  right[0] = vshrq_n_s32(vaddq_s32(right[0], one), 2);
  right[1] = vshrq_n_s32(vaddq_s32(right[1], one), 2);
  right[2] = vshrq_n_s32(vaddq_s32(right[2], one), 2);
  right[3] = vshrq_n_s32(vaddq_s32(right[3], one), 2);
  right[4] = vshrq_n_s32(vaddq_s32(right[4], one), 2);
  right[5] = vshrq_n_s32(vaddq_s32(right[5], one), 2);
  right[6] = vshrq_n_s32(vaddq_s32(right[6], one), 2);
  right[7] = vshrq_n_s32(vaddq_s32(right[7], one), 2);
  right[8] = vshrq_n_s32(vaddq_s32(right[8], one), 2);
  right[9] = vshrq_n_s32(vaddq_s32(right[9], one), 2);
  right[10] = vshrq_n_s32(vaddq_s32(right[10], one), 2);
  right[11] = vshrq_n_s32(vaddq_s32(right[11], one), 2);
  right[12] = vshrq_n_s32(vaddq_s32(right[12], one), 2);
  right[13] = vshrq_n_s32(vaddq_s32(right[13], one), 2);
  right[14] = vshrq_n_s32(vaddq_s32(right[14], one), 2);
  right[15] = vshrq_n_s32(vaddq_s32(right[15], one), 2);
}

// Store 16 32x4 vectors, assuming stride == 16.
static INLINE void store16_s32(tran_low_t *a, const int32x4_t *b /*[32]*/) {
  vst1q_s32(a, b[0]);
  a += 16;
  vst1q_s32(a, b[1]);
  a += 16;
  vst1q_s32(a, b[2]);
  a += 16;
  vst1q_s32(a, b[3]);
  a += 16;
  vst1q_s32(a, b[4]);
  a += 16;
  vst1q_s32(a, b[5]);
  a += 16;
  vst1q_s32(a, b[6]);
  a += 16;
  vst1q_s32(a, b[7]);
  a += 16;
  vst1q_s32(a, b[8]);
  a += 16;
  vst1q_s32(a, b[9]);
  a += 16;
  vst1q_s32(a, b[10]);
  a += 16;
  vst1q_s32(a, b[11]);
  a += 16;
  vst1q_s32(a, b[12]);
  a += 16;
  vst1q_s32(a, b[13]);
  a += 16;
  vst1q_s32(a, b[14]);
  a += 16;
  vst1q_s32(a, b[15]);
}

// Main body of fdct8x16 column
static void vpx_highbd_fdct8x16_body(int32x4_t *left /*[16]*/,
                                     int32x4_t *right /* [16] */) {
  int32x4_t sl[8];
  int32x4_t sr[8];
  int32x4_t xl[4];
  int32x4_t xr[4];
  int32x4_t inl[8];
  int32x4_t inr[8];
  int32x4_t stepl[8];
  int32x4_t stepr[8];

  // stage 1
  // From fwd_txfm.c: Work on the first eight values; fdct8(input,
  // even_results);"
  sl[0] = vaddq_s32(left[0], left[7]);
  sr[0] = vaddq_s32(right[0], right[7]);
  sl[1] = vaddq_s32(left[1], left[6]);
  sr[1] = vaddq_s32(right[1], right[6]);
  sl[2] = vaddq_s32(left[2], left[5]);
  sr[2] = vaddq_s32(right[2], right[5]);
  sl[3] = vaddq_s32(left[3], left[4]);
  sr[3] = vaddq_s32(right[3], right[4]);
  sl[4] = vsubq_s32(left[3], left[4]);
  sr[4] = vsubq_s32(right[3], right[4]);
  sl[5] = vsubq_s32(left[2], left[5]);
  sr[5] = vsubq_s32(right[2], right[5]);
  sl[6] = vsubq_s32(left[1], left[6]);
  sr[6] = vsubq_s32(right[1], right[6]);
  sl[7] = vsubq_s32(left[0], left[7]);
  sr[7] = vsubq_s32(right[0], right[7]);

  // Copy values 8-15 as we're storing in-place
  inl[0] = left[8];
  inr[0] = right[8];
  inl[1] = left[9];
  inr[1] = right[9];
  inl[2] = left[10];
  inr[2] = right[10];
  inl[3] = left[11];
  inr[3] = right[11];
  inl[4] = left[12];
  inr[4] = right[12];
  inl[5] = left[13];
  inr[5] = right[13];
  inl[6] = left[14];
  inr[6] = right[14];
  inl[7] = left[15];
  inr[7] = right[15];

  // fdct4(step, step);
  xl[0] = vaddq_s32(sl[0], sl[3]);
  xr[0] = vaddq_s32(sr[0], sr[3]);
  xl[1] = vaddq_s32(sl[1], sl[2]);
  xr[1] = vaddq_s32(sr[1], sr[2]);
  xl[2] = vsubq_s32(sl[1], sl[2]);
  xr[2] = vsubq_s32(sr[1], sr[2]);
  xl[3] = vsubq_s32(sl[0], sl[3]);
  xr[3] = vsubq_s32(sr[0], sr[3]);

  // out[0] = fdct_round_shift((x0 + x1) * cospi_16_64)
  // out[8] = fdct_round_shift((x0 - x1) * cospi_16_64)
  butterfly_one_coeff_s32_fast(xl[0], xr[0], xl[1], xr[1], cospi_16_64,
                               &left[0], &right[0], &left[8], &right[8]);

  // out[4]  = fdct_round_shift(x3 * cospi_8_64  + x2 * cospi_24_64);
  // out[12] = fdct_round_shift(x3 * cospi_24_64 - x2 * cospi_8_64);
  butterfly_two_coeff_s32_s64_narrow(xl[3], xr[3], xl[2], xr[2], cospi_8_64,
                                     cospi_24_64, &left[4], &right[4],
                                     &left[12], &right[12]);

  //  Stage 2
  // Re-using source s5/s6
  // s5 = fdct_round_shift((s6 - s5) * cospi_16_64)
  // s6 = fdct_round_shift((s6 + s5) * cospi_16_64)
  butterfly_one_coeff_s32_fast(sl[6], sr[6], sl[5], sr[5], cospi_16_64, &sl[6],
                               &sr[6], &sl[5], &sr[5]);

  //  Stage 3
  xl[0] = vaddq_s32(sl[4], sl[5]);
  xr[0] = vaddq_s32(sr[4], sr[5]);
  xl[1] = vsubq_s32(sl[4], sl[5]);
  xr[1] = vsubq_s32(sr[4], sr[5]);
  xl[2] = vsubq_s32(sl[7], sl[6]);
  xr[2] = vsubq_s32(sr[7], sr[6]);
  xl[3] = vaddq_s32(sl[7], sl[6]);
  xr[3] = vaddq_s32(sr[7], sr[6]);

  // Stage 4
  // out[2]  = fdct_round_shift(x3 * cospi_4_64  + x0 * cospi_28_64)
  // out[14] = fdct_round_shift(x3 * cospi_28_64 - x0 * cospi_4_64)
  butterfly_two_coeff_s32_s64_narrow(xl[3], xr[3], xl[0], xr[0], cospi_4_64,
                                     cospi_28_64, &left[2], &right[2],
                                     &left[14], &right[14]);
  // out[6]  = fdct_round_shift(x2 * cospi_20_64 + x1 * cospi_12_64)
  // out[10] = fdct_round_shift(x2 * cospi_12_64 - x1 * cospi_20_64)
  butterfly_two_coeff_s32_s64_narrow(xl[2], xr[2], xl[1], xr[1], cospi_20_64,
                                     cospi_12_64, &left[10], &right[10],
                                     &left[6], &right[6]);

  // step 2
  // From fwd_txfm.c: Work on the next eight values; step1 -> odd_results"
  // That file distinguished between "in_high" and "step1" but the only
  // difference is that "in_high" is the first 8 values and "step 1" is the
  // second. Here, since they are all in one array, "step1" values are += 8.

  // step2[2] = fdct_round_shift((step1[5] - step1[2]) * cospi_16_64)
  // step2[3] = fdct_round_shift((step1[4] - step1[3]) * cospi_16_64)
  // step2[4] = fdct_round_shift((step1[4] + step1[3]) * cospi_16_64)
  // step2[5] = fdct_round_shift((step1[5] + step1[2]) * cospi_16_64)
  butterfly_one_coeff_s32_fast(inl[5], inr[5], inl[2], inr[2], cospi_16_64,
                               &sl[5], &sr[5], &sl[2], &sr[2]);
  butterfly_one_coeff_s32_fast(inl[4], inr[4], inl[3], inr[3], cospi_16_64,
                               &sl[4], &sr[4], &sl[3], &sr[3]);

  // step 3
  sl[0] = vaddq_s32(inl[0], sl[3]);
  sr[0] = vaddq_s32(inr[0], sr[3]);
  sl[1] = vaddq_s32(inl[1], sl[2]);
  sr[1] = vaddq_s32(inr[1], sr[2]);
  xl[0] = vsubq_s32(inl[1], sl[2]);
  xr[0] = vsubq_s32(inr[1], sr[2]);
  xl[1] = vsubq_s32(inl[0], sl[3]);
  xr[1] = vsubq_s32(inr[0], sr[3]);
  xl[2] = vsubq_s32(inl[7], sl[4]);
  xr[2] = vsubq_s32(inr[7], sr[4]);
  xl[3] = vsubq_s32(inl[6], sl[5]);
  xr[3] = vsubq_s32(inr[6], sr[5]);
  sl[6] = vaddq_s32(inl[6], sl[5]);
  sr[6] = vaddq_s32(inr[6], sr[5]);
  sl[7] = vaddq_s32(inl[7], sl[4]);
  sr[7] = vaddq_s32(inr[7], sr[4]);

  // step 4
  // step2[6] = fdct_round_shift(step3[6] * cospi_8_64  + step3[1] *
  // cospi_24_64) step2[1] = fdct_round_shift(step3[6] * cospi_24_64 - step3[1]
  // * cospi_8_64)
  butterfly_two_coeff_s32_s64_narrow(sl[6], sr[6], sl[1], sr[1], cospi_8_64,
                                     cospi_24_64, &sl[6], &sr[6], &sl[1],
                                     &sr[1]);
  // step2[2] = fdct_round_shift(step3[2] * cospi_24_64 + step3[5] * cospi_8_64)
  // step2[5] = fdct_round_shift(step3[2] * cospi_8_64  - step3[5] *
  // cospi_24_64)
  butterfly_two_coeff_s32_s64_narrow(xl[0], xr[0], xl[3], xr[3], cospi_24_64,
                                     cospi_8_64, &sl[2], &sr[2], &sl[5],
                                     &sr[5]);

  // step 5
  stepl[0] = vaddq_s32(sl[0], sl[1]);
  stepr[0] = vaddq_s32(sr[0], sr[1]);
  stepl[1] = vsubq_s32(sl[0], sl[1]);
  stepr[1] = vsubq_s32(sr[0], sr[1]);
  stepl[2] = vaddq_s32(xl[1], sl[2]);
  stepr[2] = vaddq_s32(xr[1], sr[2]);
  stepl[3] = vsubq_s32(xl[1], sl[2]);
  stepr[3] = vsubq_s32(xr[1], sr[2]);
  stepl[4] = vsubq_s32(xl[2], sl[5]);
  stepr[4] = vsubq_s32(xr[2], sr[5]);
  stepl[5] = vaddq_s32(xl[2], sl[5]);
  stepr[5] = vaddq_s32(xr[2], sr[5]);
  stepl[6] = vsubq_s32(sl[7], sl[6]);
  stepr[6] = vsubq_s32(sr[7], sr[6]);
  stepl[7] = vaddq_s32(sl[7], sl[6]);
  stepr[7] = vaddq_s32(sr[7], sr[6]);

  // step 6
  // out[9] = fdct_round_shift(step1[6] * cospi_18_64 + step1[1] * cospi_14_64)
  // out[7] = fdct_round_shift(step1[6] * cospi_14_64 - step1[1] * cospi_18_64)
  butterfly_two_coeff_s32_s64_narrow(stepl[6], stepr[6], stepl[1], stepr[1],
                                     cospi_18_64, cospi_14_64, &left[9],
                                     &right[9], &left[7], &right[7]);
  // out[1]  = fdct_round_shift(step1[7] * cospi_2_64  + step1[0] * cospi_30_64)
  // out[15] = fdct_round_shift(step1[7] * cospi_30_64 - step1[0] * cospi_2_64)
  butterfly_two_coeff_s32_s64_narrow(stepl[7], stepr[7], stepl[0], stepr[0],
                                     cospi_2_64, cospi_30_64, &left[1],
                                     &right[1], &left[15], &right[15]);
  // out[13] = fdct_round_shift(step1[4] * cospi_26_64 + step1[3] * cospi_6_64)
  // out[3]  = fdct_round_shift(step1[4] * cospi_6_64  - step1[3] * cospi_26_64)
  butterfly_two_coeff_s32_s64_narrow(stepl[4], stepr[4], stepl[3], stepr[3],
                                     cospi_26_64, cospi_6_64, &left[13],
                                     &right[13], &left[3], &right[3]);
  // out[5]  = fdct_round_shift(step1[5] * cospi_10_64 + step1[2] * cospi_22_64)
  // out[11] = fdct_round_shift(step1[5] * cospi_22_64 - step1[2] * cospi_10_64)
  butterfly_two_coeff_s32_s64_narrow(stepl[5], stepr[5], stepl[2], stepr[2],
                                     cospi_10_64, cospi_22_64, &left[5],
                                     &right[5], &left[11], &right[11]);
}

#endif  // CONFIG_VP9_HIGHBITDEPTH

#endif  // VPX_VPX_DSP_ARM_FDCT16X16_NEON_H_
