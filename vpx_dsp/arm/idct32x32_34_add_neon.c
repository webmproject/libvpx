/*
 *  Copyright (c) 2016 The WebM project authors. All Rights Reserved.
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
#include "vpx_dsp/arm/transpose_neon.h"
#include "vpx_dsp/txfm_common.h"

// Multiply a by a_const. Saturate, shift and narrow by 14.
static int16x8_t multiply_shift_and_narrow(const int16x8_t a,
                                           const int16_t a_const) {
  // Shift by 14 + rounding will be within 16 bits for well formed streams.
  // See WRAPLOW and dct_const_round_shift for details.
  // This instruction doubles the result and returns the high half, essentially
  // resulting in a right shift by 15. By multiplying the constant first that
  // becomes a right shift by 14.
  // The largest possible value used here is
  // vpx_dsp/txfm_common.h:cospi_1_64 = 16364 (* 2 = 32728) a which falls *just*
  // within the range of int16_t (+32767 / -32768) even when negated.
  return vqrdmulhq_n_s16(a, a_const * 2);
}

// Add a and b, then multiply by ab_const. Shift and narrow by 14.
static int16x8_t add_multiply_shift_and_narrow(const int16x8_t a,
                                               const int16x8_t b,
                                               const int16_t ab_const) {
  // In both add_ and its pair, sub_, the input for well-formed streams will be
  // well within 16 bits (input to the idct is the difference between two frames
  // and will be within -255 to 255, or 9 bits)
  // However, for inputs over about 25,000 (valid for int16_t, but not for idct
  // input) this function can not use vaddq_s16.
  // In order to match existing behavior and intentionally out of range tests,
  // expand the addition up to 32 bits to prevent truncation.
  int32x4_t temp_low = vaddl_s16(vget_low_s16(a), vget_low_s16(b));
  int32x4_t temp_high = vaddl_s16(vget_high_s16(a), vget_high_s16(b));
  temp_low = vmulq_n_s32(temp_low, ab_const);
  temp_high = vmulq_n_s32(temp_high, ab_const);
  return vcombine_s16(vrshrn_n_s32(temp_low, 14), vrshrn_n_s32(temp_high, 14));
}

// Subtract b from a, then multiply by ab_const. Shift and narrow by 14.
static int16x8_t sub_multiply_shift_and_narrow(const int16x8_t a,
                                               const int16x8_t b,
                                               const int16_t ab_const) {
  int32x4_t temp_low = vsubl_s16(vget_low_s16(a), vget_low_s16(b));
  int32x4_t temp_high = vsubl_s16(vget_high_s16(a), vget_high_s16(b));
  temp_low = vmulq_n_s32(temp_low, ab_const);
  temp_high = vmulq_n_s32(temp_high, ab_const);
  return vcombine_s16(vrshrn_n_s32(temp_low, 14), vrshrn_n_s32(temp_high, 14));
}

// Multiply a by a_const and b by b_const, then accumulate. Shift and narrow by
// 14.
static int16x8_t multiply_accumulate_shift_and_narrow(const int16x8_t a,
                                                      const int16_t a_const,
                                                      const int16x8_t b,
                                                      const int16_t b_const) {
  int32x4_t temp_low = vmull_n_s16(vget_low_s16(a), a_const);
  int32x4_t temp_high = vmull_n_s16(vget_high_s16(a), a_const);
  temp_low = vmlal_n_s16(temp_low, vget_low_s16(b), b_const);
  temp_high = vmlal_n_s16(temp_high, vget_high_s16(b), b_const);
  return vcombine_s16(vrshrn_n_s32(temp_low, 14), vrshrn_n_s32(temp_high, 14));
}

// Shift the output down by 6 and add it to the destination buffer.
static void add_and_store(const int16x8_t a0, const int16x8_t a1,
                          const int16x8_t a2, const int16x8_t a3,
                          const int16x8_t a4, const int16x8_t a5,
                          const int16x8_t a6, const int16x8_t a7, uint8_t *b,
                          const int b_stride) {
  uint8x8_t b0, b1, b2, b3, b4, b5, b6, b7;
  int16x8_t c0, c1, c2, c3, c4, c5, c6, c7;
  b0 = vld1_u8(b);
  b += b_stride;
  b1 = vld1_u8(b);
  b += b_stride;
  b2 = vld1_u8(b);
  b += b_stride;
  b3 = vld1_u8(b);
  b += b_stride;
  b4 = vld1_u8(b);
  b += b_stride;
  b5 = vld1_u8(b);
  b += b_stride;
  b6 = vld1_u8(b);
  b += b_stride;
  b7 = vld1_u8(b);
  b -= (7 * b_stride);

  // c = b + (a >> 6)
  c0 = vrsraq_n_s16(vreinterpretq_s16_u16(vmovl_u8(b0)), a0, 6);
  c1 = vrsraq_n_s16(vreinterpretq_s16_u16(vmovl_u8(b1)), a1, 6);
  c2 = vrsraq_n_s16(vreinterpretq_s16_u16(vmovl_u8(b2)), a2, 6);
  c3 = vrsraq_n_s16(vreinterpretq_s16_u16(vmovl_u8(b3)), a3, 6);
  c4 = vrsraq_n_s16(vreinterpretq_s16_u16(vmovl_u8(b4)), a4, 6);
  c5 = vrsraq_n_s16(vreinterpretq_s16_u16(vmovl_u8(b5)), a5, 6);
  c6 = vrsraq_n_s16(vreinterpretq_s16_u16(vmovl_u8(b6)), a6, 6);
  c7 = vrsraq_n_s16(vreinterpretq_s16_u16(vmovl_u8(b7)), a7, 6);

  b0 = vqmovun_s16(c0);
  b1 = vqmovun_s16(c1);
  b2 = vqmovun_s16(c2);
  b3 = vqmovun_s16(c3);
  b4 = vqmovun_s16(c4);
  b5 = vqmovun_s16(c5);
  b6 = vqmovun_s16(c6);
  b7 = vqmovun_s16(c7);

  vst1_u8(b, b0);
  b += b_stride;
  vst1_u8(b, b1);
  b += b_stride;
  vst1_u8(b, b2);
  b += b_stride;
  vst1_u8(b, b3);
  b += b_stride;
  vst1_u8(b, b4);
  b += b_stride;
  vst1_u8(b, b5);
  b += b_stride;
  vst1_u8(b, b6);
  b += b_stride;
  vst1_u8(b, b7);
}

// Only for the first pass of the  _34_ variant. Since it only uses values from
// the top left 8x8 it can safely assume all the remaining values are 0 and skip
// an awful lot of calculations. In fact, only the first 6 columns make the cut.
// None of the elements in the 7th or 8th column are used so it skips any calls
// to input[67] too.
// In C this does a single row of 32 for each call. Here it transposes the top
// left 8x8 to allow using SIMD.

// vp9/common/vp9_scan.c:vp9_default_iscan_32x32 arranges the first 34 non-zero
// coefficients as follows:
//    0  1  2  3  4  5  6  7
// 0  0  2  5 10 17 25
// 1  1  4  8 15 22 30
// 2  3  7 12 18 28
// 3  6 11 16 23 31
// 4  9 14 19 29
// 5 13 20 26
// 6 21 27 33
// 7 24 32
static void idct32_6_neon(const int16_t *input, int16_t *output) {
  int16x8_t in0, in1, in2, in3, in4, in5, in6, in7;
  int16x8_t s1_0, s1_1, s1_2, s1_3, s1_4, s1_5, s1_6, s1_7, s1_8, s1_9, s1_10,
      s1_11, s1_12, s1_13, s1_14, s1_15, s1_16, s1_17, s1_18, s1_19, s1_20,
      s1_21, s1_22, s1_23, s1_24, s1_25, s1_26, s1_27, s1_28, s1_29, s1_30,
      s1_31;
  int16x8_t s2_0, s2_1, s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10,
      s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20,
      s2_21, s2_22, s2_23, s2_24, s2_25, s2_26, s2_27, s2_28, s2_29, s2_30,
      s2_31;
  int16x8_t s3_24, s3_25, s3_26, s3_27;

  in0 = vld1q_s16(input);
  input += 32;
  in1 = vld1q_s16(input);
  input += 32;
  in2 = vld1q_s16(input);
  input += 32;
  in3 = vld1q_s16(input);
  input += 32;
  in4 = vld1q_s16(input);
  input += 32;
  in5 = vld1q_s16(input);
  input += 32;
  in6 = vld1q_s16(input);
  input += 32;
  in7 = vld1q_s16(input);

  transpose_s16_8x8(&in0, &in1, &in2, &in3, &in4, &in5, &in6, &in7);

  // stage 1
  // input[1] * cospi_31_64 - input[31] * cospi_1_64 (but input[31] == 0)
  s1_16 = multiply_shift_and_narrow(in1, cospi_31_64);
  // input[1] * cospi_1_64 + input[31] * cospi_31_64 (but input[31] == 0)
  s1_31 = multiply_shift_and_narrow(in1, cospi_1_64);

  s1_20 = multiply_shift_and_narrow(in5, cospi_27_64);
  s1_27 = multiply_shift_and_narrow(in5, cospi_5_64);

  s1_23 = multiply_shift_and_narrow(in3, -cospi_29_64);
  s1_24 = multiply_shift_and_narrow(in3, cospi_3_64);

  // stage 2
  s2_8 = multiply_shift_and_narrow(in2, cospi_30_64);
  s2_15 = multiply_shift_and_narrow(in2, cospi_2_64);

  // stage 3
  s1_4 = multiply_shift_and_narrow(in4, cospi_28_64);
  s1_7 = multiply_shift_and_narrow(in4, cospi_4_64);

  s1_17 = multiply_accumulate_shift_and_narrow(s1_16, -cospi_4_64, s1_31,
                                               cospi_28_64);
  s1_30 = multiply_accumulate_shift_and_narrow(s1_16, cospi_28_64, s1_31,
                                               cospi_4_64);

  s1_21 = multiply_accumulate_shift_and_narrow(s1_20, -cospi_20_64, s1_27,
                                               cospi_12_64);
  s1_26 = multiply_accumulate_shift_and_narrow(s1_20, cospi_12_64, s1_27,
                                               cospi_20_64);

  s1_22 = multiply_accumulate_shift_and_narrow(s1_23, -cospi_12_64, s1_24,
                                               -cospi_20_64);
  s1_25 = multiply_accumulate_shift_and_narrow(s1_23, -cospi_20_64, s1_24,
                                               cospi_12_64);

  // stage 4
  s1_0 = multiply_shift_and_narrow(in0, cospi_16_64);

  s2_9 = multiply_accumulate_shift_and_narrow(s2_8, -cospi_8_64, s2_15,
                                              cospi_24_64);
  s2_14 = multiply_accumulate_shift_and_narrow(s2_8, cospi_24_64, s2_15,
                                               cospi_8_64);

  s2_20 = vsubq_s16(s1_23, s1_20);
  s2_21 = vsubq_s16(s1_22, s1_21);
  s2_22 = vaddq_s16(s1_21, s1_22);
  s2_23 = vaddq_s16(s1_20, s1_23);
  s2_24 = vaddq_s16(s1_24, s1_27);
  s2_25 = vaddq_s16(s1_25, s1_26);
  s2_26 = vsubq_s16(s1_25, s1_26);
  s2_27 = vsubq_s16(s1_24, s1_27);

  // stage 5
  s1_5 = sub_multiply_shift_and_narrow(s1_7, s1_4, cospi_16_64);
  s1_6 = add_multiply_shift_and_narrow(s1_4, s1_7, cospi_16_64);

  s1_18 = multiply_accumulate_shift_and_narrow(s1_17, -cospi_8_64, s1_30,
                                               cospi_24_64);
  s1_29 = multiply_accumulate_shift_and_narrow(s1_17, cospi_24_64, s1_30,
                                               cospi_8_64);

  s1_19 = multiply_accumulate_shift_and_narrow(s1_16, -cospi_8_64, s1_31,
                                               cospi_24_64);
  s1_28 = multiply_accumulate_shift_and_narrow(s1_16, cospi_24_64, s1_31,
                                               cospi_8_64);

  s1_20 = multiply_accumulate_shift_and_narrow(s2_20, -cospi_24_64, s2_27,
                                               -cospi_8_64);
  s1_27 = multiply_accumulate_shift_and_narrow(s2_20, -cospi_8_64, s2_27,
                                               cospi_24_64);

  s1_21 = multiply_accumulate_shift_and_narrow(s2_21, -cospi_24_64, s2_26,
                                               -cospi_8_64);
  s1_26 = multiply_accumulate_shift_and_narrow(s2_21, -cospi_8_64, s2_26,
                                               cospi_24_64);

  // stage 6
  s2_0 = vaddq_s16(s1_0, s1_7);
  s2_1 = vaddq_s16(s1_0, s1_6);
  s2_2 = vaddq_s16(s1_0, s1_5);
  s2_3 = vaddq_s16(s1_0, s1_4);
  s2_4 = vsubq_s16(s1_0, s1_4);
  s2_5 = vsubq_s16(s1_0, s1_5);
  s2_6 = vsubq_s16(s1_0, s1_6);
  s2_7 = vsubq_s16(s1_0, s1_7);

  s2_10 = sub_multiply_shift_and_narrow(s2_14, s2_9, cospi_16_64);
  s2_13 = add_multiply_shift_and_narrow(s2_9, s2_14, cospi_16_64);

  s2_11 = sub_multiply_shift_and_narrow(s2_15, s2_8, cospi_16_64);
  s2_12 = add_multiply_shift_and_narrow(s2_8, s2_15, cospi_16_64);

  s2_16 = vaddq_s16(s1_16, s2_23);
  s2_17 = vaddq_s16(s1_17, s2_22);
  s2_18 = vaddq_s16(s1_18, s1_21);
  s2_19 = vaddq_s16(s1_19, s1_20);
  s2_20 = vsubq_s16(s1_19, s1_20);
  s2_21 = vsubq_s16(s1_18, s1_21);
  s2_22 = vsubq_s16(s1_17, s2_22);
  s2_23 = vsubq_s16(s1_16, s2_23);

  s3_24 = vsubq_s16(s1_31, s2_24);
  s3_25 = vsubq_s16(s1_30, s2_25);
  s3_26 = vsubq_s16(s1_29, s1_26);
  s3_27 = vsubq_s16(s1_28, s1_27);
  s2_28 = vaddq_s16(s1_27, s1_28);
  s2_29 = vaddq_s16(s1_26, s1_29);
  s2_30 = vaddq_s16(s2_25, s1_30);
  s2_31 = vaddq_s16(s2_24, s1_31);

  // stage 7
  s1_0 = vaddq_s16(s2_0, s2_15);
  s1_1 = vaddq_s16(s2_1, s2_14);
  s1_2 = vaddq_s16(s2_2, s2_13);
  s1_3 = vaddq_s16(s2_3, s2_12);
  s1_4 = vaddq_s16(s2_4, s2_11);
  s1_5 = vaddq_s16(s2_5, s2_10);
  s1_6 = vaddq_s16(s2_6, s2_9);
  s1_7 = vaddq_s16(s2_7, s2_8);
  s1_8 = vsubq_s16(s2_7, s2_8);
  s1_9 = vsubq_s16(s2_6, s2_9);
  s1_10 = vsubq_s16(s2_5, s2_10);
  s1_11 = vsubq_s16(s2_4, s2_11);
  s1_12 = vsubq_s16(s2_3, s2_12);
  s1_13 = vsubq_s16(s2_2, s2_13);
  s1_14 = vsubq_s16(s2_1, s2_14);
  s1_15 = vsubq_s16(s2_0, s2_15);

  s1_20 = sub_multiply_shift_and_narrow(s3_27, s2_20, cospi_16_64);
  s1_27 = add_multiply_shift_and_narrow(s2_20, s3_27, cospi_16_64);

  s1_21 = sub_multiply_shift_and_narrow(s3_26, s2_21, cospi_16_64);
  s1_26 = add_multiply_shift_and_narrow(s2_21, s3_26, cospi_16_64);

  s1_22 = sub_multiply_shift_and_narrow(s3_25, s2_22, cospi_16_64);
  s1_25 = add_multiply_shift_and_narrow(s2_22, s3_25, cospi_16_64);

  s1_23 = sub_multiply_shift_and_narrow(s3_24, s2_23, cospi_16_64);
  s1_24 = add_multiply_shift_and_narrow(s2_23, s3_24, cospi_16_64);

  // final stage
  vst1q_s16(output, vaddq_s16(s1_0, s2_31));
  output += 8;
  vst1q_s16(output, vaddq_s16(s1_1, s2_30));
  output += 8;
  vst1q_s16(output, vaddq_s16(s1_2, s2_29));
  output += 8;
  vst1q_s16(output, vaddq_s16(s1_3, s2_28));
  output += 8;
  vst1q_s16(output, vaddq_s16(s1_4, s1_27));
  output += 8;
  vst1q_s16(output, vaddq_s16(s1_5, s1_26));
  output += 8;
  vst1q_s16(output, vaddq_s16(s1_6, s1_25));
  output += 8;
  vst1q_s16(output, vaddq_s16(s1_7, s1_24));
  output += 8;

  vst1q_s16(output, vaddq_s16(s1_8, s1_23));
  output += 8;
  vst1q_s16(output, vaddq_s16(s1_9, s1_22));
  output += 8;
  vst1q_s16(output, vaddq_s16(s1_10, s1_21));
  output += 8;
  vst1q_s16(output, vaddq_s16(s1_11, s1_20));
  output += 8;
  vst1q_s16(output, vaddq_s16(s1_12, s2_19));
  output += 8;
  vst1q_s16(output, vaddq_s16(s1_13, s2_18));
  output += 8;
  vst1q_s16(output, vaddq_s16(s1_14, s2_17));
  output += 8;
  vst1q_s16(output, vaddq_s16(s1_15, s2_16));
  output += 8;

  vst1q_s16(output, vsubq_s16(s1_15, s2_16));
  output += 8;
  vst1q_s16(output, vsubq_s16(s1_14, s2_17));
  output += 8;
  vst1q_s16(output, vsubq_s16(s1_13, s2_18));
  output += 8;
  vst1q_s16(output, vsubq_s16(s1_12, s2_19));
  output += 8;
  vst1q_s16(output, vsubq_s16(s1_11, s1_20));
  output += 8;
  vst1q_s16(output, vsubq_s16(s1_10, s1_21));
  output += 8;
  vst1q_s16(output, vsubq_s16(s1_9, s1_22));
  output += 8;
  vst1q_s16(output, vsubq_s16(s1_8, s1_23));
  output += 8;

  vst1q_s16(output, vsubq_s16(s1_7, s1_24));
  output += 8;
  vst1q_s16(output, vsubq_s16(s1_6, s1_25));
  output += 8;
  vst1q_s16(output, vsubq_s16(s1_5, s1_26));
  output += 8;
  vst1q_s16(output, vsubq_s16(s1_4, s1_27));
  output += 8;
  vst1q_s16(output, vsubq_s16(s1_3, s2_28));
  output += 8;
  vst1q_s16(output, vsubq_s16(s1_2, s2_29));
  output += 8;
  vst1q_s16(output, vsubq_s16(s1_1, s2_30));
  output += 8;
  vst1q_s16(output, vsubq_s16(s1_0, s2_31));
}

static void idct32_8_neon(const int16_t *input, uint8_t *output, int stride) {
  int16x8_t in0, in1, in2, in3, in4, in5, in6, in7;
  int16x8_t out0, out1, out2, out3, out4, out5, out6, out7;
  int16x8_t s1_0, s1_1, s1_2, s1_3, s1_4, s1_5, s1_6, s1_7, s1_8, s1_9, s1_10,
      s1_11, s1_12, s1_13, s1_14, s1_15, s1_16, s1_17, s1_18, s1_19, s1_20,
      s1_21, s1_22, s1_23, s1_24, s1_25, s1_26, s1_27, s1_28, s1_29, s1_30,
      s1_31;
  int16x8_t s2_0, s2_1, s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10,
      s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20,
      s2_21, s2_22, s2_23, s2_24, s2_25, s2_26, s2_27, s2_28, s2_29, s2_30,
      s2_31;
  int16x8_t s3_24, s3_25, s3_26, s3_27;

  in0 = vld1q_s16(input);
  input += 8;
  in1 = vld1q_s16(input);
  input += 8;
  in2 = vld1q_s16(input);
  input += 8;
  in3 = vld1q_s16(input);
  input += 8;
  in4 = vld1q_s16(input);
  input += 8;
  in5 = vld1q_s16(input);
  input += 8;
  in6 = vld1q_s16(input);
  input += 8;
  in7 = vld1q_s16(input);

  transpose_s16_8x8(&in0, &in1, &in2, &in3, &in4, &in5, &in6, &in7);

  // stage 1
  s1_16 = multiply_shift_and_narrow(in1, cospi_31_64);
  s1_31 = multiply_shift_and_narrow(in1, cospi_1_64);

  // Different for _8_
  s1_19 = multiply_shift_and_narrow(in7, -cospi_25_64);
  s1_28 = multiply_shift_and_narrow(in7, cospi_7_64);

  s1_20 = multiply_shift_and_narrow(in5, cospi_27_64);
  s1_27 = multiply_shift_and_narrow(in5, cospi_5_64);

  s1_23 = multiply_shift_and_narrow(in3, -cospi_29_64);
  s1_24 = multiply_shift_and_narrow(in3, cospi_3_64);

  // stage 2
  s2_8 = multiply_shift_and_narrow(in2, cospi_30_64);
  s2_15 = multiply_shift_and_narrow(in2, cospi_2_64);

  s2_11 = multiply_shift_and_narrow(in6, -cospi_26_64);
  s2_12 = multiply_shift_and_narrow(in6, cospi_6_64);

  // stage 3
  s1_4 = multiply_shift_and_narrow(in4, cospi_28_64);
  s1_7 = multiply_shift_and_narrow(in4, cospi_4_64);

  s1_17 = multiply_accumulate_shift_and_narrow(s1_16, -cospi_4_64, s1_31,
                                               cospi_28_64);
  s1_30 = multiply_accumulate_shift_and_narrow(s1_16, cospi_28_64, s1_31,
                                               cospi_4_64);

  // Different for _8_
  s1_18 = multiply_accumulate_shift_and_narrow(s1_19, -cospi_28_64, s1_28,
                                               -cospi_4_64);
  s1_29 = multiply_accumulate_shift_and_narrow(s1_19, -cospi_4_64, s1_28,
                                               cospi_28_64);

  s1_21 = multiply_accumulate_shift_and_narrow(s1_20, -cospi_20_64, s1_27,
                                               cospi_12_64);
  s1_26 = multiply_accumulate_shift_and_narrow(s1_20, cospi_12_64, s1_27,
                                               cospi_20_64);

  s1_22 = multiply_accumulate_shift_and_narrow(s1_23, -cospi_12_64, s1_24,
                                               -cospi_20_64);
  s1_25 = multiply_accumulate_shift_and_narrow(s1_23, -cospi_20_64, s1_24,
                                               cospi_12_64);

  // stage 4
  s1_0 = multiply_shift_and_narrow(in0, cospi_16_64);

  s2_9 = multiply_accumulate_shift_and_narrow(s2_8, -cospi_8_64, s2_15,
                                              cospi_24_64);
  s2_14 = multiply_accumulate_shift_and_narrow(s2_8, cospi_24_64, s2_15,
                                               cospi_8_64);

  s2_10 = multiply_accumulate_shift_and_narrow(s2_11, -cospi_24_64, s2_12,
                                               -cospi_8_64);
  s2_13 = multiply_accumulate_shift_and_narrow(s2_11, -cospi_8_64, s2_12,
                                               cospi_24_64);

  s2_16 = vaddq_s16(s1_16, s1_19);

  s2_17 = vaddq_s16(s1_17, s1_18);
  s2_18 = vsubq_s16(s1_17, s1_18);

  s2_19 = vsubq_s16(s1_16, s1_19);

  s2_20 = vsubq_s16(s1_23, s1_20);
  s2_21 = vsubq_s16(s1_22, s1_21);

  s2_22 = vaddq_s16(s1_21, s1_22);
  s2_23 = vaddq_s16(s1_20, s1_23);

  s2_24 = vaddq_s16(s1_24, s1_27);
  s2_25 = vaddq_s16(s1_25, s1_26);
  s2_26 = vsubq_s16(s1_25, s1_26);
  s2_27 = vsubq_s16(s1_24, s1_27);

  s2_28 = vsubq_s16(s1_31, s1_28);
  s2_29 = vsubq_s16(s1_30, s1_29);
  s2_30 = vaddq_s16(s1_29, s1_30);
  s2_31 = vaddq_s16(s1_28, s1_31);

  // stage 5
  s1_5 = sub_multiply_shift_and_narrow(s1_7, s1_4, cospi_16_64);
  s1_6 = add_multiply_shift_and_narrow(s1_4, s1_7, cospi_16_64);

  s1_8 = vaddq_s16(s2_8, s2_11);
  s1_9 = vaddq_s16(s2_9, s2_10);
  s1_10 = vsubq_s16(s2_9, s2_10);
  s1_11 = vsubq_s16(s2_8, s2_11);
  s1_12 = vsubq_s16(s2_15, s2_12);
  s1_13 = vsubq_s16(s2_14, s2_13);
  s1_14 = vaddq_s16(s2_13, s2_14);
  s1_15 = vaddq_s16(s2_12, s2_15);

  s1_18 = multiply_accumulate_shift_and_narrow(s2_18, -cospi_8_64, s2_29,
                                               cospi_24_64);
  s1_29 = multiply_accumulate_shift_and_narrow(s2_18, cospi_24_64, s2_29,
                                               cospi_8_64);

  s1_19 = multiply_accumulate_shift_and_narrow(s2_19, -cospi_8_64, s2_28,
                                               cospi_24_64);
  s1_28 = multiply_accumulate_shift_and_narrow(s2_19, cospi_24_64, s2_28,
                                               cospi_8_64);

  s1_20 = multiply_accumulate_shift_and_narrow(s2_20, -cospi_24_64, s2_27,
                                               -cospi_8_64);
  s1_27 = multiply_accumulate_shift_and_narrow(s2_20, -cospi_8_64, s2_27,
                                               cospi_24_64);

  s1_21 = multiply_accumulate_shift_and_narrow(s2_21, -cospi_24_64, s2_26,
                                               -cospi_8_64);
  s1_26 = multiply_accumulate_shift_and_narrow(s2_21, -cospi_8_64, s2_26,
                                               cospi_24_64);

  // stage 6
  s2_0 = vaddq_s16(s1_0, s1_7);
  s2_1 = vaddq_s16(s1_0, s1_6);
  s2_2 = vaddq_s16(s1_0, s1_5);
  s2_3 = vaddq_s16(s1_0, s1_4);
  s2_4 = vsubq_s16(s1_0, s1_4);
  s2_5 = vsubq_s16(s1_0, s1_5);
  s2_6 = vsubq_s16(s1_0, s1_6);
  s2_7 = vsubq_s16(s1_0, s1_7);

  s2_10 = sub_multiply_shift_and_narrow(s1_13, s1_10, cospi_16_64);
  s2_13 = add_multiply_shift_and_narrow(s1_10, s1_13, cospi_16_64);

  s2_11 = sub_multiply_shift_and_narrow(s1_12, s1_11, cospi_16_64);
  s2_12 = add_multiply_shift_and_narrow(s1_11, s1_12, cospi_16_64);

  s1_16 = vaddq_s16(s2_16, s2_23);
  s1_17 = vaddq_s16(s2_17, s2_22);
  s2_18 = vaddq_s16(s1_18, s1_21);
  s2_19 = vaddq_s16(s1_19, s1_20);
  s2_20 = vsubq_s16(s1_19, s1_20);
  s2_21 = vsubq_s16(s1_18, s1_21);
  s1_22 = vsubq_s16(s2_17, s2_22);
  s1_23 = vsubq_s16(s2_16, s2_23);

  s3_24 = vsubq_s16(s2_31, s2_24);
  s3_25 = vsubq_s16(s2_30, s2_25);
  s3_26 = vsubq_s16(s1_29, s1_26);
  s3_27 = vsubq_s16(s1_28, s1_27);
  s2_28 = vaddq_s16(s1_27, s1_28);
  s2_29 = vaddq_s16(s1_26, s1_29);
  s2_30 = vaddq_s16(s2_25, s2_30);
  s2_31 = vaddq_s16(s2_24, s2_31);

  // stage 7
  s1_0 = vaddq_s16(s2_0, s1_15);
  s1_1 = vaddq_s16(s2_1, s1_14);
  s1_2 = vaddq_s16(s2_2, s2_13);
  s1_3 = vaddq_s16(s2_3, s2_12);
  s1_4 = vaddq_s16(s2_4, s2_11);
  s1_5 = vaddq_s16(s2_5, s2_10);
  s1_6 = vaddq_s16(s2_6, s1_9);
  s1_7 = vaddq_s16(s2_7, s1_8);
  s1_8 = vsubq_s16(s2_7, s1_8);
  s1_9 = vsubq_s16(s2_6, s1_9);
  s1_10 = vsubq_s16(s2_5, s2_10);
  s1_11 = vsubq_s16(s2_4, s2_11);
  s1_12 = vsubq_s16(s2_3, s2_12);
  s1_13 = vsubq_s16(s2_2, s2_13);
  s1_14 = vsubq_s16(s2_1, s1_14);
  s1_15 = vsubq_s16(s2_0, s1_15);

  s1_20 = sub_multiply_shift_and_narrow(s3_27, s2_20, cospi_16_64);
  s1_27 = add_multiply_shift_and_narrow(s2_20, s3_27, cospi_16_64);

  s1_21 = sub_multiply_shift_and_narrow(s3_26, s2_21, cospi_16_64);
  s1_26 = add_multiply_shift_and_narrow(s2_21, s3_26, cospi_16_64);

  s2_22 = sub_multiply_shift_and_narrow(s3_25, s1_22, cospi_16_64);
  s1_25 = add_multiply_shift_and_narrow(s1_22, s3_25, cospi_16_64);

  s2_23 = sub_multiply_shift_and_narrow(s3_24, s1_23, cospi_16_64);
  s1_24 = add_multiply_shift_and_narrow(s1_23, s3_24, cospi_16_64);

  // final stage
  out0 = vaddq_s16(s1_0, s2_31);
  out1 = vaddq_s16(s1_1, s2_30);
  out2 = vaddq_s16(s1_2, s2_29);
  out3 = vaddq_s16(s1_3, s2_28);
  out4 = vaddq_s16(s1_4, s1_27);
  out5 = vaddq_s16(s1_5, s1_26);
  out6 = vaddq_s16(s1_6, s1_25);
  out7 = vaddq_s16(s1_7, s1_24);

  add_and_store(out0, out1, out2, out3, out4, out5, out6, out7, output, stride);

  out0 = vaddq_s16(s1_8, s2_23);
  out1 = vaddq_s16(s1_9, s2_22);
  out2 = vaddq_s16(s1_10, s1_21);
  out3 = vaddq_s16(s1_11, s1_20);
  out4 = vaddq_s16(s1_12, s2_19);
  out5 = vaddq_s16(s1_13, s2_18);
  out6 = vaddq_s16(s1_14, s1_17);
  out7 = vaddq_s16(s1_15, s1_16);

  add_and_store(out0, out1, out2, out3, out4, out5, out6, out7,
                output + (8 * stride), stride);

  out0 = vsubq_s16(s1_15, s1_16);
  out1 = vsubq_s16(s1_14, s1_17);
  out2 = vsubq_s16(s1_13, s2_18);
  out3 = vsubq_s16(s1_12, s2_19);
  out4 = vsubq_s16(s1_11, s1_20);
  out5 = vsubq_s16(s1_10, s1_21);
  out6 = vsubq_s16(s1_9, s2_22);
  out7 = vsubq_s16(s1_8, s2_23);

  add_and_store(out0, out1, out2, out3, out4, out5, out6, out7,
                output + (16 * stride), stride);

  out0 = vsubq_s16(s1_7, s1_24);
  out1 = vsubq_s16(s1_6, s1_25);
  out2 = vsubq_s16(s1_5, s1_26);
  out3 = vsubq_s16(s1_4, s1_27);
  out4 = vsubq_s16(s1_3, s2_28);
  out5 = vsubq_s16(s1_2, s2_29);
  out6 = vsubq_s16(s1_1, s2_30);
  out7 = vsubq_s16(s1_0, s2_31);

  add_and_store(out0, out1, out2, out3, out4, out5, out6, out7,
                output + (24 * stride), stride);
}

void vpx_idct32x32_34_add_neon(const int16_t *input, uint8_t *dest,
                               int stride) {
  int i;
  int16_t temp[32 * 8];
  int16_t *t = temp;

  idct32_6_neon(input, t);

  for (i = 0; i < 32; i += 8) {
    idct32_8_neon(t, dest, stride);
    t += (8 * 8);
    dest += 8;
  }
}
