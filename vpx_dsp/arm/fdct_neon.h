/*
 *  Copyright (c) 2022 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#ifndef VPX_VPX_DSP_ARM_FDCT_NEON_H_
#define VPX_VPX_DSP_ARM_FDCT_NEON_H_

#include <arm_neon.h>

// fdct_round_shift((a +/- b) * c)
static INLINE void butterfly_one_coeff(const int16x8_t a, const int16x8_t b,
                                       const tran_high_t constant,
                                       int16x8_t *add, int16x8_t *sub) {
  const int32x4_t a0 = vmull_n_s16(vget_low_s16(a), constant);
  const int32x4_t a1 = vmull_n_s16(vget_high_s16(a), constant);
  const int32x4_t sum0 = vmlal_n_s16(a0, vget_low_s16(b), constant);
  const int32x4_t sum1 = vmlal_n_s16(a1, vget_high_s16(b), constant);
  const int32x4_t diff0 = vmlsl_n_s16(a0, vget_low_s16(b), constant);
  const int32x4_t diff1 = vmlsl_n_s16(a1, vget_high_s16(b), constant);
  const int16x4_t rounded0 = vqrshrn_n_s32(sum0, DCT_CONST_BITS);
  const int16x4_t rounded1 = vqrshrn_n_s32(sum1, DCT_CONST_BITS);
  const int16x4_t rounded2 = vqrshrn_n_s32(diff0, DCT_CONST_BITS);
  const int16x4_t rounded3 = vqrshrn_n_s32(diff1, DCT_CONST_BITS);
  *add = vcombine_s16(rounded0, rounded1);
  *sub = vcombine_s16(rounded2, rounded3);
}

// fdct_round_shift(a * c0 +/- b * c1)
static INLINE void butterfly_two_coeff(const int16x8_t a, const int16x8_t b,
                                       const tran_coef_t constant0,
                                       const tran_coef_t constant1,
                                       int16x8_t *add, int16x8_t *sub) {
  const int32x4_t a0 = vmull_n_s16(vget_low_s16(a), constant0);
  const int32x4_t a1 = vmull_n_s16(vget_high_s16(a), constant0);
  const int32x4_t a2 = vmull_n_s16(vget_low_s16(a), constant1);
  const int32x4_t a3 = vmull_n_s16(vget_high_s16(a), constant1);
  const int32x4_t sum0 = vmlal_n_s16(a2, vget_low_s16(b), constant0);
  const int32x4_t sum1 = vmlal_n_s16(a3, vget_high_s16(b), constant0);
  const int32x4_t diff0 = vmlsl_n_s16(a0, vget_low_s16(b), constant1);
  const int32x4_t diff1 = vmlsl_n_s16(a1, vget_high_s16(b), constant1);
  const int16x4_t rounded0 = vqrshrn_n_s32(sum0, DCT_CONST_BITS);
  const int16x4_t rounded1 = vqrshrn_n_s32(sum1, DCT_CONST_BITS);
  const int16x4_t rounded2 = vqrshrn_n_s32(diff0, DCT_CONST_BITS);
  const int16x4_t rounded3 = vqrshrn_n_s32(diff1, DCT_CONST_BITS);
  *add = vcombine_s16(rounded0, rounded1);
  *sub = vcombine_s16(rounded2, rounded3);
}

// Add 2 if positive, 1 if negative, and shift by 2.
// In practice, subtract the sign bit, then shift with rounding.
static INLINE int16x8_t sub_round_shift(const int16x8_t a) {
  const uint16x8_t a_u16 = vreinterpretq_u16_s16(a);
  const uint16x8_t a_sign_u16 = vshrq_n_u16(a_u16, 15);
  const int16x8_t a_sign_s16 = vreinterpretq_s16_u16(a_sign_u16);
  return vrshrq_n_s16(vsubq_s16(a, a_sign_s16), 2);
}

// Like butterfly_one_coeff, but don't narrow results.
static INLINE void butterfly_one_coeff_s16_s32(
    const int16x8_t a, const int16x8_t b, const tran_high_t constant,
    int32x4_t *add_lo, int32x4_t *add_hi, int32x4_t *sub_lo,
    int32x4_t *sub_hi) {
  const int32x4_t a0 = vmull_n_s16(vget_low_s16(a), constant);
  const int32x4_t a1 = vmull_n_s16(vget_high_s16(a), constant);
  const int32x4_t sum0 = vmlal_n_s16(a0, vget_low_s16(b), constant);
  const int32x4_t sum1 = vmlal_n_s16(a1, vget_high_s16(b), constant);
  const int32x4_t diff0 = vmlsl_n_s16(a0, vget_low_s16(b), constant);
  const int32x4_t diff1 = vmlsl_n_s16(a1, vget_high_s16(b), constant);
  *add_lo = vrshrq_n_s32(sum0, DCT_CONST_BITS);
  *add_hi = vrshrq_n_s32(sum1, DCT_CONST_BITS);
  *sub_lo = vrshrq_n_s32(diff0, DCT_CONST_BITS);
  *sub_hi = vrshrq_n_s32(diff1, DCT_CONST_BITS);
}

// Like butterfly_one_coeff, but with s32.
static INLINE void butterfly_one_coeff_s32(
    const int32x4_t a_lo, const int32x4_t a_hi, const int32x4_t b_lo,
    const int32x4_t b_hi, const int32_t constant, int32x4_t *add_lo,
    int32x4_t *add_hi, int32x4_t *sub_lo, int32x4_t *sub_hi) {
  const int32x4_t a_lo_0 = vmulq_n_s32(a_lo, constant);
  const int32x4_t a_hi_0 = vmulq_n_s32(a_hi, constant);
  const int32x4_t sum0 = vmlaq_n_s32(a_lo_0, b_lo, constant);
  const int32x4_t sum1 = vmlaq_n_s32(a_hi_0, b_hi, constant);
  const int32x4_t diff0 = vmlsq_n_s32(a_lo_0, b_lo, constant);
  const int32x4_t diff1 = vmlsq_n_s32(a_hi_0, b_hi, constant);
  *add_lo = vrshrq_n_s32(sum0, DCT_CONST_BITS);
  *add_hi = vrshrq_n_s32(sum1, DCT_CONST_BITS);
  *sub_lo = vrshrq_n_s32(diff0, DCT_CONST_BITS);
  *sub_hi = vrshrq_n_s32(diff1, DCT_CONST_BITS);
}

// Like butterfly_two_coeff, but with s32.
static INLINE void butterfly_two_coeff_s32(
    const int32x4_t a_lo, const int32x4_t a_hi, const int32x4_t b_lo,
    const int32x4_t b_hi, const int32_t constant0, const int32_t constant1,
    int32x4_t *add_lo, int32x4_t *add_hi, int32x4_t *sub_lo,
    int32x4_t *sub_hi) {
  const int32x4_t a0 = vmulq_n_s32(a_lo, constant0);
  const int32x4_t a1 = vmulq_n_s32(a_hi, constant0);
  const int32x4_t a2 = vmulq_n_s32(a_lo, constant1);
  const int32x4_t a3 = vmulq_n_s32(a_hi, constant1);
  const int32x4_t sum0 = vmlaq_n_s32(a2, b_lo, constant0);
  const int32x4_t sum1 = vmlaq_n_s32(a3, b_hi, constant0);
  const int32x4_t diff0 = vmlsq_n_s32(a0, b_lo, constant1);
  const int32x4_t diff1 = vmlsq_n_s32(a1, b_hi, constant1);
  *add_lo = vrshrq_n_s32(sum0, DCT_CONST_BITS);
  *add_hi = vrshrq_n_s32(sum1, DCT_CONST_BITS);
  *sub_lo = vrshrq_n_s32(diff0, DCT_CONST_BITS);
  *sub_hi = vrshrq_n_s32(diff1, DCT_CONST_BITS);
}

// Add 1 if positive, 2 if negative, and shift by 2.
// In practice, add 1, then add the sign bit, then shift without rounding.
static INLINE int16x8_t add_round_shift_s16(const int16x8_t a) {
  const int16x8_t one = vdupq_n_s16(1);
  const uint16x8_t a_u16 = vreinterpretq_u16_s16(a);
  const uint16x8_t a_sign_u16 = vshrq_n_u16(a_u16, 15);
  const int16x8_t a_sign_s16 = vreinterpretq_s16_u16(a_sign_u16);
  return vshrq_n_s16(vaddq_s16(vaddq_s16(a, a_sign_s16), one), 2);
}

// Add 1 if positive, 2 if negative, and shift by 2.
// In practice, add 1, then add the sign bit, then shift without rounding.
static INLINE int16x8_t add_round_shift_s32(const int32x4_t a_lo,
                                            const int32x4_t a_hi) {
  const int32x4_t one = vdupq_n_s32(1);
  const uint32x4_t a_lo_u32 = vreinterpretq_u32_s32(a_lo);
  const uint32x4_t a_lo_sign_u32 = vshrq_n_u32(a_lo_u32, 31);
  const int32x4_t a_lo_sign_s32 = vreinterpretq_s32_u32(a_lo_sign_u32);
  const int16x4_t b_lo =
      vshrn_n_s32(vqaddq_s32(vqaddq_s32(a_lo, a_lo_sign_s32), one), 2);
  const uint32x4_t a_hi_u32 = vreinterpretq_u32_s32(a_hi);
  const uint32x4_t a_hi_sign_u32 = vshrq_n_u32(a_hi_u32, 31);
  const int32x4_t a_hi_sign_s32 = vreinterpretq_s32_u32(a_hi_sign_u32);
  const int16x4_t b_hi =
      vshrn_n_s32(vqaddq_s32(vqaddq_s32(a_hi, a_hi_sign_s32), one), 2);
  return vcombine_s16(b_lo, b_hi);
}

static INLINE void vpx_fdct4x4_pass1_neon(int16x4_t *in) {
  const int16x8_t input_01 = vcombine_s16(in[0], in[1]);
  const int16x8_t input_32 = vcombine_s16(in[3], in[2]);

  // in_0 +/- in_3, in_1 +/- in_2
  const int16x8_t s_01 = vaddq_s16(input_01, input_32);
  const int16x8_t s_32 = vsubq_s16(input_01, input_32);

  // step_0 +/- step_1, step_2 +/- step_3
  const int16x4_t s_0 = vget_low_s16(s_01);
  const int16x4_t s_1 = vget_high_s16(s_01);
  const int16x4_t s_2 = vget_high_s16(s_32);
  const int16x4_t s_3 = vget_low_s16(s_32);

  // (s_0 +/- s_1) * cospi_16_64
  // Must expand all elements to s32. See 'needs32' comment in fwd_txfm.c.
  const int32x4_t s_0_p_s_1 = vaddl_s16(s_0, s_1);
  const int32x4_t s_0_m_s_1 = vsubl_s16(s_0, s_1);
  const int32x4_t temp1 = vmulq_n_s32(s_0_p_s_1, cospi_16_64);
  const int32x4_t temp2 = vmulq_n_s32(s_0_m_s_1, cospi_16_64);

  // fdct_round_shift
  int16x4_t out_0 = vrshrn_n_s32(temp1, DCT_CONST_BITS);
  int16x4_t out_2 = vrshrn_n_s32(temp2, DCT_CONST_BITS);

  // s_3 * cospi_8_64 + s_2 * cospi_24_64
  // s_3 * cospi_24_64 - s_2 * cospi_8_64
  const int32x4_t s_3_cospi_8_64 = vmull_n_s16(s_3, cospi_8_64);
  const int32x4_t s_3_cospi_24_64 = vmull_n_s16(s_3, cospi_24_64);

  const int32x4_t temp3 = vmlal_n_s16(s_3_cospi_8_64, s_2, cospi_24_64);
  const int32x4_t temp4 = vmlsl_n_s16(s_3_cospi_24_64, s_2, cospi_8_64);

  // fdct_round_shift
  int16x4_t out_1 = vrshrn_n_s32(temp3, DCT_CONST_BITS);
  int16x4_t out_3 = vrshrn_n_s32(temp4, DCT_CONST_BITS);

  transpose_s16_4x4d(&out_0, &out_1, &out_2, &out_3);

  in[0] = out_0;
  in[1] = out_1;
  in[2] = out_2;
  in[3] = out_3;
}

static INLINE void vpx_fdct8x8_pass1_notranspose_neon(int16x8_t *in,
                                                      int16x8_t *out) {
  const int16x8_t v_s0 = vaddq_s16(in[0], in[7]);
  const int16x8_t v_s1 = vaddq_s16(in[1], in[6]);
  const int16x8_t v_s2 = vaddq_s16(in[2], in[5]);
  const int16x8_t v_s3 = vaddq_s16(in[3], in[4]);
  const int16x8_t v_s4 = vsubq_s16(in[3], in[4]);
  const int16x8_t v_s5 = vsubq_s16(in[2], in[5]);
  const int16x8_t v_s6 = vsubq_s16(in[1], in[6]);
  const int16x8_t v_s7 = vsubq_s16(in[0], in[7]);
  // fdct4(step, step);
  int16x8_t v_x0 = vaddq_s16(v_s0, v_s3);
  int16x8_t v_x1 = vaddq_s16(v_s1, v_s2);
  int16x8_t v_x2 = vsubq_s16(v_s1, v_s2);
  int16x8_t v_x3 = vsubq_s16(v_s0, v_s3);
  // fdct4(step, step);
  int32x4_t v_t0_lo = vaddl_s16(vget_low_s16(v_x0), vget_low_s16(v_x1));
  int32x4_t v_t0_hi = vaddl_s16(vget_high_s16(v_x0), vget_high_s16(v_x1));
  int32x4_t v_t1_lo = vsubl_s16(vget_low_s16(v_x0), vget_low_s16(v_x1));
  int32x4_t v_t1_hi = vsubl_s16(vget_high_s16(v_x0), vget_high_s16(v_x1));
  int32x4_t v_t2_lo = vmull_n_s16(vget_low_s16(v_x2), cospi_24_64);
  int32x4_t v_t2_hi = vmull_n_s16(vget_high_s16(v_x2), cospi_24_64);
  int32x4_t v_t3_lo = vmull_n_s16(vget_low_s16(v_x3), cospi_24_64);
  int32x4_t v_t3_hi = vmull_n_s16(vget_high_s16(v_x3), cospi_24_64);
  v_t2_lo = vmlal_n_s16(v_t2_lo, vget_low_s16(v_x3), cospi_8_64);
  v_t2_hi = vmlal_n_s16(v_t2_hi, vget_high_s16(v_x3), cospi_8_64);
  v_t3_lo = vmlsl_n_s16(v_t3_lo, vget_low_s16(v_x2), cospi_8_64);
  v_t3_hi = vmlsl_n_s16(v_t3_hi, vget_high_s16(v_x2), cospi_8_64);
  v_t0_lo = vmulq_n_s32(v_t0_lo, cospi_16_64);
  v_t0_hi = vmulq_n_s32(v_t0_hi, cospi_16_64);
  v_t1_lo = vmulq_n_s32(v_t1_lo, cospi_16_64);
  v_t1_hi = vmulq_n_s32(v_t1_hi, cospi_16_64);
  {
    const int16x4_t a = vrshrn_n_s32(v_t0_lo, DCT_CONST_BITS);
    const int16x4_t b = vrshrn_n_s32(v_t0_hi, DCT_CONST_BITS);
    const int16x4_t c = vrshrn_n_s32(v_t1_lo, DCT_CONST_BITS);
    const int16x4_t d = vrshrn_n_s32(v_t1_hi, DCT_CONST_BITS);
    const int16x4_t e = vrshrn_n_s32(v_t2_lo, DCT_CONST_BITS);
    const int16x4_t f = vrshrn_n_s32(v_t2_hi, DCT_CONST_BITS);
    const int16x4_t g = vrshrn_n_s32(v_t3_lo, DCT_CONST_BITS);
    const int16x4_t h = vrshrn_n_s32(v_t3_hi, DCT_CONST_BITS);
    out[0] = vcombine_s16(a, c);  // 00 01 02 03 40 41 42 43
    out[2] = vcombine_s16(e, g);  // 20 21 22 23 60 61 62 63
    out[4] = vcombine_s16(b, d);  // 04 05 06 07 44 45 46 47
    out[6] = vcombine_s16(f, h);  // 24 25 26 27 64 65 66 67
  }
  // Stage 2
  v_x0 = vsubq_s16(v_s6, v_s5);
  v_x1 = vaddq_s16(v_s6, v_s5);
  v_t0_lo = vmull_n_s16(vget_low_s16(v_x0), cospi_16_64);
  v_t0_hi = vmull_n_s16(vget_high_s16(v_x0), cospi_16_64);
  v_t1_lo = vmull_n_s16(vget_low_s16(v_x1), cospi_16_64);
  v_t1_hi = vmull_n_s16(vget_high_s16(v_x1), cospi_16_64);
  {
    const int16x4_t a = vrshrn_n_s32(v_t0_lo, DCT_CONST_BITS);
    const int16x4_t b = vrshrn_n_s32(v_t0_hi, DCT_CONST_BITS);
    const int16x4_t c = vrshrn_n_s32(v_t1_lo, DCT_CONST_BITS);
    const int16x4_t d = vrshrn_n_s32(v_t1_hi, DCT_CONST_BITS);
    const int16x8_t ab = vcombine_s16(a, b);
    const int16x8_t cd = vcombine_s16(c, d);
    // Stage 3
    v_x0 = vaddq_s16(v_s4, ab);
    v_x1 = vsubq_s16(v_s4, ab);
    v_x2 = vsubq_s16(v_s7, cd);
    v_x3 = vaddq_s16(v_s7, cd);
  }
  // Stage 4
  v_t0_lo = vmull_n_s16(vget_low_s16(v_x3), cospi_4_64);
  v_t0_hi = vmull_n_s16(vget_high_s16(v_x3), cospi_4_64);
  v_t0_lo = vmlal_n_s16(v_t0_lo, vget_low_s16(v_x0), cospi_28_64);
  v_t0_hi = vmlal_n_s16(v_t0_hi, vget_high_s16(v_x0), cospi_28_64);
  v_t1_lo = vmull_n_s16(vget_low_s16(v_x1), cospi_12_64);
  v_t1_hi = vmull_n_s16(vget_high_s16(v_x1), cospi_12_64);
  v_t1_lo = vmlal_n_s16(v_t1_lo, vget_low_s16(v_x2), cospi_20_64);
  v_t1_hi = vmlal_n_s16(v_t1_hi, vget_high_s16(v_x2), cospi_20_64);
  v_t2_lo = vmull_n_s16(vget_low_s16(v_x2), cospi_12_64);
  v_t2_hi = vmull_n_s16(vget_high_s16(v_x2), cospi_12_64);
  v_t2_lo = vmlsl_n_s16(v_t2_lo, vget_low_s16(v_x1), cospi_20_64);
  v_t2_hi = vmlsl_n_s16(v_t2_hi, vget_high_s16(v_x1), cospi_20_64);
  v_t3_lo = vmull_n_s16(vget_low_s16(v_x3), cospi_28_64);
  v_t3_hi = vmull_n_s16(vget_high_s16(v_x3), cospi_28_64);
  v_t3_lo = vmlsl_n_s16(v_t3_lo, vget_low_s16(v_x0), cospi_4_64);
  v_t3_hi = vmlsl_n_s16(v_t3_hi, vget_high_s16(v_x0), cospi_4_64);
  {
    const int16x4_t a = vrshrn_n_s32(v_t0_lo, DCT_CONST_BITS);
    const int16x4_t b = vrshrn_n_s32(v_t0_hi, DCT_CONST_BITS);
    const int16x4_t c = vrshrn_n_s32(v_t1_lo, DCT_CONST_BITS);
    const int16x4_t d = vrshrn_n_s32(v_t1_hi, DCT_CONST_BITS);
    const int16x4_t e = vrshrn_n_s32(v_t2_lo, DCT_CONST_BITS);
    const int16x4_t f = vrshrn_n_s32(v_t2_hi, DCT_CONST_BITS);
    const int16x4_t g = vrshrn_n_s32(v_t3_lo, DCT_CONST_BITS);
    const int16x4_t h = vrshrn_n_s32(v_t3_hi, DCT_CONST_BITS);
    out[1] = vcombine_s16(a, c);  // 10 11 12 13 50 51 52 53
    out[3] = vcombine_s16(e, g);  // 30 31 32 33 70 71 72 73
    out[5] = vcombine_s16(b, d);  // 14 15 16 17 54 55 56 57
    out[7] = vcombine_s16(f, h);  // 34 35 36 37 74 75 76 77
  }
}

static INLINE void vpx_fdct8x8_pass1_neon(int16x8_t *in) {
  int16x8_t out[8];
  vpx_fdct8x8_pass1_notranspose_neon(in, out);
  // transpose 8x8
  // Can't use transpose_s16_8x8() because the values are arranged in two 4x8
  // columns.
  {
    // 00 01 02 03 40 41 42 43
    // 10 11 12 13 50 51 52 53
    // 20 21 22 23 60 61 62 63
    // 30 31 32 33 70 71 72 73
    // 04 05 06 07 44 45 46 47
    // 14 15 16 17 54 55 56 57
    // 24 25 26 27 64 65 66 67
    // 34 35 36 37 74 75 76 77
    const int32x4x2_t r02_s32 =
        vtrnq_s32(vreinterpretq_s32_s16(out[0]), vreinterpretq_s32_s16(out[2]));
    const int32x4x2_t r13_s32 =
        vtrnq_s32(vreinterpretq_s32_s16(out[1]), vreinterpretq_s32_s16(out[3]));
    const int32x4x2_t r46_s32 =
        vtrnq_s32(vreinterpretq_s32_s16(out[4]), vreinterpretq_s32_s16(out[6]));
    const int32x4x2_t r57_s32 =
        vtrnq_s32(vreinterpretq_s32_s16(out[5]), vreinterpretq_s32_s16(out[7]));
    const int16x8x2_t r01_s16 =
        vtrnq_s16(vreinterpretq_s16_s32(r02_s32.val[0]),
                  vreinterpretq_s16_s32(r13_s32.val[0]));
    const int16x8x2_t r23_s16 =
        vtrnq_s16(vreinterpretq_s16_s32(r02_s32.val[1]),
                  vreinterpretq_s16_s32(r13_s32.val[1]));
    const int16x8x2_t r45_s16 =
        vtrnq_s16(vreinterpretq_s16_s32(r46_s32.val[0]),
                  vreinterpretq_s16_s32(r57_s32.val[0]));
    const int16x8x2_t r67_s16 =
        vtrnq_s16(vreinterpretq_s16_s32(r46_s32.val[1]),
                  vreinterpretq_s16_s32(r57_s32.val[1]));
    in[0] = r01_s16.val[0];
    in[1] = r01_s16.val[1];
    in[2] = r23_s16.val[0];
    in[3] = r23_s16.val[1];
    in[4] = r45_s16.val[0];
    in[5] = r45_s16.val[1];
    in[6] = r67_s16.val[0];
    in[7] = r67_s16.val[1];
    // 00 10 20 30 40 50 60 70
    // 01 11 21 31 41 51 61 71
    // 02 12 22 32 42 52 62 72
    // 03 13 23 33 43 53 63 73
    // 04 14 24 34 44 54 64 74
    // 05 15 25 35 45 55 65 75
    // 06 16 26 36 46 56 66 76
    // 07 17 27 37 47 57 67 77
  }
}
#endif  // VPX_VPX_DSP_ARM_FDCT_NEON_H_
