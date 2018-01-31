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

#include "./vp9_rtcd.h"
#include "./vpx_config.h"
#include "vp9/common/vp9_common.h"
#include "vpx_dsp/arm/idct_neon.h"
#include "vpx_dsp/arm/mem_neon.h"
#include "vpx_dsp/arm/transpose_neon.h"

static INLINE void iadst_half_butterfly_neon(int16x8_t *const x) {
  const int16x4_t c = vdup_n_s16(cospi_16_64);
  const int16x8_t sum = vaddq_s16(x[0], x[1]);
  const int16x8_t sub = vsubq_s16(x[0], x[1]);
  int32x4_t t0[2], t1[2];

  t0[0] = vmull_s16(c, vget_low_s16(sum));
  t0[1] = vmull_s16(c, vget_high_s16(sum));
  t1[0] = vmull_s16(c, vget_low_s16(sub));
  t1[1] = vmull_s16(c, vget_high_s16(sub));
  x[0] = dct_const_round_shift_low_8(t0);
  x[1] = dct_const_round_shift_low_8(t1);
}

static INLINE void iadst_butterfly_neon(const int16x8_t in0,
                                        const int16x8_t in1, const int c0,
                                        const int c1, int32x4_t *const s0,
                                        int32x4_t *const s1) {
  const int16x4_t cst0 = vdup_n_s16(c0);
  const int16x4_t cst1 = vdup_n_s16(c1);
  int32x4_t t0[2], t1[2];

  t0[0] = vmull_s16(cst0, vget_low_s16(in0));
  t0[1] = vmull_s16(cst0, vget_high_s16(in0));
  t1[0] = vmull_s16(cst1, vget_low_s16(in0));
  t1[1] = vmull_s16(cst1, vget_high_s16(in0));

  s0[0] = vmlal_s16(t0[0], cst1, vget_low_s16(in1));
  s0[1] = vmlal_s16(t0[1], cst1, vget_high_s16(in1));
  s1[0] = vmlsl_s16(t1[0], cst0, vget_low_s16(in1));
  s1[1] = vmlsl_s16(t1[1], cst0, vget_high_s16(in1));
}

static INLINE int16x8_t add_dct_const_round_shift_low_8(
    const int32x4_t *const in0, const int32x4_t *const in1) {
  int32x4_t sum[2];

  sum[0] = vaddq_s32(in0[0], in1[0]);
  sum[1] = vaddq_s32(in0[1], in1[1]);
  return dct_const_round_shift_low_8(sum);
}

static INLINE int16x8_t sub_dct_const_round_shift_low_8(
    const int32x4_t *const in0, const int32x4_t *const in1) {
  int32x4_t sum[2];

  sum[0] = vsubq_s32(in0[0], in1[0]);
  sum[1] = vsubq_s32(in0[1], in1[1]);
  return dct_const_round_shift_low_8(sum);
}

static INLINE void iadst8(int16x8_t *const io) {
  int16x8_t x[8], t[4];
  int32x4_t s0[2], s1[2], s2[2], s3[2], s4[2], s5[2], s6[2], s7[2];

  x[0] = io[7];
  x[1] = io[0];
  x[2] = io[5];
  x[3] = io[2];
  x[4] = io[3];
  x[5] = io[4];
  x[6] = io[1];
  x[7] = io[6];

  // stage 1
  iadst_butterfly_neon(x[0], x[1], cospi_2_64, cospi_30_64, s0, s1);
  iadst_butterfly_neon(x[2], x[3], cospi_10_64, cospi_22_64, s2, s3);
  iadst_butterfly_neon(x[4], x[5], cospi_18_64, cospi_14_64, s4, s5);
  iadst_butterfly_neon(x[6], x[7], cospi_26_64, cospi_6_64, s6, s7);

  x[0] = add_dct_const_round_shift_low_8(s0, s4);
  x[1] = add_dct_const_round_shift_low_8(s1, s5);
  x[2] = add_dct_const_round_shift_low_8(s2, s6);
  x[3] = add_dct_const_round_shift_low_8(s3, s7);
  x[4] = sub_dct_const_round_shift_low_8(s0, s4);
  x[5] = sub_dct_const_round_shift_low_8(s1, s5);
  x[6] = sub_dct_const_round_shift_low_8(s2, s6);
  x[7] = sub_dct_const_round_shift_low_8(s3, s7);

  // stage 2
  t[0] = x[0];
  t[1] = x[1];
  t[2] = x[2];
  t[3] = x[3];
  iadst_butterfly_neon(x[4], x[5], cospi_8_64, cospi_24_64, s4, s5);
  iadst_butterfly_neon(x[7], x[6], cospi_24_64, cospi_8_64, s7, s6);

  x[0] = vaddq_s16(t[0], t[2]);
  x[1] = vaddq_s16(t[1], t[3]);
  x[2] = vsubq_s16(t[0], t[2]);
  x[3] = vsubq_s16(t[1], t[3]);
  x[4] = add_dct_const_round_shift_low_8(s4, s6);
  x[5] = add_dct_const_round_shift_low_8(s5, s7);
  x[6] = sub_dct_const_round_shift_low_8(s4, s6);
  x[7] = sub_dct_const_round_shift_low_8(s5, s7);

  // stage 3
  iadst_half_butterfly_neon(x + 2);
  iadst_half_butterfly_neon(x + 6);

  io[0] = x[0];
  io[1] = vnegq_s16(x[4]);
  io[2] = x[6];
  io[3] = vnegq_s16(x[2]);
  io[4] = x[3];
  io[5] = vnegq_s16(x[7]);
  io[6] = x[5];
  io[7] = vnegq_s16(x[1]);
}

void vp9_iht8x8_64_add_neon(const tran_low_t *input, uint8_t *dest, int stride,
                            int tx_type) {
  const int16x8_t cospis = vld1q_s16(kCospi);
  const int16x4_t cospis0 = vget_low_s16(cospis);   // cospi 0, 8, 16, 24
  const int16x4_t cospis1 = vget_high_s16(cospis);  // cospi 4, 12, 20, 28
  int16x8_t a[8];

  a[0] = load_tran_low_to_s16q(input + 0 * 8);
  a[1] = load_tran_low_to_s16q(input + 1 * 8);
  a[2] = load_tran_low_to_s16q(input + 2 * 8);
  a[3] = load_tran_low_to_s16q(input + 3 * 8);
  a[4] = load_tran_low_to_s16q(input + 4 * 8);
  a[5] = load_tran_low_to_s16q(input + 5 * 8);
  a[6] = load_tran_low_to_s16q(input + 6 * 8);
  a[7] = load_tran_low_to_s16q(input + 7 * 8);

  transpose_s16_8x8(&a[0], &a[1], &a[2], &a[3], &a[4], &a[5], &a[6], &a[7]);

  switch (tx_type) {
    case DCT_DCT:
      idct8x8_64_1d_bd8_kernel(cospis0, cospis1, a);
      transpose_s16_8x8(&a[0], &a[1], &a[2], &a[3], &a[4], &a[5], &a[6], &a[7]);
      idct8x8_64_1d_bd8_kernel(cospis0, cospis1, a);
      break;

    case ADST_DCT:
      idct8x8_64_1d_bd8_kernel(cospis0, cospis1, a);
      transpose_s16_8x8(&a[0], &a[1], &a[2], &a[3], &a[4], &a[5], &a[6], &a[7]);
      iadst8(a);
      break;

    case DCT_ADST:
      iadst8(a);
      transpose_s16_8x8(&a[0], &a[1], &a[2], &a[3], &a[4], &a[5], &a[6], &a[7]);
      idct8x8_64_1d_bd8_kernel(cospis0, cospis1, a);
      break;

    default:
      assert(tx_type == ADST_ADST);
      iadst8(a);
      transpose_s16_8x8(&a[0], &a[1], &a[2], &a[3], &a[4], &a[5], &a[6], &a[7]);
      iadst8(a);
      break;
  }

  idct8x8_add8x8_neon(a, dest, stride);
}
