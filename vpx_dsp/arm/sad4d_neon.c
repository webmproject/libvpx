/*
 *  Copyright (c) 2015 The WebM project authors. All Rights Reserved.
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
#include "vpx/vpx_integer.h"
#include "vpx_dsp/arm/mem_neon.h"
#include "vpx_dsp/arm/sum_neon.h"

void vpx_sad4x4x4d_neon(const uint8_t *src, int src_stride,
                        const uint8_t *const ref[4], int ref_stride,
                        uint32_t *res) {
  int i;
  const uint8x16_t src_u8 = load_unaligned_u8q(src, src_stride);
  for (i = 0; i < 4; ++i) {
    const uint8x16_t ref_u8 = load_unaligned_u8q(ref[i], ref_stride);
    uint16x8_t abs = vabdl_u8(vget_low_u8(src_u8), vget_low_u8(ref_u8));
    abs = vabal_u8(abs, vget_high_u8(src_u8), vget_high_u8(ref_u8));
    res[i] = vget_lane_u32(horizontal_add_uint16x8(abs), 0);
  }
}

void vpx_sad4x8x4d_neon(const uint8_t *src, int src_stride,
                        const uint8_t *const ref[4], int ref_stride,
                        uint32_t *res) {
  int i;
  const uint8x16_t src_0 = load_unaligned_u8q(src, src_stride);
  const uint8x16_t src_1 = load_unaligned_u8q(src + 4 * src_stride, src_stride);
  for (i = 0; i < 4; ++i) {
    const uint8x16_t ref_0 = load_unaligned_u8q(ref[i], ref_stride);
    const uint8x16_t ref_1 =
        load_unaligned_u8q(ref[i] + 4 * ref_stride, ref_stride);
    uint16x8_t abs = vabdl_u8(vget_low_u8(src_0), vget_low_u8(ref_0));
    abs = vabal_u8(abs, vget_high_u8(src_0), vget_high_u8(ref_0));
    abs = vabal_u8(abs, vget_low_u8(src_1), vget_low_u8(ref_1));
    abs = vabal_u8(abs, vget_high_u8(src_1), vget_high_u8(ref_1));
    res[i] = vget_lane_u32(horizontal_add_uint16x8(abs), 0);
  }
}

static INLINE void sad8x_4d(const uint8_t *a, int a_stride,
                            const uint8_t *const b[4], int b_stride,
                            uint32_t *result, const int height) {
  int i, j;
  uint16x8_t sum[4] = { vdupq_n_u16(0), vdupq_n_u16(0), vdupq_n_u16(0),
                        vdupq_n_u16(0) };
  const uint8_t *b_loop[4] = { b[0], b[1], b[2], b[3] };

  for (i = 0; i < height; ++i) {
    const uint8x8_t a_u8 = vld1_u8(a);
    a += a_stride;
    for (j = 0; j < 4; ++j) {
      const uint8x8_t b_u8 = vld1_u8(b_loop[j]);
      b_loop[j] += b_stride;
      sum[j] = vabal_u8(sum[j], a_u8, b_u8);
    }
  }

  for (j = 0; j < 4; ++j) {
    result[j] = vget_lane_u32(horizontal_add_uint16x8(sum[j]), 0);
  }
}

void vpx_sad8x4x4d_neon(const uint8_t *src, int src_stride,
                        const uint8_t *const ref[4], int ref_stride,
                        uint32_t *res) {
  sad8x_4d(src, src_stride, ref, ref_stride, res, 4);
}

void vpx_sad8x8x4d_neon(const uint8_t *src, int src_stride,
                        const uint8_t *const ref[4], int ref_stride,
                        uint32_t *res) {
  sad8x_4d(src, src_stride, ref, ref_stride, res, 8);
}

void vpx_sad8x16x4d_neon(const uint8_t *src, int src_stride,
                         const uint8_t *const ref[4], int ref_stride,
                         uint32_t *res) {
  sad8x_4d(src, src_stride, ref, ref_stride, res, 16);
}

static INLINE unsigned int horizontal_long_add_16x8(const uint16x8_t vec_lo,
                                                    const uint16x8_t vec_hi) {
  const uint32x4_t vec_l_lo =
      vaddl_u16(vget_low_u16(vec_lo), vget_high_u16(vec_lo));
  const uint32x4_t vec_l_hi =
      vaddl_u16(vget_low_u16(vec_hi), vget_high_u16(vec_hi));
  const uint32x4_t a = vaddq_u32(vec_l_lo, vec_l_hi);
  const uint64x2_t b = vpaddlq_u32(a);
  const uint32x2_t c = vadd_u32(vreinterpret_u32_u64(vget_low_u64(b)),
                                vreinterpret_u32_u64(vget_high_u64(b)));
  return vget_lane_u32(c, 0);
}

// Calculate the absolute difference of 64 bytes from vec_src_00, vec_src_16,
// vec_src_32, vec_src_48 and ref. Accumulate partial sums in vec_sum_ref_lo
// and vec_sum_ref_hi.
static void sad_neon_64(const uint8x16_t vec_src_00,
                        const uint8x16_t vec_src_16,
                        const uint8x16_t vec_src_32,
                        const uint8x16_t vec_src_48, const uint8_t *ref,
                        uint16x8_t *vec_sum_ref_lo,
                        uint16x8_t *vec_sum_ref_hi) {
  const uint8x16_t vec_ref_00 = vld1q_u8(ref);
  const uint8x16_t vec_ref_16 = vld1q_u8(ref + 16);
  const uint8x16_t vec_ref_32 = vld1q_u8(ref + 32);
  const uint8x16_t vec_ref_48 = vld1q_u8(ref + 48);

  *vec_sum_ref_lo = vabal_u8(*vec_sum_ref_lo, vget_low_u8(vec_src_00),
                             vget_low_u8(vec_ref_00));
  *vec_sum_ref_hi = vabal_u8(*vec_sum_ref_hi, vget_high_u8(vec_src_00),
                             vget_high_u8(vec_ref_00));
  *vec_sum_ref_lo = vabal_u8(*vec_sum_ref_lo, vget_low_u8(vec_src_16),
                             vget_low_u8(vec_ref_16));
  *vec_sum_ref_hi = vabal_u8(*vec_sum_ref_hi, vget_high_u8(vec_src_16),
                             vget_high_u8(vec_ref_16));
  *vec_sum_ref_lo = vabal_u8(*vec_sum_ref_lo, vget_low_u8(vec_src_32),
                             vget_low_u8(vec_ref_32));
  *vec_sum_ref_hi = vabal_u8(*vec_sum_ref_hi, vget_high_u8(vec_src_32),
                             vget_high_u8(vec_ref_32));
  *vec_sum_ref_lo = vabal_u8(*vec_sum_ref_lo, vget_low_u8(vec_src_48),
                             vget_low_u8(vec_ref_48));
  *vec_sum_ref_hi = vabal_u8(*vec_sum_ref_hi, vget_high_u8(vec_src_48),
                             vget_high_u8(vec_ref_48));
}

// Calculate the absolute difference of 32 bytes from vec_src_00, vec_src_16,
// and ref. Accumulate partial sums in vec_sum_ref_lo and vec_sum_ref_hi.
static void sad_neon_32(const uint8x16_t vec_src_00,
                        const uint8x16_t vec_src_16, const uint8_t *ref,
                        uint16x8_t *vec_sum_ref_lo,
                        uint16x8_t *vec_sum_ref_hi) {
  const uint8x16_t vec_ref_00 = vld1q_u8(ref);
  const uint8x16_t vec_ref_16 = vld1q_u8(ref + 16);

  *vec_sum_ref_lo = vabal_u8(*vec_sum_ref_lo, vget_low_u8(vec_src_00),
                             vget_low_u8(vec_ref_00));
  *vec_sum_ref_hi = vabal_u8(*vec_sum_ref_hi, vget_high_u8(vec_src_00),
                             vget_high_u8(vec_ref_00));
  *vec_sum_ref_lo = vabal_u8(*vec_sum_ref_lo, vget_low_u8(vec_src_16),
                             vget_low_u8(vec_ref_16));
  *vec_sum_ref_hi = vabal_u8(*vec_sum_ref_hi, vget_high_u8(vec_src_16),
                             vget_high_u8(vec_ref_16));
}

void vpx_sad64x64x4d_neon(const uint8_t *src, int src_stride,
                          const uint8_t *const ref[4], int ref_stride,
                          uint32_t *res) {
  int i;
  uint16x8_t vec_sum_ref0_lo = vdupq_n_u16(0);
  uint16x8_t vec_sum_ref0_hi = vdupq_n_u16(0);
  uint16x8_t vec_sum_ref1_lo = vdupq_n_u16(0);
  uint16x8_t vec_sum_ref1_hi = vdupq_n_u16(0);
  uint16x8_t vec_sum_ref2_lo = vdupq_n_u16(0);
  uint16x8_t vec_sum_ref2_hi = vdupq_n_u16(0);
  uint16x8_t vec_sum_ref3_lo = vdupq_n_u16(0);
  uint16x8_t vec_sum_ref3_hi = vdupq_n_u16(0);
  const uint8_t *ref0, *ref1, *ref2, *ref3;
  ref0 = ref[0];
  ref1 = ref[1];
  ref2 = ref[2];
  ref3 = ref[3];

  for (i = 0; i < 64; ++i) {
    const uint8x16_t vec_src_00 = vld1q_u8(src);
    const uint8x16_t vec_src_16 = vld1q_u8(src + 16);
    const uint8x16_t vec_src_32 = vld1q_u8(src + 32);
    const uint8x16_t vec_src_48 = vld1q_u8(src + 48);

    sad_neon_64(vec_src_00, vec_src_16, vec_src_32, vec_src_48, ref0,
                &vec_sum_ref0_lo, &vec_sum_ref0_hi);
    sad_neon_64(vec_src_00, vec_src_16, vec_src_32, vec_src_48, ref1,
                &vec_sum_ref1_lo, &vec_sum_ref1_hi);
    sad_neon_64(vec_src_00, vec_src_16, vec_src_32, vec_src_48, ref2,
                &vec_sum_ref2_lo, &vec_sum_ref2_hi);
    sad_neon_64(vec_src_00, vec_src_16, vec_src_32, vec_src_48, ref3,
                &vec_sum_ref3_lo, &vec_sum_ref3_hi);

    src += src_stride;
    ref0 += ref_stride;
    ref1 += ref_stride;
    ref2 += ref_stride;
    ref3 += ref_stride;
  }

  res[0] = horizontal_long_add_16x8(vec_sum_ref0_lo, vec_sum_ref0_hi);
  res[1] = horizontal_long_add_16x8(vec_sum_ref1_lo, vec_sum_ref1_hi);
  res[2] = horizontal_long_add_16x8(vec_sum_ref2_lo, vec_sum_ref2_hi);
  res[3] = horizontal_long_add_16x8(vec_sum_ref3_lo, vec_sum_ref3_hi);
}

void vpx_sad32x32x4d_neon(const uint8_t *src, int src_stride,
                          const uint8_t *const ref[4], int ref_stride,
                          uint32_t *res) {
  int i;
  uint16x8_t vec_sum_ref0_lo = vdupq_n_u16(0);
  uint16x8_t vec_sum_ref0_hi = vdupq_n_u16(0);
  uint16x8_t vec_sum_ref1_lo = vdupq_n_u16(0);
  uint16x8_t vec_sum_ref1_hi = vdupq_n_u16(0);
  uint16x8_t vec_sum_ref2_lo = vdupq_n_u16(0);
  uint16x8_t vec_sum_ref2_hi = vdupq_n_u16(0);
  uint16x8_t vec_sum_ref3_lo = vdupq_n_u16(0);
  uint16x8_t vec_sum_ref3_hi = vdupq_n_u16(0);
  const uint8_t *ref0, *ref1, *ref2, *ref3;
  ref0 = ref[0];
  ref1 = ref[1];
  ref2 = ref[2];
  ref3 = ref[3];

  for (i = 0; i < 32; ++i) {
    const uint8x16_t vec_src_00 = vld1q_u8(src);
    const uint8x16_t vec_src_16 = vld1q_u8(src + 16);

    sad_neon_32(vec_src_00, vec_src_16, ref0, &vec_sum_ref0_lo,
                &vec_sum_ref0_hi);
    sad_neon_32(vec_src_00, vec_src_16, ref1, &vec_sum_ref1_lo,
                &vec_sum_ref1_hi);
    sad_neon_32(vec_src_00, vec_src_16, ref2, &vec_sum_ref2_lo,
                &vec_sum_ref2_hi);
    sad_neon_32(vec_src_00, vec_src_16, ref3, &vec_sum_ref3_lo,
                &vec_sum_ref3_hi);

    src += src_stride;
    ref0 += ref_stride;
    ref1 += ref_stride;
    ref2 += ref_stride;
    ref3 += ref_stride;
  }

  res[0] = horizontal_long_add_16x8(vec_sum_ref0_lo, vec_sum_ref0_hi);
  res[1] = horizontal_long_add_16x8(vec_sum_ref1_lo, vec_sum_ref1_hi);
  res[2] = horizontal_long_add_16x8(vec_sum_ref2_lo, vec_sum_ref2_hi);
  res[3] = horizontal_long_add_16x8(vec_sum_ref3_lo, vec_sum_ref3_hi);
}

void vpx_sad16x16x4d_neon(const uint8_t *src, int src_stride,
                          const uint8_t *const ref[4], int ref_stride,
                          uint32_t *res) {
  int i;
  uint16x8_t vec_sum_ref0_lo = vdupq_n_u16(0);
  uint16x8_t vec_sum_ref0_hi = vdupq_n_u16(0);
  uint16x8_t vec_sum_ref1_lo = vdupq_n_u16(0);
  uint16x8_t vec_sum_ref1_hi = vdupq_n_u16(0);
  uint16x8_t vec_sum_ref2_lo = vdupq_n_u16(0);
  uint16x8_t vec_sum_ref2_hi = vdupq_n_u16(0);
  uint16x8_t vec_sum_ref3_lo = vdupq_n_u16(0);
  uint16x8_t vec_sum_ref3_hi = vdupq_n_u16(0);
  const uint8_t *ref0, *ref1, *ref2, *ref3;
  ref0 = ref[0];
  ref1 = ref[1];
  ref2 = ref[2];
  ref3 = ref[3];

  for (i = 0; i < 16; ++i) {
    const uint8x16_t vec_src = vld1q_u8(src);
    const uint8x16_t vec_ref0 = vld1q_u8(ref0);
    const uint8x16_t vec_ref1 = vld1q_u8(ref1);
    const uint8x16_t vec_ref2 = vld1q_u8(ref2);
    const uint8x16_t vec_ref3 = vld1q_u8(ref3);

    vec_sum_ref0_lo =
        vabal_u8(vec_sum_ref0_lo, vget_low_u8(vec_src), vget_low_u8(vec_ref0));
    vec_sum_ref0_hi = vabal_u8(vec_sum_ref0_hi, vget_high_u8(vec_src),
                               vget_high_u8(vec_ref0));
    vec_sum_ref1_lo =
        vabal_u8(vec_sum_ref1_lo, vget_low_u8(vec_src), vget_low_u8(vec_ref1));
    vec_sum_ref1_hi = vabal_u8(vec_sum_ref1_hi, vget_high_u8(vec_src),
                               vget_high_u8(vec_ref1));
    vec_sum_ref2_lo =
        vabal_u8(vec_sum_ref2_lo, vget_low_u8(vec_src), vget_low_u8(vec_ref2));
    vec_sum_ref2_hi = vabal_u8(vec_sum_ref2_hi, vget_high_u8(vec_src),
                               vget_high_u8(vec_ref2));
    vec_sum_ref3_lo =
        vabal_u8(vec_sum_ref3_lo, vget_low_u8(vec_src), vget_low_u8(vec_ref3));
    vec_sum_ref3_hi = vabal_u8(vec_sum_ref3_hi, vget_high_u8(vec_src),
                               vget_high_u8(vec_ref3));

    src += src_stride;
    ref0 += ref_stride;
    ref1 += ref_stride;
    ref2 += ref_stride;
    ref3 += ref_stride;
  }

  res[0] = horizontal_long_add_16x8(vec_sum_ref0_lo, vec_sum_ref0_hi);
  res[1] = horizontal_long_add_16x8(vec_sum_ref1_lo, vec_sum_ref1_hi);
  res[2] = horizontal_long_add_16x8(vec_sum_ref2_lo, vec_sum_ref2_hi);
  res[3] = horizontal_long_add_16x8(vec_sum_ref3_lo, vec_sum_ref3_hi);
}
