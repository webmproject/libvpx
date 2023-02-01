/*
 *  Copyright (c) 2022 The WebM project authors. All Rights Reserved.
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

static INLINE uint32_t highbd_sad4xh_neon(const uint8_t *src_ptr,
                                          int src_stride,
                                          const uint8_t *ref_ptr,
                                          int ref_stride, int h) {
  const uint16_t *src16_ptr = CONVERT_TO_SHORTPTR(src_ptr);
  const uint16_t *ref16_ptr = CONVERT_TO_SHORTPTR(ref_ptr);
  uint32x4_t sum = vdupq_n_u32(0);

  int i = h;
  do {
    uint16x4_t s = vld1_u16(src16_ptr);
    uint16x4_t r = vld1_u16(ref16_ptr);
    sum = vabal_u16(sum, s, r);

    src16_ptr += src_stride;
    ref16_ptr += ref_stride;
  } while (--i != 0);

  return horizontal_add_uint32x4(sum);
}

static INLINE uint32_t highbd_sad8xh_neon(const uint8_t *src_ptr,
                                          int src_stride,
                                          const uint8_t *ref_ptr,
                                          int ref_stride, int h) {
  const uint16_t *src16_ptr = CONVERT_TO_SHORTPTR(src_ptr);
  const uint16_t *ref16_ptr = CONVERT_TO_SHORTPTR(ref_ptr);
  uint32x4_t sum = vdupq_n_u32(0);

  int i = h;
  do {
    uint16x8_t s = vld1q_u16(src16_ptr);
    uint16x8_t r = vld1q_u16(ref16_ptr);
    uint16x8_t diff = vabdq_u16(s, r);
    sum = vpadalq_u16(sum, diff);

    src16_ptr += src_stride;
    ref16_ptr += ref_stride;
  } while (--i != 0);

  return horizontal_add_uint32x4(sum);
}

static INLINE uint32_t highbd_sad16xh_neon(const uint8_t *src_ptr,
                                           int src_stride,
                                           const uint8_t *ref_ptr,
                                           int ref_stride, int h) {
  const uint16_t *src16_ptr = CONVERT_TO_SHORTPTR(src_ptr);
  const uint16_t *ref16_ptr = CONVERT_TO_SHORTPTR(ref_ptr);
  uint32x4_t sum[2] = { vdupq_n_u32(0), vdupq_n_u32(0) };

  int i = h;
  do {
    uint16x8_t s0, s1, r0, r1;
    uint16x8_t diff0, diff1;

    s0 = vld1q_u16(src16_ptr);
    r0 = vld1q_u16(ref16_ptr);
    diff0 = vabdq_u16(s0, r0);
    sum[0] = vpadalq_u16(sum[0], diff0);

    s1 = vld1q_u16(src16_ptr + 8);
    r1 = vld1q_u16(ref16_ptr + 8);
    diff1 = vabdq_u16(s1, r1);
    sum[1] = vpadalq_u16(sum[1], diff1);

    src16_ptr += src_stride;
    ref16_ptr += ref_stride;
  } while (--i != 0);

  sum[0] = vaddq_u32(sum[0], sum[1]);
  return horizontal_add_uint32x4(sum[0]);
}

static INLINE uint32_t highbd_sadwxh_neon(const uint8_t *src_ptr,
                                          int src_stride,
                                          const uint8_t *ref_ptr,
                                          int ref_stride, int w, int h) {
  const uint16_t *src16_ptr = CONVERT_TO_SHORTPTR(src_ptr);
  const uint16_t *ref16_ptr = CONVERT_TO_SHORTPTR(ref_ptr);
  uint32x4_t sum[4] = { vdupq_n_u32(0), vdupq_n_u32(0), vdupq_n_u32(0),
                        vdupq_n_u32(0) };

  int i = h;
  do {
    int j = 0;
    do {
      uint16x8_t s0, s1, s2, s3, r0, r1, r2, r3;
      uint16x8_t diff0, diff1, diff2, diff3;

      s0 = vld1q_u16(src16_ptr + j);
      r0 = vld1q_u16(ref16_ptr + j);
      diff0 = vabdq_u16(s0, r0);
      sum[0] = vpadalq_u16(sum[0], diff0);

      s1 = vld1q_u16(src16_ptr + j + 8);
      r1 = vld1q_u16(ref16_ptr + j + 8);
      diff1 = vabdq_u16(s1, r1);
      sum[1] = vpadalq_u16(sum[1], diff1);

      s2 = vld1q_u16(src16_ptr + j + 16);
      r2 = vld1q_u16(ref16_ptr + j + 16);
      diff2 = vabdq_u16(s2, r2);
      sum[2] = vpadalq_u16(sum[2], diff2);

      s3 = vld1q_u16(src16_ptr + j + 24);
      r3 = vld1q_u16(ref16_ptr + j + 24);
      diff3 = vabdq_u16(s3, r3);
      sum[3] = vpadalq_u16(sum[3], diff3);

      j += 32;
    } while (j < w);

    src16_ptr += src_stride;
    ref16_ptr += ref_stride;
  } while (--i != 0);

  sum[0] = vaddq_u32(sum[0], sum[1]);
  sum[2] = vaddq_u32(sum[2], sum[3]);
  sum[0] = vaddq_u32(sum[0], sum[2]);

  return horizontal_add_uint32x4(sum[0]);
}

static INLINE unsigned int highbd_sad64xh_neon(const uint8_t *src_ptr,
                                               int src_stride,
                                               const uint8_t *ref_ptr,
                                               int ref_stride, int h) {
  return highbd_sadwxh_neon(src_ptr, src_stride, ref_ptr, ref_stride, 64, h);
}

static INLINE unsigned int highbd_sad32xh_neon(const uint8_t *src_ptr,
                                               int src_stride,
                                               const uint8_t *ref_ptr,
                                               int ref_stride, int h) {
  return highbd_sadwxh_neon(src_ptr, src_stride, ref_ptr, ref_stride, 32, h);
}

#define HBD_SAD_WXH_NEON(w, h)                                            \
  unsigned int vpx_highbd_sad##w##x##h##_neon(                            \
      const uint8_t *src, int src_stride, const uint8_t *ref,             \
      int ref_stride) {                                                   \
    return highbd_sad##w##xh_neon(src, src_stride, ref, ref_stride, (h)); \
  }

HBD_SAD_WXH_NEON(4, 4)
HBD_SAD_WXH_NEON(4, 8)

HBD_SAD_WXH_NEON(8, 4)
HBD_SAD_WXH_NEON(8, 8)
HBD_SAD_WXH_NEON(8, 16)

HBD_SAD_WXH_NEON(16, 8)
HBD_SAD_WXH_NEON(16, 16)
HBD_SAD_WXH_NEON(16, 32)

HBD_SAD_WXH_NEON(32, 16)
HBD_SAD_WXH_NEON(32, 32)
HBD_SAD_WXH_NEON(32, 64)

HBD_SAD_WXH_NEON(64, 32)
HBD_SAD_WXH_NEON(64, 64)

static VPX_FORCE_INLINE uint32_t highbd_sad4_avg_neon(
    const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr,
    int ref_stride, const uint8_t *second_pred, int width, int height) {
  int i, j;
  uint32x4_t sum_abs_diff = vdupq_n_u32(0);
  const uint16_t *src16_ptr = CONVERT_TO_SHORTPTR(src_ptr);
  const uint16_t *ref16_ptr = CONVERT_TO_SHORTPTR(ref_ptr);
  const uint16_t *pred_ptr = CONVERT_TO_SHORTPTR(second_pred);
  for (i = 0; i < height; i++) {
    for (j = 0; j < width; j += 4) {
      const uint16x4_t a_u16 = vld1_u16(src16_ptr + j);
      const uint16x4_t b_u16 = vld1_u16(ref16_ptr + j);
      const uint16x4_t c_u16 = vld1_u16(pred_ptr + j);
      const uint16x4_t avg = vrhadd_u16(b_u16, c_u16);
      sum_abs_diff = vabal_u16(sum_abs_diff, a_u16, avg);
    }
    src16_ptr += src_stride;
    ref16_ptr += ref_stride;
    pred_ptr += width;
  }

  return horizontal_add_uint32x4(sum_abs_diff);
}

static VPX_FORCE_INLINE uint32_t highbd_sad8_avg_neon(
    const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr,
    int ref_stride, const uint8_t *second_pred, int width, int height) {
  int i, j;
  uint32x4_t sum_abs_diff = vdupq_n_u32(0);
  const uint16_t *src16_ptr = CONVERT_TO_SHORTPTR(src_ptr);
  const uint16_t *ref16_ptr = CONVERT_TO_SHORTPTR(ref_ptr);
  const uint16_t *pred_ptr = CONVERT_TO_SHORTPTR(second_pred);
  for (i = 0; i < height; i++) {
    for (j = 0; j < width; j += 8) {
      const uint16x8_t a_u16 = vld1q_u16(src16_ptr + j);
      const uint16x8_t b_u16 = vld1q_u16(ref16_ptr + j);
      const uint16x8_t c_u16 = vld1q_u16(pred_ptr + j);
      const uint16x8_t avg = vrhaddq_u16(b_u16, c_u16);
      sum_abs_diff =
          vabal_u16(sum_abs_diff, vget_low_u16(a_u16), vget_low_u16(avg));
      sum_abs_diff =
          vabal_u16(sum_abs_diff, vget_high_u16(a_u16), vget_high_u16(avg));
    }
    src16_ptr += src_stride;
    ref16_ptr += ref_stride;
    pred_ptr += width;
  }

  return horizontal_add_uint32x4(sum_abs_diff);
}

#define highbd_sad4MxN_avg(m, n)                                          \
  unsigned int vpx_highbd_sad##m##x##n##_avg_neon(                        \
      const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr,     \
      int ref_stride, const uint8_t *second_pred) {                       \
    return highbd_sad4_avg_neon(src_ptr, src_stride, ref_ptr, ref_stride, \
                                second_pred, m, n);                       \
  }

#define highbd_sadMxN_avg(m, n)                                           \
  unsigned int vpx_highbd_sad##m##x##n##_avg_neon(                        \
      const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr,     \
      int ref_stride, const uint8_t *second_pred) {                       \
    return highbd_sad8_avg_neon(src_ptr, src_stride, ref_ptr, ref_stride, \
                                second_pred, m, n);                       \
  }

#define highbd_sadMxNx4D(m, n)                                                 \
  void vpx_highbd_sad##m##x##n##x4d_neon(                                      \
      const uint8_t *src_ptr, int src_stride,                                  \
      const uint8_t *const ref_array[4], int ref_stride,                       \
      uint32_t sad_array[4]) {                                                 \
    int i;                                                                     \
    for (i = 0; i < 4; ++i) {                                                  \
      sad_array[i] = vpx_highbd_sad##m##x##n##_neon(src_ptr, src_stride,       \
                                                    ref_array[i], ref_stride); \
    }                                                                          \
  }

/* clang-format off */
// 4x4
highbd_sad4MxN_avg(4, 4)
highbd_sadMxNx4D(4, 4)

// 4x8
highbd_sad4MxN_avg(4, 8)
highbd_sadMxNx4D(4, 8)

// 8x4
highbd_sadMxN_avg(8, 4)
highbd_sadMxNx4D(8, 4)

// 8x8
highbd_sadMxN_avg(8, 8)
highbd_sadMxNx4D(8, 8)

// 8x16
highbd_sadMxN_avg(8, 16)
highbd_sadMxNx4D(8, 16)

// 16x8
highbd_sadMxN_avg(16, 8)
highbd_sadMxNx4D(16, 8)

// 16x16
highbd_sadMxN_avg(16, 16)
highbd_sadMxNx4D(16, 16)

// 16x32
highbd_sadMxN_avg(16, 32)
highbd_sadMxNx4D(16, 32)

// 32x16
highbd_sadMxN_avg(32, 16)
highbd_sadMxNx4D(32, 16)

// 32x32
highbd_sadMxN_avg(32, 32)
highbd_sadMxNx4D(32, 32)

// 32x64
highbd_sadMxN_avg(32, 64)
highbd_sadMxNx4D(32, 64)

// 64x32
highbd_sadMxN_avg(64, 32)
highbd_sadMxNx4D(64, 32)

// 64x64
highbd_sadMxN_avg(64, 64)
highbd_sadMxNx4D(64, 64)
    /* clang-format on */
