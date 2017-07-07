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

#include "vpx/vpx_integer.h"
#include "vpx_dsp/arm/mem_neon.h"

// TODO(johannkoenig): combine with avg_neon.h version.
static INLINE uint32_t horizontal_add_16x8(const uint16x8_t vec_16x8) {
  const uint32x4_t a = vpaddlq_u16(vec_16x8);
  const uint64x2_t b = vpaddlq_u32(a);
  const uint32x2_t c = vadd_u32(vreinterpret_u32_u64(vget_low_u64(b)),
                                vreinterpret_u32_u64(vget_high_u64(b)));
  return vget_lane_u32(c, 0);
}

uint32_t vpx_sad4x4_neon(const uint8_t *src_ptr, int src_stride,
                         const uint8_t *ref_ptr, int ref_stride) {
  const uint8x16_t src_u8 = load_unaligned_u8q(src_ptr, src_stride);
  const uint8x16_t ref_u8 = load_unaligned_u8q(ref_ptr, ref_stride);
  uint16x8_t abs = vabdl_u8(vget_low_u8(src_u8), vget_low_u8(ref_u8));
  abs = vabal_u8(abs, vget_high_u8(src_u8), vget_high_u8(ref_u8));
  return horizontal_add_16x8(abs);
}

uint32_t vpx_sad4x4_avg_neon(const uint8_t *src_ptr, int src_stride,
                             const uint8_t *ref_ptr, int ref_stride,
                             const uint8_t *second_pred) {
  const uint8x16_t src_u8 = load_unaligned_u8q(src_ptr, src_stride);
  const uint8x16_t ref_u8 = load_unaligned_u8q(ref_ptr, ref_stride);
  const uint8x16_t second_pred_u8 = vld1q_u8(second_pred);
  const uint8x16_t avg = vrhaddq_u8(ref_u8, second_pred_u8);
  uint16x8_t abs = vabdl_u8(vget_low_u8(src_u8), vget_low_u8(avg));
  abs = vabal_u8(abs, vget_high_u8(src_u8), vget_high_u8(avg));
  return horizontal_add_16x8(abs);
}

uint32_t vpx_sad4x8_neon(const uint8_t *src_ptr, int src_stride,
                         const uint8_t *ref_ptr, int ref_stride) {
  int i;
  uint16x8_t abs = vdupq_n_u16(0);
  for (i = 0; i < 8; i += 4) {
    const uint8x16_t src_u8 = load_unaligned_u8q(src_ptr, src_stride);
    const uint8x16_t ref_u8 = load_unaligned_u8q(ref_ptr, ref_stride);
    src_ptr += 4 * src_stride;
    ref_ptr += 4 * ref_stride;
    abs = vabal_u8(abs, vget_low_u8(src_u8), vget_low_u8(ref_u8));
    abs = vabal_u8(abs, vget_high_u8(src_u8), vget_high_u8(ref_u8));
  }

  return horizontal_add_16x8(abs);
}

uint32_t vpx_sad4x8_avg_neon(const uint8_t *src_ptr, int src_stride,
                             const uint8_t *ref_ptr, int ref_stride,
                             const uint8_t *second_pred) {
  int i;
  uint16x8_t abs = vdupq_n_u16(0);
  for (i = 0; i < 8; i += 4) {
    const uint8x16_t src_u8 = load_unaligned_u8q(src_ptr, src_stride);
    const uint8x16_t ref_u8 = load_unaligned_u8q(ref_ptr, ref_stride);
    const uint8x16_t second_pred_u8 = vld1q_u8(second_pred);
    const uint8x16_t avg = vrhaddq_u8(ref_u8, second_pred_u8);
    src_ptr += 4 * src_stride;
    ref_ptr += 4 * ref_stride;
    second_pred += 16;
    abs = vabal_u8(abs, vget_low_u8(src_u8), vget_low_u8(avg));
    abs = vabal_u8(abs, vget_high_u8(src_u8), vget_high_u8(avg));
  }

  return horizontal_add_16x8(abs);
}

static INLINE uint16x8_t sad8x(const uint8_t *a, int a_stride, const uint8_t *b,
                               int b_stride, const int height) {
  int i;
  uint16x8_t abs = vdupq_n_u16(0);

  for (i = 0; i < height; ++i) {
    const uint8x8_t a_u8 = vld1_u8(a);
    const uint8x8_t b_u8 = vld1_u8(b);
    a += a_stride;
    b += b_stride;
    abs = vabal_u8(abs, a_u8, b_u8);
  }
  return abs;
}

static INLINE uint16x8_t sad8x_avg(const uint8_t *a, int a_stride,
                                   const uint8_t *b, int b_stride,
                                   const uint8_t *c, const int height) {
  int i;
  uint16x8_t abs = vdupq_n_u16(0);

  for (i = 0; i < height; ++i) {
    const uint8x8_t a_u8 = vld1_u8(a);
    const uint8x8_t b_u8 = vld1_u8(b);
    const uint8x8_t c_u8 = vld1_u8(c);
    const uint8x8_t avg = vrhadd_u8(b_u8, c_u8);
    a += a_stride;
    b += b_stride;
    c += 8;
    abs = vabal_u8(abs, a_u8, avg);
  }
  return abs;
}

#define sad8xN(n)                                                      \
  uint32_t vpx_sad8x##n##_neon(const uint8_t *src, int src_stride,     \
                               const uint8_t *ref, int ref_stride) {   \
    const uint16x8_t abs = sad8x(src, src_stride, ref, ref_stride, n); \
    return horizontal_add_16x8(abs);                                   \
  }                                                                    \
                                                                       \
  uint32_t vpx_sad8x##n##_avg_neon(const uint8_t *src, int src_stride, \
                                   const uint8_t *ref, int ref_stride, \
                                   const uint8_t *second_pred) {       \
    const uint16x8_t abs =                                             \
        sad8x_avg(src, src_stride, ref, ref_stride, second_pred, n);   \
    return horizontal_add_16x8(abs);                                   \
  }

sad8xN(4);
sad8xN(8);
sad8xN(16);

static INLINE uint16x8_t sad16x(const uint8_t *a, int a_stride,
                                const uint8_t *b, int b_stride,
                                const int height) {
  int i;
  uint16x8_t abs = vdupq_n_u16(0);

  for (i = 0; i < height; ++i) {
    const uint8x16_t a_u8 = vld1q_u8(a);
    const uint8x16_t b_u8 = vld1q_u8(b);
    a += a_stride;
    b += b_stride;
    abs = vabal_u8(abs, vget_low_u8(a_u8), vget_low_u8(b_u8));
    abs = vabal_u8(abs, vget_high_u8(a_u8), vget_high_u8(b_u8));
  }
  return abs;
}

uint32_t vpx_sad16x8_neon(const uint8_t *src, int src_stride,
                          const uint8_t *ref, int ref_stride) {
  const uint16x8_t abs = sad16x(src, src_stride, ref, ref_stride, 8);
  return horizontal_add_16x8(abs);
}

uint32_t vpx_sad16x16_neon(const uint8_t *src, int src_stride,
                           const uint8_t *ref, int ref_stride) {
  const uint16x8_t abs = sad16x(src, src_stride, ref, ref_stride, 16);
  return horizontal_add_16x8(abs);
}

uint32_t vpx_sad16x32_neon(const uint8_t *src, int src_stride,
                           const uint8_t *ref, int ref_stride) {
  const uint16x8_t abs = sad16x(src, src_stride, ref, ref_stride, 32);
  return horizontal_add_16x8(abs);
}

static INLINE uint16x8_t sad32x(const uint8_t *a, int a_stride,
                                const uint8_t *b, int b_stride,
                                const int height) {
  int i;
  uint16x8_t abs = vdupq_n_u16(0);

  for (i = 0; i < height; ++i) {
    const uint8x16_t a_lo = vld1q_u8(a);
    const uint8x16_t a_hi = vld1q_u8(a + 16);
    const uint8x16_t b_lo = vld1q_u8(b);
    const uint8x16_t b_hi = vld1q_u8(b + 16);
    a += a_stride;
    b += b_stride;
    abs = vabal_u8(abs, vget_low_u8(a_lo), vget_low_u8(b_lo));
    abs = vabal_u8(abs, vget_high_u8(a_lo), vget_high_u8(b_lo));
    abs = vabal_u8(abs, vget_low_u8(a_hi), vget_low_u8(b_hi));
    abs = vabal_u8(abs, vget_high_u8(a_hi), vget_high_u8(b_hi));
  }
  return abs;
}

uint32_t vpx_sad32x16_neon(const uint8_t *src, int src_stride,
                           const uint8_t *ref, int ref_stride) {
  const uint16x8_t abs = sad32x(src, src_stride, ref, ref_stride, 16);
  return horizontal_add_16x8(abs);
}

uint32_t vpx_sad32x32_neon(const uint8_t *src, int src_stride,
                           const uint8_t *ref, int ref_stride) {
  const uint16x8_t abs = sad32x(src, src_stride, ref, ref_stride, 32);
  return horizontal_add_16x8(abs);
}

uint32_t vpx_sad32x64_neon(const uint8_t *src, int src_stride,
                           const uint8_t *ref, int ref_stride) {
  const uint16x8_t abs = sad32x(src, src_stride, ref, ref_stride, 64);
  return horizontal_add_16x8(abs);
}

static INLINE uint32_t horizontal_add_32x4(const uint32x4_t a) {
  const uint64x2_t b = vpaddlq_u32(a);
  const uint32x2_t c = vadd_u32(vreinterpret_u32_u64(vget_low_u64(b)),
                                vreinterpret_u32_u64(vget_high_u64(b)));
  return vget_lane_u32(c, 0);
}

static INLINE uint32x4_t sad64x(const uint8_t *a, int a_stride,
                                const uint8_t *b, int b_stride,
                                const int height) {
  int i;
  uint16x8_t abs_0 = vdupq_n_u16(0);
  uint16x8_t abs_1 = vdupq_n_u16(0);

  for (i = 0; i < height; ++i) {
    const uint8x16_t a_0 = vld1q_u8(a);
    const uint8x16_t a_1 = vld1q_u8(a + 16);
    const uint8x16_t a_2 = vld1q_u8(a + 32);
    const uint8x16_t a_3 = vld1q_u8(a + 48);
    const uint8x16_t b_0 = vld1q_u8(b);
    const uint8x16_t b_1 = vld1q_u8(b + 16);
    const uint8x16_t b_2 = vld1q_u8(b + 32);
    const uint8x16_t b_3 = vld1q_u8(b + 48);
    a += a_stride;
    b += b_stride;
    abs_0 = vabal_u8(abs_0, vget_low_u8(a_0), vget_low_u8(b_0));
    abs_0 = vabal_u8(abs_0, vget_high_u8(a_0), vget_high_u8(b_0));
    abs_0 = vabal_u8(abs_0, vget_low_u8(a_1), vget_low_u8(b_1));
    abs_0 = vabal_u8(abs_0, vget_high_u8(a_1), vget_high_u8(b_1));
    abs_1 = vabal_u8(abs_1, vget_low_u8(a_2), vget_low_u8(b_2));
    abs_1 = vabal_u8(abs_1, vget_high_u8(a_2), vget_high_u8(b_2));
    abs_1 = vabal_u8(abs_1, vget_low_u8(a_3), vget_low_u8(b_3));
    abs_1 = vabal_u8(abs_1, vget_high_u8(a_3), vget_high_u8(b_3));
  }

  {
    const uint32x4_t sum = vpaddlq_u16(abs_0);
    return vpadalq_u16(sum, abs_1);
  }
}

uint32_t vpx_sad64x32_neon(const uint8_t *src, int src_stride,
                           const uint8_t *ref, int ref_stride) {
  const uint32x4_t abs = sad64x(src, src_stride, ref, ref_stride, 32);
  return horizontal_add_32x4(abs);
}

uint32_t vpx_sad64x64_neon(const uint8_t *src, int src_stride,
                           const uint8_t *ref, int ref_stride) {
  const uint32x4_t abs = sad64x(src, src_stride, ref, ref_stride, 64);
  return horizontal_add_32x4(abs);
}
