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

uint32_t vpx_sad8x4_neon(const uint8_t *src, int src_stride, const uint8_t *ref,
                         int ref_stride) {
  const uint16x8_t abs = sad8x(src, src_stride, ref, ref_stride, 4);
  return horizontal_add_16x8(abs);
}

uint32_t vpx_sad8x8_neon(const uint8_t *src, int src_stride, const uint8_t *ref,
                         int ref_stride) {
  const uint16x8_t abs = sad8x(src, src_stride, ref, ref_stride, 8);
  return horizontal_add_16x8(abs);
}

uint32_t vpx_sad8x16_neon(const uint8_t *src, int src_stride,
                          const uint8_t *ref, int ref_stride) {
  const uint16x8_t abs = sad8x(src, src_stride, ref, ref_stride, 16);
  return horizontal_add_16x8(abs);
}

unsigned int vpx_sad16x8_neon(unsigned char *src_ptr, int src_stride,
                              unsigned char *ref_ptr, int ref_stride) {
  uint8x16_t q0, q4;
  uint16x8_t q12, q13;
  uint32x4_t q1;
  uint64x2_t q3;
  uint32x2_t d5;
  int i;

  q0 = vld1q_u8(src_ptr);
  src_ptr += src_stride;
  q4 = vld1q_u8(ref_ptr);
  ref_ptr += ref_stride;
  q12 = vabdl_u8(vget_low_u8(q0), vget_low_u8(q4));
  q13 = vabdl_u8(vget_high_u8(q0), vget_high_u8(q4));

  for (i = 0; i < 7; i++) {
    q0 = vld1q_u8(src_ptr);
    src_ptr += src_stride;
    q4 = vld1q_u8(ref_ptr);
    ref_ptr += ref_stride;
    q12 = vabal_u8(q12, vget_low_u8(q0), vget_low_u8(q4));
    q13 = vabal_u8(q13, vget_high_u8(q0), vget_high_u8(q4));
  }

  q12 = vaddq_u16(q12, q13);
  q1 = vpaddlq_u16(q12);
  q3 = vpaddlq_u32(q1);
  d5 = vadd_u32(vreinterpret_u32_u64(vget_low_u64(q3)),
                vreinterpret_u32_u64(vget_high_u64(q3)));

  return vget_lane_u32(d5, 0);
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

unsigned int vpx_sad64x64_neon(const uint8_t *src, int src_stride,
                               const uint8_t *ref, int ref_stride) {
  int i;
  uint16x8_t vec_accum_lo = vdupq_n_u16(0);
  uint16x8_t vec_accum_hi = vdupq_n_u16(0);
  for (i = 0; i < 64; ++i) {
    const uint8x16_t vec_src_00 = vld1q_u8(src);
    const uint8x16_t vec_src_16 = vld1q_u8(src + 16);
    const uint8x16_t vec_src_32 = vld1q_u8(src + 32);
    const uint8x16_t vec_src_48 = vld1q_u8(src + 48);
    const uint8x16_t vec_ref_00 = vld1q_u8(ref);
    const uint8x16_t vec_ref_16 = vld1q_u8(ref + 16);
    const uint8x16_t vec_ref_32 = vld1q_u8(ref + 32);
    const uint8x16_t vec_ref_48 = vld1q_u8(ref + 48);
    src += src_stride;
    ref += ref_stride;
    vec_accum_lo = vabal_u8(vec_accum_lo, vget_low_u8(vec_src_00),
                            vget_low_u8(vec_ref_00));
    vec_accum_hi = vabal_u8(vec_accum_hi, vget_high_u8(vec_src_00),
                            vget_high_u8(vec_ref_00));
    vec_accum_lo = vabal_u8(vec_accum_lo, vget_low_u8(vec_src_16),
                            vget_low_u8(vec_ref_16));
    vec_accum_hi = vabal_u8(vec_accum_hi, vget_high_u8(vec_src_16),
                            vget_high_u8(vec_ref_16));
    vec_accum_lo = vabal_u8(vec_accum_lo, vget_low_u8(vec_src_32),
                            vget_low_u8(vec_ref_32));
    vec_accum_hi = vabal_u8(vec_accum_hi, vget_high_u8(vec_src_32),
                            vget_high_u8(vec_ref_32));
    vec_accum_lo = vabal_u8(vec_accum_lo, vget_low_u8(vec_src_48),
                            vget_low_u8(vec_ref_48));
    vec_accum_hi = vabal_u8(vec_accum_hi, vget_high_u8(vec_src_48),
                            vget_high_u8(vec_ref_48));
  }
  return horizontal_long_add_16x8(vec_accum_lo, vec_accum_hi);
}

unsigned int vpx_sad32x32_neon(const uint8_t *src, int src_stride,
                               const uint8_t *ref, int ref_stride) {
  int i;
  uint16x8_t vec_accum_lo = vdupq_n_u16(0);
  uint16x8_t vec_accum_hi = vdupq_n_u16(0);

  for (i = 0; i < 32; ++i) {
    const uint8x16_t vec_src_00 = vld1q_u8(src);
    const uint8x16_t vec_src_16 = vld1q_u8(src + 16);
    const uint8x16_t vec_ref_00 = vld1q_u8(ref);
    const uint8x16_t vec_ref_16 = vld1q_u8(ref + 16);
    src += src_stride;
    ref += ref_stride;
    vec_accum_lo = vabal_u8(vec_accum_lo, vget_low_u8(vec_src_00),
                            vget_low_u8(vec_ref_00));
    vec_accum_hi = vabal_u8(vec_accum_hi, vget_high_u8(vec_src_00),
                            vget_high_u8(vec_ref_00));
    vec_accum_lo = vabal_u8(vec_accum_lo, vget_low_u8(vec_src_16),
                            vget_low_u8(vec_ref_16));
    vec_accum_hi = vabal_u8(vec_accum_hi, vget_high_u8(vec_src_16),
                            vget_high_u8(vec_ref_16));
  }
  return horizontal_add_16x8(vaddq_u16(vec_accum_lo, vec_accum_hi));
}

unsigned int vpx_sad16x16_neon(const uint8_t *src, int src_stride,
                               const uint8_t *ref, int ref_stride) {
  int i;
  uint16x8_t vec_accum_lo = vdupq_n_u16(0);
  uint16x8_t vec_accum_hi = vdupq_n_u16(0);

  for (i = 0; i < 16; ++i) {
    const uint8x16_t vec_src = vld1q_u8(src);
    const uint8x16_t vec_ref = vld1q_u8(ref);
    src += src_stride;
    ref += ref_stride;
    vec_accum_lo =
        vabal_u8(vec_accum_lo, vget_low_u8(vec_src), vget_low_u8(vec_ref));
    vec_accum_hi =
        vabal_u8(vec_accum_hi, vget_high_u8(vec_src), vget_high_u8(vec_ref));
  }
  return horizontal_add_16x8(vaddq_u16(vec_accum_lo, vec_accum_hi));
}
