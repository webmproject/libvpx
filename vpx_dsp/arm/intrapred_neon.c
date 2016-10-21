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
#include "vpx/vpx_integer.h"

//------------------------------------------------------------------------------
// DC 4x4

static INLINE uint16x4_t dc_sum_4(const uint8_t *ref) {
  const uint8x8_t ref_u8 = vld1_u8(ref);
  const uint16x4_t p0 = vpaddl_u8(ref_u8);
  return vpadd_u16(p0, p0);
}

static INLINE void dc_store_4x4(uint8_t *dst, ptrdiff_t stride,
                                const uint8x8_t dc) {
  const uint8x8_t dc_dup = vdup_lane_u8(dc, 0);
  int i;
  for (i = 0; i < 4; ++i, dst += stride) {
    vst1_lane_u32((uint32_t *)dst, vreinterpret_u32_u8(dc_dup), 0);
  }
}

void vpx_dc_predictor_4x4_neon(uint8_t *dst, ptrdiff_t stride,
                               const uint8_t *above, const uint8_t *left) {
  const uint8x8_t a = vld1_u8(above);
  const uint8x8_t l = vld1_u8(left);
  const uint16x8_t al = vaddl_u8(a, l);
  uint16x4_t sum;
  uint8x8_t dc;
  sum = vpadd_u16(vget_low_u16(al), vget_low_u16(al));
  sum = vpadd_u16(sum, sum);
  dc = vreinterpret_u8_u16(vrshr_n_u16(sum, 3));
  dc_store_4x4(dst, stride, dc);
}

void vpx_dc_left_predictor_4x4_neon(uint8_t *dst, ptrdiff_t stride,
                                    const uint8_t *above, const uint8_t *left) {
  const uint16x4_t sum = dc_sum_4(left);
  const uint8x8_t dc = vreinterpret_u8_u16(vrshr_n_u16(sum, 2));
  (void)above;
  dc_store_4x4(dst, stride, dc);
}

void vpx_dc_top_predictor_4x4_neon(uint8_t *dst, ptrdiff_t stride,
                                   const uint8_t *above, const uint8_t *left) {
  const uint16x4_t sum = dc_sum_4(above);
  const uint8x8_t dc = vreinterpret_u8_u16(vrshr_n_u16(sum, 2));
  (void)left;
  dc_store_4x4(dst, stride, dc);
}

void vpx_dc_128_predictor_4x4_neon(uint8_t *dst, ptrdiff_t stride,
                                   const uint8_t *above, const uint8_t *left) {
  const uint8x8_t dc = vdup_n_u8(0x80);
  (void)above;
  (void)left;
  dc_store_4x4(dst, stride, dc);
}

//------------------------------------------------------------------------------
// DC 8x8

static INLINE uint16x4_t dc_sum_8(const uint8_t *ref) {
  const uint8x8_t ref_u8 = vld1_u8(ref);
  uint16x4_t sum = vpaddl_u8(ref_u8);
  sum = vpadd_u16(sum, sum);
  return vpadd_u16(sum, sum);
}

static INLINE void dc_store_8x8(uint8_t *dst, ptrdiff_t stride,
                                const uint8x8_t dc) {
  const uint8x8_t dc_dup = vdup_lane_u8(dc, 0);
  int i;
  for (i = 0; i < 8; ++i, dst += stride) {
    vst1_u8(dst, dc_dup);
  }
}

void vpx_dc_predictor_8x8_neon(uint8_t *dst, ptrdiff_t stride,
                               const uint8_t *above, const uint8_t *left) {
  const uint8x8_t above_u8 = vld1_u8(above);
  const uint8x8_t left_u8 = vld1_u8(left);
  const uint8x16_t above_and_left = vcombine_u8(above_u8, left_u8);
  const uint16x8_t p0 = vpaddlq_u8(above_and_left);
  uint16x4_t sum = vadd_u16(vget_low_u16(p0), vget_high_u16(p0));
  uint8x8_t dc;
  sum = vpadd_u16(sum, sum);
  sum = vpadd_u16(sum, sum);
  dc = vreinterpret_u8_u16(vrshr_n_u16(sum, 4));
  dc_store_8x8(dst, stride, dc);
}

void vpx_dc_left_predictor_8x8_neon(uint8_t *dst, ptrdiff_t stride,
                                    const uint8_t *above, const uint8_t *left) {
  const uint16x4_t sum = dc_sum_8(left);
  const uint8x8_t dc = vreinterpret_u8_u16(vrshr_n_u16(sum, 3));
  (void)above;
  dc_store_8x8(dst, stride, dc);
}

void vpx_dc_top_predictor_8x8_neon(uint8_t *dst, ptrdiff_t stride,
                                   const uint8_t *above, const uint8_t *left) {
  const uint16x4_t sum = dc_sum_8(above);
  const uint8x8_t dc = vreinterpret_u8_u16(vrshr_n_u16(sum, 3));
  (void)left;
  dc_store_8x8(dst, stride, dc);
}

void vpx_dc_128_predictor_8x8_neon(uint8_t *dst, ptrdiff_t stride,
                                   const uint8_t *above, const uint8_t *left) {
  const uint8x8_t dc = vdup_n_u8(0x80);
  (void)above;
  (void)left;
  dc_store_8x8(dst, stride, dc);
}

//------------------------------------------------------------------------------
// DC 16x16

static INLINE uint16x4_t dc_sum_16(const uint8_t *ref) {
  const uint8x16_t ref_u8 = vld1q_u8(ref);
  const uint16x8_t p0 = vpaddlq_u8(ref_u8);
  uint16x4_t sum = vadd_u16(vget_low_u16(p0), vget_high_u16(p0));
  sum = vpadd_u16(sum, sum);
  return vpadd_u16(sum, sum);
}

static INLINE void dc_store_16x16(uint8_t *dst, ptrdiff_t stride,
                                  const uint8x8_t dc) {
  const uint8x16_t dc_dup = vdupq_lane_u8(dc, 0);
  int i;
  for (i = 0; i < 16; ++i, dst += stride) {
    vst1q_u8(dst, dc_dup);
  }
}

void vpx_dc_predictor_16x16_neon(uint8_t *dst, ptrdiff_t stride,
                                 const uint8_t *above, const uint8_t *left) {
  const uint8x16_t ref0 = vld1q_u8(above);
  const uint8x16_t ref1 = vld1q_u8(left);
  const uint16x8_t p0 = vpaddlq_u8(ref0);
  const uint16x8_t p1 = vpaddlq_u8(ref1);
  const uint16x8_t p2 = vaddq_u16(p0, p1);
  uint16x4_t sum = vadd_u16(vget_low_u16(p2), vget_high_u16(p2));
  uint8x8_t dc;
  sum = vpadd_u16(sum, sum);
  sum = vpadd_u16(sum, sum);
  dc = vreinterpret_u8_u16(vrshr_n_u16(sum, 5));
  dc_store_16x16(dst, stride, dc);
}

void vpx_dc_left_predictor_16x16_neon(uint8_t *dst, ptrdiff_t stride,
                                      const uint8_t *above,
                                      const uint8_t *left) {
  const uint16x4_t sum = dc_sum_16(left);
  const uint8x8_t dc = vreinterpret_u8_u16(vrshr_n_u16(sum, 4));
  (void)above;
  dc_store_16x16(dst, stride, dc);
}

void vpx_dc_top_predictor_16x16_neon(uint8_t *dst, ptrdiff_t stride,
                                     const uint8_t *above,
                                     const uint8_t *left) {
  const uint16x4_t sum = dc_sum_16(above);
  const uint8x8_t dc = vreinterpret_u8_u16(vrshr_n_u16(sum, 4));
  (void)left;
  dc_store_16x16(dst, stride, dc);
}

void vpx_dc_128_predictor_16x16_neon(uint8_t *dst, ptrdiff_t stride,
                                     const uint8_t *above,
                                     const uint8_t *left) {
  const uint8x8_t dc = vdup_n_u8(0x80);
  (void)above;
  (void)left;
  dc_store_16x16(dst, stride, dc);
}

//------------------------------------------------------------------------------
// DC 32x32

static INLINE uint16x4_t dc_sum_32(const uint8_t *ref) {
  const uint8x16x2_t r = vld2q_u8(ref);
  const uint16x8_t p0 = vpaddlq_u8(r.val[0]);
  const uint16x8_t p1 = vpaddlq_u8(r.val[1]);
  const uint16x8_t p2 = vaddq_u16(p0, p1);
  uint16x4_t sum = vadd_u16(vget_low_u16(p2), vget_high_u16(p2));
  sum = vpadd_u16(sum, sum);
  return vpadd_u16(sum, sum);
}

static INLINE void dc_store_32x32(uint8_t *dst, ptrdiff_t stride,
                                  const uint8x8_t dc) {
  uint8x16x2_t dc_dup;
  int i;
  dc_dup.val[0] = dc_dup.val[1] = vdupq_lane_u8(dc, 0);

  for (i = 0; i < 32; ++i, dst += stride) {
    vst2q_u8(dst, dc_dup);
  }
}

void vpx_dc_predictor_32x32_neon(uint8_t *dst, ptrdiff_t stride,
                                 const uint8_t *above, const uint8_t *left) {
  const uint8x16x2_t a = vld2q_u8(above);
  const uint8x16x2_t l = vld2q_u8(left);
  const uint16x8_t pa0 = vpaddlq_u8(a.val[0]);
  const uint16x8_t pl0 = vpaddlq_u8(l.val[0]);
  const uint16x8_t pa1 = vpaddlq_u8(a.val[1]);
  const uint16x8_t pl1 = vpaddlq_u8(l.val[1]);
  const uint16x8_t pa = vaddq_u16(pa0, pa1);
  const uint16x8_t pl = vaddq_u16(pl0, pl1);
  const uint16x8_t pal = vaddq_u16(pa, pl);
  uint16x4_t sum = vadd_u16(vget_low_u16(pal), vget_high_u16(pal));
  uint8x8_t dc;
  sum = vpadd_u16(sum, sum);
  sum = vpadd_u16(sum, sum);
  dc = vreinterpret_u8_u16(vrshr_n_u16(sum, 6));
  dc_store_32x32(dst, stride, dc);
}

void vpx_dc_left_predictor_32x32_neon(uint8_t *dst, ptrdiff_t stride,
                                      const uint8_t *above,
                                      const uint8_t *left) {
  const uint16x4_t sum = dc_sum_32(left);
  const uint8x8_t dc = vreinterpret_u8_u16(vrshr_n_u16(sum, 5));
  (void)above;
  dc_store_32x32(dst, stride, dc);
}

void vpx_dc_top_predictor_32x32_neon(uint8_t *dst, ptrdiff_t stride,
                                     const uint8_t *above,
                                     const uint8_t *left) {
  const uint16x4_t sum = dc_sum_32(above);
  const uint8x8_t dc = vreinterpret_u8_u16(vrshr_n_u16(sum, 5));
  (void)left;
  dc_store_32x32(dst, stride, dc);
}

void vpx_dc_128_predictor_32x32_neon(uint8_t *dst, ptrdiff_t stride,
                                     const uint8_t *above,
                                     const uint8_t *left) {
  const uint8x8_t dc = vdup_n_u8(0x80);
  (void)above;
  (void)left;
  dc_store_32x32(dst, stride, dc);
}

// -----------------------------------------------------------------------------

void vpx_d45_predictor_4x4_neon(uint8_t *dst, ptrdiff_t stride,
                                const uint8_t *above, const uint8_t *left) {
  const uint8x8_t ABCDEFGH = vld1_u8(above);
  const uint64x1_t A1 = vshr_n_u64(vreinterpret_u64_u8(ABCDEFGH), 8);
  const uint64x1_t A2 = vshr_n_u64(vreinterpret_u64_u8(ABCDEFGH), 16);
  const uint8x8_t BCDEFGH0 = vreinterpret_u8_u64(A1);
  const uint8x8_t CDEFGH00 = vreinterpret_u8_u64(A2);
  const uint8x8_t avg1 = vhadd_u8(ABCDEFGH, CDEFGH00);
  const uint8x8_t avg2 = vrhadd_u8(avg1, BCDEFGH0);
  const uint64x1_t avg2_u64 = vreinterpret_u64_u8(avg2);
  const uint32x2_t r0 = vreinterpret_u32_u8(avg2);
  const uint32x2_t r1 = vreinterpret_u32_u64(vshr_n_u64(avg2_u64, 8));
  const uint32x2_t r2 = vreinterpret_u32_u64(vshr_n_u64(avg2_u64, 16));
  const uint32x2_t r3 = vreinterpret_u32_u64(vshr_n_u64(avg2_u64, 24));
  (void)left;
  vst1_lane_u32((uint32_t *)(dst + 0 * stride), r0, 0);
  vst1_lane_u32((uint32_t *)(dst + 1 * stride), r1, 0);
  vst1_lane_u32((uint32_t *)(dst + 2 * stride), r2, 0);
  vst1_lane_u32((uint32_t *)(dst + 3 * stride), r3, 0);
  vst1_lane_u8(dst + 3 * stride + 3, ABCDEFGH, 7);
}

static INLINE void d45_store_8(uint8_t **dst, const ptrdiff_t stride,
                               const uint8x8_t above_right, uint8x8_t *row) {
  *row = vext_u8(*row, above_right, 1);
  vst1_u8(*dst, *row);
  *dst += stride;
}

void vpx_d45_predictor_8x8_neon(uint8_t *dst, ptrdiff_t stride,
                                const uint8_t *above, const uint8_t *left) {
  const uint8x8_t A0 = vld1_u8(above);
  const uint8x8_t above_right = vdup_lane_u8(A0, 7);
  const uint8x8_t A1 = vext_u8(A0, above_right, 1);
  const uint8x8_t A2 = vext_u8(A0, above_right, 2);
  const uint8x8_t avg1 = vhadd_u8(A0, A2);
  uint8x8_t row = vrhadd_u8(avg1, A1);
  (void)left;

  vst1_u8(dst, row);
  dst += stride;
  d45_store_8(&dst, stride, above_right, &row);
  d45_store_8(&dst, stride, above_right, &row);
  d45_store_8(&dst, stride, above_right, &row);
  d45_store_8(&dst, stride, above_right, &row);
  d45_store_8(&dst, stride, above_right, &row);
  d45_store_8(&dst, stride, above_right, &row);
  vst1_u8(dst, above_right);
}

static INLINE void d45_store_16(uint8_t **dst, const ptrdiff_t stride,
                                const uint8x16_t above_right, uint8x16_t *row) {
  *row = vextq_u8(*row, above_right, 1);
  vst1q_u8(*dst, *row);
  *dst += stride;
}

void vpx_d45_predictor_16x16_neon(uint8_t *dst, ptrdiff_t stride,
                                  const uint8_t *above, const uint8_t *left) {
  const uint8x16_t A0 = vld1q_u8(above);
  const uint8x16_t above_right = vdupq_lane_u8(vget_high_u8(A0), 7);
  const uint8x16_t A1 = vextq_u8(A0, above_right, 1);
  const uint8x16_t A2 = vextq_u8(A0, above_right, 2);
  const uint8x16_t avg1 = vhaddq_u8(A0, A2);
  uint8x16_t row = vrhaddq_u8(avg1, A1);
  (void)left;

  vst1q_u8(dst, row);
  dst += stride;
  d45_store_16(&dst, stride, above_right, &row);
  d45_store_16(&dst, stride, above_right, &row);
  d45_store_16(&dst, stride, above_right, &row);
  d45_store_16(&dst, stride, above_right, &row);
  d45_store_16(&dst, stride, above_right, &row);
  d45_store_16(&dst, stride, above_right, &row);
  d45_store_16(&dst, stride, above_right, &row);
  d45_store_16(&dst, stride, above_right, &row);
  d45_store_16(&dst, stride, above_right, &row);
  d45_store_16(&dst, stride, above_right, &row);
  d45_store_16(&dst, stride, above_right, &row);
  d45_store_16(&dst, stride, above_right, &row);
  d45_store_16(&dst, stride, above_right, &row);
  d45_store_16(&dst, stride, above_right, &row);
  vst1q_u8(dst, above_right);
}

// -----------------------------------------------------------------------------

void vpx_d135_predictor_4x4_neon(uint8_t *dst, ptrdiff_t stride,
                                 const uint8_t *above, const uint8_t *left) {
  const uint8x8_t XABCD = vld1_u8(above - 1);
  const uint32x2_t zero = vdup_n_u32(0);
  const uint32x2_t IJKL = vld1_lane_u32((const uint32_t *)left, zero, 0);
  const uint8x8_t LKJI = vrev64_u8(vreinterpret_u8_u32(IJKL));
  const uint8x8_t LKJIXABC = vext_u8(LKJI, XABCD, 4);
  const uint8x8_t KJIXABCD = vext_u8(LKJI, XABCD, 5);
  const uint8x8_t JIXABCD0 =
      vreinterpret_u8_u64(vshr_n_u64(vreinterpret_u64_u8(KJIXABCD), 8));
  const uint8x8_t avg1 = vhadd_u8(JIXABCD0, LKJIXABC);
  const uint8x8_t avg2 = vrhadd_u8(avg1, KJIXABCD);
  const uint64x1_t avg2_u64 = vreinterpret_u64_u8(avg2);
  const uint32x2_t r3 = vreinterpret_u32_u8(avg2);
  const uint32x2_t r2 = vreinterpret_u32_u64(vshr_n_u64(avg2_u64, 8));
  const uint32x2_t r1 = vreinterpret_u32_u64(vshr_n_u64(avg2_u64, 16));
  const uint32x2_t r0 = vreinterpret_u32_u64(vshr_n_u64(avg2_u64, 24));
  vst1_lane_u32((uint32_t *)dst, r0, 0);
  dst += stride;
  vst1_lane_u32((uint32_t *)dst, r1, 0);
  dst += stride;
  vst1_lane_u32((uint32_t *)dst, r2, 0);
  dst += stride;
  vst1_lane_u32((uint32_t *)dst, r3, 0);
}

// -----------------------------------------------------------------------------

#if !HAVE_NEON_ASM

void vpx_v_predictor_4x4_neon(uint8_t *dst, ptrdiff_t stride,
                              const uint8_t *above, const uint8_t *left) {
  const uint32_t d = *(const uint32_t *)above;
  int i;
  (void)left;

  for (i = 0; i < 4; i++, dst += stride) {
    *(uint32_t *)dst = d;
  }
}

void vpx_v_predictor_8x8_neon(uint8_t *dst, ptrdiff_t stride,
                              const uint8_t *above, const uint8_t *left) {
  const uint8x8_t d = vld1_u8(above);
  int i;
  (void)left;

  for (i = 0; i < 8; i++, dst += stride) {
    vst1_u8(dst, d);
  }
}

void vpx_v_predictor_16x16_neon(uint8_t *dst, ptrdiff_t stride,
                                const uint8_t *above, const uint8_t *left) {
  const uint8x16_t d = vld1q_u8(above);
  int i;
  (void)left;

  for (i = 0; i < 16; i++, dst += stride) {
    vst1q_u8(dst, d);
  }
}

void vpx_v_predictor_32x32_neon(uint8_t *dst, ptrdiff_t stride,
                                const uint8_t *above, const uint8_t *left) {
  const uint8x16_t d0 = vld1q_u8(above);
  const uint8x16_t d1 = vld1q_u8(above + 16);
  int i;
  (void)left;

  for (i = 0; i < 32; i++) {
    // Note: performance was worse using vst2q_u8 under gcc-4.9 & clang-3.8.
    // clang-3.8 unrolled the loop fully with no filler so the cause is likely
    // the latency of the instruction.
    vst1q_u8(dst, d0);
    dst += 16;
    vst1q_u8(dst, d1);
    dst += stride - 16;
  }
}

// -----------------------------------------------------------------------------

void vpx_h_predictor_4x4_neon(uint8_t *dst, ptrdiff_t stride,
                              const uint8_t *above, const uint8_t *left) {
  const uint32x2_t zero = vdup_n_u32(0);
  const uint8x8_t left_u8 =
      vreinterpret_u8_u32(vld1_lane_u32((const uint32_t *)left, zero, 0));
  uint8x8_t d;
  (void)above;

  d = vdup_lane_u8(left_u8, 0);
  vst1_lane_u32((uint32_t *)dst, vreinterpret_u32_u8(d), 0);
  dst += stride;
  d = vdup_lane_u8(left_u8, 1);
  vst1_lane_u32((uint32_t *)dst, vreinterpret_u32_u8(d), 0);
  dst += stride;
  d = vdup_lane_u8(left_u8, 2);
  vst1_lane_u32((uint32_t *)dst, vreinterpret_u32_u8(d), 0);
  dst += stride;
  d = vdup_lane_u8(left_u8, 3);
  vst1_lane_u32((uint32_t *)dst, vreinterpret_u32_u8(d), 0);
}

void vpx_h_predictor_8x8_neon(uint8_t *dst, ptrdiff_t stride,
                              const uint8_t *above, const uint8_t *left) {
  const uint8x8_t left_u8 = vld1_u8(left);
  uint8x8_t d;
  (void)above;

  d = vdup_lane_u8(left_u8, 0);
  vst1_u8(dst, d);
  dst += stride;
  d = vdup_lane_u8(left_u8, 1);
  vst1_u8(dst, d);
  dst += stride;
  d = vdup_lane_u8(left_u8, 2);
  vst1_u8(dst, d);
  dst += stride;
  d = vdup_lane_u8(left_u8, 3);
  vst1_u8(dst, d);
  dst += stride;
  d = vdup_lane_u8(left_u8, 4);
  vst1_u8(dst, d);
  dst += stride;
  d = vdup_lane_u8(left_u8, 5);
  vst1_u8(dst, d);
  dst += stride;
  d = vdup_lane_u8(left_u8, 6);
  vst1_u8(dst, d);
  dst += stride;
  d = vdup_lane_u8(left_u8, 7);
  vst1_u8(dst, d);
}

void vpx_h_predictor_16x16_neon(uint8_t *dst, ptrdiff_t stride,
                                const uint8_t *above, const uint8_t *left) {
  const uint8x16_t left_u8q = vld1q_u8(left);
  uint8x8_t left_u8d = vget_low_u8(left_u8q);
  uint8x16_t d;
  int i;
  (void)above;

  for (i = 0; i < 2; i++, left_u8d = vget_high_u8(left_u8q)) {
    d = vdupq_lane_u8(left_u8d, 0);
    vst1q_u8(dst, d);
    dst += stride;
    d = vdupq_lane_u8(left_u8d, 1);
    vst1q_u8(dst, d);
    dst += stride;
    d = vdupq_lane_u8(left_u8d, 2);
    vst1q_u8(dst, d);
    dst += stride;
    d = vdupq_lane_u8(left_u8d, 3);
    vst1q_u8(dst, d);
    dst += stride;
    d = vdupq_lane_u8(left_u8d, 4);
    vst1q_u8(dst, d);
    dst += stride;
    d = vdupq_lane_u8(left_u8d, 5);
    vst1q_u8(dst, d);
    dst += stride;
    d = vdupq_lane_u8(left_u8d, 6);
    vst1q_u8(dst, d);
    dst += stride;
    d = vdupq_lane_u8(left_u8d, 7);
    vst1q_u8(dst, d);
    dst += stride;
  }
}

void vpx_h_predictor_32x32_neon(uint8_t *dst, ptrdiff_t stride,
                                const uint8_t *above, const uint8_t *left) {
  uint8x16_t d;
  int i;
  (void)above;

  for (i = 0; i < 2; i++, left += 16) {
    const uint8x16_t left_u8 = vld1q_u8(left);
    const uint8x8_t left_low = vget_low_u8(left_u8);
    const uint8x8_t left_high = vget_high_u8(left_u8);
    d = vdupq_lane_u8(left_low, 0);
    vst1q_u8(dst, d);  // Note clang-3.8 produced poor code w/vst2q_u8
    dst += 16;
    vst1q_u8(dst, d);
    dst += stride - 16;
    d = vdupq_lane_u8(left_low, 1);
    vst1q_u8(dst, d);
    dst += 16;
    vst1q_u8(dst, d);
    dst += stride - 16;
    d = vdupq_lane_u8(left_low, 2);
    vst1q_u8(dst, d);
    dst += 16;
    vst1q_u8(dst, d);
    dst += stride - 16;
    d = vdupq_lane_u8(left_low, 3);
    vst1q_u8(dst, d);
    dst += 16;
    vst1q_u8(dst, d);
    dst += stride - 16;
    d = vdupq_lane_u8(left_low, 4);
    vst1q_u8(dst, d);
    dst += 16;
    vst1q_u8(dst, d);
    dst += stride - 16;
    d = vdupq_lane_u8(left_low, 5);
    vst1q_u8(dst, d);
    dst += 16;
    vst1q_u8(dst, d);
    dst += stride - 16;
    d = vdupq_lane_u8(left_low, 6);
    vst1q_u8(dst, d);
    dst += 16;
    vst1q_u8(dst, d);
    dst += stride - 16;
    d = vdupq_lane_u8(left_low, 7);
    vst1q_u8(dst, d);
    dst += 16;
    vst1q_u8(dst, d);
    dst += stride - 16;

    d = vdupq_lane_u8(left_high, 0);
    vst1q_u8(dst, d);
    dst += 16;
    vst1q_u8(dst, d);
    dst += stride - 16;
    d = vdupq_lane_u8(left_high, 1);
    vst1q_u8(dst, d);
    dst += 16;
    vst1q_u8(dst, d);
    dst += stride - 16;
    d = vdupq_lane_u8(left_high, 2);
    vst1q_u8(dst, d);
    dst += 16;
    vst1q_u8(dst, d);
    dst += stride - 16;
    d = vdupq_lane_u8(left_high, 3);
    vst1q_u8(dst, d);
    dst += 16;
    vst1q_u8(dst, d);
    dst += stride - 16;
    d = vdupq_lane_u8(left_high, 4);
    vst1q_u8(dst, d);
    dst += 16;
    vst1q_u8(dst, d);
    dst += stride - 16;
    d = vdupq_lane_u8(left_high, 5);
    vst1q_u8(dst, d);
    dst += 16;
    vst1q_u8(dst, d);
    dst += stride - 16;
    d = vdupq_lane_u8(left_high, 6);
    vst1q_u8(dst, d);
    dst += 16;
    vst1q_u8(dst, d);
    dst += stride - 16;
    d = vdupq_lane_u8(left_high, 7);
    vst1q_u8(dst, d);
    dst += 16;
    vst1q_u8(dst, d);
    dst += stride - 16;
  }
}

// -----------------------------------------------------------------------------

static INLINE int16x8_t convert_u8_to_s16(uint8x8_t v) {
  return vreinterpretq_s16_u16(vmovl_u8(v));
}

void vpx_tm_predictor_4x4_neon(uint8_t *dst, ptrdiff_t stride,
                               const uint8_t *above, const uint8_t *left) {
  const uint8x8_t top_left = vld1_dup_u8(above - 1);
  const uint8x8_t left_u8 = vld1_u8(left);
  const uint8x8_t above_u8 = vld1_u8(above);
  const int16x4_t left_s16 = vget_low_s16(convert_u8_to_s16(left_u8));
  int16x8_t sub, sum;
  uint32x2_t d;

  sub = vreinterpretq_s16_u16(vsubl_u8(above_u8, top_left));
  // Avoid vcombine_s16() which generates lots of redundant code with clang-3.8.
  sub = vreinterpretq_s16_s64(
      vdupq_lane_s64(vreinterpret_s64_s16(vget_low_s16(sub)), 0));

  sum = vcombine_s16(vdup_lane_s16(left_s16, 0), vdup_lane_s16(left_s16, 1));
  sum = vaddq_s16(sum, sub);
  d = vreinterpret_u32_u8(vqmovun_s16(sum));
  vst1_lane_u32((uint32_t *)dst, d, 0);
  dst += stride;
  vst1_lane_u32((uint32_t *)dst, d, 1);
  dst += stride;

  sum = vcombine_s16(vdup_lane_s16(left_s16, 2), vdup_lane_s16(left_s16, 3));
  sum = vaddq_s16(sum, sub);
  d = vreinterpret_u32_u8(vqmovun_s16(sum));
  vst1_lane_u32((uint32_t *)dst, d, 0);
  dst += stride;
  vst1_lane_u32((uint32_t *)dst, d, 1);
}

static INLINE void tm_8_kernel(uint8_t **dst, const ptrdiff_t stride,
                               const int16x8_t left_dup, const int16x8_t sub) {
  const int16x8_t sum = vaddq_s16(left_dup, sub);
  const uint8x8_t d = vqmovun_s16(sum);
  vst1_u8(*dst, d);
  *dst += stride;
}

void vpx_tm_predictor_8x8_neon(uint8_t *dst, ptrdiff_t stride,
                               const uint8_t *above, const uint8_t *left) {
  const uint8x8_t top_left = vld1_dup_u8(above - 1);
  const uint8x8_t above_u8 = vld1_u8(above);
  const uint8x8_t left_u8 = vld1_u8(left);
  const int16x8_t left_s16q = convert_u8_to_s16(left_u8);
  const int16x8_t sub = vreinterpretq_s16_u16(vsubl_u8(above_u8, top_left));
  int16x4_t left_s16d = vget_low_s16(left_s16q);
  int i;

  for (i = 0; i < 2; i++, left_s16d = vget_high_s16(left_s16q)) {
    int16x8_t left_dup;

    left_dup = vdupq_lane_s16(left_s16d, 0);
    tm_8_kernel(&dst, stride, left_dup, sub);
    left_dup = vdupq_lane_s16(left_s16d, 1);
    tm_8_kernel(&dst, stride, left_dup, sub);
    left_dup = vdupq_lane_s16(left_s16d, 2);
    tm_8_kernel(&dst, stride, left_dup, sub);
    left_dup = vdupq_lane_s16(left_s16d, 3);
    tm_8_kernel(&dst, stride, left_dup, sub);
  }
}

static INLINE void tm_16_kernel(uint8_t **dst, const ptrdiff_t stride,
                                const int16x8_t left_dup, const int16x8_t sub0,
                                const int16x8_t sub1) {
  const int16x8_t sum0 = vaddq_s16(left_dup, sub0);
  const int16x8_t sum1 = vaddq_s16(left_dup, sub1);
  const uint8x8_t d0 = vqmovun_s16(sum0);
  const uint8x8_t d1 = vqmovun_s16(sum1);
  vst1_u8(*dst, d0);
  *dst += 8;
  vst1_u8(*dst, d1);
  *dst += stride - 8;
}

void vpx_tm_predictor_16x16_neon(uint8_t *dst, ptrdiff_t stride,
                                 const uint8_t *above, const uint8_t *left) {
  const uint8x16_t top_left = vld1q_dup_u8(above - 1);
  const uint8x16_t above_u8 = vld1q_u8(above);
  const int16x8_t sub0 = vreinterpretq_s16_u16(
      vsubl_u8(vget_low_u8(above_u8), vget_low_u8(top_left)));
  const int16x8_t sub1 = vreinterpretq_s16_u16(
      vsubl_u8(vget_high_u8(above_u8), vget_high_u8(top_left)));
  int16x8_t left_dup;
  int i;

  for (i = 0; i < 2; i++, left += 8) {
    const uint8x8_t left_u8 = vld1_u8(left);
    const int16x8_t left_s16q = convert_u8_to_s16(left_u8);
    const int16x4_t left_low = vget_low_s16(left_s16q);
    const int16x4_t left_high = vget_high_s16(left_s16q);

    left_dup = vdupq_lane_s16(left_low, 0);
    tm_16_kernel(&dst, stride, left_dup, sub0, sub1);
    left_dup = vdupq_lane_s16(left_low, 1);
    tm_16_kernel(&dst, stride, left_dup, sub0, sub1);
    left_dup = vdupq_lane_s16(left_low, 2);
    tm_16_kernel(&dst, stride, left_dup, sub0, sub1);
    left_dup = vdupq_lane_s16(left_low, 3);
    tm_16_kernel(&dst, stride, left_dup, sub0, sub1);

    left_dup = vdupq_lane_s16(left_high, 0);
    tm_16_kernel(&dst, stride, left_dup, sub0, sub1);
    left_dup = vdupq_lane_s16(left_high, 1);
    tm_16_kernel(&dst, stride, left_dup, sub0, sub1);
    left_dup = vdupq_lane_s16(left_high, 2);
    tm_16_kernel(&dst, stride, left_dup, sub0, sub1);
    left_dup = vdupq_lane_s16(left_high, 3);
    tm_16_kernel(&dst, stride, left_dup, sub0, sub1);
  }
}

static INLINE void tm_32_kernel(uint8_t **dst, const ptrdiff_t stride,
                                const int16x8_t left_dup, const int16x8_t sub0,
                                const int16x8_t sub1, const int16x8_t sub2,
                                const int16x8_t sub3) {
  const int16x8_t sum0 = vaddq_s16(left_dup, sub0);
  const int16x8_t sum1 = vaddq_s16(left_dup, sub1);
  const int16x8_t sum2 = vaddq_s16(left_dup, sub2);
  const int16x8_t sum3 = vaddq_s16(left_dup, sub3);
  const uint8x8_t d0 = vqmovun_s16(sum0);
  const uint8x8_t d1 = vqmovun_s16(sum1);
  const uint8x8_t d2 = vqmovun_s16(sum2);
  const uint8x8_t d3 = vqmovun_s16(sum3);

  vst1q_u8(*dst, vcombine_u8(d0, d1));
  *dst += 16;
  vst1q_u8(*dst, vcombine_u8(d2, d3));
  *dst += stride - 16;
}

void vpx_tm_predictor_32x32_neon(uint8_t *dst, ptrdiff_t stride,
                                 const uint8_t *above, const uint8_t *left) {
  const uint8x16_t top_left = vld1q_dup_u8(above - 1);
  const uint8x16_t above_low = vld1q_u8(above);
  const uint8x16_t above_high = vld1q_u8(above + 16);
  const int16x8_t sub0 = vreinterpretq_s16_u16(
      vsubl_u8(vget_low_u8(above_low), vget_low_u8(top_left)));
  const int16x8_t sub1 = vreinterpretq_s16_u16(
      vsubl_u8(vget_high_u8(above_low), vget_high_u8(top_left)));
  const int16x8_t sub2 = vreinterpretq_s16_u16(
      vsubl_u8(vget_low_u8(above_high), vget_low_u8(top_left)));
  const int16x8_t sub3 = vreinterpretq_s16_u16(
      vsubl_u8(vget_high_u8(above_high), vget_high_u8(top_left)));
  int16x8_t left_dup;
  int i, j;

  for (j = 0; j < 4; j++, left += 8) {
    const uint8x8_t left_u8 = vld1_u8(left);
    const int16x8_t left_s16q = convert_u8_to_s16(left_u8);
    int16x4_t left_s16d = vget_low_s16(left_s16q);
    for (i = 0; i < 2; i++, left_s16d = vget_high_s16(left_s16q)) {
      left_dup = vdupq_lane_s16(left_s16d, 0);
      tm_32_kernel(&dst, stride, left_dup, sub0, sub1, sub2, sub3);
      left_dup = vdupq_lane_s16(left_s16d, 1);
      tm_32_kernel(&dst, stride, left_dup, sub0, sub1, sub2, sub3);
      left_dup = vdupq_lane_s16(left_s16d, 2);
      tm_32_kernel(&dst, stride, left_dup, sub0, sub1, sub2, sub3);
      left_dup = vdupq_lane_s16(left_s16d, 3);
      tm_32_kernel(&dst, stride, left_dup, sub0, sub1, sub2, sub3);
    }
  }
}
#endif  // !HAVE_NEON_ASM
