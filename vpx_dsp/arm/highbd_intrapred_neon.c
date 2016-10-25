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
#include "vpx/vpx_integer.h"

//------------------------------------------------------------------------------
// DC 4x4

static INLINE uint16x4_t dc_sum_4(const uint16_t *ref) {
  const uint16x4_t ref_u16 = vld1_u16(ref);
  const uint16x4_t p0 = vpadd_u16(ref_u16, ref_u16);
  return vpadd_u16(p0, p0);
}

static INLINE void dc_store_4x4(uint16_t *dst, ptrdiff_t stride,
                                const uint16x4_t dc) {
  const uint16x4_t dc_dup = vdup_lane_u16(dc, 0);
  int i;
  for (i = 0; i < 4; ++i, dst += stride) {
    vst1_u16(dst, dc_dup);
  }
}

void vpx_highbd_dc_predictor_4x4_neon(uint16_t *dst, ptrdiff_t stride,
                                      const uint16_t *above,
                                      const uint16_t *left, int bd) {
  const uint16x4_t a = vld1_u16(above);
  const uint16x4_t l = vld1_u16(left);
  uint16x4_t sum;
  uint16x4_t dc;
  (void)bd;
  sum = vadd_u16(a, l);
  sum = vpadd_u16(sum, sum);
  sum = vpadd_u16(sum, sum);
  dc = vrshr_n_u16(sum, 3);
  dc_store_4x4(dst, stride, dc);
}

void vpx_highbd_dc_left_predictor_4x4_neon(uint16_t *dst, ptrdiff_t stride,
                                           const uint16_t *above,
                                           const uint16_t *left, int bd) {
  const uint16x4_t sum = dc_sum_4(left);
  const uint16x4_t dc = vrshr_n_u16(sum, 2);
  (void)above;
  (void)bd;
  dc_store_4x4(dst, stride, dc);
}

void vpx_highbd_dc_top_predictor_4x4_neon(uint16_t *dst, ptrdiff_t stride,
                                          const uint16_t *above,
                                          const uint16_t *left, int bd) {
  const uint16x4_t sum = dc_sum_4(above);
  const uint16x4_t dc = vrshr_n_u16(sum, 2);
  (void)left;
  (void)bd;
  dc_store_4x4(dst, stride, dc);
}

void vpx_highbd_dc_128_predictor_4x4_neon(uint16_t *dst, ptrdiff_t stride,
                                          const uint16_t *above,
                                          const uint16_t *left, int bd) {
  const uint16x4_t dc = vdup_n_u16(1 << (bd - 1));
  (void)above;
  (void)left;
  dc_store_4x4(dst, stride, dc);
}

//------------------------------------------------------------------------------
// DC 8x8

static INLINE uint16x4_t dc_sum_8(const uint16_t *ref) {
  const uint16x8_t ref_u16 = vld1q_u16(ref);
  uint16x4_t sum = vadd_u16(vget_low_u16(ref_u16), vget_high_u16(ref_u16));
  sum = vpadd_u16(sum, sum);
  return vpadd_u16(sum, sum);
}

static INLINE void dc_store_8x8(uint16_t *dst, ptrdiff_t stride,
                                const uint16x4_t dc) {
  const uint16x8_t dc_dup = vdupq_lane_u16(dc, 0);
  int i;
  for (i = 0; i < 8; ++i, dst += stride) {
    vst1q_u16(dst, dc_dup);
  }
}

void vpx_highbd_dc_predictor_8x8_neon(uint16_t *dst, ptrdiff_t stride,
                                      const uint16_t *above,
                                      const uint16_t *left, int bd) {
  const uint16x8_t above_u16 = vld1q_u16(above);
  const uint16x8_t left_u16 = vld1q_u16(left);
  const uint16x8_t p0 = vaddq_u16(above_u16, left_u16);
  uint16x4_t sum = vadd_u16(vget_low_u16(p0), vget_high_u16(p0));
  uint16x4_t dc;
  (void)bd;
  sum = vpadd_u16(sum, sum);
  sum = vpadd_u16(sum, sum);
  dc = vrshr_n_u16(sum, 4);
  dc_store_8x8(dst, stride, dc);
}

void vpx_highbd_dc_left_predictor_8x8_neon(uint16_t *dst, ptrdiff_t stride,
                                           const uint16_t *above,
                                           const uint16_t *left, int bd) {
  const uint16x4_t sum = dc_sum_8(left);
  const uint16x4_t dc = vrshr_n_u16(sum, 3);
  (void)above;
  (void)bd;
  dc_store_8x8(dst, stride, dc);
}

void vpx_highbd_dc_top_predictor_8x8_neon(uint16_t *dst, ptrdiff_t stride,
                                          const uint16_t *above,
                                          const uint16_t *left, int bd) {
  const uint16x4_t sum = dc_sum_8(above);
  const uint16x4_t dc = vrshr_n_u16(sum, 3);
  (void)left;
  (void)bd;
  dc_store_8x8(dst, stride, dc);
}

void vpx_highbd_dc_128_predictor_8x8_neon(uint16_t *dst, ptrdiff_t stride,
                                          const uint16_t *above,
                                          const uint16_t *left, int bd) {
  const uint16x4_t dc = vdup_n_u16(1 << (bd - 1));
  (void)above;
  (void)left;
  dc_store_8x8(dst, stride, dc);
}

//------------------------------------------------------------------------------
// DC 16x16

static INLINE uint16x4_t dc_sum_16(const uint16_t *ref) {
  const uint16x8x2_t ref_u16 = vld2q_u16(ref);
  const uint16x8_t p0 = vaddq_u16(ref_u16.val[0], ref_u16.val[1]);
  uint16x4_t sum = vadd_u16(vget_low_u16(p0), vget_high_u16(p0));
  sum = vpadd_u16(sum, sum);
  return vpadd_u16(sum, sum);
}

static INLINE void dc_store_16x16(uint16_t *dst, ptrdiff_t stride,
                                  const uint16x4_t dc) {
  uint16x8x2_t dc_dup;
  int i;
  dc_dup.val[0] = dc_dup.val[1] = vdupq_lane_u16(dc, 0);
  for (i = 0; i < 16; ++i, dst += stride) {
    vst2q_u16(dst, dc_dup);
  }
}

void vpx_highbd_dc_predictor_16x16_neon(uint16_t *dst, ptrdiff_t stride,
                                        const uint16_t *above,
                                        const uint16_t *left, int bd) {
  const uint16x8x2_t a = vld2q_u16(above);
  const uint16x8x2_t l = vld2q_u16(left);
  const uint16x8_t pa = vaddq_u16(a.val[0], a.val[1]);
  const uint16x8_t pl = vaddq_u16(l.val[0], l.val[1]);
  const uint16x8_t pal0 = vaddq_u16(pa, pl);
  uint16x4_t pal1 = vadd_u16(vget_low_u16(pal0), vget_high_u16(pal0));
  uint32x2_t sum;
  uint16x4_t dc;
  (void)bd;
  pal1 = vpadd_u16(pal1, pal1);
  sum = vpaddl_u16(pal1);
  dc = vreinterpret_u16_u32(vrshr_n_u32(sum, 5));
  dc_store_16x16(dst, stride, dc);
}

void vpx_highbd_dc_left_predictor_16x16_neon(uint16_t *dst, ptrdiff_t stride,
                                             const uint16_t *above,
                                             const uint16_t *left, int bd) {
  const uint16x4_t sum = dc_sum_16(left);
  const uint16x4_t dc = vrshr_n_u16(sum, 4);
  (void)above;
  (void)bd;
  dc_store_16x16(dst, stride, dc);
}

void vpx_highbd_dc_top_predictor_16x16_neon(uint16_t *dst, ptrdiff_t stride,
                                            const uint16_t *above,
                                            const uint16_t *left, int bd) {
  const uint16x4_t sum = dc_sum_16(above);
  const uint16x4_t dc = vrshr_n_u16(sum, 4);
  (void)left;
  (void)bd;
  dc_store_16x16(dst, stride, dc);
}

void vpx_highbd_dc_128_predictor_16x16_neon(uint16_t *dst, ptrdiff_t stride,
                                            const uint16_t *above,
                                            const uint16_t *left, int bd) {
  const uint16x4_t dc = vdup_n_u16(1 << (bd - 1));
  (void)above;
  (void)left;
  dc_store_16x16(dst, stride, dc);
}

//------------------------------------------------------------------------------
// DC 32x32

static INLINE uint32x2_t dc_sum_32(const uint16_t *ref) {
  const uint16x8x4_t r = vld4q_u16(ref);
  const uint16x8_t p0 = vaddq_u16(r.val[0], r.val[1]);
  const uint16x8_t p1 = vaddq_u16(r.val[2], r.val[3]);
  const uint16x8_t p2 = vaddq_u16(p0, p1);
  uint16x4_t sum = vadd_u16(vget_low_u16(p2), vget_high_u16(p2));
  sum = vpadd_u16(sum, sum);
  return vpaddl_u16(sum);
}

static INLINE void dc_store_32x32(uint16_t *dst, ptrdiff_t stride,
                                  const uint16x4_t dc) {
  uint16x8x2_t dc_dup;
  int i;
  dc_dup.val[0] = dc_dup.val[1] = vdupq_lane_u16(dc, 0);

  for (i = 0; i < 32; ++i) {
    vst2q_u16(dst, dc_dup);
    dst += 16;
    vst2q_u16(dst, dc_dup);
    dst += stride - 16;
  }
}

void vpx_highbd_dc_predictor_32x32_neon(uint16_t *dst, ptrdiff_t stride,
                                        const uint16_t *above,
                                        const uint16_t *left, int bd) {
  const uint16x8x4_t a = vld4q_u16(above);
  const uint16x8x4_t l = vld4q_u16(left);
  const uint16x8_t pa0 = vaddq_u16(a.val[0], a.val[1]);
  const uint16x8_t pa1 = vaddq_u16(a.val[2], a.val[3]);
  const uint16x8_t pl0 = vaddq_u16(l.val[0], l.val[1]);
  const uint16x8_t pl1 = vaddq_u16(l.val[2], l.val[3]);
  const uint16x8_t pa = vaddq_u16(pa0, pa1);
  const uint16x8_t pl = vaddq_u16(pl0, pl1);
  const uint16x8_t pal0 = vaddq_u16(pa, pl);
  const uint16x4_t pal1 = vadd_u16(vget_low_u16(pal0), vget_high_u16(pal0));
  uint32x2_t sum = vpaddl_u16(pal1);
  uint16x4_t dc;
  (void)bd;
  sum = vpadd_u32(sum, sum);
  dc = vreinterpret_u16_u32(vrshr_n_u32(sum, 6));
  dc_store_32x32(dst, stride, dc);
}

void vpx_highbd_dc_left_predictor_32x32_neon(uint16_t *dst, ptrdiff_t stride,
                                             const uint16_t *above,
                                             const uint16_t *left, int bd) {
  const uint32x2_t sum = dc_sum_32(left);
  const uint16x4_t dc = vreinterpret_u16_u32(vrshr_n_u32(sum, 5));
  (void)above;
  (void)bd;
  dc_store_32x32(dst, stride, dc);
}

void vpx_highbd_dc_top_predictor_32x32_neon(uint16_t *dst, ptrdiff_t stride,
                                            const uint16_t *above,
                                            const uint16_t *left, int bd) {
  const uint32x2_t sum = dc_sum_32(above);
  const uint16x4_t dc = vreinterpret_u16_u32(vrshr_n_u32(sum, 5));
  (void)left;
  (void)bd;
  dc_store_32x32(dst, stride, dc);
}

void vpx_highbd_dc_128_predictor_32x32_neon(uint16_t *dst, ptrdiff_t stride,
                                            const uint16_t *above,
                                            const uint16_t *left, int bd) {
  const uint16x4_t dc = vdup_n_u16(1 << (bd - 1));
  (void)above;
  (void)left;
  dc_store_32x32(dst, stride, dc);
}
