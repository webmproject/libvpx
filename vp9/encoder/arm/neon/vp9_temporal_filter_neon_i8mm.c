/*
 *  Copyright (c) 2025 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include <assert.h>
#include <arm_neon.h>

#include "./vp9_rtcd.h"
#include "./vpx_config.h"
#include "vpx/vpx_integer.h"
#include "vpx_dsp/arm/mem_neon.h"
#include "vp9/encoder/vp9_temporal_filter.h"

DECLARE_ALIGNED(16, static const uint8_t, kMatMulPermuteTbl[32]) = {
  // clang-format off
  0,  1,  2,  3,  4,  5,  6,  7,  2,  3,  4,  5,  6,  7,  8,  9,
  4,  5,  6,  7,  8,  9, 10, 11,  6,  7,  8,  9, 10, 11, 12, 13
  // clang-format on
};

static INLINE uint8x8_t convolve12_8_h(uint8x16_t samples[2],
                                       const int8x16_t filter[2],
                                       const uint8x16x2_t perm_tbl) {
  // Permute samples ready for matrix multiply.
  // {  0,  1,  2,  3,  4,  5,  6,  7,  2,  3,  4,  5,  6,  7,  8,  9 }
  // {  4,  5,  6,  7,  8,  9, 10, 11,  6,  7,  8,  9, 10, 11, 12, 13 }
  // {  6,  7,  8,  9, 10, 11, 12, 13,  8,  9, 10, 11, 12, 13, 14, 15 }
  // { 10, 11, 12, 13, 14, 15, 16, 17, 12, 13, 14, 15, 16, 17, 18, 19 }
  uint8x16_t perm_samples[4] = { vqtbl1q_u8(samples[0], perm_tbl.val[0]),
                                 vqtbl1q_u8(samples[0], perm_tbl.val[1]),
                                 vqtbl1q_u8(samples[1], perm_tbl.val[0]),
                                 vqtbl1q_u8(samples[1], perm_tbl.val[1]) };

  // These instructions multiply a 2x8 matrix (samples) by an 8x2 matrix
  // (filter), destructively accumulating into the destination register.
  int32x4_t sum0123 = vusmmlaq_s32(vdupq_n_s32(0), perm_samples[0], filter[0]);
  int32x4_t sum4567 = vusmmlaq_s32(vdupq_n_s32(0), perm_samples[1], filter[0]);
  sum0123 = vusmmlaq_s32(sum0123, perm_samples[2], filter[1]);
  sum4567 = vusmmlaq_s32(sum4567, perm_samples[3], filter[1]);

  // Narrow and re-pack.
  int16x8_t sum_s16 = vcombine_s16(vqrshrn_n_s32(sum0123, FILTER_BITS),
                                   vqrshrn_n_s32(sum4567, FILTER_BITS));
  return vqmovun_s16(sum_s16);
}

void vpx_convolve12_horiz_neon_i8mm(const uint8_t *src, ptrdiff_t src_stride,
                                    uint8_t *dst, ptrdiff_t dst_stride,
                                    const InterpKernel12 *filter, int x0_q4,
                                    int x_step_q4, int y0_q4, int y_step_q4,
                                    int w, int h) {
  // Scaling not supported by Neon implementation.
  if (x_step_q4 != 16) {
    vpx_convolve12_horiz_c(src, src_stride, dst, dst_stride, filter, x0_q4,
                           x_step_q4, y0_q4, y_step_q4, w, h);
    return;
  }

  assert(w == 32 || w == 16 || w == 8);
  assert(h == 32 || h == 16 || h == 8);

  // Split 12-tap filter into two 6-tap filters, masking the top two elements.
  // { 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0, 0 }
  const int8x8_t mask = vcreate_s8(0x0000ffffffffffff);
  const int8x8_t filter_0 = vand_s8(vmovn_s16(vld1q_s16(filter[x0_q4])), mask);
  const int8x8_t filter_1 =
      vext_s8(vmovn_s16(vld1q_s16(filter[x0_q4] + 4)), vdup_n_s8(0), 2);

  // Stagger each 6-tap filter to enable use of matrix multiply instructions.
  // { f0, f1, f2, f3, f4, f5,  0,  0,  0, f0, f1, f2, f3, f4, f5,  0 }
  const int8x16_t x_filter[2] = {
    vcombine_s8(filter_0, vext_s8(filter_0, filter_0, 7)),
    vcombine_s8(filter_1, vext_s8(filter_1, filter_1, 7))
  };

  const uint8x16x2_t permute_tbl = vld1q_u8_x2(kMatMulPermuteTbl);

  src -= MAX_FILTER_TAP / 2 - 1;

  do {
    const uint8_t *s = src;
    uint8_t *d = dst;
    int width = w;

    do {
      uint8x16_t s0[2], s1[2], s2[2], s3[2];
      load_u8_16x4(s, src_stride, &s0[0], &s1[0], &s2[0], &s3[0]);
      load_u8_16x4(s + 6, src_stride, &s0[1], &s1[1], &s2[1], &s3[1]);

      uint8x8_t d0 = convolve12_8_h(s0, x_filter, permute_tbl);
      uint8x8_t d1 = convolve12_8_h(s1, x_filter, permute_tbl);
      uint8x8_t d2 = convolve12_8_h(s2, x_filter, permute_tbl);
      uint8x8_t d3 = convolve12_8_h(s3, x_filter, permute_tbl);

      store_u8_8x4(d, dst_stride, d0, d1, d2, d3);

      s += 8;
      d += 8;
      width -= 8;
    } while (width != 0);
    src += 4 * src_stride;
    dst += 4 * dst_stride;
    h -= 4;
  } while (h != 0);
}
