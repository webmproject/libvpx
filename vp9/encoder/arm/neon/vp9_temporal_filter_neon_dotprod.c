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

DECLARE_ALIGNED(16, static const uint8_t, kDotProdPermuteTbl[48]) = {
  // clang-format off
  0,  1,  2,  3,  1,  2,  3,  4,  2,  3,  4,  5,  3,  4,  5,  6,
  4,  5,  6,  7,  5,  6,  7,  8,  6,  7,  8,  9,  7,  8,  9, 10,
  8,  9, 10, 11,  9, 10, 11, 12, 10, 11, 12, 13, 11, 12, 13, 14
  // clang-format on
};

static INLINE uint8x8_t convolve12_8_h(uint8x16_t samples[2],
                                       const int8x16_t filter,
                                       const uint8x16x3_t perm_tbl) {
  // Transform sample range to [-128, 127] for 8-bit signed dot product.
  int8x16_t samples_128[2] = {
    vreinterpretq_s8_u8(vsubq_u8(samples[0], vdupq_n_u8(128))),
    vreinterpretq_s8_u8(vsubq_u8(samples[1], vdupq_n_u8(128)))
  };

  // Permute samples ready for dot product.
  // {  0,  1,  2,  3,  1,  2,  3,  4,  2,  3,  4,  5,  3,  4,  5,  6 }
  // {  4,  5,  6,  7,  5,  6,  7,  8,  6,  7,  8,  9,  7,  8,  9, 10 }
  // {  8,  9, 10, 11,  9, 10, 11, 12, 10, 11, 12, 13, 11, 12, 13, 14 }
  // { 12, 13, 14, 15, 13, 14, 15, 16, 14, 15, 16, 17, 15, 16, 17, 18 }
  int8x16_t perm_samples[4] = { vqtbl1q_s8(samples_128[0], perm_tbl.val[0]),
                                vqtbl1q_s8(samples_128[0], perm_tbl.val[1]),
                                vqtbl1q_s8(samples_128[0], perm_tbl.val[2]),
                                vqtbl1q_s8(samples_128[1], perm_tbl.val[2]) };

  // Accumulate into 128 << FILTER_BITS to account for range transform.
  int32x4_t acc = vdupq_n_s32(128 << FILTER_BITS);

  int32x4_t sum0123 = vdotq_laneq_s32(acc, perm_samples[0], filter, 0);
  sum0123 = vdotq_laneq_s32(sum0123, perm_samples[1], filter, 1);
  sum0123 = vdotq_laneq_s32(sum0123, perm_samples[2], filter, 2);

  int32x4_t sum4567 = vdotq_laneq_s32(acc, perm_samples[1], filter, 0);
  sum4567 = vdotq_laneq_s32(sum4567, perm_samples[2], filter, 1);
  sum4567 = vdotq_laneq_s32(sum4567, perm_samples[3], filter, 2);

  // Narrow and re-pack.
  int16x8_t sum_s16 = vcombine_s16(vqrshrn_n_s32(sum0123, FILTER_BITS),
                                   vqrshrn_n_s32(sum4567, FILTER_BITS));
  return vqmovun_s16(sum_s16);
}

void vpx_convolve12_horiz_neon_dotprod(const uint8_t *src, ptrdiff_t src_stride,
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

  const int16x8_t x_filter_0_7 = vld1q_s16(filter[x0_q4]);
  const int16x4_t x_filter_8_11 = vld1_s16(filter[x0_q4] + 8);
  const int16x8_t x_filter_8_15 = vcombine_s16(x_filter_8_11, vdup_n_s16(0));
  const int8x16_t x_filter =
      vcombine_s8(vmovn_s16(x_filter_0_7), vmovn_s16(x_filter_8_15));

  const uint8x16x3_t permute_tbl = vld1q_u8_x3(kDotProdPermuteTbl);

  src -= MAX_FILTER_TAP / 2 - 1;

  do {
    const uint8_t *s = src;
    uint8_t *d = dst;
    int width = w;

    do {
      uint8x16_t s0[2], s1[2], s2[2], s3[2];
      load_u8_16x4(s, src_stride, &s0[0], &s1[0], &s2[0], &s3[0]);
      load_u8_16x4(s + 4, src_stride, &s0[1], &s1[1], &s2[1], &s3[1]);

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
