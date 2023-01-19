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
#include "./vpx_dsp_rtcd.h"
#include "./vpx_config.h"

#include "vpx/vpx_integer.h"

#include "vpx_dsp/variance.h"
#include "vpx_dsp/arm/mem_neon.h"

// Process a block exactly 4 wide and a multiple of 2 high.
static void var_filter_block2d_bil_w4(const uint8_t *src_ptr, uint8_t *dst_ptr,
                                      int src_stride, int pixel_step,
                                      int dst_height, int filter_offset) {
  const uint8x8_t f0 = vdup_n_u8(8 - filter_offset);
  const uint8x8_t f1 = vdup_n_u8(filter_offset);

  int i = dst_height;
  do {
    uint8x8_t s0 = load_unaligned_u8(src_ptr, src_stride);
    uint8x8_t s1 = load_unaligned_u8(src_ptr + pixel_step, src_stride);
    uint16x8_t blend = vmlal_u8(vmull_u8(s0, f0), s1, f1);
    uint8x8_t blend_u8 = vrshrn_n_u16(blend, 3);
    vst1_u8(dst_ptr, blend_u8);

    src_ptr += 2 * src_stride;
    dst_ptr += 2 * 4;
    i -= 2;
  } while (i != 0);
}

// Process a block exactly 8 wide and any height.
static void var_filter_block2d_bil_w8(const uint8_t *src_ptr, uint8_t *dst_ptr,
                                      int src_stride, int pixel_step,
                                      int dst_height, int filter_offset) {
  const uint8x8_t f0 = vdup_n_u8(8 - filter_offset);
  const uint8x8_t f1 = vdup_n_u8(filter_offset);

  int i = dst_height;
  do {
    uint8x8_t s0 = vld1_u8(src_ptr);
    uint8x8_t s1 = vld1_u8(src_ptr + pixel_step);
    uint16x8_t blend = vmlal_u8(vmull_u8(s0, f0), s1, f1);
    uint8x8_t blend_u8 = vrshrn_n_u16(blend, 3);
    vst1_u8(dst_ptr, blend_u8);

    src_ptr += src_stride;
    dst_ptr += 8;
  } while (--i != 0);
}

// Process a block which is a mutiple of 16 wide and any height.
static void var_filter_block2d_bil_large(const uint8_t *src_ptr,
                                         uint8_t *dst_ptr, int src_stride,
                                         int pixel_step, int dst_width,
                                         int dst_height, int filter_offset) {
  const uint8x8_t f0 = vdup_n_u8(8 - filter_offset);
  const uint8x8_t f1 = vdup_n_u8(filter_offset);

  int i = dst_height;
  do {
    int j = 0;
    do {
      uint8x16_t s0 = vld1q_u8(src_ptr + j);
      uint8x16_t s1 = vld1q_u8(src_ptr + j + pixel_step);
      uint16x8_t blend_l =
          vmlal_u8(vmull_u8(vget_low_u8(s0), f0), vget_low_u8(s1), f1);
      uint16x8_t blend_h =
          vmlal_u8(vmull_u8(vget_high_u8(s0), f0), vget_high_u8(s1), f1);
      uint8x8_t out_lo = vrshrn_n_u16(blend_l, 3);
      uint8x8_t out_hi = vrshrn_n_u16(blend_h, 3);
      vst1q_u8(dst_ptr + j, vcombine_u8(out_lo, out_hi));

      j += 16;
    } while (j < dst_width);

    src_ptr += src_stride;
    dst_ptr += dst_width;
  } while (--i != 0);
}

static void var_filter_block2d_bil_w16(const uint8_t *src_ptr, uint8_t *dst_ptr,
                                       int src_stride, int pixel_step,
                                       int dst_height, int filter_offset) {
  var_filter_block2d_bil_large(src_ptr, dst_ptr, src_stride, pixel_step, 16,
                               dst_height, filter_offset);
}
static void var_filter_block2d_bil_w32(const uint8_t *src_ptr, uint8_t *dst_ptr,
                                       int src_stride, int pixel_step,
                                       int dst_height, int filter_offset) {
  var_filter_block2d_bil_large(src_ptr, dst_ptr, src_stride, pixel_step, 32,
                               dst_height, filter_offset);
}
static void var_filter_block2d_bil_w64(const uint8_t *src_ptr, uint8_t *dst_ptr,
                                       int src_stride, int pixel_step,
                                       int dst_height, int filter_offset) {
  var_filter_block2d_bil_large(src_ptr, dst_ptr, src_stride, pixel_step, 64,
                               dst_height, filter_offset);
}

#define SUBPEL_VARIANCE_WXH_NEON(w, h, padding)                          \
  unsigned int vpx_sub_pixel_variance##w##x##h##_neon(                   \
      const uint8_t *src, int src_stride, int xoffset, int yoffset,      \
      const uint8_t *ref, int ref_stride, uint32_t *sse) {               \
    uint8_t tmp0[w * (h + padding)];                                     \
    uint8_t tmp1[w * h];                                                 \
    var_filter_block2d_bil_w##w(src, tmp0, src_stride, 1, (h + padding), \
                                xoffset);                                \
    var_filter_block2d_bil_w##w(tmp0, tmp1, w, w, h, yoffset);           \
    return vpx_variance##w##x##h(tmp1, w, ref, ref_stride, sse);         \
  }

// 4x<h> blocks are processed two rows at a time, so require an extra row of
// padding.
SUBPEL_VARIANCE_WXH_NEON(4, 4, 2)
SUBPEL_VARIANCE_WXH_NEON(4, 8, 2)

SUBPEL_VARIANCE_WXH_NEON(8, 4, 1)
SUBPEL_VARIANCE_WXH_NEON(8, 8, 1)
SUBPEL_VARIANCE_WXH_NEON(8, 16, 1)

SUBPEL_VARIANCE_WXH_NEON(16, 8, 1)
SUBPEL_VARIANCE_WXH_NEON(16, 16, 1)
SUBPEL_VARIANCE_WXH_NEON(16, 32, 1)

SUBPEL_VARIANCE_WXH_NEON(32, 16, 1)
SUBPEL_VARIANCE_WXH_NEON(32, 32, 1)
SUBPEL_VARIANCE_WXH_NEON(32, 64, 1)

SUBPEL_VARIANCE_WXH_NEON(64, 32, 1)
SUBPEL_VARIANCE_WXH_NEON(64, 64, 1)

// 4xM filter writes an extra row to fdata because it processes two rows at a
// time.
#define SUB_PIXEL_AVG_VARIANCENXM(n, m)                                       \
  uint32_t vpx_sub_pixel_avg_variance##n##x##m##_neon(                        \
      const uint8_t *src_ptr, int src_stride, int x_offset, int y_offset,     \
      const uint8_t *ref_ptr, int ref_stride, uint32_t *sse,                  \
      const uint8_t *second_pred) {                                           \
    uint8_t temp0[n * (m + (n == 4 ? 2 : 1))];                                \
    uint8_t temp1[n * m];                                                     \
                                                                              \
    if (n == 4) {                                                             \
      var_filter_block2d_bil_w4(src_ptr, temp0, src_stride, 1, (m + 2),       \
                                x_offset);                                    \
      var_filter_block2d_bil_w4(temp0, temp1, n, n, m, y_offset);             \
    } else if (n == 8) {                                                      \
      var_filter_block2d_bil_w8(src_ptr, temp0, src_stride, 1, (m + 1),       \
                                x_offset);                                    \
      var_filter_block2d_bil_w8(temp0, temp1, n, n, m, y_offset);             \
    } else {                                                                  \
      var_filter_block2d_bil_large(src_ptr, temp0, src_stride, 1, n, (m + 1), \
                                   x_offset);                                 \
      var_filter_block2d_bil_large(temp0, temp1, n, n, n, m, y_offset);       \
    }                                                                         \
                                                                              \
    vpx_comp_avg_pred(temp0, second_pred, n, m, temp1, n);                    \
                                                                              \
    return vpx_variance##n##x##m(temp0, n, ref_ptr, ref_stride, sse);         \
  }

SUB_PIXEL_AVG_VARIANCENXM(4, 4)
SUB_PIXEL_AVG_VARIANCENXM(4, 8)
SUB_PIXEL_AVG_VARIANCENXM(8, 4)
SUB_PIXEL_AVG_VARIANCENXM(8, 8)
SUB_PIXEL_AVG_VARIANCENXM(8, 16)
SUB_PIXEL_AVG_VARIANCENXM(16, 8)
SUB_PIXEL_AVG_VARIANCENXM(16, 16)
SUB_PIXEL_AVG_VARIANCENXM(16, 32)
SUB_PIXEL_AVG_VARIANCENXM(32, 16)
SUB_PIXEL_AVG_VARIANCENXM(32, 32)
SUB_PIXEL_AVG_VARIANCENXM(32, 64)
SUB_PIXEL_AVG_VARIANCENXM(64, 32)
SUB_PIXEL_AVG_VARIANCENXM(64, 64)
