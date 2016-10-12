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

#include "./vpx_dsp_rtcd.h"
#include "vpx/vpx_integer.h"

void vpx_highbd_convolve_copy_neon(const uint8_t *src8, ptrdiff_t src_stride,
                                   uint8_t *dst8, ptrdiff_t dst_stride,
                                   const int16_t *filter_x, int filter_x_stride,
                                   const int16_t *filter_y, int filter_y_stride,
                                   int w, int h, int bd) {
  const uint16_t *src = CONVERT_TO_SHORTPTR(src8);
  uint16_t *dst = CONVERT_TO_SHORTPTR(dst8);

  (void)filter_x;
  (void)filter_x_stride;
  (void)filter_y;
  (void)filter_y_stride;
  (void)bd;

  if (w < 8) {  // copy4
    do {
      vst1_u16(dst, vld1_u16(src));
      src += src_stride;
      dst += dst_stride;
      vst1_u16(dst, vld1_u16(src));
      src += src_stride;
      dst += dst_stride;
      h -= 2;
    } while (h > 0);
  } else if (w == 8) {  // copy8
    do {
      vst1q_u16(dst, vld1q_u16(src));
      src += src_stride;
      dst += dst_stride;
      vst1q_u16(dst, vld1q_u16(src));
      src += src_stride;
      dst += dst_stride;
      h -= 2;
    } while (h > 0);
  } else if (w < 32) {  // copy16
    do {
      vst1q_u16(dst, vld1q_u16(src));
      vst1q_u16(dst + 8, vld1q_u16(src + 8));
      src += src_stride;
      dst += dst_stride;
      vst1q_u16(dst, vld1q_u16(src));
      vst1q_u16(dst + 8, vld1q_u16(src + 8));
      src += src_stride;
      dst += dst_stride;
      h -= 2;
    } while (h > 0);
  } else if (w == 32) {  // copy32
    do {
      vst1q_u16(dst, vld1q_u16(src));
      vst1q_u16(dst + 8, vld1q_u16(src + 8));
      vst1q_u16(dst + 16, vld1q_u16(src + 16));
      vst1q_u16(dst + 24, vld1q_u16(src + 24));
      src += src_stride;
      dst += dst_stride;
      vst1q_u16(dst, vld1q_u16(src));
      vst1q_u16(dst + 8, vld1q_u16(src + 8));
      vst1q_u16(dst + 16, vld1q_u16(src + 16));
      vst1q_u16(dst + 24, vld1q_u16(src + 24));
      src += src_stride;
      dst += dst_stride;
      h -= 2;
    } while (h > 0);
  } else {  // copy64
    do {
      vst1q_u16(dst, vld1q_u16(src));
      vst1q_u16(dst + 8, vld1q_u16(src + 8));
      vst1q_u16(dst + 16, vld1q_u16(src + 16));
      vst1q_u16(dst + 24, vld1q_u16(src + 24));
      vst1q_u16(dst + 32, vld1q_u16(src + 32));
      vst1q_u16(dst + 40, vld1q_u16(src + 40));
      vst1q_u16(dst + 48, vld1q_u16(src + 48));
      vst1q_u16(dst + 56, vld1q_u16(src + 56));
      src += src_stride;
      dst += dst_stride;
    } while (--h);
  }
}
