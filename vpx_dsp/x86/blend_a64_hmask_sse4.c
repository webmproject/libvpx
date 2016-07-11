/*
 *  Copyright (c) 2016 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include "vpx/vpx_integer.h"

#include "./vpx_dsp_rtcd.h"

// To start out, just dispatch to the function using the 2D mask and
// pass mask stride as 0. This can be improved upon if necessary.

void vpx_blend_a64_hmask_sse4_1(
    uint8_t *dst, uint32_t dst_stride,
    const uint8_t *src0, uint32_t src0_stride,
    const uint8_t *src1, uint32_t src1_stride,
    const uint8_t *mask, int h, int w) {
  vpx_blend_a64_mask_sse4_1(dst, dst_stride,
                            src0, src0_stride,
                            src1, src1_stride,
                            mask, 0, h, w, 0, 0);
}

#if CONFIG_VP9_HIGHBITDEPTH
void vpx_highbd_blend_a64_hmask_sse4_1(
    uint8_t *dst_8, uint32_t dst_stride,
    const uint8_t *src0_8, uint32_t src0_stride,
    const uint8_t *src1_8, uint32_t src1_stride,
    const uint8_t *mask, int h, int w,
    int bd) {
  vpx_highbd_blend_a64_mask_sse4_1(dst_8, dst_stride,
                                   src0_8, src0_stride,
                                   src1_8, src1_stride,
                                   mask, 0, h, w, 0, 0, bd);
}
#endif  // CONFIG_VP9_HIGHBITDEPTH
