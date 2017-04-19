/*
 *  Copyright (c) 2017 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include "./vp9_rtcd.h"
#include "./vpx_dsp_rtcd.h"
#include "./vpx_scale_rtcd.h"
#include "vp9/common/vp9_blockd.h"
#include "vpx_dsp/vpx_filter.h"
#include "vpx_scale/yv12config.h"

void vp9_scale_and_extend_frame_c(const YV12_BUFFER_CONFIG *src,
                                  YV12_BUFFER_CONFIG *dst, int phase_scaler) {
  const int src_w = src->y_crop_width;
  const int src_h = src->y_crop_height;
  const int dst_w = dst->y_crop_width;
  const int dst_h = dst->y_crop_height;
  const uint8_t *const srcs[3] = { src->y_buffer, src->u_buffer,
                                   src->v_buffer };
  const int src_strides[3] = { src->y_stride, src->uv_stride, src->uv_stride };
  uint8_t *const dsts[3] = { dst->y_buffer, dst->u_buffer, dst->v_buffer };
  const int dst_strides[3] = { dst->y_stride, dst->uv_stride, dst->uv_stride };
  const InterpKernel *const kernel = vp9_filter_kernels[EIGHTTAP];
  int x, y, i;

  for (i = 0; i < MAX_MB_PLANE; ++i) {
    const int factor = (i == 0 || i == 3 ? 1 : 2);
    const int src_stride = src_strides[i];
    const int dst_stride = dst_strides[i];
    for (y = 0; y < dst_h; y += 16) {
      const int y_q4 = y * (16 / factor) * src_h / dst_h + phase_scaler;
      for (x = 0; x < dst_w; x += 16) {
        const int x_q4 = x * (16 / factor) * src_w / dst_w + phase_scaler;
        const uint8_t *src_ptr = srcs[i] +
                                 (y / factor) * src_h / dst_h * src_stride +
                                 (x / factor) * src_w / dst_w;
        uint8_t *dst_ptr = dsts[i] + (y / factor) * dst_stride + (x / factor);

        vpx_scaled_2d(src_ptr, src_stride, dst_ptr, dst_stride,
                      kernel[x_q4 & 0xf], 16 * src_w / dst_w,
                      kernel[y_q4 & 0xf], 16 * src_h / dst_h, 16 / factor,
                      16 / factor);
      }
    }
  }

  vpx_extend_frame_borders(dst);
}
