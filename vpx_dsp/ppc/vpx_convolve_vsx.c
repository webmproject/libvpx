/*
 *  Copyright (c) 2017 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */
#include <string.h>
#include "./vpx_dsp_rtcd.h"
#include "vpx_dsp/ppc/types_vsx.h"

// TODO(lu_zero): unroll
static inline void copy_w16(const uint8_t *src, ptrdiff_t src_stride,
                            uint8_t *dst, ptrdiff_t dst_stride, int32_t h) {
  int i;

  for (i = h; i--;) {
    vec_vsx_st(vec_vsx_ld(0, src), 0, dst);
    src += src_stride;
    dst += dst_stride;
  }
}

static inline void copy_w32(const uint8_t *src, ptrdiff_t src_stride,
                            uint8_t *dst, ptrdiff_t dst_stride, int32_t h) {
  int i;

  for (i = h; i--;) {
    vec_vsx_st(vec_vsx_ld(0, src), 0, dst);
    vec_vsx_st(vec_vsx_ld(16, src), 16, dst);
    src += src_stride;
    dst += dst_stride;
  }
}

static inline void copy_w64(const uint8_t *src, ptrdiff_t src_stride,
                            uint8_t *dst, ptrdiff_t dst_stride, int32_t h) {
  int i;

  for (i = h; i--;) {
    vec_vsx_st(vec_vsx_ld(0, src), 0, dst);
    vec_vsx_st(vec_vsx_ld(16, src), 16, dst);
    vec_vsx_st(vec_vsx_ld(32, src), 32, dst);
    vec_vsx_st(vec_vsx_ld(48, src), 48, dst);
    src += src_stride;
    dst += dst_stride;
  }
}

void vpx_convolve_copy_vsx(const uint8_t *src, ptrdiff_t src_stride,
                           uint8_t *dst, ptrdiff_t dst_stride,
                           const int16_t *filter_x, int32_t filter_x_stride,
                           const int16_t *filter_y, int32_t filter_y_stride,
                           int32_t w, int32_t h) {
  (void)filter_x;
  (void)filter_y;
  (void)filter_x_stride;
  (void)filter_y_stride;

  switch (w) {
    case 16: {
      copy_w16(src, src_stride, dst, dst_stride, h);
      break;
    }
    case 32: {
      copy_w32(src, src_stride, dst, dst_stride, h);
      break;
    }
    case 64: {
      copy_w64(src, src_stride, dst, dst_stride, h);
      break;
    }
    default: {
      int i;
      for (i = h; i--;) {
        memcpy(dst, src, w);
        src += src_stride;
        dst += dst_stride;
      }
      break;
    }
  }
}
