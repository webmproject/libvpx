/*
 *  Copyright (c) 2013 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */
#ifndef VPX_DSP_VPX_CONVOLVE_H_
#define VPX_DSP_VPX_CONVOLVE_H_

#include "./vpx_config.h"
#include "vpx/vpx_integer.h"

#ifdef __cplusplus
extern "C" {
#endif

// Note: Fixed size intermediate buffers, place limits on parameters
// of some functions. 2d filtering proceeds in 2 steps:
//   (1) Interpolate horizontally into an intermediate buffer, temp.
//   (2) Interpolate temp vertically to derive the sub-pixel result.
// Deriving the maximum number of rows in the temp buffer (135):
// --Smallest scaling factor is x1/2 ==> y_step_q4 = 32 (Normative).
// --Largest block size is 64x64 pixels.
// --64 rows in the downscaled frame span a distance of (64 - 1) * 32 in the
//   original frame (in 1/16th pixel units).
// --Must round-up because block may be located at sub-pixel position.
// --Require an additional SUBPEL_TAPS rows for the 8-tap filter tails.
// --((64 - 1) * 32 + 15) >> 4 + 8 = 135.
#if CONFIG_VP10 && CONFIG_EXT_PARTITION
# define MAX_EXT_SIZE 263
#else
# define MAX_EXT_SIZE 135
#endif  // CONFIG_VP10 && CONFIG_EXT_PARTITION

typedef void (*convolve_fn_t)(const uint8_t *src, ptrdiff_t src_stride,
                              uint8_t *dst, ptrdiff_t dst_stride,
                              const int16_t *filter_x, int x_step_q4,
                              const int16_t *filter_y, int y_step_q4,
                              int w, int h);

#if CONFIG_VPX_HIGHBITDEPTH
typedef void (*highbd_convolve_fn_t)(const uint8_t *src, ptrdiff_t src_stride,
                                     uint8_t *dst, ptrdiff_t dst_stride,
                                     const int16_t *filter_x, int x_step_q4,
                                     const int16_t *filter_y, int y_step_q4,
                                     int w, int h, int bd);
#endif

#ifdef __cplusplus
}  // extern "C"
#endif

#endif  // VPX_DSP_VPX_CONVOLVE_H_
