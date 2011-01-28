/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include "vpx_config.h"
#include "variance.h"

#if HAVE_ARMV7

unsigned int vp8_sub_pixel_variance16x16_neon
(
    const unsigned char  *src_ptr,
    int  src_pixels_per_line,
    int  xoffset,
    int  yoffset,
    const unsigned char *dst_ptr,
    int dst_pixels_per_line,
    unsigned int *sse
)
{
  if (xoffset == 4 && yoffset == 0)
    return vp8_variance_halfpixvar16x16_h_neon(src_ptr, src_pixels_per_line, dst_ptr, dst_pixels_per_line, sse);
  else if (xoffset == 0 && yoffset == 4)
    return vp8_variance_halfpixvar16x16_v_neon(src_ptr, src_pixels_per_line, dst_ptr, dst_pixels_per_line, sse);
  else if (xoffset == 4 && yoffset == 4)
    return vp8_variance_halfpixvar16x16_hv_neon(src_ptr, src_pixels_per_line, dst_ptr, dst_pixels_per_line, sse);
  else
    return vp8_sub_pixel_variance16x16_neon_func(src_ptr, src_pixels_per_line, xoffset, yoffset, dst_ptr, dst_pixels_per_line, sse);
}

#endif
