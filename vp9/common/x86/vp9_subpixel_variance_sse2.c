/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#define HALFNDX 8

void vp9_half_horiz_variance16x_h_sse2(const unsigned char *ref_ptr,
                                       int ref_pixels_per_line,
                                       const unsigned char *src_ptr,
                                       int src_pixels_per_line,
                                       unsigned int Height,
                                       int *sum,
                                       unsigned int *sumsquared);

void vp9_half_vert_variance16x_h_sse2(const unsigned char *ref_ptr,
                                      int ref_pixels_per_line,
                                      const unsigned char *src_ptr,
                                      int src_pixels_per_line,
                                      unsigned int Height,
                                      int *sum,
                                      unsigned int *sumsquared);

void vp9_half_horiz_vert_variance16x_h_sse2(const unsigned char *ref_ptr,
                                            int ref_pixels_per_line,
                                            const unsigned char *src_ptr,
                                            int src_pixels_per_line,
                                            unsigned int Height,
                                            int *sum,
                                            unsigned int *sumsquared);

void vp9_filter_block2d_bil_var_sse2(const unsigned char *ref_ptr,
                                     int ref_pixels_per_line,
                                     const unsigned char *src_ptr,
                                     int src_pixels_per_line,
                                     unsigned int Height,
                                     int  xoffset,
                                     int  yoffset,
                                     int *sum,
                                     unsigned int *sumsquared);
