/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */


#include "vp9/common/vp9_filter.h"



static void variance(const unsigned char *src_ptr,
                     int  source_stride,
                     const unsigned char *ref_ptr,
                     int  recon_stride,
                     int  w,
                     int  h,
                     unsigned int *sse,
                     int *sum) {
  int i, j;
  int diff;

  *sum = 0;
  *sse = 0;

  for (i = 0; i < h; i++) {
    for (j = 0; j < w; j++) {
      diff = src_ptr[j] - ref_ptr[j];
      *sum += diff;
      *sse += diff * diff;
    }

    src_ptr += source_stride;
    ref_ptr += recon_stride;
  }
}

/****************************************************************************
 *
 *  ROUTINE       : filter_block2d_bil_first_pass
 *
 *  INPUTS        : UINT8  *src_ptr          : Pointer to source block.
 *                  UINT32 src_pixels_per_line : Stride of input block.
 *                  UINT32 pixel_step        : Offset between filter input samples (see notes).
 *                  UINT32 output_height     : Input block height.
 *                  UINT32 output_width      : Input block width.
 *                  INT32  *vp9_filter          : Array of 2 bi-linear filter taps.
 *
 *  OUTPUTS       : INT32 *output_ptr        : Pointer to filtered block.
 *
 *  RETURNS       : void
 *
 *  FUNCTION      : Applies a 1-D 2-tap bi-linear filter to the source block in
 *                  either horizontal or vertical direction to produce the
 *                  filtered output block. Used to implement first-pass
 *                  of 2-D separable filter.
 *
 *  SPECIAL NOTES : Produces INT32 output to retain precision for next pass.
 *                  Two filter taps should sum to VP9_FILTER_WEIGHT.
 *                  pixel_step defines whether the filter is applied
 *                  horizontally (pixel_step=1) or vertically (pixel_step=stride).
 *                  It defines the offset required to move from one input
 *                  to the next.
 *
 ****************************************************************************/
static void var_filter_block2d_bil_first_pass(const unsigned char *src_ptr,
                                              unsigned short *output_ptr,
                                              unsigned int src_pixels_per_line,
                                              int pixel_step,
                                              unsigned int output_height,
                                              unsigned int output_width,
                                              const short *vp9_filter) {
  unsigned int i, j;

  for (i = 0; i < output_height; i++) {
    for (j = 0; j < output_width; j++) {
      // Apply bilinear filter
      output_ptr[j] = (((int)src_ptr[0]          * vp9_filter[0]) +
                       ((int)src_ptr[pixel_step] * vp9_filter[1]) +
                       (VP9_FILTER_WEIGHT / 2)) >> VP9_FILTER_SHIFT;
      src_ptr++;
    }

    // Next row...
    src_ptr    += src_pixels_per_line - output_width;
    output_ptr += output_width;
  }
}

/****************************************************************************
 *
 *  ROUTINE       : filter_block2d_bil_second_pass
 *
 *  INPUTS        : INT32  *src_ptr          : Pointer to source block.
 *                  UINT32 src_pixels_per_line : Stride of input block.
 *                  UINT32 pixel_step        : Offset between filter input samples (see notes).
 *                  UINT32 output_height     : Input block height.
 *                  UINT32 output_width      : Input block width.
 *                  INT32  *vp9_filter          : Array of 2 bi-linear filter taps.
 *
 *  OUTPUTS       : UINT16 *output_ptr       : Pointer to filtered block.
 *
 *  RETURNS       : void
 *
 *  FUNCTION      : Applies a 1-D 2-tap bi-linear filter to the source block in
 *                  either horizontal or vertical direction to produce the
 *                  filtered output block. Used to implement second-pass
 *                  of 2-D separable filter.
 *
 *  SPECIAL NOTES : Requires 32-bit input as produced by filter_block2d_bil_first_pass.
 *                  Two filter taps should sum to VP9_FILTER_WEIGHT.
 *                  pixel_step defines whether the filter is applied
 *                  horizontally (pixel_step=1) or vertically (pixel_step=stride).
 *                  It defines the offset required to move from one input
 *                  to the next.
 *
 ****************************************************************************/
static void var_filter_block2d_bil_second_pass(const unsigned short *src_ptr,
                                               unsigned char *output_ptr,
                                               unsigned int src_pixels_per_line,
                                               unsigned int pixel_step,
                                               unsigned int output_height,
                                               unsigned int output_width,
                                               const short *vp9_filter) {
  unsigned int  i, j;
  int  Temp;

  for (i = 0; i < output_height; i++) {
    for (j = 0; j < output_width; j++) {
      // Apply filter
      Temp = ((int)src_ptr[0]         * vp9_filter[0]) +
             ((int)src_ptr[pixel_step] * vp9_filter[1]) +
             (VP9_FILTER_WEIGHT / 2);
      output_ptr[j] = (unsigned int)(Temp >> VP9_FILTER_SHIFT);
      src_ptr++;
    }

    // Next row...
    src_ptr    += src_pixels_per_line - output_width;
    output_ptr += output_width;
  }
}

