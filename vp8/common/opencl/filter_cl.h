/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#ifndef FILTER_CL_H_
#define FILTER_CL_H_

#ifdef	__cplusplus
extern "C" {
#endif

#include "vp8_opencl.h"

#define VP8_FILTER_WEIGHT 128
#define VP8_FILTER_SHIFT  7

#define REGISTER_FILTER 1
#define CLAMP(x,min,max) if (x < min) x = min; else if ( x > max ) x = max;
#define PRE_CALC_PIXEL_STEPS 1
#define PRE_CALC_SRC_INCREMENT 1

#if PRE_CALC_PIXEL_STEPS
#define PS2 two_pixel_steps
#define PS3 three_pixel_steps
#else
#define PS2 2*(int)pixel_step
#define PS3 3*(int)pixel_step
#endif

#if REGISTER_FILTER
#define FILTER0 filter0
#define FILTER1 filter1
#define FILTER2 filter2
#define FILTER3 filter3
#define FILTER4 filter4
#define FILTER5 filter5
#else
#define FILTER0 vp8_filter[0]
#define FILTER1 vp8_filter[1]
#define FILTER2 vp8_filter[2]
#define FILTER3 vp8_filter[3]
#define FILTER4 vp8_filter[4]
#define FILTER5 vp8_filter[5]
#endif

#if PRE_CALC_SRC_INCREMENT
#define SRC_INCREMENT src_increment
#else
#define SRC_INCREMENT (src_pixels_per_line - output_width)
#endif

#define FILTER_OFFSET //Filter data stored as CL constant memory
#define FILTER_REF sub_pel_filters[filter_offset]

extern const char *filterCompileOptions;
extern const char *filter_cl_file_name;

//Copy the -2*pixel_step (and ps*3) bytes because the filter algorithm
//accesses negative indexes
#define SIXTAP_SRC_LEN(out_width,out_height,src_px) ((out_width)*(out_height) + (((out_width)*(out_height)-1)/(out_width))*(src_px - out_width) + 5)
#define BIL_SRC_LEN(out_width,out_height,src_px) ((out_height) * src_px + out_width)
#define DST_LEN(dst_pitch,dst_height,dst_width) (dst_pitch * (dst_height) + (dst_width))

#ifdef	__cplusplus
}
#endif

#endif /* FILTER_CL_H_ */
