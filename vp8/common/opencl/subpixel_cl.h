/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */


#ifndef SUBPIXEL_CL_H
#define SUBPIXEL_CL_H

#include "../blockd.h"

/* Note:
 *
 * This platform is commonly built for runtime CPU detection. If you modify
 * any of the function mappings present in this file, be sure to also update
 * them in the function pointer initialization code
 */

#define prototype_subpixel_predict_cl(sym) \
    void sym(cl_command_queue cq, unsigned char *src_base, cl_mem src_mem, int src_offset, \
            int src_pitch, int xofst, int yofst, \
             unsigned char *dst_base, cl_mem dst_mem, int dst_offset, int dst_pitch)

extern prototype_subpixel_predict_cl(vp8_sixtap_predict16x16_cl);
extern prototype_subpixel_predict_cl(vp8_sixtap_predict8x8_cl);
extern prototype_subpixel_predict_cl(vp8_sixtap_predict8x4_cl);
extern prototype_subpixel_predict_cl(vp8_sixtap_predict4x4_cl);
extern prototype_subpixel_predict_cl(vp8_bilinear_predict16x16_cl);
extern prototype_subpixel_predict_cl(vp8_bilinear_predict8x8_cl);
extern prototype_subpixel_predict_cl(vp8_bilinear_predict8x4_cl);
extern prototype_subpixel_predict_cl(vp8_bilinear_predict4x4_cl);

typedef prototype_subpixel_predict_cl((*vp8_subpix_cl_fn_t));

//typedef enum
//{
//    SIXTAP = 0,
//    BILINEAR = 1
//} SUBPIX_TYPE;

#endif
