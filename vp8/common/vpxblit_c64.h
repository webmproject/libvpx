/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */


#ifndef _VPX_BLIT_C64_h
#define _VPX_BLIT_C64_h

/****************************************************************************
*  Typedefs
****************************************************************************/

typedef struct  // YUV buffer configuration structure
{
    int     y_width;
    int     y_height;
    int     y_stride;

    int     uv_width;
    int     uv_height;
    int     uv_stride;

    unsigned char *y_buffer;
    unsigned char *u_buffer;
    unsigned char *v_buffer;

    unsigned char *y_ptr_scrn;
    unsigned char *u_ptr_scrn;
    unsigned char *v_ptr_scrn;

} DXV_YUV_BUFFER_CONFIG;

typedef struct
{
    unsigned char *rgbptr_scrn;
    unsigned char *y_ptr_scrn;
    unsigned char *u_ptr_scrn;
    unsigned char *v_ptr_scrn;
    unsigned char *rgbptr_scrn2;
} DXV_FINAL_VIDEO;

#endif /* include guards */
