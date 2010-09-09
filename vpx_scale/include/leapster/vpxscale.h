/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */


/****************************************************************************
*
*   Module Title :     postp.h
*
*   Description  :     Post processor interface
*
****************************************************************************/
#ifndef VPXSCALE_H
#define VPXSCALE_H


// fwg 2004-10-14
typedef void (*vpxvertical_band_4_5_scale_lf)(unsigned char *dest, unsigned int dest_pitch, unsigned int dest_width);
typedef void (*vpxlast_vertical_band_4_5_scale_lf)(unsigned char *dest, unsigned int dest_pitch, unsigned int dest_width);
typedef void (*vpxvertical_band_3_5_scale_lf)(unsigned char *dest, unsigned int dest_pitch, unsigned int dest_width);
typedef void (*vpxlast_vertical_band_3_5_scale_lf)(unsigned char *dest, unsigned int dest_pitch, unsigned int dest_width);
typedef void (*vpxhorizontal_line_1_2_scale_lf)(const unsigned char *source, unsigned int source_width, unsigned char *dest, unsigned int dest_width);
typedef void (*vpxhorizontal_line_3_5_scale_lf)(const unsigned char *source, unsigned int source_width, unsigned char *dest, unsigned int dest_width);
typedef void (*vpxhorizontal_line_4_5_scale_lf)(const unsigned char *source, unsigned int source_width, unsigned char *dest, unsigned int dest_width);
typedef void (*vpxvertical_band_1_2_scale_lf)(unsigned char *dest, unsigned int dest_pitch, unsigned int dest_width);
typedef void (*vpxlast_vertical_band_1_2_scale_lf)(unsigned char *dest, unsigned int dest_pitch, unsigned int dest_width);


typedef struct vpxglobal_scalling_ptrs_t
{
    vpxvertical_band_4_5_scale_lf        vpxvertical_band_4_5_scale_t;
    vpxlast_vertical_band_4_5_scale_lf    vpxlast_vertical_band_4_5_scale_t;
    vpxvertical_band_3_5_scale_lf        vpxvertical_band_3_5_scale_t;
    vpxlast_vertical_band_3_5_scale_lf    vpxlast_vertical_band_3_5_scale_t;
    vpxhorizontal_line_1_2_scale_lf      vpxhorizontal_line_1_2_scale_t;
    vpxhorizontal_line_3_5_scale_lf      vpxhorizontal_line_3_5_scale_t;
    vpxhorizontal_line_4_5_scale_lf      vpxhorizontal_line_4_5_scale_t;
    vpxvertical_band_1_2_scale_lf        vpxvertical_band_1_2_scale_t;
    vpxlast_vertical_band_1_2_scale_lf    vpxlast_vertical_band_1_2_scale_t;
} vpxglobal_scalling_ptrs;

extern struct vpxglobal_scalling_ptrs_t *g_scaling_ptrs;

/*
extern void  (*vp8_vertical_band_4_5_scale)(unsigned char * dest,unsigned int dest_pitch,unsigned int dest_width);
extern void  (*vp8_last_vertical_band_4_5_scale)(unsigned char * dest,unsigned int dest_pitch,unsigned int dest_width);
extern void  (*vp8_vertical_band_3_5_scale)(unsigned char * dest,unsigned int dest_pitch,unsigned int dest_width);
extern void  (*vp8_last_vertical_band_3_5_scale)(unsigned char * dest,unsigned int dest_pitch,unsigned int dest_width);
extern void  (*vp8_horizontal_line_1_2_scale)(const unsigned char * source,unsigned int source_width,unsigned char * dest,unsigned int dest_width);
extern void  (*vp8_horizontal_line_3_5_scale)(const unsigned char * source,unsigned int source_width,unsigned char * dest,unsigned int dest_width);
extern void  (*vp8_horizontal_line_4_5_scale)(const unsigned char * source,unsigned int source_width,unsigned char * dest,unsigned int dest_width);
extern void  (*vp8_vertical_band_1_2_scale)(unsigned char * dest,unsigned int dest_pitch,unsigned int dest_width);
extern void  (*vp8_last_vertical_band_1_2_scale)(unsigned char * dest,unsigned int dest_pitch,unsigned int dest_width);
*/

#endif
