/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */


void  vp8cx_vertical_band_4_5_scale_c(unsigned char *dest, unsigned int dest_pitch, unsigned int dest_width);
void  vp8cx_last_vertical_band_4_5_scale_c(unsigned char *dest, unsigned int dest_pitch, unsigned int dest_width);
void  vp8cx_vertical_band_2_3_scale_c(unsigned char *dest, unsigned int dest_pitch, unsigned int dest_width);
void  vp8cx_last_vertical_band_2_3_scale_c(unsigned char *dest, unsigned int dest_pitch, unsigned int dest_width);
void  vp8cx_vertical_band_3_5_scale_c(unsigned char *dest, unsigned int dest_pitch, unsigned int dest_width);
void  vp8cx_last_vertical_band_3_5_scale_c(unsigned char *dest, unsigned int dest_pitch, unsigned int dest_width);
void  vp8cx_horizontal_line_1_2_scale_c(const unsigned char *source, unsigned int source_width, unsigned char *dest, unsigned int dest_width);
void  vp8cx_horizontal_line_3_5_scale_c(const unsigned char *source, unsigned int source_width, unsigned char *dest, unsigned int dest_width);
void  vp8cx_horizontal_line_2_3_scale_c(const unsigned char *source, unsigned int source_width, unsigned char *dest, unsigned int dest_width);
void  vp8cx_horizontal_line_4_5_scale_c(const unsigned char *source, unsigned int source_width, unsigned char *dest, unsigned int dest_width);
void  vp8cx_vertical_band_1_2_scale_c(unsigned char *dest, unsigned int dest_pitch, unsigned int dest_width);
void  vp8cx_last_vertical_band_1_2_scale_c(unsigned char *dest, unsigned int dest_pitch, unsigned int dest_width);

void  vp8cx_vertical_band_5_4_scale_c(unsigned char *source, unsigned int src_pitch, unsigned char *dest, unsigned int dest_pitch, unsigned int dest_width);
void  vp8cx_vertical_band_5_3_scale_c(unsigned char *source, unsigned int src_pitch, unsigned char *dest, unsigned int dest_pitch, unsigned int dest_width);
void  vp8cx_vertical_band_2_1_scale_c(unsigned char *source, unsigned int src_pitch, unsigned char *dest, unsigned int dest_pitch, unsigned int dest_width);
void  vp8cx_vertical_band_2_1_scale_i_c(unsigned char *source, unsigned int src_pitch, unsigned char *dest, unsigned int dest_pitch, unsigned int dest_width);
void  vp8cx_horizontal_line_2_1_scale_c(const unsigned char *source, unsigned int source_width, unsigned char *dest, unsigned int dest_width);
void  vp8cx_horizontal_line_5_3_scale_c(const unsigned char *source, unsigned int source_width, unsigned char *dest, unsigned int dest_width);
void  vp8cx_horizontal_line_5_4_scale_c(const unsigned char *source, unsigned int source_width, unsigned char *dest, unsigned int dest_width);

#define vp8_vertical_band_4_5_scale     vp8cx_vertical_band_4_5_scale_c
#define vp8_last_vertical_band_4_5_scale vp8cx_last_vertical_band_4_5_scale_c
#define vp8_vertical_band_2_3_scale     vp8cx_vertical_band_2_3_scale_c
#define vp8_last_vertical_band_2_3_scale vp8cx_last_vertical_band_2_3_scale_c
#define vp8_vertical_band_3_5_scale     vp8cx_vertical_band_3_5_scale_c
#define vp8_last_vertical_band_3_5_scale vp8cx_last_vertical_band_3_5_scale_c
#define vp8_horizontal_line_1_2_scale   vp8cx_horizontal_line_1_2_scale_c
#define vp8_horizontal_line_3_5_scale   vp8cx_horizontal_line_3_5_scale_c
#define vp8_horizontal_line_4_5_scale   vp8cx_horizontal_line_4_5_scale_c
#define vp8_horizontal_line_2_3_scale   vp8cx_horizontal_line_2_3_scale_c
#define vp8_vertical_band_1_2_scale     vp8cx_vertical_band_1_2_scale_c
#define vp8_last_vertical_band_1_2_scale vp8cx_last_vertical_band_1_2_scale_c
#define vp8_vertical_band_5_4_scale     vp8cx_vertical_band_5_4_scale_c
#define vp8_vertical_band_5_3_scale     vp8cx_vertical_band_5_3_scale_c
#define vp8_vertical_band_2_1_scale     vp8cx_vertical_band_2_1_scale_c
#define vp8_vertical_band_2_1_scale_i   vp8cx_vertical_band_2_1_scale_i_c
#define vp8_horizontal_line_2_1_scale   vp8cx_horizontal_line_2_1_scale_c
#define vp8_horizontal_line_5_3_scale   vp8cx_horizontal_line_5_3_scale_c
#define vp8_horizontal_line_5_4_scale   vp8cx_horizontal_line_5_4_scale_c
