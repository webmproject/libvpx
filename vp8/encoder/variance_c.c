/*
 *  Copyright (c) 2010 The VP8 project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license and patent
 *  grant that can be found in the LICENSE file in the root of the source
 *  tree. All contributing project authors may be found in the AUTHORS
 *  file in the root of the source tree.
 */


#include "variance.h"

const int vp8_six_tap[8][6] =
{
    { 0,  0,  128,    0,   0,  0 },         // note that 1/8 pel positions are just as per alpha -0.5 bicubic
    { 0, -6,  123,   12,  -1,  0 },
    { 2, -11, 108,   36,  -8,  1 },         // New 1/4 pel 6 tap filter
    { 0, -9,   93,   50,  -6,  0 },
    { 3, -16,  77,   77, -16,  3 },         // New 1/2 pel 6 tap filter
    { 0, -6,   50,   93,  -9,  0 },
    { 1, -8,   36,  108, -11,  2 },         // New 1/4 pel 6 tap filter
    { 0, -1,   12,  123,  -6,  0 }
};


#ifdef USEBILINEAR
const int VP8_FILTER_WEIGHT = 128;
const int VP8_FILTER_SHIFT  =   7;
const int vp8_bilinear_taps[8][2] =
{
    { 128,   0 },
    { 112,  16 },
    {  96,  32 },
    {  80,  48 },
    {  64,  64 },
    {  48,  80 },
    {  32,  96 },
    {  16, 112 }
};

unsigned int vp8_get_mb_ss_c
(
    short *src_ptr
)
{
    unsigned int i = 0, sum = 0;

    do
    {
        sum += (src_ptr[i] * src_ptr[i]);
        i++;
    }
    while (i < 256);

    return sum;
}


void  vp8_variance(
    unsigned char *src_ptr,
    int  source_stride,
    unsigned char *ref_ptr,
    int  recon_stride,
    int  w,
    int  h,
    unsigned int *sse,
    int *sum)
{
    int i, j;
    int diff;

    *sum = 0;
    *sse = 0;

    for (i = 0; i < h; i++)
    {
        for (j = 0; j < w; j++)
        {
            diff = src_ptr[j] - ref_ptr[j];
            *sum += diff;
            *sse += diff * diff;
        }

        src_ptr += source_stride;
        ref_ptr += recon_stride;
    }
}

unsigned int
vp8_get8x8var_c
(
    unsigned char *src_ptr,
    int  source_stride,
    unsigned char *ref_ptr,
    int  recon_stride,
    unsigned int *SSE,
    int *Sum
)
{

    vp8_variance(src_ptr, source_stride, ref_ptr, recon_stride, 8, 8, SSE, Sum);
    return (*SSE - (((*Sum) * (*Sum)) >> 6));
}

unsigned int
vp8_get16x16var_c
(
    unsigned char *src_ptr,
    int  source_stride,
    unsigned char *ref_ptr,
    int  recon_stride,
    unsigned int *SSE,
    int *Sum
)
{

    vp8_variance(src_ptr, source_stride, ref_ptr, recon_stride, 16, 16, SSE, Sum);
    return (*SSE - (((*Sum) * (*Sum)) >> 8));

}



unsigned int vp8_variance16x16_c(
    unsigned char *src_ptr,
    int  source_stride,
    unsigned char *ref_ptr,
    int  recon_stride,
    unsigned int *sse)
{
    unsigned int var;
    int avg;


    vp8_variance(src_ptr, source_stride, ref_ptr, recon_stride, 16, 16, &var, &avg);
    *sse = var;
    return (var - ((avg * avg) >> 8));
}

unsigned int vp8_variance8x16_c(
    unsigned char *src_ptr,
    int  source_stride,
    unsigned char *ref_ptr,
    int  recon_stride,
    unsigned int *sse)
{
    unsigned int var;
    int avg;


    vp8_variance(src_ptr, source_stride, ref_ptr, recon_stride, 8, 16, &var, &avg);
    *sse = var;
    return (var - ((avg * avg) >> 7));
}

unsigned int vp8_variance16x8_c(
    unsigned char *src_ptr,
    int  source_stride,
    unsigned char *ref_ptr,
    int  recon_stride,
    unsigned int *sse)
{
    unsigned int var;
    int avg;


    vp8_variance(src_ptr, source_stride, ref_ptr, recon_stride, 16, 8, &var, &avg);
    *sse = var;
    return (var - ((avg * avg) >> 7));
}


unsigned int vp8_variance8x8_c(
    unsigned char *src_ptr,
    int  source_stride,
    unsigned char *ref_ptr,
    int  recon_stride,
    unsigned int *sse)
{
    unsigned int var;
    int avg;


    vp8_variance(src_ptr, source_stride, ref_ptr, recon_stride, 8, 8, &var, &avg);
    *sse = var;
    return (var - ((avg * avg) >> 6));
}

unsigned int vp8_variance4x4_c(
    unsigned char *src_ptr,
    int  source_stride,
    unsigned char *ref_ptr,
    int  recon_stride,
    unsigned int *sse)
{
    unsigned int var;
    int avg;


    vp8_variance(src_ptr, source_stride, ref_ptr, recon_stride, 4, 4, &var, &avg);
    *sse = var;
    return (var - ((avg * avg) >> 4));
}


unsigned int vp8_mse16x16_c(
    unsigned char *src_ptr,
    int  source_stride,
    unsigned char *ref_ptr,
    int  recon_stride,
    unsigned int *sse)
{
    unsigned int var;
    int avg;

    vp8_variance(src_ptr, source_stride, ref_ptr, recon_stride, 16, 16, &var, &avg);
    *sse = var;
    return var;
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
 *                  INT32  *vp8_filter          : Array of 2 bi-linear filter taps.
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
 *                  Two filter taps should sum to VP8_FILTER_WEIGHT.
 *                  pixel_step defines whether the filter is applied
 *                  horizontally (pixel_step=1) or vertically (pixel_step=stride).
 *                  It defines the offset required to move from one input
 *                  to the next.
 *
 ****************************************************************************/
void vp8e_filter_block2d_bil_first_pass
(
    unsigned char *src_ptr,
    unsigned short *output_ptr,
    unsigned int src_pixels_per_line,
    int pixel_step,
    unsigned int output_height,
    unsigned int output_width,
    const int *vp8_filter
)
{
    unsigned int i, j;

    for (i = 0; i < output_height; i++)
    {
        for (j = 0; j < output_width; j++)
        {
            // Apply bilinear filter
            output_ptr[j] = (((int)src_ptr[0]          * vp8_filter[0]) +
                             ((int)src_ptr[pixel_step] * vp8_filter[1]) +
                             (VP8_FILTER_WEIGHT / 2)) >> VP8_FILTER_SHIFT;
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
 *                  INT32  *vp8_filter          : Array of 2 bi-linear filter taps.
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
 *                  Two filter taps should sum to VP8_FILTER_WEIGHT.
 *                  pixel_step defines whether the filter is applied
 *                  horizontally (pixel_step=1) or vertically (pixel_step=stride).
 *                  It defines the offset required to move from one input
 *                  to the next.
 *
 ****************************************************************************/
void vp8e_filter_block2d_bil_second_pass
(
    unsigned short *src_ptr,
    unsigned char  *output_ptr,
    unsigned int  src_pixels_per_line,
    unsigned int  pixel_step,
    unsigned int  output_height,
    unsigned int  output_width,
    const int *vp8_filter
)
{
    unsigned int  i, j;
    int  Temp;

    for (i = 0; i < output_height; i++)
    {
        for (j = 0; j < output_width; j++)
        {
            // Apply filter
            Temp = ((int)src_ptr[0]         * vp8_filter[0]) +
                   ((int)src_ptr[pixel_step] * vp8_filter[1]) +
                   (VP8_FILTER_WEIGHT / 2);
            output_ptr[j] = (unsigned int)(Temp >> VP8_FILTER_SHIFT);
            src_ptr++;
        }

        // Next row...
        src_ptr    += src_pixels_per_line - output_width;
        output_ptr += output_width;
    }
}


/****************************************************************************
 *
 *  ROUTINE       : filter_block2d_bil
 *
 *  INPUTS        : UINT8  *src_ptr          : Pointer to source block.
 *                  UINT32 src_pixels_per_line : Stride of input block.
 *                  INT32  *HFilter         : Array of 2 horizontal filter taps.
 *                  INT32  *VFilter         : Array of 2 vertical filter taps.
 *
 *  OUTPUTS       : UINT16 *output_ptr       : Pointer to filtered block.
 *
 *  RETURNS       : void
 *
 *  FUNCTION      : 2-D filters an 8x8 input block by applying a 2-tap
 *                  bi-linear filter horizontally followed by a 2-tap
 *                  bi-linear filter vertically on the result.
 *
 *  SPECIAL NOTES : The intermediate horizontally filtered block must produce
 *                  1 more point than the input block in each column. This
 *                  is to ensure that the 2-tap filter has one extra data-point
 *                  at the top of each column so filter taps do not extend
 *                  beyond data. Thus the output of the first stage filter
 *                  is an 8x9 (hx_v) block.
 *
 ****************************************************************************/
void vp8e_filter_block2d_bil
(
    unsigned char  *src_ptr,
    unsigned char *output_ptr,
    unsigned int src_pixels_per_line,
    int  *HFilter,
    int  *VFilter
)
{

    unsigned short FData[20*16];    // Temp data bufffer used in filtering

    // First filter 1-D horizontally...
    vp8e_filter_block2d_bil_first_pass(src_ptr, FData, src_pixels_per_line, 1, 9, 8, HFilter);

    // then 1-D vertically...
    vp8e_filter_block2d_bil_second_pass(FData, output_ptr, 8, 8, 8, 8, VFilter);
}



unsigned int vp8_sub_pixel_variance4x4_c
(
    unsigned char  *src_ptr,
    int  src_pixels_per_line,
    int  xoffset,
    int  yoffset,
    unsigned char *dst_ptr,
    int dst_pixels_per_line,
    unsigned int *sse
)
{
    unsigned char  temp2[20*16];
    const int *HFilter, *VFilter;
    unsigned short FData3[5*4]; // Temp data bufffer used in filtering

    HFilter = vp8_bilinear_taps[xoffset];
    VFilter = vp8_bilinear_taps[yoffset];

    // First filter 1d Horizontal
    vp8e_filter_block2d_bil_first_pass(src_ptr, FData3, src_pixels_per_line, 1, 5, 4, HFilter);

    // Now filter Verticaly
    vp8e_filter_block2d_bil_second_pass(FData3, temp2, 4,  4,  4,  4, VFilter);

    return vp8_variance4x4_c(temp2, 4, dst_ptr, dst_pixels_per_line, sse);
}


unsigned int vp8_sub_pixel_variance8x8_c
(
    unsigned char  *src_ptr,
    int  src_pixels_per_line,
    int  xoffset,
    int  yoffset,
    unsigned char *dst_ptr,
    int dst_pixels_per_line,
    unsigned int *sse
)
{
    unsigned short FData3[9*8]; // Temp data bufffer used in filtering
    unsigned char  temp2[20*16];
    const int *HFilter, *VFilter;

    HFilter = vp8_bilinear_taps[xoffset];
    VFilter = vp8_bilinear_taps[yoffset];

    vp8e_filter_block2d_bil_first_pass(src_ptr, FData3, src_pixels_per_line, 1, 9, 8, HFilter);
    vp8e_filter_block2d_bil_second_pass(FData3, temp2, 8, 8, 8, 8, VFilter);

    return vp8_variance8x8_c(temp2, 8, dst_ptr, dst_pixels_per_line, sse);
}

unsigned int vp8_sub_pixel_variance16x16_c
(
    unsigned char  *src_ptr,
    int  src_pixels_per_line,
    int  xoffset,
    int  yoffset,
    unsigned char *dst_ptr,
    int dst_pixels_per_line,
    unsigned int *sse
)
{
    unsigned short FData3[17*16];   // Temp data bufffer used in filtering
    unsigned char  temp2[20*16];
    const int *HFilter, *VFilter;

    HFilter = vp8_bilinear_taps[xoffset];
    VFilter = vp8_bilinear_taps[yoffset];

    vp8e_filter_block2d_bil_first_pass(src_ptr, FData3, src_pixels_per_line, 1, 17, 16, HFilter);
    vp8e_filter_block2d_bil_second_pass(FData3, temp2, 16, 16, 16, 16, VFilter);

    return vp8_variance16x16_c(temp2, 16, dst_ptr, dst_pixels_per_line, sse);
}

unsigned int vp8_sub_pixel_mse16x16_c
(
    unsigned char  *src_ptr,
    int  src_pixels_per_line,
    int  xoffset,
    int  yoffset,
    unsigned char *dst_ptr,
    int dst_pixels_per_line,
    unsigned int *sse
)
{
    vp8_sub_pixel_variance16x16_c(src_ptr, src_pixels_per_line, xoffset, yoffset, dst_ptr, dst_pixels_per_line, sse);
    return *sse;
}

unsigned int vp8_sub_pixel_variance16x8_c
(
    unsigned char  *src_ptr,
    int  src_pixels_per_line,
    int  xoffset,
    int  yoffset,
    unsigned char *dst_ptr,
    int dst_pixels_per_line,
    unsigned int *sse
)
{
    unsigned short FData3[16*9];    // Temp data bufffer used in filtering
    unsigned char  temp2[20*16];
    const int *HFilter, *VFilter;

    HFilter = vp8_bilinear_taps[xoffset];
    VFilter = vp8_bilinear_taps[yoffset];

    vp8e_filter_block2d_bil_first_pass(src_ptr, FData3, src_pixels_per_line, 1, 9, 16, HFilter);
    vp8e_filter_block2d_bil_second_pass(FData3, temp2, 16, 16, 8, 16, VFilter);

    return vp8_variance16x8_c(temp2, 16, dst_ptr, dst_pixels_per_line, sse);
}

unsigned int vp8_sub_pixel_variance8x16_c
(
    unsigned char  *src_ptr,
    int  src_pixels_per_line,
    int  xoffset,
    int  yoffset,
    unsigned char *dst_ptr,
    int dst_pixels_per_line,
    unsigned int *sse
)
{
    unsigned short FData3[9*16];    // Temp data bufffer used in filtering
    unsigned char  temp2[20*16];
    const int *HFilter, *VFilter;


    HFilter = vp8_bilinear_taps[xoffset];
    VFilter = vp8_bilinear_taps[yoffset];


    vp8e_filter_block2d_bil_first_pass(src_ptr, FData3, src_pixels_per_line, 1, 17, 8, HFilter);
    vp8e_filter_block2d_bil_second_pass(FData3, temp2, 8, 8, 16, 8, VFilter);

    return vp8_variance8x16_c(temp2, 8, dst_ptr, dst_pixels_per_line, sse);
}
#endif
