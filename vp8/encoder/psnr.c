/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */


#include "vpx_scale/yv12config.h"
#include "math.h"
#include "vp8/common/systemdependent.h" /* for vp8_clear_system_state() */

#define MAX_PSNR 60

double vp8_mse2psnr(double Samples, double Peak, double Mse)
{
    double psnr;

    if ((double)Mse > 0.0)
        psnr = 10.0 * log10(Peak * Peak * Samples / Mse);
    else
        psnr = MAX_PSNR;      // Limit to prevent / 0

    if (psnr > MAX_PSNR)
        psnr = MAX_PSNR;

    return psnr;
}

double vp8_calc_psnr(YV12_BUFFER_CONFIG *source, YV12_BUFFER_CONFIG *dest, double *YPsnr, double *UPsnr, double *VPsnr, double *sq_error)
{
    int i, j;
    int Diff;
    double frame_psnr;
    double Total;
    double grand_total;
    unsigned char *src = source->y_buffer;
    unsigned char *dst = dest->y_buffer;

    Total = 0.0;
    grand_total = 0.0;

    // Loop throught the Y plane raw and reconstruction data summing (square differences)
    for (i = 0; i < source->y_height; i++)
    {

        for (j = 0; j < source->y_width; j++)
        {
            Diff        = (int)(src[j]) - (int)(dst[j]);
            Total      += Diff * Diff;
        }

        src += source->y_stride;
        dst += dest->y_stride;
    }

    // Work out Y PSNR
    *YPsnr = vp8_mse2psnr(source->y_height * source->y_width, 255.0, Total);
    grand_total += Total;
    Total = 0;


    // Loop through the U plane
    src = source->u_buffer;
    dst = dest->u_buffer;

    for (i = 0; i < source->uv_height; i++)
    {

        for (j = 0; j < source->uv_width; j++)
        {
            Diff        = (int)(src[j]) - (int)(dst[j]);
            Total      += Diff * Diff;
        }

        src += source->uv_stride;
        dst += dest->uv_stride;
    }

    // Work out U PSNR
    *UPsnr = vp8_mse2psnr(source->uv_height * source->uv_width, 255.0, Total);
    grand_total += Total;
    Total = 0;


    // V PSNR
    src = source->v_buffer;
    dst = dest->v_buffer;

    for (i = 0; i < source->uv_height; i++)
    {

        for (j = 0; j < source->uv_width; j++)
        {
            Diff        = (int)(src[j]) - (int)(dst[j]);
            Total      += Diff * Diff;
        }

        src += source->uv_stride;
        dst += dest->uv_stride;
    }

    // Work out UV PSNR
    *VPsnr = vp8_mse2psnr(source->uv_height * source->uv_width, 255.0, Total);
    grand_total += Total;
    Total = 0;

    // Work out total PSNR
    frame_psnr = vp8_mse2psnr(source->y_height * source->y_width * 3 / 2 , 255.0, grand_total);

    *sq_error = 1.0 * grand_total;

    return frame_psnr;
}
