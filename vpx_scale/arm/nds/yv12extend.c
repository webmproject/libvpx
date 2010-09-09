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
*   Module Title :     yv12extend.c
*
*   Description  :
*
***************************************************************************/

/****************************************************************************
*  Header Files
****************************************************************************/
#include "vpx_scale/yv12config.h"
#include "vpx_mem/vpx_mem.h"
#include <nitro.h>
#include <nitro/mi.h>
#include <nitro/itcm_begin.h>

//---- DMA Number
#define DMA_NO  3

/****************************************************************************
*  Exports
****************************************************************************/

/****************************************************************************
*
****************************************************************************/
void
vp8_yv12_extend_frame_borders(YV12_BUFFER_CONFIG *ybf)
{
    int i;
    unsigned char *src_ptr1, *src_ptr2;
    unsigned char *dest_ptr1, *dest_ptr2;

    unsigned int Border;
    int plane_stride;
    int plane_height;
    int plane_width;

    /***********/
    /* Y Plane */
    /***********/
    Border = ybf->border;
    plane_stride = ybf->y_stride;
    plane_height = ybf->y_height;
    plane_width = ybf->y_width;

    // copy the left and right most columns out
    src_ptr1 = ybf->y_buffer;
    src_ptr2 = src_ptr1 + plane_width - 1;
    dest_ptr1 = src_ptr1 - Border;
    dest_ptr2 = src_ptr2 + 1;

    for (i = 0; i < plane_height; i++)
    {
        mi_cpu_fill8(dest_ptr1, src_ptr1[0], Border);
        mi_cpu_fill8(dest_ptr2, src_ptr2[0], Border);
        src_ptr1  += plane_stride;
        src_ptr2  += plane_stride;
        dest_ptr1 += plane_stride;
        dest_ptr2 += plane_stride;
    }

    // Now copy the top and bottom source lines into each line of the respective borders
    src_ptr1 = ybf->y_buffer - Border;
    src_ptr2 = src_ptr1 + (plane_height * plane_stride) - plane_stride;
    dest_ptr1 = src_ptr1 - (Border * plane_stride);
    dest_ptr2 = src_ptr2 + plane_stride;

    for (i = 0; i < (int)Border; i++)
    {
        mi_cpu_copy_fast(src_ptr1, dest_ptr1, plane_stride);
        mi_cpu_copy_fast(src_ptr2, dest_ptr2, plane_stride);
        dest_ptr1 += plane_stride;
        dest_ptr2 += plane_stride;
    }

    plane_stride /= 2;
    plane_height /= 2;
    plane_width /= 2;
    Border /= 2;

    /***********/
    /* U Plane */
    /***********/

    // copy the left and right most columns out
    src_ptr1 = ybf->u_buffer;
    src_ptr2 = src_ptr1 + plane_width - 1;
    dest_ptr1 = src_ptr1 - Border;
    dest_ptr2 = src_ptr2 + 1;

    for (i = 0; i < plane_height; i++)
    {
        mi_cpu_fill8(dest_ptr1, src_ptr1[0], Border);
        mi_cpu_fill8(dest_ptr2, src_ptr2[0], Border);
        src_ptr1  += plane_stride;
        src_ptr2  += plane_stride;
        dest_ptr1 += plane_stride;
        dest_ptr2 += plane_stride;
    }

    // Now copy the top and bottom source lines into each line of the respective borders
    src_ptr1 = ybf->u_buffer - Border;
    src_ptr2 = src_ptr1 + (plane_height * plane_stride) - plane_stride;
    dest_ptr1 = src_ptr1 - (Border * plane_stride);
    dest_ptr2 = src_ptr2 + plane_stride;

    for (i = 0; i < (int)(Border); i++)
    {
        mi_cpu_copy_fast(src_ptr1, dest_ptr1, plane_stride);
        mi_cpu_copy_fast(src_ptr2, dest_ptr2, plane_stride);
        dest_ptr1 += plane_stride;
        dest_ptr2 += plane_stride;
    }

    /***********/
    /* V Plane */
    /***********/

    // copy the left and right most columns out
    src_ptr1 = ybf->v_buffer;
    src_ptr2 = src_ptr1 + plane_width - 1;
    dest_ptr1 = src_ptr1 - Border;
    dest_ptr2 = src_ptr2 + 1;

    for (i = 0; i < plane_height; i++)
    {
        mi_cpu_fill8(dest_ptr1, src_ptr1[0], Border);
        mi_cpu_fill8(dest_ptr2, src_ptr2[0], Border);
        src_ptr1  += plane_stride;
        src_ptr2  += plane_stride;
        dest_ptr1 += plane_stride;
        dest_ptr2 += plane_stride;
    }

    // Now copy the top and bottom source lines into each line of the respective borders
    src_ptr1 = ybf->v_buffer - Border;
    src_ptr2 = src_ptr1 + (plane_height * plane_stride) - plane_stride;
    dest_ptr1 = src_ptr1 - (Border * plane_stride);
    dest_ptr2 = src_ptr2 + plane_stride;

    for (i = 0; i < (int)(Border); i++)
    {
        mi_cpu_copy_fast(src_ptr1, dest_ptr1, plane_stride);
        mi_cpu_copy_fast(src_ptr2, dest_ptr2, plane_stride);
        dest_ptr1 += plane_stride;
        dest_ptr2 += plane_stride;
    }
}



/****************************************************************************
*
*  ROUTINE       : vp8_yv12_copy_frame
*
*  INPUTS        :
*
*  OUTPUTS       : None.
*
*  RETURNS       : void
*
*  FUNCTION      : Copies the source image into the destination image and
*                  updates the destination's UMV borders.
*
*  SPECIAL NOTES : The frames are assumed to be identical in size.
*
****************************************************************************/
void
vp8_yv12_copy_frame(YV12_BUFFER_CONFIG *src_ybc, YV12_BUFFER_CONFIG *dst_ybc)
{
    int yplane_size = (src_ybc->y_height + 2 * src_ybc->border) * (src_ybc->y_stride);
    int mem_size = (yplane_size * 3 / 2) + (src_ybc->y_stride * 2);

    mi_cpu_copy_fast(src_ybc->buffer_alloc, dst_ybc->buffer_alloc, mem_size);

    /*  unsigned char *src_y, *dst_y;
        unsigned char *src_u, *dst_u;
        unsigned char *src_v, *dst_v;

        int yheight, uv_height;
        int ystride, uv_stride;
        int border;
        int yoffset, uvoffset;

        border   = src_ybc->border;
        yheight  = src_ybc->y_height;
        uv_height = src_ybc->uv_height;

        ystride  = src_ybc->y_stride;
        uv_stride = src_ybc->uv_stride;

        yoffset  = border * (ystride + 1);
        uvoffset = border/2 * (uv_stride + 1);

        src_y = src_ybc->y_buffer - yoffset;
        dst_y = dst_ybc->y_buffer - yoffset;
        src_u = src_ybc->u_buffer - uvoffset;
        dst_u = dst_ybc->u_buffer - uvoffset;
        src_v = src_ybc->v_buffer - uvoffset;
        dst_v = dst_ybc->v_buffer - uvoffset;

        mi_cpu_copy_fast (src_y, dst_y, ystride *  (yheight + 2 * border));
        mi_cpu_copy_fast (src_u, dst_u, uv_stride * (uv_height + border));
        mi_cpu_copy_fast (src_v, dst_v, uv_stride * (uv_height + border));
    */
}

#include <nitro/itcm_end.h>
