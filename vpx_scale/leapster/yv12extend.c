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
//#include <stdlib.h>
#include "vpx_scale/yv12config.h"
#include "vpx_mem/vpx_mem.h"

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
    char *src_ptr1, *src_ptr2;
    char *dest_ptr1, *dest_ptr2;

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
        memset(dest_ptr1, src_ptr1[0], Border);
        memset(dest_ptr2, src_ptr2[0], Border);
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
        memcpy(dest_ptr1, src_ptr1, plane_stride);
        memcpy(dest_ptr2, src_ptr2, plane_stride);
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
        memset(dest_ptr1, src_ptr1[0], Border);
        memset(dest_ptr2, src_ptr2[0], Border);
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
        memcpy(dest_ptr1, src_ptr1, plane_stride);
        memcpy(dest_ptr2, src_ptr2, plane_stride);
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
        memset(dest_ptr1, src_ptr1[0], Border);
        memset(dest_ptr2, src_ptr2[0], Border);
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
        memcpy(dest_ptr1, src_ptr1, plane_stride);
        memcpy(dest_ptr2, src_ptr2, plane_stride);
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
    int row;
    int i;
    unsigned int *source;
    _Uncached unsigned int *dest;
    int height;
    int width;

    height = src_ybc->y_height + (src_ybc->border * 2);
    width =  src_ybc->y_width + (src_ybc->border * 2);
    width /= 4;
    source = (unsigned int *)(src_ybc->y_buffer - (src_ybc->border * src_ybc->y_stride) - src_ybc->border);
    dest = (_Uncached unsigned int *)(dst_ybc->y_buffer - (dst_ybc->border * dst_ybc->y_stride) - dst_ybc->border);

    for (row = 0; row < height; row++)
    {
        for (i = 0; i < width; i++)
        {
            dest[i] = source[i];
        }

        source += width;
        dest   += width;
    }

    height = src_ybc->uv_height + (src_ybc->border);
    width =  src_ybc->uv_width + (src_ybc->border);
    width /= 4;

    source = (unsigned int *)(src_ybc->u_buffer - (src_ybc->border / 2 * src_ybc->uv_stride) - src_ybc->border / 2);
    dest = (_Uncached unsigned int *)(dst_ybc->u_buffer - (dst_ybc->border / 2 * dst_ybc->uv_stride) - dst_ybc->border / 2);

    for (row = 0; row < height; row++)
    {
        for (i = 0; i < width; i++)
        {
            dest[i] = source[i];
        }

        source += width;
        dest   += width;
    }

    source = (unsigned int *)(src_ybc->v_buffer - (src_ybc->border / 2 * src_ybc->uv_stride) - src_ybc->border / 2);
    dest = (_Uncached unsigned int *)(dst_ybc->v_buffer - (dst_ybc->border / 2 * dst_ybc->uv_stride) - dst_ybc->border / 2);

    for (row = 0; row < height; row++)
    {
        for (i = 0; i < width; i++)
        {
            dest[i] = source[i];
        }

        source += width;
        dest   += width;
    }

}
