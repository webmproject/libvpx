/*
 *  Copyright (c) 2010 The VP8 project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license and patent
 *  grant that can be found in the LICENSE file in the root of the source
 *  tree. All contributing project authors may be found in the AUTHORS
 *  file in the root of the source tree.
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
#include "csl_dat.h"
#include "vpx_scale/yv12config.h"
#include "vpx_mem/vpx_mem.h"

/****************************************************************************
*  Exports
****************************************************************************/
#define UINT8 unsigned char
#define UINT32 unsigned int


static inline
void copy_yleft_right_border(
    UINT8 *restrict src_ptr1,
    UINT8 *restrict src_ptr2,
    UINT8 *restrict dest_ptr1,
    UINT8 *restrict dest_ptr2,
    UINT32  plane_height,
    UINT32  plane_stride
)
{
    UINT32 left, right, left2, left4, right2, right4;
    double dl, dr;
    int i;

#pragma MUST_ITERATE(16,16,16)

    for (i = 0; i < plane_height; i++)
    {
        left  = src_ptr1[0];
        right = src_ptr2[0];

        left2 = _pack2(left, left);
        left4 = _packl4(left2, left2);

        right2 = _pack2(right, right);
        right4 = _packl4(right2, right2);

        dl = _itod(left4, left4);
        dr = _itod(right4, right4);

        _amemd8(&dest_ptr1[ 0]) = dl;
        _amemd8(&dest_ptr2[ 0]) = dr;

        _amemd8(&dest_ptr1[ 8]) = dl;
        _amemd8(&dest_ptr2[ 8]) = dr;

        _amemd8(&dest_ptr1[16]) = dl;
        _amemd8(&dest_ptr2[16]) = dr;

        _amemd8(&dest_ptr1[24]) = dl;
        _amemd8(&dest_ptr2[24]) = dr;

        _amemd8(&dest_ptr1[32]) = dl;
        _amemd8(&dest_ptr2[32]) = dr;

        _amemd8(&dest_ptr1[40]) = dl;
        _amemd8(&dest_ptr2[40]) = dr;


        src_ptr1 += plane_stride;
        src_ptr2 += plane_stride;
        dest_ptr1 += plane_stride;
        dest_ptr2 += plane_stride;
    }
}
/****************************************************************************
 *
 *
 ****************************************************************************/
static
void copy_uvleft_right_border(
    UINT8 *restrict src_ptr1,
    UINT8 *restrict src_ptr2,
    UINT8 *restrict dest_ptr1,
    UINT8 *restrict dest_ptr2,
    UINT32  plane_height,
    UINT32  plane_stride
)
{
    UINT32 left, right, left2, left4, right2, right4;
    double dl, dr;
    int i;

#pragma MUST_ITERATE(8,8 ,8)

    for (i = 0; i < plane_height; i++)
    {
        left  = src_ptr1[0];
        right = src_ptr2[0];

        left2 = _pack2(left, left);
        left4 = _packl4(left2, left2);

        right2 = _pack2(right, right);
        right4 = _packl4(right2, right2);

        dl = _itod(left4, left4);
        dr = _itod(right4, right4);

        _amemd8(&dest_ptr1[ 0]) = dl;
        _amemd8(&dest_ptr2[ 0]) = dr;

        _amemd8(&dest_ptr1[ 8]) = dl;
        _amemd8(&dest_ptr2[ 8]) = dr;

        _amemd8(&dest_ptr1[16]) = dl;
        _amemd8(&dest_ptr2[16]) = dr;


        src_ptr1 += plane_stride;
        src_ptr2 += plane_stride;
        dest_ptr1 += plane_stride;
        dest_ptr2 += plane_stride;
    }
}
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

#if 1
    // copy the left and right most columns out
    src_ptr1 = ybf->y_buffer;
    src_ptr2 = src_ptr1 + plane_width - 1;
    dest_ptr1 = src_ptr1 - Border;
    dest_ptr2 = src_ptr2 + 1;
    copy_yleft_right_border(src_ptr1, src_ptr2, dest_ptr1, dest_ptr2, plane_height, plane_stride);
#endif

    // Now copy the top and bottom source lines into each line of the respective borders
    src_ptr1 = ybf->y_buffer - Border;
    src_ptr2 = src_ptr1 + (plane_height * plane_stride) - plane_stride;
    dest_ptr1 = src_ptr1 - (Border * plane_stride);
    dest_ptr2 = src_ptr2 + plane_stride;

    for (i = 0; i < (int)Border; i++)
    {
        vpx_memcpy(dest_ptr1, src_ptr1, plane_stride);
        vpx_memcpy(dest_ptr2, src_ptr2, plane_stride);
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
#if 1
    // copy the left and right most columns out
    src_ptr1 = ybf->u_buffer;
    src_ptr2 = src_ptr1 + plane_width - 1;
    dest_ptr1 = src_ptr1 - Border;
    dest_ptr2 = src_ptr2 + 1;

    copy_uvleft_right_border(src_ptr1, src_ptr2, dest_ptr1, dest_ptr2, plane_height, plane_stride);


#endif

    // Now copy the top and bottom source lines into each line of the respective borders
    src_ptr1 = ybf->u_buffer - Border;
    src_ptr2 = src_ptr1 + (plane_height * plane_stride) - plane_stride;
    dest_ptr1 = src_ptr1 - (Border * plane_stride);
    dest_ptr2 = src_ptr2 + plane_stride;

    for (i = 0; i < (int)(Border); i++)
    {
        vpx_memcpy(dest_ptr1, src_ptr1, plane_stride);
        vpx_memcpy(dest_ptr2, src_ptr2, plane_stride);
        dest_ptr1 += plane_stride;
        dest_ptr2 += plane_stride;
    }

    /***********/
    /* V Plane */
    /***********/
#if 1
    // copy the left and right most columns out
    src_ptr1 = ybf->v_buffer;
    src_ptr2 = src_ptr1 + plane_width - 1;
    dest_ptr1 = src_ptr1 - Border;
    dest_ptr2 = src_ptr2 + 1;

    copy_uvleft_right_border(src_ptr1, src_ptr2, dest_ptr1, dest_ptr2, plane_height, plane_stride);

#endif

    // Now copy the top and bottom source lines into each line of the respective borders
    src_ptr1 = ybf->v_buffer - Border;
    src_ptr2 = src_ptr1 + (plane_height * plane_stride) - plane_stride;
    dest_ptr1 = src_ptr1 - (Border * plane_stride);
    dest_ptr2 = src_ptr2 + plane_stride;

    for (i = 0; i < (int)(Border); i++)
    {
        vpx_memcpy(dest_ptr1, src_ptr1, plane_stride);
        vpx_memcpy(dest_ptr2, src_ptr2, plane_stride);
        dest_ptr1 += plane_stride;
        dest_ptr2 += plane_stride;
    }
}
/****************************************************************************
 *
 ****************************************************************************/
void
vpxyv12_extend_frame_tbborders(YV12_BUFFER_CONFIG *ybf)
{
    int i;
    unsigned char *src_ptr1, *src_ptr2;
    unsigned char *dest_ptr1, *dest_ptr2;
    int tid1, tid2;

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


    // Now copy the top and bottom source lines into each line of the respective borders
    src_ptr1 = ybf->y_buffer - Border;
    src_ptr2 = src_ptr1 + (plane_height * plane_stride) - plane_stride;
    dest_ptr1 = src_ptr1 - (Border * plane_stride);
    dest_ptr2 = src_ptr2 + plane_stride;


    for (i = 0; i < (int)Border; i++)
    {
        dat_copy(src_ptr1, dest_ptr1, plane_stride);
        dat_copy(src_ptr2, dest_ptr2, plane_stride);
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
    // Now copy the top and bottom source lines into each line of the respective borders
    src_ptr1 = ybf->u_buffer - Border;
    src_ptr2 = src_ptr1 + (plane_height * plane_stride) - plane_stride;
    dest_ptr1 = src_ptr1 - (Border * plane_stride);
    dest_ptr2 = src_ptr2 + plane_stride;

    for (i = 0; i < (int)(Border); i++)
    {
        dat_copy(src_ptr1, dest_ptr1, plane_stride);
        dat_copy(src_ptr2, dest_ptr2, plane_stride);
        dest_ptr1 += plane_stride;
        dest_ptr2 += plane_stride;
    }

    /***********/
    /* V Plane */
    /***********/
    // Now copy the top and bottom source lines into each line of the respective borders
    src_ptr1 = ybf->v_buffer - Border;
    src_ptr2 = src_ptr1 + (plane_height * plane_stride) - plane_stride;
    dest_ptr1 = src_ptr1 - (Border * plane_stride);
    dest_ptr2 = src_ptr2 + plane_stride;

    for (i = 0; i < (int)(Border); i++)
    {
        tid1 = dat_copy(src_ptr1, dest_ptr1, plane_stride);
        tid2 = dat_copy(src_ptr2, dest_ptr2, plane_stride);
        dest_ptr1 += plane_stride;
        dest_ptr2 += plane_stride;
    }

    dat_wait(tid1);
    dat_wait(tid2);
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
 *                  updates the destination's UMV borders.  Because the
 *                  borders have been update prior to this so the whole frame
 *                  is copied, borders and all.  This is also to circumvent
 *                  using copy_left_right Border functions when copying data
 *                  between L2 and main memory.  When that occurs a cache
 *                  clean needs to be done, which would require invalidating
 *                  an entire frame.
 *
 *  SPECIAL NOTES : The frames are assumed to be identical in size.
 *
 ****************************************************************************/
void
vpxyv12_copy_frame_dma(YV12_BUFFER_CONFIG *src_ybc, YV12_BUFFER_CONFIG *dst_ybc)
{
    int yheight, uv_height;
    int ystride, uv_stride;
    int border;
    int yoffset, uvoffset;

    border = src_ybc->border;
    yheight = src_ybc->y_height;
    uv_height = src_ybc->uv_height;

    ystride = src_ybc->y_stride;
    uv_stride = src_ybc->uv_stride;

    yoffset = border * (ystride + 1);
    uvoffset = border / 2 * (uv_stride + 1);

    dat_copy2d(DAT_2D2D,
               src_ybc->y_buffer - yoffset,
               dst_ybc->y_buffer - yoffset,
               ystride,
               yheight + 2 * border,
               ystride);
    dat_copy2d(DAT_2D2D,
               src_ybc->u_buffer - uvoffset,
               dst_ybc->u_buffer - uvoffset,
               uv_stride,
               uv_height + border,
               uv_stride);
    dat_copy2d(DAT_2D2D,
               src_ybc->v_buffer - uvoffset,
               dst_ybc->v_buffer - uvoffset,
               uv_stride,
               uv_height + border,
               uv_stride);

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
    unsigned char *source, *dest;

    source = src_ybc->y_buffer;
    dest = dst_ybc->y_buffer;

    for (row = 0; row < src_ybc->y_height; row++)
    {
        vpx_memcpy(dest, source, src_ybc->y_width);
        source += src_ybc->y_stride;
        dest   += dst_ybc->y_stride;
    }

    source = src_ybc->u_buffer;
    dest = dst_ybc->u_buffer;

    for (row = 0; row < src_ybc->uv_height; row++)
    {
        vpx_memcpy(dest, source, src_ybc->uv_width);
        source += src_ybc->uv_stride;
        dest   += dst_ybc->uv_stride;
    }

    source = src_ybc->v_buffer;
    dest = dst_ybc->v_buffer;

    for (row = 0; row < src_ybc->uv_height; row++)
    {
        vpx_memcpy(dest, source, src_ybc->uv_width);
        source += src_ybc->uv_stride;
        dest   += dst_ybc->uv_stride;
    }

    vp8_yv12_extend_frame_borders(dst_ybc);
}
