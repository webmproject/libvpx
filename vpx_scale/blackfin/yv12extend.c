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
#include <cdef_bf533.h>

#include "vpx_scale/yv12config.h"
#include "vpx_mem/vpx_mem.h"

/****************************************************************************
*
****************************************************************************/


/****************************************************************************
*
****************************************************************************/
void
extend_memset(void *dst, unsigned char value, unsigned int size)
{
#if 0
    unsigned int quad_value;

    quad_value = (unsigned int) value;
    quad_value |= (unsigned int) value << 8;
    quad_value |= (unsigned int) value << 16;
    quad_value |= (unsigned int) value << 24;
#else
    unsigned short quad_value;

    quad_value = (unsigned int) value;
    quad_value |= (unsigned int) value << 8;
#endif


    if (size / 2 >= 64 * 1024)
        printf("_Extend_memset__________ dma memset is broken\n");

    *p_mdma_s1_start_addr = &quad_value;
    *p_mdma_s1_x_count = size / 2;
    *p_mdma_s1_x_modify = 0x0;
    *p_mdma_d1_start_addr = dst;
    *p_mdma_d1_x_count = size / 2;
    *p_mdma_d1_x_modify = 2;

    *p_mdma_s1_config = DMAEN | WDSIZE_16;
    asm("ssync;");

    *p_mdma_d1_config = DI_EN | DMAEN | WNR | WDSIZE_16;
    asm("ssync;");

    while ((*p_mdma_d1_irq_status & DMA_DONE) == 0);

    *p_mdma_d1_irq_status |= DMA_DONE;
}

/****************************************************************************
*
****************************************************************************/
void
extend_memcpy(void *dst, void *src, unsigned int size)
{
    if (size / 2 >= 64 * 1024)
        printf("_Extend_memcpy__________ dma memcpy is broken\n");


    if ((size & 0x3))
        printf("_)__________ size not a multiple of 4\n");

//32 bit dma here caused some data to be corrupted --- WHY ??????

    *p_mdma_s1_start_addr = src;
    *p_mdma_s1_x_count = size / 2;
    *p_mdma_s1_x_modify = 2;
    *p_mdma_d1_start_addr = dst;
    *p_mdma_d1_x_count = size / 2;
    *p_mdma_d1_x_modify = 2;

    *p_mdma_s1_config = DMAEN | WDSIZE_16;
    asm("ssync;");

    *p_mdma_d1_config = DI_EN | DMAEN | WNR | WDSIZE_16;
    asm("ssync;");

    while ((*p_mdma_d1_irq_status & DMA_DONE) == 0);

    *p_mdma_d1_irq_status |= DMA_DONE;
}

/****************************************************************************
 *
 ****************************************************************************/
void
vp8_yv12_extend_frame_borders(YV12_BUFFER_CONFIG *ybf)
{
#if 1
    int i;
    unsigned char *src_ptr1, *src_ptr2;
    unsigned char *dest_ptr1, *dest_ptr2;

    unsigned int Border;
    int plane_stride;
    int plane_height;
    int plane_width;

    unsigned int quad_sample;
    unsigned int sample;

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
        extend_memset(dest_ptr1, src_ptr1[0], Border);
        extend_memset(dest_ptr2, src_ptr2[0], Border);
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
        extend_memcpy(dest_ptr1, src_ptr1, plane_stride);
        dest_ptr1 += plane_stride;
    }

    for (i = 0; i < (int)Border; i++)
    {
        extend_memcpy(dest_ptr2, src_ptr2, plane_stride);
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
        extend_memset(dest_ptr1, src_ptr1[0], Border);
        extend_memset(dest_ptr2, src_ptr2[0], Border);
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
        extend_memcpy(dest_ptr1, src_ptr1, plane_stride);
        dest_ptr1 += plane_stride;
    }

    for (i = 0; i < (int)(Border); i++)
    {
        extend_memcpy(dest_ptr2, src_ptr2, plane_stride);
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
        extend_memset(dest_ptr1, src_ptr1[0], Border);
        extend_memset(dest_ptr2, src_ptr2[0], Border);
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
        extend_memcpy(dest_ptr1, src_ptr1, plane_stride);
        dest_ptr1 += plane_stride;
    }

    for (i = 0; i < (int)(Border); i++)
    {
        extend_memcpy(dest_ptr2, src_ptr2, plane_stride);
        dest_ptr2 += plane_stride;
    }

#endif
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
#if 1
    int row;
    unsigned char *source, *dest;

    source = src_ybc->y_buffer;
    dest = dst_ybc->y_buffer;

    for (row = 0; row < src_ybc->y_height; row++)
    {
        extend_memcpy(dest, source, src_ybc->y_width);
        source += src_ybc->y_stride;
        dest   += dst_ybc->y_stride;
    }

    source = src_ybc->u_buffer;
    dest = dst_ybc->u_buffer;

    for (row = 0; row < src_ybc->uv_height; row++)
    {
        extend_memcpy(dest, source, src_ybc->uv_width);
        source += src_ybc->uv_stride;
        dest   += dst_ybc->uv_stride;
    }

    source = src_ybc->v_buffer;
    dest = dst_ybc->v_buffer;

    for (row = 0; row < src_ybc->uv_height; row++)
    {
        extend_memcpy(dest, source, src_ybc->uv_width);
        source += src_ybc->uv_stride;
        dest   += dst_ybc->uv_stride;
    }

    vp8_yv12_extend_frame_borders(dst_ybc);

#else
    int row;
    char *source, *dest;
    int height;
    int width;

    height = src_ybc->y_height + (src_ybc->border * 2);
    width =  src_ybc->y_width + (src_ybc->border * 2);
    source = src_ybc->y_buffer;
    dest = dst_ybc->y_buffer;

    for (row = 0; row < height; row++)
    {
        extend_memcpy(dest, source, width);
        source += src_ybc->y_stride;
        dest   += dst_ybc->y_stride;
    }

    height = src_ybc->uv_height + (src_ybc->border);
    width =  src_ybc->uv_width + (src_ybc->border);

    source = src_ybc->u_buffer;
    dest = dst_ybc->u_buffer;

    for (row = 0; row < height; row++)
    {
        extend_memcpy(dest, source, width);
        source += src_ybc->uv_stride;
        dest   += dst_ybc->uv_stride;
    }

    source = src_ybc->v_buffer;
    dest = dst_ybc->v_buffer;

    for (row = 0; row < height; row++)
    {
        extend_memcpy(dest, source, width);
        source += src_ybc->uv_stride;
        dest   += dst_ybc->uv_stride;
    }

#endif

}
