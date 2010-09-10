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
 *   Module Title :     yv12config.c
 *
 *   Description  :
 *
 ***************************************************************************/

/****************************************************************************
*  Header Files
****************************************************************************/
#include "vpx_scale/yv12config.h"
#include "vpx_mem/vpx_mem.h"

#include <cdef_bf533.h>

/****************************************************************************
*  Imports
****************************************************************************/
void
extend_memset(void *dst, unsigned char value, unsigned int size);

/****************************************************************************
 *
 ****************************************************************************/
int
vp8_yv12_de_alloc_frame_buffer(YV12_BUFFER_CONFIG *ybf)
{
    if (ybf)
    {
        if (ybf->buffer_alloc)
        {
            duck_free(ybf->buffer_alloc);
        }

        ybf->buffer_alloc = 0;
    }
    else
    {
        return -1;
    }

    return 0;
}

/****************************************************************************
 *
 ****************************************************************************/
int
vp8_yv12_alloc_frame_buffer(YV12_BUFFER_CONFIG *ybf, int width, int height, int border)
{
//NOTE:

    int yplane_size = (height + 2 * border) * (width + 2 * border);
    int uvplane_size = (height / 2 + border) * (width / 2 + border);

    if (ybf)
    {
        vp8_yv12_de_alloc_frame_buffer(ybf);

        ybf->y_width  = width;
        ybf->y_height = height;
        ybf->y_stride = width + 2 * border;

        ybf->uv_width = width / 2;
        ybf->uv_height = height / 2;
        ybf->uv_stride = ybf->uv_width + border;

        ybf->border = border;

        // Added 2 extra lines to framebuffer so that copy12x12 doesn't fail
        // when we have a large motion vector in V on the last v block.
        // Note : We never use these pixels anyway so this doesn't hurt.
        ybf->buffer_alloc = (unsigned char *) duck_memalign(32, (yplane_size * 3 / 2) +  ybf->y_stride , 0);

        if (ybf->buffer_alloc == NULL)
            return -1;

        ybf->y_buffer = ybf->buffer_alloc + border * ybf->y_stride + border;
        ybf->u_buffer = ybf->buffer_alloc + yplane_size + border / 2  * ybf->uv_stride + border / 2;
        ybf->v_buffer = ybf->buffer_alloc + yplane_size + uvplane_size + border / 2  * ybf->uv_stride + border / 2;
    }
    else
    {
        return -2;
    }

    return 0;
}
/****************************************************************************
 *
 ****************************************************************************/
int
vp8_yv12_black_frame_buffer(YV12_BUFFER_CONFIG *ybf)
{
    if (ybf)
    {
        if (ybf->buffer_alloc)
        {
            extend_memset(ybf->y_buffer, 0x0, ybf->y_stride *(ybf->y_height + 2 * ybf->border));
            extend_memset(ybf->u_buffer, 0x80, ybf->uv_stride *(ybf->uv_height + ybf->border));
            extend_memset(ybf->v_buffer, 0x80, ybf->uv_stride *(ybf->uv_height + ybf->border));
        }

        return 0;
    }

    return -1;
}
