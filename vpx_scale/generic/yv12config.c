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
#include "vpx_mem/vpx_mem.h"

/****************************************************************************
*  Exports
****************************************************************************/

/****************************************************************************
 *
 ****************************************************************************/
int
vp8_yv12_de_alloc_frame_buffer(YV12_BUFFER_CONFIG *ybf)
{
    if (ybf)
    {
        vpx_free(ybf->buffer_alloc);

        /* buffer_alloc isn't accessed by most functions.  Rather y_buffer,
          u_buffer and v_buffer point to buffer_alloc and are used.  Clear out
          all of this so that a freed pointer isn't inadvertently used */
        vpx_memset (ybf, 0, sizeof (YV12_BUFFER_CONFIG));
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
/*NOTE:*/

    if (ybf)
    {
        int uv_width = width >> 1;
        int uv_height = height >> 1;
        int yplane_size = (height + 2 * border) * (width + 2 * border);
        int uvplane_size = (uv_height + border) * (uv_width + border);

        vp8_yv12_de_alloc_frame_buffer(ybf);

        /* only support allocating buffers that have
          a height and width that are multiples of 16 */
        if ((width & 0xf) | (height & 0xf))
            return -3;

        ybf->y_width  = width;
        ybf->y_height = height;
        ybf->y_stride = width + 2 * border;

        ybf->uv_width = uv_width;
        ybf->uv_height = uv_height;
        ybf->uv_stride = uv_width + border;

        ybf->border = border;
        ybf->frame_size = yplane_size + 2 * uvplane_size;

        ybf->buffer_alloc = (unsigned char *) vpx_memalign(32, ybf->frame_size);

        if (ybf->buffer_alloc == NULL)
            return -1;

        ybf->y_buffer = ybf->buffer_alloc + (border * ybf->y_stride) + border;
        ybf->u_buffer = ybf->buffer_alloc + yplane_size + (border / 2  * ybf->uv_stride) + border / 2;
        ybf->v_buffer = ybf->buffer_alloc + yplane_size + uvplane_size + (border / 2  * ybf->uv_stride) + border / 2;

        ybf->corrupted = 0; /* assume not currupted by errors */
    }
    else
    {
        return -2;
    }

    return 0;
}
