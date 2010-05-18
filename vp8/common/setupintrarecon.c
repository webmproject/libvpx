/*
 *  Copyright (c) 2010 The VP8 project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license and patent
 *  grant that can be found in the LICENSE file in the root of the source
 *  tree. All contributing project authors may be found in the AUTHORS
 *  file in the root of the source tree.
 */


#include "setupintrarecon.h"
#include "vpx_mem/vpx_mem.h"

void vp8_setup_intra_recon(YV12_BUFFER_CONFIG *ybf)
{
    int i;

    // set up frame new frame for intra coded blocks
    vpx_memset(ybf->y_buffer - 1 - 2 * ybf->y_stride, 127, ybf->y_width + 5);
    vpx_memset(ybf->y_buffer - 1 - ybf->y_stride, 127, ybf->y_width + 5);

    for (i = 0; i < ybf->y_height; i++)
        ybf->y_buffer[ybf->y_stride *i - 1] = (unsigned char) 129;

    vpx_memset(ybf->u_buffer - 1 - 2 * ybf->uv_stride, 127, ybf->uv_width + 5);
    vpx_memset(ybf->u_buffer - 1 - ybf->uv_stride, 127, ybf->uv_width + 5);

    for (i = 0; i < ybf->uv_height; i++)
        ybf->u_buffer[ybf->uv_stride *i - 1] = (unsigned char) 129;

    vpx_memset(ybf->v_buffer - 1 - 2 * ybf->uv_stride, 127, ybf->uv_width + 5);
    vpx_memset(ybf->v_buffer - 1 - ybf->uv_stride, 127, ybf->uv_width + 5);

    for (i = 0; i < ybf->uv_height; i++)
        ybf->v_buffer[ybf->uv_stride *i - 1] = (unsigned char) 129;

}
