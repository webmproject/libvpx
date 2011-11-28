/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */


#include "invtrans.h"


void vp8_inverse_transform_b(const vp8_idct_rtcd_vtable_t *rtcd, BLOCKD *b,
                             int pitch)
{
    if (*b->eob > 1)
    {
        IDCT_INVOKE(rtcd, idct16)(b->dqcoeff, b->predictor, pitch,
              *(b->base_dst) + b->dst, b->dst_stride);
    }
    else
    {
        IDCT_INVOKE(rtcd, idct1_scalar_add)(b->dqcoeff[0], b->predictor, pitch,
                         *(b->base_dst) + b->dst, b->dst_stride);
    }

}

void vp8_inverse_transform_mby(const vp8_idct_rtcd_vtable_t *rtcd, MACROBLOCKD *x)
{
    int i;

    if(x->mode_info_context->mbmi.mode != SPLITMV)
    {
        /* do 2nd order transform on the dc block */
        IDCT_INVOKE(rtcd, iwalsh16)(x->block[24].dqcoeff, x->dqcoeff);
    }

    for (i = 0; i < 16; i++)
    {
        vp8_inverse_transform_b(rtcd, &x->block[i], 16);
    }

}
void vp8_inverse_transform_mbuv(const vp8_idct_rtcd_vtable_t *rtcd, MACROBLOCKD *x)
{
    int i;

    for (i = 16; i < 24; i++)
    {
        vp8_inverse_transform_b(rtcd, &x->block[i], 8);
    }

}
