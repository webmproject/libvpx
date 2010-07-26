/*
 *  Copyright (c) 2010 The VP8 project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */


#include "vpx_ports/config.h"
#include "dequantize.h"
#include "predictdc.h"
#include "idct.h"
#include "vpx_mem/vpx_mem.h"

extern void vp8_short_idct4x4llm_c(short *input, short *output, int pitch) ;
extern void vp8_short_idct4x4llm_1_c(short *input, short *output, int pitch);


void vp8_dequantize_b_c(BLOCKD *d)
{
    int i;
    short *DQ  = d->dqcoeff;
    short *Q   = d->qcoeff;
    short *DQC = &d->dequant[0][0];

    for (i = 0; i < 16; i++)
    {
        DQ[i] = Q[i] * DQC[i];
    }
}

void vp8_dequant_idct_add_c(short *input, short *dq, unsigned char *pred, unsigned char *dest, int pitch, int stride)
{
    // output needs to be at least pitch * 4 for vp8_short_idct4x4llm_c to work properly
    short output[16*4];
    short *diff_ptr = output;
    int r, c;
    int i;

    for (i = 0; i < 16; i++)
    {
        input[i] = dq[i] * input[i];
    }

    vp8_short_idct4x4llm_c(input, output, pitch*2);

    vpx_memset(input, 0, 32);

    for (r = 0; r < 4; r++)
    {
        for (c = 0; c < 4; c++)
        {
            int a = diff_ptr[c] + pred[c];

            if (a < 0)
                a = 0;

            if (a > 255)
                a = 255;

            dest[c] = (unsigned char) a;
        }

        dest += stride;
        diff_ptr += pitch;
        pred += pitch;
    }
}

void vp8_dequant_dc_idct_add_c(short *input, short *dq, unsigned char *pred, unsigned char *dest, int pitch, int stride, int Dc)
{
    int i;
    // output needs to be at least pitch * 4 for vp8_short_idct4x4llm_c to work properly
    short output[16*4];
    short *diff_ptr = output;
    int r, c;

    input[0] = (short)Dc;

    for (i = 1; i < 16; i++)
    {
        input[i] = dq[i] * input[i];
    }

    vp8_short_idct4x4llm_c(input, output, pitch*2);

    vpx_memset(input, 0, 32);

    for (r = 0; r < 4; r++)
    {
        for (c = 0; c < 4; c++)
        {
            int a = diff_ptr[c] + pred[c];

            if (a < 0)
                a = 0;

            if (a > 255)
                a = 255;

            dest[c] = (unsigned char) a;
        }

        dest += stride;
        diff_ptr += pitch;
        pred += pitch;
    }
}
