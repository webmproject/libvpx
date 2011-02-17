/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */


#include <math.h>
#include "vpx_mem/vpx_mem.h"

#include "vp8/encoder/quantize.h"
#include "vp8/common/entropy.h"

DECLARE_ALIGNED(16, const short, vp8_rvsplus1_default_zig_zag1d[16]) =
{
    1,  2,  6,  7,
    3,  5,  8,  13,
    4,  9,  12, 14,
    10, 11, 15, 16,
};


extern int vp8_fast_quantize_b_neon_func(short *coeff_ptr, short *zbin_ptr, short *qcoeff_ptr, short *dqcoeff_ptr, short *dequant_ptr, const short *scan_mask, short *round_ptr, short *quant_ptr);

void vp8_fast_quantize_b_neon(BLOCK *b, BLOCKD *d)
{
    d->eob = vp8_fast_quantize_b_neon_func(b->coeff, b->zbin, d->qcoeff, d->dqcoeff, d->dequant, vp8_rvsplus1_default_zig_zag1d, b->round, b->quant_fast);
}

/*
//neon code is written according to the following rewritten c code
void vp8_fast_quantize_b_neon(BLOCK *b,BLOCKD *d)
{
    int i, rc, eob;
    int zbin;
    int x, x1, y, z, sz;
    short *coeff_ptr  = &b->Coeff[0];
    short *zbin_ptr   = &b->Zbin[0][0];
    short *round_ptr  = &b->Round[0][0];
    short *quant_ptr  = &b->Quant[0][0];
    short *qcoeff_ptr = d->qcoeff;
    short *dqcoeff_ptr= d->dqcoeff;
    short *dequant_ptr= &d->Dequant[0][0];

    eob = 0;

    for(i=0;i<16;i++)
    {
        z    = coeff_ptr[i];
        zbin = zbin_ptr[i] ;
        x  = abs(z);                                    // x = abs(z)

        if(x>=zbin)
        {
            sz = (z>>31);                               // sign of z
            y  = ((x+round_ptr[i])*quant_ptr[i])>>16;     // quantize (x)
            x1  = (y^sz) - sz;                          // get the sign back

            qcoeff_ptr[i] = x1;                          // write to destination
            dqcoeff_ptr[i] = x1 * dequant_ptr[i];         // dequantized value

            if(y)
            {
                if(eob<vp8_rvsplus1_default_zig_zag1d[i])
                    eob=(int)vp8_rvsplus1_default_zig_zag1d[i];         // last nonzero coeffs
            }
        }else
        {
            qcoeff_ptr[i] = 0;                          // write to destination
            dqcoeff_ptr[i] = 0;         // dequantized value
        }
    }
        d->eob = eob;
}
*/
