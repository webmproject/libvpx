/*
 *  Copyright (c) 2010 The VP8 project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license and patent
 *  grant that can be found in the LICENSE file in the root of the source
 *  tree. All contributing project authors may be found in the AUTHORS
 *  file in the root of the source tree.
 */


#include <math.h>
#include "vpx_mem/vpx_mem.h"

#include "quantize.h"
#include "entropy.h"
#include "predictdc.h"

void vp8_fast_quantize_b_c(BLOCK *b, BLOCKD *d)
{
    int i, rc, eob;
    int zbin;
    int x, y, z, sz;
    short *coeff_ptr  = &b->coeff[0];
    short *zbin_ptr   = &b->zbin[0][0];
    short *round_ptr  = &b->round[0][0];
    short *quant_ptr  = &b->quant[0][0];
    short *qcoeff_ptr = d->qcoeff;
    short *dqcoeff_ptr = d->dqcoeff;
    short *dequant_ptr = &d->dequant[0][0];

    vpx_memset(qcoeff_ptr, 0, 32);
    vpx_memset(dqcoeff_ptr, 0, 32);

    eob = -1;

    for (i = 0; i < 16; i++)
    {
        rc   = vp8_default_zig_zag1d[i];
        z    = coeff_ptr[rc];
        zbin = zbin_ptr[rc] ;

        sz = (z >> 31);                                 // sign of z
        x  = (z ^ sz) - sz;                             // x = abs(z)

        if (x >= zbin)
        {
            y  = ((x + round_ptr[rc]) * quant_ptr[rc]) >> 16; // quantize (x)
            x  = (y ^ sz) - sz;                         // get the sign back
            qcoeff_ptr[rc] = x;                          // write to destination
            dqcoeff_ptr[rc] = x * dequant_ptr[rc];        // dequantized value

            if (y)
            {
                eob = i;                                // last nonzero coeffs
            }
        }
    }

    d->eob = eob + 1;

}

void vp8_regular_quantize_b(BLOCK *b, BLOCKD *d)
{
    int i, rc, eob;
    int zbin;
    int x, y, z, sz;
    short *zbin_boost_ptr = &b->zrun_zbin_boost[0];
    short *coeff_ptr  = &b->coeff[0];
    short *zbin_ptr   = &b->zbin[0][0];
    short *round_ptr  = &b->round[0][0];
    short *quant_ptr  = &b->quant[0][0];
    short *qcoeff_ptr = d->qcoeff;
    short *dqcoeff_ptr = d->dqcoeff;
    short *dequant_ptr = &d->dequant[0][0];
    short zbin_oq_value = b->zbin_extra;

    vpx_memset(qcoeff_ptr, 0, 32);
    vpx_memset(dqcoeff_ptr, 0, 32);

    eob = -1;

    for (i = 0; i < 16; i++)
    {
        rc   = vp8_default_zig_zag1d[i];
        z    = coeff_ptr[rc];

        //if ( i == 0 )
        //    zbin = zbin_ptr[rc] + *zbin_boost_ptr + zbin_oq_value/2;
        //else
        zbin = zbin_ptr[rc] + *zbin_boost_ptr + zbin_oq_value;

        zbin_boost_ptr ++;
        sz = (z >> 31);                                 // sign of z
        x  = (z ^ sz) - sz;                             // x = abs(z)

        if (x >= zbin)
        {
            y  = ((x + round_ptr[rc]) * quant_ptr[rc]) >> 16; // quantize (x)
            x  = (y ^ sz) - sz;                         // get the sign back
            qcoeff_ptr[rc]  = x;                         // write to destination
            dqcoeff_ptr[rc] = x * dequant_ptr[rc];        // dequantized value

            if (y)
            {
                eob = i;                                // last nonzero coeffs
                zbin_boost_ptr = &b->zrun_zbin_boost[0];    // reset zero runlength
            }
        }
    }

    d->eob = eob + 1;
}
void vp8_quantize_mby(MACROBLOCK *x)
{
    int i;

    if (x->e_mbd.mbmi.mode != B_PRED && x->e_mbd.mbmi.mode != SPLITMV)
    {
        for (i = 0; i < 16; i++)
        {
            x->quantize_b(&x->block[i], &x->e_mbd.block[i]);
            x->e_mbd.mbmi.mb_skip_coeff &= (x->e_mbd.block[i].eob < 2);
        }

        x->quantize_b(&x->block[24], &x->e_mbd.block[24]);
        x->e_mbd.mbmi.mb_skip_coeff &= (!x->e_mbd.block[24].eob);

    }
    else
    {
        for (i = 0; i < 16; i++)
        {
            x->quantize_b(&x->block[i], &x->e_mbd.block[i]);
            x->e_mbd.mbmi.mb_skip_coeff &= (!x->e_mbd.block[i].eob);
        }
    }
}

void vp8_quantize_mb(MACROBLOCK *x)
{
    int i;

    x->e_mbd.mbmi.mb_skip_coeff = 1;

    if (x->e_mbd.mbmi.mode != B_PRED && x->e_mbd.mbmi.mode != SPLITMV)
    {
        for (i = 0; i < 16; i++)
        {
            x->quantize_b(&x->block[i], &x->e_mbd.block[i]);
            x->e_mbd.mbmi.mb_skip_coeff &= (x->e_mbd.block[i].eob < 2);
        }

        for (i = 16; i < 25; i++)
        {
            x->quantize_b(&x->block[i], &x->e_mbd.block[i]);
            x->e_mbd.mbmi.mb_skip_coeff &= (!x->e_mbd.block[i].eob);
        }
    }
    else
    {
        for (i = 0; i < 24; i++)
        {
            x->quantize_b(&x->block[i], &x->e_mbd.block[i]);
            x->e_mbd.mbmi.mb_skip_coeff &= (!x->e_mbd.block[i].eob);
        }
    }

}


void vp8_quantize_mbuv(MACROBLOCK *x)
{
    int i;

    for (i = 16; i < 24; i++)
    {
        x->quantize_b(&x->block[i], &x->e_mbd.block[i]);
        x->e_mbd.mbmi.mb_skip_coeff &= (!x->e_mbd.block[i].eob);
    }
}

// This function is not currently called
void vp8_quantize_mbrd(MACROBLOCK *x)
{
    int i;

    x->e_mbd.mbmi.mb_skip_coeff = 1;

    if (x->e_mbd.mbmi.mode != B_PRED && x->e_mbd.mbmi.mode != SPLITMV)
    {
        for (i = 0; i < 16; i++)
        {
            x->quantize_brd(&x->block[i], &x->e_mbd.block[i]);
            x->e_mbd.mbmi.mb_skip_coeff &= (x->e_mbd.block[i].eob < 2);
        }

        for (i = 16; i < 25; i++)
        {
            x->quantize_brd(&x->block[i], &x->e_mbd.block[i]);
            x->e_mbd.mbmi.mb_skip_coeff &= (!x->e_mbd.block[i].eob);
        }
    }
    else
    {
        for (i = 0; i < 24; i++)
        {
            x->quantize_brd(&x->block[i], &x->e_mbd.block[i]);
            x->e_mbd.mbmi.mb_skip_coeff &= (!x->e_mbd.block[i].eob);
        }
    }
}

void vp8_quantize_mbuvrd(MACROBLOCK *x)
{
    int i;

    for (i = 16; i < 24; i++)
    {
        x->quantize_brd(&x->block[i], &x->e_mbd.block[i]);
        x->e_mbd.mbmi.mb_skip_coeff &= (!x->e_mbd.block[i].eob);
    }
}

void vp8_quantize_mbyrd(MACROBLOCK *x)
{
    int i;

    if (x->e_mbd.mbmi.mode != B_PRED && x->e_mbd.mbmi.mode != SPLITMV)
    {
        for (i = 0; i < 16; i++)
        {
            x->quantize_brd(&x->block[i], &x->e_mbd.block[i]);
            x->e_mbd.mbmi.mb_skip_coeff &= (x->e_mbd.block[i].eob < 2);
        }

        x->quantize_brd(&x->block[24], &x->e_mbd.block[24]);
        x->e_mbd.mbmi.mb_skip_coeff &= (!x->e_mbd.block[24].eob);

    }
    else
    {
        for (i = 0; i < 16; i++)
        {
            x->quantize_brd(&x->block[i], &x->e_mbd.block[i]);
            x->e_mbd.mbmi.mb_skip_coeff &= (!x->e_mbd.block[i].eob);
        }
    }
}
