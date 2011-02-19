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

#include "quantize.h"
#include "vp8/common/entropy.h"

#define EXACT_QUANT

#ifdef EXACT_FASTQUANT
void vp8_fast_quantize_b_c(BLOCK *b, BLOCKD *d)
{
    int i, rc, eob;
    int zbin;
    int x, y, z, sz;
    short *coeff_ptr       = b->coeff;
    short *zbin_ptr        = b->zbin;
    short *round_ptr       = b->round;
    short *quant_ptr       = b->quant_fast;
    short *quant_shift_ptr = b->quant_shift;
    short *qcoeff_ptr      = d->qcoeff;
    short *dqcoeff_ptr     = d->dqcoeff;
    short *dequant_ptr     = d->dequant;

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
            x += round_ptr[rc];
            y  = (((x * quant_ptr[rc]) >> 16) + x)
                 >> quant_shift_ptr[rc];                // quantize (x)
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

#else

void vp8_fast_quantize_b_c(BLOCK *b, BLOCKD *d)
{
    int i, rc, eob;
    int x, y, z, sz;
    short *coeff_ptr   = b->coeff;
    short *round_ptr   = b->round;
    short *quant_ptr   = b->quant_fast;
    short *qcoeff_ptr  = d->qcoeff;
    short *dqcoeff_ptr = d->dqcoeff;
    short *dequant_ptr = d->dequant;

    eob = -1;
    for (i = 0; i < 16; i++)
    {
        rc   = vp8_default_zig_zag1d[i];
        z    = coeff_ptr[rc];

        sz = (z >> 31);                                 // sign of z
        x  = (z ^ sz) - sz;                             // x = abs(z)

        y  = ((x + round_ptr[rc]) * quant_ptr[rc]) >> 16; // quantize (x)
        x  = (y ^ sz) - sz;                         // get the sign back
        qcoeff_ptr[rc] = x;                          // write to destination
        dqcoeff_ptr[rc] = x * dequant_ptr[rc];        // dequantized value

        if (y)
        {
            eob = i;                                // last nonzero coeffs
        }
    }
    d->eob = eob + 1;
}

#endif

#ifdef EXACT_QUANT
void vp8_regular_quantize_b(BLOCK *b, BLOCKD *d)
{
    int i, rc, eob;
    int zbin;
    int x, y, z, sz;
    short *zbin_boost_ptr  = b->zrun_zbin_boost;
    short *coeff_ptr       = b->coeff;
    short *zbin_ptr        = b->zbin;
    short *round_ptr       = b->round;
    short *quant_ptr       = b->quant;
    short *quant_shift_ptr = b->quant_shift;
    short *qcoeff_ptr      = d->qcoeff;
    short *dqcoeff_ptr     = d->dqcoeff;
    short *dequant_ptr     = d->dequant;
    short zbin_oq_value    = b->zbin_extra;

    vpx_memset(qcoeff_ptr, 0, 32);
    vpx_memset(dqcoeff_ptr, 0, 32);

    eob = -1;

    for (i = 0; i < 16; i++)
    {
        rc   = vp8_default_zig_zag1d[i];
        z    = coeff_ptr[rc];

        zbin = zbin_ptr[rc] + *zbin_boost_ptr + zbin_oq_value;

        zbin_boost_ptr ++;
        sz = (z >> 31);                                 // sign of z
        x  = (z ^ sz) - sz;                             // x = abs(z)

        if (x >= zbin)
        {
            x += round_ptr[rc];
            y  = (((x * quant_ptr[rc]) >> 16) + x)
                 >> quant_shift_ptr[rc];                // quantize (x)
            x  = (y ^ sz) - sz;                         // get the sign back
            qcoeff_ptr[rc]  = x;                        // write to destination
            dqcoeff_ptr[rc] = x * dequant_ptr[rc];      // dequantized value

            if (y)
            {
                eob = i;                                // last nonzero coeffs
                zbin_boost_ptr = b->zrun_zbin_boost;    // reset zero runlength
            }
        }
    }

    d->eob = eob + 1;
}

/* Perform regular quantization, with unbiased rounding and no zero bin. */
void vp8_strict_quantize_b(BLOCK *b, BLOCKD *d)
{
    int i;
    int rc;
    int eob;
    int x;
    int y;
    int z;
    int sz;
    short *coeff_ptr;
    short *quant_ptr;
    short *quant_shift_ptr;
    short *qcoeff_ptr;
    short *dqcoeff_ptr;
    short *dequant_ptr;

    coeff_ptr       = b->coeff;
    quant_ptr       = b->quant;
    quant_shift_ptr = b->quant_shift;
    qcoeff_ptr      = d->qcoeff;
    dqcoeff_ptr     = d->dqcoeff;
    dequant_ptr     = d->dequant;
    eob = - 1;
    vpx_memset(qcoeff_ptr, 0, 32);
    vpx_memset(dqcoeff_ptr, 0, 32);
    for (i = 0; i < 16; i++)
    {
        int dq;
        int round;

        /*TODO: These arrays should be stored in zig-zag order.*/
        rc = vp8_default_zig_zag1d[i];
        z = coeff_ptr[rc];
        dq = dequant_ptr[rc];
        round = dq >> 1;
        /* Sign of z. */
        sz = -(z < 0);
        x = (z + sz) ^ sz;
        x += round;
        if (x >= dq)
        {
            /* Quantize x. */
            y  = (((x * quant_ptr[rc]) >> 16) + x) >> quant_shift_ptr[rc];
            /* Put the sign back. */
            x = (y + sz) ^ sz;
            /* Save the coefficient and its dequantized value. */
            qcoeff_ptr[rc] = x;
            dqcoeff_ptr[rc] = x * dq;
            /* Remember the last non-zero coefficient. */
            if (y)
                eob = i;
        }
    }

    d->eob = eob + 1;
}

#else

void vp8_regular_quantize_b(BLOCK *b, BLOCKD *d)
{
    int i, rc, eob;
    int zbin;
    int x, y, z, sz;
    short *zbin_boost_ptr = b->zrun_zbin_boost;
    short *coeff_ptr      = b->coeff;
    short *zbin_ptr       = b->zbin;
    short *round_ptr      = b->round;
    short *quant_ptr      = b->quant;
    short *qcoeff_ptr     = d->qcoeff;
    short *dqcoeff_ptr    = d->dqcoeff;
    short *dequant_ptr    = d->dequant;
    short zbin_oq_value   = b->zbin_extra;

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

#endif

void vp8_quantize_mby(MACROBLOCK *x)
{
    int i;
    int has_2nd_order = (x->e_mbd.mode_info_context->mbmi.mode != B_PRED
        && x->e_mbd.mode_info_context->mbmi.mode != SPLITMV);

    for (i = 0; i < 16; i++)
        x->quantize_b(&x->block[i], &x->e_mbd.block[i]);

    if(has_2nd_order)
        x->quantize_b(&x->block[24], &x->e_mbd.block[24]);
}

void vp8_quantize_mb(MACROBLOCK *x)
{
    int i;
    int has_2nd_order=(x->e_mbd.mode_info_context->mbmi.mode != B_PRED
        && x->e_mbd.mode_info_context->mbmi.mode != SPLITMV);

    for (i = 0; i < 24+has_2nd_order; i++)
        x->quantize_b(&x->block[i], &x->e_mbd.block[i]);
}


void vp8_quantize_mbuv(MACROBLOCK *x)
{
    int i;

    for (i = 16; i < 24; i++)
        x->quantize_b(&x->block[i], &x->e_mbd.block[i]);
}
