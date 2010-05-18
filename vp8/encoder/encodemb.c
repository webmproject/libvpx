/*
 *  Copyright (c) 2010 The VP8 project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license and patent
 *  grant that can be found in the LICENSE file in the root of the source
 *  tree. All contributing project authors may be found in the AUTHORS
 *  file in the root of the source tree.
 */


#include "vpx_ports/config.h"
#include "encodemb.h"
#include "reconinter.h"
#include "quantize.h"
#include "invtrans.h"
#include "recon.h"
#include "reconintra.h"
#include "dct.h"
#include "vpx_mem/vpx_mem.h"

#if CONFIG_RUNTIME_CPU_DETECT
#define IF_RTCD(x) (x)
#else
#define IF_RTCD(x) NULL
#endif
void vp8_subtract_b_c(BLOCK *be, BLOCKD *bd, int pitch)
{
    unsigned char *src_ptr = (*(be->base_src) + be->src);
    short *diff_ptr = be->src_diff;
    unsigned char *pred_ptr = bd->predictor;
    int src_stride = be->src_stride;

    int r, c;

    for (r = 0; r < 4; r++)
    {
        for (c = 0; c < 4; c++)
        {
            diff_ptr[c] = src_ptr[c] - pred_ptr[c];
        }

        diff_ptr += pitch;
        pred_ptr += pitch;
        src_ptr  += src_stride;
    }
}

void vp8_subtract_mbuv_c(short *diff, unsigned char *usrc, unsigned char *vsrc, unsigned char *pred, int stride)
{
    short *udiff = diff + 256;
    short *vdiff = diff + 320;
    unsigned char *upred = pred + 256;
    unsigned char *vpred = pred + 320;

    int r, c;

    for (r = 0; r < 8; r++)
    {
        for (c = 0; c < 8; c++)
        {
            udiff[c] = usrc[c] - upred[c];
        }

        udiff += 8;
        upred += 8;
        usrc  += stride;
    }

    for (r = 0; r < 8; r++)
    {
        for (c = 0; c < 8; c++)
        {
            vdiff[c] = vsrc[c] - vpred[c];
        }

        vdiff += 8;
        vpred += 8;
        vsrc  += stride;
    }
}

void vp8_subtract_mby_c(short *diff, unsigned char *src, unsigned char *pred, int stride)
{
    int r, c;

    for (r = 0; r < 16; r++)
    {
        for (c = 0; c < 16; c++)
        {
            diff[c] = src[c] - pred[c];
        }

        diff += 16;
        pred += 16;
        src  += stride;
    }
}

static void vp8_subtract_mb(const VP8_ENCODER_RTCD *rtcd, MACROBLOCK *x)
{
    ENCODEMB_INVOKE(&rtcd->encodemb, submby)(x->src_diff, x->src.y_buffer, x->e_mbd.predictor, x->src.y_stride);
    ENCODEMB_INVOKE(&rtcd->encodemb, submbuv)(x->src_diff, x->src.u_buffer, x->src.v_buffer, x->e_mbd.predictor, x->src.uv_stride);
}

void vp8_build_dcblock(MACROBLOCK *x)
{
    short *src_diff_ptr = &x->src_diff[384];
    int i;

    for (i = 0; i < 16; i++)
    {
        src_diff_ptr[i] = x->coeff[i * 16];
    }
}

void vp8_transform_mbuv(MACROBLOCK *x)
{
    int i;

    for (i = 16; i < 24; i += 2)
    {
        x->vp8_short_fdct8x4(&x->block[i].src_diff[0], &x->block[i].coeff[0], 16);
    }
}

void vp8_transform_mbuvrd(MACROBLOCK *x)
{
    int i;

    for (i = 16; i < 24; i += 2)
    {
        x->short_fdct8x4rd(&x->block[i].src_diff[0], &x->block[i].coeff[0], 16);
    }
}

void vp8_transform_intra_mby(MACROBLOCK *x)
{
    int i;

    for (i = 0; i < 16; i += 2)
    {
        x->vp8_short_fdct8x4(&x->block[i].src_diff[0], &x->block[i].coeff[0], 32);
    }

    // build dc block from 16 y dc values
    vp8_build_dcblock(x);

    // do 2nd order transform on the dc block
    x->short_walsh4x4(&x->block[24].src_diff[0], &x->block[24].coeff[0], 8);

}

void vp8_transform_intra_mbyrd(MACROBLOCK *x)
{
    int i;

    for (i = 0; i < 16; i += 2)
    {
        x->short_fdct8x4rd(&x->block[i].src_diff[0], &x->block[i].coeff[0], 32);
    }

    // build dc block from 16 y dc values
    vp8_build_dcblock(x);

    // do 2nd order transform on the dc block
    x->short_walsh4x4(&x->block[24].src_diff[0], &x->block[24].coeff[0], 8);
}

void vp8_transform_mb(MACROBLOCK *x)
{
    int i;

    for (i = 0; i < 16; i += 2)
    {
        x->vp8_short_fdct8x4(&x->block[i].src_diff[0], &x->block[i].coeff[0], 32);
    }

    // build dc block from 16 y dc values
    if (x->e_mbd.mbmi.mode != SPLITMV)
        vp8_build_dcblock(x);

    for (i = 16; i < 24; i += 2)
    {
        x->vp8_short_fdct8x4(&x->block[i].src_diff[0], &x->block[i].coeff[0], 16);
    }

    // do 2nd order transform on the dc block
    if (x->e_mbd.mbmi.mode != SPLITMV)
        x->short_walsh4x4(&x->block[24].src_diff[0], &x->block[24].coeff[0], 8);

}

void vp8_transform_mby(MACROBLOCK *x)
{
    int i;

    for (i = 0; i < 16; i += 2)
    {
        x->vp8_short_fdct8x4(&x->block[i].src_diff[0], &x->block[i].coeff[0], 32);
    }

    // build dc block from 16 y dc values
    if (x->e_mbd.mbmi.mode != SPLITMV)
    {
        vp8_build_dcblock(x);
        x->short_walsh4x4(&x->block[24].src_diff[0], &x->block[24].coeff[0], 8);
    }
}

void vp8_transform_mbrd(MACROBLOCK *x)
{
    int i;

    for (i = 0; i < 16; i += 2)
    {
        x->short_fdct8x4rd(&x->block[i].src_diff[0], &x->block[i].coeff[0], 32);
    }

    // build dc block from 16 y dc values
    if (x->e_mbd.mbmi.mode != SPLITMV)
        vp8_build_dcblock(x);

    for (i = 16; i < 24; i += 2)
    {
        x->short_fdct8x4rd(&x->block[i].src_diff[0], &x->block[i].coeff[0], 16);
    }

    // do 2nd order transform on the dc block
    if (x->e_mbd.mbmi.mode != SPLITMV)
        x->short_walsh4x4(&x->block[24].src_diff[0], &x->block[24].coeff[0], 8);
}

void vp8_stuff_inter16x16(MACROBLOCK *x)
{
    vp8_build_inter_predictors_mb_s(&x->e_mbd);
    /*
        // recon = copy from predictors to destination
        {
            BLOCKD *b = &x->e_mbd.block[0];
            unsigned char *pred_ptr = b->predictor;
            unsigned char *dst_ptr = *(b->base_dst) + b->dst;
            int stride = b->dst_stride;

            int i;
            for(i=0;i<16;i++)
                vpx_memcpy(dst_ptr+i*stride,pred_ptr+16*i,16);

            b = &x->e_mbd.block[16];
            pred_ptr = b->predictor;
            dst_ptr = *(b->base_dst) + b->dst;
            stride = b->dst_stride;

            for(i=0;i<8;i++)
                vpx_memcpy(dst_ptr+i*stride,pred_ptr+8*i,8);

            b = &x->e_mbd.block[20];
            pred_ptr = b->predictor;
            dst_ptr = *(b->base_dst) + b->dst;
            stride = b->dst_stride;

            for(i=0;i<8;i++)
                vpx_memcpy(dst_ptr+i*stride,pred_ptr+8*i,8);
        }
    */
}

#if !(CONFIG_REALTIME_ONLY)
extern const TOKENEXTRA vp8_dct_value_tokens[DCT_MAX_VALUE*2];
extern const TOKENEXTRA *vp8_dct_value_tokens_ptr;
extern int vp8_dct_value_cost[DCT_MAX_VALUE*2];
extern int *vp8_dct_value_cost_ptr;

static int cost_coeffs(MACROBLOCK *mb, BLOCKD *b, int type, ENTROPY_CONTEXT *a, ENTROPY_CONTEXT *l)
{
    int c = !type;              /* start at coef 0, unless Y with Y2 */
    int eob = b->eob;
    int pt ;    /* surrounding block/prev coef predictor */
    int cost = 0;
    short *qcoeff_ptr = b->qcoeff;

    VP8_COMBINEENTROPYCONTEXTS(pt, *a, *l);

# define QC( I)  ( qcoeff_ptr [vp8_default_zig_zag1d[I]] )

    for (; c < eob; c++)
    {
        int v = QC(c);
        int t = vp8_dct_value_tokens_ptr[v].Token;
        cost += mb->token_costs [type] [vp8_coef_bands[c]] [pt] [t];
        cost += vp8_dct_value_cost_ptr[v];
        pt = vp8_prev_token_class[t];
    }

# undef QC

    if (c < 16)
        cost += mb->token_costs [type] [vp8_coef_bands[c]] [pt] [DCT_EOB_TOKEN];

    return cost;
}

static int mbycost_coeffs(MACROBLOCK *mb)
{
    int cost = 0;
    int b;
    TEMP_CONTEXT t;
    int type = 0;

    MACROBLOCKD *x = &mb->e_mbd;

    vp8_setup_temp_context(&t, x->above_context[Y1CONTEXT], x->left_context[Y1CONTEXT], 4);

    if (x->mbmi.mode == SPLITMV)
        type = 3;

    for (b = 0; b < 16; b++)
        cost += cost_coeffs(mb, x->block + b, type,
                            t.a + vp8_block2above[b], t.l + vp8_block2left[b]);

    return cost;
}

#define RDFUNC(RM,DM,R,D,target_rd) ( ((128+(R)*(RM)) >> 8) + (DM)*(D) )

void vp8_optimize_b(MACROBLOCK *x, int i, int type, ENTROPY_CONTEXT *a, ENTROPY_CONTEXT *l, const VP8_ENCODER_RTCD *rtcd)
{
    BLOCK *b = &x->block[i];
    BLOCKD *bd = &x->e_mbd.block[i];
    short *dequant_ptr = &bd->dequant[0][0];
    int nzpos[16] = {0};
    short saved_qcoefs[16];
    short saved_dqcoefs[16];
    int baserate, baseerror, baserd;
    int rate, error, thisrd;
    int k;
    int nzcoefcount = 0;
    int nc, bestnc = 0;
    int besteob;

    // count potential coefficient to be optimized
    for (k = !type; k < 16; k++)
    {
        int qcoef = abs(bd->qcoeff[k]);
        int coef = abs(b->coeff[k]);
        int dq   = dequant_ptr[k];

        if (qcoef && (qcoef * dq > coef) && (qcoef * dq < coef + dq))
        {
            nzpos[nzcoefcount] = k;
            nzcoefcount++;
        }
    }

    // if nothing here, do nothing for this block.
    if (!nzcoefcount)
    {
        *a = *l = (bd->eob != !type);
        return;
    }

    // save a copy of quantized coefficients
    vpx_memcpy(saved_qcoefs, bd->qcoeff, 32);
    vpx_memcpy(saved_dqcoefs, bd->dqcoeff, 32);

    besteob   = bd->eob;
    baserate  = cost_coeffs(x, bd, type, a, l);
    baseerror = ENCODEMB_INVOKE(&rtcd->encodemb, berr)(b->coeff, bd->dqcoeff) >> 2;
    baserd    = RDFUNC(x->rdmult, x->rddiv, baserate, baseerror, 100);

    for (nc = 1; nc < (1 << nzcoefcount); nc++)
    {
        //reset coefficients
        vpx_memcpy(bd->qcoeff,  saved_qcoefs,  32);
        vpx_memcpy(bd->dqcoeff, saved_dqcoefs, 32);

        for (k = 0; k < nzcoefcount; k++)
        {
            int pos = nzpos[k];

            if ((nc & (1 << k)))
            {
                int cur_qcoef = bd->qcoeff[pos];

                if (cur_qcoef < 0)
                {
                    bd->qcoeff[pos]++;
                    bd->dqcoeff[pos] = bd->qcoeff[pos] * dequant_ptr[pos];
                }
                else
                {
                    bd->qcoeff[pos]--;
                    bd->dqcoeff[pos] = bd->qcoeff[pos] * dequant_ptr[pos];
                }
            }
        }

        {
            int eob = -1;
            int rc;
            int m;

            for (m = 0; m < 16; m++)
            {
                rc   = vp8_default_zig_zag1d[m];

                if (bd->qcoeff[rc])
                    eob = m;
            }

            bd->eob = eob + 1;
        }

        rate  = cost_coeffs(x, bd, type, a, l);
        error = ENCODEMB_INVOKE(&rtcd->encodemb, berr)(b->coeff, bd->dqcoeff) >> 2;
        thisrd = RDFUNC(x->rdmult, x->rddiv, rate, error, 100);

        if (thisrd < baserd)
        {
            baserd = thisrd;
            bestnc = nc;
            besteob = bd->eob;
        }
    }

    //reset coefficients
    vpx_memcpy(bd->qcoeff,  saved_qcoefs, 32);
    vpx_memcpy(bd->dqcoeff, saved_dqcoefs, 32);

    if (bestnc)
    {
        for (k = 0; k < nzcoefcount; k++)
        {
            int pos = nzpos[k];

            if (bestnc & (1 << k))
            {
                int cur_qcoef = bd->qcoeff[pos];

                if (cur_qcoef < 0)
                {
                    bd->qcoeff[pos]++;
                    bd->dqcoeff[pos] = bd->qcoeff[pos] * dequant_ptr[pos];
                }
                else
                {
                    bd->qcoeff[pos]--;
                    bd->dqcoeff[pos] = bd->qcoeff[pos] * dequant_ptr[pos];
                }
            }
        }

#if 0
        {
            int eob = -1;
            int rc;
            int m;

            for (m = 0; m < 16; m++)
            {
                rc   = vp8_default_zig_zag1d[m];

                if (bd->qcoeff[rc])
                    eob = m;
            }

            bd->eob = eob + 1;
        }
#endif
    }

#if 1
    bd->eob = besteob;
#endif
#if 0
    {
        int eob = -1;
        int rc;
        int m;

        for (m = 0; m < 16; m++)
        {
            rc   = vp8_default_zig_zag1d[m];

            if (bd->qcoeff[rc])
                eob = m;
        }

        bd->eob = eob + 1;
    }

#endif
    *a = *l = (bd->eob != !type);
    return;
}

void vp8_optimize_bplus(MACROBLOCK *x, int i, int type, ENTROPY_CONTEXT *a, ENTROPY_CONTEXT *l, const VP8_ENCODER_RTCD *rtcd)
{
    BLOCK *b = &x->block[i];
    BLOCKD *bd = &x->e_mbd.block[i];
    short *dequant_ptr = &bd->dequant[0][0];
    int nzpos[16] = {0};
    short saved_qcoefs[16];
    short saved_dqcoefs[16];
    int baserate, baseerror, baserd;
    int rate, error, thisrd;
    int k;
    int nzcoefcount = 0;
    int nc, bestnc = 0;
    int besteob;

    // count potential coefficient to be optimized
    for (k = !type; k < 16; k++)
    {
        int qcoef = abs(bd->qcoeff[k]);
        int coef = abs(b->coeff[k]);
        int dq   = dequant_ptr[k];

        if (qcoef && (qcoef * dq < coef) && (coef < (qcoef * dq + dq)))
        {
            nzpos[nzcoefcount] = k;
            nzcoefcount++;
        }
    }

    // if nothing here, do nothing for this block.
    if (!nzcoefcount)
    {
        //do not update context, we need do the other half.
        //*a = *l = (bd->eob != !type);
        return;
    }

    // save a copy of quantized coefficients
    vpx_memcpy(saved_qcoefs, bd->qcoeff, 32);
    vpx_memcpy(saved_dqcoefs, bd->dqcoeff, 32);

    besteob   = bd->eob;
    baserate  = cost_coeffs(x, bd, type, a, l);
    baseerror = ENCODEMB_INVOKE(&rtcd->encodemb, berr)(b->coeff, bd->dqcoeff) >> 2;
    baserd    = RDFUNC(x->rdmult, x->rddiv, baserate, baseerror, 100);

    for (nc = 1; nc < (1 << nzcoefcount); nc++)
    {
        //reset coefficients
        vpx_memcpy(bd->qcoeff, saved_qcoefs, 32);
        vpx_memcpy(bd->dqcoeff, saved_dqcoefs, 32);

        for (k = 0; k < nzcoefcount; k++)
        {
            int pos = nzpos[k];

            if ((nc & (1 << k)))
            {
                int cur_qcoef = bd->qcoeff[pos];

                if (cur_qcoef < 0)
                {
                    bd->qcoeff[pos]--;
                    bd->dqcoeff[pos] = bd->qcoeff[pos] * dequant_ptr[pos];
                }
                else
                {
                    bd->qcoeff[pos]++;
                    bd->dqcoeff[pos] = bd->qcoeff[pos] * dequant_ptr[pos];
                }
            }
        }

        {
            int eob = -1;
            int rc;
            int m;

            for (m = 0; m < 16; m++)
            {
                rc   = vp8_default_zig_zag1d[m];

                if (bd->qcoeff[rc])
                    eob = m;
            }

            bd->eob = eob + 1;
        }

        rate  = cost_coeffs(x, bd, type, a, l);
        error = ENCODEMB_INVOKE(&rtcd->encodemb, berr)(b->coeff, bd->dqcoeff) >> 2;
        thisrd = RDFUNC(x->rdmult, x->rddiv, rate, error, 100);

        if (thisrd < baserd)
        {
            baserd = thisrd;
            bestnc = nc;
            besteob = bd->eob;
        }
    }

    //reset coefficients
    vpx_memcpy(bd->qcoeff,  saved_qcoefs, 32);
    vpx_memcpy(bd->dqcoeff, saved_dqcoefs, 32);

    if (bestnc)
    {
        for (k = 0; k < nzcoefcount; k++)
        {
            int pos = nzpos[k];

            if (bestnc & (1 << k))
            {
                int cur_qcoef = bd->qcoeff[pos];

                if (cur_qcoef < 0)
                {
                    bd->qcoeff[pos]++;
                    bd->dqcoeff[pos] = bd->qcoeff[pos] * dequant_ptr[pos];
                }
                else
                {
                    bd->qcoeff[pos]--;
                    bd->dqcoeff[pos] = bd->qcoeff[pos] * dequant_ptr[pos];
                }
            }
        }
    }

    bd->eob = besteob;
    //do not update context, we need do the other half.
    //*a = *l = (bd->eob != !type);
    return;
}

void vp8_optimize_y2b(MACROBLOCK *x, int i, int type, ENTROPY_CONTEXT *a, ENTROPY_CONTEXT *l, const VP8_ENCODER_RTCD *rtcd)
{

    BLOCK *b = &x->block[i];
    BLOCKD *bd = &x->e_mbd.block[i];
    short *dequant_ptr = &bd->dequant[0][0];

    int baserate, baseerror, baserd;
    int rate, error, thisrd;
    int k;

    if (bd->eob == 0)
        return;

    baserate  = cost_coeffs(x, bd, type, a, l);
    baseerror = ENCODEMB_INVOKE(&rtcd->encodemb, berr)(b->coeff, bd->dqcoeff) >> 4;
    baserd = RDFUNC(x->rdmult, x->rddiv, baserate, baseerror, 100);

    for (k = 0; k < 16; k++)
    {
        int cur_qcoef = bd->qcoeff[k];

        if (!cur_qcoef)
            continue;

        if (cur_qcoef < 0)
        {
            bd->qcoeff[k]++;
            bd->dqcoeff[k] = bd->qcoeff[k] * dequant_ptr[k];
        }
        else
        {
            bd->qcoeff[k]--;
            bd->dqcoeff[k] = bd->qcoeff[k] * dequant_ptr[k];
        }

        if (bd->qcoeff[k] == 0)
        {
            int eob = -1;
            int rc;
            int l;

            for (l = 0; l < 16; l++)
            {
                rc   = vp8_default_zig_zag1d[l];

                if (bd->qcoeff[rc])
                    eob = l;
            }

            bd->eob = eob + 1;
        }

        rate  =   cost_coeffs(x, bd, type, a, l);
        error = ENCODEMB_INVOKE(&rtcd->encodemb, berr)(b->coeff, bd->dqcoeff) >> 4;
        thisrd = RDFUNC(x->rdmult, x->rddiv, rate, error, 100);

        if (thisrd > baserd)
        {
            bd->qcoeff[k] = cur_qcoef;
            bd->dqcoeff[k] = cur_qcoef * dequant_ptr[k];
        }
        else
        {
            baserd = thisrd;
        }

    }

    {
        int eob = -1;
        int rc;

        for (k = 0; k < 16; k++)
        {
            rc   = vp8_default_zig_zag1d[k];

            if (bd->qcoeff[rc])
                eob = k;
        }

        bd->eob = eob + 1;
    }

    return;
}


void vp8_optimize_mb(MACROBLOCK *x, const VP8_ENCODER_RTCD *rtcd)
{
    int cost = 0;
    int b;
    TEMP_CONTEXT t, t2;
    int type = 0;

    vp8_setup_temp_context(&t, x->e_mbd.above_context[Y1CONTEXT], x->e_mbd.left_context[Y1CONTEXT], 4);

    if (x->e_mbd.mbmi.mode == SPLITMV || x->e_mbd.mbmi.mode == B_PRED)
        type = 3;

    for (b = 0; b < 16; b++)
    {
        //vp8_optimize_bplus(x, b, type, t.a + vp8_block2above[b], t.l + vp8_block2left[b]);
        vp8_optimize_b(x, b, type, t.a + vp8_block2above[b], t.l + vp8_block2left[b], rtcd);
    }

    vp8_setup_temp_context(&t, x->e_mbd.above_context[UCONTEXT], x->e_mbd.left_context[UCONTEXT], 2);
    vp8_setup_temp_context(&t2, x->e_mbd.above_context[VCONTEXT], x->e_mbd.left_context[VCONTEXT], 2);

    for (b = 16; b < 20; b++)
    {
        //vp8_optimize_bplus(x, b, vp8_block2type[b], t.a + vp8_block2above[b], t.l + vp8_block2left[b]);
        vp8_optimize_b(x, b, vp8_block2type[b], t.a + vp8_block2above[b], t.l + vp8_block2left[b], rtcd);
    }

    for (b = 20; b < 24; b++)
    {
        //vp8_optimize_bplus(x, b, vp8_block2type[b], t2.a + vp8_block2above[b], t2.l + vp8_block2left[b]);
        vp8_optimize_b(x, b, vp8_block2type[b], t2.a + vp8_block2above[b], t2.l + vp8_block2left[b], rtcd);
    }
}



void vp8_super_slow_yquant_optimization(MACROBLOCK *x, int type, const VP8_ENCODER_RTCD *rtcd)
{
    BLOCK  *b = &x->block[0];
    BLOCKD *bd = &x->e_mbd.block[0];
    short *dequant_ptr = &bd->dequant[0][0];
    struct
    {
        int block;
        int pos;
    } nzpos[256];
    short saved_qcoefs[256];
    short saved_dqcoefs[256];
    short *coef_ptr   = x->coeff;
    short *qcoef_ptr  = x->e_mbd.qcoeff;
    short *dqcoef_ptr = x->e_mbd.dqcoeff;

    int baserate, baseerror, baserd;
    int rate, error, thisrd;
    int i, k;
    int nzcoefcount = 0;
    int nc, bestnc = 0;
    int besteob;

    //this code has assumption in macroblock coeff buffer layout
    for (i = 0; i < 16; i++)
    {
        // count potential coefficient to be optimized
        for (k = !type; k < 16; k++)
        {
            int qcoef = abs(qcoef_ptr[i*16 + k]);
            int coef = abs(coef_ptr[i*16 + k]);
            int dq   = dequant_ptr[k];

            if (qcoef && (qcoef * dq > coef) && (qcoef * dq < coef + dq))
            {
                nzpos[nzcoefcount].block = i;
                nzpos[nzcoefcount].pos   = k;
                nzcoefcount++;
            }
        }
    }

    // if nothing here, do nothing for this macro_block.
    if (!nzcoefcount || nzcoefcount > 15)
    {
        return;
    }

    /******************************************************************************
    looking from each coeffient's perspective, each identifed coefficent above could
    have 2 values:roundeddown(x) and roundedup(x). Therefore the total number of
    different states is less than 2**nzcoefcount.
    ******************************************************************************/
    // save the qunatized coefficents and dequantized coefficicents
    vpx_memcpy(saved_qcoefs, x->e_mbd.qcoeff,  256);
    vpx_memcpy(saved_dqcoefs, x->e_mbd.dqcoeff, 256);

    baserate    = mbycost_coeffs(x);
    baseerror   = ENCODEMB_INVOKE(&rtcd->encodemb, mberr)(x, !type);
    baserd      = RDFUNC(x->rdmult, x->rddiv, baserate, baseerror, 100);

    for (nc = 1; nc < (1 << nzcoefcount); nc++)
    {
        //reset coefficients
        vpx_memcpy(x->e_mbd.qcoeff,  saved_qcoefs, 256);
        vpx_memcpy(x->e_mbd.dqcoeff, saved_dqcoefs, 256);

        for (k = 0; k < nzcoefcount; k++)
        {
            int bk  = nzpos[k].block;
            int pos = nzpos[k].pos;
            int mbkpos  = bk * 16 + pos;

            if ((nc & (1 << k)))
            {
                int cur_qcoef = x->e_mbd.qcoeff[mbkpos];

                if (cur_qcoef < 0)
                {
                    x->e_mbd.qcoeff[mbkpos]++;
                    x->e_mbd.dqcoeff[mbkpos] = x->e_mbd.qcoeff[mbkpos] * dequant_ptr[pos];
                }
                else
                {
                    x->e_mbd.qcoeff[mbkpos]--;
                    x->e_mbd.dqcoeff[mbkpos] = x->e_mbd.qcoeff[mbkpos] * dequant_ptr[pos];
                }
            }
        }

        for (i = 0; i < 16; i++)
        {
            BLOCKD *bd = &x->e_mbd.block[i];
            {
                int eob = -1;
                int rc;
                int l;

                for (l = 0; l < 16; l++)
                {
                    rc   = vp8_default_zig_zag1d[l];

                    if (bd->qcoeff[rc])
                        eob = l;
                }

                bd->eob = eob + 1;
            }
        }

        rate  = mbycost_coeffs(x);
        error = ENCODEMB_INVOKE(&rtcd->encodemb, mberr)(x, !type);;
        thisrd = RDFUNC(x->rdmult, x->rddiv, rate, error, 100);

        if (thisrd < baserd)
        {
            baserd = thisrd;
            bestnc = nc;
            besteob = bd->eob;
        }
    }

    //reset coefficients
    vpx_memcpy(x->e_mbd.qcoeff,  saved_qcoefs, 256);
    vpx_memcpy(x->e_mbd.dqcoeff, saved_dqcoefs, 256);

    if (bestnc)
    {
        for (k = 0; k < nzcoefcount; k++)
        {
            int bk  = nzpos[k].block;
            int pos = nzpos[k].pos;
            int mbkpos  = bk * 16 + pos;

            if ((nc & (1 << k)))
            {
                int cur_qcoef = x->e_mbd.qcoeff[mbkpos];

                if (cur_qcoef < 0)
                {
                    x->e_mbd.qcoeff[mbkpos]++;
                    x->e_mbd.dqcoeff[mbkpos] = x->e_mbd.qcoeff[mbkpos] * dequant_ptr[pos];
                }
                else
                {
                    x->e_mbd.qcoeff[mbkpos]--;
                    x->e_mbd.dqcoeff[mbkpos] = x->e_mbd.qcoeff[mbkpos] * dequant_ptr[pos];
                }
            }
        }
    }

    for (i = 0; i < 16; i++)
    {
        BLOCKD *bd = &x->e_mbd.block[i];
        {
            int eob = -1;
            int rc;
            int l;

            for (l = 0; l < 16; l++)
            {
                rc   = vp8_default_zig_zag1d[l];

                if (bd->qcoeff[rc])
                    eob = l;
            }

            bd->eob = eob + 1;
        }
    }

    return;
}

static void vp8_find_mb_skip_coef(MACROBLOCK *x)
{
    int i;

    x->e_mbd.mbmi.mb_skip_coeff = 1;

    if (x->e_mbd.mbmi.mode != B_PRED && x->e_mbd.mbmi.mode != SPLITMV)
    {
        for (i = 0; i < 16; i++)
        {
            x->e_mbd.mbmi.mb_skip_coeff &= (x->e_mbd.block[i].eob < 2);
        }

        for (i = 16; i < 25; i++)
        {
            x->e_mbd.mbmi.mb_skip_coeff &= (!x->e_mbd.block[i].eob);
        }
    }
    else
    {
        for (i = 0; i < 24; i++)
        {
            x->e_mbd.mbmi.mb_skip_coeff &= (!x->e_mbd.block[i].eob);
        }
    }
}


void vp8_optimize_mb_slow(MACROBLOCK *x, const VP8_ENCODER_RTCD *rtcd)
{
    int cost = 0;
    int b;
    TEMP_CONTEXT t, t2;
    int type = 0;


    vp8_setup_temp_context(&t, x->e_mbd.above_context[Y1CONTEXT], x->e_mbd.left_context[Y1CONTEXT], 4);

    if (x->e_mbd.mbmi.mode == SPLITMV || x->e_mbd.mbmi.mode == B_PRED)
        type = 3;

    vp8_super_slow_yquant_optimization(x, type, rtcd);
    /*
    for(b=0;b<16;b++)
    {
        vp8_optimize_b(x, b, type, t.a + vp8_block2above[b], t.l + vp8_block2left[b]);
    }
    */

    vp8_setup_temp_context(&t, x->e_mbd.above_context[UCONTEXT], x->e_mbd.left_context[UCONTEXT], 2);

    for (b = 16; b < 20; b++)
    {
        vp8_optimize_b(x, b, vp8_block2type[b], t.a + vp8_block2above[b], t.l + vp8_block2left[b], rtcd);
    }

    vp8_setup_temp_context(&t2, x->e_mbd.above_context[VCONTEXT], x->e_mbd.left_context[VCONTEXT], 2);

    for (b = 20; b < 24; b++)
    {
        vp8_optimize_b(x, b, vp8_block2type[b], t2.a + vp8_block2above[b], t2.l + vp8_block2left[b], rtcd);
    }
}


void vp8_optimize_mby(MACROBLOCK *x, const VP8_ENCODER_RTCD *rtcd)
{
    int cost = 0;
    int b;
    TEMP_CONTEXT t;
    int type = 0;

    if (!x->e_mbd.above_context[Y1CONTEXT])
        return;

    if (!x->e_mbd.left_context[Y1CONTEXT])
        return;

    vp8_setup_temp_context(&t, x->e_mbd.above_context[Y1CONTEXT], x->e_mbd.left_context[Y1CONTEXT], 4);

    if (x->e_mbd.mbmi.mode == SPLITMV || x->e_mbd.mbmi.mode == B_PRED)
        type = 3;

    for (b = 0; b < 16; b++)
    {
        vp8_optimize_b(x, b, type, t.a + vp8_block2above[b], t.l + vp8_block2left[b], rtcd);
    }

}

void vp8_optimize_mbuv(MACROBLOCK *x, const VP8_ENCODER_RTCD *rtcd)
{
    int cost = 0;
    int b;
    TEMP_CONTEXT t, t2;
    int type = 0;

    if (!x->e_mbd.above_context[UCONTEXT])
        return;

    if (!x->e_mbd.left_context[UCONTEXT])
        return;

    if (!x->e_mbd.above_context[VCONTEXT])
        return;

    if (!x->e_mbd.left_context[VCONTEXT])
        return;


    vp8_setup_temp_context(&t, x->e_mbd.above_context[UCONTEXT], x->e_mbd.left_context[UCONTEXT], 2);
    vp8_setup_temp_context(&t2, x->e_mbd.above_context[VCONTEXT], x->e_mbd.left_context[VCONTEXT], 2);

    for (b = 16; b < 20; b++)
    {
        vp8_optimize_b(x, b, vp8_block2type[b],
                       t.a + vp8_block2above[b], t.l + vp8_block2left[b], rtcd);

    }

    for (b = 20; b < 24; b++)
    {
        vp8_optimize_b(x, b, vp8_block2type[b],
                       t2.a + vp8_block2above[b], t2.l + vp8_block2left[b], rtcd);
    }

}
#endif

void vp8_encode_inter16x16(const VP8_ENCODER_RTCD *rtcd, MACROBLOCK *x)
{
    vp8_build_inter_predictors_mb(&x->e_mbd);

    vp8_subtract_mb(rtcd, x);

    vp8_transform_mb(x);

    vp8_quantize_mb(x);

#if !(CONFIG_REALTIME_ONLY)
#if 1

    if (x->optimize && x->rddiv > 1)
    {
        vp8_optimize_mb(x, rtcd);
        vp8_find_mb_skip_coef(x);
    }

#endif
#endif

    vp8_inverse_transform_mb(IF_RTCD(&rtcd->common->idct), &x->e_mbd);

    vp8_recon16x16mb(IF_RTCD(&rtcd->common->recon), &x->e_mbd);
}


/* this funciton is used by first pass only */
void vp8_encode_inter16x16y(const VP8_ENCODER_RTCD *rtcd, MACROBLOCK *x)
{
    vp8_build_inter_predictors_mby(&x->e_mbd);

    ENCODEMB_INVOKE(&rtcd->encodemb, submby)(x->src_diff, x->src.y_buffer, x->e_mbd.predictor, x->src.y_stride);

    vp8_transform_mby(x);

    vp8_quantize_mby(x);

    vp8_inverse_transform_mby(IF_RTCD(&rtcd->common->idct), &x->e_mbd);

    vp8_recon16x16mby(IF_RTCD(&rtcd->common->recon), &x->e_mbd);
}


void vp8_encode_inter16x16uv(const VP8_ENCODER_RTCD *rtcd, MACROBLOCK *x)
{
    vp8_build_inter_predictors_mbuv(&x->e_mbd);

    ENCODEMB_INVOKE(&rtcd->encodemb, submbuv)(x->src_diff, x->src.u_buffer, x->src.v_buffer, x->e_mbd.predictor, x->src.uv_stride);

    vp8_transform_mbuv(x);

    vp8_quantize_mbuv(x);

    vp8_inverse_transform_mbuv(IF_RTCD(&rtcd->common->idct), &x->e_mbd);

    vp8_recon_intra_mbuv(IF_RTCD(&rtcd->common->recon), &x->e_mbd);
}


void vp8_encode_inter16x16uvrd(const VP8_ENCODER_RTCD *rtcd, MACROBLOCK *x)
{
    vp8_build_inter_predictors_mbuv(&x->e_mbd);
    ENCODEMB_INVOKE(&rtcd->encodemb, submbuv)(x->src_diff, x->src.u_buffer, x->src.v_buffer, x->e_mbd.predictor, x->src.uv_stride);

    vp8_transform_mbuvrd(x);

    vp8_quantize_mbuvrd(x);

}
