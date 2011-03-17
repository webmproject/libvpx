/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */


#include "vpx_ports/config.h"
#include "vp8/common/idct.h"
#include "quantize.h"
#include "vp8/common/reconintra.h"
#include "vp8/common/reconintra4x4.h"
#include "encodemb.h"
#include "vp8/common/invtrans.h"
#include "vp8/common/recon.h"
#include "dct.h"
#include "vp8/common/g_common.h"
#include "encodeintra.h"

#define intra4x4ibias_rate    128
#define intra4x4pbias_rate    256


#if CONFIG_RUNTIME_CPU_DETECT
#define IF_RTCD(x) (x)
#else
#define IF_RTCD(x) NULL
#endif
void vp8_encode_intra4x4block(const VP8_ENCODER_RTCD *rtcd, MACROBLOCK *x, BLOCK *be, BLOCKD *b, int best_mode)
{
    vp8_predict_intra4x4(b, best_mode, b->predictor);

    ENCODEMB_INVOKE(&rtcd->encodemb, subb)(be, b, 16);

    x->vp8_short_fdct4x4(be->src_diff, be->coeff, 32);

    x->quantize_b(be, b);

    vp8_inverse_transform_b(IF_RTCD(&rtcd->common->idct), b, 32);

    RECON_INVOKE(&rtcd->common->recon, recon)(b->predictor, b->diff, *(b->base_dst) + b->dst, b->dst_stride);
}

void vp8_encode_intra4x4mby(const VP8_ENCODER_RTCD *rtcd, MACROBLOCK *mb)
{
    int i;

    MACROBLOCKD *x = &mb->e_mbd;
    vp8_intra_prediction_down_copy(x);

    for (i = 0; i < 16; i++)
    {
        BLOCK *be = &mb->block[i];
        BLOCKD *b = &x->block[i];

        vp8_encode_intra4x4block(rtcd, mb, be, b, b->bmi.mode);
    }

    return;
}

void vp8_encode_intra16x16mby(const VP8_ENCODER_RTCD *rtcd, MACROBLOCK *x)
{
    int b;

    RECON_INVOKE(&rtcd->common->recon, build_intra_predictors_mby)(&x->e_mbd);

    ENCODEMB_INVOKE(&rtcd->encodemb, submby)(x->src_diff, x->src.y_buffer, x->e_mbd.predictor, x->src.y_stride);

    vp8_transform_intra_mby(x);

    vp8_quantize_mby(x);

#if !(CONFIG_REALTIME_ONLY)
#if 1
    if (x->optimize)
        vp8_optimize_mby(x, rtcd);

#endif
#endif

    vp8_inverse_transform_mby(IF_RTCD(&rtcd->common->idct), &x->e_mbd);

    RECON_INVOKE(&rtcd->common->recon, recon_mby)
        (IF_RTCD(&rtcd->common->recon), &x->e_mbd);

    // make sure block modes are set the way we want them for context updates
    for (b = 0; b < 16; b++)
    {
        BLOCKD *d = &x->e_mbd.block[b];

        switch (x->e_mbd.mode_info_context->mbmi.mode)
        {

        case DC_PRED:
            d->bmi.mode = B_DC_PRED;
            break;
        case V_PRED:
            d->bmi.mode = B_VE_PRED;
            break;
        case H_PRED:
            d->bmi.mode = B_HE_PRED;
            break;
        case TM_PRED:
            d->bmi.mode = B_TM_PRED;
            break;
        default:
            d->bmi.mode = B_DC_PRED;
            break;

        }
    }
}

void vp8_encode_intra16x16mbuv(const VP8_ENCODER_RTCD *rtcd, MACROBLOCK *x)
{
    vp8_build_intra_predictors_mbuv(&x->e_mbd);

    ENCODEMB_INVOKE(&rtcd->encodemb, submbuv)(x->src_diff, x->src.u_buffer, x->src.v_buffer, x->e_mbd.predictor, x->src.uv_stride);

    vp8_transform_mbuv(x);

    vp8_quantize_mbuv(x);

#if !(CONFIG_REALTIME_ONLY)
#if 1

    if (x->optimize==2 ||(x->optimize && x->rddiv > 1))
        vp8_optimize_mbuv(x, rtcd);

#endif
#endif

    vp8_inverse_transform_mbuv(IF_RTCD(&rtcd->common->idct), &x->e_mbd);

    vp8_recon_intra_mbuv(IF_RTCD(&rtcd->common->recon), &x->e_mbd);
}

