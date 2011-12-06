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
#include "vpx/vpx_integer.h"
#include "recon.h"
#include "subpixel.h"
#include "blockd.h"
#include "reconinter.h"
#if CONFIG_RUNTIME_CPU_DETECT
#include "onyxc_int.h"
#endif

void vp8_copy_mem16x16_c(
    unsigned char *src,
    int src_stride,
    unsigned char *dst,
    int dst_stride)
{

    int r;

    for (r = 0; r < 16; r++)
    {
#if !(CONFIG_FAST_UNALIGNED)
        dst[0] = src[0];
        dst[1] = src[1];
        dst[2] = src[2];
        dst[3] = src[3];
        dst[4] = src[4];
        dst[5] = src[5];
        dst[6] = src[6];
        dst[7] = src[7];
        dst[8] = src[8];
        dst[9] = src[9];
        dst[10] = src[10];
        dst[11] = src[11];
        dst[12] = src[12];
        dst[13] = src[13];
        dst[14] = src[14];
        dst[15] = src[15];

#else
        ((uint32_t *)dst)[0] = ((uint32_t *)src)[0] ;
        ((uint32_t *)dst)[1] = ((uint32_t *)src)[1] ;
        ((uint32_t *)dst)[2] = ((uint32_t *)src)[2] ;
        ((uint32_t *)dst)[3] = ((uint32_t *)src)[3] ;

#endif
        src += src_stride;
        dst += dst_stride;

    }

}

#if CONFIG_DUALPRED
void vp8_avg_mem16x16_c(
    unsigned char *src,
    int src_stride,
    unsigned char *dst,
    int dst_stride)
{
    int r;

    for (r = 0; r < 16; r++)
    {
        int n;

        for (n = 0; n < 16; n++)
        {
            dst[n] = (dst[n] + src[n] + 1) >> 1;
        }

        src += src_stride;
        dst += dst_stride;
    }
}
#endif /* CONFIG_DUALPRED */

void vp8_copy_mem8x8_c(
    unsigned char *src,
    int src_stride,
    unsigned char *dst,
    int dst_stride)
{
    int r;

    for (r = 0; r < 8; r++)
    {
#if !(CONFIG_FAST_UNALIGNED)
        dst[0] = src[0];
        dst[1] = src[1];
        dst[2] = src[2];
        dst[3] = src[3];
        dst[4] = src[4];
        dst[5] = src[5];
        dst[6] = src[6];
        dst[7] = src[7];
#else
        ((uint32_t *)dst)[0] = ((uint32_t *)src)[0] ;
        ((uint32_t *)dst)[1] = ((uint32_t *)src)[1] ;
#endif
        src += src_stride;
        dst += dst_stride;

    }

}

#if CONFIG_DUALPRED
void vp8_avg_mem8x8_c(
    unsigned char *src,
    int src_stride,
    unsigned char *dst,
    int dst_stride)
{
    int r;

    for (r = 0; r < 8; r++)
    {
        int n;

        for (n = 0; n < 8; n++)
        {
            dst[n] = (dst[n] + src[n] + 1) >> 1;
        }

        src += src_stride;
        dst += dst_stride;
    }
}
#endif /* CONFIG_DUALPRED */

void vp8_copy_mem8x4_c(
    unsigned char *src,
    int src_stride,
    unsigned char *dst,
    int dst_stride)
{
    int r;

    for (r = 0; r < 4; r++)
    {
#if !(CONFIG_FAST_UNALIGNED)
        dst[0] = src[0];
        dst[1] = src[1];
        dst[2] = src[2];
        dst[3] = src[3];
        dst[4] = src[4];
        dst[5] = src[5];
        dst[6] = src[6];
        dst[7] = src[7];
#else
        ((uint32_t *)dst)[0] = ((uint32_t *)src)[0] ;
        ((uint32_t *)dst)[1] = ((uint32_t *)src)[1] ;
#endif
        src += src_stride;
        dst += dst_stride;

    }

}



void vp8_build_inter_predictors_b(BLOCKD *d, int pitch, vp8_subpix_fn_t sppf)
{
    int r;
    unsigned char *ptr_base;
    unsigned char *ptr;
    unsigned char *pred_ptr = d->predictor;

    ptr_base = *(d->base_pre);

    if (d->bmi.mv.as_mv.row & 7 || d->bmi.mv.as_mv.col & 7)
    {
        ptr = ptr_base + d->pre + (d->bmi.mv.as_mv.row >> 3) * d->pre_stride + (d->bmi.mv.as_mv.col >> 3);
        sppf(ptr, d->pre_stride, d->bmi.mv.as_mv.col & 7, d->bmi.mv.as_mv.row & 7, pred_ptr, pitch);
    }
    else
    {
        ptr_base += d->pre + (d->bmi.mv.as_mv.row >> 3) * d->pre_stride + (d->bmi.mv.as_mv.col >> 3);
        ptr = ptr_base;

        for (r = 0; r < 4; r++)
        {
#if !(CONFIG_FAST_UNALIGNED)
            pred_ptr[0]  = ptr[0];
            pred_ptr[1]  = ptr[1];
            pred_ptr[2]  = ptr[2];
            pred_ptr[3]  = ptr[3];
#else
            *(uint32_t *)pred_ptr = *(uint32_t *)ptr ;
#endif
            pred_ptr     += pitch;
            ptr         += d->pre_stride;
        }
    }
}

static void build_inter_predictors4b(MACROBLOCKD *x, BLOCKD *d, int pitch)
{
    unsigned char *ptr_base;
    unsigned char *ptr;
    unsigned char *pred_ptr = d->predictor;

    ptr_base = *(d->base_pre);
    ptr = ptr_base + d->pre + (d->bmi.mv.as_mv.row >> 3) * d->pre_stride + (d->bmi.mv.as_mv.col >> 3);

    if (d->bmi.mv.as_mv.row & 7 || d->bmi.mv.as_mv.col & 7)
    {
        x->subpixel_predict8x8(ptr, d->pre_stride, d->bmi.mv.as_mv.col & 7, d->bmi.mv.as_mv.row & 7, pred_ptr, pitch);
    }
    else
    {
        RECON_INVOKE(&x->rtcd->recon, copy8x8)(ptr, d->pre_stride, pred_ptr, pitch);
    }
}

static void build_inter_predictors2b(MACROBLOCKD *x, BLOCKD *d, int pitch)
{
    unsigned char *ptr_base;
    unsigned char *ptr;
    unsigned char *pred_ptr = d->predictor;

    ptr_base = *(d->base_pre);
    ptr = ptr_base + d->pre + (d->bmi.mv.as_mv.row >> 3) * d->pre_stride + (d->bmi.mv.as_mv.col >> 3);

    if (d->bmi.mv.as_mv.row & 7 || d->bmi.mv.as_mv.col & 7)
    {
        x->subpixel_predict8x4(ptr, d->pre_stride, d->bmi.mv.as_mv.col & 7, d->bmi.mv.as_mv.row & 7, pred_ptr, pitch);
    }
    else
    {
        RECON_INVOKE(&x->rtcd->recon, copy8x4)(ptr, d->pre_stride, pred_ptr, pitch);
    }
}


/*encoder only*/
void vp8_build_inter16x16_predictors_mbuv(MACROBLOCKD *x)
{
    unsigned char *uptr, *vptr;
    unsigned char *upred_ptr = &x->predictor[256];
    unsigned char *vpred_ptr = &x->predictor[320];

    int mv_row = x->mode_info_context->mbmi.mv.as_mv.row;
    int mv_col = x->mode_info_context->mbmi.mv.as_mv.col;
    int offset;
    int pre_stride = x->block[16].pre_stride;

    /* calc uv motion vectors */
    if (mv_row < 0)
        mv_row -= 1;
    else
        mv_row += 1;

    if (mv_col < 0)
        mv_col -= 1;
    else
        mv_col += 1;

    mv_row /= 2;
    mv_col /= 2;

    mv_row &= x->fullpixel_mask;
    mv_col &= x->fullpixel_mask;

    offset = (mv_row >> 3) * pre_stride + (mv_col >> 3);
    uptr = x->pre.u_buffer + offset;
    vptr = x->pre.v_buffer + offset;

    if ((mv_row | mv_col) & 7)
    {
        x->subpixel_predict8x8(uptr, pre_stride, mv_col & 7, mv_row & 7, upred_ptr, 8);
        x->subpixel_predict8x8(vptr, pre_stride, mv_col & 7, mv_row & 7, vpred_ptr, 8);
    }
    else
    {
        RECON_INVOKE(&x->rtcd->recon, copy8x8)(uptr, pre_stride, upred_ptr, 8);
        RECON_INVOKE(&x->rtcd->recon, copy8x8)(vptr, pre_stride, vpred_ptr, 8);
    }
}

/*encoder only*/
void vp8_build_inter4x4_predictors_mbuv(MACROBLOCKD *x)
{
    int i, j;

    /* build uv mvs */
    for (i = 0; i < 2; i++)
    {
        for (j = 0; j < 2; j++)
        {
            int yoffset = i * 8 + j * 2;
            int uoffset = 16 + i * 2 + j;
            int voffset = 20 + i * 2 + j;

            int temp;

            temp = x->block[yoffset  ].bmi.mv.as_mv.row
                   + x->block[yoffset+1].bmi.mv.as_mv.row
                   + x->block[yoffset+4].bmi.mv.as_mv.row
                   + x->block[yoffset+5].bmi.mv.as_mv.row;

            if (temp < 0) temp -= 4;
            else temp += 4;

            x->block[uoffset].bmi.mv.as_mv.row = (temp / 8) & x->fullpixel_mask;

            temp = x->block[yoffset  ].bmi.mv.as_mv.col
                   + x->block[yoffset+1].bmi.mv.as_mv.col
                   + x->block[yoffset+4].bmi.mv.as_mv.col
                   + x->block[yoffset+5].bmi.mv.as_mv.col;

            if (temp < 0) temp -= 4;
            else temp += 4;

            x->block[uoffset].bmi.mv.as_mv.col = (temp / 8) & x->fullpixel_mask;

            x->block[voffset].bmi.mv.as_mv.row =
                x->block[uoffset].bmi.mv.as_mv.row ;
            x->block[voffset].bmi.mv.as_mv.col =
                x->block[uoffset].bmi.mv.as_mv.col ;
        }
    }

    for (i = 16; i < 24; i += 2)
    {
        BLOCKD *d0 = &x->block[i];
        BLOCKD *d1 = &x->block[i+1];

        if (d0->bmi.mv.as_int == d1->bmi.mv.as_int)
            build_inter_predictors2b(x, d0, 8);
        else
        {
            vp8_build_inter_predictors_b(d0, 8, x->subpixel_predict);
            vp8_build_inter_predictors_b(d1, 8, x->subpixel_predict);
        }
    }
}


/*encoder only*/
void vp8_build_inter16x16_predictors_mby(MACROBLOCKD *x)
{
    unsigned char *ptr_base;
    unsigned char *ptr;
    unsigned char *pred_ptr = x->predictor;
    int mv_row = x->mode_info_context->mbmi.mv.as_mv.row;
    int mv_col = x->mode_info_context->mbmi.mv.as_mv.col;
    int pre_stride = x->block[0].pre_stride;

    ptr_base = x->pre.y_buffer;
    ptr = ptr_base + (mv_row >> 3) * pre_stride + (mv_col >> 3);

    if ((mv_row | mv_col) & 7)
    {
        x->subpixel_predict16x16(ptr, pre_stride, mv_col & 7, mv_row & 7, pred_ptr, 16);
    }
    else
    {
        RECON_INVOKE(&x->rtcd->recon, copy16x16)(ptr, pre_stride, pred_ptr, 16);
    }
}

void vp8_build_inter16x16_predictors_mb(MACROBLOCKD *x,
                                        unsigned char *dst_y,
                                        unsigned char *dst_u,
                                        unsigned char *dst_v,
                                        int dst_ystride,
                                        int dst_uvstride)
{
    int offset;
    unsigned char *ptr;
    unsigned char *uptr, *vptr;

    int mv_row = x->mode_info_context->mbmi.mv.as_mv.row;
    int mv_col = x->mode_info_context->mbmi.mv.as_mv.col;

    unsigned char *ptr_base = x->pre.y_buffer;
    int pre_stride = x->block[0].pre_stride;

    ptr = ptr_base + (mv_row >> 3) * pre_stride + (mv_col >> 3);

    if ((mv_row | mv_col) & 7)
    {
        x->subpixel_predict16x16(ptr, pre_stride, mv_col & 7, mv_row & 7, dst_y, dst_ystride);
    }
    else
    {
        RECON_INVOKE(&x->rtcd->recon, copy16x16)(ptr, pre_stride, dst_y, dst_ystride);
    }

    /* calc uv motion vectors */
    if (mv_row < 0)
        mv_row -= 1;
    else
        mv_row += 1;

    if (mv_col < 0)
        mv_col -= 1;
    else
        mv_col += 1;

    mv_row /= 2;
    mv_col /= 2;

    mv_row &= x->fullpixel_mask;
    mv_col &= x->fullpixel_mask;

    pre_stride >>= 1;
    offset = (mv_row >> 3) * pre_stride + (mv_col >> 3);
    uptr = x->pre.u_buffer + offset;
    vptr = x->pre.v_buffer + offset;

    if ((mv_row | mv_col) & 7)
    {
        x->subpixel_predict8x8(uptr, pre_stride, mv_col & 7, mv_row & 7, dst_u, dst_uvstride);
        x->subpixel_predict8x8(vptr, pre_stride, mv_col & 7, mv_row & 7, dst_v, dst_uvstride);
    }
    else
    {
        RECON_INVOKE(&x->rtcd->recon, copy8x8)(uptr, pre_stride, dst_u, dst_uvstride);
        RECON_INVOKE(&x->rtcd->recon, copy8x8)(vptr, pre_stride, dst_v, dst_uvstride);
    }

}

#if CONFIG_DUALPRED
/*
 * This function should be called after an initial call to
 * vp8_build_inter16x16_predictors_mb() or _mby()/_mbuv().
 * It will run a second sixtap filter on a (different) ref
 * frame and average the result with the output of the
 * first sixtap filter. The second reference frame is stored
 * in x->second_pre (the reference frame index is in
 * x->mode_info_context->mbmi.second_ref_frame). The second
 * motion vector is x->mode_info_context->mbmi.second_mv.
 *
 * This allows blending prediction from two reference frames
 * which sometimes leads to better prediction than from a
 * single reference framer.
 */
void vp8_build_2nd_inter16x16_predictors_mb(MACROBLOCKD *x,
                                            unsigned char *dst_y,
                                            unsigned char *dst_u,
                                            unsigned char *dst_v,
                                            int dst_ystride,
                                            int dst_uvstride)
{
    int offset;
    unsigned char *ptr;
    unsigned char *uptr, *vptr;

    int mv_row = x->mode_info_context->mbmi.second_mv.as_mv.row;
    int mv_col = x->mode_info_context->mbmi.second_mv.as_mv.col;

    unsigned char *ptr_base = x->second_pre.y_buffer;
    int pre_stride = x->block[0].pre_stride;

    ptr = ptr_base + (mv_row >> 3) * pre_stride + (mv_col >> 3);

    if ((mv_row | mv_col) & 7)
    {
        x->subpixel_predict_avg16x16(ptr, pre_stride, mv_col & 7, mv_row & 7, dst_y, dst_ystride);
    }
    else
    {
        RECON_INVOKE(&x->rtcd->recon, avg16x16)(ptr, pre_stride, dst_y, dst_ystride);
    }

    /* calc uv motion vectors */
    mv_row = (mv_row + (mv_row > 0)) >> 1;
    mv_col = (mv_col + (mv_col > 0)) >> 1;

    mv_row &= x->fullpixel_mask;
    mv_col &= x->fullpixel_mask;

    pre_stride >>= 1;
    offset = (mv_row >> 3) * pre_stride + (mv_col >> 3);
    uptr = x->second_pre.u_buffer + offset;
    vptr = x->second_pre.v_buffer + offset;

    if ((mv_row | mv_col) & 7)
    {
        x->subpixel_predict_avg8x8(uptr, pre_stride, mv_col & 7, mv_row & 7, dst_u, dst_uvstride);
        x->subpixel_predict_avg8x8(vptr, pre_stride, mv_col & 7, mv_row & 7, dst_v, dst_uvstride);
    }
    else
    {
        RECON_INVOKE(&x->rtcd->recon, avg8x8)(uptr, pre_stride, dst_u, dst_uvstride);
        RECON_INVOKE(&x->rtcd->recon, avg8x8)(vptr, pre_stride, dst_v, dst_uvstride);
    }
}
#endif /* CONFIG_DUALPRED */

static void build_inter4x4_predictors_mb(MACROBLOCKD *x)
{
    int i;

    if (x->mode_info_context->mbmi.partitioning < 3)
    {
        x->block[ 0].bmi = x->mode_info_context->bmi[ 0];
        x->block[ 2].bmi = x->mode_info_context->bmi[ 2];
        x->block[ 8].bmi = x->mode_info_context->bmi[ 8];
        x->block[10].bmi = x->mode_info_context->bmi[10];

        build_inter_predictors4b(x, &x->block[ 0], 16);
        build_inter_predictors4b(x, &x->block[ 2], 16);
        build_inter_predictors4b(x, &x->block[ 8], 16);
        build_inter_predictors4b(x, &x->block[10], 16);
    }
    else
    {
        for (i = 0; i < 16; i += 2)
        {
            BLOCKD *d0 = &x->block[i];
            BLOCKD *d1 = &x->block[i+1];

            x->block[i+0].bmi = x->mode_info_context->bmi[i+0];
            x->block[i+1].bmi = x->mode_info_context->bmi[i+1];

            if (d0->bmi.mv.as_int == d1->bmi.mv.as_int)
                build_inter_predictors2b(x, d0, 16);
            else
            {
                vp8_build_inter_predictors_b(d0, 16, x->subpixel_predict);
                vp8_build_inter_predictors_b(d1, 16, x->subpixel_predict);
            }

        }

    }

    for (i = 16; i < 24; i += 2)
    {
        BLOCKD *d0 = &x->block[i];
        BLOCKD *d1 = &x->block[i+1];

        if (d0->bmi.mv.as_int == d1->bmi.mv.as_int)
            build_inter_predictors2b(x, d0, 8);
        else
        {
            vp8_build_inter_predictors_b(d0, 8, x->subpixel_predict);
            vp8_build_inter_predictors_b(d1, 8, x->subpixel_predict);
        }
    }
}

static
void build_4x4uvmvs(MACROBLOCKD *x)
{
    int i, j;

    for (i = 0; i < 2; i++)
    {
        for (j = 0; j < 2; j++)
        {
            int yoffset = i * 8 + j * 2;
            int uoffset = 16 + i * 2 + j;
            int voffset = 20 + i * 2 + j;

            int temp;

            temp = x->mode_info_context->bmi[yoffset + 0].mv.as_mv.row
                 + x->mode_info_context->bmi[yoffset + 1].mv.as_mv.row
                 + x->mode_info_context->bmi[yoffset + 4].mv.as_mv.row
                 + x->mode_info_context->bmi[yoffset + 5].mv.as_mv.row;

            if (temp < 0) temp -= 4;
            else temp += 4;

            x->block[uoffset].bmi.mv.as_mv.row = (temp / 8) & x->fullpixel_mask;

            temp = x->mode_info_context->bmi[yoffset + 0].mv.as_mv.col
                 + x->mode_info_context->bmi[yoffset + 1].mv.as_mv.col
                 + x->mode_info_context->bmi[yoffset + 4].mv.as_mv.col
                 + x->mode_info_context->bmi[yoffset + 5].mv.as_mv.col;

            if (temp < 0) temp -= 4;
            else temp += 4;

            x->block[uoffset].bmi.mv.as_mv.col = (temp / 8) & x->fullpixel_mask;

            x->block[voffset].bmi.mv.as_mv.row =
                x->block[uoffset].bmi.mv.as_mv.row ;
            x->block[voffset].bmi.mv.as_mv.col =
                x->block[uoffset].bmi.mv.as_mv.col ;
        }
    }
}

void vp8_build_inter_predictors_mb(MACROBLOCKD *x)
{
    if (x->mode_info_context->mbmi.mode != SPLITMV)
    {
        vp8_build_inter16x16_predictors_mb(x, x->predictor, &x->predictor[256],
                                           &x->predictor[320], 16, 8);
#if CONFIG_DUALPRED
        if (x->mode_info_context->mbmi.second_ref_frame)
        {
            /* 256 = offset of U plane in Y+U+V buffer;
             * 320 = offset of V plane in Y+U+V buffer.
             * (256=16x16, 320=16x16+8x8). */
            vp8_build_2nd_inter16x16_predictors_mb(x, x->predictor,
                                                   &x->predictor[256],
                                                   &x->predictor[320], 16, 8);
        }
#endif /* CONFIG_DUALPRED */
    }
    else
    {
        build_4x4uvmvs(x);
        build_inter4x4_predictors_mb(x);
    }
}

