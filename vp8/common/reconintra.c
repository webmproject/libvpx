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
#include "recon.h"
#include "reconintra.h"
#include "vpx_mem/vpx_mem.h"

/* For skip_recon_mb(), add vp8_build_intra_predictors_mby_s(MACROBLOCKD *x) and
 * vp8_build_intra_predictors_mbuv_s(MACROBLOCKD *x).
 */
void vp8_recon_intra_mbuv(const vp8_recon_rtcd_vtable_t *rtcd, MACROBLOCKD *x)
{
    int i;

    for (i = 16; i < 24; i += 2)
    {
        BLOCKD *b = &x->block[i];
        RECON_INVOKE(rtcd, recon2)(b->predictor, b->diff, *(b->base_dst) + b->dst, b->dst_stride);
    }
}

void vp8_build_intra_predictors_mby_internal(MACROBLOCKD *x, unsigned char *ypred_ptr, int y_stride, int mode)
{

    unsigned char *yabove_row = x->dst.y_buffer - x->dst.y_stride;
    unsigned char yleft_col[16];
    unsigned char ytop_left = yabove_row[-1];
    int r, c, i;

    for (i = 0; i < 16; i++)
    {
        yleft_col[i] = x->dst.y_buffer [i* x->dst.y_stride -1];
    }

    /* for Y */
    switch (mode)
    {
    case DC_PRED:
    {
        int expected_dc;
        int i;
        int shift;
        int average = 0;


        if (x->up_available || x->left_available)
        {
            if (x->up_available)
            {
                for (i = 0; i < 16; i++)
                {
                    average += yabove_row[i];
                }
            }

            if (x->left_available)
            {

                for (i = 0; i < 16; i++)
                {
                    average += yleft_col[i];
                }
            }
            shift = 3 + x->up_available + x->left_available;
            expected_dc = (average + (1 << (shift - 1))) >> shift;
        }
        else
        {
            expected_dc = 128;
        }

        for (r = 0; r < 16; r++)
        {
            vpx_memset(ypred_ptr, expected_dc, 16);
            ypred_ptr += y_stride; /*16;*/
        }
    }
    break;
    case V_PRED:
    {

        for (r = 0; r < 16; r++)
        {

            ((int *)ypred_ptr)[0] = ((int *)yabove_row)[0];
            ((int *)ypred_ptr)[1] = ((int *)yabove_row)[1];
            ((int *)ypred_ptr)[2] = ((int *)yabove_row)[2];
            ((int *)ypred_ptr)[3] = ((int *)yabove_row)[3];
            ypred_ptr += y_stride;
        }
    }
    break;
    case H_PRED:
    {

        for (r = 0; r < 16; r++)
        {

            vpx_memset(ypred_ptr, yleft_col[r], 16);
            ypred_ptr += y_stride;
        }

    }
    break;
    case TM_PRED:
    {

        for (r = 0; r < 16; r++)
        {
            for (c = 0; c < 16; c++)
            {
                int pred =  yleft_col[r] + yabove_row[ c] - ytop_left;

                if (pred < 0)
                    pred = 0;

                if (pred > 255)
                    pred = 255;

                ypred_ptr[c] = pred;
            }

            ypred_ptr += y_stride;
        }

    }
    break;
#if CONIFG_I8X8
    case I8X8_PRED:
#endif
    case B_PRED:
    case NEARESTMV:
    case NEARMV:
    case ZEROMV:
    case NEWMV:
    case SPLITMV:
    case MB_MODE_COUNT:
        break;
    }
}

void vp8_build_intra_predictors_mby(MACROBLOCKD *x)
{
    vp8_build_intra_predictors_mby_internal(x, x->predictor, 16,
                                            x->mode_info_context->mbmi.mode);
}

void vp8_build_intra_predictors_mby_s(MACROBLOCKD *x)
{
    vp8_build_intra_predictors_mby_internal(x, x->dst.y_buffer, x->dst.y_stride,
                                            x->mode_info_context->mbmi.mode);
}

#if CONFIG_COMP_INTRA_PRED
void vp8_build_comp_intra_predictors_mby(MACROBLOCKD *x)
{
    unsigned char predictor[2][256];
    int i;

    vp8_build_intra_predictors_mby_internal(x, predictor[0], 16,
                                            x->mode_info_context->mbmi.mode);
    vp8_build_intra_predictors_mby_internal(x, predictor[1], 16,
                                            x->mode_info_context->mbmi.second_mode);

    for (i = 0; i < 256; i++)
    {
        x->predictor[i] = (predictor[0][i] + predictor[1][i] + 1) >> 1;
    }
}
#endif

void vp8_build_intra_predictors_mbuv_internal(MACROBLOCKD *x,
                                              unsigned char *upred_ptr,
                                              unsigned char *vpred_ptr,
                                              int uv_stride,
                                              int mode)
{
    unsigned char *uabove_row = x->dst.u_buffer - x->dst.uv_stride;
    unsigned char uleft_col[16];
    unsigned char utop_left = uabove_row[-1];
    unsigned char *vabove_row = x->dst.v_buffer - x->dst.uv_stride;
    unsigned char vleft_col[20];
    unsigned char vtop_left = vabove_row[-1];

    int i, j;

    for (i = 0; i < 8; i++)
    {
        uleft_col[i] = x->dst.u_buffer [i* x->dst.uv_stride -1];
        vleft_col[i] = x->dst.v_buffer [i* x->dst.uv_stride -1];
    }

    switch (mode)
    {
    case DC_PRED:
    {
        int expected_udc;
        int expected_vdc;
        int i;
        int shift;
        int Uaverage = 0;
        int Vaverage = 0;

        if (x->up_available)
        {
            for (i = 0; i < 8; i++)
            {
                Uaverage += uabove_row[i];
                Vaverage += vabove_row[i];
            }
        }

        if (x->left_available)
        {
            for (i = 0; i < 8; i++)
            {
                Uaverage += uleft_col[i];
                Vaverage += vleft_col[i];
            }
        }

        if (!x->up_available && !x->left_available)
        {
            expected_udc = 128;
            expected_vdc = 128;
        }
        else
        {
            shift = 2 + x->up_available + x->left_available;
            expected_udc = (Uaverage + (1 << (shift - 1))) >> shift;
            expected_vdc = (Vaverage + (1 << (shift - 1))) >> shift;
        }


        /*vpx_memset(upred_ptr,expected_udc,64);*/
        /*vpx_memset(vpred_ptr,expected_vdc,64);*/
        for (i = 0; i < 8; i++)
        {
            vpx_memset(upred_ptr, expected_udc, 8);
            vpx_memset(vpred_ptr, expected_vdc, 8);
            upred_ptr += uv_stride; /*8;*/
            vpred_ptr += uv_stride; /*8;*/
        }
    }
    break;
    case V_PRED:
    {
        int i;

        for (i = 0; i < 8; i++)
        {
            vpx_memcpy(upred_ptr, uabove_row, 8);
            vpx_memcpy(vpred_ptr, vabove_row, 8);
            upred_ptr += uv_stride; /*8;*/
            vpred_ptr += uv_stride; /*8;*/
        }

    }
    break;
    case H_PRED:
    {
        int i;

        for (i = 0; i < 8; i++)
        {
            vpx_memset(upred_ptr, uleft_col[i], 8);
            vpx_memset(vpred_ptr, vleft_col[i], 8);
            upred_ptr += uv_stride; /*8;*/
            vpred_ptr += uv_stride; /*8;*/
        }
    }

    break;
    case TM_PRED:
    {
        int i;

        for (i = 0; i < 8; i++)
        {
            for (j = 0; j < 8; j++)
            {
                int predu = uleft_col[i] + uabove_row[j] - utop_left;
                int predv = vleft_col[i] + vabove_row[j] - vtop_left;

                if (predu < 0)
                    predu = 0;

                if (predu > 255)
                    predu = 255;

                if (predv < 0)
                    predv = 0;

                if (predv > 255)
                    predv = 255;

                upred_ptr[j] = predu;
                vpred_ptr[j] = predv;
            }

            upred_ptr += uv_stride; /*8;*/
            vpred_ptr += uv_stride; /*8;*/
        }

    }
    break;
    case B_PRED:
    case NEARESTMV:
    case NEARMV:
    case ZEROMV:
    case NEWMV:
    case SPLITMV:
    case MB_MODE_COUNT:
        break;
    }
}

void vp8_build_intra_predictors_mbuv(MACROBLOCKD *x)
{
    vp8_build_intra_predictors_mbuv_internal(x,
                                             &x->predictor[256],
                                             &x->predictor[320],
                                             8,
                                             x->mode_info_context->mbmi.uv_mode);
}

void vp8_build_intra_predictors_mbuv_s(MACROBLOCKD *x)
{
    vp8_build_intra_predictors_mbuv_internal(x,
                                             x->dst.u_buffer,
                                             x->dst.v_buffer,
                                             x->dst.uv_stride,
                                             x->mode_info_context->mbmi.uv_mode);
}

#if CONFIG_COMP_INTRA_PRED
void vp8_build_comp_intra_predictors_mbuv(MACROBLOCKD *x)
{
    unsigned char predictor[2][2][64];
    int i;

    vp8_build_intra_predictors_mbuv_internal(x, predictor[0][0], predictor[1][0], 8,
                                             x->mode_info_context->mbmi.uv_mode);
    vp8_build_intra_predictors_mbuv_internal(x, predictor[0][1], predictor[1][1], 8,
                                             x->mode_info_context->mbmi.second_uv_mode);
    for (i = 0; i < 64; i++)
    {
        x->predictor[256 + i] = (predictor[0][0][i] + predictor[0][1][i] + 1) >> 1;
        x->predictor[256 + 64 + i] = (predictor[1][0][i] + predictor[1][1][i] + 1) >> 1;
    }
}
#endif

void vp8_intra8x8_predict(BLOCKD *x,
                          int mode,
                          unsigned char *predictor)
{

    unsigned char *yabove_row = *(x->base_dst) + x->dst - x->dst_stride;
    unsigned char yleft_col[8];
    unsigned char ytop_left = yabove_row[-1];
    int r, c, i;

    for (i = 0; i < 8; i++)
    {
        yleft_col[i] = (*(x->base_dst))[x->dst - 1 + i * x->dst_stride];
    }
    switch (mode)
    {
    case DC_PRED:
        {
            int expected_dc = 0;

            for (i = 0; i < 8; i++)
            {
                expected_dc += yabove_row[i];
                expected_dc += yleft_col[i];
            }
            expected_dc = (expected_dc + 8) >> 4;

            for (r = 0; r < 8; r++)
            {
                for (c = 0; c < 8; c++)
                {
                    predictor[c] = expected_dc;
                }
                predictor += 16;
            }
        }
        break;
    case V_PRED:
        {
            for (r = 0; r < 8; r++)
            {
                for (c = 0; c < 8; c++)
                {

                    predictor[c] = yabove_row[c];
                }
                predictor += 16;
            }

        }
        break;
    case H_PRED:
        {

            for (r = 0; r < 8; r++)
            {
                for (c = 0; c < 8; c++)
                {
                    predictor[c] = yleft_col[r];
                }
                predictor += 16;
            }
        }
        break;
    case TM_PRED:
        {
            /* prediction similar to true_motion prediction */
            for (r = 0; r < 8; r++)
            {
                for (c = 0; c < 8; c++)
                {
                    int pred = yabove_row[c] - ytop_left + yleft_col[r];
                    if (pred < 0)
                        pred = 0;

                    if (pred > 255)
                        pred = 255;
                    predictor[c] = pred;
                }

                predictor += 16;
            }
        }
        break;
    }
}

#if CONFIG_COMP_INTRA_PRED
void vp8_comp_intra8x8_predict(BLOCKD *x,
                               int mode, int second_mode,
                               unsigned char *out_predictor)
{
    
    unsigned char predictor[2][8*16];
    int i, j;

    vp8_intra8x8_predict(x, mode, predictor[0]);
    vp8_intra8x8_predict(x, second_mode, predictor[1]);

    for (i = 0; i < 8*16; i += 16)
    {
        for (j = i; j < i + 8; j++)
        {
            out_predictor[j] = (predictor[0][j] + predictor[1][j] + 1) >> 1;
        }
    }
}
#endif

void vp8_intra_uv4x4_predict(BLOCKD *x,
                             int mode,
                             unsigned char *predictor)
{

    unsigned char *above_row = *(x->base_dst) + x->dst - x->dst_stride;
    unsigned char left_col[4];
    unsigned char top_left = above_row[-1];
    int r, c, i;

    for (i = 0; i < 4; i++)
    {
        left_col[i] = (*(x->base_dst))[x->dst - 1 + i * x->dst_stride];
    }
    switch (mode)
    {
    case DC_PRED:
        {
            int expected_dc = 0;

            for (i = 0; i < 4; i++)
            {
                expected_dc += above_row[i];
                expected_dc += left_col[i];
            }
            expected_dc = (expected_dc + 4) >> 3;

            for (r = 0; r < 4; r++)
            {
                for (c = 0; c < 4; c++)
                {
                    predictor[c] = expected_dc;
                }
                predictor += 8;
            }
        }
        break;
    case V_PRED:
        {
            for (r = 0; r < 4; r++)
            {
                for (c = 0; c < 4; c++)
                {

                    predictor[c] = above_row[c];
                }
                predictor += 8;
            }

        }
        break;
    case H_PRED:
        {

            for (r = 0; r < 4; r++)
            {
                for (c = 0; c < 4; c++)
                {
                    predictor[c] = left_col[r];
                }
                predictor += 8;
            }
        }
        break;
    case TM_PRED:
        {
            /* prediction similar to true_motion prediction */
            for (r = 0; r < 4; r++)
            {
                for (c = 0; c < 4; c++)
                {
                    int pred = above_row[c] - top_left + left_col[r];
                    if (pred < 0)
                        pred = 0;

                    if (pred > 255)
                        pred = 255;
                    predictor[c] = pred;
                }

                predictor += 8;
            }
        }
        break;
    }
}

#if CONFIG_COMP_INTRA_PRED
void vp8_comp_intra_uv4x4_predict(BLOCKD *x,
                                  int mode, int mode2,
                                  unsigned char *out_predictor)
{
    unsigned char predictor[2][8*4];
    int i, j;

    vp8_intra_uv4x4_predict(x, mode, predictor[0]);
    vp8_intra_uv4x4_predict(x, mode2, predictor[1]);

    for (i = 0; i < 4*8; i += 8)
    {
        for (j = i; j < i + 4; j++)
        {
            out_predictor[j] = (predictor[0][j] + predictor[1][j] + 1) >> 1;
        }
    }
}
#endif

/* TODO: try different ways of use Y-UV mode correlation
 Current code assumes that a uv 4x4 block use same mode
 as corresponding Y 8x8 area
 */
