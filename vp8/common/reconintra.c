/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */


#include "./vpx_config.h"
#include "./vpx_dsp_rtcd.h"
#include "./vp8_rtcd.h"
#include "vpx_mem/vpx_mem.h"
#include "vpx_ports/vpx_once.h"
#include "blockd.h"
#include "vp8/common/reconintra.h"

typedef void (*intra_pred_fn)(uint8_t *dst, ptrdiff_t stride,
                              const uint8_t *above, const uint8_t *left);

static intra_pred_fn pred[4];
static intra_pred_fn dc_pred[2][2];

static void vp8_init_intra_predictors_internal(void)
{
    pred[V_PRED] = vpx_v_predictor_16x16;
    pred[H_PRED] = vpx_h_predictor_16x16;
    pred[TM_PRED] = vpx_tm_predictor_16x16;

    dc_pred[0][0] = vpx_dc_128_predictor_16x16;
    dc_pred[0][1] = vpx_dc_top_predictor_16x16;
    dc_pred[1][0] = vpx_dc_left_predictor_16x16;
    dc_pred[1][1] = vpx_dc_predictor_16x16;
}

void vp8_build_intra_predictors_mby_s(MACROBLOCKD *x,
                                      unsigned char * yabove_row,
                                      unsigned char * yleft,
                                      int left_stride,
                                      unsigned char * ypred_ptr,
                                      int y_stride)
{
    MB_PREDICTION_MODE mode = x->mode_info_context->mbmi.mode;
    unsigned char yleft_col[16];
    int i;

    for (i = 0; i < 16; i++)
    {
        yleft_col[i] = yleft[i* left_stride];
    }

    if (mode == DC_PRED)
    {
        dc_pred[x->left_available][x->up_available](ypred_ptr, y_stride,
                                                    yabove_row, yleft_col);
    }
    else
    {
        pred[mode](ypred_ptr, y_stride, yabove_row, yleft_col);
    }
}

void vp8_build_intra_predictors_mbuv_s_c(MACROBLOCKD *x,
                                         unsigned char * uabove_row,
                                         unsigned char * vabove_row,
                                         unsigned char * uleft,
                                         unsigned char * vleft,
                                         int left_stride,
                                         unsigned char * upred_ptr,
                                         unsigned char * vpred_ptr,
                                         int pred_stride)
{
    unsigned char uleft_col[8];
    unsigned char utop_left = uabove_row[-1];
    unsigned char vleft_col[8];
    unsigned char vtop_left = vabove_row[-1];

    int i, j;

    for (i = 0; i < 8; i++)
    {
        uleft_col[i] = uleft [i* left_stride];
        vleft_col[i] = vleft [i* left_stride];
    }

    switch (x->mode_info_context->mbmi.uv_mode)
    {
    case DC_PRED:
    {
        int expected_udc;
        int expected_vdc;
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


        /*memset(upred_ptr,expected_udc,64);*/
        /*memset(vpred_ptr,expected_vdc,64);*/
        for (i = 0; i < 8; i++)
        {
            memset(upred_ptr, expected_udc, 8);
            memset(vpred_ptr, expected_vdc, 8);
            upred_ptr += pred_stride;
            vpred_ptr += pred_stride;
        }
    }
    break;
    case V_PRED:
    {
        for (i = 0; i < 8; i++)
        {
            memcpy(upred_ptr, uabove_row, 8);
            memcpy(vpred_ptr, vabove_row, 8);
            upred_ptr += pred_stride;
            vpred_ptr += pred_stride;
        }

    }
    break;
    case H_PRED:
    {
        for (i = 0; i < 8; i++)
        {
            memset(upred_ptr, uleft_col[i], 8);
            memset(vpred_ptr, vleft_col[i], 8);
            upred_ptr += pred_stride;
            vpred_ptr += pred_stride;
        }
    }

    break;
    case TM_PRED:
    {
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

            upred_ptr += pred_stride;
            vpred_ptr += pred_stride;
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

void vp8_init_intra_predictors(void)
{
    once(vp8_init_intra_predictors_internal);
}
