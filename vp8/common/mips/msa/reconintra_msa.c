/*
 *  Copyright (c) 2015 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include "./vp8_rtcd.h"
#include "vp8/common/blockd.h"
#include "vp8/common/mips/msa/vp8_macros_msa.h"

static void intra_predict_vert_8x8_msa(uint8_t *src, uint8_t *dst,
                                       int32_t dst_stride)
{
    uint64_t out = LD(src);

    SD4(out, out, out, out, dst, dst_stride);
    dst += (4 * dst_stride);
    SD4(out, out, out, out, dst, dst_stride);
}

static void intra_predict_vert_16x16_msa(uint8_t *src, uint8_t *dst,
                                         int32_t dst_stride)
{
    v16u8 out = LD_UB(src);

    ST_UB8(out, out, out, out, out, out, out, out, dst, dst_stride);
    dst += (8 * dst_stride);
    ST_UB8(out, out, out, out, out, out, out, out, dst, dst_stride);
}

static void intra_predict_horiz_8x8_msa(uint8_t *src, int32_t src_stride,
                                        uint8_t *dst, int32_t dst_stride)
{
    uint64_t out0, out1, out2, out3, out4, out5, out6, out7;

    out0 = src[0 * src_stride] * 0x0101010101010101ull;
    out1 = src[1 * src_stride] * 0x0101010101010101ull;
    out2 = src[2 * src_stride] * 0x0101010101010101ull;
    out3 = src[3 * src_stride] * 0x0101010101010101ull;
    out4 = src[4 * src_stride] * 0x0101010101010101ull;
    out5 = src[5 * src_stride] * 0x0101010101010101ull;
    out6 = src[6 * src_stride] * 0x0101010101010101ull;
    out7 = src[7 * src_stride] * 0x0101010101010101ull;

    SD4(out0, out1, out2, out3, dst, dst_stride);
    dst += (4 * dst_stride);
    SD4(out4, out5, out6, out7, dst, dst_stride);
}

static void intra_predict_horiz_16x16_msa(uint8_t *src, int32_t src_stride,
                                          uint8_t *dst, int32_t dst_stride)
{
    uint32_t row;
    uint8_t inp0, inp1, inp2, inp3;
    v16u8 src0, src1, src2, src3;

    for (row = 4; row--;)
    {
        inp0 = src[0];
        src += src_stride;
        inp1 = src[0];
        src += src_stride;
        inp2 = src[0];
        src += src_stride;
        inp3 = src[0];
        src += src_stride;

        src0 = (v16u8)__msa_fill_b(inp0);
        src1 = (v16u8)__msa_fill_b(inp1);
        src2 = (v16u8)__msa_fill_b(inp2);
        src3 = (v16u8)__msa_fill_b(inp3);

        ST_UB4(src0, src1, src2, src3, dst, dst_stride);
        dst += (4 * dst_stride);
    }
}

static void intra_predict_dc_8x8_msa(uint8_t *src_top, uint8_t *src_left,
                                     int32_t src_stride_left,
                                     uint8_t *dst, int32_t dst_stride,
                                     uint8_t is_above, uint8_t is_left)
{
    uint32_t row, addition = 0;
    uint64_t out;
    v16u8 src_above, store;
    v8u16 sum_above;
    v4u32 sum_top;
    v2u64 sum;

    if (is_left && is_above)
    {
        src_above = LD_UB(src_top);

        sum_above = __msa_hadd_u_h(src_above, src_above);
        sum_top = __msa_hadd_u_w(sum_above, sum_above);
        sum = __msa_hadd_u_d(sum_top, sum_top);
        addition = __msa_copy_u_w((v4i32)sum, 0);

        for (row = 0; row < 8; ++row)
        {
            addition += src_left[row * src_stride_left];
        }

        addition = (addition + 8) >> 4;
        store = (v16u8)__msa_fill_b(addition);
    }
    else if (is_left)
    {
        for (row = 0; row < 8; ++row)
        {
            addition += src_left[row * src_stride_left];
        }

        addition = (addition + 4) >> 3;
        store = (v16u8)__msa_fill_b(addition);
    }
    else if (is_above)
    {
        src_above = LD_UB(src_top);

        sum_above = __msa_hadd_u_h(src_above, src_above);
        sum_top = __msa_hadd_u_w(sum_above, sum_above);
        sum = __msa_hadd_u_d(sum_top, sum_top);
        sum = (v2u64)__msa_srari_d((v2i64)sum, 3);
        store = (v16u8)__msa_splati_b((v16i8)sum, 0);
    }
    else
    {
        store = (v16u8)__msa_ldi_b(128);
    }

    out = __msa_copy_u_d((v2i64)store, 0);

    SD4(out, out, out, out, dst, dst_stride);
    dst += (4 * dst_stride);
    SD4(out, out, out, out, dst, dst_stride);
}

static void intra_predict_dc_16x16_msa(uint8_t *src_top, uint8_t *src_left,
                                       int32_t src_stride_left,
                                       uint8_t *dst, int32_t dst_stride,
                                       uint8_t is_above, uint8_t is_left)
{
    uint32_t row;
    uint32_t addition = 0;
    v16u8 src_above, out;
    v8u16 sum_above;
    v4u32 sum_top;
    v2u64 sum;

    if (is_left && is_above)
    {
        src_above = LD_UB(src_top);

        sum_above = __msa_hadd_u_h(src_above, src_above);
        sum_top = __msa_hadd_u_w(sum_above, sum_above);
        sum = __msa_hadd_u_d(sum_top, sum_top);
        sum_top = (v4u32)__msa_pckev_w((v4i32)sum, (v4i32)sum);
        sum = __msa_hadd_u_d(sum_top, sum_top);
        addition = __msa_copy_u_w((v4i32)sum, 0);

        for (row = 0; row < 16; ++row)
        {
            addition += src_left[row * src_stride_left];
        }

        addition = (addition + 16) >> 5;
        out = (v16u8)__msa_fill_b(addition);
    }
    else if (is_left)
    {
        for (row = 0; row < 16; ++row)
        {
            addition += src_left[row * src_stride_left];
        }

        addition = (addition + 8) >> 4;
        out = (v16u8)__msa_fill_b(addition);
    }
    else if (is_above)
    {
        src_above = LD_UB(src_top);

        sum_above = __msa_hadd_u_h(src_above, src_above);
        sum_top = __msa_hadd_u_w(sum_above, sum_above);
        sum = __msa_hadd_u_d(sum_top, sum_top);
        sum_top = (v4u32)__msa_pckev_w((v4i32)sum, (v4i32)sum);
        sum = __msa_hadd_u_d(sum_top, sum_top);
        sum = (v2u64)__msa_srari_d((v2i64)sum, 4);
        out = (v16u8)__msa_splati_b((v16i8)sum, 0);
    }
    else
    {
        out = (v16u8)__msa_ldi_b(128);
    }

    ST_UB8(out, out, out, out, out, out, out, out, dst, dst_stride);
    dst += (8 * dst_stride);
    ST_UB8(out, out, out, out, out, out, out, out, dst, dst_stride);
}

void vp8_build_intra_predictors_mby_s_msa(struct macroblockd *x,
                                          unsigned char *yabove_row,
                                          unsigned char *yleft,
                                          int left_stride,
                                          unsigned char *ypred_ptr,
                                          int y_stride)
{
    uint32_t row, col;
    uint8_t ytop_left = yabove_row[-1];

    switch (x->mode_info_context->mbmi.mode)
    {
        case DC_PRED:
            intra_predict_dc_16x16_msa(yabove_row, yleft, left_stride,
                                       ypred_ptr, y_stride,
                                       x->up_available, x->left_available);
            break;

        case V_PRED:
            intra_predict_vert_16x16_msa(yabove_row, ypred_ptr, y_stride);
            break;

        case H_PRED:
            intra_predict_horiz_16x16_msa(yleft, left_stride, ypred_ptr,
                                          y_stride);
            break;

        case TM_PRED:
            for (row = 0; row < 16; ++row)
            {
                for (col = 0; col < 16; ++col)
                {
                    int pred = yleft[row * left_stride] + yabove_row[col] -
                               ytop_left;

                    if (pred < 0)
                        pred = 0;

                    if (pred > 255)
                        pred = 255;

                    ypred_ptr[col] = pred;
                }

                ypred_ptr += y_stride;
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

void vp8_build_intra_predictors_mbuv_s_msa(struct macroblockd *x,
                                           unsigned char *uabove_row,
                                           unsigned char *vabove_row,
                                           unsigned char *uleft,
                                           unsigned char *vleft,
                                           int left_stride,
                                           unsigned char *upred_ptr,
                                           unsigned char *vpred_ptr,
                                           int pred_stride)
{
    uint32_t row, col;
    uint8_t utop_left = uabove_row[-1];
    uint8_t vtop_left = vabove_row[-1];

    switch (x->mode_info_context->mbmi.uv_mode)
    {
        case DC_PRED:
            intra_predict_dc_8x8_msa(uabove_row, uleft, left_stride,
                                     upred_ptr, pred_stride,
                                     x->up_available, x->left_available);
            intra_predict_dc_8x8_msa(vabove_row, vleft, left_stride,
                                     vpred_ptr, pred_stride,
                                     x->up_available, x->left_available);
            break;

        case V_PRED:
            intra_predict_vert_8x8_msa(uabove_row, upred_ptr, pred_stride);
            intra_predict_vert_8x8_msa(vabove_row, vpred_ptr, pred_stride);
            break;

        case H_PRED:
            intra_predict_horiz_8x8_msa(uleft, left_stride, upred_ptr,
                                        pred_stride);
            intra_predict_horiz_8x8_msa(vleft, left_stride, vpred_ptr,
                                        pred_stride);
            break;

        case TM_PRED:
            for (row = 0; row < 8; ++row)
            {
                for (col = 0; col < 8; ++col)
                {
                    int predu = uleft[row * left_stride] + uabove_row[col] -
                                utop_left;
                    int predv = vleft[row * left_stride] + vabove_row[col] -
                                vtop_left;

                    if (predu < 0)
                        predu = 0;

                    if (predu > 255)
                        predu = 255;

                    if (predv < 0)
                        predv = 0;

                    if (predv > 255)
                        predv = 255;

                    upred_ptr[col] = predu;
                    vpred_ptr[col] = predv;
                }

                upred_ptr += pred_stride;
                vpred_ptr += pred_stride;
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
