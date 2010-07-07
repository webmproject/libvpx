/*
 *  Copyright (c) 2010 The VP8 project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */
#include "dixie.h"
#include "idct_add.h"


void
vp8_dixie_walsh(const short *input, short *output)
{
    int i;
    int a1, b1, c1, d1;
    int a2, b2, c2, d2;
    const short *ip = input;
    short *op = output;

    for (i = 0; i < 4; i++)
    {
        a1 = ip[0] + ip[12];
        b1 = ip[4] + ip[8];
        c1 = ip[4] - ip[8];
        d1 = ip[0] - ip[12];

        op[0] = a1 + b1;
        op[4] = c1 + d1;
        op[8] = a1 - b1;
        op[12] = d1 - c1;
        ip++;
        op++;
    }

    ip = output;
    op = output;

    for (i = 0; i < 4; i++)
    {
        a1 = ip[0] + ip[3];
        b1 = ip[1] + ip[2];
        c1 = ip[1] - ip[2];
        d1 = ip[0] - ip[3];

        a2 = a1 + b1;
        b2 = c1 + d1;
        c2 = a1 - b1;
        d2 = d1 - c1;

        op[0] = (a2 + 3) >> 3;
        op[1] = (b2 + 3) >> 3;
        op[2] = (c2 + 3) >> 3;
        op[3] = (d2 + 3) >> 3;

        ip += 4;
        op += 4;
    }
}


static void
walsh1(short *input, short *output)
{
    int i;
    int a1;
    short *op = output;

    a1 = ((input[0] + 3) >> 3);

    for (i = 0; i < 4; i++)
    {
        op[0] = a1;
        op[1] = a1;
        op[2] = a1;
        op[3] = a1;
        op += 4;
    }
}


#define cospi8sqrt2minus1 20091
#define sinpi8sqrt2       35468
#define rounding          0
static void
idct_columns(const short *input, short *output)
{
    int i;
    int a1, b1, c1, d1;

    const short *ip = input;
    short *op = output;
    int temp1, temp2;
    int shortpitch = 4;

    for (i = 0; i < 4; i++)
    {
        a1 = ip[0] + ip[8];
        b1 = ip[0] - ip[8];

        temp1 = (ip[4] * sinpi8sqrt2 + rounding) >> 16;
        temp2 = ip[12] + ((ip[12] * cospi8sqrt2minus1 + rounding) >> 16);
        c1 = temp1 - temp2;

        temp1 = ip[4] + ((ip[4] * cospi8sqrt2minus1 + rounding) >> 16);
        temp2 = (ip[12] * sinpi8sqrt2 + rounding) >> 16;
        d1 = temp1 + temp2;

        op[shortpitch*0] = a1 + d1;
        op[shortpitch*3] = a1 - d1;

        op[shortpitch*1] = b1 + c1;
        op[shortpitch*2] = b1 - c1;

        ip++;
        op++;
    }
}


static void
idct_dc(unsigned char        *predict,
        int                   stride,
        const short          *coeffs,
        int                   pixel_offset)
{
    int i;
    int a1, b1, c1, d1, temp1, temp2;
    short tmp[16];
    unsigned char dc = predict[pixel_offset];

    idct_columns(coeffs, tmp);
    coeffs = tmp;

    for (i = 0; i < 4; i++)
    {
        a1 = coeffs[0] + coeffs[2];
        b1 = coeffs[0] - coeffs[2];

        temp1 = (coeffs[1] * sinpi8sqrt2 + rounding) >> 16;
        temp2 = coeffs[3] + ((coeffs[3] * cospi8sqrt2minus1 + rounding) >> 16);
        c1 = temp1 - temp2;

        temp1 = coeffs[1] + ((coeffs[1] * cospi8sqrt2minus1 + rounding) >> 16);
        temp2 = (coeffs[3] * sinpi8sqrt2 + rounding) >> 16;
        d1 = temp1 + temp2;

        predict[0] = CLAMP_255(dc + ((a1 + d1 + 4) >> 3));
        predict[3] = CLAMP_255(dc + ((a1 - d1 + 4) >> 3));
        predict[1] = CLAMP_255(dc + ((b1 + c1 + 4) >> 3));
        predict[2] = CLAMP_255(dc + ((b1 - c1 + 4) >> 3));

        coeffs += 4;
        predict += stride;
    }
}


static void
idct_row(unsigned char        *predict,
         int                   stride,
         const short          *coeffs,
         int                   pixel_offset)
{
    int i;
    int a1, b1, c1, d1, temp1, temp2;
    short tmp[16];
    unsigned char *recon;

    idct_columns(coeffs, tmp);
    coeffs = tmp;
    recon = predict;
    predict += pixel_offset * stride;

    for (i = 0; i < 4; i++)
    {
        a1 = coeffs[0] + coeffs[2];
        b1 = coeffs[0] - coeffs[2];

        temp1 = (coeffs[1] * sinpi8sqrt2 + rounding) >> 16;
        temp2 = coeffs[3] + ((coeffs[3] * cospi8sqrt2minus1 + rounding) >> 16);
        c1 = temp1 - temp2;

        temp1 = coeffs[1] + ((coeffs[1] * cospi8sqrt2minus1 + rounding) >> 16);
        temp2 = (coeffs[3] * sinpi8sqrt2 + rounding) >> 16;
        d1 = temp1 + temp2;

        recon[0] = CLAMP_255(predict[0] + ((a1 + d1 + 4) >> 3));
        recon[3] = CLAMP_255(predict[3] + ((a1 - d1 + 4) >> 3));
        recon[1] = CLAMP_255(predict[1] + ((b1 + c1 + 4) >> 3));
        recon[2] = CLAMP_255(predict[2] + ((b1 - c1 + 4) >> 3));

        coeffs += 4;
        recon += stride;
    }
}


static void
idct_col(unsigned char        *predict,
         int                   stride,
         const short          *coeffs,
         int                   pixel_offset)
{
    int i;
    int a1, b1, c1, d1, temp1, temp2;
    short tmp[16];

    idct_columns(coeffs, tmp);
    coeffs = tmp;

    for (i = 0; i < 4; i++)
    {
        unsigned char predictor = predict[pixel_offset];
        a1 = coeffs[0] + coeffs[2];
        b1 = coeffs[0] - coeffs[2];

        temp1 = (coeffs[1] * sinpi8sqrt2 + rounding) >> 16;
        temp2 = coeffs[3] + ((coeffs[3] * cospi8sqrt2minus1 + rounding) >> 16);
        c1 = temp1 - temp2;

        temp1 = coeffs[1] + ((coeffs[1] * cospi8sqrt2minus1 + rounding) >> 16);
        temp2 = (coeffs[3] * sinpi8sqrt2 + rounding) >> 16;
        d1 = temp1 + temp2;

        predict[0] = CLAMP_255(predictor + ((a1 + d1 + 4) >> 3));
        predict[1] = CLAMP_255(predictor + ((b1 + c1 + 4) >> 3));
        predict[2] = CLAMP_255(predictor + ((b1 - c1 + 4) >> 3));
        predict[3] = CLAMP_255(predictor + ((a1 - d1 + 4) >> 3));

        coeffs += 4;
        predict += stride;
    }
}


static void
idct_block(unsigned char        *predict,
           int                   stride,
           const short          *coeffs)
{
    int i;
    int a1, b1, c1, d1, temp1, temp2;
    short tmp[16];

    idct_columns(coeffs, tmp);
    coeffs = tmp;

    for (i = 0; i < 4; i++)
    {
        a1 = coeffs[0] + coeffs[2];
        b1 = coeffs[0] - coeffs[2];

        temp1 = (coeffs[1] * sinpi8sqrt2 + rounding) >> 16;
        temp2 = coeffs[3] + ((coeffs[3] * cospi8sqrt2minus1 + rounding) >> 16);
        c1 = temp1 - temp2;

        temp1 = coeffs[1] + ((coeffs[1] * cospi8sqrt2minus1 + rounding) >> 16);
        temp2 = (coeffs[3] * sinpi8sqrt2 + rounding) >> 16;
        d1 = temp1 + temp2;

        predict[0] = CLAMP_255(predict[0] + ((a1 + d1 + 4) >> 3));
        predict[3] = CLAMP_255(predict[3] + ((a1 - d1 + 4) >> 3));
        predict[1] = CLAMP_255(predict[1] + ((b1 + c1 + 4) >> 3));
        predict[2] = CLAMP_255(predict[2] + ((b1 - c1 + 4) >> 3));

        coeffs += 4;
        predict += stride;
    }
}


static void
fill_dc_dc(unsigned char        *predict,
           int                   stride,
           int                   dc,
           int                   pixel_offset)
{
    int i;

    dc = CLAMP_255(dc + predict[pixel_offset]);

    for (i = 0; i < 4; i++)
    {
        predict[0] = predict[1] = predict[2] = predict[3] = dc;
        predict += stride;
    }
}


static void
fill_dc_row(unsigned char        *predict,
            int                   stride,
            int                   dc,
            int                   pixel_offset)
{
    unsigned char *recon;
    int            i;

    recon = predict;
    predict += pixel_offset * stride;
    recon[0] = CLAMP_255(predict[0] + dc);
    recon[1] = CLAMP_255(predict[1] + dc);
    recon[2] = CLAMP_255(predict[2] + dc);
    recon[3] = CLAMP_255(predict[3] + dc);
    predict = recon;
    recon += stride;

    for (i = 1; i < 4; i++)
    {
        recon[0] = predict[0];
        recon[1] = predict[1];
        recon[2] = predict[2];
        recon[3] = predict[3];
        recon += stride;
    }
}


static void
fill_dc_col(unsigned char        *predict,
            int                   stride,
            int                   dc,
            int                   pixel_offset)
{
    int i;

    for (i = 0; i < 4; i++)
    {
        predict[0] = predict[1] = predict[2] = predict[3] =
                CLAMP_255(predict[pixel_offset] + dc);
        predict += stride;
    }
}


static void
fill_dc_block(unsigned char        *predict,
              int                   stride,
              int                   dc)
{
    int i;

    for (i = 0; i < 4; i++)
    {
        predict[0] = CLAMP_255(predict[0] + dc);
        predict[1] = CLAMP_255(predict[1] + dc);
        predict[2] = CLAMP_255(predict[2] + dc);
        predict[3] = CLAMP_255(predict[3] + dc);
        predict += stride;
    }
}


/* TODO: can this map be eliminated by reorering the B_ modes? */
static const enum prediction_mode map_prediction_mode[MB_MODE_COUNT] =
{
    B_DC_PRED, B_VE_PRED, B_HE_PRED, B_TM_PRED, B_TM_PRED,
    B_TM_PRED, B_TM_PRED, B_TM_PRED, B_TM_PRED, B_TM_PRED

};


void
vp8_dixie_idct_add(unsigned char        *predict,
                   int                   stride,
                   const short          *coeffs,
                   struct mb_info       *mbi,
                   int                   block)
{
    int                  do_whole_block;
    int                  row_offset, col_offset;
    enum prediction_mode mode;

    do_whole_block = mbi->base.eob_mask & (1 << block);

    if (mbi->base.y_mode == B_PRED && block < 16)
    {
        mode = mbi->split.modes[block];
        row_offset = 3;
        col_offset = 3;
    }
    else if (block < 16)
    {
        mode = map_prediction_mode[mbi->base.y_mode];
        row_offset = 15 - (block & 0x0c);
        col_offset = 15 - 4 * (block & 0x03);
    }
    else
    {
        mode = map_prediction_mode[mbi->base.uv_mode];
        row_offset = 7 - 2 * (block & 0x02);
        col_offset = 7 - 4 * (block & 0x01);
    }

    if (!do_whole_block)
    {
        int dc;

        dc = ((coeffs[0] + 4) >> 3);

        switch (mode)
        {
            /* Modes that predict to a single value */
        case B_DC_PRED:
            fill_dc_dc(predict, stride, dc, col_offset);
            break;

            /* Modes that predict to a single row */
        case B_VE_PRED:
            fill_dc_row(predict, stride, dc, row_offset);
            break;

            /* Modes that predict to a single column */
        case B_HE_PRED:
            fill_dc_col(predict, stride, dc, col_offset);
            break;

            /* Modes that predict to a whole block */
        default:

            if (dc)
                fill_dc_block(predict, stride, dc);
        }
    }
    else
    {
        switch (mode)
        {
            /* Modes that predict to a single value */
        case B_DC_PRED:
            idct_dc(predict, stride, coeffs, col_offset);
            break;

            /* Modes that predict to a single row */
        case B_VE_PRED:
            idct_row(predict, stride, coeffs, row_offset);
            break;

            /* Modes that predict to a single column */
        case B_HE_PRED:
            idct_col(predict, stride, coeffs, col_offset);
            break;

            /* Modes that predict to a whole block */
        default:
            idct_block(predict, stride, coeffs);
        }
    }
}

void
vp8_dixie_idct_add_process_row(struct vp8_decoder_ctx *ctx,
                               short                  *coeffs,
                               unsigned int            row,
                               unsigned int            start_col,
                               unsigned int            num_cols)
{
    unsigned char  *mb_y, *mb_u, *mb_v;
    int             stride, uv_stride;
    struct mb_info *mbi;
    unsigned int    col;
    int             i;

    /* Adjust pointers based on row, start_col */
    stride    = ctx->ref_frames[CURRENT_FRAME]->img.stride[PLANE_Y];
    uv_stride = ctx->ref_frames[CURRENT_FRAME]->img.stride[PLANE_U];
    mb_y = ctx->ref_frames[CURRENT_FRAME]->img.planes[PLANE_Y];
    mb_u = ctx->ref_frames[CURRENT_FRAME]->img.planes[PLANE_U];
    mb_v = ctx->ref_frames[CURRENT_FRAME]->img.planes[PLANE_V];
    mb_y += (stride * row + start_col) * 16;
    mb_u += (uv_stride * row + start_col) * 8;
    mb_v += (uv_stride * row + start_col) * 8;
    mbi = ctx->mb_info_rows[row] + start_col;
    coeffs += 25 * 16 * start_col;

    for (col = start_col; col < start_col + num_cols; col++)
    {
        short *b_coeffs = coeffs;
        unsigned char *y = mb_y, *u = mb_u, *v = mb_v;

        if (mbi->base.y_mode == B_PRED || mbi->base.y_mode == SPLITMV)
        {
            for (i = 0; i < 16; i++)
            {
                vp8_dixie_idct_add(y, stride, b_coeffs, mbi, i);
                b_coeffs += 16;
                y += 4;

                if ((i & 3) == 3)
                    y += 4 * stride - 16;
            }
        }
        else
        {
            short y2[16];

            if (mbi->base.eob_mask & (1 << 24))
                vp8_dixie_walsh(coeffs + 24 * 16, y2);
            else
                walsh1(coeffs + 24 * 16, y2);

            for (i = 0; i < 16; i++)
            {
                b_coeffs[0] = y2[i];
                vp8_dixie_idct_add(y, stride, b_coeffs, mbi, i);
                b_coeffs += 16;
                y += 4;

                if ((i & 3) == 3)
                    y += 4 * stride - 16;
            }
        }

        for (; i < 20; i++)
        {
            vp8_dixie_idct_add(u, uv_stride, b_coeffs, mbi, i);
            b_coeffs += 16;
            u += 4;

            if (i & 1)
                u += 4 * uv_stride - 8;
        }

        for (; i < 24; i++)
        {
            vp8_dixie_idct_add(v, uv_stride, b_coeffs, mbi, i);
            b_coeffs += 16;
            v += 4;

            if (i & 1)
                v += 4 * uv_stride - 8;
        }

        coeffs += 25 * 16;
        mb_y += 16;
        mb_u += 8;
        mb_v += 8;
    }
}
