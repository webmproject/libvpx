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
#include "predict.h"
#include "idct_add.h"
#include <assert.h>
#include <string.h>

enum
{
    BORDER_PIXELS     = 16,
};


static void
predict_h_nxn(unsigned char *predict,
              int            stride,
              int            n)
{
    unsigned char *left = predict - 1;
    int            i;

    for (i = 0; i < n; i++)
    {
        predict[n-1] = *left;
        predict += stride;
        left += stride;
    }
}


static void
predict_v_nxn(unsigned char *predict,
              int            stride,
              int            n)
{
    unsigned char *above = predict - stride;
    int            i;

    predict += (n - 1) * stride;

    for (i = 0; i < n; i++)
        predict[i] = above[i];
}


static void
predict_tm_nxn(unsigned char *predict,
               int            stride,
               int            n)
{
    /* Transposes the left column to the top row for later consumption
     * by the idct/recon stage
     */
    unsigned char *left = predict - 1;
    unsigned char *above = predict - stride;
    unsigned char  p = above[-1];
    int            i, j;

    for (j = 0; j < n; j++)
    {
        for (i = 0; i < n; i++)
            predict[i] = CLAMP_255(*left + above[i] - p);

        predict += stride;
        left += stride;
    }
}


static void
predict_dc_16x16(unsigned char *predict,
                 int            stride)
{
    unsigned char *left = predict - 1;
    unsigned char *above = predict - stride;
    int            i, dc = 0;

    for (i = 0; i < 16; i++)
    {
        dc += *left + above[i];
        left += stride;
    }

    dc = (dc + 16) >> 5;
    predict[15] = dc;
    predict[15+4*stride] = dc;
    predict[15+8*stride] = dc;
    predict[15+12*stride] = dc;
}


static void
predict_dc_8x8(unsigned char *predict,
               int            stride)
{
    unsigned char *left = predict - 1;
    unsigned char *above = predict - stride;
    int            i, dc = 0;

    for (i = 0; i < 8; i++)
    {
        dc += *left + above[i];
        left += stride;
    }

    dc = (dc + 8) >> 4;
    predict[7] = dc;
    predict[7+4*stride] = dc;
}


static void
predict_dc_4x4(unsigned char *predict,
               int            stride)
{
    unsigned char *left = predict - 1;
    unsigned char *above = predict - stride;
    int            i, dc = 0;

    for (i = 0; i < 4; i++)
    {
        dc += *left + above[i];
        left += stride;
    }

    dc = (dc + 4) >> 3;
    predict[3] = dc;
}


static void
predict_ve_4x4(unsigned char *predict,
               int            stride)
{
    unsigned char *above = predict - stride;

    predict += 3 * stride;
    predict[0] = (above[-1] + 2 * above[0] + above[1] + 2) >> 2;
    predict[1] = (above[ 0] + 2 * above[1] + above[2] + 2) >> 2;
    predict[2] = (above[ 1] + 2 * above[2] + above[3] + 2) >> 2;
    predict[3] = (above[ 2] + 2 * above[3] + above[4] + 2) >> 2;
}


static void
predict_he_4x4(unsigned char *predict,
               int            stride)
{
    predict = predict - 1;
    predict[4] = (predict[-stride] + 2 * predict[0] + predict[stride] + 2) >> 2;
    predict += stride;
    predict[4] = (predict[-stride] + 2 * predict[0] + predict[stride] + 2) >> 2;
    predict += stride;
    predict[4] = (predict[-stride] + 2 * predict[0] + predict[stride] + 2) >> 2;
    predict += stride;
    predict[4] = (predict[-stride] + 2 * predict[0] + predict[0] + 2) >> 2;
}


static void
predict_ld_4x4(unsigned char *predict,
               int            stride)
{
    unsigned char *above = predict - stride;
    int            pred0, pred1, pred2, pred3, pred4, pred5, pred6;

    predict[0] = pred0 = (above[0] + 2 * above[1] + above[2] + 2) >> 2;
    predict[1] = pred1 = (above[1] + 2 * above[2] + above[3] + 2) >> 2;
    predict[2] = pred2 = (above[2] + 2 * above[3] + above[4] + 2) >> 2;
    predict[3] = pred3 = (above[3] + 2 * above[4] + above[5] + 2) >> 2;
    predict += stride;

    predict[0] = pred1;
    predict[1] = pred2;
    predict[2] = pred3;
    predict[3] = pred4 = (above[4] + 2 * above[5] + above[6] + 2) >> 2;
    predict += stride;

    predict[0] = pred2;
    predict[1] = pred3;
    predict[2] = pred4;
    predict[3] = pred5 = (above[5] + 2 * above[6] + above[7] + 2) >> 2;
    predict += stride;

    predict[0] = pred3;
    predict[1] = pred4;
    predict[2] = pred5;
    predict[3] = pred6 = (above[6] + 2 * above[7] + above[7] + 2) >> 2;
}


static void
predict_rd_4x4(unsigned char *predict,
               int            stride)
{
    unsigned char *left = predict - 1;
    unsigned char *above = predict - stride;
    int            pred0, pred1, pred2, pred3, pred4, pred5, pred6;

    predict[0] = pred0 = (left[ 0] + 2 * above[-1] + above[0] + 2) >> 2;
    predict[1] = pred1 = (above[-1] + 2 * above[ 0] + above[1] + 2) >> 2;
    predict[2] = pred2 = (above[ 0] + 2 * above[ 1] + above[2] + 2) >> 2;
    predict[3] = pred3 = (above[ 1] + 2 * above[ 2] + above[3] + 2) >> 2;
    predict += stride;

    predict[0] = pred4 = (left[stride] + 2 * left[0] + above[-1] + 2) >> 2;
    predict[1] = pred0;
    predict[2] = pred1;
    predict[3] = pred2;
    predict += stride;

    predict[0] = pred5 = (left[stride*2] + 2 * left[stride] + left[0] + 2) >> 2;
    predict[1] = pred4;
    predict[2] = pred0;
    predict[3] = pred1;
    predict += stride;

    predict[0] = pred6 = (left[stride*3] + 2 * left[stride*2] + left[stride] + 2) >> 2;
    predict[1] = pred5;
    predict[2] = pred4;
    predict[3] = pred0;
}


static void
predict_vr_4x4(unsigned char *predict,
               int            stride)
{
    unsigned char *left = predict - 1;
    unsigned char *above = predict - stride;
    int            pred0, pred1, pred2, pred3, pred4, pred5, pred6, pred7,
      pred8, pred9;

    predict[0] = pred0 = (above[-1] + above[0] + 1) >> 1;
    predict[1] = pred1 = (above[ 0] + above[1] + 1) >> 1;
    predict[2] = pred2 = (above[ 1] + above[2] + 1) >> 1;
    predict[3] = pred3 = (above[ 2] + above[3] + 1) >> 1;
    predict += stride;

    predict[0] = pred4 = (left[ 0] + 2 * above[-1] + above[0] + 2) >> 2;
    predict[1] = pred5 = (above[-1] + 2 * above[ 0] + above[1] + 2) >> 2;
    predict[2] = pred6 = (above[ 0] + 2 * above[ 1] + above[2] + 2) >> 2;
    predict[3] = pred7 = (above[ 1] + 2 * above[ 2] + above[3] + 2) >> 2;
    predict += stride;

    predict[0] = pred8 = (left[stride] + 2 * left[0] + above[-1] + 2) >> 2;
    predict[1] = pred0;
    predict[2] = pred1;
    predict[3] = pred2;
    predict += stride;

    predict[0] = pred9 = (left[stride*2] + 2 * left[stride] + left[0] + 2) >> 2;
    predict[1] = pred4;
    predict[2] = pred5;
    predict[3] = pred6;
}


static void
predict_vl_4x4(unsigned char *predict,
               int            stride)
{
    unsigned char *above = predict - stride;
    int            pred0, pred1, pred2, pred3, pred4, pred5, pred6, pred7,
      pred8, pred9;

    predict[0] = pred0 = (above[0] + above[1] + 1) >> 1;
    predict[1] = pred1 = (above[1] + above[2] + 1) >> 1;
    predict[2] = pred2 = (above[2] + above[3] + 1) >> 1;
    predict[3] = pred3 = (above[3] + above[4] + 1) >> 1;
    predict += stride;

    predict[0] = pred4 = (above[0] + 2 * above[1] + above[2] + 2) >> 2;
    predict[1] = pred5 = (above[1] + 2 * above[2] + above[3] + 2) >> 2;
    predict[2] = pred6 = (above[2] + 2 * above[3] + above[4] + 2) >> 2;
    predict[3] = pred7 = (above[3] + 2 * above[4] + above[5] + 2) >> 2;
    predict += stride;

    predict[0] = pred1;
    predict[1] = pred2;
    predict[2] = pred3;
    predict[3] = pred8 = (above[4] + 2 * above[5] + above[6] + 2) >> 2;
    predict += stride;

    predict[0] = pred5;
    predict[1] = pred6;
    predict[2] = pred7;
    predict[3] = pred9 = (above[5] + 2 * above[6] + above[7] + 2) >> 2;
}


static void
predict_hd_4x4(unsigned char *predict,
               int            stride)
{
    unsigned char *left = predict - 1;
    unsigned char *above = predict - stride;
    int            pred0, pred1, pred2, pred3, pred4, pred5, pred6, pred7,
      pred8, pred9;

    predict[0] = pred0 = (left[ 0] + above[-1] + 1) >> 1;
    predict[1] = pred1 = (left[ 0] + 2 * above[-1] + above[0] + 2) >> 2;
    predict[2] = pred2 = (above[-1] + 2 * above[ 0] + above[1] + 2) >> 2;
    predict[3] = pred3 = (above[ 0] + 2 * above[ 1] + above[2] + 2) >> 2;
    predict += stride;

    predict[0] = pred4 = (left[stride] + left[0] + 1) >> 1;
    predict[1] = pred5 = (left[stride] + 2 * left[0] + above[-1] + 2) >> 2;
    predict[2] = pred0;
    predict[3] = pred1;
    predict += stride;

    predict[0] = pred6 = (left[stride*2] +   left[stride] + 1) >> 1;
    predict[1] = pred7 = (left[stride*2] + 2 * left[stride] + left[0] + 2) >> 2;
    predict[2] = pred4;
    predict[3] = pred5;
    predict += stride;

    predict[0] = pred8 = (left[stride*3] +   left[stride*2] + 1) >> 1;
    predict[1] = pred9 = (left[stride*3] + 2 * left[stride*2] + left[stride] + 2) >> 2;
    predict[2] = pred6;
    predict[3] = pred7;
}


static void
predict_hu_4x4(unsigned char *predict,
               int            stride)
{
    unsigned char *left = predict - 1;
    int            pred0, pred1, pred2, pred3, pred4, pred5, pred6;

    predict[0] = pred0 = (left[stride*0] +   left[stride*1] + 1) >> 1;
    predict[1] = pred1 = (left[stride*0] + 2 * left[stride*1] + left[stride*2] + 2) >> 2;
    predict[2] = pred2 = (left[stride*1] +   left[stride*2] + 1) >> 1;
    predict[3] = pred3 = (left[stride*1] + 2 * left[stride*2] + left[stride*3] + 2) >> 2;
    predict += stride;

    predict[0] = pred2;
    predict[1] = pred3;
    predict[2] = pred4 = (left[stride*2] + left[stride*3] + 1) >> 1;
    predict[3] = pred5 = (left[stride*2] + 2 * left[stride*3] + left[stride*3] + 2) >> 2;
    predict += stride;

    predict[0] = pred4;
    predict[1] = pred5;
    predict[2] = pred6 = left[stride*3];
    predict[3] = pred6;
    predict += stride;

    predict[0] = pred6;
    predict[1] = pred6;
    predict[2] = pred6;
    predict[3] = pred6;
}


static void
predict_h_16x16(unsigned char *predict, int stride)
{
    predict_h_nxn(predict, stride, 16);
}


static void
predict_v_16x16(unsigned char *predict, int stride)
{
    predict_v_nxn(predict, stride, 16);
}


static void
predict_tm_16x16(unsigned char *predict, int stride)
{
    predict_tm_nxn(predict, stride, 16);
}


static void
predict_h_8x8(unsigned char *predict, int stride)
{
    predict_h_nxn(predict, stride, 8);
}


static void
predict_v_8x8(unsigned char *predict, int stride)
{
    predict_v_nxn(predict, stride, 8);
}


static void
predict_tm_8x8(unsigned char *predict, int stride)
{
    predict_tm_nxn(predict, stride, 8);
}


static void
predict_h_4x4(unsigned char *predict, int stride)
{
    predict_h_nxn(predict, stride, 4);
}


static void
predict_v_4x4(unsigned char *predict, int stride)
{
    predict_v_nxn(predict, stride, 4);
}


static void
predict_tm_4x4(unsigned char *predict, int stride)
{
    predict_tm_nxn(predict, stride, 4);
}


static void
copy_down(unsigned char           *recon,
          int                      stride)
{
    /* Copy the four pixels above-right of subblock 3 to
     * above-right of subblocks 7, 11, and 15
     */
    uint32_t tmp, *copy = (void *)(recon + 16 - stride);

    stride = stride / sizeof(unsigned int);
    tmp = *copy;
    copy += stride * 4;
    *copy = tmp;
    copy += stride * 4;
    *copy = tmp;
    copy += stride * 4;
    *copy = tmp;
}


static void
b_pred(unsigned char  *predict,
       int             stride,
       struct mb_info *mbi,
       short          *coeffs)
{
    int i;

    copy_down(predict, stride);

    for (i = 0; i < 16; i++)
    {
        unsigned char *b_predict = predict + (i & 3) * 4;

        switch (mbi->split.modes[i])
        {
        case B_DC_PRED:
            predict_dc_4x4(b_predict, stride);
            break;
        case B_TM_PRED:
            predict_tm_4x4(b_predict, stride);
            break;
        case B_VE_PRED:
            predict_ve_4x4(b_predict, stride);
            break;
        case B_HE_PRED:
            predict_he_4x4(b_predict, stride);
            break;
        case B_LD_PRED:
            predict_ld_4x4(b_predict, stride);
            break;
        case B_RD_PRED:
            predict_rd_4x4(b_predict, stride);
            break;
        case B_VR_PRED:
            predict_vr_4x4(b_predict, stride);
            break;
        case B_VL_PRED:
            predict_vl_4x4(b_predict, stride);
            break;
        case B_HD_PRED:
            predict_hd_4x4(b_predict, stride);
            break;
        case B_HU_PRED:
            predict_hu_4x4(b_predict, stride);
            break;
        default:
            assert(0);
        }

        vp8_dixie_idct_add(b_predict, stride, coeffs, mbi, i);
        coeffs += 16;

        if ((i & 3) == 3)
        {
            predict += stride * 4;
        }
    }
}

static void
predict_intra_luma(unsigned char   *predict,
                   int              stride,
                   struct mb_info  *mbi,
                   const ptrdiff_t  reference_offsets[4],
                   short           *coeffs)
{
    if (mbi->base.y_mode == B_PRED)
        b_pred(predict, stride, mbi, coeffs);
    else
    {
        short y2[16];
        int i;

        switch (mbi->base.y_mode)
        {
        case DC_PRED:
            predict_dc_16x16(predict, stride);
            break;
        case V_PRED:
            predict_v_16x16(predict, stride);
            break;
        case H_PRED:
            predict_h_16x16(predict, stride);
            break;
        case TM_PRED:
            predict_tm_16x16(predict, stride);
            break;
        default:
            assert(0);
        }

        /* TODO: clean up, this will be dup'd for inter */
        if (mbi->base.eob_mask & (1 << 24))
        {
            vp8_dixie_walsh(coeffs + 24 * 16, y2);

            for (i = 0; i < 16; i++)
                coeffs[i*16] = y2[i];
        }
        else
        {
            int dc = ((coeffs[24*16] + 3) >> 3);

            for (i = 0; i < 16; i++)
                coeffs[i*16] = dc;
        }

        for (i = 0; i < 16; i++)
        {
            vp8_dixie_idct_add(predict, stride, coeffs, mbi, i);
            coeffs += 16;
            predict += 4;

            if ((i & 3) == 3)
                predict += stride * 4 - 16;
        }

    }
}


static void
predict_intra_chroma(unsigned char   *predict_u,
                     unsigned char   *predict_v,
                     int              stride,
                     struct mb_info  *mbi,
                     const ptrdiff_t  reference_offsets[4],
                     short           *coeffs)
{
    int i;

    switch (mbi->base.uv_mode)
    {
    case DC_PRED:
        predict_dc_8x8(predict_u, stride);
        predict_dc_8x8(predict_v, stride);
        break;
    case V_PRED:
        predict_v_8x8(predict_u, stride);
        predict_v_8x8(predict_v, stride);
        break;
    case H_PRED:
        predict_h_8x8(predict_u, stride);
        predict_h_8x8(predict_v, stride);
        break;
    case TM_PRED:
        predict_tm_8x8(predict_u, stride);
        predict_tm_8x8(predict_v, stride);
        break;
    default:
        assert(0);
    }

    coeffs += 16 * 16;

    for (i = 16; i < 20; i++)
    {
        vp8_dixie_idct_add(predict_u, stride, coeffs, mbi, i);
        coeffs += 16;
        predict_u += 4;

        if (i & 1)
            predict_u += stride * 4 - 8;
    }

    for (i = 20; i < 24; i++)
    {
        vp8_dixie_idct_add(predict_v, stride, coeffs, mbi, i);
        coeffs += 16;
        predict_v += 4;

        if (i & 1)
            predict_v += stride * 4 - 8;
    }

}

static void
release_ref_frame(struct ref_cnt_img *rcimg)
{
    if (rcimg)
    {
        assert(rcimg->ref_cnt);
        rcimg->ref_cnt--;
    }
}


static struct ref_cnt_img *
ref_frame(struct ref_cnt_img *rcimg)
{
    rcimg->ref_cnt++;
    return rcimg;
}


static struct ref_cnt_img *
find_free_ref_frame(struct ref_cnt_img *frames)
{
    int i;

    for (i = 0; i < NUM_REF_FRAMES; i++)
        if (frames[i].ref_cnt == 0)
        {
            frames[i].ref_cnt = 1;
            return &frames[i];
        }

    assert(0);
    return NULL;
}


static void
fixup_left(unsigned char        *predict,
           int                   width,
           int                   stride,
           unsigned int          row,
           enum prediction_mode  mode)
{
    /* The left column of out-of-frame pixels is taken to be 129,
     * unless we're doing DC_PRED, in which case we duplicate the
     * above row, unless this is also row 0, in which case we use
     * 129.
     */
    unsigned char *left = predict - 1;
    int i;

    if (mode == DC_PRED && row)
    {
        unsigned char *above = predict - stride;

        for (i = 0; i < width; i++)
        {
            *left = above[i];
            left += stride;
        }
    }
    else
    {
        /* Need to re-set the above row, in case the above MB was
         * DC_PRED.
         */
        left -= stride;

        for (i = -1; i < width; i++)
        {
            *left = 129;
            left += stride;
        }
    }
}


static void
fixup_above(unsigned char        *predict,
            int                   width,
            int                   stride,
            unsigned int          col,
            enum prediction_mode  mode)
{
    /* The above row of out-of-frame pixels is taken to be 127,
     * unless we're doing DC_PRED, in which case we duplicate the
     * left col, unless this is also col 0, in which case we use
     * 127.
     */
    unsigned char *above = predict - stride;
    int i;

    if (mode == DC_PRED && col)
    {
        unsigned char *left = predict - 1;

        for (i = 0; i < width; i++)
        {
            above[i] = *left;
            left += stride;
        }
    }
    else
        /* Need to re-set the left col, in case the last MB was
         * DC_PRED.
         */
        memset(above - 1, 127, width + 1);

    memset(above + width, 127, 4); // for above-right subblock modes
}


void
vp8_dixie_predict_init(struct vp8_decoder_ctx *ctx)
{

    int i;
    unsigned char *this_frame_base;

    if (ctx->frame_hdr.frame_size_updated)
    {
        for (i = 0; i < NUM_REF_FRAMES; i++)
        {
            unsigned int w = ctx->mb_cols * 16 + BORDER_PIXELS * 2;
            unsigned int h = ctx->mb_rows * 16 + BORDER_PIXELS * 2;

            vpx_img_free(&ctx->frame_strg[i].img);
            ctx->frame_strg[i].ref_cnt = 0;

            if (!vpx_img_alloc(&ctx->frame_strg[i].img,
                               IMG_FMT_I420, w, h, 16))
                vpx_internal_error(&ctx->error, VPX_CODEC_MEM_ERROR,
                                   "Failed to allocate %dx%d framebuffer",
                                   w, h);

            vpx_img_set_rect(&ctx->frame_strg[i].img,
                             BORDER_PIXELS, BORDER_PIXELS,
                             ctx->frame_hdr.kf.w, ctx->frame_hdr.kf.h);

        }
    }

    /* Find a free framebuffer to predict into */
    if (ctx->ref_frames[CURRENT_FRAME])
        release_ref_frame(ctx->ref_frames[CURRENT_FRAME]);

    ctx->ref_frames[CURRENT_FRAME] = find_free_ref_frame(ctx->frame_strg);
    this_frame_base = ctx->ref_frames[CURRENT_FRAME]->img.img_data;

    /* Calculate offsets to the other reference frames */
    for (i = 0; i < NUM_REF_FRAMES; i++)
    {
        struct ref_cnt_img  *ref = ctx->ref_frames[i];

        ctx->ref_frame_offsets[i] = ref ? ref->img.img_data - this_frame_base
                                    : 0;
    }
}


void
vp8_dixie_predict_process_row(struct vp8_decoder_ctx *ctx,
                              unsigned int            row,
                              unsigned int            start_col,
                              unsigned int            num_cols)
{
    unsigned char  *y, *u, *v;
    int             stride, uv_stride;
    struct mb_info *mbi;
    unsigned int    col;
    short          *coeffs;

    /* Adjust pointers based on row, start_col */
    stride    = ctx->ref_frames[CURRENT_FRAME]->img.stride[PLANE_Y];
    uv_stride = ctx->ref_frames[CURRENT_FRAME]->img.stride[PLANE_U];
    y = ctx->ref_frames[CURRENT_FRAME]->img.planes[PLANE_Y];
    u = ctx->ref_frames[CURRENT_FRAME]->img.planes[PLANE_U];
    v = ctx->ref_frames[CURRENT_FRAME]->img.planes[PLANE_V];
    y += (stride * row + start_col) * 16;
    u += (uv_stride * row + start_col) * 8;
    v += (uv_stride * row + start_col) * 8;
    mbi = ctx->mb_info_rows[row] + start_col;
    coeffs = ctx->tokens[row & (ctx->token_hdr.partitions - 1)].coeffs;

    /* Fix up the out-of-frame pixels */
    if (start_col == 0)
    {
        fixup_left(y, 16, stride, row, mbi->base.y_mode);
        fixup_left(u, 8, uv_stride, row, mbi->base.uv_mode);
        fixup_left(v, 8, uv_stride, row, mbi->base.uv_mode);

        if (row == 0)
            *(y - stride - 1) = 127;
    }

    for (col = start_col; col < start_col + num_cols; col++)
    {
        if (row == 0)
        {
            fixup_above(y, 16, stride, col, mbi->base.y_mode);
            fixup_above(u, 8, uv_stride, col, mbi->base.uv_mode);
            fixup_above(v, 8, uv_stride, col, mbi->base.uv_mode);
        }

        if (mbi->base.y_mode <= B_PRED)
        {
            predict_intra_luma(y, stride, mbi, ctx->ref_frame_offsets, coeffs);
            predict_intra_chroma(u, v, uv_stride, mbi, ctx->ref_frame_offsets,
                                 coeffs);
        }
        else
            assert(0);

        /* Advance to the next macroblock */
        mbi++;
        y += 16;
        u += 8;
        v += 8;
        coeffs += 25 * 16;
    }

    if (col == ctx->mb_cols)
    {
        /* Extend the last row by four pixels for intra prediction. This will
         * be propagated later by copy_down.
         */
        uint32_t *extend = (uint32_t *)(y + 15 * stride);
        uint32_t  val = 0x01010101 * y[-1 + 15 * stride];
        *extend = val;
    }
}
