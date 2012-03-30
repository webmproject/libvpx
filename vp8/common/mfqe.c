/*
 *  Copyright (c) 2012 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */


/* MFQE: Multiframe Quality Enhancement
 * In rate limited situations keyframes may cause significant visual artifacts
 * commonly referred to as "popping." This file implements a postproccesing
 * algorithm which blends data from the preceeding frame when there is no
 * motion and the q from the previous frame is lower which indicates that it is
 * higher quality.
 */

#include "postproc.h"
#include "variance.h"
#include "vpx_mem/vpx_mem.h"
#include "vpx_rtcd.h"
#include "vpx_scale/yv12config.h"

#include <limits.h>
#include <stdlib.h>


static void filter_by_weight(unsigned char *src, int src_stride,
                             unsigned char *dst, int dst_stride,
                             int block_size, int src_weight)
{
    int dst_weight = (1 << MFQE_PRECISION) - src_weight;
    int rounding_bit = 1 << (MFQE_PRECISION - 1);
    int r, c;

    for (r = 0; r < block_size; r++)
    {
        for (c = 0; c < block_size; c++)
        {
            dst[c] = (src[c] * src_weight +
                      dst[c] * dst_weight +
                      rounding_bit) >> MFQE_PRECISION;
        }
        src += src_stride;
        dst += dst_stride;
    }
}

void vp8_filter_by_weight16x16_c(unsigned char *src, int src_stride,
                                 unsigned char *dst, int dst_stride,
                                 int src_weight)
{
    filter_by_weight(src, src_stride, dst, dst_stride, 16, src_weight);
}

void vp8_filter_by_weight8x8_c(unsigned char *src, int src_stride,
                               unsigned char *dst, int dst_stride,
                               int src_weight)
{
    filter_by_weight(src, src_stride, dst, dst_stride, 8, src_weight);
}

void vp8_filter_by_weight4x4_c(unsigned char *src, int src_stride,
                               unsigned char *dst, int dst_stride,
                               int src_weight)
{
    filter_by_weight(src, src_stride, dst, dst_stride, 4, src_weight);
}

static void apply_ifactor(unsigned char *y_src,
                          int y_src_stride,
                          unsigned char *y_dst,
                          int y_dst_stride,
                          unsigned char *u_src,
                          unsigned char *v_src,
                          int uv_src_stride,
                          unsigned char *u_dst,
                          unsigned char *v_dst,
                          int uv_dst_stride,
                          int block_size,
                          int src_weight)
{
    if (block_size == 16)
    {
        vp8_filter_by_weight16x16(y_src, y_src_stride, y_dst, y_dst_stride, src_weight);
        vp8_filter_by_weight8x8(u_src, uv_src_stride, u_dst, uv_dst_stride, src_weight);
        vp8_filter_by_weight8x8(v_src, uv_src_stride, v_dst, uv_dst_stride, src_weight);
    }
    else /* if (block_size == 8) */
    {
        vp8_filter_by_weight8x8(y_src, y_src_stride, y_dst, y_dst_stride, src_weight);
        vp8_filter_by_weight4x4(u_src, uv_src_stride, u_dst, uv_dst_stride, src_weight);
        vp8_filter_by_weight4x4(v_src, uv_src_stride, v_dst, uv_dst_stride, src_weight);
    }
}

static void multiframe_quality_enhance_block
(
    int blksize, /* Currently only values supported are 16, 8, 4 */
    int qcurr,
    int qprev,
    unsigned char *y,
    unsigned char *u,
    unsigned char *v,
    int y_stride,
    int uv_stride,
    unsigned char *yd,
    unsigned char *ud,
    unsigned char *vd,
    int yd_stride,
    int uvd_stride
)
{
    static const unsigned char VP8_ZEROS[16]=
    {
         0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
    };
    int blksizeby2 = blksize >> 1;
    int qdiff = qcurr - qprev;

    int i, j;
    unsigned char *yp;
    unsigned char *ydp;
    unsigned char *up;
    unsigned char *udp;
    unsigned char *vp;
    unsigned char *vdp;

    unsigned int act_sum = 0, sse, sad_sum = 0, thr, uvsad = UINT_MAX;
    unsigned int act[4], sad[4];

    if (blksize == 16)
    {
        act[0] = (vp8_variance8x8(yd, yd_stride,
                                  VP8_ZEROS, 0, &sse) + 32) >> 6;
        act_sum = act[0];
        sad[0] = (vp8_sad8x8(y, y_stride,
                             yd, yd_stride, INT_MAX) + 32) >> 6;
        sad_sum = sad[0];

        act[1] = (vp8_variance8x8(yd + 8, yd_stride,
                                  VP8_ZEROS, 0, &sse) + 32) >> 6;
        act_sum += act[1];
        sad[1] = (vp8_sad8x8(y + 8, y_stride,
                             yd + 8, yd_stride, INT_MAX) + 32) >> 6;
        sad_sum = sad[1];

        act[2] = (vp8_variance8x8(yd + 8 * yd_stride, yd_stride,
                                  VP8_ZEROS, 0, &sse) + 32) >> 6;
        act_sum += act[2];
        sad[2] = (vp8_sad8x8(y + 8 * y_stride, y_stride,
                             yd + 8 * yd_stride, yd_stride, INT_MAX) + 32) >> 6;
        sad_sum = sad[2];

        act[3] = (vp8_variance8x8(yd + 8 * yd_stride + 8, yd_stride,
                                  VP8_ZEROS, 0, &sse) + 32) >> 6;
        act_sum += act[3];
        sad[3] = (vp8_sad8x8(y + 8 * y_stride + 8, y_stride,
                             yd + 8 * y_stride + 8, yd_stride, INT_MAX)
                  + 32) >> 6;
        sad_sum = sad[3];
    }
    else if (blksize == 8)
    {
        act[0] = (vp8_variance4x4(yd, yd_stride,
                                  VP8_ZEROS, 0, &sse) + 4) >> 4;
        act_sum = act[0];
        sad[0] = (vp8_sad4x4(y, y_stride,
                             yd, yd_stride, INT_MAX) + 4) >> 4;
        sad_sum = sad[0];

        act[1] = (vp8_variance4x4(yd + 4, yd_stride,
                                  VP8_ZEROS, 0, &sse) + 4) >> 4;
        act_sum += act[1];
        sad[1] = (vp8_sad4x4(y + 4, y_stride,
                             yd + 4, yd_stride, INT_MAX) + 4) >> 4;
        sad_sum = sad[1];

        act[2] = (vp8_variance4x4(yd + 4 * yd_stride, yd_stride,
                                  VP8_ZEROS, 0, &sse) + 4) >> 4;
        act_sum += act[2];
        sad[2] = (vp8_sad4x4(y + 4 * y_stride, y_stride,
                             yd + 4 * yd_stride, yd_stride, INT_MAX) + 4) >> 4;
        sad_sum = sad[2];

        act[3] = (vp8_variance4x4(yd + 4 * yd_stride + 4, yd_stride,
                                  VP8_ZEROS, 0, &sse) + 4) >> 4;
        act_sum += act[3];
        sad[3] = (vp8_sad4x4(y + 4 * y_stride + 4, y_stride,
                             yd + 4 * y_stride + 4, yd_stride, INT_MAX)
                  + 4) >> 4;
        sad_sum = sad[3];
    }
    else
    {
        act_sum = (vp8_variance4x4(yd, yd_stride, VP8_ZEROS, 0, &sse) + 8) >> 4;
        sad_sum = (vp8_sad4x4(y, y_stride, yd, yd_stride, INT_MAX) + 8) >> 4;
    }

    /* thr = qdiff/8 + log2(act) + log4(qprev) */
    thr = (qdiff >> 3);

    while (qprev >>= 2) thr++;

    if (blksize > 4)
    {
        unsigned int base_thr = thr, this_thr, this_act;
        int i;

        for (i = 0; i < 4; i++)
        {
            this_thr = base_thr;
            this_act = act[i];

            while (this_act >>= 1) this_thr++;

            if (sad[i] >= this_thr || act[i] < 16)
            {
                sad_sum = UINT_MAX;
                break;
            }
        }
    }

    while (act_sum >>= 1) thr++;

    if (sad_sum < thr)
    {
        if (blksize == 16)
        {
            uvsad = (vp8_sad8x8(u, uv_stride, ud, uvd_stride, INT_MAX) + 32)
                    >> 6;

            if (uvsad < thr)
                uvsad = (vp8_sad8x8(v, uv_stride, vd, uvd_stride, INT_MAX) + 32)
                        >> 6;
        }
        else
        {
            uvsad = (vp8_sad4x4(u, uv_stride, ud, uvd_stride, INT_MAX) + 8)
                    >> 4;

            if (uvsad < thr)
                uvsad = (vp8_sad4x4(v, uv_stride, vd, uvd_stride, INT_MAX) + 8)
                        >> 4;
        }
    }

    if (uvsad < thr)
    {
        static const int roundoff = (1 << (MFQE_PRECISION - 1));
        int ifactor = (sad_sum << MFQE_PRECISION) / thr;
        ifactor >>= (qdiff >> 5);
        // TODO: SIMD optimize this section
        if (ifactor)
        {
            int icfactor = (1 << MFQE_PRECISION) - ifactor;
            for (yp = y, ydp = yd, i = 0; i < blksize; ++i, yp += y_stride, ydp += yd_stride)
            {
                for (j = 0; j < blksize; ++j)
                    ydp[j] = (int)((yp[j] * ifactor + ydp[j] * icfactor + roundoff) >> MFQE_PRECISION);
            }
            for (up = u, udp = ud, i = 0; i < blksizeby2; ++i, up += uv_stride, udp += uvd_stride)
            {
                for (j = 0; j < blksizeby2; ++j)
                    udp[j] = (int)((up[j] * ifactor + udp[j] * icfactor + roundoff) >> MFQE_PRECISION);
            }
            for (vp = v, vdp = vd, i = 0; i < blksizeby2; ++i, vp += uv_stride, vdp += uvd_stride)
            {
                for (j = 0; j < blksizeby2; ++j)
                    vdp[j] = (int)((vp[j] * ifactor + vdp[j] * icfactor + roundoff) >> MFQE_PRECISION);
            }
        }
    }
    else
    {
        if (blksize == 16)
        {
            vp8_copy_mem16x16(y, y_stride, yd, yd_stride);
            vp8_copy_mem8x8(u, uv_stride, ud, uvd_stride);
            vp8_copy_mem8x8(v, uv_stride, vd, uvd_stride);
        }
        else if (blksize == 8)
        {
            vp8_copy_mem8x8(y, y_stride, yd, yd_stride);
            for (up = u, udp = ud, i = 0; i < blksizeby2; ++i, up += uv_stride, udp += uvd_stride)
                vpx_memcpy(udp, up, blksizeby2);
            for (vp = v, vdp = vd, i = 0; i < blksizeby2; ++i, vp += uv_stride, vdp += uvd_stride)
                vpx_memcpy(vdp, vp, blksizeby2);
        }
        else
        {
            for (yp = y, ydp = yd, i = 0; i < blksize; ++i, yp += y_stride, ydp += yd_stride)
                vpx_memcpy(ydp, yp, blksize);
            for (up = u, udp = ud, i = 0; i < blksizeby2; ++i, up += uv_stride, udp += uvd_stride)
                vpx_memcpy(udp, up, blksizeby2);
            for (vp = v, vdp = vd, i = 0; i < blksizeby2; ++i, vp += uv_stride, vdp += uvd_stride)
                vpx_memcpy(vdp, vp, blksizeby2);
        }
    }
}

void vp8_multiframe_quality_enhance
(
    VP8_COMMON *cm
)
{
    YV12_BUFFER_CONFIG *show = cm->frame_to_show;
    YV12_BUFFER_CONFIG *dest = &cm->post_proc_buffer;

    FRAME_TYPE frame_type = cm->frame_type;
    /* Point at base of Mb MODE_INFO list has motion vectors etc */
    const MODE_INFO *mode_info_context = cm->mi;
    int mb_row;
    int mb_col;
    int qcurr = cm->base_qindex;
    int qprev = cm->postproc_state.last_base_qindex;

    unsigned char *y_ptr, *u_ptr, *v_ptr;
    unsigned char *yd_ptr, *ud_ptr, *vd_ptr;

    /* Set up the buffer pointers */
    y_ptr = show->y_buffer;
    u_ptr = show->u_buffer;
    v_ptr = show->v_buffer;
    yd_ptr = dest->y_buffer;
    ud_ptr = dest->u_buffer;
    vd_ptr = dest->v_buffer;

    /* postprocess each macro block */
    for (mb_row = 0; mb_row < cm->mb_rows; mb_row++)
    {
        for (mb_col = 0; mb_col < cm->mb_cols; mb_col++)
        {
            /* if motion is high there will likely be no benefit */
            if (((frame_type == INTER_FRAME &&
                  abs(mode_info_context->mbmi.mv.as_mv.row) <= 10 &&
                  abs(mode_info_context->mbmi.mv.as_mv.col) <= 10) ||
                 (frame_type == KEY_FRAME)))
            {
                if (mode_info_context->mbmi.mode == B_PRED || mode_info_context->mbmi.mode == SPLITMV)
                {
                    int i, j;
                    for (i=0; i<2; ++i)
                        for (j=0; j<2; ++j)
                            multiframe_quality_enhance_block(8, qcurr, qprev,
                                                             y_ptr + 8*(i*show->y_stride+j),
                                                             u_ptr + 4*(i*show->uv_stride+j),
                                                             v_ptr + 4*(i*show->uv_stride+j),
                                                             show->y_stride,
                                                             show->uv_stride,
                                                             yd_ptr + 8*(i*dest->y_stride+j),
                                                             ud_ptr + 4*(i*dest->uv_stride+j),
                                                             vd_ptr + 4*(i*dest->uv_stride+j),
                                                             dest->y_stride,
                                                             dest->uv_stride);
                }
                else
                {
                    multiframe_quality_enhance_block(16, qcurr, qprev, y_ptr,
                                                     u_ptr, v_ptr,
                                                     show->y_stride,
                                                     show->uv_stride,
                                                     yd_ptr, ud_ptr, vd_ptr,
                                                     dest->y_stride,
                                                     dest->uv_stride);
                }
            }
            else
            {
                vp8_copy_mem16x16(y_ptr, show->y_stride, yd_ptr, dest->y_stride);
                vp8_copy_mem8x8(u_ptr, show->uv_stride, ud_ptr, dest->uv_stride);
                vp8_copy_mem8x8(v_ptr, show->uv_stride, vd_ptr, dest->uv_stride);
            }
            y_ptr += 16;
            u_ptr += 8;
            v_ptr += 8;
            yd_ptr += 16;
            ud_ptr += 8;
            vd_ptr += 8;
            mode_info_context++;     /* step to next MB */
        }

        y_ptr += show->y_stride  * 16 - 16 * cm->mb_cols;
        u_ptr += show->uv_stride *  8 - 8 * cm->mb_cols;
        v_ptr += show->uv_stride *  8 - 8 * cm->mb_cols;
        yd_ptr += dest->y_stride  * 16 - 16 * cm->mb_cols;
        ud_ptr += dest->uv_stride *  8 - 8 * cm->mb_cols;
        vd_ptr += dest->uv_stride *  8 - 8 * cm->mb_cols;

        mode_info_context++;         /* Skip border mb */
    }
}
