/*
 *  Copyright (c) 2011 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include "error_concealment.h"
#include "onyxd_int.h"
#include "vpx_mem/vpx_mem.h"

void vp8_estimate_missing_mvs(VP8D_COMP *pbi)
{
    VP8_COMMON *const pc = &pbi->common;
    unsigned int first_corrupt = pbi->mvs_corrupt_from_mb;
    const unsigned int num_mbs = pc->mb_rows * pc->mb_cols;
    if (first_corrupt < num_mbs)
    {
        MODE_INFO *mi = pc->mi;
        MODE_INFO *correct_mi;
        const int num_corrupt = num_mbs - first_corrupt;
        int i;
        int mb_row, mb_col;
        if (first_corrupt == 0)
        {
            /* if the first MB is corrupt we just copy from it
               the previous frame */
            mi[0].mbmi.mv.as_int = 0;
            mi[0].mbmi.mode = ZEROMV;
            mi[0].mbmi.uv_mode = ZEROMV;
            mi[0].mbmi.ref_frame = LAST_FRAME;
            first_corrupt = 1;
            correct_mi = mi + 1;
        }
        for (mb_row = 0; mb_row < pc->mb_rows; ++mb_row)
        {
            for (mb_col = 0; mb_col < pc->mb_cols; ++mb_col)
            {
                int mb_idx = mb_row * pc->mb_cols + mb_col;
                if (mb_idx >= first_corrupt)
                {
                    *mi = *correct_mi;
                }
                else if (mb_idx == first_corrupt - 1)
                {
                    correct_mi = mi;
                }
                ++mi;
            }
            ++mi;
        }
    }
}

void vp8_conceal_corrupt_block(MACROBLOCKD *xd)
{
    /* this macroblock has corrupt residual, use the motion compensated
       image for concealment */
    int i;
    for (i=0; i < 16; i++)
        vpx_memcpy(xd->dst.y_buffer + i*xd->dst.y_stride,
                   xd->predictor + i*16, 16);
    for (i=0; i < 8; i++)
        vpx_memcpy(xd->dst.u_buffer + i*xd->dst.uv_stride,
                   xd->predictor + 256 + i*8, 8);
    for (i=0; i < 8; i++)
        vpx_memcpy(xd->dst.v_buffer + i*xd->dst.uv_stride,
                   xd->predictor + 320 + i*8, 8);
}
