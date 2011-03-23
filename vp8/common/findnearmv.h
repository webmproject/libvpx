/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */


#ifndef __INC_FINDNEARMV_H
#define __INC_FINDNEARMV_H

#include "mv.h"
#include "blockd.h"
#include "modecont.h"
#include "treecoder.h"

typedef union
{
    unsigned int as_int;
    MV           as_mv;
} int_mv;        /* facilitates rapid equality tests */

static void mv_bias(int refmb_ref_frame_sign_bias, int refframe, int_mv *mvp, const int *ref_frame_sign_bias)
{
    MV xmv;
    xmv = mvp->as_mv;

    if (refmb_ref_frame_sign_bias != ref_frame_sign_bias[refframe])
    {
        xmv.row *= -1;
        xmv.col *= -1;
    }

    mvp->as_mv = xmv;
}

#define LEFT_TOP_MARGIN (16 << 3)
#define RIGHT_BOTTOM_MARGIN (16 << 3)
static void vp8_clamp_mv(MV *mv, const MACROBLOCKD *xd)
{
    if (mv->col < (xd->mb_to_left_edge - LEFT_TOP_MARGIN))
        mv->col = xd->mb_to_left_edge - LEFT_TOP_MARGIN;
    else if (mv->col > xd->mb_to_right_edge + RIGHT_BOTTOM_MARGIN)
        mv->col = xd->mb_to_right_edge + RIGHT_BOTTOM_MARGIN;

    if (mv->row < (xd->mb_to_top_edge - LEFT_TOP_MARGIN))
        mv->row = xd->mb_to_top_edge - LEFT_TOP_MARGIN;
    else if (mv->row > xd->mb_to_bottom_edge + RIGHT_BOTTOM_MARGIN)
        mv->row = xd->mb_to_bottom_edge + RIGHT_BOTTOM_MARGIN;
}

void vp8_find_near_mvs
(
    MACROBLOCKD *xd,
    const MODE_INFO *here,
    MV *nearest, MV *nearby, MV *best,
    int near_mv_ref_cts[4],
    int refframe,
    int *ref_frame_sign_bias
);

vp8_prob *vp8_mv_ref_probs(
    vp8_prob p[VP8_MVREFS-1], const int near_mv_ref_ct[4]
);

const B_MODE_INFO *vp8_left_bmi(const MODE_INFO *cur_mb, int b);

const B_MODE_INFO *vp8_above_bmi(const MODE_INFO *cur_mb, int b, int mi_stride);

extern const unsigned char vp8_mbsplit_offset[4][16];

#endif
