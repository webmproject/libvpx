/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */


#include "vp8/common/blockd.h"
#include "onyx_int.h"
#include "treewriter.h"
#include "vp8/common/entropymode.h"


void vp8_init_mode_costs(VP8_COMP *c)
{
    VP8_COMMON *x = &c->common;
    {
        const vp8_tree_p T = vp8_bmode_tree;

        int i = 0;

        do
        {
            int j = 0;

            do
            {
                vp8_cost_tokens((int *)c->mb.bmode_costs[i][j], x->kf_bmode_prob[i][j], T);
            }
            while (++j < VP8_BINTRAMODES);
        }
        while (++i < VP8_BINTRAMODES);

        vp8_cost_tokens((int *)c->mb.inter_bmode_costs, x->fc.bmode_prob, T);
    }
#if CONFIG_ADAPTIVE_ENTROPY
    vp8_cost_tokens((int *)c->mb.inter_bmode_costs,
        vp8_sub_mv_ref_prob, vp8_sub_mv_ref_tree);
#else
    vp8_cost_tokens(c->mb.inter_bmode_costs,
        vp8_sub_mv_ref_prob, vp8_sub_mv_ref_tree);
#endif

    vp8_cost_tokens(c->mb.mbmode_cost[1], x->fc.ymode_prob, vp8_ymode_tree);
    vp8_cost_tokens(c->mb.mbmode_cost[0],
                    x->kf_ymode_prob[c->common.kf_ymode_probs_index],
                    vp8_kf_ymode_tree);
    vp8_cost_tokens(c->mb.intra_uv_mode_cost[1],
                    x->fc.uv_mode_prob[VP8_YMODES-1], vp8_uv_mode_tree);
    vp8_cost_tokens(c->mb.intra_uv_mode_cost[0],
                    x->kf_uv_mode_prob[VP8_YMODES-1], vp8_uv_mode_tree);
    vp8_cost_tokens(c->mb.i8x8_mode_costs,
                    x->fc.i8x8_mode_prob,vp8_i8x8_mode_tree);

}
