/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */


#include "segmentation.h"
#include "vpx_mem/vpx_mem.h"

void vp8_update_gf_useage_maps(VP8_COMP *cpi, VP8_COMMON *cm, MACROBLOCK *x)
{
    int mb_row, mb_col;

    MODE_INFO *this_mb_mode_info = cm->mi;

    x->gf_active_ptr = (signed char *)cpi->gf_active_flags;

    if ((cm->frame_type == KEY_FRAME) || (cm->refresh_golden_frame))
    {
        // Reset Gf useage monitors
        vpx_memset(cpi->gf_active_flags, 1, (cm->mb_rows * cm->mb_cols));
        cpi->gf_active_count = cm->mb_rows * cm->mb_cols;
    }
    else
    {
        // for each macroblock row in image
        for (mb_row = 0; mb_row < cm->mb_rows; mb_row++)
        {
            // for each macroblock col in image
            for (mb_col = 0; mb_col < cm->mb_cols; mb_col++)
            {

                // If using golden then set GF active flag if not already set.
                // If using last frame 0,0 mode then leave flag as it is
                // else if using non 0,0 motion or intra modes then clear
                // flag if it is currently set
                if ((this_mb_mode_info->mbmi.ref_frame == GOLDEN_FRAME) ||
                    (this_mb_mode_info->mbmi.ref_frame == ALTREF_FRAME))
                {
                    if (*(x->gf_active_ptr) == 0)
                    {
                        *(x->gf_active_ptr) = 1;
                        cpi->gf_active_count ++;
                    }
                }
                else if ((this_mb_mode_info->mbmi.mode != ZEROMV) &&
                         *(x->gf_active_ptr))
                {
                    *(x->gf_active_ptr) = 0;
                    cpi->gf_active_count--;
                }

                x->gf_active_ptr++;          // Step onto next entry
                this_mb_mode_info++;           // skip to next mb

            }

            // this is to account for the border
            this_mb_mode_info++;
        }
    }
}

void vp8_enable_segmentation(VP8_PTR ptr)
{
    VP8_COMP *cpi = (VP8_COMP *)(ptr);

    // Set the appropriate feature bit
    cpi->mb.e_mbd.segmentation_enabled = 1;
    cpi->mb.e_mbd.update_mb_segmentation_map = 1;
    cpi->mb.e_mbd.update_mb_segmentation_data = 1;
}

void vp8_disable_segmentation(VP8_PTR ptr)
{
    VP8_COMP *cpi = (VP8_COMP *)(ptr);

    // Clear the appropriate feature bit
    cpi->mb.e_mbd.segmentation_enabled = 0;
}

void vp8_set_segmentation_map(VP8_PTR ptr,
                              unsigned char *segmentation_map)
{
    VP8_COMP *cpi = (VP8_COMP *)(ptr);

    // Copy in the new segmentation map
    vpx_memcpy( cpi->segmentation_map, segmentation_map,
                (cpi->common.mb_rows * cpi->common.mb_cols) );

    // Signal that the map should be updated.
    cpi->mb.e_mbd.update_mb_segmentation_map = 1;
    cpi->mb.e_mbd.update_mb_segmentation_data = 1;
}

void vp8_set_segment_data(VP8_PTR ptr,
                          signed char *feature_data,
                          unsigned char abs_delta)
{
    VP8_COMP *cpi = (VP8_COMP *)(ptr);

    cpi->mb.e_mbd.mb_segement_abs_delta = abs_delta;

    vpx_memcpy(cpi->mb.e_mbd.segment_feature_data, feature_data,
               sizeof(cpi->mb.e_mbd.segment_feature_data));

//#if CONFIG_SEGFEATURES
    // TBD ?? Set the feature mask
    // vpx_memcpy(cpi->mb.e_mbd.segment_feature_mask, 0,
    //            sizeof(cpi->mb.e_mbd.segment_feature_mask));
}

#if CONFIG_SEGMENTATION
void choose_segmap_coding_method( VP8_COMP *cpi,
                                  int * segment_counts )
{
    VP8_COMMON *const cm = & cpi->common;
    MACROBLOCKD *const xd = & cpi->mb.e_mbd;

    int tot_count;
    int i;
    int count1,count2,count3,count4;
    int prob[3];
    int new_cost, original_cost;

    // Select the coding strategy for the segment map (temporal or spatial)
    tot_count = segment_counts[12] + segment_counts[13] +
                segment_counts[14] + segment_counts[15];
    count1 = segment_counts[12] + segment_counts[13];
    count2 = segment_counts[14] + segment_counts[15];

    if (tot_count)
        prob[0] = (count1 * 255) / tot_count;

    if (count1 > 0)
        prob[1] = (segment_counts[12] * 255) /count1;

    if (count2 > 0)
        prob[2] = (segment_counts[14] * 255) /count2;

    if (cm->frame_type != KEY_FRAME)
    {
        tot_count = segment_counts[4] + segment_counts[7];
        if (tot_count)
            xd->mb_segment_tree_probs[3] = (segment_counts[4] * 255)/tot_count;

        tot_count = segment_counts[5] + segment_counts[8];
        if (tot_count)
            xd->mb_segment_tree_probs[4] = (segment_counts[5] * 255)/tot_count;

        tot_count = segment_counts[6] + segment_counts[9];
        if (tot_count)
            xd->mb_segment_tree_probs[5] = (segment_counts[6] * 255)/tot_count;
    }

    tot_count = segment_counts[0] + segment_counts[1] +
                segment_counts[2] + segment_counts[3];
    count3 = segment_counts[0] + segment_counts[1];
    count4 = segment_counts[2] + segment_counts[3];

    if (tot_count)
        xd->mb_segment_tree_probs[0] = (count3 * 255) / tot_count;

    if (count3 > 0)
        xd->mb_segment_tree_probs[1] = (segment_counts[0] * 255) /count3;

    if (count4 > 0)
        xd->mb_segment_tree_probs[2] = (segment_counts[2] * 255) /count4;

    for (i = 0; i < MB_FEATURE_TREE_PROBS+3; i++)
    {
        if (xd->mb_segment_tree_probs[i] == 0)
            xd->mb_segment_tree_probs[i] = 1;
    }

    original_cost = count1 * vp8_cost_zero(prob[0]) +
                    count2 * vp8_cost_one(prob[0]);

    if (count1 > 0)
        original_cost += segment_counts[12] * vp8_cost_zero(prob[1]) +
                         segment_counts[13] * vp8_cost_one(prob[1]);

    if (count2 > 0)
        original_cost += segment_counts[14] * vp8_cost_zero(prob[2]) +
                         segment_counts[15] * vp8_cost_one(prob[2]) ;

    new_cost = 0;

    if (cm->frame_type != KEY_FRAME)
    {
        new_cost = segment_counts[4] *
                        vp8_cost_zero(xd->mb_segment_tree_probs[3]) +
                   segment_counts[7] *
                        vp8_cost_one(xd->mb_segment_tree_probs[3]);

        new_cost += segment_counts[5] *
                        vp8_cost_zero(xd->mb_segment_tree_probs[4]) +
                    segment_counts[8] *
                        vp8_cost_one(xd->mb_segment_tree_probs[4]);

        new_cost += segment_counts[6] *
                        vp8_cost_zero(xd->mb_segment_tree_probs[5]) +
                    segment_counts[9] *
                        vp8_cost_one (xd->mb_segment_tree_probs[5]);
    }

    if (tot_count > 0)
        new_cost += count3 * vp8_cost_zero(xd->mb_segment_tree_probs[0]) +
                    count4 * vp8_cost_one(xd->mb_segment_tree_probs[0]);

    if (count3 > 0)
        new_cost += segment_counts[0] *
                        vp8_cost_zero(xd->mb_segment_tree_probs[1]) +
                    segment_counts[1] *
                        vp8_cost_one(xd->mb_segment_tree_probs[1]);

    if (count4 > 0)
        new_cost += segment_counts[2] *
                        vp8_cost_zero(xd->mb_segment_tree_probs[2]) +
                    segment_counts[3] *
                        vp8_cost_one(xd->mb_segment_tree_probs[2]) ;

    if (new_cost < original_cost)
        xd->temporal_update = 1;
    else
    {
        xd->temporal_update = 0;
        xd->mb_segment_tree_probs[0] = prob[0];
        xd->mb_segment_tree_probs[1] = prob[1];
        xd->mb_segment_tree_probs[2] = prob[2];
    }

    // ***** TODO
    // PGW temp test code fix value as spatial
    xd->temporal_update = 0;
}
#endif
