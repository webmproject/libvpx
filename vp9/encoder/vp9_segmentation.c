/*
 *  Copyright (c) 2012 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */


#include "limits.h"
#include "vpx_mem/vpx_mem.h"
#include "vp9/encoder/vp9_segmentation.h"
#include "vp9/common/vp9_pred_common.h"

void vp9_update_gf_useage_maps(VP9_COMP *cpi, VP9_COMMON *cm, MACROBLOCK *x) {
  int mb_row, mb_col;

  MODE_INFO *this_mb_mode_info = cm->mi;

  x->gf_active_ptr = (signed char *)cpi->gf_active_flags;

  if ((cm->frame_type == KEY_FRAME) || (cm->refresh_golden_frame)) {
    // Reset Gf useage monitors
    vpx_memset(cpi->gf_active_flags, 1, (cm->mb_rows * cm->mb_cols));
    cpi->gf_active_count = cm->mb_rows * cm->mb_cols;
  } else {
    // for each macroblock row in image
    for (mb_row = 0; mb_row < cm->mb_rows; mb_row++) {
      // for each macroblock col in image
      for (mb_col = 0; mb_col < cm->mb_cols; mb_col++) {

        // If using golden then set GF active flag if not already set.
        // If using last frame 0,0 mode then leave flag as it is
        // else if using non 0,0 motion or intra modes then clear
        // flag if it is currently set
        if ((this_mb_mode_info->mbmi.ref_frame == GOLDEN_FRAME) ||
            (this_mb_mode_info->mbmi.ref_frame == ALTREF_FRAME)) {
          if (*(x->gf_active_ptr) == 0) {
            *(x->gf_active_ptr) = 1;
            cpi->gf_active_count++;
          }
        } else if ((this_mb_mode_info->mbmi.mode != ZEROMV) &&
                   *(x->gf_active_ptr)) {
          *(x->gf_active_ptr) = 0;
          cpi->gf_active_count--;
        }

        x->gf_active_ptr++;          // Step onto next entry
        this_mb_mode_info++;         // skip to next mb

      }

      // this is to account for the border
      this_mb_mode_info++;
    }
  }
}

void vp9_enable_segmentation(VP9_PTR ptr) {
  VP9_COMP *cpi = (VP9_COMP *)(ptr);

  // Set the appropriate feature bit
  cpi->mb.e_mbd.segmentation_enabled = 1;
  cpi->mb.e_mbd.update_mb_segmentation_map = 1;
  cpi->mb.e_mbd.update_mb_segmentation_data = 1;
}

void vp9_disable_segmentation(VP9_PTR ptr) {
  VP9_COMP *cpi = (VP9_COMP *)(ptr);

  // Clear the appropriate feature bit
  cpi->mb.e_mbd.segmentation_enabled = 0;
}

void vp9_set_segmentation_map(VP9_PTR ptr,
                              unsigned char *segmentation_map) {
  VP9_COMP *cpi = (VP9_COMP *)(ptr);

  // Copy in the new segmentation map
  vpx_memcpy(cpi->segmentation_map, segmentation_map,
             (cpi->common.mb_rows * cpi->common.mb_cols));

  // Signal that the map should be updated.
  cpi->mb.e_mbd.update_mb_segmentation_map = 1;
  cpi->mb.e_mbd.update_mb_segmentation_data = 1;
}

void vp9_set_segment_data(VP9_PTR ptr,
                          signed char *feature_data,
                          unsigned char abs_delta) {
  VP9_COMP *cpi = (VP9_COMP *)(ptr);

  cpi->mb.e_mbd.mb_segment_abs_delta = abs_delta;

  vpx_memcpy(cpi->mb.e_mbd.segment_feature_data, feature_data,
             sizeof(cpi->mb.e_mbd.segment_feature_data));

  // TBD ?? Set the feature mask
  // vpx_memcpy(cpi->mb.e_mbd.segment_feature_mask, 0,
  //            sizeof(cpi->mb.e_mbd.segment_feature_mask));
}

// Based on set of segment counts calculate a probability tree
static void calc_segtree_probs(MACROBLOCKD *xd,
                               int *segcounts,
                               vp9_prob *segment_tree_probs) {
  int count1, count2;
  int tot_count;
  int i;

  // Blank the strtucture to start with
  vpx_memset(segment_tree_probs, 0,
             MB_FEATURE_TREE_PROBS * sizeof(*segment_tree_probs));

  // Total count for all segments
  count1 = segcounts[0] + segcounts[1];
  count2 = segcounts[2] + segcounts[3];
  tot_count = count1 + count2;

  // Work out probabilities of each segment
  if (tot_count)
    segment_tree_probs[0] = (count1 * 255) / tot_count;
  if (count1 > 0)
    segment_tree_probs[1] = (segcounts[0] * 255) / count1;
  if (count2 > 0)
    segment_tree_probs[2] = (segcounts[2] * 255) / count2;

  // Clamp probabilities to minimum allowed value
  for (i = 0; i < MB_FEATURE_TREE_PROBS; i++) {
    if (segment_tree_probs[i] == 0)
      segment_tree_probs[i] = 1;
  }
}

// Based on set of segment counts and probabilities calculate a cost estimate
static int cost_segmap(MACROBLOCKD *xd,
                       int *segcounts,
                       vp9_prob *probs) {
  int cost;
  int count1, count2;

  // Cost the top node of the tree
  count1 = segcounts[0] + segcounts[1];
  count2 = segcounts[2] + segcounts[3];
  cost = count1 * vp9_cost_zero(probs[0]) +
         count2 * vp9_cost_one(probs[0]);

  // Now add the cost of each individual segment branch
  if (count1 > 0)
    cost += segcounts[0] * vp9_cost_zero(probs[1]) +
            segcounts[1] * vp9_cost_one(probs[1]);

  if (count2 > 0)
    cost += segcounts[2] * vp9_cost_zero(probs[2]) +
            segcounts[3] * vp9_cost_one(probs[2]);

  return cost;

}

void vp9_choose_segmap_coding_method(VP9_COMP *cpi) {
  VP9_COMMON *const cm = &cpi->common;
  MACROBLOCKD *const xd = &cpi->mb.e_mbd;

  int i;
  int tot_count;
  int no_pred_cost;
  int t_pred_cost = INT_MAX;
  int pred_context;

  int mb_row, mb_col;
  int segmap_index = 0;
  unsigned char segment_id;

  int temporal_predictor_count[PREDICTION_PROBS][2];
  int no_pred_segcounts[MAX_MB_SEGMENTS];
  int t_unpred_seg_counts[MAX_MB_SEGMENTS];

  vp9_prob no_pred_tree[MB_FEATURE_TREE_PROBS];
  vp9_prob t_pred_tree[MB_FEATURE_TREE_PROBS];
  vp9_prob t_nopred_prob[PREDICTION_PROBS];

#if CONFIG_SUPERBLOCKS
  const int mis = cm->mode_info_stride;
#endif

  // Set default state for the segment tree probabilities and the
  // temporal coding probabilities
  vpx_memset(xd->mb_segment_tree_probs, 255,
             sizeof(xd->mb_segment_tree_probs));
  vpx_memset(cm->segment_pred_probs, 255,
             sizeof(cm->segment_pred_probs));

  vpx_memset(no_pred_segcounts, 0, sizeof(no_pred_segcounts));
  vpx_memset(t_unpred_seg_counts, 0, sizeof(t_unpred_seg_counts));
  vpx_memset(temporal_predictor_count, 0, sizeof(temporal_predictor_count));

  // First of all generate stats regarding how well the last segment map
  // predicts this one

  // Initialize macroblock decoder mode info context for the first mb
  // in the frame
  xd->mode_info_context = cm->mi;

  for (mb_row = 0; mb_row < cm->mb_rows; mb_row += 2) {
    for (mb_col = 0; mb_col < cm->mb_cols; mb_col += 2) {
      for (i = 0; i < 4; i++) {
        static const int dx[4] = { +1, -1, +1, +1 };
        static const int dy[4] = {  0, +1,  0, -1 };
        int x_idx = i & 1, y_idx = i >> 1;

        if (mb_col + x_idx >= cm->mb_cols ||
            mb_row + y_idx >= cm->mb_rows) {
          goto end;
        }

        xd->mb_to_top_edge = -((mb_row * 16) << 3);
        xd->mb_to_left_edge = -((mb_col * 16) << 3);

        segmap_index = (mb_row + y_idx) * cm->mb_cols + mb_col + x_idx;
        segment_id = xd->mode_info_context->mbmi.segment_id;
#if CONFIG_SUPERBLOCKS
        if (xd->mode_info_context->mbmi.encoded_as_sb) {
          if (mb_col + 1 < cm->mb_cols)
            segment_id = segment_id &&
                         xd->mode_info_context[1].mbmi.segment_id;
          if (mb_row + 1 < cm->mb_rows) {
            segment_id = segment_id &&
                         xd->mode_info_context[mis].mbmi.segment_id;
            if (mb_col + 1 < cm->mb_cols)
              segment_id = segment_id &&
                           xd->mode_info_context[mis + 1].mbmi.segment_id;
          }
          xd->mb_to_bottom_edge = ((cm->mb_rows - 2 - mb_row) * 16) << 3;
          xd->mb_to_right_edge  = ((cm->mb_cols - 2 - mb_col) * 16) << 3;
        } else {
#endif
          xd->mb_to_bottom_edge = ((cm->mb_rows - 1 - mb_row) * 16) << 3;
          xd->mb_to_right_edge  = ((cm->mb_cols - 1 - mb_col) * 16) << 3;
#if CONFIG_SUPERBLOCKS
        }
#endif

        // Count the number of hits on each segment with no prediction
        no_pred_segcounts[segment_id]++;

        // Temporal prediction not allowed on key frames
        if (cm->frame_type != KEY_FRAME) {
          // Test to see if the segment id matches the predicted value.
          int seg_predicted =
            (segment_id == vp9_get_pred_mb_segid(cm, xd, segmap_index));

          // Get the segment id prediction context
          pred_context =
            vp9_get_pred_context(cm, xd, PRED_SEG_ID);

          // Store the prediction status for this mb and update counts
          // as appropriate
          vp9_set_pred_flag(xd, PRED_SEG_ID, seg_predicted);
          temporal_predictor_count[pred_context][seg_predicted]++;

          if (!seg_predicted)
            // Update the "unpredicted" segment count
            t_unpred_seg_counts[segment_id]++;
        }

#if CONFIG_SUPERBLOCKS
        if (xd->mode_info_context->mbmi.encoded_as_sb) {
          assert(!i);
          xd->mode_info_context += 2;
          break;
        }
#endif
      end:
        xd->mode_info_context += dx[i] + dy[i] * cm->mode_info_stride;
      }
    }

    // this is to account for the border in mode_info_context
    xd->mode_info_context -= mb_col;
    xd->mode_info_context += cm->mode_info_stride * 2;
  }

  // Work out probability tree for coding segments without prediction
  // and the cost.
  calc_segtree_probs(xd, no_pred_segcounts, no_pred_tree);
  no_pred_cost = cost_segmap(xd, no_pred_segcounts, no_pred_tree);

  // Key frames cannot use temporal prediction
  if (cm->frame_type != KEY_FRAME) {
    // Work out probability tree for coding those segments not
    // predicted using the temporal method and the cost.
    calc_segtree_probs(xd, t_unpred_seg_counts, t_pred_tree);
    t_pred_cost = cost_segmap(xd, t_unpred_seg_counts, t_pred_tree);

    // Add in the cost of the signalling for each prediction context
    for (i = 0; i < PREDICTION_PROBS; i++) {
      tot_count = temporal_predictor_count[i][0] +
                  temporal_predictor_count[i][1];

      // Work out the context probabilities for the segment
      // prediction flag
      if (tot_count) {
        t_nopred_prob[i] = (temporal_predictor_count[i][0] * 255) /
                           tot_count;

        // Clamp to minimum allowed value
        if (t_nopred_prob[i] < 1)
          t_nopred_prob[i] = 1;
      } else
        t_nopred_prob[i] = 1;

      // Add in the predictor signaling cost
      t_pred_cost += (temporal_predictor_count[i][0] *
                      vp9_cost_zero(t_nopred_prob[i])) +
                     (temporal_predictor_count[i][1] *
                      vp9_cost_one(t_nopred_prob[i]));
    }
  }

  // Now choose which coding method to use.
  if (t_pred_cost < no_pred_cost) {
    cm->temporal_update = 1;
    vpx_memcpy(xd->mb_segment_tree_probs,
               t_pred_tree, sizeof(t_pred_tree));
    vpx_memcpy(&cm->segment_pred_probs,
               t_nopred_prob, sizeof(t_nopred_prob));
  } else {
    cm->temporal_update = 0;
    vpx_memcpy(xd->mb_segment_tree_probs,
               no_pred_tree, sizeof(no_pred_tree));
  }
}
