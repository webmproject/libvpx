/*
 *  Copyright (c) 2012 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */


#include <limits.h>
#include "vpx_mem/vpx_mem.h"
#include "vp9/encoder/vp9_segmentation.h"
#include "vp9/common/vp9_pred_common.h"
#include "vp9/common/vp9_tile_common.h"

void vp9_enable_segmentation(VP9_PTR ptr) {
  VP9_COMP *cpi = (VP9_COMP *)ptr;
  struct segmentation *const seg =  &cpi->common.seg;

  seg->enabled = 1;
  seg->update_map = 1;
  seg->update_data = 1;
}

void vp9_disable_segmentation(VP9_PTR ptr) {
  VP9_COMP *cpi = (VP9_COMP *)ptr;
  struct segmentation *const seg =  &cpi->common.seg;
  seg->enabled = 0;
}

void vp9_set_segmentation_map(VP9_PTR ptr,
                              unsigned char *segmentation_map) {
  VP9_COMP *cpi = (VP9_COMP *)ptr;
  struct segmentation *const seg = &cpi->common.seg;

  // Copy in the new segmentation map
  vpx_memcpy(cpi->segmentation_map, segmentation_map,
             (cpi->common.mi_rows * cpi->common.mi_cols));

  // Signal that the map should be updated.
  seg->update_map = 1;
  seg->update_data = 1;
}

void vp9_set_segment_data(VP9_PTR ptr,
                          signed char *feature_data,
                          unsigned char abs_delta) {
  VP9_COMP *cpi = (VP9_COMP *)ptr;
  struct segmentation *const seg = &cpi->common.seg;

  seg->abs_delta = abs_delta;

  vpx_memcpy(seg->feature_data, feature_data, sizeof(seg->feature_data));

  // TBD ?? Set the feature mask
  // vpx_memcpy(cpi->mb.e_mbd.segment_feature_mask, 0,
  //            sizeof(cpi->mb.e_mbd.segment_feature_mask));
}

// Based on set of segment counts calculate a probability tree
static void calc_segtree_probs(int *segcounts, vp9_prob *segment_tree_probs) {
  // Work out probabilities of each segment
  const int c01 = segcounts[0] + segcounts[1];
  const int c23 = segcounts[2] + segcounts[3];
  const int c45 = segcounts[4] + segcounts[5];
  const int c67 = segcounts[6] + segcounts[7];

  segment_tree_probs[0] = get_binary_prob(c01 + c23, c45 + c67);
  segment_tree_probs[1] = get_binary_prob(c01, c23);
  segment_tree_probs[2] = get_binary_prob(c45, c67);
  segment_tree_probs[3] = get_binary_prob(segcounts[0], segcounts[1]);
  segment_tree_probs[4] = get_binary_prob(segcounts[2], segcounts[3]);
  segment_tree_probs[5] = get_binary_prob(segcounts[4], segcounts[5]);
  segment_tree_probs[6] = get_binary_prob(segcounts[6], segcounts[7]);
}

// Based on set of segment counts and probabilities calculate a cost estimate
static int cost_segmap(int *segcounts, vp9_prob *probs) {
  const int c01 = segcounts[0] + segcounts[1];
  const int c23 = segcounts[2] + segcounts[3];
  const int c45 = segcounts[4] + segcounts[5];
  const int c67 = segcounts[6] + segcounts[7];
  const int c0123 = c01 + c23;
  const int c4567 = c45 + c67;

  // Cost the top node of the tree
  int cost = c0123 * vp9_cost_zero(probs[0]) +
             c4567 * vp9_cost_one(probs[0]);

  // Cost subsequent levels
  if (c0123 > 0) {
    cost += c01 * vp9_cost_zero(probs[1]) +
            c23 * vp9_cost_one(probs[1]);

    if (c01 > 0)
      cost += segcounts[0] * vp9_cost_zero(probs[3]) +
              segcounts[1] * vp9_cost_one(probs[3]);
    if (c23 > 0)
      cost += segcounts[2] * vp9_cost_zero(probs[4]) +
              segcounts[3] * vp9_cost_one(probs[4]);
  }

  if (c4567 > 0) {
    cost += c45 * vp9_cost_zero(probs[2]) +
            c67 * vp9_cost_one(probs[2]);

    if (c45 > 0)
      cost += segcounts[4] * vp9_cost_zero(probs[5]) +
              segcounts[5] * vp9_cost_one(probs[5]);
    if (c67 > 0)
      cost += segcounts[6] * vp9_cost_zero(probs[6]) +
              segcounts[7] * vp9_cost_one(probs[6]);
  }

  return cost;
}

static void count_segs(VP9_COMP *cpi, MODE_INFO *mi,
                       int *no_pred_segcounts,
                       int (*temporal_predictor_count)[2],
                       int *t_unpred_seg_counts,
                       int bw, int bh, int mi_row, int mi_col) {
  VP9_COMMON *const cm = &cpi->common;
  MACROBLOCKD *const xd = &cpi->mb.e_mbd;
  int segment_id;

  if (mi_row >= cm->mi_rows || mi_col >= cm->mi_cols)
    return;

  segment_id = mi->mbmi.segment_id;
  xd->mode_info_context = mi;
  set_mi_row_col(cm, xd, mi_row, bh, mi_col, bw);

  // Count the number of hits on each segment with no prediction
  no_pred_segcounts[segment_id]++;

  // Temporal prediction not allowed on key frames
  if (cm->frame_type != KEY_FRAME) {
    const BLOCK_SIZE_TYPE bsize = mi->mbmi.sb_type;
    // Test to see if the segment id matches the predicted value.
    const int pred_segment_id = vp9_get_segment_id(cm, cm->last_frame_seg_map,
                                                   bsize, mi_row, mi_col);
    const int pred_flag = pred_segment_id == segment_id;
    const int pred_context = vp9_get_pred_context_seg_id(xd);

    // Store the prediction status for this mb and update counts
    // as appropriate
    vp9_set_pred_flag_seg_id(cm, bsize, mi_row, mi_col, pred_flag);
    temporal_predictor_count[pred_context][pred_flag]++;

    if (!pred_flag)
      // Update the "unpredicted" segment count
      t_unpred_seg_counts[segment_id]++;
  }
}

static void count_segs_sb(VP9_COMP *cpi, MODE_INFO *mi,
                          int *no_pred_segcounts,
                          int (*temporal_predictor_count)[2],
                          int *t_unpred_seg_counts,
                          int mi_row, int mi_col,
                          BLOCK_SIZE_TYPE bsize) {
  VP9_COMMON *const cm = &cpi->common;
  const int mis = cm->mode_info_stride;
  int bwl, bhl;
  const int bsl = mi_width_log2(bsize), bs = 1 << (bsl - 1);

  if (mi_row >= cm->mi_rows || mi_col >= cm->mi_cols)
    return;

  bwl = mi_width_log2(mi->mbmi.sb_type);
  bhl = mi_height_log2(mi->mbmi.sb_type);

  if (bwl == bsl && bhl == bsl) {
    count_segs(cpi, mi, no_pred_segcounts, temporal_predictor_count,
               t_unpred_seg_counts, 1 << bsl, 1 << bsl, mi_row, mi_col);
  } else if (bwl == bsl && bhl < bsl) {
    count_segs(cpi, mi, no_pred_segcounts, temporal_predictor_count,
               t_unpred_seg_counts, 1 << bsl, bs, mi_row, mi_col);
    count_segs(cpi, mi + bs * mis, no_pred_segcounts, temporal_predictor_count,
               t_unpred_seg_counts, 1 << bsl, bs, mi_row + bs, mi_col);
  } else if (bwl < bsl && bhl == bsl) {
    count_segs(cpi, mi, no_pred_segcounts, temporal_predictor_count,
               t_unpred_seg_counts, bs, 1 << bsl, mi_row, mi_col);
    count_segs(cpi, mi + bs, no_pred_segcounts, temporal_predictor_count,
               t_unpred_seg_counts, bs, 1 << bsl, mi_row, mi_col + bs);
  } else {
    BLOCK_SIZE_TYPE subsize;
    int n;

    assert(bwl < bsl && bhl < bsl);
    if (bsize == BLOCK_64X64) {
      subsize = BLOCK_32X32;
    } else if (bsize == BLOCK_32X32) {
      subsize = BLOCK_16X16;
    } else {
      assert(bsize == BLOCK_16X16);
      subsize = BLOCK_8X8;
    }

    for (n = 0; n < 4; n++) {
      const int y_idx = n >> 1, x_idx = n & 0x01;

      count_segs_sb(cpi, mi + y_idx * bs * mis + x_idx * bs,
                    no_pred_segcounts, temporal_predictor_count,
                    t_unpred_seg_counts,
                    mi_row + y_idx * bs, mi_col + x_idx * bs, subsize);
    }
  }
}

void vp9_choose_segmap_coding_method(VP9_COMP *cpi) {
  VP9_COMMON *const cm = &cpi->common;
  struct segmentation *seg = &cm->seg;

  int no_pred_cost;
  int t_pred_cost = INT_MAX;

  int i, tile_col, mi_row, mi_col;

  int temporal_predictor_count[PREDICTION_PROBS][2];
  int no_pred_segcounts[MAX_SEGMENTS];
  int t_unpred_seg_counts[MAX_SEGMENTS];

  vp9_prob no_pred_tree[SEG_TREE_PROBS];
  vp9_prob t_pred_tree[SEG_TREE_PROBS];
  vp9_prob t_nopred_prob[PREDICTION_PROBS];

  const int mis = cm->mode_info_stride;
  MODE_INFO *mi_ptr, *mi;

  // Set default state for the segment tree probabilities and the
  // temporal coding probabilities
  vpx_memset(seg->tree_probs, 255, sizeof(seg->tree_probs));
  vpx_memset(seg->pred_probs, 255, sizeof(seg->pred_probs));

  vpx_memset(no_pred_segcounts, 0, sizeof(no_pred_segcounts));
  vpx_memset(t_unpred_seg_counts, 0, sizeof(t_unpred_seg_counts));
  vpx_memset(temporal_predictor_count, 0, sizeof(temporal_predictor_count));

  // First of all generate stats regarding how well the last segment map
  // predicts this one
  for (tile_col = 0; tile_col < 1 << cm->log2_tile_cols; tile_col++) {
    vp9_get_tile_col_offsets(cm, tile_col);
    mi_ptr = cm->mi + cm->cur_tile_mi_col_start;
    for (mi_row = 0; mi_row < cm->mi_rows;
         mi_row += 8, mi_ptr += 8 * mis) {
      mi = mi_ptr;
      for (mi_col = cm->cur_tile_mi_col_start; mi_col < cm->cur_tile_mi_col_end;
           mi_col += 8, mi += 8)
        count_segs_sb(cpi, mi, no_pred_segcounts, temporal_predictor_count,
                      t_unpred_seg_counts, mi_row, mi_col, BLOCK_64X64);
    }
  }

  // Work out probability tree for coding segments without prediction
  // and the cost.
  calc_segtree_probs(no_pred_segcounts, no_pred_tree);
  no_pred_cost = cost_segmap(no_pred_segcounts, no_pred_tree);

  // Key frames cannot use temporal prediction
  if (cm->frame_type != KEY_FRAME) {
    // Work out probability tree for coding those segments not
    // predicted using the temporal method and the cost.
    calc_segtree_probs(t_unpred_seg_counts, t_pred_tree);
    t_pred_cost = cost_segmap(t_unpred_seg_counts, t_pred_tree);

    // Add in the cost of the signalling for each prediction context
    for (i = 0; i < PREDICTION_PROBS; i++) {
      const int count0 = temporal_predictor_count[i][0];
      const int count1 = temporal_predictor_count[i][1];

      t_nopred_prob[i] = get_binary_prob(count0, count1);

      // Add in the predictor signaling cost
      t_pred_cost += count0 * vp9_cost_zero(t_nopred_prob[i]) +
                     count1 * vp9_cost_one(t_nopred_prob[i]);
    }
  }

  // Now choose which coding method to use.
  if (t_pred_cost < no_pred_cost) {
    seg->temporal_update = 1;
    vpx_memcpy(seg->tree_probs, t_pred_tree, sizeof(t_pred_tree));
    vpx_memcpy(seg->pred_probs, t_nopred_prob, sizeof(t_nopred_prob));
  } else {
    seg->temporal_update = 0;
    vpx_memcpy(seg->tree_probs, no_pred_tree, sizeof(no_pred_tree));
  }
}
