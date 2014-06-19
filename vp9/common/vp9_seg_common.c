/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include <assert.h>

#include "vp9/common/vp9_blockd.h"
#include "vp9/common/vp9_loopfilter.h"
#include "vp9/common/vp9_seg_common.h"
#include "vp9/common/vp9_quant_common.h"


static const int seg_feature_data_signed[SEG_LVL_MAX] = { 1, 1, 0, 0 };

static const int seg_feature_data_max[SEG_LVL_MAX] = {
  MAXQ, MAX_LOOP_FILTER, 3, 0 };

#if CONFIG_VP9_HIGH && CONFIG_HIGH_TRANSFORMS && CONFIG_HIGH_QUANT
static const int seg_feature_data_max_10[SEG_LVL_MAX] = {
  MAXQ_10, MAX_LOOP_FILTER, 3, 0 };

static const int seg_feature_data_max_12[SEG_LVL_MAX] = {
  MAXQ_12, MAX_LOOP_FILTER, 3, 0 };
#endif

// These functions provide access to new segment level features.
// Eventually these function may be "optimized out" but for the moment,
// the coding mechanism is still subject to change so these provide a
// convenient single point of change.

int vp9_segfeature_active(const struct segmentation *seg, int segment_id,
                          SEG_LVL_FEATURES feature_id) {
  return seg->enabled &&
         (seg->feature_mask[segment_id] & (1 << feature_id));
}

void vp9_clearall_segfeatures(struct segmentation *seg) {
  vp9_zero(seg->feature_data);
  vp9_zero(seg->feature_mask);
}

void vp9_enable_segfeature(struct segmentation *seg, int segment_id,
                           SEG_LVL_FEATURES feature_id) {
  seg->feature_mask[segment_id] |= 1 << feature_id;
}

int vp9_seg_feature_data_max(SEG_LVL_FEATURES feature_id,
                             vpx_bit_depth_t bit_depth) {
#if CONFIG_VP9_HIGH && CONFIG_HIGH_TRANSFORMS && CONFIG_HIGH_QUANT
  switch (bit_depth) {
    case VPX_BITS_8:
      return seg_feature_data_max[feature_id];
    case VPX_BITS_10:
      return seg_feature_data_max_10[feature_id];
    case VPX_BITS_12:
      return seg_feature_data_max_12[feature_id];
    default:
      assert(0 && "bit_depth should be VPX_BITS_8, VPX_BITS_10 or VPX_BITS_12");
  }
#else
  (void) bit_depth;
  return seg_feature_data_max[feature_id];
#endif
}

int vp9_is_segfeature_signed(SEG_LVL_FEATURES feature_id) {
  return seg_feature_data_signed[feature_id];
}

void vp9_set_segdata(struct segmentation *seg, int segment_id,
                     SEG_LVL_FEATURES feature_id, int seg_data,
                     vpx_bit_depth_t bit_depth) {
  const int data_max = vp9_seg_feature_data_max(feature_id, bit_depth);
  assert(seg_data <= data_max);
  if (seg_data < 0) {
    assert(seg_feature_data_signed[feature_id]);
    assert(-seg_data <= data_max);
  }

  seg->feature_data[segment_id][feature_id] = seg_data;
}

int vp9_get_segdata(const struct segmentation *seg, int segment_id,
                    SEG_LVL_FEATURES feature_id) {
  return seg->feature_data[segment_id][feature_id];
}


const vp9_tree_index vp9_segment_tree[TREE_SIZE(MAX_SEGMENTS)] = {
  2,  4,  6,  8, 10, 12,
  0, -1, -2, -3, -4, -5, -6, -7
};


// TBD? Functions to read and write segment data with range / validity checking
