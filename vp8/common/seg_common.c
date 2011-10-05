/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include "vp8/common/seg_common.h"

// These functions provide access to new segment level features.
// Eventually these function may be "optimized out" but for the moment,
// the coding mechanism is still subject to change so these provide a
// convenient single point of change.

int segfeature_active( MACROBLOCKD *xd,
                       int segment_id,
                       SEG_LVL_FEATURES feature_id )
{
    // Return true if mask bit set and segmentation enabled.
    return ( xd->segmentation_enabled &&
             ( xd->segment_feature_mask[segment_id] &
               (0x01 << feature_id) ) );
}

void enable_segfeature( MACROBLOCKD *xd,
                        int segment_id,
                        SEG_LVL_FEATURES feature_id )
{
     xd->segment_feature_mask[segment_id] |= (0x01 << feature_id);
}
void disable_segfeature( MACROBLOCKD *xd,
                         int segment_id,
                         SEG_LVL_FEATURES feature_id )
{
     xd->segment_feature_mask[segment_id] &= ~(1 << feature_id);
}

// TBD? Functions to read and write segment data with range / validity checking