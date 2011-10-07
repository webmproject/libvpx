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

const int segfeaturedata_signed[SEG_LVL_MAX] = {1, 1, 0, 0, 0, 0};


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

void clearall_segfeatures( MACROBLOCKD *xd )
{
     vpx_memset(xd->segment_feature_data, 0, sizeof(xd->segment_feature_data));
     vpx_memset(xd->segment_feature_mask, 0, sizeof(xd->segment_feature_mask));
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

int is_segfeature_signed( SEG_LVL_FEATURES feature_id )
{
    return ( segfeaturedata_signed[feature_id] );
}

// TBD? Functions to read and write segment data with range / validity checking