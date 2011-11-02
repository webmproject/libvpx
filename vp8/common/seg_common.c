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

#if CONFIG_SEGFEATURES
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

void clear_segdata( MACROBLOCKD *xd,
                    int segment_id,
                    SEG_LVL_FEATURES feature_id)
{
    xd->segment_feature_data[segment_id][feature_id] = 0;
}

void set_segdata( MACROBLOCKD *xd,
                  int segment_id,
                  SEG_LVL_FEATURES feature_id,
                  int seg_data )
{
    xd->segment_feature_data[segment_id][feature_id] = seg_data;
}

int get_segdata( MACROBLOCKD *xd,
                 int segment_id,
                 SEG_LVL_FEATURES feature_id )
{
    return xd->segment_feature_data[segment_id][feature_id];
}

void clear_segref( MACROBLOCKD *xd, int segment_id )
{
    xd->segment_feature_data[segment_id][SEG_LVL_REF_FRAME] = 0;
}

void set_segref( MACROBLOCKD *xd,
                 int segment_id,
                 MV_REFERENCE_FRAME ref_frame )
{
    xd->segment_feature_data[segment_id][SEG_LVL_REF_FRAME] |=
        (1 << ref_frame);
}

int check_segref( MACROBLOCKD *xd,
                  int segment_id,
                  MV_REFERENCE_FRAME ref_frame )
{
    return ( xd->segment_feature_data[segment_id][SEG_LVL_REF_FRAME] &
             (1 << ref_frame) ) ? 1 : 0;
}

int check_segref_inter(MACROBLOCKD *xd, int segment_id)
{
    return ( xd->segment_feature_data[segment_id][SEG_LVL_REF_FRAME] &
             ~(1 << INTRA_FRAME) ) ? 1 : 0;
}

// TBD? Functions to read and write segment data with range / validity checking

#endif
