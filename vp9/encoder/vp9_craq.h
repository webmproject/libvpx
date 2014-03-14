/*
 *  Copyright (c) 2014 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */


#ifndef VP9_ENCODER_VP9_CRAQ_H_
#define VP9_ENCODER_VP9_CRAQ_H_

#include "vp9/encoder/vp9_onyx_int.h"

#ifdef __cplusplus
extern "C" {
#endif

// Check if we should turn off cyclic refresh based on bitrate condition.
static int apply_cyclic_refresh_bitrate(VP9_COMP *const cpi);

// Check if this coding block, of size bsize, should be considered for refresh
// (lower-qp coding).
static int candidate_refresh_aq(VP9_COMP *const cpi,
                                MODE_INFO *const mi,
                                int bsize,
                                int use_rd);

// Prior to coding a given prediction block, of size bsize at (mi_row, mi_col),
// check if we should reset the segment_id, and update the cyclic_refresh map
// and segmentation map.
void vp9_update_segment_aq(VP9_COMP *const cpi,
                           MODE_INFO *const mi,
                           int mi_row,
                           int mi_col,
                           int bsize,
                           int use_rd);

// Setup cyclic background refresh: set delta q and segmentation map.
void vp9_setup_cyclic_refresh_aq(VP9_COMP *const cpi);

#ifdef __cplusplus
}  // extern "C"
#endif

#endif  // VP9_ENCODER_VP9_CRAQ_H_
