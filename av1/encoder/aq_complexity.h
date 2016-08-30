/*
 *  Copyright (c) 2014 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#ifndef AV1_ENCODER_AQ_COMPLEXITY_H_
#define AV1_ENCODER_AQ_COMPLEXITY_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "av1/common/enums.h"

struct AV1_COMP;
struct macroblock;

// Select a segment for the current Block.
void av1_caq_select_segment(struct AV1_COMP *cpi, struct macroblock *,
                            BLOCK_SIZE bs, int mi_row, int mi_col,
                            int projected_rate);

// This function sets up a set of segments with delta Q values around
// the baseline frame quantizer.
void av1_setup_in_frame_q_adj(struct AV1_COMP *cpi);

#ifdef __cplusplus
}  // extern "C"
#endif

#endif  // AV1_ENCODER_AQ_COMPLEXITY_H_
