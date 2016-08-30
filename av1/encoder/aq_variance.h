/*
 *  Copyright (c) 2013 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#ifndef AV1_ENCODER_AQ_VARIANCE_H_
#define AV1_ENCODER_AQ_VARIANCE_H_

#include "av1/encoder/encoder.h"

#ifdef __cplusplus
extern "C" {
#endif

unsigned int av1_vaq_segment_id(int energy);
void av1_vaq_frame_setup(AV1_COMP *cpi);

int av1_block_energy(AV1_COMP *cpi, MACROBLOCK *x, BLOCK_SIZE bs);
double av1_log_block_var(AV1_COMP *cpi, MACROBLOCK *x, BLOCK_SIZE bs);

#ifdef __cplusplus
}  // extern "C"
#endif

#endif  // AV1_ENCODER_AQ_VARIANCE_H_
