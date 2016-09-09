/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#ifndef AV1_ENCODER_GLOBAL_MOTION_H_
#define AV1_ENCODER_GLOBAL_MOTION_H_

#include "aom/aom_integer.h"

#ifdef __cplusplus
extern "C" {
#endif

/*
  Computes global motion parameters between two frames. The array
  "params" should be length 9, where the first 2 slots are translation
  parameters in (row, col) order, and the remaining slots correspond
  to values in the transformation matrix of the corresponding motion
  model. They are arranged in "params" such that values on the tx-matrix
  diagonal have odd numbered indices so the folowing matrix:
  A | B
  C | D
  would produce params = [trans row, trans col, B, A, C, D]
*/
int compute_global_motion_feature_based(TransformationType type,
                                        YV12_BUFFER_CONFIG *frm,
                                        YV12_BUFFER_CONFIG *ref,
                                        double *params);
#ifdef __cplusplus
}  // extern "C"
#endif
#endif  // AV1_ENCODER_GLOBAL_MOTION_H_
