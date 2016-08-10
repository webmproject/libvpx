/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */


#ifndef VP10_ENCODER_GLOBAL_MOTION_H_
#define VP10_ENCODER_GLOBAL_MOTION_H_

#include "vpx/vpx_integer.h"

#ifdef __cplusplus
extern "C" {
#endif

int compute_global_motion_feature_based(struct VP10_COMP *cpi,
                                        TransformationType type,
                                        YV12_BUFFER_CONFIG *frm,
                                        YV12_BUFFER_CONFIG *ref,
                                        double inlier_prob, double *H);
#ifdef __cplusplus
}  // extern "C"
#endif
#endif  // VP10_ENCODER_GLOBAL_MOTION_H_

