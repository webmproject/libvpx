/*
 *  Copyright (c) 2015 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be
 *  found  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#ifndef VP9_ENCODER_VP9_GLOBAL_MOTION_H
#define VP9_ENCODER_VP9_GLOBAL_MOTION_H

#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <math.h>
#include <assert.h>

#include "vp9_corner_detect.h"
#include "vp9_corner_match.h"
#include "vp9_ransac.h"
#include "vp9/common/vp9_motion_model.h"

struct VP9_COMP;

static const int CONFIDENCE_THRESHOLD = 1.0;

INLINE ransacType get_ransacType(TransformationType type);

// Searches around each parameter and seeks to minimize MSE  between
// the warped frame produced from the set of parameters and the frame being
// approximated.
void refine_quant_param(Global_Motion_Params *gm, TransformationType type,
                        unsigned char *ref, int ref_width,
                        int ref_height, int ref_stride,
                        unsigned char *frm, int frm_width, int frm_height,
                        int frm_stride, int n_refinements);

// Returns number of models actually returned: 1 - if success, 0 - if failure
int vp9_compute_global_motion_single_feature_based(struct VP9_COMP *cpi,
                                                   TransformationType type,
                                                   YV12_BUFFER_CONFIG *frm,
                                                   YV12_BUFFER_CONFIG *ref,
                                                   double *H);

// Returns number of models actually returned: 1+ - #models, 0 - if failure
// max_models is the maximum number of models returned
// inlier_prob is the probability of being inlier over all the models
// combined, beyond which no more models are computed.
// Ex. if max_models = 4, and inlier_prob = 0.8, and during the
// process three models together already cover more than 80% of the
// matching points, then only three models are returned.
int vp9_compute_global_motion_multiple_feature_based(struct VP9_COMP *cpi,
                                                     TransformationType type,
                                                     YV12_BUFFER_CONFIG *frm,
                                                     YV12_BUFFER_CONFIG *ref,
                                                     int max_models,
                                                     double inlier_prob,
                                                     double *H);

// Returns number of models actually returned: 1 - if success, 0 - if failure
int vp9_compute_global_motion_single_block_based(struct VP9_COMP *cpi,
                                                 TransformationType type,
                                                 YV12_BUFFER_CONFIG *frm,
                                                 YV12_BUFFER_CONFIG *ref,
                                                 BLOCK_SIZE bsize,
                                                 double *H);

int vp9_compute_global_motion_multiple_optical_flow(struct VP9_COMP *cpi,
                                                    TransformationType type,
                                                    YV12_BUFFER_CONFIG *frm,
                                                    YV12_BUFFER_CONFIG *ref,
                                                    int max_models,
                                                    double inlier_prob,
                                                    double *H);


// Returns number of models actually returned: 1+ - #models, 0 - if failure
// max_models is the maximum number of models returned
// inlier_prob is the probability of being inlier over all the models
// combined, beyond which no more models are computed.
// Ex. if max_models = 4, and inlier_prob = 0.8, and during the
// process three models together already cover more than 80% of the
// matching points, then only three models are returned.
int vp9_compute_global_motion_multiple_block_based(struct VP9_COMP *cpi,
                                                   TransformationType type,
                                                   YV12_BUFFER_CONFIG *frm,
                                                   YV12_BUFFER_CONFIG *ref,
                                                   BLOCK_SIZE bsize,
                                                   int max_models,
                                                   double inlier_prob,
                                                   double *H);

#endif  // VP9_ENCODER_VP9_GLOBAL_MOTION_H
