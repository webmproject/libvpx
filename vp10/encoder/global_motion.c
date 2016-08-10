/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <math.h>
#include <assert.h>

#include "vp10/common/warped_motion.h"

#include "vp10/encoder/segmentation.h"
#include "vp10/encoder/global_motion.h"

int compute_global_motion_feature_based(struct VP10_COMP *cpi,
                                        TransformationType type,
                                        YV12_BUFFER_CONFIG *frm,
                                        YV12_BUFFER_CONFIG *ref,
                                        double inlier_prob, double *H) {
  (void) cpi;
  (void) type;
  (void) frm;
  (void) ref;
  (void) inlier_prob;
  (void) H;
  return 0;
}
