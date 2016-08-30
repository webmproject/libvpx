/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#ifndef AV1_ENCODER_RANSAC_H_
#define AV1_ENCODER_RANSAC_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <memory.h>

#include "av1/common/warped_motion.h"

typedef int (*RansacType)(double *matched_points, int npoints,
                          int *number_of_inliers, int *best_inlier_mask,
                          double *bestH);

/* Each of these functions fits a motion model from a set of
corresponding points in 2 frames using RANSAC.*/
int ransacHomography(double *matched_points, int npoints,
                     int *number_of_inliers, int *best_inlier_indices,
                     double *bestH);
int ransacAffine(double *matched_points, int npoints, int *number_of_inliers,
                 int *best_inlier_indices, double *bestH);
int ransacRotZoom(double *matched_points, int npoints, int *number_of_inliers,
                  int *best_inlier_indices, double *bestH);
int ransacTranslation(double *matched_points, int npoints,
                      int *number_of_inliers, int *best_inlier_indices,
                      double *bestH);
#endif  // AV1_ENCODER_RANSAC_H
