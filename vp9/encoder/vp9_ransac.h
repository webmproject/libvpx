/*
 *  Copyright (c) 2015 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be
 *  found  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#ifndef VP9_ENCODER_VP9_RANSAC_H
#define VP9_ENCODER_VP9_RANSAC_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <memory.h>

#include "vp9/common/vp9_motion_model.h"

typedef int (*ransacType)(double *matched_points, int npoints,
                          int *number_of_inliers, int *best_inlier_mask,
                          double *bestH);

int ransacHomography(double *matched_points, int npoints,
                     int *number_of_inliers, int *best_inlier_indices,
                     double *bestH);
int ransacAffine(double *matched_points, int npoints,
                 int *number_of_inliers, int *best_inlier_indices,
                 double *bestH);
int ransacRotZoom(double *matched_points, int npoints,
                  int *number_of_inliers, int *best_inlier_indices,
                  double *bestH);
int ransacTranslation(double *matched_points, int npoints,
                      int *number_of_inliers, int *best_inlier_indices,
                      double *bestH);

void MultiplyMat(double *m1, double *m2, double *res,
                 const int M1, const int N1, const int N2);
int PseudoInverse(double *inv, double *matx, const int M, const int N);

#endif  // VP9_ENCODER_VP9_RANSAC_H
