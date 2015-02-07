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

typedef int (*ransacType)(double *matched_points, int npoints,
                          int *number_of_inliers, int *best_inlier_mask,
                          double *bestH);
typedef void (*projectPointsType)(double *mat, double *points, double *proj,
                                  const int n, const int stride_points,
                                  const int stride_proj);

int ransacHomography(double *matched_points, int npoints,
                     int *number_of_inliers, int *best_inlier_indices,
                     double *bestH);
int ransacAffine(double *matched_points, int npoints,
                 int *number_of_inliers, int *best_inlier_indices,
                 double *bestH);
int ransacRotZoom(double *matched_points, int npoints,
                  int *number_of_inliers, int *best_inlier_indices,
                  double *bestH);

void projectPointsHomography(double *mat, double *points, double *proj,
                             const int n, const int stride_points, const int stride_proj);
void projectPointsAffine(double *mat, double *points, double *proj,
                         const int n, const int stride_points, const int stride_proj);
void projectPointsRotZoom(double *mat, double *points, double *proj,
                          const int n, const int stride_points, const int stride_proj);

#endif  // VP9_ENCODER_VP9_RANSAC_H
