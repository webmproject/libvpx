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

#include "av1/common/warped_motion.h"

#include "av1/encoder/segmentation.h"
#include "av1/encoder/global_motion.h"
#include "av1/encoder/corner_detect.h"
#include "av1/encoder/corner_match.h"
#include "av1/encoder/ransac.h"

#define MAX_CORNERS 4096
#define MIN_INLIER_PROB 0.1

static INLINE RansacFunc get_ransac_type(TransformationType type) {
  switch (type) {
    case HOMOGRAPHY: return ransac_homography;
    case AFFINE: return ransac_affine;
    case ROTZOOM: return ransac_rotzoom;
    case TRANSLATION: return ransac_translation;
    default: assert(0); return NULL;
  }
}

// computes global motion parameters by fitting a model using RANSAC
static int compute_global_motion_params(TransformationType type,
                                        double *correspondences,
                                        int num_correspondences, double *params,
                                        int *inlier_map) {
  int result;
  int num_inliers = 0;
  RansacFunc ransac = get_ransac_type(type);
  if (ransac == NULL) return 0;

  result = ransac(correspondences, num_correspondences, &num_inliers,
                  inlier_map, params);
  if (!result && num_inliers < MIN_INLIER_PROB * num_correspondences) {
    result = 1;
    num_inliers = 0;
  }
  return num_inliers;
}

int compute_global_motion_feature_based(TransformationType type,
                                        YV12_BUFFER_CONFIG *frm,
                                        YV12_BUFFER_CONFIG *ref,
                                        double *params) {
  int num_frm_corners, num_ref_corners;
  int num_correspondences;
  double *correspondences;
  int num_inliers;
  int frm_corners[2 * MAX_CORNERS], ref_corners[2 * MAX_CORNERS];
  int *inlier_map = NULL;

  // compute interest points in images using FAST features
  num_frm_corners =
      fast_corner_detect(frm->y_buffer, frm->y_width, frm->y_height,
                         frm->y_stride, frm_corners, MAX_CORNERS);
  num_ref_corners =
      fast_corner_detect(ref->y_buffer, ref->y_width, ref->y_height,
                         ref->y_stride, ref_corners, MAX_CORNERS);

  // find correspondences between the two images
  correspondences =
      (double *)malloc(num_frm_corners * 4 * sizeof(*correspondences));
  num_correspondences = determine_correspondence(
      frm->y_buffer, (int *)frm_corners, num_frm_corners, ref->y_buffer,
      (int *)ref_corners, num_ref_corners, frm->y_width, frm->y_height,
      frm->y_stride, ref->y_stride, correspondences);

  inlier_map = (int *)malloc(num_correspondences * sizeof(*inlier_map));
  num_inliers = compute_global_motion_params(
      type, correspondences, num_correspondences, params, inlier_map);
  free(correspondences);
  free(inlier_map);
  return (num_inliers > 0);
}
