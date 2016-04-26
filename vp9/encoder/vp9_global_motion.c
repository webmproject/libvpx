/*
 *  Copyright (c) 2015 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be
 *  found  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <math.h>
#include <assert.h>


#include "vpx_mem/vpx_mem.h"
#include "vp9/encoder/vp9_segmentation.h"
#include "vp9/encoder/vp9_mcomp.h"
#include "vp9/common/vp9_blockd.h"
#include "vp9/common/vp9_reconinter.h"
#include "vp9/common/vp9_reconintra.h"
#include "vp9/common/vp9_systemdependent.h"
#include "vp9/common/vp9_motion_model.h"
#include "vp9/encoder/vp9_corner_detect.h"
#include "vp9/encoder/vp9_corner_match.h"
#include "vp9/encoder/vp9_opticalflow.h"
#include "vp9/encoder/vp9_ransac.h"
#include "vp9/encoder/vp9_global_motion.h"
#include "vp9/encoder/vp9_motion_field.h"

// #define VERBOSE

// Default is Harris
#define USE_FAST_CORNER

#define MIN_INLIER_PROB 0.1

// #define PARAM_SEARCH

#define MAX_CORNERS 4096

INLINE ransacType get_ransacType(TransformationType type) {
  switch (type) {
    case HOMOGRAPHY:
      return ransacHomography;
    case AFFINE:
      return ransacAffine;
    case ROTZOOM:
      return ransacRotZoom;
    case TRANSLATION:
      return ransacTranslation;
    default:
      assert(0);
      return NULL;
  }
}

static double compute_error_score(TransformationType type,
                                  int *points1, int stride1,
                                  int *points2, int stride2,
                                  int npoints, double *H,
                                  int *map) {
  int i, n = 0;
  double ot[2], pt[2];
  int *mp1 = points1;
  int *mp2 = points2;
  double sqerr = 0.0;
  const int numparams = get_numparams(type);

  projectPointsType projectPoints = get_projectPointsType(type);
  if (projectPoints == NULL) return -1.0;
  if (map) {
    for (i = 0; i < npoints; ++i, mp1+=stride1, mp2+=stride2) {
      if (map[i] != -1) {
        ot[0] = mp1[0];
        ot[1] = mp1[1];
        projectPoints(&H[map[i] * numparams], ot, pt, 1, stride1, stride2);
        sqerr += (pt[0] - mp2[0]) * (pt[0] - mp2[0]) +
                 (pt[1] - mp2[1]) * (pt[1] - mp2[1]);
        n++;
      }
    }
  } else {
    for (i = 0; i < npoints; ++i, mp1+=stride1, mp2+=stride2) {
      ot[0] = mp1[0];
      ot[1] = mp1[1];
      projectPoints(H, ot, pt, 1, stride1, stride2);
      sqerr += (pt[0] - mp2[0]) * (pt[0] - mp2[0]) +
               (pt[1] - mp2[1]) * (pt[1] - mp2[1]);
      n++;
    }
  }
  return sqrt(sqerr / n);
}

static int16_t* get_gm_param_trans(Global_Motion_Params *gm, int index) {
  switch (index) {
    case 0 :
      return &gm->mv.as_mv.row;
      break;
    case 1 :
      return &gm->mv.as_mv.col;
      break;
  }
  assert(0);
  return NULL;
}

static int* get_gm_param_rotzoom(Global_Motion_Params *gm, int index) {
  switch (index) {
    case 0 :
      return &gm->rotation;
      break;
    case 1 :
      return &gm->zoom;
      break;
  }
  assert(0);
  return NULL;
}

static void refine_quant_param_trans(Global_Motion_Params *gm,
                                     TransformationType type,
                                     unsigned char *ref, int ref_width,
                                     int ref_height, int ref_stride,
                                     unsigned char *frm, int frm_width,
                                     int frm_height, int frm_stride,
                                     int n_refinements) {
  int i = 0, p;
  double step_mse;
  int step;
  int16_t *param;
  int16_t curr_param;
  int16_t best_param;
  projectPointsType projectPoints = get_projectPointsType(type);
  double best_mse = compute_warp_and_error(gm, projectPoints, ref, ref_width,
                                           ref_height, ref_stride, frm, 0, 0,
                                           frm_width, frm_height, frm_stride,
                                           0, 0, 16, 16);
  for (p = 0; p < 2; ++p) {
    param = get_gm_param_trans(gm, p);
    step = 1 << (n_refinements + 1);
    curr_param = *param;
    best_param = curr_param;
    for (i = 0; i < n_refinements; i++) {
      // look to the left
      *param = curr_param - step;
      step_mse = compute_warp_and_error(gm, projectPoints, ref, ref_width,
                                        ref_height, ref_stride, frm, 0, 0,
                                        frm_width, frm_height, frm_stride,
                                        0, 0, 16, 16);
      if (step_mse < best_mse) {
        step >>= 1;
        best_mse = step_mse;
        best_param = *param;
        curr_param = best_param;
        continue;
      }

      // look to the right
      *param = curr_param + step;
      step_mse = compute_warp_and_error(gm, projectPoints, ref, ref_width,
                                        ref_height, ref_stride, frm, 0, 0,
                                        frm_width, frm_height, frm_stride,
                                        0, 0, 16, 16);
      if (step_mse < best_mse) {
        step >>= 1;
        best_mse = step_mse;
        best_param = *param;
        curr_param = best_param;
        continue;
      }

      // no improvement found-> means we're either already at a minimum or
      // step is too wide
      step >>= 1;
    }

    *param = best_param;
  }
}

static void refine_quant_param_rotzoom(Global_Motion_Params *gm,
                                       TransformationType type,
                                       unsigned char *ref, int ref_width,
                                       int ref_height, int ref_stride,
                                       unsigned char *frm, int frm_width,
                                       int frm_height, int frm_stride,
                                       int n_refinements) {
  int i = 0, p;
  double step_mse;
  int step;
  int *param;
  int curr_param;
  int best_param;
  projectPointsType projectPoints = get_projectPointsType(type);
  double best_mse = compute_warp_and_error(gm, projectPoints, ref, ref_width,
                                           ref_height, ref_stride, frm, 0, 0,
                                           frm_width, frm_height, frm_stride,
                                           0, 0, 16, 16);
  for (p = 0; p < 2; ++p) {
    param = get_gm_param_rotzoom(gm, p);
    step = 1 << (n_refinements + 1);
    curr_param = *param;
    best_param = curr_param;
    for (i = 0; i < n_refinements; i++) {
      // look to the left
      *param = curr_param - step;
      step_mse = compute_warp_and_error(gm, projectPoints, ref, ref_width,
                                        ref_height, ref_stride, frm, 0, 0,
                                        frm_width, frm_height, frm_stride,
                                        0, 0, 16, 16);
      if (step_mse < best_mse) {
        step >>= 1;
        best_mse = step_mse;
        best_param = *param;
        curr_param = best_param;
        continue;
      }

      // look to the right
      *param = curr_param + step;
      step_mse = compute_warp_and_error(gm, projectPoints, ref, ref_width,
                                        ref_height, ref_stride, frm, 0, 0,
                                        frm_width, frm_height, frm_stride,
                                        0, 0, 16, 16);
      if (step_mse < best_mse) {
        step >>= 1;
        best_mse = step_mse;
        best_param = *param;
        curr_param = best_param;
        continue;
      }

      // no improvement found-> means we're either already at a minimum or
      // step is too wide
      step >>= 1;
    }

    *param = best_param;
  }
}

void refine_quant_param(Global_Motion_Params *gm,
                        TransformationType type,
                        unsigned char *ref, int ref_width,
                        int ref_height, int ref_stride,
                        unsigned char *frm, int frm_width, int frm_height,
                        int frm_stride, int n_refinements) {
  switch (gm->gmtype) {
    case GLOBAL_TRANSLATION :
      refine_quant_param_trans(gm, type, ref, ref_width, ref_height,
                                ref_stride, frm, frm_width, frm_height,
                                frm_stride, n_refinements);
      break;
    case GLOBAL_ROTZOOM :
      refine_quant_param_rotzoom(gm, type, ref, ref_width, ref_height,
                                 ref_stride, frm, frm_width, frm_height,
                                 frm_stride, n_refinements);
      refine_quant_param_trans(gm, type, ref, ref_width, ref_height,
                               ref_stride, frm, frm_width, frm_height,
                               frm_stride, n_refinements);
      break;
    default :
      break;
  }
}

static int compute_global_motion_single(TransformationType type,
                                        double *correspondences,
                                        int num_correspondences,
                                        double *H,
                                        int *inlier_map) {
  double *mp, *matched_points;
  double *cp = correspondences;
  int i, result;
  int num_inliers = 0;
  ransacType ransac = get_ransacType(type);
  if (ransac == NULL)
    return 0;
  matched_points =
      (double *)malloc(4 * num_correspondences * sizeof(double));

  for (mp = matched_points, cp = correspondences, i = 0;
       i < num_correspondences; ++i) {
    *mp++ = *cp++;
    *mp++ = *cp++;
    *mp++ = *cp++;
    *mp++ = *cp++;
  }
  result = ransac(matched_points, num_correspondences,
                  &num_inliers, inlier_map, H);
  if (!result && num_inliers < MIN_INLIER_PROB * num_correspondences) {
    result = 1;
    num_inliers = 0;
  }
  if (!result) {
    for (i = 0; i < num_correspondences; ++i) {
      inlier_map[i] = inlier_map[i] - 1;
    }
  }
  free(matched_points);
  return num_inliers;
}

// Returns number of models actually returned: 1 - if success, 0 - if failure
int vp9_compute_global_motion_single_feature_based(
    struct VP9_COMP *cpi,
    TransformationType type,
    YV12_BUFFER_CONFIG *frm,
    YV12_BUFFER_CONFIG *ref,
    double *H) {
  int num_frm_corners, num_ref_corners;
  int num_correspondences;
  double *correspondences;
  int num_inliers;
  int *inlier_map = NULL;
  int frm_corners[2 * MAX_CORNERS], ref_corners[2 * MAX_CORNERS];
  (void) cpi;


#ifdef USE_FAST_CORNER
  num_frm_corners = FastCornerDetect(frm->y_buffer, frm->y_width,
                                     frm->y_height, frm->y_stride,
                                     frm_corners, MAX_CORNERS);
  num_ref_corners = FastCornerDetect(ref->y_buffer, ref->y_width,
                                     ref->y_height, ref->y_stride,
                                     ref_corners, MAX_CORNERS);
#else
  num_frm_corners = HarrisCornerDetect(frmbuf, width, height, frm_stride,
                                       frm_corners, MAX_CORNERS);
  num_ref_corners = HarrisCornerDetect(refbuf, width, height, ref_stride,
                                       ref_corners, MAX_CORNERS);
#endif
#ifdef VERBOSE
  printf("Frame corners     = %d\n", num_frm_corners);
  printf("Reference corners = %d\n", num_ref_corners);
#endif

  correspondences = (double *) malloc(num_frm_corners * 4 *
                                   sizeof(*correspondences));

  num_correspondences = determine_correspondence(frm->y_buffer,
                                                (int*)frm_corners,
                                                 num_frm_corners,
                                                 ref->y_buffer,
                                                 (int*)ref_corners,
                                                 num_ref_corners,
                                                 frm->y_width, frm->y_height,
                                                 frm->y_stride, ref->y_stride,
                                                 correspondences);
#ifdef VERBOSE
  printf("Number of correspondences = %d\n", num_correspondences);
#endif
  inlier_map = (int *)malloc(num_correspondences * sizeof(int));
  num_inliers = compute_global_motion_single(type, correspondences,
                                             num_correspondences, H,
                                             inlier_map);
#ifdef VERBOSE
  printf("Inliers = %d\n", num_inliers);
  if (num_inliers) {
    printf("Error Score (inliers) = %g\n",
           compute_error_score(type, correspondences, 4, correspondences + 2, 4,
                               num_correspondences, H, inlier_map));
    printf("Warp error score = %g\n",
           vp9_warp_error_unq(ROTZOOM, H, ref->y_buffer, ref->y_crop_width,
                              ref->y_crop_height, ref->y_stride,
                              frm->y_buffer, 0, 0,
                              frm->y_crop_width, frm->y_crop_height,
                              frm->y_stride, 0, 0, 16, 16));
  }
#endif
  free(correspondences);
  free(inlier_map);
  return (num_inliers > 0);
}

static int compute_global_motion_multiple(TransformationType type,
                                          double *correspondences,
                                          int num_correspondences,
                                          double *H,
                                          int max_models,
                                          double inlier_prob,
                                          int *num_models,
                                          int *processed_mask) {
  double *cp = correspondences;
  double *mp, *matched_points;
  int *best_inlier_mask;
  int i, result;
  int num_points = 0;
  int num_inliers = 0;
  int num_inliers_sum = 0;
  const int numparams = get_numparams(type);
  ransacType ransac = get_ransacType(type);
  if (ransac == NULL)
    return 0;
  matched_points =
      (double *)malloc(4 * num_correspondences * sizeof(double));
  best_inlier_mask =
      (int *)malloc(num_correspondences * sizeof(int));
  for (i = 0; i < num_correspondences; ++i)
    processed_mask[i] = -1;
  *num_models = 0;

  while ((double)num_inliers_sum / (double)num_correspondences < inlier_prob &&
         *num_models < max_models) {
    num_points = 0;
    for (mp = matched_points, cp = correspondences, i = 0;
         i < num_correspondences; ++i) {
      if (processed_mask[i] == -1) {
        *mp++ = *cp++;
        *mp++ = *cp++;
        *mp++ = *cp++;
        *mp++ = *cp++;
        num_points++;
      } else {
        cp += 4;
      }
    }
    result = ransac(matched_points, num_points,
                    &num_inliers, best_inlier_mask,
                    &H[(*num_models) * numparams]);
    if (!result && num_inliers < MIN_INLIER_PROB * num_correspondences) {
      result = 1;
      num_inliers = 0;
    }
    if (!result) {
      num_points = 0;
      for (i = 0; i < num_correspondences; ++i) {
        if (processed_mask[i] == -1) {
          if (best_inlier_mask[num_points]) processed_mask[i] = *num_models;
          num_points++;
        }
      }
      num_inliers_sum += num_inliers;
      (*num_models)++;
    } else {
      break;
    }
  }
  free(matched_points);
  free(best_inlier_mask);
  return num_inliers_sum;
}

int vp9_compute_global_motion_multiple_optical_flow(struct VP9_COMP *cpi,
                                                    TransformationType type,
                                                    YV12_BUFFER_CONFIG *frm,
                                                    YV12_BUFFER_CONFIG *ref,
                                                    int max_models,
                                                    double inlier_prob,
                                                    double *H) {
  int num_correspondences = 0;
  double *correspondence_pts = (double *)malloc(frm->y_width * frm->y_height *
                                                sizeof(correspondence));
  correspondence *correspondences = (correspondence *)correspondence_pts;
  int num_inliers;
  int i, j;
  int num_models = 0;
  int *inlier_map = NULL;
  double *u = (double *)malloc(frm->y_width * frm->y_height * sizeof(u));
  double *v = (double *)malloc(frm->y_width * frm->y_height * sizeof(v));
  double *confidence = (double *)malloc(frm->y_width * frm->y_height *
                                        sizeof(confidence));

  (void) cpi;
  compute_flow(frm->y_buffer, ref->y_buffer, u, v, confidence,
               frm->y_width, frm->y_height, frm->y_stride,
               ref->y_stride, 1);

  // avoid the boundaries because of artifacts
  for (i = 10; i < frm->y_width - 10; ++i)
    for (j = 10; j < frm->y_height - 10; ++j) {
      // Note that this threshold could still benefit from tuning. Confidence
      // is based on the inverse of the harris corner score.
      if (confidence[i + j * frm->y_width] < 0.01) {
        correspondences[num_correspondences].x = i;
        correspondences[num_correspondences].y = j;
        correspondences[num_correspondences].rx = i - u[i + j * frm->y_width];
        correspondences[num_correspondences].ry = j - v[i + j * frm->y_width];
        num_correspondences++;
      }
    }
  printf("Number of correspondences = %d\n", num_correspondences);
#ifdef VERBOSE
  printf("Number of correspondences = %d\n", num_correspondences);
#endif
  inlier_map = (int *)malloc(num_correspondences * sizeof(*inlier_map));
  num_inliers = compute_global_motion_multiple(type, correspondence_pts,
                                               num_correspondences, H,
                                               max_models, inlier_prob,
                                               &num_models, inlier_map);
#ifdef VERBOSE
  printf("Models = %d, Inliers = %d\n", num_models, num_inliers);
  if (num_models)
    printf("Error Score (inliers) = %g\n",
           compute_error_score(type, correspondences, 4, correspondences + 2, 4,
                               num_correspondences, H, inlier_map));
    printf("Warp error score = %g\n",
           vp9_warp_error_unq(ROTZOOM, H, ref->y_buffer, ref->y_crop_width,
                              ref->y_crop_height, ref->y_stride,
                              frm->y_buffer, 0, 0,
                              frm->y_crop_width, frm->y_crop_height,
                              frm->y_stride, 0, 0, 16, 16));
#endif
  (void) num_inliers;
  free(correspondences);
  free(inlier_map);
  free(u);
  free(v);
  free(confidence);
  return num_models;
}

// Returns number of models actually returned
int vp9_compute_global_motion_multiple_feature_based(
    struct VP9_COMP *cpi,
    TransformationType type,
    YV12_BUFFER_CONFIG *frm,
    YV12_BUFFER_CONFIG *ref,
    int max_models,
    double inlier_prob,
    double *H) {
  int num_frm_corners, num_ref_corners;
  int num_correspondences;
  double *correspondences;
  int num_inliers;
  int i, j;
  int frm_corners[2 * MAX_CORNERS], ref_corners[2 * MAX_CORNERS];
  int num_models = 0;
  int *inlier_map = NULL;
  (void) cpi;
  (void) i;
  (void) j;

#ifdef USE_FAST_CORNER
  num_frm_corners = FastCornerDetect(frm->y_buffer, frm->y_width,
                                     frm->y_height, frm->y_stride,
                                     frm_corners, MAX_CORNERS);
  num_ref_corners = FastCornerDetect(ref->y_buffer, ref->y_width,
                                     ref->y_height, ref->y_stride,
                                     ref_corners, MAX_CORNERS);
#else
  num_frm_corners = HarrisCornerDetect(frmbuf, width, height, frm_stride,
                                       frm_corners, MAX_CORNERS);
  num_ref_corners = HarrisCornerDetect(refbuf, width, height, ref_stride,
                                       ref_corners, MAX_CORNERS);
#endif
#ifdef VERBOSE
  printf("Frame corners     = %d\n", num_frm_corners);
  printf("Reference corners = %d\n", num_ref_corners);
#endif

  correspondences = (double *) malloc(num_frm_corners * 4 *
                                   sizeof(*correspondences));

  num_correspondences = determine_correspondence(frm->y_buffer,
                                                 (int*)frm_corners,
                                                 num_frm_corners,
                                                 ref->y_buffer,
                                                 (int*)ref_corners,
                                                 num_ref_corners,
                                                 frm->y_width, frm->y_height,
                                                 frm->y_stride, ref->y_stride,
                                                 correspondences);
  printf("Number of correspondences = %d\n", num_correspondences);
#ifdef VERBOSE
  printf("Number of correspondences = %d\n", num_correspondences);
#endif
  inlier_map = (int *)malloc(num_correspondences * sizeof(int));
  num_inliers = compute_global_motion_multiple(type, correspondences,
                                               num_correspondences, H,
                                               max_models, inlier_prob,
                                               &num_models, inlier_map);
#ifdef VERBOSE
  printf("Models = %d, Inliers = %d\n", num_models, num_inliers);
  if (num_models)
    printf("Error Score (inliers) = %g\n",
           compute_error_score(type, correspondences, 4, correspondences + 2, 4,
                               num_correspondences, H, inlier_map));
    printf("Warp error score = %g\n",
           vp9_warp_error_unq(ROTZOOM, H, ref->y_buffer, ref->y_crop_width,
                              ref->y_crop_height, ref->y_stride,
                              frm->y_buffer, 0, 0,
                              frm->y_crop_width, frm->y_crop_height,
                              frm->y_stride, 0, 0, 16, 16));
#endif
  (void) num_inliers;
  free(correspondences);
  free(inlier_map);
  return num_models;
}

// Returns number of models actually returned: 1 - if success, 0 - if failure
int vp9_compute_global_motion_single_block_based(struct VP9_COMP *cpi,
                                                 TransformationType type,
                                                 YV12_BUFFER_CONFIG *frm,
                                                 YV12_BUFFER_CONFIG *ref,
                                                 BLOCK_SIZE bsize,
                                                 double *H) {
  VP9_COMMON *const cm = &cpi->common;
  int num_correspondences = 0;
  double *correspondences;
  int num_inliers;
  int *inlier_map = NULL;
  int bwidth = num_4x4_blocks_wide_lookup[bsize] << 2;
  int bheight = num_4x4_blocks_high_lookup[bsize] << 2;
  int i;
  MV motionfield[4096];
  double confidence[4096];

  vp9_get_frame_motionfield(cpi, frm, ref, bsize, motionfield, confidence);

  correspondences = (double *)malloc(4 * cm->mb_rows * cm->mb_cols *
                                  sizeof(*correspondences));

  for (i = 0; i < cm->mb_rows * cm->mb_cols; i ++) {
      int x = (i % cm->mb_cols) * bwidth + bwidth / 2;
      int y = (i / cm->mb_cols) * bheight + bheight / 2;
      if (confidence[i] > CONFIDENCE_THRESHOLD) {
        correspondences[num_correspondences * 4]   = x;
        correspondences[num_correspondences * 4 + 1] = y;
        correspondences[num_correspondences * 4 + 2] =
            (double)motionfield[i].col / 8 + x;
        correspondences[num_correspondences * 4 + 3] =
            (double)motionfield[i].row / 8 + y;
        num_correspondences++;
      }
  }

  inlier_map = (int *)malloc(num_correspondences * sizeof(*inlier_map));

  num_inliers = compute_global_motion_single(type, correspondences,
                                             num_correspondences, H,
                                             inlier_map);
#ifdef VERBOSE
  printf("Inliers = %d\n", num_inliers);
  printf("Error Score (inliers) = %g\n",
         compute_error_score(type, correspondences, 4, correspondences + 2, 4,
                             num_correspondences, H, inlier_map));
#endif
  free(correspondences);
  free(inlier_map);
  return (num_inliers > 0);
}

// Returns number of models actually returned
int vp9_compute_global_motion_multiple_block_based(struct VP9_COMP *cpi,
                                                   TransformationType type,
                                                   YV12_BUFFER_CONFIG *frm,
                                                   YV12_BUFFER_CONFIG *ref,
                                                   BLOCK_SIZE bsize,
                                                   int max_models,
                                                   double inlier_prob,
                                                   double *H) {
  VP9_COMMON *const cm = &cpi->common;
  int num_correspondences = 0;
  double *correspondences;
  int num_inliers;
  int num_models = 0;
  int *inlier_map = NULL;
  int bwidth = num_4x4_blocks_wide_lookup[bsize] << 2;
  int bheight = num_4x4_blocks_high_lookup[bsize] << 2;

  int i;
  MV motionfield[4096];
  double confidence[4096];
  vp9_get_frame_motionfield(cpi, frm, ref, bsize, motionfield, confidence);

  correspondences = (double *)malloc(4 * cm->mb_rows * cm->mb_cols *
                                  sizeof(*correspondences));

  for (i = 0; i < cm->mb_rows * cm->mb_cols; i ++) {
      int x = (i % cm->mb_cols) * bwidth + bwidth / 2;
      int y = (i / cm->mb_cols) * bheight + bheight / 2;
      if (confidence[i] > CONFIDENCE_THRESHOLD) {
        correspondences[num_correspondences * 4]   = x;
        correspondences[num_correspondences * 4 + 1] = y;
        correspondences[num_correspondences * 4 + 2] =
            (double)motionfield[i].col / 8 + x;
        correspondences[num_correspondences * 4 + 3] =
            (double)motionfield[i].row / 8 + y;
        num_correspondences++;
      }
  }

  inlier_map = (int *)malloc(num_correspondences * sizeof(*inlier_map));
  num_inliers = compute_global_motion_multiple(type, correspondences,
                                               num_correspondences, H,
                                               max_models, inlier_prob,
                                               &num_models, inlier_map);
#ifdef VERBOSE
  printf("Models = %d, Inliers = %d\n", num_models, num_inliers);
  if (num_models)
    printf("Error Score (inliers) = %g\n",
           compute_error_score(type, correspondences, 4, correspondences + 2, 4,
                               num_correspondences, H, inlier_map));
#endif
  (void) num_inliers;
  free(correspondences);
  free(inlier_map);
  return num_models;
}
