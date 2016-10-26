/*
 *   (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include <memory.h>
#include <math.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "av1/encoder/ransac.h"

#define MAX_MINPTS 4

#define MAX_DEGENERATE_ITER 10
#define MINPTS_MULTIPLIER 5

////////////////////////////////////////////////////////////////////////////////
// ransac
typedef int (*IsDegenerateFunc)(double *p);
typedef void (*NormalizeFunc)(double *p, int np, double *T);
typedef void (*DenormalizeFunc)(double *params, double *T1, double *T2);
typedef int (*FindTransformationFunc)(int points, double *points1,
                                      double *points2, double *params);
typedef void (*ProjectPointsDoubleFunc)(double *mat, double *points,
                                        double *proj, const int n,
                                        const int stride_points,
                                        const int stride_proj);

static void project_points_double_translation(double *mat, double *points,
                                              double *proj, const int n,
                                              const int stride_points,
                                              const int stride_proj) {
  int i;
  for (i = 0; i < n; ++i) {
    const double x = *(points++), y = *(points++);
    *(proj++) = x + mat[1];
    *(proj++) = y + mat[0];
    points += stride_points - 2;
    proj += stride_proj - 2;
  }
}

static void project_points_double_rotzoom(double *mat, double *points,
                                          double *proj, const int n,
                                          const int stride_points,
                                          const int stride_proj) {
  int i;
  for (i = 0; i < n; ++i) {
    const double x = *(points++), y = *(points++);
    *(proj++) = mat[3] * x + mat[2] * y + mat[1];
    *(proj++) = -mat[2] * x + mat[3] * y + mat[0];
    points += stride_points - 2;
    proj += stride_proj - 2;
  }
}

static void project_points_double_affine(double *mat, double *points,
                                         double *proj, const int n,
                                         const int stride_points,
                                         const int stride_proj) {
  int i;
  for (i = 0; i < n; ++i) {
    const double x = *(points++), y = *(points++);
    *(proj++) = mat[3] * x + mat[2] * y + mat[1];
    *(proj++) = mat[4] * x + mat[5] * y + mat[0];
    points += stride_points - 2;
    proj += stride_proj - 2;
  }
}

static void project_points_double_homography(double *mat, double *points,
                                             double *proj, const int n,
                                             const int stride_points,
                                             const int stride_proj) {
  int i;
  double x, y, Z, Z_inv;
  for (i = 0; i < n; ++i) {
    x = *(points++), y = *(points++);
    Z_inv = mat[7] * x + mat[6] * y + 1;
    assert(fabs(Z_inv) > 0.00001);
    Z = 1. / Z_inv;
    *(proj++) = (mat[1] * x + mat[0] * y + mat[3]) * Z;
    *(proj++) = (mat[2] * x + mat[4] * y + mat[4]) * Z;
    points += stride_points - 2;
    proj += stride_proj - 2;
  }
}

static int get_rand_indices(int npoints, int minpts, int *indices,
                            unsigned int *seed) {
  int i, j;
  int ptr = rand_r(seed) % npoints;
  if (minpts > npoints) return 0;
  indices[0] = ptr;
  ptr = (ptr == npoints - 1 ? 0 : ptr + 1);
  i = 1;
  while (i < minpts) {
    int index = rand_r(seed) % npoints;
    while (index) {
      ptr = (ptr == npoints - 1 ? 0 : ptr + 1);
      for (j = 0; j < i; ++j) {
        if (indices[j] == ptr) break;
      }
      if (j == i) index--;
    }
    indices[i++] = ptr;
  }
  return 1;
}

static int ransac(double *matched_points, int npoints, int *number_of_inliers,
                  int *best_inlier_mask, double *best_params, const int minpts,
                  const int paramdim, IsDegenerateFunc is_degenerate,
                  NormalizeFunc normalize, DenormalizeFunc denormalize,
                  FindTransformationFunc find_transformation,
                  ProjectPointsDoubleFunc projectpoints) {
  static const double INLIER_THRESHOLD_NORMALIZED = 0.1;
  static const double INLIER_THRESHOLD_UNNORMALIZED = 1.0;
  static const double PROBABILITY_REQUIRED = 0.9;
  static const double EPS = 1e-12;
  static const int MIN_TRIALS = 20;

  const double inlier_threshold =
      (normalize && denormalize ? INLIER_THRESHOLD_NORMALIZED
                                : INLIER_THRESHOLD_UNNORMALIZED);
  int N = 10000, trial_count = 0;
  int i;
  int ret_val = 0;
  unsigned int seed = (unsigned int)npoints;

  int max_inliers = 0;
  double best_variance = 0.0;
  double params[MAX_PARAMDIM];
  WarpedMotionParams wm;
  double points1[2 * MAX_MINPTS];
  double points2[2 * MAX_MINPTS];
  int indices[MAX_MINPTS] = { 0 };

  double *best_inlier_set1;
  double *best_inlier_set2;
  double *inlier_set1;
  double *inlier_set2;
  double *corners1;
  double *corners2;
  double *image1_coord;
  int *inlier_mask;

  double *cnp1, *cnp2;
  double T1[9], T2[9];

  *number_of_inliers = 0;
  if (npoints < minpts * MINPTS_MULTIPLIER || npoints == 0) {
    printf("Cannot find motion with %d matches\n", npoints);
    return 1;
  }

  memset(&wm, 0, sizeof(wm));
  best_inlier_set1 =
      (double *)aom_malloc(sizeof(*best_inlier_set1) * npoints * 2);
  best_inlier_set2 =
      (double *)aom_malloc(sizeof(*best_inlier_set2) * npoints * 2);
  inlier_set1 = (double *)aom_malloc(sizeof(*inlier_set1) * npoints * 2);
  inlier_set2 = (double *)aom_malloc(sizeof(*inlier_set2) * npoints * 2);
  corners1 = (double *)aom_malloc(sizeof(*corners1) * npoints * 2);
  corners2 = (double *)aom_malloc(sizeof(*corners2) * npoints * 2);
  image1_coord = (double *)aom_malloc(sizeof(*image1_coord) * npoints * 2);
  inlier_mask = (int *)aom_malloc(sizeof(*inlier_mask) * npoints);

  if (!(best_inlier_set1 && best_inlier_set2 && inlier_set1 && inlier_set2 &&
        corners1 && corners2 && image1_coord && inlier_mask)) {
    ret_val = 1;
    goto finish_ransac;
  }

  for (cnp1 = corners1, cnp2 = corners2, i = 0; i < npoints; ++i) {
    *(cnp1++) = *(matched_points++);
    *(cnp1++) = *(matched_points++);
    *(cnp2++) = *(matched_points++);
    *(cnp2++) = *(matched_points++);
  }
  matched_points -= 4 * npoints;

  if (normalize && denormalize) {
    normalize(corners1, npoints, T1);
    normalize(corners2, npoints, T2);
  }

  while (N > trial_count) {
    int num_inliers = 0;
    double sum_distance = 0.0;
    double sum_distance_squared = 0.0;

    int degenerate = 1;
    int num_degenerate_iter = 0;
    while (degenerate) {
      num_degenerate_iter++;
      if (!get_rand_indices(npoints, minpts, indices, &seed)) {
        ret_val = 1;
        goto finish_ransac;
      }
      i = 0;
      while (i < minpts) {
        int index = indices[i];
        // add to list
        points1[i * 2] = corners1[index * 2];
        points1[i * 2 + 1] = corners1[index * 2 + 1];
        points2[i * 2] = corners2[index * 2];
        points2[i * 2 + 1] = corners2[index * 2 + 1];
        i++;
      }
      degenerate = is_degenerate(points1);
      if (num_degenerate_iter > MAX_DEGENERATE_ITER) {
        ret_val = 1;
        goto finish_ransac;
      }
    }

    if (find_transformation(minpts, points1, points2, params)) {
      trial_count++;
      continue;
    }

    projectpoints(params, corners1, image1_coord, npoints, 2, 2);

    for (i = 0; i < npoints; ++i) {
      double dx = image1_coord[i * 2] - corners2[i * 2];
      double dy = image1_coord[i * 2 + 1] - corners2[i * 2 + 1];
      double distance = sqrt(dx * dx + dy * dy);

      inlier_mask[i] = distance < inlier_threshold;
      if (inlier_mask[i]) {
        inlier_set1[num_inliers * 2] = corners1[i * 2];
        inlier_set1[num_inliers * 2 + 1] = corners1[i * 2 + 1];
        inlier_set2[num_inliers * 2] = corners2[i * 2];
        inlier_set2[num_inliers * 2 + 1] = corners2[i * 2 + 1];
        num_inliers++;
        sum_distance += distance;
        sum_distance_squared += distance * distance;
      }
    }

    if (num_inliers >= max_inliers && num_inliers > 1) {
      int temp;
      double fracinliers, pNoOutliers, mean_distance, variance;

      assert(num_inliers > 1);
      mean_distance = sum_distance / ((double)num_inliers);
      variance = sum_distance_squared / ((double)num_inliers - 1.0) -
                 mean_distance * mean_distance * ((double)num_inliers) /
                     ((double)num_inliers - 1.0);
      if ((num_inliers > max_inliers) ||
          (num_inliers == max_inliers && variance < best_variance)) {
        best_variance = variance;
        max_inliers = num_inliers;
        memcpy(best_params, params, paramdim * sizeof(*best_params));
        memcpy(best_inlier_set1, inlier_set1,
               num_inliers * 2 * sizeof(*best_inlier_set1));
        memcpy(best_inlier_set2, inlier_set2,
               num_inliers * 2 * sizeof(*best_inlier_set2));
        memcpy(best_inlier_mask, inlier_mask,
               npoints * sizeof(*best_inlier_mask));

        assert(npoints > 0);
        fracinliers = (double)num_inliers / (double)npoints;
        pNoOutliers = 1 - pow(fracinliers, minpts);
        pNoOutliers = fmax(EPS, pNoOutliers);
        pNoOutliers = fmin(1 - EPS, pNoOutliers);
        assert(fabs(1.0 - pNoOutliers) > 0.00001);
        temp = (int)(log(1.0 - PROBABILITY_REQUIRED) / log(pNoOutliers));
        if (temp > 0 && temp < N) {
          N = AOMMAX(temp, MIN_TRIALS);
        }
      }
    }
    trial_count++;
  }
  find_transformation(max_inliers, best_inlier_set1, best_inlier_set2,
                      best_params);
  if (normalize && denormalize) {
    denormalize(best_params, T1, T2);
  }
  *number_of_inliers = max_inliers;
finish_ransac:
  aom_free(best_inlier_set1);
  aom_free(best_inlier_set2);
  aom_free(inlier_set1);
  aom_free(inlier_set2);
  aom_free(corners1);
  aom_free(corners2);
  aom_free(image1_coord);
  aom_free(inlier_mask);
  return ret_val;
}

static int is_collinear3(double *p1, double *p2, double *p3) {
  static const double collinear_eps = 1e-3;
  const double v =
      (p2[0] - p1[0]) * (p3[1] - p1[1]) - (p2[1] - p1[1]) * (p3[0] - p1[0]);
  return fabs(v) < collinear_eps;
}

static int is_degenerate_translation(double *p) {
  return (p[0] - p[2]) * (p[0] - p[2]) + (p[1] - p[3]) * (p[1] - p[3]) <= 2;
}

static int is_degenerate_affine(double *p) {
  return is_collinear3(p, p + 2, p + 4);
}

static int is_degenerate_homography(double *p) {
  return is_collinear3(p, p + 2, p + 4) || is_collinear3(p, p + 2, p + 6) ||
         is_collinear3(p, p + 4, p + 6) || is_collinear3(p + 2, p + 4, p + 6);
}

int ransac_translation(double *matched_points, int npoints,
                       int *number_of_inliers, int *best_inlier_mask,
                       double *best_params) {
  return ransac(matched_points, npoints, number_of_inliers, best_inlier_mask,
                best_params, 3, 2, is_degenerate_translation,
                NULL,  // normalize_homography,
                NULL,  // denormalize_rotzoom,
                find_translation, project_points_double_translation);
}

int ransac_rotzoom(double *matched_points, int npoints, int *number_of_inliers,
                   int *best_inlier_mask, double *best_params) {
  return ransac(matched_points, npoints, number_of_inliers, best_inlier_mask,
                best_params, 3, 4, is_degenerate_affine,
                NULL,  // normalize_homography,
                NULL,  // denormalize_rotzoom,
                find_rotzoom, project_points_double_rotzoom);
}

int ransac_affine(double *matched_points, int npoints, int *number_of_inliers,
                  int *best_inlier_mask, double *best_params) {
  return ransac(matched_points, npoints, number_of_inliers, best_inlier_mask,
                best_params, 3, 6, is_degenerate_affine,
                NULL,  // normalize_homography,
                NULL,  // denormalize_affine,
                find_affine, project_points_double_affine);
}

int ransac_homography(double *matched_points, int npoints,
                      int *number_of_inliers, int *best_inlier_mask,
                      double *best_params) {
  const int result =
      ransac(matched_points, npoints, number_of_inliers, best_inlier_mask,
             best_params, 4, 8, is_degenerate_homography,
             NULL,  // normalize_homography,
             NULL,  // denormalize_homography,
             find_homography, project_points_double_homography);
  if (!result) {
    // normalize so that H33 = 1
    int i;
    const double m = 1.0 / best_params[8];
    assert(fabs(best_params[8]) > 0.00001);
    for (i = 0; i < 8; ++i) best_params[i] *= m;
    best_params[8] = 1.0;
  }
  return result;
}
