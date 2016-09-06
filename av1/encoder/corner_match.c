/*
 *  Copyright (c) 2016 The WebM project authors. All Rights Reserved.
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

#include "av1/encoder/corner_match.h"

#define MATCH_SZ 15
#define MATCH_SZ_BY2 ((MATCH_SZ - 1) / 2)
#define MATCH_SZ_SQ (MATCH_SZ * MATCH_SZ)
#define SEARCH_SZ 9
#define SEARCH_SZ_BY2 ((SEARCH_SZ - 1) / 2)

#define THRESHOLD_NCC 0.80

static double compute_variance(unsigned char *im, int stride, int x, int y,
                               double *mean) {
  double sum = 0.0;
  double sumsq = 0.0;
  double var;
  int i, j;
  for (i = 0; i < MATCH_SZ; ++i)
    for (j = 0; j < MATCH_SZ; ++j) {
      sum += im[(i + y - MATCH_SZ_BY2) * stride + (j + x - MATCH_SZ_BY2)];
      sumsq += im[(i + y - MATCH_SZ_BY2) * stride + (j + x - MATCH_SZ_BY2)] *
               im[(i + y - MATCH_SZ_BY2) * stride + (j + x - MATCH_SZ_BY2)];
    }
  var = (sumsq * MATCH_SZ_SQ - sum * sum) / (MATCH_SZ_SQ * MATCH_SZ_SQ);
  if (mean) *mean = sum / MATCH_SZ_SQ;
  return var;
}

static double compute_cross_correlation(unsigned char *im1, int stride1, int x1,
                                        int y1, unsigned char *im2, int stride2,
                                        int x2, int y2) {
  double sum1 = 0;
  double sum2 = 0;
  double cross = 0;
  double corr;
  int i, j;
  for (i = 0; i < MATCH_SZ; ++i)
    for (j = 0; j < MATCH_SZ; ++j) {
      sum1 += im1[(i + y1 - MATCH_SZ_BY2) * stride1 + (j + x1 - MATCH_SZ_BY2)];
      sum2 += im2[(i + y2 - MATCH_SZ_BY2) * stride2 + (j + x2 - MATCH_SZ_BY2)];
      cross +=
          im1[(i + y1 - MATCH_SZ_BY2) * stride1 + (j + x1 - MATCH_SZ_BY2)] *
          im2[(i + y2 - MATCH_SZ_BY2) * stride2 + (j + x2 - MATCH_SZ_BY2)];
    }
  corr = (cross * MATCH_SZ_SQ - sum1 * sum2) / (MATCH_SZ_SQ * MATCH_SZ_SQ);
  return corr;
}

static int is_eligible_point(double pointx, double pointy, int width,
                             int height) {
  return (pointx >= MATCH_SZ_BY2 && pointy >= MATCH_SZ_BY2 &&
          pointx + MATCH_SZ_BY2 < width && pointy + MATCH_SZ_BY2 < height);
}

static int is_eligible_distance(double point1x, double point1y, double point2x,
                                double point2y, int width, int height) {
  const int thresh = (width < height ? height : width) >> 4;
  return ((point1x - point2x) * (point1x - point2x) +
          (point1y - point2y) * (point1y - point2y)) <= thresh * thresh;
}

static void improve_correspondence(unsigned char *frm, unsigned char *ref,
                                   int width, int height, int frm_stride,
                                   int ref_stride,
                                   Correspondence *correspondences,
                                   int num_correspondences) {
  int i;
  for (i = 0; i < num_correspondences; ++i) {
    double template_norm =
        compute_variance(frm, frm_stride, (int)correspondences[i].x,
                         (int)correspondences[i].y, NULL);
    int x, y, best_x = 0, best_y = 0;
    double best_match_ncc = 0.0;
    for (y = -SEARCH_SZ_BY2; y <= SEARCH_SZ_BY2; ++y) {
      for (x = -SEARCH_SZ_BY2; x <= SEARCH_SZ_BY2; ++x) {
        double match_ncc;
        double subimage_norm;
        if (!is_eligible_point((int)correspondences[i].rx + x,
                               (int)correspondences[i].ry + y, width, height))
          continue;
        if (!is_eligible_distance(
                (int)correspondences[i].x, (int)correspondences[i].y,
                (int)correspondences[i].rx + x, (int)correspondences[i].ry + y,
                width, height))
          continue;
        subimage_norm =
            compute_variance(ref, ref_stride, (int)correspondences[i].rx + x,
                             (int)correspondences[i].ry + y, NULL);
        match_ncc = compute_cross_correlation(
                        frm, frm_stride, (int)correspondences[i].x,
                        (int)correspondences[i].y, ref, ref_stride,
                        (int)correspondences[i].rx + x,
                        (int)correspondences[i].ry + y) /
                    sqrt(template_norm * subimage_norm);
        if (match_ncc > best_match_ncc) {
          best_match_ncc = match_ncc;
          best_y = y;
          best_x = x;
        }
      }
    }
    correspondences[i].rx += (double)best_x;
    correspondences[i].ry += (double)best_y;
  }
  for (i = 0; i < num_correspondences; ++i) {
    double template_norm =
        compute_variance(ref, ref_stride, (int)correspondences[i].rx,
                         (int)correspondences[i].ry, NULL);
    int x, y, best_x = 0, best_y = 0;
    double best_match_ncc = 0.0;
    for (y = -SEARCH_SZ_BY2; y <= SEARCH_SZ_BY2; ++y)
      for (x = -SEARCH_SZ_BY2; x <= SEARCH_SZ_BY2; ++x) {
        double match_ncc;
        double subimage_norm;
        if (!is_eligible_point((int)correspondences[i].x + x,
                               (int)correspondences[i].y + y, width, height))
          continue;
        if (!is_eligible_distance((int)correspondences[i].x + x,
                                  (int)correspondences[i].y + y,
                                  (int)correspondences[i].rx,
                                  (int)correspondences[i].ry, width, height))
          continue;
        subimage_norm =
            compute_variance(frm, frm_stride, (int)correspondences[i].x + x,
                             (int)correspondences[i].y + y, NULL);
        match_ncc =
            compute_cross_correlation(
                frm, frm_stride, (int)correspondences[i].x + x,
                (int)correspondences[i].y + y, ref, ref_stride,
                (int)correspondences[i].rx, (int)correspondences[i].ry) /
            sqrt(template_norm * subimage_norm);
        if (match_ncc > best_match_ncc) {
          best_match_ncc = match_ncc;
          best_y = y;
          best_x = x;
        }
      }
    correspondences[i].x += best_x;
    correspondences[i].y += best_y;
  }
}

int determine_correspondence(unsigned char *frm, int *frm_corners,
                             int num_frm_corners, unsigned char *ref,
                             int *ref_corners, int num_ref_corners, int width,
                             int height, int frm_stride, int ref_stride,
                             double *correspondence_pts) {
  // TODO(sarahparker) Improve this to include 2-way match
  int i, j;
  Correspondence *correspondences = (Correspondence *)correspondence_pts;
  int num_correspondences = 0;
  for (i = 0; i < num_frm_corners; ++i) {
    double best_match_ncc = 0.0;
    double template_norm;
    int best_match_j = -1;
    if (!is_eligible_point(frm_corners[2 * i], frm_corners[2 * i + 1], width,
                           height))
      continue;
    template_norm = compute_variance(frm, frm_stride, frm_corners[2 * i],
                                     frm_corners[2 * i + 1], NULL);
    for (j = 0; j < num_ref_corners; ++j) {
      double match_ncc;
      double subimage_norm;
      if (!is_eligible_point(ref_corners[2 * j], ref_corners[2 * j + 1], width,
                             height))
        continue;
      if (!is_eligible_distance(frm_corners[2 * i], frm_corners[2 * i + 1],
                                ref_corners[2 * j], ref_corners[2 * j + 1],
                                width, height))
        continue;
      subimage_norm = compute_variance(ref, ref_stride, ref_corners[2 * j],
                                       ref_corners[2 * j + 1], NULL);
      match_ncc = compute_cross_correlation(frm, frm_stride, frm_corners[2 * i],
                                            frm_corners[2 * i + 1], ref,
                                            ref_stride, ref_corners[2 * j],
                                            ref_corners[2 * j + 1]) /
                  sqrt(template_norm * subimage_norm);
      if (match_ncc > best_match_ncc) {
        best_match_ncc = match_ncc;
        best_match_j = j;
      }
    }
    if (best_match_ncc > THRESHOLD_NCC) {
      correspondences[num_correspondences].x = (double)frm_corners[2 * i];
      correspondences[num_correspondences].y = (double)frm_corners[2 * i + 1];
      correspondences[num_correspondences].rx =
          (double)ref_corners[2 * best_match_j];
      correspondences[num_correspondences].ry =
          (double)ref_corners[2 * best_match_j + 1];
      num_correspondences++;
    }
  }
  improve_correspondence(frm, ref, width, height, frm_stride, ref_stride,
                         correspondences, num_correspondences);
  return num_correspondences;
}
