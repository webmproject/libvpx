/*
 *  Copyright (c) 2015 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be
 *  found  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#ifndef VP9_COMMON_VP9_MOTION_MODEL_H
#define VP9_COMMON_VP9_MOTION_MODEL_H

#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <math.h>
#include <assert.h>

#include "./vpx_config.h"
#include "vpx_ports/mem.h"
#include "vp9/common/vp9_enums.h"
#include "vp9/common/vp9_mv.h"

typedef void (*projectPointsType)(double *mat, double *points, double *proj,
                                  const int n, const int stride_points,
                                  const int stride_proj);
typedef enum {
  UNKNOWN_TRANSFORM = -1,
  HOMOGRAPHY,  // homography, 8-parameter
  AFFINE,      // affine, 6-parameter
  ROTZOOM,     // simplified affine with rotation and zoom only, 4-parameter
  TRANSLATION  // translational motion 2-parameter
} TransformationType;

static INLINE int get_numparams(TransformationType type) {
  switch (type) {
    case HOMOGRAPHY:
      return 9;
    case AFFINE:
      return 6;
    case ROTZOOM:
      return 4;
    case TRANSLATION:
      return 2;
    default:
      assert(0);
      return 0;
  }
}

void projectPointsHomography(double *mat, double *points, double *proj,
                             const int n, const int stride_points,
                             const int stride_proj);
void projectPointsAffine(double *mat, double *points, double *proj,
                         const int n, const int stride_points,
                         const int stride_proj);
void projectPointsRotZoom(double *mat, double *points, double *proj,
                          const int n, const int stride_points,
                          const int stride_proj);
void projectPointsTranslation(double *mat, double *points, double *proj,
                              const int n, const int stride_points,
                              const int stride_proj);

projectPointsType get_projectPointsType(TransformationType type);

void vp9_convert_params_to_rotzoom(Global_Motion_Params *model, double *H);

void vp9_warp_plane(Global_Motion_Params *gm,
                    unsigned char *ref,
                    int width, int height, int stride,
                    unsigned char *pred,
                    int p_col, int p_row,
                    int p_width, int p_height, int p_stride,
                    int subsampling_col, int subsampling_row,
                    int x_scale, int y_scale);

double vp9_warp_erroradv(Global_Motion_Params *gm,
                         unsigned char *ref,
                         int width, int height, int stride,
                         unsigned char *src,
                         int p_col, int p_row,
                         int p_width, int p_height, int p_stride,
                         int subsampling_col, int subsampling_row,
                         int x_scale, int y_scale);
double vp9_warp_erroradv_unq(TransformationType type, double *H,
                             unsigned char *ref,
                             int width, int height, int stride,
                             unsigned char *src,
                             int p_col, int p_row,
                             int p_width, int p_height, int p_stride,
                             int subsampling_col, int subsampling_row,
                             int x_scale, int y_scale);

double compute_warp_and_error(Global_Motion_Params *gm,
                              projectPointsType projectPoints,
                              unsigned char *ref,
                              int width, int height, int stride,
                              unsigned char *src,
                              int p_col, int p_row,
                              int p_width, int p_height, int p_stride,
                              int subsampling_col, int subsampling_row,
                              int x_scale, int y_scale);

unsigned char interpolate(unsigned char *ref, double x, double y,
                          int width, int height, int stride);


int_mv vp9_get_global_sb_center_mv(int col, int row, int bw, int bh,
                                   Global_Motion_Params *model);
int_mv vp9_get_global_sub8x8_center_mv(int col, int row, int block,
                                       Global_Motion_Params *model);

#endif  // VP9_COMMON_VP9_MOTION_MODEL_H
