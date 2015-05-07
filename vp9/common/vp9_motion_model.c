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

#include "vp9/common/vp9_common_data.h"
#include "vp9/common/vp9_mv.h"
#include "vp9/common/vp9_motion_model.h"

INLINE projectPointsType get_projectPointsType(TransformationType type) {
  switch (type) {
    case HOMOGRAPHY:
      return projectPointsHomography;
    case AFFINE:
      return projectPointsAffine;
    case ROTZOOM:
      return projectPointsRotZoom;
    case TRANSLATION:
      return projectPointsTranslation;
    default:
      assert(0);
      return NULL;
  }
}

void projectPointsTranslation(double *mat, double *points, double *proj,
                              const int n,
                              const int stride_points,
                              const int stride_proj) {
  int i;
  for (i = 0; i < n; ++i) {
    const double x = *(points++), y = *(points++);
    *(proj++) = x + mat[0];
    *(proj++) = y + mat[1];
    points += stride_points - 2;
    proj += stride_proj - 2;
  }
}

void projectPointsRotZoom(double *mat, double *points,
                          double *proj, const int n,
                          const int stride_points, const int stride_proj) {
  int i;
  for (i = 0; i < n; ++i) {
    const double x = *(points++), y = *(points++);
    *(proj++) =  mat[0] * x + mat[1] * y + mat[2];
    *(proj++) = -mat[1] * x + mat[0] * y + mat[3];
    points += stride_points - 2;
    proj += stride_proj - 2;
  }
}

void projectPointsAffine(double *mat, double *points,
                         double *proj, const int n,
                         const int stride_points, const int stride_proj) {
  int i;
  for (i = 0; i < n; ++i) {
    const double x = *(points++), y = *(points++);
    *(proj++) = mat[0] * x + mat[1] * y + mat[4];
    *(proj++) = mat[2] * x + mat[3] * y + mat[5];
    points += stride_points - 2;
    proj += stride_proj - 2;
  }
}

void projectPointsHomography(double *mat, double *points,
                             double *proj, const int n,
                             const int stride_points, const int stride_proj) {
  int i;
  double x, y, Z;
  for (i = 0; i < n; ++i) {
    x = *(points++), y = *(points++);
    Z = 1. / (mat[6] * x + mat[7] * y + mat[8]);
    *(proj++) = (mat[0] * x + mat[1] * y + mat[2]) * Z;
    *(proj++) = (mat[3] * x + mat[4] * y + mat[5]) * Z;
    points += stride_points - 2;
    proj += stride_proj - 2;
  }
}

#define clip_pixel(v) ((v) < 0 ? 0 : ((v) > 255 ? 255 : (v)))

static unsigned char bilinear(unsigned char *ref, double x, double y,
                              int width, int height, int stride) {
  if (x < 0 && y < 0) return ref[0];
  else if (x < 0 && y > height - 1)
    return ref[(height - 1) * stride];
  else if (x > width - 1 && y < 0)
    return ref[width - 1];
  else if (x > width - 1 && y > height - 1)
    return ref[(height - 1) * stride + (width - 1)];
  else if (x < 0) {
    int i = (int) y;
    double a = y - i;
    int v = (int)(ref[i * stride] * (1 - a) + ref[(i + 1) * stride] * a + 0.5);
    return clip_pixel(v);
  } else if (y < 0) {
    int j = (int) x;
    double b = x - j;
    int v = (int)(ref[j] * (1 - b) + ref[j + 1] * b + 0.5);
    return clip_pixel(v);
  } else if (x > width - 1) {
    int i = (int) y;
    double a = y - i;
    int v = (int)(ref[i * stride + width - 1] * (1 - a) +
                  ref[(i + 1) * stride + width - 1] * a + 0.5);
    return clip_pixel(v);
  } else if (y > height - 1) {
    int j = (int) x;
    double b = x - j;
    int v = (int)(ref[(height - 1) * stride + j] * (1 - b) +
                  ref[(height - 1) * stride + j + 1] * b + 0.5);
    return clip_pixel(v);
  } else {
    int i = (int) y;
    int j = (int) x;
    double a = y - i;
    double b = x - j;
    int v = (int)(ref[i * stride + j] * (1 - a) * (1 - b) +
                  ref[i * stride + j + 1] * (1 - a) * b +
                  ref[(i + 1) * stride + j] * a * (1 - b) +
                  ref[(i + 1) * stride + j + 1] * a * b);
    return clip_pixel(v);
  }
}

static void WarpImage(TransformationType type, double *H,
                      unsigned char *ref,
                      int width, int height, int stride,
                      unsigned char *pred,
                      int p_col, int p_row,
                      int p_width, int p_height, int p_stride,
                      int subsampling_col, int subsampling_row,
                      int x_scale, int y_scale) {
  int i, j;
  projectPointsType projectPoints = get_projectPointsType(type);
  if (projectPoints == NULL)
    return;
  for (i = p_row; i < p_row + p_height; ++i) {
    for (j = p_col; j < p_col + p_width; ++j) {
      double in[2], out[2];
      in[0] = subsampling_col ? 2 * j + 0.5 : j;
      in[1] = subsampling_row ? 2 * i + 0.5 : i;
      projectPoints(H, in, out, 1, 2, 2);
      out[0] = subsampling_col ? (out[0] - 0.5) / 2.0 : out[0];
      out[1] = subsampling_row ? (out[1] - 0.5) / 2.0 : out[1];
      out[0] *= x_scale / 16.0;
      out[1] *= y_scale / 16.0;
      pred[(j - p_col) + (i - p_row) * p_stride] =
          bilinear(ref, out[0], out[1], width, height, stride);
    }
  }
}

static void convert_params_to_rotzoom(Global_Motion_Params *model,
                                      double *H) {
  double z = 1.0 + (double) model->zoom / (1 << ZOOM_PRECISION_BITS);
  double r = (double) model->rotation / (1 << ROTATION_PRECISION_BITS);
  H[0] =  (1 + z) * cos(r * M_PI / 180.0);
  H[1] = -(1 + z) * sin(r * M_PI / 180.0);
  H[2] = (double) model->mv.as_mv.col / 8.0;
  H[3] = (double) model->mv.as_mv.row / 8.0;
}

void vp9_warp_plane(Global_Motion_Params *gm,
                    unsigned char *ref,
                    int width, int height, int stride,
                    unsigned char *pred,
                    int p_col, int p_row,
                    int p_width, int p_height, int p_stride,
                    int subsampling_col, int subsampling_row,
                    int x_scale, int y_scale) {
  double H[9];
  convert_params_to_rotzoom(gm, H);
  WarpImage(ROTZOOM, H,
            ref, width, height, stride,
            pred, p_col, p_row, p_width, p_height, p_stride,
            subsampling_col,  subsampling_row,
            x_scale, y_scale);
}

static int_mv vp9_get_global_mv(int col, int row, Global_Motion_Params *model) {
  int_mv mv;
  double H[4];
  double x, y;
  convert_params_to_rotzoom(model, H);
  x =  H[0] * col + H[1] * row + H[2];
  y = -H[1] * col + H[0] * row + H[3];
  mv.as_mv.col = (int)floor(x * 8 + 0.5) - col;
  mv.as_mv.row = (int)floor(y * 8 + 0.5) - row;
  return mv;
}

int_mv vp9_get_global_sb_center_mv(int col, int row, int bw, int bh,
                                   Global_Motion_Params *model) {
  col += bw / 2;
  row += bh / 2;
  return vp9_get_global_mv(col, row, model);
}

int_mv vp9_get_global_sub8x8_center_mv(int col, int row, int block,
                                       Global_Motion_Params *model) {
  if (block == 0 || block == 2)
    col += 2;
  else
    col += 6;
  if (block == 0 || block == 1)
    row += 2;
  else
    row += 6;
  return vp9_get_global_mv(col, row, model);
}
