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

#include "vp10/common/warped_motion.h"


typedef void (*projectPointsType)(int *mat,
                                  int *points,
                                  int *proj,
                                  const int n,
                                  const int stride_points,
                                  const int stride_proj,
                                  const int subsampling_x,
                                  const int subsampling_y);

static void projectPointsHomography(int *mat,
                                    int *points,
                                    int *proj,
                                    const int n,
                                    const int stride_points,
                                    const int stride_proj,
                                    const int subsampling_x,
                                    const int subsampling_y);
static void projectPointsAffine(int *mat,
                                int *points,
                                int *proj,
                                const int n,
                                const int stride_points,
                                const int stride_proj,
                                const int subsampling_x,
                                const int subsampling_y);
static void projectPointsRotZoom(int *mat,
                                 int *points,
                                 int *proj,
                                 const int n,
                                 const int stride_points,
                                 const int stride_proj,
                                 const int subsampling_x,
                                 const int subsampling_y);
static void projectPointsTranslation(int *mat,
                                     int *points,
                                     int *proj,
                                     const int n,
                                     const int stride_points,
                                     const int stride_proj,
                                     const int subsampling_x,
                                     const int subsampling_y);

static projectPointsType get_projectPointsType(TransformationType type) {
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

static void projectPointsTranslation(int *mat, int *points, int *proj,
                                     const int n,
                                     const int stride_points,
                                     const int stride_proj,
                                     const int subsampling_x,
                                     const int subsampling_y) {
  int i;
  for (i = 0; i < n; ++i) {
    const int x = *(points++), y = *(points++);
    if (subsampling_x)
      *(proj++) = ROUND_POWER_OF_TWO_SIGNED(
          ((x << (WARPEDMODEL_PREC_BITS + 1)) + mat[0]),
          WARPEDPIXEL_PREC_BITS + 1);
    else
      *(proj++) = ROUND_POWER_OF_TWO_SIGNED(
          ((x << WARPEDMODEL_PREC_BITS)) + mat[0],
          WARPEDPIXEL_PREC_BITS);
    if (subsampling_y)
      *(proj++) = ROUND_POWER_OF_TWO_SIGNED(
          ((y << (WARPEDMODEL_PREC_BITS + 1)) + mat[1]),
          WARPEDPIXEL_PREC_BITS + 1);
    else
      *(proj++) = ROUND_POWER_OF_TWO_SIGNED(
          ((y << WARPEDMODEL_PREC_BITS)) + mat[1],
          WARPEDPIXEL_PREC_BITS);
    points += stride_points - 2;
    proj += stride_proj - 2;
  }
}

void projectPointsRotZoom(int *mat, int *points, int *proj,
                          const int n,
                          const int stride_points,
                          const int stride_proj,
                          const int subsampling_x,
                          const int subsampling_y) {
  int i;
  for (i = 0; i < n; ++i) {
    const int x = *(points++), y = *(points++);
    if (subsampling_x)
      *(proj++) = ROUND_POWER_OF_TWO_SIGNED(
          mat[2] * 2 * x + mat[3] * 2 * y + mat[0] +
          (mat[2] + mat[3] - (1 << WARPEDMODEL_PREC_BITS)) / 2,
          WARPEDDIFF_PREC_BITS + 1);
    else
      *(proj++) = ROUND_POWER_OF_TWO_SIGNED(mat[2] * x + mat[3] * y + mat[0],
                                            WARPEDDIFF_PREC_BITS);
    if (subsampling_y)
      *(proj++) = ROUND_POWER_OF_TWO_SIGNED(
          -mat[3] * 2 * x + mat[2] * 2 * y + mat[1] +
          (-mat[3] + mat[2] - (1 << WARPEDMODEL_PREC_BITS)) / 2,
          WARPEDDIFF_PREC_BITS + 1);
    else
      *(proj++) = ROUND_POWER_OF_TWO_SIGNED(-mat[3] * x + mat[2] * y + mat[1],
                                            WARPEDDIFF_PREC_BITS);
    points += stride_points - 2;
    proj += stride_proj - 2;
  }
}

static void projectPointsAffine(int *mat, int *points, int *proj,
                                const int n,
                                const int stride_points,
                                const int stride_proj,
                                const int subsampling_x,
                                const int subsampling_y) {
  int i;
  for (i = 0; i < n; ++i) {
    const int x = *(points++), y = *(points++);
    if (subsampling_x)
      *(proj++) = ROUND_POWER_OF_TWO_SIGNED(
          mat[2] * 2 * x + mat[3] * 2 * y + mat[0] +
          (mat[2] + mat[3] - (1 << WARPEDMODEL_PREC_BITS)) / 2,
          WARPEDDIFF_PREC_BITS + 1);
    else
      *(proj++) = ROUND_POWER_OF_TWO_SIGNED(mat[2] * x + mat[3] * y + mat[0],
                                            WARPEDDIFF_PREC_BITS);
    if (subsampling_y)
      *(proj++) = ROUND_POWER_OF_TWO_SIGNED(
          mat[4] * 2 * x + mat[5] * 2 * y + mat[1] +
          (mat[4] + mat[5] - (1 << WARPEDMODEL_PREC_BITS)) / 2,
          WARPEDDIFF_PREC_BITS + 1);
    else
      *(proj++) = ROUND_POWER_OF_TWO_SIGNED(mat[4] * x + mat[5] * y + mat[1],
                                            WARPEDDIFF_PREC_BITS);
    points += stride_points - 2;
    proj += stride_proj - 2;
  }
}

static void projectPointsHomography(int *mat, int *points, int *proj,
                                    const int n,
                                    const int stride_points,
                                    const int stride_proj,
                                    const int subsampling_x,
                                    const int subsampling_y) {
  int i;
  int64_t x, y, Z;
  int64_t xp, yp;
  for (i = 0; i < n; ++i) {
    x = *(points++), y = *(points++);
    x = (subsampling_x ? 4 * x + 1 : 2 * x);
    y = (subsampling_y ? 4 * y + 1 : 2 * y);

    Z = (mat[6] * x + mat[7] * y + (1 << (WARPEDMODEL_ROW3HOMO_PREC_BITS + 1)));
    xp = (mat[0] * x + mat[1] * y + 2 * mat[2])
        << (WARPEDPIXEL_PREC_BITS +
        WARPEDMODEL_ROW3HOMO_PREC_BITS - WARPEDMODEL_PREC_BITS);
    yp = (mat[3] * x + mat[4] * y + 2 * mat[5])
        << (WARPEDPIXEL_PREC_BITS +
        WARPEDMODEL_ROW3HOMO_PREC_BITS - WARPEDMODEL_PREC_BITS);

    xp = xp > 0 ? (xp + Z / 2) / Z : (xp - Z / 2) / Z;
    yp = yp > 0 ? (yp + Z / 2) / Z : (yp - Z / 2) / Z;

    if (subsampling_x)
      xp = (xp - (1 << (WARPEDPIXEL_PREC_BITS - 1))) / 2;
    if (subsampling_y)
      yp = (yp - (1 << (WARPEDPIXEL_PREC_BITS - 1))) / 2;
    *(proj++) = xp;
    *(proj++) = yp;

    points += stride_points - 2;
    proj += stride_proj - 2;
  }
}

static const int16_t
filter_4tap[WARPEDPIXEL_PREC_SHIFTS][4] = {
  {0, 128,   0, 0},
  {-1, 127,   2, 0},
  {-2, 127,   4, -1},
  {-3, 126,   6, -1},
  {-3, 125,   8, -2},
  {-4, 124,  11, -3},
  {-5, 123,  13, -3},
  {-5, 121,  15, -3},
  {-6, 120,  18, -4},
  {-7, 119,  20, -4},
  {-7, 118,  22, -5},
  {-8, 116,  25, -5},
  {-8, 115,  27, -6},
  {-9, 113,  30, -6},
  {-9, 112,  32, -7},
  {-9, 110,  34, -7},
  {-10, 108,  37, -7},
  {-10, 107,  39, -8},
  {-10, 105,  41, -8},
  {-11, 103,  44, -8},
  {-11, 101,  47, -9},
  {-11,  99,  49, -9},
  {-11,  97,  51, -9},
  {-11,  95,  54, -10},
  {-11,  93,  56, -10},
  {-12,  91,  59, -10},
  {-12,  89,  61, -10},
  {-12,  87,  64, -11},
  {-12,  85,  66, -11},
  {-12,  82,  69, -11},
  {-12,  80,  71, -11},
  {-12,  78,  73, -11},
  {-11,  75,  75, -11},
  {-11,  73,  78, -12},
  {-11,  71,  80, -12},
  {-11,  69,  82, -12},
  {-11,  66,  85, -12},
  {-11,  64,  87, -12},
  {-10,  61,  89, -12},
  {-10,  59,  91, -12},
  {-10,  56,  93, -11},
  {-10,  54,  95, -11},
  {-9, 51, 97, -11},
  {-9, 49, 99, -11},
  {-9,  47, 101, -11},
  {-8,  44, 103, -11},
  {-8,  41, 105, -10},
  {-8,  39, 107, -10},
  {-7,  37, 108, -10},
  {-7,  34, 110, -9},
  {-7,  32, 112, -9},
  {-6,  30, 113, -9},
  {-6,  27, 115, -8},
  {-5,  25, 116, -8},
  {-5,  22, 118, -7},
  {-4,  20, 119, -7},
  {-4,  18, 120, -6},
  {-3,  15, 121, -5},
  {-3,  13, 123, -5},
  {-3,  11, 124, -4},
  {-2,   8, 125, -3},
  {-1,   6, 126, -3},
  {-1,   4, 127, -2},
  {0,   2, 127, -1},
};

static const int16_t
filter_ntap[WARPEDPIXEL_PREC_SHIFTS][WARPEDPIXEL_FILTER_TAPS] = {
  {0,   0, 128,   0,   0, 0},
  {0,  -1, 128,   2,  -1, 0},
  {1,  -3, 127,   4,  -1, 0},
  {1,  -4, 126,   6,  -2, 1},
  {1,  -5, 126,   8,  -3, 1},
  {1,  -6, 125,  11,  -4, 1},
  {1,  -7, 124,  13,  -4, 1},
  {2,  -8, 123,  15,  -5, 1},
  {2,  -9, 122,  18,  -6, 1},
  {2, -10, 121,  20,  -6, 1},
  {2, -11, 120,  22,  -7, 2},
  {2, -12, 119,  25,  -8, 2},
  {3, -13, 117,  27,  -8, 2},
  {3, -13, 116,  29,  -9, 2},
  {3, -14, 114,  32, -10, 3},
  {3, -15, 113,  35, -10, 2},
  {3, -15, 111,  37, -11, 3},
  {3, -16, 109,  40, -11, 3},
  {3, -16, 108,  42, -12, 3},
  {4, -17, 106,  45, -13, 3},
  {4, -17, 104,  47, -13, 3},
  {4, -17, 102,  50, -14, 3},
  {4, -17, 100,  52, -14, 3},
  {4, -18,  98,  55, -15, 4},
  {4, -18,  96,  58, -15, 3},
  {4, -18,  94,  60, -16, 4},
  {4, -18,  91,  63, -16, 4},
  {4, -18,  89,  65, -16, 4},
  {4, -18,  87,  68, -17, 4},
  {4, -18,  85,  70, -17, 4},
  {4, -18,  82,  73, -17, 4},
  {4, -18,  80,  75, -17, 4},
  {4, -18,  78,  78, -18, 4},
  {4, -17,  75,  80, -18, 4},
  {4, -17,  73,  82, -18, 4},
  {4, -17,  70,  85, -18, 4},
  {4, -17,  68,  87, -18, 4},
  {4, -16,  65,  89, -18, 4},
  {4, -16,  63,  91, -18, 4},
  {4, -16,  60,  94, -18, 4},
  {3, -15,  58,  96, -18, 4},
  {4, -15,  55,  98, -18, 4},
  {3, -14,  52, 100, -17, 4},
  {3, -14,  50, 102, -17, 4},
  {3, -13,  47, 104, -17, 4},
  {3, -13,  45, 106, -17, 4},
  {3, -12,  42, 108, -16, 3},
  {3, -11,  40, 109, -16, 3},
  {3, -11,  37, 111, -15, 3},
  {2, -10,  35, 113, -15, 3},
  {3, -10,  32, 114, -14, 3},
  {2,  -9,  29, 116, -13, 3},
  {2,  -8,  27, 117, -13, 3},
  {2,  -8,  25, 119, -12, 2},
  {2,  -7,  22, 120, -11, 2},
  {1,  -6,  20, 121, -10, 2},
  {1,  -6,  18, 122,  -9, 2},
  {1,  -5,  15, 123,  -8, 2},
  {1,  -4,  13, 124,  -7, 1},
  {1,  -4,  11, 125,  -6, 1},
  {1,  -3,   8, 126,  -5, 1},
  {1,  -2,   6, 126,  -4, 1},
  {0,  -1,   4, 127,  -3, 1},
  {0,  -1,   2, 128,  -1, 0},
};

static int32_t do_ntap_filter(int32_t *p, int x) {
  int i;
  int32_t sum = 0;
  for (i = 0; i < WARPEDPIXEL_FILTER_TAPS; ++i) {
    sum += p[i - WARPEDPIXEL_FILTER_TAPS / 2 + 1] * filter_ntap[x][i];
  }
  return sum;
}

static int32_t do_cubic_filter(int32_t *p, int x) {
  if (x == 0)  {
    return p[0];
  } else if (x == (1 << WARPEDPIXEL_PREC_BITS)) {
    return p[1];
  } else {
    const int64_t v1 = x * x * x * (3 * (p[0] - p[1]) + p[2] - p[-1]);
    const int64_t v2 = x * x * (2 * p[-1] - 5 * p[0] + 4 * p[1] - p[2]);
    const int64_t v3 = x * (p[1] - p[-1]);
    const int64_t v4 = 2 * p[0];
    return (int32_t)ROUND_POWER_OF_TWO_SIGNED(
        (v4 << (3 * WARPEDPIXEL_PREC_BITS)) +
        (v3 << (2 * WARPEDPIXEL_PREC_BITS)) +
        (v2 << WARPEDPIXEL_PREC_BITS) + v1,
        3 * WARPEDPIXEL_PREC_BITS + 1 - WARPEDPIXEL_FILTER_BITS);
  }
}

/*
static int32_t do_linear_filter(int32_t *p, int x) {
  int32_t sum = 0;
  sum = p[0] * (WARPEDPIXEL_PREC_SHIFTS - x) + p[1] * x;
  sum <<= (WARPEDPIXEL_FILTER_BITS - WARPEDPIXEL_PREC_BITS);
  return sum;
}

static int32_t do_4tap_filter(int32_t *p, int x) {
  int i;
  int32_t sum = 0;
  for (i = 0; i < 4; ++i) {
    sum += p[i - 1] * filter_4tap[x][i];
  }
  return sum;
}
*/

static INLINE void get_subcolumn(int taps, uint8_t *ref, int32_t *col,
                                 int stride, int x, int y_start) {
  int i;
  for (i = 0; i < taps; ++i) {
    col[i] = ref[(i + y_start) * stride + x];
  }
}

static uint8_t bi_ntap_filter(uint8_t *ref, int x, int y, int stride) {
  int32_t val, arr[WARPEDPIXEL_FILTER_TAPS];
  int k;
  int i = (int) x >> WARPEDPIXEL_PREC_BITS;
  int j = (int) y >> WARPEDPIXEL_PREC_BITS;
  for (k = 0; k < WARPEDPIXEL_FILTER_TAPS; ++k) {
    int32_t arr_temp[WARPEDPIXEL_FILTER_TAPS];
    get_subcolumn(WARPEDPIXEL_FILTER_TAPS, ref, arr_temp, stride,
                  i + k + 1 - WARPEDPIXEL_FILTER_TAPS / 2,
                  j + 1 - WARPEDPIXEL_FILTER_TAPS / 2);
    arr[k] = do_ntap_filter(arr_temp + WARPEDPIXEL_FILTER_TAPS / 2 - 1,
                            y - (j << WARPEDPIXEL_PREC_BITS));
  }
  val = do_ntap_filter(arr + WARPEDPIXEL_FILTER_TAPS / 2 - 1,
                       x - (i << WARPEDPIXEL_PREC_BITS));
  val = ROUND_POWER_OF_TWO_SIGNED(val, WARPEDPIXEL_FILTER_BITS * 2);
  return (uint8_t)clip_pixel(val);
}

static uint8_t bi_cubic_filter(uint8_t *ref, int x, int y, int stride) {
  int32_t val, arr[4];
  int k;
  int i = (int) x >> WARPEDPIXEL_PREC_BITS;
  int j = (int) y >> WARPEDPIXEL_PREC_BITS;
  for (k = 0; k < 4; ++k) {
    int32_t arr_temp[4];
    get_subcolumn(4, ref, arr_temp, stride,
                  i + k - 1, j - 1);
    arr[k] = do_cubic_filter(arr_temp + 1, y - (j << WARPEDPIXEL_PREC_BITS));
  }
  val = do_cubic_filter(arr + 1, x - (i << WARPEDPIXEL_PREC_BITS));
  val = ROUND_POWER_OF_TWO_SIGNED(val, WARPEDPIXEL_FILTER_BITS * 2);
  return (uint8_t)clip_pixel(val);
}

static uint8_t bi_linear_filter(uint8_t *ref, int x, int y, int stride) {
  const int ix = x >> WARPEDPIXEL_PREC_BITS;
  const int iy = y >> WARPEDPIXEL_PREC_BITS;
  const int sx = x - (ix << WARPEDPIXEL_PREC_BITS);
  const int sy = y - (iy << WARPEDPIXEL_PREC_BITS);
  int32_t val;
  val = ROUND_POWER_OF_TWO_SIGNED(
      ref[iy * stride + ix] * (WARPEDPIXEL_PREC_SHIFTS - sy) *
                              (WARPEDPIXEL_PREC_SHIFTS - sx) +
      ref[iy * stride + ix + 1] * (WARPEDPIXEL_PREC_SHIFTS - sy) * sx +
      ref[(iy + 1) * stride + ix] * sy * (WARPEDPIXEL_PREC_SHIFTS - sx) +
      ref[(iy + 1) * stride + ix + 1] * sy * sx,
      WARPEDPIXEL_PREC_BITS * 2);
  return (uint8_t)clip_pixel(val);
}

static uint8_t warp_interpolate(uint8_t *ref, int x, int y,
                                int width, int height, int stride) {
  int ix = x >> WARPEDPIXEL_PREC_BITS;
  int iy = y >> WARPEDPIXEL_PREC_BITS;
  int sx = x - (ix << WARPEDPIXEL_PREC_BITS);
  int sy = y - (iy << WARPEDPIXEL_PREC_BITS);
  int32_t v;

  if (ix < 0 && iy < 0) return ref[0];
  else if (ix < 0 && iy > height - 1)
    return ref[(height - 1) * stride];
  else if (ix > width - 1 && iy < 0)
    return ref[width - 1];
  else if (ix > width - 1 && iy > height - 1)
    return ref[(height - 1) * stride + (width - 1)];
  else if (ix < 0) {
    v = ROUND_POWER_OF_TWO_SIGNED(
        ref[iy * stride] * (WARPEDPIXEL_PREC_SHIFTS - sy) +
        ref[(iy + 1) * stride] * sy,
        WARPEDPIXEL_PREC_BITS);
    return clip_pixel(v);
  } else if (iy < 0) {
    v = ROUND_POWER_OF_TWO_SIGNED(
        ref[ix] * (WARPEDPIXEL_PREC_SHIFTS - sx) +
        ref[ix + 1] * sx,
        WARPEDPIXEL_PREC_BITS);
    return clip_pixel(v);
  } else if (ix > width - 1) {
    v = ROUND_POWER_OF_TWO_SIGNED(
        ref[iy * stride + width - 1] * (WARPEDPIXEL_PREC_SHIFTS - sy) +
        ref[(iy + 1) * stride + width - 1] * sy,
        WARPEDPIXEL_PREC_BITS);
    return clip_pixel(v);
  } else if (iy > height - 1) {
    v = ROUND_POWER_OF_TWO_SIGNED(
        ref[(height - 1) * stride + ix] * (WARPEDPIXEL_PREC_SHIFTS - sx) +
        ref[(height - 1) * stride + ix + 1] * sx,
        WARPEDPIXEL_PREC_BITS);
    return clip_pixel(v);
  } else if (ix >= WARPEDPIXEL_FILTER_TAPS / 2 - 1 &&
             iy >= WARPEDPIXEL_FILTER_TAPS / 2 - 1 &&
             ix < width - WARPEDPIXEL_FILTER_TAPS / 2 &&
             iy < height - WARPEDPIXEL_FILTER_TAPS / 2) {
    return bi_ntap_filter(ref, x, y, stride);
  } else if (ix >= 1 && iy >= 1 &&
             ix < width - 2 && iy < height - 2) {
    return bi_cubic_filter(ref, x, y, stride);
  } else {
    return bi_linear_filter(ref, x, y, stride);
  }
}

void vp10_warp_plane(WarpedMotionParams *wm,
                     uint8_t *ref,
                     int width, int height, int stride,
                     uint8_t *pred,
                     int p_col, int p_row,
                     int p_width, int p_height, int p_stride,
                     int subsampling_x, int subsampling_y,
                     int x_scale, int y_scale) {
  int i, j;
  projectPointsType projectPoints = get_projectPointsType(wm->wmtype);
  if (projectPoints == NULL)
    return;
  for (i = p_row; i < p_row + p_height; ++i) {
    for (j = p_col; j < p_col + p_width; ++j) {
      int in[2], out[2];
      projectPoints(wm->wmmat, in, out, 1, 2, 2, subsampling_x, subsampling_y);
      out[0] = ROUND_POWER_OF_TWO_SIGNED(out[0] * x_scale, 4);
      out[1] = ROUND_POWER_OF_TWO_SIGNED(out[1] * y_scale, 4);
      pred[(j - p_col) + (i - p_row) * p_stride] =
          warp_interpolate(ref, out[0], out[1], width, height, stride);
    }
  }
}

#if CONFIG_VP9_HIGHBITDEPTH
static INLINE void highbd_get_subcolumn(int taps, uint16_t *ref, int32_t *col,
                                        int stride, int x, int y_start) {
  int i;
  for (i = 0; i < taps; ++i) {
    col[i] = ref[(i + y_start) * stride + x];
  }
}

static uint16_t highbd_bi_ntap_filter(uint16_t *ref,
                                      int x, int y, int stride,
                                      int bd) {
  int32_t val, arr[WARPEDPIXEL_FILTER_TAPS];
  int k;
  int i = (int) x >> WARPEDPIXEL_PREC_BITS;
  int j = (int) y >> WARPEDPIXEL_PREC_BITS;
  for (k = 0; k < WARPEDPIXEL_FILTER_TAPS; ++k) {
    int32_t arr_temp[WARPEDPIXEL_FILTER_TAPS];
    highbd_get_subcolumn(WARPEDPIXEL_FILTER_TAPS, ref, arr_temp, stride,
                         i + k + 1 - WARPEDPIXEL_FILTER_TAPS / 2,
                         j + 1 - WARPEDPIXEL_FILTER_TAPS / 2);
    arr[k] = do_ntap_filter(arr_temp + WARPEDPIXEL_FILTER_TAPS / 2 - 1,
                            y - (j << WARPEDPIXEL_PREC_BITS));
  }
  val = do_ntap_filter(arr + WARPEDPIXEL_FILTER_TAPS / 2 - 1,
                       x - (i << WARPEDPIXEL_PREC_BITS));
  val = ROUND_POWER_OF_TWO_SIGNED(val, WARPEDPIXEL_FILTER_BITS * 2);
  return (uint16_t)highbd_clip_pixel(val, bd);
}

static uint16_t highbd_bi_cubic_filter(uint16_t *ref,
                                       int x, int y, int stride,
                                       int bd) {
  int32_t val, arr[4];
  int k;
  int i = (int) x >> WARPEDPIXEL_PREC_BITS;
  int j = (int) y >> WARPEDPIXEL_PREC_BITS;
  for (k = 0; k < 4; ++k) {
    int32_t arr_temp[4];
    highbd_get_subcolumn(4, ref, arr_temp, stride,
                         i + k - 1, j - 1);
    arr[k] = do_cubic_filter(arr_temp + 1, y - (j << WARPEDPIXEL_PREC_BITS));
  }
  val = do_cubic_filter(arr + 1, x - (i << WARPEDPIXEL_PREC_BITS));
  val = ROUND_POWER_OF_TWO_SIGNED(val, WARPEDPIXEL_FILTER_BITS * 2);
  return (uint16_t)highbd_clip_pixel(val, bd);
}

static uint16_t highbd_bi_linear_filter(uint16_t *ref,
                                        int x, int y, int stride,
                                        int bd) {
  const int ix = x >> WARPEDPIXEL_PREC_BITS;
  const int iy = y >> WARPEDPIXEL_PREC_BITS;
  const int sx = x - (ix << WARPEDPIXEL_PREC_BITS);
  const int sy = y - (iy << WARPEDPIXEL_PREC_BITS);
  int32_t val;
  val = ROUND_POWER_OF_TWO_SIGNED(
      ref[iy * stride + ix] * (WARPEDPIXEL_PREC_SHIFTS - sy) *
                              (WARPEDPIXEL_PREC_SHIFTS - sx) +
      ref[iy * stride + ix + 1] * (WARPEDPIXEL_PREC_SHIFTS - sy) * sx +
      ref[(iy + 1) * stride + ix] * sy * (WARPEDPIXEL_PREC_SHIFTS - sx) +
      ref[(iy + 1) * stride + ix + 1] * sy * sx,
      WARPEDPIXEL_PREC_BITS * 2);
  return (uint16_t)highbd_clip_pixel(val, bd);
}

static uint16_t highbd_warp_interpolate(uint16_t *ref,
                                        int x, int y,
                                        int width, int height, int stride,
                                        int bd) {
  int ix = x >> WARPEDPIXEL_PREC_BITS;
  int iy = y >> WARPEDPIXEL_PREC_BITS;
  int sx = x - (ix << WARPEDPIXEL_PREC_BITS);
  int sy = y - (iy << WARPEDPIXEL_PREC_BITS);
  int32_t v;

  if (ix < 0 && iy < 0) return ref[0];
  else if (ix < 0 && iy > height - 1)
    return ref[(height - 1) * stride];
  else if (ix > width - 1 && iy < 0)
    return ref[width - 1];
  else if (ix > width - 1 && iy > height - 1)
    return ref[(height - 1) * stride + (width - 1)];
  else if (ix < 0) {
    v = ROUND_POWER_OF_TWO_SIGNED(
        ref[iy * stride] * (WARPEDPIXEL_PREC_SHIFTS - sy) +
        ref[(iy + 1) * stride] * sy,
        WARPEDPIXEL_PREC_BITS);
    return highbd_clip_pixel(v, bd);
  } else if (iy < 0) {
    v = ROUND_POWER_OF_TWO_SIGNED(
        ref[ix] * (WARPEDPIXEL_PREC_SHIFTS - sx) +
        ref[ix + 1] * sx,
        WARPEDPIXEL_PREC_BITS);
    return highbd_clip_pixel(v, bd);
  } else if (ix > width - 1) {
    v = ROUND_POWER_OF_TWO_SIGNED(
        ref[iy * stride + width - 1] * (WARPEDPIXEL_PREC_SHIFTS - sy) +
        ref[(iy + 1) * stride + width - 1] * sy,
        WARPEDPIXEL_PREC_BITS);
    return highbd_clip_pixel(v, bd);
  } else if (iy > height - 1) {
    v = ROUND_POWER_OF_TWO_SIGNED(
        ref[(height - 1) * stride + ix] * (WARPEDPIXEL_PREC_SHIFTS - sx) +
        ref[(height - 1) * stride + ix + 1] * sx,
        WARPEDPIXEL_PREC_BITS);
    return highbd_clip_pixel(v, bd);
  } else if (ix >= WARPEDPIXEL_FILTER_TAPS / 2 - 1 &&
             iy >= WARPEDPIXEL_FILTER_TAPS / 2 - 1 &&
             ix < width - WARPEDPIXEL_FILTER_TAPS / 2 &&
             iy < height - WARPEDPIXEL_FILTER_TAPS / 2) {
    return highbd_bi_ntap_filter(ref, x, y, stride, bd);
  } else if (ix >= 1 && iy >= 1 &&
             ix < width - 2 && iy < height - 2) {
    return highbd_bi_cubic_filter(ref, x, y, stride, bd);
  } else {
    return highbd_bi_linear_filter(ref, x, y, stride, bd);
  }
}

void vp10_highbd_warp_plane(WarpedMotionParams *wm,
                            uint8_t *ref8,
                            int width, int height, int stride,
                            uint8_t *pred8,
                            int p_col, int p_row,
                            int p_width, int p_height, int p_stride,
                            int subsampling_col, int subsampling_row,
                            int x_scale, int y_scale,
                            int bd) {
  int i, j;
  projectPointsType projectPoints = get_projectPointsType(wm->wmtype);
  uint16_t *pred = CONVERT_TO_SHORTPTR(pred8);
  uint16_t *ref = CONVERT_TO_SHORTPTR(ref8);
  if (projectPoints == NULL)
    return;
  for (i = p_row; i < p_row + p_height; ++i) {
    for (j = p_col; j < p_col + p_width; ++j) {
      int in[2], out[2];
      projectPoints(wm->wmmat, in, out, 1, 2, 2, subsampling_x, subsampling_y);
      out[0] = ROUND_POWER_OF_TWO_SIGNED(out[0] * x_scale, 4);
      out[1] = ROUND_POWER_OF_TWO_SIGNED(out[1] * y_scale, 4);
      pred[(j - p_col) + (i - p_row) * p_stride] =
          highbd_warp_interpolate(
              ref, out[0], out[1], width, height, stride, bd);
    }
  }
}
#endif  // CONFIG_VP9_HIGHBITDEPTH

void vp10_integerize_model(double *H, TransformationType wmtype,
                           WarpedMotionParams *wm) {
  wm->wmtype = wmtype;
  switch (wmtype) {
    case HOMOGRAPHY:
      assert(fabs(H[8] - 1.0) < 1e-12);
      wm->wmmat[7] = rint(H[7] * (1 << WARPEDMODEL_ROW3HOMO_PREC_BITS));
      wm->wmmat[6] = rint(H[6] * (1 << WARPEDMODEL_ROW3HOMO_PREC_BITS));
    case AFFINE:
      wm->wmmat[5] = rint(H[5] * (1 << WARPEDMODEL_PREC_BITS));
      wm->wmmat[4] = rint(H[4] * (1 << WARPEDMODEL_PREC_BITS));
    case ROTZOOM:
      wm->wmmat[3] = rint(H[3] * (1 << WARPEDMODEL_PREC_BITS));
      wm->wmmat[2] = rint(H[2] * (1 << WARPEDMODEL_PREC_BITS));
    case TRANSLATION:
      wm->wmmat[1] = rint(H[1] * (1 << WARPEDMODEL_PREC_BITS));
      wm->wmmat[0] = rint(H[0] * (1 << WARPEDMODEL_PREC_BITS));
      break;
    default:
      assert(0);
  };
  return;
}
