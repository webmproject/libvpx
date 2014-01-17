/*
 *  Copyright (c) 2014 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include <assert.h>
#include <inttypes.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "vp9/common/vp9_common.h"
#include "vp9/encoder/vp9_resize.h"

#define FILTER_BITS               7

#define INTERP_TAPS               8
#define SUBPEL_BITS               5
#define SUBPEL_MASK               ((1 << SUBPEL_BITS) - 1)
#define INTERP_PRECISION_BITS     32

#define ROUND_POWER_OF_TWO(value, n) \
    (((value) + (1 << ((n) - 1))) >> (n))

typedef int16_t interp_kernel[INTERP_TAPS];

// Filters for interpolation - note this also filters integer pels.
const interp_kernel vp9_filteredinterp_filters[(1 << SUBPEL_BITS)] = {
  {-1, -8, 33, 80, 33, -8, -1, 0},
  {-1, -8, 30, 80, 35, -8, -1, 1},
  {-1, -8, 28, 80, 37, -7, -2, 1},
  {0, -8, 26, 79, 39, -7, -2, 1},
  {0, -8, 24, 79, 41, -7, -2, 1},
  {0, -8, 22, 78, 43, -6, -2, 1},
  {0, -8, 20, 78, 45, -5, -3, 1},
  {0, -8, 18, 77, 48, -5, -3, 1},
  {0, -8, 16, 76, 50, -4, -3, 1},
  {0, -8, 15, 75, 52, -3, -4, 1},
  {0, -7, 13, 74, 54, -3, -4, 1},
  {0, -7, 11, 73, 56, -2, -4, 1},
  {0, -7, 10, 71, 58, -1, -4, 1},
  {1, -7,  8, 70, 60,  0, -5, 1},
  {1, -6,  6, 68, 62,  1, -5, 1},
  {1, -6,  5, 67, 63,  2, -5, 1},
  {1, -6,  4, 65, 65,  4, -6, 1},
  {1, -5,  2, 63, 67,  5, -6, 1},
  {1, -5,  1, 62, 68,  6, -6, 1},
  {1, -5,  0, 60, 70,  8, -7, 1},
  {1, -4, -1, 58, 71, 10, -7, 0},
  {1, -4, -2, 56, 73, 11, -7, 0},
  {1, -4, -3, 54, 74, 13, -7, 0},
  {1, -4, -3, 52, 75, 15, -8, 0},
  {1, -3, -4, 50, 76, 16, -8, 0},
  {1, -3, -5, 48, 77, 18, -8, 0},
  {1, -3, -5, 45, 78, 20, -8, 0},
  {1, -2, -6, 43, 78, 22, -8, 0},
  {1, -2, -7, 41, 79, 24, -8, 0},
  {1, -2, -7, 39, 79, 26, -8, 0},
  {1, -2, -7, 37, 80, 28, -8, -1},
  {1, -1, -8, 35, 80, 30, -8, -1},
};

// Filters for factor of 2 downsampling.
static const int16_t vp9_down2_symeven_half_filter[] = {56, 12, -3, -1};
static const int16_t vp9_down2_symodd_half_filter[] = {64, 35, 0, -3};

static void interpolate(const uint8_t *const input, int inlength,
                        uint8_t *output, int outlength) {
  const int64_t delta = (((uint64_t)inlength << 32) + outlength / 2) /
      outlength;
  const int64_t offset = inlength > outlength ?
      (((int64_t)(inlength - outlength) << 31) + outlength / 2) / outlength :
      -(((int64_t)(outlength - inlength) << 31) + outlength / 2) / outlength;
  uint8_t *optr = output;
  int x, x1, x2, sum, k, int_pel, sub_pel;
  int64_t y;

  x = 0;
  y = offset;
  while ((y >> INTERP_PRECISION_BITS) < (INTERP_TAPS / 2 - 1)) {
    x++;
    y += delta;
  }
  x1 = x;
  x = outlength - 1;
  y = delta * x + offset;
  while ((y >> INTERP_PRECISION_BITS) +
         (int64_t)(INTERP_TAPS / 2) >= inlength) {
    x--;
    y -= delta;
  }
  x2 = x;
  if (x1 > x2) {
    for (x = 0, y = offset; x < outlength; ++x, y += delta) {
      const int16_t *filter;
      int_pel = y >> INTERP_PRECISION_BITS;
      sub_pel = (y >> (INTERP_PRECISION_BITS - SUBPEL_BITS)) & SUBPEL_MASK;
      filter = vp9_filteredinterp_filters[sub_pel];
      sum = 0;
      for (k = 0; k < INTERP_TAPS; ++k) {
        const int pk = int_pel - INTERP_TAPS / 2 + 1 + k;
        sum += filter[k] * input[(pk < 0 ? 0 :
                                  (pk >= inlength ? inlength - 1 : pk))];
      }
      *optr++ = clip_pixel(ROUND_POWER_OF_TWO(sum, FILTER_BITS));
    }
  } else {
    // Initial part.
    for (x = 0, y = offset; x < x1; ++x, y += delta) {
      const int16_t *filter;
      int_pel = y >> INTERP_PRECISION_BITS;
      sub_pel = (y >> (INTERP_PRECISION_BITS - SUBPEL_BITS)) & SUBPEL_MASK;
      filter = vp9_filteredinterp_filters[sub_pel];
      sum = 0;
      for (k = 0; k < INTERP_TAPS; ++k)
        sum += filter[k] * input[(int_pel - INTERP_TAPS / 2 + 1 + k < 0 ?
                                  0 :
                                  int_pel - INTERP_TAPS / 2 + 1 + k)];
      *optr++ = clip_pixel(ROUND_POWER_OF_TWO(sum, FILTER_BITS));
    }
    // Middle part.
    for (; x <= x2; ++x, y += delta) {
      const int16_t *filter;
      int_pel = y >> INTERP_PRECISION_BITS;
      sub_pel = (y >> (INTERP_PRECISION_BITS - SUBPEL_BITS)) & SUBPEL_MASK;
      filter = vp9_filteredinterp_filters[sub_pel];
      sum = 0;
      for (k = 0; k < INTERP_TAPS; ++k)
        sum += filter[k] * input[int_pel - INTERP_TAPS / 2 + 1 + k];
      *optr++ = clip_pixel(ROUND_POWER_OF_TWO(sum, FILTER_BITS));
    }
    // End part.
    for (; x < outlength; ++x, y += delta) {
      const int16_t *filter;
      int_pel = y >> INTERP_PRECISION_BITS;
      sub_pel = (y >> (INTERP_PRECISION_BITS - SUBPEL_BITS)) & SUBPEL_MASK;
      filter = vp9_filteredinterp_filters[sub_pel];
      sum = 0;
      for (k = 0; k < INTERP_TAPS; ++k)
        sum += filter[k] * input[(int_pel - INTERP_TAPS / 2 + 1 + k >=
                                  inlength ?  inlength - 1 :
                                  int_pel - INTERP_TAPS / 2 + 1 + k)];
      *optr++ = clip_pixel(ROUND_POWER_OF_TWO(sum, FILTER_BITS));
    }
  }
}

static void down2_symeven(const uint8_t *const input, int length,
                          uint8_t *output) {
  // Actual filter len = 2 * filter_len_half.
  static const int16_t *filter = vp9_down2_symeven_half_filter;
  const int filter_len_half = sizeof(vp9_down2_symeven_half_filter) / 2;
  int i, j;
  uint8_t *optr = output;
  int l1 = filter_len_half;
  int l2 = (length - filter_len_half);
  l1 += (l1 & 1);
  l2 += (l2 & 1);
  if (l1 > l2) {
    // Short input length.
    for (i = 0; i < length; i += 2) {
      int sum = (1 << (FILTER_BITS - 1));
      for (j = 0; j < filter_len_half; ++j) {
        sum += (input[(i - j < 0 ? 0 : i - j)] +
                input[(i + 1 + j >= length ? length - 1 : i + 1 + j)]) *
            filter[j];
      }
      sum >>= FILTER_BITS;
      *optr++ = clip_pixel(sum);
    }
  } else {
    // Initial part.
    for (i = 0; i < l1; i += 2) {
      int sum = (1 << (FILTER_BITS - 1));
      for (j = 0; j < filter_len_half; ++j) {
        sum += (input[(i - j < 0 ? 0 : i - j)] + input[i + 1 + j]) * filter[j];
      }
      sum >>= FILTER_BITS;
      *optr++ = clip_pixel(sum);
    }
    // Middle part.
    for (; i < l2; i += 2) {
      int sum = (1 << (FILTER_BITS - 1));
      for (j = 0; j < filter_len_half; ++j) {
        sum += (input[i - j] + input[i + 1 + j]) * filter[j];
      }
      sum >>= FILTER_BITS;
      *optr++ = clip_pixel(sum);
    }
    // End part.
    for (; i < length; i += 2) {
      int sum = (1 << (FILTER_BITS - 1));
      for (j = 0; j < filter_len_half; ++j) {
        sum += (input[i - j] +
                input[(i + 1 + j >= length ? length - 1 : i + 1 + j)]) *
            filter[j];
      }
      sum >>= FILTER_BITS;
      *optr++ = clip_pixel(sum);
    }
  }
}

static void down2_symodd(const uint8_t *const input, int length,
                         uint8_t *output) {
  // Actual filter len = 2 * filter_len_half - 1.
  static const int16_t *filter = vp9_down2_symodd_half_filter;
  const int filter_len_half = sizeof(vp9_down2_symodd_half_filter) / 2;
  int i, j;
  uint8_t *optr = output;
  int l1 = filter_len_half - 1;
  int l2 = (length - filter_len_half + 1);
  l1 += (l1 & 1);
  l2 += (l2 & 1);
  if (l1 > l2) {
    // Short input length.
    for (i = 0; i < length; i += 2) {
      int sum = (1 << (FILTER_BITS - 1)) + input[i] * filter[0];
      for (j = 1; j < filter_len_half; ++j) {
        sum += (input[(i - j < 0 ? 0 : i - j)] +
                input[(i + j >= length ? length - 1 : i + j)]) *
            filter[j];
      }
      sum >>= FILTER_BITS;
      *optr++ = clip_pixel(sum);
    }
  } else {
    // Initial part.
    for (i = 0; i < l1; i += 2) {
      int sum = (1 << (FILTER_BITS - 1)) + input[i] * filter[0];
      for (j = 1; j < filter_len_half; ++j) {
        sum += (input[(i - j < 0 ? 0 : i - j)] + input[i + j]) * filter[j];
      }
      sum >>= FILTER_BITS;
      *optr++ = clip_pixel(sum);
    }
    // Middle part.
    for (; i < l2; i += 2) {
      int sum = (1 << (FILTER_BITS - 1)) + input[i] * filter[0];
      for (j = 1; j < filter_len_half; ++j) {
        sum += (input[i - j] + input[i + j]) * filter[j];
      }
      sum >>= FILTER_BITS;
      *optr++ = clip_pixel(sum);
    }
    // End part.
    for (; i < length; i += 2) {
      int sum = (1 << (FILTER_BITS - 1)) + input[i] * filter[0];
      for (j = 1; j < filter_len_half; ++j) {
        sum += (input[i - j] + input[(i + j >= length ? length - 1 : i + j)]) *
            filter[j];
      }
      sum >>= FILTER_BITS;
      *optr++ = clip_pixel(sum);
    }
  }
}

static int get_down2_length(int length, int steps) {
  int s;
  for (s = 0; s < steps; ++s)
    length = (length + 1) >> 1;
  return length;
}

int get_down2_steps(int in_length, int out_length) {
  int steps = 0;
  int proj_in_length;
  while ((proj_in_length = get_down2_length(in_length, 1)) >= out_length) {
    ++steps;
    in_length = proj_in_length;
  }
  return steps;
}

static void resize_multistep(const uint8_t *const input,
                             int length,
                             uint8_t *output,
                             int olength,
                             uint8_t *buf) {
  int steps;
  if (length == olength) {
    memcpy(output, input, sizeof(uint8_t) * length);
    return;
  }
  steps = get_down2_steps(length, olength);

  if (steps > 0) {
    int s;
    uint8_t *out = NULL;
    uint8_t *tmpbuf = NULL;
    uint8_t *otmp, *otmp2;
    int filteredlength = length;
    if (!tmpbuf) {
      tmpbuf = (uint8_t *)malloc(sizeof(uint8_t) * length);
      otmp = tmpbuf;
    } else {
      otmp = buf;
    }
    otmp2 = otmp + get_down2_length(length, 1);
    for (s = 0; s < steps; ++s) {
      const int proj_filteredlength = get_down2_length(filteredlength, 1);
      const uint8_t *const in = (s == 0 ? input : out);
      if (s == steps - 1 && proj_filteredlength == olength)
        out = output;
      else
        out = (s & 1 ? otmp2 : otmp);
      if (filteredlength & 1)
        down2_symodd(in, filteredlength, out);
      else
        down2_symeven(in, filteredlength, out);
      filteredlength = proj_filteredlength;
    }
    if (filteredlength != olength) {
      interpolate(out, filteredlength, output, olength);
    }
    if (tmpbuf)
      free(tmpbuf);
  } else {
    interpolate(input, length, output, olength);
  }
}

static void fill_col_to_arr(uint8_t *img, int stride, int len, uint8_t *arr) {
  int i;
  uint8_t *iptr = img;
  uint8_t *aptr = arr;
  for (i = 0; i < len; ++i, iptr += stride) {
    *aptr++ = *iptr;
  }
}

static void fill_arr_to_col(uint8_t *img, int stride, int len, uint8_t *arr) {
  int i;
  uint8_t *iptr = img;
  uint8_t *aptr = arr;
  for (i = 0; i < len; ++i, iptr += stride) {
    *iptr = *aptr++;
  }
}

void vp9_resize_plane(const uint8_t *const input,
                      int height,
                      int width,
                      int in_stride,
                      uint8_t *output,
                      int height2,
                      int width2,
                      int out_stride) {
  int i;
  uint8_t *intbuf = (uint8_t *)malloc(sizeof(uint8_t) * width2 * height);
  uint8_t *tmpbuf = (uint8_t *)malloc(sizeof(uint8_t) *
                                      (width < height ? height : width));
  uint8_t *arrbuf = (uint8_t *)malloc(sizeof(uint8_t) * (height + height2));
  for (i = 0; i < height; ++i)
    resize_multistep(input + in_stride * i, width,
                        intbuf + width2 * i, width2, tmpbuf);
  for (i = 0; i < width2; ++i) {
    fill_col_to_arr(intbuf + i, width2, height, arrbuf);
    resize_multistep(arrbuf, height, arrbuf + height, height2, tmpbuf);
    fill_arr_to_col(output + i, out_stride, height2, arrbuf + height);
  }
  free(intbuf);
  free(tmpbuf);
  free(arrbuf);
}

void vp9_resize_frame420(const uint8_t *const y,
                         int y_stride,
                         const uint8_t *const u, const uint8_t *const v,
                         int uv_stride,
                         int height, int width,
                         uint8_t *oy, int oy_stride,
                         uint8_t *ou, uint8_t *ov, int ouv_stride,
                         int oheight, int owidth) {
  vp9_resize_plane(y, height, width, y_stride,
                   oy, oheight, owidth, oy_stride);
  vp9_resize_plane(u, height / 2, width / 2, uv_stride,
                   ou, oheight / 2, owidth / 2, ouv_stride);
  vp9_resize_plane(v, height / 2, width / 2, uv_stride,
                   ov, oheight / 2, owidth / 2, ouv_stride);
}

void vp9_resize_frame422(const uint8_t *const y, int y_stride,
                         const uint8_t *const u, const uint8_t *const v,
                         int uv_stride,
                         int height, int width,
                         uint8_t *oy, int oy_stride,
                         uint8_t *ou, uint8_t *ov, int ouv_stride,
                         int oheight, int owidth) {
  vp9_resize_plane(y, height, width, y_stride,
                   oy, oheight, owidth, oy_stride);
  vp9_resize_plane(u, height, width / 2, uv_stride,
                   ou, oheight, owidth / 2, ouv_stride);
  vp9_resize_plane(v, height, width / 2, uv_stride,
                   ov, oheight, owidth / 2, ouv_stride);
}

void vp9_resize_frame444(const uint8_t *const y, int y_stride,
                         const uint8_t *const u, const uint8_t *const v,
                         int uv_stride,
                         int height, int width,
                         uint8_t *oy, int oy_stride,
                         uint8_t *ou, uint8_t *ov, int ouv_stride,
                         int oheight, int owidth) {
  vp9_resize_plane(y, height, width, y_stride,
                   oy, oheight, owidth, oy_stride);
  vp9_resize_plane(u, height, width, uv_stride,
                   ou, oheight, owidth, ouv_stride);
  vp9_resize_plane(v, height, width, uv_stride,
                   ov, oheight, owidth, ouv_stride);
}
