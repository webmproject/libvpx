/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include <assert.h>
#include <math.h>

#include "./vpx_config.h"
#include "./vp9_rtcd.h"

#include "vp9/encoder/vp9_dct.h"
#include "vp9/encoder/vp9_dwt.h"

// Note: block length must be even for this implementation
static void analysis_53_row(int length, tran_low_t *x,
                            tran_low_t *lowpass, tran_low_t *highpass) {
  int n;
  tran_low_t r, *a, *b;

  n = length >> 1;
  b = highpass;
  a = lowpass;
  while (--n) {
    *a++ = (r = *x++) << 1;
    *b++ = *x - ((r + x[1] + 1) >> 1);
    x++;
  }
  *a = (r = *x++) << 1;
  *b = *x - r;

  n = length >> 1;
  b = highpass;
  a = lowpass;
  r = *highpass;
  while (n--) {
    *a++ += (r + (*b) + 1) >> 1;
    r = *b++;
  }
}

static void analysis_53_col(int length, tran_low_t *x,
                            tran_low_t *lowpass, tran_low_t *highpass) {
  int n;
  tran_low_t r, *a, *b;

  n = length >> 1;
  b = highpass;
  a = lowpass;
  while (--n) {
    *a++ = (r = *x++);
    *b++ = (((*x) << 1) - (r + x[1]) + 2) >> 2;
    x++;
  }
  *a = (r = *x++);
  *b = (*x - r + 1) >> 1;

  n = length >> 1;
  b = highpass;
  a = lowpass;
  r = *highpass;
  while (n--) {
    *a++ += (r + (*b) + 1) >> 1;
    r = *b++;
  }
}

static void dyadic_analyze_53(int levels, int width, int height,
                              const int16_t *x, int pitch_x,
                              tran_low_t *c, int pitch_c,
                              int dwt_scale_bits) {
  int lv, i, j, nh, nw, hh = height, hw = width;
  tran_low_t buffer[2 * DWT_MAX_LENGTH];
  for (i = 0; i < height; i++) {
    for (j = 0; j < width; j++) {
      c[i * pitch_c + j] = x[i * pitch_x + j] << dwt_scale_bits;
    }
  }
  for (lv = 0; lv < levels; lv++) {
    nh = hh;
    hh = (hh + 1) >> 1;
    nw = hw;
    hw = (hw + 1) >> 1;
    if ((nh < 2) || (nw < 2)) return;
    for (i = 0; i < nh; i++) {
      memcpy(buffer, &c[i * pitch_c], nw * sizeof(tran_low_t));
      analysis_53_row(nw, buffer, &c[i * pitch_c], &c[i * pitch_c] + hw);
    }
    for (j = 0; j < nw; j++) {
      for (i = 0; i < nh; i++)
        buffer[i + nh] = c[i * pitch_c + j];
      analysis_53_col(nh, buffer + nh, buffer, buffer + hh);
      for (i = 0; i < nh; i++)
        c[i * pitch_c + j] = buffer[i];
    }
  }
}

static void analysis_26_row(int length, tran_low_t *x,
                            tran_low_t *lowpass, tran_low_t *highpass) {
  int i, n;
  tran_low_t r, s, *a, *b;
  a = lowpass;
  b = highpass;
  for (i = length >> 1; i; i--) {
    r = *x++;
    s = *x++;
    *a++ = r + s;
    *b++ = r - s;
  }
  n = length >> 1;
  if (n >= 4) {
    a = lowpass;
    b = highpass;
    r = *lowpass;
    while (--n) {
      *b++ -= (r - a[1] + 4) >> 3;
      r = *a++;
    }
    *b -= (r - *a + 4) >> 3;
  }
}

static void analysis_26_col(int length, tran_low_t *x,
                            tran_low_t *lowpass, tran_low_t *highpass) {
  int i, n;
  tran_low_t r, s, *a, *b;
  a = lowpass;
  b = highpass;
  for (i = length >> 1; i; i--) {
    r = *x++;
    s = *x++;
    *a++ = (r + s + 1) >> 1;
    *b++ = (r - s + 1) >> 1;
  }
  n = length >> 1;
  if (n >= 4) {
    a = lowpass;
    b = highpass;
    r = *lowpass;
    while (--n) {
      *b++ -= (r - a[1] + 4) >> 3;
      r = *a++;
    }
    *b -= (r - *a + 4) >> 3;
  }
}

static void dyadic_analyze_26(int levels, int width, int height,
                              const int16_t *x, int pitch_x,
                              tran_low_t *c, int pitch_c,
                              int dwt_scale_bits) {
  int lv, i, j, nh, nw, hh = height, hw = width;
  tran_low_t buffer[2 * DWT_MAX_LENGTH];
  for (i = 0; i < height; i++) {
    for (j = 0; j < width; j++) {
      c[i * pitch_c + j] = x[i * pitch_x + j] << dwt_scale_bits;
    }
  }
  for (lv = 0; lv < levels; lv++) {
    nh = hh;
    hh = (hh + 1) >> 1;
    nw = hw;
    hw = (hw + 1) >> 1;
    if ((nh < 2) || (nw < 2)) return;
    for (i = 0; i < nh; i++) {
      memcpy(buffer, &c[i * pitch_c], nw * sizeof(tran_low_t));
      analysis_26_row(nw, buffer, &c[i * pitch_c], &c[i * pitch_c] + hw);
    }
    for (j = 0; j < nw; j++) {
      for (i = 0; i < nh; i++)
        buffer[i + nh] = c[i * pitch_c + j];
      analysis_26_col(nh, buffer + nh, buffer, buffer + hh);
      for (i = 0; i < nh; i++)
        c[i * pitch_c + j] = buffer[i];
    }
  }
}

static void analysis_97(int length, double *x,
                        double *lowpass, double *highpass) {
  static const double a_predict1 = -1.586134342;
  static const double a_update1 = -0.05298011854;
  static const double a_predict2 = 0.8829110762;
  static const double a_update2 = 0.4435068522;
  static const double s_low = 1.149604398;
  static const double s_high = 1/1.149604398;
  int i;
  double y[DWT_MAX_LENGTH];
  // Predict 1
  for (i = 1; i < length - 2; i += 2) {
    x[i] += a_predict1 * (x[i - 1] + x[i + 1]);
  }
  x[length - 1] += 2 * a_predict1 * x[length - 2];
  // Update 1
  for (i = 2; i < length; i += 2) {
    x[i] += a_update1 * (x[i - 1] + x[i + 1]);
  }
  x[0] += 2 * a_update1 * x[1];
  // Predict 2
  for (i = 1; i < length - 2; i += 2) {
    x[i] += a_predict2 * (x[i - 1] + x[i + 1]);
  }
  x[length - 1] += 2 * a_predict2 * x[length - 2];
  // Update 2
  for (i = 2; i < length; i += 2) {
    x[i] += a_update2 * (x[i - 1] + x[i + 1]);
  }
  x[0] += 2 * a_update2 * x[1];
  memcpy(y, x, sizeof(*y) * length);
  // Scale and pack
  for (i = 0; i < length / 2; i++) {
    lowpass[i] = y[2 * i] * s_low;
    highpass[i] = y[2 * i + 1] * s_high;
  }
}

static void dyadic_analyze_97(int levels, int width, int height,
                              const int16_t *x, int pitch_x,
                              tran_low_t *c, int pitch_c,
                              int dwt_scale_bits) {
  int lv, i, j, nh, nw, hh = height, hw = width;
  double buffer[2 * DWT_MAX_LENGTH];
  double y[DWT_MAX_LENGTH * DWT_MAX_LENGTH];
  for (i = 0; i < height; i++) {
    for (j = 0; j < width; j++) {
      y[i * DWT_MAX_LENGTH + j] = x[i * pitch_x + j] << dwt_scale_bits;
    }
  }
  for (lv = 0; lv < levels; lv++) {
    nh = hh;
    hh = (hh + 1) >> 1;
    nw = hw;
    hw = (hw + 1) >> 1;
    if ((nh < 2) || (nw < 2)) return;
    for (i = 0; i < nh; i++) {
      memcpy(buffer, &y[i * DWT_MAX_LENGTH], nw * sizeof(*buffer));
      analysis_97(nw, buffer, &y[i * DWT_MAX_LENGTH],
                  &y[i * DWT_MAX_LENGTH] + hw);
    }
    for (j = 0; j < nw; j++) {
      for (i = 0; i < nh; i++)
        buffer[i + nh] = y[i * DWT_MAX_LENGTH + j];
      analysis_97(nh, buffer + nh, buffer, buffer + hh);
      for (i = 0; i < nh; i++)
        y[i * DWT_MAX_LENGTH + j] = buffer[i];
    }
  }
  for (i = 0; i < height; i++) {
    for (j = 0; j < width; j++) {
      c[i * pitch_c + j] = round(y[i * DWT_MAX_LENGTH + j]);
    }
  }
}

void vp9_fdwt32x32_c(const int16_t *input, tran_low_t *output, int stride) {
#if DWT_TYPE == 26
  dyadic_analyze_26(4, 32, 32, input, stride, output, 32, 2);
#elif DWT_TYPE == 97
  dyadic_analyze_97(4, 32, 32, input, stride, output, 32, 2);
#elif DWT_TYPE == 53
  dyadic_analyze_53(4, 32, 32, input, stride, output, 32, 2);
#endif
}

void vp9_fdwtdct32x32_c(const int16_t *input, tran_low_t *output,
                        int stride) {
  const int dwt_levels = 1;
  tran_low_t buffer[16 * 16];
  int i;
  // Scales up by 2-bit from unitary
#if DWT_TYPE == 26
  dyadic_analyze_26(dwt_levels, 32, 32, input, stride, output, 32, 2);
#elif DWT_TYPE == 97
  dyadic_analyze_97(dwt_levels, 32, 32, input, stride, output, 32, 2);
#elif DWT_TYPE == 53
  dyadic_analyze_53(dwt_levels, 32, 32, input, stride, output, 32, 2);
#endif
  // 16x16 dct in LL band that is unitary
  vp9_fdct16x16_noscale(output, buffer, 32);
  // Note that the transform overall is 2-bit scaled up from unitary
  for (i = 0; i < 16; ++i) {
    memcpy(&output[i * 32], &buffer[i * 16], sizeof(buffer[0]) * 16);
  }
}

#if CONFIG_TX64X64
void vp9_fdwt64x64_c(const int16_t *input, tran_low_t *output, int stride) {
#if DWT_TYPE == 26
  dyadic_analyze_26(4, 64, 64, input, stride, output, 64, 1);
#elif DWT_TYPE == 97
  dyadic_analyze_97(4, 64, 64, input, stride, output, 64, 1);
#elif DWT_TYPE == 53
  dyadic_analyze_53(4, 64, 64, input, stride, output, 64, 1);
#endif
}

void vp9_fdwtdct64x64_c(const int16_t *input, tran_low_t *output, int stride) {
  const int dwt_levels = 1;
  tran_low_t buffer[32 * 32];
  int i;
  // Scales up by 1-bit from unitary
#if DWT_TYPE == 26
  dyadic_analyze_26(dwt_levels, 64, 64, input, stride, output, 64, 1);
#elif DWT_TYPE == 97
  dyadic_analyze_97(dwt_levels, 64, 64, input, stride, output, 64, 1);
#elif DWT_TYPE == 53
  dyadic_analyze_53(dwt_levels, 64, 64, input, stride, output, 64, 1);
#endif
  // 32x32 dct in LL band that is unitary
  vp9_fdct32x32_noscale(output, buffer, 64);
  // Note that the transform overall is 1-bit scaled up from unitary
  for (i = 0; i < 32; ++i) {
    memcpy(&output[i * 64], &buffer[i * 32], sizeof(buffer[0]) * 32);
  }
}
#endif  // CONFIG_TX64X64
