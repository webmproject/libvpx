/*
 *  Copyright (c) 2015 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include <math.h>

#include "./vp9_rtcd.h"
#include "vp9/common/vp9_systemdependent.h"
#include "vp9/common/vp9_blockd.h"
#include "vp9/common/vp9_idwt.h"

// Note: block length must be even for this implementation
static void synthesis_53_row(int length,
                             tran_low_t *lowpass, tran_low_t *highpass,
                             tran_low_t *x) {
  tran_low_t r, *a, *b;
  int n;

  n = length >> 1;
  b = highpass;
  a = lowpass;
  r = *highpass;
  while (n--) {
    *a++ -= (r + (*b) + 1) >> 1;
    r = *b++;
  }

  n = length >> 1;
  b = highpass;
  a = lowpass;
  while (--n) {
    *x++ = ((r = *a++) + 1) >> 1;
    *x++ = *b++ + ((r + (*a) + 2) >> 2);
  }
  *x++ = ((r = *a) + 1) >> 1;
  *x++ = *b + ((r + 1) >> 1);
}

static void synthesis_53_col(int length,
                             tran_low_t *lowpass, tran_low_t *highpass,
                             tran_low_t *x) {
  tran_low_t r, *a, *b;
  int n;

  n = length >> 1;
  b = highpass;
  a = lowpass;
  r = *highpass;
  while (n--) {
    *a++ -= (r + (*b) + 1) >> 1;
    r = *b++;
  }

  n = length >> 1;
  b = highpass;
  a = lowpass;
  while (--n) {
    r = *a++;
    *x++ = r;
    *x++ = ((*b++) << 1) + ((r + (*a) + 1) >> 1);
  }
  *x++ = *a;
  *x++ = ((*b) << 1) + *a;
}

static void dyadic_synthesize_53(int levels, int width, int height,
                                 tran_low_t *c, int pitch_c,
                                 int16_t *x, int pitch_x,
                                 int dwt_scale_bits) {
  int th[16], tw[16], lv, i, j, nh, nw, hh = height, hw = width;
  tran_low_t buffer[2 * DWT_MAX_LENGTH];
  const int dwt_scale_rnd = 1 << (dwt_scale_bits - 1);

  th[0] = hh;
  tw[0] = hw;
  for (i = 1; i <= levels; i++) {
    th[i] = (th[i - 1] + 1) >> 1;
    tw[i] = (tw[i - 1] + 1) >> 1;
  }
  for (lv = levels - 1; lv >= 0; lv--) {
    nh = th[lv];
    nw = tw[lv];
    hh = th[lv + 1];
    hw = tw[lv + 1];
    if ((nh < 2) || (nw < 2)) continue;
    for (j = 0; j < nw; j++) {
      for (i = 0; i < nh; i++)
        buffer[i] = c[i * pitch_c + j];
      synthesis_53_col(nh, buffer, buffer + hh, buffer + nh);
      for (i = 0; i < nh; i++)
        c[i * pitch_c + j] = buffer[i + nh];
    }
    for (i = 0; i < nh; i++) {
      memcpy(buffer, &c[i * pitch_c], nw * sizeof(*buffer));
      synthesis_53_row(nw, buffer, buffer + hw, &c[i * pitch_c]);
    }
  }
  for (i = 0; i < height; i++) {
    for (j = 0; j < width; j++) {
      x[i * pitch_x + j] = c[i * pitch_c + j] >= 0 ?
          ((c[i * pitch_c + j] + dwt_scale_rnd) >> dwt_scale_bits) :
          -((-c[i * pitch_c + j] + dwt_scale_rnd) >> dwt_scale_bits);
    }
  }
}

// Note: block length must be even for this implementation
static void synthesis_26_row(int length,
                             tran_low_t *lowpass, tran_low_t *highpass,
                             tran_low_t *x) {
  tran_low_t r, s, *a, *b;
  int i, n = length >> 1;

  if (n >= 4) {
    a = lowpass;
    b = highpass;
    r = *lowpass;
    while (--n) {
      *b++ += (r - a[1] + 4) >> 3;
      r = *a++;
    }
    *b += (r - *a + 4) >> 3;
  }
  a = lowpass;
  b = highpass;
  for (i = length >> 1; i; i--) {
    s = *b++;
    r = *a++;
    *x++ = (r + s + 1) >> 1;
    *x++ = (r - s + 1) >> 1;
  }
}

static void synthesis_26_col(int length,
                             tran_low_t *lowpass, tran_low_t *highpass,
                             tran_low_t *x) {
  tran_low_t r, s, *a, *b;
  int i, n = length >> 1;

  if (n >= 4) {
    a = lowpass;
    b = highpass;
    r = *lowpass;
    while (--n) {
      *b++ += (r - a[1] + 4) >> 3;
      r = *a++;
    }
    *b += (r - *a + 4) >> 3;
  }
  a = lowpass;
  b = highpass;
  for (i = length >> 1; i; i--) {
    s = *b++;
    r = *a++;
    *x++ = r + s;
    *x++ = r - s;
  }
}

static void dyadic_synthesize_26(int levels, int width, int height,
                                 tran_low_t *c, int pitch_c,
                                 int16_t *x, int pitch_x,
                                 int dwt_scale_bits) {
  int th[16], tw[16], lv, i, j, nh, nw, hh = height, hw = width;
  tran_low_t buffer[2 * DWT_MAX_LENGTH];
  const int dwt_scale_rnd = 1 << (dwt_scale_bits - 1);

  th[0] = hh;
  tw[0] = hw;
  for (i = 1; i <= levels; i++) {
    th[i] = (th[i - 1] + 1) >> 1;
    tw[i] = (tw[i - 1] + 1) >> 1;
  }
  for (lv = levels - 1; lv >= 0; lv--) {
    nh = th[lv];
    nw = tw[lv];
    hh = th[lv + 1];
    hw = tw[lv + 1];
    if ((nh < 2) || (nw < 2)) continue;
    for (j = 0; j < nw; j++) {
      for (i = 0; i < nh; i++)
        buffer[i] = c[i * pitch_c + j];
      synthesis_26_col(nh, buffer, buffer + hh, buffer + nh);
      for (i = 0; i < nh; i++)
        c[i * pitch_c + j] = buffer[i + nh];
    }
    for (i = 0; i < nh; i++) {
      memcpy(buffer, &c[i * pitch_c], nw * sizeof(*buffer));
      synthesis_26_row(nw, buffer, buffer + hw, &c[i * pitch_c]);
    }
  }
  for (i = 0; i < height; i++) {
    for (j = 0; j < width; j++) {
      x[i * pitch_x + j] = c[i * pitch_c + j] >= 0 ?
          ((c[i * pitch_c + j] + dwt_scale_rnd) >> dwt_scale_bits) :
          -((-c[i * pitch_c + j] + dwt_scale_rnd) >> dwt_scale_bits);
    }
  }
}

static void synthesis_97(int length, double *lowpass, double *highpass,
                         double *x) {
  const double a_predict1 = -1.586134342;
  const double a_update1 = -0.05298011854;
  const double a_predict2 = 0.8829110762;
  const double a_update2 = 0.4435068522;
  const double s_low = 1.149604398;
  const double s_high = 1/1.149604398;
  const double inv_s_low = 1 / s_low;
  const double inv_s_high = 1 / s_high;
  int i;
  double y[DWT_MAX_LENGTH];
  // Undo pack and scale
  for (i = 0; i < length / 2; i++) {
    y[i * 2] = lowpass[i] * inv_s_low;
    y[i * 2 + 1] = highpass[i] * inv_s_high;
  }
  memcpy(x, y, sizeof(*y) * length);
  // Undo update 2
  for (i = 2; i < length; i += 2) {
    x[i] -= a_update2 * (x[i - 1] + x[i + 1]);
  }
  x[0] -= 2 * a_update2 * x[1];
  // Undo predict 2
  for (i = 1; i < length - 2; i += 2) {
    x[i] -= a_predict2 * (x[i - 1] + x[i + 1]);
  }
  x[length - 1] -= 2 * a_predict2 * x[length - 2];
  // Undo update 1
  for (i = 2; i < length; i += 2) {
    x[i] -= a_update1 * (x[i - 1] + x[i + 1]);
  }
  x[0] -= 2 * a_update1 * x[1];
  // Undo predict 1
  for (i = 1; i < length - 2; i += 2) {
    x[i] -= a_predict1 * (x[i - 1] + x[i + 1]);
  }
  x[length - 1] -= 2 * a_predict1 * x[length - 2];
}

static void dyadic_synthesize_97(int levels, int width, int height,
                                 tran_low_t *c, int pitch_c,
                                 int16_t *x, int pitch_x,
                                 int dwt_scale_bits) {
  int th[16], tw[16], lv, i, j, nh, nw, hh = height, hw = width;
  double buffer[2 * DWT_MAX_LENGTH];
  double y[DWT_MAX_LENGTH * DWT_MAX_LENGTH];

  for (i = 0; i < height; i++)
    for (j = 0; j < width; j++)
      y[i * DWT_MAX_LENGTH + j] = c[i * pitch_c + j];
  th[0] = hh;
  tw[0] = hw;
  for (i = 1; i <= levels; i++) {
    th[i] = (th[i - 1] + 1) >> 1;
    tw[i] = (tw[i - 1] + 1) >> 1;
  }
  for (lv = levels - 1; lv >= 0; lv--) {
    nh = th[lv];
    nw = tw[lv];
    hh = th[lv + 1];
    hw = tw[lv + 1];
    if ((nh < 2) || (nw < 2)) continue;
    for (j = 0; j < nw; j++) {
      for (i = 0; i < nh; i++)
        buffer[i] = y[i * DWT_MAX_LENGTH + j];
      synthesis_97(nh, buffer, buffer + hh, buffer + nh);
      for (i = 0; i < nh; i++)
        y[i * DWT_MAX_LENGTH + j] = buffer[i + nh];
    }
    for (i = 0; i < nh; i++) {
      memcpy(buffer, &y[i * DWT_MAX_LENGTH], nw * sizeof(*buffer));
      synthesis_97(nw, buffer, buffer + hw, &y[i * DWT_MAX_LENGTH]);
    }
  }
  for (i = 0; i < height; i++)
    for (j = 0; j < width; j++)
      x[i * pitch_x + j] = round(y[i * DWT_MAX_LENGTH + j] /
                                 (1 << dwt_scale_bits));
}

void vp9_idwt32x32_c(const tran_low_t *input, tran_low_t *output, int stride) {
  tran_low_t in[32 * 32];
  vpx_memcpy(in, input, sizeof(in));
#if DWT_TYPE == 26
  dyadic_synthesize_26(4, 32, 32, in, 32, output, stride, 2);
#elif DWT_TYPE == 97
  dyadic_synthesize_97(4, 32, 32, in, 32, output, stride, 2);
#elif DWT_TYPE == 53
  dyadic_synthesize_53(4, 32, 32, in, 32, output, stride, 2);
#endif
}

void vp9_idwtdct32x32_c(const tran_low_t *input, tran_low_t *output,
                        int stride) {
  const int dwt_levels = 1;
  tran_low_t buffer[16 * 16];
  tran_low_t buffer2[32 * 32];
  int i;
  for (i = 0; i < 32; ++i) {
    memcpy(&buffer2[i * 32], &input[i * 32], sizeof(buffer2[0]) * 32);
  }
  for (i = 0; i < 16; ++i) {
    memcpy(&buffer[i * 16], &input[i * 32], sizeof(buffer[0]) * 16);
  }
  vp9_idct16x16_noscale(buffer, buffer2, 32);

#if DWT_TYPE == 26
  dyadic_synthesize_26(dwt_levels, 32, 32, buffer2, 32, output, stride, 2);
#elif DWT_TYPE == 97
  dyadic_synthesize_97(dwt_levels, 32, 32, buffer2, 32, output, stride, 2);
#elif DWT_TYPE == 53
  dyadic_synthesize_53(dwt_levels, 32, 32, buffer2, 32, output, stride, 2);
#endif
}

void vp9_idwt32x32_add_c(const tran_low_t *input, uint8_t *dest, int stride) {
  int i, j;
  tran_low_t output[32 * 32];
  vp9_idwt32x32_c(input, output, 32);
  for (i = 0; i < 32; ++i) {
    for (j = 0; j < 32; ++j) {
      dest[j * stride + i] =
          clip_pixel_add(dest[j * stride + i], output[j * 32 + i]);
    }
  }
}

void vp9_idwtdct32x32_add_c(const tran_low_t *input, uint8_t *dest,
                            int stride) {
  int i, j;
  tran_low_t output[32 * 32];
  vp9_idwtdct32x32_c(input, output, 32);
  for (i = 0; i < 32; ++i) {
    for (j = 0; j < 32; ++j) {
      dest[j * stride + i] =
          clip_pixel_add(dest[j * stride + i], output[j * 32 + i]);
    }
  }
}

#if CONFIG_TX64X64
void vp9_idwt64x64_c(const tran_low_t *input, tran_low_t *output, int stride) {
  tran_low_t in[64 * 64];
  vpx_memcpy(in, input, sizeof(in));
#if DWT_TYPE == 26
  dyadic_synthesize_26(4, 64, 64, in, 64, output, stride, 1);
#elif DWT_TYPE == 97
  dyadic_synthesize_97(4, 64, 64, in, 64, output, stride, 1);
#elif DWT_TYPE == 53
  dyadic_synthesize_53(4, 64, 64, in, 64, output, stride, 1);
#endif
}

void vp9_idwtdct64x64_c(const tran_low_t *input, tran_low_t *output,
                        int stride) {
  const int dwt_levels = 1;
  tran_low_t buffer[32 * 32];
  tran_low_t buffer2[64 * 64];
  int i;
  for (i = 0; i < 64; ++i) {
    memcpy(&buffer2[i * 64], &input[i * 64], sizeof(buffer2[0]) * 64);
  }
  for (i = 0; i < 32; ++i) {
    memcpy(&buffer[i * 32], &input[i * 64], sizeof(buffer[0]) * 32);
  }
  vp9_idct32x32_noscale(buffer, buffer2, 64);
#if DWT_TYPE == 26
  dyadic_synthesize_26(dwt_levels, 64, 64, buffer2, 64, output, stride, 1);
#elif DWT_TYPE == 97
  dyadic_synthesize_97(dwt_levels, 64, 64, buffer2, 64, output, stride, 1);
#elif DWT_TYPE == 53
  dyadic_synthesize_53(dwt_levels, 64, 64, buffer2, 64, output, stride, 1);
#endif
}

void vp9_idwt64x64_add_c(const tran_low_t *input, uint8_t *dest, int stride) {
  int i, j;
  tran_low_t output[64 * 64];
  vp9_idwt64x64_c(input, output, 64);
  for (i = 0; i < 64; ++i) {
    for (j = 0; j < 64; ++j) {
      dest[j * stride + i] =
          clip_pixel_add(dest[j * stride + i], output[j * 64 + i]);
    }
  }
}

void vp9_idwtdct64x64_add_c(const tran_low_t *input, uint8_t *dest,
                            int stride) {
  int i, j;
  tran_low_t output[64 * 64];
  vp9_idwtdct64x64_c(input, output, 64);
  for (i = 0; i < 64; ++i) {
    for (j = 0; j < 64; ++j) {
      dest[j * stride + i] =
          clip_pixel_add(dest[j * stride + i], output[j * 64 + i]);
    }
  }
}
#endif  // CONFIG_TX64X64
