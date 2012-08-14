/*
 *  Copyright (c) 2012 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include "vpx_config.h"

#if CONFIG_ROTATION
typedef struct {
  int y;
  int x;
  unsigned long t;
} tap;

typedef struct {
  tap pt[4];
} point_taps;

typedef struct {
  point_taps pt[256];
} mb_taps;

mb_taps mt_8x8[] = {
#include "rotate2.h"
};

mb_taps mt[] = {
#include "rotate.h"
};

void predict_rotated_16x16(int rotation_index, unsigned char *src, int sp,
                           unsigned char *dst, int dp) {
  int i, j, k, p = 0;

  for (i = 0; i < 16; i++, dst += dp) {
    for (j = 0; j < 16; j++, p++) {
      unsigned int sum = 32768;

      for (k = 0; k < 4; k++) {
        tap *tp = &mt[rotation_index].pt[p].pt[k];
        sum += src[tp->y * sp + tp->x] * tp->t;
      }
      sum >>= 16;
      dst[j] = sum;
    }
  }
}
void predict_rotated_8x8(int rotation_index, unsigned char *src, int sp,
                         unsigned char *dst, int dp) {
  int i, j, k, p = 0;

  for (i = 0; i < 8; i++, dst += dp) {
    for (j = 0; j < 8; j++, p++) {
      unsigned int sum = 32768;

      for (k = 0; k < 4; k++) {
        tap *tp = &mt_8x8[rotation_index].pt[p].pt[k];
        sum += src[tp->y * sp + tp->x] * tp->t;
      }
      sum >>= 16;
      dst[j] = sum;
    }
  }
}
#endif




