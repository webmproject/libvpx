/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include "./vpx_config.h"
#include "vp9/encoder/vp9_variance.h"
#include "vpx_ports/mem.h"

unsigned int vp9_get8x8var_mmx(const uint8_t *src, int src_stride,
                               const uint8_t *ref, int ref_stride,
                               unsigned int *sse, int *sum);

unsigned int vp9_get4x4var_mmx(const uint8_t *src, int src_stride,
                               const uint8_t *ref, int ref_stride,
                               unsigned int *SSE, int *sum);

unsigned int vp9_variance4x4_mmx(const uint8_t *src, int src_stride,
                                 const uint8_t *ref, int ref_stride,
                                 unsigned int *sse) {
  int sum;
  vp9_get4x4var_mmx(src, src_stride, ref, ref_stride, sse, &sum);
  return *sse - (((unsigned int)sum * sum) >> 4);
}

unsigned int vp9_variance8x8_mmx(const uint8_t *src, int src_stride,
                                 const uint8_t *ref, int ref_stride,
                                 unsigned int *sse) {
  int sum;
  vp9_get8x8var_mmx(src, src_stride, ref, ref_stride, sse, &sum);
  return *sse - (((unsigned int)sum * sum) >> 6);
}

unsigned int vp9_mse16x16_mmx(const uint8_t *src, int src_stride,
                              const uint8_t *ref, int ref_stride,
                              unsigned int *sse) {
  unsigned int sse0, sse1, sse2, sse3;
  int sum0, sum1, sum2, sum3;

  vp9_get8x8var_mmx(src, src_stride, ref, ref_stride, &sse0, &sum0);
  vp9_get8x8var_mmx(src + 8, src_stride, ref + 8, ref_stride, &sse1, &sum1);
  vp9_get8x8var_mmx(src + 8 * src_stride, src_stride,
                    ref + 8 * ref_stride, ref_stride, &sse2, &sum2);
  vp9_get8x8var_mmx(src + 8 * src_stride + 8, src_stride,
                    ref + 8 * ref_stride + 8, ref_stride, &sse3, &sum3);

  *sse = sse0 + sse1 + sse2 + sse3;
  return *sse;
}


unsigned int vp9_variance16x16_mmx(const uint8_t *src, int src_stride,
                                   const uint8_t *ref, int ref_stride,
                                   unsigned int *sse) {
  unsigned int sse0, sse1, sse2, sse3;
  int sum0, sum1, sum2, sum3, sum;

  vp9_get8x8var_mmx(src, src_stride, ref, ref_stride, &sse0, &sum0);
  vp9_get8x8var_mmx(src + 8, src_stride, ref + 8, ref_stride, &sse1, &sum1);
  vp9_get8x8var_mmx(src + 8 * src_stride, src_stride,
                    ref + 8 * ref_stride, ref_stride, &sse2, &sum2);
  vp9_get8x8var_mmx(src + 8 * src_stride + 8, src_stride,
                    ref + 8 * ref_stride + 8, ref_stride, &sse3, &sum3);

  *sse = sse0 + sse1 + sse2 + sse3;
  sum = sum0 + sum1 + sum2 + sum3;
  return *sse - (((unsigned int)sum * sum) >> 8);
}

unsigned int vp9_variance16x8_mmx(const uint8_t *src, int src_stride,
                                  const uint8_t *ref, int ref_stride,
                                  unsigned int *sse) {
  unsigned int sse0, sse1;
  int sum0, sum1, sum;

  vp9_get8x8var_mmx(src, src_stride, ref, ref_stride, &sse0, &sum0);
  vp9_get8x8var_mmx(src + 8, src_stride, ref + 8, ref_stride, &sse1, &sum1);

  *sse = sse0 + sse1;
  sum = sum0 + sum1;
  return *sse - (((unsigned int)sum * sum) >> 7);
}


unsigned int vp9_variance8x16_mmx(const uint8_t *src, int src_stride,
                                  const uint8_t *ref, int ref_stride,
                                  unsigned int *sse) {
  unsigned int sse0, sse1;
  int sum0, sum1, sum;

  vp9_get8x8var_mmx(src, src_stride, ref, ref_stride, &sse0, &sum0);
  vp9_get8x8var_mmx(src + 8 * src_stride, src_stride,
                    ref + 8 * ref_stride, ref_stride, &sse1, &sum1);

  *sse = sse0 + sse1;
  sum = sum0 + sum1;
  return *sse - (((unsigned int)sum * sum) >> 7);
}
