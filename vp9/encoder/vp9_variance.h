/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#ifndef VP9_ENCODER_VP9_VARIANCE_H_
#define VP9_ENCODER_VP9_VARIANCE_H_

#include "vpx/vpx_integer.h"

#ifdef __cplusplus
extern "C" {
#endif

void variance(const uint8_t *a, int a_stride,
              const uint8_t *b, int b_stride,
              int  w, int  h,
              unsigned int *sse, int *sum);

#if CONFIG_VP9_HIGH
void high_variance(const uint8_t *a8, int a_stride,
                   const uint8_t *b8, int b_stride,
                   int  w, int  h,
                   unsigned int *sse, int *sum);
#endif


typedef unsigned int(*vp9_sad_fn_t)(const uint8_t *src_ptr,
                                    int source_stride,
                                    const uint8_t *ref_ptr,
                                    int ref_stride,
                                    unsigned int max_sad);

typedef unsigned int(*vp9_sad_avg_fn_t)(const uint8_t *src_ptr,
                                        int source_stride,
                                        const uint8_t *ref_ptr,
                                        int ref_stride,
                                        const uint8_t *second_pred,
                                        unsigned int max_sad);

typedef void (*vp9_sad_multi_fn_t)(const uint8_t *src_ptr,
                                   int source_stride,
                                   const uint8_t *ref_ptr,
                                   int  ref_stride,
                                   unsigned int *sad_array);

typedef void (*vp9_sad_multi_d_fn_t)(const uint8_t *src_ptr,
                                     int source_stride,
                                     const uint8_t* const ref_ptr[],
                                     int  ref_stride, unsigned int *sad_array);

typedef unsigned int (*vp9_variance_fn_t)(const uint8_t *src_ptr,
                                          int source_stride,
                                          const uint8_t *ref_ptr,
                                          int ref_stride,
                                          unsigned int *sse);

typedef unsigned int (*vp9_subpixvariance_fn_t)(const uint8_t *src_ptr,
                                                int source_stride,
                                                int xoffset,
                                                int yoffset,
                                                const uint8_t *ref_ptr,
                                                int Refstride,
                                                unsigned int *sse);

typedef unsigned int (*vp9_subp_avg_variance_fn_t)(const uint8_t *src_ptr,
                                                   int source_stride,
                                                   int xoffset,
                                                   int yoffset,
                                                   const uint8_t *ref_ptr,
                                                   int Refstride,
                                                   unsigned int *sse,
                                                   const uint8_t *second_pred);

typedef struct vp9_variance_vtable {
  vp9_sad_fn_t               sdf;
  vp9_sad_avg_fn_t           sdaf;
  vp9_variance_fn_t          vf;
  vp9_subpixvariance_fn_t    svf;
  vp9_subp_avg_variance_fn_t svaf;
  vp9_sad_multi_fn_t         sdx3f;
  vp9_sad_multi_fn_t         sdx8f;
  vp9_sad_multi_d_fn_t       sdx4df;
} vp9_variance_fn_ptr_t;

void vp9_comp_avg_pred(uint8_t *comp_pred, const uint8_t *pred, int width,
                       int height, const uint8_t *ref, int ref_stride);

#if defined(CONVERT_TO_SHORTPTR) && CONFIG_VP9_HIGH
static void high_comp_avg_pred(uint16_t *comp_pred, const uint8_t *pred8,
                               int width, int height, const uint8_t *ref8,
                               int ref_stride) {
  int i, j;
  uint16_t *pred = CONVERT_TO_SHORTPTR(pred8);
  uint16_t *ref = CONVERT_TO_SHORTPTR(ref8);
  for (i = 0; i < height; i++) {
    for (j = 0; j < width; j++) {
      const int tmp = pred[j] + ref[j];
      comp_pred[j] = ROUND_POWER_OF_TWO(tmp, 1);
    }
    comp_pred += width;
    pred += width;
    ref += ref_stride;
  }
}
#endif

#ifdef __cplusplus
}  // extern "C"
#endif

#endif  // VP9_ENCODER_VP9_VARIANCE_H_
