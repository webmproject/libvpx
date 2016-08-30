/*
*  Copyright (c) 2016 The WebM project authors. All Rights Reserved.
*
*  Use of this source code is governed by a BSD-style license
*  that can be found in the LICENSE file in the root of the source
*  tree. An additional intellectual property rights grant can be found
*  in the file PATENTS.  All contributing project authors may
*  be found in the AUTHORS file in the root of the source tree.
*/

#ifndef AOM_DSP_PSNR_H_
#define AOM_DSP_PSNR_H_

#include "aom_scale/yv12config.h"

#define MAX_PSNR 100.0

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
  double psnr[4];       // total/y/u/v
  uint64_t sse[4];      // total/y/u/v
  uint32_t samples[4];  // total/y/u/v
} PSNR_STATS;

// TODO(dkovalev) change aom_sse_to_psnr signature: double -> int64_t

/*!\brief Converts SSE to PSNR
*
* Converts sum of squared errros (SSE) to peak signal-to-noise ratio (PNSR).
*
* \param[in]    samples       Number of samples
* \param[in]    peak          Max sample value
* \param[in]    sse           Sum of squared errors
*/
double aom_sse_to_psnr(double samples, double peak, double sse);
int64_t aom_get_y_sse(const YV12_BUFFER_CONFIG *a, const YV12_BUFFER_CONFIG *b);
int64_t aom_get_u_sse(const YV12_BUFFER_CONFIG *a, const YV12_BUFFER_CONFIG *b);
int64_t aom_get_v_sse(const YV12_BUFFER_CONFIG *a, const YV12_BUFFER_CONFIG *b);
#if CONFIG_AOM_HIGHBITDEPTH
int64_t aom_highbd_get_y_sse(const YV12_BUFFER_CONFIG *a,
                             const YV12_BUFFER_CONFIG *b);
int64_t v_highbd_get_u_sse(const YV12_BUFFER_CONFIG *a,
                           const YV12_BUFFER_CONFIG *b);
int64_t aom_highbd_get_v_sse(const YV12_BUFFER_CONFIG *a,
                             const YV12_BUFFER_CONFIG *b);
void aom_calc_highbd_psnr(const YV12_BUFFER_CONFIG *a,
                          const YV12_BUFFER_CONFIG *b, PSNR_STATS *psnr,
                          unsigned int bit_depth, unsigned int in_bit_depth);
#endif
void aom_calc_psnr(const YV12_BUFFER_CONFIG *a, const YV12_BUFFER_CONFIG *b,
                   PSNR_STATS *psnr);

double aom_psnrhvs(const YV12_BUFFER_CONFIG *source,
                   const YV12_BUFFER_CONFIG *dest, double *phvs_y,
                   double *phvs_u, double *phvs_v, uint32_t bd, uint32_t in_bd);

#ifdef __cplusplus
}  // extern "C"
#endif
#endif  // AOM_DSP_PSNR_H_
