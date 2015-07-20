/*
 *  Copyright (c) 2014 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#ifndef VP9_ENCODER_VP9_DCT_H_
#define VP9_ENCODER_VP9_DCT_H_

#include "vp9/common/vp9_idct.h"

#ifdef __cplusplus
extern "C" {
#endif

void vp9_fdct32(const tran_high_t *input, tran_high_t *output, int round);

void vp9_highbd_fdct4x4_c(const int16_t *input, tran_low_t *output, int stride);
void vp9_highbd_fdct8x8_c(const int16_t *input, tran_low_t *output, int stride);
void vp9_highbd_fdct16x16_c(const int16_t *input, tran_low_t *output,
                            int stride);
void vp9_highbd_fdct32x32_c(const int16_t *input, tran_low_t *out, int stride);
void vp9_highbd_fdct32x32_rd_c(const int16_t *input, tran_low_t *out,
                               int stride);

#ifdef __cplusplus
}  // extern "C"
#endif

#endif  // VP9_ENCODER_VP9_DCT_H_
