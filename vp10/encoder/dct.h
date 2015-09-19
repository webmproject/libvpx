/*
 *  Copyright (c) 2015 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */
#ifndef VP10_ENCODER_DCT_H_
#define VP10_ENCODER_DCT_H_

#include "vpx_dsp/vpx_dsp_common.h"

#ifdef __cplusplus
extern "C" {
#endif

void vp10_fdct4(const tran_low_t *input, tran_low_t *output);
void vp10_fdct8(const tran_low_t *input, tran_low_t *output);
void vp10_fdct16(const tran_low_t *input, tran_low_t *output);
void vp10_fdct32_local(const tran_low_t *input, tran_low_t *output);

#ifdef __cplusplus
}  // extern "C"
#endif

#endif  // VP10_ENCODER_DCT_H_
