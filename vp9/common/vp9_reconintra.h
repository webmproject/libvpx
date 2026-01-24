/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#ifndef VPX_VP9_COMMON_VP9_RECONINTRA_H_
#define VPX_VP9_COMMON_VP9_RECONINTRA_H_

#include "vpx/vpx_integer.h"
#include "vp9/common/vp9_blockd.h"

#ifdef __cplusplus
extern "C" {
#endif

#if VPX_ARCH_X86 || VPX_ARCH_X86_64
// max block 32x32 pixels, 64 bytes alignment for AVX512
#define MAX_PRED_ALIGNMENT_HBD 64
// 32 bytes alignment for AVX2
#define MAX_PRED_ALIGNMENT 32
#define EXTRA_ABOVE_DATA 32
#else
#define MAX_PRED_ALIGNMENT_HBD 16
#define MAX_PRED_ALIGNMENT 16
#define EXTRA_ABOVE_DATA 16
#endif

void vp9_init_intra_predictors(void);

void vp9_predict_intra_block(const MACROBLOCKD *xd, int bwl_in, TX_SIZE tx_size,
                             PREDICTION_MODE mode, const uint8_t *ref,
                             int ref_stride, uint8_t *dst, int dst_stride,
                             int aoff, int loff, int plane);
#ifdef __cplusplus
}  // extern "C"
#endif

#endif  // VPX_VP9_COMMON_VP9_RECONINTRA_H_
