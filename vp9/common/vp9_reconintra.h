/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#ifndef VP9_COMMON_VP9_RECONINTRA_H_
#define VP9_COMMON_VP9_RECONINTRA_H_

#include "vpx/vpx_integer.h"
#include "vp9/common/vp9_blockd.h"

#ifdef __cplusplus
extern "C" {
#endif

void vp9_predict_intra_block(const MACROBLOCKD *xd, int block_idx, int bwl_in,
                             TX_SIZE tx_size, PREDICTION_MODE mode,
#if CONFIG_FILTERINTRA
                             int filterbit,
#endif
                             const uint8_t *ref, int ref_stride,
                             uint8_t *dst, int dst_stride,
                             int aoff, int loff, int plane);
#if CONFIG_INTERINTRA
void vp9_build_interintra_predictors(MACROBLOCKD *xd,
                                     uint8_t *ypred,
                                     uint8_t *upred,
                                     uint8_t *vpred,
                                     int ystride,
                                     int ustride,
                                     int vstride,
                                     BLOCK_SIZE bsize);
void vp9_build_interintra_predictors_sby(MACROBLOCKD *xd,
                                         uint8_t *ypred,
                                         int ystride,
                                         BLOCK_SIZE bsize);
void vp9_build_interintra_predictors_sbuv(MACROBLOCKD *xd,
                                          uint8_t *upred,
                                          uint8_t *vpred,
                                          int ustride, int vstride,
                                          BLOCK_SIZE bsize);
#if CONFIG_MASKED_INTERINTRA
void vp9_generate_masked_weight_interintra(int mask_index,
                                           BLOCK_SIZE sb_type,
                                           int h, int w,
                                           uint8_t *mask, int stride);
#endif
#endif
#ifdef __cplusplus
}  // extern "C"
#endif

#endif  // VP9_COMMON_VP9_RECONINTRA_H_
