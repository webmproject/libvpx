/*
 *  Copyright (c) 2015 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#ifndef VP10_INV_TXFM2D_C_H_
#define VP10_INV_TXFM2D_C_H_

#include "vp10/common/vp10_inv_txfm2d_cfg.h"
#ifdef __cplusplus
extern "C" {
#endif
void vp10_inv_txfm2d_add_4x4(const int32_t *input, uint16_t *output,
                             const int stride, const TXFM_2D_CFG *cfg,
                             const int bd);
void vp10_inv_txfm2d_add_8x8(const int32_t *input, uint16_t *output,
                             const int stride, const TXFM_2D_CFG *cfg,
                             const int bd);
void vp10_inv_txfm2d_add_16x16(const int32_t *input, uint16_t *output,
                               const int stride, const TXFM_2D_CFG *cfg,
                               const int bd);
void vp10_inv_txfm2d_add_32x32(const int32_t *input, uint16_t *output,
                               const int stride, const TXFM_2D_CFG *cfg,
                               const int bd);
#ifdef __cplusplus
}
#endif
#endif  // VP10_INV_TXFM2D_C_H_
