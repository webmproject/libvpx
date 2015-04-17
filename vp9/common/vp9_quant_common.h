/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#ifndef VP9_COMMON_VP9_QUANT_COMMON_H_
#define VP9_COMMON_VP9_QUANT_COMMON_H_

#include <stdio.h>

#include "vpx/vpx_codec.h"
#include "vp9/common/vp9_seg_common.h"

#ifdef __cplusplus
extern "C" {
#endif

#define MINQ 0
#define MAXQ 255
#define QINDEX_RANGE (MAXQ - MINQ + 1)
#define QINDEX_BITS 8
#if CONFIG_TX_SKIP
#define TX_SKIP_Q_THRESH_INTER 64
#define TX_SKIP_Q_THRESH_INTRA 255
#define TX_SKIP_SHIFT_THRESH 0
#define PXD_QUANT_INDEX 0
#endif  // CONFIG_TX_SKIP

int16_t vp9_dc_quant(int qindex, int delta, vpx_bit_depth_t bit_depth);
int16_t vp9_ac_quant(int qindex, int delta, vpx_bit_depth_t bit_depth);

int vp9_get_qindex(const struct segmentation *seg, int segment_id,
                   int base_qindex);

static INLINE int16_t vp9_round_factor_to_round(int16_t quant,
                                                int16_t round_factor) {
  return (round_factor * quant) >> 7;
}

#if CONFIG_NEW_QUANT
#define NUQ_KNOTES 5
typedef tran_low_t dequant_val_type_nuq[NUQ_KNOTES + 1];
typedef tran_low_t cumbins_type_nuq[NUQ_KNOTES];
void vp9_get_dequant_val_nuq(int q, int band, int bd,
                             tran_low_t *dq, tran_low_t *cumbins);
tran_low_t vp9_dequant_abscoeff_nuq(int v, int q, const tran_low_t *dq);
tran_low_t vp9_dequant_coeff_nuq(int v, int q, const tran_low_t *dq);
#endif  // CONFIG_NEW_QUANT

#ifdef __cplusplus
}  // extern "C"
#endif

#endif  // VP9_COMMON_VP9_QUANT_COMMON_H_
