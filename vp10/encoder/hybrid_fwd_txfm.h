/*
 *  Copyright (c) 2015 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#ifndef VP10_ENCODER_HYBRID_FWD_TXFM_H_
#define VP10_ENCODER_HYBRID_FWD_TXFM_H_

#include "./vpx_config.h"

typedef enum FWD_TXFM_OPT { FWD_TXFM_OPT_NORMAL, FWD_TXFM_OPT_DC } FWD_TXFM_OPT;

typedef struct FWD_TXFM_PARAM {
  TX_TYPE tx_type;
  TX_SIZE tx_size;
  FWD_TXFM_OPT fwd_txfm_opt;
  int rd_transform;
  int lossless;
#if CONFIG_VP9_HIGHBITDEPTH
  int bd;
#endif  // CONFIG_VP9_HIGHBITDEPTH
} FWD_TXFM_PARAM;

#ifdef __cplusplus
extern "C" {
#endif

void fwd_txfm(const int16_t *src_diff, tran_low_t *coeff, int diff_stride,
              FWD_TXFM_PARAM *fwd_txfm_param);

#if CONFIG_VP9_HIGHBITDEPTH
void highbd_fwd_txfm(const int16_t *src_diff, tran_low_t *coeff,
                     int diff_stride, FWD_TXFM_PARAM *fwd_txfm_param);
#endif  // CONFIG_VP9_HIGHBITDEPTH

#ifdef __cplusplus
}  // extern "C"
#endif

#endif  // VP10_ENCODER_HYBRID_FWD_TXFM_H_
