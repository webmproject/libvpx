/*
 * Copyright (c) 2016, Alliance for Open Media. All rights reserved
 *
 * This source code is subject to the terms of the BSD 2 Clause License and
 * the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
 * was not distributed with this source code in the LICENSE file, you can
 * obtain it at www.aomedia.org/license/software. If the Alliance for Open
 * Media Patent License 1.0 was not distributed with this source code in the
 * PATENTS file, you can obtain it at www.aomedia.org/license/patent.
 */

#ifndef AV1_COMMON_RECONINTRA_H_
#define AV1_COMMON_RECONINTRA_H_

#include "aom/aom_integer.h"
#include "av1/common/blockd.h"

#ifdef __cplusplus
extern "C" {
#endif

void av1_init_intra_predictors(void);

void av1_predict_intra_block(const MACROBLOCKD *xd, int bwl_in, int bhl_in,
                             TX_SIZE tx_size, PREDICTION_MODE mode,
                             const uint8_t *ref, int ref_stride, uint8_t *dst,
                             int dst_stride, int aoff, int loff, int plane);
#if CONFIG_EXT_INTRA
int av1_is_intra_filter_switchable(int angle);
#endif  // CONFIG_EXT_INTRA
#ifdef __cplusplus
}  // extern "C"
#endif

#endif  // AV1_COMMON_RECONINTRA_H_
