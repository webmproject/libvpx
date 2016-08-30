/*
 *  Copyright (c) 2015 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#ifndef AV1_COMMON_AV1_FWD_TXFM_H_
#define AV1_COMMON_AV1_FWD_TXFM_H_

#include "aom_dsp/txfm_common.h"
#include "aom_dsp/fwd_txfm.h"

void av1_fdct32(const tran_high_t *input, tran_high_t *output, int round);
#endif  // AV1_COMMON_AV1_FWD_TXFM_H_
