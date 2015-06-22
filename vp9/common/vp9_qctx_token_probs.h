/*
 *  Copyright (c) 2015 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#ifndef VP9_COMMON_VP9_QCTX_TOKEN_PROBS_H_
#define VP9_COMMON_VP9_QCTX_TOKEN_PROBS_H_

#include "vp9/common/vp9_entropymode.h"

#if CONFIG_QCTX_TPROBS
#define QCTX_BINS_BITS 2
extern const vp9_coeff_probs_model
default_qctx_coef_probs[1 << QCTX_BINS_BITS][TX_SIZES][PLANE_TYPES];
#endif  // CONFIG_QCTX_TPROBS

#endif  // VP9_COMMON_VP9_QCTX_TOKEN_PROBS_H_
