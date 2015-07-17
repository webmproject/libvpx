/*
 *  Copyright (c) 2015 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#ifndef VP9_COMMON_VP9_IDWT_H_
#define VP9_COMMON_VP9_IDWT_H_

#include <assert.h>

#include "./vpx_config.h"
#include "vp9/common/vp9_common.h"
#include "vp9/common/vp9_enums.h"
#include "vp9/common/vp9_idct.h"

#define DWT_MAX_LENGTH   64
#define DWT_TYPE         53    // 26/53/97

#ifdef __cplusplus
extern "C" {
#endif

#if CONFIG_TX64X64
void vp9_idwt64x64(tran_low_t *input, tran_low_t *output, int stride);
void vp9_idwtdct64x64(tran_low_t *input, tran_low_t *output, int stride);
#endif  // CONFIG_TX64X64
void vp9_idwt32x32(tran_low_t *input, tran_low_t *output, int stride);
void vp9_idwtdct32x32(tran_low_t *input, tran_low_t *output, int stride);

#ifdef __cplusplus
}  // extern "C"
#endif

#endif  // VP9_COMMON_VP9_IDWT_H_
