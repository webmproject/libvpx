/*
 *  Copyright (c) 2015 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#ifndef AV1_INV_TXFM1D_H_
#define AV1_INV_TXFM1D_H_

#include "av1/common/av1_txfm.h"

#ifdef __cplusplus
extern "C" {
#endif

void av1_idct4_new(const int32_t *input, int32_t *output, const int8_t *cos_bit,
                   const int8_t *stage_range);
void av1_idct8_new(const int32_t *input, int32_t *output, const int8_t *cos_bit,
                   const int8_t *stage_range);
void av1_idct16_new(const int32_t *input, int32_t *output,
                    const int8_t *cos_bit, const int8_t *stage_range);
void av1_idct32_new(const int32_t *input, int32_t *output,
                    const int8_t *cos_bit, const int8_t *stage_range);
void av1_idct64_new(const int32_t *input, int32_t *output,
                    const int8_t *cos_bit, const int8_t *stage_range);

void av1_iadst4_new(const int32_t *input, int32_t *output,
                    const int8_t *cos_bit, const int8_t *stage_range);
void av1_iadst8_new(const int32_t *input, int32_t *output,
                    const int8_t *cos_bit, const int8_t *stage_range);
void av1_iadst16_new(const int32_t *input, int32_t *output,
                     const int8_t *cos_bit, const int8_t *stage_range);
void av1_iadst32_new(const int32_t *input, int32_t *output,
                     const int8_t *cos_bit, const int8_t *stage_range);

#ifdef __cplusplus
}
#endif

#endif  // AV1_INV_TXFM1D_H_
