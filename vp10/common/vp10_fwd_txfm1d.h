/*
 *  Copyright (c) 2015 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#ifndef VP10_FWD_TXFM1D_H_
#define VP10_FWD_TXFM1D_H_

#include "vp10/common/vp10_txfm.h"

#ifdef __cplusplus
extern "C" {
#endif

void vp10_fdct4_new(const int32_t *input, int32_t *output,
                    const int8_t *cos_bit, const int8_t *stage_range);
void vp10_fdct8_new(const int32_t *input, int32_t *output,
                    const int8_t *cos_bit, const int8_t *stage_range);
void vp10_fdct16_new(const int32_t *input, int32_t *output,
                     const int8_t *cos_bit, const int8_t *stage_range);
void vp10_fdct32_new(const int32_t *input, int32_t *output,
                     const int8_t *cos_bit, const int8_t *stage_range);
void vp10_fdct64_new(const int32_t *input, int32_t *output,
                     const int8_t *cos_bit, const int8_t *stage_range);

void vp10_fadst4_new(const int32_t *input, int32_t *output,
                     const int8_t *cos_bit, const int8_t *stage_range);
void vp10_fadst8_new(const int32_t *input, int32_t *output,
                     const int8_t *cos_bit, const int8_t *stage_range);
void vp10_fadst16_new(const int32_t *input, int32_t *output,
                      const int8_t *cos_bit, const int8_t *stage_range);
void vp10_fadst32_new(const int32_t *input, int32_t *output,
                      const int8_t *cos_bit, const int8_t *stage_range);

#ifdef __cplusplus
}
#endif

#endif  // VP10_FWD_TXFM1D_H_
