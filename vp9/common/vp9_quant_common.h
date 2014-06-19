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

#include "vpx/vpx_codec.h"
#include "vp9/common/vp9_blockd.h"

#ifdef __cplusplus
extern "C" {
#endif

#define MINQ 0
#define MAXQ 255
#define MAXQ_10 327
#define MAXQ_12 398
#define QINDEX_RANGE (MAXQ - MINQ + 1)
#define QINDEX_RANGE_10 (MAXQ_10 - MINQ + 1)
#define QINDEX_RANGE_12 (MAXQ_12 - MINQ + 1)
#define QINDEX_BITS 8
#define QINDEX_BITS_10 9
#define QINDEX_BITS_12 9

#if CONFIG_VP9_HIGH && CONFIG_HIGH_TRANSFORMS && CONFIG_HIGH_QUANT
#define QINDEX_RANGE_MAX QINDEX_RANGE_12
#else
#define QINDEX_RANGE_MAX QINDEX_RANGE
#endif

void vp9_init_quant_tables();

int16_t vp9_dc_quant(int qindex, int delta, vpx_bit_depth_t bit_depth);
int16_t vp9_ac_quant(int qindex, int delta, vpx_bit_depth_t bit_depth);

int vp9_get_qindex(const struct segmentation *seg, int segment_id,
                   int base_qindex, vpx_bit_depth_t bit_depth);

int vp9_get_maxq(vpx_bit_depth_t bit_depth);

int vp9_get_qindex_range(vpx_bit_depth_t bit_depth);

int vp9_get_qindex_bits(vpx_bit_depth_t bit_depth);

#ifdef __cplusplus
}  // extern "C"
#endif

#endif  // VP9_COMMON_VP9_QUANT_COMMON_H_
