/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#ifndef VP9_ENCODER_VP9_ENCODEMB_H_
#define VP9_ENCODER_VP9_ENCODEMB_H_

#include "./vpx_config.h"
#include "vp9/encoder/vp9_block.h"
#include "vp9/encoder/vp9_onyx_int.h"
#include "vp9/common/vp9_onyxc_int.h"

typedef struct {
  MB_PREDICTION_MODE mode;
  MV_REFERENCE_FRAME ref_frame;
  MV_REFERENCE_FRAME second_ref_frame;
} MODE_DEFINITION;


struct VP9_ENCODER_RTCD;
#if !CONFIG_SB8X8
void vp9_encode_inter16x16(VP9_COMMON *const cm, MACROBLOCK *x,
                           int mb_row, int mb_col);
#endif

void vp9_encode_inter16x16y(MACROBLOCK *x, int mb_row, int mb_col);

void vp9_transform_sby_32x32(MACROBLOCK *x, BLOCK_SIZE_TYPE bsize);
void vp9_transform_sby_16x16(MACROBLOCK *x, BLOCK_SIZE_TYPE bsize);
void vp9_transform_sby_8x8(MACROBLOCK *x, BLOCK_SIZE_TYPE bsize);
void vp9_transform_sby_4x4(MACROBLOCK *x, BLOCK_SIZE_TYPE bsize);
void vp9_transform_sbuv_32x32(MACROBLOCK *x, BLOCK_SIZE_TYPE bsize);
void vp9_transform_sbuv_16x16(MACROBLOCK *x, BLOCK_SIZE_TYPE bsize);
void vp9_transform_sbuv_8x8(MACROBLOCK *x, BLOCK_SIZE_TYPE bsize);
void vp9_transform_sbuv_4x4(MACROBLOCK *x, BLOCK_SIZE_TYPE bsize);

void vp9_optimize_sby(VP9_COMMON *const cm, MACROBLOCK *x,
                      BLOCK_SIZE_TYPE bsize);
void vp9_optimize_sbuv(VP9_COMMON *const cm, MACROBLOCK *x,
                       BLOCK_SIZE_TYPE bsize);

#if !CONFIG_SB8X8
void vp9_fidct_mb(VP9_COMMON *const cm, MACROBLOCK *x);
#endif

void vp9_subtract_block(int rows, int cols,
                        int16_t *diff_ptr, int diff_stride,
                        const uint8_t *src_ptr, int src_stride,
                        const uint8_t *pred_ptr, int pred_stride);
void vp9_subtract_sby(MACROBLOCK *x, BLOCK_SIZE_TYPE bsize);
void vp9_subtract_sbuv(MACROBLOCK *x, BLOCK_SIZE_TYPE bsize);
void vp9_subtract_sb(MACROBLOCK *xd, BLOCK_SIZE_TYPE bsize);

#endif  // VP9_ENCODER_VP9_ENCODEMB_H_
