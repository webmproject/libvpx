/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#ifndef AV1_ENCODER_ENCODEMB_H_
#define AV1_ENCODER_ENCODEMB_H_

#include "./aom_config.h"
#include "av1/encoder/block.h"

#ifdef __cplusplus
extern "C" {
#endif

struct optimize_ctx {
  ENTROPY_CONTEXT ta[MAX_MB_PLANE][2 * MAX_MIB_SIZE];
  ENTROPY_CONTEXT tl[MAX_MB_PLANE][2 * MAX_MIB_SIZE];
};

struct encode_b_args {
  MACROBLOCK *x;
  struct optimize_ctx *ctx;
  int8_t *skip;
  ENTROPY_CONTEXT *ta;
  ENTROPY_CONTEXT *tl;
  int8_t enable_optimize_b;
};

typedef enum AV1_XFORM_QUANT {
  AV1_XFORM_QUANT_FP = 0,
  AV1_XFORM_QUANT_B = 1,
  AV1_XFORM_QUANT_DC = 2,
  AV1_XFORM_QUANT_SKIP_QUANT = 3,
  AV1_XFORM_QUANT_LAST = 4
} AV1_XFORM_QUANT;

void av1_encode_sb(MACROBLOCK *x, BLOCK_SIZE bsize);
#if CONFIG_SUPERTX
void av1_encode_sb_supertx(MACROBLOCK *x, BLOCK_SIZE bsize);
#endif  // CONFIG_SUPERTX
void av1_encode_sby_pass1(MACROBLOCK *x, BLOCK_SIZE bsize);
void av1_xform_quant(MACROBLOCK *x, int plane, int block, int blk_row,
                     int blk_col, BLOCK_SIZE plane_bsize, TX_SIZE tx_size,
                     AV1_XFORM_QUANT xform_quant_idx);
#if CONFIG_NEW_QUANT
void av1_xform_quant_nuq(MACROBLOCK *x, int plane, int block, int blk_row,
                         int blk_col, BLOCK_SIZE plane_bsize, TX_SIZE tx_size,
                         int ctx);
void av1_xform_quant_dc_nuq(MACROBLOCK *x, int plane, int block, int blk_row,
                            int blk_col, BLOCK_SIZE plane_bsize,
                            TX_SIZE tx_size, int ctx);
void av1_xform_quant_fp_nuq(MACROBLOCK *x, int plane, int block, int blk_row,
                            int blk_col, BLOCK_SIZE plane_bsize,
                            TX_SIZE tx_size, int ctx);
void av1_xform_quant_dc_fp_nuq(MACROBLOCK *x, int plane, int block, int blk_row,
                               int blk_col, BLOCK_SIZE plane_bsize,
                               TX_SIZE tx_size, int ctx);
#endif

int av1_optimize_b(MACROBLOCK *mb, int plane, int block, TX_SIZE tx_size,
                   int ctx);

void av1_subtract_plane(MACROBLOCK *x, BLOCK_SIZE bsize, int plane);

void av1_encode_block_intra(int plane, int block, int blk_row, int blk_col,
                            BLOCK_SIZE plane_bsize, TX_SIZE tx_size, void *arg);

void av1_encode_intra_block_plane(MACROBLOCK *x, BLOCK_SIZE bsize, int plane,
                                  int enable_optimize_b);

#ifdef __cplusplus
}  // extern "C"
#endif

#endif  // AV1_ENCODER_ENCODEMB_H_
