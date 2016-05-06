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

#ifdef __cplusplus
extern "C" {
#endif

void vp9_encode_sb(MACROBLOCK *x, BLOCK_SIZE bsize);
#if CONFIG_SUPERTX
void vp9_encode_sb_supertx(MACROBLOCK *x, BLOCK_SIZE bsize);
#endif
void vp9_encode_sby_pass1(MACROBLOCK *x, BLOCK_SIZE bsize);
void vp9_xform_quant_fp(MACROBLOCK *x, int plane, int block,
                        BLOCK_SIZE plane_bsize, TX_SIZE tx_size);
void vp9_xform_quant_dc(MACROBLOCK *x, int plane, int block,
                        BLOCK_SIZE plane_bsize, TX_SIZE tx_size);
void vp9_xform_quant(MACROBLOCK *x, int plane, int block,
                     BLOCK_SIZE plane_bsize, TX_SIZE tx_size);
#if CONFIG_NEW_QUANT
void vp9_xform_quant_nuq(MACROBLOCK *x, int plane, int block,
                         BLOCK_SIZE plane_bsize, TX_SIZE tx_size);
void vp9_xform_quant_dc_nuq(MACROBLOCK *x, int plane, int block,
                            BLOCK_SIZE plane_bsize, TX_SIZE tx_size);
void vp9_xform_quant_fp_nuq(MACROBLOCK *x, int plane, int block,
                            BLOCK_SIZE plane_bsize, TX_SIZE tx_size);
void vp9_xform_quant_dc_fp_nuq(MACROBLOCK *x, int plane, int block,
                               BLOCK_SIZE plane_bsize, TX_SIZE tx_size);
#endif

void vp9_subtract_plane(MACROBLOCK *x, BLOCK_SIZE bsize, int plane);

void vp9_encode_block_intra(MACROBLOCK *x, int plane, int block,
                            BLOCK_SIZE plane_bsize, TX_SIZE tx_size,
                            int8_t *skip);

void vp9_encode_intra_block_plane(MACROBLOCK *x, BLOCK_SIZE bsize, int plane,
                                  int enable_optimize_b);

#if CONFIG_SR_MODE
void inv_trfm_sr(MACROBLOCK * x, TX_SIZE tx_size,
                 int plane, int block, uint8_t *dst, int dst_stride);
void sr_downsample(int16_t *src, int src_stride, int16_t *dst, int dst_stride,
                   int w, int h);
#endif  // CONFIG_SR_MODE
#if CONFIG_TX_SKIP
void vp9_tx_identity_rect(const int16_t *input, tran_low_t *out,
                          int row, int col,
                          int stride_in, int stride_out, int shift);
void vp9_tx_identity(const int16_t *input, tran_low_t *out, int stride,
                     int bs, int shift);
#endif
#ifdef __cplusplus
}  // extern "C"
#endif

#endif  // VP9_ENCODER_VP9_ENCODEMB_H_
