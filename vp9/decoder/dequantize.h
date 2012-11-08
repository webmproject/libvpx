/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */


#ifndef DEQUANTIZE_H
#define DEQUANTIZE_H
#include "vp9/common/blockd.h"

#if CONFIG_LOSSLESS
extern void vp9_dequant_idct_add_lossless_c(short *input, short *dq,
                                            unsigned char *pred,
                                            unsigned char *output,
                                            int pitch, int stride);
extern void vp9_dequant_dc_idct_add_lossless_c(short *input, short *dq,
                                               unsigned char *pred,
                                               unsigned char *output,
                                               int pitch, int stride, int dc);
extern void vp9_dequant_dc_idct_add_y_block_lossless_c(short *q, short *dq,
                                                       unsigned char *pre,
                                                       unsigned char *dst,
                                                       int stride,
                                                       unsigned short *eobs,
                                                       short *dc);
extern void vp9_dequant_idct_add_y_block_lossless_c(short *q, short *dq,
                                                    unsigned char *pre,
                                                    unsigned char *dst,
                                                    int stride,
                                                    unsigned short *eobs);
extern void vp9_dequant_idct_add_uv_block_lossless_c(short *q, short *dq,
                                                     unsigned char *pre,
                                                     unsigned char *dst_u,
                                                     unsigned char *dst_v,
                                                     int stride,
                                                     unsigned short *eobs);
#endif

typedef void (*vp9_dequant_idct_add_fn_t)(short *input, short *dq,
    unsigned char *pred, unsigned char *output, int pitch, int stride);
typedef void(*vp9_dequant_dc_idct_add_fn_t)(short *input, short *dq,
    unsigned char *pred, unsigned char *output, int pitch, int stride, int dc);

typedef void(*vp9_dequant_dc_idct_add_y_block_fn_t)(short *q, short *dq,
    unsigned char *pre, unsigned char *dst, int stride, unsigned short *eobs,
    short *dc);
typedef void(*vp9_dequant_idct_add_y_block_fn_t)(short *q, short *dq,
    unsigned char *pre, unsigned char *dst, int stride, unsigned short *eobs);
typedef void(*vp9_dequant_idct_add_uv_block_fn_t)(short *q, short *dq,
    unsigned char *pre, unsigned char *dst_u, unsigned char *dst_v, int stride,
    unsigned short *eobs);

void vp9_ht_dequant_idct_add_c(TX_TYPE tx_type, short *input, short *dq,
                                    unsigned char *pred, unsigned char *dest,
                                    int pitch, int stride);

void vp9_ht_dequant_idct_add_8x8_c(TX_TYPE tx_type, short *input, short *dq,
                                   unsigned char *pred, unsigned char *dest,
                                   int pitch, int stride);

void vp9_ht_dequant_idct_add_16x16_c(TX_TYPE tx_type, short *input, short *dq,
                                     unsigned char *pred, unsigned char *dest,
                                     int pitch, int stride);

#if CONFIG_SUPERBLOCKS
void vp9_dequant_dc_idct_add_y_block_8x8_inplace_c(short *q, short *dq,
                                                   unsigned char *dst,
                                                   int stride,
                                                   unsigned short *eobs,
                                                   short *dc, MACROBLOCKD *xd);

void vp9_dequant_dc_idct_add_y_block_4x4_inplace_c(short *q, short *dq,
                                                   unsigned char *dst,
                                                   int stride, char *eobs,
                                                   short *dc, MACROBLOCKD *xd);

void vp9_dequant_idct_add_uv_block_8x8_inplace_c(short *q, short *dq,
                                                 unsigned char *dstu,
                                                 unsigned char *dstv,
                                                 int stride,
                                                 unsigned short *eobs,
                                                 MACROBLOCKD *xd);

void vp9_dequant_idct_add_uv_block_4x4_inplace_c(short *q, short *dq,
                                                 unsigned char *dstu,
                                                 unsigned char *dstv,
                                                 int stride, char *eobs,
                                                 MACROBLOCKD *xd);
#endif

#endif
