/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */


#ifndef VP9_DECODER_VP9_IDCT_BLK_H_
#define VP9_DECODER_VP9_IDCT_BLK_H_

#include "vp9/common/vp9_blockd.h"


void vp9_idct_add_lossless_c(int16_t *input, unsigned char *dest, int stride,
                             int eob);

void vp9_dc_idct_add_lossless_c(int16_t *input, unsigned char *output,
                                int stride, int dc);

void vp9_dc_idct_add_y_block_lossless_c(int16_t *q, unsigned char *pre,
                                        unsigned char *dst, int stride,
                                        const int16_t *dc);

void vp9_idct_add_y_block_lossless_c(int16_t *q, unsigned char *dst, int stride,
                                     struct macroblockd *xd);

void vp9_idct_add_uv_block_lossless_c(int16_t *q, unsigned char *dst,
                                      int stride, uint16_t *eobs);

void vp9_iht_add_c(TX_TYPE tx_type, int16_t *input, unsigned char *dest,
                   int stride, int eob);

void vp9_iht_add_8x8_c(TX_TYPE tx_type, int16_t *input, unsigned char *dest,
                       int stride, int eob);

void vp9_iht_add_16x16_c(TX_TYPE tx_type, int16_t *input, unsigned char *dest,
                         int stride, int eob);

#endif  // VP9_DECODER_VP9_IDCT_BLK_H_
