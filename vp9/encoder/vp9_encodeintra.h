/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#ifndef VP9_ENCODER_VP9_ENCODEINTRA_H_
#define VP9_ENCODER_VP9_ENCODEINTRA_H_

#include "vp9/encoder/vp9_onyx_int.h"

int vp9_encode_intra(VP9_COMP *cpi, MACROBLOCK *x, int use_16x16_pred);
void encode_block_intra(int plane, int block, BLOCK_SIZE_TYPE bsize,
                               int ss_txfrm_size, void *arg);
void vp9_encode_intra_block_y(VP9_COMMON *const cm, MACROBLOCK *mb,
                              BLOCK_SIZE_TYPE bs);
void vp9_encode_intra_block_uv(VP9_COMMON *const cm, MACROBLOCK *mb,
                               BLOCK_SIZE_TYPE bs);

#endif  // VP9_ENCODER_VP9_ENCODEINTRA_H_
