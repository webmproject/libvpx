/*
 *  Copyright (c) 2014 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#ifndef VP10_ENCODER_PICKMODE_H_
#define VP10_ENCODER_PICKMODE_H_

#include "vp10/encoder/encoder.h"

#ifdef __cplusplus
extern "C" {
#endif

void vp10_pick_intra_mode(VP10_COMP *cpi, MACROBLOCK *x, RD_COST *rd_cost,
                         BLOCK_SIZE bsize, PICK_MODE_CONTEXT *ctx);

void vp10_pick_inter_mode(VP10_COMP *cpi, MACROBLOCK *x,
                         TileDataEnc *tile_data,
                         int mi_row, int mi_col, RD_COST *rd_cost,
                         BLOCK_SIZE bsize,
                         PICK_MODE_CONTEXT *ctx);

void vp10_pick_inter_mode_sub8x8(VP10_COMP *cpi, MACROBLOCK *x,
                                int mi_row, int mi_col, RD_COST *rd_cost,
                                BLOCK_SIZE bsize,
                                PICK_MODE_CONTEXT *ctx);

#ifdef __cplusplus
}  // extern "C"
#endif

#endif  // VP10_ENCODER_PICKMODE_H_
