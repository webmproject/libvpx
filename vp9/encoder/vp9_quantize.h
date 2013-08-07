/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#ifndef VP9_ENCODER_VP9_QUANTIZE_H_
#define VP9_ENCODER_VP9_QUANTIZE_H_

#include "vp9/encoder/vp9_block.h"

#define prototype_quantize_block(sym) \
  void (sym)(MACROBLOCK *mb, int b_idx)

#define prototype_quantize_block_pair(sym) \
  void (sym)(MACROBLOCK *mb, int b_idx1, int b_idx2)

#define prototype_quantize_mb(sym) \
  void (sym)(MACROBLOCK *x)

void vp9_regular_quantize_b_4x4_pair(MACROBLOCK *mb, int b_idx1, int b_idx2,
                                     int y_blocks);
void vp9_regular_quantize_b_4x4(MACROBLOCK *mb, int b_idx, TX_TYPE tx_type,
                                int y_blocks);
void vp9_regular_quantize_b_8x8(MACROBLOCK *mb, int b_idx, TX_TYPE tx_type,
                                int y_blocks);
struct VP9_COMP;

extern void vp9_set_quantizer(struct VP9_COMP *cpi, int Q);

extern void vp9_frame_init_quantizer(struct VP9_COMP *cpi);

extern void vp9_update_zbin_extra(struct VP9_COMP *cpi, MACROBLOCK *x);

extern void vp9_mb_init_quantizer(struct VP9_COMP *cpi, MACROBLOCK *x);

extern void vp9_init_quantizer(struct VP9_COMP *cpi);

#endif  // VP9_ENCODER_VP9_QUANTIZE_H_
