/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */


#ifndef VP10_DECODER_DETOKENIZE_H_
#define VP10_DECODER_DETOKENIZE_H_

#include "vp10/decoder/decoder.h"
#include "vp10/common/ans.h"
#include "vp10/common/scan.h"

#ifdef __cplusplus
extern "C" {
#endif

void vp10_decode_palette_tokens(MACROBLOCKD *const xd, int plane,
                                vpx_reader *r);
int vp10_decode_block_tokens(MACROBLOCKD *const xd,
#if CONFIG_ANS
                             const rans_dec_lut *const token_tab,
#endif  // CONFIG_ANS
                             int plane, const scan_order *sc,
                             int x, int y,
                             TX_SIZE tx_size,
#if CONFIG_ANS
                             struct AnsDecoder *const r,
#else
                             vpx_reader *r,
#endif  // CONFIG_ANS
                             int seg_id);

#ifdef __cplusplus
}  // extern "C"
#endif

#endif  // VP10_DECODER_DETOKENIZE_H_
