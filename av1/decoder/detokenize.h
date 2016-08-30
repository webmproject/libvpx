/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#ifndef AV1_DECODER_DETOKENIZE_H_
#define AV1_DECODER_DETOKENIZE_H_

#include "av1/decoder/decoder.h"
#include "av1/common/ans.h"
#include "av1/common/scan.h"

#ifdef __cplusplus
extern "C" {
#endif

void av1_decode_palette_tokens(MACROBLOCKD *const xd, int plane, aom_reader *r);
int av1_decode_block_tokens(MACROBLOCKD *const xd, int plane,
                            const scan_order *sc, int x, int y, TX_SIZE tx_size,
                            TX_TYPE tx_type,
#if CONFIG_ANS
                            struct AnsDecoder *const r,
#else
                            aom_reader *r,
#endif  // CONFIG_ANS
                            int seg_id);

#ifdef __cplusplus
}  // extern "C"
#endif

#endif  // AV1_DECODER_DETOKENIZE_H_
