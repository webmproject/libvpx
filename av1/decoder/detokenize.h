/*
 * Copyright (c) 2016, Alliance for Open Media. All rights reserved
 *
 * This source code is subject to the terms of the BSD 2 Clause License and
 * the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
 * was not distributed with this source code in the LICENSE file, you can
 * obtain it at www.aomedia.org/license/software. If the Alliance for Open
 * Media Patent License 1.0 was not distributed with this source code in the
 * PATENTS file, you can obtain it at www.aomedia.org/license/patent.
 */

#if !CONFIG_PVQ
#ifndef AV1_DECODER_DETOKENIZE_H_
#define AV1_DECODER_DETOKENIZE_H_

#include "av1/decoder/decoder.h"
#if CONFIG_ANS
#include "aom_dsp/ans.h"
#endif  // CONFIG_ANS
#include "av1/common/scan.h"

#ifdef __cplusplus
extern "C" {
#endif

#if CONFIG_PALETTE
void av1_decode_palette_tokens(MACROBLOCKD *const xd, int plane, aom_reader *r);
#endif  // CONFIG_PALETTE

int av1_decode_block_tokens(MACROBLOCKD *const xd, int plane,
                            const SCAN_ORDER *sc, int x, int y, TX_SIZE tx_size,
                            TX_TYPE tx_type, int16_t *max_scan_line,
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
#endif
