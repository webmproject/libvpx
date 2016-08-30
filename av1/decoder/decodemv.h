/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#ifndef AV1_DECODER_DECODEMV_H_
#define AV1_DECODER_DECODEMV_H_

#include "av1/decoder/bitreader.h"

#include "av1/decoder/decoder.h"

#ifdef __cplusplus
extern "C" {
#endif

void av1_read_mode_info(AV1Decoder *const pbi, MACROBLOCKD *xd,
#if CONFIG_SUPERTX
                        int supertx_enabled,
#endif

                        int mi_row, int mi_col, aom_reader *r, int x_mis,
                        int y_mis);

#ifdef __cplusplus
}  // extern "C"
#endif

#endif  // AV1_DECODER_DECODEMV_H_
