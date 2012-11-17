/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */


#ifndef DETOKENIZE_H
#define DETOKENIZE_H

#include "onyxd_int.h"

void vp9_reset_mb_tokens_context(MACROBLOCKD* const);

int vp9_decode_mb_tokens_4x4(VP9D_COMP* const, MACROBLOCKD* const,
                             BOOL_DECODER* const);

int vp9_decode_mb_tokens_8x8(VP9D_COMP* const, MACROBLOCKD* const,
                             BOOL_DECODER* const);

int vp9_decode_mb_tokens_16x16(VP9D_COMP* const, MACROBLOCKD* const,
                               BOOL_DECODER* const);

#endif /* DETOKENIZE_H */
