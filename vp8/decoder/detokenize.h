/*
 *  Copyright (c) 2010 The VP8 project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license and patent
 *  grant that can be found in the LICENSE file in the root of the source
 *  tree. All contributing project authors may be found in the AUTHORS
 *  file in the root of the source tree.
 */


#ifndef detokenize_h
#define detokenize_h 1

#include "onyxd_int.h"

void vp8_reset_mb_tokens_context(MACROBLOCKD *x);
int vp8_decode_mb_tokens(VP8D_COMP *, MACROBLOCKD *);

#endif /* detokenize_h */
