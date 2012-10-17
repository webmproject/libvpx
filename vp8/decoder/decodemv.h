/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */


#include "onyxd_int.h"

void vpx_decode_mb_mode_mv(VP8D_COMP* const pbi,
                           MACROBLOCKD* const xd,
                           int mb_row,
                           int mb_col,
                           BOOL_DECODER* const bc);
void vpx_decode_mode_mvs_init(VP8D_COMP* const pbi, BOOL_DECODER* const bc);
