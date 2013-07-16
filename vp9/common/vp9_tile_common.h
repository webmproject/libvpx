/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#ifndef VP9_COMMON_VP9_TILE_COMMON_H_
#define VP9_COMMON_VP9_TILE_COMMON_H_

#include "vp9/common/vp9_onyxc_int.h"

void vp9_get_tile_col_offsets(VP9_COMMON *cm, int tile_col_idx);

void vp9_get_tile_row_offsets(VP9_COMMON *cm, int tile_row_idx);

void vp9_get_tile_n_bits(int mi_cols,
                         int *min_log2_tile_cols, int *max_log2_tile_cols);

#endif  // VP9_COMMON_VP9_TILE_COMMON_H_
