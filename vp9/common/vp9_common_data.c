/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include "vp9/common/vp9_common_data.h"

// Log 2 conversion lookup tables for block width and height
const int b_width_log2_lookup[BLOCK_SIZE_TYPES] =
  {0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4};
const int b_height_log2_lookup[BLOCK_SIZE_TYPES] =
  {0, 1, 0, 1, 2, 1, 2, 3, 2, 3, 4, 3, 4};
