/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#ifndef VP9_COMMON_VP9_ENUMS_H_
#define VP9_COMMON_VP9_ENUMS_H_

#include "./vpx_config.h"

#define LOG2_MI_SIZE 3

#define MI_SIZE (1 << LOG2_MI_SIZE)
#define MI_MASK ((64 >> LOG2_MI_SIZE) - 1)

typedef enum BLOCK_SIZE_TYPE {
  BLOCK_SIZE_AB4X4,
  BLOCK_SIZE_SB4X8,
  BLOCK_SIZE_SB8X4,
  BLOCK_SIZE_SB8X8,
  BLOCK_SIZE_SB8X16,
  BLOCK_SIZE_SB16X8,
  BLOCK_SIZE_MB16X16,
  BLOCK_SIZE_SB16X32,
  BLOCK_SIZE_SB32X16,
  BLOCK_SIZE_SB32X32,
  BLOCK_SIZE_SB32X64,
  BLOCK_SIZE_SB64X32,
  BLOCK_SIZE_SB64X64,
  BLOCK_SIZE_TYPES
} BLOCK_SIZE_TYPE;

typedef enum PARTITION_TYPE {
  PARTITION_NONE,
  PARTITION_HORZ,
  PARTITION_VERT,
  PARTITION_SPLIT,
  PARTITION_TYPES
} PARTITION_TYPE;

#define PARTITION_PLOFFSET   4  // number of probability models per block size
#define NUM_PARTITION_CONTEXTS (4 * PARTITION_PLOFFSET)

#endif  // VP9_COMMON_VP9_ENUMS_H_
