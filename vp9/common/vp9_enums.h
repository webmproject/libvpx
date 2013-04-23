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

typedef enum BLOCK_SIZE_TYPE {
#if CONFIG_SB8X8
  BLOCK_SIZE_SB8X8,
#if CONFIG_SBSEGMENT
  BLOCK_SIZE_SB8X16,
  BLOCK_SIZE_SB16X8,
#endif
#endif
  BLOCK_SIZE_MB16X16,
#if CONFIG_SBSEGMENT
  BLOCK_SIZE_SB16X32,
  BLOCK_SIZE_SB32X16,
#endif
  BLOCK_SIZE_SB32X32,
#if CONFIG_SBSEGMENT
  BLOCK_SIZE_SB32X64,
  BLOCK_SIZE_SB64X32,
#endif
  BLOCK_SIZE_SB64X64,
} BLOCK_SIZE_TYPE;

typedef enum PARTITION_TYPE {
  PARTITION_NONE,
#if CONFIG_SBSEGMENT
  PARTITION_HORZ,
  PARTITION_VERT,
#endif
  PARTITION_SPLIT,
  PARTITION_TYPES
} PARTITION_TYPE;

#define PARTITION_PLOFFSET   4  // number of probability models per block size
#define NUM_PARTITION_CONTEXTS (2 * PARTITION_PLOFFSET)

#endif  // VP9_COMMON_VP9_ENUMS_H_
