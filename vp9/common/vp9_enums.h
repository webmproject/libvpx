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

#ifdef __cplusplus
extern "C" {
#endif

#if CONFIG_EXT_CODING_UNIT_SIZE
#define CODING_UNIT_SIZE_LOG2 7
#else
#define CODING_UNIT_SIZE_LOG2 6
#endif

#define CODING_UNIT_SIZE (1 << CODING_UNIT_SIZE_LOG2)

#define MI_SIZE_LOG2 3
#define MI_BLOCK_SIZE_LOG2 (CODING_UNIT_SIZE_LOG2 - MI_SIZE_LOG2)

#define MI_SIZE (1 << MI_SIZE_LOG2)  // pixels per mi-unit
#define MI_BLOCK_SIZE (1 << MI_BLOCK_SIZE_LOG2)  // mi-units per max block

#define MI_MASK (MI_BLOCK_SIZE - 1)
#define MI_MASK_2 (MI_BLOCK_SIZE * 2 - 1)

// Bitstream profiles indicated by 2-3 bits in the uncompressed header.
// 00: Profile 0.  8-bit 4:2:0 only.
// 10: Profile 1.  8-bit 4:4:4, 4:2:2, and 4:4:0.
// 01: Profile 2.  10-bit and 12-bit color only, with 4:2:0 sampling.
// 110: Profile 3. 10-bit and 12-bit color only, with 4:2:2/4:4:4/4:4:0
//                 sampling.
// 111: Undefined profile.
typedef enum BITSTREAM_PROFILE {
  PROFILE_0,
  PROFILE_1,
  PROFILE_2,
  PROFILE_3,
  MAX_PROFILES
} BITSTREAM_PROFILE;

typedef enum BLOCK_SIZE {
  BLOCK_4X4,
  BLOCK_4X8,
  BLOCK_8X4,
  BLOCK_8X8,
  BLOCK_8X16,
  BLOCK_16X8,
  BLOCK_16X16,
  BLOCK_16X32,
  BLOCK_32X16,
  BLOCK_32X32,
  BLOCK_32X64,
  BLOCK_64X32,
  BLOCK_64X64,
#if CONFIG_EXT_CODING_UNIT_SIZE
  BLOCK_64X128,
  BLOCK_128X64,
  BLOCK_128X128,
#endif  // CONFIG_EXT_CODING_UNIT_SIZE
  BLOCK_SIZES,
  BLOCK_INVALID = BLOCK_SIZES,
  BLOCK_LARGEST = BLOCK_SIZES - 1
} BLOCK_SIZE;

#if CONFIG_EXT_PARTITION
typedef enum PARTITION_TYPE {
  PARTITION_NONE,
  PARTITION_HORZ,
  PARTITION_VERT,
  PARTITION_SPLIT,
  PARTITION_HORZ_A,  // HORZ split and the left partition is split again
  PARTITION_HORZ_B,  // HORZ split and the right partition is split again
  PARTITION_VERT_A,  // VERT split and the top partition is split again
  PARTITION_VERT_B,  // VERT split and the bottom partition is split again
  EXT_PARTITION_TYPES,
  PARTITION_TYPES = PARTITION_SPLIT + 1,
  PARTITION_INVALID = EXT_PARTITION_TYPES
} PARTITION_TYPE;
#else
typedef enum PARTITION_TYPE {
  PARTITION_NONE,
  PARTITION_HORZ,
  PARTITION_VERT,
  PARTITION_SPLIT,
  PARTITION_TYPES,
  PARTITION_INVALID = PARTITION_TYPES
} PARTITION_TYPE;
#endif

typedef char PARTITION_CONTEXT;
#define PARTITION_PLOFFSET   4  // number of probability models per block size
#if CONFIG_EXT_CODING_UNIT_SIZE
#define PARTITION_CONTEXTS (5 * PARTITION_PLOFFSET)
#else
#define PARTITION_CONTEXTS (4 * PARTITION_PLOFFSET)
#endif

// block transform size
typedef enum {
  TX_4X4 = 0,                      // 4x4 transform
  TX_8X8 = 1,                      // 8x8 transform
  TX_16X16 = 2,                    // 16x16 transform
  TX_32X32 = 3,                    // 32x32 transform
#if CONFIG_TX64X64
  TX_64X64 = 4,                    // 64x64 transform
#endif
  TX_SIZES
} TX_SIZE;

#define MAX_TX_SIZE_LOG2 (TX_SIZES + 1)
#define MAX_MIN_TX_IN_BLOCK_LOG2 MAX((CODING_UNIT_SIZE_LOG2 - \
                                      MAX_TX_SIZE_LOG2), 1)
#define MAX_MIN_TX_IN_BLOCK (1 << MAX_MIN_TX_IN_BLOCK_LOG2)

// frame transform mode
typedef enum {
  ONLY_4X4            = 0,        // only 4x4 transform used
  ALLOW_8X8           = 1,        // allow block transform size up to 8x8
  ALLOW_16X16         = 2,        // allow block transform size up to 16x16
  ALLOW_32X32         = 3,        // allow block transform size up to 32x32
#if CONFIG_TX64X64
  ALLOW_64X64         = 4,        // allow block transform size up to 32x32
#endif
  TX_MODE_SELECT,                 // transform specified for each block
  TX_MODES,
} TX_MODE;

typedef enum {
  DCT_DCT   = 0,                      // DCT  in both horizontal and vertical
  ADST_DCT  = 1,                      // ADST in vertical, DCT in horizontal
  DCT_ADST  = 2,                      // DCT  in vertical, ADST in horizontal
  ADST_ADST = 3,                      // ADST in both directions
  TX_TYPES,
#if CONFIG_EXT_TX
  FLIPADST_DCT = 4,
  DCT_FLIPADST = 5,
  FLIPADST_FLIPADST = 6,
  ADST_FLIPADST = 7,
  FLIPADST_ADST = 8,
  DST_DST = 9,
  DST_DCT = 10,
  DCT_DST = 11,
  DST_ADST = 12,
  ADST_DST = 13,
  DST_FLIPADST = 14,
  FLIPADST_DST = 15,
#if CONFIG_WAVELETS
  WAVELET1_DCT_DCT,
#endif  // CONFIG_WAVELETS
  TOTAL_TX_TYPES,
#endif  // CONFIG_EXT_TX
} TX_TYPE;

#if CONFIG_EXT_TX
typedef enum {
  NORM = 0,
  ALT1 = 1,
#if CONFIG_WAVELETS
  EXT_TX_TYPES_LARGE = 2,
#endif  // CONFIG_WAVELETS
  ALT2 = 2,
  ALT3 = 3,
  ALT4 = 4,
  ALT5 = 5,
  ALT6 = 6,
  ALT7 = 7,
  ALT8 = 8,
  ALT9 = 9,
  ALT10 = 10,
  ALT11 = 11,
  ALT12 = 12,
  ALT13 = 13,
  ALT14 = 14,
  ALT15 = 15,
  EXT_TX_TYPES
} EXT_TX_TYPE;
#endif  // CONFIG_EXT_TX

#if CONFIG_PALETTE
typedef enum {
  TWO_COLORS,
  THREE_COLORS,
  FOUR_COLORS,
  FIVE_COLORS,
  SIX_COLORS,
  SEVEN_COLORS,
  EIGHT_COLORS,
  PALETTE_SIZES
} PALETTE_SIZE;

typedef enum {
  PALETTE_COLOR_ONE,
  PALETTE_COLOR_TWO,
  PALETTE_COLOR_THREE,
  PALETTE_COLOR_FOUR,
  PALETTE_COLOR_FIVE,
  PALETTE_COLOR_SIX,
  PALETTE_COLOR_SEVEN,
  PALETTE_COLOR_EIGHT,
  PALETTE_COLORS
} PALETTE_COLOR;
#endif  // CONFIG_PALETTE

typedef enum {
  VP9_LAST_FLAG = 1 << 0,
#if CONFIG_MULTI_REF
  VP9_LAST2_FLAG = 1 << 1,
  VP9_LAST3_FLAG = 1 << 2,
  VP9_LAST4_FLAG = 1 << 3,
  VP9_GOLD_FLAG = 1 << 4,
  VP9_ALT_FLAG = 1 << 5,
#else  // CONFIG_MULTI_REF
  VP9_GOLD_FLAG = 1 << 1,
  VP9_ALT_FLAG = 1 << 2,
#endif  // CONFIG_MULTI_REF
} VP9_REFFRAME;

#ifdef __cplusplus
}  // extern "C"
#endif

#endif  // VP9_COMMON_VP9_ENUMS_H_
