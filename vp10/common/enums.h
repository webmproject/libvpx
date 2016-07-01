/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#ifndef VP10_COMMON_ENUMS_H_
#define VP10_COMMON_ENUMS_H_

#include "./vpx_config.h"
#include "vpx/vpx_integer.h"

#ifdef __cplusplus
extern "C" {
#endif

#undef MAX_SB_SIZE

// Max superblock size
#if CONFIG_EXT_PARTITION
# define MAX_SB_SIZE_LOG2 7
#else
# define MAX_SB_SIZE_LOG2 6
#endif  // CONFIG_EXT_PARTITION
#define MAX_SB_SIZE   (1 << MAX_SB_SIZE_LOG2)
#define MAX_SB_SQUARE (MAX_SB_SIZE * MAX_SB_SIZE)

// Min superblock size
#define MIN_SB_SIZE_LOG2 6

// Pixels per Mode Info (MI) unit
#define MI_SIZE_LOG2  3
#define MI_SIZE       (1 << MI_SIZE_LOG2)

// MI-units per max superblock (MI Block - MIB)
#define MAX_MIB_SIZE_LOG2 (MAX_SB_SIZE_LOG2 - MI_SIZE_LOG2)
#define MAX_MIB_SIZE      (1 << MAX_MIB_SIZE_LOG2)

// MI-units per min superblock
#define MIN_MIB_SIZE_LOG2 (MIN_SB_SIZE_LOG2 - MI_SIZE_LOG2)

// Mask to extract MI offset within max MIB
#define MAX_MIB_MASK    (MAX_MIB_SIZE - 1)
#define MAX_MIB_MASK_2  (MAX_MIB_SIZE * 2 - 1)

// Maximum number of tile rows and tile columns
#if CONFIG_EXT_TILE
# define  MAX_TILE_ROWS 1024
# define  MAX_TILE_COLS 1024
#else
# define  MAX_TILE_ROWS 4
# define  MAX_TILE_COLS 64
#endif  // CONFIG_EXT_TILE

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

#define BLOCK_4X4       0
#define BLOCK_4X8       1
#define BLOCK_8X4       2
#define BLOCK_8X8       3
#define BLOCK_8X16      4
#define BLOCK_16X8      5
#define BLOCK_16X16     6
#define BLOCK_16X32     7
#define BLOCK_32X16     8
#define BLOCK_32X32     9
#define BLOCK_32X64    10
#define BLOCK_64X32    11
#define BLOCK_64X64    12
#if !CONFIG_EXT_PARTITION
# define BLOCK_SIZES   13
#else
# define BLOCK_64X128  13
# define BLOCK_128X64  14
# define BLOCK_128X128 15
# define BLOCK_SIZES   16
#endif  // !CONFIG_EXT_PARTITION
#define BLOCK_INVALID BLOCK_SIZES
#define BLOCK_LARGEST (BLOCK_SIZES - 1)
typedef uint8_t BLOCK_SIZE;

#if CONFIG_EXT_PARTITION_TYPES
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
#endif  // CONFIG_EXT_PARTITION_TYPES

typedef char PARTITION_CONTEXT;
#define PARTITION_PLOFFSET   4  // number of probability models per block size
#if CONFIG_EXT_PARTITION
# define PARTITION_CONTEXTS  (5 * PARTITION_PLOFFSET)
#else
# define PARTITION_CONTEXTS  (4 * PARTITION_PLOFFSET)
#endif  // CONFIG_EXT_PARTITION

// block transform size
typedef uint8_t TX_SIZE;
#define TX_4X4   ((TX_SIZE)0)   // 4x4 transform
#define TX_8X8   ((TX_SIZE)1)   // 8x8 transform
#define TX_16X16 ((TX_SIZE)2)   // 16x16 transform
#define TX_32X32 ((TX_SIZE)3)   // 32x32 transform
#define TX_SIZES ((TX_SIZE)4)

#if CONFIG_EXT_TX
#define TX_4X8   ((TX_SIZE)4)      // 4x8 transform
#define TX_8X4   ((TX_SIZE)5)      // 8x4 transform
#define TX_SIZES_ALL ((TX_SIZE)6)  // Includes rectangular transforms
#else
#define TX_SIZES_ALL ((TX_SIZE)4)
#endif  // CONFIG_EXT_TX

#define MAX_TX_SIZE_LOG2  5
#define MAX_TX_SIZE       (1 << MAX_TX_SIZE_LOG2)
#define MIN_TX_SIZE_LOG2  2
#define MIN_TX_SIZE       (1 << MIN_TX_SIZE_LOG2)
#define MAX_TX_SQUARE     (MAX_TX_SIZE * MAX_TX_SIZE)

// Number of maxium size transform blocks in the maximum size superblock
#define MAX_TX_BLOCKS_IN_MAX_SB_LOG2 \
  ((MAX_SB_SIZE_LOG2 - MAX_TX_SIZE_LOG2) * 2)
#define MAX_TX_BLOCKS_IN_MAX_SB (1 << MAX_TX_BLOCKS_IN_MAX_SB_LOG2)

#define MAX_NUM_TXB  (1 << (MAX_SB_SIZE_LOG2 - MIN_TX_SIZE_LOG2))

// frame transform mode
typedef enum {
  ONLY_4X4            = 0,        // only 4x4 transform used
  ALLOW_8X8           = 1,        // allow block transform size up to 8x8
  ALLOW_16X16         = 2,        // allow block transform size up to 16x16
  ALLOW_32X32         = 3,        // allow block transform size up to 32x32
  TX_MODE_SELECT      = 4,        // transform specified for each block
  TX_MODES            = 5,
} TX_MODE;

// 1D tx types
typedef enum {
  DCT_1D = 0,
  ADST_1D = 1,
  FLIPADST_1D = 2,
  IDTX_1D = 3,
  TX_TYPES_1D = 4,
} TX_TYPE_1D;

typedef enum {
  DCT_DCT   = 0,                  // DCT  in both horizontal and vertical
  ADST_DCT  = 1,                  // ADST in vertical, DCT in horizontal
  DCT_ADST  = 2,                  // DCT  in vertical, ADST in horizontal
  ADST_ADST = 3,                  // ADST in both directions
#if CONFIG_EXT_TX
  FLIPADST_DCT = 4,
  DCT_FLIPADST = 5,
  FLIPADST_FLIPADST = 6,
  ADST_FLIPADST = 7,
  FLIPADST_ADST = 8,
  IDTX = 9,
  V_DCT = 10,
  H_DCT = 11,
  V_ADST = 12,
  H_ADST = 13,
  V_FLIPADST = 14,
  H_FLIPADST = 15,
#endif  // CONFIG_EXT_TX
  TX_TYPES,
} TX_TYPE;

#if CONFIG_EXT_TX
#define EXT_TX_SIZES       4  // number of sizes that use extended transforms
#define EXT_TX_SETS_INTER  4  // Sets of transform selections for INTER
#define EXT_TX_SETS_INTRA  3  // Sets of transform selections for INTRA
#else
#define EXT_TX_SIZES       3  // number of sizes that use extended transforms
#endif  // CONFIG_EXT_TX

typedef enum {
  VP9_LAST_FLAG = 1 << 0,
#if CONFIG_EXT_REFS
  VP9_LAST2_FLAG = 1 << 1,
  VP9_LAST3_FLAG = 1 << 2,
  VP9_GOLD_FLAG = 1 << 3,
  VP9_BWD_FLAG = 1 << 4,
  VP9_ALT_FLAG = 1 << 5,
  VP9_REFFRAME_ALL = (1 << 6) - 1
#else
  VP9_GOLD_FLAG = 1 << 1,
  VP9_ALT_FLAG = 1 << 2,
  VP9_REFFRAME_ALL = (1 << 3) - 1
#endif  // CONFIG_EXT_REFS
} VP9_REFFRAME;

typedef enum {
  PLANE_TYPE_Y  = 0,
  PLANE_TYPE_UV = 1,
  PLANE_TYPES
} PLANE_TYPE;

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

#define DC_PRED    0       // Average of above and left pixels
#define V_PRED     1       // Vertical
#define H_PRED     2       // Horizontal
#define D45_PRED   3       // Directional 45  deg = round(arctan(1/1) * 180/pi)
#define D135_PRED  4       // Directional 135 deg = 180 - 45
#define D117_PRED  5       // Directional 117 deg = 180 - 63
#define D153_PRED  6       // Directional 153 deg = 180 - 27
#define D207_PRED  7       // Directional 207 deg = 180 + 27
#define D63_PRED   8       // Directional 63  deg = round(arctan(2/1) * 180/pi)
#define TM_PRED    9       // True-motion
#define NEARESTMV 10
#define NEARMV    11
#define ZEROMV    12
#define NEWMV     13
#if CONFIG_EXT_INTER
#define NEWFROMNEARMV     14
#define NEAREST_NEARESTMV 15
#define NEAREST_NEARMV    16
#define NEAR_NEARESTMV    17
#define NEAR_NEARMV       18
#define NEAREST_NEWMV     19
#define NEW_NEARESTMV     20
#define NEAR_NEWMV        21
#define NEW_NEARMV        22
#define ZERO_ZEROMV       23
#define NEW_NEWMV         24
#define MB_MODE_COUNT     25
#else
#define MB_MODE_COUNT     14
#endif  // CONFIG_EXT_INTER
typedef uint8_t PREDICTION_MODE;

#define INTRA_MODES (TM_PRED + 1)

typedef enum {
  SIMPLE_TRANSLATION = 0,
#if CONFIG_OBMC
  OBMC_CAUSAL,    // 2-sided OBMC
#endif  // CONFIG_OBMC
#if CONFIG_WARPED_MOTION
  WARPED_CAUSAL,  // 2-sided WARPED
#endif  // CONFIG_WARPED_MOTION
  MOTION_VARIATIONS
} MOTION_VARIATION;

#if CONFIG_EXT_INTER
typedef enum {
  II_DC_PRED = 0,
  II_V_PRED,
  II_H_PRED,
  II_D45_PRED,
  II_D135_PRED,
  II_D117_PRED,
  II_D153_PRED,
  II_D207_PRED,
  II_D63_PRED,
  II_TM_PRED,
  INTERINTRA_MODES
} INTERINTRA_MODE;

#endif  // CONFIG_EXT_INTER

#if CONFIG_EXT_INTRA
typedef enum {
  FILTER_DC_PRED,
  FILTER_V_PRED,
  FILTER_H_PRED,
  FILTER_D45_PRED,
  FILTER_D135_PRED,
  FILTER_D117_PRED,
  FILTER_D153_PRED,
  FILTER_D207_PRED,
  FILTER_D63_PRED,
  FILTER_TM_PRED,
  EXT_INTRA_MODES,
} EXT_INTRA_MODE;

#define FILTER_INTRA_MODES (FILTER_TM_PRED + 1)
#define DIRECTIONAL_MODES (INTRA_MODES - 2)
#endif  // CONFIG_EXT_INTRA

#if CONFIG_EXT_INTER
#define INTER_MODES (1 + NEWFROMNEARMV - NEARESTMV)
#else
#define INTER_MODES (1 + NEWMV - NEARESTMV)
#endif  // CONFIG_EXT_INTER

#if CONFIG_EXT_INTER
#define INTER_COMPOUND_MODES (1 + NEW_NEWMV - NEAREST_NEARESTMV)
#endif  // CONFIG_EXT_INTER

#define SKIP_CONTEXTS 3

#if CONFIG_REF_MV
#define NMV_CONTEXTS 3

#define NEWMV_MODE_CONTEXTS  7
#define ZEROMV_MODE_CONTEXTS 2
#define REFMV_MODE_CONTEXTS  9
#define DRL_MODE_CONTEXTS    5

#define ZEROMV_OFFSET 3
#define REFMV_OFFSET  4

#define NEWMV_CTX_MASK ((1 << ZEROMV_OFFSET) - 1)
#define ZEROMV_CTX_MASK ((1 << (REFMV_OFFSET - ZEROMV_OFFSET)) - 1)
#define REFMV_CTX_MASK ((1 << (8 - REFMV_OFFSET)) - 1)

#define ALL_ZERO_FLAG_OFFSET   8
#define SKIP_NEARESTMV_OFFSET  9
#define SKIP_NEARMV_OFFSET    10
#define SKIP_NEARESTMV_SUB8X8_OFFSET 11
#endif

#define INTER_MODE_CONTEXTS 7

/* Segment Feature Masks */
#define MAX_MV_REF_CANDIDATES 2

#if CONFIG_REF_MV
#define MAX_REF_MV_STACK_SIZE 16
#if CONFIG_EXT_PARTITION
#define REF_CAT_LEVEL 640
#else
#define REF_CAT_LEVEL 160
#endif  // CONFIG_EXT_PARTITION
#endif  // CONFIG_REF_MV

#define INTRA_INTER_CONTEXTS 4
#define COMP_INTER_CONTEXTS 5
#define REF_CONTEXTS 5

#if CONFIG_VAR_TX
#define TXFM_PARTITION_CONTEXTS 9
typedef TX_SIZE TXFM_CONTEXT;
#endif

#define NONE           -1
#define INTRA_FRAME     0
#define LAST_FRAME      1

#if CONFIG_EXT_REFS

#define LAST2_FRAME     2
#define LAST3_FRAME     3
#define GOLDEN_FRAME    4
#define BWDREF_FRAME    5
#define ALTREF_FRAME    6
#define MAX_REF_FRAMES  7
#define LAST_REF_FRAMES (LAST3_FRAME - LAST_FRAME + 1)

#else

#define GOLDEN_FRAME    2
#define ALTREF_FRAME    3
#define MAX_REF_FRAMES  4
#endif  // CONFIG_EXT_REFS

#define FWD_REFS (GOLDEN_FRAME - LAST_FRAME + 1)
#define FWD_RF_OFFSET(ref) (ref - LAST_FRAME)
#if CONFIG_EXT_REFS
#define BWD_REFS (ALTREF_FRAME - BWDREF_FRAME + 1)
#define BWD_RF_OFFSET(ref) (ref - BWDREF_FRAME)
#else
#define BWD_REFS 1
#define BWD_RF_OFFSET(ref) (ref - ALTREF_FRAME)
#endif

#define SINGLE_REFS (FWD_REFS + BWD_REFS)
#define COMP_REFS   (FWD_REFS * BWD_REFS)

#if CONFIG_REF_MV
#define MODE_CTX_REF_FRAMES (MAX_REF_FRAMES + COMP_REFS)
#else
#define MODE_CTX_REF_FRAMES MAX_REF_FRAMES
#endif

#if CONFIG_SUPERTX
#define PARTITION_SUPERTX_CONTEXTS 2
#define MAX_SUPERTX_BLOCK_SIZE BLOCK_32X32
#endif  // CONFIG_SUPERTX

#ifdef __cplusplus
}  // extern "C"
#endif

#endif  // VP10_COMMON_ENUMS_H_
