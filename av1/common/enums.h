/*
 * Copyright (c) 2016, Alliance for Open Media. All rights reserved
 *
 * This source code is subject to the terms of the BSD 2 Clause License and
 * the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
 * was not distributed with this source code in the LICENSE file, you can
 * obtain it at www.aomedia.org/license/software. If the Alliance for Open
 * Media Patent License 1.0 was not distributed with this source code in the
 * PATENTS file, you can obtain it at www.aomedia.org/license/patent.
 */

#ifndef AV1_COMMON_ENUMS_H_
#define AV1_COMMON_ENUMS_H_

#include "./aom_config.h"
#include "aom/aom_codec.h"
#include "aom/aom_integer.h"

#ifdef __cplusplus
extern "C" {
#endif

#undef MAX_SB_SIZE

// Max superblock size
#if CONFIG_EXT_PARTITION
#define MAX_SB_SIZE_LOG2 7
#else
#define MAX_SB_SIZE_LOG2 6
#endif  // CONFIG_EXT_PARTITION
#define MAX_SB_SIZE (1 << MAX_SB_SIZE_LOG2)
#define MAX_SB_SQUARE (MAX_SB_SIZE * MAX_SB_SIZE)

// Min superblock size
#define MIN_SB_SIZE_LOG2 6

// Pixels per Mode Info (MI) unit
#define MI_SIZE_LOG2 3
#define MI_SIZE (1 << MI_SIZE_LOG2)

// MI-units per max superblock (MI Block - MIB)
#define MAX_MIB_SIZE_LOG2 (MAX_SB_SIZE_LOG2 - MI_SIZE_LOG2)
#define MAX_MIB_SIZE (1 << MAX_MIB_SIZE_LOG2)

// MI-units per min superblock
#define MIN_MIB_SIZE_LOG2 (MIN_SB_SIZE_LOG2 - MI_SIZE_LOG2)

// Mask to extract MI offset within max MIB
#define MAX_MIB_MASK (MAX_MIB_SIZE - 1)
#define MAX_MIB_MASK_2 (MAX_MIB_SIZE * 2 - 1)

// Maximum number of tile rows and tile columns
#if CONFIG_EXT_TILE
#define MAX_TILE_ROWS 1024
#define MAX_TILE_COLS 1024
#else
#define MAX_TILE_ROWS 4
#define MAX_TILE_COLS 64
#endif  // CONFIG_EXT_TILE

#if CONFIG_VAR_TX
#define MAX_VARTX_DEPTH 2
#endif

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

// Note: Some enums use the attribute 'packed' to use smallest possible integer
// type, so that we can save memory when they are used in structs/arrays.

typedef enum ATTRIBUTE_PACKED {
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
#if CONFIG_EXT_PARTITION
  BLOCK_64X128,
  BLOCK_128X64,
  BLOCK_128X128,
#endif  // CONFIG_EXT_PARTITION

  BLOCK_SIZES,
  BLOCK_INVALID = BLOCK_SIZES,
  BLOCK_LARGEST = (BLOCK_SIZES - 1)
} BLOCK_SIZE;

typedef enum {
  PARTITION_NONE,
  PARTITION_HORZ,
  PARTITION_VERT,
  PARTITION_SPLIT,
#if CONFIG_EXT_PARTITION_TYPES
  PARTITION_HORZ_A,  // HORZ split and the left partition is split again
  PARTITION_HORZ_B,  // HORZ split and the right partition is split again
  PARTITION_VERT_A,  // VERT split and the top partition is split again
  PARTITION_VERT_B,  // VERT split and the bottom partition is split again
  EXT_PARTITION_TYPES,
#endif  // CONFIG_EXT_PARTITION_TYPES
  PARTITION_TYPES = PARTITION_SPLIT + 1,
  PARTITION_INVALID = 255
} PARTITION_TYPE;

typedef char PARTITION_CONTEXT;
#define PARTITION_PLOFFSET 4  // number of probability models per block size
#if CONFIG_EXT_PARTITION
#define PARTITION_CONTEXTS (5 * PARTITION_PLOFFSET)
#else
#define PARTITION_CONTEXTS (4 * PARTITION_PLOFFSET)
#endif  // CONFIG_EXT_PARTITION

// block transform size
typedef enum ATTRIBUTE_PACKED {
#if CONFIG_CB4X4
  TX_2X2,  // 2x2 transform
#endif
  TX_4X4,    // 4x4 transform
  TX_8X8,    // 8x8 transform
  TX_16X16,  // 16x16 transform
  TX_32X32,  // 32x32 transform
#if CONFIG_TX64X64
  TX_64X64,           // 64x64 transform
#endif                // CONFIG_TX64X64
  TX_4X8,             // 4x8 transform
  TX_8X4,             // 8x4 transform
  TX_8X16,            // 8x16 transform
  TX_16X8,            // 16x8 transform
  TX_16X32,           // 16x32 transform
  TX_32X16,           // 32x16 transform
#if 0                 // CONFIG_TX64X64
  // TODO(debargha): To be enabled later
  TX_32X64,                 // 32x64 transform
  TX_64X32,                 // 64x32 transform
#endif                // CONFIG_TX64X64
  TX_SIZES_ALL,       // Includes rectangular transforms
  TX_SIZES = TX_4X8,  // Does NOT include rectangular transforms
  TX_INVALID = 255    // Invalid transform size
} TX_SIZE;

#define MAX_TX_DEPTH (TX_32X32 - TX_4X4)

#define MAX_TX_SIZE_LOG2 (5 + CONFIG_TX64X64)
#define MAX_TX_SIZE (1 << MAX_TX_SIZE_LOG2)
#define MIN_TX_SIZE_LOG2 2
#define MIN_TX_SIZE (1 << MIN_TX_SIZE_LOG2)
#define MAX_TX_SQUARE (MAX_TX_SIZE * MAX_TX_SIZE)

// Number of maxium size transform blocks in the maximum size superblock
#define MAX_TX_BLOCKS_IN_MAX_SB_LOG2 ((MAX_SB_SIZE_LOG2 - MAX_TX_SIZE_LOG2) * 2)
#define MAX_TX_BLOCKS_IN_MAX_SB (1 << MAX_TX_BLOCKS_IN_MAX_SB_LOG2)

#define MAX_NUM_TXB (1 << (MAX_SB_SIZE_LOG2 - MIN_TX_SIZE_LOG2))

// frame transform mode
typedef enum {
  ONLY_4X4 = 0,        // only 4x4 transform used
  ALLOW_8X8 = 1,       // allow block transform size up to 8x8
  ALLOW_16X16 = 2,     // allow block transform size up to 16x16
  ALLOW_32X32 = 3,     // allow block transform size up to 32x32
  TX_MODE_SELECT = 4,  // transform specified for each block
  TX_MODES = 5,
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
  DCT_DCT = 0,    // DCT  in both horizontal and vertical
  ADST_DCT = 1,   // ADST in vertical, DCT in horizontal
  DCT_ADST = 2,   // DCT  in vertical, ADST in horizontal
  ADST_ADST = 3,  // ADST in both directions
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
#define EXT_TX_SIZES 4       // number of sizes that use extended transforms
#define EXT_TX_SETS_INTER 4  // Sets of transform selections for INTER
#define EXT_TX_SETS_INTRA 3  // Sets of transform selections for INTRA
#else
#if CONFIG_CB4X4
#define EXT_TX_SIZES 4  // number of sizes that use extended transforms
#else
#define EXT_TX_SIZES 3  // number of sizes that use extended transforms
#endif
#endif  // CONFIG_EXT_TX

typedef enum {
  AOM_LAST_FLAG = 1 << 0,
#if CONFIG_EXT_REFS
  AOM_LAST2_FLAG = 1 << 1,
  AOM_LAST3_FLAG = 1 << 2,
  AOM_GOLD_FLAG = 1 << 3,
  AOM_BWD_FLAG = 1 << 4,
  AOM_ALT_FLAG = 1 << 5,
  AOM_REFFRAME_ALL = (1 << 6) - 1
#else
  AOM_GOLD_FLAG = 1 << 1,
  AOM_ALT_FLAG = 1 << 2,
  AOM_REFFRAME_ALL = (1 << 3) - 1
#endif  // CONFIG_EXT_REFS
} AOM_REFFRAME;

typedef enum { PLANE_TYPE_Y = 0, PLANE_TYPE_UV = 1, PLANE_TYPES } PLANE_TYPE;

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

#ifdef CONFIG_CLPF
#define CLPF_NOFLAG -1
typedef enum {
  CLPF_NOSIZE = 0,
  CLPF_32X32 = 1,
  CLPF_64X64 = 2,
  CLPF_128X128 = 3
} CLPF_BLOCK_SIZE;
#endif

typedef enum ATTRIBUTE_PACKED {
  DC_PRED,    // Average of above and left pixels
  V_PRED,     // Vertical
  H_PRED,     // Horizontal
  D45_PRED,   // Directional 45  deg = round(arctan(1/1) * 180/pi)
  D135_PRED,  // Directional 135 deg = 180 - 45
  D117_PRED,  // Directional 117 deg = 180 - 63
  D153_PRED,  // Directional 153 deg = 180 - 27
  D207_PRED,  // Directional 207 deg = 180 + 27
  D63_PRED,   // Directional 63  deg = round(arctan(2/1) * 180/pi)
  TM_PRED,    // True-motion
  NEARESTMV,
  NEARMV,
  ZEROMV,
  NEWMV,
#if CONFIG_EXT_INTER
  NEWFROMNEARMV,
  NEAREST_NEARESTMV,
  NEAREST_NEARMV,
  NEAR_NEARESTMV,
  NEAR_NEARMV,
  NEAREST_NEWMV,
  NEW_NEARESTMV,
  NEAR_NEWMV,
  NEW_NEARMV,
  ZERO_ZEROMV,
  NEW_NEWMV,
#endif  // CONFIG_EXT_INTER
  MB_MODE_COUNT,
  INTRA_MODES = TM_PRED + 1
} PREDICTION_MODE;

typedef enum {
  SIMPLE_TRANSLATION = 0,
#if CONFIG_MOTION_VAR
  OBMC_CAUSAL,  // 2-sided OBMC
#endif          // CONFIG_MOTION_VAR
#if CONFIG_WARPED_MOTION
  WARPED_CAUSAL,  // 2-sided WARPED
#endif            // CONFIG_WARPED_MOTION
  MOTION_MODES
} MOTION_MODE;

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

#if CONFIG_FILTER_INTRA
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
  FILTER_INTRA_MODES,
} FILTER_INTRA_MODE;
#endif  // CONFIG_FILTER_INTRA

#if CONFIG_EXT_INTRA
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

#define NEWMV_MODE_CONTEXTS 7
#define ZEROMV_MODE_CONTEXTS 2
#define REFMV_MODE_CONTEXTS 9
#define DRL_MODE_CONTEXTS 5

#define ZEROMV_OFFSET 3
#define REFMV_OFFSET 4

#define NEWMV_CTX_MASK ((1 << ZEROMV_OFFSET) - 1)
#define ZEROMV_CTX_MASK ((1 << (REFMV_OFFSET - ZEROMV_OFFSET)) - 1)
#define REFMV_CTX_MASK ((1 << (8 - REFMV_OFFSET)) - 1)

#define ALL_ZERO_FLAG_OFFSET 8
#define SKIP_NEARESTMV_OFFSET 9
#define SKIP_NEARMV_OFFSET 10
#define SKIP_NEARESTMV_SUB8X8_OFFSET 11
#endif

#define INTER_MODE_CONTEXTS 7
#if CONFIG_DELTA_Q
#define DELTA_Q_SMALL 3
#define DELTA_Q_CONTEXTS (DELTA_Q_SMALL)
#define DEFAULT_DELTA_Q_RES 4
#endif

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
#define TXFM_PARTITION_CONTEXTS 16
typedef uint8_t TXFM_CONTEXT;
#endif

#define NONE -1
#define INTRA_FRAME 0
#define LAST_FRAME 1

#if CONFIG_EXT_REFS
#define LAST2_FRAME 2
#define LAST3_FRAME 3
#define GOLDEN_FRAME 4
#define BWDREF_FRAME 5
#define ALTREF_FRAME 6
#define LAST_REF_FRAMES (LAST3_FRAME - LAST_FRAME + 1)
#else
#define GOLDEN_FRAME 2
#define ALTREF_FRAME 3
#endif  // CONFIG_EXT_REFS

#define INTER_REFS_PER_FRAME (ALTREF_FRAME - LAST_FRAME + 1)
#define TOTAL_REFS_PER_FRAME (ALTREF_FRAME - INTRA_FRAME + 1)

#define FWD_REFS (GOLDEN_FRAME - LAST_FRAME + 1)
#define FWD_RF_OFFSET(ref) (ref - LAST_FRAME)
#if CONFIG_EXT_REFS
#define BWD_REFS (ALTREF_FRAME - BWDREF_FRAME + 1)
#define BWD_RF_OFFSET(ref) (ref - BWDREF_FRAME)
#else
#define BWD_REFS 1
#define BWD_RF_OFFSET(ref) (ref - ALTREF_FRAME)
#endif  // CONFIG_EXT_REFS

#define SINGLE_REFS (FWD_REFS + BWD_REFS)
#define COMP_REFS (FWD_REFS * BWD_REFS)

#if CONFIG_REF_MV
#define MODE_CTX_REF_FRAMES (TOTAL_REFS_PER_FRAME + COMP_REFS)
#else
#define MODE_CTX_REF_FRAMES TOTAL_REFS_PER_FRAME
#endif

#if CONFIG_SUPERTX
#define PARTITION_SUPERTX_CONTEXTS 2
#define MAX_SUPERTX_BLOCK_SIZE BLOCK_32X32
#endif  // CONFIG_SUPERTX

#if CONFIG_LOOP_RESTORATION
typedef enum {
  RESTORE_NONE,
  RESTORE_BILATERAL,
  RESTORE_WIENER,
  RESTORE_SWITCHABLE,
  RESTORE_SWITCHABLE_TYPES = RESTORE_SWITCHABLE,
  RESTORE_TYPES,
} RestorationType;
#endif  // CONFIG_LOOP_RESTORATION
#ifdef __cplusplus
}  // extern "C"
#endif

#endif  // AV1_COMMON_ENUMS_H_
