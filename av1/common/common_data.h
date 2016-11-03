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

#ifndef AV1_COMMON_COMMON_DATA_H_
#define AV1_COMMON_COMMON_DATA_H_

#include "av1/common/enums.h"
#include "aom/aom_integer.h"
#include "aom_dsp/aom_dsp_common.h"

#ifdef __cplusplus
extern "C" {
#endif

#if CONFIG_EXT_PARTITION
#define IF_EXT_PARTITION(...) __VA_ARGS__
#else
#define IF_EXT_PARTITION(...)
#endif

// Log 2 conversion lookup tables for block width and height
static const uint8_t b_width_log2_lookup[BLOCK_SIZES] = {
  0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, IF_EXT_PARTITION(4, 5, 5)
};
static const uint8_t b_height_log2_lookup[BLOCK_SIZES] = {
  0, 1, 0, 1, 2, 1, 2, 3, 2, 3, 4, 3, 4, IF_EXT_PARTITION(5, 4, 5)
};
// Log 2 conversion lookup tables for modeinfo width and height
static const uint8_t mi_width_log2_lookup[BLOCK_SIZES] = {
  0, 0, 0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3, IF_EXT_PARTITION(3, 4, 4)
};
static const uint8_t mi_height_log2_lookup[BLOCK_SIZES] = {
  0, 0, 0, 0, 1, 0, 1, 2, 1, 2, 3, 2, 3, IF_EXT_PARTITION(4, 3, 4)
};

// Width/height lookup tables in units of various block sizes
static const uint8_t block_size_wide[BLOCK_SIZES] = {
  4, 4, 8, 8, 8, 16, 16, 16, 32, 32, 32, 64, 64, IF_EXT_PARTITION(64, 128, 128)
};
static const uint8_t block_size_high[BLOCK_SIZES] = {
  4, 8, 4, 8, 16, 8, 16, 32, 16, 32, 64, 32, 64, IF_EXT_PARTITION(128, 64, 128)
};
static const uint8_t num_4x4_blocks_wide_lookup[BLOCK_SIZES] = {
  1, 1, 2, 2, 2, 4, 4, 4, 8, 8, 8, 16, 16, IF_EXT_PARTITION(16, 32, 32)
};
static const uint8_t num_4x4_blocks_high_lookup[BLOCK_SIZES] = {
  1, 2, 1, 2, 4, 2, 4, 8, 4, 8, 16, 8, 16, IF_EXT_PARTITION(32, 16, 32)
};
static const uint8_t num_8x8_blocks_wide_lookup[BLOCK_SIZES] = {
  1, 1, 1, 1, 1, 2, 2, 2, 4, 4, 4, 8, 8, IF_EXT_PARTITION(8, 16, 16)
};
static const uint8_t num_8x8_blocks_high_lookup[BLOCK_SIZES] = {
  1, 1, 1, 1, 2, 1, 2, 4, 2, 4, 8, 4, 8, IF_EXT_PARTITION(16, 8, 16)
};
static const uint8_t num_16x16_blocks_wide_lookup[BLOCK_SIZES] = {
  1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 4, 4, IF_EXT_PARTITION(4, 8, 8)
};
static const uint8_t num_16x16_blocks_high_lookup[BLOCK_SIZES] = {
  1, 1, 1, 1, 1, 1, 1, 2, 1, 2, 4, 2, 4, IF_EXT_PARTITION(8, 4, 8)
};

// AOMMIN(3, AOMMIN(b_width_log2(bsize), b_height_log2(bsize)))
static const uint8_t size_group_lookup[BLOCK_SIZES] = {
  0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3, 3, IF_EXT_PARTITION(3, 3, 3)
};

static const uint8_t num_pels_log2_lookup[BLOCK_SIZES] = {
  4, 5, 5, 6, 7, 7, 8, 9, 9, 10, 11, 11, 12, IF_EXT_PARTITION(13, 13, 14)
};

/* clang-format off */
static const PARTITION_TYPE
  partition_lookup[MAX_SB_SIZE_LOG2 - 1][BLOCK_SIZES] = {
  {     // 4X4 ->
    //                                    4X4
                                          PARTITION_NONE,
    // 4X8,            8X4,               8X8
    PARTITION_INVALID, PARTITION_INVALID, PARTITION_INVALID,
    // 8X16,           16X8,              16X16
    PARTITION_INVALID, PARTITION_INVALID, PARTITION_INVALID,
    // 16X32,          32X16,             32X32
    PARTITION_INVALID, PARTITION_INVALID, PARTITION_INVALID,
    // 32X64,          64X32,             64X64
    PARTITION_INVALID, PARTITION_INVALID, PARTITION_INVALID,
#if CONFIG_EXT_PARTITION
    PARTITION_INVALID, PARTITION_INVALID, PARTITION_INVALID,
#endif  // CONFIG_EXT_PARTITION
  }, {  // 8X8 ->
    //                                    4X4
                                          PARTITION_SPLIT,
    // 4X8,            8X4,               8X8
    PARTITION_VERT,    PARTITION_HORZ,    PARTITION_NONE,
    // 8X16,           16X8,              16X16
    PARTITION_INVALID, PARTITION_INVALID, PARTITION_INVALID,
    // 16X32,          32X16,             32X32
    PARTITION_INVALID, PARTITION_INVALID, PARTITION_INVALID,
    // 32X64,          64X32,             64X64
    PARTITION_INVALID, PARTITION_INVALID, PARTITION_INVALID,
#if CONFIG_EXT_PARTITION
    // 64x128,         128x64,            128x128
    PARTITION_INVALID, PARTITION_INVALID, PARTITION_INVALID,
#endif  // CONFIG_EXT_PARTITION
  }, {  // 16X16 ->
    //                                    4X4
                                          PARTITION_SPLIT,
    // 4X8,            8X4,               8X8
    PARTITION_SPLIT,   PARTITION_SPLIT,   PARTITION_SPLIT,
    // 8X16,           16X8,              16X16
    PARTITION_VERT,    PARTITION_HORZ,    PARTITION_NONE,
    // 16X32,          32X16,             32X32
    PARTITION_INVALID, PARTITION_INVALID, PARTITION_INVALID,
    // 32X64,          64X32,             64X64
    PARTITION_INVALID, PARTITION_INVALID, PARTITION_INVALID,
#if CONFIG_EXT_PARTITION
    // 64x128,         128x64,            128x128
    PARTITION_INVALID, PARTITION_INVALID, PARTITION_INVALID,
#endif  // CONFIG_EXT_PARTITION
  }, {  // 32X32 ->
    //                                    4X4
                                          PARTITION_SPLIT,
    // 4X8,            8X4,               8X8
    PARTITION_SPLIT,   PARTITION_SPLIT,   PARTITION_SPLIT,
    // 8X16,           16X8,              16X16
    PARTITION_SPLIT,   PARTITION_SPLIT,   PARTITION_SPLIT,
    // 16X32,          32X16,             32X32
    PARTITION_VERT,    PARTITION_HORZ,    PARTITION_NONE,
    // 32X64,          64X32,             64X64
    PARTITION_INVALID, PARTITION_INVALID, PARTITION_INVALID,
#if CONFIG_EXT_PARTITION
    // 64x128,         128x64,            128x128
    PARTITION_INVALID, PARTITION_INVALID, PARTITION_INVALID,
#endif  // CONFIG_EXT_PARTITION
  }, {  // 64X64 ->
    //                                    4X4
                                          PARTITION_SPLIT,
    // 4X8,            8X4,               8X8
    PARTITION_SPLIT,   PARTITION_SPLIT,   PARTITION_SPLIT,
    // 8X16,           16X8,              16X16
    PARTITION_SPLIT,   PARTITION_SPLIT,   PARTITION_SPLIT,
    // 16X32,          32X16,             32X32
    PARTITION_SPLIT,   PARTITION_SPLIT,   PARTITION_SPLIT,
    // 32X64,          64X32,             64X64
    PARTITION_VERT,    PARTITION_HORZ,    PARTITION_NONE,
#if CONFIG_EXT_PARTITION
    // 64x128,         128x64,            128x128
    PARTITION_INVALID, PARTITION_INVALID, PARTITION_INVALID,
  }, {  // 128x128 ->
    //                                    4X4
                                          PARTITION_SPLIT,
    // 4X8,            8X4,               8X8
    PARTITION_SPLIT,   PARTITION_SPLIT,   PARTITION_SPLIT,
    // 8X16,           16X8,              16X16
    PARTITION_SPLIT,   PARTITION_SPLIT,   PARTITION_SPLIT,
    // 16X32,          32X16,             32X32
    PARTITION_SPLIT,   PARTITION_SPLIT,   PARTITION_SPLIT,
    // 32X64,          64X32,             64X64
    PARTITION_SPLIT,   PARTITION_SPLIT,   PARTITION_SPLIT,
    // 64x128,         128x64,            128x128
    PARTITION_VERT,    PARTITION_HORZ,    PARTITION_NONE,
#endif  // CONFIG_EXT_PARTITION
  }
};

#if CONFIG_EXT_PARTITION_TYPES
static const BLOCK_SIZE subsize_lookup[EXT_PARTITION_TYPES][BLOCK_SIZES] =
#else
static const BLOCK_SIZE subsize_lookup[PARTITION_TYPES][BLOCK_SIZES] =
#endif  // CONFIG_EXT_PARTITION_TYPES
{
  {     // PARTITION_NONE
    //                            4X4
                                  BLOCK_4X4,
    // 4X8,        8X4,           8X8
    BLOCK_4X8,     BLOCK_8X4,     BLOCK_8X8,
    // 8X16,       16X8,          16X16
    BLOCK_8X16,    BLOCK_16X8,    BLOCK_16X16,
    // 16X32,      32X16,         32X32
    BLOCK_16X32,   BLOCK_32X16,   BLOCK_32X32,
    // 32X64,      64X32,         64X64
    BLOCK_32X64,   BLOCK_64X32,   BLOCK_64X64,
#if CONFIG_EXT_PARTITION
    // 64x128,     128x64,        128x128
    BLOCK_64X128,  BLOCK_128X64,  BLOCK_128X128,
#endif  // CONFIG_EXT_PARTITION
  }, {  // PARTITION_HORZ
    //                            4X4
                                  BLOCK_INVALID,
    // 4X8,        8X4,           8X8
    BLOCK_INVALID, BLOCK_INVALID, BLOCK_8X4,
    // 8X16,       16X8,          16X16
    BLOCK_INVALID, BLOCK_INVALID, BLOCK_16X8,
    // 16X32,      32X16,         32X32
    BLOCK_INVALID, BLOCK_INVALID, BLOCK_32X16,
    // 32X64,      64X32,         64X64
    BLOCK_INVALID, BLOCK_INVALID, BLOCK_64X32,
#if CONFIG_EXT_PARTITION
    // 64x128,     128x64,        128x128
    BLOCK_INVALID, BLOCK_INVALID, BLOCK_128X64,
#endif  // CONFIG_EXT_PARTITION
  }, {  // PARTITION_VERT
    //                            4X4
                                  BLOCK_INVALID,
    // 4X8,        8X4,           8X8
    BLOCK_INVALID, BLOCK_INVALID, BLOCK_4X8,
    // 8X16,       16X8,          16X16
    BLOCK_INVALID, BLOCK_INVALID, BLOCK_8X16,
    // 16X32,      32X16,         32X32
    BLOCK_INVALID, BLOCK_INVALID, BLOCK_16X32,
    // 32X64,      64X32,         64X64
    BLOCK_INVALID, BLOCK_INVALID, BLOCK_32X64,
#if CONFIG_EXT_PARTITION
    // 64x128,     128x64,        128x128
    BLOCK_INVALID, BLOCK_INVALID, BLOCK_64X128,
#endif  // CONFIG_EXT_PARTITION
  }, {  // PARTITION_SPLIT
    //                            4X4
                                  BLOCK_INVALID,
    // 4X8,        8X4,           8X8
    BLOCK_INVALID, BLOCK_INVALID, BLOCK_4X4,
    // 8X16,       16X8,          16X16
    BLOCK_INVALID, BLOCK_INVALID, BLOCK_8X8,
    // 16X32,      32X16,         32X32
    BLOCK_INVALID, BLOCK_INVALID, BLOCK_16X16,
    // 32X64,      64X32,         64X64
    BLOCK_INVALID, BLOCK_INVALID, BLOCK_32X32,
#if CONFIG_EXT_PARTITION
    // 64x128,     128x64,        128x128
    BLOCK_INVALID, BLOCK_INVALID, BLOCK_64X64,
#endif  // CONFIG_EXT_PARTITION
#if CONFIG_EXT_PARTITION_TYPES
  }, {  // PARTITION_HORZ_A
    //                            4X4
                                  BLOCK_INVALID,
    // 4X8,        8X4,           8X8
    BLOCK_INVALID, BLOCK_INVALID, BLOCK_8X4,
    // 8X16,       16X8,          16X16
    BLOCK_INVALID, BLOCK_INVALID, BLOCK_16X8,
    // 16X32,      32X16,         32X32
    BLOCK_INVALID, BLOCK_INVALID, BLOCK_32X16,
    // 32X64,      64X32,         64X64
    BLOCK_INVALID, BLOCK_INVALID, BLOCK_64X32,
#if CONFIG_EXT_PARTITION
    // 64x128,     128x64,        128x128
    BLOCK_INVALID, BLOCK_INVALID, BLOCK_128X64,
#endif  // CONFIG_EXT_PARTITION
  }, {  // PARTITION_HORZ_B
    //                            4X4
                                  BLOCK_INVALID,
    // 4X8,        8X4,           8X8
    BLOCK_INVALID, BLOCK_INVALID, BLOCK_8X4,
    // 8X16,       16X8,          16X16
    BLOCK_INVALID, BLOCK_INVALID, BLOCK_16X8,
    // 16X32,      32X16,         32X32
    BLOCK_INVALID, BLOCK_INVALID, BLOCK_32X16,
    // 32X64,      64X32,         64X64
    BLOCK_INVALID, BLOCK_INVALID, BLOCK_64X32,
#if CONFIG_EXT_PARTITION
    // 64x128,     128x64,        128x128
    BLOCK_INVALID, BLOCK_INVALID, BLOCK_128X64,
#endif  // CONFIG_EXT_PARTITION
  }, {  // PARTITION_VERT_A
    //                            4X4
                                  BLOCK_INVALID,
    // 4X8,        8X4,           8X8
    BLOCK_INVALID, BLOCK_INVALID, BLOCK_4X8,
    // 8X16,       16X8,          16X16
    BLOCK_INVALID, BLOCK_INVALID, BLOCK_8X16,
    // 16X32,      32X16,         32X32
    BLOCK_INVALID, BLOCK_INVALID, BLOCK_16X32,
    // 32X64,      64X32,         64X64
    BLOCK_INVALID, BLOCK_INVALID, BLOCK_32X64,
#if CONFIG_EXT_PARTITION
    // 64x128,     128x64,        128x128
    BLOCK_INVALID, BLOCK_INVALID, BLOCK_64X128,
#endif  // CONFIG_EXT_PARTITION
  }, {  // PARTITION_VERT_B
    //                            4X4
                                  BLOCK_INVALID,
    // 4X8,        8X4,           8X8
    BLOCK_INVALID, BLOCK_INVALID, BLOCK_4X8,
    // 8X16,       16X8,          16X16
    BLOCK_INVALID, BLOCK_INVALID, BLOCK_8X16,
    // 16X32,      32X16,         32X32
    BLOCK_INVALID, BLOCK_INVALID, BLOCK_16X32,
    // 32X64,      64X32,         64X64
    BLOCK_INVALID, BLOCK_INVALID, BLOCK_32X64,
#if CONFIG_EXT_PARTITION
    // 64x128,     128x64,        128x128
    BLOCK_INVALID, BLOCK_INVALID, BLOCK_64X128,
#endif  // CONFIG_EXT_PARTITION
#endif  // CONFIG_EXT_PARTITION_TYPES
  }
};

static const TX_SIZE max_txsize_lookup[BLOCK_SIZES] = {
  //                   4X4
                       TX_4X4,
  // 4X8,    8X4,      8X8
  TX_4X4,    TX_4X4,   TX_8X8,
  // 8X16,   16X8,     16X16
  TX_8X8,    TX_8X8,   TX_16X16,
  // 16X32,  32X16,    32X32
  TX_16X16,  TX_16X16, TX_32X32,
  // 32X64,  64X32,    64X64
  TX_32X32,  TX_32X32, TX_32X32,
#if CONFIG_EXT_PARTITION
  // 64x128, 128x64,   128x128
  TX_32X32,  TX_32X32, TX_32X32,
#endif  // CONFIG_EXT_PARTITION
};

static const TX_SIZE max_txsize_rect_lookup[BLOCK_SIZES] = {
  //                   4X4
                       TX_4X4,
  // 4X8,    8X4,      8X8
  TX_4X8,    TX_8X4,   TX_8X8,
  // 8X16,   16X8,     16X16
  TX_8X16,   TX_16X8,  TX_16X16,
  // 16X32,  32X16,    32X32
  TX_16X32,  TX_32X16, TX_32X32,
  // 32X64,  64X32,    64X64
  TX_32X32,  TX_32X32, TX_32X32,
#if CONFIG_EXT_PARTITION
  // 64x128, 128x64,   128x128
  TX_32X32,  TX_32X32, TX_32X32,
#endif  // CONFIG_EXT_PARTITION
};

// Same as "max_txsize_lookup[bsize] - TX_8X8", invalid for bsize < 8X8
static const int32_t intra_tx_size_cat_lookup[BLOCK_SIZES] = {
  //                                      4X4
                                          INT32_MIN,
  // 4X8,             8X4,                8X8
  INT32_MIN,          INT32_MIN,          TX_8X8 - TX_8X8,
  // 8X16,            16X8,               16X16
  TX_8X8 - TX_8X8,    TX_8X8 - TX_8X8,    TX_16X16 - TX_8X8,
  // 16X32,           32X16,              32X32
  TX_16X16 - TX_8X8,  TX_16X16 - TX_8X8,  TX_32X32 - TX_8X8,
  // 32X64,           64X32,              64X64
  TX_32X32 - TX_8X8,  TX_32X32 - TX_8X8,  TX_32X32 - TX_8X8,
#if CONFIG_EXT_PARTITION
  // 64x128,          128x64,             128x128
  TX_32X32 - TX_8X8,  TX_32X32 - TX_8X8,  TX_32X32 - TX_8X8,
#endif  // CONFIG_EXT_PARTITION
};

#if CONFIG_EXT_TX && CONFIG_RECT_TX
// Same as "max_txsize_lookup[bsize] - TX_8X8", except for rectangular
// block which may use a rectangular transform, in which  case it is
// "(max_txsize_lookup[bsize] + 1) - TX_8X8", invalid for bsize < 8X8
static const int32_t inter_tx_size_cat_lookup[BLOCK_SIZES] = {
  //                                      4X4
                                          INT32_MIN,
  // 4X8,             8X4,                8X8
  INT32_MIN,          INT32_MIN,           TX_8X8 - TX_8X8,
  // 8X16,            16X8,               16X16
  TX_16X16 - TX_8X8,  TX_16X16 - TX_8X8,  TX_16X16 - TX_8X8,
  // 16X32,           32X16,              32X32
  TX_32X32 - TX_8X8,  TX_32X32 - TX_8X8,  TX_32X32 - TX_8X8,
  // 32X64,           64X32,              64X64
  TX_32X32 - TX_8X8,  TX_32X32 - TX_8X8,  TX_32X32 - TX_8X8,
#if CONFIG_EXT_PARTITION
  // 64x128,          128x64,             128x128
  TX_32X32 - TX_8X8,  TX_32X32 - TX_8X8,  TX_32X32 - TX_8X8,
#endif  // CONFIG_EXT_PARTITION
};
#else
#define inter_tx_size_cat_lookup intra_tx_size_cat_lookup
#endif  // CONFIG_EXT_TX && CONFIG_RECT_TX

/* clang-format on */

static const TX_SIZE sub_tx_size_map[TX_SIZES_ALL] = {
#if CONFIG_CB4X4
  TX_2X2,  // TX_2X2
#endif
  TX_4X4,    // TX_4X4
  TX_4X4,    // TX_8X8
  TX_8X8,    // TX_16X16
  TX_16X16,  // TX_32X32
  TX_4X4,    // TX_4X8
  TX_4X4,    // TX_8X4
  TX_8X8,    // TX_8X16
  TX_8X8,    // TX_16X8
  TX_16X16,  // TX_16X32
  TX_16X16   // TX_32X16
};

static const TX_SIZE txsize_horz_map[TX_SIZES_ALL] = {
#if CONFIG_CB4X4
  TX_2X2,  // TX_2X2
#endif
  TX_4X4,    // TX_4X4
  TX_8X8,    // TX_8X8
  TX_16X16,  // TX_16X16
  TX_32X32,  // TX_32X32
  TX_4X4,    // TX_4X8
  TX_8X8,    // TX_8X4
  TX_8X8,    // TX_8X16
  TX_16X16,  // TX_16X8
  TX_16X16,  // TX_16X32
  TX_32X32   // TX_32X16
};

static const TX_SIZE txsize_vert_map[TX_SIZES_ALL] = {
#if CONFIG_CB4X4
  TX_2X2,  // TX_2X2
#endif
  TX_4X4,    // TX_4X4
  TX_8X8,    // TX_8X8
  TX_16X16,  // TX_16X16
  TX_32X32,  // TX_32X32
  TX_8X8,    // TX_4X8
  TX_4X4,    // TX_8X4
  TX_16X16,  // TX_8X16
  TX_8X8,    // TX_16X8
  TX_32X32,  // TX_16X32
  TX_16X16   // TX_32X16
};

// Transform block width in pixels
static const int tx_size_wide[TX_SIZES_ALL] = {
#if CONFIG_CB4X4
  2,
#endif
  4, 8, 16, 32, 4, 8, 8, 16, 16, 32,
};

// Transform block height in pixels
static const int tx_size_high[TX_SIZES_ALL] = {
#if CONFIG_CB4X4
  2,
#endif
  4, 8, 16, 32, 8, 4, 16, 8, 32, 16,
};

// Transform block width in unit
static const int tx_size_wide_unit[TX_SIZES_ALL] = {
#if CONFIG_CB4X4
  1,
#endif
  1, 2, 4, 8, 1, 2, 2, 4, 4, 8,
};

// Transform block height in unit
static const int tx_size_high_unit[TX_SIZES_ALL] = {
#if CONFIG_CB4X4
  1,
#endif
  1, 2, 4, 8, 2, 1, 4, 2, 8, 4,
};

// Transform block width in log2
static const int tx_size_wide_log2[TX_SIZES_ALL] = {
#if CONFIG_CB4X4
  2,
#endif
  2, 3, 4, 5, 2, 3, 3, 4, 4, 5,
};

// Transform block height in log2
static const int tx_size_high_log2[TX_SIZES_ALL] = {
#if CONFIG_CB4X4
  2,
#endif
  2, 3, 4, 5, 3, 2, 4, 3, 5, 4,
};

static const int tx_size_2d[TX_SIZES_ALL] = {
#if CONFIG_CB4X4
  4,
#endif
  16, 64, 256, 1024, 32, 32, 128, 128, 512, 512,
};

static const uint8_t tx_size_1d_log2[TX_SIZES] = { 2, 3, 4, 5 };

// TODO(jingning): Temporary table during the construction.
static const int tx_size_1d_in_unit_log2[TX_SIZES] = { 0, 1, 2, 3 };

static const BLOCK_SIZE txsize_to_bsize[TX_SIZES_ALL] = {
#if CONFIG_CB4X4
  BLOCK_4X4,  // TX_2X2
#endif
  BLOCK_4X4,    // TX_4X4
  BLOCK_8X8,    // TX_8X8
  BLOCK_16X16,  // TX_16X16
  BLOCK_32X32,  // TX_32X32
  BLOCK_4X8,    // TX_4X8
  BLOCK_8X4,    // TX_8X4
  BLOCK_8X16,   // TX_8X16
  BLOCK_16X8,   // TX_16X8
  BLOCK_16X32,  // TX_16X32
  BLOCK_32X16,  // TX_32X16
};

static const TX_SIZE txsize_sqr_map[TX_SIZES_ALL] = {
#if CONFIG_CB4X4
  TX_2X2,  // TX_2X2
#endif
  TX_4X4,    // TX_4X4
  TX_8X8,    // TX_8X8
  TX_16X16,  // TX_16X16
  TX_32X32,  // TX_32X32
  TX_4X4,    // TX_4X8
  TX_4X4,    // TX_8X4
  TX_8X8,    // TX_8X16
  TX_8X8,    // TX_16X8
  TX_16X16,  // TX_16X32
  TX_16X16,  // TX_32X16
};

static const TX_SIZE txsize_sqr_up_map[TX_SIZES_ALL] = {
#if CONFIG_CB4X4
  TX_2X2,  // TX_2X2
#endif
  TX_4X4,    // TX_4X4
  TX_8X8,    // TX_8X8
  TX_16X16,  // TX_16X16
  TX_32X32,  // TX_32X32
  TX_8X8,    // TX_4X8
  TX_8X8,    // TX_8X4
  TX_16X16,  // TX_8X16
  TX_16X16,  // TX_16X8
  TX_32X32,  // TX_16X32
  TX_32X32,  // TX_32X16
};

static const TX_SIZE tx_mode_to_biggest_tx_size[TX_MODES] = {
  TX_4X4,    // ONLY_4X4
  TX_8X8,    // ALLOW_8X8
  TX_16X16,  // ALLOW_16X16
  TX_32X32,  // ALLOW_32X32
  TX_32X32,  // TX_MODE_SELECT
};

static const BLOCK_SIZE ss_size_lookup[BLOCK_SIZES][2][2] = {
  //  ss_x == 0    ss_x == 0        ss_x == 1      ss_x == 1
  //  ss_y == 0    ss_y == 1        ss_y == 0      ss_y == 1
  { { BLOCK_4X4, BLOCK_INVALID }, { BLOCK_INVALID, BLOCK_INVALID } },
  { { BLOCK_4X8, BLOCK_4X4 }, { BLOCK_INVALID, BLOCK_INVALID } },
  { { BLOCK_8X4, BLOCK_INVALID }, { BLOCK_4X4, BLOCK_INVALID } },
  { { BLOCK_8X8, BLOCK_8X4 }, { BLOCK_4X8, BLOCK_4X4 } },
  { { BLOCK_8X16, BLOCK_8X8 }, { BLOCK_INVALID, BLOCK_4X8 } },
  { { BLOCK_16X8, BLOCK_INVALID }, { BLOCK_8X8, BLOCK_8X4 } },
  { { BLOCK_16X16, BLOCK_16X8 }, { BLOCK_8X16, BLOCK_8X8 } },
  { { BLOCK_16X32, BLOCK_16X16 }, { BLOCK_INVALID, BLOCK_8X16 } },
  { { BLOCK_32X16, BLOCK_INVALID }, { BLOCK_16X16, BLOCK_16X8 } },
  { { BLOCK_32X32, BLOCK_32X16 }, { BLOCK_16X32, BLOCK_16X16 } },
  { { BLOCK_32X64, BLOCK_32X32 }, { BLOCK_INVALID, BLOCK_16X32 } },
  { { BLOCK_64X32, BLOCK_INVALID }, { BLOCK_32X32, BLOCK_32X16 } },
  { { BLOCK_64X64, BLOCK_64X32 }, { BLOCK_32X64, BLOCK_32X32 } },
#if CONFIG_EXT_PARTITION
  { { BLOCK_64X128, BLOCK_64X64 }, { BLOCK_INVALID, BLOCK_32X64 } },
  { { BLOCK_128X64, BLOCK_INVALID }, { BLOCK_64X64, BLOCK_64X32 } },
  { { BLOCK_128X128, BLOCK_128X64 }, { BLOCK_64X128, BLOCK_64X64 } },
#endif  // CONFIG_EXT_PARTITION
};

#define USE_UV_RECT_TX 1
static const TX_SIZE uv_txsize_lookup[BLOCK_SIZES][TX_SIZES_ALL][2][2] = {
  //  ss_x == 0    ss_x == 0        ss_x == 1      ss_x == 1
  //  ss_y == 0    ss_y == 1        ss_y == 0      ss_y == 1
  {
// BLOCK_4X4
#if CONFIG_CB4X4
      { { TX_4X4, TX_4X4 }, { TX_4X4, TX_4X4 } },
#endif
      { { TX_4X4, TX_4X4 }, { TX_4X4, TX_4X4 } },
      { { TX_4X4, TX_4X4 }, { TX_4X4, TX_4X4 } },
      { { TX_4X4, TX_4X4 }, { TX_4X4, TX_4X4 } },
      { { TX_4X4, TX_4X4 }, { TX_4X4, TX_4X4 } },
      { { TX_4X4, TX_4X4 }, { TX_4X4, TX_4X4 } },
      { { TX_4X4, TX_4X4 }, { TX_4X4, TX_4X4 } },
      { { TX_4X4, TX_4X4 }, { TX_4X4, TX_4X4 } },
      { { TX_4X4, TX_4X4 }, { TX_4X4, TX_4X4 } },
      { { TX_4X4, TX_4X4 }, { TX_4X4, TX_4X4 } },
      { { TX_4X4, TX_4X4 }, { TX_4X4, TX_4X4 } },
  },
  {
// BLOCK_4X8
#if CONFIG_CB4X4
      { { TX_4X4, TX_4X4 }, { TX_4X4, TX_4X4 } },
#endif
      { { TX_4X4, TX_4X4 }, { TX_4X4, TX_4X4 } },
      { { TX_4X4, TX_4X4 }, { TX_4X4, TX_4X4 } },
      { { TX_4X4, TX_4X4 }, { TX_4X4, TX_4X4 } },
      { { TX_4X4, TX_4X4 }, { TX_4X4, TX_4X4 } },
#if CONFIG_RECT_TX && USE_UV_RECT_TX
      { { TX_4X8, TX_4X4 }, { TX_4X4, TX_4X4 } },  // used
#else
      { { TX_4X4, TX_4X4 }, { TX_4X4, TX_4X4 } },      // used
#endif  // CONFIG_RECT_TX && USE_UV_RECT_TX
      { { TX_INVALID, TX_INVALID }, { TX_INVALID, TX_INVALID } },
      { { TX_INVALID, TX_INVALID }, { TX_INVALID, TX_INVALID } },
      { { TX_INVALID, TX_INVALID }, { TX_INVALID, TX_INVALID } },
      { { TX_INVALID, TX_INVALID }, { TX_INVALID, TX_INVALID } },
      { { TX_INVALID, TX_INVALID }, { TX_INVALID, TX_INVALID } },
  },
  {
// BLOCK_8X4
#if CONFIG_CB4X4
      { { TX_4X4, TX_4X4 }, { TX_4X4, TX_4X4 } },
#endif
      { { TX_4X4, TX_4X4 }, { TX_4X4, TX_4X4 } },
      { { TX_4X4, TX_4X4 }, { TX_4X4, TX_4X4 } },
      { { TX_4X4, TX_4X4 }, { TX_4X4, TX_4X4 } },
      { { TX_4X4, TX_4X4 }, { TX_4X4, TX_4X4 } },
      { { TX_INVALID, TX_INVALID }, { TX_INVALID, TX_INVALID } },
#if CONFIG_RECT_TX && USE_UV_RECT_TX
      { { TX_8X4, TX_4X4 }, { TX_4X4, TX_4X4 } },  // used
#else
      { { TX_4X4, TX_4X4 }, { TX_4X4, TX_4X4 } },      // used
#endif  // CONFIG_RECT_TX && USE_UV_RECT_TX
      { { TX_INVALID, TX_INVALID }, { TX_INVALID, TX_INVALID } },
      { { TX_INVALID, TX_INVALID }, { TX_INVALID, TX_INVALID } },
      { { TX_INVALID, TX_INVALID }, { TX_INVALID, TX_INVALID } },
      { { TX_INVALID, TX_INVALID }, { TX_INVALID, TX_INVALID } },
  },
  {
// BLOCK_8X8
#if CONFIG_CB4X4
      { { TX_4X4, TX_4X4 }, { TX_4X4, TX_4X4 } },
#endif
      { { TX_4X4, TX_4X4 }, { TX_4X4, TX_4X4 } },
      { { TX_8X8, TX_4X4 }, { TX_4X4, TX_4X4 } },
      { { TX_8X8, TX_8X8 }, { TX_4X4, TX_4X4 } },
      { { TX_8X8, TX_8X8 }, { TX_4X4, TX_4X4 } },
      { { TX_4X8, TX_4X4 }, { TX_4X8, TX_4X4 } },
      { { TX_8X4, TX_8X4 }, { TX_4X4, TX_4X4 } },
      { { TX_8X8, TX_8X4 }, { TX_4X8, TX_4X4 } },
      { { TX_8X8, TX_8X4 }, { TX_4X8, TX_4X4 } },
      { { TX_8X8, TX_8X4 }, { TX_4X8, TX_4X4 } },
      { { TX_8X8, TX_8X4 }, { TX_4X8, TX_4X4 } },
  },
  {
// BLOCK_8X16
#if CONFIG_CB4X4
      { { TX_4X4, TX_4X4 }, { TX_4X4, TX_4X4 } },
#endif
      { { TX_4X4, TX_4X4 }, { TX_4X4, TX_4X4 } },
      { { TX_8X8, TX_8X8 }, { TX_4X4, TX_4X4 } },
      { { TX_8X8, TX_8X8 }, { TX_4X4, TX_4X4 } },
      { { TX_8X8, TX_8X8 }, { TX_4X4, TX_4X4 } },
      { { TX_4X8, TX_4X8 }, { TX_4X8, TX_4X8 } },
      { { TX_8X4, TX_8X4 }, { TX_4X4, TX_4X4 } },
#if CONFIG_RECT_TX && USE_UV_RECT_TX
      { { TX_8X16, TX_8X8 }, { TX_4X8, TX_4X8 } },  // used
#else
      { { TX_8X8, TX_8X8 }, { TX_4X4, TX_4X4 } },      // used
#endif  // CONFIG_RECT_TX && USE_UV_RECT_TX
      { { TX_INVALID, TX_INVALID }, { TX_INVALID, TX_INVALID } },
      { { TX_INVALID, TX_INVALID }, { TX_INVALID, TX_INVALID } },
      { { TX_INVALID, TX_INVALID }, { TX_INVALID, TX_INVALID } },
  },
  {
// BLOCK_16X8
#if CONFIG_CB4X4
      { { TX_4X4, TX_4X4 }, { TX_4X4, TX_4X4 } },
#endif
      { { TX_4X4, TX_4X4 }, { TX_4X4, TX_4X4 } },
      { { TX_8X8, TX_4X4 }, { TX_8X8, TX_4X4 } },
      { { TX_8X8, TX_4X4 }, { TX_8X8, TX_8X8 } },
      { { TX_8X8, TX_4X4 }, { TX_8X8, TX_8X8 } },
      { { TX_4X8, TX_4X4 }, { TX_4X8, TX_4X4 } },
      { { TX_8X4, TX_8X4 }, { TX_8X4, TX_8X4 } },
      { { TX_INVALID, TX_INVALID }, { TX_INVALID, TX_INVALID } },
#if CONFIG_RECT_TX && USE_UV_RECT_TX
      { { TX_16X8, TX_8X4 }, { TX_8X8, TX_8X4 } },  // used
#else
      { { TX_8X8, TX_4X4 }, { TX_8X8, TX_4X4 } },      // used
#endif  // CONFIG_RECT_TX && USE_UV_RECT_TX
      { { TX_INVALID, TX_INVALID }, { TX_INVALID, TX_INVALID } },
      { { TX_INVALID, TX_INVALID }, { TX_INVALID, TX_INVALID } },
  },
  {
// BLOCK_16X16
#if CONFIG_CB4X4
      { { TX_4X4, TX_4X4 }, { TX_4X4, TX_4X4 } },
#endif
      { { TX_4X4, TX_4X4 }, { TX_4X4, TX_4X4 } },
      { { TX_8X8, TX_8X8 }, { TX_8X8, TX_8X8 } },
      { { TX_16X16, TX_8X8 }, { TX_8X8, TX_8X8 } },
      { { TX_16X16, TX_8X8 }, { TX_8X8, TX_8X8 } },
      { { TX_4X8, TX_4X8 }, { TX_4X8, TX_4X8 } },
      { { TX_8X4, TX_8X4 }, { TX_8X4, TX_8X4 } },
      { { TX_8X16, TX_8X8 }, { TX_8X16, TX_8X8 } },
      { { TX_16X8, TX_16X8 }, { TX_8X8, TX_8X8 } },
      { { TX_16X16, TX_16X8 }, { TX_8X16, TX_8X8 } },
      { { TX_16X16, TX_16X8 }, { TX_8X16, TX_8X8 } },
  },
  {
// BLOCK_16X32
#if CONFIG_CB4X4
      { { TX_4X4, TX_4X4 }, { TX_4X4, TX_4X4 } },
#endif
      { { TX_4X4, TX_4X4 }, { TX_4X4, TX_4X4 } },
      { { TX_8X8, TX_8X8 }, { TX_8X8, TX_8X8 } },
      { { TX_16X16, TX_16X16 }, { TX_8X8, TX_8X8 } },
      { { TX_16X16, TX_16X16 }, { TX_8X8, TX_8X8 } },
      { { TX_4X8, TX_4X8 }, { TX_4X8, TX_4X8 } },
      { { TX_8X4, TX_8X4 }, { TX_8X4, TX_8X4 } },
      { { TX_8X16, TX_8X16 }, { TX_8X16, TX_8X16 } },
      { { TX_16X8, TX_16X8 }, { TX_8X8, TX_8X8 } },
#if CONFIG_RECT_TX && USE_UV_RECT_TX
      { { TX_16X32, TX_16X16 }, { TX_8X16, TX_8X16 } },  // used
#else
      { { TX_16X16, TX_16X16 }, { TX_8X8, TX_8X8 } },  // used
#endif  // CONFIG_RECT_TX && USE_UV_RECT_TX
      { { TX_INVALID, TX_INVALID }, { TX_INVALID, TX_INVALID } },
  },
  {
// BLOCK_32X16
#if CONFIG_CB4X4
      { { TX_4X4, TX_4X4 }, { TX_4X4, TX_4X4 } },
#endif
      { { TX_4X4, TX_4X4 }, { TX_4X4, TX_4X4 } },
      { { TX_8X8, TX_8X8 }, { TX_8X8, TX_8X8 } },
      { { TX_16X16, TX_8X8 }, { TX_16X16, TX_8X8 } },
      { { TX_16X16, TX_8X8 }, { TX_16X16, TX_8X8 } },
      { { TX_4X8, TX_4X8 }, { TX_4X8, TX_4X8 } },
      { { TX_8X4, TX_8X4 }, { TX_8X4, TX_8X4 } },
      { { TX_8X16, TX_8X8 }, { TX_8X16, TX_8X8 } },
      { { TX_16X8, TX_16X8 }, { TX_16X8, TX_16X8 } },
      { { TX_INVALID, TX_INVALID }, { TX_INVALID, TX_INVALID } },
#if CONFIG_RECT_TX && USE_UV_RECT_TX
      { { TX_32X16, TX_16X8 }, { TX_16X16, TX_16X8 } },  // used
#else
      { { TX_16X16, TX_8X8 }, { TX_16X16, TX_8X8 } },  // used
#endif  // CONFIG_RECT_TX && USE_UV_RECT_TX
  },
  {
// BLOCK_32X32
#if CONFIG_CB4X4
      { { TX_4X4, TX_4X4 }, { TX_4X4, TX_4X4 } },
#endif
      { { TX_4X4, TX_4X4 }, { TX_4X4, TX_4X4 } },
      { { TX_8X8, TX_8X8 }, { TX_8X8, TX_8X8 } },
      { { TX_16X16, TX_16X16 }, { TX_16X16, TX_16X16 } },
      { { TX_32X32, TX_16X16 }, { TX_16X16, TX_16X16 } },
      { { TX_4X8, TX_4X8 }, { TX_4X8, TX_4X8 } },
      { { TX_8X4, TX_8X4 }, { TX_8X4, TX_8X4 } },
      { { TX_8X16, TX_8X16 }, { TX_8X16, TX_8X16 } },
      { { TX_16X8, TX_16X8 }, { TX_16X8, TX_16X8 } },
      { { TX_16X32, TX_16X16 }, { TX_16X32, TX_16X16 } },
      { { TX_32X16, TX_32X16 }, { TX_16X16, TX_16X16 } },
  },
  {
// BLOCK_32X64
#if CONFIG_CB4X4
      { { TX_4X4, TX_4X4 }, { TX_4X4, TX_4X4 } },
#endif
      { { TX_4X4, TX_4X4 }, { TX_4X4, TX_4X4 } },
      { { TX_8X8, TX_8X8 }, { TX_8X8, TX_8X8 } },
      { { TX_16X16, TX_16X16 }, { TX_16X16, TX_16X16 } },
      { { TX_32X32, TX_32X32 }, { TX_16X16, TX_16X16 } },
      { { TX_4X8, TX_4X8 }, { TX_4X8, TX_4X8 } },
      { { TX_8X4, TX_8X4 }, { TX_8X4, TX_8X4 } },
      { { TX_8X16, TX_8X16 }, { TX_8X16, TX_8X16 } },
      { { TX_16X8, TX_16X8 }, { TX_16X8, TX_16X8 } },
      { { TX_16X32, TX_16X32 }, { TX_16X16, TX_16X16 } },
      { { TX_32X16, TX_32X16 }, { TX_16X16, TX_16X16 } },
  },
  {
// BLOCK_64X32
#if CONFIG_CB4X4
      { { TX_4X4, TX_4X4 }, { TX_4X4, TX_4X4 } },
#endif
      { { TX_4X4, TX_4X4 }, { TX_4X4, TX_4X4 } },
      { { TX_8X8, TX_8X8 }, { TX_8X8, TX_8X8 } },
      { { TX_16X16, TX_16X16 }, { TX_16X16, TX_16X16 } },
      { { TX_32X32, TX_16X16 }, { TX_32X32, TX_16X16 } },
      { { TX_4X8, TX_4X8 }, { TX_4X8, TX_4X8 } },
      { { TX_8X4, TX_8X4 }, { TX_8X4, TX_8X4 } },
      { { TX_8X16, TX_8X16 }, { TX_8X16, TX_8X16 } },
      { { TX_16X8, TX_16X8 }, { TX_16X8, TX_16X8 } },
      { { TX_16X32, TX_16X16 }, { TX_16X32, TX_16X16 } },
      { { TX_32X16, TX_16X16 }, { TX_32X16, TX_16X16 } },
  },
  {
// BLOCK_64X64
#if CONFIG_CB4X4
      { { TX_4X4, TX_4X4 }, { TX_4X4, TX_4X4 } },
#endif
      { { TX_4X4, TX_4X4 }, { TX_4X4, TX_4X4 } },
      { { TX_8X8, TX_8X8 }, { TX_8X8, TX_8X8 } },
      { { TX_16X16, TX_16X16 }, { TX_16X16, TX_16X16 } },
      { { TX_32X32, TX_32X32 }, { TX_32X32, TX_32X32 } },
#if CONFIG_EXT_PARTITION
  },
  {
// BLOCK_64X128
#if CONFIG_CB4X4
      { { TX_4X4, TX_4X4 }, { TX_4X4, TX_4X4 } },
#endif
      { { TX_4X4, TX_4X4 }, { TX_4X4, TX_4X4 } },
      { { TX_8X8, TX_8X8 }, { TX_8X8, TX_8X8 } },
      { { TX_16X16, TX_16X16 }, { TX_16X16, TX_16X16 } },
      { { TX_32X32, TX_32X32 }, { TX_32X32, TX_32X32 } },
      { { TX_4X8, TX_4X8 }, { TX_4X8, TX_4X8 } },
      { { TX_8X4, TX_8X4 }, { TX_8X4, TX_8X4 } },
      { { TX_8X16, TX_8X16 }, { TX_8X16, TX_8X16 } },
      { { TX_16X8, TX_16X8 }, { TX_16X8, TX_16X8 } },
      { { TX_16X32, TX_16X32 }, { TX_16X32, TX_16X32 } },
      { { TX_32X16, TX_32X16 }, { TX_32X16, TX_32X16 } },
  },
  {
// BLOCK_128X64
#if CONFIG_CB4X4
      { { TX_4X4, TX_4X4 }, { TX_4X4, TX_4X4 } },
#endif
      { { TX_4X4, TX_4X4 }, { TX_4X4, TX_4X4 } },
      { { TX_8X8, TX_8X8 }, { TX_8X8, TX_8X8 } },
      { { TX_16X16, TX_16X16 }, { TX_16X16, TX_16X16 } },
      { { TX_32X32, TX_32X32 }, { TX_32X32, TX_32X32 } },
      { { TX_4X8, TX_4X8 }, { TX_4X8, TX_4X8 } },
      { { TX_8X4, TX_8X4 }, { TX_8X4, TX_8X4 } },
      { { TX_8X16, TX_8X16 }, { TX_8X16, TX_8X16 } },
      { { TX_16X8, TX_16X8 }, { TX_16X8, TX_16X8 } },
      { { TX_16X32, TX_16X32 }, { TX_16X32, TX_16X32 } },
      { { TX_32X16, TX_32X16 }, { TX_32X16, TX_32X16 } },
  },
  {
// BLOCK_128X128
#if CONFIG_CB4X4
      { { TX_4X4, TX_4X4 }, { TX_4X4, TX_4X4 } },
#endif
      { { TX_4X4, TX_4X4 }, { TX_4X4, TX_4X4 } },
      { { TX_8X8, TX_8X8 }, { TX_8X8, TX_8X8 } },
      { { TX_16X16, TX_16X16 }, { TX_16X16, TX_16X16 } },
      { { TX_32X32, TX_32X32 }, { TX_32X32, TX_32X32 } },
      { { TX_4X8, TX_4X8 }, { TX_4X8, TX_4X8 } },
      { { TX_8X4, TX_8X4 }, { TX_8X4, TX_8X4 } },
      { { TX_8X16, TX_8X16 }, { TX_8X16, TX_8X16 } },
      { { TX_16X8, TX_16X8 }, { TX_16X8, TX_16X8 } },
      { { TX_16X32, TX_16X32 }, { TX_16X32, TX_16X32 } },
      { { TX_32X16, TX_32X16 }, { TX_32X16, TX_32X16 } },
#endif  // CONFIG_EXT_PARTITION
  },
};

// Generates 4 bit field in which each bit set to 1 represents
// a blocksize partition  1111 means we split 64x64, 32x32, 16x16
// and 8x8.  1000 means we just split the 64x64 to 32x32
static const struct {
  PARTITION_CONTEXT above;
  PARTITION_CONTEXT left;
} partition_context_lookup[BLOCK_SIZES] = {
#if CONFIG_EXT_PARTITION
  { 31, 31 },  // 4X4   - {0b11111, 0b11111}
  { 31, 30 },  // 4X8   - {0b11111, 0b11110}
  { 30, 31 },  // 8X4   - {0b11110, 0b11111}
  { 30, 30 },  // 8X8   - {0b11110, 0b11110}
  { 30, 28 },  // 8X16  - {0b11110, 0b11100}
  { 28, 30 },  // 16X8  - {0b11100, 0b11110}
  { 28, 28 },  // 16X16 - {0b11100, 0b11100}
  { 28, 24 },  // 16X32 - {0b11100, 0b11000}
  { 24, 28 },  // 32X16 - {0b11000, 0b11100}
  { 24, 24 },  // 32X32 - {0b11000, 0b11000}
  { 24, 16 },  // 32X64 - {0b11000, 0b10000}
  { 16, 24 },  // 64X32 - {0b10000, 0b11000}
  { 16, 16 },  // 64X64 - {0b10000, 0b10000}
  { 16, 0 },   // 64X128- {0b10000, 0b00000}
  { 0, 16 },   // 128X64- {0b00000, 0b10000}
  { 0, 0 },    // 128X128-{0b00000, 0b00000}
#else
  { 15, 15 },  // 4X4   - {0b1111, 0b1111}
  { 15, 14 },  // 4X8   - {0b1111, 0b1110}
  { 14, 15 },  // 8X4   - {0b1110, 0b1111}
  { 14, 14 },  // 8X8   - {0b1110, 0b1110}
  { 14, 12 },  // 8X16  - {0b1110, 0b1100}
  { 12, 14 },  // 16X8  - {0b1100, 0b1110}
  { 12, 12 },  // 16X16 - {0b1100, 0b1100}
  { 12, 8 },   // 16X32 - {0b1100, 0b1000}
  { 8, 12 },   // 32X16 - {0b1000, 0b1100}
  { 8, 8 },    // 32X32 - {0b1000, 0b1000}
  { 8, 0 },    // 32X64 - {0b1000, 0b0000}
  { 0, 8 },    // 64X32 - {0b0000, 0b1000}
  { 0, 0 },    // 64X64 - {0b0000, 0b0000}
#endif  // CONFIG_EXT_PARTITION
};

#if CONFIG_SUPERTX
static const TX_SIZE uvsupertx_size_lookup[TX_SIZES][2][2] = {
//  ss_x == 0 ss_x == 0   ss_x == 1 ss_x == 1
//  ss_y == 0 ss_y == 1   ss_y == 0 ss_y == 1
#if CONFIG_CB4X4
  { { TX_4X4, TX_4X4 }, { TX_4X4, TX_4X4 } },
#endif
  { { TX_4X4, TX_4X4 }, { TX_4X4, TX_4X4 } },
  { { TX_8X8, TX_4X4 }, { TX_4X4, TX_4X4 } },
  { { TX_16X16, TX_8X8 }, { TX_8X8, TX_8X8 } },
  { { TX_32X32, TX_16X16 }, { TX_16X16, TX_16X16 } },
};

#if CONFIG_EXT_PARTITION_TYPES
static const int partition_supertx_context_lookup[EXT_PARTITION_TYPES] = {
  -1, 0, 0, 1, 0, 0, 0, 0
};

#else
static const int partition_supertx_context_lookup[PARTITION_TYPES] = { -1, 0, 0,
                                                                       1 };
#endif  // CONFIG_EXT_PARTITION_TYPES
#endif  // CONFIG_SUPERTX

#ifdef __cplusplus
}  // extern "C"
#endif

#endif  // AV1_COMMON_COMMON_DATA_H_
