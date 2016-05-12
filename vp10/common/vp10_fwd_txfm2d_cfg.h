/*
 *  Copyright (c) 2015 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#ifndef VP10_FWD_TXFM2D_CFG_H_
#define VP10_FWD_TXFM2D_CFG_H_
#include "vp10/common/enums.h"
#include "vp10/common/vp10_fwd_txfm1d.h"
//  ---------------- config fwd_dct_dct_4 ----------------
static const int8_t fwd_shift_dct_dct_4[3] = {2, 0, 0};
static const int8_t fwd_stage_range_col_dct_dct_4[4] = {15, 16, 17, 17};
static const int8_t fwd_stage_range_row_dct_dct_4[4] = {17, 18, 18, 18};
static const int8_t fwd_cos_bit_col_dct_dct_4[4] = {13, 13, 13, 13};
static const int8_t fwd_cos_bit_row_dct_dct_4[4] = {13, 13, 13, 13};

static const TXFM_2D_CFG fwd_txfm_2d_cfg_dct_dct_4 = {
    4,  // .txfm_size
    4,  // .stage_num_col
    4,  // .stage_num_row
    // 0,  // .log_scale
    fwd_shift_dct_dct_4,            // .shift
    fwd_stage_range_col_dct_dct_4,  // .stage_range_col
    fwd_stage_range_row_dct_dct_4,  // .stage_range_row
    fwd_cos_bit_col_dct_dct_4,      // .cos_bit_col
    fwd_cos_bit_row_dct_dct_4,      // .cos_bit_row
    TXFM_TYPE_DCT4,                 // .txfm_type_col
    TXFM_TYPE_DCT4};                // .txfm_type_row

//  ---------------- config fwd_dct_dct_8 ----------------
static const int8_t fwd_shift_dct_dct_8[3] = {2, -1, 0};
static const int8_t fwd_stage_range_col_dct_dct_8[6] = {15, 16, 17, 18, 18, 18};
static const int8_t fwd_stage_range_row_dct_dct_8[6] = {17, 18, 19, 19, 19, 19};
static const int8_t fwd_cos_bit_col_dct_dct_8[6] = {13, 13, 13, 13, 13, 13};
static const int8_t fwd_cos_bit_row_dct_dct_8[6] = {13, 13, 13, 13, 13, 13};

static const TXFM_2D_CFG fwd_txfm_2d_cfg_dct_dct_8 = {
    8,  // .txfm_size
    6,  // .stage_num_col
    6,  // .stage_num_row
    // 0,  // .log_scale
    fwd_shift_dct_dct_8,            // .shift
    fwd_stage_range_col_dct_dct_8,  // .stage_range_col
    fwd_stage_range_row_dct_dct_8,  // .stage_range_row
    fwd_cos_bit_col_dct_dct_8,      // .cos_bit_col
    fwd_cos_bit_row_dct_dct_8,      // .cos_bit_row
    TXFM_TYPE_DCT8,                 // .txfm_type_col
    TXFM_TYPE_DCT8};                // .txfm_type_row

//  ---------------- config fwd_dct_dct_16 ----------------
static const int8_t fwd_shift_dct_dct_16[3] = {2, -2, 0};
static const int8_t fwd_stage_range_col_dct_dct_16[8] = {15, 16, 17, 18,
                                                         19, 19, 19, 19};
static const int8_t fwd_stage_range_row_dct_dct_16[8] = {17, 18, 19, 20,
                                                         20, 20, 20, 20};
static const int8_t fwd_cos_bit_col_dct_dct_16[8] = {13, 13, 13, 13,
                                                     13, 13, 13, 13};
static const int8_t fwd_cos_bit_row_dct_dct_16[8] = {12, 12, 12, 12,
                                                     12, 12, 12, 12};

static const TXFM_2D_CFG fwd_txfm_2d_cfg_dct_dct_16 = {
    16,  // .txfm_size
    8,   // .stage_num_col
    8,   // .stage_num_row
    // 0,  // .log_scale
    fwd_shift_dct_dct_16,            // .shift
    fwd_stage_range_col_dct_dct_16,  // .stage_range_col
    fwd_stage_range_row_dct_dct_16,  // .stage_range_row
    fwd_cos_bit_col_dct_dct_16,      // .cos_bit_col
    fwd_cos_bit_row_dct_dct_16,      // .cos_bit_row
    TXFM_TYPE_DCT16,                 // .txfm_type_col
    TXFM_TYPE_DCT16};                // .txfm_type_row

//  ---------------- config fwd_dct_dct_32 ----------------
static const int8_t fwd_shift_dct_dct_32[3] = {2, -4, 0};
static const int8_t fwd_stage_range_col_dct_dct_32[10] = {15, 16, 17, 18, 19,
                                                          20, 20, 20, 20, 20};
static const int8_t fwd_stage_range_row_dct_dct_32[10] = {16, 17, 18, 19, 20,
                                                          20, 20, 20, 20, 20};
static const int8_t fwd_cos_bit_col_dct_dct_32[10] = {12, 12, 12, 12, 12,
                                                      12, 12, 12, 12, 12};
static const int8_t fwd_cos_bit_row_dct_dct_32[10] = {12, 12, 12, 12, 12,
                                                      12, 12, 12, 12, 12};

static const TXFM_2D_CFG fwd_txfm_2d_cfg_dct_dct_32 = {
    32,  // .txfm_size
    10,  // .stage_num_col
    10,  // .stage_num_row
    // 1,  // .log_scale
    fwd_shift_dct_dct_32,            // .shift
    fwd_stage_range_col_dct_dct_32,  // .stage_range_col
    fwd_stage_range_row_dct_dct_32,  // .stage_range_row
    fwd_cos_bit_col_dct_dct_32,      // .cos_bit_col
    fwd_cos_bit_row_dct_dct_32,      // .cos_bit_row
    TXFM_TYPE_DCT32,                 // .txfm_type_col
    TXFM_TYPE_DCT32};                // .txfm_type_row

//  ---------------- config fwd_dct_dct_64 ----------------
static const int8_t fwd_shift_dct_dct_64[3] = {2, -2, -2};
static const int8_t fwd_stage_range_col_dct_dct_64[12] = {
    13, 14, 15, 16, 17, 18, 19, 19, 19, 19, 19, 19};
static const int8_t fwd_stage_range_row_dct_dct_64[12] = {
    17, 18, 19, 20, 21, 22, 22, 22, 22, 22, 22, 22};
static const int8_t fwd_cos_bit_col_dct_dct_64[12] = {15, 15, 15, 15, 15, 14,
                                                      13, 13, 13, 13, 13, 13};
static const int8_t fwd_cos_bit_row_dct_dct_64[12] = {15, 14, 13, 12, 11, 10,
                                                      10, 10, 10, 10, 10, 10};

static const TXFM_2D_CFG fwd_txfm_2d_cfg_dct_dct_64 = {
    64,                              // .txfm_size
    12,                              // .stage_num_col
    12,                              // .stage_num_row
    fwd_shift_dct_dct_64,            // .shift
    fwd_stage_range_col_dct_dct_64,  // .stage_range_col
    fwd_stage_range_row_dct_dct_64,  // .stage_range_row
    fwd_cos_bit_col_dct_dct_64,      // .cos_bit_col
    fwd_cos_bit_row_dct_dct_64,      // .cos_bit_row
    TXFM_TYPE_DCT64,                 // .txfm_type_col
    TXFM_TYPE_DCT64};                // .txfm_type_row

//  ---------------- config fwd_dct_adst_4 ----------------
static const int8_t fwd_shift_dct_adst_4[3] = {2, 0, 0};
static const int8_t fwd_stage_range_col_dct_adst_4[4] = {15, 16, 17, 17};
static const int8_t fwd_stage_range_row_dct_adst_4[6] = {17, 17, 17,
                                                         18, 18, 18};
static const int8_t fwd_cos_bit_col_dct_adst_4[4] = {13, 13, 13, 13};
static const int8_t fwd_cos_bit_row_dct_adst_4[6] = {13, 13, 13, 13, 13, 13};

static const TXFM_2D_CFG fwd_txfm_2d_cfg_dct_adst_4 = {
    4,  // .txfm_size
    4,  // .stage_num_col
    6,  // .stage_num_row
    // 0,  // .log_scale
    fwd_shift_dct_adst_4,            // .shift
    fwd_stage_range_col_dct_adst_4,  // .stage_range_col
    fwd_stage_range_row_dct_adst_4,  // .stage_range_row
    fwd_cos_bit_col_dct_adst_4,      // .cos_bit_col
    fwd_cos_bit_row_dct_adst_4,      // .cos_bit_row
    TXFM_TYPE_DCT4,                  // .txfm_type_col
    TXFM_TYPE_ADST4};                // .txfm_type_row

//  ---------------- config fwd_dct_adst_8 ----------------
static const int8_t fwd_shift_dct_adst_8[3] = {2, -1, 0};
static const int8_t fwd_stage_range_col_dct_adst_8[6] = {15, 16, 17,
                                                         18, 18, 18};
static const int8_t fwd_stage_range_row_dct_adst_8[8] = {17, 17, 17, 18,
                                                         18, 19, 19, 19};
static const int8_t fwd_cos_bit_col_dct_adst_8[6] = {13, 13, 13, 13, 13, 13};
static const int8_t fwd_cos_bit_row_dct_adst_8[8] = {13, 13, 13, 13,
                                                     13, 13, 13, 13};

static const TXFM_2D_CFG fwd_txfm_2d_cfg_dct_adst_8 = {
    8,  // .txfm_size
    6,  // .stage_num_col
    8,  // .stage_num_row
    // 0,  // .log_scale
    fwd_shift_dct_adst_8,            // .shift
    fwd_stage_range_col_dct_adst_8,  // .stage_range_col
    fwd_stage_range_row_dct_adst_8,  // .stage_range_row
    fwd_cos_bit_col_dct_adst_8,      // .cos_bit_col
    fwd_cos_bit_row_dct_adst_8,      // .cos_bit_row
    TXFM_TYPE_DCT8,                  // .txfm_type_col
    TXFM_TYPE_ADST8};                // .txfm_type_row

//  ---------------- config fwd_dct_adst_16 ----------------
static const int8_t fwd_shift_dct_adst_16[3] = {2, -2, 0};
static const int8_t fwd_stage_range_col_dct_adst_16[8] = {15, 16, 17, 18,
                                                          19, 19, 19, 19};
static const int8_t fwd_stage_range_row_dct_adst_16[10] = {17, 17, 17, 18, 18,
                                                           19, 19, 20, 20, 20};
static const int8_t fwd_cos_bit_col_dct_adst_16[8] = {13, 13, 13, 13,
                                                      13, 13, 13, 13};
static const int8_t fwd_cos_bit_row_dct_adst_16[10] = {12, 12, 12, 12, 12,
                                                       12, 12, 12, 12, 12};

static const TXFM_2D_CFG fwd_txfm_2d_cfg_dct_adst_16 = {
    16,  // .txfm_size
    8,   // .stage_num_col
    10,  // .stage_num_row
    // 0,  // .log_scale
    fwd_shift_dct_adst_16,            // .shift
    fwd_stage_range_col_dct_adst_16,  // .stage_range_col
    fwd_stage_range_row_dct_adst_16,  // .stage_range_row
    fwd_cos_bit_col_dct_adst_16,      // .cos_bit_col
    fwd_cos_bit_row_dct_adst_16,      // .cos_bit_row
    TXFM_TYPE_DCT16,                  // .txfm_type_col
    TXFM_TYPE_ADST16};                // .txfm_type_row

//  ---------------- config fwd_dct_adst_32 ----------------
static const int8_t fwd_shift_dct_adst_32[3] = {2, -4, 0};
static const int8_t fwd_stage_range_col_dct_adst_32[10] = {15, 16, 17, 18, 19,
                                                           20, 20, 20, 20, 20};
static const int8_t fwd_stage_range_row_dct_adst_32[12] = {
    16, 16, 16, 17, 17, 18, 18, 19, 19, 20, 20, 20};
static const int8_t fwd_cos_bit_col_dct_adst_32[10] = {12, 12, 12, 12, 12,
                                                       12, 12, 12, 12, 12};
static const int8_t fwd_cos_bit_row_dct_adst_32[12] = {12, 12, 12, 12, 12, 12,
                                                       12, 12, 12, 12, 12, 12};

static const TXFM_2D_CFG fwd_txfm_2d_cfg_dct_adst_32 = {
    32,  // .txfm_size
    10,  // .stage_num_col
    12,  // .stage_num_row
    // 1,  // .log_scale
    fwd_shift_dct_adst_32,            // .shift
    fwd_stage_range_col_dct_adst_32,  // .stage_range_col
    fwd_stage_range_row_dct_adst_32,  // .stage_range_row
    fwd_cos_bit_col_dct_adst_32,      // .cos_bit_col
    fwd_cos_bit_row_dct_adst_32,      // .cos_bit_row
    TXFM_TYPE_DCT32,                  // .txfm_type_col
    TXFM_TYPE_ADST32};                // .txfm_type_row
//  ---------------- config fwd_adst_adst_4 ----------------
static const int8_t fwd_shift_adst_adst_4[3] = {2, 0, 0};
static const int8_t fwd_stage_range_col_adst_adst_4[6] = {15, 15, 16,
                                                          17, 17, 17};
static const int8_t fwd_stage_range_row_adst_adst_4[6] = {17, 17, 17,
                                                          18, 18, 18};
static const int8_t fwd_cos_bit_col_adst_adst_4[6] = {13, 13, 13, 13, 13, 13};
static const int8_t fwd_cos_bit_row_adst_adst_4[6] = {13, 13, 13, 13, 13, 13};

static const TXFM_2D_CFG fwd_txfm_2d_cfg_adst_adst_4 = {
    4,  // .txfm_size
    6,  // .stage_num_col
    6,  // .stage_num_row
    // 0,  // .log_scale
    fwd_shift_adst_adst_4,            // .shift
    fwd_stage_range_col_adst_adst_4,  // .stage_range_col
    fwd_stage_range_row_adst_adst_4,  // .stage_range_row
    fwd_cos_bit_col_adst_adst_4,      // .cos_bit_col
    fwd_cos_bit_row_adst_adst_4,      // .cos_bit_row
    TXFM_TYPE_ADST4,                  // .txfm_type_col
    TXFM_TYPE_ADST4};                 // .txfm_type_row

//  ---------------- config fwd_adst_adst_8 ----------------
static const int8_t fwd_shift_adst_adst_8[3] = {2, -1, 0};
static const int8_t fwd_stage_range_col_adst_adst_8[8] = {15, 15, 16, 17,
                                                          17, 18, 18, 18};
static const int8_t fwd_stage_range_row_adst_adst_8[8] = {17, 17, 17, 18,
                                                          18, 19, 19, 19};
static const int8_t fwd_cos_bit_col_adst_adst_8[8] = {13, 13, 13, 13,
                                                      13, 13, 13, 13};
static const int8_t fwd_cos_bit_row_adst_adst_8[8] = {13, 13, 13, 13,
                                                      13, 13, 13, 13};

static const TXFM_2D_CFG fwd_txfm_2d_cfg_adst_adst_8 = {
    8,  // .txfm_size
    8,  // .stage_num_col
    8,  // .stage_num_row
    // 0,  // .log_scale
    fwd_shift_adst_adst_8,            // .shift
    fwd_stage_range_col_adst_adst_8,  // .stage_range_col
    fwd_stage_range_row_adst_adst_8,  // .stage_range_row
    fwd_cos_bit_col_adst_adst_8,      // .cos_bit_col
    fwd_cos_bit_row_adst_adst_8,      // .cos_bit_row
    TXFM_TYPE_ADST8,                  // .txfm_type_col
    TXFM_TYPE_ADST8};                 // .txfm_type_row

//  ---------------- config fwd_adst_adst_16 ----------------
static const int8_t fwd_shift_adst_adst_16[3] = {2, -2, 0};
static const int8_t fwd_stage_range_col_adst_adst_16[10] = {15, 15, 16, 17, 17,
                                                            18, 18, 19, 19, 19};
static const int8_t fwd_stage_range_row_adst_adst_16[10] = {17, 17, 17, 18, 18,
                                                            19, 19, 20, 20, 20};
static const int8_t fwd_cos_bit_col_adst_adst_16[10] = {13, 13, 13, 13, 13,
                                                        13, 13, 13, 13, 13};
static const int8_t fwd_cos_bit_row_adst_adst_16[10] = {12, 12, 12, 12, 12,
                                                        12, 12, 12, 12, 12};

static const TXFM_2D_CFG fwd_txfm_2d_cfg_adst_adst_16 = {
    16,  // .txfm_size
    10,  // .stage_num_col
    10,  // .stage_num_row
    // 0,  // .log_scale
    fwd_shift_adst_adst_16,            // .shift
    fwd_stage_range_col_adst_adst_16,  // .stage_range_col
    fwd_stage_range_row_adst_adst_16,  // .stage_range_row
    fwd_cos_bit_col_adst_adst_16,      // .cos_bit_col
    fwd_cos_bit_row_adst_adst_16,      // .cos_bit_row
    TXFM_TYPE_ADST16,                  // .txfm_type_col
    TXFM_TYPE_ADST16};                 // .txfm_type_row

//  ---------------- config fwd_adst_adst_32 ----------------
static const int8_t fwd_shift_adst_adst_32[3] = {2, -4, 0};
static const int8_t fwd_stage_range_col_adst_adst_32[12] = {
    15, 15, 16, 17, 17, 18, 18, 19, 19, 20, 20, 20};
static const int8_t fwd_stage_range_row_adst_adst_32[12] = {
    16, 16, 16, 17, 17, 18, 18, 19, 19, 20, 20, 20};
static const int8_t fwd_cos_bit_col_adst_adst_32[12] = {12, 12, 12, 12, 12, 12,
                                                        12, 12, 12, 12, 12, 12};
static const int8_t fwd_cos_bit_row_adst_adst_32[12] = {12, 12, 12, 12, 12, 12,
                                                        12, 12, 12, 12, 12, 12};

static const TXFM_2D_CFG fwd_txfm_2d_cfg_adst_adst_32 = {
    32,  // .txfm_size
    12,  // .stage_num_col
    12,  // .stage_num_row
    // 1,  // .log_scale
    fwd_shift_adst_adst_32,            // .shift
    fwd_stage_range_col_adst_adst_32,  // .stage_range_col
    fwd_stage_range_row_adst_adst_32,  // .stage_range_row
    fwd_cos_bit_col_adst_adst_32,      // .cos_bit_col
    fwd_cos_bit_row_adst_adst_32,      // .cos_bit_row
    TXFM_TYPE_ADST32,                  // .txfm_type_col
    TXFM_TYPE_ADST32};                 // .txfm_type_row

//  ---------------- config fwd_adst_dct_4 ----------------
static const int8_t fwd_shift_adst_dct_4[3] = {2, 0, 0};
static const int8_t fwd_stage_range_col_adst_dct_4[6] = {15, 15, 16,
                                                         17, 17, 17};
static const int8_t fwd_stage_range_row_adst_dct_4[4] = {17, 18, 18, 18};
static const int8_t fwd_cos_bit_col_adst_dct_4[6] = {13, 13, 13, 13, 13, 13};
static const int8_t fwd_cos_bit_row_adst_dct_4[4] = {13, 13, 13, 13};

static const TXFM_2D_CFG fwd_txfm_2d_cfg_adst_dct_4 = {
    4,  // .txfm_size
    6,  // .stage_num_col
    4,  // .stage_num_row
    // 0,  // .log_scale
    fwd_shift_adst_dct_4,            // .shift
    fwd_stage_range_col_adst_dct_4,  // .stage_range_col
    fwd_stage_range_row_adst_dct_4,  // .stage_range_row
    fwd_cos_bit_col_adst_dct_4,      // .cos_bit_col
    fwd_cos_bit_row_adst_dct_4,      // .cos_bit_row
    TXFM_TYPE_ADST4,                 // .txfm_type_col
    TXFM_TYPE_DCT4};                 // .txfm_type_row

//  ---------------- config fwd_adst_dct_8 ----------------
static const int8_t fwd_shift_adst_dct_8[3] = {2, -1, 0};
static const int8_t fwd_stage_range_col_adst_dct_8[8] = {15, 15, 16, 17,
                                                         17, 18, 18, 18};
static const int8_t fwd_stage_range_row_adst_dct_8[6] = {17, 18, 19,
                                                         19, 19, 19};
static const int8_t fwd_cos_bit_col_adst_dct_8[8] = {13, 13, 13, 13,
                                                     13, 13, 13, 13};
static const int8_t fwd_cos_bit_row_adst_dct_8[6] = {13, 13, 13, 13, 13, 13};

static const TXFM_2D_CFG fwd_txfm_2d_cfg_adst_dct_8 = {
    8,  // .txfm_size
    8,  // .stage_num_col
    6,  // .stage_num_row
    // 0,  // .log_scale
    fwd_shift_adst_dct_8,            // .shift
    fwd_stage_range_col_adst_dct_8,  // .stage_range_col
    fwd_stage_range_row_adst_dct_8,  // .stage_range_row
    fwd_cos_bit_col_adst_dct_8,      // .cos_bit_col
    fwd_cos_bit_row_adst_dct_8,      // .cos_bit_row
    TXFM_TYPE_ADST8,                 // .txfm_type_col
    TXFM_TYPE_DCT8};                 // .txfm_type_row

//  ---------------- config fwd_adst_dct_16 ----------------
static const int8_t fwd_shift_adst_dct_16[3] = {2, -2, 0};
static const int8_t fwd_stage_range_col_adst_dct_16[10] = {15, 15, 16, 17, 17,
                                                           18, 18, 19, 19, 19};
static const int8_t fwd_stage_range_row_adst_dct_16[8] = {17, 18, 19, 20,
                                                          20, 20, 20, 20};
static const int8_t fwd_cos_bit_col_adst_dct_16[10] = {13, 13, 13, 13, 13,
                                                       13, 13, 13, 13, 13};
static const int8_t fwd_cos_bit_row_adst_dct_16[8] = {12, 12, 12, 12,
                                                      12, 12, 12, 12};

static const TXFM_2D_CFG fwd_txfm_2d_cfg_adst_dct_16 = {
    16,  // .txfm_size
    10,  // .stage_num_col
    8,   // .stage_num_row
    // 0,  // .log_scale
    fwd_shift_adst_dct_16,            // .shift
    fwd_stage_range_col_adst_dct_16,  // .stage_range_col
    fwd_stage_range_row_adst_dct_16,  // .stage_range_row
    fwd_cos_bit_col_adst_dct_16,      // .cos_bit_col
    fwd_cos_bit_row_adst_dct_16,      // .cos_bit_row
    TXFM_TYPE_ADST16,                 // .txfm_type_col
    TXFM_TYPE_DCT16};                 // .txfm_type_row

//  ---------------- config fwd_adst_dct_32 ----------------
static const int8_t fwd_shift_adst_dct_32[3] = {2, -4, 0};
static const int8_t fwd_stage_range_col_adst_dct_32[12] = {
    15, 15, 16, 17, 17, 18, 18, 19, 19, 20, 20, 20};
static const int8_t fwd_stage_range_row_adst_dct_32[10] = {16, 17, 18, 19, 20,
                                                           20, 20, 20, 20, 20};
static const int8_t fwd_cos_bit_col_adst_dct_32[12] = {12, 12, 12, 12, 12, 12,
                                                       12, 12, 12, 12, 12, 12};
static const int8_t fwd_cos_bit_row_adst_dct_32[10] = {12, 12, 12, 12, 12,
                                                       12, 12, 12, 12, 12};

static const TXFM_2D_CFG fwd_txfm_2d_cfg_adst_dct_32 = {
    32,  // .txfm_size
    12,  // .stage_num_col
    10,  // .stage_num_row
    // 1,  // .log_scale
    fwd_shift_adst_dct_32,            // .shift
    fwd_stage_range_col_adst_dct_32,  // .stage_range_col
    fwd_stage_range_row_adst_dct_32,  // .stage_range_row
    fwd_cos_bit_col_adst_dct_32,      // .cos_bit_col
    fwd_cos_bit_row_adst_dct_32,      // .cos_bit_row
    TXFM_TYPE_ADST32,                 // .txfm_type_col
    TXFM_TYPE_DCT32};                 // .txfm_type_row
#endif  // VP10_FWD_TXFM2D_CFG_H_
