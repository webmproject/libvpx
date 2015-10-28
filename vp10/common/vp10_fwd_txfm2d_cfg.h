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
#include "vp10/common/vp10_fwd_txfm1d.h"

//  ---------------- config fwd_dct_dct_4 ----------------
static int8_t fwd_shift_dct_dct_4[3] = {4, 0, -2};
static int8_t fwd_stage_range_col_dct_dct_4[4] = {15, 16, 17, 17};
static int8_t fwd_stage_range_row_dct_dct_4[4] = {17, 18, 18, 18};
static int8_t fwd_cos_bit_col_dct_dct_4[4] = {15, 15, 15, 15};
static int8_t fwd_cos_bit_row_dct_dct_4[4] = {15, 14, 14, 14};

static const TXFM_2D_CFG fwd_txfm_2d_cfg_dct_dct_4 = {
    .txfm_size = 4,
    .stage_num_col = 4,
    .stage_num_row = 4,

    .shift = fwd_shift_dct_dct_4,
    .stage_range_col = fwd_stage_range_col_dct_dct_4,
    .stage_range_row = fwd_stage_range_row_dct_dct_4,
    .cos_bit_col = fwd_cos_bit_col_dct_dct_4,
    .cos_bit_row = fwd_cos_bit_row_dct_dct_4,
    .txfm_func_col = vp10_fdct4_new,
    .txfm_func_row = vp10_fdct4_new};

//  ---------------- config fwd_dct_dct_8 ----------------
static int8_t fwd_shift_dct_dct_8[3] = {5, -3, -1};
static int8_t fwd_stage_range_col_dct_dct_8[6] = {16, 17, 18, 19, 19, 19};
static int8_t fwd_stage_range_row_dct_dct_8[6] = {16, 17, 18, 18, 18, 18};
static int8_t fwd_cos_bit_col_dct_dct_8[6] = {15, 15, 14, 13, 13, 13};
static int8_t fwd_cos_bit_row_dct_dct_8[6] = {15, 15, 14, 14, 14, 14};

static const TXFM_2D_CFG fwd_txfm_2d_cfg_dct_dct_8 = {
    .txfm_size = 8,
    .stage_num_col = 6,
    .stage_num_row = 6,

    .shift = fwd_shift_dct_dct_8,
    .stage_range_col = fwd_stage_range_col_dct_dct_8,
    .stage_range_row = fwd_stage_range_row_dct_dct_8,
    .cos_bit_col = fwd_cos_bit_col_dct_dct_8,
    .cos_bit_row = fwd_cos_bit_row_dct_dct_8,
    .txfm_func_col = vp10_fdct8_new,
    .txfm_func_row = vp10_fdct8_new};

//  ---------------- config fwd_dct_dct_16 ----------------
static int8_t fwd_shift_dct_dct_16[3] = {4, -3, -1};
static int8_t fwd_stage_range_col_dct_dct_16[8] = {15, 16, 17, 18,
                                                   19, 19, 19, 19};
static int8_t fwd_stage_range_row_dct_dct_16[8] = {16, 17, 18, 19,
                                                   19, 19, 19, 19};
static int8_t fwd_cos_bit_col_dct_dct_16[8] = {15, 15, 15, 14, 13, 13, 13, 13};
static int8_t fwd_cos_bit_row_dct_dct_16[8] = {15, 15, 14, 13, 13, 13, 13, 13};

static const TXFM_2D_CFG fwd_txfm_2d_cfg_dct_dct_16 = {
    .txfm_size = 16,
    .stage_num_col = 8,
    .stage_num_row = 8,

    .shift = fwd_shift_dct_dct_16,
    .stage_range_col = fwd_stage_range_col_dct_dct_16,
    .stage_range_row = fwd_stage_range_row_dct_dct_16,
    .cos_bit_col = fwd_cos_bit_col_dct_dct_16,
    .cos_bit_row = fwd_cos_bit_row_dct_dct_16,
    .txfm_func_col = vp10_fdct16_new,
    .txfm_func_row = vp10_fdct16_new};

//  ---------------- config fwd_dct_dct_32 ----------------
static int8_t fwd_shift_dct_dct_32[3] = {3, -3, -1};
static int8_t fwd_stage_range_col_dct_dct_32[10] = {14, 15, 16, 17, 18,
                                                    19, 19, 19, 19, 19};
static int8_t fwd_stage_range_row_dct_dct_32[10] = {16, 17, 18, 19, 20,
                                                    20, 20, 20, 20, 20};
static int8_t fwd_cos_bit_col_dct_dct_32[10] = {15, 15, 15, 15, 14,
                                                13, 13, 13, 13, 13};
static int8_t fwd_cos_bit_row_dct_dct_32[10] = {15, 15, 14, 13, 12,
                                                12, 12, 12, 12, 12};

static const TXFM_2D_CFG fwd_txfm_2d_cfg_dct_dct_32 = {
    .txfm_size = 32,
    .stage_num_col = 10,
    .stage_num_row = 10,

    .shift = fwd_shift_dct_dct_32,
    .stage_range_col = fwd_stage_range_col_dct_dct_32,
    .stage_range_row = fwd_stage_range_row_dct_dct_32,
    .cos_bit_col = fwd_cos_bit_col_dct_dct_32,
    .cos_bit_row = fwd_cos_bit_row_dct_dct_32,
    .txfm_func_col = vp10_fdct32_new,
    .txfm_func_row = vp10_fdct32_new};

#endif  // VP10_FWD_TXFM2D_CFG_H_
