/*
 *  Copyright (c) 2016 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include <assert.h>

#include "./vp10_rtcd.h"
#include "vp10/common/enums.h"
#include "vp10/common/vp10_fwd_txfm1d.h"
#include "vp10/common/vp10_fwd_txfm2d_cfg.h"
#include "vp10/common/vp10_txfm.h"

static INLINE TxfmFunc fwd_txfm_type_to_func(TXFM_TYPE txfm_type) {
  switch (txfm_type) {
    case TXFM_TYPE_DCT4:
      return vp10_fdct4_new;
    case TXFM_TYPE_DCT8:
      return vp10_fdct8_new;
    case TXFM_TYPE_DCT16:
      return vp10_fdct16_new;
    case TXFM_TYPE_DCT32:
      return vp10_fdct32_new;
    case TXFM_TYPE_DCT64:
      return vp10_fdct64_new;
    case TXFM_TYPE_ADST4:
      return vp10_fadst4_new;
    case TXFM_TYPE_ADST8:
      return vp10_fadst8_new;
    case TXFM_TYPE_ADST16:
      return vp10_fadst16_new;
    case TXFM_TYPE_ADST32:
      return vp10_fadst32_new;
    default:
      assert(0);
      return NULL;
  }
}

static INLINE void fwd_txfm2d_c(const int16_t *input, int32_t *output,
                                const int stride, const TXFM_2D_FLIP_CFG *cfg,
                                int32_t *buf) {
  int c, r;
  const int txfm_size = cfg->cfg->txfm_size;
  const int8_t *shift = cfg->cfg->shift;
  const int8_t *stage_range_col = cfg->cfg->stage_range_col;
  const int8_t *stage_range_row = cfg->cfg->stage_range_row;
  const int8_t *cos_bit_col = cfg->cfg->cos_bit_col;
  const int8_t *cos_bit_row = cfg->cfg->cos_bit_row;
  const TxfmFunc txfm_func_col = fwd_txfm_type_to_func(cfg->cfg->txfm_type_col);
  const TxfmFunc txfm_func_row = fwd_txfm_type_to_func(cfg->cfg->txfm_type_row);

  // use output buffer as temp buffer
  int32_t* temp_in = output;
  int32_t* temp_out = output + txfm_size;

  // Columns
  for (c = 0; c < txfm_size; ++c) {
    if (cfg->ud_flip == 0) {
      for (r = 0; r < txfm_size; ++r)
        temp_in[r] = input[r * stride + c];
    } else {
      for (r = 0; r < txfm_size; ++r)
        // flip upside down
        temp_in[r] = input[(txfm_size - r - 1) * stride + c];
    }
    round_shift_array(temp_in, txfm_size, -shift[0]);
    txfm_func_col(temp_in, temp_out, cos_bit_col, stage_range_col);
    round_shift_array(temp_out, txfm_size, -shift[1]);
    if (cfg->lr_flip == 0) {
      for (r = 0; r < txfm_size; ++r)
        buf[r * txfm_size + c] = temp_out[r];
    } else {
      for (r = 0; r < txfm_size; ++r)
        // flip from left to right
        buf[r * txfm_size + (txfm_size - c - 1)] = temp_out[r];
    }
  }

  // Rows
  for (r = 0; r < txfm_size; ++r) {
    txfm_func_row(buf + r * txfm_size, output + r * txfm_size, cos_bit_row,
                  stage_range_row);
    round_shift_array(output + r * txfm_size, txfm_size, -shift[2]);
  }
}

void vp10_fwd_txfm2d_4x4_c(const int16_t *input, int32_t *output,
                           const int stride, int tx_type,
                           const int bd) {
  int32_t txfm_buf[4 * 4];
  TXFM_2D_FLIP_CFG cfg = vp10_get_fwd_txfm_cfg(tx_type, TX_4X4);
  (void)bd;
  fwd_txfm2d_c(input, output, stride, &cfg, txfm_buf);
}

void vp10_fwd_txfm2d_8x8_c(const int16_t *input, int32_t *output,
                           const int stride, int tx_type,
                           const int bd) {
  int32_t txfm_buf[8 * 8];
  TXFM_2D_FLIP_CFG cfg = vp10_get_fwd_txfm_cfg(tx_type, TX_8X8);
  (void)bd;
  fwd_txfm2d_c(input, output, stride, &cfg, txfm_buf);
}

void vp10_fwd_txfm2d_16x16_c(const int16_t *input, int32_t *output,
                             const int stride, int tx_type,
                             const int bd) {
  int32_t txfm_buf[16 * 16];
  TXFM_2D_FLIP_CFG cfg = vp10_get_fwd_txfm_cfg(tx_type, TX_16X16);
  (void)bd;
  fwd_txfm2d_c(input, output, stride, &cfg, txfm_buf);
}

void vp10_fwd_txfm2d_32x32_c(const int16_t *input, int32_t *output,
                             const int stride, int tx_type,
                             const int bd) {
  int32_t txfm_buf[32 * 32];
  TXFM_2D_FLIP_CFG cfg = vp10_get_fwd_txfm_cfg(tx_type, TX_32X32);
  (void)bd;
  fwd_txfm2d_c(input, output, stride, &cfg, txfm_buf);
}

void vp10_fwd_txfm2d_64x64_c(const int16_t *input, int32_t *output,
                             const int stride, int tx_type,
                             const int bd) {
  int32_t txfm_buf[64 * 64];
  TXFM_2D_FLIP_CFG cfg = vp10_get_fwd_txfm_64x64_cfg(tx_type);
  (void)bd;
  fwd_txfm2d_c(input, output, stride, &cfg, txfm_buf);
}

#if CONFIG_EXT_TX
static const TXFM_2D_CFG* fwd_txfm_cfg_ls[FLIPADST_ADST + 1][TX_SIZES] = {
    {&fwd_txfm_2d_cfg_dct_dct_4  , &fwd_txfm_2d_cfg_dct_dct_8,
     &fwd_txfm_2d_cfg_dct_dct_16  , &fwd_txfm_2d_cfg_dct_dct_32},
    {&fwd_txfm_2d_cfg_adst_dct_4 , &fwd_txfm_2d_cfg_adst_dct_8,
     &fwd_txfm_2d_cfg_adst_dct_16 , &fwd_txfm_2d_cfg_adst_dct_32},
    {&fwd_txfm_2d_cfg_dct_adst_4 , &fwd_txfm_2d_cfg_dct_adst_8,
     &fwd_txfm_2d_cfg_dct_adst_16 , &fwd_txfm_2d_cfg_dct_adst_32},
    {&fwd_txfm_2d_cfg_adst_adst_4, &fwd_txfm_2d_cfg_adst_adst_8,
     &fwd_txfm_2d_cfg_adst_adst_16, &fwd_txfm_2d_cfg_adst_adst_32},
    {&fwd_txfm_2d_cfg_adst_dct_4 , &fwd_txfm_2d_cfg_adst_dct_8,
     &fwd_txfm_2d_cfg_adst_dct_16 , &fwd_txfm_2d_cfg_adst_dct_32},
    {&fwd_txfm_2d_cfg_dct_adst_4 , &fwd_txfm_2d_cfg_dct_adst_8,
     &fwd_txfm_2d_cfg_dct_adst_16 , &fwd_txfm_2d_cfg_dct_adst_32},
    {&fwd_txfm_2d_cfg_adst_adst_4, &fwd_txfm_2d_cfg_adst_adst_8,
     &fwd_txfm_2d_cfg_adst_adst_16, &fwd_txfm_2d_cfg_adst_adst_32},
    {&fwd_txfm_2d_cfg_adst_adst_4, &fwd_txfm_2d_cfg_adst_adst_8,
     &fwd_txfm_2d_cfg_adst_adst_16, &fwd_txfm_2d_cfg_adst_adst_32},
    {&fwd_txfm_2d_cfg_adst_adst_4, &fwd_txfm_2d_cfg_adst_adst_8,
     &fwd_txfm_2d_cfg_adst_adst_16, &fwd_txfm_2d_cfg_adst_adst_32},
};
#else  // CONFIG_EXT_TX
static const TXFM_2D_CFG* fwd_txfm_cfg_ls[TX_TYPES][TX_SIZES] = {
    {&fwd_txfm_2d_cfg_dct_dct_4  , &fwd_txfm_2d_cfg_dct_dct_8,
     &fwd_txfm_2d_cfg_dct_dct_16  , &fwd_txfm_2d_cfg_dct_dct_32},
    {&fwd_txfm_2d_cfg_adst_dct_4 , &fwd_txfm_2d_cfg_adst_dct_8,
     &fwd_txfm_2d_cfg_adst_dct_16 , &fwd_txfm_2d_cfg_adst_dct_32},
    {&fwd_txfm_2d_cfg_dct_adst_4 , &fwd_txfm_2d_cfg_dct_adst_8,
     &fwd_txfm_2d_cfg_dct_adst_16 , &fwd_txfm_2d_cfg_dct_adst_32},
    {&fwd_txfm_2d_cfg_adst_adst_4, &fwd_txfm_2d_cfg_adst_adst_8,
     &fwd_txfm_2d_cfg_adst_adst_16, &fwd_txfm_2d_cfg_adst_adst_32},
};
#endif  // CONFIG_EXT_TX

TXFM_2D_FLIP_CFG vp10_get_fwd_txfm_cfg(int tx_type, int tx_size) {
  TXFM_2D_FLIP_CFG cfg;
  set_flip_cfg(tx_type, &cfg);
  cfg.cfg = fwd_txfm_cfg_ls[tx_type][tx_size];
  return cfg;
}

TXFM_2D_FLIP_CFG vp10_get_fwd_txfm_64x64_cfg(int tx_type) {
  TXFM_2D_FLIP_CFG cfg;
  switch (tx_type) {
    case DCT_DCT:
      cfg.cfg = &fwd_txfm_2d_cfg_dct_dct_64;
      cfg.ud_flip = 0;
      cfg.lr_flip = 0;
      break;
    case ADST_DCT:
    case DCT_ADST:
    case ADST_ADST:
    default:
      cfg.ud_flip = 0;
      cfg.lr_flip = 0;
      assert(0);
  }
  return cfg;
}
