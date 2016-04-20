/*
 *  Copyright (c) 2015 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include <assert.h>

#include "vp10/common/enums.h"
#include "vp10/common/vp10_fwd_txfm1d.h"
#include "vp10/common/vp10_fwd_txfm2d_cfg.h"
#include "vp10/common/vp10_txfm.h"

static inline TxfmFunc fwd_txfm_type_to_func(TXFM_TYPE txfm_type) {
  switch (txfm_type) {
    case TXFM_TYPE_DCT4:
      return vp10_fdct4_new;
      break;
    case TXFM_TYPE_DCT8:
      return vp10_fdct8_new;
      break;
    case TXFM_TYPE_DCT16:
      return vp10_fdct16_new;
      break;
    case TXFM_TYPE_DCT32:
      return vp10_fdct32_new;
      break;
    case TXFM_TYPE_DCT64:
      return vp10_fdct64_new;
      break;
    case TXFM_TYPE_ADST4:
      return vp10_fadst4_new;
      break;
    case TXFM_TYPE_ADST8:
      return vp10_fadst8_new;
      break;
    case TXFM_TYPE_ADST16:
      return vp10_fadst16_new;
      break;
    case TXFM_TYPE_ADST32:
      return vp10_fadst32_new;
      break;
    default:
      assert(0);
      return NULL;
  }
}

static inline void fwd_txfm2d_c(const int16_t *input, int32_t *output,
                                const int stride, const TXFM_2D_CFG *cfg,
                                int32_t *buf) {
  int i, j;
  const int txfm_size = cfg->txfm_size;
  const int8_t *shift = cfg->shift;
  const int8_t *stage_range_col = cfg->stage_range_col;
  const int8_t *stage_range_row = cfg->stage_range_row;
  const int8_t *cos_bit_col = cfg->cos_bit_col;
  const int8_t *cos_bit_row = cfg->cos_bit_row;
  const TxfmFunc txfm_func_col = fwd_txfm_type_to_func(cfg->txfm_type_col);
  const TxfmFunc txfm_func_row = fwd_txfm_type_to_func(cfg->txfm_type_row);

  // use output buffer as temp buffer
  int32_t* temp_in = output;
  int32_t* temp_out = output + txfm_size;

  // Columns
  for (i = 0; i < txfm_size; ++i) {
    for (j = 0; j < txfm_size; ++j)
      temp_in[j] = input[j * stride + i];
    round_shift_array(temp_in, txfm_size, -shift[0]);
    txfm_func_col(temp_in, temp_out, cos_bit_col, stage_range_col);
    round_shift_array(temp_out, txfm_size, -shift[1]);
    for (j = 0; j < txfm_size; ++j)
      buf[j * txfm_size + i] = temp_out[j];
  }

  // Rows
  for (i = 0; i < txfm_size; ++i) {
    txfm_func_row(buf + i * txfm_size, output + i * txfm_size, cos_bit_row,
                  stage_range_row);
    round_shift_array(output + i * txfm_size, txfm_size, -shift[2]);
  }
}

void vp10_fwd_txfm2d_4x4_c(const int16_t *input, int32_t *output,
                         const int stride, int tx_type,
                         const int bd) {
  int32_t txfm_buf[4 * 4];
  const TXFM_2D_CFG* cfg = vp10_get_txfm_4x4_cfg(tx_type);
  (void)bd;
  fwd_txfm2d_c(input, output, stride, cfg, txfm_buf);
}

void vp10_fwd_txfm2d_8x8_c(const int16_t *input, int32_t *output,
                         const int stride, int tx_type,
                         const int bd) {
  int32_t txfm_buf[8 * 8];
  const TXFM_2D_CFG* cfg = vp10_get_txfm_8x8_cfg(tx_type);
  (void)bd;
  fwd_txfm2d_c(input, output, stride, cfg, txfm_buf);
}

void vp10_fwd_txfm2d_16x16_c(const int16_t *input, int32_t *output,
                           const int stride, int tx_type,
                           const int bd) {
  int32_t txfm_buf[16 * 16];
  const TXFM_2D_CFG* cfg = vp10_get_txfm_16x16_cfg(tx_type);
  (void)bd;
  fwd_txfm2d_c(input, output, stride, cfg, txfm_buf);
}

void vp10_fwd_txfm2d_32x32_c(const int16_t *input, int32_t *output,
                           const int stride, int tx_type,
                           const int bd) {
  int32_t txfm_buf[32 * 32];
  const TXFM_2D_CFG* cfg = vp10_get_txfm_32x32_cfg(tx_type);
  (void)bd;
  fwd_txfm2d_c(input, output, stride, cfg, txfm_buf);
}

void vp10_fwd_txfm2d_64x64_c(const int16_t *input, int32_t *output,
                           const int stride, int tx_type,
                           const int bd) {
  int32_t txfm_buf[64 * 64];
  const TXFM_2D_CFG* cfg = vp10_get_txfm_64x64_cfg(tx_type);
  (void)bd;
  fwd_txfm2d_c(input, output, stride, cfg, txfm_buf);
}

const TXFM_2D_CFG* vp10_get_txfm_4x4_cfg(int tx_type) {
  const TXFM_2D_CFG* cfg = NULL;
  switch (tx_type) {
    case DCT_DCT:
      cfg = &fwd_txfm_2d_cfg_dct_dct_4;
      break;
    case ADST_DCT:
      cfg = &fwd_txfm_2d_cfg_adst_dct_4;
      break;
    case DCT_ADST:
      cfg = &fwd_txfm_2d_cfg_dct_adst_4;
      break;
    case ADST_ADST:
      cfg = &fwd_txfm_2d_cfg_adst_adst_4;
      break;
    default:
      assert(0);
  }
  return cfg;
}

const TXFM_2D_CFG* vp10_get_txfm_8x8_cfg(int tx_type) {
  const TXFM_2D_CFG* cfg = NULL;
  switch (tx_type) {
    case DCT_DCT:
      cfg = &fwd_txfm_2d_cfg_dct_dct_8;
      break;
    case ADST_DCT:
      cfg = &fwd_txfm_2d_cfg_adst_dct_8;
      break;
    case DCT_ADST:
      cfg = &fwd_txfm_2d_cfg_dct_adst_8;
      break;
    case ADST_ADST:
      cfg = &fwd_txfm_2d_cfg_adst_adst_8;
      break;
    default:
      assert(0);
  }
  return cfg;
}

const TXFM_2D_CFG* vp10_get_txfm_16x16_cfg(int tx_type) {
  const TXFM_2D_CFG* cfg = NULL;
  switch (tx_type) {
    case DCT_DCT:
      cfg = &fwd_txfm_2d_cfg_dct_dct_16;
      break;
    case ADST_DCT:
      cfg = &fwd_txfm_2d_cfg_adst_dct_16;
      break;
    case DCT_ADST:
      cfg = &fwd_txfm_2d_cfg_dct_adst_16;
      break;
    case ADST_ADST:
      cfg = &fwd_txfm_2d_cfg_adst_adst_16;
      break;
    default:
      assert(0);
  }
  return cfg;
}

const TXFM_2D_CFG* vp10_get_txfm_32x32_cfg(int tx_type) {
  const TXFM_2D_CFG* cfg = NULL;
  switch (tx_type) {
    case DCT_DCT:
      cfg = &fwd_txfm_2d_cfg_dct_dct_32;
      break;
    case ADST_DCT:
      cfg = &fwd_txfm_2d_cfg_adst_dct_32;
      break;
    case DCT_ADST:
      cfg = &fwd_txfm_2d_cfg_dct_adst_32;
      break;
    case ADST_ADST:
      cfg = &fwd_txfm_2d_cfg_adst_adst_32;
      break;
    default:
      assert(0);
  }
  return cfg;
}

const TXFM_2D_CFG* vp10_get_txfm_64x64_cfg(int tx_type) {
  const TXFM_2D_CFG* cfg = NULL;
  switch (tx_type) {
    case DCT_DCT:
      cfg = &fwd_txfm_2d_cfg_dct_dct_64;
      break;
    case ADST_DCT:
    case DCT_ADST:
    case ADST_ADST:
    default:
      assert(0);
  }
  return cfg;
}
