/*
 *  Copyright (c) 2015 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include "vp10/common/vp10_txfm.h"
#include "vp10/common/vp10_inv_txfm1d.h"

static INLINE TxfmFunc inv_txfm_type_to_func(TXFM_TYPE txfm_type) {
  switch (txfm_type) {
    case TXFM_TYPE_DCT4:
      return vp10_idct4_new;
      break;
    case TXFM_TYPE_DCT8:
      return vp10_idct8_new;
      break;
    case TXFM_TYPE_DCT16:
      return vp10_idct16_new;
      break;
    case TXFM_TYPE_DCT32:
      return vp10_idct32_new;
      break;
    case TXFM_TYPE_DCT64:
      return vp10_idct64_new;
      break;
    case TXFM_TYPE_ADST4:
      return vp10_iadst4_new;
      break;
    case TXFM_TYPE_ADST8:
      return vp10_iadst8_new;
      break;
    case TXFM_TYPE_ADST16:
      return vp10_iadst16_new;
      break;
    case TXFM_TYPE_ADST32:
      return vp10_iadst32_new;
      break;
    default:
      assert(0);
      return NULL;
  }
}

static INLINE void inv_txfm2d_add_c(const int32_t *input, int16_t *output,
                                    int stride, const TXFM_2D_CFG *cfg,
                                    int32_t *txfm_buf) {
  const int txfm_size = cfg->txfm_size;
  const int8_t *shift = cfg->shift;
  const int8_t *stage_range_col = cfg->stage_range_col;
  const int8_t *stage_range_row = cfg->stage_range_row;
  const int8_t *cos_bit_col = cfg->cos_bit_col;
  const int8_t *cos_bit_row = cfg->cos_bit_row;
  const TxfmFunc txfm_func_col = inv_txfm_type_to_func(cfg->txfm_type_col);
  const TxfmFunc txfm_func_row = inv_txfm_type_to_func(cfg->txfm_type_row);

  // txfm_buf's length is  txfm_size * txfm_size + 2 * txfm_size
  // it is used for intermediate data buffering
  int32_t *temp_in = txfm_buf;
  int32_t *temp_out = temp_in + txfm_size;
  int32_t *buf = temp_out + txfm_size;
  int32_t *buf_ptr = buf;
  int i, j;

  // Rows
  for (i = 0; i < txfm_size; ++i) {
    txfm_func_row(input, buf_ptr, cos_bit_row, stage_range_row);
    round_shift_array(buf_ptr, txfm_size, -shift[0]);
    input += txfm_size;
    buf_ptr += txfm_size;
  }

  // Columns
  for (i = 0; i < txfm_size; ++i) {
    for (j = 0; j < txfm_size; ++j)
      temp_in[j] = buf[j * txfm_size + i];
    txfm_func_col(temp_in, temp_out, cos_bit_col, stage_range_col);
    round_shift_array(temp_out, txfm_size, -shift[1]);
    for (j = 0; j < txfm_size; ++j)
      output[j * stride + i] += temp_out[j];
  }
}

void vp10_inv_txfm2d_add_4x4_c(const int32_t *input, uint16_t *output,
                             const int stride, const TXFM_2D_CFG *cfg,
                             const int bd) {
  int txfm_buf[4 * 4 + 4 + 4];
  // output contains the prediction signal which is always positive and smaller
  // than (1 << bd) - 1
  // since bd < 16-1, therefore we can treat the uint16_t* output buffer as an
  // int16_t*
  inv_txfm2d_add_c(input, (int16_t *)output, stride, cfg, txfm_buf);
  clamp_block((int16_t *)output, 4, stride, 0, (1 << bd) - 1);
}

void vp10_inv_txfm2d_add_8x8_c(const int32_t *input, uint16_t *output,
                             const int stride, const TXFM_2D_CFG *cfg,
                             const int bd) {
  int txfm_buf[8 * 8 + 8 + 8];
  // output contains the prediction signal which is always positive and smaller
  // than (1 << bd) - 1
  // since bd < 16-1, therefore we can treat the uint16_t* output buffer as an
  // int16_t*
  inv_txfm2d_add_c(input, (int16_t *)output, stride, cfg, txfm_buf);
  clamp_block((int16_t *)output, 8, stride, 0, (1 << bd) - 1);
}

void vp10_inv_txfm2d_add_16x16_c(const int32_t *input, uint16_t *output,
                               const int stride, const TXFM_2D_CFG *cfg,
                               const int bd) {
  int txfm_buf[16 * 16 + 16 + 16];
  // output contains the prediction signal which is always positive and smaller
  // than (1 << bd) - 1
  // since bd < 16-1, therefore we can treat the uint16_t* output buffer as an
  // int16_t*
  inv_txfm2d_add_c(input, (int16_t *)output, stride, cfg, txfm_buf);
  clamp_block((int16_t *)output, 16, stride, 0, (1 << bd) - 1);
}

void vp10_inv_txfm2d_add_32x32_c(const int32_t *input, uint16_t *output,
                               const int stride, const TXFM_2D_CFG *cfg,
                               const int bd) {
  int txfm_buf[32 * 32 + 32 + 32];
  // output contains the prediction signal which is always positive and smaller
  // than (1 << bd) - 1
  // since bd < 16-1, therefore we can treat the uint16_t* output buffer as an
  // int16_t*
  inv_txfm2d_add_c(input, (int16_t *)output, stride, cfg, txfm_buf);
  clamp_block((int16_t *)output, 32, stride, 0, (1 << bd) - 1);
}

void vp10_inv_txfm2d_add_64x64_c(const int32_t *input, uint16_t *output,
                               const int stride, const TXFM_2D_CFG *cfg,
                               const int bd) {
  int txfm_buf[64 * 64 + 64 + 64];
  // output contains the prediction signal which is always positive and smaller
  // than (1 << bd) - 1
  // since bd < 16-1, therefore we can treat the uint16_t* output buffer as an
  // int16_t*
  inv_txfm2d_add_c(input, (int16_t *)output, stride, cfg, txfm_buf);
  clamp_block((int16_t *)output, 64, stride, 0, (1 << bd) - 1);
}
