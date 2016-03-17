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

static INLINE void fwd_txfm2d_c(const int16_t *input, int32_t *output,
                                const int stride, const TXFM_2D_CFG *cfg,
                                int32_t *txfm_buf) {
  int i, j;
  const int txfm_size = cfg->txfm_size;
  const int8_t *shift = cfg->shift;
  const int8_t *stage_range_col = cfg->stage_range_col;
  const int8_t *stage_range_row = cfg->stage_range_row;
  const int8_t *cos_bit_col = cfg->cos_bit_col;
  const int8_t *cos_bit_row = cfg->cos_bit_row;
  const TxfmFunc txfm_func_col = cfg->txfm_func_col;
  const TxfmFunc txfm_func_row = cfg->txfm_func_row;

  // txfm_buf's length is  txfm_size * txfm_size + 2 * txfm_size
  // it is used for intermediate data buffering
  int32_t *temp_in = txfm_buf;
  int32_t *temp_out = temp_in + txfm_size;
  int32_t *buf = temp_out + txfm_size;

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
    for (j = 0; j < txfm_size; ++j)
      temp_in[j] = buf[j + i * txfm_size];
    txfm_func_row(temp_in, temp_out, cos_bit_row, stage_range_row);
    round_shift_array(temp_out, txfm_size, -shift[2]);
    for (j = 0; j < txfm_size; ++j)
      output[j + i * txfm_size] = (int32_t)temp_out[j];
  }
}

void vp10_fwd_txfm2d_4x4(const int16_t *input, int32_t *output,
                         const int stride, const TXFM_2D_CFG *cfg,
                         const int bd) {
  int txfm_buf[4 * 4 + 4 + 4];
  (void)bd;
  fwd_txfm2d_c(input, output, stride, cfg, txfm_buf);
}

void vp10_fwd_txfm2d_8x8(const int16_t *input, int32_t *output,
                         const int stride, const TXFM_2D_CFG *cfg,
                         const int bd) {
  int txfm_buf[8 * 8 + 8 + 8];
  (void)bd;
  fwd_txfm2d_c(input, output, stride, cfg, txfm_buf);
}

void vp10_fwd_txfm2d_16x16(const int16_t *input, int32_t *output,
                           const int stride, const TXFM_2D_CFG *cfg,
                           const int bd) {
  int txfm_buf[16 * 16 + 16 + 16];
  (void)bd;
  fwd_txfm2d_c(input, output, stride, cfg, txfm_buf);
}

void vp10_fwd_txfm2d_32x32(const int16_t *input, int32_t *output,
                           const int stride, const TXFM_2D_CFG *cfg,
                           const int bd) {
  int txfm_buf[32 * 32 + 32 + 32];
  (void)bd;
  fwd_txfm2d_c(input, output, stride, cfg, txfm_buf);
}

void vp10_fwd_txfm2d_64x64(const int16_t *input, int32_t *output,
                           const int stride, const TXFM_2D_CFG *cfg,
                           const int bd) {
  int txfm_buf[64 * 64 + 64 + 64];
  (void)bd;
  fwd_txfm2d_c(input, output, stride, cfg, txfm_buf);
}
