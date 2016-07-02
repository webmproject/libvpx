/*
 *  Copyright (c) 2016 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include "./vp10_rtcd.h"
#include "vp10/common/enums.h"
#include "vp10/common/vp10_txfm.h"
#include "vp10/common/x86/vp10_txfm1d_sse4.h"

static INLINE void int16_array_with_stride_to_int32_array_without_stride(
    const int16_t *input, int stride, int32_t *output, int txfm1d_size) {
  int r, c;
  for (r = 0; r < txfm1d_size; r++) {
    for (c = 0; c < txfm1d_size; c++) {
      output[r * txfm1d_size + c] = (int32_t)input[r * stride + c];
    }
  }
}

typedef void (*TxfmFuncSSE2)(const __m128i *input, __m128i *output,
                             const int8_t *cos_bit, const int8_t *stage_range);

static INLINE TxfmFuncSSE2 fwd_txfm_type_to_func(TXFM_TYPE txfm_type) {
  switch (txfm_type) {
    case TXFM_TYPE_DCT4:
      return vp10_fdct4_new_sse4_1;
      break;
    case TXFM_TYPE_DCT8:
      return vp10_fdct8_new_sse4_1;
      break;
    case TXFM_TYPE_DCT16:
      return vp10_fdct16_new_sse4_1;
      break;
    case TXFM_TYPE_DCT32:
      return vp10_fdct32_new_sse4_1;
      break;
    case TXFM_TYPE_DCT64:
      return vp10_fdct64_new_sse4_1;
      break;
    case TXFM_TYPE_ADST4:
      return vp10_fadst4_new_sse4_1;
      break;
    case TXFM_TYPE_ADST8:
      return vp10_fadst8_new_sse4_1;
      break;
    case TXFM_TYPE_ADST16:
      return vp10_fadst16_new_sse4_1;
      break;
    case TXFM_TYPE_ADST32:
      return vp10_fadst32_new_sse4_1;
      break;
    default:
      assert(0);
  }
  return NULL;
}

static INLINE void fwd_txfm2d_sse4_1(const int16_t *input, int32_t *output,
                                     const int stride, const TXFM_2D_CFG *cfg,
                                     int32_t *txfm_buf) {
  const int txfm_size = cfg->txfm_size;
  const int8_t *shift = cfg->shift;
  const int8_t *stage_range_col = cfg->stage_range_col;
  const int8_t *stage_range_row = cfg->stage_range_row;
  const int8_t *cos_bit_col = cfg->cos_bit_col;
  const int8_t *cos_bit_row = cfg->cos_bit_row;
  const TxfmFuncSSE2 txfm_func_col = fwd_txfm_type_to_func(cfg->txfm_type_col);
  const TxfmFuncSSE2 txfm_func_row = fwd_txfm_type_to_func(cfg->txfm_type_row);

  __m128i *buf_128 = (__m128i *)txfm_buf;
  __m128i *out_128 = (__m128i *)output;
  int num_per_128 = 4;
  int txfm2d_size_128 = txfm_size * txfm_size / num_per_128;

  int16_array_with_stride_to_int32_array_without_stride(input, stride, txfm_buf,
                                                        txfm_size);
  round_shift_array_32_sse4_1(buf_128, out_128, txfm2d_size_128, -shift[0]);
  txfm_func_col(out_128, buf_128, cos_bit_col, stage_range_col);
  round_shift_array_32_sse4_1(buf_128, out_128, txfm2d_size_128, -shift[1]);
  transpose_32(txfm_size, out_128, buf_128);
  txfm_func_row(buf_128, out_128, cos_bit_row, stage_range_row);
  round_shift_array_32_sse4_1(out_128, buf_128, txfm2d_size_128, -shift[2]);
  transpose_32(txfm_size, buf_128, out_128);
}

void vp10_fwd_txfm2d_32x32_sse4_1(const int16_t *input, int32_t *output,
                                  const int stride, int tx_type,
                                  const int bd) {
  int32_t txfm_buf[1024];
  TXFM_2D_FLIP_CFG cfg = vp10_get_fwd_txfm_cfg(tx_type, TX_32X32);
  (void)bd;
  fwd_txfm2d_sse4_1(input, output, stride, cfg.cfg, txfm_buf);
}

void vp10_fwd_txfm2d_64x64_sse4_1(const int16_t *input, int32_t *output,
                                  const int stride, int tx_type,
                                  const int bd) {
  int32_t txfm_buf[4096];
  TXFM_2D_FLIP_CFG cfg = vp10_get_fwd_txfm_64x64_cfg(tx_type);
  (void)bd;
  fwd_txfm2d_sse4_1(input, output, stride, cfg.cfg, txfm_buf);
}
