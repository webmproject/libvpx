/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include <math.h>

#include "./vp10_rtcd.h"
#include "./vpx_dsp_rtcd.h"
#include "vp10/common/blockd.h"
#include "vp10/common/enums.h"
#include "vp10/common/idct.h"
#include "vp10/common/vp10_inv_txfm2d_cfg.h"
#include "vpx_dsp/inv_txfm.h"
#include "vpx_ports/mem.h"

int get_tx_scale(const MACROBLOCKD *const xd, const TX_TYPE tx_type,
                 const TX_SIZE tx_size) {
  (void) tx_type;
#if CONFIG_VP9_HIGHBITDEPTH
  if (xd->cur_buf->flags & YV12_FLAG_HIGHBITDEPTH) {
    return tx_size == TX_32X32;
  }
#else
  (void)xd;
#endif
  return tx_size == TX_32X32;
}

#if CONFIG_EXT_TX
static void iidtx4_c(const tran_low_t *input, tran_low_t *output) {
  int i;
  for (i = 0; i < 4; ++i)
    output[i] = (tran_low_t)dct_const_round_shift(input[i] * Sqrt2);
}

static void iidtx8_c(const tran_low_t *input, tran_low_t *output) {
  int i;
  for (i = 0; i < 8; ++i)
    output[i] = input[i] * 2;
}

static void iidtx16_c(const tran_low_t *input, tran_low_t *output) {
  int i;
  for (i = 0; i < 16; ++i)
    output[i] = (tran_low_t)dct_const_round_shift(input[i] * 2 * Sqrt2);
}

static void iidtx32_c(const tran_low_t *input, tran_low_t *output) {
  int i;
  for (i = 0; i < 32; ++i)
    output[i] = input[i] * 4;
}

// For use in lieu of DST
static void ihalfright32_c(const tran_low_t *input, tran_low_t *output) {
  int i;
  tran_low_t inputhalf[16];
  for (i = 0; i < 16; ++i) {
    output[i] = input[16 + i] * 4;
  }
  // Multiply input by sqrt(2)
  for (i = 0; i < 16; ++i) {
    inputhalf[i] = (tran_low_t)dct_const_round_shift(input[i] * Sqrt2);
  }
  idct16_c(inputhalf, output + 16);
  // Note overall scaling factor is 4 times orthogonal
}

#if CONFIG_VP9_HIGHBITDEPTH
static void highbd_iidtx4_c(const tran_low_t *input, tran_low_t *output,
                            int bd) {
  int i;
  for (i = 0; i < 4; ++i)
    output[i] = HIGHBD_WRAPLOW(
        highbd_dct_const_round_shift(input[i] * Sqrt2), bd);
}

static void highbd_iidtx8_c(const tran_low_t *input, tran_low_t *output,
                            int bd) {
  int i;
  (void) bd;
  for (i = 0; i < 8; ++i)
    output[i] = input[i] * 2;
}

static void highbd_iidtx16_c(const tran_low_t *input, tran_low_t *output,
                            int bd) {
  int i;
  for (i = 0; i < 16; ++i)
    output[i] = HIGHBD_WRAPLOW(
        highbd_dct_const_round_shift(input[i] * 2 * Sqrt2), bd);
}

static void highbd_iidtx32_c(const tran_low_t *input, tran_low_t *output,
                             int bd) {
  int i;
  (void) bd;
  for (i = 0; i < 32; ++i)
    output[i] = input[i] * 4;
}

static void highbd_ihalfright32_c(const tran_low_t *input, tran_low_t *output,
                                  int bd) {
  int i;
  tran_low_t inputhalf[16];
  for (i = 0; i < 16; ++i) {
    output[i] = input[16 + i] * 4;
  }
  // Multiply input by sqrt(2)
  for (i = 0; i < 16; ++i) {
    inputhalf[i] = HIGHBD_WRAPLOW(
        highbd_dct_const_round_shift(input[i] * Sqrt2), bd);
  }
  vpx_highbd_idct16_c(inputhalf, output + 16, bd);
  // Note overall scaling factor is 4 times orthogonal
}
#endif  // CONFIG_VP9_HIGHBITDEPTH

// Inverse identity transform and add.
static void inv_idtx_add_c(const tran_low_t *input, uint8_t *dest, int stride,
                           int bs, int tx_type) {
  int r, c;
  const int shift = bs < 32 ? 3 : 2;
  if (tx_type == IDTX) {
    for (r = 0; r < bs; ++r) {
      for (c = 0; c < bs; ++c)
        dest[c] = clip_pixel_add(dest[c], input[c] >> shift);
      dest += stride;
      input += bs;
    }
  }
}

#define FLIPUD_PTR(dest, stride, size) do {     \
  (dest) = (dest) + ((size) - 1) * (stride);  \
  (stride) = - (stride);                      \
} while (0)

static void maybe_flip_strides(uint8_t **dst, int *dstride,
                               tran_low_t **src, int *sstride,
                               int tx_type, int size) {
  // Note that the transpose of src will be added to dst. In order to LR
  // flip the addends (in dst coordinates), we UD flip the src. To UD flip
  // the addends, we UD flip the dst.
  switch (tx_type) {
    case DCT_DCT:
    case ADST_DCT:
    case DCT_ADST:
    case ADST_ADST:
    case IDTX:
    case V_DCT:
    case H_DCT:
    case V_ADST:
    case H_ADST:
      break;
    case FLIPADST_DCT:
    case FLIPADST_ADST:
    case V_FLIPADST:
      // flip UD
      FLIPUD_PTR(*dst, *dstride, size);
      break;
    case DCT_FLIPADST:
    case ADST_FLIPADST:
    case H_FLIPADST:
      // flip LR
      FLIPUD_PTR(*src, *sstride, size);
      break;
    case FLIPADST_FLIPADST:
      // flip UD
      FLIPUD_PTR(*dst, *dstride, size);
      // flip LR
      FLIPUD_PTR(*src, *sstride, size);
      break;
    default:
      assert(0);
      break;
  }
}

#if CONFIG_VP9_HIGHBITDEPTH
void highbd_idst4_c(const tran_low_t *input, tran_low_t *output, int bd) {
  tran_low_t step[4];
  tran_high_t temp1, temp2;
  (void) bd;
  // stage 1
  temp1 = (input[3] + input[1]) * cospi_16_64;
  temp2 = (input[3] - input[1]) * cospi_16_64;
  step[0] = HIGHBD_WRAPLOW(dct_const_round_shift(temp1), bd);
  step[1] = HIGHBD_WRAPLOW(dct_const_round_shift(temp2), bd);
  temp1 = input[2] * cospi_24_64 - input[0] * cospi_8_64;
  temp2 = input[2] * cospi_8_64 + input[0] * cospi_24_64;
  step[2] = HIGHBD_WRAPLOW(dct_const_round_shift(temp1), bd);
  step[3] = HIGHBD_WRAPLOW(dct_const_round_shift(temp2), bd);

  // stage 2
  output[0] = HIGHBD_WRAPLOW(step[0] + step[3], bd);
  output[1] = HIGHBD_WRAPLOW(-step[1] - step[2], bd);
  output[2] = HIGHBD_WRAPLOW(step[1] - step[2], bd);
  output[3] = HIGHBD_WRAPLOW(step[3] - step[0], bd);
}

void highbd_idst8_c(const tran_low_t *input, tran_low_t *output, int bd) {
  tran_low_t step1[8], step2[8];
  tran_high_t temp1, temp2;
  (void) bd;
  // stage 1
  step1[0] = input[7];
  step1[2] = input[3];
  step1[1] = input[5];
  step1[3] = input[1];
  temp1 = input[6] * cospi_28_64 - input[0] * cospi_4_64;
  temp2 = input[6] * cospi_4_64 + input[0] * cospi_28_64;
  step1[4] = HIGHBD_WRAPLOW(dct_const_round_shift(temp1), bd);
  step1[7] = HIGHBD_WRAPLOW(dct_const_round_shift(temp2), bd);
  temp1 = input[2] * cospi_12_64 - input[4] * cospi_20_64;
  temp2 = input[2] * cospi_20_64 + input[4] * cospi_12_64;
  step1[5] = HIGHBD_WRAPLOW(dct_const_round_shift(temp1), bd);
  step1[6] = HIGHBD_WRAPLOW(dct_const_round_shift(temp2), bd);

  // stage 2
  temp1 = (step1[0] + step1[2]) * cospi_16_64;
  temp2 = (step1[0] - step1[2]) * cospi_16_64;
  step2[0] = HIGHBD_WRAPLOW(dct_const_round_shift(temp1), bd);
  step2[1] = HIGHBD_WRAPLOW(dct_const_round_shift(temp2), bd);
  temp1 = step1[1] * cospi_24_64 - step1[3] * cospi_8_64;
  temp2 = step1[1] * cospi_8_64 + step1[3] * cospi_24_64;
  step2[2] = HIGHBD_WRAPLOW(dct_const_round_shift(temp1), bd);
  step2[3] = HIGHBD_WRAPLOW(dct_const_round_shift(temp2), bd);
  step2[4] = HIGHBD_WRAPLOW(step1[4] + step1[5], bd);
  step2[5] = HIGHBD_WRAPLOW(step1[4] - step1[5], bd);
  step2[6] = HIGHBD_WRAPLOW(-step1[6] + step1[7], bd);
  step2[7] = HIGHBD_WRAPLOW(step1[6] + step1[7], bd);

  // stage 3
  step1[0] = HIGHBD_WRAPLOW(step2[0] + step2[3], bd);
  step1[1] = HIGHBD_WRAPLOW(step2[1] + step2[2], bd);
  step1[2] = HIGHBD_WRAPLOW(step2[1] - step2[2], bd);
  step1[3] = HIGHBD_WRAPLOW(step2[0] - step2[3], bd);
  step1[4] = step2[4];
  temp1 = (step2[6] - step2[5]) * cospi_16_64;
  temp2 = (step2[5] + step2[6]) * cospi_16_64;
  step1[5] = HIGHBD_WRAPLOW(dct_const_round_shift(temp1), bd);
  step1[6] = HIGHBD_WRAPLOW(dct_const_round_shift(temp2), bd);
  step1[7] = step2[7];

  // stage 4
  output[0] = HIGHBD_WRAPLOW(step1[0] + step1[7], bd);
  output[1] = HIGHBD_WRAPLOW(-step1[1] - step1[6], bd);
  output[2] = HIGHBD_WRAPLOW(step1[2] + step1[5], bd);
  output[3] = HIGHBD_WRAPLOW(-step1[3] - step1[4], bd);
  output[4] = HIGHBD_WRAPLOW(step1[3] - step1[4], bd);
  output[5] = HIGHBD_WRAPLOW(-step1[2] + step1[5], bd);
  output[6] = HIGHBD_WRAPLOW(step1[1] - step1[6], bd);
  output[7] = HIGHBD_WRAPLOW(-step1[0] + step1[7], bd);
}

void highbd_idst16_c(const tran_low_t *input, tran_low_t *output, int bd) {
  // vp9_highbd_igentx16(input, output, bd, Tx16);
  tran_low_t step1[16], step2[16];
  tran_high_t temp1, temp2;
  (void) bd;

  // stage 1
  step1[0] = input[15];
  step1[1] = input[7];
  step1[2] = input[11];
  step1[3] = input[3];
  step1[4] = input[13];
  step1[5] = input[5];
  step1[6] = input[9];
  step1[7] = input[1];
  step1[8] = input[14];
  step1[9] = input[6];
  step1[10] = input[10];
  step1[11] = input[2];
  step1[12] = input[12];
  step1[13] = input[4];
  step1[14] = input[8];
  step1[15] = input[0];

  // stage 2
  step2[0] = step1[0];
  step2[1] = step1[1];
  step2[2] = step1[2];
  step2[3] = step1[3];
  step2[4] = step1[4];
  step2[5] = step1[5];
  step2[6] = step1[6];
  step2[7] = step1[7];

  temp1 = step1[8] * cospi_30_64 - step1[15] * cospi_2_64;
  temp2 = step1[8] * cospi_2_64 + step1[15] * cospi_30_64;
  step2[8] = HIGHBD_WRAPLOW(dct_const_round_shift(temp1), bd);
  step2[15] = HIGHBD_WRAPLOW(dct_const_round_shift(temp2), bd);

  temp1 = step1[9] * cospi_14_64 - step1[14] * cospi_18_64;
  temp2 = step1[9] * cospi_18_64 + step1[14] * cospi_14_64;
  step2[9] = HIGHBD_WRAPLOW(dct_const_round_shift(temp1), bd);
  step2[14] = HIGHBD_WRAPLOW(dct_const_round_shift(temp2), bd);

  temp1 = step1[10] * cospi_22_64 - step1[13] * cospi_10_64;
  temp2 = step1[10] * cospi_10_64 + step1[13] * cospi_22_64;
  step2[10] = HIGHBD_WRAPLOW(dct_const_round_shift(temp1), bd);
  step2[13] = HIGHBD_WRAPLOW(dct_const_round_shift(temp2), bd);

  temp1 = step1[11] * cospi_6_64 - step1[12] * cospi_26_64;
  temp2 = step1[11] * cospi_26_64 + step1[12] * cospi_6_64;
  step2[11] = HIGHBD_WRAPLOW(dct_const_round_shift(temp1), bd);
  step2[12] = HIGHBD_WRAPLOW(dct_const_round_shift(temp2), bd);

  // stage 3
  step1[0] = step2[0];
  step1[1] = step2[1];
  step1[2] = step2[2];
  step1[3] = step2[3];

  temp1 = step2[4] * cospi_28_64 - step2[7] * cospi_4_64;
  temp2 = step2[4] * cospi_4_64 + step2[7] * cospi_28_64;
  step1[4] = HIGHBD_WRAPLOW(dct_const_round_shift(temp1), bd);
  step1[7] = HIGHBD_WRAPLOW(dct_const_round_shift(temp2), bd);
  temp1 = step2[5] * cospi_12_64 - step2[6] * cospi_20_64;
  temp2 = step2[5] * cospi_20_64 + step2[6] * cospi_12_64;
  step1[5] = HIGHBD_WRAPLOW(dct_const_round_shift(temp1), bd);
  step1[6] = HIGHBD_WRAPLOW(dct_const_round_shift(temp2), bd);

  step1[8] = HIGHBD_WRAPLOW(step2[8] + step2[9], bd);
  step1[9] = HIGHBD_WRAPLOW(step2[8] - step2[9], bd);
  step1[10] = HIGHBD_WRAPLOW(-step2[10] + step2[11], bd);
  step1[11] = HIGHBD_WRAPLOW(step2[10] + step2[11], bd);
  step1[12] = HIGHBD_WRAPLOW(step2[12] + step2[13], bd);
  step1[13] = HIGHBD_WRAPLOW(step2[12] - step2[13], bd);
  step1[14] = HIGHBD_WRAPLOW(-step2[14] + step2[15], bd);
  step1[15] = HIGHBD_WRAPLOW(step2[14] + step2[15], bd);

  // stage 4
  temp1 = (step1[0] + step1[1]) * cospi_16_64;
  temp2 = (step1[0] - step1[1]) * cospi_16_64;
  step2[0] = HIGHBD_WRAPLOW(dct_const_round_shift(temp1), bd);
  step2[1] = HIGHBD_WRAPLOW(dct_const_round_shift(temp2), bd);
  temp1 = step1[2] * cospi_24_64 - step1[3] * cospi_8_64;
  temp2 = step1[2] * cospi_8_64 + step1[3] * cospi_24_64;
  step2[2] = HIGHBD_WRAPLOW(dct_const_round_shift(temp1), bd);
  step2[3] = HIGHBD_WRAPLOW(dct_const_round_shift(temp2), bd);
  step2[4] = HIGHBD_WRAPLOW(step1[4] + step1[5], bd);
  step2[5] = HIGHBD_WRAPLOW(step1[4] - step1[5], bd);
  step2[6] = HIGHBD_WRAPLOW(-step1[6] + step1[7], bd);
  step2[7] = HIGHBD_WRAPLOW(step1[6] + step1[7], bd);

  step2[8] = step1[8];
  step2[15] = step1[15];
  temp1 = -step1[9] * cospi_8_64 + step1[14] * cospi_24_64;
  temp2 = step1[9] * cospi_24_64 + step1[14] * cospi_8_64;
  step2[9] = HIGHBD_WRAPLOW(dct_const_round_shift(temp1), bd);
  step2[14] = HIGHBD_WRAPLOW(dct_const_round_shift(temp2), bd);
  temp1 = -step1[10] * cospi_24_64 - step1[13] * cospi_8_64;
  temp2 = -step1[10] * cospi_8_64 + step1[13] * cospi_24_64;
  step2[10] = HIGHBD_WRAPLOW(dct_const_round_shift(temp1), bd);
  step2[13] = HIGHBD_WRAPLOW(dct_const_round_shift(temp2), bd);
  step2[11] = step1[11];
  step2[12] = step1[12];

  // stage 5
  step1[0] = HIGHBD_WRAPLOW(step2[0] + step2[3], bd);
  step1[1] = HIGHBD_WRAPLOW(step2[1] + step2[2], bd);
  step1[2] = HIGHBD_WRAPLOW(step2[1] - step2[2], bd);
  step1[3] = HIGHBD_WRAPLOW(step2[0] - step2[3], bd);
  step1[4] = step2[4];
  temp1 = (step2[6] - step2[5]) * cospi_16_64;
  temp2 = (step2[5] + step2[6]) * cospi_16_64;
  step1[5] = HIGHBD_WRAPLOW(dct_const_round_shift(temp1), bd);
  step1[6] = HIGHBD_WRAPLOW(dct_const_round_shift(temp2), bd);
  step1[7] = step2[7];

  step1[8] = HIGHBD_WRAPLOW(step2[8] + step2[11], bd);
  step1[9] = HIGHBD_WRAPLOW(step2[9] + step2[10], bd);
  step1[10] = HIGHBD_WRAPLOW(step2[9] - step2[10], bd);
  step1[11] = HIGHBD_WRAPLOW(step2[8] - step2[11], bd);
  step1[12] = HIGHBD_WRAPLOW(-step2[12] + step2[15], bd);
  step1[13] = HIGHBD_WRAPLOW(-step2[13] + step2[14], bd);
  step1[14] = HIGHBD_WRAPLOW(step2[13] + step2[14], bd);
  step1[15] = HIGHBD_WRAPLOW(step2[12] + step2[15], bd);

  // stage 6
  step2[0] = HIGHBD_WRAPLOW(step1[0] + step1[7], bd);
  step2[1] = HIGHBD_WRAPLOW(step1[1] + step1[6], bd);
  step2[2] = HIGHBD_WRAPLOW(step1[2] + step1[5], bd);
  step2[3] = HIGHBD_WRAPLOW(step1[3] + step1[4], bd);
  step2[4] = HIGHBD_WRAPLOW(step1[3] - step1[4], bd);
  step2[5] = HIGHBD_WRAPLOW(step1[2] - step1[5], bd);
  step2[6] = HIGHBD_WRAPLOW(step1[1] - step1[6], bd);
  step2[7] = HIGHBD_WRAPLOW(step1[0] - step1[7], bd);
  step2[8] = step1[8];
  step2[9] = step1[9];
  temp1 = (-step1[10] + step1[13]) * cospi_16_64;
  temp2 = (step1[10] + step1[13]) * cospi_16_64;
  step2[10] = HIGHBD_WRAPLOW(dct_const_round_shift(temp1), bd);
  step2[13] = HIGHBD_WRAPLOW(dct_const_round_shift(temp2), bd);
  temp1 = (-step1[11] + step1[12]) * cospi_16_64;
  temp2 = (step1[11] + step1[12]) * cospi_16_64;
  step2[11] = HIGHBD_WRAPLOW(dct_const_round_shift(temp1), bd);
  step2[12] = HIGHBD_WRAPLOW(dct_const_round_shift(temp2), bd);
  step2[14] = step1[14];
  step2[15] = step1[15];

  // stage 7
  output[0] = HIGHBD_WRAPLOW(step2[0] + step2[15], bd);
  output[1] = HIGHBD_WRAPLOW(-step2[1] - step2[14], bd);
  output[2] = HIGHBD_WRAPLOW(step2[2] + step2[13], bd);
  output[3] = HIGHBD_WRAPLOW(-step2[3] - step2[12], bd);
  output[4] = HIGHBD_WRAPLOW(step2[4] + step2[11], bd);
  output[5] = HIGHBD_WRAPLOW(-step2[5] - step2[10], bd);
  output[6] = HIGHBD_WRAPLOW(step2[6] + step2[9], bd);
  output[7] = HIGHBD_WRAPLOW(-step2[7] - step2[8], bd);
  output[8] = HIGHBD_WRAPLOW(step2[7] - step2[8], bd);
  output[9] = HIGHBD_WRAPLOW(-step2[6] + step2[9], bd);
  output[10] = HIGHBD_WRAPLOW(step2[5] - step2[10], bd);
  output[11] = HIGHBD_WRAPLOW(-step2[4] + step2[11], bd);
  output[12] = HIGHBD_WRAPLOW(step2[3] - step2[12], bd);
  output[13] = HIGHBD_WRAPLOW(-step2[2] + step2[13], bd);
  output[14] = HIGHBD_WRAPLOW(step2[1] - step2[14], bd);
  output[15] = HIGHBD_WRAPLOW(-step2[0] + step2[15], bd);
}

static void highbd_inv_idtx_add_c(const tran_low_t *input, uint8_t *dest8,
                                  int stride, int bs, int tx_type, int bd) {
  int r, c;
  const int shift = bs < 32 ? 3 : 2;
  uint16_t *dest = CONVERT_TO_SHORTPTR(dest8);

  if (tx_type == IDTX) {
    for (r = 0; r < bs; ++r) {
      for (c = 0; c < bs; ++c)
        dest[c] = highbd_clip_pixel_add(dest[c], input[c] >> shift, bd);
      dest += stride;
      input += bs;
    }
  }
}

static void maybe_flip_strides16(uint16_t **dst, int *dstride,
                                 tran_low_t **src, int *sstride,
                                 int tx_type, int size) {
  // Note that the transpose of src will be added to dst. In order to LR
  // flip the addends (in dst coordinates), we UD flip the src. To UD flip
  // the addends, we UD flip the dst.
  switch (tx_type) {
    case DCT_DCT:
    case ADST_DCT:
    case DCT_ADST:
    case ADST_ADST:
    case IDTX:
    case V_DCT:
    case H_DCT:
    case V_ADST:
    case H_ADST:
      break;
    case FLIPADST_DCT:
    case FLIPADST_ADST:
    case V_FLIPADST:
      // flip UD
      FLIPUD_PTR(*dst, *dstride, size);
      break;
    case DCT_FLIPADST:
    case ADST_FLIPADST:
    case H_FLIPADST:
      // flip LR
      FLIPUD_PTR(*src, *sstride, size);
      break;
    case FLIPADST_FLIPADST:
      // flip UD
      FLIPUD_PTR(*dst, *dstride, size);
      // flip LR
      FLIPUD_PTR(*src, *sstride, size);
      break;
    default:
      assert(0);
      break;
  }
}
#endif  // CONFIG_VP9_HIGHBITDEPTH
#endif  // CONFIG_EXT_TX

void vp10_iht4x4_16_add_c(const tran_low_t *input, uint8_t *dest, int stride,
                          int tx_type) {
  static const transform_2d IHT_4[] = {
    { idct4_c,  idct4_c  },  // DCT_DCT
    { iadst4_c, idct4_c  },  // ADST_DCT
    { idct4_c,  iadst4_c },  // DCT_ADST
    { iadst4_c, iadst4_c },  // ADST_ADST
#if CONFIG_EXT_TX
    { iadst4_c, idct4_c  },  // FLIPADST_DCT
    { idct4_c,  iadst4_c },  // DCT_FLIPADST
    { iadst4_c, iadst4_c },  // FLIPADST_FLIPADST
    { iadst4_c, iadst4_c },  // ADST_FLIPADST
    { iadst4_c, iadst4_c },  // FLIPADST_ADST
    { iidtx4_c, iidtx4_c },  // IDTX
    { idct4_c,  iidtx4_c },  // V_DCT
    { iidtx4_c, idct4_c  },  // H_DCT
    { iadst4_c, iidtx4_c },  // V_ADST
    { iidtx4_c, iadst4_c },  // H_ADST
    { iadst4_c, iidtx4_c },  // V_FLIPADST
    { iidtx4_c, iadst4_c },  // H_FLIPADST
#endif  // CONFIG_EXT_TX
  };

  int i, j;
  tran_low_t tmp;
  tran_low_t out[4][4];
  tran_low_t *outp = &out[0][0];
  int outstride = 4;

  // inverse transform row vectors
  for (i = 0; i < 4; ++i) {
    IHT_4[tx_type].rows(input, out[i]);
    input  += 4;
  }

  // transpose
  for (i = 1 ; i < 4; i++) {
    for (j = 0; j < i; j++) {
            tmp = out[i][j];
      out[i][j] = out[j][i];
      out[j][i] = tmp;
    }
  }

  // inverse transform column vectors
  for (i = 0; i < 4; ++i) {
    IHT_4[tx_type].cols(out[i], out[i]);
  }

#if CONFIG_EXT_TX
  maybe_flip_strides(&dest, &stride, &outp, &outstride, tx_type, 4);
#endif

  // Sum with the destination
  for (i = 0; i < 4; ++i) {
    for (j = 0; j < 4; ++j) {
      int d = i * stride + j;
      int s = j * outstride + i;
      dest[d] = clip_pixel_add(dest[d], ROUND_POWER_OF_TWO(outp[s], 4));
    }
  }
}

void vp10_iht8x8_64_add_c(const tran_low_t *input, uint8_t *dest, int stride,
                         int tx_type) {
  static const transform_2d IHT_8[] = {
    { idct8_c,  idct8_c  },  // DCT_DCT
    { iadst8_c, idct8_c  },  // ADST_DCT
    { idct8_c,  iadst8_c },  // DCT_ADST
    { iadst8_c, iadst8_c },  // ADST_ADST
#if CONFIG_EXT_TX
    { iadst8_c, idct8_c  },  // FLIPADST_DCT
    { idct8_c,  iadst8_c },  // DCT_FLIPADST
    { iadst8_c, iadst8_c },  // FLIPADST_FLIPADST
    { iadst8_c, iadst8_c },  // ADST_FLIPADST
    { iadst8_c, iadst8_c },  // FLIPADST_ADST
    { iidtx8_c, iidtx8_c },  // IDTX
    { idct8_c,  iidtx8_c },  // V_DCT
    { iidtx8_c, idct8_c  },  // H_DCT
    { iadst8_c, iidtx8_c },  // V_ADST
    { iidtx8_c, iadst8_c },  // H_ADST
    { iadst8_c, iidtx8_c },  // V_FLIPADST
    { iidtx8_c, iadst8_c },  // H_FLIPADST
#endif  // CONFIG_EXT_TX
  };

  int i, j;
  tran_low_t tmp;
  tran_low_t out[8][8];
  tran_low_t *outp = &out[0][0];
  int outstride = 8;

  // inverse transform row vectors
  for (i = 0; i < 8; ++i) {
    IHT_8[tx_type].rows(input, out[i]);
    input  += 8;
  }

  // transpose
  for (i = 1 ; i < 8; i++) {
    for (j = 0; j < i; j++) {
            tmp = out[i][j];
      out[i][j] = out[j][i];
      out[j][i] = tmp;
    }
  }

  // inverse transform column vectors
  for (i = 0; i < 8; ++i) {
    IHT_8[tx_type].cols(out[i], out[i]);
  }

#if CONFIG_EXT_TX
  maybe_flip_strides(&dest, &stride, &outp, &outstride, tx_type, 8);
#endif

  // Sum with the destination
  for (i = 0; i < 8; ++i) {
    for (j = 0; j < 8; ++j) {
      int d = i * stride + j;
      int s = j * outstride + i;
      dest[d] = clip_pixel_add(dest[d], ROUND_POWER_OF_TWO(outp[s], 5));
    }
  }
}

void vp10_iht16x16_256_add_c(const tran_low_t *input, uint8_t *dest, int stride,
                            int tx_type) {
  static const transform_2d IHT_16[] = {
    { idct16_c,  idct16_c  },  // DCT_DCT
    { iadst16_c, idct16_c  },  // ADST_DCT
    { idct16_c,  iadst16_c },  // DCT_ADST
    { iadst16_c, iadst16_c },  // ADST_ADST
#if CONFIG_EXT_TX
    { iadst16_c, idct16_c  },  // FLIPADST_DCT
    { idct16_c,  iadst16_c },  // DCT_FLIPADST
    { iadst16_c, iadst16_c },  // FLIPADST_FLIPADST
    { iadst16_c, iadst16_c },  // ADST_FLIPADST
    { iadst16_c, iadst16_c },  // FLIPADST_ADST
    { iidtx16_c, iidtx16_c },  // IDTX
    { idct16_c,  iidtx16_c },  // V_DCT
    { iidtx16_c, idct16_c  },  // H_DCT
    { iadst16_c, iidtx16_c },  // V_ADST
    { iidtx16_c, iadst16_c },  // H_ADST
    { iadst16_c, iidtx16_c },  // V_FLIPADST
    { iidtx16_c, iadst16_c },  // H_FLIPADST
#endif  // CONFIG_EXT_TX
  };

  int i, j;
  tran_low_t tmp;
  tran_low_t out[16][16];
  tran_low_t *outp = &out[0][0];
  int outstride = 16;

  // inverse transform row vectors
  for (i = 0; i < 16; ++i) {
    IHT_16[tx_type].rows(input, out[i]);
    input  += 16;
  }

  // transpose
  for (i = 1 ; i < 16; i++) {
    for (j = 0; j < i; j++) {
            tmp = out[i][j];
      out[i][j] = out[j][i];
      out[j][i] = tmp;
    }
  }

  // inverse transform column vectors
  for (i = 0; i < 16; ++i) {
    IHT_16[tx_type].cols(out[i], out[i]);
  }

#if CONFIG_EXT_TX
  maybe_flip_strides(&dest, &stride, &outp, &outstride, tx_type, 16);
#endif

  // Sum with the destination
  for (i = 0; i < 16; ++i) {
    for (j = 0; j < 16; ++j) {
      int d = i * stride + j;
      int s = j * outstride + i;
      dest[d] = clip_pixel_add(dest[d], ROUND_POWER_OF_TWO(outp[s], 6));
    }
  }
}

#if CONFIG_EXT_TX
void vp10_iht32x32_1024_add_c(const tran_low_t *input, uint8_t *dest,
                              int stride, int tx_type) {
  static const transform_2d IHT_32[] = {
    { idct32_c,  idct32_c  },                // DCT_DCT
    { ihalfright32_c, idct32_c  },           // ADST_DCT
    { idct32_c,  ihalfright32_c },           // DCT_ADST
    { ihalfright32_c, ihalfright32_c },      // ADST_ADST
    { ihalfright32_c, idct32_c  },           // FLIPADST_DCT
    { idct32_c,  ihalfright32_c },           // DCT_FLIPADST
    { ihalfright32_c, ihalfright32_c },      // FLIPADST_FLIPADST
    { ihalfright32_c, ihalfright32_c },      // ADST_FLIPADST
    { ihalfright32_c, ihalfright32_c },      // FLIPADST_ADST
    { iidtx32_c, iidtx32_c },                // IDTX
    { idct32_c,  iidtx32_c },                // V_DCT
    { iidtx32_c, idct32_c  },                // H_DCT
    { ihalfright32_c, iidtx16_c },           // V_ADST
    { iidtx16_c, ihalfright32_c },           // H_ADST
    { ihalfright32_c, iidtx16_c },           // V_FLIPADST
    { iidtx16_c, ihalfright32_c },           // H_FLIPADST
  };

  int i, j;
  tran_low_t tmp;
  tran_low_t out[32][32];
  tran_low_t *outp = &out[0][0];
  int outstride = 32;

  // inverse transform row vectors
  for (i = 0; i < 32; ++i) {
    IHT_32[tx_type].rows(input, out[i]);
    input  += 32;
  }

  // transpose
  for (i = 1 ; i < 32; i++) {
    for (j = 0; j < i; j++) {
            tmp = out[i][j];
      out[i][j] = out[j][i];
      out[j][i] = tmp;
    }
  }

  // inverse transform column vectors
  for (i = 0; i < 32; ++i) {
    IHT_32[tx_type].cols(out[i], out[i]);
  }

  maybe_flip_strides(&dest, &stride, &outp, &outstride, tx_type, 32);

  // Sum with the destination
  for (i = 0; i < 32; ++i) {
    for (j = 0; j < 32; ++j) {
      int d = i * stride + j;
      int s = j * outstride + i;
      dest[d] = clip_pixel_add(dest[d], ROUND_POWER_OF_TWO(outp[s], 6));
    }
  }
}
#endif  // CONFIG_EXT_TX

// idct
void vp10_idct4x4_add(const tran_low_t *input, uint8_t *dest, int stride,
                     int eob) {
  if (eob > 1)
    vpx_idct4x4_16_add(input, dest, stride);
  else
    vpx_idct4x4_1_add(input, dest, stride);
}


void vp10_iwht4x4_add(const tran_low_t *input, uint8_t *dest, int stride,
                     int eob) {
  if (eob > 1)
    vpx_iwht4x4_16_add(input, dest, stride);
  else
    vpx_iwht4x4_1_add(input, dest, stride);
}

void vp10_idct8x8_add(const tran_low_t *input, uint8_t *dest, int stride,
                     int eob) {
  // If dc is 1, then input[0] is the reconstructed value, do not need
  // dequantization. Also, when dc is 1, dc is counted in eobs, namely eobs >=1.

  // The calculation can be simplified if there are not many non-zero dct
  // coefficients. Use eobs to decide what to do.
  // TODO(yunqingwang): "eobs = 1" case is also handled in vp10_short_idct8x8_c.
  // Combine that with code here.
  if (eob == 1)
    // DC only DCT coefficient
    vpx_idct8x8_1_add(input, dest, stride);
  else if (eob <= 12)
    vpx_idct8x8_12_add(input, dest, stride);
  else
    vpx_idct8x8_64_add(input, dest, stride);
}

void vp10_idct16x16_add(const tran_low_t *input, uint8_t *dest, int stride,
                       int eob) {
  /* The calculation can be simplified if there are not many non-zero dct
   * coefficients. Use eobs to separate different cases. */
  if (eob == 1)
    /* DC only DCT coefficient. */
    vpx_idct16x16_1_add(input, dest, stride);
  else if (eob <= 10)
    vpx_idct16x16_10_add(input, dest, stride);
  else
    vpx_idct16x16_256_add(input, dest, stride);
}

void vp10_idct32x32_add(const tran_low_t *input, uint8_t *dest, int stride,
                       int eob) {
  if (eob == 1)
    vpx_idct32x32_1_add(input, dest, stride);
  else if (eob <= 34)
    // non-zero coeff only in upper-left 8x8
    vpx_idct32x32_34_add(input, dest, stride);
  else
    vpx_idct32x32_1024_add(input, dest, stride);
}

void vp10_inv_txfm_add_4x4(const tran_low_t *input, uint8_t *dest,
                           int stride, int eob, TX_TYPE tx_type, int lossless) {
  if (lossless) {
    assert(tx_type == DCT_DCT);
    vp10_iwht4x4_add(input, dest, stride, eob);
    return;
  }

  switch (tx_type) {
    case DCT_DCT:
      vp10_idct4x4_add(input, dest, stride, eob);
      break;
    case ADST_DCT:
    case DCT_ADST:
    case ADST_ADST:
      vp10_iht4x4_16_add(input, dest, stride, tx_type);
      break;
#if CONFIG_EXT_TX
    case FLIPADST_DCT:
    case DCT_FLIPADST:
    case FLIPADST_FLIPADST:
    case ADST_FLIPADST:
    case FLIPADST_ADST:
      vp10_iht4x4_16_add(input, dest, stride, tx_type);
      break;
    case V_DCT:
    case H_DCT:
    case V_ADST:
    case H_ADST:
    case V_FLIPADST:
    case H_FLIPADST:
      // Use C version since DST only exists in C code
      vp10_iht4x4_16_add_c(input, dest, stride, tx_type);
      break;
    case IDTX:
      inv_idtx_add_c(input, dest, stride, 4, tx_type);
      break;
#endif  // CONFIG_EXT_TX
    default:
      assert(0);
      break;
  }
}

void vp10_inv_txfm_add_8x8(const tran_low_t *input, uint8_t *dest,
                           int stride, int eob, TX_TYPE tx_type) {
  switch (tx_type) {
    case DCT_DCT:
      vp10_idct8x8_add(input, dest, stride, eob);
      break;
    case ADST_DCT:
    case DCT_ADST:
    case ADST_ADST:
      vp10_iht8x8_64_add(input, dest, stride, tx_type);
      break;
#if CONFIG_EXT_TX
    case FLIPADST_DCT:
    case DCT_FLIPADST:
    case FLIPADST_FLIPADST:
    case ADST_FLIPADST:
    case FLIPADST_ADST:
      vp10_iht8x8_64_add(input, dest, stride, tx_type);
      break;
    case V_DCT:
    case H_DCT:
    case V_ADST:
    case H_ADST:
    case V_FLIPADST:
    case H_FLIPADST:
      // Use C version since DST only exists in C code
      vp10_iht8x8_64_add_c(input, dest, stride, tx_type);
      break;
    case IDTX:
      inv_idtx_add_c(input, dest, stride, 8, tx_type);
      break;
#endif  // CONFIG_EXT_TX
    default:
      assert(0);
      break;
  }
}

void vp10_inv_txfm_add_16x16(const tran_low_t *input, uint8_t *dest,
                             int stride, int eob, TX_TYPE tx_type) {
  switch (tx_type) {
    case DCT_DCT:
      vp10_idct16x16_add(input, dest, stride, eob);
      break;
    case ADST_DCT:
    case DCT_ADST:
    case ADST_ADST:
      vp10_iht16x16_256_add(input, dest, stride, tx_type);
      break;
#if CONFIG_EXT_TX
    case FLIPADST_DCT:
    case DCT_FLIPADST:
    case FLIPADST_FLIPADST:
    case ADST_FLIPADST:
    case FLIPADST_ADST:
      vp10_iht16x16_256_add(input, dest, stride, tx_type);
      break;
    case V_DCT:
    case H_DCT:
    case V_ADST:
    case H_ADST:
    case V_FLIPADST:
    case H_FLIPADST:
      // Use C version since DST only exists in C code
      vp10_iht16x16_256_add_c(input, dest, stride, tx_type);
      break;
    case IDTX:
      inv_idtx_add_c(input, dest, stride, 16, tx_type);
      break;
#endif  // CONFIG_EXT_TX
    default:
      assert(0);
      break;
  }
}

void vp10_inv_txfm_add_32x32(const tran_low_t *input, uint8_t *dest,
                             int stride, int eob, TX_TYPE tx_type) {
  switch (tx_type) {
    case DCT_DCT:
      vp10_idct32x32_add(input, dest, stride, eob);
      break;
#if CONFIG_EXT_TX
    case ADST_DCT:
    case DCT_ADST:
    case ADST_ADST:
    case FLIPADST_DCT:
    case DCT_FLIPADST:
    case FLIPADST_FLIPADST:
    case ADST_FLIPADST:
    case FLIPADST_ADST:
    case V_DCT:
    case H_DCT:
    case V_ADST:
    case H_ADST:
    case V_FLIPADST:
    case H_FLIPADST:
      vp10_iht32x32_1024_add_c(input, dest, stride, tx_type);
      break;
    case IDTX:
      inv_idtx_add_c(input, dest, stride, 32, tx_type);
      break;
#endif  // CONFIG_EXT_TX
    default:
      assert(0);
      break;
  }
}

#if CONFIG_VP9_HIGHBITDEPTH
void vp10_highbd_iht4x4_16_add_c(const tran_low_t *input, uint8_t *dest8,
                                int stride, int tx_type, int bd) {
  static const highbd_transform_2d HIGH_IHT_4[] = {
    { vpx_highbd_idct4_c,  vpx_highbd_idct4_c  },  // DCT_DCT
    { vpx_highbd_iadst4_c, vpx_highbd_idct4_c  },  // ADST_DCT
    { vpx_highbd_idct4_c,  vpx_highbd_iadst4_c },  // DCT_ADST
    { vpx_highbd_iadst4_c, vpx_highbd_iadst4_c },  // ADST_ADST
#if CONFIG_EXT_TX
    { vpx_highbd_iadst4_c, vpx_highbd_idct4_c  },  // FLIPADST_DCT
    { vpx_highbd_idct4_c,  vpx_highbd_iadst4_c },  // DCT_FLIPADST
    { vpx_highbd_iadst4_c, vpx_highbd_iadst4_c },  // FLIPADST_FLIPADST
    { vpx_highbd_iadst4_c, vpx_highbd_iadst4_c },  // ADST_FLIPADST
    { vpx_highbd_iadst4_c, vpx_highbd_iadst4_c },  // FLIPADST_ADST
    {     highbd_iidtx4_c,     highbd_iidtx4_c },  // IDTX
    { vpx_highbd_idct4_c,      highbd_iidtx4_c },  // V_DCT
    {     highbd_iidtx4_c, vpx_highbd_idct4_c  },  // H_DCT
    { vpx_highbd_iadst4_c,     highbd_iidtx4_c },  // V_ADST
    {     highbd_iidtx4_c, vpx_highbd_iadst4_c },  // H_ADST
    { vpx_highbd_iadst4_c,     highbd_iidtx4_c },  // V_FLIPADST
    {     highbd_iidtx4_c, vpx_highbd_iadst4_c },  // H_FLIPADST
#endif  // CONFIG_EXT_TX
  };

  uint16_t *dest = CONVERT_TO_SHORTPTR(dest8);

  int i, j;
  tran_low_t tmp;
  tran_low_t out[4][4];
  tran_low_t *outp = &out[0][0];
  int outstride = 4;

  // inverse transform row vectors
  for (i = 0; i < 4; ++i) {
    HIGH_IHT_4[tx_type].rows(input, out[i], bd);
    input  += 4;
  }

  // transpose
  for (i = 1 ; i < 4; i++) {
    for (j = 0; j < i; j++) {
            tmp = out[i][j];
      out[i][j] = out[j][i];
      out[j][i] = tmp;
    }
  }

  // inverse transform column vectors
  for (i = 0; i < 4; ++i) {
    HIGH_IHT_4[tx_type].cols(out[i], out[i], bd);
  }

#if CONFIG_EXT_TX
  maybe_flip_strides16(&dest, &stride, &outp, &outstride, tx_type, 4);
#endif

  // Sum with the destination
  for (i = 0; i < 4; ++i) {
    for (j = 0; j < 4; ++j) {
      int d = i * stride + j;
      int s = j * outstride + i;
      dest[d] = highbd_clip_pixel_add(dest[d],
                                      ROUND_POWER_OF_TWO(outp[s], 4), bd);
    }
  }
}

void vp10_highbd_iht8x8_64_add_c(const tran_low_t *input, uint8_t *dest8,
                                int stride, int tx_type, int bd) {
  static const highbd_transform_2d HIGH_IHT_8[] = {
    { vpx_highbd_idct8_c,  vpx_highbd_idct8_c  },  // DCT_DCT
    { vpx_highbd_iadst8_c, vpx_highbd_idct8_c  },  // ADST_DCT
    { vpx_highbd_idct8_c,  vpx_highbd_iadst8_c },  // DCT_ADST
    { vpx_highbd_iadst8_c, vpx_highbd_iadst8_c },  // ADST_ADST
#if CONFIG_EXT_TX
    { vpx_highbd_iadst8_c, vpx_highbd_idct8_c  },  // FLIPADST_DCT
    { vpx_highbd_idct8_c,  vpx_highbd_iadst8_c },  // DCT_FLIPADST
    { vpx_highbd_iadst8_c, vpx_highbd_iadst8_c },  // FLIPADST_FLIPADST
    { vpx_highbd_iadst8_c, vpx_highbd_iadst8_c },  // ADST_FLIPADST
    { vpx_highbd_iadst8_c, vpx_highbd_iadst8_c },  // FLIPADST_ADST
    {     highbd_iidtx8_c,     highbd_iidtx8_c },  // IDTX
    { vpx_highbd_idct8_c,      highbd_iidtx8_c },  // V_DCT
    {     highbd_iidtx8_c, vpx_highbd_idct8_c  },  // H_DCT
    { vpx_highbd_iadst8_c,     highbd_iidtx8_c },  // V_ADST
    {     highbd_iidtx8_c, vpx_highbd_iadst8_c },  // H_ADST
    { vpx_highbd_iadst8_c,     highbd_iidtx8_c },  // V_FLIPADST
    {     highbd_iidtx8_c, vpx_highbd_iadst8_c },  // H_FLIPADST
#endif  // CONFIG_EXT_TX
  };

  uint16_t *dest = CONVERT_TO_SHORTPTR(dest8);

  int i, j;
  tran_low_t tmp;
  tran_low_t out[8][8];
  tran_low_t *outp = &out[0][0];
  int outstride = 8;

  // inverse transform row vectors
  for (i = 0; i < 8; ++i) {
    HIGH_IHT_8[tx_type].rows(input, out[i], bd);
    input  += 8;
  }

  // transpose
  for (i = 1 ; i < 8; i++) {
    for (j = 0; j < i; j++) {
            tmp = out[i][j];
      out[i][j] = out[j][i];
      out[j][i] = tmp;
    }
  }

  // inverse transform column vectors
  for (i = 0; i < 8; ++i) {
    HIGH_IHT_8[tx_type].cols(out[i], out[i], bd);
  }

#if CONFIG_EXT_TX
  maybe_flip_strides16(&dest, &stride, &outp, &outstride, tx_type, 8);
#endif

  // Sum with the destination
  for (i = 0; i < 8; ++i) {
    for (j = 0; j < 8; ++j) {
      int d = i * stride + j;
      int s = j * outstride + i;
      dest[d] = highbd_clip_pixel_add(dest[d],
                                      ROUND_POWER_OF_TWO(outp[s], 5), bd);
    }
  }
}

void vp10_highbd_iht16x16_256_add_c(const tran_low_t *input, uint8_t *dest8,
                                   int stride, int tx_type, int bd) {
  static const highbd_transform_2d HIGH_IHT_16[] = {
    { vpx_highbd_idct16_c,  vpx_highbd_idct16_c  },  // DCT_DCT
    { vpx_highbd_iadst16_c, vpx_highbd_idct16_c  },  // ADST_DCT
    { vpx_highbd_idct16_c,  vpx_highbd_iadst16_c },  // DCT_ADST
    { vpx_highbd_iadst16_c, vpx_highbd_iadst16_c },  // ADST_ADST
#if CONFIG_EXT_TX
    { vpx_highbd_iadst16_c, vpx_highbd_idct16_c  },  // FLIPADST_DCT
    { vpx_highbd_idct16_c,  vpx_highbd_iadst16_c },  // DCT_FLIPADST
    { vpx_highbd_iadst16_c, vpx_highbd_iadst16_c },  // FLIPADST_FLIPADST
    { vpx_highbd_iadst16_c, vpx_highbd_iadst16_c },  // ADST_FLIPADST
    { vpx_highbd_iadst16_c, vpx_highbd_iadst16_c },  // FLIPADST_ADST
    {     highbd_iidtx16_c,     highbd_iidtx16_c },  // IDTX
    { vpx_highbd_idct16_c,      highbd_iidtx16_c },  // V_DCT
    {     highbd_iidtx16_c, vpx_highbd_idct16_c  },  // H_DCT
    { vpx_highbd_iadst16_c,     highbd_iidtx16_c },  // V_ADST
    {     highbd_iidtx16_c, vpx_highbd_iadst16_c },  // H_ADST
    { vpx_highbd_iadst16_c,     highbd_iidtx16_c },  // V_FLIPADST
    {     highbd_iidtx16_c, vpx_highbd_iadst16_c },  // H_FLIPADST
#endif  // CONFIG_EXT_TX
  };

  uint16_t *dest = CONVERT_TO_SHORTPTR(dest8);

  int i, j;
  tran_low_t tmp;
  tran_low_t out[16][16];
  tran_low_t *outp = &out[0][0];
  int outstride = 16;

  // inverse transform row vectors
  for (i = 0; i < 16; ++i) {
    HIGH_IHT_16[tx_type].rows(input, out[i], bd);
    input  += 16;
  }

  // transpose
  for (i = 1 ; i < 16; i++) {
    for (j = 0; j < i; j++) {
            tmp = out[i][j];
      out[i][j] = out[j][i];
      out[j][i] = tmp;
    }
  }

  // inverse transform column vectors
  for (i = 0; i < 16; ++i) {
    HIGH_IHT_16[tx_type].cols(out[i], out[i], bd);
  }

#if CONFIG_EXT_TX
  maybe_flip_strides16(&dest, &stride, &outp, &outstride, tx_type, 16);
#endif

  // Sum with the destination
  for (i = 0; i < 16; ++i) {
    for (j = 0; j < 16; ++j) {
      int d = i * stride + j;
      int s = j * outstride + i;
      dest[d] = highbd_clip_pixel_add(dest[d],
                                      ROUND_POWER_OF_TWO(outp[s], 6), bd);
    }
  }
}

#if CONFIG_EXT_TX
void vp10_highbd_iht32x32_1024_add_c(const tran_low_t *input, uint8_t *dest8,
                                     int stride, int tx_type, int bd) {
  static const highbd_transform_2d HIGH_IHT_32[] = {
    { vpx_highbd_idct32_c,    vpx_highbd_idct32_c    },  // DCT_DCT
    { highbd_ihalfright32_c,  vpx_highbd_idct32_c    },  // ADST_DCT
    { vpx_highbd_idct32_c,    highbd_ihalfright32_c  },  // DCT_ADST
    { highbd_ihalfright32_c,  highbd_ihalfright32_c  },  // ADST_ADST
    { highbd_ihalfright32_c,  vpx_highbd_idct32_c    },  // FLIPADST_DCT
    { vpx_highbd_idct32_c,    highbd_ihalfright32_c  },  // DCT_FLIPADST
    { highbd_ihalfright32_c,  highbd_ihalfright32_c  },  // FLIPADST_FLIPADST
    { highbd_ihalfright32_c,  highbd_ihalfright32_c  },  // ADST_FLIPADST
    { highbd_ihalfright32_c,  highbd_ihalfright32_c  },  // FLIPADST_ADST
    { highbd_iidtx32_c,       highbd_iidtx32_c       },  // IDTX
    { vpx_highbd_idct32_c,    highbd_iidtx32_c       },  // V_DCT
    { highbd_iidtx32_c,       vpx_highbd_idct32_c    },  // H_DCT
    { highbd_ihalfright32_c,  highbd_iidtx32_c       },  // V_ADST
    { highbd_iidtx32_c,       highbd_ihalfright32_c  },  // H_ADST
    { highbd_ihalfright32_c,  highbd_iidtx32_c       },  // V_FLIPADST
    { highbd_iidtx32_c,       highbd_ihalfright32_c  },  // H_FLIPADST
  };

  uint16_t *dest = CONVERT_TO_SHORTPTR(dest8);

  int i, j;
  tran_low_t tmp;
  tran_low_t out[32][32];
  tran_low_t *outp = &out[0][0];
  int outstride = 32;

  // inverse transform row vectors
  for (i = 0; i < 32; ++i) {
    HIGH_IHT_32[tx_type].rows(input, out[i], bd);
    input  += 32;
  }

  // transpose
  for (i = 1 ; i < 32; i++) {
    for (j = 0; j < i; j++) {
            tmp = out[i][j];
      out[i][j] = out[j][i];
      out[j][i] = tmp;
    }
  }

  // inverse transform column vectors
  for (i = 0; i < 32; ++i) {
    HIGH_IHT_32[tx_type].cols(out[i], out[i], bd);
  }

  maybe_flip_strides16(&dest, &stride, &outp, &outstride, tx_type, 32);

  // Sum with the destination
  for (i = 0; i < 32; ++i) {
    for (j = 0; j < 32; ++j) {
      int d = i * stride + j;
      int s = j * outstride + i;
      dest[d] = highbd_clip_pixel_add(dest[d],
                                      ROUND_POWER_OF_TWO(outp[s], 6), bd);
    }
  }
}
#endif  // CONFIG_EXT_TX

// idct
void vp10_highbd_idct4x4_add(const tran_low_t *input, uint8_t *dest, int stride,
                            int eob, int bd) {
  if (eob > 1)
    vpx_highbd_idct4x4_16_add(input, dest, stride, bd);
  else
    vpx_highbd_idct4x4_1_add(input, dest, stride, bd);
}


void vp10_highbd_iwht4x4_add(const tran_low_t *input, uint8_t *dest, int stride,
                            int eob, int bd) {
  if (eob > 1)
    vpx_highbd_iwht4x4_16_add(input, dest, stride, bd);
  else
    vpx_highbd_iwht4x4_1_add(input, dest, stride, bd);
}

void vp10_highbd_idct8x8_add(const tran_low_t *input, uint8_t *dest, int stride,
                            int eob, int bd) {
  // If dc is 1, then input[0] is the reconstructed value, do not need
  // dequantization. Also, when dc is 1, dc is counted in eobs, namely eobs >=1.

  // The calculation can be simplified if there are not many non-zero dct
  // coefficients. Use eobs to decide what to do.
  // TODO(yunqingwang): "eobs = 1" case is also handled in vp10_short_idct8x8_c.
  // Combine that with code here.
  // DC only DCT coefficient
  if (eob == 1) {
    vpx_highbd_idct8x8_1_add(input, dest, stride, bd);
  } else if (eob <= 10) {
    vpx_highbd_idct8x8_10_add(input, dest, stride, bd);
  } else {
    vpx_highbd_idct8x8_64_add(input, dest, stride, bd);
  }
}

void vp10_highbd_idct16x16_add(const tran_low_t *input, uint8_t *dest,
                              int stride, int eob, int bd) {
  // The calculation can be simplified if there are not many non-zero dct
  // coefficients. Use eobs to separate different cases.
  // DC only DCT coefficient.
  if (eob == 1) {
    vpx_highbd_idct16x16_1_add(input, dest, stride, bd);
  } else if (eob <= 10) {
    vpx_highbd_idct16x16_10_add(input, dest, stride, bd);
  } else {
    vpx_highbd_idct16x16_256_add(input, dest, stride, bd);
  }
}

void vp10_highbd_idct32x32_add(const tran_low_t *input, uint8_t *dest,
                              int stride, int eob, int bd) {
  // Non-zero coeff only in upper-left 8x8
  if (eob == 1) {
    vpx_highbd_idct32x32_1_add(input, dest, stride, bd);
  } else if (eob <= 34) {
    vpx_highbd_idct32x32_34_add(input, dest, stride, bd);
  } else {
    vpx_highbd_idct32x32_1024_add(input, dest, stride, bd);
  }
}

void vp10_highbd_inv_txfm_add_4x4(const tran_low_t *input, uint8_t *dest,
                                  int stride, int eob, int bd, TX_TYPE tx_type,
                                  int lossless) {
  if (lossless) {
    assert(tx_type == DCT_DCT);
    vp10_highbd_iwht4x4_add(input, dest, stride, eob, bd);
    return;
  }

  switch (tx_type) {
    case DCT_DCT:
    case ADST_DCT:
    case DCT_ADST:
    case ADST_ADST:
      vp10_inv_txfm2d_add_4x4(input, CONVERT_TO_SHORTPTR(dest), stride,
                              tx_type, bd);
      break;
#if CONFIG_EXT_TX
    case FLIPADST_DCT:
    case DCT_FLIPADST:
    case FLIPADST_FLIPADST:
    case ADST_FLIPADST:
    case FLIPADST_ADST:
      vp10_inv_txfm2d_add_4x4(input, CONVERT_TO_SHORTPTR(dest), stride,
                              tx_type, bd);
      break;
    case V_DCT:
    case H_DCT:
    case V_ADST:
    case H_ADST:
    case V_FLIPADST:
    case H_FLIPADST:
      // Use C version since DST only exists in C code
      vp10_highbd_iht4x4_16_add_c(input, dest, stride, tx_type, bd);
      break;
    case IDTX:
      highbd_inv_idtx_add_c(input, dest, stride, 4, tx_type, bd);
      break;
#endif  // CONFIG_EXT_TX
    default:
      assert(0);
      break;
  }
}

void vp10_highbd_inv_txfm_add_8x8(const tran_low_t *input, uint8_t *dest,
                                  int stride, int eob, int bd,
                                  TX_TYPE tx_type) {
  (void)eob;
  switch (tx_type) {
    case DCT_DCT:
    case ADST_DCT:
    case DCT_ADST:
    case ADST_ADST:
      vp10_inv_txfm2d_add_8x8(input, CONVERT_TO_SHORTPTR(dest), stride,
                              tx_type, bd);
      break;
#if CONFIG_EXT_TX
    case FLIPADST_DCT:
    case DCT_FLIPADST:
    case FLIPADST_FLIPADST:
    case ADST_FLIPADST:
    case FLIPADST_ADST:
      vp10_inv_txfm2d_add_8x8(input, CONVERT_TO_SHORTPTR(dest), stride,
                              tx_type, bd);
      break;
    case V_DCT:
    case H_DCT:
    case V_ADST:
    case H_ADST:
    case V_FLIPADST:
    case H_FLIPADST:
      // Use C version since DST only exists in C code
      vp10_highbd_iht8x8_64_add_c(input, dest, stride, tx_type, bd);
      break;
    case IDTX:
      highbd_inv_idtx_add_c(input, dest, stride, 8, tx_type, bd);
      break;
#endif  // CONFIG_EXT_TX
    default:
      assert(0);
      break;
  }
}

void vp10_highbd_inv_txfm_add_16x16(const tran_low_t *input, uint8_t *dest,
                                    int stride, int eob, int bd,
                                    TX_TYPE tx_type) {
  (void)eob;
  switch (tx_type) {
    case DCT_DCT:
    case ADST_DCT:
    case DCT_ADST:
    case ADST_ADST:
      vp10_inv_txfm2d_add_16x16(input, CONVERT_TO_SHORTPTR(dest), stride,
                                tx_type, bd);
      break;
#if CONFIG_EXT_TX
    case FLIPADST_DCT:
    case DCT_FLIPADST:
    case FLIPADST_FLIPADST:
    case ADST_FLIPADST:
    case FLIPADST_ADST:
      vp10_inv_txfm2d_add_16x16(input, CONVERT_TO_SHORTPTR(dest), stride,
                                tx_type, bd);
      break;
    case V_DCT:
    case H_DCT:
    case V_ADST:
    case H_ADST:
    case V_FLIPADST:
    case H_FLIPADST:
      // Use C version since DST only exists in C code
      vp10_highbd_iht16x16_256_add_c(input, dest, stride, tx_type, bd);
      break;
    case IDTX:
      highbd_inv_idtx_add_c(input, dest, stride, 16, tx_type, bd);
      break;
#endif  // CONFIG_EXT_TX
    default:
      assert(0);
      break;
  }
}

void vp10_highbd_inv_txfm_add_32x32(const tran_low_t *input, uint8_t *dest,
                                    int stride, int eob, int bd,
                                    TX_TYPE tx_type) {
  (void)eob;
  switch (tx_type) {
    case DCT_DCT:
      vp10_inv_txfm2d_add_32x32(input, CONVERT_TO_SHORTPTR(dest), stride,
                                DCT_DCT, bd);
      break;
#if CONFIG_EXT_TX
    case ADST_DCT:
    case DCT_ADST:
    case ADST_ADST:
    case FLIPADST_DCT:
    case DCT_FLIPADST:
    case FLIPADST_FLIPADST:
    case ADST_FLIPADST:
    case FLIPADST_ADST:
    case V_DCT:
    case H_DCT:
    case V_ADST:
    case H_ADST:
    case V_FLIPADST:
    case H_FLIPADST:
      vp10_highbd_iht32x32_1024_add_c(input, dest, stride, tx_type, bd);
      break;
    case IDTX:
      highbd_inv_idtx_add_c(input, dest, stride, 32, tx_type, bd);
      break;
#endif  // CONFIG_EXT_TX
    default:
      assert(0);
      break;
  }
}
#endif  // CONFIG_VP9_HIGHBITDEPTH

void inv_txfm_add(const tran_low_t *input, uint8_t *dest, int stride,
                  INV_TXFM_PARAM *inv_txfm_param) {
  const TX_TYPE tx_type = inv_txfm_param->tx_type;
  const TX_SIZE tx_size = inv_txfm_param->tx_size;
  const int eob = inv_txfm_param->eob;
  const int lossless = inv_txfm_param->lossless;

  switch (tx_size) {
    case TX_32X32:
      vp10_inv_txfm_add_32x32(input, dest, stride, eob, tx_type);
      break;
    case TX_16X16:
      vp10_inv_txfm_add_16x16(input, dest, stride, eob, tx_type);
      break;
    case TX_8X8:
      vp10_inv_txfm_add_8x8(input, dest, stride, eob, tx_type);
      break;
    case TX_4X4:
      // this is like vp10_short_idct4x4 but has a special case around eob<=1
      // which is significant (not just an optimization) for the lossless
      // case.
      vp10_inv_txfm_add_4x4(input, dest, stride, eob, tx_type,
                            lossless);
      break;
    default:
      assert(0 && "Invalid transform size");
      break;
  }
}

#if CONFIG_VP9_HIGHBITDEPTH
void highbd_inv_txfm_add(const tran_low_t *input, uint8_t *dest, int stride,
                         INV_TXFM_PARAM *inv_txfm_param) {
  const TX_TYPE tx_type = inv_txfm_param->tx_type;
  const TX_SIZE tx_size = inv_txfm_param->tx_size;
  const int eob = inv_txfm_param->eob;
  const int bd = inv_txfm_param->bd;
  const int lossless = inv_txfm_param->lossless;

  switch (tx_size) {
    case TX_32X32:
      vp10_highbd_inv_txfm_add_32x32(input, dest, stride, eob, bd, tx_type);
      break;
    case TX_16X16:
      vp10_highbd_inv_txfm_add_16x16(input, dest, stride, eob, bd, tx_type);
      break;
    case TX_8X8:
      vp10_highbd_inv_txfm_add_8x8(input, dest, stride, eob, bd, tx_type);
      break;
    case TX_4X4:
      // this is like vp10_short_idct4x4 but has a special case around eob<=1
      // which is significant (not just an optimization) for the lossless
      // case.
      vp10_highbd_inv_txfm_add_4x4(input, dest, stride, eob, bd, tx_type,
                                   lossless);
      break;
    default:
      assert(0 && "Invalid transform size");
      break;
  }
}
#endif  // CONFIG_VP9_HIGHBITDEPTH
