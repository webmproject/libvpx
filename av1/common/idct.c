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

#include "./aom_dsp_rtcd.h"
#include "./av1_rtcd.h"
#include "aom_dsp/inv_txfm.h"
#include "aom_ports/mem.h"
#include "av1/common/av1_inv_txfm2d_cfg.h"
#include "av1/common/blockd.h"
#include "av1/common/enums.h"
#include "av1/common/idct.h"

int get_tx_scale(const MACROBLOCKD *const xd, const TX_TYPE tx_type,
                 const TX_SIZE tx_size) {
  (void)tx_type;
#if CONFIG_AOM_HIGHBITDEPTH
  if (xd->cur_buf->flags & YV12_FLAG_HIGHBITDEPTH) {
    return txsize_sqr_up_map[tx_size] == TX_32X32;
  }
#else
  (void)xd;
#endif
  return txsize_sqr_up_map[tx_size] == TX_32X32;
}

#if CONFIG_EXT_TX
static void iidtx4_c(const tran_low_t *input, tran_low_t *output) {
  int i;
  for (i = 0; i < 4; ++i)
    output[i] = (tran_low_t)dct_const_round_shift(input[i] * Sqrt2);
}

static void iidtx8_c(const tran_low_t *input, tran_low_t *output) {
  int i;
  for (i = 0; i < 8; ++i) output[i] = input[i] * 2;
}

static void iidtx16_c(const tran_low_t *input, tran_low_t *output) {
  int i;
  for (i = 0; i < 16; ++i)
    output[i] = (tran_low_t)dct_const_round_shift(input[i] * 2 * Sqrt2);
}

static void iidtx32_c(const tran_low_t *input, tran_low_t *output) {
  int i;
  for (i = 0; i < 32; ++i) output[i] = input[i] * 4;
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

#if CONFIG_AOM_HIGHBITDEPTH
static void highbd_iidtx4_c(const tran_low_t *input, tran_low_t *output,
                            int bd) {
  int i;
  for (i = 0; i < 4; ++i)
    output[i] =
        HIGHBD_WRAPLOW(highbd_dct_const_round_shift(input[i] * Sqrt2), bd);
}

static void highbd_iidtx8_c(const tran_low_t *input, tran_low_t *output,
                            int bd) {
  int i;
  (void)bd;
  for (i = 0; i < 8; ++i) output[i] = input[i] * 2;
}

static void highbd_iidtx16_c(const tran_low_t *input, tran_low_t *output,
                             int bd) {
  int i;
  for (i = 0; i < 16; ++i)
    output[i] =
        HIGHBD_WRAPLOW(highbd_dct_const_round_shift(input[i] * 2 * Sqrt2), bd);
}

static void highbd_iidtx32_c(const tran_low_t *input, tran_low_t *output,
                             int bd) {
  int i;
  (void)bd;
  for (i = 0; i < 32; ++i) output[i] = input[i] * 4;
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
    inputhalf[i] =
        HIGHBD_WRAPLOW(highbd_dct_const_round_shift(input[i] * Sqrt2), bd);
  }
  aom_highbd_idct16_c(inputhalf, output + 16, bd);
  // Note overall scaling factor is 4 times orthogonal
}
#endif  // CONFIG_AOM_HIGHBITDEPTH

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

#define FLIPUD_PTR(dest, stride, size)       \
  do {                                       \
    (dest) = (dest) + ((size)-1) * (stride); \
    (stride) = -(stride);                    \
  } while (0)

static void maybe_flip_strides(uint8_t **dst, int *dstride, tran_low_t **src,
                               int *sstride, int tx_type, int sizey,
                               int sizex) {
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
    case H_ADST: break;
    case FLIPADST_DCT:
    case FLIPADST_ADST:
    case V_FLIPADST:
      // flip UD
      FLIPUD_PTR(*dst, *dstride, sizey);
      break;
    case DCT_FLIPADST:
    case ADST_FLIPADST:
    case H_FLIPADST:
      // flip LR
      FLIPUD_PTR(*src, *sstride, sizex);
      break;
    case FLIPADST_FLIPADST:
      // flip UD
      FLIPUD_PTR(*dst, *dstride, sizey);
      // flip LR
      FLIPUD_PTR(*src, *sstride, sizex);
      break;
    default: assert(0); break;
  }
}

#if CONFIG_AOM_HIGHBITDEPTH
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

static void maybe_flip_strides16(uint16_t **dst, int *dstride, tran_low_t **src,
                                 int *sstride, int tx_type, int sizey,
                                 int sizex) {
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
    case H_ADST: break;
    case FLIPADST_DCT:
    case FLIPADST_ADST:
    case V_FLIPADST:
      // flip UD
      FLIPUD_PTR(*dst, *dstride, sizey);
      break;
    case DCT_FLIPADST:
    case ADST_FLIPADST:
    case H_FLIPADST:
      // flip LR
      FLIPUD_PTR(*src, *sstride, sizex);
      break;
    case FLIPADST_FLIPADST:
      // flip UD
      FLIPUD_PTR(*dst, *dstride, sizey);
      // flip LR
      FLIPUD_PTR(*src, *sstride, sizex);
      break;
    default: assert(0); break;
  }
}
#endif  // CONFIG_AOM_HIGHBITDEPTH
#endif  // CONFIG_EXT_TX

void av1_iht4x4_16_add_c(const tran_low_t *input, uint8_t *dest, int stride,
                         int tx_type) {
  static const transform_2d IHT_4[] = {
    { idct4_c, idct4_c },    // DCT_DCT
    { iadst4_c, idct4_c },   // ADST_DCT
    { idct4_c, iadst4_c },   // DCT_ADST
    { iadst4_c, iadst4_c },  // ADST_ADST
#if CONFIG_EXT_TX
    { iadst4_c, idct4_c },   // FLIPADST_DCT
    { idct4_c, iadst4_c },   // DCT_FLIPADST
    { iadst4_c, iadst4_c },  // FLIPADST_FLIPADST
    { iadst4_c, iadst4_c },  // ADST_FLIPADST
    { iadst4_c, iadst4_c },  // FLIPADST_ADST
    { iidtx4_c, iidtx4_c },  // IDTX
    { idct4_c, iidtx4_c },   // V_DCT
    { iidtx4_c, idct4_c },   // H_DCT
    { iadst4_c, iidtx4_c },  // V_ADST
    { iidtx4_c, iadst4_c },  // H_ADST
    { iadst4_c, iidtx4_c },  // V_FLIPADST
    { iidtx4_c, iadst4_c },  // H_FLIPADST
#endif                       // CONFIG_EXT_TX
  };

  int i, j;
  tran_low_t tmp;
  tran_low_t out[4][4];
  tran_low_t *outp = &out[0][0];
  int outstride = 4;

  // inverse transform row vectors
  for (i = 0; i < 4; ++i) {
    IHT_4[tx_type].rows(input, out[i]);
    input += 4;
  }

  // transpose
  for (i = 1; i < 4; i++) {
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
  maybe_flip_strides(&dest, &stride, &outp, &outstride, tx_type, 4, 4);
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

#if CONFIG_EXT_TX
void av1_iht4x8_32_add_c(const tran_low_t *input, uint8_t *dest, int stride,
                         int tx_type) {
  static const transform_2d IHT_4x8[] = {
    { idct8_c, idct4_c },    // DCT_DCT
    { iadst8_c, idct4_c },   // ADST_DCT
    { idct8_c, iadst4_c },   // DCT_ADST
    { iadst8_c, iadst4_c },  // ADST_ADST
    { iadst8_c, idct4_c },   // FLIPADST_DCT
    { idct8_c, iadst4_c },   // DCT_FLIPADST
    { iadst8_c, iadst4_c },  // FLIPADST_FLIPADST
    { iadst8_c, iadst4_c },  // ADST_FLIPADST
    { iadst8_c, iadst4_c },  // FLIPADST_ADST
    { iidtx8_c, iidtx4_c },  // IDTX
    { idct8_c, iidtx4_c },   // V_DCT
    { iidtx8_c, idct4_c },   // H_DCT
    { iadst8_c, iidtx4_c },  // V_ADST
    { iidtx8_c, iadst4_c },  // H_ADST
    { iadst8_c, iidtx4_c },  // V_FLIPADST
    { iidtx8_c, iadst4_c },  // H_FLIPADST
  };

  const int n = 4;
  const int n2 = 8;
  int i, j;
  tran_low_t out[4][8], outtmp[4];
  tran_low_t *outp = &out[0][0];
  int outstride = n2;

  // inverse transform row vectors and transpose
  for (i = 0; i < n2; ++i) {
    IHT_4x8[tx_type].rows(input, outtmp);
    for (j = 0; j < n; ++j)
      out[j][i] = (tran_low_t)dct_const_round_shift(outtmp[j] * Sqrt2);
    input += n;
  }

  // inverse transform column vectors
  for (i = 0; i < n; ++i) {
    IHT_4x8[tx_type].cols(out[i], out[i]);
  }

  maybe_flip_strides(&dest, &stride, &outp, &outstride, tx_type, n2, n);

  // Sum with the destination
  for (i = 0; i < n2; ++i) {
    for (j = 0; j < n; ++j) {
      int d = i * stride + j;
      int s = j * outstride + i;
      dest[d] = clip_pixel_add(dest[d], ROUND_POWER_OF_TWO(outp[s], 5));
    }
  }
}

void av1_iht8x4_32_add_c(const tran_low_t *input, uint8_t *dest, int stride,
                         int tx_type) {
  static const transform_2d IHT_8x4[] = {
    { idct4_c, idct8_c },    // DCT_DCT
    { iadst4_c, idct8_c },   // ADST_DCT
    { idct4_c, iadst8_c },   // DCT_ADST
    { iadst4_c, iadst8_c },  // ADST_ADST
    { iadst4_c, idct8_c },   // FLIPADST_DCT
    { idct4_c, iadst8_c },   // DCT_FLIPADST
    { iadst4_c, iadst8_c },  // FLIPADST_FLIPADST
    { iadst4_c, iadst8_c },  // ADST_FLIPADST
    { iadst4_c, iadst8_c },  // FLIPADST_ADST
    { iidtx4_c, iidtx8_c },  // IDTX
    { idct4_c, iidtx8_c },   // V_DCT
    { iidtx4_c, idct8_c },   // H_DCT
    { iadst4_c, iidtx8_c },  // V_ADST
    { iidtx4_c, iadst8_c },  // H_ADST
    { iadst4_c, iidtx8_c },  // V_FLIPADST
    { iidtx4_c, iadst8_c },  // H_FLIPADST
  };
  const int n = 4;
  const int n2 = 8;

  int i, j;
  tran_low_t out[8][4], outtmp[8];
  tran_low_t *outp = &out[0][0];
  int outstride = n;

  // inverse transform row vectors and transpose
  for (i = 0; i < n; ++i) {
    IHT_8x4[tx_type].rows(input, outtmp);
    for (j = 0; j < n2; ++j)
      out[j][i] = (tran_low_t)dct_const_round_shift(outtmp[j] * Sqrt2);
    input += n2;
  }

  // inverse transform column vectors
  for (i = 0; i < n2; ++i) {
    IHT_8x4[tx_type].cols(out[i], out[i]);
  }

  maybe_flip_strides(&dest, &stride, &outp, &outstride, tx_type, n, n2);

  // Sum with the destination
  for (i = 0; i < n; ++i) {
    for (j = 0; j < n2; ++j) {
      int d = i * stride + j;
      int s = j * outstride + i;
      dest[d] = clip_pixel_add(dest[d], ROUND_POWER_OF_TWO(outp[s], 5));
    }
  }
}

void av1_iht8x16_128_add_c(const tran_low_t *input, uint8_t *dest, int stride,
                           int tx_type) {
  static const transform_2d IHT_8x16[] = {
    { idct16_c, idct8_c },    // DCT_DCT
    { iadst16_c, idct8_c },   // ADST_DCT
    { idct16_c, iadst8_c },   // DCT_ADST
    { iadst16_c, iadst8_c },  // ADST_ADST
    { iadst16_c, idct8_c },   // FLIPADST_DCT
    { idct16_c, iadst8_c },   // DCT_FLIPADST
    { iadst16_c, iadst8_c },  // FLIPADST_FLIPADST
    { iadst16_c, iadst8_c },  // ADST_FLIPADST
    { iadst16_c, iadst8_c },  // FLIPADST_ADST
    { iidtx16_c, iidtx8_c },  // IDTX
    { idct16_c, iidtx8_c },   // V_DCT
    { iidtx16_c, idct8_c },   // H_DCT
    { iadst16_c, iidtx8_c },  // V_ADST
    { iidtx16_c, iadst8_c },  // H_ADST
    { iadst16_c, iidtx8_c },  // V_FLIPADST
    { iidtx16_c, iadst8_c },  // H_FLIPADST
  };

  const int n = 8;
  const int n2 = 16;
  int i, j;
  tran_low_t out[8][16], outtmp[8];
  tran_low_t *outp = &out[0][0];
  int outstride = n2;

  // inverse transform row vectors and transpose
  for (i = 0; i < n2; ++i) {
    IHT_8x16[tx_type].rows(input, outtmp);
    for (j = 0; j < n; ++j)
      out[j][i] = (tran_low_t)dct_const_round_shift(outtmp[j] * Sqrt2);
    input += n;
  }

  // inverse transform column vectors
  for (i = 0; i < n; ++i) {
    IHT_8x16[tx_type].cols(out[i], out[i]);
  }

  maybe_flip_strides(&dest, &stride, &outp, &outstride, tx_type, n2, n);

  // Sum with the destination
  for (i = 0; i < n2; ++i) {
    for (j = 0; j < n; ++j) {
      int d = i * stride + j;
      int s = j * outstride + i;
      dest[d] = clip_pixel_add(dest[d], ROUND_POWER_OF_TWO(outp[s], 6));
    }
  }
}

void av1_iht16x8_128_add_c(const tran_low_t *input, uint8_t *dest, int stride,
                           int tx_type) {
  static const transform_2d IHT_16x8[] = {
    { idct8_c, idct16_c },    // DCT_DCT
    { iadst8_c, idct16_c },   // ADST_DCT
    { idct8_c, iadst16_c },   // DCT_ADST
    { iadst8_c, iadst16_c },  // ADST_ADST
    { iadst8_c, idct16_c },   // FLIPADST_DCT
    { idct8_c, iadst16_c },   // DCT_FLIPADST
    { iadst8_c, iadst16_c },  // FLIPADST_FLIPADST
    { iadst8_c, iadst16_c },  // ADST_FLIPADST
    { iadst8_c, iadst16_c },  // FLIPADST_ADST
    { iidtx8_c, iidtx16_c },  // IDTX
    { idct8_c, iidtx16_c },   // V_DCT
    { iidtx8_c, idct16_c },   // H_DCT
    { iadst8_c, iidtx16_c },  // V_ADST
    { iidtx8_c, iadst16_c },  // H_ADST
    { iadst8_c, iidtx16_c },  // V_FLIPADST
    { iidtx8_c, iadst16_c },  // H_FLIPADST
  };
  const int n = 8;
  const int n2 = 16;

  int i, j;
  tran_low_t out[16][8], outtmp[16];
  tran_low_t *outp = &out[0][0];
  int outstride = n;

  // inverse transform row vectors and transpose
  for (i = 0; i < n; ++i) {
    IHT_16x8[tx_type].rows(input, outtmp);
    for (j = 0; j < n2; ++j)
      out[j][i] = (tran_low_t)dct_const_round_shift(outtmp[j] * Sqrt2);
    input += n2;
  }

  // inverse transform column vectors
  for (i = 0; i < n2; ++i) {
    IHT_16x8[tx_type].cols(out[i], out[i]);
  }

  maybe_flip_strides(&dest, &stride, &outp, &outstride, tx_type, n, n2);

  // Sum with the destination
  for (i = 0; i < n; ++i) {
    for (j = 0; j < n2; ++j) {
      int d = i * stride + j;
      int s = j * outstride + i;
      dest[d] = clip_pixel_add(dest[d], ROUND_POWER_OF_TWO(outp[s], 6));
    }
  }
}

void av1_iht16x32_512_add_c(const tran_low_t *input, uint8_t *dest, int stride,
                            int tx_type) {
  static const transform_2d IHT_16x32[] = {
    { idct32_c, idct16_c },         // DCT_DCT
    { ihalfright32_c, idct16_c },   // ADST_DCT
    { idct32_c, iadst16_c },        // DCT_ADST
    { ihalfright32_c, iadst16_c },  // ADST_ADST
    { ihalfright32_c, idct16_c },   // FLIPADST_DCT
    { idct32_c, iadst16_c },        // DCT_FLIPADST
    { ihalfright32_c, iadst16_c },  // FLIPADST_FLIPADST
    { ihalfright32_c, iadst16_c },  // ADST_FLIPADST
    { ihalfright32_c, iadst16_c },  // FLIPADST_ADST
    { iidtx32_c, iidtx16_c },       // IDTX
    { idct32_c, iidtx16_c },        // V_DCT
    { iidtx32_c, idct16_c },        // H_DCT
    { ihalfright32_c, iidtx16_c },  // V_ADST
    { iidtx32_c, iadst16_c },       // H_ADST
    { ihalfright32_c, iidtx16_c },  // V_FLIPADST
    { iidtx32_c, iadst16_c },       // H_FLIPADST
  };

  const int n = 16;
  const int n2 = 32;
  int i, j;
  tran_low_t out[16][32], outtmp[16];
  tran_low_t *outp = &out[0][0];
  int outstride = n2;

  // inverse transform row vectors and transpose
  for (i = 0; i < n2; ++i) {
    IHT_16x32[tx_type].rows(input, outtmp);
    for (j = 0; j < n; ++j)
      out[j][i] = (tran_low_t)dct_const_round_shift(outtmp[j] * Sqrt2);
    input += n;
  }

  // inverse transform column vectors
  for (i = 0; i < n; ++i) {
    IHT_16x32[tx_type].cols(out[i], out[i]);
  }

  maybe_flip_strides(&dest, &stride, &outp, &outstride, tx_type, n2, n);

  // Sum with the destination
  for (i = 0; i < n2; ++i) {
    for (j = 0; j < n; ++j) {
      int d = i * stride + j;
      int s = j * outstride + i;
      dest[d] = clip_pixel_add(dest[d], ROUND_POWER_OF_TWO(outp[s], 6));
    }
  }
}

void av1_iht32x16_512_add_c(const tran_low_t *input, uint8_t *dest, int stride,
                            int tx_type) {
  static const transform_2d IHT_32x16[] = {
    { idct16_c, idct32_c },         // DCT_DCT
    { iadst16_c, idct32_c },        // ADST_DCT
    { idct16_c, ihalfright32_c },   // DCT_ADST
    { iadst16_c, ihalfright32_c },  // ADST_ADST
    { iadst16_c, idct32_c },        // FLIPADST_DCT
    { idct16_c, ihalfright32_c },   // DCT_FLIPADST
    { iadst16_c, ihalfright32_c },  // FLIPADST_FLIPADST
    { iadst16_c, ihalfright32_c },  // ADST_FLIPADST
    { iadst16_c, ihalfright32_c },  // FLIPADST_ADST
    { iidtx16_c, iidtx32_c },       // IDTX
    { idct16_c, iidtx32_c },        // V_DCT
    { iidtx16_c, idct32_c },        // H_DCT
    { iadst16_c, iidtx32_c },       // V_ADST
    { iidtx16_c, ihalfright32_c },  // H_ADST
    { iadst16_c, iidtx32_c },       // V_FLIPADST
    { iidtx16_c, ihalfright32_c },  // H_FLIPADST
  };
  const int n = 16;
  const int n2 = 32;

  int i, j;
  tran_low_t out[32][16], outtmp[32];
  tran_low_t *outp = &out[0][0];
  int outstride = n;

  // inverse transform row vectors and transpose
  for (i = 0; i < n; ++i) {
    IHT_32x16[tx_type].rows(input, outtmp);
    for (j = 0; j < n2; ++j)
      out[j][i] = (tran_low_t)dct_const_round_shift(outtmp[j] * Sqrt2);
    input += n2;
  }

  // inverse transform column vectors
  for (i = 0; i < n2; ++i) {
    IHT_32x16[tx_type].cols(out[i], out[i]);
  }

  maybe_flip_strides(&dest, &stride, &outp, &outstride, tx_type, n, n2);

  // Sum with the destination
  for (i = 0; i < n; ++i) {
    for (j = 0; j < n2; ++j) {
      int d = i * stride + j;
      int s = j * outstride + i;
      dest[d] = clip_pixel_add(dest[d], ROUND_POWER_OF_TWO(outp[s], 6));
    }
  }
}
#endif  // CONFIG_EXT_TX

void av1_iht8x8_64_add_c(const tran_low_t *input, uint8_t *dest, int stride,
                         int tx_type) {
  static const transform_2d IHT_8[] = {
    { idct8_c, idct8_c },    // DCT_DCT
    { iadst8_c, idct8_c },   // ADST_DCT
    { idct8_c, iadst8_c },   // DCT_ADST
    { iadst8_c, iadst8_c },  // ADST_ADST
#if CONFIG_EXT_TX
    { iadst8_c, idct8_c },   // FLIPADST_DCT
    { idct8_c, iadst8_c },   // DCT_FLIPADST
    { iadst8_c, iadst8_c },  // FLIPADST_FLIPADST
    { iadst8_c, iadst8_c },  // ADST_FLIPADST
    { iadst8_c, iadst8_c },  // FLIPADST_ADST
    { iidtx8_c, iidtx8_c },  // IDTX
    { idct8_c, iidtx8_c },   // V_DCT
    { iidtx8_c, idct8_c },   // H_DCT
    { iadst8_c, iidtx8_c },  // V_ADST
    { iidtx8_c, iadst8_c },  // H_ADST
    { iadst8_c, iidtx8_c },  // V_FLIPADST
    { iidtx8_c, iadst8_c },  // H_FLIPADST
#endif                       // CONFIG_EXT_TX
  };

  int i, j;
  tran_low_t tmp;
  tran_low_t out[8][8];
  tran_low_t *outp = &out[0][0];
  int outstride = 8;

  // inverse transform row vectors
  for (i = 0; i < 8; ++i) {
    IHT_8[tx_type].rows(input, out[i]);
    input += 8;
  }

  // transpose
  for (i = 1; i < 8; i++) {
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
  maybe_flip_strides(&dest, &stride, &outp, &outstride, tx_type, 8, 8);
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

void av1_iht16x16_256_add_c(const tran_low_t *input, uint8_t *dest, int stride,
                            int tx_type) {
  static const transform_2d IHT_16[] = {
    { idct16_c, idct16_c },    // DCT_DCT
    { iadst16_c, idct16_c },   // ADST_DCT
    { idct16_c, iadst16_c },   // DCT_ADST
    { iadst16_c, iadst16_c },  // ADST_ADST
#if CONFIG_EXT_TX
    { iadst16_c, idct16_c },   // FLIPADST_DCT
    { idct16_c, iadst16_c },   // DCT_FLIPADST
    { iadst16_c, iadst16_c },  // FLIPADST_FLIPADST
    { iadst16_c, iadst16_c },  // ADST_FLIPADST
    { iadst16_c, iadst16_c },  // FLIPADST_ADST
    { iidtx16_c, iidtx16_c },  // IDTX
    { idct16_c, iidtx16_c },   // V_DCT
    { iidtx16_c, idct16_c },   // H_DCT
    { iadst16_c, iidtx16_c },  // V_ADST
    { iidtx16_c, iadst16_c },  // H_ADST
    { iadst16_c, iidtx16_c },  // V_FLIPADST
    { iidtx16_c, iadst16_c },  // H_FLIPADST
#endif                         // CONFIG_EXT_TX
  };

  int i, j;
  tran_low_t tmp;
  tran_low_t out[16][16];
  tran_low_t *outp = &out[0][0];
  int outstride = 16;

  // inverse transform row vectors
  for (i = 0; i < 16; ++i) {
    IHT_16[tx_type].rows(input, out[i]);
    input += 16;
  }

  // transpose
  for (i = 1; i < 16; i++) {
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
  maybe_flip_strides(&dest, &stride, &outp, &outstride, tx_type, 16, 16);
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
void av1_iht32x32_1024_add_c(const tran_low_t *input, uint8_t *dest, int stride,
                             int tx_type) {
  static const transform_2d IHT_32[] = {
    { idct32_c, idct32_c },              // DCT_DCT
    { ihalfright32_c, idct32_c },        // ADST_DCT
    { idct32_c, ihalfright32_c },        // DCT_ADST
    { ihalfright32_c, ihalfright32_c },  // ADST_ADST
    { ihalfright32_c, idct32_c },        // FLIPADST_DCT
    { idct32_c, ihalfright32_c },        // DCT_FLIPADST
    { ihalfright32_c, ihalfright32_c },  // FLIPADST_FLIPADST
    { ihalfright32_c, ihalfright32_c },  // ADST_FLIPADST
    { ihalfright32_c, ihalfright32_c },  // FLIPADST_ADST
    { iidtx32_c, iidtx32_c },            // IDTX
    { idct32_c, iidtx32_c },             // V_DCT
    { iidtx32_c, idct32_c },             // H_DCT
    { ihalfright32_c, iidtx16_c },       // V_ADST
    { iidtx16_c, ihalfright32_c },       // H_ADST
    { ihalfright32_c, iidtx16_c },       // V_FLIPADST
    { iidtx16_c, ihalfright32_c },       // H_FLIPADST
  };

  int i, j;
  tran_low_t tmp;
  tran_low_t out[32][32];
  tran_low_t *outp = &out[0][0];
  int outstride = 32;

  // inverse transform row vectors
  for (i = 0; i < 32; ++i) {
    IHT_32[tx_type].rows(input, out[i]);
    input += 32;
  }

  // transpose
  for (i = 1; i < 32; i++) {
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

  maybe_flip_strides(&dest, &stride, &outp, &outstride, tx_type, 32, 32);

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
void av1_idct4x4_add(const tran_low_t *input, uint8_t *dest, int stride,
                     int eob) {
  if (eob > 1)
    aom_idct4x4_16_add(input, dest, stride);
  else
    aom_idct4x4_1_add(input, dest, stride);
}

void av1_iwht4x4_add(const tran_low_t *input, uint8_t *dest, int stride,
                     int eob) {
  if (eob > 1)
    aom_iwht4x4_16_add(input, dest, stride);
  else
    aom_iwht4x4_1_add(input, dest, stride);
}

void av1_idct8x8_add(const tran_low_t *input, uint8_t *dest, int stride,
                     int eob) {
  // If dc is 1, then input[0] is the reconstructed value, do not need
  // dequantization. Also, when dc is 1, dc is counted in eobs, namely eobs >=1.

  // The calculation can be simplified if there are not many non-zero dct
  // coefficients. Use eobs to decide what to do.
  // TODO(yunqingwang): "eobs = 1" case is also handled in av1_short_idct8x8_c.
  // Combine that with code here.
  if (eob == 1)
    // DC only DCT coefficient
    aom_idct8x8_1_add(input, dest, stride);
  else if (eob <= 12)
    aom_idct8x8_12_add(input, dest, stride);
  else
    aom_idct8x8_64_add(input, dest, stride);
}

void av1_idct16x16_add(const tran_low_t *input, uint8_t *dest, int stride,
                       int eob) {
  /* The calculation can be simplified if there are not many non-zero dct
   * coefficients. Use eobs to separate different cases. */
  if (eob == 1) /* DC only DCT coefficient. */
    aom_idct16x16_1_add(input, dest, stride);
  else if (eob <= 10)
    aom_idct16x16_10_add(input, dest, stride);
  else
    aom_idct16x16_256_add(input, dest, stride);
}

void av1_idct32x32_add(const tran_low_t *input, uint8_t *dest, int stride,
                       int eob) {
  if (eob == 1)
    aom_idct32x32_1_add(input, dest, stride);
  else if (eob <= 34)
    // non-zero coeff only in upper-left 8x8
    aom_idct32x32_34_add(input, dest, stride);
  else
    aom_idct32x32_1024_add(input, dest, stride);
}

void av1_inv_txfm_add_4x4(const tran_low_t *input, uint8_t *dest, int stride,
                          int eob, TX_TYPE tx_type, int lossless) {
  if (lossless) {
    assert(tx_type == DCT_DCT);
    av1_iwht4x4_add(input, dest, stride, eob);
    return;
  }

  switch (tx_type) {
    case DCT_DCT: av1_idct4x4_add(input, dest, stride, eob); break;
    case ADST_DCT:
    case DCT_ADST:
    case ADST_ADST: av1_iht4x4_16_add(input, dest, stride, tx_type); break;
#if CONFIG_EXT_TX
    case FLIPADST_DCT:
    case DCT_FLIPADST:
    case FLIPADST_FLIPADST:
    case ADST_FLIPADST:
    case FLIPADST_ADST: av1_iht4x4_16_add(input, dest, stride, tx_type); break;
    case V_DCT:
    case H_DCT:
    case V_ADST:
    case H_ADST:
    case V_FLIPADST:
    case H_FLIPADST:
      // Use C version since DST only exists in C code
      av1_iht4x4_16_add_c(input, dest, stride, tx_type);
      break;
    case IDTX: inv_idtx_add_c(input, dest, stride, 4, tx_type); break;
#endif  // CONFIG_EXT_TX
    default: assert(0); break;
  }
}

#if CONFIG_EXT_TX
void av1_inv_txfm_add_4x8(const tran_low_t *input, uint8_t *dest, int stride,
                          int eob, TX_TYPE tx_type) {
  (void)eob;
  av1_iht4x8_32_add(input, dest, stride, tx_type);
}

void av1_inv_txfm_add_8x4(const tran_low_t *input, uint8_t *dest, int stride,
                          int eob, TX_TYPE tx_type) {
  (void)eob;
  av1_iht8x4_32_add(input, dest, stride, tx_type);
}

void av1_inv_txfm_add_8x16(const tran_low_t *input, uint8_t *dest, int stride,
                           int eob, TX_TYPE tx_type) {
  (void)eob;
  av1_iht8x16_128_add(input, dest, stride, tx_type);
}

void av1_inv_txfm_add_16x8(const tran_low_t *input, uint8_t *dest, int stride,
                           int eob, TX_TYPE tx_type) {
  (void)eob;
  av1_iht16x8_128_add(input, dest, stride, tx_type);
}

void av1_inv_txfm_add_16x32(const tran_low_t *input, uint8_t *dest, int stride,
                            int eob, TX_TYPE tx_type) {
  (void)eob;
  av1_iht16x32_512_add(input, dest, stride, tx_type);
}

void av1_inv_txfm_add_32x16(const tran_low_t *input, uint8_t *dest, int stride,
                            int eob, TX_TYPE tx_type) {
  (void)eob;
  av1_iht32x16_512_add(input, dest, stride, tx_type);
}
#endif  // CONFIG_EXT_TX

void av1_inv_txfm_add_8x8(const tran_low_t *input, uint8_t *dest, int stride,
                          int eob, TX_TYPE tx_type) {
  switch (tx_type) {
    case DCT_DCT: av1_idct8x8_add(input, dest, stride, eob); break;
    case ADST_DCT:
    case DCT_ADST:
    case ADST_ADST: av1_iht8x8_64_add(input, dest, stride, tx_type); break;
#if CONFIG_EXT_TX
    case FLIPADST_DCT:
    case DCT_FLIPADST:
    case FLIPADST_FLIPADST:
    case ADST_FLIPADST:
    case FLIPADST_ADST: av1_iht8x8_64_add(input, dest, stride, tx_type); break;
    case V_DCT:
    case H_DCT:
    case V_ADST:
    case H_ADST:
    case V_FLIPADST:
    case H_FLIPADST:
      // Use C version since DST only exists in C code
      av1_iht8x8_64_add_c(input, dest, stride, tx_type);
      break;
    case IDTX: inv_idtx_add_c(input, dest, stride, 8, tx_type); break;
#endif  // CONFIG_EXT_TX
    default: assert(0); break;
  }
}

void av1_inv_txfm_add_16x16(const tran_low_t *input, uint8_t *dest, int stride,
                            int eob, TX_TYPE tx_type) {
  switch (tx_type) {
    case DCT_DCT: av1_idct16x16_add(input, dest, stride, eob); break;
    case ADST_DCT:
    case DCT_ADST:
    case ADST_ADST: av1_iht16x16_256_add(input, dest, stride, tx_type); break;
#if CONFIG_EXT_TX
    case FLIPADST_DCT:
    case DCT_FLIPADST:
    case FLIPADST_FLIPADST:
    case ADST_FLIPADST:
    case FLIPADST_ADST:
      av1_iht16x16_256_add(input, dest, stride, tx_type);
      break;
    case V_DCT:
    case H_DCT:
    case V_ADST:
    case H_ADST:
    case V_FLIPADST:
    case H_FLIPADST:
      // Use C version since DST only exists in C code
      av1_iht16x16_256_add_c(input, dest, stride, tx_type);
      break;
    case IDTX: inv_idtx_add_c(input, dest, stride, 16, tx_type); break;
#endif  // CONFIG_EXT_TX
    default: assert(0); break;
  }
}

void av1_inv_txfm_add_32x32(const tran_low_t *input, uint8_t *dest, int stride,
                            int eob, TX_TYPE tx_type) {
  switch (tx_type) {
    case DCT_DCT: av1_idct32x32_add(input, dest, stride, eob); break;
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
      av1_iht32x32_1024_add_c(input, dest, stride, tx_type);
      break;
    case IDTX: inv_idtx_add_c(input, dest, stride, 32, tx_type); break;
#endif  // CONFIG_EXT_TX
    default: assert(0); break;
  }
}

#if CONFIG_AOM_HIGHBITDEPTH
void av1_highbd_iht4x4_16_add_c(const tran_low_t *input, uint8_t *dest8,
                                int stride, int tx_type, int bd) {
  static const highbd_transform_2d HIGH_IHT_4[] = {
    { aom_highbd_idct4_c, aom_highbd_idct4_c },    // DCT_DCT
    { aom_highbd_iadst4_c, aom_highbd_idct4_c },   // ADST_DCT
    { aom_highbd_idct4_c, aom_highbd_iadst4_c },   // DCT_ADST
    { aom_highbd_iadst4_c, aom_highbd_iadst4_c },  // ADST_ADST
#if CONFIG_EXT_TX
    { aom_highbd_iadst4_c, aom_highbd_idct4_c },   // FLIPADST_DCT
    { aom_highbd_idct4_c, aom_highbd_iadst4_c },   // DCT_FLIPADST
    { aom_highbd_iadst4_c, aom_highbd_iadst4_c },  // FLIPADST_FLIPADST
    { aom_highbd_iadst4_c, aom_highbd_iadst4_c },  // ADST_FLIPADST
    { aom_highbd_iadst4_c, aom_highbd_iadst4_c },  // FLIPADST_ADST
    { highbd_iidtx4_c, highbd_iidtx4_c },          // IDTX
    { aom_highbd_idct4_c, highbd_iidtx4_c },       // V_DCT
    { highbd_iidtx4_c, aom_highbd_idct4_c },       // H_DCT
    { aom_highbd_iadst4_c, highbd_iidtx4_c },      // V_ADST
    { highbd_iidtx4_c, aom_highbd_iadst4_c },      // H_ADST
    { aom_highbd_iadst4_c, highbd_iidtx4_c },      // V_FLIPADST
    { highbd_iidtx4_c, aom_highbd_iadst4_c },      // H_FLIPADST
#endif                                             // CONFIG_EXT_TX
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
    input += 4;
  }

  // transpose
  for (i = 1; i < 4; i++) {
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
  maybe_flip_strides16(&dest, &stride, &outp, &outstride, tx_type, 4, 4);
#endif

  // Sum with the destination
  for (i = 0; i < 4; ++i) {
    for (j = 0; j < 4; ++j) {
      int d = i * stride + j;
      int s = j * outstride + i;
      dest[d] =
          highbd_clip_pixel_add(dest[d], ROUND_POWER_OF_TWO(outp[s], 4), bd);
    }
  }
}

#if CONFIG_EXT_TX
void av1_highbd_iht4x8_32_add_c(const tran_low_t *input, uint8_t *dest8,
                                int stride, int tx_type, int bd) {
  static const highbd_transform_2d HIGH_IHT_4x8[] = {
    { aom_highbd_idct8_c, aom_highbd_idct4_c },    // DCT_DCT
    { aom_highbd_iadst8_c, aom_highbd_idct4_c },   // ADST_DCT
    { aom_highbd_idct8_c, aom_highbd_iadst4_c },   // DCT_ADST
    { aom_highbd_iadst8_c, aom_highbd_iadst4_c },  // ADST_ADST
    { aom_highbd_iadst8_c, aom_highbd_idct4_c },   // FLIPADST_DCT
    { aom_highbd_idct8_c, aom_highbd_iadst4_c },   // DCT_FLIPADST
    { aom_highbd_iadst8_c, aom_highbd_iadst4_c },  // FLIPADST_FLIPADST
    { aom_highbd_iadst8_c, aom_highbd_iadst4_c },  // ADST_FLIPADST
    { aom_highbd_iadst8_c, aom_highbd_iadst4_c },  // FLIPADST_ADST
    { highbd_iidtx8_c, highbd_iidtx4_c },          // IDTX
    { aom_highbd_idct8_c, highbd_iidtx4_c },       // V_DCT
    { highbd_iidtx8_c, aom_highbd_idct4_c },       // H_DCT
    { aom_highbd_iadst8_c, highbd_iidtx4_c },      // V_ADST
    { highbd_iidtx8_c, aom_highbd_iadst4_c },      // H_ADST
    { aom_highbd_iadst8_c, highbd_iidtx4_c },      // V_FLIPADST
    { highbd_iidtx8_c, aom_highbd_iadst4_c },      // H_FLIPADST
  };
  const int n = 4;
  const int n2 = 8;

  uint16_t *dest = CONVERT_TO_SHORTPTR(dest8);

  int i, j;
  tran_low_t out[4][8], outtmp[4];
  tran_low_t *outp = &out[0][0];
  int outstride = n2;

  // inverse transform row vectors, and transpose
  for (i = 0; i < n2; ++i) {
    HIGH_IHT_4x8[tx_type].rows(input, outtmp, bd);
    for (j = 0; j < n; ++j) {
      out[j][i] =
          HIGHBD_WRAPLOW(highbd_dct_const_round_shift(outtmp[j] * Sqrt2), bd);
    }
    input += n;
  }

  // inverse transform column vectors
  for (i = 0; i < n; ++i) {
    HIGH_IHT_4x8[tx_type].cols(out[i], out[i], bd);
  }

  maybe_flip_strides16(&dest, &stride, &outp, &outstride, tx_type, n2, n);

  // Sum with the destination
  for (i = 0; i < n2; ++i) {
    for (j = 0; j < n; ++j) {
      int d = i * stride + j;
      int s = j * outstride + i;
      dest[d] =
          highbd_clip_pixel_add(dest[d], ROUND_POWER_OF_TWO(outp[s], 5), bd);
    }
  }
}

void av1_highbd_iht8x4_32_add_c(const tran_low_t *input, uint8_t *dest8,
                                int stride, int tx_type, int bd) {
  static const highbd_transform_2d HIGH_IHT_8x4[] = {
    { aom_highbd_idct4_c, aom_highbd_idct8_c },    // DCT_DCT
    { aom_highbd_iadst4_c, aom_highbd_idct8_c },   // ADST_DCT
    { aom_highbd_idct4_c, aom_highbd_iadst8_c },   // DCT_ADST
    { aom_highbd_iadst4_c, aom_highbd_iadst8_c },  // ADST_ADST
    { aom_highbd_iadst4_c, aom_highbd_idct8_c },   // FLIPADST_DCT
    { aom_highbd_idct4_c, aom_highbd_iadst8_c },   // DCT_FLIPADST
    { aom_highbd_iadst4_c, aom_highbd_iadst8_c },  // FLIPADST_FLIPADST
    { aom_highbd_iadst4_c, aom_highbd_iadst8_c },  // ADST_FLIPADST
    { aom_highbd_iadst4_c, aom_highbd_iadst8_c },  // FLIPADST_ADST
    { highbd_iidtx4_c, highbd_iidtx8_c },          // IDTX
    { aom_highbd_idct4_c, highbd_iidtx8_c },       // V_DCT
    { highbd_iidtx4_c, aom_highbd_idct8_c },       // H_DCT
    { aom_highbd_iadst4_c, highbd_iidtx8_c },      // V_ADST
    { highbd_iidtx4_c, aom_highbd_iadst8_c },      // H_ADST
    { aom_highbd_iadst4_c, highbd_iidtx8_c },      // V_FLIPADST
    { highbd_iidtx4_c, aom_highbd_iadst8_c },      // H_FLIPADST
  };
  const int n = 4;
  const int n2 = 8;

  uint16_t *dest = CONVERT_TO_SHORTPTR(dest8);

  int i, j;
  tran_low_t out[8][4], outtmp[8];
  tran_low_t *outp = &out[0][0];
  int outstride = n;

  // inverse transform row vectors, and transpose
  for (i = 0; i < n; ++i) {
    HIGH_IHT_8x4[tx_type].rows(input, outtmp, bd);
    for (j = 0; j < n2; ++j) {
      out[j][i] =
          HIGHBD_WRAPLOW(highbd_dct_const_round_shift(outtmp[j] * Sqrt2), bd);
    }
    input += n2;
  }

  // inverse transform column vectors
  for (i = 0; i < n2; ++i) {
    HIGH_IHT_8x4[tx_type].cols(out[i], out[i], bd);
  }

  maybe_flip_strides16(&dest, &stride, &outp, &outstride, tx_type, n, n2);

  // Sum with the destination
  for (i = 0; i < n; ++i) {
    for (j = 0; j < n2; ++j) {
      int d = i * stride + j;
      int s = j * outstride + i;
      dest[d] =
          highbd_clip_pixel_add(dest[d], ROUND_POWER_OF_TWO(outp[s], 5), bd);
    }
  }
}

void av1_highbd_iht8x16_128_add_c(const tran_low_t *input, uint8_t *dest8,
                                  int stride, int tx_type, int bd) {
  static const highbd_transform_2d HIGH_IHT_8x16[] = {
    { aom_highbd_idct16_c, aom_highbd_idct8_c },    // DCT_DCT
    { aom_highbd_iadst16_c, aom_highbd_idct8_c },   // ADST_DCT
    { aom_highbd_idct16_c, aom_highbd_iadst8_c },   // DCT_ADST
    { aom_highbd_iadst16_c, aom_highbd_iadst8_c },  // ADST_ADST
    { aom_highbd_iadst16_c, aom_highbd_idct8_c },   // FLIPADST_DCT
    { aom_highbd_idct16_c, aom_highbd_iadst8_c },   // DCT_FLIPADST
    { aom_highbd_iadst16_c, aom_highbd_iadst8_c },  // FLIPADST_FLIPADST
    { aom_highbd_iadst16_c, aom_highbd_iadst8_c },  // ADST_FLIPADST
    { aom_highbd_iadst16_c, aom_highbd_iadst8_c },  // FLIPADST_ADST
    { highbd_iidtx16_c, highbd_iidtx8_c },          // IDTX
    { aom_highbd_idct16_c, highbd_iidtx8_c },       // V_DCT
    { highbd_iidtx16_c, aom_highbd_idct8_c },       // H_DCT
    { aom_highbd_iadst16_c, highbd_iidtx8_c },      // V_ADST
    { highbd_iidtx16_c, aom_highbd_iadst8_c },      // H_ADST
    { aom_highbd_iadst16_c, highbd_iidtx8_c },      // V_FLIPADST
    { highbd_iidtx16_c, aom_highbd_iadst8_c },      // H_FLIPADST
  };
  const int n = 8;
  const int n2 = 16;

  uint16_t *dest = CONVERT_TO_SHORTPTR(dest8);

  int i, j;
  tran_low_t out[8][16], outtmp[8];
  tran_low_t *outp = &out[0][0];
  int outstride = n2;

  // inverse transform row vectors, and transpose
  for (i = 0; i < n2; ++i) {
    HIGH_IHT_8x16[tx_type].rows(input, outtmp, bd);
    for (j = 0; j < n; ++j)
      out[j][i] =
          HIGHBD_WRAPLOW(highbd_dct_const_round_shift(outtmp[j] * Sqrt2), bd);
    input += n;
  }

  // inverse transform column vectors
  for (i = 0; i < n; ++i) {
    HIGH_IHT_8x16[tx_type].cols(out[i], out[i], bd);
  }

  maybe_flip_strides16(&dest, &stride, &outp, &outstride, tx_type, n2, n);

  // Sum with the destination
  for (i = 0; i < n2; ++i) {
    for (j = 0; j < n; ++j) {
      int d = i * stride + j;
      int s = j * outstride + i;
      dest[d] =
          highbd_clip_pixel_add(dest[d], ROUND_POWER_OF_TWO(outp[s], 6), bd);
    }
  }
}

void av1_highbd_iht16x8_128_add_c(const tran_low_t *input, uint8_t *dest8,
                                  int stride, int tx_type, int bd) {
  static const highbd_transform_2d HIGH_IHT_16x8[] = {
    { aom_highbd_idct8_c, aom_highbd_idct16_c },    // DCT_DCT
    { aom_highbd_iadst8_c, aom_highbd_idct16_c },   // ADST_DCT
    { aom_highbd_idct8_c, aom_highbd_iadst16_c },   // DCT_ADST
    { aom_highbd_iadst8_c, aom_highbd_iadst16_c },  // ADST_ADST
    { aom_highbd_iadst8_c, aom_highbd_idct16_c },   // FLIPADST_DCT
    { aom_highbd_idct8_c, aom_highbd_iadst16_c },   // DCT_FLIPADST
    { aom_highbd_iadst8_c, aom_highbd_iadst16_c },  // FLIPADST_FLIPADST
    { aom_highbd_iadst8_c, aom_highbd_iadst16_c },  // ADST_FLIPADST
    { aom_highbd_iadst8_c, aom_highbd_iadst16_c },  // FLIPADST_ADST
    { highbd_iidtx8_c, highbd_iidtx16_c },          // IDTX
    { aom_highbd_idct8_c, highbd_iidtx16_c },       // V_DCT
    { highbd_iidtx8_c, aom_highbd_idct16_c },       // H_DCT
    { aom_highbd_iadst8_c, highbd_iidtx16_c },      // V_ADST
    { highbd_iidtx8_c, aom_highbd_iadst16_c },      // H_ADST
    { aom_highbd_iadst8_c, highbd_iidtx16_c },      // V_FLIPADST
    { highbd_iidtx8_c, aom_highbd_iadst16_c },      // H_FLIPADST
  };
  const int n = 8;
  const int n2 = 16;

  uint16_t *dest = CONVERT_TO_SHORTPTR(dest8);

  int i, j;
  tran_low_t out[16][8], outtmp[16];
  tran_low_t *outp = &out[0][0];
  int outstride = n;

  // inverse transform row vectors, and transpose
  for (i = 0; i < n; ++i) {
    HIGH_IHT_16x8[tx_type].rows(input, outtmp, bd);
    for (j = 0; j < n2; ++j)
      out[j][i] =
          HIGHBD_WRAPLOW(highbd_dct_const_round_shift(outtmp[j] * Sqrt2), bd);
    input += n2;
  }

  // inverse transform column vectors
  for (i = 0; i < n2; ++i) {
    HIGH_IHT_16x8[tx_type].cols(out[i], out[i], bd);
  }

  maybe_flip_strides16(&dest, &stride, &outp, &outstride, tx_type, n, n2);

  // Sum with the destination
  for (i = 0; i < n; ++i) {
    for (j = 0; j < n2; ++j) {
      int d = i * stride + j;
      int s = j * outstride + i;
      dest[d] =
          highbd_clip_pixel_add(dest[d], ROUND_POWER_OF_TWO(outp[s], 6), bd);
    }
  }
}

void av1_highbd_iht16x32_512_add_c(const tran_low_t *input, uint8_t *dest8,
                                   int stride, int tx_type, int bd) {
  static const highbd_transform_2d HIGH_IHT_16x32[] = {
    { aom_highbd_idct32_c, aom_highbd_idct16_c },     // DCT_DCT
    { highbd_ihalfright32_c, aom_highbd_idct16_c },   // ADST_DCT
    { aom_highbd_idct32_c, aom_highbd_iadst16_c },    // DCT_ADST
    { highbd_ihalfright32_c, aom_highbd_iadst16_c },  // ADST_ADST
    { highbd_ihalfright32_c, aom_highbd_idct16_c },   // FLIPADST_DCT
    { aom_highbd_idct32_c, aom_highbd_iadst16_c },    // DCT_FLIPADST
    { highbd_ihalfright32_c, aom_highbd_iadst16_c },  // FLIPADST_FLIPADST
    { highbd_ihalfright32_c, aom_highbd_iadst16_c },  // ADST_FLIPADST
    { highbd_ihalfright32_c, aom_highbd_iadst16_c },  // FLIPADST_ADST
    { highbd_iidtx32_c, highbd_iidtx16_c },           // IDTX
    { aom_highbd_idct32_c, highbd_iidtx16_c },        // V_DCT
    { highbd_iidtx32_c, aom_highbd_idct16_c },        // H_DCT
    { highbd_ihalfright32_c, highbd_iidtx16_c },      // V_ADST
    { highbd_iidtx32_c, aom_highbd_iadst16_c },       // H_ADST
    { highbd_ihalfright32_c, highbd_iidtx16_c },      // V_FLIPADST
    { highbd_iidtx32_c, aom_highbd_iadst16_c },       // H_FLIPADST
  };
  const int n = 16;
  const int n2 = 32;

  uint16_t *dest = CONVERT_TO_SHORTPTR(dest8);

  int i, j;
  tran_low_t out[16][32], outtmp[16];
  tran_low_t *outp = &out[0][0];
  int outstride = n2;

  // inverse transform row vectors, and transpose
  for (i = 0; i < n2; ++i) {
    HIGH_IHT_16x32[tx_type].rows(input, outtmp, bd);
    for (j = 0; j < n; ++j)
      out[j][i] =
          HIGHBD_WRAPLOW(highbd_dct_const_round_shift(outtmp[j] * Sqrt2), bd);
    input += n;
  }

  // inverse transform column vectors
  for (i = 0; i < n; ++i) {
    HIGH_IHT_16x32[tx_type].cols(out[i], out[i], bd);
  }

  maybe_flip_strides16(&dest, &stride, &outp, &outstride, tx_type, n2, n);

  // Sum with the destination
  for (i = 0; i < n2; ++i) {
    for (j = 0; j < n; ++j) {
      int d = i * stride + j;
      int s = j * outstride + i;
      dest[d] =
          highbd_clip_pixel_add(dest[d], ROUND_POWER_OF_TWO(outp[s], 6), bd);
    }
  }
}

void av1_highbd_iht32x16_512_add_c(const tran_low_t *input, uint8_t *dest8,
                                   int stride, int tx_type, int bd) {
  static const highbd_transform_2d HIGH_IHT_32x16[] = {
    { aom_highbd_idct16_c, aom_highbd_idct32_c },     // DCT_DCT
    { aom_highbd_iadst16_c, aom_highbd_idct32_c },    // ADST_DCT
    { aom_highbd_idct16_c, highbd_ihalfright32_c },   // DCT_ADST
    { aom_highbd_iadst16_c, highbd_ihalfright32_c },  // ADST_ADST
    { aom_highbd_iadst16_c, aom_highbd_idct32_c },    // FLIPADST_DCT
    { aom_highbd_idct16_c, highbd_ihalfright32_c },   // DCT_FLIPADST
    { aom_highbd_iadst16_c, highbd_ihalfright32_c },  // FLIPADST_FLIPADST
    { aom_highbd_iadst16_c, highbd_ihalfright32_c },  // ADST_FLIPADST
    { aom_highbd_iadst16_c, highbd_ihalfright32_c },  // FLIPADST_ADST
    { highbd_iidtx16_c, highbd_iidtx32_c },           // IDTX
    { aom_highbd_idct16_c, highbd_iidtx32_c },        // V_DCT
    { highbd_iidtx16_c, aom_highbd_idct32_c },        // H_DCT
    { aom_highbd_iadst16_c, highbd_iidtx32_c },       // V_ADST
    { highbd_iidtx16_c, highbd_ihalfright32_c },      // H_ADST
    { aom_highbd_iadst16_c, highbd_iidtx32_c },       // V_FLIPADST
    { highbd_iidtx16_c, highbd_ihalfright32_c },      // H_FLIPADST
  };
  const int n = 16;
  const int n2 = 32;

  uint16_t *dest = CONVERT_TO_SHORTPTR(dest8);

  int i, j;
  tran_low_t out[32][16], outtmp[32];
  tran_low_t *outp = &out[0][0];
  int outstride = n;

  // inverse transform row vectors, and transpose
  for (i = 0; i < n; ++i) {
    HIGH_IHT_32x16[tx_type].rows(input, outtmp, bd);
    for (j = 0; j < n2; ++j)
      out[j][i] =
          HIGHBD_WRAPLOW(highbd_dct_const_round_shift(outtmp[j] * Sqrt2), bd);
    input += n2;
  }

  // inverse transform column vectors
  for (i = 0; i < n2; ++i) {
    HIGH_IHT_32x16[tx_type].cols(out[i], out[i], bd);
  }

  maybe_flip_strides16(&dest, &stride, &outp, &outstride, tx_type, n, n2);

  // Sum with the destination
  for (i = 0; i < n; ++i) {
    for (j = 0; j < n2; ++j) {
      int d = i * stride + j;
      int s = j * outstride + i;
      dest[d] =
          highbd_clip_pixel_add(dest[d], ROUND_POWER_OF_TWO(outp[s], 6), bd);
    }
  }
}
#endif  // CONFIG_EXT_TX

void av1_highbd_iht8x8_64_add_c(const tran_low_t *input, uint8_t *dest8,
                                int stride, int tx_type, int bd) {
  static const highbd_transform_2d HIGH_IHT_8[] = {
    { aom_highbd_idct8_c, aom_highbd_idct8_c },    // DCT_DCT
    { aom_highbd_iadst8_c, aom_highbd_idct8_c },   // ADST_DCT
    { aom_highbd_idct8_c, aom_highbd_iadst8_c },   // DCT_ADST
    { aom_highbd_iadst8_c, aom_highbd_iadst8_c },  // ADST_ADST
#if CONFIG_EXT_TX
    { aom_highbd_iadst8_c, aom_highbd_idct8_c },   // FLIPADST_DCT
    { aom_highbd_idct8_c, aom_highbd_iadst8_c },   // DCT_FLIPADST
    { aom_highbd_iadst8_c, aom_highbd_iadst8_c },  // FLIPADST_FLIPADST
    { aom_highbd_iadst8_c, aom_highbd_iadst8_c },  // ADST_FLIPADST
    { aom_highbd_iadst8_c, aom_highbd_iadst8_c },  // FLIPADST_ADST
    { highbd_iidtx8_c, highbd_iidtx8_c },          // IDTX
    { aom_highbd_idct8_c, highbd_iidtx8_c },       // V_DCT
    { highbd_iidtx8_c, aom_highbd_idct8_c },       // H_DCT
    { aom_highbd_iadst8_c, highbd_iidtx8_c },      // V_ADST
    { highbd_iidtx8_c, aom_highbd_iadst8_c },      // H_ADST
    { aom_highbd_iadst8_c, highbd_iidtx8_c },      // V_FLIPADST
    { highbd_iidtx8_c, aom_highbd_iadst8_c },      // H_FLIPADST
#endif                                             // CONFIG_EXT_TX
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
    input += 8;
  }

  // transpose
  for (i = 1; i < 8; i++) {
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
  maybe_flip_strides16(&dest, &stride, &outp, &outstride, tx_type, 8, 8);
#endif

  // Sum with the destination
  for (i = 0; i < 8; ++i) {
    for (j = 0; j < 8; ++j) {
      int d = i * stride + j;
      int s = j * outstride + i;
      dest[d] =
          highbd_clip_pixel_add(dest[d], ROUND_POWER_OF_TWO(outp[s], 5), bd);
    }
  }
}

void av1_highbd_iht16x16_256_add_c(const tran_low_t *input, uint8_t *dest8,
                                   int stride, int tx_type, int bd) {
  static const highbd_transform_2d HIGH_IHT_16[] = {
    { aom_highbd_idct16_c, aom_highbd_idct16_c },    // DCT_DCT
    { aom_highbd_iadst16_c, aom_highbd_idct16_c },   // ADST_DCT
    { aom_highbd_idct16_c, aom_highbd_iadst16_c },   // DCT_ADST
    { aom_highbd_iadst16_c, aom_highbd_iadst16_c },  // ADST_ADST
#if CONFIG_EXT_TX
    { aom_highbd_iadst16_c, aom_highbd_idct16_c },   // FLIPADST_DCT
    { aom_highbd_idct16_c, aom_highbd_iadst16_c },   // DCT_FLIPADST
    { aom_highbd_iadst16_c, aom_highbd_iadst16_c },  // FLIPADST_FLIPADST
    { aom_highbd_iadst16_c, aom_highbd_iadst16_c },  // ADST_FLIPADST
    { aom_highbd_iadst16_c, aom_highbd_iadst16_c },  // FLIPADST_ADST
    { highbd_iidtx16_c, highbd_iidtx16_c },          // IDTX
    { aom_highbd_idct16_c, highbd_iidtx16_c },       // V_DCT
    { highbd_iidtx16_c, aom_highbd_idct16_c },       // H_DCT
    { aom_highbd_iadst16_c, highbd_iidtx16_c },      // V_ADST
    { highbd_iidtx16_c, aom_highbd_iadst16_c },      // H_ADST
    { aom_highbd_iadst16_c, highbd_iidtx16_c },      // V_FLIPADST
    { highbd_iidtx16_c, aom_highbd_iadst16_c },      // H_FLIPADST
#endif                                               // CONFIG_EXT_TX
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
    input += 16;
  }

  // transpose
  for (i = 1; i < 16; i++) {
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
  maybe_flip_strides16(&dest, &stride, &outp, &outstride, tx_type, 16, 16);
#endif

  // Sum with the destination
  for (i = 0; i < 16; ++i) {
    for (j = 0; j < 16; ++j) {
      int d = i * stride + j;
      int s = j * outstride + i;
      dest[d] =
          highbd_clip_pixel_add(dest[d], ROUND_POWER_OF_TWO(outp[s], 6), bd);
    }
  }
}

#if CONFIG_EXT_TX
void av1_highbd_iht32x32_1024_add_c(const tran_low_t *input, uint8_t *dest8,
                                    int stride, int tx_type, int bd) {
  static const highbd_transform_2d HIGH_IHT_32[] = {
    { aom_highbd_idct32_c, aom_highbd_idct32_c },      // DCT_DCT
    { highbd_ihalfright32_c, aom_highbd_idct32_c },    // ADST_DCT
    { aom_highbd_idct32_c, highbd_ihalfright32_c },    // DCT_ADST
    { highbd_ihalfright32_c, highbd_ihalfright32_c },  // ADST_ADST
    { highbd_ihalfright32_c, aom_highbd_idct32_c },    // FLIPADST_DCT
    { aom_highbd_idct32_c, highbd_ihalfright32_c },    // DCT_FLIPADST
    { highbd_ihalfright32_c, highbd_ihalfright32_c },  // FLIPADST_FLIPADST
    { highbd_ihalfright32_c, highbd_ihalfright32_c },  // ADST_FLIPADST
    { highbd_ihalfright32_c, highbd_ihalfright32_c },  // FLIPADST_ADST
    { highbd_iidtx32_c, highbd_iidtx32_c },            // IDTX
    { aom_highbd_idct32_c, highbd_iidtx32_c },         // V_DCT
    { highbd_iidtx32_c, aom_highbd_idct32_c },         // H_DCT
    { highbd_ihalfright32_c, highbd_iidtx32_c },       // V_ADST
    { highbd_iidtx32_c, highbd_ihalfright32_c },       // H_ADST
    { highbd_ihalfright32_c, highbd_iidtx32_c },       // V_FLIPADST
    { highbd_iidtx32_c, highbd_ihalfright32_c },       // H_FLIPADST
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
    input += 32;
  }

  // transpose
  for (i = 1; i < 32; i++) {
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

  maybe_flip_strides16(&dest, &stride, &outp, &outstride, tx_type, 32, 32);

  // Sum with the destination
  for (i = 0; i < 32; ++i) {
    for (j = 0; j < 32; ++j) {
      int d = i * stride + j;
      int s = j * outstride + i;
      dest[d] =
          highbd_clip_pixel_add(dest[d], ROUND_POWER_OF_TWO(outp[s], 6), bd);
    }
  }
}
#endif  // CONFIG_EXT_TX

// idct
void av1_highbd_idct4x4_add(const tran_low_t *input, uint8_t *dest, int stride,
                            int eob, int bd) {
  if (eob > 1)
    aom_highbd_idct4x4_16_add(input, dest, stride, bd);
  else
    aom_highbd_idct4x4_1_add(input, dest, stride, bd);
}

void av1_highbd_iwht4x4_add(const tran_low_t *input, uint8_t *dest, int stride,
                            int eob, int bd) {
  if (eob > 1)
    aom_highbd_iwht4x4_16_add(input, dest, stride, bd);
  else
    aom_highbd_iwht4x4_1_add(input, dest, stride, bd);
}

void av1_highbd_idct8x8_add(const tran_low_t *input, uint8_t *dest, int stride,
                            int eob, int bd) {
  // If dc is 1, then input[0] is the reconstructed value, do not need
  // dequantization. Also, when dc is 1, dc is counted in eobs, namely eobs >=1.

  // The calculation can be simplified if there are not many non-zero dct
  // coefficients. Use eobs to decide what to do.
  // TODO(yunqingwang): "eobs = 1" case is also handled in av1_short_idct8x8_c.
  // Combine that with code here.
  // DC only DCT coefficient
  if (eob == 1) {
    aom_highbd_idct8x8_1_add(input, dest, stride, bd);
  } else if (eob <= 10) {
    aom_highbd_idct8x8_10_add(input, dest, stride, bd);
  } else {
    aom_highbd_idct8x8_64_add(input, dest, stride, bd);
  }
}

void av1_highbd_idct16x16_add(const tran_low_t *input, uint8_t *dest,
                              int stride, int eob, int bd) {
  // The calculation can be simplified if there are not many non-zero dct
  // coefficients. Use eobs to separate different cases.
  // DC only DCT coefficient.
  if (eob == 1) {
    aom_highbd_idct16x16_1_add(input, dest, stride, bd);
  } else if (eob <= 10) {
    aom_highbd_idct16x16_10_add(input, dest, stride, bd);
  } else {
    aom_highbd_idct16x16_256_add(input, dest, stride, bd);
  }
}

void av1_highbd_idct32x32_add(const tran_low_t *input, uint8_t *dest,
                              int stride, int eob, int bd) {
  // Non-zero coeff only in upper-left 8x8
  if (eob == 1) {
    aom_highbd_idct32x32_1_add(input, dest, stride, bd);
  } else if (eob <= 34) {
    aom_highbd_idct32x32_34_add(input, dest, stride, bd);
  } else {
    aom_highbd_idct32x32_1024_add(input, dest, stride, bd);
  }
}

void av1_highbd_inv_txfm_add_4x4(const tran_low_t *input, uint8_t *dest,
                                 int stride, int eob, int bd, TX_TYPE tx_type,
                                 int lossless) {
  if (lossless) {
    assert(tx_type == DCT_DCT);
    av1_highbd_iwht4x4_add(input, dest, stride, eob, bd);
    return;
  }

  switch (tx_type) {
    case DCT_DCT:
    case ADST_DCT:
    case DCT_ADST:
    case ADST_ADST:
      av1_inv_txfm2d_add_4x4(input, CONVERT_TO_SHORTPTR(dest), stride, tx_type,
                             bd);
      break;
#if CONFIG_EXT_TX
    case FLIPADST_DCT:
    case DCT_FLIPADST:
    case FLIPADST_FLIPADST:
    case ADST_FLIPADST:
    case FLIPADST_ADST:
      av1_inv_txfm2d_add_4x4(input, CONVERT_TO_SHORTPTR(dest), stride, tx_type,
                             bd);
      break;
    case V_DCT:
    case H_DCT:
    case V_ADST:
    case H_ADST:
    case V_FLIPADST:
    case H_FLIPADST:
      // Use C version since DST only exists in C code
      av1_highbd_iht4x4_16_add_c(input, dest, stride, tx_type, bd);
      break;
    case IDTX:
      highbd_inv_idtx_add_c(input, dest, stride, 4, tx_type, bd);
      break;
#endif  // CONFIG_EXT_TX
    default: assert(0); break;
  }
}

#if CONFIG_EXT_TX
void av1_highbd_inv_txfm_add_4x8(const tran_low_t *input, uint8_t *dest,
                                 int stride, int eob, int bd, TX_TYPE tx_type) {
  (void)eob;
  av1_highbd_iht4x8_32_add_c(input, dest, stride, tx_type, bd);
}

void av1_highbd_inv_txfm_add_8x4(const tran_low_t *input, uint8_t *dest,
                                 int stride, int eob, int bd, TX_TYPE tx_type) {
  (void)eob;
  av1_highbd_iht8x4_32_add_c(input, dest, stride, tx_type, bd);
}

void av1_highbd_inv_txfm_add_8x16(const tran_low_t *input, uint8_t *dest,
                                  int stride, int eob, int bd,
                                  TX_TYPE tx_type) {
  (void)eob;
  av1_highbd_iht8x16_128_add_c(input, dest, stride, tx_type, bd);
}

void av1_highbd_inv_txfm_add_16x8(const tran_low_t *input, uint8_t *dest,
                                  int stride, int eob, int bd,
                                  TX_TYPE tx_type) {
  (void)eob;
  av1_highbd_iht16x8_128_add_c(input, dest, stride, tx_type, bd);
}

void av1_highbd_inv_txfm_add_16x32(const tran_low_t *input, uint8_t *dest,
                                   int stride, int eob, int bd,
                                   TX_TYPE tx_type) {
  (void)eob;
  av1_highbd_iht16x32_512_add_c(input, dest, stride, tx_type, bd);
}

void av1_highbd_inv_txfm_add_32x16(const tran_low_t *input, uint8_t *dest,
                                   int stride, int eob, int bd,
                                   TX_TYPE tx_type) {
  (void)eob;
  av1_highbd_iht32x16_512_add_c(input, dest, stride, tx_type, bd);
}
#endif  // CONFIG_EXT_TX

void av1_highbd_inv_txfm_add_8x8(const tran_low_t *input, uint8_t *dest,
                                 int stride, int eob, int bd, TX_TYPE tx_type) {
  (void)eob;
  switch (tx_type) {
    case DCT_DCT:
    case ADST_DCT:
    case DCT_ADST:
    case ADST_ADST:
      av1_inv_txfm2d_add_8x8(input, CONVERT_TO_SHORTPTR(dest), stride, tx_type,
                             bd);
      break;
#if CONFIG_EXT_TX
    case FLIPADST_DCT:
    case DCT_FLIPADST:
    case FLIPADST_FLIPADST:
    case ADST_FLIPADST:
    case FLIPADST_ADST:
      av1_inv_txfm2d_add_8x8(input, CONVERT_TO_SHORTPTR(dest), stride, tx_type,
                             bd);
      break;
    case V_DCT:
    case H_DCT:
    case V_ADST:
    case H_ADST:
    case V_FLIPADST:
    case H_FLIPADST:
      // Use C version since DST only exists in C code
      av1_highbd_iht8x8_64_add_c(input, dest, stride, tx_type, bd);
      break;
    case IDTX:
      highbd_inv_idtx_add_c(input, dest, stride, 8, tx_type, bd);
      break;
#endif  // CONFIG_EXT_TX
    default: assert(0); break;
  }
}

void av1_highbd_inv_txfm_add_16x16(const tran_low_t *input, uint8_t *dest,
                                   int stride, int eob, int bd,
                                   TX_TYPE tx_type) {
  (void)eob;
  switch (tx_type) {
    case DCT_DCT:
    case ADST_DCT:
    case DCT_ADST:
    case ADST_ADST:
      av1_inv_txfm2d_add_16x16(input, CONVERT_TO_SHORTPTR(dest), stride,
                               tx_type, bd);
      break;
#if CONFIG_EXT_TX
    case FLIPADST_DCT:
    case DCT_FLIPADST:
    case FLIPADST_FLIPADST:
    case ADST_FLIPADST:
    case FLIPADST_ADST:
      av1_inv_txfm2d_add_16x16(input, CONVERT_TO_SHORTPTR(dest), stride,
                               tx_type, bd);
      break;
    case V_DCT:
    case H_DCT:
    case V_ADST:
    case H_ADST:
    case V_FLIPADST:
    case H_FLIPADST:
      // Use C version since DST only exists in C code
      av1_highbd_iht16x16_256_add_c(input, dest, stride, tx_type, bd);
      break;
    case IDTX:
      highbd_inv_idtx_add_c(input, dest, stride, 16, tx_type, bd);
      break;
#endif  // CONFIG_EXT_TX
    default: assert(0); break;
  }
}

void av1_highbd_inv_txfm_add_32x32(const tran_low_t *input, uint8_t *dest,
                                   int stride, int eob, int bd,
                                   TX_TYPE tx_type) {
  (void)eob;
  switch (tx_type) {
    case DCT_DCT:
      av1_inv_txfm2d_add_32x32(input, CONVERT_TO_SHORTPTR(dest), stride,
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
      av1_highbd_iht32x32_1024_add_c(input, dest, stride, tx_type, bd);
      break;
    case IDTX:
      highbd_inv_idtx_add_c(input, dest, stride, 32, tx_type, bd);
      break;
#endif  // CONFIG_EXT_TX
    default: assert(0); break;
  }
}
#endif  // CONFIG_AOM_HIGHBITDEPTH

void inv_txfm_add(const tran_low_t *input, uint8_t *dest, int stride,
                  INV_TXFM_PARAM *inv_txfm_param) {
  const TX_TYPE tx_type = inv_txfm_param->tx_type;
  const TX_SIZE tx_size = inv_txfm_param->tx_size;
  const int eob = inv_txfm_param->eob;
  const int lossless = inv_txfm_param->lossless;

  switch (tx_size) {
    case TX_32X32:
      av1_inv_txfm_add_32x32(input, dest, stride, eob, tx_type);
      break;
    case TX_16X16:
      av1_inv_txfm_add_16x16(input, dest, stride, eob, tx_type);
      break;
    case TX_8X8: av1_inv_txfm_add_8x8(input, dest, stride, eob, tx_type); break;
#if CONFIG_EXT_TX
    case TX_4X8: av1_inv_txfm_add_4x8(input, dest, stride, eob, tx_type); break;
    case TX_8X4: av1_inv_txfm_add_8x4(input, dest, stride, eob, tx_type); break;
    case TX_8X16:
      av1_inv_txfm_add_8x16(input, dest, stride, eob, tx_type);
      break;
    case TX_16X8:
      av1_inv_txfm_add_16x8(input, dest, stride, eob, tx_type);
      break;
    case TX_16X32:
      av1_inv_txfm_add_16x32(input, dest, stride, eob, tx_type);
      break;
    case TX_32X16:
      av1_inv_txfm_add_32x16(input, dest, stride, eob, tx_type);
      break;
#endif  // CONFIG_EXT_TX
    case TX_4X4:
      // this is like av1_short_idct4x4 but has a special case around eob<=1
      // which is significant (not just an optimization) for the lossless
      // case.
      av1_inv_txfm_add_4x4(input, dest, stride, eob, tx_type, lossless);
      break;
    default: assert(0 && "Invalid transform size"); break;
  }
}

#if CONFIG_AOM_HIGHBITDEPTH
void highbd_inv_txfm_add(const tran_low_t *input, uint8_t *dest, int stride,
                         INV_TXFM_PARAM *inv_txfm_param) {
  const TX_TYPE tx_type = inv_txfm_param->tx_type;
  const TX_SIZE tx_size = inv_txfm_param->tx_size;
  const int eob = inv_txfm_param->eob;
  const int bd = inv_txfm_param->bd;
  const int lossless = inv_txfm_param->lossless;

  switch (tx_size) {
    case TX_32X32:
      av1_highbd_inv_txfm_add_32x32(input, dest, stride, eob, bd, tx_type);
      break;
    case TX_16X16:
      av1_highbd_inv_txfm_add_16x16(input, dest, stride, eob, bd, tx_type);
      break;
    case TX_8X8:
      av1_highbd_inv_txfm_add_8x8(input, dest, stride, eob, bd, tx_type);
      break;
#if CONFIG_EXT_TX
    case TX_4X8:
      av1_highbd_inv_txfm_add_4x8(input, dest, stride, eob, bd, tx_type);
      break;
    case TX_8X4:
      av1_highbd_inv_txfm_add_8x4(input, dest, stride, eob, bd, tx_type);
      break;
    case TX_8X16:
      av1_highbd_inv_txfm_add_8x16(input, dest, stride, eob, bd, tx_type);
      break;
    case TX_16X8:
      av1_highbd_inv_txfm_add_16x8(input, dest, stride, eob, bd, tx_type);
      break;
    case TX_16X32:
      av1_highbd_inv_txfm_add_16x32(input, dest, stride, eob, bd, tx_type);
      break;
    case TX_32X16:
      av1_highbd_inv_txfm_add_32x16(input, dest, stride, eob, bd, tx_type);
      break;
#endif  // CONFIG_EXT_TX
    case TX_4X4:
      // this is like av1_short_idct4x4 but has a special case around eob<=1
      // which is significant (not just an optimization) for the lossless
      // case.
      av1_highbd_inv_txfm_add_4x4(input, dest, stride, eob, bd, tx_type,
                                  lossless);
      break;
    default: assert(0 && "Invalid transform size"); break;
  }
}
#endif  // CONFIG_AOM_HIGHBITDEPTH
