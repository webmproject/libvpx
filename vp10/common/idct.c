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
#include "vp10/common/idct.h"
#include "vpx_dsp/inv_txfm.h"
#include "vpx_ports/mem.h"

#if CONFIG_EXT_TX
void idst4_c(const tran_low_t *input, tran_low_t *output) {
  static const int N = 4;
  static const int sinvalue_lookup_table[] = {
    9630, 15582
  };
  static const int mult = 14654;  // sqrt(4/5)
  int i, j;
  for (i = 0; i < N; i++) {
    int64_t sum = 0;
    for (j = 0; j < N; j++) {
      int idx = (i + 1) * (j + 1);
      int sign = 0;
      if (idx > N + 1) {
        sign = (idx / (N + 1)) & 1;
        idx %= (N + 1);
      }
      idx = idx > N + 1 - idx ? N + 1 - idx : idx;
      if (idx == 0) continue;
      idx--;
      sum += (int64_t)input[j] * sinvalue_lookup_table[idx] * (sign ? -1 : 1);
    }
    sum = (sum * mult) >> (2 * DCT_CONST_BITS);
    output[i] = WRAPLOW(sum, 8);
  }
}

void idst8_c(const tran_low_t *input, tran_low_t *output) {
  static const int N = 8;
  static const int sinvalue_lookup_table[] = {
    5604, 10531, 14189, 16135
  };
  static const int mult = 15447;  // 2*sqrt(2/9)
  int i, j;
  for (i = 0; i < N; i++) {
    int64_t sum = 0;
    for (j = 0; j < N; j++) {
      int idx = (i + 1) * (j + 1);
      int sign = 0;
      if (idx > N + 1) {
        sign = (idx / (N + 1)) & 1;
        idx %= (N + 1);
      }
      idx = idx > N + 1 - idx ? N + 1 - idx : idx;
      if (idx == 0) continue;
      idx--;
      sum += (int64_t)input[j] * sinvalue_lookup_table[idx] * (sign ? -1 : 1);
    }
    sum = (sum * mult) >> (2 * DCT_CONST_BITS);
    output[i] = WRAPLOW(sum, 8);
  }
}

void idst16_c(const tran_low_t *input, tran_low_t *output) {
  static const int N = 16;
  static const int sinvalue_lookup_table[] = {
    3011,  5919,  8625, 11038,
    13075, 14666, 15759, 16314
  };
  static const int mult = 15895;  // 2*sqrt(4/17)
  int i, j;
  for (i = 0; i < N; i++) {
    int64_t sum = 0;
    for (j = 0; j < N; j++) {
      int idx = (i + 1) * (j + 1);
      int sign = 0;
      if (idx > N + 1) {
        sign = (idx / (N + 1)) & 1;
        idx %= (N + 1);
      }
      idx = idx > N + 1 - idx ? N + 1 - idx : idx;
      if (idx == 0) continue;
      idx--;
      sum += (int64_t)input[j] * sinvalue_lookup_table[idx] * (sign ? -1 : 1);
    }
    sum = (sum * mult) >> (2 * DCT_CONST_BITS);
    output[i] = WRAPLOW(sum, 8);
  }
}

#if CONFIG_VP9_HIGHBITDEPTH
void highbd_idst4_c(const tran_low_t *input, tran_low_t *output, int bd) {
  static const int N = 4;
  static const int sinvalue_lookup_table[] = {
    9630, 15582
  };
  static const int mult = 14654;  // sqrt(4/5)
  int i, j;
  (void) bd;
  for (i = 0; i < N; i++) {
    int64_t sum = 0;
    for (j = 0; j < N; j++) {
      int idx = (i + 1) * (j + 1);
      int sign = 0;
      if (idx > N + 1) {
        sign = (idx / (N + 1)) & 1;
        idx %= (N + 1);
      }
      idx = idx > N + 1 - idx ? N + 1 - idx : idx;
      if (idx == 0) continue;
      idx--;
      sum += (int64_t)input[j] * sinvalue_lookup_table[idx] * (sign ? -1 : 1);
    }
    sum = (sum * mult) >> (2 * DCT_CONST_BITS);
    output[i] = WRAPLOW(sum, bd);
  }
}

void highbd_idst8_c(const tran_low_t *input, tran_low_t *output, int bd) {
  static const int N = 8;
  static const int sinvalue_lookup_table[] = {
    5604, 10531, 14189, 16135
  };
  static const int mult = 15447;  // 2*sqrt(2/9)
  int i, j;
  (void) bd;
  for (i = 0; i < N; i++) {
    int64_t sum = 0;
    for (j = 0; j < N; j++) {
      int idx = (i + 1) * (j + 1);
      int sign = 0;
      if (idx > N + 1) {
        sign = (idx / (N + 1)) & 1;
        idx %= (N + 1);
      }
      idx = idx > N + 1 - idx ? N + 1 - idx : idx;
      if (idx == 0) continue;
      idx--;
      sum += (int64_t)input[j] * sinvalue_lookup_table[idx] * (sign ? -1 : 1);
    }
    sum = (sum * mult) >> (2 * DCT_CONST_BITS);
    output[i] = WRAPLOW(sum, bd);
  }
}

void highbd_idst16_c(const tran_low_t *input, tran_low_t *output, int bd) {
  static const int N = 16;
  static const int sinvalue_lookup_table[] = {
    3011,  5919,  8625, 11038,
    13075, 14666, 15759, 16314
  };
  static const int mult = 15895;  // 2*sqrt(4/17)
  int i, j;
  (void) bd;
  for (i = 0; i < N; i++) {
    int64_t sum = 0;
    for (j = 0; j < N; j++) {
      int idx = (i + 1) * (j + 1);
      int sign = 0;
      if (idx > N + 1) {
        sign = (idx / (N + 1)) & 1;
        idx %= (N + 1);
      }
      idx = idx > N + 1 - idx ? N + 1 - idx : idx;
      if (idx == 0) continue;
      idx--;
      sum += (int64_t)input[j] * sinvalue_lookup_table[idx] * (sign ? -1 : 1);
    }
    sum = (sum * mult) >> (2 * DCT_CONST_BITS);
    output[i] = WRAPLOW(sum, bd);
  }
}
#endif  // CONFIG_VP9_HIGHBITDEPTH
#endif  // CONFIG_EXT_TX

#if CONFIG_EXT_TX
void fliplr(uint8_t *dest, int stride, int l) {
  int i, j;
  for (i = 0; i < l; ++i) {
    for (j = 0; j < l / 2; ++j) {
      const uint8_t tmp = dest[i * stride + j];
      dest[i * stride + j] = dest[i * stride + l - 1 - j];
      dest[i * stride + l - 1 - j] = tmp;
    }
  }
}

void flipud(uint8_t *dest, int stride, int l) {
  int i, j;
  for (j = 0; j < l; ++j) {
    for (i = 0; i < l / 2; ++i) {
      const uint8_t tmp = dest[i * stride + j];
      dest[i * stride + j] = dest[(l - 1 - i) * stride + j];
      dest[(l - 1 - i) * stride + j] = tmp;
    }
  }
}

void fliplrud(uint8_t *dest, int stride, int l) {
  int i, j;
  for (i = 0; i < l / 2; ++i) {
    for (j = 0; j < l; ++j) {
      const uint8_t tmp = dest[i * stride + j];
      dest[i * stride + j] = dest[(l - 1 - i) * stride + l - 1 - j];
      dest[(l - 1 - i) * stride + l - 1 - j] = tmp;
    }
  }
}

void fliplr16(uint16_t *dest, int stride, int l) {
  int i, j;
  for (i = 0; i < l; ++i) {
    for (j = 0; j < l / 2; ++j) {
      const uint16_t tmp = dest[i * stride + j];
      dest[i * stride + j] = dest[i * stride + l - 1 - j];
      dest[i * stride + l - 1 - j] = tmp;
    }
  }
}

void flipud16(uint16_t *dest, int stride, int l) {
  int i, j;
  for (j = 0; j < l; ++j) {
    for (i = 0; i < l / 2; ++i) {
      const uint16_t tmp = dest[i * stride + j];
      dest[i * stride + j] = dest[(l - 1 - i) * stride + j];
      dest[(l - 1 - i) * stride + j] = tmp;
    }
  }
}

void fliplrud16(uint16_t *dest, int stride, int l) {
  int i, j;
  for (i = 0; i < l / 2; ++i) {
    for (j = 0; j < l; ++j) {
      const uint16_t tmp = dest[i * stride + j];
      dest[i * stride + j] = dest[(l - 1 - i) * stride + l - 1 - j];
      dest[(l - 1 - i) * stride + l - 1 - j] = tmp;
    }
  }
}
#endif  // CONFIG_EXT_TX

void vp10_iht4x4_16_add_c(const tran_low_t *input, uint8_t *dest, int stride,
                         int tx_type) {
  const transform_2d IHT_4[] = {
    { idct4_c, idct4_c  },   // DCT_DCT  = 0
    { iadst4_c, idct4_c  },  // ADST_DCT = 1
    { idct4_c, iadst4_c },   // DCT_ADST = 2
    { iadst4_c, iadst4_c },  // ADST_ADST = 3
#if CONFIG_EXT_TX
    { iadst4_c, idct4_c },   // FLIPADST_DCT = 4
    { idct4_c,  iadst4_c },  // DCT_FLIPADST = 5
    { iadst4_c, iadst4_c },  // FLIPADST_FLIPADST = 6
    { iadst4_c, iadst4_c },  // ADST_FLIPADST = 7
    { iadst4_c, iadst4_c },  // FLIPADST_ADST = 8
    { idst4_c,  idst4_c },   // DST_DST = 9
    { idst4_c,  idct4_c  },  // DST_DCT = 10
    { idct4_c,  idst4_c  },  // DCT_DST = 11
    { idst4_c,  iadst4_c },  // DST_ADST = 12
    { iadst4_c, idst4_c  },  // ADST_DST = 13
    { idst4_c,  iadst4_c },  // DST_FLIPADST = 14
    { iadst4_c, idst4_c  },  // FLIPADST_DST = 15
#endif  // CONFIG_EXT_TX
  };

  int i, j;
  tran_low_t out[4 * 4];
  tran_low_t *outptr = out;
  tran_low_t temp_in[4], temp_out[4];

  // inverse transform row vectors
  for (i = 0; i < 4; ++i) {
    IHT_4[tx_type].rows(input, outptr);
    input  += 4;
    outptr += 4;
  }

  // inverse transform column vectors
  for (i = 0; i < 4; ++i) {
    for (j = 0; j < 4; ++j)
      temp_in[j] = out[j * 4 + i];
    IHT_4[tx_type].cols(temp_in, temp_out);
    for (j = 0; j < 4; ++j) {
      dest[j * stride + i] = clip_pixel_add(dest[j * stride + i],
                                            ROUND_POWER_OF_TWO(temp_out[j], 4));
    }
  }
}

void vp10_iht8x8_64_add_c(const tran_low_t *input, uint8_t *dest, int stride,
                         int tx_type) {
  static const transform_2d IHT_8[] = {
    { idct8_c,  idct8_c  },  // DCT_DCT  = 0
    { iadst8_c, idct8_c  },  // ADST_DCT = 1
    { idct8_c,  iadst8_c },  // DCT_ADST = 2
    { iadst8_c, iadst8_c },  // ADST_ADST = 3
#if CONFIG_EXT_TX
    { iadst8_c, idct8_c },   // FLIPADST_DCT = 4
    { idct8_c,  iadst8_c },  // DCT_FLIPADST = 5
    { iadst8_c, iadst8_c },  // FLIPADST_FLIPADST = 6
    { iadst8_c, iadst8_c },  // ADST_FLIPADST = 7
    { iadst8_c, iadst8_c },  // FLIPADST_ADST = 8
    { idst8_c,  idst8_c },   // DST_DST = 9
    { idst8_c,  idct8_c  },  // DST_DCT = 10
    { idct8_c,  idst8_c  },  // DCT_DST = 11
    { idst8_c,  iadst8_c },  // DST_ADST = 12
    { iadst8_c, idst8_c  },  // ADST_DST = 13
    { idst8_c,  iadst8_c },  // DST_FLIPADST = 14
    { iadst8_c, idst8_c  },  // FLIPADST_DST = 15
#endif  // CONFIG_EXT_TX
  };

  int i, j;
  tran_low_t out[8 * 8];
  tran_low_t *outptr = out;
  tran_low_t temp_in[8], temp_out[8];
  const transform_2d ht = IHT_8[tx_type];

  // inverse transform row vectors
  for (i = 0; i < 8; ++i) {
    ht.rows(input, outptr);
    input += 8;
    outptr += 8;
  }

  // inverse transform column vectors
  for (i = 0; i < 8; ++i) {
    for (j = 0; j < 8; ++j)
      temp_in[j] = out[j * 8 + i];
    ht.cols(temp_in, temp_out);
    for (j = 0; j < 8; ++j) {
      dest[j * stride + i] = clip_pixel_add(dest[j * stride + i],
                                            ROUND_POWER_OF_TWO(temp_out[j], 5));
    }
  }
}

void vp10_iht16x16_256_add_c(const tran_low_t *input, uint8_t *dest, int stride,
                            int tx_type) {
  static const transform_2d IHT_16[] = {
    { idct16_c,  idct16_c  },  // DCT_DCT  = 0
    { iadst16_c, idct16_c  },  // ADST_DCT = 1
    { idct16_c,  iadst16_c },  // DCT_ADST = 2
    { iadst16_c, iadst16_c },  // ADST_ADST = 3
#if CONFIG_EXT_TX
    { iadst16_c, idct16_c  },  // FLIPADST_DCT = 4
    { idct16_c,  iadst16_c },  // DCT_FLIPADST = 5
    { iadst16_c, iadst16_c },  // FLIPADST_FLIPADST = 6
    { iadst16_c, iadst16_c },  // ADST_FLIPADST = 7
    { iadst16_c, iadst16_c },  // FLIPADST_ADST = 8
    { idst16_c,  idst16_c  },  // DST_DST = 9
    { idst16_c,  idct16_c  },  // DST_DCT = 10
    { idct16_c,  idst16_c  },  // DCT_DST = 11
    { idst16_c,  iadst16_c },  // DST_ADST = 12
    { iadst16_c, idst16_c  },  // ADST_DST = 13
    { idst16_c,  iadst16_c },  // DST_FLIPADST = 14
    { iadst16_c, idst16_c  },  // FLIPADST_DST = 15
#endif  // CONFIG_EXT_TX
  };

  int i, j;
  tran_low_t out[16 * 16];
  tran_low_t *outptr = out;
  tran_low_t temp_in[16], temp_out[16];
  const transform_2d ht = IHT_16[tx_type];

  // Rows
  for (i = 0; i < 16; ++i) {
    ht.rows(input, outptr);
    input += 16;
    outptr += 16;
  }

  // Columns
  for (i = 0; i < 16; ++i) {
    for (j = 0; j < 16; ++j)
      temp_in[j] = out[j * 16 + i];
    ht.cols(temp_in, temp_out);
    for (j = 0; j < 16; ++j) {
      dest[j * stride + i] = clip_pixel_add(dest[j * stride + i],
                                            ROUND_POWER_OF_TWO(temp_out[j], 6));
    }
  }
}

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

void vp10_inv_txfm_add_4x4(
    const tran_low_t *input, uint8_t *dest,
    int stride, int eob, TX_TYPE tx_type,
    void (*itxm_add_4x4)(const tran_low_t *input,
                         uint8_t *dest, int stride, int eob)) {
  switch (tx_type) {
    case DCT_DCT:
      itxm_add_4x4(input, dest, stride, eob);
      break;
    case ADST_DCT:
    case DCT_ADST:
    case ADST_ADST:
      vp10_iht4x4_16_add(input, dest, stride, tx_type);
      break;
#if CONFIG_EXT_TX
    case FLIPADST_DCT:
      flipud(dest, stride, 4);
      vp10_iht4x4_16_add(input, dest, stride, ADST_DCT);
      flipud(dest, stride, 4);
      break;
    case DCT_FLIPADST:
      fliplr(dest, stride, 4);
      vp10_iht4x4_16_add(input, dest, stride, DCT_ADST);
      fliplr(dest, stride, 4);
      break;
    case FLIPADST_FLIPADST:
      fliplrud(dest, stride, 4);
      vp10_iht4x4_16_add(input, dest, stride, ADST_ADST);
      fliplrud(dest, stride, 4);
      break;
    case ADST_FLIPADST:
      fliplr(dest, stride, 4);
      vp10_iht4x4_16_add(input, dest, stride, ADST_ADST);
      fliplr(dest, stride, 4);
      break;
    case FLIPADST_ADST:
      flipud(dest, stride, 4);
      vp10_iht4x4_16_add(input, dest, stride, ADST_ADST);
      flipud(dest, stride, 4);
      break;
    case DST_DST:
    case DST_DCT:
    case DCT_DST:
    case DST_ADST:
    case ADST_DST:
      // Use C version since DST only exists in C code
      vp10_iht4x4_16_add_c(input, dest, stride, tx_type);
      break;
    case FLIPADST_DST:
      flipud(dest, stride, 4);
      vp10_iht4x4_16_add_c(input, dest, stride, ADST_DST);
      flipud(dest, stride, 4);
      break;
    case DST_FLIPADST:
      fliplr(dest, stride, 4);
      vp10_iht4x4_16_add_c(input, dest, stride, DST_ADST);
      fliplr(dest, stride, 4);
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
      flipud(dest, stride, 8);
      vp10_iht8x8_64_add(input, dest, stride, ADST_DCT);
      flipud(dest, stride, 8);
      break;
    case DCT_FLIPADST:
      fliplr(dest, stride, 8);
      vp10_iht8x8_64_add(input, dest, stride, DCT_ADST);
      fliplr(dest, stride, 8);
      break;
    case FLIPADST_FLIPADST:
      fliplrud(dest, stride, 8);
      vp10_iht8x8_64_add(input, dest, stride, ADST_ADST);
      fliplrud(dest, stride, 8);
      break;
    case ADST_FLIPADST:
      fliplr(dest, stride, 8);
      vp10_iht8x8_64_add(input, dest, stride, ADST_ADST);
      fliplr(dest, stride, 8);
      break;
    case FLIPADST_ADST:
      flipud(dest, stride, 8);
      vp10_iht8x8_64_add(input, dest, stride, ADST_ADST);
      flipud(dest, stride, 8);
      break;
    case DST_DST:
    case DST_DCT:
    case DCT_DST:
    case DST_ADST:
    case ADST_DST:
      // Use C version since DST only exists in C code
      vp10_iht8x8_64_add_c(input, dest, stride, tx_type);
      break;
    case FLIPADST_DST:
      flipud(dest, stride, 8);
      vp10_iht8x8_64_add_c(input, dest, stride, ADST_DST);
      flipud(dest, stride, 8);
      break;
    case DST_FLIPADST:
      fliplr(dest, stride, 8);
      vp10_iht8x8_64_add_c(input, dest, stride, DST_ADST);
      fliplr(dest, stride, 8);
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
      flipud(dest, stride, 16);
      vp10_iht16x16_256_add(input, dest, stride, ADST_DCT);
      flipud(dest, stride, 16);
      break;
    case DCT_FLIPADST:
      fliplr(dest, stride, 16);
      vp10_iht16x16_256_add(input, dest, stride, DCT_ADST);
      fliplr(dest, stride, 16);
      break;
    case FLIPADST_FLIPADST:
      fliplrud(dest, stride, 16);
      vp10_iht16x16_256_add(input, dest, stride, ADST_ADST);
      fliplrud(dest, stride, 16);
      break;
    case ADST_FLIPADST:
      fliplr(dest, stride, 16);
      vp10_iht16x16_256_add(input, dest, stride, ADST_ADST);
      fliplr(dest, stride, 16);
      break;
    case FLIPADST_ADST:
      flipud(dest, stride, 16);
      vp10_iht16x16_256_add(input, dest, stride, ADST_ADST);
      flipud(dest, stride, 16);
      break;
    case DST_DST:
    case DST_DCT:
    case DCT_DST:
    case DST_ADST:
    case ADST_DST:
      // Use C version since DST only exists in C code
      vp10_iht16x16_256_add_c(input, dest, stride, tx_type);
      break;
    case FLIPADST_DST:
      flipud(dest, stride, 16);
      vp10_iht16x16_256_add_c(input, dest, stride, ADST_DST);
      flipud(dest, stride, 16);
      break;
    case DST_FLIPADST:
      fliplr(dest, stride, 16);
      vp10_iht16x16_256_add_c(input, dest, stride, DST_ADST);
      fliplr(dest, stride, 16);
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
    case ADST_DCT:
    case DCT_ADST:
    case ADST_ADST:
      assert(0);
      break;
    default:
      assert(0);
      break;
  }
}

#if CONFIG_VP9_HIGHBITDEPTH
void vp10_highbd_iht4x4_16_add_c(const tran_low_t *input, uint8_t *dest8,
                                int stride, int tx_type, int bd) {
  const highbd_transform_2d IHT_4[] = {
    { vpx_highbd_idct4_c,  vpx_highbd_idct4_c  },  // DCT_DCT  = 0
    { vpx_highbd_iadst4_c, vpx_highbd_idct4_c  },  // ADST_DCT = 1
    { vpx_highbd_idct4_c,  vpx_highbd_iadst4_c },  // DCT_ADST = 2
    { vpx_highbd_iadst4_c, vpx_highbd_iadst4_c },  // ADST_ADST = 3
#if CONFIG_EXT_TX
    { vpx_highbd_iadst4_c, vpx_highbd_idct4_c  },  // FLIPADST_DCT = 4
    { vpx_highbd_idct4_c,  vpx_highbd_iadst4_c },  // DCT_FLIPADST = 5
    { vpx_highbd_iadst4_c, vpx_highbd_iadst4_c },  // FLIPADST_FLIPADST = 6
    { vpx_highbd_iadst4_c, vpx_highbd_iadst4_c },  // ADST_FLIPADST = 7
    { vpx_highbd_iadst4_c, vpx_highbd_iadst4_c },  // FLIPADST_ADST = 8
    { highbd_idst4_c,      highbd_idst4_c      },  // DST_DST = 9
    { highbd_idst4_c,      vpx_highbd_idct4_c  },  // DST_DCT = 10
    { vpx_highbd_idct4_c,  highbd_idst4_c      },  // DCT_DST = 11
    { highbd_idst4_c,      vpx_highbd_iadst4_c },  // DST_ADST = 12
    { vpx_highbd_iadst4_c, highbd_idst4_c      },  // ADST_DST = 13
    { highbd_idst4_c,      vpx_highbd_iadst4_c },  // DST_FLIPADST = 14
    { vpx_highbd_iadst4_c, highbd_idst4_c      },  // FLIPADST_DST = 15
#endif  // CONFIG_EXT_TX
  };
  uint16_t *dest = CONVERT_TO_SHORTPTR(dest8);

  int i, j;
  tran_low_t out[4 * 4];
  tran_low_t *outptr = out;
  tran_low_t temp_in[4], temp_out[4];

  // Inverse transform row vectors.
  for (i = 0; i < 4; ++i) {
    IHT_4[tx_type].rows(input, outptr, bd);
    input  += 4;
    outptr += 4;
  }

  // Inverse transform column vectors.
  for (i = 0; i < 4; ++i) {
    for (j = 0; j < 4; ++j)
      temp_in[j] = out[j * 4 + i];
    IHT_4[tx_type].cols(temp_in, temp_out, bd);
    for (j = 0; j < 4; ++j) {
      dest[j * stride + i] = highbd_clip_pixel_add(
          dest[j * stride + i], ROUND_POWER_OF_TWO(temp_out[j], 4), bd);
    }
  }
}

void vp10_highbd_iht8x8_64_add_c(const tran_low_t *input, uint8_t *dest8,
                                int stride, int tx_type, int bd) {
  static const highbd_transform_2d HIGH_IHT_8[] = {
    { vpx_highbd_idct8_c,  vpx_highbd_idct8_c  },  // DCT_DCT  = 0
    { vpx_highbd_iadst8_c, vpx_highbd_idct8_c  },  // ADST_DCT = 1
    { vpx_highbd_idct8_c,  vpx_highbd_iadst8_c },  // DCT_ADST = 2
    { vpx_highbd_iadst8_c, vpx_highbd_iadst8_c },  // ADST_ADST = 3
#if CONFIG_EXT_TX
    { vpx_highbd_iadst8_c, vpx_highbd_idct8_c  },  // FLIPADST_DCT = 4
    { vpx_highbd_idct8_c,  vpx_highbd_iadst8_c },  // DCT_FLIPADST = 5
    { vpx_highbd_iadst8_c, vpx_highbd_iadst8_c },  // FLIPADST_FLIPADST = 6
    { vpx_highbd_iadst8_c, vpx_highbd_iadst8_c },  // ADST_FLIPADST = 7
    { vpx_highbd_iadst8_c, vpx_highbd_iadst8_c },  // FLIPADST_ADST = 8
    { highbd_idst8_c,      highbd_idst8_c      },  // DST_DST = 9
    { highbd_idst8_c,      vpx_highbd_idct8_c  },  // DST_DCT = 10
    { vpx_highbd_idct8_c,  highbd_idst8_c      },  // DCT_DST = 11
    { highbd_idst8_c,      vpx_highbd_iadst8_c },  // DST_ADST = 12
    { vpx_highbd_iadst8_c, highbd_idst8_c      },  // ADST_DST = 13
    { highbd_idst8_c,      vpx_highbd_iadst8_c },  // DST_FLIPADST = 14
    { vpx_highbd_iadst8_c, highbd_idst8_c      },  // FLIPADST_DST = 15
#endif  // CONFIG_EXT_TX
  };

  int i, j;
  tran_low_t out[8 * 8];
  tran_low_t *outptr = out;
  tran_low_t temp_in[8], temp_out[8];
  const highbd_transform_2d ht = HIGH_IHT_8[tx_type];
  uint16_t *dest = CONVERT_TO_SHORTPTR(dest8);

  // Inverse transform row vectors.
  for (i = 0; i < 8; ++i) {
    ht.rows(input, outptr, bd);
    input += 8;
    outptr += 8;
  }

  // Inverse transform column vectors.
  for (i = 0; i < 8; ++i) {
    for (j = 0; j < 8; ++j)
      temp_in[j] = out[j * 8 + i];
    ht.cols(temp_in, temp_out, bd);
    for (j = 0; j < 8; ++j) {
      dest[j * stride + i] = highbd_clip_pixel_add(
          dest[j * stride + i], ROUND_POWER_OF_TWO(temp_out[j], 5), bd);
    }
  }
}

void vp10_highbd_iht16x16_256_add_c(const tran_low_t *input, uint8_t *dest8,
                                   int stride, int tx_type, int bd) {
  static const highbd_transform_2d HIGH_IHT_16[] = {
    { vpx_highbd_idct16_c,  vpx_highbd_idct16_c  },  // DCT_DCT  = 0
    { vpx_highbd_iadst16_c, vpx_highbd_idct16_c  },  // ADST_DCT = 1
    { vpx_highbd_idct16_c,  vpx_highbd_iadst16_c },  // DCT_ADST = 2
    { vpx_highbd_iadst16_c, vpx_highbd_iadst16_c },  // ADST_ADST = 3
#if CONFIG_EXT_TX
    { vpx_highbd_iadst16_c, vpx_highbd_idct16_c  },  // FLIPADST_DCT = 4
    { vpx_highbd_idct16_c,  vpx_highbd_iadst16_c },  // DCT_FLIPADST = 5
    { vpx_highbd_iadst16_c, vpx_highbd_iadst16_c },  // FLIPADST_FLIPADST = 6
    { vpx_highbd_iadst16_c, vpx_highbd_iadst16_c },  // ADST_FLIPADST = 7
    { vpx_highbd_iadst16_c, vpx_highbd_iadst16_c },  // FLIPADST_ADST = 8
    { highbd_idst16_c,      highbd_idst16_c      },  // DST_DST = 9
    { highbd_idst16_c,      vpx_highbd_idct16_c  },  // DST_DCT = 10
    { vpx_highbd_idct16_c,  highbd_idst16_c      },  // DCT_DST = 11
    { highbd_idst16_c,      vpx_highbd_iadst16_c },  // DST_ADST = 12
    { vpx_highbd_iadst16_c, highbd_idst16_c      },  // ADST_DST = 13
    { highbd_idst16_c,      vpx_highbd_iadst16_c },  // DST_FLIPADST = 14
    { vpx_highbd_iadst16_c, highbd_idst16_c      },  // FLIPADST_DST = 15
#endif  // CONFIG_EXT_TX
  };

  int i, j;
  tran_low_t out[16 * 16];
  tran_low_t *outptr = out;
  tran_low_t temp_in[16], temp_out[16];
  const highbd_transform_2d ht = HIGH_IHT_16[tx_type];
  uint16_t *dest = CONVERT_TO_SHORTPTR(dest8);

  // Rows
  for (i = 0; i < 16; ++i) {
    ht.rows(input, outptr, bd);
    input += 16;
    outptr += 16;
  }

  // Columns
  for (i = 0; i < 16; ++i) {
    for (j = 0; j < 16; ++j)
      temp_in[j] = out[j * 16 + i];
    ht.cols(temp_in, temp_out, bd);
    for (j = 0; j < 16; ++j) {
      dest[j * stride + i] = highbd_clip_pixel_add(
          dest[j * stride + i], ROUND_POWER_OF_TWO(temp_out[j], 6), bd);
    }
  }
}

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
                                  void (*highbd_itxm_add_4x4)
                                  (const tran_low_t *input, uint8_t *dest,
                                      int stride, int eob, int bd)) {
  switch (tx_type) {
    case DCT_DCT:
      highbd_itxm_add_4x4(input, dest, stride, eob, bd);
      break;
    case ADST_DCT:
    case DCT_ADST:
    case ADST_ADST:
      vp10_highbd_iht4x4_16_add(input, dest, stride, tx_type, bd);
      break;
#if CONFIG_EXT_TX
    case FLIPADST_DCT:
      flipud16(CONVERT_TO_SHORTPTR(dest), stride, 4);
      vp10_highbd_iht4x4_16_add(input, dest, stride, ADST_DCT, bd);
      flipud16(CONVERT_TO_SHORTPTR(dest), stride, 4);
      break;
    case DCT_FLIPADST:
      fliplr16(CONVERT_TO_SHORTPTR(dest), stride, 4);
      vp10_highbd_iht4x4_16_add(input, dest, stride, DCT_ADST, bd);
      fliplr16(CONVERT_TO_SHORTPTR(dest), stride, 4);
      break;
    case FLIPADST_FLIPADST:
      fliplrud16(CONVERT_TO_SHORTPTR(dest), stride, 4);
      vp10_highbd_iht4x4_16_add(input, dest, stride, ADST_ADST, bd);
      fliplrud16(CONVERT_TO_SHORTPTR(dest), stride, 4);
      break;
    case ADST_FLIPADST:
      fliplr16(CONVERT_TO_SHORTPTR(dest), stride, 4);
      vp10_highbd_iht4x4_16_add(input, dest, stride, ADST_ADST, bd);
      fliplr16(CONVERT_TO_SHORTPTR(dest), stride, 4);
      break;
    case FLIPADST_ADST:
      flipud16(CONVERT_TO_SHORTPTR(dest), stride, 4);
      vp10_highbd_iht4x4_16_add(input, dest, stride, ADST_ADST, bd);
      flipud16(CONVERT_TO_SHORTPTR(dest), stride, 4);
      break;
    case DST_DST:
    case DST_DCT:
    case DCT_DST:
    case DST_ADST:
    case ADST_DST:
      // Use C version since DST only exists in C code
      vp10_highbd_iht4x4_16_add_c(input, dest, stride, tx_type, bd);
      break;
    case FLIPADST_DST:
      flipud16(CONVERT_TO_SHORTPTR(dest), stride, 4);
      vp10_highbd_iht4x4_16_add_c(input, dest, stride, ADST_DST, bd);
      flipud16(CONVERT_TO_SHORTPTR(dest), stride, 4);
      break;
    case DST_FLIPADST:
      fliplr16(CONVERT_TO_SHORTPTR(dest), stride, 4);
      vp10_highbd_iht4x4_16_add_c(input, dest, stride, DST_ADST, bd);
      fliplr16(CONVERT_TO_SHORTPTR(dest), stride, 4);
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
  switch (tx_type) {
    case DCT_DCT:
      vp10_highbd_idct8x8_add(input, dest, stride, eob, bd);
      break;
    case ADST_DCT:
    case DCT_ADST:
    case ADST_ADST:
      vp10_highbd_iht8x8_64_add(input, dest, stride, tx_type, bd);
      break;
#if CONFIG_EXT_TX
    case FLIPADST_DCT:
      flipud16(CONVERT_TO_SHORTPTR(dest), stride, 8);
      vp10_highbd_iht8x8_64_add(input, dest, stride, ADST_DCT, bd);
      flipud16(CONVERT_TO_SHORTPTR(dest), stride, 8);
      break;
    case DCT_FLIPADST:
      fliplr16(CONVERT_TO_SHORTPTR(dest), stride, 8);
      vp10_highbd_iht8x8_64_add(input, dest, stride, DCT_ADST, bd);
      fliplr16(CONVERT_TO_SHORTPTR(dest), stride, 8);
      break;
    case FLIPADST_FLIPADST:
      fliplrud16(CONVERT_TO_SHORTPTR(dest), stride, 8);
      vp10_highbd_iht8x8_64_add(input, dest, stride, ADST_ADST, bd);
      fliplrud16(CONVERT_TO_SHORTPTR(dest), stride, 8);
      break;
    case ADST_FLIPADST:
      fliplr16(CONVERT_TO_SHORTPTR(dest), stride, 8);
      vp10_highbd_iht8x8_64_add(input, dest, stride, ADST_ADST, bd);
      fliplr16(CONVERT_TO_SHORTPTR(dest), stride, 8);
      break;
    case FLIPADST_ADST:
      flipud16(CONVERT_TO_SHORTPTR(dest), stride, 8);
      vp10_highbd_iht8x8_64_add(input, dest, stride, ADST_ADST, bd);
      flipud16(CONVERT_TO_SHORTPTR(dest), stride, 8);
      break;
    case DST_DST:
    case DST_DCT:
    case DCT_DST:
    case DST_ADST:
    case ADST_DST:
      // Use C version since DST only exists in C code
      vp10_highbd_iht8x8_64_add_c(input, dest, stride, tx_type, bd);
      break;
    case FLIPADST_DST:
      flipud16(CONVERT_TO_SHORTPTR(dest), stride, 8);
      vp10_highbd_iht8x8_64_add_c(input, dest, stride, ADST_DST, bd);
      flipud16(CONVERT_TO_SHORTPTR(dest), stride, 8);
      break;
    case DST_FLIPADST:
      fliplr16(CONVERT_TO_SHORTPTR(dest), stride, 8);
      vp10_highbd_iht8x8_64_add_c(input, dest, stride, DST_ADST, bd);
      fliplr16(CONVERT_TO_SHORTPTR(dest), stride, 8);
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
  switch (tx_type) {
    case DCT_DCT:
      vp10_highbd_idct16x16_add(input, dest, stride, eob, bd);
      break;
    case ADST_DCT:
    case DCT_ADST:
    case ADST_ADST:
      vp10_highbd_iht16x16_256_add(input, dest, stride, tx_type, bd);
      break;
#if CONFIG_EXT_TX
    case FLIPADST_DCT:
      flipud16(CONVERT_TO_SHORTPTR(dest), stride, 16);
      vp10_highbd_iht16x16_256_add(input, dest, stride, ADST_DCT, bd);
      flipud16(CONVERT_TO_SHORTPTR(dest), stride, 16);
      break;
    case DCT_FLIPADST:
      fliplr16(CONVERT_TO_SHORTPTR(dest), stride, 16);
      vp10_highbd_iht16x16_256_add(input, dest, stride, DCT_ADST, bd);
      fliplr16(CONVERT_TO_SHORTPTR(dest), stride, 16);
      break;
    case FLIPADST_FLIPADST:
      fliplrud16(CONVERT_TO_SHORTPTR(dest), stride, 16);
      vp10_highbd_iht16x16_256_add(input, dest, stride, ADST_ADST, bd);
      fliplrud16(CONVERT_TO_SHORTPTR(dest), stride, 16);
      break;
    case ADST_FLIPADST:
      fliplr16(CONVERT_TO_SHORTPTR(dest), stride, 16);
      vp10_highbd_iht16x16_256_add(input, dest, stride, ADST_ADST, bd);
      fliplr16(CONVERT_TO_SHORTPTR(dest), stride, 16);
      break;
    case FLIPADST_ADST:
      flipud16(CONVERT_TO_SHORTPTR(dest), stride, 16);
      vp10_highbd_iht16x16_256_add(input, dest, stride, ADST_ADST, bd);
      flipud16(CONVERT_TO_SHORTPTR(dest), stride, 16);
      break;
    case DST_DST:
    case DST_DCT:
    case DCT_DST:
    case DST_ADST:
    case ADST_DST:
      // Use C version since DST only exists in C code
      vp10_highbd_iht16x16_256_add_c(input, dest, stride, tx_type, bd);
      break;
    case FLIPADST_DST:
      flipud16(CONVERT_TO_SHORTPTR(dest), stride, 16);
      vp10_highbd_iht16x16_256_add_c(input, dest, stride, ADST_DST, bd);
      flipud16(CONVERT_TO_SHORTPTR(dest), stride, 16);
      break;
    case DST_FLIPADST:
      fliplr16(CONVERT_TO_SHORTPTR(dest), stride, 16);
      vp10_highbd_iht16x16_256_add_c(input, dest, stride, DST_ADST, bd);
      fliplr16(CONVERT_TO_SHORTPTR(dest), stride, 16);
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
  switch (tx_type) {
    case DCT_DCT:
      vp10_highbd_idct32x32_add(input, dest, stride, eob, bd);
      break;
    case ADST_DCT:
    case DCT_ADST:
    case ADST_ADST:
      assert(0);
      break;
    default:
      assert(0);
      break;
  }
}
#endif  // CONFIG_VP9_HIGHBITDEPTH
