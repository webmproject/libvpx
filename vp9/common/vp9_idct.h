/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#ifndef VP9_COMMON_VP9_IDCT_H_
#define VP9_COMMON_VP9_IDCT_H_

#include <assert.h>

#include "./vpx_config.h"
#include "vp9/common/vp9_common.h"
#include "vp9/common/vp9_enums.h"

#ifdef __cplusplus
extern "C" {
#endif

// Constants and Macros used by all idct/dct functions
#define DCT_CONST_BITS 14
#define DCT_CONST_ROUNDING  (1 << (DCT_CONST_BITS - 1))

#define UNIT_QUANT_SHIFT 2
#define UNIT_QUANT_FACTOR (1 << UNIT_QUANT_SHIFT)

#define pair_set_epi16(a, b) \
  _mm_set_epi16(b, a, b, a, b, a, b, a)

#define dual_set_epi16(a, b) \
  _mm_set_epi16(b, b, b, b, a, a, a, a)

#if CONFIG_TX_SKIP
#define TX_SKIP_SHIFT_LQ 2
#define TX_SKIP_SHIFT_HQ 3
#endif

// Constants:
//  for (int i = 1; i< 32; ++i)
//    printf("static const int cospi_%d_64 = %.0f;\n", i,
//           round(16384 * cos(i*M_PI/64)));
// Note: sin(k*Pi/64) = cos((32-k)*Pi/64)
static const tran_high_t cospi_1_64  = 16364;
static const tran_high_t cospi_2_64  = 16305;
static const tran_high_t cospi_3_64  = 16207;
static const tran_high_t cospi_4_64  = 16069;
static const tran_high_t cospi_5_64  = 15893;
static const tran_high_t cospi_6_64  = 15679;
static const tran_high_t cospi_7_64  = 15426;
static const tran_high_t cospi_8_64  = 15137;
static const tran_high_t cospi_9_64  = 14811;
static const tran_high_t cospi_10_64 = 14449;
static const tran_high_t cospi_11_64 = 14053;
static const tran_high_t cospi_12_64 = 13623;
static const tran_high_t cospi_13_64 = 13160;
static const tran_high_t cospi_14_64 = 12665;
static const tran_high_t cospi_15_64 = 12140;
static const tran_high_t cospi_16_64 = 11585;
static const tran_high_t cospi_17_64 = 11003;
static const tran_high_t cospi_18_64 = 10394;
static const tran_high_t cospi_19_64 = 9760;
static const tran_high_t cospi_20_64 = 9102;
static const tran_high_t cospi_21_64 = 8423;
static const tran_high_t cospi_22_64 = 7723;
static const tran_high_t cospi_23_64 = 7005;
static const tran_high_t cospi_24_64 = 6270;
static const tran_high_t cospi_25_64 = 5520;
static const tran_high_t cospi_26_64 = 4756;
static const tran_high_t cospi_27_64 = 3981;
static const tran_high_t cospi_28_64 = 3196;
static const tran_high_t cospi_29_64 = 2404;
static const tran_high_t cospi_30_64 = 1606;
static const tran_high_t cospi_31_64 = 804;

//  16384 * sqrt(2) * sin(kPi/9) * 2 / 3
static const tran_high_t sinpi_1_9 = 5283;
static const tran_high_t sinpi_2_9 = 9929;
static const tran_high_t sinpi_3_9 = 13377;
static const tran_high_t sinpi_4_9 = 15212;

static INLINE tran_low_t check_range(tran_high_t input) {
#if CONFIG_VP9_HIGHBITDEPTH
  // For valid highbitdepth VP9 streams, intermediate stage coefficients will
  // stay within the ranges:
  // - 8 bit: signed 16 bit integer
  // - 10 bit: signed 18 bit integer
  // - 12 bit: signed 20 bit integer
#elif CONFIG_COEFFICIENT_RANGE_CHECKING
  // For valid VP9 input streams, intermediate stage coefficients should always
  // stay within the range of a signed 16 bit integer. Coefficients can go out
  // of this range for invalid/corrupt VP9 streams. However, strictly checking
  // this range for every intermediate coefficient can burdensome for a decoder,
  // therefore the following assertion is only enabled when configured with
  // --enable-coefficient-range-checking.
  assert(INT16_MIN <= input);
  assert(input <= INT16_MAX);
#endif
  return (tran_low_t)input;
}

static INLINE tran_low_t dct_const_round_shift(tran_high_t input) {
  tran_high_t rv = ROUND_POWER_OF_TWO(input, DCT_CONST_BITS);
  return check_range(rv);
}

typedef void (*transform_1d)(const tran_low_t*, tran_low_t*);

typedef struct {
  transform_1d cols, rows;  // vertical and horizontal
} transform_2d;

#if CONFIG_VP9_HIGHBITDEPTH
typedef void (*highbd_transform_1d)(const tran_low_t*, tran_low_t*, int bd);

typedef struct {
  highbd_transform_1d cols, rows;  // vertical and horizontal
} highbd_transform_2d;
#endif  // CONFIG_VP9_HIGHBITDEPTH

#if CONFIG_EMULATE_HARDWARE
// When CONFIG_EMULATE_HARDWARE is 1 the transform performs a
// non-normative method to handle overflows. A stream that causes
// overflows  in the inverse transform is considered invalid in VP9,
// and a hardware implementer is free to choose any reasonable
// method to handle overflows. However to aid in hardware
// verification they can use a specific implementation of the
// WRAPLOW() macro below that is identical to their intended
// hardware implementation (and also use configure options to trigger
// the C-implementation of the transform).
//
// The particular WRAPLOW implementation below performs strict
// overflow wrapping to match common hardware implementations.
// bd of 8 uses trans_low with 16bits, need to remove 16bits
// bd of 10 uses trans_low with 18bits, need to remove 14bits
// bd of 12 uses trans_low with 20bits, need to remove 12bits
// bd of x uses trans_low with 8+x bits, need to remove 24-x bits
#define WRAPLOW(x, bd) ((((int32_t)(x)) << (24 - bd)) >> (24 - bd))
#else
#define WRAPLOW(x, bd) (x)
#endif  // CONFIG_EMULATE_HARDWARE

static INLINE uint8_t clip_pixel_add(uint8_t dest, tran_high_t trans) {
  trans = WRAPLOW(trans, 8);
  return clip_pixel(WRAPLOW(dest + trans, 8));
}

#if CONFIG_EXT_TX
#define USE_DST2 1

#if USE_DST2
static const tran_high_t Tx4[4 * 4] = {
  // dst2
  6270,  15137,  15137,   6270,
  11585,  11585, -11585, -11585,
  15137,  -6270,  -6270,  15137,
  11585, -11585,  11585, -11585,
};

static const tran_high_t Tx8[8 * 8] = {
  // dst2
  3196,   9102,  13623,  16069,  16069,  13623,   9102,   3196,
  6270,  15137,  15137,   6270,  -6270, -15137, -15137,  -6270,
  9102,  16069,   3196, -13623, -13623,   3196,  16069,   9102,
  11585,  11585, -11585, -11585,  11585,  11585, -11585, -11585,
  13623,   3196, -16069,   9102,   9102, -16069,   3196,  13623,
  15137,  -6270,  -6270,  15137, -15137,   6270,   6270, -15137,
  16069, -13623,   9102,  -3196,  -3196,   9102, -13623,  16069,
  11585, -11585,  11585, -11585,  11585, -11585,  11585, -11585,
};

static const tran_high_t Tx16[16 * 16] = {
  // dst2
  1606,   4756,   7723,  10394,  12665,  14449,  15679,  16305,
  16305,  15679,  14449,  12665,  10394,   7723,   4756,   1606,
  3196,   9102,  13623,  16069,  16069,  13623,   9102,   3196,
  -3196,  -9102, -13623, -16069, -16069, -13623,  -9102,  -3196,
  4756,  12665,  16305,  14449,   7723,  -1606, -10394, -15679,
  -15679, -10394,  -1606,   7723,  14449,  16305,  12665,   4756,
  6270,  15137,  15137,   6270,  -6270, -15137, -15137,  -6270,
  6270,  15137,  15137,   6270,  -6270, -15137, -15137,  -6270,
  7723,  16305,  10394,  -4756, -15679, -12665,   1606,  14449,
  14449,   1606, -12665, -15679,  -4756,  10394,  16305,   7723,
  9102,  16069,   3196, -13623, -13623,   3196,  16069,   9102,
  -9102, -16069,  -3196,  13623,  13623,  -3196, -16069,  -9102,
  10394,  14449,  -4756, -16305,  -1606,  15679,   7723, -12665,
  -12665,   7723,  15679,  -1606, -16305,  -4756,  14449,  10394,
  11585,  11585, -11585, -11585,  11585,  11585, -11585, -11585,
  11585,  11585, -11585, -11585,  11585,  11585, -11585, -11585,
  12665,   7723, -15679,  -1606,  16305,  -4756, -14449,  10394,
  10394, -14449,  -4756,  16305,  -1606, -15679,   7723,  12665,
  13623,   3196, -16069,   9102,   9102, -16069,   3196,  13623,
  -13623,  -3196,  16069,  -9102,  -9102,  16069,  -3196, -13623,
  14449,  -1606, -12665,  15679,  -4756, -10394,  16305,  -7723,
  -7723,  16305, -10394,  -4756,  15679, -12665,  -1606,  14449,
  15137,  -6270,  -6270,  15137, -15137,   6270,   6270, -15137,
  15137,  -6270,  -6270,  15137, -15137,   6270,   6270, -15137,
  15679, -10394,   1606,   7723, -14449,  16305, -12665,   4756,
  4756, -12665,  16305, -14449,   7723,   1606, -10394,  15679,
  16069, -13623,   9102,  -3196,  -3196,   9102, -13623,  16069,
  -16069,  13623,  -9102,   3196,   3196,  -9102,  13623, -16069,
  16305, -15679,  14449, -12665,  10394,  -7723,   4756,  -1606,
  -1606,   4756,  -7723,  10394, -12665,  14449, -15679,  16305,
  11585, -11585,  11585, -11585,  11585, -11585,  11585, -11585,
  11585, -11585,  11585, -11585,  11585, -11585,  11585, -11585,
};
#endif  // USE_DST2

static INLINE void vp9_fgentx4(const tran_low_t *input, tran_low_t *output,
                               const tran_high_t *T) {
  tran_high_t sum;
  int i;
  const tran_high_t *Tx = T;
  for (i = 0; i < 4; ++i, Tx += 4) {
    sum = Tx[0] * input[0] + Tx[1] * input[1] +
          Tx[2] * input[2] + Tx[3] * input[3];
    output[i] = ROUND_POWER_OF_TWO(sum, DCT_CONST_BITS);
  }
}

static INLINE void vp9_fgentx8(const tran_low_t *input, tran_low_t *output,
                               const tran_high_t *T) {
  tran_high_t sum;
  int i;
  const tran_high_t *Tx = T;
  for (i = 0; i < 8; ++i, Tx += 8) {
    sum = Tx[0] * input[0] + Tx[1] * input[1] +
          Tx[2] * input[2] + Tx[3] * input[3] +
          Tx[4] * input[4] + Tx[5] * input[5] +
          Tx[6] * input[6] + Tx[7] * input[7];
    output[i] = ROUND_POWER_OF_TWO(sum, DCT_CONST_BITS);
  }
}

static INLINE void vp9_fgentx16(const tran_low_t *input, tran_low_t *output,
                                const tran_high_t *T) {
  tran_high_t sum;
  int i;
  const tran_high_t *Tx = T;
  for (i = 0; i < 16; ++i, Tx += 16) {
    sum = Tx[0] * input[0] + Tx[1] * input[1] +
          Tx[2] * input[2] + Tx[3] * input[3] +
          Tx[4] * input[4] + Tx[5] * input[5] +
          Tx[6] * input[6] + Tx[7] * input[7] +
          Tx[8] * input[8] + Tx[9] * input[9] +
          Tx[10] * input[10] + Tx[11] * input[11] +
          Tx[12] * input[12] + Tx[13] * input[13] +
          Tx[14] * input[14] + Tx[15] * input[15];
    output[i] = ROUND_POWER_OF_TWO(sum, DCT_CONST_BITS);
  }
}

static INLINE void vp9_igentx4(const tran_low_t *input, tran_low_t *output,
                               const tran_high_t *T) {
  tran_high_t sum[4];
  int i;
  const tran_high_t *Tx = T;
  for (i = 0; i < 4; ++i, ++Tx) {
    sum[i] = Tx[0] * input[0] + Tx[4] * input[1] +
             Tx[8] * input[2] + Tx[12] * input[3];
  }
  for (i = 0; i < 4; ++i) {
    output[i] = WRAPLOW(ROUND_POWER_OF_TWO(sum[i], DCT_CONST_BITS), 8);
  }
}

static INLINE void vp9_igentx8(const tran_low_t *input, tran_low_t *output,
                               const tran_high_t *T) {
  tran_high_t sum[8];
  int i;
  const tran_high_t *Tx = T;
  for (i = 0; i < 8; ++i, ++Tx) {
    sum[i] = Tx[0] * input[0] + Tx[8] * input[1] +
             Tx[16] * input[2] + Tx[24] * input[3] +
             Tx[32] * input[4] + Tx[40] * input[5] +
             Tx[48] * input[6] + Tx[56] * input[7];
  }
  for (i = 0; i < 8; ++i) {
    output[i] = WRAPLOW(ROUND_POWER_OF_TWO(sum[i], DCT_CONST_BITS), 8);
  }
}

static INLINE void vp9_igentx16(const tran_low_t *input, tran_low_t *output,
                                const tran_high_t *T) {
  tran_high_t sum[16];
  int i;
  const tran_high_t *Tx = T;
  for (i = 0; i < 16; ++i, ++Tx) {
    sum[i] = Tx[0] * input[0] + Tx[16] * input[1] +
             Tx[32] * input[2] + Tx[48] * input[3] +
             Tx[64] * input[4] + Tx[80] * input[5] +
             Tx[96] * input[6] + Tx[112] * input[7] +
             Tx[128] * input[8] + Tx[144] * input[9] +
             Tx[160] * input[10] + Tx[176] * input[11] +
             Tx[192] * input[12] + Tx[208] * input[13] +
             Tx[224] * input[14] + Tx[240] * input[15];
  }
  for (i = 0; i < 16; ++i) {
    output[i] = WRAPLOW(ROUND_POWER_OF_TWO(sum[i], DCT_CONST_BITS), 8);
  }
}

#if CONFIG_VP9_HIGHBITDEPTH
static INLINE void vp9_highbd_igentx4(const tran_low_t *input,
                                      tran_low_t *output,
                                      int bd, const tran_high_t *T) {
  tran_high_t sum[4];
  int i;
  const tran_high_t *Tx = T;
  (void) bd;
  for (i = 0; i < 4; ++i, Tx += 1) {
    sum[i] = Tx[0] * input[0] + Tx[4] * input[1] +
             Tx[8] * input[2] + Tx[12] * input[3];
  }
  for (i = 0; i < 4; ++i) {
    output[i] = WRAPLOW(ROUND_POWER_OF_TWO(sum[i], DCT_CONST_BITS), bd);
  }
}

static INLINE void vp9_highbd_igentx8(const tran_low_t *input,
                                      tran_low_t *output,
                                      int bd, const tran_high_t *T) {
  tran_high_t sum[8];
  int i;
  const tran_high_t *Tx = T;
  (void) bd;
  for (i = 0; i < 8; ++i, Tx += 1) {
    sum[i] = Tx[0] * input[0] + Tx[8] * input[1] +
             Tx[16] * input[2] + Tx[24] * input[3] +
             Tx[32] * input[4] + Tx[40] * input[5] +
             Tx[48] * input[6] + Tx[56] * input[7];
  }
  for (i = 0; i < 8; ++i) {
    output[i] = WRAPLOW(ROUND_POWER_OF_TWO(sum[i], DCT_CONST_BITS), bd);
  }
}

static INLINE void vp9_highbd_igentx16(const tran_low_t *input,
                                       tran_low_t *output,
                                       int bd, const tran_high_t *T) {
  tran_high_t sum[16];
  int i;
  const tran_high_t *Tx = T;
  (void) bd;
  for (i = 0; i < 16; ++i, Tx += 1) {
    sum[i] = Tx[0] * input[0] + Tx[16] * input[1] +
             Tx[32] * input[2] + Tx[48] * input[3] +
             Tx[64] * input[4] + Tx[80] * input[5] +
             Tx[96] * input[6] + Tx[112] * input[7] +
             Tx[128] * input[8] + Tx[144] * input[9] +
             Tx[160] * input[10] + Tx[176] * input[11] +
             Tx[192] * input[12] + Tx[208] * input[13] +
             Tx[224] * input[14] + Tx[240] * input[15];
  }
  for (i = 0; i < 16; ++i) {
    output[i] = WRAPLOW(ROUND_POWER_OF_TWO(sum[i], DCT_CONST_BITS), bd);
  }
}
#endif  // CONFIG_VP9_HIGHBITDEPTH
#endif  // CONFIG_EXT_TX

void vp9_iwht4x4_add(const tran_low_t *input, uint8_t *dest, int stride,
                     int eob);
void vp9_idct4x4_add(const tran_low_t *input, uint8_t *dest, int stride,
                     int eob);
void vp9_idct8x8_add(const tran_low_t *input, uint8_t *dest, int stride,
                     int eob);
void vp9_idct16x16_add(const tran_low_t *input, uint8_t *dest, int stride,
                       int eob);
void vp9_idct32x32_add(const tran_low_t *input, uint8_t *dest, int stride,
                       int eob);
#if CONFIG_TX64X64
void vp9_idct64x64_add(const tran_low_t *input, uint8_t *dest, int stride,
                       int eob);
#endif  // CONFIG_TX64X64
void vp9_iht4x4_add(TX_TYPE tx_type, const tran_low_t *input, uint8_t *dest,
                    int stride, int eob);
void vp9_iht8x8_add(TX_TYPE tx_type, const tran_low_t *input, uint8_t *dest,
                    int stride, int eob);
void vp9_iht16x16_add(TX_TYPE tx_type, const tran_low_t *input, uint8_t *dest,
                      int stride, int eob);

#if CONFIG_VP9_HIGHBITDEPTH
void vp9_highbd_idct4(const tran_low_t *input, tran_low_t *output, int bd);
void vp9_highbd_idct8(const tran_low_t *input, tran_low_t *output, int bd);
void vp9_highbd_idct16(const tran_low_t *input, tran_low_t *output, int bd);
void vp9_highbd_iwht4x4_add(const tran_low_t *input, uint8_t *dest, int stride,
                            int eob, int bd);
void vp9_highbd_idct4x4_add(const tran_low_t *input, uint8_t *dest, int stride,
                            int eob, int bd);
void vp9_highbd_idct8x8_add(const tran_low_t *input, uint8_t *dest, int stride,
                            int eob, int bd);
void vp9_highbd_idct16x16_add(const tran_low_t *input, uint8_t *dest,
                              int stride, int eob, int bd);
void vp9_highbd_idct32x32_add(const tran_low_t *input, uint8_t *dest,
                              int stride, int eob, int bd);
#if CONFIG_TX64X64
void vp9_highbd_idct64x64_add(const tran_low_t *input, uint8_t *dest,
                              int stride, int eob, int bd);
#endif
void vp9_highbd_iht4x4_add(TX_TYPE tx_type, const tran_low_t *input,
                           uint8_t *dest, int stride, int eob, int bd);
void vp9_highbd_iht8x8_add(TX_TYPE tx_type, const tran_low_t *input,
                           uint8_t *dest, int stride, int eob, int bd);
void vp9_highbd_iht16x16_add(TX_TYPE tx_type, const tran_low_t *input,
                             uint8_t *dest, int stride, int eob, int bd);
static INLINE uint16_t highbd_clip_pixel_add(uint16_t dest, tran_high_t trans,
                                             int bd) {
  trans = WRAPLOW(trans, bd);
  return clip_pixel_highbd(WRAPLOW(dest + trans, bd), bd);
}
#endif  // CONFIG_VP9_HIGHBITDEPTH

#if CONFIG_TX_SKIP
void vp9_tx_identity_add_rect(const tran_low_t *input, uint8_t *dest,
                              int row, int col, int stride_in,
                              int stride_out, int shift);
void vp9_tx_identity_add(const tran_low_t *input, uint8_t *dest,
                         int stride, int bs, int shift);
#if CONFIG_VP9_HIGHBITDEPTH
void vp9_highbd_tx_identity_add_rect(const tran_low_t *input, uint8_t *dest,
                                     int row, int col, int stride_in,
                                     int stride_out, int shift, int bd);
void vp9_highbd_tx_identity_add(const tran_low_t *input, uint8_t *dest,
                                int stride, int bs, int shift, int bd);
#endif  // CONFIG_VP9_HIGHBITDEPTH
#endif  // CONFIG_TX_SKIP

void vp9_dst1d_type1(int64_t *in, int64_t *out, int N);
void vp9_idst4x4_add(const tran_low_t *input, uint8_t *dest, int stride);
void vp9_idst8x8_add(const tran_low_t *input, uint8_t *dest, int stride);
void vp9_idst16x16_add(const tran_low_t *input, uint8_t *dest, int stride);

#ifdef __cplusplus
}  // extern "C"
#endif

#endif  // VP9_COMMON_VP9_IDCT_H_
#if CONFIG_SR_MODE
void vp9_iwht4x4(const tran_low_t *input, int16_t *dest, int stride,
                 int eob);
void vp9_idct4x4(const tran_low_t *input, int16_t *dest, int stride,
                 int eob);
void vp9_idct8x8(const tran_low_t *input, int16_t *dest, int stride,
                 int eob);
void vp9_idct16x16(const tran_low_t *input, int16_t *dest, int stride,
                   int eob);
void vp9_idct32x32(const tran_low_t *input, int16_t *dest, int stride,
                   int eob);
void vp9_iht4x4(TX_TYPE tx_type, const tran_low_t *input, int16_t *dest,
                int stride, int eob);
void vp9_iht8x8(TX_TYPE tx_type, const tran_low_t *input, int16_t *dest,
                int stride, int eob);
void vp9_iht16x16(TX_TYPE tx_type, const tran_low_t *input, int16_t *dest,
                  int stride, int eob);
#if CONFIG_TX64X64
void vp9_idct64x64(const tran_low_t *input, int16_t *dest, int stride,
                   int eob);
#endif  // CONFIG_TX64X64
#endif  // CONFIG_SR_MODE
