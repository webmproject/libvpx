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

#include "./vp9_rtcd.h"
#include "vp9/common/vp9_systemdependent.h"
#include "vp9/common/vp9_blockd.h"
#include "vp9/common/vp9_idct.h"

#if CONFIG_DST1
// Integers to represent double.
// The sine transform formula is: X{i} = Sum_{0<=j<N}( x_j *
// sin((i+1)*(j+1)/(N+1) * PI) ) * sqrt(2/(N+1))
// e.g. when N == 4, it is a series of sin(PI*i/5) and sqrt(2/5). Similar for
// N = 8 and 16.
// For integer calculation, we multiply 2^14.
// sin_pi_5 = sin(PI/5)*pow(2,14)
// sqrt_2_5 = sqrt(2/5)*pow(2,14)

// {sin(pi/5), sin(pi*2/5)}
int sinvalue_lookup_table_4[2] = { 9630, 15582 };
// {sin(pi/9), sin(pi*2/9), ..., sin(pi*4/9)}
int sinvalue_lookup_table_8[4] = { 5604, 10531, 14189, 16135 };
// {sin(pi/17), ...
int sinvalue_lookup_table_16[] = { 3011,  5919,  8625, 11038,
                                  13075, 14666, 15759, 16314 };

void vp9_dst1d_type1(int64_t *in, int64_t *out, int N) {
  int i, j;
  for (i = 0; i < N; i++) {
    int64_t sum = 0;
    for (j = 0; j < N; j++) {
      int64_t sinvalue = 0;
      int idx = (i + 1) * (j + 1);
      int sign = 0;
      if (idx > N + 1) {
        sign = idx / (N + 1);
        sign = sign % 2 ? 1 : 0;
        idx %= (N + 1);
      }
      idx = idx > N + 1 - idx ? N + 1 - idx : idx;
      if (idx == 0) continue;
      idx--;

      if (N == 4)
        sinvalue = sinvalue_lookup_table_4[idx];
      else if (N == 8)
        sinvalue = sinvalue_lookup_table_8[idx];
      else if (N == 16)
        sinvalue = sinvalue_lookup_table_16[idx];
      else
        assert(0 && "Invalid transform size.");
      if (sign) sinvalue = -sinvalue;

      sum += in[j] * sinvalue;
    }
    out[i] = sum;
  }
}

static void idstNxN_add(const tran_low_t *input, uint8_t *output,
                        int stride, int N) {
  const int val_2_5 = 6554;
  const int val_2_9 = 3641;
  const int val_2_17 = 1928;
  int i, j;
  int64_t *in = (int64_t *) malloc (N * sizeof(int64_t));
  int64_t *inter = (int64_t *) malloc (N * sizeof(int64_t));
  int64_t *mat = (int64_t *) malloc (N * N * sizeof(int64_t));
  int64_t *mat2 = (int64_t *) malloc (N * N * sizeof(int64_t));
  int64_t val;

  // 1d dst: transform columns
  for (j = 0; j < N; j++) {
    for (i = 0; i < N; i++) {
      in[i] = input[i * N + j];
    }
    vp9_dst1d_type1(in, inter, N);
    for (i = 0; i < N; i++) {
      mat2[i * N + j] = inter[i];
    }
  }

  // transpose
  for (i = 0; i < N; i++)
    for (j = 0; j < N; j++)
      mat[i * N + j] = mat2[i + j * N];

  switch (N) {
    case 4:
      val = val_2_5;
      break;
    case 8:
      val = val_2_9;
      break;
    case 16:
      val = val_2_17;
      break;
    default:
      assert(0 && "Invalid transform size.");
      return;
  }

  // 1d dst: transform rows
  for (j = 0; j < N; j++) {
    for (i = 0; i < N; i++) {
      in[i] = mat[i * N + j];
    }
    vp9_dst1d_type1(in, inter, N);
    for (i = 0; i < N; i++) {
      int64_t tmp = inter[i];
      tmp = tmp >> DCT_CONST_BITS;
      tmp *= val;
      mat[i*N + j] = tmp >> (2 * DCT_CONST_BITS);
    }
  }
  for (i = 0; i < N; i++) {
    for (j = 0; j < N; j++) {
      tran_high_t tmp = mat[i * N + j];
      tmp = WRAPLOW(tmp, 8);
      output[i * stride + j] = clip_pixel_add(output[i * stride + j],
                                              ROUND_POWER_OF_TWO(tmp, 3));
    }
  }
  free(in);
  free(inter);
  free(mat);
  free(mat2);
}

void vp9_idst4x4_add(const tran_low_t *input, uint8_t *dest, int stride) {
  idstNxN_add(input, dest, stride, 4);
}

void vp9_idst8x8_add(const tran_low_t *input, uint8_t *dest, int stride) {
  idstNxN_add(input, dest, stride, 8);
}

void vp9_idst16x16_add(const tran_low_t *input, uint8_t *dest, int stride) {
  idstNxN_add(input, dest, stride, 16);
}
#endif  // CONFIG_DST1

void vp9_iwht4x4_16_add_c(const tran_low_t *input, uint8_t *dest, int stride) {
  /* 4-point reversible, orthonormal inverse Walsh-Hadamard in 3.5 adds,
     0.5 shifts per pixel. */
  int i;
    tran_low_t output[16];
    tran_high_t a1, b1, c1, d1, e1;
    const tran_low_t *ip = input;
    tran_low_t *op = output;

    for (i = 0; i < 4; i++) {
      a1 = ip[0] >> UNIT_QUANT_SHIFT;
      c1 = ip[1] >> UNIT_QUANT_SHIFT;
      d1 = ip[2] >> UNIT_QUANT_SHIFT;
      b1 = ip[3] >> UNIT_QUANT_SHIFT;
      a1 += c1;
      d1 -= b1;
      e1 = (a1 - d1) >> 1;
      b1 = e1 - b1;
      c1 = e1 - c1;
      a1 -= b1;
      d1 += c1;
      op[0] = WRAPLOW(a1, 8);
      op[1] = WRAPLOW(b1, 8);
      op[2] = WRAPLOW(c1, 8);
      op[3] = WRAPLOW(d1, 8);
      ip += 4;
      op += 4;
    }

    ip = output;
    for (i = 0; i < 4; i++) {
      a1 = ip[4 * 0];
      c1 = ip[4 * 1];
    d1 = ip[4 * 2];
    b1 = ip[4 * 3];
    a1 += c1;
    d1 -= b1;
    e1 = (a1 - d1) >> 1;
    b1 = e1 - b1;
    c1 = e1 - c1;
    a1 -= b1;
    d1 += c1;
    dest[stride * 0] = clip_pixel_add(dest[stride * 0], a1);
    dest[stride * 1] = clip_pixel_add(dest[stride * 1], b1);
    dest[stride * 2] = clip_pixel_add(dest[stride * 2], c1);
    dest[stride * 3] = clip_pixel_add(dest[stride * 3], d1);

    ip++;
    dest++;
  }
}

void vp9_iwht4x4_1_add_c(const tran_low_t *in, uint8_t *dest, int dest_stride) {
  int i;
  tran_high_t a1, e1;
  tran_low_t tmp[4];
  const tran_low_t *ip = in;
  tran_low_t *op = tmp;

  a1 = ip[0] >> UNIT_QUANT_SHIFT;
  e1 = a1 >> 1;
  a1 -= e1;
  op[0] = WRAPLOW(a1, 8);
  op[1] = op[2] = op[3] = WRAPLOW(e1, 8);

  ip = tmp;
  for (i = 0; i < 4; i++) {
    e1 = ip[0] >> 1;
    a1 = ip[0] - e1;
    dest[dest_stride * 0] = clip_pixel_add(dest[dest_stride * 0], a1);
    dest[dest_stride * 1] = clip_pixel_add(dest[dest_stride * 1], e1);
    dest[dest_stride * 2] = clip_pixel_add(dest[dest_stride * 2], e1);
    dest[dest_stride * 3] = clip_pixel_add(dest[dest_stride * 3], e1);
    ip++;
    dest++;
  }
}

static void idct4(const tran_low_t *input, tran_low_t *output) {
  tran_low_t step[4];
  tran_high_t temp1, temp2;
  // stage 1
  temp1 = (input[0] + input[2]) * cospi_16_64;
  temp2 = (input[0] - input[2]) * cospi_16_64;
  step[0] = WRAPLOW(dct_const_round_shift(temp1), 8);
  step[1] = WRAPLOW(dct_const_round_shift(temp2), 8);
  temp1 = input[1] * cospi_24_64 - input[3] * cospi_8_64;
  temp2 = input[1] * cospi_8_64 + input[3] * cospi_24_64;
  step[2] = WRAPLOW(dct_const_round_shift(temp1), 8);
  step[3] = WRAPLOW(dct_const_round_shift(temp2), 8);

  // stage 2
  output[0] = WRAPLOW(step[0] + step[3], 8);
  output[1] = WRAPLOW(step[1] + step[2], 8);
  output[2] = WRAPLOW(step[1] - step[2], 8);
  output[3] = WRAPLOW(step[0] - step[3], 8);
}

void vp9_idct4x4_16_add_c(const tran_low_t *input, uint8_t *dest, int stride) {
  tran_low_t out[4 * 4];
  tran_low_t *outptr = out;
  int i, j;
  tran_low_t temp_in[4], temp_out[4];

  // Rows
  for (i = 0; i < 4; ++i) {
    idct4(input, outptr);
    input += 4;
    outptr += 4;
  }

  // Columns
  for (i = 0; i < 4; ++i) {
    for (j = 0; j < 4; ++j)
      temp_in[j] = out[j * 4 + i];
    idct4(temp_in, temp_out);
    for (j = 0; j < 4; ++j) {
      dest[j * stride + i] = clip_pixel_add(dest[j * stride + i],
                                            ROUND_POWER_OF_TWO(temp_out[j], 4));
    }
  }
}

void vp9_idct4x4_1_add_c(const tran_low_t *input, uint8_t *dest,
                         int dest_stride) {
  int i;
  tran_high_t a1;
  tran_low_t out = WRAPLOW(dct_const_round_shift(input[0] * cospi_16_64), 8);
  out = WRAPLOW(dct_const_round_shift(out * cospi_16_64), 8);
  a1 = ROUND_POWER_OF_TWO(out, 4);

  for (i = 0; i < 4; i++) {
    dest[0] = clip_pixel_add(dest[0], a1);
    dest[1] = clip_pixel_add(dest[1], a1);
    dest[2] = clip_pixel_add(dest[2], a1);
    dest[3] = clip_pixel_add(dest[3], a1);
    dest += dest_stride;
  }
}

static void idct8(const tran_low_t *input, tran_low_t *output) {
  tran_low_t step1[8], step2[8];
  tran_high_t temp1, temp2;
  // stage 1
  step1[0] = input[0];
  step1[2] = input[4];
  step1[1] = input[2];
  step1[3] = input[6];
  temp1 = input[1] * cospi_28_64 - input[7] * cospi_4_64;
  temp2 = input[1] * cospi_4_64 + input[7] * cospi_28_64;
  step1[4] = WRAPLOW(dct_const_round_shift(temp1), 8);
  step1[7] = WRAPLOW(dct_const_round_shift(temp2), 8);
  temp1 = input[5] * cospi_12_64 - input[3] * cospi_20_64;
  temp2 = input[5] * cospi_20_64 + input[3] * cospi_12_64;
  step1[5] = WRAPLOW(dct_const_round_shift(temp1), 8);
  step1[6] = WRAPLOW(dct_const_round_shift(temp2), 8);

  // stage 2 & stage 3 - even half
  idct4(step1, step1);

  // stage 2 - odd half
  step2[4] = WRAPLOW(step1[4] + step1[5], 8);
  step2[5] = WRAPLOW(step1[4] - step1[5], 8);
  step2[6] = WRAPLOW(-step1[6] + step1[7], 8);
  step2[7] = WRAPLOW(step1[6] + step1[7], 8);

  // stage 3 -odd half
  step1[4] = step2[4];
  temp1 = (step2[6] - step2[5]) * cospi_16_64;
  temp2 = (step2[5] + step2[6]) * cospi_16_64;
  step1[5] = WRAPLOW(dct_const_round_shift(temp1), 8);
  step1[6] = WRAPLOW(dct_const_round_shift(temp2), 8);
  step1[7] = step2[7];

  // stage 4
  output[0] = WRAPLOW(step1[0] + step1[7], 8);
  output[1] = WRAPLOW(step1[1] + step1[6], 8);
  output[2] = WRAPLOW(step1[2] + step1[5], 8);
  output[3] = WRAPLOW(step1[3] + step1[4], 8);
  output[4] = WRAPLOW(step1[3] - step1[4], 8);
  output[5] = WRAPLOW(step1[2] - step1[5], 8);
  output[6] = WRAPLOW(step1[1] - step1[6], 8);
  output[7] = WRAPLOW(step1[0] - step1[7], 8);
}

void vp9_idct8x8_64_add_c(const tran_low_t *input, uint8_t *dest, int stride) {
  tran_low_t out[8 * 8];
  tran_low_t *outptr = out;
  int i, j;
  tran_low_t temp_in[8], temp_out[8];

  // First transform rows
  for (i = 0; i < 8; ++i) {
    idct8(input, outptr);
    input += 8;
    outptr += 8;
  }

  // Then transform columns
  for (i = 0; i < 8; ++i) {
    for (j = 0; j < 8; ++j)
      temp_in[j] = out[j * 8 + i];
    idct8(temp_in, temp_out);
    for (j = 0; j < 8; ++j) {
      dest[j * stride + i] = clip_pixel_add(dest[j * stride + i],
                                            ROUND_POWER_OF_TWO(temp_out[j], 5));
    }
  }
}

void vp9_idct8x8_1_add_c(const tran_low_t *input, uint8_t *dest, int stride) {
  int i, j;
  tran_high_t a1;
  tran_low_t out = WRAPLOW(dct_const_round_shift(input[0] * cospi_16_64), 8);
  out = WRAPLOW(dct_const_round_shift(out * cospi_16_64), 8);
  a1 = ROUND_POWER_OF_TWO(out, 5);
  for (j = 0; j < 8; ++j) {
    for (i = 0; i < 8; ++i)
      dest[i] = clip_pixel_add(dest[i], a1);
    dest += stride;
  }
}

static void iadst4(const tran_low_t *input, tran_low_t *output) {
  tran_high_t s0, s1, s2, s3, s4, s5, s6, s7;

  tran_low_t x0 = input[0];
  tran_low_t x1 = input[1];
  tran_low_t x2 = input[2];
  tran_low_t x3 = input[3];

  if (!(x0 | x1 | x2 | x3)) {
    output[0] = output[1] = output[2] = output[3] = 0;
    return;
  }

  s0 = sinpi_1_9 * x0;
  s1 = sinpi_2_9 * x0;
  s2 = sinpi_3_9 * x1;
  s3 = sinpi_4_9 * x2;
  s4 = sinpi_1_9 * x2;
  s5 = sinpi_2_9 * x3;
  s6 = sinpi_4_9 * x3;
  s7 = x0 - x2 + x3;

  s0 = s0 + s3 + s5;
  s1 = s1 - s4 - s6;
  s3 = s2;
  s2 = sinpi_3_9 * s7;

  // 1-D transform scaling factor is sqrt(2).
  // The overall dynamic range is 14b (input) + 14b (multiplication scaling)
  // + 1b (addition) = 29b.
  // Hence the output bit depth is 15b.
  output[0] = WRAPLOW(dct_const_round_shift(s0 + s3), 8);
  output[1] = WRAPLOW(dct_const_round_shift(s1 + s3), 8);
  output[2] = WRAPLOW(dct_const_round_shift(s2), 8);
  output[3] = WRAPLOW(dct_const_round_shift(s0 + s1 - s3), 8);
}

void vp9_iht4x4_16_add_c(const tran_low_t *input, uint8_t *dest, int stride,
                         int tx_type) {
  const transform_2d IHT_4[] = {
    { idct4, idct4  },  // DCT_DCT  = 0
    { iadst4, idct4  },   // ADST_DCT = 1
    { idct4, iadst4 },    // DCT_ADST = 2
    { iadst4, iadst4 }      // ADST_ADST = 3
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

static void iadst8(const tran_low_t *input, tran_low_t *output) {
  int s0, s1, s2, s3, s4, s5, s6, s7;

  tran_high_t x0 = input[7];
  tran_high_t x1 = input[0];
  tran_high_t x2 = input[5];
  tran_high_t x3 = input[2];
  tran_high_t x4 = input[3];
  tran_high_t x5 = input[4];
  tran_high_t x6 = input[1];
  tran_high_t x7 = input[6];

  if (!(x0 | x1 | x2 | x3 | x4 | x5 | x6 | x7)) {
    output[0] = output[1] = output[2] = output[3] = output[4]
        = output[5] = output[6] = output[7] = 0;
    return;
  }

  // stage 1
  s0 = cospi_2_64  * x0 + cospi_30_64 * x1;
  s1 = cospi_30_64 * x0 - cospi_2_64  * x1;
  s2 = cospi_10_64 * x2 + cospi_22_64 * x3;
  s3 = cospi_22_64 * x2 - cospi_10_64 * x3;
  s4 = cospi_18_64 * x4 + cospi_14_64 * x5;
  s5 = cospi_14_64 * x4 - cospi_18_64 * x5;
  s6 = cospi_26_64 * x6 + cospi_6_64  * x7;
  s7 = cospi_6_64  * x6 - cospi_26_64 * x7;

  x0 = WRAPLOW(dct_const_round_shift(s0 + s4), 8);
  x1 = WRAPLOW(dct_const_round_shift(s1 + s5), 8);
  x2 = WRAPLOW(dct_const_round_shift(s2 + s6), 8);
  x3 = WRAPLOW(dct_const_round_shift(s3 + s7), 8);
  x4 = WRAPLOW(dct_const_round_shift(s0 - s4), 8);
  x5 = WRAPLOW(dct_const_round_shift(s1 - s5), 8);
  x6 = WRAPLOW(dct_const_round_shift(s2 - s6), 8);
  x7 = WRAPLOW(dct_const_round_shift(s3 - s7), 8);

  // stage 2
  s0 = x0;
  s1 = x1;
  s2 = x2;
  s3 = x3;
  s4 =  cospi_8_64  * x4 + cospi_24_64 * x5;
  s5 =  cospi_24_64 * x4 - cospi_8_64  * x5;
  s6 = -cospi_24_64 * x6 + cospi_8_64  * x7;
  s7 =  cospi_8_64  * x6 + cospi_24_64 * x7;

  x0 = WRAPLOW(s0 + s2, 8);
  x1 = WRAPLOW(s1 + s3, 8);
  x2 = WRAPLOW(s0 - s2, 8);
  x3 = WRAPLOW(s1 - s3, 8);
  x4 = WRAPLOW(dct_const_round_shift(s4 + s6), 8);
  x5 = WRAPLOW(dct_const_round_shift(s5 + s7), 8);
  x6 = WRAPLOW(dct_const_round_shift(s4 - s6), 8);
  x7 = WRAPLOW(dct_const_round_shift(s5 - s7), 8);

  // stage 3
  s2 = cospi_16_64 * (x2 + x3);
  s3 = cospi_16_64 * (x2 - x3);
  s6 = cospi_16_64 * (x6 + x7);
  s7 = cospi_16_64 * (x6 - x7);

  x2 = WRAPLOW(dct_const_round_shift(s2), 8);
  x3 = WRAPLOW(dct_const_round_shift(s3), 8);
  x6 = WRAPLOW(dct_const_round_shift(s6), 8);
  x7 = WRAPLOW(dct_const_round_shift(s7), 8);

  output[0] = WRAPLOW(x0, 8);
  output[1] = WRAPLOW(-x4, 8);
  output[2] = WRAPLOW(x6, 8);
  output[3] = WRAPLOW(-x2, 8);
  output[4] = WRAPLOW(x3, 8);
  output[5] = WRAPLOW(-x7, 8);
  output[6] = WRAPLOW(x5, 8);
  output[7] = WRAPLOW(-x1, 8);
}

static const transform_2d IHT_8[] = {
  { idct8,  idct8  },  // DCT_DCT  = 0
  { iadst8, idct8  },  // ADST_DCT = 1
  { idct8,  iadst8 },  // DCT_ADST = 2
  { iadst8, iadst8 }   // ADST_ADST = 3
};

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

void vp9_iht8x8_64_add_c(const tran_low_t *input, uint8_t *dest, int stride,
                         int tx_type) {
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

void vp9_idct8x8_12_add_c(const tran_low_t *input, uint8_t *dest, int stride) {
  tran_low_t out[8 * 8] = { 0 };
  tran_low_t *outptr = out;
  int i, j;
  tran_low_t temp_in[8], temp_out[8];

  // First transform rows
  // only first 4 row has non-zero coefs
  for (i = 0; i < 4; ++i) {
    idct8(input, outptr);
    input += 8;
    outptr += 8;
  }

  // Then transform columns
  for (i = 0; i < 8; ++i) {
    for (j = 0; j < 8; ++j)
      temp_in[j] = out[j * 8 + i];
    idct8(temp_in, temp_out);
    for (j = 0; j < 8; ++j) {
      dest[j * stride + i] = clip_pixel_add(dest[j * stride + i],
                                            ROUND_POWER_OF_TWO(temp_out[j], 5));
    }
  }
}

static void idct16(const tran_low_t *input, tran_low_t *output) {
  tran_low_t step1[16], step2[16];
  tran_high_t temp1, temp2;

  // stage 1
  step1[0] = input[0/2];
  step1[1] = input[16/2];
  step1[2] = input[8/2];
  step1[3] = input[24/2];
  step1[4] = input[4/2];
  step1[5] = input[20/2];
  step1[6] = input[12/2];
  step1[7] = input[28/2];
  step1[8] = input[2/2];
  step1[9] = input[18/2];
  step1[10] = input[10/2];
  step1[11] = input[26/2];
  step1[12] = input[6/2];
  step1[13] = input[22/2];
  step1[14] = input[14/2];
  step1[15] = input[30/2];

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
  step2[8] = WRAPLOW(dct_const_round_shift(temp1), 8);
  step2[15] = WRAPLOW(dct_const_round_shift(temp2), 8);

  temp1 = step1[9] * cospi_14_64 - step1[14] * cospi_18_64;
  temp2 = step1[9] * cospi_18_64 + step1[14] * cospi_14_64;
  step2[9] = WRAPLOW(dct_const_round_shift(temp1), 8);
  step2[14] = WRAPLOW(dct_const_round_shift(temp2), 8);

  temp1 = step1[10] * cospi_22_64 - step1[13] * cospi_10_64;
  temp2 = step1[10] * cospi_10_64 + step1[13] * cospi_22_64;
  step2[10] = WRAPLOW(dct_const_round_shift(temp1), 8);
  step2[13] = WRAPLOW(dct_const_round_shift(temp2), 8);

  temp1 = step1[11] * cospi_6_64 - step1[12] * cospi_26_64;
  temp2 = step1[11] * cospi_26_64 + step1[12] * cospi_6_64;
  step2[11] = WRAPLOW(dct_const_round_shift(temp1), 8);
  step2[12] = WRAPLOW(dct_const_round_shift(temp2), 8);

  // stage 3
  step1[0] = step2[0];
  step1[1] = step2[1];
  step1[2] = step2[2];
  step1[3] = step2[3];

  temp1 = step2[4] * cospi_28_64 - step2[7] * cospi_4_64;
  temp2 = step2[4] * cospi_4_64 + step2[7] * cospi_28_64;
  step1[4] = WRAPLOW(dct_const_round_shift(temp1), 8);
  step1[7] = WRAPLOW(dct_const_round_shift(temp2), 8);
  temp1 = step2[5] * cospi_12_64 - step2[6] * cospi_20_64;
  temp2 = step2[5] * cospi_20_64 + step2[6] * cospi_12_64;
  step1[5] = WRAPLOW(dct_const_round_shift(temp1), 8);
  step1[6] = WRAPLOW(dct_const_round_shift(temp2), 8);

  step1[8] = WRAPLOW(step2[8] + step2[9], 8);
  step1[9] = WRAPLOW(step2[8] - step2[9], 8);
  step1[10] = WRAPLOW(-step2[10] + step2[11], 8);
  step1[11] = WRAPLOW(step2[10] + step2[11], 8);
  step1[12] = WRAPLOW(step2[12] + step2[13], 8);
  step1[13] = WRAPLOW(step2[12] - step2[13], 8);
  step1[14] = WRAPLOW(-step2[14] + step2[15], 8);
  step1[15] = WRAPLOW(step2[14] + step2[15], 8);

  // stage 4
  temp1 = (step1[0] + step1[1]) * cospi_16_64;
  temp2 = (step1[0] - step1[1]) * cospi_16_64;
  step2[0] = WRAPLOW(dct_const_round_shift(temp1), 8);
  step2[1] = WRAPLOW(dct_const_round_shift(temp2), 8);
  temp1 = step1[2] * cospi_24_64 - step1[3] * cospi_8_64;
  temp2 = step1[2] * cospi_8_64 + step1[3] * cospi_24_64;
  step2[2] = WRAPLOW(dct_const_round_shift(temp1), 8);
  step2[3] = WRAPLOW(dct_const_round_shift(temp2), 8);
  step2[4] = WRAPLOW(step1[4] + step1[5], 8);
  step2[5] = WRAPLOW(step1[4] - step1[5], 8);
  step2[6] = WRAPLOW(-step1[6] + step1[7], 8);
  step2[7] = WRAPLOW(step1[6] + step1[7], 8);

  step2[8] = step1[8];
  step2[15] = step1[15];
  temp1 = -step1[9] * cospi_8_64 + step1[14] * cospi_24_64;
  temp2 = step1[9] * cospi_24_64 + step1[14] * cospi_8_64;
  step2[9] = WRAPLOW(dct_const_round_shift(temp1), 8);
  step2[14] = WRAPLOW(dct_const_round_shift(temp2), 8);
  temp1 = -step1[10] * cospi_24_64 - step1[13] * cospi_8_64;
  temp2 = -step1[10] * cospi_8_64 + step1[13] * cospi_24_64;
  step2[10] = WRAPLOW(dct_const_round_shift(temp1), 8);
  step2[13] = WRAPLOW(dct_const_round_shift(temp2), 8);
  step2[11] = step1[11];
  step2[12] = step1[12];

  // stage 5
  step1[0] = WRAPLOW(step2[0] + step2[3], 8);
  step1[1] = WRAPLOW(step2[1] + step2[2], 8);
  step1[2] = WRAPLOW(step2[1] - step2[2], 8);
  step1[3] = WRAPLOW(step2[0] - step2[3], 8);
  step1[4] = step2[4];
  temp1 = (step2[6] - step2[5]) * cospi_16_64;
  temp2 = (step2[5] + step2[6]) * cospi_16_64;
  step1[5] = WRAPLOW(dct_const_round_shift(temp1), 8);
  step1[6] = WRAPLOW(dct_const_round_shift(temp2), 8);
  step1[7] = step2[7];

  step1[8] = WRAPLOW(step2[8] + step2[11], 8);
  step1[9] = WRAPLOW(step2[9] + step2[10], 8);
  step1[10] = WRAPLOW(step2[9] - step2[10], 8);
  step1[11] = WRAPLOW(step2[8] - step2[11], 8);
  step1[12] = WRAPLOW(-step2[12] + step2[15], 8);
  step1[13] = WRAPLOW(-step2[13] + step2[14], 8);
  step1[14] = WRAPLOW(step2[13] + step2[14], 8);
  step1[15] = WRAPLOW(step2[12] + step2[15], 8);

  // stage 6
  step2[0] = WRAPLOW(step1[0] + step1[7], 8);
  step2[1] = WRAPLOW(step1[1] + step1[6], 8);
  step2[2] = WRAPLOW(step1[2] + step1[5], 8);
  step2[3] = WRAPLOW(step1[3] + step1[4], 8);
  step2[4] = WRAPLOW(step1[3] - step1[4], 8);
  step2[5] = WRAPLOW(step1[2] - step1[5], 8);
  step2[6] = WRAPLOW(step1[1] - step1[6], 8);
  step2[7] = WRAPLOW(step1[0] - step1[7], 8);
  step2[8] = step1[8];
  step2[9] = step1[9];
  temp1 = (-step1[10] + step1[13]) * cospi_16_64;
  temp2 = (step1[10] + step1[13]) * cospi_16_64;
  step2[10] = WRAPLOW(dct_const_round_shift(temp1), 8);
  step2[13] = WRAPLOW(dct_const_round_shift(temp2), 8);
  temp1 = (-step1[11] + step1[12]) * cospi_16_64;
  temp2 = (step1[11] + step1[12]) * cospi_16_64;
  step2[11] = WRAPLOW(dct_const_round_shift(temp1), 8);
  step2[12] = WRAPLOW(dct_const_round_shift(temp2), 8);
  step2[14] = step1[14];
  step2[15] = step1[15];

  // stage 7
  output[0] = WRAPLOW(step2[0] + step2[15], 8);
  output[1] = WRAPLOW(step2[1] + step2[14], 8);
  output[2] = WRAPLOW(step2[2] + step2[13], 8);
  output[3] = WRAPLOW(step2[3] + step2[12], 8);
  output[4] = WRAPLOW(step2[4] + step2[11], 8);
  output[5] = WRAPLOW(step2[5] + step2[10], 8);
  output[6] = WRAPLOW(step2[6] + step2[9], 8);
  output[7] = WRAPLOW(step2[7] + step2[8], 8);
  output[8] = WRAPLOW(step2[7] - step2[8], 8);
  output[9] = WRAPLOW(step2[6] - step2[9], 8);
  output[10] = WRAPLOW(step2[5] - step2[10], 8);
  output[11] = WRAPLOW(step2[4] - step2[11], 8);
  output[12] = WRAPLOW(step2[3] - step2[12], 8);
  output[13] = WRAPLOW(step2[2] - step2[13], 8);
  output[14] = WRAPLOW(step2[1] - step2[14], 8);
  output[15] = WRAPLOW(step2[0] - step2[15], 8);
}

void vp9_idct16x16_256_add_c(const tran_low_t *input, uint8_t *dest,
                             int stride) {
  tran_low_t out[16 * 16];
  tran_low_t *outptr = out;
  int i, j;
  tran_low_t temp_in[16], temp_out[16];

  // First transform rows
  for (i = 0; i < 16; ++i) {
    idct16(input, outptr);
    input += 16;
    outptr += 16;
  }

  // Then transform columns
  for (i = 0; i < 16; ++i) {
    for (j = 0; j < 16; ++j)
      temp_in[j] = out[j * 16 + i];
    idct16(temp_in, temp_out);
    for (j = 0; j < 16; ++j) {
      dest[j * stride + i] = clip_pixel_add(dest[j * stride + i],
                                            ROUND_POWER_OF_TWO(temp_out[j], 6));
    }
  }
}

#if CONFIG_WAVELETS
void vp9_idct16x16_noscale_c(const tran_low_t *input, int16_t *dest,
                             int stride) {
  tran_low_t out[16 * 16];
  tran_low_t *outptr = out;
  int i, j;
  tran_low_t temp_in[16], temp_out[16];

  // First transform rows
  for (i = 0; i < 16; ++i) {
    idct16(input, outptr);
    input += 16;
    outptr += 16;
  }

  // Then transform columns
  for (i = 0; i < 16; ++i) {
    for (j = 0; j < 16; ++j)
      temp_in[j] = out[j * 16 + i];
    idct16(temp_in, temp_out);
    for (j = 0; j < 16; ++j) {
      dest[j * stride + i] = ROUND_POWER_OF_TWO(temp_out[j], 3);
    }
  }
}
#endif  // CONFIG_WAVELETS

static void iadst16(const tran_low_t *input, tran_low_t *output) {
  tran_high_t s0, s1, s2, s3, s4, s5, s6, s7, s8;
  tran_high_t s9, s10, s11, s12, s13, s14, s15;

  tran_high_t x0 = input[15];
  tran_high_t x1 = input[0];
  tran_high_t x2 = input[13];
  tran_high_t x3 = input[2];
  tran_high_t x4 = input[11];
  tran_high_t x5 = input[4];
  tran_high_t x6 = input[9];
  tran_high_t x7 = input[6];
  tran_high_t x8 = input[7];
  tran_high_t x9 = input[8];
  tran_high_t x10 = input[5];
  tran_high_t x11 = input[10];
  tran_high_t x12 = input[3];
  tran_high_t x13 = input[12];
  tran_high_t x14 = input[1];
  tran_high_t x15 = input[14];

  if (!(x0 | x1 | x2 | x3 | x4 | x5 | x6 | x7 | x8
           | x9 | x10 | x11 | x12 | x13 | x14 | x15)) {
    output[0] = output[1] = output[2] = output[3] = output[4]
              = output[5] = output[6] = output[7] = output[8]
              = output[9] = output[10] = output[11] = output[12]
              = output[13] = output[14] = output[15] = 0;
    return;
  }

  // stage 1
  s0 = x0 * cospi_1_64  + x1 * cospi_31_64;
  s1 = x0 * cospi_31_64 - x1 * cospi_1_64;
  s2 = x2 * cospi_5_64  + x3 * cospi_27_64;
  s3 = x2 * cospi_27_64 - x3 * cospi_5_64;
  s4 = x4 * cospi_9_64  + x5 * cospi_23_64;
  s5 = x4 * cospi_23_64 - x5 * cospi_9_64;
  s6 = x6 * cospi_13_64 + x7 * cospi_19_64;
  s7 = x6 * cospi_19_64 - x7 * cospi_13_64;
  s8 = x8 * cospi_17_64 + x9 * cospi_15_64;
  s9 = x8 * cospi_15_64 - x9 * cospi_17_64;
  s10 = x10 * cospi_21_64 + x11 * cospi_11_64;
  s11 = x10 * cospi_11_64 - x11 * cospi_21_64;
  s12 = x12 * cospi_25_64 + x13 * cospi_7_64;
  s13 = x12 * cospi_7_64  - x13 * cospi_25_64;
  s14 = x14 * cospi_29_64 + x15 * cospi_3_64;
  s15 = x14 * cospi_3_64  - x15 * cospi_29_64;

  x0 = WRAPLOW(dct_const_round_shift(s0 + s8), 8);
  x1 = WRAPLOW(dct_const_round_shift(s1 + s9), 8);
  x2 = WRAPLOW(dct_const_round_shift(s2 + s10), 8);
  x3 = WRAPLOW(dct_const_round_shift(s3 + s11), 8);
  x4 = WRAPLOW(dct_const_round_shift(s4 + s12), 8);
  x5 = WRAPLOW(dct_const_round_shift(s5 + s13), 8);
  x6 = WRAPLOW(dct_const_round_shift(s6 + s14), 8);
  x7 = WRAPLOW(dct_const_round_shift(s7 + s15), 8);
  x8 = WRAPLOW(dct_const_round_shift(s0 - s8), 8);
  x9 = WRAPLOW(dct_const_round_shift(s1 - s9), 8);
  x10 = WRAPLOW(dct_const_round_shift(s2 - s10), 8);
  x11 = WRAPLOW(dct_const_round_shift(s3 - s11), 8);
  x12 = WRAPLOW(dct_const_round_shift(s4 - s12), 8);
  x13 = WRAPLOW(dct_const_round_shift(s5 - s13), 8);
  x14 = WRAPLOW(dct_const_round_shift(s6 - s14), 8);
  x15 = WRAPLOW(dct_const_round_shift(s7 - s15), 8);

  // stage 2
  s0 = x0;
  s1 = x1;
  s2 = x2;
  s3 = x3;
  s4 = x4;
  s5 = x5;
  s6 = x6;
  s7 = x7;
  s8 =    x8 * cospi_4_64   + x9 * cospi_28_64;
  s9 =    x8 * cospi_28_64  - x9 * cospi_4_64;
  s10 =   x10 * cospi_20_64 + x11 * cospi_12_64;
  s11 =   x10 * cospi_12_64 - x11 * cospi_20_64;
  s12 = - x12 * cospi_28_64 + x13 * cospi_4_64;
  s13 =   x12 * cospi_4_64  + x13 * cospi_28_64;
  s14 = - x14 * cospi_12_64 + x15 * cospi_20_64;
  s15 =   x14 * cospi_20_64 + x15 * cospi_12_64;

  x0 = WRAPLOW(s0 + s4, 8);
  x1 = WRAPLOW(s1 + s5, 8);
  x2 = WRAPLOW(s2 + s6, 8);
  x3 = WRAPLOW(s3 + s7, 8);
  x4 = WRAPLOW(s0 - s4, 8);
  x5 = WRAPLOW(s1 - s5, 8);
  x6 = WRAPLOW(s2 - s6, 8);
  x7 = WRAPLOW(s3 - s7, 8);
  x8 = WRAPLOW(dct_const_round_shift(s8 + s12), 8);
  x9 = WRAPLOW(dct_const_round_shift(s9 + s13), 8);
  x10 = WRAPLOW(dct_const_round_shift(s10 + s14), 8);
  x11 = WRAPLOW(dct_const_round_shift(s11 + s15), 8);
  x12 = WRAPLOW(dct_const_round_shift(s8 - s12), 8);
  x13 = WRAPLOW(dct_const_round_shift(s9 - s13), 8);
  x14 = WRAPLOW(dct_const_round_shift(s10 - s14), 8);
  x15 = WRAPLOW(dct_const_round_shift(s11 - s15), 8);

  // stage 3
  s0 = x0;
  s1 = x1;
  s2 = x2;
  s3 = x3;
  s4 = x4 * cospi_8_64  + x5 * cospi_24_64;
  s5 = x4 * cospi_24_64 - x5 * cospi_8_64;
  s6 = - x6 * cospi_24_64 + x7 * cospi_8_64;
  s7 =   x6 * cospi_8_64  + x7 * cospi_24_64;
  s8 = x8;
  s9 = x9;
  s10 = x10;
  s11 = x11;
  s12 = x12 * cospi_8_64  + x13 * cospi_24_64;
  s13 = x12 * cospi_24_64 - x13 * cospi_8_64;
  s14 = - x14 * cospi_24_64 + x15 * cospi_8_64;
  s15 =   x14 * cospi_8_64  + x15 * cospi_24_64;

  x0 = WRAPLOW(check_range(s0 + s2), 8);
  x1 = WRAPLOW(check_range(s1 + s3), 8);
  x2 = WRAPLOW(check_range(s0 - s2), 8);
  x3 = WRAPLOW(check_range(s1 - s3), 8);
  x4 = WRAPLOW(dct_const_round_shift(s4 + s6), 8);
  x5 = WRAPLOW(dct_const_round_shift(s5 + s7), 8);
  x6 = WRAPLOW(dct_const_round_shift(s4 - s6), 8);
  x7 = WRAPLOW(dct_const_round_shift(s5 - s7), 8);
  x8 = WRAPLOW(check_range(s8 + s10), 8);
  x9 = WRAPLOW(check_range(s9 + s11), 8);
  x10 = WRAPLOW(check_range(s8 - s10), 8);
  x11 = WRAPLOW(check_range(s9 - s11), 8);
  x12 = WRAPLOW(dct_const_round_shift(s12 + s14), 8);
  x13 = WRAPLOW(dct_const_round_shift(s13 + s15), 8);
  x14 = WRAPLOW(dct_const_round_shift(s12 - s14), 8);
  x15 = WRAPLOW(dct_const_round_shift(s13 - s15), 8);

  // stage 4
  s2 = (- cospi_16_64) * (x2 + x3);
  s3 = cospi_16_64 * (x2 - x3);
  s6 = cospi_16_64 * (x6 + x7);
  s7 = cospi_16_64 * (- x6 + x7);
  s10 = cospi_16_64 * (x10 + x11);
  s11 = cospi_16_64 * (- x10 + x11);
  s14 = (- cospi_16_64) * (x14 + x15);
  s15 = cospi_16_64 * (x14 - x15);

  x2 = WRAPLOW(dct_const_round_shift(s2), 8);
  x3 = WRAPLOW(dct_const_round_shift(s3), 8);
  x6 = WRAPLOW(dct_const_round_shift(s6), 8);
  x7 = WRAPLOW(dct_const_round_shift(s7), 8);
  x10 = WRAPLOW(dct_const_round_shift(s10), 8);
  x11 = WRAPLOW(dct_const_round_shift(s11), 8);
  x14 = WRAPLOW(dct_const_round_shift(s14), 8);
  x15 = WRAPLOW(dct_const_round_shift(s15), 8);

  output[0] = WRAPLOW(x0, 8);
  output[1] = WRAPLOW(-x8, 8);
  output[2] = WRAPLOW(x12, 8);
  output[3] = WRAPLOW(-x4, 8);
  output[4] = WRAPLOW(x6, 8);
  output[5] = WRAPLOW(x14, 8);
  output[6] = WRAPLOW(x10, 8);
  output[7] = WRAPLOW(x2, 8);
  output[8] = WRAPLOW(x3, 8);
  output[9] = WRAPLOW(x11, 8);
  output[10] = WRAPLOW(x15, 8);
  output[11] = WRAPLOW(x7, 8);
  output[12] = WRAPLOW(x5, 8);
  output[13] = WRAPLOW(-x13, 8);
  output[14] = WRAPLOW(x9, 8);
  output[15] = WRAPLOW(-x1, 8);
}

static const transform_2d IHT_16[] = {
  { idct16,  idct16  },  // DCT_DCT  = 0
  { iadst16, idct16  },  // ADST_DCT = 1
  { idct16,  iadst16 },  // DCT_ADST = 2
  { iadst16, iadst16 }   // ADST_ADST = 3
};

void vp9_iht16x16_256_add_c(const tran_low_t *input, uint8_t *dest, int stride,
                            int tx_type) {
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

void vp9_idct16x16_10_add_c(const tran_low_t *input, uint8_t *dest,
                            int stride) {
  tran_low_t out[16 * 16] = { 0 };
  tran_low_t *outptr = out;
  int i, j;
  tran_low_t temp_in[16], temp_out[16];

  // First transform rows. Since all non-zero dct coefficients are in
  // upper-left 4x4 area, we only need to calculate first 4 rows here.
  for (i = 0; i < 4; ++i) {
    idct16(input, outptr);
    input += 16;
    outptr += 16;
  }

  // Then transform columns
  for (i = 0; i < 16; ++i) {
    for (j = 0; j < 16; ++j)
      temp_in[j] = out[j*16 + i];
    idct16(temp_in, temp_out);
    for (j = 0; j < 16; ++j) {
      dest[j * stride + i] = clip_pixel_add(dest[j * stride + i],
                                            ROUND_POWER_OF_TWO(temp_out[j], 6));
    }
  }
}

void vp9_idct16x16_1_add_c(const tran_low_t *input, uint8_t *dest, int stride) {
  int i, j;
  tran_high_t a1;
  tran_low_t out = WRAPLOW(dct_const_round_shift(input[0] * cospi_16_64), 8);
  out = WRAPLOW(dct_const_round_shift(out * cospi_16_64), 8);
  a1 = ROUND_POWER_OF_TWO(out, 6);
  for (j = 0; j < 16; ++j) {
    for (i = 0; i < 16; ++i)
      dest[i] = clip_pixel_add(dest[i], a1);
    dest += stride;
  }
}

static void idct32(const tran_low_t *input, tran_low_t *output) {
  tran_low_t step1[32], step2[32];
  tran_high_t temp1, temp2;

  // stage 1
  step1[0] = input[0];
  step1[1] = input[16];
  step1[2] = input[8];
  step1[3] = input[24];
  step1[4] = input[4];
  step1[5] = input[20];
  step1[6] = input[12];
  step1[7] = input[28];
  step1[8] = input[2];
  step1[9] = input[18];
  step1[10] = input[10];
  step1[11] = input[26];
  step1[12] = input[6];
  step1[13] = input[22];
  step1[14] = input[14];
  step1[15] = input[30];

  temp1 = input[1] * cospi_31_64 - input[31] * cospi_1_64;
  temp2 = input[1] * cospi_1_64 + input[31] * cospi_31_64;
  step1[16] = WRAPLOW(dct_const_round_shift(temp1), 8);
  step1[31] = WRAPLOW(dct_const_round_shift(temp2), 8);

  temp1 = input[17] * cospi_15_64 - input[15] * cospi_17_64;
  temp2 = input[17] * cospi_17_64 + input[15] * cospi_15_64;
  step1[17] = WRAPLOW(dct_const_round_shift(temp1), 8);
  step1[30] = WRAPLOW(dct_const_round_shift(temp2), 8);

  temp1 = input[9] * cospi_23_64 - input[23] * cospi_9_64;
  temp2 = input[9] * cospi_9_64 + input[23] * cospi_23_64;
  step1[18] = WRAPLOW(dct_const_round_shift(temp1), 8);
  step1[29] = WRAPLOW(dct_const_round_shift(temp2), 8);

  temp1 = input[25] * cospi_7_64 - input[7] * cospi_25_64;
  temp2 = input[25] * cospi_25_64 + input[7] * cospi_7_64;
  step1[19] = WRAPLOW(dct_const_round_shift(temp1), 8);
  step1[28] = WRAPLOW(dct_const_round_shift(temp2), 8);

  temp1 = input[5] * cospi_27_64 - input[27] * cospi_5_64;
  temp2 = input[5] * cospi_5_64 + input[27] * cospi_27_64;
  step1[20] = WRAPLOW(dct_const_round_shift(temp1), 8);
  step1[27] = WRAPLOW(dct_const_round_shift(temp2), 8);

  temp1 = input[21] * cospi_11_64 - input[11] * cospi_21_64;
  temp2 = input[21] * cospi_21_64 + input[11] * cospi_11_64;
  step1[21] = WRAPLOW(dct_const_round_shift(temp1), 8);
  step1[26] = WRAPLOW(dct_const_round_shift(temp2), 8);

  temp1 = input[13] * cospi_19_64 - input[19] * cospi_13_64;
  temp2 = input[13] * cospi_13_64 + input[19] * cospi_19_64;
  step1[22] = WRAPLOW(dct_const_round_shift(temp1), 8);
  step1[25] = WRAPLOW(dct_const_round_shift(temp2), 8);

  temp1 = input[29] * cospi_3_64 - input[3] * cospi_29_64;
  temp2 = input[29] * cospi_29_64 + input[3] * cospi_3_64;
  step1[23] = WRAPLOW(dct_const_round_shift(temp1), 8);
  step1[24] = WRAPLOW(dct_const_round_shift(temp2), 8);

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
  step2[8] = WRAPLOW(dct_const_round_shift(temp1), 8);
  step2[15] = WRAPLOW(dct_const_round_shift(temp2), 8);

  temp1 = step1[9] * cospi_14_64 - step1[14] * cospi_18_64;
  temp2 = step1[9] * cospi_18_64 + step1[14] * cospi_14_64;
  step2[9] = WRAPLOW(dct_const_round_shift(temp1), 8);
  step2[14] = WRAPLOW(dct_const_round_shift(temp2), 8);

  temp1 = step1[10] * cospi_22_64 - step1[13] * cospi_10_64;
  temp2 = step1[10] * cospi_10_64 + step1[13] * cospi_22_64;
  step2[10] = WRAPLOW(dct_const_round_shift(temp1), 8);
  step2[13] = WRAPLOW(dct_const_round_shift(temp2), 8);

  temp1 = step1[11] * cospi_6_64 - step1[12] * cospi_26_64;
  temp2 = step1[11] * cospi_26_64 + step1[12] * cospi_6_64;
  step2[11] = WRAPLOW(dct_const_round_shift(temp1), 8);
  step2[12] = WRAPLOW(dct_const_round_shift(temp2), 8);

  step2[16] = WRAPLOW(step1[16] + step1[17], 8);
  step2[17] = WRAPLOW(step1[16] - step1[17], 8);
  step2[18] = WRAPLOW(-step1[18] + step1[19], 8);
  step2[19] = WRAPLOW(step1[18] + step1[19], 8);
  step2[20] = WRAPLOW(step1[20] + step1[21], 8);
  step2[21] = WRAPLOW(step1[20] - step1[21], 8);
  step2[22] = WRAPLOW(-step1[22] + step1[23], 8);
  step2[23] = WRAPLOW(step1[22] + step1[23], 8);
  step2[24] = WRAPLOW(step1[24] + step1[25], 8);
  step2[25] = WRAPLOW(step1[24] - step1[25], 8);
  step2[26] = WRAPLOW(-step1[26] + step1[27], 8);
  step2[27] = WRAPLOW(step1[26] + step1[27], 8);
  step2[28] = WRAPLOW(step1[28] + step1[29], 8);
  step2[29] = WRAPLOW(step1[28] - step1[29], 8);
  step2[30] = WRAPLOW(-step1[30] + step1[31], 8);
  step2[31] = WRAPLOW(step1[30] + step1[31], 8);

  // stage 3
  step1[0] = step2[0];
  step1[1] = step2[1];
  step1[2] = step2[2];
  step1[3] = step2[3];

  temp1 = step2[4] * cospi_28_64 - step2[7] * cospi_4_64;
  temp2 = step2[4] * cospi_4_64 + step2[7] * cospi_28_64;
  step1[4] = WRAPLOW(dct_const_round_shift(temp1), 8);
  step1[7] = WRAPLOW(dct_const_round_shift(temp2), 8);
  temp1 = step2[5] * cospi_12_64 - step2[6] * cospi_20_64;
  temp2 = step2[5] * cospi_20_64 + step2[6] * cospi_12_64;
  step1[5] = WRAPLOW(dct_const_round_shift(temp1), 8);
  step1[6] = WRAPLOW(dct_const_round_shift(temp2), 8);

  step1[8] = WRAPLOW(step2[8] + step2[9], 8);
  step1[9] = WRAPLOW(step2[8] - step2[9], 8);
  step1[10] = WRAPLOW(-step2[10] + step2[11], 8);
  step1[11] = WRAPLOW(step2[10] + step2[11], 8);
  step1[12] = WRAPLOW(step2[12] + step2[13], 8);
  step1[13] = WRAPLOW(step2[12] - step2[13], 8);
  step1[14] = WRAPLOW(-step2[14] + step2[15], 8);
  step1[15] = WRAPLOW(step2[14] + step2[15], 8);

  step1[16] = step2[16];
  step1[31] = step2[31];
  temp1 = -step2[17] * cospi_4_64 + step2[30] * cospi_28_64;
  temp2 = step2[17] * cospi_28_64 + step2[30] * cospi_4_64;
  step1[17] = WRAPLOW(dct_const_round_shift(temp1), 8);
  step1[30] = WRAPLOW(dct_const_round_shift(temp2), 8);
  temp1 = -step2[18] * cospi_28_64 - step2[29] * cospi_4_64;
  temp2 = -step2[18] * cospi_4_64 + step2[29] * cospi_28_64;
  step1[18] = WRAPLOW(dct_const_round_shift(temp1), 8);
  step1[29] = WRAPLOW(dct_const_round_shift(temp2), 8);
  step1[19] = step2[19];
  step1[20] = step2[20];
  temp1 = -step2[21] * cospi_20_64 + step2[26] * cospi_12_64;
  temp2 = step2[21] * cospi_12_64 + step2[26] * cospi_20_64;
  step1[21] = WRAPLOW(dct_const_round_shift(temp1), 8);
  step1[26] = WRAPLOW(dct_const_round_shift(temp2), 8);
  temp1 = -step2[22] * cospi_12_64 - step2[25] * cospi_20_64;
  temp2 = -step2[22] * cospi_20_64 + step2[25] * cospi_12_64;
  step1[22] = WRAPLOW(dct_const_round_shift(temp1), 8);
  step1[25] = WRAPLOW(dct_const_round_shift(temp2), 8);
  step1[23] = step2[23];
  step1[24] = step2[24];
  step1[27] = step2[27];
  step1[28] = step2[28];

  // stage 4
  temp1 = (step1[0] + step1[1]) * cospi_16_64;
  temp2 = (step1[0] - step1[1]) * cospi_16_64;
  step2[0] = WRAPLOW(dct_const_round_shift(temp1), 8);
  step2[1] = WRAPLOW(dct_const_round_shift(temp2), 8);
  temp1 = step1[2] * cospi_24_64 - step1[3] * cospi_8_64;
  temp2 = step1[2] * cospi_8_64 + step1[3] * cospi_24_64;
  step2[2] = WRAPLOW(dct_const_round_shift(temp1), 8);
  step2[3] = WRAPLOW(dct_const_round_shift(temp2), 8);
  step2[4] = WRAPLOW(step1[4] + step1[5], 8);
  step2[5] = WRAPLOW(step1[4] - step1[5], 8);
  step2[6] = WRAPLOW(-step1[6] + step1[7], 8);
  step2[7] = WRAPLOW(step1[6] + step1[7], 8);

  step2[8] = step1[8];
  step2[15] = step1[15];
  temp1 = -step1[9] * cospi_8_64 + step1[14] * cospi_24_64;
  temp2 = step1[9] * cospi_24_64 + step1[14] * cospi_8_64;
  step2[9] = WRAPLOW(dct_const_round_shift(temp1), 8);
  step2[14] = WRAPLOW(dct_const_round_shift(temp2), 8);
  temp1 = -step1[10] * cospi_24_64 - step1[13] * cospi_8_64;
  temp2 = -step1[10] * cospi_8_64 + step1[13] * cospi_24_64;
  step2[10] = WRAPLOW(dct_const_round_shift(temp1), 8);
  step2[13] = WRAPLOW(dct_const_round_shift(temp2), 8);
  step2[11] = step1[11];
  step2[12] = step1[12];

  step2[16] = WRAPLOW(step1[16] + step1[19], 8);
  step2[17] = WRAPLOW(step1[17] + step1[18], 8);
  step2[18] = WRAPLOW(step1[17] - step1[18], 8);
  step2[19] = WRAPLOW(step1[16] - step1[19], 8);
  step2[20] = WRAPLOW(-step1[20] + step1[23], 8);
  step2[21] = WRAPLOW(-step1[21] + step1[22], 8);
  step2[22] = WRAPLOW(step1[21] + step1[22], 8);
  step2[23] = WRAPLOW(step1[20] + step1[23], 8);

  step2[24] = WRAPLOW(step1[24] + step1[27], 8);
  step2[25] = WRAPLOW(step1[25] + step1[26], 8);
  step2[26] = WRAPLOW(step1[25] - step1[26], 8);
  step2[27] = WRAPLOW(step1[24] - step1[27], 8);
  step2[28] = WRAPLOW(-step1[28] + step1[31], 8);
  step2[29] = WRAPLOW(-step1[29] + step1[30], 8);
  step2[30] = WRAPLOW(step1[29] + step1[30], 8);
  step2[31] = WRAPLOW(step1[28] + step1[31], 8);

  // stage 5
  step1[0] = WRAPLOW(step2[0] + step2[3], 8);
  step1[1] = WRAPLOW(step2[1] + step2[2], 8);
  step1[2] = WRAPLOW(step2[1] - step2[2], 8);
  step1[3] = WRAPLOW(step2[0] - step2[3], 8);
  step1[4] = step2[4];
  temp1 = (step2[6] - step2[5]) * cospi_16_64;
  temp2 = (step2[5] + step2[6]) * cospi_16_64;
  step1[5] = WRAPLOW(dct_const_round_shift(temp1), 8);
  step1[6] = WRAPLOW(dct_const_round_shift(temp2), 8);
  step1[7] = step2[7];

  step1[8] = WRAPLOW(step2[8] + step2[11], 8);
  step1[9] = WRAPLOW(step2[9] + step2[10], 8);
  step1[10] = WRAPLOW(step2[9] - step2[10], 8);
  step1[11] = WRAPLOW(step2[8] - step2[11], 8);
  step1[12] = WRAPLOW(-step2[12] + step2[15], 8);
  step1[13] = WRAPLOW(-step2[13] + step2[14], 8);
  step1[14] = WRAPLOW(step2[13] + step2[14], 8);
  step1[15] = WRAPLOW(step2[12] + step2[15], 8);

  step1[16] = step2[16];
  step1[17] = step2[17];
  temp1 = -step2[18] * cospi_8_64 + step2[29] * cospi_24_64;
  temp2 = step2[18] * cospi_24_64 + step2[29] * cospi_8_64;
  step1[18] = WRAPLOW(dct_const_round_shift(temp1), 8);
  step1[29] = WRAPLOW(dct_const_round_shift(temp2), 8);
  temp1 = -step2[19] * cospi_8_64 + step2[28] * cospi_24_64;
  temp2 = step2[19] * cospi_24_64 + step2[28] * cospi_8_64;
  step1[19] = WRAPLOW(dct_const_round_shift(temp1), 8);
  step1[28] = WRAPLOW(dct_const_round_shift(temp2), 8);
  temp1 = -step2[20] * cospi_24_64 - step2[27] * cospi_8_64;
  temp2 = -step2[20] * cospi_8_64 + step2[27] * cospi_24_64;
  step1[20] = WRAPLOW(dct_const_round_shift(temp1), 8);
  step1[27] = WRAPLOW(dct_const_round_shift(temp2), 8);
  temp1 = -step2[21] * cospi_24_64 - step2[26] * cospi_8_64;
  temp2 = -step2[21] * cospi_8_64 + step2[26] * cospi_24_64;
  step1[21] = WRAPLOW(dct_const_round_shift(temp1), 8);
  step1[26] = WRAPLOW(dct_const_round_shift(temp2), 8);
  step1[22] = step2[22];
  step1[23] = step2[23];
  step1[24] = step2[24];
  step1[25] = step2[25];
  step1[30] = step2[30];
  step1[31] = step2[31];

  // stage 6
  step2[0] = WRAPLOW(step1[0] + step1[7], 8);
  step2[1] = WRAPLOW(step1[1] + step1[6], 8);
  step2[2] = WRAPLOW(step1[2] + step1[5], 8);
  step2[3] = WRAPLOW(step1[3] + step1[4], 8);
  step2[4] = WRAPLOW(step1[3] - step1[4], 8);
  step2[5] = WRAPLOW(step1[2] - step1[5], 8);
  step2[6] = WRAPLOW(step1[1] - step1[6], 8);
  step2[7] = WRAPLOW(step1[0] - step1[7], 8);
  step2[8] = step1[8];
  step2[9] = step1[9];
  temp1 = (-step1[10] + step1[13]) * cospi_16_64;
  temp2 = (step1[10] + step1[13]) * cospi_16_64;
  step2[10] = WRAPLOW(dct_const_round_shift(temp1), 8);
  step2[13] = WRAPLOW(dct_const_round_shift(temp2), 8);
  temp1 = (-step1[11] + step1[12]) * cospi_16_64;
  temp2 = (step1[11] + step1[12]) * cospi_16_64;
  step2[11] = WRAPLOW(dct_const_round_shift(temp1), 8);
  step2[12] = WRAPLOW(dct_const_round_shift(temp2), 8);
  step2[14] = step1[14];
  step2[15] = step1[15];

  step2[16] = WRAPLOW(step1[16] + step1[23], 8);
  step2[17] = WRAPLOW(step1[17] + step1[22], 8);
  step2[18] = WRAPLOW(step1[18] + step1[21], 8);
  step2[19] = WRAPLOW(step1[19] + step1[20], 8);
  step2[20] = WRAPLOW(step1[19] - step1[20], 8);
  step2[21] = WRAPLOW(step1[18] - step1[21], 8);
  step2[22] = WRAPLOW(step1[17] - step1[22], 8);
  step2[23] = WRAPLOW(step1[16] - step1[23], 8);

  step2[24] = WRAPLOW(-step1[24] + step1[31], 8);
  step2[25] = WRAPLOW(-step1[25] + step1[30], 8);
  step2[26] = WRAPLOW(-step1[26] + step1[29], 8);
  step2[27] = WRAPLOW(-step1[27] + step1[28], 8);
  step2[28] = WRAPLOW(step1[27] + step1[28], 8);
  step2[29] = WRAPLOW(step1[26] + step1[29], 8);
  step2[30] = WRAPLOW(step1[25] + step1[30], 8);
  step2[31] = WRAPLOW(step1[24] + step1[31], 8);

  // stage 7
  step1[0] = WRAPLOW(step2[0] + step2[15], 8);
  step1[1] = WRAPLOW(step2[1] + step2[14], 8);
  step1[2] = WRAPLOW(step2[2] + step2[13], 8);
  step1[3] = WRAPLOW(step2[3] + step2[12], 8);
  step1[4] = WRAPLOW(step2[4] + step2[11], 8);
  step1[5] = WRAPLOW(step2[5] + step2[10], 8);
  step1[6] = WRAPLOW(step2[6] + step2[9], 8);
  step1[7] = WRAPLOW(step2[7] + step2[8], 8);
  step1[8] = WRAPLOW(step2[7] - step2[8], 8);
  step1[9] = WRAPLOW(step2[6] - step2[9], 8);
  step1[10] = WRAPLOW(step2[5] - step2[10], 8);
  step1[11] = WRAPLOW(step2[4] - step2[11], 8);
  step1[12] = WRAPLOW(step2[3] - step2[12], 8);
  step1[13] = WRAPLOW(step2[2] - step2[13], 8);
  step1[14] = WRAPLOW(step2[1] - step2[14], 8);
  step1[15] = WRAPLOW(step2[0] - step2[15], 8);

  step1[16] = step2[16];
  step1[17] = step2[17];
  step1[18] = step2[18];
  step1[19] = step2[19];
  temp1 = (-step2[20] + step2[27]) * cospi_16_64;
  temp2 = (step2[20] + step2[27]) * cospi_16_64;
  step1[20] = WRAPLOW(dct_const_round_shift(temp1), 8);
  step1[27] = WRAPLOW(dct_const_round_shift(temp2), 8);
  temp1 = (-step2[21] + step2[26]) * cospi_16_64;
  temp2 = (step2[21] + step2[26]) * cospi_16_64;
  step1[21] = WRAPLOW(dct_const_round_shift(temp1), 8);
  step1[26] = WRAPLOW(dct_const_round_shift(temp2), 8);
  temp1 = (-step2[22] + step2[25]) * cospi_16_64;
  temp2 = (step2[22] + step2[25]) * cospi_16_64;
  step1[22] = WRAPLOW(dct_const_round_shift(temp1), 8);
  step1[25] = WRAPLOW(dct_const_round_shift(temp2), 8);
  temp1 = (-step2[23] + step2[24]) * cospi_16_64;
  temp2 = (step2[23] + step2[24]) * cospi_16_64;
  step1[23] = WRAPLOW(dct_const_round_shift(temp1), 8);
  step1[24] = WRAPLOW(dct_const_round_shift(temp2), 8);
  step1[28] = step2[28];
  step1[29] = step2[29];
  step1[30] = step2[30];
  step1[31] = step2[31];

  // final stage
  output[0] = WRAPLOW(step1[0] + step1[31], 8);
  output[1] = WRAPLOW(step1[1] + step1[30], 8);
  output[2] = WRAPLOW(step1[2] + step1[29], 8);
  output[3] = WRAPLOW(step1[3] + step1[28], 8);
  output[4] = WRAPLOW(step1[4] + step1[27], 8);
  output[5] = WRAPLOW(step1[5] + step1[26], 8);
  output[6] = WRAPLOW(step1[6] + step1[25], 8);
  output[7] = WRAPLOW(step1[7] + step1[24], 8);
  output[8] = WRAPLOW(step1[8] + step1[23], 8);
  output[9] = WRAPLOW(step1[9] + step1[22], 8);
  output[10] = WRAPLOW(step1[10] + step1[21], 8);
  output[11] = WRAPLOW(step1[11] + step1[20], 8);
  output[12] = WRAPLOW(step1[12] + step1[19], 8);
  output[13] = WRAPLOW(step1[13] + step1[18], 8);
  output[14] = WRAPLOW(step1[14] + step1[17], 8);
  output[15] = WRAPLOW(step1[15] + step1[16], 8);
  output[16] = WRAPLOW(step1[15] - step1[16], 8);
  output[17] = WRAPLOW(step1[14] - step1[17], 8);
  output[18] = WRAPLOW(step1[13] - step1[18], 8);
  output[19] = WRAPLOW(step1[12] - step1[19], 8);
  output[20] = WRAPLOW(step1[11] - step1[20], 8);
  output[21] = WRAPLOW(step1[10] - step1[21], 8);
  output[22] = WRAPLOW(step1[9] - step1[22], 8);
  output[23] = WRAPLOW(step1[8] - step1[23], 8);
  output[24] = WRAPLOW(step1[7] - step1[24], 8);
  output[25] = WRAPLOW(step1[6] - step1[25], 8);
  output[26] = WRAPLOW(step1[5] - step1[26], 8);
  output[27] = WRAPLOW(step1[4] - step1[27], 8);
  output[28] = WRAPLOW(step1[3] - step1[28], 8);
  output[29] = WRAPLOW(step1[2] - step1[29], 8);
  output[30] = WRAPLOW(step1[1] - step1[30], 8);
  output[31] = WRAPLOW(step1[0] - step1[31], 8);
}

void vp9_idct32x32_1024_add_c(const tran_low_t *input, uint8_t *dest,
                              int stride) {
  tran_low_t out[32 * 32];
  tran_low_t *outptr = out;
  int i, j;
  tran_low_t temp_in[32], temp_out[32];

  // Rows
  for (i = 0; i < 32; ++i) {
    int16_t zero_coeff[16];
    for (j = 0; j < 16; ++j)
      zero_coeff[j] = input[2 * j] | input[2 * j + 1];
    for (j = 0; j < 8; ++j)
      zero_coeff[j] = zero_coeff[2 * j] | zero_coeff[2 * j + 1];
    for (j = 0; j < 4; ++j)
      zero_coeff[j] = zero_coeff[2 * j] | zero_coeff[2 * j + 1];
    for (j = 0; j < 2; ++j)
      zero_coeff[j] = zero_coeff[2 * j] | zero_coeff[2 * j + 1];

    if (zero_coeff[0] | zero_coeff[1])
      idct32(input, outptr);
    else
      vpx_memset(outptr, 0, sizeof(tran_low_t) * 32);
    input += 32;
    outptr += 32;
  }

  // Columns
  for (i = 0; i < 32; ++i) {
    for (j = 0; j < 32; ++j)
      temp_in[j] = out[j * 32 + i];
    idct32(temp_in, temp_out);
    for (j = 0; j < 32; ++j) {
      dest[j * stride + i] = clip_pixel_add(dest[j * stride + i],
                                            ROUND_POWER_OF_TWO(temp_out[j], 6));
    }
  }
}

#if CONFIG_WAVELETS
void vp9_idct32x32_noscale_c(const tran_low_t *input, int16_t *dest,
                             int stride) {
  tran_low_t out[32 * 32];
  tran_low_t *outptr = out;
  int i, j;
  tran_low_t temp_in[32], temp_out[32];

  // Rows
  for (i = 0; i < 32; ++i) {
    int16_t zero_coeff[16];
    for (j = 0; j < 16; ++j)
      zero_coeff[j] = input[2 * j] | input[2 * j + 1];
    for (j = 0; j < 8; ++j)
      zero_coeff[j] = zero_coeff[2 * j] | zero_coeff[2 * j + 1];
    for (j = 0; j < 4; ++j)
      zero_coeff[j] = zero_coeff[2 * j] | zero_coeff[2 * j + 1];
    for (j = 0; j < 2; ++j)
      zero_coeff[j] = zero_coeff[2 * j] | zero_coeff[2 * j + 1];

    if (zero_coeff[0] | zero_coeff[1])
      idct32(input, outptr);
    else
      vpx_memset(outptr, 0, sizeof(tran_low_t) * 32);
    input += 32;
    outptr += 32;
  }

  // Columns
  for (i = 0; i < 32; ++i) {
    for (j = 0; j < 32; ++j)
      temp_in[j] = out[j * 32 + i];
    idct32(temp_in, temp_out);
    for (j = 0; j < 32; ++j) {
      dest[j * stride + i] = ROUND_POWER_OF_TWO(temp_out[j], 4);
    }
  }
}
#endif  // CONFIG_WAVELETS

void vp9_idct32x32_34_add_c(const tran_low_t *input, uint8_t *dest,
                            int stride) {
  tran_low_t out[32 * 32] = {0};
  tran_low_t *outptr = out;
  int i, j;
  tran_low_t temp_in[32], temp_out[32];

  // Rows
  // only upper-left 8x8 has non-zero coeff
  for (i = 0; i < 8; ++i) {
    idct32(input, outptr);
    input += 32;
    outptr += 32;
  }

  // Columns
  for (i = 0; i < 32; ++i) {
    for (j = 0; j < 32; ++j)
      temp_in[j] = out[j * 32 + i];
    idct32(temp_in, temp_out);
    for (j = 0; j < 32; ++j) {
      dest[j * stride + i] = clip_pixel_add(dest[j * stride + i],
                                            ROUND_POWER_OF_TWO(temp_out[j], 6));
    }
  }
}

void vp9_idct32x32_1_add_c(const tran_low_t *input, uint8_t *dest, int stride) {
  int i, j;
  tran_high_t a1;

  tran_low_t out = WRAPLOW(dct_const_round_shift(input[0] * cospi_16_64), 8);
  out = WRAPLOW(dct_const_round_shift(out * cospi_16_64), 8);
  a1 = ROUND_POWER_OF_TWO(out, 6);

  for (j = 0; j < 32; ++j) {
    for (i = 0; i < 32; ++i)
      dest[i] = clip_pixel_add(dest[i], a1);
    dest += stride;
  }
}

// idct
void vp9_idct4x4_add(const tran_low_t *input, uint8_t *dest, int stride,
                     int eob) {
  if (eob > 1)
    vp9_idct4x4_16_add(input, dest, stride);
  else
    vp9_idct4x4_1_add(input, dest, stride);
}


void vp9_iwht4x4_add(const tran_low_t *input, uint8_t *dest, int stride,
                     int eob) {
  if (eob > 1)
    vp9_iwht4x4_16_add(input, dest, stride);
  else
    vp9_iwht4x4_1_add(input, dest, stride);
}

void vp9_idct8x8_add(const tran_low_t *input, uint8_t *dest, int stride,
                     int eob) {
  // If dc is 1, then input[0] is the reconstructed value, do not need
  // dequantization. Also, when dc is 1, dc is counted in eobs, namely eobs >=1.

  // The calculation can be simplified if there are not many non-zero dct
  // coefficients. Use eobs to decide what to do.
  // TODO(yunqingwang): "eobs = 1" case is also handled in vp9_short_idct8x8_c.
  // Combine that with code here.
  if (eob == 1)
    // DC only DCT coefficient
    vp9_idct8x8_1_add(input, dest, stride);
  else if (eob <= 12)
    vp9_idct8x8_12_add(input, dest, stride);
  else
    vp9_idct8x8_64_add(input, dest, stride);
}

void vp9_idct16x16_add(const tran_low_t *input, uint8_t *dest, int stride,
                       int eob) {
  /* The calculation can be simplified if there are not many non-zero dct
   * coefficients. Use eobs to separate different cases. */
  if (eob == 1)
    /* DC only DCT coefficient. */
    vp9_idct16x16_1_add(input, dest, stride);
  else if (eob <= 10)
    vp9_idct16x16_10_add(input, dest, stride);
  else
    vp9_idct16x16_256_add(input, dest, stride);
}

void vp9_idct32x32_add(const tran_low_t *input, uint8_t *dest, int stride,
                       int eob) {
  if (eob == 1)
    vp9_idct32x32_1_add(input, dest, stride);
  else if (eob <= 34)
    // non-zero coeff only in upper-left 8x8
    vp9_idct32x32_34_add(input, dest, stride);
  else
    vp9_idct32x32_1024_add(input, dest, stride);
}

// iht
void vp9_iht4x4_add(TX_TYPE tx_type, const tran_low_t *input, uint8_t *dest,
                    int stride, int eob) {
  if (tx_type == DCT_DCT) {
    vp9_idct4x4_add(input, dest, stride, eob);
#if CONFIG_EXT_TX
#if CONFIG_DST1
  } else if (tx_type == DST_DST) {
    vp9_idst4x4_add(input, dest, stride);
#endif  // CONFIG_DST1
  } else if (tx_type == FLIPADST_DCT) {
    flipud(dest, stride, 4);
    vp9_iht4x4_16_add(input, dest, stride, ADST_DCT);
    flipud(dest, stride, 4);
  } else if (tx_type == DCT_FLIPADST) {
    fliplr(dest, stride, 4);
    vp9_iht4x4_16_add(input, dest, stride, DCT_ADST);
    fliplr(dest, stride, 4);
  } else if (tx_type == FLIPADST_FLIPADST) {
    fliplrud(dest, stride, 4);
    vp9_iht4x4_16_add(input, dest, stride, ADST_ADST);
    fliplrud(dest, stride, 4);
  } else if (tx_type == ADST_FLIPADST) {
    fliplr(dest, stride, 4);
    vp9_iht4x4_16_add(input, dest, stride, ADST_ADST);
    fliplr(dest, stride, 4);
  } else if (tx_type == FLIPADST_ADST) {
    flipud(dest, stride, 4);
    vp9_iht4x4_16_add(input, dest, stride, ADST_ADST);
    flipud(dest, stride, 4);
#endif  // CONFIG_EXT_TX
  } else {
    vp9_iht4x4_16_add(input, dest, stride, tx_type);
  }
}

void vp9_iht8x8_add(TX_TYPE tx_type, const tran_low_t *input, uint8_t *dest,
                    int stride, int eob) {
  if (tx_type == DCT_DCT) {
    vp9_idct8x8_add(input, dest, stride, eob);
#if CONFIG_EXT_TX
#if CONFIG_DST1
  } else if (tx_type == DST_DST) {
    vp9_idst8x8_add(input, dest, stride);
#endif  // CONFIG_DST1
  } else if (tx_type == FLIPADST_DCT) {
    flipud(dest, stride, 8);
    vp9_iht8x8_64_add(input, dest, stride, ADST_DCT);
    flipud(dest, stride, 8);
  } else if (tx_type == DCT_FLIPADST) {
    fliplr(dest, stride, 8);
    vp9_iht8x8_64_add(input, dest, stride, DCT_ADST);
    fliplr(dest, stride, 8);
  } else if (tx_type == FLIPADST_FLIPADST) {
    fliplrud(dest, stride, 8);
    vp9_iht8x8_64_add(input, dest, stride, ADST_ADST);
    fliplrud(dest, stride, 8);
  } else if (tx_type == ADST_FLIPADST) {
    fliplr(dest, stride, 8);
    vp9_iht8x8_64_add(input, dest, stride, ADST_ADST);
    fliplr(dest, stride, 8);
  } else if (tx_type == FLIPADST_ADST) {
    flipud(dest, stride, 8);
    vp9_iht8x8_64_add(input, dest, stride, ADST_ADST);
    flipud(dest, stride, 8);
#endif  // CONFIG_EXT_TX
  } else {
    vp9_iht8x8_64_add(input, dest, stride, tx_type);
  }
}

void vp9_iht16x16_add(TX_TYPE tx_type, const tran_low_t *input, uint8_t *dest,
                      int stride, int eob) {
  if (tx_type == DCT_DCT) {
    vp9_idct16x16_add(input, dest, stride, eob);
#if CONFIG_EXT_TX
#if CONFIG_DST1
  } else if (tx_type == DST_DST) {
    vp9_idst16x16_add(input, dest, stride);
#endif  // CONFIG_DST1
  } else if (tx_type == FLIPADST_DCT) {
    flipud(dest, stride, 16);
    vp9_iht16x16_256_add(input, dest, stride, ADST_DCT);
    flipud(dest, stride, 16);
  } else if (tx_type == DCT_FLIPADST) {
    fliplr(dest, stride, 16);
    vp9_iht16x16_256_add(input, dest, stride, DCT_ADST);
    fliplr(dest, stride, 16);
  } else if (tx_type == FLIPADST_FLIPADST) {
    fliplrud(dest, stride, 16);
    vp9_iht16x16_256_add(input, dest, stride, ADST_ADST);
    fliplrud(dest, stride, 16);
  } else if (tx_type == ADST_FLIPADST) {
    fliplr(dest, stride, 16);
    vp9_iht16x16_256_add(input, dest, stride, ADST_ADST);
    fliplr(dest, stride, 16);
  } else if (tx_type == FLIPADST_ADST) {
    flipud(dest, stride, 16);
    vp9_iht16x16_256_add(input, dest, stride, ADST_ADST);
    flipud(dest, stride, 16);
#endif  // CONFIG_EXT_TX
  } else {
    vp9_iht16x16_256_add(input, dest, stride, tx_type);
  }
}

#if CONFIG_TX_SKIP
void vp9_tx_identity_add_rect(const tran_low_t *input, uint8_t *dest,
                              int row, int col,
                              int stride_in, int stride_out, int shift) {
  int r, c, temp;
  for (r = 0; r < row; r++)
    for (c = 0; c < col; c++) {
      temp = dest[r * stride_out + c] + (input[r * stride_in + c] >> shift);
      dest[r * stride_out + c] = clip_pixel(temp);
    }
}

void vp9_tx_identity_add(const tran_low_t *input, uint8_t *dest,
                         int stride, int bs, int shift) {
  vp9_tx_identity_add_rect(input, dest, bs, bs, bs, stride, shift);
}

#if CONFIG_VP9_HIGHBITDEPTH
void vp9_highbd_tx_identity_add_rect(const tran_low_t *input, uint8_t *dest8,
                                     int row, int col, int stride_in,
                                     int stride_out, int shift, int bd) {
  int r, c;
  uint16_t *dest = CONVERT_TO_SHORTPTR(dest8);
  for (r = 0; r < row; r++)
    for (c = 0; c < col; c++) {
      dest[r * stride_out + c] =
          highbd_clip_pixel_add(dest[r * stride_out + c],
                                input[r * stride_in + c] >> shift, bd);
    }
}

void vp9_highbd_tx_identity_add(const tran_low_t *input, uint8_t *dest8,
                                int stride, int bs, int shift, int bd) {
  vp9_highbd_tx_identity_add_rect(input, dest8, bs, bs, bs, stride, shift, bd);
}
#endif  // CONFIG_VP9_HIGHBITDEPTH
#endif  // CONFIG_TX_SKIP

#if CONFIG_TX64X64
#define DownshiftMultiplyBy2(x) x * 2
#define DownshiftMultiply(x) x

static void idct16f(double *input, double *output, int stride) {
  static const double C1 = 0.995184726672197;
  static const double C2 = 0.98078528040323;
  static const double C3 = 0.956940335732209;
  static const double C4 = 0.923879532511287;
  static const double C5 = 0.881921264348355;
  static const double C6 = 0.831469612302545;
  static const double C7 = 0.773010453362737;
  static const double C8 = 0.707106781186548;
  static const double C9 = 0.634393284163646;
  static const double C10 = 0.555570233019602;
  static const double C11 = 0.471396736825998;
  static const double C12 = 0.38268343236509;
  static const double C13 = 0.290284677254462;
  static const double C14 = 0.195090322016128;
  static const double C15 = 0.098017140329561;

  double step[16];
  double intermediate[16];
  double temp1, temp2;

  // step 1 and 2
  step[ 0] = input[stride*0] + input[stride*8];
  step[ 1] = input[stride*0] - input[stride*8];

  temp1 = input[stride*4]*C12;
  temp2 = input[stride*12]*C4;

  temp1 -= temp2;
  temp1 = DownshiftMultiply(temp1);
  temp1 *= C8;

  step[ 2] = DownshiftMultiplyBy2(temp1);

  temp1 = input[stride*4]*C4;
  temp2 = input[stride*12]*C12;
  temp1 += temp2;
  temp1 = DownshiftMultiply(temp1);
  temp1 *= C8;
  step[ 3] = DownshiftMultiplyBy2(temp1);

  temp1 = input[stride*2]*C8;
  temp1 = DownshiftMultiplyBy2(temp1);
  temp2 = input[stride*6] + input[stride*10];

  step[ 4] = temp1 + temp2;
  step[ 5] = temp1 - temp2;

  temp1 = input[stride*14]*C8;
  temp1 = DownshiftMultiplyBy2(temp1);
  temp2 = input[stride*6] - input[stride*10];

  step[ 6] = temp2 - temp1;
  step[ 7] = temp2 + temp1;

  // for odd input
  temp1 = input[stride*3]*C12;
  temp2 = input[stride*13]*C4;
  temp1 += temp2;
  temp1 = DownshiftMultiply(temp1);
  temp1 *= C8;
  intermediate[ 8] = DownshiftMultiplyBy2(temp1);

  temp1 = input[stride*3]*C4;
  temp2 = input[stride*13]*C12;
  temp2 -= temp1;
  temp2 = DownshiftMultiply(temp2);
  temp2 *= C8;
  intermediate[ 9] = DownshiftMultiplyBy2(temp2);

  intermediate[10] = DownshiftMultiplyBy2(input[stride*9]*C8);
  intermediate[11] = input[stride*15] - input[stride*1];
  intermediate[12] = input[stride*15] + input[stride*1];
  intermediate[13] = DownshiftMultiplyBy2((input[stride*7]*C8));

  temp1 = input[stride*11]*C12;
  temp2 = input[stride*5]*C4;
  temp2 -= temp1;
  temp2 = DownshiftMultiply(temp2);
  temp2 *= C8;
  intermediate[14] = DownshiftMultiplyBy2(temp2);

  temp1 = input[stride*11]*C4;
  temp2 = input[stride*5]*C12;
  temp1 += temp2;
  temp1 = DownshiftMultiply(temp1);
  temp1 *= C8;
  intermediate[15] = DownshiftMultiplyBy2(temp1);

  step[ 8] = intermediate[ 8] + intermediate[14];
  step[ 9] = intermediate[ 9] + intermediate[15];
  step[10] = intermediate[10] + intermediate[11];
  step[11] = intermediate[10] - intermediate[11];
  step[12] = intermediate[12] + intermediate[13];
  step[13] = intermediate[12] - intermediate[13];
  step[14] = intermediate[ 8] - intermediate[14];
  step[15] = intermediate[ 9] - intermediate[15];

  // step 3
  output[stride*0] = step[ 0] + step[ 3];
  output[stride*1] = step[ 1] + step[ 2];
  output[stride*2] = step[ 1] - step[ 2];
  output[stride*3] = step[ 0] - step[ 3];

  temp1 = step[ 4]*C14;
  temp2 = step[ 7]*C2;
  temp1 -= temp2;
  output[stride*4] =  DownshiftMultiply(temp1);

  temp1 = step[ 4]*C2;
  temp2 = step[ 7]*C14;
  temp1 += temp2;
  output[stride*7] =  DownshiftMultiply(temp1);

  temp1 = step[ 5]*C10;
  temp2 = step[ 6]*C6;
  temp1 -= temp2;
  output[stride*5] =  DownshiftMultiply(temp1);

  temp1 = step[ 5]*C6;
  temp2 = step[ 6]*C10;
  temp1 += temp2;
  output[stride*6] =  DownshiftMultiply(temp1);

  output[stride*8] = step[ 8] + step[11];
  output[stride*9] = step[ 9] + step[10];
  output[stride*10] = step[ 9] - step[10];
  output[stride*11] = step[ 8] - step[11];
  output[stride*12] = step[12] + step[15];
  output[stride*13] = step[13] + step[14];
  output[stride*14] = step[13] - step[14];
  output[stride*15] = step[12] - step[15];

  // output 4
  step[ 0] = output[stride*0] + output[stride*7];
  step[ 1] = output[stride*1] + output[stride*6];
  step[ 2] = output[stride*2] + output[stride*5];
  step[ 3] = output[stride*3] + output[stride*4];
  step[ 4] = output[stride*3] - output[stride*4];
  step[ 5] = output[stride*2] - output[stride*5];
  step[ 6] = output[stride*1] - output[stride*6];
  step[ 7] = output[stride*0] - output[stride*7];

  temp1 = output[stride*8]*C7;
  temp2 = output[stride*15]*C9;
  temp1 -= temp2;
  step[ 8] = DownshiftMultiply(temp1);

  temp1 = output[stride*9]*C11;
  temp2 = output[stride*14]*C5;
  temp1 += temp2;
  step[ 9] = DownshiftMultiply(temp1);

  temp1 = output[stride*10]*C3;
  temp2 = output[stride*13]*C13;
  temp1 -= temp2;
  step[10] = DownshiftMultiply(temp1);

  temp1 = output[stride*11]*C15;
  temp2 = output[stride*12]*C1;
  temp1 += temp2;
  step[11] = DownshiftMultiply(temp1);

  temp1 = output[stride*11]*C1;
  temp2 = output[stride*12]*C15;
  temp2 -= temp1;
  step[12] = DownshiftMultiply(temp2);

  temp1 = output[stride*10]*C13;
  temp2 = output[stride*13]*C3;
  temp1 += temp2;
  step[13] = DownshiftMultiply(temp1);

  temp1 = output[stride*9]*C5;
  temp2 = output[stride*14]*C11;
  temp2 -= temp1;
  step[14] = DownshiftMultiply(temp2);

  temp1 = output[stride*8]*C9;
  temp2 = output[stride*15]*C7;
  temp1 += temp2;
  step[15] = DownshiftMultiply(temp1);

  // step 5
  output[stride*0] = step[0] + step[15];
  output[stride*1] = step[1] + step[14];
  output[stride*2] = step[2] + step[13];
  output[stride*3] = step[3] + step[12];
  output[stride*4] = step[4] + step[11];
  output[stride*5] = step[5] + step[10];
  output[stride*6] = step[6] + step[ 9];
  output[stride*7] = step[7] + step[ 8];

  output[stride*15] = step[0] - step[15];
  output[stride*14] = step[1] - step[14];
  output[stride*13] = step[2] - step[13];
  output[stride*12] = step[3] - step[12];
  output[stride*11] = step[4] - step[11];
  output[stride*10] = step[5] - step[10];
  output[stride*9] = step[6] - step[ 9];
  output[stride*8] = step[7] - step[ 8];
}

static void butterfly_32_idct_1d(double *input, double *output, int stride) {
  static const double C1 = 0.998795456205;  // cos(pi * 1 / 64)
  static const double C3 = 0.989176509965;  // cos(pi * 3 / 64)
  static const double C5 = 0.970031253195;  // cos(pi * 5 / 64)
  static const double C7 = 0.941544065183;  // cos(pi * 7 / 64)
  static const double C9 = 0.903989293123;  // cos(pi * 9 / 64)
  static const double C11 = 0.857728610000;  // cos(pi * 11 / 64)
  static const double C13 = 0.803207531481;  // cos(pi * 13 / 64)
  static const double C15 = 0.740951125355;  // cos(pi * 15 / 64)
  static const double C16 = 0.707106781187;  // cos(pi * 16 / 64)
  static const double C17 = 0.671558954847;  // cos(pi * 17 / 64)
  static const double C19 = 0.595699304492;  // cos(pi * 19 / 64)
  static const double C21 = 0.514102744193;  // cos(pi * 21 / 64)
  static const double C23 = 0.427555093430;  // cos(pi * 23 / 64)
  static const double C25 = 0.336889853392;  // cos(pi * 25 / 64)
  static const double C27 = 0.242980179903;  // cos(pi * 27 / 64)
  static const double C29 = 0.146730474455;  // cos(pi * 29 / 64)
  static const double C31 = 0.049067674327;  // cos(pi * 31 / 64)

  double step1[32];
  double step2[32];

  step1[ 0] = input[stride*0];
  step1[ 1] = input[stride*2];
  step1[ 2] = input[stride*4];
  step1[ 3] = input[stride*6];
  step1[ 4] = input[stride*8];
  step1[ 5] = input[stride*10];
  step1[ 6] = input[stride*12];
  step1[ 7] = input[stride*14];
  step1[ 8] = input[stride*16];
  step1[ 9] = input[stride*18];
  step1[10] = input[stride*20];
  step1[11] = input[stride*22];
  step1[12] = input[stride*24];
  step1[13] = input[stride*26];
  step1[14] = input[stride*28];
  step1[15] = input[stride*30];

  step1[16] = DownshiftMultiplyBy2(input[stride*1]*C16);
  step1[17] = (input[stride*3] + input[stride*1]);
  step1[18] = (input[stride*5] + input[stride*3]);
  step1[19] = (input[stride*7] + input[stride*5]);
  step1[20] = (input[stride*9] + input[stride*7]);
  step1[21] = (input[stride*11] + input[stride*9]);
  step1[22] = (input[stride*13] + input[stride*11]);
  step1[23] = (input[stride*15] + input[stride*13]);
  step1[24] = (input[stride*17] + input[stride*15]);
  step1[25] = (input[stride*19] + input[stride*17]);
  step1[26] = (input[stride*21] + input[stride*19]);
  step1[27] = (input[stride*23] + input[stride*21]);
  step1[28] = (input[stride*25] + input[stride*23]);
  step1[29] = (input[stride*27] + input[stride*25]);
  step1[30] = (input[stride*29] + input[stride*27]);
  step1[31] = (input[stride*31] + input[stride*29]);

  idct16f(step1, step2, 1);
  idct16f(step1 + 16, step2 + 16, 1);

  step2[16] = DownshiftMultiply(step2[16] / (2*C1));
  step2[17] = DownshiftMultiply(step2[17] / (2*C3));
  step2[18] = DownshiftMultiply(step2[18] / (2*C5));
  step2[19] = DownshiftMultiply(step2[19] / (2*C7));
  step2[20] = DownshiftMultiply(step2[20] / (2*C9));
  step2[21] = DownshiftMultiply(step2[21] / (2*C11));
  step2[22] = DownshiftMultiply(step2[22] / (2*C13));
  step2[23] = DownshiftMultiply(step2[23] / (2*C15));
  step2[24] = DownshiftMultiply(step2[24] / (2*C17));
  step2[25] = DownshiftMultiply(step2[25] / (2*C19));
  step2[26] = DownshiftMultiply(step2[26] / (2*C21));
  step2[27] = DownshiftMultiply(step2[27] / (2*C23));
  step2[28] = DownshiftMultiply(step2[28] / (2*C25));
  step2[29] = DownshiftMultiply(step2[29] / (2*C27));
  step2[30] = DownshiftMultiply(step2[30] / (2*C29));
  step2[31] = DownshiftMultiply(step2[31] / (2*C31));

  output[stride* 0] = step2[ 0] + step2[16];
  output[stride* 1] = step2[ 1] + step2[17];
  output[stride* 2] = step2[ 2] + step2[18];
  output[stride* 3] = step2[ 3] + step2[19];
  output[stride* 4] = step2[ 4] + step2[20];
  output[stride* 5] = step2[ 5] + step2[21];
  output[stride* 6] = step2[ 6] + step2[22];
  output[stride* 7] = step2[ 7] + step2[23];
  output[stride* 8] = step2[ 8] + step2[24];
  output[stride* 9] = step2[ 9] + step2[25];
  output[stride*10] = step2[10] + step2[26];
  output[stride*11] = step2[11] + step2[27];
  output[stride*12] = step2[12] + step2[28];
  output[stride*13] = step2[13] + step2[29];
  output[stride*14] = step2[14] + step2[30];
  output[stride*15] = step2[15] + step2[31];
  output[stride*16] = step2[15] - step2[(31 - 0)];
  output[stride*17] = step2[14] - step2[(31 - 1)];
  output[stride*18] = step2[13] - step2[(31 - 2)];
  output[stride*19] = step2[12] - step2[(31 - 3)];
  output[stride*20] = step2[11] - step2[(31 - 4)];
  output[stride*21] = step2[10] - step2[(31 - 5)];
  output[stride*22] = step2[ 9] - step2[(31 - 6)];
  output[stride*23] = step2[ 8] - step2[(31 - 7)];
  output[stride*24] = step2[ 7] - step2[(31 - 8)];
  output[stride*25] = step2[ 6] - step2[(31 - 9)];
  output[stride*26] = step2[ 5] - step2[(31 - 10)];
  output[stride*27] = step2[ 4] - step2[(31 - 11)];
  output[stride*28] = step2[ 3] - step2[(31 - 12)];
  output[stride*29] = step2[ 2] - step2[(31 - 13)];
  output[stride*30] = step2[ 1] - step2[(31 - 14)];
  output[stride*31] = step2[ 0] - step2[(31 - 15)];
}

static void butterfly_64_idct_1d(double *input, double *output, int stride) {
  double step1[64], step2[64];
  int i;
  static const double C[64] = {
    1.00000000000000000000,  // cos(0 * pi / 128)
    0.99969881869620424997,  // cos(1 * pi / 128)
    0.99879545620517240501,  // cos(2 * pi / 128)
    0.99729045667869020697,  // cos(3 * pi / 128)
    0.99518472667219692873,  // cos(4 * pi / 128)
    0.99247953459870996706,  // cos(5 * pi / 128)
    0.98917650996478101444,  // cos(6 * pi / 128)
    0.98527764238894122162,  // cos(7 * pi / 128)
    0.98078528040323043058,  // cos(8 * pi / 128)
    0.97570213003852857003,  // cos(9 * pi / 128)
    0.97003125319454397424,  // cos(10 * pi / 128)
    0.96377606579543984022,  // cos(11 * pi / 128)
    0.95694033573220882438,  // cos(12 * pi / 128)
    0.94952818059303667475,  // cos(13 * pi / 128)
    0.94154406518302080631,  // cos(14 * pi / 128)
    0.93299279883473895669,  // cos(15 * pi / 128)
    0.92387953251128673848,  // cos(16 * pi / 128)
    0.91420975570353069095,  // cos(17 * pi / 128)
    0.90398929312344333820,  // cos(18 * pi / 128)
    0.89322430119551532446,  // cos(19 * pi / 128)
    0.88192126434835504956,  // cos(20 * pi / 128)
    0.87008699110871146054,  // cos(21 * pi / 128)
    0.85772861000027211809,  // cos(22 * pi / 128)
    0.84485356524970711689,  // cos(23 * pi / 128)
    0.83146961230254523567,  // cos(24 * pi / 128)
    0.81758481315158371139,  // cos(25 * pi / 128)
    0.80320753148064494287,  // cos(26 * pi / 128)
    0.78834642762660633863,  // cos(27 * pi / 128)
    0.77301045336273699338,  // cos(28 * pi / 128)
    0.75720884650648456748,  // cos(29 * pi / 128)
    0.74095112535495921691,  // cos(30 * pi / 128)
    0.72424708295146700276,  // cos(31 * pi / 128)
    0.70710678118654757274,  // cos(32 * pi / 128)
    0.68954054473706694051,  // cos(33 * pi / 128)
    0.67155895484701844111,  // cos(34 * pi / 128)
    0.65317284295377686654,  // cos(35 * pi / 128)
    0.63439328416364559882,  // cos(36 * pi / 128)
    0.61523159058062693028,  // cos(37 * pi / 128)
    0.59569930449243346793,  // cos(38 * pi / 128)
    0.57580819141784544968,  // cos(39 * pi / 128)
    0.55557023301960228867,  // cos(40 * pi / 128)
    0.53499761988709737537,  // cos(41 * pi / 128)
    0.51410274419322177231,  // cos(42 * pi / 128)
    0.49289819222978414892,  // cos(43 * pi / 128)
    0.47139673682599780857,  // cos(44 * pi / 128)
    0.44961132965460659516,  // cos(45 * pi / 128)
    0.42755509343028219593,  // cos(46 * pi / 128)
    0.40524131400498980549,  // cos(47 * pi / 128)
    0.38268343236508983729,  // cos(48 * pi / 128)
    0.35989503653498827740,  // cos(49 * pi / 128)
    0.33688985339222005111,  // cos(50 * pi / 128)
    0.31368174039889151761,  // cos(51 * pi / 128)
    0.29028467725446227554,  // cos(52 * pi / 128)
    0.26671275747489842090,  // cos(53 * pi / 128)
    0.24298017990326398197,  // cos(54 * pi / 128)
    0.21910124015686976984,  // cos(55 * pi / 128)
    0.19509032201612830359,  // cos(56 * pi / 128)
    0.17096188876030135595,  // cos(57 * pi / 128)
    0.14673047445536174793,  // cos(58 * pi / 128)
    0.12241067519921627893,  // cos(59 * pi / 128)
    0.09801714032956077016,  // cos(60 * pi / 128)
    0.07356456359966745406,  // cos(61 * pi / 128)
    0.04906767432741813290,  // cos(62 * pi / 128)
    0.02454122852291226731,  // cos(63 * pi / 128)
  };

  for (i = 0; i < 64; i += 2) {
    step1[i / 2] = input[stride * i];
  }
  step1[32] = DownshiftMultiplyBy2(input[stride*1] * C[32]);
  for (i = 3; i < 64; i+=2) {
    step1[32 + i/2] = (input[stride * i] + input[stride * (i - 2)]);
  }

  butterfly_32_idct_1d(step1, step2, 1);
  butterfly_32_idct_1d(step1 + 32, step2 + 32, 1);

  for (i = 32; i < 64; ++i) {
    step2[i] = DownshiftMultiply(step2[i] / (2 * C[(i - 32) * 2 + 1]));
  }

  for (i = 0; i < 32; ++i) {
    output[stride * i] = step2[i] + step2[32 + i];
  }

  for (i = 0; i < 32; ++i) {
    output[stride * (i + 32)] = step2[31 - i] - step2[63 - i];
  }
}

void vp9_idct64x64_4096_add_c(const tran_low_t *input, uint8_t *dest,
                              int stride) {
  // vp9_clear_system_state();  // Make it simd safe : __asm emms;
  {
    double out[64 * 64], out2[64 * 64];
    int i, j;
    // First transform rows
    for (i = 0; i < 64; ++i) {
      double temp_in[64], temp_out[64];
      for (j = 0; j < 64; ++j)
        temp_in[j] = input[j + i * 64];
      butterfly_64_idct_1d(temp_in, temp_out, 1);
      for (j = 0; j < 64; ++j)
        out[j + i * 64] = temp_out[j];
    }
    // Then transform columns
    for (i = 0; i < 64; ++i) {
      double temp_in[64], temp_out[64];
      for (j = 0; j < 64; ++j)
        temp_in[j] = out[j * 64 + i];
      butterfly_64_idct_1d(temp_in, temp_out, 1);
      for (j = 0; j < 64; ++j)
        out2[j * 64 + i] = temp_out[j];
    }

    for (j = 0; j < 64; ++j) {
      for (i = 0; i < 64; ++i)
        dest[i] = clip_pixel_add(dest[i], round(out2[j * 64 + i] / 128));
      dest += stride;
    }
  }
  // vp9_clear_system_state();  // Make it simd safe : __asm emms;
}

void vp9_idct64x64_add(const tran_low_t *input, uint8_t *dest,
                       int stride, int eob) {
  (void) eob;
  vp9_idct64x64_4096_add_c(input, dest, stride);
}
#endif

#if CONFIG_VP9_HIGHBITDEPTH
void vp9_highbd_iwht4x4_16_add_c(const tran_low_t *input, uint8_t *dest8,
                                 int stride, int bd) {
  /* 4-point reversible, orthonormal inverse Walsh-Hadamard in 3.5 adds,
     0.5 shifts per pixel. */
  int i;
  tran_low_t output[16];
  tran_high_t a1, b1, c1, d1, e1;
  const tran_low_t *ip = input;
  tran_low_t *op = output;
  uint16_t *dest = CONVERT_TO_SHORTPTR(dest8);

  for (i = 0; i < 4; i++) {
    a1 = ip[0] >> UNIT_QUANT_SHIFT;
    c1 = ip[1] >> UNIT_QUANT_SHIFT;
    d1 = ip[2] >> UNIT_QUANT_SHIFT;
    b1 = ip[3] >> UNIT_QUANT_SHIFT;
    a1 += c1;
    d1 -= b1;
    e1 = (a1 - d1) >> 1;
    b1 = e1 - b1;
    c1 = e1 - c1;
    a1 -= b1;
    d1 += c1;
    op[0] = WRAPLOW(a1, bd);
    op[1] = WRAPLOW(b1, bd);
    op[2] = WRAPLOW(c1, bd);
    op[3] = WRAPLOW(d1, bd);
    ip += 4;
    op += 4;
  }

  ip = output;
  for (i = 0; i < 4; i++) {
    a1 = ip[4 * 0];
    c1 = ip[4 * 1];
    d1 = ip[4 * 2];
    b1 = ip[4 * 3];
    a1 += c1;
    d1 -= b1;
    e1 = (a1 - d1) >> 1;
    b1 = e1 - b1;
    c1 = e1 - c1;
    a1 -= b1;
    d1 += c1;
    dest[stride * 0] = highbd_clip_pixel_add(dest[stride * 0], a1, bd);
    dest[stride * 1] = highbd_clip_pixel_add(dest[stride * 1], b1, bd);
    dest[stride * 2] = highbd_clip_pixel_add(dest[stride * 2], c1, bd);
    dest[stride * 3] = highbd_clip_pixel_add(dest[stride * 3], d1, bd);

    ip++;
    dest++;
  }
}

void vp9_highbd_iwht4x4_1_add_c(const tran_low_t *in, uint8_t *dest8,
                                int dest_stride, int bd) {
  int i;
  tran_high_t a1, e1;
  tran_low_t tmp[4];
  const tran_low_t *ip = in;
  tran_low_t *op = tmp;
  uint16_t *dest = CONVERT_TO_SHORTPTR(dest8);
  (void) bd;

  a1 = ip[0] >> UNIT_QUANT_SHIFT;
  e1 = a1 >> 1;
  a1 -= e1;
  op[0] = WRAPLOW(a1, bd);
  op[1] = op[2] = op[3] = WRAPLOW(e1, bd);

  ip = tmp;
  for (i = 0; i < 4; i++) {
    e1 = ip[0] >> 1;
    a1 = ip[0] - e1;
    dest[dest_stride * 0] = highbd_clip_pixel_add(
        dest[dest_stride * 0], a1, bd);
    dest[dest_stride * 1] = highbd_clip_pixel_add(
        dest[dest_stride * 1], e1, bd);
    dest[dest_stride * 2] = highbd_clip_pixel_add(
        dest[dest_stride * 2], e1, bd);
    dest[dest_stride * 3] = highbd_clip_pixel_add(
        dest[dest_stride * 3], e1, bd);
    ip++;
    dest++;
  }
}

void vp9_highbd_idct4(const tran_low_t *input, tran_low_t *output, int bd) {
  tran_low_t step[4];
  tran_high_t temp1, temp2;
  (void) bd;
  // stage 1
  temp1 = (input[0] + input[2]) * cospi_16_64;
  temp2 = (input[0] - input[2]) * cospi_16_64;
  step[0] = WRAPLOW(dct_const_round_shift(temp1), bd);
  step[1] = WRAPLOW(dct_const_round_shift(temp2), bd);
  temp1 = input[1] * cospi_24_64 - input[3] * cospi_8_64;
  temp2 = input[1] * cospi_8_64 + input[3] * cospi_24_64;
  step[2] = WRAPLOW(dct_const_round_shift(temp1), bd);
  step[3] = WRAPLOW(dct_const_round_shift(temp2), bd);

  // stage 2
  output[0] = WRAPLOW(step[0] + step[3], bd);
  output[1] = WRAPLOW(step[1] + step[2], bd);
  output[2] = WRAPLOW(step[1] - step[2], bd);
  output[3] = WRAPLOW(step[0] - step[3], bd);
}

void vp9_highbd_idct4x4_16_add_c(const tran_low_t *input, uint8_t *dest8,
                                 int stride, int bd) {
  tran_low_t out[4 * 4];
  tran_low_t *outptr = out;
  int i, j;
  tran_low_t temp_in[4], temp_out[4];
  uint16_t *dest = CONVERT_TO_SHORTPTR(dest8);

  // Rows
  for (i = 0; i < 4; ++i) {
    vp9_highbd_idct4(input, outptr, bd);
    input += 4;
    outptr += 4;
  }

  // Columns
  for (i = 0; i < 4; ++i) {
    for (j = 0; j < 4; ++j)
      temp_in[j] = out[j * 4 + i];
    vp9_highbd_idct4(temp_in, temp_out, bd);
    for (j = 0; j < 4; ++j) {
      dest[j * stride + i] = highbd_clip_pixel_add(
          dest[j * stride + i], ROUND_POWER_OF_TWO(temp_out[j], 4), bd);
    }
  }
}

void vp9_highbd_idct4x4_1_add_c(const tran_low_t *input, uint8_t *dest8,
                                int dest_stride, int bd) {
  int i;
  tran_high_t a1;
  tran_low_t out = WRAPLOW(dct_const_round_shift(input[0] * cospi_16_64), bd);
  uint16_t *dest = CONVERT_TO_SHORTPTR(dest8);

  out = WRAPLOW(dct_const_round_shift(out * cospi_16_64), bd);
  a1 = ROUND_POWER_OF_TWO(out, 4);

  for (i = 0; i < 4; i++) {
    dest[0] = highbd_clip_pixel_add(dest[0], a1, bd);
    dest[1] = highbd_clip_pixel_add(dest[1], a1, bd);
    dest[2] = highbd_clip_pixel_add(dest[2], a1, bd);
    dest[3] = highbd_clip_pixel_add(dest[3], a1, bd);
    dest += dest_stride;
  }
}

void vp9_highbd_idct8(const tran_low_t *input, tran_low_t *output, int bd) {
  tran_low_t step1[8], step2[8];
  tran_high_t temp1, temp2;
  // stage 1
  step1[0] = input[0];
  step1[2] = input[4];
  step1[1] = input[2];
  step1[3] = input[6];
  temp1 = input[1] * cospi_28_64 - input[7] * cospi_4_64;
  temp2 = input[1] * cospi_4_64 + input[7] * cospi_28_64;
  step1[4] = WRAPLOW(dct_const_round_shift(temp1), bd);
  step1[7] = WRAPLOW(dct_const_round_shift(temp2), bd);
  temp1 = input[5] * cospi_12_64 - input[3] * cospi_20_64;
  temp2 = input[5] * cospi_20_64 + input[3] * cospi_12_64;
  step1[5] = WRAPLOW(dct_const_round_shift(temp1), bd);
  step1[6] = WRAPLOW(dct_const_round_shift(temp2), bd);

  // stage 2 & stage 3 - even half
  vp9_highbd_idct4(step1, step1, bd);

  // stage 2 - odd half
  step2[4] = WRAPLOW(step1[4] + step1[5], bd);
  step2[5] = WRAPLOW(step1[4] - step1[5], bd);
  step2[6] = WRAPLOW(-step1[6] + step1[7], bd);
  step2[7] = WRAPLOW(step1[6] + step1[7], bd);

  // stage 3 - odd half
  step1[4] = step2[4];
  temp1 = (step2[6] - step2[5]) * cospi_16_64;
  temp2 = (step2[5] + step2[6]) * cospi_16_64;
  step1[5] = WRAPLOW(dct_const_round_shift(temp1), bd);
  step1[6] = WRAPLOW(dct_const_round_shift(temp2), bd);
  step1[7] = step2[7];

  // stage 4
  output[0] = WRAPLOW(step1[0] + step1[7], bd);
  output[1] = WRAPLOW(step1[1] + step1[6], bd);
  output[2] = WRAPLOW(step1[2] + step1[5], bd);
  output[3] = WRAPLOW(step1[3] + step1[4], bd);
  output[4] = WRAPLOW(step1[3] - step1[4], bd);
  output[5] = WRAPLOW(step1[2] - step1[5], bd);
  output[6] = WRAPLOW(step1[1] - step1[6], bd);
  output[7] = WRAPLOW(step1[0] - step1[7], bd);
}

void vp9_highbd_idct8x8_64_add_c(const tran_low_t *input, uint8_t *dest8,
                                 int stride, int bd) {
  tran_low_t out[8 * 8];
  tran_low_t *outptr = out;
  int i, j;
  tran_low_t temp_in[8], temp_out[8];
  uint16_t *dest = CONVERT_TO_SHORTPTR(dest8);

  // First transform rows.
  for (i = 0; i < 8; ++i) {
    vp9_highbd_idct8(input, outptr, bd);
    input += 8;
    outptr += 8;
  }

  // Then transform columns.
  for (i = 0; i < 8; ++i) {
    for (j = 0; j < 8; ++j)
      temp_in[j] = out[j * 8 + i];
    vp9_highbd_idct8(temp_in, temp_out, bd);
    for (j = 0; j < 8; ++j) {
      dest[j * stride + i] = highbd_clip_pixel_add(
          dest[j * stride + i], ROUND_POWER_OF_TWO(temp_out[j], 5), bd);
    }
  }
}

void vp9_highbd_idct8x8_1_add_c(const tran_low_t *input, uint8_t *dest8,
                                int stride, int bd) {
  int i, j;
  tran_high_t a1;
  tran_low_t out = WRAPLOW(dct_const_round_shift(input[0] * cospi_16_64), bd);
  uint16_t *dest = CONVERT_TO_SHORTPTR(dest8);
  out = WRAPLOW(dct_const_round_shift(out * cospi_16_64), bd);
  a1 = ROUND_POWER_OF_TWO(out, 5);
  for (j = 0; j < 8; ++j) {
    for (i = 0; i < 8; ++i)
      dest[i] = highbd_clip_pixel_add(dest[i], a1, bd);
    dest += stride;
  }
}

static void highbd_iadst4(const tran_low_t *input, tran_low_t *output, int bd) {
  tran_high_t s0, s1, s2, s3, s4, s5, s6, s7;

  tran_low_t x0 = input[0];
  tran_low_t x1 = input[1];
  tran_low_t x2 = input[2];
  tran_low_t x3 = input[3];
  (void) bd;

  if (!(x0 | x1 | x2 | x3)) {
    vpx_memset(output, 0, 4 * sizeof(*output));
    return;
  }

  s0 = sinpi_1_9 * x0;
  s1 = sinpi_2_9 * x0;
  s2 = sinpi_3_9 * x1;
  s3 = sinpi_4_9 * x2;
  s4 = sinpi_1_9 * x2;
  s5 = sinpi_2_9 * x3;
  s6 = sinpi_4_9 * x3;
  s7 = (tran_high_t)(x0 - x2 + x3);

  s0 = s0 + s3 + s5;
  s1 = s1 - s4 - s6;
  s3 = s2;
  s2 = sinpi_3_9 * s7;

  // 1-D transform scaling factor is sqrt(2).
  // The overall dynamic range is 14b (input) + 14b (multiplication scaling)
  // + 1b (addition) = 29b.
  // Hence the output bit depth is 15b.
  output[0] = WRAPLOW(dct_const_round_shift(s0 + s3), bd);
  output[1] = WRAPLOW(dct_const_round_shift(s1 + s3), bd);
  output[2] = WRAPLOW(dct_const_round_shift(s2), bd);
  output[3] = WRAPLOW(dct_const_round_shift(s0 + s1 - s3), bd);
}

void vp9_highbd_iht4x4_16_add_c(const tran_low_t *input, uint8_t *dest8,
                                int stride, int tx_type, int bd) {
  const highbd_transform_2d IHT_4[] = {
    { vp9_highbd_idct4, vp9_highbd_idct4  },    // DCT_DCT  = 0
    { highbd_iadst4, vp9_highbd_idct4 },    // ADST_DCT = 1
    { vp9_highbd_idct4, highbd_iadst4 },    // DCT_ADST = 2
    { highbd_iadst4, highbd_iadst4 }    // ADST_ADST = 3
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

static void highbd_iadst8(const tran_low_t *input, tran_low_t *output, int bd) {
  tran_high_t s0, s1, s2, s3, s4, s5, s6, s7;

  tran_low_t x0 = input[7];
  tran_low_t x1 = input[0];
  tran_low_t x2 = input[5];
  tran_low_t x3 = input[2];
  tran_low_t x4 = input[3];
  tran_low_t x5 = input[4];
  tran_low_t x6 = input[1];
  tran_low_t x7 = input[6];
  (void) bd;

  if (!(x0 | x1 | x2 | x3 | x4 | x5 | x6 | x7)) {
    vpx_memset(output, 0, 8 * sizeof(*output));
    return;
  }

  // stage 1
  s0 = cospi_2_64  * x0 + cospi_30_64 * x1;
  s1 = cospi_30_64 * x0 - cospi_2_64  * x1;
  s2 = cospi_10_64 * x2 + cospi_22_64 * x3;
  s3 = cospi_22_64 * x2 - cospi_10_64 * x3;
  s4 = cospi_18_64 * x4 + cospi_14_64 * x5;
  s5 = cospi_14_64 * x4 - cospi_18_64 * x5;
  s6 = cospi_26_64 * x6 + cospi_6_64  * x7;
  s7 = cospi_6_64  * x6 - cospi_26_64 * x7;

  x0 = WRAPLOW(dct_const_round_shift(s0 + s4), bd);
  x1 = WRAPLOW(dct_const_round_shift(s1 + s5), bd);
  x2 = WRAPLOW(dct_const_round_shift(s2 + s6), bd);
  x3 = WRAPLOW(dct_const_round_shift(s3 + s7), bd);
  x4 = WRAPLOW(dct_const_round_shift(s0 - s4), bd);
  x5 = WRAPLOW(dct_const_round_shift(s1 - s5), bd);
  x6 = WRAPLOW(dct_const_round_shift(s2 - s6), bd);
  x7 = WRAPLOW(dct_const_round_shift(s3 - s7), bd);

  // stage 2
  s0 = x0;
  s1 = x1;
  s2 = x2;
  s3 = x3;
  s4 =  cospi_8_64  * x4 + cospi_24_64 * x5;
  s5 =  cospi_24_64 * x4 - cospi_8_64  * x5;
  s6 = -cospi_24_64 * x6 + cospi_8_64  * x7;
  s7 =  cospi_8_64  * x6 + cospi_24_64 * x7;

  x0 = WRAPLOW(s0 + s2, bd);
  x1 = WRAPLOW(s1 + s3, bd);
  x2 = WRAPLOW(s0 - s2, bd);
  x3 = WRAPLOW(s1 - s3, bd);
  x4 = WRAPLOW(dct_const_round_shift(s4 + s6), bd);
  x5 = WRAPLOW(dct_const_round_shift(s5 + s7), bd);
  x6 = WRAPLOW(dct_const_round_shift(s4 - s6), bd);
  x7 = WRAPLOW(dct_const_round_shift(s5 - s7), bd);

  // stage 3
  s2 = cospi_16_64 * (x2 + x3);
  s3 = cospi_16_64 * (x2 - x3);
  s6 = cospi_16_64 * (x6 + x7);
  s7 = cospi_16_64 * (x6 - x7);

  x2 = WRAPLOW(dct_const_round_shift(s2), bd);
  x3 = WRAPLOW(dct_const_round_shift(s3), bd);
  x6 = WRAPLOW(dct_const_round_shift(s6), bd);
  x7 = WRAPLOW(dct_const_round_shift(s7), bd);

  output[0] = WRAPLOW(x0, bd);
  output[1] = WRAPLOW(-x4, bd);
  output[2] = WRAPLOW(x6, bd);
  output[3] = WRAPLOW(-x2, bd);
  output[4] = WRAPLOW(x3, bd);
  output[5] = WRAPLOW(-x7, bd);
  output[6] = WRAPLOW(x5, bd);
  output[7] = WRAPLOW(-x1, bd);
}

static const highbd_transform_2d HIGH_IHT_8[] = {
  { vp9_highbd_idct8,  vp9_highbd_idct8  },  // DCT_DCT  = 0
  { highbd_iadst8, vp9_highbd_idct8  },  // ADST_DCT = 1
  { vp9_highbd_idct8,  highbd_iadst8 },  // DCT_ADST = 2
  { highbd_iadst8, highbd_iadst8 }   // ADST_ADST = 3
};

void vp9_highbd_iht8x8_64_add_c(const tran_low_t *input, uint8_t *dest8,
                                int stride, int tx_type, int bd) {
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

void vp9_highbd_idct8x8_10_add_c(const tran_low_t *input, uint8_t *dest8,
                                 int stride, int bd) {
  tran_low_t out[8 * 8] = { 0 };
  tran_low_t *outptr = out;
  int i, j;
  tran_low_t temp_in[8], temp_out[8];
  uint16_t *dest = CONVERT_TO_SHORTPTR(dest8);

  // First transform rows.
  // Only first 4 row has non-zero coefs.
  for (i = 0; i < 4; ++i) {
    vp9_highbd_idct8(input, outptr, bd);
    input += 8;
    outptr += 8;
  }
  // Then transform columns.
  for (i = 0; i < 8; ++i) {
    for (j = 0; j < 8; ++j)
      temp_in[j] = out[j * 8 + i];
    vp9_highbd_idct8(temp_in, temp_out, bd);
    for (j = 0; j < 8; ++j) {
      dest[j * stride + i] = highbd_clip_pixel_add(
          dest[j * stride + i], ROUND_POWER_OF_TWO(temp_out[j], 5), bd);
    }
  }
}

void vp9_highbd_idct16(const tran_low_t *input, tran_low_t *output, int bd) {
  tran_low_t step1[16], step2[16];
  tran_high_t temp1, temp2;
  (void) bd;

  // stage 1
  step1[0] = input[0/2];
  step1[1] = input[16/2];
  step1[2] = input[8/2];
  step1[3] = input[24/2];
  step1[4] = input[4/2];
  step1[5] = input[20/2];
  step1[6] = input[12/2];
  step1[7] = input[28/2];
  step1[8] = input[2/2];
  step1[9] = input[18/2];
  step1[10] = input[10/2];
  step1[11] = input[26/2];
  step1[12] = input[6/2];
  step1[13] = input[22/2];
  step1[14] = input[14/2];
  step1[15] = input[30/2];

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
  step2[8] = WRAPLOW(dct_const_round_shift(temp1), bd);
  step2[15] = WRAPLOW(dct_const_round_shift(temp2), bd);

  temp1 = step1[9] * cospi_14_64 - step1[14] * cospi_18_64;
  temp2 = step1[9] * cospi_18_64 + step1[14] * cospi_14_64;
  step2[9] = WRAPLOW(dct_const_round_shift(temp1), bd);
  step2[14] = WRAPLOW(dct_const_round_shift(temp2), bd);

  temp1 = step1[10] * cospi_22_64 - step1[13] * cospi_10_64;
  temp2 = step1[10] * cospi_10_64 + step1[13] * cospi_22_64;
  step2[10] = WRAPLOW(dct_const_round_shift(temp1), bd);
  step2[13] = WRAPLOW(dct_const_round_shift(temp2), bd);

  temp1 = step1[11] * cospi_6_64 - step1[12] * cospi_26_64;
  temp2 = step1[11] * cospi_26_64 + step1[12] * cospi_6_64;
  step2[11] = WRAPLOW(dct_const_round_shift(temp1), bd);
  step2[12] = WRAPLOW(dct_const_round_shift(temp2), bd);

  // stage 3
  step1[0] = step2[0];
  step1[1] = step2[1];
  step1[2] = step2[2];
  step1[3] = step2[3];

  temp1 = step2[4] * cospi_28_64 - step2[7] * cospi_4_64;
  temp2 = step2[4] * cospi_4_64 + step2[7] * cospi_28_64;
  step1[4] = WRAPLOW(dct_const_round_shift(temp1), bd);
  step1[7] = WRAPLOW(dct_const_round_shift(temp2), bd);
  temp1 = step2[5] * cospi_12_64 - step2[6] * cospi_20_64;
  temp2 = step2[5] * cospi_20_64 + step2[6] * cospi_12_64;
  step1[5] = WRAPLOW(dct_const_round_shift(temp1), bd);
  step1[6] = WRAPLOW(dct_const_round_shift(temp2), bd);

  step1[8] = WRAPLOW(step2[8] + step2[9], bd);
  step1[9] = WRAPLOW(step2[8] - step2[9], bd);
  step1[10] = WRAPLOW(-step2[10] + step2[11], bd);
  step1[11] = WRAPLOW(step2[10] + step2[11], bd);
  step1[12] = WRAPLOW(step2[12] + step2[13], bd);
  step1[13] = WRAPLOW(step2[12] - step2[13], bd);
  step1[14] = WRAPLOW(-step2[14] + step2[15], bd);
  step1[15] = WRAPLOW(step2[14] + step2[15], bd);

  // stage 4
  temp1 = (step1[0] + step1[1]) * cospi_16_64;
  temp2 = (step1[0] - step1[1]) * cospi_16_64;
  step2[0] = WRAPLOW(dct_const_round_shift(temp1), bd);
  step2[1] = WRAPLOW(dct_const_round_shift(temp2), bd);
  temp1 = step1[2] * cospi_24_64 - step1[3] * cospi_8_64;
  temp2 = step1[2] * cospi_8_64 + step1[3] * cospi_24_64;
  step2[2] = WRAPLOW(dct_const_round_shift(temp1), bd);
  step2[3] = WRAPLOW(dct_const_round_shift(temp2), bd);
  step2[4] = WRAPLOW(step1[4] + step1[5], bd);
  step2[5] = WRAPLOW(step1[4] - step1[5], bd);
  step2[6] = WRAPLOW(-step1[6] + step1[7], bd);
  step2[7] = WRAPLOW(step1[6] + step1[7], bd);

  step2[8] = step1[8];
  step2[15] = step1[15];
  temp1 = -step1[9] * cospi_8_64 + step1[14] * cospi_24_64;
  temp2 = step1[9] * cospi_24_64 + step1[14] * cospi_8_64;
  step2[9] = WRAPLOW(dct_const_round_shift(temp1), bd);
  step2[14] = WRAPLOW(dct_const_round_shift(temp2), bd);
  temp1 = -step1[10] * cospi_24_64 - step1[13] * cospi_8_64;
  temp2 = -step1[10] * cospi_8_64 + step1[13] * cospi_24_64;
  step2[10] = WRAPLOW(dct_const_round_shift(temp1), bd);
  step2[13] = WRAPLOW(dct_const_round_shift(temp2), bd);
  step2[11] = step1[11];
  step2[12] = step1[12];

  // stage 5
  step1[0] = WRAPLOW(step2[0] + step2[3], bd);
  step1[1] = WRAPLOW(step2[1] + step2[2], bd);
  step1[2] = WRAPLOW(step2[1] - step2[2], bd);
  step1[3] = WRAPLOW(step2[0] - step2[3], bd);
  step1[4] = step2[4];
  temp1 = (step2[6] - step2[5]) * cospi_16_64;
  temp2 = (step2[5] + step2[6]) * cospi_16_64;
  step1[5] = WRAPLOW(dct_const_round_shift(temp1), bd);
  step1[6] = WRAPLOW(dct_const_round_shift(temp2), bd);
  step1[7] = step2[7];

  step1[8] = WRAPLOW(step2[8] + step2[11], bd);
  step1[9] = WRAPLOW(step2[9] + step2[10], bd);
  step1[10] = WRAPLOW(step2[9] - step2[10], bd);
  step1[11] = WRAPLOW(step2[8] - step2[11], bd);
  step1[12] = WRAPLOW(-step2[12] + step2[15], bd);
  step1[13] = WRAPLOW(-step2[13] + step2[14], bd);
  step1[14] = WRAPLOW(step2[13] + step2[14], bd);
  step1[15] = WRAPLOW(step2[12] + step2[15], bd);

  // stage 6
  step2[0] = WRAPLOW(step1[0] + step1[7], bd);
  step2[1] = WRAPLOW(step1[1] + step1[6], bd);
  step2[2] = WRAPLOW(step1[2] + step1[5], bd);
  step2[3] = WRAPLOW(step1[3] + step1[4], bd);
  step2[4] = WRAPLOW(step1[3] - step1[4], bd);
  step2[5] = WRAPLOW(step1[2] - step1[5], bd);
  step2[6] = WRAPLOW(step1[1] - step1[6], bd);
  step2[7] = WRAPLOW(step1[0] - step1[7], bd);
  step2[8] = step1[8];
  step2[9] = step1[9];
  temp1 = (-step1[10] + step1[13]) * cospi_16_64;
  temp2 = (step1[10] + step1[13]) * cospi_16_64;
  step2[10] = WRAPLOW(dct_const_round_shift(temp1), bd);
  step2[13] = WRAPLOW(dct_const_round_shift(temp2), bd);
  temp1 = (-step1[11] + step1[12]) * cospi_16_64;
  temp2 = (step1[11] + step1[12]) * cospi_16_64;
  step2[11] = WRAPLOW(dct_const_round_shift(temp1), bd);
  step2[12] = WRAPLOW(dct_const_round_shift(temp2), bd);
  step2[14] = step1[14];
  step2[15] = step1[15];

  // stage 7
  output[0] = WRAPLOW(step2[0] + step2[15], bd);
  output[1] = WRAPLOW(step2[1] + step2[14], bd);
  output[2] = WRAPLOW(step2[2] + step2[13], bd);
  output[3] = WRAPLOW(step2[3] + step2[12], bd);
  output[4] = WRAPLOW(step2[4] + step2[11], bd);
  output[5] = WRAPLOW(step2[5] + step2[10], bd);
  output[6] = WRAPLOW(step2[6] + step2[9], bd);
  output[7] = WRAPLOW(step2[7] + step2[8], bd);
  output[8] = WRAPLOW(step2[7] - step2[8], bd);
  output[9] = WRAPLOW(step2[6] - step2[9], bd);
  output[10] = WRAPLOW(step2[5] - step2[10], bd);
  output[11] = WRAPLOW(step2[4] - step2[11], bd);
  output[12] = WRAPLOW(step2[3] - step2[12], bd);
  output[13] = WRAPLOW(step2[2] - step2[13], bd);
  output[14] = WRAPLOW(step2[1] - step2[14], bd);
  output[15] = WRAPLOW(step2[0] - step2[15], bd);
}

void vp9_highbd_idct16x16_256_add_c(const tran_low_t *input, uint8_t *dest8,
                                    int stride, int bd) {
  tran_low_t out[16 * 16];
  tran_low_t *outptr = out;
  int i, j;
  tran_low_t temp_in[16], temp_out[16];
  uint16_t *dest = CONVERT_TO_SHORTPTR(dest8);

  // First transform rows.
  for (i = 0; i < 16; ++i) {
    vp9_highbd_idct16(input, outptr, bd);
    input += 16;
    outptr += 16;
  }

  // Then transform columns.
  for (i = 0; i < 16; ++i) {
    for (j = 0; j < 16; ++j)
      temp_in[j] = out[j * 16 + i];
    vp9_highbd_idct16(temp_in, temp_out, bd);
    for (j = 0; j < 16; ++j) {
      dest[j * stride + i] = highbd_clip_pixel_add(
          dest[j * stride + i], ROUND_POWER_OF_TWO(temp_out[j], 6), bd);
    }
  }
}

static void highbd_iadst16(const tran_low_t *input, tran_low_t *output,
                           int bd) {
  tran_high_t s0, s1, s2, s3, s4, s5, s6, s7, s8;
  tran_high_t s9, s10, s11, s12, s13, s14, s15;

  tran_low_t x0 = input[15];
  tran_low_t x1 = input[0];
  tran_low_t x2 = input[13];
  tran_low_t x3 = input[2];
  tran_low_t x4 = input[11];
  tran_low_t x5 = input[4];
  tran_low_t x6 = input[9];
  tran_low_t x7 = input[6];
  tran_low_t x8 = input[7];
  tran_low_t x9 = input[8];
  tran_low_t x10 = input[5];
  tran_low_t x11 = input[10];
  tran_low_t x12 = input[3];
  tran_low_t x13 = input[12];
  tran_low_t x14 = input[1];
  tran_low_t x15 = input[14];
  (void) bd;

  if (!(x0 | x1 | x2 | x3 | x4 | x5 | x6 | x7 | x8
           | x9 | x10 | x11 | x12 | x13 | x14 | x15)) {
    vpx_memset(output, 0, 16 * sizeof(*output));
    return;
  }

  // stage 1
  s0 = x0 * cospi_1_64  + x1 * cospi_31_64;
  s1 = x0 * cospi_31_64 - x1 * cospi_1_64;
  s2 = x2 * cospi_5_64  + x3 * cospi_27_64;
  s3 = x2 * cospi_27_64 - x3 * cospi_5_64;
  s4 = x4 * cospi_9_64  + x5 * cospi_23_64;
  s5 = x4 * cospi_23_64 - x5 * cospi_9_64;
  s6 = x6 * cospi_13_64 + x7 * cospi_19_64;
  s7 = x6 * cospi_19_64 - x7 * cospi_13_64;
  s8 = x8 * cospi_17_64 + x9 * cospi_15_64;
  s9 = x8 * cospi_15_64 - x9 * cospi_17_64;
  s10 = x10 * cospi_21_64 + x11 * cospi_11_64;
  s11 = x10 * cospi_11_64 - x11 * cospi_21_64;
  s12 = x12 * cospi_25_64 + x13 * cospi_7_64;
  s13 = x12 * cospi_7_64  - x13 * cospi_25_64;
  s14 = x14 * cospi_29_64 + x15 * cospi_3_64;
  s15 = x14 * cospi_3_64  - x15 * cospi_29_64;

  x0 = WRAPLOW(dct_const_round_shift(s0 + s8), bd);
  x1 = WRAPLOW(dct_const_round_shift(s1 + s9), bd);
  x2 = WRAPLOW(dct_const_round_shift(s2 + s10), bd);
  x3 = WRAPLOW(dct_const_round_shift(s3 + s11), bd);
  x4 = WRAPLOW(dct_const_round_shift(s4 + s12), bd);
  x5 = WRAPLOW(dct_const_round_shift(s5 + s13), bd);
  x6 = WRAPLOW(dct_const_round_shift(s6 + s14), bd);
  x7 = WRAPLOW(dct_const_round_shift(s7 + s15), bd);
  x8  = WRAPLOW(dct_const_round_shift(s0 - s8), bd);
  x9  = WRAPLOW(dct_const_round_shift(s1 - s9), bd);
  x10 = WRAPLOW(dct_const_round_shift(s2 - s10), bd);
  x11 = WRAPLOW(dct_const_round_shift(s3 - s11), bd);
  x12 = WRAPLOW(dct_const_round_shift(s4 - s12), bd);
  x13 = WRAPLOW(dct_const_round_shift(s5 - s13), bd);
  x14 = WRAPLOW(dct_const_round_shift(s6 - s14), bd);
  x15 = WRAPLOW(dct_const_round_shift(s7 - s15), bd);

  // stage 2
  s0 = x0;
  s1 = x1;
  s2 = x2;
  s3 = x3;
  s4 = x4;
  s5 = x5;
  s6 = x6;
  s7 = x7;
  s8 = x8 * cospi_4_64 + x9 * cospi_28_64;
  s9 = x8 * cospi_28_64 - x9 * cospi_4_64;
  s10 = x10 * cospi_20_64 + x11 * cospi_12_64;
  s11 = x10 * cospi_12_64 - x11 * cospi_20_64;
  s12 = -x12 * cospi_28_64 + x13 * cospi_4_64;
  s13 = x12 * cospi_4_64 + x13 * cospi_28_64;
  s14 = -x14 * cospi_12_64 + x15 * cospi_20_64;
  s15 = x14 * cospi_20_64 + x15 * cospi_12_64;

  x0 = WRAPLOW(s0 + s4, bd);
  x1 = WRAPLOW(s1 + s5, bd);
  x2 = WRAPLOW(s2 + s6, bd);
  x3 = WRAPLOW(s3 + s7, bd);
  x4 = WRAPLOW(s0 - s4, bd);
  x5 = WRAPLOW(s1 - s5, bd);
  x6 = WRAPLOW(s2 - s6, bd);
  x7 = WRAPLOW(s3 - s7, bd);
  x8 = WRAPLOW(dct_const_round_shift(s8 + s12), bd);
  x9 = WRAPLOW(dct_const_round_shift(s9 + s13), bd);
  x10 = WRAPLOW(dct_const_round_shift(s10 + s14), bd);
  x11 = WRAPLOW(dct_const_round_shift(s11 + s15), bd);
  x12 = WRAPLOW(dct_const_round_shift(s8 - s12), bd);
  x13 = WRAPLOW(dct_const_round_shift(s9 - s13), bd);
  x14 = WRAPLOW(dct_const_round_shift(s10 - s14), bd);
  x15 = WRAPLOW(dct_const_round_shift(s11 - s15), bd);

  // stage 3
  s0 = x0;
  s1 = x1;
  s2 = x2;
  s3 = x3;
  s4 = x4 * cospi_8_64 + x5 * cospi_24_64;
  s5 = x4 * cospi_24_64 - x5 * cospi_8_64;
  s6 = -x6 * cospi_24_64 + x7 * cospi_8_64;
  s7 = x6 * cospi_8_64 + x7 * cospi_24_64;
  s8 = x8;
  s9 = x9;
  s10 = x10;
  s11 = x11;
  s12 = x12 * cospi_8_64 + x13 * cospi_24_64;
  s13 = x12 * cospi_24_64 - x13 * cospi_8_64;
  s14 = -x14 * cospi_24_64 + x15 * cospi_8_64;
  s15 = x14 * cospi_8_64 + x15 * cospi_24_64;

  x0 = WRAPLOW(s0 + s2, bd);
  x1 = WRAPLOW(s1 + s3, bd);
  x2 = WRAPLOW(s0 - s2, bd);
  x3 = WRAPLOW(s1 - s3, bd);
  x4 = WRAPLOW(dct_const_round_shift(s4 + s6), bd);
  x5 = WRAPLOW(dct_const_round_shift(s5 + s7), bd);
  x6 = WRAPLOW(dct_const_round_shift(s4 - s6), bd);
  x7 = WRAPLOW(dct_const_round_shift(s5 - s7), bd);
  x8 = WRAPLOW(s8 + s10, bd);
  x9 = WRAPLOW(s9 + s11, bd);
  x10 = WRAPLOW(s8 - s10, bd);
  x11 = WRAPLOW(s9 - s11, bd);
  x12 = WRAPLOW(dct_const_round_shift(s12 + s14), bd);
  x13 = WRAPLOW(dct_const_round_shift(s13 + s15), bd);
  x14 = WRAPLOW(dct_const_round_shift(s12 - s14), bd);
  x15 = WRAPLOW(dct_const_round_shift(s13 - s15), bd);

  // stage 4
  s2 = (- cospi_16_64) * (x2 + x3);
  s3 = cospi_16_64 * (x2 - x3);
  s6 = cospi_16_64 * (x6 + x7);
  s7 = cospi_16_64 * (-x6 + x7);
  s10 = cospi_16_64 * (x10 + x11);
  s11 = cospi_16_64 * (-x10 + x11);
  s14 = (- cospi_16_64) * (x14 + x15);
  s15 = cospi_16_64 * (x14 - x15);

  x2 = WRAPLOW(dct_const_round_shift(s2), bd);
  x3 = WRAPLOW(dct_const_round_shift(s3), bd);
  x6 = WRAPLOW(dct_const_round_shift(s6), bd);
  x7 = WRAPLOW(dct_const_round_shift(s7), bd);
  x10 = WRAPLOW(dct_const_round_shift(s10), bd);
  x11 = WRAPLOW(dct_const_round_shift(s11), bd);
  x14 = WRAPLOW(dct_const_round_shift(s14), bd);
  x15 = WRAPLOW(dct_const_round_shift(s15), bd);

  output[0] = WRAPLOW(x0, bd);
  output[1] = WRAPLOW(-x8, bd);
  output[2] = WRAPLOW(x12, bd);
  output[3] = WRAPLOW(-x4, bd);
  output[4] = WRAPLOW(x6, bd);
  output[5] = WRAPLOW(x14, bd);
  output[6] = WRAPLOW(x10, bd);
  output[7] = WRAPLOW(x2, bd);
  output[8] = WRAPLOW(x3, bd);
  output[9] = WRAPLOW(x11, bd);
  output[10] = WRAPLOW(x15, bd);
  output[11] = WRAPLOW(x7, bd);
  output[12] = WRAPLOW(x5, bd);
  output[13] = WRAPLOW(-x13, bd);
  output[14] = WRAPLOW(x9, bd);
  output[15] = WRAPLOW(-x1, bd);
}

static const highbd_transform_2d HIGH_IHT_16[] = {
  { vp9_highbd_idct16,  vp9_highbd_idct16  },  // DCT_DCT  = 0
  { highbd_iadst16, vp9_highbd_idct16  },  // ADST_DCT = 1
  { vp9_highbd_idct16,  highbd_iadst16 },  // DCT_ADST = 2
  { highbd_iadst16, highbd_iadst16 }   // ADST_ADST = 3
};

void vp9_highbd_iht16x16_256_add_c(const tran_low_t *input, uint8_t *dest8,
                                   int stride, int tx_type, int bd) {
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

void vp9_highbd_idct16x16_10_add_c(const tran_low_t *input, uint8_t *dest8,
                                   int stride, int bd) {
  tran_low_t out[16 * 16] = { 0 };
  tran_low_t *outptr = out;
  int i, j;
  tran_low_t temp_in[16], temp_out[16];
  uint16_t *dest = CONVERT_TO_SHORTPTR(dest8);

  // First transform rows. Since all non-zero dct coefficients are in
  // upper-left 4x4 area, we only need to calculate first 4 rows here.
  for (i = 0; i < 4; ++i) {
    vp9_highbd_idct16(input, outptr, bd);
    input += 16;
    outptr += 16;
  }

  // Then transform columns.
  for (i = 0; i < 16; ++i) {
    for (j = 0; j < 16; ++j)
      temp_in[j] = out[j*16 + i];
    vp9_highbd_idct16(temp_in, temp_out, bd);
    for (j = 0; j < 16; ++j) {
      dest[j * stride + i] = highbd_clip_pixel_add(
          dest[j * stride + i], ROUND_POWER_OF_TWO(temp_out[j], 6), bd);
    }
  }
}

void vp9_highbd_idct16x16_1_add_c(const tran_low_t *input, uint8_t *dest8,
                                  int stride, int bd) {
  int i, j;
  tran_high_t a1;
  tran_low_t out = WRAPLOW(dct_const_round_shift(input[0] * cospi_16_64), bd);
  uint16_t *dest = CONVERT_TO_SHORTPTR(dest8);

  out = WRAPLOW(dct_const_round_shift(out * cospi_16_64), bd);
  a1 = ROUND_POWER_OF_TWO(out, 6);
  for (j = 0; j < 16; ++j) {
    for (i = 0; i < 16; ++i)
      dest[i] = highbd_clip_pixel_add(dest[i], a1, bd);
    dest += stride;
  }
}

static void highbd_idct32(const tran_low_t *input, tran_low_t *output, int bd) {
  tran_low_t step1[32], step2[32];
  tran_high_t temp1, temp2;
  (void) bd;

  // stage 1
  step1[0] = input[0];
  step1[1] = input[16];
  step1[2] = input[8];
  step1[3] = input[24];
  step1[4] = input[4];
  step1[5] = input[20];
  step1[6] = input[12];
  step1[7] = input[28];
  step1[8] = input[2];
  step1[9] = input[18];
  step1[10] = input[10];
  step1[11] = input[26];
  step1[12] = input[6];
  step1[13] = input[22];
  step1[14] = input[14];
  step1[15] = input[30];

  temp1 = input[1] * cospi_31_64 - input[31] * cospi_1_64;
  temp2 = input[1] * cospi_1_64 + input[31] * cospi_31_64;
  step1[16] = WRAPLOW(dct_const_round_shift(temp1), bd);
  step1[31] = WRAPLOW(dct_const_round_shift(temp2), bd);

  temp1 = input[17] * cospi_15_64 - input[15] * cospi_17_64;
  temp2 = input[17] * cospi_17_64 + input[15] * cospi_15_64;
  step1[17] = WRAPLOW(dct_const_round_shift(temp1), bd);
  step1[30] = WRAPLOW(dct_const_round_shift(temp2), bd);

  temp1 = input[9] * cospi_23_64 - input[23] * cospi_9_64;
  temp2 = input[9] * cospi_9_64 + input[23] * cospi_23_64;
  step1[18] = WRAPLOW(dct_const_round_shift(temp1), bd);
  step1[29] = WRAPLOW(dct_const_round_shift(temp2), bd);

  temp1 = input[25] * cospi_7_64 - input[7] * cospi_25_64;
  temp2 = input[25] * cospi_25_64 + input[7] * cospi_7_64;
  step1[19] = WRAPLOW(dct_const_round_shift(temp1), bd);
  step1[28] = WRAPLOW(dct_const_round_shift(temp2), bd);

  temp1 = input[5] * cospi_27_64 - input[27] * cospi_5_64;
  temp2 = input[5] * cospi_5_64 + input[27] * cospi_27_64;
  step1[20] = WRAPLOW(dct_const_round_shift(temp1), bd);
  step1[27] = WRAPLOW(dct_const_round_shift(temp2), bd);

  temp1 = input[21] * cospi_11_64 - input[11] * cospi_21_64;
  temp2 = input[21] * cospi_21_64 + input[11] * cospi_11_64;
  step1[21] = WRAPLOW(dct_const_round_shift(temp1), bd);
  step1[26] = WRAPLOW(dct_const_round_shift(temp2), bd);

  temp1 = input[13] * cospi_19_64 - input[19] * cospi_13_64;
  temp2 = input[13] * cospi_13_64 + input[19] * cospi_19_64;
  step1[22] = WRAPLOW(dct_const_round_shift(temp1), bd);
  step1[25] = WRAPLOW(dct_const_round_shift(temp2), bd);

  temp1 = input[29] * cospi_3_64 - input[3] * cospi_29_64;
  temp2 = input[29] * cospi_29_64 + input[3] * cospi_3_64;
  step1[23] = WRAPLOW(dct_const_round_shift(temp1), bd);
  step1[24] = WRAPLOW(dct_const_round_shift(temp2), bd);

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
  step2[8] = WRAPLOW(dct_const_round_shift(temp1), bd);
  step2[15] = WRAPLOW(dct_const_round_shift(temp2), bd);

  temp1 = step1[9] * cospi_14_64 - step1[14] * cospi_18_64;
  temp2 = step1[9] * cospi_18_64 + step1[14] * cospi_14_64;
  step2[9] = WRAPLOW(dct_const_round_shift(temp1), bd);
  step2[14] = WRAPLOW(dct_const_round_shift(temp2), bd);

  temp1 = step1[10] * cospi_22_64 - step1[13] * cospi_10_64;
  temp2 = step1[10] * cospi_10_64 + step1[13] * cospi_22_64;
  step2[10] = WRAPLOW(dct_const_round_shift(temp1), bd);
  step2[13] = WRAPLOW(dct_const_round_shift(temp2), bd);

  temp1 = step1[11] * cospi_6_64 - step1[12] * cospi_26_64;
  temp2 = step1[11] * cospi_26_64 + step1[12] * cospi_6_64;
  step2[11] = WRAPLOW(dct_const_round_shift(temp1), bd);
  step2[12] = WRAPLOW(dct_const_round_shift(temp2), bd);

  step2[16] = WRAPLOW(step1[16] + step1[17], bd);
  step2[17] = WRAPLOW(step1[16] - step1[17], bd);
  step2[18] = WRAPLOW(-step1[18] + step1[19], bd);
  step2[19] = WRAPLOW(step1[18] + step1[19], bd);
  step2[20] = WRAPLOW(step1[20] + step1[21], bd);
  step2[21] = WRAPLOW(step1[20] - step1[21], bd);
  step2[22] = WRAPLOW(-step1[22] + step1[23], bd);
  step2[23] = WRAPLOW(step1[22] + step1[23], bd);
  step2[24] = WRAPLOW(step1[24] + step1[25], bd);
  step2[25] = WRAPLOW(step1[24] - step1[25], bd);
  step2[26] = WRAPLOW(-step1[26] + step1[27], bd);
  step2[27] = WRAPLOW(step1[26] + step1[27], bd);
  step2[28] = WRAPLOW(step1[28] + step1[29], bd);
  step2[29] = WRAPLOW(step1[28] - step1[29], bd);
  step2[30] = WRAPLOW(-step1[30] + step1[31], bd);
  step2[31] = WRAPLOW(step1[30] + step1[31], bd);

  // stage 3
  step1[0] = step2[0];
  step1[1] = step2[1];
  step1[2] = step2[2];
  step1[3] = step2[3];

  temp1 = step2[4] * cospi_28_64 - step2[7] * cospi_4_64;
  temp2 = step2[4] * cospi_4_64 + step2[7] * cospi_28_64;
  step1[4] = WRAPLOW(dct_const_round_shift(temp1), bd);
  step1[7] = WRAPLOW(dct_const_round_shift(temp2), bd);
  temp1 = step2[5] * cospi_12_64 - step2[6] * cospi_20_64;
  temp2 = step2[5] * cospi_20_64 + step2[6] * cospi_12_64;
  step1[5] = WRAPLOW(dct_const_round_shift(temp1), bd);
  step1[6] = WRAPLOW(dct_const_round_shift(temp2), bd);

  step1[8] = WRAPLOW(step2[8] + step2[9], bd);
  step1[9] = WRAPLOW(step2[8] - step2[9], bd);
  step1[10] = WRAPLOW(-step2[10] + step2[11], bd);
  step1[11] = WRAPLOW(step2[10] + step2[11], bd);
  step1[12] = WRAPLOW(step2[12] + step2[13], bd);
  step1[13] = WRAPLOW(step2[12] - step2[13], bd);
  step1[14] = WRAPLOW(-step2[14] + step2[15], bd);
  step1[15] = WRAPLOW(step2[14] + step2[15], bd);

  step1[16] = step2[16];
  step1[31] = step2[31];
  temp1 = -step2[17] * cospi_4_64 + step2[30] * cospi_28_64;
  temp2 = step2[17] * cospi_28_64 + step2[30] * cospi_4_64;
  step1[17] = WRAPLOW(dct_const_round_shift(temp1), bd);
  step1[30] = WRAPLOW(dct_const_round_shift(temp2), bd);
  temp1 = -step2[18] * cospi_28_64 - step2[29] * cospi_4_64;
  temp2 = -step2[18] * cospi_4_64 + step2[29] * cospi_28_64;
  step1[18] = WRAPLOW(dct_const_round_shift(temp1), bd);
  step1[29] = WRAPLOW(dct_const_round_shift(temp2), bd);
  step1[19] = step2[19];
  step1[20] = step2[20];
  temp1 = -step2[21] * cospi_20_64 + step2[26] * cospi_12_64;
  temp2 = step2[21] * cospi_12_64 + step2[26] * cospi_20_64;
  step1[21] = WRAPLOW(dct_const_round_shift(temp1), bd);
  step1[26] = WRAPLOW(dct_const_round_shift(temp2), bd);
  temp1 = -step2[22] * cospi_12_64 - step2[25] * cospi_20_64;
  temp2 = -step2[22] * cospi_20_64 + step2[25] * cospi_12_64;
  step1[22] = WRAPLOW(dct_const_round_shift(temp1), bd);
  step1[25] = WRAPLOW(dct_const_round_shift(temp2), bd);
  step1[23] = step2[23];
  step1[24] = step2[24];
  step1[27] = step2[27];
  step1[28] = step2[28];

  // stage 4
  temp1 = (step1[0] + step1[1]) * cospi_16_64;
  temp2 = (step1[0] - step1[1]) * cospi_16_64;
  step2[0] = WRAPLOW(dct_const_round_shift(temp1), bd);
  step2[1] = WRAPLOW(dct_const_round_shift(temp2), bd);
  temp1 = step1[2] * cospi_24_64 - step1[3] * cospi_8_64;
  temp2 = step1[2] * cospi_8_64 + step1[3] * cospi_24_64;
  step2[2] = WRAPLOW(dct_const_round_shift(temp1), bd);
  step2[3] = WRAPLOW(dct_const_round_shift(temp2), bd);
  step2[4] = WRAPLOW(step1[4] + step1[5], bd);
  step2[5] = WRAPLOW(step1[4] - step1[5], bd);
  step2[6] = WRAPLOW(-step1[6] + step1[7], bd);
  step2[7] = WRAPLOW(step1[6] + step1[7], bd);

  step2[8] = step1[8];
  step2[15] = step1[15];
  temp1 = -step1[9] * cospi_8_64 + step1[14] * cospi_24_64;
  temp2 = step1[9] * cospi_24_64 + step1[14] * cospi_8_64;
  step2[9] = WRAPLOW(dct_const_round_shift(temp1), bd);
  step2[14] = WRAPLOW(dct_const_round_shift(temp2), bd);
  temp1 = -step1[10] * cospi_24_64 - step1[13] * cospi_8_64;
  temp2 = -step1[10] * cospi_8_64 + step1[13] * cospi_24_64;
  step2[10] = WRAPLOW(dct_const_round_shift(temp1), bd);
  step2[13] = WRAPLOW(dct_const_round_shift(temp2), bd);
  step2[11] = step1[11];
  step2[12] = step1[12];

  step2[16] = WRAPLOW(step1[16] + step1[19], bd);
  step2[17] = WRAPLOW(step1[17] + step1[18], bd);
  step2[18] = WRAPLOW(step1[17] - step1[18], bd);
  step2[19] = WRAPLOW(step1[16] - step1[19], bd);
  step2[20] = WRAPLOW(-step1[20] + step1[23], bd);
  step2[21] = WRAPLOW(-step1[21] + step1[22], bd);
  step2[22] = WRAPLOW(step1[21] + step1[22], bd);
  step2[23] = WRAPLOW(step1[20] + step1[23], bd);

  step2[24] = WRAPLOW(step1[24] + step1[27], bd);
  step2[25] = WRAPLOW(step1[25] + step1[26], bd);
  step2[26] = WRAPLOW(step1[25] - step1[26], bd);
  step2[27] = WRAPLOW(step1[24] - step1[27], bd);
  step2[28] = WRAPLOW(-step1[28] + step1[31], bd);
  step2[29] = WRAPLOW(-step1[29] + step1[30], bd);
  step2[30] = WRAPLOW(step1[29] + step1[30], bd);
  step2[31] = WRAPLOW(step1[28] + step1[31], bd);

  // stage 5
  step1[0] = WRAPLOW(step2[0] + step2[3], bd);
  step1[1] = WRAPLOW(step2[1] + step2[2], bd);
  step1[2] = WRAPLOW(step2[1] - step2[2], bd);
  step1[3] = WRAPLOW(step2[0] - step2[3], bd);
  step1[4] = step2[4];
  temp1 = (step2[6] - step2[5]) * cospi_16_64;
  temp2 = (step2[5] + step2[6]) * cospi_16_64;
  step1[5] = WRAPLOW(dct_const_round_shift(temp1), bd);
  step1[6] = WRAPLOW(dct_const_round_shift(temp2), bd);
  step1[7] = step2[7];

  step1[8] = WRAPLOW(step2[8] + step2[11], bd);
  step1[9] = WRAPLOW(step2[9] + step2[10], bd);
  step1[10] = WRAPLOW(step2[9] - step2[10], bd);
  step1[11] = WRAPLOW(step2[8] - step2[11], bd);
  step1[12] = WRAPLOW(-step2[12] + step2[15], bd);
  step1[13] = WRAPLOW(-step2[13] + step2[14], bd);
  step1[14] = WRAPLOW(step2[13] + step2[14], bd);
  step1[15] = WRAPLOW(step2[12] + step2[15], bd);

  step1[16] = step2[16];
  step1[17] = step2[17];
  temp1 = -step2[18] * cospi_8_64 + step2[29] * cospi_24_64;
  temp2 = step2[18] * cospi_24_64 + step2[29] * cospi_8_64;
  step1[18] = WRAPLOW(dct_const_round_shift(temp1), bd);
  step1[29] = WRAPLOW(dct_const_round_shift(temp2), bd);
  temp1 = -step2[19] * cospi_8_64 + step2[28] * cospi_24_64;
  temp2 = step2[19] * cospi_24_64 + step2[28] * cospi_8_64;
  step1[19] = WRAPLOW(dct_const_round_shift(temp1), bd);
  step1[28] = WRAPLOW(dct_const_round_shift(temp2), bd);
  temp1 = -step2[20] * cospi_24_64 - step2[27] * cospi_8_64;
  temp2 = -step2[20] * cospi_8_64 + step2[27] * cospi_24_64;
  step1[20] = WRAPLOW(dct_const_round_shift(temp1), bd);
  step1[27] = WRAPLOW(dct_const_round_shift(temp2), bd);
  temp1 = -step2[21] * cospi_24_64 - step2[26] * cospi_8_64;
  temp2 = -step2[21] * cospi_8_64 + step2[26] * cospi_24_64;
  step1[21] = WRAPLOW(dct_const_round_shift(temp1), bd);
  step1[26] = WRAPLOW(dct_const_round_shift(temp2), bd);
  step1[22] = step2[22];
  step1[23] = step2[23];
  step1[24] = step2[24];
  step1[25] = step2[25];
  step1[30] = step2[30];
  step1[31] = step2[31];

  // stage 6
  step2[0] = WRAPLOW(step1[0] + step1[7], bd);
  step2[1] = WRAPLOW(step1[1] + step1[6], bd);
  step2[2] = WRAPLOW(step1[2] + step1[5], bd);
  step2[3] = WRAPLOW(step1[3] + step1[4], bd);
  step2[4] = WRAPLOW(step1[3] - step1[4], bd);
  step2[5] = WRAPLOW(step1[2] - step1[5], bd);
  step2[6] = WRAPLOW(step1[1] - step1[6], bd);
  step2[7] = WRAPLOW(step1[0] - step1[7], bd);
  step2[8] = step1[8];
  step2[9] = step1[9];
  temp1 = (-step1[10] + step1[13]) * cospi_16_64;
  temp2 = (step1[10] + step1[13]) * cospi_16_64;
  step2[10] = WRAPLOW(dct_const_round_shift(temp1), bd);
  step2[13] = WRAPLOW(dct_const_round_shift(temp2), bd);
  temp1 = (-step1[11] + step1[12]) * cospi_16_64;
  temp2 = (step1[11] + step1[12]) * cospi_16_64;
  step2[11] = WRAPLOW(dct_const_round_shift(temp1), bd);
  step2[12] = WRAPLOW(dct_const_round_shift(temp2), bd);
  step2[14] = step1[14];
  step2[15] = step1[15];

  step2[16] = WRAPLOW(step1[16] + step1[23], bd);
  step2[17] = WRAPLOW(step1[17] + step1[22], bd);
  step2[18] = WRAPLOW(step1[18] + step1[21], bd);
  step2[19] = WRAPLOW(step1[19] + step1[20], bd);
  step2[20] = WRAPLOW(step1[19] - step1[20], bd);
  step2[21] = WRAPLOW(step1[18] - step1[21], bd);
  step2[22] = WRAPLOW(step1[17] - step1[22], bd);
  step2[23] = WRAPLOW(step1[16] - step1[23], bd);

  step2[24] = WRAPLOW(-step1[24] + step1[31], bd);
  step2[25] = WRAPLOW(-step1[25] + step1[30], bd);
  step2[26] = WRAPLOW(-step1[26] + step1[29], bd);
  step2[27] = WRAPLOW(-step1[27] + step1[28], bd);
  step2[28] = WRAPLOW(step1[27] + step1[28], bd);
  step2[29] = WRAPLOW(step1[26] + step1[29], bd);
  step2[30] = WRAPLOW(step1[25] + step1[30], bd);
  step2[31] = WRAPLOW(step1[24] + step1[31], bd);

  // stage 7
  step1[0] = WRAPLOW(step2[0] + step2[15], bd);
  step1[1] = WRAPLOW(step2[1] + step2[14], bd);
  step1[2] = WRAPLOW(step2[2] + step2[13], bd);
  step1[3] = WRAPLOW(step2[3] + step2[12], bd);
  step1[4] = WRAPLOW(step2[4] + step2[11], bd);
  step1[5] = WRAPLOW(step2[5] + step2[10], bd);
  step1[6] = WRAPLOW(step2[6] + step2[9], bd);
  step1[7] = WRAPLOW(step2[7] + step2[8], bd);
  step1[8] = WRAPLOW(step2[7] - step2[8], bd);
  step1[9] = WRAPLOW(step2[6] - step2[9], bd);
  step1[10] = WRAPLOW(step2[5] - step2[10], bd);
  step1[11] = WRAPLOW(step2[4] - step2[11], bd);
  step1[12] = WRAPLOW(step2[3] - step2[12], bd);
  step1[13] = WRAPLOW(step2[2] - step2[13], bd);
  step1[14] = WRAPLOW(step2[1] - step2[14], bd);
  step1[15] = WRAPLOW(step2[0] - step2[15], bd);

  step1[16] = step2[16];
  step1[17] = step2[17];
  step1[18] = step2[18];
  step1[19] = step2[19];
  temp1 = (-step2[20] + step2[27]) * cospi_16_64;
  temp2 = (step2[20] + step2[27]) * cospi_16_64;
  step1[20] = WRAPLOW(dct_const_round_shift(temp1), bd);
  step1[27] = WRAPLOW(dct_const_round_shift(temp2), bd);
  temp1 = (-step2[21] + step2[26]) * cospi_16_64;
  temp2 = (step2[21] + step2[26]) * cospi_16_64;
  step1[21] = WRAPLOW(dct_const_round_shift(temp1), bd);
  step1[26] = WRAPLOW(dct_const_round_shift(temp2), bd);
  temp1 = (-step2[22] + step2[25]) * cospi_16_64;
  temp2 = (step2[22] + step2[25]) * cospi_16_64;
  step1[22] = WRAPLOW(dct_const_round_shift(temp1), bd);
  step1[25] = WRAPLOW(dct_const_round_shift(temp2), bd);
  temp1 = (-step2[23] + step2[24]) * cospi_16_64;
  temp2 = (step2[23] + step2[24]) * cospi_16_64;
  step1[23] = WRAPLOW(dct_const_round_shift(temp1), bd);
  step1[24] = WRAPLOW(dct_const_round_shift(temp2), bd);
  step1[28] = step2[28];
  step1[29] = step2[29];
  step1[30] = step2[30];
  step1[31] = step2[31];

  // final stage
  output[0] = WRAPLOW(step1[0] + step1[31], bd);
  output[1] = WRAPLOW(step1[1] + step1[30], bd);
  output[2] = WRAPLOW(step1[2] + step1[29], bd);
  output[3] = WRAPLOW(step1[3] + step1[28], bd);
  output[4] = WRAPLOW(step1[4] + step1[27], bd);
  output[5] = WRAPLOW(step1[5] + step1[26], bd);
  output[6] = WRAPLOW(step1[6] + step1[25], bd);
  output[7] = WRAPLOW(step1[7] + step1[24], bd);
  output[8] = WRAPLOW(step1[8] + step1[23], bd);
  output[9] = WRAPLOW(step1[9] + step1[22], bd);
  output[10] = WRAPLOW(step1[10] + step1[21], bd);
  output[11] = WRAPLOW(step1[11] + step1[20], bd);
  output[12] = WRAPLOW(step1[12] + step1[19], bd);
  output[13] = WRAPLOW(step1[13] + step1[18], bd);
  output[14] = WRAPLOW(step1[14] + step1[17], bd);
  output[15] = WRAPLOW(step1[15] + step1[16], bd);
  output[16] = WRAPLOW(step1[15] - step1[16], bd);
  output[17] = WRAPLOW(step1[14] - step1[17], bd);
  output[18] = WRAPLOW(step1[13] - step1[18], bd);
  output[19] = WRAPLOW(step1[12] - step1[19], bd);
  output[20] = WRAPLOW(step1[11] - step1[20], bd);
  output[21] = WRAPLOW(step1[10] - step1[21], bd);
  output[22] = WRAPLOW(step1[9] - step1[22], bd);
  output[23] = WRAPLOW(step1[8] - step1[23], bd);
  output[24] = WRAPLOW(step1[7] - step1[24], bd);
  output[25] = WRAPLOW(step1[6] - step1[25], bd);
  output[26] = WRAPLOW(step1[5] - step1[26], bd);
  output[27] = WRAPLOW(step1[4] - step1[27], bd);
  output[28] = WRAPLOW(step1[3] - step1[28], bd);
  output[29] = WRAPLOW(step1[2] - step1[29], bd);
  output[30] = WRAPLOW(step1[1] - step1[30], bd);
  output[31] = WRAPLOW(step1[0] - step1[31], bd);
}

void vp9_highbd_idct32x32_1024_add_c(const tran_low_t *input, uint8_t *dest8,
                                     int stride, int bd) {
  tran_low_t out[32 * 32];
  tran_low_t *outptr = out;
  int i, j;
  tran_low_t temp_in[32], temp_out[32];
  uint16_t *dest = CONVERT_TO_SHORTPTR(dest8);

  // Rows
  for (i = 0; i < 32; ++i) {
    tran_low_t zero_coeff[16];
    for (j = 0; j < 16; ++j)
      zero_coeff[j] = input[2 * j] | input[2 * j + 1];
    for (j = 0; j < 8; ++j)
      zero_coeff[j] = zero_coeff[2 * j] | zero_coeff[2 * j + 1];
    for (j = 0; j < 4; ++j)
      zero_coeff[j] = zero_coeff[2 * j] | zero_coeff[2 * j + 1];
    for (j = 0; j < 2; ++j)
      zero_coeff[j] = zero_coeff[2 * j] | zero_coeff[2 * j + 1];

    if (zero_coeff[0] | zero_coeff[1])
      highbd_idct32(input, outptr, bd);
    else
      vpx_memset(outptr, 0, sizeof(tran_low_t) * 32);
    input += 32;
    outptr += 32;
  }

  // Columns
  for (i = 0; i < 32; ++i) {
    for (j = 0; j < 32; ++j)
      temp_in[j] = out[j * 32 + i];
    highbd_idct32(temp_in, temp_out, bd);
    for (j = 0; j < 32; ++j) {
      dest[j * stride + i] = highbd_clip_pixel_add(
          dest[j * stride + i], ROUND_POWER_OF_TWO(temp_out[j], 6), bd);
    }
  }
}

void vp9_highbd_idct32x32_34_add_c(const tran_low_t *input, uint8_t *dest8,
                                   int stride, int bd) {
  tran_low_t out[32 * 32] = {0};
  tran_low_t *outptr = out;
  int i, j;
  tran_low_t temp_in[32], temp_out[32];
  uint16_t *dest = CONVERT_TO_SHORTPTR(dest8);

  // Rows
  // Only upper-left 8x8 has non-zero coeff.
  for (i = 0; i < 8; ++i) {
    highbd_idct32(input, outptr, bd);
    input += 32;
    outptr += 32;
  }
  // Columns
  for (i = 0; i < 32; ++i) {
    for (j = 0; j < 32; ++j)
      temp_in[j] = out[j * 32 + i];
    highbd_idct32(temp_in, temp_out, bd);
    for (j = 0; j < 32; ++j) {
      dest[j * stride + i] = highbd_clip_pixel_add(
          dest[j * stride + i], ROUND_POWER_OF_TWO(temp_out[j], 6), bd);
    }
  }
}

void vp9_highbd_idct32x32_1_add_c(const tran_low_t *input, uint8_t *dest8,
                                  int stride, int bd) {
  int i, j;
  int a1;
  uint16_t *dest = CONVERT_TO_SHORTPTR(dest8);

  tran_low_t out = WRAPLOW(dct_const_round_shift(input[0] * cospi_16_64), bd);
  out = WRAPLOW(dct_const_round_shift(out * cospi_16_64), bd);
  a1 = ROUND_POWER_OF_TWO(out, 6);

  for (j = 0; j < 32; ++j) {
    for (i = 0; i < 32; ++i)
      dest[i] = highbd_clip_pixel_add(dest[i], a1, bd);
    dest += stride;
  }
}

// idct
void vp9_highbd_idct4x4_add(const tran_low_t *input, uint8_t *dest, int stride,
                            int eob, int bd) {
  if (eob > 1)
    vp9_highbd_idct4x4_16_add(input, dest, stride, bd);
  else
    vp9_highbd_idct4x4_1_add(input, dest, stride, bd);
}


void vp9_highbd_iwht4x4_add(const tran_low_t *input, uint8_t *dest, int stride,
                            int eob, int bd) {
  if (eob > 1)
    vp9_highbd_iwht4x4_16_add(input, dest, stride, bd);
  else
    vp9_highbd_iwht4x4_1_add(input, dest, stride, bd);
}

void vp9_highbd_idct8x8_add(const tran_low_t *input, uint8_t *dest, int stride,
                            int eob, int bd) {
  // If dc is 1, then input[0] is the reconstructed value, do not need
  // dequantization. Also, when dc is 1, dc is counted in eobs, namely eobs >=1.

  // The calculation can be simplified if there are not many non-zero dct
  // coefficients. Use eobs to decide what to do.
  // TODO(yunqingwang): "eobs = 1" case is also handled in vp9_short_idct8x8_c.
  // Combine that with code here.
  // DC only DCT coefficient
  if (eob == 1) {
    vp9_highbd_idct8x8_1_add(input, dest, stride, bd);
  } else if (eob <= 10) {
    vp9_highbd_idct8x8_10_add(input, dest, stride, bd);
  } else {
    vp9_highbd_idct8x8_64_add(input, dest, stride, bd);
  }
}

void vp9_highbd_idct16x16_add(const tran_low_t *input, uint8_t *dest,
                              int stride, int eob, int bd) {
  // The calculation can be simplified if there are not many non-zero dct
  // coefficients. Use eobs to separate different cases.
  // DC only DCT coefficient.
  if (eob == 1) {
    vp9_highbd_idct16x16_1_add(input, dest, stride, bd);
  } else if (eob <= 10) {
    vp9_highbd_idct16x16_10_add(input, dest, stride, bd);
  } else {
    vp9_highbd_idct16x16_256_add(input, dest, stride, bd);
  }
}

void vp9_highbd_idct32x32_add(const tran_low_t *input, uint8_t *dest,
                              int stride, int eob, int bd) {
  // Non-zero coeff only in upper-left 8x8
  if (eob == 1) {
    vp9_highbd_idct32x32_1_add(input, dest, stride, bd);
  } else if (eob <= 34) {
    vp9_highbd_idct32x32_34_add(input, dest, stride, bd);
  } else {
    vp9_highbd_idct32x32_1024_add(input, dest, stride, bd);
  }
}

// iht
void vp9_highbd_iht4x4_add(TX_TYPE tx_type, const tran_low_t *input,
                           uint8_t *dest, int stride, int eob, int bd) {
  if (tx_type == DCT_DCT) {
    vp9_highbd_idct4x4_add(input, dest, stride, eob, bd);
#if CONFIG_EXT_TX
  } else if (tx_type == FLIPADST_DCT) {
    flipud16(CONVERT_TO_SHORTPTR(dest), stride, 4);
    vp9_highbd_iht4x4_16_add(input, dest, stride, ADST_DCT, bd);
    flipud16(CONVERT_TO_SHORTPTR(dest), stride, 4);
  } else if (tx_type == DCT_FLIPADST) {
    fliplr16(CONVERT_TO_SHORTPTR(dest), stride, 4);
    vp9_highbd_iht4x4_16_add(input, dest, stride, DCT_ADST, bd);
    fliplr16(CONVERT_TO_SHORTPTR(dest), stride, 4);
  } else if (tx_type == FLIPADST_FLIPADST) {
    fliplrud16(CONVERT_TO_SHORTPTR(dest), stride, 4);
    vp9_highbd_iht4x4_16_add(input, dest, stride, ADST_ADST, bd);
    fliplrud16(CONVERT_TO_SHORTPTR(dest), stride, 4);
  } else if (tx_type == ADST_FLIPADST) {
    fliplr16(CONVERT_TO_SHORTPTR(dest), stride, 4);
    vp9_highbd_iht4x4_16_add(input, dest, stride, ADST_ADST, bd);
    fliplr16(CONVERT_TO_SHORTPTR(dest), stride, 4);
  } else if (tx_type == FLIPADST_ADST) {
    flipud16(CONVERT_TO_SHORTPTR(dest), stride, 4);
    vp9_highbd_iht4x4_16_add(input, dest, stride, ADST_ADST, bd);
    flipud16(CONVERT_TO_SHORTPTR(dest), stride, 4);
#endif  // CONFIG_EXT_TX
  } else {
    vp9_highbd_iht4x4_16_add(input, dest, stride, tx_type, bd);
  }
}

void vp9_highbd_iht8x8_add(TX_TYPE tx_type, const tran_low_t *input,
                           uint8_t *dest, int stride, int eob, int bd) {
  if (tx_type == DCT_DCT) {
    vp9_highbd_idct8x8_add(input, dest, stride, eob, bd);
#if CONFIG_EXT_TX
  } else if (tx_type == FLIPADST_DCT) {
    flipud16(CONVERT_TO_SHORTPTR(dest), stride, 8);
    vp9_highbd_iht8x8_64_add(input, dest, stride, ADST_DCT, bd);
    flipud16(CONVERT_TO_SHORTPTR(dest), stride, 8);
  } else if (tx_type == DCT_FLIPADST) {
    fliplr16(CONVERT_TO_SHORTPTR(dest), stride, 8);
    vp9_highbd_iht8x8_64_add(input, dest, stride, DCT_ADST, bd);
    fliplr16(CONVERT_TO_SHORTPTR(dest), stride, 8);
  } else if (tx_type == FLIPADST_FLIPADST) {
    fliplrud16(CONVERT_TO_SHORTPTR(dest), stride, 8);
    vp9_highbd_iht8x8_64_add(input, dest, stride, ADST_ADST, bd);
    fliplrud16(CONVERT_TO_SHORTPTR(dest), stride, 8);
  } else if (tx_type == ADST_FLIPADST) {
    fliplr16(CONVERT_TO_SHORTPTR(dest), stride, 8);
    vp9_highbd_iht8x8_64_add(input, dest, stride, ADST_ADST, bd);
    fliplr16(CONVERT_TO_SHORTPTR(dest), stride, 8);
  } else if (tx_type == FLIPADST_ADST) {
    flipud16(CONVERT_TO_SHORTPTR(dest), stride, 8);
    vp9_highbd_iht8x8_64_add(input, dest, stride, ADST_ADST, bd);
    flipud16(CONVERT_TO_SHORTPTR(dest), stride, 8);
#endif  // CONFIG_EXT_TX
  } else {
    vp9_highbd_iht8x8_64_add(input, dest, stride, tx_type, bd);
  }
}

void vp9_highbd_iht16x16_add(TX_TYPE tx_type, const tran_low_t *input,
                           uint8_t *dest, int stride, int eob, int bd) {
  if (tx_type == DCT_DCT) {
    vp9_highbd_idct16x16_add(input, dest, stride, eob, bd);
#if CONFIG_EXT_TX
  } else if (tx_type == FLIPADST_DCT) {
    flipud16(CONVERT_TO_SHORTPTR(dest), stride, 16);
    vp9_highbd_iht16x16_256_add(input, dest, stride, ADST_DCT, bd);
    flipud16(CONVERT_TO_SHORTPTR(dest), stride, 16);
  } else if (tx_type == DCT_FLIPADST) {
    fliplr16(CONVERT_TO_SHORTPTR(dest), stride, 16);
    vp9_highbd_iht16x16_256_add(input, dest, stride, DCT_ADST, bd);
    fliplr16(CONVERT_TO_SHORTPTR(dest), stride, 16);
  } else if (tx_type == FLIPADST_FLIPADST) {
    fliplrud16(CONVERT_TO_SHORTPTR(dest), stride, 16);
    vp9_highbd_iht16x16_256_add(input, dest, stride, ADST_ADST, bd);
    fliplrud16(CONVERT_TO_SHORTPTR(dest), stride, 16);
  } else if (tx_type == ADST_FLIPADST) {
    fliplr16(CONVERT_TO_SHORTPTR(dest), stride, 16);
    vp9_highbd_iht16x16_256_add(input, dest, stride, ADST_ADST, bd);
    fliplr16(CONVERT_TO_SHORTPTR(dest), stride, 16);
  } else if (tx_type == FLIPADST_ADST) {
    flipud16(CONVERT_TO_SHORTPTR(dest), stride, 16);
    vp9_highbd_iht16x16_256_add(input, dest, stride, ADST_ADST, bd);
    flipud16(CONVERT_TO_SHORTPTR(dest), stride, 16);
#endif  // CONFIG_EXT_TX
  } else {
    vp9_highbd_iht16x16_256_add(input, dest, stride, tx_type, bd);
  }
}

#if CONFIG_TX64X64
void vp9_highbd_idct64x64_4096_add_c(const tran_low_t *input, uint8_t *dest,
                                     int stride, int bd) {
  // vp9_clear_system_state();  // Make it simd safe : __asm emms;
  {
    double out[64 * 64], out2[64 * 64];
    int i, j;
    // First transform rows
    for (i = 0; i < 64; ++i) {
      double temp_in[64], temp_out[64];
      for (j = 0; j < 64; ++j)
        temp_in[j] = input[j + i * 64];
      butterfly_64_idct_1d(temp_in, temp_out, 1);
      for (j = 0; j < 64; ++j)
        out[j + i * 64] = temp_out[j];
    }
    // Then transform columns
    for (i = 0; i < 64; ++i) {
      double temp_in[64], temp_out[64];
      for (j = 0; j < 64; ++j)
        temp_in[j] = out[j * 64 + i];
      butterfly_64_idct_1d(temp_in, temp_out, 1);
      for (j = 0; j < 64; ++j)
        out2[j * 64 + i] = temp_out[j];
    }

    for (j = 0; j < 64; ++j) {
      for (i = 0; i < 64; ++i)
        dest[j * stride + i] = highbd_clip_pixel_add(
            dest[j * 64 + i], round(out2[j * 64 + i] / 128), bd);
      dest += stride;
    }
  }
  // vp9_clear_system_state();  // Make it simd safe : __asm emms;
}

void vp9_highbd_idct64x64_add(const tran_low_t *input, uint8_t *dest,
                              int stride, int eob, int bd) {
  (void) eob;
  vp9_highbd_idct64x64_4096_add_c(input, dest, stride, bd);
}
#endif  // CONFIG_TX64X64
#endif  // CONFIG_VP9_HIGHBITDEPTH
