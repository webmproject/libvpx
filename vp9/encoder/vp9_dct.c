/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */


#include <assert.h>
#include <math.h>
#include "./vpx_config.h"
#include "vp9/common/vp9_systemdependent.h"

#include "vp9/common/vp9_blockd.h"
#include "vp9/common/vp9_idct.h"

static void fdct4_1d(int16_t *input, int16_t *output) {
  int16_t step[4];
  int temp1, temp2;

  step[0] = input[0] + input[3];
  step[1] = input[1] + input[2];
  step[2] = input[1] - input[2];
  step[3] = input[0] - input[3];

  temp1 = (step[0] + step[1]) * cospi_16_64;
  temp2 = (step[0] - step[1]) * cospi_16_64;
  output[0] = dct_const_round_shift(temp1);
  output[2] = dct_const_round_shift(temp2);
  temp1 = step[2] * cospi_24_64 + step[3] * cospi_8_64;
  temp2 = -step[2] * cospi_8_64 + step[3] * cospi_24_64;
  output[1] = dct_const_round_shift(temp1);
  output[3] = dct_const_round_shift(temp2);
}

void vp9_short_fdct4x4_c(short *input, short *output, int pitch) {
  int16_t out[4 * 4];
  int16_t *outptr = &out[0];
  const int short_pitch = pitch >> 1;
  int i, j;
  int16_t temp_in[4], temp_out[4];
  // First transform cols
  for (i = 0; i < 4; ++i) {
    for (j = 0; j < 4; ++j)
      temp_in[j] = input[j * short_pitch + i] << 4;
    if (i == 0 && temp_in[0])
      temp_in[0] += 1;
    fdct4_1d(temp_in, temp_out);
    for (j = 0; j < 4; ++j)
      outptr[j * 4 + i] = temp_out[j];
  }
  for (i = 0; i < 4; ++i) {
    for (j = 0; j < 4; ++j)
      temp_in[j] = out[j + i * 4];
    fdct4_1d(temp_in, temp_out);
    for (j = 0; j < 4; ++j)
        output[j + i * 4] = (temp_out[j] + 1) >> 2;
  }
}

static void fadst4_1d(int16_t *input, int16_t *output) {
  int x0, x1, x2, x3;
  int s0, s1, s2, s3, s4, s5, s6, s7;

  x0 = input[0];
  x1 = input[1];
  x2 = input[2];
  x3 = input[3];

  if (!(x0 | x1 | x2 | x3)) {
    output[0] = output[1] = output[2] = output[3] = 0;
    return;
  }

  s0 = sinpi_1_9 * x0;
  s1 = sinpi_4_9 * x0;
  s2 = sinpi_2_9 * x1;
  s3 = sinpi_1_9 * x1;
  s4 = sinpi_3_9 * x2;
  s5 = sinpi_4_9 * x3;
  s6 = sinpi_2_9 * x3;
  s7 = x0 + x1 - x3;

  x0 = s0 + s2 + s5;
  x1 = sinpi_3_9 * s7;
  x2 = s1 - s3 + s6;
  x3 = s4;

  s0 = x0 + x3;
  s1 = x1;
  s2 = x2 - x3;
  s3 = x2 - x0 + x3;

  // 1-D transform scaling factor is sqrt(2).
  output[0] = dct_const_round_shift(s0);
  output[1] = dct_const_round_shift(s1);
  output[2] = dct_const_round_shift(s2);
  output[3] = dct_const_round_shift(s3);
}

void vp9_short_fht4x4_c(int16_t *input, int16_t *output,
                        int pitch, TX_TYPE tx_type) {
  int16_t out[4 * 4];
  int16_t *outptr = &out[0];
  const int short_pitch = pitch >> 1;
  int i, j;
  int16_t temp_in[4], temp_out[4];

  void (*fwdr)(int16_t*, int16_t*);
  void (*fwdc)(int16_t*, int16_t*);

  switch (tx_type) {
    case ADST_ADST:
      fwdc = &fadst4_1d;
      fwdr = &fadst4_1d;
      break;
    case ADST_DCT:
      fwdc = &fadst4_1d;
      fwdr = &fdct4_1d;
      break;
    case DCT_ADST:
      fwdc = &fdct4_1d;
      fwdr = &fadst4_1d;
      break;
    case DCT_DCT:
      fwdc = &fdct4_1d;
      fwdr = &fdct4_1d;
      break;
    default:
      assert(0);
  }


  // column transform
  for (i = 0; i < 4; ++i) {
    for (j = 0; j < 4; ++j)
      temp_in[j] = input[j * short_pitch + i] << 4;
    if (i == 0 && temp_in[0])
      temp_in[0] += 1;
    fwdc(temp_in, temp_out);
    for (j = 0; j < 4; ++j)
      outptr[j * 4 + i] = temp_out[j];
  }

  // row transform
  for (i = 0; i < 4; ++i) {
    for (j = 0; j < 4; ++j)
      temp_in[j] = out[j + i * 4];
    fwdr(temp_in, temp_out);
    for (j = 0; j < 4; ++j)
      output[j + i * 4] = (temp_out[j] + 1) >> 2;
  }
}

void vp9_short_fdct8x4_c(short *input, short *output, int pitch)
{
    vp9_short_fdct4x4_c(input,   output,    pitch);
    vp9_short_fdct4x4_c(input + 4, output + 16, pitch);
}

static void fdct8_1d(int16_t *input, int16_t *output) {
  int16_t step[8];
  int temp1, temp2;

  // stage 1
  step[0] = input[0] + input[7];
  step[1] = input[1] + input[6];
  step[2] = input[2] + input[5];
  step[3] = input[3] + input[4];
  step[4] = input[3] - input[4];
  step[5] = input[2] - input[5];
  step[6] = input[1] - input[6];
  step[7] = input[0] - input[7];

  fdct4_1d(step, step);

  // Stage 2
  output[4] = step[4];
  temp1 = (-step[5] + step[6]) * cospi_16_64;
  temp2 = (step[6] + step[5]) * cospi_16_64;
  output[5] = dct_const_round_shift(temp1);
  output[6] = dct_const_round_shift(temp2);
  output[7] = step[7];

  // Stage 3
  step[4] = output[4] + output[5];
  step[5] = -output[5] + output[4];
  step[6] = -output[6] + output[7];
  step[7] = output[7] + output[6];

  // Stage 4
  output[0] = step[0];
  output[4] = step[2];
  output[2] = step[1];
  output[6] = step[3];

  temp1 = step[4] * cospi_28_64 + step[7] * cospi_4_64;
  temp2 = step[5] * cospi_12_64 + step[6] * cospi_20_64;
  output[1] = dct_const_round_shift(temp1);
  output[5] = dct_const_round_shift(temp2);
  temp1 = step[6] * cospi_12_64 + step[5] * -cospi_20_64;
  temp2 = step[7] * cospi_28_64 + step[4] * -cospi_4_64;
  output[3] = dct_const_round_shift(temp1);
  output[7] = dct_const_round_shift(temp2);
}

void vp9_short_fdct8x8_c(int16_t *input, int16_t *output, int pitch) {
  int shortpitch = pitch >> 1;
  int i, j;
  int16_t out[64];
  int16_t temp_in[8], temp_out[8];

  // First transform columns
  for (i = 0; i < 8; i++) {
    for (j = 0; j < 8; j++)
      temp_in[j] = input[j * shortpitch + i] << 2;
    fdct8_1d(temp_in, temp_out);
    for (j = 0; j < 8; j++)
      out[j * 8 + i] = temp_out[j];
  }

  // Then transform rows
  for (i = 0; i < 8; ++i) {
    for (j = 0; j < 8; ++j)
      temp_in[j] = out[j + i * 8];
    fdct8_1d(temp_in, temp_out);
    for (j = 0; j < 8; ++j)
      output[j + i * 8] = temp_out[j] / 2;
  }
}

static void fadst8_1d(int16_t *input, int16_t *output) {
  int x0, x1, x2, x3, x4, x5, x6, x7;
  int s0, s1, s2, s3, s4, s5, s6, s7;

  x0 = input[7];
  x1 = input[0];
  x2 = input[5];
  x3 = input[2];
  x4 = input[3];
  x5 = input[4];
  x6 = input[1];
  x7 = input[6];

  // stage 1
  s0 = cospi_2_64  * x0 + cospi_30_64 * x1;
  s1 = cospi_30_64 * x0 - cospi_2_64  * x1;
  s2 = cospi_10_64 * x2 + cospi_22_64 * x3;
  s3 = cospi_22_64 * x2 - cospi_10_64 * x3;
  s4 = cospi_18_64 * x4 + cospi_14_64 * x5;
  s5 = cospi_14_64 * x4 - cospi_18_64 * x5;
  s6 = cospi_26_64 * x6 + cospi_6_64  * x7;
  s7 = cospi_6_64  * x6 - cospi_26_64 * x7;

  x0 = dct_const_round_shift(s0 + s4);
  x1 = dct_const_round_shift(s1 + s5);
  x2 = dct_const_round_shift(s2 + s6);
  x3 = dct_const_round_shift(s3 + s7);
  x4 = dct_const_round_shift(s0 - s4);
  x5 = dct_const_round_shift(s1 - s5);
  x6 = dct_const_round_shift(s2 - s6);
  x7 = dct_const_round_shift(s3 - s7);

  // stage 2
  s0 = x0;
  s1 = x1;
  s2 = x2;
  s3 = x3;
  s4 = cospi_8_64  * x4 + cospi_24_64 * x5;
  s5 = cospi_24_64 * x4 - cospi_8_64  * x5;
  s6 = - cospi_24_64 * x6 + cospi_8_64  * x7;
  s7 =   cospi_8_64  * x6 + cospi_24_64 * x7;

  x0 = s0 + s2;
  x1 = s1 + s3;
  x2 = s0 - s2;
  x3 = s1 - s3;
  x4 = dct_const_round_shift(s4 + s6);
  x5 = dct_const_round_shift(s5 + s7);
  x6 = dct_const_round_shift(s4 - s6);
  x7 = dct_const_round_shift(s5 - s7);

  // stage 3
  s2 = cospi_16_64 * (x2 + x3);
  s3 = cospi_16_64 * (x2 - x3);
  s6 = cospi_16_64 * (x6 + x7);
  s7 = cospi_16_64 * (x6 - x7);

  x2 = dct_const_round_shift(s2);
  x3 = dct_const_round_shift(s3);
  x6 = dct_const_round_shift(s6);
  x7 = dct_const_round_shift(s7);

  output[0] =   x0;
  output[1] = - x4;
  output[2] =   x6;
  output[3] = - x2;
  output[4] =   x3;
  output[5] = - x7;
  output[6] =   x5;
  output[7] = - x1;
}

void vp9_short_fht8x8_c(int16_t *input, int16_t *output,
                        int pitch, TX_TYPE tx_type) {
  int16_t out[64];
  int16_t *outptr = &out[0];
  const int short_pitch = pitch >> 1;
  int i, j;
  int16_t temp_in[8], temp_out[8];

  void (*fwdr)(int16_t*, int16_t*);
  void (*fwdc)(int16_t*, int16_t*);

  switch (tx_type) {
    case ADST_ADST:
      fwdc = &fadst8_1d;
      fwdr = &fadst8_1d;
      break;
    case ADST_DCT:
      fwdc = &fadst8_1d;
      fwdr = &fdct8_1d;
      break;
    case DCT_ADST:
      fwdc = &fdct8_1d;
      fwdr = &fadst8_1d;
      break;
    case DCT_DCT:
      fwdc = &fdct8_1d;
      fwdr = &fdct8_1d;
      break;
    default:
      assert(0);
  }

  // column transform
  for (i = 0; i < 8; ++i) {
    for (j = 0; j < 8; ++j)
      temp_in[j] = input[j * short_pitch + i] << 2;
    fwdc(temp_in, temp_out);
    for (j = 0; j < 8; ++j)
      outptr[j * 8 + i] = temp_out[j];
  }

  // row transform
  for (i = 0; i < 8; ++i) {
    for (j = 0; j < 8; ++j)
      temp_in[j] = out[j + i * 8];
    fwdr(temp_in, temp_out);
    for (j = 0; j < 8; ++j)
      output[j + i * 8] = temp_out[j] >> 1;
  }
}

void vp9_short_walsh4x4_x8_c(short *input, short *output, int pitch) {
  int i;
  int a1, b1, c1, d1;
  short *ip = input;
  short *op = output;
  int pitch_short = pitch >> 1;

  for (i = 0; i < 4; i++) {
    a1 = ip[0 * pitch_short] + ip[3 * pitch_short];
    b1 = ip[1 * pitch_short] + ip[2 * pitch_short];
    c1 = ip[1 * pitch_short] - ip[2 * pitch_short];
    d1 = ip[0 * pitch_short] - ip[3 * pitch_short];

    op[0] = (a1 + b1 + 1) >> 1;
    op[4] = (c1 + d1) >> 1;
    op[8] = (a1 - b1) >> 1;
    op[12] = (d1 - c1) >> 1;

    ip++;
    op++;
  }
  ip = output;
  op = output;

  for (i = 0; i < 4; i++) {
    a1 = ip[0] + ip[3];
    b1 = ip[1] + ip[2];
    c1 = ip[1] - ip[2];
    d1 = ip[0] - ip[3];

    op[0] = ((a1 + b1 + 1) >> 1) << WHT_UPSCALE_FACTOR;
    op[1] = ((c1 + d1) >> 1) << WHT_UPSCALE_FACTOR;
    op[2] = ((a1 - b1) >> 1) << WHT_UPSCALE_FACTOR;
    op[3] = ((d1 - c1) >> 1) << WHT_UPSCALE_FACTOR;

    ip += 4;
    op += 4;
  }
}

void vp9_short_walsh8x4_x8_c(short *input, short *output, int pitch) {
  vp9_short_walsh4x4_x8_c(input,   output,    pitch);
  vp9_short_walsh4x4_x8_c(input + 4, output + 16, pitch);
}


// Rewrote to use same algorithm as others.
static void fdct16_1d(int16_t input[16], int16_t output[16]) {
  int16_t step[16];
  int temp1, temp2;

  // step 1
  step[ 0] = input[0] + input[15];
  step[ 1] = input[1] + input[14];
  step[ 2] = input[2] + input[13];
  step[ 3] = input[3] + input[12];
  step[ 4] = input[4] + input[11];
  step[ 5] = input[5] + input[10];
  step[ 6] = input[6] + input[ 9];
  step[ 7] = input[7] + input[ 8];
  step[ 8] = input[7] - input[ 8];
  step[ 9] = input[6] - input[ 9];
  step[10] = input[5] - input[10];
  step[11] = input[4] - input[11];
  step[12] = input[3] - input[12];
  step[13] = input[2] - input[13];
  step[14] = input[1] - input[14];
  step[15] = input[0] - input[15];

  fdct8_1d(step, step);

  // step 2
  output[8] = step[8];
  output[9] = step[9];
  temp1 = (-step[10] + step[13]) * cospi_16_64;
  temp2 = (-step[11] + step[12]) * cospi_16_64;
  output[10] = dct_const_round_shift(temp1);
  output[11] = dct_const_round_shift(temp2);
  temp1 = (step[11] + step[12]) * cospi_16_64;
  temp2 = (step[10] + step[13]) * cospi_16_64;
  output[12] = dct_const_round_shift(temp1);
  output[13] = dct_const_round_shift(temp2);
  output[14] = step[14];
  output[15] = step[15];

  // step 3
  step[ 8] = output[8] + output[11];
  step[ 9] = output[9] + output[10];
  step[ 10] = output[9] - output[10];
  step[ 11] = output[8] - output[11];
  step[ 12] = -output[12] + output[15];
  step[ 13] = -output[13] + output[14];
  step[ 14] = output[13] + output[14];
  step[ 15] = output[12] + output[15];

  // step 4
  output[8] = step[8];
  temp1 = -step[9] * cospi_8_64 + step[14] * cospi_24_64;
  temp2 = -step[10] * cospi_24_64 - step[13] * cospi_8_64;
  output[9] = dct_const_round_shift(temp1);
  output[10] = dct_const_round_shift(temp2);
  output[11] = step[11];
  output[12] = step[12];
  temp1 = -step[10] * cospi_8_64 + step[13] * cospi_24_64;
  temp2 = step[9] * cospi_24_64 + step[14] * cospi_8_64;
  output[13] = dct_const_round_shift(temp1);
  output[14] = dct_const_round_shift(temp2);
  output[15] = step[15];

  // step 5
  step[8] = output[8] + output[9];
  step[9] = output[8] - output[9];
  step[10] = -output[10] + output[11];
  step[11] = output[10] + output[11];
  step[12] = output[12] + output[13];
  step[13] = output[12] - output[13];
  step[14] = -output[14] + output[15];
  step[15] = output[14] + output[15];

  // step 6
  output[0] = step[0];
  output[8] = step[4];
  output[4] = step[2];
  output[12] = step[6];
  output[2] = step[1];
  output[10] = step[5];
  output[6] = step[3];
  output[14] = step[7];

  temp1 = step[8] * cospi_30_64 + step[15] * cospi_2_64;
  temp2 = step[9] * cospi_14_64 + step[14] * cospi_18_64;
  output[1] = dct_const_round_shift(temp1);
  output[9] = dct_const_round_shift(temp2);

  temp1 = step[10] * cospi_22_64 + step[13] * cospi_10_64;
  temp2 = step[11] * cospi_6_64 + step[12] * cospi_26_64;
  output[5] = dct_const_round_shift(temp1);
  output[13] = dct_const_round_shift(temp2);

  temp1 = -step[11] * cospi_26_64 + step[12] * cospi_6_64;
  temp2 = -step[10] * cospi_10_64 + step[13] * cospi_22_64;
  output[3] = dct_const_round_shift(temp1);
  output[11] = dct_const_round_shift(temp2);

  temp1 = -step[9] * cospi_18_64 + step[14] * cospi_14_64;
  temp2 = -step[8] * cospi_2_64 + step[15] * cospi_30_64;
  output[7] = dct_const_round_shift(temp1);
  output[15] = dct_const_round_shift(temp2);
}

void vp9_short_fdct16x16_c(int16_t *input, int16_t *out, int pitch) {
  int shortpitch = pitch >> 1;
  int i, j;
  int16_t output[256];
  int16_t temp_in[16], temp_out[16];

  // First transform columns
  for (i = 0; i < 16; i++) {
    for (j = 0; j < 16; j++)
      temp_in[j] = input[j * shortpitch + i] << 2;
    fdct16_1d(temp_in, temp_out);
    for (j = 0; j < 16; j++)
      output[j * 16 + i] = (temp_out[j] + 1) >> 2;
  }

  // Then transform rows
  for (i = 0; i < 16; ++i) {
    for (j = 0; j < 16; ++j)
      temp_in[j] = output[j + i * 16];
    fdct16_1d(temp_in, temp_out);
    for (j = 0; j < 16; ++j)
      out[j + i * 16] = temp_out[j];
  }
}

void fadst16_1d(int16_t *input, int16_t *output) {
  int x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14, x15;
  int s0, s1, s2, s3, s4, s5, s6, s7, s8, s9, s10, s11, s12, s13, s14, s15;

  x0 = input[15];
  x1 = input[0];
  x2 = input[13];
  x3 = input[2];
  x4 = input[11];
  x5 = input[4];
  x6 = input[9];
  x7 = input[6];
  x8 = input[7];
  x9 = input[8];
  x10 = input[5];
  x11 = input[10];
  x12 = input[3];
  x13 = input[12];
  x14 = input[1];
  x15 = input[14];

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

  x0 = dct_const_round_shift(s0 + s8);
  x1 = dct_const_round_shift(s1 + s9);
  x2 = dct_const_round_shift(s2 + s10);
  x3 = dct_const_round_shift(s3 + s11);
  x4 = dct_const_round_shift(s4 + s12);
  x5 = dct_const_round_shift(s5 + s13);
  x6 = dct_const_round_shift(s6 + s14);
  x7 = dct_const_round_shift(s7 + s15);
  x8  = dct_const_round_shift(s0 - s8);
  x9  = dct_const_round_shift(s1 - s9);
  x10 = dct_const_round_shift(s2 - s10);
  x11 = dct_const_round_shift(s3 - s11);
  x12 = dct_const_round_shift(s4 - s12);
  x13 = dct_const_round_shift(s5 - s13);
  x14 = dct_const_round_shift(s6 - s14);
  x15 = dct_const_round_shift(s7 - s15);

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

  x0 = s0 + s4;
  x1 = s1 + s5;
  x2 = s2 + s6;
  x3 = s3 + s7;
  x4 = s0 - s4;
  x5 = s1 - s5;
  x6 = s2 - s6;
  x7 = s3 - s7;
  x8 = dct_const_round_shift(s8 + s12);
  x9 = dct_const_round_shift(s9 + s13);
  x10 = dct_const_round_shift(s10 + s14);
  x11 = dct_const_round_shift(s11 + s15);
  x12 = dct_const_round_shift(s8 - s12);
  x13 = dct_const_round_shift(s9 - s13);
  x14 = dct_const_round_shift(s10 - s14);
  x15 = dct_const_round_shift(s11 - s15);

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

  x0 = s0 + s2;
  x1 = s1 + s3;
  x2 = s0 - s2;
  x3 = s1 - s3;
  x4 = dct_const_round_shift(s4 + s6);
  x5 = dct_const_round_shift(s5 + s7);
  x6 = dct_const_round_shift(s4 - s6);
  x7 = dct_const_round_shift(s5 - s7);
  x8 = s8 + s10;
  x9 = s9 + s11;
  x10 = s8 - s10;
  x11 = s9 - s11;
  x12 = dct_const_round_shift(s12 + s14);
  x13 = dct_const_round_shift(s13 + s15);
  x14 = dct_const_round_shift(s12 - s14);
  x15 = dct_const_round_shift(s13 - s15);

  // stage 4
  s2 = (- cospi_16_64) * (x2 + x3);
  s3 = cospi_16_64 * (x2 - x3);
  s6 = cospi_16_64 * (x6 + x7);
  s7 = cospi_16_64 * (- x6 + x7);
  s10 = cospi_16_64 * (x10 + x11);
  s11 = cospi_16_64 * (- x10 + x11);
  s14 = (- cospi_16_64) * (x14 + x15);
  s15 = cospi_16_64 * (x14 - x15);

  x2 = dct_const_round_shift(s2);
  x3 = dct_const_round_shift(s3);
  x6 = dct_const_round_shift(s6);
  x7 = dct_const_round_shift(s7);
  x10 = dct_const_round_shift(s10);
  x11 = dct_const_round_shift(s11);
  x14 = dct_const_round_shift(s14);
  x15 = dct_const_round_shift(s15);

  output[0] = x0;
  output[1] = - x8;
  output[2] = x12;
  output[3] = - x4;
  output[4] = x6;
  output[5] = x14;
  output[6] = x10;
  output[7] = x2;
  output[8] = x3;
  output[9] =  x11;
  output[10] = x15;
  output[11] = x7;
  output[12] = x5;
  output[13] = - x13;
  output[14] = x9;
  output[15] = - x1;
}

void vp9_short_fht16x16_c(int16_t *input, int16_t *output,
                          int pitch, TX_TYPE tx_type) {
  int16_t out[256];
  int16_t *outptr = &out[0];
  const int short_pitch = pitch >> 1;
  int i, j;
  int16_t temp_in[16], temp_out[16];

  void (*fwdr)(int16_t*, int16_t*);
  void (*fwdc)(int16_t*, int16_t*);

  switch (tx_type) {
    case ADST_ADST:
      fwdc = &fadst16_1d;
      fwdr = &fadst16_1d;
      break;
    case ADST_DCT:
      fwdc = &fadst16_1d;
      fwdr = &fdct16_1d;
      break;
    case DCT_ADST:
      fwdc = &fdct16_1d;
      fwdr = &fadst16_1d;
      break;
    case DCT_DCT:
      fwdc = &fdct16_1d;
      fwdr = &fdct16_1d;
      break;
    default:
      assert(0);
  }

  // column transform
  for (i = 0; i < 16; ++i) {
    for (j = 0; j < 16; ++j)
      temp_in[j] = input[j * short_pitch + i] << 2;
    fwdc(temp_in, temp_out);
    for (j = 0; j < 16; ++j)
      outptr[j * 16 + i] = (temp_out[j] + 1 + (temp_out[j] > 0)) >> 2;
  }

  // row transform
  for (i = 0; i < 16; ++i) {
    for (j = 0; j < 16; ++j)
      temp_in[j] = out[j + i * 16];
    fwdr(temp_in, temp_out);
    for (j = 0; j < 16; ++j)
      output[j + i * 16] = temp_out[j];
  }
}


static void dct32_1d(int *input, int *output) {
  int step[32];
  // Stage 1
  step[0] = input[0] + input[(32 - 1)];
  step[1] = input[1] + input[(32 - 2)];
  step[2] = input[2] + input[(32 - 3)];
  step[3] = input[3] + input[(32 - 4)];
  step[4] = input[4] + input[(32 - 5)];
  step[5] = input[5] + input[(32 - 6)];
  step[6] = input[6] + input[(32 - 7)];
  step[7] = input[7] + input[(32 - 8)];
  step[8] = input[8] + input[(32 - 9)];
  step[9] = input[9] + input[(32 - 10)];
  step[10] = input[10] + input[(32 - 11)];
  step[11] = input[11] + input[(32 - 12)];
  step[12] = input[12] + input[(32 - 13)];
  step[13] = input[13] + input[(32 - 14)];
  step[14] = input[14] + input[(32 - 15)];
  step[15] = input[15] + input[(32 - 16)];
  step[16] = -input[16] + input[(32 - 17)];
  step[17] = -input[17] + input[(32 - 18)];
  step[18] = -input[18] + input[(32 - 19)];
  step[19] = -input[19] + input[(32 - 20)];
  step[20] = -input[20] + input[(32 - 21)];
  step[21] = -input[21] + input[(32 - 22)];
  step[22] = -input[22] + input[(32 - 23)];
  step[23] = -input[23] + input[(32 - 24)];
  step[24] = -input[24] + input[(32 - 25)];
  step[25] = -input[25] + input[(32 - 26)];
  step[26] = -input[26] + input[(32 - 27)];
  step[27] = -input[27] + input[(32 - 28)];
  step[28] = -input[28] + input[(32 - 29)];
  step[29] = -input[29] + input[(32 - 30)];
  step[30] = -input[30] + input[(32 - 31)];
  step[31] = -input[31] + input[(32 - 32)];

  // Stage 2
  output[0] = step[0] + step[16 - 1];
  output[1] = step[1] + step[16 - 2];
  output[2] = step[2] + step[16 - 3];
  output[3] = step[3] + step[16 - 4];
  output[4] = step[4] + step[16 - 5];
  output[5] = step[5] + step[16 - 6];
  output[6] = step[6] + step[16 - 7];
  output[7] = step[7] + step[16 - 8];
  output[8] = -step[8] + step[16 - 9];
  output[9] = -step[9] + step[16 - 10];
  output[10] = -step[10] + step[16 - 11];
  output[11] = -step[11] + step[16 - 12];
  output[12] = -step[12] + step[16 - 13];
  output[13] = -step[13] + step[16 - 14];
  output[14] = -step[14] + step[16 - 15];
  output[15] = -step[15] + step[16 - 16];

  output[16] = step[16];
  output[17] = step[17];
  output[18] = step[18];
  output[19] = step[19];

  output[20] = dct_32_round((-step[20] + step[27]) * cospi_16_64);
  output[21] = dct_32_round((-step[21] + step[26]) * cospi_16_64);
  output[22] = dct_32_round((-step[22] + step[25]) * cospi_16_64);
  output[23] = dct_32_round((-step[23] + step[24]) * cospi_16_64);

  output[24] = dct_32_round((step[24] + step[23]) * cospi_16_64);
  output[25] = dct_32_round((step[25] + step[22]) * cospi_16_64);
  output[26] = dct_32_round((step[26] + step[21]) * cospi_16_64);
  output[27] = dct_32_round((step[27] + step[20]) * cospi_16_64);

  output[28] = step[28];
  output[29] = step[29];
  output[30] = step[30];
  output[31] = step[31];

  // Stage 3
  step[0] = output[0] + output[(8 - 1)];
  step[1] = output[1] + output[(8 - 2)];
  step[2] = output[2] + output[(8 - 3)];
  step[3] = output[3] + output[(8 - 4)];
  step[4] = -output[4] + output[(8 - 5)];
  step[5] = -output[5] + output[(8 - 6)];
  step[6] = -output[6] + output[(8 - 7)];
  step[7] = -output[7] + output[(8 - 8)];
  step[8] = output[8];
  step[9] = output[9];
  step[10] = dct_32_round((-output[10] + output[13]) * cospi_16_64);
  step[11] = dct_32_round((-output[11] + output[12]) * cospi_16_64);
  step[12] = dct_32_round((output[12] + output[11]) * cospi_16_64);
  step[13] = dct_32_round((output[13] + output[10]) * cospi_16_64);
  step[14] = output[14];
  step[15] = output[15];

  step[16] = output[16] + output[23];
  step[17] = output[17] + output[22];
  step[18] = output[18] + output[21];
  step[19] = output[19] + output[20];
  step[20] = -output[20] + output[19];
  step[21] = -output[21] + output[18];
  step[22] = -output[22] + output[17];
  step[23] = -output[23] + output[16];
  step[24] = -output[24] + output[31];
  step[25] = -output[25] + output[30];
  step[26] = -output[26] + output[29];
  step[27] = -output[27] + output[28];
  step[28] = output[28] + output[27];
  step[29] = output[29] + output[26];
  step[30] = output[30] + output[25];
  step[31] = output[31] + output[24];

  // Stage 4
  output[0] = step[0] + step[3];
  output[1] = step[1] + step[2];
  output[2] = -step[2] + step[1];
  output[3] = -step[3] + step[0];
  output[4] = step[4];
  output[5] = dct_32_round((-step[5] + step[6]) * cospi_16_64);
  output[6] = dct_32_round((step[6] + step[5]) * cospi_16_64);
  output[7] = step[7];
  output[8] = step[8] + step[11];
  output[9] = step[9] + step[10];
  output[10] = -step[10] + step[9];
  output[11] = -step[11] + step[8];
  output[12] = -step[12] + step[15];
  output[13] = -step[13] + step[14];
  output[14] = step[14] + step[13];
  output[15] = step[15] + step[12];

  output[16] = step[16];
  output[17] = step[17];
  output[18] = dct_32_round(step[18] * -cospi_8_64 + step[29] * cospi_24_64);
  output[19] = dct_32_round(step[19] * -cospi_8_64 + step[28] * cospi_24_64);
  output[20] = dct_32_round(step[20] * -cospi_24_64 + step[27] * -cospi_8_64);
  output[21] = dct_32_round(step[21] * -cospi_24_64 + step[26] * -cospi_8_64);
  output[22] = step[22];
  output[23] = step[23];
  output[24] = step[24];
  output[25] = step[25];
  output[26] = dct_32_round(step[26] * cospi_24_64 + step[21] * -cospi_8_64);
  output[27] = dct_32_round(step[27] * cospi_24_64 + step[20] * -cospi_8_64);
  output[28] = dct_32_round(step[28] * cospi_8_64 + step[19] * cospi_24_64);
  output[29] = dct_32_round(step[29] * cospi_8_64 + step[18] * cospi_24_64);
  output[30] = step[30];
  output[31] = step[31];

  // Stage 5
  step[0] = dct_32_round((output[0] + output[1]) * cospi_16_64);
  step[1] = dct_32_round((-output[1] + output[0]) * cospi_16_64);
  step[2] = dct_32_round(output[2] * cospi_24_64 + output[3] * cospi_8_64);
  step[3] = dct_32_round(output[3] * cospi_24_64 - output[2] * cospi_8_64);
  step[4] = output[4] + output[5];
  step[5] = -output[5] + output[4];
  step[6] = -output[6] + output[7];
  step[7] = output[7] + output[6];
  step[8] = output[8];
  step[9] = dct_32_round(output[9] * -cospi_8_64 + output[14] * cospi_24_64);
  step[10] = dct_32_round(output[10] * -cospi_24_64 + output[13] * -cospi_8_64);
  step[11] = output[11];
  step[12] = output[12];
  step[13] = dct_32_round(output[13] * cospi_24_64 + output[10] * -cospi_8_64);
  step[14] = dct_32_round(output[14] * cospi_8_64 + output[9] * cospi_24_64);
  step[15] = output[15];

  step[16] = output[16] + output[19];
  step[17] = output[17] + output[18];
  step[18] = -output[18] + output[17];
  step[19] = -output[19] + output[16];
  step[20] = -output[20] + output[23];
  step[21] = -output[21] + output[22];
  step[22] = output[22] + output[21];
  step[23] = output[23] + output[20];
  step[24] = output[24] + output[27];
  step[25] = output[25] + output[26];
  step[26] = -output[26] + output[25];
  step[27] = -output[27] + output[24];
  step[28] = -output[28] + output[31];
  step[29] = -output[29] + output[30];
  step[30] = output[30] + output[29];
  step[31] = output[31] + output[28];

  // Stage 6
  output[0] = step[0];
  output[1] = step[1];
  output[2] = step[2];
  output[3] = step[3];
  output[4] = dct_32_round(step[4] * cospi_28_64 + step[7] * cospi_4_64);
  output[5] = dct_32_round(step[5] * cospi_12_64 + step[6] * cospi_20_64);
  output[6] = dct_32_round(step[6] * cospi_12_64 + step[5] * -cospi_20_64);
  output[7] = dct_32_round(step[7] * cospi_28_64 + step[4] * -cospi_4_64);
  output[8] = step[8] + step[9];
  output[9] = -step[9] + step[8];
  output[10] = -step[10] + step[11];
  output[11] = step[11] + step[10];
  output[12] = step[12] + step[13];
  output[13] = -step[13] + step[12];
  output[14] = -step[14] + step[15];
  output[15] = step[15] + step[14];

  output[16] = step[16];
  output[17] = dct_32_round(step[17] * -cospi_4_64 + step[30] * cospi_28_64);
  output[18] = dct_32_round(step[18] * -cospi_28_64 + step[29] * -cospi_4_64);
  output[19] = step[19];
  output[20] = step[20];
  output[21] = dct_32_round(step[21] * -cospi_20_64 + step[26] * cospi_12_64);
  output[22] = dct_32_round(step[22] * -cospi_12_64 + step[25] * -cospi_20_64);
  output[23] = step[23];
  output[24] = step[24];
  output[25] = dct_32_round(step[25] * cospi_12_64 + step[22] * -cospi_20_64);
  output[26] = dct_32_round(step[26] * cospi_20_64 + step[21] * cospi_12_64);
  output[27] = step[27];
  output[28] = step[28];
  output[29] = dct_32_round(step[29] * cospi_28_64 + step[18] * -cospi_4_64);
  output[30] = dct_32_round(step[30] * cospi_4_64 + step[17] * cospi_28_64);
  output[31] = step[31];

  // Stage 7
  step[0] = output[0];
  step[1] = output[1];
  step[2] = output[2];
  step[3] = output[3];
  step[4] = output[4];
  step[5] = output[5];
  step[6] = output[6];
  step[7] = output[7];
  step[8] = dct_32_round(output[8] * cospi_30_64 + output[15] * cospi_2_64);
  step[9] = dct_32_round(output[9] * cospi_14_64 + output[14] * cospi_18_64);
  step[10] = dct_32_round(output[10] * cospi_22_64 + output[13] * cospi_10_64);
  step[11] = dct_32_round(output[11] * cospi_6_64 + output[12] * cospi_26_64);
  step[12] = dct_32_round(output[12] * cospi_6_64 + output[11] * -cospi_26_64);
  step[13] = dct_32_round(output[13] * cospi_22_64 + output[10] * -cospi_10_64);
  step[14] = dct_32_round(output[14] * cospi_14_64 + output[9] * -cospi_18_64);
  step[15] = dct_32_round(output[15] * cospi_30_64 + output[8] * -cospi_2_64);

  step[16] = output[16] + output[17];
  step[17] = -output[17] + output[16];
  step[18] = -output[18] + output[19];
  step[19] = output[19] + output[18];
  step[20] = output[20] + output[21];
  step[21] = -output[21] + output[20];
  step[22] = -output[22] + output[23];
  step[23] = output[23] + output[22];
  step[24] = output[24] + output[25];
  step[25] = -output[25] + output[24];
  step[26] = -output[26] + output[27];
  step[27] = output[27] + output[26];
  step[28] = output[28] + output[29];
  step[29] = -output[29] + output[28];
  step[30] = -output[30] + output[31];
  step[31] = output[31] + output[30];

  // Final stage --- outputs indices are bit-reversed.
  output[0]  = step[0];
  output[16] = step[1];
  output[8]  = step[2];
  output[24] = step[3];
  output[4]  = step[4];
  output[20] = step[5];
  output[12] = step[6];
  output[28] = step[7];
  output[2]  = step[8];
  output[18] = step[9];
  output[10] = step[10];
  output[26] = step[11];
  output[6]  = step[12];
  output[22] = step[13];
  output[14] = step[14];
  output[30] = step[15];

  output[1]  = dct_32_round(step[16] * cospi_31_64 + step[31] * cospi_1_64);
  output[17] = dct_32_round(step[17] * cospi_15_64 + step[30] * cospi_17_64);
  output[9]  = dct_32_round(step[18] * cospi_23_64 + step[29] * cospi_9_64);
  output[25] = dct_32_round(step[19] * cospi_7_64 + step[28] * cospi_25_64);
  output[5]  = dct_32_round(step[20] * cospi_27_64 + step[27] * cospi_5_64);
  output[21] = dct_32_round(step[21] * cospi_11_64 + step[26] * cospi_21_64);
  output[13] = dct_32_round(step[22] * cospi_19_64 + step[25] * cospi_13_64);
  output[29] = dct_32_round(step[23] * cospi_3_64 + step[24] * cospi_29_64);
  output[3]  = dct_32_round(step[24] * cospi_3_64 + step[23] * -cospi_29_64);
  output[19] = dct_32_round(step[25] * cospi_19_64 + step[22] * -cospi_13_64);
  output[11] = dct_32_round(step[26] * cospi_11_64 + step[21] * -cospi_21_64);
  output[27] = dct_32_round(step[27] * cospi_27_64 + step[20] * -cospi_5_64);
  output[7]  = dct_32_round(step[28] * cospi_7_64 + step[19] * -cospi_25_64);
  output[23] = dct_32_round(step[29] * cospi_23_64 + step[18] * -cospi_9_64);
  output[15] = dct_32_round(step[30] * cospi_15_64 + step[17] * -cospi_17_64);
  output[31] = dct_32_round(step[31] * cospi_31_64 + step[16] * -cospi_1_64);
}

void vp9_short_fdct32x32_c(int16_t *input, int16_t *out, int pitch) {
  int shortpitch = pitch >> 1;
  int i, j;
  int output[1024];
  // First transform columns
  for (i = 0; i < 32; i++) {
    int temp_in[32], temp_out[32];
    for (j = 0; j < 32; j++)
      temp_in[j] = input[j * shortpitch + i] << 2;
    dct32_1d(temp_in, temp_out);
    for (j = 0; j < 32; j++)
      output[j * 32 + i] = (temp_out[j] + 1 + (temp_out[j] > 0)) >> 2;
  }

  // Then transform rows
  for (i = 0; i < 32; ++i) {
    int temp_in[32], temp_out[32];
    for (j = 0; j < 32; ++j)
      temp_in[j] = output[j + i * 32];
    dct32_1d(temp_in, temp_out);
    for (j = 0; j < 32; ++j)
      out[j + i * 32] = (temp_out[j] + 1 + (temp_out[j] < 0)) >> 2;
  }
}

