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

#define TEST_INT_32x32_DCT 1

#if !TEST_INT_32x32_DCT

static void dct32_1d(double *input, double *output, int stride) {
  static const double C1 = 0.998795456205;  // cos(pi * 1 / 64)
  static const double C2 = 0.995184726672;  // cos(pi * 2 / 64)
  static const double C3 = 0.989176509965;  // cos(pi * 3 / 64)
  static const double C4 = 0.980785280403;  // cos(pi * 4 / 64)
  static const double C5 = 0.970031253195;  // cos(pi * 5 / 64)
  static const double C6 = 0.956940335732;  // cos(pi * 6 / 64)
  static const double C7 = 0.941544065183;  // cos(pi * 7 / 64)
  static const double C8 = 0.923879532511;  // cos(pi * 8 / 64)
  static const double C9 = 0.903989293123;  // cos(pi * 9 / 64)
  static const double C10 = 0.881921264348;  // cos(pi * 10 / 64)
  static const double C11 = 0.857728610000;  // cos(pi * 11 / 64)
  static const double C12 = 0.831469612303;  // cos(pi * 12 / 64)
  static const double C13 = 0.803207531481;  // cos(pi * 13 / 64)
  static const double C14 = 0.773010453363;  // cos(pi * 14 / 64)
  static const double C15 = 0.740951125355;  // cos(pi * 15 / 64)
  static const double C16 = 0.707106781187;  // cos(pi * 16 / 64)
  static const double C17 = 0.671558954847;  // cos(pi * 17 / 64)
  static const double C18 = 0.634393284164;  // cos(pi * 18 / 64)
  static const double C19 = 0.595699304492;  // cos(pi * 19 / 64)
  static const double C20 = 0.555570233020;  // cos(pi * 20 / 64)
  static const double C21 = 0.514102744193;  // cos(pi * 21 / 64)
  static const double C22 = 0.471396736826;  // cos(pi * 22 / 64)
  static const double C23 = 0.427555093430;  // cos(pi * 23 / 64)
  static const double C24 = 0.382683432365;  // cos(pi * 24 / 64)
  static const double C25 = 0.336889853392;  // cos(pi * 25 / 64)
  static const double C26 = 0.290284677254;  // cos(pi * 26 / 64)
  static const double C27 = 0.242980179903;  // cos(pi * 27 / 64)
  static const double C28 = 0.195090322016;  // cos(pi * 28 / 64)
  static const double C29 = 0.146730474455;  // cos(pi * 29 / 64)
  static const double C30 = 0.098017140330;  // cos(pi * 30 / 64)
  static const double C31 = 0.049067674327;  // cos(pi * 31 / 64)

  double step[32];

  // Stage 1
  step[0] = input[stride*0] + input[stride*(32 - 1)];
  step[1] = input[stride*1] + input[stride*(32 - 2)];
  step[2] = input[stride*2] + input[stride*(32 - 3)];
  step[3] = input[stride*3] + input[stride*(32 - 4)];
  step[4] = input[stride*4] + input[stride*(32 - 5)];
  step[5] = input[stride*5] + input[stride*(32 - 6)];
  step[6] = input[stride*6] + input[stride*(32 - 7)];
  step[7] = input[stride*7] + input[stride*(32 - 8)];
  step[8] = input[stride*8] + input[stride*(32 - 9)];
  step[9] = input[stride*9] + input[stride*(32 - 10)];
  step[10] = input[stride*10] + input[stride*(32 - 11)];
  step[11] = input[stride*11] + input[stride*(32 - 12)];
  step[12] = input[stride*12] + input[stride*(32 - 13)];
  step[13] = input[stride*13] + input[stride*(32 - 14)];
  step[14] = input[stride*14] + input[stride*(32 - 15)];
  step[15] = input[stride*15] + input[stride*(32 - 16)];
  step[16] = -input[stride*16] + input[stride*(32 - 17)];
  step[17] = -input[stride*17] + input[stride*(32 - 18)];
  step[18] = -input[stride*18] + input[stride*(32 - 19)];
  step[19] = -input[stride*19] + input[stride*(32 - 20)];
  step[20] = -input[stride*20] + input[stride*(32 - 21)];
  step[21] = -input[stride*21] + input[stride*(32 - 22)];
  step[22] = -input[stride*22] + input[stride*(32 - 23)];
  step[23] = -input[stride*23] + input[stride*(32 - 24)];
  step[24] = -input[stride*24] + input[stride*(32 - 25)];
  step[25] = -input[stride*25] + input[stride*(32 - 26)];
  step[26] = -input[stride*26] + input[stride*(32 - 27)];
  step[27] = -input[stride*27] + input[stride*(32 - 28)];
  step[28] = -input[stride*28] + input[stride*(32 - 29)];
  step[29] = -input[stride*29] + input[stride*(32 - 30)];
  step[30] = -input[stride*30] + input[stride*(32 - 31)];
  step[31] = -input[stride*31] + input[stride*(32 - 32)];

  // Stage 2
  output[stride*0] = step[0] + step[16 - 1];
  output[stride*1] = step[1] + step[16 - 2];
  output[stride*2] = step[2] + step[16 - 3];
  output[stride*3] = step[3] + step[16 - 4];
  output[stride*4] = step[4] + step[16 - 5];
  output[stride*5] = step[5] + step[16 - 6];
  output[stride*6] = step[6] + step[16 - 7];
  output[stride*7] = step[7] + step[16 - 8];
  output[stride*8] = -step[8] + step[16 - 9];
  output[stride*9] = -step[9] + step[16 - 10];
  output[stride*10] = -step[10] + step[16 - 11];
  output[stride*11] = -step[11] + step[16 - 12];
  output[stride*12] = -step[12] + step[16 - 13];
  output[stride*13] = -step[13] + step[16 - 14];
  output[stride*14] = -step[14] + step[16 - 15];
  output[stride*15] = -step[15] + step[16 - 16];

  output[stride*16] = step[16];
  output[stride*17] = step[17];
  output[stride*18] = step[18];
  output[stride*19] = step[19];

  output[stride*20] = (-step[20] + step[27])*C16;
  output[stride*21] = (-step[21] + step[26])*C16;
  output[stride*22] = (-step[22] + step[25])*C16;
  output[stride*23] = (-step[23] + step[24])*C16;

  output[stride*24] = (step[24] + step[23])*C16;
  output[stride*25] = (step[25] + step[22])*C16;
  output[stride*26] = (step[26] + step[21])*C16;
  output[stride*27] = (step[27] + step[20])*C16;

  output[stride*28] = step[28];
  output[stride*29] = step[29];
  output[stride*30] = step[30];
  output[stride*31] = step[31];

  // Stage 3
  step[0] = output[stride*0] + output[stride*(8 - 1)];
  step[1] = output[stride*1] + output[stride*(8 - 2)];
  step[2] = output[stride*2] + output[stride*(8 - 3)];
  step[3] = output[stride*3] + output[stride*(8 - 4)];
  step[4] = -output[stride*4] + output[stride*(8 - 5)];
  step[5] = -output[stride*5] + output[stride*(8 - 6)];
  step[6] = -output[stride*6] + output[stride*(8 - 7)];
  step[7] = -output[stride*7] + output[stride*(8 - 8)];
  step[8] = output[stride*8];
  step[9] = output[stride*9];
  step[10] = (-output[stride*10] + output[stride*13])*C16;
  step[11] = (-output[stride*11] + output[stride*12])*C16;
  step[12] = (output[stride*12] + output[stride*11])*C16;
  step[13] = (output[stride*13] + output[stride*10])*C16;
  step[14] = output[stride*14];
  step[15] = output[stride*15];

  step[16] = output[stride*16] + output[stride*23];
  step[17] = output[stride*17] + output[stride*22];
  step[18] = output[stride*18] + output[stride*21];
  step[19] = output[stride*19] + output[stride*20];
  step[20] = -output[stride*20] + output[stride*19];
  step[21] = -output[stride*21] + output[stride*18];
  step[22] = -output[stride*22] + output[stride*17];
  step[23] = -output[stride*23] + output[stride*16];
  step[24] = -output[stride*24] + output[stride*31];
  step[25] = -output[stride*25] + output[stride*30];
  step[26] = -output[stride*26] + output[stride*29];
  step[27] = -output[stride*27] + output[stride*28];
  step[28] = output[stride*28] + output[stride*27];
  step[29] = output[stride*29] + output[stride*26];
  step[30] = output[stride*30] + output[stride*25];
  step[31] = output[stride*31] + output[stride*24];

  // Stage 4
  output[stride*0] = step[0] + step[3];
  output[stride*1] = step[1] + step[2];
  output[stride*2] = -step[2] + step[1];
  output[stride*3] = -step[3] + step[0];
  output[stride*4] = step[4];
  output[stride*5] = (-step[5] + step[6])*C16;
  output[stride*6] = (step[6] + step[5])*C16;
  output[stride*7] = step[7];
  output[stride*8] = step[8] + step[11];
  output[stride*9] = step[9] + step[10];
  output[stride*10] = -step[10] + step[9];
  output[stride*11] = -step[11] + step[8];
  output[stride*12] = -step[12] + step[15];
  output[stride*13] = -step[13] + step[14];
  output[stride*14] = step[14] + step[13];
  output[stride*15] = step[15] + step[12];

  output[stride*16] = step[16];
  output[stride*17] = step[17];
  output[stride*18] = step[18]*-C8 + step[29]*C24;
  output[stride*19] = step[19]*-C8 + step[28]*C24;
  output[stride*20] = step[20]*-C24 + step[27]*-C8;
  output[stride*21] = step[21]*-C24 + step[26]*-C8;
  output[stride*22] = step[22];
  output[stride*23] = step[23];
  output[stride*24] = step[24];
  output[stride*25] = step[25];
  output[stride*26] = step[26]*C24 + step[21]*-C8;
  output[stride*27] = step[27]*C24 + step[20]*-C8;
  output[stride*28] = step[28]*C8 + step[19]*C24;
  output[stride*29] = step[29]*C8 + step[18]*C24;
  output[stride*30] = step[30];
  output[stride*31] = step[31];

  // Stage 5
  step[0] = (output[stride*0] + output[stride*1]) * C16;
  step[1] = (-output[stride*1] + output[stride*0]) * C16;
  step[2] = output[stride*2]*C24 + output[stride*3] * C8;
  step[3] = output[stride*3]*C24 - output[stride*2] * C8;
  step[4] = output[stride*4] + output[stride*5];
  step[5] = -output[stride*5] + output[stride*4];
  step[6] = -output[stride*6] + output[stride*7];
  step[7] = output[stride*7] + output[stride*6];
  step[8] = output[stride*8];
  step[9] = output[stride*9]*-C8 + output[stride*14]*C24;
  step[10] = output[stride*10]*-C24 + output[stride*13]*-C8;
  step[11] = output[stride*11];
  step[12] = output[stride*12];
  step[13] = output[stride*13]*C24 + output[stride*10]*-C8;
  step[14] = output[stride*14]*C8 + output[stride*9]*C24;
  step[15] = output[stride*15];

  step[16] = output[stride*16] + output[stride*19];
  step[17] = output[stride*17] + output[stride*18];
  step[18] = -output[stride*18] + output[stride*17];
  step[19] = -output[stride*19] + output[stride*16];
  step[20] = -output[stride*20] + output[stride*23];
  step[21] = -output[stride*21] + output[stride*22];
  step[22] = output[stride*22] + output[stride*21];
  step[23] = output[stride*23] + output[stride*20];
  step[24] = output[stride*24] + output[stride*27];
  step[25] = output[stride*25] + output[stride*26];
  step[26] = -output[stride*26] + output[stride*25];
  step[27] = -output[stride*27] + output[stride*24];
  step[28] = -output[stride*28] + output[stride*31];
  step[29] = -output[stride*29] + output[stride*30];
  step[30] = output[stride*30] + output[stride*29];
  step[31] = output[stride*31] + output[stride*28];

  // Stage 6
  output[stride*0] = step[0];
  output[stride*1] = step[1];
  output[stride*2] = step[2];
  output[stride*3] = step[3];
  output[stride*4] = step[4]*C28 + step[7]*C4;
  output[stride*5] = step[5]*C12 + step[6]*C20;
  output[stride*6] = step[6]*C12 + step[5]*-C20;
  output[stride*7] = step[7]*C28 + step[4]*-C4;
  output[stride*8] = step[8] + step[9];
  output[stride*9] = -step[9] + step[8];
  output[stride*10] = -step[10] + step[11];
  output[stride*11] = step[11] + step[10];
  output[stride*12] = step[12] + step[13];
  output[stride*13] = -step[13] + step[12];
  output[stride*14] = -step[14] + step[15];
  output[stride*15] = step[15] + step[14];

  output[stride*16] = step[16];
  output[stride*17] = step[17]*-C4 + step[30]*C28;
  output[stride*18] = step[18]*-C28 + step[29]*-C4;
  output[stride*19] = step[19];
  output[stride*20] = step[20];
  output[stride*21] = step[21]*-C20 + step[26]*C12;
  output[stride*22] = step[22]*-C12 + step[25]*-C20;
  output[stride*23] = step[23];
  output[stride*24] = step[24];
  output[stride*25] = step[25]*C12 + step[22]*-C20;
  output[stride*26] = step[26]*C20 + step[21]*C12;
  output[stride*27] = step[27];
  output[stride*28] = step[28];
  output[stride*29] = step[29]*C28 + step[18]*-C4;
  output[stride*30] = step[30]*C4 + step[17]*C28;
  output[stride*31] = step[31];

  // Stage 7
  step[0] = output[stride*0];
  step[1] = output[stride*1];
  step[2] = output[stride*2];
  step[3] = output[stride*3];
  step[4] = output[stride*4];
  step[5] = output[stride*5];
  step[6] = output[stride*6];
  step[7] = output[stride*7];
  step[8] = output[stride*8]*C30 + output[stride*15]*C2;
  step[9] = output[stride*9]*C14 + output[stride*14]*C18;
  step[10] = output[stride*10]*C22 + output[stride*13]*C10;
  step[11] = output[stride*11]*C6 + output[stride*12]*C26;
  step[12] = output[stride*12]*C6 + output[stride*11]*-C26;
  step[13] = output[stride*13]*C22 + output[stride*10]*-C10;
  step[14] = output[stride*14]*C14 + output[stride*9]*-C18;
  step[15] = output[stride*15]*C30 + output[stride*8]*-C2;

  step[16] = output[stride*16] + output[stride*17];
  step[17] = -output[stride*17] + output[stride*16];
  step[18] = -output[stride*18] + output[stride*19];
  step[19] = output[stride*19] + output[stride*18];
  step[20] = output[stride*20] + output[stride*21];
  step[21] = -output[stride*21] + output[stride*20];
  step[22] = -output[stride*22] + output[stride*23];
  step[23] = output[stride*23] + output[stride*22];
  step[24] = output[stride*24] + output[stride*25];
  step[25] = -output[stride*25] + output[stride*24];
  step[26] = -output[stride*26] + output[stride*27];
  step[27] = output[stride*27] + output[stride*26];
  step[28] = output[stride*28] + output[stride*29];
  step[29] = -output[stride*29] + output[stride*28];
  step[30] = -output[stride*30] + output[stride*31];
  step[31] = output[stride*31] + output[stride*30];

  // Final stage --- outputs indices are bit-reversed.
  output[stride*0] = step[0];
  output[stride*16] = step[1];
  output[stride*8] = step[2];
  output[stride*24] = step[3];
  output[stride*4] = step[4];
  output[stride*20] = step[5];
  output[stride*12] = step[6];
  output[stride*28] = step[7];
  output[stride*2] = step[8];
  output[stride*18] = step[9];
  output[stride*10] = step[10];
  output[stride*26] = step[11];
  output[stride*6] = step[12];
  output[stride*22] = step[13];
  output[stride*14] = step[14];
  output[stride*30] = step[15];

  output[stride*1] = step[16]*C31 + step[31]*C1;
  output[stride*17] = step[17]*C15 + step[30]*C17;
  output[stride*9] = step[18]*C23 + step[29]*C9;
  output[stride*25] = step[19]*C7 + step[28]*C25;
  output[stride*5] = step[20]*C27 + step[27]*C5;
  output[stride*21] = step[21]*C11 + step[26]*C21;
  output[stride*13] = step[22]*C19 + step[25]*C13;
  output[stride*29] = step[23]*C3 + step[24]*C29;
  output[stride*3] = step[24]*C3 + step[23]*-C29;
  output[stride*19] = step[25]*C19 + step[22]*-C13;
  output[stride*11] = step[26]*C11 + step[21]*-C21;
  output[stride*27] = step[27]*C27 + step[20]*-C5;
  output[stride*7] = step[28]*C7 + step[19]*-C25;
  output[stride*23] = step[29]*C23 + step[18]*-C9;
  output[stride*15] = step[30]*C15 + step[17]*-C17;
  output[stride*31] = step[31]*C31 + step[16]*-C1;
}

void vp9_short_fdct32x32_c(int16_t *input, int16_t *out, int pitch) {
  vp9_clear_system_state();  // Make it simd safe : __asm emms;
  {
    int shortpitch = pitch >> 1;
    int i, j;
    double output[1024];
    // First transform columns
    for (i = 0; i < 32; i++) {
      double temp_in[32], temp_out[32];
      for (j = 0; j < 32; j++)
        temp_in[j] = input[j*shortpitch + i];
      dct32_1d(temp_in, temp_out, 1);
      for (j = 0; j < 32; j++)
        output[j*32 + i] = temp_out[j];
    }
    // Then transform rows
    for (i = 0; i < 32; ++i) {
      double temp_in[32], temp_out[32];
      for (j = 0; j < 32; ++j)
        temp_in[j] = output[j + i*32];
      dct32_1d(temp_in, temp_out, 1);
      for (j = 0; j < 32; ++j)
        output[j + i*32] = temp_out[j];
    }
    // Scale by some magic number
    for (i = 0; i < 1024; i++) {
      out[i] = (short)round(output[i]/4);
    }
  }

  vp9_clear_system_state();  // Make it simd safe : __asm emms;
}

#else

#define RIGHT_SHIFT 13
#define ROUNDING (1 << (RIGHT_SHIFT - 1))

static void dct32_1d(int *input, int *output, int last_shift_bits) {
  static const int16_t C1 = 8182;    // 2^13
  static const int16_t C2 = 8153;
  static const int16_t C3 = 8103;
  static const int16_t C4 = 8035;
  static const int16_t C5 = 7946;
  static const int16_t C6 = 7839;
  static const int16_t C7 = 7713;
  static const int16_t C8 = 7568;
  static const int16_t C9 = 7405;
  static const int16_t C10 = 7225;
  static const int16_t C11 = 7027;
  static const int16_t C12 = 6811;
  static const int16_t C13 = 6580;
  static const int16_t C14 = 6333;
  static const int16_t C15 = 6070;
  static const int16_t C16 = 5793;
  static const int16_t C17 = 5501;
  static const int16_t C18 = 5197;
  static const int16_t C19 = 4880;
  static const int16_t C20 = 4551;
  static const int16_t C21 = 4212;
  static const int16_t C22 = 3862;
  static const int16_t C23 = 3503;
  static const int16_t C24 = 3135;
  static const int16_t C25 = 2760;
  static const int16_t C26 = 2378;
  static const int16_t C27 = 1990;
  static const int16_t C28 = 1598;
  static const int16_t C29 = 1202;
  static const int16_t C30 = 803;
  static const int16_t C31 = 402;

  int step[32];

  int last_rounding = 0;
  int final_shift = RIGHT_SHIFT;
  int final_rounding = 0;

  if (last_shift_bits > 0)
    last_rounding = 1 << (last_shift_bits - 1);

  final_shift += last_shift_bits;
  if (final_shift > 0)
    final_rounding = 1 << (final_shift - 1);

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

  output[20] = ((-step[20] + step[27]) * C16 + ROUNDING) >> RIGHT_SHIFT;
  output[21] = ((-step[21] + step[26]) * C16 + ROUNDING) >> RIGHT_SHIFT;
  output[22] = ((-step[22] + step[25]) * C16 + ROUNDING) >> RIGHT_SHIFT;
  output[23] = ((-step[23] + step[24]) * C16 + ROUNDING) >> RIGHT_SHIFT;

  output[24] = ((step[24] + step[23]) * C16 + ROUNDING) >> RIGHT_SHIFT;
  output[25] = ((step[25] + step[22]) * C16 + ROUNDING) >> RIGHT_SHIFT;
  output[26] = ((step[26] + step[21]) * C16 + ROUNDING) >> RIGHT_SHIFT;
  output[27] = ((step[27] + step[20]) * C16 + ROUNDING) >> RIGHT_SHIFT;

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
  step[10] = ((-output[10] + output[13]) * C16 + ROUNDING) >> RIGHT_SHIFT;
  step[11] = ((-output[11] + output[12]) * C16 + ROUNDING) >> RIGHT_SHIFT;
  step[12] = ((output[12] + output[11]) * C16 + ROUNDING) >> RIGHT_SHIFT;
  step[13] = ((output[13] + output[10]) * C16 + ROUNDING) >> RIGHT_SHIFT;
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
  output[5] = ((-step[5] + step[6]) * C16 + ROUNDING) >> RIGHT_SHIFT;
  output[6] = ((step[6] + step[5]) * C16 + ROUNDING) >> RIGHT_SHIFT;
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
  output[18] = (step[18] * -C8 + step[29] * C24 + ROUNDING) >> RIGHT_SHIFT;
  output[19] = (step[19] * -C8 + step[28] * C24 + ROUNDING) >> RIGHT_SHIFT;
  output[20] = (step[20] * -C24 + step[27] * -C8 + ROUNDING) >> RIGHT_SHIFT;
  output[21] = (step[21] * -C24 + step[26] * -C8 + ROUNDING) >> RIGHT_SHIFT;
  output[22] = step[22];
  output[23] = step[23];
  output[24] = step[24];
  output[25] = step[25];
  output[26] = (step[26] * C24 + step[21] * -C8 + ROUNDING) >> RIGHT_SHIFT;
  output[27] = (step[27] * C24 + step[20] * -C8 + ROUNDING) >> RIGHT_SHIFT;
  output[28] = (step[28] * C8 + step[19] * C24 + ROUNDING) >> RIGHT_SHIFT;
  output[29] = (step[29] * C8 + step[18] * C24 + ROUNDING) >> RIGHT_SHIFT;
  output[30] = step[30];
  output[31] = step[31];

  // Stage 5
  step[0] = ((output[0] + output[1]) * C16 + ROUNDING) >> RIGHT_SHIFT;
  step[1] = ((-output[1] + output[0]) * C16 + ROUNDING) >> RIGHT_SHIFT;
  step[2] = (output[2] * C24 + output[3] * C8 + ROUNDING) >> RIGHT_SHIFT;
  step[3] = (output[3] * C24 - output[2] * C8 + ROUNDING) >> RIGHT_SHIFT;
  step[4] = output[4] + output[5];
  step[5] = -output[5] + output[4];
  step[6] = -output[6] + output[7];
  step[7] = output[7] + output[6];
  step[8] = output[8];
  step[9] = (output[9] * -C8 + output[14] * C24 + ROUNDING) >> RIGHT_SHIFT;
  step[10] = (output[10] * -C24 + output[13] * -C8 + ROUNDING) >> RIGHT_SHIFT;
  step[11] = output[11];
  step[12] = output[12];
  step[13] = (output[13] * C24 + output[10] * -C8 + ROUNDING) >> RIGHT_SHIFT;
  step[14] = (output[14] * C8 + output[9] * C24 + ROUNDING) >> RIGHT_SHIFT;
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
  output[4] = (step[4] * C28 + step[7] * C4 + ROUNDING) >> RIGHT_SHIFT;
  output[5] = (step[5] * C12 + step[6] * C20 + ROUNDING) >> RIGHT_SHIFT;
  output[6] = (step[6] * C12 + step[5] * -C20 + ROUNDING) >> RIGHT_SHIFT;
  output[7] = (step[7] * C28 + step[4] * -C4 + ROUNDING) >> RIGHT_SHIFT;
  output[8] = step[8] + step[9];
  output[9] = -step[9] + step[8];
  output[10] = -step[10] + step[11];
  output[11] = step[11] + step[10];
  output[12] = step[12] + step[13];
  output[13] = -step[13] + step[12];
  output[14] = -step[14] + step[15];
  output[15] = step[15] + step[14];

  output[16] = step[16];
  output[17] = (step[17] * -C4 + step[30] * C28 + ROUNDING) >> RIGHT_SHIFT;
  output[18] = (step[18] * -C28 + step[29] * -C4 + ROUNDING) >> RIGHT_SHIFT;
  output[19] = step[19];
  output[20] = step[20];
  output[21] = (step[21] * -C20 + step[26] * C12 + ROUNDING) >> RIGHT_SHIFT;
  output[22] = (step[22] * -C12 + step[25] * -C20 + ROUNDING) >> RIGHT_SHIFT;
  output[23] = step[23];
  output[24] = step[24];
  output[25] = (step[25] * C12 + step[22] * -C20 + ROUNDING) >> RIGHT_SHIFT;
  output[26] = (step[26] * C20 + step[21] * C12 + ROUNDING) >> RIGHT_SHIFT;
  output[27] = step[27];
  output[28] = step[28];
  output[29] = (step[29] * C28 + step[18] * -C4 + ROUNDING) >> RIGHT_SHIFT;
  output[30] = (step[30] * C4 + step[17] * C28 + ROUNDING) >> RIGHT_SHIFT;
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
  step[8] = (output[8] * C30 + output[15] * C2 + ROUNDING) >> RIGHT_SHIFT;
  step[9] = (output[9] * C14 + output[14] * C18 + ROUNDING) >> RIGHT_SHIFT;
  step[10] = (output[10] * C22 + output[13] * C10 + ROUNDING) >> RIGHT_SHIFT;
  step[11] = (output[11] * C6 + output[12] * C26 + ROUNDING) >> RIGHT_SHIFT;
  step[12] = (output[12] * C6 + output[11] * -C26 + ROUNDING) >> RIGHT_SHIFT;
  step[13] = (output[13] * C22 + output[10] * -C10 + ROUNDING) >> RIGHT_SHIFT;
  step[14] = (output[14] * C14 + output[9] * -C18 + ROUNDING) >> RIGHT_SHIFT;
  step[15] = (output[15] * C30 + output[8] * -C2 + ROUNDING) >> RIGHT_SHIFT;

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
  output[0] = (step[0] + last_rounding) >> last_shift_bits;
  output[16] = (step[1] + last_rounding) >> last_shift_bits;
  output[8] = (step[2] + last_rounding) >> last_shift_bits;
  output[24] = (step[3] + last_rounding) >> last_shift_bits;
  output[4] = (step[4] + last_rounding) >> last_shift_bits;
  output[20] = (step[5] + last_rounding) >> last_shift_bits;
  output[12] = (step[6] + last_rounding) >> last_shift_bits;
  output[28] = (step[7] + last_rounding) >> last_shift_bits;
  output[2] = (step[8] + last_rounding) >> last_shift_bits;
  output[18] = (step[9] + last_rounding) >> last_shift_bits;
  output[10] = (step[10] + last_rounding) >> last_shift_bits;
  output[26] = (step[11] + last_rounding) >> last_shift_bits;
  output[6] = (step[12] + last_rounding) >> last_shift_bits;
  output[22] = (step[13] + last_rounding) >> last_shift_bits;
  output[14] = (step[14] + last_rounding) >> last_shift_bits;
  output[30] = (step[15] + last_rounding) >> last_shift_bits;

  output[1] = (step[16] * C31 + step[31] * C1 + final_rounding) >> final_shift;
  output[17] = (step[17] * C15 + step[30] * C17 + final_rounding)
      >> final_shift;
  output[9] = (step[18] * C23 + step[29] * C9 + final_rounding) >> final_shift;
  output[25] = (step[19] * C7 + step[28] * C25 + final_rounding) >> final_shift;
  output[5] = (step[20] * C27 + step[27] * C5 + final_rounding) >> final_shift;
  output[21] = (step[21] * C11 + step[26] * C21 + final_rounding)
      >> final_shift;
  output[13] = (step[22] * C19 + step[25] * C13 + final_rounding)
      >> final_shift;
  output[29] = (step[23] * C3 + step[24] * C29 + final_rounding) >> final_shift;
  output[3] = (step[24] * C3 + step[23] * -C29 + final_rounding) >> final_shift;
  output[19] = (step[25] * C19 + step[22] * -C13 + final_rounding)
      >> final_shift;
  output[11] = (step[26] * C11 + step[21] * -C21 + final_rounding)
      >> final_shift;
  output[27] = (step[27] * C27 + step[20] * -C5 + final_rounding)
      >> final_shift;
  output[7] = (step[28] * C7 + step[19] * -C25 + final_rounding) >> final_shift;
  output[23] = (step[29] * C23 + step[18] * -C9 + final_rounding)
      >> final_shift;
  output[15] = (step[30] * C15 + step[17] * -C17 + final_rounding)
      >> final_shift;
  output[31] = (step[31] * C31 + step[16] * -C1 + final_rounding)
      >> final_shift;

  // Clamp to fit 16-bit.
  if (last_shift_bits > 0) {
    int i;

    for (i = 0; i < 32; i++)
      if (output[i] < -32768)
        output[i] = -32768;
      else if (output[i] > 32767)
        output[i] = 32767;
  }
}
#undef RIGHT_SHIFT
#undef ROUNDING

void vp9_short_fdct32x32_c(int16_t *input, int16_t *out, int pitch) {
  int shortpitch = pitch >> 1;
  int i, j;
  int output[1024];
  // First transform columns
  for (i = 0; i < 32; i++) {
    int temp_in[32], temp_out[32];
    for (j = 0; j < 32; j++)
      temp_in[j] = input[j * shortpitch + i];
    dct32_1d(temp_in, temp_out, 0);
    for (j = 0; j < 32; j++)
      output[j * 32 + i] = temp_out[j];
  }

  // Then transform rows
  for (i = 0; i < 32; ++i) {
    int temp_in[32], temp_out[32];
    for (j = 0; j < 32; ++j)
      temp_in[j] = output[j + i * 32];
    dct32_1d(temp_in, temp_out, 2);
    for (j = 0; j < 32; ++j)
      out[j + i * 32] = temp_out[j];
  }
}

#endif
