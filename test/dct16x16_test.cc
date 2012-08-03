/*
 *  Copyright (c) 2012 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "third_party/googletest/src/include/gtest/gtest.h"

extern "C" {
#include "vp8/common/entropy.h"
#include "vp8/common/idct.h"
#include "vp8/encoder/dct.h"
}

#include "acm_random.h"
#include "vpx/vpx_integer.h"

using libvpx_test::ACMRandom;

namespace {

const double PI = 3.1415926535898;
void reference2_16x16_idct_2d(double *input, double *output) {
  double x;
  for (int l = 0; l < 16; ++l) {
    for (int k = 0; k < 16; ++k) {
      double s = 0;
      for (int i = 0; i < 16; ++i) {
        for (int j = 0; j < 16; ++j) {
          x=cos(PI*j*(l+0.5)/16.0)*cos(PI*i*(k+0.5)/16.0)*input[i*16+j]/256;
          if (i != 0)
            x *= sqrt(2.0);
          if (j != 0)
            x *= sqrt(2.0);
          s += x;
        }
      }
      output[k*16+l] = s;
    }
  }
}

static void butterfly_16x16_dct_1d(double input[16], double output[16]) {
  double step[16];
  double intermediate[16];
  double temp1, temp2;

  const double C1 = cos(1*PI/(double)32);
  const double C2 = cos(2*PI/(double)32);
  const double C3 = cos(3*PI/(double)32);
  const double C4 = cos(4*PI/(double)32);
  const double C5 = cos(5*PI/(double)32);
  const double C6 = cos(6*PI/(double)32);
  const double C7 = cos(7*PI/(double)32);
  const double C8 = cos(8*PI/(double)32);
  const double C9 = cos(9*PI/(double)32);
  const double C10 = cos(10*PI/(double)32);
  const double C11 = cos(11*PI/(double)32);
  const double C12 = cos(12*PI/(double)32);
  const double C13 = cos(13*PI/(double)32);
  const double C14 = cos(14*PI/(double)32);
  const double C15 = cos(15*PI/(double)32);

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

  // step 2
  output[0] = step[0] + step[7];
  output[1] = step[1] + step[6];
  output[2] = step[2] + step[5];
  output[3] = step[3] + step[4];
  output[4] = step[3] - step[4];
  output[5] = step[2] - step[5];
  output[6] = step[1] - step[6];
  output[7] = step[0] - step[7];

  temp1 = step[ 8]*C7;
  temp2 = step[15]*C9;
  output[ 8] = temp1 + temp2;

  temp1 = step[ 9]*C11;
  temp2 = step[14]*C5;
  output[ 9] = temp1 - temp2;

  temp1 = step[10]*C3;
  temp2 = step[13]*C13;
  output[10] = temp1 + temp2;

  temp1 = step[11]*C15;
  temp2 = step[12]*C1;
  output[11] = temp1 - temp2;

  temp1 = step[11]*C1;
  temp2 = step[12]*C15;
  output[12] = temp2 + temp1;

  temp1 = step[10]*C13;
  temp2 = step[13]*C3;
  output[13] = temp2 - temp1;

  temp1 = step[ 9]*C5;
  temp2 = step[14]*C11;
  output[14] = temp2 + temp1;

  temp1 = step[ 8]*C9;
  temp2 = step[15]*C7;
  output[15] = temp2 - temp1;

  // step 3
  step[ 0] = output[0] + output[3];
  step[ 1] = output[1] + output[2];
  step[ 2] = output[1] - output[2];
  step[ 3] = output[0] - output[3];

  temp1 = output[4]*C14;
  temp2 = output[7]*C2;
  step[ 4] = temp1 + temp2;

  temp1 = output[5]*C10;
  temp2 = output[6]*C6;
  step[ 5] = temp1 + temp2;

  temp1 = output[5]*C6;
  temp2 = output[6]*C10;
  step[ 6] = temp2 - temp1;

  temp1 = output[4]*C2;
  temp2 = output[7]*C14;
  step[ 7] = temp2 - temp1;

  step[ 8] = output[ 8] + output[11];
  step[ 9] = output[ 9] + output[10];
  step[10] = output[ 9] - output[10];
  step[11] = output[ 8] - output[11];

  step[12] = output[12] + output[15];
  step[13] = output[13] + output[14];
  step[14] = output[13] - output[14];
  step[15] = output[12] - output[15];

  // step 4
  output[ 0] = (step[ 0] + step[ 1]);
  output[ 8] = (step[ 0] - step[ 1]);

  temp1 = step[2]*C12;
  temp2 = step[3]*C4;
  temp1 = temp1 + temp2;
  output[ 4] = 2*(temp1*C8);

  temp1 = step[2]*C4;
  temp2 = step[3]*C12;
  temp1 = temp2 - temp1;
  output[12] = 2*(temp1*C8);

  output[ 2] = 2*((step[4] + step[ 5])*C8);
  output[14] = 2*((step[7] - step[ 6])*C8);

  temp1 = step[4] - step[5];
  temp2 = step[6] + step[7];
  output[ 6] = (temp1 + temp2);
  output[10] = (temp1 - temp2);

  intermediate[8] = step[8] + step[14];
  intermediate[9] = step[9] + step[15];

  temp1 = intermediate[8]*C12;
  temp2 = intermediate[9]*C4;
  temp1 = temp1 - temp2;
  output[3] = 2*(temp1*C8);

  temp1 = intermediate[8]*C4;
  temp2 = intermediate[9]*C12;
  temp1 = temp2 + temp1;
  output[13] = 2*(temp1*C8);

  output[ 9] = 2*((step[10] + step[11])*C8);

  intermediate[11] = step[10] - step[11];
  intermediate[12] = step[12] + step[13];
  intermediate[13] = step[12] - step[13];
  intermediate[14] = step[ 8] - step[14];
  intermediate[15] = step[ 9] - step[15];

  output[15] = (intermediate[11] + intermediate[12]);
  output[ 1] = -(intermediate[11] - intermediate[12]);

  output[ 7] = 2*(intermediate[13]*C8);

  temp1 = intermediate[14]*C12;
  temp2 = intermediate[15]*C4;
  temp1 = temp1 - temp2;
  output[11] = -2*(temp1*C8);

  temp1 = intermediate[14]*C4;
  temp2 = intermediate[15]*C12;
  temp1 = temp2 + temp1;
  output[ 5] = 2*(temp1*C8);
}

static void reference_16x16_dct_1d(double in[16], double out[16]) {
  const double kPi = 3.141592653589793238462643383279502884;
  const double kInvSqrt2 = 0.707106781186547524400844362104;
  for (int k = 0; k < 16; k++) {
    out[k] = 0.0;
    for (int n = 0; n < 16; n++)
      out[k] += in[n]*cos(kPi*(2*n+1)*k/32.0);
    if (k == 0)
      out[k] = out[k]*kInvSqrt2;
  }
}

void reference_16x16_dct_2d(int16_t input[16*16], double output[16*16]) {
  // First transform columns
  for (int i = 0; i < 16; ++i) {
    double temp_in[16], temp_out[16];
    for (int j = 0; j < 16; ++j)
      temp_in[j] = input[j*16 + i];
    butterfly_16x16_dct_1d(temp_in, temp_out);
    for (int j = 0; j < 16; ++j)
      output[j*16 + i] = temp_out[j];
  }
  // Then transform rows
  for (int i = 0; i < 16; ++i) {
    double temp_in[16], temp_out[16];
    for (int j = 0; j < 16; ++j)
      temp_in[j] = output[j + i*16];
    butterfly_16x16_dct_1d(temp_in, temp_out);
    // Scale by some magic number
    for (int j = 0; j < 16; ++j)
      output[j + i*16] = temp_out[j]/2;
  }
}


TEST(VP8Idct16x16Test, AccuracyCheck) {
  ACMRandom rnd(ACMRandom::DeterministicSeed());
  const int count_test_block = 1000;
  for (int i = 0; i < count_test_block; ++i) {
    int16_t in[256], coeff[256];
    int16_t out_c[256];
    double out_r[256];

    // Initialize a test block with input range [-255, 255].
    for (int j = 0; j < 256; ++j)
      in[j] = rnd.Rand8() - rnd.Rand8();

    reference_16x16_dct_2d(in, out_r);
    for (int j = 0; j < 256; j++)
      coeff[j] = round(out_r[j]);
    vp8_short_idct16x16_c(coeff, out_c, 32);
    for (int j = 0; j < 256; ++j) {
      const int diff = out_c[j] - in[j];
      const int error = diff * diff;
      EXPECT_GE(1, error)
          << "Error: 16x16 IDCT has error " << error
          << " at index " << j;
    }

    vp8_short_fdct16x16_c(in, out_c, 32);
    for (int j = 0; j < 256; ++j) {
      const double diff = coeff[j] - out_c[j];
      const double error = diff * diff;
      EXPECT_GE(1.0, error)
          << "Error: 16x16 FDCT has error " << error
          << " at index " << j;
    }
  }
}

TEST(VP8Fdct16x16Test, AccuracyCheck) {
  ACMRandom rnd(ACMRandom::DeterministicSeed());
  int max_error = 0;
  double total_error = 0;
  const int count_test_block = 1000;
  for (int i = 0; i < count_test_block; ++i) {
    int16_t test_input_block[256];
    int16_t test_temp_block[256];
    int16_t test_output_block[256];

    // Initialize a test block with input range [-255, 255].
    for (int j = 0; j < 256; ++j)
      test_input_block[j] = rnd.Rand8() - rnd.Rand8();

    const int pitch = 32;
    vp8_short_fdct16x16_c(test_input_block, test_temp_block, pitch);
    vp8_short_idct16x16_c(test_temp_block, test_output_block, pitch);

    for (int j = 0; j < 256; ++j) {
      const int diff = test_input_block[j] - test_output_block[j];
      const int error = diff * diff;
      if (max_error < error)
        max_error = error;
      total_error += error;
    }
  }

  EXPECT_GE(1, max_error)
      << "Error: 16x16 FDCT/IDCT has an individual roundtrip error > 1";

  EXPECT_GE(count_test_block/10, total_error)
      << "Error: 16x16 FDCT/IDCT has average roundtrip error > 1/10 per block";
}

TEST(VP8Fdct16x16Test, CoeffSizeCheck) {
  ACMRandom rnd(ACMRandom::DeterministicSeed());
  const int count_test_block = 1000;
  for (int i = 0; i < count_test_block; ++i) {
    int16_t input_block[256], input_extreme_block[256];
    int16_t output_block[256], output_extreme_block[256];

    // Initialize a test block with input range [-255, 255].
    for (int j = 0; j < 256; ++j) {
      input_block[j] = rnd.Rand8() - rnd.Rand8();
      input_extreme_block[j] = rnd.Rand8() % 2 ? 255 : -255;
    }
    if (i == 0)
      for (int j = 0; j < 256; ++j)
        input_extreme_block[j] = 255;

    const int pitch = 32;
    vp8_short_fdct16x16_c(input_block, output_block, pitch);
    vp8_short_fdct16x16_c(input_extreme_block, output_extreme_block, pitch);

    // The minimum quant value is 4.
    for (int j = 0; j < 256; ++j) {
      EXPECT_GE(4*DCT_MAX_VALUE, abs(output_block[j]))
          << "Error: 16x16 FDCT has coefficient larger than 4*DCT_MAX_VALUE";
      EXPECT_GE(4*DCT_MAX_VALUE, abs(output_extreme_block[j]))
          << "Error: 16x16 FDCT extreme has coefficient larger than 4*DCT_MAX_VALUE";
    }
  }
}
}  // namespace
