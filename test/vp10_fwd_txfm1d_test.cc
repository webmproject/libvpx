/*
 *  Copyright (c) 2015 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include "test/vp10_txfm_test.h"
#include "vp10/common/vp10_fwd_txfm1d.h"

using libvpx_test::ACMRandom;

namespace {
static int txfm_type_num = 2;
static TYPE_TXFM txfm_type_ls[2] = {TYPE_DCT, TYPE_ADST};

static int txfm_size_num = 4;
static int txfm_size_ls[4] = {4, 8, 16, 32};

static TxfmFunc fwd_txfm_func_ls[2][4] = {
    {vp10_fdct4_new, vp10_fdct8_new, vp10_fdct16_new, vp10_fdct32_new},
    {vp10_fadst4_new, vp10_fadst8_new, vp10_fadst16_new, vp10_fadst32_new}};

// the maximum stage number of fwd/inv 1d dct/adst txfm is 12
static int8_t cos_bit[12] = {14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14};
static int8_t range_bit[12] = {32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32};

TEST(vp10_fwd_txfm1d, round_shift) {
  EXPECT_EQ(round_shift(7, 1), 3);
  EXPECT_EQ(round_shift(-7, 1), -3);

  EXPECT_EQ(round_shift(7, 2), 2);
  EXPECT_EQ(round_shift(-7, 2), -2);

  EXPECT_EQ(round_shift(8, 2), 2);
  EXPECT_EQ(round_shift(-8, 2), -2);
}

TEST(vp10_fwd_txfm1d, get_max_bit) {
  int max_bit = get_max_bit(8);
  EXPECT_EQ(max_bit, 3);
}

TEST(vp10_fwd_txfm1d, half_btf) {
  int32_t max = (1 << 15) - 1;
  int32_t w0 = max;
  int32_t in0 = max;
  int32_t w1 = max;
  int32_t in1 = max;
  int32_t result_32 = half_btf(w0, in0, w1, in1, 0);
  int64_t result_64 = (int64_t)w0 * (int64_t)in0 + (int64_t)w1 * (int64_t)in1;
  EXPECT_EQ(result_32, result_64);
}

TEST(vp10_fwd_txfm1d, cospi_arr) {
  for (int i = 0; i < 7; i++) {
    for (int j = 0; j < 64; j++) {
      EXPECT_EQ(cospi_arr[i][j],
                (int32_t)round(cos(M_PI * j / 128) * (1 << (cos_bit_min + i))));
    }
  }
}

TEST(vp10_fwd_txfm1d, clamp_block) {
  int16_t block[5][5] = {{7, -5, 6, -3, 9},
                         {7, -5, 6, -3, 9},
                         {7, -5, 6, -3, 9},
                         {7, -5, 6, -3, 9},
                         {7, -5, 6, -3, 9}};

  int16_t ref_block[5][5] = {{7, -5, 6, -3, 9},
                             {7, -5, 6, -3, 9},
                             {7, -4, 2, -3, 9},
                             {7, -4, 2, -3, 9},
                             {7, -4, 2, -3, 9}};

  int row = 2;
  int col = 1;
  int block_size = 3;
  int stride = 5;
  clamp_block(block[row] + col, block_size, stride, -4, 2);
  for (int r = 0; r < stride; r++) {
    for (int c = 0; c < stride; c++) {
      EXPECT_EQ(block[r][c], ref_block[r][c]);
    }
  }
}

TEST(vp10_fwd_txfm1d, accuracy) {
  ACMRandom rnd(ACMRandom::DeterministicSeed());
  for (int si = 0; si < txfm_size_num; ++si) {
    int txfm_size = txfm_size_ls[si];
    int32_t *input = new int32_t[txfm_size];
    int32_t *output = new int32_t[txfm_size];
    double *ref_input = new double[txfm_size];
    double *ref_output = new double[txfm_size];

    for (int ti = 0; ti < txfm_type_num; ++ti) {
      TYPE_TXFM txfm_type = txfm_type_ls[ti];
      TxfmFunc fwd_txfm_func = fwd_txfm_func_ls[ti][si];
      int max_error = 7;

      const int count_test_block = 5000;
      for (int ti = 0; ti < count_test_block; ++ti) {
        for (int ni = 0; ni < txfm_size; ++ni) {
          input[ni] = rnd.Rand16() % base - rnd.Rand16() % base;
          ref_input[ni] = static_cast<double>(input[ni]);
        }

        fwd_txfm_func(input, output, cos_bit, range_bit);
        reference_hybrid_1d(ref_input, ref_output, txfm_size, txfm_type);

        for (int ni = 0; ni < txfm_size; ++ni) {
          EXPECT_LE(
              abs(output[ni] - static_cast<int32_t>(round(ref_output[ni]))),
              max_error);
        }
      }
    }

    delete[] input;
    delete[] output;
    delete[] ref_input;
    delete[] ref_output;
  }
}
}  // namespace
