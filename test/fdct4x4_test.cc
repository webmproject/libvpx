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
#include "./vp9_rtcd.h"
}

#include "test/acm_random.h"
#include "vpx/vpx_integer.h"
#include "vpx_ports/mem.h"

using libvpx_test::ACMRandom;

namespace {
void fdct4x4(int16_t *in, int16_t *out, uint8_t* /*dst*/,
             int stride, int /*tx_type*/) {
  vp9_short_fdct4x4_c(in, out, stride);
}
void idct4x4_add(int16_t* /*in*/, int16_t *out, uint8_t *dst,
                 int stride, int /*tx_type*/) {
  vp9_short_idct4x4_add_c(out, dst, stride >> 1);
}
void fht4x4(int16_t *in, int16_t *out, uint8_t* /*dst*/,
            int stride, int tx_type) {
  vp9_short_fht4x4_c(in, out, stride >> 1, tx_type);
}
void iht4x4_add(int16_t* /*in*/, int16_t *out, uint8_t *dst,
                int stride, int tx_type) {
  vp9_short_iht4x4_add_c(out, dst, stride >> 1, tx_type);
}

class FwdTrans4x4Test : public ::testing::TestWithParam<int> {
 public:
  virtual ~FwdTrans4x4Test() {}
  virtual void SetUp() {
    tx_type_ = GetParam();
    if (tx_type_ == 0) {
      fwd_txfm_ = fdct4x4;
      inv_txfm_ = idct4x4_add;
    } else {
      fwd_txfm_ = fht4x4;
      inv_txfm_ = iht4x4_add;
    }
  }

 protected:
  void RunFwdTxfm(int16_t *in, int16_t *out, uint8_t *dst,
                  int stride, int tx_type) {
    (*fwd_txfm_)(in, out, dst, stride, tx_type);
  }

  void RunInvTxfm(int16_t *in, int16_t *out, uint8_t *dst,
                  int stride, int tx_type) {
    (*inv_txfm_)(in, out, dst, stride, tx_type);
  }

  int tx_type_;
  void (*fwd_txfm_)(int16_t *in, int16_t *out, uint8_t *dst,
                   int stride, int tx_type);
  void (*inv_txfm_)(int16_t *in, int16_t *out, uint8_t *dst,
                   int stride, int tx_type);
};

TEST_P(FwdTrans4x4Test, SignBiasCheck) {
  ACMRandom rnd(ACMRandom::DeterministicSeed());
  DECLARE_ALIGNED_ARRAY(16, int16_t, test_input_block, 16);
  DECLARE_ALIGNED_ARRAY(16, int16_t, test_output_block, 16);
  const int pitch = 8;
  int count_sign_block[16][2];
  const int count_test_block = 1000000;

  memset(count_sign_block, 0, sizeof(count_sign_block));
  for (int i = 0; i < count_test_block; ++i) {
    // Initialize a test block with input range [-255, 255].
    for (int j = 0; j < 16; ++j)
      test_input_block[j] = rnd.Rand8() - rnd.Rand8();

    RunFwdTxfm(test_input_block, test_output_block, NULL, pitch, tx_type_);

    for (int j = 0; j < 16; ++j) {
      if (test_output_block[j] < 0)
        ++count_sign_block[j][0];
      else if (test_output_block[j] > 0)
        ++count_sign_block[j][1];
    }
  }

  for (int j = 0; j < 16; ++j) {
    const bool bias_acceptable = (abs(count_sign_block[j][0] -
                                      count_sign_block[j][1]) < 10000);
    EXPECT_TRUE(bias_acceptable)
        << "Error: 4x4 FDCT/FHT has a sign bias > 1%"
        << " for input range [-255, 255] at index " << j
        << " tx_type " << tx_type_;
  }

  memset(count_sign_block, 0, sizeof(count_sign_block));
  for (int i = 0; i < count_test_block; ++i) {
    // Initialize a test block with input range [-15, 15].
    for (int j = 0; j < 16; ++j)
      test_input_block[j] = (rnd.Rand8() >> 4) - (rnd.Rand8() >> 4);

    RunFwdTxfm(test_input_block, test_output_block, NULL, pitch, tx_type_);

    for (int j = 0; j < 16; ++j) {
      if (test_output_block[j] < 0)
        ++count_sign_block[j][0];
      else if (test_output_block[j] > 0)
        ++count_sign_block[j][1];
    }
  }

  for (int j = 0; j < 16; ++j) {
    const bool bias_acceptable = (abs(count_sign_block[j][0] -
                                      count_sign_block[j][1]) < 100000);
    EXPECT_TRUE(bias_acceptable)
        << "Error: 4x4 FDCT/FHT has a sign bias > 10%"
        << " for input range [-15, 15] at index " << j;
  }
}

TEST_P(FwdTrans4x4Test, RoundTripErrorCheck) {
  ACMRandom rnd(ACMRandom::DeterministicSeed());

  int max_error = 0;
  int total_error = 0;
  const int count_test_block = 1000000;
  for (int i = 0; i < count_test_block; ++i) {
    DECLARE_ALIGNED_ARRAY(16, int16_t, test_input_block, 16);
    DECLARE_ALIGNED_ARRAY(16, int16_t, test_temp_block, 16);
    DECLARE_ALIGNED_ARRAY(16, uint8_t, dst, 16);
    DECLARE_ALIGNED_ARRAY(16, uint8_t, src, 16);

    for (int j = 0; j < 16; ++j) {
      src[j] = rnd.Rand8();
      dst[j] = rnd.Rand8();
    }
    // Initialize a test block with input range [-255, 255].
    for (int j = 0; j < 16; ++j)
      test_input_block[j] = src[j] - dst[j];

    const int pitch = 8;
    RunFwdTxfm(test_input_block, test_temp_block, dst, pitch, tx_type_);

    for (int j = 0; j < 16; ++j) {
        if (test_temp_block[j] > 0) {
          test_temp_block[j] += 2;
          test_temp_block[j] /= 4;
          test_temp_block[j] *= 4;
        } else {
          test_temp_block[j] -= 2;
          test_temp_block[j] /= 4;
          test_temp_block[j] *= 4;
        }
    }

    // inverse transform and reconstruct the pixel block
    RunInvTxfm(test_input_block, test_temp_block, dst, pitch, tx_type_);

    for (int j = 0; j < 16; ++j) {
      const int diff = dst[j] - src[j];
      const int error = diff * diff;
      if (max_error < error)
        max_error = error;
      total_error += error;
    }
  }
  EXPECT_GE(1, max_error)
      << "Error: FDCT/IDCT or FHT/IHT has an individual roundtrip error > 1";

  EXPECT_GE(count_test_block, total_error)
      << "Error: FDCT/IDCT or FHT/IHT has average "
      << "roundtrip error > 1 per block";
}

INSTANTIATE_TEST_CASE_P(VP9, FwdTrans4x4Test, ::testing::Range(0, 4));
}  // namespace
