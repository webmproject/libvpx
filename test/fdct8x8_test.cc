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
#include "test/clear_system_state.h"
#include "test/register_state_check.h"
#include "vpx_ports/mem.h"

extern "C" {
#include "./vp9_rtcd.h"
void vp9_short_idct8x8_add_c(int16_t *input, uint8_t *output, int pitch);
}

#include "test/acm_random.h"
#include "vpx/vpx_integer.h"

using libvpx_test::ACMRandom;

namespace {
void fdct8x8(int16_t *in, int16_t *out, uint8_t* /*dst*/,
             int stride, int /*tx_type*/) {
  vp9_short_fdct8x8_c(in, out, stride);
}
void idct8x8_add(int16_t* /*in*/, int16_t *out, uint8_t *dst,
                 int stride, int /*tx_type*/) {
  vp9_short_idct8x8_add_c(out, dst, stride >> 1);
}
void fht8x8(int16_t *in, int16_t *out, uint8_t* /*dst*/,
            int stride, int tx_type) {
  // TODO(jingning): need to refactor this to test both _c and _sse2 functions,
  // when we have all inverse dct functions done sse2.
#if HAVE_SSE2
  vp9_short_fht8x8_sse2(in, out, stride >> 1, tx_type);
#else
  vp9_short_fht8x8_c(in, out, stride >> 1, tx_type);
#endif
}
void iht8x8_add(int16_t* /*in*/, int16_t *out, uint8_t *dst,
                int stride, int tx_type) {
  vp9_short_iht8x8_add_c(out, dst, stride >> 1, tx_type);
}

class FwdTrans8x8Test : public ::testing::TestWithParam<int> {
 public:
  virtual ~FwdTrans8x8Test() {}
  virtual void SetUp() {
    tx_type_ = GetParam();
    if (tx_type_ == 0) {
      fwd_txfm = fdct8x8;
      inv_txfm = idct8x8_add;
    } else {
      fwd_txfm = fht8x8;
      inv_txfm = iht8x8_add;
    }
  }
  virtual void TearDown() { libvpx_test::ClearSystemState(); }

 protected:
  void RunFwdTxfm(int16_t *in, int16_t *out, uint8_t *dst,
                  int stride, int tx_type) {
    (*fwd_txfm)(in, out, dst, stride, tx_type);
  }
  void RunInvTxfm(int16_t *in, int16_t *out, uint8_t *dst,
                  int stride, int tx_type) {
    (*inv_txfm)(in, out, dst, stride, tx_type);
  }

  int tx_type_;
  void (*fwd_txfm)(int16_t*, int16_t*, uint8_t*, int, int);
  void (*inv_txfm)(int16_t*, int16_t*, uint8_t*, int, int);
};

TEST_P(FwdTrans8x8Test, SignBiasCheck) {
  ACMRandom rnd(ACMRandom::DeterministicSeed());
  DECLARE_ALIGNED_ARRAY(16, int16_t, test_input_block, 64);
  DECLARE_ALIGNED_ARRAY(16, int16_t, test_output_block, 64);
  const int pitch = 16;
  int count_sign_block[64][2];
  const int count_test_block = 100000;

  memset(count_sign_block, 0, sizeof(count_sign_block));

  for (int i = 0; i < count_test_block; ++i) {
    // Initialize a test block with input range [-255, 255].
    for (int j = 0; j < 64; ++j)
      test_input_block[j] = rnd.Rand8() - rnd.Rand8();
    REGISTER_STATE_CHECK(
        RunFwdTxfm(test_input_block, test_output_block,
                   NULL, pitch, tx_type_));

    for (int j = 0; j < 64; ++j) {
      if (test_output_block[j] < 0)
        ++count_sign_block[j][0];
      else if (test_output_block[j] > 0)
        ++count_sign_block[j][1];
    }
  }

  for (int j = 0; j < 64; ++j) {
    const int diff = abs(count_sign_block[j][0] - count_sign_block[j][1]);
    const int max_diff = 1125;
    EXPECT_LT(diff, max_diff)
        << "Error: 8x8 FDCT/FHT has a sign bias > "
        << 1. * max_diff / count_test_block * 100 << "%"
        << " for input range [-255, 255] at index " << j
        << " count0: " << count_sign_block[j][0]
        << " count1: " << count_sign_block[j][1]
        << " diff: " << diff;
  }

  memset(count_sign_block, 0, sizeof(count_sign_block));

  for (int i = 0; i < count_test_block; ++i) {
    // Initialize a test block with input range [-15, 15].
    for (int j = 0; j < 64; ++j)
      test_input_block[j] = (rnd.Rand8() >> 4) - (rnd.Rand8() >> 4);
    REGISTER_STATE_CHECK(
        RunFwdTxfm(test_input_block, test_output_block,
                   NULL, pitch, tx_type_));

    for (int j = 0; j < 64; ++j) {
      if (test_output_block[j] < 0)
        ++count_sign_block[j][0];
      else if (test_output_block[j] > 0)
        ++count_sign_block[j][1];
    }
  }

  for (int j = 0; j < 64; ++j) {
    const int diff = abs(count_sign_block[j][0] - count_sign_block[j][1]);
    const int max_diff = 10000;
    EXPECT_LT(diff, max_diff)
        << "Error: 4x4 FDCT/FHT has a sign bias > "
        << 1. * max_diff / count_test_block * 100 << "%"
        << " for input range [-15, 15] at index " << j
        << " count0: " << count_sign_block[j][0]
        << " count1: " << count_sign_block[j][1]
        << " diff: " << diff;
  }
}

TEST_P(FwdTrans8x8Test, RoundTripErrorCheck) {
  ACMRandom rnd(ACMRandom::DeterministicSeed());
  int max_error = 0;
  int total_error = 0;
  const int count_test_block = 100000;
  for (int i = 0; i < count_test_block; ++i) {
    DECLARE_ALIGNED_ARRAY(16, int16_t, test_input_block, 64);
    DECLARE_ALIGNED_ARRAY(16, int16_t, test_temp_block, 64);
    DECLARE_ALIGNED_ARRAY(16, uint8_t, dst, 64);
    DECLARE_ALIGNED_ARRAY(16, uint8_t, src, 64);

    for (int j = 0; j < 64; ++j) {
      src[j] = rnd.Rand8();
      dst[j] = rnd.Rand8();
    }
    // Initialize a test block with input range [-255, 255].
    for (int j = 0; j < 64; ++j)
      test_input_block[j] = src[j] - dst[j];

    const int pitch = 16;
    REGISTER_STATE_CHECK(
        RunFwdTxfm(test_input_block, test_temp_block,
                   dst, pitch, tx_type_));
    for (int j = 0; j < 64; ++j) {
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
    REGISTER_STATE_CHECK(
        RunInvTxfm(test_input_block, test_temp_block,
                   dst, pitch, tx_type_));

    for (int j = 0; j < 64; ++j) {
      const int diff = dst[j] - src[j];
      const int error = diff * diff;
      if (max_error < error)
        max_error = error;
      total_error += error;
    }
  }

  EXPECT_GE(1, max_error)
    << "Error: 8x8 FDCT/IDCT or FHT/IHT has an individual roundtrip error > 1";

  EXPECT_GE(count_test_block/5, total_error)
    << "Error: 8x8 FDCT/IDCT or FHT/IHT has average roundtrip "
        "error > 1/5 per block";
}

TEST_P(FwdTrans8x8Test, ExtremalCheck) {
  ACMRandom rnd(ACMRandom::DeterministicSeed());
  int max_error = 0;
  int total_error = 0;
  const int count_test_block = 100000;
  for (int i = 0; i < count_test_block; ++i) {
    DECLARE_ALIGNED_ARRAY(16, int16_t, test_input_block, 64);
    DECLARE_ALIGNED_ARRAY(16, int16_t, test_temp_block, 64);
    DECLARE_ALIGNED_ARRAY(16, uint8_t, dst, 64);
    DECLARE_ALIGNED_ARRAY(16, uint8_t, src, 64);

    for (int j = 0; j < 64; ++j) {
      src[j] = rnd.Rand8() % 2 ? 255 : 0;
      dst[j] = src[j] > 0 ? 0 : 255;
    }
    // Initialize a test block with input range [-255, 255].
    for (int j = 0; j < 64; ++j)
      test_input_block[j] = src[j] - dst[j];

    const int pitch = 16;
    REGISTER_STATE_CHECK(
        RunFwdTxfm(test_input_block, test_temp_block,
                   dst, pitch, tx_type_));
    REGISTER_STATE_CHECK(
        RunInvTxfm(test_input_block, test_temp_block,
                   dst, pitch, tx_type_));

    for (int j = 0; j < 64; ++j) {
      const int diff = dst[j] - src[j];
      const int error = diff * diff;
      if (max_error < error)
        max_error = error;
      total_error += error;
    }

    EXPECT_GE(1, max_error)
        << "Error: Extremal 8x8 FDCT/IDCT or FHT/IHT has an"
        << " individual roundtrip error > 1";

    EXPECT_GE(count_test_block/5, total_error)
        << "Error: Extremal 8x8 FDCT/IDCT or FHT/IHT has average"
        << " roundtrip error > 1/5 per block";
  }
}

INSTANTIATE_TEST_CASE_P(VP9, FwdTrans8x8Test, ::testing::Range(0, 4));
}  // namespace
