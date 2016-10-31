/*
 *  Copyright (c) 2013 The WebM project authors. All Rights Reserved.
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

#include "./vp9_rtcd.h"
#include "./vpx_dsp_rtcd.h"
#include "test/acm_random.h"
#include "test/clear_system_state.h"
#include "test/register_state_check.h"
#include "test/util.h"
#include "vp9/common/vp9_blockd.h"
#include "vp9/common/vp9_scan.h"
#include "vpx/vpx_integer.h"

using libvpx_test::ACMRandom;

namespace {
typedef void (*FwdTxfmFunc)(const int16_t *in, tran_low_t *out, int stride);
typedef void (*InvTxfmFunc)(const tran_low_t *in, uint8_t *out, int stride);
typedef std::tr1::tuple<FwdTxfmFunc, InvTxfmFunc, InvTxfmFunc, TX_SIZE, int>
    PartialInvTxfmParam;
const int kMaxNumCoeffs = 1024;
class PartialIDctTest : public ::testing::TestWithParam<PartialInvTxfmParam> {
 public:
  virtual ~PartialIDctTest() {}
  virtual void SetUp() {
    ftxfm_ = GET_PARAM(0);
    full_itxfm_ = GET_PARAM(1);
    partial_itxfm_ = GET_PARAM(2);
    tx_size_ = GET_PARAM(3);
    last_nonzero_ = GET_PARAM(4);

    switch (tx_size_) {
      case TX_4X4: size_ = 4; break;
      case TX_8X8: size_ = 8; break;
      case TX_16X16: size_ = 16; break;
      case TX_32X32: size_ = 32; break;
      default: FAIL() << "Wrong Size!"; break;
    }
    block_size_ = size_ * size_;

    input_block_ = reinterpret_cast<tran_low_t *>(
        vpx_memalign(16, sizeof(*input_block_) * block_size_));
    output_block_ = reinterpret_cast<uint8_t *>(
        vpx_memalign(16, sizeof(*output_block_) * block_size_));
    output_block_ref_ = reinterpret_cast<uint8_t *>(
        vpx_memalign(16, sizeof(*output_block_ref_) * block_size_));
  }

  virtual void TearDown() {
    vpx_free(input_block_);
    input_block_ = NULL;
    vpx_free(output_block_);
    output_block_ = NULL;
    vpx_free(output_block_ref_);
    output_block_ref_ = NULL;
    libvpx_test::ClearSystemState();
  }

 protected:
  int last_nonzero_;
  TX_SIZE tx_size_;
  tran_low_t *input_block_;
  uint8_t *output_block_;
  uint8_t *output_block_ref_;
  int size_;
  int block_size_;
  FwdTxfmFunc ftxfm_;
  InvTxfmFunc full_itxfm_;
  InvTxfmFunc partial_itxfm_;
};

TEST_P(PartialIDctTest, RunQuantCheck) {
  ACMRandom rnd(ACMRandom::DeterministicSeed());

  const int count_test_block = 1000;

  DECLARE_ALIGNED(16, int16_t, input_extreme_block[kMaxNumCoeffs]);
  DECLARE_ALIGNED(16, tran_low_t, output_ref_block[kMaxNumCoeffs]);

  for (int i = 0; i < count_test_block; ++i) {
    // clear out destination buffer
    memset(input_block_, 0, sizeof(*input_block_) * block_size_);
    memset(output_block_, 0, sizeof(*output_block_) * block_size_);
    memset(output_block_ref_, 0, sizeof(*output_block_ref_) * block_size_);

    ACMRandom rnd(ACMRandom::DeterministicSeed());

    for (int i = 0; i < count_test_block; ++i) {
      // Initialize a test block with input range [-255, 255].
      if (i == 0) {
        for (int j = 0; j < block_size_; ++j) input_extreme_block[j] = 255;
      } else if (i == 1) {
        for (int j = 0; j < block_size_; ++j) input_extreme_block[j] = -255;
      } else {
        for (int j = 0; j < block_size_; ++j) {
          input_extreme_block[j] = rnd.Rand8() % 2 ? 255 : -255;
        }
      }

      ftxfm_(input_extreme_block, output_ref_block, size_);

      // quantization with maximum allowed step sizes
      input_block_[0] = (output_ref_block[0] / 1336) * 1336;
      for (int j = 1; j < last_nonzero_; ++j) {
        input_block_[vp9_default_scan_orders[tx_size_].scan[j]] =
            (output_ref_block[j] / 1828) * 1828;
      }
    }

    ASM_REGISTER_STATE_CHECK(
        full_itxfm_(input_block_, output_block_ref_, size_));
    ASM_REGISTER_STATE_CHECK(
        partial_itxfm_(input_block_, output_block_, size_));

    ASSERT_EQ(0, memcmp(output_block_ref_, output_block_,
                        sizeof(*output_block_) * block_size_))
        << "Error: partial inverse transform produces different results";
  }
}

TEST_P(PartialIDctTest, ResultsMatch) {
  ACMRandom rnd(ACMRandom::DeterministicSeed());
  const int count_test_block = 1000;
  const int max_coeff = 32766 / 4;
  for (int i = 0; i < count_test_block; ++i) {
    // clear out destination buffer
    memset(input_block_, 0, sizeof(*input_block_) * block_size_);
    memset(output_block_, 0, sizeof(*output_block_) * block_size_);
    memset(output_block_ref_, 0, sizeof(*output_block_ref_) * block_size_);
    int max_energy_leftover = max_coeff * max_coeff;
    for (int j = 0; j < last_nonzero_; ++j) {
      int16_t coef = static_cast<int16_t>(sqrt(1.0 * max_energy_leftover) *
                                          (rnd.Rand16() - 32768) / 65536);
      max_energy_leftover -= coef * coef;
      if (max_energy_leftover < 0) {
        max_energy_leftover = 0;
        coef = 0;
      }
      input_block_[vp9_default_scan_orders[tx_size_].scan[j]] = coef;
    }

    ASM_REGISTER_STATE_CHECK(
        full_itxfm_(input_block_, output_block_ref_, size_));
    ASM_REGISTER_STATE_CHECK(
        partial_itxfm_(input_block_, output_block_, size_));

    ASSERT_EQ(0, memcmp(output_block_ref_, output_block_,
                        sizeof(*output_block_) * block_size_))
        << "Error: partial inverse transform produces different results";
  }
}

TEST_P(PartialIDctTest, AddOutputBlock) {
  ACMRandom rnd(ACMRandom::DeterministicSeed());
  const int count_test_block = 10;
  for (int i = 0; i < count_test_block; ++i) {
    memset(input_block_, 0, sizeof(*input_block_) * block_size_);
    for (int j = 0; j < last_nonzero_; ++j) {
      input_block_[vp9_default_scan_orders[tx_size_].scan[j]] = 10;
    }

    for (int j = 0; j < block_size_; ++j) {
      output_block_[j] = output_block_ref_[j] = rnd.Rand8();
    }

    ASM_REGISTER_STATE_CHECK(
        full_itxfm_(input_block_, output_block_ref_, size_));
    ASM_REGISTER_STATE_CHECK(
        partial_itxfm_(input_block_, output_block_, size_));

    ASSERT_EQ(0, memcmp(output_block_ref_, output_block_,
                        sizeof(*output_block_) * block_size_))
        << "Error: Transform results are not correctly added to output.";
  }
}
using std::tr1::make_tuple;

INSTANTIATE_TEST_CASE_P(
    C, PartialIDctTest,
    ::testing::Values(make_tuple(&vpx_fdct32x32_c, &vpx_idct32x32_1024_add_c,
                                 &vpx_idct32x32_1024_add_c, TX_32X32, 1024),
                      make_tuple(&vpx_fdct32x32_c, &vpx_idct32x32_1024_add_c,
                                 &vpx_idct32x32_135_add_c, TX_32X32, 135),
                      make_tuple(&vpx_fdct32x32_c, &vpx_idct32x32_1024_add_c,
                                 &vpx_idct32x32_34_add_c, TX_32X32, 34),
                      make_tuple(&vpx_fdct32x32_c, &vpx_idct32x32_1024_add_c,
                                 &vpx_idct32x32_1_add_c, TX_32X32, 1),
                      make_tuple(&vpx_fdct16x16_c, &vpx_idct16x16_256_add_c,
                                 &vpx_idct16x16_256_add_c, TX_16X16, 256),
                      make_tuple(&vpx_fdct16x16_c, &vpx_idct16x16_256_add_c,
                                 &vpx_idct16x16_10_add_c, TX_16X16, 10),
                      make_tuple(&vpx_fdct16x16_c, &vpx_idct16x16_256_add_c,
                                 &vpx_idct16x16_1_add_c, TX_16X16, 1),
                      make_tuple(&vpx_fdct8x8_c, &vpx_idct8x8_64_add_c,
                                 &vpx_idct8x8_64_add_c, TX_8X8, 64),
                      make_tuple(&vpx_fdct8x8_c, &vpx_idct8x8_64_add_c,
                                 &vpx_idct8x8_12_add_c, TX_8X8, 12),
                      make_tuple(&vpx_fdct8x8_c, &vpx_idct8x8_64_add_c,
                                 &vpx_idct8x8_1_add_c, TX_8X8, 1),
                      make_tuple(&vpx_fdct4x4_c, &vpx_idct4x4_16_add_c,
                                 &vpx_idct4x4_16_add_c, TX_4X4, 16),
                      make_tuple(&vpx_fdct4x4_c, &vpx_idct4x4_16_add_c,
                                 &vpx_idct4x4_1_add_c, TX_4X4, 1)));

#if HAVE_NEON && !CONFIG_EMULATE_HARDWARE
#if CONFIG_VP9_HIGHBITDEPTH
INSTANTIATE_TEST_CASE_P(
    NEON, PartialIDctTest,
    ::testing::Values(make_tuple(&vpx_fdct32x32_c, &vpx_idct32x32_1024_add_c,
                                 &vpx_idct32x32_1_add_neon, TX_32X32, 1),
                      make_tuple(&vpx_fdct16x16_c, &vpx_idct16x16_256_add_c,
                                 &vpx_idct16x16_1_add_neon, TX_16X16, 1),
                      make_tuple(&vpx_fdct8x8_c, &vpx_idct8x8_64_add_c,
                                 &vpx_idct8x8_1_add_neon, TX_8X8, 1),
                      make_tuple(&vpx_fdct4x4_c, &vpx_idct4x4_16_add_c,
                                 &vpx_idct4x4_1_add_neon, TX_4X4, 1)));
#else   // !CONFIG_VP9_HIGHBITDEPTH
// 32x32_135_ is implemented using the 1024 version.
INSTANTIATE_TEST_CASE_P(
    NEON, PartialIDctTest,
    ::testing::Values(make_tuple(&vpx_fdct32x32_c, &vpx_idct32x32_1024_add_c,
                                 &vpx_idct32x32_1024_add_neon, TX_32X32, 1024),
                      make_tuple(&vpx_fdct32x32_c, &vpx_idct32x32_1024_add_c,
                                 &vpx_idct32x32_1024_add_neon, TX_32X32, 135),
                      make_tuple(&vpx_fdct32x32_c, &vpx_idct32x32_1024_add_c,
                                 &vpx_idct32x32_34_add_neon, TX_32X32, 34),
                      make_tuple(&vpx_fdct32x32_c, &vpx_idct32x32_1024_add_c,
                                 &vpx_idct32x32_1_add_neon, TX_32X32, 1),
                      make_tuple(&vpx_fdct16x16_c, &vpx_idct16x16_256_add_c,
                                 &vpx_idct16x16_256_add_neon, TX_16X16, 256),
                      make_tuple(&vpx_fdct16x16_c, &vpx_idct16x16_256_add_c,
                                 &vpx_idct16x16_10_add_neon, TX_16X16, 10),
                      make_tuple(&vpx_fdct16x16_c, &vpx_idct16x16_256_add_c,
                                 &vpx_idct16x16_1_add_neon, TX_16X16, 1),
                      make_tuple(&vpx_fdct8x8_c, &vpx_idct8x8_64_add_c,
                                 &vpx_idct8x8_64_add_neon, TX_8X8, 64),
                      make_tuple(&vpx_fdct8x8_c, &vpx_idct8x8_64_add_c,
                                 &vpx_idct8x8_12_add_neon, TX_8X8, 12),
                      make_tuple(&vpx_fdct8x8_c, &vpx_idct8x8_64_add_c,
                                 &vpx_idct8x8_1_add_neon, TX_8X8, 1),
                      make_tuple(&vpx_fdct4x4_c, &vpx_idct4x4_16_add_c,
                                 &vpx_idct4x4_16_add_neon, TX_4X4, 16),
                      make_tuple(&vpx_fdct4x4_c, &vpx_idct4x4_16_add_c,
                                 &vpx_idct4x4_1_add_neon, TX_4X4, 1)));
#endif  // CONFIG_VP9_HIGHBITDEPTH
#endif  // HAVE_NEON && !CONFIG_EMULATE_HARDWARE

#if HAVE_SSE2 && !CONFIG_VP9_HIGHBITDEPTH && !CONFIG_EMULATE_HARDWARE
// 32x32_135_ is implemented using the 1024 version.
INSTANTIATE_TEST_CASE_P(
    SSE2, PartialIDctTest,
    ::testing::Values(make_tuple(&vpx_fdct32x32_c, &vpx_idct32x32_1024_add_c,
                                 &vpx_idct32x32_1024_add_sse2, TX_32X32, 1024),
                      make_tuple(&vpx_fdct32x32_c, &vpx_idct32x32_1024_add_c,
                                 &vpx_idct32x32_1024_add_sse2, TX_32X32, 135),
                      make_tuple(&vpx_fdct32x32_c, &vpx_idct32x32_1024_add_c,
                                 &vpx_idct32x32_34_add_sse2, TX_32X32, 34),
                      make_tuple(&vpx_fdct32x32_c, &vpx_idct32x32_1024_add_c,
                                 &vpx_idct32x32_1_add_sse2, TX_32X32, 1),
                      make_tuple(&vpx_fdct16x16_c, &vpx_idct16x16_256_add_c,
                                 &vpx_idct16x16_256_add_sse2, TX_16X16, 256),
                      make_tuple(&vpx_fdct16x16_c, &vpx_idct16x16_256_add_c,
                                 &vpx_idct16x16_10_add_sse2, TX_16X16, 10),
                      make_tuple(&vpx_fdct16x16_c, &vpx_idct16x16_256_add_c,
                                 &vpx_idct16x16_1_add_sse2, TX_16X16, 1),
                      make_tuple(&vpx_fdct8x8_c, &vpx_idct8x8_64_add_c,
                                 &vpx_idct8x8_64_add_sse2, TX_8X8, 64),
                      make_tuple(&vpx_fdct8x8_c, &vpx_idct8x8_64_add_c,
                                 &vpx_idct8x8_12_add_sse2, TX_8X8, 12),
                      make_tuple(&vpx_fdct8x8_c, &vpx_idct8x8_64_add_c,
                                 &vpx_idct8x8_1_add_sse2, TX_8X8, 1),
                      make_tuple(&vpx_fdct4x4_c, &vpx_idct4x4_16_add_c,
                                 &vpx_idct4x4_16_add_sse2, TX_4X4, 16),
                      make_tuple(&vpx_fdct4x4_c, &vpx_idct4x4_16_add_c,
                                 &vpx_idct4x4_1_add_sse2, TX_4X4, 1)));
#endif  // HAVE_SSE2 && !CONFIG_VP9_HIGHBITDEPTH && !CONFIG_EMULATE_HARDWARE

#if HAVE_SSSE3 && ARCH_X86_64 && !CONFIG_VP9_HIGHBITDEPTH && \
    !CONFIG_EMULATE_HARDWARE
INSTANTIATE_TEST_CASE_P(
    SSSE3_64, PartialIDctTest,
    ::testing::Values(make_tuple(&vpx_fdct32x32_c, &vpx_idct32x32_1024_add_c,
                                 &vpx_idct32x32_1024_add_ssse3, TX_32X32, 1024),
                      make_tuple(&vpx_fdct32x32_c, &vpx_idct32x32_1024_add_c,
                                 &vpx_idct32x32_135_add_ssse3, TX_32X32, 135),
                      make_tuple(&vpx_fdct32x32_c, &vpx_idct32x32_1024_add_c,
                                 &vpx_idct32x32_34_add_ssse3, TX_32X32, 34),
                      make_tuple(&vpx_fdct8x8_c, &vpx_idct8x8_64_add_c,
                                 &vpx_idct8x8_64_add_ssse3, TX_8X8, 64),
                      make_tuple(&vpx_fdct8x8_c, &vpx_idct8x8_64_add_c,
                                 &vpx_idct8x8_12_add_ssse3, TX_8X8, 12)));
#endif  // HAVE_SSSE3 && ARCH_X86_64 && !CONFIG_VP9_HIGHBITDEPTH &&
        // !CONFIG_EMULATE_HARDWARE

#if HAVE_MSA && !CONFIG_VP9_HIGHBITDEPTH && !CONFIG_EMULATE_HARDWARE
// 32x32_135_ is implemented using the 1024 version.
INSTANTIATE_TEST_CASE_P(
    MSA, PartialIDctTest,
    ::testing::Values(make_tuple(&vpx_fdct32x32_c, &vpx_idct32x32_1024_add_c,
                                 &vpx_idct32x32_1024_add_msa, TX_32X32, 1024),
                      make_tuple(&vpx_fdct32x32_c, &vpx_idct32x32_1024_add_c,
                                 &vpx_idct32x32_1024_add_msa, TX_32X32, 135),
                      make_tuple(&vpx_fdct32x32_c, &vpx_idct32x32_1024_add_c,
                                 &vpx_idct32x32_34_add_msa, TX_32X32, 34),
                      make_tuple(&vpx_fdct32x32_c, &vpx_idct32x32_1024_add_c,
                                 &vpx_idct32x32_1_add_msa, TX_32X32, 1),
                      make_tuple(&vpx_fdct16x16_c, &vpx_idct16x16_256_add_c,
                                 &vpx_idct16x16_256_add_msa, TX_16X16, 256),
                      make_tuple(&vpx_fdct16x16_c, &vpx_idct16x16_256_add_c,
                                 &vpx_idct16x16_10_add_msa, TX_16X16, 10),
                      make_tuple(&vpx_fdct16x16_c, &vpx_idct16x16_256_add_c,
                                 &vpx_idct16x16_1_add_msa, TX_16X16, 1),
                      make_tuple(&vpx_fdct8x8_c, &vpx_idct8x8_64_add_c,
                                 &vpx_idct8x8_64_add_msa, TX_8X8, 64),
                      make_tuple(&vpx_fdct8x8_c, &vpx_idct8x8_64_add_c,
                                 &vpx_idct8x8_12_add_msa, TX_8X8, 10),
                      make_tuple(&vpx_fdct8x8_c, &vpx_idct8x8_64_add_c,
                                 &vpx_idct8x8_1_add_msa, TX_8X8, 1),
                      make_tuple(&vpx_fdct4x4_c, &vpx_idct4x4_16_add_c,
                                 &vpx_idct4x4_16_add_msa, TX_4X4, 16),
                      make_tuple(&vpx_fdct4x4_c, &vpx_idct4x4_16_add_c,
                                 &vpx_idct4x4_1_add_msa, TX_4X4, 1)));
#endif  // HAVE_MSA && !CONFIG_VP9_HIGHBITDEPTH && !CONFIG_EMULATE_HARDWARE

}  // namespace
