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

#include <limits>

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
#include "vpx_ports/vpx_timer.h"

using libvpx_test::ACMRandom;

namespace {

typedef void (*FwdTxfmFunc)(const int16_t *in, tran_low_t *out, int stride);

#if CONFIG_VP9_HIGHBITDEPTH
typedef uint16_t Pixel;
typedef void (*InvTxfmFunc)(const tran_low_t *in, uint8_t *out, int stride,
                            int bd);
#else   // !CONFIG_VP9_HIGHBITDEPTH
typedef uint8_t Pixel;
typedef void (*InvTxfmFunc)(const tran_low_t *in, uint8_t *out, int stride);
#endif  // CONFIG_VP9_HIGHBITDEPTH

typedef std::tr1::tuple<FwdTxfmFunc, InvTxfmFunc, InvTxfmFunc, TX_SIZE, int,
                        int>
    PartialInvTxfmParam;
const int kMaxNumCoeffs = 1024;
const int kCountTestBlock = 1000;

// https://bugs.chromium.org/p/webm/issues/detail?id=1332
// The functions specified do not pass with INT16_MIN/MAX. They fail at the
// value specified, but pass when 1 is added/subtracted.
int16_t MaxSupportedCoeff(InvTxfmFunc a) {
#if HAVE_SSSE3 && ARCH_X86_64 && !CONFIG_VP9_HIGHBITDEPTH && \
    !CONFIG_EMULATE_HARDWARE
  if (a == vpx_idct8x8_64_add_ssse3 || a == vpx_idct8x8_12_add_ssse3) {
    return 23625 - 1;
  }
#else
  (void)a;
#endif
  return std::numeric_limits<int16_t>::max();
}

int16_t MinSupportedCoeff(InvTxfmFunc a) {
  (void)a;
#if HAVE_SSSE3 && ARCH_X86_64 && !CONFIG_VP9_HIGHBITDEPTH && \
    !CONFIG_EMULATE_HARDWARE
  if (a == vpx_idct8x8_64_add_ssse3 || a == vpx_idct8x8_12_add_ssse3) {
    return -23625 + 1;
  }
#endif  // !CONFIG_EMULATE_HARDWARE
  return std::numeric_limits<int16_t>::min();
}

class PartialIDctTest : public ::testing::TestWithParam<PartialInvTxfmParam> {
 public:
  virtual ~PartialIDctTest() {}
  virtual void SetUp() {
    rnd_.Reset(ACMRandom::DeterministicSeed());
    ftxfm_ = GET_PARAM(0);
    full_itxfm_ = GET_PARAM(1);
    partial_itxfm_ = GET_PARAM(2);
    tx_size_ = GET_PARAM(3);
    last_nonzero_ = GET_PARAM(4);
    bit_depth_ = GET_PARAM(5);
    mask_ = (1 << bit_depth_) - 1;

    switch (tx_size_) {
      case TX_4X4: size_ = 4; break;
      case TX_8X8: size_ = 8; break;
      case TX_16X16: size_ = 16; break;
      case TX_32X32: size_ = 32; break;
      default: FAIL() << "Wrong Size!"; break;
    }

    // Randomize stride_ to a value less than or equal to 1024
    stride_ = rnd_(1024) + 1;
    if (stride_ < size_) {
      stride_ = size_;
    }
    // Align stride_ to 16 if it's bigger than 16.
    if (stride_ > 16) {
      stride_ &= ~15;
    }

    input_block_size_ = size_ * size_;
    output_block_size_ = size_ * stride_;

    input_block_ = reinterpret_cast<tran_low_t *>(
        vpx_memalign(16, sizeof(*input_block_) * input_block_size_));
    output_block_ = reinterpret_cast<Pixel *>(
        vpx_memalign(16, sizeof(*output_block_) * output_block_size_));
    output_block_ref_ = reinterpret_cast<Pixel *>(
        vpx_memalign(16, sizeof(*output_block_ref_) * output_block_size_));
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

  void InitMem() {
    memset(input_block_, 0, sizeof(*input_block_) * input_block_size_);
    for (int j = 0; j < output_block_size_; ++j) {
      output_block_[j] = output_block_ref_[j] = rnd_.Rand16() & mask_;
    }
  }

  void InitInput() {
    const int max_coeff = 32766 / 4;
    int max_energy_leftover = max_coeff * max_coeff;
    for (int j = 0; j < last_nonzero_; ++j) {
      int16_t coeff = static_cast<int16_t>(sqrt(1.0 * max_energy_leftover) *
                                           (rnd_.Rand16() - 32768) / 65536);
      max_energy_leftover -= coeff * coeff;
      if (max_energy_leftover < 0) {
        max_energy_leftover = 0;
        coeff = 0;
      }
      input_block_[vp9_default_scan_orders[tx_size_].scan[j]] = coeff;
    }
  }

  void Exec(InvTxfmFunc func, void *out) {
#if CONFIG_VP9_HIGHBITDEPTH
    func(input_block_, CONVERT_TO_BYTEPTR(out), stride_, bit_depth_);
#else
    func(input_block_, reinterpret_cast<uint8_t *>(out), stride_);
#endif
  }

 protected:
  int last_nonzero_;
  TX_SIZE tx_size_;
  tran_low_t *input_block_;
  Pixel *output_block_;
  Pixel *output_block_ref_;
  int size_;
  int stride_;
  int input_block_size_;
  int output_block_size_;
  int bit_depth_;
  int mask_;
  FwdTxfmFunc ftxfm_;
  InvTxfmFunc full_itxfm_;
  InvTxfmFunc partial_itxfm_;
  ACMRandom rnd_;
};

TEST_P(PartialIDctTest, RunQuantCheck) {
  DECLARE_ALIGNED(16, int16_t, input_extreme_block[kMaxNumCoeffs]);
  DECLARE_ALIGNED(16, tran_low_t, output_ref_block[kMaxNumCoeffs]);

  InitMem();
  for (int i = 0; i < kCountTestBlock * kCountTestBlock; ++i) {
    // Initialize a test block with input range [-mask_, mask_].
    if (i == 0) {
      for (int k = 0; k < input_block_size_; ++k) {
        input_extreme_block[k] = mask_;
      }
    } else if (i == 1) {
      for (int k = 0; k < input_block_size_; ++k) {
        input_extreme_block[k] = -mask_;
      }
    } else {
      for (int k = 0; k < input_block_size_; ++k) {
        input_extreme_block[k] = rnd_.Rand8() % 2 ? mask_ : -mask_;
      }
    }

    ftxfm_(input_extreme_block, output_ref_block, size_);

    // quantization with minimum allowed step sizes
    input_block_[0] = (output_ref_block[0] / 4) * 4;
    for (int k = 1; k < last_nonzero_; ++k) {
      const int pos = vp9_default_scan_orders[tx_size_].scan[k];
      input_block_[pos] = (output_ref_block[pos] / 4) * 4;
    }

    ASM_REGISTER_STATE_CHECK(Exec(full_itxfm_, output_block_ref_));
    ASM_REGISTER_STATE_CHECK(Exec(partial_itxfm_, output_block_));
    ASSERT_EQ(0, memcmp(output_block_ref_, output_block_,
                        sizeof(*output_block_) * output_block_size_))
        << "Error: partial inverse transform produces different results";
  }
}

TEST_P(PartialIDctTest, ResultsMatch) {
  for (int i = 0; i < kCountTestBlock; ++i) {
    InitMem();
    InitInput();

    ASM_REGISTER_STATE_CHECK(Exec(full_itxfm_, output_block_ref_));
    ASM_REGISTER_STATE_CHECK(Exec(partial_itxfm_, output_block_));

    ASSERT_EQ(0, memcmp(output_block_ref_, output_block_,
                        sizeof(*output_block_) * output_block_size_))
        << "Error: partial inverse transform produces different results";
  }
}

TEST_P(PartialIDctTest, AddOutputBlock) {
  for (int i = 0; i < kCountTestBlock; ++i) {
    InitMem();
    for (int j = 0; j < last_nonzero_; ++j) {
      input_block_[vp9_default_scan_orders[tx_size_].scan[j]] = 10;
    }

    ASM_REGISTER_STATE_CHECK(Exec(full_itxfm_, output_block_ref_));
    ASM_REGISTER_STATE_CHECK(Exec(partial_itxfm_, output_block_));

    ASSERT_EQ(0, memcmp(output_block_ref_, output_block_,
                        sizeof(*output_block_) * output_block_size_))
        << "Error: Transform results are not correctly added to output.";
  }
}

TEST_P(PartialIDctTest, SingleExtremeCoeff) {
  const int16_t max_coeff = MaxSupportedCoeff(partial_itxfm_);
  const int16_t min_coeff = MinSupportedCoeff(partial_itxfm_);
  for (int i = 0; i < last_nonzero_; ++i) {
    memset(input_block_, 0, sizeof(*input_block_) * input_block_size_);
    // Run once for min and once for max.
    for (int j = 0; j < 2; ++j) {
      const int coeff = j ? min_coeff : max_coeff;

      memset(output_block_, 0, sizeof(*output_block_) * output_block_size_);
      memset(output_block_ref_, 0,
             sizeof(*output_block_ref_) * output_block_size_);
      input_block_[vp9_default_scan_orders[tx_size_].scan[i]] = coeff;

      ASM_REGISTER_STATE_CHECK(Exec(full_itxfm_, output_block_ref_));
      ASM_REGISTER_STATE_CHECK(Exec(partial_itxfm_, output_block_));

      ASSERT_EQ(0, memcmp(output_block_ref_, output_block_,
                          sizeof(*output_block_) * output_block_size_))
          << "Error: Fails with single coeff of " << coeff << " at " << i
          << ".";
    }
  }
}

TEST_P(PartialIDctTest, DISABLED_Speed) {
  // Keep runtime stable with transform size.
  const int kCountSpeedTestBlock = 500000000 / input_block_size_;
  InitMem();
  InitInput();

  for (int i = 0; i < kCountSpeedTestBlock; ++i) {
    ASM_REGISTER_STATE_CHECK(Exec(full_itxfm_, output_block_ref_));
  }
  vpx_usec_timer timer;
  vpx_usec_timer_start(&timer);
  for (int i = 0; i < kCountSpeedTestBlock; ++i) {
    Exec(partial_itxfm_, output_block_);
  }
  libvpx_test::ClearSystemState();
  vpx_usec_timer_mark(&timer);
  const int elapsed_time =
      static_cast<int>(vpx_usec_timer_elapsed(&timer) / 1000);
  printf("idct%dx%d_%d (bitdepth %d) time: %5d ms ", size_, size_,
         last_nonzero_, bit_depth_, elapsed_time);

  ASSERT_EQ(0, memcmp(output_block_ref_, output_block_,
                      sizeof(*output_block_) * output_block_size_))
      << "Error: partial inverse transform produces different results";
}

using std::tr1::make_tuple;

#if CONFIG_VP9_HIGHBITDEPTH

INSTANTIATE_TEST_CASE_P(
    C, PartialIDctTest,
    ::testing::Values(
        make_tuple(&vpx_highbd_fdct32x32_c, &vpx_highbd_idct32x32_1024_add_c,
                   &vpx_highbd_idct32x32_1024_add_c, TX_32X32, 1024, 8),
        make_tuple(&vpx_highbd_fdct32x32_c, &vpx_highbd_idct32x32_1024_add_c,
                   &vpx_highbd_idct32x32_1024_add_c, TX_32X32, 1024, 10),
        make_tuple(&vpx_highbd_fdct32x32_c, &vpx_highbd_idct32x32_1024_add_c,
                   &vpx_highbd_idct32x32_1024_add_c, TX_32X32, 1024, 12),
        make_tuple(&vpx_highbd_fdct32x32_c, &vpx_highbd_idct32x32_1024_add_c,
                   &vpx_highbd_idct32x32_34_add_c, TX_32X32, 34, 8),
        make_tuple(&vpx_highbd_fdct32x32_c, &vpx_highbd_idct32x32_1024_add_c,
                   &vpx_highbd_idct32x32_34_add_c, TX_32X32, 34, 10),
        make_tuple(&vpx_highbd_fdct32x32_c, &vpx_highbd_idct32x32_1024_add_c,
                   &vpx_highbd_idct32x32_34_add_c, TX_32X32, 34, 12),
        make_tuple(&vpx_highbd_fdct32x32_c, &vpx_highbd_idct32x32_1024_add_c,
                   &vpx_highbd_idct32x32_1_add_c, TX_32X32, 1, 8),
        make_tuple(&vpx_highbd_fdct32x32_c, &vpx_highbd_idct32x32_1024_add_c,
                   &vpx_highbd_idct32x32_1_add_c, TX_32X32, 1, 10),
        make_tuple(&vpx_highbd_fdct32x32_c, &vpx_highbd_idct32x32_1024_add_c,
                   &vpx_highbd_idct32x32_1_add_c, TX_32X32, 1, 12),
        make_tuple(&vpx_highbd_fdct16x16_c, &vpx_highbd_idct16x16_256_add_c,
                   &vpx_highbd_idct16x16_256_add_c, TX_16X16, 256, 8),
        make_tuple(&vpx_highbd_fdct16x16_c, &vpx_highbd_idct16x16_256_add_c,
                   &vpx_highbd_idct16x16_256_add_c, TX_16X16, 256, 10),
        make_tuple(&vpx_highbd_fdct16x16_c, &vpx_highbd_idct16x16_256_add_c,
                   &vpx_highbd_idct16x16_256_add_c, TX_16X16, 256, 12),
        make_tuple(&vpx_highbd_fdct16x16_c, &vpx_highbd_idct16x16_256_add_c,
                   &vpx_highbd_idct16x16_10_add_c, TX_16X16, 10, 8),
        make_tuple(&vpx_highbd_fdct16x16_c, &vpx_highbd_idct16x16_256_add_c,
                   &vpx_highbd_idct16x16_10_add_c, TX_16X16, 10, 10),
        make_tuple(&vpx_highbd_fdct16x16_c, &vpx_highbd_idct16x16_256_add_c,
                   &vpx_highbd_idct16x16_10_add_c, TX_16X16, 10, 12),
        make_tuple(&vpx_highbd_fdct16x16_c, &vpx_highbd_idct16x16_256_add_c,
                   &vpx_highbd_idct16x16_1_add_c, TX_16X16, 1, 8),
        make_tuple(&vpx_highbd_fdct16x16_c, &vpx_highbd_idct16x16_256_add_c,
                   &vpx_highbd_idct16x16_1_add_c, TX_16X16, 1, 10),
        make_tuple(&vpx_highbd_fdct16x16_c, &vpx_highbd_idct16x16_256_add_c,
                   &vpx_highbd_idct16x16_1_add_c, TX_16X16, 1, 12),
        make_tuple(&vpx_highbd_fdct8x8_c, &vpx_highbd_idct8x8_64_add_c,
                   &vpx_highbd_idct8x8_64_add_c, TX_8X8, 64, 8),
        make_tuple(&vpx_highbd_fdct8x8_c, &vpx_highbd_idct8x8_64_add_c,
                   &vpx_highbd_idct8x8_64_add_c, TX_8X8, 64, 10),
        make_tuple(&vpx_highbd_fdct8x8_c, &vpx_highbd_idct8x8_64_add_c,
                   &vpx_highbd_idct8x8_64_add_c, TX_8X8, 64, 12),
        make_tuple(&vpx_highbd_fdct8x8_c, &vpx_highbd_idct8x8_64_add_c,
                   &vpx_highbd_idct8x8_12_add_c, TX_8X8, 12, 8),
        make_tuple(&vpx_highbd_fdct8x8_c, &vpx_highbd_idct8x8_64_add_c,
                   &vpx_highbd_idct8x8_12_add_c, TX_8X8, 12, 10),
        make_tuple(&vpx_highbd_fdct8x8_c, &vpx_highbd_idct8x8_64_add_c,
                   &vpx_highbd_idct8x8_12_add_c, TX_8X8, 12, 12),
        make_tuple(&vpx_highbd_fdct8x8_c, &vpx_highbd_idct8x8_64_add_c,
                   &vpx_highbd_idct8x8_1_add_c, TX_8X8, 1, 8),
        make_tuple(&vpx_highbd_fdct8x8_c, &vpx_highbd_idct8x8_64_add_c,
                   &vpx_highbd_idct8x8_1_add_c, TX_8X8, 1, 10),
        make_tuple(&vpx_highbd_fdct8x8_c, &vpx_highbd_idct8x8_64_add_c,
                   &vpx_highbd_idct8x8_1_add_c, TX_8X8, 1, 12),
        make_tuple(&vpx_highbd_fdct4x4_c, &vpx_highbd_idct4x4_16_add_c,
                   &vpx_highbd_idct4x4_16_add_c, TX_4X4, 16, 8),
        make_tuple(&vpx_highbd_fdct4x4_c, &vpx_highbd_idct4x4_16_add_c,
                   &vpx_highbd_idct4x4_16_add_c, TX_4X4, 16, 10),
        make_tuple(&vpx_highbd_fdct4x4_c, &vpx_highbd_idct4x4_16_add_c,
                   &vpx_highbd_idct4x4_16_add_c, TX_4X4, 16, 12),
        make_tuple(&vpx_highbd_fdct4x4_c, &vpx_highbd_idct4x4_16_add_c,
                   &vpx_highbd_idct4x4_1_add_c, TX_4X4, 1, 8),
        make_tuple(&vpx_highbd_fdct4x4_c, &vpx_highbd_idct4x4_16_add_c,
                   &vpx_highbd_idct4x4_1_add_c, TX_4X4, 1, 10),
        make_tuple(&vpx_highbd_fdct4x4_c, &vpx_highbd_idct4x4_16_add_c,
                   &vpx_highbd_idct4x4_1_add_c, TX_4X4, 1, 12)));

#if HAVE_SSE2 && !CONFIG_EMULATE_HARDWARE
INSTANTIATE_TEST_CASE_P(
    SSE2, PartialIDctTest,
    ::testing::Values(
        make_tuple(&vpx_highbd_fdct32x32_c, &vpx_highbd_idct32x32_1024_add_c,
                   &vpx_highbd_idct32x32_1_add_sse2, TX_32X32, 1, 8),
        make_tuple(&vpx_highbd_fdct32x32_c, &vpx_highbd_idct32x32_1024_add_c,
                   &vpx_highbd_idct32x32_1_add_sse2, TX_32X32, 1, 10),
        make_tuple(&vpx_highbd_fdct32x32_c, &vpx_highbd_idct32x32_1024_add_c,
                   &vpx_highbd_idct32x32_1_add_sse2, TX_32X32, 1, 12),
        make_tuple(&vpx_highbd_fdct16x16_c, &vpx_highbd_idct16x16_256_add_c,
                   &vpx_highbd_idct16x16_256_add_sse2, TX_16X16, 256, 8),
        make_tuple(&vpx_highbd_fdct16x16_c, &vpx_highbd_idct16x16_256_add_c,
                   &vpx_highbd_idct16x16_256_add_sse2, TX_16X16, 256, 10),
        make_tuple(&vpx_highbd_fdct16x16_c, &vpx_highbd_idct16x16_256_add_c,
                   &vpx_highbd_idct16x16_256_add_sse2, TX_16X16, 256, 12),
        make_tuple(&vpx_highbd_fdct16x16_c, &vpx_highbd_idct16x16_256_add_c,
                   &vpx_highbd_idct16x16_10_add_sse2, TX_16X16, 10, 8),
        make_tuple(&vpx_highbd_fdct16x16_c, &vpx_highbd_idct16x16_256_add_c,
                   &vpx_highbd_idct16x16_10_add_sse2, TX_16X16, 10, 10),
        make_tuple(&vpx_highbd_fdct16x16_c, &vpx_highbd_idct16x16_256_add_c,
                   &vpx_highbd_idct16x16_10_add_sse2, TX_16X16, 10, 12),
        make_tuple(&vpx_highbd_fdct8x8_c, &vpx_highbd_idct8x8_64_add_c,
                   &vpx_highbd_idct8x8_64_add_sse2, TX_8X8, 64, 8),
        make_tuple(&vpx_highbd_fdct8x8_c, &vpx_highbd_idct8x8_64_add_c,
                   &vpx_highbd_idct8x8_64_add_sse2, TX_8X8, 64, 10),
        make_tuple(&vpx_highbd_fdct8x8_c, &vpx_highbd_idct8x8_64_add_c,
                   &vpx_highbd_idct8x8_64_add_sse2, TX_8X8, 64, 12),
        make_tuple(&vpx_highbd_fdct8x8_c, &vpx_highbd_idct8x8_64_add_c,
                   &vpx_highbd_idct8x8_12_add_sse2, TX_8X8, 12, 8),
        make_tuple(&vpx_highbd_fdct8x8_c, &vpx_highbd_idct8x8_64_add_c,
                   &vpx_highbd_idct8x8_12_add_sse2, TX_8X8, 12, 10),
        make_tuple(&vpx_highbd_fdct8x8_c, &vpx_highbd_idct8x8_64_add_c,
                   &vpx_highbd_idct8x8_12_add_sse2, TX_8X8, 12, 12),
        make_tuple(&vpx_highbd_fdct4x4_c, &vpx_highbd_idct4x4_16_add_c,
                   &vpx_highbd_idct4x4_16_add_sse2, TX_4X4, 1, 8),
        make_tuple(&vpx_highbd_fdct4x4_c, &vpx_highbd_idct4x4_16_add_c,
                   &vpx_highbd_idct4x4_16_add_sse2, TX_4X4, 1, 10),
        make_tuple(&vpx_highbd_fdct4x4_c, &vpx_highbd_idct4x4_16_add_c,
                   &vpx_highbd_idct4x4_16_add_sse2, TX_4X4, 1, 12)));
#endif  // HAVE_SSE2 && !CONFIG_EMULATE_HARDWARE

#if HAVE_NEON && !CONFIG_EMULATE_HARDWARE
INSTANTIATE_TEST_CASE_P(
    NEON, PartialIDctTest,
    ::testing::Values(
        make_tuple(&vpx_highbd_fdct4x4_c, &vpx_highbd_idct4x4_16_add_c,
                   &vpx_highbd_idct4x4_16_add_neon, TX_4X4, 16, 8),
        make_tuple(&vpx_highbd_fdct4x4_c, &vpx_highbd_idct4x4_16_add_c,
                   &vpx_highbd_idct4x4_16_add_neon, TX_4X4, 16, 10),
        make_tuple(&vpx_highbd_fdct4x4_c, &vpx_highbd_idct4x4_16_add_c,
                   &vpx_highbd_idct4x4_16_add_neon, TX_4X4, 16, 12),
        make_tuple(&vpx_highbd_fdct4x4_c, &vpx_highbd_idct4x4_1_add_c,
                   &vpx_highbd_idct4x4_1_add_neon, TX_4X4, 1, 8),
        make_tuple(&vpx_highbd_fdct4x4_c, &vpx_highbd_idct4x4_1_add_c,
                   &vpx_highbd_idct4x4_1_add_neon, TX_4X4, 1, 10),
        make_tuple(&vpx_highbd_fdct4x4_c, &vpx_highbd_idct4x4_1_add_c,
                   &vpx_highbd_idct4x4_1_add_neon, TX_4X4, 1, 12)));
#endif  // HAVE_NEON && !CONFIG_EMULATE_HARDWARE

#else  // !CONFIG_VP9_HIGHBITDEPTH

INSTANTIATE_TEST_CASE_P(
    C, PartialIDctTest,
    ::testing::Values(make_tuple(&vpx_fdct32x32_c, &vpx_idct32x32_1024_add_c,
                                 &vpx_idct32x32_1024_add_c, TX_32X32, 1024, 8),
                      make_tuple(&vpx_fdct32x32_c, &vpx_idct32x32_1024_add_c,
                                 &vpx_idct32x32_135_add_c, TX_32X32, 135, 8),
                      make_tuple(&vpx_fdct32x32_c, &vpx_idct32x32_1024_add_c,
                                 &vpx_idct32x32_34_add_c, TX_32X32, 34, 8),
                      make_tuple(&vpx_fdct32x32_c, &vpx_idct32x32_1024_add_c,
                                 &vpx_idct32x32_1_add_c, TX_32X32, 1, 8),
                      make_tuple(&vpx_fdct16x16_c, &vpx_idct16x16_256_add_c,
                                 &vpx_idct16x16_256_add_c, TX_16X16, 256, 8),
                      make_tuple(&vpx_fdct16x16_c, &vpx_idct16x16_256_add_c,
                                 &vpx_idct16x16_10_add_c, TX_16X16, 10, 8),
                      make_tuple(&vpx_fdct16x16_c, &vpx_idct16x16_256_add_c,
                                 &vpx_idct16x16_1_add_c, TX_16X16, 1, 8),
                      make_tuple(&vpx_fdct8x8_c, &vpx_idct8x8_64_add_c,
                                 &vpx_idct8x8_64_add_c, TX_8X8, 64, 8),
                      make_tuple(&vpx_fdct8x8_c, &vpx_idct8x8_64_add_c,
                                 &vpx_idct8x8_12_add_c, TX_8X8, 12, 8),
                      make_tuple(&vpx_fdct8x8_c, &vpx_idct8x8_64_add_c,
                                 &vpx_idct8x8_1_add_c, TX_8X8, 1, 8),
                      make_tuple(&vpx_fdct4x4_c, &vpx_idct4x4_16_add_c,
                                 &vpx_idct4x4_16_add_c, TX_4X4, 16, 8),
                      make_tuple(&vpx_fdct4x4_c, &vpx_idct4x4_16_add_c,
                                 &vpx_idct4x4_1_add_c, TX_4X4, 1, 8)));

#if HAVE_NEON && !CONFIG_EMULATE_HARDWARE
INSTANTIATE_TEST_CASE_P(
    NEON, PartialIDctTest,
    ::testing::Values(make_tuple(&vpx_fdct32x32_c, &vpx_idct32x32_1024_add_c,
                                 &vpx_idct32x32_1024_add_neon, TX_32X32, 1024,
                                 8),
                      make_tuple(&vpx_fdct32x32_c, &vpx_idct32x32_1024_add_c,
                                 &vpx_idct32x32_135_add_neon, TX_32X32, 135, 8),
                      make_tuple(&vpx_fdct32x32_c, &vpx_idct32x32_1024_add_c,
                                 &vpx_idct32x32_34_add_neon, TX_32X32, 34, 8),
                      make_tuple(&vpx_fdct32x32_c, &vpx_idct32x32_1024_add_c,
                                 &vpx_idct32x32_1_add_neon, TX_32X32, 1, 8),
                      make_tuple(&vpx_fdct16x16_c, &vpx_idct16x16_256_add_c,
                                 &vpx_idct16x16_256_add_neon, TX_16X16, 256, 8),
                      make_tuple(&vpx_fdct16x16_c, &vpx_idct16x16_256_add_c,
                                 &vpx_idct16x16_10_add_neon, TX_16X16, 10, 8),
                      make_tuple(&vpx_fdct16x16_c, &vpx_idct16x16_256_add_c,
                                 &vpx_idct16x16_1_add_neon, TX_16X16, 1, 8),
                      make_tuple(&vpx_fdct8x8_c, &vpx_idct8x8_64_add_c,
                                 &vpx_idct8x8_64_add_neon, TX_8X8, 64, 8),
                      make_tuple(&vpx_fdct8x8_c, &vpx_idct8x8_64_add_c,
                                 &vpx_idct8x8_12_add_neon, TX_8X8, 12, 8),
                      make_tuple(&vpx_fdct8x8_c, &vpx_idct8x8_64_add_c,
                                 &vpx_idct8x8_1_add_neon, TX_8X8, 1, 8),
                      make_tuple(&vpx_fdct4x4_c, &vpx_idct4x4_16_add_c,
                                 &vpx_idct4x4_16_add_neon, TX_4X4, 16, 8),
                      make_tuple(&vpx_fdct4x4_c, &vpx_idct4x4_16_add_c,
                                 &vpx_idct4x4_1_add_neon, TX_4X4, 1, 8)));
#endif  // HAVE_NEON && !CONFIG_EMULATE_HARDWARE

#if HAVE_SSE2 && !CONFIG_EMULATE_HARDWARE
// 32x32_135_ is implemented using the 1024 version.
INSTANTIATE_TEST_CASE_P(
    SSE2, PartialIDctTest,
    ::testing::Values(make_tuple(&vpx_fdct32x32_c, &vpx_idct32x32_1024_add_c,
                                 &vpx_idct32x32_1024_add_sse2, TX_32X32, 1024,
                                 8),
                      make_tuple(&vpx_fdct32x32_c, &vpx_idct32x32_1024_add_c,
                                 &vpx_idct32x32_1024_add_sse2, TX_32X32, 135,
                                 8),
                      make_tuple(&vpx_fdct32x32_c, &vpx_idct32x32_1024_add_c,
                                 &vpx_idct32x32_34_add_sse2, TX_32X32, 34, 8),
                      make_tuple(&vpx_fdct32x32_c, &vpx_idct32x32_1024_add_c,
                                 &vpx_idct32x32_1_add_sse2, TX_32X32, 1, 8),
                      make_tuple(&vpx_fdct16x16_c, &vpx_idct16x16_256_add_c,
                                 &vpx_idct16x16_256_add_sse2, TX_16X16, 256, 8),
                      make_tuple(&vpx_fdct16x16_c, &vpx_idct16x16_256_add_c,
                                 &vpx_idct16x16_10_add_sse2, TX_16X16, 10, 8),
                      make_tuple(&vpx_fdct16x16_c, &vpx_idct16x16_256_add_c,
                                 &vpx_idct16x16_1_add_sse2, TX_16X16, 1, 8),
                      make_tuple(&vpx_fdct8x8_c, &vpx_idct8x8_64_add_c,
                                 &vpx_idct8x8_64_add_sse2, TX_8X8, 64, 8),
                      make_tuple(&vpx_fdct8x8_c, &vpx_idct8x8_64_add_c,
                                 &vpx_idct8x8_12_add_sse2, TX_8X8, 12, 8),
                      make_tuple(&vpx_fdct8x8_c, &vpx_idct8x8_64_add_c,
                                 &vpx_idct8x8_1_add_sse2, TX_8X8, 1, 8),
                      make_tuple(&vpx_fdct4x4_c, &vpx_idct4x4_16_add_c,
                                 &vpx_idct4x4_16_add_sse2, TX_4X4, 16, 8),
                      make_tuple(&vpx_fdct4x4_c, &vpx_idct4x4_16_add_c,
                                 &vpx_idct4x4_1_add_sse2, TX_4X4, 1, 8)));
#endif  // HAVE_SSE2 && !CONFIG_EMULATE_HARDWARE

#if HAVE_SSSE3 && ARCH_X86_64 && !CONFIG_EMULATE_HARDWARE
INSTANTIATE_TEST_CASE_P(
    SSSE3_64, PartialIDctTest,
    ::testing::Values(make_tuple(&vpx_fdct32x32_c, &vpx_idct32x32_1024_add_c,
                                 &vpx_idct32x32_1024_add_ssse3, TX_32X32, 1024,
                                 8),
                      make_tuple(&vpx_fdct32x32_c, &vpx_idct32x32_1024_add_c,
                                 &vpx_idct32x32_135_add_ssse3, TX_32X32, 135,
                                 8),
                      make_tuple(&vpx_fdct32x32_c, &vpx_idct32x32_1024_add_c,
                                 &vpx_idct32x32_34_add_ssse3, TX_32X32, 34, 8),
                      make_tuple(&vpx_fdct8x8_c, &vpx_idct8x8_64_add_c,
                                 &vpx_idct8x8_64_add_ssse3, TX_8X8, 64, 8),
                      make_tuple(&vpx_fdct8x8_c, &vpx_idct8x8_64_add_c,
                                 &vpx_idct8x8_12_add_ssse3, TX_8X8, 12, 8)));
#endif  // HAVE_SSSE3 && ARCH_X86_64 && !CONFIG_EMULATE_HARDWARE

#if HAVE_MSA && !CONFIG_EMULATE_HARDWARE
// 32x32_135_ is implemented using the 1024 version.
INSTANTIATE_TEST_CASE_P(
    MSA, PartialIDctTest,
    ::testing::Values(make_tuple(&vpx_fdct32x32_c, &vpx_idct32x32_1024_add_c,
                                 &vpx_idct32x32_1024_add_msa, TX_32X32, 1024,
                                 8),
                      make_tuple(&vpx_fdct32x32_c, &vpx_idct32x32_1024_add_c,
                                 &vpx_idct32x32_1024_add_msa, TX_32X32, 135, 8),
                      make_tuple(&vpx_fdct32x32_c, &vpx_idct32x32_1024_add_c,
                                 &vpx_idct32x32_34_add_msa, TX_32X32, 34, 8),
                      make_tuple(&vpx_fdct32x32_c, &vpx_idct32x32_1024_add_c,
                                 &vpx_idct32x32_1_add_msa, TX_32X32, 1, 8),
                      make_tuple(&vpx_fdct16x16_c, &vpx_idct16x16_256_add_c,
                                 &vpx_idct16x16_256_add_msa, TX_16X16, 256, 8),
                      make_tuple(&vpx_fdct16x16_c, &vpx_idct16x16_256_add_c,
                                 &vpx_idct16x16_10_add_msa, TX_16X16, 10, 8),
                      make_tuple(&vpx_fdct16x16_c, &vpx_idct16x16_256_add_c,
                                 &vpx_idct16x16_1_add_msa, TX_16X16, 1, 8),
                      make_tuple(&vpx_fdct8x8_c, &vpx_idct8x8_64_add_c,
                                 &vpx_idct8x8_64_add_msa, TX_8X8, 64, 8),
                      make_tuple(&vpx_fdct8x8_c, &vpx_idct8x8_64_add_c,
                                 &vpx_idct8x8_12_add_msa, TX_8X8, 10, 8),
                      make_tuple(&vpx_fdct8x8_c, &vpx_idct8x8_64_add_c,
                                 &vpx_idct8x8_1_add_msa, TX_8X8, 1, 8),
                      make_tuple(&vpx_fdct4x4_c, &vpx_idct4x4_16_add_c,
                                 &vpx_idct4x4_16_add_msa, TX_4X4, 16, 8),
                      make_tuple(&vpx_fdct4x4_c, &vpx_idct4x4_16_add_c,
                                 &vpx_idct4x4_1_add_msa, TX_4X4, 1, 8)));
#endif  // HAVE_MSA && !CONFIG_EMULATE_HARDWARE

#endif  // CONFIG_VP9_HIGHBITDEPTH

}  // namespace
