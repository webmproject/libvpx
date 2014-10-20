/*
 *  Copyright (c) 2014 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include <string.h>

#include "third_party/googletest/src/include/gtest/gtest.h"
#include "test/acm_random.h"
#include "test/clear_system_state.h"
#include "test/register_state_check.h"
#include "test/util.h"

#include "./vpx_config.h"
#include "./vp8_rtcd.h"
#include "vp8/common/blockd.h"
#include "vp8/common/onyx.h"
#include "vp8/encoder/block.h"
#include "vp8/encoder/onyx_int.h"
#include "vp8/encoder/quantize.h"
#include "vpx/vpx_integer.h"
#include "vpx_mem/vpx_mem.h"

namespace {

const int kNumBlocks = 25;
const int kNumBlockEntries = 16;

typedef void (*VP8Quantize)(BLOCK *b, BLOCKD *d);
typedef void (*VP8QuantizePair)(BLOCK *b0, BLOCK *b1, BLOCKD *d0, BLOCKD *d1);

typedef std::tr1::tuple<VP8Quantize, VP8Quantize> VP8QuantizeParam;
typedef std::tr1::tuple<VP8QuantizePair, VP8QuantizePair> VP8QuantizePairParam;

using libvpx_test::ACMRandom;
using std::tr1::make_tuple;

// Create and populate a VP8_COMP instance which has a complete set of
// quantization inputs as well as a second MACROBLOCKD for output.
class QuantizeTestBase {
 public:
  virtual ~QuantizeTestBase() {
    vp8_remove_compressor(&vp8_comp_);
    vp8_comp_ = NULL;
    vpx_free(macroblockd_dst_);
    macroblockd_dst_ = NULL;
    libvpx_test::ClearSystemState();
  }

 protected:
  void SetupCompressor() {
    rnd_.Reset(ACMRandom::DeterministicSeed());

    // The full configuration is necessary to generate the quantization tables.
    VP8_CONFIG *const vp8_config =
        reinterpret_cast<VP8_CONFIG *>(vpx_calloc(sizeof(*vp8_config), 1));

    vp8_comp_ = vp8_create_compressor(vp8_config);

    // Set the tables based on a quantizer of 0.
    vp8_set_quantizer(vp8_comp_, 0);

    // Set up all the block/blockd pointers for the mb in vp8_comp_.
    vp8cx_frame_init_quantizer(vp8_comp_);

    // Copy macroblockd from the reference to get pre-set-up dequant values.
    macroblockd_dst_ =
        reinterpret_cast<MACROBLOCKD *>(
            vpx_calloc(sizeof(*macroblockd_dst_), 1));
    vpx_memcpy(macroblockd_dst_, &vp8_comp_->mb.e_mbd,
               sizeof(*macroblockd_dst_));
    // Fix block pointers - currently they point to the blocks in the reference
    // structure.
    vp8_setup_block_dptrs(macroblockd_dst_);
  }

  void UpdateQuantizer(int q) {
    vp8_set_quantizer(vp8_comp_, q);

    vpx_memcpy(macroblockd_dst_, &vp8_comp_->mb.e_mbd,
               sizeof(*macroblockd_dst_));
    vp8_setup_block_dptrs(macroblockd_dst_);
  }

  void FillCoeffConstant(int16_t c) {
    for (int i = 0; i < kNumBlocks * kNumBlockEntries; ++i) {
      vp8_comp_->mb.coeff[i] = c;
    }
  }

  void FillCoeffRandom() {
    for (int i = 0; i < kNumBlocks * kNumBlockEntries; ++i) {
      vp8_comp_->mb.coeff[i] = rnd_.Rand8();
    }
  }

  void CheckOutput() {
    EXPECT_EQ(0, memcmp(vp8_comp_->mb.e_mbd.qcoeff, macroblockd_dst_->qcoeff,
                        sizeof(*macroblockd_dst_->qcoeff) * kNumBlocks *
                            kNumBlockEntries))
        << "qcoeff mismatch";
    EXPECT_EQ(0, memcmp(vp8_comp_->mb.e_mbd.dqcoeff, macroblockd_dst_->dqcoeff,
                        sizeof(*macroblockd_dst_->dqcoeff) * kNumBlocks *
                            kNumBlockEntries))
        << "dqcoeff mismatch";
    EXPECT_EQ(0, memcmp(vp8_comp_->mb.e_mbd.eobs, macroblockd_dst_->eobs,
                        sizeof(*macroblockd_dst_->eobs) * kNumBlocks))
        << "eobs mismatch";
  }

  VP8_COMP *vp8_comp_;
  MACROBLOCKD *macroblockd_dst_;

 private:
  ACMRandom rnd_;
};

class QuantizeTest : public QuantizeTestBase,
                     public ::testing::TestWithParam<VP8QuantizeParam> {
 protected:
  virtual void SetUp() {
    SetupCompressor();
    asm_quant_ = GET_PARAM(0);
    c_quant_ = GET_PARAM(1);
  }

  void RunComparison() {
    for (int i = 0; i < kNumBlocks; ++i) {
      ASM_REGISTER_STATE_CHECK(
          c_quant_(&vp8_comp_->mb.block[i], &vp8_comp_->mb.e_mbd.block[i]));
      ASM_REGISTER_STATE_CHECK(
          asm_quant_(&vp8_comp_->mb.block[i], &macroblockd_dst_->block[i]));
    }

    CheckOutput();
  }

 private:
  VP8Quantize asm_quant_;
  VP8Quantize c_quant_;
};

class QuantizeTestPair : public QuantizeTestBase,
                         public ::testing::TestWithParam<VP8QuantizePairParam> {
 protected:
  virtual void SetUp() {
    SetupCompressor();
    asm_quant_pair_ = GET_PARAM(0);
    c_quant_pair_ = GET_PARAM(1);
  }

  void RunComparison() {
    // Skip the last, unpaired, block.
    for (int i = 0; i < kNumBlocks - 1; i += 2) {
      ASM_REGISTER_STATE_CHECK(c_quant_pair_(
          &vp8_comp_->mb.block[i], &vp8_comp_->mb.block[i + 1],
          &vp8_comp_->mb.e_mbd.block[i], &vp8_comp_->mb.e_mbd.block[i + 1]));
      ASM_REGISTER_STATE_CHECK(asm_quant_pair_(
          &vp8_comp_->mb.block[i], &vp8_comp_->mb.block[i + 1],
          &macroblockd_dst_->block[i], &macroblockd_dst_->block[i + 1]));
    }

    CheckOutput();
  }

 private:
  VP8QuantizePair asm_quant_pair_;
  VP8QuantizePair c_quant_pair_;
};

TEST_P(QuantizeTest, TestZeroInput) {
  FillCoeffConstant(0);
  RunComparison();
}

TEST_P(QuantizeTest, TestRandomInput) {
  FillCoeffRandom();
  RunComparison();
}

TEST_P(QuantizeTest, TestMultipleQ) {
  for (int q = 0; q < QINDEX_RANGE; ++q) {
    UpdateQuantizer(q);
    FillCoeffRandom();
    RunComparison();
  }
}

TEST_P(QuantizeTestPair, TestZeroInput) {
  FillCoeffConstant(0);
  RunComparison();
}

TEST_P(QuantizeTestPair, TestRandomInput) {
  FillCoeffRandom();
  RunComparison();
}

TEST_P(QuantizeTestPair, TestMultipleQ) {
  for (int q = 0; q < QINDEX_RANGE; ++q) {
    UpdateQuantizer(q);
    FillCoeffRandom();
    RunComparison();
  }
}

#if HAVE_SSE2
INSTANTIATE_TEST_CASE_P(
    SSE2, QuantizeTest,
    ::testing::Values(
        make_tuple(vp8_fast_quantize_b_sse2, vp8_fast_quantize_b_c),
        make_tuple(vp8_regular_quantize_b_sse2, vp8_regular_quantize_b_c)));
#endif  // HAVE_SSE2

#if HAVE_SSSE3
INSTANTIATE_TEST_CASE_P(SSSE3, QuantizeTest,
                        ::testing::Values(make_tuple(vp8_fast_quantize_b_ssse3,
                                                     vp8_fast_quantize_b_c)));
#endif  // HAVE_SSSE3

#if HAVE_SSE4_1
INSTANTIATE_TEST_CASE_P(SSE4_1, QuantizeTest, ::testing::Values(make_tuple(
                                                  vp8_regular_quantize_b_sse4_1,
                                                  vp8_regular_quantize_b_c)));
#endif  // HAVE_SSE4_1

#if HAVE_MEDIA
INSTANTIATE_TEST_CASE_P(MEDIA, QuantizeTest,
                        ::testing::Values(make_tuple(vp8_fast_quantize_b_armv6,
                                                     vp8_fast_quantize_b_c)));
#endif  // HAVE_MEDIA

#if HAVE_NEON_ASM
INSTANTIATE_TEST_CASE_P(NEON, QuantizeTest,
                        ::testing::Values(make_tuple(vp8_fast_quantize_b_neon,
                                                     vp8_fast_quantize_b_c)));

INSTANTIATE_TEST_CASE_P(
    NEON, QuantizeTestPair,
    ::testing::Values(make_tuple(vp8_fast_quantize_b_pair_neon,
                                 vp8_fast_quantize_b_pair_c)));
#endif  // HAVE_NEON_ASM
}  // namespace
