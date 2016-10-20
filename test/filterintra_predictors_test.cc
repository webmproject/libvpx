/*
 *  Copyright (c) 2016 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include "third_party/googletest/src/include/gtest/gtest.h"

#include "./av1_rtcd.h"
#include "test/acm_random.h"
#include "test/clear_system_state.h"
#include "test/register_state_check.h"
#include "test/util.h"
#include "av1/common/enums.h"

namespace {

using std::tr1::tuple;
using libaom_test::ACMRandom;

typedef void (*Predictor)(uint8_t *dst, ptrdiff_t stride, int bs,
                          const uint8_t *above, const uint8_t *left);

// Note:
//  Test parameter list:
//  Reference predictor, optimized predictor, prediction mode, block size
//
typedef tuple<Predictor, Predictor, int> PredFuncMode;
typedef tuple<PredFuncMode, int> PredParams;

#if CONFIG_AOM_HIGHBITDEPTH
typedef void (*HbdPredictor)(uint16_t *dst, ptrdiff_t stride, int bs,
                             const uint16_t *above, const uint16_t *left,
                             int bd);

// Note:
//  Test parameter list:
//  Reference predictor, optimized predictor, prediction mode, block size,
//  bit depth
//
typedef tuple<HbdPredictor, HbdPredictor, int> HbdPredFuncMode;
typedef tuple<HbdPredFuncMode, int, int> HbdPredParams;
#endif

const int MaxBlkSize = 32;

// By default, disable speed test
#define PREDICTORS_SPEED_TEST (0)

#if PREDICTORS_SPEED_TEST
const int MaxTestNum = 100000;
#else
const int MaxTestNum = 100;
#endif

class AV1FilterIntraPredOptimzTest
    : public ::testing::TestWithParam<PredParams> {
 public:
  virtual ~AV1FilterIntraPredOptimzTest() {}
  virtual void SetUp() {
    PredFuncMode funcMode = GET_PARAM(0);
    predFuncRef_ = std::tr1::get<0>(funcMode);
    predFunc_ = std::tr1::get<1>(funcMode);
    mode_ = std::tr1::get<2>(funcMode);
    blockSize_ = GET_PARAM(1);

    alloc_ = new uint8_t[3 * MaxBlkSize + 2];
    predRef_ = new uint8_t[MaxBlkSize * MaxBlkSize];
    pred_ = new uint8_t[MaxBlkSize * MaxBlkSize];
  }

  virtual void TearDown() {
    delete[] alloc_;
    delete[] predRef_;
    delete[] pred_;
    libaom_test::ClearSystemState();
  }

 protected:
  void RunTest() const {
    int tstIndex = 0;
    int stride = blockSize_;
    uint8_t *left = alloc_;
    uint8_t *above = alloc_ + MaxBlkSize + 1;
    while (tstIndex < MaxTestNum) {
      PrepareBuffer();
      predFuncRef_(predRef_, stride, blockSize_, &above[1], left);
      ASM_REGISTER_STATE_CHECK(
          predFunc_(pred_, stride, blockSize_, &above[1], left));
      DiffPred(tstIndex);
      tstIndex += 1;
    }
  }

  void RunSpeedTestC() const {
    int tstIndex = 0;
    int stride = blockSize_;
    uint8_t *left = alloc_;
    uint8_t *above = alloc_ + MaxBlkSize + 1;
    PrepareBuffer();
    while (tstIndex < MaxTestNum) {
      predFuncRef_(predRef_, stride, blockSize_, &above[1], left);
      tstIndex += 1;
    }
  }

  void RunSpeedTestSSE() const {
    int tstIndex = 0;
    int stride = blockSize_;
    uint8_t *left = alloc_;
    uint8_t *above = alloc_ + MaxBlkSize + 1;
    PrepareBuffer();
    while (tstIndex < MaxTestNum) {
      predFunc_(predRef_, stride, blockSize_, &above[1], left);
      tstIndex += 1;
    }
  }

 private:
  void PrepareBuffer() const {
    ACMRandom rnd(ACMRandom::DeterministicSeed());
    int i = 0;
    while (i < (3 * MaxBlkSize + 2)) {
      alloc_[i] = rnd.Rand8();
      i += 1;
    }
  }

  void DiffPred(int testNum) const {
    int i = 0;
    while (i < blockSize_ * blockSize_) {
      EXPECT_EQ(predRef_[i], pred_[i]) << "Error at position: " << i << " "
                                       << "Block size: " << blockSize_ << " "
                                       << "Test number: " << testNum;
      i += 1;
    }
  }

  Predictor predFunc_;
  Predictor predFuncRef_;
  int mode_;
  int blockSize_;
  uint8_t *alloc_;
  uint8_t *pred_;
  uint8_t *predRef_;
};

#if CONFIG_AOM_HIGHBITDEPTH
class AV1HbdFilterIntraPredOptimzTest
    : public ::testing::TestWithParam<HbdPredParams> {
 public:
  virtual ~AV1HbdFilterIntraPredOptimzTest() {}
  virtual void SetUp() {
    HbdPredFuncMode funcMode = GET_PARAM(0);
    predFuncRef_ = std::tr1::get<0>(funcMode);
    predFunc_ = std::tr1::get<1>(funcMode);
    mode_ = std::tr1::get<2>(funcMode);
    blockSize_ = GET_PARAM(1);
    bd_ = GET_PARAM(2);

    alloc_ = new uint16_t[3 * MaxBlkSize + 2];
    predRef_ = new uint16_t[MaxBlkSize * MaxBlkSize];
    pred_ = new uint16_t[MaxBlkSize * MaxBlkSize];
  }

  virtual void TearDown() {
    delete[] alloc_;
    delete[] predRef_;
    delete[] pred_;
    libaom_test::ClearSystemState();
  }

 protected:
  void RunTest() const {
    int tstIndex = 0;
    int stride = blockSize_;
    uint16_t *left = alloc_;
    uint16_t *above = alloc_ + MaxBlkSize + 1;
    while (tstIndex < MaxTestNum) {
      PrepareBuffer();
      predFuncRef_(predRef_, stride, blockSize_, &above[1], left, bd_);
      ASM_REGISTER_STATE_CHECK(
          predFunc_(pred_, stride, blockSize_, &above[1], left, bd_));
      DiffPred(tstIndex);
      tstIndex += 1;
    }
  }

  void RunSpeedTestC() const {
    int tstIndex = 0;
    int stride = blockSize_;
    uint16_t *left = alloc_;
    uint16_t *above = alloc_ + MaxBlkSize + 1;
    PrepareBuffer();
    while (tstIndex < MaxTestNum) {
      predFuncRef_(predRef_, stride, blockSize_, &above[1], left, bd_);
      tstIndex += 1;
    }
  }

  void RunSpeedTestSSE() const {
    int tstIndex = 0;
    int stride = blockSize_;
    uint16_t *left = alloc_;
    uint16_t *above = alloc_ + MaxBlkSize + 1;
    PrepareBuffer();
    while (tstIndex < MaxTestNum) {
      predFunc_(predRef_, stride, blockSize_, &above[1], left, bd_);
      tstIndex += 1;
    }
  }

 private:
  void PrepareBuffer() const {
    ACMRandom rnd(ACMRandom::DeterministicSeed());
    int i = 0;
    while (i < (3 * MaxBlkSize + 2)) {
      alloc_[i] = rnd.Rand16() & ((1 << bd_) - 1);
      i += 1;
    }
  }

  void DiffPred(int testNum) const {
    int i = 0;
    while (i < blockSize_ * blockSize_) {
      EXPECT_EQ(predRef_[i], pred_[i]) << "Error at position: " << i << " "
                                       << "Block size: " << blockSize_ << " "
                                       << "Bit depth: " << bd_ << " "
                                       << "Test number: " << testNum;
      i += 1;
    }
  }

  HbdPredictor predFunc_;
  HbdPredictor predFuncRef_;
  int mode_;
  int blockSize_;
  int bd_;
  uint16_t *alloc_;
  uint16_t *pred_;
  uint16_t *predRef_;
};
#endif  // CONFIG_AOM_HIGHBITDEPTH

TEST_P(AV1FilterIntraPredOptimzTest, BitExactCheck) { RunTest(); }

#if PREDICTORS_SPEED_TEST
TEST_P(AV1FilterIntraPredOptimzTest, SpeedCheckC) { RunSpeedTestC(); }

TEST_P(AV1FilterIntraPredOptimzTest, SpeedCheckSSE) { RunSpeedTestSSE(); }
#endif

#if CONFIG_AOM_HIGHBITDEPTH
TEST_P(AV1HbdFilterIntraPredOptimzTest, BitExactCheck) { RunTest(); }

#if PREDICTORS_SPEED_TEST
TEST_P(AV1HbdFilterIntraPredOptimzTest, SpeedCheckC) { RunSpeedTestC(); }

TEST_P(AV1HbdFilterIntraPredOptimzTest, SpeedCheckSSE) { RunSpeedTestSSE(); }
#endif  // PREDICTORS_SPEED_TEST
#endif  // CONFIG_AOM_HIGHBITDEPTH

using std::tr1::make_tuple;

const PredFuncMode kPredFuncMdArray[] = {
  make_tuple(av1_dc_filter_predictor_c, av1_dc_filter_predictor_sse4_1,
             DC_PRED),
  make_tuple(av1_v_filter_predictor_c, av1_v_filter_predictor_sse4_1, V_PRED),
  make_tuple(av1_h_filter_predictor_c, av1_h_filter_predictor_sse4_1, H_PRED),
  make_tuple(av1_d45_filter_predictor_c, av1_d45_filter_predictor_sse4_1,
             D45_PRED),
  make_tuple(av1_d135_filter_predictor_c, av1_d135_filter_predictor_sse4_1,
             D135_PRED),
  make_tuple(av1_d117_filter_predictor_c, av1_d117_filter_predictor_sse4_1,
             D117_PRED),
  make_tuple(av1_d153_filter_predictor_c, av1_d153_filter_predictor_sse4_1,
             D153_PRED),
  make_tuple(av1_d207_filter_predictor_c, av1_d207_filter_predictor_sse4_1,
             D207_PRED),
  make_tuple(av1_d63_filter_predictor_c, av1_d63_filter_predictor_sse4_1,
             D63_PRED),
  make_tuple(av1_tm_filter_predictor_c, av1_tm_filter_predictor_sse4_1,
             TM_PRED),
};

const int kBlkSize[] = { 4, 8, 16, 32 };

INSTANTIATE_TEST_CASE_P(
    SSE4_1, AV1FilterIntraPredOptimzTest,
    ::testing::Combine(::testing::ValuesIn(kPredFuncMdArray),
                       ::testing::ValuesIn(kBlkSize)));

#if CONFIG_AOM_HIGHBITDEPTH
const HbdPredFuncMode kHbdPredFuncMdArray[] = {
  make_tuple(av1_highbd_dc_filter_predictor_c,
             av1_highbd_dc_filter_predictor_sse4_1, DC_PRED),
  make_tuple(av1_highbd_v_filter_predictor_c,
             av1_highbd_v_filter_predictor_sse4_1, V_PRED),
  make_tuple(av1_highbd_h_filter_predictor_c,
             av1_highbd_h_filter_predictor_sse4_1, H_PRED),
  make_tuple(av1_highbd_d45_filter_predictor_c,
             av1_highbd_d45_filter_predictor_sse4_1, D45_PRED),
  make_tuple(av1_highbd_d135_filter_predictor_c,
             av1_highbd_d135_filter_predictor_sse4_1, D135_PRED),
  make_tuple(av1_highbd_d117_filter_predictor_c,
             av1_highbd_d117_filter_predictor_sse4_1, D117_PRED),
  make_tuple(av1_highbd_d153_filter_predictor_c,
             av1_highbd_d153_filter_predictor_sse4_1, D153_PRED),
  make_tuple(av1_highbd_d207_filter_predictor_c,
             av1_highbd_d207_filter_predictor_sse4_1, D207_PRED),
  make_tuple(av1_highbd_d63_filter_predictor_c,
             av1_highbd_d63_filter_predictor_sse4_1, D63_PRED),
  make_tuple(av1_highbd_tm_filter_predictor_c,
             av1_highbd_tm_filter_predictor_sse4_1, TM_PRED),
};

const int kBd[] = { 10, 12 };

INSTANTIATE_TEST_CASE_P(
    SSE4_1, AV1HbdFilterIntraPredOptimzTest,
    ::testing::Combine(::testing::ValuesIn(kHbdPredFuncMdArray),
                       ::testing::ValuesIn(kBlkSize),
                       ::testing::ValuesIn(kBd)));
#endif  // CONFIG_AOM_HIGHBITDEPTH

}  // namespace
