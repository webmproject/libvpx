/*
 *  Copyright (c) 2014 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include <cmath>
#include <cstdlib>
#include <string>

#include "third_party/googletest/src/include/gtest/gtest.h"

#include "./vpx_config.h"
#include "./vpx_dsp_rtcd.h"
#include "vpx_ports/mem.h"
#include "test/acm_random.h"
#include "test/clear_system_state.h"
#include "test/register_state_check.h"
#include "test/util.h"
#include "test/function_equivalence_test.h"

using libvpx_test::ACMRandom;
using libvpx_test::FunctionEquivalenceTest;

namespace {
const int kNumIterations = 10000;

static const int16_t kInt13Max = (1 << 12) - 1;

typedef uint64_t (*SSI16Func)(const int16_t *src,
                             int stride ,
                             int size);

typedef std::tr1::tuple<SSI16Func, SSI16Func> SumSquaresParam;

class SumSquaresTest : public ::testing::TestWithParam<SumSquaresParam> {
 public:
  virtual ~SumSquaresTest() {}
  virtual void SetUp() {
    ref_func_ = GET_PARAM(0);
    tst_func_ = GET_PARAM(1);
  }

  virtual void TearDown() { libvpx_test::ClearSystemState(); }

 protected:
  SSI16Func ref_func_;
  SSI16Func tst_func_;
};

TEST_P(SumSquaresTest, OperationCheck) {
  ACMRandom rnd(ACMRandom::DeterministicSeed());
  DECLARE_ALIGNED(16, int16_t, src[256*256]);

  int failed = 0;

  const int msb = 11;   // Up to 12 bit input
  const int limit = 1 << (msb+1);

  for (int k = 0; k < kNumIterations; k++) {
    int size = 4 << rnd(6);     // Up to 128x128
    int stride = 4 << rnd(7);   // Up to 256 stride
    while (stride < size) {     // Make sure it's valid
      stride = 4 << rnd(7);
    }

    for (int ii = 0 ; ii < size; ii++) {
      for (int jj = 0; jj < size; jj++) {
        src[ii*stride+jj] = rnd(2) ? rnd(limit) : -rnd(limit);
      }
    }

    uint64_t res_ref = ref_func_(src, stride, size);
    uint64_t res_tst;
    ASM_REGISTER_STATE_CHECK(res_tst = tst_func_(src, stride, size));

    if (!failed) {
      failed = res_ref != res_tst;
      EXPECT_EQ(res_ref, res_tst)
        << "Error: Sum Squares Test"
        << " C output does not match optimized output.";
    }
  }
}

TEST_P(SumSquaresTest, ExtremeValues) {
  ACMRandom rnd(ACMRandom::DeterministicSeed());
  DECLARE_ALIGNED(16, int16_t, src[256*256]);

  int failed = 0;

  const int msb = 11;   // Up to 12 bit input
  const int limit = 1 << (msb+1);

  for (int k = 0; k < kNumIterations; k++) {
    int size = 4 << rnd(6);     // Up to 128x128
    int stride = 4 << rnd(7);   // Up to 256 stride
    while (stride < size) {     // Make sure it's valid
      stride = 4 << rnd(7);
    }

    int val = rnd(2) ? limit-1 : -(limit-1);
    for (int ii = 0 ; ii < size; ii++) {
      for (int jj = 0; jj < size; jj++) {
        src[ii*stride+jj] = val;
      }
    }

    uint64_t res_ref = ref_func_(src, stride, size);
    uint64_t res_tst;
    ASM_REGISTER_STATE_CHECK(res_tst = tst_func_(src, stride, size));

    if (!failed) {
      failed = res_ref != res_tst;
      EXPECT_EQ(res_ref, res_tst)
        << "Error: Sum Squares Test"
        << " C output does not match optimized output.";
    }
  }
}
using std::tr1::make_tuple;

#if HAVE_SSE2

INSTANTIATE_TEST_CASE_P(
    SSE2, SumSquaresTest,
    ::testing::Values(
        make_tuple(&vpx_sum_squares_2d_i16_c, &vpx_sum_squares_2d_i16_sse2)
    )
);
#endif  // HAVE_SSE2

//////////////////////////////////////////////////////////////////////////////
// 1D version
//////////////////////////////////////////////////////////////////////////////

typedef uint64_t (*F1D)(const int16_t *src, uint32_t N);

class SumSquares1DTest : public FunctionEquivalenceTest<F1D> {
 protected:
  SumSquares1DTest() : rng_(ACMRandom::DeterministicSeed()) {}

  static const int kIterations = 1000;
  static const int kMaxSize = 256;

  ACMRandom rng_;
};

TEST_P(SumSquares1DTest, RandomValues) {
  DECLARE_ALIGNED(16, int16_t, src[kMaxSize * kMaxSize]);

  for (int iter = 0 ; iter < kIterations && !HasFatalFailure(); ++iter) {
    for (int i = 0 ; i < kMaxSize * kMaxSize ; ++i)
      src[i] = rng_(kInt13Max * 2 + 1) - kInt13Max;

    const int N = rng_(2) ? rng_(kMaxSize * kMaxSize + 1 - kMaxSize) + kMaxSize
                          : rng_(kMaxSize) + 1;

    const uint64_t ref_res = ref_func_(src, N);
    const uint64_t tst_res = tst_func_(src, N);

    ASSERT_EQ(ref_res, tst_res);
  }
}

TEST_P(SumSquares1DTest, ExtremeValues) {
  DECLARE_ALIGNED(16, int16_t, src[kMaxSize * kMaxSize]);

  for (int iter = 0 ; iter < kIterations && !HasFatalFailure(); ++iter) {
    if (rng_(2)) {
      for (int i = 0 ; i < kMaxSize * kMaxSize ; ++i)
        src[i] = kInt13Max;
    } else {
      for (int i = 0 ; i < kMaxSize * kMaxSize ; ++i)
        src[i] = -kInt13Max;
    }

    const int N = rng_(2) ? rng_(kMaxSize * kMaxSize + 1 - kMaxSize) + kMaxSize
                          : rng_(kMaxSize) + 1;

    const uint64_t ref_res = ref_func_(src, N);
    const uint64_t tst_res = tst_func_(src, N);

    ASSERT_EQ(ref_res, tst_res);
  }
}

using std::tr1::make_tuple;

#if HAVE_SSE2
INSTANTIATE_TEST_CASE_P(
    SSE2, SumSquares1DTest,
    ::testing::Values(
        make_tuple(&vpx_sum_squares_i16_c, &vpx_sum_squares_i16_sse2)
    )
);
#endif  // HAVE_SSE2
}  // namespace
