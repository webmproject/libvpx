/*
 *  Copyright (c) 2014 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include "third_party/googletest/src/include/gtest/gtest.h"

#include "./vpx_config.h"
#include "./vpx_dsp_rtcd.h"
#include "vpx_ports/mem.h"

#include "test/array_utils.h"
#include "test/assertion_helpers.h"
#include "test/function_equivalence_test.h"
#include "test/randomise.h"
#include "test/register_state_check.h"
#include "test/snapshot.h"

using libvpx_test::FunctionEquivalenceTest;
using libvpx_test::Snapshot;
using libvpx_test::Randomise;
using libvpx_test::array_utils::arraySet;
using libvpx_test::assertion_helpers::ArraysEq;

namespace {

static const int16_t int13_max = (1<<12) - 1;

//////////////////////////////////////////////////////////////////////////////
// 2D version
//////////////////////////////////////////////////////////////////////////////

typedef uint64_t (*F2D)(const int16_t *src, int stride, uint32_t size);

class SumSquares2DTest : public FunctionEquivalenceTest<F2D> {
 protected:
  void Common() {
    const int sizelog2 = randomise.uniform<int>(2, 8);

    const uint32_t size = 1 << sizelog2;
    const int stride = 1 << randomise.uniform<int>(sizelog2, 9);

    snapshot(src);

    uint64_t ref_res, tst_res;

    ref_res = ref_func_(src, stride, size);
    ASM_REGISTER_STATE_CHECK(tst_res = tst_func_(src, stride, size));

    ASSERT_EQ(ref_res, tst_res);

    ASSERT_TRUE(ArraysEq(snapshot.get(src), src));
  }

  Snapshot snapshot;
  Randomise randomise;

  DECLARE_ALIGNED(16, int16_t, src[256*256]);
};

TEST_P(SumSquares2DTest, RandomValues) {
  for (int i = 0 ; i < 10000 && !HasFatalFailure(); i++) {
    randomise(src, -int13_max, int13_max + 1);

    Common();
  }
}

TEST_P(SumSquares2DTest, ExtremeValues) {
  for (int i = 0 ; i < 10000 && !HasFatalFailure(); i++) {
    if (randomise.uniform<bool>())
      arraySet(src, int13_max);
    else
      arraySet(src, -int13_max);

    Common();
  }
}
using std::tr1::make_tuple;

#if HAVE_SSE2
INSTANTIATE_TEST_CASE_P(
    SSE2, SumSquares2DTest,
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
  void Common() {
    const int N = randomise.uniform<int>(1, 256*256-1);

    snapshot(src);

    uint64_t ref_res, tst_res;

    ref_res = ref_func_(src, N);
    ASM_REGISTER_STATE_CHECK(tst_res = tst_func_(src, N));

    ASSERT_EQ(ref_res, tst_res);

    ASSERT_TRUE(ArraysEq(snapshot.get(src), src));
  }

  Snapshot snapshot;
  Randomise randomise;

  DECLARE_ALIGNED(16, int16_t, src[256*256]);
};

TEST_P(SumSquares1DTest, RandomValues) {
  for (int i = 0 ; i < 10000 && !HasFatalFailure(); i++) {
    randomise(src, -int13_max, int13_max+1);

    Common();
  }
}

TEST_P(SumSquares1DTest, ExtremeValues) {
  for (int i = 0 ; i < 10000 && !HasFatalFailure(); i++) {
    if (randomise.uniform<bool>())
      arraySet(src, int13_max);
    else
      arraySet(src, -int13_max);

    Common();
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
