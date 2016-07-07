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
#include "test/acm_random.h"

#include "test/function_equivalence_test.h"

#include "./vpx_config.h"
#include "./vpx_dsp_rtcd.h"
#include "vpx/vpx_integer.h"

#define MAX_SB_SQUARE (MAX_SB_SIZE * MAX_SB_SIZE)

using std::tr1::make_tuple;

using libvpx_test::ACMRandom;
using libvpx_test::FunctionEquivalenceTest;

namespace {

static const int kIterations = 1000;
static const int kMaskMax = 64;

typedef unsigned int (*ObmcSadF)(const uint8_t *ref, int ref_stride,
                                 const int32_t *wsrc, const int32_t *mask);

////////////////////////////////////////////////////////////////////////////////
// 8 bit
////////////////////////////////////////////////////////////////////////////////

class ObmcSadTest : public FunctionEquivalenceTest<ObmcSadF> {
 public:
  ObmcSadTest() : rng_(ACMRandom::DeterministicSeed()) {}

 protected:
  ACMRandom rng_;
};

TEST_P(ObmcSadTest, RandomValues) {
  DECLARE_ALIGNED(32, uint8_t, ref[MAX_SB_SQUARE]);
  DECLARE_ALIGNED(32, int32_t, wsrc[MAX_SB_SQUARE]);
  DECLARE_ALIGNED(32, int32_t, mask[MAX_SB_SQUARE]);

  for (int iter = 0 ; iter < kIterations && !HasFatalFailure() ; ++iter) {
    const int ref_stride = rng_(MAX_SB_SIZE + 1);

    for (int i = 0 ; i < MAX_SB_SQUARE ; ++i) {
      ref[i] = rng_.Rand8();
      wsrc[i] = rng_.Rand8() * rng_(kMaskMax * kMaskMax + 1);
      mask[i] = rng_(kMaskMax * kMaskMax + 1);
    }

    const unsigned int ref_res = ref_func_(ref, ref_stride, wsrc, mask);
    const unsigned int tst_res = tst_func_(ref, ref_stride, wsrc, mask);

    ASSERT_EQ(ref_res, tst_res);
  }
}

TEST_P(ObmcSadTest, ExtremeValues) {
  DECLARE_ALIGNED(32, uint8_t, ref[MAX_SB_SQUARE]);
  DECLARE_ALIGNED(32, int32_t, wsrc[MAX_SB_SQUARE]);
  DECLARE_ALIGNED(32, int32_t, mask[MAX_SB_SQUARE]);

  for (int iter = 0 ; iter < MAX_SB_SIZE && !HasFatalFailure() ; ++iter) {
    const int ref_stride = iter;

    for (int i = 0 ; i < MAX_SB_SQUARE ; ++i) {
      ref[i] = UINT8_MAX;
      wsrc[i] = UINT8_MAX * kMaskMax * kMaskMax;
      mask[i] = kMaskMax * kMaskMax;
    }

    const unsigned int ref_res = ref_func_(ref, ref_stride, wsrc, mask);
    const unsigned int tst_res = tst_func_(ref, ref_stride, wsrc, mask);

    ASSERT_EQ(ref_res, tst_res);
  }
}

#if HAVE_SSE4_1
const ObmcSadTest::ParamType sse4_functions[] = {
#if CONFIG_EXT_PARTITION
  make_tuple(vpx_obmc_sad128x128_c, vpx_obmc_sad128x128_sse4_1),
  make_tuple(vpx_obmc_sad128x64_c, vpx_obmc_sad128x64_sse4_1),
  make_tuple(vpx_obmc_sad64x128_c, vpx_obmc_sad64x128_sse4_1),
#endif  // CONFIG_EXT_PARTITION
  make_tuple(vpx_obmc_sad64x64_c, vpx_obmc_sad64x64_sse4_1),
  make_tuple(vpx_obmc_sad64x32_c, vpx_obmc_sad64x32_sse4_1),
  make_tuple(vpx_obmc_sad32x64_c, vpx_obmc_sad32x64_sse4_1),
  make_tuple(vpx_obmc_sad32x32_c, vpx_obmc_sad32x32_sse4_1),
  make_tuple(vpx_obmc_sad32x16_c, vpx_obmc_sad32x16_sse4_1),
  make_tuple(vpx_obmc_sad16x32_c, vpx_obmc_sad16x32_sse4_1),
  make_tuple(vpx_obmc_sad16x16_c, vpx_obmc_sad16x16_sse4_1),
  make_tuple(vpx_obmc_sad16x8_c, vpx_obmc_sad16x8_sse4_1),
  make_tuple(vpx_obmc_sad8x16_c, vpx_obmc_sad8x16_sse4_1),
  make_tuple(vpx_obmc_sad8x8_c, vpx_obmc_sad8x8_sse4_1),
  make_tuple(vpx_obmc_sad8x4_c, vpx_obmc_sad8x4_sse4_1),
  make_tuple(vpx_obmc_sad4x8_c, vpx_obmc_sad4x8_sse4_1),
  make_tuple(vpx_obmc_sad4x4_c, vpx_obmc_sad4x4_sse4_1)
};

INSTANTIATE_TEST_CASE_P(SSE4_1_C_COMPARE, ObmcSadTest,
                        ::testing::ValuesIn(sse4_functions));
#endif  // HAVE_SSE4_1

////////////////////////////////////////////////////////////////////////////////
// High bit-depth
////////////////////////////////////////////////////////////////////////////////

#if CONFIG_VP9_HIGHBITDEPTH
class ObmcSadHBDTest : public FunctionEquivalenceTest<ObmcSadF> {
 public:
  ObmcSadHBDTest() : rng_(ACMRandom::DeterministicSeed()) {}

 protected:
  ACMRandom rng_;
};

TEST_P(ObmcSadHBDTest, RandomValues) {
  DECLARE_ALIGNED(32, uint16_t, ref[MAX_SB_SQUARE]);
  DECLARE_ALIGNED(32, int32_t, wsrc[MAX_SB_SQUARE]);
  DECLARE_ALIGNED(32, int32_t, mask[MAX_SB_SQUARE]);

  for (int iter = 0 ; iter < kIterations && !HasFatalFailure() ; ++iter) {
    const int ref_stride = rng_(MAX_SB_SIZE + 1);

    for (int i = 0 ; i < MAX_SB_SQUARE ; ++i) {
      ref[i] = rng_(1<<12);
      wsrc[i] = rng_(1<<12) * rng_(kMaskMax * kMaskMax + 1);
      mask[i] = rng_(kMaskMax * kMaskMax + 1);
    }

    const unsigned int ref_res = ref_func_(CONVERT_TO_BYTEPTR(ref), ref_stride,
                                           wsrc, mask);
    const unsigned int tst_res = tst_func_(CONVERT_TO_BYTEPTR(ref), ref_stride,
                                           wsrc, mask);

    ASSERT_EQ(ref_res, tst_res);
  }
}

TEST_P(ObmcSadHBDTest, ExtremeValues) {
  DECLARE_ALIGNED(32, uint16_t, ref[MAX_SB_SQUARE]);
  DECLARE_ALIGNED(32, int32_t, wsrc[MAX_SB_SQUARE]);
  DECLARE_ALIGNED(32, int32_t, mask[MAX_SB_SQUARE]);

  for (int iter = 0 ; iter < MAX_SB_SIZE && !HasFatalFailure() ; ++iter) {
    const int ref_stride = iter;

    for (int i = 0 ; i < MAX_SB_SQUARE ; ++i) {
      ref[i] = (1 << 12) - 1;
      wsrc[i] = ((1 << 12) - 1) * kMaskMax * kMaskMax;
      mask[i] = kMaskMax * kMaskMax;
    }

    const unsigned int ref_res = ref_func_(CONVERT_TO_BYTEPTR(ref), ref_stride,
                                           wsrc, mask);
    const unsigned int tst_res = tst_func_(CONVERT_TO_BYTEPTR(ref), ref_stride,
                                           wsrc, mask);

    ASSERT_EQ(ref_res, tst_res);
  }
}

#if HAVE_SSE4_1
ObmcSadHBDTest::ParamType sse4_functions_hbd[] = {
#if CONFIG_EXT_PARTITION
  make_tuple(vpx_highbd_obmc_sad128x128_c, vpx_highbd_obmc_sad128x128_sse4_1),
  make_tuple(vpx_highbd_obmc_sad128x64_c, vpx_highbd_obmc_sad128x64_sse4_1),
  make_tuple(vpx_highbd_obmc_sad64x128_c, vpx_highbd_obmc_sad64x128_sse4_1),
#endif  // CONFIG_EXT_PARTITION
  make_tuple(vpx_highbd_obmc_sad64x64_c, vpx_highbd_obmc_sad64x64_sse4_1),
  make_tuple(vpx_highbd_obmc_sad64x32_c, vpx_highbd_obmc_sad64x32_sse4_1),
  make_tuple(vpx_highbd_obmc_sad32x64_c, vpx_highbd_obmc_sad32x64_sse4_1),
  make_tuple(vpx_highbd_obmc_sad32x32_c, vpx_highbd_obmc_sad32x32_sse4_1),
  make_tuple(vpx_highbd_obmc_sad32x16_c, vpx_highbd_obmc_sad32x16_sse4_1),
  make_tuple(vpx_highbd_obmc_sad16x32_c, vpx_highbd_obmc_sad16x32_sse4_1),
  make_tuple(vpx_highbd_obmc_sad16x16_c, vpx_highbd_obmc_sad16x16_sse4_1),
  make_tuple(vpx_highbd_obmc_sad16x8_c, vpx_highbd_obmc_sad16x8_sse4_1),
  make_tuple(vpx_highbd_obmc_sad8x16_c, vpx_highbd_obmc_sad8x16_sse4_1),
  make_tuple(vpx_highbd_obmc_sad8x8_c, vpx_highbd_obmc_sad8x8_sse4_1),
  make_tuple(vpx_highbd_obmc_sad8x4_c, vpx_highbd_obmc_sad8x4_sse4_1),
  make_tuple(vpx_highbd_obmc_sad4x8_c, vpx_highbd_obmc_sad4x8_sse4_1),
  make_tuple(vpx_highbd_obmc_sad4x4_c, vpx_highbd_obmc_sad4x4_sse4_1)
};

INSTANTIATE_TEST_CASE_P(SSE4_1_C_COMPARE, ObmcSadHBDTest,
                        ::testing::ValuesIn(sse4_functions_hbd));
#endif  // HAVE_SSE4_1
#endif  // CONFIG_VP9_HIGHBITDEPTH
}  // namespace
