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
#include "test/register_state_check.h"

#include "./vpx_config.h"
#include "./vpx_dsp_rtcd.h"
#include "vpx/vpx_integer.h"

#define MAX_SB_SQUARE (MAX_SB_SIZE * MAX_SB_SIZE)

using libvpx_test::ACMRandom;
using libvpx_test::FunctionEquivalenceTest;

namespace {

static const int kIterations = 1000;
static const int kMaskMax = 64;

typedef unsigned int (*ObmcVarF)(const uint8_t *pre, int pre_stride,
                                 const int32_t *wsrc, const int32_t *mask,
                                 unsigned int *sse);
typedef libvpx_test::FuncParam<ObmcVarF> TestFuncs;

////////////////////////////////////////////////////////////////////////////////
// 8 bit
////////////////////////////////////////////////////////////////////////////////

class ObmcVarianceTest : public FunctionEquivalenceTest<ObmcVarF> {};

TEST_P(ObmcVarianceTest, RandomValues) {
  DECLARE_ALIGNED(32, uint8_t, pre[MAX_SB_SQUARE]);
  DECLARE_ALIGNED(32, int32_t, wsrc[MAX_SB_SQUARE]);
  DECLARE_ALIGNED(32, int32_t, mask[MAX_SB_SQUARE]);

  for (int iter = 0 ; iter < kIterations && !HasFatalFailure() ; ++iter) {
    const int pre_stride = this->rng_(MAX_SB_SIZE + 1);

    for (int i = 0 ; i < MAX_SB_SQUARE ; ++i) {
      pre[i] = this->rng_.Rand8();
      wsrc[i] = this->rng_.Rand8() * this->rng_(kMaskMax * kMaskMax + 1);
      mask[i] = this->rng_(kMaskMax * kMaskMax + 1);
    }

    unsigned int ref_sse, tst_sse;
    const unsigned int ref_res = params_.ref_func(pre, pre_stride, wsrc, mask,
                                                  &ref_sse);
    unsigned int tst_res;
    ASM_REGISTER_STATE_CHECK(
        tst_res = params_.tst_func(pre, pre_stride, wsrc, mask, &tst_sse));

    ASSERT_EQ(ref_res, tst_res);
    ASSERT_EQ(ref_sse, tst_sse);
  }
}

TEST_P(ObmcVarianceTest, ExtremeValues) {
  DECLARE_ALIGNED(32, uint8_t, pre[MAX_SB_SQUARE]);
  DECLARE_ALIGNED(32, int32_t, wsrc[MAX_SB_SQUARE]);
  DECLARE_ALIGNED(32, int32_t, mask[MAX_SB_SQUARE]);

  for (int iter = 0 ; iter < MAX_SB_SIZE && !HasFatalFailure() ; ++iter) {
    const int pre_stride = iter;

    for (int i = 0 ; i < MAX_SB_SQUARE ; ++i) {
      pre[i] = UINT8_MAX;
      wsrc[i] = UINT8_MAX * kMaskMax * kMaskMax;
      mask[i] = kMaskMax * kMaskMax;
    }

    unsigned int ref_sse, tst_sse;
    const unsigned int ref_res = params_.ref_func(pre, pre_stride, wsrc, mask,
                                                  &ref_sse);
    unsigned int tst_res;
    ASM_REGISTER_STATE_CHECK(
        tst_res = params_.tst_func(pre, pre_stride, wsrc, mask, &tst_sse));

    ASSERT_EQ(ref_res, tst_res);
    ASSERT_EQ(ref_sse, tst_sse);
  }
}

#if HAVE_SSE4_1
const ObmcVarianceTest::ParamType sse4_functions[] = {
#if CONFIG_EXT_PARTITION
  TestFuncs(vpx_obmc_variance128x128_c, vpx_obmc_variance128x128_sse4_1),
  TestFuncs(vpx_obmc_variance128x64_c, vpx_obmc_variance128x64_sse4_1),
  TestFuncs(vpx_obmc_variance64x128_c, vpx_obmc_variance64x128_sse4_1),
#endif  // CONFIG_EXT_PARTITION
  TestFuncs(vpx_obmc_variance64x64_c, vpx_obmc_variance64x64_sse4_1),
  TestFuncs(vpx_obmc_variance64x32_c, vpx_obmc_variance64x32_sse4_1),
  TestFuncs(vpx_obmc_variance32x64_c, vpx_obmc_variance32x64_sse4_1),
  TestFuncs(vpx_obmc_variance32x32_c, vpx_obmc_variance32x32_sse4_1),
  TestFuncs(vpx_obmc_variance32x16_c, vpx_obmc_variance32x16_sse4_1),
  TestFuncs(vpx_obmc_variance16x32_c, vpx_obmc_variance16x32_sse4_1),
  TestFuncs(vpx_obmc_variance16x16_c, vpx_obmc_variance16x16_sse4_1),
  TestFuncs(vpx_obmc_variance16x8_c, vpx_obmc_variance16x8_sse4_1),
  TestFuncs(vpx_obmc_variance8x16_c, vpx_obmc_variance8x16_sse4_1),
  TestFuncs(vpx_obmc_variance8x8_c, vpx_obmc_variance8x8_sse4_1),
  TestFuncs(vpx_obmc_variance8x4_c, vpx_obmc_variance8x4_sse4_1),
  TestFuncs(vpx_obmc_variance4x8_c, vpx_obmc_variance4x8_sse4_1),
  TestFuncs(vpx_obmc_variance4x4_c, vpx_obmc_variance4x4_sse4_1)
};

INSTANTIATE_TEST_CASE_P(SSE4_1_C_COMPARE, ObmcVarianceTest,
                        ::testing::ValuesIn(sse4_functions));
#endif  // HAVE_SSE4_1

////////////////////////////////////////////////////////////////////////////////
// High bit-depth
////////////////////////////////////////////////////////////////////////////////

#if CONFIG_VPX_HIGHBITDEPTH
class ObmcVarianceHBDTest : public FunctionEquivalenceTest<ObmcVarF> {};

TEST_P(ObmcVarianceHBDTest, RandomValues) {
  DECLARE_ALIGNED(32, uint16_t, pre[MAX_SB_SQUARE]);
  DECLARE_ALIGNED(32, int32_t, wsrc[MAX_SB_SQUARE]);
  DECLARE_ALIGNED(32, int32_t, mask[MAX_SB_SQUARE]);

  for (int iter = 0 ; iter < kIterations && !HasFatalFailure() ; ++iter) {
    const int pre_stride = this->rng_(MAX_SB_SIZE + 1);

    for (int i = 0 ; i < MAX_SB_SQUARE ; ++i) {
      pre[i] = this->rng_(1 << params_.bit_depth);
      wsrc[i] = this->rng_(1 << params_.bit_depth) *
                this->rng_(kMaskMax * kMaskMax + 1);
      mask[i] = this->rng_(kMaskMax * kMaskMax + 1);
    }

    unsigned int ref_sse, tst_sse;
    const unsigned int ref_res = params_.ref_func(CONVERT_TO_BYTEPTR(pre),
                                                  pre_stride,
                                                  wsrc, mask, &ref_sse);
    unsigned int tst_res;
    ASM_REGISTER_STATE_CHECK(
        tst_res = params_.tst_func(CONVERT_TO_BYTEPTR(pre),
                                   pre_stride, wsrc, mask, &tst_sse));

    ASSERT_EQ(ref_res, tst_res);
    ASSERT_EQ(ref_sse, tst_sse);
  }
}

TEST_P(ObmcVarianceHBDTest, ExtremeValues) {
  DECLARE_ALIGNED(32, uint16_t, pre[MAX_SB_SQUARE]);
  DECLARE_ALIGNED(32, int32_t, wsrc[MAX_SB_SQUARE]);
  DECLARE_ALIGNED(32, int32_t, mask[MAX_SB_SQUARE]);

  for (int iter = 0 ; iter < MAX_SB_SIZE && !HasFatalFailure() ; ++iter) {
    const int pre_stride = iter;

    for (int i = 0 ; i < MAX_SB_SQUARE ; ++i) {
      pre[i] = (1 << params_.bit_depth) - 1;
      wsrc[i] = ((1 << params_.bit_depth) - 1) * kMaskMax * kMaskMax;
      mask[i] = kMaskMax * kMaskMax;
    }

    unsigned int ref_sse, tst_sse;
    const unsigned int ref_res = params_.ref_func(CONVERT_TO_BYTEPTR(pre),
                                                  pre_stride,
                                                  wsrc, mask, &ref_sse);
    unsigned int tst_res;
    ASM_REGISTER_STATE_CHECK(
        tst_res = params_.tst_func(CONVERT_TO_BYTEPTR(pre), pre_stride,
                                   wsrc, mask, &tst_sse));

    ASSERT_EQ(ref_res, tst_res);
    ASSERT_EQ(ref_sse, tst_sse);
  }
}

#if HAVE_SSE4_1
ObmcVarianceHBDTest::ParamType sse4_functions_hbd[] = {
#if CONFIG_EXT_PARTITION
  TestFuncs(vpx_highbd_obmc_variance128x128_c,
            vpx_highbd_obmc_variance128x128_sse4_1, 8),
  TestFuncs(vpx_highbd_obmc_variance128x64_c,
            vpx_highbd_obmc_variance128x64_sse4_1, 8),
  TestFuncs(vpx_highbd_obmc_variance64x128_c,
            vpx_highbd_obmc_variance64x128_sse4_1, 8),
#endif  // CONFIG_EXT_PARTITION
  TestFuncs(vpx_highbd_obmc_variance64x64_c,
            vpx_highbd_obmc_variance64x64_sse4_1, 8),
  TestFuncs(vpx_highbd_obmc_variance64x32_c,
            vpx_highbd_obmc_variance64x32_sse4_1, 8),
  TestFuncs(vpx_highbd_obmc_variance32x64_c,
            vpx_highbd_obmc_variance32x64_sse4_1, 8),
  TestFuncs(vpx_highbd_obmc_variance32x32_c,
            vpx_highbd_obmc_variance32x32_sse4_1, 8),
  TestFuncs(vpx_highbd_obmc_variance32x16_c,
            vpx_highbd_obmc_variance32x16_sse4_1, 8),
  TestFuncs(vpx_highbd_obmc_variance16x32_c,
            vpx_highbd_obmc_variance16x32_sse4_1, 8),
  TestFuncs(vpx_highbd_obmc_variance16x16_c,
            vpx_highbd_obmc_variance16x16_sse4_1, 8),
  TestFuncs(vpx_highbd_obmc_variance16x8_c,
            vpx_highbd_obmc_variance16x8_sse4_1, 8),
  TestFuncs(vpx_highbd_obmc_variance8x16_c,
            vpx_highbd_obmc_variance8x16_sse4_1, 8),
  TestFuncs(vpx_highbd_obmc_variance8x8_c,
            vpx_highbd_obmc_variance8x8_sse4_1, 8),
  TestFuncs(vpx_highbd_obmc_variance8x4_c,
            vpx_highbd_obmc_variance8x4_sse4_1, 8),
  TestFuncs(vpx_highbd_obmc_variance4x8_c,
            vpx_highbd_obmc_variance4x8_sse4_1, 8),
  TestFuncs(vpx_highbd_obmc_variance4x4_c,
            vpx_highbd_obmc_variance4x4_sse4_1, 8),
#if CONFIG_EXT_PARTITION
  TestFuncs(vpx_highbd_10_obmc_variance128x128_c,
            vpx_highbd_10_obmc_variance128x128_sse4_1, 10),
  TestFuncs(vpx_highbd_10_obmc_variance128x64_c,
            vpx_highbd_10_obmc_variance128x64_sse4_1, 10),
  TestFuncs(vpx_highbd_10_obmc_variance64x128_c,
            vpx_highbd_10_obmc_variance64x128_sse4_1, 10),
#endif  // CONFIG_EXT_PARTITION
  TestFuncs(vpx_highbd_10_obmc_variance64x64_c,
            vpx_highbd_10_obmc_variance64x64_sse4_1, 10),
  TestFuncs(vpx_highbd_10_obmc_variance64x32_c,
            vpx_highbd_10_obmc_variance64x32_sse4_1, 10),
  TestFuncs(vpx_highbd_10_obmc_variance32x64_c,
            vpx_highbd_10_obmc_variance32x64_sse4_1, 10),
  TestFuncs(vpx_highbd_10_obmc_variance32x32_c,
            vpx_highbd_10_obmc_variance32x32_sse4_1, 10),
  TestFuncs(vpx_highbd_10_obmc_variance32x16_c,
            vpx_highbd_10_obmc_variance32x16_sse4_1, 10),
  TestFuncs(vpx_highbd_10_obmc_variance16x32_c,
            vpx_highbd_10_obmc_variance16x32_sse4_1, 10),
  TestFuncs(vpx_highbd_10_obmc_variance16x16_c,
            vpx_highbd_10_obmc_variance16x16_sse4_1, 10),
  TestFuncs(vpx_highbd_10_obmc_variance16x8_c,
            vpx_highbd_10_obmc_variance16x8_sse4_1, 10),
  TestFuncs(vpx_highbd_10_obmc_variance8x16_c,
            vpx_highbd_10_obmc_variance8x16_sse4_1, 10),
  TestFuncs(vpx_highbd_10_obmc_variance8x8_c,
            vpx_highbd_10_obmc_variance8x8_sse4_1, 10),
  TestFuncs(vpx_highbd_10_obmc_variance8x4_c,
            vpx_highbd_10_obmc_variance8x4_sse4_1, 10),
  TestFuncs(vpx_highbd_10_obmc_variance4x8_c,
            vpx_highbd_10_obmc_variance4x8_sse4_1, 10),
  TestFuncs(vpx_highbd_10_obmc_variance4x4_c,
            vpx_highbd_10_obmc_variance4x4_sse4_1, 10),
#if CONFIG_EXT_PARTITION
  TestFuncs(vpx_highbd_12_obmc_variance128x128_c,
            vpx_highbd_12_obmc_variance128x128_sse4_1, 12),
  TestFuncs(vpx_highbd_12_obmc_variance128x64_c,
            vpx_highbd_12_obmc_variance128x64_sse4_1, 12),
  TestFuncs(vpx_highbd_12_obmc_variance64x128_c,
            vpx_highbd_12_obmc_variance64x128_sse4_1, 12),
#endif  // CONFIG_EXT_PARTITION
  TestFuncs(vpx_highbd_12_obmc_variance64x64_c,
            vpx_highbd_12_obmc_variance64x64_sse4_1, 12),
  TestFuncs(vpx_highbd_12_obmc_variance64x32_c,
            vpx_highbd_12_obmc_variance64x32_sse4_1, 12),
  TestFuncs(vpx_highbd_12_obmc_variance32x64_c,
            vpx_highbd_12_obmc_variance32x64_sse4_1, 12),
  TestFuncs(vpx_highbd_12_obmc_variance32x32_c,
            vpx_highbd_12_obmc_variance32x32_sse4_1, 12),
  TestFuncs(vpx_highbd_12_obmc_variance32x16_c,
            vpx_highbd_12_obmc_variance32x16_sse4_1, 12),
  TestFuncs(vpx_highbd_12_obmc_variance16x32_c,
            vpx_highbd_12_obmc_variance16x32_sse4_1, 12),
  TestFuncs(vpx_highbd_12_obmc_variance16x16_c,
            vpx_highbd_12_obmc_variance16x16_sse4_1, 12),
  TestFuncs(vpx_highbd_12_obmc_variance16x8_c,
            vpx_highbd_12_obmc_variance16x8_sse4_1, 12),
  TestFuncs(vpx_highbd_12_obmc_variance8x16_c,
            vpx_highbd_12_obmc_variance8x16_sse4_1, 12),
  TestFuncs(vpx_highbd_12_obmc_variance8x8_c,
            vpx_highbd_12_obmc_variance8x8_sse4_1, 12),
  TestFuncs(vpx_highbd_12_obmc_variance8x4_c,
            vpx_highbd_12_obmc_variance8x4_sse4_1, 12),
  TestFuncs(vpx_highbd_12_obmc_variance4x8_c,
            vpx_highbd_12_obmc_variance4x8_sse4_1, 12),
  TestFuncs(vpx_highbd_12_obmc_variance4x4_c,
            vpx_highbd_12_obmc_variance4x4_sse4_1, 12)
};

INSTANTIATE_TEST_CASE_P(SSE4_1_C_COMPARE, ObmcVarianceHBDTest,
                        ::testing::ValuesIn(sse4_functions_hbd));
#endif  // HAVE_SSE4_1
#endif  // CONFIG_VPX_HIGHBITDEPTH
}  // namespace
