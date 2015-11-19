/*
 *  Copyright (c) 2015 The WebM project authors. All Rights Reserved.
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
#include "test/acm_random.h"
#include "test/clear_system_state.h"
#include "test/register_state_check.h"
#include "test/util.h"

#include "./vpx_config.h"
#include "./vp9_rtcd.h"
#include "vp9/common/vp9_entropy.h"
#include "vpx/vpx_integer.h"
#include "vp9/common/vp9_filter.h"

#define MAX_SIZE 64

using libvpx_test::ACMRandom;

namespace {
const int number_of_iterations = 500;

typedef unsigned int (*MaskedVarianceFunc)(const uint8_t *a, int a_stride,
                                           const uint8_t *b, int b_stride,
                                           const uint8_t *m, int m_stride,
                                           unsigned int *sse);

typedef std::tr1::tuple<MaskedVarianceFunc,
                        MaskedVarianceFunc> MaskedVarianceParam;

class MaskedVarianceTest :
  public ::testing::TestWithParam<MaskedVarianceParam> {
 public:
  virtual ~MaskedVarianceTest() {}
  virtual void SetUp() {
    opt_func_ = GET_PARAM(0);
    ref_func_ = GET_PARAM(1);
  }

  virtual void TearDown() { libvpx_test::ClearSystemState(); }

 protected:
  MaskedVarianceFunc opt_func_;
  MaskedVarianceFunc ref_func_;
};

TEST_P(MaskedVarianceTest, OperationCheck) {
  unsigned int ref_ret, opt_ret;
  unsigned int ref_sse, opt_sse;
  ACMRandom rnd(ACMRandom::DeterministicSeed());
  DECLARE_ALIGNED_ARRAY(16, uint8_t,  src_ptr, MAX_SIZE*MAX_SIZE);
  DECLARE_ALIGNED_ARRAY(16, uint8_t,  ref_ptr, MAX_SIZE*MAX_SIZE);
  DECLARE_ALIGNED_ARRAY(16, uint8_t,  msk_ptr, MAX_SIZE*MAX_SIZE);
  int err_count = 0;
  int first_failure = -1;
  int src_stride = MAX_SIZE;
  int ref_stride = MAX_SIZE;
  int msk_stride = MAX_SIZE;

  for (int i = 0; i < number_of_iterations; ++i) {
    for (int j = 0; j < MAX_SIZE*MAX_SIZE; j++) {
      src_ptr[j] = rnd.Rand8();
      ref_ptr[j] = rnd.Rand8();
      msk_ptr[j] = rnd(65);
    }

    ref_ret = ref_func_(src_ptr, src_stride,
                        ref_ptr, ref_stride,
                        msk_ptr, msk_stride,
                        &ref_sse);
    ASM_REGISTER_STATE_CHECK(opt_ret = opt_func_(src_ptr, src_stride,
                                                 ref_ptr, ref_stride,
                                                 msk_ptr, msk_stride,
                                                 &opt_sse));

    if (opt_ret != ref_ret || opt_sse != ref_sse) {
      err_count++;
      if (first_failure == -1)
        first_failure = i;
    }
  }

  EXPECT_EQ(0, err_count)
  << "Error: Masked Variance Test OperationCheck,"
  << "C output doesn't match SSSE3 output. "
  << "First failed at test case " << first_failure;
}

TEST_P(MaskedVarianceTest, ExtremeValues) {
  unsigned int ref_ret, opt_ret;
  unsigned int ref_sse, opt_sse;
  ACMRandom rnd(ACMRandom::DeterministicSeed());
  DECLARE_ALIGNED_ARRAY(16, uint8_t,  src_ptr, MAX_SIZE*MAX_SIZE);
  DECLARE_ALIGNED_ARRAY(16, uint8_t,  ref_ptr, MAX_SIZE*MAX_SIZE);
  DECLARE_ALIGNED_ARRAY(16, uint8_t,  msk_ptr, MAX_SIZE*MAX_SIZE);
  int err_count = 0;
  int first_failure = -1;
  int src_stride = MAX_SIZE;
  int ref_stride = MAX_SIZE;
  int msk_stride = MAX_SIZE;

  for (int i = 0; i < 8; ++i) {
    memset(src_ptr, (i & 0x1) ? 255 : 0, MAX_SIZE*MAX_SIZE);
    memset(ref_ptr, (i & 0x2) ? 255 : 0, MAX_SIZE*MAX_SIZE);
    memset(msk_ptr, (i & 0x4) ?  64 : 0, MAX_SIZE*MAX_SIZE);

    ref_ret = ref_func_(src_ptr, src_stride,
                        ref_ptr, ref_stride,
                        msk_ptr, msk_stride,
                        &ref_sse);
    ASM_REGISTER_STATE_CHECK(opt_ret = opt_func_(src_ptr, src_stride,
                                                 ref_ptr, ref_stride,
                                                 msk_ptr, msk_stride,
                                                 &opt_sse));

    if (opt_ret != ref_ret || opt_sse != ref_sse) {
      err_count++;
      if (first_failure == -1)
        first_failure = i;
    }
  }

  EXPECT_EQ(0, err_count)
  << "Error: Masked Variance Test ExtremeValues,"
  << "C output doesn't match SSSE3 output. "
  << "First failed at test case " << first_failure;
}

typedef unsigned int (*MaskedSubPixelVarianceFunc)(
    const uint8_t *a, int a_stride,
    int xoffset, int  yoffset,
    const uint8_t *b, int b_stride,
    const uint8_t *m, int m_stride,
    unsigned int *sse);

typedef std::tr1::tuple<MaskedSubPixelVarianceFunc,
                        MaskedSubPixelVarianceFunc> MaskedSubPixelVarianceParam;

class MaskedSubPixelVarianceTest :
  public ::testing::TestWithParam<MaskedSubPixelVarianceParam> {
 public:
  virtual ~MaskedSubPixelVarianceTest() {}
  virtual void SetUp() {
    opt_func_ = GET_PARAM(0);
    ref_func_ = GET_PARAM(1);
  }

  virtual void TearDown() { libvpx_test::ClearSystemState(); }

 protected:
  MaskedSubPixelVarianceFunc opt_func_;
  MaskedSubPixelVarianceFunc ref_func_;
};

TEST_P(MaskedSubPixelVarianceTest, OperationCheck) {
  unsigned int ref_ret, opt_ret;
  unsigned int ref_sse, opt_sse;
  ACMRandom rnd(ACMRandom::DeterministicSeed());
  DECLARE_ALIGNED_ARRAY(16, uint8_t,  src_ptr, (MAX_SIZE+1)*(MAX_SIZE+1));
  DECLARE_ALIGNED_ARRAY(16, uint8_t,  ref_ptr, (MAX_SIZE+1)*(MAX_SIZE+1));
  DECLARE_ALIGNED_ARRAY(16, uint8_t,  msk_ptr, (MAX_SIZE+1)*(MAX_SIZE+1));
  int err_count = 0;
  int first_failure = -1;
  int src_stride = (MAX_SIZE+1);
  int ref_stride = (MAX_SIZE+1);
  int msk_stride = (MAX_SIZE+1);
  int xoffset;
  int yoffset;

  for (int i = 0; i < number_of_iterations; ++i) {
    int xoffsets[] = {0, 8, rnd(SUBPEL_SHIFTS)};
    int yoffsets[] = {0, 8, rnd(SUBPEL_SHIFTS)};
    for (int j = 0; j < (MAX_SIZE+1)*(MAX_SIZE+1); j++) {
      src_ptr[j] = rnd.Rand8();
      ref_ptr[j] = rnd.Rand8();
      msk_ptr[j] = rnd(65);
    }
    for (int k = 0; k < 3; k++) {
      xoffset = xoffsets[k];
      for (int l = 0; l < 3; l++) {
        xoffset = xoffsets[k];
        yoffset = yoffsets[l];

        ref_ret = ref_func_(src_ptr, src_stride,
                            xoffset, yoffset,
                            ref_ptr, ref_stride,
                            msk_ptr, msk_stride,
                            &ref_sse);
        ASM_REGISTER_STATE_CHECK(opt_ret = opt_func_(src_ptr, src_stride,
                                                    xoffset, yoffset,
                                                    ref_ptr, ref_stride,
                                                    msk_ptr, msk_stride,
                                                    &opt_sse));

        if (opt_ret != ref_ret || opt_sse != ref_sse) {
        err_count++;
        if (first_failure == -1)
            first_failure = i;
        }
      }
    }
  }

  EXPECT_EQ(0, err_count)
    << "Error: Masked Sub Pixel Variance Test OperationCheck,"
    << "C output doesn't match SSSE3 output. "
    << "First failed at test case " << first_failure;
}

TEST_P(MaskedSubPixelVarianceTest, ExtremeValues) {
  unsigned int ref_ret, opt_ret;
  unsigned int ref_sse, opt_sse;
  ACMRandom rnd(ACMRandom::DeterministicSeed());
  DECLARE_ALIGNED_ARRAY(16, uint8_t,  src_ptr, (MAX_SIZE+1)*(MAX_SIZE+1));
  DECLARE_ALIGNED_ARRAY(16, uint8_t,  ref_ptr, (MAX_SIZE+1)*(MAX_SIZE+1));
  DECLARE_ALIGNED_ARRAY(16, uint8_t,  msk_ptr, (MAX_SIZE+1)*(MAX_SIZE+1));
  int first_failure_x = -1;
  int first_failure_y = -1;
  int err_count = 0;
  int first_failure = -1;
  int src_stride = (MAX_SIZE+1);
  int ref_stride = (MAX_SIZE+1);
  int msk_stride = (MAX_SIZE+1);

  for (int xoffset = 0 ; xoffset < SUBPEL_SHIFTS ; xoffset++) {
    for (int yoffset = 0 ; yoffset < SUBPEL_SHIFTS ; yoffset++) {
      for (int i = 0; i < 8; ++i) {
        memset(src_ptr, (i & 0x1) ? 255 : 0, (MAX_SIZE+1)*(MAX_SIZE+1));
        memset(ref_ptr, (i & 0x2) ? 255 : 0, (MAX_SIZE+1)*(MAX_SIZE+1));
        memset(msk_ptr, (i & 0x4) ?  64 : 0, (MAX_SIZE+1)*(MAX_SIZE+1));

        ref_ret = ref_func_(src_ptr, src_stride,
                            xoffset, yoffset,
                            ref_ptr, ref_stride,
                            msk_ptr, msk_stride,
                            &ref_sse);
        ASM_REGISTER_STATE_CHECK(opt_ret = opt_func_(src_ptr, src_stride,
                                                     xoffset, yoffset,
                                                     ref_ptr, ref_stride,
                                                     msk_ptr, msk_stride,
                                                     &opt_sse));

        if (opt_ret != ref_ret || opt_sse != ref_sse) {
          err_count++;
          if (first_failure == -1) {
            first_failure = i;
            first_failure_x = xoffset;
            first_failure_y = yoffset;
          }
        }
      }
    }
  }

  EXPECT_EQ(0, err_count)
  << "Error: Masked Variance Test ExtremeValues,"
  << "C output doesn't match SSSE3 output. "
  << "First failed at test case " << first_failure
  << " x_offset = " << first_failure_x
  << " y_offset = " << first_failure_y;
}

using std::tr1::make_tuple;

#if HAVE_SSSE3
INSTANTIATE_TEST_CASE_P(
  SSSE3_C_COMPARE, MaskedVarianceTest,
  ::testing::Values(
    make_tuple(&vp9_masked_variance64x64_ssse3,
               &vp9_masked_variance64x64_c),
    make_tuple(&vp9_masked_variance64x32_ssse3,
               &vp9_masked_variance64x32_c),
    make_tuple(&vp9_masked_variance32x64_ssse3,
               &vp9_masked_variance32x64_c),
    make_tuple(&vp9_masked_variance32x32_ssse3,
               &vp9_masked_variance32x32_c),
    make_tuple(&vp9_masked_variance32x16_ssse3,
               &vp9_masked_variance32x16_c),
    make_tuple(&vp9_masked_variance16x32_ssse3,
               &vp9_masked_variance16x32_c),
    make_tuple(&vp9_masked_variance16x16_ssse3,
               &vp9_masked_variance16x16_c),
    make_tuple(&vp9_masked_variance16x8_ssse3,
               &vp9_masked_variance16x8_c),
    make_tuple(&vp9_masked_variance8x16_ssse3,
               &vp9_masked_variance8x16_c),
    make_tuple(&vp9_masked_variance8x8_ssse3,
               &vp9_masked_variance8x8_c),
    make_tuple(&vp9_masked_variance8x4_ssse3,
               &vp9_masked_variance8x4_c),
    make_tuple(&vp9_masked_variance4x8_ssse3,
               &vp9_masked_variance4x8_c),
    make_tuple(&vp9_masked_variance4x4_ssse3,
               &vp9_masked_variance4x4_c)));

INSTANTIATE_TEST_CASE_P(
  SSSE3_C_COMPARE, MaskedSubPixelVarianceTest,
  ::testing::Values(
    make_tuple(&vp9_masked_sub_pixel_variance64x64_ssse3,
              &vp9_masked_sub_pixel_variance64x64_c),
    make_tuple(&vp9_masked_sub_pixel_variance64x32_ssse3,
              &vp9_masked_sub_pixel_variance64x32_c),
    make_tuple(&vp9_masked_sub_pixel_variance32x64_ssse3,
              &vp9_masked_sub_pixel_variance32x64_c),
    make_tuple(&vp9_masked_sub_pixel_variance32x32_ssse3,
              &vp9_masked_sub_pixel_variance32x32_c),
    make_tuple(&vp9_masked_sub_pixel_variance32x16_ssse3,
              &vp9_masked_sub_pixel_variance32x16_c),
    make_tuple(&vp9_masked_sub_pixel_variance16x32_ssse3,
              &vp9_masked_sub_pixel_variance16x32_c),
    make_tuple(&vp9_masked_sub_pixel_variance16x16_ssse3,
              &vp9_masked_sub_pixel_variance16x16_c),
    make_tuple(&vp9_masked_sub_pixel_variance16x8_ssse3,
              &vp9_masked_sub_pixel_variance16x8_c),
    make_tuple(&vp9_masked_sub_pixel_variance8x16_ssse3,
              &vp9_masked_sub_pixel_variance8x16_c),
    make_tuple(&vp9_masked_sub_pixel_variance8x8_ssse3,
              &vp9_masked_sub_pixel_variance8x8_c),
    make_tuple(&vp9_masked_sub_pixel_variance8x4_ssse3,
              &vp9_masked_sub_pixel_variance8x4_c),
    make_tuple(&vp9_masked_sub_pixel_variance4x8_ssse3,
              &vp9_masked_sub_pixel_variance4x8_c),
    make_tuple(&vp9_masked_sub_pixel_variance4x4_ssse3,
              &vp9_masked_sub_pixel_variance4x4_c)));
#endif  // HAVE_SSSE3
}  // namespace
