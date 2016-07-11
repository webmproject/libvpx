/*
 *  Copyright (c) 2016 The WebM project authors. All Rights Reserved.
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
#include "test/register_state_check.h"

#include "test/function_equivalence_test.h"

#include "./vpx_config.h"
#include "./vpx_dsp_rtcd.h"
#include "vpx/vpx_integer.h"

#include "./vp10_rtcd.h"

#include "test/acm_random.h"
#include "vp10/common/enums.h"

#include "vpx_dsp/blend.h"

using libvpx_test::ACMRandom;
using libvpx_test::FunctionEquivalenceTest;
using std::tr1::make_tuple;

namespace {

template<typename F, typename T>
class BlendA64Mask1DTest : public FunctionEquivalenceTest<F> {
 public:
  static const int kIterations = 10000;
  static const int kMaxWidth = MAX_SB_SIZE * 5;  // * 5 to cover longer strides
  static const int kMaxHeight = MAX_SB_SIZE;
  static const int kBufSize = kMaxWidth * kMaxHeight;
  static const int kMaxMaskWidth = 2 * MAX_SB_SIZE;
  static const int kMaxMaskSize = kMaxMaskWidth;

  BlendA64Mask1DTest() : rng_(ACMRandom::DeterministicSeed()) {}

  virtual ~BlendA64Mask1DTest() {}

  virtual void Execute(T *p_src0, T *p_src1) = 0;

  void Common() {
    w_ = 1 << rng_(MAX_SB_SIZE_LOG2 + 1);
    h_ = 1 << rng_(MAX_SB_SIZE_LOG2 + 1);

    dst_offset_ = rng_(33);
    dst_stride_ = rng_(kMaxWidth + 1 - w_) + w_;

    src0_offset_ = rng_(33);
    src0_stride_ = rng_(kMaxWidth + 1 - w_) + w_;

    src1_offset_ = rng_(33);
    src1_stride_ = rng_(kMaxWidth + 1 - w_) + w_;

    T *p_src0;
    T *p_src1;

    switch (rng_(3)) {
      case 0:   // Separate sources
        p_src0 = src0_;
        p_src1 = src1_;
        break;
      case 1:   // src0 == dst
        p_src0 = dst_tst_;
        src0_stride_ = dst_stride_;
        src0_offset_ = dst_offset_;
        p_src1 = src1_;
        break;
      case 2:   // src1 == dst
        p_src0 = src0_;
        p_src1 = dst_tst_;
        src1_stride_ = dst_stride_;
        src1_offset_ = dst_offset_;
        break;
      default:
        FAIL();
    }

    Execute(p_src0, p_src1);

    for (int r = 0 ; r < h_ ; ++r) {
      for (int c = 0 ; c < w_ ; ++c) {
        ASSERT_EQ(dst_ref_[dst_offset_ + r * dst_stride_ + c],
                  dst_tst_[dst_offset_ + r * dst_stride_ + c]);
      }
    }
  }

  ACMRandom rng_;

  T dst_ref_[kBufSize];
  T dst_tst_[kBufSize];
  size_t dst_stride_;
  size_t dst_offset_;

  T src0_[kBufSize];
  size_t src0_stride_;
  size_t src0_offset_;

  T src1_[kBufSize];
  size_t src1_stride_;
  size_t src1_offset_;

  uint8_t mask_[kMaxMaskSize];

  int w_;
  int h_;
};

//////////////////////////////////////////////////////////////////////////////
// 8 bit version
//////////////////////////////////////////////////////////////////////////////

typedef void (*F8B)(uint8_t *dst, uint32_t dst_stride,
                    const uint8_t *src0, uint32_t src0_stride,
                    const uint8_t *src1, uint32_t src1_stride,
                    const uint8_t *mask, int h, int w);

class BlendA64Mask1DTest8B : public BlendA64Mask1DTest<F8B, uint8_t> {
 protected:
  void Execute(uint8_t *p_src0, uint8_t *p_src1) {
    ref_func_(dst_ref_ + dst_offset_, dst_stride_,
              p_src0 + src0_offset_, src0_stride_,
              p_src1 + src1_offset_, src1_stride_,
              mask_, h_, w_);

    tst_func_(dst_tst_ + dst_offset_, dst_stride_,
              p_src0 + src0_offset_, src0_stride_,
              p_src1 + src1_offset_, src1_stride_,
              mask_, h_, w_);
  }
};

TEST_P(BlendA64Mask1DTest8B, RandomValues) {
  for (int iter = 0 ; iter < kIterations && !HasFatalFailure(); ++iter) {
    for (int i = 0 ; i < kBufSize ; ++i) {
      dst_ref_[i] = rng_.Rand8();
      dst_tst_[i] = rng_.Rand8();

      src0_[i] = rng_.Rand8();
      src1_[i] = rng_.Rand8();
    }

    for (int i = 0 ; i < kMaxMaskSize ; ++i)
      mask_[i] = rng_(VPX_BLEND_A64_MAX_ALPHA + 1);

    Common();
  }
}

TEST_P(BlendA64Mask1DTest8B, ExtremeValues) {
  for (int iter = 0 ; iter < kIterations && !HasFatalFailure(); ++iter) {
    for (int i = 0 ; i < kBufSize ; ++i) {
      dst_ref_[i] = rng_(2) + 254;
      dst_tst_[i] = rng_(2) + 254;
      src0_[i] = rng_(2) + 254;
      src1_[i] = rng_(2) + 254;
    }

    for (int i = 0 ; i < kMaxMaskSize ; ++i)
      mask_[i] = rng_(2) + VPX_BLEND_A64_MAX_ALPHA - 1;

    Common();
  }
}

static void blend_a64_hmask_ref(
    uint8_t *dst, uint32_t dst_stride,
    const uint8_t *src0, uint32_t src0_stride,
    const uint8_t *src1, uint32_t src1_stride,
    const uint8_t *mask, int h, int w) {
  uint8_t mask2d[BlendA64Mask1DTest8B::kMaxMaskSize]
                [BlendA64Mask1DTest8B::kMaxMaskSize];

  for (int row = 0 ; row < h ; ++row)
    for (int col = 0 ; col < w ; ++col)
      mask2d[row][col] = mask[col];

  vpx_blend_a64_mask_c(dst, dst_stride,
                       src0, src0_stride,
                       src1, src1_stride,
                       &mask2d[0][0], BlendA64Mask1DTest8B::kMaxMaskSize,
                       h, w, 0, 0);
}

static void blend_a64_vmask_ref(
    uint8_t *dst, uint32_t dst_stride,
    const uint8_t *src0, uint32_t src0_stride,
    const uint8_t *src1, uint32_t src1_stride,
    const uint8_t *mask, int h, int w) {
  uint8_t mask2d[BlendA64Mask1DTest8B::kMaxMaskSize]
                [BlendA64Mask1DTest8B::kMaxMaskSize];

  for (int row = 0 ; row < h ; ++row)
    for (int col = 0 ; col < w ; ++col)
      mask2d[row][col] = mask[row];

  vpx_blend_a64_mask_c(dst, dst_stride,
                       src0, src0_stride,
                       src1, src1_stride,
                       &mask2d[0][0], BlendA64Mask1DTest8B::kMaxMaskSize,
                       h, w, 0, 0);
}

INSTANTIATE_TEST_CASE_P(
  C, BlendA64Mask1DTest8B,
  ::testing::Values(
    make_tuple(blend_a64_hmask_ref, vpx_blend_a64_hmask_c),
    make_tuple(blend_a64_vmask_ref, vpx_blend_a64_vmask_c)));

#if HAVE_SSE4_1
INSTANTIATE_TEST_CASE_P(
  SSE4_1, BlendA64Mask1DTest8B,
  ::testing::Values(
    make_tuple(blend_a64_hmask_ref, vpx_blend_a64_hmask_sse4_1),
    make_tuple(blend_a64_vmask_ref, vpx_blend_a64_vmask_sse4_1)));
#endif  // HAVE_SSE4_1

#if CONFIG_VP9_HIGHBITDEPTH
//////////////////////////////////////////////////////////////////////////////
// High bit-depth version
//////////////////////////////////////////////////////////////////////////////

typedef void (*FHBD)(uint8_t *dst, uint32_t dst_stride,
                     const uint8_t *src0, uint32_t src0_stride,
                     const uint8_t *src1, uint32_t src1_stride,
                     const uint8_t *mask, int h, int w, int bd);

class BlendA64Mask1DTestHBD : public BlendA64Mask1DTest<FHBD, uint16_t> {
 protected:
  void Execute(uint16_t *p_src0, uint16_t *p_src1) {
    ref_func_(CONVERT_TO_BYTEPTR(dst_ref_ + dst_offset_), dst_stride_,
              CONVERT_TO_BYTEPTR(p_src0 + src0_offset_), src0_stride_,
              CONVERT_TO_BYTEPTR(p_src1 + src1_offset_), src1_stride_,
              mask_, h_, w_, bit_depth_);

    ASM_REGISTER_STATE_CHECK(
      tst_func_(CONVERT_TO_BYTEPTR(dst_tst_ + dst_offset_), dst_stride_,
                CONVERT_TO_BYTEPTR(p_src0 + src0_offset_), src0_stride_,
                CONVERT_TO_BYTEPTR(p_src1 + src1_offset_), src1_stride_,
                mask_, h_, w_, bit_depth_));
  }

  int bit_depth_;
};

TEST_P(BlendA64Mask1DTestHBD, RandomValues) {
  for (int iter = 0 ; iter < kIterations && !HasFatalFailure(); ++iter) {
    switch (rng_(3)) {
    case 0:
      bit_depth_ = 8;
      break;
    case 1:
      bit_depth_ = 10;
      break;
    default:
      bit_depth_ = 12;
      break;
    }

    const int hi = 1 << bit_depth_;

    for (int i = 0 ; i < kBufSize ; ++i) {
      dst_ref_[i] = rng_(hi);
      dst_tst_[i] = rng_(hi);
      src0_[i] = rng_(hi);
      src1_[i] = rng_(hi);
    }

    for (int i = 0 ; i < kMaxMaskSize ; ++i)
      mask_[i] = rng_(VPX_BLEND_A64_MAX_ALPHA + 1);

    Common();
  }
}

TEST_P(BlendA64Mask1DTestHBD, ExtremeValues) {
  for (int iter = 0 ; iter < 1000 && !HasFatalFailure(); ++iter) {
    switch (rng_(3)) {
    case 0:
      bit_depth_ = 8;
      break;
    case 1:
      bit_depth_ = 10;
      break;
    default:
      bit_depth_ = 12;
      break;
    }

    const int hi = 1 << bit_depth_;
    const int lo = hi - 2;

    for (int i = 0 ; i < kBufSize ; ++i) {
      dst_ref_[i] = rng_(hi - lo) + lo;
      dst_tst_[i] = rng_(hi - lo) + lo;
      src0_[i] = rng_(hi - lo) + lo;
      src1_[i] = rng_(hi - lo) + lo;
    }

    for (int i = 0 ; i < kMaxMaskSize ; ++i)
      mask_[i] = rng_(2) + VPX_BLEND_A64_MAX_ALPHA - 1;

    Common();
  }
}

static void highbd_blend_a64_hmask_ref(
    uint8_t *dst, uint32_t dst_stride,
    const uint8_t *src0, uint32_t src0_stride,
    const uint8_t *src1, uint32_t src1_stride,
    const uint8_t *mask, int h, int w, int bd) {
  uint8_t mask2d[BlendA64Mask1DTestHBD::kMaxMaskSize]
                [BlendA64Mask1DTestHBD::kMaxMaskSize];

  for (int row = 0 ; row < h ; ++row)
    for (int col = 0 ; col < w ; ++col)
      mask2d[row][col] = mask[col];

  vpx_highbd_blend_a64_mask_c(dst, dst_stride,
                              src0, src0_stride,
                              src1, src1_stride,
                              &mask2d[0][0],
                              BlendA64Mask1DTestHBD::kMaxMaskSize,
                              h, w, 0, 0, bd);
}

static void highbd_blend_a64_vmask_ref(
    uint8_t *dst, uint32_t dst_stride,
    const uint8_t *src0, uint32_t src0_stride,
    const uint8_t *src1, uint32_t src1_stride,
    const uint8_t *mask, int h, int w, int bd) {
  uint8_t mask2d[BlendA64Mask1DTestHBD::kMaxMaskSize]
                [BlendA64Mask1DTestHBD::kMaxMaskSize];

  for (int row = 0 ; row < h ; ++row)
    for (int col = 0 ; col < w ; ++col)
      mask2d[row][col] = mask[row];

  vpx_highbd_blend_a64_mask_c(dst, dst_stride,
                              src0, src0_stride,
                              src1, src1_stride,
                              &mask2d[0][0],
                              BlendA64Mask1DTestHBD::kMaxMaskSize,
                              h, w, 0, 0, bd);
}

INSTANTIATE_TEST_CASE_P(
  C, BlendA64Mask1DTestHBD,
  ::testing::Values(
    make_tuple(highbd_blend_a64_hmask_ref, vpx_highbd_blend_a64_hmask_c),
    make_tuple(highbd_blend_a64_vmask_ref, vpx_highbd_blend_a64_vmask_c)));

#if HAVE_SSE4_1
INSTANTIATE_TEST_CASE_P(
  SSE4_1, BlendA64Mask1DTestHBD,
  ::testing::Values(
    make_tuple(highbd_blend_a64_hmask_ref, vpx_highbd_blend_a64_hmask_sse4_1),
    make_tuple(highbd_blend_a64_vmask_ref, vpx_highbd_blend_a64_vmask_sse4_1)));
#endif  // HAVE_SSE4_1

#endif  // CONFIG_VP9_HIGHBITDEPTH
}  // namespace
