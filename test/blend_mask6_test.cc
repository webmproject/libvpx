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
#include "test/randomise.h"
#include "test/snapshot.h"

#include "./vpx_config.h"
#include "./vpx_dsp_rtcd.h"
#include "vpx/vpx_integer.h"

#include "./vp10_rtcd.h"

#include "test/assertion_helpers.h"
#include "vp10/common/enums.h"

using libvpx_test::assertion_helpers::BuffersEqWithin;
using libvpx_test::assertion_helpers::BuffersEqOutside;
using libvpx_test::assertion_helpers::ArraysEq;
using libvpx_test::FunctionEquivalenceTest;
using libvpx_test::Snapshot;
using libvpx_test::Randomise;
using std::tr1::make_tuple;

namespace {

template<typename F, typename T>
class BlendMask6Test : public FunctionEquivalenceTest<F> {
 protected:
  virtual ~BlendMask6Test() {}

  virtual void Execute(T *p_src0, T *p_src1) = 0;

  void Common() {
    w = 1 << randomise.uniform<int>(2, MAX_SB_SIZE_LOG2 + 1);
    h = 1 << randomise.uniform<int>(2, MAX_SB_SIZE_LOG2 + 1);

    randomise(subx);
    randomise(suby);

    randomise(dst_offset, 0, 32);
    randomise(dst_stride, w, MAX_SB_SIZE * 5 + 1);

    randomise(src0_offset, 0, 32);
    randomise(src0_stride, w, MAX_SB_SIZE * 5 + 1);

    randomise(src1_offset, 0, 32);
    randomise(src1_stride, w, MAX_SB_SIZE * 5 + 1);

    randomise(mask_stride, w * (subx ? 2: 1), 2 * MAX_SB_SIZE + 1);

    T *p_src0;
    T *p_src1;

    switch (randomise.uniform<int>(3)) {
      case 0:   // Separate sources
        p_src0 = &src0[0][0];
        p_src1 = &src1[0][0];
        break;
      case 1:   // src0 == dst
        p_src0 = &dst_tst[0][0];
        src0_stride = dst_stride;
        src0_offset = dst_offset;
        p_src1 = &src1[0][0];
        break;
      case 2:   // src1 == dst
        p_src0 = &src0[0][0];
        p_src1 = &dst_tst[0][0];
        src1_stride = dst_stride;
        src1_offset = dst_offset;
        break;
      default:
        FAIL();
    }

    //////////////////////////////////////////////////////////////////////////
    // Prepare
    //////////////////////////////////////////////////////////////////////////

    snapshot(dst_ref);
    snapshot(dst_tst);

    snapshot(src0);
    snapshot(src1);

    snapshot(mask);

    //////////////////////////////////////////////////////////////////////////
    // Execute
    //////////////////////////////////////////////////////////////////////////

    Execute(p_src0, p_src1);

    //////////////////////////////////////////////////////////////////////////
    // Check
    //////////////////////////////////////////////////////////////////////////

    ASSERT_TRUE(BuffersEqWithin(dst_ref, dst_tst,
                                dst_stride, dst_stride,
                                dst_offset, dst_offset,
                                h, w));

    ASSERT_TRUE(ArraysEq(snapshot.get(src0), src0));
    ASSERT_TRUE(ArraysEq(snapshot.get(src1), src1));
    ASSERT_TRUE(ArraysEq(snapshot.get(mask), mask));

    ASSERT_TRUE(BuffersEqOutside(snapshot.get(dst_ref), dst_ref,
                                 dst_stride,
                                 dst_offset,
                                 h, w));

    ASSERT_TRUE(BuffersEqOutside(snapshot.get(dst_tst), dst_tst,
                                 dst_stride,
                                 dst_offset,
                                 h, w));
  }

  Snapshot snapshot;
  Randomise randomise;

  T dst_ref[MAX_SB_SIZE][MAX_SB_SIZE * 5];
  T dst_tst[MAX_SB_SIZE][MAX_SB_SIZE * 5];
  size_t dst_stride;
  size_t dst_offset;

  T src0[MAX_SB_SIZE][MAX_SB_SIZE * 5];
  size_t src0_stride;
  size_t src0_offset;

  T src1[MAX_SB_SIZE][MAX_SB_SIZE * 5];
  size_t src1_stride;
  size_t src1_offset;

  uint8_t mask[2 * MAX_SB_SIZE][2 * MAX_SB_SIZE];
  size_t mask_stride;

  int w;
  int h;

  bool suby;
  bool subx;
};

//////////////////////////////////////////////////////////////////////////////
// 8 bit version
//////////////////////////////////////////////////////////////////////////////

typedef void (*F8B)(uint8_t *dst, uint32_t dst_stride,
                      uint8_t *src0, uint32_t src0_stride,
                      uint8_t *src1, uint32_t src1_stride,
                      const uint8_t *mask, uint32_t mask_stride,
                      int h, int w, int suby, int subx);

class BlendMask6Test8B : public BlendMask6Test<F8B, uint8_t> {
 protected:
  void Execute(uint8_t *p_src0, uint8_t *p_src1) {
    ref_func_(&dst_ref[0][dst_offset], dst_stride,
              p_src0 + src0_offset, src0_stride,
              p_src1 + src1_offset, src1_stride,
              &mask[0][0], sizeof(mask[0]),
              h, w, suby, subx);

    ASM_REGISTER_STATE_CHECK(
      tst_func_(&dst_tst[0][dst_offset], dst_stride,
                p_src0 + src0_offset, src0_stride,
                p_src1 + src1_offset, src1_stride,
                &mask[0][0], sizeof(mask[0]),
                h, w, suby, subx));
  }
};

TEST_P(BlendMask6Test8B, RandomValues) {
  for (int i = 0 ; i < 10000 && !HasFatalFailure(); i++) {
    //////////////////////////////////////////////////////////////////////////
    // Randomise
    //////////////////////////////////////////////////////////////////////////

    randomise(dst_ref);
    randomise(dst_tst);

    randomise(src0);
    randomise(src1);

    randomise(mask, 65);

    Common();
  }
}

TEST_P(BlendMask6Test8B, ExtremeValues) {
  for (int i = 0 ; i < 1000 && !HasFatalFailure(); i++) {
    //////////////////////////////////////////////////////////////////////////
    // Randomise
    //////////////////////////////////////////////////////////////////////////

    randomise(dst_ref, 254, 256);
    randomise(dst_tst, 254, 256);

    randomise(src0, 254, 256);
    randomise(src1, 254, 256);

    randomise(mask, 63, 65);

    Common();
  }
}

#if HAVE_SSE4_1
INSTANTIATE_TEST_CASE_P(
  SSE4_1_C_COMPARE, BlendMask6Test8B,
  ::testing::Values(make_tuple(&vpx_blend_mask6_c, &vpx_blend_mask6_sse4_1)));
#endif  // HAVE_SSE4_1

#if CONFIG_VP9_HIGHBITDEPTH
//////////////////////////////////////////////////////////////////////////////
// High bit-depth version
//////////////////////////////////////////////////////////////////////////////

typedef void (*FHBD)(uint8_t *dst, uint32_t dst_stride,
                       uint8_t *src0, uint32_t src0_stride,
                       uint8_t *src1, uint32_t src1_stride,
                       const uint8_t *mask, uint32_t mask_stride,
                       int h, int w, int suby, int subx, int bd);

class BlendMask6TestHBD : public BlendMask6Test<FHBD, uint16_t> {
 protected:
  void Execute(uint16_t *p_src0, uint16_t *p_src1) {
    ref_func_(CONVERT_TO_BYTEPTR(&dst_ref[0][dst_offset]), dst_stride,
              CONVERT_TO_BYTEPTR(p_src0 + src0_offset), src0_stride,
              CONVERT_TO_BYTEPTR(p_src1 + src1_offset), src1_stride,
              &mask[0][0], sizeof(mask[0]),
              h, w, suby, subx, bit_depth);

    ASM_REGISTER_STATE_CHECK(
      tst_func_(CONVERT_TO_BYTEPTR(&dst_tst[0][dst_offset]), dst_stride,
                CONVERT_TO_BYTEPTR(p_src0 + src0_offset), src0_stride,
                CONVERT_TO_BYTEPTR(p_src1 + src1_offset), src1_stride,
                &mask[0][0], sizeof(mask[0]),
                h, w, suby, subx, bit_depth));
  }

  int bit_depth;
};

TEST_P(BlendMask6TestHBD, RandomValues) {
  for (int i = 0 ; i < 10000 && !HasFatalFailure(); i++) {
    //////////////////////////////////////////////////////////////////////////
    // Randomise
    //////////////////////////////////////////////////////////////////////////

    bit_depth = randomise.choice(8, 10, 12);

    const int hi = 1 << bit_depth;

    randomise(dst_ref, hi);
    randomise(dst_tst, hi);

    randomise(src0, hi);
    randomise(src1, hi);

    randomise(mask, 65);

    Common();
  }
}

TEST_P(BlendMask6TestHBD, ExtremeValues) {
  for (int i = 0 ; i < 1000 && !HasFatalFailure(); i++) {
    //////////////////////////////////////////////////////////////////////////
    // Randomise
    //////////////////////////////////////////////////////////////////////////

    bit_depth = randomise.choice(8, 10, 12);

    const int hi = 1 << bit_depth;
    const int lo = hi - 2;

    randomise(dst_ref, lo, hi);
    randomise(dst_tst, lo, hi);

    randomise(src0, lo, hi);
    randomise(src1, lo, hi);

    randomise(mask, 63, 65);

    Common();
  }
}

#if HAVE_SSE4_1
INSTANTIATE_TEST_CASE_P(
  SSE4_1_C_COMPARE, BlendMask6TestHBD,
  ::testing::Values(make_tuple(&vpx_highbd_blend_mask6_c,
                               &vpx_highbd_blend_mask6_sse4_1)));
#endif  // HAVE_SSE4_1
#endif  // CONFIG_VP9_HIGHBITDEPTH
}  // namespace
