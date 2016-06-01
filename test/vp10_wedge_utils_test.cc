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
#include "vpx_ports/mem.h"

#include "./vpx_dsp_rtcd.h"
#include "./vp10_rtcd.h"

#include "vpx_dsp/vpx_dsp_common.h"

#include "vp10/common/enums.h"

#include "test/array_utils.h"
#include "test/assertion_helpers.h"
#include "test/function_equivalence_test.h"
#include "test/randomise.h"
#include "test/register_state_check.h"
#include "test/snapshot.h"

#define WEDGE_WEIGHT_BITS 6
#define MAX_MASK_VALUE  (1 << (WEDGE_WEIGHT_BITS))

using std::tr1::make_tuple;
using libvpx_test::FunctionEquivalenceTest;
using libvpx_test::Snapshot;
using libvpx_test::Randomise;
using libvpx_test::array_utils::arraySet;
using libvpx_test::assertion_helpers::ArraysEq;
using libvpx_test::assertion_helpers::ArraysEqWithin;

namespace {

static const int16_t int13_max = (1<<12) - 1;

//////////////////////////////////////////////////////////////////////////////
// vp10_wedge_sse_from_residuals - functionality
//////////////////////////////////////////////////////////////////////////////

class WedgeUtilsSSEFuncTest : public testing::Test {
 protected:
  Snapshot snapshot;
  Randomise randomise;
};

static void equiv_blend_residuals(int16_t *r,
                                  const int16_t *r0,
                                  const int16_t *r1,
                                  const uint8_t *m,
                                  int N) {
  for (int i = 0 ; i < N ; i++) {
    const int32_t m0 = m[i];
    const int32_t m1 = MAX_MASK_VALUE - m0;
    const int16_t R = m0 * r0[i] + m1 * r1[i];
    // Note that this rounding is designed to match the result
    // you would get when actually blending the 2 predictors and computing
    // the residuals.
    r[i] = ROUND_POWER_OF_TWO(R - 1, WEDGE_WEIGHT_BITS);
  }
}

static uint64_t equiv_sse_from_residuals(const int16_t *r0,
                                         const int16_t *r1,
                                         const uint8_t *m,
                                         int N) {
  uint64_t acc = 0;
  for (int i = 0 ; i < N ; i++) {
    const int32_t m0 = m[i];
    const int32_t m1 = MAX_MASK_VALUE - m0;
    const int16_t R = m0 * r0[i] + m1 * r1[i];
    const int32_t r = ROUND_POWER_OF_TWO(R - 1, WEDGE_WEIGHT_BITS);
    acc += r * r;
  }
  return acc;
}

TEST_F(WedgeUtilsSSEFuncTest, ResidualBlendingEquiv) {
  for (int i = 0 ; i < 1000 && !HasFatalFailure(); i++) {
    uint8_t s[MAX_SB_SQUARE];
    uint8_t p0[MAX_SB_SQUARE];
    uint8_t p1[MAX_SB_SQUARE];
    uint8_t p[MAX_SB_SQUARE];

    int16_t r0[MAX_SB_SQUARE];
    int16_t r1[MAX_SB_SQUARE];
    int16_t r_ref[MAX_SB_SQUARE];
    int16_t r_tst[MAX_SB_SQUARE];
    uint8_t m[MAX_SB_SQUARE];

    randomise(s);
    randomise(m, 0, MAX_MASK_VALUE + 1);

    const int w = 1 << randomise.uniform<uint32_t>(3, MAX_SB_SIZE_LOG2);
    const int h = 1 << randomise.uniform<uint32_t>(3, MAX_SB_SIZE_LOG2);
    const int N = w * h;

    for (int j = 0 ; j < N ; j++) {
      p0[j] = clamp(s[j] + randomise.uniform<int>(-16, 17), 0, UINT8_MAX);
      p1[j] = clamp(s[j] + randomise.uniform<int>(-16, 17), 0, UINT8_MAX);
    }
    vpx_blend_mask6(p, w, p0, w, p1, w, m, w, h, w, 0, 0);

    vpx_subtract_block(h, w, r0, w, s, w, p0, w);
    vpx_subtract_block(h, w, r1, w, s, w, p1, w);

    vpx_subtract_block(h, w, r_ref, w, s, w, p, w);
    equiv_blend_residuals(r_tst, r0, r1, m, N);

    ASSERT_TRUE(ArraysEqWithin(r_ref, r_tst, 0, N));

    uint64_t ref_sse = vpx_sum_squares_i16(r_ref, N);
    uint64_t tst_sse = equiv_sse_from_residuals(r0, r1, m, N);

    ASSERT_EQ(ref_sse, tst_sse);
  }
}

static uint64_t sse_from_residuals(const int16_t *r0,
                                   const int16_t *r1,
                                   const uint8_t *m,
                                   int N) {
  uint64_t acc = 0;
  for (int i = 0 ; i < N ; i++) {
    const int32_t m0 = m[i];
    const int32_t m1 = MAX_MASK_VALUE - m0;
    const int32_t r = m0 * r0[i] + m1 * r1[i];
    acc += r * r;
  }
  return ROUND_POWER_OF_TWO(acc, 2 * WEDGE_WEIGHT_BITS);
}

TEST_F(WedgeUtilsSSEFuncTest, ResidualBlendingMethod) {
  for (int i = 0 ; i < 1000 && !HasFatalFailure(); i++) {
    int16_t r0[MAX_SB_SQUARE];
    int16_t r1[MAX_SB_SQUARE];
    int16_t d[MAX_SB_SQUARE];
    uint8_t m[MAX_SB_SQUARE];

    randomise(r1, 2 * INT8_MIN, 2 * INT8_MAX + 1);
    randomise(d, 2 * INT8_MIN, 2 * INT8_MAX + 1);
    randomise(m, 0, MAX_MASK_VALUE + 1);

    const int N = 64 * randomise.uniform<uint32_t>(1, MAX_SB_SQUARE/64);

    for (int j = 0 ; j < N ; j++)
      r0[j] = r1[j] + d[j];

    uint64_t ref_res, tst_res;

    ref_res = sse_from_residuals(r0, r1, m, N);
    tst_res = vp10_wedge_sse_from_residuals(r1, d, m, N);

    ASSERT_EQ(ref_res, tst_res);
  }
}

//////////////////////////////////////////////////////////////////////////////
// vp10_wedge_sse_from_residuals - optimizations
//////////////////////////////////////////////////////////////////////////////

typedef uint64_t (*FSSE)(const int16_t *r1,
                         const int16_t *d,
                         const uint8_t *m,
                         int N);

class WedgeUtilsSSEOptTest : public FunctionEquivalenceTest<FSSE> {
 protected:
  void Common() {
    const int N = 64 * randomise.uniform<uint32_t>(1, MAX_SB_SQUARE/64);

    snapshot(r1);
    snapshot(d);
    snapshot(m);

    uint64_t ref_res, tst_res;

    ref_res = ref_func_(r1, d, m, N);
    ASM_REGISTER_STATE_CHECK(tst_res = tst_func_(r1, d, m, N));

    ASSERT_EQ(ref_res, tst_res);

    ASSERT_TRUE(ArraysEq(snapshot.get(r1), r1));
    ASSERT_TRUE(ArraysEq(snapshot.get(d), d));
    ASSERT_TRUE(ArraysEq(snapshot.get(m), m));
  }

  Snapshot snapshot;
  Randomise randomise;

  DECLARE_ALIGNED(16, int16_t, r1[MAX_SB_SQUARE]);
  DECLARE_ALIGNED(16, int16_t, d[MAX_SB_SQUARE]);
  DECLARE_ALIGNED(16, uint8_t, m[MAX_SB_SQUARE]);
};

TEST_P(WedgeUtilsSSEOptTest, RandomValues) {
  for (int i = 0 ; i < 10000 && !HasFatalFailure(); i++) {
    randomise(r1, -int13_max, int13_max + 1);
    randomise(d, -int13_max, int13_max + 1);
    randomise(m, 0, 65);

    Common();
  }
}

TEST_P(WedgeUtilsSSEOptTest, ExtremeValues) {
  for (int i = 0 ; i < 10000 && !HasFatalFailure(); i++) {
    if (randomise.uniform<bool>())
      arraySet(r1, int13_max);
    else
      arraySet(r1, -int13_max);

    if (randomise.uniform<bool>())
      arraySet(d, int13_max);
    else
      arraySet(d, -int13_max);

    arraySet(m, MAX_MASK_VALUE);

    Common();
  }
}

#if HAVE_SSE2
INSTANTIATE_TEST_CASE_P(
    SSE2, WedgeUtilsSSEOptTest,
    ::testing::Values(
        make_tuple(&vp10_wedge_sse_from_residuals_c,
                   &vp10_wedge_sse_from_residuals_sse2)
    )
);
#endif  // HAVE_SSE2

//////////////////////////////////////////////////////////////////////////////
// vp10_wedge_sign_from_residuals
//////////////////////////////////////////////////////////////////////////////

typedef int (*FSign)(const int16_t *ds,
                     const uint8_t *m,
                     int N,
                     int64_t limit);

class WedgeUtilsSignOptTest : public FunctionEquivalenceTest<FSign> {
 protected:
  static const int maxSize = 8196;  // Size limited by SIMD implementation.

  void Common() {
    const int maxN = VPXMIN(maxSize, MAX_SB_SQUARE);
    const int N = 64 * randomise.uniform<uint32_t>(1, maxN/64);

    int64_t limit;
    limit = (int64_t)vpx_sum_squares_i16(r0, N);
    limit -= (int64_t)vpx_sum_squares_i16(r1, N);
    limit *= (1 << WEDGE_WEIGHT_BITS) / 2;

    for (int i = 0 ; i < N ; i++)
      ds[i] = clamp(r0[i]*r0[i] - r1[i]*r1[i], INT16_MIN, INT16_MAX);

    snapshot(r0);
    snapshot(r1);
    snapshot(ds);
    snapshot(m);

    int ref_res, tst_res;

    ref_res = ref_func_(ds, m, N, limit);
    ASM_REGISTER_STATE_CHECK(tst_res = tst_func_(ds, m, N, limit));

    ASSERT_EQ(ref_res, tst_res);

    ASSERT_TRUE(ArraysEq(snapshot.get(r0), r0));
    ASSERT_TRUE(ArraysEq(snapshot.get(r1), r1));
    ASSERT_TRUE(ArraysEq(snapshot.get(ds), ds));
    ASSERT_TRUE(ArraysEq(snapshot.get(m), m));
  }

  Snapshot snapshot;
  Randomise randomise;

  DECLARE_ALIGNED(16, int16_t, r0[MAX_SB_SQUARE]);
  DECLARE_ALIGNED(16, int16_t, r1[MAX_SB_SQUARE]);
  DECLARE_ALIGNED(16, int16_t, ds[MAX_SB_SQUARE]);
  DECLARE_ALIGNED(16, uint8_t, m[MAX_SB_SQUARE]);
};

TEST_P(WedgeUtilsSignOptTest, RandomValues) {
  for (int i = 0 ; i < 10000 && !HasFatalFailure(); i++) {
    randomise(r0, -int13_max, int13_max+1);
    randomise(r1, -int13_max, int13_max+1);
    randomise(m, 0, MAX_MASK_VALUE + 1);

    Common();
  }
}

TEST_P(WedgeUtilsSignOptTest, ExtremeValues) {
  for (int i = 0 ; i < 10000 && !HasFatalFailure(); i++) {
    switch (randomise.uniform<int>(4)) {
    case 0:
      arraySet(r0, 0);
      arraySet(r1, int13_max);
      break;
    case 1:
      arraySet(r0, int13_max);
      arraySet(r1, 0);
      break;
    case 2:
      arraySet(r0, 0);
      arraySet(r1, -int13_max);
      break;
    default:
      arraySet(r0, -int13_max);
      arraySet(r1, 0);
      break;
    }

    arraySet(m, MAX_MASK_VALUE);

    Common();
  }
}

#if HAVE_SSE2
INSTANTIATE_TEST_CASE_P(
    SSE2, WedgeUtilsSignOptTest,
    ::testing::Values(
        make_tuple(&vp10_wedge_sign_from_residuals_c,
                   &vp10_wedge_sign_from_residuals_sse2)
    )
);
#endif  // HAVE_SSE2

//////////////////////////////////////////////////////////////////////////////
// vp10_wedge_compute_delta_squares
//////////////////////////////////////////////////////////////////////////////

typedef void (*FDS)(int16_t *d,
                    const int16_t *a,
                    const int16_t *b,
                    int N);

class WedgeUtilsDeltaSquaresOptTest : public FunctionEquivalenceTest<FDS> {
 protected:
  void Common() {
    const int N = 64 * randomise.uniform<uint32_t>(1, MAX_SB_SQUARE/64);

    randomise(d_ref);
    randomise(d_tst);

    snapshot(a);
    snapshot(b);

    ref_func_(d_ref, a, b, N);
    ASM_REGISTER_STATE_CHECK(tst_func_(d_tst, a, b, N));

    ASSERT_TRUE(ArraysEqWithin(d_ref, d_tst, 0, N));

    ASSERT_TRUE(ArraysEq(snapshot.get(a), a));
    ASSERT_TRUE(ArraysEq(snapshot.get(b), b));
  }

  Snapshot snapshot;
  Randomise randomise;

  DECLARE_ALIGNED(16, int16_t, a[MAX_SB_SQUARE]);
  DECLARE_ALIGNED(16, int16_t, b[MAX_SB_SQUARE]);
  DECLARE_ALIGNED(16, int16_t, d_ref[MAX_SB_SQUARE]);
  DECLARE_ALIGNED(16, int16_t, d_tst[MAX_SB_SQUARE]);
};

TEST_P(WedgeUtilsDeltaSquaresOptTest, RandomValues) {
  for (int i = 0 ; i < 10000 && !HasFatalFailure(); i++) {
    randomise(a);
    randomise(b, -INT16_MAX, INT16_MAX + 1);

    Common();
  }
}

#if HAVE_SSE2
INSTANTIATE_TEST_CASE_P(
    SSE2, WedgeUtilsDeltaSquaresOptTest,
    ::testing::Values(
        make_tuple(&vp10_wedge_compute_delta_squares_c,
                   &vp10_wedge_compute_delta_squares_sse2)
    )
);
#endif  // HAVE_SSE2


}  // namespace
