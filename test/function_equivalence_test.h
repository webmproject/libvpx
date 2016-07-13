/*
 *  Copyright (c) 2016 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#ifndef TEST_FUNCTION_EQUIVALENCE_TEST_H_
#define TEST_FUNCTION_EQUIVALENCE_TEST_H_

#include "third_party/googletest/src/include/gtest/gtest.h"
#include "test/acm_random.h"
#include "test/clear_system_state.h"
#include "test/util.h"

using libvpx_test::ACMRandom;

namespace libvpx_test {
// Base class for tests that compare 2 implementations of the same function
// for equivalence. The template parameter should be pointer to a function
// that is being tested.
//
// The test takes a 3-parameters encapsulating struct 'FuncParam', containing:
//   - Pointer to reference function
//   - Pointer to tested function
//   - Integer bit depth (default to 0).
//
// These values are then accessible in the tests as member of params_:
// params_.ref_func, params_.tst_func, and params_.bit_depth.
//

template <typename T>
struct FuncParam {
  FuncParam(T ref = NULL, T tst = NULL, int bit_depth = 0)
      : ref_func(ref), tst_func(tst), bit_depth(bit_depth) {}
  T ref_func;
  T tst_func;
  int bit_depth;
};

template <typename T>
class FunctionEquivalenceTest : public ::testing::TestWithParam<FuncParam<T> > {
 public:
  FunctionEquivalenceTest() : rng_(ACMRandom::DeterministicSeed()) {}

  virtual ~FunctionEquivalenceTest() {}

  virtual void SetUp() {
    params_ = this->GetParam();
  }

  virtual void TearDown() {
    libvpx_test::ClearSystemState();
  }

 protected:
  ACMRandom rng_;
  FuncParam<T> params_;
};

}   // namespace libvpx_test
#endif  // TEST_FUNCTION_EQUIVALENCE_TEST_H_
