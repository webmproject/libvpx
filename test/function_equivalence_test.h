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
#include "test/clear_system_state.h"
#include "test/util.h"

namespace libvpx_test {
template <typename T>
class FunctionEquivalenceTest :
  public ::testing::TestWithParam< std::tr1::tuple< T, T > > {
 public:
  virtual ~FunctionEquivalenceTest() {}

  virtual void SetUp() {
    ref_func_ = std::tr1::get<0>(this->GetParam());
    tst_func_ = std::tr1::get<1>(this->GetParam());
  }

  virtual void TearDown() {
    libvpx_test::ClearSystemState();
  }

 protected:
  T ref_func_;
  T tst_func_;
};

}   // namespace libvpx_test
#endif  // TEST_FUNCTION_EQUIVALENCE_TEST_H_
