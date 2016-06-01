/*
 *  Copyright (c) 2016 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */


#ifndef TEST_ARRAY_UTILS_H_
#define TEST_ARRAY_UTILS_H_

#include "third_party/googletest/src/include/gtest/gtest.h"

namespace libvpx_test {
namespace array_utils {

template<typename T, size_t n, typename V>
void arraySet(T (&arr)[n], const V &v) {
  for (size_t i = 0; i < n ; i++) {
    arr[i] = v;
  }
}

template<typename T, size_t n, size_t m, typename V>
void arraySet(T (&arr)[n][m], const V &v) {
  for (size_t i = 0; i < n ; i++) {
    for (size_t j = 0; j < m ; j++) {
      arr[i][j] = v;
    }
  }
}

}   // namespace array_utils
}   // namespace libvpx_test

#endif  // TEST_ARRAY_UTILS_H_
