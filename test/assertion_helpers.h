/*
 *  Copyright (c) 2016 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */


#ifndef TEST_ASSERTION_HELPERS_H_
#define TEST_ASSERTION_HELPERS_H_

#include "third_party/googletest/src/include/gtest/gtest.h"

namespace libvpx_test {
namespace assertion_helpers {

// Arrays (1D) are element-wise equal
template<typename E, size_t n>
::testing::AssertionResult ArraysEq(const E (&a)[n],
                                    const E (&b)[n]) {
  for (size_t i = 0; i < n; i++) {
    const E &va = a[i];
    const E &vb = b[i];
    if (va != vb) {
      return ::testing::AssertionFailure()
        << "Arrays do not equal at index "
        << "[" << i << "]"
        << " values are: " << va << " vs " << vb;
    }
  }

  return ::testing::AssertionSuccess();
}

// Arrays (1D) are element-wise equal
// within the index interval [lo, hi)
template<typename E, size_t n>
::testing::AssertionResult ArraysEqWithin(const E (&a)[n],
                                          const E (&b)[n],
                                          const size_t lo,
                                          const size_t hi) {
  assert(hi > lo);
  assert(hi <= n);

  for (size_t i = lo; i < hi; i++) {
    const E &va = a[i];
    const E &vb = b[i];
    if (va != vb) {
      return ::testing::AssertionFailure()
        << "Arrays do not equal at index "
        << "[" << i << "]"
        << " values are: " << va << " vs " << vb;
    }
  }

  return ::testing::AssertionSuccess();
}

// Arrays (1D) are element-wise equal
// outside the index interval [lo, hi)
template<typename E, size_t n>
::testing::AssertionResult ArraysEqOutside(const E (&a)[n],
                                           const E (&b)[n],
                                           const size_t lo,
                                           const size_t hi) {
  assert(hi > lo);
  assert(hi <= n);

  for (size_t i = 0; i < n; i++) {
    if (lo <= i && i < hi)
      continue;

    const E &va = a[i];
    const E &vb = b[i];
    if (va != vb) {
      return ::testing::AssertionFailure()
        << "Arrays do not equal at index "
        << "[" << i << "]"
        << " values are: " << va << " vs " << vb;
    }
  }

  return ::testing::AssertionSuccess();
}

// Arrays (2D) are element-wise equal
template<typename E, size_t n, size_t m>
::testing::AssertionResult ArraysEq(const E (&a)[n][m],
                                    const E (&b)[n][m]) {
  for (size_t i = 0; i < n; i++) {
    for (size_t j = 0; j < m; j++) {
      const E &va = a[i][j];
      const E &vb = b[i][j];
      if (va != vb) {
        return ::testing::AssertionFailure()
          << "Arrays do not equal at index "
          << "[" << i << "][" << j << "]"
          << " values are: " << va << " vs " << vb;
      }
    }
  }

  return ::testing::AssertionSuccess();
}

// Arrays (2D) are element-wise equal
// within the index interval [lo0, hi0) x [lo1, hi1) (Cartesian product)
template<typename E, size_t n, size_t m>
::testing::AssertionResult ArraysEqWithin(const E (&a)[n][m],
                                          const E (&b)[n][m],
                                          const size_t lo0,
                                          const size_t hi0,
                                          const size_t lo1,
                                          const size_t hi1) {
  assert(hi0 > lo0);
  assert(hi0 <= n);
  assert(hi1 > lo1);
  assert(hi1 <= m);

  for (size_t i = lo0; i < hi0; i++) {
    for (size_t j = lo1; j < hi1; j++) {
      const E &va = a[i][j];
      const E &vb = b[i][j];
      if (va != vb) {
        return ::testing::AssertionFailure()
          << "Arrays do not equal at index "
          << "[" << i << "][" << j << "]"
          << " values are: " << va << " vs " << vb;
      }
    }
  }

  return ::testing::AssertionSuccess();
}

// Arrays (2D) are element-wise equal
// outside the index interval [lo0, hi0) x [lo1, hi1) (Cartesian product)
template<typename E, size_t n, size_t m>
::testing::AssertionResult ArraysEqOutside(const E (&a)[n][m],
                                           const E (&b)[n][m],
                                           const size_t lo0,
                                           const size_t hi0,
                                           const size_t lo1,
                                           const size_t hi1) {
  assert(hi0 > lo0);
  assert(hi0 <= n);
  assert(hi1 > lo1);
  assert(hi1 <= m);

  for (size_t i = 0; i < n; i++) {
    if (lo0 <= i && i < hi0)
      continue;

    for (size_t j = 0; j < m; j++) {
      if (lo1 <= j && j < hi1)
        continue;

      const E &va = a[i][j];
      const E &vb = b[i][j];
      if (va != vb) {
        return ::testing::AssertionFailure()
          << "Arrays do not equal at index "
          << "[" << i << "][" << j << "]"
          << " values are: " << va << " vs " << vb;
      }
    }
  }

  return ::testing::AssertionSuccess();
}

// Non contiguous 2D array buffers are element-wise equal
// at corresponding linear indices specified by rows/cols/stride/offset
template<typename E, size_t n, size_t m>
::testing::AssertionResult BuffersEqWithin(const E (&a)[n][m],
                                           const E (&b)[n][m],
                                           const size_t stridea,
                                           const size_t strideb,
                                           const size_t offseta,
                                           const size_t offsetb,
                                           const size_t rows,
                                           const size_t cols) {
  assert(rows <= n);
  assert(cols <= m);
  assert(stridea <= m);
  assert(strideb <= m);
  assert(cols <= stridea);
  assert(cols <= strideb);
  assert(offseta < n * m);
  assert(offsetb < n * m);
  assert(offseta + (rows - 1) * stridea + (cols - 1) < n * m);
  assert(offsetb + (rows - 1) * strideb + (cols - 1) < n * m);

  const E *pa = &a[0][0] + offseta;
  const E *pb = &b[0][0] + offsetb;

  for (size_t r = 0 ; r < rows ; r++) {
    for (size_t c = 0 ; c < cols ; c++) {
      const E &va = pa[c];
      const E &vb = pb[c];
      if (va != vb) {
        return ::testing::AssertionFailure()
          << "Arrays do not equal at linear index "
          << "[" << pa - &a[0][0]  << "] vs [" << pb - &b[0][0]  << "]"
          << " row=" << r << " col=" << c
          << " values are: " << va << " vs " << vb;
      }
    }
    pa += stridea;
    pb += strideb;
  }

  return ::testing::AssertionSuccess();
}

// Non contiguous 2D array buffers are element-wise equal
// except at corresponding linear indices specified by
// rows/cols/stride/offset.
template<typename E, size_t n, size_t m>
::testing::AssertionResult BuffersEqOutside(const E (&a)[n][m],
                                            const E (&b)[n][m],
                                            const size_t stride,
                                            const size_t offset,
                                            const size_t rows,
                                            const size_t cols ) {
  assert(rows <= n);
  assert(cols <= m);
  assert(stride <= m);
  assert(cols <= stride);
  assert(offset < n * m);
  assert(offset + (rows - 1) * stride + (cols - 1) < n * m);

  const E *const pa = &a[0][0];
  const E *const pb = &b[0][0];

  size_t idx = 0;
  size_t r = 0;
  size_t end = offset;  // beginning of first row

  while (idx < n * m) {
    while (idx < end) {   // until beginning of row or end of buffer
      const E &va = pa[idx];
      const E &vb = pb[idx];
      if (va != vb) {
        return ::testing::AssertionFailure()
          << "Arrays do not equal at index "
          << "[" << idx / m << "][" << idx % m << "]"
          << " values are: " << va << " vs " << vb;
      }

      idx++;
    }

    // Move past row end
    idx += cols;

    if (++r < rows) {
      // Move to next row
      end += stride;
    } else {
      // Move to end of buffer
      end = n * m;
    }
  }

  // Sanity check
  assert(idx == n * m + cols);

  return ::testing::AssertionSuccess();
}

}   // namespace assertion_helpers
}   // namespace libvpx_test

#endif  // TEST_ASSERTION_HELPERS_H_
