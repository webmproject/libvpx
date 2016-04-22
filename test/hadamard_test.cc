/*
 *  Copyright (c) 2016 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include <algorithm>

#include "third_party/googletest/src/include/gtest/gtest.h"

#include "./vpx_dsp_rtcd.h"

#include "test/acm_random.h"
#include "test/register_state_check.h"

namespace {

using ::libvpx_test::ACMRandom;

typedef void (*Hadamard8x8Func)(const int16_t *a, int a_stride,
                                int16_t *b);

class HadamardTest : public ::testing::TestWithParam<Hadamard8x8Func> {
 public:
  virtual void SetUp() {
    h_func_ = GetParam();
    rnd_.Reset(ACMRandom::DeterministicSeed());
  }

 protected:
  Hadamard8x8Func h_func_;
  ACMRandom rnd_;
};

void hadamard_loop(const int16_t *a, int a_stride, int16_t *out) {
  int16_t b[8];
  for (int i = 0; i < 8; i += 2) {
    b[i + 0] = a[i * a_stride] + a[(i + 1) * a_stride];
    b[i + 1] = a[i * a_stride] - a[(i + 1) * a_stride];
  }
  int16_t c[8];
  for (int i = 0; i < 8; i += 4) {
    c[i + 0] = b[i + 0] + b[i + 2];
    c[i + 1] = b[i + 1] + b[i + 3];
    c[i + 2] = b[i + 0] - b[i + 2];
    c[i + 3] = b[i + 1] - b[i + 3];
  }
  out[0] = c[0] + c[4];
  out[7] = c[1] + c[5];
  out[3] = c[2] + c[6];
  out[4] = c[3] + c[7];
  out[2] = c[0] - c[4];
  out[6] = c[1] - c[5];
  out[1] = c[2] - c[6];
  out[5] = c[3] - c[7];
}

void reference_hadamard(const int16_t *a, int a_stride, int16_t *b) {
  int16_t buf[64];
  for (int i = 0; i < 8; i++) {
    hadamard_loop(a + i, a_stride, buf + i * 8);
  }

  for (int i = 0; i < 8; i++) {
    hadamard_loop(buf + i, 8, b + i * 8);
  }
}

TEST_P(HadamardTest, CompareReferenceRandom) {
  DECLARE_ALIGNED(16, int16_t, a[64]);
  DECLARE_ALIGNED(16, int16_t, b[64]);
  int16_t b_ref[64];
  for (int i = 0; i < 64; i++) {
    a[i] = rnd_.Rand9Signed();
  }
  memset(b, 0, sizeof(b));
  memset(b_ref, 0, sizeof(b_ref));

  reference_hadamard(a, 8, b_ref);
  ASM_REGISTER_STATE_CHECK(h_func_(a, 8, b));

  // The order of the output is not important. Sort before checking.
  std::sort(b, b + 64);
  std::sort(b_ref, b_ref + 64);
  EXPECT_EQ(0, memcmp(b, b_ref, sizeof(b)));
}

TEST_P(HadamardTest, VaryStride) {
  DECLARE_ALIGNED(16, int16_t, a[64 * 8]);
  DECLARE_ALIGNED(16, int16_t, b[64]);
  int16_t b_ref[64];
  for (int i = 0; i < 64 * 8; i++) {
    a[i] = rnd_.Rand9Signed();
  }

  for (int i = 8; i < 64; i += 8) {
    memset(b, 0, sizeof(b));
    memset(b_ref, 0, sizeof(b_ref));

    reference_hadamard(a, i, b_ref);
    ASM_REGISTER_STATE_CHECK(h_func_(a, i, b));

    // The order of the output is not important. Sort before checking.
    std::sort(b, b + 64);
    std::sort(b_ref, b_ref + 64);
    EXPECT_EQ(0, memcmp(b, b_ref, sizeof(b)));
  }
}

INSTANTIATE_TEST_CASE_P(C, HadamardTest,
                        ::testing::Values(&vpx_hadamard_8x8_c));

#if HAVE_SSE2
INSTANTIATE_TEST_CASE_P(SSE2, HadamardTest,
                        ::testing::Values(&vpx_hadamard_8x8_sse2));
#endif  // HAVE_SSE2

#if HAVE_SSSE3 && CONFIG_USE_X86INC && ARCH_X86_64
INSTANTIATE_TEST_CASE_P(SSSE3, HadamardTest,
                        ::testing::Values(&vpx_hadamard_8x8_ssse3));
#endif  // HAVE_SSSE3 && CONFIG_USE_X86INC && ARCH_X86_64
}  // namespace
