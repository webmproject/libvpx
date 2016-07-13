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
#include "test/clear_system_state.h"
#include "test/register_state_check.h"
#include "third_party/googletest/src/include/gtest/gtest.h"
#include "./vpx_dsp_rtcd.h"
#include "vpx/vpx_integer.h"
#include "vpx_dsp/postproc.h"
#include "vpx_mem/vpx_mem.h"

namespace {

// TODO(jimbankoski): make width and height integers not unsigned.
typedef void (*AddNoiseFunc)(unsigned char *start, char *noise,
                             char blackclamp[16], char whiteclamp[16],
                             char bothclamp[16], unsigned int width,
                             unsigned int height, int pitch);

class AddNoiseTest
    : public ::testing::TestWithParam<AddNoiseFunc> {
 public:
  virtual void TearDown() {
    libvpx_test::ClearSystemState();
  }
  virtual ~AddNoiseTest() {}
};

double stddev6(char a, char b, char c, char d, char e, char f) {
  const double n = (a + b + c + d + e + f) / 6.0;
  const double v = ((a - n) * (a - n) + (b - n) * (b - n) + (c - n) * (c - n) +
                    (d - n) * (d - n) + (e - n) * (e - n) + (f - n) * (f - n)) /
                   6.0;
  return sqrt(v);
}

TEST_P(AddNoiseTest, CheckNoiseAdded) {
  DECLARE_ALIGNED(16, char, blackclamp[16]);
  DECLARE_ALIGNED(16, char, whiteclamp[16]);
  DECLARE_ALIGNED(16, char, bothclamp[16]);
  const int width  = 64;
  const int height = 64;
  const int image_size = width * height;
  char noise[3072];
  const int clamp = vpx_setup_noise(sizeof(noise), 4.4, noise);

  for (int i = 0; i < 16; i++) {
    blackclamp[i] = clamp;
    whiteclamp[i] = clamp;
    bothclamp[i] = 2 * clamp;
  }

  uint8_t *const s = reinterpret_cast<uint8_t *>(vpx_calloc(image_size, 1));
  memset(s, 99, image_size);

  ASM_REGISTER_STATE_CHECK(GetParam()(s, noise, blackclamp, whiteclamp,
                                      bothclamp, width, height, width));

  // Check to make sure we don't end up having either the same or no added
  // noise either vertically or horizontally.
  for (int i = 0; i < image_size - 6 * width - 6; ++i) {
    const double hd = stddev6(s[i] - 99, s[i + 1] - 99, s[i + 2] - 99,
                              s[i + 3] - 99, s[i + 4] - 99, s[i + 5] - 99);
    const double vd = stddev6(s[i] - 99, s[i + width] - 99,
                              s[i + 2 * width] - 99, s[i + 3 * width] - 99,
                              s[i + 4 * width] - 99, s[i + 5 * width] - 99);

    EXPECT_NE(hd, 0);
    EXPECT_NE(vd, 0);
  }

  // Initialize pixels in the image to 255 and check for roll over.
  memset(s, 255, image_size);

  ASM_REGISTER_STATE_CHECK(GetParam()(s, noise, blackclamp, whiteclamp,
                                      bothclamp, width, height, width));

  // Check to make sure don't roll over.
  for (int i = 0; i < image_size; ++i) {
    EXPECT_GT((int)s[i], clamp) << "i = " << i;
  }

  // Initialize pixels in the image to 0 and check for roll under.
  memset(s, 0, image_size);

  ASM_REGISTER_STATE_CHECK(GetParam()(s, noise, blackclamp, whiteclamp,
                                      bothclamp, width, height, width));

  // Check to make sure don't roll under.
  for (int i = 0; i < image_size; ++i) {
    EXPECT_LT((int)s[i], 255 - clamp) << "i = " << i;
  }

  vpx_free(s);
}

TEST_P(AddNoiseTest, CheckCvsAssembly) {
  DECLARE_ALIGNED(16, char, blackclamp[16]);
  DECLARE_ALIGNED(16, char, whiteclamp[16]);
  DECLARE_ALIGNED(16, char, bothclamp[16]);
  const int width  = 64;
  const int height = 64;
  const int image_size = width * height;
  char noise[3072];

  const int clamp = vpx_setup_noise(sizeof(noise), 4.4, noise);

  for (int i = 0; i < 16; i++) {
    blackclamp[i] = clamp;
    whiteclamp[i] = clamp;
    bothclamp[i] = 2 * clamp;
  }

  uint8_t *const s = reinterpret_cast<uint8_t *>(vpx_calloc(image_size, 1));
  uint8_t *const d = reinterpret_cast<uint8_t *>(vpx_calloc(image_size, 1));

  memset(s, 99, image_size);
  memset(d, 99, image_size);

  srand(0);
  ASM_REGISTER_STATE_CHECK(GetParam()(s, noise, blackclamp, whiteclamp,
                                      bothclamp, width, height, width));
  srand(0);
  ASM_REGISTER_STATE_CHECK(vpx_plane_add_noise_c(d, noise, blackclamp,
                                                 whiteclamp, bothclamp,
                                                 width, height, width));

  for (int i = 0; i < image_size; ++i) {
    EXPECT_EQ((int)s[i], (int)d[i]) << "i = " << i;
  }

  vpx_free(d);
  vpx_free(s);
}

INSTANTIATE_TEST_CASE_P(C, AddNoiseTest,
                        ::testing::Values(vpx_plane_add_noise_c));

#if HAVE_SSE2
INSTANTIATE_TEST_CASE_P(SSE2, AddNoiseTest,
                        ::testing::Values(vpx_plane_add_noise_sse2));
#endif

#if HAVE_MSA
INSTANTIATE_TEST_CASE_P(MSA, AddNoiseTest,
                        ::testing::Values(vpx_plane_add_noise_msa));
#endif
}  // namespace
