/*
 * Copyright (c) 2016, Alliance for Open Media. All rights reserved
 *
 * This source code is subject to the terms of the BSD 2 Clause License and
 * the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
 * was not distributed with this source code in the LICENSE file, you can
 * obtain it at www.aomedia.org/license/software. If the Alliance for Open
 * Media Patent License 1.0 was not distributed with this source code in the
 * PATENTS file, you can obtain it at www.aomedia.org/license/patent.
*/

#include <cstdlib>
#include <string>

#include "third_party/googletest/src/include/gtest/gtest.h"

#include "./aom_config.h"
#include "./aom_dsp_rtcd.h"
#include "aom_ports/aom_timer.h"
#include "test/acm_random.h"
#include "test/clear_system_state.h"
#include "test/register_state_check.h"
#include "test/util.h"

using libaom_test::ACMRandom;

namespace {

typedef void (*clpf_block_t)(const uint8_t *src, uint8_t *dst, int stride,
                             int x0, int y0, int sizex, int sizey, int width,
                             int height, unsigned int strength);

typedef std::tr1::tuple<clpf_block_t, clpf_block_t, int, int>
    clpf_block_param_t;

class ClpfBlockTest : public ::testing::TestWithParam<clpf_block_param_t> {
 public:
  virtual ~ClpfBlockTest() {}
  virtual void SetUp() {
    clpf = GET_PARAM(0);
    ref_clpf = GET_PARAM(1);
    sizex = GET_PARAM(2);
    sizey = GET_PARAM(3);
  }

  virtual void TearDown() { libaom_test::ClearSystemState(); }

 protected:
  int sizex;
  int sizey;
  clpf_block_t clpf;
  clpf_block_t ref_clpf;
};

typedef ClpfBlockTest ClpfSpeedTest;

TEST_P(ClpfBlockTest, TestSIMDNoMismatch) {
  int w = sizex;
  int h = sizey;
  const int size = 32;
  ACMRandom rnd(ACMRandom::DeterministicSeed());
  DECLARE_ALIGNED(16, uint8_t, s[size * size]);
  DECLARE_ALIGNED(16, uint8_t, d[size * size]);
  DECLARE_ALIGNED(16, uint8_t, ref_d[size * size]);
  memset(ref_d, 0, size * size);
  memset(d, 0, size * size);

  int error = 0;
  int pos = 0;
  int strength = 0;
  int xpos = 0, ypos = 0;
  int bits;
  int level;

  // Test every combination of:
  // * Input with 1-8 bits of noise
  // * Noise level around every value from 0 to 255
  // * Blocks anywhere in the frame (along all egdes and also fully inside)
  // * All strengths
  for (level = 0; level < 256 && !error; level++) {
    for (bits = 1; bits < 9 && !error; bits++) {
      for (int i = 0; i < size * size; i++)
        s[i] = clamp((rnd.Rand8() & ((1 << bits) - 1)) + level, 0, 255);

      for (ypos = 0; ypos < size && !error; ypos += h * !error) {
        for (xpos = 0; xpos < size && !error; xpos += w * !error) {
          for (strength = 0; strength < 3 && !error; strength += !error) {
            ref_clpf(s, ref_d, size, xpos, ypos, w, h, size, size,
                     1 << strength);
            ASM_REGISTER_STATE_CHECK(
                clpf(s, d, size, xpos, ypos, w, h, size, size, 1 << strength));

            for (pos = 0; pos < size * size && !error; pos++) {
              error = ref_d[pos] != d[pos];
            }
          }
        }
      }
    }
  }

  EXPECT_EQ(0, error)
      << "Error: ClpfBlockTest, SIMD and C mismatch." << std::endl
      << "First error at " << pos % size << "," << pos / size << " ("
      << (int16_t)ref_d[pos] << " != " << (int16_t)d[pos] << ") " << std::endl
      << "strength: " << (1 << strength) << std::endl
      << "xpos: " << xpos << std::endl
      << "ypos: " << ypos << std::endl
      << "A=" << (pos > size ? (int16_t)s[pos - size] : -1) << std::endl
      << "B=" << (pos % size - 2 >= 0 ? (int16_t)s[pos - 2] : -1) << std::endl
      << "C=" << (pos % size - 1 >= 0 ? (int16_t)s[pos - 1] : -1) << std::endl
      << "X=" << (int16_t)s[pos] << std::endl
      << "D=" << (pos % size + 1 < size ? (int16_t)s[pos + 1] : -1) << std::endl
      << "E=" << (pos % size + 2 < size ? (int16_t)s[pos + 2] : -1) << std::endl
      << "F=" << (pos + size < size * size ? (int16_t)s[pos + size] : -1)
      << std::endl;
}

TEST_P(ClpfSpeedTest, TestSpeed) {
  int w = sizex;
  int h = sizey;
  const int size = 32;
  ACMRandom rnd(ACMRandom::DeterministicSeed());
  DECLARE_ALIGNED(16, uint8_t, s[size * size]);
  DECLARE_ALIGNED(16, uint8_t, d[size * size]);

  int strength;
  int xpos, ypos;

  for (int i = 0; i < size * size; i++) s[i] = rnd.Rand8();

  aom_usec_timer ref_timer;
  aom_usec_timer timer;

  aom_usec_timer_start(&ref_timer);
  for (int c = 0; c < 65536; c++) {
    for (ypos = 0; ypos < size; ypos += h) {
      for (xpos = 0; xpos < size; xpos += w) {
        for (strength = 0; strength < 3; strength++) {
          ref_clpf(s, d, size, xpos, ypos, w, h, size, size, 1 << strength);
        }
      }
    }
  }
  aom_usec_timer_mark(&ref_timer);
  int ref_elapsed_time = aom_usec_timer_elapsed(&ref_timer);

  aom_usec_timer_start(&timer);
  for (int c = 0; c < 65536; c++) {
    for (ypos = 0; ypos < size; ypos += h) {
      for (xpos = 0; xpos < size; xpos += w) {
        for (strength = 0; strength < 3; strength++) {
          clpf(s, d, size, xpos, ypos, w, h, size, size, 1 << strength);
        }
      }
    }
  }
  aom_usec_timer_mark(&timer);
  int elapsed_time = aom_usec_timer_elapsed(&timer);

#if 0
  std::cout << "[          ] C time = " << ref_elapsed_time / 1000
            << " ms, SIMD time = " << elapsed_time / 1000 << " ms" << std::endl;
#endif

  EXPECT_GT(ref_elapsed_time, elapsed_time)
      << "Error: ClpfSpeedTest, SIMD slower than C." << std::endl
      << "C time: " << ref_elapsed_time << "ms" << std::endl
      << "SIMD time: " << elapsed_time << "ms" << std::endl;
}

using std::tr1::make_tuple;

// Test all supported architectures and block sizes
#if HAVE_SSE2
INSTANTIATE_TEST_CASE_P(
    SSE2, ClpfBlockTest,
    ::testing::Values(make_tuple(&aom_clpf_block_sse2, &aom_clpf_block_c, 8, 8),
                      make_tuple(&aom_clpf_block_sse2, &aom_clpf_block_c, 8, 4),
                      make_tuple(&aom_clpf_block_sse2, &aom_clpf_block_c, 4, 8),
                      make_tuple(&aom_clpf_block_sse2, &aom_clpf_block_c, 4,
                                 4)));
#endif

#if HAVE_SSSE3
INSTANTIATE_TEST_CASE_P(
    SSSE3, ClpfBlockTest,
    ::testing::Values(
        make_tuple(&aom_clpf_block_ssse3, &aom_clpf_block_c, 8, 8),
        make_tuple(&aom_clpf_block_ssse3, &aom_clpf_block_c, 8, 4),
        make_tuple(&aom_clpf_block_ssse3, &aom_clpf_block_c, 4, 8),
        make_tuple(&aom_clpf_block_ssse3, &aom_clpf_block_c, 4, 4)));
#endif

#if HAVE_SSE4_1
INSTANTIATE_TEST_CASE_P(
    SSSE4_1, ClpfBlockTest,
    ::testing::Values(
        make_tuple(&aom_clpf_block_sse4_1, &aom_clpf_block_c, 8, 8),
        make_tuple(&aom_clpf_block_sse4_1, &aom_clpf_block_c, 8, 4),
        make_tuple(&aom_clpf_block_sse4_1, &aom_clpf_block_c, 4, 8),
        make_tuple(&aom_clpf_block_sse4_1, &aom_clpf_block_c, 4, 4)));
#endif

#if HAVE_NEON
INSTANTIATE_TEST_CASE_P(
    NEON, ClpfBlockTest,
    ::testing::Values(make_tuple(&aom_clpf_block_neon, &aom_clpf_block_c, 8, 8),
                      make_tuple(&aom_clpf_block_neon, &aom_clpf_block_c, 8, 4),
                      make_tuple(&aom_clpf_block_neon, &aom_clpf_block_c, 4, 8),
                      make_tuple(&aom_clpf_block_neon, &aom_clpf_block_c, 4,
                                 4)));
#endif

// Test speed for all supported architectures
#if HAVE_SSE2
INSTANTIATE_TEST_CASE_P(SSE2, ClpfSpeedTest,
                        ::testing::Values(make_tuple(&aom_clpf_block_sse2,
                                                     &aom_clpf_block_c, 8, 8)));
#endif

#if HAVE_SSSE3
INSTANTIATE_TEST_CASE_P(SSSE3, ClpfSpeedTest,
                        ::testing::Values(make_tuple(&aom_clpf_block_ssse3,
                                                     &aom_clpf_block_c, 8, 8)));
#endif

#if HAVE_SSE4_1
INSTANTIATE_TEST_CASE_P(SSSE4_1, ClpfSpeedTest,
                        ::testing::Values(make_tuple(&aom_clpf_block_ssse3,
                                                     &aom_clpf_block_c, 8, 8)));
#endif

#if HAVE_NEON
INSTANTIATE_TEST_CASE_P(NEON, ClpfSpeedTest,
                        ::testing::Values(make_tuple(&aom_clpf_block_neon,
                                                     &aom_clpf_block_c, 8, 8)));
#endif
}  // namespace
