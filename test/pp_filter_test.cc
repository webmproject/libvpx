/*
 *  Copyright (c) 2012 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */
#include "./vpx_config.h"
#include "./vpx_dsp_rtcd.h"
#include "test/acm_random.h"
#include "test/clear_system_state.h"
#include "test/register_state_check.h"
#include "third_party/googletest/src/include/gtest/gtest.h"
#include "vpx/vpx_integer.h"
#include "vpx_mem/vpx_mem.h"

using libvpx_test::ACMRandom;

typedef void (*VpxPostProcDownAndAcrossMbRowFunc)(
    unsigned char *src_ptr, unsigned char *dst_ptr, int src_pixels_per_line,
    int dst_pixels_per_line, int cols, unsigned char *flimit, int size);

typedef void (*VpxMbPostProcAcrossIpFunc)(unsigned char *src, int pitch,
                                          int rows, int cols, int flimit);

typedef void (*VpxMbPostProcDownFunc)(unsigned char *dst, int pitch, int rows,
                                      int cols, int flimit);

namespace {

// Compute the filter level used in post proc from the loop filter strength
int q2mbl(int x) {
  if (x < 20) x = 20;

  x = 50 + (x - 50) * 10 / 8;
  return x * x / 3;
}

class VpxPostProcDownAndAcrossMbRowTest
    : public ::testing::TestWithParam<VpxPostProcDownAndAcrossMbRowFunc> {
 public:
  virtual void TearDown() { libvpx_test::ClearSystemState(); }
};

// Test routine for the VPx post-processing function
// vpx_post_proc_down_and_across_mb_row_c.

TEST_P(VpxPostProcDownAndAcrossMbRowTest, CheckFilterOutput) {
  // Size of the underlying data block that will be filtered.
  const int block_width = 16;
  const int block_height = 16;

  // 5-tap filter needs 2 padding rows above and below the block in the input.
  const int input_width = block_width;
  const int input_height = block_height + 4;
  const int input_stride = input_width;
  const int input_size = input_width * input_height;

  // Filter extends output block by 8 samples at left and right edges.
  const int output_width = block_width + 16;
  const int output_height = block_height;
  const int output_stride = output_width;
  const int output_size = output_width * output_height;

  uint8_t *const src_image = new uint8_t[input_size];
  ASSERT_TRUE(src_image != NULL);

  // Though the left padding is only 8 bytes, the assembly code tries to
  // read 16 bytes before the pointer.
  uint8_t *const dst_image = new uint8_t[output_size + 8];
  ASSERT_TRUE(dst_image != NULL);

  // Pointers to top-left pixel of block in the input and output images.
  uint8_t *const src_image_ptr = src_image + (input_stride << 1);

  // The assembly works in increments of 16. The first read may be offset by
  // this amount.
  uint8_t *const dst_image_ptr = dst_image + 16;
  uint8_t *const flimits =
      reinterpret_cast<uint8_t *>(vpx_memalign(16, block_width));
  (void)memset(flimits, 255, block_width);

  // Initialize pixels in the input:
  //   block pixels to value 1,
  //   border pixels to value 10.
  (void)memset(src_image, 10, input_size);
  uint8_t *pixel_ptr = src_image_ptr;
  for (int i = 0; i < block_height; ++i) {
    for (int j = 0; j < block_width; ++j) {
      pixel_ptr[j] = 1;
    }
    pixel_ptr += input_stride;
  }

  // Initialize pixels in the output to 99.
  (void)memset(dst_image, 99, output_size);

  ASM_REGISTER_STATE_CHECK(GetParam()(src_image_ptr, dst_image_ptr,
                                      input_stride, output_stride, block_width,
                                      flimits, 16));

  static const uint8_t kExpectedOutput[block_height] = {
    4, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 3, 4
  };

  pixel_ptr = dst_image_ptr;
  for (int i = 0; i < block_height; ++i) {
    for (int j = 0; j < block_width; ++j) {
      ASSERT_EQ(kExpectedOutput[i], pixel_ptr[j]) << "at (" << i << ", " << j
                                                  << ")";
    }
    pixel_ptr += output_stride;
  }

  delete[] src_image;
  delete[] dst_image;
  vpx_free(flimits);
};

class VpxMbPostProcAcrossIpTest
    : public ::testing::TestWithParam<VpxMbPostProcAcrossIpFunc> {
 public:
  virtual void TearDown() { libvpx_test::ClearSystemState(); }
};

TEST_P(VpxMbPostProcAcrossIpTest, CheckFilterOutput) {
  const int rows = 16;
  const int cols = 16;
  const int src_left_padding = 8;
  const int src_right_padding = 17;
  const int src_width = cols + src_left_padding + src_right_padding;
  const int src_size = rows * src_width;

  unsigned char *const src = new unsigned char[src_size];
  ASSERT_TRUE(src != NULL);
  memset(src, 10, src_size);
  unsigned char *s = src + src_left_padding;
  for (int r = 0; r < rows; r++) {
    for (int c = 0; c < cols; c++) {
      s[c] = c;
    }
    s += src_width;
  }

  s = src + src_left_padding;

  ASM_REGISTER_STATE_CHECK(GetParam()(s, src_width, rows, cols, q2mbl(100)));

  static const unsigned char kExpectedOutput[cols] = {
    2, 2, 3, 4, 4, 5, 6, 7, 8, 9, 10, 11, 11, 12, 13, 13
  };
  s = src + src_left_padding;
  for (int r = 0; r < rows; r++) {
    for (int c = 0; c < cols; c++) {
      ASSERT_EQ(kExpectedOutput[c], s[c]) << "at (" << r << ", " << c << ")";
    }
    s += src_width;
  }

  delete[] src;
}

class VpxMbPostProcDownTest
    : public ::testing::TestWithParam<VpxMbPostProcDownFunc> {
 public:
  virtual void TearDown() { libvpx_test::ClearSystemState(); }

 protected:
  virtual void SetRows(unsigned char *src_c, int rows, int cols) {
    for (int r = 0; r < rows; r++) {
      memset(src_c, r, cols);
      src_c += cols;
    }
  }

  virtual void SetRandom(unsigned char *src_c, unsigned char *src_asm, int rows,
                         int cols, int src_pitch) {
    ACMRandom rnd;
    rnd.Reset(ACMRandom::DeterministicSeed());

    // Add some random noise to the input
    for (int r = 0; r < rows; r++) {
      for (int c = 0; c < cols; c++) {
        const int noise = rnd(4);
        src_c[c] = r + noise;
        src_asm[c] = r + noise;
      }
      src_c += src_pitch;
      src_asm += src_pitch;
    }
  }

  virtual void SetRandomSaturation(unsigned char *src_c, unsigned char *src_asm,
                                   int rows, int cols, int src_pitch) {
    ACMRandom rnd;
    rnd.Reset(ACMRandom::DeterministicSeed());

    // Add some random noise to the input
    for (int r = 0; r < rows; r++) {
      for (int c = 0; c < cols; c++) {
        const int noise = 3 * rnd(2);
        src_c[c] = r + noise;
        src_asm[c] = r + noise;
      }
      src_c += src_pitch;
      src_asm += src_pitch;
    }
  }

  virtual void RunComparison(const unsigned char *kExpectedOutput,
                             unsigned char *src_c, int rows, int cols,
                             int src_pitch) {
    for (int r = 0; r < rows; r++) {
      for (int c = 0; c < cols; c++) {
        ASSERT_EQ(kExpectedOutput[r * rows + c], src_c[c]) << "at (" << r
                                                           << ", " << c << ")";
      }
      src_c += src_pitch;
    }
  }

  virtual void RunComparison(unsigned char *src_c, unsigned char *src_asm,
                             int rows, int cols, int src_pitch) {
    for (int r = 0; r < rows; r++) {
      for (int c = 0; c < cols; c++) {
        ASSERT_EQ(src_c[c], src_asm[c]) << "at (" << r << ", " << c << ")";
      }
      src_c += src_pitch;
      src_asm += src_pitch;
    }
  }
};

TEST_P(VpxMbPostProcDownTest, CheckFilterOutput) {
  const int rows = 16;
  const int cols = 16;
  const int src_pitch = cols;
  const int src_top_padding = 8;
  const int src_btm_padding = 17;

  const int src_size = cols * (rows + src_top_padding + src_btm_padding);
  unsigned char *c_mem = new unsigned char[src_size];
  ASSERT_TRUE(c_mem != NULL);
  memset(c_mem, 10, src_size);
  unsigned char *const src_c = c_mem + src_top_padding * src_pitch;

  SetRows(src_c, rows, cols);
  ASM_REGISTER_STATE_CHECK(
      GetParam()(src_c, src_pitch, rows, cols, q2mbl(100)));

  static const unsigned char kExpectedOutput[rows * cols] = {
    2,  2,  1,  1,  2,  2,  2,  2,  2,  2,  1,  1,  2,  2,  2,  2,  2,  2,  2,
    2,  3,  2,  2,  2,  2,  2,  2,  2,  3,  2,  2,  2,  3,  3,  3,  3,  3,  3,
    3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  4,  4,  3,  4,  4,  3,  3,  3,
    4,  4,  3,  4,  4,  3,  3,  4,  5,  4,  4,  4,  4,  4,  4,  4,  5,  4,  4,
    4,  4,  4,  4,  5,  5,  5,  5,  5,  5,  5,  5,  5,  5,  5,  5,  5,  5,  5,
    5,  6,  6,  6,  6,  6,  6,  6,  6,  6,  6,  6,  6,  6,  6,  6,  6,  7,  7,
    7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  8,  8,  8,  8,  8,
    8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  9,  8,  9,  9,  8,  8,  8,  9,
    9,  8,  9,  9,  8,  8,  8,  9,  9,  10, 10, 9,  9,  9,  10, 10, 9,  10, 10,
    9,  9,  9,  10, 10, 10, 11, 10, 10, 10, 11, 10, 11, 10, 11, 10, 10, 10, 11,
    10, 11, 11, 11, 11, 11, 11, 11, 12, 11, 11, 11, 11, 11, 11, 11, 12, 11, 12,
    12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 13, 12,
    13, 12, 13, 12, 12, 12, 13, 12, 13, 12, 13, 12, 13, 13, 13, 14, 13, 13, 13,
    13, 13, 13, 13, 14, 13, 13, 13, 13
  };

  RunComparison(kExpectedOutput, src_c, rows, cols, src_pitch);

  delete[] c_mem;
}

TEST_P(VpxMbPostProcDownTest, CheckCvsAssembly) {
  const int rows = 16;
  const int cols = 16;
  const int src_pitch = cols;
  const int src_top_padding = 8;
  const int src_btm_padding = 17;
  const int src_size = cols * (rows + src_top_padding + src_btm_padding);
  unsigned char *c_mem = new unsigned char[src_size];
  unsigned char *asm_mem = new unsigned char[src_size];
  memset(c_mem, 10, src_size);
  memset(asm_mem, 10, src_size);
  unsigned char *const src_c = c_mem + src_top_padding * src_pitch;
  unsigned char *const src_asm = asm_mem + src_top_padding * src_pitch;

  SetRandom(src_c, src_asm, rows, cols, src_pitch);
  vpx_mbpost_proc_down_c(src_c, src_pitch, rows, cols, q2mbl(100));
  ASM_REGISTER_STATE_CHECK(
      GetParam()(src_asm, src_pitch, rows, cols, q2mbl(100)));
  RunComparison(src_c, src_asm, rows, cols, src_pitch);

  SetRandomSaturation(src_c, src_asm, rows, cols, src_pitch);
  vpx_mbpost_proc_down_c(src_c, src_pitch, rows, cols, q2mbl(100));
  ASM_REGISTER_STATE_CHECK(
      GetParam()(src_asm, src_pitch, rows, cols, q2mbl(100)));
  RunComparison(src_c, src_asm, rows, cols, src_pitch);

  delete[] c_mem;
  delete[] asm_mem;
}

INSTANTIATE_TEST_CASE_P(
    C, VpxPostProcDownAndAcrossMbRowTest,
    ::testing::Values(vpx_post_proc_down_and_across_mb_row_c));

INSTANTIATE_TEST_CASE_P(C, VpxMbPostProcAcrossIpTest,
                        ::testing::Values(vpx_mbpost_proc_across_ip_c));

INSTANTIATE_TEST_CASE_P(C, VpxMbPostProcDownTest,
                        ::testing::Values(vpx_mbpost_proc_down_c));

#if HAVE_SSE2
INSTANTIATE_TEST_CASE_P(
    SSE2, VpxPostProcDownAndAcrossMbRowTest,
    ::testing::Values(vpx_post_proc_down_and_across_mb_row_sse2));

INSTANTIATE_TEST_CASE_P(SSE2, VpxMbPostProcAcrossIpTest,
                        ::testing::Values(vpx_mbpost_proc_across_ip_xmm));

INSTANTIATE_TEST_CASE_P(SSE2, VpxMbPostProcDownTest,
                        ::testing::Values(vpx_mbpost_proc_down_xmm));
#endif

#if HAVE_MSA
INSTANTIATE_TEST_CASE_P(
    MSA, VpxPostProcDownAndAcrossMbRowTest,
    ::testing::Values(vpx_post_proc_down_and_across_mb_row_msa));

INSTANTIATE_TEST_CASE_P(MSA, VpxMbPostProcAcrossIpTest,
                        ::testing::Values(vpx_mbpost_proc_across_ip_msa));

INSTANTIATE_TEST_CASE_P(MSA, VpxMbPostProcDownTest,
                        ::testing::Values(vpx_mbpost_proc_down_msa));
#endif

}  // namespace
