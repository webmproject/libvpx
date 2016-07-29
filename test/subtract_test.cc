/*
 *  Copyright (c) 2012 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include "third_party/googletest/src/include/gtest/gtest.h"

#include "./vpx_config.h"
#include "./vpx_dsp_rtcd.h"
#include "test/acm_random.h"
#include "test/clear_system_state.h"
#include "test/register_state_check.h"
#include "test/util.h"
#if CONFIG_VP10
#include "vp10/common/blockd.h"
#endif
#include "vpx_mem/vpx_mem.h"
#include "vpx_ports/mem.h"

#define USE_SPEED_TEST (0)

typedef void (*SubtractFunc)(int rows, int cols,
                             int16_t *diff_ptr, ptrdiff_t diff_stride,
                             const uint8_t *src_ptr, ptrdiff_t src_stride,
                             const uint8_t *pred_ptr, ptrdiff_t pred_stride);

namespace {

class VP9SubtractBlockTest : public ::testing::TestWithParam<SubtractFunc> {
 public:
  virtual void TearDown() {
    libvpx_test::ClearSystemState();
  }
};

using libvpx_test::ACMRandom;

TEST_P(VP9SubtractBlockTest, SimpleSubtract) {
  ACMRandom rnd(ACMRandom::DeterministicSeed());

  // FIXME(rbultje) split in its own file
  for (BLOCK_SIZE bsize = BLOCK_4X4; bsize < BLOCK_SIZES;
       bsize = static_cast<BLOCK_SIZE>(static_cast<int>(bsize) + 1)) {
    const int block_width = 4 * num_4x4_blocks_wide_lookup[bsize];
    const int block_height = 4 * num_4x4_blocks_high_lookup[bsize];
    int16_t *diff = reinterpret_cast<int16_t *>(
        vpx_memalign(16, sizeof(*diff) * block_width * block_height * 2));
    uint8_t *pred = reinterpret_cast<uint8_t *>(
        vpx_memalign(16, block_width * block_height * 2));
    uint8_t *src  = reinterpret_cast<uint8_t *>(
        vpx_memalign(16, block_width * block_height * 2));

    for (int n = 0; n < 100; n++) {
      for (int r = 0; r < block_height; ++r) {
        for (int c = 0; c < block_width * 2; ++c) {
          src[r * block_width * 2 + c] = rnd.Rand8();
          pred[r * block_width * 2 + c] = rnd.Rand8();
        }
      }

      GetParam()(block_height, block_width, diff, block_width,
                 src, block_width, pred, block_width);

      for (int r = 0; r < block_height; ++r) {
        for (int c = 0; c < block_width; ++c) {
          EXPECT_EQ(diff[r * block_width + c],
                    (src[r * block_width + c] -
                     pred[r * block_width + c])) << "r = " << r
                                                 << ", c = " << c
                                                 << ", bs = " << bsize;
        }
      }

      GetParam()(block_height, block_width, diff, block_width * 2,
                 src, block_width * 2, pred, block_width * 2);

      for (int r = 0; r < block_height; ++r) {
        for (int c = 0; c < block_width; ++c) {
          EXPECT_EQ(diff[r * block_width * 2 + c],
                    (src[r * block_width * 2 + c] -
                     pred[r * block_width * 2 + c])) << "r = " << r
                                                     << ", c = " << c
                                                     << ", bs = " << bsize;
        }
      }
    }
    vpx_free(diff);
    vpx_free(pred);
    vpx_free(src);
  }
}

INSTANTIATE_TEST_CASE_P(C, VP9SubtractBlockTest,
                        ::testing::Values(vpx_subtract_block_c));

#if HAVE_SSE2
INSTANTIATE_TEST_CASE_P(SSE2, VP9SubtractBlockTest,
                        ::testing::Values(vpx_subtract_block_sse2));
#endif
#if HAVE_NEON
INSTANTIATE_TEST_CASE_P(NEON, VP9SubtractBlockTest,
                        ::testing::Values(vpx_subtract_block_neon));
#endif
#if HAVE_MSA
INSTANTIATE_TEST_CASE_P(MSA, VP9SubtractBlockTest,
                        ::testing::Values(vpx_subtract_block_msa));
#endif

typedef void (*HBDSubtractFunc)(int rows, int cols,
                                int16_t *diff_ptr, ptrdiff_t diff_stride,
                                const uint8_t *src_ptr, ptrdiff_t src_stride,
                                const uint8_t *pred_ptr, ptrdiff_t pred_stride,
                                int bd);

using ::std::tr1::get;
using ::std::tr1::make_tuple;
using ::std::tr1::tuple;

// <width, height, bit_dpeth, subtract>
typedef tuple<int, int, int, HBDSubtractFunc> Params;

#if CONFIG_VPX_HIGHBITDEPTH
class VP10HBDSubtractBlockTest : public ::testing::TestWithParam<Params> {
 public:
  virtual void SetUp() {
    block_width_ = GET_PARAM(0);
    block_height_ = GET_PARAM(1);
    bit_depth_ = static_cast<vpx_bit_depth_t>(GET_PARAM(2));
    func_ = GET_PARAM(3);

    rnd_.Reset(ACMRandom::DeterministicSeed());

    const size_t max_width = 128;
    const size_t max_block_size = max_width * max_width;
    src_ = CONVERT_TO_BYTEPTR(reinterpret_cast<uint16_t *>(
        vpx_memalign(16, max_block_size * sizeof(uint16_t))));
    pred_ = CONVERT_TO_BYTEPTR(reinterpret_cast<uint16_t *>(
        vpx_memalign(16, max_block_size * sizeof(uint16_t))));
    diff_ = reinterpret_cast<int16_t *>(
        vpx_memalign(16, max_block_size * sizeof(int16_t)));
  }

  virtual void TearDown() {
    vpx_free(CONVERT_TO_SHORTPTR(src_));
    vpx_free(CONVERT_TO_SHORTPTR(pred_));
    vpx_free(diff_);
  }

 protected:
  void RunForSpeed();
  void CheckResult();

 private:
  ACMRandom rnd_;
  int block_height_;
  int block_width_;
  vpx_bit_depth_t bit_depth_;
  HBDSubtractFunc func_;
  uint8_t *src_;
  uint8_t *pred_;
  int16_t *diff_;
};

void VP10HBDSubtractBlockTest::RunForSpeed() {
  const int test_num = 200000;
  const int max_width = 128;
  const int max_block_size = max_width * max_width;
  const int mask = (1 << bit_depth_) - 1;
  int i, j;

  for (j = 0; j < max_block_size; ++j) {
    CONVERT_TO_SHORTPTR(src_)[j] = rnd_.Rand16() & mask;
    CONVERT_TO_SHORTPTR(pred_)[j] = rnd_.Rand16() & mask;
  }

  for (i = 0; i < test_num; ++i) {
    func_(block_height_, block_width_, diff_, block_width_,
          src_, block_width_, pred_, block_width_, bit_depth_);
  }
}

void VP10HBDSubtractBlockTest::CheckResult() {
  const int test_num = 100;
  const int max_width = 128;
  const int max_block_size = max_width * max_width;
  const int mask = (1 << bit_depth_) - 1;
  int i, j;

  for (i = 0; i < test_num; ++i) {
    for (j = 0; j < max_block_size; ++j) {
      CONVERT_TO_SHORTPTR(src_)[j] = rnd_.Rand16() & mask;
      CONVERT_TO_SHORTPTR(pred_)[j] = rnd_.Rand16() & mask;
    }

    func_(block_height_, block_width_, diff_, block_width_,
          src_, block_width_, pred_, block_width_, bit_depth_);

    for (int r = 0; r < block_height_; ++r) {
      for (int c = 0; c < block_width_; ++c) {
        EXPECT_EQ(diff_[r * block_width_ + c],
                  (CONVERT_TO_SHORTPTR(src_)[r * block_width_ + c] -
                   CONVERT_TO_SHORTPTR(pred_)[r * block_width_ + c]))
            << "r = " << r << ", c = " << c << ", test: " << i;
      }
    }
  }
}

TEST_P(VP10HBDSubtractBlockTest, CheckResult) {
  CheckResult();
}

#if USE_SPEED_TEST
TEST_P(VP10HBDSubtractBlockTest, CheckSpeed) {
  RunForSpeed();
}
#endif  // USE_SPEED_TEST

#if HAVE_SSE2
INSTANTIATE_TEST_CASE_P(SSE2, VP10HBDSubtractBlockTest, ::testing::Values(
    make_tuple(4, 4, 12, vpx_highbd_subtract_block_sse2),
    make_tuple(4, 4, 12, vpx_highbd_subtract_block_c),
    make_tuple(4, 8, 12, vpx_highbd_subtract_block_sse2),
    make_tuple(4, 8, 12, vpx_highbd_subtract_block_c),
    make_tuple(8, 4, 12, vpx_highbd_subtract_block_sse2),
    make_tuple(8, 4, 12, vpx_highbd_subtract_block_c),
    make_tuple(8, 8, 12, vpx_highbd_subtract_block_sse2),
    make_tuple(8, 8, 12, vpx_highbd_subtract_block_c),
    make_tuple(8, 16, 12, vpx_highbd_subtract_block_sse2),
    make_tuple(8, 16, 12, vpx_highbd_subtract_block_c),
    make_tuple(16, 8, 12, vpx_highbd_subtract_block_sse2),
    make_tuple(16, 8, 12, vpx_highbd_subtract_block_c),
    make_tuple(16, 16, 12, vpx_highbd_subtract_block_sse2),
    make_tuple(16, 16, 12, vpx_highbd_subtract_block_c),
    make_tuple(16, 32, 12, vpx_highbd_subtract_block_sse2),
    make_tuple(16, 32, 12, vpx_highbd_subtract_block_c),
    make_tuple(32, 16, 12, vpx_highbd_subtract_block_sse2),
    make_tuple(32, 16, 12, vpx_highbd_subtract_block_c),
    make_tuple(32, 32, 12, vpx_highbd_subtract_block_sse2),
    make_tuple(32, 32, 12, vpx_highbd_subtract_block_c),
    make_tuple(32, 64, 12, vpx_highbd_subtract_block_sse2),
    make_tuple(32, 64, 12, vpx_highbd_subtract_block_c),
    make_tuple(64, 32, 12, vpx_highbd_subtract_block_sse2),
    make_tuple(64, 32, 12, vpx_highbd_subtract_block_c),
    make_tuple(64, 64, 12, vpx_highbd_subtract_block_sse2),
    make_tuple(64, 64, 12, vpx_highbd_subtract_block_c),
    make_tuple(64, 128, 12, vpx_highbd_subtract_block_sse2),
    make_tuple(64, 128, 12, vpx_highbd_subtract_block_c),
    make_tuple(128, 64, 12, vpx_highbd_subtract_block_sse2),
    make_tuple(128, 64, 12, vpx_highbd_subtract_block_c),
    make_tuple(128, 128, 12, vpx_highbd_subtract_block_sse2),
    make_tuple(128, 128, 12, vpx_highbd_subtract_block_c)));
#endif  // HAVE_SSE2
#endif  // CONFIG_VPX_HIGHBITDEPTH
}  // namespace
