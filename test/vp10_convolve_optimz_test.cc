/*
 *  Copyright (c) 2016 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include "third_party/googletest/src/include/gtest/gtest.h"

#include "./vp10_rtcd.h"
#include "test/acm_random.h"
#include "test/clear_system_state.h"
#include "test/register_state_check.h"
#include "test/util.h"

namespace {

using std::tr1::tuple;
using libvpx_test::ACMRandom;

typedef void (*conv_filter_t)(const uint8_t*, int, uint8_t*, int,
                              int, int, const InterpFilterParams,
                              const int, int, int);
// Test parameter list:
//  <convolve_horiz_func, convolve_vert_func,
//  <width, height>, filter_params, subpel_x_q4, avg>
typedef tuple<int, int> BlockDimension;
typedef tuple<conv_filter_t, conv_filter_t, BlockDimension, INTERP_FILTER,
              int, int> ConvParams;

// Note:
//  src_ and src_ref_ have special boundary requirement
//  dst_ and dst_ref_ don't
const size_t maxWidth = 256;
const size_t maxHeight = 256;
const size_t maxBlockSize = maxWidth * maxHeight;
const int horizOffset = 32;
const int vertiOffset = 32;
const size_t testMaxBlk = 128;
const int stride = 128;
const int x_step_q4 = 16;

class VP10ConvolveOptimzTest : public ::testing::TestWithParam<ConvParams> {
 public:
  virtual ~VP10ConvolveOptimzTest() {}
  virtual void SetUp() {
    conv_horiz_ = GET_PARAM(0);
    conv_vert_ = GET_PARAM(1);
    BlockDimension block = GET_PARAM(2);
    width_ = std::tr1::get<0>(block);
    height_ = std::tr1::get<1>(block);
    filter_ = GET_PARAM(3);
    subpel_ = GET_PARAM(4);
    avg_ = GET_PARAM(5);

    alloc_ = new uint8_t[maxBlockSize * 4];
    src_ = alloc_ + (vertiOffset * maxWidth);
    src_ += horizOffset;
    src_ref_ = src_ + maxBlockSize;

    dst_ = alloc_ + 2 * maxBlockSize;
    dst_ref_ = alloc_ + 3 * maxBlockSize;
  }

  virtual void TearDown() {
    delete[] alloc_;
    libvpx_test::ClearSystemState();
  }

 protected:
  void RunHorizFilterBitExactCheck();
  void RunVertFilterBitExactCheck();

 private:
  void PrepFilterBuffer(uint8_t *src, uint8_t *src_ref,
                        uint8_t *dst, uint8_t *dst_ref,
                        int w, int h);
  void DiffFilterBuffer(const uint8_t *buf, const uint8_t *buf_ref,
                        int w, int h, int fgroup, int findex);
  conv_filter_t conv_horiz_;
  conv_filter_t conv_vert_;
  uint8_t *alloc_;
  uint8_t *src_;
  uint8_t *dst_;
  uint8_t *src_ref_;
  uint8_t *dst_ref_;
  int width_;
  int height_;
  int filter_;
  int subpel_;
  int avg_;
};

void VP10ConvolveOptimzTest::PrepFilterBuffer(uint8_t *src, uint8_t *src_ref,
                                              uint8_t *dst, uint8_t *dst_ref,
                                              int w, int h) {
  int r, c;
  ACMRandom rnd(ACMRandom::DeterministicSeed());

  memset(alloc_, 0, 4 * maxBlockSize * sizeof(alloc_[0]));

  uint8_t *src_ptr = src;
  uint8_t *dst_ptr = dst;
  uint8_t *src_ref_ptr = src_ref;
  uint8_t *dst_ref_ptr = dst_ref;

  for (r = 0; r < height_; ++r) {
    for (c = 0; c < width_; ++c) {
      src_ptr[c] = rnd.Rand8();
      src_ref_ptr[c] = src_ptr[c];
      dst_ptr[c] = rnd.Rand8();
      dst_ref_ptr[c] = dst_ptr[c];
    }
    src_ptr += stride;
    src_ref_ptr += stride;
    dst_ptr += stride;
    dst_ref_ptr += stride;
  }
}

void VP10ConvolveOptimzTest::DiffFilterBuffer(const uint8_t *buf,
                                              const uint8_t *buf_ref,
                                              int w, int h,
                                              int filter_group,
                                              int filter_index) {
  int r, c;
  const uint8_t *dst_ptr = buf;
  const uint8_t *dst_ref_ptr = buf_ref;
  for (r = 0; r < h; ++r) {
    for (c = 0; c < w; ++c) {
      EXPECT_EQ((uint8_t)dst_ref_ptr[c], (uint8_t)dst_ptr[c])
      << "Error at row: " << r << " col: " << c << " "
      << "w = " << w << " " << "h = " << h << " "
      << "filter group index = " << filter_group << " "
      << "filter index = " << filter_index;
    }
    dst_ptr += stride;
    dst_ref_ptr += stride;
  }
}

void VP10ConvolveOptimzTest::RunHorizFilterBitExactCheck() {
  PrepFilterBuffer(src_, src_ref_, dst_, dst_ref_, testMaxBlk, testMaxBlk);

  InterpFilterParams filter_params = vp10_get_interp_filter_params(filter_);

  vp10_convolve_horiz_c(src_ref_, stride, dst_ref_, stride, width_, height_,
                        filter_params, subpel_, x_step_q4, avg_);

  conv_horiz_(src_, stride, dst_, stride, width_, height_,
              filter_params, subpel_, x_step_q4, avg_);

  DiffFilterBuffer(dst_, dst_ref_, width_, height_, filter_, subpel_);

  // Note:
  // Here we need calculate a height which is different from the specified one
  // and test again.
  int intermediate_height =
      (((height_ - 1) * 16 + subpel_) >> SUBPEL_BITS) + filter_params.taps;
  PrepFilterBuffer(src_, src_ref_, dst_, dst_ref_, testMaxBlk, testMaxBlk);

  vp10_convolve_horiz_c(src_ref_, stride, dst_ref_, stride, width_,
                        intermediate_height, filter_params, subpel_, x_step_q4,
                        avg_);

  conv_horiz_(src_, stride, dst_, stride, width_,
              intermediate_height, filter_params, subpel_, x_step_q4,
              avg_);

  DiffFilterBuffer(dst_, dst_ref_, width_, intermediate_height, filter_,
                   subpel_);
}

void VP10ConvolveOptimzTest::RunVertFilterBitExactCheck() {
  PrepFilterBuffer(src_, src_ref_, dst_, dst_ref_, testMaxBlk, testMaxBlk);

  InterpFilterParams filter_params = vp10_get_interp_filter_params(filter_);

  vp10_convolve_vert_c(src_ref_, stride, dst_ref_, stride, width_, height_,
                       filter_params, subpel_, x_step_q4, avg_);

  conv_vert_(src_, stride, dst_, stride, width_, height_,
             filter_params, subpel_, x_step_q4, avg_);

  DiffFilterBuffer(dst_, dst_ref_, width_, height_, filter_, subpel_);
}

TEST_P(VP10ConvolveOptimzTest, HorizBitExactCheck) {
  RunHorizFilterBitExactCheck();
}
TEST_P(VP10ConvolveOptimzTest, VerticalBitExactCheck) {
  RunVertFilterBitExactCheck();
}

using std::tr1::make_tuple;

const BlockDimension kBlockDim[] = {
  make_tuple(2, 2),
  make_tuple(2, 4),
  make_tuple(4, 4),
  make_tuple(4, 8),
  make_tuple(8, 4),
  make_tuple(8, 8),
  make_tuple(8, 16),
  make_tuple(16, 8),
  make_tuple(16, 16),
  make_tuple(16, 32),
  make_tuple(32, 16),
  make_tuple(32, 32),
  make_tuple(32, 64),
  make_tuple(64, 32),
  make_tuple(64, 64),
  make_tuple(64, 128),
  make_tuple(128, 64),
  make_tuple(128, 128),
};
// 10/12-tap filters
const INTERP_FILTER kFilter[] = {6, 4, 2};

const int kSubpelQ4[] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};

const int kAvg[] = {0, 1};

#if HAVE_SSSE3 && CONFIG_EXT_INTERP
INSTANTIATE_TEST_CASE_P(
    SSSE3, VP10ConvolveOptimzTest,
    ::testing::Combine(
         ::testing::Values(vp10_convolve_horiz_ssse3),
         ::testing::Values(vp10_convolve_vert_ssse3),
         ::testing::ValuesIn(kBlockDim),
         ::testing::ValuesIn(kFilter),
         ::testing::ValuesIn(kSubpelQ4),
         ::testing::ValuesIn(kAvg)));
#endif  // HAVE_SSSE3 && CONFIG_EXT_INTERP
}  // namespace
