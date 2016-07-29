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
#if CONFIG_VPX_HIGHBITDEPTH
typedef void (*hbd_conv_filter_t)(const uint16_t*, int, uint16_t*, int,
                                  int, int, const InterpFilterParams,
                                  const int, int, int, int);
#endif

// Test parameter list:
//  <convolve_horiz_func, convolve_vert_func,
//  <width, height>, filter_params, subpel_x_q4, avg>
typedef tuple<int, int> BlockDimension;
typedef tuple<conv_filter_t, conv_filter_t, BlockDimension, INTERP_FILTER,
              int, int> ConvParams;
#if CONFIG_VPX_HIGHBITDEPTH
// Test parameter list:
//  <convolve_horiz_func, convolve_vert_func,
//  <width, height>, filter_params, subpel_x_q4, avg, bit_dpeth>
typedef tuple<hbd_conv_filter_t, hbd_conv_filter_t, BlockDimension,
              INTERP_FILTER, int, int, int> HbdConvParams;
#endif

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
  void PrepFilterBuffer(int w, int h);
  void DiffFilterBuffer();
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

void VP10ConvolveOptimzTest::PrepFilterBuffer(int w, int h) {
  int r, c;
  ACMRandom rnd(ACMRandom::DeterministicSeed());

  memset(alloc_, 0, 4 * maxBlockSize * sizeof(alloc_[0]));

  uint8_t *src_ptr = src_;
  uint8_t *dst_ptr = dst_;
  uint8_t *src_ref_ptr = src_ref_;
  uint8_t *dst_ref_ptr = dst_ref_;

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

void VP10ConvolveOptimzTest::DiffFilterBuffer() {
  int r, c;
  const uint8_t *dst_ptr = dst_;
  const uint8_t *dst_ref_ptr = dst_ref_;
  for (r = 0; r < height_; ++r) {
    for (c = 0; c < width_; ++c) {
      EXPECT_EQ((uint8_t)dst_ref_ptr[c], (uint8_t)dst_ptr[c])
      << "Error at row: " << r << " col: " << c << " "
      << "w = " << width_ << " " << "h = " << height_ << " "
      << "filter group index = " << filter_ << " "
      << "filter index = " << subpel_;
    }
    dst_ptr += stride;
    dst_ref_ptr += stride;
  }
}

void VP10ConvolveOptimzTest::RunHorizFilterBitExactCheck() {
  PrepFilterBuffer(testMaxBlk, testMaxBlk);

  InterpFilterParams filter_params = vp10_get_interp_filter_params(filter_);

  vp10_convolve_horiz_c(src_ref_, stride, dst_ref_, stride, width_, height_,
                        filter_params, subpel_, x_step_q4, avg_);

  conv_horiz_(src_, stride, dst_, stride, width_, height_,
              filter_params, subpel_, x_step_q4, avg_);

  DiffFilterBuffer();

  // Note:
  // Here we need calculate a height which is different from the specified one
  // and test again.
  int intermediate_height =
      (((height_ - 1) * 16 + subpel_) >> SUBPEL_BITS) + filter_params.taps;
  PrepFilterBuffer(testMaxBlk, testMaxBlk);

  vp10_convolve_horiz_c(src_ref_, stride, dst_ref_, stride, width_,
                        intermediate_height, filter_params, subpel_, x_step_q4,
                        avg_);

  conv_horiz_(src_, stride, dst_, stride, width_,
              intermediate_height, filter_params, subpel_, x_step_q4,
              avg_);

  DiffFilterBuffer();
}

void VP10ConvolveOptimzTest::RunVertFilterBitExactCheck() {
  PrepFilterBuffer(testMaxBlk, testMaxBlk);

  InterpFilterParams filter_params = vp10_get_interp_filter_params(filter_);

  vp10_convolve_vert_c(src_ref_, stride, dst_ref_, stride, width_, height_,
                       filter_params, subpel_, x_step_q4, avg_);

  conv_vert_(src_, stride, dst_, stride, width_, height_,
             filter_params, subpel_, x_step_q4, avg_);

  DiffFilterBuffer();
}

TEST_P(VP10ConvolveOptimzTest, HorizBitExactCheck) {
  RunHorizFilterBitExactCheck();
}
TEST_P(VP10ConvolveOptimzTest, VerticalBitExactCheck) {
  RunVertFilterBitExactCheck();
}

using std::tr1::make_tuple;

#if (HAVE_SSSE3 || HAVE_SSE4_1) && CONFIG_EXT_INTERP
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
#endif

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

#if CONFIG_VPX_HIGHBITDEPTH
typedef ::testing::TestWithParam<HbdConvParams> TestWithHbdConvParams;
class VP10HbdConvolveOptimzTest : public TestWithHbdConvParams {
 public:
  virtual ~VP10HbdConvolveOptimzTest() {}
  virtual void SetUp() {
    conv_horiz_ = GET_PARAM(0);
    conv_vert_ = GET_PARAM(1);
    BlockDimension block = GET_PARAM(2);
    width_ = std::tr1::get<0>(block);
    height_ = std::tr1::get<1>(block);
    filter_ = GET_PARAM(3);
    subpel_ = GET_PARAM(4);
    avg_ = GET_PARAM(5);
    bit_depth_ = GET_PARAM(6);

    alloc_ = new uint16_t[maxBlockSize * 4];
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
  void PrepFilterBuffer(int w, int h);
  void DiffFilterBuffer();
  hbd_conv_filter_t conv_horiz_;
  hbd_conv_filter_t conv_vert_;
  uint16_t *alloc_;
  uint16_t *src_;
  uint16_t *dst_;
  uint16_t *src_ref_;
  uint16_t *dst_ref_;
  int width_;
  int height_;
  int filter_;
  int subpel_;
  int avg_;
  int bit_depth_;
};

void VP10HbdConvolveOptimzTest::PrepFilterBuffer(int w, int h) {
  int r, c;
  ACMRandom rnd(ACMRandom::DeterministicSeed());

  memset(alloc_, 0, 4 * maxBlockSize * sizeof(alloc_[0]));

  uint16_t *src_ptr = src_;
  uint16_t *dst_ptr = dst_;
  uint16_t *dst_ref_ptr = dst_ref_;
  uint16_t hbd_mask = (1 << bit_depth_) - 1;

  for (r = 0; r < height_; ++r) {
    for (c = 0; c < width_; ++c) {
      src_ptr[c] = rnd.Rand16() & hbd_mask;
      dst_ptr[c] = rnd.Rand16() & hbd_mask;
      dst_ref_ptr[c] = dst_ptr[c];
    }
    src_ptr += stride;
    dst_ptr += stride;
    dst_ref_ptr += stride;
  }
}

void VP10HbdConvolveOptimzTest::DiffFilterBuffer() {
  int r, c;
  const uint16_t *dst_ptr = dst_;
  const uint16_t *dst_ref_ptr = dst_ref_;
  for (r = 0; r < height_; ++r) {
    for (c = 0; c < width_; ++c) {
      EXPECT_EQ((uint16_t)dst_ref_ptr[c], (uint16_t)dst_ptr[c])
      << "Error at row: " << r << " col: " << c << " "
      << "w = " << width_ << " " << "h = " << height_ << " "
      << "filter group index = " << filter_ << " "
      << "filter index = " << subpel_ << " "
      << "bit depth = " << bit_depth_;
    }
    dst_ptr += stride;
    dst_ref_ptr += stride;
  }
}

void VP10HbdConvolveOptimzTest::RunHorizFilterBitExactCheck() {
  PrepFilterBuffer(testMaxBlk, testMaxBlk);

  InterpFilterParams filter_params = vp10_get_interp_filter_params(filter_);

  vp10_highbd_convolve_horiz_c(src_, stride, dst_ref_, stride, width_,
                               height_, filter_params, subpel_, x_step_q4,
                               avg_, bit_depth_);

  conv_horiz_(src_, stride, dst_, stride, width_, height_,
              filter_params, subpel_, x_step_q4, avg_, bit_depth_);

  DiffFilterBuffer();

  // Note:
  // Here we need calculate a height which is different from the specified one
  // and test again.
  int intermediate_height =
      (((height_ - 1) * 16 + subpel_) >> SUBPEL_BITS) + filter_params.taps;
  PrepFilterBuffer(testMaxBlk, testMaxBlk);

  vp10_highbd_convolve_horiz_c(src_, stride, dst_ref_, stride, width_,
                               intermediate_height, filter_params, subpel_,
                               x_step_q4, avg_, bit_depth_);

  conv_horiz_(src_, stride, dst_, stride, width_, intermediate_height,
              filter_params, subpel_, x_step_q4, avg_, bit_depth_);

  DiffFilterBuffer();
}

void VP10HbdConvolveOptimzTest::RunVertFilterBitExactCheck() {
  PrepFilterBuffer(testMaxBlk, testMaxBlk);

  InterpFilterParams filter_params = vp10_get_interp_filter_params(filter_);

  vp10_highbd_convolve_vert_c(src_, stride, dst_ref_, stride, width_, height_,
                              filter_params, subpel_, x_step_q4, avg_,
                              bit_depth_);

  conv_vert_(src_, stride, dst_, stride, width_, height_,
             filter_params, subpel_, x_step_q4, avg_, bit_depth_);

  DiffFilterBuffer();
}

TEST_P(VP10HbdConvolveOptimzTest, HorizBitExactCheck) {
  RunHorizFilterBitExactCheck();
}
TEST_P(VP10HbdConvolveOptimzTest, VertBitExactCheck) {
  RunVertFilterBitExactCheck();
}

#if HAVE_SSE4_1 && CONFIG_EXT_INTERP

const int kBitdepth[] = {10, 12};

INSTANTIATE_TEST_CASE_P(
    SSE4_1, VP10HbdConvolveOptimzTest,
    ::testing::Combine(
         ::testing::Values(vp10_highbd_convolve_horiz_sse4_1),
         ::testing::Values(vp10_highbd_convolve_vert_sse4_1),
         ::testing::ValuesIn(kBlockDim),
         ::testing::ValuesIn(kFilter),
         ::testing::ValuesIn(kSubpelQ4),
         ::testing::ValuesIn(kAvg),
         ::testing::ValuesIn(kBitdepth)));
#endif  // HAVE_SSE4_1 && CONFIG_EXT_INTERP
#endif  // CONFIG_VPX_HIGHBITDEPTH
}  // namespace
