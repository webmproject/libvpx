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

#include <climits>
#include "third_party/googletest/src/include/gtest/gtest.h"
#include "test/codec_factory.h"
#include "test/encode_test_driver.h"
#include "test/i420_video_source.h"
#include "test/util.h"

namespace {

const int kTestMode = 0;
const int kSuperframeSyntax = 1;
const int kTileCols = 2;
const int kTileRows = 3;

typedef std::tr1::tuple<libaom_test::TestMode, int, int, int>
    SuperframeTestParam;

class SuperframeTest
    : public ::libaom_test::EncoderTest,
      public ::libaom_test::CodecTestWithParam<SuperframeTestParam> {
 protected:
  SuperframeTest()
      : EncoderTest(GET_PARAM(0)), modified_buf_(NULL), last_sf_pts_(0) {}
  virtual ~SuperframeTest() {}

  virtual void SetUp() {
    InitializeConfig();
    const SuperframeTestParam input = GET_PARAM(1);
    const libaom_test::TestMode mode = std::tr1::get<kTestMode>(input);
    const int syntax = std::tr1::get<kSuperframeSyntax>(input);
    SetMode(mode);
    sf_count_ = 0;
    sf_count_max_ = INT_MAX;
    is_av1_style_superframe_ = syntax;
    n_tile_cols_ = std::tr1::get<kTileCols>(input);
    n_tile_rows_ = std::tr1::get<kTileRows>(input);
  }

  virtual void TearDown() { delete[] modified_buf_; }

  virtual void PreEncodeFrameHook(libaom_test::VideoSource *video,
                                  libaom_test::Encoder *encoder) {
    if (video->frame() == 1) {
      encoder->Control(AOME_SET_ENABLEAUTOALTREF, 1);
      encoder->Control(AOME_SET_CPUUSED, 2);
      encoder->Control(AV1E_SET_TILE_COLUMNS, n_tile_cols_);
      encoder->Control(AV1E_SET_TILE_ROWS, n_tile_rows_);
    }
  }

  virtual const aom_codec_cx_pkt_t *MutateEncoderOutputHook(
      const aom_codec_cx_pkt_t *pkt) {
    if (pkt->kind != AOM_CODEC_CX_FRAME_PKT) return pkt;

    const uint8_t *buffer = reinterpret_cast<uint8_t *>(pkt->data.frame.buf);
    const uint8_t marker = buffer[pkt->data.frame.sz - 1];
    const int frames = (marker & 0x7) + 1;
    const int mag = ((marker >> 3) & 3) + 1;
    const unsigned int index_sz = 2 + mag * (frames - is_av1_style_superframe_);
    if ((marker & 0xe0) == 0xc0 && pkt->data.frame.sz >= index_sz &&
        buffer[pkt->data.frame.sz - index_sz] == marker) {
      // frame is a superframe. strip off the index.
      if (modified_buf_) delete[] modified_buf_;
      modified_buf_ = new uint8_t[pkt->data.frame.sz - index_sz];
      memcpy(modified_buf_, pkt->data.frame.buf, pkt->data.frame.sz - index_sz);
      modified_pkt_ = *pkt;
      modified_pkt_.data.frame.buf = modified_buf_;
      modified_pkt_.data.frame.sz -= index_sz;

      sf_count_++;
      last_sf_pts_ = pkt->data.frame.pts;
      return &modified_pkt_;
    }

    // Make sure we do a few frames after the last SF
    abort_ |=
        sf_count_ > sf_count_max_ && pkt->data.frame.pts - last_sf_pts_ >= 5;
    return pkt;
  }

  int is_av1_style_superframe_;
  int sf_count_;
  int sf_count_max_;
  aom_codec_cx_pkt_t modified_pkt_;
  uint8_t *modified_buf_;
  aom_codec_pts_t last_sf_pts_;

 private:
  int n_tile_cols_;
  int n_tile_rows_;
};

TEST_P(SuperframeTest, TestSuperframeIndexIsOptional) {
  sf_count_max_ = 0;  // early exit on successful test.
  cfg_.g_lag_in_frames = 25;

  ::libaom_test::I420VideoSource video("hantro_collage_w352h288.yuv", 352, 288,
                                       30, 1, 0, 40);
  ASSERT_NO_FATAL_FAILURE(RunLoop(&video));
#if CONFIG_EXT_REFS
  // NOTE: The use of BWDREF_FRAME will enable the coding of more non-show
  //       frames besides ALTREF_FRAME.
  EXPECT_GE(sf_count_, 1);
#else
  EXPECT_EQ(sf_count_, 1);
#endif  // CONFIG_EXT_REFS
}

// The superframe index is currently mandatory with ANS due to the decoder
// starting at the end of the buffer.
#if CONFIG_EXT_TILE
// Single tile does not work with ANS (see comment above).
#if CONFIG_ANS
const int tile_col_values[] = { 1, 2 };
#else
const int tile_col_values[] = { 1, 2, 32 };
#endif
const int tile_row_values[] = { 1, 2, 32 };
AV1_INSTANTIATE_TEST_CASE(
    SuperframeTest,
    ::testing::Combine(::testing::Values(::libaom_test::kTwoPassGood),
                       ::testing::Values(1),
                       ::testing::ValuesIn(tile_col_values),
                       ::testing::ValuesIn(tile_row_values)));
#else
#if !CONFIG_ANS && !CONFIG_DAALA_EC
AV1_INSTANTIATE_TEST_CASE(
    SuperframeTest,
    ::testing::Combine(::testing::Values(::libaom_test::kTwoPassGood),
                       ::testing::Values(1), ::testing::Values(0),
                       ::testing::Values(0)));
#endif  // !CONFIG_ANS
#endif  // CONFIG_EXT_TILE
}  // namespace
