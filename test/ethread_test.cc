/*
 *  Copyright (c) 2014 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include <string>
#include <vector>
#include "third_party/googletest/src/include/gtest/gtest.h"
#include "test/codec_factory.h"
#include "test/encode_test_driver.h"
#include "test/md5_helper.h"
#include "test/util.h"
#include "test/y4m_video_source.h"

namespace {
class VPxEncoderThreadTest
    : public ::libvpx_test::EncoderTest,
      public ::libvpx_test::CodecTestWith2Params<libvpx_test::TestMode, int> {
 protected:
  VPxEncoderThreadTest()
      : EncoderTest(GET_PARAM(0)),
        encoder_initialized_(false),
        encoding_mode_(GET_PARAM(1)),
        set_cpu_used_(GET_PARAM(2)) {
    init_flags_ = VPX_CODEC_USE_PSNR;
    vpx_codec_dec_cfg_t cfg = vpx_codec_dec_cfg_t();
    cfg.w = 1280;
    cfg.h = 720;
    decoder_ = codec_->CreateDecoder(cfg, 0);
#if CONFIG_VP10 && CONFIG_EXT_TILE
    if (decoder_->IsVP10()) {
      decoder_->Control(VP10_SET_DECODE_TILE_ROW, -1);
      decoder_->Control(VP10_SET_DECODE_TILE_COL, -1);
    }
#endif

    size_enc_.clear();
    md5_dec_.clear();
    md5_enc_.clear();
  }
  virtual ~VPxEncoderThreadTest() {
    delete decoder_;
  }

  virtual void SetUp() {
    InitializeConfig();
    SetMode(encoding_mode_);

    if (encoding_mode_ != ::libvpx_test::kRealTime) {
      cfg_.g_lag_in_frames = 3;
      cfg_.rc_end_usage = VPX_VBR;
      cfg_.rc_2pass_vbr_minsection_pct = 5;
      cfg_.rc_2pass_vbr_maxsection_pct = 2000;
    } else {
      cfg_.g_lag_in_frames = 0;
      cfg_.rc_end_usage = VPX_CBR;
      cfg_.g_error_resilient = 1;
    }
    cfg_.rc_max_quantizer = 56;
    cfg_.rc_min_quantizer = 0;
  }

  virtual void BeginPassHook(unsigned int /*pass*/) {
    encoder_initialized_ = false;
  }

  virtual void PreEncodeFrameHook(::libvpx_test::VideoSource * /*video*/,
                                  ::libvpx_test::Encoder *encoder) {
    if (!encoder_initialized_) {
#if CONFIG_VP10 && CONFIG_EXT_TILE
      encoder->Control(VP9E_SET_TILE_COLUMNS, 1);
      if (codec_ == &libvpx_test::kVP10) {
        // TODO(geza): Start using multiple tile rows when the multi-threaded
        // encoder can handle them
        encoder->Control(VP9E_SET_TILE_ROWS, 32);
      } else {
        encoder->Control(VP9E_SET_TILE_ROWS, 0);
      }
#else
      // Encode 4 tile columns.
      encoder->Control(VP9E_SET_TILE_COLUMNS, 2);
      encoder->Control(VP9E_SET_TILE_ROWS, 0);
#endif  // CONFIG_VP10 && CONFIG_EXT_TILE
      encoder->Control(VP8E_SET_CPUUSED, set_cpu_used_);
      if (encoding_mode_ != ::libvpx_test::kRealTime) {
        encoder->Control(VP8E_SET_ENABLEAUTOALTREF, 1);
        encoder->Control(VP8E_SET_ARNR_MAXFRAMES, 7);
        encoder->Control(VP8E_SET_ARNR_STRENGTH, 5);
        encoder->Control(VP8E_SET_ARNR_TYPE, 3);
      } else {
        encoder->Control(VP8E_SET_ENABLEAUTOALTREF, 0);
        encoder->Control(VP9E_SET_AQ_MODE, 3);
      }
      encoder_initialized_ = true;
    }
  }

  virtual void FramePktHook(const vpx_codec_cx_pkt_t *pkt) {
    size_enc_.push_back(pkt->data.frame.sz);

    ::libvpx_test::MD5 md5_enc;
    md5_enc.Add(reinterpret_cast<uint8_t*>(pkt->data.frame.buf),
                pkt->data.frame.sz);
    md5_enc_.push_back(md5_enc.Get());

    const vpx_codec_err_t res = decoder_->DecodeFrame(
        reinterpret_cast<uint8_t*>(pkt->data.frame.buf), pkt->data.frame.sz);
    if (res != VPX_CODEC_OK) {
      abort_ = true;
      ASSERT_EQ(VPX_CODEC_OK, res);
    }
    const vpx_image_t *img = decoder_->GetDxData().Next();

    if (img) {
      ::libvpx_test::MD5 md5_res;
      md5_res.Add(img);
      md5_dec_.push_back(md5_res.Get());
    }
  }

  void DoTest() {
    ::libvpx_test::Y4mVideoSource video("niklas_1280_720_30.y4m", 15, 18);
    cfg_.rc_target_bitrate = 1000;

    // Encode using single thread.
    cfg_.g_threads = 1;
    init_flags_ = VPX_CODEC_USE_PSNR;
    ASSERT_NO_FATAL_FAILURE(RunLoop(&video));
    std::vector<size_t> single_thr_size_enc;
    std::vector<std::string> single_thr_md5_enc;
    std::vector<std::string> single_thr_md5_dec;
    single_thr_size_enc = size_enc_;
    single_thr_md5_enc = md5_enc_;
    single_thr_md5_dec = md5_dec_;
    size_enc_.clear();
    md5_enc_.clear();
    md5_dec_.clear();

    // Encode using multiple threads.
    cfg_.g_threads = 4;
    ASSERT_NO_FATAL_FAILURE(RunLoop(&video));
    std::vector<size_t> multi_thr_size_enc;
    std::vector<std::string> multi_thr_md5_enc;
    std::vector<std::string> multi_thr_md5_dec;
    multi_thr_size_enc = size_enc_;
    multi_thr_md5_enc = md5_enc_;
    multi_thr_md5_dec = md5_dec_;
    size_enc_.clear();
    md5_enc_.clear();
    md5_dec_.clear();

    // Check that the vectors are equal.
    ASSERT_EQ(single_thr_size_enc, multi_thr_size_enc);
    ASSERT_EQ(single_thr_md5_enc, multi_thr_md5_enc);
    ASSERT_EQ(single_thr_md5_dec, multi_thr_md5_dec);
  }

  bool encoder_initialized_;
  ::libvpx_test::TestMode encoding_mode_;
  int set_cpu_used_;
  ::libvpx_test::Decoder *decoder_;
  std::vector<size_t> size_enc_;
  std::vector<std::string> md5_enc_;
  std::vector<std::string> md5_dec_;
};

TEST_P(VPxEncoderThreadTest, EncoderResultTest) {
  DoTest();
}

class VPxEncoderThreadTestLarge : public VPxEncoderThreadTest {};

TEST_P(VPxEncoderThreadTestLarge, EncoderResultTest) {
  DoTest();
}


VP10_INSTANTIATE_TEST_CASE(
    VPxEncoderThreadTest,
    ::testing::Values(::libvpx_test::kTwoPassGood, ::libvpx_test::kOnePassGood),
    ::testing::Range(3, 9));

VP10_INSTANTIATE_TEST_CASE(
    VPxEncoderThreadTestLarge,
    ::testing::Values(::libvpx_test::kTwoPassGood, ::libvpx_test::kOnePassGood),
    ::testing::Range(1, 3));
}  // namespace
