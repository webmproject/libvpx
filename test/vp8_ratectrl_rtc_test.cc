/*
 *  Copyright (c) 2021 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include <fstream>  // NOLINT
#include <string>

#include "./vpx_config.h"
#include "third_party/googletest/src/include/gtest/gtest.h"
#include "test/codec_factory.h"
#include "test/encode_test_driver.h"
#include "test/i420_video_source.h"
#include "test/util.h"
#include "test/video_source.h"
#include "vp8/vp8_ratectrl_rtc.h"
#include "vpx/vpx_codec.h"
#include "vpx_ports/bitops.h"

namespace {

struct Vp8RCTestVideo {
  Vp8RCTestVideo() {}
  Vp8RCTestVideo(const char *name_, int width_, int height_,
                 unsigned int frames_)
      : name(name_), width(width_), height(height_), frames(frames_) {}

  friend std::ostream &operator<<(std::ostream &os,
                                  const Vp8RCTestVideo &video) {
    os << video.name << " " << video.width << " " << video.height << " "
       << video.frames;
    return os;
  }
  const char *name;
  int width;
  int height;
  unsigned int frames;
};

const Vp8RCTestVideo kVp8RCTestVectors[] = {
  Vp8RCTestVideo("niklas_640_480_30.yuv", 640, 480, 470),
  Vp8RCTestVideo("desktop_office1.1280_720-020.yuv", 1280, 720, 300),
};

class Vp8RcInterfaceTest
    : public ::libvpx_test::EncoderTest,
      public ::libvpx_test::CodecTestWith2Params<int, Vp8RCTestVideo> {
 public:
  Vp8RcInterfaceTest()
      : EncoderTest(GET_PARAM(0)), key_interval_(3000), encoder_exit_(false) {}
  virtual ~Vp8RcInterfaceTest() {}

 protected:
  virtual void SetUp() {
    InitializeConfig();
    SetMode(::libvpx_test::kRealTime);
  }

  virtual void PreEncodeFrameHook(::libvpx_test::VideoSource *video,
                                  ::libvpx_test::Encoder *encoder) {
    if (video->frame() == 0) {
      encoder->Control(VP8E_SET_CPUUSED, -6);
      encoder->Control(VP8E_SET_RTC_EXTERNAL_RATECTRL, 1);
      encoder->Control(VP8E_SET_MAX_INTRA_BITRATE_PCT, 1000);
    }
    frame_params_.frame_type =
        video->frame() % key_interval_ == 0 ? KEY_FRAME : INTER_FRAME;
    if (frame_params_.frame_type == INTER_FRAME) {
      // Disable golden frame update.
      frame_flags_ |= VP8_EFLAG_NO_UPD_GF;
      frame_flags_ |= VP8_EFLAG_NO_UPD_ARF;
    }
    encoder_exit_ = video->frame() == test_video_.frames;
  }

  virtual void PostEncodeFrameHook(::libvpx_test::Encoder *encoder) {
    if (encoder_exit_) {
      return;
    }
    int qp;
    encoder->Control(VP8E_GET_LAST_QUANTIZER, &qp);
    rc_api_->ComputeQP(frame_params_);
    ASSERT_EQ(rc_api_->GetQP(), qp);
  }

  virtual void FramePktHook(const vpx_codec_cx_pkt_t *pkt) {
    rc_api_->PostEncodeUpdate(pkt->data.frame.sz);
  }

  void RunOneLayer() {
    test_video_ = GET_PARAM(2);
    target_bitrate_ = GET_PARAM(1);
    if (test_video_.width == 1280 && target_bitrate_ == 200) return;
    if (test_video_.width == 640 && target_bitrate_ == 1000) return;
    SetConfig();
    rc_api_ = libvpx::VP8RateControlRTC::Create(rc_cfg_);
    rc_api_->UpdateRateControl(rc_cfg_);

    ::libvpx_test::I420VideoSource video(test_video_.name, test_video_.width,
                                         test_video_.height, 30, 1, 0,
                                         test_video_.frames);

    ASSERT_NO_FATAL_FAILURE(RunLoop(&video));
  }

  void RunPeriodicKey() {
    test_video_ = GET_PARAM(2);
    target_bitrate_ = GET_PARAM(1);
    if (test_video_.width == 1280 && target_bitrate_ == 200) return;
    if (test_video_.width == 640 && target_bitrate_ == 1000) return;
    key_interval_ = 100;
    SetConfig();
    rc_api_ = libvpx::VP8RateControlRTC::Create(rc_cfg_);
    rc_api_->UpdateRateControl(rc_cfg_);

    ::libvpx_test::I420VideoSource video(test_video_.name, test_video_.width,
                                         test_video_.height, 30, 1, 0,
                                         test_video_.frames);

    ASSERT_NO_FATAL_FAILURE(RunLoop(&video));
  }

 private:
  void SetConfig() {
    rc_cfg_.width = test_video_.width;
    rc_cfg_.height = test_video_.height;
    rc_cfg_.max_quantizer = 60;
    rc_cfg_.min_quantizer = 2;
    rc_cfg_.target_bandwidth = target_bitrate_;
    rc_cfg_.buf_initial_sz = 600;
    rc_cfg_.buf_optimal_sz = 600;
    rc_cfg_.buf_sz = target_bitrate_;
    rc_cfg_.undershoot_pct = 50;
    rc_cfg_.overshoot_pct = 50;
    rc_cfg_.max_intra_bitrate_pct = 1000;
    rc_cfg_.framerate = 30.0;
    rc_cfg_.layer_target_bitrate[0] = target_bitrate_;

    // Encoder settings for ground truth.
    cfg_.g_w = test_video_.width;
    cfg_.g_h = test_video_.height;
    cfg_.rc_undershoot_pct = 50;
    cfg_.rc_overshoot_pct = 50;
    cfg_.rc_buf_initial_sz = 600;
    cfg_.rc_buf_optimal_sz = 600;
    cfg_.rc_buf_sz = target_bitrate_;
    cfg_.rc_dropframe_thresh = 0;
    cfg_.rc_min_quantizer = 2;
    cfg_.rc_max_quantizer = 60;
    cfg_.rc_end_usage = VPX_CBR;
    cfg_.g_lag_in_frames = 0;
    cfg_.g_error_resilient = 1;
    cfg_.rc_target_bitrate = target_bitrate_;
    cfg_.kf_min_dist = key_interval_;
    cfg_.kf_max_dist = key_interval_;
  }

  std::unique_ptr<libvpx::VP8RateControlRTC> rc_api_;
  libvpx::VP8RateControlRtcConfig rc_cfg_;
  int key_interval_;
  int target_bitrate_;
  Vp8RCTestVideo test_video_;
  libvpx::VP8FrameParamsQpRTC frame_params_;
  bool encoder_exit_;
};

TEST_P(Vp8RcInterfaceTest, OneLayer) { RunOneLayer(); }

TEST_P(Vp8RcInterfaceTest, OneLayerPeriodicKey) { RunPeriodicKey(); }

VP8_INSTANTIATE_TEST_SUITE(Vp8RcInterfaceTest,
                           ::testing::Values(200, 400, 1000),
                           ::testing::ValuesIn(kVp8RCTestVectors));

}  // namespace
