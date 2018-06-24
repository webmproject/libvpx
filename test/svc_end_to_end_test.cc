/*
 *  Copyright (c) 2018 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */
#include "./vpx_config.h"
#include "third_party/googletest/src/include/gtest/gtest.h"
#include "test/codec_factory.h"
#include "test/encode_test_driver.h"
#include "test/i420_video_source.h"
#include "test/svc_test.h"
#include "test/util.h"
#include "test/y4m_video_source.h"
#include "vpx/vpx_codec.h"
#include "vpx_ports/bitops.h"

namespace {

// Params: Inter layer prediction modes.
class SyncFrameOnePassCbrSvc : public ::svc_test::OnePassCbrSvc,
                               public ::libvpx_test::CodecTestWithParam<int> {
 public:
  SyncFrameOnePassCbrSvc()
      : OnePassCbrSvc(GET_PARAM(0)), current_video_frame_(0),
        frame_to_start_decode_(0), mismatch_nframes_(0), num_nonref_frames_(0),
        inter_layer_pred_mode_(GET_PARAM(1)) {
    SetMode(::libvpx_test::kRealTime);
  }

 protected:
  virtual ~SyncFrameOnePassCbrSvc() {}

  virtual void SetUp() {
    InitializeConfig();
    speed_setting_ = 7;
  }

  virtual bool DoDecode() const {
    return current_video_frame_ >= frame_to_start_decode_;
  }

  virtual void PreEncodeFrameHook(::libvpx_test::VideoSource *video,
                                  ::libvpx_test::Encoder *encoder) {
    current_video_frame_ = video->frame();
    PreEncodeFrameHookSetup(video, encoder);
    if (video->frame() == 0)
      encoder->Control(VP9E_SET_SVC_INTER_LAYER_PRED, inter_layer_pred_mode_);
    if (video->frame() == frame_to_start_decode_) {
      vpx_svc_spatial_layer_sync_t svc_layer_sync;
      for (size_t i = 0; i < cfg_.ss_number_layers; i++)
        svc_layer_sync.spatial_layer_sync[i] = 0;
      svc_layer_sync.base_layer_intra_only = 0;
      svc_layer_sync.spatial_layer_sync[0] = 1;
      encoder->Control(VP9E_SET_SVC_SPATIAL_LAYER_SYNC, &svc_layer_sync);
    }
  }

  virtual void FramePktHook(const vpx_codec_cx_pkt_t *pkt) {
    // Keep track of number of non-reference frames, needed for mismatch check.
    // Non-reference frames are top spatial and temporal layer frames,
    // for TL > 0.
    if (temporal_layer_id_ == number_temporal_layers_ - 1 &&
        temporal_layer_id_ > 0 &&
        pkt->data.frame.spatial_layer_encoded[number_spatial_layers_ - 1] &&
        current_video_frame_ >= frame_to_start_decode_)
      num_nonref_frames_++;
  }

  virtual void MismatchHook(const vpx_image_t * /*img1*/,
                            const vpx_image_t * /*img2*/) {
    ++mismatch_nframes_;
  }

  unsigned int GetMismatchFrames() const { return mismatch_nframes_; }

  unsigned int current_video_frame_;
  unsigned int frame_to_start_decode_;
  unsigned int mismatch_nframes_;
  unsigned int num_nonref_frames_;
  int inter_layer_pred_mode_;
};

// Test for sync layer for 1 pass CBR SVC: 2 spatial layers and
// 3 temporal layers. Only start decoding on the sync layer.
TEST_P(SyncFrameOnePassCbrSvc, OnePassCbrSvc2SL3TLSyncFrame) {
  cfg_.rc_buf_initial_sz = 500;
  cfg_.rc_buf_optimal_sz = 500;
  cfg_.rc_buf_sz = 1000;
  cfg_.rc_min_quantizer = 0;
  cfg_.rc_max_quantizer = 63;
  cfg_.rc_end_usage = VPX_CBR;
  cfg_.g_lag_in_frames = 0;
  cfg_.ss_number_layers = 2;
  cfg_.ts_number_layers = 3;
  cfg_.ts_rate_decimator[0] = 4;
  cfg_.ts_rate_decimator[1] = 2;
  cfg_.ts_rate_decimator[2] = 1;
  cfg_.g_error_resilient = 1;
  cfg_.g_threads = 1;
  cfg_.temporal_layering_mode = 3;
  svc_params_.scaling_factor_num[0] = 144;
  svc_params_.scaling_factor_den[0] = 288;
  svc_params_.scaling_factor_num[1] = 288;
  svc_params_.scaling_factor_den[1] = 288;
  cfg_.rc_dropframe_thresh = 30;
  cfg_.kf_max_dist = 9999;
  number_spatial_layers_ = cfg_.ss_number_layers;
  number_temporal_layers_ = cfg_.ts_number_layers;
  frame_to_start_decode_ = 100;
  ::libvpx_test::I420VideoSource video("niklas_640_480_30.yuv", 640, 480, 30, 1,
                                       0, 400);
  cfg_.rc_target_bitrate = 400;
  AssignLayerBitrates();
  ASSERT_NO_FATAL_FAILURE(RunLoop(&video));
#if CONFIG_VP9_DECODER
  // The non-reference frames are expected to be mismatched frames as the
  // encoder will avoid loopfilter on these frames.
  EXPECT_EQ(num_nonref_frames_, GetMismatchFrames());
#endif
}

VP9_INSTANTIATE_TEST_CASE(SyncFrameOnePassCbrSvc, ::testing::Range(0, 3));

}  // namespace
