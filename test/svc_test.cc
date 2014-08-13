/*
 *  Copyright (c) 2013 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include <string>
#include "third_party/googletest/src/include/gtest/gtest.h"
#include "test/codec_factory.h"
#include "test/decode_test_driver.h"
#include "test/i420_video_source.h"

#include "vp9/decoder/vp9_decoder.h"

#include "vpx/svc_context.h"
#include "vpx/vp8cx.h"
#include "vpx/vpx_encoder.h"

namespace {

using libvpx_test::CodecFactory;
using libvpx_test::Decoder;
using libvpx_test::DxDataIterator;
using libvpx_test::VP9CodecFactory;

class SvcTest : public ::testing::Test {
 protected:
  static const uint32_t kWidth = 352;
  static const uint32_t kHeight = 288;

  SvcTest()
      : codec_iface_(0),
        test_file_name_("hantro_collage_w352h288.yuv"),
        codec_initialized_(false),
        decoder_(0) {
    memset(&svc_, 0, sizeof(svc_));
    memset(&codec_, 0, sizeof(codec_));
    memset(&codec_enc_, 0, sizeof(codec_enc_));
  }

  virtual ~SvcTest() {}

  virtual void SetUp() {
    svc_.log_level = SVC_LOG_DEBUG;
    svc_.log_print = 0;

    codec_iface_ = vpx_codec_vp9_cx();
    const vpx_codec_err_t res =
        vpx_codec_enc_config_default(codec_iface_, &codec_enc_, 0);
    EXPECT_EQ(VPX_CODEC_OK, res);

    codec_enc_.g_w = kWidth;
    codec_enc_.g_h = kHeight;
    codec_enc_.g_timebase.num = 1;
    codec_enc_.g_timebase.den = 60;
    codec_enc_.kf_min_dist = 100;
    codec_enc_.kf_max_dist = 100;

    vpx_codec_dec_cfg_t dec_cfg = {0};
    VP9CodecFactory codec_factory;
    decoder_ = codec_factory.CreateDecoder(dec_cfg, 0);
  }

  virtual void TearDown() {
    ReleaseEncoder();
    delete(decoder_);
  }

  void InitializeEncoder() {
    const vpx_codec_err_t res =
        vpx_svc_init(&svc_, &codec_, vpx_codec_vp9_cx(), &codec_enc_);
    EXPECT_EQ(VPX_CODEC_OK, res);
    codec_initialized_ = true;
  }

  void ReleaseEncoder() {
    vpx_svc_release(&svc_);
    if (codec_initialized_) vpx_codec_destroy(&codec_);
    codec_initialized_ = false;
  }

  void Pass1EncodeNFrames(const int n, const int layers,
                          std::string *const stats_buf) {
    vpx_codec_err_t res;
    size_t stats_size = 0;
    const char *stats_data = NULL;

    ASSERT_GT(n, 0);
    ASSERT_GT(layers, 0);
    svc_.spatial_layers = layers;
    codec_enc_.g_pass = VPX_RC_FIRST_PASS;
    InitializeEncoder();

    libvpx_test::I420VideoSource video(test_file_name_, kWidth, kHeight,
                                       codec_enc_.g_timebase.den,
                                       codec_enc_.g_timebase.num, 0, 30);
    video.Begin();

    for (int i = 0; i < n; ++i) {
      res = vpx_svc_encode(&svc_, &codec_, video.img(), video.pts(),
                           video.duration(), VPX_DL_GOOD_QUALITY);
      ASSERT_EQ(VPX_CODEC_OK, res);
      stats_size = vpx_svc_get_rc_stats_buffer_size(&svc_);
      EXPECT_GT(stats_size, 0U);
      stats_data = vpx_svc_get_rc_stats_buffer(&svc_);
      ASSERT_TRUE(stats_data != NULL);
      stats_buf->append(stats_data, stats_size);
      video.Next();
    }

    // Flush encoder and test EOS packet
    res = vpx_svc_encode(&svc_, &codec_, NULL, video.pts(),
                         video.duration(), VPX_DL_GOOD_QUALITY);
    stats_size = vpx_svc_get_rc_stats_buffer_size(&svc_);
    EXPECT_GT(stats_size, 0U);
    stats_data = vpx_svc_get_rc_stats_buffer(&svc_);
    ASSERT_TRUE(stats_data != NULL);
    stats_buf->append(stats_data, stats_size);

    ReleaseEncoder();
  }

  void StoreFrames(const size_t max_frame_received,
                   struct vpx_fixed_buf *const outputs,
                   size_t *const frame_received) {
    size_t frame_size;
    while ((frame_size = vpx_svc_get_frame_size(&svc_)) > 0) {
      ASSERT_LT(*frame_received, max_frame_received);

      if (*frame_received == 0) {
        EXPECT_EQ(1, vpx_svc_is_keyframe(&svc_));
      }

      outputs[*frame_received].buf = malloc(frame_size);
      ASSERT_TRUE(outputs[*frame_received].buf != NULL);
      memcpy(outputs[*frame_received].buf, vpx_svc_get_buffer(&svc_),
             frame_size);
      outputs[*frame_received].sz = frame_size;
      ++(*frame_received);
    }
  }

  void Pass2EncodeNFrames(std::string *const stats_buf,
                          const int n, const int layers,
                          struct vpx_fixed_buf *const outputs) {
    vpx_codec_err_t res;
    size_t frame_received = 0;

    ASSERT_TRUE(outputs != NULL);
    ASSERT_GT(n, 0);
    ASSERT_GT(layers, 0);
    svc_.spatial_layers = layers;
    codec_enc_.rc_target_bitrate = 500;
    if (codec_enc_.g_pass == VPX_RC_LAST_PASS) {
      ASSERT_TRUE(stats_buf != NULL);
      ASSERT_GT(stats_buf->size(), 0U);
      codec_enc_.rc_twopass_stats_in.buf = &(*stats_buf)[0];
      codec_enc_.rc_twopass_stats_in.sz = stats_buf->size();
    }
    InitializeEncoder();

    libvpx_test::I420VideoSource video(test_file_name_, kWidth, kHeight,
                                       codec_enc_.g_timebase.den,
                                       codec_enc_.g_timebase.num, 0, 30);
    video.Begin();

    for (int i = 0; i < n; ++i) {
      res = vpx_svc_encode(&svc_, &codec_, video.img(), video.pts(),
                           video.duration(), VPX_DL_GOOD_QUALITY);
      ASSERT_EQ(VPX_CODEC_OK, res);
      StoreFrames(n, outputs, &frame_received);
      video.Next();
    }

    // Flush Encoder
    res = vpx_svc_encode(&svc_, &codec_, NULL, 0,
                         video.duration(), VPX_DL_GOOD_QUALITY);
    EXPECT_EQ(VPX_CODEC_OK, res);
    StoreFrames(n, outputs, &frame_received);

    EXPECT_EQ(frame_received, (size_t)n);

    ReleaseEncoder();
  }

  void DecodeNFrames(const struct vpx_fixed_buf *const inputs, const int n) {
    int decoded_frames = 0;
    int received_frames = 0;

    ASSERT_TRUE(inputs != NULL);
    ASSERT_GT(n, 0);

    for (int i = 0; i < n; ++i) {
      ASSERT_TRUE(inputs[i].buf != NULL);
      ASSERT_GT(inputs[i].sz, 0U);
      const vpx_codec_err_t res_dec =
          decoder_->DecodeFrame(static_cast<const uint8_t *>(inputs[i].buf),
                                inputs[i].sz);
      ASSERT_EQ(VPX_CODEC_OK, res_dec) << decoder_->DecodeError();
      ++decoded_frames;

      DxDataIterator dec_iter = decoder_->GetDxData();
      while (dec_iter.Next()) {
        ++received_frames;
      }
    }
    EXPECT_EQ(decoded_frames, n);
    EXPECT_EQ(received_frames, n);
  }

  void DropEnhancementLayers(struct vpx_fixed_buf *const inputs,
                             const int num_super_frames,
                             const int remained_layers) {
    ASSERT_TRUE(inputs != NULL);
    ASSERT_GT(num_super_frames, 0);
    ASSERT_GT(remained_layers, 0);

    for (int i = 0; i < num_super_frames; ++i) {
      uint32_t frame_sizes[8] = {0};
      int frame_count = 0;
      int frames_found = 0;
      int frame;
      ASSERT_TRUE(inputs[i].buf != NULL);
      ASSERT_GT(inputs[i].sz, 0U);

      vpx_codec_err_t res =
          vp9_parse_superframe_index(static_cast<const uint8_t*>(inputs[i].buf),
                                     inputs[i].sz, frame_sizes, &frame_count,
                                     NULL, NULL);
      ASSERT_EQ(VPX_CODEC_OK, res);

      uint8_t *frame_data = static_cast<uint8_t *>(inputs[i].buf);
      uint8_t *frame_start = frame_data;
      for (frame = 0; frame < frame_count; ++frame) {
        // Looking for a visible frame
        if (frame_data[0] & 0x02) {
          ++frames_found;
          if (frames_found == remained_layers)
            break;
        }
        frame_data += frame_sizes[frame];
      }
      ASSERT_LT(frame, frame_count);
      if (frame == frame_count - 1)
        continue;

      frame_data += frame_sizes[frame];
      uint8_t marker =
          static_cast<const uint8_t *>(inputs[i].buf)[inputs[i].sz - 1];
      const uint32_t mag = ((marker >> 3) & 0x3) + 1;
      const size_t index_sz = 2 + mag * frame_count;
      const size_t new_index_sz = 2 + mag * (frame + 1);
      marker &= 0x0f8;
      marker |= frame;
      frame_data[0] = marker;
      memcpy(frame_data + 1, frame_start + inputs[i].sz - index_sz + 1,
             new_index_sz - 2);
      frame_data[new_index_sz - 1] = marker;
      inputs[i].sz = frame_data - frame_start + new_index_sz;
    }
  }

  void FreeBitstreamBuffers(struct vpx_fixed_buf *const inputs, const int n) {
    ASSERT_TRUE(inputs != NULL);
    ASSERT_GT(n, 0);

    for (int i = 0; i < n; ++i) {
      free(inputs[i].buf);
      inputs[i].buf = NULL;
      inputs[i].sz = 0;
    }
  }

  SvcContext svc_;
  vpx_codec_ctx_t codec_;
  struct vpx_codec_enc_cfg codec_enc_;
  vpx_codec_iface_t *codec_iface_;
  std::string test_file_name_;
  bool codec_initialized_;
  Decoder *decoder_;
};

TEST_F(SvcTest, SvcInit) {
  // test missing parameters
  vpx_codec_err_t res = vpx_svc_init(NULL, &codec_, codec_iface_, &codec_enc_);
  EXPECT_EQ(VPX_CODEC_INVALID_PARAM, res);
  res = vpx_svc_init(&svc_, NULL, codec_iface_, &codec_enc_);
  EXPECT_EQ(VPX_CODEC_INVALID_PARAM, res);
  res = vpx_svc_init(&svc_, &codec_, NULL, &codec_enc_);
  EXPECT_EQ(VPX_CODEC_INVALID_PARAM, res);

  res = vpx_svc_init(&svc_, &codec_, codec_iface_, NULL);
  EXPECT_EQ(VPX_CODEC_INVALID_PARAM, res);

  svc_.spatial_layers = 6;  // too many layers
  res = vpx_svc_init(&svc_, &codec_, codec_iface_, &codec_enc_);
  EXPECT_EQ(VPX_CODEC_INVALID_PARAM, res);

  svc_.spatial_layers = 0;  // use default layers
  InitializeEncoder();
  EXPECT_EQ(VPX_SS_DEFAULT_LAYERS, svc_.spatial_layers);
}

TEST_F(SvcTest, InitTwoLayers) {
  svc_.spatial_layers = 2;
  vpx_svc_set_scale_factors(&svc_, "4/16,16*16");  // invalid scale values
  vpx_codec_err_t res = vpx_svc_init(&svc_, &codec_, codec_iface_, &codec_enc_);
  EXPECT_EQ(VPX_CODEC_INVALID_PARAM, res);

  vpx_svc_set_scale_factors(&svc_, "4/16,16/16");  // valid scale values
  InitializeEncoder();
}

TEST_F(SvcTest, InvalidOptions) {
  vpx_codec_err_t res = vpx_svc_set_options(&svc_, NULL);
  EXPECT_EQ(VPX_CODEC_INVALID_PARAM, res);

  res = vpx_svc_set_options(&svc_, "not-an-option=1");
  EXPECT_EQ(VPX_CODEC_OK, res);
  res = vpx_svc_init(&svc_, &codec_, vpx_codec_vp9_cx(), &codec_enc_);
  EXPECT_EQ(VPX_CODEC_INVALID_PARAM, res);
}

TEST_F(SvcTest, SetLayersOption) {
  vpx_codec_err_t res = vpx_svc_set_options(&svc_, "layers=3");
  EXPECT_EQ(VPX_CODEC_OK, res);
  InitializeEncoder();
  EXPECT_EQ(3, svc_.spatial_layers);
}

TEST_F(SvcTest, SetMultipleOptions) {
  vpx_codec_err_t res =
      vpx_svc_set_options(&svc_, "layers=2 scale-factors=1/3,2/3");
  EXPECT_EQ(VPX_CODEC_OK, res);
  InitializeEncoder();
  EXPECT_EQ(2, svc_.spatial_layers);
}

TEST_F(SvcTest, SetScaleFactorsOption) {
  svc_.spatial_layers = 2;
  vpx_codec_err_t res =
      vpx_svc_set_options(&svc_, "scale-factors=not-scale-factors");
  EXPECT_EQ(VPX_CODEC_OK, res);
  res = vpx_svc_init(&svc_, &codec_, vpx_codec_vp9_cx(), &codec_enc_);
  EXPECT_EQ(VPX_CODEC_INVALID_PARAM, res);

  res = vpx_svc_set_options(&svc_, "scale-factors=1/3,2/3");
  EXPECT_EQ(VPX_CODEC_OK, res);
  InitializeEncoder();
}

TEST_F(SvcTest, SetQuantizersOption) {
  svc_.spatial_layers = 2;
  vpx_codec_err_t res = vpx_svc_set_options(&svc_, "quantizers=not-quantizers");
  EXPECT_EQ(VPX_CODEC_OK, res);
  res = vpx_svc_init(&svc_, &codec_, vpx_codec_vp9_cx(), &codec_enc_);
  EXPECT_EQ(VPX_CODEC_INVALID_PARAM, res);

  vpx_svc_set_options(&svc_, "quantizers=40,45");
  InitializeEncoder();
}

TEST_F(SvcTest, SetAutoAltRefOption) {
  svc_.spatial_layers = 5;
  vpx_codec_err_t res = vpx_svc_set_options(&svc_, "auto-alt-refs=none");
  EXPECT_EQ(VPX_CODEC_OK, res);
  res = vpx_svc_init(&svc_, &codec_, vpx_codec_vp9_cx(), &codec_enc_);
  EXPECT_EQ(VPX_CODEC_INVALID_PARAM, res);

  res = vpx_svc_set_options(&svc_, "auto-alt-refs=1,1,1,1,0");
  EXPECT_EQ(VPX_CODEC_OK, res);
  res = vpx_svc_init(&svc_, &codec_, vpx_codec_vp9_cx(), &codec_enc_);
  EXPECT_EQ(VPX_CODEC_INVALID_PARAM, res);

  vpx_svc_set_options(&svc_, "auto-alt-refs=0,1,1,1,0");
  InitializeEncoder();
}

TEST_F(SvcTest, SetQuantizers) {
  vpx_codec_err_t res = vpx_svc_set_quantizers(NULL, "40,30");
  EXPECT_EQ(VPX_CODEC_INVALID_PARAM, res);

  res = vpx_svc_set_quantizers(&svc_, NULL);
  EXPECT_EQ(VPX_CODEC_INVALID_PARAM, res);

  svc_.spatial_layers = 2;
  res = vpx_svc_set_quantizers(&svc_, "40");
  EXPECT_EQ(VPX_CODEC_OK, res);
  res = vpx_svc_init(&svc_, &codec_, vpx_codec_vp9_cx(), &codec_enc_);
  EXPECT_EQ(VPX_CODEC_INVALID_PARAM, res);

  res = vpx_svc_set_quantizers(&svc_, "40,30");
  EXPECT_EQ(VPX_CODEC_OK, res);
  InitializeEncoder();
}

TEST_F(SvcTest, SetScaleFactors) {
  vpx_codec_err_t res = vpx_svc_set_scale_factors(NULL, "4/16,16/16");
  EXPECT_EQ(VPX_CODEC_INVALID_PARAM, res);

  res = vpx_svc_set_scale_factors(&svc_, NULL);
  EXPECT_EQ(VPX_CODEC_INVALID_PARAM, res);

  svc_.spatial_layers = 2;
  res = vpx_svc_set_scale_factors(&svc_, "4/16");
  EXPECT_EQ(VPX_CODEC_OK, res);
  res = vpx_svc_init(&svc_, &codec_, vpx_codec_vp9_cx(), &codec_enc_);
  EXPECT_EQ(VPX_CODEC_INVALID_PARAM, res);

  res = vpx_svc_set_scale_factors(&svc_, "4/16,16/16");
  EXPECT_EQ(VPX_CODEC_OK, res);
  InitializeEncoder();
}

// Test that decoder can handle an SVC frame as the first frame in a sequence.
TEST_F(SvcTest, OnePassEncodeOneFrame) {
  codec_enc_.g_pass = VPX_RC_ONE_PASS;
  vpx_fixed_buf output = {0};
  Pass2EncodeNFrames(NULL, 1, 2, &output);
  DecodeNFrames(&output, 1);
  FreeBitstreamBuffers(&output, 1);
}

TEST_F(SvcTest, OnePassEncodeThreeFrames) {
  codec_enc_.g_pass = VPX_RC_ONE_PASS;
  vpx_fixed_buf outputs[3];
  memset(&outputs[0], 0, sizeof(outputs));
  Pass2EncodeNFrames(NULL, 3, 2, &outputs[0]);
  DecodeNFrames(&outputs[0], 3);
  FreeBitstreamBuffers(&outputs[0], 3);
}

TEST_F(SvcTest, GetLayerResolution) {
  svc_.spatial_layers = 2;
  vpx_svc_set_scale_factors(&svc_, "4/16,8/16");
  vpx_svc_set_quantizers(&svc_, "40,30");

  InitializeEncoder();

  // ensure that requested layer is a valid layer
  uint32_t layer_width, layer_height;
  vpx_codec_err_t res = vpx_svc_get_layer_resolution(&svc_, svc_.spatial_layers,
                                     &layer_width, &layer_height);
  EXPECT_EQ(VPX_CODEC_INVALID_PARAM, res);

  res = vpx_svc_get_layer_resolution(NULL, 0, &layer_width, &layer_height);
  EXPECT_EQ(VPX_CODEC_INVALID_PARAM, res);

  res = vpx_svc_get_layer_resolution(&svc_, 0, NULL, &layer_height);
  EXPECT_EQ(VPX_CODEC_INVALID_PARAM, res);

  res = vpx_svc_get_layer_resolution(&svc_, 0, &layer_width, NULL);
  EXPECT_EQ(VPX_CODEC_INVALID_PARAM, res);

  res = vpx_svc_get_layer_resolution(&svc_, 0, &layer_width, &layer_height);
  EXPECT_EQ(VPX_CODEC_OK, res);
  EXPECT_EQ(kWidth * 4 / 16, layer_width);
  EXPECT_EQ(kHeight * 4 / 16, layer_height);

  res = vpx_svc_get_layer_resolution(&svc_, 1, &layer_width, &layer_height);
  EXPECT_EQ(VPX_CODEC_OK, res);
  EXPECT_EQ(kWidth * 8 / 16, layer_width);
  EXPECT_EQ(kHeight * 8 / 16, layer_height);
}

TEST_F(SvcTest, TwoPassEncode10Frames) {
  // First pass encode
  std::string stats_buf;
  Pass1EncodeNFrames(10, 2, &stats_buf);

  // Second pass encode
  codec_enc_.g_pass = VPX_RC_LAST_PASS;
  vpx_fixed_buf outputs[10];
  memset(&outputs[0], 0, sizeof(outputs));
  Pass2EncodeNFrames(&stats_buf, 10, 2, &outputs[0]);
  DecodeNFrames(&outputs[0], 10);
  FreeBitstreamBuffers(&outputs[0], 10);
}

TEST_F(SvcTest, TwoPassEncode20FramesWithAltRef) {
  // First pass encode
  std::string stats_buf;
  Pass1EncodeNFrames(20, 2, &stats_buf);

  // Second pass encode
  codec_enc_.g_pass = VPX_RC_LAST_PASS;
  vpx_svc_set_options(&svc_, "auto-alt-refs=1,1");
  vpx_fixed_buf outputs[20];
  memset(&outputs[0], 0, sizeof(outputs));
  Pass2EncodeNFrames(&stats_buf, 20, 2, &outputs[0]);
  DecodeNFrames(&outputs[0], 20);
  FreeBitstreamBuffers(&outputs[0], 20);
}

TEST_F(SvcTest, TwoPassEncode2LayersDecodeBaseLayerOnly) {
  // First pass encode
  std::string stats_buf;
  Pass1EncodeNFrames(10, 2, &stats_buf);

  // Second pass encode
  codec_enc_.g_pass = VPX_RC_LAST_PASS;
  vpx_svc_set_options(&svc_, "auto-alt-refs=1,1");
  vpx_fixed_buf outputs[10];
  memset(&outputs[0], 0, sizeof(outputs));
  Pass2EncodeNFrames(&stats_buf, 10, 2, &outputs[0]);
  DropEnhancementLayers(&outputs[0], 10, 1);
  DecodeNFrames(&outputs[0], 10);
  FreeBitstreamBuffers(&outputs[0], 10);
}

TEST_F(SvcTest, TwoPassEncode5LayersDecode54321Layers) {
  // First pass encode
  std::string stats_buf;
  Pass1EncodeNFrames(10, 5, &stats_buf);

  // Second pass encode
  codec_enc_.g_pass = VPX_RC_LAST_PASS;
  vpx_svc_set_options(&svc_, "auto-alt-refs=0,1,1,1,0");
  vpx_fixed_buf outputs[10];
  memset(&outputs[0], 0, sizeof(outputs));
  Pass2EncodeNFrames(&stats_buf, 10, 5, &outputs[0]);

  DecodeNFrames(&outputs[0], 10);
  DropEnhancementLayers(&outputs[0], 10, 4);
  DecodeNFrames(&outputs[0], 10);
  DropEnhancementLayers(&outputs[0], 10, 3);
  DecodeNFrames(&outputs[0], 10);
  DropEnhancementLayers(&outputs[0], 10, 2);
  DecodeNFrames(&outputs[0], 10);
  DropEnhancementLayers(&outputs[0], 10, 1);
  DecodeNFrames(&outputs[0], 10);

  FreeBitstreamBuffers(&outputs[0], 10);
}

TEST_F(SvcTest, TwoPassEncode2SNRLayers) {
  // First pass encode
  std::string stats_buf;
  vpx_svc_set_options(&svc_, "scale-factors=1/1,1/1");
  Pass1EncodeNFrames(20, 2, &stats_buf);

  // Second pass encode
  codec_enc_.g_pass = VPX_RC_LAST_PASS;
  vpx_svc_set_options(&svc_,
                      "auto-alt-refs=1,1 scale-factors=1/1,1/1");
  vpx_fixed_buf outputs[20];
  memset(&outputs[0], 0, sizeof(outputs));
  Pass2EncodeNFrames(&stats_buf, 20, 2, &outputs[0]);
  DecodeNFrames(&outputs[0], 20);
  FreeBitstreamBuffers(&outputs[0], 20);
}

TEST_F(SvcTest, TwoPassEncode3SNRLayersDecode321Layers) {
  // First pass encode
  std::string stats_buf;
  vpx_svc_set_options(&svc_, "scale-factors=1/1,1/1,1/1");
  Pass1EncodeNFrames(20, 3, &stats_buf);

  // Second pass encode
  codec_enc_.g_pass = VPX_RC_LAST_PASS;
  vpx_svc_set_options(&svc_,
                      "auto-alt-refs=1,1,1 scale-factors=1/1,1/1,1/1");
  vpx_fixed_buf outputs[20];
  memset(&outputs[0], 0, sizeof(outputs));
  Pass2EncodeNFrames(&stats_buf, 20, 3, &outputs[0]);
  DecodeNFrames(&outputs[0], 20);
  DropEnhancementLayers(&outputs[0], 20, 2);
  DecodeNFrames(&outputs[0], 20);
  DropEnhancementLayers(&outputs[0], 20, 1);
  DecodeNFrames(&outputs[0], 20);

  FreeBitstreamBuffers(&outputs[0], 20);
}

}  // namespace
