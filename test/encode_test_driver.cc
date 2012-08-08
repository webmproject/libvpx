/*
 *  Copyright (c) 2012 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */
#include "test/encode_test_driver.h"
#include "test/video_source.h"
#include "third_party/googletest/src/include/gtest/gtest.h"

namespace libvpx_test {

void Encoder::EncodeFrame(VideoSource *video, unsigned long flags) {
  if (video->img())
    EncodeFrameInternal(*video, flags);
  else
    Flush();

  // Handle twopass stats
  CxDataIterator iter = GetCxData();

  while (const vpx_codec_cx_pkt_t *pkt = iter.Next()) {
    if (pkt->kind != VPX_CODEC_STATS_PKT)
      continue;

    stats_->Append(*pkt);
  }
}

void Encoder::EncodeFrameInternal(const VideoSource &video,
                                  unsigned long flags) {
  vpx_codec_err_t res;
  const vpx_image_t *img = video.img();

  // Handle first frame initialization
  if (!encoder_.priv) {
    cfg_.g_w = img->d_w;
    cfg_.g_h = img->d_h;
    cfg_.g_timebase = video.timebase();
    cfg_.rc_twopass_stats_in = stats_->buf();
    res = vpx_codec_enc_init(&encoder_, &vpx_codec_vp8_cx_algo, &cfg_, 0);
    ASSERT_EQ(VPX_CODEC_OK, res) << EncoderError();
  }

  // Handle frame resizing
  if (cfg_.g_w != img->d_w || cfg_.g_h != img->d_h) {
    cfg_.g_w = img->d_w;
    cfg_.g_h = img->d_h;
    res = vpx_codec_enc_config_set(&encoder_, &cfg_);
    ASSERT_EQ(VPX_CODEC_OK, res) << EncoderError();
  }

  // Encode the frame
  res = vpx_codec_encode(&encoder_,
                         video.img(), video.pts(), video.duration(),
                         flags, deadline_);
  ASSERT_EQ(VPX_CODEC_OK, res) << EncoderError();
}

void Encoder::Flush() {
  const vpx_codec_err_t res = vpx_codec_encode(&encoder_, NULL, 0, 0, 0,
                                               deadline_);
  ASSERT_EQ(VPX_CODEC_OK, res) << EncoderError();
}

void EncoderTest::SetMode(TestMode mode) {
  switch (mode) {
    case kRealTime:
      deadline_ = VPX_DL_REALTIME;
      break;

    case kOnePassGood:
    case kTwoPassGood:
      deadline_ = VPX_DL_GOOD_QUALITY;
      break;

    case kOnePassBest:
    case kTwoPassBest:
      deadline_ = VPX_DL_BEST_QUALITY;
      break;

    default:
      ASSERT_TRUE(false) << "Unexpected mode " << mode;
  }

  if (mode == kTwoPassGood || mode == kTwoPassBest)
    passes_ = 2;
  else
    passes_ = 1;
}

void EncoderTest::RunLoop(VideoSource *video) {
  for (unsigned int pass = 0; pass < passes_; pass++) {
    last_pts_ = 0;

    if (passes_ == 1)
      cfg_.g_pass = VPX_RC_ONE_PASS;
    else if (pass == 0)
      cfg_.g_pass = VPX_RC_FIRST_PASS;
    else
      cfg_.g_pass = VPX_RC_LAST_PASS;

    BeginPassHook(pass);
    Encoder encoder(cfg_, deadline_, &stats_);

    bool again;
    for (again = true, video->Begin(); again; video->Next()) {
      again = video->img() != NULL;

      PreEncodeFrameHook(video);
      PreEncodeFrameHook(video, &encoder);
      encoder.EncodeFrame(video, flags_);

      CxDataIterator iter = encoder.GetCxData();

      while (const vpx_codec_cx_pkt_t *pkt = iter.Next()) {
        again = true;

        if (pkt->kind != VPX_CODEC_CX_FRAME_PKT)
          continue;

        ASSERT_GE(pkt->data.frame.pts, last_pts_);
        last_pts_ = pkt->data.frame.pts;
        FramePktHook(pkt);
      }

      if (!Continue())
        break;
    }

    EndPassHook();

    if (!Continue())
      break;
  }
}

}  // namespace libvpx_test
