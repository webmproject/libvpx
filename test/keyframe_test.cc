/*
 *  Copyright (c) 2012 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */
#include <climits>
#include <vector>
#include "test/encode_test_driver.h"
#include "test/video_source.h"
#include "third_party/googletest/src/include/gtest/gtest.h"

namespace {

class KeyframeTest : public ::libvpx_test::EncoderTest,
  public ::testing::TestWithParam<enum libvpx_test::TestMode> {
 protected:
  virtual void SetUp() {
    InitializeConfig();
    SetMode(GetParam());
    kf_count_ = 0;
    kf_count_max_ = INT_MAX;
    kf_do_force_kf_ = false;
  }

  virtual bool Continue() {
    return !HasFatalFailure() && !abort_;
  }

  virtual void PreEncodeFrameHook(::libvpx_test::VideoSource *video) {
    if (kf_do_force_kf_)
      flags_ = (video->frame() % 3) ? 0 : VPX_EFLAG_FORCE_KF;
  }

  virtual void FramePktHook(const vpx_codec_cx_pkt_t *pkt) {
    if (pkt->data.frame.flags & VPX_FRAME_IS_KEY) {
      kf_pts_list_.push_back(pkt->data.frame.pts);
      kf_count_++;
      abort_ |= kf_count_ > kf_count_max_;
    }
  }

  bool kf_do_force_kf_;
  int kf_count_;
  int kf_count_max_;
  std::vector< vpx_codec_pts_t > kf_pts_list_;
};

TEST_P(KeyframeTest, TestRandomVideoSource) {
  // Validate that encoding the RandomVideoSource produces multiple keyframes.
  // This validates the results of the TestDisableKeyframes test.
  kf_count_max_ = 2;  // early exit successful tests.

  ::libvpx_test::RandomVideoSource video;
  ASSERT_NO_FATAL_FAILURE(RunLoop(&video));

  EXPECT_GT(kf_count_, 1);
}

TEST_P(KeyframeTest, TestDisableKeyframes) {
  cfg_.kf_mode = VPX_KF_DISABLED;
  kf_count_max_ = 1;  // early exit failed tests.

  ::libvpx_test::RandomVideoSource video;
  ASSERT_NO_FATAL_FAILURE(RunLoop(&video));

  EXPECT_EQ(1, kf_count_);
}

TEST_P(KeyframeTest, TestForceKeyframe) {
  cfg_.kf_mode = VPX_KF_DISABLED;
  kf_do_force_kf_ = true;

  ::libvpx_test::DummyVideoSource video;
  ASSERT_NO_FATAL_FAILURE(RunLoop(&video));

  // verify that every third frame is a keyframe.
  for (std::vector<vpx_codec_pts_t>::iterator iter = kf_pts_list_.begin();
       iter != kf_pts_list_.end();
       ++iter) {
    ASSERT_EQ(0, *iter % 3) << "Unexpected keyframe at frame " << *iter;
  }
}

TEST_P(KeyframeTest, TestKeyframeMaxDistance) {
  cfg_.kf_max_dist = 25;

  ::libvpx_test::DummyVideoSource video;
  ASSERT_NO_FATAL_FAILURE(RunLoop(&video));

  // verify that keyframe interval matches kf_max_dist
  for (std::vector<vpx_codec_pts_t>::iterator iter = kf_pts_list_.begin();
       iter != kf_pts_list_.end();
       iter++) {
    ASSERT_EQ(0, *iter % 25) << "Unexpected keyframe at frame " << *iter;
  }
}

INSTANTIATE_TEST_CASE_P(AllModes, KeyframeTest, ALL_TEST_MODES);
}  // namespace
