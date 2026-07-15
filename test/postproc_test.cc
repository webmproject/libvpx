/*
 *  Copyright (c) 2026 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include <algorithm>
#include <array>
#include <memory>
#include <new>
#include <ostream>
#include <string>
#include <tuple>
#include <vector>

#include "gtest/gtest.h"

#include "test/codec_factory.h"
#include "test/decode_test_driver.h"
#include "test/ivf_video_source.h"
#include "test/test_vectors.h"
#include "test/util.h"
#if CONFIG_WEBM_IO
#include "test/webm_video_source.h"
#endif
#include "vpx/vp8dx.h"
#include "vpx/vpx_decoder.h"

namespace {

struct TestParam {
  int flags;
  const char *filename;
};

std::ostream &operator<<(std::ostream &os, const TestParam &param) {
  std::string flags_str;
  if (param.flags == VP8_NOFILTERING) {
    flags_str = "default flags";
  } else {
    if (param.flags & VP8_DEBLOCK) flags_str += "deblock";
    if (param.flags & VP8_DEMACROBLOCK) {
      if (!flags_str.empty()) flags_str += ", ";
      flags_str += "demacroblock";
    }
    if (param.flags & VP8_ADDNOISE) {
      if (!flags_str.empty()) flags_str += ", ";
      flags_str += "addnoise";
    }
    if (param.flags & VP8_MFQE) {
      if (!flags_str.empty()) flags_str += ", ";
      flags_str += "mfqe";
    }
  }

  return os << "(" << flags_str << ", " << param.filename << ")";
}

class PostProcTest : public ::libvpx_test::DecoderTest,
                     public ::libvpx_test::CodecTestWithParam<TestParam> {
 protected:
  PostProcTest()
      : DecoderTest(GET_PARAM(0)), postproc_flags_(GET_PARAM(1).flags) {}

  ~PostProcTest() override = default;

  void PreDecodeFrameHook(const libvpx_test::CompressedVideoSource &video,
                          libvpx_test::Decoder *decoder) override {
    // VP8_NOFILTERING is used to test the default flags set when only
    // VPX_CODEC_USE_POSTPROC is used.
    if (video.frame_number() == 0 && postproc_flags_ != VP8_NOFILTERING) {
      const vp8_postproc_cfg_t vp8_pp_cfg = { postproc_flags_,
                                              /*deblocking_level=*/2,
                                              /*noise_level=*/10 };
      decoder->Control(VP8_SET_POSTPROC, &vp8_pp_cfg);
    }
  }

  void DecodeTest();

 private:
  int postproc_flags_;
};

void PostProcTest::DecodeTest() {
  const std::string filename = GET_PARAM(1).filename;
  vpx_codec_dec_cfg_t cfg = vpx_codec_dec_cfg_t();
  cfg.threads = 2;

  // Open compressed video file.
  std::unique_ptr<libvpx_test::CompressedVideoSource> video;
  if (filename.substr(filename.length() - 3, 3) == "ivf") {
    video.reset(new (std::nothrow) libvpx_test::IVFVideoSource(filename));
  } else if (filename.substr(filename.length() - 4, 4) == "webm") {
#if CONFIG_WEBM_IO
    video.reset(new (std::nothrow) libvpx_test::WebMVideoSource(filename));
#else
    GTEST_SKIP() << "WebM IO is disabled, skipping test vector " << filename;
#endif
  }
  ASSERT_NE(video, nullptr);
  video->Init();

  // Set decode config and flags.
  set_cfg(cfg);
  set_flags(VPX_CODEC_USE_POSTPROC);

  // Decode frame.
  ASSERT_NO_FATAL_FAILURE(RunLoop(video.get(), cfg));
}

TEST_P(PostProcTest, Decode) { DecodeTest(); }

std::vector<int> GeneratePostProcFlags() {
  std::vector<int> flags;

  flags.push_back(VP8_NOFILTERING);
  // Add each flag (or flags) on their own, and with add noise enabled.
  static constexpr std::array<int, 7> kFlags = { VP8_DEBLOCK,
                                                 VP8_DEMACROBLOCK,
                                                 VP8_ADDNOISE,
                                                 VP8_MFQE,
                                                 VP8_MFQE | VP8_DEBLOCK,
                                                 VP8_MFQE | VP8_DEMACROBLOCK,
                                                 VP8_MFQE | VP8_DEBLOCK |
                                                     VP8_DEMACROBLOCK };
  for (const int flag : kFlags) {
    flags.push_back(flag);
    if (flag != VP8_ADDNOISE) {
      flags.push_back(flag | VP8_ADDNOISE);
    }
  }

  std::sort(flags.begin(), flags.end());
  return flags;
}

template <typename Iterator>
std::vector<TestParam> GenerateTestParams(const std::vector<int> &flags,
                                          Iterator begin, Iterator end) {
  std::vector<TestParam> params;
  for (const int flag : flags) {
    for (Iterator it = begin; it != end; ++it) {
      params.push_back({ flag, *it });
    }
  }
  return params;
}

template <size_t N>
std::vector<TestParam> GenerateTestParams(
    const std::vector<int> &flags,
    const std::array<const char *, N> &filenames) {
  return GenerateTestParams(flags, filenames.cbegin(), filenames.cend());
}

class PostProcTestInvalidFiles : public PostProcTest {
 protected:
  void HandlePeekResult(libvpx_test::Decoder *const /*decoder*/,
                        libvpx_test::CompressedVideoSource * /*video*/,
                        const vpx_codec_err_t /*res_peek*/) override {}

  bool HandleDecodeResult(const vpx_codec_err_t /*res_dec*/,
                          const libvpx_test::CompressedVideoSource & /*video*/,
                          libvpx_test::Decoder * /*decoder*/) override {
    // Ignore decode errors when decoding corrupted bitstreams.
    return true;
  }
};

TEST_P(PostProcTestInvalidFiles, Decode) { DecodeTest(); }

#if CONFIG_VP8_DECODER && CONFIG_POSTPROC
// TODO(issue 499602810): remove this test suite after valgrind errors are
// fixed, and the files are migrated to the enabled test suite.
constexpr std::array<const char *, 5> kVP8FailingTestVectors = {
  "invalid-bug-1443.ivf", "vp80-00-comprehensive-008.ivf",
  "vp80-02-inter-1418.ivf", "vp80-03-segmentation-1425.ivf",
  "vp80-03-segmentation-1436.ivf"
};

std::vector<TestParam> GenerateVP8PassingTestParams() {
  std::vector<TestParam> params;
  const std::vector<int> flags = GeneratePostProcFlags();
  for (const int flag : flags) {
    for (const char *const *it = libvpx_test::kVP8TestVectors;
         it != libvpx_test::kVP8TestVectors + libvpx_test::kNumVP8TestVectors;
         ++it) {
      if (std::find_if(kVP8FailingTestVectors.cbegin(),
                       kVP8FailingTestVectors.cend(), [it](const char *f) {
                         return std::string(f) == *it;
                       }) == kVP8FailingTestVectors.cend()) {
        params.push_back({ flag, *it });
      }
    }
  }
  return params;
}

VP8_INSTANTIATE_TEST_SUITE(PostProcTest,
                           ::testing::ValuesIn(GenerateVP8PassingTestParams()));

// TODO(issue 499602810): remove this test suite after asserts, out of
// bounds accesses, and valgrind errors are fixed, and the files are migrated
// to the enabled test suite.
INSTANTIATE_TEST_SUITE_P(
    DISABLED_VP8, PostProcTest,
    ::testing::Combine(
        ::testing::Values(
            static_cast<const libvpx_test::CodecFactory *>(&libvpx_test::kVP8)),
        ::testing::ValuesIn(GenerateTestParams(GeneratePostProcFlags(),
                                               kVP8FailingTestVectors))));

constexpr std::array<const char *, 4> kVP8InvalidFiles = {
  "invalid-bug-148271109.ivf", "invalid-token-partition.ivf",
  "invalid-vp80-00-comprehensive-018.ivf.2kf_0x6.ivf",
  "invalid-vp80-00-comprehensive-s17661_r01-05_b6-.ivf"
};

VP8_INSTANTIATE_TEST_SUITE(PostProcTestInvalidFiles,
                           ::testing::ValuesIn(GenerateTestParams(
                               GeneratePostProcFlags(), kVP8InvalidFiles)));
#endif  // CONFIG_VP8_DECODER && CONFIG_POSTPROC

#if CONFIG_VP9_DECODER && CONFIG_VP9_POSTPROC
constexpr std::array<const char *, 104> kVP9FailingTestVectors = {
  "invalid-vp90-02-v2.webm",
  "vp90-2-02-size-08x08.webm",
  "vp90-2-02-size-08x10.webm",
  "vp90-2-02-size-08x16.webm",
  "vp90-2-02-size-08x18.webm",
  "vp90-2-02-size-08x32.webm",
  "vp90-2-02-size-08x34.webm",
  "vp90-2-02-size-08x64.webm",
  "vp90-2-02-size-08x66.webm",
  "vp90-2-02-size-10x08.webm",
  "vp90-2-02-size-10x10.webm",
  "vp90-2-02-size-10x16.webm",
  "vp90-2-02-size-10x18.webm",
  "vp90-2-02-size-10x32.webm",
  "vp90-2-02-size-10x34.webm",
  "vp90-2-02-size-10x64.webm",
  "vp90-2-02-size-10x66.webm",
  "vp90-2-02-size-130x132.webm",
  "vp90-2-02-size-132x130.webm",
  "vp90-2-02-size-132x132.webm",
  "vp90-2-02-size-178x180.webm",
  "vp90-2-02-size-180x178.webm",
  "vp90-2-02-size-180x180.webm",
  "vp90-2-02-size-18x08.webm",
  "vp90-2-02-size-18x10.webm",
  "vp90-2-02-size-18x16.webm",
  "vp90-2-02-size-18x18.webm",
  "vp90-2-02-size-18x32.webm",
  "vp90-2-02-size-18x34.webm",
  "vp90-2-02-size-18x64.webm",
  "vp90-2-02-size-18x66.webm",
  "vp90-2-02-size-34x08.webm",
  "vp90-2-02-size-34x10.webm",
  "vp90-2-02-size-34x16.webm",
  "vp90-2-02-size-34x18.webm",
  "vp90-2-02-size-34x32.webm",
  "vp90-2-02-size-34x34.webm",
  "vp90-2-02-size-34x64.webm",
  "vp90-2-02-size-34x66.webm",
  "vp90-2-02-size-66x08.webm",
  "vp90-2-02-size-66x10.webm",
  "vp90-2-02-size-66x16.webm",
  "vp90-2-02-size-66x18.webm",
  "vp90-2-02-size-66x32.webm",
  "vp90-2-02-size-66x34.webm",
  "vp90-2-02-size-66x64.webm",
  "vp90-2-02-size-66x66.webm",
  "vp90-2-03-size-196x196.webm",
  "vp90-2-03-size-196x198.webm",
  "vp90-2-03-size-196x200.webm",
  "vp90-2-03-size-196x202.webm",
  "vp90-2-03-size-196x208.webm",
  "vp90-2-03-size-196x210.webm",
  "vp90-2-03-size-196x224.webm",
  "vp90-2-03-size-196x226.webm",
  "vp90-2-03-size-198x196.webm",
  "vp90-2-03-size-198x198.webm",
  "vp90-2-03-size-198x200.webm",
  "vp90-2-03-size-198x202.webm",
  "vp90-2-03-size-198x208.webm",
  "vp90-2-03-size-198x210.webm",
  "vp90-2-03-size-198x224.webm",
  "vp90-2-03-size-198x226.webm",
  "vp90-2-03-size-200x196.webm",
  "vp90-2-03-size-200x198.webm",
  "vp90-2-03-size-200x200.webm",
  "vp90-2-03-size-200x202.webm",
  "vp90-2-03-size-200x208.webm",
  "vp90-2-03-size-200x210.webm",
  "vp90-2-03-size-200x224.webm",
  "vp90-2-03-size-200x226.webm",
  "vp90-2-03-size-202x196.webm",
  "vp90-2-03-size-202x198.webm",
  "vp90-2-03-size-202x200.webm",
  "vp90-2-03-size-202x202.webm",
  "vp90-2-03-size-202x208.webm",
  "vp90-2-03-size-202x210.webm",
  "vp90-2-03-size-202x224.webm",
  "vp90-2-03-size-202x226.webm",
  "vp90-2-03-size-210x196.webm",
  "vp90-2-03-size-210x198.webm",
  "vp90-2-03-size-210x200.webm",
  "vp90-2-03-size-210x202.webm",
  "vp90-2-03-size-210x208.webm",
  "vp90-2-03-size-210x210.webm",
  "vp90-2-03-size-210x224.webm",
  "vp90-2-03-size-210x226.webm",
  "vp90-2-03-size-226x196.webm",
  "vp90-2-03-size-226x198.webm",
  "vp90-2-03-size-226x200.webm",
  "vp90-2-03-size-226x202.webm",
  "vp90-2-03-size-226x208.webm",
  "vp90-2-03-size-226x210.webm",
  "vp90-2-03-size-226x224.webm",
  "vp90-2-03-size-226x226.webm",
  "vp90-2-14-resize-10frames-fp-tiles-1-4.webm",
  "vp90-2-14-resize-10frames-fp-tiles-1-8.webm",
  "vp90-2-14-resize-10frames-fp-tiles-8-4-2-1.webm",
  "vp90-2-15-segkey_adpq.webm",
  "vp90-2-18-resize.ivf",
  "vp90-2-21-resize_inter_1920x1080_5_3-4.webm",
  "vp90-2-21-resize_inter_1920x1080_7_3-4.webm",
  "vp90-2-21-resize_inter_320x240_5_3-4.webm",
  "vp90-2-21-resize_inter_640x480_7_3-4.webm"
};

std::vector<TestParam> GenerateVP9PassingTestParams() {
  std::vector<TestParam> params;
  const std::vector<int> flags = GeneratePostProcFlags();
  for (const int flag : flags) {
    for (const char *const *it = libvpx_test::kVP9TestVectors;
         it != libvpx_test::kVP9TestVectors + libvpx_test::kNumVP9TestVectors;
         ++it) {
      if (std::find_if(kVP9FailingTestVectors.cbegin(),
                       kVP9FailingTestVectors.cend(), [it](const char *f) {
                         return std::string(f) == *it;
                       }) == kVP9FailingTestVectors.cend()) {
        params.push_back({ flag, *it });
      }
    }
  }
  return params;
}

VP9_INSTANTIATE_TEST_SUITE(PostProcTest,
                           ::testing::ValuesIn(GenerateVP9PassingTestParams()));

// TODO(issue 499602810): remove this test suite after asserts and out of
// bounds accesses are fixed, and the files are migrated to the enabled test
// suite.
INSTANTIATE_TEST_SUITE_P(
    DISABLED_VP9, PostProcTest,
    ::testing::Combine(
        ::testing::Values(
            static_cast<const libvpx_test::CodecFactory *>(&libvpx_test::kVP9)),
        ::testing::ValuesIn(GenerateTestParams(GeneratePostProcFlags(),
                                               kVP9FailingTestVectors))));

constexpr std::array<const char *, 26> kVP9InvalidFiles = {
  "invalid-crbug-1558.ivf", "invalid-crbug-1562.ivf",
  "invalid-crbug-629481.webm", "invalid-crbug-667044.webm",
  "invalid-vp90-01-v3.webm", "invalid-vp90-03-v3.webm",
  "invalid-vp90-2-00-quantizer-00.webm.ivf.s5861_r01-05_b6-.v2.ivf",
  "invalid-vp90-2-00-quantizer-11.webm.ivf.s52984_r01-05_b6-.ivf",
  "invalid-vp90-2-00-quantizer-11.webm.ivf.s52984_r01-05_b6-z.ivf",
  // This file is disabled due to its high memory usage.
  // "invalid-vp90-2-00-quantizer-63.ivf.kf_65527x61446.ivf",
  "invalid-vp90-2-03-size-202x210.webm.ivf.s113306_r01-05_b6-.ivf",
  "invalid-vp90-2-03-size-224x196.webm.ivf.s44156_r01-05_b6-.ivf",
  "invalid-vp90-2-05-resize.ivf.s59293_r01-05_b6-.ivf",
  "invalid-vp90-2-07-frame_parallel-1.webm",
  "invalid-vp90-2-07-frame_parallel-2.webm",
  "invalid-vp90-2-07-frame_parallel-3.webm",
  "invalid-vp90-2-08-tile_1x2_frame_parallel.webm.ivf.s47039_r01-05_b6-.ivf",
  "invalid-vp90-2-08-tile_1x4_frame_parallel_all_key.webm",
  "invalid-vp90-2-08-tile_1x8_frame_parallel.webm.ivf.s288_r01-05_b6-.ivf",
  "invalid-vp90-2-09-aq2.webm.ivf.s3984_r01-05_b6-.v2.ivf",
  "invalid-vp90-2-09-subpixel-00.ivf.s19552_r01-05_b6-.v2.ivf",
  "invalid-vp90-2-09-subpixel-00.ivf.s20492_r01-05_b6-.v2.ivf",
  "invalid-vp90-2-10-show-existing-frame.webm.ivf.s180315_r01-05_b6-.ivf",
  "invalid-vp90-2-12-droppable_1.ivf.s3676_r01-05_b6-.ivf",
  "invalid-vp90-2-12-droppable_1.ivf.s73804_r01-05_b6-.ivf",
  "invalid-vp90-2-21-resize_inter_320x180_5_3-4.webm.ivf.s45551_r01-05_b6-.ivf",
  "invalid-vp91-2-mixedrefcsp-444to420.ivf"
};

VP9_INSTANTIATE_TEST_SUITE(PostProcTestInvalidFiles,
                           ::testing::ValuesIn(GenerateTestParams(
                               GeneratePostProcFlags(), kVP9InvalidFiles)));

#endif  // CONFIG_VP9_DECODER && CONFIG_VP9_POSTPROC

}  // namespace
