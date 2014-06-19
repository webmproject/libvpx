/*
 *  Copyright (c) 2014 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>
#include "third_party/googletest/src/include/gtest/gtest.h"
#include "./vpx_config.h"
#include "test/codec_factory.h"
#include "test/decode_test_driver.h"
#include "test/ivf_video_source.h"
#include "test/util.h"
#if CONFIG_WEBM_IO
#include "test/webm_video_source.h"
#endif
#include "vpx_mem/vpx_mem.h"

namespace {

class InvalidFileTest
    : public ::libvpx_test::DecoderTest,
      public ::libvpx_test::CodecTestWithParam<const char*> {
 protected:
  InvalidFileTest() : DecoderTest(GET_PARAM(0)), res_file_(NULL) {}

  virtual ~InvalidFileTest() {
    if (res_file_ != NULL)
      fclose(res_file_);
  }

  void OpenResFile(const std::string &res_file_name_) {
    res_file_ = libvpx_test::OpenTestDataFile(res_file_name_);
    ASSERT_TRUE(res_file_ != NULL) << "Result file open failed. Filename: "
        << res_file_name_;
  }

  virtual bool HandleDecodeResult(
      const vpx_codec_err_t res_dec,
      const libvpx_test::CompressedVideoSource &video,
      libvpx_test::Decoder *decoder) {
    EXPECT_TRUE(res_file_ != NULL);
    int expected_res_dec;

    // Read integer result.
    const int res = fscanf(res_file_, "%d", &expected_res_dec);
    EXPECT_NE(res, EOF) << "Read result data failed";

    // Check results match.
    EXPECT_EQ(expected_res_dec, res_dec)
        << "Results don't match: frame number = " << video.frame_number();

    return !HasFailure();
  }

 private:
  FILE *res_file_;
};

TEST_P(InvalidFileTest, ReturnCode) {
  const std::string filename = GET_PARAM(1);
  libvpx_test::CompressedVideoSource *video = NULL;

  // Open compressed video file.
  if (filename.substr(filename.length() - 3, 3) == "ivf") {
    video = new libvpx_test::IVFVideoSource(filename);
  } else if (filename.substr(filename.length() - 4, 4) == "webm") {
#if CONFIG_WEBM_IO
    video = new libvpx_test::WebMVideoSource(filename);
#else
    fprintf(stderr, "WebM IO is disabled, skipping test vector %s\n",
            filename.c_str());
    return;
#endif
  }
  video->Init();

  // Construct result file name. The file holds a list of expected integer
  // results, one for each decoded frame.  Any result that doesn't match
  // the files list will cause a test failure.
  const std::string res_filename = filename + ".res";
  OpenResFile(res_filename);

  // Decode frame, and check the md5 matching.
  ASSERT_NO_FATAL_FAILURE(RunLoop(video));
  delete video;
}

const char *const kVP9InvalidFileTests[] = {
  "invalid-vp90-01.webm"
};

#define NELEMENTS(x) static_cast<int>(sizeof(x) / sizeof(x[0]))

VP9_INSTANTIATE_TEST_CASE(InvalidFileTest,
                          ::testing::ValuesIn(kVP9InvalidFileTests,
                                              kVP9InvalidFileTests +
                                              NELEMENTS(kVP9InvalidFileTests)));

}  // namespace
