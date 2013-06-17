/*
 Copyright (c) 2012 The WebM project authors. All Rights Reserved.

 Use of this source code is governed by a BSD-style license
 that can be found in the LICENSE file in the root of the source
 tree. An additional intellectual property rights grant can be found
 in the file PATENTS.  All contributing project authors may
 be found in the AUTHORS file in the root of the source tree.
 */

#include <cstdio>
#include <cstdlib>
#include <string>
#include "third_party/googletest/src/include/gtest/gtest.h"
#include "test/codec_factory.h"
#include "test/decode_test_driver.h"
#include "test/ivf_video_source.h"
#include "test/webm_video_source.h"
#include "test/util.h"
#include "test/md5_helper.h"
extern "C" {
#include "vpx_mem/vpx_mem.h"
}

namespace {
#if CONFIG_VP8_DECODER
const char *kVP8TestVectors[] = {
  "vp80-00-comprehensive-001.ivf",
  "vp80-00-comprehensive-002.ivf", "vp80-00-comprehensive-003.ivf",
  "vp80-00-comprehensive-004.ivf", "vp80-00-comprehensive-005.ivf",
  "vp80-00-comprehensive-006.ivf", "vp80-00-comprehensive-007.ivf",
  "vp80-00-comprehensive-008.ivf", "vp80-00-comprehensive-009.ivf",
  "vp80-00-comprehensive-010.ivf", "vp80-00-comprehensive-011.ivf",
  "vp80-00-comprehensive-012.ivf", "vp80-00-comprehensive-013.ivf",
  "vp80-00-comprehensive-014.ivf", "vp80-00-comprehensive-015.ivf",
  "vp80-00-comprehensive-016.ivf", "vp80-00-comprehensive-017.ivf",
  "vp80-00-comprehensive-018.ivf", "vp80-01-intra-1400.ivf",
  "vp80-01-intra-1411.ivf", "vp80-01-intra-1416.ivf",
  "vp80-01-intra-1417.ivf", "vp80-02-inter-1402.ivf",
  "vp80-02-inter-1412.ivf", "vp80-02-inter-1418.ivf",
  "vp80-02-inter-1424.ivf", "vp80-03-segmentation-01.ivf",
  "vp80-03-segmentation-02.ivf", "vp80-03-segmentation-03.ivf",
  "vp80-03-segmentation-04.ivf", "vp80-03-segmentation-1401.ivf",
  "vp80-03-segmentation-1403.ivf", "vp80-03-segmentation-1407.ivf",
  "vp80-03-segmentation-1408.ivf", "vp80-03-segmentation-1409.ivf",
  "vp80-03-segmentation-1410.ivf", "vp80-03-segmentation-1413.ivf",
  "vp80-03-segmentation-1414.ivf", "vp80-03-segmentation-1415.ivf",
  "vp80-03-segmentation-1425.ivf", "vp80-03-segmentation-1426.ivf",
  "vp80-03-segmentation-1427.ivf", "vp80-03-segmentation-1432.ivf",
  "vp80-03-segmentation-1435.ivf", "vp80-03-segmentation-1436.ivf",
  "vp80-03-segmentation-1437.ivf", "vp80-03-segmentation-1441.ivf",
  "vp80-03-segmentation-1442.ivf", "vp80-04-partitions-1404.ivf",
  "vp80-04-partitions-1405.ivf", "vp80-04-partitions-1406.ivf",
  "vp80-05-sharpness-1428.ivf", "vp80-05-sharpness-1429.ivf",
  "vp80-05-sharpness-1430.ivf", "vp80-05-sharpness-1431.ivf",
  "vp80-05-sharpness-1433.ivf", "vp80-05-sharpness-1434.ivf",
  "vp80-05-sharpness-1438.ivf", "vp80-05-sharpness-1439.ivf",
  "vp80-05-sharpness-1440.ivf", "vp80-05-sharpness-1443.ivf"
};
#endif
#if CONFIG_VP9_DECODER
const char *kVP9TestVectors[] = {
  "vp90-00-akiyo-200.webm", "vp90-00-akiyo-300.webm",
  "vp90-00-akiyo-50.webm", "vp90-00-bowing-150.webm",
  "vp90-00-bowing-25.webm", "vp90-00-bowing-400.webm",
  "vp90-00-bus-100.webm", "vp90-00-bus-2000.webm",
  "vp90-00-bus-300.webm", "vp90-00-bus-4400.webm",
  "vp90-00-bus-800.webm", "vp90-00-cheer-1600.webm",
  "vp90-00-cheer-2800.webm", "vp90-00-cheer-400.webm",
  "vp90-00-cheer-600.webm", "vp90-00-city-1200.webm",
  "vp90-00-city-2000.webm", "vp90-00-city-300.webm",
  "vp90-00-city-600.webm", "vp90-00-coastguard-1200.webm",
  "vp90-00-coastguard-200.webm", "vp90-00-coastguard-3600.webm",
  "vp90-00-coastguard-5200.webm", "vp90-00-container-1000.webm",
  "vp90-00-container-200.webm", "vp90-00-container-50.webm",
  "vp90-00-deadline-1000.webm", "vp90-00-deadline-200.webm",
  "vp90-00-deadline-50.webm", "vp90-00-flower-100.webm",
  "vp90-00-flower-2000.webm", "vp90-00-flower-300.webm",
  "vp90-00-flower-4400.webm", "vp90-00-flower-800.webm",
  "vp90-00-football-1600.webm", "vp90-00-football-2800.webm",
  "vp90-00-football-400.webm", "vp90-00-football-600.webm",
  "vp90-00-foreman-1200.webm", "vp90-00-foreman-2000.webm",
  "vp90-00-foreman-300.webm", "vp90-00-foreman-600.webm",
  "vp90-00-hallmonitor-1200.webm", "vp90-00-hallmonitor-2000.webm",
  "vp90-00-hallmonitor-300.webm", "vp90-00-hallmonitor-600.webm",
  "vp90-00-harbour-1200.webm", "vp90-00-harbour-200.webm",
  "vp90-00-harbour-3600.webm", "vp90-00-harbour-5200.webm",
  "vp90-00-highway-100.webm", "vp90-00-highway-1600.webm",
  "vp90-00-highway-2800.webm", "vp90-00-highway-50.webm",
  "vp90-00-husky-100.webm", "vp90-00-husky-2000.webm",
  "vp90-00-husky-300.webm", "vp90-00-husky-4400.webm",
  "vp90-00-husky-800.webm", "vp90-00-ice-150.webm",
  "vp90-00-ice-400.webm", "vp90-00-ice-800.webm",
  "vp90-00-mobile-1600.webm", "vp90-00-mobile-2800.webm",
  "vp90-00-mobile-400.webm", "vp90-00-mobile-600.webm",
  "vp90-00-motherdaughter-100.webm", "vp90-00-motherdaughter-300.webm",
  "vp90-00-motherdaughter-600.webm", "vp90-00-news-100.webm",
  "vp90-00-news-300.webm", "vp90-00-news-600.webm",
  "vp90-00-pamphlet-150.webm", "vp90-00-pamphlet-25.webm",
  "vp90-00-pamphlet-400.webm", "vp90-00-paris-1000.webm",
  "vp90-00-paris-200.webm", "vp90-00-paris-50.webm",
  "vp90-00-signirene-1000.webm", "vp90-00-signirene-200.webm",
  "vp90-00-signirene-50.webm", "vp90-00-silent-1000.webm",
  "vp90-00-silent-200.webm", "vp90-00-silent-50.webm",
  "vp90-00-soccer-100.webm", "vp90-00-soccer-2000.webm",
  "vp90-00-soccer-300.webm", "vp90-00-soccer-4400.webm",
  "vp90-00-soccer-800.webm", "vp90-00-stefan-1600.webm",
  "vp90-00-stefan-2800.webm", "vp90-00-stefan-400.webm",
  "vp90-00-stefan-600.webm", "vp90-00-students-100.webm",
  "vp90-00-students-300.webm", "vp90-00-students-600.webm",
  "vp90-00-tempete-1200.webm", "vp90-00-tempete-200.webm",
  "vp90-00-tempete-3600.webm", "vp90-00-tempete-5200.webm",
  "vp90-00-tennis-100.webm", "vp90-00-tennis-2000.webm",
  "vp90-00-tennis-300.webm", "vp90-00-tennis-4400.webm",
  "vp90-00-tennis-800.webm", "vp90-00-waterfall-150.webm",
  "vp90-00-waterfall-200.webm", "vp90-00-waterfall-400.webm",
  "vp90-00-waterfall-800.webm",
};
#endif

class TestVectorTest : public ::libvpx_test::DecoderTest,
    public ::libvpx_test::CodecTestWithParam<const char*> {
 protected:
  TestVectorTest() : DecoderTest(GET_PARAM(0)), md5_file_(NULL) {}

  virtual ~TestVectorTest() {
    if (md5_file_)
      fclose(md5_file_);
  }

  void OpenMD5File(const std::string& md5_file_name_) {
    md5_file_ = libvpx_test::OpenTestDataFile(md5_file_name_);
    ASSERT_TRUE(md5_file_) << "Md5 file open failed. Filename: "
        << md5_file_name_;
  }

  virtual void DecompressedFrameHook(const vpx_image_t& img,
                                     const unsigned int frame_number) {
    char expected_md5[33];
    char junk[128];

    // Read correct md5 checksums.
    const int res = fscanf(md5_file_, "%s  %s", expected_md5, junk);
    ASSERT_NE(res, EOF) << "Read md5 data failed";
    expected_md5[32] = '\0';

    ::libvpx_test::MD5 md5_res;
    md5_res.Add(&img);
    const char *actual_md5 = md5_res.Get();

    // Check md5 match.
    ASSERT_STREQ(expected_md5, actual_md5)
        << "Md5 checksums don't match: frame number = " << frame_number;
  }

 private:
  FILE *md5_file_;
};

// This test runs through the whole set of test vectors, and decodes them.
// The md5 checksums are computed for each frame in the video file. If md5
// checksums match the correct md5 data, then the test is passed. Otherwise,
// the test failed.
TEST_P(TestVectorTest, MD5Match) {
  const std::string filename = GET_PARAM(1);
  libvpx_test::CompressedVideoSource *video = NULL;

  // Open compressed video file.
  if (filename.substr(filename.length() - 3, 3) == "ivf") {
    video = new libvpx_test::IVFVideoSource(filename);
  } else if (filename.substr(filename.length() - 4, 4) == "webm") {
    video = new libvpx_test::WebMVideoSource(filename);
  }
  video->Init();

  // Construct md5 file name.
  const std::string md5_filename = filename + ".md5";
  OpenMD5File(md5_filename);

  // Decode frame, and check the md5 matching.
  ASSERT_NO_FATAL_FAILURE(RunLoop(video));
  delete video;
}

VP8_INSTANTIATE_TEST_CASE(TestVectorTest,
                          ::testing::ValuesIn(kVP8TestVectors));
VP9_INSTANTIATE_TEST_CASE(TestVectorTest,
                          ::testing::ValuesIn(kVP9TestVectors));

}  // namespace
