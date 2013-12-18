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

#include "test/codec_factory.h"
#include "test/decode_test_driver.h"
#include "test/ivf_video_source.h"
#include "test/md5_helper.h"
#include "test/test_vectors.h"
#include "test/util.h"
#include "test/webm_video_source.h"

namespace {

const int kVideoNameParam = 1;
const char kVP9TestFile[] = "vp90-2-02-size-lf-1920x1080.webm";

// Callback used by libvpx to request the application to allocate a frame
// buffer of at least |new_size| in bytes.
int realloc_vp9_frame_buffer(void *user_priv, size_t new_size,
                             vpx_codec_frame_buffer_t *fb) {
  (void)user_priv;
  if (fb == NULL)
    return -1;

  delete [] fb->data;
  fb->data = new uint8_t[new_size];
  fb->size = new_size;
  return 0;
}

// Callback will not allocate data for frame buffer.
int zero_realloc_vp9_frame_buffer(void *user_priv, size_t new_size,
                                  vpx_codec_frame_buffer_t *fb) {
  (void)user_priv;
  if (fb == NULL)
    return -1;

  delete [] fb->data;
  fb->data = NULL;
  fb->size = new_size;
  return 0;
}

// Callback will allocate one less byte.
int one_less_byte_realloc_vp9_frame_buffer(void *user_priv, size_t new_size,
                                           vpx_codec_frame_buffer_t *fb) {
  (void)user_priv;
  if (fb == NULL)
    return -1;

  delete [] fb->data;

  const size_t error_size = new_size - 1;
  fb->data = new uint8_t[error_size];
  fb->size = error_size;
  return 0;
}

// Class for testing passing in external frame buffers to libvpx.
class ExternalFrameBufferMD5Test
    : public ::libvpx_test::DecoderTest,
      public ::libvpx_test::CodecTestWithParam<const char*> {
 protected:
  ExternalFrameBufferMD5Test()
      : DecoderTest(GET_PARAM(::libvpx_test::kCodecFactoryParam)),
        md5_file_(NULL),
        num_buffers_(0),
        frame_buffers_(NULL) {}

  virtual ~ExternalFrameBufferMD5Test() {
    for (int i = 0; i < num_buffers_; ++i) {
      delete [] frame_buffers_[i].data;
    }
    delete [] frame_buffers_;

    if (md5_file_ != NULL)
      fclose(md5_file_);
  }

  virtual void PreDecodeFrameHook(
      const libvpx_test::CompressedVideoSource &video,
      libvpx_test::Decoder *decoder) {
    if (num_buffers_ > 0 && video.frame_number() == 0) {
      // Have libvpx use frame buffers we create.
      frame_buffers_ = new vpx_codec_frame_buffer_t[num_buffers_];
      memset(frame_buffers_, 0, sizeof(frame_buffers_[0]) * num_buffers_);

      ASSERT_EQ(VPX_CODEC_OK,
                decoder->SetExternalFrameBuffers(
                    frame_buffers_, num_buffers_,
                    realloc_vp9_frame_buffer, NULL));
    }
  }

  void OpenMD5File(const std::string &md5_file_name_) {
    md5_file_ = libvpx_test::OpenTestDataFile(md5_file_name_);
    ASSERT_TRUE(md5_file_ != NULL) << "Md5 file open failed. Filename: "
        << md5_file_name_;
  }

  virtual void DecompressedFrameHook(const vpx_image_t &img,
                                     const unsigned int frame_number) {
    ASSERT_TRUE(md5_file_ != NULL);
    char expected_md5[33];
    char junk[128];

    // Read correct md5 checksums.
    const int res = fscanf(md5_file_, "%s  %s", expected_md5, junk);
    ASSERT_NE(EOF, res) << "Read md5 data failed";
    expected_md5[32] = '\0';

    ::libvpx_test::MD5 md5_res;
    md5_res.Add(&img);
    const char *const actual_md5 = md5_res.Get();

    // Check md5 match.
    ASSERT_STREQ(expected_md5, actual_md5)
        << "Md5 checksums don't match: frame number = " << frame_number;
  }

  void set_num_buffers(int num_buffers) { num_buffers_ = num_buffers; }
  int num_buffers() const { return num_buffers_; }

 private:
  FILE *md5_file_;
  int num_buffers_;
  vpx_codec_frame_buffer_t *frame_buffers_;
};

class ExternalFrameBufferTest : public ::testing::Test {
 protected:
  ExternalFrameBufferTest()
      : video_(NULL),
        decoder_(NULL),
        num_buffers_(0),
        frame_buffers_(NULL) {}

  virtual void SetUp() {
    video_ = new libvpx_test::WebMVideoSource(kVP9TestFile);
    video_->Init();
    video_->Begin();

    vpx_codec_dec_cfg_t cfg = {0};
    decoder_ = new libvpx_test::VP9Decoder(cfg, 0);
  }

  virtual void TearDown() {
    for (int i = 0; i < num_buffers_; ++i) {
      delete [] frame_buffers_[i].data;
    }
    delete [] frame_buffers_;
    delete decoder_;
    delete video_;
  }

  // Passes the external frame buffer information to libvpx.
  vpx_codec_err_t SetExternalFrameBuffers(
      int num_buffers,
      vpx_realloc_frame_buffer_cb_fn_t cb) {
    if (num_buffers > 0) {
      num_buffers_ = num_buffers;

      // Have libvpx use frame buffers we create.
      frame_buffers_ = new vpx_codec_frame_buffer_t[num_buffers_];
      memset(frame_buffers_, 0, sizeof(frame_buffers_[0]) * num_buffers_);
    }

    return decoder_->SetExternalFrameBuffers(frame_buffers_, num_buffers_,
                                             cb, NULL);
  }

  // Pass Null frame buffer list to libvpx.
  vpx_codec_err_t SetNullFrameBuffers(
      int num_buffers,
      vpx_realloc_frame_buffer_cb_fn_t cb) {
    return decoder_->SetExternalFrameBuffers(NULL, num_buffers,
                                             cb, NULL);
  }

  vpx_codec_err_t DecodeOneFrame() {
    const vpx_codec_err_t res =
        decoder_->DecodeFrame(video_->cxdata(), video_->frame_size());
    if (res == VPX_CODEC_OK)
      video_->Next();
    return res;
  }

  vpx_codec_err_t DecodeRemainingFrames() {
    for (; video_->cxdata(); video_->Next()) {
      const vpx_codec_err_t res =
          decoder_->DecodeFrame(video_->cxdata(), video_->frame_size());
      if (res != VPX_CODEC_OK)
        return res;

      libvpx_test::DxDataIterator dec_iter = decoder_->GetDxData();
      const vpx_image_t *img = NULL;

      // Get decompressed data
      while ((img = dec_iter.Next())) {
      }
    }
    return VPX_CODEC_OK;
  }

  libvpx_test::WebMVideoSource *video_;
  libvpx_test::VP9Decoder *decoder_;
  int num_buffers_;
  vpx_codec_frame_buffer_t *frame_buffers_;
};


// This test runs through the set of test vectors, and decodes them.
// Libvpx will call into the application to allocate a frame buffer when
// needed. The md5 checksums are computed for each frame in the video file.
// If md5 checksums match the correct md5 data, then the test is passed.
// Otherwise, the test failed.
TEST_P(ExternalFrameBufferMD5Test, ExtFBMD5Match) {
  const std::string filename = GET_PARAM(kVideoNameParam);
  libvpx_test::CompressedVideoSource *video = NULL;

  // Number of buffers equals number of possible reference buffers(8), plus
  // one working buffer, plus four jitter buffers.
  const int num_buffers = 13;
  set_num_buffers(num_buffers);

#if CONFIG_VP8_DECODER
  // Tell compiler we are not using kVP8TestVectors.
  (void)libvpx_test::kVP8TestVectors;
#endif

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

TEST_F(ExternalFrameBufferTest, NineFrameBuffers) {
  // Minimum number of external frame buffers for VP9 is
  // #VP9_MAXIMUM_REF_BUFFERS + #VPX_MAXIMUM_WORK_BUFFERS.
  const int num_buffers = VP9_MAXIMUM_REF_BUFFERS + VPX_MAXIMUM_WORK_BUFFERS;
  ASSERT_EQ(VPX_CODEC_OK,
            SetExternalFrameBuffers(num_buffers, realloc_vp9_frame_buffer));
  ASSERT_EQ(VPX_CODEC_OK, DecodeRemainingFrames());
}

TEST_F(ExternalFrameBufferTest, EightJitterBuffers) {
  // Number of buffers equals #VP9_MAXIMUM_REF_BUFFERS +
  // #VPX_MAXIMUM_WORK_BUFFERS + eight jitter buffers.
  const int jitter_buffers = 8;
  const int num_buffers =
      VP9_MAXIMUM_REF_BUFFERS + VPX_MAXIMUM_WORK_BUFFERS + jitter_buffers;
  ASSERT_EQ(VPX_CODEC_OK,
            SetExternalFrameBuffers(num_buffers, realloc_vp9_frame_buffer));
  ASSERT_EQ(VPX_CODEC_OK, DecodeRemainingFrames());
}

TEST_F(ExternalFrameBufferTest, NotEnoughBuffers) {
  // Minimum number of external frame buffers for VP9 is
  // #VP9_MAXIMUM_REF_BUFFERS + #VPX_MAXIMUM_WORK_BUFFERS. Set one less.
  const int num_buffers =
      VP9_MAXIMUM_REF_BUFFERS + VPX_MAXIMUM_WORK_BUFFERS - 1;
  ASSERT_EQ(VPX_CODEC_INVALID_PARAM,
            SetExternalFrameBuffers(num_buffers, realloc_vp9_frame_buffer));
}

TEST_F(ExternalFrameBufferTest, NullFrameBufferList) {
  // Number of buffers equals #VP9_MAXIMUM_REF_BUFFERS +
  // #VPX_MAXIMUM_WORK_BUFFERS + four jitter buffers.
  const int jitter_buffers = 4;
  const int num_buffers =
      VP9_MAXIMUM_REF_BUFFERS + VPX_MAXIMUM_WORK_BUFFERS + jitter_buffers;
  ASSERT_EQ(VPX_CODEC_INVALID_PARAM,
            SetNullFrameBuffers(num_buffers, realloc_vp9_frame_buffer));
}

TEST_F(ExternalFrameBufferTest, NullRealloc) {
  // Number of buffers equals #VP9_MAXIMUM_REF_BUFFERS +
  // #VPX_MAXIMUM_WORK_BUFFERS + four jitter buffers.
  const int jitter_buffers = 4;
  const int num_buffers =
      VP9_MAXIMUM_REF_BUFFERS + VPX_MAXIMUM_WORK_BUFFERS + jitter_buffers;
  ASSERT_EQ(VPX_CODEC_OK,
            SetExternalFrameBuffers(num_buffers,
                                    zero_realloc_vp9_frame_buffer));
  ASSERT_EQ(VPX_CODEC_MEM_ERROR, DecodeOneFrame());
}

TEST_F(ExternalFrameBufferTest, ReallocOneLessByte) {
  // Number of buffers equals #VP9_MAXIMUM_REF_BUFFERS +
  // #VPX_MAXIMUM_WORK_BUFFERS + four jitter buffers.
  const int jitter_buffers = 4;
  const int num_buffers =
      VP9_MAXIMUM_REF_BUFFERS + VPX_MAXIMUM_WORK_BUFFERS + jitter_buffers;
  ASSERT_EQ(VPX_CODEC_OK,
            SetExternalFrameBuffers(num_buffers,
                                    one_less_byte_realloc_vp9_frame_buffer));
  ASSERT_EQ(VPX_CODEC_MEM_ERROR, DecodeOneFrame());
}

VP9_INSTANTIATE_TEST_CASE(ExternalFrameBufferMD5Test,
                          ::testing::ValuesIn(libvpx_test::kVP9TestVectors));
}  // namespace
