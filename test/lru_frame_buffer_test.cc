/*
 *  Copyright (c) 2013 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include <queue>
#include <string>

#include "test/codec_factory.h"
#include "test/decode_test_driver.h"
#include "test/ivf_video_source.h"
#include "test/md5_helper.h"
#include "test/util.h"
#include "test/webm_video_source.h"

namespace {

const int kVideoNameParam = 1;

const char *kLRUTestVectors[] = {
  "vp90-2-02-size-lf-1920x1080.webm",
  "vp90-2-05-resize.ivf",
};

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

// Class for testing libvpx is using the least recently
// used frame buffer when a new buffer is requested.
class LRUFrameBufferTest
    : public ::libvpx_test::DecoderTest,
      public ::libvpx_test::CodecTestWithParam<const char*> {
 protected:
  struct FrameBufferMD5Sum {
    int frame_buffer_index;
    vpx_image_t img;
    std::string md5;
  };

  LRUFrameBufferTest()
      : DecoderTest(GET_PARAM(::libvpx_test::kCodecFactoryParam)),
        num_buffers_(0),
        num_jitter_buffers_(0),
        frame_buffers_(NULL) {}

  virtual ~LRUFrameBufferTest() {
    for (int i = 0; i < num_buffers_; ++i) {
      delete [] frame_buffers_[i].data;
    }
    delete [] frame_buffers_;
  }

  virtual void PreDecodeFrameHook(
      const libvpx_test::CompressedVideoSource &video,
      libvpx_test::Decoder *decoder) {
    // Use external buffers for testing jitter buffers.
    if (num_jitter_buffers_ > 0 && video.frame_number() == 0) {
      const int max_reference_buffers = 8;

      // Add 1 for a work buffer.
      num_buffers_ = max_reference_buffers + 1 + num_jitter_buffers_;

      // Have libvpx use frame buffers we create.
      frame_buffers_ = new vpx_codec_frame_buffer_t[num_buffers_];
      memset(frame_buffers_, 0, sizeof(frame_buffers_[0]) * num_buffers_);

      decoder->SetExternalFrameBuffers(frame_buffers_, num_buffers_,
                                       realloc_vp9_frame_buffer, NULL);
    }

    // Turn on frame buffer LRU cache.
    decoder->Control(VP9D_SET_FRAME_BUFFER_LRU_CACHE, 1);
  }

  virtual void DecompressedFrameHook(const vpx_image_t &img,
                                     const unsigned int frame_number) {
    const uint32_t ximg_y_plane = 0;
    const uint8_t *const y_buffer = img.planes[ximg_y_plane];

    // Find which external buffer contains the y_buffer.
    int i = 0;
    for (i = 0; i < num_buffers_; ++i) {
      if (y_buffer >= frame_buffers_[i].data &&
          y_buffer < (frame_buffers_[i].data + frame_buffers_[i].size)) {
        break;
      }
    }

    FrameBufferMD5Sum fb_md5;
    fb_md5.frame_buffer_index = i;
    fb_md5.img = img;

    libvpx_test::MD5 md5;
    md5.Add(&img);
    fb_md5.md5 = md5.Get();
    jitter_buffer_md5_sums_.push(fb_md5);

    // Check to see if any of the reconstructed image changed.
    if (jitter_buffer_md5_sums_.size() >
        static_cast<size_t>(num_jitter_buffers_)) {
      fb_md5 = jitter_buffer_md5_sums_.front();

      libvpx_test::MD5 md5;
      md5.Add(&fb_md5.img);
      const std::string check_str = md5.Get();

      ASSERT_EQ(fb_md5.md5, check_str);
      jitter_buffer_md5_sums_.pop();
    }
  }

  libvpx_test::CompressedVideoSource *OpenCompressedFile(
      const std::string &filename) {
    if (filename.substr(filename.length() - 3, 3) == "ivf") {
      return new libvpx_test::IVFVideoSource(filename);
    } else if (filename.substr(filename.length() - 4, 4) == "webm") {
      return new libvpx_test::WebMVideoSource(filename);
    }
    return NULL;
  }

  void set_num_jitter_buffers(int num_buffers) {
    num_jitter_buffers_ = num_buffers;
  }

 private:
  // Total number of external frame buffers.
  int num_buffers_;
  int num_jitter_buffers_;

  // External frame buffers used by libvpx.
  vpx_codec_frame_buffer_t *frame_buffers_;

  // Save the md5 checksums for later comparison.
  std::queue<FrameBufferMD5Sum> jitter_buffer_md5_sums_;
};

// This test runs through a set of test vectors, and decodes them.
// Libvpx will call into the application to allocate a frame buffer when
// needed. The md5 checksums are computed for each frame after it is
// decoded and stored to be checked later. After a jitter frame buffer
// has expired, the md5 checksum is computed again for the expired jitter
// buffer frame and checked against the md5 checksum after the frame was
// decoded. If md5 checksums match, then the test is passed. Otherwise,
// the test failed.
TEST_P(LRUFrameBufferTest, CheckLRUOneJitterBuffer) {
  const std::string filename = GET_PARAM(kVideoNameParam);

  set_num_jitter_buffers(1);

  libvpx_test::CompressedVideoSource *const video =
      OpenCompressedFile(filename);
  video->Init();

  // Decode frame, and check the md5 matching.
  ASSERT_NO_FATAL_FAILURE(RunLoop(video));
  delete video;
}

TEST_P(LRUFrameBufferTest, CheckLRUFourJitterBuffers) {
  const std::string filename = GET_PARAM(kVideoNameParam);

  set_num_jitter_buffers(4);

  libvpx_test::CompressedVideoSource *const video =
      OpenCompressedFile(filename);
  video->Init();

  // Decode frame, and check the md5 matching.
  ASSERT_NO_FATAL_FAILURE(RunLoop(video));
  delete video;
}

TEST_P(LRUFrameBufferTest, CheckLRUEightJitterBuffers) {
  const std::string filename = GET_PARAM(kVideoNameParam);

  set_num_jitter_buffers(8);

  libvpx_test::CompressedVideoSource *const video =
      OpenCompressedFile(filename);
  video->Init();

  // Decode frame, and check the md5 matching.
  ASSERT_NO_FATAL_FAILURE(RunLoop(video));
  delete video;
}

VP9_INSTANTIATE_TEST_CASE(LRUFrameBufferTest,
                          ::testing::ValuesIn(kLRUTestVectors));
}  // namespace
