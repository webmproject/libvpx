/*
 *  Copyright (c) 2016 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include <assert.h>
#include <string>
#include <vector>
#include "third_party/googletest/src/include/gtest/gtest.h"
#include "test/codec_factory.h"
#include "test/encode_test_driver.h"
#include "test/i420_video_source.h"
#include "test/md5_helper.h"
#include "test/util.h"

namespace {
// The number of frames to be encoded/decoded
const int kLimit = 8;
// Skip 1 frame to check the frame decoding independency.
const int kSkip = 5;
const int kTileSize = 1;
const int kTIleSizeInPixels = (kTileSize << 6);
// Fake width and height so that they can be multiples of the tile size.
const int kImgWidth = 704;
const int kImgHeight = 576;

class VP10ExtTileTest
    : public ::libvpx_test::EncoderTest,
      public ::libvpx_test::CodecTestWith2Params<libvpx_test::TestMode, int> {
 protected:
  VP10ExtTileTest()
      : EncoderTest(GET_PARAM(0)),
        encoding_mode_(GET_PARAM(1)),
        set_cpu_used_(GET_PARAM(2)) {
    init_flags_ = VPX_CODEC_USE_PSNR;
    vpx_codec_dec_cfg_t cfg = vpx_codec_dec_cfg_t();
    cfg.w = kImgWidth;
    cfg.h = kImgHeight;

    decoder_ = codec_->CreateDecoder(cfg, 0);
    decoder_->Control(VP10_SET_DECODE_TILE_ROW, -1);
    decoder_->Control(VP10_SET_DECODE_TILE_COL, -1);

    // Allocate buffer to store tile image.
    vpx_img_alloc(&tile_img_, VPX_IMG_FMT_I420, kImgWidth, kImgHeight, 32);

    md5_.clear();
    tile_md5_.clear();
  }

  virtual ~VP10ExtTileTest() {
    vpx_img_free(&tile_img_);
    delete decoder_;
  }

  virtual void SetUp() {
    InitializeConfig();
    SetMode(encoding_mode_);

    cfg_.g_lag_in_frames = 0;
    cfg_.rc_end_usage = VPX_VBR;
    cfg_.g_error_resilient = 1;

    cfg_.rc_max_quantizer = 56;
    cfg_.rc_min_quantizer = 0;
  }

  virtual void PreEncodeFrameHook(::libvpx_test::VideoSource * video,
                                  ::libvpx_test::Encoder *encoder) {
    if (video->frame() == 0) {
      // Encode setting
      encoder->Control(VP8E_SET_CPUUSED, set_cpu_used_);
      encoder->Control(VP8E_SET_ENABLEAUTOALTREF, 0);
      encoder->Control(VP9E_SET_FRAME_PARALLEL_DECODING, 1);

      // The tile size is 64x64.
      encoder->Control(VP9E_SET_TILE_COLUMNS, kTileSize);
      encoder->Control(VP9E_SET_TILE_ROWS, kTileSize);
#if CONFIG_EXT_PARTITION
      // Always use 64x64 max partition.
      encoder->Control(VP10E_SET_SUPERBLOCK_SIZE, VPX_SUPERBLOCK_SIZE_64X64);
#endif
    }

    if (video->frame() == 1) {
      frame_flags_ = VP8_EFLAG_NO_UPD_LAST | VP8_EFLAG_NO_UPD_GF |
          VP8_EFLAG_NO_UPD_ARF;
    }
  }

  virtual void DecompressedFrameHook(const vpx_image_t &img,
                                     vpx_codec_pts_t pts) {
    // Skip 1 already decoded frame to be consistent with the decoder in this
    // test.
    if (pts == (vpx_codec_pts_t)kSkip)
      return;

    // Calculate MD5 as the reference.
    ::libvpx_test::MD5 md5_res;
    md5_res.Add(&img);
    md5_.push_back(md5_res.Get());
  }

  virtual void FramePktHook(const vpx_codec_cx_pkt_t *pkt) {
    // Skip decoding 1 frame.
    if (pkt->data.frame.pts == (vpx_codec_pts_t)kSkip)
      return;

    bool IsLastFrame = (pkt->data.frame.pts == (vpx_codec_pts_t)(kLimit - 1));

    // Decode the first (kLimit - 1) frames as whole frame, and decode the last
    // frame in single tiles.
    for (int r = 0; r < kImgHeight / kTIleSizeInPixels; ++r) {
      for (int c = 0; c < kImgWidth / kTIleSizeInPixels; ++c) {
        if (!IsLastFrame) {
          decoder_->Control(VP10_SET_DECODE_TILE_ROW, -1);
          decoder_->Control(VP10_SET_DECODE_TILE_COL, -1);
        } else {
          decoder_->Control(VP10_SET_DECODE_TILE_ROW, r);
          decoder_->Control(VP10_SET_DECODE_TILE_COL, c);
        }

        const vpx_codec_err_t res = decoder_->DecodeFrame(
            reinterpret_cast<uint8_t*>(pkt->data.frame.buf),
            pkt->data.frame.sz);
        if (res != VPX_CODEC_OK) {
          abort_ = true;
          ASSERT_EQ(VPX_CODEC_OK, res);
        }
        const vpx_image_t *img = decoder_->GetDxData().Next();

        if (!IsLastFrame) {
          if (img) {
            ::libvpx_test::MD5 md5_res;
            md5_res.Add(img);
            tile_md5_.push_back(md5_res.Get());
          }
          break;
        }

        const int kMaxMBPlane = 3;
        for (int plane = 0; plane < kMaxMBPlane; ++plane) {
          const int shift = (plane == 0) ? 0 : 1;
          int tile_height = kTIleSizeInPixels >> shift;
          int tile_width = kTIleSizeInPixels >> shift;

          for (int tr = 0; tr < tile_height; ++tr) {
            memcpy(tile_img_.planes[plane] +
                   tile_img_.stride[plane] * (r * tile_height + tr) +
                   c * tile_width,
                   img->planes[plane] + img->stride[plane] * tr, tile_width);
          }
        }
      }

      if (!IsLastFrame)
        break;
    }

    if (IsLastFrame) {
      ::libvpx_test::MD5 md5_res;
      md5_res.Add(&tile_img_);
      tile_md5_.push_back(md5_res.Get());
    }
  }

  ::libvpx_test::TestMode encoding_mode_;
  int set_cpu_used_;
  ::libvpx_test::Decoder *decoder_;
  vpx_image_t tile_img_;
  std::vector<std::string> md5_;
  std::vector<std::string> tile_md5_;
};

TEST_P(VP10ExtTileTest, DecoderResultTest) {
  ::libvpx_test::I420VideoSource video("hantro_collage_w352h288.yuv",
                                       kImgWidth, kImgHeight, 30, 1, 0, kLimit);
  cfg_.rc_target_bitrate = 500;
  cfg_.g_error_resilient = VPX_ERROR_RESILIENT_DEFAULT;
  cfg_.g_lag_in_frames = 0;
  cfg_.g_threads = 1;

  // Tile encoding
  init_flags_ = VPX_CODEC_USE_PSNR;
  ASSERT_NO_FATAL_FAILURE(RunLoop(&video));

  // Compare to check if two vectors are equal.
  ASSERT_EQ(md5_, tile_md5_);
}

VP10_INSTANTIATE_TEST_CASE(
    // Now only test 2-pass mode.
    VP10ExtTileTest,
    ::testing::Values(::libvpx_test::kTwoPassGood),
    ::testing::Range(0, 4));
}  // namespace
