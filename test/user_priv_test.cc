/*
 *  Copyright (c) 2013 The WebM project authors. All Rights Reserved.
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
#include "third_party/googletest/src/include/gtest/gtest.h"
#include "./vpx_config.h"
#include "test/codec_factory.h"
#include "test/decode_test_driver.h"
#include "test/ivf_video_source.h"
#include "test/md5_helper.h"
#include "test/util.h"
#if CONFIG_WEBM_IO
#include "test/webm_video_source.h"
#endif
#include "vpx_mem/vpx_mem.h"

namespace {

using std::string;

#if CONFIG_WEBM_IO
// Decodes |filename|. Passes in user_priv data when calling DecodeFrame and
// compares the user_priv from return img with the original user_priv to see if
// they match. Both the pointer values and the values inside the addresses
// should match.
string DecodeFile(const string &filename) {
  libvpx_test::WebMVideoSource video(filename);
  video.Init();

  vpx_codec_dec_cfg_t cfg = {0};
  libvpx_test::VP9Decoder decoder(cfg, 0);

  libvpx_test::MD5 md5;
  int frame_num = 0;
  for (video.Begin(); video.cxdata(); video.Next()) {
    void *user_priv = reinterpret_cast<void *>(&frame_num);
    const vpx_codec_err_t res =
        decoder.DecodeFrame(video.cxdata(), video.frame_size(),
                            (frame_num == 0) ? NULL : user_priv);
    if (res != VPX_CODEC_OK) {
      EXPECT_EQ(VPX_CODEC_OK, res) << decoder.DecodeError();
      break;
    }
    libvpx_test::DxDataIterator dec_iter = decoder.GetDxData();
    const vpx_image_t *img = NULL;

    // Get decompressed data.
    while ((img = dec_iter.Next())) {
      if (frame_num == 0) {
        // user_priv pointer value should be the same.
        EXPECT_EQ(img->user_priv, reinterpret_cast<void *>(NULL)) <<
            "user_priv pointer value does not match.";
      } else {
        // user_priv pointer value should be the same.
        EXPECT_EQ(img->user_priv, reinterpret_cast<void *>(&frame_num)) <<
            "user_priv pointer value does not match.";
        // value in user_priv pointer should also be the same.
        EXPECT_EQ(*reinterpret_cast<int *>(img->user_priv), frame_num) <<
            "Value in user_priv does not match.";
      }
      md5.Add(img);
    }

    frame_num++;
  }
  return string(md5.Get());
}

TEST(UserPrivTest, VideoDecode) {
  // no tiles or frame parallel; this exercises the decoding to test the
  // user_priv.
  EXPECT_STREQ("b35a1b707b28e82be025d960aba039bc",
               DecodeFile("vp90-2-03-size-226x226.webm").c_str());
}

#endif  // CONFIG_WEBM_IO

}  // namespace
