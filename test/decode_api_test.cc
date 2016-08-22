/*
 *  Copyright (c) 2014 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include "third_party/googletest/src/include/gtest/gtest.h"

#include "./vpx_config.h"
#include "test/ivf_video_source.h"
#include "aom/vp8dx.h"
#include "aom/vpx_decoder.h"

namespace {

#define NELEMENTS(x) static_cast<int>(sizeof(x) / sizeof(x[0]))

TEST(DecodeAPI, InvalidParams) {
  static const vpx_codec_iface_t *kCodecs[] = {
#if CONFIG_VP10_DECODER
    &vpx_codec_vp10_dx_algo,
#endif
  };
  uint8_t buf[1] = { 0 };
  vpx_codec_ctx_t dec;

  EXPECT_EQ(VPX_CODEC_INVALID_PARAM, vpx_codec_dec_init(NULL, NULL, NULL, 0));
  EXPECT_EQ(VPX_CODEC_INVALID_PARAM, vpx_codec_dec_init(&dec, NULL, NULL, 0));
  EXPECT_EQ(VPX_CODEC_INVALID_PARAM, vpx_codec_decode(NULL, NULL, 0, NULL, 0));
  EXPECT_EQ(VPX_CODEC_INVALID_PARAM, vpx_codec_decode(NULL, buf, 0, NULL, 0));
  EXPECT_EQ(VPX_CODEC_INVALID_PARAM,
            vpx_codec_decode(NULL, buf, NELEMENTS(buf), NULL, 0));
  EXPECT_EQ(VPX_CODEC_INVALID_PARAM,
            vpx_codec_decode(NULL, NULL, NELEMENTS(buf), NULL, 0));
  EXPECT_EQ(VPX_CODEC_INVALID_PARAM, vpx_codec_destroy(NULL));
  EXPECT_TRUE(vpx_codec_error(NULL) != NULL);

  for (int i = 0; i < NELEMENTS(kCodecs); ++i) {
    EXPECT_EQ(VPX_CODEC_INVALID_PARAM,
              vpx_codec_dec_init(NULL, kCodecs[i], NULL, 0));

    EXPECT_EQ(VPX_CODEC_OK, vpx_codec_dec_init(&dec, kCodecs[i], NULL, 0));
    EXPECT_EQ(VPX_CODEC_UNSUP_BITSTREAM,
              vpx_codec_decode(&dec, buf, NELEMENTS(buf), NULL, 0));
    EXPECT_EQ(VPX_CODEC_INVALID_PARAM,
              vpx_codec_decode(&dec, NULL, NELEMENTS(buf), NULL, 0));
    EXPECT_EQ(VPX_CODEC_INVALID_PARAM, vpx_codec_decode(&dec, buf, 0, NULL, 0));

    EXPECT_EQ(VPX_CODEC_OK, vpx_codec_destroy(&dec));
  }
}

}  // namespace
