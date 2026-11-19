/*
 *  Copyright (c) 2026 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

// Tests for the YV12 frame buffer allocation helpers in vpx_scale. The frame
// sizes here are computed independently in 64-bit so they cannot wrap.

#include <stdint.h>

#include <cstring>

#include "gtest/gtest.h"
#include "vpx_scale/yv12config.h"

namespace {

struct FrameLayout {
  uint64_t frame_size;
  uint64_t y_offset;
  uint64_t u_offset;
  uint64_t v_offset;
};

// The buffer layout vp8_yv12_realloc_frame_buffer() must produce, evaluated
// in 64-bit as ground truth.
FrameLayout ComputeLayout(int width, int height, int border) {
  const uint64_t aligned_width = ((uint64_t)width + 15) & ~(uint64_t)15;
  const uint64_t aligned_height = ((uint64_t)height + 15) & ~(uint64_t)15;
  const uint64_t y_stride =
      ((aligned_width + 2 * (uint64_t)border) + 31) & ~(uint64_t)31;
  const uint64_t yplane_size = (aligned_height + 2 * (uint64_t)border) * y_stride;
  const uint64_t uv_stride = y_stride >> 1;
  const uint64_t uvplane_size =
      (aligned_height / 2 + (uint64_t)border) * uv_stride;
  FrameLayout layout;
  layout.frame_size = yplane_size + 2 * uvplane_size;
  layout.y_offset = (uint64_t)border * y_stride + (uint64_t)border;
  layout.u_offset =
      yplane_size + (uint64_t)(border / 2) * uv_stride + (uint64_t)(border / 2);
  layout.v_offset = yplane_size + uvplane_size +
                    (uint64_t)(border / 2) * uv_stride +
                    (uint64_t)(border / 2);
  return layout;
}

void ExpectExactLayout(const YV12_BUFFER_CONFIG &ybf, const FrameLayout &layout) {
  EXPECT_EQ(layout.frame_size, (uint64_t)ybf.frame_size);
  EXPECT_EQ(ybf.buffer_alloc + layout.y_offset, ybf.y_buffer);
  EXPECT_EQ(ybf.buffer_alloc + layout.u_offset, ybf.u_buffer);
  EXPECT_EQ(ybf.buffer_alloc + layout.v_offset, ybf.v_buffer);
}

// Normal frames must still allocate and report the exact frame size.
TEST(YV12ConfigTest, Vp8NormalFrames) {
  static const int kDims[][2] = { { 16, 16 }, { 352, 288 }, { 16383, 16383 } };
  for (size_t i = 0; i < sizeof(kDims) / sizeof(kDims[0]); ++i) {
    const int width = kDims[i][0];
    const int height = kDims[i][1];
    YV12_BUFFER_CONFIG ybf;
    memset(&ybf, 0, sizeof(ybf));
    ASSERT_EQ(0, vp8_yv12_alloc_frame_buffer(&ybf, width, height,
                                             VP8BORDERINPIXELS));
    ExpectExactLayout(ybf, ComputeLayout(width, height, VP8BORDERINPIXELS));
    vp8_yv12_de_alloc_frame_buffer(&ybf);
  }
}

// Dimensions whose true plane sizes overflow 32-bit arithmetic must never be
// accepted with a wrapped (undersized) allocation: either the call fails, or
// the buffer is the full 64-bit size.
TEST(YV12ConfigTest, Vp8WrappingSizesNotUndersized) {
  static const int kWrap[][3] = {
    // width, height, border.
    { 65536, 65536, 32 },    // yplane_size wraps to 8392704, true 4303360000.
    { 46341, 46341, 32 },    // yplane_size wraps to -2139779584.
    { 100000, 100000, 32 },  // yplane_size wraps to 1422869504.
  };
  for (size_t i = 0; i < sizeof(kWrap) / sizeof(kWrap[0]); ++i) {
    const int width = kWrap[i][0];
    const int height = kWrap[i][1];
    const int border = kWrap[i][2];
    YV12_BUFFER_CONFIG ybf;
    memset(&ybf, 0, sizeof(ybf));
    const int ret = vp8_yv12_alloc_frame_buffer(&ybf, width, height, border);
    if (ret != 0) {
      EXPECT_EQ(NULL, ybf.buffer_alloc);
      continue;
    }
    ExpectExactLayout(ybf, ComputeLayout(width, height, border));
    vp8_yv12_de_alloc_frame_buffer(&ybf);
  }
}

// A huge border passes the "multiple of 32" rule but wraps the plane size
// arithmetic to a tiny allocation; it must be rejected outright.
TEST(YV12ConfigTest, Vp8OversizedBorderRejected) {
  YV12_BUFFER_CONFIG ybf;
  memset(&ybf, 0, sizeof(ybf));
  EXPECT_NE(0, vp8_yv12_alloc_frame_buffer(&ybf, 32, 32, 1 << 29));
  EXPECT_EQ(NULL, ybf.buffer_alloc);
}

// Negative dimensions and borders must be rejected.
TEST(YV12ConfigTest, Vp8NegativeArgsRejected) {
  YV12_BUFFER_CONFIG ybf;
  memset(&ybf, 0, sizeof(ybf));
  EXPECT_NE(0, vp8_yv12_alloc_frame_buffer(&ybf, -16, 16, 32));
  EXPECT_NE(0, vp8_yv12_alloc_frame_buffer(&ybf, 16, -16, 32));
  EXPECT_NE(0, vp8_yv12_alloc_frame_buffer(&ybf, 16, 16, -32));
  EXPECT_NE(0, vp8_yv12_alloc_frame_buffer(&ybf, 16, 16, 33));
  EXPECT_EQ(NULL, ybf.buffer_alloc);
}

}  // namespace
