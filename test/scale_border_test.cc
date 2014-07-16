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

#include "test/clear_system_state.h"
#include "test/register_state_check.h"

#include "./vpx_config.h"
#include "./vpx_scale_rtcd.h"
#include "vpx_mem/vpx_mem.h"
#include "vpx_scale/yv12config.h"

namespace {

typedef void (*ExtendFrameBorderFunc)(YV12_BUFFER_CONFIG *ybf);

class ExtendBorderTest
    : public ::testing::TestWithParam<ExtendFrameBorderFunc> {
 public:
  virtual ~ExtendBorderTest() {
    libvpx_test::ClearSystemState();
  }

  void ResetImage(int width, int height) {
    width_ = width;
    height_ = height;
    vpx_memset(&img_, 0, sizeof(img_));
    ASSERT_EQ(0, vp8_yv12_alloc_frame_buffer(&img_, width_, height_,
                                             VP8BORDERINPIXELS));

    vpx_memset(img_.buffer_alloc, kBufFiller, img_.frame_size);
    FillPlane(img_.y_buffer, img_.y_crop_width, img_.y_crop_height,
              img_.y_stride);
    FillPlane(img_.u_buffer, img_.uv_crop_width, img_.uv_crop_height,
              img_.uv_stride);
    FillPlane(img_.v_buffer, img_.uv_crop_width, img_.uv_crop_height,
              img_.uv_stride);

    vpx_memset(&ref_img_, 0, sizeof(ref_img_));
    ASSERT_EQ(0, vp8_yv12_alloc_frame_buffer(&ref_img_, width_, height_,
                                             VP8BORDERINPIXELS));

    vpx_memset(ref_img_.buffer_alloc, kBufFiller, ref_img_.frame_size);
    FillPlane(ref_img_.y_buffer, ref_img_.y_crop_width, ref_img_.y_crop_height,
              ref_img_.y_stride);
    FillPlane(ref_img_.u_buffer,
              ref_img_.uv_crop_width, ref_img_.uv_crop_height,
              ref_img_.uv_stride);
    FillPlane(ref_img_.v_buffer,
              ref_img_.uv_crop_width, ref_img_.uv_crop_height,
              ref_img_.uv_stride);
  }

  void DeallocImage() {
    vp8_yv12_de_alloc_frame_buffer(&img_);
    vp8_yv12_de_alloc_frame_buffer(&ref_img_);
  }

 private:
  static const int kBufFiller = 123;
  static const int kBufMax = kBufFiller - 1;

  virtual void SetUp() {
    extend_fn_ = GetParam();
  }

  static void FillPlane(uint8_t *buf, int width, int height, int stride) {
    for (int y = 0; y < height; ++y) {
      for (int x = 0; x < width; ++x) {
        buf[x + (y * stride)] = (x + (width * y)) % kBufMax;
      }
    }
  }

  void ReferenceExtendBorder() {
    ExtendPlane(ref_img_.y_buffer,
                ref_img_.y_crop_width, ref_img_.y_crop_height,
                ref_img_.y_width, ref_img_.y_height,
                ref_img_.y_stride,
                ref_img_.border);
    ExtendPlane(ref_img_.u_buffer,
                ref_img_.uv_crop_width, ref_img_.uv_crop_height,
                ref_img_.uv_width, ref_img_.uv_height,
                ref_img_.uv_stride,
                ref_img_.border / 2);
    ExtendPlane(ref_img_.v_buffer,
                ref_img_.uv_crop_width, ref_img_.uv_crop_height,
                ref_img_.uv_width, ref_img_.uv_height,
                ref_img_.uv_stride,
                ref_img_.border / 2);
  }

  static void ExtendPlane(uint8_t *buf, int crop_width, int crop_height,
                          int width, int height, int stride, int padding) {
    // Copy the outermost visible pixel to a distance of at least 'padding.'
    // The buffers are allocated such that there may be excess space outside the
    // padding. As long as the minimum amount of padding is achieved it is not
    // necessary to fill this space as well.
    uint8_t *left = buf - padding;
    uint8_t *right = buf + crop_width;
    const int right_extend = padding + (width - crop_width);
    const int bottom_extend = padding + (height - crop_height);

    // Fill the border pixels from the nearest image pixel.
    for (int y = 0; y < crop_height; ++y) {
      vpx_memset(left, left[padding], padding);
      vpx_memset(right, right[-1], right_extend);
      left += stride;
      right += stride;
    }

    left = buf - padding;
    uint8_t *top = left - (stride * padding);
    // The buffer does not always extend as far as the stride.
    // Equivalent to padding + width + padding.
    const int extend_width = padding + crop_width + right_extend;

    // The first row was already extended to the left and right. Copy it up.
    for (int y = 0; y < padding; ++y) {
      vpx_memcpy(top, left, extend_width);
      top += stride;
    }

    uint8_t *bottom = left + (crop_height * stride);
    for (int y = 0; y <  bottom_extend; ++y) {
      vpx_memcpy(bottom, left + (crop_height - 1) * stride, extend_width);
      bottom += stride;
    }
  }

  void ExtendBorder() {
    ASM_REGISTER_STATE_CHECK(extend_fn_(&img_));
  }

  void CompareImages() {
    EXPECT_EQ(ref_img_.frame_size, img_.frame_size);
    EXPECT_EQ(0, memcmp(ref_img_.buffer_alloc, img_.buffer_alloc,
                        ref_img_.frame_size));
  }

 protected:
  void RunTest() {
#if ARCH_ARM
    // Some arm devices OOM when trying to allocate the largest buffers.
    static const int kNumSizesToTest = 6;
#else
    static const int kNumSizesToTest = 7;
#endif
    static const int kSizesToTest[] = {1, 15, 33, 145, 512, 1025, 16383};
    for (int h = 0; h < kNumSizesToTest; ++h) {
      for (int w = 0; w < kNumSizesToTest; ++w) {
        ResetImage(kSizesToTest[w], kSizesToTest[h]);
        ExtendBorder();
        ReferenceExtendBorder();
        CompareImages();
        DeallocImage();
      }
    }
  }

  YV12_BUFFER_CONFIG img_;
  YV12_BUFFER_CONFIG ref_img_;
  ExtendFrameBorderFunc extend_fn_;
  int width_;
  int height_;
};

TEST_P(ExtendBorderTest, ExtendBorder) {
  ASSERT_NO_FATAL_FAILURE(RunTest());
}

INSTANTIATE_TEST_CASE_P(C, ExtendBorderTest,
                        ::testing::Values(vp8_yv12_extend_frame_borders_c));

#if HAVE_NEON
INSTANTIATE_TEST_CASE_P(DISABLED_NEON, ExtendBorderTest,
                        ::testing::Values(vp8_yv12_extend_frame_borders_neon));
#endif

}  // namespace
