/*
 *  Copyright (c) 2016 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#ifndef TEST_BUFFER_H_
#define TEST_BUFFER_H_

#include <stdio.h>

#include <limits>

#include "third_party/googletest/src/include/gtest/gtest.h"

#include "vpx/vpx_integer.h"

namespace libvpx_test {

template <typename T>
class Buffer {
 public:
  Buffer(int width, int height, int top_padding, int left_padding,
         int right_padding, int bottom_padding)
      : width_(width), height_(height), top_padding_(top_padding),
        left_padding_(left_padding), right_padding_(right_padding),
        bottom_padding_(bottom_padding) {
    Init();
  }

  Buffer(int width, int height, int padding)
      : width_(width), height_(height), top_padding_(padding),
        left_padding_(padding), right_padding_(padding),
        bottom_padding_(padding) {
    Init();
  }

  ~Buffer() { delete[] raw_buffer_; }

  T *TopLeftPixel() const;

  int stride() const { return stride_; }

  // Set the buffer (excluding padding) to 'value'.
  void Set(const int value);

  void DumpBuffer() const;

  bool HasPadding() const;

  // Sets all the values in the buffer to 'padding_value'.
  void SetPadding(const int padding_value);

  // Checks if all the values (excluding padding) are equal to 'value'.
  bool CheckValues(const int value) const;

  // Check that padding matches the expected value or there is no padding.
  bool CheckPadding() const;

 private:
  void Init() {
    ASSERT_GT(width_, 0);
    ASSERT_GT(height_, 0);
    ASSERT_GE(top_padding_, 0);
    ASSERT_GE(left_padding_, 0);
    ASSERT_GE(right_padding_, 0);
    ASSERT_GE(bottom_padding_, 0);
    stride_ = left_padding_ + width_ + right_padding_;
    raw_size_ = stride_ * (top_padding_ + height_ + bottom_padding_);
    raw_buffer_ = new (std::nothrow) T[raw_size_];
    ASSERT_TRUE(raw_buffer_ != NULL);
    SetPadding(std::numeric_limits<T>::max());
  }

  const int width_;
  const int height_;
  const int top_padding_;
  const int left_padding_;
  const int right_padding_;
  const int bottom_padding_;
  int padding_value_;
  int stride_;
  int raw_size_;
  T *raw_buffer_;
};

template <typename T>
T *Buffer<T>::TopLeftPixel() const {
  return raw_buffer_ + (top_padding_ * stride()) + left_padding_;
}

template <typename T>
void Buffer<T>::Set(const int value) {
  T *src = TopLeftPixel();
  for (int height = 0; height < height_; ++height) {
    for (int width = 0; width < width_; ++width) {
      src[width] = value;
    }
    src += stride();
  }
}

template <typename T>
void Buffer<T>::DumpBuffer() const {
  for (int height = 0; height < height_ + top_padding_ + bottom_padding_;
       ++height) {
    for (int width = 0; width < stride(); ++width) {
      printf("%4d", raw_buffer_[height + width * stride()]);
    }
    printf("\n");
  }
}

template <typename T>
bool Buffer<T>::HasPadding() const {
  return top_padding_ || left_padding_ || right_padding_ || bottom_padding_;
}

template <typename T>
void Buffer<T>::SetPadding(const int padding_value) {
  padding_value_ = padding_value;

  T *src = raw_buffer_;
  for (int i = 0; i < raw_size_; ++i) {
    src[i] = padding_value;
  }
}

template <typename T>
bool Buffer<T>::CheckValues(const int value) const {
  T *src = TopLeftPixel();
  for (int height = 0; height < height_; ++height) {
    for (int width = 0; width < width_; ++width) {
      if (value != src[width]) {
        return false;
      }
    }
    src += stride();
  }
  return true;
}

template <typename T>
bool Buffer<T>::CheckPadding() const {
  if (!HasPadding()) {
    return true;
  }

  // Top padding.
  T const *top = raw_buffer_;
  for (int i = 0; i < stride() * top_padding_; ++i) {
    if (padding_value_ != top[i]) {
      return false;
    }
  }

  // Left padding.
  T const *left = TopLeftPixel() - left_padding_;
  for (int height = 0; height < height_; ++height) {
    for (int width = 0; width < left_padding_; ++width) {
      if (padding_value_ != left[width]) {
        return false;
      }
    }
    left += stride();
  }

  // Right padding.
  T const *right = TopLeftPixel() + width_;
  for (int height = 0; height < height_; ++height) {
    for (int width = 0; width < right_padding_; ++width) {
      if (padding_value_ != right[width]) {
        return false;
      }
    }
    right += stride();
  }

  // Bottom padding
  T const *bottom = raw_buffer_ + (top_padding_ + height_) * stride();
  for (int i = 0; i < stride() * bottom_padding_; ++i) {
    if (padding_value_ != bottom[i]) {
      return false;
    }
  }

  return true;
}
}  // namespace libvpx_test
#endif  // TEST_BUFFER_H_
