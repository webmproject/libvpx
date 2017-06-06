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

#include "test/acm_random.h"
#include "vpx/vpx_integer.h"

namespace libvpx_test {

template <typename T>
class Buffer {
 public:
  Buffer(size_t width, size_t height, size_t top_padding, size_t left_padding,
         size_t right_padding, size_t bottom_padding)
      : width_(width), height_(height), top_padding_(top_padding),
        left_padding_(left_padding), right_padding_(right_padding),
        bottom_padding_(bottom_padding), raw_buffer_(NULL) {}

  Buffer(size_t width, size_t height, size_t padding)
      : width_(width), height_(height), top_padding_(padding),
        left_padding_(padding), right_padding_(padding),
        bottom_padding_(padding), raw_buffer_(NULL) {}

  ~Buffer() { delete[] raw_buffer_; }

  T *TopLeftPixel() const;

  size_t stride() const { return stride_; }

  // Set the buffer (excluding padding) to 'value'.
  void Set(const T value);

  // Set the buffer (excluding padding) to the output of ACMRandom function 'b'.
  void Set(ACMRandom *rand_class, T (ACMRandom::*rand_func)());

  // Copy the contents of Buffer 'a' (excluding padding).
  void CopyFrom(const Buffer<T> &a);

  void DumpBuffer() const;

  // Highlight the differences between two buffers if they are the same size.
  void PrintDifference(const Buffer<T> &a) const;

  bool HasPadding() const;

  // Sets all the values in the buffer to 'padding_value'.
  void SetPadding(const T padding_value);

  // Checks if all the values (excluding padding) are equal to 'value' if the
  // Buffers are the same size.
  bool CheckValues(const T value) const;

  // Check that padding matches the expected value or there is no padding.
  bool CheckPadding() const;

  // Compare the non-padding portion of two buffers if they are the same size.
  bool CheckValues(const Buffer<T> &a) const;

  bool Init() {
    EXPECT_GT(width_, 0u);
    EXPECT_GT(height_, 0u);
    stride_ = left_padding_ + width_ + right_padding_;
    raw_size_ = stride_ * (top_padding_ + height_ + bottom_padding_);
    raw_buffer_ = new (std::nothrow) T[raw_size_];
    EXPECT_TRUE(raw_buffer_ != NULL);
    SetPadding(std::numeric_limits<T>::max());
    return !::testing::Test::HasFailure();
  }

 private:
  bool BufferSizesMatch(const Buffer<T> &a) const;

  const size_t width_;
  const size_t height_;
  const size_t top_padding_;
  const size_t left_padding_;
  const size_t right_padding_;
  const size_t bottom_padding_;
  T padding_value_;
  size_t stride_;
  size_t raw_size_;
  T *raw_buffer_;
};

template <typename T>
T *Buffer<T>::TopLeftPixel() const {
  if (!raw_buffer_) return NULL;
  return raw_buffer_ + (top_padding_ * stride()) + left_padding_;
}

template <typename T>
void Buffer<T>::Set(const T value) {
  if (!raw_buffer_) return;
  T *src = TopLeftPixel();
  for (size_t height = 0; height < height_; ++height) {
    for (size_t width = 0; width < width_; ++width) {
      src[width] = value;
    }
    src += stride();
  }
}

template <typename T>
void Buffer<T>::Set(ACMRandom *rand_class, T (ACMRandom::*rand_func)()) {
  if (!raw_buffer_) return;
  T *src = TopLeftPixel();
  for (size_t height = 0; height < height_; ++height) {
    for (size_t width = 0; width < width_; ++width) {
      src[width] = (*rand_class.*rand_func)();
    }
    src += stride();
  }
}

template <typename T>
void Buffer<T>::CopyFrom(const Buffer<T> &a) {
  if (!raw_buffer_) return;
  if (!BufferSizesMatch(a)) return;

  T *a_src = a.TopLeftPixel();
  T *b_src = this->TopLeftPixel();
  for (size_t height = 0; height < height_; ++height) {
    for (size_t width = 0; width < width_; ++width) {
      b_src[width] = a_src[width];
    }
    a_src += a.stride();
    b_src += this->stride();
  }
}

template <typename T>
void Buffer<T>::DumpBuffer() const {
  if (!raw_buffer_) return;
  for (size_t height = 0; height < height_ + top_padding_ + bottom_padding_;
       ++height) {
    for (size_t width = 0; width < stride(); ++width) {
      printf("%4d", raw_buffer_[height + width * stride()]);
    }
    printf("\n");
  }
}

template <typename T>
bool Buffer<T>::HasPadding() const {
  if (!raw_buffer_) return false;
  return top_padding_ || left_padding_ || right_padding_ || bottom_padding_;
}

template <typename T>
void Buffer<T>::PrintDifference(const Buffer<T> &a) const {
  if (!raw_buffer_) return;
  if (!BufferSizesMatch(a)) return;

  T *a_src = a.TopLeftPixel();
  T *b_src = TopLeftPixel();

  printf("This buffer:\n");
  for (size_t height = 0; height < height_; ++height) {
    for (size_t width = 0; width < width_; ++width) {
      if (a_src[width] != b_src[width]) {
        printf("*%3d", b_src[width]);
      } else {
        printf("%4d", b_src[width]);
      }
    }
    printf("\n");
    a_src += a.stride();
    b_src += this->stride();
  }

  a_src = a.TopLeftPixel();
  b_src = TopLeftPixel();

  printf("Reference buffer:\n");
  for (size_t height = 0; height < height_; ++height) {
    for (size_t width = 0; width < width_; ++width) {
      if (a_src[width] != b_src[width]) {
        printf("*%3d", a_src[width]);
      } else {
        printf("%4d", a_src[width]);
      }
    }
    printf("\n");
    a_src += a.stride();
    b_src += this->stride();
  }
}

template <typename T>
void Buffer<T>::SetPadding(const T padding_value) {
  if (!raw_buffer_) return;
  padding_value_ = padding_value;

  T *src = raw_buffer_;
  for (size_t i = 0; i < raw_size_; ++i) {
    src[i] = padding_value;
  }
}

template <typename T>
bool Buffer<T>::CheckValues(const T value) const {
  if (!raw_buffer_) return false;
  T *src = TopLeftPixel();
  for (size_t height = 0; height < height_; ++height) {
    for (size_t width = 0; width < width_; ++width) {
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
  if (!raw_buffer_) return false;
  if (!HasPadding()) return true;

  // Top padding.
  T const *top = raw_buffer_;
  for (size_t i = 0; i < stride() * top_padding_; ++i) {
    if (padding_value_ != top[i]) {
      return false;
    }
  }

  // Left padding.
  T const *left = TopLeftPixel() - left_padding_;
  for (size_t height = 0; height < height_; ++height) {
    for (size_t width = 0; width < left_padding_; ++width) {
      if (padding_value_ != left[width]) {
        return false;
      }
    }
    left += stride();
  }

  // Right padding.
  T const *right = TopLeftPixel() + width_;
  for (size_t height = 0; height < height_; ++height) {
    for (size_t width = 0; width < right_padding_; ++width) {
      if (padding_value_ != right[width]) {
        return false;
      }
    }
    right += stride();
  }

  // Bottom padding
  T const *bottom = raw_buffer_ + (top_padding_ + height_) * stride();
  for (size_t i = 0; i < stride() * bottom_padding_; ++i) {
    if (padding_value_ != bottom[i]) {
      return false;
    }
  }

  return true;
}

template <typename T>
bool Buffer<T>::CheckValues(const Buffer<T> &a) const {
  if (!raw_buffer_) return false;
  if (!BufferSizesMatch(a)) return false;

  T *a_src = a.TopLeftPixel();
  T *b_src = this->TopLeftPixel();
  for (size_t height = 0; height < height_; ++height) {
    for (size_t width = 0; width < width_; ++width) {
      if (a_src[width] != b_src[width]) {
        return false;
      }
    }
    a_src += a.stride();
    b_src += this->stride();
  }
  return true;
}

template <typename T>
bool Buffer<T>::BufferSizesMatch(const Buffer<T> &a) const {
  if (!raw_buffer_) return false;
  if (a.width_ != this->width_ || a.height_ != this->height_) {
    printf(
        "Reference buffer of size %zux%zu does not match this buffer which is "
        "size %zux%zu\n",
        a.width_, a.height_, this->width_, this->height_);
    return false;
  }

  return true;
}
}  // namespace libvpx_test
#endif  // TEST_BUFFER_H_
