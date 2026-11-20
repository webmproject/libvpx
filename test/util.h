/*
 *  Copyright (c) 2012 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#ifndef VPX_TEST_UTIL_H_
#define VPX_TEST_UTIL_H_

#include <stdio.h>
#include <math.h>
#include <tuple>

#include "gtest/gtest.h"
#include "vpx/vpx_image.h"
#include "vpx_ports/compiler_attributes.h"

#if defined(__linux__) && !VPX_WITH_ASAN &&                                 \
    !__has_feature(memory_sanitizer) && !__has_feature(thread_sanitizer) && \
    !__has_feature(hwaddress_sanitizer) && !defined(__SANITIZE_THREAD__)
#define VPX_TEST_CAN_EXHAUST_ADDRESS_SPACE 1
#else
#define VPX_TEST_CAN_EXHAUST_ADDRESS_SPACE 0
#endif

#if VPX_TEST_CAN_EXHAUST_ADDRESS_SPACE
#include <cstdint>
#include <cstdlib>
#include <vector>
#include <sys/resource.h>
#ifdef __GLIBC__
#include <malloc.h>
#endif

namespace libvpx_test {

// RAII helper that sets RLIMIT_AS to 0 and drains residual allocator free
// lists so subsequent heap allocations return nullptr. Restore() (or
// destruction) restores the original limit and frees drained blocks. Must call
// Restore() before any EXPECT_*/ASSERT_* macros.
class ScopedAddressSpaceExhaustion {
 public:
  ScopedAddressSpaceExhaustion() {
    filler_.reserve(65536);
#ifdef __GLIBC__
    malloc_trim(0);
#endif
    if (getrlimit(RLIMIT_AS, &old_lim_) != 0) return;
    const struct rlimit zero_lim = { 0, old_lim_.rlim_max };
    if (setrlimit(RLIMIT_AS, &zero_lim) != 0) return;
    armed_ = true;

    for (size_t sz : { size_t{ 1 } << 24, size_t{ 1 } << 20, size_t{ 1 } << 16,
                       size_t{ 1 } << 12, size_t{ 256 } }) {
      while (filler_.size() < filler_.capacity()) {
        void *p = malloc(sz);
        if (p == nullptr) break;
        *static_cast<volatile uint8_t *>(p) = 0;
        filler_.push_back(p);
      }
    }
  }

  ScopedAddressSpaceExhaustion(const ScopedAddressSpaceExhaustion &) = delete;
  ScopedAddressSpaceExhaustion &operator=(
      const ScopedAddressSpaceExhaustion &) = delete;

  ~ScopedAddressSpaceExhaustion() { Restore(); }

  bool IsArmed() const { return armed_; }

  bool Restore() {
    if (!armed_) return true;
    const int restore_res = setrlimit(RLIMIT_AS, &old_lim_);
    for (void *p : filler_) free(p);
    filler_.clear();
    armed_ = false;
    return restore_res == 0;
  }

 private:
  struct rlimit old_lim_ = {};
  std::vector<void *> filler_;
  bool armed_ = false;
};

}  // namespace libvpx_test
#endif  // VPX_TEST_CAN_EXHAUST_ADDRESS_SPACE

// Macros
#define GET_PARAM(k) std::get<k>(GetParam())

inline double compute_psnr(const vpx_image_t *img1, const vpx_image_t *img2) {
  assert((img1->fmt == img2->fmt) && (img1->d_w == img2->d_w) &&
         (img1->d_h == img2->d_h));

  const unsigned int width_y = img1->d_w;
  const unsigned int height_y = img1->d_h;
  unsigned int i, j;

  int64_t sqrerr = 0;
  for (i = 0; i < height_y; ++i) {
    for (j = 0; j < width_y; ++j) {
      int64_t d = img1->planes[VPX_PLANE_Y][i * img1->stride[VPX_PLANE_Y] + j] -
                  img2->planes[VPX_PLANE_Y][i * img2->stride[VPX_PLANE_Y] + j];
      sqrerr += d * d;
    }
  }
  double mse = static_cast<double>(sqrerr) / (width_y * height_y);
  double psnr = 100.0;
  if (mse > 0.0) {
    psnr = 10 * log10(255.0 * 255.0 / mse);
  }
  return psnr;
}

#endif  // VPX_TEST_UTIL_H_
