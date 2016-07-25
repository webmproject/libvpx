/*
 *  Copyright (c) 2012 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include <limits.h>
#include <stdio.h>
#include <string.h>

#include "third_party/googletest/src/include/gtest/gtest.h"

#include "./vpx_config.h"

#include "test/acm_random.h"
#include "test/clear_system_state.h"
#include "test/register_state_check.h"
#include "test/util.h"
#include "vpx_dsp/ssim.h"
#include "vpx_mem/vpx_mem.h"

extern "C"
double vpx_get_ssim_metrics(uint8_t *img1, int img1_pitch,
                            uint8_t *img2, int img2_pitch,
                            int width, int height,
                            Ssimv *sv2, Metrics *m,
                            int do_inconsistency);

using libvpx_test::ACMRandom;

namespace {
class ConsistencyTestBase : public ::testing::Test {
 public:
  ConsistencyTestBase(int width, int height) : width_(width), height_(height) {}

  static void SetUpTestCase() {
    source_data_[0] = reinterpret_cast<uint8_t*>(
        vpx_memalign(kDataAlignment, kDataBufferSize));
    reference_data_[0] = reinterpret_cast<uint8_t*>(
        vpx_memalign(kDataAlignment, kDataBufferSize));
    source_data_[1] = reinterpret_cast<uint8_t*>(
        vpx_memalign(kDataAlignment, kDataBufferSize));
    reference_data_[1] = reinterpret_cast<uint8_t*>(
        vpx_memalign(kDataAlignment, kDataBufferSize));
    ssim_array_ = new Ssimv[kDataBufferSize / 16];
  }

  static void ClearSsim() {
    memset(ssim_array_, 0, kDataBufferSize / 16);
  }
  static void TearDownTestCase() {
    vpx_free(source_data_[0]);
    source_data_[0] = NULL;
    vpx_free(reference_data_[0]);
    reference_data_[0] = NULL;
    vpx_free(source_data_[1]);
    source_data_[1] = NULL;
    vpx_free(reference_data_[1]);
    reference_data_[1] = NULL;

    delete[] ssim_array_;
  }

  virtual void TearDown() {
    libvpx_test::ClearSystemState();
  }

 protected:
  // Handle frames up to 640x480
  static const int kDataAlignment = 16;
  static const int kDataBufferSize = 640*480;

  virtual void SetUp() {
    source_stride_ = (width_ + 31) & ~31;
    reference_stride_ = width_ * 2;
    rnd_.Reset(ACMRandom::DeterministicSeed());
  }

  void FillRandom(uint8_t *data, int stride, int width, int height) {
    for (int h = 0; h < height; ++h) {
      for (int w = 0; w < width; ++w) {
        data[h * stride + w] = rnd_.Rand8();
      }
    }
  }

  void FillRandom(uint8_t *data, int stride) {
    FillRandom(data, stride, width_, height_);
  }

  void Copy(uint8_t *reference, uint8_t *source) {
    memcpy(reference, source, kDataBufferSize);
  }

  void Blur(uint8_t *data, int stride, int taps) {
    int sum = 0;
    int half_taps = taps / 2;
    for (int h = 0; h < height_; ++h) {
      for (int w = 0; w < taps; ++w) {
        sum += data[w + h * stride];
      }
      for (int w = taps; w < width_; ++w) {
        sum += data[w + h * stride] - data[w - taps + h * stride];
        data[w - half_taps + h * stride] = (sum + half_taps) / taps;
      }
    }
    for (int w = 0; w < width_; ++w) {
      for (int h = 0; h < taps; ++h) {
        sum += data[h + w * stride];
      }
      for (int h = taps; h < height_; ++h) {
        sum += data[w + h * stride] - data[(h - taps) * stride + w];
        data[(h - half_taps) * stride + w] = (sum + half_taps) / taps;
      }
    }
  }
  int width_, height_;
  static uint8_t* source_data_[2];
  int source_stride_;
  static uint8_t* reference_data_[2];
  int reference_stride_;
  static Ssimv *ssim_array_;
  Metrics metrics_;

  ACMRandom rnd_;
};

uint8_t* ConsistencyTestBase::source_data_[2] = {NULL, NULL};
uint8_t* ConsistencyTestBase::reference_data_[2] = {NULL, NULL};
Ssimv* ConsistencyTestBase::ssim_array_ = NULL;



using std::tr1::make_tuple;

//------------------------------------------------------------------------------
// C functions

}  // namespace
