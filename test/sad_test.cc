/*
 *  Copyright (c) 2012 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */


#include <string.h>
#include <limits.h>
#include <stdio.h>

extern "C" {
#include "./vpx_config.h"
#include "./vpx_rtcd.h"
#include "vp8/common/blockd.h"
}

#include "test/acm_random.h"
#include "test/util.h"
#include "third_party/googletest/src/include/gtest/gtest.h"


typedef unsigned int (*sad_m_by_n_fn_t)(const unsigned char *source_ptr,
                                        int source_stride,
                                        const unsigned char *reference_ptr,
                                        int reference_stride,
                                        unsigned int max_sad);

using libvpx_test::ACMRandom;

namespace {
class SADTest : public PARAMS(int, int, sad_m_by_n_fn_t) {
  protected:
  virtual void SetUp() {
    sad_fn_ = GET_PARAM(2);
    height_ = GET_PARAM(1);
    width_ = GET_PARAM(0);
    source_stride_ = width_ * 2;
    reference_stride_ = width_ * 2;
    rnd_.Reset(ACMRandom::DeterministicSeed());
  }

  sad_m_by_n_fn_t sad_fn_;
  virtual unsigned int SAD(unsigned int max_sad) {
    return sad_fn_(source_data_, source_stride_,
                   reference_data_, reference_stride_,
                   max_sad);
  }

  // Sum of Absolute Differences. Given two blocks, calculate the absolute
  // difference between two pixels in the same relative location; accumulate.
  unsigned int ReferenceSAD(unsigned int max_sad) {
    unsigned int sad = 0;

    for (int h = 0; h < height_; ++h) {
      for (int w = 0; w < width_; ++w) {
        sad += abs(source_data_[h * source_stride_ + w]
               - reference_data_[h * reference_stride_ + w]);
      }
      if (sad > max_sad) {
        break;
      }
    }
    return sad;
  }

  void FillConstant(uint8_t *data, int stride, uint8_t fill_constant) {
    for (int h = 0; h < height_; ++h) {
      for (int w = 0; w < width_; ++w) {
        data[h * stride + w] = fill_constant;
      }
    }
  }

  void FillRandom(uint8_t *data, int stride) {
    for (int h = 0; h < height_; ++h) {
      for (int w = 0; w < width_; ++w) {
        data[h * stride + w] = rnd_.Rand8();
      }
    }
  }

  void CheckSad(unsigned int max_sad) {
    unsigned int reference_sad, exp_sad;

    reference_sad = ReferenceSAD(max_sad);
    exp_sad = SAD(max_sad);

    if (reference_sad <= max_sad) {
      ASSERT_EQ(exp_sad, reference_sad);
    } else {
      // Alternative implementations are not required to check max_sad
      ASSERT_GE(exp_sad, reference_sad);
    }
  }

  // Handle blocks up to 16x16 with stride up to 32
  int height_, width_;
  DECLARE_ALIGNED(16, uint8_t, source_data_[16*32]);
  int source_stride_;
  DECLARE_ALIGNED(16, uint8_t, reference_data_[16*32]);
  int reference_stride_;

  ACMRandom rnd_;
};

TEST_P(SADTest, MaxRef) {
  FillConstant(source_data_, source_stride_, 0);
  FillConstant(reference_data_, reference_stride_, 255);
  CheckSad(UINT_MAX);
}

TEST_P(SADTest, MaxSrc) {
  FillConstant(source_data_, source_stride_, 255);
  FillConstant(reference_data_, reference_stride_, 0);
  CheckSad(UINT_MAX);
}

TEST_P(SADTest, ShortRef) {
  int tmp_stride = reference_stride_;
  reference_stride_ >>= 1;
  FillRandom(source_data_, source_stride_);
  FillRandom(reference_data_, reference_stride_);
  CheckSad(UINT_MAX);
  reference_stride_ = tmp_stride;
}

TEST_P(SADTest, UnalignedRef) {
  // The reference frame, but not the source frame, may be unaligned for
  // certain types of searches.
  int tmp_stride = reference_stride_;
  reference_stride_ -= 1;
  FillRandom(source_data_, source_stride_);
  FillRandom(reference_data_, reference_stride_);
  CheckSad(UINT_MAX);
  reference_stride_ = tmp_stride;
}

TEST_P(SADTest, ShortSrc) {
  int tmp_stride = source_stride_;
  source_stride_ >>= 1;
  FillRandom(source_data_, source_stride_);
  FillRandom(reference_data_, reference_stride_);
  CheckSad(UINT_MAX);
  source_stride_ = tmp_stride;
}

TEST_P(SADTest, MaxSAD) {
  // Verify that, when max_sad is set, the implementation does not return a
  // value lower than the reference.
  FillConstant(source_data_, source_stride_, 255);
  FillConstant(reference_data_, reference_stride_, 0);
  CheckSad(128);
}

INSTANTIATE_TEST_CASE_P(C, SADTest, ::testing::Values(
                        std::tr1::make_tuple(16, 16, vp8_sad16x16_c),
                        std::tr1::make_tuple(8, 16, vp8_sad8x16_c),
                        std::tr1::make_tuple(16, 8, vp8_sad16x8_c),
                        std::tr1::make_tuple(8, 8, vp8_sad8x8_c),
                        std::tr1::make_tuple(4, 4, vp8_sad4x4_c)));

// ARM tests
#if HAVE_MEDIA
INSTANTIATE_TEST_CASE_P(MEDIA, SADTest, ::testing::Values(
                        std::tr1::make_tuple(16, 16, vp8_sad16x16_armv6)));

#endif
#if HAVE_NEON
INSTANTIATE_TEST_CASE_P(NEON, SADTest, ::testing::Values(
                        std::tr1::make_tuple(16, 16, vp8_sad16x16_neon),
                        std::tr1::make_tuple(8, 16, vp8_sad8x16_neon),
                        std::tr1::make_tuple(16, 8, vp8_sad16x8_neon),
                        std::tr1::make_tuple(8, 8, vp8_sad8x8_neon),
                        std::tr1::make_tuple(4, 4, vp8_sad4x4_neon)));
#endif

// X86 tests
#if HAVE_MMX
INSTANTIATE_TEST_CASE_P(MMX, SADTest, ::testing::Values(
                        std::tr1::make_tuple(16, 16, vp8_sad16x16_mmx),
                        std::tr1::make_tuple(8, 16, vp8_sad8x16_mmx),
                        std::tr1::make_tuple(16, 8, vp8_sad16x8_mmx),
                        std::tr1::make_tuple(8, 8, vp8_sad8x8_mmx),
                        std::tr1::make_tuple(4, 4, vp8_sad4x4_mmx)));
#endif
#if HAVE_SSE2
INSTANTIATE_TEST_CASE_P(SSE2, SADTest, ::testing::Values(
                        std::tr1::make_tuple(16, 16, vp8_sad16x16_wmt),
                        std::tr1::make_tuple(8, 16, vp8_sad8x16_wmt),
                        std::tr1::make_tuple(16, 8, vp8_sad16x8_wmt),
                        std::tr1::make_tuple(8, 8, vp8_sad8x8_wmt),
                        std::tr1::make_tuple(4, 4, vp8_sad4x4_wmt)));
#endif
#if HAVE_SSSE3
INSTANTIATE_TEST_CASE_P(SSE3, SADTest, ::testing::Values(
                        std::tr1::make_tuple(16, 16, vp8_sad16x16_sse3)));
#endif

}  // namespace
