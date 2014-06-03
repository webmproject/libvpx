/*
 *  Copyright (c) 2012 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "third_party/googletest/src/include/gtest/gtest.h"
#include "test/acm_random.h"
#include "test/clear_system_state.h"
#include "test/register_state_check.h"
#include "test/util.h"

#include "./vp9_rtcd.h"
#include "vp9/common/vp9_entropy.h"
#include "vpx/vpx_integer.h"

extern "C" {
void vp9_idct8x8_64_add_c(const tran_low_t *input, uint8_t *output, int pitch);
}

using libvpx_test::ACMRandom;

namespace {
typedef void (*fdct_t)(const int16_t *in, tran_low_t *out, int stride);
typedef void (*idct_t)(const tran_low_t *in, uint8_t *out, int stride);
typedef void (*fht_t) (const int16_t *in, tran_low_t *out, int stride,
                       int tx_type);
typedef void (*iht_t) (const tran_low_t *in, uint8_t *out, int stride,
                       int tx_type);

typedef std::tr1::tuple<fdct_t, idct_t, int, int> dct_8x8_param_t;
typedef std::tr1::tuple<fht_t, iht_t, int, int> ht_8x8_param_t;

void fdct8x8_ref(const int16_t *in, tran_low_t *out, int stride, int tx_type) {
  vp9_fdct8x8_c(in, out, stride);
}

void fht8x8_ref(const int16_t *in, tran_low_t *out, int stride, int tx_type) {
  vp9_fht8x8_c(in, out, stride, tx_type);
}

#if CONFIG_VP9_HIGH

void idct8x8_10(const tran_low_t *in, uint8_t *out, int stride) {
  vp9_high_idct8x8_64_add_c(in, out, stride, 10);
}

void idct8x8_12(const tran_low_t *in, uint8_t *out, int stride) {
  vp9_high_idct8x8_64_add_c(in, out, stride, 12);
}

void iht8x8_10(const tran_low_t *in, uint8_t *out, int stride, int tx_type) {
  vp9_high_iht8x8_64_add_c(in, out, stride, tx_type, 10);
}

void iht8x8_12(const tran_low_t *in, uint8_t *out, int stride, int tx_type) {
  vp9_high_iht8x8_64_add_c(in, out, stride, tx_type, 12);
}

#endif

class FwdTrans8x8TestBase {
 public:
  virtual ~FwdTrans8x8TestBase() {}

 protected:
  virtual void RunFwdTxfm(int16_t *in, tran_low_t *out, int stride) = 0;
  virtual void RunInvTxfm(tran_low_t *out, uint8_t *dst, int stride) = 0;

  void RunSignBiasCheck() {
    ACMRandom rnd(ACMRandom::DeterministicSeed());
    DECLARE_ALIGNED_ARRAY(16, int16_t, test_input_block, 64);
    DECLARE_ALIGNED_ARRAY(16, tran_low_t, test_output_block, 64);
    int count_sign_block[64][2];
    const int count_test_block = 100000;

    memset(count_sign_block, 0, sizeof(count_sign_block));

    for (int i = 0; i < count_test_block; ++i) {
      // Initialize a test block with input range [-255, 255].
      for (int j = 0; j < 64; ++j)
        test_input_block[j] = ((rnd.Rand16() >> (16 - bit_depth_)) & mask_) -
                              ((rnd.Rand16() >> (16 - bit_depth_)) & mask_);
      REGISTER_STATE_CHECK(
          RunFwdTxfm(test_input_block, test_output_block, pitch_));

      for (int j = 0; j < 64; ++j) {
        if (test_output_block[j] < 0)
          ++count_sign_block[j][0];
        else if (test_output_block[j] > 0)
          ++count_sign_block[j][1];
      }
    }

    for (int j = 0; j < 64; ++j) {
      const int diff = abs(count_sign_block[j][0] - count_sign_block[j][1]);
      const int max_diff = 1125;
      EXPECT_LT(diff, max_diff)
          << "Error: 8x8 FDCT/FHT has a sign bias > "
          << 1. * max_diff / count_test_block * 100 << "%"
          << " for input range [-255, 255] at index " << j
          << " count0: " << count_sign_block[j][0]
          << " count1: " << count_sign_block[j][1]
          << " diff: " << diff;
    }

    memset(count_sign_block, 0, sizeof(count_sign_block));

    for (int i = 0; i < count_test_block; ++i) {
      // Initialize a test block with input range [-15, 15].
      for (int j = 0; j < 64; ++j)
        test_input_block[j] = ((rnd.Rand16() & mask_) >> 4) -
                              ((rnd.Rand16() & mask_) >> 4);
      REGISTER_STATE_CHECK(
          RunFwdTxfm(test_input_block, test_output_block, pitch_));

      for (int j = 0; j < 64; ++j) {
        if (test_output_block[j] < 0)
          ++count_sign_block[j][0];
        else if (test_output_block[j] > 0)
          ++count_sign_block[j][1];
      }
    }

    for (int j = 0; j < 64; ++j) {
      const int diff = abs(count_sign_block[j][0] - count_sign_block[j][1]);
      const int max_diff = 10000;
      EXPECT_LT(diff, max_diff)
          << "Error: 4x4 FDCT/FHT has a sign bias > "
          << 1. * max_diff / count_test_block * 100 << "%"
          << " for input range [-15, 15] at index " << j
          << " count0: " << count_sign_block[j][0]
          << " count1: " << count_sign_block[j][1]
          << " diff: " << diff;
    }
  }

  void RunRoundTripErrorCheck() {
    ACMRandom rnd(ACMRandom::DeterministicSeed());
    int max_error = 0;
    int total_error = 0;
    const int count_test_block = 100000;
    DECLARE_ALIGNED_ARRAY(16, int16_t, test_input_block, 64);
    DECLARE_ALIGNED_ARRAY(16, tran_low_t, test_temp_block, 64);
    DECLARE_ALIGNED_ARRAY(16, uint8_t, dst, 64);
    DECLARE_ALIGNED_ARRAY(16, uint16_t, dst16, 64);
    DECLARE_ALIGNED_ARRAY(16, uint8_t, src, 64);
    DECLARE_ALIGNED_ARRAY(16, uint16_t, src16, 64);

    for (int i = 0; i < count_test_block; ++i) {
      // Initialize a test block with input range [-mask_, mask_].
      for (int j = 0; j < 64; ++j) {
        if (bit_depth_ == 8) {
          src[j] = rnd.Rand8();
          dst[j] = rnd.Rand8();
          test_input_block[j] = src[j] - dst[j];
        } else {
          src16[j] = rnd.Rand16() & mask_;
          dst16[j] = rnd.Rand16() & mask_;
          test_input_block[j] = src16[j] - dst16[j];
        }
      }

      REGISTER_STATE_CHECK(
          RunFwdTxfm(test_input_block, test_temp_block, pitch_));
      for (int j = 0; j < 64; ++j) {
          if (test_temp_block[j] > 0) {
            test_temp_block[j] += 2;
            test_temp_block[j] /= 4;
            test_temp_block[j] *= 4;
          } else {
            test_temp_block[j] -= 2;
            test_temp_block[j] /= 4;
            test_temp_block[j] *= 4;
          }
      }
      if (bit_depth_ == 8)
        REGISTER_STATE_CHECK(
            RunInvTxfm(test_temp_block, dst, pitch_));
#if CONFIG_VP9_HIGH
      else
        REGISTER_STATE_CHECK(
            RunInvTxfm(test_temp_block, CONVERT_TO_BYTEPTR(dst16), pitch_));
#endif

      for (int j = 0; j < 64; ++j) {
        const int diff = bit_depth_ == 8 ? dst[j] - src[j] :
                                           dst16[j] - src16[j];
        const int error = diff * diff;
        if (max_error < error)
          max_error = error;
        total_error += error;
      }
    }

    EXPECT_GE(1 << (bit_depth_ - 8), max_error)
      << "Error: 8x8 FDCT/IDCT or FHT/IHT has an individual"
      << " roundtrip error > 1";

    EXPECT_GE((count_test_block << (bit_depth_ - 8))/5, total_error)
      << "Error: 8x8 FDCT/IDCT or FHT/IHT has average roundtrip "
      << "error > 1/5 per block";
  }

  void RunExtremalCheck() {
    ACMRandom rnd(ACMRandom::DeterministicSeed());
    int max_error = 0;
    int total_error = 0;
    const int count_test_block = 100000;
    DECLARE_ALIGNED_ARRAY(16, int16_t, test_input_block, 64);
    DECLARE_ALIGNED_ARRAY(16, tran_low_t, test_temp_block, 64);
    DECLARE_ALIGNED_ARRAY(16, uint8_t, dst, 64);
    DECLARE_ALIGNED_ARRAY(16, uint16_t, dst16, 64);
    DECLARE_ALIGNED_ARRAY(16, uint8_t, src, 64);
    DECLARE_ALIGNED_ARRAY(16, uint16_t, src16, 64);

    for (int i = 0; i < count_test_block; ++i) {
      // Initialize a test block with input range [-mask_, mask_].
      for (int j = 0; j < 64; ++j) {
        if (bit_depth_ == 8) {
          src[j] = rnd.Rand8() % 2 ? 255 : 0;
          dst[j] = src[j] > 0 ? 0 : 255;
          test_input_block[j] = src[j] - dst[j];
        } else {
          src16[j] = rnd.Rand8() % 2 ? mask_ : 0;
          dst16[j] = src16[j] > 0 ? 0 : mask_;
          test_input_block[j] = src16[j] - dst16[j];
        }
      }

      REGISTER_STATE_CHECK(
          RunFwdTxfm(test_input_block, test_temp_block, pitch_));
      if (bit_depth_ == 8)
        REGISTER_STATE_CHECK(
            RunInvTxfm(test_temp_block, dst, pitch_));
#if CONFIG_VP9_HIGH
      else
        REGISTER_STATE_CHECK(
            RunInvTxfm(test_temp_block, CONVERT_TO_BYTEPTR(dst16), pitch_));
#endif

      for (int j = 0; j < 64; ++j) {
        const int diff = bit_depth_ == 8 ? dst[j] - src[j] :
                                           dst16[j] - src16[j];
        const int error = diff * diff;
        if (max_error < error)
          max_error = error;
        total_error += error;
      }

      EXPECT_GE(1 << (bit_depth_ - 8), max_error)
          << "Error: Extremal 8x8 FDCT/IDCT or FHT/IHT has"
          << "an individual roundtrip error > 1";

      EXPECT_GE((count_test_block << (bit_depth_ - 8))/5, total_error)
          << "Error: Extremal 8x8 FDCT/IDCT or FHT/IHT has average"
          << " roundtrip error > 1/5 per block";
    }
  }

  int pitch_;
  int tx_type_;
  fht_t fwd_txfm_ref;
  int bit_depth_;
  int mask_;
};

class FwdTrans8x8DCT
    : public FwdTrans8x8TestBase,
      public ::testing::TestWithParam<dct_8x8_param_t> {
 public:
  virtual ~FwdTrans8x8DCT() {}

  virtual void SetUp() {
    fwd_txfm_  = GET_PARAM(0);
    inv_txfm_  = GET_PARAM(1);
    tx_type_   = GET_PARAM(2);
    bit_depth_ = GET_PARAM(3);
    pitch_     = 8;
    fwd_txfm_ref = fdct8x8_ref;
    mask_ = (1 << bit_depth_) - 1;
  }

  virtual void TearDown() { libvpx_test::ClearSystemState(); }

 protected:
  void RunFwdTxfm(int16_t *in, tran_low_t *out, int stride) {
    fwd_txfm_(in, out, stride);
  }
  void RunInvTxfm(tran_low_t *out, uint8_t *dst, int stride) {
    inv_txfm_(out, dst, stride);
  }

  fdct_t fwd_txfm_;
  idct_t inv_txfm_;
};

TEST_P(FwdTrans8x8DCT, SignBiasCheck) {
  RunSignBiasCheck();
}

TEST_P(FwdTrans8x8DCT, RoundTripErrorCheck) {
  RunRoundTripErrorCheck();
}

TEST_P(FwdTrans8x8DCT, ExtremalCheck) {
  RunExtremalCheck();
}

class FwdTrans8x8HT
    : public FwdTrans8x8TestBase,
      public ::testing::TestWithParam<ht_8x8_param_t> {
 public:
  virtual ~FwdTrans8x8HT() {}

  virtual void SetUp() {
    fwd_txfm_  = GET_PARAM(0);
    inv_txfm_  = GET_PARAM(1);
    tx_type_   = GET_PARAM(2);
    bit_depth_ = GET_PARAM(3);
    pitch_     = 8;
    fwd_txfm_ref = fht8x8_ref;
    mask_ = (1 << bit_depth_) - 1;
  }

  virtual void TearDown() { libvpx_test::ClearSystemState(); }

 protected:
  void RunFwdTxfm(int16_t *in, tran_low_t *out, int stride) {
    fwd_txfm_(in, out, stride, tx_type_);
  }
  void RunInvTxfm(tran_low_t *out, uint8_t *dst, int stride) {
    inv_txfm_(out, dst, stride, tx_type_);
  }

  fht_t fwd_txfm_;
  iht_t inv_txfm_;
};

TEST_P(FwdTrans8x8HT, SignBiasCheck) {
  RunSignBiasCheck();
}

TEST_P(FwdTrans8x8HT, RoundTripErrorCheck) {
  RunRoundTripErrorCheck();
}

TEST_P(FwdTrans8x8HT, ExtremalCheck) {
  RunExtremalCheck();
}

using std::tr1::make_tuple;

INSTANTIATE_TEST_CASE_P(
    C, FwdTrans8x8DCT,
    ::testing::Values(
#if CONFIG_VP9_HIGH
        make_tuple(&vp9_fdct8x8_c, &vp9_idct8x8_64_add_c, 0, 8),
        make_tuple(&vp9_high_fdct8x8_c, &idct8x8_10, 0, 10),
        make_tuple(&vp9_high_fdct8x8_c, &idct8x8_12, 0, 12)));
#else
        make_tuple(&vp9_fdct8x8_c, &vp9_idct8x8_64_add_c, 0, 8)));
#endif
INSTANTIATE_TEST_CASE_P(
    C, FwdTrans8x8HT,
    ::testing::Values(
#if CONFIG_VP9_HIGH
        make_tuple(&vp9_fht8x8_c, &vp9_iht8x8_64_add_c, 0, 8),
        make_tuple(&vp9_fht8x8_c, &vp9_iht8x8_64_add_c, 1, 8),
        make_tuple(&vp9_fht8x8_c, &vp9_iht8x8_64_add_c, 2, 8),
        make_tuple(&vp9_fht8x8_c, &vp9_iht8x8_64_add_c, 3, 8),
        make_tuple(&vp9_high_fht8x8_c, &iht8x8_10, 0, 10),
        make_tuple(&vp9_high_fht8x8_c, &iht8x8_10, 1, 10),
        make_tuple(&vp9_high_fht8x8_c, &iht8x8_10, 2, 10),
        make_tuple(&vp9_high_fht8x8_c, &iht8x8_10, 3, 10),
        make_tuple(&vp9_high_fht8x8_c, &iht8x8_12, 0, 12),
        make_tuple(&vp9_high_fht8x8_c, &iht8x8_12, 1, 12),
        make_tuple(&vp9_high_fht8x8_c, &iht8x8_12, 2, 12),
        make_tuple(&vp9_high_fht8x8_c, &iht8x8_12, 3, 12)));
#else
        make_tuple(&vp9_fht8x8_c, &vp9_iht8x8_64_add_c, 0, 8),
        make_tuple(&vp9_fht8x8_c, &vp9_iht8x8_64_add_c, 1, 8),
        make_tuple(&vp9_fht8x8_c, &vp9_iht8x8_64_add_c, 2, 8),
        make_tuple(&vp9_fht8x8_c, &vp9_iht8x8_64_add_c, 3, 8)));
#endif

#if HAVE_NEON_ASM && !CONFIG_HIGH_TRANSFORMS
INSTANTIATE_TEST_CASE_P(
    NEON, FwdTrans8x8DCT,
    ::testing::Values(
        make_tuple(&vp9_fdct8x8_c, &vp9_idct8x8_64_add_neon, 0, 8)));
INSTANTIATE_TEST_CASE_P(
    DISABLED_NEON, FwdTrans8x8HT,
    ::testing::Values(
        make_tuple(&vp9_fht8x8_c, &vp9_iht8x8_64_add_neon, 0, 8),
        make_tuple(&vp9_fht8x8_c, &vp9_iht8x8_64_add_neon, 1, 8),
        make_tuple(&vp9_fht8x8_c, &vp9_iht8x8_64_add_neon, 2, 8),
        make_tuple(&vp9_fht8x8_c, &vp9_iht8x8_64_add_neon, 3, 8)));
#endif

#if HAVE_SSE2 && !CONFIG_HIGH_TRANSFORMS
INSTANTIATE_TEST_CASE_P(
    SSE2, FwdTrans8x8DCT,
    ::testing::Values(
        make_tuple(&vp9_fdct8x8_sse2, &vp9_idct8x8_64_add_sse2, 0, 8)));
INSTANTIATE_TEST_CASE_P(
    SSE2, FwdTrans8x8HT,
    ::testing::Values(
        make_tuple(&vp9_fht8x8_sse2, &vp9_iht8x8_64_add_sse2, 0, 8),
        make_tuple(&vp9_fht8x8_sse2, &vp9_iht8x8_64_add_sse2, 1, 8),
        make_tuple(&vp9_fht8x8_sse2, &vp9_iht8x8_64_add_sse2, 2, 8),
        make_tuple(&vp9_fht8x8_sse2, &vp9_iht8x8_64_add_sse2, 3, 8)));
#endif

#if HAVE_SSSE3 && ARCH_X86_64 && !CONFIG_HIGH_TRANSFORMS
INSTANTIATE_TEST_CASE_P(
    SSSE3, FwdTrans8x8DCT,
    ::testing::Values(
        make_tuple(&vp9_fdct8x8_ssse3, &vp9_idct8x8_64_add_ssse3, 0, 8)));
#endif
}  // namespace
