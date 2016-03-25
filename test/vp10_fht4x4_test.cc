/*
 *  Copyright (c) 2016 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include "third_party/googletest/src/include/gtest/gtest.h"

#include "./vp10_rtcd.h"
#include "./vpx_dsp_rtcd.h"

#include "test/acm_random.h"
#include "test/clear_system_state.h"
#include "test/register_state_check.h"
#include "test/transform_test_base.h"
#include "test/util.h"
#include "vpx_ports/mem.h"

using libvpx_test::ACMRandom;

namespace {
typedef void (*IhtFunc)(const tran_low_t *in, uint8_t *out, int stride,
                        int tx_type);

using libvpx_test::FhtFunc;
typedef std::tr1::tuple<FhtFunc, IhtFunc, int, vpx_bit_depth_t, int> Ht4x4Param;

void fht4x4_ref(const int16_t *in, tran_low_t *out, int stride,
                int tx_type) {
  vp10_fht4x4_c(in, out, stride, tx_type);
}

#if CONFIG_VP9_HIGHBITDEPTH
typedef void (*IhighbdHtFunc)(const tran_low_t *in, uint8_t *out, int stride,
                              int tx_type, int bd);

typedef std::tr1::tuple<FhtFunc, IhighbdHtFunc, int, vpx_bit_depth_t, int>
HighbdHt4x4Param;

void highbe_fht4x4_ref(const int16_t *in, tran_low_t *out, int stride,
                       int tx_type) {
  vp10_highbd_fht4x4_c(in, out, stride, tx_type);
}
#endif  // CONFIG_VP9_HIGHBITDEPTH

class VP10Trans4x4HT
    : public libvpx_test::TransformTestBase,
      public ::testing::TestWithParam<Ht4x4Param> {
 public:
  virtual ~VP10Trans4x4HT() {}

  virtual void SetUp() {
    fwd_txfm_ = GET_PARAM(0);
    inv_txfm_ = GET_PARAM(1);
    tx_type_  = GET_PARAM(2);
    pitch_    = 4;
    fwd_txfm_ref = fht4x4_ref;
    bit_depth_ = GET_PARAM(3);
    mask_ = (1 << bit_depth_) - 1;
    num_coeffs_ = GET_PARAM(4);
  }
  virtual void TearDown() { libvpx_test::ClearSystemState(); }

 protected:
  void RunFwdTxfm(const int16_t *in, tran_low_t *out, int stride) {
    fwd_txfm_(in, out, stride, tx_type_);
  }

  void RunInvTxfm(const tran_low_t *out, uint8_t *dst, int stride) {
    inv_txfm_(out, dst, stride, tx_type_);
  }

  FhtFunc fwd_txfm_;
  IhtFunc inv_txfm_;
};

TEST_P(VP10Trans4x4HT, CoeffCheck) {
  RunCoeffCheck();
}

#if CONFIG_VP9_HIGHBITDEPTH
class VP10HighbdTrans4x4HT
    : public libvpx_test::TransformTestBase,
      public ::testing::TestWithParam<HighbdHt4x4Param> {
 public:
  virtual ~VP10HighbdTrans4x4HT() {}

  virtual void SetUp() {
    fwd_txfm_ = GET_PARAM(0);
    inv_txfm_ = GET_PARAM(1);
    tx_type_  = GET_PARAM(2);
    pitch_    = 4;
    fwd_txfm_ref = highbe_fht4x4_ref;
    bit_depth_ = GET_PARAM(3);
    mask_ = (1 << bit_depth_) - 1;
    num_coeffs_ = GET_PARAM(4);
  }
  virtual void TearDown() { libvpx_test::ClearSystemState(); }

 protected:
  void RunFwdTxfm(const int16_t *in, tran_low_t *out, int stride) {
    fwd_txfm_(in, out, stride, tx_type_);
  }

  void RunInvTxfm(const tran_low_t *out, uint8_t *dst, int stride) {
    inv_txfm_(out, dst, stride, tx_type_, bit_depth_);
  }

  FhtFunc fwd_txfm_;
  IhighbdHtFunc inv_txfm_;
};

TEST_P(VP10HighbdTrans4x4HT, HighbdCoeffCheck) {
  RunCoeffCheck();
}
#endif  // CONFIG_VP9_HIGHBITDEPTH

#if CONFIG_EXT_TX
TEST(VP10Trans4x4HTSpeedTest, C_version) {
    ACMRandom rnd(ACMRandom::DeterministicSeed());
    const int count_test_block = 20000;
    int bit_depth = 8;
    int mask = (1 << bit_depth) - 1;
    const int num_coeffs = 16;
    int16_t *input = new int16_t[num_coeffs];
    tran_low_t *output = new tran_low_t[num_coeffs];
    const int stride = 4;
    int tx_type;

    for (int i = 0; i < count_test_block; ++i) {
      for (int j = 0; j < num_coeffs; ++j) {
        input[j] = (rnd.Rand8() & mask) - (rnd.Rand8() & mask);
      }
      for (tx_type = V_DCT; tx_type <= H_FLIPADST; ++tx_type) {
        vp10_fht4x4_c(input, output, stride, tx_type);
      }
    }

    delete[] input;
    delete[] output;
}
#endif  // CONFIG_EXT_TX

#if HAVE_SSE2 && CONFIG_EXT_TX
TEST(VP10Trans4x4HTSpeedTest, SSE2_version) {
    ACMRandom rnd(ACMRandom::DeterministicSeed());
    const int count_test_block = 20000;
    int bit_depth = 8;
    int mask = (1 << bit_depth) - 1;
    const int num_coeffs = 16;
    int16_t *input = new int16_t[num_coeffs];
    tran_low_t *output = new tran_low_t[num_coeffs];
    const int stride = 4;
    int tx_type;

    for (int i = 0; i < count_test_block; ++i) {
      for (int j = 0; j < num_coeffs; ++j) {
        input[j] = (rnd.Rand8() & mask) - (rnd.Rand8() & mask);
      }
      for (tx_type = V_DCT; tx_type <= H_FLIPADST; ++tx_type) {
        vp10_fht4x4_sse2(input, output, stride, tx_type);
      }
    }

    delete[] input;
    delete[] output;
}
#endif  // HAVE_SSE2 && CONFIG_EXT_TX

using std::tr1::make_tuple;

#if HAVE_SSE2
INSTANTIATE_TEST_CASE_P(
    SSE2, VP10Trans4x4HT,
    ::testing::Values(
#if !CONFIG_EXT_TX
      make_tuple(&vp10_fht4x4_sse2, &vp10_iht4x4_16_add_sse2, 0,
                 VPX_BITS_8, 16),
      make_tuple(&vp10_fht4x4_sse2, &vp10_iht4x4_16_add_sse2, 1,
                 VPX_BITS_8, 16),
      make_tuple(&vp10_fht4x4_sse2, &vp10_iht4x4_16_add_sse2, 2,
                 VPX_BITS_8, 16),
      make_tuple(&vp10_fht4x4_sse2, &vp10_iht4x4_16_add_sse2, 3,
                 VPX_BITS_8, 16)));
#else
      make_tuple(&vp10_fht4x4_sse2, &vp10_iht4x4_16_add_sse2, 0,
                 VPX_BITS_8, 16),
      make_tuple(&vp10_fht4x4_sse2, &vp10_iht4x4_16_add_sse2, 1,
                 VPX_BITS_8, 16),
      make_tuple(&vp10_fht4x4_sse2, &vp10_iht4x4_16_add_sse2, 2,
                 VPX_BITS_8, 16),
      make_tuple(&vp10_fht4x4_sse2, &vp10_iht4x4_16_add_sse2, 3,
                 VPX_BITS_8, 16),
      make_tuple(&vp10_fht4x4_sse2, &vp10_iht4x4_16_add_sse2, 4,
                 VPX_BITS_8, 16),
      make_tuple(&vp10_fht4x4_sse2, &vp10_iht4x4_16_add_sse2, 5,
                 VPX_BITS_8, 16),
      make_tuple(&vp10_fht4x4_sse2, &vp10_iht4x4_16_add_sse2, 6,
                 VPX_BITS_8, 16),
      make_tuple(&vp10_fht4x4_sse2, &vp10_iht4x4_16_add_sse2, 7,
                 VPX_BITS_8, 16),
      make_tuple(&vp10_fht4x4_sse2, &vp10_iht4x4_16_add_sse2, 8,
                 VPX_BITS_8, 16),
      make_tuple(&vp10_fht4x4_sse2, &vp10_iht4x4_16_add_sse2, 10,
                 VPX_BITS_8, 16),
      make_tuple(&vp10_fht4x4_sse2, &vp10_iht4x4_16_add_sse2, 11,
                 VPX_BITS_8, 16),
      make_tuple(&vp10_fht4x4_sse2, &vp10_iht4x4_16_add_sse2, 12,
                 VPX_BITS_8, 16),
      make_tuple(&vp10_fht4x4_sse2, &vp10_iht4x4_16_add_sse2, 13,
                 VPX_BITS_8, 16),
      make_tuple(&vp10_fht4x4_sse2, &vp10_iht4x4_16_add_sse2, 14,
                 VPX_BITS_8, 16),
      make_tuple(&vp10_fht4x4_sse2, &vp10_iht4x4_16_add_sse2, 15,
                 VPX_BITS_8, 16)));
#endif  // !CONFIG_EXT_TX
#endif  // HAVE_SSE2

#if HAVE_SSE4_1 && CONFIG_VP9_HIGHBITDEPTH
INSTANTIATE_TEST_CASE_P(
    SSE4_1, VP10HighbdTrans4x4HT,
    ::testing::Values(
#if !CONFIG_EXT_TX
      make_tuple(&vp10_highbd_fht4x4_sse4_1, &vp10_highbd_iht4x4_16_add_c, 0,
                 VPX_BITS_10, 16),
      make_tuple(&vp10_highbd_fht4x4_sse4_1, &vp10_highbd_iht4x4_16_add_c, 1,
                 VPX_BITS_10, 16),
      make_tuple(&vp10_highbd_fht4x4_sse4_1, &vp10_highbd_iht4x4_16_add_c, 2,
                 VPX_BITS_10, 16),
      make_tuple(&vp10_highbd_fht4x4_sse4_1, &vp10_highbd_iht4x4_16_add_c, 3,
                 VPX_BITS_10, 16),
      make_tuple(&vp10_highbd_fht4x4_sse4_1, &vp10_highbd_iht4x4_16_add_c, 0,
                 VPX_BITS_12, 16),
      make_tuple(&vp10_highbd_fht4x4_sse4_1, &vp10_highbd_iht4x4_16_add_c, 1,
                 VPX_BITS_12, 16),
      make_tuple(&vp10_highbd_fht4x4_sse4_1, &vp10_highbd_iht4x4_16_add_c, 2,
                 VPX_BITS_12, 16),
      make_tuple(&vp10_highbd_fht4x4_sse4_1, &vp10_highbd_iht4x4_16_add_c, 3,
                 VPX_BITS_12, 16)));
#else
      make_tuple(&vp10_highbd_fht4x4_sse4_1, &vp10_highbd_iht4x4_16_add_c, 0,
                 VPX_BITS_10, 16),
      make_tuple(&vp10_highbd_fht4x4_sse4_1, &vp10_highbd_iht4x4_16_add_c, 1,
                 VPX_BITS_10, 16),
      make_tuple(&vp10_highbd_fht4x4_sse4_1, &vp10_highbd_iht4x4_16_add_c, 2,
                 VPX_BITS_10, 16),
      make_tuple(&vp10_highbd_fht4x4_sse4_1, &vp10_highbd_iht4x4_16_add_c, 3,
                 VPX_BITS_10, 16),
      make_tuple(&vp10_highbd_fht4x4_sse4_1, &vp10_highbd_iht4x4_16_add_c, 4,
                 VPX_BITS_10, 16),
      make_tuple(&vp10_highbd_fht4x4_sse4_1, &vp10_highbd_iht4x4_16_add_c, 5,
                 VPX_BITS_10, 16),
      make_tuple(&vp10_highbd_fht4x4_sse4_1, &vp10_highbd_iht4x4_16_add_c, 6,
                 VPX_BITS_10, 16),
      make_tuple(&vp10_highbd_fht4x4_sse4_1, &vp10_highbd_iht4x4_16_add_c, 7,
                 VPX_BITS_10, 16),
      make_tuple(&vp10_highbd_fht4x4_sse4_1, &vp10_highbd_iht4x4_16_add_c, 8,
                 VPX_BITS_10, 16),
      make_tuple(&vp10_highbd_fht4x4_sse4_1, &vp10_highbd_iht4x4_16_add_c, 10,
                 VPX_BITS_10, 16),
      make_tuple(&vp10_highbd_fht4x4_sse4_1, &vp10_highbd_iht4x4_16_add_c, 11,
                 VPX_BITS_10, 16),
      make_tuple(&vp10_highbd_fht4x4_sse4_1, &vp10_highbd_iht4x4_16_add_c, 12,
                 VPX_BITS_10, 16),
      make_tuple(&vp10_highbd_fht4x4_sse4_1, &vp10_highbd_iht4x4_16_add_c, 13,
                 VPX_BITS_10, 16),
      make_tuple(&vp10_highbd_fht4x4_sse4_1, &vp10_highbd_iht4x4_16_add_c, 14,
                 VPX_BITS_10, 16),
      make_tuple(&vp10_highbd_fht4x4_sse4_1, &vp10_highbd_iht4x4_16_add_c, 15,
                 VPX_BITS_10, 16),
      make_tuple(&vp10_highbd_fht4x4_sse4_1, &vp10_highbd_iht4x4_16_add_c, 0,
                 VPX_BITS_12, 16),
      make_tuple(&vp10_highbd_fht4x4_sse4_1, &vp10_highbd_iht4x4_16_add_c, 1,
                 VPX_BITS_12, 16),
      make_tuple(&vp10_highbd_fht4x4_sse4_1, &vp10_highbd_iht4x4_16_add_c, 2,
                 VPX_BITS_12, 16),
      make_tuple(&vp10_highbd_fht4x4_sse4_1, &vp10_highbd_iht4x4_16_add_c, 3,
                 VPX_BITS_12, 16),
      make_tuple(&vp10_highbd_fht4x4_sse4_1, &vp10_highbd_iht4x4_16_add_c, 4,
                 VPX_BITS_12, 16),
      make_tuple(&vp10_highbd_fht4x4_sse4_1, &vp10_highbd_iht4x4_16_add_c, 5,
                 VPX_BITS_12, 16),
      make_tuple(&vp10_highbd_fht4x4_sse4_1, &vp10_highbd_iht4x4_16_add_c, 6,
                 VPX_BITS_12, 16),
      make_tuple(&vp10_highbd_fht4x4_sse4_1, &vp10_highbd_iht4x4_16_add_c, 7,
                 VPX_BITS_12, 16),
      make_tuple(&vp10_highbd_fht4x4_sse4_1, &vp10_highbd_iht4x4_16_add_c, 8,
                 VPX_BITS_12, 16),
      make_tuple(&vp10_highbd_fht4x4_sse4_1, &vp10_highbd_iht4x4_16_add_c, 10,
                 VPX_BITS_12, 16),
      make_tuple(&vp10_highbd_fht4x4_sse4_1, &vp10_highbd_iht4x4_16_add_c, 11,
                 VPX_BITS_12, 16),
      make_tuple(&vp10_highbd_fht4x4_sse4_1, &vp10_highbd_iht4x4_16_add_c, 12,
                 VPX_BITS_12, 16),
      make_tuple(&vp10_highbd_fht4x4_sse4_1, &vp10_highbd_iht4x4_16_add_c, 13,
                 VPX_BITS_12, 16),
      make_tuple(&vp10_highbd_fht4x4_sse4_1, &vp10_highbd_iht4x4_16_add_c, 14,
                 VPX_BITS_12, 16),
      make_tuple(&vp10_highbd_fht4x4_sse4_1, &vp10_highbd_iht4x4_16_add_c, 15,
                 VPX_BITS_12, 16)));
#endif  // !CONFIG_EXT_TX
#endif  // HAVE_SSE4_1 && CONFIG_VP9_HIGHBITDEPTH

}  // namespace
