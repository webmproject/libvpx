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
using std::tr1::tuple;
using libvpx_test::FhtFunc;
typedef tuple<FhtFunc, IhtFunc, int, vpx_bit_depth_t, int> Ht4x4Param;

void fht4x4_ref(const int16_t *in, tran_low_t *out, int stride,
                int tx_type) {
  vp10_fht4x4_c(in, out, stride, tx_type);
}

#if CONFIG_VPX_HIGHBITDEPTH
typedef void (*IhighbdHtFunc)(const tran_low_t *in, uint8_t *out, int stride,
                              int tx_type, int bd);
typedef void (*HBDFhtFunc)(const int16_t *input, int32_t *output, int stride,
                           int tx_type, int bd);

// HighbdHt4x4Param argument list:
// <Target optimized function, tx_type, bit depth>
typedef tuple<HBDFhtFunc, int, int> HighbdHt4x4Param;

void highbe_fht4x4_ref(const int16_t *in, int32_t *out, int stride,
                       int tx_type, int bd) {
  vp10_fwd_txfm2d_4x4_c(in, out, stride, tx_type, bd);
}
#endif  // CONFIG_VPX_HIGHBITDEPTH

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

#if CONFIG_VPX_HIGHBITDEPTH
class VP10HighbdTrans4x4HT : public ::testing::TestWithParam<HighbdHt4x4Param> {
 public:
  virtual ~VP10HighbdTrans4x4HT() {}

  virtual void SetUp() {
    fwd_txfm_ = GET_PARAM(0);
    fwd_txfm_ref_ = highbe_fht4x4_ref;
    tx_type_  = GET_PARAM(1);
    bit_depth_ = GET_PARAM(2);
    mask_ = (1 << bit_depth_) - 1;
    num_coeffs_ = 16;

    input_ = reinterpret_cast<int16_t *>(
        vpx_memalign(16, sizeof(int16_t) * num_coeffs_));
    output_ = reinterpret_cast<int32_t *>(
        vpx_memalign(16, sizeof(int32_t) * num_coeffs_));
    output_ref_ = reinterpret_cast<int32_t *>(
        vpx_memalign(16, sizeof(int32_t) * num_coeffs_));
  }

  virtual void TearDown() {
    vpx_free(input_);
    vpx_free(output_);
    vpx_free(output_ref_);
    libvpx_test::ClearSystemState();
  }

 protected:
  void RunBitexactCheck();

 private:
  HBDFhtFunc fwd_txfm_;
  HBDFhtFunc fwd_txfm_ref_;
  int tx_type_;
  int bit_depth_;
  int mask_;
  int num_coeffs_;
  int16_t *input_;
  int32_t *output_;
  int32_t *output_ref_;
};

void VP10HighbdTrans4x4HT::RunBitexactCheck() {
  ACMRandom rnd(ACMRandom::DeterministicSeed());
  int i, j;
  const int stride = 4;
  const int num_tests = 1000;
  const int num_coeffs = 16;

  for (i = 0; i < num_tests; ++i) {
    for (j = 0; j < num_coeffs; ++j) {
      input_[j] = (rnd.Rand16() & mask_) - (rnd.Rand16() & mask_);
    }

    fwd_txfm_ref_(input_, output_ref_, stride, tx_type_, bit_depth_);
    fwd_txfm_(input_, output_, stride, tx_type_, bit_depth_);

    for (j = 0; j < num_coeffs; ++j) {
      EXPECT_EQ(output_[j], output_ref_[j])
          << "Not bit-exact result at index: " << j
          << " at test block: " << i;
    }
  }
}

TEST_P(VP10HighbdTrans4x4HT, HighbdCoeffCheck) {
  RunBitexactCheck();
}
#endif  // CONFIG_VPX_HIGHBITDEPTH

using std::tr1::make_tuple;

#if HAVE_SSE2
const Ht4x4Param kArrayHt4x4Param_sse2[] = {
      make_tuple(&vp10_fht4x4_sse2, &vp10_iht4x4_16_add_sse2, 0,
                 VPX_BITS_8, 16),
      make_tuple(&vp10_fht4x4_sse2, &vp10_iht4x4_16_add_sse2, 1,
                 VPX_BITS_8, 16),
      make_tuple(&vp10_fht4x4_sse2, &vp10_iht4x4_16_add_sse2, 2,
                 VPX_BITS_8, 16),
      make_tuple(&vp10_fht4x4_sse2, &vp10_iht4x4_16_add_sse2, 3,
                 VPX_BITS_8, 16),
#if CONFIG_EXT_TX
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
                 VPX_BITS_8, 16)
#endif  // CONFIG_EXT_TX
};
INSTANTIATE_TEST_CASE_P(
    SSE2, VP10Trans4x4HT,
    ::testing::ValuesIn(kArrayHt4x4Param_sse2));
#endif  // HAVE_SSE2

#if HAVE_SSE4_1 && CONFIG_VPX_HIGHBITDEPTH
const HighbdHt4x4Param kArrayHighbdHt4x4Param[] = {
         make_tuple(&vp10_fwd_txfm2d_4x4_sse4_1, 0, 10),
         make_tuple(&vp10_fwd_txfm2d_4x4_sse4_1, 0, 12),
         make_tuple(&vp10_fwd_txfm2d_4x4_sse4_1, 1, 10),
         make_tuple(&vp10_fwd_txfm2d_4x4_sse4_1, 1, 12),
         make_tuple(&vp10_fwd_txfm2d_4x4_sse4_1, 2, 10),
         make_tuple(&vp10_fwd_txfm2d_4x4_sse4_1, 2, 12),
         make_tuple(&vp10_fwd_txfm2d_4x4_sse4_1, 3, 10),
         make_tuple(&vp10_fwd_txfm2d_4x4_sse4_1, 3, 12),
#if CONFIG_EXT_TX
         make_tuple(&vp10_fwd_txfm2d_4x4_sse4_1, 4, 10),
         make_tuple(&vp10_fwd_txfm2d_4x4_sse4_1, 4, 12),
         make_tuple(&vp10_fwd_txfm2d_4x4_sse4_1, 5, 10),
         make_tuple(&vp10_fwd_txfm2d_4x4_sse4_1, 5, 12),
         make_tuple(&vp10_fwd_txfm2d_4x4_sse4_1, 6, 10),
         make_tuple(&vp10_fwd_txfm2d_4x4_sse4_1, 6, 12),
         make_tuple(&vp10_fwd_txfm2d_4x4_sse4_1, 7, 10),
         make_tuple(&vp10_fwd_txfm2d_4x4_sse4_1, 7, 12),
         make_tuple(&vp10_fwd_txfm2d_4x4_sse4_1, 8, 10),
         make_tuple(&vp10_fwd_txfm2d_4x4_sse4_1, 8, 12),
#endif  // CONFIG_EXT_TX
};

INSTANTIATE_TEST_CASE_P(
    SSE4_1, VP10HighbdTrans4x4HT,
      ::testing::ValuesIn(kArrayHighbdHt4x4Param));

#endif  // HAVE_SSE4_1 && CONFIG_VPX_HIGHBITDEPTH

}  // namespace
