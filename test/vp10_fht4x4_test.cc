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
                 VPX_BITS_8, 16)));
#endif  // !CONFIG_EXT_TX
#endif  // HAVE_SSE2

}  // namespace
