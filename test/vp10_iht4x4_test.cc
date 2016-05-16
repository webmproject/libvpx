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
#include "test/acm_random.h"
#include "test/clear_system_state.h"
#include "test/register_state_check.h"
#include "test/util.h"
#include "vpx_ports/mem.h"

namespace {

using std::tr1::tuple;
using libvpx_test::ACMRandom;

void iht4x4_ref(const int32_t *coeff, uint16_t *output, int stride,
                int tx_type, int bd) {
  vp10_inv_txfm2d_add_4x4_c(coeff, output, stride, tx_type, bd);
}

typedef void (*IHbdHtFunc)(const int32_t *coeff, uint16_t *output, int stride,
                           int tx_type, int bd);

// IhbdHt4x4Param argument list:
//   <target optimization function, tx_type, bit_depth>
typedef tuple<IHbdHtFunc, int, int> IHbdHt4x4Param;

class VP10HighbdInvTrans4x4HT :
      public ::testing::TestWithParam<IHbdHt4x4Param> {
 public:
  virtual ~VP10HighbdInvTrans4x4HT() {}

  virtual void SetUp() {
    inv_txfm_ = GET_PARAM(0);
    inv_txfm_ref_ = iht4x4_ref;
    tx_type_ = GET_PARAM(1);
    bit_depth_ = GET_PARAM(2);
    num_coeffs_ = 4 * 4;

    coeffs_ = reinterpret_cast<int32_t *>(
        vpx_memalign(16, sizeof(int32_t) * num_coeffs_));
    output_ = reinterpret_cast<uint16_t *>(
        vpx_memalign(16, sizeof(uint16_t) * num_coeffs_));
    output_ref_ = reinterpret_cast<uint16_t *>(
        vpx_memalign(16, sizeof(uint16_t) * num_coeffs_));
  }

  virtual void TearDown() {
    vpx_free(coeffs_);
    vpx_free(output_);
    vpx_free(output_ref_);
    libvpx_test::ClearSystemState();
  }

 protected:
  void RunBitexactCheck();

 private:
  IHbdHtFunc inv_txfm_;
  IHbdHtFunc inv_txfm_ref_;
  int tx_type_;
  int bit_depth_;
  int num_coeffs_;
  int32_t *coeffs_;
  uint16_t *output_;
  uint16_t *output_ref_;

  int32_t clamp(int32_t number, int bit) {
    int32_t ret = number;
    const int32_t max = (int32_t)(1 << bit) - 1;
    const int32_t min = -max;

    if (number > max) {
      ret = max;
    } else if (number < min) {
      ret = min;
    }
    return ret;
  }
};

void VP10HighbdInvTrans4x4HT::RunBitexactCheck() {
  ACMRandom rnd(ACMRandom::DeterministicSeed());
  const int stride = 4;
  const int num_tests = 2000000;
  int i;
  int j;
  const uint16_t mask = (1 << bit_depth_) - 1;

  for (i = 0; i < num_tests; ++i) {
    for (j = 0; j < num_coeffs_; ++j) {
      coeffs_[j] = clamp((rnd.Rand16() - rnd.Rand16()) << 2, 18);
      output_ref_[j] = rnd.Rand16() & mask;
      output_[j] = output_ref_[j];
    }

    inv_txfm_ref_(coeffs_, output_ref_, stride, tx_type_, bit_depth_);
    ASM_REGISTER_STATE_CHECK(inv_txfm_(coeffs_, output_, stride, tx_type_,
                                       bit_depth_));

    for (j = 0; j < num_coeffs_; ++j) {
      EXPECT_EQ(output_ref_[j], output_[j])
          << "Not bit-exact result at index: " << j
          << "At test block: " << i;
    }
  }
}

TEST_P(VP10HighbdInvTrans4x4HT, InvTransResultCheck) {
  RunBitexactCheck();
}

using std::tr1::make_tuple;

#if HAVE_SSE4_1 && CONFIG_VP9_HIGHBITDEPTH
const IHbdHt4x4Param kArrayIht4x4Param[] = {
  make_tuple(&vp10_inv_txfm2d_add_4x4_sse4_1, 0, 10),
  make_tuple(&vp10_inv_txfm2d_add_4x4_sse4_1, 0, 12),
  make_tuple(&vp10_inv_txfm2d_add_4x4_sse4_1, 1, 10),
  make_tuple(&vp10_inv_txfm2d_add_4x4_sse4_1, 1, 12),
  make_tuple(&vp10_inv_txfm2d_add_4x4_sse4_1, 2, 10),
  make_tuple(&vp10_inv_txfm2d_add_4x4_sse4_1, 2, 12),
  make_tuple(&vp10_inv_txfm2d_add_4x4_sse4_1, 3, 10),
  make_tuple(&vp10_inv_txfm2d_add_4x4_sse4_1, 3, 12)
};

INSTANTIATE_TEST_CASE_P(
    SSE4_1, VP10HighbdInvTrans4x4HT,
    ::testing::ValuesIn(kArrayIht4x4Param));
#endif  // HAVE_SSE4_1 && CONFIG_VP9_HIGHBITDEPTH

}  // namespace
