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
#include "vp10/common/enums.h"
#include "vpx_dsp/vpx_dsp_common.h"
#include "vpx_ports/mem.h"

namespace {

using std::tr1::tuple;
using libvpx_test::ACMRandom;

typedef void (*HbdHtFunc)(const int16_t *input, int32_t *output, int stride,
                          int tx_type, int bd);

typedef void (*IHbdHtFunc)(const int32_t *coeff, uint16_t *output, int stride,
                           int tx_type, int bd);

// Test parameter argument list:
//   <transform reference function,
//    optimized inverse transform function,
//    inverse transform reference function,
//    num_coeffs,
//    tx_type,
//    bit_depth>
typedef tuple<HbdHtFunc, IHbdHtFunc, IHbdHtFunc, int, int, int> IHbdHtParam;

class VP10HighbdInvHTNxN : public ::testing::TestWithParam<IHbdHtParam> {
 public:
  virtual ~VP10HighbdInvHTNxN() {}

  virtual void SetUp() {
    txfm_ref_ = GET_PARAM(0);
    inv_txfm_ = GET_PARAM(1);
    inv_txfm_ref_ = GET_PARAM(2);
    num_coeffs_ = GET_PARAM(3);
    tx_type_ = GET_PARAM(4);
    bit_depth_ = GET_PARAM(5);

    input_ = reinterpret_cast<int16_t *>(
        vpx_memalign(16, sizeof(input_[0]) * num_coeffs_));

    // Note:
    // Inverse transform input buffer is 32-byte aligned
    // Refer to <root>/vp10/encoder/context_tree.c, function,
    // void alloc_mode_context().
    coeffs_ = reinterpret_cast<int32_t *>(
        vpx_memalign(32, sizeof(coeffs_[0]) * num_coeffs_));
    output_ = reinterpret_cast<uint16_t *>(
        vpx_memalign(32, sizeof(output_[0]) * num_coeffs_));
    output_ref_ = reinterpret_cast<uint16_t *>(
        vpx_memalign(32, sizeof(output_ref_[0]) * num_coeffs_));
  }

  virtual void TearDown() {
    vpx_free(input_);
    vpx_free(coeffs_);
    vpx_free(output_);
    vpx_free(output_ref_);
    libvpx_test::ClearSystemState();
  }

 protected:
  void RunBitexactCheck();

 private:
  int GetStride() const {
    if (16 == num_coeffs_) {
      return 4;
    } else if (64 == num_coeffs_) {
      return 8;
    } else if (256 == num_coeffs_) {
      return 16;
    } else {
      return 0;
    }
  }

  HbdHtFunc txfm_ref_;
  IHbdHtFunc inv_txfm_;
  IHbdHtFunc inv_txfm_ref_;
  int num_coeffs_;
  int tx_type_;
  int bit_depth_;

  int16_t *input_;
  int32_t *coeffs_;
  uint16_t *output_;
  uint16_t *output_ref_;
};

void VP10HighbdInvHTNxN::RunBitexactCheck() {
  ACMRandom rnd(ACMRandom::DeterministicSeed());
  const int stride = GetStride();
  const int num_tests = 20000;
  const uint16_t mask = (1 << bit_depth_) - 1;

  for (int i = 0; i < num_tests; ++i) {
    for (int j = 0; j < num_coeffs_; ++j) {
      input_[j] = (rnd.Rand16() & mask) - (rnd.Rand16() & mask);
      output_ref_[j] = rnd.Rand16() & mask;
      output_[j] = output_ref_[j];
    }

    txfm_ref_(input_, coeffs_, stride, tx_type_, bit_depth_);
    inv_txfm_ref_(coeffs_, output_ref_, stride, tx_type_, bit_depth_);
    ASM_REGISTER_STATE_CHECK(inv_txfm_(coeffs_, output_, stride, tx_type_,
                                       bit_depth_));

    for (int j = 0; j < num_coeffs_; ++j) {
      EXPECT_EQ(output_ref_[j], output_[j])
          << "Not bit-exact result at index: " << j
          << " At test block: " << i;
    }
  }
}

TEST_P(VP10HighbdInvHTNxN, InvTransResultCheck) {
  RunBitexactCheck();
}

using std::tr1::make_tuple;

#if HAVE_SSE4_1 && CONFIG_VPX_HIGHBITDEPTH
#define PARAM_LIST_4X4 &vp10_fwd_txfm2d_4x4_c, \
             &vp10_inv_txfm2d_add_4x4_sse4_1,  \
             &vp10_inv_txfm2d_add_4x4_c, 16

#define PARAM_LIST_8X8 &vp10_fwd_txfm2d_8x8_c, \
             &vp10_inv_txfm2d_add_8x8_sse4_1,  \
             &vp10_inv_txfm2d_add_8x8_c, 64

#define PARAM_LIST_16X16 &vp10_fwd_txfm2d_16x16_c, \
             &vp10_inv_txfm2d_add_16x16_sse4_1,    \
             &vp10_inv_txfm2d_add_16x16_c, 256

const IHbdHtParam kArrayIhtParam[] = {
  // 16x16
  make_tuple(PARAM_LIST_16X16, DCT_DCT, 10),
  make_tuple(PARAM_LIST_16X16, DCT_DCT, 12),
  make_tuple(PARAM_LIST_16X16, ADST_DCT, 10),
  make_tuple(PARAM_LIST_16X16, ADST_DCT, 12),
  make_tuple(PARAM_LIST_16X16, DCT_ADST, 10),
  make_tuple(PARAM_LIST_16X16, DCT_ADST, 12),
  make_tuple(PARAM_LIST_16X16, ADST_ADST, 10),
  make_tuple(PARAM_LIST_16X16, ADST_ADST, 12),
#if CONFIG_EXT_TX
  make_tuple(PARAM_LIST_16X16, FLIPADST_DCT, 10),
  make_tuple(PARAM_LIST_16X16, FLIPADST_DCT, 12),
  make_tuple(PARAM_LIST_16X16, DCT_FLIPADST, 10),
  make_tuple(PARAM_LIST_16X16, DCT_FLIPADST, 12),
  make_tuple(PARAM_LIST_16X16, FLIPADST_FLIPADST, 10),
  make_tuple(PARAM_LIST_16X16, FLIPADST_FLIPADST, 12),
  make_tuple(PARAM_LIST_16X16, ADST_FLIPADST, 10),
  make_tuple(PARAM_LIST_16X16, ADST_FLIPADST, 12),
  make_tuple(PARAM_LIST_16X16, FLIPADST_ADST, 10),
  make_tuple(PARAM_LIST_16X16, FLIPADST_ADST, 12),
#endif
  // 8x8
  make_tuple(PARAM_LIST_8X8, DCT_DCT, 10),
  make_tuple(PARAM_LIST_8X8, DCT_DCT, 12),
  make_tuple(PARAM_LIST_8X8, ADST_DCT, 10),
  make_tuple(PARAM_LIST_8X8, ADST_DCT, 12),
  make_tuple(PARAM_LIST_8X8, DCT_ADST, 10),
  make_tuple(PARAM_LIST_8X8, DCT_ADST, 12),
  make_tuple(PARAM_LIST_8X8, ADST_ADST, 10),
  make_tuple(PARAM_LIST_8X8, ADST_ADST, 12),
#if CONFIG_EXT_TX
  make_tuple(PARAM_LIST_8X8, FLIPADST_DCT, 10),
  make_tuple(PARAM_LIST_8X8, FLIPADST_DCT, 12),
  make_tuple(PARAM_LIST_8X8, DCT_FLIPADST, 10),
  make_tuple(PARAM_LIST_8X8, DCT_FLIPADST, 12),
  make_tuple(PARAM_LIST_8X8, FLIPADST_FLIPADST, 10),
  make_tuple(PARAM_LIST_8X8, FLIPADST_FLIPADST, 12),
  make_tuple(PARAM_LIST_8X8, ADST_FLIPADST, 10),
  make_tuple(PARAM_LIST_8X8, ADST_FLIPADST, 12),
  make_tuple(PARAM_LIST_8X8, FLIPADST_ADST, 10),
  make_tuple(PARAM_LIST_8X8, FLIPADST_ADST, 12),
#endif
  // 4x4
  make_tuple(PARAM_LIST_4X4, DCT_DCT, 10),
  make_tuple(PARAM_LIST_4X4, DCT_DCT, 12),
  make_tuple(PARAM_LIST_4X4, ADST_DCT, 10),
  make_tuple(PARAM_LIST_4X4, ADST_DCT, 12),
  make_tuple(PARAM_LIST_4X4, DCT_ADST, 10),
  make_tuple(PARAM_LIST_4X4, DCT_ADST, 12),
  make_tuple(PARAM_LIST_4X4, ADST_ADST, 10),
  make_tuple(PARAM_LIST_4X4, ADST_ADST, 12),
#if CONFIG_EXT_TX
  make_tuple(PARAM_LIST_4X4, FLIPADST_DCT, 10),
  make_tuple(PARAM_LIST_4X4, FLIPADST_DCT, 12),
  make_tuple(PARAM_LIST_4X4, DCT_FLIPADST, 10),
  make_tuple(PARAM_LIST_4X4, DCT_FLIPADST, 12),
  make_tuple(PARAM_LIST_4X4, FLIPADST_FLIPADST, 10),
  make_tuple(PARAM_LIST_4X4, FLIPADST_FLIPADST, 12),
  make_tuple(PARAM_LIST_4X4, ADST_FLIPADST, 10),
  make_tuple(PARAM_LIST_4X4, ADST_FLIPADST, 12),
  make_tuple(PARAM_LIST_4X4, FLIPADST_ADST, 10),
  make_tuple(PARAM_LIST_4X4, FLIPADST_ADST, 12),
#endif
};

INSTANTIATE_TEST_CASE_P(
    SSE4_1, VP10HighbdInvHTNxN,
    ::testing::ValuesIn(kArrayIhtParam));
#endif  // HAVE_SSE4_1 && CONFIG_VPX_HIGHBITDEPTH

}  // namespace
