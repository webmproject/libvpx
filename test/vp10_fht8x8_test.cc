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
typedef std::tr1::tuple<FhtFunc, IhtFunc, int, vpx_bit_depth_t, int> Ht8x8Param;

void fht8x8_ref(const int16_t *in, tran_low_t *out, int stride,
                int tx_type) {
  vp10_fht8x8_c(in, out, stride, tx_type);
}

class VP10Trans8x8HT
    : public libvpx_test::TransformTestBase,
      public ::testing::TestWithParam<Ht8x8Param> {
 public:
  virtual ~VP10Trans8x8HT() {}

  virtual void SetUp() {
    fwd_txfm_ = GET_PARAM(0);
    inv_txfm_ = GET_PARAM(1);
    tx_type_  = GET_PARAM(2);
    pitch_    = 8;
    fwd_txfm_ref = fht8x8_ref;
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

TEST_P(VP10Trans8x8HT, CoeffCheck) {
  RunCoeffCheck();
}

#if CONFIG_EXT_TX && !CONFIG_VP9_HIGHBITDEPTH
TEST(VP10Trans8x8HTSpeedTest, C_version) {
    ACMRandom rnd(ACMRandom::DeterministicSeed());
    const int count_test_block = 20000;
    int bit_depth = 8;
    int mask = (1 << bit_depth) - 1;
    const int num_coeffs = 64;
    int16_t *input = new int16_t[num_coeffs];
    tran_low_t *output = new tran_low_t[num_coeffs];
    const int stride = 8;
    int tx_type;

    for (int i = 0; i < count_test_block; ++i) {
      for (int j = 0; j < num_coeffs; ++j) {
        input[j] = (rnd.Rand8() & mask) - (rnd.Rand8() & mask);
      }
      for (tx_type = V_DCT; tx_type <= H_FLIPADST; ++tx_type) {
        vp10_fht8x8_c(input, output, stride, tx_type);
      }
    }

    delete[] input;
    delete[] output;
}
#endif  // CONFIG_EXT_TX && !CONFIG_VP9_HIGHBITDEPTH

#if HAVE_SSE2 && CONFIG_EXT_TX && !CONFIG_VP9_HIGHBITDEPTH
TEST(VP10Trans8x8HTSpeedTest, SSE2_version) {
    ACMRandom rnd(ACMRandom::DeterministicSeed());
    const int count_test_block = 20000;
    int bit_depth = 8;
    int mask = (1 << bit_depth) - 1;
    const int num_coeffs = 64;
    int16_t *input = reinterpret_cast<int16_t *>
        (vpx_memalign(16, sizeof(int16_t) * num_coeffs));
    tran_low_t *output = reinterpret_cast<tran_low_t *>
        (vpx_memalign(16, sizeof(tran_low_t) * num_coeffs));

    const int stride = 8;
    int tx_type;

    for (int i = 0; i < count_test_block; ++i) {
      for (int j = 0; j < num_coeffs; ++j) {
        input[j] = (rnd.Rand8() & mask) - (rnd.Rand8() & mask);
      }
      for (tx_type = V_DCT; tx_type <= H_FLIPADST; ++tx_type) {
        vp10_fht8x8_sse2(input, output, stride, tx_type);
      }
    }

    vpx_free(input);
    vpx_free(output);
}
#endif  // HAVE_SSE2 && CONFIG_EXT_TX && !CONFIG_VP9_HIGHBITDEPTH

using std::tr1::make_tuple;

#if HAVE_SSE2
const Ht8x8Param kArrayHt8x8Param_sse2[] = {
      make_tuple(&vp10_fht8x8_sse2, &vp10_iht8x8_64_add_sse2, 0,
                 VPX_BITS_8, 64),
      make_tuple(&vp10_fht8x8_sse2, &vp10_iht8x8_64_add_sse2, 1,
                 VPX_BITS_8, 64),
      make_tuple(&vp10_fht8x8_sse2, &vp10_iht8x8_64_add_sse2, 2,
                 VPX_BITS_8, 64),
      make_tuple(&vp10_fht8x8_sse2, &vp10_iht8x8_64_add_sse2, 3,
                 VPX_BITS_8, 64),
#if CONFIG_EXT_TX
      make_tuple(&vp10_fht8x8_sse2, &vp10_iht8x8_64_add_sse2, 4,
                 VPX_BITS_8, 64),
      make_tuple(&vp10_fht8x8_sse2, &vp10_iht8x8_64_add_sse2, 5,
                 VPX_BITS_8, 64),
      make_tuple(&vp10_fht8x8_sse2, &vp10_iht8x8_64_add_sse2, 6,
                 VPX_BITS_8, 64),
      make_tuple(&vp10_fht8x8_sse2, &vp10_iht8x8_64_add_sse2, 7,
                 VPX_BITS_8, 64),
      make_tuple(&vp10_fht8x8_sse2, &vp10_iht8x8_64_add_sse2, 8,
                 VPX_BITS_8, 64),
      make_tuple(&vp10_fht8x8_sse2, &vp10_iht8x8_64_add_sse2, 10,
                 VPX_BITS_8, 64),
      make_tuple(&vp10_fht8x8_sse2, &vp10_iht8x8_64_add_sse2, 11,
                 VPX_BITS_8, 64),
      make_tuple(&vp10_fht8x8_sse2, &vp10_iht8x8_64_add_sse2, 12,
                 VPX_BITS_8, 64),
      make_tuple(&vp10_fht8x8_sse2, &vp10_iht8x8_64_add_sse2, 13,
                 VPX_BITS_8, 64),
      make_tuple(&vp10_fht8x8_sse2, &vp10_iht8x8_64_add_sse2, 14,
                 VPX_BITS_8, 64),
      make_tuple(&vp10_fht8x8_sse2, &vp10_iht8x8_64_add_sse2, 15,
                 VPX_BITS_8, 64)
#endif  // CONFIG_EXT_TX
};
INSTANTIATE_TEST_CASE_P(
    SSE2, VP10Trans8x8HT,
    ::testing::ValuesIn(kArrayHt8x8Param_sse2));
#endif  // HAVE_SSE2

}  // namespace
