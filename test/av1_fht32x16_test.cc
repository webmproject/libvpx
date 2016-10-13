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

#include "./aom_dsp_rtcd.h"
#include "./av1_rtcd.h"

#include "aom_ports/mem.h"
#include "test/acm_random.h"
#include "test/clear_system_state.h"
#include "test/register_state_check.h"
#include "test/transform_test_base.h"
#include "test/util.h"

using libaom_test::ACMRandom;

namespace {
typedef void (*IhtFunc)(const tran_low_t *in, uint8_t *out, int stride,
                        int tx_type);
using std::tr1::tuple;
using libaom_test::FhtFunc;
typedef tuple<FhtFunc, IhtFunc, int, aom_bit_depth_t, int> Ht32x16Param;

void fht32x16_ref(const int16_t *in, tran_low_t *out, int stride, int tx_type) {
  av1_fht32x16_c(in, out, stride, tx_type);
}

void iht32x16_ref(const tran_low_t *in, uint8_t *out, int stride, int tx_type) {
  av1_iht32x16_512_add_c(in, out, stride, tx_type);
}

class AV1Trans32x16HT : public libaom_test::TransformTestBase,
                        public ::testing::TestWithParam<Ht32x16Param> {
 public:
  virtual ~AV1Trans32x16HT() {}

  virtual void SetUp() {
    fwd_txfm_ = GET_PARAM(0);
    inv_txfm_ = GET_PARAM(1);
    tx_type_ = GET_PARAM(2);
    pitch_ = 32;
    height_ = 16;
    fwd_txfm_ref = fht32x16_ref;
    inv_txfm_ref = iht32x16_ref;
    bit_depth_ = GET_PARAM(3);
    mask_ = (1 << bit_depth_) - 1;
    num_coeffs_ = GET_PARAM(4);
  }
  virtual void TearDown() { libaom_test::ClearSystemState(); }

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

TEST_P(AV1Trans32x16HT, CoeffCheck) { RunCoeffCheck(); }
TEST_P(AV1Trans32x16HT, InvCoeffCheck) { RunInvCoeffCheck(); }

using std::tr1::make_tuple;

#if HAVE_SSE2
const Ht32x16Param kArrayHt32x16Param_sse2[] = {
  make_tuple(&av1_fht32x16_sse2, &av1_iht32x16_512_add_sse2, 0, AOM_BITS_8,
             512),
  make_tuple(&av1_fht32x16_sse2, &av1_iht32x16_512_add_sse2, 1, AOM_BITS_8,
             512),
  make_tuple(&av1_fht32x16_sse2, &av1_iht32x16_512_add_sse2, 2, AOM_BITS_8,
             512),
  make_tuple(&av1_fht32x16_sse2, &av1_iht32x16_512_add_sse2, 3, AOM_BITS_8,
             512),
#if CONFIG_EXT_TX
  make_tuple(&av1_fht32x16_sse2, &av1_iht32x16_512_add_sse2, 4, AOM_BITS_8,
             512),
  make_tuple(&av1_fht32x16_sse2, &av1_iht32x16_512_add_sse2, 5, AOM_BITS_8,
             512),
  make_tuple(&av1_fht32x16_sse2, &av1_iht32x16_512_add_sse2, 6, AOM_BITS_8,
             512),
  make_tuple(&av1_fht32x16_sse2, &av1_iht32x16_512_add_sse2, 7, AOM_BITS_8,
             512),
  make_tuple(&av1_fht32x16_sse2, &av1_iht32x16_512_add_sse2, 8, AOM_BITS_8,
             512),
  make_tuple(&av1_fht32x16_sse2, &av1_iht32x16_512_add_sse2, 9, AOM_BITS_8,
             512),
  make_tuple(&av1_fht32x16_sse2, &av1_iht32x16_512_add_sse2, 10, AOM_BITS_8,
             512),
  make_tuple(&av1_fht32x16_sse2, &av1_iht32x16_512_add_sse2, 11, AOM_BITS_8,
             512),
  make_tuple(&av1_fht32x16_sse2, &av1_iht32x16_512_add_sse2, 12, AOM_BITS_8,
             512),
  make_tuple(&av1_fht32x16_sse2, &av1_iht32x16_512_add_sse2, 13, AOM_BITS_8,
             512),
  make_tuple(&av1_fht32x16_sse2, &av1_iht32x16_512_add_sse2, 14, AOM_BITS_8,
             512),
  make_tuple(&av1_fht32x16_sse2, &av1_iht32x16_512_add_sse2, 15, AOM_BITS_8,
             512)
#endif  // CONFIG_EXT_TX
};
INSTANTIATE_TEST_CASE_P(SSE2, AV1Trans32x16HT,
                        ::testing::ValuesIn(kArrayHt32x16Param_sse2));
#endif  // HAVE_SSE2

}  // namespace
