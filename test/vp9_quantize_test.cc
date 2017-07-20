/*
 *  Copyright (c) 2014 The WebM project authors. All Rights Reserved.
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

#include "./vpx_config.h"
#include "./vpx_dsp_rtcd.h"
#include "test/acm_random.h"
#include "test/buffer.h"
#include "test/clear_system_state.h"
#include "test/register_state_check.h"
#include "test/util.h"
#include "vp9/common/vp9_entropy.h"
#include "vp9/common/vp9_scan.h"
#include "vpx/vpx_codec.h"
#include "vpx/vpx_integer.h"

using libvpx_test::ACMRandom;
using libvpx_test::Buffer;

namespace {
const int number_of_iterations = 100;

typedef void (*QuantizeFunc)(const tran_low_t *coeff, intptr_t count,
                             int skip_block, const int16_t *zbin,
                             const int16_t *round, const int16_t *quant,
                             const int16_t *quant_shift, tran_low_t *qcoeff,
                             tran_low_t *dqcoeff, const int16_t *dequant,
                             uint16_t *eob, const int16_t *scan,
                             const int16_t *iscan);
typedef std::tr1::tuple<QuantizeFunc, QuantizeFunc, vpx_bit_depth_t>
    QuantizeParam;

class VP9QuantizeTest : public ::testing::TestWithParam<QuantizeParam> {
 public:
  virtual ~VP9QuantizeTest() {}
  virtual void SetUp() {
    quantize_op_ = GET_PARAM(0);
    ref_quantize_op_ = GET_PARAM(1);
    bit_depth_ = GET_PARAM(2);
    max_value_ = (1 << bit_depth_) - 1;
  }

  virtual void TearDown() { libvpx_test::ClearSystemState(); }

 protected:
  vpx_bit_depth_t bit_depth_;
  int max_value_;
  QuantizeFunc quantize_op_;
  QuantizeFunc ref_quantize_op_;
};

class VP9Quantize32Test : public ::testing::TestWithParam<QuantizeParam> {
 public:
  virtual ~VP9Quantize32Test() {}
  virtual void SetUp() {
    quantize_op_ = GET_PARAM(0);
    ref_quantize_op_ = GET_PARAM(1);
    bit_depth_ = GET_PARAM(2);
    max_value_ = (1 << bit_depth_) - 1;
  }

  virtual void TearDown() { libvpx_test::ClearSystemState(); }

 protected:
  vpx_bit_depth_t bit_depth_;
  int max_value_;
  QuantizeFunc quantize_op_;
  QuantizeFunc ref_quantize_op_;
};

TEST_P(VP9QuantizeTest, OperationCheck) {
  ACMRandom rnd(ACMRandom::DeterministicSeed());
  Buffer<tran_low_t> coeff = Buffer<tran_low_t>(16, 16, 0, 16);
  ASSERT_TRUE(coeff.Init());
  DECLARE_ALIGNED(16, int16_t, zbin_ptr[8]);
  DECLARE_ALIGNED(16, int16_t, round_ptr[8]);
  DECLARE_ALIGNED(16, int16_t, quant_ptr[8]);
  DECLARE_ALIGNED(16, int16_t, quant_shift_ptr[8]);
  DECLARE_ALIGNED(16, int16_t, dequant_ptr[8]);
  Buffer<tran_low_t> qcoeff = Buffer<tran_low_t>(16, 16, 0, 32);
  ASSERT_TRUE(qcoeff.Init());
  Buffer<tran_low_t> dqcoeff = Buffer<tran_low_t>(16, 16, 0, 32);
  ASSERT_TRUE(dqcoeff.Init());
  Buffer<tran_low_t> ref_qcoeff = Buffer<tran_low_t>(16, 16, 0);
  ASSERT_TRUE(ref_qcoeff.Init());
  Buffer<tran_low_t> ref_dqcoeff = Buffer<tran_low_t>(16, 16, 0);
  ASSERT_TRUE(ref_dqcoeff.Init());
  uint16_t eob, ref_eob;

  for (int i = 0; i < number_of_iterations; ++i) {
    const int skip_block = i == 0;
    const TX_SIZE sz = (TX_SIZE)(i % 3);  // TX_4X4, TX_8X8 TX_16X16
    const TX_TYPE tx_type = (TX_TYPE)((i >> 2) % 3);
    const scan_order *scan_order = &vp9_scan_orders[sz][tx_type];
    const int count = (4 << sz) * (4 << sz);  // 16, 64, 256
    coeff.Set(&rnd, 0, max_value_);
    for (int j = 0; j < 2; j++) {
      // Values determined by deconstructing vp9_init_quantizer().
      // zbin may be up to 1143 for 8 and 10 bit Y values, or 1200 for 12 bit Y
      // values or U/V values of any bit depth. This is because y_delta is not
      // factored into the vp9_ac_quant() call.
      zbin_ptr[j] = rnd.RandRange(1200);
      // round may be up to 685 for Y values or 914 for U/V.
      round_ptr[j] = rnd.RandRange(914);
      // quant ranges from 1 to -32703
      quant_ptr[j] = static_cast<int>(rnd.RandRange(32704)) - 32703;
      // quant_shift goes up to 1 << 16.
      quant_shift_ptr[j] = rnd.RandRange(16384);
      // dequant maxes out at 1828 for all cases.
      dequant_ptr[j] = rnd.RandRange(1828);
    }
    for (int j = 2; j < 8; j++) {
      zbin_ptr[j] = zbin_ptr[1];
      round_ptr[j] = round_ptr[1];
      quant_ptr[j] = quant_ptr[1];
      quant_shift_ptr[j] = quant_shift_ptr[1];
      dequant_ptr[j] = dequant_ptr[1];
    }

    ref_quantize_op_(
        coeff.TopLeftPixel(), count, skip_block, zbin_ptr, round_ptr, quant_ptr,
        quant_shift_ptr, ref_qcoeff.TopLeftPixel(), ref_dqcoeff.TopLeftPixel(),
        dequant_ptr, &ref_eob, scan_order->scan, scan_order->iscan);
    ASM_REGISTER_STATE_CHECK(quantize_op_(
        coeff.TopLeftPixel(), count, skip_block, zbin_ptr, round_ptr, quant_ptr,
        quant_shift_ptr, qcoeff.TopLeftPixel(), dqcoeff.TopLeftPixel(),
        dequant_ptr, &eob, scan_order->scan, scan_order->iscan));

    EXPECT_TRUE(qcoeff.CheckValues(ref_qcoeff));
    EXPECT_TRUE(dqcoeff.CheckValues(ref_dqcoeff));

    EXPECT_EQ(eob, ref_eob);

    if (HasFailure()) {
      printf("Failure on iteration %d.\n", i);
      qcoeff.PrintDifference(ref_qcoeff);
      dqcoeff.PrintDifference(ref_dqcoeff);
      return;
    }
  }
}

TEST_P(VP9Quantize32Test, OperationCheck) {
  ACMRandom rnd(ACMRandom::DeterministicSeed());
  Buffer<tran_low_t> coeff = Buffer<tran_low_t>(32, 32, 0, 16);
  ASSERT_TRUE(coeff.Init());
  DECLARE_ALIGNED(16, int16_t, zbin_ptr[8]);
  DECLARE_ALIGNED(16, int16_t, round_ptr[8]);
  DECLARE_ALIGNED(16, int16_t, quant_ptr[8]);
  DECLARE_ALIGNED(16, int16_t, quant_shift_ptr[8]);
  DECLARE_ALIGNED(16, int16_t, dequant_ptr[8]);
  Buffer<tran_low_t> qcoeff = Buffer<tran_low_t>(32, 32, 0, 32);
  ASSERT_TRUE(qcoeff.Init());
  Buffer<tran_low_t> dqcoeff = Buffer<tran_low_t>(32, 32, 0, 32);
  ASSERT_TRUE(dqcoeff.Init());
  Buffer<tran_low_t> ref_qcoeff = Buffer<tran_low_t>(32, 32, 0);
  ASSERT_TRUE(ref_qcoeff.Init());
  Buffer<tran_low_t> ref_dqcoeff = Buffer<tran_low_t>(32, 32, 0);
  ASSERT_TRUE(ref_dqcoeff.Init());
  uint16_t eob, ref_eob;

  for (int i = 0; i < number_of_iterations; ++i) {
    const int skip_block = i == 0;
    const TX_SIZE sz = TX_32X32;
    const TX_TYPE tx_type = (TX_TYPE)(i % 4);
    const scan_order *scan_order = &vp9_scan_orders[sz][tx_type];
    const int count = (4 << sz) * (4 << sz);  // 1024
    coeff.Set(&rnd, 0, max_value_);
    for (int j = 0; j < 2; j++) {
      zbin_ptr[j] = rnd.RandRange(1200);
      round_ptr[j] = rnd.RandRange(914);
      quant_ptr[j] = static_cast<int>(rnd.RandRange(32704)) - 32703;
      quant_shift_ptr[j] = rnd.RandRange(16384);
      dequant_ptr[j] = rnd.RandRange(1828);
    }
    for (int j = 2; j < 8; j++) {
      zbin_ptr[j] = zbin_ptr[1];
      round_ptr[j] = round_ptr[1];
      quant_ptr[j] = quant_ptr[1];
      quant_shift_ptr[j] = quant_shift_ptr[1];
      dequant_ptr[j] = dequant_ptr[1];
    }

    ref_quantize_op_(
        coeff.TopLeftPixel(), count, skip_block, zbin_ptr, round_ptr, quant_ptr,
        quant_shift_ptr, ref_qcoeff.TopLeftPixel(), ref_dqcoeff.TopLeftPixel(),
        dequant_ptr, &ref_eob, scan_order->scan, scan_order->iscan);
    ASM_REGISTER_STATE_CHECK(quantize_op_(
        coeff.TopLeftPixel(), count, skip_block, zbin_ptr, round_ptr, quant_ptr,
        quant_shift_ptr, qcoeff.TopLeftPixel(), dqcoeff.TopLeftPixel(),
        dequant_ptr, &eob, scan_order->scan, scan_order->iscan));

    EXPECT_TRUE(qcoeff.CheckValues(ref_qcoeff));
    EXPECT_TRUE(dqcoeff.CheckValues(ref_dqcoeff));

    EXPECT_EQ(eob, ref_eob);

    if (HasFailure()) {
      printf("Failure on iteration %d.\n", i);
      qcoeff.PrintDifference(ref_qcoeff);
      dqcoeff.PrintDifference(ref_dqcoeff);
      return;
    }
  }
}

TEST_P(VP9QuantizeTest, EOBCheck) {
  ACMRandom rnd(ACMRandom::DeterministicSeed());
  Buffer<tran_low_t> coeff = Buffer<tran_low_t>(16, 16, 0, 16);
  ASSERT_TRUE(coeff.Init());
  DECLARE_ALIGNED(16, int16_t, zbin_ptr[8]);
  DECLARE_ALIGNED(16, int16_t, round_ptr[8]);
  DECLARE_ALIGNED(16, int16_t, quant_ptr[8]);
  DECLARE_ALIGNED(16, int16_t, quant_shift_ptr[8]);
  DECLARE_ALIGNED(16, int16_t, dequant_ptr[8]);
  Buffer<tran_low_t> qcoeff = Buffer<tran_low_t>(16, 16, 0, 32);
  ASSERT_TRUE(qcoeff.Init());
  Buffer<tran_low_t> dqcoeff = Buffer<tran_low_t>(16, 16, 0, 32);
  ASSERT_TRUE(dqcoeff.Init());
  Buffer<tran_low_t> ref_qcoeff = Buffer<tran_low_t>(16, 16, 0);
  ASSERT_TRUE(ref_qcoeff.Init());
  Buffer<tran_low_t> ref_dqcoeff = Buffer<tran_low_t>(16, 16, 0);
  ASSERT_TRUE(ref_dqcoeff.Init());
  uint16_t eob, ref_eob;

  for (int i = 0; i < number_of_iterations; ++i) {
    int skip_block = i == 0;
    TX_SIZE sz = (TX_SIZE)(i % 3);  // TX_4X4, TX_8X8 TX_16X16
    TX_TYPE tx_type = (TX_TYPE)((i >> 2) % 3);
    const scan_order *scan_order = &vp9_scan_orders[sz][tx_type];
    int count = (4 << sz) * (4 << sz);  // 16, 64, 256
    // Two random entries
    coeff.Set(0);
    coeff.TopLeftPixel()[rnd(count)] = rnd.RandRange(max_value_);
    coeff.TopLeftPixel()[rnd(count)] = rnd.RandRange(max_value_);
    for (int j = 0; j < 2; j++) {
      zbin_ptr[j] = rnd.RandRange(1200);
      round_ptr[j] = rnd.RandRange(914);
      quant_ptr[j] = static_cast<int>(rnd.RandRange(32704)) - 32703;
      quant_shift_ptr[j] = rnd.RandRange(16384);
      dequant_ptr[j] = rnd.RandRange(1828);
    }
    for (int j = 2; j < 8; j++) {
      zbin_ptr[j] = zbin_ptr[1];
      round_ptr[j] = round_ptr[1];
      quant_ptr[j] = quant_ptr[1];
      quant_shift_ptr[j] = quant_shift_ptr[1];
      dequant_ptr[j] = dequant_ptr[1];
    }

    ref_quantize_op_(
        coeff.TopLeftPixel(), count, skip_block, zbin_ptr, round_ptr, quant_ptr,
        quant_shift_ptr, ref_qcoeff.TopLeftPixel(), ref_dqcoeff.TopLeftPixel(),
        dequant_ptr, &ref_eob, scan_order->scan, scan_order->iscan);
    ASM_REGISTER_STATE_CHECK(quantize_op_(
        coeff.TopLeftPixel(), count, skip_block, zbin_ptr, round_ptr, quant_ptr,
        quant_shift_ptr, qcoeff.TopLeftPixel(), dqcoeff.TopLeftPixel(),
        dequant_ptr, &eob, scan_order->scan, scan_order->iscan));

    EXPECT_TRUE(qcoeff.CheckValues(ref_qcoeff));
    EXPECT_TRUE(dqcoeff.CheckValues(ref_dqcoeff));

    EXPECT_EQ(eob, ref_eob);

    if (HasFailure()) {
      printf("Failure on iteration %d.\n", i);
      qcoeff.PrintDifference(ref_qcoeff);
      dqcoeff.PrintDifference(ref_dqcoeff);
      return;
    }
  }
}

TEST_P(VP9Quantize32Test, EOBCheck) {
  ACMRandom rnd(ACMRandom::DeterministicSeed());
  Buffer<tran_low_t> coeff = Buffer<tran_low_t>(32, 32, 0, 16);
  ASSERT_TRUE(coeff.Init());
  DECLARE_ALIGNED(16, int16_t, zbin_ptr[8]);
  DECLARE_ALIGNED(16, int16_t, round_ptr[8]);
  DECLARE_ALIGNED(16, int16_t, quant_ptr[8]);
  DECLARE_ALIGNED(16, int16_t, quant_shift_ptr[8]);
  DECLARE_ALIGNED(16, int16_t, dequant_ptr[8]);
  Buffer<tran_low_t> qcoeff = Buffer<tran_low_t>(32, 32, 0, 32);
  ASSERT_TRUE(qcoeff.Init());
  Buffer<tran_low_t> dqcoeff = Buffer<tran_low_t>(32, 32, 0, 32);
  ASSERT_TRUE(dqcoeff.Init());
  Buffer<tran_low_t> ref_qcoeff = Buffer<tran_low_t>(32, 32, 0);
  ASSERT_TRUE(ref_qcoeff.Init());
  Buffer<tran_low_t> ref_dqcoeff = Buffer<tran_low_t>(32, 32, 0);
  ASSERT_TRUE(ref_dqcoeff.Init());
  uint16_t eob, ref_eob;

  for (int i = 0; i < number_of_iterations; ++i) {
    int skip_block = i == 0;
    TX_SIZE sz = TX_32X32;
    TX_TYPE tx_type = (TX_TYPE)(i % 4);
    const scan_order *scan_order = &vp9_scan_orders[sz][tx_type];
    int count = (4 << sz) * (4 << sz);  // 1024
    coeff.Set(0);
    // Two random entries
    coeff.TopLeftPixel()[rnd(count)] = rnd.RandRange(max_value_);
    coeff.TopLeftPixel()[rnd(count)] = rnd.RandRange(max_value_);
    for (int j = 0; j < 2; j++) {
      zbin_ptr[j] = rnd.RandRange(1200);
      round_ptr[j] = rnd.RandRange(914);
      quant_ptr[j] = static_cast<int>(rnd.RandRange(32704)) - 32703;
      quant_shift_ptr[j] = rnd.RandRange(16384);
      dequant_ptr[j] = rnd.RandRange(1828);
    }
    for (int j = 2; j < 8; j++) {
      zbin_ptr[j] = zbin_ptr[1];
      round_ptr[j] = round_ptr[1];
      quant_ptr[j] = quant_ptr[1];
      quant_shift_ptr[j] = quant_shift_ptr[1];
      dequant_ptr[j] = dequant_ptr[1];
    }

    ref_quantize_op_(
        coeff.TopLeftPixel(), count, skip_block, zbin_ptr, round_ptr, quant_ptr,
        quant_shift_ptr, ref_qcoeff.TopLeftPixel(), ref_dqcoeff.TopLeftPixel(),
        dequant_ptr, &ref_eob, scan_order->scan, scan_order->iscan);
    ASM_REGISTER_STATE_CHECK(quantize_op_(
        coeff.TopLeftPixel(), count, skip_block, zbin_ptr, round_ptr, quant_ptr,
        quant_shift_ptr, qcoeff.TopLeftPixel(), dqcoeff.TopLeftPixel(),
        dequant_ptr, &eob, scan_order->scan, scan_order->iscan));

    EXPECT_TRUE(qcoeff.CheckValues(ref_qcoeff));
    EXPECT_TRUE(dqcoeff.CheckValues(ref_dqcoeff));

    EXPECT_EQ(eob, ref_eob);

    if (HasFailure()) {
      printf("Failure on iteration %d.\n", i);
      qcoeff.PrintDifference(ref_qcoeff);
      dqcoeff.PrintDifference(ref_dqcoeff);
      return;
    }
  }
}
using std::tr1::make_tuple;

#if HAVE_SSE2
#if CONFIG_VP9_HIGHBITDEPTH
// TODO(johannkoenig): Fix vpx_quantize_b_sse2 in highbitdepth builds.
// make_tuple(&vpx_quantize_b_sse2, &vpx_highbd_quantize_b_c, VPX_BITS_8),
INSTANTIATE_TEST_CASE_P(
    SSE2, VP9QuantizeTest,
    ::testing::Values(make_tuple(&vpx_highbd_quantize_b_sse2,
                                 &vpx_highbd_quantize_b_c, VPX_BITS_8),
                      make_tuple(&vpx_highbd_quantize_b_sse2,
                                 &vpx_highbd_quantize_b_c, VPX_BITS_10),
                      make_tuple(&vpx_highbd_quantize_b_sse2,
                                 &vpx_highbd_quantize_b_c, VPX_BITS_12)));
INSTANTIATE_TEST_CASE_P(
    SSE2, VP9Quantize32Test,
    ::testing::Values(make_tuple(&vpx_highbd_quantize_b_32x32_sse2,
                                 &vpx_highbd_quantize_b_32x32_c, VPX_BITS_8),
                      make_tuple(&vpx_highbd_quantize_b_32x32_sse2,
                                 &vpx_highbd_quantize_b_32x32_c, VPX_BITS_10),
                      make_tuple(&vpx_highbd_quantize_b_32x32_sse2,
                                 &vpx_highbd_quantize_b_32x32_c, VPX_BITS_12)));
#else
INSTANTIATE_TEST_CASE_P(SSE2, VP9QuantizeTest,
                        ::testing::Values(make_tuple(&vpx_quantize_b_sse2,
                                                     &vpx_quantize_b_c,
                                                     VPX_BITS_8)));
#endif  // CONFIG_VP9_HIGHBITDEPTH
#endif  // HAVE_SSE2

// TODO(johannkoenig): SSSE3 optimizations do not yet pass these tests.
#if HAVE_SSSE3 && ARCH_X86_64
INSTANTIATE_TEST_CASE_P(DISABLED_SSSE3, VP9QuantizeTest,
                        ::testing::Values(make_tuple(&vpx_quantize_b_ssse3,
                                                     &vpx_quantize_b_c,
                                                     VPX_BITS_8)));

INSTANTIATE_TEST_CASE_P(
    DISABLED_SSSE3, VP9Quantize32Test,
    ::testing::Values(make_tuple(&vpx_quantize_b_32x32_ssse3,
                                 &vpx_quantize_b_32x32_c, VPX_BITS_8)));
#endif  // HAVE_SSSE3 && ARCH_X86_64

// TODO(johannkoenig): AVX optimizations do not yet pass the 32x32 test or
// highbitdepth configurations.
#if HAVE_AVX && ARCH_X86_64 && !CONFIG_VP9_HIGHBITDEPTH
INSTANTIATE_TEST_CASE_P(AVX, VP9QuantizeTest,
                        ::testing::Values(make_tuple(&vpx_quantize_b_avx,
                                                     &vpx_quantize_b_c,
                                                     VPX_BITS_8)));

INSTANTIATE_TEST_CASE_P(DISABLED_AVX, VP9Quantize32Test,
                        ::testing::Values(make_tuple(&vpx_quantize_b_32x32_avx,
                                                     &vpx_quantize_b_32x32_c,
                                                     VPX_BITS_8)));
#endif  // HAVE_AVX && ARCH_X86_64 && !CONFIG_VP9_HIGHBITDEPTH
}  // namespace
