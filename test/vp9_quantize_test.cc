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
#include "./vp9_rtcd.h"
#include "test/acm_random.h"
#include "test/buffer.h"
#include "test/clear_system_state.h"
#include "test/register_state_check.h"
#include "test/util.h"
#include "vp9/common/vp9_entropy.h"
#include "vp9/common/vp9_scan.h"
#include "vpx/vpx_codec.h"
#include "vpx/vpx_integer.h"
#include "vpx_ports/vpx_timer.h"

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
typedef std::tr1::tuple<QuantizeFunc, QuantizeFunc, vpx_bit_depth_t,
                        int /*max_size*/>
    QuantizeParam;

class VP9QuantizeBase {
 public:
  VP9QuantizeBase(vpx_bit_depth_t bit_depth, int max_size)
      : bit_depth_(bit_depth), max_size_(max_size) {
    max_value_ = (1 << bit_depth_) - 1;
    zbin_ptr_ =
        reinterpret_cast<int16_t *>(vpx_memalign(16, 8 * sizeof(*zbin_ptr_)));
    round_ptr_ =
        reinterpret_cast<int16_t *>(vpx_memalign(16, 8 * sizeof(*round_ptr_)));
    quant_ptr_ =
        reinterpret_cast<int16_t *>(vpx_memalign(16, 8 * sizeof(*quant_ptr_)));
    quant_shift_ptr_ = reinterpret_cast<int16_t *>(
        vpx_memalign(16, 8 * sizeof(*quant_shift_ptr_)));
    dequant_ptr_ = reinterpret_cast<int16_t *>(
        vpx_memalign(16, 8 * sizeof(*dequant_ptr_)));
  }

  ~VP9QuantizeBase() {
    vpx_free(zbin_ptr_);
    vpx_free(round_ptr_);
    vpx_free(quant_ptr_);
    vpx_free(quant_shift_ptr_);
    vpx_free(dequant_ptr_);
    zbin_ptr_ = NULL;
    round_ptr_ = NULL;
    quant_ptr_ = NULL;
    quant_shift_ptr_ = NULL;
    dequant_ptr_ = NULL;
    libvpx_test::ClearSystemState();
  }

 protected:
  int16_t *zbin_ptr_;
  int16_t *round_ptr_;
  int16_t *quant_ptr_;
  int16_t *quant_shift_ptr_;
  int16_t *dequant_ptr_;
  const vpx_bit_depth_t bit_depth_;
  int max_value_;
  const int max_size_;
};

class VP9QuantizeTest : public VP9QuantizeBase,
                        public ::testing::TestWithParam<QuantizeParam> {
 public:
  VP9QuantizeTest()
      : VP9QuantizeBase(GET_PARAM(2), GET_PARAM(3)), quantize_op_(GET_PARAM(0)),
        ref_quantize_op_(GET_PARAM(1)) {}

 protected:
  const QuantizeFunc quantize_op_;
  const QuantizeFunc ref_quantize_op_;
};

void GenerateHelperArrays(ACMRandom *rnd, int16_t *zbin, int16_t *round,
                          int16_t *quant, int16_t *quant_shift,
                          int16_t *dequant) {
  for (int j = 0; j < 2; j++) {
    // Values determined by deconstructing vp9_init_quantizer().
    // zbin may be up to 1143 for 8 and 10 bit Y values, or 1200 for 12 bit Y
    // values or U/V values of any bit depth. This is because y_delta is not
    // factored into the vp9_ac_quant() call.
    zbin[j] = rnd->RandRange(1200);
    // round may be up to 685 for Y values or 914 for U/V.
    round[j] = rnd->RandRange(914);
    // quant ranges from 1 to -32703
    quant[j] = static_cast<int>(rnd->RandRange(32704)) - 32703;
    // quant_shift goes up to 1 << 16.
    quant_shift[j] = rnd->RandRange(16384);
    // dequant maxes out at 1828 for all cases.
    dequant[j] = rnd->RandRange(1828);
  }
  for (int j = 2; j < 8; j++) {
    zbin[j] = zbin[1];
    round[j] = round[1];
    quant[j] = quant[1];
    quant_shift[j] = quant_shift[1];
    dequant[j] = dequant[1];
  }
}

TEST_P(VP9QuantizeTest, OperationCheck) {
  ACMRandom rnd(ACMRandom::DeterministicSeed());
  Buffer<tran_low_t> coeff = Buffer<tran_low_t>(max_size_, max_size_, 0, 16);
  ASSERT_TRUE(coeff.Init());
  Buffer<tran_low_t> qcoeff = Buffer<tran_low_t>(max_size_, max_size_, 0, 32);
  ASSERT_TRUE(qcoeff.Init());
  Buffer<tran_low_t> dqcoeff = Buffer<tran_low_t>(max_size_, max_size_, 0, 32);
  ASSERT_TRUE(dqcoeff.Init());
  Buffer<tran_low_t> ref_qcoeff = Buffer<tran_low_t>(max_size_, max_size_, 0);
  ASSERT_TRUE(ref_qcoeff.Init());
  Buffer<tran_low_t> ref_dqcoeff = Buffer<tran_low_t>(max_size_, max_size_, 0);
  ASSERT_TRUE(ref_dqcoeff.Init());
  uint16_t eob, ref_eob;

  for (int i = 0; i < number_of_iterations; ++i) {
    const int skip_block = i == 0;
    TX_SIZE sz;
    if (max_size_ == 16) {
      sz = (TX_SIZE)(i % 3);  // TX_4X4, TX_8X8 TX_16X16
    } else {
      sz = TX_32X32;
    }
    const TX_TYPE tx_type = (TX_TYPE)((i >> 2) % 3);
    const scan_order *scan_order = &vp9_scan_orders[sz][tx_type];
    const int count = (4 << sz) * (4 << sz);  // 16, 64, 256
    coeff.Set(&rnd, 0, max_value_);
    GenerateHelperArrays(&rnd, zbin_ptr_, round_ptr_, quant_ptr_,
                         quant_shift_ptr_, dequant_ptr_);

    ref_quantize_op_(coeff.TopLeftPixel(), count, skip_block, zbin_ptr_,
                     round_ptr_, quant_ptr_, quant_shift_ptr_,
                     ref_qcoeff.TopLeftPixel(), ref_dqcoeff.TopLeftPixel(),
                     dequant_ptr_, &ref_eob, scan_order->scan,
                     scan_order->iscan);
    ASM_REGISTER_STATE_CHECK(
        quantize_op_(coeff.TopLeftPixel(), count, skip_block, zbin_ptr_,
                     round_ptr_, quant_ptr_, quant_shift_ptr_,
                     qcoeff.TopLeftPixel(), dqcoeff.TopLeftPixel(),
                     dequant_ptr_, &eob, scan_order->scan, scan_order->iscan));

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
  Buffer<tran_low_t> coeff = Buffer<tran_low_t>(max_size_, max_size_, 0, 16);
  ASSERT_TRUE(coeff.Init());
  Buffer<tran_low_t> qcoeff = Buffer<tran_low_t>(max_size_, max_size_, 0, 32);
  ASSERT_TRUE(qcoeff.Init());
  Buffer<tran_low_t> dqcoeff = Buffer<tran_low_t>(max_size_, max_size_, 0, 32);
  ASSERT_TRUE(dqcoeff.Init());
  Buffer<tran_low_t> ref_qcoeff = Buffer<tran_low_t>(max_size_, max_size_, 0);
  ASSERT_TRUE(ref_qcoeff.Init());
  Buffer<tran_low_t> ref_dqcoeff = Buffer<tran_low_t>(max_size_, max_size_, 0);
  ASSERT_TRUE(ref_dqcoeff.Init());
  uint16_t eob, ref_eob;

  for (int i = 0; i < number_of_iterations; ++i) {
    int skip_block = i == 0;
    TX_SIZE sz;
    if (max_size_ == 16) {
      sz = (TX_SIZE)(i % 3);  // TX_4X4, TX_8X8 TX_16X16
    } else {
      sz = TX_32X32;
    }
    TX_TYPE tx_type = (TX_TYPE)((i >> 2) % 3);
    const scan_order *scan_order = &vp9_scan_orders[sz][tx_type];
    int count = (4 << sz) * (4 << sz);  // 16, 64, 256
    // Two random entries
    coeff.Set(0);
    coeff.TopLeftPixel()[rnd(count)] = rnd.RandRange(max_value_);
    coeff.TopLeftPixel()[rnd(count)] = rnd.RandRange(max_value_);
    GenerateHelperArrays(&rnd, zbin_ptr_, round_ptr_, quant_ptr_,
                         quant_shift_ptr_, dequant_ptr_);

    ref_quantize_op_(coeff.TopLeftPixel(), count, skip_block, zbin_ptr_,
                     round_ptr_, quant_ptr_, quant_shift_ptr_,
                     ref_qcoeff.TopLeftPixel(), ref_dqcoeff.TopLeftPixel(),
                     dequant_ptr_, &ref_eob, scan_order->scan,
                     scan_order->iscan);
    ASM_REGISTER_STATE_CHECK(
        quantize_op_(coeff.TopLeftPixel(), count, skip_block, zbin_ptr_,
                     round_ptr_, quant_ptr_, quant_shift_ptr_,
                     qcoeff.TopLeftPixel(), dqcoeff.TopLeftPixel(),
                     dequant_ptr_, &eob, scan_order->scan, scan_order->iscan));

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

TEST_P(VP9QuantizeTest, DISABLED_Speed) {
  ACMRandom rnd(ACMRandom::DeterministicSeed());
  Buffer<tran_low_t> coeff = Buffer<tran_low_t>(max_size_, max_size_, 0, 16);
  ASSERT_TRUE(coeff.Init());
  Buffer<tran_low_t> qcoeff = Buffer<tran_low_t>(max_size_, max_size_, 0, 32);
  ASSERT_TRUE(qcoeff.Init());
  Buffer<tran_low_t> dqcoeff = Buffer<tran_low_t>(max_size_, max_size_, 0, 32);
  ASSERT_TRUE(dqcoeff.Init());
  uint16_t eob;
  int starting_sz, ending_sz;

  if (max_size_ == 16) {
    // TX_4X4, TX_8X8 TX_16X16
    starting_sz = 0;
    ending_sz = 2;
  } else {
    // TX_32X32
    starting_sz = 3;
    ending_sz = 3;
  }

  for (TX_SIZE sz = starting_sz; sz <= ending_sz; ++sz) {
    // skip_block, zbin > coeff, zbin < coeff.
    for (int i = 0; i < 3; ++i) {
      const int skip_block = i == 0;
      // TX_TYPE defines the scan order. That is not relevant to the speed test.
      // Pick the first one.
      const TX_TYPE tx_type = DCT_DCT;
      const scan_order *scan_order = &vp9_scan_orders[sz][tx_type];
      const int count = (4 << sz) * (4 << sz);  // 16, 64, 256

      GenerateHelperArrays(&rnd, zbin_ptr_, round_ptr_, quant_ptr_,
                           quant_shift_ptr_, dequant_ptr_);

      if (i == 0) {
        // zbin values are unused when skip_block == 1.
        zbin_ptr_[0] = zbin_ptr_[1] = 0;
        coeff.Set(0);
      } else if (i == 1) {
        // When |coeff values| are less than zbin the results are 0.
        zbin_ptr_[0] = zbin_ptr_[1] = 100;
        coeff.Set(&rnd, -99, 99);
      } else if (i == 2) {
        zbin_ptr_[0] = zbin_ptr_[1] = 50;
        coeff.Set(&rnd, -500, 500);
      }

      vpx_usec_timer timer;
      vpx_usec_timer_start(&timer);
      for (int j = 0; j < 100000000 / count; ++j) {
        quantize_op_(coeff.TopLeftPixel(), count, skip_block, zbin_ptr_,
                     round_ptr_, quant_ptr_, quant_shift_ptr_,
                     qcoeff.TopLeftPixel(), dqcoeff.TopLeftPixel(),
                     dequant_ptr_, &eob, scan_order->scan, scan_order->iscan);
      }
      vpx_usec_timer_mark(&timer);
      const int elapsed_time = static_cast<int>(vpx_usec_timer_elapsed(&timer));
      if (i == 0) printf("Skip block.\n");
      if (i == 1) printf("Bypass calculations.\n");
      if (i == 2) printf("Full calculations.\n");
      printf("Quantize %dx%d time: %5d ms\n", 4 << sz, 4 << sz,
             elapsed_time / 1000);
    }
    printf("\n");
  }
}

using std::tr1::make_tuple;

#if HAVE_SSE2
#if CONFIG_VP9_HIGHBITDEPTH
// TODO(johannkoenig): Fix vpx_quantize_b_sse2 in highbitdepth builds.
// make_tuple(&vpx_quantize_b_sse2, &vpx_highbd_quantize_b_c, VPX_BITS_8),
INSTANTIATE_TEST_CASE_P(
    SSE2, VP9QuantizeTest,
    ::testing::Values(
        make_tuple(&vpx_highbd_quantize_b_sse2, &vpx_highbd_quantize_b_c,
                   VPX_BITS_8, 16),
        make_tuple(&vpx_highbd_quantize_b_sse2, &vpx_highbd_quantize_b_c,
                   VPX_BITS_10, 16),
        make_tuple(&vpx_highbd_quantize_b_sse2, &vpx_highbd_quantize_b_c,
                   VPX_BITS_12, 16),
        make_tuple(&vpx_highbd_quantize_b_32x32_sse2,
                   &vpx_highbd_quantize_b_32x32_c, VPX_BITS_8, 32),
        make_tuple(&vpx_highbd_quantize_b_32x32_sse2,
                   &vpx_highbd_quantize_b_32x32_c, VPX_BITS_10, 32),
        make_tuple(&vpx_highbd_quantize_b_32x32_sse2,
                   &vpx_highbd_quantize_b_32x32_c, VPX_BITS_12, 32)));
#else
INSTANTIATE_TEST_CASE_P(SSE2, VP9QuantizeTest,
                        ::testing::Values(make_tuple(&vpx_quantize_b_sse2,
                                                     &vpx_quantize_b_c,
                                                     VPX_BITS_8, 16)));
#endif  // CONFIG_VP9_HIGHBITDEPTH
#endif  // HAVE_SSE2

// TODO(johannkoenig): SSSE3 optimizations do not yet pass these tests.
#if HAVE_SSSE3 && ARCH_X86_64
INSTANTIATE_TEST_CASE_P(
    DISABLED_SSSE3, VP9QuantizeTest,
    ::testing::Values(make_tuple(&vpx_quantize_b_ssse3, &vpx_quantize_b_c,
                                 VPX_BITS_8, 16),
                      make_tuple(&vpx_quantize_b_32x32_ssse3,
                                 &vpx_quantize_b_32x32_c, VPX_BITS_8, 32)));
#endif  // HAVE_SSSE3 && ARCH_X86_64

// TODO(johannkoenig): AVX optimizations do not yet pass the 32x32 test or
// highbitdepth configurations.
#if HAVE_AVX && ARCH_X86_64 && !CONFIG_VP9_HIGHBITDEPTH
INSTANTIATE_TEST_CASE_P(AVX, VP9QuantizeTest,
                        ::testing::Values(make_tuple(&vpx_quantize_b_avx,
                                                     &vpx_quantize_b_c,
                                                     VPX_BITS_8, 16)));
INSTANTIATE_TEST_CASE_P(DISABLED_AVX, VP9QuantizeTest,
                        ::testing::Values(make_tuple(&vpx_quantize_b_32x32_avx,
                                                     &vpx_quantize_b_32x32_c,
                                                     VPX_BITS_8, 32)));
#endif  // HAVE_AVX && ARCH_X86_64 && !CONFIG_VP9_HIGHBITDEPTH

// TODO(webm:1448): dqcoeff is not handled correctly in HBD builds.
#if HAVE_NEON && !CONFIG_VP9_HIGHBITDEPTH
INSTANTIATE_TEST_CASE_P(NEON, VP9QuantizeTest,
                        ::testing::Values(make_tuple(&vpx_quantize_b_neon,
                                                     &vpx_quantize_b_c,
                                                     VPX_BITS_8, 16)));
#endif  // HAVE_NEON && !CONFIG_VP9_HIGHBITDEPTH

// Only useful to compare "Speed" test results.
INSTANTIATE_TEST_CASE_P(DISABLED_C, VP9QuantizeTest,
                        ::testing::Values(make_tuple(&vpx_quantize_b_c,
                                                     &vpx_quantize_b_c,
                                                     VPX_BITS_8, 16)));
}  // namespace
