/*
 *  Copyright (c) 2016 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include <stdlib.h>

#include "third_party/googletest/src/include/gtest/gtest.h"

#include "./vpx_config.h"
#include "./vp10_rtcd.h"
#include "test/acm_random.h"
#include "test/clear_system_state.h"
#include "test/register_state_check.h"
#include "vp10/common/scan.h"

namespace {

typedef void (*QuantizeFpFunc)(const tran_low_t *coeff_ptr, intptr_t count,
                               int skip_block, const int16_t *zbin_ptr,
                               const int16_t *round_ptr,
                               const int16_t *quant_ptr,
                               const int16_t *quant_shift_ptr,
                               tran_low_t *qcoeff_ptr, tran_low_t *dqcoeff_ptr,
                               const int16_t *dequant_ptr, uint16_t *eob_ptr,
                               const int16_t *scan, const int16_t *iscan,
                               const int log_scale);

struct QuantizeFuncParams {
  QuantizeFuncParams(QuantizeFpFunc qF = NULL, QuantizeFpFunc qRefF = NULL,
                     int count = 16) : qFunc(qF), qFuncRef(qRefF),
                                       coeffCount(count) {}
  QuantizeFpFunc qFunc;
  QuantizeFpFunc qFuncRef;
  int coeffCount;
};

using libvpx_test::ACMRandom;

const int numTests = 1000;
const int maxSize = 1024;
const int roundFactorRange = 127;
const int dequantRange = 32768;
const int coeffRange = (1 << 20) - 1;

class VP10QuantizeTest : public ::testing::TestWithParam<QuantizeFuncParams> {
 public:
  void RunQuantizeTest() {
    ACMRandom rnd(ACMRandom::DeterministicSeed());
    DECLARE_ALIGNED(16, tran_low_t, coeff_ptr[maxSize]);
    DECLARE_ALIGNED(16, int16_t, zbin_ptr[2]);
    DECLARE_ALIGNED(16, int16_t, round_ptr[2]);
    DECLARE_ALIGNED(16, int16_t, quant_ptr[2]);
    DECLARE_ALIGNED(16, int16_t, quant_shift_ptr[2]);
    DECLARE_ALIGNED(16, tran_low_t, qcoeff_ptr[maxSize]);
    DECLARE_ALIGNED(16, tran_low_t, dqcoeff_ptr[maxSize]);
    DECLARE_ALIGNED(16, tran_low_t, ref_qcoeff_ptr[maxSize]);
    DECLARE_ALIGNED(16, tran_low_t, ref_dqcoeff_ptr[maxSize]);
    DECLARE_ALIGNED(16, int16_t, dequant_ptr[2]);
    uint16_t eob;
    uint16_t ref_eob;
    int err_count_total = 0;
    int first_failure = -1;
    int skip_block = 0;
    int count = params_.coeffCount;
    const TX_SIZE txSize = getTxSize(count);
    int log_scale = (txSize == TX_32X32);
    QuantizeFpFunc quanFunc = params_.qFunc;
    QuantizeFpFunc quanFuncRef = params_.qFuncRef;

    const scan_order scanOrder = vp10_default_scan_orders[txSize];
    for (int i = 0; i < numTests; i++) {
      int err_count = 0;
      ref_eob = eob = -1;
      for (int j = 0; j < count; j++) {
        coeff_ptr[j] = rnd(coeffRange);
      }

      for (int j = 0; j < 2; j++) {
        zbin_ptr[j] = rnd.Rand16();
        quant_shift_ptr[j] = rnd.Rand16();
        // int16_t positive
        dequant_ptr[j] = abs(rnd(dequantRange));
        quant_ptr[j] = (1 << 16) / dequant_ptr[j];
        round_ptr[j] = (abs(rnd(roundFactorRange)) * dequant_ptr[j]) >> 7;
      }

      quanFuncRef(coeff_ptr, count, skip_block, zbin_ptr,
                  round_ptr, quant_ptr, quant_shift_ptr,
                  ref_qcoeff_ptr, ref_dqcoeff_ptr, dequant_ptr,
                  &ref_eob, scanOrder.scan, scanOrder.iscan,
                  log_scale);

      ASM_REGISTER_STATE_CHECK(quanFunc(coeff_ptr, count, skip_block, zbin_ptr,
                                        round_ptr, quant_ptr, quant_shift_ptr,
                                        qcoeff_ptr, dqcoeff_ptr, dequant_ptr,
                                        &eob, scanOrder.scan, scanOrder.iscan,
                                        log_scale));

      for (int j = 0; j < count; ++j) {
        err_count += (ref_qcoeff_ptr[j]  != qcoeff_ptr[j]) |
            (ref_dqcoeff_ptr[j] != dqcoeff_ptr[j]);
        EXPECT_EQ(ref_qcoeff_ptr[j], qcoeff_ptr[j])
            << "qcoeff error: i = " << i << " j = " << j << "\n";
        EXPECT_EQ(ref_dqcoeff_ptr[j], dqcoeff_ptr[j])
            << "dqcoeff error: i = " << i << " j = " << j << "\n";
      }
      EXPECT_EQ(ref_eob, eob)
          << "eob error: " << "i = " << i << "\n";
      err_count += (ref_eob != eob);
      if (err_count && !err_count_total) {
        first_failure = i;
      }
      err_count_total += err_count;
    }
    EXPECT_EQ(0, err_count_total)
        << "Error: Quantization Test, C output doesn't match SSE2 output. "
        << "First failed at test case " << first_failure;
  }

  void RunEobTest() {
    ACMRandom rnd(ACMRandom::DeterministicSeed());
    DECLARE_ALIGNED(16, tran_low_t, coeff_ptr[maxSize]);
    DECLARE_ALIGNED(16, int16_t, zbin_ptr[2]);
    DECLARE_ALIGNED(16, int16_t, round_ptr[2]);
    DECLARE_ALIGNED(16, int16_t, quant_ptr[2]);
    DECLARE_ALIGNED(16, int16_t, quant_shift_ptr[2]);
    DECLARE_ALIGNED(16, tran_low_t, qcoeff_ptr[maxSize]);
    DECLARE_ALIGNED(16, tran_low_t, dqcoeff_ptr[maxSize]);
    DECLARE_ALIGNED(16, tran_low_t, ref_qcoeff_ptr[maxSize]);
    DECLARE_ALIGNED(16, tran_low_t, ref_dqcoeff_ptr[maxSize]);
    DECLARE_ALIGNED(16, int16_t, dequant_ptr[2]);
    uint16_t eob;
    uint16_t ref_eob;
    int skip_block = 0;
    int count = params_.coeffCount;
    const TX_SIZE txSize = getTxSize(count);
    int log_scale = (txSize == TX_32X32);
    QuantizeFpFunc quanFunc = params_.qFunc;
    QuantizeFpFunc quanFuncRef = params_.qFuncRef;
    const scan_order scanOrder = vp10_default_scan_orders[txSize];

    for (int i = 0; i < numTests; i++) {
      ref_eob = eob = -1;
      for (int j = 0; j < count; j++) {
        coeff_ptr[j] = 0;
      }

      coeff_ptr[rnd(count)] = rnd(coeffRange);
      coeff_ptr[rnd(count)] = rnd(coeffRange);
      coeff_ptr[rnd(count)] = rnd(coeffRange);

      for (int j = 0; j < 2; j++) {
        zbin_ptr[j] = rnd.Rand16();
        quant_shift_ptr[j] = rnd.Rand16();
        // int16_t positive
        dequant_ptr[j] = abs(rnd(dequantRange));
        quant_ptr[j] = (1 << 16) / dequant_ptr[j];
        round_ptr[j] = (abs(rnd(roundFactorRange)) * dequant_ptr[j]) >> 7;
      }

      quanFuncRef(coeff_ptr, count, skip_block, zbin_ptr,
                  round_ptr, quant_ptr, quant_shift_ptr,
                  ref_qcoeff_ptr, ref_dqcoeff_ptr, dequant_ptr,
                  &ref_eob, scanOrder.scan, scanOrder.iscan,
                  log_scale);

      ASM_REGISTER_STATE_CHECK(quanFunc(coeff_ptr, count, skip_block, zbin_ptr,
                                        round_ptr, quant_ptr, quant_shift_ptr,
                                        qcoeff_ptr, dqcoeff_ptr, dequant_ptr,
                                        &eob, scanOrder.scan, scanOrder.iscan,
                                        log_scale));
      EXPECT_EQ(ref_eob, eob)
          << "eob error: " << "i = " << i << "\n";
    }
  }

  virtual void SetUp() {
    params_ = GetParam();
  }

  virtual void TearDown() {
    libvpx_test::ClearSystemState();
  }

  virtual ~VP10QuantizeTest() {}

 private:
  TX_SIZE getTxSize(int count) {
    TX_SIZE txSize = 0;
    if (16 == count) {
      txSize = 0;
    } else if (64 == count) {
      txSize = 1;
    } else if (256 == count) {
      txSize = 2;
    } else if (1024 == count) {
      txSize = 3;
    }
    return txSize;
  }

  QuantizeFuncParams params_;
};

TEST_P(VP10QuantizeTest, BitExactCheck) {
  RunQuantizeTest();
}
TEST_P(VP10QuantizeTest, EobVerify) {
  RunEobTest();
}

#if HAVE_SSE4_1
INSTANTIATE_TEST_CASE_P(
    SSE4_1, VP10QuantizeTest,
    ::testing::Values(QuantizeFuncParams(&vp10_highbd_quantize_fp_sse4_1,
                                         &vp10_highbd_quantize_fp_c, 16),
                      QuantizeFuncParams(&vp10_highbd_quantize_fp_sse4_1,
                                         &vp10_highbd_quantize_fp_c, 64),
                      QuantizeFuncParams(&vp10_highbd_quantize_fp_sse4_1,
                                         &vp10_highbd_quantize_fp_c, 256),
                      QuantizeFuncParams(&vp10_highbd_quantize_fp_sse4_1,
                                         &vp10_highbd_quantize_fp_c, 1024)));
#endif  // HAVE_SSE4_1
}  // namespace
