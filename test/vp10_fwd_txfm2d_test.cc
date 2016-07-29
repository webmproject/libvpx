/*
 *  Copyright (c) 2015 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "test/acm_random.h"
#include "test/util.h"
#include "test/vp10_txfm_test.h"
#include "vp10/common/vp10_txfm.h"
#include "./vp10_rtcd.h"

using libvpx_test::ACMRandom;
using libvpx_test::input_base;
using libvpx_test::bd;
using libvpx_test::compute_avg_abs_error;
using libvpx_test::Fwd_Txfm2d_Func;
using libvpx_test::TYPE_TXFM;

namespace {
#if CONFIG_VPX_HIGHBITDEPTH
// tx_type_, tx_size_, max_error_, max_avg_error_
typedef std::tr1::tuple<TX_TYPE, TX_SIZE, double, double> VP10FwdTxfm2dParam;

class VP10FwdTxfm2d : public ::testing::TestWithParam<VP10FwdTxfm2dParam> {
 public:
  virtual void SetUp() {
    tx_type_ = GET_PARAM(0);
    tx_size_ = GET_PARAM(1);
    max_error_ = GET_PARAM(2);
    max_avg_error_ = GET_PARAM(3);
    count_ = 500;
    TXFM_2D_FLIP_CFG fwd_txfm_flip_cfg =
        vp10_get_fwd_txfm_cfg(tx_type_, tx_size_);
    const TXFM_2D_CFG *fwd_txfm_cfg = fwd_txfm_flip_cfg.cfg;
    int amplify_bit = fwd_txfm_cfg->shift[0] + fwd_txfm_cfg->shift[1] +
                      fwd_txfm_cfg->shift[2];
    ud_flip_ = fwd_txfm_flip_cfg.ud_flip;
    lr_flip_ = fwd_txfm_flip_cfg.lr_flip;
    amplify_factor_ =
        amplify_bit >= 0 ? (1 << amplify_bit) : (1.0 / (1 << -amplify_bit));

    fwd_txfm_ = libvpx_test::fwd_txfm_func_ls[tx_size_];
    txfm1d_size_ = libvpx_test::get_txfm1d_size(tx_size_);
    txfm2d_size_ = txfm1d_size_ * txfm1d_size_;
    get_txfm1d_type(tx_type_, &type0_, &type1_);
    input_ = reinterpret_cast<int16_t *>
       (vpx_memalign(16, sizeof(int16_t) * txfm2d_size_));
    output_ = reinterpret_cast<int32_t *>
        (vpx_memalign(16, sizeof(int32_t) * txfm2d_size_));
    ref_input_ = reinterpret_cast<double *>
        (vpx_memalign(16, sizeof(double) * txfm2d_size_));
    ref_output_ = reinterpret_cast<double *>
        (vpx_memalign(16, sizeof(double) * txfm2d_size_));
  }

  void RunFwdAccuracyCheck() {
    ACMRandom rnd(ACMRandom::DeterministicSeed());
    double avg_abs_error = 0;
    for (int ci = 0; ci < count_; ci++) {
      for (int ni = 0; ni < txfm2d_size_; ++ni) {
        input_[ni] = rnd.Rand16() % input_base;
        ref_input_[ni] = static_cast<double>(input_[ni]);
        output_[ni] = 0;
        ref_output_[ni] = 0;
      }

      fwd_txfm_(input_, output_, txfm1d_size_, tx_type_, bd);

      if (lr_flip_ && ud_flip_)
        libvpx_test::fliplrud(ref_input_, txfm1d_size_, txfm1d_size_);
      else if (lr_flip_)
        libvpx_test::fliplr(ref_input_, txfm1d_size_, txfm1d_size_);
      else if (ud_flip_)
        libvpx_test::flipud(ref_input_, txfm1d_size_, txfm1d_size_);

      reference_hybrid_2d(ref_input_, ref_output_, txfm1d_size_,
                          type0_, type1_);

      for (int ni = 0; ni < txfm2d_size_; ++ni) {
        ref_output_[ni] = round(ref_output_[ni] * amplify_factor_);
        EXPECT_GE(max_error_,
                  fabs(output_[ni] - ref_output_[ni]) / amplify_factor_);
      }
      avg_abs_error += compute_avg_abs_error<int32_t, double>(
          output_, ref_output_, txfm2d_size_);
    }

    avg_abs_error /= amplify_factor_;
    avg_abs_error /= count_;
    // max_abs_avg_error comes from upper bound of avg_abs_error
    // printf("type0: %d type1: %d txfm_size: %d accuracy_avg_abs_error:
    // %f\n", type0_, type1_, txfm1d_size_, avg_abs_error);
    EXPECT_GE(max_avg_error_, avg_abs_error);
  }

  virtual void TearDown() {
    vpx_free(input_);
    vpx_free(output_);
    vpx_free(ref_input_);
    vpx_free(ref_output_);
  }

 private:
  double max_error_;
  double max_avg_error_;
  int count_;
  double amplify_factor_;
  TX_TYPE tx_type_;
  TX_SIZE tx_size_;
  int txfm1d_size_;
  int txfm2d_size_;
  Fwd_Txfm2d_Func fwd_txfm_;
  TYPE_TXFM type0_;
  TYPE_TXFM type1_;
  int16_t* input_;
  int32_t* output_;
  double* ref_input_;
  double* ref_output_;
  int ud_flip_;  // flip upside down
  int lr_flip_;  // flip left to right
};

TEST_P(VP10FwdTxfm2d, RunFwdAccuracyCheck) {
  RunFwdAccuracyCheck();
}
const VP10FwdTxfm2dParam vp10_fwd_txfm2d_param_c[] = {
#if CONFIG_EXT_TX
  VP10FwdTxfm2dParam(FLIPADST_DCT,  TX_4X4, 2, 0.2),
  VP10FwdTxfm2dParam(DCT_FLIPADST,  TX_4X4, 2, 0.2),
  VP10FwdTxfm2dParam(FLIPADST_FLIPADST, TX_4X4, 2, 0.2),
  VP10FwdTxfm2dParam(ADST_FLIPADST, TX_4X4, 2, 0.2),
  VP10FwdTxfm2dParam(FLIPADST_ADST, TX_4X4, 2, 0.2),
  VP10FwdTxfm2dParam(FLIPADST_DCT,  TX_8X8, 5, 0.6),
  VP10FwdTxfm2dParam(DCT_FLIPADST,  TX_8X8, 5, 0.6),
  VP10FwdTxfm2dParam(FLIPADST_FLIPADST, TX_8X8, 5, 0.6),
  VP10FwdTxfm2dParam(ADST_FLIPADST, TX_8X8, 5, 0.6),
  VP10FwdTxfm2dParam(FLIPADST_ADST, TX_8X8, 5, 0.6),
  VP10FwdTxfm2dParam(FLIPADST_DCT,  TX_16X16, 11, 1.5),
  VP10FwdTxfm2dParam(DCT_FLIPADST,  TX_16X16, 11, 1.5),
  VP10FwdTxfm2dParam(FLIPADST_FLIPADST, TX_16X16, 11, 1.5),
  VP10FwdTxfm2dParam(ADST_FLIPADST, TX_16X16, 11, 1.5),
  VP10FwdTxfm2dParam(FLIPADST_ADST, TX_16X16, 11, 1.5),
  VP10FwdTxfm2dParam(FLIPADST_DCT,  TX_32X32, 70, 7),
  VP10FwdTxfm2dParam(DCT_FLIPADST,  TX_32X32, 70, 7),
  VP10FwdTxfm2dParam(FLIPADST_FLIPADST, TX_32X32, 70, 7),
  VP10FwdTxfm2dParam(ADST_FLIPADST, TX_32X32, 70, 7),
  VP10FwdTxfm2dParam(FLIPADST_ADST, TX_32X32, 70, 7),
#endif
  VP10FwdTxfm2dParam(DCT_DCT,   TX_4X4, 2, 0.2),
  VP10FwdTxfm2dParam(ADST_DCT,  TX_4X4, 2, 0.2),
  VP10FwdTxfm2dParam(DCT_ADST,  TX_4X4, 2, 0.2),
  VP10FwdTxfm2dParam(ADST_ADST, TX_4X4, 2, 0.2),
  VP10FwdTxfm2dParam(DCT_DCT,   TX_8X8, 5, 0.6),
  VP10FwdTxfm2dParam(ADST_DCT,  TX_8X8, 5, 0.6),
  VP10FwdTxfm2dParam(DCT_ADST,  TX_8X8, 5, 0.6),
  VP10FwdTxfm2dParam(ADST_ADST, TX_8X8, 5, 0.6),
  VP10FwdTxfm2dParam(DCT_DCT,   TX_16X16, 11, 1.5),
  VP10FwdTxfm2dParam(ADST_DCT,  TX_16X16, 11, 1.5),
  VP10FwdTxfm2dParam(DCT_ADST,  TX_16X16, 11, 1.5),
  VP10FwdTxfm2dParam(ADST_ADST, TX_16X16, 11, 1.5),
  VP10FwdTxfm2dParam(DCT_DCT,   TX_32X32, 70, 7),
  VP10FwdTxfm2dParam(ADST_DCT,  TX_32X32, 70, 7),
  VP10FwdTxfm2dParam(DCT_ADST,  TX_32X32, 70, 7),
  VP10FwdTxfm2dParam(ADST_ADST, TX_32X32, 70, 7)
};

INSTANTIATE_TEST_CASE_P(
    C, VP10FwdTxfm2d,
    ::testing::ValuesIn(vp10_fwd_txfm2d_param_c));

#endif  // CONFIG_VPX_HIGHBITDEPTH
}  // namespace
