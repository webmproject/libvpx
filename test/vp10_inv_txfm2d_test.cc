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

#include "third_party/googletest/src/include/gtest/gtest.h"

#include "test/acm_random.h"
#include "test/vp10_txfm_test.h"
#include "vp10/common/vp10_fwd_txfm2d.h"
#include "vp10/common/vp10_fwd_txfm2d_cfg.h"
#include "vp10/common/vp10_inv_txfm2d.h"
#include "vp10/common/vp10_inv_txfm2d_cfg.h"

using libvpx_test::ACMRandom;

namespace {

const int txfm_size_num = 4;
const int txfm_size_ls[4] = {4, 8, 16, 32};
const TXFM_2D_CFG fwd_txfm_cfg_ls[4][4] = {
    {fwd_txfm_2d_cfg_dct_dct_4, fwd_txfm_2d_cfg_dct_adst_4,
     fwd_txfm_2d_cfg_adst_adst_4, fwd_txfm_2d_cfg_adst_dct_4},
    {fwd_txfm_2d_cfg_dct_dct_8, fwd_txfm_2d_cfg_dct_adst_8,
     fwd_txfm_2d_cfg_adst_adst_8, fwd_txfm_2d_cfg_adst_dct_8},
    {fwd_txfm_2d_cfg_dct_dct_16, fwd_txfm_2d_cfg_dct_adst_16,
     fwd_txfm_2d_cfg_adst_adst_16, fwd_txfm_2d_cfg_adst_dct_16},
    {fwd_txfm_2d_cfg_dct_dct_32, fwd_txfm_2d_cfg_dct_adst_32,
     fwd_txfm_2d_cfg_adst_adst_32, fwd_txfm_2d_cfg_adst_dct_32}};

const TXFM_2D_CFG inv_txfm_cfg_ls[4][4] = {
    {inv_txfm_2d_cfg_dct_dct_4, inv_txfm_2d_cfg_dct_adst_4,
     inv_txfm_2d_cfg_adst_adst_4, inv_txfm_2d_cfg_adst_dct_4},
    {inv_txfm_2d_cfg_dct_dct_8, inv_txfm_2d_cfg_dct_adst_8,
     inv_txfm_2d_cfg_adst_adst_8, inv_txfm_2d_cfg_adst_dct_8},
    {inv_txfm_2d_cfg_dct_dct_16, inv_txfm_2d_cfg_dct_adst_16,
     inv_txfm_2d_cfg_adst_adst_16, inv_txfm_2d_cfg_adst_dct_16},
    {inv_txfm_2d_cfg_dct_dct_32, inv_txfm_2d_cfg_dct_adst_32,
     inv_txfm_2d_cfg_adst_adst_32, inv_txfm_2d_cfg_adst_dct_32}};

const Fwd_Txfm2d_Func fwd_txfm_func_ls[4] = {
    vp10_fwd_txfm2d_4x4, vp10_fwd_txfm2d_8x8, vp10_fwd_txfm2d_16x16,
    vp10_fwd_txfm2d_32x32};
const Inv_Txfm2d_Func inv_txfm_func_ls[4] = {
    vp10_inv_txfm2d_add_4x4, vp10_inv_txfm2d_add_8x8, vp10_inv_txfm2d_add_16x16,
    vp10_inv_txfm2d_add_32x32};

const int txfm_type_num = 4;

TEST(vp10_inv_txfm2d, round_trip) {
  for (int txfm_size_idx = 0; txfm_size_idx < txfm_size_num; ++txfm_size_idx) {
    const int txfm_size = txfm_size_ls[txfm_size_idx];
    const int sqr_txfm_size = txfm_size * txfm_size;
    int16_t* input = new int16_t[sqr_txfm_size];
    uint16_t* ref_input = new uint16_t[sqr_txfm_size];
    int32_t* output = new int32_t[sqr_txfm_size];

    for (int txfm_type_idx = 0; txfm_type_idx < txfm_type_num;
         ++txfm_type_idx) {
      const TXFM_2D_CFG fwd_txfm_cfg =
          fwd_txfm_cfg_ls[txfm_size_idx][txfm_type_idx];
      const TXFM_2D_CFG inv_txfm_cfg =
          inv_txfm_cfg_ls[txfm_size_idx][txfm_type_idx];
      const Fwd_Txfm2d_Func fwd_txfm_func = fwd_txfm_func_ls[txfm_size_idx];
      const Inv_Txfm2d_Func inv_txfm_func = inv_txfm_func_ls[txfm_size_idx];
      const int count = 5000;
      double avg_abs_error = 0;
      ACMRandom rnd(ACMRandom::DeterministicSeed());
      for (int ci = 0; ci < count; ci++) {
        for (int ni = 0; ni < sqr_txfm_size; ++ni) {
          if (ci == 0) {
            int extreme_input = base - 1;
            input[ni] = extreme_input;  // extreme case
            ref_input[ni] = 0;
          } else {
            input[ni] = rnd.Rand16() % base;
            ref_input[ni] = 0;
          }
        }

        fwd_txfm_func(input, output, txfm_size, &fwd_txfm_cfg, bd);
        inv_txfm_func(output, ref_input, txfm_size, &inv_txfm_cfg, bd);

        for (int ni = 0; ni < sqr_txfm_size; ++ni) {
          EXPECT_LE(abs(input[ni] - ref_input[ni]), 2);
        }
        avg_abs_error += compute_avg_abs_error<int16_t, uint16_t>(
            input, ref_input, sqr_txfm_size);
      }

      avg_abs_error /= count;
      // max_abs_avg_error comes from upper bound of
      // printf("txfm_size: %d accuracy_avg_abs_error: %f\n", txfm_size,
      // avg_abs_error);
      // TODO(angiebird): this upper bound is from adst_adst_8
      const double max_abs_avg_error = 0.024;
      EXPECT_LE(avg_abs_error, max_abs_avg_error);
    }

    delete[] input;
    delete[] ref_input;
    delete[] output;
  }
}

}  // anonymous namespace
