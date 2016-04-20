#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "./vp10_rtcd.h"
#include "test/acm_random.h"
#include "test/vp10_txfm_test.h"
#include "vp10/common/vp10_fwd_txfm2d_cfg.h"

using libvpx_test::ACMRandom;
using libvpx_test::Fwd_Txfm2d_Func;
using libvpx_test::input_base;
using libvpx_test::bd;

namespace {

#if CONFIG_VP9_HIGHBITDEPTH
TEST(vp10_fwd_txfm2d_sse4_1, accuracy) {
  int16_t input[4096] = {0};
  int32_t output_sse4_1[4096] = {0};
  int32_t output_c[4096] = {0};

  int txfm_num = 17;

  TXFM_2D_CFG cfg_list[] = {
      fwd_txfm_2d_cfg_dct_dct_4,    fwd_txfm_2d_cfg_dct_dct_8,
      fwd_txfm_2d_cfg_dct_dct_16,   fwd_txfm_2d_cfg_dct_dct_32,
      fwd_txfm_2d_cfg_dct_dct_64,   fwd_txfm_2d_cfg_dct_adst_4,
      fwd_txfm_2d_cfg_dct_adst_8,   fwd_txfm_2d_cfg_dct_adst_16,
      fwd_txfm_2d_cfg_dct_adst_32,  fwd_txfm_2d_cfg_adst_dct_4,
      fwd_txfm_2d_cfg_adst_dct_8,   fwd_txfm_2d_cfg_adst_dct_16,
      fwd_txfm_2d_cfg_adst_dct_32,  fwd_txfm_2d_cfg_adst_adst_4,
      fwd_txfm_2d_cfg_adst_adst_8,  fwd_txfm_2d_cfg_adst_adst_16,
      fwd_txfm_2d_cfg_adst_adst_32,
  };

  Fwd_Txfm2d_Func txfm2d_func_c_list[] = {
      vp10_fwd_txfm2d_4x4_c,   vp10_fwd_txfm2d_8x8_c,   vp10_fwd_txfm2d_16x16_c,
      vp10_fwd_txfm2d_32x32_c, vp10_fwd_txfm2d_64x64_c,
  };

  Fwd_Txfm2d_Func txfm2d_func_sse4_1_list[] = {
      vp10_fwd_txfm2d_4x4_sse4_1,   vp10_fwd_txfm2d_8x8_sse4_1,
      vp10_fwd_txfm2d_16x16_sse4_1, vp10_fwd_txfm2d_32x32_sse4_1,
      vp10_fwd_txfm2d_64x64_sse4_1,
  };

  for (int i = 0; i < txfm_num; i++) {
    TXFM_2D_CFG cfg = cfg_list[i];
    int txfm_size = cfg.txfm_size;
    int func_idx = get_max_bit(txfm_size) - 2;
    Fwd_Txfm2d_Func txfm2d_func_c = txfm2d_func_c_list[func_idx];
    Fwd_Txfm2d_Func txfm2d_func_sse4_1 = txfm2d_func_sse4_1_list[func_idx];
    int tx_type = libvpx_test::get_tx_type(&cfg);

    ACMRandom rnd(ACMRandom::DeterministicSeed());

    // init input
    for (int r = 0; r < txfm_size; r++) {
      for (int c = 0; c < txfm_size; c++) {
        input[r * txfm_size + c] = rnd.Rand16() % input_base;
      }
    }

    txfm2d_func_c(input, output_c, cfg.txfm_size, tx_type, bd);
    txfm2d_func_sse4_1(input, output_sse4_1, cfg.txfm_size, tx_type, bd);
    for (int r = 0; r < txfm_size; r++) {
      for (int c = 0; c < txfm_size; c++) {
        EXPECT_EQ(output_c[r * txfm_size + c],
                  output_sse4_1[r * txfm_size + c]);
      }
    }
  }
}
#endif  // CONFIG_VP9_HIGHBITDEPTH
}  // anonymous namespace
