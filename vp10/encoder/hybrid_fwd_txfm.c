/*
 *  Copyright (c) 2015 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include "./vp10_rtcd.h"
#include "./vpx_config.h"
#include "./vpx_dsp_rtcd.h"

#include "vp10/common/idct.h"
#include "vp10/encoder/hybrid_fwd_txfm.h"

static INLINE void fdct32x32(int rd_transform, const int16_t *src,
                             tran_low_t *dst, int src_stride) {
  if (rd_transform)
    vpx_fdct32x32_rd(src, dst, src_stride);
  else
    vpx_fdct32x32(src, dst, src_stride);
}

static void fwd_txfm_4x4(const int16_t *src_diff, tran_low_t *coeff,
                         int diff_stride, TX_TYPE tx_type, int lossless) {
  if (lossless) {
    assert(tx_type == DCT_DCT);
    vp10_fwht4x4(src_diff, coeff, diff_stride);
    return;
  }

  switch (tx_type) {
    case DCT_DCT:
    case ADST_DCT:
    case DCT_ADST:
    case ADST_ADST:
      vp10_fht4x4(src_diff, coeff, diff_stride, tx_type);
      break;
#if CONFIG_EXT_TX
    case FLIPADST_DCT:
    case DCT_FLIPADST:
    case FLIPADST_FLIPADST:
    case ADST_FLIPADST:
    case FLIPADST_ADST:
    case V_DCT:
    case H_DCT:
    case V_ADST:
    case H_ADST:
    case V_FLIPADST:
    case H_FLIPADST:
      vp10_fht4x4(src_diff, coeff, diff_stride, tx_type);
      break;
    case IDTX:
      vp10_fwd_idtx_c(src_diff, coeff, diff_stride, 4, tx_type);
      break;
#endif  // CONFIG_EXT_TX
    default:
      assert(0);
  }
}

#if CONFIG_EXT_TX
static void fwd_txfm_8x4(const int16_t *src_diff, tran_low_t *coeff,
                         int diff_stride, TX_TYPE tx_type,
                         FWD_TXFM_OPT fwd_txfm_opt) {
  (void) fwd_txfm_opt;
  vp10_fht8x4(src_diff, coeff, diff_stride, tx_type);
}

static void fwd_txfm_4x8(const int16_t *src_diff, tran_low_t *coeff,
                         int diff_stride, TX_TYPE tx_type,
                         FWD_TXFM_OPT fwd_txfm_opt) {
  (void) fwd_txfm_opt;
  vp10_fht4x8(src_diff, coeff, diff_stride, tx_type);
}
#endif  // CONFIG_EXT_TX

static void fwd_txfm_8x8(const int16_t *src_diff, tran_low_t *coeff,
                         int diff_stride, TX_TYPE tx_type,
                         FWD_TXFM_OPT fwd_txfm_opt) {
  switch (tx_type) {
    case DCT_DCT:
    case ADST_DCT:
    case DCT_ADST:
    case ADST_ADST:
      if (fwd_txfm_opt == FWD_TXFM_OPT_NORMAL)
        vp10_fht8x8(src_diff, coeff, diff_stride, tx_type);
      else  // FWD_TXFM_OPT_DC
        vpx_fdct8x8_1(src_diff, coeff, diff_stride);
      break;
#if CONFIG_EXT_TX
    case FLIPADST_DCT:
    case DCT_FLIPADST:
    case FLIPADST_FLIPADST:
    case ADST_FLIPADST:
    case FLIPADST_ADST:
    case V_DCT:
    case H_DCT:
    case V_ADST:
    case H_ADST:
    case V_FLIPADST:
    case H_FLIPADST:
      vp10_fht8x8(src_diff, coeff, diff_stride, tx_type);
      break;
    case IDTX:
      vp10_fwd_idtx_c(src_diff, coeff, diff_stride, 8, tx_type);
      break;
#endif  // CONFIG_EXT_TX
    default:
      assert(0);
  }
}

static void fwd_txfm_16x16(const int16_t *src_diff, tran_low_t *coeff,
                           int diff_stride, TX_TYPE tx_type,
                           FWD_TXFM_OPT fwd_txfm_opt) {
  switch (tx_type) {
    case DCT_DCT:
    case ADST_DCT:
    case DCT_ADST:
    case ADST_ADST:
      if (fwd_txfm_opt == FWD_TXFM_OPT_NORMAL)
        vp10_fht16x16(src_diff, coeff, diff_stride, tx_type);
      else  // FWD_TXFM_OPT_DC
        vpx_fdct16x16_1(src_diff, coeff, diff_stride);
      break;
#if CONFIG_EXT_TX
    case FLIPADST_DCT:
    case DCT_FLIPADST:
    case FLIPADST_FLIPADST:
    case ADST_FLIPADST:
    case FLIPADST_ADST:
    case V_DCT:
    case H_DCT:
    case V_ADST:
    case H_ADST:
    case V_FLIPADST:
    case H_FLIPADST:
      vp10_fht16x16(src_diff, coeff, diff_stride, tx_type);
      break;
    case IDTX:
      vp10_fwd_idtx_c(src_diff, coeff, diff_stride, 16, tx_type);
      break;
#endif  // CONFIG_EXT_TX
    default:
      assert(0);
  }
}

static void fwd_txfm_32x32(int rd_transform, const int16_t *src_diff,
                           tran_low_t *coeff, int diff_stride, TX_TYPE tx_type,
                           FWD_TXFM_OPT fwd_txfm_opt) {
  switch (tx_type) {
    case DCT_DCT:
      if (fwd_txfm_opt == FWD_TXFM_OPT_NORMAL)
        fdct32x32(rd_transform, src_diff, coeff, diff_stride);
      else  // FWD_TXFM_OPT_DC
        vpx_fdct32x32_1(src_diff, coeff, diff_stride);
      break;
#if CONFIG_EXT_TX
    case ADST_DCT:
    case DCT_ADST:
    case ADST_ADST:
    case FLIPADST_DCT:
    case DCT_FLIPADST:
    case FLIPADST_FLIPADST:
    case ADST_FLIPADST:
    case FLIPADST_ADST:
      vp10_fht32x32_c(src_diff, coeff, diff_stride, tx_type);
      break;
    case V_DCT:
    case H_DCT:
    case V_ADST:
    case H_ADST:
    case V_FLIPADST:
    case H_FLIPADST:
      vp10_fht32x32_c(src_diff, coeff, diff_stride, tx_type);
      break;
    case IDTX:
      vp10_fwd_idtx_c(src_diff, coeff, diff_stride, 32, tx_type);
      break;
#endif  // CONFIG_EXT_TX
    default:
      assert(0);
      break;
  }
}

#if CONFIG_VP9_HIGHBITDEPTH
static void highbd_fwd_txfm_4x4(const int16_t *src_diff, tran_low_t *coeff,
                                int diff_stride, TX_TYPE tx_type, int lossless,
                                const int bd) {
  if (lossless) {
    assert(tx_type == DCT_DCT);
    vp10_highbd_fwht4x4(src_diff, coeff, diff_stride);
    return;
  }

  switch (tx_type) {
    case DCT_DCT:
    case ADST_DCT:
    case DCT_ADST:
    case ADST_ADST:
      vp10_fwd_txfm2d_4x4(src_diff, coeff, diff_stride, tx_type, bd);
      break;
#if CONFIG_EXT_TX
    case FLIPADST_DCT:
    case DCT_FLIPADST:
    case FLIPADST_FLIPADST:
    case ADST_FLIPADST:
    case FLIPADST_ADST:
      vp10_fwd_txfm2d_4x4(src_diff, coeff, diff_stride, tx_type, bd);
      break;
    case V_DCT:
    case H_DCT:
    case V_ADST:
    case H_ADST:
    case V_FLIPADST:
    case H_FLIPADST:
      vp10_highbd_fht4x4_c(src_diff, coeff, diff_stride, tx_type);
      break;
    case IDTX:
      vp10_fwd_idtx_c(src_diff, coeff, diff_stride, 4, tx_type);
      break;
#endif  // CONFIG_EXT_TX
    default:
      assert(0);
  }
}

#if CONFIG_EXT_TX
static void highbd_fwd_txfm_8x4(const int16_t *src_diff, tran_low_t *coeff,
                                int diff_stride, TX_TYPE tx_type,
                                FWD_TXFM_OPT fwd_txfm_opt, const int bd) {
  (void) fwd_txfm_opt;
  (void) bd;
  vp10_highbd_fht8x4(src_diff, coeff, diff_stride, tx_type);
}

static void highbd_fwd_txfm_4x8(const int16_t *src_diff, tran_low_t *coeff,
                                int diff_stride, TX_TYPE tx_type,
                                FWD_TXFM_OPT fwd_txfm_opt, const int bd) {
  (void) fwd_txfm_opt;
  (void) bd;
  vp10_highbd_fht4x8(src_diff, coeff, diff_stride, tx_type);
}
#endif  // CONFIG_EXT_TX

static void highbd_fwd_txfm_8x8(const int16_t *src_diff, tran_low_t *coeff,
                                int diff_stride, TX_TYPE tx_type,
                                FWD_TXFM_OPT fwd_txfm_opt, const int bd) {
  (void)fwd_txfm_opt;
  switch (tx_type) {
    case DCT_DCT:
    case ADST_DCT:
    case DCT_ADST:
    case ADST_ADST:
      vp10_fwd_txfm2d_8x8(src_diff, coeff, diff_stride, tx_type, bd);
      break;
#if CONFIG_EXT_TX
    case FLIPADST_DCT:
    case DCT_FLIPADST:
    case FLIPADST_FLIPADST:
    case ADST_FLIPADST:
    case FLIPADST_ADST:
      vp10_fwd_txfm2d_8x8(src_diff, coeff, diff_stride, tx_type, bd);
      break;
    case V_DCT:
    case H_DCT:
    case V_ADST:
    case H_ADST:
    case V_FLIPADST:
    case H_FLIPADST:
      // Use C version since DST exists only in C
      vp10_highbd_fht8x8_c(src_diff, coeff, diff_stride, tx_type);
      break;
    case IDTX:
      vp10_fwd_idtx_c(src_diff, coeff, diff_stride, 8, tx_type);
      break;
#endif  // CONFIG_EXT_TX
    default:
      assert(0);
  }
}

static void highbd_fwd_txfm_16x16(const int16_t *src_diff, tran_low_t *coeff,
                                  int diff_stride, TX_TYPE tx_type,
                                  FWD_TXFM_OPT fwd_txfm_opt, const int bd) {
  (void)fwd_txfm_opt;
  switch (tx_type) {
    case DCT_DCT:
    case ADST_DCT:
    case DCT_ADST:
    case ADST_ADST:
      vp10_fwd_txfm2d_16x16(src_diff, coeff, diff_stride, tx_type, bd);
      break;
#if CONFIG_EXT_TX
    case FLIPADST_DCT:
    case DCT_FLIPADST:
    case FLIPADST_FLIPADST:
    case ADST_FLIPADST:
    case FLIPADST_ADST:
      vp10_fwd_txfm2d_16x16(src_diff, coeff, diff_stride, tx_type, bd);
      break;
    case V_DCT:
    case H_DCT:
    case V_ADST:
    case H_ADST:
    case V_FLIPADST:
    case H_FLIPADST:
      // Use C version since DST exists only in C
      vp10_highbd_fht16x16_c(src_diff, coeff, diff_stride, tx_type);
      break;
    case IDTX:
      vp10_fwd_idtx_c(src_diff, coeff, diff_stride, 16, tx_type);
      break;
#endif  // CONFIG_EXT_TX
    default:
      assert(0);
  }
}

static void highbd_fwd_txfm_32x32(int rd_transform, const int16_t *src_diff,
                                  tran_low_t *coeff, int diff_stride,
                                  TX_TYPE tx_type, FWD_TXFM_OPT fwd_txfm_opt,
                                  const int bd) {
  (void)rd_transform;
  (void)fwd_txfm_opt;
  switch (tx_type) {
    case DCT_DCT:
      vp10_fwd_txfm2d_32x32(src_diff, coeff, diff_stride, tx_type, bd);
      break;
#if CONFIG_EXT_TX
    case ADST_DCT:
    case DCT_ADST:
    case ADST_ADST:
    case FLIPADST_DCT:
    case DCT_FLIPADST:
    case FLIPADST_FLIPADST:
    case ADST_FLIPADST:
    case FLIPADST_ADST:
    case V_DCT:
    case H_DCT:
    case V_ADST:
    case H_ADST:
    case V_FLIPADST:
    case H_FLIPADST:
      vp10_highbd_fht32x32_c(src_diff, coeff, diff_stride, tx_type);
      break;
    case IDTX:
      vp10_fwd_idtx_c(src_diff, coeff, diff_stride, 32, tx_type);
      break;
#endif  // CONFIG_EXT_TX
    default:
      assert(0);
      break;
  }
}
#endif  // CONFIG_VP9_HIGHBITDEPTH

void fwd_txfm(const int16_t *src_diff, tran_low_t *coeff, int diff_stride,
              FWD_TXFM_PARAM *fwd_txfm_param) {
  const int fwd_txfm_opt = fwd_txfm_param->fwd_txfm_opt;
  const TX_TYPE tx_type = fwd_txfm_param->tx_type;
  const TX_SIZE tx_size = fwd_txfm_param->tx_size;
  const int rd_transform = fwd_txfm_param->rd_transform;
  const int lossless = fwd_txfm_param->lossless;
  switch (tx_size) {
    case TX_32X32:
      fwd_txfm_32x32(rd_transform, src_diff, coeff, diff_stride, tx_type,
                     fwd_txfm_opt);
      break;
    case TX_16X16:
      fwd_txfm_16x16(src_diff, coeff, diff_stride, tx_type, fwd_txfm_opt);
      break;
    case TX_8X8:
      fwd_txfm_8x8(src_diff, coeff, diff_stride, tx_type, fwd_txfm_opt);
      break;
#if CONFIG_EXT_TX
    case TX_4X8:
      fwd_txfm_4x8(src_diff, coeff, diff_stride, tx_type, fwd_txfm_opt);
      break;
    case TX_8X4:
      fwd_txfm_8x4(src_diff, coeff, diff_stride, tx_type, fwd_txfm_opt);
      break;
#endif  // CONFIG_EXT_TX
    case TX_4X4:
      fwd_txfm_4x4(src_diff, coeff, diff_stride, tx_type, lossless);
      break;
    default:
      assert(0);
      break;
  }
}

#if CONFIG_VP9_HIGHBITDEPTH
void highbd_fwd_txfm(const int16_t *src_diff, tran_low_t *coeff,
                     int diff_stride, FWD_TXFM_PARAM *fwd_txfm_param) {
  const int fwd_txfm_opt = fwd_txfm_param->fwd_txfm_opt;
  const TX_TYPE tx_type = fwd_txfm_param->tx_type;
  const TX_SIZE tx_size = fwd_txfm_param->tx_size;
  const int rd_transform = fwd_txfm_param->rd_transform;
  const int lossless = fwd_txfm_param->lossless;
  const int bd = fwd_txfm_param->bd;
  switch (tx_size) {
    case TX_32X32:
      highbd_fwd_txfm_32x32(rd_transform, src_diff, coeff, diff_stride, tx_type,
                            fwd_txfm_opt, bd);
      break;
    case TX_16X16:
      highbd_fwd_txfm_16x16(src_diff, coeff, diff_stride, tx_type,
                            fwd_txfm_opt, bd);
      break;
    case TX_8X8:
      highbd_fwd_txfm_8x8(src_diff, coeff, diff_stride, tx_type,
                          fwd_txfm_opt, bd);
      break;
#if CONFIG_EXT_TX
    case TX_4X8:
      highbd_fwd_txfm_4x8(src_diff, coeff, diff_stride, tx_type,
                          fwd_txfm_opt, bd);
      break;
    case TX_8X4:
      highbd_fwd_txfm_8x4(src_diff, coeff, diff_stride, tx_type,
                          fwd_txfm_opt, bd);
      break;
#endif  // CONFIG_EXT_TX
    case TX_4X4:
      highbd_fwd_txfm_4x4(src_diff, coeff, diff_stride, tx_type,
                               lossless, bd);
      break;
    default:
      assert(0);
      break;
  }
}
#endif  // CONFIG_VP9_HIGHBITDEPTH
