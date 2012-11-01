/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include "vpx_ports/config.h"
#include "vpx_rtcd.h"
#include "vp9/common/idct.h"
#include "quantize.h"
#include "vp9/common/reconintra.h"
#include "vp9/common/reconintra4x4.h"
#include "encodemb.h"
#include "vp9/common/invtrans.h"
#include "encodeintra.h"

#if CONFIG_RUNTIME_CPU_DETECT
#define IF_RTCD(x) (x)
#else
#define IF_RTCD(x) NULL
#endif

int vp9_encode_intra(VP9_COMP *cpi, MACROBLOCK *x, int use_16x16_pred) {
  int i;
  int intra_pred_var = 0;
  MB_MODE_INFO * mbmi = &x->e_mbd.mode_info_context->mbmi;
  (void) cpi;

  if (use_16x16_pred) {
    mbmi->mode = DC_PRED;
#if CONFIG_COMP_INTRA_PRED
    mbmi->second_mode = (MB_PREDICTION_MODE)(DC_PRED - 1);
#endif
    mbmi->uv_mode = DC_PRED;
    mbmi->ref_frame = INTRA_FRAME;

    vp9_encode_intra16x16mby(IF_RTCD(&cpi->rtcd), x);
  } else {
    for (i = 0; i < 16; i++) {
      x->e_mbd.block[i].bmi.as_mode.first = B_DC_PRED;
      vp9_encode_intra4x4block(IF_RTCD(&cpi->rtcd), x, i);
    }
  }

  intra_pred_var = vp9_get_mb_ss(x->src_diff);

  return intra_pred_var;
}

void vp9_encode_intra4x4block(const VP9_ENCODER_RTCD *rtcd,
                              MACROBLOCK *x, int ib) {
  BLOCKD *b = &x->e_mbd.block[ib];
  BLOCK *be = &x->block[ib];
  TX_TYPE tx_type;

#if CONFIG_COMP_INTRA_PRED
  if (b->bmi.as_mode.second == (B_PREDICTION_MODE)(B_DC_PRED - 1)) {
#endif
    vp9_intra4x4_predict(b, b->bmi.as_mode.first, b->predictor);
#if CONFIG_COMP_INTRA_PRED
  } else {
    vp9_comp_intra4x4_predict(b, b->bmi.as_mode.first, b->bmi.as_mode.second,
                              b->predictor);
  }
#endif

  vp9_subtract_b(be, b, 16);

  tx_type = get_tx_type(&x->e_mbd, b);
  if (tx_type != DCT_DCT) {
    vp9_fht(be->src_diff, 32, be->coeff, tx_type, 4);
    vp9_ht_quantize_b_4x4(be, b, tx_type);
    vp9_ihtllm_c(b->dqcoeff, b->diff, 32, tx_type, 4);
  } else {
    x->vp9_short_fdct4x4(be->src_diff, be->coeff, 32);
    x->quantize_b_4x4(be, b) ;
    vp9_inverse_transform_b_4x4(IF_RTCD(&rtcd->common->idct), b, 32);
  }

  vp9_recon_b(b->predictor, b->diff, *(b->base_dst) + b->dst, b->dst_stride);
}

void vp9_encode_intra4x4mby(const VP9_ENCODER_RTCD *rtcd, MACROBLOCK *mb) {
  int i;

  for (i = 0; i < 16; i++)
    vp9_encode_intra4x4block(rtcd, mb, i);
  return;
}

void vp9_encode_intra16x16mby(const VP9_ENCODER_RTCD *rtcd, MACROBLOCK *x) {
  MACROBLOCKD *xd = &x->e_mbd;
  BLOCK *b = &x->block[0];
  TX_SIZE tx_size = xd->mode_info_context->mbmi.txfm_size;
  TX_TYPE tx_type;

#if CONFIG_COMP_INTRA_PRED
  if (xd->mode_info_context->mbmi.second_mode == (MB_PREDICTION_MODE)(DC_PRED - 1))
#endif
    vp9_build_intra_predictors_mby(xd);
#if CONFIG_COMP_INTRA_PRED
  else
    vp9_build_comp_intra_predictors_mby(xd);
#endif

  vp9_subtract_mby(x->src_diff, *(b->base_src), xd->predictor, b->src_stride);

  if (tx_size == TX_16X16) {
    BLOCKD  *bd = &xd->block[0];
    tx_type = get_tx_type(xd, bd);
    if (tx_type != DCT_DCT) {
      vp9_fht(b->src_diff, 32, b->coeff, tx_type, 16);
      vp9_quantize_mby_16x16(x);
      if (x->optimize)
        vp9_optimize_mby_16x16(x, rtcd);
      vp9_ihtllm_c(bd->dqcoeff, bd->diff, 32, tx_type, 16);
    } else {
      vp9_transform_mby_16x16(x);
      vp9_quantize_mby_16x16(x);
      if (x->optimize)
        vp9_optimize_mby_16x16(x, rtcd);
      vp9_inverse_transform_mby_16x16(IF_RTCD(&rtcd->common->idct), xd);
    }
  } else if (tx_size == TX_8X8) {
    vp9_transform_mby_8x8(x);
    vp9_quantize_mby_8x8(x);
    if (x->optimize)
      vp9_optimize_mby_8x8(x, rtcd);
    vp9_inverse_transform_mby_8x8(IF_RTCD(&rtcd->common->idct), xd);
  } else {
    vp9_transform_mby_4x4(x);
    vp9_quantize_mby_4x4(x);
    if (x->optimize)
      vp9_optimize_mby_4x4(x, rtcd);
    vp9_inverse_transform_mby_4x4(IF_RTCD(&rtcd->common->idct), xd);
  }

  vp9_recon_mby(xd);
}

void vp9_encode_intra16x16mbuv(const VP9_ENCODER_RTCD *rtcd, MACROBLOCK *x) {
  MACROBLOCKD *xd = &x->e_mbd;
  TX_SIZE tx_size = xd->mode_info_context->mbmi.txfm_size;

#if CONFIG_COMP_INTRA_PRED
  if (xd->mode_info_context->mbmi.second_uv_mode == (MB_PREDICTION_MODE)(DC_PRED - 1)) {
#endif
    vp9_build_intra_predictors_mbuv(xd);
#if CONFIG_COMP_INTRA_PRED
  } else {
    vp9_build_comp_intra_predictors_mbuv(xd);
  }
#endif

  vp9_subtract_mbuv(x->src_diff, x->src.u_buffer, x->src.v_buffer,
                    xd->predictor, x->src.uv_stride);

  if (tx_size == TX_4X4) {
    vp9_transform_mbuv_4x4(x);
    vp9_quantize_mbuv_4x4(x);
    if (x->optimize)
      vp9_optimize_mbuv_4x4(x, rtcd);
    vp9_inverse_transform_mbuv_4x4(IF_RTCD(&rtcd->common->idct), xd);
  } else /* 16x16 or 8x8 */ {
    vp9_transform_mbuv_8x8(x);
    vp9_quantize_mbuv_8x8(x);
    if (x->optimize)
      vp9_optimize_mbuv_8x8(x, rtcd);
    vp9_inverse_transform_mbuv_8x8(IF_RTCD(&rtcd->common->idct), xd);
  }

  vp9_recon_intra_mbuv(xd);
}

void vp9_encode_intra8x8(const VP9_ENCODER_RTCD *rtcd,
                         MACROBLOCK *x, int ib) {
  MACROBLOCKD *xd = &x->e_mbd;
  BLOCKD *b = &xd->block[ib];
  BLOCK *be = &x->block[ib];
  const int iblock[4] = {0, 1, 4, 5};
  int i;
  TX_TYPE tx_type;

#if CONFIG_COMP_INTRA_PRED
  if (b->bmi.as_mode.second == (MB_PREDICTION_MODE)(DC_PRED - 1)) {
#endif
    vp9_intra8x8_predict(b, b->bmi.as_mode.first, b->predictor);
#if CONFIG_COMP_INTRA_PRED
  } else {
    vp9_comp_intra8x8_predict(b, b->bmi.as_mode.first, b->bmi.as_mode.second,
                              b->predictor);
  }
#endif

  if (xd->mode_info_context->mbmi.txfm_size == TX_8X8) {
    int idx = (ib & 0x02) ? (ib + 2) : ib;

    // generate residual blocks
    vp9_subtract_4b_c(be, b, 16);

    tx_type = get_tx_type(xd, xd->block + idx);
    if (tx_type != DCT_DCT) {
      vp9_fht(be->src_diff, 32, (x->block + idx)->coeff,
                tx_type, 8);
      x->quantize_b_8x8(x->block + idx, xd->block + idx);
      vp9_ihtllm_c(xd->block[idx].dqcoeff, xd->block[ib].diff, 32,
                   tx_type, 8);
    } else {
      x->vp9_short_fdct8x8(be->src_diff, (x->block + idx)->coeff, 32);
      x->quantize_b_8x8(x->block + idx, xd->block + idx);
      vp9_idct_idct8(xd->block[idx].dqcoeff, xd->block[ib].diff, 32);
    }
  } else {
    for (i = 0; i < 4; i++) {
      b = &xd->block[ib + iblock[i]];
      be = &x->block[ib + iblock[i]];
      vp9_subtract_b(be, b, 16);
      x->vp9_short_fdct4x4(be->src_diff, be->coeff, 32);
      x->quantize_b_4x4(be, b);
      vp9_inverse_transform_b_4x4(IF_RTCD(&rtcd->common->idct), b, 32);
    }
  }

  // reconstruct submacroblock
  for (i = 0; i < 4; i++) {
    b = &xd->block[ib + iblock[i]];
    vp9_recon_b_c(b->predictor, b->diff, *(b->base_dst) + b->dst,
                  b->dst_stride);
  }
}

void vp9_encode_intra8x8mby(const VP9_ENCODER_RTCD *rtcd, MACROBLOCK *x) {
  int i, ib;

  for (i = 0; i < 4; i++) {
    ib = vp9_i8x8_block[i];
    vp9_encode_intra8x8(rtcd, x, ib);
  }
}

void vp9_encode_intra_uv4x4(const VP9_ENCODER_RTCD *rtcd,
                            MACROBLOCK *x, int ib,
                            int mode, int second) {
  BLOCKD *b = &x->e_mbd.block[ib];
  BLOCK *be = &x->block[ib];

#if CONFIG_COMP_INTRA_PRED
  if (second == -1) {
#endif
    vp9_intra_uv4x4_predict(b, mode, b->predictor);
#if CONFIG_COMP_INTRA_PRED
  } else {
    vp9_comp_intra_uv4x4_predict(b, mode, second, b->predictor);
  }
#endif

  vp9_subtract_b(be, b, 8);

  x->vp9_short_fdct4x4(be->src_diff, be->coeff, 16);
  x->quantize_b_4x4(be, b);
  vp9_inverse_transform_b_4x4(IF_RTCD(&rtcd->common->idct), b, 16);

  vp9_recon_uv_b_c(b->predictor, b->diff, *(b->base_dst) + b->dst,
                   b->dst_stride);
}

void vp9_encode_intra8x8mbuv(const VP9_ENCODER_RTCD *rtcd, MACROBLOCK *x) {
  int i, ib, mode, second;
  BLOCKD *b;

  for (i = 0; i < 4; i++) {
    ib = vp9_i8x8_block[i];
    b = &x->e_mbd.block[ib];
    mode = b->bmi.as_mode.first;
#if CONFIG_COMP_INTRA_PRED
    second = b->bmi.as_mode.second;
#else
    second = -1;
#endif
    /*u */
    vp9_encode_intra_uv4x4(rtcd, x, i + 16, mode, second);
    /*v */
    vp9_encode_intra_uv4x4(rtcd, x, i + 20, mode, second);
  }
}
