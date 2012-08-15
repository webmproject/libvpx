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
#include "vp8/common/idct.h"
#include "quantize.h"
#include "vp8/common/reconintra.h"
#include "vp8/common/reconintra4x4.h"
#include "encodemb.h"
#include "vp8/common/invtrans.h"
#include "vp8/common/recon.h"
#include "dct.h"
#include "vp8/common/g_common.h"
#include "encodeintra.h"


#ifdef ENC_DEBUG
extern int enc_debug;
#endif

#if CONFIG_RUNTIME_CPU_DETECT
#define IF_RTCD(x) (x)
#else
#define IF_RTCD(x) NULL
#endif

#if CONFIG_HYBRIDTRANSFORM
extern void vp8_ht_quantize_b(BLOCK *b, BLOCKD *d);
#endif

int vp8_encode_intra(VP8_COMP *cpi, MACROBLOCK *x, int use_16x16_pred) {
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

    vp8_encode_intra16x16mby(IF_RTCD(&cpi->rtcd), x);
  } else {
    for (i = 0; i < 16; i++) {
      x->e_mbd.block[i].bmi.as_mode.first = B_DC_PRED;
      vp8_encode_intra4x4block(IF_RTCD(&cpi->rtcd), x, i);
    }
  }

  intra_pred_var = VARIANCE_INVOKE(&cpi->rtcd.variance, getmbss)(x->src_diff);

  return intra_pred_var;
}

void vp8_encode_intra4x4block(const VP8_ENCODER_RTCD *rtcd,
                              MACROBLOCK *x, int ib) {
  BLOCKD *b = &x->e_mbd.block[ib];
  BLOCK *be = &x->block[ib];

#if CONFIG_HYBRIDTRANSFORM
    int QIndex = x->q_index;
    int active_ht = (QIndex < ACTIVE_HT);
#endif


#if CONFIG_COMP_INTRA_PRED
  if (b->bmi.as_mode.second == (B_PREDICTION_MODE)(B_DC_PRED - 1)) {
#endif
    RECON_INVOKE(&rtcd->common->recon, intra4x4_predict)
    (b, b->bmi.as_mode.first, b->predictor);
#if CONFIG_COMP_INTRA_PRED
  } else {
    RECON_INVOKE(&rtcd->common->recon, comp_intra4x4_predict)
    (b, b->bmi.as_mode.first, b->bmi.as_mode.second, b->predictor);
  }
#endif

  ENCODEMB_INVOKE(&rtcd->encodemb, subb)(be, b, 16);

#if CONFIG_HYBRIDTRANSFORM
    if(active_ht) {
      b->bmi.as_mode.test = b->bmi.as_mode.first;
      txfm_map(b, b->bmi.as_mode.first);
      vp8_fht_c(be->src_diff, be->coeff, 32, b->bmi.as_mode.tx_type, 4);
      vp8_ht_quantize_b(be, b);
      vp8_inverse_htransform_b(IF_RTCD(&rtcd->common->idct), b, 32) ;
    } else {
      x->vp8_short_fdct4x4(be->src_diff, be->coeff, 32) ;
      x->quantize_b(be, b) ;
      vp8_inverse_transform_b(IF_RTCD(&rtcd->common->idct), b, 32) ;
    }
#else
    x->vp8_short_fdct4x4(be->src_diff, be->coeff, 32);
    x->quantize_b(be, b);
    vp8_inverse_transform_b(IF_RTCD(&rtcd->common->idct), b, 32);
#endif

  RECON_INVOKE(&rtcd->common->recon, recon)(b->predictor, b->diff, *(b->base_dst) + b->dst, b->dst_stride);
}

void vp8_encode_intra4x4mby(const VP8_ENCODER_RTCD *rtcd, MACROBLOCK *mb) {
  int i;

#if 0
  MACROBLOCKD *xd = &mb->e_mbd;
  // Intra modes requiring top-right MB reconstructed data have been disabled
  vp8_intra_prediction_down_copy(xd);
#endif

  for (i = 0; i < 16; i++)
    vp8_encode_intra4x4block(rtcd, mb, i);
  return;
}

void vp8_encode_intra16x16mby(const VP8_ENCODER_RTCD *rtcd, MACROBLOCK *x) {
  BLOCK *b = &x->block[0];

  int tx_type = x->e_mbd.mode_info_context->mbmi.txfm_size;

#if CONFIG_COMP_INTRA_PRED
  if (x->e_mbd.mode_info_context->mbmi.second_mode == (MB_PREDICTION_MODE)(DC_PRED - 1))
#endif
    RECON_INVOKE(&rtcd->common->recon, build_intra_predictors_mby)(&x->e_mbd);
#if CONFIG_COMP_INTRA_PRED
  else
    RECON_INVOKE(&rtcd->common->recon, build_comp_intra_predictors_mby)(&x->e_mbd);
#endif

  ENCODEMB_INVOKE(&rtcd->encodemb, submby)(x->src_diff, *(b->base_src), x->e_mbd.predictor, b->src_stride);

#if CONFIG_TX16X16
  if (tx_type == TX_16X16)
    vp8_transform_intra_mby_16x16(x);
  else
#endif
  if (tx_type == TX_8X8)
    vp8_transform_intra_mby_8x8(x);
  else
    vp8_transform_intra_mby(x);

#if CONFIG_TX16X16
  if (tx_type == TX_16X16)
    vp8_quantize_mby_16x16(x);
  else
#endif
  if (tx_type == TX_8X8)
    vp8_quantize_mby_8x8(x);
  else
    vp8_quantize_mby(x);

  if (x->optimize) {
#if CONFIG_TX16X16
    if (tx_type == TX_16X16)
      vp8_optimize_mby_16x16(x, rtcd);
    else
#endif
    if (tx_type == TX_8X8)
      vp8_optimize_mby_8x8(x, rtcd);
    else
      vp8_optimize_mby(x, rtcd);
  }

#if CONFIG_TX16X16
  if (tx_type == TX_16X16)
    vp8_inverse_transform_mby_16x16(IF_RTCD(&rtcd->common->idct), &x->e_mbd);
  else
#endif
  if (tx_type == TX_8X8)
    vp8_inverse_transform_mby_8x8(IF_RTCD(&rtcd->common->idct), &x->e_mbd);
  else
    vp8_inverse_transform_mby(IF_RTCD(&rtcd->common->idct), &x->e_mbd);

#ifdef ENC_DEBUG
  if (enc_debug) {
    int i;
    printf("Intra qcoeff:\n");
    printf("%d %d:\n", x->e_mbd.mb_to_left_edge, x->e_mbd.mb_to_top_edge);
    for (i = 0; i < 400; i++) {
      printf("%3d ", x->e_mbd.qcoeff[i]);
      if (i % 16 == 15) printf("\n");
    }
    printf("Intra dqcoeff:\n");
    for (i = 0; i < 400; i++) {
      printf("%3d ", x->e_mbd.dqcoeff[i]);
      if (i % 16 == 15) printf("\n");
    }
    printf("Intra diff:\n");
    for (i = 0; i < 400; i++) {
      printf("%3d ", x->e_mbd.diff[i]);
      if (i % 16 == 15) printf("\n");
    }
    printf("Intra predictor:\n");
    for (i = 0; i < 400; i++) {
      printf("%3d ", x->e_mbd.predictor[i]);
      if (i % 16 == 15) printf("\n");
    }
    printf("eobs:\n");
    for (i = 0; i < 25; i++)
      printf("%d ", x->e_mbd.block[i].eob);
    printf("\n");
  }
#endif

  RECON_INVOKE(&rtcd->common->recon, recon_mby)
  (IF_RTCD(&rtcd->common->recon), &x->e_mbd);

}

void vp8_encode_intra16x16mbuv(const VP8_ENCODER_RTCD *rtcd, MACROBLOCK *x) {
  int tx_type = x->e_mbd.mode_info_context->mbmi.txfm_size;
#if CONFIG_TX16X16
  if (tx_type == TX_16X16) tx_type = TX_8X8; // 16x16 for U and V should default to 8x8 behavior.
#endif
#if CONFIG_COMP_INTRA_PRED
  if (x->e_mbd.mode_info_context->mbmi.second_uv_mode == (MB_PREDICTION_MODE)(DC_PRED - 1)) {
#endif
    RECON_INVOKE(&rtcd->common->recon, build_intra_predictors_mbuv)(&x->e_mbd);
#if CONFIG_COMP_INTRA_PRED
  } else {
    RECON_INVOKE(&rtcd->common->recon, build_comp_intra_predictors_mbuv)(&x->e_mbd);
  }
#endif

  ENCODEMB_INVOKE(&rtcd->encodemb, submbuv)(x->src_diff, x->src.u_buffer, x->src.v_buffer, x->e_mbd.predictor, x->src.uv_stride);
  if (tx_type == TX_8X8)
    vp8_transform_mbuv_8x8(x);
  else
    vp8_transform_mbuv(x);

  if (tx_type == TX_8X8)
    vp8_quantize_mbuv_8x8(x);
  else
    vp8_quantize_mbuv(x);

#ifdef ENC_DEBUG
  if (enc_debug) {
    int i;
    printf("vp8_encode_intra16x16mbuv\n");
    printf("%d %d:\n", x->e_mbd.mb_to_left_edge, x->e_mbd.mb_to_top_edge);
    printf("qcoeff:\n");
    for (i = 0; i < 400; i++) {
      printf("%3d ", x->e_mbd.qcoeff[i]);
      if (i % 16 == 15) printf("\n");
    }
    printf("dqcoeff:\n");
    for (i = 0; i < 400; i++) {
      printf("%3d ", x->e_mbd.dqcoeff[i]);
      if (i % 16 == 15) printf("\n");
    }
    printf("diff:\n");
    for (i = 0; i < 400; i++) {
      printf("%3d ", x->e_mbd.diff[i]);
      if (i % 16 == 15) printf("\n");
    }
    printf("predictor:\n");
    for (i = 0; i < 400; i++) {
      printf("%3d ", x->e_mbd.predictor[i]);
      if (i % 16 == 15) printf("\n");
    }
    printf("eobs:\n");
    for (i = 0; i < 25; i++)
      printf("%d ", x->e_mbd.block[i].eob);
    printf("\n");
  }
#endif
  if (x->optimize) {
    if (tx_type == TX_8X8)
      vp8_optimize_mbuv_8x8(x, rtcd);
    else
      vp8_optimize_mbuv(x, rtcd);
  }

  if (tx_type == TX_8X8)
    vp8_inverse_transform_mbuv_8x8(IF_RTCD(&rtcd->common->idct), &x->e_mbd);
  else
    vp8_inverse_transform_mbuv(IF_RTCD(&rtcd->common->idct), &x->e_mbd);

  vp8_recon_intra_mbuv(IF_RTCD(&rtcd->common->recon), &x->e_mbd);
}

void vp8_encode_intra8x8(const VP8_ENCODER_RTCD *rtcd,
                         MACROBLOCK *x, int ib) {
  BLOCKD *b = &x->e_mbd.block[ib];
  BLOCK *be = &x->block[ib];
  const int iblock[4] = {0, 1, 4, 5};
  int i;

#if CONFIG_COMP_INTRA_PRED
  if (b->bmi.as_mode.second == (MB_PREDICTION_MODE)(DC_PRED - 1)) {
#endif
    RECON_INVOKE(&rtcd->common->recon, intra8x8_predict)
    (b, b->bmi.as_mode.first, b->predictor);
#if CONFIG_COMP_INTRA_PRED
  } else {
    RECON_INVOKE(&rtcd->common->recon, comp_intra8x8_predict)
    (b, b->bmi.as_mode.first, b->bmi.as_mode.second, b->predictor);
  }
#endif

#if CONFIG_HYBRIDTRANSFORM8X8
  {
    MACROBLOCKD *xd = &x->e_mbd;
    int idx = (ib & 0x02) ? (ib + 2) : ib;

    // generate residual blocks
    vp8_subtract_4b_c(be, b, 16);

    txfm_map(b, pred_mode_conv(b->bmi.as_mode.first));
    vp8_fht_c(be->src_diff, (x->block + idx)->coeff, 32,
              b->bmi.as_mode.tx_type, 8);
    x->quantize_b_8x8(x->block + idx, xd->block + idx);
    vp8_ihtllm_c(xd->block[idx].dqcoeff, xd->block[ib].diff, 32,
                 b->bmi.as_mode.tx_type, 8);

    // reconstruct submacroblock
    for (i = 0; i < 4; i++) {
      b = &xd->block[ib + iblock[i]];
      vp8_recon_b_c(b->predictor, b->diff, *(b->base_dst) + b->dst,
                    b->dst_stride);
    }
  }
#else
  for (i = 0; i < 4; i++) {
    b = &x->e_mbd.block[ib + iblock[i]];
    be = &x->block[ib + iblock[i]];
    ENCODEMB_INVOKE(&rtcd->encodemb, subb)(be, b, 16);
    x->vp8_short_fdct4x4(be->src_diff, be->coeff, 32);
    x->quantize_b(be, b);
    vp8_inverse_transform_b(IF_RTCD(&rtcd->common->idct), b, 32);
    RECON_INVOKE(&rtcd->common->recon, recon)(b->predictor,
                                              b->diff, *(b->base_dst) + b->dst,
                                              b->dst_stride);
  }
#endif
}

extern const int vp8_i8x8_block[4];
void vp8_encode_intra8x8mby(const VP8_ENCODER_RTCD *rtcd, MACROBLOCK *x) {
  int i, ib;

  for (i = 0; i < 4; i++) {
    ib = vp8_i8x8_block[i];
    vp8_encode_intra8x8(rtcd, x, ib);
  }

}

void vp8_encode_intra_uv4x4(const VP8_ENCODER_RTCD *rtcd,
                            MACROBLOCK *x, int ib,
                            int mode, int second) {
  BLOCKD *b = &x->e_mbd.block[ib];
  BLOCK *be = &x->block[ib];

#if CONFIG_COMP_INTRA_PRED
  if (second == -1) {
#endif
    RECON_INVOKE(&rtcd->common->recon, intra_uv4x4_predict)
    (b, mode, b->predictor);
#if CONFIG_COMP_INTRA_PRED
  } else {
    RECON_INVOKE(&rtcd->common->recon, comp_intra_uv4x4_predict)
    (b, mode, second, b->predictor);
  }
#endif

  ENCODEMB_INVOKE(&rtcd->encodemb, subb)(be, b, 8);

  x->vp8_short_fdct4x4(be->src_diff, be->coeff, 16);

  x->quantize_b(be, b);

  vp8_inverse_transform_b(IF_RTCD(&rtcd->common->idct), b, 16);

  RECON_INVOKE(&rtcd->common->recon, recon_uv)(b->predictor,
                                               b->diff, *(b->base_dst) + b->dst, b->dst_stride);
}



void vp8_encode_intra8x8mbuv(const VP8_ENCODER_RTCD *rtcd, MACROBLOCK *x) {
  int i, ib, mode, second;
  BLOCKD *b;
  for (i = 0; i < 4; i++) {
    ib = vp8_i8x8_block[i];
    b = &x->e_mbd.block[ib];
    mode = b->bmi.as_mode.first;
#if CONFIG_COMP_INTRA_PRED
    second = b->bmi.as_mode.second;
#else
    second = -1;
#endif
    /*u */
    vp8_encode_intra_uv4x4(rtcd, x, i + 16, mode, second);
    /*v */
    vp8_encode_intra_uv4x4(rtcd, x, i + 20, mode, second);
  }
}
