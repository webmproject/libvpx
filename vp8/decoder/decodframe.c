/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */


#include "onyxd_int.h"
#include "vp8/common/header.h"
#include "vp8/common/reconintra.h"
#include "vp8/common/reconintra4x4.h"
#include "vp8/common/reconinter.h"
#include "detokenize.h"
#include "vp8/common/invtrans.h"
#include "vp8/common/alloccommon.h"
#include "vp8/common/entropymode.h"
#include "vp8/common/quant_common.h"
#include "vpx_scale/vpxscale.h"
#include "vpx_scale/yv12extend.h"
#include "vp8/common/setupintrarecon.h"

#include "decodemv.h"
#include "vp8/common/extend.h"
#include "vp8/common/modecont.h"
#include "vpx_mem/vpx_mem.h"
#include "vp8/common/idct.h"
#include "dboolhuff.h"

#include "vp8/common/seg_common.h"
#include "vp8/common/entropy.h"
#include "vpx_rtcd.h"

#include <assert.h>
#include <stdio.h>


#define COEFCOUNT_TESTING

static int merge_index(int v, int n, int modulus) {
  int max1 = (n - 1 - modulus / 2) / modulus + 1;
  if (v < max1) v = v * modulus + modulus / 2;
  else {
    int w;
    v -= max1;
    w = v;
    v += (v + modulus - modulus / 2) / modulus;
    while (v % modulus == modulus / 2 ||
           w != v - (v + modulus - modulus / 2) / modulus) v++;
  }
  return v;
}

static int inv_remap_prob(int v, int m) {
  const int n = 256;
  const int modulus = MODULUS_PARAM;
  int i;
  v = merge_index(v, n - 1, modulus);
  if ((m << 1) <= n) {
    i = vp9_inv_recenter_nonneg(v + 1, m);
  } else {
    i = n - 1 - vp9_inv_recenter_nonneg(v + 1, n - 1 - m);
  }
  return i;
}

static vp9_prob read_prob_diff_update(vp9_reader *const bc, int oldp) {
  int delp = vp9_decode_term_subexp(bc, SUBEXP_PARAM, 255);
  return (vp9_prob)inv_remap_prob(delp, oldp);
}

void vp9_init_de_quantizer(VP9D_COMP *pbi) {
  int i;
  int Q;
  VP9_COMMON *const pc = &pbi->common;

  for (Q = 0; Q < QINDEX_RANGE; Q++) {
    pc->Y1dequant[Q][0] = (short)vp9_dc_quant(Q, pc->y1dc_delta_q);
    pc->Y2dequant[Q][0] = (short)vp9_dc2quant(Q, pc->y2dc_delta_q);
    pc->UVdequant[Q][0] = (short)vp9_dc_uv_quant(Q, pc->uvdc_delta_q);

    /* all the ac values =; */
    for (i = 1; i < 16; i++) {
      int rc = vp9_default_zig_zag1d[i];

      pc->Y1dequant[Q][rc] = (short)vp9_ac_yquant(Q);
      pc->Y2dequant[Q][rc] = (short)vp9_ac2quant(Q, pc->y2ac_delta_q);
      pc->UVdequant[Q][rc] = (short)vp9_ac_uv_quant(Q, pc->uvac_delta_q);
    }
  }
}

static void mb_init_dequantizer(VP9D_COMP *pbi, MACROBLOCKD *xd) {
  int i;
  int QIndex;
  VP9_COMMON *const pc = &pbi->common;
  int segment_id = xd->mode_info_context->mbmi.segment_id;

  // Set the Q baseline allowing for any segment level adjustment
  if (vp9_segfeature_active(xd, segment_id, SEG_LVL_ALT_Q)) {
    /* Abs Value */
    if (xd->mb_segment_abs_delta == SEGMENT_ABSDATA)
      QIndex = vp9_get_segdata(xd, segment_id, SEG_LVL_ALT_Q);

    /* Delta Value */
    else {
      QIndex = pc->base_qindex +
               vp9_get_segdata(xd, segment_id, SEG_LVL_ALT_Q);
      QIndex = (QIndex >= 0) ? ((QIndex <= MAXQ) ? QIndex : MAXQ) : 0;    /* Clamp to valid range */
    }
  } else
    QIndex = pc->base_qindex;
  xd->q_index = QIndex;

  /* Set up the block level dequant pointers */
  for (i = 0; i < 16; i++) {
    xd->block[i].dequant = pc->Y1dequant[QIndex];
  }

#if CONFIG_LOSSLESS
  if (!QIndex) {
    pbi->common.rtcd.idct.idct1        = vp9_short_inv_walsh4x4_1_x8_c;
    pbi->common.rtcd.idct.idct16       = vp9_short_inv_walsh4x4_x8_c;
    pbi->common.rtcd.idct.idct1_scalar_add  = vp9_dc_only_inv_walsh_add_c;
    pbi->common.rtcd.idct.iwalsh1      = vp9_short_inv_walsh4x4_1_lossless_c;
    pbi->common.rtcd.idct.iwalsh16     = vp9_short_inv_walsh4x4_lossless_c;
    pbi->idct_add            = vp9_dequant_idct_add_lossless_c;
    pbi->dc_idct_add         = vp9_dequant_dc_idct_add_lossless_c;
    pbi->dc_idct_add_y_block = vp9_dequant_dc_idct_add_y_block_lossless_c;
    pbi->idct_add_y_block    = vp9_dequant_idct_add_y_block_lossless_c;
    pbi->idct_add_uv_block   = vp9_dequant_idct_add_uv_block_lossless_c;
  } else {
    pbi->common.rtcd.idct.idct1        = vp9_short_idct4x4llm_1_c;
    pbi->common.rtcd.idct.idct16       = vp9_short_idct4x4llm_c;
    pbi->common.rtcd.idct.idct1_scalar_add  = vp9_dc_only_idct_add_c;
    pbi->common.rtcd.idct.iwalsh1      = vp9_short_inv_walsh4x4_1_c;
    pbi->common.rtcd.idct.iwalsh16     = vp9_short_inv_walsh4x4_c;
    pbi->idct_add            = vp9_dequant_idct_add;
    pbi->dc_idct_add         = vp9_dequant_dc_idct_add;
    pbi->dc_idct_add_y_block = vp9_dequant_dc_idct_add_y_block;
    pbi->idct_add_y_block    = vp9_dequant_idct_add_y_block;
    pbi->idct_add_uv_block   = vp9_dequant_idct_add_uv_block;
  }
#else
  pbi->idct_add            = vp9_dequant_idct_add;
  pbi->dc_idct_add         = vp9_dequant_dc_idct_add;
  pbi->dc_idct_add_y_block = vp9_dequant_dc_idct_add_y_block;
  pbi->idct_add_y_block    = vp9_dequant_idct_add_y_block;
  pbi->idct_add_uv_block   = vp9_dequant_idct_add_uv_block;
#endif

  for (i = 16; i < 24; i++) {
    xd->block[i].dequant = pc->UVdequant[QIndex];
  }

  xd->block[24].dequant = pc->Y2dequant[QIndex];

}

#if CONFIG_RUNTIME_CPU_DETECT
#define RTCD_VTABLE(x) (&(pbi)->common.rtcd.x)
#else
#define RTCD_VTABLE(x) NULL
#endif

/* skip_recon_mb() is Modified: Instead of writing the result to predictor buffer and then copying it
 *  to dst buffer, we can write the result directly to dst buffer. This eliminates unnecessary copy.
 */
static void skip_recon_mb(VP9D_COMP *pbi, MACROBLOCKD *xd) {
  if (xd->mode_info_context->mbmi.ref_frame == INTRA_FRAME) {
#if CONFIG_SUPERBLOCKS
    if (xd->mode_info_context->mbmi.encoded_as_sb) {
      vp9_build_intra_predictors_sbuv_s(xd);
      vp9_build_intra_predictors_sby_s(xd);
    } else {
#endif
    vp9_build_intra_predictors_mbuv_s(xd);
    vp9_build_intra_predictors_mby_s(xd);
#if CONFIG_SUPERBLOCKS
    }
#endif
  } else {
#if CONFIG_SUPERBLOCKS
    if (xd->mode_info_context->mbmi.encoded_as_sb) {
      vp9_build_inter32x32_predictors_sb(xd, xd->dst.y_buffer,
                                         xd->dst.u_buffer, xd->dst.v_buffer,
                                         xd->dst.y_stride, xd->dst.uv_stride);
    } else {
#endif
    vp9_build_1st_inter16x16_predictors_mb(xd, xd->dst.y_buffer,
                                           xd->dst.u_buffer, xd->dst.v_buffer,
                                           xd->dst.y_stride, xd->dst.uv_stride);

    if (xd->mode_info_context->mbmi.second_ref_frame) {
      vp9_build_2nd_inter16x16_predictors_mb(xd, xd->dst.y_buffer,
                                             xd->dst.u_buffer, xd->dst.v_buffer,
                                             xd->dst.y_stride, xd->dst.uv_stride);
    }
#if CONFIG_SUPERBLOCKS
    }
#endif
  }
}

static void decode_macroblock(VP9D_COMP *pbi, MACROBLOCKD *xd,
                              int mb_row, unsigned int mb_col,
                              BOOL_DECODER* const bc) {
  int eobtotal = 0;
  MB_PREDICTION_MODE mode;
  int i;
  int tx_size;
  TX_TYPE tx_type;
  VP9_COMMON *pc = &pbi->common;
#if CONFIG_SUPERBLOCKS
  int orig_skip_flag = xd->mode_info_context->mbmi.mb_skip_coeff;
#endif

  // re-initialize macroblock dequantizer before detokenization
  if (xd->segmentation_enabled)
    mb_init_dequantizer(pbi, xd);

  tx_size = xd->mode_info_context->mbmi.txfm_size;
  mode = xd->mode_info_context->mbmi.mode;

  if (xd->mode_info_context->mbmi.mb_skip_coeff) {
    vp9_reset_mb_tokens_context(xd);
#if CONFIG_SUPERBLOCKS
    if (xd->mode_info_context->mbmi.encoded_as_sb &&
        (mb_col < pc->mb_cols - 1 || mb_row < pc->mb_rows - 1)) {
      if (mb_col < pc->mb_cols - 1)
        xd->above_context++;
      if (mb_row < pc->mb_rows - 1)
        xd->left_context++;
      vp9_reset_mb_tokens_context(xd);
      if (mb_col < pc->mb_cols - 1)
        xd->above_context--;
      if (mb_row < pc->mb_rows - 1)
        xd->left_context--;
    }
#endif
  } else if (!bool_error(bc)) {
    for (i = 0; i < 25; i++) {
      xd->block[i].eob = 0;
      xd->eobs[i] = 0;
    }
    if (tx_size == TX_16X16) {
      eobtotal = vp9_decode_mb_tokens_16x16(pbi, xd, bc);
    } else if (tx_size == TX_8X8) {
      eobtotal = vp9_decode_mb_tokens_8x8(pbi, xd, bc);
    } else {
      eobtotal = vp9_decode_mb_tokens(pbi, xd, bc);
    }
  }

  //mode = xd->mode_info_context->mbmi.mode;
  if (pbi->common.frame_type != KEY_FRAME)
    vp9_setup_interp_filters(xd, xd->mode_info_context->mbmi.interp_filter,
                             &pbi->common);

  if (eobtotal == 0 && mode != B_PRED && mode != SPLITMV
      && mode != I8X8_PRED
      && !bool_error(bc)) {
    /* Special case:  Force the loopfilter to skip when eobtotal and
     * mb_skip_coeff are zero.
     * */
    xd->mode_info_context->mbmi.mb_skip_coeff = 1;

#if CONFIG_SUPERBLOCKS
    if (!xd->mode_info_context->mbmi.encoded_as_sb || orig_skip_flag)
#endif
    {
      skip_recon_mb(pbi, xd);
      return;
    }
  }

  // moved to be performed before detokenization
//  if (xd->segmentation_enabled)
//    mb_init_dequantizer(pbi, xd);

  /* do prediction */
  if (xd->mode_info_context->mbmi.ref_frame == INTRA_FRAME) {
#if CONFIG_SUPERBLOCKS
    if (xd->mode_info_context->mbmi.encoded_as_sb) {
      vp9_build_intra_predictors_sby_s(xd);
      vp9_build_intra_predictors_sbuv_s(xd);
    } else
#endif
    if (mode != I8X8_PRED) {
      vp9_build_intra_predictors_mbuv(xd);
      if (mode != B_PRED) {
        vp9_build_intra_predictors_mby(xd);
      }
    }
  } else {
#if CONFIG_SUPERBLOCKS
    if (xd->mode_info_context->mbmi.encoded_as_sb) {
      vp9_build_inter32x32_predictors_sb(xd, xd->dst.y_buffer,
                                         xd->dst.u_buffer, xd->dst.v_buffer,
                                         xd->dst.y_stride, xd->dst.uv_stride);
    } else
#endif
    vp9_build_inter_predictors_mb(xd);
  }

  /* dequantization and idct */
  if (mode == I8X8_PRED) {
    for (i = 0; i < 4; i++) {
      int ib = vp9_i8x8_block[i];
      const int iblock[4] = {0, 1, 4, 5};
      int j;
      int i8x8mode;
      BLOCKD *b;

      int idx = (ib & 0x02) ? (ib + 2) : ib;

      short *q  = xd->block[idx].qcoeff;
      short *dq = xd->block[0].dequant;
      unsigned char *pre = xd->block[ib].predictor;
      unsigned char *dst = *(xd->block[ib].base_dst) + xd->block[ib].dst;
      int stride = xd->dst.y_stride;

      b = &xd->block[ib];
      i8x8mode = b->bmi.as_mode.first;
      vp9_intra8x8_predict(b, i8x8mode, b->predictor);

      if (xd->mode_info_context->mbmi.txfm_size == TX_8X8) {
        tx_type = get_tx_type(xd, &xd->block[idx]);
        if (tx_type != DCT_DCT) {
          vp9_ht_dequant_idct_add_8x8_c(tx_type,
                                        q, dq, pre, dst, 16, stride);
        } else {
          vp9_dequant_idct_add_8x8_c(q, dq, pre, dst, 16, stride);
        }
        q += 64;
      } else {
        for (j = 0; j < 4; j++) {
          b = &xd->block[ib + iblock[j]];
          vp9_dequant_idct_add(b->qcoeff, b->dequant, b->predictor,
                                 *(b->base_dst) + b->dst, 16, b->dst_stride);
        }
      }
      b = &xd->block[16 + i];
      vp9_intra_uv4x4_predict(b, i8x8mode, b->predictor);
      pbi->idct_add(b->qcoeff, b->dequant, b->predictor,
                    *(b->base_dst) + b->dst, 8, b->dst_stride);
      b = &xd->block[20 + i];
      vp9_intra_uv4x4_predict(b, i8x8mode, b->predictor);
      pbi->idct_add(b->qcoeff, b->dequant, b->predictor,
                    *(b->base_dst) + b->dst, 8, b->dst_stride);
    }
  } else if (mode == B_PRED) {
    for (i = 0; i < 16; i++) {
      BLOCKD *b = &xd->block[i];
      int b_mode = xd->mode_info_context->bmi[i].as_mode.first;
#if CONFIG_COMP_INTRA_PRED
      int b_mode2 = xd->mode_info_context->bmi[i].as_mode.second;

      if (b_mode2 == (B_PREDICTION_MODE)(B_DC_PRED - 1)) {
#endif
        vp9_intra4x4_predict(b, b_mode, b->predictor);
#if CONFIG_COMP_INTRA_PRED
      } else {
        vp9_comp_intra4x4_predict(b, b_mode, b_mode2, b->predictor);
      }
#endif

      tx_type = get_tx_type(xd, b);
      if (tx_type != DCT_DCT) {
        vp9_ht_dequant_idct_add_c(tx_type, b->qcoeff,
                                  b->dequant, b->predictor,
                                  *(b->base_dst) + b->dst, 16, b->dst_stride);
      } else {
        vp9_dequant_idct_add(b->qcoeff, b->dequant, b->predictor,
                               *(b->base_dst) + b->dst, 16, b->dst_stride);
      }
    }
  } else if (mode == SPLITMV) {
    if (tx_size == TX_8X8) {
      vp9_dequant_idct_add_y_block_8x8(xd->qcoeff, xd->block[0].dequant,
                                         xd->predictor, xd->dst.y_buffer,
                                         xd->dst.y_stride, xd->eobs, xd);
    } else {
      pbi->idct_add_y_block(xd->qcoeff, xd->block[0].dequant,
                                       xd->predictor, xd->dst.y_buffer,
                                       xd->dst.y_stride, xd->eobs);
    }
  } else {
    BLOCKD *b = &xd->block[24];

    if (tx_size == TX_16X16) {
      BLOCKD *bd = &xd->block[0];
      tx_type = get_tx_type(xd, bd);
      if (tx_type != DCT_DCT) {
        vp9_ht_dequant_idct_add_16x16_c(tx_type, xd->qcoeff,
                                        xd->block[0].dequant, xd->predictor,
                                        xd->dst.y_buffer, 16, xd->dst.y_stride);
      } else {
        vp9_dequant_idct_add_16x16(xd->qcoeff, xd->block[0].dequant,
                                     xd->predictor, xd->dst.y_buffer,
                                     16, xd->dst.y_stride);
      }
    } else if (tx_size == TX_8X8) {
#if CONFIG_SUPERBLOCKS
      void *orig = xd->mode_info_context;
      int n, num = xd->mode_info_context->mbmi.encoded_as_sb ? 4 : 1;
      for (n = 0; n < num; n++) {
        int x_idx = n & 1, y_idx = n >> 1;
        if (num == 4 && (mb_col + x_idx >= pc->mb_cols ||
                         mb_row + y_idx >= pc->mb_rows))
          continue;

        if (n != 0) {
          for (i = 0; i < 25; i++) {
            xd->block[i].eob = 0;
            xd->eobs[i] = 0;
          }
          xd->above_context = pc->above_context + mb_col + (n & 1);
          xd->left_context = pc->left_context + (n >> 1);
          xd->mode_info_context = orig;
          xd->mode_info_context += (n & 1);
          xd->mode_info_context += (n >> 1) * pc->mode_info_stride;
          if (!orig_skip_flag) {
            eobtotal = vp9_decode_mb_tokens_8x8(pbi, xd, bc);
            if (eobtotal == 0) // skip loopfilter
              xd->mode_info_context->mbmi.mb_skip_coeff = 1;
          } else {
            vp9_reset_mb_tokens_context(xd);
          }
        }

        if (xd->mode_info_context->mbmi.mb_skip_coeff)
          continue; // only happens for SBs, which are already in dest buffer
#endif
      vp9_dequantize_b_2x2(b);
      IDCT_INVOKE(RTCD_VTABLE(idct), ihaar2)(&b->dqcoeff[0], b->diff, 8);
      ((int *)b->qcoeff)[0] = 0;// 2nd order block are set to 0 after inverse transform
      ((int *)b->qcoeff)[1] = 0;
      ((int *)b->qcoeff)[2] = 0;
      ((int *)b->qcoeff)[3] = 0;
      ((int *)b->qcoeff)[4] = 0;
      ((int *)b->qcoeff)[5] = 0;
      ((int *)b->qcoeff)[6] = 0;
      ((int *)b->qcoeff)[7] = 0;
#if CONFIG_SUPERBLOCKS
      if (xd->mode_info_context->mbmi.encoded_as_sb) {
        vp9_dequant_dc_idct_add_y_block_8x8_inplace_c(xd->qcoeff,
          xd->block[0].dequant,
          xd->dst.y_buffer + (n >> 1) * 16 * xd->dst.y_stride + (n & 1) * 16,
          xd->dst.y_stride, xd->eobs, xd->block[24].diff, xd);
        // do UV inline also
        vp9_dequant_idct_add_uv_block_8x8_inplace_c(xd->qcoeff + 16 * 16,
          xd->block[16].dequant,
          xd->dst.u_buffer + (n >> 1) * 8 * xd->dst.uv_stride + (n & 1) * 8,
          xd->dst.v_buffer + (n >> 1) * 8 * xd->dst.uv_stride + (n & 1) * 8,
          xd->dst.uv_stride, xd->eobs + 16, xd);
      } else
#endif
        vp9_dequant_dc_idct_add_y_block_8x8(xd->qcoeff,
          xd->block[0].dequant, xd->predictor, xd->dst.y_buffer,
          xd->dst.y_stride, xd->eobs, xd->block[24].diff, xd);
#if CONFIG_SUPERBLOCKS
      }
      xd->mode_info_context = orig;
#endif
    } else {
      vp9_dequantize_b(b);
      if (xd->eobs[24] > 1) {
        IDCT_INVOKE(RTCD_VTABLE(idct), iwalsh16)(&b->dqcoeff[0], b->diff);
        ((int *)b->qcoeff)[0] = 0;
        ((int *)b->qcoeff)[1] = 0;
        ((int *)b->qcoeff)[2] = 0;
        ((int *)b->qcoeff)[3] = 0;
        ((int *)b->qcoeff)[4] = 0;
        ((int *)b->qcoeff)[5] = 0;
        ((int *)b->qcoeff)[6] = 0;
        ((int *)b->qcoeff)[7] = 0;
      } else {
        IDCT_INVOKE(RTCD_VTABLE(idct), iwalsh1)(&b->dqcoeff[0], b->diff);
        ((int *)b->qcoeff)[0] = 0;
      }

      pbi->dc_idct_add_y_block(xd->qcoeff, xd->block[0].dequant, xd->predictor,
                               xd->dst.y_buffer, xd->dst.y_stride, xd->eobs,
                               xd->block[24].diff);
    }
  }

#if CONFIG_SUPERBLOCKS
  if (!xd->mode_info_context->mbmi.encoded_as_sb) {
#endif
    if ((tx_size == TX_8X8 &&
         xd->mode_info_context->mbmi.mode != I8X8_PRED &&
         xd->mode_info_context->mbmi.mode != SPLITMV)
        || tx_size == TX_16X16
       )
      vp9_dequant_idct_add_uv_block_8x8
          (xd->qcoeff + 16 * 16, xd->block[16].dequant,
           xd->predictor + 16 * 16, xd->dst.u_buffer, xd->dst.v_buffer,
           xd->dst.uv_stride, xd->eobs + 16, xd); //
    else if (xd->mode_info_context->mbmi.mode != I8X8_PRED)
      pbi->idct_add_uv_block(xd->qcoeff + 16 * 16, xd->block[16].dequant,
           xd->predictor + 16 * 16, xd->dst.u_buffer, xd->dst.v_buffer,
           xd->dst.uv_stride, xd->eobs + 16);
#if CONFIG_SUPERBLOCKS
  }
#endif
}


static int get_delta_q(vp9_reader *bc, int prev, int *q_update) {
  int ret_val = 0;

  if (vp9_read_bit(bc)) {
    ret_val = vp9_read_literal(bc, 4);

    if (vp9_read_bit(bc))
      ret_val = -ret_val;
  }

  /* Trigger a quantizer update if the delta-q value has changed */
  if (ret_val != prev)
    *q_update = 1;

  return ret_val;
}

#ifdef PACKET_TESTING
#include <stdio.h>
FILE *vpxlog = 0;
#endif

/* Decode a row of Superblocks (2x2 region of MBs) */
static void
decode_sb_row(VP9D_COMP *pbi, VP9_COMMON *pc, int mbrow, MACROBLOCKD *xd,
              BOOL_DECODER* const bc) {
  int i;
  int sb_col;
  int mb_row, mb_col;
  int recon_yoffset, recon_uvoffset;
  int ref_fb_idx = pc->lst_fb_idx;
  int dst_fb_idx = pc->new_fb_idx;
  int recon_y_stride = pc->yv12_fb[ref_fb_idx].y_stride;
  int recon_uv_stride = pc->yv12_fb[ref_fb_idx].uv_stride;
  int row_delta[4] = { 0, +1,  0, -1};
  int col_delta[4] = { +1, -1, +1, +1};
  int sb_cols = (pc->mb_cols + 1) >> 1;

  // For a SB there are 2 left contexts, each pertaining to a MB row within
  vpx_memset(pc->left_context, 0, sizeof(pc->left_context));

  mb_row = mbrow;
  mb_col = 0;

  for (sb_col = 0; sb_col < sb_cols; sb_col++) {
    MODE_INFO *mi = xd->mode_info_context;

#if CONFIG_SUPERBLOCKS
    mi->mbmi.encoded_as_sb = vp9_read(bc, pc->sb_coded);
#endif

    // Process the 4 MBs within the SB in the order:
    // top-left, top-right, bottom-left, bottom-right
    for (i = 0; i < 4; i++) {
      int dy = row_delta[i];
      int dx = col_delta[i];
      int offset_extended = dy * xd->mode_info_stride + dx;

      xd->mb_index = i;

      mi = xd->mode_info_context;
      if ((mb_row >= pc->mb_rows) || (mb_col >= pc->mb_cols)) {
        // MB lies outside frame, skip on to next
        mb_row += dy;
        mb_col += dx;
        xd->mode_info_context += offset_extended;
        xd->prev_mode_info_context += offset_extended;
        continue;
      }

      // Set above context pointer
      xd->above_context = pc->above_context + mb_col;
      xd->left_context = pc->left_context + (i >> 1);

      /* Distance of Mb to the various image edges.
       * These are specified to 8th pel as they are always compared to
       * values that are in 1/8th pel units
       */
      xd->mb_to_top_edge = -((mb_row * 16)) << 3;
      xd->mb_to_bottom_edge = ((pc->mb_rows - 1 - mb_row) * 16) << 3;

      xd->mb_to_left_edge = -((mb_col * 16) << 3);
      xd->mb_to_right_edge = ((pc->mb_cols - 1 - mb_col) * 16) << 3;

      xd->up_available = (mb_row != 0);
      xd->left_available = (mb_col != 0);


      recon_yoffset = (mb_row * recon_y_stride * 16) + (mb_col * 16);
      recon_uvoffset = (mb_row * recon_uv_stride * 8) + (mb_col * 8);

      xd->dst.y_buffer = pc->yv12_fb[dst_fb_idx].y_buffer + recon_yoffset;
      xd->dst.u_buffer = pc->yv12_fb[dst_fb_idx].u_buffer + recon_uvoffset;
      xd->dst.v_buffer = pc->yv12_fb[dst_fb_idx].v_buffer + recon_uvoffset;

#if CONFIG_SUPERBLOCKS
      if (i)
        mi->mbmi.encoded_as_sb = 0;
#endif
      vp9_decode_mb_mode_mv(pbi, xd, mb_row, mb_col, bc);

      update_blockd_bmi(xd);

      /* Select the appropriate reference frame for this MB */
      if (xd->mode_info_context->mbmi.ref_frame == LAST_FRAME)
        ref_fb_idx = pc->lst_fb_idx;
      else if (xd->mode_info_context->mbmi.ref_frame == GOLDEN_FRAME)
        ref_fb_idx = pc->gld_fb_idx;
      else
        ref_fb_idx = pc->alt_fb_idx;

      xd->pre.y_buffer = pc->yv12_fb[ref_fb_idx].y_buffer + recon_yoffset;
      xd->pre.u_buffer = pc->yv12_fb[ref_fb_idx].u_buffer + recon_uvoffset;
      xd->pre.v_buffer = pc->yv12_fb[ref_fb_idx].v_buffer + recon_uvoffset;

      if (xd->mode_info_context->mbmi.second_ref_frame) {
        int second_ref_fb_idx;

        /* Select the appropriate reference frame for this MB */
        if (xd->mode_info_context->mbmi.second_ref_frame == LAST_FRAME)
          second_ref_fb_idx = pc->lst_fb_idx;
        else if (xd->mode_info_context->mbmi.second_ref_frame ==
                 GOLDEN_FRAME)
          second_ref_fb_idx = pc->gld_fb_idx;
        else
          second_ref_fb_idx = pc->alt_fb_idx;

        xd->second_pre.y_buffer =
          pc->yv12_fb[second_ref_fb_idx].y_buffer + recon_yoffset;
        xd->second_pre.u_buffer =
          pc->yv12_fb[second_ref_fb_idx].u_buffer + recon_uvoffset;
        xd->second_pre.v_buffer =
          pc->yv12_fb[second_ref_fb_idx].v_buffer + recon_uvoffset;
      }

      if (xd->mode_info_context->mbmi.ref_frame != INTRA_FRAME) {
        /* propagate errors from reference frames */
        xd->corrupted |= pc->yv12_fb[ref_fb_idx].corrupted;
      }

#if CONFIG_SUPERBLOCKS
      if (xd->mode_info_context->mbmi.encoded_as_sb) {
        if (mb_col < pc->mb_cols - 1)
          mi[1] = mi[0];
        if (mb_row < pc->mb_rows - 1) {
          mi[pc->mode_info_stride] = mi[0];
          if (mb_col < pc->mb_cols - 1)
            mi[pc->mode_info_stride + 1] = mi[0];
        }
      }
#endif
      vp9_intra_prediction_down_copy(xd);
      decode_macroblock(pbi, xd, mb_row, mb_col, bc);

      /* check if the boolean decoder has suffered an error */
      xd->corrupted |= bool_error(bc);

#if CONFIG_SUPERBLOCKS
      if (mi->mbmi.encoded_as_sb) {
        assert(!i);
        mb_col += 2;
        xd->mode_info_context += 2;
        xd->prev_mode_info_context += 2;
        break;
      }
#endif

      // skip to next MB
      xd->mode_info_context += offset_extended;
      xd->prev_mode_info_context += offset_extended;
      mb_row += dy;
      mb_col += dx;
    }
  }

  /* skip prediction column */
  xd->mode_info_context += 1 - (pc->mb_cols & 0x1) + xd->mode_info_stride;
  xd->prev_mode_info_context += 1 - (pc->mb_cols & 0x1) + xd->mode_info_stride;
}

static unsigned int read_partition_size(const unsigned char *cx_size) {
  const unsigned int size =
    cx_size[0] + (cx_size[1] << 8) + (cx_size[2] << 16);
  return size;
}

static int read_is_valid(const unsigned char *start,
                         size_t               len,
                         const unsigned char *end) {
  return (start + len > start && start + len <= end);
}


static void setup_token_decoder(VP9D_COMP *pbi,
                                const unsigned char *cx_data,
                                BOOL_DECODER* const bool_decoder) {
  VP9_COMMON          *pc = &pbi->common;
  const unsigned char *user_data_end = pbi->Source + pbi->source_sz;
  const unsigned char *partition;

  ptrdiff_t            partition_size;
  ptrdiff_t            bytes_left;

  // Set up pointers to token partition
  partition = cx_data;
  bytes_left = user_data_end - partition;
  partition_size = bytes_left;

  /* Validate the calculated partition length. If the buffer
   * described by the partition can't be fully read, then restrict
   * it to the portion that can be (for EC mode) or throw an error.
   */
  if (!read_is_valid(partition, partition_size, user_data_end)) {
    vpx_internal_error(&pc->error, VPX_CODEC_CORRUPT_FRAME,
                       "Truncated packet or corrupt partition "
                       "%d length", 1);
  }

  if (vp9_start_decode(bool_decoder, partition, partition_size))
    vpx_internal_error(&pc->error, VPX_CODEC_MEM_ERROR,
                       "Failed to allocate bool decoder %d", 1);
}

static void init_frame(VP9D_COMP *pbi) {
  VP9_COMMON *const pc = &pbi->common;
  MACROBLOCKD *const xd  = &pbi->mb;

  if (pc->frame_type == KEY_FRAME) {
    /* Various keyframe initializations */
    vp9_init_mv_probs(pc);

    vp9_init_mbmode_probs(pc);
    vp9_default_bmode_probs(pc->fc.bmode_prob);

    vp9_default_coef_probs(pc);
    vp9_kf_default_bmode_probs(pc->kf_bmode_prob);

    // Reset the segment feature data to the default stats:
    // Features disabled, 0, with delta coding (Default state).
    vp9_clearall_segfeatures(xd);

    xd->mb_segment_abs_delta = SEGMENT_DELTADATA;

    /* reset the mode ref deltasa for loop filter */
    vpx_memset(xd->ref_lf_deltas, 0, sizeof(xd->ref_lf_deltas));
    vpx_memset(xd->mode_lf_deltas, 0, sizeof(xd->mode_lf_deltas));

    /* All buffers are implicitly updated on key frames. */
    pc->refresh_golden_frame = 1;
    pc->refresh_alt_ref_frame = 1;
    pc->copy_buffer_to_gf = 0;
    pc->copy_buffer_to_arf = 0;

    /* Note that Golden and Altref modes cannot be used on a key frame so
     * ref_frame_sign_bias[] is undefined and meaningless
     */
    pc->ref_frame_sign_bias[GOLDEN_FRAME] = 0;
    pc->ref_frame_sign_bias[ALTREF_FRAME] = 0;

    vp9_init_mode_contexts(&pbi->common);
    vpx_memcpy(&pc->lfc, &pc->fc, sizeof(pc->fc));
    vpx_memcpy(&pc->lfc_a, &pc->fc, sizeof(pc->fc));

    vpx_memcpy(pbi->common.fc.vp8_mode_contexts,
               pbi->common.fc.mode_context,
               sizeof(pbi->common.fc.mode_context));
    vpx_memset(pc->prev_mip, 0,
               (pc->mb_cols + 1) * (pc->mb_rows + 1)* sizeof(MODE_INFO));
    vpx_memset(pc->mip, 0,
               (pc->mb_cols + 1) * (pc->mb_rows + 1)* sizeof(MODE_INFO));

    vp9_update_mode_info_border(pc, pc->mip);
    vp9_update_mode_info_in_image(pc, pc->mi);

  } else {

    if (!pc->use_bilinear_mc_filter)
      pc->mcomp_filter_type = EIGHTTAP;
    else
      pc->mcomp_filter_type = BILINEAR;

    /* To enable choice of different interpolation filters */
    vp9_setup_interp_filters(xd, pc->mcomp_filter_type, pc);
  }

  xd->mode_info_context = pc->mi;
  xd->prev_mode_info_context = pc->prev_mi;
  xd->frame_type = pc->frame_type;
  xd->mode_info_context->mbmi.mode = DC_PRED;
  xd->mode_info_stride = pc->mode_info_stride;
  xd->corrupted = 0; /* init without corruption */

  xd->fullpixel_mask = 0xffffffff;
  if (pc->full_pixel)
    xd->fullpixel_mask = 0xfffffff8;

}

#if 0
static void read_coef_probs2(VP9D_COMP *pbi) {
  const vp9_prob grpupd = 192;
  int i, j, k, l;
  vp9_reader *const bc = &pbi->bc;
  VP9_COMMON *const pc = &pbi->common;
  for (l = 0; l < ENTROPY_NODES; l++) {
    if (vp9_read(bc, grpupd)) {
      // printf("Decoding %d\n", l);
      for (i = 0; i < BLOCK_TYPES; i++)
        for (j = !i; j < COEF_BANDS; j++)
          for (k = 0; k < PREV_COEF_CONTEXTS; k++) {
            if (k >= 3 && ((i == 0 && j == 1) ||
                           (i > 0 && j == 0)))
              continue;
            {
              vp9_prob *const p = pc->fc.coef_probs [i][j][k] + l;
              int u = vp9_read(bc, COEF_UPDATE_PROB);
              if (u) *p = read_prob_diff_update(bc, *p);
            }
          }
    }
  }
  if (pbi->common.txfm_mode == ALLOW_8X8) {
    for (l = 0; l < ENTROPY_NODES; l++) {
      if (vp9_read(bc, grpupd)) {
        for (i = 0; i < BLOCK_TYPES_8X8; i++)
          for (j = !i; j < COEF_BANDS; j++)
            for (k = 0; k < PREV_COEF_CONTEXTS; k++) {
              if (k >= 3 && ((i == 0 && j == 1) ||
                             (i > 0 && j == 0)))
                continue;
              {
                vp9_prob *const p = pc->fc.coef_probs_8x8 [i][j][k] + l;

                int u = vp9_read(bc, COEF_UPDATE_PROB_8X8);
                if (u) *p = read_prob_diff_update(bc, *p);
              }
            }
      }
    }
  }
}
#endif

static void read_coef_probs_common(
    BOOL_DECODER* const bc,
    vp9_prob coef_probs[BLOCK_TYPES][COEF_BANDS]
                       [PREV_COEF_CONTEXTS][ENTROPY_NODES]) {
  int i, j, k, l;

  if (vp9_read_bit(bc)) {
    for (i = 0; i < BLOCK_TYPES; i++) {
      for (j = !i; j < COEF_BANDS; j++) {
        /* NB: This j loop starts from 1 on block type i == 0 */
        for (k = 0; k < PREV_COEF_CONTEXTS; k++) {
          if (k >= 3 && ((i == 0 && j == 1) ||
                         (i > 0 && j == 0)))
            continue;
          for (l = 0; l < ENTROPY_NODES; l++) {
            vp9_prob *const p = coef_probs[i][j][k] + l;

            if (vp9_read(bc, COEF_UPDATE_PROB)) {
              *p = read_prob_diff_update(bc, *p);
            }
          }
        }
      }
    }
  }
}

static void read_coef_probs(VP9D_COMP *pbi, BOOL_DECODER* const bc) {
  VP9_COMMON *const pc = &pbi->common;

  read_coef_probs_common(bc, pc->fc.coef_probs);
  read_coef_probs_common(bc, pc->fc.hybrid_coef_probs);

  if (pbi->common.txfm_mode != ONLY_4X4) {
    read_coef_probs_common(bc, pc->fc.coef_probs_8x8);
    read_coef_probs_common(bc, pc->fc.hybrid_coef_probs_8x8);
  }
  if (pbi->common.txfm_mode > ALLOW_8X8) {
    read_coef_probs_common(bc, pc->fc.coef_probs_16x16);
    read_coef_probs_common(bc, pc->fc.hybrid_coef_probs_16x16);
  }
}

int vp9_decode_frame(VP9D_COMP *pbi) {
  BOOL_DECODER header_bc, residual_bc;
  VP9_COMMON *const pc = &pbi->common;
  MACROBLOCKD *const xd  = &pbi->mb;
  const unsigned char *data = (const unsigned char *)pbi->Source;
  const unsigned char *data_end = data + pbi->source_sz;
  ptrdiff_t first_partition_length_in_bytes = 0;

  int mb_row;
  int i, j;
  int corrupt_tokens = 0;

  /* start with no corruption of current frame */
  xd->corrupted = 0;
  pc->yv12_fb[pc->new_fb_idx].corrupted = 0;

  if (data_end - data < 3) {
    vpx_internal_error(&pc->error, VPX_CODEC_CORRUPT_FRAME,
                       "Truncated packet");
  } else {
    pc->last_frame_type = pc->frame_type;
    pc->frame_type = (FRAME_TYPE)(data[0] & 1);
    pc->version = (data[0] >> 1) & 7;
    pc->show_frame = (data[0] >> 4) & 1;
    first_partition_length_in_bytes =
      (data[0] | (data[1] << 8) | (data[2] << 16)) >> 5;

    if ((data + first_partition_length_in_bytes > data_end
         || data + first_partition_length_in_bytes < data))
      vpx_internal_error(&pc->error, VPX_CODEC_CORRUPT_FRAME,
                         "Truncated packet or corrupt partition 0 length");

    data += 3;

    vp9_setup_version(pc);

    if (pc->frame_type == KEY_FRAME) {
      const int Width = pc->Width;
      const int Height = pc->Height;

      /* vet via sync code */
      /* When error concealment is enabled we should only check the sync
       * code if we have enough bits available
       */
      if (data + 3 < data_end) {
        if (data[0] != 0x9d || data[1] != 0x01 || data[2] != 0x2a)
          vpx_internal_error(&pc->error, VPX_CODEC_UNSUP_BITSTREAM,
                             "Invalid frame sync code");
      }

      /* If error concealment is enabled we should only parse the new size
       * if we have enough data. Otherwise we will end up with the wrong
       * size.
       */
      if (data + 6 < data_end) {
        pc->Width = (data[3] | (data[4] << 8)) & 0x3fff;
        pc->horiz_scale = data[4] >> 6;
        pc->Height = (data[5] | (data[6] << 8)) & 0x3fff;
        pc->vert_scale = data[6] >> 6;
      }
      data += 7;

      if (Width != pc->Width  ||  Height != pc->Height) {
        if (pc->Width <= 0) {
          pc->Width = Width;
          vpx_internal_error(&pc->error, VPX_CODEC_CORRUPT_FRAME,
                             "Invalid frame width");
        }

        if (pc->Height <= 0) {
          pc->Height = Height;
          vpx_internal_error(&pc->error, VPX_CODEC_CORRUPT_FRAME,
                             "Invalid frame height");
        }

        if (vp9_alloc_frame_buffers(pc, pc->Width, pc->Height))
          vpx_internal_error(&pc->error, VPX_CODEC_MEM_ERROR,
                             "Failed to allocate frame buffers");
      }
    }
  }

  if ((!pbi->decoded_key_frame && pc->frame_type != KEY_FRAME) ||
      pc->Width == 0 || pc->Height == 0) {
    return -1;
  }

  init_frame(pbi);

  if (vp9_start_decode(&header_bc, data, first_partition_length_in_bytes))
    vpx_internal_error(&pc->error, VPX_CODEC_MEM_ERROR,
                       "Failed to allocate bool decoder 0");
  if (pc->frame_type == KEY_FRAME) {
    pc->clr_type    = (YUV_TYPE)vp9_read_bit(&header_bc);
    pc->clamp_type  = (CLAMP_TYPE)vp9_read_bit(&header_bc);
  }

  /* Is segmentation enabled */
  xd->segmentation_enabled = (unsigned char)vp9_read_bit(&header_bc);

  if (xd->segmentation_enabled) {
    // Read whether or not the segmentation map is being explicitly
    // updated this frame.
    xd->update_mb_segmentation_map = (unsigned char)vp9_read_bit(&header_bc);

    // If so what method will be used.
    if (xd->update_mb_segmentation_map) {
      // Which macro block level features are enabled

      // Read the probs used to decode the segment id for each macro
      // block.
      for (i = 0; i < MB_FEATURE_TREE_PROBS; i++) {
          xd->mb_segment_tree_probs[i] = vp9_read_bit(&header_bc) ?
              (vp9_prob)vp9_read_literal(&header_bc, 8) : 255;
      }

      // Read the prediction probs needed to decode the segment id
      pc->temporal_update = (unsigned char)vp9_read_bit(&header_bc);
      for (i = 0; i < PREDICTION_PROBS; i++) {
        if (pc->temporal_update) {
          pc->segment_pred_probs[i] = vp9_read_bit(&header_bc) ?
              (vp9_prob)vp9_read_literal(&header_bc, 8) : 255;
        } else {
          pc->segment_pred_probs[i] = 255;
        }
      }
    }
    // Is the segment data being updated
    xd->update_mb_segmentation_data = (unsigned char)vp9_read_bit(&header_bc);

    if (xd->update_mb_segmentation_data) {
      int data;

      xd->mb_segment_abs_delta = (unsigned char)vp9_read_bit(&header_bc);

      vp9_clearall_segfeatures(xd);

      // For each segmentation...
      for (i = 0; i < MAX_MB_SEGMENTS; i++) {
        // For each of the segments features...
        for (j = 0; j < SEG_LVL_MAX; j++) {
          // Is the feature enabled
          if (vp9_read_bit(&header_bc)) {
            // Update the feature data and mask
            vp9_enable_segfeature(xd, i, j);

            data = (signed char)vp9_read_literal(
                     &header_bc, vp9_seg_feature_data_bits(j));

            // Is the segment data signed..
            if (vp9_is_segfeature_signed(j)) {
              if (vp9_read_bit(&header_bc))
                data = - data;
            }
          } else
            data = 0;

          vp9_set_segdata(xd, i, j, data);
        }
      }
    }
  }

  // Read common prediction model status flag probability updates for the
  // reference frame
  if (pc->frame_type == KEY_FRAME) {
    // Set the prediction probabilities to defaults
    pc->ref_pred_probs[0] = 120;
    pc->ref_pred_probs[1] = 80;
    pc->ref_pred_probs[2] = 40;
  } else {
    for (i = 0; i < PREDICTION_PROBS; i++) {
      if (vp9_read_bit(&header_bc))
        pc->ref_pred_probs[i] = (vp9_prob)vp9_read_literal(&header_bc, 8);
    }
  }

#if CONFIG_SUPERBLOCKS
  pc->sb_coded = vp9_read_literal(&header_bc, 8);
#endif

  /* Read the loop filter level and type */
  pc->txfm_mode = vp9_read_literal(&header_bc, 2);
  if (pc->txfm_mode == TX_MODE_SELECT) {
    pc->prob_tx[0] = vp9_read_literal(&header_bc, 8);
    pc->prob_tx[1] = vp9_read_literal(&header_bc, 8);
  }

  pc->filter_type = (LOOPFILTERTYPE) vp9_read_bit(&header_bc);
  pc->filter_level = vp9_read_literal(&header_bc, 6);
  pc->sharpness_level = vp9_read_literal(&header_bc, 3);

  /* Read in loop filter deltas applied at the MB level based on mode or ref frame. */
  xd->mode_ref_lf_delta_update = 0;
  xd->mode_ref_lf_delta_enabled = (unsigned char)vp9_read_bit(&header_bc);

  if (xd->mode_ref_lf_delta_enabled) {
    /* Do the deltas need to be updated */
    xd->mode_ref_lf_delta_update = (unsigned char)vp9_read_bit(&header_bc);

    if (xd->mode_ref_lf_delta_update) {
      /* Send update */
      for (i = 0; i < MAX_REF_LF_DELTAS; i++) {
        if (vp9_read_bit(&header_bc)) {
          /*sign = vp9_read_bit( &header_bc );*/
          xd->ref_lf_deltas[i] = (signed char)vp9_read_literal(&header_bc, 6);

          if (vp9_read_bit(&header_bc))        /* Apply sign */
            xd->ref_lf_deltas[i] = xd->ref_lf_deltas[i] * -1;
        }
      }

      /* Send update */
      for (i = 0; i < MAX_MODE_LF_DELTAS; i++) {
        if (vp9_read_bit(&header_bc)) {
          /*sign = vp9_read_bit( &header_bc );*/
          xd->mode_lf_deltas[i] = (signed char)vp9_read_literal(&header_bc, 6);

          if (vp9_read_bit(&header_bc))        /* Apply sign */
            xd->mode_lf_deltas[i] = xd->mode_lf_deltas[i] * -1;
        }
      }
    }
  }

  // Dummy read for now
  vp9_read_literal(&header_bc, 2);

  setup_token_decoder(pbi, data + first_partition_length_in_bytes,
                      &residual_bc);

  /* Read the default quantizers. */
  {
    int Q, q_update;

    Q = vp9_read_literal(&header_bc, QINDEX_BITS);
    pc->base_qindex = Q;
    q_update = 0;
    /* AC 1st order Q = default */
    pc->y1dc_delta_q = get_delta_q(&header_bc, pc->y1dc_delta_q, &q_update);
    pc->y2dc_delta_q = get_delta_q(&header_bc, pc->y2dc_delta_q, &q_update);
    pc->y2ac_delta_q = get_delta_q(&header_bc, pc->y2ac_delta_q, &q_update);
    pc->uvdc_delta_q = get_delta_q(&header_bc, pc->uvdc_delta_q, &q_update);
    pc->uvac_delta_q = get_delta_q(&header_bc, pc->uvac_delta_q, &q_update);

    if (q_update)
      vp9_init_de_quantizer(pbi);

    /* MB level dequantizer setup */
    mb_init_dequantizer(pbi, &pbi->mb);
  }

  /* Determine if the golden frame or ARF buffer should be updated and how.
   * For all non key frames the GF and ARF refresh flags and sign bias
   * flags must be set explicitly.
   */
  if (pc->frame_type != KEY_FRAME) {
    /* Should the GF or ARF be updated from the current frame */
    pc->refresh_golden_frame = vp9_read_bit(&header_bc);
    pc->refresh_alt_ref_frame = vp9_read_bit(&header_bc);

    if (pc->refresh_alt_ref_frame) {
      vpx_memcpy(&pc->fc, &pc->lfc_a, sizeof(pc->fc));
      vpx_memcpy(pc->fc.vp8_mode_contexts,
                 pc->fc.mode_context_a,
                 sizeof(pc->fc.vp8_mode_contexts));
    } else {
      vpx_memcpy(&pc->fc, &pc->lfc, sizeof(pc->fc));
      vpx_memcpy(pc->fc.vp8_mode_contexts,
                 pc->fc.mode_context,
                 sizeof(pc->fc.vp8_mode_contexts));
    }

    /* Buffer to buffer copy flags. */
    pc->copy_buffer_to_gf = 0;

    if (!pc->refresh_golden_frame)
      pc->copy_buffer_to_gf = vp9_read_literal(&header_bc, 2);

    pc->copy_buffer_to_arf = 0;

    if (!pc->refresh_alt_ref_frame)
      pc->copy_buffer_to_arf = vp9_read_literal(&header_bc, 2);

    pc->ref_frame_sign_bias[GOLDEN_FRAME] = vp9_read_bit(&header_bc);
    pc->ref_frame_sign_bias[ALTREF_FRAME] = vp9_read_bit(&header_bc);

    /* Is high precision mv allowed */
    xd->allow_high_precision_mv = (unsigned char)vp9_read_bit(&header_bc);
    // Read the type of subpel filter to use
    if (vp9_read_bit(&header_bc)) {
      pc->mcomp_filter_type = SWITCHABLE;
    } else {
      pc->mcomp_filter_type = vp9_read_literal(&header_bc, 2);
    }
    /* To enable choice of different interploation filters */
    vp9_setup_interp_filters(xd, pc->mcomp_filter_type, pc);
  }

  pc->refresh_entropy_probs = vp9_read_bit(&header_bc);
  if (pc->refresh_entropy_probs == 0) {
    vpx_memcpy(&pc->lfc, &pc->fc, sizeof(pc->fc));
  }

  pc->refresh_last_frame = (pc->frame_type == KEY_FRAME)
                           || vp9_read_bit(&header_bc);

  if (0) {
    FILE *z = fopen("decodestats.stt", "a");
    fprintf(z, "%6d F:%d,G:%d,A:%d,L:%d,Q:%d\n",
            pc->current_video_frame,
            pc->frame_type,
            pc->refresh_golden_frame,
            pc->refresh_alt_ref_frame,
            pc->refresh_last_frame,
            pc->base_qindex);
    fclose(z);
  }

  vp9_copy(pbi->common.fc.pre_coef_probs,
           pbi->common.fc.coef_probs);
  vp9_copy(pbi->common.fc.pre_hybrid_coef_probs,
           pbi->common.fc.hybrid_coef_probs);
  vp9_copy(pbi->common.fc.pre_coef_probs_8x8,
           pbi->common.fc.coef_probs_8x8);
  vp9_copy(pbi->common.fc.pre_hybrid_coef_probs_8x8,
           pbi->common.fc.hybrid_coef_probs_8x8);
  vp9_copy(pbi->common.fc.pre_coef_probs_16x16,
           pbi->common.fc.coef_probs_16x16);
  vp9_copy(pbi->common.fc.pre_hybrid_coef_probs_16x16,
           pbi->common.fc.hybrid_coef_probs_16x16);
  vp9_copy(pbi->common.fc.pre_ymode_prob, pbi->common.fc.ymode_prob);
  vp9_copy(pbi->common.fc.pre_uv_mode_prob, pbi->common.fc.uv_mode_prob);
  vp9_copy(pbi->common.fc.pre_bmode_prob, pbi->common.fc.bmode_prob);
  vp9_copy(pbi->common.fc.pre_i8x8_mode_prob, pbi->common.fc.i8x8_mode_prob);
  vp9_copy(pbi->common.fc.pre_sub_mv_ref_prob, pbi->common.fc.sub_mv_ref_prob);
  vp9_copy(pbi->common.fc.pre_mbsplit_prob, pbi->common.fc.mbsplit_prob);
  pbi->common.fc.pre_nmvc = pbi->common.fc.nmvc;
  vp9_zero(pbi->common.fc.coef_counts);
  vp9_zero(pbi->common.fc.hybrid_coef_counts);
  vp9_zero(pbi->common.fc.coef_counts_8x8);
  vp9_zero(pbi->common.fc.hybrid_coef_counts_8x8);
  vp9_zero(pbi->common.fc.coef_counts_16x16);
  vp9_zero(pbi->common.fc.hybrid_coef_counts_16x16);
  vp9_zero(pbi->common.fc.ymode_counts);
  vp9_zero(pbi->common.fc.uv_mode_counts);
  vp9_zero(pbi->common.fc.bmode_counts);
  vp9_zero(pbi->common.fc.i8x8_mode_counts);
  vp9_zero(pbi->common.fc.sub_mv_ref_counts);
  vp9_zero(pbi->common.fc.mbsplit_counts);
  vp9_zero(pbi->common.fc.NMVcount);
  vp9_zero(pbi->common.fc.mv_ref_ct);
  vp9_zero(pbi->common.fc.mv_ref_ct_a);

  read_coef_probs(pbi, &header_bc);

  vpx_memcpy(&xd->pre, &pc->yv12_fb[pc->lst_fb_idx], sizeof(YV12_BUFFER_CONFIG));
  vpx_memcpy(&xd->dst, &pc->yv12_fb[pc->new_fb_idx], sizeof(YV12_BUFFER_CONFIG));

  // Create the segmentation map structure and set to 0
  if (!pc->last_frame_seg_map)
    CHECK_MEM_ERROR(pc->last_frame_seg_map,
                    vpx_calloc((pc->mb_rows * pc->mb_cols), 1));

  /* set up frame new frame for intra coded blocks */
  vp9_setup_intra_recon(&pc->yv12_fb[pc->new_fb_idx]);

  vp9_setup_block_dptrs(xd);

  vp9_build_block_doffsets(xd);

  /* clear out the coeff buffer */
  vpx_memset(xd->qcoeff, 0, sizeof(xd->qcoeff));

  /* Read the mb_no_coeff_skip flag */
  pc->mb_no_coeff_skip = (int)vp9_read_bit(&header_bc);

  vp9_decode_mode_mvs_init(pbi, &header_bc);

  vpx_memset(pc->above_context, 0, sizeof(ENTROPY_CONTEXT_PLANES) * pc->mb_cols);

  // Resset the macroblock mode info context to the start of the list
  xd->mode_info_context = pc->mi;
  xd->prev_mode_info_context = pc->prev_mi;

  /* Decode a row of superblocks */
  for (mb_row = 0; mb_row < pc->mb_rows; mb_row += 2) {
    decode_sb_row(pbi, pc, mb_row, xd, &residual_bc);
  }
  corrupt_tokens |= xd->corrupted;

  /* Collect information about decoder corruption. */
  /* 1. Check first boolean decoder for errors. */
  pc->yv12_fb[pc->new_fb_idx].corrupted = bool_error(&header_bc);
  /* 2. Check the macroblock information */
  pc->yv12_fb[pc->new_fb_idx].corrupted |= corrupt_tokens;

  if (!pbi->decoded_key_frame) {
    if (pc->frame_type == KEY_FRAME &&
        !pc->yv12_fb[pc->new_fb_idx].corrupted)
      pbi->decoded_key_frame = 1;
    else
      vpx_internal_error(&pbi->common.error, VPX_CODEC_CORRUPT_FRAME,
                         "A stream must start with a complete key frame");
  }

  vp9_adapt_coef_probs(pc);
  if (pc->frame_type != KEY_FRAME) {
    vp9_adapt_mode_probs(pc);
    vp9_adapt_nmv_probs(pc, xd->allow_high_precision_mv);
    vp9_update_mode_context(&pbi->common);
  }

  /* If this was a kf or Gf note the Q used */
  if ((pc->frame_type == KEY_FRAME) ||
      pc->refresh_golden_frame || pc->refresh_alt_ref_frame) {
    pc->last_kf_gf_q = pc->base_qindex;
  }
  if (pc->refresh_entropy_probs) {
    if (pc->refresh_alt_ref_frame)
      vpx_memcpy(&pc->lfc_a, &pc->fc, sizeof(pc->fc));
    else
      vpx_memcpy(&pc->lfc, &pc->fc, sizeof(pc->fc));
  }

#ifdef PACKET_TESTING
  {
    FILE *f = fopen("decompressor.VP8", "ab");
    unsigned int size = residual_bc.pos + header_bc.pos + 8;
    fwrite((void *) &size, 4, 1, f);
    fwrite((void *) pbi->Source, size, 1, f);
    fclose(f);
  }
#endif
  // printf("Frame %d Done\n", frame_count++);

  return 0;
}
