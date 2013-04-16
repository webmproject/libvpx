/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */


#include "vp9/decoder/vp9_onyxd_int.h"
#include "vp9/common/vp9_common.h"
#include "vp9/common/vp9_header.h"
#include "vp9/common/vp9_reconintra.h"
#include "vp9/common/vp9_reconinter.h"
#include "vp9/common/vp9_entropy.h"
#include "vp9/decoder/vp9_decodframe.h"
#include "vp9/decoder/vp9_detokenize.h"
#include "vp9/common/vp9_invtrans.h"
#include "vp9/common/vp9_alloccommon.h"
#include "vp9/common/vp9_entropymode.h"
#include "vp9/common/vp9_quant_common.h"
#include "vpx_scale/vpx_scale.h"
#include "vp9/common/vp9_setupintrarecon.h"

#include "vp9/decoder/vp9_decodemv.h"
#include "vp9/common/vp9_extend.h"
#include "vp9/common/vp9_modecont.h"
#include "vpx_mem/vpx_mem.h"
#include "vp9/decoder/vp9_dboolhuff.h"

#include "vp9/common/vp9_seg_common.h"
#include "vp9/common/vp9_tile_common.h"
#include "vp9_rtcd.h"

#include <assert.h>
#include <stdio.h>

#define COEFCOUNT_TESTING

// #define DEC_DEBUG
#ifdef DEC_DEBUG
int dec_debug = 0;
#endif

static int read_le16(const uint8_t *p) {
  return (p[1] << 8) | p[0];
}

static int read_le32(const uint8_t *p) {
  return (p[3] << 24) | (p[2] << 16) | (p[1] << 8) | p[0];
}

// len == 0 is not allowed
static int read_is_valid(const uint8_t *start, size_t len,
                         const uint8_t *end) {
  return start + len > start && start + len <= end;
}

static TXFM_MODE read_txfm_mode(vp9_reader *r) {
  TXFM_MODE mode = vp9_read_literal(r, 2);
  if (mode == ALLOW_32X32)
    mode += vp9_read_bit(r);
  return mode;
}

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

  v = merge_index(v, n - 1, modulus);
  if ((m << 1) <= n) {
    return vp9_inv_recenter_nonneg(v + 1, m);
  } else {
    return n - 1 - vp9_inv_recenter_nonneg(v + 1, n - 1 - m);
  }
}

static vp9_prob read_prob_diff_update(vp9_reader *const bc, int oldp) {
  int delp = vp9_decode_term_subexp(bc, SUBEXP_PARAM, 255);
  return (vp9_prob)inv_remap_prob(delp, oldp);
}

void vp9_init_de_quantizer(VP9D_COMP *pbi) {
  int i;
  int q;
  VP9_COMMON *const pc = &pbi->common;

  for (q = 0; q < QINDEX_RANGE; q++) {
    pc->y_dequant[q][0] = (int16_t)vp9_dc_quant(q, pc->y_dc_delta_q);
    pc->uv_dequant[q][0] = (int16_t)vp9_dc_uv_quant(q, pc->uv_dc_delta_q);

    /* all the ac values =; */
    for (i = 1; i < 16; i++) {
      const int rc = vp9_default_zig_zag1d_4x4[i];

      pc->y_dequant[q][rc] = (int16_t)vp9_ac_yquant(q);
      pc->uv_dequant[q][rc] = (int16_t)vp9_ac_uv_quant(q, pc->uv_ac_delta_q);
    }
  }
}

static int get_qindex(MACROBLOCKD *mb, int segment_id, int base_qindex) {
  // Set the Q baseline allowing for any segment level adjustment
  if (vp9_segfeature_active(mb, segment_id, SEG_LVL_ALT_Q)) {
    const int data = vp9_get_segdata(mb, segment_id, SEG_LVL_ALT_Q);
    return mb->mb_segment_abs_delta == SEGMENT_ABSDATA ?
               data :  // Abs value
               clamp(base_qindex + data, 0, MAXQ);  // Delta value
  } else {
    return base_qindex;
  }
}

static void mb_init_dequantizer(VP9D_COMP *pbi, MACROBLOCKD *mb) {
  int i;

  VP9_COMMON *const pc = &pbi->common;
  const int segment_id = mb->mode_info_context->mbmi.segment_id;
  const int qindex = get_qindex(mb, segment_id, pc->base_qindex);
  mb->q_index = qindex;

  for (i = 0; i < 16; i++)
    mb->block[i].dequant = pc->y_dequant[qindex];

  for (i = 16; i < 24; i++)
    mb->block[i].dequant = pc->uv_dequant[qindex];

  if (mb->lossless) {
    assert(qindex == 0);
    mb->inv_txm4x4_1      = vp9_short_iwalsh4x4_1;
    mb->inv_txm4x4        = vp9_short_iwalsh4x4;
    mb->itxm_add          = vp9_dequant_idct_add_lossless_c;
    mb->itxm_add_y_block  = vp9_dequant_idct_add_y_block_lossless_c;
    mb->itxm_add_uv_block = vp9_dequant_idct_add_uv_block_lossless_c;
  } else {
    mb->inv_txm4x4_1      = vp9_short_idct4x4_1;
    mb->inv_txm4x4        = vp9_short_idct4x4;
    mb->itxm_add          = vp9_dequant_idct_add;
    mb->itxm_add_y_block  = vp9_dequant_idct_add_y_block;
    mb->itxm_add_uv_block = vp9_dequant_idct_add_uv_block;
  }
}

#if CONFIG_CODE_NONZEROCOUNT
static void propagate_nzcs(VP9_COMMON *cm, MACROBLOCKD *xd) {
  MODE_INFO *m = xd->mode_info_context;
  BLOCK_SIZE_TYPE sb_type = m->mbmi.sb_type;
  const int mis = cm->mode_info_stride;
  int n;
  if (sb_type == BLOCK_SIZE_SB64X64) {
    for (n = 0; n < 16; ++n) {
      int i = n >> 2;
      int j = n & 3;
      if (i == 0 && j == 0) continue;
      vpx_memcpy((m + j + mis * i)->mbmi.nzcs, m->mbmi.nzcs,
                 384 * sizeof(m->mbmi.nzcs[0]));
    }
  } else if (sb_type == BLOCK_SIZE_SB32X32) {
    for (n = 0; n < 4; ++n) {
      int i = n >> 1;
      int j = n & 1;
      if (i == 0 && j == 0) continue;
      vpx_memcpy((m + j + mis * i)->mbmi.nzcs, m->mbmi.nzcs,
                 384 * sizeof(m->mbmi.nzcs[0]));
    }
  }
}
#endif

static void skip_recon_sb(VP9D_COMP *pbi, MACROBLOCKD *xd,
                          int mb_row, int mb_col,
                          BLOCK_SIZE_TYPE bsize) {
  MODE_INFO *m = xd->mode_info_context;

  if (m->mbmi.ref_frame == INTRA_FRAME) {
    vp9_build_intra_predictors_sbuv_s(xd, bsize);
    vp9_build_intra_predictors_sby_s(xd, bsize);
  } else {
    vp9_build_inter_predictors_sb(xd, mb_row, mb_col, bsize);
  }
#if CONFIG_CODE_NONZEROCOUNT
  vpx_memset(m->mbmi.nzcs, 0, 384 * sizeof(m->mbmi.nzcs[0]));
  propagate_nzcs(&pbi->common, xd);
#endif
}

static void decode_16x16(VP9D_COMP *pbi, MACROBLOCKD *xd,
                         BOOL_DECODER* const bc) {
  const TX_TYPE tx_type = get_tx_type_16x16(xd, 0);

  vp9_dequant_iht_add_16x16_c(tx_type, xd->plane[0].qcoeff,
                              xd->block[0].dequant, xd->dst.y_buffer,
                              xd->dst.y_stride, xd->plane[0].eobs[0]);

  vp9_dequant_idct_add_8x8(xd->plane[1].qcoeff, xd->block[16].dequant,
                           xd->dst.u_buffer, xd->dst.uv_stride,
                           xd->plane[1].eobs[0]);

  vp9_dequant_idct_add_8x8(xd->plane[2].qcoeff, xd->block[20].dequant,
                           xd->dst.v_buffer, xd->dst.uv_stride,
                           xd->plane[2].eobs[0]);
}

static void decode_8x8(VP9D_COMP *pbi, MACROBLOCKD *xd,
                       BOOL_DECODER* const bc) {
  const MB_PREDICTION_MODE mode = xd->mode_info_context->mbmi.mode;
  // luma
  // if the first one is DCT_DCT assume all the rest are as well
  TX_TYPE tx_type = get_tx_type_8x8(xd, 0);
  if (tx_type != DCT_DCT || mode == I8X8_PRED) {
    int i;
    for (i = 0; i < 4; i++) {
      int ib = vp9_i8x8_block[i];
      int idx = (ib & 0x02) ? (ib + 2) : ib;
      int16_t *q  = BLOCK_OFFSET(xd->plane[0].qcoeff, idx, 16);
      int16_t *dq = xd->block[0].dequant;
      uint8_t *dst = *(xd->block[ib].base_dst) + xd->block[ib].dst;
      int stride = xd->dst.y_stride;
      if (mode == I8X8_PRED) {
        BLOCKD *b = &xd->block[ib];
        int i8x8mode = b->bmi.as_mode.first;
        vp9_intra8x8_predict(xd, b, i8x8mode, dst, stride);
      }
      tx_type = get_tx_type_8x8(xd, ib);
      vp9_dequant_iht_add_8x8_c(tx_type, q, dq, dst, stride,
                                xd->plane[0].eobs[idx]);
    }
  } else {
    vp9_dequant_idct_add_y_block_8x8(xd->plane[0].qcoeff,
                                     xd->block[0].dequant, xd->dst.y_buffer,
                                     xd->dst.y_stride, xd);
  }

  // chroma
  if (mode == I8X8_PRED) {
    int i;
    for (i = 0; i < 4; i++) {
      int ib = vp9_i8x8_block[i];
      BLOCKD *b = &xd->block[ib];
      int i8x8mode = b->bmi.as_mode.first;

      b = &xd->block[16 + i];
      vp9_intra_uv4x4_predict(xd, b, i8x8mode, *(b->base_dst) + b->dst,
                              b->dst_stride);
      xd->itxm_add(BLOCK_OFFSET(xd->plane[1].qcoeff, i, 16),
                   b->dequant, *(b->base_dst) + b->dst, b->dst_stride,
                   xd->plane[1].eobs[i]);

      b = &xd->block[20 + i];
      vp9_intra_uv4x4_predict(xd, b, i8x8mode, *(b->base_dst) + b->dst,
                              b->dst_stride);
      xd->itxm_add(BLOCK_OFFSET(xd->plane[2].qcoeff, i, 16),
                   b->dequant, *(b->base_dst) + b->dst, b->dst_stride,
                   xd->plane[2].eobs[i]);
    }
  } else if (mode == SPLITMV) {
    xd->itxm_add_uv_block(xd->plane[1].qcoeff, xd->block[16].dequant,
         xd->dst.u_buffer, xd->dst.uv_stride, xd->plane[1].eobs);
    xd->itxm_add_uv_block(xd->plane[2].qcoeff, xd->block[16].dequant,
         xd->dst.v_buffer, xd->dst.uv_stride, xd->plane[2].eobs);
  } else {
    vp9_dequant_idct_add_8x8(xd->plane[1].qcoeff, xd->block[16].dequant,
                             xd->dst.u_buffer, xd->dst.uv_stride,
                             xd->plane[1].eobs[0]);

    vp9_dequant_idct_add_8x8(xd->plane[2].qcoeff, xd->block[16].dequant,
                             xd->dst.v_buffer, xd->dst.uv_stride,
                             xd->plane[2].eobs[0]);
  }
}

static INLINE void dequant_add_y(MACROBLOCKD *xd, TX_TYPE tx_type, int idx) {
  BLOCKD *const b = &xd->block[idx];
  struct mb_plane *const y = &xd->plane[0];
  if (tx_type != DCT_DCT) {
    vp9_dequant_iht_add_c(tx_type,
                          BLOCK_OFFSET(y->qcoeff, idx, 16),
                          b->dequant, *(b->base_dst) + b->dst,
                          b->dst_stride, y->eobs[idx]);
  } else {
    xd->itxm_add(BLOCK_OFFSET(y->qcoeff, idx, 16),
                 b->dequant, *(b->base_dst) + b->dst,
                 b->dst_stride, y->eobs[idx]);
  }
}


static void decode_4x4(VP9D_COMP *pbi, MACROBLOCKD *xd,
                       BOOL_DECODER* const bc) {
  TX_TYPE tx_type;
  int i = 0;
  const MB_PREDICTION_MODE mode = xd->mode_info_context->mbmi.mode;
  if (mode == I8X8_PRED) {
    for (i = 0; i < 4; i++) {
      int ib = vp9_i8x8_block[i];
      const int iblock[4] = {0, 1, 4, 5};
      int j;
      BLOCKD *b = &xd->block[ib];
      int i8x8mode = b->bmi.as_mode.first;
      vp9_intra8x8_predict(xd, b, i8x8mode, *(b->base_dst) + b->dst,
                           b->dst_stride);
      for (j = 0; j < 4; j++) {
        tx_type = get_tx_type_4x4(xd, ib + iblock[j]);
        dequant_add_y(xd, tx_type, ib + iblock[j]);
      }
      b = &xd->block[16 + i];
      vp9_intra_uv4x4_predict(xd, b, i8x8mode, *(b->base_dst) + b->dst,
                              b->dst_stride);
      xd->itxm_add(BLOCK_OFFSET(xd->plane[1].qcoeff, i, 16),
                   b->dequant, *(b->base_dst) + b->dst, b->dst_stride,
                   xd->plane[1].eobs[i]);
      b = &xd->block[20 + i];
      vp9_intra_uv4x4_predict(xd, b, i8x8mode, *(b->base_dst) + b->dst,
                              b->dst_stride);
      xd->itxm_add(BLOCK_OFFSET(xd->plane[2].qcoeff, i, 16),
                   b->dequant, *(b->base_dst) + b->dst, b->dst_stride,
                   xd->plane[2].eobs[i]);
    }
  } else if (mode == I4X4_PRED) {
    for (i = 0; i < 16; i++) {
      BLOCKD *b = &xd->block[i];
      int b_mode = xd->mode_info_context->bmi[i].as_mode.first;
#if CONFIG_NEWBINTRAMODES
      xd->mode_info_context->bmi[i].as_mode.context = b->bmi.as_mode.context =
          vp9_find_bpred_context(xd, b);
      if (!xd->mode_info_context->mbmi.mb_skip_coeff)
        vp9_decode_coefs_4x4(pbi, xd, bc, PLANE_TYPE_Y_WITH_DC, i);
#endif
      vp9_intra4x4_predict(xd, b, b_mode, *(b->base_dst) + b->dst,
                           b->dst_stride);
      tx_type = get_tx_type_4x4(xd, i);
      dequant_add_y(xd, tx_type, i);
    }
#if CONFIG_NEWBINTRAMODES
    if (!xd->mode_info_context->mbmi.mb_skip_coeff)
      vp9_decode_mb_tokens_4x4_uv(pbi, xd, bc);
#endif
    vp9_build_intra_predictors_sbuv_s(xd, BLOCK_SIZE_MB16X16);
    xd->itxm_add_uv_block(xd->plane[1].qcoeff, xd->block[16].dequant,
         xd->dst.u_buffer, xd->dst.uv_stride, xd->plane[1].eobs);
    xd->itxm_add_uv_block(xd->plane[2].qcoeff, xd->block[16].dequant,
         xd->dst.v_buffer, xd->dst.uv_stride, xd->plane[2].eobs);
  } else if (mode == SPLITMV || get_tx_type_4x4(xd, 0) == DCT_DCT) {
    xd->itxm_add_y_block(xd->plane[0].qcoeff,
                          xd->block[0].dequant,
                          xd->dst.y_buffer, xd->dst.y_stride, xd);
    xd->itxm_add_uv_block(xd->plane[1].qcoeff, xd->block[16].dequant,
         xd->dst.u_buffer, xd->dst.uv_stride, xd->plane[1].eobs);
    xd->itxm_add_uv_block(xd->plane[2].qcoeff, xd->block[16].dequant,
         xd->dst.v_buffer, xd->dst.uv_stride, xd->plane[2].eobs);
  } else {
    for (i = 0; i < 16; i++) {
      tx_type = get_tx_type_4x4(xd, i);
      dequant_add_y(xd, tx_type, i);
    }
    xd->itxm_add_uv_block(xd->plane[1].qcoeff, xd->block[16].dequant,
                          xd->dst.u_buffer, xd->dst.uv_stride,
                          xd->plane[1].eobs);
    xd->itxm_add_uv_block(xd->plane[2].qcoeff, xd->block[16].dequant,
                          xd->dst.v_buffer, xd->dst.uv_stride,
                          xd->plane[2].eobs);
  }
}

static INLINE void decode_sby_32x32(MACROBLOCKD *mb, BLOCK_SIZE_TYPE bsize) {
  const int bwl = mb_width_log2(bsize) - 1, bw = 1 << bwl;
  const int bhl = mb_height_log2(bsize) - 1, bh = 1 << bhl;
  const int y_count = bw * bh;
  int n;

  for (n = 0; n < y_count; n++) {
    const int x_idx = n & (bw - 1);
    const int y_idx = n >> bwl;
    const int y_offset = (y_idx * 32) * mb->dst.y_stride + (x_idx * 32);
    vp9_dequant_idct_add_32x32(BLOCK_OFFSET(mb->plane[0].qcoeff, n, 1024),
                               mb->block[0].dequant ,
                               mb->dst.y_buffer + y_offset, mb->dst.y_stride,
                               mb->plane[0].eobs[n * 64]);
  }
}

static INLINE void decode_sbuv_32x32(MACROBLOCKD *mb, BLOCK_SIZE_TYPE bsize) {
  const int bwl = mb_width_log2(bsize) - 1, bw = (1 << bwl) / 2;
  const int bhl = mb_height_log2(bsize) - 1, bh = (1 << bhl) / 2;
  const int uv_count = bw * bh;
  int n;
  for (n = 0; n < uv_count; n++) {
     const int x_idx = n & (bw - 1);
     const int y_idx = n >> (bwl - 1);
     const int uv_offset = (y_idx * 32) * mb->dst.uv_stride + (x_idx * 32);
     vp9_dequant_idct_add_32x32(BLOCK_OFFSET(mb->plane[1].qcoeff, n, 1024),
                                mb->block[16].dequant,
                                mb->dst.u_buffer + uv_offset,
                                mb->dst.uv_stride, mb->plane[1].eobs[n * 64]);
     vp9_dequant_idct_add_32x32(BLOCK_OFFSET(mb->plane[2].qcoeff, n, 1024),
                                mb->block[20].dequant,
                                mb->dst.v_buffer + uv_offset,
                                mb->dst.uv_stride, mb->plane[2].eobs[n * 64]);
  }
}

static INLINE void decode_sby_16x16(MACROBLOCKD *mb, BLOCK_SIZE_TYPE bsize) {
  const int bwl = mb_width_log2(bsize), bw = 1 << bwl;
  const int bhl = mb_height_log2(bsize), bh = 1 << bhl;
  const int y_count = bw * bh;
  int n;

  for (n = 0; n < y_count; n++) {
    const int x_idx = n & (bw - 1);
    const int y_idx = n >> bwl;
    const int y_offset = (y_idx * 16) * mb->dst.y_stride + (x_idx * 16);
    const TX_TYPE tx_type = get_tx_type_16x16(mb,
                                (y_idx * (4 * bw) + x_idx) * 4);
    vp9_dequant_iht_add_16x16_c(tx_type,
                                BLOCK_OFFSET(mb->plane[0].qcoeff, n, 256),
                                mb->block[0].dequant,
                                mb->dst.y_buffer + y_offset,
                                mb->dst.y_stride,
                                mb->plane[0].eobs[n * 16]);
  }
}

static INLINE void decode_sbuv_16x16(MACROBLOCKD *mb, BLOCK_SIZE_TYPE bsize) {
  const int bwl = mb_width_log2(bsize), bw = (1 << bwl) / 2;
  const int bhl = mb_height_log2(bsize), bh = (1 << bhl) / 2;
  const int uv_count = bw * bh;
  int n;

  assert(bsize >= BLOCK_SIZE_SB32X32);

  for (n = 0; n < uv_count; n++) {
    const int x_idx = n & (bw - 1);
    const int y_idx = n >> (bwl - 1);
    const int uv_offset = (y_idx * 16) * mb->dst.uv_stride + (x_idx * 16);
    vp9_dequant_idct_add_16x16(BLOCK_OFFSET(mb->plane[1].qcoeff, n, 256),
                               mb->block[16].dequant,
                               mb->dst.u_buffer + uv_offset, mb->dst.uv_stride,
                               mb->plane[1].eobs[n * 16]);
    vp9_dequant_idct_add_16x16(BLOCK_OFFSET(mb->plane[2].qcoeff, n, 256),
                               mb->block[20].dequant,
                               mb->dst.v_buffer + uv_offset, mb->dst.uv_stride,
                               mb->plane[2].eobs[n * 16]);
  }
}

static INLINE void decode_sby_8x8(MACROBLOCKD *xd, BLOCK_SIZE_TYPE bsize) {
  const int bwl = mb_width_log2(bsize)  + 1, bw = 1 << bwl;
  const int bhl = mb_height_log2(bsize) + 1, bh = 1 << bhl;
  const int y_count = bw * bh;
  int n;

  // luma
  for (n = 0; n < y_count; n++) {
    const int x_idx = n & (bw - 1);
    const int y_idx = n >> bwl;
    const int y_offset = (y_idx * 8) * xd->dst.y_stride + (x_idx * 8);
    const TX_TYPE tx_type = get_tx_type_8x8(xd,
                                            (y_idx * (2 * bw) + x_idx) * 2);

    vp9_dequant_iht_add_8x8_c(tx_type,
                              BLOCK_OFFSET(xd->plane[0].qcoeff, n, 64),
                              xd->block[0].dequant,
                              xd->dst.y_buffer + y_offset, xd->dst.y_stride,
                              xd->plane[0].eobs[n * 4]);
  }
}

static INLINE void decode_sbuv_8x8(MACROBLOCKD *xd, BLOCK_SIZE_TYPE bsize) {
  const int bwl = mb_width_log2(bsize)  + 1, bw = 1 << (bwl - 1);
  const int bhl = mb_height_log2(bsize) + 1, bh = 1 << (bhl - 1);
  const int uv_count = bw * bh;
  int n;

  // chroma
  for (n = 0; n < uv_count; n++) {
    const int x_idx = n & (bw - 1);
    const int y_idx = n >> (bwl - 1);
    const int uv_offset = (y_idx * 8) * xd->dst.uv_stride + (x_idx * 8);
    vp9_dequant_idct_add_8x8(BLOCK_OFFSET(xd->plane[1].qcoeff, n, 64),
                             xd->block[16].dequant,
                             xd->dst.u_buffer + uv_offset, xd->dst.uv_stride,
                             xd->plane[1].eobs[n * 4]);
    vp9_dequant_idct_add_8x8(BLOCK_OFFSET(xd->plane[2].qcoeff, n, 64),
                             xd->block[20].dequant,
                             xd->dst.v_buffer + uv_offset, xd->dst.uv_stride,
                             xd->plane[2].eobs[n * 4]);
  }
}

static INLINE void decode_sby_4x4(MACROBLOCKD *xd, BLOCK_SIZE_TYPE bsize) {
  const int bwl = mb_width_log2(bsize)  + 2, bw = 1 << bwl;
  const int bhl = mb_height_log2(bsize) + 2, bh = 1 << bhl;
  const int y_count = bw * bh;
  int n;

  for (n = 0; n < y_count; n++) {
    const int x_idx = n & (bw - 1);
    const int y_idx = n >> bwl;
    const int y_offset = (y_idx * 4) * xd->dst.y_stride + (x_idx * 4);
    const TX_TYPE tx_type = get_tx_type_4x4(xd, n);
    if (tx_type == DCT_DCT) {
      xd->itxm_add(BLOCK_OFFSET(xd->plane[0].qcoeff, n, 16),
                   xd->block[0].dequant,
                   xd->dst.y_buffer + y_offset, xd->dst.y_stride,
                   xd->plane[0].eobs[n]);
    } else {
      vp9_dequant_iht_add_c(tx_type,
                            BLOCK_OFFSET(xd->plane[0].qcoeff, n, 16),
                            xd->block[0].dequant, xd->dst.y_buffer + y_offset,
                            xd->dst.y_stride, xd->plane[0].eobs[n]);
    }
  }
}

static INLINE void decode_sbuv_4x4(MACROBLOCKD *xd, BLOCK_SIZE_TYPE bsize) {
  const int bwl = mb_width_log2(bsize)  + 2, bw = 1 << (bwl - 1);
  const int bhl = mb_height_log2(bsize) + 2, bh = 1 << (bhl - 1);
  const int uv_count = bw * bh;
  int n;

  for (n = 0; n < uv_count; n++) {
    const int x_idx = n & (bw - 1);
    const int y_idx = n >> (bwl - 1);
    const int uv_offset = (y_idx * 4) * xd->dst.uv_stride + (x_idx * 4);
    xd->itxm_add(BLOCK_OFFSET(xd->plane[1].qcoeff, n, 16),
        xd->block[16].dequant,
        xd->dst.u_buffer + uv_offset, xd->dst.uv_stride, xd->plane[1].eobs[n]);
    xd->itxm_add(BLOCK_OFFSET(xd->plane[2].qcoeff, n, 16),
        xd->block[20].dequant,
        xd->dst.v_buffer + uv_offset, xd->dst.uv_stride, xd->plane[2].eobs[n]);
  }
}

// TODO(jingning): combine luma and chroma dequantization and inverse
// transform into a single function looping over planes.
static void decode_sb_32x32(MACROBLOCKD *mb, BLOCK_SIZE_TYPE bsize) {
  decode_sby_32x32(mb, bsize);
  if (bsize == BLOCK_SIZE_SB64X64)
    decode_sbuv_32x32(mb, bsize);
  else
    decode_sbuv_16x16(mb, bsize);
}

static void decode_sb_16x16(MACROBLOCKD *mb, BLOCK_SIZE_TYPE bsize) {
  decode_sby_16x16(mb, bsize);
  if (bsize >= BLOCK_SIZE_SB32X32)
    decode_sbuv_16x16(mb, bsize);
  else
    decode_sbuv_8x8(mb, bsize);
}

static void decode_sb(VP9D_COMP *pbi, MACROBLOCKD *xd, int mb_row, int mb_col,
                      BOOL_DECODER* const bc, BLOCK_SIZE_TYPE bsize) {
  const int bwl = mb_width_log2(bsize), bhl = mb_height_log2(bsize);
  const int bw = 1 << bwl, bh = 1 << bhl;
  int n, eobtotal;
  VP9_COMMON *const pc = &pbi->common;
  MODE_INFO *mi = xd->mode_info_context;
  const int mis = pc->mode_info_stride;

  assert(mi->mbmi.sb_type == bsize);

  if (pbi->common.frame_type != KEY_FRAME)
    vp9_setup_interp_filters(xd, mi->mbmi.interp_filter, pc);

  // re-initialize macroblock dequantizer before detokenization
  if (xd->segmentation_enabled)
    mb_init_dequantizer(pbi, xd);

  if (mi->mbmi.mb_skip_coeff) {
    vp9_reset_sb_tokens_context(xd, bsize);

    // Special case:  Force the loopfilter to skip when eobtotal and
    // mb_skip_coeff are zero.
    skip_recon_sb(pbi, xd, mb_row, mb_col, bsize);
    return;
  }

  // generate prediction
  if (xd->mode_info_context->mbmi.ref_frame == INTRA_FRAME) {
    vp9_build_intra_predictors_sby_s(xd, bsize);
    vp9_build_intra_predictors_sbuv_s(xd, bsize);
  } else {
    vp9_build_inter_predictors_sb(xd, mb_row, mb_col, bsize);
  }

  // dequantization and idct
  eobtotal = vp9_decode_tokens(pbi, xd, bc, bsize);
  if (eobtotal == 0) {  // skip loopfilter
    for (n = 0; n < bw * bh; n++) {
      const int x_idx = n & (bw - 1), y_idx = n >> bwl;

      if (mb_col + x_idx < pc->mb_cols && mb_row + y_idx < pc->mb_rows)
        mi[y_idx * mis + x_idx].mbmi.mb_skip_coeff = mi->mbmi.mb_skip_coeff;
    }
  } else {
    switch (xd->mode_info_context->mbmi.txfm_size) {
      case TX_32X32:
        decode_sb_32x32(xd, bsize);
        break;
      case TX_16X16:
        decode_sb_16x16(xd, bsize);
        break;
      case TX_8X8:
        decode_sby_8x8(xd, bsize);
        decode_sbuv_8x8(xd, bsize);
        break;
      case TX_4X4:
        decode_sby_4x4(xd, bsize);
        decode_sbuv_4x4(xd, bsize);
        break;
      default: assert(0);
    }
  }
#if CONFIG_CODE_NONZEROCOUNT
  propagate_nzcs(&pbi->common, xd);
#endif
}

// TODO(jingning): Need to merge SB and MB decoding. The MB decoding currently
// couples special handles on I8x8, B_PRED, and splitmv modes.
static void decode_mb(VP9D_COMP *pbi, MACROBLOCKD *xd,
                     int mb_row, int mb_col,
                     BOOL_DECODER* const bc) {
  int eobtotal = 0;
  const MB_PREDICTION_MODE mode = xd->mode_info_context->mbmi.mode;
  const int tx_size = xd->mode_info_context->mbmi.txfm_size;

  assert(!xd->mode_info_context->mbmi.sb_type);

  // re-initialize macroblock dequantizer before detokenization
  if (xd->segmentation_enabled)
    mb_init_dequantizer(pbi, xd);

  if (xd->mode_info_context->mbmi.mb_skip_coeff) {
    vp9_reset_sb_tokens_context(xd, BLOCK_SIZE_MB16X16);
  } else if (!bool_error(bc)) {
#if CONFIG_NEWBINTRAMODES
    if (mode != I4X4_PRED)
#endif
      eobtotal = vp9_decode_tokens(pbi, xd, bc, BLOCK_SIZE_MB16X16);
  }

  //mode = xd->mode_info_context->mbmi.mode;
  if (pbi->common.frame_type != KEY_FRAME)
    vp9_setup_interp_filters(xd, xd->mode_info_context->mbmi.interp_filter,
                             &pbi->common);

  if (eobtotal == 0 &&
      mode != I4X4_PRED &&
      mode != SPLITMV &&
      mode != I8X8_PRED &&
      !bool_error(bc)) {
    // Special case:  Force the loopfilter to skip when eobtotal and
    // mb_skip_coeff are zero.
    xd->mode_info_context->mbmi.mb_skip_coeff = 1;
    skip_recon_sb(pbi, xd, mb_row, mb_col, BLOCK_SIZE_MB16X16);
    return;
  }
#if 0  // def DEC_DEBUG
  if (dec_debug)
    printf("Decoding mb:  %d %d\n", xd->mode_info_context->mbmi.mode, tx_size);
#endif

  // moved to be performed before detokenization
  //  if (xd->segmentation_enabled)
  //    mb_init_dequantizer(pbi, xd);

  // do prediction
  if (xd->mode_info_context->mbmi.ref_frame == INTRA_FRAME) {
    if (mode != I8X8_PRED) {
      vp9_build_intra_predictors_sbuv_s(xd, BLOCK_SIZE_MB16X16);
      if (mode != I4X4_PRED)
        vp9_build_intra_predictors_sby_s(xd, BLOCK_SIZE_MB16X16);
    }
  } else {
#if 0  // def DEC_DEBUG
  if (dec_debug)
    printf("Decoding mb:  %d %d interp %d\n",
           xd->mode_info_context->mbmi.mode, tx_size,
           xd->mode_info_context->mbmi.interp_filter);
#endif
    vp9_build_inter_predictors_mb_s(xd, mb_row, mb_col);
  }

  if (tx_size == TX_16X16) {
    decode_16x16(pbi, xd, bc);
  } else if (tx_size == TX_8X8) {
    decode_8x8(pbi, xd, bc);
  } else {
    decode_4x4(pbi, xd, bc);
  }
#ifdef DEC_DEBUG
  if (dec_debug) {
    int i, j;
    printf("\n");
    printf("predictor y\n");
    for (i = 0; i < 16; i++) {
      for (j = 0; j < 16; j++)
        printf("%3d ", xd->predictor[i * 16 + j]);
      printf("\n");
    }
    printf("\n");
    printf("final y\n");
    for (i = 0; i < 16; i++) {
      for (j = 0; j < 16; j++)
        printf("%3d ", xd->dst.y_buffer[i * xd->dst.y_stride + j]);
      printf("\n");
    }
    printf("\n");
    printf("final u\n");
    for (i = 0; i < 8; i++) {
      for (j = 0; j < 8; j++)
        printf("%3d ", xd->dst.u_buffer[i * xd->dst.uv_stride + j]);
      printf("\n");
    }
    printf("\n");
    printf("final v\n");
    for (i = 0; i < 8; i++) {
      for (j = 0; j < 8; j++)
        printf("%3d ", xd->dst.v_buffer[i * xd->dst.uv_stride + j]);
      printf("\n");
    }
    fflush(stdout);
  }
#endif
}

static int get_delta_q(vp9_reader *r, int *dq) {
  const int old_value = *dq;

  if (vp9_read_bit(r)) {  // Update bit
    const int value = vp9_read_literal(r, 4);
    *dq = vp9_read_and_apply_sign(r, value);
  }

  // Trigger a quantizer update if the delta-q value has changed
  return old_value != *dq;
}

#ifdef PACKET_TESTING
#include <stdio.h>
FILE *vpxlog = 0;
#endif

static void set_offsets(VP9D_COMP *pbi, BLOCK_SIZE_TYPE bsize,
                        int mb_row, int mb_col) {
  const int bh = 1 << mb_height_log2(bsize);
  const int bw = 1 << mb_width_log2(bsize);
  VP9_COMMON *const cm = &pbi->common;
  MACROBLOCKD *const xd = &pbi->mb;

  const int mb_idx = mb_row * cm->mode_info_stride + mb_col;
  const YV12_BUFFER_CONFIG *dst_fb = &cm->yv12_fb[cm->new_fb_idx];
  const int recon_yoffset = (16 * mb_row) * dst_fb->y_stride + (16 * mb_col);
  const int recon_uvoffset = (8 * mb_row) * dst_fb->uv_stride + (8 * mb_col);

  xd->mode_info_context = cm->mi + mb_idx;
  xd->mode_info_context->mbmi.sb_type = bsize;
  xd->prev_mode_info_context = cm->prev_mi + mb_idx;
  xd->above_context = cm->above_context + mb_col;
  xd->left_context = cm->left_context + mb_row % 4;

  // Distance of Mb to the various image edges. These are specified to 8th pel
  // as they are always compared to values that are in 1/8th pel units
  set_mb_row(cm, xd, mb_row, bh);
  set_mb_col(cm, xd, mb_col, bw);

  xd->dst.y_buffer = dst_fb->y_buffer + recon_yoffset;
  xd->dst.u_buffer = dst_fb->u_buffer + recon_uvoffset;
  xd->dst.v_buffer = dst_fb->v_buffer + recon_uvoffset;
}

static void set_refs(VP9D_COMP *pbi, int mb_row, int mb_col) {
  VP9_COMMON *const cm = &pbi->common;
  MACROBLOCKD *const xd = &pbi->mb;
  MB_MODE_INFO *const mbmi = &xd->mode_info_context->mbmi;

  if (mbmi->ref_frame > INTRA_FRAME) {
    // Select the appropriate reference frame for this MB
    const int fb_idx = cm->active_ref_idx[mbmi->ref_frame - 1];
    const YV12_BUFFER_CONFIG *cfg = &cm->yv12_fb[fb_idx];
    xd->scale_factor[0]    = cm->active_ref_scale[mbmi->ref_frame - 1];
    xd->scale_factor_uv[0] = cm->active_ref_scale[mbmi->ref_frame - 1];
    setup_pred_block(&xd->pre, cfg, mb_row, mb_col,
                     &xd->scale_factor[0], &xd->scale_factor_uv[0]);
    xd->corrupted |= cfg->corrupted;

    if (mbmi->second_ref_frame > INTRA_FRAME) {
      // Select the appropriate reference frame for this MB
      const int second_fb_idx = cm->active_ref_idx[mbmi->second_ref_frame - 1];
      const YV12_BUFFER_CONFIG *second_cfg = &cm->yv12_fb[second_fb_idx];
      xd->scale_factor[1]    = cm->active_ref_scale[mbmi->second_ref_frame - 1];
      xd->scale_factor_uv[1] = cm->active_ref_scale[mbmi->second_ref_frame - 1];
      setup_pred_block(&xd->second_pre, second_cfg, mb_row, mb_col,
                       &xd->scale_factor[1], &xd->scale_factor_uv[1]);
      xd->corrupted |= second_cfg->corrupted;
    }
  }
}

/* Decode a row of Superblocks (2x2 region of MBs) */
static void decode_sb_row(VP9D_COMP *pbi, int mb_row, vp9_reader* r) {
  VP9_COMMON *const pc = &pbi->common;
  MACROBLOCKD *const xd = &pbi->mb;
  int mb_col;

  // For a SB there are 2 left contexts, each pertaining to a MB row within
  vpx_memset(pc->left_context, 0, sizeof(pc->left_context));

  for (mb_col = pc->cur_tile_mb_col_start;
       mb_col < pc->cur_tile_mb_col_end; mb_col += 4) {
    if (vp9_read(r, pc->prob_sb64_coded)) {
      // SB64 decoding
      set_offsets(pbi, BLOCK_SIZE_SB64X64, mb_row, mb_col);
      vp9_decode_mb_mode_mv(pbi, xd, mb_row, mb_col, r);
      set_refs(pbi, mb_row, mb_col);
      decode_sb(pbi, xd, mb_row, mb_col, r, BLOCK_SIZE_SB64X64);
      xd->corrupted |= bool_error(r);
    } else {
      // not SB64
      int j;
      for (j = 0; j < 4; j++) {
        const int x_idx_sb = mb_col + 2 * (j % 2);
        const int y_idx_sb = mb_row + 2 * (j / 2);

        if (y_idx_sb >= pc->mb_rows || x_idx_sb >= pc->mb_cols)
          continue;  // MB lies outside frame, skip on to next

        xd->sb_index = j;

        if (vp9_read(r, pc->prob_sb32_coded)) {
          // SB32 decoding
          set_offsets(pbi, BLOCK_SIZE_SB32X32, y_idx_sb, x_idx_sb);
          vp9_decode_mb_mode_mv(pbi, xd, y_idx_sb, x_idx_sb, r);
          set_refs(pbi, y_idx_sb, x_idx_sb);
          decode_sb(pbi, xd, y_idx_sb, x_idx_sb, r, BLOCK_SIZE_SB32X32);
          xd->corrupted |= bool_error(r);
        } else {
          // not SB32
          // Process the 4 MBs within the SB in the order:
          // top-left, top-right, bottom-left, bottom-right
          int i;
          for (i = 0; i < 4; i++) {
            const int x_idx_mb = x_idx_sb + (i % 2);
            const int y_idx_mb = y_idx_sb + (i / 2);

            if (y_idx_mb >= pc->mb_rows || x_idx_mb >= pc->mb_cols)
              continue;  // MB lies outside frame, skip on to next

            xd->mb_index = i;

            // MB decoding
            set_offsets(pbi, BLOCK_SIZE_MB16X16, y_idx_mb, x_idx_mb);
            vp9_decode_mb_mode_mv(pbi, xd, y_idx_mb, x_idx_mb, r);
            set_refs(pbi, y_idx_mb, x_idx_mb);
            decode_mb(pbi, xd, y_idx_mb, x_idx_mb, r);
            xd->corrupted |= bool_error(r);
          }
        }
      }
    }
  }
}


static void setup_token_decoder(VP9D_COMP *pbi,
                                const uint8_t *data,
                                vp9_reader *r) {
  VP9_COMMON *pc = &pbi->common;
  const uint8_t *data_end = pbi->source + pbi->source_sz;
  const size_t partition_size = data_end - data;

  // Validate the calculated partition length. If the buffer
  // described by the partition can't be fully read, then restrict
  // it to the portion that can be (for EC mode) or throw an error.
  if (!read_is_valid(data, partition_size, data_end))
    vpx_internal_error(&pc->error, VPX_CODEC_CORRUPT_FRAME,
                       "Truncated packet or corrupt partition "
                       "%d length", 1);

  if (vp9_start_decode(r, data, partition_size))
    vpx_internal_error(&pc->error, VPX_CODEC_MEM_ERROR,
                       "Failed to allocate bool decoder %d", 1);
}

static void init_frame(VP9D_COMP *pbi) {
  VP9_COMMON *const pc = &pbi->common;
  MACROBLOCKD *const xd = &pbi->mb;

  if (pc->frame_type == KEY_FRAME) {
    vp9_setup_past_independence(pc, xd);
    // All buffers are implicitly updated on key frames.
    pbi->refresh_frame_flags = (1 << NUM_REF_FRAMES) - 1;
  } else if (pc->error_resilient_mode) {
    vp9_setup_past_independence(pc, xd);
  }

  xd->mode_info_context = pc->mi;
  xd->prev_mode_info_context = pc->prev_mi;
  xd->frame_type = pc->frame_type;
  xd->mode_info_context->mbmi.mode = DC_PRED;
  xd->mode_info_stride = pc->mode_info_stride;
  xd->corrupted = 0;
}

#if CONFIG_CODE_NONZEROCOUNT
static void read_nzc_probs_common(VP9_COMMON *cm,
                                  BOOL_DECODER* const bc,
                                  TX_SIZE tx_size) {
  int c, r, b, t;
  int tokens, nodes;
  vp9_prob *nzc_probs;
  vp9_prob upd;

  if (!get_nzc_used(tx_size)) return;
  if (!vp9_read_bit(bc)) return;

  if (tx_size == TX_32X32) {
    tokens = NZC32X32_TOKENS;
    nzc_probs = cm->fc.nzc_probs_32x32[0][0][0];
    upd = NZC_UPDATE_PROB_32X32;
  } else if (tx_size == TX_16X16) {
    tokens = NZC16X16_TOKENS;
    nzc_probs = cm->fc.nzc_probs_16x16[0][0][0];
    upd = NZC_UPDATE_PROB_16X16;
  } else if (tx_size == TX_8X8) {
    tokens = NZC8X8_TOKENS;
    nzc_probs = cm->fc.nzc_probs_8x8[0][0][0];
    upd = NZC_UPDATE_PROB_8X8;
  } else {
    tokens = NZC4X4_TOKENS;
    nzc_probs = cm->fc.nzc_probs_4x4[0][0][0];
    upd = NZC_UPDATE_PROB_4X4;
  }
  nodes = tokens - 1;
  for (c = 0; c < MAX_NZC_CONTEXTS; ++c) {
    for (r = 0; r < REF_TYPES; ++r) {
      for (b = 0; b < BLOCK_TYPES; ++b) {
        int offset = c * REF_TYPES * BLOCK_TYPES + r * BLOCK_TYPES + b;
        int offset_nodes = offset * nodes;
        for (t = 0; t < nodes; ++t) {
          vp9_prob *p = &nzc_probs[offset_nodes + t];
          if (vp9_read(bc, upd)) {
            *p = read_prob_diff_update(bc, *p);
          }
        }
      }
    }
  }
}

static void read_nzc_pcat_probs(VP9_COMMON *cm, BOOL_DECODER* const bc) {
  int c, t, b;
  vp9_prob upd = NZC_UPDATE_PROB_PCAT;
  if (!(get_nzc_used(TX_4X4) || get_nzc_used(TX_8X8) ||
        get_nzc_used(TX_16X16) || get_nzc_used(TX_32X32)))
    return;
  if (!vp9_read_bit(bc)) {
    return;
  }
  for (c = 0; c < MAX_NZC_CONTEXTS; ++c) {
    for (t = 0; t < NZC_TOKENS_EXTRA; ++t) {
      int bits = vp9_extranzcbits[t + NZC_TOKENS_NOEXTRA];
      for (b = 0; b < bits; ++b) {
        vp9_prob *p = &cm->fc.nzc_pcat_probs[c][t][b];
        if (vp9_read(bc, upd)) {
          *p = read_prob_diff_update(bc, *p);
        }
      }
    }
  }
}

static void read_nzc_probs(VP9_COMMON *cm,
                           BOOL_DECODER* const bc) {
  read_nzc_probs_common(cm, bc, TX_4X4);
  if (cm->txfm_mode != ONLY_4X4)
    read_nzc_probs_common(cm, bc, TX_8X8);
  if (cm->txfm_mode > ALLOW_8X8)
    read_nzc_probs_common(cm, bc, TX_16X16);
  if (cm->txfm_mode > ALLOW_16X16)
    read_nzc_probs_common(cm, bc, TX_32X32);
#ifdef NZC_PCAT_UPDATE
  read_nzc_pcat_probs(cm, bc);
#endif
}
#endif  // CONFIG_CODE_NONZEROCOUNT

static void read_coef_probs_common(VP9D_COMP *pbi,
                                   BOOL_DECODER* const bc,
                                   vp9_coeff_probs *coef_probs,
                                   TX_SIZE tx_size) {
#if CONFIG_MODELCOEFPROB && MODEL_BASED_UPDATE
  const int entropy_nodes_update = UNCONSTRAINED_UPDATE_NODES;
#else
  const int entropy_nodes_update = ENTROPY_NODES;
#endif

  int i, j, k, l, m;

  if (vp9_read_bit(bc)) {
    for (i = 0; i < BLOCK_TYPES; i++) {
      for (j = 0; j < REF_TYPES; j++) {
        for (k = 0; k < COEF_BANDS; k++) {
          for (l = 0; l < PREV_COEF_CONTEXTS; l++) {
#if CONFIG_CODE_NONZEROCOUNT
            const int mstart = get_nzc_used(tx_size);
#else
            const int mstart = 0;
#endif
            if (l >= 3 && k == 0)
              continue;

            for (m = mstart; m < entropy_nodes_update; m++) {
              vp9_prob *const p = coef_probs[i][j][k][l] + m;

              if (vp9_read(bc, vp9_coef_update_prob[m])) {
                *p = read_prob_diff_update(bc, *p);
#if CONFIG_MODELCOEFPROB && MODEL_BASED_UPDATE
                if (m == UNCONSTRAINED_NODES - 1)
                  vp9_get_model_distribution(*p, coef_probs[i][j][k][l], i, j);
#endif
              }
            }
          }
        }
      }
    }
  }
}

static void read_coef_probs(VP9D_COMP *pbi, BOOL_DECODER* const bc) {
  VP9_COMMON *const pc = &pbi->common;

  read_coef_probs_common(pbi, bc, pc->fc.coef_probs_4x4, TX_4X4);

  if (pbi->common.txfm_mode != ONLY_4X4)
    read_coef_probs_common(pbi, bc, pc->fc.coef_probs_8x8, TX_8X8);

  if (pbi->common.txfm_mode > ALLOW_8X8)
    read_coef_probs_common(pbi, bc, pc->fc.coef_probs_16x16, TX_16X16);

  if (pbi->common.txfm_mode > ALLOW_16X16)
    read_coef_probs_common(pbi, bc, pc->fc.coef_probs_32x32, TX_32X32);
}

static void update_frame_size(VP9D_COMP *pbi) {
  VP9_COMMON *cm = &pbi->common;

  const int width = multiple16(cm->width);
  const int height = multiple16(cm->height);

  cm->mb_rows = height / 16;
  cm->mb_cols = width / 16;
  cm->MBs = cm->mb_rows * cm->mb_cols;
  cm->mode_info_stride = cm->mb_cols + 1;
  memset(cm->mip, 0,
        (cm->mb_cols + 1) * (cm->mb_rows + 1) * sizeof(MODE_INFO));
  vp9_update_mode_info_border(cm, cm->mip);

  cm->mi = cm->mip + cm->mode_info_stride + 1;
  cm->prev_mi = cm->prev_mip + cm->mode_info_stride + 1;
  vp9_update_mode_info_in_image(cm, cm->mi);
}

static void setup_segmentation(VP9_COMMON *pc, MACROBLOCKD *xd, vp9_reader *r) {
  int i, j;

  xd->segmentation_enabled = vp9_read_bit(r);
  if (xd->segmentation_enabled) {
    // Read whether or not the segmentation map is being explicitly updated
    // this frame.
    xd->update_mb_segmentation_map = vp9_read_bit(r);

    if (xd->update_mb_segmentation_map) {
      // Which macro block level features are enabled. Read the probs used to
      // decode the segment id for each macro block.
      for (i = 0; i < MB_FEATURE_TREE_PROBS; i++)
        xd->mb_segment_tree_probs[i] = vp9_read_bit(r) ? vp9_read_prob(r) : 255;

      // Read the prediction probs needed to decode the segment id
      pc->temporal_update = vp9_read_bit(r);
      if (pc->temporal_update) {
        const vp9_prob *p = xd->mb_segment_tree_probs;
        vp9_prob *p_mod = xd->mb_segment_mispred_tree_probs;

        const int c0 =        p[0]  *        p[1];
        const int c1 =        p[0]  * (256 - p[1]);
        const int c2 = (256 - p[0]) *        p[2];
        const int c3 = (256 - p[0]) * (256 - p[2]);

        p_mod[0] = get_binary_prob(c1, c2 + c3);
        p_mod[1] = get_binary_prob(c0, c2 + c3);
        p_mod[2] = get_binary_prob(c0 + c1, c3);
        p_mod[3] = get_binary_prob(c0 + c1, c2);

        for (i = 0; i < PREDICTION_PROBS; i++)
          pc->segment_pred_probs[i] = vp9_read_bit(r) ? vp9_read_prob(r) : 255;
      } else {
        for (i = 0; i < PREDICTION_PROBS; i++)
          pc->segment_pred_probs[i] = 255;
      }
    }

    xd->update_mb_segmentation_data = vp9_read_bit(r);
    if (xd->update_mb_segmentation_data) {
      xd->mb_segment_abs_delta = vp9_read_bit(r);

      vp9_clearall_segfeatures(xd);

      for (i = 0; i < MAX_MB_SEGMENTS; i++) {
        for (j = 0; j < SEG_LVL_MAX; j++) {
          int data = 0;
          const int feature_enabled = vp9_read_bit(r);
          if (feature_enabled) {
            vp9_enable_segfeature(xd, i, j);
            data = vp9_decode_unsigned_max(r, vp9_seg_feature_data_max(j));
            if (vp9_is_segfeature_signed(j))
              data = vp9_read_and_apply_sign(r, data);
          }
          vp9_set_segdata(xd, i, j, data);
        }
      }
    }
  }
}

static void setup_pred_probs(VP9_COMMON *pc, vp9_reader *r) {
  // Read common prediction model status flag probability updates for the
  // reference frame
  if (pc->frame_type == KEY_FRAME) {
    // Set the prediction probabilities to defaults
    pc->ref_pred_probs[0] = 120;
    pc->ref_pred_probs[1] = 80;
    pc->ref_pred_probs[2] = 40;
  } else {
    int i;
    for (i = 0; i < PREDICTION_PROBS; ++i)
      if (vp9_read_bit(r))
        pc->ref_pred_probs[i] = vp9_read_prob(r);
  }
}

static void setup_loopfilter(VP9_COMMON *pc, MACROBLOCKD *xd, vp9_reader *r) {
  pc->filter_type = (LOOPFILTER_TYPE) vp9_read_bit(r);
  pc->filter_level = vp9_read_literal(r, 6);
  pc->sharpness_level = vp9_read_literal(r, 3);

#if CONFIG_LOOP_DERING
  if (vp9_read_bit(r))
    pc->dering_enabled = 1 + vp9_read_literal(r, 4);
  else
    pc->dering_enabled = 0;
#endif

  // Read in loop filter deltas applied at the MB level based on mode or ref
  // frame.
  xd->mode_ref_lf_delta_update = 0;
  xd->mode_ref_lf_delta_enabled = vp9_read_bit(r);

  if (xd->mode_ref_lf_delta_enabled) {
    xd->mode_ref_lf_delta_update = vp9_read_bit(r);
    if (xd->mode_ref_lf_delta_update) {
      int i;

      for (i = 0; i < MAX_REF_LF_DELTAS; i++) {
        if (vp9_read_bit(r)) {
          const int value = vp9_read_literal(r, 6);
          xd->ref_lf_deltas[i] = vp9_read_and_apply_sign(r, value);
        }
      }

      for (i = 0; i < MAX_MODE_LF_DELTAS; i++) {
        if (vp9_read_bit(r)) {
          const int value = vp9_read_literal(r, 6);
          xd->mode_lf_deltas[i] = vp9_read_and_apply_sign(r, value);
        }
      }
    }
  }
}

static void setup_quantization(VP9D_COMP *pbi, vp9_reader *r) {
  // Read the default quantizers
  VP9_COMMON *const pc = &pbi->common;

  pc->base_qindex = vp9_read_literal(r, QINDEX_BITS);
  if (get_delta_q(r, &pc->y_dc_delta_q) |
      get_delta_q(r, &pc->uv_dc_delta_q) |
      get_delta_q(r, &pc->uv_ac_delta_q))
    vp9_init_de_quantizer(pbi);

  mb_init_dequantizer(pbi, &pbi->mb);  // MB level dequantizer setup
}

static const uint8_t *read_frame_size(VP9_COMMON *const pc, const uint8_t *data,
                                      const uint8_t *data_end,
                                      int *width, int *height) {
  if (data + 4 < data_end) {
    *width = read_le16(data);
    *height = read_le16(data + 2);
    data += 4;
  } else {
    vpx_internal_error(&pc->error, VPX_CODEC_CORRUPT_FRAME,
                       "Failed to read frame size");
  }
  return data;
}

static const uint8_t *setup_frame_size(VP9D_COMP *pbi, int scaling_active,
                                      const uint8_t *data,
                                      const uint8_t *data_end) {
  // If error concealment is enabled we should only parse the new size
  // if we have enough data. Otherwise we will end up with the wrong size.
  VP9_COMMON *const pc = &pbi->common;
  int display_width = pc->display_width;
  int display_height = pc->display_height;
  int width = pc->width;
  int height = pc->height;

  if (scaling_active)
    data = read_frame_size(pc, data, data_end, &display_width, &display_height);

  data = read_frame_size(pc, data, data_end, &width, &height);

  if (pc->width != width || pc->height != height) {
    if (width <= 0)
      vpx_internal_error(&pc->error, VPX_CODEC_CORRUPT_FRAME,
                         "Invalid frame width");

    if (height <= 0)
      vpx_internal_error(&pc->error, VPX_CODEC_CORRUPT_FRAME,
                         "Invalid frame height");

    if (!pbi->initial_width || !pbi->initial_height) {
      if (vp9_alloc_frame_buffers(pc, width, height))
        vpx_internal_error(&pc->error, VPX_CODEC_MEM_ERROR,
                           "Failed to allocate frame buffers");
        pbi->initial_width = width;
        pbi->initial_height = height;
    } else {
      if (width > pbi->initial_width)
        vpx_internal_error(&pc->error, VPX_CODEC_CORRUPT_FRAME,
                           "Frame width too large");

      if (height > pbi->initial_height)
        vpx_internal_error(&pc->error, VPX_CODEC_CORRUPT_FRAME,
                           "Frame height too large");
    }

    pc->width = width;
    pc->height = height;
    pc->display_width = scaling_active ? display_width : width;
    pc->display_height = scaling_active ? display_height : height;

    update_frame_size(pbi);
  }

  return data;
}

static void update_frame_context(VP9D_COMP *pbi, vp9_reader *r) {
  FRAME_CONTEXT *const fc = &pbi->common.fc;

  vp9_copy(fc->pre_coef_probs_4x4, fc->coef_probs_4x4);
  vp9_copy(fc->pre_coef_probs_8x8, fc->coef_probs_8x8);
  vp9_copy(fc->pre_coef_probs_16x16, fc->coef_probs_16x16);
  vp9_copy(fc->pre_coef_probs_32x32, fc->coef_probs_32x32);
  vp9_copy(fc->pre_ymode_prob, fc->ymode_prob);
  vp9_copy(fc->pre_sb_ymode_prob, fc->sb_ymode_prob);
  vp9_copy(fc->pre_uv_mode_prob, fc->uv_mode_prob);
  vp9_copy(fc->pre_bmode_prob, fc->bmode_prob);
  vp9_copy(fc->pre_i8x8_mode_prob, fc->i8x8_mode_prob);
  vp9_copy(fc->pre_sub_mv_ref_prob, fc->sub_mv_ref_prob);
  vp9_copy(fc->pre_mbsplit_prob, fc->mbsplit_prob);
  fc->pre_nmvc = fc->nmvc;

  vp9_zero(fc->coef_counts_4x4);
  vp9_zero(fc->coef_counts_8x8);
  vp9_zero(fc->coef_counts_16x16);
  vp9_zero(fc->coef_counts_32x32);
  vp9_zero(fc->eob_branch_counts);
  vp9_zero(fc->ymode_counts);
  vp9_zero(fc->sb_ymode_counts);
  vp9_zero(fc->uv_mode_counts);
  vp9_zero(fc->bmode_counts);
  vp9_zero(fc->i8x8_mode_counts);
  vp9_zero(fc->sub_mv_ref_counts);
  vp9_zero(fc->mbsplit_counts);
  vp9_zero(fc->NMVcount);
  vp9_zero(fc->mv_ref_ct);

#if CONFIG_COMP_INTERINTRA_PRED
  fc->pre_interintra_prob = fc->interintra_prob;
  vp9_zero(fc->interintra_counts);
#endif

#if CONFIG_CODE_NONZEROCOUNT
  vp9_copy(fc->pre_nzc_probs_4x4, fc->nzc_probs_4x4);
  vp9_copy(fc->pre_nzc_probs_8x8, fc->nzc_probs_8x8);
  vp9_copy(fc->pre_nzc_probs_16x16, fc->nzc_probs_16x16);
  vp9_copy(fc->pre_nzc_probs_32x32, fc->nzc_probs_32x32);
  vp9_copy(fc->pre_nzc_pcat_probs, fc->nzc_pcat_probs);

  vp9_zero(fc->nzc_counts_4x4);
  vp9_zero(fc->nzc_counts_8x8);
  vp9_zero(fc->nzc_counts_16x16);
  vp9_zero(fc->nzc_counts_32x32);
  vp9_zero(fc->nzc_pcat_counts);
#endif

  read_coef_probs(pbi, r);
#if CONFIG_CODE_NONZEROCOUNT
  read_nzc_probs(&pbi->common, r);
#endif
}

static void decode_tiles(VP9D_COMP *pbi,
                         const uint8_t *data, int first_partition_size,
                         BOOL_DECODER *header_bc, BOOL_DECODER *residual_bc) {
  VP9_COMMON *const pc = &pbi->common;

  const uint8_t *data_ptr = data + first_partition_size;
  int tile_row, tile_col, delta_log2_tiles;
  int mb_row;

  vp9_get_tile_n_bits(pc, &pc->log2_tile_columns, &delta_log2_tiles);
  while (delta_log2_tiles--) {
    if (vp9_read_bit(header_bc)) {
      pc->log2_tile_columns++;
    } else {
      break;
    }
  }
  pc->log2_tile_rows = vp9_read_bit(header_bc);
  if (pc->log2_tile_rows)
    pc->log2_tile_rows += vp9_read_bit(header_bc);
  pc->tile_columns = 1 << pc->log2_tile_columns;
  pc->tile_rows    = 1 << pc->log2_tile_rows;

  vpx_memset(pc->above_context, 0,
             sizeof(ENTROPY_CONTEXT_PLANES) * pc->mb_cols);

  if (pbi->oxcf.inv_tile_order) {
    const int n_cols = pc->tile_columns;
    const uint8_t *data_ptr2[4][1 << 6];
    BOOL_DECODER UNINITIALIZED_IS_SAFE(bc_bak);

    // pre-initialize the offsets, we're going to read in inverse order
    data_ptr2[0][0] = data_ptr;
    for (tile_row = 0; tile_row < pc->tile_rows; tile_row++) {
      if (tile_row) {
        const int size = read_le32(data_ptr2[tile_row - 1][n_cols - 1]);
        data_ptr2[tile_row - 1][n_cols - 1] += 4;
        data_ptr2[tile_row][0] = data_ptr2[tile_row - 1][n_cols - 1] + size;
      }

      for (tile_col = 1; tile_col < n_cols; tile_col++) {
        const int size = read_le32(data_ptr2[tile_row][tile_col - 1]);
        data_ptr2[tile_row][tile_col - 1] += 4;
        data_ptr2[tile_row][tile_col] =
            data_ptr2[tile_row][tile_col - 1] + size;
      }
    }

    for (tile_row = 0; tile_row < pc->tile_rows; tile_row++) {
      vp9_get_tile_row_offsets(pc, tile_row);
      for (tile_col = n_cols - 1; tile_col >= 0; tile_col--) {
        vp9_get_tile_col_offsets(pc, tile_col);
        setup_token_decoder(pbi, data_ptr2[tile_row][tile_col], residual_bc);

        // Decode a row of superblocks
        for (mb_row = pc->cur_tile_mb_row_start;
             mb_row < pc->cur_tile_mb_row_end; mb_row += 4) {
          decode_sb_row(pbi, mb_row, residual_bc);
        }

        if (tile_row == pc->tile_rows - 1 && tile_col == n_cols - 1)
          bc_bak = *residual_bc;
      }
    }
    *residual_bc = bc_bak;
  } else {
    int has_more;

    for (tile_row = 0; tile_row < pc->tile_rows; tile_row++) {
      vp9_get_tile_row_offsets(pc, tile_row);
      for (tile_col = 0; tile_col < pc->tile_columns; tile_col++) {
        vp9_get_tile_col_offsets(pc, tile_col);

        has_more = tile_col < pc->tile_columns - 1 ||
                   tile_row < pc->tile_rows - 1;

        // Setup decoder
        setup_token_decoder(pbi, data_ptr + (has_more ? 4 : 0), residual_bc);

        // Decode a row of superblocks
        for (mb_row = pc->cur_tile_mb_row_start;
             mb_row < pc->cur_tile_mb_row_end; mb_row += 4) {
          decode_sb_row(pbi, mb_row, residual_bc);
        }

        if (has_more) {
          const int size = read_le32(data_ptr);
          data_ptr += 4 + size;
        }
      }
    }
  }
}

int vp9_decode_frame(VP9D_COMP *pbi, const uint8_t **p_data_end) {
  BOOL_DECODER header_bc, residual_bc;
  VP9_COMMON *const pc = &pbi->common;
  MACROBLOCKD *const xd  = &pbi->mb;
  const uint8_t *data = pbi->source;
  const uint8_t *data_end = data + pbi->source_sz;
  size_t first_partition_size = 0;
  int i, corrupt_tokens = 0;

  // printf("Decoding frame %d\n", pc->current_video_frame);

  xd->corrupted = 0;  // start with no corruption of current frame
  pc->yv12_fb[pc->new_fb_idx].corrupted = 0;

  if (data_end - data < 3) {
    vpx_internal_error(&pc->error, VPX_CODEC_CORRUPT_FRAME, "Truncated packet");
  } else {
    int scaling_active;
    pc->last_frame_type = pc->frame_type;
    pc->frame_type = (FRAME_TYPE)(data[0] & 1);
    pc->version = (data[0] >> 1) & 7;
    pc->show_frame = (data[0] >> 4) & 1;
    scaling_active = (data[0] >> 5) & 1;
    first_partition_size = read_le16(data + 1);

    if (!read_is_valid(data, first_partition_size, data_end))
      vpx_internal_error(&pc->error, VPX_CODEC_CORRUPT_FRAME,
                         "Truncated packet or corrupt partition 0 length");

    data += 3;

    vp9_setup_version(pc);

    if (pc->frame_type == KEY_FRAME) {
      // When error concealment is enabled we should only check the sync
      // code if we have enough bits available
      if (data + 3 < data_end) {
        if (data[0] != 0x9d || data[1] != 0x01 || data[2] != 0x2a)
          vpx_internal_error(&pc->error, VPX_CODEC_UNSUP_BITSTREAM,
                             "Invalid frame sync code");
      }
      data += 3;
    }

    data = setup_frame_size(pbi, scaling_active, data, data_end);
  }

  if ((!pbi->decoded_key_frame && pc->frame_type != KEY_FRAME) ||
      pc->width == 0 || pc->height == 0) {
    return -1;
  }

  init_frame(pbi);

  // Reset the frame pointers to the current frame size
  vp8_yv12_realloc_frame_buffer(&pc->yv12_fb[pc->new_fb_idx],
                                pc->width, pc->height,
                                VP9BORDERINPIXELS);

  if (vp9_start_decode(&header_bc, data, first_partition_size))
    vpx_internal_error(&pc->error, VPX_CODEC_MEM_ERROR,
                       "Failed to allocate bool decoder 0");

  pc->clr_type = (YUV_TYPE)vp9_read_bit(&header_bc);
  pc->clamp_type = (CLAMP_TYPE)vp9_read_bit(&header_bc);
  pc->error_resilient_mode = vp9_read_bit(&header_bc);

  setup_segmentation(pc, xd, &header_bc);

  setup_pred_probs(pc, &header_bc);

  pc->prob_sb64_coded = vp9_read_prob(&header_bc);
  pc->prob_sb32_coded = vp9_read_prob(&header_bc);
  xd->lossless = vp9_read_bit(&header_bc);
  pc->txfm_mode = xd->lossless ? ONLY_4X4 : read_txfm_mode(&header_bc);
  if (pc->txfm_mode == TX_MODE_SELECT) {
    pc->prob_tx[0] = vp9_read_prob(&header_bc);
    pc->prob_tx[1] = vp9_read_prob(&header_bc);
    pc->prob_tx[2] = vp9_read_prob(&header_bc);
  }

  setup_loopfilter(pc, xd, &header_bc);

  // Dummy read for now
  vp9_read_literal(&header_bc, 2);

  setup_quantization(pbi, &header_bc);

  // Determine if the golden frame or ARF buffer should be updated and how.
  // For all non key frames the GF and ARF refresh flags and sign bias
  // flags must be set explicitly.
  if (pc->frame_type == KEY_FRAME) {
    for (i = 0; i < ALLOWED_REFS_PER_FRAME; ++i)
      pc->active_ref_idx[i] = pc->new_fb_idx;
  } else {
    // Should the GF or ARF be updated from the current frame
    pbi->refresh_frame_flags = vp9_read_literal(&header_bc, NUM_REF_FRAMES);

    // Select active reference frames
    for (i = 0; i < ALLOWED_REFS_PER_FRAME; ++i) {
      int ref_frame_num = vp9_read_literal(&header_bc, NUM_REF_FRAMES_LG2);
      pc->active_ref_idx[i] = pc->ref_frame_map[ref_frame_num];
    }

    pc->ref_frame_sign_bias[GOLDEN_FRAME] = vp9_read_bit(&header_bc);
    pc->ref_frame_sign_bias[ALTREF_FRAME] = vp9_read_bit(&header_bc);

    // Is high precision mv allowed
    xd->allow_high_precision_mv = vp9_read_bit(&header_bc);

    // Read the type of subpel filter to use
    pc->mcomp_filter_type = vp9_read_bit(&header_bc)
                                ? SWITCHABLE
                                : vp9_read_literal(&header_bc, 2);

#if CONFIG_COMP_INTERINTRA_PRED
    pc->use_interintra = vp9_read_bit(&header_bc);
#endif

    // Calculate scaling factors for each of the 3 available references
    for (i = 0; i < ALLOWED_REFS_PER_FRAME; ++i) {
      const int idx = pc->active_ref_idx[i];
      struct scale_factors *sf = &pc->active_ref_scale[i];
      if (idx >= NUM_YV12_BUFFERS)
        memset(sf, 0, sizeof(*sf));
      else
        vp9_setup_scale_factors_for_frame(sf, &pc->yv12_fb[idx],
                                          pc->width, pc->height);
    }

    // To enable choice of different interpolation filters
    vp9_setup_interp_filters(xd, pc->mcomp_filter_type, pc);
  }

  if (!pc->error_resilient_mode) {
    pc->refresh_entropy_probs = vp9_read_bit(&header_bc);
    pc->frame_parallel_decoding_mode = vp9_read_bit(&header_bc);
  } else {
    pc->refresh_entropy_probs = 0;
    pc->frame_parallel_decoding_mode = 1;
  }
  pc->frame_context_idx = vp9_read_literal(&header_bc, NUM_FRAME_CONTEXTS_LG2);
  vpx_memcpy(&pc->fc, &pc->frame_contexts[pc->frame_context_idx],
             sizeof(pc->fc));

  // Read inter mode probability context updates
  if (pc->frame_type != KEY_FRAME) {
    int i, j;
    for (i = 0; i < INTER_MODE_CONTEXTS; ++i)
      for (j = 0; j < 4; ++j)
        if (vp9_read(&header_bc, 252))
          pc->fc.vp9_mode_contexts[i][j] = vp9_read_prob(&header_bc);
  }
#if CONFIG_MODELCOEFPROB
  if (pc->frame_type == KEY_FRAME) {
    vp9_default_coef_probs(pc);
  }
#endif

#if CONFIG_NEW_MVREF
  // If Key frame reset mv ref id probabilities to defaults
  if (pc->frame_type != KEY_FRAME) {
    // Read any mv_ref index probability updates
    int i, j;

    for (i = 0; i < MAX_REF_FRAMES; ++i) {
      // Skip the dummy entry for intra ref frame.
      if (i == INTRA_FRAME)
        continue;

      // Read any updates to probabilities
      for (j = 0; j < MAX_MV_REF_CANDIDATES - 1; ++j)
        if (vp9_read(&header_bc, VP9_MVREF_UPDATE_PROB))
          xd->mb_mv_ref_probs[i][j] = vp9_read_prob(&header_bc);
    }
  }
#endif

  if (0) {
    FILE *z = fopen("decodestats.stt", "a");
    fprintf(z, "%6d F:%d,R:%d,Q:%d\n",
            pc->current_video_frame,
            pc->frame_type,
            pbi->refresh_frame_flags,
            pc->base_qindex);
    fclose(z);
  }

  update_frame_context(pbi, &header_bc);

  // Initialize xd pointers. Any reference should do for xd->pre, so use 0.
  vpx_memcpy(&xd->pre, &pc->yv12_fb[pc->active_ref_idx[0]],
             sizeof(YV12_BUFFER_CONFIG));
  vpx_memcpy(&xd->dst, &pc->yv12_fb[pc->new_fb_idx],
             sizeof(YV12_BUFFER_CONFIG));

  // Create the segmentation map structure and set to 0
  if (!pc->last_frame_seg_map)
    CHECK_MEM_ERROR(pc->last_frame_seg_map,
                    vpx_calloc((pc->mb_rows * pc->mb_cols), 1));

  // set up frame new frame for intra coded blocks
  vp9_setup_intra_recon(&pc->yv12_fb[pc->new_fb_idx]);

  vp9_setup_block_dptrs(xd);
  vp9_build_block_doffsets(xd);

  // clear out the coeff buffer
  vpx_memset(xd->plane[0].qcoeff, 0, sizeof(xd->plane[0].qcoeff));
  vpx_memset(xd->plane[1].qcoeff, 0, sizeof(xd->plane[1].qcoeff));
  vpx_memset(xd->plane[2].qcoeff, 0, sizeof(xd->plane[2].qcoeff));

  vp9_read_bit(&header_bc);  // unused

  vp9_decode_mode_mvs_init(pbi, &header_bc);

  decode_tiles(pbi, data, first_partition_size, &header_bc, &residual_bc);
  corrupt_tokens |= xd->corrupted;

  // keep track of the last coded dimensions
  pc->last_width = pc->width;
  pc->last_height = pc->height;

  // Collect information about decoder corruption.
  // 1. Check first boolean decoder for errors.
  // 2. Check the macroblock information
  pc->yv12_fb[pc->new_fb_idx].corrupted = bool_error(&header_bc) |
                                          corrupt_tokens;

  if (!pbi->decoded_key_frame) {
    if (pc->frame_type == KEY_FRAME && !pc->yv12_fb[pc->new_fb_idx].corrupted)
      pbi->decoded_key_frame = 1;
    else
      vpx_internal_error(&pbi->common.error, VPX_CODEC_CORRUPT_FRAME,
                         "A stream must start with a complete key frame");
  }

  if (!pc->error_resilient_mode && !pc->frame_parallel_decoding_mode) {
    vp9_adapt_coef_probs(pc);
#if CONFIG_CODE_NONZEROCOUNT
    vp9_adapt_nzc_probs(pc);
#endif
  }

  if (pc->frame_type != KEY_FRAME) {
    if (!pc->error_resilient_mode && !pc->frame_parallel_decoding_mode) {
      vp9_adapt_mode_probs(pc);
      vp9_adapt_nmv_probs(pc, xd->allow_high_precision_mv);
      vp9_adapt_mode_context(&pbi->common);
    }
  }

  if (pc->refresh_entropy_probs) {
    vpx_memcpy(&pc->frame_contexts[pc->frame_context_idx], &pc->fc,
               sizeof(pc->fc));
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

  *p_data_end = vp9_reader_find_end(&residual_bc);
  return 0;
}
