/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include "./vpx_config.h"
#include "vp9_rtcd.h"
#include "vp9/encoder/vp9_quantize.h"
#include "vp9/common/vp9_reconintra.h"
#include "vp9/encoder/vp9_encodemb.h"
#include "vp9/common/vp9_invtrans.h"
#include "vp9/encoder/vp9_encodeintra.h"

static void encode_intra4x4block(MACROBLOCK *x, int ib, BLOCK_SIZE_TYPE bs);

int vp9_encode_intra(VP9_COMP *cpi, MACROBLOCK *x, int use_16x16_pred) {
  MB_MODE_INFO * mbmi = &x->e_mbd.mode_info_context->mbmi;
  (void) cpi;

  if (use_16x16_pred) {
    mbmi->mode = DC_PRED;
    mbmi->uv_mode = DC_PRED;
    mbmi->ref_frame = INTRA_FRAME;

    vp9_encode_intra16x16mby(&cpi->common, x);
  } else {
    int i;

    for (i = 0; i < 16; i++) {
      x->e_mbd.mode_info_context->bmi[i].as_mode.first = B_DC_PRED;
      encode_intra4x4block(x, i, BLOCK_SIZE_MB16X16);
    }
  }

  return vp9_get_mb_ss(x->plane[0].src_diff);
}

static void encode_intra4x4block(MACROBLOCK *x, int ib,
                                 BLOCK_SIZE_TYPE bsize) {
  MACROBLOCKD * const xd = &x->e_mbd;
  TX_TYPE tx_type;
  uint8_t* const src =
      raster_block_offset_uint8(xd, bsize, 0, ib,
                                x->plane[0].src.buf, x->plane[0].src.stride);
  uint8_t* const dst =
      raster_block_offset_uint8(xd, bsize, 0, ib,
                                xd->plane[0].dst.buf, xd->plane[0].dst.stride);
  int16_t* const src_diff =
      raster_block_offset_int16(xd, bsize, 0, ib,
                                x->plane[0].src_diff);
  int16_t* const diff =
      raster_block_offset_int16(xd, bsize, 0, ib,
                                xd->plane[0].diff);
  int16_t* const coeff = BLOCK_OFFSET(x->plane[0].coeff, ib, 16);

  assert(ib < 16);

#if CONFIG_NEWBINTRAMODES
  xd->mode_info_context->bmi[ib].as_mode.context =
    vp9_find_bpred_context(&x->e_mbd, ib, dst, xd->plane[0].dst.stride);
#endif

  vp9_intra4x4_predict(&x->e_mbd, ib,
                       xd->mode_info_context->bmi[ib].as_mode.first,
                       dst, xd->plane[0].dst.stride);
  vp9_subtract_block(4, 4, src_diff, 16,
                     src, x->plane[0].src.stride,
                     dst, xd->plane[0].dst.stride);

  tx_type = get_tx_type_4x4(&x->e_mbd, ib);
  if (tx_type != DCT_DCT) {
    vp9_short_fht4x4(src_diff, coeff, 16, tx_type);
    x->quantize_b_4x4(x, ib, tx_type, 16);
    vp9_short_iht4x4(BLOCK_OFFSET(xd->plane[0].dqcoeff, ib, 16),
                     diff, 16, tx_type);
  } else {
    x->fwd_txm4x4(src_diff, coeff, 32);
    x->quantize_b_4x4(x, ib, tx_type, 16);
    vp9_inverse_transform_b_4x4(&x->e_mbd, xd->plane[0].eobs[ib],
                                BLOCK_OFFSET(xd->plane[0].dqcoeff, ib, 16),
                                diff, 32);
  }

  vp9_recon_b(dst, diff, dst, xd->plane[0].dst.stride);
}

void vp9_encode_intra4x4mby(MACROBLOCK *mb, BLOCK_SIZE_TYPE bsize) {
  int i;
  int bwl = b_width_log2(bsize), bhl = b_height_log2(bsize);
  int bc = 1 << (bwl + bhl);

  for (i = 0; i < bc; i++)
    encode_intra4x4block(mb, i, bsize);
}

void vp9_encode_intra16x16mby(VP9_COMMON *const cm, MACROBLOCK *x) {
  MACROBLOCKD *xd = &x->e_mbd;
  TX_SIZE tx_size = xd->mode_info_context->mbmi.txfm_size;

  vp9_build_intra_predictors_sby_s(xd, BLOCK_SIZE_MB16X16);
  vp9_subtract_sby(x, BLOCK_SIZE_MB16X16);

  switch (tx_size) {
    case TX_16X16:
      vp9_transform_sby_16x16(x, BLOCK_SIZE_MB16X16);
      vp9_quantize_sby_16x16(x, BLOCK_SIZE_MB16X16);
      if (x->optimize)
        vp9_optimize_sby_16x16(cm, x, BLOCK_SIZE_MB16X16);
      vp9_inverse_transform_sby_16x16(xd, BLOCK_SIZE_MB16X16);
      break;
    case TX_8X8:
      vp9_transform_sby_8x8(x, BLOCK_SIZE_MB16X16);
      vp9_quantize_sby_8x8(x, BLOCK_SIZE_MB16X16);
      if (x->optimize)
        vp9_optimize_sby_8x8(cm, x, BLOCK_SIZE_MB16X16);
      vp9_inverse_transform_sby_8x8(xd, BLOCK_SIZE_MB16X16);
      break;
    default:
      vp9_transform_sby_4x4(x, BLOCK_SIZE_MB16X16);
      vp9_quantize_sby_4x4(x, BLOCK_SIZE_MB16X16);
      if (x->optimize)
        vp9_optimize_sby_4x4(cm, x, BLOCK_SIZE_MB16X16);
      vp9_inverse_transform_sby_4x4(xd, BLOCK_SIZE_MB16X16);
      break;
  }

  vp9_recon_sby(xd, BLOCK_SIZE_MB16X16);
}

void vp9_encode_intra16x16mbuv(VP9_COMMON *const cm, MACROBLOCK *x) {
  MACROBLOCKD *xd = &x->e_mbd;
  TX_SIZE tx_size = xd->mode_info_context->mbmi.txfm_size;

  vp9_build_intra_predictors_sbuv_s(xd, BLOCK_SIZE_MB16X16);
  vp9_subtract_sbuv(x, BLOCK_SIZE_MB16X16);

  switch (tx_size) {
    case TX_4X4:
      vp9_transform_sbuv_4x4(x, BLOCK_SIZE_MB16X16);
      vp9_quantize_sbuv_4x4(x, BLOCK_SIZE_MB16X16);
      if (x->optimize)
        vp9_optimize_sbuv_4x4(cm, x, BLOCK_SIZE_MB16X16);
      vp9_inverse_transform_sbuv_4x4(xd, BLOCK_SIZE_MB16X16);
      break;
    default:  // 16x16 or 8x8
      vp9_transform_sbuv_8x8(x, BLOCK_SIZE_MB16X16);
      vp9_quantize_sbuv_8x8(x, BLOCK_SIZE_MB16X16);
      if (x->optimize)
        vp9_optimize_sbuv_8x8(cm, x, BLOCK_SIZE_MB16X16);
      vp9_inverse_transform_sbuv_8x8(xd, BLOCK_SIZE_MB16X16);
      break;
    }

  vp9_recon_sbuv(xd, BLOCK_SIZE_MB16X16);
}

void vp9_encode_intra8x8(MACROBLOCK *x, int ib) {
  MACROBLOCKD *xd = &x->e_mbd;
  uint8_t* const src =
      raster_block_offset_uint8(xd, BLOCK_SIZE_MB16X16, 0, ib,
                                x->plane[0].src.buf, x->plane[0].src.stride);
  int16_t* const src_diff =
      raster_block_offset_int16(xd, BLOCK_SIZE_MB16X16, 0, ib,
                                x->plane[0].src_diff);
  int16_t* const diff =
      raster_block_offset_int16(xd, BLOCK_SIZE_MB16X16, 0, ib,
                                xd->plane[0].diff);
  uint8_t* const dst =
      raster_block_offset_uint8(xd, BLOCK_SIZE_MB16X16, 0, ib,
                                xd->plane[0].dst.buf, xd->plane[0].dst.stride);
  const int iblock[4] = {0, 1, 4, 5};
  int i;
  TX_TYPE tx_type;

  vp9_intra8x8_predict(xd, ib, xd->mode_info_context->bmi[ib].as_mode.first,
                       dst, xd->plane[0].dst.stride);
  // generate residual blocks
  vp9_subtract_block(8, 8, src_diff, 16,
                     src, x->plane[0].src.stride,
                     dst, xd->plane[0].dst.stride);

  if (xd->mode_info_context->mbmi.txfm_size == TX_8X8) {
    int idx = (ib & 0x02) ? (ib + 2) : ib;
    int16_t* const dqcoeff = BLOCK_OFFSET(xd->plane[0].dqcoeff, idx, 16);
    int16_t* const coeff = BLOCK_OFFSET(x->plane[0].coeff, idx, 16);

    assert(idx < 16);
    tx_type = get_tx_type_8x8(xd, ib);
    if (tx_type != DCT_DCT) {
      vp9_short_fht8x8(src_diff, coeff, 16, tx_type);
      x->quantize_b_8x8(x, idx, tx_type, 16);
      vp9_short_iht8x8(dqcoeff, diff, 16, tx_type);
    } else {
      x->fwd_txm8x8(src_diff, coeff, 32);
      x->quantize_b_8x8(x, idx, DCT_DCT, 16);
      vp9_short_idct8x8(dqcoeff, diff, 32);
    }
  } else {
    for (i = 0; i < 4; i++) {
      int idx = ib + iblock[i];
      int16_t* const dqcoeff = BLOCK_OFFSET(xd->plane[0].dqcoeff, idx, 16);
      int16_t* const coeff = BLOCK_OFFSET(x->plane[0].coeff, idx, 16);
      int16_t* const src_diff =
          raster_block_offset_int16(xd, BLOCK_SIZE_MB16X16, 0, idx,
                                    x->plane[0].src_diff);
      int16_t* const diff =
          raster_block_offset_int16(xd, BLOCK_SIZE_MB16X16, 0, idx,
                                    xd->plane[0].diff);

      assert(idx < 16);
      tx_type = get_tx_type_4x4(xd, ib + iblock[i]);
      if (tx_type != DCT_DCT) {
        vp9_short_fht4x4(src_diff, coeff, 16, tx_type);
        x->quantize_b_4x4(x, ib + iblock[i], tx_type, 16);
        vp9_short_iht4x4(dqcoeff, diff, 16, tx_type);
      } else if (!(i & 1) &&
                 get_tx_type_4x4(xd, ib + iblock[i] + 1) == DCT_DCT) {
        x->fwd_txm8x4(src_diff, coeff, 32);
        x->quantize_b_4x4_pair(x, ib + iblock[i], ib + iblock[i] + 1, 16);
        vp9_inverse_transform_b_4x4(xd, xd->plane[0].eobs[ib + iblock[i]],
                                    dqcoeff, diff, 32);
        vp9_inverse_transform_b_4x4(xd, xd->plane[0].eobs[ib + iblock[i] + 1],
                                    dqcoeff + 16, diff + 4, 32);
        i++;
      } else {
        x->fwd_txm4x4(src_diff, coeff, 32);
        x->quantize_b_4x4(x, ib + iblock[i], tx_type, 16);
        vp9_inverse_transform_b_4x4(xd, xd->plane[0].eobs[ib + iblock[i]],
                                    dqcoeff, diff, 32);
      }
    }
  }

  // reconstruct submacroblock
  for (i = 0; i < 4; i++) {
    int16_t* const diff =
        raster_block_offset_int16(xd, BLOCK_SIZE_MB16X16, 0, ib + iblock[i],
                                  xd->plane[0].diff);
    uint8_t* const dst =
        raster_block_offset_uint8(xd, BLOCK_SIZE_MB16X16, 0, ib + iblock[i],
                                  xd->plane[0].dst.buf,
                                  xd->plane[0].dst.stride);
    vp9_recon_b_c(dst, diff, dst, xd->plane[0].dst.stride);
  }
}

void vp9_encode_intra8x8mby(MACROBLOCK *x) {
  int i;

  for (i = 0; i < 4; i++)
    vp9_encode_intra8x8(x, vp9_i8x8_block[i]);
}

static void encode_intra_uv4x4(MACROBLOCK *x, int ib, int mode) {
  MACROBLOCKD * const xd = &x->e_mbd;
  int16_t * const dqcoeff = MB_SUBBLOCK_FIELD(xd, dqcoeff, ib);
  int16_t* const coeff = MB_SUBBLOCK_FIELD(x, coeff, ib);
  const int plane = ib < 20 ? 1 : 2;
  const int block = ib < 20 ? ib - 16 : ib - 20;
  uint8_t* const src =
      raster_block_offset_uint8(xd, BLOCK_SIZE_MB16X16, plane, block,
                                x->plane[plane].src.buf,
                                x->plane[plane].src.stride);
  int16_t* const src_diff =
      raster_block_offset_int16(xd, BLOCK_SIZE_MB16X16, plane, block,
                                x->plane[plane].src_diff);
  int16_t* const diff =
      raster_block_offset_int16(xd, BLOCK_SIZE_MB16X16, plane, block,
                                xd->plane[plane].diff);
  uint8_t* const dst =
      raster_block_offset_uint8(xd, BLOCK_SIZE_MB16X16, plane, block,
                                xd->plane[plane].dst.buf,
                                xd->plane[plane].dst.stride);

  assert(ib >= 16 && ib < 24);
  vp9_intra_uv4x4_predict(&x->e_mbd, ib, mode,
                          dst, xd->plane[plane].dst.stride);

  assert(xd->plane[1].subsampling_x == 1);
  vp9_subtract_block(4, 4, src_diff, 8,
                     src, x->plane[plane].src.stride,
                     dst, xd->plane[plane].dst.stride);

  x->fwd_txm4x4(src_diff, coeff, 16);
  x->quantize_b_4x4(x, ib, DCT_DCT, 16);
  vp9_inverse_transform_b_4x4(&x->e_mbd, xd->plane[plane].eobs[block],
                              dqcoeff, diff, 16);

  vp9_recon_uv_b_c(dst, diff, dst, xd->plane[plane].dst.stride);
}

void vp9_encode_intra8x8mbuv(MACROBLOCK *x) {
  int i;

  for (i = 0; i < 4; i++) {
    int mode = x->e_mbd.mode_info_context->bmi[vp9_i8x8_block[i]].as_mode.first;

    encode_intra_uv4x4(x, i + 16, mode);  // u
    encode_intra_uv4x4(x, i + 20, mode);  // v
  }
}
