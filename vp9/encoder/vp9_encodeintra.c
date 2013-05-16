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
  int16_t* const coeff = BLOCK_OFFSET(x->plane[0].coeff, ib, 16);
  const int bwl = b_width_log2(bsize), bhl = b_height_log2(bsize);

  assert(ib < (1 << (bwl + bhl)));

  vp9_intra4x4_predict(&x->e_mbd, ib, bsize,
                       xd->mode_info_context->bmi[ib].as_mode.first,
                       dst, xd->plane[0].dst.stride);
  vp9_subtract_block(4, 4, src_diff, 4 << bwl,
                     src, x->plane[0].src.stride,
                     dst, xd->plane[0].dst.stride);

  tx_type = get_tx_type_4x4(&x->e_mbd, ib);
  if (tx_type != DCT_DCT) {
    vp9_short_fht4x4(src_diff, coeff, 4 << bwl, tx_type);
    x->quantize_b_4x4(x, ib, tx_type, 16);
    vp9_short_iht4x4_add(BLOCK_OFFSET(xd->plane[0].dqcoeff, ib, 16), dst,
                         xd->plane[0].dst.stride, tx_type);
  } else {
    x->fwd_txm4x4(src_diff, coeff, 8 << bwl);
    x->quantize_b_4x4(x, ib, tx_type, 16);
    vp9_inverse_transform_b_4x4_add(&x->e_mbd, xd->plane[0].eobs[ib],
                                BLOCK_OFFSET(xd->plane[0].dqcoeff, ib, 16),
                                dst, xd->plane[0].dst.stride);
  }
}

void vp9_encode_intra16x16mby(VP9_COMMON *const cm, MACROBLOCK *x) {
  MACROBLOCKD *xd = &x->e_mbd;

  vp9_build_intra_predictors_sby_s(xd, BLOCK_SIZE_MB16X16);
  vp9_encode_sby(cm, x, BLOCK_SIZE_MB16X16);
}

void vp9_encode_intra16x16mbuv(VP9_COMMON *const cm, MACROBLOCK *x) {
  MACROBLOCKD *xd = &x->e_mbd;

  vp9_build_intra_predictors_sbuv_s(xd, BLOCK_SIZE_MB16X16);
  vp9_encode_sbuv(cm, x, BLOCK_SIZE_MB16X16);
}


