/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include "vp9_rtcd.h"
#include "vp9/common/vp9_blockd.h"
#include "vp9/decoder/vp9_dequantize.h"

void vp9_dequant_idct_add_y_block_c(int16_t *q, const int16_t *dq,
                                    uint8_t *pre,
                                    uint8_t *dst,
                                    int stride, MACROBLOCKD *xd) {
  int i, j;

  for (i = 0; i < 4; i++) {
    for (j = 0; j < 4; j++) {
      vp9_dequant_idct_add(q, dq, pre, dst, 16, stride,
                           xd->plane[0].eobs[i * 4  + j]);
      q   += 16;
      pre += 4;
      dst += 4;
    }

    pre += 64 - 16;
    dst += 4 * stride - 16;
  }
}

void vp9_dequant_idct_add_uv_block_c(int16_t *q, const int16_t *dq,
                                     uint8_t *pre, uint8_t *dst,
                                     int stride, uint16_t *eobs) {
  int i, j;

  for (i = 0; i < 2; i++) {
    for (j = 0; j < 2; j++) {
      vp9_dequant_idct_add(q, dq, pre, dst, 8, stride, eobs[i * 2 + j]);
      q   += 16;
      pre += 4;
      dst += 4;
    }

    pre  += 32 - 8;
    dst += 4 * stride - 8;
  }
}

void vp9_dequant_idct_add_y_block_8x8_c(int16_t *q, const int16_t *dq,
                                        uint8_t *pre,
                                        uint8_t *dst,
                                        int stride, MACROBLOCKD *xd) {
  uint8_t *origdest = dst;
  uint8_t *origpred = pre;

  vp9_dequant_idct_add_8x8_c(q, dq, pre, dst, 16, stride,
                             xd->plane[0].eobs[0]);
  vp9_dequant_idct_add_8x8_c(&q[64], dq, origpred + 8,
                             origdest + 8, 16, stride,
                             xd->plane[0].eobs[4]);
  vp9_dequant_idct_add_8x8_c(&q[128], dq, origpred + 8 * 16,
                             origdest + 8 * stride, 16, stride,
                             xd->plane[0].eobs[8]);
  vp9_dequant_idct_add_8x8_c(&q[192], dq, origpred + 8 * 16 + 8,
                             origdest + 8 * stride + 8, 16, stride,
                             xd->plane[0].eobs[12]);
}

void vp9_dequant_idct_add_y_block_lossless_c(int16_t *q, const int16_t *dq,
                                             uint8_t *pre,
                                             uint8_t *dst,
                                             int stride, MACROBLOCKD *xd) {
  int i, j;

  for (i = 0; i < 4; i++) {
    for (j = 0; j < 4; j++) {
      vp9_dequant_idct_add_lossless_c(q, dq, pre, dst, 16, stride,
                                      xd->plane[0].eobs[i * 4 + j]);
      q   += 16;
      pre += 4;
      dst += 4;
    }

    pre += 64 - 16;
    dst += 4 * stride - 16;
  }
}

void vp9_dequant_idct_add_uv_block_lossless_c(int16_t *q, const int16_t *dq,
                                              uint8_t *pre,
                                              uint8_t *dst,
                                              int stride,
                                              uint16_t *eobs) {
  int i, j;

  for (i = 0; i < 2; i++) {
    for (j = 0; j < 2; j++) {
      vp9_dequant_idct_add_lossless_c(q, dq, pre, dst, 8, stride,
                                      eobs[i * 2 + j]);
      q   += 16;
      pre += 4;
      dst += 4;
    }

    pre += 32 - 8;
    dst += 4 * stride - 8;
  }
}

