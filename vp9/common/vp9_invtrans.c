/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include "vp9/common/vp9_invtrans.h"
#include "./vp9_rtcd.h"

void vp9_inverse_transform_b_4x4(MACROBLOCKD *xd, int eob,
                                 int16_t *dqcoeff, int16_t *diff,
                                 int pitch) {
  if (eob <= 1)
    xd->inv_txm4x4_1(dqcoeff, diff, pitch);
  else
    xd->inv_txm4x4(dqcoeff, diff, pitch);
}

void vp9_inverse_transform_mby_4x4(MACROBLOCKD *xd) {
  int i;

  for (i = 0; i < 16; i++) {
    TX_TYPE tx_type = get_tx_type_4x4(xd, i);
    const int x = i & 3, y = i >> 2;
    if (tx_type != DCT_DCT) {
      vp9_short_iht4x4(BLOCK_OFFSET(xd->plane[0].dqcoeff, i, 16),
                       xd->diff + 64 * y + 4 * x, 16, tx_type);
    } else {
      vp9_inverse_transform_b_4x4(xd,
                                  xd->plane[0].eobs[i],
                                  BLOCK_OFFSET(xd->plane[0].dqcoeff, i, 16),
                                  xd->diff + 64 * y + 4 * x, 32);
    }
  }
}

void vp9_inverse_transform_mbuv_4x4(MACROBLOCKD *xd) {
  int i;

  for (i = 0; i < 4; i++) {
    const int y = i >> 1, x = i & 1;
    vp9_inverse_transform_b_4x4(xd, xd->plane[1].eobs[i],
                                BLOCK_OFFSET(xd->plane[1].dqcoeff, i, 16),
                                xd->diff + 256 + y * 32 + x * 4, 16);
    vp9_inverse_transform_b_4x4(xd, xd->plane[2].eobs[i],
                                BLOCK_OFFSET(xd->plane[2].dqcoeff, i, 16),
                                xd->diff + 320 + y * 32 + x * 4, 16);
  }
}

void vp9_inverse_transform_mb_4x4(MACROBLOCKD *xd) {
  vp9_inverse_transform_mby_4x4(xd);
  vp9_inverse_transform_mbuv_4x4(xd);
}

void vp9_inverse_transform_b_8x8(int16_t *input_dqcoeff, int16_t *output_coeff,
                                 int pitch) {
  vp9_short_idct8x8(input_dqcoeff, output_coeff, pitch);
}

void vp9_inverse_transform_mby_8x8(MACROBLOCKD *xd) {
  int i;

  for (i = 0; i < 4; i++) {
    const int y = i >> 1, x = i & 1;
    TX_TYPE tx_type = get_tx_type_8x8(xd, x * 2 + y * 8);
    if (tx_type != DCT_DCT) {
      vp9_short_iht8x8(BLOCK_OFFSET(xd->plane[0].dqcoeff, i * 4, 16),
                       xd->diff + y * 128 + x * 8, 16, tx_type);
    } else {
      vp9_inverse_transform_b_8x8(BLOCK_OFFSET(xd->plane[0].dqcoeff, i * 4, 16),
                                  xd->diff + y * 128 + x * 8, 32);
    }
  }
}

void vp9_inverse_transform_mbuv_8x8(MACROBLOCKD *xd) {
  vp9_inverse_transform_b_8x8(BLOCK_OFFSET(xd->plane[1].dqcoeff, 0, 16),
                              xd->diff + 256, 16);
  vp9_inverse_transform_b_8x8(BLOCK_OFFSET(xd->plane[2].dqcoeff, 0, 16),
                              xd->diff + 320, 16);
}

void vp9_inverse_transform_mb_8x8(MACROBLOCKD *xd) {
  vp9_inverse_transform_mby_8x8(xd);
  vp9_inverse_transform_mbuv_8x8(xd);
}

void vp9_inverse_transform_b_16x16(int16_t *input_dqcoeff,
                                   int16_t *output_coeff, int pitch) {
  vp9_short_idct16x16(input_dqcoeff, output_coeff, pitch);
}

void vp9_inverse_transform_mby_16x16(MACROBLOCKD *xd) {
  TX_TYPE tx_type = get_tx_type_16x16(xd, 0);
  if (tx_type != DCT_DCT) {
    vp9_short_iht16x16(BLOCK_OFFSET(xd->plane[0].dqcoeff, 0, 16),
                       xd->diff, 16, tx_type);
  } else {
    vp9_inverse_transform_b_16x16(BLOCK_OFFSET(xd->plane[0].dqcoeff, 0, 16),
                                  xd->diff, 32);
  }
}

void vp9_inverse_transform_mb_16x16(MACROBLOCKD *xd) {
  vp9_inverse_transform_mby_16x16(xd);
  vp9_inverse_transform_mbuv_8x8(xd);
}

void vp9_inverse_transform_sby_32x32(MACROBLOCKD *xd, BLOCK_SIZE_TYPE bsize) {
  const int bwl = mb_width_log2(bsize) - 1, bw = 1 << bwl;
  const int bh = 1 << (mb_height_log2(bsize) - 1);
  const int stride = 32 << bwl;
  int n;

  for (n = 0; n < bw * bh; n++) {
    const int x_idx = n & (bw - 1), y_idx = n >> bwl;

    vp9_short_idct32x32(BLOCK_OFFSET(xd->plane[0].dqcoeff, n, 1024),
                        xd->diff + x_idx * 32 + y_idx * 32 * stride,
                        stride * 2);
  }
}

void vp9_inverse_transform_sby_16x16(MACROBLOCKD *xd, BLOCK_SIZE_TYPE bsize) {
  const int bwl = mb_width_log2(bsize), bw = 1 << bwl;
  const int bh = 1 << mb_height_log2(bsize);
  const int stride = 16 << bwl, bstride = 4 << bwl;
  int n;

  for (n = 0; n < bw * bh; n++) {
    const int x_idx = n & (bw - 1), y_idx = n >> bwl;
    const TX_TYPE tx_type = get_tx_type_16x16(xd,
                                              (y_idx * bstride + x_idx) * 4);

    if (tx_type == DCT_DCT) {
      vp9_inverse_transform_b_16x16(BLOCK_OFFSET(xd->plane[0].dqcoeff, n, 256),
                                    xd->diff + x_idx * 16 + y_idx * stride * 16,
                                    stride * 2);
    } else {
      vp9_short_iht16x16(BLOCK_OFFSET(xd->plane[0].dqcoeff, n, 256),
                         xd->diff + x_idx * 16 + y_idx * stride * 16,
                         stride, tx_type);
    }
  }
}

void vp9_inverse_transform_sby_8x8(MACROBLOCKD *xd, BLOCK_SIZE_TYPE bsize) {
  const int bwl = mb_width_log2(bsize) + 1, bw = 1 << bwl;
  const int bh = 1 << (mb_height_log2(bsize) + 1);
  const int stride = 8 << bwl, bstride = 2 << bwl;
  int n;

  for (n = 0; n < bw * bh; n++) {
    const int x_idx = n & (bw - 1), y_idx = n >> bwl;
    const TX_TYPE tx_type = get_tx_type_8x8(xd, (y_idx * bstride + x_idx) * 2);

    if (tx_type == DCT_DCT) {
      vp9_inverse_transform_b_8x8(BLOCK_OFFSET(xd->plane[0].dqcoeff, n, 64),
                                  xd->diff + x_idx * 8 + y_idx * stride * 8,
                                  stride * 2);
    } else {
      vp9_short_iht8x8(BLOCK_OFFSET(xd->plane[0].dqcoeff, n, 64),
                       xd->diff + x_idx * 8 + y_idx * stride * 8,
                       stride, tx_type);
    }
  }
}

void vp9_inverse_transform_sby_4x4(MACROBLOCKD *xd, BLOCK_SIZE_TYPE bsize) {
  const int bwl = mb_width_log2(bsize) + 2, bw = 1 << bwl;
  const int bh = 1 << (mb_height_log2(bsize) + 2);
  const int stride = 4 << bwl, bstride = 1 << bwl;
  int n;

  for (n = 0; n < bw * bh; n++) {
    const int x_idx = n & (bw - 1), y_idx = n >> bwl;
    const TX_TYPE tx_type = get_tx_type_4x4(xd, y_idx * bstride + x_idx);

    if (tx_type == DCT_DCT) {
      vp9_inverse_transform_b_4x4(xd, xd->plane[0].eobs[n],
                                  BLOCK_OFFSET(xd->plane[0].dqcoeff, n, 16),
                                  xd->diff + x_idx * 4 + y_idx * 4 * stride,
                                  stride * 2);
    } else {
      vp9_short_iht4x4(BLOCK_OFFSET(xd->plane[0].dqcoeff, n, 16),
                       xd->diff + x_idx * 4 + y_idx * 4 * stride,
                       stride, tx_type);
    }
  }
}

void vp9_inverse_transform_sbuv_32x32(MACROBLOCKD *xd, BLOCK_SIZE_TYPE bsize) {
  assert(bsize == BLOCK_SIZE_SB64X64);

  vp9_short_idct32x32(xd->plane[1].dqcoeff,
                      xd->diff + 4096, 64);
  vp9_short_idct32x32(xd->plane[2].dqcoeff,
                      xd->diff + 4096 + 1024, 64);
}

void vp9_inverse_transform_sbuv_16x16(MACROBLOCKD *xd, BLOCK_SIZE_TYPE bsize) {
  const int bwl = mb_width_log2(bsize), bhl = mb_height_log2(bsize);
  const int uoff = (16 * 16) << (bwl + bhl), voff = (uoff * 5) >> 2;
  const int bw = 1 << (bwl - 1), bh = 1 << (bhl - 1);
  const int stride = 16 << (bwl - 1);
  int n;

  for (n = 0; n < bw * bh; n++) {
    const int x_idx = n & (bw - 1), y_idx = n >> (bwl - 1);
    const int off = x_idx * 16 + y_idx * stride * 16;

    vp9_inverse_transform_b_16x16(BLOCK_OFFSET(xd->plane[1].dqcoeff, n, 256),
                                  xd->diff + uoff + off, stride * 2);
    vp9_inverse_transform_b_16x16(BLOCK_OFFSET(xd->plane[2].dqcoeff, n, 256),
                                  xd->diff + voff + off, stride * 2);
  }
}

void vp9_inverse_transform_sbuv_8x8(MACROBLOCKD *xd, BLOCK_SIZE_TYPE bsize) {
  const int bwl = mb_width_log2(bsize) + 1, bhl = mb_height_log2(bsize) + 1;
  const int uoff = (8 * 8) << (bwl + bhl), voff = (uoff * 5) >> 2;
  const int bw = 1 << (bwl - 1), bh = 1 << (bhl - 1);
  const int stride = 8 << (bwl - 1);
  int n;

  for (n = 0; n < bw * bh; n++) {
    const int x_idx = n & (bw - 1), y_idx = n >> (bwl - 1);
    const int off = x_idx * 8 + y_idx * stride * 8;

    vp9_inverse_transform_b_8x8(BLOCK_OFFSET(xd->plane[1].dqcoeff, n, 64),
                                xd->diff + uoff + off, stride * 2);
    vp9_inverse_transform_b_8x8(BLOCK_OFFSET(xd->plane[2].dqcoeff, n, 64),
                                xd->diff + voff + off, stride * 2);
  }
}

void vp9_inverse_transform_sbuv_4x4(MACROBLOCKD *xd, BLOCK_SIZE_TYPE bsize) {
  const int bwl = mb_width_log2(bsize) + 2, bhl = mb_height_log2(bsize) + 2;
  const int uoff = (4 * 4) << (bwl + bhl), voff = (uoff * 5) >> 2;
  const int bw = 1 << (bwl - 1), bh = 1 << (bhl - 1);
  const int stride = 4 << (bwl - 1);
  int n;

  for (n = 0; n < bw * bh; n++) {
    const int x_idx = n & (bw - 1), y_idx = n >> (bwl - 1);
    const int off = x_idx * 4 + y_idx * stride * 4;

    vp9_inverse_transform_b_4x4(xd, xd->plane[1].eobs[n],
                                BLOCK_OFFSET(xd->plane[1].dqcoeff, n, 16),
                                xd->diff + uoff + off, stride * 2);
    vp9_inverse_transform_b_4x4(xd, xd->plane[2].eobs[n],
                                BLOCK_OFFSET(xd->plane[2].dqcoeff, n, 16),
                                xd->diff + voff + off, stride * 2);
  }
}
