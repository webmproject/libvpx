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

void vp9_inverse_transform_sby_32x32(MACROBLOCKD *xd) {
  vp9_short_idct32x32(BLOCK_OFFSET(xd->plane[0].dqcoeff, 0, 16), xd->diff, 64);
}

void vp9_inverse_transform_sby_16x16(MACROBLOCKD *xd) {
  int n;

  for (n = 0; n < 4; n++) {
    const int x_idx = n & 1, y_idx = n >> 1;
    const TX_TYPE tx_type = get_tx_type_16x16(xd, (y_idx * 8 + x_idx) * 4);

    if (tx_type == DCT_DCT) {
      vp9_inverse_transform_b_16x16(BLOCK_OFFSET(xd->plane[0].dqcoeff, n, 256),
                                    xd->diff + x_idx * 16 + y_idx * 32 * 16,
                                    64);
    } else {
      vp9_short_iht16x16(BLOCK_OFFSET(xd->plane[0].dqcoeff, n, 256),
                         xd->diff + x_idx * 16 + y_idx * 32 * 16, 32, tx_type);
    }
  }
}

void vp9_inverse_transform_sby_8x8(MACROBLOCKD *xd) {
  int n;

  for (n = 0; n < 16; n++) {
    const int x_idx = n & 3, y_idx = n >> 2;
    const TX_TYPE tx_type = get_tx_type_8x8(xd, (y_idx * 8 + x_idx) * 2);

    if (tx_type == DCT_DCT) {
      vp9_inverse_transform_b_8x8(BLOCK_OFFSET(xd->plane[0].dqcoeff, n, 64),
                                  xd->diff + x_idx * 8 + y_idx * 32 * 8, 64);
    } else {
      vp9_short_iht8x8(BLOCK_OFFSET(xd->plane[0].dqcoeff, n, 64),
                       xd->diff + x_idx * 8 + y_idx * 32 * 8, 32, tx_type);
    }
  }
}

void vp9_inverse_transform_sby_4x4(MACROBLOCKD *xd) {
  int n;

  for (n = 0; n < 64; n++) {
    const int x_idx = n & 7, y_idx = n >> 3;
    const TX_TYPE tx_type = get_tx_type_4x4(xd, y_idx * 8 + x_idx);

    if (tx_type == DCT_DCT) {
      vp9_inverse_transform_b_4x4(xd, xd->plane[0].eobs[n],
                                  BLOCK_OFFSET(xd->plane[0].dqcoeff, n, 16),
                                  xd->diff + x_idx * 4 + y_idx * 4 * 32, 64);
    } else {
      vp9_short_iht4x4(BLOCK_OFFSET(xd->plane[0].dqcoeff, n, 16),
                       xd->diff + x_idx * 4 + y_idx * 4 * 32, 32, tx_type);
    }
  }
}

void vp9_inverse_transform_sbuv_16x16(MACROBLOCKD *xd) {
  vp9_inverse_transform_b_16x16(xd->plane[1].dqcoeff,
                                xd->diff + 1024, 32);
  vp9_inverse_transform_b_16x16(xd->plane[2].dqcoeff,
                                xd->diff + 1280, 32);
}

void vp9_inverse_transform_sbuv_8x8(MACROBLOCKD *xd) {
  int n;

  for (n = 0; n < 4; n++) {
    const int x_idx = n & 1, y_idx = n >> 1;

    vp9_inverse_transform_b_8x8(BLOCK_OFFSET(xd->plane[1].dqcoeff, n, 64),
                                xd->diff + 1024 + x_idx * 8 + y_idx * 16 * 8,
                                32);
    vp9_inverse_transform_b_8x8(BLOCK_OFFSET(xd->plane[2].dqcoeff, n, 64),
                                xd->diff + 1280 + x_idx * 8 + y_idx * 16 * 8,
                                32);
  }
}

void vp9_inverse_transform_sbuv_4x4(MACROBLOCKD *xd) {
  int n;

  for (n = 0; n < 16; n++) {
    const int x_idx = n & 3, y_idx = n >> 2;

    vp9_inverse_transform_b_4x4(xd, xd->plane[1].eobs[n],
                                BLOCK_OFFSET(xd->plane[1].dqcoeff, n, 16),
                                xd->diff + 1024 + x_idx * 4 + y_idx * 16 * 4,
                                32);
    vp9_inverse_transform_b_4x4(xd, xd->plane[2].eobs[n],
                                BLOCK_OFFSET(xd->plane[2].dqcoeff, n, 16),
                                xd->diff + 1280 + x_idx * 4 + y_idx * 16 * 4,
                                32);
  }
}

void vp9_inverse_transform_sb64y_32x32(MACROBLOCKD *xd) {
  int n;

  for (n = 0; n < 4; n++) {
    const int x_idx = n & 1, y_idx = n >> 1;

    vp9_short_idct32x32(BLOCK_OFFSET(xd->plane[0].dqcoeff, n, 1024),
                        xd->diff + x_idx * 32 + y_idx * 32 * 64, 128);
  }
}

void vp9_inverse_transform_sb64y_16x16(MACROBLOCKD *xd) {
  int n;

  for (n = 0; n < 16; n++) {
    const int x_idx = n & 3, y_idx = n >> 2;
    const TX_TYPE tx_type = get_tx_type_16x16(xd, (y_idx * 16 + x_idx) * 4);

    if (tx_type == DCT_DCT) {
      vp9_inverse_transform_b_16x16(BLOCK_OFFSET(xd->plane[0].dqcoeff, n, 256),
                                    xd->diff + x_idx * 16 + y_idx * 64 * 16,
                                    128);
    } else {
      vp9_short_iht16x16(BLOCK_OFFSET(xd->plane[0].dqcoeff, n, 256),
                         xd->diff + x_idx * 16 + y_idx * 64 * 16, 64, tx_type);
    }
  }
}

void vp9_inverse_transform_sb64y_8x8(MACROBLOCKD *xd) {
  int n;

  for (n = 0; n < 64; n++) {
    const int x_idx = n & 7, y_idx = n >> 3;
    const TX_TYPE tx_type = get_tx_type_8x8(xd, (y_idx * 16 + x_idx) * 2);

    if (tx_type == DCT_DCT) {
      vp9_inverse_transform_b_8x8(BLOCK_OFFSET(xd->plane[0].dqcoeff, n, 64),
                                  xd->diff + x_idx * 8 + y_idx * 64 * 8, 128);
    } else {
      vp9_short_iht8x8(BLOCK_OFFSET(xd->plane[0].dqcoeff, n, 64),
                       xd->diff + x_idx * 8 + y_idx * 64 * 8, 64, tx_type);
    }
  }
}

void vp9_inverse_transform_sb64y_4x4(MACROBLOCKD *xd) {
  int n;

  for (n = 0; n < 256; n++) {
    const int x_idx = n & 15, y_idx = n >> 4;
    const TX_TYPE tx_type = get_tx_type_4x4(xd, y_idx * 16 + x_idx);

    if (tx_type == DCT_DCT) {
      vp9_inverse_transform_b_4x4(xd, xd->plane[0].eobs[n],
                                  BLOCK_OFFSET(xd->plane[0].dqcoeff, n, 16),
                                  xd->diff + x_idx * 4 + y_idx * 4 * 64, 128);
    } else {
      vp9_short_iht4x4(BLOCK_OFFSET(xd->plane[0].dqcoeff, n, 16),
                       xd->diff + x_idx * 4 + y_idx * 4 * 64, 64, tx_type);
    }
  }
}

void vp9_inverse_transform_sb64uv_32x32(MACROBLOCKD *xd) {
  vp9_short_idct32x32(xd->plane[1].dqcoeff,
                      xd->diff + 4096, 64);
  vp9_short_idct32x32(xd->plane[2].dqcoeff,
                      xd->diff + 4096 + 1024, 64);
}

void vp9_inverse_transform_sb64uv_16x16(MACROBLOCKD *xd) {
  int n;

  for (n = 0; n < 4; n++) {
    const int x_idx = n & 1, y_idx = n >> 1, off = x_idx * 16 + y_idx * 32 * 16;

    vp9_inverse_transform_b_16x16(BLOCK_OFFSET(xd->plane[1].dqcoeff, n, 256),
                                  xd->diff + 4096 + off, 64);
    vp9_inverse_transform_b_16x16(BLOCK_OFFSET(xd->plane[2].dqcoeff, n, 256),
                                  xd->diff + 4096 + 1024 + off, 64);
  }
}

void vp9_inverse_transform_sb64uv_8x8(MACROBLOCKD *xd) {
  int n;

  for (n = 0; n < 16; n++) {
    const int x_idx = n & 3, y_idx = n >> 2, off = x_idx * 8 + y_idx * 32 * 8;

    vp9_inverse_transform_b_8x8(BLOCK_OFFSET(xd->plane[1].dqcoeff, n, 64),
                                xd->diff + 4096 + off, 64);
    vp9_inverse_transform_b_8x8(BLOCK_OFFSET(xd->plane[2].dqcoeff, n, 64),
                                xd->diff + 4096 + 1024 + off, 64);
  }
}

void vp9_inverse_transform_sb64uv_4x4(MACROBLOCKD *xd) {
  int n;

  for (n = 0; n < 64; n++) {
    const int x_idx = n & 7, y_idx = n >> 3, off = x_idx * 4 + y_idx * 32 * 4;

    vp9_inverse_transform_b_4x4(xd, xd->plane[1].eobs[n],
                                BLOCK_OFFSET(xd->plane[1].dqcoeff, n, 16),
                                xd->diff + 4096 + off, 64);
    vp9_inverse_transform_b_4x4(xd, xd->plane[2].eobs[n],
                                BLOCK_OFFSET(xd->plane[2].dqcoeff, n, 16),
                                xd->diff + 4096 + 1024 + off, 64);
  }
}
