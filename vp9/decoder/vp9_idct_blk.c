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
#include "vp9/decoder/vp9_idct_blk.h"

void vp9_idct_add_y_block_c(int16_t *q, uint8_t *dst, int stride,
                            MACROBLOCKD *xd) {
  int i, j;

  for (i = 0; i < 4; i++) {
    for (j = 0; j < 4; j++) {
      vp9_idct_add(q, dst, stride, xd->plane[0].eobs[i * 4  + j]);
      q   += 16;
      dst += 4;
    }

    dst += 4 * stride - 16;
  }
}

void vp9_idct_add_uv_block_c(int16_t *q, uint8_t *dst, int stride,
                             uint16_t *eobs) {
  int i, j;

  for (i = 0; i < 2; i++) {
    for (j = 0; j < 2; j++) {
      vp9_idct_add(q, dst, stride, eobs[i * 2 + j]);
      q   += 16;
      dst += 4;
    }

    dst += 4 * stride - 8;
  }
}

void vp9_idct_add_y_block_8x8_c(int16_t *q, uint8_t *dst, int stride,
                                MACROBLOCKD *xd) {
  uint8_t *origdest = dst;

  vp9_idct_add_8x8_c(q, dst, stride, xd->plane[0].eobs[0]);
  vp9_idct_add_8x8_c(&q[64], origdest + 8, stride, xd->plane[0].eobs[4]);
  vp9_idct_add_8x8_c(&q[128], origdest + 8 * stride, stride,
                     xd->plane[0].eobs[8]);
  vp9_idct_add_8x8_c(&q[192], origdest + 8 * stride + 8, stride,
                     xd->plane[0].eobs[12]);
}

void vp9_idct_add_y_block_lossless_c(int16_t *q, uint8_t *dst, int stride,
                                     MACROBLOCKD *xd) {
  int i, j;

  for (i = 0; i < 4; i++) {
    for (j = 0; j < 4; j++) {
      vp9_idct_add_lossless_c(q, dst, stride, xd->plane[0].eobs[i * 4 + j]);
      q   += 16;
      dst += 4;
    }

    dst += 4 * stride - 16;
  }
}

void vp9_idct_add_uv_block_lossless_c(int16_t *q, uint8_t *dst, int stride,
                                      uint16_t *eobs) {
  int i, j;

  for (i = 0; i < 2; i++) {
    for (j = 0; j < 2; j++) {
      vp9_idct_add_lossless_c(q, dst, stride, eobs[i * 2 + j]);
      q   += 16;
      dst += 4;
    }

    dst += 4 * stride - 8;
  }
}

static void add_constant_residual(const int16_t diff, uint8_t *dest, int stride,
                                  int width, int height) {
  int r, c;

  for (r = 0; r < height; r++) {
    for (c = 0; c < width; c++)
      dest[c] = clip_pixel(diff + dest[c]);

    dest += stride;
  }
}

void vp9_add_constant_residual_8x8_c(const int16_t diff, uint8_t *dest,
                                     int stride) {
  add_constant_residual(diff, dest, stride, 8, 8);
}

void vp9_add_constant_residual_16x16_c(const int16_t diff, uint8_t *dest,
                                       int stride) {
  add_constant_residual(diff, dest, stride, 16, 16);
}

void vp9_add_constant_residual_32x32_c(const int16_t diff,  uint8_t *dest,
                                       int stride) {
  add_constant_residual(diff, dest, stride, 32, 32);
}

void vp9_iht_add_c(TX_TYPE tx_type, int16_t *input, uint8_t *dest, int stride,
                   int eob) {
  if (tx_type == DCT_DCT) {
    vp9_idct_add(input, dest, stride, eob);
  } else {
    vp9_short_iht4x4_add(input, dest, stride, tx_type);
    vpx_memset(input, 0, 32);
  }
}

void vp9_iht_add_8x8_c(TX_TYPE tx_type, int16_t *input, uint8_t *dest,
                       int stride, int eob) {
  if (tx_type == DCT_DCT) {
    vp9_idct_add_8x8(input, dest, stride, eob);
  } else {
    if (eob > 0) {
      vp9_short_iht8x8_add(input, dest, stride, tx_type);
      vpx_memset(input, 0, 128);
    }
  }
}

void vp9_idct_add_c(int16_t *input, uint8_t *dest, int stride, int eob) {
  if (eob > 1) {
    vp9_short_idct4x4_add(input, dest, stride);
    vpx_memset(input, 0, 32);
  } else {
    vp9_dc_only_idct_add(input[0], dest, dest, stride, stride);
    ((int *)input)[0] = 0;
  }
}

void vp9_dc_idct_add_c(int16_t *input, uint8_t *dest, int stride, int dc) {
  input[0] = dc;
  vp9_short_idct4x4_add(input, dest, stride);
  vpx_memset(input, 0, 32);
}

void vp9_idct_add_lossless_c(int16_t *input, uint8_t *dest, int stride,
                             int eob) {
  if (eob > 1) {
    vp9_short_iwalsh4x4_add(input, dest, stride);
    vpx_memset(input, 0, 32);
  } else {
    vp9_short_iwalsh4x4_1_add_c(input, dest, stride);
    ((int *)input)[0] = 0;
  }
}

void vp9_dc_idct_add_lossless_c(int16_t *input, uint8_t *dest,
                                int stride, int dc) {
  input[0] = dc;
  vp9_short_iwalsh4x4_add(input, dest, stride);
  vpx_memset(input, 0, 32);
}

void vp9_idct_add_8x8_c(int16_t *input, uint8_t *dest, int stride, int eob) {
  // If dc is 1, then input[0] is the reconstructed value, do not need
  // dequantization. Also, when dc is 1, dc is counted in eobs, namely eobs >=1.

  // The calculation can be simplified if there are not many non-zero dct
  // coefficients. Use eobs to decide what to do.
  // TODO(yunqingwang): "eobs = 1" case is also handled in vp9_short_idct8x8_c.
  // Combine that with code here.
  if (eob) {
    if (eob == 1) {
      // DC only DCT coefficient
      int16_t in = input[0];
      int16_t out;

      // Note: the idct1 will need to be modified accordingly whenever
      // vp9_short_idct8x8_c() is modified.
      vp9_short_idct1_8x8_c(&in, &out);
      input[0] = 0;

      vp9_add_constant_residual_8x8(out, dest, stride);
#if !CONFIG_SCATTERSCAN
    } else if (eob <= 10) {
      vp9_short_idct10_8x8_add(input, dest, stride);
      input[0] = input[1] = input[2] = input[3] = 0;
      input[8] = input[9] = input[10] = 0;
      input[16] = input[17] = 0;
      input[24] = 0;
#endif
    } else {
      vp9_short_idct8x8_add(input, dest, stride);
      vpx_memset(input, 0, 128);
    }
  }
}

void vp9_iht_add_16x16_c(TX_TYPE tx_type, int16_t *input, uint8_t *dest,
                         int stride, int eob) {
  if (tx_type == DCT_DCT) {
    vp9_idct_add_16x16(input, dest, stride, eob);
  } else {
    if (eob > 0) {
      vp9_short_iht16x16_add(input, dest, stride, tx_type);
      vpx_memset(input, 0, 512);
    }
  }
}

void vp9_idct_add_16x16_c(int16_t *input, uint8_t *dest, int stride, int eob) {
  /* The calculation can be simplified if there are not many non-zero dct
   * coefficients. Use eobs to separate different cases. */
  if (eob) {
    if (eob == 1) {
      /* DC only DCT coefficient. */
      int16_t in = input[0];
      int16_t out;
      /* Note: the idct1 will need to be modified accordingly whenever
       * vp9_short_idct16x16() is modified. */
      vp9_short_idct1_16x16_c(&in, &out);
      input[0] = 0;

      vp9_add_constant_residual_16x16(out, dest, stride);
#if !CONFIG_SCATTERSCAN
    } else if (eob <= 10) {
      vp9_short_idct10_16x16_add(input, dest, stride);
      input[0] = input[1] = input[2] = input[3] = 0;
      input[16] = input[17] = input[18] = 0;
      input[32] = input[33] = 0;
      input[48] = 0;
#endif
    } else {
      vp9_short_idct16x16_add(input, dest, stride);
      vpx_memset(input, 0, 512);
    }
  }
}

void vp9_idct_add_32x32_c(int16_t *input, uint8_t *dest, int stride, int eob) {
  DECLARE_ALIGNED_ARRAY(16, int16_t, output, 1024);

  if (eob) {
    if (eob == 1) {
      vp9_short_idct1_32x32(input, output);
      vp9_add_constant_residual_32x32(output[0], dest, stride);
      input[0] = 0;
#if !CONFIG_SCATTERSCAN
    } else if (eob <= 10) {
      vp9_short_idct10_32x32_add_c(input, dest, stride);
      input[0] = input[1] = input[2] = input[3] = 0;
      input[32] = input[33] = input[34] = 0;
      input[64] = input[65] = 0;
      input[96] = 0;

#endif
    } else {
      vp9_short_idct32x32_add(input, dest, stride);
      vpx_memset(input, 0, 2048);
    }
  }
}

