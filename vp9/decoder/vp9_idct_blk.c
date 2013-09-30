/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include "./vp9_rtcd.h"
#include "vp9/common/vp9_blockd.h"
#include "vp9/decoder/vp9_idct_blk.h"

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
    vp9_short_idct4x4_1_add(input, dest, stride);
    ((int *)input)[0] = 0;
  }
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
      vp9_short_idct8x8_1_add(input, dest, stride);
      input[0] = 0;
    } else if (eob <= 10) {
      vp9_short_idct8x8_10_add(input, dest, stride);
      vpx_memset(input, 0, 128);
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
      vp9_short_idct16x16_1_add(input, dest, stride);
      input[0] = 0;
    } else if (eob <= 10) {
      vp9_short_idct16x16_10_add(input, dest, stride);
      vpx_memset(input, 0, 512);
    } else {
      vp9_short_idct16x16_add(input, dest, stride);
      vpx_memset(input, 0, 512);
    }
  }
}

void vp9_idct_add_32x32_c(int16_t *input, uint8_t *dest, int stride, int eob) {
  if (eob) {
    if (eob == 1) {
      vp9_short_idct32x32_1_add(input, dest, stride);
      input[0] = 0;
    } else {
      vp9_short_idct32x32_add(input, dest, stride);
      vpx_memset(input, 0, 2048);
    }
  }
}

