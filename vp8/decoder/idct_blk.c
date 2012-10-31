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
#include "dequantize.h"

void vp9_dequant_dc_idct_add_c(short *input, short *dq, unsigned char *pred,
                               unsigned char *dest, int pitch, int stride,
                               int Dc);
void vp9_dequant_idct_add_c(short *input, short *dq, unsigned char *pred,
                            unsigned char *dest, int pitch, int stride);
void vp9_dc_only_idct_add_c(short input_dc, unsigned char *pred_ptr,
                            unsigned char *dst_ptr, int pitch, int stride);
#if CONFIG_LOSSLESS
void vp9_dequant_idct_add_lossless_c(short *input, short *dq,
                                     unsigned char *pred, unsigned char *dest,
                                     int pitch, int stride);
void vp9_dc_only_idct_add_lossless_c(short input_dc, unsigned char *pred_ptr,
                                     unsigned char *dst_ptr,
                                     int pitch, int stride);
#endif

void vp9_dequant_dc_idct_add_y_block_c(short *q, short *dq,
                                       unsigned char *pre,
                                       unsigned char *dst,
                                       int stride, char *eobs,
                                       short *dc) {
  int i, j;

  for (i = 0; i < 4; i++) {
    for (j = 0; j < 4; j++) {
      if (*eobs++ > 1)
        vp9_dequant_dc_idct_add_c(q, dq, pre, dst, 16, stride, dc[0]);
      else
        vp9_dc_only_idct_add_c(dc[0], pre, dst, 16, stride);

      q   += 16;
      pre += 4;
      dst += 4;
      dc++;
    }

    pre += 64 - 16;
    dst += 4 * stride - 16;
  }
}

void vp9_dequant_idct_add_y_block_c(short *q, short *dq,
                                    unsigned char *pre,
                                    unsigned char *dst,
                                    int stride, char *eobs) {
  int i, j;

  for (i = 0; i < 4; i++) {
    for (j = 0; j < 4; j++) {
      if (*eobs++ > 1)
        vp9_dequant_idct_add_c(q, dq, pre, dst, 16, stride);
      else {
        vp9_dc_only_idct_add_c(q[0]*dq[0], pre, dst, 16, stride);
        ((int *)q)[0] = 0;
      }

      q   += 16;
      pre += 4;
      dst += 4;
    }

    pre += 64 - 16;
    dst += 4 * stride - 16;
  }
}

void vp9_dequant_idct_add_uv_block_c(short *q, short *dq, unsigned char *pre,
                                     unsigned char *dstu, unsigned char *dstv,
                                     int stride, char *eobs) {
  int i, j;

  for (i = 0; i < 2; i++) {
    for (j = 0; j < 2; j++) {
      if (*eobs++ > 1)
        vp9_dequant_idct_add_c(q, dq, pre, dstu, 8, stride);
      else {
        vp9_dc_only_idct_add_c(q[0]*dq[0], pre, dstu, 8, stride);
        ((int *)q)[0] = 0;
      }

      q    += 16;
      pre  += 4;
      dstu += 4;
    }

    pre  += 32 - 8;
    dstu += 4 * stride - 8;
  }

  for (i = 0; i < 2; i++) {
    for (j = 0; j < 2; j++) {
      if (*eobs++ > 1)
        vp9_dequant_idct_add_c(q, dq, pre, dstv, 8, stride);
      else {
        vp9_dc_only_idct_add_c(q[0]*dq[0], pre, dstv, 8, stride);
        ((int *)q)[0] = 0;
      }

      q    += 16;
      pre  += 4;
      dstv += 4;
    }

    pre  += 32 - 8;
    dstv += 4 * stride - 8;
  }
}


void vp9_dequant_dc_idct_add_y_block_8x8_c(short *q, short *dq,
                                           unsigned char *pre,
                                           unsigned char *dst,
                                           int stride, char *eobs, short *dc,
                                           MACROBLOCKD *xd) {
  vp9_dequant_dc_idct_add_8x8_c(q, dq, pre, dst, 16, stride, dc[0]);
  vp9_dequant_dc_idct_add_8x8_c(&q[64], dq, pre + 8, dst + 8, 16, stride, dc[1]);
  vp9_dequant_dc_idct_add_8x8_c(&q[128], dq, pre + 8 * 16,
                                dst + 8 * stride, 16, stride, dc[4]);
  vp9_dequant_dc_idct_add_8x8_c(&q[192], dq, pre + 8 * 16 + 8,
                                dst + 8 * stride + 8, 16, stride, dc[8]);
}

#if CONFIG_SUPERBLOCKS
void vp9_dequant_dc_idct_add_y_block_8x8_inplace_c(short *q, short *dq,
                                                   unsigned char *dst,
                                                   int stride, char *eobs,
                                                   short *dc, MACROBLOCKD *xd) {
  vp9_dequant_dc_idct_add_8x8_c(q, dq, dst, dst, stride, stride, dc[0]);
  vp9_dequant_dc_idct_add_8x8_c(&q[64], dq, dst + 8,
                                dst + 8, stride, stride, dc[1]);
  vp9_dequant_dc_idct_add_8x8_c(&q[128], dq, dst + 8 * stride,
                                dst + 8 * stride, stride, stride, dc[4]);
  vp9_dequant_dc_idct_add_8x8_c(&q[192], dq, dst + 8 * stride + 8,
                                dst + 8 * stride + 8, stride, stride, dc[8]);
}
#endif

void vp9_dequant_idct_add_y_block_8x8_c(short *q, short *dq,
                                        unsigned char *pre,
                                        unsigned char *dst,
                                        int stride, char *eobs,
                                        MACROBLOCKD *xd) {
  unsigned char *origdest = dst;
  unsigned char *origpred = pre;

  vp9_dequant_idct_add_8x8_c(q, dq, pre, dst, 16, stride);
  vp9_dequant_idct_add_8x8_c(&q[64], dq, origpred + 8,
                             origdest + 8, 16, stride);
  vp9_dequant_idct_add_8x8_c(&q[128], dq, origpred + 8 * 16,
                             origdest + 8 * stride, 16, stride);
  vp9_dequant_idct_add_8x8_c(&q[192], dq, origpred + 8 * 16 + 8,
                             origdest + 8 * stride + 8, 16, stride);
}

void vp9_dequant_idct_add_uv_block_8x8_c(short *q, short *dq,
                                         unsigned char *pre,
                                         unsigned char *dstu,
                                         unsigned char *dstv,
                                         int stride, char *eobs,
                                         MACROBLOCKD *xd) {
  vp9_dequant_idct_add_8x8_c(q, dq, pre, dstu, 8, stride);

  q    += 64;
  pre  += 64;

  vp9_dequant_idct_add_8x8_c(q, dq, pre, dstv, 8, stride);
}

#if CONFIG_SUPERBLOCKS
void vp9_dequant_idct_add_uv_block_8x8_inplace_c(short *q, short *dq,
                                                 unsigned char *dstu,
                                                 unsigned char *dstv,
                                                 int stride, char *eobs,
                                                 MACROBLOCKD *xd) {
  vp9_dequant_idct_add_8x8_c(q, dq, dstu, dstu, stride, stride);

  q    += 64;

  vp9_dequant_idct_add_8x8_c(q, dq, dstv, dstv, stride, stride);
}
#endif

#if CONFIG_LOSSLESS
void vp9_dequant_dc_idct_add_y_block_lossless_c(short *q, short *dq,
                                                unsigned char *pre,
                                                unsigned char *dst,
                                                int stride, char *eobs,
                                                short *dc) {
  int i, j;

  for (i = 0; i < 4; i++) {
    for (j = 0; j < 4; j++) {
      if (*eobs++ > 1)
        vp9_dequant_dc_idct_add_lossless_c(q, dq, pre, dst, 16, stride, dc[0]);
      else
        vp9_dc_only_inv_walsh_add_c(dc[0], pre, dst, 16, stride);

      q   += 16;
      pre += 4;
      dst += 4;
      dc++;
    }

    pre += 64 - 16;
    dst += 4 * stride - 16;
  }
}

void vp9_dequant_idct_add_y_block_lossless_c(short *q, short *dq,
                                             unsigned char *pre,
                                             unsigned char *dst,
                                             int stride, char *eobs) {
  int i, j;

  for (i = 0; i < 4; i++) {
    for (j = 0; j < 4; j++) {
      if (*eobs++ > 1)
        vp9_dequant_idct_add_lossless_c(q, dq, pre, dst, 16, stride);
      else {
        vp9_dc_only_inv_walsh_add_c(q[0]*dq[0], pre, dst, 16, stride);
        ((int *)q)[0] = 0;
      }

      q   += 16;
      pre += 4;
      dst += 4;
    }

    pre += 64 - 16;
    dst += 4 * stride - 16;
  }
}

void vp9_dequant_idct_add_uv_block_lossless_c(short *q, short *dq,
                                              unsigned char *pre,
                                              unsigned char *dstu,
                                              unsigned char *dstv,
                                              int stride, char *eobs) {
  int i, j;

  for (i = 0; i < 2; i++) {
    for (j = 0; j < 2; j++) {
      if (*eobs++ > 1)
        vp9_dequant_idct_add_lossless_c(q, dq, pre, dstu, 8, stride);
      else {
        vp9_dc_only_inv_walsh_add_c(q[0]*dq[0], pre, dstu, 8, stride);
        ((int *)q)[0] = 0;
      }

      q    += 16;
      pre  += 4;
      dstu += 4;
    }

    pre  += 32 - 8;
    dstu += 4 * stride - 8;
  }

  for (i = 0; i < 2; i++) {
    for (j = 0; j < 2; j++) {
      if (*eobs++ > 1)
        vp9_dequant_idct_add_lossless_c(q, dq, pre, dstv, 8, stride);
      else {
        vp9_dc_only_inv_walsh_add_c(q[0]*dq[0], pre, dstv, 8, stride);
        ((int *)q)[0] = 0;
      }

      q    += 16;
      pre  += 4;
      dstv += 4;
    }

    pre  += 32 - 8;
    dstv += 4 * stride - 8;
  }
}
#endif

