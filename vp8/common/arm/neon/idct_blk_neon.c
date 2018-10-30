/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include "vpx_config.h"
#include "vp8_rtcd.h"

/* place these declarations here because we don't want to maintain them
 * outside of this scope
 */
void idct_dequant_full_2x_neon(short *q, short *dq, unsigned char *dst,
                               int stride);
void idct_dequant_0_2x_neon(short *q, short dq, unsigned char *dst, int stride);

void vp8_dequant_idct_add_y_block_neon(short *q, short *dq, unsigned char *dst,
                                       int stride, char *eobs) {
  int i;

  for (i = 0; i < 4; ++i) {
    if (((short *)(eobs))[0]) {
      if (((short *)eobs)[0] & 0xfefe)
        idct_dequant_full_2x_neon(q, dq, dst, stride);
      else
        idct_dequant_0_2x_neon(q, dq[0], dst, stride);
    }

    if (((short *)(eobs))[1]) {
      if (((short *)eobs)[1] & 0xfefe)
        idct_dequant_full_2x_neon(q + 32, dq, dst + 8, stride);
      else
        idct_dequant_0_2x_neon(q + 32, dq[0], dst + 8, stride);
    }
    q += 64;
    dst += 4 * stride;
    eobs += 4;
  }
}

void vp8_dequant_idct_add_uv_block_neon(short *q, short *dq,
                                        unsigned char *dst_u,
                                        unsigned char *dst_v, int stride,
                                        char *eobs) {
  if (((short *)(eobs))[0]) {
    if (((short *)eobs)[0] & 0xfefe)
      idct_dequant_full_2x_neon(q, dq, dst_u, stride);
    else
      idct_dequant_0_2x_neon(q, dq[0], dst_u, stride);
  }

  q += 32;
  dst_u += 4 * stride;

  if (((short *)(eobs))[1]) {
    if (((short *)eobs)[1] & 0xfefe)
      idct_dequant_full_2x_neon(q, dq, dst_u, stride);
    else
      idct_dequant_0_2x_neon(q, dq[0], dst_u, stride);
  }

  q += 32;

  if (((short *)(eobs))[2]) {
    if (((short *)eobs)[2] & 0xfefe)
      idct_dequant_full_2x_neon(q, dq, dst_v, stride);
    else
      idct_dequant_0_2x_neon(q, dq[0], dst_v, stride);
  }

  q += 32;
  dst_v += 4 * stride;

  if (((short *)(eobs))[3]) {
    if (((short *)eobs)[3] & 0xfefe)
      idct_dequant_full_2x_neon(q, dq, dst_v, stride);
    else
      idct_dequant_0_2x_neon(q, dq[0], dst_v, stride);
  }
}
