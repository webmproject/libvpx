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
#include "vp9_rtcd.h"
#include "vp9/common/vp9_blockd.h"

void vp9_recon_b_c
(
  unsigned char *pred_ptr,
  short *diff_ptr,
  unsigned char *dst_ptr,
  int stride
) {
  int r, c;

  for (r = 0; r < 4; r++) {
    for (c = 0; c < 4; c++) {
      int a = diff_ptr[c] + pred_ptr[c];

      if (a < 0)
        a = 0;

      if (a > 255)
        a = 255;

      dst_ptr[c] = (unsigned char) a;
    }

    dst_ptr += stride;
    diff_ptr += 16;
    pred_ptr += 16;
  }
}

void vp9_recon_uv_b_c
(
  unsigned char *pred_ptr,
  short *diff_ptr,
  unsigned char *dst_ptr,
  int stride
) {
  int r, c;

  for (r = 0; r < 4; r++) {
    for (c = 0; c < 4; c++) {
      int a = diff_ptr[c] + pred_ptr[c];

      if (a < 0)
        a = 0;

      if (a > 255)
        a = 255;

      dst_ptr[c] = (unsigned char) a;
    }

    dst_ptr += stride;
    diff_ptr += 8;
    pred_ptr += 8;
  }
}
void vp9_recon4b_c
(
  unsigned char *pred_ptr,
  short *diff_ptr,
  unsigned char *dst_ptr,
  int stride
) {
  int r, c;

  for (r = 0; r < 4; r++) {
    for (c = 0; c < 16; c++) {
      int a = diff_ptr[c] + pred_ptr[c];

      if (a < 0)
        a = 0;

      if (a > 255)
        a = 255;

      dst_ptr[c] = (unsigned char) a;
    }

    dst_ptr += stride;
    diff_ptr += 16;
    pred_ptr += 16;
  }
}

void vp9_recon2b_c
(
  unsigned char *pred_ptr,
  short *diff_ptr,
  unsigned char *dst_ptr,
  int stride
) {
  int r, c;

  for (r = 0; r < 4; r++) {
    for (c = 0; c < 8; c++) {
      int a = diff_ptr[c] + pred_ptr[c];

      if (a < 0)
        a = 0;

      if (a > 255)
        a = 255;

      dst_ptr[c] = (unsigned char) a;
    }

    dst_ptr += stride;
    diff_ptr += 8;
    pred_ptr += 8;
  }
}

#if CONFIG_SUPERBLOCKS
void vp9_recon_mby_s_c(MACROBLOCKD *xd, uint8_t *dst) {
  int x, y;
  BLOCKD *b = &xd->block[0];
  int stride = b->dst_stride;
  short *diff = b->diff;

  for (y = 0; y < 16; y++) {
    for (x = 0; x < 16; x++) {
      int a = dst[x] + diff[x];
      if (a < 0)
        a = 0;
      else if (a > 255)
        a = 255;
      dst[x] = a;
    }
    dst += stride;
    diff += 16;
  }
}

void vp9_recon_mbuv_s_c(MACROBLOCKD *xd, uint8_t *udst, uint8_t *vdst) {
  int x, y, i;
  uint8_t *dst = udst;

  for (i = 0; i < 2; i++, dst = vdst) {
    BLOCKD *b = &xd->block[16 + 4 * i];
    int stride = b->dst_stride;
    short *diff = b->diff;

    for (y = 0; y < 8; y++) {
      for (x = 0; x < 8; x++) {
        int a = dst[x] + diff[x];
        if (a < 0)
          a = 0;
        else if (a > 255)
          a = 255;
        dst[x] = a;
      }
      dst += stride;
      diff += 8;
    }
  }
}
#endif

void vp9_recon_mby_c(MACROBLOCKD *xd) {
  int i;

  for (i = 0; i < 16; i += 4) {
    BLOCKD *b = &xd->block[i];

    vp9_recon4b(b->predictor, b->diff, *(b->base_dst) + b->dst, b->dst_stride);
  }
}

void vp9_recon_mb_c(MACROBLOCKD *xd) {
  int i;

  for (i = 0; i < 16; i += 4) {
    BLOCKD *b = &xd->block[i];

    vp9_recon4b(b->predictor, b->diff, *(b->base_dst) + b->dst, b->dst_stride);
  }

  for (i = 16; i < 24; i += 2) {
    BLOCKD *b = &xd->block[i];

    vp9_recon2b(b->predictor, b->diff, *(b->base_dst) + b->dst, b->dst_stride);
  }
}
