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
#include "vp9/common/vp9_blockd.h"

static INLINE void recon(int rows, int cols,
                         const int16_t *diff_ptr, int diff_stride,
                         uint8_t *dst_ptr, int dst_stride) {
  int r, c;

  for (r = 0; r < rows; r++) {
    for (c = 0; c < cols; c++)
      dst_ptr[c] = clip_pixel(diff_ptr[c] + dst_ptr[c]);

    dst_ptr += dst_stride;
    diff_ptr += diff_stride;
  }
}


void vp9_recon_b_c(uint8_t *pred_ptr, int16_t *diff_ptr, int diff_stride,
                   uint8_t *dst_ptr, int stride) {
  assert(pred_ptr == dst_ptr);
  recon(4, 4, diff_ptr, diff_stride, dst_ptr, stride);
}

#if !CONFIG_SB8X8
void vp9_recon_uv_b_c(uint8_t *pred_ptr, int16_t *diff_ptr, uint8_t *dst_ptr,
                      int stride) {
  assert(pred_ptr == dst_ptr);
  recon(4, 4, diff_ptr, 8, dst_ptr, stride);
}

void vp9_recon4b_c(uint8_t *pred_ptr, int16_t *diff_ptr, uint8_t *dst_ptr,
                   int stride) {
  assert(pred_ptr == dst_ptr);
  recon(4, 16, diff_ptr, 16, dst_ptr, stride);
}

void vp9_recon2b_c(uint8_t *pred_ptr, int16_t *diff_ptr, uint8_t *dst_ptr,
                   int stride) {
  assert(pred_ptr == dst_ptr);
  recon(4, 8, diff_ptr, 8, dst_ptr, stride);
}
#endif

static void recon_plane(MACROBLOCKD *xd, BLOCK_SIZE_TYPE bsize, int plane) {
  const int bw = 4 << (b_width_log2(bsize) - xd->plane[plane].subsampling_x);
  const int bh = 4 << (b_height_log2(bsize) - xd->plane[plane].subsampling_y);
  recon(bh, bw,
        xd->plane[plane].diff, bw,
        xd->plane[plane].dst.buf, xd->plane[plane].dst.stride);
}

void vp9_recon_sby_c(MACROBLOCKD *mb, BLOCK_SIZE_TYPE bsize) {
  recon_plane(mb, bsize, 0);
}

void vp9_recon_sbuv_c(MACROBLOCKD *mb, BLOCK_SIZE_TYPE bsize) {
  int i;

  for (i = 1; i < MAX_MB_PLANE; i++)
    recon_plane(mb, bsize, i);
}

void vp9_recon_sb_c(MACROBLOCKD *xd, BLOCK_SIZE_TYPE bsize) {
  vp9_recon_sby(xd, bsize);
  vp9_recon_sbuv(xd, bsize);
}
