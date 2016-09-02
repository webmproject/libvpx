/*
 * Copyright (c) 2016, Alliance for Open Media. All rights reserved
 *
 * This source code is subject to the terms of the BSD 2 Clause License and
 * the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
 * was not distributed with this source code in the LICENSE file, you can
 * obtain it at www.aomedia.org/license/software. If the Alliance for Open
 * Media Patent License 1.0 was not distributed with this source code in the
 * PATENTS file, you can obtain it at www.aomedia.org/license/patent.
 */
#include "av1/common/clpf.h"

// Apply the filter on a single block
static void clpf_block(const uint8_t *src, uint8_t *dst, int sstride,
                       int dstride, int has_top, int has_left, int has_bottom,
                       int has_right, int width, int height) {
  int x, y;

  for (y = 0; y < height; y++) {
    for (x = 0; x < width; x++) {
      int X = src[(y + 0) * sstride + x + 0];
      int A = has_top ? src[(y - 1) * sstride + x + 0] : X;
      int B = has_left ? src[(y + 0) * sstride + x - 1] : X;
      int C = has_right ? src[(y + 0) * sstride + x + 1] : X;
      int D = has_bottom ? src[(y + 1) * sstride + x + 0] : X;
      int delta = ((A > X) + (B > X) + (C > X) + (D > X) > 2) -
                  ((A < X) + (B < X) + (C < X) + (D < X) > 2);
      dst[y * dstride + x] = X + delta;
    }
  }
}

#define BS (MI_SIZE * MAX_MIB_SIZE)

// Iterate over blocks within a superblock
static void av1_clpf_sb(const YV12_BUFFER_CONFIG *frame_buffer,
                        const AV1_COMMON *cm, MACROBLOCKD *xd,
                        MODE_INFO *const *mi_8x8, int xpos, int ypos) {
  // Temporary buffer (to allow SIMD parallelism)
  uint8_t buf_unaligned[BS * BS + 15];
  uint8_t *buf = (uint8_t *)(((intptr_t)buf_unaligned + 15) & ~15);
  int x, y, p;

  for (p = 0; p < (CLPF_FILTER_ALL_PLANES ? MAX_MB_PLANE : 1); p++) {
    for (y = 0; y < MAX_MIB_SIZE && ypos + y < cm->mi_rows; y++) {
      for (x = 0; x < MAX_MIB_SIZE && xpos + x < cm->mi_cols; x++) {
        const MB_MODE_INFO *mbmi =
            &mi_8x8[(ypos + y) * cm->mi_stride + xpos + x]->mbmi;

        // Do not filter if there is no residual
        if (!mbmi->skip) {
          // Do not filter frame edges
          int has_top = ypos + y > 0;
          int has_left = xpos + x > 0;
          int has_bottom = ypos + y < cm->mi_rows - 1;
          int has_right = xpos + x < cm->mi_cols - 1;
#if CLPF_ALLOW_BLOCK_PARALLELISM
          // Do not filter superblock edges
          has_top &= !!y;
          has_left &= !!x;
          has_bottom &= y != MAX_MIB_SIZE - 1;
          has_right &= x != MAX_MIB_SIZE - 1;
#endif
          av1_setup_dst_planes(xd->plane, frame_buffer, ypos + y, xpos + x);
          clpf_block(
              xd->plane[p].dst.buf, CLPF_ALLOW_PIXEL_PARALLELISM
                                        ? buf + y * MI_SIZE * BS + x * MI_SIZE
                                        : xd->plane[p].dst.buf,
              xd->plane[p].dst.stride,
              CLPF_ALLOW_PIXEL_PARALLELISM ? BS : xd->plane[p].dst.stride,
              has_top, has_left, has_bottom, has_right,
              MI_SIZE >> xd->plane[p].subsampling_x,
              MI_SIZE >> xd->plane[p].subsampling_y);
        }
      }
    }
#if CLPF_ALLOW_PIXEL_PARALLELISM
    for (y = 0; y < MAX_MIB_SIZE && ypos + y < cm->mi_rows; y++) {
      for (x = 0; x < MAX_MIB_SIZE && xpos + x < cm->mi_cols; x++) {
        const MB_MODE_INFO *mbmi =
            &mi_8x8[(ypos + y) * cm->mi_stride + xpos + x]->mbmi;
        av1_setup_dst_planes(xd->plane, frame_buffer, ypos + y, xpos + x);
        if (!mbmi->skip) {
          int i = 0;
          for (i = 0; i<MI_SIZE>> xd->plane[p].subsampling_y; i++)
            memcpy(xd->plane[p].dst.buf + i * xd->plane[p].dst.stride,
                   buf + (y * MI_SIZE + i) * BS + x * MI_SIZE,
                   MI_SIZE >> xd->plane[p].subsampling_x);
        }
      }
    }
#endif
  }
}

// Iterate over the superblocks of an entire frame
void av1_clpf_frame(const YV12_BUFFER_CONFIG *frame, const AV1_COMMON *cm,
                    MACROBLOCKD *xd) {
  int x, y;

  for (y = 0; y < cm->mi_rows; y += MAX_MIB_SIZE)
    for (x = 0; x < cm->mi_cols; x += MAX_MIB_SIZE)
      av1_clpf_sb(frame, cm, xd, cm->mi_grid_visible, x, y);
}
