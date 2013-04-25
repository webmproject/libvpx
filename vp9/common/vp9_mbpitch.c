/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */


#include "vp9/common/vp9_blockd.h"

typedef enum {
  PRED = 0,
  DEST = 1
} BLOCKSET;

static void setup_block(BLOCKD *b, uint8_t **base, uint8_t **base2,
                        int stride, int offset, BLOCKSET bs) {
  if (bs == DEST) {
    b->dst_stride = stride;
    b->dst = offset;
    b->base_dst = base;
  }
}

static void setup_macroblock(MACROBLOCKD *mb, BLOCKSET bs) {
  BLOCKD *blockd = mb->block;
  uint8_t **y, **u, **v, **y2, **u2, **v2;
  int i, stride;

  if (bs == DEST) {
    y = &mb->plane[0].dst.buf;
    u = &mb->plane[1].dst.buf;
    v = &mb->plane[2].dst.buf;

    y2 = NULL;
    u2 = NULL;
    v2 = NULL;
  }

  // luma
  stride = mb->plane[0].dst.stride;
  for (i = 0; i < 16; ++i) {
    const int offset = (i >> 2) * 4 * stride + (i & 3) * 4;
    setup_block(&blockd[i], y, y2, stride, offset, bs);
  }

  // chroma
  stride = mb->plane[1].dst.stride;
  for (i = 16; i < 20; i++) {
    const int offset = ((i - 16) >> 1) * 4 * stride + (i & 1) * 4;
    setup_block(&blockd[i],     u, u2, stride, offset, bs);
    setup_block(&blockd[i + 4], v, v2, stride, offset, bs);
  }
}

void vp9_setup_block_dptrs(MACROBLOCKD *mb) {
  int i;

  for (i = 0; i < MAX_MB_PLANE; i++) {
    mb->plane[i].plane_type = i ? PLANE_TYPE_UV : PLANE_TYPE_Y_WITH_DC;
    mb->plane[i].subsampling_x = !!i;
    mb->plane[i].subsampling_y = !!i;
  }
}

void vp9_build_block_doffsets(MACROBLOCKD *mb) {
  // handle the destination pitch features
  setup_macroblock(mb, DEST);
}
