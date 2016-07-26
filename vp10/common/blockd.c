/*
 *  Copyright (c) 2014 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include <math.h>

#include "vpx_ports/system_state.h"

#include "vp10/common/blockd.h"

PREDICTION_MODE vp10_left_block_mode(const MODE_INFO *cur_mi,
                                    const MODE_INFO *left_mi, int b) {
  if (b == 0 || b == 2) {
    if (!left_mi || is_inter_block(&left_mi->mbmi))
      return DC_PRED;

    return get_y_mode(left_mi, b + 1);
  } else {
    assert(b == 1 || b == 3);
    return cur_mi->bmi[b - 1].as_mode;
  }
}

PREDICTION_MODE vp10_above_block_mode(const MODE_INFO *cur_mi,
                                     const MODE_INFO *above_mi, int b) {
  if (b == 0 || b == 1) {
    if (!above_mi || is_inter_block(&above_mi->mbmi))
      return DC_PRED;

    return get_y_mode(above_mi, b + 2);
  } else {
    assert(b == 2 || b == 3);
    return cur_mi->bmi[b - 2].as_mode;
  }
}

void vp10_foreach_transformed_block_in_plane(
    const MACROBLOCKD *const xd, BLOCK_SIZE bsize, int plane,
    foreach_transformed_block_visitor visit, void *arg) {
  const struct macroblockd_plane *const pd = &xd->plane[plane];
  const MB_MODE_INFO* mbmi = &xd->mi[0]->mbmi;
  // block and transform sizes, in number of 4x4 blocks log 2 ("*_b")
  // 4x4=0, 8x8=2, 16x16=4, 32x32=6, 64x64=8
  // transform size varies per plane, look it up in a common way.
  const TX_SIZE tx_size = plane ? get_uv_tx_size(mbmi, pd)
                                : mbmi->tx_size;
  const BLOCK_SIZE plane_bsize = get_plane_block_size(bsize, pd);
  const int num_4x4_w = num_4x4_blocks_wide_lookup[plane_bsize];
  const int num_4x4_h = num_4x4_blocks_high_lookup[plane_bsize];
  const uint8_t num_4x4_tw = num_4x4_blocks_wide_txsize_lookup[tx_size];
  const uint8_t num_4x4_th = num_4x4_blocks_high_txsize_lookup[tx_size];
  const int step = num_4x4_tw * num_4x4_th;
  int i = 0, r, c;

  // If mb_to_right_edge is < 0 we are in a situation in which
  // the current block size extends into the UMV and we won't
  // visit the sub blocks that are wholly within the UMV.
  const int max_blocks_wide = num_4x4_w + (xd->mb_to_right_edge >= 0 ? 0 :
      xd->mb_to_right_edge >> (5 + pd->subsampling_x));
  const int max_blocks_high = num_4x4_h + (xd->mb_to_bottom_edge >= 0 ? 0 :
      xd->mb_to_bottom_edge >> (5 + pd->subsampling_y));
  const int extra_step =
      ((num_4x4_w - max_blocks_wide) >>
       num_4x4_blocks_wide_txsize_log2_lookup[tx_size]) * step;

  // Keep track of the row and column of the blocks we use so that we know
  // if we are in the unrestricted motion border.
  for (r = 0; r < max_blocks_high; r += num_4x4_th) {
    // Skip visiting the sub blocks that are wholly within the UMV.
    for (c = 0; c < max_blocks_wide; c += num_4x4_tw) {
      visit(plane, i, r, c, plane_bsize, tx_size, arg);
      i += step;
    }
    i += extra_step;
  }
}

void vp10_foreach_transformed_block(const MACROBLOCKD* const xd,
                                   BLOCK_SIZE bsize,
                                   foreach_transformed_block_visitor visit,
                                   void *arg) {
  int plane;
  for (plane = 0; plane < MAX_MB_PLANE; ++plane)
    vp10_foreach_transformed_block_in_plane(xd, bsize, plane, visit, arg);
}

void vp10_set_contexts(const MACROBLOCKD *xd, struct macroblockd_plane *pd,
                       BLOCK_SIZE plane_bsize, TX_SIZE tx_size, int has_eob,
                       int aoff, int loff) {
  ENTROPY_CONTEXT *const a = pd->above_context + aoff;
  ENTROPY_CONTEXT *const l = pd->left_context + loff;
  const int tx_w_in_blocks = num_4x4_blocks_wide_txsize_lookup[tx_size];
  const int tx_h_in_blocks = num_4x4_blocks_high_txsize_lookup[tx_size];

  // above
  if (has_eob && xd->mb_to_right_edge < 0) {
    int i;
    const int blocks_wide = num_4x4_blocks_wide_lookup[plane_bsize] +
                            (xd->mb_to_right_edge >> (5 + pd->subsampling_x));
    int above_contexts = tx_w_in_blocks;
    if (above_contexts + aoff > blocks_wide)
      above_contexts = blocks_wide - aoff;

    for (i = 0; i < above_contexts; ++i)
      a[i] = has_eob;
    for (i = above_contexts; i < tx_w_in_blocks; ++i)
      a[i] = 0;
  } else {
    memset(a, has_eob, sizeof(ENTROPY_CONTEXT) * tx_w_in_blocks);
  }

  // left
  if (has_eob && xd->mb_to_bottom_edge < 0) {
    int i;
    const int blocks_high = num_4x4_blocks_high_lookup[plane_bsize] +
                            (xd->mb_to_bottom_edge >> (5 + pd->subsampling_y));
    int left_contexts = tx_h_in_blocks;
    if (left_contexts + loff > blocks_high)
      left_contexts = blocks_high - loff;

    for (i = 0; i < left_contexts; ++i)
      l[i] = has_eob;
    for (i = left_contexts; i < tx_h_in_blocks; ++i)
      l[i] = 0;
  } else {
    memset(l, has_eob, sizeof(ENTROPY_CONTEXT) * tx_h_in_blocks);
  }
}

void vp10_setup_block_planes(MACROBLOCKD *xd, int ss_x, int ss_y) {
  int i;

  for (i = 0; i < MAX_MB_PLANE; i++) {
    xd->plane[i].plane_type = i ? PLANE_TYPE_UV : PLANE_TYPE_Y;
    xd->plane[i].subsampling_x = i ? ss_x : 0;
    xd->plane[i].subsampling_y = i ? ss_y : 0;
  }
}

#if CONFIG_EXT_INTRA
const int16_t dr_intra_derivative[90] = {
    1,     14666,  7330,  4884,  3660,
    2926,  2435,   2084,  1821,  1616,
    1451,  1317,   1204,  1108,  1026,
    955,   892,    837,   787,   743,
    703,   666,    633,   603,   574,
    548,   524,    502,   481,   461,
    443,   426,    409,   394,   379,
    365,   352,    339,   327,   316,
    305,   294,    284,   274,   265,
    256,   247,    238,   230,   222,
    214,   207,    200,   192,   185,
    179,   172,    166,   159,   153,
    147,   141,    136,   130,   124,
    119,   113,    108,   103,   98,
    93,    88,     83,    78,    73,
    68,    63,     59,    54,    49,
    45,    40,     35,    31,    26,
    22,    17,     13,    8,     4,
};

// Returns whether filter selection is needed for a given
// intra prediction angle.
int vp10_is_intra_filter_switchable(int angle) {
  assert(angle > 0 && angle < 270);
  if (angle % 45 == 0)
    return 0;
  if (angle > 90 && angle < 180) {
    return 1;
  } else {
    return ((angle < 90 ? dr_intra_derivative[angle] :
        dr_intra_derivative[270 - angle]) & 0xFF) > 0;
  }
}
#endif  // CONFIG_EXT_INTRA
