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
// If angle > 0 && angle < 90, dx = -((int)(256 / t)), dy = 1;
// If angle > 90 && angle < 180, dx = (int)(256 / t), dy = (int)(256 * t);
// If angle > 180 && angle < 270, dx = 1, dy = -((int)(256 * t));
const int16_t dr_intra_derivative[270][2] = {
    {     1,     1 }, { -14666,    1 }, { -7330,     1 }, { -4884,     1 },
    { -3660,     1 }, { -2926,     1 }, { -2435,     1 }, { -2084,     1 },
    { -1821,     1 }, { -1616,     1 }, { -1451,     1 }, { -1317,     1 },
    { -1204,     1 }, { -1108,     1 }, { -1026,     1 }, {  -955,     1 },
    {  -892,     1 }, {  -837,     1 }, {  -787,     1 }, {  -743,     1 },
    {  -703,     1 }, {  -666,     1 }, {  -633,     1 }, {  -603,     1 },
    {  -574,     1 }, {  -548,     1 }, {  -524,     1 }, {  -502,     1 },
    {  -481,     1 }, {  -461,     1 }, {  -443,     1 }, {  -426,     1 },
    {  -409,     1 }, {  -394,     1 }, {  -379,     1 }, {  -365,     1 },
    {  -352,     1 }, {  -339,     1 }, {  -327,     1 }, {  -316,     1 },
    {  -305,     1 }, {  -294,     1 }, {  -284,     1 }, {  -274,     1 },
    {  -265,     1 }, {  -256,     1 }, {  -247,     1 }, {  -238,     1 },
    {  -230,     1 }, {  -222,     1 }, {  -214,     1 }, {  -207,     1 },
    {  -200,     1 }, {  -192,     1 }, {  -185,     1 }, {  -179,     1 },
    {  -172,     1 }, {  -166,     1 }, {  -159,     1 }, {  -153,     1 },
    {  -147,     1 }, {  -141,     1 }, {  -136,     1 }, {  -130,     1 },
    {  -124,     1 }, {  -119,     1 }, {  -113,     1 }, {  -108,     1 },
    {  -103,     1 }, {   -98,     1 }, {   -93,     1 }, {   -88,     1 },
    {   -83,     1 }, {   -78,     1 }, {   -73,     1 }, {   -68,     1 },
    {   -63,     1 }, {   -59,     1 }, {   -54,     1 }, {   -49,     1 },
    {   -45,     1 }, {   -40,     1 }, {   -35,     1 }, {   -31,     1 },
    {   -26,     1 }, {   -22,     1 }, {   -17,     1 }, {   -13,     1 },
    {    -8,     1 }, {    -4,     1 }, {     1,     1 }, {     4, 14666 },
    {     8,  7330 }, {    13,  4884 }, {    17,  3660 }, {    22,  2926 },
    {    26,  2435 }, {    31,  2084 }, {    35,  1821 }, {    40,  1616 },
    {    45,  1451 }, {    49,  1317 }, {    54,  1204 }, {    59,  1108 },
    {    63,  1026 }, {    68,   955 }, {    73,   892 }, {    78,   837 },
    {    83,   787 }, {    88,   743 }, {    93,   703 }, {    98,   666 },
    {   103,   633 }, {   108,   603 }, {   113,   574 }, {   119,   548 },
    {   124,   524 }, {   130,   502 }, {   136,   481 }, {   141,   461 },
    {   147,   443 }, {   153,   426 }, {   159,   409 }, {   166,   394 },
    {   172,   379 }, {   179,   365 }, {   185,   352 }, {   192,   339 },
    {   200,   327 }, {   207,   316 }, {   214,   305 }, {   222,   294 },
    {   230,   284 }, {   238,   274 }, {   247,   265 }, {   255,   256 },
    {   265,   247 }, {   274,   238 }, {   284,   230 }, {   294,   222 },
    {   305,   214 }, {   316,   207 }, {   327,   200 }, {   339,   192 },
    {   352,   185 }, {   365,   179 }, {   379,   172 }, {   394,   166 },
    {   409,   159 }, {   426,   153 }, {   443,   147 }, {   461,   141 },
    {   481,   136 }, {   502,   130 }, {   524,   124 }, {   548,   119 },
    {   574,   113 }, {   603,   108 }, {   633,   103 }, {   666,    98 },
    {   703,    93 }, {   743,    88 }, {   787,    83 }, {   837,    78 },
    {   892,    73 }, {   955,    68 }, {  1026,    63 }, {  1108,    59 },
    {  1204,    54 }, {  1317,    49 }, {  1451,    45 }, {  1616,    40 },
    {  1821,    35 }, {  2084,    31 }, {  2435,    26 }, {  2926,    22 },
    {  3660,    17 }, {  4884,    13 }, {  7330,     8 }, { 14666,     4 },
    {     1,     1 }, {     1,    -4 }, {     1,    -8 }, {     1,   -13 },
    {     1,   -17 }, {     1,   -22 }, {     1,   -26 }, {     1,   -31 },
    {     1,   -35 }, {     1,   -40 }, {     1,   -45 }, {     1,   -49 },
    {     1,   -54 }, {     1,   -59 }, {     1,   -63 }, {     1,   -68 },
    {     1,   -73 }, {     1,   -78 }, {     1,   -83 }, {     1,   -88 },
    {     1,   -93 }, {     1,   -98 }, {     1,  -103 }, {     1,  -108 },
    {     1,  -113 }, {     1,  -119 }, {     1,  -124 }, {     1,  -130 },
    {     1,  -136 }, {     1,  -141 }, {     1,  -147 }, {     1,  -153 },
    {     1,  -159 }, {     1,  -166 }, {     1,  -172 }, {     1,  -179 },
    {     1,  -185 }, {     1,  -192 }, {     1,  -200 }, {     1,  -207 },
    {     1,  -214 }, {     1,  -222 }, {     1,  -230 }, {     1,  -238 },
    {     1,  -247 }, {     1,  -255 }, {     1,  -265 }, {     1,  -274 },
    {     1,  -284 }, {     1,  -294 }, {     1,  -305 }, {     1,  -316 },
    {     1,  -327 }, {     1,  -339 }, {     1,  -352 }, {     1,  -365 },
    {     1,  -379 }, {     1,  -394 }, {     1,  -409 }, {     1,  -426 },
    {     1,  -443 }, {     1,  -461 }, {     1,  -481 }, {     1,  -502 },
    {     1,  -524 }, {     1,  -548 }, {     1,  -574 }, {     1,  -603 },
    {     1,  -633 }, {     1,  -666 }, {     1,  -703 }, {     1,  -743 },
    {     1,  -787 }, {     1,  -837 }, {     1,  -892 }, {     1,  -955 },
    {     1, -1026 }, {     1, -1108 }, {     1, -1204 }, {     1, -1317 },
    {     1, -1451 }, {     1, -1616 }, {     1, -1821 }, {     1, -2084 },
    {     1, -2435 }, {     1, -2926 }, {     1, -3660 }, {     1, -4884 },
    {     1, -7330 }, {     1, -14666 },
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
    return ((-(dr_intra_derivative[angle][angle > 180])) & 0xFF) > 0;
  }
}
#endif  // CONFIG_EXT_INTRA
