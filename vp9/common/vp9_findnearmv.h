/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */


#ifndef VP9_COMMON_VP9_FINDNEARMV_H_
#define VP9_COMMON_VP9_FINDNEARMV_H_

#include "vp9/common/vp9_mv.h"
#include "vp9/common/vp9_blockd.h"
#include "vp9/common/vp9_treecoder.h"
#include "vp9/common/vp9_onyxc_int.h"

#define LEFT_TOP_MARGIN     ((VP9BORDERINPIXELS - VP9_INTERP_EXTEND) << 3)
#define RIGHT_BOTTOM_MARGIN ((VP9BORDERINPIXELS - VP9_INTERP_EXTEND) << 3)

// check a list of motion vectors by sad score using a number rows of pixels
// above and a number cols of pixels in the left to select the one with best
// score to use as ref motion vector
void vp9_find_best_ref_mvs(MACROBLOCKD *xd,
                           int_mv *mvlist,
                           int_mv *nearest,
                           int_mv *near);

// TODO(jingning): this mv clamping function should be block size dependent.
static void clamp_mv2(MV *mv, const MACROBLOCKD *xd) {
  clamp_mv(mv, xd->mb_to_left_edge - LEFT_TOP_MARGIN,
               xd->mb_to_right_edge + RIGHT_BOTTOM_MARGIN,
               xd->mb_to_top_edge - LEFT_TOP_MARGIN,
               xd->mb_to_bottom_edge + RIGHT_BOTTOM_MARGIN);
}

void vp9_append_sub8x8_mvs_for_idx(VP9_COMMON *cm,
                                   MACROBLOCKD *xd,
                                   int_mv *dst_nearest,
                                   int_mv *dst_near,
                                   int block_idx, int ref_idx,
                                   int mi_row, int mi_col);

static MB_PREDICTION_MODE left_block_mode(const MODE_INFO *cur_mb,
                                          const MODE_INFO *left_mb, int b) {
  // FIXME(rbultje, jingning): temporary hack because jenkins doesn't
  // understand this condition. This will go away soon.
  const MODE_INFO *mi = cur_mb;

  if (b == 0 || b == 2) {
    /* On L edge, get from MB to left of us */
    mi = left_mb;
    if (!mi)
      return DC_PRED;

    if (mi->mbmi.ref_frame[0] != INTRA_FRAME) {
      return DC_PRED;
    } else if (mi->mbmi.sb_type < BLOCK_8X8) {
      return ((mi->bmi + 1 + b)->as_mode);
    } else {
      return mi->mbmi.mode;
    }
  }
  assert(b == 1 || b == 3);
  return (mi->bmi + b - 1)->as_mode;
}

static MB_PREDICTION_MODE above_block_mode(const MODE_INFO *cur_mb,
                                           const MODE_INFO *above_mb, int b) {
  const MODE_INFO *mi = cur_mb;

  if (!(b >> 1)) {
    /* On top edge, get from MB above us */
    mi = above_mb;
    if (!mi)
      return DC_PRED;

    if (mi->mbmi.ref_frame[0] != INTRA_FRAME) {
      return DC_PRED;
    } else if (mi->mbmi.sb_type < BLOCK_8X8) {
      return ((mi->bmi + 2 + b)->as_mode);
    } else {
      return mi->mbmi.mode;
    }
  }

  return (mi->bmi + b - 2)->as_mode;
}

#endif  // VP9_COMMON_VP9_FINDNEARMV_H_
