/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include <limits.h>

#include "vp9/common/vp9_findnearmv.h"
#include "vp9/common/vp9_sadmxn.h"
#include "vp9/common/vp9_subpelvar.h"

const uint8_t vp9_mbsplit_offset[4][16] = {
  { 0,  8,  0,  0,  0,  0,  0,  0,  0,  0,   0,  0,  0,  0,  0,  0},
  { 0,  2,  0,  0,  0,  0,  0,  0,  0,  0,   0,  0,  0,  0,  0,  0},
  { 0,  2,  8, 10,  0,  0,  0,  0,  0,  0,   0,  0,  0,  0,  0,  0},
  { 0,  1,  2,  3,  4,  5,  6,  7,  8,  9,  10, 11, 12, 13, 14, 15}
};

static void lower_mv_precision(int_mv *mv, int usehp) {
  if (!usehp || !vp9_use_nmv_hp(&mv->as_mv)) {
    if (mv->as_mv.row & 1)
      mv->as_mv.row += (mv->as_mv.row > 0 ? -1 : 1);
    if (mv->as_mv.col & 1)
      mv->as_mv.col += (mv->as_mv.col > 0 ? -1 : 1);
  }
}

vp9_prob *vp9_mv_ref_probs(VP9_COMMON *pc, vp9_prob p[4], int context) {
  p[0] = pc->fc.vp9_mode_contexts[context][0];
  p[1] = pc->fc.vp9_mode_contexts[context][1];
  p[2] = pc->fc.vp9_mode_contexts[context][2];
  p[3] = pc->fc.vp9_mode_contexts[context][3];
  return p;
}

void vp9_find_best_ref_mvs(MACROBLOCKD *xd,
                           int_mv *mvlist,
                           int_mv *nearest,
                           int_mv *near) {
  int i;
  // Make sure all the candidates are properly clamped etc
  for (i = 0; i < MAX_MV_REF_CANDIDATES; ++i) {
    lower_mv_precision(&mvlist[i], xd->allow_high_precision_mv);
    clamp_mv2(&mvlist[i], xd);
  }
  *nearest = mvlist[0];
  *near = mvlist[1];
}
