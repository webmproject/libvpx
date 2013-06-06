/*
 *  Copyright (c) 2012 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#ifndef VP9_COMMON_VP9_PRED_COMMON_H_
#define VP9_COMMON_VP9_PRED_COMMON_H_

#include "vp9/common/vp9_blockd.h"
#include "vp9/common/vp9_onyxc_int.h"

// Predicted items
typedef enum {
  PRED_SEG_ID = 0,  // Segment identifier
  PRED_MBSKIP = 1,
  PRED_SWITCHABLE_INTERP = 2,
  PRED_INTRA_INTER = 3,
  PRED_COMP_INTER_INTER = 4,
  PRED_SINGLE_REF_P1 = 5,
  PRED_SINGLE_REF_P2 = 6,
  PRED_COMP_REF_P = 7,
  PRED_TX_SIZE = 8
} PRED_ID;

unsigned char vp9_get_pred_context(const VP9_COMMON *const cm,
                                   const MACROBLOCKD *const xd,
                                   PRED_ID pred_id);

vp9_prob vp9_get_pred_prob(const VP9_COMMON *const cm,
                           const MACROBLOCKD *const xd,
                           PRED_ID pred_id);

const vp9_prob *vp9_get_pred_probs(const VP9_COMMON *const cm,
                                   const MACROBLOCKD *const xd,
                                   PRED_ID pred_id);

unsigned char vp9_get_pred_flag(const MACROBLOCKD *const xd,
                                PRED_ID pred_id);

void vp9_set_pred_flag(MACROBLOCKD *const xd,
                       PRED_ID pred_id,
                       unsigned char pred_flag);


int vp9_get_pred_mi_segid(VP9_COMMON *cm, BLOCK_SIZE_TYPE sb_type,
                          int mi_row, int mi_col);

#endif  // VP9_COMMON_VP9_PRED_COMMON_H_
