/*
 *  Copyright (c) 2011 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */


#ifndef ERROR_CONCEALMENT_H
#define ERROR_CONCEALMENT_H

#include "onyxd_int.h"

int vp8_alloc_overlap_lists(VP8D_COMP *pbi);
void vp8_de_alloc_overlap_lists(VP8D_COMP *pbi);

void vp8_assign_overlap(OVERLAP_NODE* overlaps,
                        B_MODE_INFO *bmi,
                        MV_REFERENCE_FRAME ref_frame,
                        int overlap);
int vp8_block_overlap(int b1_row, int b1_col, int b2_row, int b2_col);
void vp8_calculate_overlaps_mb(B_OVERLAP *b_overlaps, B_MODE_INFO *bmi,
                               MV_REFERENCE_FRAME ref_frame,
                               int new_row, int new_col,
                               int first_ol_mb_row, int first_ol_mb_col,
                               int first_ol_blk_row, int first_ol_blk_col);
MV_REFERENCE_FRAME vp8_largest_overlap_type(const B_OVERLAP *block_overlaps);
void vp8_estimate_mv(const OVERLAP_NODE *overlaps, B_MODE_INFO *bmi,
                     MV_REFERENCE_FRAME type);
void vp8_estimate_mb_mvs(const B_OVERLAP *block_overlaps,
                         B_MODE_INFO *bmi,
                         MV_REFERENCE_FRAME type,
                         MV* filtered_mv);
void vp8_estimate_missing_mvs(VP8D_COMP *pbi);
void vp8_estimate_missing_mvs_ex(MB_OVERLAP *overlaps,
                                 MODE_INFO *mi, MODE_INFO *prev_mi,
                                 int mb_rows, int mb_cols,
                                 unsigned int first_corrupt);
void vp8_conceal_corrupt_block(MACROBLOCKD *);

#endif
