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
#include "ec_types.h"

/* Allocate memory for the overlap lists */
int vp8_alloc_overlap_lists(VP8D_COMP *pbi);

/* Deallocate the overlap lists */
void vp8_de_alloc_overlap_lists(VP8D_COMP *pbi);

/* Inserts a new overlap area value to the list of overlaps of a block */
void vp8_assign_overlap(OVERLAP_NODE* overlaps,
                        B_MODE_INFO *bmi,
                        MV_REFERENCE_FRAME ref_frame,
                        int overlap);

/* Calculates the overlap area between two 4x4 squares, where the first
 * square has its upper-left corner at (b1_row, b1_col) and the second
 * square has its upper-left corner at (b2_row, b2_col). Doesn't
 * properly handle squares which doesn't overlap.
 */
int vp8_block_overlap(int b1_row, int b1_col, int b2_row, int b2_col);

/* Finds the reference frame type which has the largest overlapping area. */
MV_REFERENCE_FRAME vp8_largest_overlap_type(const B_OVERLAP *block_overlaps);

/* Calculates the overlap area for all blocks in a macroblock at position
 * (mb_row, mb_col) in macroblocks, which are being overlapped by a given
 * overlapping block at position (new_row, new_col) (in pixels, Q3). The
 * first block being overlapped in the macroblock has position (first_blk_row,
 * first_blk_col) in blocks relative the upper-left corner of the image.
 */
void vp8_calculate_overlaps_mb(B_OVERLAP *b_overlaps, B_MODE_INFO *bmi,
                               MV_REFERENCE_FRAME ref_frame,
                               int new_row, int new_col,
                               int mb_row, int mb_col,
                               int first_blk_row, int first_blk_col);

/* Estimates a motion vector given the overlapping blocks' motion vectors.
 * Filters out all overlapping blocks which doesn't refer to the correct
 * reference frame type.
 */
void vp8_estimate_mv(const OVERLAP_NODE *overlaps, B_MODE_INFO *bmi,
                     MV_REFERENCE_FRAME type);

/* Estimates all motion vectors for a macroblock given the lists of
 * overlaps for each block. Decides whether or not the MVs must be clamped.
 */
void vp8_estimate_mb_mvs(const B_OVERLAP *block_overlaps,
                         MODE_INFO *mi,
                         int mb_to_left_edge,
                         int mb_to_right_edge,
                         int mb_to_top_edge,
                         int mb_to_bottom_edge);

/* Estimate all missing motion vectors.
 */
void vp8_estimate_missing_mvs(VP8D_COMP *pbi);

/* Estimate all missing motion vectors */
void vp8_estimate_missing_mvs_ex(MB_OVERLAP *overlaps,
                                 MODE_INFO *mi, MODE_INFO *prev_mi,
                                 int mb_rows, int mb_cols,
                                 unsigned int first_corrupt);

/* Functions for spatial MV interpolation */

/* Finds the neighboring blocks of a macroblocks. In the general case
 * 20 blocks are found. If a fewer number of blocks are found due to
 * image boundaries, those positions in the EC_BLOCK array are left "empty".
 * The neighbors are enumerated with the upper-left neighbor as the first
 * element, the second element refers to the neighbor to right of the previous
 * neighbor, and so on. The last element refers to the neighbor below the first
 * neighbor.
 */
void vp8_find_neighboring_blocks(MODE_INFO *mi,
                                 EC_BLOCK *neighbors,
                                 int mb_row, int mb_col,
                                 int mb_rows, int mb_cols,
                                 int mi_stride);

/* Calculates which reference frame type is dominating among the neighbors */
MV_REFERENCE_FRAME vp8_dominant_ref_frame(EC_BLOCK *neighbors);

/* Interpolates all motion vectors for a macroblock from the neighboring blocks'
 * motion vectors.
 */
void vp8_interpolate_mvs(MACROBLOCKD *mb,
                         EC_BLOCK *neighbors,
                         MV_REFERENCE_FRAME dom_ref_frame);

/* Interpolates all motion vectors for a macroblock mb at position
 * (mb_row, mb_col). */
void vp8_interpolate_motion(MACROBLOCKD *mb,
                            int mb_row, int mb_col,
                            int mb_rows, int mb_cols,
                            int mi_stride);

#endif
