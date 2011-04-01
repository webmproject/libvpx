/*
 *  Copyright (c) 2011 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include "error_concealment.h"
#include "onyxd_int.h"
#include "vpx_mem/vpx_mem.h"

#include <assert.h>

#define MIN(x,y) (((x)<(y))?(x):(y))
#define MAX(x,y) (((x)>(y))?(x):(y))

#define FLOOR(x,q) ((x) & -(1 << (q)))

int vp8_alloc_overlap_lists(VP8D_COMP *pbi)
{
    if (pbi->overlaps != NULL)
    {
        vpx_free(pbi->overlaps);
        pbi->overlaps = NULL;
    }
    pbi->overlaps = vpx_calloc(pbi->common.mb_rows * pbi->common.mb_cols,
                               sizeof(MB_OVERLAP));
    vpx_memset(pbi->overlaps, 0,
               sizeof(MB_OVERLAP) * pbi->common.mb_rows * pbi->common.mb_cols);
    if (pbi->overlaps == NULL)
        return -1;
    return 0;
}

void vp8_de_alloc_overlap_lists(VP8D_COMP *pbi)
{
    if (pbi->overlaps != NULL)
    {
        vpx_free(pbi->overlaps);
        pbi->overlaps = NULL;
    }
}

void vp8_assign_overlap(OVERLAP_NODE* overlaps,
                        B_MODE_INFO *bmi,
                        MV_REFERENCE_FRAME ref_frame,
                        int overlap)
{
    int i;
    if (overlap <= 0)
        return;
    for (i = 0; i < MAX_OVERLAPS; i++)
    {
        if (overlaps[i].bmi == NULL)
        {
            overlaps[i].bmi = bmi;
            overlaps[i].ref_frame = ref_frame;
            overlaps[i].overlap = overlap;
            break;
        }
    }
}

int vp8_block_overlap(int b1_row, int b1_col, int b2_row, int b2_col)
{
    const int int_top = MAX(b1_row, b2_row); // top
    const int int_left = MAX(b1_col, b2_col); // left
    const int int_right = MIN(b1_col + (4<<3), b2_col + (4<<3)); // right
    const int int_bottom = MIN(b1_row + (4<<3), b2_row + (4<<3)); // bottom
    return (int_bottom - int_top) * (int_right - int_left);
}

void vp8_calculate_overlaps_mb(B_OVERLAP *b_overlaps, B_MODE_INFO *bmi,
                               MV_REFERENCE_FRAME ref_frame,
                               int new_row, int new_col,
                               int first_ol_mb_row, int first_ol_mb_col,
                               int first_ol_blk_row, int first_ol_blk_col)
{
    /* find the blocks it's overlapping */
    const int rel_ol_blk_row = first_ol_blk_row - first_ol_mb_row * 4;
    const int rel_ol_blk_col = first_ol_blk_col - first_ol_mb_col * 4;
    const int blk_idx = MAX(rel_ol_blk_row,0) * 4 + MAX(rel_ol_blk_col,0);
    /* Upper left overlapping block */
    B_OVERLAP *b_ol_ul = &(b_overlaps[blk_idx]);

    /* Calculate and assign overlaps for all blocks in this MB
     * which the motion compensated block overlaps
     */
    int row, col;
    int end_row = MIN(4 + first_ol_mb_row * 4 - first_ol_blk_row, 2);
    int end_col = MIN(4 + first_ol_mb_col * 4 - first_ol_blk_col, 2);

    /* Check if new_row and new_col are evenly divisible by 4 (Q3),
     * and if so we shouldn't check neighboring blocks
     */
    if (new_row >= 0 && (new_row & 0x1F) == 0)
        end_row = 1;
    if (new_col >= 0 && (new_col & 0x1F) == 0)
        end_col = 1;

    /* Avoid calculating overlap for blocks in the previous MB */
    if (new_row < (first_ol_mb_row*16)<<3)
        end_row = 1;
    if (new_col < (first_ol_mb_col*16)<<3)
        end_col = 1;

    for (row = 0; row < end_row; ++row)
    {
        for (col = 0; col < end_col; ++col)
        {
            /* input in Q3, result in Q6 */
            const int overlap = vp8_block_overlap(new_row, new_col,
                                                  (((first_ol_blk_row + row) *
                                                      4) << 3),
                                                  (((first_ol_blk_col + col) *
                                                      4) << 3));
            vp8_assign_overlap(b_ol_ul[row * 4 + col].overlaps,
                               bmi,
                               ref_frame,
                               overlap);
        }
    }
}

void vp8_calculate_overlaps_submb(MB_OVERLAP *overlap_ul,
                                  int mb_rows, int mb_cols,
                                  B_MODE_INFO *bmi,
                                  MV_REFERENCE_FRAME ref_frame,
                                  int b_row, int b_col)
{
    MB_OVERLAP *mb_overlap;
    int row, col, rel_row, rel_col;
    int new_row, new_col;
    int new_row_pos, new_col_pos;
    int end_row, end_col;
    int overlap_b_row, overlap_b_col;
    int overlap_mb_row, overlap_mb_col;
    int i;
    B_MODE_INFO *obmi;
    int overlap;

    if (ref_frame == INTRA_FRAME)
        return;

    /* mb subpixel position */
    row = (4 * b_row) << 3; /* Q3 */
    col = (4 * b_col) << 3; /* Q3 */

    /* reverse compensate for motion */
    new_row = row - bmi->mv.as_mv.row;
    new_col = col - bmi->mv.as_mv.col;

    if (new_row >= ((16*mb_rows) << 3) || new_col >= ((16*mb_cols) << 3))
    {
        /* the new block ended up outside the frame */
        return;
    }

    if (new_row <= (-4 << 3) || new_col <= (-4 << 3))
    {
        /* outside the frame */
        return;
    }
    /* overlapping block's position in blocks */
    overlap_b_row = FLOOR(new_row / 4, 3) >> 3;
    overlap_b_col = FLOOR(new_col / 4, 3) >> 3;

    /* overlapping block's MB position in MBs
     * operations are done in Q3
     */
    overlap_mb_row = FLOOR((overlap_b_row << 3) / 4, 3) >> 3;
    overlap_mb_col = FLOOR((overlap_b_col << 3) / 4, 3) >> 3;

    end_row = MIN(mb_rows - overlap_mb_row, 2);
    end_col = MIN(mb_cols - overlap_mb_col, 2);

    /* Don't calculate overlap for MBs we don't overlap */
    /* Check if the new block row starts at the last block row of the MB */
    if (abs(new_row - ((16*overlap_mb_row) << 3)) < ((3*4) << 3))
        end_row = 1;
    /* Check if the new block col starts at the last block col of the MB */
    if (abs(new_col - ((16*overlap_mb_col) << 3)) < ((3*4) << 3))
        end_col = 1;
    /* Check if our MV is even (only overlapping one block) */
//    if (abs(bmi->mv.as_mv.row) % (4 << 3) == 0)
//        end_row = 1;
//    if (abs(bmi->mv.as_mv.col) % (4 << 3) == 0)
//        end_col = 1;

    /* find the MB(s) this block is overlapping */
    for (rel_row = 0; rel_row < end_row; ++rel_row)
    {
        for (rel_col = 0; rel_col < end_col; ++rel_col)
        {
            if (overlap_mb_row + rel_row < 0 ||
                overlap_mb_col + rel_col < 0)
                continue;
            mb_overlap = overlap_ul + (overlap_mb_row + rel_row) * mb_cols +
                 overlap_mb_col + rel_col;

            vp8_calculate_overlaps_mb(mb_overlap->overlaps, bmi, ref_frame,
                                      new_row, new_col,
                                      overlap_mb_row + rel_row,
                                      overlap_mb_col + rel_col,
                                      overlap_b_row + rel_row,
                                      overlap_b_col + rel_col);
        }
    }
}

MV_REFERENCE_FRAME vp8_largest_overlap_type(const B_OVERLAP *block_overlaps)
{
    int i, j;
    int overlap_per_type[MAX_REF_FRAMES] = {0};
    int largest_overlap = 0;
    MV_REFERENCE_FRAME largest_overlap_type = LAST_FRAME;
    for (i=0; i < 16; ++i)
    {
        const OVERLAP_NODE *overlaps = block_overlaps->overlaps;
        for (j=0; j < MAX_OVERLAPS; ++j)
        {
            if (overlaps[j].bmi != NULL)
            {
                overlap_per_type[overlaps[j].ref_frame] += overlaps[j].overlap;
                if (overlap_per_type[overlaps[j].ref_frame] > largest_overlap)
                {
                    largest_overlap = overlap_per_type[overlaps[j].ref_frame];
                    largest_overlap_type = overlaps[j].ref_frame;
                }
                assert(overlaps[j].overlap < (16*16)<<6);
            }
        }
        ++block_overlaps;
    }
    return largest_overlap_type;
}

void vp8_estimate_mv(const OVERLAP_NODE *overlaps, B_MODE_INFO *bmi,
                                   MV_REFERENCE_FRAME type)
{
    int i;
    int overlap_sum = 0;
    int row_acc = 0;
    int col_acc = 0;

    bmi->mv.as_int = 0;
    for (i=0; i < MAX_OVERLAPS; ++i)
    {
        if (overlaps[i].bmi != NULL &&
            overlaps[i].ref_frame == type)
        {
            col_acc += overlaps[i].overlap * overlaps[i].bmi->mv.as_mv.col;
            row_acc += overlaps[i].overlap * overlaps[i].bmi->mv.as_mv.row;
            overlap_sum += overlaps[i].overlap;
        }
    }
    if (overlap_sum > 0)
    {
        /* Q9 / Q6 = Q3 */
        bmi->mv.as_mv.col = col_acc / overlap_sum;
        bmi->mv.as_mv.row = row_acc / overlap_sum;
        /* TODO(holmer): Get the mode from the most overlapping block?
         * Or is the mode only used when decoding the MVs?
         */
        bmi->mode = NEW4X4;
    }
    else
    {
        bmi->mv.as_mv.col = 0;
        bmi->mv.as_mv.row = 0;
        bmi->mode = NEW4X4;
    }
}

void vp8_estimate_mb_mvs(const B_OVERLAP *block_overlaps,
                         B_MODE_INFO *bmi,
                         MV_REFERENCE_FRAME type,
                         MV* filtered_mv)
{
    int i;
    int non_zero_count = 0;
    filtered_mv->col = 0;
    filtered_mv->row = 0;
    for (i = 0; i < 16; ++i)
    {
        /* TODO(holmer): How can we be certain that all blocks refer
         * to the same frame buffer? We can't
         */
        /* Estimate vectors for all blocks which are overlapped by this
         * type
         */
        /* Interpolate/extrapolate the rest of the block's MVs */
        vp8_estimate_mv(block_overlaps[i].overlaps, bmi + i, type);
        if (bmi[i].mv.as_int != 0)
        {
            ++non_zero_count;
            filtered_mv->col += bmi[i].mv.as_mv.col;
            filtered_mv->row += bmi[i].mv.as_mv.row;
        }
    }
    if (non_zero_count > 0)
    {
        filtered_mv->col /= non_zero_count;
        filtered_mv->row /= non_zero_count;
    }
}

void vp8_estimate_missing_mvs(VP8D_COMP *pbi)
{
    VP8_COMMON * const pc = &pbi->common;
    vp8_estimate_missing_mvs_ex(pbi->overlaps,
                                pc->mi, pc->prev_mi,
                                pc->mb_rows, pc->mb_cols,
                                pbi->mvs_corrupt_from_mb);
}

void vp8_estimate_missing_mvs_ex(MB_OVERLAP *overlaps,
                                 MODE_INFO *mi, MODE_INFO *prev_mi,
                                 int mb_rows, int mb_cols,
                                 unsigned int first_corrupt)
{
    const unsigned int num_mbs = mb_rows * mb_cols;
    int mb_row, mb_col;
    vpx_memset(overlaps, 0, sizeof(MB_OVERLAP) * mb_rows * mb_cols);
    for (mb_row = 0; mb_row < mb_rows; ++mb_row)
    {
        for (mb_col = 0; mb_col < mb_cols; ++mb_col)
        {
            int sub_row;
            int sub_col;
            for (sub_row = 0; sub_row < 4; ++sub_row)
            {
                for (sub_col = 0; sub_col < 4; ++sub_col)
                {
                    vp8_calculate_overlaps_submb(
                                        overlaps, mb_rows, mb_cols,
                                        &(prev_mi->bmi[sub_row * 4 + sub_col]),
                                        prev_mi->mbmi.ref_frame,
                                        4 * mb_row + sub_row,
                                        4 * mb_col + sub_col);
                }
            }
            ++prev_mi;
        }
        ++prev_mi;
    }

    mb_row = first_corrupt / mb_cols;
    mb_col = first_corrupt - mb_row * mb_cols;
    mi += mb_row*(mb_cols + 1) + mb_col;
    for (; mb_row < mb_rows; ++mb_row)
    {
        for (; mb_col < mb_cols; ++mb_col)
        {
            int i;
            MV_REFERENCE_FRAME type = LAST_FRAME;
            int largest_overlap = 0;
            const B_OVERLAP *block_overlaps =
                    overlaps[mb_row*mb_cols + mb_col].overlaps;
            /* Find largest overlap and its type */
            mi->mbmi.ref_frame = vp8_largest_overlap_type(block_overlaps);
            vp8_estimate_mb_mvs(block_overlaps,
                                mi->bmi,
                                mi->mbmi.ref_frame,
                                &mi->mbmi.mv.as_mv);
            mi->mbmi.uv_mode = SPLITMV;

            mi->mbmi.mb_skip_coeff = 1;
            /* TODO(holmer): should this be enabled, when? */
            mi->mbmi.need_to_clamp_mvs = 1;
            ++mi;
        }
        mb_col = 0;
        ++mi;
    }
}

void vp8_conceal_corrupt_block(MACROBLOCKD *xd)
{
    /* this macroblock has corrupt residual, use the motion compensated
     image for concealment */
    int i;
    for (i = 0; i < 16; i++)
        vpx_memcpy(xd->dst.y_buffer + i * xd->dst.y_stride,
                   xd->predictor + i * 16, 16);
    for (i = 0; i < 8; i++)
        vpx_memcpy(xd->dst.u_buffer + i * xd->dst.uv_stride,
                   xd->predictor + 256 + i * 8, 8);
    for (i = 0; i < 8; i++)
        vpx_memcpy(xd->dst.v_buffer + i * xd->dst.uv_stride,
                   xd->predictor + 320 + i * 8, 8);
}
