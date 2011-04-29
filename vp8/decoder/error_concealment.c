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

#define NUM_NEIGHBORS 20

typedef struct ec_position
{
    int row;
    int col;
} EC_POS;

/*
 * Regenerate the table in Matlab with:
 * x = meshgrid((1:4), (1:4));
 * y = meshgrid((1:4), (1:4))';
 * W = round((1./(sqrt(x.^2 + y.^2))*2^7));
 * W(1,1) = 0;
 */
static const int weights_q7[5][5] = {
       {  0,   128,    64,    43,    32 },
       {128,    91,    57,    40,    31 },
       { 64,    57,    45,    36,    29 },
       { 43,    40,    36,    30,    26 },
       { 32,    31,    29,    26,    23 }
};

static int vp8_need_to_clamp_mv(MV *mv,
                         int mb_to_left_edge,
                         int mb_to_right_edge,
                         int mb_to_top_edge,
                         int mb_to_bottom_edge)
{
    return (mv->col < mb_to_left_edge) ||
           (mv->col > mb_to_right_edge) ||
           (mv->row < mb_to_top_edge) ||
           (mv->row > mb_to_bottom_edge);
}

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
                               int mb_row, int mb_col,
                               int first_blk_row, int first_blk_col)
{
    /* Find the blocks within this MB which are overlapped by bmi and
     * calculate and assign overlap for each of those blocks. */

    /* Block coordinates relative the upper-left block */
    const int rel_ol_blk_row = first_blk_row - mb_row * 4;
    const int rel_ol_blk_col = first_blk_col - mb_col * 4;
    /* If the block partly overlaps any previous MB, these coordinates
     * can be < 0. We don't want to access blocks in previous MBs.
     */
    const int blk_idx = MAX(rel_ol_blk_row,0) * 4 + MAX(rel_ol_blk_col,0);
    /* Upper left overlapping block */
    B_OVERLAP *b_ol_ul = &(b_overlaps[blk_idx]);

    /* Calculate and assign overlaps for all blocks in this MB
     * which the motion compensated block overlaps
     */
    /* Avoid calculating overlaps for blocks in later MBs */
    int end_row = MIN(4 + mb_row * 4 - first_blk_row, 2);
    int end_col = MIN(4 + mb_col * 4 - first_blk_col, 2);
    int row, col;

    /* Check if new_row and new_col are evenly divisible by 4 (Q3),
     * and if so we shouldn't check neighboring blocks
     */
    if (new_row >= 0 && (new_row & 0x1F) == 0)
        end_row = 1;
    if (new_col >= 0 && (new_col & 0x1F) == 0)
        end_col = 1;

    /* Check if the overlapping block partly overlaps a previous MB
     * and if so, we're overlapping fewer blocks in this MB.
     */
    if (new_row < (mb_row*16)<<3)
        end_row = 1;
    if (new_col < (mb_col*16)<<3)
        end_col = 1;

    for (row = 0; row < end_row; ++row)
    {
        for (col = 0; col < end_col; ++col)
        {
            /* input in Q3, result in Q6 */
            const int overlap = vp8_block_overlap(new_row, new_col,
                                                  (((first_blk_row + row) *
                                                      4) << 3),
                                                  (((first_blk_col + col) *
                                                      4) << 3));
            vp8_assign_overlap(b_ol_ul[row * 4 + col].overlaps,
                               bmi,
                               ref_frame,
                               overlap);
        }
    }
}

void vp8_calculate_overlaps(MB_OVERLAP *overlap_ul,
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
                         MODE_INFO *mi,
                         int mb_to_left_edge,
                         int mb_to_right_edge,
                         int mb_to_top_edge,
                         int mb_to_bottom_edge)
{
    int i;
    int non_zero_count = 0;
    MV * const filtered_mv = &(mi->mbmi.mv.as_mv);
    B_MODE_INFO * const bmi = mi->bmi;
    filtered_mv->col = 0;
    filtered_mv->row = 0;
    for (i = 0; i < 16; ++i)
    {
        /* Estimate vectors for all blocks which are overlapped by this
         * type
         */
        /* Interpolate/extrapolate the rest of the block's MVs */
        vp8_estimate_mv(block_overlaps[i].overlaps, bmi + i,
                mi->mbmi.ref_frame);
        mi->mbmi.need_to_clamp_mvs = vp8_need_to_clamp_mv(&(bmi[i].mv.as_mv),
                                                          mb_to_left_edge,
                                                          mb_to_right_edge,
                                                          mb_to_top_edge,
                                                          mb_to_bottom_edge);
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
    /* First calculate the overlaps for all blocks */
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
                    vp8_calculate_overlaps(
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
    /* Go through all macroblocks in the current image with missing MVs
     * and calculate new MVs using the overlaps.
     */
    for (; mb_row < mb_rows; ++mb_row)
    {
        int mb_to_top_edge = -((mb_row * 16)) << 3;
        int mb_to_bottom_edge = ((mb_rows - 1 - mb_row) * 16) << 3;
        for (; mb_col < mb_cols; ++mb_col)
        {
            int mb_to_left_edge = -((mb_col * 16) << 3);
            int mb_to_right_edge = ((mb_cols - 1 - mb_col) * 16) << 3;
            int i;
            MV_REFERENCE_FRAME type = LAST_FRAME;
            int largest_overlap = 0;
            const B_OVERLAP *block_overlaps =
                    overlaps[mb_row*mb_cols + mb_col].overlaps;
            /* Find largest overlap and its type */
            mi->mbmi.ref_frame = vp8_largest_overlap_type(block_overlaps);
            vp8_estimate_mb_mvs(block_overlaps,
                                mi,
                                mb_to_left_edge,
                                mb_to_right_edge,
                                mb_to_top_edge,
                                mb_to_bottom_edge);
            mi->mbmi.mode = SPLITMV;
            mi->mbmi.uv_mode = DC_PRED;
            mi->mbmi.partitioning = 3;
            ++mi;
        }
        mb_col = 0;
        ++mi;
    }
}

static void assign_neighbor(EC_BLOCK *neighbor, MODE_INFO *mi, int block_idx)
{
    assert(mi->mbmi.ref_frame < MAX_REF_FRAMES);
    neighbor->ref_frame = mi->mbmi.ref_frame;
    neighbor->mv = mi->bmi[block_idx].mv.as_mv;
}

void vp8_find_neighboring_blocks(MODE_INFO *mi,
                                 EC_BLOCK *neighbors,
                                 int mb_row, int mb_col,
                                 int mb_rows, int mb_cols,
                                 int mi_stride)
{
    int i = 0;
    int j;
    if (mb_row > 0)
    {
        /* upper left */
        if (mb_col > 0)
            assign_neighbor(&neighbors[i], mi - mi_stride - 1, 15);
        ++i;
        /* above */
        for (j = 12; j < 16; ++j, ++i)
            assign_neighbor(&neighbors[i], mi - mi_stride, j);
    }
    else
        i += 5;
    if (mb_col < mb_cols - 1)
    {
        /* upper right */
        if (mb_row > 0)
            assign_neighbor(&neighbors[i], mi - mi_stride + 1, 12);
        ++i;
        /* right */
        for (j = 0; j <= 12; j += 4, ++i)
            assign_neighbor(&neighbors[i], mi + 1, j);
    }
    else
        i += 5;
    if (mb_row < mb_rows - 1)
    {
        /* lower right */
        if (mb_col < mb_cols - 1)
            assign_neighbor(&neighbors[i], mi + mi_stride + 1, 0);
        ++i;
        /* below */
        for (j = 0; j < 4; ++j, ++i)
            assign_neighbor(&neighbors[i], mi + mi_stride, j);
    }
    else
        i += 5;
    if (mb_col > 0)
    {
        /* lower left */
        if (mb_row < mb_rows - 1)
            assign_neighbor(&neighbors[i], mi + mi_stride - 1, 4);
        ++i;
        /* left */
        for (j = 3; j < 16; j += 4, ++i)
        {
            assign_neighbor(&neighbors[i], mi - 1, j);
        }
    }
    else
        i += 5;
    assert(i == 20);
}

MV_REFERENCE_FRAME vp8_dominant_ref_frame(EC_BLOCK *neighbors)
{
    /* Default to referring to "skip" */
    MV_REFERENCE_FRAME dom_ref_frame = LAST_FRAME;
    int max_ref_frame_cnt = 0;
    int ref_frame_cnt[MAX_REF_FRAMES] = {0};
    int i;
    /* Count neighboring reference frames */
    for (i = 0; i < NUM_NEIGHBORS; ++i)
    {
        if (neighbors[i].ref_frame < MAX_REF_FRAMES)
            ++ref_frame_cnt[neighbors[i].ref_frame];
    }
    /* Find maximum */
    for (i = 0; i < MAX_REF_FRAMES; ++i)
    {
        if (ref_frame_cnt[i] > max_ref_frame_cnt)
        {
            dom_ref_frame = i;
            max_ref_frame_cnt = ref_frame_cnt[i];
        }
    }
    return dom_ref_frame;
}

void vp8_interpolate_mvs(MACROBLOCKD *mb,
                         EC_BLOCK *neighbors,
                         MV_REFERENCE_FRAME dom_ref_frame)
{
    int row, col, i;
    MODE_INFO * const mi = mb->mode_info_context;
    /* Table with the position of the neighboring blocks relative the position
     * of the upper left block of the current MB. Starting with the upper left
     * neighbor and going to the right.
     */
    const EC_POS neigh_pos[NUM_NEIGHBORS] = {
                                        {-1,-1}, {-1,0}, {-1,1}, {-1,2}, {-1,3},
                                        {-1,4}, {0,4}, {1,4}, {2,4}, {3,4},
                                        {4,4}, {4,3}, {4,2}, {4,1}, {4,0},
                                        {4,-1}, {3,-1}, {2,-1}, {1,-1}, {0,-1}
                                      };
    for (row = 0; row < 4; ++row)
    {
        for (col = 0; col < 4; ++col)
        {
            int w_sum = 0;
            int mv_row_sum = 0;
            int mv_col_sum = 0;
            MV * const mv = &(mi->bmi[row*4 + col].mv.as_mv);
            for (i = 0; i < NUM_NEIGHBORS; ++i)
            {
                /* Calculate the weighted sum of neighboring MVs referring
                 * to the dominant frame type.
                 */
                const int w = weights_q7[abs(row - neigh_pos[i].row)]
                                        [abs(col - neigh_pos[i].col)];
                if (neighbors[i].ref_frame != dom_ref_frame)
                    continue;
                w_sum += w;
                /* Q7 * Q3 = Q10 */
                mv_row_sum += w*neighbors[i].mv.row;
                mv_col_sum += w*neighbors[i].mv.col;
            }
            if (w_sum > 0)
            {
                /* Avoid division by zero.
                 * Normalize with the sum of the coefficients
                 * Q3 = Q10 / Q7
                 */
                mv->row = mv_row_sum / w_sum;
                mv->col = mv_col_sum / w_sum;
                mi->bmi[row*4 + col].mode = NEW4X4;
                mi->mbmi.need_to_clamp_mvs = vp8_need_to_clamp_mv(mv,
                                                       mb->mb_to_left_edge,
                                                       mb->mb_to_right_edge,
                                                       mb->mb_to_top_edge,
                                                       mb->mb_to_bottom_edge);
            }
        }
    }
}

void vp8_interpolate_motion(MACROBLOCKD *mb,
                        int mb_row, int mb_col,
                        int mb_rows, int mb_cols,
                        int mi_stride)
{
    /* Find relevant neighboring blocks */
    EC_BLOCK neighbors[NUM_NEIGHBORS];
    MV_REFERENCE_FRAME dom_ref_frame;
    int i;
    /* Initialize the array. MAX_REF_FRAMES is interpreted as "doesn't exist" */
    for (i = 0; i < NUM_NEIGHBORS; ++i)
    {
        neighbors[i].ref_frame = MAX_REF_FRAMES;
        neighbors[i].mv.row = neighbors[i].mv.col = 0;
    }
    vp8_find_neighboring_blocks(mb->mode_info_context,
                                neighbors,
                                mb_row, mb_col,
                                mb_rows, mb_cols,
                                mb->mode_info_stride);
    /* Determine the dominant block type */
    dom_ref_frame = vp8_dominant_ref_frame(neighbors);
    /* Interpolate MVs for the missing blocks
     * from the dominating MVs */
    vp8_interpolate_mvs(mb, neighbors, dom_ref_frame);

    mb->mode_info_context->mbmi.ref_frame = dom_ref_frame;
    mb->mode_info_context->mbmi.mode = SPLITMV;
    mb->mode_info_context->mbmi.uv_mode = DC_PRED;
    mb->mode_info_context->mbmi.partitioning = 3;
}
