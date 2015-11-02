/*
 *  Copyright (c) 2012 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */
#ifndef VP9_COMMON_VP9_MVREF_COMMON_H_
#define VP9_COMMON_VP9_MVREF_COMMON_H_

#include "vp9/common/vp9_onyxc_int.h"
#include "vp9/common/vp9_blockd.h"

#ifdef __cplusplus
extern "C" {
#endif

#define LEFT_TOP_MARGIN ((VP9_ENC_BORDER_IN_PIXELS - VP9_INTERP_EXTEND) << 3)
#define RIGHT_BOTTOM_MARGIN ((VP9_ENC_BORDER_IN_PIXELS -\
                                VP9_INTERP_EXTEND) << 3)

#define MVREF_NEIGHBOURS 8
#if CONFIG_NEWMVREF
#define MAX_ZONES 2
#endif  // CONFIG_NEWMVREF

typedef struct position {
  int row;
  int col;
} POSITION;

typedef enum {
  BOTH_ZERO = 0,
  ZERO_PLUS_PREDICTED = 1,
  BOTH_PREDICTED = 2,
  NEW_PLUS_NON_INTRA = 3,
  BOTH_NEW = 4,
  INTRA_PLUS_NON_INTRA = 5,
  BOTH_INTRA = 6,
  INVALID_CASE = 9
} motion_vector_context;

// This is used to figure out a context for the ref blocks. The code flattens
// an array that would have 3 possible counts (0, 1 & 2) for 3 choices by
// adding 9 for each intra block, 3 for each zero mv and 1 for each new
// motion vector. This single number is then converted into a context
// with a single lookup ( counter_to_context ).
static const int mode_2_counter[MB_MODE_COUNT] = {
  9,  // DC_PRED
  9,  // V_PRED
  9,  // H_PRED
  9,  // D45_PRED
  9,  // D135_PRED
  9,  // D117_PRED
  9,  // D153_PRED
  9,  // D207_PRED
  9,  // D63_PRED
  9,  // TM_PRED
#if CONFIG_INTRABC
  9,  // NEWDV
#endif  // CONFIG_INTRABC
  0,  // NEARESTMV
  0,  // NEARMV
  3,  // ZEROMV
  1,  // NEWMV
#if CONFIG_NEW_INTER
  1,  // NEW2MV
  0,  // NEAREST_NEARESTMV
  0,  // NEAREST_NEARMV
  0,  // NEAR_NEARESTMV
  1,  // NEAREST_NEWMV
  1,  // NEW_NEARESTMV
  1,  // NEAR_NEWMV
  1,  // NEW_NEARMV
  3,  // ZERO_ZEROMV
  1,  // NEW_NEWMV
#endif  // CONFIG_NEW_INTER
};

// There are 3^3 different combinations of 3 counts that can be either 0,1 or
// 2. However the actual count can never be greater than 2 so the highest
// counter we need is 18. 9 is an invalid counter that's never used.
static const int counter_to_context[19] = {
  BOTH_PREDICTED,  // 0
  NEW_PLUS_NON_INTRA,  // 1
  BOTH_NEW,  // 2
  ZERO_PLUS_PREDICTED,  // 3
  NEW_PLUS_NON_INTRA,  // 4
  INVALID_CASE,  // 5
  BOTH_ZERO,  // 6
  INVALID_CASE,  // 7
  INVALID_CASE,  // 8
  INTRA_PLUS_NON_INTRA,  // 9
  INTRA_PLUS_NON_INTRA,  // 10
  INVALID_CASE,  // 11
  INTRA_PLUS_NON_INTRA,  // 12
  INVALID_CASE,  // 13
  INVALID_CASE,  // 14
  INVALID_CASE,  // 15
  INVALID_CASE,  // 16
  INVALID_CASE,  // 17
  BOTH_INTRA  // 18
};

static const POSITION mv_ref_blocks[BLOCK_SIZES][MVREF_NEIGHBOURS] = {
  // 4X4
  {{-1, 0}, {0, -1}, {-1, -1}, {-2, 0}, {0, -2}, {-2, -1}, {-1, -2}, {-2, -2}},
  // 4X8
  {{-1, 0}, {0, -1}, {-1, -1}, {-2, 0}, {0, -2}, {-2, -1}, {-1, -2}, {-2, -2}},
  // 8X4
  {{-1, 0}, {0, -1}, {-1, -1}, {-2, 0}, {0, -2}, {-2, -1}, {-1, -2}, {-2, -2}},
  // 8X8
  {{-1, 0}, {0, -1}, {-1, -1}, {-2, 0}, {0, -2}, {-2, -1}, {-1, -2}, {-2, -2}},
  // 8X16
  {{0, -1}, {-1, 0}, {1, -1}, {-1, -1}, {0, -2}, {-2, 0}, {-2, -1}, {-1, -2}},
  // 16X8
  {{-1, 0}, {0, -1}, {-1, 1}, {-1, -1}, {-2, 0}, {0, -2}, {-1, -2}, {-2, -1}},
  // 16X16
  {{-1, 0}, {0, -1}, {-1, 1}, {1, -1}, {-1, -1}, {-3, 0}, {0, -3}, {-3, -3}},
  // 16X32
  {{0, -1}, {-1, 0}, {2, -1}, {-1, -1}, {-1, 1}, {0, -3}, {-3, 0}, {-3, -3}},
  // 32X16
  {{-1, 0}, {0, -1}, {-1, 2}, {-1, -1}, {1, -1}, {-3, 0}, {0, -3}, {-3, -3}},
  // 32X32
  {{-1, 1}, {1, -1}, {-1, 2}, {2, -1}, {-1, -1}, {-3, 0}, {0, -3}, {-3, -3}},
  // 32X64
  {{0, -1}, {-1, 0}, {4, -1}, {-1, 2}, {-1, -1}, {0, -3}, {-3, 0}, {2, -1}},
  // 64X32
  {{-1, 0}, {0, -1}, {-1, 4}, {2, -1}, {-1, -1}, {-3, 0}, {0, -3}, {-1, 2}},
  // 64X64
  {{-1, 3}, {3, -1}, {-1, 4}, {4, -1}, {-1, -1}, {-1, 0}, {0, -1}, {-1, 6}},
#if CONFIG_EXT_CODING_UNIT_SIZE
  // 64x128
  {{0, -1}, {-1, 0}, {4, -1}, {-1, 2}, {-1, -1}, {0, -3}, {-3, 0}, {2, -1}},
  // 128x64
  {{-1, 0}, {0, -1}, {-1, 4}, {2, -1}, {-1, -1}, {-3, 0}, {0, -3}, {-1, 2}},
  // 128x128
  {{-1, 3}, {3, -1}, {-1, 4}, {4, -1}, {-1, -1}, {-1, 0}, {0, -1}, {-1, 6}},
#endif
};

static const int idx_n_column_to_subblock[4][2] = {
  {1, 2},
  {1, 3},
  {3, 2},
  {3, 3}
};

#if CONFIG_NEWMVREF
static const POSITION mv_ref_blocks_8x8[MAX_ZONES][MVREF_NEIGHBOURS] = {
  {  // 8X8, Zone I,  where top right neighbors are available
    {-1,  0}, { 0, -1}, {-1,  1}, {-1, -1},  // nearest neighboring blocks
    {-2,  0}, { 0, -2}, {-2, -1}, {-1, -2}
  },
  {  // 8X8, Zone II, where no top right neighbor is available
    {-1,  0}, { 0, -1}, {-1, -1},            // nearest neighboring blocks
    {-2,  0}, { 0, -2}, {-2, -1}, {-1, -2}, {-2, -2}
  }
};

static const int mv_ref_topright_avail_8x8[8][8] = {
  {1, 1, 1, 1, 1, 1, 1, 1},
  {1, 0, 1, 0, 1, 0, 1, 0},
  {1, 1, 1, 0, 1, 1, 1, 0},
  {1, 0, 1, 0, 1, 0, 1, 0},
  {1, 1, 1, 1, 1, 1, 1, 0},
  {1, 0, 1, 0, 1, 0, 1, 0},
  {1, 1, 1, 0, 1, 1, 1, 0},
  {1, 0, 1, 0, 1, 0, 1, 0}
};

static const int idx_to_subblock_top_left[12][2][3] = {
  {  // 4x4 subblock 0 (current)
    {2, 0, 2},  // top:  4x4, 4x8, 8x4
    {1, 1, 0}   // left: 4x4, 4x8, 8x4
  },
  {  // 4x4 subblock 1 (current)
    {3, 1, 2},  // top:  4x4, 4x8, 8x4
    {1, 1, 0}   // left: 4x4, 4x8, 8x4
  },
  {  // 4x4 subblock 2 (current)
    {2, 0, 2},  // top:  4x4, 4x8, 8x4
    {3, 1, 2}   // left: 4x4, 4x8, 8x4
  },
  {  // 4x4 subblock 3 (current)
    {3, 1, 2},  // top:  4x4, 4x8, 8x4
    {3, 1, 2}   // left: 4x4, 4x8, 8x4
  },
  {  // 4x8 subblock 0 (current)
    { 2, 0,  2},  // top:  4x4, 4x8, 8x4
    {-1, 1, -1}   // left: 4x4, 4x8, 8x4
  },
  {  // 4x8 subblock 1 (current)
    { 3, 1,  2},  // top:  4x4, 4x8, 8x4
    {-1, 1, -1}   // left: 4x4, 4x8, 8x4
  },
  {  // 4x8 subblock 2 (current)
    { 2, 0,  2},  // top:  4x4, 4x8, 8x4
    {-1, 1, -1}   // left: 4x4, 4x8, 8x4
  },
  {  // 4x8 subblock 3 (current)
    { 3, 1,  2},  // top:  4x4, 4x8, 8x4
    {-1, 1, -1}   // left: 4x4, 4x8, 8x4
  },
  {  // 8x4 subblock 0 (current)
    {-1, -1, 2},  // top:  4x4, 4x8, 8x4
    { 1,  1, 0}   // left: 4x4, 4x8, 8x4
  },
  {  // 8x4 subblock 1 (current)
    {-1, -1, 2},  // top:  4x4, 4x8, 8x4
    { 1,  1, 0}   // left: 4x4, 4x8, 8x4
  },
  {  // 8x4 subblock 2 (current)
    {-1, -1, 2},  // top:  4x4, 4x8, 8x4
    { 3,  1, 2}   // left: 4x4, 4x8, 8x4
  },
  {  // 8x4 subblock 3 (current)
    {-1, -1, 2},  // top:  4x4, 4x8, 8x4
    { 3,  1, 2}   // left: 4x4, 4x8, 8x4
  }
};

static const int idx_to_subblock_topright_topleft[2][3] = {
  {2, 0, 2},  // top-right: 4x4, 4x8, 8x4
  {3, 1, 2}   // top-left:  4x4, 4x8, 8x4
};
#endif  // CONFIG_NEWMVREF

// clamp_mv_ref
#if CONFIG_EXT_CODING_UNIT_SIZE
#define MV_BORDER (32 << 3)  // Allow 32 pels in 1/8th pel units
#else
#define MV_BORDER (16 << 3)  // Allow 16 pels in 1/8th pel units
#endif

static INLINE void clamp_mv_ref(MV *mv, const MACROBLOCKD *xd) {
  clamp_mv(mv, xd->mb_to_left_edge - MV_BORDER,
               xd->mb_to_right_edge + MV_BORDER,
               xd->mb_to_top_edge - MV_BORDER,
               xd->mb_to_bottom_edge + MV_BORDER);
}

// This function returns either the appropriate sub block or block's mv
// on whether the block_size < 8x8 and we have check_sub_blocks set.
static INLINE int_mv get_sub_block_mv(const MODE_INFO *candidate, int which_mv,
                                      int search_col, int block_idx) {
  return block_idx >= 0 && candidate->mbmi.sb_type < BLOCK_8X8
          ? candidate->bmi[idx_n_column_to_subblock[block_idx][search_col == 0]]
              .as_mv[which_mv]
          : candidate->mbmi.mv[which_mv];
}

// Performs mv sign inversion if indicated by the reference frame combination.
static INLINE int_mv scale_mv(const MB_MODE_INFO *mbmi, int ref,
                              const MV_REFERENCE_FRAME this_ref_frame,
                              const int *ref_sign_bias) {
  int_mv mv = mbmi->mv[ref];
  if (ref_sign_bias[mbmi->ref_frame[ref]] != ref_sign_bias[this_ref_frame]) {
    mv.as_mv.row *= -1;
    mv.as_mv.col *= -1;
  }
  return mv;
}

// This macro is used to add a motion vector mv_ref list if it isn't
// already in the list.  If it's the second motion vector it will also
// skip all additional processing and jump to done!
#define ADD_MV_REF_LIST(mv) \
  do { \
    if (refmv_count) { \
      if ((mv).as_int != mv_ref_list[0].as_int) { \
        mv_ref_list[refmv_count] = (mv); \
        goto Done; \
      } \
    } else { \
      mv_ref_list[refmv_count++] = (mv); \
    } \
  } while (0)

// If either reference frame is different, not INTRA, and they
// are different from each other scale and add the mv to our list.
#define IF_DIFF_REF_FRAME_ADD_MV(mbmi) \
  do { \
    if (is_inter_block(mbmi)) { \
      if ((mbmi)->ref_frame[0] != ref_frame) \
        ADD_MV_REF_LIST(scale_mv((mbmi), 0, ref_frame, ref_sign_bias)); \
      if (has_second_ref(mbmi) && \
          (mbmi)->ref_frame[1] != ref_frame && \
          (mbmi)->mv[1].as_int != (mbmi)->mv[0].as_int) \
        ADD_MV_REF_LIST(scale_mv((mbmi), 1, ref_frame, ref_sign_bias)); \
    } \
  } while (0)


// Checks that the given mi_row, mi_col and search point
// are inside the borders of the tile.
static INLINE int is_inside(const TileInfo *const tile,
                            int mi_col, int mi_row, int mi_rows,
                            const POSITION *mi_pos) {
#if CONFIG_ROW_TILE
  (void) mi_rows;
  return !(mi_row + mi_pos->row < tile->mi_row_start ||
           mi_col + mi_pos->col < tile->mi_col_start ||
           mi_row + mi_pos->row >= tile->mi_row_end ||
           mi_col + mi_pos->col >= tile->mi_col_end);
#else
  return !(mi_row + mi_pos->row < 0 ||
           mi_col + mi_pos->col < tile->mi_col_start ||
           mi_row + mi_pos->row >= mi_rows ||
           mi_col + mi_pos->col >= tile->mi_col_end);
#endif
}

#if CONFIG_NEWMVREF
// This macro is used to add a motion vector as a candidate for the mv ref if
// it isn't already taken. If it's the third motion vector, it will also skip
// all additional processing and jump to done!
#define ADD_MV_REF_CANDIDATE(mv) \
  do { \
    if (refmv_count) { \
      if (refmv_count == 1 && \
          (mv).as_int != mv_ref_candidates[0].as_int) { \
        mv_ref_candidates[refmv_count++] = (mv); \
      } else if (refmv_count == 2 && \
                 (mv).as_int != mv_ref_candidates[0].as_int && \
                 (mv).as_int != mv_ref_candidates[1].as_int) { \
        mv_ref_candidates[refmv_count] = (mv); \
        goto Done; \
      } \
    } else { \
      mv_ref_candidates[refmv_count++] = (mv); \
    } \
  } while (0)

// If either reference frame is different, not INTRA, and they
// are different from each other scale and add the mv as candidate.
#define IF_DIFF_REF_FRAME_ADD_MV_CANDIDATE(mbmi) \
  do { \
    if (is_inter_block(mbmi)) { \
      if ((mbmi)->ref_frame[0] != ref_frame) \
        ADD_MV_REF_CANDIDATE(scale_mv((mbmi), 0, ref_frame, ref_sign_bias)); \
      if (has_second_ref(mbmi) && \
          (mbmi)->ref_frame[1] != ref_frame && \
          (mbmi)->mv[1].as_int != (mbmi)->mv[0].as_int) \
        ADD_MV_REF_CANDIDATE(scale_mv((mbmi), 1, ref_frame, ref_sign_bias)); \
    } \
  } while (0)
#endif  // CONFIG_NEWMVREF

// TODO(jingning): this mv clamping function should be block size dependent.
static INLINE void clamp_mv2(MV *mv, const MACROBLOCKD *xd) {
  clamp_mv(mv, xd->mb_to_left_edge - LEFT_TOP_MARGIN,
               xd->mb_to_right_edge + RIGHT_BOTTOM_MARGIN,
               xd->mb_to_top_edge - LEFT_TOP_MARGIN,
               xd->mb_to_bottom_edge + RIGHT_BOTTOM_MARGIN);
}

#if CONFIG_NEW_INTER
// This function keeps a mode count for a given MB/SB
void vp9_update_mv_context(const VP9_COMMON *cm, const MACROBLOCKD *xd,
                           const TileInfo *const tile,
                           MODE_INFO *mi, MV_REFERENCE_FRAME ref_frame,
                           int_mv *mv_ref_list,
                           int block, int mi_row, int mi_col);
#endif  // CONFIG_NEW_INTER

void vp9_find_mv_refs(const VP9_COMMON *cm, const MACROBLOCKD *xd,
                      const TileInfo *const tile,
                      MODE_INFO *mi, MV_REFERENCE_FRAME ref_frame,
                      int_mv *mv_ref_list, int mi_row, int mi_col);

static INLINE void vp9_lower_mv_precision(MV *mv, const int usehp) {
  if (!usehp) {
    if (mv->row & 1)
      mv->row += (mv->row > 0 ? -1 : 1);
    if (mv->col & 1)
      mv->col += (mv->col > 0 ? -1 : 1);
  }
}

// check a list of motion vectors by sad score using a number rows of pixels
// above and a number cols of pixels in the left to select the one with best
// score to use as ref motion vector
void vp9_find_best_ref_mvs(MACROBLOCKD *xd, int allow_hp,
                           int_mv *mvlist, int_mv *nearest, int_mv *near);

void vp9_append_sub8x8_mvs_for_idx(VP9_COMMON *cm, MACROBLOCKD *xd,
                                   const TileInfo *const tile,
                                   int block, int ref, int mi_row, int mi_col,
#if CONFIG_NEW_INTER
                                   int_mv *mv_list,
#endif  // CONFIG_NEW_INTER
                                   int_mv *nearest, int_mv *near);

#if CONFIG_COPY_MODE
int vp9_construct_ref_inter_list(VP9_COMMON *cm,  MACROBLOCKD *xd,
                                 const TileInfo *const tile,
                                 BLOCK_SIZE bsize,
#if CONFIG_EXT_PARTITION
                                 PARTITION_TYPE partition,
#endif
                                 int mi_row, int mi_col,
                                 MB_MODE_INFO *ref_list[18]);
#endif  // CONFIG_COPY_MODE

#if CONFIG_INTRABC
void vp9_find_ref_dv(int_mv *ref_dv, int mi_row, int mi_col);
#endif  // CONFIOG_INTRABC
#ifdef __cplusplus
}  // extern "C"
#endif

#endif  // VP9_COMMON_VP9_MVREF_COMMON_H_
