/*
 * Copyright (c) 2016, Alliance for Open Media. All rights reserved
 *
 * This source code is subject to the terms of the BSD 2 Clause License and
 * the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
 * was not distributed with this source code in the LICENSE file, you can
 * obtain it at www.aomedia.org/license/software. If the Alliance for Open
 * Media Patent License 1.0 was not distributed with this source code in the
 * PATENTS file, you can obtain it at www.aomedia.org/license/patent.
 */
#ifndef AV1_COMMON_MVREF_COMMON_H_
#define AV1_COMMON_MVREF_COMMON_H_

#include "av1/common/onyxc_int.h"
#include "av1/common/blockd.h"

#ifdef __cplusplus
extern "C" {
#endif

#if CONFIG_SIMP_MV_PRED
#define MVREF_NEIGHBOURS 9
#else
#define MVREF_NEIGHBOURS 8
#endif

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
  0,  // NEARESTMV
  0,  // NEARMV
  3,  // ZEROMV
  1,  // NEWMV
#if CONFIG_EXT_INTER
  1,    // NEWFROMNEARMV
  0,    // NEAREST_NEARESTMV
  0,    // NEAREST_NEARMV
  0,    // NEAR_NEARESTMV
  0,    // NEAR_NEARMV
  1,    // NEAREST_NEWMV
  1,    // NEW_NEARESTMV
  1,    // NEAR_NEWMV
  1,    // NEW_NEARMV
  3,    // ZERO_ZEROMV
  1,    // NEW_NEWMV
#endif  // CONFIG_EXT_INTER
};

// There are 3^3 different combinations of 3 counts that can be either 0,1 or
// 2. However the actual count can never be greater than 2 so the highest
// counter we need is 18. 9 is an invalid counter that's never used.
static const int counter_to_context[19] = {
  BOTH_PREDICTED,        // 0
  NEW_PLUS_NON_INTRA,    // 1
  BOTH_NEW,              // 2
  ZERO_PLUS_PREDICTED,   // 3
  NEW_PLUS_NON_INTRA,    // 4
  INVALID_CASE,          // 5
  BOTH_ZERO,             // 6
  INVALID_CASE,          // 7
  INVALID_CASE,          // 8
  INTRA_PLUS_NON_INTRA,  // 9
  INTRA_PLUS_NON_INTRA,  // 10
  INVALID_CASE,          // 11
  INTRA_PLUS_NON_INTRA,  // 12
  INVALID_CASE,          // 13
  INVALID_CASE,          // 14
  INVALID_CASE,          // 15
  INVALID_CASE,          // 16
  INVALID_CASE,          // 17
  BOTH_INTRA             // 18
};

#if !CONFIG_SIMP_MV_PRED
static const POSITION mv_ref_blocks[BLOCK_SIZES][MVREF_NEIGHBOURS] = {
  // 4X4
  { { -1, 0 },
    { 0, -1 },
    { -1, -1 },
    { -2, 0 },
    { 0, -2 },
    { -2, -1 },
    { -1, -2 },
    { -2, -2 } },
  // 4X8
  { { -1, 0 },
    { 0, -1 },
    { -1, -1 },
    { -2, 0 },
    { 0, -2 },
    { -2, -1 },
    { -1, -2 },
    { -2, -2 } },
  // 8X4
  { { -1, 0 },
    { 0, -1 },
    { -1, -1 },
    { -2, 0 },
    { 0, -2 },
    { -2, -1 },
    { -1, -2 },
    { -2, -2 } },
  // 8X8
  { { -1, 0 },
    { 0, -1 },
    { -1, -1 },
    { -2, 0 },
    { 0, -2 },
    { -2, -1 },
    { -1, -2 },
    { -2, -2 } },
  // 8X16
  { { 0, -1 },
    { -1, 0 },
    { 1, -1 },
    { -1, -1 },
    { 0, -2 },
    { -2, 0 },
    { -2, -1 },
    { -1, -2 } },
  // 16X8
  { { -1, 0 },
    { 0, -1 },
    { -1, 1 },
    { -1, -1 },
    { -2, 0 },
    { 0, -2 },
    { -1, -2 },
    { -2, -1 } },
  // 16X16
  { { -1, 0 },
    { 0, -1 },
    { -1, 1 },
    { 1, -1 },
    { -1, -1 },
    { -3, 0 },
    { 0, -3 },
    { -3, -3 } },
  // 16X32
  { { 0, -1 },
    { -1, 0 },
    { 2, -1 },
    { -1, -1 },
    { -1, 1 },
    { 0, -3 },
    { -3, 0 },
    { -3, -3 } },
  // 32X16
  { { -1, 0 },
    { 0, -1 },
    { -1, 2 },
    { -1, -1 },
    { 1, -1 },
    { -3, 0 },
    { 0, -3 },
    { -3, -3 } },
  // 32X32
  { { -1, 1 },
    { 1, -1 },
    { -1, 2 },
    { 2, -1 },
    { -1, -1 },
    { -3, 0 },
    { 0, -3 },
    { -3, -3 } },
  // 32X64
  { { 0, -1 },
    { -1, 0 },
    { 4, -1 },
    { -1, 2 },
    { -1, -1 },
    { 0, -3 },
    { -3, 0 },
    { 2, -1 } },
  // 64X32
  { { -1, 0 },
    { 0, -1 },
    { -1, 4 },
    { 2, -1 },
    { -1, -1 },
    { -3, 0 },
    { 0, -3 },
    { -1, 2 } },
  // 64X64
  { { -1, 3 },
    { 3, -1 },
    { -1, 4 },
    { 4, -1 },
    { -1, -1 },
    { -1, 0 },
    { 0, -1 },
    { -1, 6 } },
#if CONFIG_EXT_PARTITION
  // TODO(debargha/jingning) Making them twice the 32x64, .. ones above
  // 64x128
  { { 0, -2 },
    { -2, 0 },
    { 8, -2 },
    { -2, 4 },
    { -2, -2 },
    { 0, -6 },
    { -6, 0 },
    { 4, -2 } },
  // 128x64
  { { -2, 0 },
    { 0, -2 },
    { -2, 8 },
    { 4, -2 },
    { -2, -2 },
    { -6, 0 },
    { 0, -6 },
    { -2, 4 } },
  // 128x128
  { { -2, 6 },
    { 6, -2 },
    { -2, 8 },
    { 8, -2 },
    { -2, -2 },
    { -2, 0 },
    { 0, -2 },
    { -2, 12 } },
#endif  // CONFIG_EXT_PARTITION
};
#endif

static const int idx_n_column_to_subblock[4][2] = {
  { 1, 2 }, { 1, 3 }, { 3, 2 }, { 3, 3 }
};

// clamp_mv_ref
#if CONFIG_EXT_PARTITION
#define MV_BORDER (16 << 3)  // Allow 16 pels in 1/8th pel units
#else
#define MV_BORDER (8 << 3)  // Allow 8 pels in 1/8th pel units
#endif                      // CONFIG_EXT_PARTITION

static INLINE void clamp_mv_ref(MV *mv, int bw, int bh, const MACROBLOCKD *xd) {
  clamp_mv(mv, xd->mb_to_left_edge - bw * 8 - MV_BORDER,
           xd->mb_to_right_edge + bw * 8 + MV_BORDER,
           xd->mb_to_top_edge - bh * 8 - MV_BORDER,
           xd->mb_to_bottom_edge + bh * 8 + MV_BORDER);
}

// This function returns either the appropriate sub block or block's mv
// on whether the block_size < 8x8 and we have check_sub_blocks set.
static INLINE int_mv get_sub_block_mv(const MODE_INFO *candidate, int which_mv,
                                      int search_col, int block_idx) {
#if CONFIG_SIMP_MV_PRED
  return candidate->mbmi.mv[which_mv];
#else
  return block_idx >= 0 && candidate->mbmi.sb_type < BLOCK_8X8
             ? candidate
                   ->bmi[idx_n_column_to_subblock[block_idx][search_col == 0]]
                   .as_mv[which_mv]
             : candidate->mbmi.mv[which_mv];
#endif
}

#if CONFIG_REF_MV
static INLINE int_mv get_sub_block_pred_mv(const MODE_INFO *candidate,
                                           int which_mv, int search_col,
                                           int block_idx) {
#if CONFIG_SIMP_MV_PRED
  return candidate->mbmi.mv[which_mv];
#else
  return block_idx >= 0 && candidate->mbmi.sb_type < BLOCK_8X8
             ? candidate
                   ->bmi[idx_n_column_to_subblock[block_idx][search_col == 0]]
                   .pred_mv[which_mv]
             : candidate->mbmi.pred_mv[which_mv];
#endif
}
#endif

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

#define CLIP_IN_ADD(mv, bw, bh, xd) clamp_mv_ref(mv, bw, bh, xd)

// This macro is used to add a motion vector mv_ref list if it isn't
// already in the list.  If it's the second motion vector it will also
// skip all additional processing and jump to done!
#define ADD_MV_REF_LIST(mv, refmv_count, mv_ref_list, bw, bh, xd, Done)      \
  do {                                                                       \
    (mv_ref_list)[(refmv_count)] = (mv);                                     \
    CLIP_IN_ADD(&(mv_ref_list)[(refmv_count)].as_mv, (bw), (bh), (xd));      \
    if (refmv_count && (mv_ref_list)[1].as_int != (mv_ref_list)[0].as_int) { \
      (refmv_count) = 2;                                                     \
      goto Done;                                                             \
    }                                                                        \
    (refmv_count) = 1;                                                       \
  } while (0)

// If either reference frame is different, not INTRA, and they
// are different from each other scale and add the mv to our list.
#define IF_DIFF_REF_FRAME_ADD_MV(mbmi, ref_frame, ref_sign_bias, refmv_count, \
                                 mv_ref_list, bw, bh, xd, Done)               \
  do {                                                                        \
    if (is_inter_block(mbmi)) {                                               \
      if ((mbmi)->ref_frame[0] != ref_frame)                                  \
        ADD_MV_REF_LIST(scale_mv((mbmi), 0, ref_frame, ref_sign_bias),        \
                        refmv_count, mv_ref_list, bw, bh, xd, Done);          \
      if (has_second_ref(mbmi) && (mbmi)->ref_frame[1] != ref_frame)          \
        ADD_MV_REF_LIST(scale_mv((mbmi), 1, ref_frame, ref_sign_bias),        \
                        refmv_count, mv_ref_list, bw, bh, xd, Done);          \
    }                                                                         \
  } while (0)

// Checks that the given mi_row, mi_col and search point
// are inside the borders of the tile.
static INLINE int is_inside(const TileInfo *const tile, int mi_col, int mi_row,
                            const POSITION *mi_pos) {
  return !(mi_row + mi_pos->row < tile->mi_row_start ||
           mi_col + mi_pos->col < tile->mi_col_start ||
           mi_row + mi_pos->row >= tile->mi_row_end ||
           mi_col + mi_pos->col >= tile->mi_col_end);
}

static INLINE void lower_mv_precision(MV *mv, int allow_hp) {
  if (!allow_hp) {
    if (mv->row & 1) mv->row += (mv->row > 0 ? -1 : 1);
    if (mv->col & 1) mv->col += (mv->col > 0 ? -1 : 1);
  }
}

#if CONFIG_REF_MV
static INLINE int av1_nmv_ctx(const uint8_t ref_mv_count,
                              const CANDIDATE_MV *ref_mv_stack, int ref,
                              int ref_mv_idx) {
  int_mv this_mv = (ref == 0) ? ref_mv_stack[ref_mv_idx].this_mv
                              : ref_mv_stack[ref_mv_idx].comp_mv;
#if CONFIG_EXT_INTER
  return 0;
#endif

  if (ref_mv_stack[ref_mv_idx].weight >= REF_CAT_LEVEL && ref_mv_count > 0) {
    if (abs(this_mv.as_mv.row -
            ref_mv_stack[ref_mv_idx].pred_mv[ref].as_mv.row) <= 4 &&
        abs(this_mv.as_mv.col -
            ref_mv_stack[ref_mv_idx].pred_mv[ref].as_mv.col) <= 4)
      return 2;
    else
      return 1;
  }

  return 0;
}

static INLINE int8_t av1_ref_frame_type(const MV_REFERENCE_FRAME *const rf) {
  if (rf[1] > INTRA_FRAME) {
    return TOTAL_REFS_PER_FRAME + FWD_RF_OFFSET(rf[0]) +
           BWD_RF_OFFSET(rf[1]) * FWD_REFS;
  }

  return rf[0];
}

static MV_REFERENCE_FRAME ref_frame_map[COMP_REFS][2] = {
#if CONFIG_EXT_REFS
  { LAST_FRAME, BWDREF_FRAME },  { LAST2_FRAME, BWDREF_FRAME },
  { LAST3_FRAME, BWDREF_FRAME }, { GOLDEN_FRAME, BWDREF_FRAME },

  { LAST_FRAME, ALTREF_FRAME },  { LAST2_FRAME, ALTREF_FRAME },
  { LAST3_FRAME, ALTREF_FRAME }, { GOLDEN_FRAME, ALTREF_FRAME }
#else
  { LAST_FRAME, ALTREF_FRAME }, { GOLDEN_FRAME, ALTREF_FRAME }
#endif
};

static INLINE void av1_set_ref_frame(MV_REFERENCE_FRAME *rf,
                                     int8_t ref_frame_type) {
  if (ref_frame_type >= TOTAL_REFS_PER_FRAME) {
    rf[0] = ref_frame_map[ref_frame_type - TOTAL_REFS_PER_FRAME][0];
    rf[1] = ref_frame_map[ref_frame_type - TOTAL_REFS_PER_FRAME][1];
  } else {
    rf[0] = ref_frame_type;
    rf[1] = NONE;
    assert(ref_frame_type > INTRA_FRAME &&
           ref_frame_type < TOTAL_REFS_PER_FRAME);
  }
}

static INLINE int16_t av1_mode_context_analyzer(
    const int16_t *const mode_context, const MV_REFERENCE_FRAME *const rf,
    BLOCK_SIZE bsize, int block) {
  int16_t mode_ctx = 0;
  int8_t ref_frame_type = av1_ref_frame_type(rf);

  if (block >= 0) {
    mode_ctx = mode_context[rf[0]] & 0x00ff;

    if (block > 0 && bsize < BLOCK_8X8 && bsize > BLOCK_4X4)
      mode_ctx |= (1 << SKIP_NEARESTMV_SUB8X8_OFFSET);

    return mode_ctx;
  }

  if (rf[1] > INTRA_FRAME)
    return mode_context[rf[0]] & (mode_context[rf[1]] | 0x00ff);
  else if (rf[0] != ALTREF_FRAME)
    return mode_context[rf[0]] & ~(mode_context[ALTREF_FRAME] & 0xfe00);
  else
    return mode_context[ref_frame_type];
}

static INLINE uint8_t av1_drl_ctx(const CANDIDATE_MV *ref_mv_stack,
                                  int ref_idx) {
  if (ref_mv_stack[ref_idx].weight >= REF_CAT_LEVEL &&
      ref_mv_stack[ref_idx + 1].weight >= REF_CAT_LEVEL)
    return 0;

  if (ref_mv_stack[ref_idx].weight >= REF_CAT_LEVEL &&
      ref_mv_stack[ref_idx + 1].weight < REF_CAT_LEVEL)
    return 2;

  if (ref_mv_stack[ref_idx].weight < REF_CAT_LEVEL &&
      ref_mv_stack[ref_idx + 1].weight < REF_CAT_LEVEL)
    return 3;

  return 0;
}
#endif

typedef void (*find_mv_refs_sync)(void *const data, int mi_row);
void av1_find_mv_refs(const AV1_COMMON *cm, const MACROBLOCKD *xd,
                      MODE_INFO *mi, MV_REFERENCE_FRAME ref_frame,
#if CONFIG_REF_MV
                      uint8_t *ref_mv_count, CANDIDATE_MV *ref_mv_stack,
#if CONFIG_EXT_INTER
                      int16_t *compound_mode_context,
#endif  // CONFIG_EXT_INTER
#endif
                      int_mv *mv_ref_list, int mi_row, int mi_col,
                      find_mv_refs_sync sync, void *const data,
                      int16_t *mode_context);

// check a list of motion vectors by sad score using a number rows of pixels
// above and a number cols of pixels in the left to select the one with best
// score to use as ref motion vector
void av1_find_best_ref_mvs(int allow_hp, int_mv *mvlist, int_mv *nearest_mv,
                           int_mv *near_mv);

void av1_append_sub8x8_mvs_for_idx(const AV1_COMMON *cm, MACROBLOCKD *xd,
                                   int block, int ref, int mi_row, int mi_col,
#if CONFIG_REF_MV
                                   CANDIDATE_MV *ref_mv_stack,
                                   uint8_t *ref_mv_count,
#endif
#if CONFIG_EXT_INTER
                                   int_mv *mv_list,
#endif  // CONFIG_EXT_INTER
                                   int_mv *nearest_mv, int_mv *near_mv);

#if CONFIG_EXT_INTER
// This function keeps a mode count for a given MB/SB
void av1_update_mv_context(const MACROBLOCKD *xd, MODE_INFO *mi,
                           MV_REFERENCE_FRAME ref_frame, int_mv *mv_ref_list,
                           int block, int mi_row, int mi_col,
                           int16_t *mode_context);
#endif  // CONFIG_EXT_INTER

#ifdef __cplusplus
}  // extern "C"
#endif

#endif  // AV1_COMMON_MVREF_COMMON_H_
