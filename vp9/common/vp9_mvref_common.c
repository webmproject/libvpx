/*
 *  Copyright (c) 2012 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include "vp9/common/vp9_mvref_common.h"

#if CONFIG_NEW_INTER && CONFIG_NEWMVREF
// This function returns either the appropriate subblock or block's mv,
// depending on whether block_size < 8x8 for both current block and the
// examined candidate block.
static int_mv get_subblock_mv(const MODE_INFO *candidate,
                              const MODE_INFO *current,
                              int curr_blk_idx, int ref,
                              int search_row, int search_col) {
  int candidate_type = candidate->mbmi.sb_type;

  if (curr_blk_idx >= 0 && candidate_type < BLOCK_8X8) {
    int candidate_blk_idx = 0;
    assert(current->mbmi.sb_type < BLOCK_8X8);

    // Both current block and the candidate block are in sub8x8 mode
    if ((search_row == -1 && search_col == 0) ||  // top
        (search_row == 0 && search_col == -1)) {  // left
      int i = curr_blk_idx + current->mbmi.sb_type * 4;
      int j = (search_row == 0);  // top: 0; left: 1

      candidate_blk_idx = idx_to_subblock_top_left[i][j][candidate_type];
      return (candidate_blk_idx >= 0) ?
          candidate->bmi[candidate_blk_idx].as_mv[ref] :
          candidate->mbmi.mv[ref];
    } else if ((search_row == -1 && search_col ==  1) ||  // top_right
               (search_row == -1 && search_col == -1)) {  // top_left
      candidate_blk_idx =
          idx_to_subblock_topright_topleft[search_col == -1][candidate_type];
      return candidate->bmi[candidate_blk_idx].as_mv[ref];
    }
  }

  return candidate->mbmi.mv[ref];
}

static int get_mvref_zone_idx(const TileInfo *const tile, int bsize,
                              int mi_row, int mi_col) {
  int mvref_zone_idx = 0;
  int row_8x8 = mi_row % 8;
  int col_8x8 = mi_col % 8;

  switch (bsize) {
    case BLOCK_4X4:
    case BLOCK_4X8:
    case BLOCK_8X4:
    case BLOCK_8X8:
      mvref_zone_idx =
          (mi_col >= (tile->mi_col_end - 1) ||  // right-most column
           (mv_ref_topright_avail_8x8[row_8x8][col_8x8] == 0)) ? 1 : 0;
      break;
    default:
      // Only <= BLOCK_8X8 are supported currently
      assert(0);
      break;
  }
  return mvref_zone_idx;
}

// This function searches the neighbourhood of a given MB/SB
// to try to find candidate reference vectors.
static void find_mv_refs_idx_8x8(const VP9_COMMON *cm, const MACROBLOCKD *xd,
                                 const TileInfo *const tile,
                                 MODE_INFO *mi, MV_REFERENCE_FRAME ref_frame,
                                 int_mv *mv_ref_list,
                                 int block, int mi_row, int mi_col) {
  int_mv mv_ref_candidates[MAX_MV_REF_CANDIDATES + 1];
  const int *ref_sign_bias = cm->ref_frame_sign_bias;
  int i;
  int refmv_count = 0;
  int different_ref_found = 0;

  int zone_idx = get_mvref_zone_idx(tile, mi->mbmi.sb_type, mi_row, mi_col);
  int max_nearest_blks = (zone_idx == 0) ? 4 : 3;
  const POSITION *mv_ref_search = mv_ref_blocks_8x8[zone_idx];

  // Zero out the mv reference vector list
  vpx_memset(mv_ref_list, 0, sizeof(*mv_ref_list) * MAX_MV_REF_CANDIDATES);
  vpx_memset(mv_ref_candidates, 0,
             sizeof(*mv_ref_candidates) * (MAX_MV_REF_CANDIDATES + 1));

  // The nearest 4 (when top right is available) or 3 neighboring blocks
  // are treated differently:
  //   If their block size < 8x8, we get the mv from the bmi substructure.
  for (i = 0; i < max_nearest_blks; ++i) {
    const POSITION *const mv_ref = &mv_ref_search[i];
    if (is_inside(tile, mi_col, mi_row, cm->mi_rows, mv_ref)) {
      const MODE_INFO *const candidate_mi =
          xd->mi[mv_ref->col + mv_ref->row * xd->mi_stride].src_mi;
      const MB_MODE_INFO *const candidate = &candidate_mi->mbmi;

      different_ref_found = 1;

      if (candidate->ref_frame[0] == ref_frame) {
        ADD_MV_REF_CANDIDATE(get_subblock_mv(
            candidate_mi, mi, block, 0, mv_ref->row, mv_ref->col));
      } else if (candidate->ref_frame[1] == ref_frame) {
        ADD_MV_REF_CANDIDATE(get_subblock_mv(
            candidate_mi, mi, block, 1, mv_ref->row, mv_ref->col));
      }
    }
  }

  // Check the rest of the neighbors in much the same way as before
  // except we don't need to keep track of subblocks.
  for (; i < MVREF_NEIGHBOURS; ++i) {
    const POSITION *const mv_ref = &mv_ref_search[i];
    if (is_inside(tile, mi_col, mi_row, cm->mi_rows, mv_ref)) {
      const MB_MODE_INFO *const candidate =
          &xd->mi[mv_ref->col + mv_ref->row * xd->mi_stride].src_mi->mbmi;

      different_ref_found = 1;

      if (candidate->ref_frame[0] == ref_frame)
        ADD_MV_REF_CANDIDATE(candidate->mv[0]);
      else if (candidate->ref_frame[1] == ref_frame)
        ADD_MV_REF_CANDIDATE(candidate->mv[1]);
    }
  }

  // Since we couldn't find 3 mvs from the same reference frame,
  // go back through the neighbors and find motion vectors from
  // different reference frames.
  if (different_ref_found) {
    for (i = 0; i < MVREF_NEIGHBOURS; ++i) {
      const POSITION *mv_ref = &mv_ref_search[i];
      if (is_inside(tile, mi_col, mi_row, cm->mi_rows, mv_ref)) {
        const MB_MODE_INFO *const candidate =
            &xd->mi[mv_ref->col + mv_ref->row * xd->mi_stride].src_mi->mbmi;

        // If the candidate is INTRA we don't want to consider its mv.
        IF_DIFF_REF_FRAME_ADD_MV_CANDIDATE(candidate);
      }
    }
  }

 Done:

  if (refmv_count == 2) {
    mv_ref_list[0].as_mv.row =
        (mv_ref_candidates[0].as_mv.row + mv_ref_candidates[1].as_mv.row) >> 1;
    mv_ref_list[0].as_mv.col =
        (mv_ref_candidates[0].as_mv.col + mv_ref_candidates[1].as_mv.col) >> 1;
    mv_ref_list[1].as_int = mv_ref_candidates[2].as_int;
  } else {
    for (i = 0; i < 2; ++i) {
      mv_ref_list[i].as_int = mv_ref_candidates[i].as_int;
    }
  }

  // Clamp vectors
  for (i = 0; i < MAX_MV_REF_CANDIDATES; ++i)
    clamp_mv_ref(&mv_ref_list[i].as_mv, xd);
}

typedef enum MV_SEARCH_POS {
  TOP            = 0,
  LEFT           = 1,
  TOPLEFT        = 2,
  TOPRIGHT       = 3,
  TOPRIGHT_ALT   = 4,
  NUM_SEARCH_POS = 5
} MV_SEARCH_POS;

// Adaptive median
static int get_adaptive_median(int topright, int left, int topleft) {
  int a = topright;
  int b = left;
  int c = topright + left - topleft;

  if (a >= b) {
    if (b >= c)      return b;
    else if (a >= c) return c;
    else             return a;
  } else {
    if (b < c)       return b;
    else if (a >= c) return a;
    else             return c;
  }
}

// This function searches the neighbourhood of a given MB/SB to try
// to find the nearestmv through adaptive median filtering.
static int find_best_mvref_8x8(const VP9_COMMON *cm, const MACROBLOCKD *xd,
                               const TileInfo *const tile,
                               MODE_INFO *mi, MV_REFERENCE_FRAME ref_frame,
                               int_mv *best_mvref,
                               int block, int mi_row, int mi_col) {
  int i;
  int zone_idx = get_mvref_zone_idx(tile, mi->mbmi.sb_type, mi_row, mi_col);
  int max_nearest_blks = (zone_idx == 0) ? 4 : 3;

  const POSITION adapt_median_neighbor_pos[NUM_SEARCH_POS - 1] = {
    // TOP, LEFT, TOPLEFT, TOPRIGHT
    {-1, 0}, {0, -1}, {-1, -1}, {-1, 1}
  };
  int_mv mv_ref_mvs[NUM_SEARCH_POS];
  int is_avail[NUM_SEARCH_POS] = { 0, 0, 0, 0, 0 };

  vpx_memset(mv_ref_mvs, 0, sizeof(mv_ref_mvs[0]) * NUM_SEARCH_POS);

  // If the neighboring block size < 8x8, the mv is obtained from
  // the bmi substructure.
  for (i = 0; i < max_nearest_blks; ++i) {
    const POSITION *const mv_ref_pos = &adapt_median_neighbor_pos[i];

    if (is_inside(tile, mi_col, mi_row, cm->mi_rows, mv_ref_pos)) {
      const MODE_INFO *const candidate_mi =
          xd->mi[mv_ref_pos->col + mv_ref_pos->row * xd->mi_stride].src_mi;
      const MB_MODE_INFO *const candidate = &candidate_mi->mbmi;

      if (candidate->ref_frame[0] == ref_frame) {
        mv_ref_mvs[i] = get_subblock_mv(candidate_mi, mi, block, 0,
                                        mv_ref_pos->row, mv_ref_pos->col);
        is_avail[i] = 1;
      } else if (candidate->ref_frame[1] == ref_frame) {
        mv_ref_mvs[i] = get_subblock_mv(candidate_mi, mi, block, 1,
                                        mv_ref_pos->row, mv_ref_pos->col);
        is_avail[i] = 1;
      }
    }
  }

  if (is_avail[TOP] && is_avail[TOPRIGHT]) {
    mv_ref_mvs[TOPRIGHT_ALT].as_mv.row =
        (mv_ref_mvs[TOP].as_mv.row + mv_ref_mvs[TOPRIGHT].as_mv.row) >> 1;
    mv_ref_mvs[TOPRIGHT_ALT].as_mv.col =
        (mv_ref_mvs[TOP].as_mv.col + mv_ref_mvs[TOPRIGHT].as_mv.col) >> 1;
  } else if (is_avail[TOP]) {
    mv_ref_mvs[TOPRIGHT_ALT].as_int = mv_ref_mvs[TOP].as_int;
  } else if (is_avail[TOPRIGHT]) {
    mv_ref_mvs[TOPRIGHT_ALT].as_int = mv_ref_mvs[TOPRIGHT].as_int;
  }

  if (is_avail[TOP] || is_avail[LEFT] || is_avail[TOPLEFT] ||
      is_avail[TOPRIGHT]) {
    best_mvref->as_mv.row = get_adaptive_median(
        mv_ref_mvs[TOPRIGHT_ALT].as_mv.row,
        mv_ref_mvs[LEFT].as_mv.row,
        mv_ref_mvs[TOPLEFT].as_mv.row);
    best_mvref->as_mv.col = get_adaptive_median(
        mv_ref_mvs[TOPRIGHT_ALT].as_mv.col,
        mv_ref_mvs[LEFT].as_mv.col,
        mv_ref_mvs[TOPLEFT].as_mv.col);

    // Clamp vectors
    clamp_mv_ref(&(best_mvref->as_mv), xd);
    return 1;
  } else {
    best_mvref->as_int = 0;
    return 0;
  }
}
#endif  // CONFIG_NEW_INTER && CONFIG_NEWMVREF

// This function searches the neighbourhood of a given MB/SB
// to try and find candidate reference vectors.
static void find_mv_refs_idx(const VP9_COMMON *cm, const MACROBLOCKD *xd,
                             const TileInfo *const tile,
                             MODE_INFO *mi, MV_REFERENCE_FRAME ref_frame,
                             int_mv *mv_ref_list,
                             int block, int mi_row, int mi_col) {
  const int *ref_sign_bias = cm->ref_frame_sign_bias;
  int i, refmv_count = 0;
  const MODE_INFO *prev_mi = !cm->error_resilient_mode && cm->prev_mi
        ? cm->prev_mi[mi_row * xd->mi_stride + mi_col].src_mi
        : NULL;
  const MB_MODE_INFO *const prev_mbmi = prev_mi ? &prev_mi->mbmi : NULL;
  const POSITION *const mv_ref_search = mv_ref_blocks[mi->mbmi.sb_type];
  int different_ref_found = 0;
#if !CONFIG_NEW_INTER
  int context_counter = 0;
#endif  // !CONFIG_NEW_INTER

  // Blank the reference vector list
  vpx_memset(mv_ref_list, 0, sizeof(*mv_ref_list) * MAX_MV_REF_CANDIDATES);

  // The nearest 2 blocks are treated differently
  // if the size < 8x8 we get the mv from the bmi substructure,
  // and we also need to keep a mode count.
  for (i = 0; i < 2; ++i) {
    const POSITION *const mv_ref = &mv_ref_search[i];
    if (is_inside(tile, mi_col, mi_row, cm->mi_rows, mv_ref)) {
      const MODE_INFO *const candidate_mi = xd->mi[mv_ref->col + mv_ref->row *
                                                   xd->mi_stride].src_mi;
      const MB_MODE_INFO *const candidate = &candidate_mi->mbmi;
      // Keep counts for entropy encoding.
#if !CONFIG_NEW_INTER
      context_counter += mode_2_counter[candidate->mode];
#endif  // !CONFIG_NEW_INTER
      different_ref_found = 1;

      if (candidate->ref_frame[0] == ref_frame) {
        ADD_MV_REF_LIST(get_sub_block_mv(candidate_mi, 0, mv_ref->col, block));
      } else if (candidate->ref_frame[1] == ref_frame) {
        ADD_MV_REF_LIST(get_sub_block_mv(candidate_mi, 1, mv_ref->col, block));
      }
    }
  }

  // Check the rest of the neighbors in much the same way
  // as before except we don't need to keep track of sub blocks or
  // mode counts.
  for (; i < MVREF_NEIGHBOURS; ++i) {
    const POSITION *const mv_ref = &mv_ref_search[i];
    if (is_inside(tile, mi_col, mi_row, cm->mi_rows, mv_ref)) {
      const MB_MODE_INFO *const candidate = &xd->mi[mv_ref->col + mv_ref->row *
                                                    xd->mi_stride].src_mi->mbmi;
      different_ref_found = 1;

      if (candidate->ref_frame[0] == ref_frame)
        ADD_MV_REF_LIST(candidate->mv[0]);
      else if (candidate->ref_frame[1] == ref_frame)
        ADD_MV_REF_LIST(candidate->mv[1]);
    }
  }

  // Check the last frame's mode and mv info.
  if (prev_mbmi) {
    if (prev_mbmi->ref_frame[0] == ref_frame)
      ADD_MV_REF_LIST(prev_mbmi->mv[0]);
    else if (prev_mbmi->ref_frame[1] == ref_frame)
      ADD_MV_REF_LIST(prev_mbmi->mv[1]);
  }

  // Since we couldn't find 2 mvs from the same reference frame
  // go back through the neighbors and find motion vectors from
  // different reference frames.
  if (different_ref_found) {
    for (i = 0; i < MVREF_NEIGHBOURS; ++i) {
      const POSITION *mv_ref = &mv_ref_search[i];
      if (is_inside(tile, mi_col, mi_row, cm->mi_rows, mv_ref)) {
        const MB_MODE_INFO *const candidate = &xd->mi[mv_ref->col + mv_ref->row
                                              * xd->mi_stride].src_mi->mbmi;

        // If the candidate is INTRA we don't want to consider its mv.
        IF_DIFF_REF_FRAME_ADD_MV(candidate);
      }
    }
  }

  // Since we still don't have a candidate we'll try the last frame.
  if (prev_mbmi)
    IF_DIFF_REF_FRAME_ADD_MV(prev_mbmi);

 Done:

#if !CONFIG_NEW_INTER
  mi->mbmi.mode_context[ref_frame] = counter_to_context[context_counter];
#endif  // !CONFIG_NEW_INTER

  // Clamp vectors
  for (i = 0; i < MAX_MV_REF_CANDIDATES; ++i)
    clamp_mv_ref(&mv_ref_list[i].as_mv, xd);
}

#if CONFIG_NEW_INTER
// This function keeps a mode count for a given MB/SB
void vp9_update_mv_context(const VP9_COMMON *cm, const MACROBLOCKD *xd,
                           const TileInfo *const tile,
                           MODE_INFO *mi, MV_REFERENCE_FRAME ref_frame,
                           int_mv *mv_ref_list,
                           int block, int mi_row, int mi_col) {
  int i, refmv_count = 0;
  const POSITION *const mv_ref_search = mv_ref_blocks[mi->mbmi.sb_type];
  int context_counter = 0;

  // Blank the reference vector list
  vpx_memset(mv_ref_list, 0, sizeof(*mv_ref_list) * MAX_MV_REF_CANDIDATES);

  // The nearest 2 blocks are examined only.
  // If the size < 8x8, we get the mv from the bmi substructure;
  for (i = 0; i < 2; ++i) {
    const POSITION *const mv_ref = &mv_ref_search[i];
    if (is_inside(tile, mi_col, mi_row, cm->mi_rows, mv_ref)) {
      const MODE_INFO *const candidate_mi =
          xd->mi[mv_ref->col + mv_ref->row * xd->mi_stride].src_mi;
      const MB_MODE_INFO *const candidate = &candidate_mi->mbmi;

      // Keep counts for entropy encoding.
      context_counter += mode_2_counter[candidate->mode];

      if (candidate->ref_frame[0] == ref_frame) {
        ADD_MV_REF_LIST(get_sub_block_mv(candidate_mi, 0, mv_ref->col, block));
      } else if (candidate->ref_frame[1] == ref_frame) {
        ADD_MV_REF_LIST(get_sub_block_mv(candidate_mi, 1, mv_ref->col, block));
      }
    }
  }

 Done:

  mi->mbmi.mode_context[ref_frame] = counter_to_context[context_counter];
}
#endif  // CONFIG_NEW_INTER

void vp9_find_mv_refs(const VP9_COMMON *cm, const MACROBLOCKD *xd,
                      const TileInfo *const tile,
                      MODE_INFO *mi, MV_REFERENCE_FRAME ref_frame,
                      int_mv *mv_ref_list,
                      int mi_row, int mi_col) {
#if CONFIG_NEW_INTER
  vp9_update_mv_context(cm, xd, tile, mi, ref_frame, mv_ref_list, -1,
                        mi_row, mi_col);
#if CONFIG_NEWMVREF
  if (mi->mbmi.sb_type <= BLOCK_8X8) {
    int_mv best_mvref;
    find_best_mvref_8x8(cm, xd, tile, mi, ref_frame, &best_mvref,
                        -1, mi_row, mi_col);
    find_mv_refs_idx_8x8(cm, xd, tile, mi, ref_frame, mv_ref_list,
                         -1, mi_row, mi_col);
    if (best_mvref.as_int != 0) {
      mv_ref_list[1].as_int = mv_ref_list[0].as_int;
      mv_ref_list[0].as_int = best_mvref.as_int;
    }
  } else {
#endif  // CONFIG_NEWMVREF
#endif  // CONFIG_NEW_INTER
  find_mv_refs_idx(cm, xd, tile, mi, ref_frame, mv_ref_list, -1,
                   mi_row, mi_col);
#if CONFIG_NEW_INTER && CONFIG_NEWMVREF
  }
#endif  // CONFIG_NEW_INTER && CONFIG_NEWMVREF
}

void vp9_find_best_ref_mvs(MACROBLOCKD *xd, int allow_hp,
                           int_mv *mvlist, int_mv *nearest, int_mv *near) {
  int i;
  // Make sure all the candidates are properly clamped etc
  for (i = 0; i < MAX_MV_REF_CANDIDATES; ++i) {
    MV *mv = &mvlist[i].as_mv;
    const int usehp = allow_hp && vp9_use_mv_hp(mv);
    vp9_lower_mv_precision(mv, usehp);
    clamp_mv2(mv, xd);
  }
  *nearest = mvlist[0];
  *near = mvlist[1];
}

void vp9_append_sub8x8_mvs_for_idx(VP9_COMMON *cm, MACROBLOCKD *xd,
                                   const TileInfo *const tile,
                                   int block, int ref, int mi_row, int mi_col,
#if CONFIG_NEW_INTER
                                   int_mv *mv_list,
#endif  // CONFIG_NEW_INTER
                                   int_mv *nearest, int_mv *near) {
#if CONFIG_NEW_INTER
#if CONFIG_NEWMVREF
  int_mv best_mvref;
#endif  // CONFIG_NEWMVREF
#else
  int_mv mv_list[MAX_MV_REF_CANDIDATES];
#endif  // !CONFIG_NEW_INTER
  MODE_INFO *const mi = xd->mi[0].src_mi;
  b_mode_info *bmi = mi->bmi;
  int n;

  assert(MAX_MV_REF_CANDIDATES == 2);

#if CONFIG_NEW_INTER && CONFIG_NEWMVREF
  find_best_mvref_8x8(cm, xd, tile, mi, mi->mbmi.ref_frame[ref],
                      &best_mvref, block, mi_row, mi_col);
  find_mv_refs_idx_8x8(cm, xd, tile, mi, mi->mbmi.ref_frame[ref],
                       mv_list, block, mi_row, mi_col);
#else
  find_mv_refs_idx(cm, xd, tile, mi, mi->mbmi.ref_frame[ref],
                   mv_list, block, mi_row, mi_col);
#endif  // CONFIG_NEW_INTER && CONFIG_NEWMVREF

  near->as_int = 0;

  switch (block) {
    case 0:
#if CONFIG_NEW_INTER && CONFIG_NEWMVREF
      if (best_mvref.as_int != 0) {
        nearest->as_int = best_mvref.as_int;
        if (best_mvref.as_int != mv_list[0].as_int)
          near->as_int = mv_list[0].as_int;
        else
          near->as_int = mv_list[1].as_int;
      } else {
#endif  // CONFIG_NEW_INTER && CONFIG_NEWMVREF
        nearest->as_int = mv_list[0].as_int;
        near->as_int = mv_list[1].as_int;
#if CONFIG_NEW_INTER && CONFIG_NEWMVREF
      }
#endif  // CONFIG_NEW_INTER && CONFIG_NEWMVREF
      break;
    case 1:
#if !CONFIG_NEW_INTER
    case 2:
#endif  // !CONFIG_NEW_INTER
      nearest->as_int = bmi[0].as_mv[ref].as_int;
#if CONFIG_NEW_INTER && CONFIG_NEWMVREF
      if (best_mvref.as_int != 0 &&
          best_mvref.as_int != nearest->as_int)
        near->as_int = best_mvref.as_int;
      else
#endif  // CONFIG_NEW_INTER && CONFIG_NEWMVREF
      for (n = 0; n < MAX_MV_REF_CANDIDATES; ++n)
        if (nearest->as_int != mv_list[n].as_int) {
          near->as_int = mv_list[n].as_int;
          break;
        }
      break;
#if CONFIG_NEW_INTER
    case 2: {
#if CONFIG_NEWMVREF
      if (bmi[0].as_mv[ref].as_int !=
          bmi[1].as_mv[ref].as_int) {
        // Average of TOP and TOPRIGHT
        nearest->as_mv.row = (
            bmi[0].as_mv[ref].as_mv.row +
            bmi[1].as_mv[ref].as_mv.row) >> 1;
        nearest->as_mv.col = (
            bmi[0].as_mv[ref].as_mv.col +
            bmi[1].as_mv[ref].as_mv.col) >> 1;
        near->as_int = bmi[0].as_mv[ref].as_int;
      } else {
        nearest->as_int = bmi[0].as_mv[ref].as_int;
        if (best_mvref.as_int != 0 &&
            best_mvref.as_int != nearest->as_int) {
          near->as_int = best_mvref.as_int;
        } else {
          for (n = 0; n < MAX_MV_REF_CANDIDATES; ++n)
            if (nearest->as_int != mv_list[n].as_int) {
              near->as_int = mv_list[n].as_int;
              break;
            }
        }
      }
#else
      int_mv candidates[1 + MAX_MV_REF_CANDIDATES];
      candidates[0] = bmi[1].as_mv[ref];
      candidates[1] = mv_list[0];
      candidates[2] = mv_list[1];

      nearest->as_int = bmi[0].as_mv[ref].as_int;
      for (n = 0; n < 1 + MAX_MV_REF_CANDIDATES; ++n)
        if (nearest->as_int != candidates[n].as_int) {
          near->as_int = candidates[n].as_int;
          break;
        }
#endif  // CONFIG_NEWMVREF
      break;
    }
#endif  // CONFIG_NEW_INTER
    case 3: {
#if CONFIG_NEW_INTER && CONFIG_NEWMVREF
      if (bmi[0].as_mv[ref].as_int != bmi[1].as_mv[ref].as_int ||
          bmi[0].as_mv[ref].as_int != bmi[2].as_mv[ref].as_int ||
          bmi[1].as_mv[ref].as_int != bmi[2].as_mv[ref].as_int) {
        nearest->as_mv.row = get_adaptive_median(
            bmi[1].as_mv[ref].as_mv.row,
            bmi[2].as_mv[ref].as_mv.row,
            bmi[0].as_mv[ref].as_mv.row);
        nearest->as_mv.col = get_adaptive_median(
            bmi[1].as_mv[ref].as_mv.col,
            bmi[2].as_mv[ref].as_mv.col,
            bmi[0].as_mv[ref].as_mv.col);
        /*nearest->as_mv.row =
            (bmi[0].as_mv[ref].as_mv.row +
             bmi[1].as_mv[ref].as_mv.row +
             (bmi[2].as_mv[ref].as_mv.row << 1)) >> 2;
        nearest->as_mv.col =
            (bmi[0].as_mv[ref].as_mv.col +
             bmi[1].as_mv[ref].as_mv.col +
             (bmi[2].as_mv[ref].as_mv.col << 1)) >> 2;*/
        for (n = 2; n >= 0; --n)
          if (nearest->as_int != bmi[n].as_mv[ref].as_int) {
            near->as_int = bmi[n].as_mv[ref].as_int;
            break;
          }
      } else {
        nearest->as_int = bmi[2].as_mv[ref].as_int;
        if (best_mvref.as_int != 0 &&
            best_mvref.as_int != nearest->as_int) {
          near->as_int = best_mvref.as_int;
        } else {
          for (n = 0; n < MAX_MV_REF_CANDIDATES; ++n)
            if (nearest->as_int != mv_list[n].as_int) {
              near->as_int = mv_list[n].as_int;
              break;
          }
        }
      }
#else
      int_mv candidates[2 + MAX_MV_REF_CANDIDATES];
      candidates[0] = bmi[1].as_mv[ref];
      candidates[1] = bmi[0].as_mv[ref];
      candidates[2] = mv_list[0];
      candidates[3] = mv_list[1];

      nearest->as_int = bmi[2].as_mv[ref].as_int;
      for (n = 0; n < 2 + MAX_MV_REF_CANDIDATES; ++n)
        if (nearest->as_int != candidates[n].as_int) {
          near->as_int = candidates[n].as_int;
          break;
        }
#endif  // CONFIG_NEW_INTER && CONFIG_NEWMREF
      break;
    }
    default:
      assert(0 && "Invalid block index.");
  }
}

#if CONFIG_COPY_MODE
static int compare_interinfo(MB_MODE_INFO *mbmi, MB_MODE_INFO *ref_mbmi) {
  if (mbmi == ref_mbmi) {
    return 1;
  } else {
    int is_same;
#if CONFIG_INTERINTRA
    MV_REFERENCE_FRAME mbmi_ref1_backup = mbmi->ref_frame[1];
    MV_REFERENCE_FRAME refmbmi_ref1_backup = ref_mbmi->ref_frame[1];

    if (mbmi->ref_frame[1] == INTRA_FRAME)
      mbmi->ref_frame[1] = NONE;
    if (ref_mbmi->ref_frame[1] == INTRA_FRAME)
      ref_mbmi->ref_frame[1] = NONE;
#endif  // CONFIG_INTERINTRA
    if (mbmi->ref_frame[0] == ref_mbmi->ref_frame[0] &&
        mbmi->ref_frame[1] == ref_mbmi->ref_frame[1]) {
      if (mbmi->ref_frame[1] > INTRA_FRAME)
        is_same = mbmi->mv[0].as_int == ref_mbmi->mv[0].as_int &&
                  mbmi->mv[1].as_int == ref_mbmi->mv[1].as_int &&
                  mbmi->interp_filter == ref_mbmi->interp_filter;
      else
        is_same = mbmi->mv[0].as_int == ref_mbmi->mv[0].as_int &&
                  mbmi->interp_filter == ref_mbmi->interp_filter;
    } else {
      is_same = 0;
    }
#if CONFIG_INTERINTRA
    mbmi->ref_frame[1] = mbmi_ref1_backup;
    ref_mbmi->ref_frame[1] = refmbmi_ref1_backup;
#endif  // CONFIG_INTERINTRA

    return is_same;
  }
}

static int check_inside(const TileInfo *const tile, int mi_row, int mi_col) {
  return mi_row >= tile->mi_row_start && mi_col >= tile->mi_col_start &&
         mi_row < tile->mi_row_end && mi_col < tile->mi_col_end;
}

static int is_right_available(BLOCK_SIZE bsize,
#if CONFIG_EXT_PARTITION
                              PARTITION_TYPE partition,
#endif
                              int mi_row, int mi_col) {
  int depth, max_depth = (CODING_UNIT_SIZE_LOG2 - 2) -
          MIN(b_width_log2_lookup[bsize], b_height_log2_lookup[bsize]);
  int block[(CODING_UNIT_SIZE_LOG2 - 2)] = {0};

  if (bsize == BLOCK_LARGEST)
    return 1;
  mi_row = mi_row % MI_BLOCK_SIZE;
  mi_col = mi_col % MI_BLOCK_SIZE;
  for (depth = 1; depth <= max_depth; depth++) {
    block[depth] = (mi_row >> (MI_BLOCK_SIZE_LOG2 - depth)) * 2 +
                   (mi_col >> (MI_BLOCK_SIZE_LOG2 - depth));
    mi_row = mi_row % (MI_BLOCK_SIZE >> depth);
    mi_col = mi_col % (MI_BLOCK_SIZE >> depth);
  }

  if (b_width_log2_lookup[bsize] < b_height_log2_lookup[bsize]) {
    if (block[max_depth] == 0)
      return 1;
  } else if (b_width_log2_lookup[bsize] > b_height_log2_lookup[bsize]) {
    if (block[max_depth] > 0)
      return 0;
  } else {
#if CONFIG_EXT_PARTITION
    if (block[max_depth] == 0)
      return 1;
    if (block[max_depth] == 2)
      return partition != PARTITION_VERT_A;
#else
    if (block[max_depth] == 0 || block[max_depth] == 2)
      return 1;
#endif
    else if (block[max_depth] == 3)
      return 0;
  }

  for (depth = max_depth - 1; depth > 0; depth--) {
    if (block[depth] == 0 || block[depth] == 2)
      return 1;
    else if (block[depth] == 3)
      return 0;
  }
  return 1;
}

static int is_second_rec(int mi_row, int mi_col, BLOCK_SIZE bsize) {
  int bw = 4 << b_width_log2_lookup[bsize];
  int bh = 4 << b_height_log2_lookup[bsize];

  if (bw < bh)
    return (mi_col << 3) % (bw << 1) == 0 ? 0 : 1;
  else if (bh < bw)
    return (mi_row << 3) % (bh << 1) == 0 ? 0 : 2;
  else
    return 0;
}

int vp9_construct_ref_inter_list(VP9_COMMON *cm,  MACROBLOCKD *xd,
                                 const TileInfo *const tile,
                                 BLOCK_SIZE bsize,
#if CONFIG_EXT_PARTITION
                                 PARTITION_TYPE partition,
#endif
                                 int mi_row, int mi_col,
                                 MB_MODE_INFO *ref_list[2 *
                                                        (MI_BLOCK_SIZE + 1)]) {
  int bw = 4 << b_width_log2_lookup[bsize];
  int bh = 4 << b_height_log2_lookup[bsize];
  int row_offset, col_offset;
  int mi_offset;
  MB_MODE_INFO *ref_mbmi;
  int ref_index, ref_num = 0;
  int row_offset_cand[2 * (MI_BLOCK_SIZE + 1)];
  int col_offset_cand[2 * (MI_BLOCK_SIZE + 1)];
  int offset_num = 0, i, switchflag;
  int is_sec_rec = is_second_rec(mi_row, mi_col, bsize);

  if (is_sec_rec != 2) {
    row_offset_cand[offset_num] = -1; col_offset_cand[offset_num] = 0;
    offset_num++;
  }
  if (is_sec_rec != 1) {
    row_offset_cand[offset_num] = bh / (2 * MI_SIZE);
    col_offset_cand[offset_num] = -1;
    offset_num++;
  }

  row_offset = bh / MI_SIZE - 1;
  col_offset = 1;
  if (is_sec_rec < 2)
    switchflag = 1;
  else
    switchflag = 0;
  while ((is_sec_rec == 0 && ((row_offset >=0) ||
                              col_offset < (bw / MI_SIZE + 1))) ||
         (is_sec_rec == 1 && col_offset < (bw / MI_SIZE + 1)) ||
         (is_sec_rec == 2 && row_offset >=0)) {
    switch (switchflag) {
      case 0:
        if (row_offset >= 0) {
          if (row_offset != bh / (2 * MI_SIZE)) {
            row_offset_cand[offset_num] = row_offset;
            col_offset_cand[offset_num] = -1;
            offset_num++;
          }
          row_offset--;
        }
        break;
      case 1:
        if (col_offset < (bw / MI_SIZE + 1)) {
          row_offset_cand[offset_num] = -1;
          col_offset_cand[offset_num] = col_offset;
          offset_num++;
          col_offset++;
        }
        break;
      default:
        assert(0);
    }
    if (is_sec_rec == 0)
      switchflag = 1 - switchflag;
  }
  row_offset_cand[offset_num] = -1;
  col_offset_cand[offset_num] = -1;
  offset_num++;

  for (i = 0; i < offset_num; i++) {
    row_offset = row_offset_cand[i];
    col_offset = col_offset_cand[i];
    if ((col_offset < (bw / MI_SIZE) ||
        (col_offset == (bw / MI_SIZE) && is_right_available(bsize,
#if CONFIG_EXT_PARTITION
                                                      partition,
#endif
                                                      mi_row, mi_col)))
        && check_inside(tile, mi_row + row_offset, mi_col + col_offset)) {
      mi_offset = row_offset * cm->mi_stride + col_offset;
      ref_mbmi = &xd->mi[mi_offset].src_mi->mbmi;
      if (is_inter_block(ref_mbmi)) {
        for (ref_index = 0; ref_index < ref_num; ref_index++) {
          if (compare_interinfo(ref_mbmi, ref_list[ref_index]))
            break;
        }
        if (ref_index == ref_num) {
          ref_list[ref_num] = ref_mbmi;
          ref_num++;
        }
      }
    }
  }
  return ref_num;
}
#endif  // CONFIG_COPY_MODE

#if CONFIG_INTRABC
void vp9_find_ref_dv(int_mv *ref_dv, int mi_row, int mi_col) {
  (void) mi_col;
  if (mi_row < 8) {
    ref_dv->as_mv.row = 0;
    ref_dv->as_mv.col = -8 * 8;
  } else {
    ref_dv->as_mv.row = -8 * 8;
    ref_dv->as_mv.col = 0;
  }
}
#endif  // CONFIG_INTRABC
