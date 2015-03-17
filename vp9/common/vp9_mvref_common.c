
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
#if !CONFIG_NEWMVREF
  int context_counter = 0;
#endif  // CONFIG_NEWMVREF

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
#if !CONFIG_NEWMVREF
      context_counter += mode_2_counter[candidate->mode];
#endif  // CONFIG_NEWMVREF
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

#if !CONFIG_NEWMVREF
  mi->mbmi.mode_context[ref_frame] = counter_to_context[context_counter];
#endif  // CONFIG_NEWMVREF

  // Clamp vectors
  for (i = 0; i < MAX_MV_REF_CANDIDATES; ++i)
    clamp_mv_ref(&mv_ref_list[i].as_mv, xd);
}

#if CONFIG_NEWMVREF
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
#endif  // CONFIG_NEWMVREF

void vp9_find_mv_refs(const VP9_COMMON *cm, const MACROBLOCKD *xd,
                      const TileInfo *const tile,
                      MODE_INFO *mi, MV_REFERENCE_FRAME ref_frame,
                      int_mv *mv_ref_list,
                      int mi_row, int mi_col) {
#if CONFIG_NEWMVREF
  vp9_update_mv_context(cm, xd, tile, mi, ref_frame, mv_ref_list, -1,
                        mi_row, mi_col);
#endif  // CONFIG_NEWMVREF
  find_mv_refs_idx(cm, xd, tile, mi, ref_frame, mv_ref_list, -1,
                   mi_row, mi_col);
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
#if CONFIG_NEWMVREF
                                   int_mv *mv_list,
#endif  // CONFIG_NEWMVREF
                                   int_mv *nearest, int_mv *near) {
#if !CONFIG_NEWMVREF
  int_mv mv_list[MAX_MV_REF_CANDIDATES];
#endif  // CONFIG_NEWMVREF
  MODE_INFO *const mi = xd->mi[0].src_mi;
  b_mode_info *bmi = mi->bmi;
  int n;

  assert(MAX_MV_REF_CANDIDATES == 2);

  find_mv_refs_idx(cm, xd, tile, mi, mi->mbmi.ref_frame[ref], mv_list, block,
                   mi_row, mi_col);

  near->as_int = 0;
  switch (block) {
    case 0:
      nearest->as_int = mv_list[0].as_int;
      near->as_int = mv_list[1].as_int;
      break;
    case 1:
    case 2:
      nearest->as_int = bmi[0].as_mv[ref].as_int;
      for (n = 0; n < MAX_MV_REF_CANDIDATES; ++n)
        if (nearest->as_int != mv_list[n].as_int) {
          near->as_int = mv_list[n].as_int;
          break;
        }
      break;
    case 3: {
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
      break;
    }
    default:
      assert("Invalid block index.");
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

static int check_inside(VP9_COMMON *cm, int mi_row, int mi_col) {
  return mi_row >= 0 && mi_col >= 0 &&
         mi_row < cm->mi_rows && mi_col < cm->mi_cols;
}

static int is_right_available(BLOCK_SIZE bsize, int mi_row, int mi_col) {
  int depth, max_depth = 4 - MIN(b_width_log2_lookup[bsize],
                                 b_height_log2_lookup[bsize]);
  int block[4] = {0};

  if (bsize == BLOCK_64X64)
    return 1;
  mi_row = mi_row % 8;
  mi_col = mi_col % 8;
  for (depth = 1; depth <= max_depth; depth++) {
    block[depth] = (mi_row >> (3 - depth)) * 2 + (mi_col >> (3 - depth));
    mi_row = mi_row % (8 >> depth);
    mi_col = mi_col % (8 >> depth);
  }

  if (b_width_log2_lookup[bsize] < b_height_log2_lookup[bsize]) {
    if (block[max_depth] == 0)
      return 1;
  } else if (b_width_log2_lookup[bsize] > b_height_log2_lookup[bsize]) {
    if (block[max_depth] > 0)
      return 0;
  } else {
    if (block[max_depth] == 0 || block[max_depth] == 2)
      return 1;
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
                                 BLOCK_SIZE bsize, int mi_row, int mi_col,
                                 MB_MODE_INFO *ref_list[18]) {
  int bw = 4 << b_width_log2_lookup[bsize];
  int bh = 4 << b_height_log2_lookup[bsize];
  int row_offset, col_offset;
  int mi_offset;
  MB_MODE_INFO *ref_mbmi;
  int ref_index, ref_num = 0;
  int row_offset_cand[18], col_offset_cand[18];
  int offset_num = 0, i, switchflag;
  int is_sec_rec = is_second_rec(mi_row, mi_col, bsize);

  if (is_sec_rec != 2) {
    row_offset_cand[offset_num] = -1; col_offset_cand[offset_num] = 0;
    offset_num++;
  }
  if (is_sec_rec != 1) {
    row_offset_cand[offset_num] = bh / 16; col_offset_cand[offset_num] = -1;
    offset_num++;
  }

  row_offset = bh / 8 - 1;
  col_offset = 1;
  if (is_sec_rec < 2)
    switchflag = 1;
  else
    switchflag = 0;
  while ((is_sec_rec == 0 && ((row_offset >=0) || col_offset < (bw / 8 + 1))) ||
         (is_sec_rec == 1 && col_offset < (bw / 8 + 1)) ||
         (is_sec_rec == 2 && row_offset >=0)) {
    switch (switchflag) {
      case 0:
        if (row_offset >= 0) {
          if (row_offset != bh / 16) {
            row_offset_cand[offset_num] = row_offset;
            col_offset_cand[offset_num] = -1;
            offset_num++;
          }
          row_offset--;
        }
        break;
      case 1:
        if (col_offset < (bw / 8 + 1)) {
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
    if ((col_offset < (bw / 8) ||
        (col_offset == (bw / 8) && is_right_available(bsize, mi_row, mi_col)))
        && check_inside(cm, mi_row + row_offset, mi_col + col_offset)) {
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
