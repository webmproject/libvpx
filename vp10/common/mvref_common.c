
/*
 *  Copyright (c) 2012 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include "vp10/common/mvref_common.h"

#if CONFIG_REF_MV

static uint8_t add_ref_mv_candidate(const MODE_INFO *const candidate_mi,
                                    const MB_MODE_INFO *const candidate,
                                    const MV_REFERENCE_FRAME rf[2],
                                    uint8_t *refmv_count,
                                    CANDIDATE_MV *ref_mv_stack,
                                    const int use_hp,
                                    int len, int block, int col) {
  const int weight = len;
  int index = 0, ref;
  int newmv_count = 0;

  assert(2 * weight < REF_CAT_LEVEL);

  if (rf[1] == NONE) {
    // single reference frame
    for (ref = 0; ref < 2; ++ref) {
      if (candidate->ref_frame[ref] == rf[0]) {
        int_mv this_refmv =
            get_sub_block_mv(candidate_mi, ref, col, block);
        lower_mv_precision(&this_refmv.as_mv, use_hp);

        for (index = 0; index < *refmv_count; ++index)
          if (ref_mv_stack[index].this_mv.as_int == this_refmv.as_int)
            break;

        if (index < *refmv_count)
          ref_mv_stack[index].weight += 2 * weight;

        // Add a new item to the list.
        if (index == *refmv_count) {
          ref_mv_stack[index].this_mv = this_refmv;
          ref_mv_stack[index].pred_mv =
              get_sub_block_pred_mv(candidate_mi, ref, col, block);
          ref_mv_stack[index].weight = 2 * weight;
          ++(*refmv_count);

#if CONFIG_EXT_INTER
          if (candidate->mode == NEWMV || candidate->mode == NEWFROMNEARMV)
#else
          if (candidate->mode == NEWMV)
#endif  // CONFIG_EXT_INTER
            ++newmv_count;
        }

        if (candidate_mi->mbmi.sb_type < BLOCK_8X8 && block >= 0) {
          int alt_block = 3 - block;
          this_refmv =
              get_sub_block_mv(candidate_mi, ref, col, alt_block);
          lower_mv_precision(&this_refmv.as_mv, use_hp);

          for (index = 0; index < *refmv_count; ++index)
            if (ref_mv_stack[index].this_mv.as_int == this_refmv.as_int)
              break;

          if (index < *refmv_count)
            ref_mv_stack[index].weight += weight;

          // Add a new item to the list.
          if (index == *refmv_count) {
            ref_mv_stack[index].this_mv = this_refmv;
            ref_mv_stack[index].pred_mv =
                get_sub_block_pred_mv(candidate_mi, ref, col, alt_block);
            ref_mv_stack[index].weight = weight;
            ++(*refmv_count);

#if CONFIG_EXT_INTER
          if (candidate->mode == NEWMV || candidate->mode == NEWFROMNEARMV)
#else
            if (candidate->mode == NEWMV)
#endif  // CONFIG_EXT_INTER
              ++newmv_count;
          }
        }
      }
    }
  } else {
    // compound reference frame
    if (candidate->ref_frame[0] == rf[0] &&
        candidate->ref_frame[1] == rf[1]) {
      int_mv this_refmv[2];

      for (ref = 0; ref < 2; ++ref) {
        this_refmv[ref] = get_sub_block_mv(candidate_mi, ref, col, block);
        lower_mv_precision(&this_refmv[ref].as_mv, use_hp);
      }

      for (index = 0; index < *refmv_count; ++index)
        if ((ref_mv_stack[index].this_mv.as_int == this_refmv[0].as_int) &&
            (ref_mv_stack[index].comp_mv.as_int == this_refmv[1].as_int))
          break;

      if (index < *refmv_count)
        ref_mv_stack[index].weight += 2 * weight;

      // Add a new item to the list.
      if (index == *refmv_count) {
        ref_mv_stack[index].this_mv = this_refmv[0];
        ref_mv_stack[index].comp_mv = this_refmv[1];
        ref_mv_stack[index].weight = 2 * weight;
        ++(*refmv_count);

#if CONFIG_EXT_INTER
        if (candidate->mode == NEW_NEWMV)
#else
        if (candidate->mode == NEWMV)
#endif  // CONFIG_EXT_INTER
          ++newmv_count;
      }

      if (candidate_mi->mbmi.sb_type < BLOCK_8X8 && block >= 0) {
        int alt_block = 3 - block;
        this_refmv[0] = get_sub_block_mv(candidate_mi, 0, col, alt_block);
        this_refmv[1] = get_sub_block_mv(candidate_mi, 1, col, alt_block);

        for (ref = 0; ref < 2; ++ref)
          lower_mv_precision(&this_refmv[ref].as_mv, use_hp);

        for (index = 0; index < *refmv_count; ++index)
          if (ref_mv_stack[index].this_mv.as_int == this_refmv[0].as_int &&
              ref_mv_stack[index].comp_mv.as_int == this_refmv[1].as_int)
            break;

        if (index < *refmv_count)
          ref_mv_stack[index].weight += weight;

        // Add a new item to the list.
        if (index == *refmv_count) {
          ref_mv_stack[index].this_mv = this_refmv[0];
          ref_mv_stack[index].comp_mv = this_refmv[1];
          ref_mv_stack[index].weight = weight;
          ++(*refmv_count);

#if CONFIG_EXT_INTER
          if (candidate->mode == NEW_NEWMV)
#else
          if (candidate->mode == NEWMV)
#endif  // CONFIG_EXT_INTER
            ++newmv_count;
        }
      }
    }
  }
  return newmv_count;
}

static uint8_t scan_row_mbmi(const VP10_COMMON *cm, const MACROBLOCKD *xd,
                             const int mi_row, const int mi_col, int block,
                             const MV_REFERENCE_FRAME rf[2],
                             int row_offset,
                             CANDIDATE_MV *ref_mv_stack,
                             uint8_t *refmv_count) {
  const TileInfo *const tile = &xd->tile;
  int i;
  uint8_t newmv_count = 0;

  for (i = 0; i < xd->n8_w && *refmv_count < MAX_REF_MV_STACK_SIZE;) {
    POSITION mi_pos;
    mi_pos.row = row_offset;
    mi_pos.col = i;

    if (is_inside(tile, mi_col, mi_row, &mi_pos)) {
      const MODE_INFO *const candidate_mi =
          xd->mi[mi_pos.row * xd->mi_stride + mi_pos.col];
      const MB_MODE_INFO *const candidate = &candidate_mi->mbmi;
      const int len = VPXMIN(xd->n8_w,
                             num_8x8_blocks_wide_lookup[candidate->sb_type]);

      newmv_count += add_ref_mv_candidate(candidate_mi, candidate, rf,
                                          refmv_count, ref_mv_stack,
                                          cm->allow_high_precision_mv,
                                          len, block, mi_pos.col);
      i += len;
    } else {
      ++i;
    }
  }

  return newmv_count;
}

static uint8_t scan_col_mbmi(const VP10_COMMON *cm, const MACROBLOCKD *xd,
                             const int mi_row, const int mi_col, int block,
                             const MV_REFERENCE_FRAME rf[2],
                             int col_offset,
                             CANDIDATE_MV *ref_mv_stack,
                             uint8_t *refmv_count) {
  const TileInfo *const tile = &xd->tile;
  int i;
  uint8_t newmv_count = 0;

  for (i = 0; i < xd->n8_h && *refmv_count < MAX_REF_MV_STACK_SIZE;) {
    POSITION mi_pos;
    mi_pos.row = i;
    mi_pos.col = col_offset;

    if (is_inside(tile, mi_col, mi_row, &mi_pos)) {
      const MODE_INFO *const candidate_mi =
          xd->mi[mi_pos.row * xd->mi_stride + mi_pos.col];
      const MB_MODE_INFO *const candidate = &candidate_mi->mbmi;
      const int len = VPXMIN(xd->n8_h,
                       num_8x8_blocks_high_lookup[candidate->sb_type]);

      newmv_count += add_ref_mv_candidate(candidate_mi, candidate, rf,
                                          refmv_count, ref_mv_stack,
                                          cm->allow_high_precision_mv,
                                          len, block, mi_pos.col);
      i += len;
    } else {
      ++i;
    }
  }

  return newmv_count;
}

static uint8_t scan_blk_mbmi(const VP10_COMMON *cm, const MACROBLOCKD *xd,
                             const int mi_row, const int mi_col, int block,
                             const MV_REFERENCE_FRAME rf[2],
                             int row_offset, int col_offset,
                             CANDIDATE_MV *ref_mv_stack,
                             uint8_t *refmv_count) {
  const TileInfo *const tile = &xd->tile;
  POSITION mi_pos;
  uint8_t newmv_count = 0;

  mi_pos.row = row_offset;
  mi_pos.col = col_offset;

  if (is_inside(tile, mi_col, mi_row, &mi_pos) &&
      *refmv_count < MAX_REF_MV_STACK_SIZE) {
    const MODE_INFO *const candidate_mi =
        xd->mi[mi_pos.row * xd->mi_stride + mi_pos.col];
    const MB_MODE_INFO *const candidate = &candidate_mi->mbmi;
    const int len = 1;

    newmv_count += add_ref_mv_candidate(candidate_mi, candidate, rf,
                                        refmv_count, ref_mv_stack,
                                        cm->allow_high_precision_mv,
                                        len, block, mi_pos.col);
  }  // Analyze a single 8x8 block motion information.
  return newmv_count;
}

static int has_top_right(const MACROBLOCKD *xd,
                         int mi_row, int mi_col, int bs) {
  // In a split partition all apart from the bottom right has a top right
  int has_tr = !((mi_row & bs) && (mi_col & bs));

  // bs > 0 and bs is a power of 2
  assert(bs > 0 && !(bs & (bs - 1)));

  // For each 4x4 group of blocks, when the bottom right is decoded the blocks
  // to the right have not been decoded therefore the bottom right does
  // not have a top right
  while (bs < MAX_MIB_SIZE) {
    if (mi_col & bs) {
      if ((mi_col & (2 * bs)) && (mi_row & (2 * bs))) {
        has_tr = 0;
        break;
      }
    } else {
      break;
    }
    bs <<= 1;
  }

  // The left hand of two vertical rectangles always has a top right (as the
  // block above will have been decoded)
  if (xd->n8_w < xd->n8_h)
    if (!xd->is_sec_rect)
      has_tr = 1;

  // The bottom of two horizontal rectangles never has a top right (as the block
  // to the right won't have been decoded)
  if (xd->n8_w > xd->n8_h)
    if (xd->is_sec_rect)
      has_tr = 0;

#if CONFIG_EXT_PARTITION_TYPES
  // The bottom left square of a Vertical A does not have a top right as it is
  // decoded before the right hand rectangle of the partition
  if (xd->mi[0]->mbmi.partition == PARTITION_VERT_A)
    if ((mi_row & bs) && !(mi_col & bs))
      has_tr = 0;
#endif  // CONFIG_EXT_PARTITION_TYPES

  return has_tr;
}

static void handle_sec_rect_block(const MB_MODE_INFO * const candidate,
                                  uint8_t refmv_count,
                                  CANDIDATE_MV *ref_mv_stack,
                                  MV_REFERENCE_FRAME ref_frame,
                                  int16_t *mode_context) {
  int rf, idx;

  for (rf = 0; rf < 2; ++rf) {
    if (candidate->ref_frame[rf] == ref_frame) {
      const int list_range = VPXMIN(refmv_count, MAX_MV_REF_CANDIDATES);

      const int_mv pred_mv = candidate->mv[rf];
      for (idx = 0; idx < list_range; ++idx)
        if (pred_mv.as_int == ref_mv_stack[idx].this_mv.as_int)
          break;

      if (idx < list_range) {
        if (idx == 0)
          mode_context[ref_frame] |= (1 << SKIP_NEARESTMV_OFFSET);
        else if (idx == 1)
          mode_context[ref_frame] |= (1 << SKIP_NEARMV_OFFSET);
      }
    }
  }
}

static void setup_ref_mv_list(const VP10_COMMON *cm, const MACROBLOCKD *xd,
                              MV_REFERENCE_FRAME ref_frame,
                              uint8_t *refmv_count,
                              CANDIDATE_MV *ref_mv_stack,
                              int_mv *mv_ref_list,
                              int block, int mi_row, int mi_col,
                              int16_t *mode_context) {
  int idx, nearest_refmv_count = 0;
  uint8_t newmv_count = 0;

  CANDIDATE_MV tmp_mv;
  int len, nr_len;

  const MV_REF *const prev_frame_mvs_base = cm->use_prev_frame_mvs ?
      cm->prev_frame->mvs + mi_row * cm->mi_cols + mi_col : NULL;

  int bs = VPXMAX(xd->n8_w, xd->n8_h);
  int has_tr = has_top_right(xd, mi_row, mi_col, bs);

  MV_REFERENCE_FRAME rf[2];
  vp10_set_ref_frame(rf, ref_frame);

  mode_context[ref_frame] = 0;
  *refmv_count = 0;

  // Scan the first above row mode info.
  newmv_count = scan_row_mbmi(cm, xd, mi_row, mi_col, block, rf,
                              -1, ref_mv_stack, refmv_count);
  // Scan the first left column mode info.
  newmv_count += scan_col_mbmi(cm, xd, mi_row, mi_col, block, rf,
                               -1, ref_mv_stack, refmv_count);

  // Check top-right boundary
  if (has_tr)
    newmv_count += scan_blk_mbmi(cm, xd, mi_row, mi_col, block, rf,
                                 -1, 1, ref_mv_stack, refmv_count);

  nearest_refmv_count = *refmv_count;

  for (idx = 0; idx < nearest_refmv_count; ++idx) {
    assert(ref_mv_stack[idx].weight > 0 &&
           ref_mv_stack[idx].weight < REF_CAT_LEVEL);
    ref_mv_stack[idx].weight += REF_CAT_LEVEL;
  }

  if (prev_frame_mvs_base && cm->show_frame && cm->last_show_frame
      && rf[1] == NONE) {
    int ref;
    int blk_row, blk_col;

    for (blk_row = 0; blk_row < xd->n8_h; ++blk_row) {
      for (blk_col = 0; blk_col < xd->n8_w; ++blk_col) {
        const MV_REF *prev_frame_mvs =
            prev_frame_mvs_base + blk_row * cm->mi_cols + blk_col;

        POSITION mi_pos;
        mi_pos.row = blk_row;
        mi_pos.col = blk_col;

        if (!is_inside(&xd->tile, mi_col, mi_row, &mi_pos))
          continue;

        for (ref = 0; ref < 2; ++ref) {
          if (prev_frame_mvs->ref_frame[ref] == ref_frame) {
            int_mv this_refmv = prev_frame_mvs->mv[ref];
            lower_mv_precision(&this_refmv.as_mv,
                               cm->allow_high_precision_mv);

            for (idx = 0; idx < *refmv_count; ++idx)
              if (this_refmv.as_int == ref_mv_stack[idx].this_mv.as_int)
                break;

            if (idx < *refmv_count)
              ref_mv_stack[idx].weight += 2;

            if (idx == *refmv_count &&
                *refmv_count < MAX_REF_MV_STACK_SIZE) {
              ref_mv_stack[idx].this_mv.as_int = this_refmv.as_int;
              ref_mv_stack[idx].weight = 2;
              ++(*refmv_count);

              if (abs(ref_mv_stack[idx].this_mv.as_mv.row) >= 8 ||
                  abs(ref_mv_stack[idx].this_mv.as_mv.col) >= 8)
                mode_context[ref_frame] |= (1 << ZEROMV_OFFSET);
            }
          }
        }
      }
    }
  }

  if (*refmv_count == nearest_refmv_count)
    mode_context[ref_frame] |= (1 << ZEROMV_OFFSET);

  // Analyze the top-left corner block mode info.
//  scan_blk_mbmi(cm, xd, mi_row, mi_col, block, ref_frame,
//                -1, -1, ref_mv_stack, refmv_count);

  // Scan the second outer area.
  scan_row_mbmi(cm, xd, mi_row, mi_col, block, rf,
                -2, ref_mv_stack, refmv_count);
  scan_col_mbmi(cm, xd, mi_row, mi_col, block, rf,
                -2, ref_mv_stack, refmv_count);

  // Scan the third outer area.
  scan_row_mbmi(cm, xd, mi_row, mi_col, block, rf,
                -3, ref_mv_stack, refmv_count);
  scan_col_mbmi(cm, xd, mi_row, mi_col, block, rf,
                -3, ref_mv_stack, refmv_count);

  // Scan the fourth outer area.
  scan_row_mbmi(cm, xd, mi_row, mi_col, block, rf,
                -4, ref_mv_stack, refmv_count);
  // Scan the third left row mode info.
  scan_col_mbmi(cm, xd, mi_row, mi_col, block, rf,
                -4, ref_mv_stack, refmv_count);

  switch (nearest_refmv_count) {
    case 0:
      mode_context[ref_frame] |= 0;
      if (*refmv_count >= 1)
        mode_context[ref_frame] |= 1;

      if (*refmv_count == 1)
        mode_context[ref_frame] |= (1 << REFMV_OFFSET);
      else if (*refmv_count >= 2)
        mode_context[ref_frame] |= (2 << REFMV_OFFSET);
      break;
    case 1:
      mode_context[ref_frame] |= (newmv_count > 0) ? 2 : 3;

      if (*refmv_count == 1)
        mode_context[ref_frame] |= (3 << REFMV_OFFSET);
      else if (*refmv_count >= 2)
        mode_context[ref_frame] |= (4 << REFMV_OFFSET);
      break;

    case 2:
    default:
      if (newmv_count >= 2)
        mode_context[ref_frame] |= 4;
      else if (newmv_count == 1)
        mode_context[ref_frame] |= 5;
      else
        mode_context[ref_frame] |= 6;

      mode_context[ref_frame] |= (5 << REFMV_OFFSET);
      break;
  }

  // Rank the likelihood and assign nearest and near mvs.
  len = nearest_refmv_count;
  while (len > 0) {
    nr_len = 0;
    for (idx = 1; idx < len; ++idx) {
      if (ref_mv_stack[idx - 1].weight < ref_mv_stack[idx].weight) {
        tmp_mv = ref_mv_stack[idx - 1];
        ref_mv_stack[idx - 1] = ref_mv_stack[idx];
        ref_mv_stack[idx] = tmp_mv;
        nr_len = idx;
      }
    }
    len = nr_len;
  }

  len = *refmv_count;
  while (len > nearest_refmv_count) {
    nr_len = nearest_refmv_count;
    for (idx = nearest_refmv_count + 1; idx < len; ++idx) {
      if (ref_mv_stack[idx - 1].weight < ref_mv_stack[idx].weight) {
        tmp_mv = ref_mv_stack[idx - 1];
        ref_mv_stack[idx - 1] = ref_mv_stack[idx];
        ref_mv_stack[idx] = tmp_mv;
        nr_len = idx;
      }
    }
    len = nr_len;
  }

  // TODO(jingning): Clean-up needed.
  if (xd->is_sec_rect) {
    if (xd->n8_w < xd->n8_h) {
      const MODE_INFO *const candidate_mi = xd->mi[-1];
      const MB_MODE_INFO *const candidate = &candidate_mi->mbmi;
      handle_sec_rect_block(candidate, nearest_refmv_count, ref_mv_stack,
                            ref_frame, mode_context);
    }

    if (xd->n8_w > xd->n8_h) {
      const MODE_INFO *const candidate_mi = xd->mi[-xd->mi_stride];
      const MB_MODE_INFO *const candidate = &candidate_mi->mbmi;
      handle_sec_rect_block(candidate, nearest_refmv_count, ref_mv_stack,
                            ref_frame, mode_context);
    }
  }

  if (rf[1] > NONE) {
    for (idx = 0; idx < *refmv_count; ++idx) {
      clamp_mv_ref(&ref_mv_stack[idx].this_mv.as_mv,
                   xd->n8_w << 3 , xd->n8_h << 3, xd);
      clamp_mv_ref(&ref_mv_stack[idx].comp_mv.as_mv,
                   xd->n8_w << 3 , xd->n8_h << 3, xd);
    }
  } else {
    for (idx = 0; idx < VPXMIN(MAX_MV_REF_CANDIDATES, *refmv_count); ++idx) {
      mv_ref_list[idx].as_int = ref_mv_stack[idx].this_mv.as_int;
      clamp_mv_ref(&mv_ref_list[idx].as_mv,
                   xd->n8_w << 3, xd->n8_h << 3, xd);
    }
  }
}
#endif

// This function searches the neighbourhood of a given MB/SB
// to try and find candidate reference vectors.
static void find_mv_refs_idx(const VP10_COMMON *cm, const MACROBLOCKD *xd,
                             MODE_INFO *mi, MV_REFERENCE_FRAME ref_frame,
                             int_mv *mv_ref_list,
                             int block, int mi_row, int mi_col,
                             find_mv_refs_sync sync, void *const data,
                             int16_t *mode_context) {
  const int *ref_sign_bias = cm->ref_frame_sign_bias;
  int i, refmv_count = 0;
  const POSITION *const mv_ref_search = mv_ref_blocks[mi->mbmi.sb_type];
  int different_ref_found = 0;
  int context_counter = 0;
  const MV_REF *const  prev_frame_mvs = cm->use_prev_frame_mvs ?
      cm->prev_frame->mvs + mi_row * cm->mi_cols + mi_col : NULL;
  const TileInfo *const tile = &xd->tile;
  const int bw = num_8x8_blocks_wide_lookup[mi->mbmi.sb_type] << 3;
  const int bh = num_8x8_blocks_high_lookup[mi->mbmi.sb_type] << 3;

  // The nearest 2 blocks are treated differently
  // if the size < 8x8 we get the mv from the bmi substructure,
  // and we also need to keep a mode count.
  for (i = 0; i < 2; ++i) {
    const POSITION *const mv_ref = &mv_ref_search[i];
    if (is_inside(tile, mi_col, mi_row, mv_ref)) {
      const MODE_INFO *const candidate_mi = xd->mi[mv_ref->col + mv_ref->row *
                                                   xd->mi_stride];
      const MB_MODE_INFO *const candidate = &candidate_mi->mbmi;
      // Keep counts for entropy encoding.
      context_counter += mode_2_counter[candidate->mode];
      different_ref_found = 1;

      if (candidate->ref_frame[0] == ref_frame)
        ADD_MV_REF_LIST(get_sub_block_mv(candidate_mi, 0, mv_ref->col, block),
                        refmv_count, mv_ref_list, bw, bh, xd, Done);
      else if (candidate->ref_frame[1] == ref_frame)
        ADD_MV_REF_LIST(get_sub_block_mv(candidate_mi, 1, mv_ref->col, block),
                        refmv_count, mv_ref_list, bw, bh, xd, Done);
    }
  }

  // Check the rest of the neighbors in much the same way
  // as before except we don't need to keep track of sub blocks or
  // mode counts.
  for (; i < MVREF_NEIGHBOURS; ++i) {
    const POSITION *const mv_ref = &mv_ref_search[i];
    if (is_inside(tile, mi_col, mi_row, mv_ref)) {
      const MB_MODE_INFO *const candidate = &xd->mi[mv_ref->col + mv_ref->row *
                                                    xd->mi_stride]->mbmi;
      different_ref_found = 1;

      if (candidate->ref_frame[0] == ref_frame)
        ADD_MV_REF_LIST(candidate->mv[0], refmv_count, mv_ref_list,
                        bw, bh, xd, Done);
      else if (candidate->ref_frame[1] == ref_frame)
        ADD_MV_REF_LIST(candidate->mv[1], refmv_count, mv_ref_list,
                        bw, bh, xd, Done);
    }
  }

  // TODO(hkuang): Remove this sync after fixing pthread_cond_broadcast
  // on windows platform. The sync here is unncessary if use_perv_frame_mvs
  // is 0. But after removing it, there will be hang in the unit test on windows
  // due to several threads waiting for a thread's signal.
#if defined(_WIN32) && !HAVE_PTHREAD_H
    if (cm->frame_parallel_decode && sync != NULL) {
      sync(data, mi_row);
    }
#endif

  // Check the last frame's mode and mv info.
  if (cm->use_prev_frame_mvs) {
    // Synchronize here for frame parallel decode if sync function is provided.
    if (cm->frame_parallel_decode && sync != NULL) {
      sync(data, mi_row);
    }

    if (prev_frame_mvs->ref_frame[0] == ref_frame) {
      ADD_MV_REF_LIST(prev_frame_mvs->mv[0], refmv_count, mv_ref_list,
                      bw, bh, xd, Done);
    } else if (prev_frame_mvs->ref_frame[1] == ref_frame) {
      ADD_MV_REF_LIST(prev_frame_mvs->mv[1], refmv_count, mv_ref_list,
                      bw, bh, xd, Done);
    }
  }

  // Since we couldn't find 2 mvs from the same reference frame
  // go back through the neighbors and find motion vectors from
  // different reference frames.
  if (different_ref_found) {
    for (i = 0; i < MVREF_NEIGHBOURS; ++i) {
      const POSITION *mv_ref = &mv_ref_search[i];
      if (is_inside(tile, mi_col, mi_row, mv_ref)) {
        const MB_MODE_INFO *const candidate = &xd->mi[mv_ref->col + mv_ref->row
                                              * xd->mi_stride]->mbmi;

        // If the candidate is INTRA we don't want to consider its mv.
        IF_DIFF_REF_FRAME_ADD_MV(candidate, ref_frame, ref_sign_bias,
                                 refmv_count, mv_ref_list, bw, bh, xd, Done);
      }
    }
  }

  // Since we still don't have a candidate we'll try the last frame.
  if (cm->use_prev_frame_mvs) {
    if (prev_frame_mvs->ref_frame[0] != ref_frame &&
        prev_frame_mvs->ref_frame[0] > INTRA_FRAME) {
      int_mv mv = prev_frame_mvs->mv[0];
      if (ref_sign_bias[prev_frame_mvs->ref_frame[0]] !=
          ref_sign_bias[ref_frame]) {
        mv.as_mv.row *= -1;
        mv.as_mv.col *= -1;
      }
      ADD_MV_REF_LIST(mv, refmv_count, mv_ref_list, bw, bh, xd, Done);
    }

    if (prev_frame_mvs->ref_frame[1] > INTRA_FRAME &&
        prev_frame_mvs->ref_frame[1] != ref_frame) {
      int_mv mv = prev_frame_mvs->mv[1];
      if (ref_sign_bias[prev_frame_mvs->ref_frame[1]] !=
          ref_sign_bias[ref_frame]) {
        mv.as_mv.row *= -1;
        mv.as_mv.col *= -1;
      }
      ADD_MV_REF_LIST(mv, refmv_count, mv_ref_list, bw, bh, xd, Done);
    }
  }

Done:
  if (mode_context)
    mode_context[ref_frame] = counter_to_context[context_counter];
  for (i = refmv_count; i < MAX_MV_REF_CANDIDATES; ++i)
      mv_ref_list[i].as_int = 0;
}

#if CONFIG_EXT_INTER
// This function keeps a mode count for a given MB/SB
void vp10_update_mv_context(const MACROBLOCKD *xd,
                            MODE_INFO *mi, MV_REFERENCE_FRAME ref_frame,
                            int_mv *mv_ref_list,
                            int block, int mi_row, int mi_col,
                            int16_t *mode_context) {
  int i, refmv_count = 0;
  const POSITION *const mv_ref_search = mv_ref_blocks[mi->mbmi.sb_type];
  int context_counter = 0;
  const int bw = num_8x8_blocks_wide_lookup[mi->mbmi.sb_type] << 3;
  const int bh = num_8x8_blocks_high_lookup[mi->mbmi.sb_type] << 3;
  const TileInfo *const tile = &xd->tile;

  // Blank the reference vector list
  memset(mv_ref_list, 0, sizeof(*mv_ref_list) * MAX_MV_REF_CANDIDATES);

  // The nearest 2 blocks are examined only.
  // If the size < 8x8, we get the mv from the bmi substructure;
  for (i = 0; i < 2; ++i) {
    const POSITION *const mv_ref = &mv_ref_search[i];
    if (is_inside(tile, mi_col, mi_row, mv_ref)) {
      const MODE_INFO *const candidate_mi =
          xd->mi[mv_ref->col + mv_ref->row * xd->mi_stride];
      const MB_MODE_INFO *const candidate = &candidate_mi->mbmi;

      // Keep counts for entropy encoding.
      context_counter += mode_2_counter[candidate->mode];

      if (candidate->ref_frame[0] == ref_frame) {
        ADD_MV_REF_LIST(get_sub_block_mv(candidate_mi, 0, mv_ref->col, block),
                        refmv_count, mv_ref_list, bw, bh, xd, Done);
      } else if (candidate->ref_frame[1] == ref_frame) {
        ADD_MV_REF_LIST(get_sub_block_mv(candidate_mi, 1, mv_ref->col, block),
                        refmv_count, mv_ref_list, bw, bh, xd, Done);
      }
    }
  }

 Done:

  if (mode_context)
    mode_context[ref_frame] = counter_to_context[context_counter];
}
#endif  // CONFIG_EXT_INTER

void vp10_find_mv_refs(const VP10_COMMON *cm, const MACROBLOCKD *xd,
                      MODE_INFO *mi, MV_REFERENCE_FRAME ref_frame,
#if CONFIG_REF_MV
                      uint8_t *ref_mv_count,
                      CANDIDATE_MV *ref_mv_stack,
#if CONFIG_EXT_INTER
                      int16_t *compound_mode_context,
#endif  // CONFIG_EXT_INTER
#endif
                      int_mv *mv_ref_list,
                      int mi_row, int mi_col,
                      find_mv_refs_sync sync, void *const data,
                      int16_t *mode_context) {
#if CONFIG_REF_MV
  int idx, all_zero = 1;
#endif
#if CONFIG_EXT_INTER
  vp10_update_mv_context(xd, mi, ref_frame, mv_ref_list, -1,
                         mi_row, mi_col,
#if CONFIG_REF_MV
                         compound_mode_context);
#else
                         mode_context);
#endif  // CONFIG_REF_MV
  find_mv_refs_idx(cm, xd, mi, ref_frame, mv_ref_list, -1,
                   mi_row, mi_col, sync, data, NULL);
#else
  find_mv_refs_idx(cm, xd, mi, ref_frame, mv_ref_list, -1,
                   mi_row, mi_col, sync, data, mode_context);
#endif  // CONFIG_EXT_INTER

#if CONFIG_REF_MV
  setup_ref_mv_list(cm, xd, ref_frame, ref_mv_count, ref_mv_stack,
                    mv_ref_list, -1, mi_row, mi_col, mode_context);

  for (idx = 0; idx < MAX_MV_REF_CANDIDATES; ++idx)
    if (mv_ref_list[idx].as_int != 0)
      all_zero = 0;

  if (all_zero)
    mode_context[ref_frame] |= (1 << ALL_ZERO_FLAG_OFFSET);
#endif
}

void vp10_find_best_ref_mvs(int allow_hp,
                           int_mv *mvlist, int_mv *nearest_mv,
                           int_mv *near_mv) {
  int i;
  // Make sure all the candidates are properly clamped etc
  for (i = 0; i < MAX_MV_REF_CANDIDATES; ++i) {
    lower_mv_precision(&mvlist[i].as_mv, allow_hp);
  }
  *nearest_mv = mvlist[0];
  *near_mv = mvlist[1];
}

void vp10_append_sub8x8_mvs_for_idx(VP10_COMMON *cm, MACROBLOCKD *xd,
                                    int block, int ref, int mi_row, int mi_col,
#if CONFIG_REF_MV
                                    CANDIDATE_MV *ref_mv_stack,
                                    uint8_t *ref_mv_count,
#endif
#if CONFIG_EXT_INTER
                                    int_mv *mv_list,
#endif  // CONFIG_EXT_INTER
                                    int_mv *nearest_mv, int_mv *near_mv) {
#if !CONFIG_EXT_INTER
  int_mv mv_list[MAX_MV_REF_CANDIDATES];
#endif  // !CONFIG_EXT_INTER
  MODE_INFO *const mi = xd->mi[0];
  b_mode_info *bmi = mi->bmi;
  int n;
#if CONFIG_REF_MV
  CANDIDATE_MV tmp_mv;
  uint8_t idx;
  uint8_t above_count = 0, left_count = 0;
  MV_REFERENCE_FRAME rf[2] = { mi->mbmi.ref_frame[ref], NONE };
  *ref_mv_count = 0;
#endif

  assert(MAX_MV_REF_CANDIDATES == 2);

  find_mv_refs_idx(cm, xd, mi, mi->mbmi.ref_frame[ref], mv_list, block,
                   mi_row, mi_col, NULL, NULL, NULL);

#if CONFIG_REF_MV
  scan_blk_mbmi(cm, xd, mi_row, mi_col, block, rf,
                -1, 0, ref_mv_stack, ref_mv_count);
  above_count = *ref_mv_count;

  scan_blk_mbmi(cm, xd, mi_row, mi_col, block, rf,
                0, -1, ref_mv_stack, ref_mv_count);
  left_count = *ref_mv_count - above_count;

  if (above_count > 1 && left_count > 0) {
    tmp_mv = ref_mv_stack[1];
    ref_mv_stack[1] = ref_mv_stack[above_count];
    ref_mv_stack[above_count] = tmp_mv;
  }

  for (idx = 0; idx < *ref_mv_count; ++idx)
    clamp_mv_ref(&ref_mv_stack[idx].this_mv.as_mv,
                 xd->n8_w << 3, xd->n8_h << 3, xd);

  for (idx = 0; idx < VPXMIN(MAX_MV_REF_CANDIDATES, *ref_mv_count); ++idx)
    mv_list[idx].as_int = ref_mv_stack[idx].this_mv.as_int;
#endif

  near_mv->as_int = 0;
  switch (block) {
    case 0:
      nearest_mv->as_int = mv_list[0].as_int;
      near_mv->as_int = mv_list[1].as_int;
      break;
    case 1:
    case 2:
      nearest_mv->as_int = bmi[0].as_mv[ref].as_int;
      for (n = 0; n < MAX_MV_REF_CANDIDATES; ++n)
        if (nearest_mv->as_int != mv_list[n].as_int) {
          near_mv->as_int = mv_list[n].as_int;
          break;
        }
      break;
    case 3: {
      int_mv candidates[2 + MAX_MV_REF_CANDIDATES];
      candidates[0] = bmi[1].as_mv[ref];
      candidates[1] = bmi[0].as_mv[ref];
      candidates[2] = mv_list[0];
      candidates[3] = mv_list[1];

      nearest_mv->as_int = bmi[2].as_mv[ref].as_int;
      for (n = 0; n < 2 + MAX_MV_REF_CANDIDATES; ++n)
        if (nearest_mv->as_int != candidates[n].as_int) {
          near_mv->as_int = candidates[n].as_int;
          break;
        }
      break;
    }
    default:
      assert(0 && "Invalid block index.");
  }
}
