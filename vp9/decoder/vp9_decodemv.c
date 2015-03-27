/*
  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include <assert.h>

#include "vp9/common/vp9_common.h"
#include "vp9/common/vp9_entropy.h"
#include "vp9/common/vp9_entropymode.h"
#include "vp9/common/vp9_entropymv.h"
#include "vp9/common/vp9_mvref_common.h"
#if CONFIG_PALETTE
#include "vp9/common/vp9_palette.h"
#endif
#include "vp9/common/vp9_pred_common.h"
#include "vp9/common/vp9_reconinter.h"
#include "vp9/common/vp9_seg_common.h"

#include "vp9/decoder/vp9_decodemv.h"
#include "vp9/decoder/vp9_decodeframe.h"
#include "vp9/decoder/vp9_reader.h"

static PREDICTION_MODE read_intra_mode(vp9_reader *r, const vp9_prob *p) {
  return (PREDICTION_MODE)vp9_read_tree(r, vp9_intra_mode_tree, p);
}

static PREDICTION_MODE read_intra_mode_y(VP9_COMMON *cm, vp9_reader *r,
                                            int size_group) {
  const PREDICTION_MODE y_mode =
      read_intra_mode(r, cm->fc.y_mode_prob[size_group]);
  if (!cm->frame_parallel_decoding_mode)
    ++cm->counts.y_mode[size_group][y_mode];
  return y_mode;
}

static PREDICTION_MODE read_intra_mode_uv(VP9_COMMON *cm, vp9_reader *r,
                                          PREDICTION_MODE y_mode) {
  const PREDICTION_MODE uv_mode = read_intra_mode(r,
                                         cm->fc.uv_mode_prob[y_mode]);
  if (!cm->frame_parallel_decoding_mode)
    ++cm->counts.uv_mode[y_mode][uv_mode];
  return uv_mode;
}

#if CONFIG_COMPOUND_MODES
static PREDICTION_MODE read_inter_compound_mode(VP9_COMMON *cm, vp9_reader *r,
                                                int ctx) {
  int mode = 0;
  mode = vp9_read_tree(r, vp9_inter_compound_mode_tree,
                       cm->fc.inter_compound_mode_probs[ctx]);
  if (!cm->frame_parallel_decoding_mode) {
    ++cm->counts.inter_compound_mode[ctx][mode];
  }
  assert(is_inter_compound_mode(NEAREST_NEARESTMV + mode));
  return NEAREST_NEARESTMV + mode;
}
#endif

static PREDICTION_MODE read_inter_mode(VP9_COMMON *cm, vp9_reader *r,
                                       int ctx) {
  const int mode = vp9_read_tree(r, vp9_inter_mode_tree,
                                 cm->fc.inter_mode_probs[ctx]);
  if (!cm->frame_parallel_decoding_mode)
    ++cm->counts.inter_mode[ctx][mode];
  return NEARESTMV + mode;
}

#if CONFIG_COPY_MODE
static COPY_MODE read_copy_mode(VP9_COMMON *cm, vp9_reader *r,
                                int num_candidate, int ctx) {
  COPY_MODE mode;

  switch (num_candidate) {
    case 0:
      assert(0);
      break;
    case 1:
      mode = REF0;
      break;
    case 2:
      mode = REF0 + vp9_read_tree(r, vp9_copy_mode_tree_l2,
                                  cm->fc.copy_mode_probs_l2[ctx]);
      if (!cm->frame_parallel_decoding_mode)
          ++cm->counts.copy_mode_l2[ctx][mode - REF0];
      break;
    default:
      mode = REF0 + vp9_read_tree(r, vp9_copy_mode_tree,
                                  cm->fc.copy_mode_probs[ctx]);
      if (!cm->frame_parallel_decoding_mode)
          ++cm->counts.copy_mode[ctx][mode - REF0];
      break;
  }

  return mode;
}
#endif  // CONFIG_COPY_MODE

static int read_segment_id(vp9_reader *r, const struct segmentation *seg) {
  return vp9_read_tree(r, vp9_segment_tree, seg->tree_probs);
}

static TX_SIZE read_selected_tx_size(VP9_COMMON *cm, MACROBLOCKD *xd,
                                     TX_SIZE max_tx_size, vp9_reader *r) {
  const int ctx = vp9_get_tx_size_context(xd);
  const vp9_prob *tx_probs = get_tx_probs(max_tx_size, ctx, &cm->fc.tx_probs);
  int tx_size = vp9_read(r, tx_probs[0]);
  if (tx_size != TX_4X4 && max_tx_size >= TX_16X16) {
    tx_size += vp9_read(r, tx_probs[1]);
    if (tx_size != TX_8X8 && max_tx_size >= TX_32X32) {
      tx_size += vp9_read(r, tx_probs[2]);
#if CONFIG_TX64X64
      if (tx_size != TX_16X16 && max_tx_size >= TX_64X64) {
        tx_size += vp9_read(r, tx_probs[3]);
      }
#endif
    }
  }

  if (!cm->frame_parallel_decoding_mode)
    ++get_tx_counts(max_tx_size, ctx, &cm->counts.tx)[tx_size];
  return (TX_SIZE)tx_size;
}

static TX_SIZE read_tx_size(VP9_COMMON *cm, MACROBLOCKD *xd, TX_MODE tx_mode,
                            BLOCK_SIZE bsize, int allow_select, vp9_reader *r) {
  const TX_SIZE max_tx_size = max_txsize_lookup[bsize];
  if (allow_select && tx_mode == TX_MODE_SELECT && bsize >= BLOCK_8X8)
    return read_selected_tx_size(cm, xd, max_tx_size, r);
  else
    return MIN(max_tx_size, tx_mode_to_biggest_tx_size[tx_mode]);
}

static void set_segment_id(VP9_COMMON *cm, BLOCK_SIZE bsize,
                           int mi_row, int mi_col, int segment_id) {
  const int mi_offset = mi_row * cm->mi_cols + mi_col;
  const int bw = num_8x8_blocks_wide_lookup[bsize];
  const int bh = num_8x8_blocks_high_lookup[bsize];
  const int xmis = MIN(cm->mi_cols - mi_col, bw);
  const int ymis = MIN(cm->mi_rows - mi_row, bh);
  int x, y;

  assert(segment_id >= 0 && segment_id < MAX_SEGMENTS);

  for (y = 0; y < ymis; y++)
    for (x = 0; x < xmis; x++)
      cm->last_frame_seg_map[mi_offset + y * cm->mi_cols + x] = segment_id;
}

static int read_intra_segment_id(VP9_COMMON *const cm, MACROBLOCKD *const xd,
                                 int mi_row, int mi_col,
                                 vp9_reader *r) {
  struct segmentation *const seg = &cm->seg;
  const BLOCK_SIZE bsize = xd->mi[0].src_mi->mbmi.sb_type;
  int segment_id;

  if (!seg->enabled)
    return 0;  // Default for disabled segmentation

  if (!seg->update_map)
    return 0;

  segment_id = read_segment_id(r, seg);
  set_segment_id(cm, bsize, mi_row, mi_col, segment_id);
  return segment_id;
}

static int read_inter_segment_id(VP9_COMMON *const cm, MACROBLOCKD *const xd,
                                 int mi_row, int mi_col, vp9_reader *r) {
  struct segmentation *const seg = &cm->seg;
  MB_MODE_INFO *const mbmi = &xd->mi[0].src_mi->mbmi;
  const BLOCK_SIZE bsize = mbmi->sb_type;
  int predicted_segment_id, segment_id;

  if (!seg->enabled)
    return 0;  // Default for disabled segmentation

  predicted_segment_id = vp9_get_segment_id(cm, cm->last_frame_seg_map,
                                            bsize, mi_row, mi_col);
  if (!seg->update_map)
    return predicted_segment_id;

  if (seg->temporal_update) {
    const vp9_prob pred_prob = vp9_get_pred_prob_seg_id(seg, xd);
    mbmi->seg_id_predicted = vp9_read(r, pred_prob);
    segment_id = mbmi->seg_id_predicted ? predicted_segment_id
                                        : read_segment_id(r, seg);
  } else {
    segment_id = read_segment_id(r, seg);
  }
  set_segment_id(cm, bsize, mi_row, mi_col, segment_id);
  return segment_id;
}

static int read_skip(VP9_COMMON *cm, const MACROBLOCKD *xd,
                     int segment_id, vp9_reader *r) {
  if (vp9_segfeature_active(&cm->seg, segment_id, SEG_LVL_SKIP)) {
    return 1;
  } else {
    const int ctx = vp9_get_skip_context(xd);
    const int skip = vp9_read(r, cm->fc.skip_probs[ctx]);
    if (!cm->frame_parallel_decoding_mode)
      ++cm->counts.skip[ctx][skip];
    return skip;
  }
}

static void read_intra_frame_mode_info(VP9_COMMON *const cm,
                                       MACROBLOCKD *const xd,
                                       int mi_row, int mi_col, vp9_reader *r) {
  MODE_INFO *const mi = xd->mi[0].src_mi;
  MB_MODE_INFO *const mbmi = &mi->mbmi;
  const MODE_INFO *above_mi = xd->mi[-cm->mi_stride].src_mi;
  const MODE_INFO *left_mi  = xd->left_available ? xd->mi[-1].src_mi : NULL;
  const BLOCK_SIZE bsize = mbmi->sb_type;
  int i;

  mbmi->segment_id = read_intra_segment_id(cm, xd, mi_row, mi_col, r);
  mbmi->skip = read_skip(cm, xd, mbmi->segment_id, r);
#if CONFIG_PALETTE
  if (bsize >= BLOCK_8X8 && cm->allow_palette_mode) {
    int palette_ctx = 0;
    if (above_mi)
      palette_ctx += (above_mi->mbmi.palette_enabled[0] == 1);
    if (left_mi)
      palette_ctx += (left_mi->mbmi.palette_enabled[0] == 1);
    mbmi->palette_enabled[0] =
        vp9_read(r,
                 cm->fc.palette_enabled_prob[bsize - BLOCK_8X8][palette_ctx]);
    mbmi->palette_enabled[1] =
        vp9_read(r, cm->fc.palette_uv_enabled_prob[mbmi->palette_enabled[0]]);
  } else {
    mbmi->palette_enabled[0] = 0;
    mbmi->palette_enabled[1] = 0;
  }

  if (mbmi->palette_enabled[0]) {
    int i, j, m1, m2, val, n, color_idx, color_ctx;
    int rows = 4 * num_4x4_blocks_high_lookup[bsize];
    int cols = 4 * num_4x4_blocks_wide_lookup[bsize];
    int color_order[PALETTE_MAX_SIZE];
    uint8_t *color_map = xd->plane[0].color_index_map;

    mbmi->mode = DC_PRED;
    mbmi->palette_size[0] =
        vp9_read_tree(r, vp9_palette_size_tree,
                      cm->fc.palette_size_prob[bsize - BLOCK_8X8]);
    mbmi->palette_size[0] += 2;
    if ((xd->plane[1].subsampling_x && xd->plane[1].subsampling_y)
        || !mbmi->palette_enabled[1])
      mbmi->palette_indexed_size =
          vp9_decode_uniform(r, MIN(mbmi->palette_size[0] + 1, 8));
    else
      mbmi->palette_indexed_size = 0;
    mbmi->palette_literal_size = mbmi->palette_size[0] -
        mbmi->palette_indexed_size;

    if (PALETTE_DELTA_BIT)
      mbmi->palette_delta_bitdepth =
          vp9_read_literal(r, PALETTE_DELTA_BIT);
    else
      mbmi->palette_delta_bitdepth = 0;

    m1 = mbmi->palette_indexed_size;
    m2 = mbmi->palette_literal_size;
    n = mbmi->palette_size[0];

    if (m1 > 0) {
      for (i = 0; i < m1; i++)
        mbmi->palette_indexed_colors[i] =
            vp9_read_literal(r, vp9_get_bit_depth(cm->current_palette_size));
      if (mbmi->palette_delta_bitdepth > 0) {
        int s;
        for (i = 0; i < m1; i++) {
          s = vp9_read_bit(r);
          s = 1 - 2 * s;
          mbmi->palette_color_delta[i] =
              s * vp9_read_literal(r, mbmi->palette_delta_bitdepth);
        }
      } else {
        vpx_memset(mbmi->palette_color_delta, 0,
                   m1 * sizeof(mbmi->palette_color_delta[0]));
      }
      for (i = 0; i < m1; i++) {
        val = cm->current_palette_colors[mbmi->palette_indexed_colors[i]];
        if (mbmi->palette_color_delta[i])
          val += mbmi->palette_color_delta[i];
        mbmi->palette_colors[i] = val;
      }
    }
    if (m2 > 0) {
      for (i = 0; i < m2; i++) {
        mbmi->palette_literal_colors[i] = vp9_read_literal(r, 8);
        mbmi->palette_colors[m1 + i] = mbmi->palette_literal_colors[i];
      }
    }

    vp9_palette_color_insertion(cm->current_palette_colors,
                                &cm ->current_palette_size,
                                cm->current_palette_count, mbmi);

    color_map[0] = vp9_read_literal(r,
                                    vp9_get_bit_depth(mbmi->palette_size[0]));
    for (i = 0; i < rows; i++) {
      for (j = (i == 0 ? 1 : 0); j < cols; j++) {
        color_ctx = vp9_get_palette_color_context(color_map, cols, i, j, n,
                                                  color_order);
        color_idx = vp9_read_tree(r, vp9_palette_color_tree,
                                  cm->fc.palette_color_prob[n - 2][color_ctx]);

        color_map[i * cols + j] = color_order[color_idx];
      }
    }

    mbmi->tx_size = MIN(max_txsize_lookup[bsize],
                        tx_mode_to_biggest_tx_size[cm->tx_mode]);
  }

  if (mbmi->palette_enabled[1]) {
    int i, j;
    int rows = 4 * num_4x4_blocks_high_lookup[bsize] >>
        xd->plane[1].subsampling_y;
    int cols = 4 * num_4x4_blocks_wide_lookup[bsize] >>
        xd->plane[1].subsampling_x;

    mbmi->uv_mode = DC_PRED;
    if (xd->plane[1].subsampling_x && xd->plane[1].subsampling_y) {
      mbmi->palette_size[1] =
          vp9_read_tree(r, vp9_palette_size_tree,
                        cm->fc.palette_uv_size_prob[bsize - BLOCK_8X8]);
      mbmi->palette_size[1] += 2;
    } else {
      mbmi->palette_size[1] = mbmi->palette_size[0];
    }

    for (i = 0; i < mbmi->palette_size[1]; i++)
      mbmi->palette_colors[PALETTE_MAX_SIZE + i] = vp9_read_literal(r, 8);
    for (i = 0; i < mbmi->palette_size[1]; i++)
      mbmi->palette_colors[2 * PALETTE_MAX_SIZE + i] = vp9_read_literal(r, 8);

    if (xd->plane[1].subsampling_x && xd->plane[1].subsampling_y) {
      int color_idx = 0, color_ctx = 0;
      int n = mbmi->palette_size[1];
      int color_order[PALETTE_MAX_SIZE];
      uint8_t *color_map = xd->plane[1].color_index_map;

      color_map[0] = vp9_read_literal(r, vp9_get_bit_depth(n));
      for (i = 0; i < rows; i++) {
        for (j = (i == 0 ? 1 : 0); j < cols; j++) {
          color_ctx = vp9_get_palette_color_context(color_map, cols, i, j, n,
                                                    color_order);
          color_idx = vp9_read_tree(r, vp9_palette_color_tree,
                              cm->fc.palette_uv_color_prob[n - 2][color_ctx]);
          color_map[i * cols + j] = color_order[color_idx];
        }
      }
    }
  }

  if (!mbmi->palette_enabled[0]) {
    mbmi->tx_size = read_tx_size(cm, xd, cm->tx_mode, bsize, 1, r);
  }
#else
  mbmi->tx_size = read_tx_size(cm, xd, cm->tx_mode, bsize, 1, r);
#endif
  mbmi->ref_frame[0] = INTRA_FRAME;
  mbmi->ref_frame[1] = NONE;

#if CONFIG_TX_SKIP
  if (mbmi->sb_type >= BLOCK_8X8) {
    int q_idx = vp9_get_qindex(&cm->seg, mbmi->segment_id, cm->base_qindex);
    int try_tx_skip = q_idx <= TX_SKIP_Q_THRESH_INTRA;
    if (try_tx_skip) {
      if (xd->lossless) {
        if (mbmi->tx_size == TX_4X4)
          mbmi->tx_skip[0] = vp9_read(r, cm->fc.y_tx_skip_prob[0]);
        else
          mbmi->tx_skip[0] = 1;

        if (get_uv_tx_size(mbmi, &xd->plane[1]) == TX_4X4)
          mbmi->tx_skip[1] =
            vp9_read(r, cm->fc.uv_tx_skip_prob[mbmi->tx_skip[0]]);
        else
          mbmi->tx_skip[1] = 1;
        } else {
          mbmi->tx_skip[0] = vp9_read(r, cm->fc.y_tx_skip_prob[0]);
          mbmi->tx_skip[1] =
            vp9_read(r, cm->fc.uv_tx_skip_prob[mbmi->tx_skip[0]]);
        }
    } else {
      mbmi->tx_skip[0] = 0;
      mbmi->tx_skip[1] = 0;
    }
  } else {
    mbmi->tx_skip[0] = 0;
    mbmi->tx_skip[1] = 0;
  }
#endif

  switch (bsize) {
    case BLOCK_4X4:
#if CONFIG_FILTERINTRA
      for (i = 0; i < 4; ++i) {
#else
      for (i = 0; i < 4; ++i)
#endif
        mi->bmi[i].as_mode =
            read_intra_mode(r, get_y_mode_probs(mi, above_mi, left_mi, i));
#if CONFIG_FILTERINTRA
        if (is_filter_allowed(mi->bmi[i].as_mode))
          mi->b_filter_info[i] =
              vp9_read(r, cm->fc.filterintra_prob[0][mi->bmi[i].as_mode]);
        else
          mi->b_filter_info[i] = 0;
      }
      mbmi->filterbit = mi->b_filter_info[3];
#endif
      mbmi->mode = mi->bmi[3].as_mode;
      break;
    case BLOCK_4X8:
      mi->bmi[0].as_mode = mi->bmi[2].as_mode =
          read_intra_mode(r, get_y_mode_probs(mi, above_mi, left_mi, 0));
#if CONFIG_FILTERINTRA
      if (is_filter_allowed(mi->bmi[0].as_mode))
        mi->b_filter_info[0] = mi->b_filter_info[2] =
            vp9_read(r, cm->fc.filterintra_prob[0][mi->bmi[0].as_mode]);
      else
        mi->b_filter_info[0] = mi->b_filter_info[2] = 0;
#endif
      mi->bmi[1].as_mode = mi->bmi[3].as_mode = mbmi->mode =
          read_intra_mode(r, get_y_mode_probs(mi, above_mi, left_mi, 1));
#if CONFIG_FILTERINTRA
      if (is_filter_allowed(mi->bmi[1].as_mode))
        mi->b_filter_info[1] = mi->b_filter_info[3] = mbmi->filterbit =
            vp9_read(r, cm->fc.filterintra_prob[0][mi->bmi[1].as_mode]);
      else
        mi->b_filter_info[1] = mi->b_filter_info[3] = mbmi->filterbit = 0;
#endif
      break;
    case BLOCK_8X4:
      mi->bmi[0].as_mode = mi->bmi[1].as_mode =
          read_intra_mode(r, get_y_mode_probs(mi, above_mi, left_mi, 0));
#if CONFIG_FILTERINTRA
      if (is_filter_allowed(mi->bmi[0].as_mode))
        mi->b_filter_info[0] = mi->b_filter_info[1] =
            vp9_read(r, cm->fc.filterintra_prob[0][mi->bmi[0].as_mode]);
      else
        mi->b_filter_info[0] = mi->b_filter_info[1] = 0;
#endif
      mi->bmi[2].as_mode = mi->bmi[3].as_mode = mbmi->mode =
          read_intra_mode(r, get_y_mode_probs(mi, above_mi, left_mi, 2));
#if CONFIG_FILTERINTRA
      if (is_filter_allowed(mi->bmi[2].as_mode))
        mi->b_filter_info[2] = mi->b_filter_info[3] = mbmi->filterbit =
            vp9_read(r, cm->fc.filterintra_prob[0][mi->bmi[2].as_mode]);
      else
        mi->b_filter_info[2] = mi->b_filter_info[3] = mbmi->filterbit = 0;
#endif
      break;
    default:
#if CONFIG_PALETTE
      if (!mbmi->palette_enabled[0])
        mbmi->mode = read_intra_mode(r,
                       get_y_mode_probs(mi, above_mi, left_mi, 0));
#else
      mbmi->mode = read_intra_mode(r,
                                   get_y_mode_probs(mi, above_mi, left_mi, 0));
#endif  // CONFIG_PALETTE
#if CONFIG_FILTERINTRA
      if (is_filter_enabled(mbmi->tx_size) && is_filter_allowed(mbmi->mode)
#if CONFIG_PALETTE
            && !mbmi->palette_enabled[0]
#endif  // CONFIG_PALETTE
      )
        mbmi->filterbit = vp9_read(r,
                            cm->fc.filterintra_prob[mbmi->tx_size][mbmi->mode]);
      else
        mbmi->filterbit = 0;
#endif  // CONFIG_FILTERINTRA
  }

#if CONFIG_PALETTE
  if (!mbmi->palette_enabled[1])
    mbmi->uv_mode = read_intra_mode(r, vp9_kf_uv_mode_prob[mbmi->mode]);
#else
  mbmi->uv_mode = read_intra_mode(r, vp9_kf_uv_mode_prob[mbmi->mode]);
#endif
#if CONFIG_FILTERINTRA
  if (is_filter_enabled(get_uv_tx_size(mbmi, &xd->plane[1])) &&
      is_filter_allowed(mbmi->uv_mode)
#if CONFIG_PALETTE
            && !mbmi->palette_enabled[1]
#endif  // CONFIG_PALETTE
                                      )
    mbmi->uv_filterbit = vp9_read(r,
        cm->fc.filterintra_prob[get_uv_tx_size(mbmi, &xd->plane[1])][mbmi->uv_mode]);
  else
    mbmi->uv_filterbit = 0;
#endif
}

static int read_mv_component(vp9_reader *r,
                             const nmv_component *mvcomp, int usehp) {
  int mag, d, fr, hp;
  const int sign = vp9_read(r, mvcomp->sign);
  const int mv_class = vp9_read_tree(r, vp9_mv_class_tree, mvcomp->classes);
  const int class0 = mv_class == MV_CLASS_0;

  // Integer part
  if (class0) {
    d = vp9_read_tree(r, vp9_mv_class0_tree, mvcomp->class0);
  } else {
    int i;
    const int n = mv_class + CLASS0_BITS - 1;  // number of bits

    d = 0;
    for (i = 0; i < n; ++i)
      d |= vp9_read(r, mvcomp->bits[i]) << i;
  }

  // Fractional part
  fr = vp9_read_tree(r, vp9_mv_fp_tree, class0 ? mvcomp->class0_fp[d]
                                               : mvcomp->fp);

  // High precision part (if hp is not used, the default value of the hp is 1)
  hp = usehp ? vp9_read(r, class0 ? mvcomp->class0_hp : mvcomp->hp)
             : 1;

  // Result
  mag = vp9_get_mv_mag(mv_class, (d << 3) | (fr << 1) | hp) + 1;
  return sign ? -mag : mag;
}

static INLINE void read_mv(vp9_reader *r, MV *mv, const MV *ref,
                           const nmv_context *ctx,
                           nmv_context_counts *counts, int allow_hp) {
  const MV_JOINT_TYPE joint_type =
      (MV_JOINT_TYPE)vp9_read_tree(r, vp9_mv_joint_tree, ctx->joints);
  const int use_hp = allow_hp && vp9_use_mv_hp(ref);
  MV diff = {0, 0};

  if (mv_joint_vertical(joint_type))
    diff.row = read_mv_component(r, &ctx->comps[0], use_hp);

  if (mv_joint_horizontal(joint_type))
    diff.col = read_mv_component(r, &ctx->comps[1], use_hp);

  vp9_inc_mv(&diff, counts);

  mv->row = ref->row + diff.row;
  mv->col = ref->col + diff.col;
}

static REFERENCE_MODE read_block_reference_mode(VP9_COMMON *cm,
                                                const MACROBLOCKD *xd,
                                                vp9_reader *r) {
  if (cm->reference_mode == REFERENCE_MODE_SELECT) {
    const int ctx = vp9_get_reference_mode_context(cm, xd);
    const REFERENCE_MODE mode =
        (REFERENCE_MODE)vp9_read(r, cm->fc.comp_inter_prob[ctx]);
    if (!cm->frame_parallel_decoding_mode)
      ++cm->counts.comp_inter[ctx][mode];
    return mode;  // SINGLE_REFERENCE or COMPOUND_REFERENCE
  } else {
    return cm->reference_mode;
  }
}

// Read the referncence frame
static void read_ref_frames(VP9_COMMON *const cm, MACROBLOCKD *const xd,
                            vp9_reader *r,
                            int segment_id, MV_REFERENCE_FRAME ref_frame[2]) {
  FRAME_CONTEXT *const fc = &cm->fc;
  FRAME_COUNTS *const counts = &cm->counts;

  if (vp9_segfeature_active(&cm->seg, segment_id, SEG_LVL_REF_FRAME)) {
    ref_frame[0] = (MV_REFERENCE_FRAME)vp9_get_segdata(&cm->seg, segment_id,
                                                       SEG_LVL_REF_FRAME);
    ref_frame[1] = NONE;
  } else {
    const REFERENCE_MODE mode = read_block_reference_mode(cm, xd, r);
    // FIXME(rbultje) I'm pretty sure this breaks segmentation ref frame coding
    if (mode == COMPOUND_REFERENCE) {
      const int idx = cm->ref_frame_sign_bias[cm->comp_fixed_ref];
      const int ctx = vp9_get_pred_context_comp_ref_p(cm, xd);
      const int bit = vp9_read(r, fc->comp_ref_prob[ctx]);
      if (!cm->frame_parallel_decoding_mode)
        ++counts->comp_ref[ctx][bit];
      ref_frame[idx] = cm->comp_fixed_ref;
      ref_frame[!idx] = cm->comp_var_ref[bit];
    } else if (mode == SINGLE_REFERENCE) {
      const int ctx0 = vp9_get_pred_context_single_ref_p1(xd);
      const int bit0 = vp9_read(r, fc->single_ref_prob[ctx0][0]);
      if (!cm->frame_parallel_decoding_mode)
        ++counts->single_ref[ctx0][0][bit0];
      if (bit0) {
        const int ctx1 = vp9_get_pred_context_single_ref_p2(xd);
        const int bit1 = vp9_read(r, fc->single_ref_prob[ctx1][1]);
        if (!cm->frame_parallel_decoding_mode)
          ++counts->single_ref[ctx1][1][bit1];
        ref_frame[0] = bit1 ? ALTREF_FRAME : GOLDEN_FRAME;
      } else {
        ref_frame[0] = LAST_FRAME;
      }

      ref_frame[1] = NONE;
    } else {
      assert(0 && "Invalid prediction mode.");
    }
  }
}


static INLINE INTERP_FILTER read_switchable_interp_filter(
    VP9_COMMON *const cm, MACROBLOCKD *const xd, vp9_reader *r) {
  const int ctx = vp9_get_pred_context_switchable_interp(xd);
  const INTERP_FILTER type =
      (INTERP_FILTER)vp9_read_tree(r, vp9_switchable_interp_tree,
                                   cm->fc.switchable_interp_prob[ctx]);
  if (!cm->frame_parallel_decoding_mode)
    ++cm->counts.switchable_interp[ctx][type];
  return type;
}

static void read_intra_block_mode_info(VP9_COMMON *const cm, MODE_INFO *mi,
#if CONFIG_FILTERINTRA
                                       MACROBLOCKD *const xd,
#endif
                                       vp9_reader *r) {
  MB_MODE_INFO *const mbmi = &mi->mbmi;
  const BLOCK_SIZE bsize = mi->mbmi.sb_type;
  int i;

  mbmi->ref_frame[0] = INTRA_FRAME;
  mbmi->ref_frame[1] = NONE;

  switch (bsize) {
    case BLOCK_4X4:
#if CONFIG_FILTERINTRA
      for (i = 0; i < 4; ++i) {
#else
      for (i = 0; i < 4; ++i)
#endif
        mi->bmi[i].as_mode = read_intra_mode_y(cm, r, 0);
#if CONFIG_FILTERINTRA
        if (is_filter_allowed(mi->bmi[i].as_mode)) {
          mi->b_filter_info[i] =
              vp9_read(r, cm->fc.filterintra_prob[0][mi->bmi[i].as_mode]);
          cm->counts.filterintra[0][mi->bmi[i].as_mode]
                                   [mi->b_filter_info[i]]++;
        } else {
          mi->b_filter_info[i] = 0;
        }
      }
      mbmi->filterbit = mi->b_filter_info[3];
#endif
      mbmi->mode = mi->bmi[3].as_mode;
      break;
    case BLOCK_4X8:
      mi->bmi[0].as_mode = mi->bmi[2].as_mode = read_intra_mode_y(cm, r, 0);
#if CONFIG_FILTERINTRA
      if (is_filter_allowed(mi->bmi[0].as_mode)) {
        mi->b_filter_info[0] = mi->b_filter_info[2] =
            vp9_read(r, cm->fc.filterintra_prob[0][mi->bmi[0].as_mode]);
        cm->counts.filterintra[0][mi->bmi[0].as_mode][mi->b_filter_info[0]]++;
      } else {
        mi->b_filter_info[0] = mi->b_filter_info[2] = 0;
      }
#endif
      mi->bmi[1].as_mode = mi->bmi[3].as_mode = mbmi->mode =
          read_intra_mode_y(cm, r, 0);
#if CONFIG_FILTERINTRA
      if (is_filter_allowed(mi->bmi[1].as_mode)) {
        mi->b_filter_info[1] = mi->b_filter_info[3] = mbmi->filterbit =
            vp9_read(r, cm->fc.filterintra_prob[0][mi->bmi[1].as_mode]);
        cm->counts.filterintra[0][mi->bmi[1].as_mode][mi->b_filter_info[1]]++;
      } else {
        mi->b_filter_info[1] = mi->b_filter_info[3] = mbmi->filterbit = 0;
      }
#endif
      break;
    case BLOCK_8X4:
      mi->bmi[0].as_mode = mi->bmi[1].as_mode = read_intra_mode_y(cm, r, 0);
#if CONFIG_FILTERINTRA
      if (is_filter_allowed(mi->bmi[0].as_mode)) {
        mi->b_filter_info[0] = mi->b_filter_info[1] =
            vp9_read(r, cm->fc.filterintra_prob[0][mi->bmi[0].as_mode]);
        cm->counts.filterintra[0][mi->bmi[0].as_mode][mi->b_filter_info[0]]++;
      } else {
        mi->b_filter_info[0] = mi->b_filter_info[1] = 0;
      }
#endif
      mi->bmi[2].as_mode = mi->bmi[3].as_mode = mbmi->mode =
          read_intra_mode_y(cm, r, 0);
#if CONFIG_FILTERINTRA
      if (is_filter_allowed(mi->bmi[2].as_mode)) {
        mi->b_filter_info[2] = mi->b_filter_info[3] = mbmi->filterbit =
            vp9_read(r, cm->fc.filterintra_prob[0][mi->bmi[2].as_mode]);
        cm->counts.filterintra[0][mi->bmi[2].as_mode][mi->b_filter_info[2]]++;
      } else {
        mi->b_filter_info[2] = mi->b_filter_info[3] = mbmi->filterbit = 0;
      }
#endif
      break;
    default:
#if CONFIG_PALETTE
      if (!mbmi->palette_enabled[0]) {
        mbmi->mode = read_intra_mode_y(cm, r, size_group_lookup[bsize]);
      } else {
        mbmi->mode = DC_PRED;
        if (!cm->frame_parallel_decoding_mode)
            ++cm->counts.y_mode[size_group_lookup[bsize]][DC_PRED];
      }
#else
      mbmi->mode = read_intra_mode_y(cm, r, size_group_lookup[bsize]);
#endif  // CONFIG_PALETTE
#if CONFIG_FILTERINTRA
      if (is_filter_allowed(mbmi->mode) && is_filter_enabled(mbmi->tx_size)
#if CONFIG_PALETTE
          && !mbmi->palette_enabled[0]
#endif  // CONFIG_PALETTE
      ) {
        mbmi->filterbit = vp9_read(r,
            cm->fc.filterintra_prob[mbmi->tx_size][mbmi->mode]);
        cm->counts.filterintra[mbmi->tx_size][mbmi->mode][mbmi->filterbit]++;
      } else {
        mbmi->filterbit = 0;
#if CONFIG_PALETTE
        if (mbmi->palette_enabled[0])
          cm->counts.filterintra[mbmi->tx_size][mbmi->mode][mbmi->filterbit]++;
#endif  // CONFIG_PALETTE
      }
#endif  // CONFIG_FILTERINTRA
  }

#if CONFIG_PALETTE
  if (!mbmi->palette_enabled[1]) {
    mbmi->uv_mode = read_intra_mode_uv(cm, r, mbmi->mode);
  } else {
    mbmi->uv_mode = DC_PRED;
    if (!cm->frame_parallel_decoding_mode)
      ++cm->counts.uv_mode[mbmi->mode][DC_PRED];
  }
#else
  mbmi->uv_mode = read_intra_mode_uv(cm, r, mbmi->mode);
#endif  // CONFIG_PALETTE
#if CONFIG_FILTERINTRA
  if (is_filter_allowed(mbmi->uv_mode) &&
      is_filter_enabled(get_uv_tx_size(mbmi, &xd->plane[1]))
#if CONFIG_PALETTE
      && !mbmi->palette_enabled[1]
#endif  // CONFIG_PALETTE
  ) {
    mbmi->uv_filterbit = vp9_read(r,
        cm->fc.filterintra_prob[get_uv_tx_size(mbmi, &xd->plane[1])][mbmi->uv_mode]);
    cm->counts.filterintra[get_uv_tx_size(mbmi, &xd->plane[1])]
                           [mbmi->uv_mode][mbmi->uv_filterbit]++;
  } else {
    mbmi->uv_filterbit = 0;
#if CONFIG_PALETTE
    if (mbmi->palette_enabled[1])
      cm->counts.filterintra[get_uv_tx_size(mbmi, &xd->plane[1])]
                             [mbmi->uv_mode][mbmi->uv_filterbit]++;
#endif  // CONFIG_PALETTE
  }
#endif  // CONFIG_FILTERINTRA
}

static INLINE int is_mv_valid(const MV *mv) {
  return mv->row > MV_LOW && mv->row < MV_UPP &&
         mv->col > MV_LOW && mv->col < MV_UPP;
}

static INLINE int assign_mv(VP9_COMMON *cm, PREDICTION_MODE mode,
                            int_mv mv[2], int_mv ref_mv[2],
                            int_mv nearest_mv[2], int_mv near_mv[2],
                            int is_compound, int allow_hp, vp9_reader *r) {
  int i;
  int ret = 1;
#if CONFIG_COMPOUND_MODES
  assert(is_inter_mode(mode) || is_inter_compound_mode(mode));
#else
  assert(is_inter_mode(mode));
#endif
  switch (mode) {
    case NEWMV: {
      nmv_context_counts *const mv_counts = cm->frame_parallel_decoding_mode ?
                                            NULL : &cm->counts.mv;
      for (i = 0; i < 1 + is_compound; ++i) {
        read_mv(r, &mv[i].as_mv, &ref_mv[i].as_mv, &cm->fc.nmvc, mv_counts,
                allow_hp);
        ret = ret && is_mv_valid(&mv[i].as_mv);
        assert(ret);
      }
      break;
    }
    case NEARESTMV: {
      mv[0].as_int = nearest_mv[0].as_int;
      if (is_compound)
        mv[1].as_int = nearest_mv[1].as_int;
      break;
    }
    case NEARMV: {
      mv[0].as_int = near_mv[0].as_int;
      if (is_compound)
        mv[1].as_int = near_mv[1].as_int;
      break;
    }
    case ZEROMV: {
      mv[0].as_int = 0;
      if (is_compound)
        mv[1].as_int = 0;
      break;
    }
#if CONFIG_COMPOUND_MODES
    case NEW_NEWMV: {
      nmv_context_counts *const mv_counts = cm->frame_parallel_decoding_mode ?
                                            NULL : &cm->counts.mv;
      assert(is_compound);
      for (i = 0; i < 2; ++i) {
        read_mv(r, &mv[i].as_mv, &ref_mv[i].as_mv, &cm->fc.nmvc, mv_counts,
                allow_hp);
        ret = ret && is_mv_valid(&mv[i].as_mv);
      }
      break;
    }
    case NEAREST_NEARESTMV: {
      assert(is_compound);
      mv[0].as_int = nearest_mv[0].as_int;
      mv[1].as_int = nearest_mv[1].as_int;
      break;
    }
    case NEAREST_NEARMV: {
      assert(is_compound);
      mv[0].as_int = nearest_mv[0].as_int;
      mv[1].as_int = near_mv[1].as_int;
      break;
    }
    case NEAR_NEARESTMV: {
      assert(is_compound);
      mv[0].as_int = near_mv[0].as_int;
      mv[1].as_int = nearest_mv[1].as_int;
      break;
    }
    case NEW_NEARESTMV: {
      nmv_context_counts *const mv_counts = cm->frame_parallel_decoding_mode ?
                                            NULL : &cm->counts.mv;
      assert(is_compound);
      read_mv(r, &mv[0].as_mv, &ref_mv[0].as_mv, &cm->fc.nmvc, mv_counts,
              allow_hp);
      ret = ret && is_mv_valid(&mv[0].as_mv);
      mv[1].as_int = nearest_mv[1].as_int;
      break;
    }
    case NEAREST_NEWMV: {
      nmv_context_counts *const mv_counts = cm->frame_parallel_decoding_mode ?
                                            NULL : &cm->counts.mv;
      assert(is_compound);
      mv[0].as_int = nearest_mv[0].as_int;
      read_mv(r, &mv[1].as_mv, &ref_mv[1].as_mv, &cm->fc.nmvc, mv_counts,
              allow_hp);
      ret = ret && is_mv_valid(&mv[1].as_mv);
      break;
    }
    case NEAR_NEWMV: {
      nmv_context_counts *const mv_counts = cm->frame_parallel_decoding_mode ?
                        NULL : &cm->counts.mv;
      assert(is_compound);
      mv[0].as_int = near_mv[0].as_int;
      read_mv(r, &mv[1].as_mv, &ref_mv[1].as_mv, &cm->fc.nmvc, mv_counts,
          allow_hp);
      ret = ret && is_mv_valid(&mv[1].as_mv);
      break;
    }
    case NEW_NEARMV: {
      nmv_context_counts *const mv_counts = cm->frame_parallel_decoding_mode ?
                        NULL : &cm->counts.mv;
      assert(is_compound);
      read_mv(r, &mv[0].as_mv, &ref_mv[0].as_mv, &cm->fc.nmvc, mv_counts,
          allow_hp);
      ret = ret && is_mv_valid(&mv[0].as_mv);
      mv[1].as_int = near_mv[1].as_int;
      break;
    }
    case ZERO_ZEROMV: {
      assert(is_compound);
      mv[0].as_int = 0;
      mv[1].as_int = 0;
      break;
    }
#endif  // CONFIG_COMPOUND_MODES
    default: {
      return 0;
    }
  }

  return ret;
}

static int read_is_inter_block(VP9_COMMON *const cm, MACROBLOCKD *const xd,
                               int segment_id, vp9_reader *r) {
  if (vp9_segfeature_active(&cm->seg, segment_id, SEG_LVL_REF_FRAME)) {
    return vp9_get_segdata(&cm->seg, segment_id, SEG_LVL_REF_FRAME) !=
           INTRA_FRAME;
  } else {
    const int ctx = vp9_get_intra_inter_context(xd);
    const int is_inter = vp9_read(r, cm->fc.intra_inter_prob[ctx]);
    if (!cm->frame_parallel_decoding_mode)
      ++cm->counts.intra_inter[ctx][is_inter];
    return is_inter;
  }
}

static void read_inter_block_mode_info(VP9_COMMON *const cm,
                                       MACROBLOCKD *const xd,
                                       const TileInfo *const tile,
                                       MODE_INFO *const mi,
#if CONFIG_SUPERTX
                                       int supertx_enabled,
#endif
                                       int mi_row, int mi_col, vp9_reader *r) {
  MB_MODE_INFO *const mbmi = &mi->mbmi;
  const BLOCK_SIZE bsize = mbmi->sb_type;
  const int allow_hp = cm->allow_high_precision_mv;

  int_mv nearestmv[2], nearmv[2];
  int inter_mode_ctx, ref, is_compound;
#if CONFIG_SUPERTX
  (void) supertx_enabled;
#endif

#if CONFIG_COPY_MODE
  if (mbmi->copy_mode == NOREF)
#endif
    read_ref_frames(cm, xd, r, mbmi->segment_id, mbmi->ref_frame);
  is_compound = has_second_ref(mbmi);

  for (ref = 0; ref < 1 + is_compound; ++ref) {
    const MV_REFERENCE_FRAME frame = mbmi->ref_frame[ref];
    RefBuffer *ref_buf = &cm->frame_refs[frame - LAST_FRAME];
    xd->block_refs[ref] = ref_buf;
    if ((!vp9_is_valid_scale(&ref_buf->sf)))
      vpx_internal_error(&cm->error, VPX_CODEC_UNSUP_BITSTREAM,
                         "Reference frame has invalid dimensions");
    if (ref_buf->buf->corrupted)
      vpx_internal_error(&cm->error, VPX_CODEC_CORRUPT_FRAME,
                         "Block reference is corrupt");
    vp9_setup_pre_planes(xd, ref, ref_buf->buf, mi_row, mi_col,
                         &ref_buf->sf);
#if CONFIG_COPY_MODE
    if (mbmi->copy_mode == NOREF)
#endif
      vp9_find_mv_refs(cm, xd, tile, mi, frame, mbmi->ref_mvs[frame],
                       mi_row, mi_col);
  }
#if CONFIG_COPY_MODE
  if (mbmi->copy_mode != NOREF)
    return;
#endif

  inter_mode_ctx = mbmi->mode_context[mbmi->ref_frame[0]];
  if (vp9_segfeature_active(&cm->seg, mbmi->segment_id, SEG_LVL_SKIP)) {
    mbmi->mode = ZEROMV;
    if (bsize < BLOCK_8X8) {
        vpx_internal_error(&cm->error, VPX_CODEC_UNSUP_BITSTREAM,
                           "Invalid usage of segement feature on small blocks");
        return;
    }
  } else {
    if (bsize >= BLOCK_8X8) {
#if CONFIG_COMPOUND_MODES
      if (is_compound) {
        mbmi->mode = read_inter_compound_mode(cm, r, inter_mode_ctx);
      } else {
        mbmi->mode = read_inter_mode(cm, r, inter_mode_ctx);
      }
#else
      mbmi->mode = read_inter_mode(cm, r, inter_mode_ctx);
#endif
    }
  }

#if CONFIG_COMPOUND_MODES
  if (bsize < BLOCK_8X8 ||
      (mbmi->mode != ZEROMV && mbmi->mode != ZERO_ZEROMV)) {
#else
  if (bsize < BLOCK_8X8 || mbmi->mode != ZEROMV) {
#endif
    for (ref = 0; ref < 1 + is_compound; ++ref) {
      vp9_find_best_ref_mvs(xd, allow_hp, mbmi->ref_mvs[mbmi->ref_frame[ref]],
                            &nearestmv[ref], &nearmv[ref]);
    }
  }

  mbmi->interp_filter = (cm->interp_filter == SWITCHABLE)
                      ? read_switchable_interp_filter(cm, xd, r)
                      : cm->interp_filter;

#if CONFIG_INTERINTRA
    if (is_interintra_allowed(bsize) &&
        is_inter_mode(mbmi->mode) &&
#if CONFIG_SUPERTX
        !supertx_enabled &&
#endif
        mbmi->ref_frame[1] <= INTRA_FRAME) {
      mbmi->ref_frame[1] = vp9_read(r, cm->fc.interintra_prob[bsize]) ?
                           INTRA_FRAME : NONE;
      cm->counts.interintra[bsize][mbmi->ref_frame[1] == INTRA_FRAME]++;
#if CONFIG_WEDGE_PARTITION
      mbmi->use_wedge_interintra = 0;
#endif  // CONFIG_WEDGE_PARTITION
      if (mbmi->ref_frame[1] == INTRA_FRAME) {
        mbmi->interintra_mode =
            read_intra_mode_y(cm, r, size_group_lookup[bsize]);
        mbmi->interintra_uv_mode = mbmi->interintra_mode;
#if CONFIG_WEDGE_PARTITION
        if (get_wedge_bits(bsize)) {
          mbmi->use_wedge_interintra = vp9_read(
              r, cm->fc.wedge_interintra_prob[bsize]);
          cm->counts.wedge_interintra[bsize][mbmi->use_wedge_interintra]++;
          if (mbmi->use_wedge_interintra) {
            mbmi->interintra_wedge_index = vp9_read_literal(
                r, get_wedge_bits(bsize));
            mbmi->interintra_uv_wedge_index = mbmi->interintra_wedge_index;
          }
        }
#endif  // CONFIG_WEDGE_PARTITION
      }
    }
#endif  // CONFIG_INTERINTRA

  if (bsize < BLOCK_8X8) {
    const int num_4x4_w = num_4x4_blocks_wide_lookup[bsize];  // 1 or 2
    const int num_4x4_h = num_4x4_blocks_high_lookup[bsize];  // 1 or 2
    int idx, idy;
    PREDICTION_MODE b_mode;
    int_mv nearest_sub8x8[2], near_sub8x8[2];
#if CONFIG_NEWMVREF_SUB8X8
    int_mv ref_mv[2];
#endif  // CONFIG_NEWMVREF_SUB8X8
    for (idy = 0; idy < 2; idy += num_4x4_h) {
      for (idx = 0; idx < 2; idx += num_4x4_w) {
        int_mv block[2];
        const int j = idy * 2 + idx;
#if CONFIG_COMPOUND_MODES
        if (is_compound) {
          b_mode = read_inter_compound_mode(cm, r, inter_mode_ctx);
        } else {
          b_mode = read_inter_mode(cm, r, inter_mode_ctx);
        }
#else
        b_mode = read_inter_mode(cm, r, inter_mode_ctx);
#endif
#if CONFIG_COMPOUND_MODES
        if (b_mode == NEARESTMV || b_mode == NEARMV ||
#if CONFIG_NEWMVREF_SUB8X8
            b_mode == NEWMV || b_mode == NEW_NEWMV ||
#endif  // CONFIG_NEWMVREF_SUB8X8
            b_mode == NEAREST_NEARESTMV || b_mode == NEAREST_NEARMV ||
            b_mode == NEAR_NEARESTMV || b_mode == NEAREST_NEWMV ||
            b_mode == NEW_NEARESTMV || b_mode == NEAR_NEWMV ||
            b_mode == NEW_NEARMV)
#else
        if (b_mode == NEARESTMV || b_mode == NEARMV
#if CONFIG_NEWMVREF_SUB8X8
            || b_mode == NEWMV
#endif  // CONFIG_NEWMVREF_SUB8X8
            )
#endif  // CONFIG_COMPOUND_MODES
        {
          for (ref = 0; ref < 1 + is_compound; ++ref) {
#if CONFIG_NEWMVREF_SUB8X8
            int_mv mv_ref_list[MAX_MV_REF_CANDIDATES];
            int_mv second_ref_mv;
            vp9_update_mv_context(cm, xd, tile, mi, mbmi->ref_frame[ref],
                                  mv_ref_list, j, mi_row, mi_col);
#endif  // CONFIG_NEWMVREF_SUB8X8
            vp9_append_sub8x8_mvs_for_idx(cm, xd, tile, j, ref, mi_row, mi_col,
#if CONFIG_NEWMVREF_SUB8X8
                                          mv_ref_list,
#endif  // CONFIG_NEWMVREF_SUB8X8
                                          &nearest_sub8x8[ref],
                                          &near_sub8x8[ref]);
#if CONFIG_NEWMVREF_SUB8X8
            if (b_mode == NEWMV
#if CONFIG_COMPOUND_MODES
                || b_mode == NEW_NEWMV ||
                b_mode == NEAREST_NEWMV ||
                b_mode == NEW_NEARESTMV ||
                b_mode == NEAR_NEWMV ||
                b_mode == NEW_NEARMV
#endif  // CONFIG_COMPOUND_MODES
                ) {
              mv_ref_list[0].as_int = nearest_sub8x8[ref].as_int;
              mv_ref_list[1].as_int = near_sub8x8[ref].as_int;
              vp9_find_best_ref_mvs(xd, allow_hp, mv_ref_list,
                                    &ref_mv[ref], &second_ref_mv);
            }
#endif  // CONFIG_NEWMVREF_SUB8X8
          }
        }

        if (!assign_mv(cm, b_mode, block,
#if CONFIG_NEWMVREF_SUB8X8
                       ref_mv,
#else
                       nearestmv,
#endif  // CONFIG_NEWMVREF_SUB8X8
                       nearest_sub8x8, near_sub8x8,
                       is_compound, allow_hp, r)) {
          xd->corrupted |= 1;
          break;
        }

        mi->bmi[j].as_mv[0].as_int = block[0].as_int;
        if (is_compound)
          mi->bmi[j].as_mv[1].as_int = block[1].as_int;

        if (num_4x4_h == 2)
          mi->bmi[j + 2] = mi->bmi[j];
        if (num_4x4_w == 2)
          mi->bmi[j + 1] = mi->bmi[j];
      }
    }

    mi->mbmi.mode = b_mode;
    mbmi->mv[0].as_int = mi->bmi[3].as_mv[0].as_int;
    mbmi->mv[1].as_int = mi->bmi[3].as_mv[1].as_int;
  } else {
    xd->corrupted |= !assign_mv(cm, mbmi->mode, mbmi->mv, nearestmv,
                                nearestmv, nearmv, is_compound, allow_hp, r);
  }
#if CONFIG_TX_SKIP
  mbmi->uv_mode = mbmi->mode;
#endif
#if CONFIG_WEDGE_PARTITION
  mbmi->use_wedge_interinter = 0;
  if (cm->reference_mode != SINGLE_REFERENCE &&
#if CONFIG_COMPOUND_MODES
      is_inter_compound_mode(mbmi->mode) &&
#endif  // CONFIG_COMPOUND_MODES
      get_wedge_bits(bsize) &&
      mbmi->ref_frame[1] > INTRA_FRAME) {
    mbmi->use_wedge_interinter =
        vp9_read(r, cm->fc.wedge_interinter_prob[bsize]);
    cm->counts.wedge_interinter[bsize][mbmi->use_wedge_interinter]++;
    if (mbmi->use_wedge_interinter) {
      mbmi->interinter_wedge_index = vp9_read_literal(r, get_wedge_bits(bsize));
    }
  }
#endif  // CONFIG_WEDGE_PARTITION
}

static void read_inter_frame_mode_info(VP9_COMMON *const cm,
                                       MACROBLOCKD *const xd,
                                       const TileInfo *const tile,
#if CONFIG_SUPERTX
                                       int supertx_enabled,
#endif
                                       int mi_row, int mi_col, vp9_reader *r) {
  MODE_INFO *const mi = xd->mi[0].src_mi;
  MB_MODE_INFO *const mbmi = &mi->mbmi;
  int inter_block;
#if CONFIG_COPY_MODE
  int num_candidate = 0;
  MB_MODE_INFO *inter_ref_list[18] = {NULL};
#endif
#if CONFIG_SUPERTX
  (void) supertx_enabled;
#endif

  mbmi->mv[0].as_int = 0;
  mbmi->mv[1].as_int = 0;

#if CONFIG_COPY_MODE
  if (mbmi->sb_type >= BLOCK_8X8)
    num_candidate = vp9_construct_ref_inter_list(
        cm, xd, mbmi->sb_type, mi_row, mi_col, inter_ref_list);
  if (mbmi->sb_type >= BLOCK_8X8 && num_candidate > 0) {
    int ctx = vp9_get_copy_mode_context(xd);
    int is_copy = vp9_read(r, cm->fc.copy_noref_prob[ctx][mbmi->sb_type]);

    ++cm->counts.copy_noref[ctx][mbmi->sb_type][is_copy];
    if (!is_copy) {
      mbmi->copy_mode = NOREF;
    } else {
      mbmi->copy_mode = read_copy_mode(cm, r, num_candidate, ctx);
    }
  } else {
    mbmi->copy_mode = NOREF;
  }
  if (mbmi->copy_mode != NOREF) {
    BLOCK_SIZE bsize_backup = mbmi->sb_type;
    int skip_backup = mbmi->skip;
    COPY_MODE copy_mode_backup = mbmi->copy_mode;
#if CONFIG_SUPERTX
    TX_SIZE tx_size_backup = mbmi->tx_size;
#endif  // CONFIG_SUPERTX
#if CONFIG_EXT_TX
    EXT_TX_TYPE ext_txfrm_backup = mbmi->ext_txfrm;
#endif  // CONFIG_EXT_TX

    inter_block = 1;
    *mbmi = *inter_ref_list[mbmi->copy_mode - REF0];
#if CONFIG_SUPERTX
    mbmi->tx_size = tx_size_backup;
#endif  // CONFIG_SUPERTX
#if CONFIG_EXT_TX
    mbmi->ext_txfrm = ext_txfrm_backup;
#endif  // CONFIG_EXT_TX
#if CONFIG_INTERINTRA
    if (mbmi->ref_frame[1] == INTRA_FRAME)
      mbmi->ref_frame[1] = NONE;
#endif  // CONFIG_INTERINTRA
#if CONFIG_WEDGE_PARTITION
    mbmi->use_wedge_interinter = 0;
#endif  // CONFIG_WEDGE_PARTITION
    mbmi->sb_type = bsize_backup;
    mbmi->mode = NEARESTMV;
    mbmi->skip = skip_backup;
    mbmi->copy_mode = copy_mode_backup;
  }
#endif  // CONFIG_COPY_MODE

  mbmi->segment_id = read_inter_segment_id(cm, xd, mi_row, mi_col, r);
#if CONFIG_SUPERTX
  if (!supertx_enabled) {
#endif
    mbmi->skip = read_skip(cm, xd, mbmi->segment_id, r);
#if CONFIG_COPY_MODE
    if (mbmi->copy_mode == NOREF)
#endif
      inter_block = read_is_inter_block(cm, xd, mbmi->segment_id, r);

#if CONFIG_PALETTE
    mbmi->palette_enabled[0] = 0;
    mbmi->palette_enabled[1] = 0;

    if (!inter_block && mbmi->sb_type >= BLOCK_8X8 && cm->allow_palette_mode) {
      const MODE_INFO *above_mi = xd->up_available ?
          xd->mi[-xd->mi_stride].src_mi : NULL;
      const MODE_INFO *left_mi = xd->left_available ?
          xd->mi[-1].src_mi : NULL;
      int ctx = 0;

      if (above_mi)
        ctx += (above_mi->mbmi.palette_enabled[0] == 1);
      if (left_mi)
        ctx += (left_mi->mbmi.palette_enabled[0] == 1);
      mbmi->palette_enabled[0] =
          vp9_read(r,
                   cm->fc.palette_enabled_prob[mbmi->sb_type - BLOCK_8X8][ctx]);
      mbmi->palette_enabled[1] =
          vp9_read(r, cm->fc.palette_uv_enabled_prob[mbmi->palette_enabled[0]]);
    }

    if (mbmi->palette_enabled[0]) {
      BLOCK_SIZE bsize = mbmi->sb_type;
      int i, j, n, color_ctx, color_idx;
      int rows = 4 * num_4x4_blocks_high_lookup[bsize];
      int cols = 4 * num_4x4_blocks_wide_lookup[bsize];
      int color_order[PALETTE_MAX_SIZE];
      uint8_t *color_map = xd->plane[0].color_index_map;

      mbmi->mode = DC_PRED;
      mbmi->palette_size[0] =
          vp9_read_tree(r, vp9_palette_size_tree,
                        cm->fc.palette_size_prob[bsize - BLOCK_8X8]);
      mbmi->palette_size[0] += 2;
      n = mbmi->palette_size[0];

      for (i = 0; i < mbmi->palette_size[0]; i++)
        mbmi->palette_colors[i] = vp9_read_literal(r, 8);

      color_map[0] = vp9_read_literal(r,
                                      vp9_get_bit_depth(mbmi->palette_size[0]));
      for (i = 0; i < rows; i++) {
        for (j = (i == 0 ? 1 : 0); j < cols; j++) {
          color_ctx = vp9_get_palette_color_context(color_map, cols, i, j, n,
                                                    color_order);
          color_idx = vp9_read_tree(r, vp9_palette_color_tree,
                                    cm->fc.palette_color_prob[n - 2]
                                                              [color_ctx]);

          color_map[i * cols + j] = color_order[color_idx];
        }
      }

      mbmi->tx_size = MIN(max_txsize_lookup[bsize],
                          tx_mode_to_biggest_tx_size[cm->tx_mode]);
      if (!cm->frame_parallel_decoding_mode)
        ++get_tx_counts(max_txsize_lookup[bsize], vp9_get_tx_size_context(xd),
                        &cm->counts.tx)[mbmi->tx_size];
    }

    if (mbmi->palette_enabled[1]) {
      int i, j;
      BLOCK_SIZE bsize = mbmi->sb_type;
      int rows = 4 * num_4x4_blocks_high_lookup[bsize] >>
          xd->plane[1].subsampling_y;
      int cols = 4 * num_4x4_blocks_wide_lookup[bsize] >>
          xd->plane[1].subsampling_x;

      mbmi->uv_mode = DC_PRED;
      if (xd->plane[1].subsampling_x && xd->plane[1].subsampling_y) {
        mbmi->palette_size[1] =
            vp9_read_tree(r, vp9_palette_size_tree,
                          cm->fc.palette_uv_size_prob[bsize - BLOCK_8X8]);
        mbmi->palette_size[1] += 2;
      } else {
        mbmi->palette_size[1] = mbmi->palette_size[0];
      }

      for (i = 0; i < mbmi->palette_size[1]; i++)
        mbmi->palette_colors[PALETTE_MAX_SIZE + i] = vp9_read_literal(r, 8);
      for (i = 0; i < mbmi->palette_size[1]; i++)
        mbmi->palette_colors[2 * PALETTE_MAX_SIZE + i] = vp9_read_literal(r, 8);

      if (xd->plane[1].subsampling_x && xd->plane[1].subsampling_y) {
        int color_idx = 0, color_ctx = 0;
        int n = mbmi->palette_size[1];
        int color_order[PALETTE_MAX_SIZE];
        uint8_t *color_map = xd->plane[1].color_index_map;

        color_map[0] = vp9_read_literal(r, vp9_get_bit_depth(n));
        for (i = 0; i < rows; i++) {
          for (j = (i == 0 ? 1 : 0); j < cols; j++) {
            color_ctx = vp9_get_palette_color_context(color_map, cols, i, j, n,
                                                      color_order);
            color_idx = vp9_read_tree(r, vp9_palette_color_tree,
                                      cm->fc.palette_uv_color_prob[n - 2]
                                                                   [color_ctx]);
            color_map[i * cols + j] = color_order[color_idx];
          }
        }
      }
    }

    if (!inter_block && mbmi->sb_type >= BLOCK_8X8 && cm->allow_palette_mode) {
      BLOCK_SIZE bsize = mbmi->sb_type;
      int palette_ctx = 0;
      const MODE_INFO *above_mi = xd->up_available ?
          xd->mi[-xd->mi_stride].src_mi : NULL;
      const MODE_INFO *left_mi = xd->left_available ?
          xd->mi[-1].src_mi : NULL;

      if (above_mi)
        palette_ctx += (above_mi->mbmi.palette_enabled[0] == 1);
      if (left_mi)
        palette_ctx += (left_mi->mbmi.palette_enabled[0] == 1);
      vp9_update_palette_counts(&cm->counts, mbmi, bsize, palette_ctx);
    }

    if (!mbmi->palette_enabled[0]) {
      mbmi->tx_size = read_tx_size(cm, xd, cm->tx_mode, mbmi->sb_type,
                                   !mbmi->skip || !inter_block, r);
    }
#else
    mbmi->tx_size = read_tx_size(cm, xd, cm->tx_mode, mbmi->sb_type,
                                     !mbmi->skip || !inter_block, r);
#endif  // CONFIG_PALETTE

#if CONFIG_EXT_TX
    if (inter_block &&
        mbmi->tx_size <= TX_16X16 &&
        cm->base_qindex > 0 &&
        mbmi->sb_type >= BLOCK_8X8 &&
#if CONFIG_SUPERTX
      !supertx_enabled &&
#endif
      !vp9_segfeature_active(&cm->seg, mbmi->segment_id, SEG_LVL_SKIP) &&
      !mbmi->skip) {
      mbmi->ext_txfrm = vp9_read_tree(r, vp9_ext_tx_tree,
                                      cm->fc.ext_tx_prob[mbmi->tx_size]);
      if (!cm->frame_parallel_decoding_mode)
        ++cm->counts.ext_tx[mbmi->tx_size][mbmi->ext_txfrm];
    } else {
      mbmi->ext_txfrm = NORM;
    }
#endif  // CONFIG_EXT_TX
#if CONFIG_SUPERTX
  } else {
    const int ctx = vp9_get_intra_inter_context(xd);
    inter_block = 1;
    if (!cm->frame_parallel_decoding_mode)
#if CONFIG_COPY_MODE
      if (mbmi->copy_mode == NOREF)
#endif  // CONFIG_COPY_MODE
        ++cm->counts.intra_inter[ctx][1];
#if CONFIG_PALETTE
    mbmi->palette_enabled[0] = 0;
    mbmi->palette_enabled[1] = 0;
#endif  // CONFIG_PALETTE
  }
#endif  // CONFIG_SUPERTX

#if CONFIG_TX_SKIP
  if (mbmi->sb_type >= BLOCK_8X8) {
    int q_idx = vp9_get_qindex(&cm->seg, mbmi->segment_id, cm->base_qindex);
    int try_tx_skip = inter_block ? q_idx <= TX_SKIP_Q_THRESH_INTER :
                                    q_idx <= TX_SKIP_Q_THRESH_INTRA;
#if CONFIG_COPY_MODE
    if (mbmi->copy_mode != NOREF)
      try_tx_skip = 0;
#endif  // CONFIG_COPY_MODE

#if CONFIG_SUPERTX
    if (try_tx_skip && !supertx_enabled) {
#else
    if (try_tx_skip && (!mbmi->skip || !inter_block)) {
#endif  // CONFIG_SUPERTX
      if (xd->lossless) {
#if CONFIG_SUPERTX
        if (1)
#else
        if (mbmi->tx_size == TX_4X4)
#endif  // CONFIG_SUPERTX
          mbmi->tx_skip[0] = vp9_read(r, cm->fc.y_tx_skip_prob[inter_block]);
        else
          mbmi->tx_skip[0] = 1;

#if CONFIG_SUPERTX
        if (1)
#else
        if (get_uv_tx_size(mbmi, &xd->plane[1]) == TX_4X4)
#endif  // CONFIG_SUPERTX
          mbmi->tx_skip[1] =
            vp9_read(r, cm->fc.uv_tx_skip_prob[mbmi->tx_skip[0]]);
        else
          mbmi->tx_skip[1] = 1;
      } else {
        mbmi->tx_skip[0] = vp9_read(r, cm->fc.y_tx_skip_prob[inter_block]);
        mbmi->tx_skip[1] =
          vp9_read(r, cm->fc.uv_tx_skip_prob[mbmi->tx_skip[0]]);
      }
#if CONFIG_SUPERTX
      if (!cm->frame_parallel_decoding_mode && !supertx_enabled) {
#else
      if (!cm->frame_parallel_decoding_mode) {
#endif
        ++cm->counts.y_tx_skip[inter_block][mbmi->tx_skip[0]];
        ++cm->counts.uv_tx_skip[mbmi->tx_skip[0]][mbmi->tx_skip[1]];
      }
    } else {
      mbmi->tx_skip[0] = 0;
      mbmi->tx_skip[1] = 0;
    }
  } else {
    mbmi->tx_skip[0] = 0;
    mbmi->tx_skip[1] = 0;
  }
#endif  // CONFIG_TX_SKIP

  if (inter_block) {
    read_inter_block_mode_info(cm, xd, tile, mi,
#if CONFIG_SUPERTX
                               supertx_enabled,
#endif
                               mi_row, mi_col, r);
  } else {
    read_intra_block_mode_info(cm, mi,
#if CONFIG_FILTERINTRA
                               xd,
#endif
                               r);
  }
}

void vp9_read_mode_info(VP9_COMMON *cm, MACROBLOCKD *xd,
                        const TileInfo *const tile,
#if CONFIG_SUPERTX
                        int supertx_enabled,
#endif
                        int mi_row, int mi_col, vp9_reader *r) {
  if (frame_is_intra_only(cm))
    read_intra_frame_mode_info(cm, xd, mi_row, mi_col, r);
  else
    read_inter_frame_mode_info(cm, xd, tile,
#if CONFIG_SUPERTX
                               supertx_enabled,
#endif
                               mi_row, mi_col, r);
}
