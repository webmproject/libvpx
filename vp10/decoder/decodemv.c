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

#include "vp10/common/common.h"
#include "vp10/common/entropy.h"
#include "vp10/common/entropymode.h"
#include "vp10/common/entropymv.h"
#include "vp10/common/mvref_common.h"
#include "vp10/common/pred_common.h"
#include "vp10/common/reconinter.h"
#include "vp10/common/seg_common.h"

#include "vp10/decoder/decodemv.h"
#include "vp10/decoder/decodeframe.h"

#include "vpx_dsp/vpx_dsp_common.h"

static INLINE int read_uniform(vp10_reader *r, int n) {
  int l = get_unsigned_bits(n);
  int m = (1 << l) - n;
  int v = vp10_read_literal(r, l-1);

  assert(l != 0);

  if (v < m)
    return v;
  else
    return (v << 1) - m + vp10_read_literal(r, 1);
}

static PREDICTION_MODE read_intra_mode(vp10_reader *r, const vpx_prob *p) {
  return (PREDICTION_MODE)vp10_read_tree(r, vp10_intra_mode_tree, p);
}

static PREDICTION_MODE read_intra_mode_y(VP10_COMMON *cm, MACROBLOCKD *xd,
                                         vp10_reader *r, int size_group) {
  const PREDICTION_MODE y_mode =
      read_intra_mode(r, cm->fc->y_mode_prob[size_group]);
  FRAME_COUNTS *counts = xd->counts;
  if (counts)
    ++counts->y_mode[size_group][y_mode];
  return y_mode;
}

static PREDICTION_MODE read_intra_mode_uv(VP10_COMMON *cm, MACROBLOCKD *xd,
                                          vp10_reader *r,
                                          PREDICTION_MODE y_mode) {
  const PREDICTION_MODE uv_mode = read_intra_mode(r,
                                         cm->fc->uv_mode_prob[y_mode]);
  FRAME_COUNTS *counts = xd->counts;
  if (counts)
    ++counts->uv_mode[y_mode][uv_mode];
  return uv_mode;
}

#if CONFIG_EXT_INTER
static INTERINTRA_MODE read_interintra_mode(VP10_COMMON *cm, MACROBLOCKD *xd,
                                            vp10_reader *r, int size_group) {
  const INTERINTRA_MODE ii_mode =
      (INTERINTRA_MODE)vp10_read_tree(r, vp10_interintra_mode_tree,
                                      cm->fc->interintra_mode_prob[size_group]);
  FRAME_COUNTS *counts = xd->counts;
  if (counts)
    ++counts->interintra_mode[size_group][ii_mode];
  return ii_mode;
}
#endif  // CONFIG_EXT_INTER

static PREDICTION_MODE read_inter_mode(VP10_COMMON *cm, MACROBLOCKD *xd,
#if CONFIG_REF_MV && CONFIG_EXT_INTER
                                       MB_MODE_INFO *mbmi,
#endif
                                       vp10_reader *r, int16_t ctx) {
#if CONFIG_REF_MV
  FRAME_COUNTS *counts = xd->counts;
  int16_t mode_ctx = ctx & NEWMV_CTX_MASK;
  vpx_prob mode_prob = cm->fc->newmv_prob[mode_ctx];

  if (vp10_read(r, mode_prob) == 0) {
    if (counts)
      ++counts->newmv_mode[mode_ctx][0];

#if CONFIG_EXT_INTER
    if (has_second_ref(mbmi)) {
#endif  // CONFIG_EXT_INTER
    return NEWMV;
#if CONFIG_EXT_INTER
    } else {
      mode_prob = cm->fc->new2mv_prob;
      if (vp10_read(r, mode_prob) == 0) {
        if (counts)
          ++counts->new2mv_mode[0];
        return NEWMV;
      } else {
        if (counts)
          ++counts->new2mv_mode[1];
        return NEWFROMNEARMV;
      }
    }
#endif  // CONFIG_EXT_INTER
  }
  if (counts)
    ++counts->newmv_mode[mode_ctx][1];

  if (ctx & (1 << ALL_ZERO_FLAG_OFFSET))
    return ZEROMV;

  mode_ctx = (ctx >> ZEROMV_OFFSET) & ZEROMV_CTX_MASK;

  mode_prob = cm->fc->zeromv_prob[mode_ctx];
  if (vp10_read(r, mode_prob) == 0) {
    if (counts)
      ++counts->zeromv_mode[mode_ctx][0];
    return ZEROMV;
  }
  if (counts)
    ++counts->zeromv_mode[mode_ctx][1];

  mode_ctx = (ctx >> REFMV_OFFSET) & REFMV_CTX_MASK;

  if (ctx & (1 << SKIP_NEARESTMV_OFFSET))
    mode_ctx = 6;
  if (ctx & (1 << SKIP_NEARMV_OFFSET))
    mode_ctx = 7;
  if (ctx & (1 << SKIP_NEARESTMV_SUB8X8_OFFSET))
    mode_ctx = 8;

  mode_prob = cm->fc->refmv_prob[mode_ctx];

  if (vp10_read(r, mode_prob) == 0) {
    if (counts)
      ++counts->refmv_mode[mode_ctx][0];

    return NEARESTMV;
  } else {
    if (counts)
      ++counts->refmv_mode[mode_ctx][1];
    return NEARMV;
  }

  // Invalid prediction mode.
  assert(0);
#else
  const int mode = vp10_read_tree(r, vp10_inter_mode_tree,
                                 cm->fc->inter_mode_probs[ctx]);
  FRAME_COUNTS *counts = xd->counts;
  if (counts)
    ++counts->inter_mode[ctx][mode];

  return NEARESTMV + mode;
#endif
}

#if CONFIG_REF_MV
static void read_drl_idx(const VP10_COMMON *cm,
                         MACROBLOCKD *xd,
                         MB_MODE_INFO *mbmi,
                         vp10_reader *r) {
  uint8_t ref_frame_type = vp10_ref_frame_type(mbmi->ref_frame);
  mbmi->ref_mv_idx = 0;

  if (mbmi->mode == NEWMV) {
    int idx;
    for (idx = 0; idx < 2; ++idx) {
      if (xd->ref_mv_count[ref_frame_type] > idx + 1) {
        uint8_t drl_ctx = vp10_drl_ctx(xd->ref_mv_stack[ref_frame_type], idx);
        vpx_prob drl_prob = cm->fc->drl_prob[drl_ctx];
        if (!vp10_read(r, drl_prob)) {
          mbmi->ref_mv_idx = idx;
          if (xd->counts)
            ++xd->counts->drl_mode[drl_ctx][0];
          return;
        }
        mbmi->ref_mv_idx = idx + 1;
        if (xd->counts)
          ++xd->counts->drl_mode[drl_ctx][1];
      }
    }
  }

  if (mbmi->mode == NEARMV) {
    int idx;
    // Offset the NEARESTMV mode.
    // TODO(jingning): Unify the two syntax decoding loops after the NEARESTMV
    // mode is factored in.
    for (idx = 1; idx < 3; ++idx) {
      if (xd->ref_mv_count[ref_frame_type] > idx + 1) {
        uint8_t drl_ctx = vp10_drl_ctx(xd->ref_mv_stack[ref_frame_type], idx);
        vpx_prob drl_prob = cm->fc->drl_prob[drl_ctx];
        if (!vp10_read(r, drl_prob)) {
          mbmi->ref_mv_idx = idx - 1;
          if (xd->counts)
            ++xd->counts->drl_mode[drl_ctx][0];
          return;
        }
        mbmi->ref_mv_idx = idx;
        if (xd->counts)
          ++xd->counts->drl_mode[drl_ctx][1];
      }
    }
  }
}
#endif

#if CONFIG_EXT_INTER
static PREDICTION_MODE read_inter_compound_mode(VP10_COMMON *cm,
                                                MACROBLOCKD *xd,
                                                vp10_reader *r, int16_t ctx) {
  const int mode = vp10_read_tree(r, vp10_inter_compound_mode_tree,
                                 cm->fc->inter_compound_mode_probs[ctx]);
  FRAME_COUNTS *counts = xd->counts;

  if (counts)
    ++counts->inter_compound_mode[ctx][mode];

  assert(is_inter_compound_mode(NEAREST_NEARESTMV + mode));
  return NEAREST_NEARESTMV + mode;
}
#endif  // CONFIG_EXT_INTER

static int read_segment_id(vp10_reader *r,
    const struct segmentation_probs *segp) {
  return vp10_read_tree(r, vp10_segment_tree, segp->tree_probs);
}

#if CONFIG_VAR_TX
static void read_tx_size_vartx(VP10_COMMON *cm, MACROBLOCKD *xd,
                               MB_MODE_INFO *mbmi, FRAME_COUNTS *counts,
                               TX_SIZE tx_size, int blk_row, int blk_col,
                               vp10_reader *r) {
  int is_split = 0;
  const int tx_row = blk_row >> 1;
  const int tx_col = blk_col >> 1;
  int max_blocks_high = num_4x4_blocks_high_lookup[mbmi->sb_type];
  int max_blocks_wide = num_4x4_blocks_wide_lookup[mbmi->sb_type];
  int ctx = txfm_partition_context(xd->above_txfm_context + tx_col,
                                   xd->left_txfm_context + tx_row,
                                   tx_size);
  TX_SIZE (*const inter_tx_size)[MAX_MIB_SIZE] =
    (TX_SIZE (*)[MAX_MIB_SIZE])&mbmi->inter_tx_size[tx_row][tx_col];

  if (xd->mb_to_bottom_edge < 0)
    max_blocks_high += xd->mb_to_bottom_edge >> 5;
  if (xd->mb_to_right_edge < 0)
     max_blocks_wide += xd->mb_to_right_edge >> 5;

  if (blk_row >= max_blocks_high || blk_col >= max_blocks_wide)
     return;

  is_split = vp10_read(r, cm->fc->txfm_partition_prob[ctx]);

  if (is_split) {
    BLOCK_SIZE bsize = txsize_to_bsize[tx_size];
    int bsl = b_width_log2_lookup[bsize];
    int i;

    if (counts)
      ++counts->txfm_partition[ctx][1];

    if (tx_size == TX_8X8) {
      inter_tx_size[0][0] = TX_4X4;
      mbmi->tx_size = TX_4X4;
      txfm_partition_update(xd->above_txfm_context + tx_col,
                            xd->left_txfm_context + tx_row, TX_4X4);
      return;
    }

    assert(bsl > 0);
    --bsl;
    for (i = 0; i < 4; ++i) {
      int offsetr = blk_row + ((i >> 1) << bsl);
      int offsetc = blk_col + ((i & 0x01) << bsl);
      read_tx_size_vartx(cm, xd, mbmi, counts,
                         tx_size - 1, offsetr, offsetc, r);
    }
  } else {
    int idx, idy;
    inter_tx_size[0][0] = tx_size;
    for (idy = 0; idy < num_4x4_blocks_high_txsize_lookup[tx_size] / 2; ++idy)
      for (idx = 0; idx < num_4x4_blocks_wide_txsize_lookup[tx_size] / 2; ++idx)
        inter_tx_size[idy][idx] = tx_size;
    mbmi->tx_size = tx_size;
    if (counts)
      ++counts->txfm_partition[ctx][0];
    txfm_partition_update(xd->above_txfm_context + tx_col,
                          xd->left_txfm_context + tx_row, tx_size);
  }
}
#endif

static TX_SIZE read_selected_tx_size(VP10_COMMON *cm, MACROBLOCKD *xd,
                                     TX_SIZE max_tx_size, vp10_reader *r) {
  FRAME_COUNTS *counts = xd->counts;
  const int ctx = get_tx_size_context(xd);
  const int tx_size_cat = max_tx_size - TX_8X8;
  int tx_size = vp10_read_tree(r, vp10_tx_size_tree[tx_size_cat],
                              cm->fc->tx_size_probs[tx_size_cat][ctx]);
  if (counts)
    ++counts->tx_size[tx_size_cat][ctx][tx_size];
  return (TX_SIZE)tx_size;
}

static TX_SIZE read_tx_size_intra(VP10_COMMON *cm, MACROBLOCKD *xd,
                                  vp10_reader *r) {
  TX_MODE tx_mode = cm->tx_mode;
  BLOCK_SIZE bsize = xd->mi[0]->mbmi.sb_type;
  if (xd->lossless[xd->mi[0]->mbmi.segment_id])
    return TX_4X4;
  if (bsize >= BLOCK_8X8) {
    const TX_SIZE max_tx_size = max_txsize_lookup[bsize];
    if (tx_mode == TX_MODE_SELECT) {
      return read_selected_tx_size(cm, xd, max_tx_size, r);
    } else {
      return VPXMIN(max_tx_size, tx_mode_to_biggest_tx_size[tx_mode]);
    }
  } else {
    return TX_4X4;
  }
}

static TX_SIZE read_tx_size_inter(VP10_COMMON *cm, MACROBLOCKD *xd,
                                  int allow_select, vp10_reader *r) {
  TX_MODE tx_mode = cm->tx_mode;
  BLOCK_SIZE bsize = xd->mi[0]->mbmi.sb_type;
  if (xd->lossless[xd->mi[0]->mbmi.segment_id])
    return TX_4X4;
  if (bsize >= BLOCK_8X8) {
    const TX_SIZE max_tx_size = max_txsize_lookup[bsize];
    if (allow_select && tx_mode == TX_MODE_SELECT) {
      return read_selected_tx_size(cm, xd, max_tx_size, r);
    } else {
      return VPXMIN(max_tx_size, tx_mode_to_biggest_tx_size[tx_mode]);
    }
  } else {
#if CONFIG_EXT_TX && CONFIG_RECT_TX && !CONFIG_VAR_TX
    return max_txsize_rect_lookup[bsize];
#else
    return TX_4X4;
#endif  // CONFIG_EXT_TX && CONFIG_RECT_TX && !CONFIG_VAR_TX
  }
}

static int dec_get_segment_id(const VP10_COMMON *cm, const uint8_t *segment_ids,
                              int mi_offset, int x_mis, int y_mis) {
  int x, y, segment_id = INT_MAX;

  for (y = 0; y < y_mis; y++)
    for (x = 0; x < x_mis; x++)
      segment_id =
          VPXMIN(segment_id, segment_ids[mi_offset + y * cm->mi_cols + x]);

  assert(segment_id >= 0 && segment_id < MAX_SEGMENTS);
  return segment_id;
}

static void set_segment_id(VP10_COMMON *cm, int mi_offset,
                           int x_mis, int y_mis, int segment_id) {
  int x, y;

  assert(segment_id >= 0 && segment_id < MAX_SEGMENTS);

  for (y = 0; y < y_mis; y++)
    for (x = 0; x < x_mis; x++)
      cm->current_frame_seg_map[mi_offset + y * cm->mi_cols + x] = segment_id;
}

static int read_intra_segment_id(VP10_COMMON *const cm, MACROBLOCKD *const xd,
                                 int mi_offset, int x_mis, int y_mis,
                                 vp10_reader *r) {
  struct segmentation *const seg = &cm->seg;
  FRAME_COUNTS *counts = xd->counts;
  struct segmentation_probs *const segp = &cm->fc->seg;
  int segment_id;

  if (!seg->enabled)
    return 0;  // Default for disabled segmentation

  assert(seg->update_map && !seg->temporal_update);

  segment_id = read_segment_id(r, segp);
  if (counts)
    ++counts->seg.tree_total[segment_id];
  set_segment_id(cm, mi_offset, x_mis, y_mis, segment_id);
  return segment_id;
}

static void copy_segment_id(const VP10_COMMON *cm,
                           const uint8_t *last_segment_ids,
                           uint8_t *current_segment_ids,
                           int mi_offset, int x_mis, int y_mis) {
  int x, y;

  for (y = 0; y < y_mis; y++)
    for (x = 0; x < x_mis; x++)
      current_segment_ids[mi_offset + y * cm->mi_cols + x] =  last_segment_ids ?
          last_segment_ids[mi_offset + y * cm->mi_cols + x] : 0;
}

static int read_inter_segment_id(VP10_COMMON *const cm, MACROBLOCKD *const xd,
                                 int mi_row, int mi_col, vp10_reader *r) {
  struct segmentation *const seg = &cm->seg;
  FRAME_COUNTS *counts = xd->counts;
  struct segmentation_probs *const segp = &cm->fc->seg;
  MB_MODE_INFO *const mbmi = &xd->mi[0]->mbmi;
  int predicted_segment_id, segment_id;
  const int mi_offset = mi_row * cm->mi_cols + mi_col;
  const int bw = num_8x8_blocks_wide_lookup[mbmi->sb_type];
  const int bh = num_8x8_blocks_high_lookup[mbmi->sb_type];

  // TODO(slavarnway): move x_mis, y_mis into xd ?????
  const int x_mis = VPXMIN(cm->mi_cols - mi_col, bw);
  const int y_mis = VPXMIN(cm->mi_rows - mi_row, bh);

  if (!seg->enabled)
    return 0;  // Default for disabled segmentation

  predicted_segment_id = cm->last_frame_seg_map ?
      dec_get_segment_id(cm, cm->last_frame_seg_map, mi_offset, x_mis, y_mis) :
      0;

  if (!seg->update_map) {
    copy_segment_id(cm, cm->last_frame_seg_map, cm->current_frame_seg_map,
                    mi_offset, x_mis, y_mis);
    return predicted_segment_id;
  }

  if (seg->temporal_update) {
    const int ctx = vp10_get_pred_context_seg_id(xd);
    const vpx_prob pred_prob = segp->pred_probs[ctx];
    mbmi->seg_id_predicted = vp10_read(r, pred_prob);
    if (counts)
      ++counts->seg.pred[ctx][mbmi->seg_id_predicted];
    if (mbmi->seg_id_predicted) {
      segment_id = predicted_segment_id;
    } else {
      segment_id = read_segment_id(r, segp);
      if (counts)
        ++counts->seg.tree_mispred[segment_id];
    }
  } else {
    segment_id = read_segment_id(r, segp);
    if (counts)
      ++counts->seg.tree_total[segment_id];
  }
  set_segment_id(cm, mi_offset, x_mis, y_mis, segment_id);
  return segment_id;
}

static int read_skip(VP10_COMMON *cm, const MACROBLOCKD *xd,
                     int segment_id, vp10_reader *r) {
  if (segfeature_active(&cm->seg, segment_id, SEG_LVL_SKIP)) {
    return 1;
  } else {
    const int ctx = vp10_get_skip_context(xd);
    const int skip = vp10_read(r, cm->fc->skip_probs[ctx]);
    FRAME_COUNTS *counts = xd->counts;
    if (counts)
      ++counts->skip[ctx][skip];
    return skip;
  }
}

static void read_palette_mode_info(VP10_COMMON *const cm,
                                   MACROBLOCKD *const xd,
                                   vp10_reader *r) {
  MODE_INFO *const mi = xd->mi[0];
  MB_MODE_INFO *const mbmi = &mi->mbmi;
  const MODE_INFO *const above_mi = xd->above_mi;
  const MODE_INFO *const left_mi  = xd->left_mi;
  const BLOCK_SIZE bsize = mbmi->sb_type;
  int i, n, palette_ctx = 0;
  PALETTE_MODE_INFO *const pmi = &mbmi->palette_mode_info;

  if (mbmi->mode == DC_PRED) {
    if (above_mi)
      palette_ctx += (above_mi->mbmi.palette_mode_info.palette_size[0] > 0);
    if (left_mi)
      palette_ctx += (left_mi->mbmi.palette_mode_info.palette_size[0] > 0);
    if (vp10_read(r, vp10_default_palette_y_mode_prob[bsize - BLOCK_8X8]
                                                     [palette_ctx])) {
      pmi->palette_size[0] =
        vp10_read_tree(r, vp10_palette_size_tree,
                      vp10_default_palette_y_size_prob[bsize - BLOCK_8X8]) + 2;
      n = pmi->palette_size[0];
      for (i = 0; i < n; ++i)
        pmi->palette_colors[i] = vp10_read_literal(r, cm->bit_depth);

      xd->plane[0].color_index_map[0] = read_uniform(r, n);
      assert(xd->plane[0].color_index_map[0] < n);
    }
  }

  if (mbmi->uv_mode == DC_PRED) {
    if (vp10_read(r,
                 vp10_default_palette_uv_mode_prob[pmi->palette_size[0] > 0])) {
      pmi->palette_size[1] =
          vp10_read_tree(r, vp10_palette_size_tree,
                        vp10_default_palette_uv_size_prob[bsize - BLOCK_8X8])
                        + 2;
      n = pmi->palette_size[1];
      for (i = 0; i < n; ++i) {
        pmi->palette_colors[PALETTE_MAX_SIZE + i] =
            vp10_read_literal(r, cm->bit_depth);
        pmi->palette_colors[2 * PALETTE_MAX_SIZE + i] =
            vp10_read_literal(r, cm->bit_depth);
      }
      xd->plane[1].color_index_map[0] = read_uniform(r, n);
      assert(xd->plane[1].color_index_map[0] < n);
    }
  }
}

#if CONFIG_EXT_INTRA
static void read_ext_intra_mode_info(VP10_COMMON *const cm,
                                     MACROBLOCKD *const xd, vp10_reader *r) {
  MODE_INFO *const mi = xd->mi[0];
  MB_MODE_INFO *const mbmi = &mi->mbmi;
  FRAME_COUNTS *counts = xd->counts;

#if !ALLOW_FILTER_INTRA_MODES
  return;
#endif
  if (mbmi->mode == DC_PRED &&
      mbmi->palette_mode_info.palette_size[0] == 0) {
    mbmi->ext_intra_mode_info.use_ext_intra_mode[0] =
        vp10_read(r, cm->fc->ext_intra_probs[0]);
    if (mbmi->ext_intra_mode_info.use_ext_intra_mode[0]) {
      mbmi->ext_intra_mode_info.ext_intra_mode[0] =
          read_uniform(r, FILTER_INTRA_MODES);
    }
    if (counts)
      ++counts->ext_intra[0][mbmi->ext_intra_mode_info.use_ext_intra_mode[0]];
  }
  if (mbmi->uv_mode == DC_PRED &&
      mbmi->palette_mode_info.palette_size[1] == 0) {
    mbmi->ext_intra_mode_info.use_ext_intra_mode[1] =
        vp10_read(r, cm->fc->ext_intra_probs[1]);
    if (mbmi->ext_intra_mode_info.use_ext_intra_mode[1]) {
      mbmi->ext_intra_mode_info.ext_intra_mode[1] =
          read_uniform(r, FILTER_INTRA_MODES);
    }
    if (counts)
      ++counts->ext_intra[1][mbmi->ext_intra_mode_info.use_ext_intra_mode[1]];
  }
}

static void read_intra_angle_info(VP10_COMMON *const cm, MACROBLOCKD *const xd,
                                  vp10_reader *r) {
  MB_MODE_INFO *const mbmi = &xd->mi[0]->mbmi;
  const BLOCK_SIZE bsize = mbmi->sb_type;
  const int ctx = vp10_get_pred_context_intra_interp(xd);
  int p_angle;

  if (bsize < BLOCK_8X8)
    return;

  if (mbmi->mode != DC_PRED && mbmi->mode != TM_PRED) {
    mbmi->angle_delta[0] =
        read_uniform(r, 2 * MAX_ANGLE_DELTAS + 1) - MAX_ANGLE_DELTAS;
    p_angle = mode_to_angle_map[mbmi->mode] + mbmi->angle_delta[0] * ANGLE_STEP;
    if (vp10_is_intra_filter_switchable(p_angle)) {
      FRAME_COUNTS *counts = xd->counts;
      mbmi->intra_filter = vp10_read_tree(r, vp10_intra_filter_tree,
                                          cm->fc->intra_filter_probs[ctx]);
      if (counts)
        ++counts->intra_filter[ctx][mbmi->intra_filter];
    } else {
      mbmi->intra_filter = INTRA_FILTER_LINEAR;
    }
  }

  if (mbmi->uv_mode != DC_PRED && mbmi->uv_mode != TM_PRED) {
    mbmi->angle_delta[1] =
        read_uniform(r, 2 * MAX_ANGLE_DELTAS + 1) - MAX_ANGLE_DELTAS;
  }
}
#endif  // CONFIG_EXT_INTRA

static void read_intra_frame_mode_info(VP10_COMMON *const cm,
                                       MACROBLOCKD *const xd,
                                       int mi_row, int mi_col, vp10_reader *r) {
  MODE_INFO *const mi = xd->mi[0];
  MB_MODE_INFO *const mbmi = &mi->mbmi;
  const MODE_INFO *above_mi = xd->above_mi;
  const MODE_INFO *left_mi  = xd->left_mi;
  const BLOCK_SIZE bsize = mbmi->sb_type;
  int i;
  const int mi_offset = mi_row * cm->mi_cols + mi_col;
  const int bw = xd->plane[0].n4_w >> 1;
  const int bh = xd->plane[0].n4_h >> 1;

  // TODO(slavarnway): move x_mis, y_mis into xd ?????
  const int x_mis = VPXMIN(cm->mi_cols - mi_col, bw);
  const int y_mis = VPXMIN(cm->mi_rows - mi_row, bh);

  mbmi->segment_id = read_intra_segment_id(cm, xd, mi_offset, x_mis, y_mis, r);
  mbmi->skip = read_skip(cm, xd, mbmi->segment_id, r);
  mbmi->tx_size = read_tx_size_intra(cm, xd, r);
  mbmi->ref_frame[0] = INTRA_FRAME;
  mbmi->ref_frame[1] = NONE;

  switch (bsize) {
    case BLOCK_4X4:
      for (i = 0; i < 4; ++i)
        mi->bmi[i].as_mode =
            read_intra_mode(r, get_y_mode_probs(cm, mi, above_mi, left_mi, i));
      mbmi->mode = mi->bmi[3].as_mode;
      break;
    case BLOCK_4X8:
      mi->bmi[0].as_mode = mi->bmi[2].as_mode =
          read_intra_mode(r, get_y_mode_probs(cm, mi, above_mi, left_mi, 0));
      mi->bmi[1].as_mode = mi->bmi[3].as_mode = mbmi->mode =
          read_intra_mode(r, get_y_mode_probs(cm, mi, above_mi, left_mi, 1));
      break;
    case BLOCK_8X4:
      mi->bmi[0].as_mode = mi->bmi[1].as_mode =
          read_intra_mode(r, get_y_mode_probs(cm, mi, above_mi, left_mi, 0));
      mi->bmi[2].as_mode = mi->bmi[3].as_mode = mbmi->mode =
          read_intra_mode(r, get_y_mode_probs(cm, mi, above_mi, left_mi, 2));
      break;
    default:
      mbmi->mode = read_intra_mode(r,
          get_y_mode_probs(cm, mi, above_mi, left_mi, 0));
  }

  mbmi->uv_mode = read_intra_mode_uv(cm, xd, r, mbmi->mode);
#if CONFIG_EXT_INTRA
  read_intra_angle_info(cm, xd, r);
#endif  // CONFIG_EXT_INTRA
  mbmi->palette_mode_info.palette_size[0] = 0;
  mbmi->palette_mode_info.palette_size[1] = 0;
  if (bsize >= BLOCK_8X8 && cm->allow_screen_content_tools)
    read_palette_mode_info(cm, xd, r);
#if CONFIG_EXT_INTRA
    mbmi->ext_intra_mode_info.use_ext_intra_mode[0] = 0;
    mbmi->ext_intra_mode_info.use_ext_intra_mode[1] = 0;
    if (bsize >= BLOCK_8X8)
      read_ext_intra_mode_info(cm, xd, r);
#endif  // CONFIG_EXT_INTRA

  if (!FIXED_TX_TYPE) {
#if CONFIG_EXT_TX
    if (get_ext_tx_types(mbmi->tx_size, mbmi->sb_type, 0) > 1 &&
        cm->base_qindex > 0 && !mbmi->skip &&
        !segfeature_active(&cm->seg, mbmi->segment_id, SEG_LVL_SKIP) &&
        ALLOW_INTRA_EXT_TX) {
      FRAME_COUNTS *counts = xd->counts;
      int eset = get_ext_tx_set(mbmi->tx_size, mbmi->sb_type, 0);
      if (eset > 0) {
        mbmi->tx_type = vp10_read_tree(
            r, vp10_ext_tx_intra_tree[eset],
            cm->fc->intra_ext_tx_prob[eset][mbmi->tx_size][mbmi->mode]);
        if (counts)
          ++counts->intra_ext_tx[eset][mbmi->tx_size][mbmi->mode]
                                                     [mbmi->tx_type];
      }
    } else {
      mbmi->tx_type = DCT_DCT;
    }
#else
    if (mbmi->tx_size < TX_32X32 &&
        cm->base_qindex > 0 && !mbmi->skip &&
        !segfeature_active(&cm->seg, mbmi->segment_id, SEG_LVL_SKIP)) {
      FRAME_COUNTS *counts = xd->counts;
      TX_TYPE tx_type_nom = intra_mode_to_tx_type_context[mbmi->mode];
      mbmi->tx_type = vp10_read_tree(
          r, vp10_ext_tx_tree,
          cm->fc->intra_ext_tx_prob[mbmi->tx_size][tx_type_nom]);
      if (counts)
        ++counts->intra_ext_tx[mbmi->tx_size][tx_type_nom][mbmi->tx_type];
    } else {
      mbmi->tx_type = DCT_DCT;
    }
#endif  // CONFIG_EXT_TX
  }
}

static int read_mv_component(vp10_reader *r,
                             const nmv_component *mvcomp, int usehp) {
  int mag, d, fr, hp;
  const int sign = vp10_read(r, mvcomp->sign);
  const int mv_class = vp10_read_tree(r, vp10_mv_class_tree, mvcomp->classes);
  const int class0 = mv_class == MV_CLASS_0;

  // Integer part
  if (class0) {
    d = vp10_read_tree(r, vp10_mv_class0_tree, mvcomp->class0);
    mag = 0;
  } else {
    int i;
    const int n = mv_class + CLASS0_BITS - 1;  // number of bits

    d = 0;
    for (i = 0; i < n; ++i)
      d |= vp10_read(r, mvcomp->bits[i]) << i;
    mag = CLASS0_SIZE << (mv_class + 2);
  }

  // Fractional part
  fr = vp10_read_tree(r, vp10_mv_fp_tree, class0 ? mvcomp->class0_fp[d]
                                               : mvcomp->fp);

  // High precision part (if hp is not used, the default value of the hp is 1)
  hp = usehp ? vp10_read(r, class0 ? mvcomp->class0_hp : mvcomp->hp)
             : 1;

  // Result
  mag += ((d << 3) | (fr << 1) | hp) + 1;
  return sign ? -mag : mag;
}

static INLINE void read_mv(vp10_reader *r, MV *mv, const MV *ref,
#if CONFIG_REF_MV
                           int is_compound,
#endif
                           const nmv_context *ctx,
                           nmv_context_counts *counts, int allow_hp) {
  MV_JOINT_TYPE joint_type;
  const int use_hp = allow_hp && vp10_use_mv_hp(ref);
  MV diff = {0, 0};

#if CONFIG_REF_MV && !CONFIG_EXT_INTER
  if (is_compound) {
    int is_zero_rmv = vp10_read(r, ctx->zero_rmv);
    if (is_zero_rmv) {
      joint_type = MV_JOINT_ZERO;
    } else {
      joint_type = (MV_JOINT_TYPE)vp10_read_tree(r, vp10_mv_joint_tree,
                                                 ctx->joints);
    }
  } else {
    joint_type = (MV_JOINT_TYPE)vp10_read_tree(r, vp10_mv_joint_tree,
                                               ctx->joints);
  }
#else
  joint_type = (MV_JOINT_TYPE)vp10_read_tree(r, vp10_mv_joint_tree,
                                             ctx->joints);
#endif

#if CONFIG_REF_MV && CONFIG_EXT_INTER
  (void)is_compound;
#endif

  if (mv_joint_vertical(joint_type))
    diff.row = read_mv_component(r, &ctx->comps[0], use_hp);

  if (mv_joint_horizontal(joint_type))
    diff.col = read_mv_component(r, &ctx->comps[1], use_hp);

  vp10_inc_mv(&diff, counts, use_hp);

  mv->row = ref->row + diff.row;
  mv->col = ref->col + diff.col;
}

static REFERENCE_MODE read_block_reference_mode(VP10_COMMON *cm,
                                                const MACROBLOCKD *xd,
                                                vp10_reader *r) {
  if (cm->reference_mode == REFERENCE_MODE_SELECT) {
    const int ctx = vp10_get_reference_mode_context(cm, xd);
    const REFERENCE_MODE mode =
        (REFERENCE_MODE)vp10_read(r, cm->fc->comp_inter_prob[ctx]);
    FRAME_COUNTS *counts = xd->counts;
    if (counts)
      ++counts->comp_inter[ctx][mode];
    return mode;  // SINGLE_REFERENCE or COMPOUND_REFERENCE
  } else {
    return cm->reference_mode;
  }
}

// Read the referncence frame
static void read_ref_frames(VP10_COMMON *const cm, MACROBLOCKD *const xd,
                            vp10_reader *r,
                            int segment_id, MV_REFERENCE_FRAME ref_frame[2]) {
  FRAME_CONTEXT *const fc = cm->fc;
  FRAME_COUNTS *counts = xd->counts;

  if (segfeature_active(&cm->seg, segment_id, SEG_LVL_REF_FRAME)) {
    ref_frame[0] = (MV_REFERENCE_FRAME)get_segdata(&cm->seg, segment_id,
                                                   SEG_LVL_REF_FRAME);
    ref_frame[1] = NONE;
  } else {
    const REFERENCE_MODE mode = read_block_reference_mode(cm, xd, r);
    // FIXME(rbultje) I'm pretty sure this breaks segmentation ref frame coding
    if (mode == COMPOUND_REFERENCE) {
#if CONFIG_EXT_REFS
      const int idx = cm->ref_frame_sign_bias[cm->comp_bwd_ref[0]];
#else
      const int idx = cm->ref_frame_sign_bias[cm->comp_fixed_ref];
#endif  // CONFIG_EXT_REFS
      const int ctx = vp10_get_pred_context_comp_ref_p(cm, xd);
      const int bit = vp10_read(r, fc->comp_ref_prob[ctx][0]);

      if (counts)
        ++counts->comp_ref[ctx][0][bit];

#if CONFIG_EXT_REFS
      // Decode forward references.
      if (!bit) {
        const int ctx1 = vp10_get_pred_context_comp_ref_p1(cm, xd);
        const int bit1 = vp10_read(r, fc->comp_ref_prob[ctx1][1]);
        if (counts)
          ++counts->comp_ref[ctx1][1][bit1];
        ref_frame[!idx] = cm->comp_fwd_ref[bit1 ? 0 : 1];
      } else {
        const int ctx2 = vp10_get_pred_context_comp_ref_p2(cm, xd);
        const int bit2 = vp10_read(r, fc->comp_ref_prob[ctx2][2]);
        if (counts)
          ++counts->comp_ref[ctx2][2][bit2];
        ref_frame[!idx] = cm->comp_fwd_ref[bit2 ? 3 : 2];
      }

      // Decode backward references.
      {
        const int ctx_bwd = vp10_get_pred_context_comp_bwdref_p(cm, xd);
        const int bit_bwd = vp10_read(r, fc->comp_bwdref_prob[ctx_bwd][0]);
        if (counts)
          ++counts->comp_bwdref[ctx_bwd][0][bit_bwd];
        ref_frame[idx] = cm->comp_bwd_ref[bit_bwd];
      }
#else
      ref_frame[!idx] = cm->comp_var_ref[bit];
      ref_frame[idx] = cm->comp_fixed_ref;
#endif  // CONFIG_EXT_REFS
    } else if (mode == SINGLE_REFERENCE) {
#if CONFIG_EXT_REFS
      const int ctx0 = vp10_get_pred_context_single_ref_p1(xd);
      const int bit0 = vp10_read(r, fc->single_ref_prob[ctx0][0]);
      if (counts)
        ++counts->single_ref[ctx0][0][bit0];

      if (bit0) {
        const int ctx1 = vp10_get_pred_context_single_ref_p2(xd);
        const int bit1 = vp10_read(r, fc->single_ref_prob[ctx1][1]);
        if (counts)
          ++counts->single_ref[ctx1][1][bit1];
        ref_frame[0] = bit1 ? ALTREF_FRAME : BWDREF_FRAME;
      } else {
        const int ctx2 = vp10_get_pred_context_single_ref_p3(xd);
        const int bit2 = vp10_read(r, fc->single_ref_prob[ctx2][2]);
        if (counts)
          ++counts->single_ref[ctx2][2][bit2];
        if (bit2) {
          const int ctx4 = vp10_get_pred_context_single_ref_p5(xd);
          const int bit4 = vp10_read(r, fc->single_ref_prob[ctx4][4]);
          if (counts)
            ++counts->single_ref[ctx4][4][bit4];
          ref_frame[0] = bit4 ? GOLDEN_FRAME : LAST3_FRAME;
        } else {
          const int ctx3 = vp10_get_pred_context_single_ref_p4(xd);
          const int bit3 = vp10_read(r, fc->single_ref_prob[ctx3][3]);
          if (counts)
            ++counts->single_ref[ctx3][3][bit3];
          ref_frame[0] = bit3 ? LAST2_FRAME : LAST_FRAME;
        }
      }
#else
      const int ctx0 = vp10_get_pred_context_single_ref_p1(xd);
      const int bit0 = vp10_read(r, fc->single_ref_prob[ctx0][0]);
      if (counts)
        ++counts->single_ref[ctx0][0][bit0];

      if (bit0) {
        const int ctx1 = vp10_get_pred_context_single_ref_p2(xd);
        const int bit1 = vp10_read(r, fc->single_ref_prob[ctx1][1]);
        if (counts)
          ++counts->single_ref[ctx1][1][bit1];
        ref_frame[0] = bit1 ? ALTREF_FRAME : GOLDEN_FRAME;
      } else {
        ref_frame[0] = LAST_FRAME;
      }
#endif  // CONFIG_EXT_REFS

      ref_frame[1] = NONE;
    } else {
      assert(0 && "Invalid prediction mode.");
    }
  }
}


#if CONFIG_OBMC || CONFIG_WARPED_MOTION
static MOTION_VARIATION read_motvar_block(
    VP10_COMMON *const cm, MACROBLOCKD *const xd, vp10_reader *r) {
  BLOCK_SIZE bsize = xd->mi[0]->mbmi.sb_type;
  FRAME_COUNTS *counts = xd->counts;
  MOTION_VARIATION motvar;

  if (is_motvar_allowed(&xd->mi[0]->mbmi)) {
    motvar = (MOTION_VARIATION)
        vp10_read_tree(r, vp10_motvar_tree, cm->fc->motvar_prob[bsize]);
    if (counts)
      ++counts->motvar[bsize][motvar];
    return motvar;
  } else {
    return SIMPLE_TRANSLATION;
  }
}
#endif  // CONFIG_OBMC || CONFIG_WARPED_MOTION

static INLINE INTERP_FILTER read_interp_filter(
    VP10_COMMON *const cm, MACROBLOCKD *const xd,
#if CONFIG_DUAL_FILTER
    int dir,
#endif
    vp10_reader *r) {
#if CONFIG_EXT_INTERP
  if (!vp10_is_interp_needed(xd)) return EIGHTTAP_REGULAR;
#endif
  if (cm->interp_filter != SWITCHABLE) {
    return cm->interp_filter;
  } else {
#if CONFIG_DUAL_FILTER
    const int ctx = vp10_get_pred_context_switchable_interp(xd, dir);
#else
    const int ctx = vp10_get_pred_context_switchable_interp(xd);
#endif
    FRAME_COUNTS *counts = xd->counts;
    const INTERP_FILTER type =
      (INTERP_FILTER)vp10_read_tree(r, vp10_switchable_interp_tree,
                                    cm->fc->switchable_interp_prob[ctx]);
    if (counts)
      ++counts->switchable_interp[ctx][type];
    return type;
  }
}

static void read_intra_block_mode_info(VP10_COMMON *const cm,
                                       MACROBLOCKD *const xd, MODE_INFO *mi,
                                       vp10_reader *r) {
  MB_MODE_INFO *const mbmi = &mi->mbmi;
  const BLOCK_SIZE bsize = mi->mbmi.sb_type;
  int i;

  mbmi->ref_frame[0] = INTRA_FRAME;
  mbmi->ref_frame[1] = NONE;

  switch (bsize) {
    case BLOCK_4X4:
      for (i = 0; i < 4; ++i)
        mi->bmi[i].as_mode = read_intra_mode_y(cm, xd, r, 0);
      mbmi->mode = mi->bmi[3].as_mode;
      break;
    case BLOCK_4X8:
      mi->bmi[0].as_mode = mi->bmi[2].as_mode = read_intra_mode_y(cm, xd,
                                                                  r, 0);
      mi->bmi[1].as_mode = mi->bmi[3].as_mode = mbmi->mode =
          read_intra_mode_y(cm, xd, r, 0);
      break;
    case BLOCK_8X4:
      mi->bmi[0].as_mode = mi->bmi[1].as_mode = read_intra_mode_y(cm, xd,
                                                                  r, 0);
      mi->bmi[2].as_mode = mi->bmi[3].as_mode = mbmi->mode =
          read_intra_mode_y(cm, xd, r, 0);
      break;
    default:
      mbmi->mode = read_intra_mode_y(cm, xd, r, size_group_lookup[bsize]);
  }

  mbmi->uv_mode = read_intra_mode_uv(cm, xd, r, mbmi->mode);
#if CONFIG_EXT_INTRA
  read_intra_angle_info(cm, xd, r);
#endif  // CONFIG_EXT_INTRA
  mbmi->palette_mode_info.palette_size[0] = 0;
  mbmi->palette_mode_info.palette_size[1] = 0;
  if (bsize >= BLOCK_8X8 && cm->allow_screen_content_tools)
    read_palette_mode_info(cm, xd, r);
#if CONFIG_EXT_INTRA
  mbmi->ext_intra_mode_info.use_ext_intra_mode[0] = 0;
  mbmi->ext_intra_mode_info.use_ext_intra_mode[1] = 0;
  if (bsize >= BLOCK_8X8)
    read_ext_intra_mode_info(cm, xd, r);
#endif  // CONFIG_EXT_INTRA
}

static INLINE int is_mv_valid(const MV *mv) {
  return mv->row > MV_LOW && mv->row < MV_UPP &&
         mv->col > MV_LOW && mv->col < MV_UPP;
}

static INLINE int assign_mv(VP10_COMMON *cm, MACROBLOCKD *xd,
                            PREDICTION_MODE mode,
#if CONFIG_REF_MV
                            int block,
#endif
                            int_mv mv[2], int_mv ref_mv[2],
                            int_mv nearest_mv[2], int_mv near_mv[2],
                            int is_compound, int allow_hp, vp10_reader *r) {
  int i;
  int ret = 1;
#if CONFIG_REF_MV
  MB_MODE_INFO *mbmi = &xd->mi[0]->mbmi;
  BLOCK_SIZE bsize = mbmi->sb_type;
  int_mv *pred_mv = (bsize >= BLOCK_8X8) ?
      mbmi->pred_mv : xd->mi[0]->bmi[block].pred_mv_s8;
#endif

  switch (mode) {
#if CONFIG_EXT_INTER
    case NEWFROMNEARMV:
#endif  // CONFIG_EXT_INTER
    case NEWMV: {
      FRAME_COUNTS *counts = xd->counts;
#if !CONFIG_REF_MV
      nmv_context_counts *const mv_counts = counts ? &counts->mv : NULL;
#endif
      for (i = 0; i < 1 + is_compound; ++i) {
#if CONFIG_REF_MV
        int nmv_ctx = vp10_nmv_ctx(xd->ref_mv_count[mbmi->ref_frame[i]],
                                   xd->ref_mv_stack[mbmi->ref_frame[i]]);
        nmv_context_counts *const mv_counts =
            counts ? &counts->mv[nmv_ctx] : NULL;
        read_mv(r, &mv[i].as_mv, &ref_mv[i].as_mv,
#if CONFIG_REF_MV
                is_compound,
#endif
                &cm->fc->nmvc[nmv_ctx], mv_counts, allow_hp);
#else
        read_mv(r, &mv[i].as_mv, &ref_mv[i].as_mv, &cm->fc->nmvc, mv_counts,
                allow_hp);
#endif
        ret = ret && is_mv_valid(&mv[i].as_mv);

#if CONFIG_REF_MV
        pred_mv[i].as_int = ref_mv[i].as_int;
#endif
      }
      break;
    }
    case NEARESTMV: {
      mv[0].as_int = nearest_mv[0].as_int;
      if (is_compound)
        mv[1].as_int = nearest_mv[1].as_int;

#if CONFIG_REF_MV
      pred_mv[0].as_int = nearest_mv[0].as_int;
      if (is_compound)
        pred_mv[1].as_int = nearest_mv[1].as_int;
#endif
      break;
    }
    case NEARMV: {
      mv[0].as_int = near_mv[0].as_int;
      if (is_compound)
        mv[1].as_int = near_mv[1].as_int;

#if CONFIG_REF_MV
      pred_mv[0].as_int = near_mv[0].as_int;
      if (is_compound)
        pred_mv[1].as_int = near_mv[1].as_int;
#endif
      break;
    }
    case ZEROMV: {
      mv[0].as_int = 0;
      if (is_compound)
        mv[1].as_int = 0;

#if CONFIG_REF_MV
      pred_mv[0].as_int = 0;
      if (is_compound)
        pred_mv[1].as_int = 0;
#endif
      break;
    }
#if CONFIG_EXT_INTER
    case NEW_NEWMV: {
      FRAME_COUNTS *counts = xd->counts;
#if !CONFIG_REF_MV
      nmv_context_counts *const mv_counts = counts ? &counts->mv : NULL;
#endif
      assert(is_compound);
      for (i = 0; i < 2; ++i) {
#if CONFIG_REF_MV
        int nmv_ctx = vp10_nmv_ctx(xd->ref_mv_count[mbmi->ref_frame[i]],
                                   xd->ref_mv_stack[mbmi->ref_frame[i]]);
        nmv_context_counts *const mv_counts =
            counts ? &counts->mv[nmv_ctx] : NULL;
        read_mv(r, &mv[i].as_mv, &ref_mv[i].as_mv, is_compound,
                &cm->fc->nmvc[nmv_ctx], mv_counts,
                allow_hp);
#else
        read_mv(r, &mv[i].as_mv, &ref_mv[i].as_mv, &cm->fc->nmvc, mv_counts,
                allow_hp);
#endif
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
    case NEAR_NEARMV: {
      assert(is_compound);
      mv[0].as_int = near_mv[0].as_int;
      mv[1].as_int = near_mv[1].as_int;
      break;
    }
    case NEW_NEARESTMV: {
      FRAME_COUNTS *counts = xd->counts;
#if CONFIG_REF_MV
      int nmv_ctx = vp10_nmv_ctx(xd->ref_mv_count[mbmi->ref_frame[0]],
                                 xd->ref_mv_stack[mbmi->ref_frame[0]]);
      nmv_context_counts *const mv_counts =
          counts ? &counts->mv[nmv_ctx] : NULL;
      read_mv(r, &mv[0].as_mv, &ref_mv[0].as_mv, is_compound,
              &cm->fc->nmvc[nmv_ctx], mv_counts,
              allow_hp);
#else
      nmv_context_counts *const mv_counts = counts ? &counts->mv : NULL;
      read_mv(r, &mv[0].as_mv, &ref_mv[0].as_mv, &cm->fc->nmvc, mv_counts,
              allow_hp);
#endif
      assert(is_compound);
      ret = ret && is_mv_valid(&mv[0].as_mv);
      mv[1].as_int = nearest_mv[1].as_int;
      break;
    }
    case NEAREST_NEWMV: {
      FRAME_COUNTS *counts = xd->counts;
#if CONFIG_REF_MV
      int nmv_ctx = vp10_nmv_ctx(xd->ref_mv_count[mbmi->ref_frame[1]],
                                 xd->ref_mv_stack[mbmi->ref_frame[1]]);
      nmv_context_counts *const mv_counts =
          counts ? &counts->mv[nmv_ctx] : NULL;
      mv[0].as_int = nearest_mv[0].as_int;
      read_mv(r, &mv[1].as_mv, &ref_mv[1].as_mv, is_compound,
              &cm->fc->nmvc[nmv_ctx], mv_counts,
              allow_hp);
#else
      nmv_context_counts *const mv_counts = counts ? &counts->mv : NULL;
      mv[0].as_int = nearest_mv[0].as_int;
      read_mv(r, &mv[1].as_mv, &ref_mv[1].as_mv, &cm->fc->nmvc, mv_counts,
              allow_hp);
#endif
      assert(is_compound);
      ret = ret && is_mv_valid(&mv[1].as_mv);
      break;
    }
    case NEAR_NEWMV: {
      FRAME_COUNTS *counts = xd->counts;
#if CONFIG_REF_MV
      int nmv_ctx = vp10_nmv_ctx(xd->ref_mv_count[mbmi->ref_frame[1]],
                                 xd->ref_mv_stack[mbmi->ref_frame[1]]);
      nmv_context_counts *const mv_counts =
          counts ? &counts->mv[nmv_ctx] : NULL;
      mv[0].as_int = near_mv[0].as_int;
      read_mv(r, &mv[1].as_mv, &ref_mv[1].as_mv, is_compound,
              &cm->fc->nmvc[nmv_ctx], mv_counts,
              allow_hp);
#else
      nmv_context_counts *const mv_counts = counts ? &counts->mv : NULL;
      mv[0].as_int = near_mv[0].as_int;
      read_mv(r, &mv[1].as_mv, &ref_mv[1].as_mv, &cm->fc->nmvc, mv_counts,
              allow_hp);
#endif
      assert(is_compound);

      ret = ret && is_mv_valid(&mv[1].as_mv);
      break;
    }
    case NEW_NEARMV: {
      FRAME_COUNTS *counts = xd->counts;
#if CONFIG_REF_MV
      int nmv_ctx = vp10_nmv_ctx(xd->ref_mv_count[mbmi->ref_frame[0]],
                                 xd->ref_mv_stack[mbmi->ref_frame[0]]);
      nmv_context_counts *const mv_counts =
          counts ? &counts->mv[nmv_ctx] : NULL;
      read_mv(r, &mv[0].as_mv, &ref_mv[0].as_mv, is_compound,
              &cm->fc->nmvc[nmv_ctx], mv_counts,
              allow_hp);
#else
      nmv_context_counts *const mv_counts = counts ? &counts->mv : NULL;
      read_mv(r, &mv[0].as_mv, &ref_mv[0].as_mv, &cm->fc->nmvc, mv_counts,
              allow_hp);
#endif
      assert(is_compound);
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
#endif  // CONFIG_EXT_INTER
    default: {
      return 0;
    }
  }
  return ret;
}

static int read_is_inter_block(VP10_COMMON *const cm, MACROBLOCKD *const xd,
                               int segment_id, vp10_reader *r) {
  if (segfeature_active(&cm->seg, segment_id, SEG_LVL_REF_FRAME)) {
    return get_segdata(&cm->seg, segment_id, SEG_LVL_REF_FRAME) != INTRA_FRAME;
  } else {
    const int ctx = vp10_get_intra_inter_context(xd);
    const int is_inter = vp10_read(r, cm->fc->intra_inter_prob[ctx]);
    FRAME_COUNTS *counts = xd->counts;
    if (counts)
      ++counts->intra_inter[ctx][is_inter];
    return is_inter;
  }
}

static void fpm_sync(void *const data, int mi_row) {
  VP10Decoder *const pbi = (VP10Decoder *)data;
  vp10_frameworker_wait(pbi->frame_worker_owner, pbi->common.prev_frame,
                       mi_row << pbi->common.mib_size_log2);
}

static void read_inter_block_mode_info(VP10Decoder *const pbi,
                                       MACROBLOCKD *const xd,
                                       MODE_INFO *const mi,
#if (CONFIG_OBMC || CONFIG_EXT_INTER) && CONFIG_SUPERTX
                                       int mi_row, int mi_col, vp10_reader *r,
                                       int supertx_enabled) {
#else
                                       int mi_row, int mi_col, vp10_reader *r) {
#endif  // CONFIG_OBMC && CONFIG_SUPERTX
  VP10_COMMON *const cm = &pbi->common;
  MB_MODE_INFO *const mbmi = &mi->mbmi;
  const BLOCK_SIZE bsize = mbmi->sb_type;
  const int allow_hp = cm->allow_high_precision_mv;
  int_mv nearestmv[2], nearmv[2];
  int_mv ref_mvs[MODE_CTX_REF_FRAMES][MAX_MV_REF_CANDIDATES];
#if CONFIG_EXT_INTER
  int mv_idx;
#endif  // CONFIG_EXT_INTER
  int ref, is_compound;
  int16_t inter_mode_ctx[MODE_CTX_REF_FRAMES];
#if CONFIG_REF_MV && CONFIG_EXT_INTER
  int16_t compound_inter_mode_ctx[MODE_CTX_REF_FRAMES];
#endif  // CONFIG_REF_MV && CONFIG_EXT_INTER
  int16_t mode_ctx = 0;
  MV_REFERENCE_FRAME ref_frame;

  mbmi->palette_mode_info.palette_size[0] = 0;
  mbmi->palette_mode_info.palette_size[1] = 0;

  read_ref_frames(cm, xd, r, mbmi->segment_id, mbmi->ref_frame);
  is_compound = has_second_ref(mbmi);

  for (ref = 0; ref < 1 + is_compound; ++ref) {
    MV_REFERENCE_FRAME frame = mbmi->ref_frame[ref];
    RefBuffer *ref_buf = &cm->frame_refs[frame - LAST_FRAME];

    xd->block_refs[ref] = ref_buf;
    if ((!vp10_is_valid_scale(&ref_buf->sf)))
      vpx_internal_error(xd->error_info, VPX_CODEC_UNSUP_BITSTREAM,
                         "Reference frame has invalid dimensions");
    vp10_setup_pre_planes(xd, ref, ref_buf->buf, mi_row, mi_col,
                         &ref_buf->sf);
  }

  for (ref_frame = LAST_FRAME; ref_frame < MODE_CTX_REF_FRAMES; ++ref_frame) {
    vp10_find_mv_refs(cm, xd, mi, ref_frame,
#if CONFIG_REF_MV
                      &xd->ref_mv_count[ref_frame],
                      xd->ref_mv_stack[ref_frame],
#if CONFIG_EXT_INTER
                      compound_inter_mode_ctx,
#endif  // CONFIG_EXT_INTER
#endif
                      ref_mvs[ref_frame],
                      mi_row, mi_col, fpm_sync, (void *)pbi, inter_mode_ctx);
  }

#if CONFIG_REF_MV
#if CONFIG_EXT_INTER
  if (is_compound)
    mode_ctx = compound_inter_mode_ctx[mbmi->ref_frame[0]];
  else
#endif  // CONFIG_EXT_INTER
  mode_ctx = vp10_mode_context_analyzer(inter_mode_ctx,
                                        mbmi->ref_frame, bsize, -1);
  mbmi->ref_mv_idx = 0;
#else
  mode_ctx = inter_mode_ctx[mbmi->ref_frame[0]];
#endif

  if (segfeature_active(&cm->seg, mbmi->segment_id, SEG_LVL_SKIP)) {
    mbmi->mode = ZEROMV;
    if (bsize < BLOCK_8X8) {
        vpx_internal_error(xd->error_info, VPX_CODEC_UNSUP_BITSTREAM,
                           "Invalid usage of segement feature on small blocks");
        return;
    }
  } else {
    if (bsize >= BLOCK_8X8) {
#if CONFIG_EXT_INTER
      if (is_compound)
        mbmi->mode = read_inter_compound_mode(cm, xd, r, mode_ctx);
      else
#endif  // CONFIG_EXT_INTER
      mbmi->mode = read_inter_mode(cm, xd,
#if CONFIG_REF_MV && CONFIG_EXT_INTER
                                   mbmi,
#endif  // CONFIG_REF_MV && CONFIG_EXT_INTER
                                   r, mode_ctx);
#if CONFIG_REF_MV
      if (mbmi->mode == NEARMV || mbmi->mode == NEWMV)
        read_drl_idx(cm, xd, mbmi, r);
#endif
    }
  }

#if CONFIG_EXT_INTER
  if (bsize < BLOCK_8X8 ||
      (mbmi->mode != ZEROMV && mbmi->mode != ZERO_ZEROMV)) {
#else
  if (bsize < BLOCK_8X8 || mbmi->mode != ZEROMV) {
#endif  // CONFIG_EXT_INTER
    for (ref = 0; ref < 1 + is_compound; ++ref) {
      vp10_find_best_ref_mvs(allow_hp, ref_mvs[mbmi->ref_frame[ref]],
                             &nearestmv[ref], &nearmv[ref]);
    }
  }

#if CONFIG_REF_MV
  if (mbmi->ref_mv_idx > 0) {
    int_mv cur_mv =
        xd->ref_mv_stack[mbmi->ref_frame[0]][1 + mbmi->ref_mv_idx].this_mv;
    nearmv[0] = cur_mv;
  }

#if CONFIG_EXT_INTER
  if (is_compound && bsize >= BLOCK_8X8 && mbmi->mode != ZERO_ZEROMV) {
#else
  if (is_compound && bsize >= BLOCK_8X8 && mbmi->mode != NEWMV &&
      mbmi->mode != ZEROMV) {
#endif  // CONFIG_EXT_INTER
    uint8_t ref_frame_type = vp10_ref_frame_type(mbmi->ref_frame);

#if CONFIG_EXT_INTER
    if (xd->ref_mv_count[ref_frame_type] > 0) {
#else
    if (xd->ref_mv_count[ref_frame_type] == 1 && mbmi->mode == NEARESTMV) {
#endif  // CONFIG_EXT_INTER
#if CONFIG_EXT_INTER
      if (mbmi->mode == NEAREST_NEARESTMV) {
#endif  // CONFIG_EXT_INTER
        nearestmv[0] = xd->ref_mv_stack[ref_frame_type][0].this_mv;
        nearestmv[1] = xd->ref_mv_stack[ref_frame_type][0].comp_mv;
        lower_mv_precision(&nearestmv[0].as_mv, allow_hp);
        lower_mv_precision(&nearestmv[1].as_mv, allow_hp);
#if CONFIG_EXT_INTER
      } else if (mbmi->mode == NEAREST_NEWMV || mbmi->mode == NEAREST_NEARMV) {
        nearestmv[0] = xd->ref_mv_stack[ref_frame_type][0].this_mv;
        lower_mv_precision(&nearestmv[0].as_mv, allow_hp);
      } else if (mbmi->mode == NEW_NEARESTMV || mbmi->mode == NEAR_NEARESTMV) {
        nearestmv[1] = xd->ref_mv_stack[ref_frame_type][0].comp_mv;
        lower_mv_precision(&nearestmv[1].as_mv, allow_hp);
      }
#endif  // CONFIG_EXT_INTER
    }

#if CONFIG_EXT_INTER
    if (xd->ref_mv_count[ref_frame_type] > 1) {
      if (mbmi->mode == NEAR_NEWMV ||
          mbmi->mode == NEAR_NEARESTMV ||
          mbmi->mode == NEAR_NEARMV) {
        nearmv[0] = xd->ref_mv_stack[ref_frame_type][1].this_mv;
        lower_mv_precision(&nearmv[0].as_mv, allow_hp);
      }

      if (mbmi->mode == NEW_NEARMV ||
          mbmi->mode == NEAREST_NEARMV ||
          mbmi->mode == NEAR_NEARMV) {
        nearmv[1] = xd->ref_mv_stack[ref_frame_type][1].comp_mv;
        lower_mv_precision(&nearmv[1].as_mv, allow_hp);
      }
    }
#else
    if (xd->ref_mv_count[ref_frame_type] > 1) {
      int ref_mv_idx = 1 + mbmi->ref_mv_idx;
      nearestmv[0] = xd->ref_mv_stack[ref_frame_type][0].this_mv;
      nearestmv[1] = xd->ref_mv_stack[ref_frame_type][0].comp_mv;
      nearmv[0] = xd->ref_mv_stack[ref_frame_type][ref_mv_idx].this_mv;
      nearmv[1] = xd->ref_mv_stack[ref_frame_type][ref_mv_idx].comp_mv;
    }
#endif  // CONFIG_EXT_INTER
  }
#endif

#if !CONFIG_EXT_INTERP && !CONFIG_DUAL_FILTER
  mbmi->interp_filter = read_interp_filter(cm, xd, r);
#endif  // !CONFIG_EXT_INTERP && !CONFIG_DUAL_FILTER

  if (bsize < BLOCK_8X8) {
    const int num_4x4_w = 1 << xd->bmode_blocks_wl;
    const int num_4x4_h = 1 << xd->bmode_blocks_hl;
    int idx, idy;
    PREDICTION_MODE b_mode;
    int_mv nearest_sub8x8[2], near_sub8x8[2];
#if CONFIG_EXT_INTER
    int_mv ref_mv[2][2];
#endif  // CONFIG_EXT_INTER
    for (idy = 0; idy < 2; idy += num_4x4_h) {
      for (idx = 0; idx < 2; idx += num_4x4_w) {
        int_mv block[2];
        const int j = idy * 2 + idx;
        int_mv ref_mv_s8[2];
#if CONFIG_REF_MV
#if CONFIG_EXT_INTER
        if (!is_compound)
#endif  // CONFIG_EXT_INTER
        mode_ctx = vp10_mode_context_analyzer(inter_mode_ctx,  mbmi->ref_frame,
                                              bsize, j);
#endif
#if CONFIG_EXT_INTER
        if (is_compound)
          b_mode = read_inter_compound_mode(cm, xd, r, mode_ctx);
        else
#endif  // CONFIG_EXT_INTER
        b_mode = read_inter_mode(cm, xd,
#if CONFIG_REF_MV && CONFIG_EXT_INTER
                                 mbmi,
#endif  // CONFIG_REF_MV && CONFIG_EXT_INTER
                                 r, mode_ctx);

#if CONFIG_EXT_INTER
        mv_idx = (b_mode == NEWFROMNEARMV) ? 1 : 0;

        if (b_mode != ZEROMV && b_mode != ZERO_ZEROMV) {
#else
        if (b_mode != ZEROMV) {
#endif  // CONFIG_EXT_INTER
#if CONFIG_REF_MV
          CANDIDATE_MV ref_mv_stack[2][MAX_REF_MV_STACK_SIZE];
          uint8_t ref_mv_count[2];
#endif
          for (ref = 0; ref < 1 + is_compound; ++ref)
#if CONFIG_EXT_INTER
          {
            int_mv mv_ref_list[MAX_MV_REF_CANDIDATES];
            vp10_update_mv_context(xd, mi, mbmi->ref_frame[ref],
                                   mv_ref_list, j, mi_row, mi_col, NULL);
#endif  // CONFIG_EXT_INTER
            vp10_append_sub8x8_mvs_for_idx(cm, xd, j, ref, mi_row, mi_col,
#if CONFIG_REF_MV
                                           ref_mv_stack[ref],
                                           &ref_mv_count[ref],
#endif
#if CONFIG_EXT_INTER
                                           mv_ref_list,
#endif  // CONFIG_EXT_INTER
                                           &nearest_sub8x8[ref],
                                           &near_sub8x8[ref]);
#if CONFIG_EXT_INTER
            if (have_newmv_in_inter_mode(b_mode)) {
              mv_ref_list[0].as_int = nearest_sub8x8[ref].as_int;
              mv_ref_list[1].as_int = near_sub8x8[ref].as_int;
              vp10_find_best_ref_mvs(allow_hp, mv_ref_list,
                                     &ref_mv[0][ref], &ref_mv[1][ref]);
            }
          }
#endif  // CONFIG_EXT_INTER
        }

        for (ref = 0; ref < 1 + is_compound && b_mode != ZEROMV; ++ref) {
#if CONFIG_REF_MV
          ref_mv_s8[ref] = nearest_sub8x8[ref];
          lower_mv_precision(&ref_mv_s8[ref].as_mv, allow_hp);
#else
          ref_mv_s8[ref] = nearestmv[ref];
#endif
        }
#if CONFIG_EXT_INTER
        (void)ref_mv_s8;
#endif

        if (!assign_mv(cm, xd, b_mode,
#if CONFIG_REF_MV
                       j,
#endif
                       block,
#if CONFIG_EXT_INTER
                       ref_mv[mv_idx],
#else
                       ref_mv_s8,
#endif  // CONFIG_EXT_INTER
                       nearest_sub8x8, near_sub8x8,
                       is_compound, allow_hp, r)) {
          xd->corrupted |= 1;
          break;
        };

        mi->bmi[j].as_mv[0].as_int = block[0].as_int;
        if (is_compound)
          mi->bmi[j].as_mv[1].as_int = block[1].as_int;

        if (num_4x4_h == 2)
          mi->bmi[j + 2] = mi->bmi[j];
        if (num_4x4_w == 2)
          mi->bmi[j + 1] = mi->bmi[j];
      }
    }

#if CONFIG_REF_MV
    mbmi->pred_mv[0].as_int = mi->bmi[3].pred_mv_s8[0].as_int;
    mbmi->pred_mv[1].as_int = mi->bmi[3].pred_mv_s8[1].as_int;
#endif
    mi->mbmi.mode = b_mode;

    mbmi->mv[0].as_int = mi->bmi[3].as_mv[0].as_int;
    mbmi->mv[1].as_int = mi->bmi[3].as_mv[1].as_int;
  } else {
    int ref;
    int_mv ref_mv[2];
    ref_mv[0] = nearestmv[0];
    ref_mv[1] = nearestmv[1];

    for (ref = 0; ref < 1 + is_compound && mbmi->mode == NEWMV; ++ref) {
#if CONFIG_REF_MV
      uint8_t ref_frame_type = vp10_ref_frame_type(mbmi->ref_frame);
      if (xd->ref_mv_count[ref_frame_type] > 1) {
        ref_mv[ref] = (ref == 0) ?
            xd->ref_mv_stack[ref_frame_type][mbmi->ref_mv_idx].this_mv :
            xd->ref_mv_stack[ref_frame_type][mbmi->ref_mv_idx].comp_mv;
        clamp_mv_ref(&ref_mv[ref].as_mv, xd->n8_w << 3, xd->n8_h << 3, xd);
      }
#endif
      nearestmv[ref] = ref_mv[ref];
    }

    xd->corrupted |= !assign_mv(cm, xd, mbmi->mode,
#if CONFIG_REF_MV
                                0,
#endif
                                mbmi->mv,
#if CONFIG_EXT_INTER
                                mbmi->mode == NEWFROMNEARMV ?
                                              nearmv : nearestmv,
#else
                                ref_mv,
#endif  // CONFIG_EXT_INTER
                                nearestmv, nearmv, is_compound, allow_hp, r);
  }

#if CONFIG_EXT_INTER
  mbmi->use_wedge_interintra = 0;
  if (cm->reference_mode != COMPOUND_REFERENCE &&
#if CONFIG_SUPERTX
      !supertx_enabled &&
#endif
      is_interintra_allowed(mbmi)) {
    const int bsize_group = size_group_lookup[bsize];
    const int interintra = vp10_read(r, cm->fc->interintra_prob[bsize_group]);
    if (xd->counts)
      xd->counts->interintra[bsize_group][interintra]++;
    assert(mbmi->ref_frame[1] == NONE);
    if (interintra) {
      const INTERINTRA_MODE interintra_mode =
          read_interintra_mode(cm, xd, r, bsize_group);
      mbmi->ref_frame[1] = INTRA_FRAME;
      mbmi->interintra_mode = interintra_mode;
#if CONFIG_EXT_INTRA
      mbmi->ext_intra_mode_info.use_ext_intra_mode[0] = 0;
      mbmi->ext_intra_mode_info.use_ext_intra_mode[1] = 0;
      mbmi->angle_delta[0] = 0;
      mbmi->angle_delta[1] = 0;
      mbmi->intra_filter = INTRA_FILTER_LINEAR;
#endif  // CONFIG_EXT_INTRA
      if (is_interintra_wedge_used(bsize)) {
        mbmi->use_wedge_interintra =
            vp10_read(r, cm->fc->wedge_interintra_prob[bsize]);
        if (xd->counts)
          xd->counts->wedge_interintra[bsize][mbmi->use_wedge_interintra]++;
        if (mbmi->use_wedge_interintra) {
          mbmi->interintra_wedge_index =
              vp10_read_literal(r, get_wedge_bits_lookup(bsize));
          mbmi->interintra_wedge_sign = 0;
        }
      }
    }
  }
#endif  // CONFIG_EXT_INTER

#if CONFIG_OBMC || CONFIG_WARPED_MOTION
  mbmi->motion_variation = SIMPLE_TRANSLATION;
#if CONFIG_SUPERTX
  if (!supertx_enabled)
#endif  // CONFIG_SUPERTX
#if CONFIG_EXT_INTER
    if (mbmi->ref_frame[1] != INTRA_FRAME)
#endif  // CONFIG_EXT_INTER
    mbmi->motion_variation = read_motvar_block(cm, xd, r);
#endif  // CONFIG_OBMC || CONFIG_WARPED_MOTION

#if CONFIG_EXT_INTER
  mbmi->use_wedge_interinter = 0;
  if (cm->reference_mode != SINGLE_REFERENCE &&
      is_inter_compound_mode(mbmi->mode) &&
#if CONFIG_OBMC || CONFIG_WARPED_MOTION
      !(is_motvar_allowed(mbmi) &&
        mbmi->motion_variation != SIMPLE_TRANSLATION) &&
#endif  // CONFIG_OBMC || CONFIG_WARPED_MOTION
      is_interinter_wedge_used(bsize)) {
    mbmi->use_wedge_interinter =
        vp10_read(r, cm->fc->wedge_interinter_prob[bsize]);
    if (xd->counts)
      xd->counts->wedge_interinter[bsize][mbmi->use_wedge_interinter]++;
    if (mbmi->use_wedge_interinter) {
      mbmi->interinter_wedge_index =
          vp10_read_literal(r, get_wedge_bits_lookup(bsize));
      mbmi->interinter_wedge_sign = vp10_read_bit(r);
    }
  }
#endif  // CONFIG_EXT_INTER

#if CONFIG_DUAL_FILTER
  for (ref = 0; ref < 2; ++ref) {
    mbmi->interp_filter[ref] = (cm->interp_filter == SWITCHABLE) ?
        EIGHTTAP_REGULAR : cm->interp_filter;

    if (has_subpel_mv_component(xd->mi[0], xd, ref) ||
        (mbmi->ref_frame[1] > INTRA_FRAME &&
         has_subpel_mv_component(xd->mi[0], xd, ref + 2)))
      mbmi->interp_filter[ref] = read_interp_filter(cm, xd, ref, r);
  }
  // The index system worsk as:
  // (0, 1) -> (vertical, horizontal) filter types for the first ref frame.
  // (2, 3) -> (vertical, horizontal) filter types for the second ref frame.
  mbmi->interp_filter[2] = mbmi->interp_filter[0];
  mbmi->interp_filter[3] = mbmi->interp_filter[1];
#else
#if CONFIG_EXT_INTERP
  mbmi->interp_filter = read_interp_filter(cm, xd, r);
#endif  // CONFIG_EXT_INTERP
#endif  // CONFIG_DUAL_FILTER
}

static void read_inter_frame_mode_info(VP10Decoder *const pbi,
                                       MACROBLOCKD *const xd,
#if CONFIG_SUPERTX
                                       int supertx_enabled,
#endif  // CONFIG_SUPERTX
                                       int mi_row, int mi_col, vp10_reader *r) {
  VP10_COMMON *const cm = &pbi->common;
  MODE_INFO *const mi = xd->mi[0];
  MB_MODE_INFO *const mbmi = &mi->mbmi;
  int inter_block = 1;
#if CONFIG_VAR_TX
  BLOCK_SIZE bsize = mbmi->sb_type;
#endif  // CONFIG_VAR_TX

  mbmi->mv[0].as_int = 0;
  mbmi->mv[1].as_int = 0;
  mbmi->segment_id = read_inter_segment_id(cm, xd, mi_row, mi_col, r);
#if CONFIG_SUPERTX
  if (!supertx_enabled) {
#endif  // CONFIG_SUPERTX
    mbmi->skip = read_skip(cm, xd, mbmi->segment_id, r);
    inter_block = read_is_inter_block(cm, xd, mbmi->segment_id, r);

#if CONFIG_VAR_TX
    xd->above_txfm_context = cm->above_txfm_context + mi_col;
    xd->left_txfm_context =
      xd->left_txfm_context_buffer + (mi_row & MAX_MIB_MASK);
    if (bsize >= BLOCK_8X8 && cm->tx_mode == TX_MODE_SELECT &&
        !mbmi->skip && inter_block) {
      const TX_SIZE max_tx_size = max_txsize_lookup[bsize];
      const BLOCK_SIZE txb_size = txsize_to_bsize[max_tx_size];
      const int bs = num_4x4_blocks_wide_lookup[txb_size];
      const int width  = num_4x4_blocks_wide_lookup[bsize];
      const int height = num_4x4_blocks_high_lookup[bsize];
      int idx, idy;
      for (idy = 0; idy < height; idy += bs)
        for (idx = 0; idx < width; idx += bs)
          read_tx_size_vartx(cm, xd, mbmi, xd->counts, max_tx_size,
                             idy, idx, r);
      if (xd->counts) {
        const int ctx = get_tx_size_context(xd);
        ++xd->counts->tx_size[max_tx_size - TX_8X8][ctx][mbmi->tx_size];
      }
    } else {
      if (inter_block)
        mbmi->tx_size = read_tx_size_inter(cm, xd, !mbmi->skip, r);
      else
        mbmi->tx_size = read_tx_size_intra(cm, xd, r);

      if (inter_block) {
        const int width  = num_4x4_blocks_wide_lookup[bsize];
        const int height = num_4x4_blocks_high_lookup[bsize];
        int idx, idy;
        for (idy = 0; idy < height; ++idy)
          for (idx = 0; idx < width; ++idx)
            mbmi->inter_tx_size[idy >> 1][idx >> 1] = mbmi->tx_size;
      }

      set_txfm_ctx(xd->left_txfm_context, mbmi->tx_size, xd->n8_h);
      set_txfm_ctx(xd->above_txfm_context, mbmi->tx_size, xd->n8_w);
    }
#else
    if (inter_block)
      mbmi->tx_size = read_tx_size_inter(cm, xd, !mbmi->skip, r);
    else
      mbmi->tx_size = read_tx_size_intra(cm, xd, r);
#endif  // CONFIG_VAR_TX
#if CONFIG_SUPERTX
  }
#if CONFIG_VAR_TX
  else if (inter_block) {
    const int width  = num_4x4_blocks_wide_lookup[bsize];
    const int height = num_4x4_blocks_high_lookup[bsize];
    int idx, idy;
    xd->mi[0]->mbmi.tx_size = xd->supertx_size;
    for (idy = 0; idy < height; ++idy)
      for (idx = 0; idx < width; ++idx)
        xd->mi[0]->mbmi.inter_tx_size[idy >> 1][idx >> 1] =
            xd->supertx_size;
  }
#endif  // CONFIG_VAR_TX
#endif  // CONFIG_SUPERTX

  if (inter_block)
    read_inter_block_mode_info(pbi, xd,
#if (CONFIG_OBMC || CONFIG_EXT_INTER) && CONFIG_SUPERTX

                               mi, mi_row, mi_col, r, supertx_enabled);
#else
                               mi, mi_row, mi_col, r);
#endif  // CONFIG_OBMC && CONFIG_SUPERTX
  else
    read_intra_block_mode_info(cm, xd, mi, r);

  if (!FIXED_TX_TYPE) {
#if CONFIG_EXT_TX
    if (get_ext_tx_types(mbmi->tx_size, mbmi->sb_type, inter_block) > 1 &&
        cm->base_qindex > 0 && !mbmi->skip &&
#if CONFIG_SUPERTX
        !supertx_enabled &&
#endif  // CONFIG_SUPERTX
        !segfeature_active(&cm->seg, mbmi->segment_id, SEG_LVL_SKIP)) {
      int eset = get_ext_tx_set(mbmi->tx_size, mbmi->sb_type,
                                inter_block);
      FRAME_COUNTS *counts = xd->counts;

      if (inter_block) {
        if (eset > 0) {
          mbmi->tx_type =
              vp10_read_tree(r, vp10_ext_tx_inter_tree[eset],
                            cm->fc->inter_ext_tx_prob[eset][mbmi->tx_size]);
          if (counts)
            ++counts->inter_ext_tx[eset][mbmi->tx_size][mbmi->tx_type];
        }
      } else if (ALLOW_INTRA_EXT_TX) {
        if (eset > 0) {
          mbmi->tx_type = vp10_read_tree(r, vp10_ext_tx_intra_tree[eset],
                                        cm->fc->intra_ext_tx_prob[eset]
                                                [mbmi->tx_size][mbmi->mode]);
          if (counts)
            ++counts->intra_ext_tx[eset][mbmi->tx_size]
                                         [mbmi->mode][mbmi->tx_type];
        }
      }
    } else {
      mbmi->tx_type = DCT_DCT;
    }
#else
    if (mbmi->tx_size < TX_32X32 &&
        cm->base_qindex > 0 && !mbmi->skip &&
#if CONFIG_SUPERTX
        !supertx_enabled &&
#endif  // CONFIG_SUPERTX
        !segfeature_active(&cm->seg, mbmi->segment_id, SEG_LVL_SKIP)) {
      FRAME_COUNTS *counts = xd->counts;
      if (inter_block) {
        mbmi->tx_type = vp10_read_tree(
            r, vp10_ext_tx_tree,
            cm->fc->inter_ext_tx_prob[mbmi->tx_size]);
        if (counts)
          ++counts->inter_ext_tx[mbmi->tx_size][mbmi->tx_type];
      } else {
        const TX_TYPE tx_type_nom = intra_mode_to_tx_type_context[mbmi->mode];
        mbmi->tx_type = vp10_read_tree(
            r, vp10_ext_tx_tree,
            cm->fc->intra_ext_tx_prob[mbmi->tx_size][tx_type_nom]);
        if (counts)
          ++counts->intra_ext_tx[mbmi->tx_size][tx_type_nom][mbmi->tx_type];
      }
    } else {
      mbmi->tx_type = DCT_DCT;
    }
#endif  // CONFIG_EXT_TX
  }
}

void vp10_read_mode_info(VP10Decoder *const pbi, MACROBLOCKD *xd,
#if CONFIG_SUPERTX
                         int supertx_enabled,
#endif  // CONFIG_SUPERTX
                         int mi_row, int mi_col, vp10_reader *r,
                         int x_mis, int y_mis) {
  VP10_COMMON *const cm = &pbi->common;
  MODE_INFO *const mi = xd->mi[0];
  MV_REF* frame_mvs = cm->cur_frame->mvs + mi_row * cm->mi_cols + mi_col;
  int w, h;

  if (frame_is_intra_only(cm)) {
    read_intra_frame_mode_info(cm, xd, mi_row, mi_col, r);
#if CONFIG_REF_MV
    for (h = 0; h < y_mis; ++h) {
      MV_REF *const frame_mv = frame_mvs + h * cm->mi_cols;
      for (w = 0; w < x_mis; ++w) {
        MV_REF *const mv = frame_mv + w;
        mv->ref_frame[0] = NONE;
        mv->ref_frame[1] = NONE;
      }
    }
#endif
  } else {
    read_inter_frame_mode_info(pbi, xd,
#if CONFIG_SUPERTX
                               supertx_enabled,
#endif  // CONFIG_SUPERTX
                               mi_row, mi_col, r);
    for (h = 0; h < y_mis; ++h) {
      MV_REF *const frame_mv = frame_mvs + h * cm->mi_cols;
      for (w = 0; w < x_mis; ++w) {
        MV_REF *const mv = frame_mv + w;
        mv->ref_frame[0] = mi->mbmi.ref_frame[0];
        mv->ref_frame[1] = mi->mbmi.ref_frame[1];
        mv->mv[0].as_int = mi->mbmi.mv[0].as_int;
        mv->mv[1].as_int = mi->mbmi.mv[1].as_int;
      }
    }
  }
}
