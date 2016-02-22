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

static INLINE int read_uniform(vpx_reader *r, int n) {
  int l = get_unsigned_bits(n);
  int m = (1 << l) - n;
  int v = vpx_read_literal(r, l-1);

  assert(l != 0);

  if (v < m)
    return v;
  else
    return (v << 1) - m + vpx_read_literal(r, 1);
}

static PREDICTION_MODE read_intra_mode(vpx_reader *r, const vpx_prob *p) {
  return (PREDICTION_MODE)vpx_read_tree(r, vp10_intra_mode_tree, p);
}

static PREDICTION_MODE read_intra_mode_y(VP10_COMMON *cm, MACROBLOCKD *xd,
                                         vpx_reader *r, int size_group) {
  const PREDICTION_MODE y_mode =
      read_intra_mode(r, cm->fc->y_mode_prob[size_group]);
  FRAME_COUNTS *counts = xd->counts;
  if (counts)
    ++counts->y_mode[size_group][y_mode];
  return y_mode;
}

static PREDICTION_MODE read_intra_mode_uv(VP10_COMMON *cm, MACROBLOCKD *xd,
                                          vpx_reader *r,
                                          PREDICTION_MODE y_mode) {
  const PREDICTION_MODE uv_mode = read_intra_mode(r,
                                         cm->fc->uv_mode_prob[y_mode]);
  FRAME_COUNTS *counts = xd->counts;
  if (counts)
    ++counts->uv_mode[y_mode][uv_mode];
  return uv_mode;
}

static PREDICTION_MODE read_inter_mode(VP10_COMMON *cm, MACROBLOCKD *xd,
#if CONFIG_REF_MV && CONFIG_EXT_INTER
                                       MB_MODE_INFO *mbmi,
#endif
                                       vpx_reader *r, int16_t ctx) {
#if CONFIG_REF_MV
  FRAME_COUNTS *counts = xd->counts;
  int16_t mode_ctx = ctx & NEWMV_CTX_MASK;
  vpx_prob mode_prob = cm->fc->newmv_prob[mode_ctx];

  if (vpx_read(r, mode_prob) == 0) {
    if (counts)
      ++counts->newmv_mode[mode_ctx][0];

#if CONFIG_EXT_INTER
    if (has_second_ref(mbmi)) {
#endif  // CONFIG_EXT_INTER
    return NEWMV;
#if CONFIG_EXT_INTER
    } else {
      mode_prob = cm->fc->new2mv_prob;
      if (vpx_read(r, mode_prob) == 0) {
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
  if (vpx_read(r, mode_prob) == 0) {
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

  if (vpx_read(r, mode_prob) == 0) {
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
  const int mode = vpx_read_tree(r, vp10_inter_mode_tree,
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
                         vpx_reader *r) {
  uint8_t ref_frame_type = vp10_ref_frame_type(mbmi->ref_frame);
  mbmi->ref_mv_idx = 0;

  if (xd->ref_mv_count[ref_frame_type] > 2) {
    uint8_t drl0_ctx = vp10_drl_ctx(xd->ref_mv_stack[ref_frame_type], 1);
    vpx_prob drl0_prob = cm->fc->drl_prob0[drl0_ctx];
    if (vpx_read(r, drl0_prob)) {
      mbmi->ref_mv_idx = 1;
      if (xd->counts)
        ++xd->counts->drl_mode0[drl0_ctx][1];
      if (xd->ref_mv_count[ref_frame_type] > 3) {
        uint8_t drl1_ctx =
            vp10_drl_ctx(xd->ref_mv_stack[ref_frame_type], 2);
        vpx_prob drl1_prob = cm->fc->drl_prob1[drl1_ctx];
        if (vpx_read(r, drl1_prob)) {
          mbmi->ref_mv_idx = 2;
          if (xd->counts)
            ++xd->counts->drl_mode1[drl1_ctx][1];

          return;
        }

        if (xd->counts)
          ++xd->counts->drl_mode1[drl1_ctx][0];
      }
      return;
    }

    if (xd->counts)
      ++xd->counts->drl_mode0[drl0_ctx][0];
  }
}
#endif

#if CONFIG_EXT_INTER
static PREDICTION_MODE read_inter_compound_mode(VP10_COMMON *cm,
                                                MACROBLOCKD *xd,
                                                vpx_reader *r, int16_t ctx) {
  const int mode = vpx_read_tree(r, vp10_inter_compound_mode_tree,
                                 cm->fc->inter_compound_mode_probs[ctx]);
  FRAME_COUNTS *counts = xd->counts;

  if (counts)
    ++counts->inter_compound_mode[ctx][mode];

  assert(is_inter_compound_mode(NEAREST_NEARESTMV + mode));
  return NEAREST_NEARESTMV + mode;
}
#endif  // CONFIG_EXT_INTER

static int read_segment_id(vpx_reader *r,
    const struct segmentation_probs *segp) {
  return vpx_read_tree(r, vp10_segment_tree, segp->tree_probs);
}

#if CONFIG_VAR_TX
static void read_tx_size_inter(VP10_COMMON *cm, MACROBLOCKD *xd,
                               MB_MODE_INFO *mbmi, FRAME_COUNTS *counts,
                               TX_SIZE tx_size, int blk_row, int blk_col,
                               vpx_reader *r) {
  int is_split = 0;
  const int tx_idx = (blk_row >> 1) * 8 + (blk_col >> 1);
  int max_blocks_high = num_4x4_blocks_high_lookup[mbmi->sb_type];
  int max_blocks_wide = num_4x4_blocks_wide_lookup[mbmi->sb_type];
  int ctx = txfm_partition_context(xd->above_txfm_context + (blk_col >> 1),
                                   xd->left_txfm_context + (blk_row >> 1),
                                   tx_size);

  if (xd->mb_to_bottom_edge < 0)
    max_blocks_high += xd->mb_to_bottom_edge >> 5;
  if (xd->mb_to_right_edge < 0)
     max_blocks_wide += xd->mb_to_right_edge >> 5;

  if (blk_row >= max_blocks_high || blk_col >= max_blocks_wide)
     return;

  is_split = vpx_read(r, cm->fc->txfm_partition_prob[ctx]);

  if (is_split) {
    BLOCK_SIZE bsize = txsize_to_bsize[tx_size];
    int bsl = b_width_log2_lookup[bsize];
    int i;

    if (counts)
      ++counts->txfm_partition[ctx][1];

    if (tx_size == TX_8X8) {
      mbmi->inter_tx_size[tx_idx] = TX_4X4;
      mbmi->tx_size = mbmi->inter_tx_size[tx_idx];
      txfm_partition_update(xd->above_txfm_context + (blk_col >> 1),
                            xd->left_txfm_context + (blk_row >> 1), TX_4X4);
      return;
    }

    assert(bsl > 0);
    --bsl;
    for (i = 0; i < 4; ++i) {
      int offsetr = blk_row + ((i >> 1) << bsl);
      int offsetc = blk_col + ((i & 0x01) << bsl);
      read_tx_size_inter(cm, xd, mbmi, counts,
                         tx_size - 1, offsetr, offsetc, r);
    }
  } else {
    int idx, idy;
    mbmi->inter_tx_size[tx_idx] = tx_size;
    for (idy = 0; idy < (1 << tx_size) / 2; ++idy)
      for (idx = 0; idx < (1 << tx_size) / 2; ++idx)
        mbmi->inter_tx_size[tx_idx + (idy << 3) + idx] = tx_size;
    mbmi->tx_size = mbmi->inter_tx_size[tx_idx];
    if (counts)
      ++counts->txfm_partition[ctx][0];
    txfm_partition_update(xd->above_txfm_context + (blk_col >> 1),
                          xd->left_txfm_context + (blk_row >> 1), tx_size);
  }
}
#endif

static TX_SIZE read_selected_tx_size(VP10_COMMON *cm, MACROBLOCKD *xd,
                                     TX_SIZE max_tx_size, vpx_reader *r) {
  FRAME_COUNTS *counts = xd->counts;
  const int ctx = get_tx_size_context(xd);
  const vpx_prob *tx_probs = get_tx_probs(max_tx_size, ctx, &cm->fc->tx_probs);
  int tx_size = vpx_read(r, tx_probs[0]);
  if (tx_size != TX_4X4 && max_tx_size >= TX_16X16) {
    tx_size += vpx_read(r, tx_probs[1]);
    if (tx_size != TX_8X8 && max_tx_size >= TX_32X32)
      tx_size += vpx_read(r, tx_probs[2]);
  }

  if (counts)
    ++get_tx_counts(max_tx_size, ctx, &counts->tx)[tx_size];
  return (TX_SIZE)tx_size;
}

static TX_SIZE read_tx_size(VP10_COMMON *cm, MACROBLOCKD *xd,
                            int allow_select, vpx_reader *r) {
  TX_MODE tx_mode = cm->tx_mode;
  BLOCK_SIZE bsize = xd->mi[0]->mbmi.sb_type;
  const TX_SIZE max_tx_size = max_txsize_lookup[bsize];
  if (xd->lossless[xd->mi[0]->mbmi.segment_id])
    return TX_4X4;
  if (allow_select && tx_mode == TX_MODE_SELECT && bsize >= BLOCK_8X8)
    return read_selected_tx_size(cm, xd, max_tx_size, r);
  else
    return VPXMIN(max_tx_size, tx_mode_to_biggest_tx_size[tx_mode]);
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
                                 vpx_reader *r) {
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
                                 int mi_row, int mi_col, vpx_reader *r) {
  struct segmentation *const seg = &cm->seg;
  FRAME_COUNTS *counts = xd->counts;
  struct segmentation_probs *const segp = &cm->fc->seg;
  MB_MODE_INFO *const mbmi = &xd->mi[0]->mbmi;
  int predicted_segment_id, segment_id;
  const int mi_offset = mi_row * cm->mi_cols + mi_col;
  const int bw = xd->plane[0].n4_w >> 1;
  const int bh = xd->plane[0].n4_h >> 1;

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
    mbmi->seg_id_predicted = vpx_read(r, pred_prob);
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
                     int segment_id, vpx_reader *r) {
  if (segfeature_active(&cm->seg, segment_id, SEG_LVL_SKIP)) {
    return 1;
  } else {
    const int ctx = vp10_get_skip_context(xd);
    const int skip = vpx_read(r, cm->fc->skip_probs[ctx]);
    FRAME_COUNTS *counts = xd->counts;
    if (counts)
      ++counts->skip[ctx][skip];
    return skip;
  }
}

static void read_palette_mode_info(VP10_COMMON *const cm,
                                   MACROBLOCKD *const xd,
                                   vpx_reader *r) {
  MODE_INFO *const mi = xd->mi[0];
  MB_MODE_INFO *const mbmi = &mi->mbmi;
  const MODE_INFO *above_mi = xd->above_mi;
  const MODE_INFO *left_mi  = xd->left_mi;
  const BLOCK_SIZE bsize = mbmi->sb_type;
  int i, palette_ctx = 0;

  if (above_mi)
    palette_ctx += (above_mi->mbmi.palette_mode_info.palette_size[0] > 0);
  if (left_mi)
    palette_ctx += (left_mi->mbmi.palette_mode_info.palette_size[0] > 0);
  if (vpx_read(r, vp10_default_palette_y_mode_prob[bsize - BLOCK_8X8]
                                                   [palette_ctx])) {
    int n;
    PALETTE_MODE_INFO *pmi = &mbmi->palette_mode_info;

    pmi->palette_size[0] =
        vpx_read_tree(r, vp10_palette_size_tree,
                      vp10_default_palette_y_size_prob[bsize - BLOCK_8X8]) + 2;
    n = pmi->palette_size[0];

    for (i = 0; i < n; ++i)
      pmi->palette_colors[i] = vpx_read_literal(r, cm->bit_depth);

    xd->plane[0].color_index_map[0] = read_uniform(r, n);
    assert(xd->plane[0].color_index_map[0] < n);
  }
}

#if CONFIG_EXT_INTRA
static void read_ext_intra_mode_info(VP10_COMMON *const cm,
                                     MACROBLOCKD *const xd, vpx_reader *r) {
  MODE_INFO *const mi = xd->mi[0];
  MB_MODE_INFO *const mbmi = &mi->mbmi;
  FRAME_COUNTS *counts = xd->counts;

#if !ALLOW_FILTER_INTRA_MODES
  return;
#endif
  if (mbmi->mode == DC_PRED) {
    mbmi->ext_intra_mode_info.use_ext_intra_mode[0] =
        vpx_read(r, cm->fc->ext_intra_probs[0]);
    if (mbmi->ext_intra_mode_info.use_ext_intra_mode[0]) {
      mbmi->ext_intra_mode_info.ext_intra_mode[0] =
          read_uniform(r, FILTER_INTRA_MODES);
    }
    if (counts)
      ++counts->ext_intra[0][mbmi->ext_intra_mode_info.use_ext_intra_mode[0]];
  }
  if (mbmi->uv_mode == DC_PRED) {
    mbmi->ext_intra_mode_info.use_ext_intra_mode[1] =
        vpx_read(r, cm->fc->ext_intra_probs[1]);
    if (mbmi->ext_intra_mode_info.use_ext_intra_mode[1]) {
      mbmi->ext_intra_mode_info.ext_intra_mode[1] =
          read_uniform(r, FILTER_INTRA_MODES);
    }
    if (counts)
      ++counts->ext_intra[1][mbmi->ext_intra_mode_info.use_ext_intra_mode[1]];
  }
}
#endif  // CONFIG_EXT_INTRA

static void read_intra_frame_mode_info(VP10_COMMON *const cm,
                                       MACROBLOCKD *const xd,
                                       int mi_row, int mi_col, vpx_reader *r) {
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
  mbmi->tx_size = read_tx_size(cm, xd, 1, r);
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
#if CONFIG_EXT_INTRA
      if (mbmi->mode != DC_PRED && mbmi->mode != TM_PRED) {
        int p_angle;
        const int ctx = vp10_get_pred_context_intra_interp(xd);
        mbmi->angle_delta[0] =
            read_uniform(r, 2 * MAX_ANGLE_DELTAS + 1) - MAX_ANGLE_DELTAS;
        p_angle = mode_to_angle_map[mbmi->mode] +
            mbmi->angle_delta[0] * ANGLE_STEP;
        if (pick_intra_filter(p_angle)) {
          FRAME_COUNTS *counts = xd->counts;
          mbmi->intra_filter = vpx_read_tree(r, vp10_intra_filter_tree,
                                             cm->fc->intra_filter_probs[ctx]);
          if (counts)
            ++counts->intra_filter[ctx][mbmi->intra_filter];
        } else {
          mbmi->intra_filter = INTRA_FILTER_LINEAR;
        }
      }
#endif  // CONFIG_EXT_INTRA
  }

  mbmi->uv_mode = read_intra_mode_uv(cm, xd, r, mbmi->mode);
#if CONFIG_EXT_INTRA
  if (mbmi->uv_mode != DC_PRED && mbmi->uv_mode != TM_PRED &&
      bsize >= BLOCK_8X8)
    mbmi->angle_delta[1] =
        read_uniform(r, 2 * MAX_ANGLE_DELTAS + 1) - MAX_ANGLE_DELTAS;
#endif

  mbmi->palette_mode_info.palette_size[0] = 0;
  mbmi->palette_mode_info.palette_size[1] = 0;
  if (bsize >= BLOCK_8X8 && cm->allow_screen_content_tools &&
      mbmi->mode == DC_PRED)
    read_palette_mode_info(cm, xd, r);

  if (!FIXED_TX_TYPE) {
#if CONFIG_EXT_TX
    if (get_ext_tx_types(mbmi->tx_size, mbmi->sb_type, 0) > 1 &&
        cm->base_qindex > 0 && !mbmi->skip &&
        !segfeature_active(&cm->seg, mbmi->segment_id, SEG_LVL_SKIP) &&
        ALLOW_INTRA_EXT_TX) {
      FRAME_COUNTS *counts = xd->counts;
      int eset = get_ext_tx_set(mbmi->tx_size, mbmi->sb_type, 0);
      if (eset > 0) {
        mbmi->tx_type = vpx_read_tree(
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
      mbmi->tx_type = vpx_read_tree(
          r, vp10_ext_tx_tree,
          cm->fc->intra_ext_tx_prob[mbmi->tx_size][tx_type_nom]);
      if (counts)
        ++counts->intra_ext_tx[mbmi->tx_size][tx_type_nom][mbmi->tx_type];
    } else {
      mbmi->tx_type = DCT_DCT;
    }
#endif  // CONFIG_EXT_TX
  }

#if CONFIG_EXT_INTRA
    mbmi->ext_intra_mode_info.use_ext_intra_mode[0] = 0;
    mbmi->ext_intra_mode_info.use_ext_intra_mode[1] = 0;
    if (bsize >= BLOCK_8X8)
      read_ext_intra_mode_info(cm, xd, r);
#endif  // CONFIG_EXT_INTRA
}

static int read_mv_component(vpx_reader *r,
                             const nmv_component *mvcomp, int usehp) {
  int mag, d, fr, hp;
  const int sign = vpx_read(r, mvcomp->sign);
  const int mv_class = vpx_read_tree(r, vp10_mv_class_tree, mvcomp->classes);
  const int class0 = mv_class == MV_CLASS_0;

  // Integer part
  if (class0) {
    d = vpx_read_tree(r, vp10_mv_class0_tree, mvcomp->class0);
    mag = 0;
  } else {
    int i;
    const int n = mv_class + CLASS0_BITS - 1;  // number of bits

    d = 0;
    for (i = 0; i < n; ++i)
      d |= vpx_read(r, mvcomp->bits[i]) << i;
    mag = CLASS0_SIZE << (mv_class + 2);
  }

  // Fractional part
  fr = vpx_read_tree(r, vp10_mv_fp_tree, class0 ? mvcomp->class0_fp[d]
                                               : mvcomp->fp);

  // High precision part (if hp is not used, the default value of the hp is 1)
  hp = usehp ? vpx_read(r, class0 ? mvcomp->class0_hp : mvcomp->hp)
             : 1;

  // Result
  mag += ((d << 3) | (fr << 1) | hp) + 1;
  return sign ? -mag : mag;
}

static INLINE void read_mv(vpx_reader *r, MV *mv, const MV *ref,
                           const nmv_context *ctx,
                           nmv_context_counts *counts, int allow_hp) {
  const MV_JOINT_TYPE joint_type =
      (MV_JOINT_TYPE)vpx_read_tree(r, vp10_mv_joint_tree, ctx->joints);
  const int use_hp = allow_hp && vp10_use_mv_hp(ref);
  MV diff = {0, 0};

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
                                                vpx_reader *r) {
  if (cm->reference_mode == REFERENCE_MODE_SELECT) {
    const int ctx = vp10_get_reference_mode_context(cm, xd);
    const REFERENCE_MODE mode =
        (REFERENCE_MODE)vpx_read(r, cm->fc->comp_inter_prob[ctx]);
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
                            vpx_reader *r,
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
      const int idx = cm->ref_frame_sign_bias[cm->comp_fixed_ref];
      const int ctx = vp10_get_pred_context_comp_ref_p(cm, xd);
      const int bit = vpx_read(r, fc->comp_ref_prob[ctx][0]);
      if (counts)
        ++counts->comp_ref[ctx][0][bit];
      ref_frame[idx] = cm->comp_fixed_ref;

#if CONFIG_EXT_REFS
      if (!bit) {
        const int ctx1 = vp10_get_pred_context_comp_ref_p1(cm, xd);
        const int bit1 = vpx_read(r, fc->comp_ref_prob[ctx1][1]);
        if (counts)
          ++counts->comp_ref[ctx1][1][bit1];
        ref_frame[!idx] = cm->comp_var_ref[bit1 ? 0 : 1];
      } else {
        const int ctx2 = vp10_get_pred_context_comp_ref_p2(cm, xd);
        const int bit2 = vpx_read(r, fc->comp_ref_prob[ctx2][2]);
        if (counts)
          ++counts->comp_ref[ctx2][2][bit2];
        if (!bit2) {
          const int ctx3 = vp10_get_pred_context_comp_ref_p3(cm, xd);
          const int bit3 = vpx_read(r, fc->comp_ref_prob[ctx3][3]);
          if (counts)
            ++counts->comp_ref[ctx3][3][bit3];
          ref_frame[!idx] = cm->comp_var_ref[bit3 ? 2 : 3];
        } else {
          ref_frame[!idx] = cm->comp_var_ref[4];
        }
      }
#else
      ref_frame[!idx] = cm->comp_var_ref[bit];
#endif  // CONFIG_EXT_REFS
    } else if (mode == SINGLE_REFERENCE) {
#if CONFIG_EXT_REFS
      const int ctx0 = vp10_get_pred_context_single_ref_p1(xd);
      const int bit0 = vpx_read(r, fc->single_ref_prob[ctx0][0]);
      if (counts)
        ++counts->single_ref[ctx0][0][bit0];
      if (bit0) {
        const int ctx1 = vp10_get_pred_context_single_ref_p2(xd);
        const int bit1 = vpx_read(r, fc->single_ref_prob[ctx1][1]);
        if (counts)
          ++counts->single_ref[ctx1][1][bit1];
        ref_frame[0] = bit1 ? ALTREF_FRAME : GOLDEN_FRAME;
      } else {
        const int ctx2 = vp10_get_pred_context_single_ref_p3(xd);
        const int bit2 = vpx_read(r, fc->single_ref_prob[ctx2][2]);
        if (counts)
          ++counts->single_ref[ctx2][2][bit2];
        if (bit2) {
          const int ctx4 = vp10_get_pred_context_single_ref_p5(xd);
          const int bit4 = vpx_read(r, fc->single_ref_prob[ctx4][4]);
          if (counts)
            ++counts->single_ref[ctx4][4][bit4];
          ref_frame[0] = bit4 ? LAST4_FRAME : LAST3_FRAME;
        } else {
          const int ctx3 = vp10_get_pred_context_single_ref_p4(xd);
          const int bit3 = vpx_read(r, fc->single_ref_prob[ctx3][3]);
          if (counts)
            ++counts->single_ref[ctx3][3][bit3];
          ref_frame[0] = bit3 ? LAST2_FRAME : LAST_FRAME;
        }
      }
#else
      const int ctx0 = vp10_get_pred_context_single_ref_p1(xd);
      const int bit0 = vpx_read(r, fc->single_ref_prob[ctx0][0]);
      if (counts)
        ++counts->single_ref[ctx0][0][bit0];
      if (bit0) {
        const int ctx1 = vp10_get_pred_context_single_ref_p2(xd);
        const int bit1 = vpx_read(r, fc->single_ref_prob[ctx1][1]);
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


#if CONFIG_OBMC
static int read_is_obmc_block(VP10_COMMON *const cm, MACROBLOCKD *const xd,
                              vpx_reader *r) {
  BLOCK_SIZE bsize = xd->mi[0]->mbmi.sb_type;
  FRAME_COUNTS *counts = xd->counts;
  int is_obmc;

  if (is_obmc_allowed(&xd->mi[0]->mbmi)) {
    is_obmc = vpx_read(r, cm->fc->obmc_prob[bsize]);
    if (counts)
      ++counts->obmc[bsize][is_obmc];
    return is_obmc;
  } else {
    return 0;
  }
}
#endif  // CONFIG_OBMC

static INLINE INTERP_FILTER read_switchable_interp_filter(
    VP10_COMMON *const cm, MACROBLOCKD *const xd,
    vpx_reader *r) {
  const int ctx = vp10_get_pred_context_switchable_interp(xd);
  FRAME_COUNTS *counts = xd->counts;
  INTERP_FILTER type;
#if CONFIG_EXT_INTERP
  if (!vp10_is_interp_needed(xd)) return EIGHTTAP;
#endif
  type = (INTERP_FILTER)vpx_read_tree(r, vp10_switchable_interp_tree,
                                      cm->fc->switchable_interp_prob[ctx]);
  if (counts)
    ++counts->switchable_interp[ctx][type];
  return type;
}

static void read_intra_block_mode_info(VP10_COMMON *const cm,
                                       MACROBLOCKD *const xd, MODE_INFO *mi,
                                       vpx_reader *r) {
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
#if CONFIG_EXT_INTRA
      mbmi->angle_delta[0] = 0;
      if (mbmi->mode != DC_PRED && mbmi->mode != TM_PRED) {
        int p_angle;
        mbmi->angle_delta[0] =
            read_uniform(r, 2 * MAX_ANGLE_DELTAS + 1) - MAX_ANGLE_DELTAS;
        p_angle =
            mode_to_angle_map[mbmi->mode] + mbmi->angle_delta[0] * ANGLE_STEP;
        if (pick_intra_filter(p_angle)) {
          FRAME_COUNTS *counts = xd->counts;
          const int ctx = vp10_get_pred_context_intra_interp(xd);
          mbmi->intra_filter = vpx_read_tree(r, vp10_intra_filter_tree,
                                             cm->fc->intra_filter_probs[ctx]);
          if (counts)
            ++counts->intra_filter[ctx][mbmi->intra_filter];
        } else {
          mbmi->intra_filter = INTRA_FILTER_LINEAR;
        }
      }
#endif  // CONFIG_EXT_INTRA
  }

  mbmi->uv_mode = read_intra_mode_uv(cm, xd, r, mbmi->mode);
#if CONFIG_EXT_INTRA
  if (mbmi->uv_mode != DC_PRED && mbmi->uv_mode != TM_PRED &&
      bsize >= BLOCK_8X8)
    mbmi->angle_delta[1] =
        read_uniform(r, 2 * MAX_ANGLE_DELTAS + 1) - MAX_ANGLE_DELTAS;
#endif  // CONFIG_EXT_INTRA

  mbmi->palette_mode_info.palette_size[0] = 0;
  mbmi->palette_mode_info.palette_size[1] = 0;
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
                            int is_compound, int allow_hp, vpx_reader *r) {
  int i;
  int ret = 1;
#if CONFIG_REF_MV
  MB_MODE_INFO *mbmi = &xd->mi[0]->mbmi;
  BLOCK_SIZE bsize = mbmi->sb_type;
  int_mv *pred_mv = (bsize >= BLOCK_8X8) ?
      mbmi->pred_mv : xd->mi[0]->bmi[block].pred_mv;
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
        read_mv(r, &mv[i].as_mv, &ref_mv[i].as_mv, &cm->fc->nmvc[nmv_ctx],
                mv_counts, allow_hp);
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
        read_mv(r, &mv[i].as_mv, &ref_mv[i].as_mv,
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
    case NEW_NEARESTMV: {
      FRAME_COUNTS *counts = xd->counts;
#if CONFIG_REF_MV
      int nmv_ctx = vp10_nmv_ctx(xd->ref_mv_count[mbmi->ref_frame[0]],
                                 xd->ref_mv_stack[mbmi->ref_frame[0]]);
      nmv_context_counts *const mv_counts =
          counts ? &counts->mv[nmv_ctx] : NULL;
      read_mv(r, &mv[0].as_mv, &ref_mv[0].as_mv,
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
      read_mv(r, &mv[1].as_mv, &ref_mv[1].as_mv,
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
      read_mv(r, &mv[1].as_mv, &ref_mv[1].as_mv,
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
      read_mv(r, &mv[0].as_mv, &ref_mv[0].as_mv,
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
                               int segment_id, vpx_reader *r) {
  if (segfeature_active(&cm->seg, segment_id, SEG_LVL_REF_FRAME)) {
    return get_segdata(&cm->seg, segment_id, SEG_LVL_REF_FRAME) != INTRA_FRAME;
  } else {
    const int ctx = vp10_get_intra_inter_context(xd);
    const int is_inter = vpx_read(r, cm->fc->intra_inter_prob[ctx]);
    FRAME_COUNTS *counts = xd->counts;
    if (counts)
      ++counts->intra_inter[ctx][is_inter];
    return is_inter;
  }
}

static void fpm_sync(void *const data, int mi_row) {
  VP10Decoder *const pbi = (VP10Decoder *)data;
  vp10_frameworker_wait(pbi->frame_worker_owner, pbi->common.prev_frame,
                       mi_row << MI_BLOCK_SIZE_LOG2);
}

static void read_inter_block_mode_info(VP10Decoder *const pbi,
                                       MACROBLOCKD *const xd,
                                       MODE_INFO *const mi,
#if (CONFIG_OBMC || CONFIG_EXT_INTER) && CONFIG_SUPERTX
                                       int mi_row, int mi_col, vpx_reader *r,
                                       int supertx_enabled) {
#else
                                       int mi_row, int mi_col, vpx_reader *r) {
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

#if CONFIG_OBMC
#if CONFIG_SUPERTX
  if (!supertx_enabled)
#endif  // CONFIG_SUPERTX
  mbmi->obmc = read_is_obmc_block(cm, xd, r);
#endif  // CONFIG_OBMC

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
      if (mbmi->mode == NEARMV)
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
    lower_mv_precision(&cur_mv.as_mv, cm->allow_high_precision_mv);
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
      int i;
#if CONFIG_EXT_INTER
      if (mbmi->mode == NEAREST_NEARESTMV) {
#endif  // CONFIG_EXT_INTER
      nearestmv[0] = xd->ref_mv_stack[ref_frame_type][0].this_mv;
      nearestmv[1] = xd->ref_mv_stack[ref_frame_type][0].comp_mv;

      for (i = 0; i < MAX_MV_REF_CANDIDATES; ++i)
        lower_mv_precision(&nearestmv[i].as_mv, allow_hp);
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
      if (mbmi->mode == NEAR_NEWMV || mbmi->mode == NEAR_NEARESTMV) {
        nearmv[0] = xd->ref_mv_stack[ref_frame_type][1].this_mv;
        lower_mv_precision(&nearmv[0].as_mv, allow_hp);
      }

      if (mbmi->mode == NEW_NEARMV || mbmi->mode == NEAREST_NEARMV) {
        nearmv[1] = xd->ref_mv_stack[ref_frame_type][1].comp_mv;
        lower_mv_precision(&nearmv[1].as_mv, allow_hp);
      }
    }
#else
    if (xd->ref_mv_count[ref_frame_type] > 1) {
      int i;
      int ref_mv_idx = 1 + mbmi->ref_mv_idx;
      nearestmv[0] = xd->ref_mv_stack[ref_frame_type][0].this_mv;
      nearestmv[1] = xd->ref_mv_stack[ref_frame_type][0].comp_mv;
      nearmv[0] = xd->ref_mv_stack[ref_frame_type][ref_mv_idx].this_mv;
      nearmv[1] = xd->ref_mv_stack[ref_frame_type][ref_mv_idx].comp_mv;

      for (i = 0; i < MAX_MV_REF_CANDIDATES; ++i) {
        lower_mv_precision(&nearestmv[i].as_mv, allow_hp);
        lower_mv_precision(&nearmv[i].as_mv, allow_hp);
      }
    }
#endif  // CONFIG_EXT_INTER
  }
#endif

#if !CONFIG_EXT_INTERP
  mbmi->interp_filter = (cm->interp_filter == SWITCHABLE)
                        ? read_switchable_interp_filter(cm, xd, r)
                        : cm->interp_filter;
#endif  // !CONFIG_EXT_INTERP

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
        if (b_mode == NEARESTMV || b_mode == NEARMV) {
#endif  // CONFIG_EXT_INTER
          for (ref = 0; ref < 1 + is_compound; ++ref)
#if CONFIG_EXT_INTER
          {
            int_mv mv_ref_list[MAX_MV_REF_CANDIDATES];
            vp10_update_mv_context(cm, xd, mi, mbmi->ref_frame[ref],
                                   mv_ref_list, j, mi_row, mi_col, NULL);
#endif  // CONFIG_EXT_INTER
            vp10_append_sub8x8_mvs_for_idx(cm, xd, j, ref, mi_row, mi_col,
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

        if (!assign_mv(cm, xd, b_mode,
#if CONFIG_REF_MV
                       j,
#endif
                       block,
#if CONFIG_EXT_INTER
                       ref_mv[mv_idx],
#else
                       nearestmv,
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
    mbmi->pred_mv[0].as_int = mi->bmi[3].pred_mv[0].as_int;
    mbmi->pred_mv[1].as_int = mi->bmi[3].pred_mv[1].as_int;
#endif
    mi->mbmi.mode = b_mode;

    mbmi->mv[0].as_int = mi->bmi[3].as_mv[0].as_int;
    mbmi->mv[1].as_int = mi->bmi[3].as_mv[1].as_int;
  } else {
    xd->corrupted |= !assign_mv(cm, xd, mbmi->mode,
#if CONFIG_REF_MV
                                0,
#endif
                                mbmi->mv,
#if CONFIG_EXT_INTER
                                mbmi->mode == NEWFROMNEARMV ?
                                              nearmv : nearestmv,
#else
                                nearestmv,
#endif  // CONFIG_EXT_INTER
                                nearestmv, nearmv, is_compound, allow_hp, r);
  }

#if CONFIG_EXT_INTER
  if (cm->reference_mode != COMPOUND_REFERENCE &&
#if CONFIG_SUPERTX
      !supertx_enabled &&
#endif
      is_interintra_allowed(mbmi)) {
    const int interintra = vpx_read(r, cm->fc->interintra_prob[bsize]);
    if (xd->counts)
      xd->counts->interintra[bsize][interintra]++;
    assert(mbmi->ref_frame[1] == NONE);
    if (interintra) {
      const PREDICTION_MODE interintra_mode =
          read_intra_mode_y(cm, xd, r, size_group_lookup[bsize]);

      mbmi->ref_frame[1] = INTRA_FRAME;
      mbmi->interintra_mode = interintra_mode;
      mbmi->interintra_uv_mode = interintra_mode;
#if CONFIG_EXT_INTRA
      // TODO(debargha|geza.lore):
      // Should we use ext_intra modes for interintra?
      mbmi->ext_intra_mode_info.use_ext_intra_mode[0] = 0;
      mbmi->ext_intra_mode_info.use_ext_intra_mode[1] = 0;
      mbmi->angle_delta[0] = 0;
      mbmi->angle_delta[1] = 0;
      mbmi->intra_filter = INTRA_FILTER_LINEAR;
#endif  // CONFIG_EXT_INTRA
    }
  }
#endif  // CONFIG_EXT_INTER

#if CONFIG_EXT_INTERP
  mbmi->interp_filter = (cm->interp_filter == SWITCHABLE)
                        ? read_switchable_interp_filter(cm, xd, r)
                        : cm->interp_filter;
#endif  // CONFIG_EXT_INTERP
}

static void read_inter_frame_mode_info(VP10Decoder *const pbi,
                                       MACROBLOCKD *const xd,
#if CONFIG_SUPERTX
                                       int supertx_enabled,
#endif  // CONFIG_SUPERTX
                                       int mi_row, int mi_col, vpx_reader *r) {
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
    xd->left_txfm_context = xd->left_txfm_context_buffer + (mi_row & MI_MASK);
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
          read_tx_size_inter(cm, xd, mbmi, xd->counts, max_tx_size,
                             idy, idx, r);
      if (xd->counts) {
        const int ctx = get_tx_size_context(xd);
        ++get_tx_counts(max_tx_size, ctx, &xd->counts->tx)[mbmi->tx_size];
      }
    } else {
      mbmi->tx_size = read_tx_size(cm, xd, !mbmi->skip || !inter_block, r);
      if (inter_block) {
        const int width  = num_4x4_blocks_wide_lookup[bsize];
        const int height = num_4x4_blocks_high_lookup[bsize];
        int idx, idy;
        for (idy = 0; idy < height; ++idy)
          for (idx = 0; idx < width; ++idx)
            mbmi->inter_tx_size[(idy >> 1) * 8 + (idx >> 1)] = mbmi->tx_size;
      }

      set_txfm_ctx(xd->left_txfm_context, mbmi->tx_size, xd->n8_h);
      set_txfm_ctx(xd->above_txfm_context, mbmi->tx_size, xd->n8_w);
    }
#else
    mbmi->tx_size = read_tx_size(cm, xd, !mbmi->skip || !inter_block, r);
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
        xd->mi[0]->mbmi.inter_tx_size[(idy >> 1) * 8 + (idx >> 1)] =
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
              vpx_read_tree(r, vp10_ext_tx_inter_tree[eset],
                            cm->fc->inter_ext_tx_prob[eset][mbmi->tx_size]);
          if (counts)
            ++counts->inter_ext_tx[eset][mbmi->tx_size][mbmi->tx_type];
        }
      } else if (ALLOW_INTRA_EXT_TX) {
        if (eset > 0) {
          mbmi->tx_type = vpx_read_tree(r, vp10_ext_tx_intra_tree[eset],
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
        mbmi->tx_type = vpx_read_tree(
            r, vp10_ext_tx_tree,
            cm->fc->inter_ext_tx_prob[mbmi->tx_size]);
        if (counts)
          ++counts->inter_ext_tx[mbmi->tx_size][mbmi->tx_type];
      } else {
        const TX_TYPE tx_type_nom = intra_mode_to_tx_type_context[mbmi->mode];
        mbmi->tx_type = vpx_read_tree(
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
                         int mi_row, int mi_col, vpx_reader *r,
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
