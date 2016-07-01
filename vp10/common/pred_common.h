/*
 *  Copyright (c) 2012 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#ifndef VP10_COMMON_PRED_COMMON_H_
#define VP10_COMMON_PRED_COMMON_H_

#include "vp10/common/blockd.h"
#include "vp10/common/onyxc_int.h"
#include "vpx_dsp/vpx_dsp_common.h"

#ifdef __cplusplus
extern "C" {
#endif

static INLINE int get_segment_id(const VP10_COMMON *cm,
                                 const uint8_t *segment_ids,
                                 BLOCK_SIZE bsize, int mi_row, int mi_col) {
  const int mi_offset = mi_row * cm->mi_cols + mi_col;
  const int bw = num_8x8_blocks_wide_lookup[bsize];
  const int bh = num_8x8_blocks_high_lookup[bsize];
  const int xmis = VPXMIN(cm->mi_cols - mi_col, bw);
  const int ymis = VPXMIN(cm->mi_rows - mi_row, bh);
  int x, y, segment_id = MAX_SEGMENTS;

  for (y = 0; y < ymis; ++y)
    for (x = 0; x < xmis; ++x)
      segment_id =
          VPXMIN(segment_id, segment_ids[mi_offset + y * cm->mi_cols + x]);

  assert(segment_id >= 0 && segment_id < MAX_SEGMENTS);
  return segment_id;
}

static INLINE int vp10_get_pred_context_seg_id(const MACROBLOCKD *xd) {
  const MODE_INFO *const above_mi = xd->above_mi;
  const MODE_INFO *const left_mi = xd->left_mi;
  const int above_sip = (above_mi != NULL) ?
                        above_mi->mbmi.seg_id_predicted : 0;
  const int left_sip = (left_mi != NULL) ? left_mi->mbmi.seg_id_predicted : 0;

  return above_sip + left_sip;
}

static INLINE vpx_prob vp10_get_pred_prob_seg_id(
    const struct segmentation_probs *segp, const MACROBLOCKD *xd) {
  return segp->pred_probs[vp10_get_pred_context_seg_id(xd)];
}

static INLINE int vp10_get_skip_context(const MACROBLOCKD *xd) {
  const MODE_INFO *const above_mi = xd->above_mi;
  const MODE_INFO *const left_mi = xd->left_mi;
  const int above_skip = (above_mi != NULL) ? above_mi->mbmi.skip : 0;
  const int left_skip = (left_mi != NULL) ? left_mi->mbmi.skip : 0;
  return above_skip + left_skip;
}

static INLINE vpx_prob vp10_get_skip_prob(const VP10_COMMON *cm,
                                         const MACROBLOCKD *xd) {
  return cm->fc->skip_probs[vp10_get_skip_context(xd)];
}

#if CONFIG_DUAL_FILTER
int vp10_get_pred_context_switchable_interp(const MACROBLOCKD *xd, int dir);
#else
int vp10_get_pred_context_switchable_interp(const MACROBLOCKD *xd);
#endif

#if CONFIG_EXT_INTRA
int vp10_get_pred_context_intra_interp(const MACROBLOCKD *xd);
#endif  // CONFIG_EXT_INTRA

int vp10_get_intra_inter_context(const MACROBLOCKD *xd);

static INLINE vpx_prob vp10_get_intra_inter_prob(const VP10_COMMON *cm,
                                                const MACROBLOCKD *xd) {
  return cm->fc->intra_inter_prob[vp10_get_intra_inter_context(xd)];
}

int vp10_get_reference_mode_context(const VP10_COMMON *cm,
                                    const MACROBLOCKD *xd);

static INLINE vpx_prob vp10_get_reference_mode_prob(const VP10_COMMON *cm,
                                                    const MACROBLOCKD *xd) {
  return cm->fc->comp_inter_prob[vp10_get_reference_mode_context(cm, xd)];
}

int vp10_get_pred_context_comp_ref_p(const VP10_COMMON *cm,
                                    const MACROBLOCKD *xd);

static INLINE vpx_prob vp10_get_pred_prob_comp_ref_p(const VP10_COMMON *cm,
                                                     const MACROBLOCKD *xd) {
  const int pred_context = vp10_get_pred_context_comp_ref_p(cm, xd);
  return cm->fc->comp_ref_prob[pred_context][0];
}

#if CONFIG_EXT_REFS
int vp10_get_pred_context_comp_ref_p1(const VP10_COMMON *cm,
                                      const MACROBLOCKD *xd);

static INLINE vpx_prob vp10_get_pred_prob_comp_ref_p1(const VP10_COMMON *cm,
                                                     const MACROBLOCKD *xd) {
  const int pred_context = vp10_get_pred_context_comp_ref_p1(cm, xd);
  return cm->fc->comp_ref_prob[pred_context][1];
}

int vp10_get_pred_context_comp_ref_p2(const VP10_COMMON *cm,
                                      const MACROBLOCKD *xd);

static INLINE vpx_prob vp10_get_pred_prob_comp_ref_p2(const VP10_COMMON *cm,
                                                     const MACROBLOCKD *xd) {
  const int pred_context = vp10_get_pred_context_comp_ref_p2(cm, xd);
  return cm->fc->comp_ref_prob[pred_context][2];
}

int vp10_get_pred_context_comp_bwdref_p(const VP10_COMMON *cm,
                                        const MACROBLOCKD *xd);

static INLINE vpx_prob vp10_get_pred_prob_comp_bwdref_p(const VP10_COMMON *cm,
                                                        const MACROBLOCKD *xd) {
  const int pred_context = vp10_get_pred_context_comp_bwdref_p(cm, xd);
  return cm->fc->comp_bwdref_prob[pred_context][0];
}

#endif  // CONFIG_EXT_REFS

int vp10_get_pred_context_single_ref_p1(const MACROBLOCKD *xd);

static INLINE vpx_prob vp10_get_pred_prob_single_ref_p1(const VP10_COMMON *cm,
                                                        const MACROBLOCKD *xd) {
  return cm->fc->single_ref_prob[vp10_get_pred_context_single_ref_p1(xd)][0];
}

int vp10_get_pred_context_single_ref_p2(const MACROBLOCKD *xd);

static INLINE vpx_prob vp10_get_pred_prob_single_ref_p2(const VP10_COMMON *cm,
                                                        const MACROBLOCKD *xd) {
  return cm->fc->single_ref_prob[vp10_get_pred_context_single_ref_p2(xd)][1];
}

#if CONFIG_EXT_REFS
int vp10_get_pred_context_single_ref_p3(const MACROBLOCKD *xd);

static INLINE vpx_prob vp10_get_pred_prob_single_ref_p3(const VP10_COMMON *cm,
                                                        const MACROBLOCKD *xd) {
  return cm->fc->single_ref_prob[vp10_get_pred_context_single_ref_p3(xd)][2];
}

int vp10_get_pred_context_single_ref_p4(const MACROBLOCKD *xd);

static INLINE vpx_prob vp10_get_pred_prob_single_ref_p4(const VP10_COMMON *cm,
                                                        const MACROBLOCKD *xd) {
  return cm->fc->single_ref_prob[vp10_get_pred_context_single_ref_p4(xd)][3];
}

int vp10_get_pred_context_single_ref_p5(const MACROBLOCKD *xd);

static INLINE vpx_prob vp10_get_pred_prob_single_ref_p5(const VP10_COMMON *cm,
                                                        const MACROBLOCKD *xd) {
  return cm->fc->single_ref_prob[vp10_get_pred_context_single_ref_p5(xd)][4];
}
#endif  // CONFIG_EXT_REFS

// Returns a context number for the given MB prediction signal
// The mode info data structure has a one element border above and to the
// left of the entries corresponding to real blocks.
// The prediction flags in these dummy entries are initialized to 0.
static INLINE int get_tx_size_context(const MACROBLOCKD *xd) {
  const int max_tx_size = max_txsize_lookup[xd->mi[0]->mbmi.sb_type];
  const MB_MODE_INFO *const above_mbmi = xd->above_mbmi;
  const MB_MODE_INFO *const left_mbmi = xd->left_mbmi;
  const int has_above = xd->up_available;
  const int has_left = xd->left_available;
  int above_ctx = (has_above && !above_mbmi->skip) ?
      (int)txsize_sqr_map[above_mbmi->tx_size] : max_tx_size;
  int left_ctx = (has_left && !left_mbmi->skip) ?
      (int)txsize_sqr_map[left_mbmi->tx_size] : max_tx_size;
  assert(xd->mi[0]->mbmi.sb_type >= BLOCK_8X8);
  if (!has_left)
    left_ctx = above_ctx;

  if (!has_above)
    above_ctx = left_ctx;

  return (above_ctx + left_ctx) > max_tx_size;
}

#if CONFIG_VAR_TX
static void update_tx_counts(VP10_COMMON *cm, MACROBLOCKD *xd,
                             MB_MODE_INFO *mbmi, BLOCK_SIZE plane_bsize,
                             TX_SIZE tx_size, int blk_row, int blk_col,
                             TX_SIZE max_tx_size, int ctx) {
  const struct macroblockd_plane *const pd = &xd->plane[0];
  const BLOCK_SIZE bsize = txsize_to_bsize[tx_size];
  const int tx_row = blk_row >> (1 - pd->subsampling_y);
  const int tx_col = blk_col >> (1 - pd->subsampling_x);
  const TX_SIZE plane_tx_size = mbmi->inter_tx_size[tx_row][tx_col];
  int max_blocks_high = num_4x4_blocks_high_lookup[plane_bsize];
  int max_blocks_wide = num_4x4_blocks_wide_lookup[plane_bsize];

  if (xd->mb_to_bottom_edge < 0)
    max_blocks_high += xd->mb_to_bottom_edge >> (5 + pd->subsampling_y);
  if (xd->mb_to_right_edge < 0)
    max_blocks_wide += xd->mb_to_right_edge >> (5 + pd->subsampling_x);

  if (blk_row >= max_blocks_high || blk_col >= max_blocks_wide)
    return;

  if (tx_size == plane_tx_size) {
    ++xd->counts->tx_size[max_tx_size - TX_8X8][ctx][tx_size];
    mbmi->tx_size = tx_size;
  } else {
    int bsl = b_width_log2_lookup[bsize];
    int i;

    assert(bsl > 0);
    --bsl;

    for (i = 0; i < 4; ++i) {
      const int offsetr = blk_row + ((i >> 1) << bsl);
      const int offsetc = blk_col + ((i & 0x01) << bsl);

      if (offsetr >= max_blocks_high || offsetc >= max_blocks_wide)
        continue;
      update_tx_counts(cm, xd, mbmi, plane_bsize,
                       tx_size - 1, offsetr, offsetc, max_tx_size, ctx);
    }
  }
}

static INLINE void inter_block_tx_count_update(VP10_COMMON *cm,
                                               MACROBLOCKD *xd,
                                               MB_MODE_INFO *mbmi,
                                               BLOCK_SIZE plane_bsize,
                                               int ctx) {
  const int mi_width = num_4x4_blocks_wide_lookup[plane_bsize];
  const int mi_height = num_4x4_blocks_high_lookup[plane_bsize];
  TX_SIZE max_tx_size = max_txsize_lookup[plane_bsize];
  BLOCK_SIZE txb_size = txsize_to_bsize[max_tx_size];
  int bh = num_4x4_blocks_wide_lookup[txb_size];
  int idx, idy;

  for (idy = 0; idy < mi_height; idy += bh)
    for (idx = 0; idx < mi_width; idx += bh)
      update_tx_counts(cm, xd, mbmi, plane_bsize, max_tx_size, idy, idx,
                       max_tx_size, ctx);
}
#endif

#ifdef __cplusplus
}  // extern "C"
#endif

#endif  // VP10_COMMON_PRED_COMMON_H_
