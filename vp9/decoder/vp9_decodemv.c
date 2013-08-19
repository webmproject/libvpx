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
#include "vp9/common/vp9_findnearmv.h"
#include "vp9/common/vp9_mvref_common.h"
#include "vp9/common/vp9_pred_common.h"
#include "vp9/common/vp9_reconinter.h"
#include "vp9/common/vp9_seg_common.h"

#include "vp9/decoder/vp9_decodemv.h"
#include "vp9/decoder/vp9_decodframe.h"
#include "vp9/decoder/vp9_onyxd_int.h"
#include "vp9/decoder/vp9_dsubexp.h"
#include "vp9/decoder/vp9_treereader.h"

static MB_PREDICTION_MODE read_intra_mode(vp9_reader *r, const vp9_prob *p) {
  return (MB_PREDICTION_MODE)treed_read(r, vp9_intra_mode_tree, p);
}

static MB_PREDICTION_MODE read_inter_mode(VP9_COMMON *cm, vp9_reader *r,
                                          uint8_t context) {
  MB_PREDICTION_MODE mode = treed_read(r, vp9_inter_mode_tree,
                            cm->fc.inter_mode_probs[context]);
  ++cm->counts.inter_mode[context][inter_mode_offset(mode)];
  return mode;
}

static int read_segment_id(vp9_reader *r, const struct segmentation *seg) {
  return treed_read(r, vp9_segment_tree, seg->tree_probs);
}

static TX_SIZE read_selected_tx_size(VP9_COMMON *cm, MACROBLOCKD *xd,
                                     BLOCK_SIZE_TYPE bsize, vp9_reader *r) {
  const uint8_t context = vp9_get_pred_context_tx_size(xd);
  const vp9_prob *tx_probs = get_tx_probs(bsize, context, &cm->fc.tx_probs);
  TX_SIZE tx_size = vp9_read(r, tx_probs[0]);
  if (tx_size != TX_4X4 && bsize >= BLOCK_16X16) {
    tx_size += vp9_read(r, tx_probs[1]);
    if (tx_size != TX_8X8 && bsize >= BLOCK_32X32)
      tx_size += vp9_read(r, tx_probs[2]);
  }

  update_tx_counts(bsize, context, tx_size, &cm->counts.tx);
  return tx_size;
}

static TX_SIZE read_tx_size(VP9D_COMP *pbi, TX_MODE tx_mode,
                            BLOCK_SIZE_TYPE bsize, int allow_select,
                            vp9_reader *r) {
  VP9_COMMON *const cm = &pbi->common;
  MACROBLOCKD *const xd = &pbi->mb;

  if (allow_select && tx_mode == TX_MODE_SELECT && bsize >= BLOCK_8X8)
    return read_selected_tx_size(cm, xd, bsize, r);
  else if (tx_mode >= ALLOW_32X32 && bsize >= BLOCK_32X32)
    return TX_32X32;
  else if (tx_mode >= ALLOW_16X16 && bsize >= BLOCK_16X16)
    return TX_16X16;
  else if (tx_mode >= ALLOW_8X8 && bsize >= BLOCK_8X8)
    return TX_8X8;
  else
    return TX_4X4;
}

static void set_segment_id(VP9_COMMON *cm, BLOCK_SIZE_TYPE bsize,
                           int mi_row, int mi_col, int segment_id) {
  const int mi_offset = mi_row * cm->mi_cols + mi_col;
  const int bw = 1 << mi_width_log2(bsize);
  const int bh = 1 << mi_height_log2(bsize);
  const int xmis = MIN(cm->mi_cols - mi_col, bw);
  const int ymis = MIN(cm->mi_rows - mi_row, bh);
  int x, y;

  assert(segment_id >= 0 && segment_id < MAX_SEGMENTS);

  for (y = 0; y < ymis; y++)
    for (x = 0; x < xmis; x++)
      cm->last_frame_seg_map[mi_offset + y * cm->mi_cols + x] = segment_id;
}

static int read_intra_segment_id(VP9D_COMP *pbi, int mi_row, int mi_col,
                                 vp9_reader *r) {
  MACROBLOCKD *const xd = &pbi->mb;
  struct segmentation *const seg = &xd->seg;
  const BLOCK_SIZE_TYPE bsize = xd->mode_info_context->mbmi.sb_type;
  int segment_id;

  if (!seg->enabled)
    return 0;  // Default for disabled segmentation

  if (!seg->update_map)
    return 0;

  segment_id = read_segment_id(r, seg);
  set_segment_id(&pbi->common, bsize, mi_row, mi_col, segment_id);
  return segment_id;
}

static int read_inter_segment_id(VP9D_COMP *pbi, int mi_row, int mi_col,
                                 vp9_reader *r) {
  VP9_COMMON *const cm = &pbi->common;
  MACROBLOCKD *const xd = &pbi->mb;
  struct segmentation *const seg = &xd->seg;
  const BLOCK_SIZE_TYPE bsize = xd->mode_info_context->mbmi.sb_type;
  int pred_segment_id, segment_id;

  if (!seg->enabled)
    return 0;  // Default for disabled segmentation

  pred_segment_id = vp9_get_segment_id(cm, cm->last_frame_seg_map,
                                       bsize, mi_row, mi_col);
  if (!seg->update_map)
    return pred_segment_id;

  if (seg->temporal_update) {
    const vp9_prob pred_prob = vp9_get_pred_prob_seg_id(xd);
    const int pred_flag = vp9_read(r, pred_prob);
    vp9_set_pred_flag_seg_id(cm, bsize, mi_row, mi_col, pred_flag);
    segment_id = pred_flag ? pred_segment_id
                           : read_segment_id(r, seg);
  } else {
    segment_id = read_segment_id(r, seg);
  }
  set_segment_id(cm, bsize, mi_row, mi_col, segment_id);
  return segment_id;
}

static uint8_t read_skip_coeff(VP9D_COMP *pbi, int segment_id, vp9_reader *r) {
  VP9_COMMON *const cm = &pbi->common;
  MACROBLOCKD *const xd = &pbi->mb;
  int skip_coeff = vp9_segfeature_active(&xd->seg, segment_id, SEG_LVL_SKIP);
  if (!skip_coeff) {
    const int ctx = vp9_get_pred_context_mbskip(xd);
    skip_coeff = vp9_read(r, vp9_get_pred_prob_mbskip(cm, xd));
    cm->counts.mbskip[ctx][skip_coeff]++;
  }
  return skip_coeff;
}

static void read_intra_frame_mode_info(VP9D_COMP *pbi, MODE_INFO *m,
                                       int mi_row, int mi_col, vp9_reader *r) {
  VP9_COMMON *const cm = &pbi->common;
  MACROBLOCKD *const xd = &pbi->mb;
  MB_MODE_INFO *const mbmi = &m->mbmi;
  const BLOCK_SIZE_TYPE bsize = mbmi->sb_type;
  const int mis = cm->mode_info_stride;

  mbmi->segment_id = read_intra_segment_id(pbi, mi_row, mi_col, r);
  mbmi->mb_skip_coeff = read_skip_coeff(pbi, mbmi->segment_id, r);
  mbmi->txfm_size = read_tx_size(pbi, cm->tx_mode, bsize, 1, r);
  mbmi->ref_frame[0] = INTRA_FRAME;
  mbmi->ref_frame[1] = NONE;

  if (bsize >= BLOCK_8X8) {
    const MB_PREDICTION_MODE A = above_block_mode(m, 0, mis);
    const MB_PREDICTION_MODE L = xd->left_available ?
                                  left_block_mode(m, 0) : DC_PRED;
    mbmi->mode = read_intra_mode(r, vp9_kf_y_mode_prob[A][L]);
#if CONFIG_FILTERINTRA
    if ((mbmi->txfm_size <= TX_8X8) && is_filter_allowed(mbmi->mode))
      mbmi->filterbit = vp9_read(r,
                         cm->fc.filterintra_prob[mbmi->txfm_size][mbmi->mode]);
    else
      mbmi->filterbit = 0;
#endif
  } else {
    // Only 4x4, 4x8, 8x4 blocks
    const int num_4x4_w = num_4x4_blocks_wide_lookup[bsize];  // 1 or 2
    const int num_4x4_h = num_4x4_blocks_high_lookup[bsize];  // 1 or 2
    int idx, idy;

    for (idy = 0; idy < 2; idy += num_4x4_h) {
      for (idx = 0; idx < 2; idx += num_4x4_w) {
        const int ib = idy * 2 + idx;
        const MB_PREDICTION_MODE A = above_block_mode(m, ib, mis);
        const MB_PREDICTION_MODE L = (xd->left_available || idx) ?
                                      left_block_mode(m, ib) : DC_PRED;
        const MB_PREDICTION_MODE b_mode = read_intra_mode(r,
                                              vp9_kf_y_mode_prob[A][L]);
        m->bmi[ib].as_mode = b_mode;
        if (num_4x4_h == 2)
          m->bmi[ib + 2].as_mode = b_mode;
        if (num_4x4_w == 2)
          m->bmi[ib + 1].as_mode = b_mode;
#if CONFIG_FILTERINTRA
        if (is_filter_allowed(b_mode))
          m->b_filter_info[ib] = vp9_read(r,
                                          cm->fc.filterintra_prob[0][b_mode]);
        else
          m->b_filter_info[ib] = 0;
        if (num_4x4_h == 2)
          m->b_filter_info[ib + 2] = m->b_filter_info[ib];
        if (num_4x4_w == 2)
          m->b_filter_info[ib + 1] = m->b_filter_info[ib];
#endif
      }
    }

    mbmi->mode = m->bmi[3].as_mode;
#if CONFIG_FILTERINTRA
    mbmi->filterbit = m->b_filter_info[3];
#endif
  }

  mbmi->uv_mode = read_intra_mode(r, vp9_kf_uv_mode_prob[mbmi->mode]);
#if CONFIG_FILTERINTRA
  if ((get_uv_tx_size(mbmi) <= TX_8X8) && is_filter_allowed(mbmi->uv_mode))
    mbmi->uv_filterbit = vp9_read(r,
                  cm->fc.filterintra_prob[get_uv_tx_size(mbmi)][mbmi->uv_mode]);
  else
    mbmi->uv_filterbit = 0;
#endif
}

static int read_mv_component(vp9_reader *r,
                             const nmv_component *mvcomp, int usehp) {

  int mag, d, fr, hp;
  const int sign = vp9_read(r, mvcomp->sign);
  const int mv_class = treed_read(r, vp9_mv_class_tree, mvcomp->classes);
  const int class0 = mv_class == MV_CLASS_0;

  // Integer part
  if (class0) {
    d = treed_read(r, vp9_mv_class0_tree, mvcomp->class0);
  } else {
    int i;
    const int n = mv_class + CLASS0_BITS - 1;  // number of bits

    d = 0;
    for (i = 0; i < n; ++i)
      d |= vp9_read(r, mvcomp->bits[i]) << i;
  }

  // Fractional part
  fr = treed_read(r, vp9_mv_fp_tree,
                  class0 ? mvcomp->class0_fp[d] : mvcomp->fp);


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
  const MV_JOINT_TYPE j = treed_read(r, vp9_mv_joint_tree, ctx->joints);
  const int use_hp = allow_hp && vp9_use_mv_hp(ref);
  MV diff = {0, 0};

  if (mv_joint_vertical(j))
    diff.row = read_mv_component(r, &ctx->comps[0], use_hp);

  if (mv_joint_horizontal(j))
    diff.col = read_mv_component(r, &ctx->comps[1], use_hp);

  vp9_inc_mv(&diff, counts);

  mv->row = ref->row + diff.row;
  mv->col = ref->col + diff.col;
}

static void update_mv(vp9_reader *r, vp9_prob *p) {
  if (vp9_read(r, VP9_NMV_UPDATE_PROB))
    *p = (vp9_read_literal(r, 7) << 1) | 1;
}

static void read_mv_probs(vp9_reader *r, nmv_context *mvc, int allow_hp) {
  int i, j, k;

  for (j = 0; j < MV_JOINTS - 1; ++j)
    update_mv(r, &mvc->joints[j]);

  for (i = 0; i < 2; ++i) {
    nmv_component *const comp = &mvc->comps[i];

    update_mv(r, &comp->sign);

    for (j = 0; j < MV_CLASSES - 1; ++j)
      update_mv(r, &comp->classes[j]);

    for (j = 0; j < CLASS0_SIZE - 1; ++j)
      update_mv(r, &comp->class0[j]);

    for (j = 0; j < MV_OFFSET_BITS; ++j)
      update_mv(r, &comp->bits[j]);
  }

  for (i = 0; i < 2; ++i) {
    nmv_component *const comp = &mvc->comps[i];

    for (j = 0; j < CLASS0_SIZE; ++j)
      for (k = 0; k < 3; ++k)
        update_mv(r, &comp->class0_fp[j][k]);

    for (j = 0; j < 3; ++j)
      update_mv(r, &comp->fp[j]);
  }

  if (allow_hp) {
    for (i = 0; i < 2; ++i) {
      update_mv(r, &mvc->comps[i].class0_hp);
      update_mv(r, &mvc->comps[i].hp);
    }
  }
}

// Read the referncence frame
static void read_ref_frames(VP9D_COMP *pbi, vp9_reader *r,
                            int segment_id, MV_REFERENCE_FRAME ref_frame[2]) {
  VP9_COMMON *const cm = &pbi->common;
  MACROBLOCKD *const xd = &pbi->mb;
  FRAME_CONTEXT *const fc = &cm->fc;
  FRAME_COUNTS *const counts = &cm->counts;

  if (vp9_segfeature_active(&xd->seg, segment_id, SEG_LVL_REF_FRAME)) {
    ref_frame[0] = vp9_get_segdata(&xd->seg, segment_id, SEG_LVL_REF_FRAME);
    ref_frame[1] = NONE;
  } else {
    const int comp_ctx = vp9_get_pred_context_comp_inter_inter(cm, xd);
    int is_comp;

    if (cm->comp_pred_mode == HYBRID_PREDICTION) {
      is_comp = vp9_read(r, fc->comp_inter_prob[comp_ctx]);
      counts->comp_inter[comp_ctx][is_comp]++;
    } else {
      is_comp = cm->comp_pred_mode == COMP_PREDICTION_ONLY;
    }

    // FIXME(rbultje) I'm pretty sure this breaks segmentation ref frame coding
    if (is_comp) {
      const int fix_ref_idx = cm->ref_frame_sign_bias[cm->comp_fixed_ref];
      const int ref_ctx = vp9_get_pred_context_comp_ref_p(cm, xd);
      const int b = vp9_read(r, fc->comp_ref_prob[ref_ctx]);
      counts->comp_ref[ref_ctx][b]++;
      ref_frame[fix_ref_idx] = cm->comp_fixed_ref;
      ref_frame[!fix_ref_idx] = cm->comp_var_ref[b];
    } else {
      const int ctx0 = vp9_get_pred_context_single_ref_p1(xd);
      const int bit0 = vp9_read(r, fc->single_ref_prob[ctx0][0]);
      ++counts->single_ref[ctx0][0][bit0];
      if (bit0) {
        const int ctx1 = vp9_get_pred_context_single_ref_p2(xd);
        const int bit1 = vp9_read(r, fc->single_ref_prob[ctx1][1]);
        ref_frame[0] = bit1 ? ALTREF_FRAME : GOLDEN_FRAME;
        ++counts->single_ref[ctx1][1][bit1];
      } else {
        ref_frame[0] = LAST_FRAME;
      }

      ref_frame[1] = NONE;
    }
  }
}

static void read_switchable_interp_probs(FRAME_CONTEXT *fc, vp9_reader *r) {
  int i, j;
  for (j = 0; j < VP9_SWITCHABLE_FILTERS + 1; ++j)
    for (i = 0; i < VP9_SWITCHABLE_FILTERS - 1; ++i)
      if (vp9_read(r, VP9_MODE_UPDATE_PROB))
        vp9_diff_update_prob(r, &fc->switchable_interp_prob[j][i]);
}

static void read_inter_mode_probs(FRAME_CONTEXT *fc, vp9_reader *r) {
  int i, j;
  for (i = 0; i < INTER_MODE_CONTEXTS; ++i)
    for (j = 0; j < VP9_INTER_MODES - 1; ++j)
      if (vp9_read(r, VP9_MODE_UPDATE_PROB))
        vp9_diff_update_prob(r, &fc->inter_mode_probs[i][j]);
}

static INLINE COMPPREDMODE_TYPE read_comp_pred_mode(vp9_reader *r) {
  COMPPREDMODE_TYPE mode = vp9_read_bit(r);
  if (mode)
    mode += vp9_read_bit(r);
  return mode;
}

static INLINE INTERPOLATIONFILTERTYPE read_switchable_filter_type(
    VP9D_COMP *pbi, vp9_reader *r) {
  VP9_COMMON *const cm = &pbi->common;
  MACROBLOCKD *const xd = &pbi->mb;
  const int ctx = vp9_get_pred_context_switchable_interp(xd);
  const int type = treed_read(r, vp9_switchable_interp_tree,
                              cm->fc.switchable_interp_prob[ctx]);
  ++cm->counts.switchable_interp[ctx][type];
  return type;
}

static void read_intra_block_mode_info(VP9D_COMP *pbi, MODE_INFO *mi,
                                  vp9_reader *r) {
  VP9_COMMON *const cm = &pbi->common;
  MB_MODE_INFO *const mbmi = &mi->mbmi;
  const BLOCK_SIZE_TYPE bsize = mi->mbmi.sb_type;

  mbmi->ref_frame[0] = INTRA_FRAME;
  mbmi->ref_frame[1] = NONE;

  if (bsize >= BLOCK_8X8) {
    const int size_group = size_group_lookup[bsize];
    mbmi->mode = read_intra_mode(r, cm->fc.y_mode_prob[size_group]);
    cm->counts.y_mode[size_group][mbmi->mode]++;
#if CONFIG_FILTERINTRA
    if (is_filter_allowed(mbmi->mode) && (mbmi->txfm_size <= TX_8X8)) {
      mbmi->filterbit = vp9_read(r,
                          cm->fc.filterintra_prob[mbmi->txfm_size][mbmi->mode]);
      cm->counts.filterintra[mbmi->txfm_size][mbmi->mode][mbmi->filterbit]++;
    } else {
      mbmi->filterbit = 0;
    }
#endif
  } else {
     // Only 4x4, 4x8, 8x4 blocks
     const int num_4x4_w = num_4x4_blocks_wide_lookup[bsize];  // 1 or 2
     const int num_4x4_h = num_4x4_blocks_high_lookup[bsize];  // 1 or 2
     int idx, idy;

     for (idy = 0; idy < 2; idy += num_4x4_h) {
       for (idx = 0; idx < 2; idx += num_4x4_w) {
         const int ib = idy * 2 + idx;
         const int b_mode = read_intra_mode(r, cm->fc.y_mode_prob[0]);
         mi->bmi[ib].as_mode = b_mode;
         cm->counts.y_mode[0][b_mode]++;

         if (num_4x4_h == 2)
           mi->bmi[ib + 2].as_mode = b_mode;
         if (num_4x4_w == 2)
           mi->bmi[ib + 1].as_mode = b_mode;
#if CONFIG_FILTERINTRA
         if (is_filter_allowed(b_mode)) {
           mi->b_filter_info[ib] = vp9_read(r,
                                            cm->fc.filterintra_prob[0][b_mode]);
           cm->counts.filterintra[0][b_mode][mi->b_filter_info[ib]]++;
         } else {
           mi->b_filter_info[ib] = 0;
         }

         if (num_4x4_h == 2)
           mi->b_filter_info[ib + 2] = mi->b_filter_info[ib];
         if (num_4x4_w == 2)
           mi->b_filter_info[ib + 1] = mi->b_filter_info[ib];
#endif
      }
    }
    mbmi->mode = mi->bmi[3].as_mode;
#if CONFIG_FILTERINTRA
    mbmi->filterbit = mi->b_filter_info[3];
#endif
  }

  mbmi->uv_mode = read_intra_mode(r, cm->fc.uv_mode_prob[mbmi->mode]);
  cm->counts.uv_mode[mbmi->mode][mbmi->uv_mode]++;
#if CONFIG_FILTERINTRA
  if (is_filter_allowed(mbmi->uv_mode) && (get_uv_tx_size(mbmi) <= TX_8X8)) {
    mbmi->uv_filterbit = vp9_read(r,
                  cm->fc.filterintra_prob[get_uv_tx_size(mbmi)][mbmi->uv_mode]);
    cm->counts.filterintra[get_uv_tx_size(mbmi)]
                          [mbmi->uv_mode][mbmi->uv_filterbit]++;
  } else {
    mbmi->uv_filterbit = 0;
  }
#endif
}

static int read_is_inter_block(VP9D_COMP *pbi, int segment_id, vp9_reader *r) {
  VP9_COMMON *const cm = &pbi->common;
  MACROBLOCKD *const xd = &pbi->mb;

  if (vp9_segfeature_active(&xd->seg, segment_id, SEG_LVL_REF_FRAME)) {
    return vp9_get_segdata(&xd->seg, segment_id, SEG_LVL_REF_FRAME) !=
           INTRA_FRAME;
  } else {
    const int ctx = vp9_get_pred_context_intra_inter(xd);
    const int is_inter = vp9_read(r, vp9_get_pred_prob_intra_inter(cm, xd));
    ++cm->counts.intra_inter[ctx][is_inter];
    return is_inter;
  }
}

static void read_inter_block_mode_info(VP9D_COMP *pbi, MODE_INFO *mi,
                                       int mi_row, int mi_col, vp9_reader *r) {
  VP9_COMMON *const cm = &pbi->common;
  MACROBLOCKD *const xd = &pbi->mb;
  nmv_context *const nmvc = &cm->fc.nmvc;
  MB_MODE_INFO *const mbmi = &mi->mbmi;
  int_mv *const mv0 = &mbmi->mv[0];
  int_mv *const mv1 = &mbmi->mv[1];
  const BLOCK_SIZE_TYPE bsize = mbmi->sb_type;
  const int allow_hp = xd->allow_high_precision_mv;

  int_mv nearest, nearby, best_mv;
  int_mv nearest_second, nearby_second, best_mv_second;
  uint8_t inter_mode_ctx;
  MV_REFERENCE_FRAME ref0, ref1;

#if CONFIG_MASKED_COMPOUND
  mbmi->use_masked_compound = 0;
  mbmi->mask_index = MASK_NONE;
#endif

  read_ref_frames(pbi, r, mbmi->segment_id, mbmi->ref_frame);
  ref0 = mbmi->ref_frame[0];
  ref1 = mbmi->ref_frame[1];

  vp9_find_mv_refs(cm, xd, mi, xd->prev_mode_info_context,
                   ref0, mbmi->ref_mvs[ref0], cm->ref_frame_sign_bias,
                   mi_row, mi_col);

  inter_mode_ctx = mbmi->mb_mode_context[ref0];

  if (vp9_segfeature_active(&xd->seg, mbmi->segment_id, SEG_LVL_SKIP))
    mbmi->mode = ZEROMV;
  else if (bsize >= BLOCK_8X8)
    mbmi->mode = read_inter_mode(cm, r, inter_mode_ctx);

  mbmi->uv_mode = DC_PRED;

  // nearest, nearby
  if (bsize < BLOCK_8X8 || mbmi->mode != ZEROMV) {
    vp9_find_best_ref_mvs(xd, mbmi->ref_mvs[ref0], &nearest, &nearby);
    best_mv.as_int = mbmi->ref_mvs[ref0][0].as_int;
  }

  mbmi->interp_filter = cm->mcomp_filter_type == SWITCHABLE
                            ? read_switchable_filter_type(pbi, r)
                            : cm->mcomp_filter_type;

  if (ref1 > INTRA_FRAME) {
    vp9_find_mv_refs(cm, xd, mi, xd->prev_mode_info_context,
                     ref1, mbmi->ref_mvs[ref1], cm->ref_frame_sign_bias,
                     mi_row, mi_col);

    if (bsize < BLOCK_8X8 || mbmi->mode != ZEROMV) {
      vp9_find_best_ref_mvs(xd, mbmi->ref_mvs[ref1],
                            &nearest_second, &nearby_second);
      best_mv_second.as_int = mbmi->ref_mvs[ref1][0].as_int;
    }
  }

#if CONFIG_INTERINTRA
    if ((cm->use_interintra)
        && is_interintra_allowed(bsize)
        && is_inter_mode(mbmi->mode)
        && (mbmi->ref_frame[1] == NONE)
        ) {
      mbmi->ref_frame[1] = (vp9_read(r, cm->fc.interintra_prob[bsize]) ?
                            INTRA_FRAME : NONE);
      cm->counts.interintra[bsize][mbmi->ref_frame[1] == INTRA_FRAME]++;
#if CONFIG_MASKED_COMPOUND
      mbmi->use_masked_interintra = 0;
#endif
      if (mbmi->ref_frame[1] == INTRA_FRAME) {
        int bsg = MIN(MIN(b_width_log2(bsize), b_height_log2(bsize)), 3);
        mbmi->interintra_mode = read_intra_mode(r, cm->fc.y_mode_prob[bsg]);
        cm->counts.y_mode[bsg][mbmi->interintra_mode]++;
#if SEPARATE_INTERINTRA_UV
        mbmi->interintra_uv_mode = read_intra_mode(r,
                                   cm->fc.uv_mode_prob[mbmi->interintra_mode]);
        cm->counts.uv_mode[mbmi->interintra_mode][mbmi->interintra_uv_mode]++;
#else
        mbmi->interintra_uv_mode = mbmi->interintra_mode;
#endif
#if CONFIG_MASKED_COMPOUND
        if (cm->use_masked_interintra && get_mask_bits_interintra(bsize)) {
          mbmi->use_masked_interintra = vp9_read(r,
                                          cm->fc.masked_interintra_prob[bsize]);
          cm->counts.masked_interintra[bsize][mbmi->use_masked_interintra]++;
          if (mbmi->use_masked_interintra) {
            mbmi->interintra_mask_index = vp9_read_literal(r,
                                               get_mask_bits_interintra(bsize));
            mbmi->interintra_uv_mask_index = mbmi->interintra_mask_index;
          }
        }
#endif
      }
    }
#endif

  if (bsize < BLOCK_8X8) {
    const int num_4x4_w = num_4x4_blocks_wide_lookup[bsize];  // 1 or 2
    const int num_4x4_h = num_4x4_blocks_high_lookup[bsize];  // 1 or 2
    int idx, idy;
    for (idy = 0; idy < 2; idy += num_4x4_h) {
      for (idx = 0; idx < 2; idx += num_4x4_w) {
        int_mv blockmv, secondmv;
        const int j = idy * 2 + idx;
        const int b_mode = read_inter_mode(cm, r, inter_mode_ctx);

        if (b_mode == NEARESTMV || b_mode == NEARMV) {
          vp9_append_sub8x8_mvs_for_idx(cm, xd, &nearest, &nearby, j, 0,
                                        mi_row, mi_col);

          if (ref1 > 0)
            vp9_append_sub8x8_mvs_for_idx(cm, xd,  &nearest_second,
                                         &nearby_second, j, 1,
                                         mi_row, mi_col);
        }

        switch (b_mode) {
          case NEWMV:
            read_mv(r, &blockmv.as_mv, &best_mv.as_mv, nmvc,
                    &cm->counts.mv, allow_hp);

            if (ref1 > 0)
              read_mv(r, &secondmv.as_mv, &best_mv_second.as_mv, nmvc,
                      &cm->counts.mv, allow_hp);
            break;
          case NEARESTMV:
            blockmv.as_int = nearest.as_int;
            if (ref1 > 0)
              secondmv.as_int = nearest_second.as_int;
            break;
          case NEARMV:
            blockmv.as_int = nearby.as_int;
            if (ref1 > 0)
              secondmv.as_int = nearby_second.as_int;
            break;
          case ZEROMV:
            blockmv.as_int = 0;
            if (ref1 > 0)
              secondmv.as_int = 0;
            break;
          default:
            assert(!"Invalid inter mode value");
        }
        mi->bmi[j].as_mv[0].as_int = blockmv.as_int;
        if (ref1 > 0)
          mi->bmi[j].as_mv[1].as_int = secondmv.as_int;

        if (num_4x4_h == 2)
          mi->bmi[j + 2] = mi->bmi[j];
        if (num_4x4_w == 2)
          mi->bmi[j + 1] = mi->bmi[j];
        mi->mbmi.mode = b_mode;
      }
    }

    mv0->as_int = mi->bmi[3].as_mv[0].as_int;
    mv1->as_int = mi->bmi[3].as_mv[1].as_int;
  } else {
    switch (mbmi->mode) {
      case NEARMV:
        mv0->as_int = nearby.as_int;
        clamp_mv2(&mv0->as_mv, xd);

        if (ref1 > 0) {
          mv1->as_int = nearby_second.as_int;
          clamp_mv2(&mv1->as_mv, xd);
        }
        break;

      case NEARESTMV:
        mv0->as_int = nearest.as_int;
        clamp_mv2(&mv0->as_mv, xd);

        if (ref1 > 0) {
          mv1->as_int = nearest_second.as_int;
          clamp_mv2(&mv1->as_mv, xd);
        }
        break;

      case ZEROMV:
        mv0->as_int = 0;
        if (ref1 > 0)
          mv1->as_int = 0;
        break;

      case NEWMV:
        read_mv(r, &mv0->as_mv, &best_mv.as_mv, nmvc, &cm->counts.mv, allow_hp);
        if (ref1 > 0)
          read_mv(r, &mv1->as_mv, &best_mv_second.as_mv, nmvc, &cm->counts.mv,
                  allow_hp);
        break;
      default:
        assert(!"Invalid inter mode value");
    }
  }
#if CONFIG_MASKED_COMPOUND
    mbmi->use_masked_compound = 0;
    if (pbi->common.use_masked_compound &&
        pbi->common.comp_pred_mode != SINGLE_PREDICTION_ONLY &&
        is_inter_mode(mbmi->mode) &&
        get_mask_bits(mi->mbmi.sb_type) &&
        mbmi->ref_frame[1] > INTRA_FRAME) {
      mbmi->use_masked_compound =
          vp9_read(r, pbi->common.fc.masked_compound_prob[bsize]);
      pbi->common.counts.masked_compound[bsize][mbmi->use_masked_compound]++;
      if (mbmi->use_masked_compound) {
        mbmi->mask_index = vp9_read_literal(r, get_mask_bits(mi->mbmi.sb_type));
      }
    }
#endif
}

static void read_inter_frame_mode_info(VP9D_COMP *pbi, MODE_INFO *mi,
                                       int mi_row, int mi_col, vp9_reader *r) {
  VP9_COMMON *const cm = &pbi->common;
  MB_MODE_INFO *const mbmi = &mi->mbmi;
  int inter_block;

  mbmi->mv[0].as_int = 0;
  mbmi->mv[1].as_int = 0;
  mbmi->segment_id = read_inter_segment_id(pbi, mi_row, mi_col, r);
  mbmi->mb_skip_coeff = read_skip_coeff(pbi, mbmi->segment_id, r);
  inter_block = read_is_inter_block(pbi, mbmi->segment_id, r);
  mbmi->txfm_size = read_tx_size(pbi, cm->tx_mode, mbmi->sb_type,
                                 !mbmi->mb_skip_coeff || !inter_block, r);

  if (inter_block)
    read_inter_block_mode_info(pbi, mi, mi_row, mi_col, r);
  else
    read_intra_block_mode_info(pbi, mi, r);
}

static void read_comp_pred(VP9_COMMON *cm, vp9_reader *r) {
  int i;

  cm->comp_pred_mode = cm->allow_comp_inter_inter ? read_comp_pred_mode(r)
                                                  : SINGLE_PREDICTION_ONLY;

  if (cm->comp_pred_mode == HYBRID_PREDICTION)
    for (i = 0; i < COMP_INTER_CONTEXTS; i++)
      if (vp9_read(r, VP9_MODE_UPDATE_PROB))
        vp9_diff_update_prob(r, &cm->fc.comp_inter_prob[i]);

#if CONFIG_MASKED_COMPOUND
  if (cm->comp_pred_mode != SINGLE_PREDICTION_ONLY) {
    cm->use_masked_compound = vp9_read_bit(r);
    if (cm->use_masked_compound) {
      for (i = 0; i < BLOCK_SIZE_TYPES; ++i) {
        if (get_mask_bits(i))
          if (vp9_read(r, VP9_UPD_MASKED_COMPOUND_PROB))
            vp9_diff_update_prob(r, &cm->fc.masked_compound_prob[i]);
      }
    }
  } else {
    cm->use_masked_compound = 0;
  }
#endif

  if (cm->comp_pred_mode != COMP_PREDICTION_ONLY)
    for (i = 0; i < REF_CONTEXTS; i++) {
      if (vp9_read(r, VP9_MODE_UPDATE_PROB))
        vp9_diff_update_prob(r, &cm->fc.single_ref_prob[i][0]);
      if (vp9_read(r, VP9_MODE_UPDATE_PROB))
        vp9_diff_update_prob(r, &cm->fc.single_ref_prob[i][1]);
    }

  if (cm->comp_pred_mode != SINGLE_PREDICTION_ONLY)
    for (i = 0; i < REF_CONTEXTS; i++)
      if (vp9_read(r, VP9_MODE_UPDATE_PROB))
        vp9_diff_update_prob(r, &cm->fc.comp_ref_prob[i]);
}

void vp9_prepare_read_mode_info(VP9D_COMP* pbi, vp9_reader *r) {
  VP9_COMMON *const cm = &pbi->common;
  int k;

  // TODO(jkoleszar): does this clear more than MBSKIP_CONTEXTS? Maybe remove.
  // vpx_memset(cm->fc.mbskip_probs, 0, sizeof(cm->fc.mbskip_probs));
  for (k = 0; k < MBSKIP_CONTEXTS; ++k)
    if (vp9_read(r, VP9_MODE_UPDATE_PROB))
      vp9_diff_update_prob(r, &cm->fc.mbskip_probs[k]);

  if (cm->frame_type != KEY_FRAME && !cm->intra_only) {
    nmv_context *const nmvc = &pbi->common.fc.nmvc;
    MACROBLOCKD *const xd = &pbi->mb;
    int i, j;

    read_inter_mode_probs(&cm->fc, r);

    if (cm->mcomp_filter_type == SWITCHABLE)
      read_switchable_interp_probs(&cm->fc, r);

#if CONFIG_INTERINTRA
    if (cm->use_interintra) {
      int b;
      for (b = 0; b < BLOCK_SIZE_TYPES; ++b) {
        if (is_interintra_allowed(b))
          if (vp9_read(r, VP9_UPD_INTERINTRA_PROB))
            vp9_diff_update_prob(r, &cm->fc.interintra_prob[b]);
      }
#if CONFIG_MASKED_COMPOUND
      if (cm->use_masked_interintra) {
        int k;
        for (k = 0; k < BLOCK_SIZE_TYPES; ++k) {
          if (is_interintra_allowed(k) && get_mask_bits_interintra(k))
            if (vp9_read(r, VP9_UPD_MASKED_INTERINTRA_PROB))
              vp9_diff_update_prob(r, &cm->fc.masked_interintra_prob[k]);
        }
      }
#endif
    }
#endif

    for (i = 0; i < INTRA_INTER_CONTEXTS; i++)
      if (vp9_read(r, VP9_MODE_UPDATE_PROB))
        vp9_diff_update_prob(r, &cm->fc.intra_inter_prob[i]);

    read_comp_pred(cm, r);

    for (j = 0; j < BLOCK_SIZE_GROUPS; j++)
      for (i = 0; i < VP9_INTRA_MODES - 1; ++i)
        if (vp9_read(r, VP9_MODE_UPDATE_PROB))
          vp9_diff_update_prob(r, &cm->fc.y_mode_prob[j][i]);

    for (j = 0; j < NUM_PARTITION_CONTEXTS; ++j)
      for (i = 0; i < PARTITION_TYPES - 1; ++i)
        if (vp9_read(r, VP9_MODE_UPDATE_PROB))
          vp9_diff_update_prob(r, &cm->fc.partition_prob[INTER_FRAME][j][i]);

    read_mv_probs(r, nmvc, xd->allow_high_precision_mv);
  }
}

void vp9_read_mode_info(VP9D_COMP* pbi, int mi_row, int mi_col, vp9_reader *r) {
  VP9_COMMON *const cm = &pbi->common;
  MACROBLOCKD *const xd = &pbi->mb;
  MODE_INFO *mi = xd->mode_info_context;
  const BLOCK_SIZE_TYPE bsize = mi->mbmi.sb_type;
  const int bw = 1 << mi_width_log2(bsize);
  const int bh = 1 << mi_height_log2(bsize);
  const int y_mis = MIN(bh, cm->mi_rows - mi_row);
  const int x_mis = MIN(bw, cm->mi_cols - mi_col);
  int x, y;

  if (cm->frame_type == KEY_FRAME || cm->intra_only)
    read_intra_frame_mode_info(pbi, mi, mi_row, mi_col, r);
  else
    read_inter_frame_mode_info(pbi, mi, mi_row, mi_col, r);

  for (y = 0; y < y_mis; y++)
    for (x = !y; x < x_mis; x++)
      mi[y * cm->mode_info_stride + x] = *mi;
}
