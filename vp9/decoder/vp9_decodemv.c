/*
  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */


#include "vp9/decoder/vp9_treereader.h"
#include "vp9/common/vp9_entropymv.h"
#include "vp9/common/vp9_entropymode.h"
#include "vp9/common/vp9_reconinter.h"
#include "vp9/decoder/vp9_onyxd_int.h"
#include "vp9/common/vp9_findnearmv.h"
#include "vp9/common/vp9_common.h"
#include "vp9/common/vp9_seg_common.h"
#include "vp9/common/vp9_pred_common.h"
#include "vp9/common/vp9_entropy.h"
#include "vp9/decoder/vp9_decodemv.h"
#include "vp9/common/vp9_mvref_common.h"
#if CONFIG_DEBUG
#include <assert.h>
#endif

// #define DEBUG_DEC_MV
#ifdef DEBUG_DEC_MV
int dec_mvcount = 0;
#endif

// #define DEC_DEBUG
#ifdef DEC_DEBUG
extern int dec_debug;
#endif

static B_PREDICTION_MODE read_bmode(vp9_reader *r, const vp9_prob *p) {
  B_PREDICTION_MODE m = treed_read(r, vp9_bmode_tree, p);
  return m;
}

static B_PREDICTION_MODE read_kf_bmode(vp9_reader *r, const vp9_prob *p) {
  return (B_PREDICTION_MODE)treed_read(r, vp9_kf_bmode_tree, p);
}

static MB_PREDICTION_MODE read_ymode(vp9_reader *r, const vp9_prob *p) {
  return (MB_PREDICTION_MODE)treed_read(r, vp9_ymode_tree, p);
}

static MB_PREDICTION_MODE read_sb_ymode(vp9_reader *r, const vp9_prob *p) {
  return (MB_PREDICTION_MODE)treed_read(r, vp9_sb_ymode_tree, p);
}

static MB_PREDICTION_MODE read_kf_sb_ymode(vp9_reader *r, const vp9_prob *p) {
  return (MB_PREDICTION_MODE)treed_read(r, vp9_uv_mode_tree, p);
}

static MB_PREDICTION_MODE read_kf_mb_ymode(vp9_reader *r, const vp9_prob *p) {
  return (MB_PREDICTION_MODE)treed_read(r, vp9_kf_ymode_tree, p);
}

static MB_PREDICTION_MODE read_uv_mode(vp9_reader *r, const vp9_prob *p) {
  return (MB_PREDICTION_MODE)treed_read(r, vp9_uv_mode_tree, p);
}

static int read_mb_segid(vp9_reader *r, MACROBLOCKD *xd) {
  return treed_read(r, vp9_segment_tree, xd->mb_segment_tree_probs);
}

static void set_segment_id(VP9_COMMON *cm, MB_MODE_INFO *mbmi,
                           int mi_row, int mi_col, int segment_id) {
  const int mi_index = mi_row * cm->mi_cols + mi_col;
  const BLOCK_SIZE_TYPE sb_type = mbmi->sb_type;
  const int bw = 1 << mi_width_log2(sb_type);
  const int bh = 1 << mi_height_log2(sb_type);
  const int ymis = MIN(cm->mi_rows - mi_row, bh);
  const int xmis = MIN(cm->mi_cols - mi_col, bw);
  int x, y;

  for (y = 0; y < ymis; y++) {
    for (x = 0; x < xmis; x++) {
      const int index = mi_index + (y * cm->mi_cols + x);
      cm->last_frame_seg_map[index] = segment_id;
    }
  }
}

static TX_SIZE select_txfm_size(VP9_COMMON *cm, vp9_reader *r,
                                int allow_16x16, int allow_32x32) {
  TX_SIZE txfm_size = vp9_read(r, cm->prob_tx[0]);  // TX_4X4 or >TX_4X4
  if (txfm_size != TX_4X4 && allow_16x16) {
    txfm_size += vp9_read(r, cm->prob_tx[1]);       // TX_8X8 or >TX_8X8
    if (txfm_size != TX_8X8 && allow_32x32)
      txfm_size += vp9_read(r, cm->prob_tx[2]);     // TX_16X16 or >TX_16X16
  }
  return txfm_size;
}


static void kfread_modes(VP9D_COMP *pbi, MODE_INFO *m,
                         int mi_row, int mi_col,
                         vp9_reader *r) {
  VP9_COMMON *const cm = &pbi->common;
  MACROBLOCKD *const xd = &pbi->mb;
  m->mbmi.ref_frame = INTRA_FRAME;

  // Read segmentation map if it is being updated explicitly this frame
  m->mbmi.segment_id = 0;
  if (xd->segmentation_enabled && xd->update_mb_segmentation_map) {
    m->mbmi.segment_id = read_mb_segid(r, xd);
    set_segment_id(cm, &m->mbmi, mi_row, mi_col, m->mbmi.segment_id);
  }

  m->mbmi.mb_skip_coeff = vp9_segfeature_active(xd, m->mbmi.segment_id,
                                                SEG_LVL_SKIP);
  if (!m->mbmi.mb_skip_coeff)
    m->mbmi.mb_skip_coeff = vp9_read(r, vp9_get_pred_prob(cm, xd, PRED_MBSKIP));

  // luma mode
#if CONFIG_AB4X4
  if (m->mbmi.sb_type >= BLOCK_SIZE_SB8X8)
    m->mbmi.mode = read_kf_sb_ymode(r,
                     cm->sb_kf_ymode_prob[cm->kf_ymode_probs_index]);
  else
     m->mbmi.mode = I4X4_PRED;
#else
  m->mbmi.mode = m->mbmi.sb_type > BLOCK_SIZE_SB8X8 ?
      read_kf_sb_ymode(r, cm->sb_kf_ymode_prob[cm->kf_ymode_probs_index]):
      read_kf_mb_ymode(r, cm->kf_ymode_prob[cm->kf_ymode_probs_index]);
#endif

  m->mbmi.ref_frame = INTRA_FRAME;

#if CONFIG_AB4X4
  if (m->mbmi.sb_type < BLOCK_SIZE_SB8X8) {
#else
  if (m->mbmi.mode == I4X4_PRED) {
#endif
    int idx, idy;
    int bw = 1 << b_width_log2(m->mbmi.sb_type);
    int bh = 1 << b_height_log2(m->mbmi.sb_type);
    // FIXME(jingning): fix intra4x4 rate-distortion optimization, then
    // use bw and bh as the increment values.
#if !CONFIG_AB4X4 || CONFIG_AB4X4
    bw = 1, bh = 1;
#endif
    for (idy = 0; idy < 2; idy += bh)
      for (idx = 0; idx < 2; idx += bw)
        m->bmi[idy * 2 + idx].as_mode.first =
            read_kf_sb_ymode(r, cm->sb_kf_ymode_prob[cm->kf_ymode_probs_index]);
  }

  m->mbmi.uv_mode = read_uv_mode(r, cm->kf_uv_mode_prob[m->mbmi.mode]);

  if (cm->txfm_mode == TX_MODE_SELECT &&
      !m->mbmi.mb_skip_coeff &&
#if CONFIG_AB4X4
      m->mbmi.sb_type >= BLOCK_SIZE_SB8X8
#else
      m->mbmi.mode != I4X4_PRED
#endif
      ) {
    const int allow_16x16 = m->mbmi.sb_type >= BLOCK_SIZE_MB16X16;
    const int allow_32x32 = m->mbmi.sb_type >= BLOCK_SIZE_SB32X32;
    m->mbmi.txfm_size = select_txfm_size(cm, r, allow_16x16, allow_32x32);
  } else if (cm->txfm_mode >= ALLOW_32X32 &&
             m->mbmi.sb_type >= BLOCK_SIZE_SB32X32) {
    m->mbmi.txfm_size = TX_32X32;
  } else if (cm->txfm_mode >= ALLOW_16X16 &&
             m->mbmi.sb_type >= BLOCK_SIZE_MB16X16 &&
             m->mbmi.mode <= TM_PRED) {
    m->mbmi.txfm_size = TX_16X16;
  } else if (cm->txfm_mode >= ALLOW_8X8 &&
#if CONFIG_AB4X4
             m->mbmi.sb_type >= BLOCK_SIZE_SB8X8
#else
             m->mbmi.mode != I4X4_PRED
#endif
             ) {
    m->mbmi.txfm_size = TX_8X8;
  } else {
    m->mbmi.txfm_size = TX_4X4;
  }
}

static int read_mv_component(vp9_reader *r,
                             const nmv_component *mvcomp, int usehp) {

  int mag, d, fr, hp;
  const int sign = vp9_read(r, mvcomp->sign);
  const int mv_class = treed_read(r, vp9_mv_class_tree, mvcomp->classes);

  // Integer part
  if (mv_class == MV_CLASS_0) {
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
                  mv_class == MV_CLASS_0 ? mvcomp->class0_fp[d] : mvcomp->fp);


  // High precision part (if hp is not used, the default value of the hp is 1)
  hp = usehp ? vp9_read(r,
                        mv_class == MV_CLASS_0 ? mvcomp->class0_hp : mvcomp->hp)
             : 1;

  // result
  mag = vp9_get_mv_mag(mv_class, (d << 3) | (fr << 1) | hp) + 1;
  return sign ? -mag : mag;
}

static void update_nmv(vp9_reader *r, vp9_prob *const p,
                       const vp9_prob upd_p) {
  if (vp9_read(r, upd_p)) {
#ifdef LOW_PRECISION_MV_UPDATE
    *p = (vp9_read_literal(r, 7) << 1) | 1;
#else
    *p = (vp9_read_literal(r, 8));
#endif
  }
}

static void read_nmvprobs(vp9_reader *r, nmv_context *mvctx,
                          int usehp) {
  int i, j, k;

#ifdef MV_GROUP_UPDATE
  if (!vp9_read_bit(r))
    return;
#endif
  for (j = 0; j < MV_JOINTS - 1; ++j)
    update_nmv(r, &mvctx->joints[j], VP9_NMV_UPDATE_PROB);

  for (i = 0; i < 2; ++i) {
    update_nmv(r, &mvctx->comps[i].sign, VP9_NMV_UPDATE_PROB);
    for (j = 0; j < MV_CLASSES - 1; ++j)
      update_nmv(r, &mvctx->comps[i].classes[j], VP9_NMV_UPDATE_PROB);

    for (j = 0; j < CLASS0_SIZE - 1; ++j)
      update_nmv(r, &mvctx->comps[i].class0[j], VP9_NMV_UPDATE_PROB);

    for (j = 0; j < MV_OFFSET_BITS; ++j)
      update_nmv(r, &mvctx->comps[i].bits[j], VP9_NMV_UPDATE_PROB);
  }

  for (i = 0; i < 2; ++i) {
    for (j = 0; j < CLASS0_SIZE; ++j)
      for (k = 0; k < 3; ++k)
        update_nmv(r, &mvctx->comps[i].class0_fp[j][k], VP9_NMV_UPDATE_PROB);

    for (j = 0; j < 3; ++j)
      update_nmv(r, &mvctx->comps[i].fp[j], VP9_NMV_UPDATE_PROB);
  }

  if (usehp) {
    for (i = 0; i < 2; ++i) {
      update_nmv(r, &mvctx->comps[i].class0_hp, VP9_NMV_UPDATE_PROB);
      update_nmv(r, &mvctx->comps[i].hp, VP9_NMV_UPDATE_PROB);
    }
  }
}

// Read the referncence frame
static MV_REFERENCE_FRAME read_ref_frame(VP9D_COMP *pbi,
                                         vp9_reader *r,
                                         int segment_id) {
  MV_REFERENCE_FRAME ref_frame;
  VP9_COMMON *const cm = &pbi->common;
  MACROBLOCKD *const xd = &pbi->mb;

  int seg_ref_count = 0;
  const int seg_ref_active = vp9_segfeature_active(xd, segment_id,
                                                   SEG_LVL_REF_FRAME);

  const int intra = vp9_check_segref(xd, segment_id, INTRA_FRAME);
  const int last = vp9_check_segref(xd, segment_id, LAST_FRAME);
  const int golden = vp9_check_segref(xd, segment_id, GOLDEN_FRAME);
  const int altref = vp9_check_segref(xd, segment_id, ALTREF_FRAME);

  // If segment coding enabled does the segment allow for more than one
  // possible reference frame
  if (seg_ref_active)
    seg_ref_count = intra + last + golden + altref;

  // Segment reference frame features not available or allows for
  // multiple reference frame options
  if (!seg_ref_active || seg_ref_count > 1) {
    // Values used in prediction model coding
    MV_REFERENCE_FRAME pred_ref;

    // Get the context probability the prediction flag
    vp9_prob pred_prob = vp9_get_pred_prob(cm, xd, PRED_REF);

    // Read the prediction status flag
    unsigned char prediction_flag = vp9_read(r, pred_prob);

    // Store the prediction flag.
    vp9_set_pred_flag(xd, PRED_REF, prediction_flag);

    // Get the predicted reference frame.
    pred_ref = vp9_get_pred_ref(cm, xd);

    // If correctly predicted then use the predicted value
    if (prediction_flag) {
      ref_frame = pred_ref;
    } else {
      // decode the explicitly coded value
      vp9_prob mod_refprobs[PREDICTION_PROBS];
      vpx_memcpy(mod_refprobs, cm->mod_refprobs[pred_ref],
                 sizeof(mod_refprobs));

      // If segment coding enabled blank out options that cant occur by
      // setting the branch probability to 0.
      if (seg_ref_active) {
        mod_refprobs[INTRA_FRAME] *= intra;
        mod_refprobs[LAST_FRAME] *= last;
        mod_refprobs[GOLDEN_FRAME] *= golden * altref;
      }

      // Default to INTRA_FRAME (value 0)
      ref_frame = INTRA_FRAME;

      // Do we need to decode the Intra/Inter branch
      if (mod_refprobs[0])
        ref_frame = vp9_read(r, mod_refprobs[0]);
      else
        ref_frame++;

      if (ref_frame) {
        // Do we need to decode the Last/Gf_Arf branch
        if (mod_refprobs[1])
          ref_frame += vp9_read(r, mod_refprobs[1]);
        else
          ref_frame++;

        if (ref_frame > 1) {
          // Do we need to decode the GF/Arf branch
          if (mod_refprobs[2]) {
            ref_frame += vp9_read(r, mod_refprobs[2]);
          } else {
            if (seg_ref_active)
              ref_frame = pred_ref == GOLDEN_FRAME || !golden ? ALTREF_FRAME
                                                              : GOLDEN_FRAME;
            else
              ref_frame = pred_ref == GOLDEN_FRAME ? ALTREF_FRAME
                                                   : GOLDEN_FRAME;
          }
        }
      }
    }
  } else {
    // Segment reference frame features are enabled
    // The reference frame for the mb is considered as correclty predicted
    // if it is signaled at the segment level for the purposes of the
    // common prediction model
    vp9_set_pred_flag(xd, PRED_REF, 1);
    ref_frame = vp9_get_pred_ref(cm, xd);
  }

  return ref_frame;
}

static MB_PREDICTION_MODE read_sb_mv_ref(vp9_reader *r, const vp9_prob *p) {
  return (MB_PREDICTION_MODE) treed_read(r, vp9_sb_mv_ref_tree, p);
}

static MB_PREDICTION_MODE read_mv_ref(vp9_reader *r, const vp9_prob *p) {
  return (MB_PREDICTION_MODE) treed_read(r, vp9_mv_ref_tree, p);
}

static B_PREDICTION_MODE read_sub_mv_ref(vp9_reader *r, const vp9_prob *p) {
  return (B_PREDICTION_MODE) treed_read(r, vp9_sub_mv_ref_tree, p);
}

#ifdef VPX_MODE_COUNT
unsigned int vp9_mv_cont_count[5][4] = {
  { 0, 0, 0, 0 },
  { 0, 0, 0, 0 },
  { 0, 0, 0, 0 },
  { 0, 0, 0, 0 },
  { 0, 0, 0, 0 }
};
#endif

static void read_switchable_interp_probs(VP9D_COMP* const pbi, vp9_reader *r) {
  VP9_COMMON *const cm = &pbi->common;
  int i, j;
  for (j = 0; j < VP9_SWITCHABLE_FILTERS + 1; ++j)
    for (i = 0; i < VP9_SWITCHABLE_FILTERS - 1; ++i)
      cm->fc.switchable_interp_prob[j][i] = vp9_read_prob(r);
}

static INLINE COMPPREDMODE_TYPE read_comp_pred_mode(vp9_reader *r) {
  COMPPREDMODE_TYPE mode = vp9_read_bit(r);
  if (mode)
     mode += vp9_read_bit(r);
  return mode;
}

static void mb_mode_mv_init(VP9D_COMP *pbi, vp9_reader *r) {
  VP9_COMMON *const cm = &pbi->common;

  if (cm->frame_type == KEY_FRAME) {
    if (!cm->kf_ymode_probs_update)
      cm->kf_ymode_probs_index = vp9_read_literal(r, 3);
  } else {
    nmv_context *const nmvc = &pbi->common.fc.nmvc;
    MACROBLOCKD *const xd = &pbi->mb;
    int i, j;

    if (cm->mcomp_filter_type == SWITCHABLE)
      read_switchable_interp_probs(pbi, r);

    // Baseline probabilities for decoding reference frame
    cm->prob_intra_coded = vp9_read_prob(r);
    cm->prob_last_coded  = vp9_read_prob(r);
    cm->prob_gf_coded    = vp9_read_prob(r);

    // Computes a modified set of probabilities for use when reference
    // frame prediction fails.
    vp9_compute_mod_refprobs(cm);

    cm->comp_pred_mode = read_comp_pred_mode(r);
    if (cm->comp_pred_mode == HYBRID_PREDICTION)
      for (i = 0; i < COMP_PRED_CONTEXTS; i++)
        cm->prob_comppred[i] = vp9_read_prob(r);

    // VP9_YMODES
    if (vp9_read_bit(r))
      for (i = 0; i < VP9_YMODES - 1; ++i)
        cm->fc.ymode_prob[i] = vp9_read_prob(r);

    // VP9_I32X32_MODES
    if (vp9_read_bit(r))
      for (i = 0; i < VP9_I32X32_MODES - 1; ++i)
        cm->fc.sb_ymode_prob[i] = vp9_read_prob(r);

    for (j = 0; j < NUM_PARTITION_CONTEXTS; ++j)
      if (vp9_read_bit(r))
        for (i = 0; i < PARTITION_TYPES - 1; ++i)
          cm->fc.partition_prob[j][i] = vp9_read_prob(r);

    read_nmvprobs(r, nmvc, xd->allow_high_precision_mv);
  }
}

// This function either reads the segment id for the current macroblock from
// the bitstream or if the value is temporally predicted asserts the predicted
// value
static int read_mb_segment_id(VP9D_COMP *pbi, int mi_row, int mi_col,
                              vp9_reader *r) {
  VP9_COMMON *const cm = &pbi->common;
  MACROBLOCKD *const xd = &pbi->mb;
  MODE_INFO *const mi = xd->mode_info_context;
  MB_MODE_INFO *const mbmi = &mi->mbmi;

  if (!xd->segmentation_enabled)
    return 0;  // Default for disabled segmentation

  if (xd->update_mb_segmentation_map) {
    int segment_id;

    if (cm->temporal_update) {
      // Temporal coding of the segment id for this mb is enabled.
      // Get the context based probability for reading the
      // prediction status flag
      const vp9_prob pred_prob = vp9_get_pred_prob(cm, xd, PRED_SEG_ID);
      const int pred_flag = vp9_read(r, pred_prob);
      vp9_set_pred_flag(xd, PRED_SEG_ID, pred_flag);

      // If the value is flagged as correctly predicted
      // then use the predicted value, otherwise decode it explicitly
      segment_id = pred_flag ? vp9_get_pred_mi_segid(cm, mbmi->sb_type,
                                                     mi_row, mi_col)
                             : read_mb_segid(r, xd);
    } else {
      segment_id = read_mb_segid(r, xd);  // Normal unpredicted coding mode
    }

    set_segment_id(cm, mbmi, mi_row, mi_col, segment_id);  // Side effect
    return segment_id;
  } else {
    return vp9_get_pred_mi_segid(cm, mbmi->sb_type, mi_row, mi_col);
  }
}


static INLINE void assign_and_clamp_mv(int_mv *dst, const int_mv *src,
                                       int mb_to_left_edge,
                                       int mb_to_right_edge,
                                       int mb_to_top_edge,
                                       int mb_to_bottom_edge) {
  dst->as_int = src->as_int;
  clamp_mv(dst, mb_to_left_edge, mb_to_right_edge, mb_to_top_edge,
           mb_to_bottom_edge);
}

static INLINE void decode_mv(vp9_reader *r, MV *mv, const MV *ref,
                             const nmv_context *ctx,
                             nmv_context_counts *counts,
                             int usehp) {
  const MV_JOINT_TYPE j = treed_read(r, vp9_mv_joint_tree, ctx->joints);
  MV diff = {0, 0};

  usehp = usehp && vp9_use_nmv_hp(ref);
  if (mv_joint_vertical(j))
    diff.row = read_mv_component(r, &ctx->comps[0], usehp);

  if (mv_joint_horizontal(j))
    diff.col = read_mv_component(r, &ctx->comps[1], usehp);

  vp9_increment_nmv(&diff, ref, counts, usehp);

  mv->row = diff.row + ref->row;
  mv->col = diff.col + ref->col;
}

static INLINE INTERPOLATIONFILTERTYPE read_switchable_filter_type(
    VP9D_COMP *pbi, vp9_reader *r) {
  const int index = treed_read(r, vp9_switchable_interp_tree,
                               vp9_get_pred_probs(&pbi->common, &pbi->mb,
                                                  PRED_SWITCHABLE_INTERP));
  return vp9_switchable_interp[index];
}

static void read_mb_modes_mv(VP9D_COMP *pbi, MODE_INFO *mi, MB_MODE_INFO *mbmi,
                             int mi_row, int mi_col,
                             vp9_reader *r) {
  VP9_COMMON *const cm = &pbi->common;
  nmv_context *const nmvc = &cm->fc.nmvc;
  const int mis = cm->mode_info_stride;
  MACROBLOCKD *const xd = &pbi->mb;

  int_mv *const mv0 = &mbmi->mv[0];
  int_mv *const mv1 = &mbmi->mv[1];
  BLOCK_SIZE_TYPE bsize = mi->mbmi.sb_type;
  int bw = 1 << b_width_log2(bsize);
  int bh = 1 << b_height_log2(bsize);

  const int use_prev_in_find_mv_refs = cm->width == cm->last_width &&
                                       cm->height == cm->last_height &&
                                       !cm->error_resilient_mode &&
                                       cm->last_show_frame;

  int mb_to_left_edge, mb_to_right_edge, mb_to_top_edge, mb_to_bottom_edge;
  int j, idx, idy;

  mbmi->need_to_clamp_mvs = 0;
  mbmi->need_to_clamp_secondmv = 0;
  mbmi->second_ref_frame = NONE;

  // Make sure the MACROBLOCKD mode info pointer is pointed at the
  // correct entry for the current macroblock.
  xd->mode_info_context = mi;

  // Distance of Mb to the various image edges.
  // These specified to 8th pel as they are always compared to MV values
  // that are in 1/8th pel units
  set_mi_row_col(cm, xd, mi_row, 1 << mi_height_log2(bsize),
                         mi_col, 1 << mi_width_log2(bsize));

  mb_to_top_edge = xd->mb_to_top_edge - LEFT_TOP_MARGIN;
  mb_to_bottom_edge = xd->mb_to_bottom_edge + RIGHT_BOTTOM_MARGIN;
  mb_to_left_edge = xd->mb_to_left_edge - LEFT_TOP_MARGIN;
  mb_to_right_edge = xd->mb_to_right_edge + RIGHT_BOTTOM_MARGIN;

  // Read the macroblock segment id.
  mbmi->segment_id = read_mb_segment_id(pbi, mi_row, mi_col, r);

  mbmi->mb_skip_coeff = vp9_segfeature_active(xd, mbmi->segment_id,
                                              SEG_LVL_SKIP);
  if (!mbmi->mb_skip_coeff)
    mbmi->mb_skip_coeff = vp9_read(r, vp9_get_pred_prob(cm, xd, PRED_MBSKIP));

  // Read the reference frame
  mbmi->ref_frame = read_ref_frame(pbi, r, mbmi->segment_id);

  // If reference frame is an Inter frame
  if (mbmi->ref_frame) {
    int_mv nearest, nearby, best_mv;
    int_mv nearest_second, nearby_second, best_mv_second;
    vp9_prob mv_ref_p[VP9_MVREFS - 1];

    const MV_REFERENCE_FRAME ref_frame = mbmi->ref_frame;
    struct scale_factors *sf0 = &xd->scale_factor[0];
    *sf0 = cm->active_ref_scale[mbmi->ref_frame - 1];

    {
      // Select the appropriate reference frame for this MB
      const int ref_fb_idx = cm->active_ref_idx[ref_frame - 1];

      setup_pre_planes(xd, &cm->yv12_fb[ref_fb_idx], NULL,
                       mi_row, mi_col, xd->scale_factor, xd->scale_factor_uv);

#ifdef DEC_DEBUG
      if (dec_debug)
        printf("%d %d\n", xd->mode_info_context->mbmi.mv[0].as_mv.row,
               xd->mode_info_context->mbmi.mv[0].as_mv.col);
#endif
      vp9_find_mv_refs(cm, xd, mi, use_prev_in_find_mv_refs ?
                       xd->prev_mode_info_context : NULL,
                       ref_frame, mbmi->ref_mvs[ref_frame],
                       cm->ref_frame_sign_bias);

      vp9_mv_ref_probs(cm, mv_ref_p, mbmi->mb_mode_context[ref_frame]);

      // If the segment level skip mode enabled
      if (vp9_segfeature_active(xd, mbmi->segment_id, SEG_LVL_SKIP)) {
        mbmi->mode = ZEROMV;
      } else {
#if CONFIG_AB4X4
        if (bsize >= BLOCK_SIZE_SB8X8)
          mbmi->mode = read_sb_mv_ref(r, mv_ref_p);
        else
          mbmi->mode = SPLITMV;
#else
        mbmi->mode = bsize > BLOCK_SIZE_SB8X8 ?
                                   read_sb_mv_ref(r, mv_ref_p)
                                 : read_mv_ref(r, mv_ref_p);
#endif
        vp9_accum_mv_refs(cm, mbmi->mode, mbmi->mb_mode_context[ref_frame]);
      }

      if (mbmi->mode != ZEROMV) {
        vp9_find_best_ref_mvs(xd,
                              mbmi->ref_mvs[ref_frame],
                              &nearest, &nearby);

        best_mv.as_int = mbmi->ref_mvs[ref_frame][0].as_int;
      }

#ifdef DEC_DEBUG
      if (dec_debug)
        printf("[D %d %d] %d %d %d %d\n", ref_frame,
               mbmi->mb_mode_context[ref_frame],
               mv_ref_p[0], mv_ref_p[1], mv_ref_p[2], mv_ref_p[3]);
#endif
    }

    mbmi->interp_filter = cm->mcomp_filter_type == SWITCHABLE
                              ? read_switchable_filter_type(pbi, r)
                              : cm->mcomp_filter_type;

    if (cm->comp_pred_mode == COMP_PREDICTION_ONLY ||
        (cm->comp_pred_mode == HYBRID_PREDICTION &&
         vp9_read(r, vp9_get_pred_prob(cm, xd, PRED_COMP)))) {
      /* Since we have 3 reference frames, we can only have 3 unique
       * combinations of combinations of 2 different reference frames
       * (A-G, G-L or A-L). In the bitstream, we use this to simply
       * derive the second reference frame from the first reference
       * frame, by saying it's the next one in the enumerator, and
       * if that's > n_refs, then the second reference frame is the
       * first one in the enumerator. */
      mbmi->second_ref_frame = mbmi->ref_frame + 1;
      if (mbmi->second_ref_frame == 4)
        mbmi->second_ref_frame = 1;
      if (mbmi->second_ref_frame > 0) {
        const MV_REFERENCE_FRAME second_ref_frame = mbmi->second_ref_frame;
        struct scale_factors *sf1 = &xd->scale_factor[1];
        const int second_ref_fb_idx = cm->active_ref_idx[second_ref_frame - 1];
        *sf1 = cm->active_ref_scale[second_ref_frame - 1];

        setup_pre_planes(xd, NULL, &cm->yv12_fb[second_ref_fb_idx],
                         mi_row, mi_col, xd->scale_factor, xd->scale_factor_uv);

        vp9_find_mv_refs(cm, xd, mi,
                         use_prev_in_find_mv_refs ?
                         xd->prev_mode_info_context : NULL,
                         second_ref_frame, mbmi->ref_mvs[second_ref_frame],
                         cm->ref_frame_sign_bias);

        if (mbmi->mode != ZEROMV) {
          vp9_find_best_ref_mvs(xd,
                                mbmi->ref_mvs[second_ref_frame],
                                &nearest_second,
                                &nearby_second);
          best_mv_second.as_int = mbmi->ref_mvs[second_ref_frame][0].as_int;
        }
      }

    }

    mbmi->uv_mode = DC_PRED;
    switch (mbmi->mode) {
      case SPLITMV:
#if !CONFIG_AB4X4
        bw = 1, bh = 1;
#endif
        mbmi->need_to_clamp_mvs = 0;
        for (idy = 0; idy < 2; idy += bh) {
          for (idx = 0; idx < 2; idx += bw) {
            int_mv leftmv, abovemv, second_leftmv, second_abovemv;
            int_mv blockmv, secondmv;
            int mv_contz;
            int blockmode;
            int i, k;
            j = idy * 2 + idx;
            k = j;

            leftmv.as_int = left_block_mv(xd, mi, k);
            abovemv.as_int = above_block_mv(mi, k, mis);
            second_leftmv.as_int = 0;
            second_abovemv.as_int = 0;
            if (mbmi->second_ref_frame > 0) {
              second_leftmv.as_int = left_block_second_mv(xd, mi, k);
              second_abovemv.as_int = above_block_second_mv(mi, k, mis);
            }
            mv_contz = vp9_mv_cont(&leftmv, &abovemv);
            blockmode = read_sub_mv_ref(r, cm->fc.sub_mv_ref_prob[mv_contz]);
            cm->fc.sub_mv_ref_counts[mv_contz][blockmode - LEFT4X4]++;

            switch (blockmode) {
              case NEW4X4:
                decode_mv(r, &blockmv.as_mv, &best_mv.as_mv, nmvc,
                           &cm->fc.NMVcount, xd->allow_high_precision_mv);

                if (mbmi->second_ref_frame > 0)
                  decode_mv(r, &secondmv.as_mv, &best_mv_second.as_mv, nmvc,
                            &cm->fc.NMVcount, xd->allow_high_precision_mv);

  #ifdef VPX_MODE_COUNT
                vp9_mv_cont_count[mv_contz][3]++;
  #endif
                break;
              case LEFT4X4:
                blockmv.as_int = leftmv.as_int;
                if (mbmi->second_ref_frame > 0)
                  secondmv.as_int = second_leftmv.as_int;
  #ifdef VPX_MODE_COUNT
                vp9_mv_cont_count[mv_contz][0]++;
  #endif
                break;
              case ABOVE4X4:
                blockmv.as_int = abovemv.as_int;
                if (mbmi->second_ref_frame > 0)
                  secondmv.as_int = second_abovemv.as_int;
  #ifdef VPX_MODE_COUNT
                vp9_mv_cont_count[mv_contz][1]++;
  #endif
                break;
              case ZERO4X4:
                blockmv.as_int = 0;
                if (mbmi->second_ref_frame > 0)
                  secondmv.as_int = 0;
  #ifdef VPX_MODE_COUNT
                vp9_mv_cont_count[mv_contz][2]++;
  #endif
                break;
              default:
                break;
            }
            mi->bmi[j].as_mv[0].as_int = blockmv.as_int;
            if (mbmi->second_ref_frame > 0)
              mi->bmi[j].as_mv[1].as_int = secondmv.as_int;

            for (i = 1; i < bh; ++i)
              vpx_memcpy(&mi->bmi[j + i * 2], &mi->bmi[j], sizeof(mi->bmi[j]));
            for (i = 1; i < bw; ++i)
              vpx_memcpy(&mi->bmi[j + i], &mi->bmi[j], sizeof(mi->bmi[j]));
          }
        }

        mv0->as_int = mi->bmi[3].as_mv[0].as_int;
        mv1->as_int = mi->bmi[3].as_mv[1].as_int;
        break;  /* done with SPLITMV */

      case NEARMV:
        // Clip "next_nearest" so that it does not extend to far out of image
        assign_and_clamp_mv(mv0, &nearby, mb_to_left_edge,
                                          mb_to_right_edge,
                                          mb_to_top_edge,
                                          mb_to_bottom_edge);
        if (mbmi->second_ref_frame > 0)
          assign_and_clamp_mv(mv1, &nearby_second, mb_to_left_edge,
                                                   mb_to_right_edge,
                                                   mb_to_top_edge,
                                                   mb_to_bottom_edge);
        break;

      case NEARESTMV:
        // Clip "next_nearest" so that it does not extend to far out of image
        assign_and_clamp_mv(mv0, &nearest, mb_to_left_edge,
                                           mb_to_right_edge,
                                           mb_to_top_edge,
                                           mb_to_bottom_edge);
        if (mbmi->second_ref_frame > 0)
          assign_and_clamp_mv(mv1, &nearest_second, mb_to_left_edge,
                                                    mb_to_right_edge,
                                                    mb_to_top_edge,
                                                    mb_to_bottom_edge);
        break;

      case ZEROMV:
        mv0->as_int = 0;
        if (mbmi->second_ref_frame > 0)
          mv1->as_int = 0;
        break;

      case NEWMV:
        decode_mv(r, &mv0->as_mv, &best_mv.as_mv, nmvc, &cm->fc.NMVcount,
                  xd->allow_high_precision_mv);
        mbmi->need_to_clamp_mvs = check_mv_bounds(mv0,
                                                  mb_to_left_edge,
                                                  mb_to_right_edge,
                                                  mb_to_top_edge,
                                                  mb_to_bottom_edge);

        if (mbmi->second_ref_frame > 0) {
          decode_mv(r, &mv1->as_mv, &best_mv_second.as_mv, nmvc,
                    &cm->fc.NMVcount, xd->allow_high_precision_mv);
          mbmi->need_to_clamp_secondmv = check_mv_bounds(mv1,
                                                         mb_to_left_edge,
                                                         mb_to_right_edge,
                                                         mb_to_top_edge,
                                                         mb_to_bottom_edge);
        }
        break;
      default:
;
#if CONFIG_DEBUG
        assert(0);
#endif
    }
  } else {
    // required for left and above block mv
    mv0->as_int = 0;

#if CONFIG_AB4X4
    if (bsize >= BLOCK_SIZE_SB8X8) {
      mbmi->mode = read_sb_ymode(r, cm->fc.sb_ymode_prob);
      cm->fc.sb_ymode_counts[mbmi->mode]++;
    } else {
      mbmi->mode = I4X4_PRED;
    }
#else
    if (bsize > BLOCK_SIZE_SB8X8) {
      mbmi->mode = read_sb_ymode(r, cm->fc.sb_ymode_prob);
      cm->fc.sb_ymode_counts[mbmi->mode]++;
    } else {
      mbmi->mode = read_ymode(r, cm->fc.ymode_prob);
      cm->fc.ymode_counts[mbmi->mode]++;
    }
#endif

    // If MB mode is I4X4_PRED read the block modes
#if CONFIG_AB4X4
    if (bsize < BLOCK_SIZE_SB8X8) {
#else
    if (mbmi->mode == I4X4_PRED) {
#endif
      int idx, idy;
      // FIXME(jingning): fix intra4x4 rate-distortion optimization, then
      // use bw and bh as the increment values.
#if !CONFIG_AB4X4 || CONFIG_AB4X4
      bw = 1, bh = 1;
#endif
      for (idy = 0; idy < 2; idy += bh) {
        for (idx = 0; idx < 2; idx += bw) {
          int m = read_sb_ymode(r, cm->fc.sb_ymode_prob);
          mi->bmi[idy * 2 + idx].as_mode.first = m;
          cm->fc.sb_ymode_counts[m]++;
        }
      }
    }

    mbmi->uv_mode = read_uv_mode(r, cm->fc.uv_mode_prob[mbmi->mode]);
    cm->fc.uv_mode_counts[mbmi->mode][mbmi->uv_mode]++;
  }

#if CONFIG_AB4X4
  if (cm->txfm_mode == TX_MODE_SELECT && mbmi->mb_skip_coeff == 0 &&
      bsize >= BLOCK_SIZE_SB8X8) {
#else
  if (cm->txfm_mode == TX_MODE_SELECT && mbmi->mb_skip_coeff == 0 &&
      ((mbmi->ref_frame == INTRA_FRAME && mbmi->mode != I4X4_PRED) ||
       (mbmi->ref_frame != INTRA_FRAME && mbmi->mode != SPLITMV))) {
#endif
    const int allow_16x16 = bsize >= BLOCK_SIZE_MB16X16;
    const int allow_32x32 = bsize >= BLOCK_SIZE_SB32X32;
    mbmi->txfm_size = select_txfm_size(cm, r, allow_16x16, allow_32x32);
  } else if (bsize >= BLOCK_SIZE_SB32X32 &&
             cm->txfm_mode >= ALLOW_32X32) {
    mbmi->txfm_size = TX_32X32;
  } else if (cm->txfm_mode >= ALLOW_16X16 &&
             bsize >= BLOCK_SIZE_MB16X16
#if !CONFIG_AB4X4
      && ((mbmi->ref_frame == INTRA_FRAME && mbmi->mode <= TM_PRED) ||
       (mbmi->ref_frame != INTRA_FRAME && mbmi->mode != SPLITMV))
#endif
       ) {
    mbmi->txfm_size = TX_16X16;
  } else if (cm->txfm_mode >= ALLOW_8X8 &&
#if CONFIG_AB4X4
      (bsize >= BLOCK_SIZE_SB8X8))
#else
      (!(mbmi->ref_frame == INTRA_FRAME && mbmi->mode == I4X4_PRED) &&
       !(mbmi->ref_frame != INTRA_FRAME && mbmi->mode == SPLITMV)))
#endif
  {
    mbmi->txfm_size = TX_8X8;
  } else {
    mbmi->txfm_size = TX_4X4;
  }
}

void vp9_decode_mode_mvs_init(VP9D_COMP* const pbi, vp9_reader *r) {
  VP9_COMMON *cm = &pbi->common;
  int k;

  // TODO(jkoleszar): does this clear more than MBSKIP_CONTEXTS? Maybe remove.
  vpx_memset(cm->mbskip_pred_probs, 0, sizeof(cm->mbskip_pred_probs));
  for (k = 0; k < MBSKIP_CONTEXTS; ++k)
    cm->mbskip_pred_probs[k] = vp9_read_prob(r);

  mb_mode_mv_init(pbi, r);
}

void vp9_decode_mb_mode_mv(VP9D_COMP* const pbi,
                           MACROBLOCKD* const xd,
                           int mi_row,
                           int mi_col,
                           vp9_reader *r) {
  VP9_COMMON *const cm = &pbi->common;
  MODE_INFO *mi = xd->mode_info_context;
  MB_MODE_INFO *const mbmi = &mi->mbmi;

  if (cm->frame_type == KEY_FRAME) {
    kfread_modes(pbi, mi, mi_row, mi_col, r);
  } else {
    read_mb_modes_mv(pbi, mi, &mi->mbmi, mi_row, mi_col, r);
    set_scale_factors(xd,
                      mi->mbmi.ref_frame - 1, mi->mbmi.second_ref_frame - 1,
                      cm->active_ref_scale);
  }

  if (1) {
    const int bw = 1 << mi_width_log2(mbmi->sb_type);
    const int bh = 1 << mi_height_log2(mbmi->sb_type);
    const int y_mis = MIN(bh, cm->mi_rows - mi_row);
    const int x_mis = MIN(bw, cm->mi_cols - mi_col);
    const int mis = cm->mode_info_stride;
    int x, y;

    for (y = 0; y < y_mis; y++)
      for (x = !y; x < x_mis; x++)
        mi[y * mis + x] = *mi;
  }
}
