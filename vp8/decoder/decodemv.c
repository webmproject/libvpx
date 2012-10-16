/*
  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */


#include "treereader.h"
#include "vp8/common/entropymv.h"
#include "vp8/common/entropymode.h"
#include "onyxd_int.h"
#include "vp8/common/findnearmv.h"

#include "vp8/common/seg_common.h"
#include "vp8/common/pred_common.h"
#include "vp8/common/entropy.h"

#if CONFIG_DEBUG
#include <assert.h>
#endif

// #define DEBUG_DEC_MV
#ifdef DEBUG_DEC_MV
int dec_mvcount = 0;
#endif

static int vp8_read_bmode(vp8_reader *bc, const vp8_prob *p) {
  return vp8_treed_read(bc, vp8_bmode_tree, p);
}


static int vp8_read_ymode(vp8_reader *bc, const vp8_prob *p) {
  return vp8_treed_read(bc, vp8_ymode_tree, p);
}

#if CONFIG_SUPERBLOCKS
static int vp8_sb_kfread_ymode(vp8_reader *bc, const vp8_prob *p) {
  return vp8_treed_read(bc, vp8_uv_mode_tree, p);
}
#endif

static int vp8_kfread_ymode(vp8_reader *bc, const vp8_prob *p) {
  return vp8_treed_read(bc, vp8_kf_ymode_tree, p);
}

static int vp8_read_i8x8_mode(vp8_reader *bc, const vp8_prob *p) {
  return vp8_treed_read(bc, vp8_i8x8_mode_tree, p);
}


static int vp8_read_uv_mode(vp8_reader *bc, const vp8_prob *p) {
  return vp8_treed_read(bc, vp8_uv_mode_tree, p);
}

// This function reads the current macro block's segnent id from the bitstream
// It should only be called if a segment map update is indicated.
static void vp8_read_mb_segid(vp8_reader *r, MB_MODE_INFO *mi,
                              MACROBLOCKD *xd) {
  /* Is segmentation enabled */
  if (xd->segmentation_enabled && xd->update_mb_segmentation_map) {
    /* If so then read the segment id. */
    if (vp8_read(r, xd->mb_segment_tree_probs[0]))
      mi->segment_id =
        (unsigned char)(2 + vp8_read(r, xd->mb_segment_tree_probs[2]));
    else
      mi->segment_id =
        (unsigned char)(vp8_read(r, xd->mb_segment_tree_probs[1]));
  }
}

extern const int vp8_i8x8_block[4];
static void vp8_kfread_modes(VP8D_COMP *pbi,
                             MODE_INFO *m,
                             int mb_row,
                             int mb_col) {
  VP8_COMMON *const cm = & pbi->common;
  vp8_reader *const bc = & pbi->bc;
  const int mis = pbi->common.mode_info_stride;
  int map_index = mb_row * pbi->common.mb_cols + mb_col;
  MB_PREDICTION_MODE y_mode;

  // Read the Macroblock segmentation map if it is being updated explicitly
  // this frame (reset to 0 by default).
  m->mbmi.segment_id = 0;
  if (pbi->mb.update_mb_segmentation_map) {
    vp8_read_mb_segid(bc, &m->mbmi, &pbi->mb);
    pbi->common.last_frame_seg_map[map_index] = m->mbmi.segment_id;
  }

  m->mbmi.mb_skip_coeff = 0;
  if (pbi->common.mb_no_coeff_skip &&
      (!segfeature_active(&pbi->mb,
                          m->mbmi.segment_id, SEG_LVL_EOB) ||
       (get_segdata(&pbi->mb,
                    m->mbmi.segment_id, SEG_LVL_EOB) != 0))) {
    MACROBLOCKD *const xd  = & pbi->mb;
    m->mbmi.mb_skip_coeff = vp8_read(bc, get_pred_prob(cm, xd, PRED_MBSKIP));
  } else {
    if (segfeature_active(&pbi->mb,
                          m->mbmi.segment_id, SEG_LVL_EOB) &&
        (get_segdata(&pbi->mb,
                     m->mbmi.segment_id, SEG_LVL_EOB) == 0)) {
      m->mbmi.mb_skip_coeff = 1;
    } else
      m->mbmi.mb_skip_coeff = 0;
  }

#if CONFIG_SUPERBLOCKS
  if (m->mbmi.encoded_as_sb) {
    y_mode = (MB_PREDICTION_MODE) vp8_sb_kfread_ymode(bc,
      pbi->common.sb_kf_ymode_prob[pbi->common.kf_ymode_probs_index]);
  } else
#endif
  y_mode = (MB_PREDICTION_MODE) vp8_kfread_ymode(bc,
    pbi->common.kf_ymode_prob[pbi->common.kf_ymode_probs_index]);
#if CONFIG_COMP_INTRA_PRED
  m->mbmi.second_mode = (MB_PREDICTION_MODE)(DC_PRED - 1);
#endif

  m->mbmi.ref_frame = INTRA_FRAME;

  if ((m->mbmi.mode = y_mode) == B_PRED) {
    int i = 0;
#if CONFIG_COMP_INTRA_PRED
    int use_comp_pred = vp8_read(bc, 128);
#endif
    do {
      const B_PREDICTION_MODE A = above_block_mode(m, i, mis);
      const B_PREDICTION_MODE L = left_block_mode(m, i);

      m->bmi[i].as_mode.first =
        (B_PREDICTION_MODE) vp8_read_bmode(
          bc, pbi->common.kf_bmode_prob [A] [L]);
#if CONFIG_COMP_INTRA_PRED
      if (use_comp_pred) {
        m->bmi[i].as_mode.second =
          (B_PREDICTION_MODE) vp8_read_bmode(
            bc, pbi->common.kf_bmode_prob [A] [L]);
      } else {
        m->bmi[i].as_mode.second = (B_PREDICTION_MODE)(B_DC_PRED - 1);
      }
#endif
    } while (++i < 16);
  }
  if ((m->mbmi.mode = y_mode) == I8X8_PRED) {
    int i;
    int mode8x8;
    for (i = 0; i < 4; i++) {
      int ib = vp8_i8x8_block[i];
      mode8x8 = vp8_read_i8x8_mode(bc, pbi->common.fc.i8x8_mode_prob);
      m->bmi[ib + 0].as_mode.first = mode8x8;
      m->bmi[ib + 1].as_mode.first = mode8x8;
      m->bmi[ib + 4].as_mode.first = mode8x8;
      m->bmi[ib + 5].as_mode.first = mode8x8;
#if CONFIG_COMP_INTRA_PRED
      m->bmi[ib + 0].as_mode.second = (MB_PREDICTION_MODE)(DC_PRED - 1);
      m->bmi[ib + 1].as_mode.second = (MB_PREDICTION_MODE)(DC_PRED - 1);
      m->bmi[ib + 4].as_mode.second = (MB_PREDICTION_MODE)(DC_PRED - 1);
      m->bmi[ib + 5].as_mode.second = (MB_PREDICTION_MODE)(DC_PRED - 1);
#endif
    }
  } else
    m->mbmi.uv_mode = (MB_PREDICTION_MODE)vp8_read_uv_mode(bc,
                                                           pbi->common.kf_uv_mode_prob[m->mbmi.mode]);
#if CONFIG_COMP_INTRA_PRED
  m->mbmi.second_uv_mode = (MB_PREDICTION_MODE)(DC_PRED - 1);
#endif

#if CONFIG_TX_SELECT
  if (cm->txfm_mode == TX_MODE_SELECT && m->mbmi.mb_skip_coeff == 0 &&
      m->mbmi.mode <= I8X8_PRED) {
    // FIXME(rbultje) code ternary symbol once all experiments are merged
    m->mbmi.txfm_size = vp8_read(bc, cm->prob_tx[0]);
    if (m->mbmi.txfm_size != TX_4X4 && m->mbmi.mode != I8X8_PRED)
      m->mbmi.txfm_size += vp8_read(bc, cm->prob_tx[1]);
  } else
#endif
  if (cm->txfm_mode >= ALLOW_16X16 && m->mbmi.mode <= TM_PRED) {
    m->mbmi.txfm_size = TX_16X16;
  } else if (cm->txfm_mode >= ALLOW_8X8 && m->mbmi.mode != B_PRED) {
    m->mbmi.txfm_size = TX_8X8;
  } else {
    m->mbmi.txfm_size = TX_4X4;
  }
}

#if CONFIG_NEWMVENTROPY
static int read_nmv_component(vp8_reader *r,
                              int rv,
                              const nmv_component *mvcomp) {
  int v, s, z, c, o, d;
  s = vp8_read(r, mvcomp->sign);
  c = vp8_treed_read(r, vp8_mv_class_tree, mvcomp->classes);
  if (c == MV_CLASS_0) {
    d = vp8_treed_read(r, vp8_mv_class0_tree, mvcomp->class0);
  } else {
    int i, b;
    d = 0;
    b = c + CLASS0_BITS - 1;  /* number of bits */
    for (i = 0; i < b; ++i)
      d |= (vp8_read(r, mvcomp->bits[i]) << i);
  }
  o = d << 3;

  z = vp8_get_mv_mag(c, o);
  v = (s ? -(z + 1) : (z + 1));
  return v;
}

static int read_nmv_component_fp(vp8_reader *r,
                                 int v,
                                 int rv,
                                 const nmv_component *mvcomp,
                                 int usehp) {
  int s, z, c, o, d, e, f;
  s = v < 0;
  z = (s ? -v : v) - 1;       /* magnitude - 1 */

  c = vp8_get_mv_class(z, &o);
  d = o >> 3;

  if (c == MV_CLASS_0) {
    f = vp8_treed_read(r, vp8_mv_fp_tree, mvcomp->class0_fp[d]);
  } else {
    f = vp8_treed_read(r, vp8_mv_fp_tree, mvcomp->fp);
  }
  o += (f << 1);

  if (usehp) {
    if (c == MV_CLASS_0) {
      e = vp8_read(r, mvcomp->class0_hp);
    } else {
      e = vp8_read(r, mvcomp->hp);
    }
    o += e;
  } else {
    ++o;  /* Note if hp is not used, the default value of the hp bit is 1 */
  }
  z = vp8_get_mv_mag(c, o);
  v = (s ? -(z + 1) : (z + 1));
  return v;
}

static void read_nmv(vp8_reader *r, MV *mv, const MV *ref,
                     const nmv_context *mvctx) {
  MV_JOINT_TYPE j = vp8_treed_read(r, vp8_mv_joint_tree, mvctx->joints);
  mv->row = mv-> col = 0;
  if (j == MV_JOINT_HZVNZ || j == MV_JOINT_HNZVNZ) {
    mv->row = read_nmv_component(r, ref->row, &mvctx->comps[0]);
  }
  if (j == MV_JOINT_HNZVZ || j == MV_JOINT_HNZVNZ) {
    mv->col = read_nmv_component(r, ref->col, &mvctx->comps[1]);
  }
}

static void read_nmv_fp(vp8_reader *r, MV *mv, const MV *ref,
                        const nmv_context *mvctx, int usehp) {
  MV_JOINT_TYPE j = vp8_get_mv_joint(*mv);
  usehp = usehp && vp8_use_nmv_hp(ref);
  if (j == MV_JOINT_HZVNZ || j == MV_JOINT_HNZVNZ) {
    mv->row = read_nmv_component_fp(r, mv->row, ref->row, &mvctx->comps[0],
                                    usehp);
  }
  if (j == MV_JOINT_HNZVZ || j == MV_JOINT_HNZVNZ) {
    mv->col = read_nmv_component_fp(r, mv->col, ref->col, &mvctx->comps[1],
                                    usehp);
  }
  //printf("  %d: %d %d ref: %d %d\n", usehp, mv->row, mv-> col, ref->row, ref->col);
}

static void update_nmv(vp8_reader *bc, vp8_prob *const p,
                       const vp8_prob upd_p) {
  if (vp8_read(bc, upd_p)) {
#ifdef LOW_PRECISION_MV_UPDATE
    *p = (vp8_read_literal(bc, 7) << 1) | 1;
#else
    *p = (vp8_read_literal(bc, 8));
#endif
  }
}

static void read_nmvprobs(vp8_reader *bc, nmv_context *mvctx,
                          int usehp) {
  int i, j, k;
#ifdef MV_GROUP_UPDATE
  if (!vp8_read_bit(bc)) return;
#endif
  for (j = 0; j < MV_JOINTS - 1; ++j) {
    update_nmv(bc, &mvctx->joints[j],
               VP8_NMV_UPDATE_PROB);
  }
  for (i = 0; i < 2; ++i) {
    update_nmv(bc, &mvctx->comps[i].sign,
               VP8_NMV_UPDATE_PROB);
    for (j = 0; j < MV_CLASSES - 1; ++j) {
      update_nmv(bc, &mvctx->comps[i].classes[j],
                 VP8_NMV_UPDATE_PROB);
    }
    for (j = 0; j < CLASS0_SIZE - 1; ++j) {
      update_nmv(bc, &mvctx->comps[i].class0[j],
                 VP8_NMV_UPDATE_PROB);
    }
    for (j = 0; j < MV_OFFSET_BITS; ++j) {
      update_nmv(bc, &mvctx->comps[i].bits[j],
                 VP8_NMV_UPDATE_PROB);
    }
  }

  for (i = 0; i < 2; ++i) {
    for (j = 0; j < CLASS0_SIZE; ++j) {
      for (k = 0; k < 3; ++k)
        update_nmv(bc, &mvctx->comps[i].class0_fp[j][k],
                   VP8_NMV_UPDATE_PROB);
    }
    for (j = 0; j < 3; ++j) {
      update_nmv(bc, &mvctx->comps[i].fp[j],
                 VP8_NMV_UPDATE_PROB);
    }
  }

  if (usehp) {
    for (i = 0; i < 2; ++i) {
      update_nmv(bc, &mvctx->comps[i].class0_hp,
                 VP8_NMV_UPDATE_PROB);
      update_nmv(bc, &mvctx->comps[i].hp,
                 VP8_NMV_UPDATE_PROB);
    }
  }
}

#else

static int read_mvcomponent(vp8_reader *r, const MV_CONTEXT *mvc) {
  const vp8_prob *const p = (const vp8_prob *) mvc;
  int x = 0;

  if (vp8_read(r, p [mvpis_short])) { /* Large */
    int i = 0;

    do {
      x += vp8_read(r, p [MVPbits + i]) << i;
    } while (++i < mvnum_short_bits);

    i = mvlong_width - 1;  /* Skip bit 3, which is sometimes implicit */

    do {
      x += vp8_read(r, p [MVPbits + i]) << i;
    } while (--i > mvnum_short_bits);

    if (!(x & ~((2 << mvnum_short_bits) - 1))  ||  vp8_read(r, p [MVPbits + mvnum_short_bits]))
      x += (mvnum_short);
  } else /* small */
    x = vp8_treed_read(r, vp8_small_mvtree, p + MVPshort);

  if (x  &&  vp8_read(r, p [MVPsign]))
    x = -x;

  return x;
}

static void read_mv(vp8_reader *r, MV *mv, const MV_CONTEXT *mvc) {
  mv->row = (short)(read_mvcomponent(r,   mvc) << 1);
  mv->col = (short)(read_mvcomponent(r, ++mvc) << 1);
#ifdef DEBUG_DEC_MV
  int i;
  printf("%d (np): %d %d\n", dec_mvcount++, mv->row, mv->col);
  // for (i=0; i<MVPcount;++i) printf("  %d", (&mvc[-1])->prob[i]); printf("\n");
  // for (i=0; i<MVPcount;++i) printf("  %d", (&mvc[0])->prob[i]); printf("\n");
#endif
}

static void read_mvcontexts(vp8_reader *bc, MV_CONTEXT *mvc) {
  int i = 0;

  do {
    const vp8_prob *up = vp8_mv_update_probs[i].prob;
    vp8_prob *p = (vp8_prob *)(mvc + i);
    vp8_prob *const pstop = p + MVPcount;

    do {
      if (vp8_read(bc, *up++)) {
        const vp8_prob x = (vp8_prob)vp8_read_literal(bc, 7);

        *p = x ? x << 1 : 1;
      }
    } while (++p < pstop);
  } while (++i < 2);
}

static int read_mvcomponent_hp(vp8_reader *r, const MV_CONTEXT_HP *mvc) {
  const vp8_prob *const p = (const vp8_prob *) mvc;
  int x = 0;

  if (vp8_read(r, p [mvpis_short_hp])) { /* Large */
    int i = 0;

    do {
      x += vp8_read(r, p [MVPbits_hp + i]) << i;
    } while (++i < mvnum_short_bits_hp);

    i = mvlong_width_hp - 1;  /* Skip bit 3, which is sometimes implicit */

    do {
      x += vp8_read(r, p [MVPbits_hp + i]) << i;
    } while (--i > mvnum_short_bits_hp);

    if (!(x & ~((2 << mvnum_short_bits_hp) - 1))  ||  vp8_read(r, p [MVPbits_hp + mvnum_short_bits_hp]))
      x += (mvnum_short_hp);
  } else /* small */
    x = vp8_treed_read(r, vp8_small_mvtree_hp, p + MVPshort_hp);

  if (x  &&  vp8_read(r, p [MVPsign_hp]))
    x = -x;

  return x;
}

static void read_mv_hp(vp8_reader *r, MV *mv, const MV_CONTEXT_HP *mvc) {
  mv->row = (short)(read_mvcomponent_hp(r,   mvc));
  mv->col = (short)(read_mvcomponent_hp(r, ++mvc));
#ifdef DEBUG_DEC_MV
  int i;
  printf("%d (hp): %d %d\n", dec_mvcount++, mv->row, mv->col);
  // for (i=0; i<MVPcount_hp;++i) printf("  %d", (&mvc[-1])->prob[i]); printf("\n");
  // for (i=0; i<MVPcount_hp;++i) printf("  %d", (&mvc[0])->prob[i]); printf("\n");
#endif
}

static void read_mvcontexts_hp(vp8_reader *bc, MV_CONTEXT_HP *mvc) {
  int i = 0;

  do {
    const vp8_prob *up = vp8_mv_update_probs_hp[i].prob;
    vp8_prob *p = (vp8_prob *)(mvc + i);
    vp8_prob *const pstop = p + MVPcount_hp;

    do {
      if (vp8_read(bc, *up++)) {
        const vp8_prob x = (vp8_prob)vp8_read_literal(bc, 7);

        *p = x ? x << 1 : 1;
      }
    } while (++p < pstop);
  } while (++i < 2);
}

#endif  /* CONFIG_NEWMVENTROPY */

// Read the referncence frame
static MV_REFERENCE_FRAME read_ref_frame(VP8D_COMP *pbi,
                                         vp8_reader *const bc,
                                         unsigned char segment_id) {
  MV_REFERENCE_FRAME ref_frame;
  int seg_ref_active;
  int seg_ref_count = 0;

  VP8_COMMON *const cm = & pbi->common;
  MACROBLOCKD *const xd = &pbi->mb;

  seg_ref_active = segfeature_active(xd,
                                     segment_id,
                                     SEG_LVL_REF_FRAME);

  // If segment coding enabled does the segment allow for more than one
  // possible reference frame
  if (seg_ref_active) {
    seg_ref_count = check_segref(xd, segment_id, INTRA_FRAME) +
                    check_segref(xd, segment_id, LAST_FRAME) +
                    check_segref(xd, segment_id, GOLDEN_FRAME) +
                    check_segref(xd, segment_id, ALTREF_FRAME);
  }

  // Segment reference frame features not available or allows for
  // multiple reference frame options
  if (!seg_ref_active || (seg_ref_count > 1)) {
    // Values used in prediction model coding
    unsigned char prediction_flag;
    vp8_prob pred_prob;
    MV_REFERENCE_FRAME pred_ref;

    // Get the context probability the prediction flag
    pred_prob = get_pred_prob(cm, xd, PRED_REF);

    // Read the prediction status flag
    prediction_flag = (unsigned char)vp8_read(bc, pred_prob);

    // Store the prediction flag.
    set_pred_flag(xd, PRED_REF, prediction_flag);

    // Get the predicted reference frame.
    pred_ref = get_pred_ref(cm, xd);

    // If correctly predicted then use the predicted value
    if (prediction_flag) {
      ref_frame = pred_ref;
    }
    // else decode the explicitly coded value
    else {
      vp8_prob mod_refprobs[PREDICTION_PROBS];
      vpx_memcpy(mod_refprobs,
                 cm->mod_refprobs[pred_ref], sizeof(mod_refprobs));

      // If segment coding enabled blank out options that cant occur by
      // setting the branch probability to 0.
      if (seg_ref_active) {
        mod_refprobs[INTRA_FRAME] *=
          check_segref(xd, segment_id, INTRA_FRAME);
        mod_refprobs[LAST_FRAME] *=
          check_segref(xd, segment_id, LAST_FRAME);
        mod_refprobs[GOLDEN_FRAME] *=
          (check_segref(xd, segment_id, GOLDEN_FRAME) *
           check_segref(xd, segment_id, ALTREF_FRAME));
      }

      // Default to INTRA_FRAME (value 0)
      ref_frame = INTRA_FRAME;

      // Do we need to decode the Intra/Inter branch
      if (mod_refprobs[0])
        ref_frame = (MV_REFERENCE_FRAME) vp8_read(bc, mod_refprobs[0]);
      else
        ref_frame++;

      if (ref_frame) {
        // Do we need to decode the Last/Gf_Arf branch
        if (mod_refprobs[1])
          ref_frame += vp8_read(bc, mod_refprobs[1]);
        else
          ref_frame++;

        if (ref_frame > 1) {
          // Do we need to decode the GF/Arf branch
          if (mod_refprobs[2])
            ref_frame += vp8_read(bc, mod_refprobs[2]);
          else {
            if (seg_ref_active) {
              if ((pred_ref == GOLDEN_FRAME) ||
                  !check_segref(xd, segment_id, GOLDEN_FRAME)) {
                ref_frame = ALTREF_FRAME;
              } else
                ref_frame = GOLDEN_FRAME;
            } else
              ref_frame = (pred_ref == GOLDEN_FRAME)
                          ? ALTREF_FRAME : GOLDEN_FRAME;
          }
        }
      }
    }
  }

  // Segment reference frame features are enabled
  else {
    // The reference frame for the mb is considered as correclty predicted
    // if it is signaled at the segment level for the purposes of the
    // common prediction model
    set_pred_flag(xd, PRED_REF, 1);
    ref_frame = get_pred_ref(cm, xd);
  }

  return (MV_REFERENCE_FRAME)ref_frame;
}

#if CONFIG_SUPERBLOCKS
static MB_PREDICTION_MODE read_sb_mv_ref(vp8_reader *bc, const vp8_prob *p) {
  return (MB_PREDICTION_MODE) vp8_treed_read(bc, vp8_sb_mv_ref_tree, p);
}
#endif

static MB_PREDICTION_MODE read_mv_ref(vp8_reader *bc, const vp8_prob *p) {
  return (MB_PREDICTION_MODE) vp8_treed_read(bc, vp8_mv_ref_tree, p);
}

static B_PREDICTION_MODE sub_mv_ref(vp8_reader *bc, const vp8_prob *p) {
  return (B_PREDICTION_MODE) vp8_treed_read(bc, vp8_sub_mv_ref_tree, p);
}

#ifdef VPX_MODE_COUNT
unsigned int vp8_mv_cont_count[5][4] = {
  { 0, 0, 0, 0 },
  { 0, 0, 0, 0 },
  { 0, 0, 0, 0 },
  { 0, 0, 0, 0 },
  { 0, 0, 0, 0 }
};
#endif

static const unsigned char mbsplit_fill_count[4] = {8, 8, 4, 1};
static const unsigned char mbsplit_fill_offset[4][16] = {
  { 0,  1,  2,  3,  4,  5,  6,  7,  8,  9,  10, 11, 12, 13, 14, 15},
  { 0,  1,  4,  5,  8,  9, 12, 13,  2,  3,   6,  7, 10, 11, 14, 15},
  { 0,  1,  4,  5,  2,  3,  6,  7,  8,  9,  12, 13, 10, 11, 14, 15},
  { 0,  1,  2,  3,  4,  5,  6,  7,  8,  9,  10, 11, 12, 13, 14, 15}
};

#if CONFIG_SWITCHABLE_INTERP
static void read_switchable_interp_probs(VP8D_COMP *pbi) {
  VP8_COMMON *const cm = & pbi->common;
  vp8_reader *const bc = & pbi->bc;
  int i, j;
  for (j = 0; j <= VP8_SWITCHABLE_FILTERS; ++j) {
  //for (j = 0; j <= 0; ++j) {
    for (i = 0; i < VP8_SWITCHABLE_FILTERS - 1; ++i) {
      cm->fc.switchable_interp_prob[j][i] = vp8_read_literal(bc, 8);
    }
  }
  //printf("DECODER: %d %d\n", cm->fc.switchable_interp_prob[0],
  //cm->fc.switchable_interp_prob[1]);
}
#endif

static void mb_mode_mv_init(VP8D_COMP *pbi) {
  VP8_COMMON *const cm = & pbi->common;
  vp8_reader *const bc = & pbi->bc;
#if CONFIG_NEWMVENTROPY
  nmv_context *const nmvc = &pbi->common.fc.nmvc;
#else
  MV_CONTEXT *const mvc = pbi->common.fc.mvc;
  MV_CONTEXT_HP *const mvc_hp = pbi->common.fc.mvc_hp;
#endif
  MACROBLOCKD *const xd  = & pbi->mb;

  vpx_memset(cm->mbskip_pred_probs, 0, sizeof(cm->mbskip_pred_probs));
  if (pbi->common.mb_no_coeff_skip) {
    int k;
    for (k = 0; k < MBSKIP_CONTEXTS; ++k)
      cm->mbskip_pred_probs[k] = (vp8_prob)vp8_read_literal(bc, 8);
  }

  if (cm->frame_type != KEY_FRAME) {
#if CONFIG_PRED_FILTER
    cm->pred_filter_mode = (vp8_prob)vp8_read_literal(bc, 2);

    if (cm->pred_filter_mode == 2)
      cm->prob_pred_filter_off = (vp8_prob)vp8_read_literal(bc, 8);
#endif
#if CONFIG_SWITCHABLE_INTERP
    if (cm->mcomp_filter_type == SWITCHABLE)
      read_switchable_interp_probs(pbi);
#endif
    // Decode the baseline probabilities for decoding reference frame
    cm->prob_intra_coded = (vp8_prob)vp8_read_literal(bc, 8);
    cm->prob_last_coded  = (vp8_prob)vp8_read_literal(bc, 8);
    cm->prob_gf_coded    = (vp8_prob)vp8_read_literal(bc, 8);

    // Computes a modified set of probabilities for use when reference
    // frame prediction fails.
    compute_mod_refprobs(cm);

    pbi->common.comp_pred_mode = vp8_read(bc, 128);
    if (cm->comp_pred_mode)
      cm->comp_pred_mode += vp8_read(bc, 128);
    if (cm->comp_pred_mode == HYBRID_PREDICTION) {
      int i;
      for (i = 0; i < COMP_PRED_CONTEXTS; i++)
        cm->prob_comppred[i] = (vp8_prob)vp8_read_literal(bc, 8);
    }

    if (vp8_read_bit(bc)) {
      int i = 0;

      do {
        cm->fc.ymode_prob[i] = (vp8_prob) vp8_read_literal(bc, 8);
      } while (++i < VP8_YMODES - 1);
    }
#if CONFIG_NEWMVENTROPY
    read_nmvprobs(bc, nmvc, xd->allow_high_precision_mv);
#else
    if (xd->allow_high_precision_mv)
      read_mvcontexts_hp(bc, mvc_hp);
    else
      read_mvcontexts(bc, mvc);
#endif
  }
}

// This function either reads the segment id for the current macroblock from
// the bitstream or if the value is temporally predicted asserts the predicted
// value
static void read_mb_segment_id(VP8D_COMP *pbi,
                               int mb_row, int mb_col) {
  vp8_reader *const bc = & pbi->bc;
  VP8_COMMON *const cm = & pbi->common;
  MACROBLOCKD *const xd  = & pbi->mb;
  MODE_INFO *mi = xd->mode_info_context;
  MB_MODE_INFO *mbmi = &mi->mbmi;
  int index = mb_row * pbi->common.mb_cols + mb_col;

  if (xd->segmentation_enabled) {
    if (xd->update_mb_segmentation_map) {
      // Is temporal coding of the segment id for this mb enabled.
      if (cm->temporal_update) {
        // Get the context based probability for reading the
        // prediction status flag
        vp8_prob pred_prob =
          get_pred_prob(cm, xd, PRED_SEG_ID);

        // Read the prediction status flag
        unsigned char seg_pred_flag =
          (unsigned char)vp8_read(bc, pred_prob);

        // Store the prediction flag.
        set_pred_flag(xd, PRED_SEG_ID, seg_pred_flag);

        // If the value is flagged as correctly predicted
        // then use the predicted value
        if (seg_pred_flag) {
          mbmi->segment_id = get_pred_mb_segid(cm, index);
        }
        // Else .... decode it explicitly
        else {
          vp8_read_mb_segid(bc, mbmi, xd);
        }
      }
      // Normal unpredicted coding mode
      else {
        vp8_read_mb_segid(bc, mbmi, xd);
      }
#if CONFIG_SUPERBLOCKS
      if (mbmi->encoded_as_sb) {
        cm->last_frame_seg_map[index] =
        cm->last_frame_seg_map[index + 1] =
        cm->last_frame_seg_map[index + cm->mb_cols] =
        cm->last_frame_seg_map[index + cm->mb_cols + 1] = mbmi->segment_id;
      } else
#endif
      {
        cm->last_frame_seg_map[index] = mbmi->segment_id;
      }
    } else {
#if CONFIG_SUPERBLOCKS
      if (mbmi->encoded_as_sb) {
        mbmi->segment_id =
              cm->last_frame_seg_map[index] &&
              cm->last_frame_seg_map[index + 1] &&
              cm->last_frame_seg_map[index + cm->mb_cols] &&
              cm->last_frame_seg_map[index + cm->mb_cols + 1];
      } else
#endif
      {
        mbmi->segment_id = cm->last_frame_seg_map[index];
      }
    }
  } else {
    // The encoder explicitly sets the segment_id to 0
    // when segmentation is disabled
    mbmi->segment_id = 0;
  }
}

static void read_mb_modes_mv(VP8D_COMP *pbi, MODE_INFO *mi, MB_MODE_INFO *mbmi,
                             MODE_INFO *prev_mi,
                             int mb_row, int mb_col) {
  VP8_COMMON *const cm = & pbi->common;
  vp8_reader *const bc = & pbi->bc;
#if CONFIG_NEWMVENTROPY
  nmv_context *const nmvc = &pbi->common.fc.nmvc;
#else
  MV_CONTEXT *const mvc = pbi->common.fc.mvc;
  MV_CONTEXT_HP *const mvc_hp = pbi->common.fc.mvc_hp;
#endif
  const int mis = pbi->common.mode_info_stride;
  MACROBLOCKD *const xd  = & pbi->mb;

  int_mv *const mv = & mbmi->mv;
  int mb_to_left_edge;
  int mb_to_right_edge;
  int mb_to_top_edge;
  int mb_to_bottom_edge;

  mb_to_top_edge = xd->mb_to_top_edge;
  mb_to_bottom_edge = xd->mb_to_bottom_edge;
  mb_to_top_edge -= LEFT_TOP_MARGIN;
  mb_to_bottom_edge += RIGHT_BOTTOM_MARGIN;
  mbmi->need_to_clamp_mvs = 0;
  mbmi->need_to_clamp_secondmv = 0;
  mbmi->second_ref_frame = 0;
  /* Distance of Mb to the various image edges.
   * These specified to 8th pel as they are always compared to MV values that are in 1/8th pel units
   */
  xd->mb_to_left_edge =
    mb_to_left_edge = -((mb_col * 16) << 3);
  mb_to_left_edge -= LEFT_TOP_MARGIN;

  xd->mb_to_right_edge =
    mb_to_right_edge = ((pbi->common.mb_cols - 1 - mb_col) * 16) << 3;
  mb_to_right_edge += RIGHT_BOTTOM_MARGIN;

  // Make sure the MACROBLOCKD mode info pointer is pointed at the
  // correct entry for the current macroblock.
  xd->mode_info_context = mi;
  xd->prev_mode_info_context = prev_mi;

  // Read the macroblock segment id.
  read_mb_segment_id(pbi, mb_row, mb_col);

  if (pbi->common.mb_no_coeff_skip &&
      (!segfeature_active(xd,
                          mbmi->segment_id, SEG_LVL_EOB) ||
       (get_segdata(xd, mbmi->segment_id, SEG_LVL_EOB) != 0))) {
    // Read the macroblock coeff skip flag if this feature is in use,
    // else default to 0
    mbmi->mb_skip_coeff = vp8_read(bc, get_pred_prob(cm, xd, PRED_MBSKIP));
  } else {
    if (segfeature_active(xd,
                          mbmi->segment_id, SEG_LVL_EOB) &&
        (get_segdata(xd, mbmi->segment_id, SEG_LVL_EOB) == 0)) {
      mbmi->mb_skip_coeff = 1;
    } else
      mbmi->mb_skip_coeff = 0;
  }

  // Read the reference frame
  mbmi->ref_frame = read_ref_frame(pbi, bc, mbmi->segment_id);

  // If reference frame is an Inter frame
  if (mbmi->ref_frame) {
    int rct[4];
    int_mv nearest, nearby, best_mv;
    int_mv nearest_second, nearby_second, best_mv_second;
    vp8_prob mv_ref_p [VP8_MVREFS - 1];

#if CONFIG_NEWBESTREFMV
    int recon_y_stride, recon_yoffset;
    int recon_uv_stride, recon_uvoffset;
#endif

    vp8_find_near_mvs(xd, mi,
                      prev_mi,
                      &nearest, &nearby, &best_mv, rct,
                      mbmi->ref_frame, cm->ref_frame_sign_bias);

#if CONFIG_NEWBESTREFMV
    {
      int ref_fb_idx;
      MV_REFERENCE_FRAME ref_frame = mbmi->ref_frame;

      /* Select the appropriate reference frame for this MB */
      if (ref_frame == LAST_FRAME)
        ref_fb_idx = cm->lst_fb_idx;
      else if (ref_frame == GOLDEN_FRAME)
        ref_fb_idx = cm->gld_fb_idx;
      else
        ref_fb_idx = cm->alt_fb_idx;

      recon_y_stride = cm->yv12_fb[ref_fb_idx].y_stride  ;
      recon_uv_stride = cm->yv12_fb[ref_fb_idx].uv_stride;

      recon_yoffset = (mb_row * recon_y_stride * 16) + (mb_col * 16);
      recon_uvoffset = (mb_row * recon_uv_stride * 8) + (mb_col * 8);

      xd->pre.y_buffer = cm->yv12_fb[ref_fb_idx].y_buffer + recon_yoffset;
      xd->pre.u_buffer = cm->yv12_fb[ref_fb_idx].u_buffer + recon_uvoffset;
      xd->pre.v_buffer = cm->yv12_fb[ref_fb_idx].v_buffer + recon_uvoffset;

      // Update stats on relative distance of chosen vector to the
      // possible best reference vectors.
      {
        find_mv_refs(xd, mi, prev_mi,
                     ref_frame, mbmi->ref_mvs[ref_frame],
                     cm->ref_frame_sign_bias );
      }

      vp8_find_best_ref_mvs(xd,
                            xd->pre.y_buffer,
                            recon_y_stride,
                            mbmi->ref_mvs[ref_frame],
                            &best_mv, &nearest, &nearby);
    }
#endif

    vp8_mv_ref_probs(&pbi->common, mv_ref_p, rct);

    // Is the segment level mode feature enabled for this segment
    if (segfeature_active(xd, mbmi->segment_id, SEG_LVL_MODE)) {
      mbmi->mode =
        get_segdata(xd, mbmi->segment_id, SEG_LVL_MODE);
    } else {
#if CONFIG_SUPERBLOCKS
      if (mbmi->encoded_as_sb) {
        mbmi->mode = read_sb_mv_ref(bc, mv_ref_p);
      } else
#endif
      mbmi->mode = read_mv_ref(bc, mv_ref_p);

      vp8_accum_mv_refs(&pbi->common, mbmi->mode, rct);
    }

#if CONFIG_PRED_FILTER
    if (mbmi->mode >= NEARESTMV && mbmi->mode < SPLITMV) {
      // Is the prediction filter enabled
      if (cm->pred_filter_mode == 2)
        mbmi->pred_filter_enabled =
          vp8_read(bc, cm->prob_pred_filter_off);
      else
        mbmi->pred_filter_enabled = cm->pred_filter_mode;
    }
#endif
#if CONFIG_SWITCHABLE_INTERP
    if (mbmi->mode >= NEARESTMV && mbmi->mode <= SPLITMV)
    {
      if (cm->mcomp_filter_type == SWITCHABLE) {
        mbmi->interp_filter = vp8_switchable_interp[
            vp8_treed_read(bc, vp8_switchable_interp_tree,
                           get_pred_probs(cm, xd, PRED_SWITCHABLE_INTERP))];
        //printf("Reading: %d\n", mbmi->interp_filter);
      } else {
        mbmi->interp_filter = cm->mcomp_filter_type;
      }
    }
#endif

    if (cm->comp_pred_mode == COMP_PREDICTION_ONLY ||
        (cm->comp_pred_mode == HYBRID_PREDICTION &&
         vp8_read(bc, get_pred_prob(cm, xd, PRED_COMP)))) {
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
#if CONFIG_NEWBESTREFMV
      if (mbmi->second_ref_frame) {
        int second_ref_fb_idx;
        /* Select the appropriate reference frame for this MB */
        if (mbmi->second_ref_frame == LAST_FRAME)
          second_ref_fb_idx = cm->lst_fb_idx;
        else if (mbmi->second_ref_frame ==
          GOLDEN_FRAME)
          second_ref_fb_idx = cm->gld_fb_idx;
        else
          second_ref_fb_idx = cm->alt_fb_idx;

        xd->second_pre.y_buffer =
          cm->yv12_fb[second_ref_fb_idx].y_buffer + recon_yoffset;
        xd->second_pre.u_buffer =
          cm->yv12_fb[second_ref_fb_idx].u_buffer + recon_uvoffset;
        xd->second_pre.v_buffer =
          cm->yv12_fb[second_ref_fb_idx].v_buffer + recon_uvoffset;
        vp8_find_near_mvs(xd, mi, prev_mi,
                          &nearest_second, &nearby_second, &best_mv_second,
                          rct,
                          mbmi->second_ref_frame,
                          cm->ref_frame_sign_bias);

        // Update stats on relative distance of chosen vector to the
        // possible best reference vectors.
        {
          MV_REFERENCE_FRAME ref_frame = mbmi->second_ref_frame;

          find_mv_refs(xd, mi, prev_mi,
                       ref_frame, mbmi->ref_mvs[ref_frame],
                       cm->ref_frame_sign_bias );
        }

        vp8_find_best_ref_mvs(xd,
                              xd->second_pre.y_buffer,
                              recon_y_stride,
                              mbmi->ref_mvs[mbmi->second_ref_frame],
                              &best_mv_second,
                              &nearest_second,
                              &nearby_second);
      }
#else
      vp8_find_near_mvs(xd, mi, prev_mi,
                        &nearest_second, &nearby_second, &best_mv_second,
                        rct,
                        mbmi->second_ref_frame,
                        pbi->common.ref_frame_sign_bias);
#endif
    } else {
      mbmi->second_ref_frame = 0;
    }

    mbmi->uv_mode = DC_PRED;
    switch (mbmi->mode) {
      case SPLITMV: {
        const int s = mbmi->partitioning =
                        vp8_treed_read(bc, vp8_mbsplit_tree, cm->fc.mbsplit_prob);
        const int num_p = vp8_mbsplit_count [s];
        int j = 0;
        cm->fc.mbsplit_counts[s]++;

        mbmi->need_to_clamp_mvs = 0;
        do { /* for each subset j */
          int_mv leftmv, abovemv, second_leftmv, second_abovemv;
          int_mv blockmv, secondmv;
          int k;  /* first block in subset j */
          int mv_contz;
          int blockmode;

          k = vp8_mbsplit_offset[s][j];

          leftmv.as_int = left_block_mv(mi, k);
          abovemv.as_int = above_block_mv(mi, k, mis);
          if (mbmi->second_ref_frame) {
            second_leftmv.as_int = left_block_second_mv(mi, k);
            second_abovemv.as_int = above_block_second_mv(mi, k, mis);
          }
          mv_contz = vp8_mv_cont(&leftmv, &abovemv);
          blockmode = sub_mv_ref(bc, cm->fc.sub_mv_ref_prob [mv_contz]);
          cm->fc.sub_mv_ref_counts[mv_contz][blockmode - LEFT4X4]++;

          switch (blockmode) {
            case NEW4X4:
#if CONFIG_NEWMVENTROPY
              read_nmv(bc, &blockmv.as_mv, &best_mv.as_mv, nmvc);
              read_nmv_fp(bc, &blockmv.as_mv, &best_mv.as_mv, nmvc,
                          xd->allow_high_precision_mv);
              vp8_increment_nmv(&blockmv.as_mv, &best_mv.as_mv,
                                &cm->fc.NMVcount, xd->allow_high_precision_mv);
#else
              if (xd->allow_high_precision_mv) {
                read_mv_hp(bc, &blockmv.as_mv, (const MV_CONTEXT_HP *) mvc_hp);
                cm->fc.MVcount_hp[0][mv_max_hp + (blockmv.as_mv.row)]++;
                cm->fc.MVcount_hp[1][mv_max_hp + (blockmv.as_mv.col)]++;
              } else {
                read_mv(bc, &blockmv.as_mv, (const MV_CONTEXT *) mvc);
                cm->fc.MVcount[0][mv_max + (blockmv.as_mv.row >> 1)]++;
                cm->fc.MVcount[1][mv_max + (blockmv.as_mv.col >> 1)]++;
              }
#endif  /* CONFIG_NEWMVENTROPY */
              blockmv.as_mv.row += best_mv.as_mv.row;
              blockmv.as_mv.col += best_mv.as_mv.col;

              if (mbmi->second_ref_frame) {
#if CONFIG_NEWMVENTROPY
                read_nmv(bc, &secondmv.as_mv, &best_mv_second.as_mv, nmvc);
                read_nmv_fp(bc, &secondmv.as_mv, &best_mv_second.as_mv, nmvc,
                            xd->allow_high_precision_mv);
                vp8_increment_nmv(&secondmv.as_mv, &best_mv_second.as_mv,
                                  &cm->fc.NMVcount, xd->allow_high_precision_mv);
#else
                if (xd->allow_high_precision_mv) {
                  read_mv_hp(bc, &secondmv.as_mv, (const MV_CONTEXT_HP *) mvc_hp);
                  cm->fc.MVcount_hp[0][mv_max_hp + (secondmv.as_mv.row)]++;
                  cm->fc.MVcount_hp[1][mv_max_hp + (secondmv.as_mv.col)]++;
                } else {
                  read_mv(bc, &secondmv.as_mv, (const MV_CONTEXT *) mvc);
                  cm->fc.MVcount[0][mv_max + (secondmv.as_mv.row >> 1)]++;
                  cm->fc.MVcount[1][mv_max + (secondmv.as_mv.col >> 1)]++;
                }
#endif  /* CONFIG_NEWMVENTROPY */
                secondmv.as_mv.row += best_mv_second.as_mv.row;
                secondmv.as_mv.col += best_mv_second.as_mv.col;
              }
#ifdef VPX_MODE_COUNT
              vp8_mv_cont_count[mv_contz][3]++;
#endif
              break;
            case LEFT4X4:
              blockmv.as_int = leftmv.as_int;
              if (mbmi->second_ref_frame)
                secondmv.as_int = second_leftmv.as_int;
#ifdef VPX_MODE_COUNT
              vp8_mv_cont_count[mv_contz][0]++;
#endif
              break;
            case ABOVE4X4:
              blockmv.as_int = abovemv.as_int;
              if (mbmi->second_ref_frame)
                secondmv.as_int = second_abovemv.as_int;
#ifdef VPX_MODE_COUNT
              vp8_mv_cont_count[mv_contz][1]++;
#endif
              break;
            case ZERO4X4:
              blockmv.as_int = 0;
              if (mbmi->second_ref_frame)
                secondmv.as_int = 0;
#ifdef VPX_MODE_COUNT
              vp8_mv_cont_count[mv_contz][2]++;
#endif
              break;
            default:
              break;
          }

          mbmi->need_to_clamp_mvs |= vp8_check_mv_bounds(&blockmv,
                                                         mb_to_left_edge,
                                                         mb_to_right_edge,
                                                         mb_to_top_edge,
                                                         mb_to_bottom_edge);
          if (mbmi->second_ref_frame) {
            mbmi->need_to_clamp_mvs |= vp8_check_mv_bounds(&secondmv,
                                                           mb_to_left_edge,
                                                           mb_to_right_edge,
                                                           mb_to_top_edge,
                                                           mb_to_bottom_edge);
          }

          {
            /* Fill (uniform) modes, mvs of jth subset.
             Must do it here because ensuing subsets can
             refer back to us via "left" or "above". */
            const unsigned char *fill_offset;
            unsigned int fill_count = mbsplit_fill_count[s];

            fill_offset = &mbsplit_fill_offset[s][(unsigned char)j * mbsplit_fill_count[s]];

            do {
              mi->bmi[ *fill_offset].as_mv.first.as_int = blockmv.as_int;
              if (mbmi->second_ref_frame)
                mi->bmi[ *fill_offset].as_mv.second.as_int = secondmv.as_int;
              fill_offset++;
            } while (--fill_count);
          }

        } while (++j < num_p);
      }

      mv->as_int = mi->bmi[15].as_mv.first.as_int;
      mbmi->mv[1].as_int = mi->bmi[15].as_mv.second.as_int;

      break;  /* done with SPLITMV */

      case NEARMV:
        mv->as_int = nearby.as_int;
        /* Clip "next_nearest" so that it does not extend to far out of image */
        vp8_clamp_mv(mv, mb_to_left_edge, mb_to_right_edge,
                     mb_to_top_edge, mb_to_bottom_edge);
        if (mbmi->second_ref_frame) {
          mbmi->mv[1].as_int = nearby_second.as_int;
          vp8_clamp_mv(&mbmi->mv[1], mb_to_left_edge, mb_to_right_edge,
                       mb_to_top_edge, mb_to_bottom_edge);
        }
        break;

      case NEARESTMV:
        mv->as_int = nearest.as_int;
        /* Clip "next_nearest" so that it does not extend to far out of image */
        vp8_clamp_mv(mv, mb_to_left_edge, mb_to_right_edge,
                     mb_to_top_edge, mb_to_bottom_edge);
        if (mbmi->second_ref_frame) {
          mbmi->mv[1].as_int = nearest_second.as_int;
          vp8_clamp_mv(&mbmi->mv[1], mb_to_left_edge, mb_to_right_edge,
                       mb_to_top_edge, mb_to_bottom_edge);
        }
        break;

      case ZEROMV:
        mv->as_int = 0;
        if (mbmi->second_ref_frame)
          mbmi->mv[1].as_int = 0;
        break;

      case NEWMV:
#if CONFIG_NEWMVENTROPY
        read_nmv(bc, &mv->as_mv, &best_mv.as_mv, nmvc);
        read_nmv_fp(bc, &mv->as_mv, &best_mv.as_mv, nmvc,
                    xd->allow_high_precision_mv);
        vp8_increment_nmv(&mv->as_mv, &best_mv.as_mv, &cm->fc.NMVcount,
                          xd->allow_high_precision_mv);
#else
        if (xd->allow_high_precision_mv) {
          read_mv_hp(bc, &mv->as_mv, (const MV_CONTEXT_HP *) mvc_hp);
          cm->fc.MVcount_hp[0][mv_max_hp + (mv->as_mv.row)]++;
          cm->fc.MVcount_hp[1][mv_max_hp + (mv->as_mv.col)]++;
        } else {
          read_mv(bc, &mv->as_mv, (const MV_CONTEXT *) mvc);
          cm->fc.MVcount[0][mv_max + (mv->as_mv.row >> 1)]++;
          cm->fc.MVcount[1][mv_max + (mv->as_mv.col >> 1)]++;
        }
#endif  /* CONFIG_NEWMVENTROPY */
        mv->as_mv.row += best_mv.as_mv.row;
        mv->as_mv.col += best_mv.as_mv.col;

        /* Don't need to check this on NEARMV and NEARESTMV modes
         * since those modes clamp the MV. The NEWMV mode does not,
         * so signal to the prediction stage whether special
         * handling may be required.
         */
        mbmi->need_to_clamp_mvs = vp8_check_mv_bounds(mv,
                                                      mb_to_left_edge,
                                                      mb_to_right_edge,
                                                      mb_to_top_edge,
                                                      mb_to_bottom_edge);
        if (mbmi->second_ref_frame) {
#if CONFIG_NEWMVENTROPY
          read_nmv(bc, &mbmi->mv[1].as_mv, &best_mv_second.as_mv, nmvc);
          read_nmv_fp(bc, &mbmi->mv[1].as_mv, &best_mv_second.as_mv, nmvc,
                      xd->allow_high_precision_mv);
          vp8_increment_nmv(&mbmi->mv[1].as_mv, &best_mv_second.as_mv,
                            &cm->fc.NMVcount, xd->allow_high_precision_mv);
#else
          if (xd->allow_high_precision_mv) {
            read_mv_hp(bc, &mbmi->mv[1].as_mv, (const MV_CONTEXT_HP *) mvc_hp);
            cm->fc.MVcount_hp[0][mv_max_hp + (mbmi->mv[1].as_mv.row)]++;
            cm->fc.MVcount_hp[1][mv_max_hp + (mbmi->mv[1].as_mv.col)]++;
          } else {
            read_mv(bc, &mbmi->mv[1].as_mv, (const MV_CONTEXT *) mvc);
            cm->fc.MVcount[0][mv_max + (mbmi->mv[1].as_mv.row >> 1)]++;
            cm->fc.MVcount[1][mv_max + (mbmi->mv[1].as_mv.col >> 1)]++;
          }
#endif  /* CONFIG_NEWMVENTROPY */
          mbmi->mv[1].as_mv.row += best_mv_second.as_mv.row;
          mbmi->mv[1].as_mv.col += best_mv_second.as_mv.col;
          mbmi->need_to_clamp_secondmv |=
            vp8_check_mv_bounds(&mbmi->mv[1],
                                mb_to_left_edge, mb_to_right_edge,
                                mb_to_top_edge, mb_to_bottom_edge);
        }
        break;
      default:
;
#if CONFIG_DEBUG
        assert(0);
#endif
    }
  } else {
    /* required for left and above block mv */
    mbmi->mv[0].as_int = 0;

    if (segfeature_active(xd, mbmi->segment_id, SEG_LVL_MODE))
      mbmi->mode = (MB_PREDICTION_MODE)
                   get_segdata(xd, mbmi->segment_id, SEG_LVL_MODE);
    else {
      // FIXME write using SB mode tree
      mbmi->mode = (MB_PREDICTION_MODE)
                   vp8_read_ymode(bc, pbi->common.fc.ymode_prob);
      pbi->common.fc.ymode_counts[mbmi->mode]++;
    }
#if CONFIG_COMP_INTRA_PRED
    mbmi->second_mode = (MB_PREDICTION_MODE)(DC_PRED - 1);
#endif

    // If MB mode is BPRED read the block modes
    if (mbmi->mode == B_PRED) {
      int j = 0;
#if CONFIG_COMP_INTRA_PRED
      int use_comp_pred = vp8_read(bc, 128);
#endif
      do {
        mi->bmi[j].as_mode.first = (B_PREDICTION_MODE)vp8_read_bmode(bc, pbi->common.fc.bmode_prob);
        pbi->common.fc.bmode_counts[mi->bmi[j].as_mode.first]++;
#if CONFIG_COMP_INTRA_PRED
        if (use_comp_pred) {
          mi->bmi[j].as_mode.second = (B_PREDICTION_MODE)vp8_read_bmode(bc, pbi->common.fc.bmode_prob);
        } else {
          mi->bmi[j].as_mode.second = (B_PREDICTION_MODE)(B_DC_PRED - 1);
        }
#endif
      } while (++j < 16);
    }

    if (mbmi->mode == I8X8_PRED) {
      int i;
      int mode8x8;
      for (i = 0; i < 4; i++) {
        int ib = vp8_i8x8_block[i];
        mode8x8 = vp8_read_i8x8_mode(bc, pbi->common.fc.i8x8_mode_prob);
        mi->bmi[ib + 0].as_mode.first = mode8x8;
        mi->bmi[ib + 1].as_mode.first = mode8x8;
        mi->bmi[ib + 4].as_mode.first = mode8x8;
        mi->bmi[ib + 5].as_mode.first = mode8x8;
        pbi->common.fc.i8x8_mode_counts[mode8x8]++;
#if CONFIG_COMP_INTRA_PRED
        mi->bmi[ib + 0].as_mode.second = (MB_PREDICTION_MODE)(DC_PRED - 1);
        mi->bmi[ib + 1].as_mode.second = (MB_PREDICTION_MODE)(DC_PRED - 1);
        mi->bmi[ib + 4].as_mode.second = (MB_PREDICTION_MODE)(DC_PRED - 1);
        mi->bmi[ib + 5].as_mode.second = (MB_PREDICTION_MODE)(DC_PRED - 1);
#endif
      }
    } else {
      mbmi->uv_mode = (MB_PREDICTION_MODE)vp8_read_uv_mode(
        bc, pbi->common.fc.uv_mode_prob[mbmi->mode]);
      pbi->common.fc.uv_mode_counts[mbmi->mode][mbmi->uv_mode]++;
    }

#if CONFIG_COMP_INTRA_PRED
    mbmi->second_uv_mode = (MB_PREDICTION_MODE)(DC_PRED - 1);
#endif
  }

#if CONFIG_TX_SELECT
  if (cm->txfm_mode == TX_MODE_SELECT && mbmi->mb_skip_coeff == 0 &&
      ((mbmi->ref_frame == INTRA_FRAME && mbmi->mode <= I8X8_PRED) ||
       (mbmi->ref_frame != INTRA_FRAME && mbmi->mode != SPLITMV))) {
    // FIXME(rbultje) code ternary symbol once all experiments are merged
    mbmi->txfm_size = vp8_read(bc, cm->prob_tx[0]);
    if (mbmi->txfm_size != TX_4X4 && mbmi->mode != I8X8_PRED)
      mbmi->txfm_size += vp8_read(bc, cm->prob_tx[1]);
  } else
#endif
  if (cm->txfm_mode >= ALLOW_16X16 &&
      ((mbmi->ref_frame == INTRA_FRAME && mbmi->mode <= TM_PRED) ||
       (mbmi->ref_frame != INTRA_FRAME && mbmi->mode != SPLITMV))) {
    mbmi->txfm_size = TX_16X16;
  } else if (cm->txfm_mode >= ALLOW_8X8 &&
      ((mbmi->ref_frame == INTRA_FRAME && mbmi->mode != B_PRED) ||
       (mbmi->ref_frame != INTRA_FRAME && mbmi->mode != SPLITMV))) {
    mbmi->txfm_size = TX_8X8;
  } else {
    mbmi->txfm_size = TX_4X4;
  }
}

void vpx_decode_mode_mvs_init(VP8D_COMP *pbi){
  VP8_COMMON *cm = &pbi->common;
  mb_mode_mv_init(pbi);
  if (cm->frame_type == KEY_FRAME && !cm->kf_ymode_probs_update)
    cm->kf_ymode_probs_index = vp8_read_literal(&pbi->bc, 3);
}
void vpx_decode_mb_mode_mv(VP8D_COMP *pbi,
                           MACROBLOCKD *xd,
                           int mb_row,
                           int mb_col){
  MODE_INFO *mi = xd->mode_info_context;
  MODE_INFO *prev_mi = xd->prev_mode_info_context;

  if (pbi->common.frame_type == KEY_FRAME)
    vp8_kfread_modes(pbi, mi, mb_row, mb_col);
  else
    read_mb_modes_mv(pbi, mi, &mi->mbmi, prev_mi, mb_row, mb_col);
}
