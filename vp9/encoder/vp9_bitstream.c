/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include <assert.h>
#include <stdio.h>
#include <limits.h>

#include "vpx/vpx_encoder.h"
#include "vpx_mem/vpx_mem.h"

#include "vp9/common/vp9_entropy.h"
#include "vp9/common/vp9_entropymode.h"
#include "vp9/common/vp9_entropymv.h"
#include "vp9/common/vp9_mvref_common.h"
#if CONFIG_PALETTE
#include "vp9/common/vp9_palette.h"
#endif  // CONFIG_PALETTE
#if CONFIG_SR_MODE
#include "vp9/common/vp9_sr_txfm.h"
#endif  // CONFIG_SR_MODE
#include "vp9/common/vp9_pred_common.h"
#include "vp9/common/vp9_seg_common.h"
#include "vp9/common/vp9_systemdependent.h"
#include "vp9/common/vp9_tile_common.h"

#include "vp9/encoder/vp9_cost.h"
#include "vp9/encoder/vp9_bitstream.h"
#include "vp9/encoder/vp9_encodemv.h"
#include "vp9/encoder/vp9_mcomp.h"
#include "vp9/encoder/vp9_segmentation.h"
#include "vp9/encoder/vp9_subexp.h"
#include "vp9/encoder/vp9_tokenize.h"
#include "vp9/encoder/vp9_write_bit_buffer.h"

static struct vp9_token intra_mode_encodings[INTRA_MODES];
static struct vp9_token switchable_interp_encodings[SWITCHABLE_FILTERS];
#if CONFIG_EXT_PARTITION
static struct vp9_token ext_partition_encodings[EXT_PARTITION_TYPES];
#endif
static struct vp9_token partition_encodings[PARTITION_TYPES];
static struct vp9_token inter_mode_encodings[INTER_MODES];
#if CONFIG_SR_MODE && SR_USE_MULTI_F
static struct vp9_token sr_usfilter_encodings[SR_USFILTER_NUM];
#endif  // CONFIG_SR_MODE && SR_USE_MULTI_F
#if CONFIG_EXT_TX
static struct vp9_token ext_tx_encodings[EXT_TX_TYPES];
#if CONFIG_WAVELETS
static struct vp9_token ext_tx_large_encodings[EXT_TX_TYPES_LARGE];
#endif  // CONFIG_WAVELETS
#endif  // CONFIG_EXT_TX
#if CONFIG_PALETTE
static struct vp9_token palette_size_encodings[PALETTE_SIZES];
static struct vp9_token palette_color_encodings[PALETTE_COLORS];
#endif  // CONFIG_PALETTE
#if CONFIG_COPY_MODE
static struct vp9_token copy_mode_encodings_l2[2];
static struct vp9_token copy_mode_encodings[COPY_MODE_COUNT - 1];
#endif
#if CONFIG_NEW_INTER
static struct vp9_token inter_compound_mode_encodings[INTER_COMPOUND_MODES];
#endif  // CONFIG_NEW_INTER
#if CONFIG_GLOBAL_MOTION
static struct vp9_token global_motion_types_encodings[GLOBAL_MOTION_TYPES];
#endif  // CONFIG_GLOBAL_MOTION
#if CONFIG_NEW_QUANT && QUANT_PROFILES > 1 && !Q_CTX_BASED_PROFILES
static struct vp9_token dq_profile_encodings[QUANT_PROFILES];
#endif  // CONFIG_NEW_QUANT && QUANT_PROFILES > 1 && !Q_CTX_BASED_PROFILES

#if CONFIG_SUPERTX
static int vp9_check_supertx(VP9_COMMON *cm, int mi_row, int mi_col,
                             BLOCK_SIZE bsize) {
  MODE_INFO *mi;

  mi = cm->mi + (mi_row * cm->mi_stride + mi_col);

  return mi[0].mbmi.tx_size == bsize_to_tx_size(bsize) &&
         mi[0].mbmi.sb_type < bsize;
}
#endif  // CONFIG_SUPERTX

void vp9_entropy_mode_init() {
  vp9_tokens_from_tree(intra_mode_encodings, vp9_intra_mode_tree);
  vp9_tokens_from_tree(switchable_interp_encodings, vp9_switchable_interp_tree);
  vp9_tokens_from_tree(partition_encodings, vp9_partition_tree);
#if CONFIG_EXT_PARTITION
  vp9_tokens_from_tree(ext_partition_encodings, vp9_ext_partition_tree);
#endif
  vp9_tokens_from_tree(inter_mode_encodings, vp9_inter_mode_tree);
#if CONFIG_SR_MODE && SR_USE_MULTI_F
  vp9_tokens_from_tree(sr_usfilter_encodings, vp9_sr_usfilter_tree);
#endif  // CONFIG_SR_MODE && SR_USE_MULTI_F
#if CONFIG_EXT_TX
  vp9_tokens_from_tree(ext_tx_encodings, vp9_ext_tx_tree);
#if CONFIG_WAVELETS
  vp9_tokens_from_tree(ext_tx_large_encodings, vp9_ext_tx_large_tree);
#endif
#endif  // CONFIG_EXT_TX
#if CONFIG_PALETTE
  vp9_tokens_from_tree(palette_size_encodings, vp9_palette_size_tree);
  vp9_tokens_from_tree(palette_color_encodings, vp9_palette_color_tree);
#endif  // CONFIG_PALETTE
#if CONFIG_NEW_INTER
  vp9_tokens_from_tree(inter_compound_mode_encodings,
                       vp9_inter_compound_mode_tree);
#endif  // CONFIG_NEW_INTER
#if CONFIG_COPY_MODE
  vp9_tokens_from_tree(copy_mode_encodings_l2, vp9_copy_mode_tree_l2);
  vp9_tokens_from_tree(copy_mode_encodings, vp9_copy_mode_tree);
#endif  // CONFIG_COPY_MODE
#if CONFIG_GLOBAL_MOTION
  vp9_tokens_from_tree(global_motion_types_encodings,
                       vp9_global_motion_types_tree);
#endif  // CONFIG_GLOBAL_MOTION
#if CONFIG_NEW_QUANT && QUANT_PROFILES > 1 && !Q_CTX_BASED_PROFILES
  vp9_tokens_from_tree(dq_profile_encodings, vp9_dq_profile_tree);
#endif  // CONFIG_NEW_QUANT && QUANT_PROFILES > 1 && !Q_CTX_BASED_PROFILES
}

static void write_intra_mode(vp9_writer *w, PREDICTION_MODE mode,
                             const vp9_prob *probs) {
#if CONFIG_INTRABC
  assert(!is_intrabc_mode(mode));
#endif  // CONFIG_INTRABC
  vp9_write_token(w, vp9_intra_mode_tree, probs, &intra_mode_encodings[mode]);
}

static void write_inter_mode(vp9_writer *w, PREDICTION_MODE mode,
                             const vp9_prob *probs) {
  assert(is_inter_mode(mode));
  vp9_write_token(w, vp9_inter_mode_tree, probs,
                  &inter_mode_encodings[INTER_OFFSET(mode)]);
}

#if CONFIG_COPY_MODE
static void write_copy_mode(VP9_COMMON *cm, vp9_writer *w, COPY_MODE mode,
                            int inter_ref_count, int copy_mode_context) {
  if (inter_ref_count == 2) {
    vp9_write_token(w, vp9_copy_mode_tree_l2,
                    cm->fc.copy_mode_probs_l2[copy_mode_context],
                    &copy_mode_encodings_l2[mode - REF0]);
  } else if (inter_ref_count > 2) {
    vp9_write_token(w, vp9_copy_mode_tree,
                    cm->fc.copy_mode_probs[copy_mode_context],
                    &copy_mode_encodings[mode - REF0]);
  }
}
#endif  // CONFIG_COPY_MODE

#if CONFIG_NEW_INTER
static void write_inter_compound_mode(vp9_writer *w, PREDICTION_MODE mode,
                                      const vp9_prob *probs) {
  assert(is_inter_compound_mode(mode));
  vp9_write_token(w, vp9_inter_compound_mode_tree, probs,
                  &inter_compound_mode_encodings[INTER_COMPOUND_OFFSET(mode)]);
}
#endif  // CONFIG_NEW_INTER

static void encode_unsigned_max(struct vp9_write_bit_buffer *wb,
                                int data, int max) {
  vp9_wb_write_literal(wb, data, get_unsigned_bits(max));
}

static void prob_diff_update(const vp9_tree_index *tree,
                             vp9_prob probs[/*n - 1*/],
                             const unsigned int counts[/*n - 1*/],
                             int n, vp9_writer *w) {
  int i;
  unsigned int branch_ct[32][2];

  // Assuming max number of probabilities <= 32
  assert(n <= 32);

  vp9_tree_probs_from_distribution(tree, branch_ct, counts);
  for (i = 0; i < n - 1; ++i)
    vp9_cond_prob_diff_update(w, &probs[i], branch_ct[i]);
}

static int prob_diff_update_savings(const vp9_tree_index *tree,
                                    vp9_prob probs[/*n - 1*/],
                                    const unsigned int counts[/*n - 1*/],
                                    int n) {
  int i;
  unsigned int branch_ct[32][2];
  int savings = 0;

  // Assuming max number of probabilities <= 32
  assert(n <= 32);

  vp9_tree_probs_from_distribution(tree, branch_ct, counts);
  for (i = 0; i < n - 1; ++i) {
    savings += vp9_cond_prob_diff_update_savings(&probs[i], branch_ct[i]);
  }
  return savings;
}

static void write_selected_tx_size(const VP9_COMMON *cm,
                                   const MACROBLOCKD *xd,
                                   TX_SIZE tx_size, BLOCK_SIZE bsize,
                                   vp9_writer *w) {
  const TX_SIZE max_tx_size = max_txsize_lookup[bsize];
  const vp9_prob *const tx_probs = get_tx_probs2(max_tx_size, xd,
                                                 &cm->fc.tx_probs);
  vp9_write(w, tx_size != TX_4X4, tx_probs[0]);
  if (tx_size != TX_4X4 && max_tx_size >= TX_16X16) {
    vp9_write(w, tx_size != TX_8X8, tx_probs[1]);
    if (tx_size != TX_8X8 && max_tx_size >= TX_32X32) {
      vp9_write(w, tx_size != TX_16X16, tx_probs[2]);
#if CONFIG_TX64X64
      if (tx_size != TX_16X16 && max_tx_size >= TX_64X64)
        vp9_write(w, tx_size != TX_32X32, tx_probs[3]);
#endif
    }
  }
}
#if CONFIG_MISC_ENTROPY
static int write_skip(const VP9_COMMON *cm, const MACROBLOCKD *xd,
                      int segment_id, const MODE_INFO *mi, int is_inter,
                      vp9_writer *w) {
#else
  static int write_skip(const VP9_COMMON *cm, const MACROBLOCKD *xd,
                        int segment_id, const MODE_INFO *mi, vp9_writer *w) {
#endif
  if (vp9_segfeature_active(&cm->seg, segment_id, SEG_LVL_SKIP)) {
    return 1;
#if CONFIG_MISC_ENTROPY
  } else if (!is_inter && cm->frame_type == INTER_FRAME) {
    return 0;
#endif
  } else {
    const int skip = mi->mbmi.skip;
    vp9_write(w, skip, vp9_get_skip_prob(cm, xd));
    return skip;
  }
}

static void update_skip_probs(VP9_COMMON *cm, vp9_writer *w) {
  int k;
  for (k = 0; k < SKIP_CONTEXTS; ++k)
    vp9_cond_prob_diff_update(w, &cm->fc.skip_probs[k], cm->counts.skip[k]);
}

#if CONFIG_NEW_QUANT && QUANT_PROFILES > 1 && !Q_CTX_BASED_PROFILES
static void write_dq_profile(const VP9_COMMON *cm, int dq_profile,
                             vp9_writer *w) {
  vp9_write_token(w, vp9_dq_profile_tree, cm->fc.dq_profile_prob,
                  &dq_profile_encodings[dq_profile]);
}

static void update_dq_profile_probs(VP9_COMMON *cm, vp9_writer *w) {
  prob_diff_update(vp9_dq_profile_tree,
                   cm->fc.dq_profile_prob,
                   cm->counts.dq_profile, QUANT_PROFILES, w);
}
#endif  // CONFIG_NEW_QUANT && QUANT_PROFILES > 1 && !Q_CTX_BASED_PROFILES

#if CONFIG_SR_MODE
#if SR_USE_MULTI_F
static int write_sr_usfilter(const VP9_COMMON *cm, const MACROBLOCKD *xd,
                             const MODE_INFO *mi, vp9_writer *w) {
  const int f_idx = mi->mbmi.us_filter_idx;
  assert(f_idx >= 0 && f_idx < SR_USFILTER_NUM);
  vp9_write_token(w, vp9_sr_usfilter_tree, vp9_get_sr_usfilter_prob(cm, xd),
                  &sr_usfilter_encodings[f_idx]);
  return f_idx;
}
#endif  // SR_USE_MULTI_F

static int write_sr(const VP9_COMMON *cm, const MACROBLOCKD *xd,
                      int segment_id, const MODE_INFO *mi, vp9_writer *w) {
  if (vp9_segfeature_active(&cm->seg, segment_id, SEG_LVL_SKIP))
    return 1;
  else if (is_enable_srmode(mi->mbmi.sb_type)) {
    const int sr = mi->mbmi.sr;
    vp9_write(w, sr, vp9_get_sr_prob(cm, xd, mi->mbmi.sb_type));
#if SR_USE_MULTI_F
    if (sr) {
      write_sr_usfilter(cm, xd, mi, w);
    }
#endif  // SR_USE_MULTI_F
    return sr;
  } else {
    assert(!mi->mbmi.sr);
    return 0;
  }
}

static void update_sr_probs(VP9_COMMON *cm, vp9_writer *w) {
  int k;
  for (k = 0; k < SR_CONTEXTS; ++k)
    vp9_cond_prob_diff_update(w, &cm->fc.sr_probs[k], cm->counts.sr[k]);
}

#if SR_USE_MULTI_F
static void update_sr_usfilter_probs(VP9_COMMON *cm, vp9_writer *w) {
  int i;
  for (i = 0; i < SR_USFILTER_CONTEXTS; ++i) {
    prob_diff_update(vp9_sr_usfilter_tree, cm->fc.sr_usfilter_probs[i],
                     cm->counts.sr_usfilters[i], SR_USFILTER_NUM, w);
  }
}
#endif  // SR_USE_MULTI_F
#endif  // CONFIG_SR_MODE

static void update_switchable_interp_probs(VP9_COMMON *cm, vp9_writer *w) {
  int j;
  for (j = 0; j < SWITCHABLE_FILTER_CONTEXTS; ++j)
    prob_diff_update(vp9_switchable_interp_tree,
                     cm->fc.switchable_interp_prob[j],
                     cm->counts.switchable_interp[j], SWITCHABLE_FILTERS, w);
}

#if CONFIG_EXT_TX
static void update_ext_tx_probs(VP9_COMMON *cm, vp9_writer *w) {
  const int savings_thresh = vp9_cost_one(GROUP_DIFF_UPDATE_PROB) -
                             vp9_cost_zero(GROUP_DIFF_UPDATE_PROB);
  int i;
  int savings = 0;
  int do_update = 0;
  for (i = TX_4X4; i <= TX_16X16; ++i) {
    savings += prob_diff_update_savings(vp9_ext_tx_tree, cm->fc.ext_tx_prob[i],
                                        cm->counts.ext_tx[i], EXT_TX_TYPES);
  }
#if CONFIG_WAVELETS
  for (; i < TX_SIZES; ++i) {
    savings += prob_diff_update_savings(
        vp9_ext_tx_large_tree, cm->fc.ext_tx_prob[i], cm->counts.ext_tx[i],
        EXT_TX_TYPES_LARGE);
  }
#endif
  do_update = savings > savings_thresh;
  vp9_write(w, do_update, GROUP_DIFF_UPDATE_PROB);
  if (do_update) {
    for (i = TX_4X4; i <= TX_16X16; ++i) {
      prob_diff_update(vp9_ext_tx_tree, cm->fc.ext_tx_prob[i],
                       cm->counts.ext_tx[i], EXT_TX_TYPES, w);
    }
#if CONFIG_WAVELETS
    for (; i < TX_SIZES; ++i) {
      prob_diff_update(vp9_ext_tx_large_tree, cm->fc.ext_tx_prob[i],
                       cm->counts.ext_tx[i], EXT_TX_TYPES_LARGE, w);
    }
#endif  // CONFIG_WAVELETS
  }
}
#endif  // CONFIG_EXT_TX

#if CONFIG_SUPERTX
static void update_supertx_probs(VP9_COMMON *cm, vp9_writer *w) {
  const int savings_thresh = vp9_cost_one(GROUP_DIFF_UPDATE_PROB) -
                             vp9_cost_zero(GROUP_DIFF_UPDATE_PROB);
  int i, j;
  int savings = 0;
  int do_update = 0;
  for (i = 0; i < PARTITION_SUPERTX_CONTEXTS; ++i) {
    for (j = 1; j < TX_SIZES; ++j) {
      savings +=
          vp9_cond_prob_diff_update_savings(&cm->fc.supertx_prob[i][j],
                                            cm->counts.supertx[i][j]);
    }
  }
  do_update = savings > savings_thresh;
  vp9_write(w, do_update, GROUP_DIFF_UPDATE_PROB);
  if (do_update) {
    for (i = 0; i < PARTITION_SUPERTX_CONTEXTS; ++i) {
      for (j = 1; j < TX_SIZES; ++j) {
        vp9_cond_prob_diff_update(w, &cm->fc.supertx_prob[i][j],
                                  cm->counts.supertx[i][j]);
      }
    }
  }
}
#endif  // CONFIG_SUPERTX

#if CONFIG_NEW_INTER
static void update_inter_compound_mode_probs(VP9_COMMON *cm, vp9_writer *w) {
  const int savings_thresh = vp9_cost_one(GROUP_DIFF_UPDATE_PROB) -
                             vp9_cost_zero(GROUP_DIFF_UPDATE_PROB);
  int i;
  int savings = 0;
  int do_update = 0;
  for (i = 0; i < INTER_MODE_CONTEXTS; ++i) {
    savings += prob_diff_update_savings(vp9_inter_compound_mode_tree,
                                        cm->fc.inter_compound_mode_probs[i],
                                        cm->counts.inter_compound_mode[i],
                                        INTER_COMPOUND_MODES);
  }
  do_update = savings > savings_thresh;
  vp9_write(w, do_update, GROUP_DIFF_UPDATE_PROB);
  if (do_update) {
    for (i = 0; i < INTER_MODE_CONTEXTS; ++i) {
      prob_diff_update(vp9_inter_compound_mode_tree,
                       cm->fc.inter_compound_mode_probs[i],
                       cm->counts.inter_compound_mode[i],
                       INTER_COMPOUND_MODES, w);
    }
  }
}
#endif  // CONFIG_NEW_INTER

static void pack_mb_tokens(vp9_writer *w,
                           TOKENEXTRA **tp, const TOKENEXTRA *const stop,
                           vpx_bit_depth_t bit_depth) {
  TOKENEXTRA *p = *tp;

  while (p < stop && p->token != EOSB_TOKEN) {
    const int t = p->token;
    const struct vp9_token *const a = &vp9_coef_encodings[t];
    int i = 0;
    int v = a->value;
    int n = a->len;
#if CONFIG_VP9_HIGHBITDEPTH
    const vp9_extra_bit *b;
    if (bit_depth == VPX_BITS_12)
      b = &vp9_extra_bits_high12[t];
    else if (bit_depth == VPX_BITS_10)
      b = &vp9_extra_bits_high10[t];
    else
      b = &vp9_extra_bits[t];
#else
    const vp9_extra_bit *const b = &vp9_extra_bits[t];
    (void) bit_depth;
#endif  // CONFIG_VP9_HIGHBITDEPTH

    /* skip one or two nodes */
    if (p->skip_eob_node) {
      n -= p->skip_eob_node;
      i = 2 * p->skip_eob_node;
    }

    // TODO(jbb): expanding this can lead to big gains.  It allows
    // much better branch prediction and would enable us to avoid numerous
    // lookups and compares.

    // If we have a token that's in the constrained set, the coefficient tree
    // is split into two treed writes.  The first treed write takes care of the
    // unconstrained nodes.  The second treed write takes care of the
    // constrained nodes.
#if CONFIG_TX_SKIP
    if (p->is_pxd_token && FOR_SCREEN_CONTENT) {
      vp9_write_tree(w, vp9_coef_tree, p->context_tree, v, n, i);
    } else {
#endif  // CONFIG_TX_SKIP
      if (t >= TWO_TOKEN && t < EOB_TOKEN) {
        int len = UNCONSTRAINED_NODES - p->skip_eob_node;
        int bits = v >> (n - len);
        vp9_write_tree(w, vp9_coef_tree, p->context_tree, bits, len, i);
        vp9_write_tree(w, vp9_coef_con_tree,
                       vp9_pareto8_full[p->context_tree[PIVOT_NODE] - 1],
                       v, n - len, 0);
      } else {
        vp9_write_tree(w, vp9_coef_tree, p->context_tree, v, n, i);
      }
#if CONFIG_TX_SKIP
    }
#endif  // CONFIG_TX_SKIP

    if (b->base_val) {
      const int e = p->extra, l = b->len;

      if (l) {
        const unsigned char *pb = b->prob;
        int v = e >> 1;
        int n = l;              /* number of bits in v, assumed nonzero */
        int i = 0;

        do {
          const int bb = (v >> --n) & 1;
          vp9_write(w, bb, pb[i >> 1]);
          i = b->tree[i + bb];
        } while (n);
      }

      vp9_write_bit(w, e & 1);
    }
    ++p;
  }

  *tp = p + (p->token == EOSB_TOKEN);
}

static void write_segment_id(vp9_writer *w, const struct segmentation *seg,
                             int segment_id) {
  if (seg->enabled && seg->update_map)
    vp9_write_tree(w, vp9_segment_tree, seg->tree_probs, segment_id, 3, 0);
}

// This function encodes the reference frame
static void write_ref_frames(const VP9_COMMON *cm, const MACROBLOCKD *xd,
                             vp9_writer *w) {
  const MB_MODE_INFO *const mbmi = &xd->mi[0].src_mi->mbmi;
  const int is_compound = has_second_ref(mbmi);
  const int segment_id = mbmi->segment_id;

  // If segment level coding of this signal is disabled...
  // or the segment allows multiple reference frame options
  if (vp9_segfeature_active(&cm->seg, segment_id, SEG_LVL_REF_FRAME)) {
    assert(!is_compound);
    assert(mbmi->ref_frame[0] ==
               vp9_get_segdata(&cm->seg, segment_id, SEG_LVL_REF_FRAME));
  } else {
    // does the feature use compound prediction or not
    // (if not specified at the frame/segment level)
    if (cm->reference_mode == REFERENCE_MODE_SELECT) {
      vp9_write(w, is_compound, vp9_get_reference_mode_prob(cm, xd));
    } else {
      assert(!is_compound == (cm->reference_mode == SINGLE_REFERENCE));
    }

    if (is_compound) {
#if CONFIG_MULTI_REF
      const int bit = (mbmi->ref_frame[0] == GOLDEN_FRAME ||
                       mbmi->ref_frame[0] == LAST3_FRAME ||
                       mbmi->ref_frame[0] == LAST4_FRAME);
#else
      const int bit = mbmi->ref_frame[0] == GOLDEN_FRAME;
#endif  // CONFIG_MULTI_REF
      vp9_write(w, bit, vp9_get_pred_prob_comp_ref_p(cm, xd));

#if CONFIG_MULTI_REF
      if (!bit) {
        const int bit1 = mbmi->ref_frame[0] == LAST_FRAME;
        vp9_write(w, bit1, vp9_get_pred_prob_comp_ref_p1(cm, xd));
      } else {
        const int bit2 = mbmi->ref_frame[0] == GOLDEN_FRAME;
        vp9_write(w, bit2, vp9_get_pred_prob_comp_ref_p2(cm, xd));

        if (!bit2) {
          const int bit3 = mbmi->ref_frame[0] == LAST3_FRAME;
          vp9_write(w, bit3, vp9_get_pred_prob_comp_ref_p3(cm, xd));
        }
      }
#endif  // CONFIG_MULTI_REF
    } else {
#if CONFIG_MULTI_REF
      const int bit0 = (mbmi->ref_frame[0] == GOLDEN_FRAME ||
                        mbmi->ref_frame[0] == ALTREF_FRAME);
      vp9_write(w, bit0, vp9_get_pred_prob_single_ref_p1(cm, xd));

      if (bit0) {
        const int bit1 = mbmi->ref_frame[0] != GOLDEN_FRAME;
        vp9_write(w, bit1, vp9_get_pred_prob_single_ref_p2(cm, xd));
      } else {
        const int bit2 = (mbmi->ref_frame[0] == LAST3_FRAME ||
                          mbmi->ref_frame[0] == LAST4_FRAME);
        vp9_write(w, bit2, vp9_get_pred_prob_single_ref_p3(cm, xd));

        if (!bit2) {
          const int bit3 = mbmi->ref_frame[0] != LAST_FRAME;
          vp9_write(w, bit3, vp9_get_pred_prob_single_ref_p4(cm, xd));
        } else {
          const int bit4 = mbmi->ref_frame[0] != LAST3_FRAME;
          vp9_write(w, bit4, vp9_get_pred_prob_single_ref_p5(cm, xd));
        }
      }
#else  // CONFIG_MULTI_REF
      const int bit0 = mbmi->ref_frame[0] != LAST_FRAME;
      vp9_write(w, bit0, vp9_get_pred_prob_single_ref_p1(cm, xd));
      if (bit0) {
        const int bit1 = mbmi->ref_frame[0] != GOLDEN_FRAME;
        vp9_write(w, bit1, vp9_get_pred_prob_single_ref_p2(cm, xd));
      }
#endif  // CONFIG_MULTI_REF
    }
  }
}

static void pack_inter_mode_mvs(VP9_COMP *cpi, const MODE_INFO *mi,
#if CONFIG_SUPERTX
                                int supertx_enabled,
#endif
                                vp9_writer *w) {
  VP9_COMMON *const cm = &cpi->common;
  const nmv_context *nmvc = &cm->fc.nmvc;
  const MACROBLOCK *const x = &cpi->mb;
  const MACROBLOCKD *const xd = &x->e_mbd;
  const struct segmentation *const seg = &cm->seg;
  const MB_MODE_INFO *const mbmi = &mi->mbmi;
  const PREDICTION_MODE mode = mbmi->mode;
  const int segment_id = mbmi->segment_id;
  const BLOCK_SIZE bsize = mbmi->sb_type;
  const int allow_hp = cm->allow_high_precision_mv;
  const int is_inter = is_inter_block(mbmi);
  const int is_compound = has_second_ref(mbmi);
  int skip, ref;
#if CONFIG_COPY_MODE
  int copy_mode_context = vp9_get_copy_mode_context(xd);
  if (bsize >= BLOCK_8X8 && mbmi->inter_ref_count > 0) {
    vp9_write(w, mbmi->copy_mode != NOREF,
              cm->fc.copy_noref_prob[copy_mode_context][bsize]);
    if (mbmi->copy_mode != NOREF)
      write_copy_mode(cm, w, mbmi->copy_mode, mbmi->inter_ref_count,
                      copy_mode_context);
  }
#endif

  if (seg->update_map) {
    if (seg->temporal_update) {
      const int pred_flag = mbmi->seg_id_predicted;
      vp9_prob pred_prob = vp9_get_pred_prob_seg_id(seg, xd);
      vp9_write(w, pred_flag, pred_prob);
      if (!pred_flag)
        write_segment_id(w, seg, segment_id);
    } else {
      write_segment_id(w, seg, segment_id);
    }
  }

#if !CONFIG_MISC_ENTROPY
#if CONFIG_SUPERTX
  if (supertx_enabled)
    skip = mbmi->skip;
  else
    skip = write_skip(cm, xd, segment_id, mi, w);
#else
  skip = write_skip(cm, xd, segment_id, mi, w);
#endif  // CONFIG_SUPERTX
#endif

#if CONFIG_SUPERTX
  if (!supertx_enabled) {
#endif
#if CONFIG_COPY_MODE
    if (mbmi->copy_mode == NOREF)
#endif
      if (!vp9_segfeature_active(seg, segment_id, SEG_LVL_REF_FRAME))
        vp9_write(w, is_inter, vp9_get_intra_inter_prob(cm, xd));
#if CONFIG_MISC_ENTROPY
    skip = write_skip(cm, xd, segment_id, mi, is_inter, w);
#endif

#if CONFIG_SUPERTX
  } else {
    skip = mbmi->skip;
  }
#endif

#if CONFIG_PALETTE
  if (!is_inter && bsize >= BLOCK_8X8 && cm->allow_palette_mode) {
    int n, i, j, k, rows, cols, palette_ctx, color_ctx;
    int color_new_idx = -1, color_order[PALETTE_MAX_SIZE];
    uint8_t buffer[CODING_UNIT_SIZE * CODING_UNIT_SIZE];
    const MODE_INFO *above_mi = xd->up_available ?
        xd->mi[-xd->mi_stride].src_mi : NULL;
    const MODE_INFO *left_mi = xd->left_available ?
        xd->mi[-1].src_mi : NULL;

    palette_ctx = 0;
    if (above_mi)
      palette_ctx += (above_mi->mbmi.palette_enabled[0] == 1);
    if (left_mi)
      palette_ctx += (left_mi->mbmi.palette_enabled[0] == 1);
    vp9_write(w, mbmi->palette_enabled[0],
              cm->fc.palette_enabled_prob[bsize - BLOCK_8X8][palette_ctx]);
    vp9_write(w, mbmi->palette_enabled[1],
              cm->fc.palette_uv_enabled_prob[mbmi->palette_enabled[0]]);

    if (mbmi->palette_enabled[0]) {
      rows = 4 * num_4x4_blocks_high_lookup[bsize];
      cols = 4 * num_4x4_blocks_wide_lookup[bsize];
      n = mbmi->palette_size[0];

      vp9_write_token(w, vp9_palette_size_tree,
                      cm->fc.palette_size_prob[bsize - BLOCK_8X8],
                      &palette_size_encodings[n - 2]);
      for (i = 0; i < n; i++)
        vp9_write_literal(w, mbmi->palette_colors[i], cm->bit_depth);

      memcpy(buffer, mbmi->palette_color_map,
             rows * cols * sizeof(buffer[0]));
      vp9_write_literal(w, buffer[0], vp9_ceil_log2(n));
      for (i = 0; i < rows; i++) {
        for (j = (i == 0 ? 1 : 0); j < cols; j++) {
          color_ctx = vp9_get_palette_color_context(buffer, cols, i, j, n,
                                                    color_order);
          for (k = 0; k < n; k++)
            if (buffer[i * cols + j] == color_order[k]) {
              color_new_idx = k;
              break;
            }
          vp9_write_token(w, vp9_palette_color_tree,
                          cm->fc.palette_color_prob[n - 2][color_ctx],
                          &palette_color_encodings[color_new_idx]);
        }
      }
    }

    if (mbmi->palette_enabled[1]) {
      rows = 4 * num_4x4_blocks_high_lookup[bsize] >>
          xd->plane[1].subsampling_y;
      cols = 4 * num_4x4_blocks_wide_lookup[bsize] >>
          xd->plane[1].subsampling_x;
      n = mbmi->palette_size[1];

      if (xd->plane[1].subsampling_x && xd->plane[1].subsampling_y) {
        vp9_write_token(w, vp9_palette_size_tree,
                        cm->fc.palette_uv_size_prob[bsize - BLOCK_8X8],
                        &palette_size_encodings[n - 2]);
      }

      for (i = 0; i < n; i++)
        vp9_write_literal(w, mbmi->palette_colors[PALETTE_MAX_SIZE + i],
                          cm->bit_depth);
      for (i = 0; i < n; i++)
        vp9_write_literal(w, mbmi->palette_colors[2 * PALETTE_MAX_SIZE + i],
                          cm->bit_depth);

      if (xd->plane[1].subsampling_x && xd->plane[1].subsampling_y) {
        memcpy(buffer, mbmi->palette_uv_color_map,
               rows * cols * sizeof(buffer[0]));
        vp9_write_literal(w, buffer[0], vp9_ceil_log2(n));
        for (i = 0; i < rows; i++) {
          for (j = (i == 0 ? 1 : 0); j < cols; j++) {
            color_ctx = vp9_get_palette_color_context(buffer, cols, i, j, n,
                                                      color_order);
            for (k = 0; k < n; k++)
              if (buffer[i * cols + j] == color_order[k]) {
                color_new_idx = k;
                break;
              }
            vp9_write_token(w, vp9_palette_color_tree,
                            cm->fc.palette_uv_color_prob[n - 2][color_ctx],
                            &palette_color_encodings[color_new_idx]);
          }
        }
      }
    }
  }
#endif
#if CONFIG_SR_MODE
  if (!(is_inter && skip)) {
    write_sr(cm, xd, segment_id, mi, w);
  }
  if (mbmi->sr && !mbmi->skip) {
    assert(mbmi->tx_size == max_txsize_lookup[bsize]);
  }
#endif  // CONFIG_SR_MODE
  if (bsize >= BLOCK_8X8 && cm->tx_mode == TX_MODE_SELECT &&
#if CONFIG_SUPERTX
      !supertx_enabled &&
#endif  // CONFIG_SUPERTX
#if CONFIG_PALETTE
      !mbmi->palette_enabled[0] &&
#endif  // CONFIG_PALETTE
#if CONFIG_SR_MODE
      !mi->mbmi.sr &&
#endif  // CONFIG_SR_MODE
      !(is_inter &&
        (skip || vp9_segfeature_active(seg, segment_id, SEG_LVL_SKIP)))) {
    write_selected_tx_size(cm, xd, mbmi->tx_size, bsize, w);
  }

#if CONFIG_NEW_QUANT && QUANT_PROFILES > 1 && !Q_CTX_BASED_PROFILES
  if (cm->base_qindex > Q_THRESHOLD_MIN && cm->base_qindex < Q_THRESHOLD_MAX &&
      mbmi->send_dq_bit &&
#if CONFIG_SUPERTX
      !supertx_enabled &&
#endif  // CONFIG_SUPERTX
      !mbmi->skip &&
      !vp9_segfeature_active(&cm->seg, mbmi->segment_id, SEG_LVL_SKIP)) {
#if CONFIG_COPY_MODE
    if (mbmi->copy_mode == NOREF)
#endif  // CONFIG_COPY_MODE

      write_dq_profile(cm, mbmi->dq_off_index, w);
  } else {
    assert(mbmi->dq_off_index == 0);
  }
#endif  // CONFIG_NEW_QUANT && QUANT_PROFILES > 1 && !Q_CTX_BASED_PROFILES

#if CONFIG_EXT_TX
  if (is_inter &&
#if !CONFIG_WAVELETS
      mbmi->tx_size <= TX_16X16 &&
#endif
      cm->base_qindex > 0 &&
      bsize >= BLOCK_8X8 &&
#if CONFIG_SUPERTX
      !supertx_enabled &&
#endif
      !mbmi->skip &&
      !vp9_segfeature_active(&cm->seg, mbmi->segment_id, SEG_LVL_SKIP)) {
#if CONFIG_WAVELETS
    vp9_write_token(w, GET_EXT_TX_TREE(mbmi->tx_size),
                    cm->fc.ext_tx_prob[mbmi->tx_size],
                    &((GET_EXT_TX_ENCODINGS(mbmi->tx_size))[mbmi->ext_txfrm]));
#else
    vp9_write_token(w, vp9_ext_tx_tree, cm->fc.ext_tx_prob[mbmi->tx_size],
                    &ext_tx_encodings[mbmi->ext_txfrm]);
#endif  // CONFIG_WAVELETS
  }
#endif  // CONFIG_EXT_TX

#if CONFIG_TX_SKIP
  if (bsize >= BLOCK_8X8) {
    int q_idx = vp9_get_qindex(seg, segment_id, cm->base_qindex);
    int try_tx_skip = is_inter ? q_idx <= tx_skip_q_thresh_inter :
                                 q_idx <= tx_skip_q_thresh_intra;

#if CONFIG_COPY_MODE
    if (mbmi->copy_mode != NOREF) {
      try_tx_skip = 0;
    }
#endif  // CONFIG_COPY_MODE

#if CONFIG_SUPERTX
    if (try_tx_skip && !supertx_enabled) {
#else
    if (try_tx_skip && (!skip || !is_inter)) {
#endif  // CONFIG_SUPERTX
      if (xd->lossless) {
#if CONFIG_SUPERTX
        if (1)
#else
        if (mbmi->tx_size == TX_4X4)
#endif  // CONFIG_SUPERTX
          vp9_write(w, mbmi->tx_skip[0], cm->fc.y_tx_skip_prob[is_inter]);
#if CONFIG_SUPERTX
        if (1)
#else
        if (get_uv_tx_size(mbmi, &xd->plane[1]) == TX_4X4)
#endif  // CONFIG_SUPERTX
          vp9_write(w, mbmi->tx_skip[1],
                    cm->fc.uv_tx_skip_prob[mbmi->tx_skip[0]]);
      } else {
        vp9_write(w, mbmi->tx_skip[0], cm->fc.y_tx_skip_prob[is_inter]);
        vp9_write(w, mbmi->tx_skip[1],
                  cm->fc.uv_tx_skip_prob[mbmi->tx_skip[0]]);
      }
    }
  }
#endif  // CONFIG_TX_SKIP

  if (!is_inter) {
    if (bsize >= BLOCK_8X8) {
#if CONFIG_PALETTE
      if (!mbmi->palette_enabled[0])
        write_intra_mode(w, mode, cm->fc.y_mode_prob[size_group_lookup[bsize]]);
#else
      write_intra_mode(w, mode, cm->fc.y_mode_prob[size_group_lookup[bsize]]);
#endif  // CONFIG_PALETTE
#if CONFIG_FILTERINTRA
      if (is_filter_allowed(mode) && is_filter_enabled(mbmi->tx_size)
#if CONFIG_PALETTE
          && !mbmi->palette_enabled[0]
#endif  // CONFIG_PALETTE
      ) {
        vp9_write(w, mbmi->filterbit,
                  cm->fc.filterintra_prob[mbmi->tx_size][mode]);
      }
#endif  // CONFIG_FILTERINTRA
    } else {
      int idx, idy;
      const int num_4x4_w = num_4x4_blocks_wide_lookup[bsize];
      const int num_4x4_h = num_4x4_blocks_high_lookup[bsize];
      for (idy = 0; idy < 2; idy += num_4x4_h) {
        for (idx = 0; idx < 2; idx += num_4x4_w) {
          const PREDICTION_MODE b_mode = mi->bmi[idy * 2 + idx].as_mode;
          write_intra_mode(w, b_mode, cm->fc.y_mode_prob[0]);
#if CONFIG_FILTERINTRA
          if (is_filter_allowed(b_mode)) {
            vp9_write(w, mi->b_filter_info[idy * 2 + idx],
                      cm->fc.filterintra_prob[0][b_mode]);
          }
#endif  // CONFIG_FILTERINTRA
        }
      }
    }
#if CONFIG_PALETTE
    if (!mbmi->palette_enabled[1])
#endif  // CONFIG_PALETTE
    write_intra_mode(w, mbmi->uv_mode, cm->fc.uv_mode_prob[mode]);
#if CONFIG_FILTERINTRA
    if (is_filter_allowed(mbmi->uv_mode) &&
        is_filter_enabled(get_uv_tx_size(mbmi, &xd->plane[1]))
#if CONFIG_PALETTE
        && !mbmi->palette_enabled[1]
#endif  // CONFIG_PALETTE
    ) {
      vp9_write(w, mbmi->uv_filterbit,
                cm->fc.filterintra_prob[get_uv_tx_size(mbmi, &xd->plane[1])][mbmi->uv_mode]);
    }
#endif  // CONFIG_FILTERINTRA
#if CONFIG_COPY_MODE
  } else if (mbmi->copy_mode == NOREF) {
#else
  } else {
#endif  // CONFIG_COPY_MODE
    const int mode_ctx = mbmi->mode_context[mbmi->ref_frame[0]];
    const vp9_prob *const inter_probs = cm->fc.inter_mode_probs[mode_ctx];
#if CONFIG_NEW_INTER
    const vp9_prob *const inter_compound_probs =
        cm->fc.inter_compound_mode_probs[mode_ctx];
#endif  // CONFIG_NEW_INTER
    write_ref_frames(cm, xd, w);

    // If segment skip is not enabled code the mode.
    if (!vp9_segfeature_active(seg, segment_id, SEG_LVL_SKIP)) {
      if (bsize >= BLOCK_8X8) {
#if CONFIG_NEW_INTER
        if (is_inter_compound_mode(mode))
          write_inter_compound_mode(w, mode, inter_compound_probs);
        else if (is_inter_mode(mode))
#endif  // CONFIG_NEW_INTER
        write_inter_mode(w, mode, inter_probs);
      }
    }

    if (cm->interp_filter == SWITCHABLE) {
      const int ctx = vp9_get_pred_context_switchable_interp(xd);
      vp9_write_token(w, vp9_switchable_interp_tree,
                      cm->fc.switchable_interp_prob[ctx],
                      &switchable_interp_encodings[mbmi->interp_filter]);
      ++cpi->interp_filter_selected[0][mbmi->interp_filter];
    } else {
      assert(mbmi->interp_filter == cm->interp_filter);
    }
#if CONFIG_INTERINTRA
    if (cpi->common.reference_mode != COMPOUND_REFERENCE &&
        is_interintra_allowed(bsize) &&
        is_inter_mode(mode) &&
#if CONFIG_SUPERTX
        !supertx_enabled &&
#endif  // CONFIG_SUPERTX
        mbmi->ref_frame[1] <= INTRA_FRAME) {
      vp9_write(w, mbmi->ref_frame[1] == INTRA_FRAME,
                cm->fc.interintra_prob[bsize]);
      if (mbmi->ref_frame[1] == INTRA_FRAME) {
        write_intra_mode(w, mbmi->interintra_mode,
                         cm->fc.y_mode_prob[size_group_lookup[bsize]]);
#if CONFIG_WEDGE_PARTITION
        if (get_wedge_bits(bsize)) {
          vp9_write(w, mbmi->use_wedge_interintra,
                    cm->fc.wedge_interintra_prob[bsize]);
          if (mbmi->use_wedge_interintra) {
            vp9_write_literal(w, mbmi->interintra_wedge_index,
                              get_wedge_bits(bsize));
          }
        }
#endif  // CONFIG_WEDGE_PARTITION
      }
    }
#endif  // CONFIG_INTERINTRA

    if (bsize < BLOCK_8X8) {
      const int num_4x4_w = num_4x4_blocks_wide_lookup[bsize];
      const int num_4x4_h = num_4x4_blocks_high_lookup[bsize];
      int idx, idy;
      for (idy = 0; idy < 2; idy += num_4x4_h) {
        for (idx = 0; idx < 2; idx += num_4x4_w) {
          const int j = idy * 2 + idx;
          const PREDICTION_MODE b_mode = mi->bmi[j].as_mode;
#if CONFIG_NEW_INTER
          if (is_inter_compound_mode(b_mode))
            write_inter_compound_mode(w, b_mode, inter_compound_probs);
          else if (is_inter_mode(b_mode))
#endif  // CONFIG_NEW_INTER
          write_inter_mode(w, b_mode, inter_probs);

#if CONFIG_NEW_INTER
          if (b_mode == NEWMV || b_mode == NEW2MV ||
              b_mode == NEW_NEWMV) {
#else
          if (b_mode == NEWMV) {
#endif  // CONFIG_NEW_INTER
            for (ref = 0; ref < 1 + is_compound; ++ref) {
              vp9_encode_mv(cpi, w, &mi->bmi[j].as_mv[ref].as_mv,
#if CONFIG_NEW_INTER
                            &mi->bmi[j].ref_mv[ref].as_mv,
#else
                            &mbmi->ref_mvs[mbmi->ref_frame[ref]][0].as_mv,
#endif  // CONFIG_NEW_INTER
                            nmvc, allow_hp);
            }
          }
#if CONFIG_NEW_INTER
          else if (b_mode == NEAREST_NEWMV || b_mode == NEAR_NEWMV) {
            vp9_encode_mv(cpi, w, &mi->bmi[j].as_mv[1].as_mv,
                          &mi->bmi[j].ref_mv[1].as_mv, nmvc, allow_hp);
          } else if (b_mode == NEW_NEARESTMV || b_mode == NEW_NEARMV) {
            vp9_encode_mv(cpi, w, &mi->bmi[j].as_mv[0].as_mv,
                          &mi->bmi[j].ref_mv[0].as_mv, nmvc, allow_hp);
          }
#endif  // CONFIG_NEW_INTER
        }
      }
    } else {
#if CONFIG_NEW_INTER
      if (mode == NEWMV || mode == NEW2MV || mode == NEW_NEWMV) {
#else
      if (mode == NEWMV) {
#endif  // CONFIG_NEW_INTER
        for (ref = 0; ref < 1 + is_compound; ++ref) {
#if CONFIG_NEW_INTER
          if (mode == NEW2MV)
            vp9_encode_mv(cpi, w, &mbmi->mv[ref].as_mv,
                          &mbmi->ref_mvs[mbmi->ref_frame[ref]][1].as_mv, nmvc,
                          allow_hp);
          else
#endif  // CONFIG_NEW_INTER
          vp9_encode_mv(cpi, w, &mbmi->mv[ref].as_mv,
                        &mbmi->ref_mvs[mbmi->ref_frame[ref]][0].as_mv, nmvc,
                        allow_hp);
        }
      }
#if CONFIG_NEW_INTER
      else if (mode == NEAREST_NEWMV || mode == NEAR_NEWMV) {
        vp9_encode_mv(cpi, w, &mbmi->mv[1].as_mv,
                      &mbmi->ref_mvs[mbmi->ref_frame[1]][0].as_mv, nmvc,
                      allow_hp);
      } else if (mode == NEW_NEARESTMV || mode == NEW_NEARMV) {
        vp9_encode_mv(cpi, w, &mbmi->mv[0].as_mv,
                      &mbmi->ref_mvs[mbmi->ref_frame[0]][0].as_mv, nmvc,
                      allow_hp);
      }
#endif  // CONFIG_NEW_INTER
    }
#if CONFIG_WEDGE_PARTITION
    if (cm->reference_mode != SINGLE_REFERENCE &&
#if CONFIG_NEW_INTER
        is_inter_compound_mode(mode) &&
#endif  // CONFIG_NEW_INTER
        get_wedge_bits(bsize) &&
        mbmi->ref_frame[1] > INTRA_FRAME) {
      vp9_write(w, mbmi->use_wedge_interinter,
                cm->fc.wedge_interinter_prob[bsize]);
      if (mbmi->use_wedge_interinter)
        vp9_write_literal(w, mbmi->interinter_wedge_index,
                          get_wedge_bits(bsize));
    }
#endif  // CONFIG_NEW_INTER
  }
}

static void write_mb_modes_kf(const VP9_COMMON *cm,
#if CONFIG_PALETTE
                              MACROBLOCKD *xd,
#else
                              const MACROBLOCKD *xd,
#endif  // CONFIG_PALETTE
                              MODE_INFO *mi_8x8, vp9_writer *w) {
  const struct segmentation *const seg = &cm->seg;
  const MODE_INFO *const mi = mi_8x8;
  const MODE_INFO *const above_mi = xd->up_available ?
                                    mi_8x8[-xd->mi_stride].src_mi : NULL;
  const MODE_INFO *const left_mi =
      xd->left_available ? mi_8x8[-1].src_mi : NULL;
  const MB_MODE_INFO *const mbmi = &mi->mbmi;
  const BLOCK_SIZE bsize = mbmi->sb_type;
#if CONFIG_INTRABC
  const nmv_context *ndvc = &cm->fc.ndvc;
#endif  // CONFIG_INTRABC

  if (seg->update_map)
    write_segment_id(w, seg, mbmi->segment_id);
#if !CONFIG_MISC_ENTROPY
  write_skip(cm, xd, mbmi->segment_id, mi, w);
#endif

#if CONFIG_INTRABC
  if (bsize >= BLOCK_8X8 && cm->allow_intrabc_mode) {
    vp9_write(w, is_intrabc_mode(mbmi->mode), INTRABC_PROB);
  }
#endif  // CONFIG_INTRABC

#if CONFIG_PALETTE
  if (bsize >= BLOCK_8X8 && cm->allow_palette_mode
#if CONFIG_INTRABC
      && !is_intrabc_mode(mbmi->mode)
#endif  // CONFIG_INTRABC
      ) {
    int n, m1, m2, i, j, k, rows, cols, palette_ctx, color_ctx;
    int color_new_idx = -1, color_order[PALETTE_MAX_SIZE];
    uint8_t buffer[CODING_UNIT_SIZE * CODING_UNIT_SIZE];

    palette_ctx = 0;
    if (above_mi)
      palette_ctx += (above_mi->mbmi.palette_enabled[0] == 1);
    if (left_mi)
      palette_ctx += (left_mi->mbmi.palette_enabled[0] == 1);
    vp9_write(w, mbmi->palette_enabled[0],
              cm->fc.palette_enabled_prob[bsize - BLOCK_8X8][palette_ctx]);
    vp9_write(w, mbmi->palette_enabled[1],
              cm->fc.palette_uv_enabled_prob[mbmi->palette_enabled[0]]);

    if (mbmi->palette_enabled[0]) {
      rows = 4 * num_4x4_blocks_high_lookup[bsize];
      cols = 4 * num_4x4_blocks_wide_lookup[bsize];
      n = mbmi->palette_size[0];
      m1 = mbmi->palette_indexed_size;
      m2 = mbmi->palette_literal_size;

      vp9_write_token(w, vp9_palette_size_tree,
                      cm->fc.palette_size_prob[bsize - BLOCK_8X8],
                      &palette_size_encodings[mbmi->palette_size[0] - 2]);
      if ((xd->plane[1].subsampling_x && xd->plane[1].subsampling_y)
          || !mbmi->palette_enabled[1])
        vp9_encode_uniform(w, MIN(mbmi->palette_size[0] + 1, 8),
                           mbmi->palette_indexed_size);
      if (PALETTE_DELTA_BIT)
        vp9_write_literal(w, mbmi->palette_delta_bitdepth, PALETTE_DELTA_BIT);

      if (m1 > 0) {
        for (i = 0; i < m1; i++)
          vp9_write_literal(w, mbmi->palette_indexed_colors[i],
                            vp9_ceil_log2(mbmi->current_palette_size));
        if (mbmi->palette_delta_bitdepth > 0) {
          for (i = 0; i < m1; i++) {
            vp9_write_bit(w, mbmi->palette_color_delta[i] < 0);
            vp9_write_literal(w, abs(mbmi->palette_color_delta[i]),
                              mbmi->palette_delta_bitdepth);
          }
        }
      }
      if (m2 > 0) {
        for (i = 0; i < m2; i++)
          vp9_write_literal(w, mbmi->palette_literal_colors[i], cm->bit_depth);
      }

      memcpy(buffer, mbmi->palette_color_map,
             rows * cols * sizeof(buffer[0]));
      vp9_write_literal(w, buffer[0], vp9_ceil_log2(n));
      for (i = 0; i < rows; i++) {
        for (j = (i == 0 ? 1 : 0); j < cols; j++) {
          color_ctx = vp9_get_palette_color_context(buffer, cols, i, j, n,
                                                    color_order);
          for (k = 0; k < n; k++)
            if (buffer[i * cols + j] == color_order[k]) {
              color_new_idx = k;
              break;
            }
          vp9_write_token(w, vp9_palette_color_tree,
                          cm->fc.palette_color_prob[n - 2][color_ctx],
                          &palette_color_encodings[color_new_idx]);
        }
      }
    }

    if (mbmi->palette_enabled[1]) {
      rows = 4 * num_4x4_blocks_high_lookup[bsize] >>
          xd->plane[1].subsampling_y;
      cols = 4 * num_4x4_blocks_wide_lookup[bsize] >>
          xd->plane[1].subsampling_x;
      n = mbmi->palette_size[1];

      if (xd->plane[1].subsampling_x && xd->plane[1].subsampling_y) {
        vp9_write_token(w, vp9_palette_size_tree,
                        cm->fc.palette_uv_size_prob[bsize - BLOCK_8X8],
                        &palette_size_encodings[n - 2]);
      }

      for (i = 0; i < n; i++)
        vp9_write_literal(w, mbmi->palette_colors[PALETTE_MAX_SIZE + i],
                          cm->bit_depth);
      for (i = 0; i < n; i++)
        vp9_write_literal(w, mbmi->palette_colors[2 * PALETTE_MAX_SIZE + i],
                          cm->bit_depth);

      if (xd->plane[1].subsampling_x && xd->plane[1].subsampling_y) {
        memcpy(buffer, mbmi->palette_uv_color_map,
               rows * cols * sizeof(buffer[0]));
        vp9_write_literal(w, buffer[0], vp9_ceil_log2(n));
        for (i = 0; i < rows; i++) {
          for (j = (i == 0 ? 1 : 0); j < cols; j++) {
            color_ctx = vp9_get_palette_color_context(buffer, cols, i, j, n,
                                                      color_order);
            for (k = 0; k < n; k++)
              if (buffer[i * cols + j] == color_order[k]) {
                color_new_idx = k;
                break;
              }
            vp9_write_token(w, vp9_palette_color_tree,
                            cm->fc.palette_uv_color_prob[n - 2][color_ctx],
                            &palette_color_encodings[color_new_idx]);
          }
        }
      }
    }
  }
#endif

#if CONFIG_SR_MODE
  write_sr(cm, xd, mbmi->segment_id, mi, w);
  if (mbmi->sr && !mbmi->skip) {
    assert(mbmi->tx_size == max_txsize_lookup[bsize]);
  }
#endif  // CONFIG_SR_MODE
  if (bsize >= BLOCK_8X8 && cm->tx_mode == TX_MODE_SELECT
#if CONFIG_PALETTE
    && !mbmi->palette_enabled[0]
#endif
#if CONFIG_SR_MODE
     && !mi->mbmi.sr
#endif  // CONFIG_SR_MODE
    ) {
    write_selected_tx_size(cm, xd, mbmi->tx_size, bsize, w);
  }

#if CONFIG_NEW_QUANT && QUANT_PROFILES > 1 && !Q_CTX_BASED_PROFILES
  if (cm->base_qindex > Q_THRESHOLD_MIN && cm->base_qindex < Q_THRESHOLD_MAX &&
      mbmi->send_dq_bit &&
      !mbmi->skip &&
      !vp9_segfeature_active(&cm->seg, mbmi->segment_id, SEG_LVL_SKIP)) {
    write_dq_profile(cm, mbmi->dq_off_index, w);
  } else {
    assert(mbmi->dq_off_index == 0);
  }
#endif  // CONFIG_NEW_QUANT  && QUANT_PROFILES > 1 && !Q_CTX_BASED_PROFILES

#if CONFIG_TX_SKIP
  if (bsize >= BLOCK_8X8) {
    int q_idx = vp9_get_qindex(seg, mbmi->segment_id, cm->base_qindex);
    int try_tx_skip = q_idx <= tx_skip_q_thresh_intra;
    if (try_tx_skip) {
      if (xd->lossless) {
        if (mbmi->tx_size == TX_4X4)
          vp9_write(w, mbmi->tx_skip[0], cm->fc.y_tx_skip_prob[0]);
        if (get_uv_tx_size(mbmi, &xd->plane[1]) == TX_4X4)
          vp9_write(w, mbmi->tx_skip[1],
                    cm->fc.uv_tx_skip_prob[mbmi->tx_skip[0]]);
      } else {
        vp9_write(w, mbmi->tx_skip[0], cm->fc.y_tx_skip_prob[0]);
        vp9_write(w, mbmi->tx_skip[1],
                  cm->fc.uv_tx_skip_prob[mbmi->tx_skip[0]]);
      }
    }
  }
#endif

  if (bsize >= BLOCK_8X8) {
#if CONFIG_PALETTE
    if (!mbmi->palette_enabled[0])
#endif  // CONFIG_PALETTE
#if CONFIG_INTRABC
    if (!is_intrabc_mode(mbmi->mode))
#endif  // CONFIG_INTRABC
    write_intra_mode(w, mbmi->mode, get_y_mode_probs(mi, above_mi, left_mi, 0));
#if CONFIG_FILTERINTRA
    if (is_filter_allowed(mbmi->mode) && is_filter_enabled(mbmi->tx_size)
#if CONFIG_PALETTE
            && !mbmi->palette_enabled[0]
#endif  // CONFIG_PALETTE
    )
      vp9_write(w, mbmi->filterbit,
                cm->fc.filterintra_prob[mbmi->tx_size][mbmi->mode]);
#endif  // CONFIG_FILTERINTRA
#if CONFIG_INTRABC
    if (mbmi->mode == NEWDV) {
      int_mv dv_ref = mbmi->ref_mvs[INTRA_FRAME][0];
      vp9_encode_dv(w, &mbmi->mv[0].as_mv, &dv_ref.as_mv, ndvc);
    }
#endif  // CONFIG_INTRABC
  } else {
    const int num_4x4_w = num_4x4_blocks_wide_lookup[bsize];
    const int num_4x4_h = num_4x4_blocks_high_lookup[bsize];
    int idx, idy;

    for (idy = 0; idy < 2; idy += num_4x4_h) {
      for (idx = 0; idx < 2; idx += num_4x4_w) {
        const int block = idy * 2 + idx;
        write_intra_mode(w, mi->bmi[block].as_mode,
                         get_y_mode_probs(mi, above_mi, left_mi, block));
#if CONFIG_FILTERINTRA
        if (is_filter_allowed(mi->bmi[block].as_mode))
          vp9_write(w, mi->b_filter_info[block],
                    cm->fc.filterintra_prob[0][mi->bmi[block].as_mode]);
#endif
      }
    }
  }

#if CONFIG_PALETTE
  if (!mbmi->palette_enabled[1])
#endif  // CONFIG_PALETTE
#if CONFIG_INTRABC
  if (!is_intrabc_mode(mbmi->mode))
#endif  // CONFIG_INTRABC
  write_intra_mode(w, mbmi->uv_mode, vp9_kf_uv_mode_prob[mbmi->mode]);

#if CONFIG_FILTERINTRA
  if (is_filter_allowed(mbmi->uv_mode) &&
      is_filter_enabled(get_uv_tx_size(mbmi, &xd->plane[1]))
#if CONFIG_PALETTE
      && !mbmi->palette_enabled[1]
#endif  // CONFIG_PALETTE
  )
    vp9_write(w, mbmi->uv_filterbit,
              cm->fc.filterintra_prob[get_uv_tx_size(mbmi, &xd->plane[1])][mbmi->uv_mode]);
#endif  // CONFIG_FILTERINTRA
}

static void write_modes_b(VP9_COMP *cpi, const TileInfo *const tile,
                          vp9_writer *w, TOKENEXTRA **tok,
                          const TOKENEXTRA *const tok_end,
#if CONFIG_SUPERTX
                          int supertx_enabled,
#endif
                          int mi_row, int mi_col) {
  const VP9_COMMON *const cm = &cpi->common;
  MACROBLOCKD *const xd = &cpi->mb.e_mbd;
  MODE_INFO *m;

  xd->mi = cm->mi + (mi_row * cm->mi_stride + mi_col);
  m = xd->mi;

  set_mi_row_col(xd, tile,
                 mi_row, num_8x8_blocks_high_lookup[m->mbmi.sb_type],
                 mi_col, num_8x8_blocks_wide_lookup[m->mbmi.sb_type],
                 cm->mi_rows, cm->mi_cols);
  if (frame_is_intra_only(cm)) {
    write_mb_modes_kf(cm, xd, xd->mi, w);
  } else {
    pack_inter_mode_mvs(cpi, m,
#if CONFIG_SUPERTX
                        supertx_enabled,
#endif
                        w);
  }

#if CONFIG_SUPERTX
  if (!supertx_enabled) {
#endif
  assert(*tok < tok_end);
  pack_mb_tokens(w, tok, tok_end, cm->bit_depth);
#if CONFIG_SUPERTX
  }
#endif
}

static void write_partition(const VP9_COMMON *const cm,
                            const MACROBLOCKD *const xd,
                            int hbs, int mi_row, int mi_col,
                            PARTITION_TYPE p, BLOCK_SIZE bsize, vp9_writer *w) {
  const int ctx = partition_plane_context(xd, mi_row, mi_col, bsize);
  const vp9_prob *const probs = get_partition_probs(cm, ctx);
  const int has_rows = (mi_row + hbs) < cm->mi_rows;
  const int has_cols = (mi_col + hbs) < cm->mi_cols;

  if (has_rows && has_cols) {
#if CONFIG_EXT_PARTITION
    if (bsize <= BLOCK_8X8)
      vp9_write_token(w, vp9_partition_tree, probs, &partition_encodings[p]);
    else
      vp9_write_token(w, vp9_ext_partition_tree, probs,
                      &ext_partition_encodings[p]);
#else
    vp9_write_token(w, vp9_partition_tree, probs, &partition_encodings[p]);
#endif
  } else if (!has_rows && has_cols) {
    assert(p == PARTITION_SPLIT || p == PARTITION_HORZ);
    vp9_write(w, p == PARTITION_SPLIT, probs[1]);
  } else if (has_rows && !has_cols) {
    assert(p == PARTITION_SPLIT || p == PARTITION_VERT);
    vp9_write(w, p == PARTITION_SPLIT, probs[2]);
  } else {
    assert(p == PARTITION_SPLIT);
  }
}

static void write_modes_sb(VP9_COMP *cpi,
                           const TileInfo *const tile, vp9_writer *w,
                           TOKENEXTRA **tok, const TOKENEXTRA *const tok_end,
#if CONFIG_SUPERTX
                           int supertx_enabled,
#endif
                           int mi_row, int mi_col, BLOCK_SIZE bsize) {
  const VP9_COMMON *const cm = &cpi->common;
  MACROBLOCKD *const xd = &cpi->mb.e_mbd;

  const int bsl = b_width_log2_lookup[bsize];
  const int bs = (1 << bsl) / 4;
  PARTITION_TYPE partition;
  BLOCK_SIZE subsize;
  MODE_INFO *m = NULL;
#if CONFIG_SUPERTX
  const int pack_token = !supertx_enabled;
#endif

  if (mi_row >= cm->mi_rows || mi_col >= cm->mi_cols)
    return;

  m = cm->mi[mi_row * cm->mi_stride + mi_col].src_mi;
  partition = partition_lookup[bsl][m->mbmi.sb_type];
#if CONFIG_EXT_PARTITION
  partition = get_partition(cm->mi, cm->mi_stride, cm->mi_rows, cm->mi_cols,
                            mi_row, mi_col, bsize);
#endif
  write_partition(cm, xd, bs, mi_row, mi_col, partition, bsize, w);
  subsize = get_subsize(bsize, partition);
#if CONFIG_SUPERTX
  xd->mi = m;
  set_mi_row_col(xd, tile,
                 mi_row, num_8x8_blocks_high_lookup[bsize],
                 mi_col, num_8x8_blocks_wide_lookup[bsize],
                 cm->mi_rows, cm->mi_cols);
  if (!supertx_enabled && cm->frame_type != KEY_FRAME &&
      partition != PARTITION_NONE && bsize <= MAX_SUPERTX_BLOCK_SIZE &&
      !xd->lossless) {
    TX_SIZE supertx_size = bsize_to_tx_size(bsize);
    vp9_prob prob =
        cm->fc.supertx_prob[partition_supertx_context_lookup[partition]]
                           [supertx_size];
    supertx_enabled = (xd->mi[0].mbmi.tx_size == supertx_size);
    vp9_write(w, supertx_enabled, prob);
    if (supertx_enabled) {
      vp9_write(w, xd->mi[0].mbmi.skip, vp9_get_skip_prob(cm, xd));
#if CONFIG_EXT_TX
#if CONFIG_WAVELETS
      if (!xd->mi[0].mbmi.skip)
        vp9_write_token(w, GET_EXT_TX_TREE(supertx_size),
                        cm->fc.ext_tx_prob[supertx_size],
                        &GET_EXT_TX_ENCODINGS(supertx_size)
                            [xd->mi[0].mbmi.ext_txfrm]);
#else
      if (supertx_size <= TX_16X16 && !xd->mi[0].mbmi.skip)
        vp9_write_token(w, vp9_ext_tx_tree, cm->fc.ext_tx_prob[supertx_size],
                        &ext_tx_encodings[xd->mi[0].mbmi.ext_txfrm]);
#endif  // CONFIG_WAVELETS
#endif  // CONFIG_EXT_TX
#if CONFIG_NEW_QUANT && QUANT_PROFILES > 1 && !Q_CTX_BASED_PROFILES
      if (!xd->mi[0].mbmi.skip &&
          cm->base_qindex > Q_THRESHOLD_MIN &&
          cm->base_qindex < Q_THRESHOLD_MAX &&
          xd->mi[0].mbmi.send_dq_bit &&
          !vp9_segfeature_active(
              &cm->seg, xd->mi[0].mbmi.segment_id, SEG_LVL_SKIP)) {
        write_dq_profile(cm, xd->mi[0].mbmi.dq_off_index, w);
      }
#endif  // CONFIG_NEW_QUANT && QUANT_PROFILES > 1 && !Q_CTX_BASED_PROFILES
    }
  }
#endif  // CONFIG_SUPERTX
  if (subsize < BLOCK_8X8) {
    write_modes_b(cpi, tile, w, tok, tok_end,
#if CONFIG_SUPERTX
                  supertx_enabled,
#endif
                  mi_row, mi_col);
  } else {
    switch (partition) {
      case PARTITION_NONE:
        write_modes_b(cpi, tile, w, tok, tok_end,
#if CONFIG_SUPERTX
                      supertx_enabled,
#endif
                      mi_row, mi_col);
        break;
      case PARTITION_HORZ:
        write_modes_b(cpi, tile, w, tok, tok_end,
#if CONFIG_SUPERTX
                      supertx_enabled,
#endif
                      mi_row, mi_col);
        if (mi_row + bs < cm->mi_rows)
          write_modes_b(cpi, tile, w, tok, tok_end,
#if CONFIG_SUPERTX
                        supertx_enabled,
#endif
                        mi_row + bs, mi_col);
        break;
      case PARTITION_VERT:
        write_modes_b(cpi, tile, w, tok, tok_end,
#if CONFIG_SUPERTX
                      supertx_enabled,
#endif
                      mi_row, mi_col);
        if (mi_col + bs < cm->mi_cols)
          write_modes_b(cpi, tile, w, tok, tok_end,
#if CONFIG_SUPERTX
                        supertx_enabled,
#endif
                        mi_row, mi_col + bs);
        break;
      case PARTITION_SPLIT:
        write_modes_sb(cpi, tile, w, tok, tok_end,
#if CONFIG_SUPERTX
                       supertx_enabled,
#endif
                       mi_row, mi_col, subsize);
        write_modes_sb(cpi, tile, w, tok, tok_end,
#if CONFIG_SUPERTX
                       supertx_enabled,
#endif
                       mi_row, mi_col + bs, subsize);
        write_modes_sb(cpi, tile, w, tok, tok_end,
#if CONFIG_SUPERTX
                       supertx_enabled,
#endif
                       mi_row + bs, mi_col, subsize);
        write_modes_sb(cpi, tile, w, tok, tok_end,
#if CONFIG_SUPERTX
                       supertx_enabled,
#endif
                       mi_row + bs, mi_col + bs, subsize);
        break;
#if CONFIG_EXT_PARTITION
      case PARTITION_HORZ_A:
        write_modes_b(cpi, tile, w, tok, tok_end,
#if CONFIG_SUPERTX
                      supertx_enabled,
#endif
                      mi_row, mi_col);
        write_modes_b(cpi, tile, w, tok, tok_end,
#if CONFIG_SUPERTX
                      supertx_enabled,
#endif
                      mi_row, mi_col + bs);
        write_modes_b(cpi, tile, w, tok, tok_end,
#if CONFIG_SUPERTX
                      supertx_enabled,
#endif
                      mi_row + bs, mi_col);
        break;
      case PARTITION_HORZ_B:
        write_modes_b(cpi, tile, w, tok, tok_end,
#if CONFIG_SUPERTX
                      supertx_enabled,
#endif
                      mi_row, mi_col);
        write_modes_b(cpi, tile, w, tok, tok_end,
#if CONFIG_SUPERTX
                      supertx_enabled,
#endif
                      mi_row + bs, mi_col);
        write_modes_b(cpi, tile, w, tok, tok_end,
#if CONFIG_SUPERTX
                      supertx_enabled,
#endif
                      mi_row + bs, mi_col + bs);
        break;
      case PARTITION_VERT_A:
        write_modes_b(cpi, tile, w, tok, tok_end,
#if CONFIG_SUPERTX
                      supertx_enabled,
#endif
                      mi_row, mi_col);
        write_modes_b(cpi, tile, w, tok, tok_end,
#if CONFIG_SUPERTX
                      supertx_enabled,
#endif
                      mi_row + bs, mi_col);
        write_modes_b(cpi, tile, w, tok, tok_end,
#if CONFIG_SUPERTX
                      supertx_enabled,
#endif
                      mi_row, mi_col + bs);
        break;
      case PARTITION_VERT_B:
        write_modes_b(cpi, tile, w, tok, tok_end,
#if CONFIG_SUPERTX
                      supertx_enabled,
#endif
                      mi_row, mi_col);
        write_modes_b(cpi, tile, w, tok, tok_end,
#if CONFIG_SUPERTX
                      supertx_enabled,
#endif
                      mi_row, mi_col + bs);
        write_modes_b(cpi, tile, w, tok, tok_end,
#if CONFIG_SUPERTX
                      supertx_enabled,
#endif
                      mi_row + bs, mi_col + bs);
        break;
#endif
      default:
        assert(0);
    }
  }
#if CONFIG_SUPERTX
  if (partition != PARTITION_NONE && supertx_enabled && pack_token) {
    assert(*tok < tok_end);
    pack_mb_tokens(w, tok, tok_end, cm->bit_depth);
  }
#endif

  // update partition context
#if CONFIG_EXT_PARTITION
  update_ext_partition_context(xd, mi_row, mi_col, subsize, bsize, partition);
#else
  if (bsize >= BLOCK_8X8 &&
      (bsize == BLOCK_8X8 || partition != PARTITION_SPLIT))
    update_partition_context(xd, mi_row, mi_col, subsize, bsize);
#endif
}

static void write_modes(VP9_COMP *cpi,
                        const TileInfo *const tile, vp9_writer *w,
                        TOKENEXTRA **tok, const TOKENEXTRA *const tok_end) {
  int mi_row, mi_col;
#if CONFIG_ROW_TILE
  VP9_COMMON *cm = &cpi->common;
  vpx_memset(&cm->above_seg_context[tile->mi_col_start], 0,
             sizeof(*cm->above_seg_context) *
             mi_cols_aligned_to_sb(tile->mi_col_end - tile->mi_col_start));
#endif
  for (mi_row = tile->mi_row_start; mi_row < tile->mi_row_end;
       mi_row += MI_BLOCK_SIZE) {
    vp9_zero(cpi->mb.e_mbd.left_seg_context);
    for (mi_col = tile->mi_col_start; mi_col < tile->mi_col_end;
         mi_col += MI_BLOCK_SIZE)
      write_modes_sb(cpi, tile, w, tok, tok_end,
#if CONFIG_SUPERTX
                     0,
#endif
                     mi_row, mi_col, BLOCK_LARGEST);
  }
}

static void build_tree_distribution(VP9_COMP *cpi, TX_SIZE tx_size,
                                    vp9_coeff_stats *coef_branch_ct,
                                    vp9_coeff_probs_model *coef_probs) {
  vp9_coeff_count *coef_counts = cpi->coef_counts[tx_size];
  unsigned int (*eob_branch_ct)[REF_TYPES][COEF_BANDS][COEFF_CONTEXTS] =
      cpi->common.counts.eob_branch[tx_size];
  int i, j, k, l, m;

  for (i = 0; i < PLANE_TYPES; ++i) {
    for (j = 0; j < REF_TYPES; ++j) {
      for (k = 0; k < COEF_BANDS; ++k) {
        for (l = 0; l < BAND_COEFF_CONTEXTS(k); ++l) {
          vp9_tree_probs_from_distribution(vp9_coef_tree,
                                           coef_branch_ct[i][j][k][l],
                                           coef_counts[i][j][k][l]);
          coef_branch_ct[i][j][k][l][0][1] = eob_branch_ct[i][j][k][l] -
                                             coef_branch_ct[i][j][k][l][0][0];
          for (m = 0; m < UNCONSTRAINED_NODES; ++m)
            coef_probs[i][j][k][l][m] = get_binary_prob(
                                            coef_branch_ct[i][j][k][l][m][0],
                                            coef_branch_ct[i][j][k][l][m][1]);
        }
      }
    }
  }
}

static void update_coef_probs_common(vp9_writer* const bc, VP9_COMP *cpi,
                                     TX_SIZE tx_size,
                                     vp9_coeff_stats *frame_branch_ct,
                                     vp9_coeff_probs_model *new_coef_probs) {
  vp9_coeff_probs_model *old_coef_probs = cpi->common.fc.coef_probs[tx_size];
  const vp9_prob upd = DIFF_UPDATE_PROB;
  const int entropy_nodes_update = UNCONSTRAINED_NODES;
  int i, j, k, l, t;
  switch (cpi->sf.use_fast_coef_updates) {
    case TWO_LOOP: {
      /* dry run to see if there is any update at all needed */
      int savings = 0;
      int update[2] = {0, 0};
      for (i = 0; i < PLANE_TYPES; ++i) {
        for (j = 0; j < REF_TYPES; ++j) {
          for (k = 0; k < COEF_BANDS; ++k) {
            for (l = 0; l < BAND_COEFF_CONTEXTS(k); ++l) {
              for (t = 0; t < entropy_nodes_update; ++t) {
                vp9_prob newp = new_coef_probs[i][j][k][l][t];
                const vp9_prob oldp = old_coef_probs[i][j][k][l][t];
                int s;
                int u = 0;
                if (t == PIVOT_NODE)
                  s = vp9_prob_diff_update_savings_search_model(
                      frame_branch_ct[i][j][k][l][0],
                      old_coef_probs[i][j][k][l], &newp, upd);
                else
                  s = vp9_prob_diff_update_savings_search(
                      frame_branch_ct[i][j][k][l][t], oldp, &newp, upd);
                if (s > 0 && newp != oldp)
                  u = 1;
                if (u)
                  savings += s - (int)(vp9_cost_zero(upd));
                else
                  savings -= (int)(vp9_cost_zero(upd));
                update[u]++;
              }
            }
          }
        }
      }

      // printf("Update %d %d, savings %d\n", update[0], update[1], savings);
      /* Is coef updated at all */
      if (update[1] == 0 || savings < 0) {
        vp9_write_bit(bc, 0);
        return;
      }
      vp9_write_bit(bc, 1);
      for (i = 0; i < PLANE_TYPES; ++i) {
        for (j = 0; j < REF_TYPES; ++j) {
          for (k = 0; k < COEF_BANDS; ++k) {
            for (l = 0; l < BAND_COEFF_CONTEXTS(k); ++l) {
              // calc probs and branch cts for this frame only
              for (t = 0; t < entropy_nodes_update; ++t) {
                vp9_prob newp = new_coef_probs[i][j][k][l][t];
                vp9_prob *oldp = old_coef_probs[i][j][k][l] + t;
                const vp9_prob upd = DIFF_UPDATE_PROB;
                int s;
                int u = 0;
                if (t == PIVOT_NODE)
                  s = vp9_prob_diff_update_savings_search_model(
                      frame_branch_ct[i][j][k][l][0],
                      old_coef_probs[i][j][k][l], &newp, upd);
                else
                  s = vp9_prob_diff_update_savings_search(
                      frame_branch_ct[i][j][k][l][t],
                      *oldp, &newp, upd);
                if (s > 0 && newp != *oldp)
                  u = 1;
                vp9_write(bc, u, upd);
                if (u) {
                  /* send/use new probability */
                  vp9_write_prob_diff_update(bc, newp, *oldp);
                  *oldp = newp;
                }
              }
            }
          }
        }
      }
      return;
    }

    case ONE_LOOP:
    case ONE_LOOP_REDUCED: {
      const int prev_coef_contexts_to_update =
          cpi->sf.use_fast_coef_updates == ONE_LOOP_REDUCED ?
              COEFF_CONTEXTS >> 1 : COEFF_CONTEXTS;
      const int coef_band_to_update =
          cpi->sf.use_fast_coef_updates == ONE_LOOP_REDUCED ?
              COEF_BANDS >> 1 : COEF_BANDS;
      int updates = 0;
      int noupdates_before_first = 0;
      for (i = 0; i < PLANE_TYPES; ++i) {
        for (j = 0; j < REF_TYPES; ++j) {
          for (k = 0; k < COEF_BANDS; ++k) {
            for (l = 0; l < BAND_COEFF_CONTEXTS(k); ++l) {
              // calc probs and branch cts for this frame only
              for (t = 0; t < entropy_nodes_update; ++t) {
                vp9_prob newp = new_coef_probs[i][j][k][l][t];
                vp9_prob *oldp = old_coef_probs[i][j][k][l] + t;
                int s;
                int u = 0;
                if (l >= prev_coef_contexts_to_update ||
                    k >= coef_band_to_update) {
                  u = 0;
                } else {
                  if (t == PIVOT_NODE)
                    s = vp9_prob_diff_update_savings_search_model(
                        frame_branch_ct[i][j][k][l][0],
                        old_coef_probs[i][j][k][l], &newp, upd);
                  else
                    s = vp9_prob_diff_update_savings_search(
                        frame_branch_ct[i][j][k][l][t],
                        *oldp, &newp, upd);
                  if (s > 0 && newp != *oldp)
                    u = 1;
                }
                updates += u;
                if (u == 0 && updates == 0) {
                  noupdates_before_first++;
                  continue;
                }
                if (u == 1 && updates == 1) {
                  int v;
                  // first update
                  vp9_write_bit(bc, 1);
                  for (v = 0; v < noupdates_before_first; ++v)
                    vp9_write(bc, 0, upd);
                }
                vp9_write(bc, u, upd);
                if (u) {
                  /* send/use new probability */
                  vp9_write_prob_diff_update(bc, newp, *oldp);
                  *oldp = newp;
                }
              }
            }
          }
        }
      }
      if (updates == 0) {
        vp9_write_bit(bc, 0);  // no updates
      }
      return;
    }

    default:
      assert(0);
  }
}

#if CONFIG_TX_SKIP
static void build_tree_distribution_pxd(VP9_COMP *cpi, TX_SIZE tx_size,
                                        vp9_coeff_stats_pxd *coef_branch_ct,
                                        vp9_coeff_probs_pxd *coef_probs) {
  vp9_coeff_counts_pxd *coef_counts = cpi->common.counts.coef_pxd[tx_size];
  unsigned int (*eob_branch_ct)[REF_TYPES][COEFF_CONTEXTS] =
      cpi->common.counts.eob_branch_pxd[tx_size];
  int i, j, l, m;

  for (i = 0; i < PLANE_TYPES; ++i) {
    for (j = 0; j < REF_TYPES; ++j) {
      for (l = 0; l < COEFF_CONTEXTS; ++l) {
        vp9_tree_probs_from_distribution(vp9_coef_tree,
                                         coef_branch_ct[i][j][l],
                                         coef_counts[i][j][l]);
        coef_branch_ct[i][j][l][0][1] = eob_branch_ct[i][j][l] -
            coef_branch_ct[i][j][l][0][0];
        for (m = 0; m < ENTROPY_NODES; ++m)
          coef_probs[i][j][l][m] = get_binary_prob(
              coef_branch_ct[i][j][l][m][0],
              coef_branch_ct[i][j][l][m][1]);
      }
    }
  }
}

static void update_coef_probs_common_pxd(vp9_writer* const bc, VP9_COMP *cpi,
                                         TX_SIZE tx_size,
                                         vp9_coeff_stats_pxd *frame_branch_ct,
                                         vp9_coeff_probs_pxd *new_coef_probs) {
  vp9_coeff_probs_pxd *old_coef_probs = cpi->common.fc.coef_probs_pxd[tx_size];
  const vp9_prob upd = DIFF_UPDATE_PROB;
  const int entropy_nodes_update = ENTROPY_NODES;
  int i, j, l, t;

  switch (cpi->sf.use_fast_coef_updates) {
    case TWO_LOOP: {
      /* dry run to see if there is any update at all needed */
      int savings = 0;
      int update[2] = {0, 0};
      for (i = 0; i < PLANE_TYPES; ++i) {
        for (j = 0; j < REF_TYPES; ++j) {
          for (l = 0; l < COEFF_CONTEXTS; ++l) {
            for (t = 0; t < entropy_nodes_update; ++t) {
              vp9_prob newp = new_coef_probs[i][j][l][t];
              const vp9_prob oldp = old_coef_probs[i][j][l][t];
              int s;
              int u = 0;
              s = vp9_prob_diff_update_savings_search(
                  frame_branch_ct[i][j][l][t], oldp, &newp, upd);
              if (s > 0 && newp != oldp)
                u = 1;
              if (u)
                savings += s - (int)(vp9_cost_zero(upd));
              else
                savings -= (int)(vp9_cost_zero(upd));
              update[u]++;
            }
          }
        }
      }

      if (update[1] == 0 || savings < 0) {
        vp9_write_bit(bc, 0);
        return;
      }
      vp9_write_bit(bc, 1);
      for (i = 0; i < PLANE_TYPES; ++i) {
        for (j = 0; j < REF_TYPES; ++j) {
          for (l = 0; l < COEFF_CONTEXTS; ++l) {
            // calc probs and branch cts for this frame only
            for (t = 0; t < entropy_nodes_update; ++t) {
              vp9_prob newp = new_coef_probs[i][j][l][t];
              vp9_prob *oldp = old_coef_probs[i][j][l] + t;
              const vp9_prob upd = DIFF_UPDATE_PROB;
              int s;
              int u = 0;
              s = vp9_prob_diff_update_savings_search(
                  frame_branch_ct[i][j][l][t],
                  *oldp, &newp, upd);
              if (s > 0 && newp != *oldp)
                u = 1;
              vp9_write(bc, u, upd);
              if (u) {
                /* send/use new probability */
                vp9_write_prob_diff_update(bc, newp, *oldp);
                *oldp = newp;
              }
            }
          }
        }
      }
      return;
    }

    case ONE_LOOP:
    case ONE_LOOP_REDUCED: {
      int updates = 0;
      int noupdates_before_first = 0;
      for (i = 0; i < PLANE_TYPES; ++i) {
        for (j = 0; j < REF_TYPES; ++j) {
          for (l = 0; l < COEFF_CONTEXTS; ++l) {
            // calc probs and branch cts for this frame only
            for (t = 0; t < entropy_nodes_update; ++t) {
              vp9_prob newp = new_coef_probs[i][j][l][t];
              vp9_prob *oldp = old_coef_probs[i][j][l] + t;
              int s;
              int u = 0;

              s = vp9_prob_diff_update_savings_search(
                  frame_branch_ct[i][j][l][t],
                  *oldp, &newp, upd);
              if (s > 0 && newp != *oldp)
                u = 1;

              updates += u;
              if (u == 0 && updates == 0) {
                noupdates_before_first++;
                continue;
              }
              if (u == 1 && updates == 1) {
                int v;
                // first update
                vp9_write_bit(bc, 1);
                for (v = 0; v < noupdates_before_first; ++v)
                  vp9_write(bc, 0, upd);
              }
              vp9_write(bc, u, upd);
              if (u) {
                // send/use new probability
                vp9_write_prob_diff_update(bc, newp, *oldp);
                *oldp = newp;
              }
            }
          }
        }
      }
      if (updates == 0) {
        vp9_write_bit(bc, 0);  // no updates
      }
      return;
    }

    default:
      assert(0);
  }
}
#endif  // CONFIG_TX_SKIP

static void update_coef_probs(VP9_COMP *cpi, vp9_writer* w) {
  const TX_MODE tx_mode = cpi->common.tx_mode;
  const TX_SIZE max_tx_size = tx_mode_to_biggest_tx_size[tx_mode];
  TX_SIZE tx_size;
  vp9_coeff_stats frame_branch_ct[TX_SIZES][PLANE_TYPES];
  vp9_coeff_probs_model frame_coef_probs[TX_SIZES][PLANE_TYPES];
#if CONFIG_TX_SKIP
  vp9_coeff_stats_pxd frame_branch_ct_pxd[TX_SIZES][PLANE_TYPES];
  vp9_coeff_probs_pxd frame_coef_probs_pxd[TX_SIZES][PLANE_TYPES];
#endif  // CONFIG_TX_SKIP

  for (tx_size = TX_4X4; tx_size <= max_tx_size; ++tx_size)
    build_tree_distribution(cpi, tx_size, frame_branch_ct[tx_size],
                            frame_coef_probs[tx_size]);

  for (tx_size = TX_4X4; tx_size <= max_tx_size; ++tx_size)
    update_coef_probs_common(w, cpi, tx_size, frame_branch_ct[tx_size],
                             frame_coef_probs[tx_size]);

#if CONFIG_TX_SKIP
  if (FOR_SCREEN_CONTENT) {
    for (tx_size = TX_4X4; tx_size <= max_tx_size; ++tx_size)
      build_tree_distribution_pxd(cpi, tx_size, frame_branch_ct_pxd[tx_size],
                                  frame_coef_probs_pxd[tx_size]);

    for (tx_size = TX_4X4; tx_size <= max_tx_size; ++tx_size)
      update_coef_probs_common_pxd(w, cpi, tx_size,
                                   frame_branch_ct_pxd[tx_size],
                                   frame_coef_probs_pxd[tx_size]);
  }
#endif  // CONFIG_TX_SKIP
}

static void encode_loopfilter(VP9_COMMON *cm,
                              struct vp9_write_bit_buffer *wb) {
  int i;
  struct loopfilter *lf = &cm->lf;

  // Encode the loop filter level and type
  vp9_wb_write_literal(wb, lf->filter_level, 6);
  vp9_wb_write_literal(wb, lf->sharpness_level, 3);

  // Write out loop filter deltas applied at the MB level based on mode or
  // ref frame (if they are enabled).
  vp9_wb_write_bit(wb, lf->mode_ref_delta_enabled);

  if (lf->mode_ref_delta_enabled) {
    vp9_wb_write_bit(wb, lf->mode_ref_delta_update);
    if (lf->mode_ref_delta_update) {
      for (i = 0; i < MAX_REF_LF_DELTAS; i++) {
        const int delta = lf->ref_deltas[i];
        const int changed = delta != lf->last_ref_deltas[i];
        vp9_wb_write_bit(wb, changed);
        if (changed) {
          lf->last_ref_deltas[i] = delta;
          vp9_wb_write_literal(wb, abs(delta) & 0x3F, 6);
          vp9_wb_write_bit(wb, delta < 0);
        }
      }

      for (i = 0; i < MAX_MODE_LF_DELTAS; i++) {
        const int delta = lf->mode_deltas[i];
        const int changed = delta != lf->last_mode_deltas[i];
        vp9_wb_write_bit(wb, changed);
        if (changed) {
          lf->last_mode_deltas[i] = delta;
          vp9_wb_write_literal(wb, abs(delta) & 0x3F, 6);
          vp9_wb_write_bit(wb, delta < 0);
        }
      }
    }
  }
#if CONFIG_LOOP_POSTFILTER
  vp9_wb_write_bit(wb, lf->bilateral_level != lf->last_bilateral_level);
  if (lf->bilateral_level != lf->last_bilateral_level) {
    int level = lf->bilateral_level -
                (lf->bilateral_level > lf->last_bilateral_level);
    vp9_wb_write_literal(wb, level,
                         vp9_bilateral_level_bits(cm));
  }
#endif  // CONFIG_LOOP_POSTFILTER
}

static void write_delta_q(struct vp9_write_bit_buffer *wb, int delta_q) {
  if (delta_q != 0) {
    vp9_wb_write_bit(wb, 1);
    vp9_wb_write_literal(wb, abs(delta_q), 4);
    vp9_wb_write_bit(wb, delta_q < 0);
  } else {
    vp9_wb_write_bit(wb, 0);
  }
}

static void encode_quantization(const VP9_COMMON *const cm,
                                struct vp9_write_bit_buffer *wb) {
  vp9_wb_write_literal(wb, cm->base_qindex, QINDEX_BITS);
  write_delta_q(wb, cm->y_dc_delta_q);
  write_delta_q(wb, cm->uv_dc_delta_q);
  write_delta_q(wb, cm->uv_ac_delta_q);
}

static void encode_segmentation(VP9_COMMON *cm, MACROBLOCKD *xd,
                                struct vp9_write_bit_buffer *wb) {
  int i, j;

  const struct segmentation *seg = &cm->seg;

  vp9_wb_write_bit(wb, seg->enabled);
  if (!seg->enabled)
    return;

  // Segmentation map
  vp9_wb_write_bit(wb, seg->update_map);
  if (seg->update_map) {
    // Select the coding strategy (temporal or spatial)
    vp9_choose_segmap_coding_method(cm, xd);
    // Write out probabilities used to decode unpredicted  macro-block segments
    for (i = 0; i < SEG_TREE_PROBS; i++) {
      const int prob = seg->tree_probs[i];
      const int update = prob != MAX_PROB;
      vp9_wb_write_bit(wb, update);
      if (update)
        vp9_wb_write_literal(wb, prob, 8);
    }

    // Write out the chosen coding method.
    vp9_wb_write_bit(wb, seg->temporal_update);
    if (seg->temporal_update) {
      for (i = 0; i < PREDICTION_PROBS; i++) {
        const int prob = seg->pred_probs[i];
        const int update = prob != MAX_PROB;
        vp9_wb_write_bit(wb, update);
        if (update)
          vp9_wb_write_literal(wb, prob, 8);
      }
    }
  }

  // Segmentation data
  vp9_wb_write_bit(wb, seg->update_data);
  if (seg->update_data) {
    vp9_wb_write_bit(wb, seg->abs_delta);

    for (i = 0; i < MAX_SEGMENTS; i++) {
      for (j = 0; j < SEG_LVL_MAX; j++) {
        const int active = vp9_segfeature_active(seg, i, j);
        vp9_wb_write_bit(wb, active);
        if (active) {
          const int data = vp9_get_segdata(seg, i, j);
          const int data_max = vp9_seg_feature_data_max(j);

          if (vp9_is_segfeature_signed(j)) {
            encode_unsigned_max(wb, abs(data), data_max);
            vp9_wb_write_bit(wb, data < 0);
          } else {
            encode_unsigned_max(wb, data, data_max);
          }
        }
      }
    }
  }
}

static void encode_txfm_probs(VP9_COMMON *cm, vp9_writer *w) {
  // Mode
#if CONFIG_TX64X64
  if (cm->tx_mode == ALLOW_16X16 || cm->tx_mode == ALLOW_32X32) {
    vp9_write_literal(w, 2, 2);
    vp9_write_bit(w, cm->tx_mode == ALLOW_32X32);
  } else if (cm->tx_mode == ALLOW_64X64 || cm->tx_mode == TX_MODE_SELECT) {
    vp9_write_literal(w, 3, 2);
    vp9_write_bit(w, cm->tx_mode == TX_MODE_SELECT);
  } else {
    vp9_write_literal(w, cm->tx_mode, 2);
  }
#else
  vp9_write_literal(w, MIN(cm->tx_mode, ALLOW_32X32), 2);
  if (cm->tx_mode >= ALLOW_32X32)
    vp9_write_bit(w, cm->tx_mode == TX_MODE_SELECT);
#endif  // CONFIG_TX64X64

  // Probabilities
  if (cm->tx_mode == TX_MODE_SELECT) {
    int i, j;
    unsigned int ct_8x8p[1][2];
    unsigned int ct_16x16p[2][2];
    unsigned int ct_32x32p[3][2];
#if CONFIG_TX64X64
    unsigned int ct_64x64p[4][2];
#endif

    for (i = 0; i < TX_SIZE_CONTEXTS; i++) {
      tx_counts_to_branch_counts_8x8(cm->counts.tx.p8x8[i], ct_8x8p);
      for (j = 0; j < 1; j++)
        vp9_cond_prob_diff_update(w, &cm->fc.tx_probs.p8x8[i][j], ct_8x8p[j]);
    }

    for (i = 0; i < TX_SIZE_CONTEXTS; i++) {
      tx_counts_to_branch_counts_16x16(cm->counts.tx.p16x16[i], ct_16x16p);
      for (j = 0; j < 2; j++)
        vp9_cond_prob_diff_update(w, &cm->fc.tx_probs.p16x16[i][j],
                                  ct_16x16p[j]);
    }

    for (i = 0; i < TX_SIZE_CONTEXTS; i++) {
      tx_counts_to_branch_counts_32x32(cm->counts.tx.p32x32[i], ct_32x32p);
      for (j = 0; j < 3; j++)
        vp9_cond_prob_diff_update(w, &cm->fc.tx_probs.p32x32[i][j],
                                  ct_32x32p[j]);
    }

#if CONFIG_TX64X64
    for (i = 0; i < TX_SIZE_CONTEXTS; i++) {
      tx_counts_to_branch_counts_64x64(cm->counts.tx.p64x64[i], ct_64x64p);
      for (j = 0; j < 4; j++)
        vp9_cond_prob_diff_update(w, &cm->fc.tx_probs.p64x64[i][j],
                                  ct_64x64p[j]);
    }
#endif  // CONFIG_TX64X64
  }
}

static void write_interp_filter(INTERP_FILTER filter,
                                struct vp9_write_bit_buffer *wb) {
#if CONFIG_BITSTREAM_FIXES
  vp9_wb_write_bit(wb, filter == SWITCHABLE);
  if (filter != SWITCHABLE)
    vp9_wb_write_literal(wb, filter, 2);
#else
  const int filter_to_literal[] = { 1, 0, 2, 3 };

  vp9_wb_write_bit(wb, filter == SWITCHABLE);
  if (filter != SWITCHABLE)
    vp9_wb_write_literal(wb, filter_to_literal[filter], 2);
#endif
}

static void fix_interp_filter(VP9_COMMON *cm) {
  if (cm->interp_filter == SWITCHABLE) {
    // Check to see if only one of the filters is actually used
    int count[SWITCHABLE_FILTERS];
    int i, j, c = 0;
    for (i = 0; i < SWITCHABLE_FILTERS; ++i) {
      count[i] = 0;
      for (j = 0; j < SWITCHABLE_FILTER_CONTEXTS; ++j)
        count[i] += cm->counts.switchable_interp[j][i];
      c += (count[i] > 0);
    }
    if (c == 1) {
      // Only one filter is used. So set the filter at frame level
      for (i = 0; i < SWITCHABLE_FILTERS; ++i) {
        if (count[i]) {
          cm->interp_filter = i;
          break;
        }
      }
    }
  }
}

#if CONFIG_ROW_TILE
// Decide how many bytes used for storing tile column size and tile size based
// on the tile dimensions or the information gathered in a previous bitstream
// packing while recode is on.
static INLINE void set_tile_size_in_bytes(VP9_COMP *cpi,
                                          int tile_width, int tile_height,
                                          int final_packing) {
  VP9_COMMON *const cm = &cpi->common;

  // Default number of bytes used is 4.
  cm->tile_size_bytes = 4;
  cm->tile_col_size_bytes = 4;

  if (!final_packing ||
      (cpi->max_tile_size == UINT_MAX || cpi->max_tile_col_size == UINT_MAX)) {
    // Determine num of bytes by tile size.
    if (tile_width * tile_height <= 4)
      cm->tile_size_bytes = 2;

    if (!final_packing) {
      cpi->max_tile_size = 0;
      cpi->max_tile_col_size = 0;
    }
  } else {
    // When we have ideas on maximum tile size, use that to estimate the num of
    // bytes needed. Note: This decision is not always guaranteed to be true.
    // In later bitstream packing, need to check if the actual tile size is
    // out of the range.
    if (cpi->max_tile_size < ONE_BYTE_THRESH(0))
      cm->tile_size_bytes = 1;
    else if (cpi->max_tile_size < TWO_BYTE_THRESH(0))
      cm->tile_size_bytes = 2;
    else if (cpi->max_tile_size < THREE_BYTE_THRESH(0))
      cm->tile_size_bytes = 3;

    if (cpi->max_tile_col_size < ONE_BYTE_THRESH(1))
      cm->tile_col_size_bytes = 1;
    else if (cpi->max_tile_col_size < TWO_BYTE_THRESH(1))
      cm->tile_col_size_bytes = 2;
    else if (cpi->max_tile_col_size < THREE_BYTE_THRESH(1))
      cm->tile_col_size_bytes = 3;
  }
}

static void write_tile_info(VP9_COMP *cpi,
                            struct vp9_write_bit_buffer *wb,
                            int final_packing) {
  VP9_COMMON *const cm = &cpi->common;

  int tile_width  = mi_cols_aligned_to_sb(cm->tile_width) >> MI_BLOCK_SIZE_LOG2;
  int tile_height =
      mi_cols_aligned_to_sb(cm->tile_height) >> MI_BLOCK_SIZE_LOG2;
  vp9_wb_write_literal(wb, tile_width, 6);
  vp9_wb_write_literal(wb, tile_height, 6);

  set_tile_size_in_bytes(cpi, tile_width, tile_height, final_packing);

  assert(cm->tile_col_size_bytes > 0 && cm->tile_col_size_bytes <= 4);
  assert(cm->tile_size_bytes > 0 && cm->tile_size_bytes <= 4);

  // Write the num of bytes decision to bitstream
  vp9_wb_write_literal(wb, cm->tile_col_size_bytes - 1, 2);
  vp9_wb_write_literal(wb, cm->tile_size_bytes - 1, 2);
}
#else
static void write_tile_info(VP9_COMMON *const cm,
                            struct vp9_write_bit_buffer *wb) {
  int min_log2_tiles, max_log2_tiles, ones;
  vp9_get_tile_n_bits(cm->mi_cols, &min_log2_tiles, &max_log2_tiles);

  // columns
  ones = cm->log2_tile_cols - min_log2_tiles;
  while (ones--)
    vp9_wb_write_bit(wb, 1);

  if (cm->log2_tile_cols < max_log2_tiles)
    vp9_wb_write_bit(wb, 0);

  // rows
  vp9_wb_write_bit(wb, cm->log2_tile_rows != 0);
  if (cm->log2_tile_rows != 0)
    vp9_wb_write_bit(wb, cm->log2_tile_rows != 1);
}
#endif

static int get_refresh_mask(VP9_COMP *cpi) {
  if (vp9_preserve_existing_gf(cpi)) {
    // We have decided to preserve the previously existing golden frame as our
    // new ARF frame. However, in the short term we leave it in the GF slot and,
    // if we're updating the GF with the current decoded frame, we save it
    // instead to the ARF slot.
    // Later, in the function vp9_encoder.c:vp9_update_reference_frames() we
    // will swap gld_fb_idx and alt_fb_idx to achieve our objective. We do it
    // there so that it can be done outside of the recode loop.
    // Note: This is highly specific to the use of ARF as a forward reference,
    // and this needs to be generalized as other uses are implemented
    // (like RTC/temporal scalability).
    return (cpi->refresh_last_frame << cpi->lst_fb_idx) |
#if CONFIG_MULTI_REF
           (cpi->refresh_last2_frame << cpi->lst2_fb_idx) |
           (cpi->refresh_last3_frame << cpi->lst3_fb_idx) |
           (cpi->refresh_last4_frame << cpi->lst4_fb_idx) |
#endif  // CONFIG_MULTI_REF
           (cpi->refresh_golden_frame << cpi->alt_fb_idx);
  } else {
    int arf_idx = cpi->alt_fb_idx;
    if ((cpi->oxcf.pass == 2) && cpi->multi_arf_allowed) {
      const GF_GROUP *const gf_group = &cpi->twopass.gf_group;
      arf_idx = gf_group->arf_update_idx[gf_group->index];
    }
    return (cpi->refresh_last_frame << cpi->lst_fb_idx) |
#if CONFIG_MULTI_REF
           (cpi->refresh_last2_frame << cpi->lst2_fb_idx) |
           (cpi->refresh_last3_frame << cpi->lst3_fb_idx) |
           (cpi->refresh_last4_frame << cpi->lst4_fb_idx) |
#endif  // CONFIG_MULTI_REF
           (cpi->refresh_golden_frame << cpi->gld_fb_idx) |
           (cpi->refresh_alt_ref_frame << arf_idx);
  }
}

#if CONFIG_ROW_TILE
static INLINE void setup_size_storing(int num_bytes, MemPut *output,
                                      unsigned int *size_limit, int type) {
  *output = mem_put_be32;
  *size_limit = UINT_MAX;
  if (num_bytes == 3) {
    *output = mem_put_be24;
    *size_limit = THREE_BYTE_LIMIT(type);
  } else if (num_bytes == 2) {
    *output = mem_put_be16;
    *size_limit = TWO_BYTE_LIMIT(type);
  } else if (num_bytes == 1) {
    *output = mem_put_be8;
    *size_limit = ONE_BYTE_LIMIT(type);
  }
}
#endif

#if CONFIG_ROW_TILE
static size_t encode_tiles(VP9_COMP *cpi, uint8_t *data_ptr,
                           int final_packing) {
#else
static size_t encode_tiles(VP9_COMP *cpi, uint8_t *data_ptr) {
#endif

  VP9_COMMON *const cm = &cpi->common;
  vp9_writer residual_bc;

  int tile_row, tile_col;
#if CONFIG_ROW_TILE
  TOKENEXTRA *(*tok)[1024] = cpi->tile_tok;
  TileInfo (*tile)[1024] = cpi->tile_info;
  EncTileBuffer (*tile_buf)[1024] = cpi->tile_buffers;
#else
  TOKENEXTRA *tok[4][1 << 6];
  TileInfo tile[4][1 << 6];
#endif
  TOKENEXTRA *tok_end;
  size_t total_size = 0;
  const int tile_cols = cm->tile_cols;
  const int tile_rows = cm->tile_rows;
  TOKENEXTRA *pre_tok = cpi->tok;
  int tile_tok = 0;

#if CONFIG_ROW_TILE
  MemPut output_tile_size;
  MemPut output_tile_col_size;
  unsigned int tile_size_limit;
  unsigned int tile_col_size_limit;

  setup_size_storing(cm->tile_size_bytes, &output_tile_size, &tile_size_limit,
                     0);
  setup_size_storing(cm->tile_col_size_bytes, &output_tile_col_size,
                     &tile_col_size_limit, 1);
#endif

  vpx_memset(cm->above_seg_context, 0, sizeof(*cm->above_seg_context) *
             mi_cols_aligned_to_sb(cm->mi_cols));

  for (tile_row = 0; tile_row < tile_rows; ++tile_row) {
    for (tile_col = 0; tile_col < tile_cols; ++tile_col) {
      vp9_tile_init(&tile[tile_row][tile_col], cm, tile_row, tile_col);

      tok[tile_row][tile_col] = pre_tok + tile_tok;
      pre_tok = tok[tile_row][tile_col];
      tile_tok = allocated_tokens(tile[tile_row][tile_col]);
    }
  }

#if CONFIG_ROW_TILE
  for (tile_col = 0; tile_col < tile_cols; tile_col++) {
    int is_last_col = (tile_col == tile_cols - 1);
    size_t col_offset = total_size;

    if (!is_last_col)
      total_size += cm->tile_col_size_bytes;

    for (tile_row = 0; tile_row < tile_rows; tile_row++) {
      TileInfo * const ptile = &tile[tile_row][tile_col];
      EncTileBuffer *const ptile_buf = &tile_buf[tile_row][tile_col];
      uint8_t *source;
      int write_tile_data = 1;

      tok_end = tok[tile_row][tile_col] + cpi->tok_count[tile_row][tile_col];

      // Is CONFIG_ROW_TILE = 1, every tile in the row has a header even for
      // the last one.
      ptile_buf->data_start = data_ptr + total_size;
      source = data_ptr + total_size + cm->tile_size_bytes;

      vp9_start_encode(&residual_bc, source);

      write_modes(cpi, ptile, &residual_bc, &tok[tile_row][tile_col], tok_end);
      assert(tok[tile_row][tile_col] == tok_end);
      vp9_stop_encode(&residual_bc);

      // If it is not final packing, record the maximum tile size we see,
      // otherwise, check if the tile size is out of the range.
      if (!final_packing) {
        if (cpi->max_tile_size < residual_bc.pos)
          cpi->max_tile_size = residual_bc.pos;
      } else {
        if (residual_bc.pos > tile_size_limit)
          return 0;
      }

      ptile_buf->data_size = residual_bc.pos;

      // Check if this tile is a copy tile.
      // Very low chances to have copy tiles on the key frame. Thus, don't
      // search on key frame to reduce unnecessary search.
      if (cm->frame_type != KEY_FRAME && final_packing) {
        const MV32 candidates[1] = {{1, 0}};
        int i;

        assert(cm->tile_size_bytes >= 1);

        // (TODO: yunqingwang) For now, only above tile is checked and used.
        // More candidates such as left tile can be added later.
        for (i = 0; i < 1; i++) {
          int cand_row = tile_row - candidates[0].row;
          int cand_col = tile_col - candidates[0].col;
          uint8_t tile_hdr;
          uint8_t *ref_tile;
          unsigned int ref_tile_size;
          int identical_tile_offset = 0;

          if (tile_row == 0 )
            continue;

          tile_hdr = *(tile_buf[cand_row][cand_col].data_start);

          // Read out tcm bit
          if ((tile_hdr >> 7) == 1) {
            // The candidate is a copy tile itself
            tile_hdr &= 0x7f;
            identical_tile_offset = tile_hdr + 1;
            ref_tile = tile_buf[cand_row - tile_hdr][cand_col].data_start +
                cm->tile_size_bytes;
            ref_tile_size = tile_buf[cand_row - tile_hdr][cand_col].data_size;
          } else {
            identical_tile_offset = 1;
            ref_tile = tile_buf[cand_row][cand_col].data_start +
                cm->tile_size_bytes;
            ref_tile_size = tile_buf[cand_row][cand_col].data_size;
          }

          if (identical_tile_offset < 128 && ref_tile_size == residual_bc.pos) {
            unsigned int m;
            uint8_t *cur_tile = tile_buf[tile_row][tile_col].data_start +
                cm->tile_size_bytes;
            int match = 1;

            for (m = 0; m < residual_bc.pos; m++) {
              if (*ref_tile++ != *cur_tile++) {
                match = 0;
                break;
              }
            }

            if (match) {
              write_tile_data = 0;
              identical_tile_offset |= 0x80;
              identical_tile_offset <<= (cm->tile_size_bytes - 1) * 8;
              output_tile_size(data_ptr + total_size, identical_tile_offset);
              break;
            }
          }
        }
      }

      if (write_tile_data) {
        // size of this tile
        output_tile_size(data_ptr + total_size, residual_bc.pos);
        total_size += residual_bc.pos;
      }

      total_size += cm->tile_size_bytes;
    }

    if (!is_last_col) {
      size_t col_size = total_size - col_offset - cm->tile_col_size_bytes;
      output_tile_col_size(data_ptr + col_offset, col_size);

      // If it is not final packing, record the maximum tile column size we see,
      // otherwise, check if the tile size is out of the range.
      if (!final_packing) {
        if (cpi->max_tile_col_size < col_size)
          cpi->max_tile_col_size = col_size;
      } else {
        if (col_size > tile_col_size_limit)
          return 0;
      }
    }
  }
#else
  for (tile_row = 0; tile_row < tile_rows; tile_row++) {
    for (tile_col = 0; tile_col < tile_cols; tile_col++) {
      const TileInfo * const ptile = &tile[tile_row][tile_col];

      tok_end = tok[tile_row][tile_col] + cpi->tok_count[tile_row][tile_col];

      if (tile_col < tile_cols - 1 || tile_row < tile_rows - 1)
        vp9_start_encode(&residual_bc, data_ptr + total_size + 4);
      else
        vp9_start_encode(&residual_bc, data_ptr + total_size);

      write_modes(cpi, ptile, &residual_bc, &tok[tile_row][tile_col], tok_end);
      assert(tok[tile_row][tile_col] == tok_end);
      vp9_stop_encode(&residual_bc);
      if (tile_col < tile_cols - 1 || tile_row < tile_rows - 1) {
        // size of this tile
        mem_put_be32(data_ptr + total_size, residual_bc.pos);
        total_size += 4;
      }

      total_size += residual_bc.pos;
    }
  }
#endif

  return total_size;
}

static void write_display_size(const VP9_COMMON *cm,
                               struct vp9_write_bit_buffer *wb) {
  const int scaling_active = cm->width != cm->display_width ||
                             cm->height != cm->display_height;
  vp9_wb_write_bit(wb, scaling_active);
  if (scaling_active) {
    vp9_wb_write_literal(wb, cm->display_width - 1, 16);
    vp9_wb_write_literal(wb, cm->display_height - 1, 16);
  }
}

static void write_frame_size(const VP9_COMMON *cm,
                             struct vp9_write_bit_buffer *wb) {
  vp9_wb_write_literal(wb, cm->width - 1, 16);
  vp9_wb_write_literal(wb, cm->height - 1, 16);

  write_display_size(cm, wb);
}

static void write_frame_size_with_refs(VP9_COMP *cpi,
                                       struct vp9_write_bit_buffer *wb) {
  VP9_COMMON *const cm = &cpi->common;
  int found = 0;

  MV_REFERENCE_FRAME ref_frame;
  for (ref_frame = LAST_FRAME; ref_frame <= ALTREF_FRAME; ++ref_frame) {
    YV12_BUFFER_CONFIG *cfg = get_ref_frame_buffer(cpi, ref_frame);
    found = cm->width == cfg->y_crop_width &&
            cm->height == cfg->y_crop_height;

    vp9_wb_write_bit(wb, found);
    if (found) {
      break;
    }
  }

  if (!found) {
    vp9_wb_write_literal(wb, cm->width - 1, 16);
    vp9_wb_write_literal(wb, cm->height - 1, 16);
  }

  write_display_size(cm, wb);
}

static void write_sync_code(struct vp9_write_bit_buffer *wb) {
  vp9_wb_write_literal(wb, VP9_SYNC_CODE_0, 8);
  vp9_wb_write_literal(wb, VP9_SYNC_CODE_1, 8);
  vp9_wb_write_literal(wb, VP9_SYNC_CODE_2, 8);
}

static void write_profile(BITSTREAM_PROFILE profile,
                          struct vp9_write_bit_buffer *wb) {
  switch (profile) {
    case PROFILE_0:
      vp9_wb_write_literal(wb, 0, 2);
      break;
    case PROFILE_1:
      vp9_wb_write_literal(wb, 2, 2);
      break;
    case PROFILE_2:
      vp9_wb_write_literal(wb, 1, 2);
      break;
    case PROFILE_3:
      vp9_wb_write_literal(wb, 6, 3);
      break;
    default:
      assert(0);
  }
}

static void write_bitdepth_colorspace_sampling(
    VP9_COMMON *const cm, struct vp9_write_bit_buffer *wb) {
  if (cm->profile >= PROFILE_2) {
    assert(cm->bit_depth > VPX_BITS_8);
    vp9_wb_write_bit(wb, cm->bit_depth == VPX_BITS_10 ? 0 : 1);
  }
  vp9_wb_write_literal(wb, cm->color_space, 3);
  if (cm->color_space != VPX_CS_SRGB) {
    vp9_wb_write_bit(wb, 0);  // 0: [16, 235] (i.e. xvYCC), 1: [0, 255]
    if (cm->profile == PROFILE_1 || cm->profile == PROFILE_3) {
      assert(cm->subsampling_x != 1 || cm->subsampling_y != 1);
      vp9_wb_write_bit(wb, cm->subsampling_x);
      vp9_wb_write_bit(wb, cm->subsampling_y);
      vp9_wb_write_bit(wb, 0);  // unused
    } else {
      assert(cm->subsampling_x == 1 && cm->subsampling_y == 1);
    }
  } else {
    assert(cm->profile == PROFILE_1 || cm->profile == PROFILE_3);
    vp9_wb_write_bit(wb, 0);  // unused
  }
}

#if CONFIG_ROW_TILE
static void write_uncompressed_header(VP9_COMP *cpi,
                                      struct vp9_write_bit_buffer *wb,
                                      int final_packing) {
#else
static void write_uncompressed_header(VP9_COMP *cpi,
                                      struct vp9_write_bit_buffer *wb) {
#endif
  VP9_COMMON *const cm = &cpi->common;

  vp9_wb_write_literal(wb, VP9_FRAME_MARKER, 2);

  write_profile(cm->profile, wb);

  vp9_wb_write_bit(wb, 0);  // show_existing_frame
  vp9_wb_write_bit(wb, cm->frame_type);
  vp9_wb_write_bit(wb, cm->show_frame);
  vp9_wb_write_bit(wb, cm->error_resilient_mode);

#if CONFIG_MULTI_REF
  cpi->refresh_last2_frame =
      (cm->frame_type == KEY_FRAME || cpi->refresh_last_frame) ? 1 : 0;
  cpi->refresh_last3_frame = cpi->refresh_last2_frame ? 1 : 0;
  cpi->refresh_last4_frame = cpi->refresh_last3_frame ? 1 : 0;
#endif  // CONFIG_MULTI_REF

  if (cm->frame_type == KEY_FRAME) {
    write_sync_code(wb);
    write_bitdepth_colorspace_sampling(cm, wb);
    write_frame_size(cm, wb);
  } else {
    if (!cm->show_frame)
      vp9_wb_write_bit(wb, cm->intra_only);

    if (!cm->error_resilient_mode)
      vp9_wb_write_literal(wb, cm->reset_frame_context, 2);

    if (cm->intra_only) {
      write_sync_code(wb);

      // Note for profile 0, 420 8bpp is assumed.
      if (cm->profile > PROFILE_0) {
        write_bitdepth_colorspace_sampling(cm, wb);
      }

      vp9_wb_write_literal(wb, get_refresh_mask(cpi), REF_FRAMES);
      write_frame_size(cm, wb);
    } else {
      MV_REFERENCE_FRAME ref_frame;
      vp9_wb_write_literal(wb, get_refresh_mask(cpi), REF_FRAMES);
      for (ref_frame = LAST_FRAME; ref_frame <= ALTREF_FRAME; ++ref_frame) {
        vp9_wb_write_literal(wb, get_ref_frame_idx(cpi, ref_frame),
                             REF_FRAMES_LOG2);
        vp9_wb_write_bit(wb, cm->ref_frame_sign_bias[ref_frame]);
      }

      write_frame_size_with_refs(cpi, wb);

      vp9_wb_write_bit(wb, cm->allow_high_precision_mv);

      fix_interp_filter(cm);
      write_interp_filter(cm->interp_filter, wb);
    }
  }

  if (!cm->error_resilient_mode) {
    vp9_wb_write_bit(wb, cm->refresh_frame_context);
    vp9_wb_write_bit(wb, cm->frame_parallel_decoding_mode);
  }

  vp9_wb_write_literal(wb, cm->frame_context_idx, FRAME_CONTEXTS_LOG2);

  encode_loopfilter(cm, wb);
  encode_quantization(cm, wb);
  encode_segmentation(cm, &cpi->mb.e_mbd, wb);

#if CONFIG_ROW_TILE
  write_tile_info(cpi, wb, final_packing);
#else
  write_tile_info(cm, wb);
#endif
}

#if CONFIG_GLOBAL_MOTION
static void write_global_motion_params(Global_Motion_Params *params,
                                       vp9_prob *probs,
                                       vp9_writer *w) {
  GLOBAL_MOTION_TYPE gmtype = get_gmtype(params);
  vp9_write_token(w, vp9_global_motion_types_tree, probs,
                  &global_motion_types_encodings[gmtype]);
  switch (gmtype) {
    case GLOBAL_ZERO:
      break;
    case GLOBAL_TRANSLATION:
      vp9_write_primitive_symmetric(w, params->mv.as_mv.col,
                                    ABS_TRANSLATION_BITS);
      vp9_write_primitive_symmetric(w, params->mv.as_mv.row,
                                    ABS_TRANSLATION_BITS);
      break;
    case GLOBAL_ROTZOOM:
      vp9_write_primitive_symmetric(w, params->mv.as_mv.col,
                                    ABS_TRANSLATION_BITS);
      vp9_write_primitive_symmetric(w, params->mv.as_mv.row,
                                    ABS_TRANSLATION_BITS);
      vp9_write_primitive_symmetric(w, params->zoom, ABS_ZOOM_BITS);
      vp9_write_primitive_symmetric(w, params->rotation, ABS_ROTATION_BITS);
      break;
    default:
      assert(0);
  }
}

static void write_global_motion(VP9_COMP *cpi, vp9_writer *w) {
  VP9_COMMON *const cm = &cpi->common;
  int frame, i;
  for (frame = LAST_FRAME; frame <= ALTREF_FRAME; ++frame) {
    for (i = 0; i < cm->num_global_motion[frame]; ++i) {
      if (!cpi->global_motion_used[frame]) {
        vpx_memset(
            cm->global_motion[frame], 0,
            MAX_GLOBAL_MOTION_MODELS * sizeof(*cm->global_motion[frame]));
      }
      write_global_motion_params(
          &cm->global_motion[frame][0], cm->fc.global_motion_types_prob, w);
      /*
      printf("Enc Ref %d [%d] (used %d): %d %d %d %d\n",
             frame, cm->current_video_frame, cpi->global_motion_used[frame],
             cm->global_motion[frame][i].zoom,
             cm->global_motion[frame][i].rotation,
             cm->global_motion[frame][i].mv.as_mv.col,
             cm->global_motion[frame][i].mv.as_mv.row);
             */
    }
  }
}
#endif

static size_t write_compressed_header(VP9_COMP *cpi, uint8_t *data) {
  VP9_COMMON *const cm = &cpi->common;
#if !CONFIG_TX_SKIP || CONFIG_SUPERTX
  MACROBLOCKD *const xd = &cpi->mb.e_mbd;
#endif
  FRAME_CONTEXT *const fc = &cm->fc;
  vp9_writer header_bc;

  vp9_start_encode(&header_bc, data);

#if CONFIG_TX_SKIP
  encode_txfm_probs(cm, &header_bc);
#else
  if (xd->lossless)
    cm->tx_mode = ONLY_4X4;
  else
    encode_txfm_probs(cm, &header_bc);
#endif  // CONFIG_TX_SKIP

  update_coef_probs(cpi, &header_bc);
  update_skip_probs(cm, &header_bc);
#if CONFIG_SR_MODE
  update_sr_probs(cm, &header_bc);
#if SR_USE_MULTI_F
  update_sr_usfilter_probs(cm, &header_bc);
#endif  // SR_USE_MULTI_F
#endif  // CONFIG_SR_MODE

  if (!frame_is_intra_only(cm)) {
    int i, j;
    for (i = 0; i < INTER_MODE_CONTEXTS; ++i) {
      prob_diff_update(vp9_inter_mode_tree, cm->fc.inter_mode_probs[i],
                       cm->counts.inter_mode[i], INTER_MODES, &header_bc);
    }

#if CONFIG_NEW_INTER
    update_inter_compound_mode_probs(cm, &header_bc);
#endif  // CONFIG_NEW_INTER

#if CONFIG_NEW_QUANT && QUANT_PROFILES > 1 && !Q_CTX_BASED_PROFILES
    update_dq_profile_probs(cm, &header_bc);
#endif  // CONFIG_NEW_QUANT && QUANT_PROFILES > 1 && !Q_CTX_BASED_PROFILES

    if (cm->interp_filter == SWITCHABLE)
      update_switchable_interp_probs(cm, &header_bc);

    for (i = 0; i < INTRA_INTER_CONTEXTS; i++)
      vp9_cond_prob_diff_update(&header_bc, &fc->intra_inter_prob[i],
                                cm->counts.intra_inter[i]);

    if (cm->allow_comp_inter_inter) {
      const int use_compound_pred = cm->reference_mode != SINGLE_REFERENCE;
      const int use_hybrid_pred =
          cm->reference_mode == REFERENCE_MODE_SELECT;

      vp9_write_bit(&header_bc, use_compound_pred);
      if (use_compound_pred) {
        vp9_write_bit(&header_bc, use_hybrid_pred);

        if (use_hybrid_pred)
          for (i = 0; i < COMP_INTER_CONTEXTS; i++)
            vp9_cond_prob_diff_update(&header_bc, &fc->comp_inter_prob[i],
                                      cm->counts.comp_inter[i]);
      }
    }

    if (cm->reference_mode != COMPOUND_REFERENCE) {
      for (i = 0; i < REF_CONTEXTS; i++) {
        for (j = 0; j < (SINGLE_REFS - 1); j++) {
          vp9_cond_prob_diff_update(&header_bc, &fc->single_ref_probs[i][j],
                                    cm->counts.single_ref[i][j]);
        }
      }
    }

    if (cm->reference_mode != SINGLE_REFERENCE) {
      for (i = 0; i < REF_CONTEXTS; i++) {
        for (j = 0; j < (COMP_REFS - 1); j++) {
          vp9_cond_prob_diff_update(&header_bc, &fc->comp_ref_probs[i][j],
                                    cm->counts.comp_ref[i][j]);
        }
      }
    }

    for (i = 0; i < BLOCK_SIZE_GROUPS; ++i)
      prob_diff_update(vp9_intra_mode_tree, cm->fc.y_mode_prob[i],
                       cm->counts.y_mode[i], INTRA_MODES, &header_bc);

#if CONFIG_EXT_PARTITION
    prob_diff_update(vp9_partition_tree, fc->partition_prob[0],
                     cm->counts.partition[0], PARTITION_TYPES, &header_bc);
    for (i = 1; i < PARTITION_CONTEXTS; ++i)
      prob_diff_update(vp9_ext_partition_tree, fc->partition_prob[i],
                       cm->counts.partition[i], EXT_PARTITION_TYPES,
                       &header_bc);
#else
    for (i = 0; i < PARTITION_CONTEXTS; ++i)
      prob_diff_update(vp9_partition_tree, fc->partition_prob[i],
                       cm->counts.partition[i], PARTITION_TYPES, &header_bc);
#endif

    vp9_write_nmv_probs(cm, cm->allow_high_precision_mv, &header_bc);

#if CONFIG_EXT_TX
    update_ext_tx_probs(cm, &header_bc);
#endif
#if CONFIG_SUPERTX
    if (!xd->lossless)
      update_supertx_probs(cm, &header_bc);
#endif
#if CONFIG_TX_SKIP
    for (i = 0; i < 2; i++)
      vp9_cond_prob_diff_update(&header_bc, &fc->y_tx_skip_prob[i],
                                cm->counts.y_tx_skip[i]);
    for (i = 0; i < 2; i++)
      vp9_cond_prob_diff_update(&header_bc, &fc->uv_tx_skip_prob[i],
                                cm->counts.uv_tx_skip[i]);
#endif
#if CONFIG_COPY_MODE
    for (i = 0; i < COPY_MODE_CONTEXTS; i++) {
      prob_diff_update(vp9_copy_mode_tree_l2, cm->fc.copy_mode_probs_l2[i],
                       cm->counts.copy_mode_l2[i], 2, &header_bc);
      prob_diff_update(vp9_copy_mode_tree, cm->fc.copy_mode_probs[i],
                       cm->counts.copy_mode[i], COPY_MODE_COUNT - 1,
                       &header_bc);
    }
#endif
#if CONFIG_INTERINTRA
    if (cm->reference_mode != COMPOUND_REFERENCE) {
      for (i = 0; i < BLOCK_SIZES; i++) {
        if (is_interintra_allowed(i)) {
          vp9_cond_prob_diff_update(&header_bc,
                                    &fc->interintra_prob[i],
                                    cm->counts.interintra[i]);
        }
      }
#if CONFIG_WEDGE_PARTITION
      for (i = 0; i < BLOCK_SIZES; i++) {
        if (is_interintra_allowed(i) && get_wedge_bits(i))
          vp9_cond_prob_diff_update(&header_bc,
                                    &fc->wedge_interintra_prob[i],
                                    cm->counts.wedge_interintra[i]);
      }
#endif  // CONFIG_WEDGE_PARTITION
    }
#endif  // CONFIG_INTERINTRA
#if CONFIG_WEDGE_PARTITION
    if (cm->reference_mode != SINGLE_REFERENCE) {
      for (i = 0; i < BLOCK_SIZES; i++)
        if (get_wedge_bits(i))
          vp9_cond_prob_diff_update(&header_bc,
                                    &fc->wedge_interinter_prob[i],
                                    cm->counts.wedge_interinter[i]);
    }
#endif  // CONFIG_WEDGE_PARTITION
#if CONFIG_GLOBAL_MOTION
    write_global_motion(cpi, &header_bc);
#endif  // CONFIG_GLOBAL_MOTION
  }

#if CONFIG_INTRABC
  if (frame_is_intra_only(cm))
    vp9_write_bit(&header_bc, cm->allow_intrabc_mode);
#endif  // CONFIG_INTRABC
#if CONFIG_PALETTE
  if (frame_is_intra_only(cm))
    vp9_write_bit(&header_bc, cm->allow_palette_mode);
#endif  // CONFIG_PALETTE

  vp9_stop_encode(&header_bc);
  assert(header_bc.pos <= 0xffff);

  return header_bc.pos;
}

#if CONFIG_ROW_TILE
int vp9_pack_bitstream(VP9_COMP *cpi, uint8_t *dest, size_t *size,
                       int final_packing) {
#else
void vp9_pack_bitstream(VP9_COMP *cpi, uint8_t *dest, size_t *size) {
#endif
  uint8_t *data = dest;
  size_t first_part_size, uncompressed_hdr_size;
  struct vp9_write_bit_buffer wb = {data, 0};
  struct vp9_write_bit_buffer saved_wb;

#if CONFIG_ROW_TILE
  write_uncompressed_header(cpi, &wb, final_packing);
#else
  write_uncompressed_header(cpi, &wb);
#endif

  saved_wb = wb;
  vp9_wb_write_literal(&wb, 0, 16);  // don't know in advance first part. size

  uncompressed_hdr_size = vp9_wb_bytes_written(&wb);
  data += uncompressed_hdr_size;

  vp9_clear_system_state();

  first_part_size = write_compressed_header(cpi, data);
  data += first_part_size;
  // TODO(jbb): Figure out what to do if first_part_size > 16 bits.
  vp9_wb_write_literal(&saved_wb, (int)first_part_size, 16);

#if CONFIG_ROW_TILE
  {
    size_t total_size = encode_tiles(cpi, data, final_packing);
    // If tile size is out of the range, exit and do another packing.
    if (total_size == 0)
      return -1;

    data += total_size;
  }
#else
  data += encode_tiles(cpi, data);
#endif

  *size = data - dest;
#if CONFIG_ROW_TILE
  return 0;
#endif
}
