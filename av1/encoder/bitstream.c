/*
 * Copyright (c) 2016, Alliance for Open Media. All rights reserved
 *
 * This source code is subject to the terms of the BSD 2 Clause License and
 * the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
 * was not distributed with this source code in the LICENSE file, you can
 * obtain it at www.aomedia.org/license/software. If the Alliance for Open
 * Media Patent License 1.0 was not distributed with this source code in the
 * PATENTS file, you can obtain it at www.aomedia.org/license/patent.
 */

#include <assert.h>
#include <limits.h>
#include <stdio.h>

#include "aom/aom_encoder.h"
#include "aom_dsp/bitwriter_buffer.h"
#include "aom_dsp/aom_dsp_common.h"
#include "aom_mem/aom_mem.h"
#include "aom_ports/mem_ops.h"
#include "aom_ports/system_state.h"
#if CONFIG_BITSTREAM_DEBUG
#include "aom_util/debug_util.h"
#endif  // CONFIG_BITSTREAM_DEBUG

#if CONFIG_CLPF
#include "av1/common/clpf.h"
#endif
#if CONFIG_DERING
#include "av1/common/dering.h"
#endif  // CONFIG_DERING
#include "av1/common/entropy.h"
#include "av1/common/entropymode.h"
#include "av1/common/entropymv.h"
#include "av1/common/mvref_common.h"
#include "av1/common/odintrin.h"
#include "av1/common/pred_common.h"
#include "av1/common/reconinter.h"
#include "av1/common/seg_common.h"
#include "av1/common/tile_common.h"

#if CONFIG_ANS
#include "aom_dsp/buf_ans.h"
#endif  // CONFIG_ANS
#include "av1/encoder/bitstream.h"
#include "av1/encoder/cost.h"
#include "av1/encoder/encodemv.h"
#include "av1/encoder/mcomp.h"
#include "av1/encoder/segmentation.h"
#include "av1/encoder/subexp.h"
#include "av1/encoder/tokenize.h"
#if CONFIG_PVQ
#include "av1/encoder/pvq_encoder.h"
#endif

static struct av1_token intra_mode_encodings[INTRA_MODES];
static struct av1_token switchable_interp_encodings[SWITCHABLE_FILTERS];
#if CONFIG_EXT_PARTITION_TYPES
static const struct av1_token ext_partition_encodings[EXT_PARTITION_TYPES] = {
  { 0, 1 },  { 4, 3 },  { 12, 4 }, { 7, 3 },
  { 10, 4 }, { 11, 4 }, { 26, 5 }, { 27, 5 }
};
#endif
static struct av1_token partition_encodings[PARTITION_TYPES];
#if !CONFIG_REF_MV
static struct av1_token inter_mode_encodings[INTER_MODES];
#endif
#if CONFIG_EXT_INTER
static const struct av1_token
    inter_compound_mode_encodings[INTER_COMPOUND_MODES] = {
      { 2, 2 },  { 50, 6 }, { 51, 6 }, { 24, 5 }, { 52, 6 },
      { 53, 6 }, { 54, 6 }, { 55, 6 }, { 0, 1 },  { 7, 3 }
    };
#endif  // CONFIG_EXT_INTER
#if CONFIG_PALETTE
static struct av1_token palette_size_encodings[PALETTE_MAX_SIZE - 1];
static struct av1_token palette_color_encodings[PALETTE_MAX_SIZE - 1]
                                               [PALETTE_MAX_SIZE];
#endif  // CONFIG_PALETTE
static const struct av1_token tx_size_encodings[MAX_TX_DEPTH][TX_SIZES] = {
  { { 0, 1 }, { 1, 1 } },                      // Max tx_size is 8X8
  { { 0, 1 }, { 2, 2 }, { 3, 2 } },            // Max tx_size is 16X16
  { { 0, 1 }, { 2, 2 }, { 6, 3 }, { 7, 3 } },  // Max tx_size is 32X32
};

#if CONFIG_EXT_INTRA || CONFIG_FILTER_INTRA || CONFIG_PALETTE
static INLINE void write_uniform(aom_writer *w, int n, int v) {
  int l = get_unsigned_bits(n);
  int m = (1 << l) - n;
  if (l == 0) return;
  if (v < m) {
    aom_write_literal(w, v, l - 1);
  } else {
    aom_write_literal(w, m + ((v - m) >> 1), l - 1);
    aom_write_literal(w, (v - m) & 1, 1);
  }
}
#endif  // CONFIG_EXT_INTRA || CONFIG_FILTER_INTRA || CONFIG_PALETTE

#if CONFIG_EXT_TX
static struct av1_token ext_tx_inter_encodings[EXT_TX_SETS_INTER][TX_TYPES];
static struct av1_token ext_tx_intra_encodings[EXT_TX_SETS_INTRA][TX_TYPES];
#else
static struct av1_token ext_tx_encodings[TX_TYPES];
#endif  // CONFIG_EXT_TX
#if CONFIG_GLOBAL_MOTION
static struct av1_token global_motion_types_encodings[GLOBAL_MOTION_TYPES];
#endif  // CONFIG_GLOBAL_MOTION
#if CONFIG_EXT_INTRA
static struct av1_token intra_filter_encodings[INTRA_FILTERS];
#endif  // CONFIG_EXT_INTRA
#if CONFIG_EXT_INTER
static struct av1_token interintra_mode_encodings[INTERINTRA_MODES];
#endif  // CONFIG_EXT_INTER
#if CONFIG_MOTION_VAR || CONFIG_WARPED_MOTION
static struct av1_token motion_mode_encodings[MOTION_MODES];
#endif  // CONFIG_MOTION_VAR || CONFIG_WARPED_MOTION
#if CONFIG_LOOP_RESTORATION
static struct av1_token switchable_restore_encodings[RESTORE_SWITCHABLE_TYPES];
#endif  // CONFIG_LOOP_RESTORATION
static void write_uncompressed_header(AV1_COMP *cpi,
                                      struct aom_write_bit_buffer *wb);
static uint32_t write_compressed_header(AV1_COMP *cpi, uint8_t *data);

void av1_encode_token_init(void) {
#if CONFIG_EXT_TX || CONFIG_PALETTE
  int s;
#endif  // CONFIG_EXT_TX || CONFIG_PALETTE
#if CONFIG_EXT_TX
  for (s = 1; s < EXT_TX_SETS_INTER; ++s) {
    av1_tokens_from_tree(ext_tx_inter_encodings[s], av1_ext_tx_inter_tree[s]);
  }
  for (s = 1; s < EXT_TX_SETS_INTRA; ++s) {
    av1_tokens_from_tree(ext_tx_intra_encodings[s], av1_ext_tx_intra_tree[s]);
  }
#else
  av1_tokens_from_tree(ext_tx_encodings, av1_ext_tx_tree);
#endif  // CONFIG_EXT_TX
  av1_tokens_from_tree(intra_mode_encodings, av1_intra_mode_tree);
  av1_tokens_from_tree(switchable_interp_encodings, av1_switchable_interp_tree);
  av1_tokens_from_tree(partition_encodings, av1_partition_tree);
#if !CONFIG_REF_MV
  av1_tokens_from_tree(inter_mode_encodings, av1_inter_mode_tree);
#endif

#if CONFIG_PALETTE
  av1_tokens_from_tree(palette_size_encodings, av1_palette_size_tree);
  for (s = 0; s < PALETTE_MAX_SIZE - 1; ++s) {
    av1_tokens_from_tree(palette_color_encodings[s], av1_palette_color_tree[s]);
  }
#endif  // CONFIG_PALETTE

#if CONFIG_EXT_INTRA
  av1_tokens_from_tree(intra_filter_encodings, av1_intra_filter_tree);
#endif  // CONFIG_EXT_INTRA
#if CONFIG_EXT_INTER
  av1_tokens_from_tree(interintra_mode_encodings, av1_interintra_mode_tree);
#endif  // CONFIG_EXT_INTER
#if CONFIG_MOTION_VAR || CONFIG_WARPED_MOTION
  av1_tokens_from_tree(motion_mode_encodings, av1_motion_mode_tree);
#endif  // CONFIG_MOTION_VAR || CONFIG_WARPED_MOTION
#if CONFIG_GLOBAL_MOTION
  av1_tokens_from_tree(global_motion_types_encodings,
                       av1_global_motion_types_tree);
#endif  // CONFIG_GLOBAL_MOTION
#if CONFIG_LOOP_RESTORATION
  av1_tokens_from_tree(switchable_restore_encodings,
                       av1_switchable_restore_tree);
#endif  // CONFIG_LOOP_RESTORATION

#if CONFIG_DAALA_EC
  /* This hack is necessary when CONFIG_EXT_INTERP is enabled because the five
      SWITCHABLE_FILTERS are not consecutive, e.g., 0, 1, 2, 3, 4, when doing
      an in-order traversal of the av1_switchable_interp_tree structure. */
  av1_indices_from_tree(av1_switchable_interp_ind, av1_switchable_interp_inv,
                        SWITCHABLE_FILTERS, av1_switchable_interp_tree);
  /* This hack is necessary because the four TX_TYPES are not consecutive,
      e.g., 0, 1, 2, 3, when doing an in-order traversal of the av1_ext_tx_tree
      structure. */
  av1_indices_from_tree(av1_ext_tx_ind, av1_ext_tx_inv, TX_TYPES,
                        av1_ext_tx_tree);
  av1_indices_from_tree(av1_intra_mode_ind, av1_intra_mode_inv, INTRA_MODES,
                        av1_intra_mode_tree);
  av1_indices_from_tree(av1_inter_mode_ind, av1_inter_mode_inv, INTER_MODES,
                        av1_inter_mode_tree);
#endif
}

#if !CONFIG_DAALA_EC
static void write_intra_mode(aom_writer *w, PREDICTION_MODE mode,
                             const aom_prob *probs) {
  av1_write_token(w, av1_intra_mode_tree, probs, &intra_mode_encodings[mode]);
}
#endif

#if CONFIG_EXT_INTER
static void write_interintra_mode(aom_writer *w, INTERINTRA_MODE mode,
                                  const aom_prob *probs) {
  av1_write_token(w, av1_interintra_mode_tree, probs,
                  &interintra_mode_encodings[mode]);
}
#endif  // CONFIG_EXT_INTER

static void write_inter_mode(AV1_COMMON *cm, aom_writer *w,
                             PREDICTION_MODE mode,
#if CONFIG_REF_MV && CONFIG_EXT_INTER
                             int is_compound,
#endif  // CONFIG_REF_MV && CONFIG_EXT_INTER
                             const int16_t mode_ctx) {
#if CONFIG_REF_MV
  const int16_t newmv_ctx = mode_ctx & NEWMV_CTX_MASK;
  const aom_prob newmv_prob = cm->fc->newmv_prob[newmv_ctx];
#if CONFIG_EXT_INTER
  aom_write(w, mode != NEWMV && mode != NEWFROMNEARMV, newmv_prob);

  if (!is_compound && (mode == NEWMV || mode == NEWFROMNEARMV))
    aom_write(w, mode == NEWFROMNEARMV, cm->fc->new2mv_prob);

  if (mode != NEWMV && mode != NEWFROMNEARMV) {
#else
  aom_write(w, mode != NEWMV, newmv_prob);

  if (mode != NEWMV) {
#endif  // CONFIG_EXT_INTER
    const int16_t zeromv_ctx = (mode_ctx >> ZEROMV_OFFSET) & ZEROMV_CTX_MASK;
    const aom_prob zeromv_prob = cm->fc->zeromv_prob[zeromv_ctx];

    if (mode_ctx & (1 << ALL_ZERO_FLAG_OFFSET)) {
      assert(mode == ZEROMV);
      return;
    }

    aom_write(w, mode != ZEROMV, zeromv_prob);

    if (mode != ZEROMV) {
      int16_t refmv_ctx = (mode_ctx >> REFMV_OFFSET) & REFMV_CTX_MASK;
      aom_prob refmv_prob;

      if (mode_ctx & (1 << SKIP_NEARESTMV_OFFSET)) refmv_ctx = 6;
      if (mode_ctx & (1 << SKIP_NEARMV_OFFSET)) refmv_ctx = 7;
      if (mode_ctx & (1 << SKIP_NEARESTMV_SUB8X8_OFFSET)) refmv_ctx = 8;

      refmv_prob = cm->fc->refmv_prob[refmv_ctx];
      aom_write(w, mode != NEARESTMV, refmv_prob);
    }
  }
#else
  assert(is_inter_mode(mode));
#if CONFIG_DAALA_EC
  aom_write_symbol(w, av1_inter_mode_ind[INTER_OFFSET(mode)],
                   cm->fc->inter_mode_cdf[mode_ctx], INTER_MODES);
#else
  {
    const aom_prob *const inter_probs = cm->fc->inter_mode_probs[mode_ctx];
    av1_write_token(w, av1_inter_mode_tree, inter_probs,
                    &inter_mode_encodings[INTER_OFFSET(mode)]);
  }
#endif
#endif
}

#if CONFIG_REF_MV
static void write_drl_idx(const AV1_COMMON *cm, const MB_MODE_INFO *mbmi,
                          const MB_MODE_INFO_EXT *mbmi_ext, aom_writer *w) {
  uint8_t ref_frame_type = av1_ref_frame_type(mbmi->ref_frame);

  assert(mbmi->ref_mv_idx < 3);

  if (mbmi->mode == NEWMV) {
    int idx;
    for (idx = 0; idx < 2; ++idx) {
      if (mbmi_ext->ref_mv_count[ref_frame_type] > idx + 1) {
        uint8_t drl_ctx =
            av1_drl_ctx(mbmi_ext->ref_mv_stack[ref_frame_type], idx);
        aom_prob drl_prob = cm->fc->drl_prob[drl_ctx];

        aom_write(w, mbmi->ref_mv_idx != idx, drl_prob);
        if (mbmi->ref_mv_idx == idx) return;
      }
    }
    return;
  }

  if (mbmi->mode == NEARMV) {
    int idx;
    // TODO(jingning): Temporary solution to compensate the NEARESTMV offset.
    for (idx = 1; idx < 3; ++idx) {
      if (mbmi_ext->ref_mv_count[ref_frame_type] > idx + 1) {
        uint8_t drl_ctx =
            av1_drl_ctx(mbmi_ext->ref_mv_stack[ref_frame_type], idx);
        aom_prob drl_prob = cm->fc->drl_prob[drl_ctx];

        aom_write(w, mbmi->ref_mv_idx != (idx - 1), drl_prob);
        if (mbmi->ref_mv_idx == (idx - 1)) return;
      }
    }
    return;
  }
}
#endif

#if CONFIG_EXT_INTER
static void write_inter_compound_mode(AV1_COMMON *cm, aom_writer *w,
                                      PREDICTION_MODE mode,
                                      const int16_t mode_ctx) {
  const aom_prob *const inter_compound_probs =
      cm->fc->inter_compound_mode_probs[mode_ctx];

  assert(is_inter_compound_mode(mode));
  av1_write_token(w, av1_inter_compound_mode_tree, inter_compound_probs,
                  &inter_compound_mode_encodings[INTER_COMPOUND_OFFSET(mode)]);
}
#endif  // CONFIG_EXT_INTER

static void encode_unsigned_max(struct aom_write_bit_buffer *wb, int data,
                                int max) {
  aom_wb_write_literal(wb, data, get_unsigned_bits(max));
}

static void prob_diff_update(const aom_tree_index *tree,
                             aom_prob probs[/*n - 1*/],
                             const unsigned int counts[/*n - 1*/], int n,
                             int probwt, aom_writer *w) {
  int i;
  unsigned int branch_ct[32][2];

  // Assuming max number of probabilities <= 32
  assert(n <= 32);

  av1_tree_probs_from_distribution(tree, branch_ct, counts);
  for (i = 0; i < n - 1; ++i)
    av1_cond_prob_diff_update(w, &probs[i], branch_ct[i], probwt);
}

#if !CONFIG_EC_ADAPT
static int prob_diff_update_savings(const aom_tree_index *tree,
                                    aom_prob probs[/*n - 1*/],
                                    const unsigned int counts[/*n - 1*/], int n,
                                    int probwt) {
  int i;
  unsigned int branch_ct[32][2];
  int savings = 0;

  // Assuming max number of probabilities <= 32
  assert(n <= 32);
  av1_tree_probs_from_distribution(tree, branch_ct, counts);
  for (i = 0; i < n - 1; ++i) {
    savings +=
        av1_cond_prob_diff_update_savings(&probs[i], branch_ct[i], probwt);
  }
  return savings;
}
#endif

#if CONFIG_VAR_TX
static void write_tx_size_vartx(const AV1_COMMON *cm, const MACROBLOCKD *xd,
                                const MB_MODE_INFO *mbmi, TX_SIZE tx_size,
                                int depth, int blk_row, int blk_col,
                                aom_writer *w) {
  const int tx_row = blk_row >> 1;
  const int tx_col = blk_col >> 1;
  const int max_blocks_high = max_block_high(xd, mbmi->sb_type, 0);
  const int max_blocks_wide = max_block_wide(xd, mbmi->sb_type, 0);

  int ctx = txfm_partition_context(xd->above_txfm_context + tx_col,
                                   xd->left_txfm_context + tx_row,
                                   mbmi->sb_type, tx_size);

  if (blk_row >= max_blocks_high || blk_col >= max_blocks_wide) return;

  if (depth == MAX_VARTX_DEPTH) {
    txfm_partition_update(xd->above_txfm_context + tx_col,
                          xd->left_txfm_context + tx_row, tx_size);
    return;
  }

  if (tx_size == mbmi->inter_tx_size[tx_row][tx_col]) {
    aom_write(w, 0, cm->fc->txfm_partition_prob[ctx]);
    txfm_partition_update(xd->above_txfm_context + tx_col,
                          xd->left_txfm_context + tx_row, tx_size);
  } else {
    const TX_SIZE sub_txs = sub_tx_size_map[tx_size];
    const int bsl = tx_size_wide_unit[sub_txs];
    int i;

    aom_write(w, 1, cm->fc->txfm_partition_prob[ctx]);

    if (tx_size == TX_8X8) {
      txfm_partition_update(xd->above_txfm_context + tx_col,
                            xd->left_txfm_context + tx_row, TX_4X4);
      return;
    }

    assert(bsl > 0);
    for (i = 0; i < 4; ++i) {
      int offsetr = blk_row + (i >> 1) * bsl;
      int offsetc = blk_col + (i & 0x01) * bsl;
      write_tx_size_vartx(cm, xd, mbmi, sub_txs, depth + 1, offsetr, offsetc,
                          w);
    }
  }
}

static void update_txfm_partition_probs(AV1_COMMON *cm, aom_writer *w,
                                        FRAME_COUNTS *counts, int probwt) {
  int k;
  for (k = 0; k < TXFM_PARTITION_CONTEXTS; ++k)
    av1_cond_prob_diff_update(w, &cm->fc->txfm_partition_prob[k],
                              counts->txfm_partition[k], probwt);
}
#endif

static void write_selected_tx_size(const AV1_COMMON *cm, const MACROBLOCKD *xd,
                                   aom_writer *w) {
  const MB_MODE_INFO *const mbmi = &xd->mi[0]->mbmi;
  const BLOCK_SIZE bsize = mbmi->sb_type;
  // For sub8x8 blocks the tx_size symbol does not need to be sent
  if (bsize >= BLOCK_8X8) {
    const TX_SIZE tx_size = mbmi->tx_size;
    const int is_inter = is_inter_block(mbmi);
    const int tx_size_ctx = get_tx_size_context(xd);
    const int tx_size_cat = is_inter ? inter_tx_size_cat_lookup[bsize]
                                     : intra_tx_size_cat_lookup[bsize];
    const TX_SIZE coded_tx_size = txsize_sqr_up_map[tx_size];
    const int depth = tx_size_to_depth(coded_tx_size);

#if CONFIG_EXT_TX && CONFIG_RECT_TX
    assert(IMPLIES(is_rect_tx(tx_size), is_rect_tx_allowed(xd, mbmi)));
    assert(
        IMPLIES(is_rect_tx(tx_size), tx_size == max_txsize_rect_lookup[bsize]));
#endif  // CONFIG_EXT_TX && CONFIG_RECT_TX

    av1_write_token(w, av1_tx_size_tree[tx_size_cat],
                    cm->fc->tx_size_probs[tx_size_cat][tx_size_ctx],
                    &tx_size_encodings[tx_size_cat][depth]);
  }
}

#if CONFIG_REF_MV
static void update_inter_mode_probs(AV1_COMMON *cm, aom_writer *w,
                                    FRAME_COUNTS *counts) {
  int i;
#if CONFIG_TILE_GROUPS
  const int probwt = cm->num_tg;
#else
  const int probwt = 1;
#endif
  for (i = 0; i < NEWMV_MODE_CONTEXTS; ++i)
    av1_cond_prob_diff_update(w, &cm->fc->newmv_prob[i], counts->newmv_mode[i],
                              probwt);
  for (i = 0; i < ZEROMV_MODE_CONTEXTS; ++i)
    av1_cond_prob_diff_update(w, &cm->fc->zeromv_prob[i],
                              counts->zeromv_mode[i], probwt);
  for (i = 0; i < REFMV_MODE_CONTEXTS; ++i)
    av1_cond_prob_diff_update(w, &cm->fc->refmv_prob[i], counts->refmv_mode[i],
                              probwt);
  for (i = 0; i < DRL_MODE_CONTEXTS; ++i)
    av1_cond_prob_diff_update(w, &cm->fc->drl_prob[i], counts->drl_mode[i],
                              probwt);
#if CONFIG_EXT_INTER
  av1_cond_prob_diff_update(w, &cm->fc->new2mv_prob, counts->new2mv_mode,
                            probwt);
#endif  // CONFIG_EXT_INTER
}
#endif

#if CONFIG_EXT_INTER
static void update_inter_compound_mode_probs(AV1_COMMON *cm, int probwt,
                                             aom_writer *w) {
  const int savings_thresh = av1_cost_one(GROUP_DIFF_UPDATE_PROB) -
                             av1_cost_zero(GROUP_DIFF_UPDATE_PROB);
  int i;
  int savings = 0;
  int do_update = 0;
  for (i = 0; i < INTER_MODE_CONTEXTS; ++i) {
    savings += prob_diff_update_savings(
        av1_inter_compound_mode_tree, cm->fc->inter_compound_mode_probs[i],
        cm->counts.inter_compound_mode[i], INTER_COMPOUND_MODES, probwt);
  }
  do_update = savings > savings_thresh;
  aom_write(w, do_update, GROUP_DIFF_UPDATE_PROB);
  if (do_update) {
    for (i = 0; i < INTER_MODE_CONTEXTS; ++i) {
      prob_diff_update(
          av1_inter_compound_mode_tree, cm->fc->inter_compound_mode_probs[i],
          cm->counts.inter_compound_mode[i], INTER_COMPOUND_MODES, probwt, w);
    }
  }
}
#endif  // CONFIG_EXT_INTER

static int write_skip(const AV1_COMMON *cm, const MACROBLOCKD *xd,
                      int segment_id, const MODE_INFO *mi, aom_writer *w) {
  if (segfeature_active(&cm->seg, segment_id, SEG_LVL_SKIP)) {
    return 1;
  } else {
    const int skip = mi->mbmi.skip;
    aom_write(w, skip, av1_get_skip_prob(cm, xd));
    return skip;
  }
}

#if CONFIG_DELTA_Q
static void write_delta_qindex(const AV1_COMMON *cm, int delta_qindex,
                               aom_writer *w) {
  int sign = delta_qindex < 0;
  int abs = sign ? -delta_qindex : delta_qindex;
  int rem_bits, thr, i = 0;
  int smallval = abs < DELTA_Q_SMALL ? 1 : 0;

  while (i < DELTA_Q_SMALL && i <= abs) {
    int bit = (i < abs);
    aom_write(w, bit, cm->fc->delta_q_prob[i]);
    i++;
  }

  if (!smallval) {
    rem_bits = OD_ILOG_NZ(abs - 1) - 1;
    thr = (1 << rem_bits) + 1;
    aom_write_literal(w, rem_bits, 3);
    aom_write_literal(w, abs - thr, rem_bits);
  }
  if (abs > 0) {
    aom_write_bit(w, sign);
  }
}

static void update_delta_q_probs(AV1_COMMON *cm, aom_writer *w,
                                 FRAME_COUNTS *counts) {
  int k;
#if CONFIG_TILE_GROUPS
  const int probwt = cm->num_tg;
#else
  const int probwt = 1;
#endif
  for (k = 0; k < DELTA_Q_CONTEXTS; ++k) {
    av1_cond_prob_diff_update(w, &cm->fc->delta_q_prob[k], counts->delta_q[k],
                              probwt);
  }
}
#endif

static void update_skip_probs(AV1_COMMON *cm, aom_writer *w,
                              FRAME_COUNTS *counts) {
  int k;
#if CONFIG_TILE_GROUPS
  const int probwt = cm->num_tg;
#else
  const int probwt = 1;
#endif
  for (k = 0; k < SKIP_CONTEXTS; ++k) {
    av1_cond_prob_diff_update(w, &cm->fc->skip_probs[k], counts->skip[k],
                              probwt);
  }
}

#if !CONFIG_EC_ADAPT
static void update_switchable_interp_probs(AV1_COMMON *cm, aom_writer *w,
                                           FRAME_COUNTS *counts) {
  int j;
  for (j = 0; j < SWITCHABLE_FILTER_CONTEXTS; ++j) {
#if CONFIG_TILE_GROUPS
    const int probwt = cm->num_tg;
#else
    const int probwt = 1;
#endif
    prob_diff_update(
        av1_switchable_interp_tree, cm->fc->switchable_interp_prob[j],
        counts->switchable_interp[j], SWITCHABLE_FILTERS, probwt, w);
  }
}
#endif

#if CONFIG_EXT_TX
static void update_ext_tx_probs(AV1_COMMON *cm, aom_writer *w) {
  const int savings_thresh = av1_cost_one(GROUP_DIFF_UPDATE_PROB) -
                             av1_cost_zero(GROUP_DIFF_UPDATE_PROB);
  int i, j;
  int s;
#if CONFIG_TILE_GROUPS
  const int probwt = cm->num_tg;
#else
  const int probwt = 1;
#endif
  for (s = 1; s < EXT_TX_SETS_INTER; ++s) {
    int savings = 0;
    int do_update = 0;
    for (i = TX_4X4; i < EXT_TX_SIZES; ++i) {
      if (!use_inter_ext_tx_for_txsize[s][i]) continue;
      savings += prob_diff_update_savings(
          av1_ext_tx_inter_tree[s], cm->fc->inter_ext_tx_prob[s][i],
          cm->counts.inter_ext_tx[s][i], num_ext_tx_set_inter[s], probwt);
    }
    do_update = savings > savings_thresh;
    aom_write(w, do_update, GROUP_DIFF_UPDATE_PROB);
    if (do_update) {
      for (i = TX_4X4; i < EXT_TX_SIZES; ++i) {
        if (!use_inter_ext_tx_for_txsize[s][i]) continue;
        prob_diff_update(
            av1_ext_tx_inter_tree[s], cm->fc->inter_ext_tx_prob[s][i],
            cm->counts.inter_ext_tx[s][i], num_ext_tx_set_inter[s], probwt, w);
      }
    }
  }

  for (s = 1; s < EXT_TX_SETS_INTRA; ++s) {
    int savings = 0;
    int do_update = 0;
    for (i = TX_4X4; i < EXT_TX_SIZES; ++i) {
      if (!use_intra_ext_tx_for_txsize[s][i]) continue;
      for (j = 0; j < INTRA_MODES; ++j)
        savings += prob_diff_update_savings(
            av1_ext_tx_intra_tree[s], cm->fc->intra_ext_tx_prob[s][i][j],
            cm->counts.intra_ext_tx[s][i][j], num_ext_tx_set_intra[s], probwt);
    }
    do_update = savings > savings_thresh;
    aom_write(w, do_update, GROUP_DIFF_UPDATE_PROB);
    if (do_update) {
      for (i = TX_4X4; i < EXT_TX_SIZES; ++i) {
        if (!use_intra_ext_tx_for_txsize[s][i]) continue;
        for (j = 0; j < INTRA_MODES; ++j)
          prob_diff_update(av1_ext_tx_intra_tree[s],
                           cm->fc->intra_ext_tx_prob[s][i][j],
                           cm->counts.intra_ext_tx[s][i][j],
                           num_ext_tx_set_intra[s], probwt, w);
      }
    }
  }
}

#else
#if !CONFIG_EC_ADAPT
static void update_ext_tx_probs(AV1_COMMON *cm, aom_writer *w) {
  const int savings_thresh = av1_cost_one(GROUP_DIFF_UPDATE_PROB) -
                             av1_cost_zero(GROUP_DIFF_UPDATE_PROB);
  int i, j;

  int savings = 0;
  int do_update = 0;
#if CONFIG_TILE_GROUPS
  const int probwt = cm->num_tg;
#else
  const int probwt = 1;
#endif
  for (i = TX_4X4; i < EXT_TX_SIZES; ++i) {
    for (j = 0; j < TX_TYPES; ++j)
      savings += prob_diff_update_savings(
          av1_ext_tx_tree, cm->fc->intra_ext_tx_prob[i][j],
          cm->counts.intra_ext_tx[i][j], TX_TYPES, probwt);
  }
  do_update = savings > savings_thresh;
  aom_write(w, do_update, GROUP_DIFF_UPDATE_PROB);
  if (do_update) {
    for (i = TX_4X4; i < EXT_TX_SIZES; ++i) {
      for (j = 0; j < TX_TYPES; ++j) {
        prob_diff_update(av1_ext_tx_tree, cm->fc->intra_ext_tx_prob[i][j],
                         cm->counts.intra_ext_tx[i][j], TX_TYPES, probwt, w);
      }
    }
  }

  savings = 0;
  for (i = TX_4X4; i < EXT_TX_SIZES; ++i) {
    savings +=
        prob_diff_update_savings(av1_ext_tx_tree, cm->fc->inter_ext_tx_prob[i],
                                 cm->counts.inter_ext_tx[i], TX_TYPES, probwt);
  }
  do_update = savings > savings_thresh;
  aom_write(w, do_update, GROUP_DIFF_UPDATE_PROB);
  if (do_update) {
    for (i = TX_4X4; i < EXT_TX_SIZES; ++i) {
      prob_diff_update(av1_ext_tx_tree, cm->fc->inter_ext_tx_prob[i],
                       cm->counts.inter_ext_tx[i], TX_TYPES, probwt, w);
    }
  }
}
#endif  // CONFIG_EXT_TX
#endif
#if CONFIG_PALETTE
static void pack_palette_tokens(aom_writer *w, const TOKENEXTRA **tp, int n,
                                int num) {
  int i;
  const TOKENEXTRA *p = *tp;

  for (i = 0; i < num; ++i) {
    av1_write_token(w, av1_palette_color_tree[n - 2], p->context_tree,
                    &palette_color_encodings[n - 2][p->token]);
    ++p;
  }

  *tp = p;
}
#endif  // CONFIG_PALETTE
#if !CONFIG_PVQ
#if CONFIG_SUPERTX
static void update_supertx_probs(AV1_COMMON *cm, int probwt, aom_writer *w) {
  const int savings_thresh = av1_cost_one(GROUP_DIFF_UPDATE_PROB) -
                             av1_cost_zero(GROUP_DIFF_UPDATE_PROB);
  int i, j;
  int savings = 0;
  int do_update = 0;
  for (i = 0; i < PARTITION_SUPERTX_CONTEXTS; ++i) {
    for (j = 1; j < TX_SIZES; ++j) {
      savings += av1_cond_prob_diff_update_savings(
          &cm->fc->supertx_prob[i][j], cm->counts.supertx[i][j], probwt);
    }
  }
  do_update = savings > savings_thresh;
  aom_write(w, do_update, GROUP_DIFF_UPDATE_PROB);
  if (do_update) {
    for (i = 0; i < PARTITION_SUPERTX_CONTEXTS; ++i) {
      for (j = 1; j < TX_SIZES; ++j) {
        av1_cond_prob_diff_update(w, &cm->fc->supertx_prob[i][j],
                                  cm->counts.supertx[i][j], probwt);
      }
    }
  }
}
#endif  // CONFIG_SUPERTX

static void pack_mb_tokens(aom_writer *w, const TOKENEXTRA **tp,
                           const TOKENEXTRA *const stop,
                           aom_bit_depth_t bit_depth, const TX_SIZE tx_size,
                           TOKEN_STATS *token_stats) {
  const TOKENEXTRA *p = *tp;
#if CONFIG_VAR_TX
  int count = 0;
  const int seg_eob = tx_size_2d[tx_size];
#endif
#if CONFIG_AOM_HIGHBITDEPTH
  const av1_extra_bit *const extra_bits_table =
      (bit_depth == AOM_BITS_12)
          ? av1_extra_bits_high12
          : (bit_depth == AOM_BITS_10) ? av1_extra_bits_high10 : av1_extra_bits;
#else
  const av1_extra_bit *const extra_bits_table = av1_extra_bits;
  (void)bit_depth;
#endif  // CONFIG_AOM_HIGHBITDEPTH

  while (p < stop && p->token != EOSB_TOKEN) {
    const int token = p->token;
    aom_tree_index index = 0;
#if !CONFIG_EC_MULTISYMBOL
    const struct av1_token *const coef_encoding = &av1_coef_encodings[token];
    int coef_value = coef_encoding->value;
    int coef_length = coef_encoding->len;
#endif  // !CONFIG_EC_MULTISYMBOL
    const av1_extra_bit *const extra_bits = &extra_bits_table[token];

#if CONFIG_EC_MULTISYMBOL
    /* skip one or two nodes */
    if (!p->skip_eob_node)
      aom_write_record(w, token != EOB_TOKEN, p->context_tree[0], token_stats);

    if (token != EOB_TOKEN) {
      aom_write_record(w, token != ZERO_TOKEN, p->context_tree[1], token_stats);

      if (token != ZERO_TOKEN) {
        aom_write_symbol(w, token - ONE_TOKEN, *p->token_cdf,
                         CATEGORY6_TOKEN - ONE_TOKEN + 1);
      }
    }
#else
    /* skip one or two nodes */
    if (p->skip_eob_node)
      coef_length -= p->skip_eob_node;
    else
      aom_write_record(w, token != EOB_TOKEN, p->context_tree[0], token_stats);

    if (token != EOB_TOKEN) {
      aom_write_record(w, token != ZERO_TOKEN, p->context_tree[1], token_stats);

      if (token != ZERO_TOKEN) {
        aom_write_record(w, token != ONE_TOKEN, p->context_tree[2],
                         token_stats);

        if (token != ONE_TOKEN) {
          const int unconstrained_len = UNCONSTRAINED_NODES - p->skip_eob_node;
          aom_write_tree_record(
              w, av1_coef_con_tree,
              av1_pareto8_full[p->context_tree[PIVOT_NODE] - 1], coef_value,
              coef_length - unconstrained_len, 0, token_stats);
        }
      }
    }
#endif  // CONFIG_EC_MULTISYMBOL

    if (extra_bits->base_val) {
      const int bit_string = p->extra;
      const int bit_string_length = extra_bits->len;  // Length of extra bits to
                                                      // be written excluding
                                                      // the sign bit.
      int skip_bits = (extra_bits->base_val == CAT6_MIN_VAL)
                          ? TX_SIZES - 1 - txsize_sqr_up_map[tx_size]
                          : 0;

      if (bit_string_length > 0) {
        const unsigned char *pb = extra_bits->prob;
        const int value = bit_string >> 1;
        const int num_bits = bit_string_length;  // number of bits in value
        assert(num_bits > 0);

        for (index = 0; index < num_bits; ++index) {
          const int shift = num_bits - index - 1;
          const int bb = (value >> shift) & 1;
          if (skip_bits) {
            --skip_bits;
            assert(!bb);
          } else {
            aom_write_record(w, bb, pb[index], token_stats);
          }
        }
      }

      aom_write_bit_record(w, bit_string & 1, token_stats);
    }
    ++p;

#if CONFIG_VAR_TX
    ++count;
    if (token == EOB_TOKEN || count == seg_eob) break;
#endif
  }

  *tp = p;
}
#endif  // !CONFIG_PVG
#if CONFIG_VAR_TX
static void pack_txb_tokens(aom_writer *w, const TOKENEXTRA **tp,
                            const TOKENEXTRA *const tok_end, MACROBLOCKD *xd,
                            MB_MODE_INFO *mbmi, int plane,
                            BLOCK_SIZE plane_bsize, aom_bit_depth_t bit_depth,
                            int block, int blk_row, int blk_col,
                            TX_SIZE tx_size, TOKEN_STATS *token_stats) {
  const struct macroblockd_plane *const pd = &xd->plane[plane];
  const BLOCK_SIZE bsize = txsize_to_bsize[tx_size];
  const int tx_row = blk_row >> (1 - pd->subsampling_y);
  const int tx_col = blk_col >> (1 - pd->subsampling_x);
  TX_SIZE plane_tx_size;
  const int max_blocks_high = max_block_high(xd, plane_bsize, plane);
  const int max_blocks_wide = max_block_wide(xd, plane_bsize, plane);

  if (blk_row >= max_blocks_high || blk_col >= max_blocks_wide) return;

  plane_tx_size =
      plane ? uv_txsize_lookup[bsize][mbmi->inter_tx_size[tx_row][tx_col]][0][0]
            : mbmi->inter_tx_size[tx_row][tx_col];

  if (tx_size == plane_tx_size) {
    pack_mb_tokens(w, tp, tok_end, bit_depth, tx_size, token_stats);
  } else {
    const int bsl = block_size_wide[bsize] >> (tx_size_wide_log2[0] + 1);
    int i;

    assert(bsl > 0);

    for (i = 0; i < 4; ++i) {
      const int offsetr = blk_row + (i >> 1) * bsl;
      const int offsetc = blk_col + (i & 0x01) * bsl;
      const TX_SIZE sub_txs = tx_size - 1;
      const int step = tx_size_wide_unit[sub_txs] * tx_size_high_unit[sub_txs];

      if (offsetr >= max_blocks_high || offsetc >= max_blocks_wide) continue;

      pack_txb_tokens(w, tp, tok_end, xd, mbmi, plane, plane_bsize, bit_depth,
                      block, offsetr, offsetc, sub_txs, token_stats);
      block += step;
    }
  }
}
#endif

static void write_segment_id(aom_writer *w, const struct segmentation *seg,
                             struct segmentation_probs *segp, int segment_id) {
  if (seg->enabled && seg->update_map) {
#if CONFIG_DAALA_EC
    aom_write_symbol(w, segment_id, segp->tree_cdf, MAX_SEGMENTS);
#else
    aom_write_tree(w, av1_segment_tree, segp->tree_probs, segment_id, 3, 0);
#endif
  }
}

// This function encodes the reference frame
static void write_ref_frames(const AV1_COMMON *cm, const MACROBLOCKD *xd,
                             aom_writer *w) {
  const MB_MODE_INFO *const mbmi = &xd->mi[0]->mbmi;
  const int is_compound = has_second_ref(mbmi);
  const int segment_id = mbmi->segment_id;

  // If segment level coding of this signal is disabled...
  // or the segment allows multiple reference frame options
  if (segfeature_active(&cm->seg, segment_id, SEG_LVL_REF_FRAME)) {
    assert(!is_compound);
    assert(mbmi->ref_frame[0] ==
           get_segdata(&cm->seg, segment_id, SEG_LVL_REF_FRAME));
  } else {
    // does the feature use compound prediction or not
    // (if not specified at the frame/segment level)
    if (cm->reference_mode == REFERENCE_MODE_SELECT) {
      aom_write(w, is_compound, av1_get_reference_mode_prob(cm, xd));
    } else {
      assert((!is_compound) == (cm->reference_mode == SINGLE_REFERENCE));
    }

    if (is_compound) {
#if CONFIG_EXT_REFS
      const int bit = (mbmi->ref_frame[0] == GOLDEN_FRAME ||
                       mbmi->ref_frame[0] == LAST3_FRAME);
      const int bit_bwd = mbmi->ref_frame[1] == ALTREF_FRAME;
#else  // CONFIG_EXT_REFS
      const int bit = mbmi->ref_frame[0] == GOLDEN_FRAME;
#endif  // CONFIG_EXT_REFS

      aom_write(w, bit, av1_get_pred_prob_comp_ref_p(cm, xd));

#if CONFIG_EXT_REFS
      if (!bit) {
        const int bit1 = mbmi->ref_frame[0] == LAST_FRAME;
        aom_write(w, bit1, av1_get_pred_prob_comp_ref_p1(cm, xd));
      } else {
        const int bit2 = mbmi->ref_frame[0] == GOLDEN_FRAME;
        aom_write(w, bit2, av1_get_pred_prob_comp_ref_p2(cm, xd));
      }
      aom_write(w, bit_bwd, av1_get_pred_prob_comp_bwdref_p(cm, xd));
#endif  // CONFIG_EXT_REFS
    } else {
#if CONFIG_EXT_REFS
      const int bit0 = (mbmi->ref_frame[0] == ALTREF_FRAME ||
                        mbmi->ref_frame[0] == BWDREF_FRAME);
      aom_write(w, bit0, av1_get_pred_prob_single_ref_p1(cm, xd));

      if (bit0) {
        const int bit1 = mbmi->ref_frame[0] == ALTREF_FRAME;
        aom_write(w, bit1, av1_get_pred_prob_single_ref_p2(cm, xd));
      } else {
        const int bit2 = (mbmi->ref_frame[0] == LAST3_FRAME ||
                          mbmi->ref_frame[0] == GOLDEN_FRAME);
        aom_write(w, bit2, av1_get_pred_prob_single_ref_p3(cm, xd));

        if (!bit2) {
          const int bit3 = mbmi->ref_frame[0] != LAST_FRAME;
          aom_write(w, bit3, av1_get_pred_prob_single_ref_p4(cm, xd));
        } else {
          const int bit4 = mbmi->ref_frame[0] != LAST3_FRAME;
          aom_write(w, bit4, av1_get_pred_prob_single_ref_p5(cm, xd));
        }
      }
#else   // CONFIG_EXT_REFS
      const int bit0 = mbmi->ref_frame[0] != LAST_FRAME;
      aom_write(w, bit0, av1_get_pred_prob_single_ref_p1(cm, xd));

      if (bit0) {
        const int bit1 = mbmi->ref_frame[0] != GOLDEN_FRAME;
        aom_write(w, bit1, av1_get_pred_prob_single_ref_p2(cm, xd));
      }
#endif  // CONFIG_EXT_REFS
    }
  }
}

#if CONFIG_FILTER_INTRA
static void write_filter_intra_mode_info(const AV1_COMMON *const cm,
                                         const MB_MODE_INFO *const mbmi,
                                         aom_writer *w) {
  if (mbmi->mode == DC_PRED
#if CONFIG_PALETTE
      && mbmi->palette_mode_info.palette_size[0] == 0
#endif  // CONFIG_PALETTE
      ) {
    aom_write(w, mbmi->filter_intra_mode_info.use_filter_intra_mode[0],
              cm->fc->filter_intra_probs[0]);
    if (mbmi->filter_intra_mode_info.use_filter_intra_mode[0]) {
      const FILTER_INTRA_MODE mode =
          mbmi->filter_intra_mode_info.filter_intra_mode[0];
      write_uniform(w, FILTER_INTRA_MODES, mode);
    }
  }

  if (mbmi->uv_mode == DC_PRED
#if CONFIG_PALETTE
      && mbmi->palette_mode_info.palette_size[1] == 0
#endif  // CONFIG_PALETTE
      ) {
    aom_write(w, mbmi->filter_intra_mode_info.use_filter_intra_mode[1],
              cm->fc->filter_intra_probs[1]);
    if (mbmi->filter_intra_mode_info.use_filter_intra_mode[1]) {
      const FILTER_INTRA_MODE mode =
          mbmi->filter_intra_mode_info.filter_intra_mode[1];
      write_uniform(w, FILTER_INTRA_MODES, mode);
    }
  }
}
#endif  // CONFIG_FILTER_INTRA

#if CONFIG_EXT_INTRA
static void write_intra_angle_info(const AV1_COMMON *cm, const MACROBLOCKD *xd,
                                   aom_writer *w) {
  const MB_MODE_INFO *const mbmi = &xd->mi[0]->mbmi;
  const BLOCK_SIZE bsize = mbmi->sb_type;
  const int intra_filter_ctx = av1_get_pred_context_intra_interp(xd);
  int p_angle;

  if (bsize < BLOCK_8X8) return;

  if (mbmi->mode != DC_PRED && mbmi->mode != TM_PRED) {
    write_uniform(w, 2 * MAX_ANGLE_DELTAS + 1,
                  MAX_ANGLE_DELTAS + mbmi->angle_delta[0]);
    p_angle = mode_to_angle_map[mbmi->mode] + mbmi->angle_delta[0] * ANGLE_STEP;
    if (av1_is_intra_filter_switchable(p_angle)) {
      av1_write_token(w, av1_intra_filter_tree,
                      cm->fc->intra_filter_probs[intra_filter_ctx],
                      &intra_filter_encodings[mbmi->intra_filter]);
    }
  }

  if (mbmi->uv_mode != DC_PRED && mbmi->uv_mode != TM_PRED) {
    write_uniform(w, 2 * MAX_ANGLE_DELTAS + 1,
                  MAX_ANGLE_DELTAS + mbmi->angle_delta[1]);
  }
}
#endif  // CONFIG_EXT_INTRA

static void write_switchable_interp_filter(AV1_COMP *cpi, const MACROBLOCKD *xd,
                                           aom_writer *w) {
  AV1_COMMON *const cm = &cpi->common;
  const MB_MODE_INFO *const mbmi = &xd->mi[0]->mbmi;
#if CONFIG_DUAL_FILTER
  int dir;
#endif
  if (cm->interp_filter == SWITCHABLE) {
#if CONFIG_EXT_INTERP
#if CONFIG_DUAL_FILTER
    if (!av1_is_interp_needed(xd)) {
      assert(mbmi->interp_filter[0] == EIGHTTAP_REGULAR);
      return;
    }
#else
    if (!av1_is_interp_needed(xd)) {
#if CONFIG_DUAL_FILTER
      assert(mbmi->interp_filter[0] == EIGHTTAP_REGULAR);
      assert(mbmi->interp_filter[1] == EIGHTTAP_REGULAR);
#else
      assert(mbmi->interp_filter == EIGHTTAP_REGULAR);
#endif
      return;
    }
#endif  // CONFIG_DUAL_FILTER
#endif  // CONFIG_EXT_INTERP
#if CONFIG_DUAL_FILTER
    for (dir = 0; dir < 2; ++dir) {
      if (has_subpel_mv_component(xd->mi[0], xd, dir) ||
          (mbmi->ref_frame[1] > INTRA_FRAME &&
           has_subpel_mv_component(xd->mi[0], xd, dir + 2))) {
        const int ctx = av1_get_pred_context_switchable_interp(xd, dir);
        av1_write_token(w, av1_switchable_interp_tree,
                        cm->fc->switchable_interp_prob[ctx],
                        &switchable_interp_encodings[mbmi->interp_filter[dir]]);
        ++cpi->interp_filter_selected[0][mbmi->interp_filter[dir]];
      }
    }
#else
    {
      const int ctx = av1_get_pred_context_switchable_interp(xd);
#if CONFIG_DAALA_EC
      aom_write_symbol(w, av1_switchable_interp_ind[mbmi->interp_filter],
                       cm->fc->switchable_interp_cdf[ctx], SWITCHABLE_FILTERS);
#else
      av1_write_token(w, av1_switchable_interp_tree,
                      cm->fc->switchable_interp_prob[ctx],
                      &switchable_interp_encodings[mbmi->interp_filter]);
#endif
      ++cpi->interp_filter_selected[0][mbmi->interp_filter];
    }
#endif
  }
}

#if CONFIG_PALETTE
static void write_palette_mode_info(const AV1_COMMON *cm, const MACROBLOCKD *xd,
                                    const MODE_INFO *const mi, aom_writer *w) {
  const MB_MODE_INFO *const mbmi = &mi->mbmi;
  const MODE_INFO *const above_mi = xd->above_mi;
  const MODE_INFO *const left_mi = xd->left_mi;
  const BLOCK_SIZE bsize = mbmi->sb_type;
  const PALETTE_MODE_INFO *const pmi = &mbmi->palette_mode_info;
  int palette_ctx = 0;
  int n, i;

  if (mbmi->mode == DC_PRED) {
    n = pmi->palette_size[0];
    if (above_mi)
      palette_ctx += (above_mi->mbmi.palette_mode_info.palette_size[0] > 0);
    if (left_mi)
      palette_ctx += (left_mi->mbmi.palette_mode_info.palette_size[0] > 0);
    aom_write(w, n > 0,
              av1_default_palette_y_mode_prob[bsize - BLOCK_8X8][palette_ctx]);
    if (n > 0) {
      av1_write_token(w, av1_palette_size_tree,
                      av1_default_palette_y_size_prob[bsize - BLOCK_8X8],
                      &palette_size_encodings[n - 2]);
      for (i = 0; i < n; ++i)
        aom_write_literal(w, pmi->palette_colors[i], cm->bit_depth);
      write_uniform(w, n, pmi->palette_first_color_idx[0]);
    }
  }

  if (mbmi->uv_mode == DC_PRED) {
    n = pmi->palette_size[1];
    aom_write(w, n > 0,
              av1_default_palette_uv_mode_prob[pmi->palette_size[0] > 0]);
    if (n > 0) {
      av1_write_token(w, av1_palette_size_tree,
                      av1_default_palette_uv_size_prob[bsize - BLOCK_8X8],
                      &palette_size_encodings[n - 2]);
      for (i = 0; i < n; ++i) {
        aom_write_literal(w, pmi->palette_colors[PALETTE_MAX_SIZE + i],
                          cm->bit_depth);
        aom_write_literal(w, pmi->palette_colors[2 * PALETTE_MAX_SIZE + i],
                          cm->bit_depth);
      }
      write_uniform(w, n, pmi->palette_first_color_idx[1]);
    }
  }
}
#endif  // CONFIG_PALETTE

static void write_tx_type(const AV1_COMMON *const cm,
                          const MB_MODE_INFO *const mbmi,
#if CONFIG_SUPERTX
                          const int supertx_enabled,
#endif
                          aom_writer *w) {
  const int is_inter = is_inter_block(mbmi);
#if CONFIG_VAR_TX
  const TX_SIZE tx_size = is_inter ? mbmi->min_tx_size : mbmi->tx_size;
#else
  const TX_SIZE tx_size = mbmi->tx_size;
#endif
  if (!FIXED_TX_TYPE) {
#if CONFIG_EXT_TX
    const BLOCK_SIZE bsize = mbmi->sb_type;
    if (get_ext_tx_types(tx_size, bsize, is_inter) > 1 && cm->base_qindex > 0 &&
        !mbmi->skip &&
#if CONFIG_SUPERTX
        !supertx_enabled &&
#endif  // CONFIG_SUPERTX
        !segfeature_active(&cm->seg, mbmi->segment_id, SEG_LVL_SKIP)) {
      int eset = get_ext_tx_set(tx_size, bsize, is_inter);
      if (is_inter) {
        assert(ext_tx_used_inter[eset][mbmi->tx_type]);
        if (eset > 0)
          av1_write_token(
              w, av1_ext_tx_inter_tree[eset],
              cm->fc->inter_ext_tx_prob[eset][txsize_sqr_map[tx_size]],
              &ext_tx_inter_encodings[eset][mbmi->tx_type]);
      } else if (ALLOW_INTRA_EXT_TX) {
        if (eset > 0)
          av1_write_token(w, av1_ext_tx_intra_tree[eset],
                          cm->fc->intra_ext_tx_prob[eset][tx_size][mbmi->mode],
                          &ext_tx_intra_encodings[eset][mbmi->tx_type]);
      }
    }
#else
    if (tx_size < TX_32X32 && cm->base_qindex > 0 && !mbmi->skip &&
#if CONFIG_SUPERTX
        !supertx_enabled &&
#endif  // CONFIG_SUPERTX
        !segfeature_active(&cm->seg, mbmi->segment_id, SEG_LVL_SKIP)) {
      if (is_inter) {
#if CONFIG_DAALA_EC
        aom_write_symbol(w, av1_ext_tx_ind[mbmi->tx_type],
                         cm->fc->inter_ext_tx_cdf[tx_size], TX_TYPES);
#else
        av1_write_token(w, av1_ext_tx_tree, cm->fc->inter_ext_tx_prob[tx_size],
                        &ext_tx_encodings[mbmi->tx_type]);
#endif
      } else {
#if CONFIG_DAALA_EC
        aom_write_symbol(
            w, av1_ext_tx_ind[mbmi->tx_type],
            cm->fc->intra_ext_tx_cdf[tx_size]
                                    [intra_mode_to_tx_type_context[mbmi->mode]],
            TX_TYPES);
#else
        av1_write_token(
            w, av1_ext_tx_tree,
            cm->fc
                ->intra_ext_tx_prob[tx_size]
                                   [intra_mode_to_tx_type_context[mbmi->mode]],
            &ext_tx_encodings[mbmi->tx_type]);
#endif
      }
    }
#endif  // CONFIG_EXT_TX
  }
}

static void pack_inter_mode_mvs(AV1_COMP *cpi, const MODE_INFO *mi,
#if CONFIG_SUPERTX
                                int supertx_enabled,
#endif
                                aom_writer *w) {
  AV1_COMMON *const cm = &cpi->common;
#if !CONFIG_REF_MV
  nmv_context *nmvc = &cm->fc->nmvc;
#endif

#if CONFIG_DELTA_Q
  MACROBLOCK *const x = &cpi->td.mb;
  MACROBLOCKD *const xd = &x->e_mbd;
#else
  const MACROBLOCK *x = &cpi->td.mb;
  const MACROBLOCKD *xd = &x->e_mbd;
#endif
  const struct segmentation *const seg = &cm->seg;
  struct segmentation_probs *const segp = &cm->fc->seg;
  const MB_MODE_INFO *const mbmi = &mi->mbmi;
  const MB_MODE_INFO_EXT *const mbmi_ext = x->mbmi_ext;
  const PREDICTION_MODE mode = mbmi->mode;
  const int segment_id = mbmi->segment_id;
  const BLOCK_SIZE bsize = mbmi->sb_type;
  const int allow_hp = cm->allow_high_precision_mv;
  const int is_inter = is_inter_block(mbmi);
  const int is_compound = has_second_ref(mbmi);
  int skip, ref;

  if (seg->update_map) {
    if (seg->temporal_update) {
      const int pred_flag = mbmi->seg_id_predicted;
      aom_prob pred_prob = av1_get_pred_prob_seg_id(segp, xd);
      aom_write(w, pred_flag, pred_prob);
      if (!pred_flag) write_segment_id(w, seg, segp, segment_id);
    } else {
      write_segment_id(w, seg, segp, segment_id);
    }
  }

#if CONFIG_SUPERTX
  if (supertx_enabled)
    skip = mbmi->skip;
  else
    skip = write_skip(cm, xd, segment_id, mi, w);
#else
  skip = write_skip(cm, xd, segment_id, mi, w);
#endif  // CONFIG_SUPERTX
#if CONFIG_DELTA_Q
  if (cm->delta_q_present_flag) {
    int mi_row = (-xd->mb_to_top_edge) >> (MI_SIZE_LOG2 + 3);
    int mi_col = (-xd->mb_to_left_edge) >> (MI_SIZE_LOG2 + 3);
    int super_block_upper_left =
        ((mi_row & MAX_MIB_MASK) == 0) && ((mi_col & MAX_MIB_MASK) == 0);
    if ((bsize != BLOCK_64X64 || skip == 0) && super_block_upper_left) {
      int reduced_delta_qindex =
          (mbmi->current_q_index - xd->prev_qindex) / cm->delta_q_res;
      write_delta_qindex(cm, reduced_delta_qindex, w);
      xd->prev_qindex = mbmi->current_q_index;
    }
  }
#endif

#if CONFIG_SUPERTX
  if (!supertx_enabled)
#endif  // CONFIG_SUPERTX
    if (!segfeature_active(seg, segment_id, SEG_LVL_REF_FRAME))
      aom_write(w, is_inter, av1_get_intra_inter_prob(cm, xd));

  if (bsize >= BLOCK_8X8 && cm->tx_mode == TX_MODE_SELECT &&
#if CONFIG_SUPERTX
      !supertx_enabled &&
#endif  // CONFIG_SUPERTX
      !(is_inter && skip) && !xd->lossless[segment_id]) {
#if CONFIG_VAR_TX
    if (is_inter) {  // This implies skip flag is 0.
      const TX_SIZE max_tx_size = max_txsize_lookup[bsize];
      const int bh = tx_size_high_unit[max_tx_size];
      const int bw = tx_size_wide_unit[max_tx_size];
      const int width = num_4x4_blocks_wide_lookup[bsize];
      const int height = num_4x4_blocks_high_lookup[bsize];
      int idx, idy;

#if CONFIG_EXT_TX && CONFIG_RECT_TX
      if (is_rect_tx_allowed(xd, mbmi)) {
        int tx_size_cat = inter_tx_size_cat_lookup[bsize];

        aom_write(w, is_rect_tx(mbmi->tx_size),
                  cm->fc->rect_tx_prob[tx_size_cat]);
      }

      if (is_rect_tx(mbmi->tx_size)) {
        set_txfm_ctxs(mbmi->tx_size, xd->n8_w, xd->n8_h, xd);
      } else {
#endif  // CONFIG_EXT_TX && CONFIG_RECT_TX
        for (idy = 0; idy < height; idy += bh)
          for (idx = 0; idx < width; idx += bw)
            write_tx_size_vartx(cm, xd, mbmi, max_tx_size, height != width, idy,
                                idx, w);
#if CONFIG_EXT_TX && CONFIG_RECT_TX
      }
#endif  // CONFIG_EXT_TX && CONFIG_RECT_TX
    } else {
      set_txfm_ctxs(mbmi->tx_size, xd->n8_w, xd->n8_h, xd);
      write_selected_tx_size(cm, xd, w);
    }
  } else {
    set_txfm_ctxs(mbmi->tx_size, xd->n8_w, xd->n8_h, xd);
#else
    write_selected_tx_size(cm, xd, w);
#endif
  }

  if (!is_inter) {
    if (bsize >= BLOCK_8X8) {
#if CONFIG_DAALA_EC
      aom_write_symbol(w, av1_intra_mode_ind[mode],
                       cm->fc->y_mode_cdf[size_group_lookup[bsize]],
                       INTRA_MODES);
#else
      write_intra_mode(w, mode, cm->fc->y_mode_prob[size_group_lookup[bsize]]);
#endif
    } else {
      int idx, idy;
      const int num_4x4_w = num_4x4_blocks_wide_lookup[bsize];
      const int num_4x4_h = num_4x4_blocks_high_lookup[bsize];
      for (idy = 0; idy < 2; idy += num_4x4_h) {
        for (idx = 0; idx < 2; idx += num_4x4_w) {
          const PREDICTION_MODE b_mode = mi->bmi[idy * 2 + idx].as_mode;
#if CONFIG_DAALA_EC
          aom_write_symbol(w, av1_intra_mode_ind[b_mode], cm->fc->y_mode_cdf[0],
                           INTRA_MODES);
#else
          write_intra_mode(w, b_mode, cm->fc->y_mode_prob[0]);
#endif
        }
      }
    }
#if CONFIG_DAALA_EC
    aom_write_symbol(w, av1_intra_mode_ind[mbmi->uv_mode],
                     cm->fc->uv_mode_cdf[mode], INTRA_MODES);
#else
    write_intra_mode(w, mbmi->uv_mode, cm->fc->uv_mode_prob[mode]);
#endif
#if CONFIG_EXT_INTRA
    write_intra_angle_info(cm, xd, w);
#endif  // CONFIG_EXT_INTRA
#if CONFIG_PALETTE
    if (bsize >= BLOCK_8X8 && cm->allow_screen_content_tools)
      write_palette_mode_info(cm, xd, mi, w);
#endif  // CONFIG_PALETTE
#if CONFIG_FILTER_INTRA
    if (bsize >= BLOCK_8X8) write_filter_intra_mode_info(cm, mbmi, w);
#endif  // CONFIG_FILTER_INTRA
  } else {
    int16_t mode_ctx = mbmi_ext->mode_context[mbmi->ref_frame[0]];
    write_ref_frames(cm, xd, w);

#if CONFIG_REF_MV
#if CONFIG_EXT_INTER
    if (is_compound)
      mode_ctx = mbmi_ext->compound_mode_context[mbmi->ref_frame[0]];
    else
#endif  // CONFIG_EXT_INTER
      mode_ctx = av1_mode_context_analyzer(mbmi_ext->mode_context,
                                           mbmi->ref_frame, bsize, -1);
#endif

    // If segment skip is not enabled code the mode.
    if (!segfeature_active(seg, segment_id, SEG_LVL_SKIP)) {
      if (bsize >= BLOCK_8X8) {
#if CONFIG_EXT_INTER
        if (is_inter_compound_mode(mode))
          write_inter_compound_mode(cm, w, mode, mode_ctx);
        else if (is_inter_singleref_mode(mode))
#endif  // CONFIG_EXT_INTER
          write_inter_mode(cm, w, mode,
#if CONFIG_REF_MV && CONFIG_EXT_INTER
                           is_compound,
#endif  // CONFIG_REF_MV && CONFIG_EXT_INTER
                           mode_ctx);

#if CONFIG_REF_MV
        if (mode == NEARMV || mode == NEWMV)
          write_drl_idx(cm, mbmi, mbmi_ext, w);
#endif
      }
    }

#if !CONFIG_EXT_INTERP && !CONFIG_DUAL_FILTER
    write_switchable_interp_filter(cpi, xd, w);
#endif  // !CONFIG_EXT_INTERP

    if (bsize < BLOCK_8X8) {
      const int num_4x4_w = num_4x4_blocks_wide_lookup[bsize];
      const int num_4x4_h = num_4x4_blocks_high_lookup[bsize];
      int idx, idy;
      for (idy = 0; idy < 2; idy += num_4x4_h) {
        for (idx = 0; idx < 2; idx += num_4x4_w) {
          const int j = idy * 2 + idx;
          const PREDICTION_MODE b_mode = mi->bmi[j].as_mode;
#if CONFIG_REF_MV
#if CONFIG_EXT_INTER
          if (!is_compound)
#endif  // CONFIG_EXT_INTER
            mode_ctx = av1_mode_context_analyzer(mbmi_ext->mode_context,
                                                 mbmi->ref_frame, bsize, j);
#endif
#if CONFIG_EXT_INTER
          if (is_inter_compound_mode(b_mode))
            write_inter_compound_mode(cm, w, b_mode, mode_ctx);
          else if (is_inter_singleref_mode(b_mode))
#endif  // CONFIG_EXT_INTER
            write_inter_mode(cm, w, b_mode,
#if CONFIG_REF_MV && CONFIG_EXT_INTER
                             has_second_ref(mbmi),
#endif  // CONFIG_REF_MV && CONFIG_EXT_INTER
                             mode_ctx);

#if CONFIG_EXT_INTER
          if (b_mode == NEWMV || b_mode == NEWFROMNEARMV ||
              b_mode == NEW_NEWMV) {
#else
          if (b_mode == NEWMV) {
#endif  // CONFIG_EXT_INTER
            for (ref = 0; ref < 1 + is_compound; ++ref) {
#if CONFIG_REF_MV
              int8_t rf_type = av1_ref_frame_type(mbmi->ref_frame);
              int nmv_ctx = av1_nmv_ctx(mbmi_ext->ref_mv_count[rf_type],
                                        mbmi_ext->ref_mv_stack[rf_type], ref,
                                        mbmi->ref_mv_idx);
              nmv_context *nmvc = &cm->fc->nmvc[nmv_ctx];
#endif
              av1_encode_mv(cpi, w, &mi->bmi[j].as_mv[ref].as_mv,
#if CONFIG_EXT_INTER
                            &mi->bmi[j].ref_mv[ref].as_mv,
#if CONFIG_REF_MV
                            is_compound,
#endif
#else
#if CONFIG_REF_MV
                            &mi->bmi[j].pred_mv[ref].as_mv, is_compound,
#else
                            &mbmi_ext->ref_mvs[mbmi->ref_frame[ref]][0].as_mv,
#endif  // CONFIG_REF_MV
#endif  // CONFIG_EXT_INTER
                            nmvc, allow_hp);
            }
          }
#if CONFIG_EXT_INTER
          else if (b_mode == NEAREST_NEWMV || b_mode == NEAR_NEWMV) {
#if CONFIG_REF_MV
            int8_t rf_type = av1_ref_frame_type(mbmi->ref_frame);
            int nmv_ctx = av1_nmv_ctx(mbmi_ext->ref_mv_count[rf_type],
                                      mbmi_ext->ref_mv_stack[rf_type], 1,
                                      mbmi->ref_mv_idx);
            nmv_context *nmvc = &cm->fc->nmvc[nmv_ctx];
#endif
            av1_encode_mv(cpi, w, &mi->bmi[j].as_mv[1].as_mv,
                          &mi->bmi[j].ref_mv[1].as_mv,
#if CONFIG_REF_MV
                          is_compound,
#endif
                          nmvc, allow_hp);
          } else if (b_mode == NEW_NEARESTMV || b_mode == NEW_NEARMV) {
#if CONFIG_REF_MV
            int8_t rf_type = av1_ref_frame_type(mbmi->ref_frame);
            int nmv_ctx = av1_nmv_ctx(mbmi_ext->ref_mv_count[rf_type],
                                      mbmi_ext->ref_mv_stack[rf_type], 0,
                                      mbmi->ref_mv_idx);
            nmv_context *nmvc = &cm->fc->nmvc[nmv_ctx];
#endif
            av1_encode_mv(cpi, w, &mi->bmi[j].as_mv[0].as_mv,
                          &mi->bmi[j].ref_mv[0].as_mv,
#if CONFIG_REF_MV
                          is_compound,
#endif
                          nmvc, allow_hp);
          }
#endif  // CONFIG_EXT_INTER
        }
      }
    } else {
#if CONFIG_EXT_INTER
      if (mode == NEWMV || mode == NEWFROMNEARMV || mode == NEW_NEWMV) {
#else
      if (mode == NEWMV) {
#endif  // CONFIG_EXT_INTER
        int_mv ref_mv;
        for (ref = 0; ref < 1 + is_compound; ++ref) {
#if CONFIG_REF_MV
          int8_t rf_type = av1_ref_frame_type(mbmi->ref_frame);
          int nmv_ctx = av1_nmv_ctx(mbmi_ext->ref_mv_count[rf_type],
                                    mbmi_ext->ref_mv_stack[rf_type], ref,
                                    mbmi->ref_mv_idx);
          nmv_context *nmvc = &cm->fc->nmvc[nmv_ctx];
#endif
          ref_mv = mbmi_ext->ref_mvs[mbmi->ref_frame[ref]][0];
#if CONFIG_EXT_INTER
          if (mode == NEWFROMNEARMV)
            av1_encode_mv(cpi, w, &mbmi->mv[ref].as_mv,
                          &mbmi_ext->ref_mvs[mbmi->ref_frame[ref]][1].as_mv,
#if CONFIG_REF_MV
                          is_compound,
#endif
                          nmvc, allow_hp);
          else
#endif  // CONFIG_EXT_INTER
            av1_encode_mv(cpi, w, &mbmi->mv[ref].as_mv, &ref_mv.as_mv,
#if CONFIG_REF_MV
                          is_compound,
#endif
                          nmvc, allow_hp);
        }
#if CONFIG_EXT_INTER
      } else if (mode == NEAREST_NEWMV || mode == NEAR_NEWMV) {
#if CONFIG_REF_MV
        int8_t rf_type = av1_ref_frame_type(mbmi->ref_frame);
        int nmv_ctx =
            av1_nmv_ctx(mbmi_ext->ref_mv_count[rf_type],
                        mbmi_ext->ref_mv_stack[rf_type], 1, mbmi->ref_mv_idx);
        nmv_context *nmvc = &cm->fc->nmvc[nmv_ctx];
#endif
        av1_encode_mv(cpi, w, &mbmi->mv[1].as_mv,
                      &mbmi_ext->ref_mvs[mbmi->ref_frame[1]][0].as_mv,
#if CONFIG_REF_MV
                      is_compound,
#endif
                      nmvc, allow_hp);
      } else if (mode == NEW_NEARESTMV || mode == NEW_NEARMV) {
#if CONFIG_REF_MV
        int8_t rf_type = av1_ref_frame_type(mbmi->ref_frame);
        int nmv_ctx =
            av1_nmv_ctx(mbmi_ext->ref_mv_count[rf_type],
                        mbmi_ext->ref_mv_stack[rf_type], 0, mbmi->ref_mv_idx);
        nmv_context *nmvc = &cm->fc->nmvc[nmv_ctx];
#endif
        av1_encode_mv(cpi, w, &mbmi->mv[0].as_mv,
                      &mbmi_ext->ref_mvs[mbmi->ref_frame[0]][0].as_mv,
#if CONFIG_REF_MV
                      is_compound,
#endif
                      nmvc, allow_hp);
#endif  // CONFIG_EXT_INTER
      }
    }

#if CONFIG_EXT_INTER
    if (cpi->common.reference_mode != COMPOUND_REFERENCE &&
#if CONFIG_SUPERTX
        !supertx_enabled &&
#endif  // CONFIG_SUPERTX
        is_interintra_allowed(mbmi)) {
      const int interintra = mbmi->ref_frame[1] == INTRA_FRAME;
      const int bsize_group = size_group_lookup[bsize];
      aom_write(w, interintra, cm->fc->interintra_prob[bsize_group]);
      if (interintra) {
        write_interintra_mode(w, mbmi->interintra_mode,
                              cm->fc->interintra_mode_prob[bsize_group]);
        if (is_interintra_wedge_used(bsize)) {
          aom_write(w, mbmi->use_wedge_interintra,
                    cm->fc->wedge_interintra_prob[bsize]);
          if (mbmi->use_wedge_interintra) {
            aom_write_literal(w, mbmi->interintra_wedge_index,
                              get_wedge_bits_lookup(bsize));
            assert(mbmi->interintra_wedge_sign == 0);
          }
        }
      }
    }
#endif  // CONFIG_EXT_INTER

#if CONFIG_MOTION_VAR || CONFIG_WARPED_MOTION
#if CONFIG_SUPERTX
    if (!supertx_enabled)
#endif  // CONFIG_SUPERTX
#if CONFIG_EXT_INTER
      if (mbmi->ref_frame[1] != INTRA_FRAME)
#endif  // CONFIG_EXT_INTER
        if (is_motion_variation_allowed(mbmi)) {
          // TODO(debargha): Might want to only emit this if SEG_LVL_SKIP
          // is not active, and assume SIMPLE_TRANSLATION in the decoder if
          // it is active.
          assert(mbmi->motion_mode < MOTION_MODES);
          av1_write_token(w, av1_motion_mode_tree,
                          cm->fc->motion_mode_prob[bsize],
                          &motion_mode_encodings[mbmi->motion_mode]);
        }
#endif  // CONFIG_MOTION_VAR || CONFIG_WARPED_MOTION

#if CONFIG_EXT_INTER
    if (cpi->common.reference_mode != SINGLE_REFERENCE &&
        is_inter_compound_mode(mbmi->mode) &&
#if CONFIG_MOTION_VAR
        !(is_motion_variation_allowed(mbmi) &&
          mbmi->motion_mode != SIMPLE_TRANSLATION) &&
#endif  // CONFIG_MOTION_VAR
        is_interinter_wedge_used(bsize)) {
      aom_write(w, mbmi->use_wedge_interinter,
                cm->fc->wedge_interinter_prob[bsize]);
      if (mbmi->use_wedge_interinter) {
        aom_write_literal(w, mbmi->interinter_wedge_index,
                          get_wedge_bits_lookup(bsize));
        aom_write_bit(w, mbmi->interinter_wedge_sign);
      }
    }
#endif  // CONFIG_EXT_INTER

#if CONFIG_EXT_INTERP || CONFIG_DUAL_FILTER
    write_switchable_interp_filter(cpi, xd, w);
#endif  // CONFIG_EXT_INTERP
  }

  write_tx_type(cm, mbmi,
#if CONFIG_SUPERTX
                supertx_enabled,
#endif
                w);
}

#if CONFIG_DELTA_Q
static void write_mb_modes_kf(AV1_COMMON *cm, MACROBLOCKD *xd,
                              MODE_INFO **mi_8x8, aom_writer *w) {
  int skip;
#else
static void write_mb_modes_kf(AV1_COMMON *cm, const MACROBLOCKD *xd,
                              MODE_INFO **mi_8x8, aom_writer *w) {
#endif
  const struct segmentation *const seg = &cm->seg;
  struct segmentation_probs *const segp = &cm->fc->seg;
  const MODE_INFO *const mi = mi_8x8[0];
  const MODE_INFO *const above_mi = xd->above_mi;
  const MODE_INFO *const left_mi = xd->left_mi;
  const MB_MODE_INFO *const mbmi = &mi->mbmi;
  const BLOCK_SIZE bsize = mbmi->sb_type;

  if (seg->update_map) write_segment_id(w, seg, segp, mbmi->segment_id);

#if CONFIG_DELTA_Q
  skip = write_skip(cm, xd, mbmi->segment_id, mi, w);
  if (cm->delta_q_present_flag) {
    int mi_row = (-xd->mb_to_top_edge) >> 6;
    int mi_col = (-xd->mb_to_left_edge) >> 6;
    int super_block_upper_left = ((mi_row & 7) == 0) && ((mi_col & 7) == 0);
    if ((bsize != BLOCK_64X64 || skip == 0) && super_block_upper_left) {
      int reduced_delta_qindex =
          (mbmi->current_q_index - xd->prev_qindex) / cm->delta_q_res;
      write_delta_qindex(cm, reduced_delta_qindex, w);
      xd->prev_qindex = mbmi->current_q_index;
    }
  }
#else
  write_skip(cm, xd, mbmi->segment_id, mi, w);
#endif

  if (bsize >= BLOCK_8X8 && cm->tx_mode == TX_MODE_SELECT &&
      !xd->lossless[mbmi->segment_id])
    write_selected_tx_size(cm, xd, w);

  if (bsize >= BLOCK_8X8) {
#if CONFIG_DAALA_EC
    aom_write_symbol(w, av1_intra_mode_ind[mbmi->mode],
                     get_y_mode_cdf(cm, mi, above_mi, left_mi, 0), INTRA_MODES);
#else
    write_intra_mode(w, mbmi->mode,
                     get_y_mode_probs(cm, mi, above_mi, left_mi, 0));
#endif
  } else {
    const int num_4x4_w = num_4x4_blocks_wide_lookup[bsize];
    const int num_4x4_h = num_4x4_blocks_high_lookup[bsize];
    int idx, idy;

    for (idy = 0; idy < 2; idy += num_4x4_h) {
      for (idx = 0; idx < 2; idx += num_4x4_w) {
        const int block = idy * 2 + idx;
#if CONFIG_DAALA_EC
        aom_write_symbol(w, av1_intra_mode_ind[mi->bmi[block].as_mode],
                         get_y_mode_cdf(cm, mi, above_mi, left_mi, block),
                         INTRA_MODES);
#else
        write_intra_mode(w, mi->bmi[block].as_mode,
                         get_y_mode_probs(cm, mi, above_mi, left_mi, block));
#endif
      }
    }
  }
#if CONFIG_DAALA_EC
  aom_write_symbol(w, av1_intra_mode_ind[mbmi->uv_mode],
                   cm->fc->uv_mode_cdf[mbmi->mode], INTRA_MODES);
#else
  write_intra_mode(w, mbmi->uv_mode, cm->fc->uv_mode_prob[mbmi->mode]);
#endif
#if CONFIG_EXT_INTRA
  write_intra_angle_info(cm, xd, w);
#endif  // CONFIG_EXT_INTRA
#if CONFIG_PALETTE
  if (bsize >= BLOCK_8X8 && cm->allow_screen_content_tools)
    write_palette_mode_info(cm, xd, mi, w);
#endif  // CONFIG_PALETTE
#if CONFIG_FILTER_INTRA
  if (bsize >= BLOCK_8X8) write_filter_intra_mode_info(cm, mbmi, w);
#endif  // CONFIG_FILTER_INTRA

  write_tx_type(cm, mbmi,
#if CONFIG_SUPERTX
                0,
#endif
                w);
}

#if CONFIG_SUPERTX
#define write_modes_b_wrapper(cpi, tile, w, tok, tok_end, supertx_enabled, \
                              mi_row, mi_col)                              \
  write_modes_b(cpi, tile, w, tok, tok_end, supertx_enabled, mi_row, mi_col)
#else
#define write_modes_b_wrapper(cpi, tile, w, tok, tok_end, supertx_enabled, \
                              mi_row, mi_col)                              \
  write_modes_b(cpi, tile, w, tok, tok_end, mi_row, mi_col)
#endif  // CONFIG_SUPERTX

#if CONFIG_RD_DEBUG
static void dump_mode_info(MODE_INFO *mi) {
  printf("\nmi->mbmi.mi_row == %d\n", mi->mbmi.mi_row);
  printf("&& mi->mbmi.mi_col == %d\n", mi->mbmi.mi_col);
  printf("&& mi->mbmi.sb_type == %d\n", mi->mbmi.sb_type);
  printf("&& mi->mbmi.tx_size == %d\n", mi->mbmi.tx_size);
  if (mi->mbmi.sb_type >= BLOCK_8X8) {
    printf("&& mi->mbmi.mode == %d\n", mi->mbmi.mode);
  } else {
    printf("&& mi->bmi[0].as_mode == %d\n", mi->bmi[0].as_mode);
  }
}
#endif

#if CONFIG_PVQ
PVQ_INFO *get_pvq_block(PVQ_QUEUE *pvq_q) {
  PVQ_INFO *pvq;

  assert(pvq_q->curr_pos <= pvq_q->last_pos);
  assert(pvq_q->curr_pos < pvq_q->buf_len);

  pvq = pvq_q->buf + pvq_q->curr_pos;
  ++pvq_q->curr_pos;

  return pvq;
}
#endif

static void write_modes_b(AV1_COMP *cpi, const TileInfo *const tile,
                          aom_writer *w, const TOKENEXTRA **tok,
                          const TOKENEXTRA *const tok_end,
#if CONFIG_SUPERTX
                          int supertx_enabled,
#endif
                          int mi_row, int mi_col) {
  AV1_COMMON *const cm = &cpi->common;
  MACROBLOCKD *const xd = &cpi->td.mb.e_mbd;
  MODE_INFO *m;
  int plane;
  int bh, bw;
#if CONFIG_RD_DEBUG
  int64_t txb_coeff_cost[MAX_MB_PLANE] = { 0 };
#endif
#if CONFIG_RANS
  (void)tok;
  (void)tok_end;
  (void)plane;
#endif  // !CONFIG_RANS
#if CONFIG_PVQ
  MB_MODE_INFO *mbmi;
  BLOCK_SIZE bsize;
  od_adapt_ctx *adapt;
  (void)tok;
  (void)tok_end;
#endif
  xd->mi = cm->mi_grid_visible + (mi_row * cm->mi_stride + mi_col);
  m = xd->mi[0];

  assert(m->mbmi.sb_type <= cm->sb_size);

  bh = num_8x8_blocks_high_lookup[m->mbmi.sb_type];
  bw = num_8x8_blocks_wide_lookup[m->mbmi.sb_type];

  cpi->td.mb.mbmi_ext = cpi->mbmi_ext_base + (mi_row * cm->mi_cols + mi_col);

  set_mi_row_col(xd, tile, mi_row, bh, mi_col, bw, cm->mi_rows, cm->mi_cols);
#if CONFIG_PVQ
  mbmi = &m->mbmi;
  bsize = mbmi->sb_type;
  adapt = &cpi->td.mb.daala_enc.state.adapt;
#endif

  if (frame_is_intra_only(cm)) {
    write_mb_modes_kf(cm, xd, xd->mi, w);
  } else {
#if CONFIG_VAR_TX
    xd->above_txfm_context = cm->above_txfm_context + mi_col;
    xd->left_txfm_context =
        xd->left_txfm_context_buffer + (mi_row & MAX_MIB_MASK);
#endif
#if CONFIG_EXT_INTERP
    // av1_is_interp_needed needs the ref frame buffers set up to look
    // up if they are scaled. av1_is_interp_needed is in turn needed by
    // write_switchable_interp_filter, which is called by pack_inter_mode_mvs.
    set_ref_ptrs(cm, xd, m->mbmi.ref_frame[0], m->mbmi.ref_frame[1]);
#endif  // CONFIG_EXT_INTERP
#if 0
    // NOTE(zoeliu): For debug
    if (cm->current_video_frame == FRAME_TO_CHECK && cm->show_frame == 1) {
      const PREDICTION_MODE mode = m->mbmi.mode;
      const int segment_id = m->mbmi.segment_id;
      const BLOCK_SIZE bsize = m->mbmi.sb_type;

      // For sub8x8, simply dump out the first sub8x8 block info
      const PREDICTION_MODE b_mode =
          (bsize < BLOCK_8X8) ? m->bmi[0].as_mode : -1;
      const int mv_x = (bsize < BLOCK_8X8) ?
          m->bmi[0].as_mv[0].as_mv.row : m->mbmi.mv[0].as_mv.row;
      const int mv_y = (bsize < BLOCK_8X8) ?
          m->bmi[0].as_mv[0].as_mv.col : m->mbmi.mv[0].as_mv.col;

      printf("Before pack_inter_mode_mvs(): "
             "Frame=%d, (mi_row,mi_col)=(%d,%d), "
             "mode=%d, segment_id=%d, bsize=%d, b_mode=%d, "
             "mv[0]=(%d, %d), ref[0]=%d, ref[1]=%d\n",
             cm->current_video_frame, mi_row, mi_col,
             mode, segment_id, bsize, b_mode, mv_x, mv_y,
             m->mbmi.ref_frame[0], m->mbmi.ref_frame[1]);
    }
#endif  // 0
    pack_inter_mode_mvs(cpi, m,
#if CONFIG_SUPERTX
                        supertx_enabled,
#endif
                        w);
  }

#if CONFIG_PALETTE
  for (plane = 0; plane <= 1; ++plane) {
    if (m->mbmi.palette_mode_info.palette_size[plane] > 0) {
      const int rows = (4 * num_4x4_blocks_high_lookup[m->mbmi.sb_type]) >>
                       (xd->plane[plane].subsampling_y);
      const int cols = (4 * num_4x4_blocks_wide_lookup[m->mbmi.sb_type]) >>
                       (xd->plane[plane].subsampling_x);
      assert(*tok < tok_end);
      pack_palette_tokens(w, tok, m->mbmi.palette_mode_info.palette_size[plane],
                          rows * cols - 1);
      assert(*tok < tok_end + m->mbmi.skip);
    }
  }
#endif  // CONFIG_PALETTE
#if !CONFIG_PVQ
#if CONFIG_SUPERTX
  if (supertx_enabled) return;
#endif  // CONFIG_SUPERTX

  if (!m->mbmi.skip) {
    assert(*tok < tok_end);
    for (plane = 0; plane < MAX_MB_PLANE; ++plane) {
#if CONFIG_VAR_TX
      const struct macroblockd_plane *const pd = &xd->plane[plane];
      MB_MODE_INFO *mbmi = &m->mbmi;
      BLOCK_SIZE bsize = mbmi->sb_type;
      const BLOCK_SIZE plane_bsize =
          get_plane_block_size(AOMMAX(bsize, BLOCK_8X8), pd);

      const int num_4x4_w =
          block_size_wide[plane_bsize] >> tx_size_wide_log2[0];
      const int num_4x4_h =
          block_size_high[plane_bsize] >> tx_size_wide_log2[0];
      int row, col;
#if CONFIG_EXT_TX && CONFIG_RECT_TX
      TX_SIZE tx_size =
          plane ? get_uv_tx_size(mbmi, &xd->plane[plane]) : mbmi->tx_size;
#endif

      TOKEN_STATS token_stats;
      token_stats.cost = 0;

#if CONFIG_EXT_TX && CONFIG_RECT_TX

      if (is_inter_block(mbmi) && !is_rect_tx(tx_size))
#else
      if (is_inter_block(mbmi))
#endif
      {
        const TX_SIZE max_tx_size = max_txsize_lookup[plane_bsize];
        int block = 0;
        const int step =
            tx_size_wide_unit[max_tx_size] * tx_size_high_unit[max_tx_size];
        const int bkw = tx_size_wide_unit[max_tx_size];
        const int bkh = tx_size_high_unit[max_tx_size];
        for (row = 0; row < num_4x4_h; row += bkh) {
          for (col = 0; col < num_4x4_w; col += bkw) {
            pack_txb_tokens(w, tok, tok_end, xd, mbmi, plane, plane_bsize,
                            cm->bit_depth, block, row, col, max_tx_size,
                            &token_stats);
            block += step;
          }
        }
      } else {
        TX_SIZE tx = plane ? get_uv_tx_size(&m->mbmi, &xd->plane[plane])
                           : m->mbmi.tx_size;
        const int bkw = tx_size_wide_unit[tx];
        const int bkh = tx_size_high_unit[tx];

        for (row = 0; row < num_4x4_h; row += bkh)
          for (col = 0; col < num_4x4_w; col += bkw)
            pack_mb_tokens(w, tok, tok_end, cm->bit_depth, tx, &token_stats);
      }
#else
      TX_SIZE tx =
          plane ? get_uv_tx_size(&m->mbmi, &xd->plane[plane]) : m->mbmi.tx_size;
      TOKEN_STATS token_stats;
      token_stats.cost = 0;
      pack_mb_tokens(w, tok, tok_end, cm->bit_depth, tx, &token_stats);
#endif  // CONFIG_VAR_TX

#if CONFIG_RD_DEBUG
      txb_coeff_cost[plane] += token_stats.cost;
#else
      (void)token_stats;
#endif

      assert(*tok < tok_end && (*tok)->token == EOSB_TOKEN);
      (*tok)++;
    }
  }

#if CONFIG_RD_DEBUG
  for (plane = 0; plane < MAX_MB_PLANE; ++plane) {
    if (m->mbmi.rd_stats.txb_coeff_cost[plane] != txb_coeff_cost[plane]) {
      dump_mode_info(m);
      assert(0);
    }
  }
#endif  // CONFIG_RD_DEBUG
#else
  // PVQ writes its tokens (i.e. symbols) here.
  if (!m->mbmi.skip) {
    for (plane = 0; plane < MAX_MB_PLANE; ++plane) {
      PVQ_INFO *pvq;
      TX_SIZE tx_size =
          plane ? get_uv_tx_size(&m->mbmi, &xd->plane[plane]) : m->mbmi.tx_size;
      int idx, idy;
      const struct macroblockd_plane *const pd = &xd->plane[plane];
      int num_4x4_w;
      int num_4x4_h;
      int max_blocks_wide;
      int max_blocks_high;
      int step = (1 << tx_size);
      const int step_xy = 1 << (tx_size << 1);
      int block = 0;

      if (tx_size == TX_4X4 && bsize <= BLOCK_8X8) {
        num_4x4_w = 2 >> xd->plane[plane].subsampling_x;
        num_4x4_h = 2 >> xd->plane[plane].subsampling_y;
      } else {
        num_4x4_w =
            num_4x4_blocks_wide_lookup[bsize] >> xd->plane[plane].subsampling_x;
        num_4x4_h =
            num_4x4_blocks_high_lookup[bsize] >> xd->plane[plane].subsampling_y;
      }
      // TODO: Do we need below for 4x4,4x8,8x4 cases as well?
      max_blocks_wide =
          num_4x4_w + (xd->mb_to_right_edge >= 0
                           ? 0
                           : xd->mb_to_right_edge >> (5 + pd->subsampling_x));
      max_blocks_high =
          num_4x4_h + (xd->mb_to_bottom_edge >= 0
                           ? 0
                           : xd->mb_to_bottom_edge >> (5 + pd->subsampling_y));

      // TODO(yushin) Try to use av1_foreach_transformed_block_in_plane().
      // Logic like the mb_to_right_edge/mb_to_bottom_edge stuff should
      // really be centralized in one place.

      for (idy = 0; idy < max_blocks_high; idy += step) {
        for (idx = 0; idx < max_blocks_wide; idx += step) {
          const int is_keyframe = 0;
          const int encode_flip = 0;
          const int flip = 0;
          const int robust = 1;
          int i;
          const int has_dc_skip = 1;
          int *exg = &adapt->pvq.pvq_exg[plane][tx_size][0];
          int *ext = adapt->pvq.pvq_ext + tx_size * PVQ_MAX_PARTITIONS;
          generic_encoder *model = adapt->pvq.pvq_param_model;

          pvq = get_pvq_block(cpi->td.mb.pvq_q);

          // encode block skip info
          od_encode_cdf_adapt(&w->ec, pvq->ac_dc_coded,
                              adapt->skip_cdf[2 * tx_size + (plane != 0)], 4,
                              adapt->skip_increment);

          // AC coeffs coded?
          if (pvq->ac_dc_coded & 0x02) {
            assert(pvq->bs <= tx_size);
            for (i = 0; i < pvq->nb_bands; i++) {
              if (i == 0 || (!pvq->skip_rest &&
                             !(pvq->skip_dir & (1 << ((i - 1) % 3))))) {
                pvq_encode_partition(
                    &w->ec, pvq->qg[i], pvq->theta[i], pvq->max_theta[i],
                    pvq->y + pvq->off[i], pvq->size[i], pvq->k[i], model, adapt,
                    exg + i, ext + i, robust || is_keyframe,
                    (plane != 0) * OD_TXSIZES * PVQ_MAX_PARTITIONS +
                        pvq->bs * PVQ_MAX_PARTITIONS + i,
                    is_keyframe, i == 0 && (i < pvq->nb_bands - 1),
                    pvq->skip_rest, encode_flip, flip);
              }
              if (i == 0 && !pvq->skip_rest && pvq->bs > 0) {
                od_encode_cdf_adapt(
                    &w->ec, pvq->skip_dir,
                    &adapt->pvq
                         .pvq_skip_dir_cdf[(plane != 0) + 2 * (pvq->bs - 1)][0],
                    7, adapt->pvq.pvq_skip_dir_increment);
              }
            }
          }
          // Encode residue of DC coeff, if exist.
          if (!has_dc_skip || (pvq->ac_dc_coded & 1)) {  // DC coded?
            generic_encode(&w->ec, &adapt->model_dc[plane],
                           abs(pvq->dq_dc_residue) - has_dc_skip, -1,
                           &adapt->ex_dc[plane][pvq->bs][0], 2);
          }
          if ((pvq->ac_dc_coded & 1)) {  // DC coded?
            od_ec_enc_bits(&w->ec, pvq->dq_dc_residue < 0, 1);
          }
          block += step_xy;
        }
      }  // for (idy = 0;
    }    // for (plane =
  }      // if (!m->mbmi.skip)
#endif
}

static void write_partition(const AV1_COMMON *const cm,
                            const MACROBLOCKD *const xd, int hbs, int mi_row,
                            int mi_col, PARTITION_TYPE p, BLOCK_SIZE bsize,
                            aom_writer *w) {
  const int ctx = partition_plane_context(xd, mi_row, mi_col, bsize);
  const aom_prob *const probs = cm->fc->partition_prob[ctx];
  const int has_rows = (mi_row + hbs) < cm->mi_rows;
  const int has_cols = (mi_col + hbs) < cm->mi_cols;

  if (has_rows && has_cols) {
#if CONFIG_EXT_PARTITION_TYPES
    if (bsize <= BLOCK_8X8)
      av1_write_token(w, av1_partition_tree, probs, &partition_encodings[p]);
    else
      av1_write_token(w, av1_ext_partition_tree, probs,
                      &ext_partition_encodings[p]);
#else
#if CONFIG_DAALA_EC
    aom_write_symbol(w, p, cm->fc->partition_cdf[ctx], PARTITION_TYPES);
#else
    av1_write_token(w, av1_partition_tree, probs, &partition_encodings[p]);
#endif
#endif  // CONFIG_EXT_PARTITION_TYPES
  } else if (!has_rows && has_cols) {
    assert(p == PARTITION_SPLIT || p == PARTITION_HORZ);
    aom_write(w, p == PARTITION_SPLIT, probs[1]);
  } else if (has_rows && !has_cols) {
    assert(p == PARTITION_SPLIT || p == PARTITION_VERT);
    aom_write(w, p == PARTITION_SPLIT, probs[2]);
  } else {
    assert(p == PARTITION_SPLIT);
  }
}

#if CONFIG_SUPERTX
#define write_modes_sb_wrapper(cpi, tile, w, tok, tok_end, supertx_enabled,   \
                               mi_row, mi_col, bsize)                         \
  write_modes_sb(cpi, tile, w, tok, tok_end, supertx_enabled, mi_row, mi_col, \
                 bsize)
#else
#define write_modes_sb_wrapper(cpi, tile, w, tok, tok_end, supertx_enabled, \
                               mi_row, mi_col, bsize)                       \
  write_modes_sb(cpi, tile, w, tok, tok_end, mi_row, mi_col, bsize)
#endif  // CONFIG_SUPERTX

static void write_modes_sb(AV1_COMP *const cpi, const TileInfo *const tile,
                           aom_writer *const w, const TOKENEXTRA **tok,
                           const TOKENEXTRA *const tok_end,
#if CONFIG_SUPERTX
                           int supertx_enabled,
#endif
                           int mi_row, int mi_col, BLOCK_SIZE bsize) {
  const AV1_COMMON *const cm = &cpi->common;
  MACROBLOCKD *const xd = &cpi->td.mb.e_mbd;
  const int hbs = num_8x8_blocks_wide_lookup[bsize] / 2;
  const PARTITION_TYPE partition = get_partition(cm, mi_row, mi_col, bsize);
  const BLOCK_SIZE subsize = get_subsize(bsize, partition);
#if CONFIG_SUPERTX
  const int mi_offset = mi_row * cm->mi_stride + mi_col;
  MB_MODE_INFO *mbmi;
  const int pack_token = !supertx_enabled;
  TX_SIZE supertx_size;
  int plane;
#endif

  if (mi_row >= cm->mi_rows || mi_col >= cm->mi_cols) return;

  write_partition(cm, xd, hbs, mi_row, mi_col, partition, bsize, w);
#if CONFIG_SUPERTX
  mbmi = &cm->mi_grid_visible[mi_offset]->mbmi;
  xd->mi = cm->mi_grid_visible + mi_offset;
  set_mi_row_col(xd, tile, mi_row, num_8x8_blocks_high_lookup[bsize], mi_col,
                 num_8x8_blocks_wide_lookup[bsize], cm->mi_rows, cm->mi_cols);
  if (!supertx_enabled && !frame_is_intra_only(cm) &&
      partition != PARTITION_NONE && bsize <= MAX_SUPERTX_BLOCK_SIZE &&
      !xd->lossless[0]) {
    aom_prob prob;
    supertx_size = max_txsize_lookup[bsize];
    prob = cm->fc->supertx_prob[partition_supertx_context_lookup[partition]]
                               [supertx_size];
    supertx_enabled = (xd->mi[0]->mbmi.tx_size == supertx_size);
    aom_write(w, supertx_enabled, prob);
  }
#endif  // CONFIG_SUPERTX
  if (subsize < BLOCK_8X8) {
    write_modes_b_wrapper(cpi, tile, w, tok, tok_end, supertx_enabled, mi_row,
                          mi_col);
  } else {
    switch (partition) {
      case PARTITION_NONE:
        write_modes_b_wrapper(cpi, tile, w, tok, tok_end, supertx_enabled,
                              mi_row, mi_col);
        break;
      case PARTITION_HORZ:
        write_modes_b_wrapper(cpi, tile, w, tok, tok_end, supertx_enabled,
                              mi_row, mi_col);
        if (mi_row + hbs < cm->mi_rows)
          write_modes_b_wrapper(cpi, tile, w, tok, tok_end, supertx_enabled,
                                mi_row + hbs, mi_col);
        break;
      case PARTITION_VERT:
        write_modes_b_wrapper(cpi, tile, w, tok, tok_end, supertx_enabled,
                              mi_row, mi_col);
        if (mi_col + hbs < cm->mi_cols)
          write_modes_b_wrapper(cpi, tile, w, tok, tok_end, supertx_enabled,
                                mi_row, mi_col + hbs);
        break;
      case PARTITION_SPLIT:
        write_modes_sb_wrapper(cpi, tile, w, tok, tok_end, supertx_enabled,
                               mi_row, mi_col, subsize);
        write_modes_sb_wrapper(cpi, tile, w, tok, tok_end, supertx_enabled,
                               mi_row, mi_col + hbs, subsize);
        write_modes_sb_wrapper(cpi, tile, w, tok, tok_end, supertx_enabled,
                               mi_row + hbs, mi_col, subsize);
        write_modes_sb_wrapper(cpi, tile, w, tok, tok_end, supertx_enabled,
                               mi_row + hbs, mi_col + hbs, subsize);
        break;
#if CONFIG_EXT_PARTITION_TYPES
      case PARTITION_HORZ_A:
        write_modes_b_wrapper(cpi, tile, w, tok, tok_end, supertx_enabled,
                              mi_row, mi_col);
        write_modes_b_wrapper(cpi, tile, w, tok, tok_end, supertx_enabled,
                              mi_row, mi_col + hbs);
        write_modes_b_wrapper(cpi, tile, w, tok, tok_end, supertx_enabled,
                              mi_row + hbs, mi_col);
        break;
      case PARTITION_HORZ_B:
        write_modes_b_wrapper(cpi, tile, w, tok, tok_end, supertx_enabled,
                              mi_row, mi_col);
        write_modes_b_wrapper(cpi, tile, w, tok, tok_end, supertx_enabled,
                              mi_row + hbs, mi_col);
        write_modes_b_wrapper(cpi, tile, w, tok, tok_end, supertx_enabled,
                              mi_row + hbs, mi_col + hbs);
        break;
      case PARTITION_VERT_A:
        write_modes_b_wrapper(cpi, tile, w, tok, tok_end, supertx_enabled,
                              mi_row, mi_col);
        write_modes_b_wrapper(cpi, tile, w, tok, tok_end, supertx_enabled,
                              mi_row + hbs, mi_col);
        write_modes_b_wrapper(cpi, tile, w, tok, tok_end, supertx_enabled,
                              mi_row, mi_col + hbs);
        break;
      case PARTITION_VERT_B:
        write_modes_b_wrapper(cpi, tile, w, tok, tok_end, supertx_enabled,
                              mi_row, mi_col);
        write_modes_b_wrapper(cpi, tile, w, tok, tok_end, supertx_enabled,
                              mi_row, mi_col + hbs);
        write_modes_b_wrapper(cpi, tile, w, tok, tok_end, supertx_enabled,
                              mi_row + hbs, mi_col + hbs);
        break;
#endif  // CONFIG_EXT_PARTITION_TYPES
      default: assert(0);
    }
  }
#if CONFIG_SUPERTX
  if (partition != PARTITION_NONE && supertx_enabled && pack_token) {
    int skip;
    xd->mi = cm->mi_grid_visible + mi_offset;
    supertx_size = mbmi->tx_size;
    set_mi_row_col(xd, tile, mi_row, num_8x8_blocks_high_lookup[bsize], mi_col,
                   num_8x8_blocks_wide_lookup[bsize], cm->mi_rows, cm->mi_cols);

    assert(IMPLIES(!cm->seg.enabled, mbmi->segment_id_supertx == 0));
    assert(mbmi->segment_id_supertx < MAX_SEGMENTS);

    skip = write_skip(cm, xd, mbmi->segment_id_supertx, xd->mi[0], w);
#if CONFIG_EXT_TX
    if (get_ext_tx_types(supertx_size, bsize, 1) > 1 && !skip) {
      int eset = get_ext_tx_set(supertx_size, bsize, 1);
      if (eset > 0) {
        av1_write_token(w, av1_ext_tx_inter_tree[eset],
                        cm->fc->inter_ext_tx_prob[eset][supertx_size],
                        &ext_tx_inter_encodings[eset][mbmi->tx_type]);
      }
    }
#else
    if (supertx_size < TX_32X32 && !skip) {
      av1_write_token(w, av1_ext_tx_tree,
                      cm->fc->inter_ext_tx_prob[supertx_size],
                      &ext_tx_encodings[mbmi->tx_type]);
    }
#endif  // CONFIG_EXT_TX

    if (!skip) {
      assert(*tok < tok_end);
      for (plane = 0; plane < MAX_MB_PLANE; ++plane) {
        const int mbmi_txb_size = txsize_to_bsize[mbmi->tx_size];
        const int num_4x4_w = num_4x4_blocks_wide_lookup[mbmi_txb_size];
        const int num_4x4_h = num_4x4_blocks_high_lookup[mbmi_txb_size];
        int row, col;
        TX_SIZE tx =
            plane ? get_uv_tx_size(mbmi, &xd->plane[plane]) : mbmi->tx_size;
        BLOCK_SIZE txb_size = txsize_to_bsize[tx];
        int bw = num_4x4_blocks_wide_lookup[txb_size];

        TOKEN_STATS token_stats;
        token_stats.cost = 0;
        for (row = 0; row < num_4x4_h; row += bw)
          for (col = 0; col < num_4x4_w; col += bw)
            pack_mb_tokens(w, tok, tok_end, cm->bit_depth, tx, &token_stats);
        assert(*tok < tok_end && (*tok)->token == EOSB_TOKEN);
        (*tok)++;
      }
    }
  }
#endif  // CONFIG_SUPERTX

// update partition context
#if CONFIG_EXT_PARTITION_TYPES
  update_ext_partition_context(xd, mi_row, mi_col, subsize, bsize, partition);
#else
  if (bsize >= BLOCK_8X8 &&
      (bsize == BLOCK_8X8 || partition != PARTITION_SPLIT))
    update_partition_context(xd, mi_row, mi_col, subsize, bsize);
#endif  // CONFIG_EXT_PARTITION_TYPES

#if CONFIG_DERING
  if (bsize == BLOCK_64X64 && cm->dering_level != 0 &&
      !sb_all_skip(cm, mi_row, mi_col)) {
    aom_write_literal(
        w,
        cm->mi_grid_visible[mi_row * cm->mi_stride + mi_col]->mbmi.dering_gain,
        DERING_REFINEMENT_BITS);
  }
#endif

#if CONFIG_CLPF
  if (bsize == BLOCK_64X64 && cm->clpf_blocks && cm->clpf_strength_y &&
      cm->clpf_size != CLPF_NOSIZE) {
    const int tl = mi_row * MI_SIZE / MIN_FB_SIZE * cm->clpf_stride +
                   mi_col * MI_SIZE / MIN_FB_SIZE;
    const int tr = tl + 1;
    const int bl = tl + cm->clpf_stride;
    const int br = tr + cm->clpf_stride;

    // Up to four bits per SB.
    // When clpf_size indicates a size larger than the SB size
    // (CLPF_128X128), one bit for every fourth SB will be transmitted
    // regardless of skip blocks.
    if (cm->clpf_blocks[tl] != CLPF_NOFLAG)
      aom_write_literal(w, cm->clpf_blocks[tl], 1);

    if (mi_col + MI_SIZE / 2 < cm->mi_cols &&
        cm->clpf_blocks[tr] != CLPF_NOFLAG)
      aom_write_literal(w, cm->clpf_blocks[tr], 1);

    if (mi_row + MI_SIZE / 2 < cm->mi_rows &&
        cm->clpf_blocks[bl] != CLPF_NOFLAG)
      aom_write_literal(w, cm->clpf_blocks[bl], 1);

    if (mi_row + MI_SIZE / 2 < cm->mi_rows &&
        mi_col + MI_SIZE / 2 < cm->mi_cols &&
        cm->clpf_blocks[br] != CLPF_NOFLAG)
      aom_write_literal(w, cm->clpf_blocks[br], 1);
  }
#endif
}

static void write_modes(AV1_COMP *const cpi, const TileInfo *const tile,
                        aom_writer *const w, const TOKENEXTRA **tok,
                        const TOKENEXTRA *const tok_end) {
  AV1_COMMON *const cm = &cpi->common;
  MACROBLOCKD *const xd = &cpi->td.mb.e_mbd;
  const int mi_row_start = tile->mi_row_start;
  const int mi_row_end = tile->mi_row_end;
  const int mi_col_start = tile->mi_col_start;
  const int mi_col_end = tile->mi_col_end;
  int mi_row, mi_col;
  av1_zero_above_context(cm, mi_col_start, mi_col_end);
#if CONFIG_PVQ
  assert(cpi->td.mb.pvq_q->curr_pos == 0);
#endif
#if CONFIG_DELTA_Q
  if (cpi->common.delta_q_present_flag) {
    xd->prev_qindex = cpi->common.base_qindex;
  }
#endif

  for (mi_row = mi_row_start; mi_row < mi_row_end; mi_row += cm->mib_size) {
    av1_zero_left_context(xd);

    for (mi_col = mi_col_start; mi_col < mi_col_end; mi_col += cm->mib_size) {
      write_modes_sb_wrapper(cpi, tile, w, tok, tok_end, 0, mi_row, mi_col,
                             cm->sb_size);
    }
  }
#if CONFIG_PVQ
  // Check that the number of PVQ blocks encoded and written to the bitstream
  // are the same
  assert(cpi->td.mb.pvq_q->curr_pos == cpi->td.mb.pvq_q->last_pos);
  // Reset curr_pos in case we repack the bitstream
  cpi->td.mb.pvq_q->curr_pos = 0;
#endif
}

#if !CONFIG_PVQ
static void build_tree_distribution(AV1_COMP *cpi, TX_SIZE tx_size,
                                    av1_coeff_stats *coef_branch_ct,
                                    av1_coeff_probs_model *coef_probs) {
  av1_coeff_count *coef_counts = cpi->td.rd_counts.coef_counts[tx_size];
  unsigned int(*eob_branch_ct)[REF_TYPES][COEF_BANDS][COEFF_CONTEXTS] =
      cpi->common.counts.eob_branch[tx_size];
  int i, j, k, l, m;

  for (i = 0; i < PLANE_TYPES; ++i) {
    for (j = 0; j < REF_TYPES; ++j) {
      for (k = 0; k < COEF_BANDS; ++k) {
        for (l = 0; l < BAND_COEFF_CONTEXTS(k); ++l) {
          av1_tree_probs_from_distribution(av1_coef_tree,
                                           coef_branch_ct[i][j][k][l],
                                           coef_counts[i][j][k][l]);
          coef_branch_ct[i][j][k][l][0][1] =
              eob_branch_ct[i][j][k][l] - coef_branch_ct[i][j][k][l][0][0];
          for (m = 0; m < UNCONSTRAINED_NODES; ++m)
            coef_probs[i][j][k][l][m] =
                get_binary_prob(coef_branch_ct[i][j][k][l][m][0],
                                coef_branch_ct[i][j][k][l][m][1]);
        }
      }
    }
  }
}

static void update_coef_probs_common(aom_writer *const bc, AV1_COMP *cpi,
                                     TX_SIZE tx_size,
                                     av1_coeff_stats *frame_branch_ct,
                                     av1_coeff_probs_model *new_coef_probs) {
  av1_coeff_probs_model *old_coef_probs = cpi->common.fc->coef_probs[tx_size];
  const aom_prob upd = DIFF_UPDATE_PROB;
#if CONFIG_EC_ADAPT
  const int entropy_nodes_update = UNCONSTRAINED_NODES - 1;
#else
  const int entropy_nodes_update = UNCONSTRAINED_NODES;
#endif
  int i, j, k, l, t;
  int stepsize = cpi->sf.coeff_prob_appx_step;
#if CONFIG_TILE_GROUPS
  const int probwt = cpi->common.num_tg;
#else
  const int probwt = 1;
#endif

  switch (cpi->sf.use_fast_coef_updates) {
    case TWO_LOOP: {
      /* dry run to see if there is any update at all needed */
      int savings = 0;
      int update[2] = { 0, 0 };
      for (i = 0; i < PLANE_TYPES; ++i) {
        for (j = 0; j < REF_TYPES; ++j) {
          for (k = 0; k < COEF_BANDS; ++k) {
            for (l = 0; l < BAND_COEFF_CONTEXTS(k); ++l) {
              for (t = 0; t < entropy_nodes_update; ++t) {
                aom_prob newp = new_coef_probs[i][j][k][l][t];
                const aom_prob oldp = old_coef_probs[i][j][k][l][t];
                int s;
                int u = 0;
                if (t == PIVOT_NODE)
                  s = av1_prob_diff_update_savings_search_model(
                      frame_branch_ct[i][j][k][l][0],
                      old_coef_probs[i][j][k][l], &newp, upd, stepsize, probwt);
                else
                  s = av1_prob_diff_update_savings_search(
                      frame_branch_ct[i][j][k][l][t], oldp, &newp, upd, probwt);

                if (s > 0 && newp != oldp) u = 1;
                if (u)
                  savings += s - (int)(av1_cost_zero(upd));
                else
                  savings -= (int)(av1_cost_zero(upd));
                update[u]++;
              }
            }
          }
        }
      }

      /* Is coef updated at all */
      if (update[1] == 0 || savings < 0) {
        aom_write_bit(bc, 0);
        return;
      }
      aom_write_bit(bc, 1);
      for (i = 0; i < PLANE_TYPES; ++i) {
        for (j = 0; j < REF_TYPES; ++j) {
          for (k = 0; k < COEF_BANDS; ++k) {
            for (l = 0; l < BAND_COEFF_CONTEXTS(k); ++l) {
              // calc probs and branch cts for this frame only
              for (t = 0; t < entropy_nodes_update; ++t) {
                aom_prob newp = new_coef_probs[i][j][k][l][t];
                aom_prob *oldp = old_coef_probs[i][j][k][l] + t;
                int s;
                int u = 0;
                if (t == PIVOT_NODE)
                  s = av1_prob_diff_update_savings_search_model(
                      frame_branch_ct[i][j][k][l][0],
                      old_coef_probs[i][j][k][l], &newp, upd, stepsize, probwt);
                else
                  s = av1_prob_diff_update_savings_search(
                      frame_branch_ct[i][j][k][l][t], *oldp, &newp, upd,
                      probwt);
                if (s > 0 && newp != *oldp) u = 1;
                aom_write(bc, u, upd);
                if (u) {
                  /* send/use new probability */
                  av1_write_prob_diff_update(bc, newp, *oldp);
                  *oldp = newp;
                }
              }
            }
          }
        }
      }
      return;
    }

    case ONE_LOOP_REDUCED: {
      int updates = 0;
      int noupdates_before_first = 0;
      for (i = 0; i < PLANE_TYPES; ++i) {
        for (j = 0; j < REF_TYPES; ++j) {
          for (k = 0; k < COEF_BANDS; ++k) {
            for (l = 0; l < BAND_COEFF_CONTEXTS(k); ++l) {
              // calc probs and branch cts for this frame only
              for (t = 0; t < entropy_nodes_update; ++t) {
                aom_prob newp = new_coef_probs[i][j][k][l][t];
                aom_prob *oldp = old_coef_probs[i][j][k][l] + t;
                int s;
                int u = 0;
                if (t == PIVOT_NODE) {
                  s = av1_prob_diff_update_savings_search_model(
                      frame_branch_ct[i][j][k][l][0],
                      old_coef_probs[i][j][k][l], &newp, upd, stepsize, probwt);
                } else {
                  s = av1_prob_diff_update_savings_search(
                      frame_branch_ct[i][j][k][l][t], *oldp, &newp, upd,
                      probwt);
                }

                if (s > 0 && newp != *oldp) u = 1;
                updates += u;
                if (u == 0 && updates == 0) {
                  noupdates_before_first++;
                  continue;
                }
                if (u == 1 && updates == 1) {
                  int v;
                  // first update
                  aom_write_bit(bc, 1);
                  for (v = 0; v < noupdates_before_first; ++v)
                    aom_write(bc, 0, upd);
                }
                aom_write(bc, u, upd);
                if (u) {
                  /* send/use new probability */
                  av1_write_prob_diff_update(bc, newp, *oldp);
                  *oldp = newp;
                }
              }
            }
          }
        }
      }
      if (updates == 0) {
        aom_write_bit(bc, 0);  // no updates
      }
      return;
    }
    default: assert(0);
  }
}

#if CONFIG_ENTROPY
// Calculate the token counts between subsequent subframe updates.
static void get_coef_counts_diff(AV1_COMP *cpi, int index,
                                 av1_coeff_count coef_counts[TX_SIZES]
                                                            [PLANE_TYPES],
                                 unsigned int eob_counts[TX_SIZES][PLANE_TYPES]
                                                        [REF_TYPES][COEF_BANDS]
                                                        [COEFF_CONTEXTS]) {
  int i, j, k, l, m, tx_size, val;
  const int max_idx = cpi->common.coef_probs_update_idx;
  const TX_MODE tx_mode = cpi->common.tx_mode;
  const int max_tx_size = tx_mode_to_biggest_tx_size[tx_mode];
  const SUBFRAME_STATS *subframe_stats = &cpi->subframe_stats;

  assert(max_idx < COEF_PROBS_BUFS);

  for (tx_size = TX_4X4; tx_size <= max_tx_size; ++tx_size)
    for (i = 0; i < PLANE_TYPES; ++i)
      for (j = 0; j < REF_TYPES; ++j)
        for (k = 0; k < COEF_BANDS; ++k)
          for (l = 0; l < BAND_COEFF_CONTEXTS(k); ++l) {
            if (index == max_idx) {
              val =
                  cpi->common.counts.eob_branch[tx_size][i][j][k][l] -
                  subframe_stats->eob_counts_buf[max_idx][tx_size][i][j][k][l];
            } else {
              val = subframe_stats->eob_counts_buf[index + 1][tx_size][i][j][k]
                                                  [l] -
                    subframe_stats->eob_counts_buf[index][tx_size][i][j][k][l];
            }
            assert(val >= 0);
            eob_counts[tx_size][i][j][k][l] = val;

            for (m = 0; m < ENTROPY_TOKENS; ++m) {
              if (index == max_idx) {
                val = cpi->td.rd_counts.coef_counts[tx_size][i][j][k][l][m] -
                      subframe_stats->coef_counts_buf[max_idx][tx_size][i][j][k]
                                                     [l][m];
              } else {
                val = subframe_stats->coef_counts_buf[index + 1][tx_size][i][j]
                                                     [k][l][m] -
                      subframe_stats->coef_counts_buf[index][tx_size][i][j][k]
                                                     [l][m];
              }
              assert(val >= 0);
              coef_counts[tx_size][i][j][k][l][m] = val;
            }
          }
}

static void update_coef_probs_subframe(
    aom_writer *const bc, AV1_COMP *cpi, TX_SIZE tx_size,
    av1_coeff_stats branch_ct[COEF_PROBS_BUFS][TX_SIZES][PLANE_TYPES],
    av1_coeff_probs_model *new_coef_probs) {
  av1_coeff_probs_model *old_coef_probs = cpi->common.fc->coef_probs[tx_size];
  const aom_prob upd = DIFF_UPDATE_PROB;
  const int entropy_nodes_update = UNCONSTRAINED_NODES;
  int i, j, k, l, t;
  int stepsize = cpi->sf.coeff_prob_appx_step;
  const int max_idx = cpi->common.coef_probs_update_idx;
  int idx;
  unsigned int this_branch_ct[ENTROPY_NODES][COEF_PROBS_BUFS][2];

  switch (cpi->sf.use_fast_coef_updates) {
    case TWO_LOOP: {
      /* dry run to see if there is any update at all needed */
      int savings = 0;
      int update[2] = { 0, 0 };
      for (i = 0; i < PLANE_TYPES; ++i) {
        for (j = 0; j < REF_TYPES; ++j) {
          for (k = 0; k < COEF_BANDS; ++k) {
            for (l = 0; l < BAND_COEFF_CONTEXTS(k); ++l) {
              for (t = 0; t < ENTROPY_NODES; ++t) {
                for (idx = 0; idx <= max_idx; ++idx) {
                  memcpy(this_branch_ct[t][idx],
                         branch_ct[idx][tx_size][i][j][k][l][t],
                         2 * sizeof(this_branch_ct[t][idx][0]));
                }
              }
              for (t = 0; t < entropy_nodes_update; ++t) {
                aom_prob newp = new_coef_probs[i][j][k][l][t];
                const aom_prob oldp = old_coef_probs[i][j][k][l][t];
                int s, u = 0;

                if (t == PIVOT_NODE)
                  s = av1_prob_update_search_model_subframe(
                      this_branch_ct, old_coef_probs[i][j][k][l], &newp, upd,
                      stepsize, max_idx);
                else
                  s = av1_prob_update_search_subframe(this_branch_ct[t], oldp,
                                                      &newp, upd, max_idx);
                if (s > 0 && newp != oldp) u = 1;
                if (u)
                  savings += s - (int)(av1_cost_zero(upd));
                else
                  savings -= (int)(av1_cost_zero(upd));
                update[u]++;
              }
            }
          }
        }
      }

      /* Is coef updated at all */
      if (update[1] == 0 || savings < 0) {
        aom_write_bit(bc, 0);
        return;
      }
      aom_write_bit(bc, 1);
      for (i = 0; i < PLANE_TYPES; ++i) {
        for (j = 0; j < REF_TYPES; ++j) {
          for (k = 0; k < COEF_BANDS; ++k) {
            for (l = 0; l < BAND_COEFF_CONTEXTS(k); ++l) {
              for (t = 0; t < ENTROPY_NODES; ++t) {
                for (idx = 0; idx <= max_idx; ++idx) {
                  memcpy(this_branch_ct[t][idx],
                         branch_ct[idx][tx_size][i][j][k][l][t],
                         2 * sizeof(this_branch_ct[t][idx][0]));
                }
              }
              for (t = 0; t < entropy_nodes_update; ++t) {
                aom_prob newp = new_coef_probs[i][j][k][l][t];
                aom_prob *oldp = old_coef_probs[i][j][k][l] + t;
                int s;
                int u = 0;

                if (t == PIVOT_NODE)
                  s = av1_prob_update_search_model_subframe(
                      this_branch_ct, old_coef_probs[i][j][k][l], &newp, upd,
                      stepsize, max_idx);
                else
                  s = av1_prob_update_search_subframe(this_branch_ct[t], *oldp,
                                                      &newp, upd, max_idx);
                if (s > 0 && newp != *oldp) u = 1;
                aom_write(bc, u, upd);
                if (u) {
                  /* send/use new probability */
                  av1_write_prob_diff_update(bc, newp, *oldp);
                  *oldp = newp;
                }
              }
            }
          }
        }
      }
      return;
    }

    case ONE_LOOP_REDUCED: {
      int updates = 0;
      int noupdates_before_first = 0;
      for (i = 0; i < PLANE_TYPES; ++i) {
        for (j = 0; j < REF_TYPES; ++j) {
          for (k = 0; k < COEF_BANDS; ++k) {
            for (l = 0; l < BAND_COEFF_CONTEXTS(k); ++l) {
              for (t = 0; t < ENTROPY_NODES; ++t) {
                for (idx = 0; idx <= max_idx; ++idx) {
                  memcpy(this_branch_ct[t][idx],
                         branch_ct[idx][tx_size][i][j][k][l][t],
                         2 * sizeof(this_branch_ct[t][idx][0]));
                }
              }
              for (t = 0; t < entropy_nodes_update; ++t) {
                aom_prob newp = new_coef_probs[i][j][k][l][t];
                aom_prob *oldp = old_coef_probs[i][j][k][l] + t;
                int s;
                int u = 0;

                if (t == PIVOT_NODE)
                  s = av1_prob_update_search_model_subframe(
                      this_branch_ct, old_coef_probs[i][j][k][l], &newp, upd,
                      stepsize, max_idx);
                else
                  s = av1_prob_update_search_subframe(this_branch_ct[t], *oldp,
                                                      &newp, upd, max_idx);
                if (s > 0 && newp != *oldp) u = 1;
                updates += u;
                if (u == 0 && updates == 0) {
                  noupdates_before_first++;
                  continue;
                }
                if (u == 1 && updates == 1) {
                  int v;
                  // first update
                  aom_write_bit(bc, 1);
                  for (v = 0; v < noupdates_before_first; ++v)
                    aom_write(bc, 0, upd);
                }
                aom_write(bc, u, upd);
                if (u) {
                  /* send/use new probability */
                  av1_write_prob_diff_update(bc, newp, *oldp);
                  *oldp = newp;
                }
              }
            }
          }
        }
      }
      if (updates == 0) {
        aom_write_bit(bc, 0);  // no updates
      }
      return;
    }
    default: assert(0);
  }
}
#endif  // CONFIG_ENTROPY

static void update_coef_probs(AV1_COMP *cpi, aom_writer *w) {
  const TX_MODE tx_mode = cpi->common.tx_mode;
  const TX_SIZE max_tx_size = tx_mode_to_biggest_tx_size[tx_mode];
  TX_SIZE tx_size;
#if CONFIG_RANS
  int update = 0;
#endif  // CONFIG_RANS
#if CONFIG_ENTROPY
  AV1_COMMON *cm = &cpi->common;
  SUBFRAME_STATS *subframe_stats = &cpi->subframe_stats;
  int i;
  av1_coeff_probs_model dummy_frame_coef_probs[PLANE_TYPES];

  if (cm->do_subframe_update &&
      cm->refresh_frame_context == REFRESH_FRAME_CONTEXT_BACKWARD) {
    av1_copy(cpi->common.fc->coef_probs,
             subframe_stats->enc_starting_coef_probs);
    for (i = 0; i <= cpi->common.coef_probs_update_idx; ++i) {
      get_coef_counts_diff(cpi, i, cpi->wholeframe_stats.coef_counts_buf[i],
                           cpi->wholeframe_stats.eob_counts_buf[i]);
    }
  }
#endif  // CONFIG_ENTROPY

  for (tx_size = TX_4X4; tx_size <= max_tx_size; ++tx_size) {
    av1_coeff_stats frame_branch_ct[PLANE_TYPES];
    av1_coeff_probs_model frame_coef_probs[PLANE_TYPES];
    if (cpi->td.counts->tx_size_totals[tx_size] <= 20 ||
        (tx_size >= TX_16X16 && cpi->sf.tx_size_search_method == USE_TX_8X8)) {
      aom_write_bit(w, 0);
    } else {
#if CONFIG_ENTROPY
      if (cm->do_subframe_update &&
          cm->refresh_frame_context == REFRESH_FRAME_CONTEXT_BACKWARD) {
        unsigned int this_eob_counts_copy[PLANE_TYPES][REF_TYPES][COEF_BANDS]
                                         [COEFF_CONTEXTS];
        av1_coeff_count coef_counts_copy[PLANE_TYPES];
        av1_copy(this_eob_counts_copy, cpi->common.counts.eob_branch[tx_size]);
        av1_copy(coef_counts_copy, cpi->td.rd_counts.coef_counts[tx_size]);
        build_tree_distribution(cpi, tx_size, frame_branch_ct,
                                frame_coef_probs);
        for (i = 0; i <= cpi->common.coef_probs_update_idx; ++i) {
          av1_copy(cpi->common.counts.eob_branch[tx_size],
                   cpi->wholeframe_stats.eob_counts_buf[i][tx_size]);
          av1_copy(cpi->td.rd_counts.coef_counts[tx_size],
                   cpi->wholeframe_stats.coef_counts_buf[i][tx_size]);
          build_tree_distribution(cpi, tx_size, cpi->branch_ct_buf[i][tx_size],
                                  dummy_frame_coef_probs);
        }
        av1_copy(cpi->common.counts.eob_branch[tx_size], this_eob_counts_copy);
        av1_copy(cpi->td.rd_counts.coef_counts[tx_size], coef_counts_copy);

        update_coef_probs_subframe(w, cpi, tx_size, cpi->branch_ct_buf,
                                   frame_coef_probs);
#if CONFIG_RANS
        update = 1;
#endif  // CONFIG_RANS
      } else {
#endif  // CONFIG_ENTROPY
        build_tree_distribution(cpi, tx_size, frame_branch_ct,
                                frame_coef_probs);
        update_coef_probs_common(w, cpi, tx_size, frame_branch_ct,
                                 frame_coef_probs);
#if CONFIG_RANS
        update = 1;
#endif  // CONFIG_RANS
#if CONFIG_ENTROPY
      }
#endif  // CONFIG_ENTROPY
    }
  }

#if CONFIG_ENTROPY
  av1_copy(cm->starting_coef_probs, cm->fc->coef_probs);
  av1_copy(subframe_stats->coef_probs_buf[0], cm->fc->coef_probs);
  if (cm->do_subframe_update &&
      cm->refresh_frame_context == REFRESH_FRAME_CONTEXT_BACKWARD) {
    unsigned int eob_counts_copy[TX_SIZES][PLANE_TYPES][REF_TYPES][COEF_BANDS]
                                [COEFF_CONTEXTS];
    av1_copy(eob_counts_copy, cm->counts.eob_branch);
    for (i = 1; i <= cpi->common.coef_probs_update_idx; ++i) {
      for (tx_size = TX_4X4; tx_size <= max_tx_size; ++tx_size)
        av1_full_to_model_counts(cm->counts.coef[tx_size],
                                 subframe_stats->coef_counts_buf[i][tx_size]);
      av1_copy(cm->counts.eob_branch, subframe_stats->eob_counts_buf[i]);
      av1_partial_adapt_probs(cm, 0, 0);
      av1_copy(subframe_stats->coef_probs_buf[i], cm->fc->coef_probs);
    }
    av1_copy(cm->fc->coef_probs, subframe_stats->coef_probs_buf[0]);
    av1_copy(cm->counts.eob_branch, eob_counts_copy);
  }
#endif  // CONFIG_ENTROPY
#if CONFIG_RANS
  if (update) av1_coef_pareto_cdfs(cpi->common.fc);
#endif  // CONFIG_RANS
}
#endif

#if CONFIG_LOOP_RESTORATION
static void encode_restoration_mode(AV1_COMMON *cm,
                                    struct aom_write_bit_buffer *wb) {
  RestorationInfo *rst = &cm->rst_info;
  switch (rst->frame_restoration_type) {
    case RESTORE_NONE:
      aom_wb_write_bit(wb, 0);
      aom_wb_write_bit(wb, 0);
      break;
    case RESTORE_SWITCHABLE:
      aom_wb_write_bit(wb, 0);
      aom_wb_write_bit(wb, 1);
      break;
    case RESTORE_BILATERAL:
      aom_wb_write_bit(wb, 1);
      aom_wb_write_bit(wb, 0);
      break;
    case RESTORE_WIENER:
      aom_wb_write_bit(wb, 1);
      aom_wb_write_bit(wb, 1);
      break;
    default: assert(0);
  }
}

static void encode_restoration(AV1_COMMON *cm, aom_writer *wb) {
  int i;
  RestorationInfo *rsi = &cm->rst_info;
  if (rsi->frame_restoration_type != RESTORE_NONE) {
    if (rsi->frame_restoration_type == RESTORE_SWITCHABLE) {
      // RESTORE_SWITCHABLE
      for (i = 0; i < cm->rst_internal.ntiles; ++i) {
        av1_write_token(
            wb, av1_switchable_restore_tree, cm->fc->switchable_restore_prob,
            &switchable_restore_encodings[rsi->restoration_type[i]]);
        if (rsi->restoration_type[i] == RESTORE_BILATERAL) {
          int s;
          for (s = 0; s < BILATERAL_SUBTILES; ++s) {
#if BILATERAL_SUBTILES == 0
            aom_write_literal(wb, rsi->bilateral_info[i].level[s],
                              av1_bilateral_level_bits(cm));
#else
            aom_write(wb, rsi->bilateral_info[i].level[s] >= 0,
                      RESTORE_NONE_BILATERAL_PROB);
            if (rsi->bilateral_info[i].level[s] >= 0) {
              aom_write_literal(wb, rsi->bilateral_info[i].level[s],
                                av1_bilateral_level_bits(cm));
            }
#endif
          }
        } else if (rsi->restoration_type[i] == RESTORE_WIENER) {
          aom_write_literal(
              wb, rsi->wiener_info[i].vfilter[0] - WIENER_FILT_TAP0_MINV,
              WIENER_FILT_TAP0_BITS);
          aom_write_literal(
              wb, rsi->wiener_info[i].vfilter[1] - WIENER_FILT_TAP1_MINV,
              WIENER_FILT_TAP1_BITS);
          aom_write_literal(
              wb, rsi->wiener_info[i].vfilter[2] - WIENER_FILT_TAP2_MINV,
              WIENER_FILT_TAP2_BITS);
          aom_write_literal(
              wb, rsi->wiener_info[i].hfilter[0] - WIENER_FILT_TAP0_MINV,
              WIENER_FILT_TAP0_BITS);
          aom_write_literal(
              wb, rsi->wiener_info[i].hfilter[1] - WIENER_FILT_TAP1_MINV,
              WIENER_FILT_TAP1_BITS);
          aom_write_literal(
              wb, rsi->wiener_info[i].hfilter[2] - WIENER_FILT_TAP2_MINV,
              WIENER_FILT_TAP2_BITS);
        }
      }
    } else if (rsi->frame_restoration_type == RESTORE_BILATERAL) {
      for (i = 0; i < cm->rst_internal.ntiles; ++i) {
        int s;
        for (s = 0; s < BILATERAL_SUBTILES; ++s) {
          aom_write(wb, rsi->bilateral_info[i].level[s] >= 0,
                    RESTORE_NONE_BILATERAL_PROB);
          if (rsi->bilateral_info[i].level[s] >= 0) {
            aom_write_literal(wb, rsi->bilateral_info[i].level[s],
                              av1_bilateral_level_bits(cm));
          }
        }
      }
    } else if (rsi->frame_restoration_type == RESTORE_WIENER) {
      for (i = 0; i < cm->rst_internal.ntiles; ++i) {
        aom_write(wb, rsi->wiener_info[i].level != 0, RESTORE_NONE_WIENER_PROB);
        if (rsi->wiener_info[i].level) {
          aom_write_literal(
              wb, rsi->wiener_info[i].vfilter[0] - WIENER_FILT_TAP0_MINV,
              WIENER_FILT_TAP0_BITS);
          aom_write_literal(
              wb, rsi->wiener_info[i].vfilter[1] - WIENER_FILT_TAP1_MINV,
              WIENER_FILT_TAP1_BITS);
          aom_write_literal(
              wb, rsi->wiener_info[i].vfilter[2] - WIENER_FILT_TAP2_MINV,
              WIENER_FILT_TAP2_BITS);
          aom_write_literal(
              wb, rsi->wiener_info[i].hfilter[0] - WIENER_FILT_TAP0_MINV,
              WIENER_FILT_TAP0_BITS);
          aom_write_literal(
              wb, rsi->wiener_info[i].hfilter[1] - WIENER_FILT_TAP1_MINV,
              WIENER_FILT_TAP1_BITS);
          aom_write_literal(
              wb, rsi->wiener_info[i].hfilter[2] - WIENER_FILT_TAP2_MINV,
              WIENER_FILT_TAP2_BITS);
        }
      }
    }
  }
}
#endif  // CONFIG_LOOP_RESTORATION

static void encode_loopfilter(AV1_COMMON *cm, struct aom_write_bit_buffer *wb) {
  int i;
  struct loopfilter *lf = &cm->lf;

  // Encode the loop filter level and type
  aom_wb_write_literal(wb, lf->filter_level, 6);
  aom_wb_write_literal(wb, lf->sharpness_level, 3);

  // Write out loop filter deltas applied at the MB level based on mode or
  // ref frame (if they are enabled).
  aom_wb_write_bit(wb, lf->mode_ref_delta_enabled);

  if (lf->mode_ref_delta_enabled) {
    aom_wb_write_bit(wb, lf->mode_ref_delta_update);
    if (lf->mode_ref_delta_update) {
      for (i = 0; i < TOTAL_REFS_PER_FRAME; i++) {
        const int delta = lf->ref_deltas[i];
        const int changed = delta != lf->last_ref_deltas[i];
        aom_wb_write_bit(wb, changed);
        if (changed) {
          lf->last_ref_deltas[i] = delta;
          aom_wb_write_inv_signed_literal(wb, delta, 6);
        }
      }

      for (i = 0; i < MAX_MODE_LF_DELTAS; i++) {
        const int delta = lf->mode_deltas[i];
        const int changed = delta != lf->last_mode_deltas[i];
        aom_wb_write_bit(wb, changed);
        if (changed) {
          lf->last_mode_deltas[i] = delta;
          aom_wb_write_inv_signed_literal(wb, delta, 6);
        }
      }
    }
  }
}

#if CONFIG_CLPF
static void encode_clpf(const AV1_COMMON *cm, struct aom_write_bit_buffer *wb) {
  aom_wb_write_literal(wb, cm->clpf_strength_y, 2);
  aom_wb_write_literal(wb, cm->clpf_strength_u, 2);
  aom_wb_write_literal(wb, cm->clpf_strength_v, 2);
  if (cm->clpf_strength_y) {
    aom_wb_write_literal(wb, cm->clpf_size, 2);
  }
}
#endif

#if CONFIG_DERING
static void encode_dering(int level, struct aom_write_bit_buffer *wb) {
  aom_wb_write_literal(wb, level, DERING_LEVEL_BITS);
}
#endif  // CONFIG_DERING

static void write_delta_q(struct aom_write_bit_buffer *wb, int delta_q) {
  if (delta_q != 0) {
    aom_wb_write_bit(wb, 1);
    aom_wb_write_inv_signed_literal(wb, delta_q, 6);
  } else {
    aom_wb_write_bit(wb, 0);
  }
}

static void encode_quantization(const AV1_COMMON *const cm,
                                struct aom_write_bit_buffer *wb) {
  aom_wb_write_literal(wb, cm->base_qindex, QINDEX_BITS);
  write_delta_q(wb, cm->y_dc_delta_q);
  write_delta_q(wb, cm->uv_dc_delta_q);
  write_delta_q(wb, cm->uv_ac_delta_q);
#if CONFIG_AOM_QM
  aom_wb_write_bit(wb, cm->using_qmatrix);
  if (cm->using_qmatrix) {
    aom_wb_write_literal(wb, cm->min_qmlevel, QM_LEVEL_BITS);
    aom_wb_write_literal(wb, cm->max_qmlevel, QM_LEVEL_BITS);
  }
#endif
}

static void encode_segmentation(AV1_COMMON *cm, MACROBLOCKD *xd,
                                struct aom_write_bit_buffer *wb) {
  int i, j;
  const struct segmentation *seg = &cm->seg;

  aom_wb_write_bit(wb, seg->enabled);
  if (!seg->enabled) return;

  // Segmentation map
  if (!frame_is_intra_only(cm) && !cm->error_resilient_mode) {
    aom_wb_write_bit(wb, seg->update_map);
  } else {
    assert(seg->update_map == 1);
  }
  if (seg->update_map) {
    // Select the coding strategy (temporal or spatial)
    av1_choose_segmap_coding_method(cm, xd);

    // Write out the chosen coding method.
    if (!frame_is_intra_only(cm) && !cm->error_resilient_mode) {
      aom_wb_write_bit(wb, seg->temporal_update);
    } else {
      assert(seg->temporal_update == 0);
    }
  }

  // Segmentation data
  aom_wb_write_bit(wb, seg->update_data);
  if (seg->update_data) {
    aom_wb_write_bit(wb, seg->abs_delta);

    for (i = 0; i < MAX_SEGMENTS; i++) {
      for (j = 0; j < SEG_LVL_MAX; j++) {
        const int active = segfeature_active(seg, i, j);
        aom_wb_write_bit(wb, active);
        if (active) {
          const int data = get_segdata(seg, i, j);
          const int data_max = av1_seg_feature_data_max(j);

          if (av1_is_segfeature_signed(j)) {
            encode_unsigned_max(wb, abs(data), data_max);
            aom_wb_write_bit(wb, data < 0);
          } else {
            encode_unsigned_max(wb, data, data_max);
          }
        }
      }
    }
  }
}

#if !CONFIG_EC_ADAPT
static void update_seg_probs(AV1_COMP *cpi, aom_writer *w) {
  AV1_COMMON *cm = &cpi->common;
#if CONFIG_TILE_GROUPS
  const int probwt = cm->num_tg;
#else
  const int probwt = 1;
#endif

  if (!cm->seg.enabled || !cm->seg.update_map) return;

  if (cm->seg.temporal_update) {
    int i;

    for (i = 0; i < PREDICTION_PROBS; i++)
      av1_cond_prob_diff_update(w, &cm->fc->seg.pred_probs[i],
                                cm->counts.seg.pred[i], probwt);

    prob_diff_update(av1_segment_tree, cm->fc->seg.tree_probs,
                     cm->counts.seg.tree_mispred, MAX_SEGMENTS, probwt, w);
  } else {
    prob_diff_update(av1_segment_tree, cm->fc->seg.tree_probs,
                     cm->counts.seg.tree_total, MAX_SEGMENTS, probwt, w);
  }
}
#endif

static void write_txfm_mode(TX_MODE mode, struct aom_write_bit_buffer *wb) {
  aom_wb_write_bit(wb, mode == TX_MODE_SELECT);
  if (mode != TX_MODE_SELECT) aom_wb_write_literal(wb, mode, 2);
}

static void update_txfm_probs(AV1_COMMON *cm, aom_writer *w,
                              FRAME_COUNTS *counts) {
#if CONFIG_TILE_GROUPS
  const int probwt = cm->num_tg;
#else
  const int probwt = 1;
#endif
  if (cm->tx_mode == TX_MODE_SELECT) {
    int i, j;
    for (i = 0; i < MAX_TX_DEPTH; ++i)
      for (j = 0; j < TX_SIZE_CONTEXTS; ++j)
        prob_diff_update(av1_tx_size_tree[i], cm->fc->tx_size_probs[i][j],
                         counts->tx_size[i][j], i + 2, probwt, w);
  }
}

static void write_interp_filter(InterpFilter filter,
                                struct aom_write_bit_buffer *wb) {
  aom_wb_write_bit(wb, filter == SWITCHABLE);
  if (filter != SWITCHABLE)
    aom_wb_write_literal(wb, filter, LOG_SWITCHABLE_FILTERS);
}

static void fix_interp_filter(AV1_COMMON *cm, FRAME_COUNTS *counts) {
  if (cm->interp_filter == SWITCHABLE) {
    // Check to see if only one of the filters is actually used
    int count[SWITCHABLE_FILTERS];
    int i, j, c = 0;
    for (i = 0; i < SWITCHABLE_FILTERS; ++i) {
      count[i] = 0;
      for (j = 0; j < SWITCHABLE_FILTER_CONTEXTS; ++j)
        count[i] += counts->switchable_interp[j][i];
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

static void write_tile_info(const AV1_COMMON *const cm,
                            struct aom_write_bit_buffer *wb) {
#if CONFIG_EXT_TILE
  const int tile_width =
      ALIGN_POWER_OF_TWO(cm->tile_width, cm->mib_size_log2) >>
      cm->mib_size_log2;
  const int tile_height =
      ALIGN_POWER_OF_TWO(cm->tile_height, cm->mib_size_log2) >>
      cm->mib_size_log2;

  assert(tile_width > 0);
  assert(tile_height > 0);

// Write the tile sizes
#if CONFIG_EXT_PARTITION
  if (cm->sb_size == BLOCK_128X128) {
    assert(tile_width <= 32);
    assert(tile_height <= 32);
    aom_wb_write_literal(wb, tile_width - 1, 5);
    aom_wb_write_literal(wb, tile_height - 1, 5);
  } else
#endif  // CONFIG_EXT_PARTITION
  {
    assert(tile_width <= 64);
    assert(tile_height <= 64);
    aom_wb_write_literal(wb, tile_width - 1, 6);
    aom_wb_write_literal(wb, tile_height - 1, 6);
  }
#else
  int min_log2_tile_cols, max_log2_tile_cols, ones;
  av1_get_tile_n_bits(cm->mi_cols, &min_log2_tile_cols, &max_log2_tile_cols);

  // columns
  ones = cm->log2_tile_cols - min_log2_tile_cols;
  while (ones--) aom_wb_write_bit(wb, 1);

  if (cm->log2_tile_cols < max_log2_tile_cols) aom_wb_write_bit(wb, 0);

  // rows
  aom_wb_write_bit(wb, cm->log2_tile_rows != 0);
  if (cm->log2_tile_rows != 0) aom_wb_write_bit(wb, cm->log2_tile_rows != 1);
#endif  // CONFIG_EXT_TILE
}

static int get_refresh_mask(AV1_COMP *cpi) {
  int refresh_mask = 0;

#if CONFIG_EXT_REFS
  // NOTE(zoeliu): When LAST_FRAME is to get refreshed, the decoder will be
  // notified to get LAST3_FRAME refreshed and then the virtual indexes for all
  // the 3 LAST reference frames will be updated accordingly, i.e.:
  // (1) The original virtual index for LAST3_FRAME will become the new virtual
  //     index for LAST_FRAME; and
  // (2) The original virtual indexes for LAST_FRAME and LAST2_FRAME will be
  //     shifted and become the new virtual indexes for LAST2_FRAME and
  //     LAST3_FRAME.
  refresh_mask |=
      (cpi->refresh_last_frame << cpi->lst_fb_idxes[LAST_REF_FRAMES - 1]);
  if (cpi->rc.is_bwd_ref_frame && cpi->num_extra_arfs) {
    // We have swapped the virtual indices
    refresh_mask |= (cpi->refresh_bwd_ref_frame << cpi->arf_map[0]);
  } else {
    refresh_mask |= (cpi->refresh_bwd_ref_frame << cpi->bwd_fb_idx);
  }
#else
  refresh_mask |= (cpi->refresh_last_frame << cpi->lst_fb_idx);
#endif  // CONFIG_EXT_REFS

  if (av1_preserve_existing_gf(cpi)) {
    // We have decided to preserve the previously existing golden frame as our
    // new ARF frame. However, in the short term we leave it in the GF slot and,
    // if we're updating the GF with the current decoded frame, we save it
    // instead to the ARF slot.
    // Later, in the function av1_encoder.c:av1_update_reference_frames() we
    // will swap gld_fb_idx and alt_fb_idx to achieve our objective. We do it
    // there so that it can be done outside of the recode loop.
    // Note: This is highly specific to the use of ARF as a forward reference,
    // and this needs to be generalized as other uses are implemented
    // (like RTC/temporal scalability).
    return refresh_mask | (cpi->refresh_golden_frame << cpi->alt_fb_idx);
  } else {
    int arf_idx = cpi->alt_fb_idx;
#if CONFIG_EXT_REFS
    const GF_GROUP *const gf_group = &cpi->twopass.gf_group;
    arf_idx = cpi->arf_map[gf_group->arf_update_idx[gf_group->index]];
#else
    if ((cpi->oxcf.pass == 2) && cpi->multi_arf_allowed) {
      const GF_GROUP *const gf_group = &cpi->twopass.gf_group;
      arf_idx = gf_group->arf_update_idx[gf_group->index];
    }
#endif  // CONFIG_EXT_REFS
    return refresh_mask | (cpi->refresh_golden_frame << cpi->gld_fb_idx) |
           (cpi->refresh_alt_ref_frame << arf_idx);
  }
}

#if CONFIG_EXT_TILE
static INLINE int find_identical_tile(
    const int tile_row, const int tile_col,
    TileBufferEnc (*const tile_buffers)[1024]) {
  const MV32 candidate_offset[1] = { { 1, 0 } };
  const uint8_t *const cur_tile_data =
      tile_buffers[tile_row][tile_col].data + 4;
  const unsigned int cur_tile_size = tile_buffers[tile_row][tile_col].size;

  int i;

  if (tile_row == 0) return 0;

  // (TODO: yunqingwang) For now, only above tile is checked and used.
  // More candidates such as left tile can be added later.
  for (i = 0; i < 1; i++) {
    int row_offset = candidate_offset[0].row;
    int col_offset = candidate_offset[0].col;
    int row = tile_row - row_offset;
    int col = tile_col - col_offset;
    uint8_t tile_hdr;
    const uint8_t *tile_data;
    TileBufferEnc *candidate;

    if (row < 0 || col < 0) continue;

    tile_hdr = *(tile_buffers[row][col].data);

    // Read out tcm bit
    if ((tile_hdr >> 7) == 1) {
      // The candidate is a copy tile itself
      row_offset += tile_hdr & 0x7f;
      row = tile_row - row_offset;
    }

    candidate = &tile_buffers[row][col];

    if (row_offset >= 128 || candidate->size != cur_tile_size) continue;

    tile_data = candidate->data + 4;

    if (memcmp(tile_data, cur_tile_data, cur_tile_size) != 0) continue;

    // Identical tile found
    assert(row_offset > 0);
    return row_offset;
  }

  // No identical tile found
  return 0;
}
#endif  // CONFIG_EXT_TILE

#if CONFIG_TILE_GROUPS
static uint32_t write_tiles(AV1_COMP *const cpi,
                            struct aom_write_bit_buffer *wb,
                            unsigned int *max_tile_size,
                            unsigned int *max_tile_col_size) {
#else
static uint32_t write_tiles(AV1_COMP *const cpi, uint8_t *const dst,
                            unsigned int *max_tile_size,
                            unsigned int *max_tile_col_size) {
#endif
  const AV1_COMMON *const cm = &cpi->common;
#if CONFIG_ANS
  struct AnsCoder token_ans;
#else
  aom_writer mode_bc;
#endif  // CONFIG_ANS
  int tile_row, tile_col;
  TOKENEXTRA *(*const tok_buffers)[MAX_TILE_COLS] = cpi->tile_tok;
  TileBufferEnc(*const tile_buffers)[MAX_TILE_COLS] = cpi->tile_buffers;
  size_t total_size = 0;
  const int tile_cols = cm->tile_cols;
  const int tile_rows = cm->tile_rows;
#if CONFIG_TILE_GROUPS
  const int n_log2_tiles = cm->log2_tile_rows + cm->log2_tile_cols;
  const int have_tiles = n_log2_tiles > 0;
  size_t comp_hdr_size;
  // Fixed size tile groups for the moment
  const int num_tg_hdrs = cm->num_tg;
  const int tg_size = (tile_rows * tile_cols + num_tg_hdrs - 1) / num_tg_hdrs;
  int tile_count = 0;
  int uncompressed_hdr_size = 0;
  uint8_t *dst = NULL;
  struct aom_write_bit_buffer comp_hdr_len_wb;
  struct aom_write_bit_buffer tg_params_wb;
  int saved_offset;
#endif
#if CONFIG_EXT_TILE
  const int have_tiles = tile_cols * tile_rows > 1;
#endif  // CONFIG_EXT_TILE
#if CONFIG_ANS
  struct BufAnsCoder *buf_ans = &cpi->buf_ans;
#endif  // CONFIG_ANS

  *max_tile_size = 0;
  *max_tile_col_size = 0;

// All tile size fields are output on 4 bytes. A call to remux_tiles will
// later compact the data if smaller headers are adequate.

#if CONFIG_EXT_TILE
  for (tile_col = 0; tile_col < tile_cols; tile_col++) {
    TileInfo tile_info;
    const int is_last_col = (tile_col == tile_cols - 1);
    const size_t col_offset = total_size;

    av1_tile_set_col(&tile_info, cm, tile_col);

    // The last column does not have a column header
    if (!is_last_col) total_size += 4;

    for (tile_row = 0; tile_row < tile_rows; tile_row++) {
      TileBufferEnc *const buf = &tile_buffers[tile_row][tile_col];
      unsigned int tile_size;
      const TOKENEXTRA *tok = tok_buffers[tile_row][tile_col];
      const TOKENEXTRA *tok_end = tok + cpi->tok_count[tile_row][tile_col];
      const int data_offset = have_tiles ? 4 : 0;

      av1_tile_set_row(&tile_info, cm, tile_row);

      buf->data = dst + total_size;

      // Is CONFIG_EXT_TILE = 1, every tile in the row has a header,
      // even for the last one, unless no tiling is used at all.
      total_size += data_offset;
#if !CONFIG_ANS
      aom_start_encode(&mode_bc, buf->data + data_offset);
      write_modes(cpi, &tile_info, &mode_bc, &tok, tok_end);
      assert(tok == tok_end);
      aom_stop_encode(&mode_bc);
      tile_size = mode_bc.pos;
#else
      buf_ans_write_reset(buf_ans);
      write_modes(cpi, &tile_info, buf_ans, &tok, tok_end);
      assert(tok == tok_end);
      ans_write_init(&token_ans, buf->data + data_offset);
      buf_ans_flush(buf_ans, &token_ans);
      tile_size = ans_write_end(&token_ans);
#endif  // !CONFIG_ANS

      buf->size = tile_size;

      // Record the maximum tile size we see, so we can compact headers later.
      *max_tile_size = AOMMAX(*max_tile_size, tile_size);

      if (have_tiles) {
        // tile header: size of this tile, or copy offset
        uint32_t tile_header = tile_size;

        // Check if this tile is a copy tile.
        // Very low chances to have copy tiles on the key frames, so don't
        // search on key frames to reduce unnecessary search.
        if (cm->frame_type != KEY_FRAME) {
          const int idendical_tile_offset =
              find_identical_tile(tile_row, tile_col, tile_buffers);

          if (idendical_tile_offset > 0) {
            tile_size = 0;
            tile_header = idendical_tile_offset | 0x80;
            tile_header <<= 24;
          }
        }

        mem_put_le32(buf->data, tile_header);
      }

      total_size += tile_size;
    }

    if (!is_last_col) {
      size_t col_size = total_size - col_offset - 4;
      mem_put_le32(dst + col_offset, col_size);

      // If it is not final packing, record the maximum tile column size we see,
      // otherwise, check if the tile size is out of the range.
      *max_tile_col_size = AOMMAX(*max_tile_col_size, col_size);
    }
  }
#else
#if CONFIG_TILE_GROUPS
  write_uncompressed_header(cpi, wb);

  // Write the tile length code. Use full 32 bit length fields for the moment
  if (have_tiles) aom_wb_write_literal(wb, 3, 2);

  /* Write a placeholder for the number of tiles in each tile group */
  tg_params_wb = *wb;
  saved_offset = wb->bit_offset;
  if (n_log2_tiles) aom_wb_write_literal(wb, 0, n_log2_tiles * 2);

  /* Write a placeholder for the compressed header length */
  comp_hdr_len_wb = *wb;
  aom_wb_write_literal(wb, 0, 16);

  uncompressed_hdr_size = aom_wb_bytes_written(wb);
  dst = wb->bit_buffer;
  comp_hdr_size = write_compressed_header(cpi, dst + uncompressed_hdr_size);
  aom_wb_write_literal(&comp_hdr_len_wb, (int)(comp_hdr_size), 16);
  total_size += uncompressed_hdr_size + comp_hdr_size;
#endif

  for (tile_row = 0; tile_row < tile_rows; tile_row++) {
    TileInfo tile_info;
    const int is_last_row = (tile_row == tile_rows - 1);

    av1_tile_set_row(&tile_info, cm, tile_row);

    for (tile_col = 0; tile_col < tile_cols; tile_col++) {
      const int tile_idx = tile_row * tile_cols + tile_col;
      TileBufferEnc *const buf = &tile_buffers[tile_row][tile_col];
      const int is_last_col = (tile_col == tile_cols - 1);
      unsigned int tile_size;
#if CONFIG_PVQ
      TileDataEnc *this_tile = &cpi->tile_data[tile_idx];
#endif
      const TOKENEXTRA *tok = tok_buffers[tile_row][tile_col];
      const TOKENEXTRA *tok_end = tok + cpi->tok_count[tile_row][tile_col];
#if !CONFIG_TILE_GROUPS
      const int is_last_tile = is_last_col && is_last_row;
      (void)tile_idx;
#else
      // All tiles in a tile group have a length
      const int is_last_tile = 0;
      if (tile_count >= tg_size) {
        // Copy uncompressed header
        memcpy(dst + total_size, dst, uncompressed_hdr_size * sizeof(uint8_t));
        // Write the number of tiles in the group into the last uncompressed
        // header
        aom_wb_write_literal(&tg_params_wb, tile_idx - tile_count,
                             n_log2_tiles);
        aom_wb_write_literal(&tg_params_wb, tile_count - 1, n_log2_tiles);
        tg_params_wb.bit_offset = saved_offset + 8 * total_size;
        // Copy compressed header
        memcpy(dst + total_size + uncompressed_hdr_size,
               dst + uncompressed_hdr_size, comp_hdr_size * sizeof(uint8_t));
        total_size += uncompressed_hdr_size;
        total_size += comp_hdr_size;
        tile_count = 0;
      }
      tile_count++;
#endif
      av1_tile_set_col(&tile_info, cm, tile_col);

      buf->data = dst + total_size;

      // The last tile does not have a header.
      if (!is_last_tile) total_size += 4;

#if CONFIG_ANS
      buf_ans_write_reset(buf_ans);
      write_modes(cpi, &tile_info, buf_ans, &tok, tok_end);
      assert(tok == tok_end);
      ans_write_init(&token_ans, dst + total_size);
      buf_ans_flush(buf_ans, &token_ans);
      tile_size = ans_write_end(&token_ans);
#else
      aom_start_encode(&mode_bc, dst + total_size);
#if CONFIG_PVQ
      // NOTE: This will not work with CONFIG_ANS turned on.
      od_adapt_ctx_reset(&cpi->td.mb.daala_enc.state.adapt, 0);
      cpi->td.mb.pvq_q = &this_tile->pvq_q;
#endif
      write_modes(cpi, &tile_info, &mode_bc, &tok, tok_end);
      assert(tok == tok_end);
      aom_stop_encode(&mode_bc);
      tile_size = mode_bc.pos;
#endif
#if CONFIG_PVQ
      cpi->td.mb.pvq_q = NULL;
#endif  // !CONFIG_ANS

      assert(tile_size > 0);

      buf->size = tile_size;

      if (!is_last_tile) {
        *max_tile_size = AOMMAX(*max_tile_size, tile_size);
        // size of this tile
        mem_put_le32(buf->data, tile_size);
      }

      total_size += tile_size;
    }
  }
#if CONFIG_TILE_GROUPS
  // Write the final tile group size
  if (n_log2_tiles) {
    aom_wb_write_literal(&tg_params_wb, (1 << n_log2_tiles) - tile_count,
                         n_log2_tiles);
    aom_wb_write_literal(&tg_params_wb, tile_count - 1, n_log2_tiles);
  }
#endif
#endif  // CONFIG_EXT_TILE
  return (uint32_t)total_size;
}

static void write_render_size(const AV1_COMMON *cm,
                              struct aom_write_bit_buffer *wb) {
  const int scaling_active =
      cm->width != cm->render_width || cm->height != cm->render_height;
  aom_wb_write_bit(wb, scaling_active);
  if (scaling_active) {
    aom_wb_write_literal(wb, cm->render_width - 1, 16);
    aom_wb_write_literal(wb, cm->render_height - 1, 16);
  }
}

static void write_frame_size(const AV1_COMMON *cm,
                             struct aom_write_bit_buffer *wb) {
  aom_wb_write_literal(wb, cm->width - 1, 16);
  aom_wb_write_literal(wb, cm->height - 1, 16);

  write_render_size(cm, wb);
}

static void write_frame_size_with_refs(AV1_COMP *cpi,
                                       struct aom_write_bit_buffer *wb) {
  AV1_COMMON *const cm = &cpi->common;
  int found = 0;

  MV_REFERENCE_FRAME ref_frame;
  for (ref_frame = LAST_FRAME; ref_frame <= ALTREF_FRAME; ++ref_frame) {
    YV12_BUFFER_CONFIG *cfg = get_ref_frame_buffer(cpi, ref_frame);

    if (cfg != NULL) {
      found =
          cm->width == cfg->y_crop_width && cm->height == cfg->y_crop_height;
      found &= cm->render_width == cfg->render_width &&
               cm->render_height == cfg->render_height;
    }
    aom_wb_write_bit(wb, found);
    if (found) {
      break;
    }
  }

  if (!found) {
    aom_wb_write_literal(wb, cm->width - 1, 16);
    aom_wb_write_literal(wb, cm->height - 1, 16);
    write_render_size(cm, wb);
  }
}

static void write_sync_code(struct aom_write_bit_buffer *wb) {
  aom_wb_write_literal(wb, AV1_SYNC_CODE_0, 8);
  aom_wb_write_literal(wb, AV1_SYNC_CODE_1, 8);
  aom_wb_write_literal(wb, AV1_SYNC_CODE_2, 8);
}

static void write_profile(BITSTREAM_PROFILE profile,
                          struct aom_write_bit_buffer *wb) {
  switch (profile) {
    case PROFILE_0: aom_wb_write_literal(wb, 0, 2); break;
    case PROFILE_1: aom_wb_write_literal(wb, 2, 2); break;
    case PROFILE_2: aom_wb_write_literal(wb, 1, 2); break;
    case PROFILE_3: aom_wb_write_literal(wb, 6, 3); break;
    default: assert(0);
  }
}

static void write_bitdepth_colorspace_sampling(
    AV1_COMMON *const cm, struct aom_write_bit_buffer *wb) {
  if (cm->profile >= PROFILE_2) {
    assert(cm->bit_depth > AOM_BITS_8);
    aom_wb_write_bit(wb, cm->bit_depth == AOM_BITS_10 ? 0 : 1);
  }
  aom_wb_write_literal(wb, cm->color_space, 3);
  if (cm->color_space != AOM_CS_SRGB) {
    // 0: [16, 235] (i.e. xvYCC), 1: [0, 255]
    aom_wb_write_bit(wb, cm->color_range);
    if (cm->profile == PROFILE_1 || cm->profile == PROFILE_3) {
      assert(cm->subsampling_x != 1 || cm->subsampling_y != 1);
      aom_wb_write_bit(wb, cm->subsampling_x);
      aom_wb_write_bit(wb, cm->subsampling_y);
      aom_wb_write_bit(wb, 0);  // unused
    } else {
      assert(cm->subsampling_x == 1 && cm->subsampling_y == 1);
    }
  } else {
    assert(cm->profile == PROFILE_1 || cm->profile == PROFILE_3);
    aom_wb_write_bit(wb, 0);  // unused
  }
}

static void write_uncompressed_header(AV1_COMP *cpi,
                                      struct aom_write_bit_buffer *wb) {
  AV1_COMMON *const cm = &cpi->common;
  MACROBLOCKD *const xd = &cpi->td.mb.e_mbd;

  aom_wb_write_literal(wb, AOM_FRAME_MARKER, 2);

  write_profile(cm->profile, wb);

#if CONFIG_EXT_REFS
  // NOTE: By default all coded frames to be used as a reference
  cm->is_reference_frame = 1;

  if (cm->show_existing_frame) {
    RefCntBuffer *const frame_bufs = cm->buffer_pool->frame_bufs;
    const int frame_to_show = cm->ref_frame_map[cpi->existing_fb_idx_to_show];

    if (frame_to_show < 0 || frame_bufs[frame_to_show].ref_count < 1) {
      aom_internal_error(&cm->error, AOM_CODEC_UNSUP_BITSTREAM,
                         "Buffer %d does not contain a reconstructed frame",
                         frame_to_show);
    }
    ref_cnt_fb(frame_bufs, &cm->new_fb_idx, frame_to_show);

    aom_wb_write_bit(wb, 1);  // show_existing_frame
    aom_wb_write_literal(wb, cpi->existing_fb_idx_to_show, 3);

    return;
  } else {
#endif                        // CONFIG_EXT_REFS
    aom_wb_write_bit(wb, 0);  // show_existing_frame
#if CONFIG_EXT_REFS
  }
#endif  // CONFIG_EXT_REFS

  aom_wb_write_bit(wb, cm->frame_type);
  aom_wb_write_bit(wb, cm->show_frame);
  aom_wb_write_bit(wb, cm->error_resilient_mode);

  if (cm->frame_type == KEY_FRAME) {
    write_sync_code(wb);
    write_bitdepth_colorspace_sampling(cm, wb);
    write_frame_size(cm, wb);
#if CONFIG_PALETTE
    aom_wb_write_bit(wb, cm->allow_screen_content_tools);
#endif  // CONFIG_PALETTE
  } else {
    if (!cm->show_frame) aom_wb_write_bit(wb, cm->intra_only);
#if CONFIG_PALETTE
    if (cm->intra_only) aom_wb_write_bit(wb, cm->allow_screen_content_tools);
#endif  // CONFIG_PALETTE
    if (!cm->error_resilient_mode) {
      if (cm->intra_only) {
        aom_wb_write_bit(wb,
                         cm->reset_frame_context == RESET_FRAME_CONTEXT_ALL);
      } else {
        aom_wb_write_bit(wb,
                         cm->reset_frame_context != RESET_FRAME_CONTEXT_NONE);
        if (cm->reset_frame_context != RESET_FRAME_CONTEXT_NONE)
          aom_wb_write_bit(wb,
                           cm->reset_frame_context == RESET_FRAME_CONTEXT_ALL);
      }
    }

#if CONFIG_EXT_REFS
    cpi->refresh_frame_mask = get_refresh_mask(cpi);
#endif  // CONFIG_EXT_REFS

    if (cm->intra_only) {
      write_sync_code(wb);
      write_bitdepth_colorspace_sampling(cm, wb);

#if CONFIG_EXT_REFS
      aom_wb_write_literal(wb, cpi->refresh_frame_mask, REF_FRAMES);
#else
      aom_wb_write_literal(wb, get_refresh_mask(cpi), REF_FRAMES);
#endif  // CONFIG_EXT_REFS
      write_frame_size(cm, wb);
    } else {
      MV_REFERENCE_FRAME ref_frame;

#if CONFIG_EXT_REFS
      aom_wb_write_literal(wb, cpi->refresh_frame_mask, REF_FRAMES);
#else
      aom_wb_write_literal(wb, get_refresh_mask(cpi), REF_FRAMES);
#endif  // CONFIG_EXT_REFS

#if CONFIG_EXT_REFS
      if (!cpi->refresh_frame_mask) {
        // NOTE: "cpi->refresh_frame_mask == 0" indicates that the coded frame
        //       will not be used as a reference
        cm->is_reference_frame = 0;
      }
#endif  // CONFIG_EXT_REFS

      for (ref_frame = LAST_FRAME; ref_frame <= ALTREF_FRAME; ++ref_frame) {
        assert(get_ref_frame_map_idx(cpi, ref_frame) != INVALID_IDX);
        aom_wb_write_literal(wb, get_ref_frame_map_idx(cpi, ref_frame),
                             REF_FRAMES_LOG2);
        aom_wb_write_bit(wb, cm->ref_frame_sign_bias[ref_frame]);
      }

#if CONFIG_FRAME_SIZE
      if (cm->error_resilient_mode == 0) {
        write_frame_size_with_refs(cpi, wb);
      } else {
        write_frame_size(cm, wb);
      }
#else
      write_frame_size_with_refs(cpi, wb);
#endif

      aom_wb_write_bit(wb, cm->allow_high_precision_mv);

      fix_interp_filter(cm, cpi->td.counts);
      write_interp_filter(cm->interp_filter, wb);
    }
  }

  if (!cm->error_resilient_mode) {
    aom_wb_write_bit(
        wb, cm->refresh_frame_context == REFRESH_FRAME_CONTEXT_FORWARD);
  }

  aom_wb_write_literal(wb, cm->frame_context_idx, FRAME_CONTEXTS_LOG2);

  assert(cm->mib_size == num_8x8_blocks_wide_lookup[cm->sb_size]);
  assert(cm->mib_size == 1 << cm->mib_size_log2);
#if CONFIG_EXT_PARTITION
  assert(cm->sb_size == BLOCK_128X128 || cm->sb_size == BLOCK_64X64);
  aom_wb_write_bit(wb, cm->sb_size == BLOCK_128X128 ? 1 : 0);
#else
  assert(cm->sb_size == BLOCK_64X64);
#endif  // CONFIG_EXT_PARTITION

  encode_loopfilter(cm, wb);
#if CONFIG_DERING
  encode_dering(cm->dering_level, wb);
#endif  // CONFIG_DERING
#if CONFIG_CLPF
  encode_clpf(cm, wb);
#endif
#if CONFIG_LOOP_RESTORATION
  encode_restoration_mode(cm, wb);
#endif  // CONFIG_LOOP_RESTORATION
  encode_quantization(cm, wb);
  encode_segmentation(cm, xd, wb);
#if CONFIG_DELTA_Q
  {
    int i;
    struct segmentation *const seg = &cm->seg;
    int segment_quantizer_active = 0;
    for (i = 0; i < MAX_SEGMENTS; i++) {
      if (segfeature_active(seg, i, SEG_LVL_ALT_Q)) {
        segment_quantizer_active = 1;
      }
    }
    if (segment_quantizer_active == 0) {
      cm->delta_q_present_flag = cpi->oxcf.aq_mode == DELTA_AQ;
      aom_wb_write_bit(wb, cm->delta_q_present_flag);
      if (cm->delta_q_present_flag) {
        aom_wb_write_literal(wb, OD_ILOG_NZ(cm->delta_q_res) - 1, 2);
        xd->prev_qindex = cm->base_qindex;
      }
    }
  }
#endif

  if (!cm->seg.enabled && xd->lossless[0])
    cm->tx_mode = ONLY_4X4;
  else
    write_txfm_mode(cm->tx_mode, wb);

  if (cpi->allow_comp_inter_inter) {
    const int use_hybrid_pred = cm->reference_mode == REFERENCE_MODE_SELECT;
    const int use_compound_pred = cm->reference_mode != SINGLE_REFERENCE;

    aom_wb_write_bit(wb, use_hybrid_pred);
    if (!use_hybrid_pred) aom_wb_write_bit(wb, use_compound_pred);
  }

  write_tile_info(cm, wb);
}

#if CONFIG_GLOBAL_MOTION
static void write_global_motion_params(Global_Motion_Params *params,
                                       aom_prob *probs, aom_writer *w) {
  GLOBAL_MOTION_TYPE gmtype = get_gmtype(params);
  av1_write_token(w, av1_global_motion_types_tree, probs,
                  &global_motion_types_encodings[gmtype]);
  switch (gmtype) {
    case GLOBAL_ZERO: break;
    case GLOBAL_AFFINE:
      aom_write_primitive_symmetric(
          w, (params->motion_params.wmmat[4] >> GM_ALPHA_PREC_DIFF),
          GM_ABS_ALPHA_BITS);
      aom_write_primitive_symmetric(
          w, (params->motion_params.wmmat[5] >> GM_ALPHA_PREC_DIFF) -
                 (1 << GM_ALPHA_PREC_BITS),
          GM_ABS_ALPHA_BITS);
    // fallthrough intended
    case GLOBAL_ROTZOOM:
      aom_write_primitive_symmetric(
          w, (params->motion_params.wmmat[2] >> GM_ALPHA_PREC_DIFF),
          GM_ABS_ALPHA_BITS);
      aom_write_primitive_symmetric(
          w, (params->motion_params.wmmat[3] >> GM_ALPHA_PREC_DIFF) -
                 (1 << GM_ALPHA_PREC_BITS),
          GM_ABS_ALPHA_BITS);
    // fallthrough intended
    case GLOBAL_TRANSLATION:
      aom_write_primitive_symmetric(
          w, (params->motion_params.wmmat[0] >> GM_TRANS_PREC_DIFF),
          GM_ABS_TRANS_BITS);
      aom_write_primitive_symmetric(
          w, (params->motion_params.wmmat[1] >> GM_TRANS_PREC_DIFF),
          GM_ABS_TRANS_BITS);
      break;
    default: assert(0);
  }
}

static void write_global_motion(AV1_COMP *cpi, aom_writer *w) {
  AV1_COMMON *const cm = &cpi->common;
  int frame;
  for (frame = LAST_FRAME; frame <= ALTREF_FRAME; ++frame) {
    if (!cpi->global_motion_used[frame]) {
      memset(&cm->global_motion[frame], 0, sizeof(*cm->global_motion));
    }
    write_global_motion_params(&cm->global_motion[frame],
                               cm->fc->global_motion_types_prob, w);
    /*
          printf("Enc Ref %d [%d] (used %d): %d %d %d %d\n",
                 frame, cm->current_video_frame, cpi->global_motion_used[frame],
                 cm->global_motion[frame].motion_params.wmmat[0].as_mv.row,
                 cm->global_motion[frame].motion_params.wmmat[0].as_mv.col,
                 cm->global_motion[frame].motion_params.wmmat[1].as_mv.row,
                 cm->global_motion[frame].motion_params.wmmat[1].as_mv.col);
    */
  }
}
#endif

static uint32_t write_compressed_header(AV1_COMP *cpi, uint8_t *data) {
  AV1_COMMON *const cm = &cpi->common;
#if CONFIG_SUPERTX
  MACROBLOCKD *const xd = &cpi->td.mb.e_mbd;
#endif  // CONFIG_SUPERTX
  FRAME_CONTEXT *const fc = cm->fc;
  FRAME_COUNTS *counts = cpi->td.counts;
  aom_writer *header_bc;
  int i, j;

#if CONFIG_TILE_GROUPS
  const int probwt = cm->num_tg;
#else
  const int probwt = 1;
#endif

#if CONFIG_ANS
  struct AnsCoder header_ans;
  int header_size;
  header_bc = &cpi->buf_ans;
  buf_ans_write_reset(header_bc);
#else
  aom_writer real_header_bc;
  header_bc = &real_header_bc;
  aom_start_encode(header_bc, data);
#endif

#if CONFIG_LOOP_RESTORATION
  encode_restoration(cm, header_bc);
#endif  // CONFIG_LOOP_RESTORATION
  update_txfm_probs(cm, header_bc, counts);
#if !CONFIG_PVQ
  update_coef_probs(cpi, header_bc);
#endif
#if CONFIG_VAR_TX
  update_txfm_partition_probs(cm, header_bc, counts, probwt);
#if CONFIG_EXT_TX && CONFIG_RECT_TX
  if (cm->tx_mode == TX_MODE_SELECT) {
    for (i = 1; i < TX_SIZES - 1; ++i)
      av1_cond_prob_diff_update(header_bc, &fc->rect_tx_prob[i],
                                counts->rect_tx[i], probwt);
  }
#endif  // CONFIG_EXT_TX && CONFIG_RECT_TX
#endif

  update_skip_probs(cm, header_bc, counts);
#if CONFIG_DELTA_Q
  update_delta_q_probs(cm, header_bc, counts);
#endif
#if !CONFIG_EC_ADAPT
  update_seg_probs(cpi, header_bc);

  for (i = 0; i < INTRA_MODES; ++i) {
    prob_diff_update(av1_intra_mode_tree, fc->uv_mode_prob[i],
                     counts->uv_mode[i], INTRA_MODES, probwt, header_bc);
  }

#if CONFIG_EXT_PARTITION_TYPES
  prob_diff_update(av1_partition_tree, fc->partition_prob[0],
                   counts->partition[0], PARTITION_TYPES, probwt, header_bc);
  for (i = 1; i < PARTITION_CONTEXTS; ++i)
    prob_diff_update(av1_ext_partition_tree, fc->partition_prob[i],
                     counts->partition[i], EXT_PARTITION_TYPES, probwt,
                     header_bc);
#else
  for (i = 0; i < PARTITION_CONTEXTS; ++i) {
    prob_diff_update(av1_partition_tree, fc->partition_prob[i],
                     counts->partition[i], PARTITION_TYPES, probwt, header_bc);
  }
#endif  // CONFIG_EC_ADAPT, CONFIG_DAALA_EC

#if CONFIG_EXT_INTRA
  for (i = 0; i < INTRA_FILTERS + 1; ++i)
    prob_diff_update(av1_intra_filter_tree, fc->intra_filter_probs[i],
                     counts->intra_filter[i], INTRA_FILTERS, probwt, header_bc);
#endif  // CONFIG_EXT_INTRA
#endif  // CONFIG_EC_ADAPT, CONFIG_DAALA_EC
  if (frame_is_intra_only(cm)) {
    av1_copy(cm->kf_y_prob, av1_kf_y_mode_prob);
#if CONFIG_DAALA_EC
    av1_copy(cm->kf_y_cdf, av1_kf_y_mode_cdf);
#endif

#if !CONFIG_EC_ADAPT
    for (i = 0; i < INTRA_MODES; ++i)
      for (j = 0; j < INTRA_MODES; ++j)
        prob_diff_update(av1_intra_mode_tree, cm->kf_y_prob[i][j],
                         counts->kf_y_mode[i][j], INTRA_MODES, probwt,
                         header_bc);
#endif  // CONFIG_EC_ADAPT
  } else {
#if CONFIG_REF_MV
    update_inter_mode_probs(cm, header_bc, counts);
#else
#if !CONFIG_EC_ADAPT
    for (i = 0; i < INTER_MODE_CONTEXTS; ++i) {
      prob_diff_update(av1_inter_mode_tree, cm->fc->inter_mode_probs[i],
                       counts->inter_mode[i], INTER_MODES, probwt, header_bc);
    }
#endif
#endif
#if CONFIG_EXT_INTER
    update_inter_compound_mode_probs(cm, probwt, header_bc);

    if (cm->reference_mode != COMPOUND_REFERENCE) {
      for (i = 0; i < BLOCK_SIZE_GROUPS; i++) {
        if (is_interintra_allowed_bsize_group(i)) {
          av1_cond_prob_diff_update(header_bc, &fc->interintra_prob[i],
                                    cm->counts.interintra[i], probwt);
        }
      }
      for (i = 0; i < BLOCK_SIZE_GROUPS; i++) {
        prob_diff_update(
            av1_interintra_mode_tree, cm->fc->interintra_mode_prob[i],
            counts->interintra_mode[i], INTERINTRA_MODES, probwt, header_bc);
      }
      for (i = 0; i < BLOCK_SIZES; i++) {
        if (is_interintra_allowed_bsize(i) && is_interintra_wedge_used(i))
          av1_cond_prob_diff_update(header_bc, &fc->wedge_interintra_prob[i],
                                    cm->counts.wedge_interintra[i], probwt);
      }
    }
    if (cm->reference_mode != SINGLE_REFERENCE) {
      for (i = 0; i < BLOCK_SIZES; i++)
        if (is_interinter_wedge_used(i))
          av1_cond_prob_diff_update(header_bc, &fc->wedge_interinter_prob[i],
                                    cm->counts.wedge_interinter[i], probwt);
    }
#endif  // CONFIG_EXT_INTER

#if CONFIG_MOTION_VAR || CONFIG_WARPED_MOTION
    for (i = BLOCK_8X8; i < BLOCK_SIZES; ++i)
      prob_diff_update(av1_motion_mode_tree, fc->motion_mode_prob[i],
                       counts->motion_mode[i], MOTION_MODES, probwt, header_bc);
#endif  // CONFIG_MOTION_VAR || CONFIG_WARPED_MOTION
#if !CONFIG_EC_ADAPT
    if (cm->interp_filter == SWITCHABLE)
      update_switchable_interp_probs(cm, header_bc, counts);
#endif

    for (i = 0; i < INTRA_INTER_CONTEXTS; i++)
      av1_cond_prob_diff_update(header_bc, &fc->intra_inter_prob[i],
                                counts->intra_inter[i], probwt);

    if (cpi->allow_comp_inter_inter) {
      const int use_hybrid_pred = cm->reference_mode == REFERENCE_MODE_SELECT;
      if (use_hybrid_pred)
        for (i = 0; i < COMP_INTER_CONTEXTS; i++)
          av1_cond_prob_diff_update(header_bc, &fc->comp_inter_prob[i],
                                    counts->comp_inter[i], probwt);
    }

    if (cm->reference_mode != COMPOUND_REFERENCE) {
      for (i = 0; i < REF_CONTEXTS; i++) {
        for (j = 0; j < (SINGLE_REFS - 1); j++) {
          av1_cond_prob_diff_update(header_bc, &fc->single_ref_prob[i][j],
                                    counts->single_ref[i][j], probwt);
        }
      }
    }
    if (cm->reference_mode != SINGLE_REFERENCE) {
      for (i = 0; i < REF_CONTEXTS; i++) {
#if CONFIG_EXT_REFS
        for (j = 0; j < (FWD_REFS - 1); j++) {
          av1_cond_prob_diff_update(header_bc, &fc->comp_ref_prob[i][j],
                                    counts->comp_ref[i][j], probwt);
        }
        for (j = 0; j < (BWD_REFS - 1); j++) {
          av1_cond_prob_diff_update(header_bc, &fc->comp_bwdref_prob[i][j],
                                    counts->comp_bwdref[i][j], probwt);
        }
#else
        for (j = 0; j < (COMP_REFS - 1); j++) {
          av1_cond_prob_diff_update(header_bc, &fc->comp_ref_prob[i][j],
                                    counts->comp_ref[i][j], probwt);
        }
#endif  // CONFIG_EXT_REFS
      }
    }

#if !CONFIG_EC_ADAPT
    for (i = 0; i < BLOCK_SIZE_GROUPS; ++i) {
      prob_diff_update(av1_intra_mode_tree, cm->fc->y_mode_prob[i],
                       counts->y_mode[i], INTRA_MODES, probwt, header_bc);
    }
#endif

    av1_write_nmv_probs(cm, cm->allow_high_precision_mv, header_bc,
#if CONFIG_REF_MV
                        counts->mv);
#else
                        &counts->mv);
#endif
#if CONFIG_DAALA_EC
    av1_tree_to_cdf(av1_mv_joint_tree, cm->fc->nmvc.joints,
                    cm->fc->nmvc.joint_cdf);
#endif
#if !CONFIG_EC_ADAPT
    update_ext_tx_probs(cm, header_bc);
#endif
#if CONFIG_SUPERTX
    if (!xd->lossless[0]) update_supertx_probs(cm, probwt, header_bc);
#endif  // CONFIG_SUPERTX
#if CONFIG_GLOBAL_MOTION
    write_global_motion(cpi, header_bc);
#endif  // CONFIG_GLOBAL_MOTION
  }
#if CONFIG_EC_MULTISYMBOL
  av1_coef_pareto_cdfs(fc);
  av1_set_mv_cdfs(&fc->nmvc);
#if CONFIG_DAALA_EC
  av1_set_mode_cdfs(cm);
#endif
#endif
#if CONFIG_ANS
  ans_write_init(&header_ans, data);
  buf_ans_flush(header_bc, &header_ans);
  header_size = ans_write_end(&header_ans);
  assert(header_size <= 0xffff);
  return header_size;
#else
  aom_stop_encode(header_bc);
  assert(header_bc->pos <= 0xffff);
  return header_bc->pos;
#endif  // CONFIG_ANS
}

#if !CONFIG_TILE_GROUPS
static int choose_size_bytes(uint32_t size, int spare_msbs) {
  // Choose the number of bytes required to represent size, without
  // using the 'spare_msbs' number of most significant bits.

  // Make sure we will fit in 4 bytes to start with..
  if (spare_msbs > 0 && size >> (32 - spare_msbs) != 0) return -1;

  // Normalise to 32 bits
  size <<= spare_msbs;

  if (size >> 24 != 0)
    return 4;
  else if (size >> 16 != 0)
    return 3;
  else if (size >> 8 != 0)
    return 2;
  else
    return 1;
}

static void mem_put_varsize(uint8_t *const dst, const int sz, const int val) {
  switch (sz) {
    case 1: dst[0] = (uint8_t)(val & 0xff); break;
    case 2: mem_put_le16(dst, val); break;
    case 3: mem_put_le24(dst, val); break;
    case 4: mem_put_le32(dst, val); break;
    default: assert("Invalid size" && 0); break;
  }
}
static int remux_tiles(const AV1_COMMON *const cm, uint8_t *dst,
                       const uint32_t data_size, const uint32_t max_tile_size,
                       const uint32_t max_tile_col_size,
                       int *const tile_size_bytes,
                       int *const tile_col_size_bytes) {
// Choose the tile size bytes (tsb) and tile column size bytes (tcsb)
#if CONFIG_EXT_TILE
  // The top bit in the tile size field indicates tile copy mode, so we
  // have 1 less bit to code the tile size
  const int tsb = choose_size_bytes(max_tile_size, 1);
  const int tcsb = choose_size_bytes(max_tile_col_size, 0);
#else
  const int tsb = choose_size_bytes(max_tile_size, 0);
  const int tcsb = 4;  // This is ignored
  (void)max_tile_col_size;
#endif  // CONFIG_EXT_TILE

  assert(tsb > 0);
  assert(tcsb > 0);

  *tile_size_bytes = tsb;
  *tile_col_size_bytes = tcsb;

  if (tsb == 4 && tcsb == 4) {
    return data_size;
  } else {
    uint32_t wpos = 0;
    uint32_t rpos = 0;

#if CONFIG_EXT_TILE
    int tile_row;
    int tile_col;

    for (tile_col = 0; tile_col < cm->tile_cols; tile_col++) {
      // All but the last column has a column header
      if (tile_col < cm->tile_cols - 1) {
        uint32_t tile_col_size = mem_get_le32(dst + rpos);
        rpos += 4;

        // Adjust the tile column size by the number of bytes removed
        // from the tile size fields.
        tile_col_size -= (4 - tsb) * cm->tile_rows;

        mem_put_varsize(dst + wpos, tcsb, tile_col_size);
        wpos += tcsb;
      }

      for (tile_row = 0; tile_row < cm->tile_rows; tile_row++) {
        // All, including the last row has a header
        uint32_t tile_header = mem_get_le32(dst + rpos);
        rpos += 4;

        // If this is a copy tile, we need to shift the MSB to the
        // top bit of the new width, and there is no data to copy.
        if (tile_header >> 31 != 0) {
          if (tsb < 4) tile_header >>= 32 - 8 * tsb;
          mem_put_varsize(dst + wpos, tsb, tile_header);
          wpos += tsb;
        } else {
          mem_put_varsize(dst + wpos, tsb, tile_header);
          wpos += tsb;

          memmove(dst + wpos, dst + rpos, tile_header);
          rpos += tile_header;
          wpos += tile_header;
        }
      }
    }
#else
    const int n_tiles = cm->tile_cols * cm->tile_rows;
    int n;

    for (n = 0; n < n_tiles; n++) {
      int tile_size;

      if (n == n_tiles - 1) {
        tile_size = data_size - rpos;
      } else {
        tile_size = mem_get_le32(dst + rpos);
        rpos += 4;
        mem_put_varsize(dst + wpos, tsb, tile_size);
        wpos += tsb;
      }

      memmove(dst + wpos, dst + rpos, tile_size);

      rpos += tile_size;
      wpos += tile_size;
    }
#endif  // CONFIG_EXT_TILE

    assert(rpos > wpos);
    assert(rpos == data_size);

    return wpos;
  }
}
#endif  // CONFIG_TILE_GROUPS

void av1_pack_bitstream(AV1_COMP *const cpi, uint8_t *dst, size_t *size) {
  uint8_t *data = dst;
#if !CONFIG_TILE_GROUPS
  uint32_t compressed_header_size;
  uint32_t uncompressed_header_size;
  struct aom_write_bit_buffer saved_wb;
#endif
  uint32_t data_size;
  struct aom_write_bit_buffer wb = { data, 0 };

  unsigned int max_tile_size;
  unsigned int max_tile_col_size;

#if !CONFIG_TILE_GROUPS
  int tile_size_bytes;
  int tile_col_size_bytes;
  AV1_COMMON *const cm = &cpi->common;
  const int have_tiles = cm->tile_cols * cm->tile_rows > 1;

#if CONFIG_BITSTREAM_DEBUG
  bitstream_queue_reset_write();
#endif

  // Write the uncompressed header
  write_uncompressed_header(cpi, &wb);

#if CONFIG_EXT_REFS
  if (cm->show_existing_frame) {
    *size = aom_wb_bytes_written(&wb);
    return;
  }
#endif  // CONFIG_EXT_REFS

  // We do not know these in advance. Output placeholder bit.
  saved_wb = wb;
  // Write tile size magnitudes
  if (have_tiles) {
// Note that the last item in the uncompressed header is the data
// describing tile configuration.
#if CONFIG_EXT_TILE
    // Number of bytes in tile column size - 1
    aom_wb_write_literal(&wb, 0, 2);
#endif  // CONFIG_EXT_TILE
    // Number of bytes in tile size - 1
    aom_wb_write_literal(&wb, 0, 2);
  }
  // Size of compressed header
  aom_wb_write_literal(&wb, 0, 16);

  uncompressed_header_size = (uint32_t)aom_wb_bytes_written(&wb);
  data += uncompressed_header_size;

  aom_clear_system_state();

  // Write the compressed header
  compressed_header_size = write_compressed_header(cpi, data);
  data += compressed_header_size;

  // Write the encoded tile data
  data_size = write_tiles(cpi, data, &max_tile_size, &max_tile_col_size);
#else
  data_size = write_tiles(cpi, &wb, &max_tile_size, &max_tile_col_size);
#endif
#if !CONFIG_TILE_GROUPS
  if (have_tiles) {
    data_size =
        remux_tiles(cm, data, data_size, max_tile_size, max_tile_col_size,
                    &tile_size_bytes, &tile_col_size_bytes);
  }

  data += data_size;

  // Now fill in the gaps in the uncompressed header.
  if (have_tiles) {
#if CONFIG_EXT_TILE
    assert(tile_col_size_bytes >= 1 && tile_col_size_bytes <= 4);
    aom_wb_write_literal(&saved_wb, tile_col_size_bytes - 1, 2);
#endif  // CONFIG_EXT_TILE
    assert(tile_size_bytes >= 1 && tile_size_bytes <= 4);
    aom_wb_write_literal(&saved_wb, tile_size_bytes - 1, 2);
  }
  // TODO(jbb): Figure out what to do if compressed_header_size > 16 bits.
  assert(compressed_header_size <= 0xffff);
  aom_wb_write_literal(&saved_wb, compressed_header_size, 16);
#else
  data += data_size;
#endif
  *size = data - dst;
}
