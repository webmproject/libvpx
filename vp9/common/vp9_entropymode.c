/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include "vpx_mem/vpx_mem.h"

#include "vp9/common/vp9_onyxc_int.h"
#include "vp9/common/vp9_seg_common.h"

const vp9_prob vp9_intra_mode_prob[INTRA_MODES] = {
    227, 223, 219, 213, 204, 191, 170, 127
};

const vp9_prob vp9_intra_predictor_prob[3] = {
    170, 192, 170
};

const vp9_prob vp9_kf_uv_mode_prob[INTRA_MODES][INTRA_MODES - 1] = {
  { 144,  11,  54, 157, 195, 130,  46,  58, 108 },  // y = dc
  { 118,  15, 123, 148, 131, 101,  44,  93, 131 },  // y = v
  { 113,  12,  23, 188, 226, 142,  26,  32, 125 },  // y = h
  { 120,  11,  50, 123, 163, 135,  64,  77, 103 },  // y = d45
  { 113,   9,  36, 155, 111, 157,  32,  44, 161 },  // y = d135
  { 116,   9,  55, 176,  76,  96,  37,  61, 149 },  // y = d117
  { 115,   9,  28, 141, 161, 167,  21,  25, 193 },  // y = d153
  { 120,  12,  32, 145, 195, 142,  32,  38,  86 },  // y = d207
  { 116,  12,  64, 120, 140, 125,  49, 115, 121 },  // y = d63
  { 102,  19,  66, 162, 182, 122,  35,  59, 128 }   // y = tm
};

static const vp9_prob default_if_y_probs[BLOCK_SIZE_GROUPS][INTRA_MODES - 1] = {
  {  65,  32,  18, 144, 162, 194,  41,  51,  98 },  // block_size < 8x8
  { 132,  68,  18, 165, 217, 196,  45,  40,  78 },  // block_size < 16x16
  { 173,  80,  19, 176, 240, 193,  64,  35,  46 },  // block_size < 32x32
  { 221, 135,  38, 194, 248, 121,  96,  85,  29 }   // block_size >= 32x32
};

static const vp9_prob default_if_uv_probs[INTRA_MODES][INTRA_MODES - 1] = {
  { 120,   7,  76, 176, 208, 126,  28,  54, 103 },  // y = dc
  {  48,  12, 154, 155, 139,  90,  34, 117, 119 },  // y = v
  {  67,   6,  25, 204, 243, 158,  13,  21,  96 },  // y = h
  {  97,   5,  44, 131, 176, 139,  48,  68,  97 },  // y = d45
  {  83,   5,  42, 156, 111, 152,  26,  49, 152 },  // y = d135
  {  80,   5,  58, 178,  74,  83,  33,  62, 145 },  // y = d117
  {  86,   5,  32, 154, 192, 168,  14,  22, 163 },  // y = d153
  {  85,   5,  32, 156, 216, 148,  19,  29,  73 },  // y = d207
  {  77,   7,  64, 116, 132, 122,  37, 126, 120 },  // y = d63
  { 101,  21, 107, 181, 192, 103,  19,  67, 125 }   // y = tm
};

const vp9_prob vp9_kf_partition_probs[PARTITION_CONTEXTS]
                                     [PARTITION_TYPES - 1] = {
  // 8x8 -> 4x4
  { 158,  97,  94 },  // a/l both not split
  {  93,  24,  99 },  // a split, l not split
  {  85, 119,  44 },  // l split, a not split
  {  62,  59,  67 },  // a/l both split
  // 16x16 -> 8x8
  { 149,  53,  53 },  // a/l both not split
  {  94,  20,  48 },  // a split, l not split
  {  83,  53,  24 },  // l split, a not split
  {  52,  18,  18 },  // a/l both split
  // 32x32 -> 16x16
  { 150,  40,  39 },  // a/l both not split
  {  78,  12,  26 },  // a split, l not split
  {  67,  33,  11 },  // l split, a not split
  {  24,   7,   5 },  // a/l both split
  // 64x64 -> 32x32
  { 174,  35,  49 },  // a/l both not split
  {  68,  11,  27 },  // a split, l not split
  {  57,  15,   9 },  // l split, a not split
  {  12,   3,   3 },  // a/l both split
};

static const vp9_prob default_partition_probs[PARTITION_CONTEXTS]
                                             [PARTITION_TYPES - 1] = {
  // 8x8 -> 4x4
  { 199, 122, 141 },  // a/l both not split
  { 147,  63, 159 },  // a split, l not split
  { 148, 133, 118 },  // l split, a not split
  { 121, 104, 114 },  // a/l both split
  // 16x16 -> 8x8
  { 174,  73,  87 },  // a/l both not split
  {  92,  41,  83 },  // a split, l not split
  {  82,  99,  50 },  // l split, a not split
  {  53,  39,  39 },  // a/l both split
  // 32x32 -> 16x16
  { 177,  58,  59 },  // a/l both not split
  {  68,  26,  63 },  // a split, l not split
  {  52,  79,  25 },  // l split, a not split
  {  17,  14,  12 },  // a/l both split
  // 64x64 -> 32x32
  { 222,  34,  30 },  // a/l both not split
  {  72,  16,  44 },  // a split, l not split
  {  58,  32,  12 },  // l split, a not split
  {  10,   7,   6 },  // a/l both split
};

static const vp9_prob default_inter_mode_probs[INTER_MODE_CONTEXTS]
                                              [INTER_MODES - 1] = {
  {2,       173,   34},  // 0 = both zero mv
  {7,       145,   85},  // 1 = one zero mv + one a predicted mv
  {7,       166,   63},  // 2 = two predicted mvs
  {7,       94,    66},  // 3 = one predicted/zero and one new mv
  {8,       64,    46},  // 4 = two new mvs
  {17,      81,    31},  // 5 = one intra neighbour + x
  {25,      29,    30},  // 6 = two intra neighbours
};

/* Array indices are identical to previously-existing INTRAMODECONTEXTNODES. */
const vp9_tree_index vp9_intra_mode_tree[TREE_SIZE(INTRA_MODES)] = {
  -DC_PRED, 2,                      /* 0 = DC_NODE */
  -TM_PRED, 4,                      /* 1 = TM_NODE */
  -V_PRED, 6,                       /* 2 = V_NODE */
  8, 12,                            /* 3 = COM_NODE */
  -H_PRED, 10,                      /* 4 = H_NODE */
  -D135_PRED, -D117_PRED,           /* 5 = D135_NODE */
  -D45_PRED, 14,                    /* 6 = D45_NODE */
  -D63_PRED, 16,                    /* 7 = D63_NODE */
  -D153_PRED, -D207_PRED             /* 8 = D153_NODE */
};

const vp9_tree_index vp9_inter_mode_tree[TREE_SIZE(INTER_MODES)] = {
  -INTER_OFFSET(ZEROMV), 2,
  -INTER_OFFSET(NEARESTMV), 4,
  -INTER_OFFSET(NEARMV), -INTER_OFFSET(NEWMV)
};

const vp9_tree_index vp9_partition_tree[TREE_SIZE(PARTITION_TYPES)] = {
  -PARTITION_NONE, 2,
  -PARTITION_HORZ, 4,
  -PARTITION_VERT, -PARTITION_SPLIT
};

static const vp9_prob default_intra_inter_p[INTRA_INTER_CONTEXTS] = {
  9, 102, 187, 225
};

static const vp9_prob default_comp_inter_p[COMP_INTER_CONTEXTS] = {
  239, 183, 119,  96,  41
};

static const vp9_prob default_comp_ref_p[REF_CONTEXTS] = {
  50, 126, 123, 221, 226
};

static const vp9_prob default_single_ref_p[REF_CONTEXTS][2] = {
  {  33,  16 },
  {  77,  74 },
  { 142, 142 },
  { 172, 170 },
  { 238, 247 }
};

static const struct tx_probs default_tx_probs = {
  { { 3, 136, 37 },
    { 5, 52,  13 } },

  { { 20, 152 },
    { 15, 101 } },

  { { 100 },
    { 66  } }
};

void tx_counts_to_branch_counts_32x32(const unsigned int *tx_count_32x32p,
                                      unsigned int (*ct_32x32p)[2]) {
  ct_32x32p[0][0] = tx_count_32x32p[TX_4X4];
  ct_32x32p[0][1] = tx_count_32x32p[TX_8X8] +
                    tx_count_32x32p[TX_16X16] +
                    tx_count_32x32p[TX_32X32];
  ct_32x32p[1][0] = tx_count_32x32p[TX_8X8];
  ct_32x32p[1][1] = tx_count_32x32p[TX_16X16] +
                    tx_count_32x32p[TX_32X32];
  ct_32x32p[2][0] = tx_count_32x32p[TX_16X16];
  ct_32x32p[2][1] = tx_count_32x32p[TX_32X32];
}

void tx_counts_to_branch_counts_16x16(const unsigned int *tx_count_16x16p,
                                      unsigned int (*ct_16x16p)[2]) {
  ct_16x16p[0][0] = tx_count_16x16p[TX_4X4];
  ct_16x16p[0][1] = tx_count_16x16p[TX_8X8] + tx_count_16x16p[TX_16X16];
  ct_16x16p[1][0] = tx_count_16x16p[TX_8X8];
  ct_16x16p[1][1] = tx_count_16x16p[TX_16X16];
}

void tx_counts_to_branch_counts_8x8(const unsigned int *tx_count_8x8p,
                                    unsigned int (*ct_8x8p)[2]) {
  ct_8x8p[0][0] = tx_count_8x8p[TX_4X4];
  ct_8x8p[0][1] = tx_count_8x8p[TX_8X8];
}

static const vp9_prob default_txfm_partition_probs[TXFM_PARTITION_CONTEXTS] = {
    141, 139, 175, 87, 196, 165, 177, 75, 220, 179, 205, 197
};

static const vp9_prob default_skip_probs[SKIP_CONTEXTS] = {
  192, 128, 64
};

static const vp9_prob default_switchable_interp_prob[SWITCHABLE_FILTER_CONTEXTS]
                                                    [SWITCHABLE_FILTERS - 1] = {
  { 235, 162, },
  { 36, 255, },
  { 34, 3, },
  { 149, 144, },
};

void vp9_init_mode_probs(FRAME_CONTEXT *fc) {
  vp9_copy(fc->uv_mode_prob, default_if_uv_probs);
  vp9_copy(fc->y_mode_prob, default_if_y_probs);
  vp9_copy(fc->switchable_interp_prob, default_switchable_interp_prob);
  vp9_copy(fc->partition_prob, default_partition_probs);
  vp9_copy(fc->intra_inter_prob, default_intra_inter_p);
  vp9_copy(fc->comp_inter_prob, default_comp_inter_p);
  vp9_copy(fc->comp_ref_prob, default_comp_ref_p);
  vp9_copy(fc->single_ref_prob, default_single_ref_p);
  fc->tx_probs = default_tx_probs;
  vp9_copy(fc->txfm_partition_prob, default_txfm_partition_probs);
  vp9_copy(fc->intra_predictor_prob, vp9_intra_predictor_prob);
  vp9_copy(fc->skip_probs, default_skip_probs);
  vp9_copy(fc->inter_mode_probs, default_inter_mode_probs);
}

const vp9_tree_index vp9_switchable_interp_tree
                         [TREE_SIZE(SWITCHABLE_FILTERS)] = {
  -EIGHTTAP, 2,
  -EIGHTTAP_SMOOTH, -EIGHTTAP_SHARP
};

void vp9_adapt_mode_probs(VP9_COMMON *cm) {
  int i, j;
  FRAME_CONTEXT *fc = cm->fc;
  const FRAME_CONTEXT *pre_fc = &cm->frame_contexts[cm->frame_context_idx];
  const FRAME_COUNTS *counts = &cm->counts;

  for (i = 0; i < INTRA_INTER_CONTEXTS; i++)
    fc->intra_inter_prob[i] = mode_mv_merge_probs(pre_fc->intra_inter_prob[i],
                                                  counts->intra_inter[i]);
  for (i = 0; i < COMP_INTER_CONTEXTS; i++)
    fc->comp_inter_prob[i] = mode_mv_merge_probs(pre_fc->comp_inter_prob[i],
                                                 counts->comp_inter[i]);
  for (i = 0; i < REF_CONTEXTS; i++)
    fc->comp_ref_prob[i] = mode_mv_merge_probs(pre_fc->comp_ref_prob[i],
                                               counts->comp_ref[i]);
  for (i = 0; i < REF_CONTEXTS; i++)
    for (j = 0; j < 2; j++)
      fc->single_ref_prob[i][j] = mode_mv_merge_probs(
          pre_fc->single_ref_prob[i][j], counts->single_ref[i][j]);

  for (i = 0; i < INTER_MODE_CONTEXTS; i++)
    vp9_tree_merge_probs(vp9_inter_mode_tree, pre_fc->inter_mode_probs[i],
                counts->inter_mode[i], fc->inter_mode_probs[i]);

  for (i = 0; i < BLOCK_SIZE_GROUPS; i++)
    vp9_tree_merge_probs(vp9_intra_mode_tree, pre_fc->y_mode_prob[i],
                counts->y_mode[i], fc->y_mode_prob[i]);

  for (i = 0; i < INTRA_MODES; ++i)
    vp9_tree_merge_probs(vp9_intra_mode_tree, pre_fc->uv_mode_prob[i],
                         counts->uv_mode[i], fc->uv_mode_prob[i]);

  for (i = 0; i < PARTITION_CONTEXTS; i++)
    vp9_tree_merge_probs(vp9_partition_tree, pre_fc->partition_prob[i],
                         counts->partition[i], fc->partition_prob[i]);

  if (cm->interp_filter == SWITCHABLE) {
    for (i = 0; i < SWITCHABLE_FILTER_CONTEXTS; i++)
      vp9_tree_merge_probs(vp9_switchable_interp_tree,
                           pre_fc->switchable_interp_prob[i],
                           counts->switchable_interp[i],
                           fc->switchable_interp_prob[i]);
  }

  if (cm->tx_mode == TX_MODE_SELECT) {
    int j;
    unsigned int branch_ct_8x8p[TX_SIZES - 3][2];
    unsigned int branch_ct_16x16p[TX_SIZES - 2][2];
    unsigned int branch_ct_32x32p[TX_SIZES - 1][2];

    for (i = 0; i < TX_SIZE_CONTEXTS; ++i) {
      tx_counts_to_branch_counts_8x8(counts->tx.p8x8[i], branch_ct_8x8p);
      for (j = 0; j < TX_SIZES - 3; ++j)
        fc->tx_probs.p8x8[i][j] = mode_mv_merge_probs(
            pre_fc->tx_probs.p8x8[i][j], branch_ct_8x8p[j]);

      tx_counts_to_branch_counts_16x16(counts->tx.p16x16[i], branch_ct_16x16p);
      for (j = 0; j < TX_SIZES - 2; ++j)
        fc->tx_probs.p16x16[i][j] = mode_mv_merge_probs(
            pre_fc->tx_probs.p16x16[i][j], branch_ct_16x16p[j]);

      tx_counts_to_branch_counts_32x32(counts->tx.p32x32[i], branch_ct_32x32p);
      for (j = 0; j < TX_SIZES - 1; ++j)
        fc->tx_probs.p32x32[i][j] = mode_mv_merge_probs(
            pre_fc->tx_probs.p32x32[i][j], branch_ct_32x32p[j]);
    }
  }

  for (i = 0; i < TXFM_PARTITION_CONTEXTS; ++i)
    fc->txfm_partition_prob[i] =
        mode_mv_merge_probs(pre_fc->txfm_partition_prob[i],
                            counts->txfm_partition[i]);

  for (i = 0; i < SKIP_CONTEXTS; ++i)
    fc->skip_probs[i] = mode_mv_merge_probs(
        pre_fc->skip_probs[i], counts->skip[i]);
}

static void set_default_lf_deltas(struct loopfilter *lf) {
  lf->mode_ref_delta_enabled = 1;
  lf->mode_ref_delta_update = 1;

  lf->ref_deltas[INTRA_FRAME] = 1;
  lf->ref_deltas[LAST_FRAME] = 0;
  lf->ref_deltas[GOLDEN_FRAME] = -1;
  lf->ref_deltas[ALTREF_FRAME] = -1;

  lf->mode_deltas[0] = 0;
  lf->mode_deltas[1] = 0;
}

void vp9_setup_past_independence(VP9_COMMON *cm) {
  // Reset the segment feature data to the default stats:
  // Features disabled, 0, with delta coding (Default state).
  struct loopfilter *const lf = &cm->lf;

  int i;
  vp9_clearall_segfeatures(&cm->seg);
  cm->seg.abs_delta = SEGMENT_DELTADATA;

  if (cm->last_frame_seg_map && !cm->frame_parallel_decode)
    vpx_memset(cm->last_frame_seg_map, 0, (cm->mi_rows * cm->mi_cols));

  if (cm->current_frame_seg_map)
    vpx_memset(cm->current_frame_seg_map, 0, (cm->mi_rows * cm->mi_cols));

  // Reset the mode ref deltas for loop filter
  vp9_zero(lf->last_ref_deltas);
  vp9_zero(lf->last_mode_deltas);
  set_default_lf_deltas(lf);

  // To force update of the sharpness
  lf->last_sharpness_level = -1;

  vp9_default_coef_probs(cm);
  vp9_init_mode_probs(cm->fc);
  vp9_init_mv_probs(cm);
  cm->fc->initialized = 1;

  if (cm->frame_type == KEY_FRAME ||
      cm->error_resilient_mode || cm->reset_frame_context == 3) {
    // Reset all frame contexts.
    for (i = 0; i < FRAME_CONTEXTS; ++i)
      cm->frame_contexts[i] = *cm->fc;
  } else if (cm->reset_frame_context == 2) {
    // Reset only the frame context specified in the frame header.
    cm->frame_contexts[cm->frame_context_idx] = *cm->fc;
  }

  // prev_mip will only be allocated in encoder.
  if (frame_is_intra_only(cm) && cm->prev_mip && !cm->frame_parallel_decode)
    vpx_memset(cm->prev_mip, 0, cm->mi_stride * (cm->mi_rows + 1) *
                                    sizeof(*cm->prev_mip));

  vp9_zero(cm->ref_frame_sign_bias);

  cm->frame_context_idx = 0;
}
