/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */


#include "vp9/common/vp9_onyxc_int.h"
#include "vp9/common/vp9_modecont.h"
#include "vp9/common/vp9_seg_common.h"
#include "vp9/common/vp9_alloccommon.h"
#include "vpx_mem/vpx_mem.h"

static const vp9_prob default_kf_uv_probs[VP9_INTRA_MODES]
                                         [VP9_INTRA_MODES - 1] = {
  { 149,  13,  48, 141, 174, 131,  54,  61, 109 } /* y = dc */,
  { 120,  17, 119, 132, 103, 103,  54, 100, 130 } /* y = v */,
  { 114,  16,  19, 177, 220, 145,  31,  33, 122 } /* y = h */,
  { 119,  12,  43, 102, 133, 133,  77,  90, 102 } /* y = d45 */,
  { 110,  10,  28, 144,  78, 158,  40,  49, 161 } /* y = d135 */,
  { 114,  10,  46, 169,  50,  96,  48,  70, 150 } /* y = d117 */,
  { 116,  10,  24, 125, 134, 168,  26,  27, 193 } /* y = d153 */,
  { 121,  14,  26, 124, 175, 143,  36,  37,  79 } /* y = d27 */,
  { 116,  13,  54, 100, 105, 122,  58, 126, 122 } /* y = d63 */,
  {  98,  22,  60, 147, 159, 124,  45,  68, 128 } /* y = tm */
};

static const vp9_prob default_if_y_probs[BLOCK_SIZE_GROUPS]
                                        [VP9_INTRA_MODES - 1] = {
  {  42,  31,  23, 150, 161, 193,  32,  53, 100 } /* block_size < 8x8 */,
  { 132,  58,  30, 160, 209, 195,  52,  47,  76 } /* block_size < 16x16 */,
  { 179,  85,  24, 168, 236, 198,  87,  45,  46 } /* block_size < 32x32 */,
  { 221, 176,  63, 133, 233, 121, 125, 105,  34 } /* block_size >= 32x32 */
};

static const vp9_prob default_if_uv_probs[VP9_INTRA_MODES]
                                         [VP9_INTRA_MODES - 1] = {
  { 115,   7,  78, 180, 210, 127,  34,  57, 104 } /* y = dc */,
  {  43,   9, 165, 140, 112,  93,  45, 125, 117 } /* y = v */,
  {  68,   6,  25, 206, 241, 154,  16,  23, 102 } /* y = h */,
  {  90,   5,  48, 117, 155, 134,  61,  88,  96 } /* y = d45 */,
  {  77,   5,  43, 148, 100, 147,  37,  60, 146 } /* y = d135 */,
  {  75,   5,  57, 167,  62,  91,  45,  76, 139 } /* y = d117 */,
  {  86,   4,  34, 155, 185, 163,  22,  29, 160 } /* y = d153 */,
  {  82,   5,  34, 155, 207, 144,  26,  38,  79 } /* y = d27 */,
  {  69,   6,  65, 105, 104, 122,  48, 131, 116 } /* y = d63 */,
  {  86,  16, 114, 177, 189, 108,  28,  72, 120 } /* y = tm */
};

const vp9_prob vp9_partition_probs[NUM_FRAME_TYPES][NUM_PARTITION_CONTEXTS]
                                  [PARTITION_TYPES - 1] = {
  { /* frame_type = keyframe */
    /* 8x8 -> 4x4 */
    { 164, 121, 109 } /* a/l both not split */,
    {  69,  11, 129 } /* a split, l not split */,
    {  52, 181,  37 } /* l split, a not split */,
    {  66,  71,  93 } /* a/l both split */,
    /* 16x16 -> 8x8 */
    { 154,  48,  43 } /* a/l both not split */,
    {  81,  11,  63 } /* a split, l not split */,
    {  67,  65,  17 } /* l split, a not split */,
    {  57,  18,  24 } /* a/l both split */,
    /* 32x32 -> 16x16 */
    { 156,  42,  35 } /* a/l both not split */,
    {  74,  10,  40 } /* a split, l not split */,
    {  59,  53,  10 } /* l split, a not split */,
    {  28,  10,   9 } /* a/l both split */,
    /* 64x64 -> 32x32 */
    { 168,  32,  43 } /* a/l both not split */,
    {  59,  13,  41 } /* a split, l not split */,
    {  60,  25,  10 } /* l split, a not split */,
    {  13,   5,   4 } /* a/l both split */
  }, { /* frame_type = interframe */
    /* 8x8 -> 4x4 */
    { 192, 121, 151 } /* a/l both not split */,
    { 134,  63, 162 } /* a split, l not split */,
    { 136, 134, 127 } /* l split, a not split */,
    { 101,  97, 131 } /* a/l both split */,
    /* 16x16 -> 8x8 */
    { 167,  67,  80 } /* a/l both not split */,
    {  87,  36,  70 } /* a split, l not split */,
    {  90,  61,  45 } /* l split, a not split */,
    {  46,  31,  32 } /* a/l both split */,
    /* 32x32 -> 16x16 */
    { 167,  63,  75 } /* a/l both not split */,
    {  67,  27,  61 } /* a split, l not split */,
    {  56,  87,  31 } /* l split, a not split */,
    {  15,  13,  11 } /* a/l both split */,
    /* 64x64 -> 32x32 */
    { 222,  45,  44 } /* a/l both not split */,
    {  62,  17,  62 } /* a split, l not split */,
    {  52,  65,  16 } /* l split, a not split */,
    {   9,   7,   6 } /* a/l both split */
  }
};

/* Array indices are identical to previously-existing INTRAMODECONTEXTNODES. */
const vp9_tree_index vp9_intra_mode_tree[VP9_INTRA_MODES * 2 - 2] = {
  -DC_PRED, 2,                      /* 0 = DC_NODE */
  -TM_PRED, 4,                      /* 1 = TM_NODE */
  -V_PRED, 6,                       /* 2 = V_NODE */
  8, 12,                            /* 3 = COM_NODE */
  -H_PRED, 10,                      /* 4 = H_NODE */
  -D135_PRED, -D117_PRED,           /* 5 = D135_NODE */
  -D45_PRED, 14,                    /* 6 = D45_NODE */
  -D63_PRED, 16,                    /* 7 = D63_NODE */
  -D153_PRED, -D27_PRED             /* 8 = D153_NODE */
};

const vp9_tree_index vp9_sb_mv_ref_tree[6] = {
  -ZEROMV, 2,
  -NEARESTMV, 4,
  -NEARMV, -NEWMV
};

const vp9_tree_index vp9_partition_tree[6] = {
  -PARTITION_NONE, 2,
  -PARTITION_HORZ, 4,
  -PARTITION_VERT, -PARTITION_SPLIT
};

struct vp9_token vp9_intra_mode_encodings[VP9_INTRA_MODES];

struct vp9_token vp9_sb_mv_ref_encoding_array[VP9_INTER_MODES];

struct vp9_token vp9_partition_encodings[PARTITION_TYPES];

static const vp9_prob default_intra_inter_p[INTRA_INTER_CONTEXTS] = {
  6, 87, 165, 213
};

static const vp9_prob default_comp_inter_p[COMP_INTER_CONTEXTS] = {
  25, 66, 106, 142, 183
};

static const vp9_prob default_comp_ref_p[REF_CONTEXTS] = {
  36, 93, 136, 205, 236
};

static const vp9_prob default_single_ref_p[REF_CONTEXTS][2] = {
  { 30, 17 },
  { 80, 66 },
  { 142, 129 },
  { 192, 178 },
  { 235, 248 },
};

void tx_counts_to_branch_counts(unsigned int *tx_count_32x32p,
                                unsigned int *tx_count_16x16p,
                                unsigned int *tx_count_8x8p,
                                unsigned int (*ct)[2]) {
#if TX_SIZE_PROBS == 6
  ct[0][0] = tx_count_8x8p[TX_4X4];
  ct[0][1] = tx_count_8x8p[TX_8X8];
  ct[1][0] = tx_count_16x16p[TX_4X4];
  ct[1][1] = tx_count_16x16p[TX_8X8] + tx_count_16x16p[TX_16X16];
  ct[2][0] = tx_count_16x16p[TX_8X8];
  ct[2][1] = tx_count_16x16p[TX_16X16];
  ct[3][0] = tx_count_32x32p[TX_4X4];
  ct[3][1] = tx_count_32x32p[TX_8X8] + tx_count_32x32p[TX_16X16] +
             tx_count_32x32p[TX_32X32];
  ct[4][0] = tx_count_32x32p[TX_8X8];
  ct[4][1] = tx_count_32x32p[TX_16X16] + tx_count_32x32p[TX_32X32];
  ct[5][0] = tx_count_32x32p[TX_16X16];
  ct[5][1] = tx_count_32x32p[TX_32X32];
#else
  ct[0][0] = tx_count_32x32p[TX_4X4] +
             tx_count_16x16p[TX_4X4] +
             tx_count_8x8p[TX_4X4];
  ct[0][1] = tx_count_32x32p[TX_8X8] +
             tx_count_32x32p[TX_16X16] +
             tx_count_32x32p[TX_32X32] +
             tx_count_16x16p[TX_8X8] +
             tx_count_16x16p[TX_16X16] +
             tx_count_8x8p[TX_8X8];
  ct[1][0] = tx_count_32x32p[TX_8X8] +
             tx_count_16x16p[TX_8X8];
  ct[1][1] = tx_count_32x32p[TX_16X16] +
             tx_count_32x32p[TX_32X32] +
             tx_count_16x16p[TX_16X16];
  ct[2][0] = tx_count_32x32p[TX_16X16];
  ct[2][1] = tx_count_32x32p[TX_32X32];
#endif
}

#if TX_SIZE_PROBS == 6
const vp9_prob vp9_default_tx_probs[TX_SIZE_PROBS] = {
  96, 96, 96, 96, 96, 96
};
#else
const vp9_prob vp9_default_tx_probs[TX_SIZE_PROBS] = {
  96, 96, 96
};
#endif

void vp9_init_mbmode_probs(VP9_COMMON *x) {
  vpx_memcpy(x->fc.uv_mode_prob, default_if_uv_probs,
             sizeof(default_if_uv_probs));
  vpx_memcpy(x->kf_uv_mode_prob, default_kf_uv_probs,
             sizeof(default_kf_uv_probs));
  vpx_memcpy(x->fc.y_mode_prob, default_if_y_probs,
             sizeof(default_if_y_probs));

  vpx_memcpy(x->fc.switchable_interp_prob, vp9_switchable_interp_prob,
             sizeof(vp9_switchable_interp_prob));

  vpx_memcpy(x->fc.partition_prob, vp9_partition_probs,
             sizeof(vp9_partition_probs));

  vpx_memcpy(x->fc.intra_inter_prob, default_intra_inter_p,
             sizeof(default_intra_inter_p));
  vpx_memcpy(x->fc.comp_inter_prob, default_comp_inter_p,
             sizeof(default_comp_inter_p));
  vpx_memcpy(x->fc.comp_ref_prob, default_comp_ref_p,
             sizeof(default_comp_ref_p));
  vpx_memcpy(x->fc.single_ref_prob, default_single_ref_p,
             sizeof(default_single_ref_p));
  vpx_memcpy(x->fc.tx_probs, vp9_default_tx_probs,
             sizeof(vp9_default_tx_probs));
}

#if VP9_SWITCHABLE_FILTERS == 3
const vp9_tree_index vp9_switchable_interp_tree[VP9_SWITCHABLE_FILTERS*2-2] = {
  -0, 2,
  -1, -2
};
struct vp9_token vp9_switchable_interp_encodings[VP9_SWITCHABLE_FILTERS];
const INTERPOLATIONFILTERTYPE vp9_switchable_interp[VP9_SWITCHABLE_FILTERS] = {
  EIGHTTAP, EIGHTTAP_SMOOTH, EIGHTTAP_SHARP};
const int vp9_switchable_interp_map[SWITCHABLE+1] = {1, 0, 2, -1, -1};
const vp9_prob vp9_switchable_interp_prob [VP9_SWITCHABLE_FILTERS+1]
                                          [VP9_SWITCHABLE_FILTERS-1] = {
  {248, 192}, { 32, 248}, { 32,  32}, {192, 160}
};
#elif VP9_SWITCHABLE_FILTERS == 2
const vp9_tree_index vp9_switchable_interp_tree[VP9_SWITCHABLE_FILTERS*2-2] = {
  -0, -1,
};
struct vp9_token vp9_switchable_interp_encodings[VP9_SWITCHABLE_FILTERS];
const vp9_prob vp9_switchable_interp_prob [VP9_SWITCHABLE_FILTERS+1]
                                          [VP9_SWITCHABLE_FILTERS-1] = {
  {248},
  { 64},
  {192},
};
const INTERPOLATIONFILTERTYPE vp9_switchable_interp[VP9_SWITCHABLE_FILTERS] = {
  EIGHTTAP, EIGHTTAP_SHARP};
const int vp9_switchable_interp_map[SWITCHABLE+1] = {-1, 0, 1, -1, -1};
#endif  // VP9_SWITCHABLE_FILTERS

// Indicates if the filter is interpolating or non-interpolating
// Note currently only the EIGHTTAP_SMOOTH is non-interpolating
const int vp9_is_interpolating_filter[SWITCHABLE + 1] = {0, 1, 1, 1, -1};

void vp9_entropy_mode_init() {
  vp9_tokens_from_tree(vp9_intra_mode_encodings, vp9_intra_mode_tree);
  vp9_tokens_from_tree(vp9_switchable_interp_encodings,
                       vp9_switchable_interp_tree);
  vp9_tokens_from_tree(vp9_partition_encodings, vp9_partition_tree);

  vp9_tokens_from_tree_offset(vp9_sb_mv_ref_encoding_array,
                              vp9_sb_mv_ref_tree, NEARESTMV);
}

void vp9_init_mode_contexts(VP9_COMMON *pc) {
  vpx_memset(pc->fc.inter_mode_counts, 0, sizeof(pc->fc.inter_mode_counts));
  vpx_memcpy(pc->fc.inter_mode_probs,
             vp9_default_inter_mode_probs,
             sizeof(vp9_default_inter_mode_probs));
}

void vp9_accum_mv_refs(VP9_COMMON *pc,
                       MB_PREDICTION_MODE m,
                       const int context) {
  unsigned int (*inter_mode_counts)[VP9_INTER_MODES - 1][2] =
      pc->fc.inter_mode_counts;

  if (m == ZEROMV) {
    ++inter_mode_counts[context][0][0];
  } else {
    ++inter_mode_counts[context][0][1];
    if (m == NEARESTMV) {
      ++inter_mode_counts[context][1][0];
    } else {
      ++inter_mode_counts[context][1][1];
      if (m == NEARMV) {
        ++inter_mode_counts[context][2][0];
      } else {
        ++inter_mode_counts[context][2][1];
      }
    }
  }
}

#define MVREF_COUNT_SAT 20
#define MVREF_MAX_UPDATE_FACTOR 128
void vp9_adapt_mode_context(VP9_COMMON *pc) {
  int i, j;
  unsigned int (*inter_mode_counts)[VP9_INTER_MODES - 1][2] =
      pc->fc.inter_mode_counts;
  vp9_prob (*mode_context)[VP9_INTER_MODES - 1] = pc->fc.inter_mode_probs;

  for (j = 0; j < INTER_MODE_CONTEXTS; j++) {
    for (i = 0; i < VP9_INTER_MODES - 1; i++) {
      int count = inter_mode_counts[j][i][0] + inter_mode_counts[j][i][1];
      int factor;
      count = count > MVREF_COUNT_SAT ? MVREF_COUNT_SAT : count;
      factor = (MVREF_MAX_UPDATE_FACTOR * count / MVREF_COUNT_SAT);
      mode_context[j][i] = weighted_prob(
          pc->fc.pre_inter_mode_probs[j][i],
          get_binary_prob(inter_mode_counts[j][i][0],
                          inter_mode_counts[j][i][1]),
          factor);
    }
  }
}

#define MODE_COUNT_SAT 20
#define MODE_MAX_UPDATE_FACTOR 144
static int update_mode_ct(int pre_prob, int prob,
                          unsigned int branch_ct[2]) {
  int factor, count = branch_ct[0] + branch_ct[1];
  count = count > MODE_COUNT_SAT ? MODE_COUNT_SAT : count;
  factor = (MODE_MAX_UPDATE_FACTOR * count / MODE_COUNT_SAT);
  return weighted_prob(pre_prob, prob, factor);
}

static void update_mode_probs(int n_modes,
                              const vp9_tree_index *tree, unsigned int *cnt,
                              vp9_prob *pre_probs, vp9_prob *dst_probs,
                              unsigned int tok0_offset) {
#define MAX_PROBS 32
  vp9_prob probs[MAX_PROBS];
  unsigned int branch_ct[MAX_PROBS][2];
  int t;

  assert(n_modes - 1 < MAX_PROBS);
  vp9_tree_probs_from_distribution(tree, probs, branch_ct, cnt, tok0_offset);
  for (t = 0; t < n_modes - 1; ++t)
    dst_probs[t] = update_mode_ct(pre_probs[t], probs[t], branch_ct[t]);
}

static int update_mode_ct2(int pre_prob, unsigned int branch_ct[2]) {
  return update_mode_ct(pre_prob, get_binary_prob(branch_ct[0],
                                                  branch_ct[1]), branch_ct);
}

// #define MODE_COUNT_TESTING
void vp9_adapt_mode_probs(VP9_COMMON *cm) {
  int i, j;
  FRAME_CONTEXT *fc = &cm->fc;
#ifdef MODE_COUNT_TESTING
  int t;

  printf("static const unsigned int\nymode_counts"
         "[VP9_INTRA_MODES] = {\n");
  for (t = 0; t < VP9_INTRA_MODES; ++t)
    printf("%d, ", fc->ymode_counts[t]);
  printf("};\n");
  printf("static const unsigned int\nuv_mode_counts"
         "[VP9_INTRA_MODES] [VP9_INTRA_MODES] = {\n");
  for (i = 0; i < VP9_INTRA_MODES; ++i) {
    printf("  {");
    for (t = 0; t < VP9_INTRA_MODES; ++t)
      printf("%d, ", fc->uv_mode_counts[i][t]);
    printf("},\n");
  }
  printf("};\n");
  printf("static const unsigned int\nbmode_counts"
         "[VP9_NKF_BINTRAMODES] = {\n");
  for (t = 0; t < VP9_NKF_BINTRAMODES; ++t)
    printf("%d, ", fc->bmode_counts[t]);
  printf("};\n");
  printf("static const unsigned int\ni8x8_mode_counts"
         "[VP9_I8X8_MODES] = {\n");
  for (t = 0; t < VP9_I8X8_MODES; ++t)
    printf("%d, ", fc->i8x8_mode_counts[t]);
  printf("};\n");
  printf("static const unsigned int\nmbsplit_counts"
         "[VP9_NUMMBSPLITS] = {\n");
  for (t = 0; t < VP9_NUMMBSPLITS; ++t)
    printf("%d, ", fc->mbsplit_counts[t]);
  printf("};\n");
#endif

  for (i = 0; i < INTRA_INTER_CONTEXTS; i++)
    fc->intra_inter_prob[i] = update_mode_ct2(fc->pre_intra_inter_prob[i],
                                              fc->intra_inter_count[i]);
  for (i = 0; i < COMP_INTER_CONTEXTS; i++)
    fc->comp_inter_prob[i] = update_mode_ct2(fc->pre_comp_inter_prob[i],
                                             fc->comp_inter_count[i]);
  for (i = 0; i < REF_CONTEXTS; i++)
    fc->comp_ref_prob[i] = update_mode_ct2(fc->pre_comp_ref_prob[i],
                                           fc->comp_ref_count[i]);
  for (i = 0; i < REF_CONTEXTS; i++)
    for (j = 0; j < 2; j++)
      fc->single_ref_prob[i][j] = update_mode_ct2(fc->pre_single_ref_prob[i][j],
                                                  fc->single_ref_count[i][j]);

  for (i = 0; i < BLOCK_SIZE_GROUPS; i++)
    update_mode_probs(VP9_INTRA_MODES, vp9_intra_mode_tree,
                      fc->y_mode_counts[i], fc->pre_y_mode_prob[i],
                      fc->y_mode_prob[i], 0);

  for (i = 0; i < VP9_INTRA_MODES; ++i)
    update_mode_probs(VP9_INTRA_MODES, vp9_intra_mode_tree,
                      fc->uv_mode_counts[i], fc->pre_uv_mode_prob[i],
                      fc->uv_mode_prob[i], 0);

  for (i = 0; i < NUM_PARTITION_CONTEXTS; i++)
    update_mode_probs(PARTITION_TYPES, vp9_partition_tree,
                      fc->partition_counts[i], fc->pre_partition_prob[i],
                      fc->partition_prob[INTER_FRAME][i], 0);

  if (cm->mcomp_filter_type == SWITCHABLE) {
    for (i = 0; i <= VP9_SWITCHABLE_FILTERS; i++) {
      update_mode_probs(VP9_SWITCHABLE_FILTERS, vp9_switchable_interp_tree,
                        fc->switchable_interp_count[i],
                        fc->pre_switchable_interp_prob[i],
                        fc->switchable_interp_prob[i], 0);
    }
  }
  if (cm->txfm_mode == TX_MODE_SELECT) {
    unsigned int branch_ct[TX_SIZE_PROBS][2];
    tx_counts_to_branch_counts(cm->fc.tx_count_32x32p,
                               cm->fc.tx_count_16x16p,
                               cm->fc.tx_count_8x8p, branch_ct);
    for (i = 0; i < TX_SIZE_PROBS; ++i) {
      int factor;
      int count = branch_ct[i][0] + branch_ct[i][1];
      vp9_prob prob = get_binary_prob(branch_ct[i][0], branch_ct[i][1]);
      count = count > MODE_COUNT_SAT ? MODE_COUNT_SAT : count;
      factor = (MODE_MAX_UPDATE_FACTOR * count / MODE_COUNT_SAT);
      cm->fc.tx_probs[i] = weighted_prob(cm->fc.pre_tx_probs[i], prob, factor);
    }
  }
}

static void set_default_lf_deltas(MACROBLOCKD *xd) {
  xd->mode_ref_lf_delta_enabled = 1;
  xd->mode_ref_lf_delta_update = 1;

  xd->ref_lf_deltas[INTRA_FRAME] = 1;
  xd->ref_lf_deltas[LAST_FRAME] = 0;
  xd->ref_lf_deltas[GOLDEN_FRAME] = -1;
  xd->ref_lf_deltas[ALTREF_FRAME] = -1;

  xd->mode_lf_deltas[0] = 2;               // I4X4_PRED
  xd->mode_lf_deltas[1] = -1;              // Zero
  xd->mode_lf_deltas[2] = 1;               // New mv
  xd->mode_lf_deltas[3] = 2;               // Split mv
}

void vp9_setup_past_independence(VP9_COMMON *cm, MACROBLOCKD *xd) {
  // Reset the segment feature data to the default stats:
  // Features disabled, 0, with delta coding (Default state).
  int i;
  vp9_clearall_segfeatures(xd);
  xd->mb_segment_abs_delta = SEGMENT_DELTADATA;
  if (cm->last_frame_seg_map)
    vpx_memset(cm->last_frame_seg_map, 0, (cm->mi_rows * cm->mi_cols));

  // Reset the mode ref deltas for loop filter
  vpx_memset(xd->last_ref_lf_deltas, 0, sizeof(xd->last_ref_lf_deltas));
  vpx_memset(xd->last_mode_lf_deltas, 0, sizeof(xd->last_mode_lf_deltas));
  set_default_lf_deltas(xd);

  vp9_default_coef_probs(cm);
  vp9_init_mbmode_probs(cm);
  vpx_memcpy(cm->kf_y_mode_prob, vp9_kf_default_bmode_probs,
             sizeof(vp9_kf_default_bmode_probs));
  vp9_init_mv_probs(cm);

  // To force update of the sharpness
  cm->last_sharpness_level = -1;

  vp9_init_mode_contexts(cm);

  for (i = 0; i < NUM_FRAME_CONTEXTS; i++)
    vpx_memcpy(&cm->frame_contexts[i], &cm->fc, sizeof(cm->fc));

  vpx_memset(cm->prev_mip, 0,
             cm->mode_info_stride * (cm->mi_rows + 1) * sizeof(MODE_INFO));
  vpx_memset(cm->mip, 0,
             cm->mode_info_stride * (cm->mi_rows + 1) * sizeof(MODE_INFO));

  vp9_update_mode_info_border(cm, cm->mip);
  vp9_update_mode_info_in_image(cm, cm->mi);

  vp9_update_mode_info_border(cm, cm->prev_mip);
  vp9_update_mode_info_in_image(cm, cm->prev_mi);

  vpx_memset(cm->ref_frame_sign_bias, 0, sizeof(cm->ref_frame_sign_bias));

  cm->frame_context_idx = 0;
}
