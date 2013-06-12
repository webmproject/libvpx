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
  { 144,  11,  54, 157, 195, 130,  46,  58, 108 } /* y = dc */,
  { 118,  15, 123, 148, 131, 101,  44,  93, 131 } /* y = v */,
  { 113,  12,  23, 188, 226, 142,  26,  32, 125 } /* y = h */,
  { 120,  11,  50, 123, 163, 135,  64,  77, 103 } /* y = d45 */,
  { 113,   9,  36, 155, 111, 157,  32,  44, 161 } /* y = d135 */,
  { 116,   9,  55, 176,  76,  96,  37,  61, 149 } /* y = d117 */,
  { 115,   9,  28, 141, 161, 167,  21,  25, 193 } /* y = d153 */,
  { 120,  12,  32, 145, 195, 142,  32,  38,  86 } /* y = d27 */,
  { 116,  12,  64, 120, 140, 125,  49, 115, 121 } /* y = d63 */,
  { 102,  19,  66, 162, 182, 122,  35,  59, 128 } /* y = tm */
};

static const vp9_prob default_if_y_probs[BLOCK_SIZE_GROUPS]
                                        [VP9_INTRA_MODES - 1] = {
  {  65,  32,  18, 144, 162, 194,  41,  51,  98 } /* block_size < 8x8 */,
  { 132,  68,  18, 165, 217, 196,  45,  40,  78 } /* block_size < 16x16 */,
  { 173,  80,  19, 176, 240, 193,  64,  35,  46 } /* block_size < 32x32 */,
  { 221, 135,  38, 194, 248, 121,  96,  85,  29 } /* block_size >= 32x32 */
};

static const vp9_prob default_if_uv_probs[VP9_INTRA_MODES]
                                         [VP9_INTRA_MODES - 1] = {
  { 120,   7,  76, 176, 208, 126,  28,  54, 103 } /* y = dc */,
  {  48,  12, 154, 155, 139,  90,  34, 117, 119 } /* y = v */,
  {  67,   6,  25, 204, 243, 158,  13,  21,  96 } /* y = h */,
  {  97,   5,  44, 131, 176, 139,  48,  68,  97 } /* y = d45 */,
  {  83,   5,  42, 156, 111, 152,  26,  49, 152 } /* y = d135 */,
  {  80,   5,  58, 178,  74,  83,  33,  62, 145 } /* y = d117 */,
  {  86,   5,  32, 154, 192, 168,  14,  22, 163 } /* y = d153 */,
  {  85,   5,  32, 156, 216, 148,  19,  29,  73 } /* y = d27 */,
  {  77,   7,  64, 116, 132, 122,  37, 126, 120 } /* y = d63 */,
  { 101,  21, 107, 181, 192, 103,  19,  67, 125 } /* y = tm */
};

const vp9_prob vp9_partition_probs[NUM_FRAME_TYPES][NUM_PARTITION_CONTEXTS]
                                  [PARTITION_TYPES - 1] = {
  { /* frame_type = keyframe */
    /* 8x8 -> 4x4 */
    { 158,  97,  94 } /* a/l both not split */,
    {  93,  24,  99 } /* a split, l not split */,
    {  85, 119,  44 } /* l split, a not split */,
    {  62,  59,  67 } /* a/l both split */,
    /* 16x16 -> 8x8 */
    { 149,  53,  53 } /* a/l both not split */,
    {  94,  20,  48 } /* a split, l not split */,
    {  83,  53,  24 } /* l split, a not split */,
    {  52,  18,  18 } /* a/l both split */,
    /* 32x32 -> 16x16 */
    { 150,  40,  39 } /* a/l both not split */,
    {  78,  12,  26 } /* a split, l not split */,
    {  67,  33,  11 } /* l split, a not split */,
    {  24,   7,   5 } /* a/l both split */,
    /* 64x64 -> 32x32 */
    { 174,  35,  49 } /* a/l both not split */,
    {  68,  11,  27 } /* a split, l not split */,
    {  57,  15,   9 } /* l split, a not split */,
    {  12,   3,   3 } /* a/l both split */
  }, { /* frame_type = interframe */
    /* 8x8 -> 4x4 */
    { 199, 122, 141 } /* a/l both not split */,
    { 147,  63, 159 } /* a split, l not split */,
    { 148, 133, 118 } /* l split, a not split */,
    { 121, 104, 114 } /* a/l both split */,
    /* 16x16 -> 8x8 */
    { 174,  73,  87 } /* a/l both not split */,
    {  92,  41,  83 } /* a split, l not split */,
    {  82,  99,  50 } /* l split, a not split */,
    {  53,  39,  39 } /* a/l both split */,
    /* 32x32 -> 16x16 */
    { 177,  58,  59 } /* a/l both not split */,
    {  68,  26,  63 } /* a split, l not split */,
    {  52,  79,  25 } /* l split, a not split */,
    {  17,  14,  12 } /* a/l both split */,
    /* 64x64 -> 32x32 */
    { 222,  34,  30 } /* a/l both not split */,
    {  72,  16,  44 } /* a split, l not split */,
    {  58,  32,  12 } /* l split, a not split */,
    {  10,   7,   6 } /* a/l both split */
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

const vp9_prob vp9_default_tx_probs_32x32p[TX_SIZE_CONTEXTS]
                                          [TX_SIZE_MAX_SB - 1] = {
  { 3, 136, 37, },
  { 5, 52, 13, },
};
const vp9_prob vp9_default_tx_probs_16x16p[TX_SIZE_CONTEXTS]
                                          [TX_SIZE_MAX_SB - 2] = {
  { 20, 152, },
  { 15, 101, },
};
const vp9_prob vp9_default_tx_probs_8x8p[TX_SIZE_CONTEXTS]
                                        [TX_SIZE_MAX_SB - 3] = {
  { 100, },
  { 66, },
};

void tx_counts_to_branch_counts_32x32(unsigned int *tx_count_32x32p,
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

void tx_counts_to_branch_counts_16x16(unsigned int *tx_count_16x16p,
                                      unsigned int (*ct_16x16p)[2]) {
  ct_16x16p[0][0] = tx_count_16x16p[TX_4X4];
  ct_16x16p[0][1] = tx_count_16x16p[TX_8X8] +
                    tx_count_16x16p[TX_16X16];
  ct_16x16p[1][0] = tx_count_16x16p[TX_8X8];
  ct_16x16p[1][1] = tx_count_16x16p[TX_16X16];
}

void tx_counts_to_branch_counts_8x8(unsigned int *tx_count_8x8p,
                                    unsigned int (*ct_8x8p)[2]) {
  ct_8x8p[0][0] =   tx_count_8x8p[TX_4X4];
  ct_8x8p[0][1] =   tx_count_8x8p[TX_8X8];
}

const vp9_prob vp9_default_mbskip_probs[MBSKIP_CONTEXTS] = {
  192, 128, 64
};

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
  vpx_memcpy(x->fc.tx_probs_32x32p, vp9_default_tx_probs_32x32p,
             sizeof(vp9_default_tx_probs_32x32p));
  vpx_memcpy(x->fc.tx_probs_16x16p, vp9_default_tx_probs_16x16p,
             sizeof(vp9_default_tx_probs_16x16p));
  vpx_memcpy(x->fc.tx_probs_8x8p, vp9_default_tx_probs_8x8p,
             sizeof(vp9_default_tx_probs_8x8p));
  vpx_memcpy(x->fc.mbskip_probs, vp9_default_mbskip_probs,
             sizeof(vp9_default_mbskip_probs));
}

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
  { 235, 162, },
  { 36, 255, },
  { 34, 3, },
  { 149, 144, },
};

// Indicates if the filter is interpolating or non-interpolating
const int vp9_is_interpolating_filter[SWITCHABLE + 1] = {1, 1, 1, 1, -1};

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
#define MODE_MAX_UPDATE_FACTOR 128
static int update_mode_ct(vp9_prob pre_prob, vp9_prob prob,
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

static int update_mode_ct2(vp9_prob pre_prob, unsigned int branch_ct[2]) {
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
    int j;
    unsigned int branch_ct_8x8p[TX_SIZE_MAX_SB - 3][2];
    unsigned int branch_ct_16x16p[TX_SIZE_MAX_SB - 2][2];
    unsigned int branch_ct_32x32p[TX_SIZE_MAX_SB - 1][2];
    for (i = 0; i < TX_SIZE_CONTEXTS; ++i) {
      tx_counts_to_branch_counts_8x8(cm->fc.tx_count_8x8p[i],
                                     branch_ct_8x8p);
      for (j = 0; j < TX_SIZE_MAX_SB - 3; ++j) {
        int factor;
        int count = branch_ct_8x8p[j][0] + branch_ct_8x8p[j][1];
        vp9_prob prob = get_binary_prob(branch_ct_8x8p[j][0],
                                        branch_ct_8x8p[j][1]);
        count = count > MODE_COUNT_SAT ? MODE_COUNT_SAT : count;
        factor = (MODE_MAX_UPDATE_FACTOR * count / MODE_COUNT_SAT);
        cm->fc.tx_probs_8x8p[i][j] = weighted_prob(
            cm->fc.pre_tx_probs_8x8p[i][j], prob, factor);
      }
    }
    for (i = 0; i < TX_SIZE_CONTEXTS; ++i) {
      tx_counts_to_branch_counts_16x16(cm->fc.tx_count_16x16p[i],
                                       branch_ct_16x16p);
      for (j = 0; j < TX_SIZE_MAX_SB - 2; ++j) {
        int factor;
        int count = branch_ct_16x16p[j][0] + branch_ct_16x16p[j][1];
        vp9_prob prob = get_binary_prob(branch_ct_16x16p[j][0],
                                        branch_ct_16x16p[j][1]);
        count = count > MODE_COUNT_SAT ? MODE_COUNT_SAT : count;
        factor = (MODE_MAX_UPDATE_FACTOR * count / MODE_COUNT_SAT);
        cm->fc.tx_probs_16x16p[i][j] = weighted_prob(
            cm->fc.pre_tx_probs_16x16p[i][j], prob, factor);
      }
    }
    for (i = 0; i < TX_SIZE_CONTEXTS; ++i) {
      tx_counts_to_branch_counts_32x32(cm->fc.tx_count_32x32p[i],
                                       branch_ct_32x32p);
      for (j = 0; j < TX_SIZE_MAX_SB - 1; ++j) {
        int factor;
        int count = branch_ct_32x32p[j][0] + branch_ct_32x32p[j][1];
        vp9_prob prob = get_binary_prob(branch_ct_32x32p[j][0],
                                        branch_ct_32x32p[j][1]);
        count = count > MODE_COUNT_SAT ? MODE_COUNT_SAT : count;
        factor = (MODE_MAX_UPDATE_FACTOR * count / MODE_COUNT_SAT);
        cm->fc.tx_probs_32x32p[i][j] = weighted_prob(
            cm->fc.pre_tx_probs_32x32p[i][j], prob, factor);
      }
    }
  }
  for (i = 0; i < MBSKIP_CONTEXTS; ++i)
    fc->mbskip_probs[i] = update_mode_ct2(fc->pre_mbskip_probs[i],
                                          fc->mbskip_count[i]);
}

static void set_default_lf_deltas(MACROBLOCKD *xd) {
  xd->mode_ref_lf_delta_enabled = 1;
  xd->mode_ref_lf_delta_update = 1;

  xd->ref_lf_deltas[INTRA_FRAME] = 1;
  xd->ref_lf_deltas[LAST_FRAME] = 0;
  xd->ref_lf_deltas[GOLDEN_FRAME] = -1;
  xd->ref_lf_deltas[ALTREF_FRAME] = -1;

  xd->mode_lf_deltas[0] = 0;              // Zero
  xd->mode_lf_deltas[1] = 0;               // New mv
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

  if ((cm->frame_type == KEY_FRAME) ||
      cm->error_resilient_mode || (cm->reset_frame_context == 3)) {
    // Reset all frame contexts.
    for (i = 0; i < NUM_FRAME_CONTEXTS; ++i)
      vpx_memcpy(&cm->frame_contexts[i], &cm->fc, sizeof(cm->fc));
  } else if (cm->reset_frame_context == 2) {
    // Reset only the frame context specified in the frame header.
    vpx_memcpy(&cm->frame_contexts[cm->frame_context_idx], &cm->fc,
               sizeof(cm->fc));
  }

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
