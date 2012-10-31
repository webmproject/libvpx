/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */


#include "onyxc_int.h"
#include "modecont.h"
#include "vpx_mem/vpx_mem.h"


const unsigned int kf_y_mode_cts[8][VP8_YMODES] = {
  /* DC V   H  D45 135 117 153 D27 D63 TM i8x8 BPRED */
  {12,  6,  5,  5,  5,  5,  5,  5,  5,  2, 22, 200},
  {25, 13, 13,  7,  7,  7,  7,  7,  7,  6, 27, 160},
  {31, 17, 18,  8,  8,  8,  8,  8,  8,  9, 26, 139},
  {40, 22, 23,  8,  8,  8,  8,  8,  8, 12, 27, 116},
  {53, 26, 28,  8,  8,  8,  8,  8,  8, 13, 26,  94},
  {68, 33, 35,  8,  8,  8,  8,  8,  8, 17, 20,  68},
  {78, 38, 38,  8,  8,  8,  8,  8,  8, 19, 16,  52},
  {89, 42, 42,  8,  8,  8,  8,  8,  8, 21, 12,  34},
};

static const unsigned int y_mode_cts  [VP8_YMODES] = {
  /* DC V   H  D45 135 117 153 D27 D63 TM i8x8 BPRED */
  98, 19, 15, 14, 14, 14, 14, 12, 12, 13, 16, 70
};

static const unsigned int uv_mode_cts [VP8_YMODES] [VP8_UV_MODES] = {
  /* DC   V   H  D45 135 117 153 D27 D63 TM */
  { 200, 15, 15, 10, 10, 10, 10, 10, 10,  6}, /* DC */
  { 130, 75, 10, 10, 10, 10, 10, 10, 10,  6}, /* V */
  { 130, 10, 75, 10, 10, 10, 10, 10, 10,  6}, /* H */
  { 130, 15, 10, 75, 10, 10, 10, 10, 10,  6}, /* D45 */
  { 150, 15, 10, 10, 75, 10, 10, 10, 10,  6}, /* D135 */
  { 150, 15, 10, 10, 10, 75, 10, 10, 10,  6}, /* D117 */
  { 150, 15, 10, 10, 10, 10, 75, 10, 10,  6}, /* D153 */
  { 150, 15, 10, 10, 10, 10, 10, 75, 10,  6}, /* D27 */
  { 150, 15, 10, 10, 10, 10, 10, 10, 75,  6}, /* D63 */
  { 160, 30, 30, 10, 10, 10, 10, 10, 10, 16}, /* TM */
  { 132, 46, 40, 10, 10, 10, 10, 10, 10, 18}, /* i8x8 - never used */
  { 150, 35, 41, 10, 10, 10, 10, 10, 10, 10}, /* BPRED */
};

static const unsigned int i8x8_mode_cts  [VP8_I8X8_MODES] = {
  /* DC V   H D45 135 117 153 D27 D63  TM */
  73, 49, 61, 30, 30, 30, 30, 30, 30, 13
};

static const unsigned int kf_uv_mode_cts [VP8_YMODES] [VP8_UV_MODES] = {
  // DC   V   H  D45 135 117 153 D27 D63 TM
  { 160, 24, 24, 20, 20, 20, 20, 20, 20,  8}, /* DC */
  { 102, 64, 30, 20, 20, 20, 20, 20, 20, 10}, /* V */
  { 102, 30, 64, 20, 20, 20, 20, 20, 20, 10}, /* H */
  { 102, 33, 20, 64, 20, 20, 20, 20, 20, 14}, /* D45 */
  { 102, 33, 20, 20, 64, 20, 20, 20, 20, 14}, /* D135 */
  { 122, 33, 20, 20, 20, 64, 20, 20, 20, 14}, /* D117 */
  { 102, 33, 20, 20, 20, 20, 64, 20, 20, 14}, /* D153 */
  { 102, 33, 20, 20, 20, 20, 20, 64, 20, 14}, /* D27 */
  { 102, 33, 20, 20, 20, 20, 20, 20, 64, 14}, /* D63 */
  { 132, 36, 30, 20, 20, 20, 20, 20, 20, 18}, /* TM */
  { 122, 41, 35, 20, 20, 20, 20, 20, 20, 18}, /* i8x8 - never used */
  { 122, 41, 35, 20, 20, 20, 20, 20, 20, 18}, /* BPRED */
};

static const unsigned int bmode_cts[VP8_BINTRAMODES] = {
  /* DC    TM     VE     HE   LD    RD    VR    VL    HD    HU */
  43891, 17694, 10036, 3920, 3363, 2546, 5119, 3221, 2471, 1723
};

typedef enum {
  SUBMVREF_NORMAL,
  SUBMVREF_LEFT_ZED,
  SUBMVREF_ABOVE_ZED,
  SUBMVREF_LEFT_ABOVE_SAME,
  SUBMVREF_LEFT_ABOVE_ZED
} sumvfref_t;

int vp9_mv_cont(const int_mv *l, const int_mv *a) {
  int lez = (l->as_int == 0);
  int aez = (a->as_int == 0);
  int lea = (l->as_int == a->as_int);

  if (lea && lez)
    return SUBMVREF_LEFT_ABOVE_ZED;

  if (lea)
    return SUBMVREF_LEFT_ABOVE_SAME;

  if (aez)
    return SUBMVREF_ABOVE_ZED;

  if (lez)
    return SUBMVREF_LEFT_ZED;

  return SUBMVREF_NORMAL;
}

const vp8_prob vp8_sub_mv_ref_prob [VP8_SUBMVREFS - 1] = { 180, 162, 25};

const vp8_prob vp8_sub_mv_ref_prob2 [SUBMVREF_COUNT][VP8_SUBMVREFS - 1] = {
  { 147, 136, 18 },
  { 106, 145, 1  },
  { 179, 121, 1  },
  { 223, 1, 34 },
  { 208, 1, 1  }
};



vp8_mbsplit vp8_mbsplits [VP8_NUMMBSPLITS] = {
  {
    0,  0,  0,  0,
    0,  0,  0,  0,
    1,  1,  1,  1,
    1,  1,  1,  1,
  },
  {
    0,  0,  1,  1,
    0,  0,  1,  1,
    0,  0,  1,  1,
    0,  0,  1,  1,
  },
  {
    0,  0,  1,  1,
    0,  0,  1,  1,
    2,  2,  3,  3,
    2,  2,  3,  3,
  },
  {
    0,  1,  2,  3,
    4,  5,  6,  7,
    8,  9,  10, 11,
    12, 13, 14, 15,
  },
};

const int vp8_mbsplit_count [VP8_NUMMBSPLITS] = { 2, 2, 4, 16};

const vp8_prob vp8_mbsplit_probs [VP8_NUMMBSPLITS - 1] = { 110, 111, 150};


/* Array indices are identical to previously-existing INTRAMODECONTEXTNODES. */

const vp8_tree_index vp8_bmode_tree[VP8_BINTRAMODES * 2 - 2] = /* INTRAMODECONTEXTNODE value */
{
  -B_DC_PRED, 2,                             /* 0 = DC_NODE */
  -B_TM_PRED, 4,                            /* 1 = TM_NODE */
  -B_VE_PRED, 6,                           /* 2 = VE_NODE */
  8, 12,                                  /* 3 = COM_NODE */
  -B_HE_PRED, 10,                        /* 4 = HE_NODE */
  -B_RD_PRED, -B_VR_PRED,               /* 5 = RD_NODE */
  -B_LD_PRED, 14,                        /* 6 = LD_NODE */
  -B_VL_PRED, 16,                      /* 7 = VL_NODE */
  -B_HD_PRED, -B_HU_PRED             /* 8 = HD_NODE */
};

/* Again, these trees use the same probability indices as their
   explicitly-programmed predecessors. */
const vp8_tree_index vp8_ymode_tree[VP8_YMODES * 2 - 2] = {
  2, 14,
  -DC_PRED, 4,
  6, 8,
  -D45_PRED, -D135_PRED,
  10, 12,
  -D117_PRED, -D153_PRED,
  -D27_PRED, -D63_PRED,
  16, 18,
  -V_PRED, -H_PRED,
  -TM_PRED, 20,
  -B_PRED, -I8X8_PRED
};

const vp8_tree_index vp8_kf_ymode_tree[VP8_YMODES * 2 - 2] = {
  2, 14,
  -DC_PRED, 4,
  6, 8,
  -D45_PRED, -D135_PRED,
  10, 12,
  -D117_PRED, -D153_PRED,
  -D27_PRED, -D63_PRED,
  16, 18,
  -V_PRED, -H_PRED,
  -TM_PRED, 20,
  -B_PRED, -I8X8_PRED
};

const vp8_tree_index vp8_i8x8_mode_tree[VP8_I8X8_MODES * 2 - 2] = {
  2, 14,
  -DC_PRED, 4,
  6, 8,
  -D45_PRED, -D135_PRED,
  10, 12,
  -D117_PRED, -D153_PRED,
  -D27_PRED, -D63_PRED,
  -V_PRED, 16,
  -H_PRED, -TM_PRED
};

const vp8_tree_index vp8_uv_mode_tree[VP8_UV_MODES * 2 - 2] = {
  2, 14,
  -DC_PRED, 4,
  6, 8,
  -D45_PRED, -D135_PRED,
  10, 12,
  -D117_PRED, -D153_PRED,
  -D27_PRED, -D63_PRED,
  -V_PRED, 16,
  -H_PRED, -TM_PRED
};

const vp8_tree_index vp8_mbsplit_tree[6] = {
  -PARTITIONING_4X4,   2,
  -PARTITIONING_8X8,   4,
  -PARTITIONING_16X8, -PARTITIONING_8X16,
};

const vp8_tree_index vp8_mv_ref_tree[8] = {
  -ZEROMV, 2,
  -NEARESTMV, 4,
  -NEARMV, 6,
  -NEWMV, -SPLITMV
};

#if CONFIG_SUPERBLOCKS
const vp8_tree_index vp8_sb_mv_ref_tree[6] = {
  -ZEROMV, 2,
  -NEARESTMV, 4,
  -NEARMV, -NEWMV
};
#endif

const vp8_tree_index vp8_sub_mv_ref_tree[6] = {
  -LEFT4X4, 2,
  -ABOVE4X4, 4,
  -ZERO4X4, -NEW4X4
};


struct vp8_token_struct vp8_bmode_encodings   [VP8_BINTRAMODES];
struct vp8_token_struct vp8_ymode_encodings   [VP8_YMODES];
#if CONFIG_SUPERBLOCKS
struct vp8_token_struct vp8_sb_kf_ymode_encodings [VP8_I32X32_MODES];
#endif
struct vp8_token_struct vp8_kf_ymode_encodings [VP8_YMODES];
struct vp8_token_struct vp8_uv_mode_encodings  [VP8_UV_MODES];
struct vp8_token_struct vp8_i8x8_mode_encodings  [VP8_I8X8_MODES];
struct vp8_token_struct vp8_mbsplit_encodings [VP8_NUMMBSPLITS];

struct vp8_token_struct vp8_mv_ref_encoding_array    [VP8_MVREFS];
#if CONFIG_SUPERBLOCKS
struct vp8_token_struct vp8_sb_mv_ref_encoding_array  [VP8_MVREFS];
#endif
struct vp8_token_struct vp8_sub_mv_ref_encoding_array [VP8_SUBMVREFS];



void vp9_init_mbmode_probs(VP8_COMMON *x) {
  unsigned int bct [VP8_YMODES] [2];      /* num Ymodes > num UV modes */

  vp9_tree_probs_from_distribution(VP8_YMODES, vp8_ymode_encodings,
    vp8_ymode_tree, x->fc.ymode_prob, bct, y_mode_cts, 256, 1);
  {
    int i;
    for (i = 0; i < 8; i++) {
      vp9_tree_probs_from_distribution(
        VP8_YMODES, vp8_kf_ymode_encodings, vp8_kf_ymode_tree,
        x->kf_ymode_prob[i], bct, kf_y_mode_cts[i],
        256, 1);
#if CONFIG_SUPERBLOCKS
      vp9_tree_probs_from_distribution(
        VP8_I32X32_MODES, vp8_sb_kf_ymode_encodings, vp8_sb_ymode_tree,
        x->sb_kf_ymode_prob[i], bct, kf_y_mode_cts[i],
        256, 1);
#endif
    }
  }
  {
    int i;
    for (i = 0; i < VP8_YMODES; i++) {
      vp9_tree_probs_from_distribution(
        VP8_UV_MODES, vp8_uv_mode_encodings, vp8_uv_mode_tree,
        x->kf_uv_mode_prob[i], bct, kf_uv_mode_cts[i],
        256, 1);
      vp9_tree_probs_from_distribution(
        VP8_UV_MODES, vp8_uv_mode_encodings, vp8_uv_mode_tree,
        x->fc.uv_mode_prob[i], bct, uv_mode_cts[i],
        256, 1);
    }
  }

  vp9_tree_probs_from_distribution(
    VP8_I8X8_MODES, vp8_i8x8_mode_encodings, vp8_i8x8_mode_tree,
    x->fc.i8x8_mode_prob, bct, i8x8_mode_cts,
    256, 1);

  vpx_memcpy(x->fc.sub_mv_ref_prob, vp8_sub_mv_ref_prob2, sizeof(vp8_sub_mv_ref_prob2));
  vpx_memcpy(x->fc.mbsplit_prob, vp8_mbsplit_probs, sizeof(vp8_mbsplit_probs));
  vpx_memcpy(x->fc.switchable_interp_prob, vp8_switchable_interp_prob,
             sizeof(vp8_switchable_interp_prob));
}


static void intra_bmode_probs_from_distribution(
  vp8_prob p [VP8_BINTRAMODES - 1],
  unsigned int branch_ct [VP8_BINTRAMODES - 1] [2],
  const unsigned int events [VP8_BINTRAMODES]) {
  vp9_tree_probs_from_distribution(VP8_BINTRAMODES, vp8_bmode_encodings,
    vp8_bmode_tree, p, branch_ct, events, 256, 1);
}

void vp9_default_bmode_probs(vp8_prob p [VP8_BINTRAMODES - 1]) {
  unsigned int branch_ct [VP8_BINTRAMODES - 1] [2];
  intra_bmode_probs_from_distribution(p, branch_ct, bmode_cts);
}

void vp9_kf_default_bmode_probs(vp8_prob p [VP8_BINTRAMODES] [VP8_BINTRAMODES] [VP8_BINTRAMODES - 1]) {
  unsigned int branch_ct [VP8_BINTRAMODES - 1] [2];

  int i = 0;

  do {
    int j = 0;

    do {
      intra_bmode_probs_from_distribution(
        p[i][j], branch_ct, vp8_kf_default_bmode_counts[i][j]);

    } while (++j < VP8_BINTRAMODES);
  } while (++i < VP8_BINTRAMODES);
}

#if VP8_SWITCHABLE_FILTERS == 3
const vp8_tree_index vp8_switchable_interp_tree[VP8_SWITCHABLE_FILTERS*2-2] = {
  -0, 2,
  -1, -2
};
struct vp8_token_struct vp8_switchable_interp_encodings[VP8_SWITCHABLE_FILTERS];
const INTERPOLATIONFILTERTYPE vp8_switchable_interp[VP8_SWITCHABLE_FILTERS] = {
  EIGHTTAP, SIXTAP, EIGHTTAP_SHARP};
const int vp8_switchable_interp_map[SWITCHABLE+1] = {1, -1, 0, 2, -1};
const vp8_prob vp8_switchable_interp_prob [VP8_SWITCHABLE_FILTERS+1]
                                          [VP8_SWITCHABLE_FILTERS-1] = {
  {248, 192}, { 32, 248}, { 32,  32}, {192, 160}
};
#elif VP8_SWITCHABLE_FILTERS == 2
const vp8_tree_index vp8_switchable_interp_tree[VP8_SWITCHABLE_FILTERS*2-2] = {
  -0, -1,
};
struct vp8_token_struct vp8_switchable_interp_encodings[VP8_SWITCHABLE_FILTERS];
const vp8_prob vp8_switchable_interp_prob [VP8_SWITCHABLE_FILTERS+1]
                                          [VP8_SWITCHABLE_FILTERS-1] = {
  {248},
  { 64},
  {192},
};
const INTERPOLATIONFILTERTYPE vp8_switchable_interp[VP8_SWITCHABLE_FILTERS] = {
  EIGHTTAP, EIGHTTAP_SHARP};
const int vp8_switchable_interp_map[SWITCHABLE+1] = {-1, -1, 0, 1, -1}; //8, 8s
#endif

void vp9_entropy_mode_init() {
  vp9_tokens_from_tree(vp8_bmode_encodings,   vp8_bmode_tree);
  vp9_tokens_from_tree(vp8_ymode_encodings,   vp8_ymode_tree);
  vp9_tokens_from_tree(vp8_kf_ymode_encodings, vp8_kf_ymode_tree);
#if CONFIG_SUPERBLOCKS
  vp9_tokens_from_tree(vp8_sb_kf_ymode_encodings, vp8_sb_ymode_tree);
#endif
  vp9_tokens_from_tree(vp8_uv_mode_encodings,  vp8_uv_mode_tree);
  vp9_tokens_from_tree(vp8_i8x8_mode_encodings,  vp8_i8x8_mode_tree);
  vp9_tokens_from_tree(vp8_mbsplit_encodings, vp8_mbsplit_tree);
  vp9_tokens_from_tree(vp8_switchable_interp_encodings,
                       vp8_switchable_interp_tree);

  vp9_tokens_from_tree_offset(vp8_mv_ref_encoding_array,
                              vp8_mv_ref_tree, NEARESTMV);
#if CONFIG_SUPERBLOCKS
  vp9_tokens_from_tree_offset(vp8_sb_mv_ref_encoding_array,
                              vp8_sb_mv_ref_tree, NEARESTMV);
#endif
  vp9_tokens_from_tree_offset(vp8_sub_mv_ref_encoding_array,
                              vp8_sub_mv_ref_tree, LEFT4X4);
}

void vp9_init_mode_contexts(VP8_COMMON *pc) {
  vpx_memset(pc->fc.mv_ref_ct, 0, sizeof(pc->fc.mv_ref_ct));
  vpx_memset(pc->fc.mv_ref_ct_a, 0, sizeof(pc->fc.mv_ref_ct_a));

  vpx_memcpy(pc->fc.mode_context,
             default_vp8_mode_contexts,
             sizeof(pc->fc.mode_context));
  vpx_memcpy(pc->fc.mode_context_a,
             default_vp8_mode_contexts,
             sizeof(pc->fc.mode_context_a));

}

void vp9_accum_mv_refs(VP8_COMMON *pc,
                       MB_PREDICTION_MODE m,
                       const int ct[4]) {
  int (*mv_ref_ct)[4][2];

  if (pc->refresh_alt_ref_frame)
    mv_ref_ct = pc->fc.mv_ref_ct_a;
  else
    mv_ref_ct = pc->fc.mv_ref_ct;

  if (m == ZEROMV) {
    ++mv_ref_ct [ct[0]] [0] [0];
  } else {
    ++mv_ref_ct [ct[0]] [0] [1];
    if (m == NEARESTMV) {
      ++mv_ref_ct [ct[1]] [1] [0];
    } else {
      ++mv_ref_ct [ct[1]] [1] [1];
      if (m == NEARMV) {
        ++mv_ref_ct [ct[2]] [2] [0];
      } else {
        ++mv_ref_ct [ct[2]] [2] [1];
        if (m == NEWMV) {
          ++mv_ref_ct [ct[3]] [3] [0];
        } else {
          ++mv_ref_ct [ct[3]] [3] [1];
        }
      }
    }
  }
}

#define MVREF_COUNT_SAT 20
#define MVREF_MAX_UPDATE_FACTOR 144
void vp9_update_mode_context(VP8_COMMON *pc) {
  int i, j;
  int (*mv_ref_ct)[4][2];
  int (*mode_context)[4];

  if (pc->refresh_alt_ref_frame) {
    mv_ref_ct = pc->fc.mv_ref_ct_a;
    mode_context = pc->fc.mode_context_a;
  } else {
    mv_ref_ct = pc->fc.mv_ref_ct;
    mode_context = pc->fc.mode_context;
  }

  for (j = 0; j < 6; j++) {
    for (i = 0; i < 4; i++) {
      int this_prob;
      int count = mv_ref_ct[j][i][0] + mv_ref_ct[j][i][1];
      int factor;
      {
        this_prob = count > 0 ? 256 * mv_ref_ct[j][i][0] / count : 128;
        count = count > MVREF_COUNT_SAT ? MVREF_COUNT_SAT : count;
        factor = (MVREF_MAX_UPDATE_FACTOR * count / MVREF_COUNT_SAT);
        this_prob = (pc->fc.vp8_mode_contexts[j][i] * (256 - factor) +
                     this_prob * factor + 128) >> 8;
        this_prob = this_prob ? (this_prob < 255 ? this_prob : 255) : 1;
        mode_context[j][i] = this_prob;
      }
    }
  }
}

#ifdef MODE_STATS
#include "vp8/common/modecont.h"
void print_mode_contexts(VP8_COMMON *pc) {
  int j, i;
  printf("\n====================\n");
  for (j = 0; j < 6; j++) {
    for (i = 0; i < 4; i++) {
      printf("%4d ", pc->fc.mode_context[j][i]);
    }
    printf("\n");
  }
  printf("====================\n");
  for (j = 0; j < 6; j++) {
    for (i = 0; i < 4; i++) {
      printf("%4d ", pc->fc.mode_context_a[j][i]);
    }
    printf("\n");
  }
}
#endif

// #define MODE_COUNT_TESTING
#define MODE_COUNT_SAT 20
#define MODE_MAX_UPDATE_FACTOR 144
void vp9_adapt_mode_probs(VP8_COMMON *cm) {
  int i, t, count, factor;
  unsigned int branch_ct[32][2];
  vp8_prob ymode_probs[VP8_YMODES - 1];
  vp8_prob uvmode_probs[VP8_UV_MODES - 1];
  vp8_prob bmode_probs[VP8_BINTRAMODES - 1];
  vp8_prob i8x8_mode_probs[VP8_I8X8_MODES - 1];
  vp8_prob sub_mv_ref_probs[VP8_SUBMVREFS - 1];
  vp8_prob mbsplit_probs[VP8_NUMMBSPLITS - 1];
#ifdef MODE_COUNT_TESTING
  printf("static const unsigned int\nymode_counts"
         "[VP8_YMODES] = {\n");
  for (t = 0; t < VP8_YMODES; ++t) printf("%d, ", cm->fc.ymode_counts[t]);
  printf("};\n");
  printf("static const unsigned int\nuv_mode_counts"
         "[VP8_YMODES] [VP8_UV_MODES] = {\n");
  for (i = 0; i < VP8_YMODES; ++i) {
    printf("  {");
    for (t = 0; t < VP8_UV_MODES; ++t) printf("%d, ", cm->fc.uv_mode_counts[i][t]);
    printf("},\n");
  }
  printf("};\n");
  printf("static const unsigned int\nbmode_counts"
         "[VP8_BINTRAMODES] = {\n");
  for (t = 0; t < VP8_BINTRAMODES; ++t) printf("%d, ", cm->fc.bmode_counts[t]);
  printf("};\n");
  printf("static const unsigned int\ni8x8_mode_counts"
         "[VP8_I8X8_MODES] = {\n");
  for (t = 0; t < VP8_I8X8_MODES; ++t) printf("%d, ", cm->fc.i8x8_mode_counts[t]);
  printf("};\n");
  printf("static const unsigned int\nsub_mv_ref_counts"
         "[SUBMVREF_COUNT] [VP8_SUBMVREFS] = {\n");
  for (i = 0; i < SUBMVREF_COUNT; ++i) {
    printf("  {");
    for (t = 0; t < VP8_SUBMVREFS; ++t) printf("%d, ", cm->fc.sub_mv_ref_counts[i][t]);
    printf("},\n");
  }
  printf("};\n");
  printf("static const unsigned int\nmbsplit_counts"
         "[VP8_NUMMBSPLITS] = {\n");
  for (t = 0; t < VP8_NUMMBSPLITS; ++t) printf("%d, ", cm->fc.mbsplit_counts[t]);
  printf("};\n");
#endif
  vp9_tree_probs_from_distribution(
    VP8_YMODES, vp8_ymode_encodings, vp8_ymode_tree,
    ymode_probs, branch_ct, cm->fc.ymode_counts,
    256, 1);
  for (t = 0; t < VP8_YMODES - 1; ++t) {
    int prob;
    count = branch_ct[t][0] + branch_ct[t][1];
    count = count > MODE_COUNT_SAT ? MODE_COUNT_SAT : count;
    factor = (MODE_MAX_UPDATE_FACTOR * count / MODE_COUNT_SAT);
    prob = ((int)cm->fc.pre_ymode_prob[t] * (256 - factor) +
            (int)ymode_probs[t] * factor + 128) >> 8;
    if (prob <= 0) cm->fc.ymode_prob[t] = 1;
    else if (prob > 255) cm->fc.ymode_prob[t] = 255;
    else cm->fc.ymode_prob[t] = prob;
  }
  for (i = 0; i < VP8_YMODES; ++i) {
    vp9_tree_probs_from_distribution(
      VP8_UV_MODES, vp8_uv_mode_encodings, vp8_uv_mode_tree,
      uvmode_probs, branch_ct, cm->fc.uv_mode_counts[i],
      256, 1);
    for (t = 0; t < VP8_UV_MODES - 1; ++t) {
      int prob;
      count = branch_ct[t][0] + branch_ct[t][1];
      count = count > MODE_COUNT_SAT ? MODE_COUNT_SAT : count;
      factor = (MODE_MAX_UPDATE_FACTOR * count / MODE_COUNT_SAT);
      prob = ((int)cm->fc.pre_uv_mode_prob[i][t] * (256 - factor) +
              (int)uvmode_probs[t] * factor + 128) >> 8;
      if (prob <= 0) cm->fc.uv_mode_prob[i][t] = 1;
      else if (prob > 255) cm->fc.uv_mode_prob[i][t] = 255;
      else cm->fc.uv_mode_prob[i][t] = prob;
    }
  }
  vp9_tree_probs_from_distribution(
    VP8_BINTRAMODES, vp8_bmode_encodings, vp8_bmode_tree,
    bmode_probs, branch_ct, cm->fc.bmode_counts,
    256, 1);
  for (t = 0; t < VP8_BINTRAMODES - 1; ++t) {
    int prob;
    count = branch_ct[t][0] + branch_ct[t][1];
    count = count > MODE_COUNT_SAT ? MODE_COUNT_SAT : count;
    factor = (MODE_MAX_UPDATE_FACTOR * count / MODE_COUNT_SAT);
    prob = ((int)cm->fc.pre_bmode_prob[t] * (256 - factor) +
            (int)bmode_probs[t] * factor + 128) >> 8;
    if (prob <= 0) cm->fc.bmode_prob[t] = 1;
    else if (prob > 255) cm->fc.bmode_prob[t] = 255;
    else cm->fc.bmode_prob[t] = prob;
  }
  vp9_tree_probs_from_distribution(
    VP8_I8X8_MODES, vp8_i8x8_mode_encodings, vp8_i8x8_mode_tree,
    i8x8_mode_probs, branch_ct, cm->fc.i8x8_mode_counts,
    256, 1);
  for (t = 0; t < VP8_I8X8_MODES - 1; ++t) {
    int prob;
    count = branch_ct[t][0] + branch_ct[t][1];
    count = count > MODE_COUNT_SAT ? MODE_COUNT_SAT : count;
    factor = (MODE_MAX_UPDATE_FACTOR * count / MODE_COUNT_SAT);
    prob = ((int)cm->fc.pre_i8x8_mode_prob[t] * (256 - factor) +
            (int)i8x8_mode_probs[t] * factor + 128) >> 8;
    if (prob <= 0) cm->fc.i8x8_mode_prob[t] = 1;
    else if (prob > 255) cm->fc.i8x8_mode_prob[t] = 255;
    else cm->fc.i8x8_mode_prob[t] = prob;
  }
  for (i = 0; i < SUBMVREF_COUNT; ++i) {
    vp9_tree_probs_from_distribution(
      VP8_SUBMVREFS, vp8_sub_mv_ref_encoding_array, vp8_sub_mv_ref_tree,
      sub_mv_ref_probs, branch_ct, cm->fc.sub_mv_ref_counts[i],
      256, 1);
    for (t = 0; t < VP8_SUBMVREFS - 1; ++t) {
      int prob;
      count = branch_ct[t][0] + branch_ct[t][1];
      count = count > MODE_COUNT_SAT ? MODE_COUNT_SAT : count;
      factor = (MODE_MAX_UPDATE_FACTOR * count / MODE_COUNT_SAT);
      prob = ((int)cm->fc.pre_sub_mv_ref_prob[i][t] * (256 - factor) +
              (int)sub_mv_ref_probs[t] * factor + 128) >> 8;
      if (prob <= 0) cm->fc.sub_mv_ref_prob[i][t] = 1;
      else if (prob > 255) cm->fc.sub_mv_ref_prob[i][t] = 255;
      else cm->fc.sub_mv_ref_prob[i][t] = prob;
    }
  }
  vp9_tree_probs_from_distribution(
    VP8_NUMMBSPLITS, vp8_mbsplit_encodings, vp8_mbsplit_tree,
    mbsplit_probs, branch_ct, cm->fc.mbsplit_counts,
    256, 1);
  for (t = 0; t < VP8_NUMMBSPLITS - 1; ++t) {
    int prob;
    count = branch_ct[t][0] + branch_ct[t][1];
    count = count > MODE_COUNT_SAT ? MODE_COUNT_SAT : count;
    factor = (MODE_MAX_UPDATE_FACTOR * count / MODE_COUNT_SAT);
    prob = ((int)cm->fc.pre_mbsplit_prob[t] * (256 - factor) +
            (int)mbsplit_probs[t] * factor + 128) >> 8;
    if (prob <= 0) cm->fc.mbsplit_prob[t] = 1;
    else if (prob > 255) cm->fc.mbsplit_prob[t] = 255;
    else cm->fc.mbsplit_prob[t] = prob;
  }
}
