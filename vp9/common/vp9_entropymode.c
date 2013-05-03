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

static const unsigned int kf_y_mode_cts[8][VP9_YMODES] = {
#if CONFIG_SB8X8
  /* DC V   H  D45 135 117 153 D27 D63 TM i4X4 */
  {12,  6,  5,  5,  5,  5,  5,  5,  5,  2, 200},
  {25, 13, 13,  7,  7,  7,  7,  7,  7,  6, 160},
  {31, 17, 18,  8,  8,  8,  8,  8,  8,  9, 139},
  {40, 22, 23,  8,  8,  8,  8,  8,  8, 12, 116},
  {53, 26, 28,  8,  8,  8,  8,  8,  8, 13,  94},
  {68, 33, 35,  8,  8,  8,  8,  8,  8, 17,  68},
  {78, 38, 38,  8,  8,  8,  8,  8,  8, 19,  52},
  {89, 42, 42,  8,  8,  8,  8,  8,  8, 21,  34},
#else
  /* DC V   H  D45 135 117 153 D27 D63 TM i8x8 i4X4 */
  {12,  6,  5,  5,  5,  5,  5,  5,  5,  2, 22, 200},
  {25, 13, 13,  7,  7,  7,  7,  7,  7,  6, 27, 160},
  {31, 17, 18,  8,  8,  8,  8,  8,  8,  9, 26, 139},
  {40, 22, 23,  8,  8,  8,  8,  8,  8, 12, 27, 116},
  {53, 26, 28,  8,  8,  8,  8,  8,  8, 13, 26,  94},
  {68, 33, 35,  8,  8,  8,  8,  8,  8, 17, 20,  68},
  {78, 38, 38,  8,  8,  8,  8,  8,  8, 19, 16,  52},
  {89, 42, 42,  8,  8,  8,  8,  8,  8, 21, 12,  34},
#endif
};

static const unsigned int y_mode_cts  [VP9_YMODES] = {
#if CONFIG_SB8X8
  /* DC V   H  D45 135 117 153 D27 D63 TM i4X4 */
  98, 19, 15, 14, 14, 14, 14, 12, 12, 13, 70
#else
  /* DC V   H  D45 135 117 153 D27 D63 TM i8x8 i4X4 */
  98, 19, 15, 14, 14, 14, 14, 12, 12, 13, 16, 70
#endif
};

static const unsigned int uv_mode_cts [VP9_YMODES] [VP9_UV_MODES] = {
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
#if !CONFIG_SB8X8
  { 132, 46, 40, 10, 10, 10, 10, 10, 10, 18}, /* i8x8 - never used */
#endif
  { 150, 35, 41, 10, 10, 10, 10, 10, 10, 10}, /* i4X4 */
};

#if !CONFIG_SB8X8
static const unsigned int i8x8_mode_cts  [VP9_I8X8_MODES] = {
  /* DC V  H D45 135 117 153 D27 D63  TM */
  73, 49, 61, 30, 30, 30, 30, 30, 30, 13
};
#endif

static const unsigned int kf_uv_mode_cts [VP9_YMODES] [VP9_UV_MODES] = {
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
#if !CONFIG_SB8X8
  { 122, 41, 35, 20, 20, 20, 20, 20, 20, 18}, /* i8x8 - never used */
#endif
  { 122, 41, 35, 20, 20, 20, 20, 20, 20, 18}, /* I4X4 */
};

static const unsigned int bmode_cts[VP9_NKF_BINTRAMODES] = {
#if CONFIG_NEWBINTRAMODES
#if CONTEXT_PRED_REPLACEMENTS == 6
  /* DC    TM     V      H   CONTEXT */
  43891, 17694, 10036, 3920, 20000
#elif CONTEXT_PRED_REPLACEMENTS == 4
  /* DC    TM     V      H   D45   D135   CONTEXT */
  43891, 17694, 10036, 3920, 3363, 2546, 14000
#elif CONTEXT_PRED_REPLACEMENTS == 0
  /* DC    V     H    D45   D135  D117  D153   D27   D63   TM    CONTEXT */
  43891, 10036, 3920, 3363, 2546, 5119, 2471, 1723, 3221, 17694, 50000
#endif
#else
  /* DC    V     H    D45   D135  D117  D153   D27   D63   TM  */
  43891, 10036, 3920, 3363, 2546, 5119, 2471, 1723, 3221, 17694
#endif
};

typedef enum {
  SUBMVREF_NORMAL,
  SUBMVREF_LEFT_ZED,
  SUBMVREF_ABOVE_ZED,
  SUBMVREF_LEFT_ABOVE_SAME,
  SUBMVREF_LEFT_ABOVE_ZED
} sumvfref_t;

int vp9_mv_cont(const int_mv *l, const int_mv *a) {
  const int lez = (l->as_int == 0);
  const int aez = (a->as_int == 0);
  const int lea = (l->as_int == a->as_int);

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

const vp9_prob vp9_sub_mv_ref_prob2 [SUBMVREF_COUNT][VP9_SUBMVREFS - 1] = {
  { 147, 136, 18 },
  { 106, 145, 1  },
  { 179, 121, 1  },
  { 223, 1, 34 },
  { 208, 1, 1  }
};

#if !CONFIG_SB8X8
vp9_mbsplit vp9_mbsplits [VP9_NUMMBSPLITS] = {
  {
    0,  0,  0,  0,
    0,  0,  0,  0,
    1,  1,  1,  1,
    1,  1,  1,  1,
  }, {
    0,  0,  1,  1,
    0,  0,  1,  1,
    0,  0,  1,  1,
    0,  0,  1,  1,
  }, {
    0,  0,  1,  1,
    0,  0,  1,  1,
    2,  2,  3,  3,
    2,  2,  3,  3,
  }, {
    0,  1,  2,  3,
    4,  5,  6,  7,
    8,  9,  10, 11,
    12, 13, 14, 15,
  },
};

const int vp9_mbsplit_count [VP9_NUMMBSPLITS] = { 2, 2, 4, 16};

const vp9_prob vp9_mbsplit_probs [VP9_NUMMBSPLITS - 1] = { 110, 111, 150};
#endif

const vp9_prob vp9_partition_probs[NUM_PARTITION_CONTEXTS]
                                  [PARTITION_TYPES - 1] = {
#if CONFIG_SB8X8
  // FIXME(jingning,rbultje) put real probabilities here
  {202, 162, 107},
  {16,  2,   169},
  {3,   246,  19},
  {104, 90,  134},
#endif
  {202, 162, 107},
  {16,  2,   169},
  {3,   246,  19},
  {104, 90,  134},
  {183, 70,  109},
  {30,  14,  162},
  {67,  208,  22},
  {4,   17,   5},
};

/* Array indices are identical to previously-existing INTRAMODECONTEXTNODES. */

const vp9_tree_index vp9_kf_bmode_tree[VP9_KF_BINTRAMODES * 2 - 2] = {
  -B_DC_PRED, 2,                      /* 0 = DC_NODE */
  -B_TM_PRED, 4,                      /* 1 = TM_NODE */
  -B_V_PRED, 6,                       /* 2 = V_NODE */
  8, 12,                              /* 3 = COM_NODE */
  -B_H_PRED, 10,                      /* 4 = H_NODE */
  -B_D135_PRED, -B_D117_PRED,         /* 5 = D135_NODE */
  -B_D45_PRED, 14,                    /* 6 = D45_NODE */
  -B_D63_PRED, 16,                    /* 7 = D63_NODE */
  -B_D153_PRED, -B_D27_PRED           /* 8 = D153_NODE */
};

const vp9_tree_index vp9_bmode_tree[VP9_NKF_BINTRAMODES * 2 - 2] = {
#if CONFIG_NEWBINTRAMODES
#if CONTEXT_PRED_REPLACEMENTS == 6
  -B_DC_PRED, 2,
  -B_TM_PRED, 4,
  6, -(B_CONTEXT_PRED - CONTEXT_PRED_REPLACEMENTS),
  -B_V_PRED, -B_H_PRED
#elif CONTEXT_PRED_REPLACEMENTS == 4
  -B_DC_PRED, 2,
  -B_TM_PRED, 4,
  6, 8,
  -B_V_PRED, -B_H_PRED,
  10, -(B_CONTEXT_PRED - CONTEXT_PRED_REPLACEMENTS),
  -B_D135_PRED, -B_D45_PRED,
#elif CONTEXT_PRED_REPLACEMENTS == 0
  -B_DC_PRED, 2,                      /* 0 = DC_NODE */
  -B_TM_PRED, 4,                      /* 1 = TM_NODE */
  -B_V_PRED, 6,                       /* 2 = V_NODE */
  8, 12,                              /* 3 = COM_NODE */
  -B_H_PRED, 10,                      /* 4 = H_NODE */
  -B_D135_PRED, -B_D117_PRED,         /* 5 = D135_NODE */
  -B_D45_PRED, 14,                    /* 6 = D45_NODE */
  -B_D63_PRED, 16,                    /* 7 = D63_NODE */
  -B_D153_PRED, 18,                   /* 8 = D153_NODE */
  -B_D27_PRED, -B_CONTEXT_PRED        /* 9 = D27_NODE */
#endif
#else
  -B_DC_PRED, 2,                      /* 0 = DC_NODE */
  -B_TM_PRED, 4,                      /* 1 = TM_NODE */
  -B_V_PRED, 6,                       /* 2 = V_NODE */
  8, 12,                              /* 3 = COM_NODE */
  -B_H_PRED, 10,                      /* 4 = H_NODE */
  -B_D135_PRED, -B_D117_PRED,         /* 5 = D135_NODE */
  -B_D45_PRED, 14,                    /* 6 = D45_NODE */
  -B_D63_PRED, 16,                    /* 7 = D63_NODE */
  -B_D153_PRED, -B_D27_PRED           /* 8 = D153_NODE */
#endif
};

/* Again, these trees use the same probability indices as their
   explicitly-programmed predecessors. */
const vp9_tree_index vp9_ymode_tree[VP9_YMODES * 2 - 2] = {
  2, 14,
  -DC_PRED, 4,
  6, 8,
  -D45_PRED, -D135_PRED,
  10, 12,
  -D117_PRED, -D153_PRED,
  -D27_PRED, -D63_PRED,
  16, 18,
  -V_PRED, -H_PRED,
#if CONFIG_SB8X8
  -TM_PRED, -I4X4_PRED
#else
  -TM_PRED, 20,
  -I4X4_PRED, -I8X8_PRED
#endif
};

const vp9_tree_index vp9_kf_ymode_tree[VP9_YMODES * 2 - 2] = {
  2, 14,
  -DC_PRED, 4,
  6, 8,
  -D45_PRED, -D135_PRED,
  10, 12,
  -D117_PRED, -D153_PRED,
  -D27_PRED, -D63_PRED,
  16, 18,
  -V_PRED, -H_PRED,
#if CONFIG_SB8X8
  -TM_PRED, -I4X4_PRED
#else
  -TM_PRED, 20,
  -I4X4_PRED, -I8X8_PRED
#endif
};

#if !CONFIG_SB8X8
const vp9_tree_index vp9_i8x8_mode_tree[VP9_I8X8_MODES * 2 - 2] = {
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
#endif

const vp9_tree_index vp9_uv_mode_tree[VP9_UV_MODES * 2 - 2] = {
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

#if !CONFIG_SB8X8
const vp9_tree_index vp9_mbsplit_tree[6] = {
  -PARTITIONING_4X4,   2,
  -PARTITIONING_8X8,   4,
  -PARTITIONING_16X8, -PARTITIONING_8X16,
};
#endif

const vp9_tree_index vp9_mv_ref_tree[8] = {
  -ZEROMV, 2,
  -NEARESTMV, 4,
  -NEARMV, 6,
  -NEWMV, -SPLITMV
};

const vp9_tree_index vp9_sb_mv_ref_tree[6] = {
  -ZEROMV, 2,
  -NEARESTMV, 4,
  -NEARMV, -NEWMV
};

const vp9_tree_index vp9_sub_mv_ref_tree[6] = {
  -LEFT4X4, 2,
  -ABOVE4X4, 4,
  -ZERO4X4, -NEW4X4
};

const vp9_tree_index vp9_partition_tree[6] = {
  -PARTITION_NONE, 2,
  -PARTITION_HORZ, 4,
  -PARTITION_VERT, -PARTITION_SPLIT
};

struct vp9_token vp9_bmode_encodings[VP9_NKF_BINTRAMODES];
struct vp9_token vp9_kf_bmode_encodings[VP9_KF_BINTRAMODES];
struct vp9_token vp9_ymode_encodings[VP9_YMODES];
struct vp9_token vp9_sb_ymode_encodings[VP9_I32X32_MODES];
struct vp9_token vp9_sb_kf_ymode_encodings[VP9_I32X32_MODES];
struct vp9_token vp9_kf_ymode_encodings[VP9_YMODES];
struct vp9_token vp9_uv_mode_encodings[VP9_UV_MODES];
#if !CONFIG_SB8X8
struct vp9_token vp9_i8x8_mode_encodings[VP9_I8X8_MODES];
struct vp9_token vp9_mbsplit_encodings[VP9_NUMMBSPLITS];
#endif

struct vp9_token vp9_mv_ref_encoding_array[VP9_MVREFS];
struct vp9_token vp9_sb_mv_ref_encoding_array[VP9_MVREFS];
struct vp9_token vp9_sub_mv_ref_encoding_array[VP9_SUBMVREFS];

struct vp9_token vp9_partition_encodings[PARTITION_TYPES];

void vp9_init_mbmode_probs(VP9_COMMON *x) {
  unsigned int bct[VP9_YMODES][2];  // num Ymodes > num UV modes
  int i;

  vp9_tree_probs_from_distribution(vp9_ymode_tree, x->fc.ymode_prob,
                                   bct, y_mode_cts, 0);
  vp9_tree_probs_from_distribution(vp9_sb_ymode_tree, x->fc.sb_ymode_prob,
                                   bct, y_mode_cts, 0);
  for (i = 0; i < 8; i++) {
    vp9_tree_probs_from_distribution(vp9_kf_ymode_tree, x->kf_ymode_prob[i],
                                     bct, kf_y_mode_cts[i], 0);
    vp9_tree_probs_from_distribution(vp9_sb_kf_ymode_tree,
                                     x->sb_kf_ymode_prob[i], bct,
                                     kf_y_mode_cts[i], 0);
  }

  for (i = 0; i < VP9_YMODES; i++) {
    vp9_tree_probs_from_distribution(vp9_uv_mode_tree, x->kf_uv_mode_prob[i],
                                     bct, kf_uv_mode_cts[i], 0);
    vp9_tree_probs_from_distribution(vp9_uv_mode_tree, x->fc.uv_mode_prob[i],
                                     bct, uv_mode_cts[i], 0);
  }

#if !CONFIG_SB8X8
  vp9_tree_probs_from_distribution(vp9_i8x8_mode_tree, x->fc.i8x8_mode_prob,
                                   bct, i8x8_mode_cts, 0);
#endif

  vpx_memcpy(x->fc.sub_mv_ref_prob, vp9_sub_mv_ref_prob2,
             sizeof(vp9_sub_mv_ref_prob2));
#if !CONFIG_SB8X8
  vpx_memcpy(x->fc.mbsplit_prob, vp9_mbsplit_probs, sizeof(vp9_mbsplit_probs));
#endif
  vpx_memcpy(x->fc.switchable_interp_prob, vp9_switchable_interp_prob,
             sizeof(vp9_switchable_interp_prob));

  vpx_memcpy(x->fc.partition_prob, vp9_partition_probs,
             sizeof(vp9_partition_probs));

#if CONFIG_COMP_INTERINTRA_PRED
  x->fc.interintra_prob = VP9_DEF_INTERINTRA_PROB;
#endif
  x->ref_pred_probs[0] = DEFAULT_PRED_PROB_0;
  x->ref_pred_probs[1] = DEFAULT_PRED_PROB_1;
  x->ref_pred_probs[2] = DEFAULT_PRED_PROB_2;
}


static void intra_bmode_probs_from_distribution(
  vp9_prob p[VP9_NKF_BINTRAMODES - 1],
  unsigned int branch_ct[VP9_NKF_BINTRAMODES - 1][2],
  const unsigned int events[VP9_NKF_BINTRAMODES]) {
  vp9_tree_probs_from_distribution(vp9_bmode_tree, p, branch_ct, events, 0);
}

void vp9_default_bmode_probs(vp9_prob p[VP9_NKF_BINTRAMODES - 1]) {
  unsigned int branch_ct[VP9_NKF_BINTRAMODES - 1][2];
  intra_bmode_probs_from_distribution(p, branch_ct, bmode_cts);
}

static void intra_kf_bmode_probs_from_distribution(
  vp9_prob p[VP9_KF_BINTRAMODES - 1],
  unsigned int branch_ct[VP9_KF_BINTRAMODES - 1][2],
  const unsigned int events[VP9_KF_BINTRAMODES]) {
  vp9_tree_probs_from_distribution(vp9_kf_bmode_tree, p, branch_ct, events, 0);
}

void vp9_kf_default_bmode_probs(vp9_prob p[VP9_KF_BINTRAMODES]
                                          [VP9_KF_BINTRAMODES]
                                          [VP9_KF_BINTRAMODES - 1]) {
  unsigned int branch_ct[VP9_KF_BINTRAMODES - 1][2];
  int i, j;

  for (i = 0; i < VP9_KF_BINTRAMODES; ++i) {
    for (j = 0; j < VP9_KF_BINTRAMODES; ++j) {
      intra_kf_bmode_probs_from_distribution(
          p[i][j], branch_ct, vp9_kf_default_bmode_counts[i][j]);
    }
  }
}

#if VP9_SWITCHABLE_FILTERS == 3
const vp9_tree_index vp9_switchable_interp_tree[VP9_SWITCHABLE_FILTERS*2-2] = {
  -0, 2,
  -1, -2
};
struct vp9_token vp9_switchable_interp_encodings[VP9_SWITCHABLE_FILTERS];
#if CONFIG_ENABLE_6TAP
const INTERPOLATIONFILTERTYPE vp9_switchable_interp[VP9_SWITCHABLE_FILTERS] = {
  SIXTAP, EIGHTTAP, EIGHTTAP_SHARP};
const int vp9_switchable_interp_map[SWITCHABLE+1] = {0, -1, 1, 2, -1, -1};
#else
const INTERPOLATIONFILTERTYPE vp9_switchable_interp[VP9_SWITCHABLE_FILTERS] = {
  EIGHTTAP, EIGHTTAP_SMOOTH, EIGHTTAP_SHARP};
const int vp9_switchable_interp_map[SWITCHABLE+1] = {1, 0, 2, -1, -1};
#endif
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
#if CONFIG_ENABLE_6TAP
const int vp9_switchable_interp_map[SWITCHABLE+1] = {-1, -1, 0, 1, -1, -1};
#else
const int vp9_switchable_interp_map[SWITCHABLE+1] = {-1, 0, 1, -1, -1};
#endif
#endif  // VP9_SWITCHABLE_FILTERS

// Indicates if the filter is interpolating or non-interpolating
// Note currently only the EIGHTTAP_SMOOTH is non-interpolating
#if CONFIG_ENABLE_6TAP
const int vp9_is_interpolating_filter[SWITCHABLE + 1] = {1, 0, 1, 1, 1, -1};
#else
const int vp9_is_interpolating_filter[SWITCHABLE + 1] = {0, 1, 1, 1, -1};
#endif

void vp9_entropy_mode_init() {
  vp9_tokens_from_tree(vp9_kf_bmode_encodings,   vp9_kf_bmode_tree);
  vp9_tokens_from_tree(vp9_bmode_encodings,   vp9_bmode_tree);
  vp9_tokens_from_tree(vp9_ymode_encodings,   vp9_ymode_tree);
  vp9_tokens_from_tree(vp9_kf_ymode_encodings, vp9_kf_ymode_tree);
  vp9_tokens_from_tree(vp9_sb_ymode_encodings, vp9_sb_ymode_tree);
  vp9_tokens_from_tree(vp9_sb_kf_ymode_encodings, vp9_sb_kf_ymode_tree);
  vp9_tokens_from_tree(vp9_uv_mode_encodings,  vp9_uv_mode_tree);
#if !CONFIG_SB8X8
  vp9_tokens_from_tree(vp9_i8x8_mode_encodings,  vp9_i8x8_mode_tree);
  vp9_tokens_from_tree(vp9_mbsplit_encodings, vp9_mbsplit_tree);
#endif
  vp9_tokens_from_tree(vp9_switchable_interp_encodings,
                       vp9_switchable_interp_tree);
  vp9_tokens_from_tree(vp9_partition_encodings, vp9_partition_tree);

  vp9_tokens_from_tree_offset(vp9_mv_ref_encoding_array,
                              vp9_mv_ref_tree, NEARESTMV);
  vp9_tokens_from_tree_offset(vp9_sb_mv_ref_encoding_array,
                              vp9_sb_mv_ref_tree, NEARESTMV);
  vp9_tokens_from_tree_offset(vp9_sub_mv_ref_encoding_array,
                              vp9_sub_mv_ref_tree, LEFT4X4);
}

void vp9_init_mode_contexts(VP9_COMMON *pc) {
  vpx_memset(pc->fc.mv_ref_ct, 0, sizeof(pc->fc.mv_ref_ct));
  vpx_memcpy(pc->fc.vp9_mode_contexts,
             vp9_default_mode_contexts,
             sizeof(vp9_default_mode_contexts));
}

void vp9_accum_mv_refs(VP9_COMMON *pc,
                       MB_PREDICTION_MODE m,
                       const int context) {
  unsigned int (*mv_ref_ct)[4][2] = pc->fc.mv_ref_ct;

  if (m == ZEROMV) {
    ++mv_ref_ct[context][0][0];
  } else {
    ++mv_ref_ct[context][0][1];
    if (m == NEARESTMV) {
      ++mv_ref_ct[context][1][0];
    } else {
      ++mv_ref_ct[context][1][1];
      if (m == NEARMV) {
        ++mv_ref_ct[context][2][0];
      } else {
        ++mv_ref_ct[context][2][1];
        if (m == NEWMV) {
          ++mv_ref_ct[context][3][0];
        } else {
          ++mv_ref_ct[context][3][1];
        }
      }
    }
  }
}

#define MVREF_COUNT_SAT 20
#define MVREF_MAX_UPDATE_FACTOR 128
void vp9_adapt_mode_context(VP9_COMMON *pc) {
  int i, j;
  unsigned int (*mv_ref_ct)[4][2] = pc->fc.mv_ref_ct;
  int (*mode_context)[4] = pc->fc.vp9_mode_contexts;

  for (j = 0; j < INTER_MODE_CONTEXTS; j++) {
    for (i = 0; i < 4; i++) {
      int count = mv_ref_ct[j][i][0] + mv_ref_ct[j][i][1], factor;

      count = count > MVREF_COUNT_SAT ? MVREF_COUNT_SAT : count;
      factor = (MVREF_MAX_UPDATE_FACTOR * count / MVREF_COUNT_SAT);
      mode_context[j][i] = weighted_prob(pc->fc.vp9_mode_contexts[j][i],
                                         get_binary_prob(mv_ref_ct[j][i][0],
                                                         mv_ref_ct[j][i][1]),
                                         factor);
    }
  }
}

#ifdef MODE_STATS
#include "vp9/common/vp9_modecont.h"
void print_mode_contexts(VP9_COMMON *pc) {
  int j, i;
  printf("\n====================\n");
  for (j = 0; j < INTER_MODE_CONTEXTS; j++) {
    for (i = 0; i < 4; i++) {
      printf("%4d ", pc->fc.mode_context[j][i]);
    }
    printf("\n");
  }
  printf("====================\n");
  for (j = 0; j < INTER_MODE_CONTEXTS; j++) {
    for (i = 0; i < 4; i++) {
      printf("%4d ", pc->fc.mode_context_a[j][i]);
    }
    printf("\n");
  }
}
#endif

#define MODE_COUNT_SAT 20
#define MODE_MAX_UPDATE_FACTOR 144
static void update_mode_probs(int n_modes,
                              const vp9_tree_index *tree, unsigned int *cnt,
                              vp9_prob *pre_probs, vp9_prob *dst_probs,
                              unsigned int tok0_offset) {
#define MAX_PROBS 32
  vp9_prob probs[MAX_PROBS];
  unsigned int branch_ct[MAX_PROBS][2];
  int t, count, factor;

  assert(n_modes - 1 < MAX_PROBS);
  vp9_tree_probs_from_distribution(tree, probs, branch_ct, cnt, tok0_offset);
  for (t = 0; t < n_modes - 1; ++t) {
    count = branch_ct[t][0] + branch_ct[t][1];
    count = count > MODE_COUNT_SAT ? MODE_COUNT_SAT : count;
    factor = (MODE_MAX_UPDATE_FACTOR * count / MODE_COUNT_SAT);
    dst_probs[t] = weighted_prob(pre_probs[t], probs[t], factor);
  }
}

// #define MODE_COUNT_TESTING
void vp9_adapt_mode_probs(VP9_COMMON *cm) {
  int i;
  FRAME_CONTEXT *fc = &cm->fc;
#ifdef MODE_COUNT_TESTING
  int t;

  printf("static const unsigned int\nymode_counts"
         "[VP9_YMODES] = {\n");
  for (t = 0; t < VP9_YMODES; ++t)
    printf("%d, ", fc->ymode_counts[t]);
  printf("};\n");
  printf("static const unsigned int\nuv_mode_counts"
         "[VP9_YMODES] [VP9_UV_MODES] = {\n");
  for (i = 0; i < VP9_YMODES; ++i) {
    printf("  {");
    for (t = 0; t < VP9_UV_MODES; ++t)
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
  printf("static const unsigned int\nsub_mv_ref_counts"
         "[SUBMVREF_COUNT] [VP9_SUBMVREFS] = {\n");
  for (i = 0; i < SUBMVREF_COUNT; ++i) {
    printf("  {");
    for (t = 0; t < VP9_SUBMVREFS; ++t)
      printf("%d, ", fc->sub_mv_ref_counts[i][t]);
    printf("},\n");
  }
  printf("};\n");
  printf("static const unsigned int\nmbsplit_counts"
         "[VP9_NUMMBSPLITS] = {\n");
  for (t = 0; t < VP9_NUMMBSPLITS; ++t)
    printf("%d, ", fc->mbsplit_counts[t]);
  printf("};\n");
#if CONFIG_COMP_INTERINTRA_PRED
  printf("static const unsigned int\ninterintra_counts"
         "[2] = {\n");
  for (t = 0; t < 2; ++t)
    printf("%d, ", fc->interintra_counts[t]);
  printf("};\n");
#endif
#endif

  update_mode_probs(VP9_YMODES, vp9_ymode_tree,
                    fc->ymode_counts, fc->pre_ymode_prob,
                    fc->ymode_prob, 0);
  update_mode_probs(VP9_I32X32_MODES, vp9_sb_ymode_tree,
                    fc->sb_ymode_counts, fc->pre_sb_ymode_prob,
                    fc->sb_ymode_prob, 0);

  for (i = 0; i < VP9_YMODES; ++i)
    update_mode_probs(VP9_UV_MODES, vp9_uv_mode_tree,
                      fc->uv_mode_counts[i], fc->pre_uv_mode_prob[i],
                      fc->uv_mode_prob[i], 0);

  update_mode_probs(VP9_NKF_BINTRAMODES, vp9_bmode_tree,
                    fc->bmode_counts, fc->pre_bmode_prob,
                    fc->bmode_prob, 0);
#if !CONFIG_SB8X8
  update_mode_probs(VP9_I8X8_MODES,
                    vp9_i8x8_mode_tree, fc->i8x8_mode_counts,
                    fc->pre_i8x8_mode_prob, fc->i8x8_mode_prob, 0);
#endif

  for (i = 0; i < SUBMVREF_COUNT; ++i)
    update_mode_probs(VP9_SUBMVREFS,
                      vp9_sub_mv_ref_tree, fc->sub_mv_ref_counts[i],
                      fc->pre_sub_mv_ref_prob[i], fc->sub_mv_ref_prob[i],
                      LEFT4X4);

#if !CONFIG_SB8X8
  update_mode_probs(VP9_NUMMBSPLITS, vp9_mbsplit_tree,
                    fc->mbsplit_counts, fc->pre_mbsplit_prob,
                    fc->mbsplit_prob, 0);
#endif
#if CONFIG_COMP_INTERINTRA_PRED
  if (cm->use_interintra) {
    int factor, interintra_prob, count;

    interintra_prob = get_binary_prob(fc->interintra_counts[0],
                                      fc->interintra_counts[1]);
    count = fc->interintra_counts[0] + fc->interintra_counts[1];
    count = count > MODE_COUNT_SAT ? MODE_COUNT_SAT : count;
    factor = (MODE_MAX_UPDATE_FACTOR * count / MODE_COUNT_SAT);
    fc->interintra_prob = weighted_prob(fc->pre_interintra_prob,
                                        interintra_prob, factor);
  }
#endif
  for (i = 0; i < NUM_PARTITION_CONTEXTS; i++)
    update_mode_probs(PARTITION_TYPES, vp9_partition_tree,
                      fc->partition_counts[i], fc->pre_partition_prob[i],
                      fc->partition_prob[i], 0);
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
  vp9_default_bmode_probs(cm->fc.bmode_prob);
  vp9_kf_default_bmode_probs(cm->kf_bmode_prob);
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
