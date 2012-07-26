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
#include "entropymv.h"

//#define MV_COUNT_TESTING

#if CONFIG_NEWMVENTROPY

#define MV_COUNT_SAT 16
#define MV_MAX_UPDATE_FACTOR 160

/* Smooth or bias the mv-counts before prob computation */
/* #define SMOOTH_MV_COUNTS */

const vp8_tree_index vp8_mv_joint_tree[2 * MV_JOINTS - 2] = {
  -MV_JOINT_ZERO, 2,
  -MV_JOINT_HNZVZ, 4,
  -MV_JOINT_HZVNZ, -MV_JOINT_HNZVNZ
};
struct vp8_token_struct vp8_mv_joint_encodings[MV_JOINTS];

const vp8_tree_index vp8_mv_class_tree[2 * MV_CLASSES - 2] = {
  -MV_CLASS_0, 2,
  -MV_CLASS_1, 4,
  6, 8,
  -MV_CLASS_2, -MV_CLASS_3,
  10, 12,
  -MV_CLASS_4, -MV_CLASS_5,
  -MV_CLASS_6, -MV_CLASS_7,
};
struct vp8_token_struct vp8_mv_class_encodings[MV_CLASSES];

const vp8_tree_index vp8_mv_class0_tree [2 * CLASS0_SIZE - 2] = {
  -0, -1,
};
struct vp8_token_struct vp8_mv_class0_encodings[CLASS0_SIZE];

const vp8_tree_index vp8_mv_fp_tree [2 * 4 - 2] = {
  -0, 2,
  -1, 4,
  -2, -3
};
struct vp8_token_struct vp8_mv_fp_encodings[4];

const nmv_context vp8_default_nmv_context = {
  {32, 64, 96},
  {
    { /* vert component */
      128,                                             /* sign */
      {224, 144, 192, 168, 192, 176, 192},             /* class */
      {216},                                           /* class0 */
      {136, 140, 148, 160, 176, 192, 224},             /* bits */
      {{128, 128, 64}, {96, 112, 64}},                 /* class0_fp */
      {64, 96, 64},                                    /* fp */
      160,                                             /* class0_hp bit */
      128,                                             /* hp */
    },
    { /* hor component */
      128,                                             /* sign */
      {216, 128, 176, 160, 176, 176, 192},             /* class */
      {208},                                           /* class0 */
      {136, 140, 148, 160, 176, 192, 224},             /* bits */
      {{128, 128, 64}, {96, 112, 64}},                 /* class0_fp */
      {64, 96, 64},                                    /* fp */
      160,                                             /* class0_hp bit */
      128,                                             /* hp */
    }
  },
};

MV_JOINT_TYPE vp8_get_mv_joint(MV mv) {
  if (mv.row == 0 && mv.col == 0) return MV_JOINT_ZERO;
  else if (mv.row == 0 && mv.col != 0) return MV_JOINT_HNZVZ;
  else if (mv.row != 0 && mv.col == 0) return MV_JOINT_HZVNZ;
  else return MV_JOINT_HNZVNZ;
}

#define mv_class_base(c) ((c) ? (CLASS0_SIZE << (c + 2)) : 0)

MV_CLASS_TYPE vp8_get_mv_class(int z, int *offset) {
  MV_CLASS_TYPE c;
  if      (z < CLASS0_SIZE * 8)    c = MV_CLASS_0;
  else if (z < CLASS0_SIZE * 16)   c = MV_CLASS_1;
  else if (z < CLASS0_SIZE * 32)   c = MV_CLASS_2;
  else if (z < CLASS0_SIZE * 64)   c = MV_CLASS_3;
  else if (z < CLASS0_SIZE * 128)  c = MV_CLASS_4;
  else if (z < CLASS0_SIZE * 256)  c = MV_CLASS_5;
  else if (z < CLASS0_SIZE * 512)  c = MV_CLASS_6;
  else if (z < CLASS0_SIZE * 1024) c = MV_CLASS_7;
  else assert(0);
  if (offset)
    *offset = z - mv_class_base(c);
  return c;
}

int vp8_get_mv_mag(MV_CLASS_TYPE c, int offset) {
  return mv_class_base(c) + offset;
}

static void increment_nmv_component_count(int v,
                                          nmv_component_counts *mvcomp,
                                          int incr,
                                          int usehp) {
  assert (v != 0);            /* should not be zero */
  mvcomp->mvcount[MV_MAX + v] += incr;
}

static void increment_nmv_component(int v,
                                    nmv_component_counts *mvcomp,
                                    int incr,
                                    int usehp) {
  int s, z, c, o, d, e, f;
  assert (v != 0);            /* should not be zero */
  s = v < 0;
  mvcomp->sign[s] += incr;
  z = (s ? -v : v) - 1;       /* magnitude - 1 */

  c = vp8_get_mv_class(z, &o);
  mvcomp->classes[c] += incr;

  d = (o >> 3);               /* int mv data */
  f = (o >> 1) & 3;           /* fractional pel mv data */
  e = (o & 1);                /* high precision mv data */
  if (c == MV_CLASS_0) {
    mvcomp->class0[d] += incr;
  } else {
    int i, b;
    b = c + CLASS0_BITS - 1;  /* number of bits */
    for (i = 0; i < b; ++i)
      mvcomp->bits[i][((d >> i) & 1)] += incr;
  }

  /* Code the fractional pel bits */
  if (c == MV_CLASS_0) {
    mvcomp->class0_fp[d][f] += incr;
  } else {
    mvcomp->fp[f] += incr;
  }

  /* Code the high precision bit */
  if (usehp) {
    if (c == MV_CLASS_0) {
      mvcomp->class0_hp[e] += incr;
    } else {
      mvcomp->hp[e] += incr;
    }
  } else {  /* assume the extra bit is 1 */
    if (c == MV_CLASS_0) {
      mvcomp->class0_hp[1] += incr;
    } else {
      mvcomp->hp[1] += incr;
    }
  }
}

#ifdef SMOOTH_MV_COUNTS
static void smooth_counts(nmv_component_counts *mvcomp) {
  static const int flen = 3;  // (filter_length + 1) / 2
  static const int fval[] = {8, 3, 1};
  static const int fvalbits = 4;
  int i;
  unsigned int smvcount[MV_VALS];
  vpx_memcpy(smvcount, mvcomp->mvcount, sizeof(smvcount));
  smvcount[MV_MAX] = (smvcount[MV_MAX - 1] + smvcount[MV_MAX + 1]) >> 1;
  for (i = flen - 1; i <= MV_VALS - flen; ++i) {
    int j, s = smvcount[i] * fval[0];
    for (j = 1; j < flen; ++j)
      s += (smvcount[i - j] + smvcount[i + j]) * fval[j];
    mvcomp->mvcount[i] = (s + (1 << (fvalbits - 1))) >> fvalbits;
  }
}
#endif

static void counts_to_context(nmv_component_counts *mvcomp, int usehp) {
  int v;
  vpx_memset(mvcomp->sign, 0, sizeof(nmv_component_counts) - sizeof(mvcomp->mvcount));
  for (v = 1; v <= MV_MAX; v++) {
    increment_nmv_component(-v, mvcomp, mvcomp->mvcount[MV_MAX - v], usehp);
    increment_nmv_component( v, mvcomp, mvcomp->mvcount[MV_MAX + v], usehp);
  }
}

void vp8_increment_nmv(const MV *mv, const MV *ref, nmv_context_counts *mvctx,
                       int usehp) {
  MV_JOINT_TYPE j = vp8_get_mv_joint(*mv);
  mvctx->joints[j]++;
  if (j == MV_JOINT_HZVNZ || j == MV_JOINT_HNZVNZ) {
    increment_nmv_component_count(mv->row, &mvctx->comps[0], 1, usehp);
  }
  if (j == MV_JOINT_HNZVZ || j == MV_JOINT_HNZVNZ) {
    increment_nmv_component_count(mv->col, &mvctx->comps[1], 1, usehp);
  }
}

static void adapt_prob(vp8_prob *dest, vp8_prob prep, vp8_prob newp,
                       unsigned int ct[2]) {
  int factor;
  int prob;
  int count = ct[0] + ct[1];
  if (count) {
    count = count > MV_COUNT_SAT ? MV_COUNT_SAT : count;
    factor = (MV_MAX_UPDATE_FACTOR * count / MV_COUNT_SAT);
    prob = ((int)prep * (256 - factor) + (int)(newp) * factor + 128) >> 8;
    prob += !prob;
    prob = (prob > 255 ? 255 : prob);
    *dest = prob;
  }
}

void vp8_counts_to_nmv_context(
    nmv_context_counts *NMVcount,
    nmv_context *prob,
    int usehp,
    unsigned int (*branch_ct_joint)[2],
    unsigned int (*branch_ct_sign)[2],
    unsigned int (*branch_ct_classes)[MV_CLASSES - 1][2],
    unsigned int (*branch_ct_class0)[CLASS0_SIZE - 1][2],
    unsigned int (*branch_ct_bits)[MV_OFFSET_BITS][2],
    unsigned int (*branch_ct_class0_fp)[CLASS0_SIZE][4 - 1][2],
    unsigned int (*branch_ct_fp)[4 - 1][2],
    unsigned int (*branch_ct_class0_hp)[2],
    unsigned int (*branch_ct_hp)[2]) {
  int i, j, k;
  counts_to_context(&NMVcount->comps[0], usehp);
  counts_to_context(&NMVcount->comps[1], usehp);
  vp8_tree_probs_from_distribution(MV_JOINTS,
                                   vp8_mv_joint_encodings,
                                   vp8_mv_joint_tree,
                                   prob->joints,
                                   branch_ct_joint,
                                   NMVcount->joints,
                                   256, 1);
  for (i = 0; i < 2; ++i) {
    prob->comps[i].sign =
        vp8_bin_prob_from_distribution(NMVcount->comps[i].sign);
    branch_ct_sign[i][0] = NMVcount->comps[i].sign[0];
    branch_ct_sign[i][1] = NMVcount->comps[i].sign[1];
    vp8_tree_probs_from_distribution(MV_CLASSES,
                                     vp8_mv_class_encodings,
                                     vp8_mv_class_tree,
                                     prob->comps[i].classes,
                                     branch_ct_classes[i],
                                     NMVcount->comps[i].classes,
                                     256, 1);
    vp8_tree_probs_from_distribution(CLASS0_SIZE,
                                     vp8_mv_class0_encodings,
                                     vp8_mv_class0_tree,
                                     prob->comps[i].class0,
                                     branch_ct_class0[i],
                                     NMVcount->comps[i].class0,
                                     256, 1);
    for (j = 0; j < MV_OFFSET_BITS; ++j) {
      prob->comps[i].bits[j] = vp8_bin_prob_from_distribution(
          NMVcount->comps[i].bits[j]);
      branch_ct_bits[i][j][0] = NMVcount->comps[i].bits[j][0];
      branch_ct_bits[i][j][1] = NMVcount->comps[i].bits[j][1];
    }
  }
  for (i = 0; i < 2; ++i) {
    for (k = 0; k < CLASS0_SIZE; ++k) {
      vp8_tree_probs_from_distribution(4,
                                       vp8_mv_fp_encodings,
                                       vp8_mv_fp_tree,
                                       prob->comps[i].class0_fp[k],
                                       branch_ct_class0_fp[i][k],
                                       NMVcount->comps[i].class0_fp[k],
                                       256, 1);
    }
    vp8_tree_probs_from_distribution(4,
                                     vp8_mv_fp_encodings,
                                     vp8_mv_fp_tree,
                                     prob->comps[i].fp,
                                     branch_ct_fp[i],
                                     NMVcount->comps[i].fp,
                                     256, 1);
  }
  if (usehp) {
    for (i = 0; i < 2; ++i) {
      prob->comps[i].class0_hp = vp8_bin_prob_from_distribution(
          NMVcount->comps[i].class0_hp);
      branch_ct_class0_hp[i][0] = NMVcount->comps[i].class0_hp[0];
      branch_ct_class0_hp[i][1] = NMVcount->comps[i].class0_hp[1];

      prob->comps[i].hp =
          vp8_bin_prob_from_distribution(NMVcount->comps[i].hp);
      branch_ct_hp[i][0] = NMVcount->comps[i].hp[0];
      branch_ct_hp[i][1] = NMVcount->comps[i].hp[1];
    }
  }
}

void vp8_adapt_nmv_probs(VP8_COMMON *cm, int usehp) {
  int i, j, k;
  nmv_context prob;
  unsigned int branch_ct_joint[MV_JOINTS - 1][2];
  unsigned int branch_ct_sign[2][2];
  unsigned int branch_ct_classes[2][MV_CLASSES - 1][2];
  unsigned int branch_ct_class0[2][CLASS0_SIZE - 1][2];
  unsigned int branch_ct_bits[2][MV_OFFSET_BITS][2];
  unsigned int branch_ct_class0_fp[2][CLASS0_SIZE][4 - 1][2];
  unsigned int branch_ct_fp[2][4 - 1][2];
  unsigned int branch_ct_class0_hp[2][2];
  unsigned int branch_ct_hp[2][2];
#ifdef MV_COUNT_TESTING
  printf("joints count: ");
  for (j = 0; j < MV_JOINTS; ++j) printf("%d ", cm->fc.NMVcount.joints[j]);
  printf("\n"); fflush(stdout);
  printf("signs count:\n");
  for (i = 0; i < 2; ++i)
    printf("%d/%d ", cm->fc.NMVcount.comps[i].sign[0], cm->fc.NMVcount.comps[i].sign[1]);
  printf("\n"); fflush(stdout);
  printf("classes count:\n");
  for (i = 0; i < 2; ++i) {
    for (j = 0; j < MV_CLASSES; ++j)
      printf("%d ", cm->fc.NMVcount.comps[i].classes[j]);
    printf("\n"); fflush(stdout);
  }
  printf("class0 count:\n");
  for (i = 0; i < 2; ++i) {
    for (j = 0; j < CLASS0_SIZE; ++j)
      printf("%d ", cm->fc.NMVcount.comps[i].class0[j]);
    printf("\n"); fflush(stdout);
  }
  printf("bits count:\n");
  for (i = 0; i < 2; ++i) {
    for (j = 0; j < MV_OFFSET_BITS; ++j)
      printf("%d/%d ", cm->fc.NMVcount.comps[i].bits[j][0],
                       cm->fc.NMVcount.comps[i].bits[j][1]);
    printf("\n"); fflush(stdout);
  }
  printf("class0_fp count:\n");
  for (i = 0; i < 2; ++i) {
    for (j = 0; j < CLASS0_SIZE; ++j) {
      printf("{");
      for (k = 0; k < 4; ++k)
        printf("%d ", cm->fc.NMVcount.comps[i].class0_fp[j][k]);
      printf("}, ");
    }
    printf("\n"); fflush(stdout);
  }
  printf("fp count:\n");
  for (i = 0; i < 2; ++i) {
    for (j = 0; j < 4; ++j)
      printf("%d ", cm->fc.NMVcount.comps[i].fp[j]);
    printf("\n"); fflush(stdout);
  }
  if (usehp) {
    printf("class0_hp count:\n");
    for (i = 0; i < 2; ++i)
      printf("%d/%d ", cm->fc.NMVcount.comps[i].class0_hp[0],
                       cm->fc.NMVcount.comps[i].class0_hp[1]);
    printf("\n"); fflush(stdout);
    printf("hp count:\n");
    for (i = 0; i < 2; ++i)
      printf("%d/%d ", cm->fc.NMVcount.comps[i].hp[0],
                       cm->fc.NMVcount.comps[i].hp[1]);
    printf("\n"); fflush(stdout);
  }
#endif
#ifdef SMOOTH_MV_COUNTS
  smooth_counts(&cm->fc.NMVcount.comps[0]);
  smooth_counts(&cm->fc.NMVcount.comps[1]);
#endif
  vp8_counts_to_nmv_context(&cm->fc.NMVcount,
                            &prob,
                            usehp,
                            branch_ct_joint,
                            branch_ct_sign,
                            branch_ct_classes,
                            branch_ct_class0,
                            branch_ct_bits,
                            branch_ct_class0_fp,
                            branch_ct_fp,
                            branch_ct_class0_hp,
                            branch_ct_hp);

  for (j = 0; j < MV_JOINTS - 1; ++j) {
    adapt_prob(&cm->fc.nmvc.joints[j],
               cm->fc.pre_nmvc.joints[j],
               prob.joints[j],
               branch_ct_joint[j]);
  }
  for (i = 0; i < 2; ++i) {
    adapt_prob(&cm->fc.nmvc.comps[i].sign,
               cm->fc.pre_nmvc.comps[i].sign,
               prob.comps[i].sign,
               branch_ct_sign[i]);
    for (j = 0; j < MV_CLASSES - 1; ++j) {
      adapt_prob(&cm->fc.nmvc.comps[i].classes[j],
                 cm->fc.pre_nmvc.comps[i].classes[j],
                 prob.comps[i].classes[j],
                 branch_ct_classes[i][j]);
    }
    for (j = 0; j < CLASS0_SIZE - 1; ++j) {
      adapt_prob(&cm->fc.nmvc.comps[i].class0[j],
                 cm->fc.pre_nmvc.comps[i].class0[j],
                 prob.comps[i].class0[j],
                 branch_ct_class0[i][j]);
    }
    for (j = 0; j < MV_OFFSET_BITS; ++j) {
      adapt_prob(&cm->fc.nmvc.comps[i].bits[j],
                 cm->fc.pre_nmvc.comps[i].bits[j],
                 prob.comps[i].bits[j],
                 branch_ct_bits[i][j]);
    }
  }
  for (i = 0; i < 2; ++i) {
    for (j = 0; j < CLASS0_SIZE; ++j) {
      for (k = 0; k < 3; ++k) {
        adapt_prob(&cm->fc.nmvc.comps[i].class0_fp[j][k],
                   cm->fc.pre_nmvc.comps[i].class0_fp[j][k],
                   prob.comps[i].class0_fp[j][k],
                   branch_ct_class0_fp[i][j][k]);
      }
    }
    for (j = 0; j < 3; ++j) {
      adapt_prob(&cm->fc.nmvc.comps[i].fp[j],
                 cm->fc.pre_nmvc.comps[i].fp[j],
                 prob.comps[i].fp[j],
                 branch_ct_fp[i][j]);
    }
  }
  if (usehp) {
    for (i = 0; i < 2; ++i) {
      adapt_prob(&cm->fc.nmvc.comps[i].class0_hp,
                 cm->fc.pre_nmvc.comps[i].class0_hp,
                 prob.comps[i].class0_hp,
                 branch_ct_class0_hp[i]);
      adapt_prob(&cm->fc.nmvc.comps[i].hp,
                 cm->fc.pre_nmvc.comps[i].hp,
                 prob.comps[i].hp,
                 branch_ct_hp[i]);
    }
  }
}

#else   /* CONFIG_NEWMVENTROPY */

#define MV_COUNT_SAT 16
#define MV_MAX_UPDATE_FACTOR 128

const MV_CONTEXT_HP vp8_mv_update_probs_hp[2] = {
  {{
      237,
      246,
      253, 253, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254,
      254, 254, 254, 254, 254, 250, 250, 252, 254, 254, 254
    }
  },
  {{
      231,
      243,
      245, 253, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254,
      254, 254, 254, 254, 254, 251, 251, 254, 254, 254, 254
    }
  }
};
const MV_CONTEXT_HP vp8_default_mv_context_hp[2] = {
  {{
      /* row */
      162,                                        /* is short */
      128,                                        /* sign */
      220, 204, 180, 192, 192, 119, 192, 192, 180, 140, 192, 192, 224, 224, 224, /* short tree */
      128, 129, 132,  75, 145, 178, 206, 239, 254, 254, 254 /* long bits */
    }
  },
  {{
      /* same for column */
      164,                                        /* is short */
      128,
      220, 204, 180, 192, 192, 119, 192, 192, 180, 140, 192, 192, 224, 224, 224, /* short tree */
      128, 130, 130,  74, 148, 180, 203, 236, 254, 254, 254 /* long bits */
    }
  }
};

const MV_CONTEXT vp8_mv_update_probs[2] = {
  {{
      237,
      246,
      253, 253, 254, 254, 254, 254, 254,
      254, 254, 254, 254, 254, 250, 250, 252, 254, 254
    }
  },
  {{
      231,
      243,
      245, 253, 254, 254, 254, 254, 254,
      254, 254, 254, 254, 254, 251, 251, 254, 254, 254
    }
  }
};
const MV_CONTEXT vp8_default_mv_context[2] = {
  {{
      /* row */
      162,                                        /* is short */
      128,                                        /* sign */
      225, 146, 172, 147, 214,  39, 156,          /* short tree */
      128, 129, 132,  75, 145, 178, 206, 239, 254, 254 /* long bits */
    }
  },
  {{
      /* same for column */
      164,                                        /* is short */
      128,
      204, 170, 119, 235, 140, 230, 228,
      128, 130, 130,  74, 148, 180, 203, 236, 254, 254 /* long bits */
    }
  }
};

const vp8_tree_index vp8_small_mvtree_hp [30] = {
  2,  16,
  4,  10,
  6,   8,
  -0,  -1,
  -2,  -3,
  12,  14,
  -4,  -5,
  -6,  -7,
  18,  24,
  20,  22,
  -8,  -9,
  -10, -11,
  26,  28,
  -12, -13,
  -14, -15
};
struct vp8_token_struct vp8_small_mvencodings_hp [16];

const vp8_tree_index vp8_small_mvtree [14] = {
  2, 8,
  4, 6,
  -0, -1,
  -2, -3,
  10, 12,
  -4, -5,
  -6, -7
};
struct vp8_token_struct vp8_small_mvencodings [8];

__inline static void calc_prob(vp8_prob *p, const unsigned int ct[2], int pbits) {
  const unsigned int tot = ct[0] + ct[1];
  if (tot) {
    const vp8_prob x = ((ct[0] * 255) / tot) & -(1 << (8 - pbits));
    *p = x ? x : 1;
  } else {
    *p = 128;
  }
}

static void compute_component_probs(
  const unsigned int events [MVvals],
  vp8_prob Pnew [MVPcount],
  unsigned int is_short_ct[2],
  unsigned int sign_ct[2],
  unsigned int bit_ct [mvlong_width] [2],
  unsigned int short_ct  [mvnum_short],
  unsigned int short_bct [mvnum_short - 1] [2]
) {
  is_short_ct[0] = is_short_ct[1] = 0;
  sign_ct[0] = sign_ct[1] = 0;
  vpx_memset(bit_ct, 0, sizeof(unsigned int)*mvlong_width * 2);
  vpx_memset(short_ct, 0, sizeof(unsigned int)*mvnum_short);
  vpx_memset(short_bct, 0, sizeof(unsigned int) * (mvnum_short - 1) * 2);

  {
    const int c = events [mv_max];
    is_short_ct [0] += c;    // Short vector
    short_ct [0] += c;       // Magnitude distribution
  }
  {
    int j = 1;
    do {
      const int c1 = events [mv_max + j];  // positive
      const int c2 = events [mv_max - j];  // negative
      const int c  = c1 + c2;
      int a = j;

      sign_ct [0] += c1;
      sign_ct [1] += c2;

      if (a < mvnum_short) {
        is_short_ct [0] += c;     // Short vector
        short_ct [a] += c;       // Magnitude distribution
      } else {
        int k = mvlong_width - 1;
        is_short_ct [1] += c;     // Long vector

        do
          bit_ct [k] [(a >> k) & 1] += c;

        while (--k >= 0);
      }
    } while (++j <= mv_max);
  }
  calc_prob(Pnew + mvpis_short, is_short_ct, 8);

  calc_prob(Pnew + MVPsign, sign_ct, 8);

  {
    vp8_prob p [mvnum_short - 1];    /* actually only need branch ct */
    int j = 0;

    vp8_tree_probs_from_distribution(
      mvnum_short, vp8_small_mvencodings, vp8_small_mvtree,
      p, short_bct, short_ct,
      256, 1
    );

    do
      calc_prob(Pnew + MVPshort + j, short_bct[j], 8);
    while (++j < mvnum_short - 1);
  }

  {
    int j = 0;
    do
      calc_prob(Pnew + MVPbits + j, bit_ct[j], 8);
    while (++j < mvlong_width);
  }
}

static void compute_component_probs_hp(
  const unsigned int events [MVvals_hp],
  vp8_prob Pnew [MVPcount_hp],
  unsigned int is_short_ct[2],
  unsigned int sign_ct[2],
  unsigned int bit_ct [mvlong_width_hp] [2],
  unsigned int short_ct  [mvnum_short_hp],
  unsigned int short_bct [mvnum_short_hp - 1] [2]
) {
  is_short_ct[0] = is_short_ct[1] = 0;
  sign_ct[0] = sign_ct[1] = 0;
  vpx_memset(bit_ct, 0, sizeof(unsigned int)*mvlong_width_hp * 2);
  vpx_memset(short_ct, 0, sizeof(unsigned int)*mvnum_short_hp);
  vpx_memset(short_bct, 0, sizeof(unsigned int) * (mvnum_short_hp - 1) * 2);

  {
    const int c = events [mv_max_hp];
    is_short_ct [0] += c;    // Short vector
    short_ct [0] += c;       // Magnitude distribution
  }
  {
    int j = 1;
    do {
      const int c1 = events [mv_max_hp + j];  // positive
      const int c2 = events [mv_max_hp - j];  // negative
      const int c  = c1 + c2;
      int a = j;

      sign_ct [0] += c1;
      sign_ct [1] += c2;

      if (a < mvnum_short_hp) {
        is_short_ct [0] += c;     // Short vector
        short_ct [a] += c;       // Magnitude distribution
      } else {
        int k = mvlong_width_hp - 1;
        is_short_ct [1] += c;     // Long vector

        do
          bit_ct [k] [(a >> k) & 1] += c;

        while (--k >= 0);
      }
    } while (++j <= mv_max_hp);
  }
  calc_prob(Pnew + mvpis_short_hp, is_short_ct, 8);

  calc_prob(Pnew + MVPsign_hp, sign_ct, 8);

  {
    vp8_prob p [mvnum_short_hp - 1];    /* actually only need branch ct */
    int j = 0;

    vp8_tree_probs_from_distribution(
      mvnum_short_hp, vp8_small_mvencodings_hp, vp8_small_mvtree_hp,
      p, short_bct, short_ct,
      256, 1
    );

    do
      calc_prob(Pnew + MVPshort_hp + j, short_bct[j], 8);
    while (++j < mvnum_short_hp - 1);
  }

  {
    int j = 0;
    do
      calc_prob(Pnew + MVPbits_hp + j, bit_ct[j], 8);
    while (++j < mvlong_width_hp);
  }
}

void vp8_adapt_mv_probs(VP8_COMMON *cm) {
  int i, t, count, factor;
#ifdef MV_COUNT_TESTING
  printf("static const unsigned int\nMVcount[2][MVvals]={\n");
  for (i = 0; i < 2; ++i) {
    printf("  { ");
    for (t = 0; t < MVvals; t++) {
      printf("%d, ", cm->fc.MVcount[i][t]);
      if (t % 16 == 15 && t != MVvals - 1) printf("\n    ");
    }
    printf("},\n");
  }
  printf("};\n");
  printf("static const unsigned int\nMVcount_hp[2][MVvals_hp]={\n");
  for (i = 0; i < 2; ++i) {
    printf("  { ");
    for (t = 0; t < MVvals_hp; t++) {
      printf("%d, ", cm->fc.MVcount_hp[i][t]);
      if (t % 16 == 15 && t != MVvals_hp - 1) printf("\n    ");
    }
    printf("},\n");
  }
  printf("};\n");
#endif  /* MV_COUNT_TESTING */

  for (i = 0; i < 2; ++i) {
    int prob;
    unsigned int is_short_ct[2];
    unsigned int sign_ct[2];
    unsigned int bit_ct [mvlong_width] [2];
    unsigned int short_ct  [mvnum_short];
    unsigned int short_bct [mvnum_short - 1] [2];
    vp8_prob Pnew [MVPcount];
    compute_component_probs(cm->fc.MVcount[i], Pnew,
                            is_short_ct, sign_ct,
                            bit_ct, short_ct, short_bct);
    count = is_short_ct[0] + is_short_ct[1];
    count = count > MV_COUNT_SAT ? MV_COUNT_SAT : count;
    factor = (MV_MAX_UPDATE_FACTOR * count / MV_COUNT_SAT);
    prob = ((int)cm->fc.pre_mvc[i].prob[mvpis_short] * (256 - factor) +
            (int)Pnew[mvpis_short] * factor + 128) >> 8;
    if (prob <= 0) cm->fc.mvc[i].prob[mvpis_short] = 1;
    else if (prob > 255) cm->fc.mvc[i].prob[mvpis_short] = 255;
    else cm->fc.mvc[i].prob[mvpis_short] = prob;

    count = sign_ct[0] + sign_ct[1];
    count = count > MV_COUNT_SAT ? MV_COUNT_SAT : count;
    factor = (MV_MAX_UPDATE_FACTOR * count / MV_COUNT_SAT);
    prob = ((int)cm->fc.pre_mvc[i].prob[MVPsign] * (256 - factor) +
            (int)Pnew[MVPsign] * factor + 128) >> 8;
    if (prob <= 0) cm->fc.mvc[i].prob[MVPsign] = 1;
    else if (prob > 255) cm->fc.mvc[i].prob[MVPsign] = 255;
    else cm->fc.mvc[i].prob[MVPsign] = prob;

    for (t = 0; t < mvnum_short - 1; ++t) {
      count = short_bct[t][0] + short_bct[t][1];
      count = count > MV_COUNT_SAT ? MV_COUNT_SAT : count;
      factor = (MV_MAX_UPDATE_FACTOR * count / MV_COUNT_SAT);
      prob = ((int)cm->fc.pre_mvc[i].prob[MVPshort + t] * (256 - factor) +
              (int)Pnew[MVPshort + t] * factor + 128) >> 8;
      if (prob <= 0) cm->fc.mvc[i].prob[MVPshort + t] = 1;
      else if (prob > 255) cm->fc.mvc[i].prob[MVPshort + t] = 255;
      else cm->fc.mvc[i].prob[MVPshort + t] = prob;
    }
    for (t = 0; t < mvlong_width; ++t) {
      count = bit_ct[t][0] + bit_ct[t][1];
      count = count > MV_COUNT_SAT ? MV_COUNT_SAT : count;
      factor = (MV_MAX_UPDATE_FACTOR * count / MV_COUNT_SAT);
      prob = ((int)cm->fc.pre_mvc[i].prob[MVPbits + t] * (256 - factor) +
              (int)Pnew[MVPbits + t] * factor + 128) >> 8;
      if (prob <= 0) cm->fc.mvc[i].prob[MVPbits + t] = 1;
      else if (prob > 255) cm->fc.mvc[i].prob[MVPbits + t] = 255;
      else cm->fc.mvc[i].prob[MVPbits + t] = prob;
    }
  }
  for (i = 0; i < 2; ++i) {
    int prob;
    unsigned int is_short_ct[2];
    unsigned int sign_ct[2];
    unsigned int bit_ct [mvlong_width_hp] [2];
    unsigned int short_ct  [mvnum_short_hp];
    unsigned int short_bct [mvnum_short_hp - 1] [2];
    vp8_prob Pnew [MVPcount_hp];
    compute_component_probs_hp(cm->fc.MVcount_hp[i], Pnew,
                               is_short_ct, sign_ct,
                               bit_ct, short_ct, short_bct);
    count = is_short_ct[0] + is_short_ct[1];
    count = count > MV_COUNT_SAT ? MV_COUNT_SAT : count;
    factor = (MV_MAX_UPDATE_FACTOR * count / MV_COUNT_SAT);
    prob = ((int)cm->fc.pre_mvc_hp[i].prob[mvpis_short_hp] * (256 - factor) +
            (int)Pnew[mvpis_short_hp] * factor + 128) >> 8;
    if (prob <= 0) cm->fc.mvc_hp[i].prob[mvpis_short_hp] = 1;
    else if (prob > 255) cm->fc.mvc_hp[i].prob[mvpis_short_hp] = 255;
    else cm->fc.mvc_hp[i].prob[mvpis_short_hp] = prob;

    count = sign_ct[0] + sign_ct[1];
    count = count > MV_COUNT_SAT ? MV_COUNT_SAT : count;
    factor = (MV_MAX_UPDATE_FACTOR * count / MV_COUNT_SAT);
    prob = ((int)cm->fc.pre_mvc_hp[i].prob[MVPsign_hp] * (256 - factor) +
            (int)Pnew[MVPsign_hp] * factor + 128) >> 8;
    if (prob <= 0) cm->fc.mvc_hp[i].prob[MVPsign_hp] = 1;
    else if (prob > 255) cm->fc.mvc_hp[i].prob[MVPsign_hp] = 255;
    else cm->fc.mvc_hp[i].prob[MVPsign_hp] = prob;

    for (t = 0; t < mvnum_short_hp - 1; ++t) {
      count = short_bct[t][0] + short_bct[t][1];
      count = count > MV_COUNT_SAT ? MV_COUNT_SAT : count;
      factor = (MV_MAX_UPDATE_FACTOR * count / MV_COUNT_SAT);
      prob = ((int)cm->fc.pre_mvc_hp[i].prob[MVPshort_hp + t] * (256 - factor) +
              (int)Pnew[MVPshort_hp + t] * factor + 128) >> 8;
      if (prob <= 0) cm->fc.mvc_hp[i].prob[MVPshort_hp + t] = 1;
      else if (prob > 255) cm->fc.mvc_hp[i].prob[MVPshort_hp + t] = 255;
      else cm->fc.mvc_hp[i].prob[MVPshort_hp + t] = prob;
    }
    for (t = 0; t < mvlong_width_hp; ++t) {
      count = bit_ct[t][0] + bit_ct[t][1];
      count = count > MV_COUNT_SAT ? MV_COUNT_SAT : count;
      factor = (MV_MAX_UPDATE_FACTOR * count / MV_COUNT_SAT);
      prob = ((int)cm->fc.pre_mvc_hp[i].prob[MVPbits_hp + t] * (256 - factor) +
              (int)Pnew[MVPbits_hp + t] * factor + 128) >> 8;
      if (prob <= 0) cm->fc.mvc_hp[i].prob[MVPbits_hp + t] = 1;
      else if (prob > 255) cm->fc.mvc_hp[i].prob[MVPbits_hp + t] = 255;
      else cm->fc.mvc_hp[i].prob[MVPbits_hp + t] = prob;
    }
  }
}

#endif  /* CONFIG_NEWMVENTROPY */

void vp8_entropy_mv_init() {
#if CONFIG_NEWMVENTROPY
  vp8_tokens_from_tree(vp8_mv_joint_encodings, vp8_mv_joint_tree);
  vp8_tokens_from_tree(vp8_mv_class_encodings, vp8_mv_class_tree);
  vp8_tokens_from_tree(vp8_mv_class0_encodings, vp8_mv_class0_tree);
  vp8_tokens_from_tree(vp8_mv_fp_encodings, vp8_mv_fp_tree);
#else
  vp8_tokens_from_tree(vp8_small_mvencodings, vp8_small_mvtree);
  vp8_tokens_from_tree(vp8_small_mvencodings_hp, vp8_small_mvtree_hp);
#endif
}

void vp8_init_mv_probs(VP8_COMMON *cm) {
#if CONFIG_NEWMVENTROPY
  vpx_memcpy(&cm->fc.nmvc, &vp8_default_nmv_context, sizeof(nmv_context));
#else
  vpx_memcpy(cm->fc.mvc,
             vp8_default_mv_context, sizeof(vp8_default_mv_context));
  vpx_memcpy(cm->fc.mvc_hp,
             vp8_default_mv_context_hp, sizeof(vp8_default_mv_context_hp));
#endif
}
