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
#include "vp9/common/vp9_entropymv.h"

#define MV_COUNT_SAT 20
#define MV_MAX_UPDATE_FACTOR 128

/* Integer pel reference mv threshold for use of high-precision 1/8 mv */
#define COMPANDED_MVREF_THRESH 8

const vp9_tree_index vp9_mv_joint_tree[2 * MV_JOINTS - 2] = {
  -MV_JOINT_ZERO, 2,
  -MV_JOINT_HNZVZ, 4,
  -MV_JOINT_HZVNZ, -MV_JOINT_HNZVNZ
};
struct vp9_token vp9_mv_joint_encodings[MV_JOINTS];

const vp9_tree_index vp9_mv_class_tree[2 * MV_CLASSES - 2] = {
  -MV_CLASS_0, 2,
  -MV_CLASS_1, 4,
  6, 8,
  -MV_CLASS_2, -MV_CLASS_3,
  10, 12,
  -MV_CLASS_4, -MV_CLASS_5,
  -MV_CLASS_6, 14,
  16, 18,
  -MV_CLASS_7, -MV_CLASS_8,
  -MV_CLASS_9, -MV_CLASS_10,
};
struct vp9_token vp9_mv_class_encodings[MV_CLASSES];

const vp9_tree_index vp9_mv_class0_tree [2 * CLASS0_SIZE - 2] = {
  -0, -1,
};
struct vp9_token vp9_mv_class0_encodings[CLASS0_SIZE];

const vp9_tree_index vp9_mv_fp_tree [2 * 4 - 2] = {
  -0, 2,
  -1, 4,
  -2, -3
};
struct vp9_token vp9_mv_fp_encodings[4];

static const nmv_context default_nmv_context = {
  {32, 64, 96},
  {
    { /* vert component */
      128,                                                  /* sign */
      {224, 144, 192, 168, 192, 176, 192, 198, 198, 245},   /* class */
      {216},                                                /* class0 */
      {136, 140, 148, 160, 176, 192, 224, 234, 234, 240},   /* bits */
      {{128, 128, 64}, {96, 112, 64}},                      /* class0_fp */
      {64, 96, 64},                                         /* fp */
      160,                                                  /* class0_hp bit */
      128,                                                  /* hp */
    },
    { /* hor component */
      128,                                                  /* sign */
      {216, 128, 176, 160, 176, 176, 192, 198, 198, 208},   /* class */
      {208},                                                /* class0 */
      {136, 140, 148, 160, 176, 192, 224, 234, 234, 240},   /* bits */
      {{128, 128, 64}, {96, 112, 64}},                      /* class0_fp */
      {64, 96, 64},                                         /* fp */
      160,                                                  /* class0_hp bit */
      128,                                                  /* hp */
    }
  },
};

#define mv_class_base(c) ((c) ? (CLASS0_SIZE << (c + 2)) : 0)

static const uint8_t log_in_base_2[] = {
  0, 0, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4,
  4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
  5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6,
  6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
  6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
  6, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
  7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
  7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
  7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
  7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
  7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 8, 8, 8, 8, 8, 8, 8, 8,
  8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
  8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
  8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
  8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
  8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
  8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
  8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
  8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
  8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
  8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
  8, 8, 8, 8, 8, 8, 8, 8, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9,
  9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9,
  9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9,
  9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9,
  9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9,
  9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9,
  9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9,
  9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9,
  9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9,
  9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9,
  9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9,
  9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9,
  9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9,
  9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9,
  9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9,
  9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9,
  9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9,
  9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9,
  9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9,
  9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9,
  9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9,
  9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 10
};

MV_CLASS_TYPE vp9_get_mv_class(int z, int *offset) {
  MV_CLASS_TYPE c = MV_CLASS_0;
  if (z >= CLASS0_SIZE * 4096)
    c = MV_CLASS_10;
  else
    c = log_in_base_2[z >> 3];

  if (offset)
    *offset = z - mv_class_base(c);
  return c;
}

int vp9_use_mv_hp(const MV *ref) {
  return (abs(ref->row) >> 3) < COMPANDED_MVREF_THRESH &&
         (abs(ref->col) >> 3) < COMPANDED_MVREF_THRESH;
}

int vp9_get_mv_mag(MV_CLASS_TYPE c, int offset) {
  return mv_class_base(c) + offset;
}

static void inc_mv_component(int v, nmv_component_counts *comp_counts,
                             int incr, int usehp) {
  int s, z, c, o, d, e, f;
  if (!incr)
    return;
  assert (v != 0);            /* should not be zero */
  s = v < 0;
  comp_counts->sign[s] += incr;
  z = (s ? -v : v) - 1;       /* magnitude - 1 */

  c = vp9_get_mv_class(z, &o);
  comp_counts->classes[c] += incr;

  d = (o >> 3);               /* int mv data */
  f = (o >> 1) & 3;           /* fractional pel mv data */
  e = (o & 1);                /* high precision mv data */

  if (c == MV_CLASS_0) {
    comp_counts->class0[d] += incr;
    comp_counts->class0_fp[d][f] += incr;
    comp_counts->class0_hp[e] += usehp * incr;
  } else {
    int i;
    int b = c + CLASS0_BITS - 1;  // number of bits
    for (i = 0; i < b; ++i)
      comp_counts->bits[i][((d >> i) & 1)] += incr;
    comp_counts->fp[f] += incr;
    comp_counts->hp[e] += usehp * incr;
  }
}

static void counts_to_context(nmv_component_counts *mvcomp, int usehp) {
  int v;
  vpx_memset(mvcomp->sign, 0, sizeof(nmv_component_counts) - sizeof(mvcomp->mvcount));
  for (v = 1; v <= MV_MAX; v++) {
    inc_mv_component(-v, mvcomp, mvcomp->mvcount[MV_MAX - v], usehp);
    inc_mv_component( v, mvcomp, mvcomp->mvcount[MV_MAX + v], usehp);
  }
}

void vp9_inc_mv(const MV *mv,  nmv_context_counts *counts) {
  const MV_JOINT_TYPE j = vp9_get_mv_joint(mv);
  ++counts->joints[j];

  if (mv_joint_vertical(j))
    ++counts->comps[0].mvcount[MV_MAX + mv->row];

  if (mv_joint_horizontal(j))
    ++counts->comps[1].mvcount[MV_MAX + mv->col];
}

static vp9_prob adapt_prob(vp9_prob prep, const unsigned int ct[2]) {
  return merge_probs2(prep, ct, MV_COUNT_SAT, MV_MAX_UPDATE_FACTOR);
}

void vp9_counts_process(nmv_context_counts *nmv_count, int usehp) {
  counts_to_context(&nmv_count->comps[0], usehp);
  counts_to_context(&nmv_count->comps[1], usehp);
}

static unsigned int adapt_probs(unsigned int i,
                                vp9_tree tree,
                                vp9_prob this_probs[],
                                const vp9_prob last_probs[],
                                const unsigned int num_events[]) {


  const unsigned int left = tree[i] <= 0
          ? num_events[-tree[i]]
          : adapt_probs(tree[i], tree, this_probs, last_probs, num_events);

  const unsigned int right = tree[i + 1] <= 0
          ? num_events[-tree[i + 1]]
          : adapt_probs(tree[i + 1], tree, this_probs, last_probs, num_events);
  const unsigned int ct[2] = { left, right };
  this_probs[i >> 1] = adapt_prob(last_probs[i >> 1], ct);
  return left + right;
}


void vp9_adapt_mv_probs(VP9_COMMON *cm, int allow_hp) {
  int i, j;

  FRAME_CONTEXT *pre_fc = &cm->frame_contexts[cm->frame_context_idx];

  nmv_context *ctx = &cm->fc.nmvc;
  nmv_context *pre_ctx = &pre_fc->nmvc;
  nmv_context_counts *cts = &cm->counts.mv;

  vp9_counts_process(cts, allow_hp);

  adapt_probs(0, vp9_mv_joint_tree, ctx->joints, pre_ctx->joints, cts->joints);

  for (i = 0; i < 2; ++i) {
    ctx->comps[i].sign = adapt_prob(pre_ctx->comps[i].sign, cts->comps[i].sign);
    adapt_probs(0, vp9_mv_class_tree, ctx->comps[i].classes,
                pre_ctx->comps[i].classes, cts->comps[i].classes);
    adapt_probs(0, vp9_mv_class0_tree, ctx->comps[i].class0,
                pre_ctx->comps[i].class0, cts->comps[i].class0);

    for (j = 0; j < MV_OFFSET_BITS; ++j)
        ctx->comps[i].bits[j] = adapt_prob(pre_ctx->comps[i].bits[j],
                                           cts->comps[i].bits[j]);

    for (j = 0; j < CLASS0_SIZE; ++j)
      adapt_probs(0, vp9_mv_fp_tree, ctx->comps[i].class0_fp[j],
                  pre_ctx->comps[i].class0_fp[j], cts->comps[i].class0_fp[j]);

    adapt_probs(0, vp9_mv_fp_tree, ctx->comps[i].fp, pre_ctx->comps[i].fp,
                cts->comps[i].fp);

    if (allow_hp) {
      ctx->comps[i].class0_hp = adapt_prob(pre_ctx->comps[i].class0_hp,
                                           cts->comps[i].class0_hp);
      ctx->comps[i].hp = adapt_prob(pre_ctx->comps[i].hp, cts->comps[i].hp);
    }
  }
}

void vp9_entropy_mv_init() {
  vp9_tokens_from_tree(vp9_mv_joint_encodings, vp9_mv_joint_tree);
  vp9_tokens_from_tree(vp9_mv_class_encodings, vp9_mv_class_tree);
  vp9_tokens_from_tree(vp9_mv_class0_encodings, vp9_mv_class0_tree);
  vp9_tokens_from_tree(vp9_mv_fp_encodings, vp9_mv_fp_tree);
}

void vp9_init_mv_probs(VP9_COMMON *cm) {
  cm->fc.nmvc = default_nmv_context;
}
