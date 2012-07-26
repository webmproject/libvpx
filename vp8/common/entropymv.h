/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */


#ifndef __INC_ENTROPYMV_H
#define __INC_ENTROPYMV_H

#include "treecoder.h"
#include "vpx_config.h"
#include "blockd.h"

struct VP8Common;

void vp8_entropy_mv_init();
void vp8_init_mv_probs(struct VP8Common *cm);
void vp8_adapt_mv_probs(struct VP8Common *cm);
#if CONFIG_NEWMVENTROPY
void vp8_adapt_nmv_probs(struct VP8Common *cm, int usehp);
#endif

#if CONFIG_NEWMVENTROPY

#define VP8_NMV_UPDATE_PROB  255
//#define MV_GROUP_UPDATE

#define LOW_PRECISION_MV_UPDATE  /* Use 7 bit forward update */

/* Symbols for coding which components are zero jointly */
#define MV_JOINTS     4
typedef enum {
  MV_JOINT_ZERO = 0,             /* Zero vector */
  MV_JOINT_HNZVZ = 1,            /* Vert zero, hor nonzero */
  MV_JOINT_HZVNZ = 2,            /* Hor zero, vert nonzero */
  MV_JOINT_HNZVNZ = 3,           /* Both components nonzero */
} MV_JOINT_TYPE;

extern const vp8_tree_index vp8_mv_joint_tree[2 * MV_JOINTS - 2];
extern struct vp8_token_struct vp8_mv_joint_encodings [MV_JOINTS];

/* Symbols for coding magnitude class of nonzero components */
#define MV_CLASSES     8
typedef enum {
  MV_CLASS_0 = 0,      /* (0, 2]     integer pel */
  MV_CLASS_1 = 1,      /* (2, 4]     integer pel */
  MV_CLASS_2 = 2,      /* (4, 8]     integer pel */
  MV_CLASS_3 = 3,      /* (8, 16]    integer pel */
  MV_CLASS_4 = 4,      /* (16, 32]   integer pel */
  MV_CLASS_5 = 5,      /* (32, 64]   integer pel */
  MV_CLASS_6 = 6,      /* (64, 128]  integer pel */
  MV_CLASS_7 = 7,      /* (128, 256] integer pel */
} MV_CLASS_TYPE;

extern const vp8_tree_index vp8_mv_class_tree[2 * MV_CLASSES - 2];
extern struct vp8_token_struct vp8_mv_class_encodings [MV_CLASSES];

#define CLASS0_BITS    1  /* bits at integer precision for class 0 */
#define CLASS0_SIZE    (1 << CLASS0_BITS)
#define MV_OFFSET_BITS (MV_CLASSES + CLASS0_BITS - 2)

#define MV_MAX_BITS    (MV_CLASSES + CLASS0_BITS + 2)
#define MV_MAX         ((1 << MV_MAX_BITS) - 1)
#define MV_VALS        ((MV_MAX << 1) + 1)

extern const vp8_tree_index vp8_mv_class0_tree[2 * CLASS0_SIZE - 2];
extern struct vp8_token_struct vp8_mv_class0_encodings[CLASS0_SIZE];

extern const vp8_tree_index vp8_mv_fp_tree[2 * 4 - 2];
extern struct vp8_token_struct vp8_mv_fp_encodings[4];

typedef struct {
  vp8_prob sign;
  vp8_prob classes[MV_CLASSES - 1];
  vp8_prob class0[CLASS0_SIZE - 1];
  vp8_prob bits[MV_OFFSET_BITS];
  vp8_prob class0_fp[CLASS0_SIZE][4 - 1];
  vp8_prob fp[4 - 1];
  vp8_prob class0_hp;
  vp8_prob hp;
} nmv_component;

typedef struct {
  vp8_prob joints[MV_JOINTS - 1];
  nmv_component comps[2];
} nmv_context;

MV_JOINT_TYPE vp8_get_mv_joint(MV mv);
MV_CLASS_TYPE vp8_get_mv_class(int z, int *offset);
int vp8_get_mv_mag(MV_CLASS_TYPE c, int offset);


typedef struct {
  unsigned int mvcount[MV_VALS];
  unsigned int sign[2];
  unsigned int classes[MV_CLASSES];
  unsigned int class0[CLASS0_SIZE];
  unsigned int bits[MV_OFFSET_BITS][2];
  unsigned int class0_fp[CLASS0_SIZE][4];
  unsigned int fp[4];
  unsigned int class0_hp[2];
  unsigned int hp[2];
} nmv_component_counts;

typedef struct {
  unsigned int joints[MV_JOINTS];
  nmv_component_counts comps[2];
} nmv_context_counts;

void vp8_increment_nmv(const MV *mv, const MV *ref, nmv_context_counts *mvctx,
                       int usehp);
extern const nmv_context vp8_default_nmv_context;
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
    unsigned int (*branch_ct_hp)[2]);

#else  /* CONFIG_NEWMVENTROPY */

enum {
  mv_max  = 1023,              /* max absolute value of a MV component */
  MVvals = (2 * mv_max) + 1,   /* # possible values "" */
  mvlong_width = 10,       /* Large MVs have 9 bit magnitudes */
  mvnum_short = 8,         /* magnitudes 0 through 7 */
  mvnum_short_bits = 3,         /* number of bits for short mvs */

  mvfp_max  = 255,              /* max absolute value of a full pixel MV component */
  MVfpvals = (2 * mvfp_max) + 1, /* # possible full pixel MV values */

  /* probability offsets for coding each MV component */

  mvpis_short = 0,         /* short (<= 7) vs long (>= 8) */
  MVPsign,                /* sign for non-zero */
  MVPshort,               /* 8 short values = 7-position tree */

  MVPbits = MVPshort + mvnum_short - 1, /* mvlong_width long value bits */
  MVPcount = MVPbits + mvlong_width    /* (with independent probabilities) */
};

typedef struct mv_context {
  vp8_prob prob[MVPcount];  /* often come in row, col pairs */
} MV_CONTEXT;

extern const MV_CONTEXT vp8_mv_update_probs[2], vp8_default_mv_context[2];

enum {
  mv_max_hp  = 2047,              /* max absolute value of a MV component */
  MVvals_hp = (2 * mv_max_hp) + 1,   /* # possible values "" */
  mvlong_width_hp = 11,       /* Large MVs have 9 bit magnitudes */
  mvnum_short_hp = 16,         /* magnitudes 0 through 15 */
  mvnum_short_bits_hp = 4,         /* number of bits for short mvs */

  mvfp_max_hp  = 255,              /* max absolute value of a full pixel MV component */
  MVfpvals_hp = (2 * mvfp_max_hp) + 1, /* # possible full pixel MV values */

  /* probability offsets for coding each MV component */

  mvpis_short_hp = 0,         /* short (<= 7) vs long (>= 8) */
  MVPsign_hp,                /* sign for non-zero */
  MVPshort_hp,               /* 8 short values = 7-position tree */

  MVPbits_hp = MVPshort_hp + mvnum_short_hp - 1, /* mvlong_width long value bits */
  MVPcount_hp = MVPbits_hp + mvlong_width_hp    /* (with independent probabilities) */
};

typedef struct mv_context_hp {
  vp8_prob prob[MVPcount_hp];  /* often come in row, col pairs */
} MV_CONTEXT_HP;

extern const MV_CONTEXT_HP vp8_mv_update_probs_hp[2], vp8_default_mv_context_hp[2];

extern const vp8_tree_index vp8_small_mvtree[];
extern struct vp8_token_struct vp8_small_mvencodings [8];
extern const vp8_tree_index vp8_small_mvtree_hp[];
extern struct vp8_token_struct vp8_small_mvencodings_hp [16];

#endif  /* CONFIG_NEWMVENTROPY */

#endif
