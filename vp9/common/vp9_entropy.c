/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */


#include <stdio.h>

#include "vp9/common/vp9_entropy.h"
#include "string.h"
#include "vp9/common/vp9_blockd.h"
#include "vp9/common/vp9_onyxc_int.h"
#include "vp9/common/vp9_entropymode.h"
#include "vpx_mem/vpx_mem.h"

#define uchar unsigned char     /* typedefs can clash */
#define uint  unsigned int

typedef const uchar cuchar;
typedef const uint cuint;

typedef vp9_prob Prob;

#include "vp9/common/vp9_coefupdateprobs.h"

const int vp9_i8x8_block[4] = {0, 2, 8, 10};

DECLARE_ALIGNED(16, const unsigned char, vp9_norm[256]) = {
  0, 7, 6, 6, 5, 5, 5, 5, 4, 4, 4, 4, 4, 4, 4, 4,
  3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
  2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
  2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
  1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
  1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
  1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
  1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
};

DECLARE_ALIGNED(16, const int, vp9_coef_bands[16]) = {
  0, 1, 2, 3, 6, 4, 5, 6, 6, 6, 6, 6, 6, 7, 7, 7
};

DECLARE_ALIGNED(16, cuchar, vp9_prev_token_class[MAX_ENTROPY_TOKENS]) = {
  0, 1, 2, 2, 3, 3, 3, 3, 3, 3, 3, 0
};

DECLARE_ALIGNED(16, const int, vp9_default_zig_zag1d[16]) = {
  0,  1,  4,  8,
  5,  2,  3,  6,
  9, 12, 13, 10,
  7, 11, 14, 15,
};

DECLARE_ALIGNED(16, const int, vp9_col_scan[16]) = {
  0, 4,  8, 12,
  1, 5,  9, 13,
  2, 6, 10, 14,
  3, 7, 11, 15
};
DECLARE_ALIGNED(16, const int, vp9_row_scan[16]) = {
  0,   1,  2,  3,
  4,   5,  6,  7,
  8,   9, 10, 11,
  12, 13, 14, 15
};


DECLARE_ALIGNED(64, const int, vp9_coef_bands_8x8[64]) = { 0, 1, 2, 3, 5, 4, 4, 5,
                                                           5, 3, 6, 3, 5, 4, 6, 6,
                                                           6, 5, 5, 6, 6, 6, 6, 6,
                                                           6, 6, 6, 6, 6, 6, 6, 6,
                                                           6, 6, 6, 6, 7, 7, 7, 7,
                                                           7, 7, 7, 7, 7, 7, 7, 7,
                                                           7, 7, 7, 7, 7, 7, 7, 7,
                                                           7, 7, 7, 7, 7, 7, 7, 7
                                                         };
DECLARE_ALIGNED(64, const int, vp9_default_zig_zag1d_8x8[64]) = {
  0,  1,  8, 16,  9,  2,  3, 10, 17, 24, 32, 25, 18, 11,  4,  5,
  12, 19, 26, 33, 40, 48, 41, 34, 27, 20, 13,  6,  7, 14, 21, 28,
  35, 42, 49, 56, 57, 50, 43, 36, 29, 22, 15, 23, 30, 37, 44, 51,
  58, 59, 52, 45, 38, 31, 39, 46, 53, 60, 61, 54, 47, 55, 62, 63,
};

// Table can be optimized.
DECLARE_ALIGNED(16, const int, vp9_coef_bands_16x16[256]) = {
    0, 1, 2, 3, 5, 4, 4, 5, 5, 3, 6, 3, 5, 4, 6, 6,
    6, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
    6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
    7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
    7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
    7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
    7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
    7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
    7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
    7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
    7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
    7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
    7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
    7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
    7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
    7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
};
DECLARE_ALIGNED(16, const int, vp9_default_zig_zag1d_16x16[256]) = {
      0,   1,  16,  32,  17,   2,   3,  18,  33,  48,  64,  49,  34,  19,   4,   5,
     20,  35,  50,  65,  80,  96,  81,  66,  51,  36,  21,   6,   7,  22,  37,  52,
     67,  82,  97, 112, 128, 113,  98,  83,  68,  53,  38,  23,   8,   9,  24,  39,
     54,  69,  84,  99, 114, 129, 144, 160, 145, 130, 115, 100,  85,  70,  55,  40,
     25,  10,  11,  26,  41,  56,  71,  86, 101, 116, 131, 146, 161, 176, 192, 177,
    162, 147, 132, 117, 102,  87,  72,  57,  42,  27,  12,  13,  28,  43,  58,  73,
     88, 103, 118, 133, 148, 163, 178, 193, 208, 224, 209, 194, 179, 164, 149, 134,
    119, 104,  89,  74,  59,  44,  29,  14,  15,  30,  45,  60,  75,  90, 105, 120,
    135, 150, 165, 180, 195, 210, 225, 240, 241, 226, 211, 196, 181, 166, 151, 136,
    121, 106,  91,  76,  61,  46,  31,  47,  62,  77,  92, 107, 122, 137, 152, 167,
    182, 197, 212, 227, 242, 243, 228, 213, 198, 183, 168, 153, 138, 123, 108,  93,
     78,  63,  79,  94, 109, 124, 139, 154, 169, 184, 199, 214, 229, 244, 245, 230,
    215, 200, 185, 170, 155, 140, 125, 110,  95, 111, 126, 141, 156, 171, 186, 201,
    216, 231, 246, 247, 232, 217, 202, 187, 172, 157, 142, 127, 143, 158, 173, 188,
    203, 218, 233, 248, 249, 234, 219, 204, 189, 174, 159, 175, 190, 205, 220, 235,
    250, 251, 236, 221, 206, 191, 207, 222, 237, 252, 253, 238, 223, 239, 254, 255,
};


/* Array indices are identical to previously-existing CONTEXT_NODE indices */

const vp9_tree_index vp9_coef_tree[ 22] =     /* corresponding _CONTEXT_NODEs */
{
  -DCT_EOB_TOKEN, 2,                             /* 0 = EOB */
  -ZERO_TOKEN, 4,                               /* 1 = ZERO */
  -ONE_TOKEN, 6,                               /* 2 = ONE */
  8, 12,                                      /* 3 = LOW_VAL */
  -TWO_TOKEN, 10,                            /* 4 = TWO */
  -THREE_TOKEN, -FOUR_TOKEN,                /* 5 = THREE */
  14, 16,                                    /* 6 = HIGH_LOW */
  -DCT_VAL_CATEGORY1, -DCT_VAL_CATEGORY2,   /* 7 = CAT_ONE */
  18, 20,                                   /* 8 = CAT_THREEFOUR */
  -DCT_VAL_CATEGORY3, -DCT_VAL_CATEGORY4,  /* 9 = CAT_THREE */
  -DCT_VAL_CATEGORY5, -DCT_VAL_CATEGORY6   /* 10 = CAT_FIVE */
};

struct vp9_token_struct vp9_coef_encodings[MAX_ENTROPY_TOKENS];

/* Trees for extra bits.  Probabilities are constant and
   do not depend on previously encoded bits */

static const Prob Pcat1[] = { 159};
static const Prob Pcat2[] = { 165, 145};
static const Prob Pcat3[] = { 173, 148, 140};
static const Prob Pcat4[] = { 176, 155, 140, 135};
static const Prob Pcat5[] = { 180, 157, 141, 134, 130};
static const Prob Pcat6[] =
{ 254, 254, 252, 249, 243, 230, 196, 177, 153, 140, 133, 130, 129};

static vp9_tree_index cat1[2], cat2[4], cat3[6], cat4[8], cat5[10], cat6[26];

static void init_bit_tree(vp9_tree_index *p, int n) {
  int i = 0;

  while (++i < n) {
    p[0] = p[1] = i << 1;
    p += 2;
  }

  p[0] = p[1] = 0;
}

static void init_bit_trees() {
  init_bit_tree(cat1, 1);
  init_bit_tree(cat2, 2);
  init_bit_tree(cat3, 3);
  init_bit_tree(cat4, 4);
  init_bit_tree(cat5, 5);
  init_bit_tree(cat6, 13);
}

vp9_extra_bit_struct vp9_extra_bits[12] = {
  { 0, 0, 0, 0},
  { 0, 0, 0, 1},
  { 0, 0, 0, 2},
  { 0, 0, 0, 3},
  { 0, 0, 0, 4},
  { cat1, Pcat1, 1, 5},
  { cat2, Pcat2, 2, 7},
  { cat3, Pcat3, 3, 11},
  { cat4, Pcat4, 4, 19},
  { cat5, Pcat5, 5, 35},
  { cat6, Pcat6, 13, 67},
  { 0, 0, 0, 0}
};

#include "vp9/common/vp9_default_coef_probs.h"

void vp9_default_coef_probs(VP9_COMMON *pc) {
  vpx_memcpy(pc->fc.coef_probs, default_coef_probs,
             sizeof(pc->fc.coef_probs));
  vpx_memcpy(pc->fc.hybrid_coef_probs, default_hybrid_coef_probs,
             sizeof(pc->fc.hybrid_coef_probs));

  vpx_memcpy(pc->fc.coef_probs_8x8, default_coef_probs_8x8,
             sizeof(pc->fc.coef_probs_8x8));
  vpx_memcpy(pc->fc.hybrid_coef_probs_8x8, default_hybrid_coef_probs_8x8,
             sizeof(pc->fc.hybrid_coef_probs_8x8));

  vpx_memcpy(pc->fc.coef_probs_16x16, default_coef_probs_16x16,
             sizeof(pc->fc.coef_probs_16x16));
  vpx_memcpy(pc->fc.hybrid_coef_probs_16x16,
             default_hybrid_coef_probs_16x16,
             sizeof(pc->fc.hybrid_coef_probs_16x16));
}

void vp9_coef_tree_initialize() {
  init_bit_trees();
  vp9_tokens_from_tree(vp9_coef_encodings, vp9_coef_tree);
}

// #define COEF_COUNT_TESTING

#define COEF_COUNT_SAT 24
#define COEF_MAX_UPDATE_FACTOR 112
#define COEF_COUNT_SAT_KEY 24
#define COEF_MAX_UPDATE_FACTOR_KEY 112
#define COEF_COUNT_SAT_AFTER_KEY 24
#define COEF_MAX_UPDATE_FACTOR_AFTER_KEY 128

void vp9_adapt_coef_probs(VP9_COMMON *cm) {
  int t, i, j, k, count;
  unsigned int branch_ct[ENTROPY_NODES][2];
  vp9_prob coef_probs[ENTROPY_NODES];
  int update_factor; /* denominator 256 */
  int factor;
  int count_sat;

  // printf("Frame type: %d\n", cm->frame_type);
  if (cm->frame_type == KEY_FRAME) {
    update_factor = COEF_MAX_UPDATE_FACTOR_KEY;
    count_sat = COEF_COUNT_SAT_KEY;
  } else if (cm->last_frame_type == KEY_FRAME) {
    update_factor = COEF_MAX_UPDATE_FACTOR_AFTER_KEY;  /* adapt quickly */
    count_sat = COEF_COUNT_SAT_AFTER_KEY;
  } else {
    update_factor = COEF_MAX_UPDATE_FACTOR;
    count_sat = COEF_COUNT_SAT;
  }

#ifdef COEF_COUNT_TESTING
  {
    printf("static const unsigned int\ncoef_counts"
           "[BLOCK_TYPES] [COEF_BANDS]"
           "[PREV_COEF_CONTEXTS] [MAX_ENTROPY_TOKENS] = {\n");
    for (i = 0; i < BLOCK_TYPES; ++i) {
      printf("  {\n");
      for (j = 0; j < COEF_BANDS; ++j) {
        printf("    {\n");
        for (k = 0; k < PREV_COEF_CONTEXTS; ++k) {
          printf("      {");
          for (t = 0; t < MAX_ENTROPY_TOKENS; ++t)
            printf("%d, ", cm->fc.coef_counts[i][j][k][t]);
          printf("},\n");
        }
        printf("    },\n");
      }
      printf("  },\n");
    }
    printf("};\n");
    printf("static const unsigned int\ncoef_counts_8x8"
           "[BLOCK_TYPES_8X8] [COEF_BANDS]"
           "[PREV_COEF_CONTEXTS] [MAX_ENTROPY_TOKENS] = {\n");
    for (i = 0; i < BLOCK_TYPES_8X8; ++i) {
      printf("  {\n");
      for (j = 0; j < COEF_BANDS; ++j) {
        printf("    {\n");
        for (k = 0; k < PREV_COEF_CONTEXTS; ++k) {
          printf("      {");
          for (t = 0; t < MAX_ENTROPY_TOKENS; ++t)
            printf("%d, ", cm->fc.coef_counts_8x8[i][j][k][t]);
          printf("},\n");
        }
        printf("    },\n");
      }
      printf("  },\n");
    }
    printf("};\n");
    printf("static const unsigned int\nhybrid_coef_counts"
           "[BLOCK_TYPES] [COEF_BANDS]"
           "[PREV_COEF_CONTEXTS] [MAX_ENTROPY_TOKENS] = {\n");
    for (i = 0; i < BLOCK_TYPES; ++i) {
      printf("  {\n");
      for (j = 0; j < COEF_BANDS; ++j) {
        printf("    {\n");
        for (k = 0; k < PREV_COEF_CONTEXTS; ++k) {
          printf("      {");
          for (t = 0; t < MAX_ENTROPY_TOKENS; ++t)
            printf("%d, ", cm->fc.hybrid_coef_counts[i][j][k][t]);
          printf("},\n");
        }
        printf("    },\n");
      }
      printf("  },\n");
    }
    printf("};\n");
  }
#endif

  for (i = 0; i < BLOCK_TYPES; ++i)
    for (j = 0; j < COEF_BANDS; ++j)
      for (k = 0; k < PREV_COEF_CONTEXTS; ++k) {
        if (k >= 3 && ((i == 0 && j == 1) || (i > 0 && j == 0)))
          continue;
        vp9_tree_probs_from_distribution(
          MAX_ENTROPY_TOKENS, vp9_coef_encodings, vp9_coef_tree,
          coef_probs, branch_ct, cm->fc.coef_counts [i][j][k],
          256, 1);
        for (t = 0; t < ENTROPY_NODES; ++t) {
          int prob;
          count = branch_ct[t][0] + branch_ct[t][1];
          count = count > count_sat ? count_sat : count;
          factor = (update_factor * count / count_sat);
          prob = ((int)cm->fc.pre_coef_probs[i][j][k][t] * (256 - factor) +
                  (int)coef_probs[t] * factor + 128) >> 8;
          if (prob <= 0) cm->fc.coef_probs[i][j][k][t] = 1;
          else if (prob > 255) cm->fc.coef_probs[i][j][k][t] = 255;
          else cm->fc.coef_probs[i][j][k][t] = prob;
        }
      }

  for (i = 0; i < BLOCK_TYPES; ++i)
    for (j = 0; j < COEF_BANDS; ++j)
      for (k = 0; k < PREV_COEF_CONTEXTS; ++k) {
        if (k >= 3 && ((i == 0 && j == 1) || (i > 0 && j == 0)))
          continue;
        vp9_tree_probs_from_distribution(
          MAX_ENTROPY_TOKENS, vp9_coef_encodings, vp9_coef_tree,
          coef_probs, branch_ct, cm->fc.hybrid_coef_counts [i][j][k],
          256, 1);
        for (t = 0; t < ENTROPY_NODES; ++t) {
          int prob;
          count = branch_ct[t][0] + branch_ct[t][1];
          count = count > count_sat ? count_sat : count;
          factor = (update_factor * count / count_sat);
          prob = ((int)cm->fc.pre_hybrid_coef_probs[i][j][k][t] * (256 - factor) +
                  (int)coef_probs[t] * factor + 128) >> 8;
          if (prob <= 0) cm->fc.hybrid_coef_probs[i][j][k][t] = 1;
          else if (prob > 255) cm->fc.hybrid_coef_probs[i][j][k][t] = 255;
          else cm->fc.hybrid_coef_probs[i][j][k][t] = prob;
        }
      }

  for (i = 0; i < BLOCK_TYPES_8X8; ++i)
    for (j = 0; j < COEF_BANDS; ++j)
      for (k = 0; k < PREV_COEF_CONTEXTS; ++k) {
        if (k >= 3 && ((i == 0 && j == 1) || (i > 0 && j == 0)))
          continue;
        vp9_tree_probs_from_distribution(
          MAX_ENTROPY_TOKENS, vp9_coef_encodings, vp9_coef_tree,
          coef_probs, branch_ct, cm->fc.coef_counts_8x8 [i][j][k],
          256, 1);
        for (t = 0; t < ENTROPY_NODES; ++t) {
          int prob;
          count = branch_ct[t][0] + branch_ct[t][1];
          count = count > count_sat ? count_sat : count;
          factor = (update_factor * count / count_sat);
          prob = ((int)cm->fc.pre_coef_probs_8x8[i][j][k][t] * (256 - factor) +
                  (int)coef_probs[t] * factor + 128) >> 8;
          if (prob <= 0) cm->fc.coef_probs_8x8[i][j][k][t] = 1;
          else if (prob > 255) cm->fc.coef_probs_8x8[i][j][k][t] = 255;
          else cm->fc.coef_probs_8x8[i][j][k][t] = prob;
        }
      }

  for (i = 0; i < BLOCK_TYPES_8X8; ++i)
    for (j = 0; j < COEF_BANDS; ++j)
      for (k = 0; k < PREV_COEF_CONTEXTS; ++k) {
        if (k >= 3 && ((i == 0 && j == 1) || (i > 0 && j == 0)))
          continue;
        vp9_tree_probs_from_distribution(
          MAX_ENTROPY_TOKENS, vp9_coef_encodings, vp9_coef_tree,
          coef_probs, branch_ct, cm->fc.hybrid_coef_counts_8x8 [i][j][k],
          256, 1);
        for (t = 0; t < ENTROPY_NODES; ++t) {
          int prob;
          count = branch_ct[t][0] + branch_ct[t][1];
          count = count > count_sat ? count_sat : count;
          factor = (update_factor * count / count_sat);
          prob = ((int)cm->fc.pre_hybrid_coef_probs_8x8[i][j][k][t] *
                  (256 - factor) +
                  (int)coef_probs[t] * factor + 128) >> 8;
          if (prob <= 0) cm->fc.hybrid_coef_probs_8x8[i][j][k][t] = 1;
          else if (prob > 255) cm->fc.hybrid_coef_probs_8x8[i][j][k][t] = 255;
          else cm->fc.hybrid_coef_probs_8x8[i][j][k][t] = prob;
        }
      }

  for (i = 0; i < BLOCK_TYPES_16X16; ++i)
    for (j = 0; j < COEF_BANDS; ++j)
      for (k = 0; k < PREV_COEF_CONTEXTS; ++k) {
        if (k >= 3 && ((i == 0 && j == 1) || (i > 0 && j == 0)))
          continue;
        vp9_tree_probs_from_distribution(
          MAX_ENTROPY_TOKENS, vp9_coef_encodings, vp9_coef_tree,
          coef_probs, branch_ct, cm->fc.coef_counts_16x16[i][j][k], 256, 1);
        for (t = 0; t < ENTROPY_NODES; ++t) {
          int prob;
          count = branch_ct[t][0] + branch_ct[t][1];
          count = count > count_sat ? count_sat : count;
          factor = (update_factor * count / count_sat);
          prob = ((int)cm->fc.pre_coef_probs_16x16[i][j][k][t] *
                  (256 - factor) +
                  (int)coef_probs[t] * factor + 128) >> 8;
          if (prob <= 0) cm->fc.coef_probs_16x16[i][j][k][t] = 1;
          else if (prob > 255) cm->fc.coef_probs_16x16[i][j][k][t] = 255;
          else cm->fc.coef_probs_16x16[i][j][k][t] = prob;
        }
      }

  for (i = 0; i < BLOCK_TYPES_16X16; ++i)
    for (j = 0; j < COEF_BANDS; ++j)
      for (k = 0; k < PREV_COEF_CONTEXTS; ++k) {
        if (k >= 3 && ((i == 0 && j == 1) || (i > 0 && j == 0)))
          continue;
        vp9_tree_probs_from_distribution(
          MAX_ENTROPY_TOKENS, vp9_coef_encodings, vp9_coef_tree,
          coef_probs, branch_ct, cm->fc.hybrid_coef_counts_16x16[i][j][k], 256, 1);
        for (t = 0; t < ENTROPY_NODES; ++t) {
          int prob;
          count = branch_ct[t][0] + branch_ct[t][1];
          count = count > count_sat ? count_sat : count;
          factor = (update_factor * count / count_sat);
          prob = ((int)cm->fc.pre_hybrid_coef_probs_16x16[i][j][k][t] * (256 - factor) +
                  (int)coef_probs[t] * factor + 128) >> 8;
          if (prob <= 0) cm->fc.hybrid_coef_probs_16x16[i][j][k][t] = 1;
          else if (prob > 255) cm->fc.hybrid_coef_probs_16x16[i][j][k][t] = 255;
          else cm->fc.hybrid_coef_probs_16x16[i][j][k][t] = prob;
        }
      }
}
