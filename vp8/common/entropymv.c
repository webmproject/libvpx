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

#if CONFIG_HIGH_PRECISION_MV
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
#endif  /* CONFIG_HIGH_PRECISION_MV */

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

#if CONFIG_HIGH_PRECISION_MV
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
#endif  /* CONFIG_HIGH_PRECISION_MV */

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

#if CONFIG_HIGH_PRECISION_MV
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
#endif  /* CONFIG_HIGH_PRECISION_MV */

void vp8_entropy_mv_init() {
  vp8_tokens_from_tree(vp8_small_mvencodings, vp8_small_mvtree);
#if CONFIG_HIGH_PRECISION_MV
  vp8_tokens_from_tree(vp8_small_mvencodings_hp, vp8_small_mvtree_hp);
#endif
}

#if CONFIG_ADAPTIVE_ENTROPY
// #define MV_COUNT_TESTING
#define MV_COUNT_SAT 16
#define MV_MAX_UPDATE_FACTOR 128
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
#if CONFIG_HIGH_PRECISION_MV
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
#endif
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
#if CONFIG_HIGH_PRECISION_MV
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
#endif
}
#endif  /* CONFIG_ADAPTIVE_ENTROPY */
