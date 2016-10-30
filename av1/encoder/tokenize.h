/*
 * Copyright (c) 2016, Alliance for Open Media. All rights reserved
 *
 * This source code is subject to the terms of the BSD 2 Clause License and
 * the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
 * was not distributed with this source code in the LICENSE file, you can
 * obtain it at www.aomedia.org/license/software. If the Alliance for Open
 * Media Patent License 1.0 was not distributed with this source code in the
 * PATENTS file, you can obtain it at www.aomedia.org/license/patent.
 */

#ifndef AV1_ENCODER_TOKENIZE_H_
#define AV1_ENCODER_TOKENIZE_H_

#include "av1/common/entropy.h"

#include "av1/encoder/block.h"
#include "av1/encoder/treewriter.h"

#ifdef __cplusplus
extern "C" {
#endif

#define EOSB_TOKEN 127  // Not signalled, encoder only

#if CONFIG_AOM_HIGHBITDEPTH
typedef int32_t EXTRABIT;
#else
typedef int16_t EXTRABIT;
#endif

typedef struct {
  int16_t token;
  EXTRABIT extra;
} TOKENVALUE;

typedef struct {
  const aom_prob *context_tree;
#if CONFIG_EC_MULTISYMBOL
  aom_cdf_prob (*token_cdf)[ENTROPY_TOKENS];
#endif
  EXTRABIT extra;
  uint8_t token;
  uint8_t skip_eob_node;
} TOKENEXTRA;

extern const aom_tree_index av1_coef_tree[];
extern const aom_tree_index av1_coef_con_tree[];
#if !CONFIG_RANS
extern const struct av1_token av1_coef_encodings[];
#endif  // !CONFIG_RANS

int av1_is_skippable_in_plane(MACROBLOCK *x, BLOCK_SIZE bsize, int plane);
int av1_has_high_freq_in_plane(MACROBLOCK *x, BLOCK_SIZE bsize, int plane);

struct AV1_COMP;
struct ThreadData;

typedef enum {
  OUTPUT_ENABLED = 0,
  DRY_RUN_NORMAL,
  DRY_RUN_COSTCOEFFS,
} RUN_TYPE;

// Note in all the tokenize functions rate if non NULL is incremented
// with the coefficient token cost only if dry_run = DRY_RUN_COSTCOEFS,
// otherwise rate is not incremented.
#if CONFIG_VAR_TX
void av1_tokenize_sb_vartx(const struct AV1_COMP *cpi, struct ThreadData *td,
                           TOKENEXTRA **t, RUN_TYPE dry_run, int mi_row,
                           int mi_col, BLOCK_SIZE bsize, int *rate);
#endif
#if CONFIG_PALETTE
void av1_tokenize_palette_sb(const struct AV1_COMP *cpi,
                             const struct ThreadData *const td, int plane,
                             TOKENEXTRA **t, RUN_TYPE dry_run, BLOCK_SIZE bsize,
                             int *rate);
#endif  // CONFIG_PALETTE
void av1_tokenize_sb(const struct AV1_COMP *cpi, struct ThreadData *td,
                     TOKENEXTRA **t, RUN_TYPE dry_run, BLOCK_SIZE bsize,
                     int *rate);
#if CONFIG_SUPERTX
void av1_tokenize_sb_supertx(const struct AV1_COMP *cpi, struct ThreadData *td,
                             TOKENEXTRA **t, RUN_TYPE dry_run, BLOCK_SIZE bsize,
                             int *rate);
#endif

extern const int16_t *av1_dct_value_cost_ptr;
/* TODO: The Token field should be broken out into a separate char array to
 *  improve cache locality, since it's needed for costing when the rest of the
 *  fields are not.
 */
extern const TOKENVALUE *av1_dct_value_tokens_ptr;
extern const TOKENVALUE *av1_dct_cat_lt_10_value_tokens;
extern const int *av1_dct_cat_lt_10_value_cost;
extern const int16_t av1_cat6_low_cost[256];
extern const int av1_cat6_high_cost[64];
extern const int av1_cat6_high10_high_cost[256];
extern const int av1_cat6_high12_high_cost[1024];
static INLINE int av1_get_cost(int16_t token, EXTRABIT extrabits,
                               const int *cat6_high_table) {
  if (token != CATEGORY6_TOKEN)
    return av1_extra_bits[token].cost[extrabits >> 1];
  return av1_cat6_low_cost[(extrabits >> 1) & 0xff] +
         cat6_high_table[extrabits >> 9];
}

#if CONFIG_AOM_HIGHBITDEPTH
static INLINE const int *av1_get_high_cost_table(int bit_depth) {
  return bit_depth == 8 ? av1_cat6_high_cost
                        : (bit_depth == 10 ? av1_cat6_high10_high_cost
                                           : av1_cat6_high12_high_cost);
}
#else
static INLINE const int *av1_get_high_cost_table(int bit_depth) {
  (void)bit_depth;
  return av1_cat6_high_cost;
}
#endif  // CONFIG_AOM_HIGHBITDEPTH

static INLINE void av1_get_token_extra(int v, int16_t *token, EXTRABIT *extra) {
  if (v >= CAT6_MIN_VAL || v <= -CAT6_MIN_VAL) {
    *token = CATEGORY6_TOKEN;
    if (v >= CAT6_MIN_VAL)
      *extra = 2 * v - 2 * CAT6_MIN_VAL;
    else
      *extra = -2 * v - 2 * CAT6_MIN_VAL + 1;
    return;
  }
  *token = av1_dct_cat_lt_10_value_tokens[v].token;
  *extra = av1_dct_cat_lt_10_value_tokens[v].extra;
}
static INLINE int16_t av1_get_token(int v) {
  if (v >= CAT6_MIN_VAL || v <= -CAT6_MIN_VAL) return 10;
  return av1_dct_cat_lt_10_value_tokens[v].token;
}

static INLINE int av1_get_token_cost(int v, int16_t *token,
                                     const int *cat6_high_table) {
  if (v >= CAT6_MIN_VAL || v <= -CAT6_MIN_VAL) {
    EXTRABIT extrabits;
    *token = CATEGORY6_TOKEN;
    extrabits = abs(v) - CAT6_MIN_VAL;
    return av1_cat6_low_cost[extrabits & 0xff] +
           cat6_high_table[extrabits >> 8];
  }
  *token = av1_dct_cat_lt_10_value_tokens[v].token;
  return av1_dct_cat_lt_10_value_cost[v];
}

#ifdef __cplusplus
}  // extern "C"
#endif

#endif  // AV1_ENCODER_TOKENIZE_H_
