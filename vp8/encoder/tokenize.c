/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */


#include <math.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include "onyx_int.h"
#include "tokenize.h"
#include "vpx_mem/vpx_mem.h"

#include "vp8/common/pred_common.h"
#include "vp8/common/seg_common.h"
#include "vp8/common/entropy.h"

/* Global event counters used for accumulating statistics across several
   compressions, then generating context.c = initial stats. */

#ifdef ENTROPY_STATS
INT64 context_counters[BLOCK_TYPES] [COEF_BANDS] [PREV_COEF_CONTEXTS] [MAX_ENTROPY_TOKENS];
INT64 hybrid_context_counters[BLOCK_TYPES] [COEF_BANDS] [PREV_COEF_CONTEXTS] [MAX_ENTROPY_TOKENS];

INT64 context_counters_8x8[BLOCK_TYPES_8X8] [COEF_BANDS] [PREV_COEF_CONTEXTS] [MAX_ENTROPY_TOKENS];
INT64 hybrid_context_counters_8x8[BLOCK_TYPES_8X8] [COEF_BANDS] [PREV_COEF_CONTEXTS] [MAX_ENTROPY_TOKENS];

INT64 context_counters_16x16[BLOCK_TYPES_16X16] [COEF_BANDS] [PREV_COEF_CONTEXTS] [MAX_ENTROPY_TOKENS];
INT64 hybrid_context_counters_16x16[BLOCK_TYPES_16X16] [COEF_BANDS] [PREV_COEF_CONTEXTS] [MAX_ENTROPY_TOKENS];

extern unsigned int tree_update_hist[BLOCK_TYPES][COEF_BANDS]
                    [PREV_COEF_CONTEXTS][ENTROPY_NODES][2];
extern unsigned int hybrid_tree_update_hist[BLOCK_TYPES][COEF_BANDS]
                    [PREV_COEF_CONTEXTS][ENTROPY_NODES][2];
extern unsigned int tree_update_hist_8x8[BLOCK_TYPES_8X8][COEF_BANDS]
                    [PREV_COEF_CONTEXTS][ENTROPY_NODES] [2];
extern unsigned int hybrid_tree_update_hist_8x8[BLOCK_TYPES_8X8][COEF_BANDS]
                    [PREV_COEF_CONTEXTS][ENTROPY_NODES] [2];
extern unsigned int tree_update_hist_16x16[BLOCK_TYPES_16X16][COEF_BANDS]
                    [PREV_COEF_CONTEXTS][ENTROPY_NODES] [2];
extern unsigned int hybrid_tree_update_hist_16x16[BLOCK_TYPES_16X16][COEF_BANDS]
                    [PREV_COEF_CONTEXTS][ENTROPY_NODES] [2];
#endif  /* ENTROPY_STATS */

void vp8_stuff_mb(VP8_COMP *cpi, MACROBLOCKD *xd, TOKENEXTRA **t, int dry_run);
void vp8_fix_contexts(MACROBLOCKD *xd);

static TOKENVALUE dct_value_tokens[DCT_MAX_VALUE * 2];
const TOKENVALUE *vp8_dct_value_tokens_ptr;
static int dct_value_cost[DCT_MAX_VALUE * 2];
const int *vp8_dct_value_cost_ptr;

static void fill_value_tokens() {

  TOKENVALUE *const t = dct_value_tokens + DCT_MAX_VALUE;
  vp8_extra_bit_struct *const e = vp8_extra_bits;

  int i = -DCT_MAX_VALUE;
  int sign = 1;

  do {
    if (!i)
      sign = 0;

    {
      const int a = sign ? -i : i;
      int eb = sign;

      if (a > 4) {
        int j = 4;

        while (++j < 11  &&  e[j].base_val <= a) {}

        t[i].Token = --j;
        eb |= (a - e[j].base_val) << 1;
      } else
        t[i].Token = a;

      t[i].Extra = eb;
    }

    // initialize the cost for extra bits for all possible coefficient value.
    {
      int cost = 0;
      vp8_extra_bit_struct *p = vp8_extra_bits + t[i].Token;

      if (p->base_val) {
        const int extra = t[i].Extra;
        const int Length = p->Len;

        if (Length)
          cost += vp8_treed_cost(p->tree, p->prob, extra >> 1, Length);

        cost += vp8_cost_bit(vp8_prob_half, extra & 1); /* sign */
        dct_value_cost[i + DCT_MAX_VALUE] = cost;
      }

    }

  } while (++i < DCT_MAX_VALUE);

  vp8_dct_value_tokens_ptr = dct_value_tokens + DCT_MAX_VALUE;
  vp8_dct_value_cost_ptr   = dct_value_cost + DCT_MAX_VALUE;
}

static void tokenize1st_order_b_16x16(MACROBLOCKD *xd,
                                      const BLOCKD *const b,
                                      TOKENEXTRA **tp,
                                      PLANE_TYPE type,
                                      ENTROPY_CONTEXT *a,
                                      ENTROPY_CONTEXT *l,
                                      VP8_COMP *cpi,
                                      int dry_run) {
  int pt; /* near block/prev token context index */
  int c = (type == PLANE_TYPE_Y_NO_DC) ? 1 : 0;
  const int eob = b->eob;     /* one beyond last nonzero coeff */
  TOKENEXTRA *t = *tp;        /* store tokens starting here */
  const short *qcoeff_ptr = b->qcoeff;
  TX_TYPE tx_type = get_tx_type(xd, b);
  int seg_eob = 256;
  int segment_id = xd->mode_info_context->mbmi.segment_id;

  if (segfeature_active(xd, segment_id, SEG_LVL_EOB))
    seg_eob = get_segdata(xd, segment_id, SEG_LVL_EOB);

  VP8_COMBINEENTROPYCONTEXTS(pt, *a, *l);

  do {
    const int band = vp8_coef_bands_16x16[c];
    int x;

    if (c < eob) {
      const int rc = vp8_default_zig_zag1d_16x16[c];
      const int v = qcoeff_ptr[rc];

      assert(-DCT_MAX_VALUE <= v  &&  v < (DCT_MAX_VALUE));

      t->Extra = vp8_dct_value_tokens_ptr[v].Extra;
      x        = vp8_dct_value_tokens_ptr[v].Token;
    } else {
      x = DCT_EOB_TOKEN;
    }

    t->Token = x;
    if (tx_type != DCT_DCT)
      t->context_tree = cpi->common.fc.hybrid_coef_probs_16x16[type][band][pt];
    else
      t->context_tree = cpi->common.fc.coef_probs_16x16[type][band][pt];

    t->skip_eob_node = pt == 0 && ((band > 0 && type != PLANE_TYPE_Y_NO_DC) ||
                                   (band > 1 && type == PLANE_TYPE_Y_NO_DC));
    assert(vp8_coef_encodings[t->Token].Len - t->skip_eob_node > 0);
    if (!dry_run) {
      if (tx_type != DCT_DCT)
        ++cpi->hybrid_coef_counts_16x16[type][band][pt][x];
      else
        ++cpi->coef_counts_16x16[type][band][pt][x];
    }
    pt = vp8_prev_token_class[x];
    ++t;
  } while (c < eob  &&  ++c < seg_eob);

  *tp = t;
  pt = (c != !type); /* 0 <-> all coeff data is zero */
  *a = *l = pt;
}

static void tokenize2nd_order_b_8x8(MACROBLOCKD *xd,
                                    const BLOCKD *const b,
                                    TOKENEXTRA **tp,
                                    ENTROPY_CONTEXT *a,
                                    ENTROPY_CONTEXT *l,
                                    VP8_COMP *cpi,
                                    int dry_run) {
  int pt; /* near block/prev token context index */
  int c = 0;          /* start at DC */
  const int eob = b->eob;     /* one beyond last nonzero coeff */
  TOKENEXTRA *t = *tp;        /* store tokens starting here */
  const short *qcoeff_ptr = b->qcoeff;
  int seg_eob = 4;
  int segment_id = xd->mode_info_context->mbmi.segment_id;

  if (segfeature_active(xd, segment_id, SEG_LVL_EOB)) {
    seg_eob = get_segdata(xd, segment_id, SEG_LVL_EOB);
  }

  VP8_COMBINEENTROPYCONTEXTS(pt, *a, *l);

  assert(eob <= 4);

  do {
    const int band = vp8_coef_bands[c];
    int x;

    if (c < eob) {
      const int rc = vp8_default_zig_zag1d[c];
      const int v = qcoeff_ptr[rc];

      assert(-DCT_MAX_VALUE <= v  &&  v < (DCT_MAX_VALUE));

      t->Extra = vp8_dct_value_tokens_ptr[v].Extra;
      x        = vp8_dct_value_tokens_ptr[v].Token;
    } else {
      x = DCT_EOB_TOKEN;
    }

    t->Token = x;
    t->context_tree = cpi->common.fc.coef_probs_8x8[PLANE_TYPE_Y2][band][pt];

    t->skip_eob_node = ((pt == 0) && (band > 0));
    assert(vp8_coef_encodings[t->Token].Len - t->skip_eob_node > 0);

    if (!dry_run)
      ++cpi->coef_counts_8x8[PLANE_TYPE_Y2][band][pt][x];
    pt = vp8_prev_token_class[x];
    ++t;
  } while (c < eob && ++c < seg_eob);

  *tp = t;
  pt = (c != 0); /* 0 <-> all coeff data is zero */
  *a = *l = pt;
}

static void tokenize2nd_order_b_4x4(MACROBLOCKD *xd,
                                    TOKENEXTRA **tp,
                                    VP8_COMP *cpi,
                                    int dry_run) {
  int pt;             /* near block/prev token context index */
  int c = 0;          /* start at DC */
  TOKENEXTRA *t = *tp;/* store tokens starting here */
  const BLOCKD *b = xd->block + 24;
  const short *qcoeff_ptr = b->qcoeff;
  ENTROPY_CONTEXT *a;
  ENTROPY_CONTEXT *l;
  const int eob = b->eob;
  int seg_eob = 16;
  int segment_id = xd->mode_info_context->mbmi.segment_id;

  if (segfeature_active(xd, segment_id, SEG_LVL_EOB))
    seg_eob = get_segdata(xd, segment_id, SEG_LVL_EOB);

  a = (ENTROPY_CONTEXT *)xd->above_context + 8;
  l = (ENTROPY_CONTEXT *)xd->left_context + 8;

  VP8_COMBINEENTROPYCONTEXTS(pt, *a, *l);

  do {
    const int band = vp8_coef_bands[c];
    int token;

    if (c < eob) {
      const int rc = vp8_default_zig_zag1d[c];
      const int v = qcoeff_ptr[rc];

      t->Extra = vp8_dct_value_tokens_ptr[v].Extra;
      token    = vp8_dct_value_tokens_ptr[v].Token;
    } else
      token    = DCT_EOB_TOKEN;

    t->Token = token;
    t->context_tree = cpi->common.fc.coef_probs[PLANE_TYPE_Y2][band][pt];

    t->skip_eob_node = ((pt == 0) && (band > 0));
    assert(vp8_coef_encodings[t->Token].Len - t->skip_eob_node > 0);

    if (!dry_run)
      ++cpi->coef_counts[PLANE_TYPE_Y2][band][pt][token];
    pt = vp8_prev_token_class[token];
    ++t;
  } while (c < eob && ++c < seg_eob);

  *tp = t;
  pt = (c != 0); /* 0 <-> all coeff data is zero */
  *a = *l = pt;
}

static void tokenize1st_order_b_8x8(MACROBLOCKD *xd,
                                    const BLOCKD *const b,
                                    TOKENEXTRA **tp,
                                    PLANE_TYPE type,
                                    ENTROPY_CONTEXT *a,
                                    ENTROPY_CONTEXT *l,
                                    VP8_COMP *cpi,
                                    int dry_run) {
  int pt; /* near block/prev token context index */
  int c = (type == PLANE_TYPE_Y_NO_DC) ? 1 : 0; /* start at DC unless type 0 */
  TOKENEXTRA *t = *tp;        /* store tokens starting here */
  const short *qcoeff_ptr = b->qcoeff;
  TX_TYPE tx_type = get_tx_type(xd, b);
  const int eob = b->eob;
  int seg_eob = 64;
  int segment_id = xd->mode_info_context->mbmi.segment_id;

  if (segfeature_active(xd, segment_id, SEG_LVL_EOB))
    seg_eob = get_segdata(xd, segment_id, SEG_LVL_EOB);

  VP8_COMBINEENTROPYCONTEXTS(pt, *a, *l);

  do {
    const int band = vp8_coef_bands_8x8[c];
    int x;

    if (c < eob) {
      const int rc = vp8_default_zig_zag1d_8x8[c];
      const int v = qcoeff_ptr[rc];

      assert(-DCT_MAX_VALUE <= v  &&  v < (DCT_MAX_VALUE));

      t->Extra = vp8_dct_value_tokens_ptr[v].Extra;
      x        = vp8_dct_value_tokens_ptr[v].Token;
    } else
      x = DCT_EOB_TOKEN;

    t->Token = x;
    if (tx_type != DCT_DCT)
      t->context_tree = cpi->common.fc.hybrid_coef_probs_8x8[type][band][pt];
    else
      t->context_tree = cpi->common.fc.coef_probs_8x8[type][band][pt];

    t->skip_eob_node = pt == 0 && ((band > 0 && type != PLANE_TYPE_Y_NO_DC) ||
                                   (band > 1 && type == PLANE_TYPE_Y_NO_DC));
    assert(vp8_coef_encodings[t->Token].Len - t->skip_eob_node > 0);

    if (!dry_run) {
      if (tx_type != DCT_DCT)
        ++cpi->hybrid_coef_counts_8x8[type][band][pt][x];
      else
        ++cpi->coef_counts_8x8[type][band][pt][x];
    }
    pt = vp8_prev_token_class[x];
    ++t;
  } while (c < eob && ++c < seg_eob);

  *tp = t;
  pt = (c != !type); /* 0 <-> all coeff data is zero */
  *a = *l = pt;
}

static void tokenize1st_order_chroma_4x4(MACROBLOCKD *xd,
                                         TOKENEXTRA **tp,
                                         VP8_COMP *cpi,
                                         int dry_run) {
  unsigned int block;
  const BLOCKD *b = xd->block + 16;
  int pt;             /* near block/prev token context index */
  TOKENEXTRA *t = *tp;/* store tokens starting here */
  ENTROPY_CONTEXT *a;
  ENTROPY_CONTEXT *l;
  int seg_eob = 16;
  int segment_id = xd->mode_info_context->mbmi.segment_id;

  if (segfeature_active(xd, segment_id, SEG_LVL_EOB)) {
    seg_eob = get_segdata(xd, segment_id, SEG_LVL_EOB);
  }

  /* Chroma */
  for (block = 16; block < 24; block++, b++) {
    const int eob = b->eob;
    const int tmp1 = vp8_block2above[block];
    const int tmp2 = vp8_block2left[block];
    const int16_t *qcoeff_ptr = b->qcoeff;
    int c = 0;

    a = (ENTROPY_CONTEXT *)xd->above_context + tmp1;
    l = (ENTROPY_CONTEXT *)xd->left_context + tmp2;

    VP8_COMBINEENTROPYCONTEXTS(pt, *a, *l);

    do {
      const int band = vp8_coef_bands[c];
      int token;

      if (c < eob) {
        const int rc = vp8_default_zig_zag1d[c];
        const int v = qcoeff_ptr[rc];

        t->Extra = vp8_dct_value_tokens_ptr[v].Extra;
        token    = vp8_dct_value_tokens_ptr[v].Token;
      } else
        token = DCT_EOB_TOKEN;

      t->Token = token;
      t->context_tree = cpi->common.fc.coef_probs[PLANE_TYPE_UV][band][pt];

      t->skip_eob_node = ((pt == 0) && (band > 0));
      assert(vp8_coef_encodings[t->Token].Len - t->skip_eob_node > 0);

      if (!dry_run)
        ++cpi->coef_counts[PLANE_TYPE_UV][band][pt][token];
      pt = vp8_prev_token_class[token];
      ++t;
    } while (c < eob && ++c < seg_eob);

    *tp = t;
    pt = (c != 0); /* 0 <-> all coeff data is zero */
    *a = *l = pt;
  }
}

static void tokenize1st_order_b_4x4(MACROBLOCKD *xd,
                                    TOKENEXTRA **tp,
                                    PLANE_TYPE type,
                                    VP8_COMP *cpi,
                                    int dry_run) {
  unsigned int block;
  const BLOCKD *b = xd->block;
  int pt;             /* near block/prev token context index */
  TOKENEXTRA *t = *tp;/* store tokens starting here */
  ENTROPY_CONTEXT *a, *l;
  int seg_eob = 16;
  int segment_id = xd->mode_info_context->mbmi.segment_id;
  int const *pt_scan = vp8_default_zig_zag1d;

  if (segfeature_active(xd, segment_id, SEG_LVL_EOB)) {
    seg_eob = get_segdata(xd, segment_id, SEG_LVL_EOB);
  }

  /* Luma */
  for (block = 0; block < 16; block++, b++) {
    const int eob = b->eob;
    const int16_t *qcoeff_ptr = b->qcoeff;
    int c = (type == PLANE_TYPE_Y_NO_DC) ? 1 : 0;

    TX_TYPE tx_type = get_tx_type(xd, &xd->block[block]);
    switch (tx_type) {
      case ADST_DCT:
        pt_scan = vp8_row_scan;
        break;
      case DCT_ADST:
        pt_scan = vp8_col_scan;
        break;
      default :
        pt_scan = vp8_default_zig_zag1d;
        break;
    }
    a = (ENTROPY_CONTEXT *)xd->above_context + vp8_block2above[block];
    l = (ENTROPY_CONTEXT *)xd->left_context + vp8_block2left[block];
    VP8_COMBINEENTROPYCONTEXTS(pt, *a, *l);

    assert(b->eob <= 16);

    do {
      const int band = vp8_coef_bands[c];
      int token;

      if (c < eob) {
        const int rc = pt_scan[c];
        const int v = qcoeff_ptr[rc];

        t->Extra = vp8_dct_value_tokens_ptr[v].Extra;
        token    = vp8_dct_value_tokens_ptr[v].Token;
      } else
        token = DCT_EOB_TOKEN;

      t->Token = token;
      if (tx_type != DCT_DCT)
        t->context_tree = cpi->common.fc.hybrid_coef_probs[type][band][pt];
      else
        t->context_tree = cpi->common.fc.coef_probs[type][band][pt];

      t->skip_eob_node = pt == 0 && ((band > 0 && type != PLANE_TYPE_Y_NO_DC) ||
                                     (band > 1 && type == PLANE_TYPE_Y_NO_DC));
      assert(vp8_coef_encodings[t->Token].Len - t->skip_eob_node > 0);
      if (!dry_run) {
        if (tx_type != DCT_DCT)
          ++cpi->hybrid_coef_counts[type][band][pt][token];
        else
          ++cpi->coef_counts[type][band][pt][token];
      }
      pt = vp8_prev_token_class[token];
      ++t;
    } while (c < eob && ++c < seg_eob);

    *tp = t;
    pt = (c != !type); /* 0 <-> all coeff data is zero */
    *a = *l = pt;
  }

  tokenize1st_order_chroma_4x4(xd, tp, cpi, dry_run);
}

int mby_is_skippable_4x4(MACROBLOCKD *xd, int has_y2_block) {
  int skip = 1;
  int i = 0;

  if (has_y2_block) {
    for (i = 0; i < 16; i++)
      skip &= (xd->block[i].eob < 2);
    skip &= (!xd->block[24].eob);
  } else {
    for (i = 0; i < 16; i++)
      skip &= (!xd->block[i].eob);
  }
  return skip;
}

int mbuv_is_skippable_4x4(MACROBLOCKD *xd) {
  int skip = 1;
  int i;

  for (i = 16; i < 24; i++)
    skip &= (!xd->block[i].eob);
  return skip;
}

int mb_is_skippable_4x4(MACROBLOCKD *xd, int has_y2_block) {
  return (mby_is_skippable_4x4(xd, has_y2_block) &
          mbuv_is_skippable_4x4(xd));
}

int mby_is_skippable_8x8(MACROBLOCKD *xd, int has_y2_block) {
  int skip = 1;
  int i = 0;

  if (has_y2_block) {
    for (i = 0; i < 16; i += 4)
      skip &= (xd->block[i].eob < 2);
    skip &= (!xd->block[24].eob);
  } else {
    for (i = 0; i < 16; i += 4)
      skip &= (!xd->block[i].eob);
  }
  return skip;
}

int mbuv_is_skippable_8x8(MACROBLOCKD *xd) {
  return (!xd->block[16].eob) & (!xd->block[20].eob);
}

int mb_is_skippable_8x8(MACROBLOCKD *xd, int has_y2_block) {
  return (mby_is_skippable_8x8(xd, has_y2_block) &
          mbuv_is_skippable_8x8(xd));
}

int mb_is_skippable_8x8_4x4uv(MACROBLOCKD *xd, int has_y2_block) {
  return (mby_is_skippable_8x8(xd, has_y2_block) &
          mbuv_is_skippable_4x4(xd));
}

int mby_is_skippable_16x16(MACROBLOCKD *xd) {
  int skip = 1;
  //skip &= (xd->block[0].eob < 2); // I think this should be commented? No second order == DC must be coded
  //skip &= (xd->block[0].eob < 1);
  //skip &= (!xd->block[24].eob);
  skip &= !xd->block[0].eob;
  return skip;
}

int mb_is_skippable_16x16(MACROBLOCKD *xd) {
  return (mby_is_skippable_16x16(xd) & mbuv_is_skippable_8x8(xd));
}

void vp8_tokenize_mb(VP8_COMP *cpi,
                     MACROBLOCKD *xd,
                     TOKENEXTRA **t,
                     int dry_run) {
  PLANE_TYPE plane_type;
  int has_y2_block;
  int b;
  int tx_size = xd->mode_info_context->mbmi.txfm_size;
  int mb_skip_context = get_pred_context(&cpi->common, xd, PRED_MBSKIP);
  TOKENEXTRA *t_backup = *t;

  // If the MB is going to be skipped because of a segment level flag
  // exclude this from the skip count stats used to calculate the
  // transmitted skip probability;
  int skip_inc;
  int segment_id = xd->mode_info_context->mbmi.segment_id;

  if (!segfeature_active(xd, segment_id, SEG_LVL_EOB) ||
      (get_segdata(xd, segment_id, SEG_LVL_EOB) != 0)) {
    skip_inc = 1;
  } else
    skip_inc = 0;

  has_y2_block = (tx_size != TX_16X16
                  && xd->mode_info_context->mbmi.mode != B_PRED
                  && xd->mode_info_context->mbmi.mode != I8X8_PRED
                  && xd->mode_info_context->mbmi.mode != SPLITMV);

  switch (tx_size) {
    case TX_16X16:
      xd->mode_info_context->mbmi.mb_skip_coeff = mb_is_skippable_16x16(xd);
      break;
    case TX_8X8:
      if (xd->mode_info_context->mbmi.mode == I8X8_PRED)
        xd->mode_info_context->mbmi.mb_skip_coeff = mb_is_skippable_8x8_4x4uv(xd, 0);
      else
        xd->mode_info_context->mbmi.mb_skip_coeff = mb_is_skippable_8x8(xd, has_y2_block);
      break;

    default:
      xd->mode_info_context->mbmi.mb_skip_coeff = mb_is_skippable_4x4(xd, has_y2_block);
      break;
  }

  if (xd->mode_info_context->mbmi.mb_skip_coeff) {
    if (!dry_run)
      cpi->skip_true_count[mb_skip_context] += skip_inc;
    if (!cpi->common.mb_no_coeff_skip) {
      vp8_stuff_mb(cpi, xd, t, dry_run);
    } else {
      vp8_fix_contexts(xd);
    }
    if (dry_run)
      *t = t_backup;
    return;
  }

  if (!dry_run)
    cpi->skip_false_count[mb_skip_context] += skip_inc;

  if (has_y2_block) {
    if (tx_size == TX_8X8) {
      ENTROPY_CONTEXT *A = (ENTROPY_CONTEXT *)xd->above_context;
      ENTROPY_CONTEXT *L = (ENTROPY_CONTEXT *)xd->left_context;
      tokenize2nd_order_b_8x8(xd,
                              xd->block + 24, t,
                              A + vp8_block2above_8x8[24],
                              L + vp8_block2left_8x8[24],
                              cpi, dry_run);
    } else
      tokenize2nd_order_b_4x4(xd, t, cpi, dry_run);

    plane_type = PLANE_TYPE_Y_NO_DC;
  } else
    plane_type = PLANE_TYPE_Y_WITH_DC;

  if (tx_size == TX_16X16) {
    ENTROPY_CONTEXT * A = (ENTROPY_CONTEXT *)xd->above_context;
    ENTROPY_CONTEXT * L = (ENTROPY_CONTEXT *)xd->left_context;

    tokenize1st_order_b_16x16(xd, xd->block, t, PLANE_TYPE_Y_WITH_DC,
                              A, L, cpi, dry_run);

    for (b = 1; b < 16; b++) {
      *(A + vp8_block2above[b]) = *(A);
      *(L + vp8_block2left[b] ) = *(L);
    }
    for (b = 16; b < 24; b += 4) {
      tokenize1st_order_b_8x8(xd, xd->block + b, t, PLANE_TYPE_UV,
                              A + vp8_block2above_8x8[b],
                              L + vp8_block2left_8x8[b], cpi, dry_run);
      *(A + vp8_block2above_8x8[b]+1) = *(A + vp8_block2above_8x8[b]);
      *(L + vp8_block2left_8x8[b]+1 ) = *(L + vp8_block2left_8x8[b]);
    }
    vpx_memset(&A[8], 0, sizeof(A[8]));
    vpx_memset(&L[8], 0, sizeof(L[8]));
  }
  else if (tx_size == TX_8X8) {
    ENTROPY_CONTEXT *A = (ENTROPY_CONTEXT *)xd->above_context;
    ENTROPY_CONTEXT *L = (ENTROPY_CONTEXT *)xd->left_context;
    for (b = 0; b < 16; b += 4) {
      tokenize1st_order_b_8x8(xd,
                              xd->block + b, t, plane_type,
                              A + vp8_block2above_8x8[b],
                              L + vp8_block2left_8x8[b],
                              cpi, dry_run);
      *(A + vp8_block2above_8x8[b] + 1) = *(A + vp8_block2above_8x8[b]);
      *(L + vp8_block2left_8x8[b] + 1)  = *(L + vp8_block2left_8x8[b]);
    }
    if (xd->mode_info_context->mbmi.mode == I8X8_PRED) {
      tokenize1st_order_chroma_4x4(xd, t, cpi, dry_run);
    } else {
      for (b = 16; b < 24; b += 4) {
        tokenize1st_order_b_8x8(xd, xd->block + b, t, PLANE_TYPE_UV,
                                A + vp8_block2above_8x8[b],
                                L + vp8_block2left_8x8[b], cpi, dry_run);
        *(A + vp8_block2above_8x8[b] + 1) = *(A + vp8_block2above_8x8[b]);
        *(L + vp8_block2left_8x8[b] + 1) = *(L + vp8_block2left_8x8[b]);
      }
    }
  } else {
    tokenize1st_order_b_4x4(xd, t, plane_type, cpi, dry_run);
  }
  if (dry_run)
    *t = t_backup;
}


#ifdef ENTROPY_STATS
void init_context_counters(void) {
  FILE *f = fopen("context.bin", "rb");
  if (!f) {
    vpx_memset(context_counters, 0, sizeof(context_counters));
    vpx_memset(context_counters_8x8, 0, sizeof(context_counters_8x8));
    vpx_memset(context_counters_16x16, 0, sizeof(context_counters_16x16));
  } else {
    fread(context_counters, sizeof(context_counters), 1, f);
    fread(context_counters_8x8, sizeof(context_counters_8x8), 1, f);
    fread(context_counters_16x16, sizeof(context_counters_16x16), 1, f);
    fclose(f);
  }

  f = fopen("treeupdate.bin", "rb");
  if (!f) {
    vpx_memset(tree_update_hist, 0, sizeof(tree_update_hist));
    vpx_memset(tree_update_hist_8x8, 0, sizeof(tree_update_hist_8x8));
    vpx_memset(tree_update_hist_16x16, 0, sizeof(tree_update_hist_16x16));
  } else {
    fread(tree_update_hist, sizeof(tree_update_hist), 1, f);
    fread(tree_update_hist_8x8, sizeof(tree_update_hist_8x8), 1, f);
    fread(tree_update_hist_16x16, sizeof(tree_update_hist_16x16), 1, f);
    fclose(f);
  }
}

void print_context_counters() {
  int type, band, pt, t;
  FILE *f = fopen("context.c", "w");

  fprintf(f, "#include \"entropy.h\"\n");
  fprintf(f, "\n/* *** GENERATED FILE: DO NOT EDIT *** */\n\n");
  fprintf(f, "static const unsigned int\n"
          "vp8_default_coef_counts[BLOCK_TYPES]\n"
          "                      [COEF_BANDS]\n"
          "                      [PREV_COEF_CONTEXTS]\n"
          "                      [MAX_ENTROPY_TOKENS]={\n");

# define Comma( X) (X? ",":"")
  type = 0;
  do {
    fprintf(f, "%s\n  { /* block Type %d */", Comma(type), type);
    band = 0;
    do {
      fprintf(f, "%s\n    { /* Coeff Band %d */", Comma(band), band);
      pt = 0;
      do {
        fprintf(f, "%s\n      {", Comma(pt));

        t = 0;
        do {
          const INT64 x = context_counters [type] [band] [pt] [t];
          const int y = (int) x;
          assert(x == (INT64) y);  /* no overflow handling yet */
          fprintf(f, "%s %d", Comma(t), y);
        } while (++t < MAX_ENTROPY_TOKENS);
        fprintf(f, "}");
      } while (++pt < PREV_COEF_CONTEXTS);
      fprintf(f, "\n    }");
    } while (++band < COEF_BANDS);
    fprintf(f, "\n  }");
  } while (++type < BLOCK_TYPES);
  fprintf(f, "\n};\n");

  fprintf(f, "static const unsigned int\nvp8_default_coef_counts_8x8"
          "[BLOCK_TYPES_8X8] [COEF_BANDS]"
          "[PREV_COEF_CONTEXTS] [MAX_ENTROPY_TOKENS] = {");
  type = 0;
  do {
    fprintf(f, "%s\n  { /* block Type %d */", Comma(type), type);
    band = 0;
    do {
      fprintf(f, "%s\n    { /* Coeff Band %d */", Comma(band), band);
      pt = 0;
      do {
        fprintf(f, "%s\n      {", Comma(pt));
        t = 0;
        do {
          const INT64 x = context_counters_8x8 [type] [band] [pt] [t];
          const int y = (int) x;

          assert(x == (INT64) y);  /* no overflow handling yet */
          fprintf(f, "%s %d", Comma(t), y);

        } while (++t < MAX_ENTROPY_TOKENS);

        fprintf(f, "}");
      } while (++pt < PREV_COEF_CONTEXTS);

      fprintf(f, "\n    }");

    } while (++band < COEF_BANDS);

    fprintf(f, "\n  }");
  } while (++type < BLOCK_TYPES_8X8);
  fprintf(f, "\n};\n");

  fprintf(f, "static const unsigned int\nvp8_default_coef_counts_16x16"
          "[BLOCK_TYPES_16X16] [COEF_BANDS]"
          "[PREV_COEF_CONTEXTS] [MAX_ENTROPY_TOKENS] = {");
  type = 0;
  do {
    fprintf(f, "%s\n  { /* block Type %d */", Comma(type), type);
    band = 0;
    do {
      fprintf(f, "%s\n    { /* Coeff Band %d */", Comma(band), band);
      pt = 0;
      do {
        fprintf(f, "%s\n      {", Comma(pt));
        t = 0;
        do {
          const INT64 x = context_counters_16x16 [type] [band] [pt] [t];
          const int y = (int) x;

          assert(x == (INT64) y);  /* no overflow handling yet */
          fprintf(f, "%s %d", Comma(t), y);

        } while (++t < MAX_ENTROPY_TOKENS);

        fprintf(f, "}");
      } while (++pt < PREV_COEF_CONTEXTS);

      fprintf(f, "\n    }");

    } while (++band < COEF_BANDS);

    fprintf(f, "\n  }");
  } while (++type < BLOCK_TYPES_16X16);
  fprintf(f, "\n};\n");

  fprintf(f, "static const vp8_prob\n"
          "vp8_default_coef_probs[BLOCK_TYPES] [COEF_BANDS] \n"
          "[PREV_COEF_CONTEXTS] [ENTROPY_NODES] = {");
  type = 0;
  do {
    fprintf(f, "%s\n  { /* block Type %d */", Comma(type), type);
    band = 0;
    do {
      fprintf(f, "%s\n    { /* Coeff Band %d */", Comma(band), band);
      pt = 0;
      do {
        unsigned int branch_ct [ENTROPY_NODES] [2];
        unsigned int coef_counts[MAX_ENTROPY_TOKENS];
        vp8_prob coef_probs[ENTROPY_NODES];
        for (t = 0; t < MAX_ENTROPY_TOKENS; ++t)
          coef_counts[t] = context_counters [type] [band] [pt] [t];
        vp8_tree_probs_from_distribution(
          MAX_ENTROPY_TOKENS, vp8_coef_encodings, vp8_coef_tree,
          coef_probs, branch_ct, coef_counts, 256, 1);
        fprintf(f, "%s\n      {", Comma(pt));

        t = 0;
        do {
          fprintf(f, "%s %d", Comma(t), coef_probs[t]);

        } while (++t < ENTROPY_NODES);

        fprintf(f, "}");
      } while (++pt < PREV_COEF_CONTEXTS);
      fprintf(f, "\n    }");
    } while (++band < COEF_BANDS);
    fprintf(f, "\n  }");
  } while (++type < BLOCK_TYPES);
  fprintf(f, "\n};\n");

  fprintf(f, "static const vp8_prob\n"
          "vp8_default_coef_probs_8x8[BLOCK_TYPES_8X8] [COEF_BANDS]\n"
          "[PREV_COEF_CONTEXTS] [ENTROPY_NODES] = {");
  type = 0;
  do {
    fprintf(f, "%s\n  { /* block Type %d */", Comma(type), type);
    band = 0;
    do {
      fprintf(f, "%s\n    { /* Coeff Band %d */", Comma(band), band);
      pt = 0;
      do {
        unsigned int branch_ct [ENTROPY_NODES] [2];
        unsigned int coef_counts[MAX_ENTROPY_TOKENS];
        vp8_prob coef_probs[ENTROPY_NODES];
        for (t = 0; t < MAX_ENTROPY_TOKENS; ++t)
          coef_counts[t] = context_counters_8x8[type] [band] [pt] [t];
        vp8_tree_probs_from_distribution(
          MAX_ENTROPY_TOKENS, vp8_coef_encodings, vp8_coef_tree,
          coef_probs, branch_ct, coef_counts, 256, 1);
        fprintf(f, "%s\n      {", Comma(pt));

        t = 0;
        do {
          fprintf(f, "%s %d", Comma(t), coef_probs[t]);
        } while (++t < ENTROPY_NODES);
        fprintf(f, "}");
      } while (++pt < PREV_COEF_CONTEXTS);
      fprintf(f, "\n    }");
    } while (++band < COEF_BANDS);
    fprintf(f, "\n  }");
  } while (++type < BLOCK_TYPES_8X8);
  fprintf(f, "\n};\n");

  fprintf(f, "static const vp8_prob\n"
          "vp8_default_coef_probs_16x16[BLOCK_TYPES_16X16] [COEF_BANDS]\n"
          "[PREV_COEF_CONTEXTS] [ENTROPY_NODES] = {");
  type = 0;
  do {
    fprintf(f, "%s\n  { /* block Type %d */", Comma(type), type);
    band = 0;
    do {
      fprintf(f, "%s\n    { /* Coeff Band %d */", Comma(band), band);
      pt = 0;
      do {
        unsigned int branch_ct [ENTROPY_NODES] [2];
        unsigned int coef_counts[MAX_ENTROPY_TOKENS];
        vp8_prob coef_probs[ENTROPY_NODES];
        for (t = 0; t < MAX_ENTROPY_TOKENS; ++t)
          coef_counts[t] = context_counters_16x16[type] [band] [pt] [t];
        vp8_tree_probs_from_distribution(
          MAX_ENTROPY_TOKENS, vp8_coef_encodings, vp8_coef_tree,
          coef_probs, branch_ct, coef_counts, 256, 1);
        fprintf(f, "%s\n      {", Comma(pt));

        t = 0;
        do {
          fprintf(f, "%s %d", Comma(t), coef_probs[t]);
        } while (++t < ENTROPY_NODES);
        fprintf(f, "}");
      } while (++pt < PREV_COEF_CONTEXTS);
      fprintf(f, "\n    }");
    } while (++band < COEF_BANDS);
    fprintf(f, "\n  }");
  } while (++type < BLOCK_TYPES_16X16);
  fprintf(f, "\n};\n");

  fclose(f);

  f = fopen("context.bin", "wb");
  fwrite(context_counters, sizeof(context_counters), 1, f);
  fwrite(context_counters_8x8, sizeof(context_counters_8x8), 1, f);
  fwrite(context_counters_16x16, sizeof(context_counters_16x16), 1, f);
  fclose(f);
}
#endif

void vp8_tokenize_initialize() {
  fill_value_tokens();
}

static __inline void stuff2nd_order_b_8x8(MACROBLOCKD *xd,
                                          const BLOCKD *const b,
                                          TOKENEXTRA **tp,
                                          ENTROPY_CONTEXT *a,
                                          ENTROPY_CONTEXT *l,
                                          VP8_COMP *cpi,
                                          int dry_run) {
  int pt; /* near block/prev token context index */
  TOKENEXTRA *t = *tp;        /* store tokens starting here */
  VP8_COMBINEENTROPYCONTEXTS(pt, *a, *l);
  (void) b;

  t->Token = DCT_EOB_TOKEN;
  t->context_tree = cpi->common.fc.coef_probs_8x8[PLANE_TYPE_Y2][0][pt];
  // t->section = 11;
  t->skip_eob_node = 0;
  ++t;

  *tp = t;
  if (!dry_run)
    ++cpi->coef_counts_8x8[PLANE_TYPE_Y2][0][pt][DCT_EOB_TOKEN];
  pt = 0;
  *a = *l = pt;
}

static __inline void stuff1st_order_b_8x8(MACROBLOCKD *xd,
                                          const BLOCKD *const b,
                                          TOKENEXTRA **tp,
                                          PLANE_TYPE type,
                                          ENTROPY_CONTEXT *a,
                                          ENTROPY_CONTEXT *l,
                                          VP8_COMP *cpi,
                                          int dry_run) {
  int pt; /* near block/prev token context index */
  TOKENEXTRA *t = *tp;        /* store tokens starting here */
  TX_TYPE tx_type = get_tx_type(xd, b);
  const int band = vp8_coef_bands_8x8[(type == PLANE_TYPE_Y_NO_DC) ? 1 : 0];
  VP8_COMBINEENTROPYCONTEXTS(pt, *a, *l);
  (void) b;

  t->Token = DCT_EOB_TOKEN;
  if (tx_type != DCT_DCT)
    t->context_tree = cpi->common.fc.hybrid_coef_probs_8x8[type][band][pt];
  else
    t->context_tree = cpi->common.fc.coef_probs_8x8[type][band][pt];
  // t->section = 8;
  t->skip_eob_node = 0;
  ++t;
  *tp = t;
  if (!dry_run) {
    if (tx_type == DCT_DCT)
      ++cpi->hybrid_coef_counts_8x8[type][band][pt][DCT_EOB_TOKEN];
    else
      ++cpi->coef_counts_8x8[type][band][pt][DCT_EOB_TOKEN];
  }
  pt = 0; /* 0 <-> all coeff data is zero */
  *a = *l = pt;
}

static __inline void stuff1st_order_buv_8x8(MACROBLOCKD *xd,
                                            const BLOCKD *const b,
                                            TOKENEXTRA **tp,
                                            ENTROPY_CONTEXT *a,
                                            ENTROPY_CONTEXT *l,
                                            VP8_COMP *cpi,
                                            int dry_run) {
  int pt; /* near block/prev token context index */
  TOKENEXTRA *t = *tp;        /* store tokens starting here */
  VP8_COMBINEENTROPYCONTEXTS(pt, *a, *l);
  (void) b;

  t->Token = DCT_EOB_TOKEN;
  t->context_tree = cpi->common.fc.coef_probs_8x8[PLANE_TYPE_UV][0][pt];
  // t->section = 13;
  t->skip_eob_node = 0;
  ++t;
  *tp = t;
  if (!dry_run)
    ++cpi->coef_counts_8x8[PLANE_TYPE_UV][0][pt][DCT_EOB_TOKEN];
  pt = 0; /* 0 <-> all coeff data is zero */
  *a = *l = pt;
}

static void vp8_stuff_mb_8x8(VP8_COMP *cpi, MACROBLOCKD *xd,
                             TOKENEXTRA **t, int dry_run) {
  ENTROPY_CONTEXT *A = (ENTROPY_CONTEXT *)xd->above_context;
  ENTROPY_CONTEXT *L = (ENTROPY_CONTEXT *)xd->left_context;
  PLANE_TYPE plane_type;
  int b;
  TOKENEXTRA *t_backup = *t;
  const int has_y2_block = (xd->mode_info_context->mbmi.mode != B_PRED
                            && xd->mode_info_context->mbmi.mode != I8X8_PRED
                            && xd->mode_info_context->mbmi.mode != SPLITMV);

  if (has_y2_block) {
    stuff2nd_order_b_8x8(xd, xd->block + 24, t,
                         A + vp8_block2above_8x8[24],
                         L + vp8_block2left_8x8[24], cpi, dry_run);
    plane_type = PLANE_TYPE_Y_NO_DC;
  } else {
    plane_type = PLANE_TYPE_Y_WITH_DC;
  }

  for (b = 0; b < 16; b += 4) {
    stuff1st_order_b_8x8(xd, xd->block + b, t, plane_type,
                         A + vp8_block2above_8x8[b],
                         L + vp8_block2left_8x8[b],
                         cpi, dry_run);
    *(A + vp8_block2above_8x8[b] + 1) = *(A + vp8_block2above_8x8[b]);
    *(L + vp8_block2left_8x8[b] + 1)  = *(L + vp8_block2left_8x8[b]);
  }

  for (b = 16; b < 24; b += 4) {
    stuff1st_order_buv_8x8(xd, xd->block + b, t,
                           A + vp8_block2above[b],
                           L + vp8_block2left[b],
                           cpi, dry_run);
    *(A + vp8_block2above_8x8[b] + 1) = *(A + vp8_block2above_8x8[b]);
    *(L + vp8_block2left_8x8[b] + 1) = *(L + vp8_block2left_8x8[b]);
  }
  if (dry_run)
    *t = t_backup;
}

static __inline void stuff1st_order_b_16x16(MACROBLOCKD *xd,
                                            const BLOCKD *const b,
                                            TOKENEXTRA **tp,
                                            PLANE_TYPE type,
                                            ENTROPY_CONTEXT *a,
                                            ENTROPY_CONTEXT *l,
                                            VP8_COMP *cpi,
                                            int dry_run) {
  int pt; /* near block/prev token context index */
  TOKENEXTRA *t = *tp;        /* store tokens starting here */
  TX_TYPE tx_type = get_tx_type(xd, b);
  const int band = vp8_coef_bands_16x16[(type == PLANE_TYPE_Y_NO_DC) ? 1 : 0];
  VP8_COMBINEENTROPYCONTEXTS(pt, *a, *l);
  (void) b;

  t->Token = DCT_EOB_TOKEN;
  if (tx_type != DCT_DCT)
    t->context_tree = cpi->common.fc.hybrid_coef_probs_16x16[type][band][pt];
  else
    t->context_tree = cpi->common.fc.coef_probs_16x16[type][band][pt];
  t->skip_eob_node = 0;
  ++t;
  *tp = t;
  if (!dry_run) {
    if (tx_type != DCT_DCT)
      ++cpi->hybrid_coef_counts_16x16[type][band][pt][DCT_EOB_TOKEN];
    else
      ++cpi->coef_counts_16x16[type][band][pt][DCT_EOB_TOKEN];
  }
  pt = 0; /* 0 <-> all coeff data is zero */
  *a = *l = pt;
}

static void vp8_stuff_mb_16x16(VP8_COMP *cpi, MACROBLOCKD *xd,
                               TOKENEXTRA **t, int dry_run) {
  ENTROPY_CONTEXT * A = (ENTROPY_CONTEXT *)xd->above_context;
  ENTROPY_CONTEXT * L = (ENTROPY_CONTEXT *)xd->left_context;
  int b, i;
  TOKENEXTRA *t_backup = *t;

  stuff1st_order_b_16x16(xd, xd->block, t, PLANE_TYPE_Y_WITH_DC,
                         A, L, cpi, dry_run);
  for (i = 1; i < 16; i++) {
    *(A + vp8_block2above[i]) = *(A);
    *(L +  vp8_block2left[i]) = *(L);
  }
  for (b = 16; b < 24; b += 4) {
    stuff1st_order_buv_8x8(xd, xd->block + b, t,
        A + vp8_block2above[b],
        L + vp8_block2left[b],
        cpi, dry_run);
    *(A + vp8_block2above_8x8[b]+1) = *(A + vp8_block2above_8x8[b]);
    *(L + vp8_block2left_8x8[b]+1 ) = *(L + vp8_block2left_8x8[b]);
  }
  vpx_memset(&A[8], 0, sizeof(A[8]));
  vpx_memset(&L[8], 0, sizeof(L[8]));
  if (dry_run)
    *t = t_backup;
}

static __inline void stuff2nd_order_b_4x4(MACROBLOCKD *xd,
                                          const BLOCKD *const b,
                                          TOKENEXTRA **tp,
                                          ENTROPY_CONTEXT *a,
                                          ENTROPY_CONTEXT *l,
                                          VP8_COMP *cpi,
                                          int dry_run) {
  int pt; /* near block/prev token context index */
  TOKENEXTRA *t = *tp;        /* store tokens starting here */
  VP8_COMBINEENTROPYCONTEXTS(pt, *a, *l);

  t->Token = DCT_EOB_TOKEN;
  t->context_tree = cpi->common.fc.coef_probs[PLANE_TYPE_Y2][0][pt];
  t->skip_eob_node = 0;
  ++t;
  *tp = t;
  if (!dry_run)
    ++cpi->coef_counts[PLANE_TYPE_Y2][0][pt] [DCT_EOB_TOKEN];

  pt = 0;
  *a = *l = pt;
}

static __inline void stuff1st_order_b_4x4(MACROBLOCKD *xd,
                                          const BLOCKD *const b,
                                          TOKENEXTRA **tp,
                                          PLANE_TYPE type,
                                          ENTROPY_CONTEXT *a,
                                          ENTROPY_CONTEXT *l,
                                          VP8_COMP *cpi,
                                          int dry_run) {
  int pt; /* near block/prev token context index */
  TOKENEXTRA *t = *tp;        /* store tokens starting here */
  TX_TYPE tx_type = get_tx_type(xd, b);
  const int band = vp8_coef_bands[(type == PLANE_TYPE_Y_NO_DC) ? 1 : 0];
  VP8_COMBINEENTROPYCONTEXTS(pt, *a, *l);

  t->Token = DCT_EOB_TOKEN;
  if (tx_type != DCT_DCT)
    t->context_tree = cpi->common.fc.hybrid_coef_probs[type][band][pt];
  else
    t->context_tree = cpi->common.fc.coef_probs[type][band][pt];
  t->skip_eob_node = 0;
  ++t;
  *tp = t;
  if (!dry_run) {
    if (tx_type != DCT_DCT)
      ++cpi->hybrid_coef_counts[type][band][pt][DCT_EOB_TOKEN];
    else
      ++cpi->coef_counts[type][band][pt][DCT_EOB_TOKEN];
  }
  pt = 0; /* 0 <-> all coeff data is zero */
  *a = *l = pt;
}

static __inline void stuff1st_order_buv_4x4(MACROBLOCKD *xd,
                                            const BLOCKD *const b,
                                            TOKENEXTRA **tp,
                                            ENTROPY_CONTEXT *a,
                                            ENTROPY_CONTEXT *l,
                                            VP8_COMP *cpi,
                                            int dry_run) {
  int pt; /* near block/prev token context index */
  TOKENEXTRA *t = *tp;        /* store tokens starting here */
  VP8_COMBINEENTROPYCONTEXTS(pt, *a, *l);

  t->Token = DCT_EOB_TOKEN;
  t->context_tree = cpi->common.fc.coef_probs[PLANE_TYPE_UV][0][pt];
  t->skip_eob_node = 0;
  ++t;
  *tp = t;
  if (!dry_run)
    ++cpi->coef_counts[PLANE_TYPE_UV][0][pt][DCT_EOB_TOKEN];
  pt = 0; /* 0 <-> all coeff data is zero */
  *a = *l = pt;
}

static void vp8_stuff_mb_4x4(VP8_COMP *cpi, MACROBLOCKD *xd,
                             TOKENEXTRA **t, int dry_run) {
  ENTROPY_CONTEXT *A = (ENTROPY_CONTEXT *)xd->above_context;
  ENTROPY_CONTEXT *L = (ENTROPY_CONTEXT *)xd->left_context;
  int b;
  TOKENEXTRA *t_backup = *t;
  PLANE_TYPE plane_type;
  const int has_y2_block = (xd->mode_info_context->mbmi.mode != B_PRED
                            && xd->mode_info_context->mbmi.mode != I8X8_PRED
                            && xd->mode_info_context->mbmi.mode != SPLITMV);

  if (has_y2_block) {
    stuff2nd_order_b_4x4(xd, xd->block + 24, t,
                         A + vp8_block2above[24],
                         L + vp8_block2left[24],
                         cpi, dry_run);
    plane_type = PLANE_TYPE_Y_NO_DC;
  } else {
    plane_type = PLANE_TYPE_Y_WITH_DC;
  }

  for (b = 0; b < 16; b++)
    stuff1st_order_b_4x4(xd, xd->block + b, t, plane_type,
                         A + vp8_block2above[b],
                         L + vp8_block2left[b],
                         cpi, dry_run);

  for (b = 16; b < 24; b++)
    stuff1st_order_buv_4x4(xd, xd->block + b, t,
                           A + vp8_block2above[b],
                           L + vp8_block2left[b],
                           cpi, dry_run);

  if (dry_run)
    *t = t_backup;
}

static void vp8_stuff_mb_8x8_4x4uv(VP8_COMP *cpi, MACROBLOCKD *xd,
                                   TOKENEXTRA **t, int dry_run) {
  ENTROPY_CONTEXT *A = (ENTROPY_CONTEXT *)xd->above_context;
  ENTROPY_CONTEXT *L = (ENTROPY_CONTEXT *)xd->left_context;
  int b;
  TOKENEXTRA *t_backup = *t;

  for (b = 0; b < 16; b += 4) {
    stuff1st_order_b_8x8(xd, xd->block + b, t, PLANE_TYPE_Y_WITH_DC,
                         A + vp8_block2above_8x8[b],
                         L + vp8_block2left_8x8[b],
                         cpi, dry_run);
    *(A + vp8_block2above_8x8[b] + 1) = *(A + vp8_block2above_8x8[b]);
    *(L + vp8_block2left_8x8[b] + 1)  = *(L + vp8_block2left_8x8[b]);
  }

  for (b = 16; b < 24; b++)
    stuff1st_order_buv_4x4(xd, xd->block + b, t,
                           A + vp8_block2above[b],
                           L + vp8_block2left[b],
                           cpi, dry_run);

  if (dry_run)
    *t = t_backup;
}

void vp8_stuff_mb(VP8_COMP *cpi, MACROBLOCKD *xd, TOKENEXTRA **t, int dry_run) {
  TX_SIZE tx_size = xd->mode_info_context->mbmi.txfm_size;

  if (tx_size == TX_16X16) {
    vp8_stuff_mb_16x16(cpi, xd, t, dry_run);
  } else if (tx_size == TX_8X8) {
    if (xd->mode_info_context->mbmi.mode == I8X8_PRED) {
      vp8_stuff_mb_8x8_4x4uv(cpi, xd, t, dry_run);
    } else {
      vp8_stuff_mb_8x8(cpi, xd, t, dry_run);
    }
  } else {
    vp8_stuff_mb_4x4(cpi, xd, t, dry_run);
  }
}

void vp8_fix_contexts(MACROBLOCKD *xd) {
  /* Clear entropy contexts for Y2 blocks */
  if ((xd->mode_info_context->mbmi.mode != B_PRED
      && xd->mode_info_context->mbmi.mode != I8X8_PRED
      && xd->mode_info_context->mbmi.mode != SPLITMV)
      || xd->mode_info_context->mbmi.txfm_size == TX_16X16
      ) {
    vpx_memset(xd->above_context, 0, sizeof(ENTROPY_CONTEXT_PLANES));
    vpx_memset(xd->left_context, 0, sizeof(ENTROPY_CONTEXT_PLANES));
  } else {
    vpx_memset(xd->above_context, 0, sizeof(ENTROPY_CONTEXT_PLANES) - 1);
    vpx_memset(xd->left_context, 0, sizeof(ENTROPY_CONTEXT_PLANES) - 1);
  }
}
