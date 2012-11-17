/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */


#include "vp9/common/type_aliases.h"
#include "vp9/common/blockd.h"
#include "onyxd_int.h"
#include "vpx_mem/vpx_mem.h"
#include "vpx_ports/mem.h"
#include "detokenize.h"

#include "vp9/common/seg_common.h"

#define BOOL_DATA UINT8

#define OCB_X PREV_COEF_CONTEXTS * ENTROPY_NODES

DECLARE_ALIGNED(16, static const int, coef_bands_x[16]) = {
  0 * OCB_X, 1 * OCB_X, 2 * OCB_X, 3 * OCB_X,
  6 * OCB_X, 4 * OCB_X, 5 * OCB_X, 6 * OCB_X,
  6 * OCB_X, 6 * OCB_X, 6 * OCB_X, 6 * OCB_X,
  6 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X
};
DECLARE_ALIGNED(16, static const int, coef_bands_x_8x8[64]) = {
  0 * OCB_X, 1 * OCB_X, 2 * OCB_X, 3 * OCB_X, 5 * OCB_X, 4 * OCB_X, 4 * OCB_X, 5 * OCB_X,
  5 * OCB_X, 3 * OCB_X, 6 * OCB_X, 3 * OCB_X, 5 * OCB_X, 4 * OCB_X, 6 * OCB_X, 6 * OCB_X,
  6 * OCB_X, 5 * OCB_X, 5 * OCB_X, 6 * OCB_X, 6 * OCB_X, 6 * OCB_X, 6 * OCB_X, 6 * OCB_X,
  6 * OCB_X, 6 * OCB_X, 6 * OCB_X, 6 * OCB_X, 6 * OCB_X, 6 * OCB_X, 6 * OCB_X, 6 * OCB_X,
  6 * OCB_X, 6 * OCB_X, 6 * OCB_X, 6 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X,
  7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X,
  7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X,
  7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X,
};

DECLARE_ALIGNED(16, static const int, coef_bands_x_16x16[256]) = {
  0 * OCB_X, 1 * OCB_X, 2 * OCB_X, 3 * OCB_X, 5 * OCB_X, 4 * OCB_X, 4 * OCB_X, 5 * OCB_X, 5 * OCB_X, 3 * OCB_X, 6 * OCB_X, 3 * OCB_X, 5 * OCB_X, 4 * OCB_X, 6 * OCB_X, 6 * OCB_X,
  6 * OCB_X, 5 * OCB_X, 5 * OCB_X, 6 * OCB_X, 6 * OCB_X, 6 * OCB_X, 6 * OCB_X, 6 * OCB_X, 6 * OCB_X, 6 * OCB_X, 6 * OCB_X, 6 * OCB_X, 6 * OCB_X, 6 * OCB_X, 6 * OCB_X, 6 * OCB_X,
  6 * OCB_X, 6 * OCB_X, 6 * OCB_X, 6 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X,
  7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X,
  7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X,
  7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X,
  7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X,
  7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X,
  7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X,
  7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X,
  7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X,
  7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X,
  7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X,
  7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X,
  7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X,
  7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X
};

#define EOB_CONTEXT_NODE            0
#define ZERO_CONTEXT_NODE           1
#define ONE_CONTEXT_NODE            2
#define LOW_VAL_CONTEXT_NODE        3
#define TWO_CONTEXT_NODE            4
#define THREE_CONTEXT_NODE          5
#define HIGH_LOW_CONTEXT_NODE       6
#define CAT_ONE_CONTEXT_NODE        7
#define CAT_THREEFOUR_CONTEXT_NODE  8
#define CAT_THREE_CONTEXT_NODE      9
#define CAT_FIVE_CONTEXT_NODE       10

#define CAT1_MIN_VAL    5
#define CAT2_MIN_VAL    7
#define CAT3_MIN_VAL   11
#define CAT4_MIN_VAL   19
#define CAT5_MIN_VAL   35
#define CAT6_MIN_VAL   67
#define CAT1_PROB0    159
#define CAT2_PROB0    145
#define CAT2_PROB1    165

#define CAT3_PROB0 140
#define CAT3_PROB1 148
#define CAT3_PROB2 173

#define CAT4_PROB0 135
#define CAT4_PROB1 140
#define CAT4_PROB2 155
#define CAT4_PROB3 176

#define CAT5_PROB0 130
#define CAT5_PROB1 134
#define CAT5_PROB2 141
#define CAT5_PROB3 157
#define CAT5_PROB4 180

static const unsigned char cat6_prob[14] =
{ 254, 254, 252, 249, 243, 230, 196, 177, 153, 140, 133, 130, 129, 0 };

void vp9_reset_mb_tokens_context(MACROBLOCKD* const xd) {
  /* Clear entropy contexts for Y2 blocks */
  if ((xd->mode_info_context->mbmi.mode != B_PRED &&
      xd->mode_info_context->mbmi.mode != I8X8_PRED &&
      xd->mode_info_context->mbmi.mode != SPLITMV)
      || xd->mode_info_context->mbmi.txfm_size == TX_16X16
      ) {
    vpx_memset(xd->above_context, 0, sizeof(ENTROPY_CONTEXT_PLANES));
    vpx_memset(xd->left_context, 0, sizeof(ENTROPY_CONTEXT_PLANES));
  } else {
    vpx_memset(xd->above_context, 0, sizeof(ENTROPY_CONTEXT_PLANES) - 1);
    vpx_memset(xd->left_context, 0, sizeof(ENTROPY_CONTEXT_PLANES) - 1);
  }
}

DECLARE_ALIGNED(16, extern const unsigned char, vp9_norm[256]);

// #define PREV_CONTEXT_INC(val) (2+((val)>2))
// #define PREV_CONTEXT_INC(val) (vp9_prev_token_class[(val)])
#define PREV_CONTEXT_INC(val) (vp9_prev_token_class[(val)>10?10:(val)])

static int get_signed(BOOL_DECODER *br, int value_to_sign) {
  const int split = (br->range + 1) >> 1;
  const VP9_BD_VALUE bigsplit = (VP9_BD_VALUE)split << (VP9_BD_VALUE_SIZE - 8);
  int v;

  if (br->count < 0)
    vp9_bool_decoder_fill(br);

  if (br->value < bigsplit) {
    br->range = split;
    v = value_to_sign;
  } else {
    br->range = br->range - split;
    br->value = br->value - bigsplit;
    v = -value_to_sign;
  }
  br->range += br->range;
  br->value += br->value;
  --br->count;

  return v;
}

#define INCREMENT_COUNT(token)               \
  do {                                       \
    coef_counts[coef_bands[c]][pt][token]++; \
    pt = vp9_prev_token_class[token];        \
  } while (0)

#define WRITE_COEF_CONTINUE(val, token)                       \
  {                                                           \
    prob = coef_probs + (ENTROPY_NODES*PREV_CONTEXT_INC(val));\
    qcoeff_ptr[scan[c]] = (INT16) get_signed(br, val);        \
    INCREMENT_COUNT(token);                                   \
    c++;                                                      \
    continue;                                                 \
  }

#define ADJUST_COEF(prob, bits_count)  \
  do {                                 \
    if (vp9_read(br, prob))            \
      val += (UINT16)(1 << bits_count);\
  } while (0);

static int decode_coefs(VP9D_COMP *dx, const MACROBLOCKD *xd,
                        BOOL_DECODER* const br,
                        ENTROPY_CONTEXT *a, ENTROPY_CONTEXT *l,
                        PLANE_TYPE type,
                        TX_TYPE tx_type,
                        int seg_eob, INT16 *qcoeff_ptr, int i,
                        const int *const scan, int block_type,
                        const int *coef_bands_x,
                        const int *coef_bands) {
  FRAME_CONTEXT *const fc = &dx->common.fc;
  int pt, c = (type == PLANE_TYPE_Y_NO_DC);
  const vp9_prob *prob, *coef_probs;
  unsigned int (*coef_counts)[PREV_COEF_CONTEXTS][MAX_ENTROPY_TOKENS];

  switch (block_type) {
    default:
    case TX_4X4:
      if (tx_type == DCT_DCT) {
        coef_probs  = fc->coef_probs[type][0][0];
        coef_counts = fc->coef_counts[type];
      } else {
        coef_probs  = fc->hybrid_coef_probs[type][0][0];
        coef_counts = fc->hybrid_coef_counts[type];
      }
      break;
    case TX_8X8:
      if (tx_type == DCT_DCT) {
        coef_probs  = fc->coef_probs_8x8[type][0][0];
        coef_counts = fc->coef_counts_8x8[type];
      } else {
        coef_probs  = fc->hybrid_coef_probs_8x8[type][0][0];
        coef_counts = fc->hybrid_coef_counts_8x8[type];
      }
      break;
    case TX_16X16:
      if (tx_type == DCT_DCT) {
        coef_probs  = fc->coef_probs_16x16[type][0][0];
        coef_counts = fc->coef_counts_16x16[type];
      } else {
        coef_probs  = fc->hybrid_coef_probs_16x16[type][0][0];
        coef_counts = fc->hybrid_coef_counts_16x16[type];
      }
      break;
  }

  VP9_COMBINEENTROPYCONTEXTS(pt, *a, *l);
  prob = coef_probs + pt * ENTROPY_NODES;

  while (1) {
    int val;
    const uint8_t *cat6 = cat6_prob;
    if (c >= seg_eob) break;
    prob += coef_bands_x[c];
    if (!vp9_read(br, prob[EOB_CONTEXT_NODE]))
      break;
SKIP_START:
    if (c >= seg_eob) break;
    if (!vp9_read(br, prob[ZERO_CONTEXT_NODE])) {
      INCREMENT_COUNT(ZERO_TOKEN);
      ++c;
      prob = coef_probs + coef_bands_x[c];
      goto SKIP_START;
    }
    // ONE_CONTEXT_NODE_0_
    if (!vp9_read(br, prob[ONE_CONTEXT_NODE])) {
      prob = coef_probs + ENTROPY_NODES;
      qcoeff_ptr[scan[c]] = (INT16) get_signed(br, 1);
      INCREMENT_COUNT(ONE_TOKEN);
      ++c;
      continue;
    }
    // LOW_VAL_CONTEXT_NODE_0_
    if (!vp9_read(br, prob[LOW_VAL_CONTEXT_NODE])) {
      if (!vp9_read(br, prob[TWO_CONTEXT_NODE])) {
        WRITE_COEF_CONTINUE(2, TWO_TOKEN);
      }
      if (!vp9_read(br, prob[THREE_CONTEXT_NODE])) {
        WRITE_COEF_CONTINUE(3, THREE_TOKEN);
      }
      WRITE_COEF_CONTINUE(4, FOUR_TOKEN);
    }
    // HIGH_LOW_CONTEXT_NODE_0_
    if (!vp9_read(br, prob[HIGH_LOW_CONTEXT_NODE])) {
      if (!vp9_read(br, prob[CAT_ONE_CONTEXT_NODE])) {
        val = CAT1_MIN_VAL;
        ADJUST_COEF(CAT1_PROB0, 0);
        WRITE_COEF_CONTINUE(val, DCT_VAL_CATEGORY1);
      }
      val = CAT2_MIN_VAL;
      ADJUST_COEF(CAT2_PROB1, 1);
      ADJUST_COEF(CAT2_PROB0, 0);
      WRITE_COEF_CONTINUE(val, DCT_VAL_CATEGORY2);
    }
    // CAT_THREEFOUR_CONTEXT_NODE_0_
    if (!vp9_read(br, prob[CAT_THREEFOUR_CONTEXT_NODE])) {
      if (!vp9_read(br, prob[CAT_THREE_CONTEXT_NODE])) {
        val = CAT3_MIN_VAL;
        ADJUST_COEF(CAT3_PROB2, 2);
        ADJUST_COEF(CAT3_PROB1, 1);
        ADJUST_COEF(CAT3_PROB0, 0);
        WRITE_COEF_CONTINUE(val, DCT_VAL_CATEGORY3);
      }
      val = CAT4_MIN_VAL;
      ADJUST_COEF(CAT4_PROB3, 3);
      ADJUST_COEF(CAT4_PROB2, 2);
      ADJUST_COEF(CAT4_PROB1, 1);
      ADJUST_COEF(CAT4_PROB0, 0);
      WRITE_COEF_CONTINUE(val, DCT_VAL_CATEGORY4);
    }
    // CAT_FIVE_CONTEXT_NODE_0_:
    if (!vp9_read(br, prob[CAT_FIVE_CONTEXT_NODE])) {
      val = CAT5_MIN_VAL;
      ADJUST_COEF(CAT5_PROB4, 4);
      ADJUST_COEF(CAT5_PROB3, 3);
      ADJUST_COEF(CAT5_PROB2, 2);
      ADJUST_COEF(CAT5_PROB1, 1);
      ADJUST_COEF(CAT5_PROB0, 0);
      WRITE_COEF_CONTINUE(val, DCT_VAL_CATEGORY5);
    }
    val = 0;
    while (*cat6) {
      val = (val << 1) | vp9_read(br, *cat6++);
    }
    val += CAT6_MIN_VAL;
    WRITE_COEF_CONTINUE(val, DCT_VAL_CATEGORY6);
  }

  if (c < seg_eob)
    coef_counts[coef_bands[c]][pt][DCT_EOB_TOKEN]++;

  return c;
}


int get_eob(MACROBLOCKD* const xd, int segment_id, int eob_max) {
  int active = vp9_segfeature_active(xd, segment_id, SEG_LVL_EOB);
  int eob = vp9_get_segdata(xd, segment_id, SEG_LVL_EOB);

  if (!active || eob > eob_max)
    eob = eob_max;
  return eob;
}


int vp9_decode_mb_tokens_16x16(VP9D_COMP* const pbi,
                               MACROBLOCKD* const xd,
                               BOOL_DECODER* const bc) {
  ENTROPY_CONTEXT* const A = (ENTROPY_CONTEXT *)xd->above_context;
  ENTROPY_CONTEXT* const L = (ENTROPY_CONTEXT *)xd->left_context;

  unsigned short* const eobs = xd->eobs;
  PLANE_TYPE type;
  int c, i, eobtotal = 0, seg_eob;
  const int segment_id = xd->mode_info_context->mbmi.segment_id;
  INT16 *qcoeff_ptr = &xd->qcoeff[0];
  TX_TYPE tx_type = get_tx_type(xd, &xd->block[0]);

  type = PLANE_TYPE_Y_WITH_DC;
  seg_eob = get_eob(xd, segment_id, 256);

  // Luma block
  {
    const int* const scan = vp9_default_zig_zag1d_16x16;
    c = decode_coefs(pbi, xd, bc, A, L, type,
                     tx_type,
                     seg_eob, qcoeff_ptr,
                     0, scan, TX_16X16, coef_bands_x_16x16,
                     vp9_coef_bands_16x16);
    eobs[0] = c;
    A[0] = L[0] = (c != !type);
    A[1] = A[2] = A[3] = A[0];
    L[1] = L[2] = L[3] = L[0];
    eobtotal += c;
  }

  // 8x8 chroma blocks
  qcoeff_ptr += 256;
  type = PLANE_TYPE_UV;
  tx_type = DCT_DCT;
  seg_eob = get_eob(xd, segment_id, 64);
  for (i = 16; i < 24; i += 4) {
    ENTROPY_CONTEXT* const a = A + vp9_block2above_8x8[i];
    ENTROPY_CONTEXT* const l = L + vp9_block2left_8x8[i];
    const int* const scan = vp9_default_zig_zag1d_8x8;

    c = decode_coefs(pbi, xd, bc, a, l, type,
                     tx_type,
                     seg_eob, qcoeff_ptr,
                     i, scan, TX_8X8, coef_bands_x_8x8,
                     vp9_coef_bands_8x8);
    a[0] = l[0] = ((eobs[i] = c) != !type);
    a[1] = a[0];
    l[1] = l[0];

    eobtotal += c;
    qcoeff_ptr += 64;
  }
  vpx_memset(&A[8], 0, sizeof(A[8]));
  vpx_memset(&L[8], 0, sizeof(L[8]));
  return eobtotal;
}

int vp9_decode_mb_tokens_8x8(VP9D_COMP* const pbi,
                             MACROBLOCKD* const xd,
                             BOOL_DECODER* const bc) {
  ENTROPY_CONTEXT *const A = (ENTROPY_CONTEXT *)xd->above_context;
  ENTROPY_CONTEXT *const L = (ENTROPY_CONTEXT *)xd->left_context;

  unsigned short *const eobs = xd->eobs;
  PLANE_TYPE type;
  int c, i, eobtotal = 0, seg_eob;
  const int segment_id = xd->mode_info_context->mbmi.segment_id;
  INT16 *qcoeff_ptr = &xd->qcoeff[0];
  TX_TYPE tx_type = DCT_DCT;

  int bufthred = (xd->mode_info_context->mbmi.mode == I8X8_PRED ||
                  xd->mode_info_context->mbmi.mode == SPLITMV) ? 16 : 24;
  if (xd->mode_info_context->mbmi.mode != B_PRED &&
      xd->mode_info_context->mbmi.mode != SPLITMV &&
      xd->mode_info_context->mbmi.mode != I8X8_PRED) {
    ENTROPY_CONTEXT *const a = A + vp9_block2above_8x8[24];
    ENTROPY_CONTEXT *const l = L + vp9_block2left_8x8[24];
    const int *const scan = vp9_default_zig_zag1d;
    type = PLANE_TYPE_Y2;

    seg_eob = get_eob(xd, segment_id, 4);
    c = decode_coefs(pbi, xd, bc, a, l, type,
                     tx_type,
                     seg_eob, qcoeff_ptr + 24 * 16,
                     24, scan, TX_8X8, coef_bands_x,
                     vp9_coef_bands);
    a[0] = l[0] = ((eobs[24] = c) != !type);

    eobtotal += c - 4;

    type = PLANE_TYPE_Y_NO_DC;
  } else
    type = PLANE_TYPE_Y_WITH_DC;

  seg_eob = get_eob(xd, segment_id, 64);

  for (i = 0; i < bufthred ; i += 4) {
    ENTROPY_CONTEXT *const a = A + vp9_block2above_8x8[i];
    ENTROPY_CONTEXT *const l = L + vp9_block2left_8x8[i];
    const int *const scan = vp9_default_zig_zag1d_8x8;
    tx_type = DCT_DCT;

    if (i == 16)
      type = PLANE_TYPE_UV;
    if (type == PLANE_TYPE_Y_WITH_DC) {
      tx_type = get_tx_type(xd, xd->block + i);
    }

    c = decode_coefs(pbi, xd, bc, a, l, type,
                     tx_type,
                     seg_eob, qcoeff_ptr,
                     i, scan, TX_8X8, coef_bands_x_8x8,
                     vp9_coef_bands_8x8);
    a[0] = l[0] = ((eobs[i] = c) != !type);
    a[1] = a[0];
    l[1] = l[0];

    eobtotal += c;
    qcoeff_ptr += 64;
  }

  if (bufthred == 16) {
    type = PLANE_TYPE_UV;
    tx_type = DCT_DCT;
    seg_eob = get_eob(xd, segment_id, 16);

    // use 4x4 transform for U, V components in I8X8 prediction mode
    for (i = 16; i < 24; i++) {
      ENTROPY_CONTEXT *const a = A + vp9_block2above[i];
      ENTROPY_CONTEXT *const l = L + vp9_block2left[i];
      const int *scan = vp9_default_zig_zag1d;

      c = decode_coefs(pbi, xd, bc, a, l, type,
                       tx_type,
                       seg_eob, qcoeff_ptr,
                       i, scan, TX_4X4, coef_bands_x,
                       vp9_coef_bands);
      a[0] = l[0] = ((eobs[i] = c) != !type);

      eobtotal += c;
      qcoeff_ptr += 16;
    }
  }

  return eobtotal;
}

static int decode_coefs_4x4(VP9D_COMP *dx, MACROBLOCKD *xd,
                            BOOL_DECODER* const bc,
                            PLANE_TYPE type, int i) {
  ENTROPY_CONTEXT *const A = (ENTROPY_CONTEXT *)xd->above_context;
  ENTROPY_CONTEXT *const L = (ENTROPY_CONTEXT *)xd->left_context;
  ENTROPY_CONTEXT *const a = A + vp9_block2above[i];
  ENTROPY_CONTEXT *const l = L + vp9_block2left[i];
  INT16 *qcoeff_ptr = &xd->qcoeff[0];
  const int *scan = vp9_default_zig_zag1d;
  unsigned short *const eobs = xd->eobs;
  int segment_id = xd->mode_info_context->mbmi.segment_id;
  int c, seg_eob = get_eob(xd, segment_id, 16);
  TX_TYPE tx_type = DCT_DCT;

  if (type == PLANE_TYPE_Y_WITH_DC)
    tx_type = get_tx_type(xd, &xd->block[i]);
  switch (tx_type) {
    case ADST_DCT :
      scan = vp9_row_scan;
      break;

    case DCT_ADST :
      scan = vp9_col_scan;
      break;

    default :
      scan = vp9_default_zig_zag1d;
      break;
  }
  c = decode_coefs(dx, xd, bc, a, l, type,
                   tx_type,
                   seg_eob, qcoeff_ptr + i * 16,
                   i, scan, TX_4X4, coef_bands_x, vp9_coef_bands);
  a[0] = l[0] = ((eobs[i] = c) != !type);
  return c;
}

int vp9_decode_mb_tokens_4x4(VP9D_COMP* const dx,
                             MACROBLOCKD* const xd,
                             BOOL_DECODER* const bc) {
  int i, eobtotal = 0;
  PLANE_TYPE type;

  if (xd->mode_info_context->mbmi.mode != B_PRED &&
      xd->mode_info_context->mbmi.mode != I8X8_PRED &&
      xd->mode_info_context->mbmi.mode != SPLITMV) {
    eobtotal += decode_coefs_4x4(dx, xd, bc, PLANE_TYPE_Y2, 24) - 16;
    type = PLANE_TYPE_Y_NO_DC;
  } else {
    type = PLANE_TYPE_Y_WITH_DC;
  }

  for (i = 0; i < 16; ++i) {
    eobtotal += decode_coefs_4x4(dx, xd, bc, type, i);
  }
  do {
    eobtotal += decode_coefs_4x4(dx, xd, bc, PLANE_TYPE_UV, i);
  } while (++i < 24);
  return eobtotal;
}
