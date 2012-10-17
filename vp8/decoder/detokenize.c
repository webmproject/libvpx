/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */


#include "vp8/common/type_aliases.h"
#include "vp8/common/blockd.h"
#include "onyxd_int.h"
#include "vpx_mem/vpx_mem.h"
#include "vpx_ports/mem.h"
#include "detokenize.h"

#include "vp8/common/seg_common.h"

#define BOOL_DATA UINT8

#define OCB_X PREV_COEF_CONTEXTS * ENTROPY_NODES

DECLARE_ALIGNED(16, const int, coef_bands_x[16]) = {
  0 * OCB_X, 1 * OCB_X, 2 * OCB_X, 3 * OCB_X,
  6 * OCB_X, 4 * OCB_X, 5 * OCB_X, 6 * OCB_X,
  6 * OCB_X, 6 * OCB_X, 6 * OCB_X, 6 * OCB_X,
  6 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X
};
DECLARE_ALIGNED(16, const int, coef_bands_x_8x8[64]) = {
  0 * OCB_X, 1 * OCB_X, 2 * OCB_X, 3 * OCB_X, 5 * OCB_X, 4 * OCB_X, 4 * OCB_X, 5 * OCB_X,
  5 * OCB_X, 3 * OCB_X, 6 * OCB_X, 3 * OCB_X, 5 * OCB_X, 4 * OCB_X, 6 * OCB_X, 6 * OCB_X,
  6 * OCB_X, 5 * OCB_X, 5 * OCB_X, 6 * OCB_X, 6 * OCB_X, 6 * OCB_X, 6 * OCB_X, 6 * OCB_X,
  6 * OCB_X, 6 * OCB_X, 6 * OCB_X, 6 * OCB_X, 6 * OCB_X, 6 * OCB_X, 6 * OCB_X, 6 * OCB_X,
  6 * OCB_X, 6 * OCB_X, 6 * OCB_X, 6 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X,
  7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X,
  7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X,
  7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X, 7 * OCB_X,
};

DECLARE_ALIGNED(16, const int, coef_bands_x_16x16[256]) = {
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

void vp8_reset_mb_tokens_context(MACROBLOCKD *xd) {
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

DECLARE_ALIGNED(16, extern const unsigned char, vp8_norm[256]);

// #define PREV_CONTEXT_INC(val) (2+((val)>2))
// #define PREV_CONTEXT_INC(val) (vp8_prev_token_class[(val)])
#define PREV_CONTEXT_INC(val) (vp8_prev_token_class[(val)>10?10:(val)])


int get_token(int v) {
  if (v < 0) v = -v;
  if (v == 0) return ZERO_TOKEN;
  else if (v == 1) return ONE_TOKEN;
  else if (v == 2) return TWO_TOKEN;
  else if (v == 3) return THREE_TOKEN;
  else if (v == 4) return FOUR_TOKEN;
  else if (v <= 6) return DCT_VAL_CATEGORY1;
  else if (v <= 10) return DCT_VAL_CATEGORY2;
  else if (v <= 18) return DCT_VAL_CATEGORY3;
  else if (v <= 34) return DCT_VAL_CATEGORY4;
  else if (v <= 66) return DCT_VAL_CATEGORY5;
  else return DCT_VAL_CATEGORY6;
}

#if CONFIG_HYBRIDTRANSFORM
void static count_tokens_adaptive_scan(const MACROBLOCKD *xd, INT16 *qcoeff_ptr,
                                       int block, PLANE_TYPE type,
                                       TX_TYPE tx_type,
                                       ENTROPY_CONTEXT *a, ENTROPY_CONTEXT *l,
                                       int eob, int seg_eob,
                                       FRAME_CONTEXT *fc) {
  int c, pt, token, band;
  const int *scan;

  switch(tx_type) {
    case ADST_DCT :
      scan = vp8_row_scan;
      break;

    case DCT_ADST :
      scan = vp8_col_scan;
      break;

    default :
      scan = vp8_default_zig_zag1d;
      break;
  }

  VP8_COMBINEENTROPYCONTEXTS(pt, *a, *l);
  for (c = !type; c < eob; ++c) {
    int rc = scan[c];
    int v = qcoeff_ptr[rc];
    band = vp8_coef_bands[c];
    token = get_token(v);
    if (tx_type != DCT_DCT)
      fc->hybrid_coef_counts[type][band][pt][token]++;
    else
      fc->coef_counts[type][band][pt][token]++;
    pt = vp8_prev_token_class[token];
  }

  if (eob < seg_eob) {
    band = vp8_coef_bands[c];
    if (tx_type != DCT_DCT)
      fc->hybrid_coef_counts[type][band][pt][DCT_EOB_TOKEN]++;
    else
      fc->coef_counts[type][band][pt][DCT_EOB_TOKEN]++;
  }
}
#endif

void static count_tokens(INT16 *qcoeff_ptr, int block, PLANE_TYPE type,
                         ENTROPY_CONTEXT *a, ENTROPY_CONTEXT *l,
                         int eob, int seg_eob, FRAME_CONTEXT *const fc) {
  int c, pt, token, band;
  VP8_COMBINEENTROPYCONTEXTS(pt, *a, *l);
  for (c = !type; c < eob; ++c) {
    int rc = vp8_default_zig_zag1d[c];
    int v = qcoeff_ptr[rc];
    band = vp8_coef_bands[c];
    token = get_token(v);
    fc->coef_counts[type][band][pt][token]++;
    pt = vp8_prev_token_class[token];
  }
  if (eob < seg_eob) {
    band = vp8_coef_bands[c];
    fc->coef_counts[type][band][pt][DCT_EOB_TOKEN]++;
  }
}

void static count_tokens_8x8(INT16 *qcoeff_ptr, int block, PLANE_TYPE type,
#if CONFIG_HYBRIDTRANSFORM8X8
                             TX_TYPE tx_type,
#endif
                             ENTROPY_CONTEXT *a, ENTROPY_CONTEXT *l,
                             int eob, int seg_eob, FRAME_CONTEXT *fc) {
  int c, pt, token, band;
  VP8_COMBINEENTROPYCONTEXTS(pt, *a, *l);
  for (c = !type; c < eob; ++c) {
    int rc = (type == 1 ? vp8_default_zig_zag1d[c] : vp8_default_zig_zag1d_8x8[c]);
    int v = qcoeff_ptr[rc];
    band = (type == 1 ? vp8_coef_bands[c] : vp8_coef_bands_8x8[c]);
    token = get_token(v);
#if CONFIG_HYBRIDTRANSFORM8X8
    if (tx_type != DCT_DCT)
      fc->hybrid_coef_counts_8x8[type][band][pt][token]++;
    else
#endif
      fc->coef_counts_8x8[type][band][pt][token]++;
    pt = vp8_prev_token_class[token];
  }
  if (eob < seg_eob) {
    band = (type == 1 ? vp8_coef_bands[c] : vp8_coef_bands_8x8[c]);
#if CONFIG_HYBRIDTRANSFORM8X8
    if (tx_type != DCT_DCT)
      fc->hybrid_coef_counts_8x8[type][band][pt][DCT_EOB_TOKEN]++;
    else
#endif
      fc->coef_counts_8x8[type][band][pt][DCT_EOB_TOKEN]++;
  }
}

void static count_tokens_16x16(INT16 *qcoeff_ptr, int block, PLANE_TYPE type,
#if CONFIG_HYBRIDTRANSFORM16X16
                               TX_TYPE tx_type,
#endif
                               ENTROPY_CONTEXT *a, ENTROPY_CONTEXT *l,
                               int eob, int seg_eob, FRAME_CONTEXT *fc) {
  int c, pt, token;
  VP8_COMBINEENTROPYCONTEXTS(pt, *a, *l);
  for (c = !type; c < eob; ++c) {
    int rc = vp8_default_zig_zag1d_16x16[c];
    int v = qcoeff_ptr[rc];
    int band = vp8_coef_bands_16x16[c];
    token = get_token(v);
#if CONFIG_HYBRIDTRANSFORM16X16
    if (tx_type != DCT_DCT)
      fc->hybrid_coef_counts_16x16[type][band][pt][token]++;
    else
#endif
      fc->coef_counts_16x16[type][band][pt][token]++;
    pt = vp8_prev_token_class[token];
  }
  if (eob < seg_eob) {
    int band = vp8_coef_bands_16x16[c];
#if CONFIG_HYBRIDTRANSFORM16X16
    if (tx_type != DCT_DCT)
      fc->hybrid_coef_counts_16x16[type][band][pt][DCT_EOB_TOKEN]++;
    else
#endif
      fc->coef_counts_16x16[type][band][pt][DCT_EOB_TOKEN]++;
  }
}

static int vp8_get_signed(BOOL_DECODER *br, int value_to_sign) {
  const int split = (br->range + 1) >> 1;
  const VP8_BD_VALUE bigsplit = (VP8_BD_VALUE)split << (VP8_BD_VALUE_SIZE - 8);
  int v;

  if (br->count < 0)
    vp8dx_bool_decoder_fill(br);

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

#define WRITE_COEF_CONTINUE(val)                              \
  {                                                           \
    prob = coef_probs + (ENTROPY_NODES*PREV_CONTEXT_INC(val));\
    qcoeff_ptr[scan[c]] = (INT16) vp8_get_signed(br, val);    \
    c++;                                                      \
    continue;                                                 \
  }

#define ADJUST_COEF(prob, bits_count)  \
  do {                                 \
    if (vp8_read(br, prob))            \
      val += (UINT16)(1 << bits_count);\
  } while (0);

static int decode_coefs(VP8D_COMP *dx, const MACROBLOCKD *xd,
                        BOOL_DECODER* const br,
                        ENTROPY_CONTEXT *a, ENTROPY_CONTEXT *l,
                        PLANE_TYPE type,
#if CONFIG_HYBRIDTRANSFORM8X8 || CONFIG_HYBRIDTRANSFORM || CONFIG_HYBRIDTRANSFORM16X16
                        TX_TYPE tx_type,
#endif
                        int seg_eob, INT16 *qcoeff_ptr, int i,
                        const int *const scan, int block_type,
                        const int *coef_bands) {
  FRAME_CONTEXT *const fc = &dx->common.fc;
  int tmp, c = (type == PLANE_TYPE_Y_NO_DC);
  const vp8_prob *prob, *coef_probs;

  switch (block_type) {
    default:
    case TX_4X4:
      coef_probs =
#if CONFIG_HYBRIDTRANSFORM
        tx_type != DCT_DCT ? fc->hybrid_coef_probs[type][0][0] :
#endif
        fc->coef_probs[type][0][0];
      break;
    case TX_8X8:
      coef_probs =
#if CONFIG_HYBRIDTRANSFORM8X8
        tx_type != DCT_DCT ? fc->hybrid_coef_probs_8x8[type][0][0] :
#endif
        fc->coef_probs_8x8[type][0][0];
      break;
    case TX_16X16:
      coef_probs =
#if CONFIG_HYBRIDTRANSFORM16X16
        tx_type != DCT_DCT ? fc->hybrid_coef_probs_16x16[type][0][0] :
#endif
        fc->coef_probs_16x16[type][0][0];
      break;
  }

  VP8_COMBINEENTROPYCONTEXTS(tmp, *a, *l);
  prob = coef_probs + tmp * ENTROPY_NODES;

  while (1) {
    int val;
    const uint8_t *cat6 = cat6_prob;
    if (c == seg_eob) break;
    prob += coef_bands[c];
    if (!vp8_read(br, prob[EOB_CONTEXT_NODE]))
      break;
SKIP_START:
    if (c == seg_eob) break;
    if (!vp8_read(br, prob[ZERO_CONTEXT_NODE])) {
      ++c;
      prob = coef_probs + coef_bands[c];
      goto SKIP_START;
    }
    // ONE_CONTEXT_NODE_0_
    if (!vp8_read(br, prob[ONE_CONTEXT_NODE])) {
      prob = coef_probs + ENTROPY_NODES;
      qcoeff_ptr[scan[c]] = (INT16) vp8_get_signed(br, 1);
      ++c;
      continue;
    }
    // LOW_VAL_CONTEXT_NODE_0_
    if (!vp8_read(br, prob[LOW_VAL_CONTEXT_NODE])) {
      if (!vp8_read(br, prob[TWO_CONTEXT_NODE])) {
        WRITE_COEF_CONTINUE(2);
      }
      if (!vp8_read(br, prob[THREE_CONTEXT_NODE])) {
        WRITE_COEF_CONTINUE(3);
      }
      WRITE_COEF_CONTINUE(4);
    }
    // HIGH_LOW_CONTEXT_NODE_0_
    if (!vp8_read(br, prob[HIGH_LOW_CONTEXT_NODE])) {
      if (!vp8_read(br, prob[CAT_ONE_CONTEXT_NODE])) {
        val = CAT1_MIN_VAL;
        ADJUST_COEF(CAT1_PROB0, 0);
        WRITE_COEF_CONTINUE(val);
      }
      val = CAT2_MIN_VAL;
      ADJUST_COEF(CAT2_PROB1, 1);
      ADJUST_COEF(CAT2_PROB0, 0);
      WRITE_COEF_CONTINUE(val);
    }
    // CAT_THREEFOUR_CONTEXT_NODE_0_
    if (!vp8_read(br, prob[CAT_THREEFOUR_CONTEXT_NODE])) {
      if (!vp8_read(br, prob[CAT_THREE_CONTEXT_NODE])) {
        val = CAT3_MIN_VAL;
        ADJUST_COEF(CAT3_PROB2, 2);
        ADJUST_COEF(CAT3_PROB1, 1);
        ADJUST_COEF(CAT3_PROB0, 0);
        WRITE_COEF_CONTINUE(val);
      }
      val = CAT4_MIN_VAL;
      ADJUST_COEF(CAT4_PROB3, 3);
      ADJUST_COEF(CAT4_PROB2, 2);
      ADJUST_COEF(CAT4_PROB1, 1);
      ADJUST_COEF(CAT4_PROB0, 0);
      WRITE_COEF_CONTINUE(val);
    }
    // CAT_FIVE_CONTEXT_NODE_0_:
    if (!vp8_read(br, prob[CAT_FIVE_CONTEXT_NODE])) {
      val = CAT5_MIN_VAL;
      ADJUST_COEF(CAT5_PROB4, 4);
      ADJUST_COEF(CAT5_PROB3, 3);
      ADJUST_COEF(CAT5_PROB2, 2);
      ADJUST_COEF(CAT5_PROB1, 1);
      ADJUST_COEF(CAT5_PROB0, 0);
      WRITE_COEF_CONTINUE(val);
    }
    val = 0;
    while (*cat6) {
      val = (val << 1) | vp8_read(br, *cat6++);
    }
    val += CAT6_MIN_VAL;
    WRITE_COEF_CONTINUE(val);
  }

  if (block_type == TX_4X4) {
#if CONFIG_HYBRIDTRANSFORM
    count_tokens_adaptive_scan(xd, qcoeff_ptr, i, type,
                               tx_type,
                               a, l, c, seg_eob, fc);
#else
    count_tokens(qcoeff_ptr, i, type,
                 a, l, c, seg_eob, fc);
#endif
  }
  else if (block_type == TX_8X8)
    count_tokens_8x8(qcoeff_ptr, i, type,
#if CONFIG_HYBRIDTRANSFORM8X8
                     tx_type,
#endif
                     a, l, c, seg_eob, fc);
  else
    count_tokens_16x16(qcoeff_ptr, i, type,
#if CONFIG_HYBRIDTRANSFORM16X16
                       tx_type,
#endif
                       a, l, c, seg_eob, fc);
  return c;
}

int vp8_decode_mb_tokens_16x16(VP8D_COMP *pbi, MACROBLOCKD *xd,
                               BOOL_DECODER* const bc) {
  ENTROPY_CONTEXT* const A = (ENTROPY_CONTEXT *)xd->above_context;
  ENTROPY_CONTEXT* const L = (ENTROPY_CONTEXT *)xd->left_context;

  char* const eobs = xd->eobs;
  PLANE_TYPE type;
  int c, i, eobtotal = 0, seg_eob;
  const int segment_id = xd->mode_info_context->mbmi.segment_id;
  const int seg_active = segfeature_active(xd, segment_id, SEG_LVL_EOB);
  INT16 *qcoeff_ptr = &xd->qcoeff[0];
#if CONFIG_HYBRIDTRANSFORM8X8 || CONFIG_HYBRIDTRANSFORM || CONFIG_HYBRIDTRANSFORM16X16
  TX_TYPE tx_type = DCT_DCT;
#endif
#if CONFIG_HYBRIDTRANSFORM16X16
  tx_type = get_tx_type(xd, &xd->block[0]);
#endif

  type = PLANE_TYPE_Y_WITH_DC;

  if (seg_active)
      seg_eob = get_segdata(xd, segment_id, SEG_LVL_EOB);
  else
      seg_eob = 256;

  // Luma block
  {
    const int* const scan = vp8_default_zig_zag1d_16x16;
    c = decode_coefs(pbi, xd, bc, A, L, type,
#if CONFIG_HYBRIDTRANSFORM8X8 || CONFIG_HYBRIDTRANSFORM || CONFIG_HYBRIDTRANSFORM16X16
                     tx_type,
#endif
                     seg_eob, qcoeff_ptr,
                     0, scan, TX_16X16, coef_bands_x_16x16);
    eobs[0] = c;
    *A = *L = (c != !type);
    for (i = 1; i < 16; i++) {
      *(A + vp8_block2above[i]) = *(A);
      *(L +  vp8_block2left[i]) = *(L);
    }
    eobtotal += c;
  }

  // 8x8 chroma blocks
  qcoeff_ptr += 256;
  type = PLANE_TYPE_UV;
#if CONFIG_HYBRIDTRANSFORM8X8 || CONFIG_HYBRIDTRANSFORM || CONFIG_HYBRIDTRANSFORM16X16
  tx_type = DCT_DCT;
#endif
  if (seg_active)
    seg_eob = get_segdata(xd, segment_id, SEG_LVL_EOB);
  else
    seg_eob = 64;
  for (i = 16; i < 24; i += 4) {
    ENTROPY_CONTEXT* const a = A + vp8_block2above_8x8[i];
    ENTROPY_CONTEXT* const l = L + vp8_block2left_8x8[i];
    const int* const scan = vp8_default_zig_zag1d_8x8;

    c = decode_coefs(pbi, xd, bc, a, l, type,
#if CONFIG_HYBRIDTRANSFORM8X8 || CONFIG_HYBRIDTRANSFORM || CONFIG_HYBRIDTRANSFORM16X16
                     tx_type,
#endif
                     seg_eob, qcoeff_ptr,
                     i, scan, TX_8X8, coef_bands_x_8x8);
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

int vp8_decode_mb_tokens_8x8(VP8D_COMP *pbi, MACROBLOCKD *xd,
                             BOOL_DECODER* const bc) {
  ENTROPY_CONTEXT *const A = (ENTROPY_CONTEXT *)xd->above_context;
  ENTROPY_CONTEXT *const L = (ENTROPY_CONTEXT *)xd->left_context;

  char *const eobs = xd->eobs;
  PLANE_TYPE type;
  int c, i, eobtotal = 0, seg_eob;
  const int segment_id = xd->mode_info_context->mbmi.segment_id;
  const int seg_active = segfeature_active(xd, segment_id, SEG_LVL_EOB);
  INT16 *qcoeff_ptr = &xd->qcoeff[0];
#if CONFIG_HYBRIDTRANSFORM8X8 || CONFIG_HYBRIDTRANSFORM || CONFIG_HYBRIDTRANSFORM16X16
  TX_TYPE tx_type = DCT_DCT;
#endif

  int bufthred = (xd->mode_info_context->mbmi.mode == I8X8_PRED) ? 16 : 24;
  if (xd->mode_info_context->mbmi.mode != B_PRED &&
      xd->mode_info_context->mbmi.mode != SPLITMV &&
      xd->mode_info_context->mbmi.mode != I8X8_PRED) {
    ENTROPY_CONTEXT *const a = A + vp8_block2above_8x8[24];
    ENTROPY_CONTEXT *const l = L + vp8_block2left_8x8[24];
    const int *const scan = vp8_default_zig_zag1d;
    type = PLANE_TYPE_Y2;

    if (seg_active)
      seg_eob = get_segdata(xd, segment_id, SEG_LVL_EOB);
    else
      seg_eob = 4;
    c = decode_coefs(pbi, xd, bc, a, l, type,
#if CONFIG_HYBRIDTRANSFORM8X8 || CONFIG_HYBRIDTRANSFORM || CONFIG_HYBRIDTRANSFORM16X16
                     tx_type,
#endif
                     seg_eob, qcoeff_ptr + 24 * 16,
                     24, scan, TX_8X8, coef_bands_x);
    a[0] = l[0] = ((eobs[24] = c) != !type);

    eobtotal += c - 4;

    type = PLANE_TYPE_Y_NO_DC;
  } else
    type = PLANE_TYPE_Y_WITH_DC;

  if (seg_active)
    seg_eob = get_segdata(xd, segment_id, SEG_LVL_EOB);
  else
    seg_eob = 64;

  for (i = 0; i < bufthred ; i += 4) {
    ENTROPY_CONTEXT *const a = A + vp8_block2above_8x8[i];
    ENTROPY_CONTEXT *const l = L + vp8_block2left_8x8[i];
    const int *const scan = vp8_default_zig_zag1d_8x8;
#if CONFIG_HYBRIDTRANSFORM8X8 || CONFIG_HYBRIDTRANSFORM || CONFIG_HYBRIDTRANSFORM16X16
    tx_type = DCT_DCT;
#endif

    if (i == 16)
      type = PLANE_TYPE_UV;
#if CONFIG_HYBRIDTRANSFORM8X8
    if (type == PLANE_TYPE_Y_WITH_DC) {
      tx_type = get_tx_type(xd, xd->block + i);
    }
#endif

    c = decode_coefs(pbi, xd, bc, a, l, type,
#if CONFIG_HYBRIDTRANSFORM8X8 || CONFIG_HYBRIDTRANSFORM || CONFIG_HYBRIDTRANSFORM16X16
                     tx_type,
#endif
                     seg_eob, qcoeff_ptr,
                     i, scan, TX_8X8, coef_bands_x_8x8);
    a[0] = l[0] = ((eobs[i] = c) != !type);
    a[1] = a[0];
    l[1] = l[0];

    eobtotal += c;
    qcoeff_ptr += 64;
  }

  if (bufthred == 16) {
    type = PLANE_TYPE_UV;
#if CONFIG_HYBRIDTRANSFORM8X8 || CONFIG_HYBRIDTRANSFORM || CONFIG_HYBRIDTRANSFORM16X16
    tx_type = DCT_DCT;
#endif
    seg_eob = 16;

    // use 4x4 transform for U, V components in I8X8 prediction mode
    for (i = 16; i < 24; i++) {
      ENTROPY_CONTEXT *const a = A + vp8_block2above[i];
      ENTROPY_CONTEXT *const l = L + vp8_block2left[i];
      const int *scan = vp8_default_zig_zag1d;

      c = decode_coefs(pbi, xd, bc, a, l, type,
#if CONFIG_HYBRIDTRANSFORM8X8 || CONFIG_HYBRIDTRANSFORM || CONFIG_HYBRIDTRANSFORM16X16
                       tx_type,
#endif
                       seg_eob, qcoeff_ptr,
                       i, scan, TX_4X4, coef_bands_x);
      a[0] = l[0] = ((eobs[i] = c) != !type);

      eobtotal += c;
      qcoeff_ptr += 16;
    }
  }

  return eobtotal;
}


int vp8_decode_mb_tokens(VP8D_COMP *dx, MACROBLOCKD *xd,
                         BOOL_DECODER* const bc) {
  ENTROPY_CONTEXT *const A = (ENTROPY_CONTEXT *)xd->above_context;
  ENTROPY_CONTEXT *const L = (ENTROPY_CONTEXT *)xd->left_context;

  char *const eobs = xd->eobs;
  const int *scan = vp8_default_zig_zag1d;
  PLANE_TYPE type;
  int c, i, eobtotal = 0, seg_eob = 16;
  INT16 *qcoeff_ptr = &xd->qcoeff[0];

  int segment_id = xd->mode_info_context->mbmi.segment_id;
  if (segfeature_active(xd, segment_id, SEG_LVL_EOB))
    seg_eob = get_segdata(xd, segment_id, SEG_LVL_EOB);

  if (xd->mode_info_context->mbmi.mode != B_PRED &&
      xd->mode_info_context->mbmi.mode != I8X8_PRED &&
      xd->mode_info_context->mbmi.mode != SPLITMV) {
    ENTROPY_CONTEXT *const a = A + vp8_block2above[24];
    ENTROPY_CONTEXT *const l = L + vp8_block2left[24];
    type = PLANE_TYPE_Y2;

    c = decode_coefs(dx, xd, bc, a, l, type,
#if CONFIG_HYBRIDTRANSFORM8X8 || CONFIG_HYBRIDTRANSFORM || CONFIG_HYBRIDTRANSFORM16X16
                     DCT_DCT,
#endif
                     seg_eob, qcoeff_ptr + 24 * 16, 24,
                     scan, TX_4X4, coef_bands_x);
    a[0] = l[0] = ((eobs[24] = c) != !type);
    eobtotal += c - 16;

    type = PLANE_TYPE_Y_NO_DC;
  } else {
    type = PLANE_TYPE_Y_WITH_DC;
  }

  for (i = 0; i < 24; ++i) {
    ENTROPY_CONTEXT *const a = A + vp8_block2above[i];
    ENTROPY_CONTEXT *const l = L + vp8_block2left[i];
#if CONFIG_HYBRIDTRANSFORM8X8 || CONFIG_HYBRIDTRANSFORM || CONFIG_HYBRIDTRANSFORM16X16
    TX_TYPE tx_type = DCT_DCT;
#endif
    if (i == 16)
      type = PLANE_TYPE_UV;

#if CONFIG_HYBRIDTRANSFORM
    tx_type = get_tx_type(xd, &xd->block[i]);
    switch(tx_type) {
      case ADST_DCT :
        scan = vp8_row_scan;
        break;

      case DCT_ADST :
        scan = vp8_col_scan;
        break;

      default :
        scan = vp8_default_zig_zag1d;
        break;
    }
#endif

    c = decode_coefs(dx, xd, bc, a, l, type,
#if CONFIG_HYBRIDTRANSFORM8X8 || CONFIG_HYBRIDTRANSFORM || CONFIG_HYBRIDTRANSFORM16X16
                     tx_type,
#endif
                     seg_eob, qcoeff_ptr,
                     i, scan, TX_4X4, coef_bands_x);
    a[0] = l[0] = ((eobs[i] = c) != !type);

    eobtotal += c;
    qcoeff_ptr += 16;
  }

  return eobtotal;
}
