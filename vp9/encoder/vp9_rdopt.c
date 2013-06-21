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
#include <math.h>
#include <limits.h>
#include <assert.h>

#include "vp9/common/vp9_pragmas.h"
#include "vp9/encoder/vp9_tokenize.h"
#include "vp9/encoder/vp9_treewriter.h"
#include "vp9/encoder/vp9_onyx_int.h"
#include "vp9/encoder/vp9_modecosts.h"
#include "vp9/encoder/vp9_encodeintra.h"
#include "vp9/common/vp9_entropymode.h"
#include "vp9/common/vp9_reconinter.h"
#include "vp9/common/vp9_reconintra.h"
#include "vp9/common/vp9_findnearmv.h"
#include "vp9/common/vp9_quant_common.h"
#include "vp9/encoder/vp9_encodemb.h"
#include "vp9/encoder/vp9_quantize.h"
#include "vp9/encoder/vp9_variance.h"
#include "vp9/encoder/vp9_mcomp.h"
#include "vp9/encoder/vp9_rdopt.h"
#include "vp9/encoder/vp9_ratectrl.h"
#include "vpx_mem/vpx_mem.h"
#include "vp9/common/vp9_systemdependent.h"
#include "vp9/encoder/vp9_encodemv.h"
#include "vp9/common/vp9_seg_common.h"
#include "vp9/common/vp9_pred_common.h"
#include "vp9/common/vp9_entropy.h"
#include "vp9_rtcd.h"
#include "vp9/common/vp9_mvref_common.h"
#include "vp9/common/vp9_common.h"

#define INVALID_MV 0x80008000

/* Factor to weigh the rate for switchable interp filters */
#define SWITCHABLE_INTERP_RATE_FACTOR 1

DECLARE_ALIGNED(16, extern const uint8_t,
                vp9_pt_energy_class[MAX_ENTROPY_TOKENS]);

#define I4X4_PRED 0x8000
#define SPLITMV 0x10000

const MODE_DEFINITION vp9_mode_order[MAX_MODES] = {
  {ZEROMV,    LAST_FRAME,   NONE},
  {DC_PRED,   INTRA_FRAME,  NONE},

  {NEARESTMV, LAST_FRAME,   NONE},
  {NEARMV,    LAST_FRAME,   NONE},

  {ZEROMV,    GOLDEN_FRAME, NONE},
  {NEARESTMV, GOLDEN_FRAME, NONE},

  {ZEROMV,    ALTREF_FRAME, NONE},
  {NEARESTMV, ALTREF_FRAME, NONE},

  {NEARMV,    GOLDEN_FRAME, NONE},
  {NEARMV,    ALTREF_FRAME, NONE},

  {V_PRED,    INTRA_FRAME,  NONE},
  {H_PRED,    INTRA_FRAME,  NONE},
  {D45_PRED,  INTRA_FRAME,  NONE},
  {D135_PRED, INTRA_FRAME,  NONE},
  {D117_PRED, INTRA_FRAME,  NONE},
  {D153_PRED, INTRA_FRAME,  NONE},
  {D27_PRED,  INTRA_FRAME,  NONE},
  {D63_PRED,  INTRA_FRAME,  NONE},

  {TM_PRED,   INTRA_FRAME,  NONE},

  {NEWMV,     LAST_FRAME,   NONE},
  {NEWMV,     GOLDEN_FRAME, NONE},
  {NEWMV,     ALTREF_FRAME, NONE},

  {SPLITMV,   LAST_FRAME,   NONE},
  {SPLITMV,   GOLDEN_FRAME, NONE},
  {SPLITMV,   ALTREF_FRAME, NONE},

  {I4X4_PRED, INTRA_FRAME,  NONE},

  /* compound prediction modes */
  {ZEROMV,    LAST_FRAME,   ALTREF_FRAME},
  {NEARESTMV, LAST_FRAME,   ALTREF_FRAME},
  {NEARMV,    LAST_FRAME,   ALTREF_FRAME},

  {ZEROMV,    GOLDEN_FRAME, ALTREF_FRAME},
  {NEARESTMV, GOLDEN_FRAME, ALTREF_FRAME},
  {NEARMV,    GOLDEN_FRAME, ALTREF_FRAME},

  {NEWMV,     LAST_FRAME,   ALTREF_FRAME},
  {NEWMV,     GOLDEN_FRAME, ALTREF_FRAME},

  {SPLITMV,   LAST_FRAME,   ALTREF_FRAME},
  {SPLITMV,   GOLDEN_FRAME, ALTREF_FRAME},
};

// The baseline rd thresholds for breaking out of the rd loop for
// certain modes are assumed to be based on 8x8 blocks.
// This table is used to correct for blocks size.
// The factors here are << 2 (2 = x0.5, 32 = x8 etc).
static int rd_thresh_block_size_factor[BLOCK_SIZE_TYPES] =
  {2, 3, 3, 4, 6, 6, 8, 12, 12, 16, 24, 24, 32};

#define BASE_RD_THRESH_FREQ_FACT 16
#define MAX_RD_THRESH_FREQ_FACT 32
#define MAX_RD_THRESH_FREQ_INC 1

static void fill_token_costs(vp9_coeff_count (*c)[BLOCK_TYPES],
                             vp9_coeff_count (*cnoskip)[BLOCK_TYPES],
                             vp9_coeff_probs_model (*p)[BLOCK_TYPES]) {
  int i, j, k, l;
  TX_SIZE t;
  for (t = TX_4X4; t <= TX_32X32; t++)
    for (i = 0; i < BLOCK_TYPES; i++)
      for (j = 0; j < REF_TYPES; j++)
        for (k = 0; k < COEF_BANDS; k++)
          for (l = 0; l < PREV_COEF_CONTEXTS; l++) {
            vp9_prob probs[ENTROPY_NODES];
            vp9_model_to_full_probs(p[t][i][j][k][l], probs);
            vp9_cost_tokens((int *)cnoskip[t][i][j][k][l], probs,
                            vp9_coef_tree);
#if CONFIG_BALANCED_COEFTREE
            // Replace the eob node prob with a very small value so that the
            // cost approximately equals the cost without the eob node
            probs[1] = 1;
            vp9_cost_tokens((int *)c[t][i][j][k][l], probs, vp9_coef_tree);
#else
            vp9_cost_tokens_skip((int *)c[t][i][j][k][l], probs,
                                 vp9_coef_tree);
            assert(c[t][i][j][k][l][DCT_EOB_TOKEN] ==
                   cnoskip[t][i][j][k][l][DCT_EOB_TOKEN]);
#endif
          }
}

static int rd_iifactor[32] =  { 4, 4, 3, 2, 1, 0, 0, 0,
                                0, 0, 0, 0, 0, 0, 0, 0,
                                0, 0, 0, 0, 0, 0, 0, 0,
                                0, 0, 0, 0, 0, 0, 0, 0, };

// 3* dc_qlookup[Q]*dc_qlookup[Q];

/* values are now correlated to quantizer */
static int sad_per_bit16lut[QINDEX_RANGE];
static int sad_per_bit4lut[QINDEX_RANGE];

void vp9_init_me_luts() {
  int i;

  // Initialize the sad lut tables using a formulaic calculation for now
  // This is to make it easier to resolve the impact of experimental changes
  // to the quantizer tables.
  for (i = 0; i < QINDEX_RANGE; i++) {
    sad_per_bit16lut[i] =
      (int)((0.0418 * vp9_convert_qindex_to_q(i)) + 2.4107);
    sad_per_bit4lut[i] = (int)(0.063 * vp9_convert_qindex_to_q(i) + 2.742);
  }
}

static int compute_rd_mult(int qindex) {
  const int q = vp9_dc_quant(qindex, 0);
  return (11 * q * q) >> 2;
}

void vp9_initialize_me_consts(VP9_COMP *cpi, int qindex) {
  cpi->mb.sadperbit16 = sad_per_bit16lut[qindex];
  cpi->mb.sadperbit4 = sad_per_bit4lut[qindex];
}


void vp9_initialize_rd_consts(VP9_COMP *cpi, int qindex) {
  int q, i, bsize;

  vp9_clear_system_state();  // __asm emms;

  // Further tests required to see if optimum is different
  // for key frames, golden frames and arf frames.
  // if (cpi->common.refresh_golden_frame ||
  //     cpi->common.refresh_alt_ref_frame)
  qindex = clamp(qindex, 0, MAXQ);

  cpi->RDMULT = compute_rd_mult(qindex);
  if (cpi->pass == 2 && (cpi->common.frame_type != KEY_FRAME)) {
    if (cpi->twopass.next_iiratio > 31)
      cpi->RDMULT += (cpi->RDMULT * rd_iifactor[31]) >> 4;
    else
      cpi->RDMULT +=
          (cpi->RDMULT * rd_iifactor[cpi->twopass.next_iiratio]) >> 4;
  }
  cpi->mb.errorperbit = cpi->RDMULT >> 6;
  cpi->mb.errorperbit += (cpi->mb.errorperbit == 0);

  vp9_set_speed_features(cpi);

  q = (int)pow(vp9_dc_quant(qindex, 0) >> 2, 1.25);
  q <<= 2;
  if (q < 8)
    q = 8;

  if (cpi->RDMULT > 1000) {
    cpi->RDDIV = 1;
    cpi->RDMULT /= 100;

    for (bsize = 0; bsize < BLOCK_SIZE_TYPES; ++bsize) {
      for (i = 0; i < MAX_MODES; ++i) {
        // Threshold here seem unecessarily harsh but fine given actual
        // range of values used for cpi->sf.thresh_mult[]
        int thresh_max = INT_MAX / (q * rd_thresh_block_size_factor[bsize]);

        // *4 relates to the scaling of rd_thresh_block_size_factor[]
        if ((int64_t)cpi->sf.thresh_mult[i] < thresh_max) {
          cpi->rd_threshes[bsize][i] =
            cpi->sf.thresh_mult[i] * q *
            rd_thresh_block_size_factor[bsize] / (4 * 100);
        } else {
          cpi->rd_threshes[bsize][i] = INT_MAX;
        }
        cpi->rd_baseline_thresh[bsize][i] = cpi->rd_threshes[bsize][i];
        cpi->rd_thresh_freq_fact[bsize][i] = BASE_RD_THRESH_FREQ_FACT;
      }
    }
  } else {
    cpi->RDDIV = 100;

    for (bsize = 0; bsize < BLOCK_SIZE_TYPES; ++bsize) {
      for (i = 0; i < MAX_MODES; i++) {
        // Threshold here seem unecessarily harsh but fine given actual
        // range of values used for cpi->sf.thresh_mult[]
        int thresh_max = INT_MAX / (q * rd_thresh_block_size_factor[bsize]);

        if (cpi->sf.thresh_mult[i] < thresh_max) {
          cpi->rd_threshes[bsize][i] =
            cpi->sf.thresh_mult[i] * q *
            rd_thresh_block_size_factor[bsize] / 4;
        } else {
          cpi->rd_threshes[bsize][i] = INT_MAX;
        }
        cpi->rd_baseline_thresh[bsize][i] = cpi->rd_threshes[bsize][i];
        cpi->rd_thresh_freq_fact[bsize][i] = BASE_RD_THRESH_FREQ_FACT;
      }
    }
  }

  fill_token_costs(cpi->mb.token_costs,
                   cpi->mb.token_costs_noskip,
                   cpi->common.fc.coef_probs);

  for (i = 0; i < NUM_PARTITION_CONTEXTS; i++)
    vp9_cost_tokens(cpi->mb.partition_cost[i],
                    cpi->common.fc.partition_prob[cpi->common.frame_type][i],
                    vp9_partition_tree);

  /*rough estimate for costing*/
  vp9_init_mode_costs(cpi);

  if (cpi->common.frame_type != KEY_FRAME) {
    vp9_build_nmv_cost_table(
        cpi->mb.nmvjointcost,
        cpi->mb.e_mbd.allow_high_precision_mv ?
        cpi->mb.nmvcost_hp : cpi->mb.nmvcost,
        &cpi->common.fc.nmvc,
        cpi->mb.e_mbd.allow_high_precision_mv, 1, 1);
  }
}

int64_t vp9_block_error_c(int16_t *coeff, int16_t *dqcoeff,
                          intptr_t block_size) {
  int i;
  int64_t error = 0;

  for (i = 0; i < block_size; i++) {
    int this_diff = coeff[i] - dqcoeff[i];
    error += (unsigned)this_diff * this_diff;
  }

  return error;
}

static INLINE int cost_coeffs(VP9_COMMON *const cm, MACROBLOCK *mb,
                              int plane, int block, PLANE_TYPE type,
                              ENTROPY_CONTEXT *A,
                              ENTROPY_CONTEXT *L,
                              TX_SIZE tx_size,
                              int y_blocks) {
  MACROBLOCKD *const xd = &mb->e_mbd;
  MB_MODE_INFO *mbmi = &xd->mode_info_context->mbmi;
  int pt;
  int c = 0;
  int cost = 0, pad;
  const int *scan, *nb;
  const int eob = xd->plane[plane].eobs[block];
  const int16_t *qcoeff_ptr = BLOCK_OFFSET(xd->plane[plane].qcoeff,
                                           block, 16);
  const int ref = mbmi->ref_frame[0] != INTRA_FRAME;
  unsigned int (*token_costs)[PREV_COEF_CONTEXTS][MAX_ENTROPY_TOKENS] =
      mb->token_costs[tx_size][type][ref];
  ENTROPY_CONTEXT above_ec, left_ec;
  TX_TYPE tx_type = DCT_DCT;

  const int segment_id = xd->mode_info_context->mbmi.segment_id;
  unsigned int (*token_costs_noskip)[PREV_COEF_CONTEXTS][MAX_ENTROPY_TOKENS] =
      mb->token_costs_noskip[tx_size][type][ref];

  int seg_eob, default_eob;
  uint8_t token_cache[1024];
  const uint8_t * band_translate;

  // Check for consistency of tx_size with mode info
  assert((!type && !plane) || (type && plane));
  if (type == PLANE_TYPE_Y_WITH_DC) {
    assert(xd->mode_info_context->mbmi.txfm_size == tx_size);
  } else {
    TX_SIZE tx_size_uv = get_uv_tx_size(mbmi);
    assert(tx_size == tx_size_uv);
  }

  switch (tx_size) {
    case TX_4X4: {
      tx_type = (type == PLANE_TYPE_Y_WITH_DC) ?
          get_tx_type_4x4(xd, block) : DCT_DCT;
      above_ec = A[0] != 0;
      left_ec = L[0] != 0;
      seg_eob = 16;
      scan = get_scan_4x4(tx_type);
      band_translate = vp9_coefband_trans_4x4;
      break;
    }
    case TX_8X8: {
      const BLOCK_SIZE_TYPE sb_type = xd->mode_info_context->mbmi.sb_type;
      const int sz = 1 + b_width_log2(sb_type);
      const int x = block & ((1 << sz) - 1), y = block - x;
      TX_TYPE tx_type = (type == PLANE_TYPE_Y_WITH_DC) ?
          get_tx_type_8x8(xd, y + (x >> 1)) : DCT_DCT;
      above_ec = (A[0] + A[1]) != 0;
      left_ec = (L[0] + L[1]) != 0;
      scan = get_scan_8x8(tx_type);
      seg_eob = 64;
      band_translate = vp9_coefband_trans_8x8plus;
      break;
    }
    case TX_16X16: {
      const BLOCK_SIZE_TYPE sb_type = xd->mode_info_context->mbmi.sb_type;
      const int sz = 2 + b_width_log2(sb_type);
      const int x = block & ((1 << sz) - 1), y = block - x;
      TX_TYPE tx_type = (type == PLANE_TYPE_Y_WITH_DC) ?
          get_tx_type_16x16(xd, y + (x >> 2)) : DCT_DCT;
      scan = get_scan_16x16(tx_type);
      seg_eob = 256;
      above_ec = (A[0] + A[1] + A[2] + A[3]) != 0;
      left_ec = (L[0] + L[1] + L[2] + L[3]) != 0;
      band_translate = vp9_coefband_trans_8x8plus;
      break;
    }
    case TX_32X32:
      scan = vp9_default_scan_32x32;
      seg_eob = 1024;
      above_ec = (A[0] + A[1] + A[2] + A[3] + A[4] + A[5] + A[6] + A[7]) != 0;
      left_ec = (L[0] + L[1] + L[2] + L[3] + L[4] + L[5] + L[6] + L[7]) != 0;
      band_translate = vp9_coefband_trans_8x8plus;
      break;
    default:
      abort();
      break;
  }
  assert(eob <= seg_eob);

  pt = combine_entropy_contexts(above_ec, left_ec);
  nb = vp9_get_coef_neighbors_handle(scan, &pad);
  default_eob = seg_eob;

  if (vp9_segfeature_active(xd, segment_id, SEG_LVL_SKIP))
    seg_eob = 0;

  /* sanity check to ensure that we do not have spurious non-zero q values */
  if (eob < seg_eob)
    assert(qcoeff_ptr[scan[eob]] == 0);

  {
    for (c = 0; c < eob; c++) {
      int v = qcoeff_ptr[scan[c]];
      int t = vp9_dct_value_tokens_ptr[v].token;
      int band = get_coef_band(band_translate, c);
      if (c)
        pt = vp9_get_coef_context(scan, nb, pad, token_cache, c, default_eob);

      if (!c || token_cache[scan[c - 1]])  // do not skip eob
        cost += token_costs_noskip[band][pt][t] + vp9_dct_value_cost_ptr[v];
      else
        cost += token_costs[band][pt][t] + vp9_dct_value_cost_ptr[v];
      token_cache[scan[c]] = vp9_pt_energy_class[t];
    }
    if (c < seg_eob) {
      if (c)
        pt = vp9_get_coef_context(scan, nb, pad, token_cache, c, default_eob);
      cost += mb->token_costs_noskip[tx_size][type][ref]
          [get_coef_band(band_translate, c)]
          [pt][DCT_EOB_TOKEN];
    }
  }

  // is eob first coefficient;
  for (pt = 0; pt < (1 << tx_size); pt++) {
    A[pt] = L[pt] = c > 0;
  }

  return cost;
}

static void choose_txfm_size_from_rd(VP9_COMP *cpi, MACROBLOCK *x,
                                     int (*r)[2], int *rate,
                                     int64_t *d, int64_t *distortion,
                                     int *s, int *skip,
                                     int64_t txfm_cache[NB_TXFM_MODES],
                                     TX_SIZE max_txfm_size) {
  VP9_COMMON *const cm = &cpi->common;
  MACROBLOCKD *const xd = &x->e_mbd;
  MB_MODE_INFO *const mbmi = &xd->mode_info_context->mbmi;
  vp9_prob skip_prob = vp9_get_pred_prob(cm, xd, PRED_MBSKIP);
  int64_t rd[TX_SIZE_MAX_SB][2];
  int n, m;
  int s0, s1;

  const vp9_prob *tx_probs = vp9_get_pred_probs(cm, xd, PRED_TX_SIZE);

  for (n = TX_4X4; n <= max_txfm_size; n++) {
    r[n][1] = r[n][0];
    for (m = 0; m <= n - (n == max_txfm_size); m++) {
      if (m == n)
        r[n][1] += vp9_cost_zero(tx_probs[m]);
      else
        r[n][1] += vp9_cost_one(tx_probs[m]);
    }
  }

  assert(skip_prob > 0);
  s0 = vp9_cost_bit(skip_prob, 0);
  s1 = vp9_cost_bit(skip_prob, 1);

  for (n = TX_4X4; n <= max_txfm_size; n++) {
    if (s[n]) {
      rd[n][0] = rd[n][1] = RDCOST(x->rdmult, x->rddiv, s1, d[n]);
    } else {
      rd[n][0] = RDCOST(x->rdmult, x->rddiv, r[n][0] + s0, d[n]);
      rd[n][1] = RDCOST(x->rdmult, x->rddiv, r[n][1] + s0, d[n]);
    }
  }

  if (max_txfm_size == TX_32X32 &&
      (cm->txfm_mode == ALLOW_32X32 ||
       (cm->txfm_mode == TX_MODE_SELECT &&
        rd[TX_32X32][1] < rd[TX_16X16][1] && rd[TX_32X32][1] < rd[TX_8X8][1] &&
        rd[TX_32X32][1] < rd[TX_4X4][1]))) {
    mbmi->txfm_size = TX_32X32;
  } else if (max_txfm_size >= TX_16X16 &&
             (cm->txfm_mode == ALLOW_16X16 ||
              cm->txfm_mode == ALLOW_32X32 ||
              (cm->txfm_mode == TX_MODE_SELECT &&
               rd[TX_16X16][1] < rd[TX_8X8][1] &&
               rd[TX_16X16][1] < rd[TX_4X4][1]))) {
    mbmi->txfm_size = TX_16X16;
  } else if (cm->txfm_mode == ALLOW_8X8 ||
             cm->txfm_mode == ALLOW_16X16 ||
             cm->txfm_mode == ALLOW_32X32 ||
           (cm->txfm_mode == TX_MODE_SELECT && rd[TX_8X8][1] < rd[TX_4X4][1])) {
    mbmi->txfm_size = TX_8X8;
  } else {
    mbmi->txfm_size = TX_4X4;
  }

  *distortion = d[mbmi->txfm_size];
  *rate       = r[mbmi->txfm_size][cm->txfm_mode == TX_MODE_SELECT];
  *skip       = s[mbmi->txfm_size];

  txfm_cache[ONLY_4X4] = rd[TX_4X4][0];
  txfm_cache[ALLOW_8X8] = rd[TX_8X8][0];
  txfm_cache[ALLOW_16X16] = rd[MIN(max_txfm_size, TX_16X16)][0];
  txfm_cache[ALLOW_32X32] = rd[MIN(max_txfm_size, TX_32X32)][0];
  if (max_txfm_size == TX_32X32 &&
      rd[TX_32X32][1] < rd[TX_16X16][1] && rd[TX_32X32][1] < rd[TX_8X8][1] &&
      rd[TX_32X32][1] < rd[TX_4X4][1])
    txfm_cache[TX_MODE_SELECT] = rd[TX_32X32][1];
  else if (max_txfm_size >= TX_16X16 &&
           rd[TX_16X16][1] < rd[TX_8X8][1] && rd[TX_16X16][1] < rd[TX_4X4][1])
    txfm_cache[TX_MODE_SELECT] = rd[TX_16X16][1];
  else
    txfm_cache[TX_MODE_SELECT] = rd[TX_4X4][1] < rd[TX_8X8][1] ?
                                 rd[TX_4X4][1] : rd[TX_8X8][1];
}

static int64_t block_error_sby(MACROBLOCK *x, BLOCK_SIZE_TYPE bsize,
                               int shift) {
  const int bwl = b_width_log2(bsize), bhl = b_height_log2(bsize);
  return vp9_block_error(x->plane[0].coeff, x->e_mbd.plane[0].dqcoeff,
                         16 << (bwl + bhl)) >> shift;
}

static int64_t block_error_sbuv(MACROBLOCK *x, BLOCK_SIZE_TYPE bsize,
                                int shift) {
  const int bwl = b_width_log2(bsize), bhl = b_height_log2(bsize);
  int64_t sum = 0;
  int plane;

  for (plane = 1; plane < MAX_MB_PLANE; plane++) {
    const int subsampling = x->e_mbd.plane[plane].subsampling_x +
                            x->e_mbd.plane[plane].subsampling_y;
    sum += vp9_block_error(x->plane[plane].coeff, x->e_mbd.plane[plane].dqcoeff,
                           16 << (bwl + bhl - subsampling));
  }
  return sum >> shift;
}

struct rdcost_block_args {
  VP9_COMMON *cm;
  MACROBLOCK *x;
  ENTROPY_CONTEXT t_above[16];
  ENTROPY_CONTEXT t_left[16];
  TX_SIZE tx_size;
  int bw;
  int bh;
  int cost;
};

static void rdcost_block(int plane, int block, BLOCK_SIZE_TYPE bsize,
                         int ss_txfrm_size, void *arg) {
  struct rdcost_block_args* args = arg;
  int x_idx, y_idx;
  MACROBLOCKD * const xd = &args->x->e_mbd;

  txfrm_block_to_raster_xy(xd, bsize, plane, block, args->tx_size * 2, &x_idx,
                           &y_idx);

  args->cost += cost_coeffs(args->cm, args->x, plane, block,
                            xd->plane[plane].plane_type, args->t_above + x_idx,
                            args->t_left + y_idx, args->tx_size,
                            args->bw * args->bh);
}

static int rdcost_plane(VP9_COMMON * const cm, MACROBLOCK *x, int plane,
                        BLOCK_SIZE_TYPE bsize, TX_SIZE tx_size) {
  MACROBLOCKD * const xd = &x->e_mbd;
  const int bwl = b_width_log2(bsize) - xd->plane[plane].subsampling_x;
  const int bhl = b_height_log2(bsize) - xd->plane[plane].subsampling_y;
  const int bw = 1 << bwl, bh = 1 << bhl;
  struct rdcost_block_args args = { cm, x, { 0 }, { 0 }, tx_size, bw, bh, 0 };

  vpx_memcpy(&args.t_above, xd->plane[plane].above_context,
             sizeof(ENTROPY_CONTEXT) * bw);
  vpx_memcpy(&args.t_left, xd->plane[plane].left_context,
             sizeof(ENTROPY_CONTEXT) * bh);

  foreach_transformed_block_in_plane(xd, bsize, plane, rdcost_block, &args);

  return args.cost;
}

static int rdcost_uv(VP9_COMMON *const cm, MACROBLOCK *x,
                     BLOCK_SIZE_TYPE bsize, TX_SIZE tx_size) {
  int cost = 0, plane;

  for (plane = 1; plane < MAX_MB_PLANE; plane++) {
    cost += rdcost_plane(cm, x, plane, bsize, tx_size);
  }
  return cost;
}

static void super_block_yrd_for_txfm(VP9_COMMON *const cm, MACROBLOCK *x,
                                     int *rate, int64_t *distortion,
                                     int *skippable,
                                     BLOCK_SIZE_TYPE bsize, TX_SIZE tx_size) {
  MACROBLOCKD *const xd = &x->e_mbd;
  xd->mode_info_context->mbmi.txfm_size = tx_size;

  if (xd->mode_info_context->mbmi.ref_frame[0] == INTRA_FRAME)
    vp9_encode_intra_block_y(cm, x, bsize);
  else
    vp9_xform_quant_sby(cm, x, bsize);

  *distortion = block_error_sby(x, bsize, tx_size == TX_32X32 ? 0 : 2);
  *rate       = rdcost_plane(cm, x, 0, bsize, tx_size);
  *skippable  = vp9_sby_is_skippable(xd, bsize);
}

static void super_block_yrd(VP9_COMP *cpi,
                            MACROBLOCK *x, int *rate, int64_t *distortion,
                            int *skip, BLOCK_SIZE_TYPE bs,
                            int64_t txfm_cache[NB_TXFM_MODES]) {
  VP9_COMMON *const cm = &cpi->common;
  int r[TX_SIZE_MAX_SB][2], s[TX_SIZE_MAX_SB];
  int64_t d[TX_SIZE_MAX_SB];
  MACROBLOCKD *xd = &x->e_mbd;
  MB_MODE_INFO *const mbmi = &xd->mode_info_context->mbmi;

  assert(bs == mbmi->sb_type);
  if (mbmi->ref_frame[0] > INTRA_FRAME)
    vp9_subtract_sby(x, bs);

  if (cpi->sf.use_largest_txform) {
    if (bs >= BLOCK_SIZE_SB32X32) {
      mbmi->txfm_size = TX_32X32;
    } else if (bs >= BLOCK_SIZE_MB16X16) {
      mbmi->txfm_size = TX_16X16;
    } else if (bs >= BLOCK_SIZE_SB8X8) {
      mbmi->txfm_size = TX_8X8;
    } else {
      mbmi->txfm_size = TX_4X4;
    }
    vpx_memset(txfm_cache, 0, NB_TXFM_MODES * sizeof(int64_t));
    super_block_yrd_for_txfm(cm, x, rate, distortion, skip, bs,
                             mbmi->txfm_size);
    return;
  }
  if (bs >= BLOCK_SIZE_SB32X32)
    super_block_yrd_for_txfm(cm, x, &r[TX_32X32][0], &d[TX_32X32], &s[TX_32X32],
                             bs, TX_32X32);
  if (bs >= BLOCK_SIZE_MB16X16)
    super_block_yrd_for_txfm(cm, x, &r[TX_16X16][0], &d[TX_16X16], &s[TX_16X16],
                             bs, TX_16X16);
  super_block_yrd_for_txfm(cm, x, &r[TX_8X8][0], &d[TX_8X8], &s[TX_8X8], bs,
                           TX_8X8);
  super_block_yrd_for_txfm(cm, x, &r[TX_4X4][0], &d[TX_4X4], &s[TX_4X4], bs,
                           TX_4X4);

  choose_txfm_size_from_rd(cpi, x, r, rate, d, distortion, s,
                           skip, txfm_cache,
                           TX_32X32 - (bs < BLOCK_SIZE_SB32X32)
                           - (bs < BLOCK_SIZE_MB16X16));
}

static int64_t rd_pick_intra4x4block(VP9_COMP *cpi, MACROBLOCK *x, int ib,
                                     MB_PREDICTION_MODE *best_mode,
                                     int *bmode_costs,
                                     ENTROPY_CONTEXT *a, ENTROPY_CONTEXT *l,
                                     int *bestrate, int *bestratey,
                                     int64_t *bestdistortion,
                                     BLOCK_SIZE_TYPE bsize) {
  MB_PREDICTION_MODE mode;
  MACROBLOCKD *xd = &x->e_mbd;
  int64_t best_rd = INT64_MAX;
  int rate = 0;
  int64_t distortion;
  VP9_COMMON *const cm = &cpi->common;
  const int src_stride = x->plane[0].src.stride;
  uint8_t *src, *dst;
  int16_t *src_diff, *coeff;

  ENTROPY_CONTEXT ta[2], tempa[2];
  ENTROPY_CONTEXT tl[2], templ[2];
  TX_TYPE tx_type = DCT_DCT;
  TX_TYPE best_tx_type = DCT_DCT;
  int bw = 1 << b_width_log2(bsize);
  int bh = 1 << b_height_log2(bsize);
  int idx, idy, block;
  DECLARE_ALIGNED(16, int16_t, best_dqcoeff[4][16]);

  assert(ib < 4);

  vpx_memcpy(ta, a, sizeof(ta));
  vpx_memcpy(tl, l, sizeof(tl));
  xd->mode_info_context->mbmi.txfm_size = TX_4X4;

  for (mode = DC_PRED; mode <= TM_PRED; ++mode) {
    int64_t this_rd;
    int ratey = 0;

    rate = bmode_costs[mode];
    distortion = 0;

    vpx_memcpy(tempa, ta, sizeof(ta));
    vpx_memcpy(templ, tl, sizeof(tl));

    for (idy = 0; idy < bh; ++idy) {
      for (idx = 0; idx < bw; ++idx) {
        block = ib + idy * 2 + idx;
        xd->mode_info_context->bmi[block].as_mode.first = mode;
        src = raster_block_offset_uint8(xd, BLOCK_SIZE_SB8X8, 0, block,
                                        x->plane[0].src.buf, src_stride);
        src_diff = raster_block_offset_int16(xd, BLOCK_SIZE_SB8X8, 0, block,
                                             x->plane[0].src_diff);
        coeff = BLOCK_OFFSET(x->plane[0].coeff, block, 16);
        dst = raster_block_offset_uint8(xd, BLOCK_SIZE_SB8X8, 0, block,
                                        xd->plane[0].dst.buf,
                                        xd->plane[0].dst.stride);
        vp9_intra4x4_predict(xd, block, BLOCK_SIZE_SB8X8, mode,
                             dst, xd->plane[0].dst.stride);
        vp9_subtract_block(4, 4, src_diff, 8,
                           src, src_stride,
                           dst, xd->plane[0].dst.stride);

        tx_type = get_tx_type_4x4(xd, block);
        if (tx_type != DCT_DCT) {
          vp9_short_fht4x4(src_diff, coeff, 8, tx_type);
          x->quantize_b_4x4(x, block, tx_type, 16);
        } else {
          x->fwd_txm4x4(src_diff, coeff, 16);
          x->quantize_b_4x4(x, block, tx_type, 16);
        }

        ratey += cost_coeffs(cm, x, 0, block, PLANE_TYPE_Y_WITH_DC,
                             tempa + idx, templ + idy, TX_4X4, 16);
        distortion += vp9_block_error(coeff, BLOCK_OFFSET(xd->plane[0].dqcoeff,
                                                         block, 16), 16) >> 2;

        if (best_tx_type != DCT_DCT)
          vp9_short_iht4x4_add(BLOCK_OFFSET(xd->plane[0].dqcoeff, block, 16),
                               dst, xd->plane[0].dst.stride, best_tx_type);
        else
          xd->inv_txm4x4_add(BLOCK_OFFSET(xd->plane[0].dqcoeff, block, 16),
                             dst, xd->plane[0].dst.stride);
      }
    }

    rate += ratey;
    this_rd = RDCOST(x->rdmult, x->rddiv, rate, distortion);

    if (this_rd < best_rd) {
      *bestrate = rate;
      *bestratey = ratey;
      *bestdistortion = distortion;
      best_rd = this_rd;
      *best_mode = mode;
      best_tx_type = tx_type;
      vpx_memcpy(a, tempa, sizeof(tempa));
      vpx_memcpy(l, templ, sizeof(templ));
      for (idy = 0; idy < bh; ++idy) {
        for (idx = 0; idx < bw; ++idx) {
          block = ib + idy * 2 + idx;
          vpx_memcpy(best_dqcoeff[idy * 2 + idx],
                     BLOCK_OFFSET(xd->plane[0].dqcoeff, block, 16),
                     sizeof(best_dqcoeff[0]));
        }
      }
    }
  }

  for (idy = 0; idy < bh; ++idy) {
    for (idx = 0; idx < bw; ++idx) {
      block = ib + idy * 2 + idx;
      xd->mode_info_context->bmi[block].as_mode.first = *best_mode;
      dst = raster_block_offset_uint8(xd, BLOCK_SIZE_SB8X8, 0, block,
                                      xd->plane[0].dst.buf,
                                      xd->plane[0].dst.stride);

      vp9_intra4x4_predict(xd, block, BLOCK_SIZE_SB8X8, *best_mode,
                           dst, xd->plane[0].dst.stride);
      // inverse transform
      if (best_tx_type != DCT_DCT)
        vp9_short_iht4x4_add(best_dqcoeff[idy * 2 + idx], dst,
                             xd->plane[0].dst.stride, best_tx_type);
      else
        xd->inv_txm4x4_add(best_dqcoeff[idy * 2 + idx], dst,
                           xd->plane[0].dst.stride);
    }
  }

  return best_rd;
}

static int64_t rd_pick_intra4x4mby_modes(VP9_COMP *cpi, MACROBLOCK *mb,
                                         int *Rate, int *rate_y,
                                         int64_t *Distortion, int64_t best_rd) {
  int i, j;
  MACROBLOCKD *const xd = &mb->e_mbd;
  BLOCK_SIZE_TYPE bsize = xd->mode_info_context->mbmi.sb_type;
  int bw = 1 << b_width_log2(bsize);
  int bh = 1 << b_height_log2(bsize);
  int idx, idy;
  int cost = 0;
  int64_t distortion = 0;
  int tot_rate_y = 0;
  int64_t total_rd = 0;
  ENTROPY_CONTEXT t_above[4], t_left[4];
  int *bmode_costs;
  MODE_INFO *const mic = xd->mode_info_context;

  vpx_memcpy(t_above, xd->plane[0].above_context, sizeof(t_above));
  vpx_memcpy(t_left, xd->plane[0].left_context, sizeof(t_left));

  bmode_costs = mb->mbmode_cost;

  for (idy = 0; idy < 2; idy += bh) {
    for (idx = 0; idx < 2; idx += bw) {
      const int mis = xd->mode_info_stride;
      MB_PREDICTION_MODE UNINITIALIZED_IS_SAFE(best_mode);
      int UNINITIALIZED_IS_SAFE(r), UNINITIALIZED_IS_SAFE(ry);
      int64_t UNINITIALIZED_IS_SAFE(d);
      i = idy * 2 + idx;

      if (xd->frame_type == KEY_FRAME) {
        const MB_PREDICTION_MODE A = above_block_mode(mic, i, mis);
        const MB_PREDICTION_MODE L = (xd->left_available || idx) ?
                                     left_block_mode(mic, i) : DC_PRED;

        bmode_costs  = mb->y_mode_costs[A][L];
      }

      total_rd += rd_pick_intra4x4block(cpi, mb, i, &best_mode, bmode_costs,
                                        t_above + idx, t_left + idy,
                                        &r, &ry, &d, bsize);
      cost += r;
      distortion += d;
      tot_rate_y += ry;

      mic->bmi[i].as_mode.first = best_mode;
      for (j = 1; j < bh; ++j)
        mic->bmi[i + j * 2].as_mode.first = best_mode;
      for (j = 1; j < bw; ++j)
        mic->bmi[i + j].as_mode.first = best_mode;

      if (total_rd >= best_rd)
        break;
    }
  }

  if (total_rd >= best_rd)
    return INT64_MAX;

  *Rate = cost;
  *rate_y = tot_rate_y;
  *Distortion = distortion;
  xd->mode_info_context->mbmi.mode = mic->bmi[3].as_mode.first;

  return RDCOST(mb->rdmult, mb->rddiv, cost, distortion);
}

static int64_t rd_pick_intra_sby_mode(VP9_COMP *cpi, MACROBLOCK *x,
                                      int *rate, int *rate_tokenonly,
                                      int64_t *distortion, int *skippable,
                                      BLOCK_SIZE_TYPE bsize,
                                      int64_t txfm_cache[NB_TXFM_MODES]) {
  MB_PREDICTION_MODE mode;
  MB_PREDICTION_MODE UNINITIALIZED_IS_SAFE(mode_selected);
  MACROBLOCKD *const xd = &x->e_mbd;
  int this_rate, this_rate_tokenonly, s;
  int64_t this_distortion;
  int64_t best_rd = INT64_MAX, this_rd;
  TX_SIZE UNINITIALIZED_IS_SAFE(best_tx);
  int i;
  int *bmode_costs = x->mbmode_cost;

  if (bsize < BLOCK_SIZE_SB8X8) {
    x->e_mbd.mode_info_context->mbmi.txfm_size = TX_4X4;
    return best_rd;
  }

  for (i = 0; i < NB_TXFM_MODES; i++)
    txfm_cache[i] = INT64_MAX;

  /* Y Search for 32x32 intra prediction mode */
  for (mode = DC_PRED; mode <= TM_PRED; mode++) {
    int64_t local_txfm_cache[NB_TXFM_MODES];
    MODE_INFO *const mic = xd->mode_info_context;
    const int mis = xd->mode_info_stride;

    if (cpi->common.frame_type == KEY_FRAME) {
      const MB_PREDICTION_MODE A = above_block_mode(mic, 0, mis);
      const MB_PREDICTION_MODE L = xd->left_available ?
                                   left_block_mode(mic, 0) : DC_PRED;

      bmode_costs = x->y_mode_costs[A][L];
    }
    x->e_mbd.mode_info_context->mbmi.mode = mode;

    super_block_yrd(cpi, x, &this_rate_tokenonly, &this_distortion, &s,
                    bsize, local_txfm_cache);

    this_rate = this_rate_tokenonly + bmode_costs[mode];
    this_rd = RDCOST(x->rdmult, x->rddiv, this_rate, this_distortion);

    if (this_rd < best_rd) {
      mode_selected   = mode;
      best_rd         = this_rd;
      best_tx         = x->e_mbd.mode_info_context->mbmi.txfm_size;
      *rate           = this_rate;
      *rate_tokenonly = this_rate_tokenonly;
      *distortion     = this_distortion;
      *skippable      = s;
    }

    for (i = 0; i < NB_TXFM_MODES; i++) {
      int64_t adj_rd = this_rd + local_txfm_cache[i] -
                       local_txfm_cache[cpi->common.txfm_mode];
      if (adj_rd < txfm_cache[i]) {
        txfm_cache[i] = adj_rd;
      }
    }
  }

  x->e_mbd.mode_info_context->mbmi.mode = mode_selected;
  x->e_mbd.mode_info_context->mbmi.txfm_size = best_tx;

  return best_rd;
}

static void super_block_uvrd_for_txfm(VP9_COMMON *const cm, MACROBLOCK *x,
                                      int *rate, int64_t *distortion,
                                      int *skippable, BLOCK_SIZE_TYPE bsize,
                                      TX_SIZE uv_tx_size) {
  MACROBLOCKD *const xd = &x->e_mbd;
  if (xd->mode_info_context->mbmi.ref_frame[0] == INTRA_FRAME)
    vp9_encode_intra_block_uv(cm, x, bsize);
  else
    vp9_xform_quant_sbuv(cm, x, bsize);

  *distortion = block_error_sbuv(x, bsize, uv_tx_size == TX_32X32 ? 0 : 2);
  *rate       = rdcost_uv(cm, x, bsize, uv_tx_size);
  *skippable  = vp9_sbuv_is_skippable(xd, bsize);
}

static void super_block_uvrd(VP9_COMMON *const cm, MACROBLOCK *x,
                             int *rate, int64_t *distortion, int *skippable,
                             BLOCK_SIZE_TYPE bsize) {
  MACROBLOCKD *const xd = &x->e_mbd;
  MB_MODE_INFO *const mbmi = &xd->mode_info_context->mbmi;

  if (mbmi->ref_frame[0] > INTRA_FRAME)
    vp9_subtract_sbuv(x, bsize);

  if (mbmi->txfm_size >= TX_32X32 && bsize >= BLOCK_SIZE_SB64X64) {
    super_block_uvrd_for_txfm(cm, x, rate, distortion, skippable, bsize,
                              TX_32X32);
  } else if (mbmi->txfm_size >= TX_16X16 && bsize >= BLOCK_SIZE_SB32X32) {
    super_block_uvrd_for_txfm(cm, x, rate, distortion, skippable, bsize,
                              TX_16X16);
  } else if (mbmi->txfm_size >= TX_8X8 && bsize >= BLOCK_SIZE_MB16X16) {
    super_block_uvrd_for_txfm(cm, x, rate, distortion, skippable, bsize,
                              TX_8X8);
  } else {
    super_block_uvrd_for_txfm(cm, x, rate, distortion, skippable, bsize,
                              TX_4X4);
  }
}

static int64_t rd_pick_intra_sbuv_mode(VP9_COMP *cpi, MACROBLOCK *x,
                                       int *rate, int *rate_tokenonly,
                                       int64_t *distortion, int *skippable,
                                       BLOCK_SIZE_TYPE bsize) {
  MB_PREDICTION_MODE mode;
  MB_PREDICTION_MODE UNINITIALIZED_IS_SAFE(mode_selected);
  int64_t best_rd = INT64_MAX, this_rd;
  int this_rate_tokenonly, this_rate, s;
  int64_t this_distortion;

  for (mode = DC_PRED; mode <= TM_PRED; mode++) {
    x->e_mbd.mode_info_context->mbmi.uv_mode = mode;
    super_block_uvrd(&cpi->common, x, &this_rate_tokenonly,
                     &this_distortion, &s, bsize);
    this_rate = this_rate_tokenonly +
                x->intra_uv_mode_cost[x->e_mbd.frame_type][mode];
    this_rd = RDCOST(x->rdmult, x->rddiv, this_rate, this_distortion);

    if (this_rd < best_rd) {
      mode_selected   = mode;
      best_rd         = this_rd;
      *rate           = this_rate;
      *rate_tokenonly = this_rate_tokenonly;
      *distortion     = this_distortion;
      *skippable      = s;
    }
  }

  x->e_mbd.mode_info_context->mbmi.uv_mode = mode_selected;

  return best_rd;
}

int vp9_cost_mv_ref(VP9_COMP *cpi,
                    MB_PREDICTION_MODE m,
                    const int mode_context) {
  MACROBLOCKD *xd = &cpi->mb.e_mbd;
  int segment_id = xd->mode_info_context->mbmi.segment_id;

  // Dont account for mode here if segment skip is enabled.
  if (!vp9_segfeature_active(xd, segment_id, SEG_LVL_SKIP)) {
    VP9_COMMON *pc = &cpi->common;
    assert(NEARESTMV <= m  &&  m <= NEWMV);
    return cost_token(vp9_sb_mv_ref_tree,
                      pc->fc.inter_mode_probs[mode_context],
                      vp9_sb_mv_ref_encoding_array - NEARESTMV + m);
  } else
    return 0;
}

void vp9_set_mbmode_and_mvs(MACROBLOCK *x, MB_PREDICTION_MODE mb, int_mv *mv) {
  x->e_mbd.mode_info_context->mbmi.mode = mb;
  x->e_mbd.mode_info_context->mbmi.mv[0].as_int = mv->as_int;
}

static void joint_motion_search(VP9_COMP *cpi, MACROBLOCK *x,
                                BLOCK_SIZE_TYPE bsize,
                                int_mv *frame_mv,
                                int mi_row, int mi_col,
                                int_mv single_newmv[MAX_REF_FRAMES],
                                int *rate_mv);
static void single_motion_search(VP9_COMP *cpi, MACROBLOCK *x,
                                 BLOCK_SIZE_TYPE bsize,
                                 int mi_row, int mi_col,
                                 int_mv *tmp_mv, int *rate_mv);

static int labels2mode(MACROBLOCK *x, int i,
                       MB_PREDICTION_MODE this_mode,
                       int_mv *this_mv, int_mv *this_second_mv,
                       int_mv frame_mv[MB_MODE_COUNT][MAX_REF_FRAMES],
                       int_mv seg_mvs[MAX_REF_FRAMES],
                       int_mv *best_ref_mv,
                       int_mv *second_best_ref_mv,
                       int *mvjcost, int *mvcost[2], VP9_COMP *cpi) {
  MACROBLOCKD *const xd = &x->e_mbd;
  MODE_INFO *const mic = xd->mode_info_context;
  MB_MODE_INFO * mbmi = &mic->mbmi;
  int cost = 0, thismvcost = 0;
  int idx, idy;
  int bw = 1 << b_width_log2(mbmi->sb_type);
  int bh = 1 << b_height_log2(mbmi->sb_type);

  /* We have to be careful retrieving previously-encoded motion vectors.
   Ones from this macroblock have to be pulled from the BLOCKD array
   as they have not yet made it to the bmi array in our MB_MODE_INFO. */
  MB_PREDICTION_MODE m;

  // the only time we should do costing for new motion vector or mode
  // is when we are on a new label  (jbb May 08, 2007)
  switch (m = this_mode) {
    case NEWMV:
      this_mv->as_int = seg_mvs[mbmi->ref_frame[0]].as_int;
      thismvcost  = vp9_mv_bit_cost(this_mv, best_ref_mv, mvjcost, mvcost,
                                    102, xd->allow_high_precision_mv);
      if (mbmi->ref_frame[1] > 0) {
        this_second_mv->as_int = seg_mvs[mbmi->ref_frame[1]].as_int;
        thismvcost += vp9_mv_bit_cost(this_second_mv, second_best_ref_mv,
                                      mvjcost, mvcost, 102,
                                      xd->allow_high_precision_mv);
      }
      break;
    case NEARESTMV:
      this_mv->as_int = frame_mv[NEARESTMV][mbmi->ref_frame[0]].as_int;
      if (mbmi->ref_frame[1] > 0)
        this_second_mv->as_int =
            frame_mv[NEARESTMV][mbmi->ref_frame[1]].as_int;
      break;
    case NEARMV:
      this_mv->as_int = frame_mv[NEARMV][mbmi->ref_frame[0]].as_int;
      if (mbmi->ref_frame[1] > 0)
        this_second_mv->as_int =
            frame_mv[NEARMV][mbmi->ref_frame[1]].as_int;
      break;
    case ZEROMV:
      this_mv->as_int = 0;
      if (mbmi->ref_frame[1] > 0)
        this_second_mv->as_int = 0;
      break;
    default:
      break;
  }

  cost = vp9_cost_mv_ref(cpi, this_mode,
                         mbmi->mb_mode_context[mbmi->ref_frame[0]]);

  mic->bmi[i].as_mv[0].as_int = this_mv->as_int;
  if (mbmi->ref_frame[1] > 0)
    mic->bmi[i].as_mv[1].as_int = this_second_mv->as_int;

  x->partition_info->bmi[i].mode = m;
  x->partition_info->bmi[i].mv.as_int = this_mv->as_int;
  if (mbmi->ref_frame[1] > 0)
    x->partition_info->bmi[i].second_mv.as_int = this_second_mv->as_int;
  for (idy = 0; idy < bh; ++idy) {
    for (idx = 0; idx < bw; ++idx) {
      vpx_memcpy(&mic->bmi[i + idy * 2 + idx],
                 &mic->bmi[i], sizeof(mic->bmi[i]));
      vpx_memcpy(&x->partition_info->bmi[i + idy * 2 + idx],
                 &x->partition_info->bmi[i],
                 sizeof(x->partition_info->bmi[i]));
    }
  }

  cost += thismvcost;
  return cost;
}

static int64_t encode_inter_mb_segment(VP9_COMMON *const cm,
                                       MACROBLOCK *x,
                                       int i,
                                       int *labelyrate,
                                       int64_t *distortion,
                                       ENTROPY_CONTEXT *ta,
                                       ENTROPY_CONTEXT *tl) {
  int k;
  MACROBLOCKD *xd = &x->e_mbd;
  BLOCK_SIZE_TYPE bsize = xd->mode_info_context->mbmi.sb_type;
  int bwl = b_width_log2(bsize), bw = 1 << bwl;
  int bhl = b_height_log2(bsize), bh = 1 << bhl;
  int idx, idy;
  const int src_stride = x->plane[0].src.stride;
  uint8_t* const src =
  raster_block_offset_uint8(xd, BLOCK_SIZE_SB8X8, 0, i,
                            x->plane[0].src.buf, src_stride);
  int16_t* src_diff =
  raster_block_offset_int16(xd, BLOCK_SIZE_SB8X8, 0, i,
                            x->plane[0].src_diff);
  int16_t* coeff = BLOCK_OFFSET(x->plane[0].coeff, 16, i);
  uint8_t* const pre =
  raster_block_offset_uint8(xd, BLOCK_SIZE_SB8X8, 0, i,
                            xd->plane[0].pre[0].buf,
                            xd->plane[0].pre[0].stride);
  uint8_t* const dst =
  raster_block_offset_uint8(xd, BLOCK_SIZE_SB8X8, 0, i,
                            xd->plane[0].dst.buf,
                            xd->plane[0].dst.stride);
  int64_t thisdistortion = 0;
  int thisrate = 0;

  *labelyrate = 0;
  *distortion = 0;

  vp9_build_inter_predictor(pre,
                            xd->plane[0].pre[0].stride,
                            dst,
                            xd->plane[0].dst.stride,
                            &xd->mode_info_context->bmi[i].as_mv[0],
                            &xd->scale_factor[0],
                            4 * bw, 4 * bh, 0 /* no avg */, &xd->subpix);

  // TODO(debargha): Make this work properly with the
  // implicit-compoundinter-weight experiment when implicit
  // weighting for splitmv modes is turned on.
  if (xd->mode_info_context->mbmi.ref_frame[1] > 0) {
    uint8_t* const second_pre =
    raster_block_offset_uint8(xd, BLOCK_SIZE_SB8X8, 0, i,
                              xd->plane[0].pre[1].buf,
                              xd->plane[0].pre[1].stride);
    vp9_build_inter_predictor(second_pre, xd->plane[0].pre[1].stride,
                              dst, xd->plane[0].dst.stride,
                              &xd->mode_info_context->bmi[i].as_mv[1],
                              &xd->scale_factor[1], 4 * bw, 4 * bh, 1,
                              &xd->subpix);
  }

  vp9_subtract_block(4 * bh, 4 * bw, src_diff, 8,
                     src, src_stride,
                     dst, xd->plane[0].dst.stride);

  k = i;
  for (idy = 0; idy < bh; ++idy) {
    for (idx = 0; idx < bw; ++idx) {
      k += (idy * 2 + idx);
      src_diff = raster_block_offset_int16(xd, BLOCK_SIZE_SB8X8, 0, k,
                                           x->plane[0].src_diff);
      coeff = BLOCK_OFFSET(x->plane[0].coeff, 16, k);
      x->fwd_txm4x4(src_diff, coeff, 16);
      x->quantize_b_4x4(x, k, DCT_DCT, 16);
      thisdistortion += vp9_block_error(coeff,
                                        BLOCK_OFFSET(xd->plane[0].dqcoeff,
                                                     k, 16), 16);
      thisrate += cost_coeffs(cm, x, 0, k, PLANE_TYPE_Y_WITH_DC,
                              ta + (k & 1),
                              tl + (k >> 1), TX_4X4, 16);
    }
  }
  *distortion += thisdistortion;
  *labelyrate += thisrate;

  *distortion >>= 2;
  return RDCOST(x->rdmult, x->rddiv, *labelyrate, *distortion);
}

typedef struct {
  int_mv *ref_mv, *second_ref_mv;
  int_mv mvp;

  int64_t segment_rd;
  int r;
  int64_t d;
  int segment_yrate;
  MB_PREDICTION_MODE modes[4];
  int_mv mvs[4], second_mvs[4];
  int eobs[4];
  int mvthresh;
} BEST_SEG_INFO;

static INLINE int mv_check_bounds(MACROBLOCK *x, int_mv *mv) {
  int r = 0;
  r |= (mv->as_mv.row >> 3) < x->mv_row_min;
  r |= (mv->as_mv.row >> 3) > x->mv_row_max;
  r |= (mv->as_mv.col >> 3) < x->mv_col_min;
  r |= (mv->as_mv.col >> 3) > x->mv_col_max;
  return r;
}

static enum BlockSize get_block_size(int bw, int bh) {
  if (bw == 4 && bh == 4)
    return BLOCK_4X4;

  if (bw == 4 && bh == 8)
    return BLOCK_4X8;

  if (bw == 8 && bh == 4)
    return BLOCK_8X4;

  if (bw == 8 && bh == 8)
    return BLOCK_8X8;

  if (bw == 8 && bh == 16)
    return BLOCK_8X16;

  if (bw == 16 && bh == 8)
    return BLOCK_16X8;

  if (bw == 16 && bh == 16)
    return BLOCK_16X16;

  if (bw == 32 && bh == 32)
    return BLOCK_32X32;

  if (bw == 32 && bh == 16)
    return BLOCK_32X16;

  if (bw == 16 && bh == 32)
    return BLOCK_16X32;

  if (bw == 64 && bh == 32)
    return BLOCK_64X32;

  if (bw == 32 && bh == 64)
    return BLOCK_32X64;

  if (bw == 64 && bh == 64)
    return BLOCK_64X64;

  assert(0);
  return -1;
}

static INLINE void mi_buf_shift(MACROBLOCK *x, int i) {
  MB_MODE_INFO *mbmi = &x->e_mbd.mode_info_context->mbmi;
  x->plane[0].src.buf =
      raster_block_offset_uint8(&x->e_mbd, BLOCK_SIZE_SB8X8, 0, i,
                                x->plane[0].src.buf,
                                x->plane[0].src.stride);
  assert(((intptr_t)x->e_mbd.plane[0].pre[0].buf & 0x7) == 0);
  x->e_mbd.plane[0].pre[0].buf =
      raster_block_offset_uint8(&x->e_mbd, BLOCK_SIZE_SB8X8, 0, i,
                                x->e_mbd.plane[0].pre[0].buf,
                                x->e_mbd.plane[0].pre[0].stride);
  if (mbmi->ref_frame[1])
    x->e_mbd.plane[0].pre[1].buf =
        raster_block_offset_uint8(&x->e_mbd, BLOCK_SIZE_SB8X8, 0, i,
                                  x->e_mbd.plane[0].pre[1].buf,
                                  x->e_mbd.plane[0].pre[1].stride);
}

static INLINE void mi_buf_restore(MACROBLOCK *x, struct buf_2d orig_src,
                                  struct buf_2d orig_pre[2]) {
  MB_MODE_INFO *mbmi = &x->e_mbd.mode_info_context->mbmi;
  x->plane[0].src = orig_src;
  x->e_mbd.plane[0].pre[0] = orig_pre[0];
  if (mbmi->ref_frame[1])
    x->e_mbd.plane[0].pre[1] = orig_pre[1];
}

static void rd_check_segment_txsize(VP9_COMP *cpi, MACROBLOCK *x,
                                    BEST_SEG_INFO *bsi,
                                    int_mv seg_mvs[4][MAX_REF_FRAMES],
                                    int mi_row, int mi_col) {
  int i, j, br = 0, rate = 0, sbr = 0, idx, idy;
  int64_t bd = 0, sbd = 0;
  MB_PREDICTION_MODE this_mode;
  MB_MODE_INFO * mbmi = &x->e_mbd.mode_info_context->mbmi;
  const int label_count = 4;
  int64_t this_segment_rd = 0, other_segment_rd;
  int label_mv_thresh;
  int segmentyrate = 0;
  int best_eobs[4] = { 0 };
  BLOCK_SIZE_TYPE bsize = mbmi->sb_type;
  int bwl = b_width_log2(bsize), bw = 1 << bwl;
  int bhl = b_height_log2(bsize), bh = 1 << bhl;
  vp9_variance_fn_ptr_t *v_fn_ptr;
  ENTROPY_CONTEXT t_above[4], t_left[4];
  ENTROPY_CONTEXT t_above_b[4], t_left_b[4];

  vpx_memcpy(t_above, x->e_mbd.plane[0].above_context, sizeof(t_above));
  vpx_memcpy(t_left, x->e_mbd.plane[0].left_context, sizeof(t_left));

  v_fn_ptr = &cpi->fn_ptr[get_block_size(4 << bwl, 4 << bhl)];

  // 64 makes this threshold really big effectively
  // making it so that we very rarely check mvs on
  // segments.   setting this to 1 would make mv thresh
  // roughly equal to what it is for macroblocks
  label_mv_thresh = 1 * bsi->mvthresh / label_count;

  // Segmentation method overheads
  other_segment_rd = this_segment_rd;

  for (idy = 0; idy < 2; idy += bh) {
    for (idx = 0; idx < 2; idx += bw) {
      // TODO(jingning,rbultje): rewrite the rate-distortion optimization
      // loop for 4x4/4x8/8x4 block coding. to be replaced with new rd loop
      int_mv mode_mv[MB_MODE_COUNT], second_mode_mv[MB_MODE_COUNT];
      int_mv frame_mv[MB_MODE_COUNT][MAX_REF_FRAMES];
      int64_t best_label_rd = INT64_MAX, best_other_rd = INT64_MAX;
      MB_PREDICTION_MODE mode_selected = ZEROMV;
      int bestlabelyrate = 0;
      i = idy * 2 + idx;

      frame_mv[ZEROMV][mbmi->ref_frame[0]].as_int = 0;
      frame_mv[ZEROMV][mbmi->ref_frame[1]].as_int = 0;
      vp9_append_sub8x8_mvs_for_idx(&cpi->common, &x->e_mbd,
                                    &frame_mv[NEARESTMV][mbmi->ref_frame[0]],
                                    &frame_mv[NEARMV][mbmi->ref_frame[0]],
                                    i, 0);
      if (mbmi->ref_frame[1] > 0)
        vp9_append_sub8x8_mvs_for_idx(&cpi->common, &x->e_mbd,
                                   &frame_mv[NEARESTMV][mbmi->ref_frame[1]],
                                   &frame_mv[NEARMV][mbmi->ref_frame[1]],
                                   i, 1);

      // search for the best motion vector on this segment
      for (this_mode = NEARESTMV; this_mode <= NEWMV; ++this_mode) {
        int64_t this_rd;
        int64_t distortion;
        int labelyrate;
        ENTROPY_CONTEXT t_above_s[4], t_left_s[4];
        const struct buf_2d orig_src = x->plane[0].src;
        struct buf_2d orig_pre[2];

        vpx_memcpy(orig_pre, x->e_mbd.plane[0].pre, sizeof(orig_pre));

        vpx_memcpy(t_above_s, t_above, sizeof(t_above_s));
        vpx_memcpy(t_left_s, t_left, sizeof(t_left_s));

        // motion search for newmv (single predictor case only)
        if (mbmi->ref_frame[1] <= 0 && this_mode == NEWMV) {
          int step_param = 0;
          int further_steps;
          int thissme, bestsme = INT_MAX;
          int sadpb = x->sadperbit4;
          int_mv mvp_full;

          /* Is the best so far sufficiently good that we cant justify doing
           * and new motion search. */
          if (best_label_rd < label_mv_thresh)
            break;

          if (cpi->compressor_speed) {
            // use previous block's result as next block's MV predictor.
            if (i > 0) {
              bsi->mvp.as_int =
              x->e_mbd.mode_info_context->bmi[i - 1].as_mv[0].as_int;
              if (i == 2)
                bsi->mvp.as_int =
                x->e_mbd.mode_info_context->bmi[i - 2].as_mv[0].as_int;
              step_param = 2;
            }
          }

          further_steps = (MAX_MVSEARCH_STEPS - 1) - step_param;

          mvp_full.as_mv.row = bsi->mvp.as_mv.row >> 3;
          mvp_full.as_mv.col = bsi->mvp.as_mv.col >> 3;

          // adjust src pointer for this block
          mi_buf_shift(x, i);
          bestsme = vp9_full_pixel_diamond(cpi, x, &mvp_full, step_param,
                                           sadpb, further_steps, 0, v_fn_ptr,
                                           bsi->ref_mv, &mode_mv[NEWMV]);

          // Should we do a full search (best quality only)
          if (cpi->compressor_speed == 0) {
            /* Check if mvp_full is within the range. */
            clamp_mv(&mvp_full, x->mv_col_min, x->mv_col_max,
                     x->mv_row_min, x->mv_row_max);

            thissme = cpi->full_search_sad(x, &mvp_full,
                                           sadpb, 16, v_fn_ptr,
                                           x->nmvjointcost, x->mvcost,
                                           bsi->ref_mv, i);

            if (thissme < bestsme) {
              bestsme = thissme;
              mode_mv[NEWMV].as_int =
                  x->e_mbd.mode_info_context->bmi[i].as_mv[0].as_int;
            } else {
              /* The full search result is actually worse so re-instate the
               * previous best vector */
              x->e_mbd.mode_info_context->bmi[i].as_mv[0].as_int =
                  mode_mv[NEWMV].as_int;
            }
          }

          if (bestsme < INT_MAX) {
            int distortion;
            unsigned int sse;
            cpi->find_fractional_mv_step(x, &mode_mv[NEWMV],
                                         bsi->ref_mv, x->errorperbit, v_fn_ptr,
                                         x->nmvjointcost, x->mvcost,
                                         &distortion, &sse);

            // safe motion search result for use in compound prediction
            seg_mvs[i][mbmi->ref_frame[0]].as_int = mode_mv[NEWMV].as_int;
          }

          // restore src pointers
          mi_buf_restore(x, orig_src, orig_pre);
        } else if (mbmi->ref_frame[1] > 0 && this_mode == NEWMV) {
          if (seg_mvs[i][mbmi->ref_frame[1]].as_int == INVALID_MV ||
              seg_mvs[i][mbmi->ref_frame[0]].as_int == INVALID_MV)
            continue;

          // adjust src pointers
          mi_buf_shift(x, i);
          if (cpi->sf.comp_inter_joint_search_thresh < bsize) {
            int rate_mv;
            joint_motion_search(cpi, x, bsize, frame_mv[this_mode],
                                mi_row, mi_col, seg_mvs[i],
                                &rate_mv);
            seg_mvs[i][mbmi->ref_frame[0]].as_int =
                frame_mv[this_mode][mbmi->ref_frame[0]].as_int;
            seg_mvs[i][mbmi->ref_frame[1]].as_int =
                frame_mv[this_mode][mbmi->ref_frame[1]].as_int;
          }
          // restore src pointers
          mi_buf_restore(x, orig_src, orig_pre);
        }

        rate = labels2mode(x, i, this_mode, &mode_mv[this_mode],
                           &second_mode_mv[this_mode], frame_mv, seg_mvs[i],
                           bsi->ref_mv, bsi->second_ref_mv, x->nmvjointcost,
                           x->mvcost, cpi);

        // Trap vectors that reach beyond the UMV borders
        if (((mode_mv[this_mode].as_mv.row >> 3) < x->mv_row_min) ||
            ((mode_mv[this_mode].as_mv.row >> 3) > x->mv_row_max) ||
            ((mode_mv[this_mode].as_mv.col >> 3) < x->mv_col_min) ||
            ((mode_mv[this_mode].as_mv.col >> 3) > x->mv_col_max)) {
          continue;
        }
        if (mbmi->ref_frame[1] > 0 &&
            mv_check_bounds(x, &second_mode_mv[this_mode]))
          continue;

        this_rd = encode_inter_mb_segment(&cpi->common,
                                          x, i, &labelyrate,
                                          &distortion, t_above_s, t_left_s);
        this_rd += RDCOST(x->rdmult, x->rddiv, rate, 0);
        rate += labelyrate;

        if (this_rd < best_label_rd) {
          sbr = rate;
          sbd = distortion;
          bestlabelyrate = labelyrate;
          mode_selected = this_mode;
          best_label_rd = this_rd;
          best_eobs[i] = x->e_mbd.plane[0].eobs[i];
          vpx_memcpy(t_above_b, t_above_s, sizeof(t_above_s));
          vpx_memcpy(t_left_b, t_left_s, sizeof(t_left_s));
        }
      } /*for each 4x4 mode*/

      vpx_memcpy(t_above, t_above_b, sizeof(t_above));
      vpx_memcpy(t_left, t_left_b, sizeof(t_left));

      labels2mode(x, i, mode_selected, &mode_mv[mode_selected],
                  &second_mode_mv[mode_selected], frame_mv, seg_mvs[i],
                  bsi->ref_mv, bsi->second_ref_mv, x->nmvjointcost,
                  x->mvcost, cpi);

      br += sbr;
      bd += sbd;
      segmentyrate += bestlabelyrate;
      this_segment_rd += best_label_rd;
      other_segment_rd += best_other_rd;

      for (j = 1; j < bh; ++j)
        vpx_memcpy(&x->partition_info->bmi[i + j * 2],
                   &x->partition_info->bmi[i],
                   sizeof(x->partition_info->bmi[i]));
      for (j = 1; j < bw; ++j)
        vpx_memcpy(&x->partition_info->bmi[i + j],
                   &x->partition_info->bmi[i],
                   sizeof(x->partition_info->bmi[i]));
    }
  } /* for each label */

  if (this_segment_rd < bsi->segment_rd) {
    bsi->r = br;
    bsi->d = bd;
    bsi->segment_yrate = segmentyrate;
    bsi->segment_rd = this_segment_rd;

    // store everything needed to come back to this!!
    for (i = 0; i < 4; i++) {
      bsi->mvs[i].as_mv = x->partition_info->bmi[i].mv.as_mv;
      if (mbmi->ref_frame[1] > 0)
        bsi->second_mvs[i].as_mv = x->partition_info->bmi[i].second_mv.as_mv;
      bsi->modes[i] = x->partition_info->bmi[i].mode;
      bsi->eobs[i] = best_eobs[i];
    }
  }
}

static int rd_pick_best_mbsegmentation(VP9_COMP *cpi, MACROBLOCK *x,
                                       int_mv *best_ref_mv,
                                       int_mv *second_best_ref_mv,
                                       int64_t best_rd,
                                       int *returntotrate,
                                       int *returnyrate,
                                       int64_t *returndistortion,
                                       int *skippable, int mvthresh,
                                       int_mv seg_mvs[4][MAX_REF_FRAMES],
                                       int mi_row, int mi_col) {
  int i;
  BEST_SEG_INFO bsi;
  MB_MODE_INFO * mbmi = &x->e_mbd.mode_info_context->mbmi;

  vpx_memset(&bsi, 0, sizeof(bsi));

  bsi.segment_rd = best_rd;
  bsi.ref_mv = best_ref_mv;
  bsi.second_ref_mv = second_best_ref_mv;
  bsi.mvp.as_int = best_ref_mv->as_int;
  bsi.mvthresh = mvthresh;

  for (i = 0; i < 4; i++)
    bsi.modes[i] = ZEROMV;

  rd_check_segment_txsize(cpi, x, &bsi, seg_mvs, mi_row, mi_col);

  /* set it to the best */
  for (i = 0; i < 4; i++) {
    x->e_mbd.mode_info_context->bmi[i].as_mv[0].as_int = bsi.mvs[i].as_int;
    if (mbmi->ref_frame[1] > 0)
      x->e_mbd.mode_info_context->bmi[i].as_mv[1].as_int =
      bsi.second_mvs[i].as_int;
    x->e_mbd.plane[0].eobs[i] = bsi.eobs[i];
  }

  /* save partitions */
  x->partition_info->count = 4;

  for (i = 0; i < x->partition_info->count; i++) {
    x->partition_info->bmi[i].mode = bsi.modes[i];
    x->partition_info->bmi[i].mv.as_mv = bsi.mvs[i].as_mv;
    if (mbmi->ref_frame[1] > 0)
      x->partition_info->bmi[i].second_mv.as_mv = bsi.second_mvs[i].as_mv;
  }
  /*
   * used to set mbmi->mv.as_int
   */
  x->partition_info->bmi[3].mv.as_int = bsi.mvs[3].as_int;
  if (mbmi->ref_frame[1] > 0)
    x->partition_info->bmi[3].second_mv.as_int = bsi.second_mvs[3].as_int;

  *returntotrate = bsi.r;
  *returndistortion = bsi.d;
  *returnyrate = bsi.segment_yrate;
  *skippable = vp9_sby_is_skippable(&x->e_mbd, BLOCK_SIZE_SB8X8);
  mbmi->mode = bsi.modes[3];

  return (int)(bsi.segment_rd);
}

static void mv_pred(VP9_COMP *cpi, MACROBLOCK *x,
                    uint8_t *ref_y_buffer, int ref_y_stride,
                    int ref_frame, enum BlockSize block_size ) {
  MACROBLOCKD *xd = &x->e_mbd;
  MB_MODE_INFO *mbmi = &xd->mode_info_context->mbmi;
  int_mv this_mv;
  int i;
  int zero_seen = 0;
  int best_index = 0;
  int best_sad = INT_MAX;
  int this_sad = INT_MAX;

  uint8_t *src_y_ptr = x->plane[0].src.buf;
  uint8_t *ref_y_ptr;
  int row_offset, col_offset;

  // Get the sad for each candidate reference mv
  for (i = 0; i < MAX_MV_REF_CANDIDATES; i++) {
    this_mv.as_int = mbmi->ref_mvs[ref_frame][i].as_int;

    // The list is at an end if we see 0 for a second time.
    if (!this_mv.as_int && zero_seen)
      break;
    zero_seen = zero_seen || !this_mv.as_int;

    row_offset = this_mv.as_mv.row >> 3;
    col_offset = this_mv.as_mv.col >> 3;
    ref_y_ptr = ref_y_buffer + (ref_y_stride * row_offset) + col_offset;

    // Find sad for current vector.
    this_sad = cpi->fn_ptr[block_size].sdf(src_y_ptr, x->plane[0].src.stride,
                                           ref_y_ptr, ref_y_stride,
                                           0x7fffffff);

    // Note if it is the best so far.
    if (this_sad < best_sad) {
      best_sad = this_sad;
      best_index = i;
    }
  }

  // Note the index of the mv that worked best in the reference list.
  x->mv_best_ref_index[ref_frame] = best_index;
}

static void estimate_ref_frame_costs(VP9_COMP *cpi, int segment_id,
                                     unsigned int *ref_costs_single,
                                     unsigned int *ref_costs_comp,
                                     vp9_prob *comp_mode_p) {
  VP9_COMMON *const cm = &cpi->common;
  MACROBLOCKD *const xd = &cpi->mb.e_mbd;
  int seg_ref_active = vp9_segfeature_active(xd, segment_id,
                                             SEG_LVL_REF_FRAME);
  if (seg_ref_active) {
    vpx_memset(ref_costs_single, 0, MAX_REF_FRAMES * sizeof(*ref_costs_single));
    vpx_memset(ref_costs_comp,   0, MAX_REF_FRAMES * sizeof(*ref_costs_comp));
    *comp_mode_p = 128;
  } else {
    vp9_prob intra_inter_p = vp9_get_pred_prob(cm, xd, PRED_INTRA_INTER);
    vp9_prob comp_inter_p = 128;

    if (cm->comp_pred_mode == HYBRID_PREDICTION) {
      comp_inter_p = vp9_get_pred_prob(cm, xd, PRED_COMP_INTER_INTER);
      *comp_mode_p = comp_inter_p;
    } else {
      *comp_mode_p = 128;
    }

    ref_costs_single[INTRA_FRAME] = vp9_cost_bit(intra_inter_p, 0);

    if (cm->comp_pred_mode != COMP_PREDICTION_ONLY) {
      vp9_prob ref_single_p1 = vp9_get_pred_prob(cm, xd, PRED_SINGLE_REF_P1);
      vp9_prob ref_single_p2 = vp9_get_pred_prob(cm, xd, PRED_SINGLE_REF_P2);
      unsigned int base_cost = vp9_cost_bit(intra_inter_p, 1);

      if (cm->comp_pred_mode == HYBRID_PREDICTION)
        base_cost += vp9_cost_bit(comp_inter_p, 0);

      ref_costs_single[LAST_FRAME] = ref_costs_single[GOLDEN_FRAME] =
          ref_costs_single[ALTREF_FRAME] = base_cost;
      ref_costs_single[LAST_FRAME]   += vp9_cost_bit(ref_single_p1, 0);
      ref_costs_single[GOLDEN_FRAME] += vp9_cost_bit(ref_single_p1, 1);
      ref_costs_single[ALTREF_FRAME] += vp9_cost_bit(ref_single_p1, 1);
      ref_costs_single[GOLDEN_FRAME] += vp9_cost_bit(ref_single_p2, 0);
      ref_costs_single[ALTREF_FRAME] += vp9_cost_bit(ref_single_p2, 1);
    } else {
      ref_costs_single[LAST_FRAME]   = 512;
      ref_costs_single[GOLDEN_FRAME] = 512;
      ref_costs_single[ALTREF_FRAME] = 512;
    }
    if (cm->comp_pred_mode != SINGLE_PREDICTION_ONLY) {
      vp9_prob ref_comp_p = vp9_get_pred_prob(cm, xd, PRED_COMP_REF_P);
      unsigned int base_cost = vp9_cost_bit(intra_inter_p, 1);

      if (cm->comp_pred_mode == HYBRID_PREDICTION)
        base_cost += vp9_cost_bit(comp_inter_p, 1);

      ref_costs_comp[LAST_FRAME]   = base_cost + vp9_cost_bit(ref_comp_p, 0);
      ref_costs_comp[GOLDEN_FRAME] = base_cost + vp9_cost_bit(ref_comp_p, 1);
    } else {
      ref_costs_comp[LAST_FRAME]   = 512;
      ref_costs_comp[GOLDEN_FRAME] = 512;
    }
  }
}

static void store_coding_context(MACROBLOCK *x, PICK_MODE_CONTEXT *ctx,
                                 int mode_index,
                                 PARTITION_INFO *partition,
                                 int_mv *ref_mv,
                                 int_mv *second_ref_mv,
                                 int64_t comp_pred_diff[NB_PREDICTION_TYPES],
                                 int64_t txfm_size_diff[NB_TXFM_MODES]) {
  MACROBLOCKD *const xd = &x->e_mbd;

  // Take a snapshot of the coding context so it can be
  // restored if we decide to encode this way
  ctx->skip = x->skip;
  ctx->best_mode_index = mode_index;
  ctx->mic = *xd->mode_info_context;

  if (partition)
    ctx->partition_info = *partition;

  ctx->best_ref_mv.as_int = ref_mv->as_int;
  ctx->second_best_ref_mv.as_int = second_ref_mv->as_int;

  ctx->single_pred_diff = (int)comp_pred_diff[SINGLE_PREDICTION_ONLY];
  ctx->comp_pred_diff   = (int)comp_pred_diff[COMP_PREDICTION_ONLY];
  ctx->hybrid_pred_diff = (int)comp_pred_diff[HYBRID_PREDICTION];

  memcpy(ctx->txfm_rd_diff, txfm_size_diff, sizeof(ctx->txfm_rd_diff));
}

static void setup_pred_block(const MACROBLOCKD *xd,
                             struct buf_2d dst[MAX_MB_PLANE],
                             const YV12_BUFFER_CONFIG *src,
                             int mi_row, int mi_col,
                             const struct scale_factors *scale,
                             const struct scale_factors *scale_uv) {
  int i;

  dst[0].buf = src->y_buffer;
  dst[0].stride = src->y_stride;
  dst[1].buf = src->u_buffer;
  dst[2].buf = src->v_buffer;
  dst[1].stride = dst[2].stride = src->uv_stride;
#if CONFIG_ALPHA
  dst[3].buf = src->alpha_buffer;
  dst[3].stride = src->alpha_stride;
#endif

  // TODO(jkoleszar): Make scale factors per-plane data
  for (i = 0; i < MAX_MB_PLANE; i++) {
    setup_pred_plane(dst + i, dst[i].buf, dst[i].stride, mi_row, mi_col,
                     i ? scale_uv : scale,
                     xd->plane[i].subsampling_x, xd->plane[i].subsampling_y);
  }
}

static void setup_buffer_inter(VP9_COMP *cpi, MACROBLOCK *x,
                               int idx, MV_REFERENCE_FRAME frame_type,
                               enum BlockSize block_size,
                               int mi_row, int mi_col,
                               int_mv frame_nearest_mv[MAX_REF_FRAMES],
                               int_mv frame_near_mv[MAX_REF_FRAMES],
                               struct buf_2d yv12_mb[4][MAX_MB_PLANE],
                               struct scale_factors scale[MAX_REF_FRAMES]) {
  VP9_COMMON *cm = &cpi->common;
  YV12_BUFFER_CONFIG *yv12 = &cm->yv12_fb[cpi->common.ref_frame_map[idx]];
  MACROBLOCKD *const xd = &x->e_mbd;
  MB_MODE_INFO *const mbmi = &xd->mode_info_context->mbmi;

  // set up scaling factors
  scale[frame_type] = cpi->common.active_ref_scale[frame_type - 1];

  scale[frame_type].x_offset_q4 =
      ROUND_POWER_OF_TWO(mi_col * MI_SIZE * scale[frame_type].x_scale_fp,
       VP9_REF_SCALE_SHIFT) & 0xf;
  scale[frame_type].y_offset_q4 =
      ROUND_POWER_OF_TWO(mi_row * MI_SIZE * scale[frame_type].y_scale_fp,
       VP9_REF_SCALE_SHIFT) & 0xf;

  // TODO(jkoleszar): Is the UV buffer ever used here? If so, need to make this
  // use the UV scaling factors.
  setup_pred_block(xd, yv12_mb[frame_type], yv12, mi_row, mi_col,
                   &scale[frame_type], &scale[frame_type]);

  // Gets an initial list of candidate vectors from neighbours and orders them
  vp9_find_mv_refs(&cpi->common, xd, xd->mode_info_context,
                   xd->prev_mode_info_context,
                   frame_type,
                   mbmi->ref_mvs[frame_type],
                   cpi->common.ref_frame_sign_bias);

  // Candidate refinement carried out at encoder and decoder
  vp9_find_best_ref_mvs(xd,
                        mbmi->ref_mvs[frame_type],
                        &frame_nearest_mv[frame_type],
                        &frame_near_mv[frame_type]);

  // Further refinement that is encode side only to test the top few candidates
  // in full and choose the best as the centre point for subsequent searches.
  // The current implementation doesn't support scaling.
  if (scale[frame_type].x_scale_fp == (1 << VP9_REF_SCALE_SHIFT) &&
      scale[frame_type].y_scale_fp == (1 << VP9_REF_SCALE_SHIFT))
    mv_pred(cpi, x, yv12_mb[frame_type][0].buf, yv12->y_stride,
            frame_type, block_size);
}

static YV12_BUFFER_CONFIG *get_scaled_ref_frame(VP9_COMP *cpi, int ref_frame) {
  YV12_BUFFER_CONFIG *scaled_ref_frame = NULL;
  int fb = get_ref_frame_idx(cpi, ref_frame);
  if (cpi->scaled_ref_idx[fb] != cpi->common.ref_frame_map[fb])
    scaled_ref_frame = &cpi->common.yv12_fb[cpi->scaled_ref_idx[fb]];
  return scaled_ref_frame;
}

static double linear_interpolate(double x, int ntab, double step,
                                 const double *tab) {
  double y = x / step;
  int d = (int) y;
  double a = y - d;
  if (d >= ntab - 1)
    return tab[ntab - 1];
  else
    return tab[d] * (1 - a) + tab[d + 1] * a;
}

static double model_rate_norm(double x) {
  // Normalized rate
  // This function models the rate for a Laplacian source
  // source with given variance when quantized with a uniform quantizer
  // with given stepsize. The closed form expressions are in:
  // Hang and Chen, "Source Model for transform video coder and its
  // application - Part I: Fundamental Theory", IEEE Trans. Circ.
  // Sys. for Video Tech., April 1997.
  static const double rate_tab_step = 0.125;
  static const double rate_tab[] = {
    256.0000, 4.944453, 3.949276, 3.371593,
    2.965771, 2.654550, 2.403348, 2.193612,
    2.014208, 1.857921, 1.719813, 1.596364,
    1.484979, 1.383702, 1.291025, 1.205767,
    1.126990, 1.053937, 0.985991, 0.922644,
    0.863472, 0.808114, 0.756265, 0.707661,
    0.662070, 0.619287, 0.579129, 0.541431,
    0.506043, 0.472828, 0.441656, 0.412411,
    0.384980, 0.359260, 0.335152, 0.312563,
    0.291407, 0.271600, 0.253064, 0.235723,
    0.219508, 0.204351, 0.190189, 0.176961,
    0.164611, 0.153083, 0.142329, 0.132298,
    0.122945, 0.114228, 0.106106, 0.098541,
    0.091496, 0.084937, 0.078833, 0.073154,
    0.067872, 0.062959, 0.058392, 0.054147,
    0.050202, 0.046537, 0.043133, 0.039971,
    0.037036, 0.034312, 0.031783, 0.029436,
    0.027259, 0.025240, 0.023367, 0.021631,
    0.020021, 0.018528, 0.017145, 0.015863,
    0.014676, 0.013575, 0.012556, 0.011612,
    0.010738, 0.009929, 0.009180, 0.008487,
    0.007845, 0.007251, 0.006701, 0.006193,
    0.005722, 0.005287, 0.004884, 0.004512,
    0.004168, 0.003850, 0.003556, 0.003284,
    0.003032, 0.002800, 0.002585, 0.002386,
    0.002203, 0.002034, 0.001877, 0.001732,
    0.001599, 0.001476, 0.001362, 0.001256,
    0.001159, 0.001069, 0.000987, 0.000910,
    0.000840, 0.000774, 0.000714, 0.000659,
    0.000608, 0.000560, 0.000517, 0.000476,
    0.000439, 0.000405, 0.000373, 0.000344,
    0.000317, 0.000292, 0.000270, 0.000248,
    0.000229, 0.000211, 0.000195, 0.000179,
    0.000165, 0.000152, 0.000140, 0.000129,
    0.000119, 0.000110, 0.000101, 0.000093,
    0.000086, 0.000079, 0.000073, 0.000067,
    0.000062, 0.000057, 0.000052, 0.000048,
    0.000044, 0.000041, 0.000038, 0.000035,
    0.000032, 0.000029, 0.000027, 0.000025,
    0.000023, 0.000021, 0.000019, 0.000018,
    0.000016, 0.000015, 0.000014, 0.000013,
    0.000012, 0.000011, 0.000010, 0.000009,
    0.000008, 0.000008, 0.000007, 0.000007,
    0.000006, 0.000006, 0.000005, 0.000005,
    0.000004, 0.000004, 0.000004, 0.000003,
    0.000003, 0.000003, 0.000003, 0.000002,
    0.000002, 0.000002, 0.000002, 0.000002,
    0.000002, 0.000001, 0.000001, 0.000001,
    0.000001, 0.000001, 0.000001, 0.000001,
    0.000001, 0.000001, 0.000001, 0.000001,
    0.000001, 0.000001, 0.000000, 0.000000,
  };
  const int rate_tab_num = sizeof(rate_tab)/sizeof(rate_tab[0]);
  assert(x >= 0.0);
  return linear_interpolate(x, rate_tab_num, rate_tab_step, rate_tab);
}

static double model_dist_norm(double x) {
  // Normalized distortion
  // This function models the normalized distortion for a Laplacian source
  // source with given variance when quantized with a uniform quantizer
  // with given stepsize. The closed form expression is:
  // Dn(x) = 1 - 1/sqrt(2) * x / sinh(x/sqrt(2))
  // where x = qpstep / sqrt(variance)
  // Note the actual distortion is Dn * variance.
  static const double dist_tab_step = 0.25;
  static const double dist_tab[] = {
    0.000000, 0.005189, 0.020533, 0.045381,
    0.078716, 0.119246, 0.165508, 0.215979,
    0.269166, 0.323686, 0.378318, 0.432034,
    0.484006, 0.533607, 0.580389, 0.624063,
    0.664475, 0.701581, 0.735418, 0.766092,
    0.793751, 0.818575, 0.840761, 0.860515,
    0.878045, 0.893554, 0.907238, 0.919281,
    0.929857, 0.939124, 0.947229, 0.954306,
    0.960475, 0.965845, 0.970512, 0.974563,
    0.978076, 0.981118, 0.983750, 0.986024,
    0.987989, 0.989683, 0.991144, 0.992402,
    0.993485, 0.994417, 0.995218, 0.995905,
    0.996496, 0.997002, 0.997437, 0.997809,
    0.998128, 0.998401, 0.998635, 0.998835,
    0.999006, 0.999152, 0.999277, 0.999384,
    0.999475, 0.999553, 0.999619, 0.999676,
    0.999724, 0.999765, 0.999800, 0.999830,
    0.999855, 0.999877, 0.999895, 0.999911,
    0.999924, 0.999936, 0.999945, 0.999954,
    0.999961, 0.999967, 0.999972, 0.999976,
    0.999980, 0.999983, 0.999985, 0.999988,
    0.999989, 0.999991, 0.999992, 0.999994,
    0.999995, 0.999995, 0.999996, 0.999997,
    0.999997, 0.999998, 0.999998, 0.999998,
    0.999999, 0.999999, 0.999999, 0.999999,
    0.999999, 0.999999, 0.999999, 1.000000,
  };
  const int dist_tab_num = sizeof(dist_tab)/sizeof(dist_tab[0]);
  assert(x >= 0.0);
  return linear_interpolate(x, dist_tab_num, dist_tab_step, dist_tab);
}

static void model_rd_from_var_lapndz(int var, int n, int qstep,
                                     int *rate, int64_t *dist) {
  // This function models the rate and distortion for a Laplacian
  // source with given variance when quantized with a uniform quantizer
  // with given stepsize. The closed form expression is:
  // Rn(x) = H(sqrt(r)) + sqrt(r)*[1 + H(r)/(1 - r)],
  // where r = exp(-sqrt(2) * x) and x = qpstep / sqrt(variance)
  vp9_clear_system_state();
  if (var == 0 || n == 0) {
    *rate = 0;
    *dist = 0;
  } else {
    double D, R;
    double s2 = (double) var / n;
    double x = qstep / sqrt(s2);
    // TODO(debargha): Make the modeling functions take (qstep^2 / s2)
    // as argument rather than qstep / sqrt(s2) to obviate the need for
    // the sqrt() operation.
    D = model_dist_norm(x);
    R = model_rate_norm(x);
    if (R < 0) {
      R = 0;
      D = var;
    }
    *rate = (n * R * 256 + 0.5);
    *dist = (n * D * s2 + 0.5);
  }
  vp9_clear_system_state();
}

static enum BlockSize get_plane_block_size(BLOCK_SIZE_TYPE bsize,
                                           struct macroblockd_plane *pd) {
  return get_block_size(plane_block_width(bsize, pd),
                        plane_block_height(bsize, pd));
}

static void model_rd_for_sb(VP9_COMP *cpi, BLOCK_SIZE_TYPE bsize,
                            MACROBLOCK *x, MACROBLOCKD *xd,
                            int *out_rate_sum, int64_t *out_dist_sum) {
  // Note our transform coeffs are 8 times an orthogonal transform.
  // Hence quantizer step is also 8 times. To get effective quantizer
  // we need to divide by 8 before sending to modeling function.
  unsigned int sse;
  int i, rate_sum = 0;
  int64_t dist_sum = 0;

  for (i = 0; i < MAX_MB_PLANE; ++i) {
    struct macroblock_plane *const p = &x->plane[i];
    struct macroblockd_plane *const pd = &xd->plane[i];

    // TODO(dkovalev) the same code in get_plane_block_size
    const int bw = plane_block_width(bsize, pd);
    const int bh = plane_block_height(bsize, pd);
    const enum BlockSize bs = get_block_size(bw, bh);
    int rate;
    int64_t dist;
    cpi->fn_ptr[bs].vf(p->src.buf, p->src.stride,
                       pd->dst.buf, pd->dst.stride, &sse);
    model_rd_from_var_lapndz(sse, bw * bh, pd->dequant[1] >> 3, &rate, &dist);

    rate_sum += rate;
    dist_sum += dist;
  }

  *out_rate_sum = rate_sum;
  *out_dist_sum = dist_sum << 4;
}

static INLINE int get_switchable_rate(VP9_COMMON *cm, MACROBLOCK *x) {
  MACROBLOCKD *xd = &x->e_mbd;
  MB_MODE_INFO *const mbmi = &xd->mode_info_context->mbmi;

  const int c = vp9_get_pred_context(cm, xd, PRED_SWITCHABLE_INTERP);
  const int m = vp9_switchable_interp_map[mbmi->interp_filter];
  return SWITCHABLE_INTERP_RATE_FACTOR * x->switchable_interp_costs[c][m];
}

static void single_motion_search(VP9_COMP *cpi, MACROBLOCK *x,
                                 BLOCK_SIZE_TYPE bsize,
                                 int mi_row, int mi_col,
                                 int_mv *tmp_mv, int *rate_mv) {
  MACROBLOCKD *xd = &x->e_mbd;
  MB_MODE_INFO *mbmi = &xd->mode_info_context->mbmi;
  struct buf_2d backup_yv12[MAX_MB_PLANE] = {{0}};
  int bestsme = INT_MAX;
  int further_steps, step_param = cpi->sf.first_step;
  int sadpb = x->sadperbit16;
  int_mv mvp_full;
  int ref = mbmi->ref_frame[0];
  int_mv ref_mv = mbmi->ref_mvs[ref][0];
  int sr = 0;
  const enum BlockSize block_size = get_plane_block_size(bsize, &xd->plane[0]);

  int tmp_col_min = x->mv_col_min;
  int tmp_col_max = x->mv_col_max;
  int tmp_row_min = x->mv_row_min;
  int tmp_row_max = x->mv_row_max;

  YV12_BUFFER_CONFIG *scaled_ref_frame = get_scaled_ref_frame(cpi, ref);

  if (scaled_ref_frame) {
    int i;
    // Swap out the reference frame for a version that's been scaled to
    // match the resolution of the current frame, allowing the existing
    // motion search code to be used without additional modifications.
    for (i = 0; i < MAX_MB_PLANE; i++)
      backup_yv12[i] = xd->plane[i].pre[0];

    setup_pre_planes(xd, scaled_ref_frame, NULL, mi_row, mi_col,
                     NULL, NULL);
  }

  vp9_clamp_mv_min_max(x, &ref_mv);

  sr = vp9_init_search_range(cpi->common.width, cpi->common.height);

  // mvp_full.as_int = ref_mv[0].as_int;
  mvp_full.as_int =
      mbmi->ref_mvs[ref][x->mv_best_ref_index[ref]].as_int;

  mvp_full.as_mv.col >>= 3;
  mvp_full.as_mv.row >>= 3;

  // adjust search range according to sr from mv prediction
  step_param = MAX(step_param, sr);

  // Further step/diamond searches as necessary
  further_steps = (cpi->sf.max_step_search_steps - 1) - step_param;

  bestsme = vp9_full_pixel_diamond(cpi, x, &mvp_full, step_param,
                                   sadpb, further_steps, 1,
                                   &cpi->fn_ptr[block_size],
                                   &ref_mv, tmp_mv);

  x->mv_col_min = tmp_col_min;
  x->mv_col_max = tmp_col_max;
  x->mv_row_min = tmp_row_min;
  x->mv_row_max = tmp_row_max;

  if (bestsme < INT_MAX) {
    int dis; /* TODO: use dis in distortion calculation later. */
    unsigned int sse;
    cpi->find_fractional_mv_step(x, tmp_mv, &ref_mv,
                                 x->errorperbit,
                                 &cpi->fn_ptr[block_size],
                                 x->nmvjointcost, x->mvcost,
                                 &dis, &sse);
  }
  *rate_mv = vp9_mv_bit_cost(tmp_mv, &ref_mv,
                             x->nmvjointcost, x->mvcost,
                             96, xd->allow_high_precision_mv);
  if (scaled_ref_frame) {
    int i;
    for (i = 0; i < MAX_MB_PLANE; i++)
      xd->plane[i].pre[0] = backup_yv12[i];
  }
}

static void joint_motion_search(VP9_COMP *cpi, MACROBLOCK *x,
                                BLOCK_SIZE_TYPE bsize,
                                int_mv *frame_mv,
                                int mi_row, int mi_col,
                                int_mv single_newmv[MAX_REF_FRAMES],
                                int *rate_mv) {
  int pw = 4 << b_width_log2(bsize), ph = 4 << b_height_log2(bsize);
  MACROBLOCKD *xd = &x->e_mbd;
  MB_MODE_INFO *mbmi = &xd->mode_info_context->mbmi;
  int refs[2] = { mbmi->ref_frame[0],
    (mbmi->ref_frame[1] < 0 ? 0 : mbmi->ref_frame[1]) };
  int_mv ref_mv[2];
  const enum BlockSize block_size = get_plane_block_size(bsize, &xd->plane[0]);
  int ite;
  // Prediction buffer from second frame.
  uint8_t *second_pred = vpx_memalign(16, pw * ph * sizeof(uint8_t));

  // Do joint motion search in compound mode to get more accurate mv.
  struct buf_2d backup_yv12[MAX_MB_PLANE] = {{0}};
  struct buf_2d backup_second_yv12[MAX_MB_PLANE] = {{0}};
  struct buf_2d scaled_first_yv12;
  int last_besterr[2] = {INT_MAX, INT_MAX};
  YV12_BUFFER_CONFIG *scaled_ref_frame[2] = {NULL, NULL};
  scaled_ref_frame[0] = get_scaled_ref_frame(cpi, mbmi->ref_frame[0]);
  scaled_ref_frame[1] = get_scaled_ref_frame(cpi, mbmi->ref_frame[1]);

  ref_mv[0] = mbmi->ref_mvs[refs[0]][0];
  ref_mv[1] = mbmi->ref_mvs[refs[1]][0];

  if (scaled_ref_frame[0]) {
    int i;
    // Swap out the reference frame for a version that's been scaled to
    // match the resolution of the current frame, allowing the existing
    // motion search code to be used without additional modifications.
    for (i = 0; i < MAX_MB_PLANE; i++)
      backup_yv12[i] = xd->plane[i].pre[0];
    setup_pre_planes(xd, scaled_ref_frame[0], NULL, mi_row, mi_col,
                     NULL, NULL);
  }

  if (scaled_ref_frame[1]) {
    int i;
    for (i = 0; i < MAX_MB_PLANE; i++)
      backup_second_yv12[i] = xd->plane[i].pre[1];

    setup_pre_planes(xd, scaled_ref_frame[1], NULL, mi_row, mi_col,
                     NULL, NULL);
  }

  xd->scale_factor[0].set_scaled_offsets(&xd->scale_factor[0],
                                         mi_row, mi_col);
  xd->scale_factor[1].set_scaled_offsets(&xd->scale_factor[1],
                                         mi_row, mi_col);
  scaled_first_yv12 = xd->plane[0].pre[0];

  // Initialize mv using single prediction mode result.
  frame_mv[refs[0]].as_int = single_newmv[refs[0]].as_int;
  frame_mv[refs[1]].as_int = single_newmv[refs[1]].as_int;

  // Allow joint search multiple times iteratively for each ref frame
  // and break out the search loop if it couldn't find better mv.
  for (ite = 0; ite < 4; ite++) {
    struct buf_2d ref_yv12[2];
    int bestsme = INT_MAX;
    int sadpb = x->sadperbit16;
    int_mv tmp_mv;
    int search_range = 3;

    int tmp_col_min = x->mv_col_min;
    int tmp_col_max = x->mv_col_max;
    int tmp_row_min = x->mv_row_min;
    int tmp_row_max = x->mv_row_max;
    int id = ite % 2;

    // Initialized here because of compiler problem in Visual Studio.
    ref_yv12[0] = xd->plane[0].pre[0];
    ref_yv12[1] = xd->plane[0].pre[1];

    // Get pred block from second frame.
    vp9_build_inter_predictor(ref_yv12[!id].buf,
                              ref_yv12[!id].stride,
                              second_pred, pw,
                              &frame_mv[refs[!id]],
                              &xd->scale_factor[!id],
                              pw, ph, 0,
                              &xd->subpix);

    // Compound motion search on first ref frame.
    if (id)
      xd->plane[0].pre[0] = ref_yv12[id];
    vp9_clamp_mv_min_max(x, &ref_mv[id]);

    // Use mv result from single mode as mvp.
    tmp_mv.as_int = frame_mv[refs[id]].as_int;

    tmp_mv.as_mv.col >>= 3;
    tmp_mv.as_mv.row >>= 3;

    // Small-range full-pixel motion search
    bestsme = vp9_refining_search_8p_c(x, &tmp_mv, sadpb,
                                       search_range,
                                       &cpi->fn_ptr[block_size],
                                       x->nmvjointcost, x->mvcost,
                                       &ref_mv[id], second_pred,
                                       pw, ph);

    x->mv_col_min = tmp_col_min;
    x->mv_col_max = tmp_col_max;
    x->mv_row_min = tmp_row_min;
    x->mv_row_max = tmp_row_max;

    if (bestsme < INT_MAX) {
      int dis; /* TODO: use dis in distortion calculation later. */
      unsigned int sse;

      bestsme = vp9_find_best_sub_pixel_comp(x, &tmp_mv,
                                             &ref_mv[id],
                                             x->errorperbit,
                                             &cpi->fn_ptr[block_size],
                                             x->nmvjointcost, x->mvcost,
                                             &dis, &sse, second_pred,
                                             pw, ph);
    }

    if (id)
      xd->plane[0].pre[0] = scaled_first_yv12;

    if (bestsme < last_besterr[id]) {
      frame_mv[refs[id]].as_int = tmp_mv.as_int;
      last_besterr[id] = bestsme;
    } else {
      break;
    }
  }

  // restore the predictor
  if (scaled_ref_frame[0]) {
    int i;
    for (i = 0; i < MAX_MB_PLANE; i++)
      xd->plane[i].pre[0] = backup_yv12[i];
  }

  if (scaled_ref_frame[1]) {
    int i;
    for (i = 0; i < MAX_MB_PLANE; i++)
      xd->plane[i].pre[1] = backup_second_yv12[i];
  }
  *rate_mv  = vp9_mv_bit_cost(&frame_mv[refs[0]],
                              &mbmi->ref_mvs[refs[0]][0],
                              x->nmvjointcost, x->mvcost, 96,
                              x->e_mbd.allow_high_precision_mv);
  *rate_mv += vp9_mv_bit_cost(&frame_mv[refs[1]],
                              &mbmi->ref_mvs[refs[1]][0],
                              x->nmvjointcost, x->mvcost, 96,
                              x->e_mbd.allow_high_precision_mv);

  vpx_free(second_pred);
}

static int64_t handle_inter_mode(VP9_COMP *cpi, MACROBLOCK *x,
                                 BLOCK_SIZE_TYPE bsize,
                                 int64_t txfm_cache[],
                                 int *rate2, int64_t *distortion,
                                 int *skippable,
                                 int *rate_y, int64_t *distortion_y,
                                 int *rate_uv, int64_t *distortion_uv,
                                 int *mode_excluded, int *disable_skip,
                                 INTERPOLATIONFILTERTYPE *best_filter,
                                 int_mv *frame_mv,
                                 int mi_row, int mi_col,
                                 int_mv single_newmv[MAX_REF_FRAMES]) {
  const int bw = 1 << mi_width_log2(bsize), bh = 1 << mi_height_log2(bsize);

  VP9_COMMON *cm = &cpi->common;
  MACROBLOCKD *xd = &x->e_mbd;
  const enum BlockSize block_size = get_plane_block_size(bsize, &xd->plane[0]);
  const enum BlockSize uv_block_size = get_plane_block_size(bsize,
                                                            &xd->plane[1]);
  MB_MODE_INFO *mbmi = &xd->mode_info_context->mbmi;
  const int is_comp_pred = (mbmi->ref_frame[1] > 0);
  const int num_refs = is_comp_pred ? 2 : 1;
  const int this_mode = mbmi->mode;
  int i;
  int refs[2] = { mbmi->ref_frame[0],
    (mbmi->ref_frame[1] < 0 ? 0 : mbmi->ref_frame[1]) };
  int_mv cur_mv[2];
  int64_t this_rd = 0;
  unsigned char tmp_buf[MAX_MB_PLANE][64 * 64];
  int pred_exists = 0;
  int interpolating_intpel_seen = 0;
  int intpel_mv;
  int64_t rd, best_rd = INT64_MAX;

  switch (this_mode) {
    int rate_mv;
    case NEWMV:
      if (is_comp_pred) {
        // Initialize mv using single prediction mode result.
        frame_mv[refs[0]].as_int = single_newmv[refs[0]].as_int;
        frame_mv[refs[1]].as_int = single_newmv[refs[1]].as_int;

        if (cpi->sf.comp_inter_joint_search_thresh < bsize) {
          joint_motion_search(cpi, x, bsize, frame_mv,
                              mi_row, mi_col, single_newmv, &rate_mv);
        } else {
          rate_mv  = vp9_mv_bit_cost(&frame_mv[refs[0]],
                                     &mbmi->ref_mvs[refs[0]][0],
                                     x->nmvjointcost, x->mvcost, 96,
                                     x->e_mbd.allow_high_precision_mv);
          rate_mv += vp9_mv_bit_cost(&frame_mv[refs[1]],
                                     &mbmi->ref_mvs[refs[1]][0],
                                     x->nmvjointcost, x->mvcost, 96,
                                     x->e_mbd.allow_high_precision_mv);
        }
        if (frame_mv[refs[0]].as_int == INVALID_MV ||
            frame_mv[refs[1]].as_int == INVALID_MV)
          return INT64_MAX;
        *rate2 += rate_mv;

      } else {
        int_mv tmp_mv;
        single_motion_search(cpi, x, bsize, mi_row, mi_col,
                             &tmp_mv, &rate_mv);
        *rate2 += rate_mv;
        frame_mv[refs[0]].as_int =
            xd->mode_info_context->bmi[0].as_mv[0].as_int = tmp_mv.as_int;
        single_newmv[refs[0]].as_int = tmp_mv.as_int;
      }
      break;
    case NEARMV:
    case NEARESTMV:
    case ZEROMV:
    default:
      break;
  }
  for (i = 0; i < num_refs; ++i) {
    cur_mv[i] = frame_mv[refs[i]];
    // Clip "next_nearest" so that it does not extend to far out of image
    if (this_mode == NEWMV)
      assert(!clamp_mv2(&cur_mv[i], xd));
    else
      clamp_mv2(&cur_mv[i], xd);

    if (mv_check_bounds(x, &cur_mv[i]))
      return INT64_MAX;
    mbmi->mv[i].as_int = cur_mv[i].as_int;
  }

  /* We don't include the cost of the second reference here, because there
   * are only three options: Last/Golden, ARF/Last or Golden/ARF, or in other
   * words if you present them in that order, the second one is always known
   * if the first is known */
  *rate2 += vp9_cost_mv_ref(cpi, this_mode,
                            mbmi->mb_mode_context[mbmi->ref_frame[0]]);

  pred_exists = 0;
  interpolating_intpel_seen = 0;
  // Are all MVs integer pel for Y and UV
  intpel_mv = (mbmi->mv[0].as_mv.row & 15) == 0 &&
      (mbmi->mv[0].as_mv.col & 15) == 0;
  if (is_comp_pred)
    intpel_mv &= (mbmi->mv[1].as_mv.row & 15) == 0 &&
        (mbmi->mv[1].as_mv.col & 15) == 0;
  // Search for best switchable filter by checking the variance of
  // pred error irrespective of whether the filter will be used
  if (cpi->sf.use_8tap_always) {
    *best_filter = EIGHTTAP;
  } else {
    int i, newbest;
    int tmp_rate_sum = 0;
    int64_t tmp_dist_sum = 0;
    for (i = 0; i < VP9_SWITCHABLE_FILTERS; ++i) {
      int rs = 0;
      const INTERPOLATIONFILTERTYPE filter = vp9_switchable_interp[i];
      const int is_intpel_interp = intpel_mv &&
          vp9_is_interpolating_filter[filter];
      mbmi->interp_filter = filter;
      vp9_setup_interp_filters(xd, mbmi->interp_filter, cm);

      if (cm->mcomp_filter_type == SWITCHABLE)
        rs = get_switchable_rate(cm, x);

      if (interpolating_intpel_seen && is_intpel_interp) {
        rd = RDCOST(x->rdmult, x->rddiv, rs + tmp_rate_sum, tmp_dist_sum);
      } else {
        int rate_sum = 0;
        int64_t dist_sum = 0;
        vp9_build_inter_predictors_sb(xd, mi_row, mi_col, bsize);
        model_rd_for_sb(cpi, bsize, x, xd, &rate_sum, &dist_sum);
        rd = RDCOST(x->rdmult, x->rddiv, rs + rate_sum, dist_sum);
        if (!interpolating_intpel_seen && is_intpel_interp) {
          tmp_rate_sum = rate_sum;
          tmp_dist_sum = dist_sum;
        }
      }
      newbest = i == 0 || rd < best_rd;

      if (newbest) {
        best_rd = rd;
        *best_filter = mbmi->interp_filter;
      }

      if ((cm->mcomp_filter_type == SWITCHABLE && newbest) ||
          (cm->mcomp_filter_type != SWITCHABLE &&
           cm->mcomp_filter_type == mbmi->interp_filter)) {
        int p;

        for (p = 0; p < MAX_MB_PLANE; p++) {
          const int y = (MI_SIZE * bh) >> xd->plane[p].subsampling_y;
          const int x = (MI_SIZE * bw) >> xd->plane[p].subsampling_x;
          int i;

          for (i = 0; i < y; i++)
            vpx_memcpy(&tmp_buf[p][64 * i],
                       xd->plane[p].dst.buf + i * xd->plane[p].dst.stride, x);
        }
        pred_exists = 1;
      }
      interpolating_intpel_seen |= is_intpel_interp;
    }
  }

  // Set the appripriate filter
  mbmi->interp_filter = cm->mcomp_filter_type != SWITCHABLE ?
      cm->mcomp_filter_type : *best_filter;
  vp9_setup_interp_filters(xd, mbmi->interp_filter, cm);


  if (pred_exists) {
    int p;

    for (p = 0; p < MAX_MB_PLANE; p++) {
      const int y = (MI_SIZE * bh) >> xd->plane[p].subsampling_y;
      const int x = (MI_SIZE * bw) >> xd->plane[p].subsampling_x;
      int i;

      for (i = 0; i < y; i++)
        vpx_memcpy(xd->plane[p].dst.buf + i * xd->plane[p].dst.stride,
                   &tmp_buf[p][64 * i], x);
    }
  } else {
    // Handles the special case when a filter that is not in the
    // switchable list (ex. bilinear, 6-tap) is indicated at the frame level
    vp9_build_inter_predictors_sb(xd, mi_row, mi_col, bsize);
  }

  if (cpi->common.mcomp_filter_type == SWITCHABLE)
    *rate2 += get_switchable_rate(cm, x);

  if (cpi->active_map_enabled && x->active_ptr[0] == 0)
    x->skip = 1;
  else if (x->encode_breakout) {
    unsigned int var, sse;
    int threshold = (xd->plane[0].dequant[1]
                     * xd->plane[0].dequant[1] >> 4);

    if (threshold < x->encode_breakout)
      threshold = x->encode_breakout;

    var = cpi->fn_ptr[block_size].vf(x->plane[0].src.buf,
                                     x->plane[0].src.stride,
                                     xd->plane[0].dst.buf,
                                     xd->plane[0].dst.stride,
                                     &sse);

    if ((int)sse < threshold) {
      unsigned int q2dc = xd->plane[0].dequant[0];
      /* If there is no codeable 2nd order dc
         or a very small uniform pixel change change */
      if ((sse - var < q2dc * q2dc >> 4) ||
          (sse / 2 > var && sse - var < 64)) {
        // Check u and v to make sure skip is ok
        int sse2;
        unsigned int sse2u, sse2v;
        var = cpi->fn_ptr[uv_block_size].vf(x->plane[1].src.buf,
                                            x->plane[1].src.stride,
                                            xd->plane[1].dst.buf,
                                            xd->plane[1].dst.stride, &sse2u);
        var = cpi->fn_ptr[uv_block_size].vf(x->plane[2].src.buf,
                                            x->plane[1].src.stride,
                                            xd->plane[2].dst.buf,
                                            xd->plane[1].dst.stride, &sse2v);
        sse2 = sse2u + sse2v;

        if (sse2 * 2 < threshold) {
          x->skip = 1;
          *distortion = sse + sse2;
          *rate2 = 500;

          /* for best_yrd calculation */
          *rate_uv = 0;
          *distortion_uv = sse2;

          *disable_skip = 1;
          this_rd = RDCOST(x->rdmult, x->rddiv, *rate2, *distortion);
        }
      }
    }
  }

  if (!x->skip) {
    int skippable_y, skippable_uv;

    // Y cost and distortion
    super_block_yrd(cpi, x, rate_y, distortion_y, &skippable_y,
                    bsize, txfm_cache);

    *rate2 += *rate_y;
    *distortion += *distortion_y;

    super_block_uvrd(cm, x, rate_uv, distortion_uv,
                     &skippable_uv, bsize);

    *rate2 += *rate_uv;
    *distortion += *distortion_uv;
    *skippable = skippable_y && skippable_uv;
  }

  if (!(*mode_excluded)) {
    if (is_comp_pred) {
      *mode_excluded = (cpi->common.comp_pred_mode == SINGLE_PREDICTION_ONLY);
    } else {
      *mode_excluded = (cpi->common.comp_pred_mode == COMP_PREDICTION_ONLY);
    }
  }

  return this_rd;  // if 0, this will be re-calculated by caller
}

void vp9_rd_pick_intra_mode_sb(VP9_COMP *cpi, MACROBLOCK *x,
                               int *returnrate, int64_t *returndist,
                               BLOCK_SIZE_TYPE bsize,
                               PICK_MODE_CONTEXT *ctx) {
  VP9_COMMON *cm = &cpi->common;
  MACROBLOCKD *xd = &x->e_mbd;
  int rate_y = 0, rate_uv = 0;
  int rate_y_tokenonly = 0, rate_uv_tokenonly = 0;
  int64_t dist_y = 0, dist_uv = 0;
  int y_skip = 0, uv_skip = 0;
  int64_t txfm_cache[NB_TXFM_MODES], err;
  MB_PREDICTION_MODE mode;
  TX_SIZE txfm_size;
  int rate4x4_y, rate4x4_y_tokenonly;
  int64_t dist4x4_y;
  int64_t err4x4 = INT64_MAX;
  int i;

  vpx_memset(&txfm_cache,0,sizeof(txfm_cache));
  ctx->skip = 0;
  xd->mode_info_context->mbmi.mode = DC_PRED;
  xd->mode_info_context->mbmi.ref_frame[0] = INTRA_FRAME;
  err = rd_pick_intra_sby_mode(cpi, x, &rate_y, &rate_y_tokenonly,
                               &dist_y, &y_skip, bsize, txfm_cache);
  mode = xd->mode_info_context->mbmi.mode;
  txfm_size = xd->mode_info_context->mbmi.txfm_size;
  rd_pick_intra_sbuv_mode(cpi, x, &rate_uv, &rate_uv_tokenonly,
                          &dist_uv, &uv_skip,
                          (bsize < BLOCK_SIZE_SB8X8) ? BLOCK_SIZE_SB8X8 :
                                                       bsize);
  if (bsize < BLOCK_SIZE_SB8X8)
    err4x4 = rd_pick_intra4x4mby_modes(cpi, x, &rate4x4_y,
                                       &rate4x4_y_tokenonly,
                                       &dist4x4_y, err);

  if (y_skip && uv_skip) {
    *returnrate = rate_y + rate_uv - rate_y_tokenonly - rate_uv_tokenonly +
                  vp9_cost_bit(vp9_get_pred_prob(cm, xd, PRED_MBSKIP), 1);
    *returndist = dist_y + (dist_uv >> 2);
    memset(ctx->txfm_rd_diff, 0, sizeof(ctx->txfm_rd_diff));
    xd->mode_info_context->mbmi.mode = mode;
    xd->mode_info_context->mbmi.txfm_size = txfm_size;
  } else if (bsize < BLOCK_SIZE_SB8X8 && err4x4 < err) {
    *returnrate = rate4x4_y + rate_uv +
        vp9_cost_bit(vp9_get_pred_prob(cm, xd, PRED_MBSKIP), 0);
    *returndist = dist4x4_y + (dist_uv >> 2);
    vpx_memset(ctx->txfm_rd_diff, 0, sizeof(ctx->txfm_rd_diff));
    xd->mode_info_context->mbmi.txfm_size = TX_4X4;
  } else {
    *returnrate = rate_y + rate_uv +
        vp9_cost_bit(vp9_get_pred_prob(cm, xd, PRED_MBSKIP), 0);
    *returndist = dist_y + (dist_uv >> 2);
    for (i = 0; i < NB_TXFM_MODES; i++) {
      ctx->txfm_rd_diff[i] = txfm_cache[i] - txfm_cache[cm->txfm_mode];
    }
    xd->mode_info_context->mbmi.txfm_size = txfm_size;
    xd->mode_info_context->mbmi.mode = mode;
  }

  ctx->mic = *xd->mode_info_context;
}

int64_t vp9_rd_pick_inter_mode_sb(VP9_COMP *cpi, MACROBLOCK *x,
                                  int mi_row, int mi_col,
                                  int *returnrate,
                                  int64_t *returndistortion,
                                  BLOCK_SIZE_TYPE bsize,
                                  PICK_MODE_CONTEXT *ctx) {
  VP9_COMMON *cm = &cpi->common;
  MACROBLOCKD *xd = &x->e_mbd;
  MB_MODE_INFO *mbmi = &xd->mode_info_context->mbmi;
  const enum BlockSize block_size = get_plane_block_size(bsize, &xd->plane[0]);
  MB_PREDICTION_MODE this_mode;
  MB_PREDICTION_MODE best_mode = DC_PRED;
  MV_REFERENCE_FRAME ref_frame;
  unsigned char segment_id = xd->mode_info_context->mbmi.segment_id;
  int comp_pred, i;
  int_mv frame_mv[MB_MODE_COUNT][MAX_REF_FRAMES];
  struct buf_2d yv12_mb[4][MAX_MB_PLANE];
  int_mv single_newmv[MAX_REF_FRAMES];
  static const int flag_list[4] = { 0, VP9_LAST_FLAG, VP9_GOLD_FLAG,
                                    VP9_ALT_FLAG };
  int idx_list[4] = {0,
                     cpi->lst_fb_idx,
                     cpi->gld_fb_idx,
                     cpi->alt_fb_idx};
  int64_t best_rd = INT64_MAX;
  int64_t best_txfm_rd[NB_TXFM_MODES];
  int64_t best_txfm_diff[NB_TXFM_MODES];
  int64_t best_pred_diff[NB_PREDICTION_TYPES];
  int64_t best_pred_rd[NB_PREDICTION_TYPES];
  MB_MODE_INFO best_mbmode;
  int j;
  int mode_index, best_mode_index = 0;
  unsigned int ref_costs_single[MAX_REF_FRAMES], ref_costs_comp[MAX_REF_FRAMES];
  vp9_prob comp_mode_p;
  int64_t best_overall_rd = INT64_MAX;
  INTERPOLATIONFILTERTYPE best_filter = SWITCHABLE;
  INTERPOLATIONFILTERTYPE tmp_best_filter = SWITCHABLE;
  int rate_uv_intra[TX_SIZE_MAX_SB], rate_uv_tokenonly[TX_SIZE_MAX_SB];
  int64_t dist_uv[TX_SIZE_MAX_SB];
  int skip_uv[TX_SIZE_MAX_SB];
  MB_PREDICTION_MODE mode_uv[TX_SIZE_MAX_SB];
  struct scale_factors scale_factor[4];
  unsigned int ref_frame_mask = 0;
  unsigned int mode_mask = 0;
  int64_t mode_distortions[MB_MODE_COUNT] = {-1};
  int64_t frame_distortions[MAX_REF_FRAMES] = {-1};
  int intra_cost_penalty = 20 * vp9_dc_quant(cpi->common.base_qindex,
                                             cpi->common.y_dc_delta_q);
  int_mv seg_mvs[4][MAX_REF_FRAMES];
  union b_mode_info best_bmodes[4];
  PARTITION_INFO best_partition;
  int bwsl = b_width_log2(bsize);
  int bws = (1 << bwsl) / 4;  // mode_info step for subsize
  int bhsl = b_height_log2(bsize);
  int bhs = (1 << bhsl) / 4;  // mode_info step for subsize

  for (i = 0; i < 4; i++) {
    int j;

    for (j = 0; j < MAX_REF_FRAMES; j++)
      seg_mvs[i][j].as_int = INVALID_MV;
  }
  // Everywhere the flag is set the error is much higher than its neighbors.
  ctx->frames_with_high_error = 0;
  ctx->modes_with_high_error = 0;

  xd->mode_info_context->mbmi.segment_id = segment_id;
  estimate_ref_frame_costs(cpi, segment_id, ref_costs_single, ref_costs_comp,
                           &comp_mode_p);
  vpx_memset(&best_mbmode, 0, sizeof(best_mbmode));
  vpx_memset(&single_newmv, 0, sizeof(single_newmv));

  for (i = 0; i < NB_PREDICTION_TYPES; ++i)
    best_pred_rd[i] = INT64_MAX;
  for (i = 0; i < NB_TXFM_MODES; i++)
    best_txfm_rd[i] = INT64_MAX;

  // Create a mask set to 1 for each frame used by a smaller resolution.
  if (cpi->sf.use_avoid_tested_higherror) {
    switch (block_size) {
      case BLOCK_64X64:
        for (i = 0; i < 4; i++) {
          for (j = 0; j < 4; j++) {
            ref_frame_mask |= x->mb_context[i][j].frames_with_high_error;
            mode_mask |= x->mb_context[i][j].modes_with_high_error;
          }
        }
        for (i = 0; i < 4; i++) {
          ref_frame_mask |= x->sb32_context[i].frames_with_high_error;
          mode_mask |= x->sb32_context[i].modes_with_high_error;
        }
        break;
      case BLOCK_32X32:
        for (i = 0; i < 4; i++) {
          ref_frame_mask |=
              x->mb_context[xd->sb_index][i].frames_with_high_error;
          mode_mask |= x->mb_context[xd->sb_index][i].modes_with_high_error;
        }
        break;
      default:
        // Until we handle all block sizes set it to present;
        ref_frame_mask = 0;
        mode_mask = 0;
        break;
    }
    ref_frame_mask = ~ref_frame_mask;
    mode_mask = ~mode_mask;
  }

  for (ref_frame = LAST_FRAME; ref_frame <= ALTREF_FRAME; ref_frame++) {
    if (cpi->ref_frame_flags & flag_list[ref_frame]) {
      setup_buffer_inter(cpi, x, idx_list[ref_frame], ref_frame, block_size,
                         mi_row, mi_col, frame_mv[NEARESTMV], frame_mv[NEARMV],
                         yv12_mb, scale_factor);
    }
    frame_mv[NEWMV][ref_frame].as_int = INVALID_MV;
    frame_mv[ZEROMV][ref_frame].as_int = 0;
  }
  if (!cpi->sf.use_avoid_tested_higherror
      || (cpi->sf.use_avoid_tested_higherror
          && (ref_frame_mask & (1 << INTRA_FRAME)))) {
    mbmi->mode = DC_PRED;
    mbmi->ref_frame[0] = INTRA_FRAME;
    for (i = 0; i <= (bsize < BLOCK_SIZE_MB16X16 ? TX_4X4 :
                      (bsize < BLOCK_SIZE_SB32X32 ? TX_8X8 :
                       (bsize < BLOCK_SIZE_SB64X64 ? TX_16X16 : TX_32X32)));
         i++) {
      mbmi->txfm_size = i;
      rd_pick_intra_sbuv_mode(cpi, x, &rate_uv_intra[i], &rate_uv_tokenonly[i],
                              &dist_uv[i], &skip_uv[i],
                              (bsize < BLOCK_SIZE_SB8X8) ? BLOCK_SIZE_SB8X8 :
                                                           bsize);
      mode_uv[i] = mbmi->uv_mode;
    }
  }

  for (mode_index = 0; mode_index < MAX_MODES; ++mode_index) {
    int mode_excluded = 0;
    int64_t this_rd = INT64_MAX;
    int disable_skip = 0;
    int compmode_cost = 0;
    int rate2 = 0, rate_y = 0, rate_uv = 0;
    int64_t distortion2 = 0, distortion_y = 0, distortion_uv = 0;
    int skippable;
    int64_t txfm_cache[NB_TXFM_MODES];
    int i;

    for (i = 0; i < NB_TXFM_MODES; ++i)
      txfm_cache[i] = INT64_MAX;

    // Test best rd so far against threshold for trying this mode.
    if ((best_rd < ((cpi->rd_threshes[bsize][mode_index] *
                     cpi->rd_thresh_freq_fact[bsize][mode_index]) >> 4)) ||
        cpi->rd_threshes[bsize][mode_index] == INT_MAX)
      continue;

    // Do not allow compound prediction if the segment level reference
    // frame feature is in use as in this case there can only be one reference.
    if ((vp9_mode_order[mode_index].second_ref_frame > INTRA_FRAME) &&
         vp9_segfeature_active(xd, segment_id, SEG_LVL_REF_FRAME))
      continue;

    x->skip = 0;
    this_mode = vp9_mode_order[mode_index].mode;
    ref_frame = vp9_mode_order[mode_index].ref_frame;

    if (cpi->sf.use_avoid_tested_higherror && bsize >= BLOCK_SIZE_SB8X8) {
      if (!(ref_frame_mask & (1 << ref_frame))) {
        continue;
      }
      if (!(mode_mask & (1 << this_mode))) {
        continue;
      }
      if (vp9_mode_order[mode_index].second_ref_frame != NONE
          && !(ref_frame_mask
              & (1 << vp9_mode_order[mode_index].second_ref_frame))) {
        continue;
      }
    }

    mbmi->ref_frame[0] = ref_frame;
    mbmi->ref_frame[1] = vp9_mode_order[mode_index].second_ref_frame;

    if (!(ref_frame == INTRA_FRAME
        || (cpi->ref_frame_flags & flag_list[ref_frame]))) {
      continue;
    }
    if (!(mbmi->ref_frame[1] == NONE
        || (cpi->ref_frame_flags & flag_list[mbmi->ref_frame[1]]))) {
      continue;
    }

    // TODO(jingning, jkoleszar): scaling reference frame not supported for
    // SPLITMV.
    if (mbmi->ref_frame[0] > 0 &&
          (scale_factor[mbmi->ref_frame[0]].x_scale_fp !=
           (1 << VP9_REF_SCALE_SHIFT) ||
           scale_factor[mbmi->ref_frame[0]].y_scale_fp !=
           (1 << VP9_REF_SCALE_SHIFT)) &&
        this_mode == SPLITMV)
      continue;

    if (mbmi->ref_frame[1] > 0 &&
          (scale_factor[mbmi->ref_frame[1]].x_scale_fp !=
           (1 << VP9_REF_SCALE_SHIFT) ||
           scale_factor[mbmi->ref_frame[1]].y_scale_fp !=
           (1 << VP9_REF_SCALE_SHIFT)) &&
        this_mode == SPLITMV)
      continue;

    set_scale_factors(xd, mbmi->ref_frame[0], mbmi->ref_frame[1],
                      scale_factor);
    comp_pred = mbmi->ref_frame[1] > INTRA_FRAME;
    mbmi->mode = this_mode;
    mbmi->uv_mode = DC_PRED;

    // Evaluate all sub-pel filters irrespective of whether we can use
    // them for this frame.
    mbmi->interp_filter = cm->mcomp_filter_type;
    vp9_setup_interp_filters(xd, mbmi->interp_filter, &cpi->common);

    if (bsize >= BLOCK_SIZE_SB8X8 &&
        (this_mode == I4X4_PRED || this_mode == SPLITMV))
      continue;
    if (bsize < BLOCK_SIZE_SB8X8 &&
        !(this_mode == I4X4_PRED || this_mode == SPLITMV))
      continue;

    if (comp_pred) {
      if (!(cpi->ref_frame_flags & flag_list[mbmi->ref_frame[1]]))
        continue;
      set_scale_factors(xd, mbmi->ref_frame[0], mbmi->ref_frame[1],
                        scale_factor);

      mode_excluded =
          mode_excluded ?
              mode_excluded : cm->comp_pred_mode == SINGLE_PREDICTION_ONLY;
    } else {
      // mbmi->ref_frame[1] = vp9_mode_order[mode_index].ref_frame[1];
      if (ref_frame != INTRA_FRAME) {
        if (mbmi->ref_frame[1] != INTRA_FRAME)
          mode_excluded =
              mode_excluded ?
                  mode_excluded : cm->comp_pred_mode == COMP_PREDICTION_ONLY;
      }
    }

    // Select predictors
    for (i = 0; i < MAX_MB_PLANE; i++) {
      xd->plane[i].pre[0] = yv12_mb[ref_frame][i];
      if (comp_pred)
        xd->plane[i].pre[1] = yv12_mb[mbmi->ref_frame[1]][i];
    }

    // If the segment reference frame feature is enabled....
    // then do nothing if the current ref frame is not allowed..
    if (vp9_segfeature_active(xd, segment_id, SEG_LVL_REF_FRAME) &&
        vp9_get_segdata(xd, segment_id, SEG_LVL_REF_FRAME) != (int)ref_frame) {
      continue;
    // If the segment skip feature is enabled....
    // then do nothing if the current mode is not allowed..
    } else if (vp9_segfeature_active(xd, segment_id, SEG_LVL_SKIP) &&
               (this_mode != ZEROMV && ref_frame != INTRA_FRAME)) {
      continue;
    // Disable this drop out case if the ref frame
    // segment level feature is enabled for this segment. This is to
    // prevent the possibility that we end up unable to pick any mode.
    } else if (!vp9_segfeature_active(xd, segment_id, SEG_LVL_REF_FRAME)) {
      // Only consider ZEROMV/ALTREF_FRAME for alt ref frame,
      // unless ARNR filtering is enabled in which case we want
      // an unfiltered alternative
      if (cpi->is_src_frame_alt_ref && (cpi->oxcf.arnr_max_frames == 0)) {
        if (this_mode != ZEROMV || ref_frame != ALTREF_FRAME) {
          continue;
        }
      }
    }
    // TODO(JBB): This is to make up for the fact that we don't have sad
    // functions that work when the block size reads outside the umv.  We
    // should fix this either by making the motion search just work on
    // a representative block in the boundary ( first ) and then implement a
    // function that does sads when inside the border..
    if (((mi_row + bhs) > cm->mi_rows || (mi_col + bws) > cm->mi_cols) &&
        this_mode == NEWMV) {
      continue;
    }

    if (this_mode == I4X4_PRED) {
      int rate;

      mbmi->txfm_size = TX_4X4;
      rd_pick_intra4x4mby_modes(cpi, x, &rate, &rate_y,
                                &distortion_y, INT64_MAX);
      rate2 += rate;
      rate2 += intra_cost_penalty;
      distortion2 += distortion_y;

      rate2 += rate_uv_intra[TX_4X4];
      rate_uv = rate_uv_intra[TX_4X4];
      distortion2 += dist_uv[TX_4X4];
      distortion_uv = dist_uv[TX_4X4];
      mbmi->uv_mode = mode_uv[TX_4X4];
      txfm_cache[ONLY_4X4] = RDCOST(x->rdmult, x->rddiv, rate2, distortion2);
      for (i = 0; i < NB_TXFM_MODES; ++i)
        txfm_cache[i] = txfm_cache[ONLY_4X4];
    } else if (ref_frame == INTRA_FRAME) {
      TX_SIZE uv_tx;
      super_block_yrd(cpi, x, &rate_y, &distortion_y, &skippable,
                      bsize, txfm_cache);

      uv_tx = mbmi->txfm_size;
      if (bsize < BLOCK_SIZE_MB16X16 && uv_tx == TX_8X8)
        uv_tx = TX_4X4;
      if (bsize < BLOCK_SIZE_SB32X32 && uv_tx == TX_16X16)
        uv_tx = TX_8X8;
      else if (bsize < BLOCK_SIZE_SB64X64 && uv_tx == TX_32X32)
        uv_tx = TX_16X16;

      rate_uv = rate_uv_intra[uv_tx];
      distortion_uv = dist_uv[uv_tx];
      skippable = skippable && skip_uv[uv_tx];
      mbmi->uv_mode = mode_uv[uv_tx];

      rate2 = rate_y + x->mbmode_cost[mbmi->mode] + rate_uv;
      if (mbmi->mode != DC_PRED && mbmi->mode != TM_PRED)
        rate2 += intra_cost_penalty;
      distortion2 = distortion_y + distortion_uv;
    } else if (this_mode == SPLITMV) {
      const int is_comp_pred = mbmi->ref_frame[1] > 0;
      int rate;
      int64_t distortion;
      int64_t this_rd_thresh;
      int64_t tmp_rd, tmp_best_rd = INT64_MAX, tmp_best_rdu = INT64_MAX;
      int tmp_best_rate = INT_MAX, tmp_best_ratey = INT_MAX;
      int64_t tmp_best_distortion = INT_MAX;
      int tmp_best_skippable = 0;
      int switchable_filter_index;
      int_mv *second_ref = is_comp_pred ?
          &mbmi->ref_mvs[mbmi->ref_frame[1]][0] : NULL;
      union b_mode_info tmp_best_bmodes[16];
      MB_MODE_INFO tmp_best_mbmode;
      PARTITION_INFO tmp_best_partition;
      int pred_exists = 0;
      int uv_skippable;

      this_rd_thresh = (mbmi->ref_frame[0] == LAST_FRAME) ?
          cpi->rd_threshes[bsize][THR_NEWMV] :
          cpi->rd_threshes[bsize][THR_NEWA];
      this_rd_thresh = (mbmi->ref_frame[0] == GOLDEN_FRAME) ?
          cpi->rd_threshes[bsize][THR_NEWG] : this_rd_thresh;
      xd->mode_info_context->mbmi.txfm_size = TX_4X4;

      for (switchable_filter_index = 0;
           switchable_filter_index < VP9_SWITCHABLE_FILTERS;
           ++switchable_filter_index) {
        int newbest;
        mbmi->interp_filter =
        vp9_switchable_interp[switchable_filter_index];
        vp9_setup_interp_filters(xd, mbmi->interp_filter, &cpi->common);

        tmp_rd = rd_pick_best_mbsegmentation(cpi, x,
                     &mbmi->ref_mvs[mbmi->ref_frame[0]][0],
                     second_ref, INT64_MAX,
                     &rate, &rate_y, &distortion,
                     &skippable,
                     (int)this_rd_thresh, seg_mvs,
                     mi_row, mi_col);
        if (cpi->common.mcomp_filter_type == SWITCHABLE) {
          const int rs = get_switchable_rate(cm, x);
          tmp_rd += RDCOST(x->rdmult, x->rddiv, rs, 0);
        }
        newbest = (tmp_rd < tmp_best_rd);
        if (newbest) {
          tmp_best_filter = mbmi->interp_filter;
          tmp_best_rd = tmp_rd;
        }
        if ((newbest && cm->mcomp_filter_type == SWITCHABLE) ||
            (mbmi->interp_filter == cm->mcomp_filter_type &&
             cm->mcomp_filter_type != SWITCHABLE)) {
              tmp_best_rdu = tmp_rd;
              tmp_best_rate = rate;
              tmp_best_ratey = rate_y;
              tmp_best_distortion = distortion;
              tmp_best_skippable = skippable;
              tmp_best_mbmode = *mbmi;
              tmp_best_partition = *x->partition_info;
              for (i = 0; i < 4; i++)
                tmp_best_bmodes[i] = xd->mode_info_context->bmi[i];
              pred_exists = 1;
            }
      }  // switchable_filter_index loop

      mbmi->interp_filter = (cm->mcomp_filter_type == SWITCHABLE ?
                             tmp_best_filter : cm->mcomp_filter_type);
      vp9_setup_interp_filters(xd, mbmi->interp_filter, &cpi->common);
      if (!pred_exists) {
        // Handles the special case when a filter that is not in the
        // switchable list (bilinear, 6-tap) is indicated at the frame level
        tmp_rd = rd_pick_best_mbsegmentation(cpi, x,
                     &mbmi->ref_mvs[mbmi->ref_frame[0]][0],
                     second_ref, INT64_MAX,
                     &rate, &rate_y, &distortion,
                     &skippable,
                     (int)this_rd_thresh, seg_mvs,
                     mi_row, mi_col);
      } else {
        if (cpi->common.mcomp_filter_type == SWITCHABLE) {
          int rs = get_switchable_rate(cm, x);
          tmp_best_rdu -= RDCOST(x->rdmult, x->rddiv, rs, 0);
        }
        tmp_rd = tmp_best_rdu;
        rate = tmp_best_rate;
        rate_y = tmp_best_ratey;
        distortion = tmp_best_distortion;
        skippable = tmp_best_skippable;
        *mbmi = tmp_best_mbmode;
        *x->partition_info = tmp_best_partition;
        for (i = 0; i < 4; i++)
          xd->mode_info_context->bmi[i] = tmp_best_bmodes[i];
      }

      rate2 += rate;
      distortion2 += distortion;

      if (cpi->common.mcomp_filter_type == SWITCHABLE)
        rate2 += get_switchable_rate(cm, x);

      // If even the 'Y' rd value of split is higher than best so far
      // then dont bother looking at UV
      vp9_build_inter_predictors_sbuv(&x->e_mbd, mi_row, mi_col,
                                      BLOCK_SIZE_SB8X8);
      vp9_subtract_sbuv(x, BLOCK_SIZE_SB8X8);
      super_block_uvrd_for_txfm(cm, x, &rate_uv, &distortion_uv,
                                &uv_skippable, BLOCK_SIZE_SB8X8, TX_4X4);
      rate2 += rate_uv;
      distortion2 += distortion_uv;
      skippable = skippable && uv_skippable;

      txfm_cache[ONLY_4X4] = RDCOST(x->rdmult, x->rddiv, rate2, distortion2);
      for (i = 0; i < NB_TXFM_MODES; ++i)
        txfm_cache[i] = txfm_cache[ONLY_4X4];

      if (!mode_excluded) {
        if (is_comp_pred)
          mode_excluded = cpi->common.comp_pred_mode == SINGLE_PREDICTION_ONLY;
        else
          mode_excluded = cpi->common.comp_pred_mode == COMP_PREDICTION_ONLY;
      }

      compmode_cost = vp9_cost_bit(comp_mode_p, is_comp_pred);
    } else {
      compmode_cost = vp9_cost_bit(comp_mode_p,
                                   mbmi->ref_frame[1] > INTRA_FRAME);
      this_rd = handle_inter_mode(cpi, x, bsize,
                                  txfm_cache,
                                  &rate2, &distortion2, &skippable,
                                  &rate_y, &distortion_y,
                                  &rate_uv, &distortion_uv,
                                  &mode_excluded, &disable_skip,
                                  &tmp_best_filter, frame_mv[this_mode],
                                  mi_row, mi_col,
                                  single_newmv);
      if (this_rd == INT64_MAX)
        continue;
    }

    if (cpi->common.comp_pred_mode == HYBRID_PREDICTION) {
      rate2 += compmode_cost;
    }

    // Estimate the reference frame signaling cost and add it
    // to the rolling cost variable.
    if (mbmi->ref_frame[1] > INTRA_FRAME) {
      rate2 += ref_costs_comp[mbmi->ref_frame[0]];
    } else {
      rate2 += ref_costs_single[mbmi->ref_frame[0]];
    }

    if (!disable_skip) {
      // Test for the condition where skip block will be activated
      // because there are no non zero coefficients and make any
      // necessary adjustment for rate. Ignore if skip is coded at
      // segment level as the cost wont have been added in.
      int mb_skip_allowed;

      // Is Mb level skip allowed (i.e. not coded at segment level).
      mb_skip_allowed = !vp9_segfeature_active(xd, segment_id, SEG_LVL_SKIP);

      if (skippable && bsize >= BLOCK_SIZE_SB8X8) {
        // Back out the coefficient coding costs
        rate2 -= (rate_y + rate_uv);
        // for best_yrd calculation
        rate_uv = 0;

        if (mb_skip_allowed) {
          int prob_skip_cost;

          // Cost the skip mb case
          vp9_prob skip_prob =
            vp9_get_pred_prob(cm, xd, PRED_MBSKIP);

          if (skip_prob) {
            prob_skip_cost = vp9_cost_bit(skip_prob, 1);
            rate2 += prob_skip_cost;
          }
        }
      } else if (mb_skip_allowed) {
        // Add in the cost of the no skip flag.
        int prob_skip_cost = vp9_cost_bit(vp9_get_pred_prob(cm, xd,
                                                        PRED_MBSKIP), 0);
        rate2 += prob_skip_cost;
      }

      // Calculate the final RD estimate for this mode.
      this_rd = RDCOST(x->rdmult, x->rddiv, rate2, distortion2);
    }

#if 0
    // Keep record of best intra distortion
    if ((xd->mode_info_context->mbmi.ref_frame[0] == INTRA_FRAME) &&
        (this_rd < best_intra_rd)) {
      best_intra_rd = this_rd;
      *returnintra = distortion2;
    }
#endif

    if (!disable_skip && mbmi->ref_frame[0] == INTRA_FRAME)
      for (i = 0; i < NB_PREDICTION_TYPES; ++i)
        best_pred_rd[i] = MIN(best_pred_rd[i], this_rd);

    if (this_rd < best_overall_rd) {
      best_overall_rd = this_rd;
      best_filter = tmp_best_filter;
      best_mode = this_mode;
    }

    if (this_mode != I4X4_PRED && this_mode != SPLITMV) {
      // Store the respective mode distortions for later use.
      if (mode_distortions[this_mode] == -1
          || distortion2 < mode_distortions[this_mode]) {
        mode_distortions[this_mode] = distortion2;
      }
      if (frame_distortions[mbmi->ref_frame[0]] == -1
          || distortion2 < frame_distortions[mbmi->ref_frame[0]]) {
        frame_distortions[mbmi->ref_frame[0]] = distortion2;
      }
    }

    // Did this mode help.. i.e. is it the new best mode
    if (this_rd < best_rd || x->skip) {
      if (!mode_excluded) {
        // Note index of best mode so far
        best_mode_index = mode_index;

        if (ref_frame == INTRA_FRAME) {
          /* required for left and above block mv */
          mbmi->mv[0].as_int = 0;
        }

        *returnrate = rate2;
        *returndistortion = distortion2;
        best_rd = this_rd;
        best_mbmode = *mbmi;
        best_partition = *x->partition_info;

        if (this_mode == I4X4_PRED || this_mode == SPLITMV)
          for (i = 0; i < 4; i++)
            best_bmodes[i] = xd->mode_info_context->bmi[i];
      }
#if 0
      // Testing this mode gave rise to an improvement in best error score.
      // Lower threshold a bit for next time
      cpi->rd_thresh_mult[mode_index] =
          (cpi->rd_thresh_mult[mode_index] >= (MIN_THRESHMULT + 2)) ?
              cpi->rd_thresh_mult[mode_index] - 2 : MIN_THRESHMULT;
      cpi->rd_threshes[mode_index] =
          (cpi->rd_baseline_thresh[mode_index] >> 7)
              * cpi->rd_thresh_mult[mode_index];
#endif
    } else {
      // If the mode did not help improve the best error case then
      // raise the threshold for testing that mode next time around.
#if 0
      cpi->rd_thresh_mult[mode_index] += 4;

      if (cpi->rd_thresh_mult[mode_index] > MAX_THRESHMULT)
        cpi->rd_thresh_mult[mode_index] = MAX_THRESHMULT;

      cpi->rd_threshes[mode_index] =
          (cpi->rd_baseline_thresh[mode_index] >> 7)
              * cpi->rd_thresh_mult[mode_index];
#endif
    }

    /* keep record of best compound/single-only prediction */
    if (!disable_skip && mbmi->ref_frame[0] != INTRA_FRAME) {
      int single_rd, hybrid_rd, single_rate, hybrid_rate;

      if (cpi->common.comp_pred_mode == HYBRID_PREDICTION) {
        single_rate = rate2 - compmode_cost;
        hybrid_rate = rate2;
      } else {
        single_rate = rate2;
        hybrid_rate = rate2 + compmode_cost;
      }

      single_rd = RDCOST(x->rdmult, x->rddiv, single_rate, distortion2);
      hybrid_rd = RDCOST(x->rdmult, x->rddiv, hybrid_rate, distortion2);

      if (mbmi->ref_frame[1] <= INTRA_FRAME &&
          single_rd < best_pred_rd[SINGLE_PREDICTION_ONLY]) {
        best_pred_rd[SINGLE_PREDICTION_ONLY] = single_rd;
      } else if (mbmi->ref_frame[1] > INTRA_FRAME &&
                 single_rd < best_pred_rd[COMP_PREDICTION_ONLY]) {
        best_pred_rd[COMP_PREDICTION_ONLY] = single_rd;
      }
      if (hybrid_rd < best_pred_rd[HYBRID_PREDICTION])
        best_pred_rd[HYBRID_PREDICTION] = hybrid_rd;
    }

    /* keep record of best txfm size */
    if (bsize < BLOCK_SIZE_SB32X32) {
      if (bsize < BLOCK_SIZE_MB16X16) {
        if (this_mode == SPLITMV || this_mode == I4X4_PRED)
          txfm_cache[ALLOW_8X8] = txfm_cache[ONLY_4X4];
        txfm_cache[ALLOW_16X16] = txfm_cache[ALLOW_8X8];
      }
      txfm_cache[ALLOW_32X32] = txfm_cache[ALLOW_16X16];
    }
    if (!mode_excluded && this_rd != INT64_MAX) {
      for (i = 0; i < NB_TXFM_MODES; i++) {
        int64_t adj_rd = INT64_MAX;
        if (this_mode != I4X4_PRED) {
          adj_rd = this_rd + txfm_cache[i] - txfm_cache[cm->txfm_mode];
        } else {
          adj_rd = this_rd;
        }

        if (adj_rd < best_txfm_rd[i])
          best_txfm_rd[i] = adj_rd;
      }
    }

    if (x->skip && !mode_excluded)
      break;
  }
  // Flag all modes that have a distortion thats > 2x the best we found at
  // this level.
  for (mode_index = 0; mode_index < MB_MODE_COUNT; ++mode_index) {
    if (mode_index == NEARESTMV || mode_index == NEARMV || mode_index == NEWMV)
      continue;

    if (mode_distortions[mode_index] > 2 * *returndistortion) {
      ctx->modes_with_high_error |= (1 << mode_index);
    }
  }

  // Flag all ref frames that have a distortion thats > 2x the best we found at
  // this level.
  for (ref_frame = INTRA_FRAME; ref_frame <= ALTREF_FRAME; ref_frame++) {
    if (frame_distortions[ref_frame] > 2 * *returndistortion) {
      ctx->frames_with_high_error |= (1 << ref_frame);
    }
  }

  if (best_rd == INT64_MAX && bsize < BLOCK_SIZE_SB8X8) {
    *returnrate = INT_MAX;
    *returndistortion = INT_MAX;
    return best_rd;
  }

  assert((cm->mcomp_filter_type == SWITCHABLE) ||
         (cm->mcomp_filter_type == best_mbmode.interp_filter) ||
         (best_mbmode.ref_frame[0] == INTRA_FRAME));

  // Accumulate filter usage stats
  // TODO(agrange): Use RD criteria to select interpolation filter mode.
  if (is_inter_mode(best_mode))
    ++cpi->best_switchable_interp_count[vp9_switchable_interp_map[best_filter]];

  // Updating rd_thresh_freq_fact[] here means that the differnt
  // partition/block sizes are handled independently based on the best
  // choice for the current partition. It may well be better to keep a scaled
  // best rd so far value and update rd_thresh_freq_fact based on the mode/size
  // combination that wins out.
  if (cpi->sf.adpative_rd_thresh) {
    for (mode_index = 0; mode_index < MAX_MODES; ++mode_index) {
      if (mode_index == best_mode_index) {
        cpi->rd_thresh_freq_fact[bsize][mode_index] = BASE_RD_THRESH_FREQ_FACT;
      } else {
        cpi->rd_thresh_freq_fact[bsize][mode_index] += MAX_RD_THRESH_FREQ_INC;
        if (cpi->rd_thresh_freq_fact[bsize][mode_index] >
            (cpi->sf.adpative_rd_thresh * MAX_RD_THRESH_FREQ_FACT)) {
          cpi->rd_thresh_freq_fact[bsize][mode_index] =
            cpi->sf.adpative_rd_thresh * MAX_RD_THRESH_FREQ_FACT;
        }
      }
    }
  }

  // TODO(rbultje) integrate with RD trd_thresh_freq_facthresholding
#if 0
  // Reduce the activation RD thresholds for the best choice mode
  if ((cpi->rd_baseline_thresh[best_mode_index] > 0) &&
      (cpi->rd_baseline_thresh[best_mode_index] < (INT_MAX >> 2))) {
    int best_adjustment = (cpi->rd_thresh_mult[best_mode_index] >> 2);

    cpi->rd_thresh_mult[best_mode_index] =
      (cpi->rd_thresh_mult[best_mode_index] >= (MIN_THRESHMULT + best_adjustment)) ?
      cpi->rd_thresh_mult[best_mode_index] - best_adjustment : MIN_THRESHMULT;
    cpi->rd_threshes[best_mode_index] =
      (cpi->rd_baseline_thresh[best_mode_index] >> 7) * cpi->rd_thresh_mult[best_mode_index];
  }
#endif

  // This code forces Altref,0,0 and skip for the frame that overlays a
  // an alrtef unless Altref is filtered. However, this is unsafe if
  // segment level coding of ref frame is enabled for this segment.
  if (!vp9_segfeature_active(xd, segment_id, SEG_LVL_REF_FRAME) &&
      cpi->is_src_frame_alt_ref &&
      (cpi->oxcf.arnr_max_frames == 0) &&
      (best_mbmode.mode != ZEROMV || best_mbmode.ref_frame[0] != ALTREF_FRAME)
      && bsize >= BLOCK_SIZE_SB8X8) {
    mbmi->mode = ZEROMV;
    mbmi->ref_frame[0] = ALTREF_FRAME;
    mbmi->ref_frame[1] = NONE;
    mbmi->mv[0].as_int = 0;
    mbmi->uv_mode = DC_PRED;
    mbmi->mb_skip_coeff = 1;
    if (cm->txfm_mode == TX_MODE_SELECT) {
      if (bsize >= BLOCK_SIZE_SB32X32)
        mbmi->txfm_size = TX_32X32;
      else if (bsize >= BLOCK_SIZE_MB16X16)
        mbmi->txfm_size = TX_16X16;
      else
        mbmi->txfm_size = TX_8X8;
    }

    vpx_memset(best_txfm_diff, 0, sizeof(best_txfm_diff));
    vpx_memset(best_pred_diff, 0, sizeof(best_pred_diff));
    goto end;
  }

  // macroblock modes
  *mbmi = best_mbmode;
  if (best_mbmode.ref_frame[0] == INTRA_FRAME &&
      best_mbmode.sb_type < BLOCK_SIZE_SB8X8) {
    for (i = 0; i < 4; i++)
      xd->mode_info_context->bmi[i].as_mode = best_bmodes[i].as_mode;
  }

  if (best_mbmode.ref_frame[0] != INTRA_FRAME &&
      best_mbmode.sb_type < BLOCK_SIZE_SB8X8) {
    for (i = 0; i < 4; i++)
      xd->mode_info_context->bmi[i].as_mv[0].as_int =
          best_bmodes[i].as_mv[0].as_int;

    if (mbmi->ref_frame[1] > 0)
      for (i = 0; i < 4; i++)
        xd->mode_info_context->bmi[i].as_mv[1].as_int =
            best_bmodes[i].as_mv[1].as_int;

    *x->partition_info = best_partition;

    mbmi->mv[0].as_int = x->partition_info->bmi[3].mv.as_int;
    mbmi->mv[1].as_int = x->partition_info->bmi[3].second_mv.as_int;
  }

  for (i = 0; i < NB_PREDICTION_TYPES; ++i) {
    if (best_pred_rd[i] == INT64_MAX)
      best_pred_diff[i] = INT_MIN;
    else
      best_pred_diff[i] = best_rd - best_pred_rd[i];
  }

  if (!x->skip) {
    for (i = 0; i < NB_TXFM_MODES; i++) {
      if (best_txfm_rd[i] == INT64_MAX)
        best_txfm_diff[i] = 0;
      else
        best_txfm_diff[i] = best_rd - best_txfm_rd[i];
    }
  } else {
    vpx_memset(best_txfm_diff, 0, sizeof(best_txfm_diff));
  }

 end:
  set_scale_factors(xd, mbmi->ref_frame[0], mbmi->ref_frame[1],
                    scale_factor);
  store_coding_context(x, ctx, best_mode_index,
                       &best_partition,
                       &mbmi->ref_mvs[mbmi->ref_frame[0]][0],
                       &mbmi->ref_mvs[mbmi->ref_frame[1] < 0 ? 0 :
                                      mbmi->ref_frame[1]][0],
                       best_pred_diff, best_txfm_diff);

  return best_rd;
}
