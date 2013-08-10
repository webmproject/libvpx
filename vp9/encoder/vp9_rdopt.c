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
  {NEARESTMV, LAST_FRAME,   NONE},
  {DC_PRED,   INTRA_FRAME,  NONE},

  {NEARESTMV, ALTREF_FRAME, NONE},
  {NEARESTMV, GOLDEN_FRAME, NONE},
  {NEWMV,     LAST_FRAME,   NONE},
  {NEARESTMV, LAST_FRAME,   ALTREF_FRAME},
  {NEARMV,    LAST_FRAME,   NONE},
  {NEARESTMV, GOLDEN_FRAME, ALTREF_FRAME},

  {NEWMV,     GOLDEN_FRAME, NONE},
  {NEWMV,     ALTREF_FRAME, NONE},
  {NEARMV,    ALTREF_FRAME, NONE},

  {TM_PRED,   INTRA_FRAME,  NONE},

  {NEARMV,    LAST_FRAME,   ALTREF_FRAME},
  {NEWMV,     LAST_FRAME,   ALTREF_FRAME},
  {NEARMV,    GOLDEN_FRAME, NONE},
  {NEARMV,    GOLDEN_FRAME, ALTREF_FRAME},
  {NEWMV,     GOLDEN_FRAME, ALTREF_FRAME},

  {SPLITMV,   LAST_FRAME,   NONE},
  {SPLITMV,   GOLDEN_FRAME, NONE},
  {SPLITMV,   ALTREF_FRAME, NONE},
  {SPLITMV,   LAST_FRAME,   ALTREF_FRAME},
  {SPLITMV,   GOLDEN_FRAME, ALTREF_FRAME},

  {ZEROMV,    LAST_FRAME,   NONE},
  {ZEROMV,    GOLDEN_FRAME, NONE},
  {ZEROMV,    ALTREF_FRAME, NONE},
  {ZEROMV,    LAST_FRAME,   ALTREF_FRAME},
  {ZEROMV,    GOLDEN_FRAME, ALTREF_FRAME},

  {I4X4_PRED, INTRA_FRAME,  NONE},
  {H_PRED,    INTRA_FRAME,  NONE},
  {V_PRED,    INTRA_FRAME,  NONE},
  {D135_PRED, INTRA_FRAME,  NONE},
  {D27_PRED,  INTRA_FRAME,  NONE},
  {D153_PRED, INTRA_FRAME,  NONE},
  {D63_PRED,  INTRA_FRAME,  NONE},
  {D117_PRED, INTRA_FRAME,  NONE},
  {D45_PRED,  INTRA_FRAME,  NONE},
};

// The baseline rd thresholds for breaking out of the rd loop for
// certain modes are assumed to be based on 8x8 blocks.
// This table is used to correct for blocks size.
// The factors here are << 2 (2 = x0.5, 32 = x8 etc).
static int rd_thresh_block_size_factor[BLOCK_SIZES] =
  {2, 3, 3, 4, 6, 6, 8, 12, 12, 16, 24, 24, 32};

#define BASE_RD_THRESH_FREQ_FACT 16
#define MAX_RD_THRESH_FREQ_FACT 32
#define MAX_RD_THRESH_FREQ_INC 1

static void fill_token_costs(vp9_coeff_cost *c,
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
            vp9_cost_tokens((int *)c[t][i][j][k][0][l], probs,
                            vp9_coef_tree);
            vp9_cost_tokens_skip((int *)c[t][i][j][k][1][l], probs,
                                 vp9_coef_tree);
            assert(c[t][i][j][k][0][l][DCT_EOB_TOKEN] ==
                   c[t][i][j][k][1][l][DCT_EOB_TOKEN]);
          }
}

static const int rd_iifactor[32] = {
  4, 4, 3, 2, 1, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0,
};

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

    for (bsize = 0; bsize < BLOCK_SIZES; ++bsize) {
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

        if (cpi->sf.adaptive_rd_thresh)
          cpi->rd_thresh_freq_fact[bsize][i] = MAX_RD_THRESH_FREQ_FACT;
        else
          cpi->rd_thresh_freq_fact[bsize][i] = BASE_RD_THRESH_FREQ_FACT;
      }
    }
  } else {
    cpi->RDDIV = 100;

    for (bsize = 0; bsize < BLOCK_SIZES; ++bsize) {
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

        if (cpi->sf.adaptive_rd_thresh)
          cpi->rd_thresh_freq_fact[bsize][i] = MAX_RD_THRESH_FREQ_FACT;
        else
          cpi->rd_thresh_freq_fact[bsize][i] = BASE_RD_THRESH_FREQ_FACT;
      }
    }
  }

  fill_token_costs(cpi->mb.token_costs, cpi->common.fc.coef_probs);

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

    for (i = 0; i < INTER_MODE_CONTEXTS; i++) {
      MB_PREDICTION_MODE m;

      for (m = NEARESTMV; m < MB_MODE_COUNT; m++)
        cpi->mb.inter_mode_cost[i][m - NEARESTMV] =
            cost_token(vp9_inter_mode_tree,
                       cpi->common.fc.inter_mode_probs[i],
                       vp9_inter_mode_encodings - NEARESTMV + m);
    }
  }
}

static INLINE void linear_interpolate2(double x, int ntab, int inv_step,
                                       const double *tab1, const double *tab2,
                                       double *v1, double *v2) {
  double y = x * inv_step;
  int d = (int) y;
  if (d >= ntab - 1) {
    *v1 = tab1[ntab - 1];
    *v2 = tab2[ntab - 1];
  } else {
    double a = y - d;
    *v1 = tab1[d] * (1 - a) + tab1[d + 1] * a;
    *v2 = tab2[d] * (1 - a) + tab2[d + 1] * a;
  }
}

static void model_rd_norm(double x, double *R, double *D) {
  static const int inv_tab_step = 8;
  static const int tab_size = 120;
  // NOTE: The tables below must be of the same size
  //
  // Normalized rate
  // This table models the rate for a Laplacian source
  // source with given variance when quantized with a uniform quantizer
  // with given stepsize. The closed form expression is:
  // Rn(x) = H(sqrt(r)) + sqrt(r)*[1 + H(r)/(1 - r)],
  // where r = exp(-sqrt(2) * x) and x = qpstep / sqrt(variance),
  // and H(x) is the binary entropy function.
  static const double rate_tab[] = {
    64.00, 4.944, 3.949, 3.372, 2.966, 2.655, 2.403, 2.194,
    2.014, 1.858, 1.720, 1.596, 1.485, 1.384, 1.291, 1.206,
    1.127, 1.054, 0.986, 0.923, 0.863, 0.808, 0.756, 0.708,
    0.662, 0.619, 0.579, 0.541, 0.506, 0.473, 0.442, 0.412,
    0.385, 0.359, 0.335, 0.313, 0.291, 0.272, 0.253, 0.236,
    0.220, 0.204, 0.190, 0.177, 0.165, 0.153, 0.142, 0.132,
    0.123, 0.114, 0.106, 0.099, 0.091, 0.085, 0.079, 0.073,
    0.068, 0.063, 0.058, 0.054, 0.050, 0.047, 0.043, 0.040,
    0.037, 0.034, 0.032, 0.029, 0.027, 0.025, 0.023, 0.022,
    0.020, 0.019, 0.017, 0.016, 0.015, 0.014, 0.013, 0.012,
    0.011, 0.010, 0.009, 0.008, 0.008, 0.007, 0.007, 0.006,
    0.006, 0.005, 0.005, 0.005, 0.004, 0.004, 0.004, 0.003,
    0.003, 0.003, 0.003, 0.002, 0.002, 0.002, 0.002, 0.002,
    0.002, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001,
    0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.000,
  };
  // Normalized distortion
  // This table models the normalized distortion for a Laplacian source
  // source with given variance when quantized with a uniform quantizer
  // with given stepsize. The closed form expression is:
  // Dn(x) = 1 - 1/sqrt(2) * x / sinh(x/sqrt(2))
  // where x = qpstep / sqrt(variance)
  // Note the actual distortion is Dn * variance.
  static const double dist_tab[] = {
    0.000, 0.001, 0.005, 0.012, 0.021, 0.032, 0.045, 0.061,
    0.079, 0.098, 0.119, 0.142, 0.166, 0.190, 0.216, 0.242,
    0.269, 0.296, 0.324, 0.351, 0.378, 0.405, 0.432, 0.458,
    0.484, 0.509, 0.534, 0.557, 0.580, 0.603, 0.624, 0.645,
    0.664, 0.683, 0.702, 0.719, 0.735, 0.751, 0.766, 0.780,
    0.794, 0.807, 0.819, 0.830, 0.841, 0.851, 0.861, 0.870,
    0.878, 0.886, 0.894, 0.901, 0.907, 0.913, 0.919, 0.925,
    0.930, 0.935, 0.939, 0.943, 0.947, 0.951, 0.954, 0.957,
    0.960, 0.963, 0.966, 0.968, 0.971, 0.973, 0.975, 0.976,
    0.978, 0.980, 0.981, 0.982, 0.984, 0.985, 0.986, 0.987,
    0.988, 0.989, 0.990, 0.990, 0.991, 0.992, 0.992, 0.993,
    0.993, 0.994, 0.994, 0.995, 0.995, 0.996, 0.996, 0.996,
    0.996, 0.997, 0.997, 0.997, 0.997, 0.998, 0.998, 0.998,
    0.998, 0.998, 0.998, 0.999, 0.999, 0.999, 0.999, 0.999,
    0.999, 0.999, 0.999, 0.999, 0.999, 0.999, 0.999, 1.000,
  };
  /*
  assert(sizeof(rate_tab) == tab_size * sizeof(rate_tab[0]);
  assert(sizeof(dist_tab) == tab_size * sizeof(dist_tab[0]);
  assert(sizeof(rate_tab) == sizeof(dist_tab));
  */
  assert(x >= 0.0);
  linear_interpolate2(x, tab_size, inv_tab_step,
                      rate_tab, dist_tab, R, D);
}

static void model_rd_from_var_lapndz(int var, int n, int qstep,
                                     int *rate, int64_t *dist) {
  // This function models the rate and distortion for a Laplacian
  // source with given variance when quantized with a uniform quantizer
  // with given stepsize. The closed form expressions are in:
  // Hang and Chen, "Source Model for transform video coder and its
  // application - Part I: Fundamental Theory", IEEE Trans. Circ.
  // Sys. for Video Tech., April 1997.
  vp9_clear_system_state();
  if (var == 0 || n == 0) {
    *rate = 0;
    *dist = 0;
  } else {
    double D, R;
    double s2 = (double) var / n;
    double x = qstep / sqrt(s2);
    model_rd_norm(x, &R, &D);
    *rate = ((n << 8) * R + 0.5);
    *dist = (var * D + 0.5);
  }
  vp9_clear_system_state();
}

static void model_rd_for_sb(VP9_COMP *cpi, BLOCK_SIZE_TYPE bsize,
                            MACROBLOCK *x, MACROBLOCKD *xd,
                            int *out_rate_sum, int64_t *out_dist_sum) {
  // Note our transform coeffs are 8 times an orthogonal transform.
  // Hence quantizer step is also 8 times. To get effective quantizer
  // we need to divide by 8 before sending to modeling function.
  int i, rate_sum = 0, dist_sum = 0;

  for (i = 0; i < MAX_MB_PLANE; ++i) {
    struct macroblock_plane *const p = &x->plane[i];
    struct macroblockd_plane *const pd = &xd->plane[i];
    const BLOCK_SIZE_TYPE bs = get_plane_block_size(bsize, pd);
    unsigned int sse;
    int rate;
    int64_t dist;
    (void) cpi->fn_ptr[bs].vf(p->src.buf, p->src.stride,
                              pd->dst.buf, pd->dst.stride, &sse);
    // sse works better than var, since there is no dc prediction used
    model_rd_from_var_lapndz(sse, 1 << num_pels_log2_lookup[bs],
                             pd->dequant[1] >> 3, &rate, &dist);

    rate_sum += rate;
    dist_sum += dist;
  }

  *out_rate_sum = rate_sum;
  *out_dist_sum = dist_sum << 4;
}

static void model_rd_for_sb_y_tx(VP9_COMP *cpi, BLOCK_SIZE_TYPE bsize,
                                 TX_SIZE tx_size,
                                 MACROBLOCK *x, MACROBLOCKD *xd,
                                 int *out_rate_sum, int64_t *out_dist_sum,
                                 int *out_skip) {
  int t = 4, j, k;
  BLOCK_SIZE_TYPE bs = BLOCK_4X4;
  struct macroblock_plane *const p = &x->plane[0];
  struct macroblockd_plane *const pd = &xd->plane[0];
  const int width = plane_block_width(bsize, pd);
  const int height = plane_block_height(bsize, pd);
  int rate_sum = 0;
  int64_t dist_sum = 0;

  if (tx_size == TX_4X4) {
    bs = BLOCK_4X4;
    t = 4;
  } else if (tx_size == TX_8X8) {
    bs = BLOCK_8X8;
    t = 8;
  } else if (tx_size == TX_16X16) {
    bs = BLOCK_16X16;
    t = 16;
  } else if (tx_size == TX_32X32) {
    bs = BLOCK_32X32;
    t = 32;
  } else {
    assert(0);
  }
  *out_skip = 1;
  for (j = 0; j < height; j += t) {
    for (k = 0; k < width; k += t) {
      int rate;
      int64_t dist;
      unsigned int sse;
      (void) cpi->fn_ptr[bs].vf(p->src.buf + j * p->src.stride + k,
                                p->src.stride,
                                pd->dst.buf + j * pd->dst.stride + k,
                                pd->dst.stride, &sse);
      // sse works better than var, since there is no dc prediction used
      model_rd_from_var_lapndz(sse, t * t, pd->dequant[1] >> 3,
                               &rate, &dist);
      rate_sum += rate;
      dist_sum += dist;
      *out_skip &= (rate < 1024);
    }
  }
  *out_rate_sum = rate_sum;
  *out_dist_sum = (dist_sum << 4);
}

int64_t vp9_block_error_c(int16_t *coeff, int16_t *dqcoeff,
                          intptr_t block_size, int64_t *ssz) {
  int i;
  int64_t error = 0, sqcoeff = 0;

  for (i = 0; i < block_size; i++) {
    int this_diff = coeff[i] - dqcoeff[i];
    error += (unsigned)this_diff * this_diff;
    sqcoeff += (unsigned) coeff[i] * coeff[i];
  }

  *ssz = sqcoeff;
  return error;
}

/* The trailing '0' is a terminator which is used inside cost_coeffs() to
 * decide whether to include cost of a trailing EOB node or not (i.e. we
 * can skip this if the last coefficient in this transform block, e.g. the
 * 16th coefficient in a 4x4 block or the 64th coefficient in a 8x8 block,
 * were non-zero). */
static const int16_t band_counts[TX_SIZES][8] = {
  { 1, 2, 3, 4,  3,   16 - 13, 0 },
  { 1, 2, 3, 4, 11,   64 - 21, 0 },
  { 1, 2, 3, 4, 11,  256 - 21, 0 },
  { 1, 2, 3, 4, 11, 1024 - 21, 0 },
};

static INLINE int cost_coeffs(MACROBLOCK *mb,
                              int plane, int block, PLANE_TYPE type,
                              ENTROPY_CONTEXT *A, ENTROPY_CONTEXT *L,
                              TX_SIZE tx_size,
                              const int16_t *scan, const int16_t *nb) {
  MACROBLOCKD *const xd = &mb->e_mbd;
  MB_MODE_INFO *mbmi = &xd->mode_info_context->mbmi;
  int pt, c, cost;
  const int16_t *band_count = &band_counts[tx_size][1];
  const int eob = xd->plane[plane].eobs[block];
  const int16_t *qcoeff_ptr = BLOCK_OFFSET(xd->plane[plane].qcoeff, block);
  const int ref = mbmi->ref_frame[0] != INTRA_FRAME;
  unsigned int (*token_costs)[2][PREV_COEF_CONTEXTS]
                    [MAX_ENTROPY_TOKENS] = mb->token_costs[tx_size][type][ref];
  ENTROPY_CONTEXT above_ec = !!*A, left_ec = !!*L;
  uint8_t token_cache[1024];

  // Check for consistency of tx_size with mode info
  assert((!type && !plane) || (type && plane));
  if (type == PLANE_TYPE_Y_WITH_DC) {
    assert(xd->mode_info_context->mbmi.txfm_size == tx_size);
  } else {
    assert(tx_size == get_uv_tx_size(mbmi));
  }

  pt = combine_entropy_contexts(above_ec, left_ec);

  if (eob == 0) {
    // single eob token
    cost = token_costs[0][0][pt][DCT_EOB_TOKEN];
    c = 0;
  } else {
    int v, prev_t, band_left = *band_count++;

    // dc token
    v = qcoeff_ptr[0];
    prev_t = vp9_dct_value_tokens_ptr[v].token;
    cost = (*token_costs)[0][pt][prev_t] + vp9_dct_value_cost_ptr[v];
    token_cache[0] = vp9_pt_energy_class[prev_t];
    ++token_costs;

    // ac tokens
    for (c = 1; c < eob; c++) {
      const int rc = scan[c];
      int t;

      v = qcoeff_ptr[rc];
      t = vp9_dct_value_tokens_ptr[v].token;
      pt = get_coef_context(nb, token_cache, c);
      cost += (*token_costs)[!prev_t][pt][t] + vp9_dct_value_cost_ptr[v];
      token_cache[rc] = vp9_pt_energy_class[t];
      prev_t = t;
      if (!--band_left) {
        band_left = *band_count++;
        ++token_costs;
      }
    }

    // eob token
    if (band_left) {
      pt = get_coef_context(nb, token_cache, c);
      cost += (*token_costs)[0][pt][DCT_EOB_TOKEN];
    }
  }

  // is eob first coefficient;
  *A = *L = c > 0;

  return cost;
}

struct rdcost_block_args {
  VP9_COMMON *cm;
  MACROBLOCK *x;
  ENTROPY_CONTEXT t_above[16];
  ENTROPY_CONTEXT t_left[16];
  TX_SIZE tx_size;
  int bw;
  int bh;
  int rate;
  int64_t dist;
  int64_t sse;
  int64_t best_rd;
  int skip;
  const int16_t *scan, *nb;
};

static void dist_block(int plane, int block, BLOCK_SIZE_TYPE bsize,
                       int ss_txfrm_size, void *arg) {
  struct rdcost_block_args* args = arg;
  MACROBLOCK* const x = args->x;
  MACROBLOCKD* const xd = &x->e_mbd;
  struct macroblock_plane *const p = &x->plane[0];
  struct macroblockd_plane *const pd = &xd->plane[0];
  int64_t this_sse;
  int shift = args->tx_size == TX_32X32 ? 0 : 2;
  int16_t *const coeff = BLOCK_OFFSET(p->coeff, block);
  int16_t *const dqcoeff = BLOCK_OFFSET(pd->dqcoeff, block);
  args->dist += vp9_block_error(coeff, dqcoeff, 16 << ss_txfrm_size,
                                &this_sse) >> shift;
  args->sse += this_sse >> shift;

  if (x->skip_encode &&
      xd->mode_info_context->mbmi.ref_frame[0] == INTRA_FRAME) {
    // TODO(jingning): tune the model to better capture the distortion.
    int64_t p = (pd->dequant[1] * pd->dequant[1] *
                    (1 << ss_txfrm_size)) >> shift;
    args->dist += p;
    args->sse  += p;
  }
}

static void rate_block(int plane, int block, BLOCK_SIZE_TYPE bsize,
                       int ss_txfrm_size, void *arg) {
  struct rdcost_block_args* args = arg;
  int x_idx, y_idx;
  MACROBLOCKD * const xd = &args->x->e_mbd;

  txfrm_block_to_raster_xy(xd, bsize, plane, block, args->tx_size * 2, &x_idx,
                           &y_idx);

  args->rate += cost_coeffs(args->x, plane, block,
                            xd->plane[plane].plane_type, args->t_above + x_idx,
                            args->t_left + y_idx, args->tx_size,
                            args->scan, args->nb);
}

// FIXME(jingning): need to make the rd test of chroma components consistent
// with that of luma component. this function should be deprecated afterwards.
static int rdcost_plane(VP9_COMMON * const cm, MACROBLOCK *x, int plane,
                        BLOCK_SIZE_TYPE bsize, TX_SIZE tx_size) {
  MACROBLOCKD *const xd = &x->e_mbd;
  struct macroblockd_plane *pd = &xd->plane[plane];
  const BLOCK_SIZE_TYPE bs = get_plane_block_size(bsize, pd);
  const int num_4x4_blocks_wide = num_4x4_blocks_wide_lookup[bs];
  const int num_4x4_blocks_high = num_4x4_blocks_high_lookup[bs];
  int i;
  struct rdcost_block_args args = { cm, x, { 0 }, { 0 }, tx_size,
                                    num_4x4_blocks_wide, num_4x4_blocks_high,
                                    0, 0, 0, INT64_MAX, 0 };

  switch (tx_size) {
    case TX_4X4:
      vpx_memcpy(&args.t_above, pd->above_context,
                 sizeof(ENTROPY_CONTEXT) * num_4x4_blocks_wide);
      vpx_memcpy(&args.t_left, pd->left_context,
                 sizeof(ENTROPY_CONTEXT) * num_4x4_blocks_high);
      args.scan = vp9_default_scan_4x4;
      args.nb = vp9_default_scan_4x4_neighbors;
      break;
    case TX_8X8:
      for (i = 0; i < num_4x4_blocks_wide; i += 2)
        args.t_above[i] = !!*(uint16_t *)&pd->above_context[i];
      for (i = 0; i < num_4x4_blocks_high; i += 2)
        args.t_left[i] = !!*(uint16_t *)&pd->left_context[i];
      args.scan = vp9_default_scan_8x8;
      args.nb = vp9_default_scan_8x8_neighbors;
      break;
    case TX_16X16:
      for (i = 0; i < num_4x4_blocks_wide; i += 4)
        args.t_above[i] = !!*(uint32_t *)&pd->above_context[i];
      for (i = 0; i < num_4x4_blocks_high; i += 4)
        args.t_left[i] = !!*(uint32_t *)&pd->left_context[i];
      args.scan = vp9_default_scan_16x16;
      args.nb = vp9_default_scan_16x16_neighbors;
      break;
    case TX_32X32:
      for (i = 0; i < num_4x4_blocks_wide; i += 8)
        args.t_above[i] = !!*(uint64_t *)&pd->above_context[i];
      for (i = 0; i < num_4x4_blocks_high; i += 8)
        args.t_left[i] = !!*(uint64_t *)&pd->left_context[i];
      args.scan = vp9_default_scan_32x32;
      args.nb = vp9_default_scan_32x32_neighbors;
      break;
    default:
      assert(0);
  }

  foreach_transformed_block_in_plane(xd, bsize, plane, rate_block, &args);
  return args.rate;
}

static int rdcost_uv(VP9_COMMON *const cm, MACROBLOCK *x,
                     BLOCK_SIZE_TYPE bsize, TX_SIZE tx_size) {
  int cost = 0, plane;

  for (plane = 1; plane < MAX_MB_PLANE; plane++) {
    cost += rdcost_plane(cm, x, plane, bsize, tx_size);
  }
  return cost;
}

static int64_t block_error_sbuv(MACROBLOCK *x, BLOCK_SIZE_TYPE bsize,
                                int shift, int64_t *sse) {
  int64_t sum = 0, this_sse;
  int plane;

  *sse = 0;
  for (plane = 1; plane < MAX_MB_PLANE; plane++) {
    struct macroblockd_plane *pd = &x->e_mbd.plane[plane];
    const BLOCK_SIZE_TYPE bs = get_plane_block_size(bsize, pd);
    sum += vp9_block_error(x->plane[plane].coeff, pd->dqcoeff,
                           1 << num_pels_log2_lookup[bs], &this_sse);
    *sse += this_sse;
  }
  *sse >>= shift;
  return sum >> shift;
}

static void block_yrd_txfm(int plane, int block, BLOCK_SIZE_TYPE bsize,
                           int ss_txfrm_size, void *arg) {
  struct rdcost_block_args *args = arg;
  MACROBLOCK *const x = args->x;
  MACROBLOCKD *const xd = &x->e_mbd;
  struct encode_b_args encode_args = {args->cm, x, NULL};
  int64_t rd1, rd2, rd;

  if (args->skip)
    return;
  rd1 = RDCOST(x->rdmult, x->rddiv, args->rate, args->dist);
  rd2 = RDCOST(x->rdmult, x->rddiv, 0, args->sse);
  rd = MIN(rd1, rd2);
  if (rd > args->best_rd) {
    args->skip = 1;
    args->rate = INT_MAX;
    args->dist = INT64_MAX;
    args->sse  = INT64_MAX;
    return;
  }

  if (xd->mode_info_context->mbmi.ref_frame[0] == INTRA_FRAME)
    encode_block_intra(plane, block, bsize, ss_txfrm_size, &encode_args);
  else
    xform_quant(plane, block, bsize, ss_txfrm_size, &encode_args);

  dist_block(plane, block, bsize, ss_txfrm_size, args);
  rate_block(plane, block, bsize, ss_txfrm_size, args);
}

static void super_block_yrd_for_txfm(VP9_COMMON *const cm, MACROBLOCK *x,
                                     int *rate, int64_t *distortion,
                                     int *skippable, int64_t *sse,
                                     int64_t ref_best_rd,
                                     BLOCK_SIZE_TYPE bsize, TX_SIZE tx_size) {
  MACROBLOCKD *const xd = &x->e_mbd;
  struct macroblockd_plane *const pd = &xd->plane[0];
  const BLOCK_SIZE_TYPE bs = get_plane_block_size(bsize, pd);
  const int num_4x4_blocks_wide = num_4x4_blocks_wide_lookup[bs];
  const int num_4x4_blocks_high = num_4x4_blocks_high_lookup[bs];
  int i;
  struct rdcost_block_args args = { cm, x, { 0 }, { 0 }, tx_size,
                                    num_4x4_blocks_wide, num_4x4_blocks_high,
                                    0, 0, 0, ref_best_rd, 0 };
  xd->mode_info_context->mbmi.txfm_size = tx_size;
  switch (tx_size) {
    case TX_4X4:
      vpx_memcpy(&args.t_above, pd->above_context,
                 sizeof(ENTROPY_CONTEXT) * num_4x4_blocks_wide);
      vpx_memcpy(&args.t_left, pd->left_context,
                 sizeof(ENTROPY_CONTEXT) * num_4x4_blocks_high);
      get_scan_nb_4x4(get_tx_type_4x4(PLANE_TYPE_Y_WITH_DC, xd, 0),
                      &args.scan, &args.nb);
      break;
    case TX_8X8:
      for (i = 0; i < num_4x4_blocks_wide; i += 2)
        args.t_above[i] = !!*(uint16_t *)&pd->above_context[i];
      for (i = 0; i < num_4x4_blocks_high; i += 2)
        args.t_left[i] = !!*(uint16_t *)&pd->left_context[i];
      get_scan_nb_8x8(get_tx_type_8x8(PLANE_TYPE_Y_WITH_DC, xd),
                      &args.scan, &args.nb);
      break;
    case TX_16X16:
      for (i = 0; i < num_4x4_blocks_wide; i += 4)
        args.t_above[i] = !!*(uint32_t *)&pd->above_context[i];
      for (i = 0; i < num_4x4_blocks_high; i += 4)
        args.t_left[i] = !!*(uint32_t *)&pd->left_context[i];
      get_scan_nb_16x16(get_tx_type_16x16(PLANE_TYPE_Y_WITH_DC, xd),
                        &args.scan, &args.nb);
      break;
    case TX_32X32:
      for (i = 0; i < num_4x4_blocks_wide; i += 8)
        args.t_above[i] = !!*(uint64_t *)&pd->above_context[i];
      for (i = 0; i < num_4x4_blocks_high; i += 8)
        args.t_left[i] = !!*(uint64_t *)&pd->left_context[i];
      args.scan = vp9_default_scan_32x32;
      args.nb = vp9_default_scan_32x32_neighbors;
      break;
    default:
      assert(0);
  }

  foreach_transformed_block_in_plane(xd, bsize, 0, block_yrd_txfm, &args);
  *distortion = args.dist;
  *rate       = args.rate;
  *sse        = args.sse;
  *skippable  = vp9_sby_is_skippable(xd, bsize) && (!args.skip);
}

static void choose_largest_txfm_size(VP9_COMP *cpi, MACROBLOCK *x,
                                     int *rate, int64_t *distortion,
                                     int *skip, int64_t *sse,
                                     int64_t ref_best_rd,
                                     BLOCK_SIZE_TYPE bs) {
  const TX_SIZE max_txfm_size = TX_32X32
      - (bs < BLOCK_32X32) - (bs < BLOCK_16X16);
  VP9_COMMON *const cm = &cpi->common;
  MACROBLOCKD *const xd = &x->e_mbd;
  MB_MODE_INFO *const mbmi = &xd->mode_info_context->mbmi;
  if (max_txfm_size == TX_32X32 &&
      (cm->tx_mode == ALLOW_32X32 ||
       cm->tx_mode == TX_MODE_SELECT)) {
    mbmi->txfm_size = TX_32X32;
  } else if (max_txfm_size >= TX_16X16 &&
             (cm->tx_mode == ALLOW_16X16 ||
              cm->tx_mode == ALLOW_32X32 ||
              cm->tx_mode == TX_MODE_SELECT)) {
    mbmi->txfm_size = TX_16X16;
  } else if (cm->tx_mode != ONLY_4X4) {
    mbmi->txfm_size = TX_8X8;
  } else {
    mbmi->txfm_size = TX_4X4;
  }
  super_block_yrd_for_txfm(cm, x, rate, distortion, skip,
                           &sse[mbmi->txfm_size], ref_best_rd, bs,
                           mbmi->txfm_size);
  cpi->txfm_stepdown_count[0]++;
}

static void choose_txfm_size_from_rd(VP9_COMP *cpi, MACROBLOCK *x,
                                     int (*r)[2], int *rate,
                                     int64_t *d, int64_t *distortion,
                                     int *s, int *skip,
                                     int64_t tx_cache[TX_MODES],
                                     BLOCK_SIZE_TYPE bs) {
  const TX_SIZE max_tx_size = TX_32X32
      - (bs < BLOCK_32X32) - (bs < BLOCK_16X16);
  VP9_COMMON *const cm = &cpi->common;
  MACROBLOCKD *const xd = &x->e_mbd;
  MB_MODE_INFO *const mbmi = &xd->mode_info_context->mbmi;
  vp9_prob skip_prob = vp9_get_pred_prob_mbskip(cm, xd);
  int64_t rd[TX_SIZES][2];
  int n, m;
  int s0, s1;

  const vp9_prob *tx_probs = get_tx_probs2(xd, &cm->fc.tx_probs);

  for (n = TX_4X4; n <= max_tx_size; n++) {
    r[n][1] = r[n][0];
    if (r[n][0] == INT_MAX)
      continue;
    for (m = 0; m <= n - (n == max_tx_size); m++) {
      if (m == n)
        r[n][1] += vp9_cost_zero(tx_probs[m]);
      else
        r[n][1] += vp9_cost_one(tx_probs[m]);
    }
  }

  assert(skip_prob > 0);
  s0 = vp9_cost_bit(skip_prob, 0);
  s1 = vp9_cost_bit(skip_prob, 1);

  for (n = TX_4X4; n <= max_tx_size; n++) {
    if (d[n] == INT64_MAX) {
      rd[n][0] = rd[n][1] = INT64_MAX;
      continue;
    }
    if (s[n]) {
      rd[n][0] = rd[n][1] = RDCOST(x->rdmult, x->rddiv, s1, d[n]);
    } else {
      rd[n][0] = RDCOST(x->rdmult, x->rddiv, r[n][0] + s0, d[n]);
      rd[n][1] = RDCOST(x->rdmult, x->rddiv, r[n][1] + s0, d[n]);
    }
  }

  if (max_tx_size == TX_32X32 &&
      (cm->tx_mode == ALLOW_32X32 ||
       (cm->tx_mode == TX_MODE_SELECT &&
        rd[TX_32X32][1] < rd[TX_16X16][1] && rd[TX_32X32][1] < rd[TX_8X8][1] &&
        rd[TX_32X32][1] < rd[TX_4X4][1]))) {
    mbmi->txfm_size = TX_32X32;
  } else if (max_tx_size >= TX_16X16 &&
             (cm->tx_mode == ALLOW_16X16 ||
              cm->tx_mode == ALLOW_32X32 ||
              (cm->tx_mode == TX_MODE_SELECT &&
               rd[TX_16X16][1] < rd[TX_8X8][1] &&
               rd[TX_16X16][1] < rd[TX_4X4][1]))) {
    mbmi->txfm_size = TX_16X16;
  } else if (cm->tx_mode == ALLOW_8X8 ||
             cm->tx_mode == ALLOW_16X16 ||
             cm->tx_mode == ALLOW_32X32 ||
           (cm->tx_mode == TX_MODE_SELECT && rd[TX_8X8][1] < rd[TX_4X4][1])) {
    mbmi->txfm_size = TX_8X8;
  } else {
    mbmi->txfm_size = TX_4X4;
  }

  *distortion = d[mbmi->txfm_size];
  *rate       = r[mbmi->txfm_size][cm->tx_mode == TX_MODE_SELECT];
  *skip       = s[mbmi->txfm_size];

  tx_cache[ONLY_4X4] = rd[TX_4X4][0];
  tx_cache[ALLOW_8X8] = rd[TX_8X8][0];
  tx_cache[ALLOW_16X16] = rd[MIN(max_tx_size, TX_16X16)][0];
  tx_cache[ALLOW_32X32] = rd[MIN(max_tx_size, TX_32X32)][0];
  if (max_tx_size == TX_32X32 &&
      rd[TX_32X32][1] < rd[TX_16X16][1] && rd[TX_32X32][1] < rd[TX_8X8][1] &&
      rd[TX_32X32][1] < rd[TX_4X4][1])
    tx_cache[TX_MODE_SELECT] = rd[TX_32X32][1];
  else if (max_tx_size >= TX_16X16 &&
           rd[TX_16X16][1] < rd[TX_8X8][1] && rd[TX_16X16][1] < rd[TX_4X4][1])
    tx_cache[TX_MODE_SELECT] = rd[TX_16X16][1];
  else
    tx_cache[TX_MODE_SELECT] = rd[TX_4X4][1] < rd[TX_8X8][1] ?
                                 rd[TX_4X4][1] : rd[TX_8X8][1];

  if (max_tx_size == TX_32X32 &&
      rd[TX_32X32][1] < rd[TX_16X16][1] &&
      rd[TX_32X32][1] < rd[TX_8X8][1] &&
      rd[TX_32X32][1] < rd[TX_4X4][1]) {
    cpi->txfm_stepdown_count[0]++;
  } else if (max_tx_size >= TX_16X16 &&
             rd[TX_16X16][1] < rd[TX_8X8][1] &&
             rd[TX_16X16][1] < rd[TX_4X4][1]) {
    cpi->txfm_stepdown_count[max_tx_size - TX_16X16]++;
  } else if (rd[TX_8X8][1] < rd[TX_4X4][1]) {
    cpi->txfm_stepdown_count[max_tx_size - TX_8X8]++;
  } else {
    cpi->txfm_stepdown_count[max_tx_size - TX_4X4]++;
  }
}

static void choose_txfm_size_from_modelrd(VP9_COMP *cpi, MACROBLOCK *x,
                                          int (*r)[2], int *rate,
                                          int64_t *d, int64_t *distortion,
                                          int *s, int *skip, int64_t *sse,
                                          int64_t ref_best_rd,
                                          BLOCK_SIZE_TYPE bs,
                                          int *model_used) {
  const TX_SIZE max_txfm_size = TX_32X32
      - (bs < BLOCK_32X32) - (bs < BLOCK_16X16);
  VP9_COMMON *const cm = &cpi->common;
  MACROBLOCKD *const xd = &x->e_mbd;
  MB_MODE_INFO *const mbmi = &xd->mode_info_context->mbmi;
  vp9_prob skip_prob = vp9_get_pred_prob_mbskip(cm, xd);
  int64_t rd[TX_SIZES][2];
  int n, m;
  int s0, s1;
  double scale_rd[TX_SIZES] = {1.73, 1.44, 1.20, 1.00};
  // double scale_r[TX_SIZES] = {2.82, 2.00, 1.41, 1.00};

  const vp9_prob *tx_probs = get_tx_probs2(xd, &cm->fc.tx_probs);

  // for (n = TX_4X4; n <= max_txfm_size; n++)
  //   r[n][0] = (r[n][0] * scale_r[n]);

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
  for (n = TX_4X4; n <= max_txfm_size; n++) {
    rd[n][0] = (scale_rd[n] * rd[n][0]);
    rd[n][1] = (scale_rd[n] * rd[n][1]);
  }

  if (max_txfm_size == TX_32X32 &&
      (cm->tx_mode == ALLOW_32X32 ||
       (cm->tx_mode == TX_MODE_SELECT &&
        rd[TX_32X32][1] <= rd[TX_16X16][1] &&
        rd[TX_32X32][1] <= rd[TX_8X8][1] &&
        rd[TX_32X32][1] <= rd[TX_4X4][1]))) {
    mbmi->txfm_size = TX_32X32;
  } else if (max_txfm_size >= TX_16X16 &&
             (cm->tx_mode == ALLOW_16X16 ||
              cm->tx_mode == ALLOW_32X32 ||
              (cm->tx_mode == TX_MODE_SELECT &&
               rd[TX_16X16][1] <= rd[TX_8X8][1] &&
               rd[TX_16X16][1] <= rd[TX_4X4][1]))) {
    mbmi->txfm_size = TX_16X16;
  } else if (cm->tx_mode == ALLOW_8X8 ||
             cm->tx_mode == ALLOW_16X16 ||
             cm->tx_mode == ALLOW_32X32 ||
           (cm->tx_mode == TX_MODE_SELECT &&
            rd[TX_8X8][1] <= rd[TX_4X4][1])) {
    mbmi->txfm_size = TX_8X8;
  } else {
    mbmi->txfm_size = TX_4X4;
  }

  if (model_used[mbmi->txfm_size]) {
    // Actually encode using the chosen mode if a model was used, but do not
    // update the r, d costs
    super_block_yrd_for_txfm(cm, x, rate, distortion, skip,
                             &sse[mbmi->txfm_size], ref_best_rd,
                             bs, mbmi->txfm_size);
  } else {
    *distortion = d[mbmi->txfm_size];
    *rate       = r[mbmi->txfm_size][cm->tx_mode == TX_MODE_SELECT];
    *skip       = s[mbmi->txfm_size];
  }

  if (max_txfm_size == TX_32X32 &&
      rd[TX_32X32][1] <= rd[TX_16X16][1] &&
      rd[TX_32X32][1] <= rd[TX_8X8][1] &&
      rd[TX_32X32][1] <= rd[TX_4X4][1]) {
    cpi->txfm_stepdown_count[0]++;
  } else if (max_txfm_size >= TX_16X16 &&
             rd[TX_16X16][1] <= rd[TX_8X8][1] &&
             rd[TX_16X16][1] <= rd[TX_4X4][1]) {
    cpi->txfm_stepdown_count[max_txfm_size - TX_16X16]++;
  } else if (rd[TX_8X8][1] <= rd[TX_4X4][1]) {
    cpi->txfm_stepdown_count[max_txfm_size - TX_8X8]++;
  } else {
    cpi->txfm_stepdown_count[max_txfm_size - TX_4X4]++;
  }
}

static void super_block_yrd(VP9_COMP *cpi,
                            MACROBLOCK *x, int *rate, int64_t *distortion,
                            int *skip, int64_t *psse, BLOCK_SIZE_TYPE bs,
                            int64_t txfm_cache[TX_MODES],
                            int64_t ref_best_rd) {
  VP9_COMMON *const cm = &cpi->common;
  int r[TX_SIZES][2], s[TX_SIZES];
  int64_t d[TX_SIZES], sse[TX_SIZES];
  MACROBLOCKD *xd = &x->e_mbd;
  MB_MODE_INFO *const mbmi = &xd->mode_info_context->mbmi;

  assert(bs == mbmi->sb_type);
  if (mbmi->ref_frame[0] > INTRA_FRAME)
    vp9_subtract_sby(x, bs);

  if (cpi->sf.tx_size_search_method == USE_LARGESTALL ||
      (cpi->sf.tx_size_search_method != USE_FULL_RD &&
       mbmi->ref_frame[0] == INTRA_FRAME)) {
    vpx_memset(txfm_cache, 0, TX_MODES * sizeof(int64_t));
    choose_largest_txfm_size(cpi, x, rate, distortion, skip, sse,
                             ref_best_rd, bs);
    if (psse)
      *psse = sse[mbmi->txfm_size];
    return;
  }

  if (cpi->sf.tx_size_search_method == USE_LARGESTINTRA_MODELINTER &&
      mbmi->ref_frame[0] > INTRA_FRAME) {
    int model_used[TX_SIZES] = {1, 1, 1, 1};
    if (bs >= BLOCK_32X32) {
      if (model_used[TX_32X32])
        model_rd_for_sb_y_tx(cpi, bs, TX_32X32, x, xd,
                             &r[TX_32X32][0], &d[TX_32X32], &s[TX_32X32]);
      else
        super_block_yrd_for_txfm(cm, x, &r[TX_32X32][0], &d[TX_32X32],
                                 &s[TX_32X32], &sse[TX_32X32], INT64_MAX,
                                 bs, TX_32X32);
    }
    if (bs >= BLOCK_16X16) {
      if (model_used[TX_16X16])
        model_rd_for_sb_y_tx(cpi, bs, TX_16X16, x, xd,
                             &r[TX_16X16][0], &d[TX_16X16], &s[TX_16X16]);
      else
        super_block_yrd_for_txfm(cm, x, &r[TX_16X16][0], &d[TX_16X16],
                                 &s[TX_16X16], &sse[TX_16X16], INT64_MAX,
                                 bs, TX_16X16);
    }
    if (model_used[TX_8X8])
      model_rd_for_sb_y_tx(cpi, bs, TX_8X8, x, xd,
                           &r[TX_8X8][0], &d[TX_8X8], &s[TX_8X8]);
    else
      super_block_yrd_for_txfm(cm, x, &r[TX_8X8][0], &d[TX_8X8], &s[TX_8X8],
                               &sse[TX_8X8], INT64_MAX, bs, TX_8X8);

    if (model_used[TX_4X4])
      model_rd_for_sb_y_tx(cpi, bs, TX_4X4, x, xd,
                           &r[TX_4X4][0], &d[TX_4X4], &s[TX_4X4]);
    else
      super_block_yrd_for_txfm(cm, x, &r[TX_4X4][0], &d[TX_4X4], &s[TX_4X4],
                               &sse[TX_4X4], INT64_MAX, bs, TX_4X4);

    choose_txfm_size_from_modelrd(cpi, x, r, rate, d, distortion, s,
                                  skip, sse, ref_best_rd, bs, model_used);
  } else {
    if (bs >= BLOCK_32X32)
      super_block_yrd_for_txfm(cm, x, &r[TX_32X32][0], &d[TX_32X32],
                               &s[TX_32X32], &sse[TX_32X32], ref_best_rd,
                               bs, TX_32X32);
    if (bs >= BLOCK_16X16)
      super_block_yrd_for_txfm(cm, x, &r[TX_16X16][0], &d[TX_16X16],
                               &s[TX_16X16], &sse[TX_16X16], ref_best_rd,
                               bs, TX_16X16);
    super_block_yrd_for_txfm(cm, x, &r[TX_8X8][0], &d[TX_8X8], &s[TX_8X8],
                             &sse[TX_8X8], ref_best_rd, bs, TX_8X8);
    super_block_yrd_for_txfm(cm, x, &r[TX_4X4][0], &d[TX_4X4], &s[TX_4X4],
                             &sse[TX_4X4], ref_best_rd, bs, TX_4X4);
    choose_txfm_size_from_rd(cpi, x, r, rate, d, distortion, s,
                             skip, txfm_cache, bs);
  }
  if (psse)
    *psse = sse[mbmi->txfm_size];
}

static int conditional_skipintra(MB_PREDICTION_MODE mode,
                                 MB_PREDICTION_MODE best_intra_mode) {
  if (mode == D117_PRED &&
      best_intra_mode != V_PRED &&
      best_intra_mode != D135_PRED)
    return 1;
  if (mode == D63_PRED &&
      best_intra_mode != V_PRED &&
      best_intra_mode != D45_PRED)
    return 1;
  if (mode == D27_PRED &&
      best_intra_mode != H_PRED &&
      best_intra_mode != D45_PRED)
    return 1;
  if (mode == D153_PRED &&
      best_intra_mode != H_PRED &&
      best_intra_mode != D135_PRED)
    return 1;
  return 0;
}

static int64_t rd_pick_intra4x4block(VP9_COMP *cpi, MACROBLOCK *x, int ib,
                                     MB_PREDICTION_MODE *best_mode,
                                     int *bmode_costs,
                                     ENTROPY_CONTEXT *a, ENTROPY_CONTEXT *l,
                                     int *bestrate, int *bestratey,
                                     int64_t *bestdistortion,
                                     BLOCK_SIZE_TYPE bsize,
                                     int64_t rd_thresh) {
  MB_PREDICTION_MODE mode;
  MACROBLOCKD *xd = &x->e_mbd;
  int64_t best_rd = rd_thresh;
  int rate = 0;
  int64_t distortion;
  struct macroblock_plane *p = &x->plane[0];
  struct macroblockd_plane *pd = &xd->plane[0];
  const int src_stride = p->src.stride;
  const int dst_stride = pd->dst.stride;
  uint8_t *src_init = raster_block_offset_uint8(xd, BLOCK_8X8, 0, ib,
                                                p->src.buf, src_stride);
  uint8_t *dst_init = raster_block_offset_uint8(xd, BLOCK_8X8, 0, ib,
                                                pd->dst.buf, dst_stride);
  int16_t *src_diff, *coeff;

  ENTROPY_CONTEXT ta[2], tempa[2];
  ENTROPY_CONTEXT tl[2], templ[2];
  TX_TYPE tx_type = DCT_DCT;
  const int num_4x4_blocks_wide = num_4x4_blocks_wide_lookup[bsize];
  const int num_4x4_blocks_high = num_4x4_blocks_high_lookup[bsize];
  int idx, idy, block;
  uint8_t best_dst[8 * 8];

  assert(ib < 4);

  vpx_memcpy(ta, a, sizeof(ta));
  vpx_memcpy(tl, l, sizeof(tl));
  xd->mode_info_context->mbmi.txfm_size = TX_4X4;

  for (mode = DC_PRED; mode <= TM_PRED; ++mode) {
    int64_t this_rd;
    int ratey = 0;
    // Only do the oblique modes if the best so far is
    // one of the neighboring directional modes
    if (cpi->sf.mode_search_skip_flags & FLAG_SKIP_INTRA_DIRMISMATCH) {
      if (conditional_skipintra(mode, *best_mode))
          continue;
    }

    rate = bmode_costs[mode];
    distortion = 0;

    vpx_memcpy(tempa, ta, sizeof(ta));
    vpx_memcpy(templ, tl, sizeof(tl));

    for (idy = 0; idy < num_4x4_blocks_high; ++idy) {
      for (idx = 0; idx < num_4x4_blocks_wide; ++idx) {
        int64_t ssz;
        const int16_t *scan;
        uint8_t *src = src_init + idx * 4 + idy * 4 * src_stride;
        uint8_t *dst = dst_init + idx * 4 + idy * 4 * dst_stride;

        block = ib + idy * 2 + idx;
        xd->mode_info_context->bmi[block].as_mode = mode;
        src_diff = raster_block_offset_int16(xd, BLOCK_8X8, 0, block,
                                             p->src_diff);
        coeff = BLOCK_OFFSET(x->plane[0].coeff, block);
        vp9_predict_intra_block(xd, block, 1,
                                TX_4X4, mode,
                                x->skip_encode ? src : dst,
                                x->skip_encode ? src_stride : dst_stride,
                                dst, dst_stride);
        vp9_subtract_block(4, 4, src_diff, 8,
                           src, src_stride,
                           dst, dst_stride);

        tx_type = get_tx_type_4x4(PLANE_TYPE_Y_WITH_DC, xd, block);
        if (tx_type != DCT_DCT) {
          vp9_short_fht4x4(src_diff, coeff, 8, tx_type);
          x->quantize_b_4x4(x, block, tx_type, 16);
        } else {
          x->fwd_txm4x4(src_diff, coeff, 16);
          x->quantize_b_4x4(x, block, tx_type, 16);
        }

        scan = get_scan_4x4(get_tx_type_4x4(PLANE_TYPE_Y_WITH_DC, xd, block));
        ratey += cost_coeffs(x, 0, block, PLANE_TYPE_Y_WITH_DC,
                             tempa + idx, templ + idy, TX_4X4, scan,
                             vp9_get_coef_neighbors_handle(scan));
        distortion += vp9_block_error(coeff, BLOCK_OFFSET(pd->dqcoeff, block),
                                      16, &ssz) >> 2;
        if (RDCOST(x->rdmult, x->rddiv, ratey, distortion) >= best_rd)
          goto next;

        if (tx_type != DCT_DCT)
          vp9_short_iht4x4_add(BLOCK_OFFSET(pd->dqcoeff, block),
                               dst, pd->dst.stride, tx_type);
        else
          xd->inv_txm4x4_add(BLOCK_OFFSET(pd->dqcoeff, block),
                             dst, pd->dst.stride);
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
      vpx_memcpy(a, tempa, sizeof(tempa));
      vpx_memcpy(l, templ, sizeof(templ));
      for (idy = 0; idy < num_4x4_blocks_high * 4; ++idy)
        vpx_memcpy(best_dst + idy * 8, dst_init + idy * dst_stride,
                   num_4x4_blocks_wide * 4);
    }
  next:
    {}
  }

  if (best_rd >= rd_thresh || x->skip_encode)
    return best_rd;

  for (idy = 0; idy < num_4x4_blocks_high * 4; ++idy)
    vpx_memcpy(dst_init + idy * dst_stride, best_dst + idy * 8,
               num_4x4_blocks_wide * 4);

  return best_rd;
}

static int64_t rd_pick_intra_sub_8x8_y_mode(VP9_COMP * const cpi,
                                            MACROBLOCK * const mb,
                                            int * const rate,
                                            int * const rate_y,
                                            int64_t * const distortion,
                                            int64_t best_rd) {
  int i, j;
  MACROBLOCKD *const xd = &mb->e_mbd;
  BLOCK_SIZE_TYPE bsize = xd->mode_info_context->mbmi.sb_type;
  const int num_4x4_blocks_wide = num_4x4_blocks_wide_lookup[bsize];
  const int num_4x4_blocks_high = num_4x4_blocks_high_lookup[bsize];
  int idx, idy;
  int cost = 0;
  int64_t total_distortion = 0;
  int tot_rate_y = 0;
  int64_t total_rd = 0;
  ENTROPY_CONTEXT t_above[4], t_left[4];
  int *bmode_costs;
  MODE_INFO *const mic = xd->mode_info_context;

  vpx_memcpy(t_above, xd->plane[0].above_context, sizeof(t_above));
  vpx_memcpy(t_left, xd->plane[0].left_context, sizeof(t_left));

  bmode_costs = mb->mbmode_cost;

  // Pick modes for each sub-block (of size 4x4, 4x8, or 8x4) in an 8x8 block.
  for (idy = 0; idy < 2; idy += num_4x4_blocks_high) {
    for (idx = 0; idx < 2; idx += num_4x4_blocks_wide) {
      const int mis = xd->mode_info_stride;
      MB_PREDICTION_MODE UNINITIALIZED_IS_SAFE(best_mode);
      int UNINITIALIZED_IS_SAFE(r), UNINITIALIZED_IS_SAFE(ry);
      int64_t UNINITIALIZED_IS_SAFE(d), this_rd;
      i = idy * 2 + idx;

      if (cpi->common.frame_type == KEY_FRAME) {
        const MB_PREDICTION_MODE A = above_block_mode(mic, i, mis);
        const MB_PREDICTION_MODE L = (xd->left_available || idx) ?
                                     left_block_mode(mic, i) : DC_PRED;

        bmode_costs  = mb->y_mode_costs[A][L];
      }

      this_rd = rd_pick_intra4x4block(cpi, mb, i, &best_mode, bmode_costs,
                                      t_above + idx, t_left + idy,
                                      &r, &ry, &d, bsize,
                                      best_rd - total_rd);
      if (this_rd >= best_rd - total_rd)
        return INT64_MAX;

      total_rd += this_rd;
      cost += r;
      total_distortion += d;
      tot_rate_y += ry;

      mic->bmi[i].as_mode = best_mode;
      for (j = 1; j < num_4x4_blocks_high; ++j)
        mic->bmi[i + j * 2].as_mode = best_mode;
      for (j = 1; j < num_4x4_blocks_wide; ++j)
        mic->bmi[i + j].as_mode = best_mode;

      if (total_rd >= best_rd)
        return INT64_MAX;
    }
  }

  *rate = cost;
  *rate_y = tot_rate_y;
  *distortion = total_distortion;
  xd->mode_info_context->mbmi.mode = mic->bmi[3].as_mode;

  return RDCOST(mb->rdmult, mb->rddiv, cost, total_distortion);
}

static int64_t rd_pick_intra_sby_mode(VP9_COMP *cpi, MACROBLOCK *x,
                                      int *rate, int *rate_tokenonly,
                                      int64_t *distortion, int *skippable,
                                      BLOCK_SIZE_TYPE bsize,
                                      int64_t tx_cache[TX_MODES],
                                      int64_t best_rd) {
  MB_PREDICTION_MODE mode;
  MB_PREDICTION_MODE UNINITIALIZED_IS_SAFE(mode_selected);
  MACROBLOCKD *const xd = &x->e_mbd;
  int this_rate, this_rate_tokenonly, s;
  int64_t this_distortion, this_rd;
  TX_SIZE UNINITIALIZED_IS_SAFE(best_tx);
  int i;
  int *bmode_costs = x->mbmode_cost;

  if (cpi->sf.tx_size_search_method == USE_FULL_RD)
    for (i = 0; i < TX_MODES; i++)
      tx_cache[i] = INT64_MAX;

  /* Y Search for intra prediction mode */
  for (mode = DC_PRED; mode <= TM_PRED; mode++) {
    int64_t local_tx_cache[TX_MODES];
    MODE_INFO *const mic = xd->mode_info_context;
    const int mis = xd->mode_info_stride;

    if (cpi->common.frame_type == KEY_FRAME) {
      const MB_PREDICTION_MODE A = above_block_mode(mic, 0, mis);
      const MB_PREDICTION_MODE L = xd->left_available ?
                                   left_block_mode(mic, 0) : DC_PRED;

      bmode_costs = x->y_mode_costs[A][L];
    }
    x->e_mbd.mode_info_context->mbmi.mode = mode;

    super_block_yrd(cpi, x, &this_rate_tokenonly, &this_distortion, &s, NULL,
                    bsize, local_tx_cache, best_rd);

    if (this_rate_tokenonly == INT_MAX)
      continue;

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

    if (cpi->sf.tx_size_search_method == USE_FULL_RD && this_rd < INT64_MAX) {
      for (i = 0; i < TX_MODES; i++) {
        const int64_t adj_rd = this_rd + local_tx_cache[i] -
            local_tx_cache[cpi->common.tx_mode];
        if (adj_rd < tx_cache[i]) {
          tx_cache[i] = adj_rd;
        }
      }
    }
  }

  x->e_mbd.mode_info_context->mbmi.mode = mode_selected;
  x->e_mbd.mode_info_context->mbmi.txfm_size = best_tx;

  return best_rd;
}

static void super_block_uvrd_for_txfm(VP9_COMMON *const cm, MACROBLOCK *x,
                                      int *rate, int64_t *distortion,
                                      int *skippable, int64_t *sse,
                                      BLOCK_SIZE_TYPE bsize,
                                      TX_SIZE uv_tx_size) {
  MACROBLOCKD *const xd = &x->e_mbd;
  int64_t dummy;
  if (xd->mode_info_context->mbmi.ref_frame[0] == INTRA_FRAME)
    vp9_encode_intra_block_uv(cm, x, bsize);
  else
    vp9_xform_quant_sbuv(cm, x, bsize);

  *distortion = block_error_sbuv(x, bsize, uv_tx_size == TX_32X32 ? 0 : 2,
                                 sse ? sse : &dummy);
  *rate       = rdcost_uv(cm, x, bsize, uv_tx_size);
  *skippable  = vp9_sbuv_is_skippable(xd, bsize);
}

static void super_block_uvrd(VP9_COMMON *const cm, MACROBLOCK *x,
                             int *rate, int64_t *distortion, int *skippable,
                             int64_t *sse, BLOCK_SIZE_TYPE bsize) {
  MACROBLOCKD *const xd = &x->e_mbd;
  MB_MODE_INFO *const mbmi = &xd->mode_info_context->mbmi;
  TX_SIZE uv_txfm_size = get_uv_tx_size(mbmi);

  if (mbmi->ref_frame[0] > INTRA_FRAME)
    vp9_subtract_sbuv(x, bsize);

  super_block_uvrd_for_txfm(cm, x, rate, distortion, skippable, sse, bsize,
                            uv_txfm_size);
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

  MB_PREDICTION_MODE last_mode = bsize <= BLOCK_8X8 ?
              TM_PRED : cpi->sf.last_chroma_intra_mode;

  for (mode = DC_PRED; mode <= last_mode; mode++) {
    x->e_mbd.mode_info_context->mbmi.uv_mode = mode;
    super_block_uvrd(&cpi->common, x, &this_rate_tokenonly,
                     &this_distortion, &s, NULL, bsize);
    this_rate = this_rate_tokenonly +
                x->intra_uv_mode_cost[cpi->common.frame_type][mode];
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

static int64_t rd_sbuv_dcpred(VP9_COMP *cpi, MACROBLOCK *x,
                              int *rate, int *rate_tokenonly,
                              int64_t *distortion, int *skippable,
                              BLOCK_SIZE_TYPE bsize) {
  int64_t this_rd;

  x->e_mbd.mode_info_context->mbmi.uv_mode = DC_PRED;
  super_block_uvrd(&cpi->common, x, rate_tokenonly,
                   distortion, skippable, NULL, bsize);
  *rate = *rate_tokenonly +
          x->intra_uv_mode_cost[cpi->common.frame_type][DC_PRED];
  this_rd = RDCOST(x->rdmult, x->rddiv, *rate, *distortion);

  return this_rd;
}

static void choose_intra_uv_mode(VP9_COMP *cpi, BLOCK_SIZE_TYPE bsize,
                                 int *rate_uv, int *rate_uv_tokenonly,
                                 int64_t *dist_uv, int *skip_uv,
                                 MB_PREDICTION_MODE *mode_uv) {
  MACROBLOCK *const x = &cpi->mb;

  // Use an estimated rd for uv_intra based on DC_PRED if the
  // appropriate speed flag is set.
  if (cpi->sf.use_uv_intra_rd_estimate) {
    rd_sbuv_dcpred(cpi, x, rate_uv, rate_uv_tokenonly, dist_uv, skip_uv,
                   bsize < BLOCK_8X8 ? BLOCK_8X8 : bsize);
  // Else do a proper rd search for each possible transform size that may
  // be considered in the main rd loop.
  } else {
    rd_pick_intra_sbuv_mode(cpi, x,
                            rate_uv, rate_uv_tokenonly, dist_uv, skip_uv,
                            bsize < BLOCK_8X8 ? BLOCK_8X8 : bsize);
  }
  *mode_uv = x->e_mbd.mode_info_context->mbmi.uv_mode;
}

static int cost_mv_ref(VP9_COMP *cpi, MB_PREDICTION_MODE mode,
                       int mode_context) {
  MACROBLOCK *const x = &cpi->mb;
  MACROBLOCKD *const xd = &x->e_mbd;
  const int segment_id = xd->mode_info_context->mbmi.segment_id;

  // Don't account for mode here if segment skip is enabled.
  if (!vp9_segfeature_active(&xd->seg, segment_id, SEG_LVL_SKIP)) {
    assert(is_inter_mode(mode));
    return x->inter_mode_cost[mode_context][mode - NEARESTMV];
  } else {
    return 0;
  }
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
  const int num_4x4_blocks_wide = num_4x4_blocks_wide_lookup[mbmi->sb_type];
  const int num_4x4_blocks_high = num_4x4_blocks_high_lookup[mbmi->sb_type];

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
                                    102);
      if (mbmi->ref_frame[1] > 0) {
        this_second_mv->as_int = seg_mvs[mbmi->ref_frame[1]].as_int;
        thismvcost += vp9_mv_bit_cost(this_second_mv, second_best_ref_mv,
                                      mvjcost, mvcost, 102);
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

  cost = cost_mv_ref(cpi, this_mode,
                     mbmi->mb_mode_context[mbmi->ref_frame[0]]);

  mic->bmi[i].as_mv[0].as_int = this_mv->as_int;
  if (mbmi->ref_frame[1] > 0)
    mic->bmi[i].as_mv[1].as_int = this_second_mv->as_int;

  x->partition_info->bmi[i].mode = m;
  for (idy = 0; idy < num_4x4_blocks_high; ++idy)
    for (idx = 0; idx < num_4x4_blocks_wide; ++idx)
      vpx_memcpy(&mic->bmi[i + idy * 2 + idx],
                 &mic->bmi[i], sizeof(mic->bmi[i]));

  cost += thismvcost;
  return cost;
}

static int64_t encode_inter_mb_segment(VP9_COMP *cpi,
                                       MACROBLOCK *x,
                                       int64_t best_yrd,
                                       int i,
                                       int *labelyrate,
                                       int64_t *distortion, int64_t *sse,
                                       ENTROPY_CONTEXT *ta,
                                       ENTROPY_CONTEXT *tl) {
  int k;
  MACROBLOCKD *xd = &x->e_mbd;
  struct macroblockd_plane *const pd = &xd->plane[0];
  MODE_INFO *const mi = xd->mode_info_context;
  const BLOCK_SIZE_TYPE bsize = mi->mbmi.sb_type;
  const int width = plane_block_width(bsize, pd);
  const int height = plane_block_height(bsize, pd);
  int idx, idy;
  const int src_stride = x->plane[0].src.stride;
  uint8_t* const src = raster_block_offset_uint8(xd, BLOCK_8X8, 0, i,
                                                 x->plane[0].src.buf,
                                                 src_stride);
  int16_t* src_diff = raster_block_offset_int16(xd, BLOCK_8X8, 0, i,
                                                x->plane[0].src_diff);
  int16_t* coeff = BLOCK_OFFSET(x->plane[0].coeff, i);
  uint8_t* const pre = raster_block_offset_uint8(xd, BLOCK_8X8, 0, i,
                                                 pd->pre[0].buf,
                                                 pd->pre[0].stride);
  uint8_t* const dst = raster_block_offset_uint8(xd, BLOCK_8X8, 0, i,
                                                 pd->dst.buf,
                                                 pd->dst.stride);
  int64_t thisdistortion = 0, thissse = 0;
  int thisrate = 0;

  vp9_build_inter_predictor(pre, pd->pre[0].stride,
                            dst, pd->dst.stride,
                            &mi->bmi[i].as_mv[0].as_mv,
                            &xd->scale_factor[0],
                            width, height, 0, &xd->subpix, MV_PRECISION_Q3);

  if (mi->mbmi.ref_frame[1] > 0) {
    uint8_t* const second_pre =
    raster_block_offset_uint8(xd, BLOCK_8X8, 0, i,
                              pd->pre[1].buf, pd->pre[1].stride);
    vp9_build_inter_predictor(second_pre, pd->pre[1].stride,
                              dst, pd->dst.stride,
                              &mi->bmi[i].as_mv[1].as_mv,
                              &xd->scale_factor[1],
                              width, height, 1, &xd->subpix, MV_PRECISION_Q3);
  }

  vp9_subtract_block(height, width, src_diff, 8, src, src_stride,
                     dst, pd->dst.stride);

  k = i;
  for (idy = 0; idy < height / 4; ++idy) {
    for (idx = 0; idx < width / 4; ++idx) {
      int64_t ssz, rd, rd1, rd2;

      k += (idy * 2 + idx);
      src_diff = raster_block_offset_int16(xd, BLOCK_8X8, 0, k,
                                           x->plane[0].src_diff);
      coeff = BLOCK_OFFSET(x->plane[0].coeff, k);
      x->fwd_txm4x4(src_diff, coeff, 16);
      x->quantize_b_4x4(x, k, DCT_DCT, 16);
      thisdistortion += vp9_block_error(coeff, BLOCK_OFFSET(pd->dqcoeff, k),
                                        16, &ssz);
      thissse += ssz;
      thisrate += cost_coeffs(x, 0, k, PLANE_TYPE_Y_WITH_DC,
                              ta + (k & 1),
                              tl + (k >> 1), TX_4X4,
                              vp9_default_scan_4x4,
                              vp9_default_scan_4x4_neighbors);
      rd1 = RDCOST(x->rdmult, x->rddiv, thisrate, thisdistortion >> 2);
      rd2 = RDCOST(x->rdmult, x->rddiv, 0, thissse >> 2);
      rd = MIN(rd1, rd2);
      if (rd >= best_yrd)
        return INT64_MAX;
    }
  }
  *distortion = thisdistortion >> 2;
  *labelyrate = thisrate;
  *sse = thissse >> 2;

  return RDCOST(x->rdmult, x->rddiv, *labelyrate, *distortion);
}

typedef struct {
  int eobs;
  int brate;
  int byrate;
  int64_t bdist;
  int64_t bsse;
  int64_t brdcost;
  int_mv mvs[2];
  ENTROPY_CONTEXT ta[2];
  ENTROPY_CONTEXT tl[2];
} SEG_RDSTAT;

typedef struct {
  int_mv *ref_mv, *second_ref_mv;
  int_mv mvp;

  int64_t segment_rd;
  int r;
  int64_t d;
  int64_t sse;
  int segment_yrate;
  MB_PREDICTION_MODE modes[4];
  SEG_RDSTAT rdstat[4][VP9_INTER_MODES];
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

static INLINE void mi_buf_shift(MACROBLOCK *x, int i) {
  MB_MODE_INFO *mbmi = &x->e_mbd.mode_info_context->mbmi;
  x->plane[0].src.buf =
      raster_block_offset_uint8(&x->e_mbd, BLOCK_8X8, 0, i,
                                x->plane[0].src.buf,
                                x->plane[0].src.stride);
  assert(((intptr_t)x->e_mbd.plane[0].pre[0].buf & 0x7) == 0);
  x->e_mbd.plane[0].pre[0].buf =
      raster_block_offset_uint8(&x->e_mbd, BLOCK_8X8, 0, i,
                                x->e_mbd.plane[0].pre[0].buf,
                                x->e_mbd.plane[0].pre[0].stride);
  if (mbmi->ref_frame[1])
    x->e_mbd.plane[0].pre[1].buf =
        raster_block_offset_uint8(&x->e_mbd, BLOCK_8X8, 0, i,
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
                                    BEST_SEG_INFO *bsi_buf, int filter_idx,
                                    int_mv seg_mvs[4][MAX_REF_FRAMES],
                                    int mi_row, int mi_col) {
  int i, j, br = 0, idx, idy;
  int64_t bd = 0, block_sse = 0;
  MB_PREDICTION_MODE this_mode;
  MODE_INFO *mi = x->e_mbd.mode_info_context;
  MB_MODE_INFO *const mbmi = &mi->mbmi;
  const int label_count = 4;
  int64_t this_segment_rd = 0;
  int label_mv_thresh;
  int segmentyrate = 0;
  BLOCK_SIZE_TYPE bsize = mbmi->sb_type;
  const int num_4x4_blocks_wide = num_4x4_blocks_wide_lookup[bsize];
  const int num_4x4_blocks_high = num_4x4_blocks_high_lookup[bsize];
  vp9_variance_fn_ptr_t *v_fn_ptr;
  ENTROPY_CONTEXT t_above[2], t_left[2];
  BEST_SEG_INFO *bsi = bsi_buf + filter_idx;
  int mode_idx;
  int subpelmv = 1, have_ref = 0;

  vpx_memcpy(t_above, x->e_mbd.plane[0].above_context, sizeof(t_above));
  vpx_memcpy(t_left, x->e_mbd.plane[0].left_context, sizeof(t_left));

  v_fn_ptr = &cpi->fn_ptr[bsize];

  // 64 makes this threshold really big effectively
  // making it so that we very rarely check mvs on
  // segments.   setting this to 1 would make mv thresh
  // roughly equal to what it is for macroblocks
  label_mv_thresh = 1 * bsi->mvthresh / label_count;

  // Segmentation method overheads
  for (idy = 0; idy < 2; idy += num_4x4_blocks_high) {
    for (idx = 0; idx < 2; idx += num_4x4_blocks_wide) {
      // TODO(jingning,rbultje): rewrite the rate-distortion optimization
      // loop for 4x4/4x8/8x4 block coding. to be replaced with new rd loop
      int_mv mode_mv[MB_MODE_COUNT], second_mode_mv[MB_MODE_COUNT];
      int_mv frame_mv[MB_MODE_COUNT][MAX_REF_FRAMES];
      MB_PREDICTION_MODE mode_selected = ZEROMV;
      int64_t best_rd = INT64_MAX;
      i = idy * 2 + idx;

      frame_mv[ZEROMV][mbmi->ref_frame[0]].as_int = 0;
      frame_mv[ZEROMV][mbmi->ref_frame[1]].as_int = 0;
      vp9_append_sub8x8_mvs_for_idx(&cpi->common, &x->e_mbd,
                                    &frame_mv[NEARESTMV][mbmi->ref_frame[0]],
                                    &frame_mv[NEARMV][mbmi->ref_frame[0]],
                                    i, 0, mi_row, mi_col);
      if (mbmi->ref_frame[1] > 0)
        vp9_append_sub8x8_mvs_for_idx(&cpi->common, &x->e_mbd,
                                   &frame_mv[NEARESTMV][mbmi->ref_frame[1]],
                                   &frame_mv[NEARMV][mbmi->ref_frame[1]],
                                   i, 1, mi_row, mi_col);

      // search for the best motion vector on this segment
      for (this_mode = NEARESTMV; this_mode <= NEWMV; ++this_mode) {
        const struct buf_2d orig_src = x->plane[0].src;
        struct buf_2d orig_pre[2];

        mode_idx = inter_mode_offset(this_mode);
        bsi->rdstat[i][mode_idx].brdcost = INT64_MAX;

        // if we're near/nearest and mv == 0,0, compare to zeromv
        if ((this_mode == NEARMV || this_mode == NEARESTMV ||
             this_mode == ZEROMV) &&
            frame_mv[this_mode][mbmi->ref_frame[0]].as_int == 0 &&
            (mbmi->ref_frame[1] <= 0 ||
             frame_mv[this_mode][mbmi->ref_frame[1]].as_int == 0)) {
          int rfc = mbmi->mb_mode_context[mbmi->ref_frame[0]];
          int c1 = cost_mv_ref(cpi, NEARMV, rfc);
          int c2 = cost_mv_ref(cpi, NEARESTMV, rfc);
          int c3 = cost_mv_ref(cpi, ZEROMV, rfc);

          if (this_mode == NEARMV) {
            if (c1 > c3)
              continue;
          } else if (this_mode == NEARESTMV) {
            if (c2 > c3)
              continue;
          } else {
            assert(this_mode == ZEROMV);
            if (mbmi->ref_frame[1] <= 0) {
              if ((c3 >= c2 &&
                   frame_mv[NEARESTMV][mbmi->ref_frame[0]].as_int == 0) ||
                  (c3 >= c1 &&
                   frame_mv[NEARMV][mbmi->ref_frame[0]].as_int == 0))
                continue;
            } else {
              if ((c3 >= c2 &&
                   frame_mv[NEARESTMV][mbmi->ref_frame[0]].as_int == 0 &&
                   frame_mv[NEARESTMV][mbmi->ref_frame[1]].as_int == 0) ||
                  (c3 >= c1 &&
                   frame_mv[NEARMV][mbmi->ref_frame[0]].as_int == 0 &&
                   frame_mv[NEARMV][mbmi->ref_frame[1]].as_int == 0))
                continue;
            }
          }
        }

        vpx_memcpy(orig_pre, x->e_mbd.plane[0].pre, sizeof(orig_pre));
        vpx_memcpy(bsi->rdstat[i][mode_idx].ta, t_above,
                   sizeof(bsi->rdstat[i][mode_idx].ta));
        vpx_memcpy(bsi->rdstat[i][mode_idx].tl, t_left,
                   sizeof(bsi->rdstat[i][mode_idx].tl));

        // motion search for newmv (single predictor case only)
        if (mbmi->ref_frame[1] <= 0 && this_mode == NEWMV &&
            seg_mvs[i][mbmi->ref_frame[0]].as_int == INVALID_MV) {
          int step_param = 0;
          int further_steps;
          int thissme, bestsme = INT_MAX;
          int sadpb = x->sadperbit4;
          int_mv mvp_full;
          int max_mv;

          /* Is the best so far sufficiently good that we cant justify doing
           * and new motion search. */
          if (best_rd < label_mv_thresh)
            break;

          if (cpi->compressor_speed) {
            // use previous block's result as next block's MV predictor.
            if (i > 0) {
              bsi->mvp.as_int =
              x->e_mbd.mode_info_context->bmi[i - 1].as_mv[0].as_int;
              if (i == 2)
                bsi->mvp.as_int =
                x->e_mbd.mode_info_context->bmi[i - 2].as_mv[0].as_int;
            }
          }
          if (i == 0)
            max_mv = x->max_mv_context[mbmi->ref_frame[0]];
          else
            max_mv = MAX(abs(bsi->mvp.as_mv.row), abs(bsi->mvp.as_mv.col)) >> 3;
          if (cpi->sf.auto_mv_step_size && cpi->common.show_frame) {
            // Take wtd average of the step_params based on the last frame's
            // max mv magnitude and the best ref mvs of the current block for
            // the given reference.
            step_param = (vp9_init_search_range(cpi, max_mv) +
                          cpi->mv_step_param) >> 1;
          } else {
            step_param = cpi->mv_step_param;
          }

          further_steps = (MAX_MVSEARCH_STEPS - 1) - step_param;

          mvp_full.as_mv.row = bsi->mvp.as_mv.row >> 3;
          mvp_full.as_mv.col = bsi->mvp.as_mv.col >> 3;

          // adjust src pointer for this block
          mi_buf_shift(x, i);
          if (cpi->sf.search_method == HEX) {
            bestsme = vp9_hex_search(x, &mvp_full,
                                     step_param,
                                     sadpb, 1, v_fn_ptr, 1,
                                     bsi->ref_mv, &mode_mv[NEWMV]);
          } else if (cpi->sf.search_method == SQUARE) {
            bestsme = vp9_square_search(x, &mvp_full,
                                        step_param,
                                        sadpb, 1, v_fn_ptr, 1,
                                        bsi->ref_mv, &mode_mv[NEWMV]);
          } else if (cpi->sf.search_method == BIGDIA) {
            bestsme = vp9_bigdia_search(x, &mvp_full,
                                        step_param,
                                        sadpb, 1, v_fn_ptr, 1,
                                        bsi->ref_mv, &mode_mv[NEWMV]);
          } else {
            bestsme = vp9_full_pixel_diamond(cpi, x, &mvp_full, step_param,
                                             sadpb, further_steps, 0, v_fn_ptr,
                                             bsi->ref_mv, &mode_mv[NEWMV]);
          }

          // Should we do a full search (best quality only)
          if (cpi->compressor_speed == 0) {
            /* Check if mvp_full is within the range. */
            clamp_mv(&mvp_full.as_mv, x->mv_col_min, x->mv_col_max,
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
                                         0, cpi->sf.subpel_iters_per_step,
                                         x->nmvjointcost, x->mvcost,
                                         &distortion, &sse);

            // safe motion search result for use in compound prediction
            seg_mvs[i][mbmi->ref_frame[0]].as_int = mode_mv[NEWMV].as_int;
          }

          // restore src pointers
          mi_buf_restore(x, orig_src, orig_pre);
        }

        if (mbmi->ref_frame[1] > 0 && this_mode == NEWMV &&
            mbmi->interp_filter == EIGHTTAP) {
          if (seg_mvs[i][mbmi->ref_frame[1]].as_int == INVALID_MV ||
              seg_mvs[i][mbmi->ref_frame[0]].as_int == INVALID_MV)
            continue;

          // adjust src pointers
          mi_buf_shift(x, i);
          if (cpi->sf.comp_inter_joint_search_thresh <= bsize) {
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

        bsi->rdstat[i][mode_idx].brate =
            labels2mode(x, i, this_mode, &mode_mv[this_mode],
                        &second_mode_mv[this_mode], frame_mv, seg_mvs[i],
                        bsi->ref_mv, bsi->second_ref_mv, x->nmvjointcost,
                        x->mvcost, cpi);

        bsi->rdstat[i][mode_idx].mvs[0].as_int = mode_mv[this_mode].as_int;
        if (num_4x4_blocks_wide > 1)
          bsi->rdstat[i + 1][mode_idx].mvs[0].as_int =
              mode_mv[this_mode].as_int;
        if (num_4x4_blocks_high > 1)
          bsi->rdstat[i + 2][mode_idx].mvs[0].as_int =
              mode_mv[this_mode].as_int;
        if (mbmi->ref_frame[1] > 0) {
          bsi->rdstat[i][mode_idx].mvs[1].as_int =
              second_mode_mv[this_mode].as_int;
          if (num_4x4_blocks_wide > 1)
            bsi->rdstat[i + 1][mode_idx].mvs[1].as_int =
                second_mode_mv[this_mode].as_int;
          if (num_4x4_blocks_high > 1)
            bsi->rdstat[i + 2][mode_idx].mvs[1].as_int =
                second_mode_mv[this_mode].as_int;
        }

        // Trap vectors that reach beyond the UMV borders
        if (mv_check_bounds(x, &mode_mv[this_mode]))
          continue;
        if (mbmi->ref_frame[1] > 0 &&
            mv_check_bounds(x, &second_mode_mv[this_mode]))
          continue;

        if (filter_idx > 0) {
          BEST_SEG_INFO *ref_bsi = bsi_buf;
          subpelmv = (mode_mv[this_mode].as_mv.row & 0x0f) ||
                     (mode_mv[this_mode].as_mv.col & 0x0f);
          have_ref = mode_mv[this_mode].as_int ==
                     ref_bsi->rdstat[i][mode_idx].mvs[0].as_int;
          if (mbmi->ref_frame[1] > 0) {
            subpelmv |= (second_mode_mv[this_mode].as_mv.row & 0x0f) ||
                        (second_mode_mv[this_mode].as_mv.col & 0x0f);
            have_ref  &= second_mode_mv[this_mode].as_int ==
                         ref_bsi->rdstat[i][mode_idx].mvs[1].as_int;
          }

          if (filter_idx > 1 && !subpelmv && !have_ref) {
            ref_bsi = bsi_buf + 1;
            have_ref = mode_mv[this_mode].as_int ==
                       ref_bsi->rdstat[i][mode_idx].mvs[0].as_int;
            if (mbmi->ref_frame[1] > 0) {
              have_ref  &= second_mode_mv[this_mode].as_int ==
                           ref_bsi->rdstat[i][mode_idx].mvs[1].as_int;
            }
          }

          if (!subpelmv && have_ref &&
              ref_bsi->rdstat[i][mode_idx].brdcost < INT64_MAX) {
            vpx_memcpy(&bsi->rdstat[i][mode_idx], &ref_bsi->rdstat[i][mode_idx],
                       sizeof(SEG_RDSTAT));
            if (bsi->rdstat[i][mode_idx].brdcost < best_rd) {
              mode_selected = this_mode;
              best_rd = bsi->rdstat[i][mode_idx].brdcost;
            }
            continue;
          }
        }

        bsi->rdstat[i][mode_idx].brdcost =
            encode_inter_mb_segment(cpi, x,
                                    bsi->segment_rd - this_segment_rd, i,
                                    &bsi->rdstat[i][mode_idx].byrate,
                                    &bsi->rdstat[i][mode_idx].bdist,
                                    &bsi->rdstat[i][mode_idx].bsse,
                                    bsi->rdstat[i][mode_idx].ta,
                                    bsi->rdstat[i][mode_idx].tl);
        if (bsi->rdstat[i][mode_idx].brdcost < INT64_MAX) {
          bsi->rdstat[i][mode_idx].brdcost += RDCOST(x->rdmult, x->rddiv,
                                            bsi->rdstat[i][mode_idx].brate, 0);
          bsi->rdstat[i][mode_idx].brate += bsi->rdstat[i][mode_idx].byrate;
          bsi->rdstat[i][mode_idx].eobs = x->e_mbd.plane[0].eobs[i];
        }

        if (bsi->rdstat[i][mode_idx].brdcost < best_rd) {
          mode_selected = this_mode;
          best_rd = bsi->rdstat[i][mode_idx].brdcost;
        }
      } /*for each 4x4 mode*/

      if (best_rd == INT64_MAX) {
        int iy, midx;
        for (iy = i + 1; iy < 4; ++iy)
          for (midx = 0; midx < VP9_INTER_MODES; ++midx)
            bsi->rdstat[iy][midx].brdcost = INT64_MAX;
        bsi->segment_rd = INT64_MAX;
        return;
      }

      mode_idx = inter_mode_offset(mode_selected);
      vpx_memcpy(t_above, bsi->rdstat[i][mode_idx].ta, sizeof(t_above));
      vpx_memcpy(t_left, bsi->rdstat[i][mode_idx].tl, sizeof(t_left));

      labels2mode(x, i, mode_selected, &mode_mv[mode_selected],
                  &second_mode_mv[mode_selected], frame_mv, seg_mvs[i],
                  bsi->ref_mv, bsi->second_ref_mv, x->nmvjointcost,
                  x->mvcost, cpi);

      br += bsi->rdstat[i][mode_idx].brate;
      bd += bsi->rdstat[i][mode_idx].bdist;
      block_sse += bsi->rdstat[i][mode_idx].bsse;
      segmentyrate += bsi->rdstat[i][mode_idx].byrate;
      this_segment_rd += bsi->rdstat[i][mode_idx].brdcost;

      if (this_segment_rd > bsi->segment_rd) {
        int iy, midx;
        for (iy = i + 1; iy < 4; ++iy)
          for (midx = 0; midx < VP9_INTER_MODES; ++midx)
            bsi->rdstat[iy][midx].brdcost = INT64_MAX;
        bsi->segment_rd = INT64_MAX;
        return;
      }

      for (j = 1; j < num_4x4_blocks_high; ++j)
        vpx_memcpy(&x->partition_info->bmi[i + j * 2],
                   &x->partition_info->bmi[i],
                   sizeof(x->partition_info->bmi[i]));
      for (j = 1; j < num_4x4_blocks_wide; ++j)
        vpx_memcpy(&x->partition_info->bmi[i + j],
                   &x->partition_info->bmi[i],
                   sizeof(x->partition_info->bmi[i]));
    }
  } /* for each label */

  bsi->r = br;
  bsi->d = bd;
  bsi->segment_yrate = segmentyrate;
  bsi->segment_rd = this_segment_rd;
  bsi->sse = block_sse;

  // update the coding decisions
  for (i = 0; i < 4; ++i)
    bsi->modes[i] = x->partition_info->bmi[i].mode;
}

static int64_t rd_pick_best_mbsegmentation(VP9_COMP *cpi, MACROBLOCK *x,
                                           int_mv *best_ref_mv,
                                           int_mv *second_best_ref_mv,
                                           int64_t best_rd,
                                           int *returntotrate,
                                           int *returnyrate,
                                           int64_t *returndistortion,
                                           int *skippable, int64_t *psse,
                                           int mvthresh,
                                           int_mv seg_mvs[4][MAX_REF_FRAMES],
                                           BEST_SEG_INFO *bsi_buf,
                                           int filter_idx,
                                           int mi_row, int mi_col) {
  int i;
  BEST_SEG_INFO *bsi = bsi_buf + filter_idx;
  MACROBLOCKD *xd = &x->e_mbd;
  MODE_INFO *mi = xd->mode_info_context;
  MB_MODE_INFO *mbmi = &mi->mbmi;
  int mode_idx;

  vp9_zero(*bsi);

  bsi->segment_rd = best_rd;
  bsi->ref_mv = best_ref_mv;
  bsi->second_ref_mv = second_best_ref_mv;
  bsi->mvp.as_int = best_ref_mv->as_int;
  bsi->mvthresh = mvthresh;

  for (i = 0; i < 4; i++)
    bsi->modes[i] = ZEROMV;

  rd_check_segment_txsize(cpi, x, bsi_buf, filter_idx, seg_mvs, mi_row, mi_col);

  if (bsi->segment_rd > best_rd)
    return INT64_MAX;
  /* set it to the best */
  for (i = 0; i < 4; i++) {
    mode_idx = inter_mode_offset(bsi->modes[i]);
    mi->bmi[i].as_mv[0].as_int = bsi->rdstat[i][mode_idx].mvs[0].as_int;
    if (mbmi->ref_frame[1] > 0)
      mi->bmi[i].as_mv[1].as_int = bsi->rdstat[i][mode_idx].mvs[1].as_int;
    xd->plane[0].eobs[i] = bsi->rdstat[i][mode_idx].eobs;
    x->partition_info->bmi[i].mode = bsi->modes[i];
  }

  /*
   * used to set mbmi->mv.as_int
   */
  *returntotrate = bsi->r;
  *returndistortion = bsi->d;
  *returnyrate = bsi->segment_yrate;
  *skippable = vp9_sby_is_skippable(&x->e_mbd, BLOCK_8X8);
  *psse = bsi->sse;
  mbmi->mode = bsi->modes[3];

  return bsi->segment_rd;
}

static void mv_pred(VP9_COMP *cpi, MACROBLOCK *x,
                    uint8_t *ref_y_buffer, int ref_y_stride,
                    int ref_frame, BLOCK_SIZE_TYPE block_size ) {
  MACROBLOCKD *xd = &x->e_mbd;
  MB_MODE_INFO *mbmi = &xd->mode_info_context->mbmi;
  int_mv this_mv;
  int i;
  int zero_seen = 0;
  int best_index = 0;
  int best_sad = INT_MAX;
  int this_sad = INT_MAX;
  unsigned int max_mv = 0;

  uint8_t *src_y_ptr = x->plane[0].src.buf;
  uint8_t *ref_y_ptr;
  int row_offset, col_offset;

  // Get the sad for each candidate reference mv
  for (i = 0; i < MAX_MV_REF_CANDIDATES; i++) {
    this_mv.as_int = mbmi->ref_mvs[ref_frame][i].as_int;

    max_mv = MAX(max_mv,
                 MAX(abs(this_mv.as_mv.row), abs(this_mv.as_mv.col)) >> 3);
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
  x->max_mv_context[ref_frame] = max_mv;
}

static void estimate_ref_frame_costs(VP9_COMP *cpi, int segment_id,
                                     unsigned int *ref_costs_single,
                                     unsigned int *ref_costs_comp,
                                     vp9_prob *comp_mode_p) {
  VP9_COMMON *const cm = &cpi->common;
  MACROBLOCKD *const xd = &cpi->mb.e_mbd;
  int seg_ref_active = vp9_segfeature_active(&xd->seg, segment_id,
                                             SEG_LVL_REF_FRAME);
  if (seg_ref_active) {
    vpx_memset(ref_costs_single, 0, MAX_REF_FRAMES * sizeof(*ref_costs_single));
    vpx_memset(ref_costs_comp,   0, MAX_REF_FRAMES * sizeof(*ref_costs_comp));
    *comp_mode_p = 128;
  } else {
    vp9_prob intra_inter_p = vp9_get_pred_prob_intra_inter(cm, xd);
    vp9_prob comp_inter_p = 128;

    if (cm->comp_pred_mode == HYBRID_PREDICTION) {
      comp_inter_p = vp9_get_pred_prob_comp_inter_inter(cm, xd);
      *comp_mode_p = comp_inter_p;
    } else {
      *comp_mode_p = 128;
    }

    ref_costs_single[INTRA_FRAME] = vp9_cost_bit(intra_inter_p, 0);

    if (cm->comp_pred_mode != COMP_PREDICTION_ONLY) {
      vp9_prob ref_single_p1 = vp9_get_pred_prob_single_ref_p1(cm, xd);
      vp9_prob ref_single_p2 = vp9_get_pred_prob_single_ref_p2(cm, xd);
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
      vp9_prob ref_comp_p = vp9_get_pred_prob_comp_ref_p(cm, xd);
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
                         int64_t tx_size_diff[TX_MODES],
                         int64_t best_filter_diff[VP9_SWITCHABLE_FILTERS + 1]) {
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

  // FIXME(rbultje) does this memcpy the whole array? I believe sizeof()
  // doesn't actually work this way
  memcpy(ctx->tx_rd_diff, tx_size_diff, sizeof(ctx->tx_rd_diff));
  memcpy(ctx->best_filter_diff, best_filter_diff,
         sizeof(*best_filter_diff) * (VP9_SWITCHABLE_FILTERS + 1));
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
                               BLOCK_SIZE_TYPE block_size,
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
                   cpi->common.ref_frame_sign_bias, mi_row, mi_col);

  // Candidate refinement carried out at encoder and decoder
  vp9_find_best_ref_mvs(xd,
                        mbmi->ref_mvs[frame_type],
                        &frame_nearest_mv[frame_type],
                        &frame_near_mv[frame_type]);

  // Further refinement that is encode side only to test the top few candidates
  // in full and choose the best as the centre point for subsequent searches.
  // The current implementation doesn't support scaling.
  if (scale[frame_type].x_scale_fp == VP9_REF_NO_SCALE &&
      scale[frame_type].y_scale_fp == VP9_REF_NO_SCALE)
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

static INLINE int get_switchable_rate(const MACROBLOCK *x) {
  const MACROBLOCKD *const xd = &x->e_mbd;
  const MB_MODE_INFO *const mbmi = &xd->mode_info_context->mbmi;
  const int ctx = vp9_get_pred_context_switchable_interp(xd);
  return SWITCHABLE_INTERP_RATE_FACTOR *
             x->switchable_interp_costs[ctx][mbmi->interp_filter];
}

static void single_motion_search(VP9_COMP *cpi, MACROBLOCK *x,
                                 BLOCK_SIZE_TYPE bsize,
                                 int mi_row, int mi_col,
                                 int_mv *tmp_mv, int *rate_mv) {
  MACROBLOCKD *xd = &x->e_mbd;
  VP9_COMMON *cm = &cpi->common;
  MB_MODE_INFO *mbmi = &xd->mode_info_context->mbmi;
  struct buf_2d backup_yv12[MAX_MB_PLANE] = {{0}};
  int bestsme = INT_MAX;
  int further_steps, step_param;
  int sadpb = x->sadperbit16;
  int_mv mvp_full;
  int ref = mbmi->ref_frame[0];
  int_mv ref_mv = mbmi->ref_mvs[ref][0];
  const BLOCK_SIZE_TYPE block_size = get_plane_block_size(bsize, &xd->plane[0]);

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

    setup_pre_planes(xd, 0, scaled_ref_frame, mi_row, mi_col, NULL);
  }

  vp9_clamp_mv_min_max(x, &ref_mv);

  // Adjust search parameters based on small partitions' result.
  if (x->fast_ms) {
    // && abs(mvp_full.as_mv.row - x->pred_mv.as_mv.row) < 24 &&
    // abs(mvp_full.as_mv.col - x->pred_mv.as_mv.col) < 24) {
    // adjust search range
    step_param = 6;
    if (x->fast_ms > 1)
      step_param = 8;

    // Get prediction MV.
    mvp_full.as_int = x->pred_mv.as_int;

    // Adjust MV sign if needed.
    if (cm->ref_frame_sign_bias[ref]) {
      mvp_full.as_mv.col *= -1;
      mvp_full.as_mv.row *= -1;
    }
  } else {
    // Work out the size of the first step in the mv step search.
    // 0 here is maximum length first step. 1 is MAX >> 1 etc.
    if (cpi->sf.auto_mv_step_size && cpi->common.show_frame) {
      // Take wtd average of the step_params based on the last frame's
      // max mv magnitude and that based on the best ref mvs of the current
      // block for the given reference.
      step_param = (vp9_init_search_range(cpi, x->max_mv_context[ref]) +
                    cpi->mv_step_param) >> 1;
    } else {
      step_param = cpi->mv_step_param;
    }
    // mvp_full.as_int = ref_mv[0].as_int;
    mvp_full.as_int =
        mbmi->ref_mvs[ref][x->mv_best_ref_index[ref]].as_int;
  }

  mvp_full.as_mv.col >>= 3;
  mvp_full.as_mv.row >>= 3;

  // Further step/diamond searches as necessary
  further_steps = (cpi->sf.max_step_search_steps - 1) - step_param;

  if (cpi->sf.search_method == HEX) {
    bestsme = vp9_hex_search(x, &mvp_full,
                             step_param,
                             sadpb, 1,
                             &cpi->fn_ptr[block_size], 1,
                             &ref_mv, tmp_mv);
  } else if (cpi->sf.search_method == SQUARE) {
    bestsme = vp9_square_search(x, &mvp_full,
                                step_param,
                                sadpb, 1,
                                &cpi->fn_ptr[block_size], 1,
                                &ref_mv, tmp_mv);
  } else if (cpi->sf.search_method == BIGDIA) {
    bestsme = vp9_bigdia_search(x, &mvp_full,
                                step_param,
                                sadpb, 1,
                                &cpi->fn_ptr[block_size], 1,
                                &ref_mv, tmp_mv);
  } else {
    bestsme = vp9_full_pixel_diamond(cpi, x, &mvp_full, step_param,
                                     sadpb, further_steps, 1,
                                     &cpi->fn_ptr[block_size],
                                     &ref_mv, tmp_mv);
  }

  x->mv_col_min = tmp_col_min;
  x->mv_col_max = tmp_col_max;
  x->mv_row_min = tmp_row_min;
  x->mv_row_max = tmp_row_max;

  if (bestsme < INT_MAX) {
    int dis;  /* TODO: use dis in distortion calculation later. */
    unsigned int sse;
    cpi->find_fractional_mv_step(x, tmp_mv, &ref_mv,
                                 x->errorperbit,
                                 &cpi->fn_ptr[block_size],
                                 0, cpi->sf.subpel_iters_per_step,
                                 x->nmvjointcost, x->mvcost,
                                 &dis, &sse);
  }
  *rate_mv = vp9_mv_bit_cost(tmp_mv, &ref_mv,
                             x->nmvjointcost, x->mvcost,
                             96);
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
  const BLOCK_SIZE_TYPE block_size = get_plane_block_size(bsize, &xd->plane[0]);
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
    setup_pre_planes(xd, 0, scaled_ref_frame[0], mi_row, mi_col, NULL);
  }

  if (scaled_ref_frame[1]) {
    int i;
    for (i = 0; i < MAX_MB_PLANE; i++)
      backup_second_yv12[i] = xd->plane[i].pre[1];

    setup_pre_planes(xd, 0, scaled_ref_frame[1], mi_row, mi_col, NULL);
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
                              &frame_mv[refs[!id]].as_mv,
                              &xd->scale_factor[!id],
                              pw, ph, 0,
                              &xd->subpix, MV_PRECISION_Q3);

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

      bestsme = cpi->find_fractional_mv_step_comp(
          x, &tmp_mv,
          &ref_mv[id],
          x->errorperbit,
          &cpi->fn_ptr[block_size],
          0, cpi->sf.subpel_iters_per_step,
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
                              x->nmvjointcost, x->mvcost, 96);
  *rate_mv += vp9_mv_bit_cost(&frame_mv[refs[1]],
                              &mbmi->ref_mvs[refs[1]][0],
                              x->nmvjointcost, x->mvcost, 96);

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
                                 int_mv (*mode_mv)[MAX_REF_FRAMES],
                                 int mi_row, int mi_col,
                                 int_mv single_newmv[MAX_REF_FRAMES],
                                 int64_t *psse, int64_t ref_best_rd) {
  VP9_COMMON *cm = &cpi->common;
  MACROBLOCKD *xd = &x->e_mbd;
  MB_MODE_INFO *mbmi = &xd->mode_info_context->mbmi;
  const int is_comp_pred = (mbmi->ref_frame[1] > 0);
  const int num_refs = is_comp_pred ? 2 : 1;
  const int this_mode = mbmi->mode;
  int_mv *frame_mv = mode_mv[this_mode];
  int i;
  int refs[2] = { mbmi->ref_frame[0],
    (mbmi->ref_frame[1] < 0 ? 0 : mbmi->ref_frame[1]) };
  int_mv cur_mv[2];
  int64_t this_rd = 0;
  DECLARE_ALIGNED_ARRAY(16, uint8_t, tmp_buf, MAX_MB_PLANE * 64 * 64);
  int pred_exists = 0;
  int interpolating_intpel_seen = 0;
  int intpel_mv;
  int64_t rd, best_rd = INT64_MAX;
  int best_needs_copy = 0;
  uint8_t *orig_dst[MAX_MB_PLANE];
  int orig_dst_stride[MAX_MB_PLANE];
  int rs = 0;

  if (this_mode == NEWMV) {
    int rate_mv;
    if (is_comp_pred) {
      // Initialize mv using single prediction mode result.
      frame_mv[refs[0]].as_int = single_newmv[refs[0]].as_int;
      frame_mv[refs[1]].as_int = single_newmv[refs[1]].as_int;

      if (cpi->sf.comp_inter_joint_search_thresh <= bsize) {
        joint_motion_search(cpi, x, bsize, frame_mv,
                            mi_row, mi_col, single_newmv, &rate_mv);
      } else {
        rate_mv  = vp9_mv_bit_cost(&frame_mv[refs[0]],
                                   &mbmi->ref_mvs[refs[0]][0],
                                   x->nmvjointcost, x->mvcost, 96);
        rate_mv += vp9_mv_bit_cost(&frame_mv[refs[1]],
                                   &mbmi->ref_mvs[refs[1]][0],
                                   x->nmvjointcost, x->mvcost, 96);
      }
      if (frame_mv[refs[0]].as_int == INVALID_MV ||
          frame_mv[refs[1]].as_int == INVALID_MV)
        return INT64_MAX;
      *rate2 += rate_mv;
    } else {
      int_mv tmp_mv;
      single_motion_search(cpi, x, bsize, mi_row, mi_col, &tmp_mv, &rate_mv);
      *rate2 += rate_mv;
      frame_mv[refs[0]].as_int =
          xd->mode_info_context->bmi[0].as_mv[0].as_int = tmp_mv.as_int;
      single_newmv[refs[0]].as_int = tmp_mv.as_int;
    }
  }

  // if we're near/nearest and mv == 0,0, compare to zeromv
  if ((this_mode == NEARMV || this_mode == NEARESTMV || this_mode == ZEROMV) &&
      frame_mv[refs[0]].as_int == 0 &&
      !vp9_segfeature_active(&xd->seg, mbmi->segment_id, SEG_LVL_SKIP) &&
      (num_refs == 1 || frame_mv[refs[1]].as_int == 0)) {
    int rfc = mbmi->mb_mode_context[mbmi->ref_frame[0]];
    int c1 = cost_mv_ref(cpi, NEARMV, rfc);
    int c2 = cost_mv_ref(cpi, NEARESTMV, rfc);
    int c3 = cost_mv_ref(cpi, ZEROMV, rfc);

    if (this_mode == NEARMV) {
      if (c1 > c3)
        return INT64_MAX;
    } else if (this_mode == NEARESTMV) {
      if (c2 > c3)
        return INT64_MAX;
    } else {
      assert(this_mode == ZEROMV);
      if (num_refs == 1) {
        if ((c3 >= c2 &&
             mode_mv[NEARESTMV][mbmi->ref_frame[0]].as_int == 0) ||
            (c3 >= c1 &&
             mode_mv[NEARMV][mbmi->ref_frame[0]].as_int == 0))
          return INT64_MAX;
      } else {
        if ((c3 >= c2 &&
             mode_mv[NEARESTMV][mbmi->ref_frame[0]].as_int == 0 &&
             mode_mv[NEARESTMV][mbmi->ref_frame[1]].as_int == 0) ||
            (c3 >= c1 &&
             mode_mv[NEARMV][mbmi->ref_frame[0]].as_int == 0 &&
             mode_mv[NEARMV][mbmi->ref_frame[1]].as_int == 0))
          return INT64_MAX;
      }
    }
  }

  for (i = 0; i < num_refs; ++i) {
    cur_mv[i] = frame_mv[refs[i]];
    // Clip "next_nearest" so that it does not extend to far out of image
    if (this_mode != NEWMV)
      clamp_mv2(&cur_mv[i].as_mv, xd);

    if (mv_check_bounds(x, &cur_mv[i]))
      return INT64_MAX;
    mbmi->mv[i].as_int = cur_mv[i].as_int;
  }

  // do first prediction into the destination buffer. Do the next
  // prediction into a temporary buffer. Then keep track of which one
  // of these currently holds the best predictor, and use the other
  // one for future predictions. In the end, copy from tmp_buf to
  // dst if necessary.
  for (i = 0; i < MAX_MB_PLANE; i++) {
    orig_dst[i] = xd->plane[i].dst.buf;
    orig_dst_stride[i] = xd->plane[i].dst.stride;
  }

  /* We don't include the cost of the second reference here, because there
   * are only three options: Last/Golden, ARF/Last or Golden/ARF, or in other
   * words if you present them in that order, the second one is always known
   * if the first is known */
  *rate2 += cost_mv_ref(cpi, this_mode,
                        mbmi->mb_mode_context[mbmi->ref_frame[0]]);

  if (!(*mode_excluded)) {
    if (is_comp_pred) {
      *mode_excluded = (cpi->common.comp_pred_mode == SINGLE_PREDICTION_ONLY);
    } else {
      *mode_excluded = (cpi->common.comp_pred_mode == COMP_PREDICTION_ONLY);
    }
  }

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
  *best_filter = EIGHTTAP;
  if (cpi->sf.use_8tap_always) {
    *best_filter = EIGHTTAP;
    vp9_zero(cpi->rd_filter_cache);
  } else {
    int i, newbest;
    int tmp_rate_sum = 0;
    int64_t tmp_dist_sum = 0;

    cpi->rd_filter_cache[VP9_SWITCHABLE_FILTERS] = INT64_MAX;
    for (i = 0; i < VP9_SWITCHABLE_FILTERS; ++i) {
      int j;
      int64_t rs_rd;
      const int is_intpel_interp = intpel_mv;
      mbmi->interp_filter = i;
      vp9_setup_interp_filters(xd, mbmi->interp_filter, cm);
      rs = get_switchable_rate(x);
      rs_rd = RDCOST(x->rdmult, x->rddiv, rs, 0);

      if (interpolating_intpel_seen && is_intpel_interp) {
        cpi->rd_filter_cache[i] = RDCOST(x->rdmult, x->rddiv,
                                         tmp_rate_sum, tmp_dist_sum);
        cpi->rd_filter_cache[VP9_SWITCHABLE_FILTERS] =
            MIN(cpi->rd_filter_cache[VP9_SWITCHABLE_FILTERS],
                cpi->rd_filter_cache[i] + rs_rd);
        rd = cpi->rd_filter_cache[i];
        if (cm->mcomp_filter_type == SWITCHABLE)
          rd += rs_rd;
      } else {
        int rate_sum = 0;
        int64_t dist_sum = 0;
        if ((cm->mcomp_filter_type == SWITCHABLE &&
             (!i || best_needs_copy)) ||
            (cm->mcomp_filter_type != SWITCHABLE &&
             (cm->mcomp_filter_type == mbmi->interp_filter ||
              (!interpolating_intpel_seen && is_intpel_interp)))) {
          for (j = 0; j < MAX_MB_PLANE; j++) {
            xd->plane[j].dst.buf = orig_dst[j];
            xd->plane[j].dst.stride = orig_dst_stride[j];
          }
        } else {
          for (j = 0; j < MAX_MB_PLANE; j++) {
            xd->plane[j].dst.buf = tmp_buf + j * 64 * 64;
            xd->plane[j].dst.stride = 64;
          }
        }
        vp9_build_inter_predictors_sb(xd, mi_row, mi_col, bsize);
        model_rd_for_sb(cpi, bsize, x, xd, &rate_sum, &dist_sum);
        cpi->rd_filter_cache[i] = RDCOST(x->rdmult, x->rddiv,
                                         rate_sum, dist_sum);
        cpi->rd_filter_cache[VP9_SWITCHABLE_FILTERS] =
            MIN(cpi->rd_filter_cache[VP9_SWITCHABLE_FILTERS],
                cpi->rd_filter_cache[i] + rs_rd);
        rd = cpi->rd_filter_cache[i];
        if (cm->mcomp_filter_type == SWITCHABLE)
          rd += rs_rd;
        if (!interpolating_intpel_seen && is_intpel_interp) {
          tmp_rate_sum = rate_sum;
          tmp_dist_sum = dist_sum;
        }
      }
      if (i == 0 && cpi->sf.use_rd_breakout && ref_best_rd < INT64_MAX) {
        if (rd / 2 > ref_best_rd) {
          for (i = 0; i < MAX_MB_PLANE; i++) {
            xd->plane[i].dst.buf = orig_dst[i];
            xd->plane[i].dst.stride = orig_dst_stride[i];
          }
          return INT64_MAX;
        }
      }
      newbest = i == 0 || rd < best_rd;

      if (newbest) {
        best_rd = rd;
        *best_filter = mbmi->interp_filter;
        if (cm->mcomp_filter_type == SWITCHABLE && i &&
            !(interpolating_intpel_seen && is_intpel_interp))
          best_needs_copy = !best_needs_copy;
      }

      if ((cm->mcomp_filter_type == SWITCHABLE && newbest) ||
          (cm->mcomp_filter_type != SWITCHABLE &&
           cm->mcomp_filter_type == mbmi->interp_filter)) {
        pred_exists = 1;
      }
      interpolating_intpel_seen |= is_intpel_interp;
    }

    for (i = 0; i < MAX_MB_PLANE; i++) {
      xd->plane[i].dst.buf = orig_dst[i];
      xd->plane[i].dst.stride = orig_dst_stride[i];
    }
  }
  // Set the appropriate filter
  mbmi->interp_filter = cm->mcomp_filter_type != SWITCHABLE ?
      cm->mcomp_filter_type : *best_filter;
  vp9_setup_interp_filters(xd, mbmi->interp_filter, cm);
  rs = cm->mcomp_filter_type == SWITCHABLE ? get_switchable_rate(x) : 0;

  if (pred_exists) {
    if (best_needs_copy) {
      // again temporarily set the buffers to local memory to prevent a memcpy
      for (i = 0; i < MAX_MB_PLANE; i++) {
        xd->plane[i].dst.buf = tmp_buf + i * 64 * 64;
        xd->plane[i].dst.stride = 64;
      }
    }
  } else {
    // Handles the special case when a filter that is not in the
    // switchable list (ex. bilinear, 6-tap) is indicated at the frame level
    vp9_build_inter_predictors_sb(xd, mi_row, mi_col, bsize);
  }


  if (cpi->sf.use_rd_breakout && ref_best_rd < INT64_MAX) {
    int tmp_rate;
    int64_t tmp_dist;
    model_rd_for_sb(cpi, bsize, x, xd, &tmp_rate, &tmp_dist);
    rd = RDCOST(x->rdmult, x->rddiv, rs + tmp_rate, tmp_dist);
    // if current pred_error modeled rd is substantially more than the best
    // so far, do not bother doing full rd
    if (rd / 2 > ref_best_rd) {
      for (i = 0; i < MAX_MB_PLANE; i++) {
        xd->plane[i].dst.buf = orig_dst[i];
        xd->plane[i].dst.stride = orig_dst_stride[i];
      }
      return INT64_MAX;
    }
  }

  if (cpi->common.mcomp_filter_type == SWITCHABLE)
    *rate2 += get_switchable_rate(x);

  if (!is_comp_pred) {
    if (cpi->active_map_enabled && x->active_ptr[0] == 0)
      x->skip = 1;
    else if (x->encode_breakout) {
      const BLOCK_SIZE_TYPE y_size = get_plane_block_size(bsize, &xd->plane[0]);
      const BLOCK_SIZE_TYPE uv_size = get_plane_block_size(bsize,
                                                           &xd->plane[1]);
      unsigned int var, sse;
      // Skipping threshold for ac.
      unsigned int thresh_ac;
      // The encode_breakout input
      unsigned int encode_breakout = x->encode_breakout << 4;

      // Calculate threshold according to dequant value.
      thresh_ac = (xd->plane[0].dequant[1] * xd->plane[0].dequant[1]) / 9;

      // Set a maximum for threshold to avoid big PSNR loss in low bitrate case.
      if (thresh_ac > 36000)
        thresh_ac = 36000;

      // Use encode_breakout input if it is bigger than internal threshold.
      if (thresh_ac < encode_breakout)
        thresh_ac = encode_breakout;

      var = cpi->fn_ptr[y_size].vf(x->plane[0].src.buf, x->plane[0].src.stride,
                                   xd->plane[0].dst.buf,
                                   xd->plane[0].dst.stride, &sse);

      // Adjust threshold according to partition size.
      thresh_ac >>= 8 - (b_width_log2_lookup[bsize] +
          b_height_log2_lookup[bsize]);

      // Y skipping condition checking
      if (sse < thresh_ac || sse == 0) {
        // Skipping threshold for dc
        unsigned int thresh_dc;

        thresh_dc = (xd->plane[0].dequant[0] * xd->plane[0].dequant[0] >> 6);

        // dc skipping checking
        if ((sse - var) < thresh_dc || sse == var) {
          unsigned int sse_u, sse_v;
          unsigned int var_u, var_v;

          var_u = cpi->fn_ptr[uv_size].vf(x->plane[1].src.buf,
                                          x->plane[1].src.stride,
                                          xd->plane[1].dst.buf,
                                          xd->plane[1].dst.stride, &sse_u);

          // U skipping condition checking
          if ((sse_u * 4 < thresh_ac || sse_u == 0) &&
              (sse_u - var_u < thresh_dc || sse_u == var_u)) {
            var_v = cpi->fn_ptr[uv_size].vf(x->plane[2].src.buf,
                                            x->plane[2].src.stride,
                                            xd->plane[2].dst.buf,
                                            xd->plane[2].dst.stride, &sse_v);

            // V skipping condition checking
            if ((sse_v * 4 < thresh_ac || sse_v == 0) &&
                (sse_v - var_v < thresh_dc || sse_v == var_v)) {
              x->skip = 1;

              *rate2 = 500;
              *rate_uv = 0;

              // Scaling factor for SSE from spatial domain to frequency domain
              // is 16. Adjust distortion accordingly.
              *distortion_uv = (sse_u + sse_v) << 4;
              *distortion = (sse << 4) + *distortion_uv;

              *disable_skip = 1;
              this_rd = RDCOST(x->rdmult, x->rddiv, *rate2, *distortion);
            }
          }
        }
      }
    }
  }

  if (!x->skip) {
    int skippable_y, skippable_uv;
    int64_t sseuv = INT_MAX;

    // Y cost and distortion
    super_block_yrd(cpi, x, rate_y, distortion_y, &skippable_y, psse,
                    bsize, txfm_cache, ref_best_rd);

    if (*rate_y == INT_MAX) {
      *rate2 = INT_MAX;
      *distortion = INT64_MAX;
      for (i = 0; i < MAX_MB_PLANE; i++) {
        xd->plane[i].dst.buf = orig_dst[i];
        xd->plane[i].dst.stride = orig_dst_stride[i];
      }
      return INT64_MAX;
    }

    *rate2 += *rate_y;
    *distortion += *distortion_y;

    super_block_uvrd(cm, x, rate_uv, distortion_uv,
                     &skippable_uv, &sseuv, bsize);

    *psse += sseuv;
    *rate2 += *rate_uv;
    *distortion += *distortion_uv;
    *skippable = skippable_y && skippable_uv;
  }

  for (i = 0; i < MAX_MB_PLANE; i++) {
    xd->plane[i].dst.buf = orig_dst[i];
    xd->plane[i].dst.stride = orig_dst_stride[i];
  }

  return this_rd;  // if 0, this will be re-calculated by caller
}

void vp9_rd_pick_intra_mode_sb(VP9_COMP *cpi, MACROBLOCK *x,
                               int *returnrate, int64_t *returndist,
                               BLOCK_SIZE_TYPE bsize,
                               PICK_MODE_CONTEXT *ctx, int64_t best_rd) {
  VP9_COMMON *const cm = &cpi->common;
  MACROBLOCKD *const xd = &x->e_mbd;
  int rate_y = 0, rate_uv = 0, rate_y_tokenonly = 0, rate_uv_tokenonly = 0;
  int y_skip = 0, uv_skip;
  int64_t dist_y = 0, dist_uv = 0, tx_cache[TX_MODES] = { 0 };
  x->skip_encode = 0;
  ctx->skip = 0;
  xd->mode_info_context->mbmi.ref_frame[0] = INTRA_FRAME;
  if (bsize >= BLOCK_8X8) {
    if (rd_pick_intra_sby_mode(cpi, x, &rate_y, &rate_y_tokenonly,
                               &dist_y, &y_skip, bsize, tx_cache,
                               best_rd) >= best_rd) {
      *returnrate = INT_MAX;
      return;
    }
    rd_pick_intra_sbuv_mode(cpi, x, &rate_uv, &rate_uv_tokenonly,
                            &dist_uv, &uv_skip, bsize);
  } else {
    y_skip = 0;
    if (rd_pick_intra_sub_8x8_y_mode(cpi, x, &rate_y, &rate_y_tokenonly,
                                     &dist_y, best_rd) >= best_rd) {
      *returnrate = INT_MAX;
      return;
    }
    rd_pick_intra_sbuv_mode(cpi, x, &rate_uv, &rate_uv_tokenonly,
                            &dist_uv, &uv_skip, BLOCK_8X8);
  }

  if (y_skip && uv_skip) {
    *returnrate = rate_y + rate_uv - rate_y_tokenonly - rate_uv_tokenonly +
                  vp9_cost_bit(vp9_get_pred_prob_mbskip(cm, xd), 1);
    *returndist = dist_y + (dist_uv >> 2);
    vp9_zero(ctx->tx_rd_diff);
  } else {
    int i;
    *returnrate = rate_y + rate_uv +
        vp9_cost_bit(vp9_get_pred_prob_mbskip(cm, xd), 0);
    *returndist = dist_y + (dist_uv >> 2);
    if (cpi->sf.tx_size_search_method == USE_FULL_RD)
      for (i = 0; i < TX_MODES; i++)
        ctx->tx_rd_diff[i] = tx_cache[i] - tx_cache[cm->tx_mode];
  }

  ctx->mic = *xd->mode_info_context;
}

int64_t vp9_rd_pick_inter_mode_sb(VP9_COMP *cpi, MACROBLOCK *x,
                                  int mi_row, int mi_col,
                                  int *returnrate,
                                  int64_t *returndistortion,
                                  BLOCK_SIZE_TYPE bsize,
                                  PICK_MODE_CONTEXT *ctx,
                                  int64_t best_rd_so_far) {
  VP9_COMMON *cm = &cpi->common;
  MACROBLOCKD *xd = &x->e_mbd;
  MB_MODE_INFO *mbmi = &xd->mode_info_context->mbmi;
  const struct segmentation *seg = &xd->seg;
  const BLOCK_SIZE_TYPE block_size = get_plane_block_size(bsize, &xd->plane[0]);
  MB_PREDICTION_MODE this_mode;
  MV_REFERENCE_FRAME ref_frame, second_ref_frame;
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
  int64_t best_rd = best_rd_so_far;
  int64_t best_yrd = best_rd_so_far;  // FIXME(rbultje) more precise
  int64_t best_tx_rd[TX_MODES];
  int64_t best_tx_diff[TX_MODES];
  int64_t best_pred_diff[NB_PREDICTION_TYPES];
  int64_t best_pred_rd[NB_PREDICTION_TYPES];
  int64_t best_filter_rd[VP9_SWITCHABLE_FILTERS + 1];
  int64_t best_filter_diff[VP9_SWITCHABLE_FILTERS + 1];
  MB_MODE_INFO best_mbmode;
  int j;
  int mode_index, best_mode_index = 0;
  unsigned int ref_costs_single[MAX_REF_FRAMES], ref_costs_comp[MAX_REF_FRAMES];
  vp9_prob comp_mode_p;
  int64_t best_intra_rd = INT64_MAX;
  int64_t best_inter_rd = INT64_MAX;
  MB_PREDICTION_MODE best_intra_mode = DC_PRED;
  // MB_PREDICTION_MODE best_inter_mode = ZEROMV;
  MV_REFERENCE_FRAME best_inter_ref_frame = LAST_FRAME;
  INTERPOLATIONFILTERTYPE tmp_best_filter = SWITCHABLE;
  int rate_uv_intra[TX_SIZES], rate_uv_tokenonly[TX_SIZES];
  int64_t dist_uv[TX_SIZES];
  int skip_uv[TX_SIZES];
  MB_PREDICTION_MODE mode_uv[TX_SIZES];
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
  int best_skip2 = 0;

  x->skip_encode = (cpi->sf.skip_encode_frame &&
                    xd->q_index < QIDX_SKIP_THRESH);

  for (i = 0; i < 4; i++) {
    int j;
    for (j = 0; j < MAX_REF_FRAMES; j++)
      seg_mvs[i][j].as_int = INVALID_MV;
  }
  // Everywhere the flag is set the error is much higher than its neighbors.
  ctx->frames_with_high_error = 0;
  ctx->modes_with_high_error = 0;

  estimate_ref_frame_costs(cpi, segment_id, ref_costs_single, ref_costs_comp,
                           &comp_mode_p);
  vpx_memset(&best_mbmode, 0, sizeof(best_mbmode));
  vpx_memset(&single_newmv, 0, sizeof(single_newmv));

  for (i = 0; i < NB_PREDICTION_TYPES; ++i)
    best_pred_rd[i] = INT64_MAX;
  for (i = 0; i < TX_MODES; i++)
    best_tx_rd[i] = INT64_MAX;
  for (i = 0; i <= VP9_SWITCHABLE_FILTERS; i++)
    best_filter_rd[i] = INT64_MAX;
  for (i = 0; i < TX_SIZES; i++)
    rate_uv_intra[i] = INT_MAX;

  *returnrate = INT_MAX;

  // Create a mask set to 1 for each reference frame used by a smaller
  // resolution.
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

  for (mode_index = 0; mode_index < MAX_MODES; ++mode_index) {
    int mode_excluded = 0;
    int64_t this_rd = INT64_MAX;
    int disable_skip = 0;
    int compmode_cost = 0;
    int rate2 = 0, rate_y = 0, rate_uv = 0;
    int64_t distortion2 = 0, distortion_y = 0, distortion_uv = 0;
    int skippable = 0;
    int64_t tx_cache[TX_MODES];
    int i;
    int this_skip2 = 0;
    int64_t total_sse = INT_MAX;
    int early_term = 0;

    for (i = 0; i < TX_MODES; ++i)
      tx_cache[i] = INT64_MAX;

    x->skip = 0;
    this_mode = vp9_mode_order[mode_index].mode;
    ref_frame = vp9_mode_order[mode_index].ref_frame;
    second_ref_frame = vp9_mode_order[mode_index].second_ref_frame;

    // Skip modes that have been masked off but always consider first mode.
    if (mode_index && (bsize > cpi->sf.unused_mode_skip_lvl) &&
         (cpi->unused_mode_skip_mask & (1 << mode_index)) )
      continue;

    // Skip if the current reference frame has been masked off
    if (cpi->sf.reference_masking && !cpi->set_ref_frame_mask &&
        (cpi->ref_frame_mask & (1 << ref_frame)))
      continue;

    // Test best rd so far against threshold for trying this mode.
    if ((best_rd < ((cpi->rd_threshes[bsize][mode_index] *
                     cpi->rd_thresh_freq_fact[bsize][mode_index]) >> 4)) ||
        cpi->rd_threshes[bsize][mode_index] == INT_MAX)
      continue;

    // Do not allow compound prediction if the segment level reference
    // frame feature is in use as in this case there can only be one reference.
    if ((second_ref_frame > INTRA_FRAME) &&
         vp9_segfeature_active(seg, segment_id, SEG_LVL_REF_FRAME))
      continue;

    // Skip some checking based on small partitions' result.
    if (x->fast_ms > 1 && !ref_frame)
      continue;
    if (x->fast_ms > 2 && ref_frame != x->subblock_ref)
      continue;

    if (cpi->sf.use_avoid_tested_higherror && bsize >= BLOCK_8X8) {
      if (!(ref_frame_mask & (1 << ref_frame))) {
        continue;
      }
      if (!(mode_mask & (1 << this_mode))) {
        continue;
      }
      if (second_ref_frame != NONE
          && !(ref_frame_mask & (1 << second_ref_frame))) {
        continue;
      }
    }

    mbmi->ref_frame[0] = ref_frame;
    mbmi->ref_frame[1] = second_ref_frame;

    if (!(ref_frame == INTRA_FRAME
        || (cpi->ref_frame_flags & flag_list[ref_frame]))) {
      continue;
    }
    if (!(second_ref_frame == NONE
        || (cpi->ref_frame_flags & flag_list[second_ref_frame]))) {
      continue;
    }

    comp_pred = second_ref_frame > INTRA_FRAME;
    if (comp_pred) {
      if (cpi->sf.mode_search_skip_flags & FLAG_SKIP_COMP_BESTINTRA)
        if (vp9_mode_order[best_mode_index].ref_frame == INTRA_FRAME)
          continue;
      if (cpi->sf.mode_search_skip_flags & FLAG_SKIP_COMP_REFMISMATCH)
        if (ref_frame != best_inter_ref_frame &&
            second_ref_frame != best_inter_ref_frame)
          continue;
    }
    // TODO(jingning, jkoleszar): scaling reference frame not supported for
    // SPLITMV.
    if (ref_frame > 0 &&
        (scale_factor[ref_frame].x_scale_fp != VP9_REF_NO_SCALE ||
         scale_factor[ref_frame].y_scale_fp != VP9_REF_NO_SCALE) &&
        this_mode == SPLITMV)
      continue;

    if (second_ref_frame > 0 &&
        (scale_factor[second_ref_frame].x_scale_fp != VP9_REF_NO_SCALE ||
         scale_factor[second_ref_frame].y_scale_fp != VP9_REF_NO_SCALE) &&
        this_mode == SPLITMV)
      continue;

    set_scale_factors(xd, ref_frame, second_ref_frame, scale_factor);
    mbmi->mode = this_mode;
    mbmi->uv_mode = DC_PRED;

    // Evaluate all sub-pel filters irrespective of whether we can use
    // them for this frame.
    mbmi->interp_filter = cm->mcomp_filter_type;
    vp9_setup_interp_filters(xd, mbmi->interp_filter, &cpi->common);

    if (bsize >= BLOCK_8X8 &&
        (this_mode == I4X4_PRED || this_mode == SPLITMV))
      continue;
    if (bsize < BLOCK_8X8 &&
        !(this_mode == I4X4_PRED || this_mode == SPLITMV))
      continue;

    if (comp_pred) {
      if (!(cpi->ref_frame_flags & flag_list[second_ref_frame]))
        continue;
      set_scale_factors(xd, ref_frame, second_ref_frame, scale_factor);

      mode_excluded = mode_excluded
                         ? mode_excluded
                         : cm->comp_pred_mode == SINGLE_PREDICTION_ONLY;
    } else {
      if (ref_frame != INTRA_FRAME && second_ref_frame != INTRA_FRAME) {
        mode_excluded =
            mode_excluded ?
                mode_excluded : cm->comp_pred_mode == COMP_PREDICTION_ONLY;
      }
    }

    // Select prediction reference frames.
    for (i = 0; i < MAX_MB_PLANE; i++) {
      xd->plane[i].pre[0] = yv12_mb[ref_frame][i];
      if (comp_pred)
        xd->plane[i].pre[1] = yv12_mb[second_ref_frame][i];
    }

    // If the segment reference frame feature is enabled....
    // then do nothing if the current ref frame is not allowed..
    if (vp9_segfeature_active(seg, segment_id, SEG_LVL_REF_FRAME) &&
        vp9_get_segdata(seg, segment_id, SEG_LVL_REF_FRAME) !=
            (int)ref_frame) {
      continue;
    // If the segment skip feature is enabled....
    // then do nothing if the current mode is not allowed..
    } else if (vp9_segfeature_active(seg, segment_id, SEG_LVL_SKIP) &&
               (this_mode != ZEROMV && ref_frame != INTRA_FRAME)) {
      continue;
    // Disable this drop out case if the ref frame
    // segment level feature is enabled for this segment. This is to
    // prevent the possibility that we end up unable to pick any mode.
    } else if (!vp9_segfeature_active(seg, segment_id,
                                      SEG_LVL_REF_FRAME)) {
      // Only consider ZEROMV/ALTREF_FRAME for alt ref frame,
      // unless ARNR filtering is enabled in which case we want
      // an unfiltered alternative. We allow near/nearest as well
      // because they may result in zero-zero MVs but be cheaper.
      if (cpi->is_src_frame_alt_ref && (cpi->oxcf.arnr_max_frames == 0)) {
        if ((this_mode != ZEROMV &&
             !(this_mode == NEARMV &&
               frame_mv[NEARMV][ALTREF_FRAME].as_int == 0) &&
             !(this_mode == NEARESTMV &&
               frame_mv[NEARESTMV][ALTREF_FRAME].as_int == 0)) ||
            ref_frame != ALTREF_FRAME) {
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

      /*
      if ((cpi->sf.mode_search_skip_flags & FLAG_SKIP_INTRA_BESTINTER) &&
          (vp9_mode_order[best_mode_index].ref_frame > INTRA_FRAME))
        continue;
        */

      // I4X4_PRED is only considered for block sizes less than 8x8.
      mbmi->txfm_size = TX_4X4;
      if (rd_pick_intra_sub_8x8_y_mode(cpi, x, &rate, &rate_y,
                                       &distortion_y, best_rd) >= best_rd)
        continue;
      rate2 += rate;
      rate2 += intra_cost_penalty;
      distortion2 += distortion_y;

      if (rate_uv_intra[TX_4X4] == INT_MAX) {
        choose_intra_uv_mode(cpi, bsize, &rate_uv_intra[TX_4X4],
                             &rate_uv_tokenonly[TX_4X4],
                             &dist_uv[TX_4X4], &skip_uv[TX_4X4],
                             &mode_uv[TX_4X4]);
      }
      rate2 += rate_uv_intra[TX_4X4];
      rate_uv = rate_uv_tokenonly[TX_4X4];
      distortion2 += dist_uv[TX_4X4];
      distortion_uv = dist_uv[TX_4X4];
      mbmi->uv_mode = mode_uv[TX_4X4];
      tx_cache[ONLY_4X4] = RDCOST(x->rdmult, x->rddiv, rate2, distortion2);
      for (i = 0; i < TX_MODES; ++i)
        tx_cache[i] = tx_cache[ONLY_4X4];
    } else if (ref_frame == INTRA_FRAME) {
      TX_SIZE uv_tx;
      // Disable intra modes other than DC_PRED for blocks with low variance
      // Threshold for intra skipping based on source variance
      // TODO(debargha): Specialize the threshold for super block sizes
      static const int skip_intra_var_thresh[BLOCK_SIZES] = {
        64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64,
      };
      if ((cpi->sf.mode_search_skip_flags & FLAG_SKIP_INTRA_LOWVAR) &&
          this_mode != DC_PRED &&
          x->source_variance < skip_intra_var_thresh[mbmi->sb_type])
        continue;
      // Only search the oblique modes if the best so far is
      // one of the neighboring directional modes
      if ((cpi->sf.mode_search_skip_flags & FLAG_SKIP_INTRA_BESTINTER) &&
          (this_mode >= D45_PRED && this_mode <= TM_PRED)) {
        if (vp9_mode_order[best_mode_index].ref_frame > INTRA_FRAME)
          continue;
      }
      if (cpi->sf.mode_search_skip_flags & FLAG_SKIP_INTRA_DIRMISMATCH) {
        if (conditional_skipintra(mbmi->mode, best_intra_mode))
            continue;
      }
      super_block_yrd(cpi, x, &rate_y, &distortion_y, &skippable, NULL,
                      bsize, tx_cache, best_rd);

      if (rate_y == INT_MAX)
        continue;

      uv_tx = MIN(mbmi->txfm_size, max_uv_txsize_lookup[bsize]);
      if (rate_uv_intra[uv_tx] == INT_MAX) {
        choose_intra_uv_mode(cpi, bsize, &rate_uv_intra[uv_tx],
                             &rate_uv_tokenonly[uv_tx],
                             &dist_uv[uv_tx], &skip_uv[uv_tx],
                             &mode_uv[uv_tx]);
      }

      rate_uv = rate_uv_tokenonly[uv_tx];
      distortion_uv = dist_uv[uv_tx];
      skippable = skippable && skip_uv[uv_tx];
      mbmi->uv_mode = mode_uv[uv_tx];

      rate2 = rate_y + x->mbmode_cost[mbmi->mode] + rate_uv_intra[uv_tx];
      if (mbmi->mode != DC_PRED && mbmi->mode != TM_PRED)
        rate2 += intra_cost_penalty;
      distortion2 = distortion_y + distortion_uv;
    } else if (this_mode == SPLITMV) {
      const int is_comp_pred = second_ref_frame > 0;
      int rate;
      int64_t distortion;
      int64_t this_rd_thresh;
      int64_t tmp_rd, tmp_best_rd = INT64_MAX, tmp_best_rdu = INT64_MAX;
      int tmp_best_rate = INT_MAX, tmp_best_ratey = INT_MAX;
      int64_t tmp_best_distortion = INT_MAX, tmp_best_sse, uv_sse;
      int tmp_best_skippable = 0;
      int switchable_filter_index;
      int_mv *second_ref = is_comp_pred ?
          &mbmi->ref_mvs[second_ref_frame][0] : NULL;
      union b_mode_info tmp_best_bmodes[16];
      MB_MODE_INFO tmp_best_mbmode;
      PARTITION_INFO tmp_best_partition;
      BEST_SEG_INFO bsi[VP9_SWITCHABLE_FILTERS];
      int pred_exists = 0;
      int uv_skippable;
      if (is_comp_pred) {
        if (cpi->sf.mode_search_skip_flags & FLAG_SKIP_COMP_BESTINTRA)
          if (vp9_mode_order[best_mode_index].ref_frame == INTRA_FRAME)
            continue;
        if (cpi->sf.mode_search_skip_flags & FLAG_SKIP_COMP_REFMISMATCH)
          if (ref_frame != best_inter_ref_frame &&
              second_ref_frame != best_inter_ref_frame)
            continue;
      }

      this_rd_thresh = (ref_frame == LAST_FRAME) ?
          cpi->rd_threshes[bsize][THR_NEWMV] :
          cpi->rd_threshes[bsize][THR_NEWA];
      this_rd_thresh = (ref_frame == GOLDEN_FRAME) ?
          cpi->rd_threshes[bsize][THR_NEWG] : this_rd_thresh;
      xd->mode_info_context->mbmi.txfm_size = TX_4X4;

      cpi->rd_filter_cache[VP9_SWITCHABLE_FILTERS] = INT64_MAX;
      for (switchable_filter_index = 0;
           switchable_filter_index < VP9_SWITCHABLE_FILTERS;
           ++switchable_filter_index) {
        int newbest, rs;
        int64_t rs_rd;
        mbmi->interp_filter = switchable_filter_index;
        vp9_setup_interp_filters(xd, mbmi->interp_filter, &cpi->common);

        tmp_rd = rd_pick_best_mbsegmentation(cpi, x,
                     &mbmi->ref_mvs[ref_frame][0],
                     second_ref,
                     best_yrd,
                     &rate, &rate_y, &distortion,
                     &skippable, &total_sse,
                     (int)this_rd_thresh, seg_mvs,
                     bsi, switchable_filter_index,
                     mi_row, mi_col);

        if (tmp_rd == INT64_MAX)
          continue;
        cpi->rd_filter_cache[switchable_filter_index] = tmp_rd;
        rs = get_switchable_rate(x);
        rs_rd = RDCOST(x->rdmult, x->rddiv, rs, 0);
        cpi->rd_filter_cache[VP9_SWITCHABLE_FILTERS] =
            MIN(cpi->rd_filter_cache[VP9_SWITCHABLE_FILTERS], tmp_rd + rs_rd);
        if (cm->mcomp_filter_type == SWITCHABLE)
          tmp_rd += rs_rd;

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
          tmp_best_sse = total_sse;
          tmp_best_skippable = skippable;
          tmp_best_mbmode = *mbmi;
          tmp_best_partition = *x->partition_info;
          for (i = 0; i < 4; i++)
            tmp_best_bmodes[i] = xd->mode_info_context->bmi[i];
          pred_exists = 1;
          if (switchable_filter_index == 0 &&
              cpi->sf.use_rd_breakout &&
              best_rd < INT64_MAX) {
            if (tmp_best_rdu / 2 > best_rd) {
              // skip searching the other filters if the first is
              // already substantially larger than the best so far
              tmp_best_filter = mbmi->interp_filter;
              tmp_best_rdu = INT64_MAX;
              break;
            }
          }
        }
      }  // switchable_filter_index loop

      if (tmp_best_rdu == INT64_MAX)
        continue;

      mbmi->interp_filter = (cm->mcomp_filter_type == SWITCHABLE ?
                             tmp_best_filter : cm->mcomp_filter_type);
      vp9_setup_interp_filters(xd, mbmi->interp_filter, &cpi->common);
      if (!pred_exists) {
        // Handles the special case when a filter that is not in the
        // switchable list (bilinear, 6-tap) is indicated at the frame level
        tmp_rd = rd_pick_best_mbsegmentation(cpi, x,
                     &mbmi->ref_mvs[ref_frame][0],
                     second_ref,
                     best_yrd,
                     &rate, &rate_y, &distortion,
                     &skippable, &total_sse,
                     (int)this_rd_thresh, seg_mvs,
                     bsi, 0,
                     mi_row, mi_col);
        if (tmp_rd == INT64_MAX)
          continue;
      } else {
        if (cpi->common.mcomp_filter_type == SWITCHABLE) {
          int rs = get_switchable_rate(x);
          tmp_best_rdu -= RDCOST(x->rdmult, x->rddiv, rs, 0);
        }
        tmp_rd = tmp_best_rdu;
        total_sse = tmp_best_sse;
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
        rate2 += get_switchable_rate(x);

      if (!mode_excluded) {
        if (is_comp_pred)
          mode_excluded = cpi->common.comp_pred_mode == SINGLE_PREDICTION_ONLY;
        else
          mode_excluded = cpi->common.comp_pred_mode == COMP_PREDICTION_ONLY;
      }
      compmode_cost = vp9_cost_bit(comp_mode_p, is_comp_pred);

      if (RDCOST(x->rdmult, x->rddiv, rate2, distortion2) <
          best_rd) {
        // If even the 'Y' rd value of split is higher than best so far
        // then dont bother looking at UV
        vp9_build_inter_predictors_sbuv(&x->e_mbd, mi_row, mi_col,
                                        BLOCK_8X8);
        vp9_subtract_sbuv(x, BLOCK_8X8);
        super_block_uvrd_for_txfm(cm, x, &rate_uv, &distortion_uv,
                                  &uv_skippable, &uv_sse, BLOCK_8X8, TX_4X4);
        rate2 += rate_uv;
        distortion2 += distortion_uv;
        skippable = skippable && uv_skippable;
        total_sse += uv_sse;

        tx_cache[ONLY_4X4] = RDCOST(x->rdmult, x->rddiv, rate2, distortion2);
        for (i = 0; i < TX_MODES; ++i)
          tx_cache[i] = tx_cache[ONLY_4X4];
      }
    } else {
      compmode_cost = vp9_cost_bit(comp_mode_p, second_ref_frame > INTRA_FRAME);
      this_rd = handle_inter_mode(cpi, x, bsize,
                                  tx_cache,
                                  &rate2, &distortion2, &skippable,
                                  &rate_y, &distortion_y,
                                  &rate_uv, &distortion_uv,
                                  &mode_excluded, &disable_skip,
                                  &tmp_best_filter, frame_mv,
                                  mi_row, mi_col,
                                  single_newmv, &total_sse, best_rd);
      if (this_rd == INT64_MAX)
        continue;
    }

    if (cpi->common.comp_pred_mode == HYBRID_PREDICTION) {
      rate2 += compmode_cost;
    }

    // Estimate the reference frame signaling cost and add it
    // to the rolling cost variable.
    if (second_ref_frame > INTRA_FRAME) {
      rate2 += ref_costs_comp[ref_frame];
    } else {
      rate2 += ref_costs_single[ref_frame];
    }

    if (!disable_skip) {
      // Test for the condition where skip block will be activated
      // because there are no non zero coefficients and make any
      // necessary adjustment for rate. Ignore if skip is coded at
      // segment level as the cost wont have been added in.
      // Is Mb level skip allowed (i.e. not coded at segment level).
      const int mb_skip_allowed = !vp9_segfeature_active(seg, segment_id,
                                                         SEG_LVL_SKIP);

      if (skippable && bsize >= BLOCK_8X8) {
        // Back out the coefficient coding costs
        rate2 -= (rate_y + rate_uv);
        // for best yrd calculation
        rate_uv = 0;

        if (mb_skip_allowed) {
          int prob_skip_cost;

          // Cost the skip mb case
          vp9_prob skip_prob =
            vp9_get_pred_prob_mbskip(cm, xd);

          if (skip_prob) {
            prob_skip_cost = vp9_cost_bit(skip_prob, 1);
            rate2 += prob_skip_cost;
          }
        }
      } else if (mb_skip_allowed && ref_frame != INTRA_FRAME && !xd->lossless) {
        if (RDCOST(x->rdmult, x->rddiv, rate_y + rate_uv, distortion2) <
            RDCOST(x->rdmult, x->rddiv, 0, total_sse)) {
          // Add in the cost of the no skip flag.
          int prob_skip_cost = vp9_cost_bit(vp9_get_pred_prob_mbskip(cm, xd),
                                            0);
          rate2 += prob_skip_cost;
        } else {
          // FIXME(rbultje) make this work for splitmv also
          int prob_skip_cost = vp9_cost_bit(vp9_get_pred_prob_mbskip(cm, xd),
                                            1);
          rate2 += prob_skip_cost;
          distortion2 = total_sse;
          assert(total_sse >= 0);
          rate2 -= (rate_y + rate_uv);
          rate_y = 0;
          rate_uv = 0;
          this_skip2 = 1;
        }
      } else if (mb_skip_allowed) {
        // Add in the cost of the no skip flag.
        int prob_skip_cost = vp9_cost_bit(vp9_get_pred_prob_mbskip(cm, xd),
                                          0);
        rate2 += prob_skip_cost;
      }

      // Calculate the final RD estimate for this mode.
      this_rd = RDCOST(x->rdmult, x->rddiv, rate2, distortion2);
    }

    // Keep record of best intra rd
    if (xd->mode_info_context->mbmi.ref_frame[0] == INTRA_FRAME &&
        is_intra_mode(xd->mode_info_context->mbmi.mode) &&
        this_rd < best_intra_rd) {
      best_intra_rd = this_rd;
      best_intra_mode = xd->mode_info_context->mbmi.mode;
    }
    // Keep record of best inter rd with single reference
    if (xd->mode_info_context->mbmi.ref_frame[0] > INTRA_FRAME &&
        xd->mode_info_context->mbmi.ref_frame[1] == NONE &&
        !mode_excluded &&
        this_rd < best_inter_rd) {
      best_inter_rd = this_rd;
      best_inter_ref_frame = ref_frame;
      // best_inter_mode = xd->mode_info_context->mbmi.mode;
    }

    if (!disable_skip && ref_frame == INTRA_FRAME) {
      for (i = 0; i < NB_PREDICTION_TYPES; ++i)
        best_pred_rd[i] = MIN(best_pred_rd[i], this_rd);
      for (i = 0; i <= VP9_SWITCHABLE_FILTERS; i++)
        best_filter_rd[i] = MIN(best_filter_rd[i], this_rd);
    }

    if (this_mode != I4X4_PRED && this_mode != SPLITMV) {
      // Store the respective mode distortions for later use.
      if (mode_distortions[this_mode] == -1
          || distortion2 < mode_distortions[this_mode]) {
        mode_distortions[this_mode] = distortion2;
      }
      if (frame_distortions[ref_frame] == -1
          || distortion2 < frame_distortions[ref_frame]) {
        frame_distortions[ref_frame] = distortion2;
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
        best_yrd = best_rd -
                   RDCOST(x->rdmult, x->rddiv, rate_uv, distortion_uv);
        best_mbmode = *mbmi;
        best_skip2 = this_skip2;
        best_partition = *x->partition_info;

        if (this_mode == I4X4_PRED || this_mode == SPLITMV)
          for (i = 0; i < 4; i++)
            best_bmodes[i] = xd->mode_info_context->bmi[i];

        // TODO(debargha): enhance this test with a better distortion prediction
        // based on qp, activity mask and history
        if (cpi->sf.mode_search_skip_flags & FLAG_EARLY_TERMINATE) {
          const int qstep = xd->plane[0].dequant[1];
          // TODO(debargha): Enhance this by specializing for each mode_index
          int scale = 4;
          if (x->source_variance < UINT_MAX) {
            const int var_adjust = (x->source_variance < 16);
            scale -= var_adjust;
          }
          if (ref_frame > INTRA_FRAME &&
              distortion2 * scale < qstep * qstep) {
            early_term = 1;
          }
        }
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
    if (!disable_skip && ref_frame != INTRA_FRAME) {
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

      if (second_ref_frame <= INTRA_FRAME &&
          single_rd < best_pred_rd[SINGLE_PREDICTION_ONLY]) {
        best_pred_rd[SINGLE_PREDICTION_ONLY] = single_rd;
      } else if (second_ref_frame > INTRA_FRAME &&
                 single_rd < best_pred_rd[COMP_PREDICTION_ONLY]) {
        best_pred_rd[COMP_PREDICTION_ONLY] = single_rd;
      }
      if (hybrid_rd < best_pred_rd[HYBRID_PREDICTION])
        best_pred_rd[HYBRID_PREDICTION] = hybrid_rd;
    }

    /* keep record of best filter type */
    if (!mode_excluded && !disable_skip && ref_frame != INTRA_FRAME &&
        cm->mcomp_filter_type != BILINEAR) {
      int64_t ref = cpi->rd_filter_cache[cm->mcomp_filter_type == SWITCHABLE ?
                              VP9_SWITCHABLE_FILTERS : cm->mcomp_filter_type];
      for (i = 0; i <= VP9_SWITCHABLE_FILTERS; i++) {
        int64_t adj_rd;
        // In cases of poor prediction, filter_cache[] can contain really big
        // values, which actually are bigger than this_rd itself. This can
        // cause negative best_filter_rd[] values, which is obviously silly.
        // Therefore, if filter_cache < ref, we do an adjusted calculation.
        if (cpi->rd_filter_cache[i] >= ref)
          adj_rd = this_rd + cpi->rd_filter_cache[i] - ref;
        else  // FIXME(rbultje) do this for comppred also
          adj_rd = this_rd - (ref - cpi->rd_filter_cache[i]) * this_rd / ref;
        best_filter_rd[i] = MIN(best_filter_rd[i], adj_rd);
      }
    }

    /* keep record of best txfm size */
    if (bsize < BLOCK_32X32) {
      if (bsize < BLOCK_16X16) {
        if (this_mode == SPLITMV || this_mode == I4X4_PRED)
          tx_cache[ALLOW_8X8] = tx_cache[ONLY_4X4];
        tx_cache[ALLOW_16X16] = tx_cache[ALLOW_8X8];
      }
      tx_cache[ALLOW_32X32] = tx_cache[ALLOW_16X16];
    }
    if (!mode_excluded && this_rd != INT64_MAX) {
      for (i = 0; i < TX_MODES; i++) {
        int64_t adj_rd = INT64_MAX;
        if (this_mode != I4X4_PRED) {
          adj_rd = this_rd + tx_cache[i] - tx_cache[cm->tx_mode];
        } else {
          adj_rd = this_rd;
        }

        if (adj_rd < best_tx_rd[i])
          best_tx_rd[i] = adj_rd;
      }
    }

    if (early_term)
      break;

    if (x->skip && !comp_pred)
      break;
  }

  if (best_rd >= best_rd_so_far)
    return INT64_MAX;

  // If we used an estimate for the uv intra rd in the loop above...
  if (cpi->sf.use_uv_intra_rd_estimate) {
    // Do Intra UV best rd mode selection if best mode choice above was intra.
    if (vp9_mode_order[best_mode_index].ref_frame == INTRA_FRAME) {
      TX_SIZE uv_tx_size = get_uv_tx_size(mbmi);
      rd_pick_intra_sbuv_mode(cpi, x, &rate_uv_intra[uv_tx_size],
                              &rate_uv_tokenonly[uv_tx_size],
                              &dist_uv[uv_tx_size],
                              &skip_uv[uv_tx_size],
                              bsize < BLOCK_8X8 ? BLOCK_8X8 : bsize);
    }
  }

  // If indicated then mark the index of the chosen mode to be inspected at
  // other block sizes.
  if (bsize <= cpi->sf.unused_mode_skip_lvl) {
    cpi->unused_mode_skip_mask = cpi->unused_mode_skip_mask &
                                 (~((int64_t)1 << best_mode_index));
  }

  // If we are using reference masking and the set mask flag is set then
  // create the reference frame mask.
  if (cpi->sf.reference_masking && cpi->set_ref_frame_mask)
    cpi->ref_frame_mask = ~(1 << vp9_mode_order[best_mode_index].ref_frame);

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

  if (best_rd == INT64_MAX && bsize < BLOCK_8X8) {
    *returnrate = INT_MAX;
    *returndistortion = INT_MAX;
    return best_rd;
  }

  assert((cm->mcomp_filter_type == SWITCHABLE) ||
         (cm->mcomp_filter_type == best_mbmode.interp_filter) ||
         (best_mbmode.ref_frame[0] == INTRA_FRAME));

  // Updating rd_thresh_freq_fact[] here means that the different
  // partition/block sizes are handled independently based on the best
  // choice for the current partition. It may well be better to keep a scaled
  // best rd so far value and update rd_thresh_freq_fact based on the mode/size
  // combination that wins out.
  if (cpi->sf.adaptive_rd_thresh) {
    for (mode_index = 0; mode_index < MAX_MODES; ++mode_index) {
      if (mode_index == best_mode_index) {
        cpi->rd_thresh_freq_fact[bsize][mode_index] = BASE_RD_THRESH_FREQ_FACT;
      } else {
        cpi->rd_thresh_freq_fact[bsize][mode_index] += MAX_RD_THRESH_FREQ_INC;
        if (cpi->rd_thresh_freq_fact[bsize][mode_index] >
            (cpi->sf.adaptive_rd_thresh * MAX_RD_THRESH_FREQ_FACT)) {
          cpi->rd_thresh_freq_fact[bsize][mode_index] =
            cpi->sf.adaptive_rd_thresh * MAX_RD_THRESH_FREQ_FACT;
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

  // macroblock modes
  *mbmi = best_mbmode;
  x->skip |= best_skip2;
  if (best_mbmode.ref_frame[0] == INTRA_FRAME &&
      best_mbmode.sb_type < BLOCK_8X8) {
    for (i = 0; i < 4; i++)
      xd->mode_info_context->bmi[i].as_mode = best_bmodes[i].as_mode;
  }

  if (best_mbmode.ref_frame[0] != INTRA_FRAME &&
      best_mbmode.sb_type < BLOCK_8X8) {
    for (i = 0; i < 4; i++)
      xd->mode_info_context->bmi[i].as_mv[0].as_int =
          best_bmodes[i].as_mv[0].as_int;

    if (mbmi->ref_frame[1] > 0)
      for (i = 0; i < 4; i++)
        xd->mode_info_context->bmi[i].as_mv[1].as_int =
            best_bmodes[i].as_mv[1].as_int;

    *x->partition_info = best_partition;

    mbmi->mv[0].as_int = xd->mode_info_context->bmi[3].as_mv[0].as_int;
    mbmi->mv[1].as_int = xd->mode_info_context->bmi[3].as_mv[1].as_int;
  }

  for (i = 0; i < NB_PREDICTION_TYPES; ++i) {
    if (best_pred_rd[i] == INT64_MAX)
      best_pred_diff[i] = INT_MIN;
    else
      best_pred_diff[i] = best_rd - best_pred_rd[i];
  }

  if (!x->skip) {
    for (i = 0; i <= VP9_SWITCHABLE_FILTERS; i++) {
      if (best_filter_rd[i] == INT64_MAX)
        best_filter_diff[i] = 0;
      else
        best_filter_diff[i] = best_rd - best_filter_rd[i];
    }
    if (cm->mcomp_filter_type == SWITCHABLE)
      assert(best_filter_diff[VP9_SWITCHABLE_FILTERS] == 0);
  } else {
    vpx_memset(best_filter_diff, 0, sizeof(best_filter_diff));
  }

  if (!x->skip) {
    for (i = 0; i < TX_MODES; i++) {
      if (best_tx_rd[i] == INT64_MAX)
        best_tx_diff[i] = 0;
      else
        best_tx_diff[i] = best_rd - best_tx_rd[i];
    }
  } else {
    vpx_memset(best_tx_diff, 0, sizeof(best_tx_diff));
  }

  set_scale_factors(xd, mbmi->ref_frame[0], mbmi->ref_frame[1],
                    scale_factor);
  store_coding_context(x, ctx, best_mode_index,
                       &best_partition,
                       &mbmi->ref_mvs[mbmi->ref_frame[0]][0],
                       &mbmi->ref_mvs[mbmi->ref_frame[1] < 0 ? 0 :
                                      mbmi->ref_frame[1]][0],
                       best_pred_diff, best_tx_diff, best_filter_diff);

  return best_rd;
}
