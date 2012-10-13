/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include "vpx_ports/config.h"
#include "encodemb.h"
#include "vp8/common/reconinter.h"
#include "quantize.h"
#include "tokenize.h"
#include "vp8/common/invtrans.h"
#include "vp8/common/recon.h"
#include "vp8/common/reconintra.h"
#include "dct.h"
#include "vpx_mem/vpx_mem.h"
#include "rdopt.h"
#include "vp8/common/systemdependent.h"

#if CONFIG_RUNTIME_CPU_DETECT
#define IF_RTCD(x) (x)
#else
#define IF_RTCD(x) NULL
#endif

#ifdef ENC_DEBUG
extern int enc_debug;
#endif

void vp8_subtract_b_c(BLOCK *be, BLOCKD *bd, int pitch) {
  unsigned char *src_ptr = (*(be->base_src) + be->src);
  short *diff_ptr = be->src_diff;
  unsigned char *pred_ptr = bd->predictor;
  int src_stride = be->src_stride;

  int r, c;

  for (r = 0; r < 4; r++) {
    for (c = 0; c < 4; c++) {
      diff_ptr[c] = src_ptr[c] - pred_ptr[c];
    }

    diff_ptr += pitch;
    pred_ptr += pitch;
    src_ptr  += src_stride;
  }
}

void vp8_subtract_4b_c(BLOCK *be, BLOCKD *bd, int pitch) {
  unsigned char *src_ptr = (*(be->base_src) + be->src);
  short *diff_ptr = be->src_diff;
  unsigned char *pred_ptr = bd->predictor;
  int src_stride = be->src_stride;
  int r, c;

  for (r = 0; r < 8; r++) {
    for (c = 0; c < 8; c++) {
      diff_ptr[c] = src_ptr[c] - pred_ptr[c];
    }
    diff_ptr += pitch;
    pred_ptr += pitch;
    src_ptr  += src_stride;
  }
}

void vp8_subtract_mbuv_s_c(short *diff, const unsigned char *usrc,
                           const unsigned char *vsrc, int src_stride,
                           const unsigned char *upred,
                           const unsigned char *vpred, int dst_stride) {
  short *udiff = diff + 256;
  short *vdiff = diff + 320;
  int r, c;

  for (r = 0; r < 8; r++) {
    for (c = 0; c < 8; c++) {
      udiff[c] = usrc[c] - upred[c];
    }

    udiff += 8;
    upred += dst_stride;
    usrc  += src_stride;
  }

  for (r = 0; r < 8; r++) {
    for (c = 0; c < 8; c++) {
      vdiff[c] = vsrc[c] - vpred[c];
    }

    vdiff += 8;
    vpred += dst_stride;
    vsrc  += src_stride;
  }
}

void vp8_subtract_mbuv_c(short *diff, unsigned char *usrc,
                         unsigned char *vsrc, unsigned char *pred, int stride) {
  unsigned char *upred = pred + 256;
  unsigned char *vpred = pred + 320;

  vp8_subtract_mbuv_s_c(diff, usrc, vsrc, stride, upred, vpred, 8);
}

void vp8_subtract_mby_s_c(short *diff, const unsigned char *src, int src_stride,
                          const unsigned char *pred, int dst_stride) {
  int r, c;

  for (r = 0; r < 16; r++) {
    for (c = 0; c < 16; c++) {
      diff[c] = src[c] - pred[c];
    }

    diff += 16;
    pred += dst_stride;
    src  += src_stride;
  }
}

void vp8_subtract_mby_c(short *diff, unsigned char *src,
                        unsigned char *pred, int stride) {
  vp8_subtract_mby_s_c(diff, src, stride, pred, 16);
}

static void vp8_subtract_mb(const VP8_ENCODER_RTCD *rtcd, MACROBLOCK *x) {
  BLOCK *b = &x->block[0];

  ENCODEMB_INVOKE(&rtcd->encodemb, submby)(x->src_diff, *(b->base_src), x->e_mbd.predictor, b->src_stride);
  ENCODEMB_INVOKE(&rtcd->encodemb, submbuv)(x->src_diff, x->src.u_buffer, x->src.v_buffer, x->e_mbd.predictor, x->src.uv_stride);
}

static void build_dcblock_4x4(MACROBLOCK *x) {
  short *src_diff_ptr = &x->src_diff[384];
  int i;

  for (i = 0; i < 16; i++) {
    src_diff_ptr[i] = x->coeff[i * 16];
  }
}

void vp8_transform_mby_4x4(MACROBLOCK *x) {
  int i;

  for (i = 0; i < 16; i += 2) {
    x->vp8_short_fdct8x4(&x->block[i].src_diff[0],
                         &x->block[i].coeff[0], 32);
  }

  if (x->e_mbd.mode_info_context->mbmi.mode != SPLITMV) {
    // build dc block from 16 y dc values
    build_dcblock_4x4(x);

    // do 2nd order transform on the dc block
    x->short_walsh4x4(&x->block[24].src_diff[0],
                      &x->block[24].coeff[0], 8);
  }
}

void vp8_transform_mbuv_4x4(MACROBLOCK *x) {
  int i;

  for (i = 16; i < 24; i += 2) {
    x->vp8_short_fdct8x4(&x->block[i].src_diff[0],
                         &x->block[i].coeff[0], 16);
  }
}

static void transform_mb_4x4(MACROBLOCK *x) {
  vp8_transform_mby_4x4(x);
  vp8_transform_mbuv_4x4(x);
}

void vp8_build_dcblock_8x8(MACROBLOCK *x) {
  int16_t *src_diff_ptr = x->block[24].src_diff;
  int i;

  for (i = 0; i < 16; i++) {
    src_diff_ptr[i] = 0;
  }
  src_diff_ptr[0] = x->coeff[0 * 16];
  src_diff_ptr[1] = x->coeff[4 * 16];
  src_diff_ptr[4] = x->coeff[8 * 16];
  src_diff_ptr[8] = x->coeff[12 * 16];
}

void vp8_transform_mby_8x8(MACROBLOCK *x) {
  int i;

  for (i = 0; i < 9; i += 8) {
    x->vp8_short_fdct8x8(&x->block[i].src_diff[0],
                         &x->block[i].coeff[0], 32);
  }
  for (i = 2; i < 11; i += 8) {
    x->vp8_short_fdct8x8(&x->block[i].src_diff[0],
                         &x->block[i + 2].coeff[0], 32);
  }

  if (x->e_mbd.mode_info_context->mbmi.mode != SPLITMV) {
    // build dc block from 2x2 y dc values
    vp8_build_dcblock_8x8(x);

    // do 2nd order transform on the dc block
    x->short_fhaar2x2(&x->block[24].src_diff[0],
                      &x->block[24].coeff[0], 8);
  }
}

void vp8_transform_mbuv_8x8(MACROBLOCK *x) {
  int i;

  for (i = 16; i < 24; i += 4) {
    x->vp8_short_fdct8x8(&x->block[i].src_diff[0],
                         &x->block[i].coeff[0], 16);
  }
}

void vp8_transform_mb_8x8(MACROBLOCK *x) {
  vp8_transform_mby_8x8(x);
  vp8_transform_mbuv_8x8(x);
}

void vp8_transform_mby_16x16(MACROBLOCK *x) {
  vp8_clear_system_state();
  x->vp8_short_fdct16x16(&x->block[0].src_diff[0],
                         &x->block[0].coeff[0], 32);
}

void vp8_transform_mb_16x16(MACROBLOCK *x) {
  vp8_transform_mby_16x16(x);
  vp8_transform_mbuv_8x8(x);
}

#define RDTRUNC(RM,DM,R,D) ( (128+(R)*(RM)) & 0xFF )
#define RDTRUNC_8x8(RM,DM,R,D) ( (128+(R)*(RM)) & 0xFF )
typedef struct vp8_token_state vp8_token_state;

struct vp8_token_state {
  int           rate;
  int           error;
  int           next;
  signed char   token;
  short         qc;
};

// TODO: experiments to find optimal multiple numbers
#define Y1_RD_MULT 4
#define UV_RD_MULT 2
#define Y2_RD_MULT 4

static const int plane_rd_mult[4] = {
  Y1_RD_MULT,
  Y2_RD_MULT,
  UV_RD_MULT,
  Y1_RD_MULT
};

#define UPDATE_RD_COST()\
{\
  rd_cost0 = RDCOST(rdmult, rddiv, rate0, error0);\
  rd_cost1 = RDCOST(rdmult, rddiv, rate1, error1);\
  if (rd_cost0 == rd_cost1) {\
    rd_cost0 = RDTRUNC(rdmult, rddiv, rate0, error0);\
    rd_cost1 = RDTRUNC(rdmult, rddiv, rate1, error1);\
  }\
}

void optimize_b(MACROBLOCK *mb, int i, int type,
                ENTROPY_CONTEXT *a, ENTROPY_CONTEXT *l,
                const VP8_ENCODER_RTCD *rtcd, int tx_type) {
  BLOCK *b;
  BLOCKD *d;
  vp8_token_state tokens[65][2];
  uint64_t best_mask[2];
  const short *dequant_ptr;
  const short *coeff_ptr;
  short *qcoeff_ptr;
  short *dqcoeff_ptr;
  int eob;
  int i0;
  int rc;
  int x;
  int sz = 0;
  int next;
  int rdmult;
  int rddiv;
  int final_eob;
  int64_t rd_cost0, rd_cost1;
  int rate0, rate1;
  int error0, error1;
  int t0, t1;
  int best;
  int band;
  int pt;
  int err_mult = plane_rd_mult[type];
  int default_eob;
  int const *scan, *bands;

  b = &mb->block[i];
  d = &mb->e_mbd.block[i];
  switch (tx_type) {
    default:
    case TX_4X4:
      scan = vp8_default_zig_zag1d;
      bands = vp8_coef_bands;
      default_eob = 16;
#if CONFIG_HYBRIDTRANSFORM
      // TODO: this isn't called (for intra4x4 modes), but will be left in
      // since it could be used later
      {
        int active_ht = (mb->q_index < ACTIVE_HT) &&
                        (mb->e_mbd.mode_info_context->mbmi.mode == B_PRED);

        if((type == PLANE_TYPE_Y_WITH_DC) && active_ht) {
          switch (d->bmi.as_mode.tx_type) {
            case ADST_DCT:
              scan = vp8_row_scan;
              break;

            case DCT_ADST:
              scan = vp8_col_scan;
              break;

            default:
              scan = vp8_default_zig_zag1d;
              break;
          }

        } else
          scan = vp8_default_zig_zag1d;
      }
#endif
      break;
    case TX_8X8:
      scan = vp8_default_zig_zag1d_8x8;
      bands = vp8_coef_bands_8x8;
      default_eob = 64;
      break;
  }

  dequant_ptr = d->dequant;
  coeff_ptr = b->coeff;
  qcoeff_ptr = d->qcoeff;
  dqcoeff_ptr = d->dqcoeff;
  i0 = !type;
  eob = d->eob;

  /* Now set up a Viterbi trellis to evaluate alternative roundings. */
  rdmult = mb->rdmult * err_mult;
  if (mb->e_mbd.mode_info_context->mbmi.ref_frame == INTRA_FRAME)
    rdmult = (rdmult * 9) >> 4;
  rddiv = mb->rddiv;
  best_mask[0] = best_mask[1] = 0;
  /* Initialize the sentinel node of the trellis. */
  tokens[eob][0].rate = 0;
  tokens[eob][0].error = 0;
  tokens[eob][0].next = default_eob;
  tokens[eob][0].token = DCT_EOB_TOKEN;
  tokens[eob][0].qc = 0;
  *(tokens[eob] + 1) = *(tokens[eob] + 0);
  next = eob;
  for (i = eob; i-- > i0;) {
    int base_bits;
    int d2;
    int dx;

    rc = scan[i];
    x = qcoeff_ptr[rc];
    /* Only add a trellis state for non-zero coefficients. */
    if (x) {
      int shortcut = 0;
      error0 = tokens[next][0].error;
      error1 = tokens[next][1].error;
      /* Evaluate the first possibility for this state. */
      rate0 = tokens[next][0].rate;
      rate1 = tokens[next][1].rate;
      t0 = (vp8_dct_value_tokens_ptr + x)->Token;
      /* Consider both possible successor states. */
      if (next < default_eob) {
        band = bands[i + 1];
        pt = vp8_prev_token_class[t0];
        rate0 +=
          mb->token_costs[tx_type][type][band][pt][tokens[next][0].token];
        rate1 +=
          mb->token_costs[tx_type][type][band][pt][tokens[next][1].token];
      }
      UPDATE_RD_COST();
      /* And pick the best. */
      best = rd_cost1 < rd_cost0;
      base_bits = *(vp8_dct_value_cost_ptr + x);
      dx = dqcoeff_ptr[rc] - coeff_ptr[rc];
      d2 = dx * dx;
      tokens[i][0].rate = base_bits + (best ? rate1 : rate0);
      tokens[i][0].error = d2 + (best ? error1 : error0);
      tokens[i][0].next = next;
      tokens[i][0].token = t0;
      tokens[i][0].qc = x;
      best_mask[0] |= best << i;
      /* Evaluate the second possibility for this state. */
      rate0 = tokens[next][0].rate;
      rate1 = tokens[next][1].rate;

      if ((abs(x)*dequant_ptr[rc != 0] > abs(coeff_ptr[rc])) &&
          (abs(x)*dequant_ptr[rc != 0] < abs(coeff_ptr[rc]) + dequant_ptr[rc != 0]))
        shortcut = 1;
      else
        shortcut = 0;

      if (shortcut) {
        sz = -(x < 0);
        x -= 2 * sz + 1;
      }

      /* Consider both possible successor states. */
      if (!x) {
        /* If we reduced this coefficient to zero, check to see if
         *  we need to move the EOB back here.
         */
        t0 = tokens[next][0].token == DCT_EOB_TOKEN ?
             DCT_EOB_TOKEN : ZERO_TOKEN;
        t1 = tokens[next][1].token == DCT_EOB_TOKEN ?
             DCT_EOB_TOKEN : ZERO_TOKEN;
      } else {
        t0 = t1 = (vp8_dct_value_tokens_ptr + x)->Token;
      }
      if (next < default_eob) {
        band = bands[i + 1];
        if (t0 != DCT_EOB_TOKEN) {
          pt = vp8_prev_token_class[t0];
          rate0 += mb->token_costs[tx_type][type][band][pt][
              tokens[next][0].token];
        }
        if (t1 != DCT_EOB_TOKEN) {
          pt = vp8_prev_token_class[t1];
          rate1 += mb->token_costs[tx_type][type][band][pt][
              tokens[next][1].token];
        }
      }

      UPDATE_RD_COST();
      /* And pick the best. */
      best = rd_cost1 < rd_cost0;
      base_bits = *(vp8_dct_value_cost_ptr + x);

      if (shortcut) {
        dx -= (dequant_ptr[rc != 0] + sz) ^ sz;
        d2 = dx * dx;
      }
      tokens[i][1].rate = base_bits + (best ? rate1 : rate0);
      tokens[i][1].error = d2 + (best ? error1 : error0);
      tokens[i][1].next = next;
      tokens[i][1].token = best ? t1 : t0;
      tokens[i][1].qc = x;
      best_mask[1] |= best << i;
      /* Finally, make this the new head of the trellis. */
      next = i;
    }
    /* There's no choice to make for a zero coefficient, so we don't
     *  add a new trellis node, but we do need to update the costs.
     */
    else {
      band = bands[i + 1];
      t0 = tokens[next][0].token;
      t1 = tokens[next][1].token;
      /* Update the cost of each path if we're past the EOB token. */
      if (t0 != DCT_EOB_TOKEN) {
        tokens[next][0].rate += mb->token_costs[tx_type][type][band][0][t0];
        tokens[next][0].token = ZERO_TOKEN;
      }
      if (t1 != DCT_EOB_TOKEN) {
        tokens[next][1].rate += mb->token_costs[tx_type][type][band][0][t1];
        tokens[next][1].token = ZERO_TOKEN;
      }
      /* Don't update next, because we didn't add a new node. */
    }
  }

  /* Now pick the best path through the whole trellis. */
  band = bands[i + 1];
  VP8_COMBINEENTROPYCONTEXTS(pt, *a, *l);
  rate0 = tokens[next][0].rate;
  rate1 = tokens[next][1].rate;
  error0 = tokens[next][0].error;
  error1 = tokens[next][1].error;
  t0 = tokens[next][0].token;
  t1 = tokens[next][1].token;
  rate0 += mb->token_costs[tx_type][type][band][pt][t0];
  rate1 += mb->token_costs[tx_type][type][band][pt][t1];
  UPDATE_RD_COST();
  best = rd_cost1 < rd_cost0;
  final_eob = i0 - 1;
  for (i = next; i < eob; i = next) {
    x = tokens[i][best].qc;
    if (x)
      final_eob = i;
    rc = scan[i];
    qcoeff_ptr[rc] = x;
    dqcoeff_ptr[rc] = (x * dequant_ptr[rc != 0]);

    next = tokens[i][best].next;
    best = (best_mask[best] >> i) & 1;
  }
  final_eob++;

  d->eob = final_eob;
  *a = *l = (d->eob != !type);
}

/**************************************************************************
our inverse hadamard transform effectively is weighted sum of all 16 inputs
with weight either 1 or -1. It has a last stage scaling of (sum+1)>>2. And
dc only idct is (dc+16)>>5. So if all the sums are between -65 and 63 the
output after inverse wht and idct will be all zero. A sum of absolute value
smaller than 65 guarantees all 16 different (+1/-1) weighted sums in wht
fall between -65 and +65.
**************************************************************************/
#define SUM_2ND_COEFF_THRESH 65

static void check_reset_2nd_coeffs(MACROBLOCKD *xd, int type,
                                   ENTROPY_CONTEXT *a, ENTROPY_CONTEXT *l) {
  int sum = 0;
  int i;
  BLOCKD *bd = &xd->block[24];
  if (bd->dequant[0] >= SUM_2ND_COEFF_THRESH
      && bd->dequant[1] >= SUM_2ND_COEFF_THRESH)
    return;

  for (i = 0; i < bd->eob; i++) {
    int coef = bd->dqcoeff[vp8_default_zig_zag1d[i]];
    sum += (coef >= 0) ? coef : -coef;
    if (sum >= SUM_2ND_COEFF_THRESH)
      return;
  }

  if (sum < SUM_2ND_COEFF_THRESH) {
    for (i = 0; i < bd->eob; i++) {
      int rc = vp8_default_zig_zag1d[i];
      bd->qcoeff[rc] = 0;
      bd->dqcoeff[rc] = 0;
    }
    bd->eob = 0;
    *a = *l = (bd->eob != !type);
  }
}

#define SUM_2ND_COEFF_THRESH_8X8 32
static void check_reset_8x8_2nd_coeffs(MACROBLOCKD *xd, int type,
                                       ENTROPY_CONTEXT *a, ENTROPY_CONTEXT *l) {
  int sum = 0;
  BLOCKD *bd = &xd->block[24];
  int coef;

  coef = bd->dqcoeff[0];
  sum += (coef >= 0) ? coef : -coef;
  coef = bd->dqcoeff[1];
  sum += (coef >= 0) ? coef : -coef;
  coef = bd->dqcoeff[4];
  sum += (coef >= 0) ? coef : -coef;
  coef = bd->dqcoeff[8];
  sum += (coef >= 0) ? coef : -coef;

  if (sum < SUM_2ND_COEFF_THRESH_8X8) {
    bd->qcoeff[0] = 0;
    bd->dqcoeff[0] = 0;
    bd->qcoeff[1] = 0;
    bd->dqcoeff[1] = 0;
    bd->qcoeff[4] = 0;
    bd->dqcoeff[4] = 0;
    bd->qcoeff[8] = 0;
    bd->dqcoeff[8] = 0;
    bd->eob = 0;
    *a = *l = (bd->eob != !type);
  }
}

static void optimize_mb_4x4(MACROBLOCK *x, const VP8_ENCODER_RTCD *rtcd) {
  int b;
  int type;
  int has_2nd_order;
  ENTROPY_CONTEXT_PLANES t_above, t_left;
  ENTROPY_CONTEXT *ta;
  ENTROPY_CONTEXT *tl;
  MB_PREDICTION_MODE mode = x->e_mbd.mode_info_context->mbmi.mode;

  vpx_memcpy(&t_above, x->e_mbd.above_context, sizeof(ENTROPY_CONTEXT_PLANES));
  vpx_memcpy(&t_left, x->e_mbd.left_context, sizeof(ENTROPY_CONTEXT_PLANES));

  ta = (ENTROPY_CONTEXT *)&t_above;
  tl = (ENTROPY_CONTEXT *)&t_left;

  has_2nd_order = (mode != B_PRED && mode != I8X8_PRED && mode != SPLITMV);
  type = has_2nd_order ? PLANE_TYPE_Y_NO_DC : PLANE_TYPE_Y_WITH_DC;

  for (b = 0; b < 16; b++) {
    optimize_b(x, b, type,
               ta + vp8_block2above[b], tl + vp8_block2left[b], rtcd, TX_4X4);
  }

  for (b = 16; b < 24; b++) {
    optimize_b(x, b, PLANE_TYPE_UV,
               ta + vp8_block2above[b], tl + vp8_block2left[b], rtcd, TX_4X4);
  }

  if (has_2nd_order) {
    b = 24;
    optimize_b(x, b, PLANE_TYPE_Y2,
               ta + vp8_block2above[b], tl + vp8_block2left[b], rtcd, TX_4X4);
    check_reset_2nd_coeffs(&x->e_mbd, PLANE_TYPE_Y2,
                           ta + vp8_block2above[b], tl + vp8_block2left[b]);
  }
}

void vp8_optimize_mby_4x4(MACROBLOCK *x, const VP8_ENCODER_RTCD *rtcd) {
  int b;
  int type;
  int has_2nd_order;
  ENTROPY_CONTEXT_PLANES t_above, t_left;
  ENTROPY_CONTEXT *ta;
  ENTROPY_CONTEXT *tl;
  MB_PREDICTION_MODE mode = x->e_mbd.mode_info_context->mbmi.mode;

  if (!x->e_mbd.above_context || !x->e_mbd.left_context)
    return;

  vpx_memcpy(&t_above, x->e_mbd.above_context, sizeof(ENTROPY_CONTEXT_PLANES));
  vpx_memcpy(&t_left, x->e_mbd.left_context, sizeof(ENTROPY_CONTEXT_PLANES));

  ta = (ENTROPY_CONTEXT *)&t_above;
  tl = (ENTROPY_CONTEXT *)&t_left;

  has_2nd_order = (mode != B_PRED && mode != I8X8_PRED && mode != SPLITMV);
  type = has_2nd_order ? PLANE_TYPE_Y_NO_DC : PLANE_TYPE_Y_WITH_DC;

  for (b = 0; b < 16; b++) {
    optimize_b(x, b, type,
               ta + vp8_block2above[b], tl + vp8_block2left[b], rtcd, TX_4X4);
  }

  if (has_2nd_order) {
    b = 24;
    optimize_b(x, b, PLANE_TYPE_Y2,
               ta + vp8_block2above[b], tl + vp8_block2left[b], rtcd, TX_4X4);
    check_reset_2nd_coeffs(&x->e_mbd, PLANE_TYPE_Y2,
                           ta + vp8_block2above[b], tl + vp8_block2left[b]);
  }
}

void vp8_optimize_mbuv_4x4(MACROBLOCK *x, const VP8_ENCODER_RTCD *rtcd) {
  int b;
  ENTROPY_CONTEXT_PLANES t_above, t_left;
  ENTROPY_CONTEXT *ta;
  ENTROPY_CONTEXT *tl;

  if (!x->e_mbd.above_context || !x->e_mbd.left_context)
    return;

  vpx_memcpy(&t_above, x->e_mbd.above_context, sizeof(ENTROPY_CONTEXT_PLANES));
  vpx_memcpy(&t_left, x->e_mbd.left_context, sizeof(ENTROPY_CONTEXT_PLANES));

  ta = (ENTROPY_CONTEXT *)&t_above;
  tl = (ENTROPY_CONTEXT *)&t_left;

  for (b = 16; b < 24; b++) {
    optimize_b(x, b, PLANE_TYPE_UV,
               ta + vp8_block2above[b], tl + vp8_block2left[b], rtcd, TX_4X4);
  }
}

void optimize_mb_8x8(MACROBLOCK *x, const VP8_ENCODER_RTCD *rtcd) {
  int b;
  int type;
  ENTROPY_CONTEXT_PLANES t_above, t_left;
  ENTROPY_CONTEXT *ta;
  ENTROPY_CONTEXT *tl;

  vpx_memcpy(&t_above, x->e_mbd.above_context, sizeof(ENTROPY_CONTEXT_PLANES));
  vpx_memcpy(&t_left, x->e_mbd.left_context, sizeof(ENTROPY_CONTEXT_PLANES));

  ta = (ENTROPY_CONTEXT *)&t_above;
  tl = (ENTROPY_CONTEXT *)&t_left;

  type = 0;
  for (b = 0; b < 16; b += 4) {
    optimize_b(x, b, type,
               ta + vp8_block2above_8x8[b], tl + vp8_block2left_8x8[b],
               rtcd, TX_8X8);
    *(ta + vp8_block2above_8x8[b] + 1) = *(ta + vp8_block2above_8x8[b]);
    *(tl + vp8_block2left_8x8[b] + 1)  = *(tl + vp8_block2left_8x8[b]);
  }

  for (b = 16; b < 24; b += 4) {
    optimize_b(x, b, PLANE_TYPE_UV,
               ta + vp8_block2above_8x8[b], tl + vp8_block2left_8x8[b],
               rtcd, TX_8X8);
    *(ta + vp8_block2above_8x8[b] + 1) = *(ta + vp8_block2above_8x8[b]);
    *(tl + vp8_block2left_8x8[b] + 1) = *(tl + vp8_block2left_8x8[b]);
  }

  // 8x8 always have 2nd roder haar block
  check_reset_8x8_2nd_coeffs(&x->e_mbd, PLANE_TYPE_Y2,
                             ta + vp8_block2above_8x8[24], tl + vp8_block2left_8x8[24]);
}

void vp8_optimize_mby_8x8(MACROBLOCK *x, const VP8_ENCODER_RTCD *rtcd) {
  int b;
  int type;
  ENTROPY_CONTEXT_PLANES t_above, t_left;
  ENTROPY_CONTEXT *ta;
  ENTROPY_CONTEXT *tl;

  if (!x->e_mbd.above_context || !x->e_mbd.left_context)
    return;

  vpx_memcpy(&t_above, x->e_mbd.above_context, sizeof(ENTROPY_CONTEXT_PLANES));
  vpx_memcpy(&t_left, x->e_mbd.left_context, sizeof(ENTROPY_CONTEXT_PLANES));

  ta = (ENTROPY_CONTEXT *)&t_above;
  tl = (ENTROPY_CONTEXT *)&t_left;
  type = 0;
  for (b = 0; b < 16; b += 4) {
    optimize_b(x, b, type,
               ta + vp8_block2above[b], tl + vp8_block2left[b],
               rtcd, TX_8X8);
    *(ta + vp8_block2above_8x8[b] + 1) = *(ta + vp8_block2above_8x8[b]);
    *(tl + vp8_block2left_8x8[b] + 1)  = *(tl + vp8_block2left_8x8[b]);
  }

  // 8x8 always have 2nd roder haar block
  check_reset_8x8_2nd_coeffs(&x->e_mbd, PLANE_TYPE_Y2,
                             ta + vp8_block2above_8x8[24], tl + vp8_block2left_8x8[24]);
}

void vp8_optimize_mbuv_8x8(MACROBLOCK *x, const VP8_ENCODER_RTCD *rtcd) {
  int b;
  ENTROPY_CONTEXT_PLANES t_above, t_left;
  ENTROPY_CONTEXT *ta;
  ENTROPY_CONTEXT *tl;

  if (!x->e_mbd.above_context || !x->e_mbd.left_context)
    return;

  vpx_memcpy(&t_above, x->e_mbd.above_context, sizeof(ENTROPY_CONTEXT_PLANES));
  vpx_memcpy(&t_left, x->e_mbd.left_context, sizeof(ENTROPY_CONTEXT_PLANES));

  ta = (ENTROPY_CONTEXT *)&t_above;
  tl = (ENTROPY_CONTEXT *)&t_left;

  for (b = 16; b < 24; b += 4) {
    optimize_b(x, b, PLANE_TYPE_UV,
               ta + vp8_block2above_8x8[b], tl + vp8_block2left_8x8[b],
               rtcd, TX_8X8);
    *(ta + vp8_block2above_8x8[b] + 1) = *(ta + vp8_block2above_8x8[b]);
    *(tl + vp8_block2left_8x8[b] + 1) = *(tl + vp8_block2left_8x8[b]);
  }
}

void optimize_b_16x16(MACROBLOCK *mb, int i, int type,
                      ENTROPY_CONTEXT *a, ENTROPY_CONTEXT *l,
                      const VP8_ENCODER_RTCD *rtcd) {
  BLOCK *b = &mb->block[i];
  BLOCKD *d = &mb->e_mbd.block[i];
  vp8_token_state tokens[257][2];
  unsigned best_index[257][2];
  const short *dequant_ptr = d->dequant, *coeff_ptr = b->coeff;
  short *qcoeff_ptr = qcoeff_ptr = d->qcoeff;
  short *dqcoeff_ptr = dqcoeff_ptr = d->dqcoeff;
  int eob = d->eob, final_eob, sz = 0;
  int rc, x, next;
  int64_t rdmult, rddiv, rd_cost0, rd_cost1;
  int rate0, rate1, error0, error1, t0, t1;
  int best, band, pt;
  int err_mult = plane_rd_mult[type];

  /* Now set up a Viterbi trellis to evaluate alternative roundings. */
  rdmult = mb->rdmult * err_mult;
  if (mb->e_mbd.mode_info_context->mbmi.ref_frame == INTRA_FRAME)
      rdmult = (rdmult * 9)>>4;
  rddiv = mb->rddiv;
  memset(best_index, 0, sizeof(best_index));
  /* Initialize the sentinel node of the trellis. */
  tokens[eob][0].rate = 0;
  tokens[eob][0].error = 0;
  tokens[eob][0].next = 256;
  tokens[eob][0].token = DCT_EOB_TOKEN;
  tokens[eob][0].qc = 0;
  *(tokens[eob] + 1) = *(tokens[eob] + 0);
  next = eob;
  for (i = eob; i-- > 0;) {
    int base_bits, d2, dx;

    rc = vp8_default_zig_zag1d_16x16[i];
    x = qcoeff_ptr[rc];
    /* Only add a trellis state for non-zero coefficients. */
    if (x) {
      int shortcut = 0;
      error0 = tokens[next][0].error;
      error1 = tokens[next][1].error;
      /* Evaluate the first possibility for this state. */
      rate0 = tokens[next][0].rate;
      rate1 = tokens[next][1].rate;
      t0 = (vp8_dct_value_tokens_ptr + x)->Token;
      /* Consider both possible successor states. */
      if (next < 256) {
        band = vp8_coef_bands_16x16[i + 1];
        pt = vp8_prev_token_class[t0];
        rate0 += mb->token_costs[TX_16X16][type][band][pt][tokens[next][0].token];
        rate1 += mb->token_costs[TX_16X16][type][band][pt][tokens[next][1].token];
      }
      UPDATE_RD_COST();
      /* And pick the best. */
      best = rd_cost1 < rd_cost0;
      base_bits = *(vp8_dct_value_cost_ptr + x);
      dx = dqcoeff_ptr[rc] - coeff_ptr[rc];
      d2 = dx*dx;
      tokens[i][0].rate = base_bits + (best ? rate1 : rate0);
      tokens[i][0].error = d2 + (best ? error1 : error0);
      tokens[i][0].next = next;
      tokens[i][0].token = t0;
      tokens[i][0].qc = x;
      best_index[i][0] = best;
      /* Evaluate the second possibility for this state. */
      rate0 = tokens[next][0].rate;
      rate1 = tokens[next][1].rate;

      if((abs(x)*dequant_ptr[rc!=0]>abs(coeff_ptr[rc])) &&
         (abs(x)*dequant_ptr[rc!=0]<abs(coeff_ptr[rc])+dequant_ptr[rc!=0]))
        shortcut = 1;
      else
        shortcut = 0;

      if (shortcut) {
        sz = -(x < 0);
        x -= 2*sz + 1;
      }

      /* Consider both possible successor states. */
      if (!x) {
        /* If we reduced this coefficient to zero, check to see if
         *  we need to move the EOB back here.
         */
        t0 = tokens[next][0].token == DCT_EOB_TOKEN ?
             DCT_EOB_TOKEN : ZERO_TOKEN;
        t1 = tokens[next][1].token == DCT_EOB_TOKEN ?
             DCT_EOB_TOKEN : ZERO_TOKEN;
      }
      else
        t0=t1 = (vp8_dct_value_tokens_ptr + x)->Token;
      if (next < 256) {
        band = vp8_coef_bands_16x16[i + 1];
        if (t0 != DCT_EOB_TOKEN) {
            pt = vp8_prev_token_class[t0];
            rate0 += mb->token_costs[TX_16X16][type][band][pt]
                [tokens[next][0].token];
        }
        if (t1!=DCT_EOB_TOKEN) {
            pt = vp8_prev_token_class[t1];
            rate1 += mb->token_costs[TX_16X16][type][band][pt]
                [tokens[next][1].token];
        }
      }
      UPDATE_RD_COST();
      /* And pick the best. */
      best = rd_cost1 < rd_cost0;
      base_bits = *(vp8_dct_value_cost_ptr + x);

      if(shortcut) {
        dx -= (dequant_ptr[rc!=0] + sz) ^ sz;
        d2 = dx*dx;
      }
      tokens[i][1].rate = base_bits + (best ? rate1 : rate0);
      tokens[i][1].error = d2 + (best ? error1 : error0);
      tokens[i][1].next = next;
      tokens[i][1].token = best ? t1 : t0;
      tokens[i][1].qc = x;
      best_index[i][1] = best;
      /* Finally, make this the new head of the trellis. */
      next = i;
    }
    /* There's no choice to make for a zero coefficient, so we don't
     *  add a new trellis node, but we do need to update the costs.
     */
    else {
      band = vp8_coef_bands_16x16[i + 1];
      t0 = tokens[next][0].token;
      t1 = tokens[next][1].token;
      /* Update the cost of each path if we're past the EOB token. */
      if (t0 != DCT_EOB_TOKEN) {
        tokens[next][0].rate += mb->token_costs[TX_16X16][type][band][0][t0];
        tokens[next][0].token = ZERO_TOKEN;
      }
      if (t1 != DCT_EOB_TOKEN) {
        tokens[next][1].rate += mb->token_costs[TX_16X16][type][band][0][t1];
        tokens[next][1].token = ZERO_TOKEN;
      }
      /* Don't update next, because we didn't add a new node. */
    }
  }

  /* Now pick the best path through the whole trellis. */
  band = vp8_coef_bands_16x16[i + 1];
  VP8_COMBINEENTROPYCONTEXTS(pt, *a, *l);
  rate0 = tokens[next][0].rate;
  rate1 = tokens[next][1].rate;
  error0 = tokens[next][0].error;
  error1 = tokens[next][1].error;
  t0 = tokens[next][0].token;
  t1 = tokens[next][1].token;
  rate0 += mb->token_costs[TX_16X16][type][band][pt][t0];
  rate1 += mb->token_costs[TX_16X16][type][band][pt][t1];
  UPDATE_RD_COST();
  best = rd_cost1 < rd_cost0;
  final_eob = -1;

  for (i = next; i < eob; i = next) {
    x = tokens[i][best].qc;
    if (x)
      final_eob = i;
    rc = vp8_default_zig_zag1d_16x16[i];
    qcoeff_ptr[rc] = x;
    dqcoeff_ptr[rc] = (x * dequant_ptr[rc!=0]);

    next = tokens[i][best].next;
    best = best_index[i][best];
  }
  final_eob++;

  d->eob = final_eob;
  *a = *l = (d->eob != !type);
}

void vp8_optimize_mby_16x16(MACROBLOCK *x, const VP8_ENCODER_RTCD *rtcd) {
  ENTROPY_CONTEXT_PLANES t_above, t_left;
  ENTROPY_CONTEXT *ta, *tl;

  if (!x->e_mbd.above_context || !x->e_mbd.left_context)
    return;

  vpx_memcpy(&t_above, x->e_mbd.above_context, sizeof(ENTROPY_CONTEXT_PLANES));
  vpx_memcpy(&t_left, x->e_mbd.left_context, sizeof(ENTROPY_CONTEXT_PLANES));

  ta = (ENTROPY_CONTEXT *)&t_above;
  tl = (ENTROPY_CONTEXT *)&t_left;
  optimize_b_16x16(x, 0, PLANE_TYPE_Y_WITH_DC, ta, tl, rtcd);
  *(ta + 1) = *ta;
  *(tl + 1) = *tl;
}

void optimize_mb_16x16(MACROBLOCK *x, const VP8_ENCODER_RTCD *rtcd) {
  int b;
  ENTROPY_CONTEXT_PLANES t_above, t_left;
  ENTROPY_CONTEXT *ta, *tl;

  vpx_memcpy(&t_above, x->e_mbd.above_context, sizeof(ENTROPY_CONTEXT_PLANES));
  vpx_memcpy(&t_left, x->e_mbd.left_context, sizeof(ENTROPY_CONTEXT_PLANES));

  ta = (ENTROPY_CONTEXT *)&t_above;
  tl = (ENTROPY_CONTEXT *)&t_left;

  optimize_b_16x16(x, 0, PLANE_TYPE_Y_WITH_DC, ta, tl, rtcd);
  *(ta + 1) = *ta;
  *(tl + 1) = *tl;

  for (b = 16; b < 24; b += 4) {
    optimize_b(x, b, PLANE_TYPE_UV,
               ta + vp8_block2above_8x8[b], tl + vp8_block2left_8x8[b],
               rtcd, TX_8X8);
    *(ta + vp8_block2above_8x8[b] + 1) = *(ta + vp8_block2above_8x8[b]);
    *(tl + vp8_block2left_8x8[b] + 1) = *(tl + vp8_block2left_8x8[b]);
  }
}

void vp8_encode_inter16x16(const VP8_ENCODER_RTCD *rtcd, MACROBLOCK *x) {
  MACROBLOCKD *xd = &x->e_mbd;
  TX_SIZE tx_size = xd->mode_info_context->mbmi.txfm_size;

  vp8_build_inter_predictors_mb(xd);
  vp8_subtract_mb(rtcd, x);

  if (tx_size == TX_16X16) {
    vp8_transform_mb_16x16(x);
    vp8_quantize_mb_16x16(x);
    if (x->optimize)
      optimize_mb_16x16(x, rtcd);
    vp8_inverse_transform_mb_16x16(IF_RTCD(&rtcd->common->idct), xd);
  } else if (tx_size == TX_8X8) {
    vp8_transform_mb_8x8(x);
    vp8_quantize_mb_8x8(x);
    if (x->optimize)
      optimize_mb_8x8(x, rtcd);
    vp8_inverse_transform_mb_8x8(IF_RTCD(&rtcd->common->idct), xd);
  } else {
    transform_mb_4x4(x);
    vp8_quantize_mb_4x4(x);
    if (x->optimize)
      optimize_mb_4x4(x, rtcd);
    vp8_inverse_transform_mb_4x4(IF_RTCD(&rtcd->common->idct), xd);
  }

  RECON_INVOKE(&rtcd->common->recon, recon_mb)(IF_RTCD(&rtcd->common->recon),
                                               xd);
}

/* this function is used by first pass only */
void vp8_encode_inter16x16y(const VP8_ENCODER_RTCD *rtcd, MACROBLOCK *x) {
  MACROBLOCKD *xd = &x->e_mbd;
  BLOCK *b = &x->block[0];

#if CONFIG_PRED_FILTER
  // Disable the prediction filter for firstpass
  xd->mode_info_context->mbmi.pred_filter_enabled = 0;
#endif

  vp8_build_1st_inter16x16_predictors_mby(xd, xd->predictor, 16, 0);

  ENCODEMB_INVOKE(&rtcd->encodemb, submby)(x->src_diff, *(b->base_src),
                                           xd->predictor, b->src_stride);

  vp8_transform_mby_4x4(x);
  vp8_quantize_mby_4x4(x);
  vp8_inverse_transform_mby_4x4(IF_RTCD(&rtcd->common->idct), xd);

  RECON_INVOKE(&rtcd->common->recon, recon_mby)(IF_RTCD(&rtcd->common->recon),
                                                xd);
}
