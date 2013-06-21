/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include "./vpx_config.h"
#include "vp9/encoder/vp9_encodemb.h"
#include "vp9/common/vp9_reconinter.h"
#include "vp9/encoder/vp9_quantize.h"
#include "vp9/encoder/vp9_tokenize.h"
#include "vp9/common/vp9_reconintra.h"
#include "vpx_mem/vpx_mem.h"
#include "vp9/encoder/vp9_rdopt.h"
#include "vp9/common/vp9_systemdependent.h"
#include "vp9_rtcd.h"

DECLARE_ALIGNED(16, extern const uint8_t,
                vp9_pt_energy_class[MAX_ENTROPY_TOKENS]);

void vp9_subtract_block_c(int rows, int cols,
                          int16_t *diff_ptr, ptrdiff_t diff_stride,
                          const uint8_t *src_ptr, ptrdiff_t src_stride,
                          const uint8_t *pred_ptr, ptrdiff_t pred_stride) {
  int r, c;

  for (r = 0; r < rows; r++) {
    for (c = 0; c < cols; c++)
      diff_ptr[c] = src_ptr[c] - pred_ptr[c];

    diff_ptr += diff_stride;
    pred_ptr += pred_stride;
    src_ptr  += src_stride;
  }
}

static void inverse_transform_b_4x4_add(MACROBLOCKD *xd, int eob,
                                        int16_t *dqcoeff, uint8_t *dest,
                                        int stride) {
  if (eob <= 1)
    xd->inv_txm4x4_1_add(dqcoeff, dest, stride);
  else
    xd->inv_txm4x4_add(dqcoeff, dest, stride);
}


static void subtract_plane(MACROBLOCK *x, BLOCK_SIZE_TYPE bsize, int plane) {
  struct macroblock_plane *const p = &x->plane[plane];
  const MACROBLOCKD *const xd = &x->e_mbd;
  const struct macroblockd_plane *const pd = &xd->plane[plane];
  const int bw = plane_block_width(bsize, pd);
  const int bh = plane_block_height(bsize, pd);

  vp9_subtract_block(bh, bw, p->src_diff, bw,
                     p->src.buf, p->src.stride,
                     pd->dst.buf, pd->dst.stride);
}

void vp9_subtract_sby(MACROBLOCK *x, BLOCK_SIZE_TYPE bsize) {
  subtract_plane(x, bsize, 0);
}

void vp9_subtract_sbuv(MACROBLOCK *x, BLOCK_SIZE_TYPE bsize) {
  int i;

  for (i = 1; i < MAX_MB_PLANE; i++)
    subtract_plane(x, bsize, i);
}

void vp9_subtract_sb(MACROBLOCK *x, BLOCK_SIZE_TYPE bsize) {
  vp9_subtract_sby(x, bsize);
  vp9_subtract_sbuv(x, bsize);
}


#define RDTRUNC(RM,DM,R,D) ( (128+(R)*(RM)) & 0xFF )
#define RDTRUNC_8x8(RM,DM,R,D) ( (128+(R)*(RM)) & 0xFF )
typedef struct vp9_token_state vp9_token_state;

struct vp9_token_state {
  int           rate;
  int           error;
  int           next;
  signed char   token;
  short         qc;
};

// TODO: experiments to find optimal multiple numbers
#define Y1_RD_MULT 4
#define UV_RD_MULT 2

static const int plane_rd_mult[4] = {
  Y1_RD_MULT,
  UV_RD_MULT,
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

// This function is a place holder for now but may ultimately need
// to scan previous tokens to work out the correct context.
static int trellis_get_coeff_context(const int *scan,
                                     const int *nb,
                                     int idx, int token,
                                     uint8_t *token_cache,
                                     int pad, int l) {
  int bak = token_cache[scan[idx]], pt;
  token_cache[scan[idx]] = vp9_pt_energy_class[token];
  pt = vp9_get_coef_context(scan, nb, pad, token_cache, idx + 1, l);
  token_cache[scan[idx]] = bak;
  return pt;
}

static void optimize_b(VP9_COMMON *const cm, MACROBLOCK *mb,
                       int plane, int block, BLOCK_SIZE_TYPE bsize,
                       ENTROPY_CONTEXT *a, ENTROPY_CONTEXT *l,
                       TX_SIZE tx_size) {
  const int ref = mb->e_mbd.mode_info_context->mbmi.ref_frame[0] != INTRA_FRAME;
  MACROBLOCKD *const xd = &mb->e_mbd;
  vp9_token_state tokens[1025][2];
  unsigned best_index[1025][2];
  const int16_t *coeff_ptr = BLOCK_OFFSET(mb->plane[plane].coeff,
                                          block, 16);
  int16_t *qcoeff_ptr;
  int16_t *dqcoeff_ptr;
  int eob = xd->plane[plane].eobs[block], final_eob, sz = 0;
  const int i0 = 0;
  int rc, x, next, i;
  int64_t rdmult, rddiv, rd_cost0, rd_cost1;
  int rate0, rate1, error0, error1, t0, t1;
  int best, band, pt;
  PLANE_TYPE type = xd->plane[plane].plane_type;
  int err_mult = plane_rd_mult[type];
  int default_eob, pad;
  int const *scan, *nb;
  const int mul = 1 + (tx_size == TX_32X32);
  uint8_t token_cache[1024];
  const int ib = txfrm_block_to_raster_block(xd, bsize, plane,
                                             block, 2 * tx_size);
  const int16_t *dequant_ptr = xd->plane[plane].dequant;
  const uint8_t * band_translate;

  assert((!type && !plane) || (type && plane));
  dqcoeff_ptr = BLOCK_OFFSET(xd->plane[plane].dqcoeff, block, 16);
  qcoeff_ptr = BLOCK_OFFSET(xd->plane[plane].qcoeff, block, 16);
  switch (tx_size) {
    default:
    case TX_4X4: {
      const TX_TYPE tx_type = plane == 0 ? get_tx_type_4x4(xd, ib) : DCT_DCT;
      default_eob = 16;
      scan = get_scan_4x4(tx_type);
      band_translate = vp9_coefband_trans_4x4;
      break;
    }
    case TX_8X8: {
      const TX_TYPE tx_type = plane == 0 ? get_tx_type_8x8(xd, ib) : DCT_DCT;
      scan = get_scan_8x8(tx_type);
      default_eob = 64;
      band_translate = vp9_coefband_trans_8x8plus;
      break;
    }
    case TX_16X16: {
      const TX_TYPE tx_type = plane == 0 ? get_tx_type_16x16(xd, ib) : DCT_DCT;
      scan = get_scan_16x16(tx_type);
      default_eob = 256;
      band_translate = vp9_coefband_trans_8x8plus;
      break;
    }
    case TX_32X32:
      scan = vp9_default_scan_32x32;
      default_eob = 1024;
      band_translate = vp9_coefband_trans_8x8plus;
      break;
  }
  assert(eob <= default_eob);

  /* Now set up a Viterbi trellis to evaluate alternative roundings. */
  rdmult = mb->rdmult * err_mult;
  if (mb->e_mbd.mode_info_context->mbmi.ref_frame[0] == INTRA_FRAME)
    rdmult = (rdmult * 9) >> 4;
  rddiv = mb->rddiv;
  memset(best_index, 0, sizeof(best_index));
  /* Initialize the sentinel node of the trellis. */
  tokens[eob][0].rate = 0;
  tokens[eob][0].error = 0;
  tokens[eob][0].next = default_eob;
  tokens[eob][0].token = DCT_EOB_TOKEN;
  tokens[eob][0].qc = 0;
  *(tokens[eob] + 1) = *(tokens[eob] + 0);
  next = eob;
  for (i = 0; i < eob; i++)
    token_cache[scan[i]] = vp9_pt_energy_class[vp9_dct_value_tokens_ptr[
        qcoeff_ptr[scan[i]]].token];
  nb = vp9_get_coef_neighbors_handle(scan, &pad);

  for (i = eob; i-- > i0;) {
    int base_bits, d2, dx;

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
      t0 = (vp9_dct_value_tokens_ptr + x)->token;
      /* Consider both possible successor states. */
      if (next < default_eob) {
        band = get_coef_band(band_translate, i + 1);
        pt = trellis_get_coeff_context(scan, nb, i, t0, token_cache,
                                       pad, default_eob);
        rate0 +=
          mb->token_costs_noskip[tx_size][type][ref][band][pt]
                                [tokens[next][0].token];
        rate1 +=
          mb->token_costs_noskip[tx_size][type][ref][band][pt]
                                [tokens[next][1].token];
      }
      UPDATE_RD_COST();
      /* And pick the best. */
      best = rd_cost1 < rd_cost0;
      base_bits = *(vp9_dct_value_cost_ptr + x);
      dx = mul * (dqcoeff_ptr[rc] - coeff_ptr[rc]);
      d2 = dx * dx;
      tokens[i][0].rate = base_bits + (best ? rate1 : rate0);
      tokens[i][0].error = d2 + (best ? error1 : error0);
      tokens[i][0].next = next;
      tokens[i][0].token = t0;
      tokens[i][0].qc = x;
      best_index[i][0] = best;

      /* Evaluate the second possibility for this state. */
      rate0 = tokens[next][0].rate;
      rate1 = tokens[next][1].rate;

      if ((abs(x)*dequant_ptr[rc != 0] > abs(coeff_ptr[rc]) * mul) &&
          (abs(x)*dequant_ptr[rc != 0] < abs(coeff_ptr[rc]) * mul +
                                         dequant_ptr[rc != 0]))
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
        t0 = t1 = (vp9_dct_value_tokens_ptr + x)->token;
      }
      if (next < default_eob) {
        band = get_coef_band(band_translate, i + 1);
        if (t0 != DCT_EOB_TOKEN) {
          pt = trellis_get_coeff_context(scan, nb, i, t0, token_cache,
                                         pad, default_eob);
          if (!x)
            rate0 += mb->token_costs[tx_size][type][ref][band][pt][
                tokens[next][0].token];
          else
            rate0 += mb->token_costs_noskip[tx_size][type][ref][band][pt][
                tokens[next][0].token];
        }
        if (t1 != DCT_EOB_TOKEN) {
          pt = trellis_get_coeff_context(scan, nb, i, t1, token_cache,
                                         pad, default_eob);
          if (!x)
            rate1 += mb->token_costs[tx_size][type][ref][band][pt][
                tokens[next][1].token];
          else
            rate1 += mb->token_costs_noskip[tx_size][type][ref][band][pt][
                tokens[next][1].token];
        }
      }

      UPDATE_RD_COST();
      /* And pick the best. */
      best = rd_cost1 < rd_cost0;
      base_bits = *(vp9_dct_value_cost_ptr + x);

      if (shortcut) {
        dx -= (dequant_ptr[rc != 0] + sz) ^ sz;
        d2 = dx * dx;
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
      band = get_coef_band(band_translate, i + 1);
      t0 = tokens[next][0].token;
      t1 = tokens[next][1].token;
      /* Update the cost of each path if we're past the EOB token. */
      if (t0 != DCT_EOB_TOKEN) {
        tokens[next][0].rate +=
            mb->token_costs[tx_size][type][ref][band][0][t0];
        tokens[next][0].token = ZERO_TOKEN;
      }
      if (t1 != DCT_EOB_TOKEN) {
        tokens[next][1].rate +=
            mb->token_costs[tx_size][type][ref][band][0][t1];
        tokens[next][1].token = ZERO_TOKEN;
      }
      /* Don't update next, because we didn't add a new node. */
    }
  }

  /* Now pick the best path through the whole trellis. */
  band = get_coef_band(band_translate, i + 1);
  pt = combine_entropy_contexts(*a, *l);
  rate0 = tokens[next][0].rate;
  rate1 = tokens[next][1].rate;
  error0 = tokens[next][0].error;
  error1 = tokens[next][1].error;
  t0 = tokens[next][0].token;
  t1 = tokens[next][1].token;
  rate0 += mb->token_costs_noskip[tx_size][type][ref][band][pt][t0];
  rate1 += mb->token_costs_noskip[tx_size][type][ref][band][pt][t1];
  UPDATE_RD_COST();
  best = rd_cost1 < rd_cost0;
  final_eob = i0 - 1;
  vpx_memset(qcoeff_ptr, 0, sizeof(*qcoeff_ptr) * (16 << (tx_size * 2)));
  vpx_memset(dqcoeff_ptr, 0, sizeof(*dqcoeff_ptr) * (16 << (tx_size * 2)));
  for (i = next; i < eob; i = next) {
    x = tokens[i][best].qc;
    if (x) {
      final_eob = i;
    }
    rc = scan[i];
    qcoeff_ptr[rc] = x;
    dqcoeff_ptr[rc] = (x * dequant_ptr[rc != 0]) / mul;

    next = tokens[i][best].next;
    best = best_index[i][best];
  }
  final_eob++;

  xd->plane[plane].eobs[block] = final_eob;
  *a = *l = (final_eob > 0);
}

struct optimize_block_args {
  VP9_COMMON *cm;
  MACROBLOCK *x;
  struct optimize_ctx *ctx;
};

void vp9_optimize_b(int plane, int block, BLOCK_SIZE_TYPE bsize,
                    int ss_txfrm_size, VP9_COMMON *cm, MACROBLOCK *mb,
                    struct optimize_ctx *ctx) {
  MACROBLOCKD *const xd = &mb->e_mbd;
  int x, y;

  // find current entropy context
  txfrm_block_to_raster_xy(xd, bsize, plane, block, ss_txfrm_size, &x, &y);

  optimize_b(cm, mb, plane, block, bsize,
             &ctx->ta[plane][x], &ctx->tl[plane][y], ss_txfrm_size / 2);
}

static void optimize_block(int plane, int block, BLOCK_SIZE_TYPE bsize,
                           int ss_txfrm_size, void *arg) {
  const struct optimize_block_args* const args = arg;
  vp9_optimize_b(plane, block, bsize, ss_txfrm_size, args->cm, args->x,
                 args->ctx);
}

void vp9_optimize_init(MACROBLOCKD *xd, BLOCK_SIZE_TYPE bsize,
                       struct optimize_ctx *ctx) {
  int p;

  for (p = 0; p < MAX_MB_PLANE; p++) {
    const struct macroblockd_plane* const plane = &xd->plane[p];
    const int bwl = b_width_log2(bsize) - plane->subsampling_x;
    const int bhl = b_height_log2(bsize) - plane->subsampling_y;
    const MB_MODE_INFO *mbmi = &xd->mode_info_context->mbmi;
    const TX_SIZE tx_size = p ? get_uv_tx_size(mbmi)
                              : mbmi->txfm_size;
    int i, j;

    for (i = 0; i < 1 << bwl; i += 1 << tx_size) {
      int c = 0;
      ctx->ta[p][i] = 0;
      for (j = 0; j < 1 << tx_size && !c; j++) {
        c = ctx->ta[p][i] |= plane->above_context[i + j];
      }
    }
    for (i = 0; i < 1 << bhl; i += 1 << tx_size) {
      int c = 0;
      ctx->tl[p][i] = 0;
      for (j = 0; j < 1 << tx_size && !c; j++) {
        c = ctx->tl[p][i] |= plane->left_context[i + j];
      }
    }
  }
}

void vp9_optimize_sby(VP9_COMMON *cm, MACROBLOCK *x, BLOCK_SIZE_TYPE bsize) {
  struct optimize_ctx ctx;
  struct optimize_block_args arg = {cm, x, &ctx};
  vp9_optimize_init(&x->e_mbd, bsize, &ctx);
  foreach_transformed_block_in_plane(&x->e_mbd, bsize, 0, optimize_block, &arg);
}

void vp9_optimize_sbuv(VP9_COMMON *const cm, MACROBLOCK *x,
                       BLOCK_SIZE_TYPE bsize) {
  struct optimize_ctx ctx;
  struct optimize_block_args arg = {cm, x, &ctx};
  vp9_optimize_init(&x->e_mbd, bsize, &ctx);
  foreach_transformed_block_uv(&x->e_mbd, bsize, optimize_block, &arg);
}

struct encode_b_args {
  VP9_COMMON *cm;
  MACROBLOCK *x;
  struct optimize_ctx *ctx;
};

static void xform_quant(int plane, int block, BLOCK_SIZE_TYPE bsize,
                         int ss_txfrm_size, void *arg) {
  struct encode_b_args* const args = arg;
  MACROBLOCK* const x = args->x;
  MACROBLOCKD* const xd = &x->e_mbd;
  const int bw = plane_block_width(bsize, &xd->plane[plane]);
  const int raster_block = txfrm_block_to_raster_block(xd, bsize, plane,
                                                       block, ss_txfrm_size);
  int16_t *const coeff = BLOCK_OFFSET(x->plane[plane].coeff, block, 16);
  int16_t *const src_diff = raster_block_offset_int16(xd, bsize, plane,
                                                      raster_block,
                                                      x->plane[plane].src_diff);
  TX_TYPE tx_type = DCT_DCT;

  switch (ss_txfrm_size / 2) {
    case TX_32X32:
      if (x->rd_search)
        vp9_short_fdct32x32_rd(src_diff, coeff, bw * 2);
      else
        vp9_short_fdct32x32(src_diff, coeff, bw * 2);
      break;
    case TX_16X16:
      tx_type = plane == 0 ? get_tx_type_16x16(xd, raster_block) : DCT_DCT;
      if (tx_type != DCT_DCT)
        vp9_short_fht16x16(src_diff, coeff, bw, tx_type);
      else
        x->fwd_txm16x16(src_diff, coeff, bw * 2);
      break;
    case TX_8X8:
      tx_type = plane == 0 ? get_tx_type_8x8(xd, raster_block) : DCT_DCT;
      if (tx_type != DCT_DCT)
        vp9_short_fht8x8(src_diff, coeff, bw, tx_type);
      else
        x->fwd_txm8x8(src_diff, coeff, bw * 2);
      break;
    case TX_4X4:
      tx_type = plane == 0 ? get_tx_type_4x4(xd, raster_block) : DCT_DCT;
      if (tx_type != DCT_DCT)
        vp9_short_fht4x4(src_diff, coeff, bw, tx_type);
      else
        x->fwd_txm4x4(src_diff, coeff, bw * 2);
      break;
    default:
      assert(0);
  }

  vp9_quantize(x, plane, block, 16 << ss_txfrm_size, tx_type);
}

static void encode_block(int plane, int block, BLOCK_SIZE_TYPE bsize,
                         int ss_txfrm_size, void *arg) {
  struct encode_b_args *const args = arg;
  MACROBLOCK *const x = args->x;
  MACROBLOCKD *const xd = &x->e_mbd;
  const int raster_block = txfrm_block_to_raster_block(xd, bsize, plane,
                                                       block, ss_txfrm_size);
  struct macroblockd_plane *const pd = &xd->plane[plane];
  int16_t *const dqcoeff = BLOCK_OFFSET(pd->dqcoeff, block, 16);
  uint8_t *const dst = raster_block_offset_uint8(xd, bsize, plane,
                                                 raster_block,
                                                 pd->dst.buf, pd->dst.stride);
  TX_TYPE tx_type = DCT_DCT;

  xform_quant(plane, block, bsize, ss_txfrm_size, arg);

  if (x->optimize)
    vp9_optimize_b(plane, block, bsize, ss_txfrm_size, args->cm, x, args->ctx);

  switch (ss_txfrm_size / 2) {
    case TX_32X32:
      vp9_short_idct32x32_add(dqcoeff, dst, pd->dst.stride);
      break;
    case TX_16X16:
      tx_type = plane == 0 ? get_tx_type_16x16(xd, raster_block) : DCT_DCT;
      if (tx_type == DCT_DCT)
        vp9_short_idct16x16_add(dqcoeff, dst, pd->dst.stride);
      else
        vp9_short_iht16x16_add(dqcoeff, dst, pd->dst.stride, tx_type);
      break;
    case TX_8X8:
      tx_type = plane == 0 ? get_tx_type_8x8(xd, raster_block) : DCT_DCT;
      if (tx_type == DCT_DCT)
        vp9_short_idct8x8_add(dqcoeff, dst, pd->dst.stride);
      else
        vp9_short_iht8x8_add(dqcoeff, dst, pd->dst.stride, tx_type);
      break;
    case TX_4X4:
      tx_type = plane == 0 ? get_tx_type_4x4(xd, raster_block) : DCT_DCT;
      if (tx_type == DCT_DCT)
        // this is like vp9_short_idct4x4 but has a special case around eob<=1
        // which is significant (not just an optimization) for the lossless
        // case.
        inverse_transform_b_4x4_add(xd, pd->eobs[block], dqcoeff,
                                    dst, pd->dst.stride);
      else
        vp9_short_iht4x4_add(dqcoeff, dst, pd->dst.stride, tx_type);
      break;
  }
}

void vp9_xform_quant_sby(VP9_COMMON *cm, MACROBLOCK *x, BLOCK_SIZE_TYPE bsize) {
  MACROBLOCKD* const xd = &x->e_mbd;
  struct encode_b_args arg = {cm, x, NULL};

  foreach_transformed_block_in_plane(xd, bsize, 0, xform_quant, &arg);
}

void vp9_xform_quant_sbuv(VP9_COMMON *cm, MACROBLOCK *x,
                          BLOCK_SIZE_TYPE bsize) {
  MACROBLOCKD* const xd = &x->e_mbd;
  struct encode_b_args arg = {cm, x, NULL};

  foreach_transformed_block_uv(xd, bsize, xform_quant, &arg);
}

void vp9_encode_sby(VP9_COMMON *cm, MACROBLOCK *x, BLOCK_SIZE_TYPE bsize) {
  MACROBLOCKD *const xd = &x->e_mbd;
  struct optimize_ctx ctx;
  struct encode_b_args arg = {cm, x, &ctx};

  vp9_subtract_sby(x, bsize);
  if (x->optimize)
    vp9_optimize_init(xd, bsize, &ctx);

  foreach_transformed_block_in_plane(xd, bsize, 0, encode_block, &arg);
}

void vp9_encode_sbuv(VP9_COMMON *cm, MACROBLOCK *x, BLOCK_SIZE_TYPE bsize) {
  MACROBLOCKD *const xd = &x->e_mbd;
  struct optimize_ctx ctx;
  struct encode_b_args arg = {cm, x, &ctx};

  vp9_subtract_sbuv(x, bsize);
  if (x->optimize)
    vp9_optimize_init(xd, bsize, &ctx);

  foreach_transformed_block_uv(xd, bsize, encode_block, &arg);
}

void vp9_encode_sb(VP9_COMMON *cm, MACROBLOCK *x, BLOCK_SIZE_TYPE bsize) {
  MACROBLOCKD *const xd = &x->e_mbd;
  struct optimize_ctx ctx;
  struct encode_b_args arg = {cm, x, &ctx};

  vp9_subtract_sb(x, bsize);
  if (x->optimize)
    vp9_optimize_init(xd, bsize, &ctx);

  foreach_transformed_block(xd, bsize, encode_block, &arg);
}

static void encode_block_intra(int plane, int block, BLOCK_SIZE_TYPE bsize,
                               int ss_txfrm_size, void *arg) {
  struct encode_b_args* const args = arg;
  MACROBLOCK *const x = args->x;
  MACROBLOCKD *const xd = &x->e_mbd;
  MB_MODE_INFO *const mbmi = &xd->mode_info_context->mbmi;
  const TX_SIZE tx_size = (TX_SIZE)(ss_txfrm_size / 2);
  struct macroblock_plane *const p = &x->plane[plane];
  struct macroblockd_plane *const pd = &xd->plane[plane];
  int16_t *const dqcoeff = BLOCK_OFFSET(pd->dqcoeff, block, 16);
  const int bw = plane_block_width(bsize, pd);
  const int raster_block = txfrm_block_to_raster_block(xd, bsize, plane,
                                                       block, ss_txfrm_size);

  uint8_t *const src = raster_block_offset_uint8(xd, bsize, plane, raster_block,
                                                 p->src.buf, p->src.stride);
  uint8_t *const dst = raster_block_offset_uint8(xd, bsize, plane, raster_block,
                                                 pd->dst.buf, pd->dst.stride);
  int16_t *const src_diff = raster_block_offset_int16(xd, bsize, plane,
                                                      raster_block,
                                                      p->src_diff);

  const int txfm_b_size = 4 << tx_size;
  int ib = raster_block;
  int tx_ib = ib >> tx_size;
  int plane_b_size;

  TX_TYPE tx_type;
  int mode, b_mode;

  if (xd->mb_to_right_edge < 0 || xd->mb_to_bottom_edge < 0) {
    extend_for_intra(xd, plane, block, bsize, ss_txfrm_size);
  }

  mode = plane == 0? mbmi->mode: mbmi->uv_mode;
  if (plane == 0 &&
      mbmi->sb_type < BLOCK_SIZE_SB8X8 &&
      mbmi->ref_frame[0] == INTRA_FRAME)
    b_mode = xd->mode_info_context->bmi[ib].as_mode.first;
  else
    b_mode = mode;

  assert(b_mode >= DC_PRED && b_mode <= TM_PRED);

  plane_b_size = b_width_log2(bsize) - pd->subsampling_x;
  vp9_predict_intra_block(xd, tx_ib, plane_b_size, tx_size, b_mode,
                          dst, pd->dst.stride);
  vp9_subtract_block(txfm_b_size, txfm_b_size, src_diff, bw,
                     src, p->src.stride, dst, pd->dst.stride);

  xform_quant(plane, block, bsize, ss_txfrm_size, arg);


  // if (x->optimize)
  // vp9_optimize_b(plane, block, bsize, ss_txfrm_size,
  //                args->cm, x, args->ctx);

  switch (ss_txfrm_size / 2) {
    case TX_32X32:
        vp9_short_idct32x32_add(dqcoeff, dst, pd->dst.stride);
      break;
    case TX_16X16:
      tx_type = plane == 0 ? get_tx_type_16x16(xd, raster_block) : DCT_DCT;
      if (tx_type == DCT_DCT)
        vp9_short_idct16x16_add(dqcoeff, dst, pd->dst.stride);
      else
        vp9_short_iht16x16_add(dqcoeff, dst, pd->dst.stride, tx_type);
      break;
    case TX_8X8:
      tx_type = plane == 0 ? get_tx_type_8x8(xd, raster_block) : DCT_DCT;
      if (tx_type == DCT_DCT)
        vp9_short_idct8x8_add(dqcoeff, dst, pd->dst.stride);
      else
        vp9_short_iht8x8_add(dqcoeff, dst, pd->dst.stride, tx_type);
      break;
    case TX_4X4:
      tx_type = plane == 0 ? get_tx_type_4x4(xd, raster_block) : DCT_DCT;
      if (tx_type == DCT_DCT)
        // this is like vp9_short_idct4x4 but has a special case around eob<=1
        // which is significant (not just an optimization) for the lossless
        // case.
        inverse_transform_b_4x4_add(xd, pd->eobs[block], dqcoeff,
                                    dst, pd->dst.stride);
      else
        vp9_short_iht4x4_add(dqcoeff, dst, pd->dst.stride, tx_type);
      break;
  }
}

void vp9_encode_intra_block_y(VP9_COMMON *cm, MACROBLOCK *x,
                              BLOCK_SIZE_TYPE bsize) {
  MACROBLOCKD* const xd = &x->e_mbd;
  struct optimize_ctx ctx;
  struct encode_b_args arg = {cm, x, &ctx};

  foreach_transformed_block_in_plane(xd, bsize, 0,
                                     encode_block_intra, &arg);
}
void vp9_encode_intra_block_uv(VP9_COMMON *cm, MACROBLOCK *x,
                              BLOCK_SIZE_TYPE bsize) {
  MACROBLOCKD* const xd = &x->e_mbd;
  struct optimize_ctx ctx;
  struct encode_b_args arg = {cm, x, &ctx};
  foreach_transformed_block_uv(xd, bsize, encode_block_intra, &arg);
}

