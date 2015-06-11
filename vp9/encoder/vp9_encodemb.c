/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */


#include "./vp9_rtcd.h"
#include "./vpx_config.h"

#include "vpx_mem/vpx_mem.h"

#include "vp9/common/vp9_idct.h"
#include "vp9/common/vp9_reconinter.h"
#include "vp9/common/vp9_reconintra.h"
#include "vp9/common/vp9_systemdependent.h"

#include "vp9/encoder/vp9_encodemb.h"
#include "vp9/encoder/vp9_quantize.h"
#include "vp9/encoder/vp9_rd.h"
#include "vp9/encoder/vp9_tokenize.h"

struct optimize_ctx {
  ENTROPY_CONTEXT ta[MAX_MB_PLANE][16];
  ENTROPY_CONTEXT tl[MAX_MB_PLANE][16];
};

void vp9_subtract_block_c(int rows, int cols,
                          int16_t *diff, ptrdiff_t diff_stride,
                          const uint8_t *src, ptrdiff_t src_stride,
                          const uint8_t *pred, ptrdiff_t pred_stride) {
  int r, c;

  for (r = 0; r < rows; r++) {
    for (c = 0; c < cols; c++)
      diff[c] = src[c] - pred[c];

    diff += diff_stride;
    pred += pred_stride;
    src  += src_stride;
  }
}

#if CONFIG_VP9_HIGHBITDEPTH
void vp9_highbd_subtract_block_c(int rows, int cols,
                                 int16_t *diff, ptrdiff_t diff_stride,
                                 const uint8_t *src8, ptrdiff_t src_stride,
                                 const uint8_t *pred8, ptrdiff_t pred_stride,
                                 int bd) {
  int r, c;
  uint16_t *src = CONVERT_TO_SHORTPTR(src8);
  uint16_t *pred = CONVERT_TO_SHORTPTR(pred8);
  (void) bd;

  for (r = 0; r < rows; r++) {
    for (c = 0; c < cols; c++) {
      diff[c] = src[c] - pred[c];
    }

    diff += diff_stride;
    pred += pred_stride;
    src  += src_stride;
  }
}
#endif  // CONFIG_VP9_HIGHBITDEPTH

void vp9_subtract_plane(MACROBLOCK *x, BLOCK_SIZE bsize, int plane) {
  struct macroblock_plane *const p = &x->plane[plane];
  const struct macroblockd_plane *const pd = &x->e_mbd.plane[plane];
  const BLOCK_SIZE plane_bsize = get_plane_block_size(bsize, pd);
  const int bw = 4 * num_4x4_blocks_wide_lookup[plane_bsize];
  const int bh = 4 * num_4x4_blocks_high_lookup[plane_bsize];

#if CONFIG_VP9_HIGHBITDEPTH
  if (x->e_mbd.cur_buf->flags & YV12_FLAG_HIGHBITDEPTH) {
    vp9_highbd_subtract_block(bh, bw, p->src_diff, bw, p->src.buf,
                              p->src.stride, pd->dst.buf, pd->dst.stride,
                              x->e_mbd.bd);
    return;
  }
#endif  // CONFIG_VP9_HIGHBITDEPTH
  vp9_subtract_block(bh, bw, p->src_diff, bw, p->src.buf, p->src.stride,
                     pd->dst.buf, pd->dst.stride);
}

#define RDTRUNC(RM, DM, R, D) ((128 + (R) * (RM)) & 0xFF)

typedef struct vp9_token_state {
  int           rate;
  int           error;
  int           next;
  int16_t       token;
  short         qc;
} vp9_token_state;

// TODO(jimbankoski): experiment to find optimal RD numbers.
static const int plane_rd_mult[PLANE_TYPES] = { 4, 2 };

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
static int trellis_get_coeff_context(const int16_t *scan,
                                     const int16_t *nb,
                                     int idx, int token,
                                     uint8_t *token_cache) {
  int bak = token_cache[scan[idx]], pt;
  token_cache[scan[idx]] = vp9_pt_energy_class[token];
  pt = get_coef_context(nb, token_cache, idx + 1);
  token_cache[scan[idx]] = bak;
  return pt;
}

static int optimize_b(MACROBLOCK *mb, int plane, int block,
                      TX_SIZE tx_size, int ctx) {
  MACROBLOCKD *const xd = &mb->e_mbd;
  struct macroblock_plane *const p = &mb->plane[plane];
  struct macroblockd_plane *const pd = &xd->plane[plane];
  const int ref = is_inter_block(&xd->mi[0].src_mi->mbmi);
  vp9_token_state tokens[1025][2];
  unsigned best_index[1025][2];
  uint8_t token_cache[1024];
  const tran_low_t *const coeff = BLOCK_OFFSET(mb->plane[plane].coeff, block);
  tran_low_t *const qcoeff = BLOCK_OFFSET(p->qcoeff, block);
  tran_low_t *const dqcoeff = BLOCK_OFFSET(pd->dqcoeff, block);
  const int eob = p->eobs[block];
  const PLANE_TYPE type = pd->plane_type;
  const int default_eob = 16 << (tx_size << 1);
  const int mul = 1 + (tx_size == TX_32X32);
  const int16_t *dequant_ptr = pd->dequant;
  const uint8_t *const band_translate = get_band_translate(tx_size);
  const scan_order *const so = get_scan(xd, tx_size, type, block);
  const int16_t *const scan = so->scan;
  const int16_t *const nb = so->neighbors;
  int next = eob, sz = 0;
  int64_t rdmult = mb->rdmult * plane_rd_mult[type], rddiv = mb->rddiv;
  int64_t rd_cost0, rd_cost1;
  int rate0, rate1, error0, error1;
  int16_t t0, t1;
  EXTRABIT e0;
  int best, band, pt, i, final_eob;
#if CONFIG_VP9_HIGHBITDEPTH
  const int16_t *cat6_high_cost = vp9_get_high_cost_table(xd->bd);
#else
  const int16_t *cat6_high_cost = vp9_get_high_cost_table(8);
#endif

  assert((!type && !plane) || (type && plane));
  assert(eob <= default_eob);

  /* Now set up a Viterbi trellis to evaluate alternative roundings. */
  if (!ref)
    rdmult = (rdmult * 9) >> 4;

  /* Initialize the sentinel node of the trellis. */
  tokens[eob][0].rate = 0;
  tokens[eob][0].error = 0;
  tokens[eob][0].next = default_eob;
  tokens[eob][0].token = EOB_TOKEN;
  tokens[eob][0].qc = 0;
  tokens[eob][1] = tokens[eob][0];

  for (i = 0; i < eob; i++)
    token_cache[scan[i]] =
        vp9_pt_energy_class[vp9_get_token(qcoeff[scan[i]])];

  for (i = eob; i-- > 0;) {
    int base_bits, d2, dx;
    const int rc = scan[i];
    int x = qcoeff[rc];
    /* Only add a trellis state for non-zero coefficients. */
    if (x) {
      int shortcut = 0;
      error0 = tokens[next][0].error;
      error1 = tokens[next][1].error;
      /* Evaluate the first possibility for this state. */
      rate0 = tokens[next][0].rate;
      rate1 = tokens[next][1].rate;
      vp9_get_token_extra(x, &t0, &e0);
      /* Consider both possible successor states. */
      if (next < default_eob) {
        band = band_translate[i + 1];
        pt = trellis_get_coeff_context(scan, nb, i, t0, token_cache);
        rate0 += mb->token_costs[tx_size][type][ref][band][0][pt]
                                [tokens[next][0].token];
        rate1 += mb->token_costs[tx_size][type][ref][band][0][pt]
                                [tokens[next][1].token];
      }
      UPDATE_RD_COST();
      /* And pick the best. */
      best = rd_cost1 < rd_cost0;
      base_bits = vp9_get_cost(t0, e0, cat6_high_cost);
      dx = mul * (dqcoeff[rc] - coeff[rc]);
#if CONFIG_VP9_HIGHBITDEPTH
      if (xd->cur_buf->flags & YV12_FLAG_HIGHBITDEPTH) {
        dx >>= xd->bd - 8;
      }
#endif  // CONFIG_VP9_HIGHBITDEPTH
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

      if ((abs(x) * dequant_ptr[rc != 0] > abs(coeff[rc]) * mul) &&
          (abs(x) * dequant_ptr[rc != 0] < abs(coeff[rc]) * mul +
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
        t0 = tokens[next][0].token == EOB_TOKEN ? EOB_TOKEN : ZERO_TOKEN;
        t1 = tokens[next][1].token == EOB_TOKEN ? EOB_TOKEN : ZERO_TOKEN;
        e0 = 0;
      } else {
        vp9_get_token_extra(x, &t0, &e0);
        t1 = t0;
      }
      if (next < default_eob) {
        band = band_translate[i + 1];
        if (t0 != EOB_TOKEN) {
          pt = trellis_get_coeff_context(scan, nb, i, t0, token_cache);
          rate0 += mb->token_costs[tx_size][type][ref][band][!x][pt]
                                  [tokens[next][0].token];
        }
        if (t1 != EOB_TOKEN) {
          pt = trellis_get_coeff_context(scan, nb, i, t1, token_cache);
          rate1 += mb->token_costs[tx_size][type][ref][band][!x][pt]
                                  [tokens[next][1].token];
        }
      }

      UPDATE_RD_COST();
      /* And pick the best. */
      best = rd_cost1 < rd_cost0;
      base_bits = vp9_get_cost(t0, e0, cat6_high_cost);

      if (shortcut) {
#if CONFIG_VP9_HIGHBITDEPTH
        if (xd->cur_buf->flags & YV12_FLAG_HIGHBITDEPTH) {
          dx -= ((dequant_ptr[rc != 0] >> (xd->bd - 8)) + sz) ^ sz;
        } else {
          dx -= (dequant_ptr[rc != 0] + sz) ^ sz;
        }
#else
        dx -= (dequant_ptr[rc != 0] + sz) ^ sz;
#endif  // CONFIG_VP9_HIGHBITDEPTH
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
    } else {
      /* There's no choice to make for a zero coefficient, so we don't
       *  add a new trellis node, but we do need to update the costs.
       */
      band = band_translate[i + 1];
      t0 = tokens[next][0].token;
      t1 = tokens[next][1].token;
      /* Update the cost of each path if we're past the EOB token. */
      if (t0 != EOB_TOKEN) {
        tokens[next][0].rate +=
            mb->token_costs[tx_size][type][ref][band][1][0][t0];
        tokens[next][0].token = ZERO_TOKEN;
      }
      if (t1 != EOB_TOKEN) {
        tokens[next][1].rate +=
            mb->token_costs[tx_size][type][ref][band][1][0][t1];
        tokens[next][1].token = ZERO_TOKEN;
      }
      best_index[i][0] = best_index[i][1] = 0;
      /* Don't update next, because we didn't add a new node. */
    }
  }

  /* Now pick the best path through the whole trellis. */
  band = band_translate[i + 1];
  rate0 = tokens[next][0].rate;
  rate1 = tokens[next][1].rate;
  error0 = tokens[next][0].error;
  error1 = tokens[next][1].error;
  t0 = tokens[next][0].token;
  t1 = tokens[next][1].token;
  rate0 += mb->token_costs[tx_size][type][ref][band][0][ctx][t0];
  rate1 += mb->token_costs[tx_size][type][ref][band][0][ctx][t1];
  UPDATE_RD_COST();
  best = rd_cost1 < rd_cost0;
  final_eob = -1;
  vpx_memset(qcoeff, 0, sizeof(*qcoeff) * (16 << (tx_size * 2)));
  vpx_memset(dqcoeff, 0, sizeof(*dqcoeff) * (16 << (tx_size * 2)));
  for (i = next; i < eob; i = next) {
    const int x = tokens[i][best].qc;
    const int rc = scan[i];
    if (x) {
      final_eob = i;
    }

    qcoeff[rc] = x;
    dqcoeff[rc] = (x * dequant_ptr[rc != 0]) / mul;

    next = tokens[i][best].next;
    best = best_index[i][best];
  }
  final_eob++;

  mb->plane[plane].eobs[block] = final_eob;
  return final_eob;
}

static INLINE void fdct32x32(int rd_transform,
                             const int16_t *src, tran_low_t *dst,
                             int src_stride) {
  if (rd_transform)
    vp9_fdct32x32_rd(src, dst, src_stride);
  else
    vp9_fdct32x32(src, dst, src_stride);
}

#if CONFIG_VP9_HIGHBITDEPTH
static INLINE void highbd_fdct32x32(int rd_transform, const int16_t *src,
                                    tran_low_t *dst, int src_stride) {
  if (rd_transform)
    vp9_highbd_fdct32x32_rd(src, dst, src_stride);
  else
    vp9_highbd_fdct32x32(src, dst, src_stride);
}
#endif  // CONFIG_VP9_HIGHBITDEPTH

void vp9_xform_quant_fp(MACROBLOCK *x, int plane, int block,
                        BLOCK_SIZE plane_bsize, TX_SIZE tx_size) {
  MACROBLOCKD *const xd = &x->e_mbd;
  const struct macroblock_plane *const p = &x->plane[plane];
  const struct macroblockd_plane *const pd = &xd->plane[plane];
  const scan_order *const scan_order = &vp9_default_scan_orders[tx_size];
  tran_low_t *const coeff = BLOCK_OFFSET(p->coeff, block);
  tran_low_t *const qcoeff = BLOCK_OFFSET(p->qcoeff, block);
  tran_low_t *const dqcoeff = BLOCK_OFFSET(pd->dqcoeff, block);
  uint16_t *const eob = &p->eobs[block];
  const int diff_stride = 4 * num_4x4_blocks_wide_lookup[plane_bsize];
  int i, j;
  const int16_t *src_diff;
  txfrm_block_to_raster_xy(plane_bsize, tx_size, block, &i, &j);
  src_diff = &p->src_diff[4 * (j * diff_stride + i)];

#if CONFIG_VP9_HIGHBITDEPTH
  if (xd->cur_buf->flags & YV12_FLAG_HIGHBITDEPTH) {
    switch (tx_size) {
      case TX_32X32:
        highbd_fdct32x32(x->use_lp32x32fdct, src_diff, coeff, diff_stride);
        vp9_highbd_quantize_fp_32x32(coeff, 1024, x->skip_block, p->zbin,
                                     p->round_fp, p->quant_fp, p->quant_shift,
                                     qcoeff, dqcoeff, pd->dequant,
                                     eob, scan_order->scan,
                                     scan_order->iscan);
        break;
      case TX_16X16:
        vp9_highbd_fdct16x16(src_diff, coeff, diff_stride);
        vp9_highbd_quantize_fp(coeff, 256, x->skip_block, p->zbin, p->round_fp,
                               p->quant_fp, p->quant_shift, qcoeff, dqcoeff,
                               pd->dequant, eob,
                               scan_order->scan, scan_order->iscan);
        break;
      case TX_8X8:
        vp9_highbd_fdct8x8(src_diff, coeff, diff_stride);
        vp9_highbd_quantize_fp(coeff, 64, x->skip_block, p->zbin, p->round_fp,
                               p->quant_fp, p->quant_shift, qcoeff, dqcoeff,
                               pd->dequant, eob,
                               scan_order->scan, scan_order->iscan);
        break;
      case TX_4X4:
        x->fwd_txm4x4(src_diff, coeff, diff_stride);
        vp9_highbd_quantize_fp(coeff, 16, x->skip_block, p->zbin, p->round_fp,
                               p->quant_fp, p->quant_shift, qcoeff, dqcoeff,
                               pd->dequant, eob,
                               scan_order->scan, scan_order->iscan);
        break;
      default:
        assert(0);
    }
    return;
  }
#endif  // CONFIG_VP9_HIGHBITDEPTH

  switch (tx_size) {
    case TX_32X32:
      fdct32x32(x->use_lp32x32fdct, src_diff, coeff, diff_stride);
      vp9_quantize_fp_32x32(coeff, 1024, x->skip_block, p->zbin, p->round_fp,
                            p->quant_fp, p->quant_shift, qcoeff, dqcoeff,
                            pd->dequant, eob, scan_order->scan,
                            scan_order->iscan);
      break;
    case TX_16X16:
      vp9_fdct16x16(src_diff, coeff, diff_stride);
      vp9_quantize_fp(coeff, 256, x->skip_block, p->zbin, p->round_fp,
                      p->quant_fp, p->quant_shift, qcoeff, dqcoeff,
                      pd->dequant, eob,
                      scan_order->scan, scan_order->iscan);
      break;
    case TX_8X8:
      vp9_fdct8x8_quant(src_diff, diff_stride, coeff, 64,
                        x->skip_block, p->zbin, p->round_fp,
                        p->quant_fp, p->quant_shift, qcoeff, dqcoeff,
                        pd->dequant, eob,
                        scan_order->scan, scan_order->iscan);
      break;
    case TX_4X4:
      x->fwd_txm4x4(src_diff, coeff, diff_stride);
      vp9_quantize_fp(coeff, 16, x->skip_block, p->zbin, p->round_fp,
                      p->quant_fp, p->quant_shift, qcoeff, dqcoeff,
                      pd->dequant, eob,
                      scan_order->scan, scan_order->iscan);
      break;
    default:
      assert(0);
      break;
  }
}

void vp9_xform_quant_dc(MACROBLOCK *x, int plane, int block,
                        BLOCK_SIZE plane_bsize, TX_SIZE tx_size) {
  MACROBLOCKD *const xd = &x->e_mbd;
  const struct macroblock_plane *const p = &x->plane[plane];
  const struct macroblockd_plane *const pd = &xd->plane[plane];
  tran_low_t *const coeff = BLOCK_OFFSET(p->coeff, block);
  tran_low_t *const qcoeff = BLOCK_OFFSET(p->qcoeff, block);
  tran_low_t *const dqcoeff = BLOCK_OFFSET(pd->dqcoeff, block);
  uint16_t *const eob = &p->eobs[block];
  const int diff_stride = 4 * num_4x4_blocks_wide_lookup[plane_bsize];
  int i, j;
  const int16_t *src_diff;

  txfrm_block_to_raster_xy(plane_bsize, tx_size, block, &i, &j);
  src_diff = &p->src_diff[4 * (j * diff_stride + i)];

#if CONFIG_VP9_HIGHBITDEPTH
  if (xd->cur_buf->flags & YV12_FLAG_HIGHBITDEPTH) {
    switch (tx_size) {
      case TX_32X32:
        vp9_highbd_fdct32x32_1(src_diff, coeff, diff_stride);
        vp9_highbd_quantize_dc_32x32(coeff, x->skip_block, p->round,
                                     p->quant_fp[0], qcoeff, dqcoeff,
                                     pd->dequant[0], eob);
        break;
      case TX_16X16:
        vp9_highbd_fdct16x16_1(src_diff, coeff, diff_stride);
        vp9_highbd_quantize_dc(coeff, 256, x->skip_block, p->round,
                               p->quant_fp[0], qcoeff, dqcoeff,
                               pd->dequant[0], eob);
        break;
      case TX_8X8:
        vp9_highbd_fdct8x8_1(src_diff, coeff, diff_stride);
        vp9_highbd_quantize_dc(coeff, 64, x->skip_block, p->round,
                               p->quant_fp[0], qcoeff, dqcoeff,
                               pd->dequant[0], eob);
        break;
      case TX_4X4:
        x->fwd_txm4x4(src_diff, coeff, diff_stride);
        vp9_highbd_quantize_dc(coeff, 16, x->skip_block, p->round,
                               p->quant_fp[0], qcoeff, dqcoeff,
                               pd->dequant[0], eob);
        break;
      default:
        assert(0);
    }
    return;
  }
#endif  // CONFIG_VP9_HIGHBITDEPTH

  switch (tx_size) {
    case TX_32X32:
      vp9_fdct32x32_1(src_diff, coeff, diff_stride);
      vp9_quantize_dc_32x32(coeff, x->skip_block, p->round,
                            p->quant_fp[0], qcoeff, dqcoeff,
                            pd->dequant[0], eob);
      break;
    case TX_16X16:
      vp9_fdct16x16_1(src_diff, coeff, diff_stride);
      vp9_quantize_dc(coeff, 256, x->skip_block, p->round,
                     p->quant_fp[0], qcoeff, dqcoeff,
                     pd->dequant[0], eob);
      break;
    case TX_8X8:
      vp9_fdct8x8_1(src_diff, coeff, diff_stride);
      vp9_quantize_dc(coeff, 64, x->skip_block, p->round,
                      p->quant_fp[0], qcoeff, dqcoeff,
                      pd->dequant[0], eob);
      break;
    case TX_4X4:
      x->fwd_txm4x4(src_diff, coeff, diff_stride);
      vp9_quantize_dc(coeff, 16, x->skip_block, p->round,
                      p->quant_fp[0], qcoeff, dqcoeff,
                      pd->dequant[0], eob);
      break;
    default:
      assert(0);
      break;
  }
}

void vp9_xform_quant_inter(MACROBLOCK *x, int plane, int block,
                           int blk_row, int blk_col,
                           BLOCK_SIZE plane_bsize, TX_SIZE tx_size) {
  MACROBLOCKD *const xd = &x->e_mbd;
  const struct macroblock_plane *const p = &x->plane[plane];
  const struct macroblockd_plane *const pd = &xd->plane[plane];
  const scan_order *const scan_order = &vp9_default_scan_orders[tx_size];
  tran_low_t *const coeff = BLOCK_OFFSET(p->coeff, block);
  tran_low_t *const qcoeff = BLOCK_OFFSET(p->qcoeff, block);
  tran_low_t *const dqcoeff = BLOCK_OFFSET(pd->dqcoeff, block);
  uint16_t *const eob = &p->eobs[block];
  const int diff_stride = 4 * num_4x4_blocks_wide_lookup[plane_bsize];
  const int16_t *src_diff;
  src_diff = &p->src_diff[4 * (blk_row * diff_stride + blk_col)];

#if CONFIG_VP9_HIGHBITDEPTH
  if (xd->cur_buf->flags & YV12_FLAG_HIGHBITDEPTH) {
     switch (tx_size) {
      case TX_32X32:
        highbd_fdct32x32(x->use_lp32x32fdct, src_diff, coeff, diff_stride);
        vp9_highbd_quantize_b_32x32(coeff, 1024, x->skip_block, p->zbin,
                                    p->round, p->quant, p->quant_shift, qcoeff,
                                    dqcoeff, pd->dequant, eob,
                                    scan_order->scan, scan_order->iscan);
        break;
      case TX_16X16:
        vp9_highbd_fdct16x16(src_diff, coeff, diff_stride);
        vp9_highbd_quantize_b(coeff, 256, x->skip_block, p->zbin, p->round,
                              p->quant, p->quant_shift, qcoeff, dqcoeff,
                              pd->dequant, eob,
                              scan_order->scan, scan_order->iscan);
        break;
      case TX_8X8:
        vp9_highbd_fdct8x8(src_diff, coeff, diff_stride);
        vp9_highbd_quantize_b(coeff, 64, x->skip_block, p->zbin, p->round,
                              p->quant, p->quant_shift, qcoeff, dqcoeff,
                              pd->dequant, eob,
                              scan_order->scan, scan_order->iscan);
        break;
      case TX_4X4:
        x->fwd_txm4x4(src_diff, coeff, diff_stride);
        vp9_highbd_quantize_b(coeff, 16, x->skip_block, p->zbin, p->round,
                              p->quant, p->quant_shift, qcoeff, dqcoeff,
                              pd->dequant, eob,
                              scan_order->scan, scan_order->iscan);
        break;
      default:
        assert(0);
    }
    return;
  }
#endif  // CONFIG_VP9_HIGHBITDEPTH

  switch (tx_size) {
    case TX_32X32:
      fdct32x32(x->use_lp32x32fdct, src_diff, coeff, diff_stride);
      vp9_quantize_b_32x32(coeff, 1024, x->skip_block, p->zbin, p->round,
                           p->quant, p->quant_shift, qcoeff, dqcoeff,
                           pd->dequant, eob, scan_order->scan,
                           scan_order->iscan);
      break;
    case TX_16X16:
      vp9_fdct16x16(src_diff, coeff, diff_stride);
      vp9_quantize_b(coeff, 256, x->skip_block, p->zbin, p->round,
                     p->quant, p->quant_shift, qcoeff, dqcoeff,
                     pd->dequant, eob,
                     scan_order->scan, scan_order->iscan);
      break;
    case TX_8X8:
      vp9_fdct8x8(src_diff, coeff, diff_stride);
      vp9_quantize_b(coeff, 64, x->skip_block, p->zbin, p->round,
                     p->quant, p->quant_shift, qcoeff, dqcoeff,
                     pd->dequant, eob,
                     scan_order->scan, scan_order->iscan);
      break;
    case TX_4X4:
      x->fwd_txm4x4(src_diff, coeff, diff_stride);
      vp9_quantize_b(coeff, 16, x->skip_block, p->zbin, p->round,
                     p->quant, p->quant_shift, qcoeff, dqcoeff,
                     pd->dequant, eob,
                     scan_order->scan, scan_order->iscan);
      break;
    default:
      assert(0);
      break;
  }
}

void vp9_xform_quant(MACROBLOCK *x, int plane, int block,
                     BLOCK_SIZE plane_bsize, TX_SIZE tx_size) {
  MACROBLOCKD *const xd = &x->e_mbd;
  const struct macroblock_plane *const p = &x->plane[plane];
  const struct macroblockd_plane *const pd = &xd->plane[plane];
  const scan_order *const scan_order = &vp9_default_scan_orders[tx_size];
  tran_low_t *const coeff = BLOCK_OFFSET(p->coeff, block);
  tran_low_t *const qcoeff = BLOCK_OFFSET(p->qcoeff, block);
  tran_low_t *const dqcoeff = BLOCK_OFFSET(pd->dqcoeff, block);
  uint16_t *const eob = &p->eobs[block];
  const int diff_stride = 4 * num_4x4_blocks_wide_lookup[plane_bsize];
  int i, j;
  const int16_t *src_diff;
  txfrm_block_to_raster_xy(plane_bsize, tx_size, block, &i, &j);
  src_diff = &p->src_diff[4 * (j * diff_stride + i)];

#if CONFIG_VP9_HIGHBITDEPTH
  if (xd->cur_buf->flags & YV12_FLAG_HIGHBITDEPTH) {
     switch (tx_size) {
      case TX_32X32:
        highbd_fdct32x32(x->use_lp32x32fdct, src_diff, coeff, diff_stride);
        vp9_highbd_quantize_b_32x32(coeff, 1024, x->skip_block, p->zbin,
                                    p->round, p->quant, p->quant_shift, qcoeff,
                                    dqcoeff, pd->dequant, eob,
                                    scan_order->scan, scan_order->iscan);
        break;
      case TX_16X16:
        vp9_highbd_fdct16x16(src_diff, coeff, diff_stride);
        vp9_highbd_quantize_b(coeff, 256, x->skip_block, p->zbin, p->round,
                              p->quant, p->quant_shift, qcoeff, dqcoeff,
                              pd->dequant, eob,
                              scan_order->scan, scan_order->iscan);
        break;
      case TX_8X8:
        vp9_highbd_fdct8x8(src_diff, coeff, diff_stride);
        vp9_highbd_quantize_b(coeff, 64, x->skip_block, p->zbin, p->round,
                              p->quant, p->quant_shift, qcoeff, dqcoeff,
                              pd->dequant, eob,
                              scan_order->scan, scan_order->iscan);
        break;
      case TX_4X4:
        x->fwd_txm4x4(src_diff, coeff, diff_stride);
        vp9_highbd_quantize_b(coeff, 16, x->skip_block, p->zbin, p->round,
                              p->quant, p->quant_shift, qcoeff, dqcoeff,
                              pd->dequant, eob,
                              scan_order->scan, scan_order->iscan);
        break;
      default:
        assert(0);
    }
    return;
  }
#endif  // CONFIG_VP9_HIGHBITDEPTH

  switch (tx_size) {
    case TX_32X32:
      fdct32x32(x->use_lp32x32fdct, src_diff, coeff, diff_stride);
      vp9_quantize_b_32x32(coeff, 1024, x->skip_block, p->zbin, p->round,
                           p->quant, p->quant_shift, qcoeff, dqcoeff,
                           pd->dequant, eob, scan_order->scan,
                           scan_order->iscan);
      break;
    case TX_16X16:
      vp9_fdct16x16(src_diff, coeff, diff_stride);
      vp9_quantize_b(coeff, 256, x->skip_block, p->zbin, p->round,
                     p->quant, p->quant_shift, qcoeff, dqcoeff,
                     pd->dequant, eob,
                     scan_order->scan, scan_order->iscan);
      break;
    case TX_8X8:
      vp9_fdct8x8(src_diff, coeff, diff_stride);
      vp9_quantize_b(coeff, 64, x->skip_block, p->zbin, p->round,
                     p->quant, p->quant_shift, qcoeff, dqcoeff,
                     pd->dequant, eob,
                     scan_order->scan, scan_order->iscan);
      break;
    case TX_4X4:
      x->fwd_txm4x4(src_diff, coeff, diff_stride);
      vp9_quantize_b(coeff, 16, x->skip_block, p->zbin, p->round,
                     p->quant, p->quant_shift, qcoeff, dqcoeff,
                     pd->dequant, eob,
                     scan_order->scan, scan_order->iscan);
      break;
    default:
      assert(0);
      break;
  }
}

static void encode_block_b(int blk_row, int blk_col, int plane,
                           int block, BLOCK_SIZE plane_bsize,
                           TX_SIZE tx_size, void *arg) {
  struct encode_b_args *const args = arg;
  MACROBLOCK *const x = args->x;
  MACROBLOCKD *const xd = &x->e_mbd;
  struct optimize_ctx *const ctx = args->ctx;
  struct macroblock_plane *const p = &x->plane[plane];
  struct macroblockd_plane *const pd = &xd->plane[plane];
  tran_low_t *const dqcoeff = BLOCK_OFFSET(pd->dqcoeff, block);
  uint8_t *dst;
  ENTROPY_CONTEXT *a, *l;
  const int block_stride = num_4x4_blocks_wide_lookup[plane_bsize];
  int i;
  dst = &pd->dst.buf[4 * blk_row * pd->dst.stride + 4 * blk_col];
  a = &ctx->ta[plane][blk_col];
  l = &ctx->tl[plane][blk_row];

  // TODO(jingning): per transformed block zero forcing only enabled for
  // luma component. will integrate chroma components as well.
//  if (x->zcoeff_blk[tx_size][block] && plane == 0) {
//    p->eobs[block] = 0;
//    *a = *l = 0;
//    return;
//  }

  if (x->blk_skip[plane][blk_row * block_stride + blk_col] == 0)
    vp9_xform_quant_inter(x, plane, block, blk_row, blk_col,
                          plane_bsize, tx_size);
  else
    p->eobs[block] = 0;

  if (x->optimize) {
    int context;
    switch (tx_size) {
      case TX_4X4:
        break;
      case TX_8X8:
        a[0] = !!*(const uint16_t *)&a[0];
        l[0] = !!*(const uint16_t *)&l[0];
        break;
      case TX_16X16:
        a[0] = !!*(const uint32_t *)&a[0];
        l[0] = !!*(const uint32_t *)&l[0];
        break;
      case TX_32X32:
        a[0] = !!*(const uint64_t *)&a[0];
        l[0] = !!*(const uint64_t *)&l[0];
        break;
      default:
        assert(0 && "Invalid transform size.");
        break;
    }
    context = combine_entropy_contexts(*a, *l);
    *a = *l = optimize_b(x, plane, block, tx_size, context) > 0;
  } else {
    *a = *l = p->eobs[block] > 0;
  }

  for (i = 0; i < (1 << tx_size); ++i) {
    a[i] = a[0];
    l[i] = l[0];
  }

  if (p->eobs[block])
    *(args->skip) = 0;

  if (x->skip_encode || p->eobs[block] == 0)
    return;
#if CONFIG_VP9_HIGHBITDEPTH
  if (xd->cur_buf->flags & YV12_FLAG_HIGHBITDEPTH) {
    switch (tx_size) {
      case TX_32X32:
        vp9_highbd_idct32x32_add(dqcoeff, dst, pd->dst.stride,
                                 p->eobs[block], xd->bd);
        break;
      case TX_16X16:
        vp9_highbd_idct16x16_add(dqcoeff, dst, pd->dst.stride,
                                 p->eobs[block], xd->bd);
        break;
      case TX_8X8:
        vp9_highbd_idct8x8_add(dqcoeff, dst, pd->dst.stride,
                               p->eobs[block], xd->bd);
        break;
      case TX_4X4:
        // this is like vp9_short_idct4x4 but has a special case around eob<=1
        // which is significant (not just an optimization) for the lossless
        // case.
        x->highbd_itxm_add(dqcoeff, dst, pd->dst.stride,
                           p->eobs[block], xd->bd);
        break;
      default:
        assert(0 && "Invalid transform size");
    }
    return;
  }
#endif  // CONFIG_VP9_HIGHBITDEPTH

  switch (tx_size) {
    case TX_32X32:
      vp9_idct32x32_add(dqcoeff, dst, pd->dst.stride, p->eobs[block]);
      break;
    case TX_16X16:
      vp9_idct16x16_add(dqcoeff, dst, pd->dst.stride, p->eobs[block]);
      break;
    case TX_8X8:
      vp9_idct8x8_add(dqcoeff, dst, pd->dst.stride, p->eobs[block]);
      break;
    case TX_4X4:
      // this is like vp9_short_idct4x4 but has a special case around eob<=1
      // which is significant (not just an optimization) for the lossless
      // case.
      x->itxm_add(dqcoeff, dst, pd->dst.stride, p->eobs[block]);
      break;
    default:
      assert(0 && "Invalid transform size");
      break;
  }
}

static void encode_block_inter(int blk_row, int blk_col,
                               int plane, int block, BLOCK_SIZE plane_bsize,
                               TX_SIZE tx_size, void *arg) {
  struct encode_b_args *const args = arg;
  MACROBLOCK *const x = args->x;
  MACROBLOCKD *const xd = &x->e_mbd;
  MB_MODE_INFO *const mbmi = &xd->mi[0].src_mi->mbmi;
  struct macroblockd_plane *const pd = &xd->plane[plane];
  int tx_idx = (blk_row >> (1 - pd->subsampling_y)) * 8 +
               (blk_col >> (1 - pd->subsampling_x));
  TX_SIZE plane_tx_size = plane ?
      get_uv_tx_size_impl(mbmi->inter_tx_size[tx_idx], plane_bsize, 0, 0) :
      mbmi->inter_tx_size[tx_idx];
  int max_blocks_high = num_4x4_blocks_high_lookup[plane_bsize];
  int max_blocks_wide = num_4x4_blocks_wide_lookup[plane_bsize];
  if (xd->mb_to_bottom_edge < 0)
    max_blocks_high += xd->mb_to_bottom_edge >> (5 + pd->subsampling_y);
  if (xd->mb_to_right_edge < 0)
    max_blocks_wide += xd->mb_to_right_edge >> (5 + pd->subsampling_x);

  if (blk_row >= max_blocks_high || blk_col >= max_blocks_wide)
    return;

  if (tx_size == plane_tx_size) {
    encode_block_b(blk_row, blk_col, plane, block, plane_bsize, tx_size, arg);
  } else {
    BLOCK_SIZE bsize = txsize_to_bsize[tx_size];
    int bh = num_4x4_blocks_high_lookup[bsize];
    int step = 1 << (2 *(tx_size - 1));
    int i;
    for (i = 0; i < 4; ++i) {
      int offsetr = (i >> 1) * bh / 2;
      int offsetc = (i & 0x01) * bh / 2;
      encode_block_inter(blk_row + offsetr, blk_col + offsetc,
                         plane, block + i * step, plane_bsize,
                         tx_size - 1, arg);
    }
  }
}

static void encode_block_pass1(int plane, int block, BLOCK_SIZE plane_bsize,
                               TX_SIZE tx_size, void *arg) {
  MACROBLOCK *const x = (MACROBLOCK *)arg;
  MACROBLOCKD *const xd = &x->e_mbd;
  struct macroblock_plane *const p = &x->plane[plane];
  struct macroblockd_plane *const pd = &xd->plane[plane];
  tran_low_t *const dqcoeff = BLOCK_OFFSET(pd->dqcoeff, block);
  int i, j;
  uint8_t *dst;
  txfrm_block_to_raster_xy(plane_bsize, tx_size, block, &i, &j);
  dst = &pd->dst.buf[4 * j * pd->dst.stride + 4 * i];

  vp9_xform_quant(x, plane, block, plane_bsize, tx_size);

  if (p->eobs[block] > 0) {
#if CONFIG_VP9_HIGHBITDEPTH
    if (xd->cur_buf->flags & YV12_FLAG_HIGHBITDEPTH) {
       x->highbd_itxm_add(dqcoeff, dst, pd->dst.stride, p->eobs[block], xd->bd);
       return;
    }
#endif  // CONFIG_VP9_HIGHBITDEPTH
    x->itxm_add(dqcoeff, dst, pd->dst.stride, p->eobs[block]);
  }
}

void vp9_encode_sby_pass1(MACROBLOCK *x, BLOCK_SIZE bsize) {
  vp9_subtract_plane(x, bsize, 0);
  vp9_foreach_transformed_block_in_plane(&x->e_mbd, bsize, 0,
                                         encode_block_pass1, x);
}

void vp9_encode_sb(MACROBLOCK *x, BLOCK_SIZE bsize) {
  MACROBLOCKD *const xd = &x->e_mbd;
  struct optimize_ctx ctx;
  MB_MODE_INFO *mbmi = &xd->mi[0].src_mi->mbmi;
  struct encode_b_args arg = {x, &ctx, &mbmi->skip};
  int plane;

  mbmi->skip = 1;

  if (x->skip)
    return;

  for (plane = 0; plane < MAX_MB_PLANE; ++plane) {
    const struct macroblockd_plane *const pd = &xd->plane[plane];
    const BLOCK_SIZE plane_bsize = get_plane_block_size(bsize, pd);
    const int mi_width = num_4x4_blocks_wide_lookup[plane_bsize];
    const int mi_height = num_4x4_blocks_high_lookup[plane_bsize];
    BLOCK_SIZE txb_size = txsize_to_bsize[max_txsize_lookup[plane_bsize]];
    int bh = num_4x4_blocks_wide_lookup[txb_size];
    int idx, idy;
    int block = 0;
    int step = 1 << (max_txsize_lookup[plane_bsize] * 2);

    vp9_subtract_plane(x, bsize, plane);

    vp9_get_entropy_contexts(bsize, TX_4X4, pd, ctx.ta[plane], ctx.tl[plane]);

    for (idy = 0; idy < mi_height; idy += bh) {
      for (idx = 0; idx < mi_width; idx += bh) {
        encode_block_inter(idy, idx, plane, block, plane_bsize,
                           max_txsize_lookup[plane_bsize], &arg);
        block += step;
      }
    }
  }
}

void vp9_encode_block_intra(int plane, int block, BLOCK_SIZE plane_bsize,
                               TX_SIZE tx_size, void *arg) {
  struct encode_b_args* const args = arg;
  MACROBLOCK *const x = args->x;
  MACROBLOCKD *const xd = &x->e_mbd;
  MB_MODE_INFO *mbmi = &xd->mi[0].src_mi->mbmi;
  struct macroblock_plane *const p = &x->plane[plane];
  struct macroblockd_plane *const pd = &xd->plane[plane];
  tran_low_t *coeff = BLOCK_OFFSET(p->coeff, block);
  tran_low_t *qcoeff = BLOCK_OFFSET(p->qcoeff, block);
  tran_low_t *dqcoeff = BLOCK_OFFSET(pd->dqcoeff, block);
  const scan_order *scan_order;
  TX_TYPE tx_type;
  PREDICTION_MODE mode;
  const int bwl = b_width_log2_lookup[plane_bsize];
  const int diff_stride = 4 * (1 << bwl);
  uint8_t *src, *dst;
  int16_t *src_diff;
  uint16_t *eob = &p->eobs[block];
  const int src_stride = p->src.stride;
  const int dst_stride = pd->dst.stride;
  int i, j;

  txfrm_block_to_raster_xy(plane_bsize, tx_size, block, &i, &j);
  dst = &pd->dst.buf[4 * (j * dst_stride + i)];
  src = &p->src.buf[4 * (j * src_stride + i)];
  src_diff = &p->src_diff[4 * (j * diff_stride + i)];

#if CONFIG_VP9_HIGHBITDEPTH
  if (xd->cur_buf->flags & YV12_FLAG_HIGHBITDEPTH) {
    switch (tx_size) {
      case TX_32X32:
        scan_order = &vp9_default_scan_orders[TX_32X32];
        mode = plane == 0 ? mbmi->mode : mbmi->uv_mode;
        vp9_predict_intra_block(xd, block >> 6, bwl, TX_32X32, mode,
                                x->skip_encode ? src : dst,
                                x->skip_encode ? src_stride : dst_stride,
                                dst, dst_stride, i, j, plane);
        if (!x->skip_recode) {
          vp9_highbd_subtract_block(32, 32, src_diff, diff_stride,
                                    src, src_stride, dst, dst_stride, xd->bd);
          highbd_fdct32x32(x->use_lp32x32fdct, src_diff, coeff, diff_stride);
          vp9_highbd_quantize_b_32x32(coeff, 1024, x->skip_block, p->zbin,
                                      p->round, p->quant, p->quant_shift,
                                      qcoeff, dqcoeff, pd->dequant, eob,
                                      scan_order->scan, scan_order->iscan);
        }
        if (!x->skip_encode && *eob) {
          vp9_highbd_idct32x32_add(dqcoeff, dst, dst_stride, *eob, xd->bd);
        }
        break;
      case TX_16X16:
        tx_type = get_tx_type(pd->plane_type, xd);
        scan_order = &vp9_scan_orders[TX_16X16][tx_type];
        mode = plane == 0 ? mbmi->mode : mbmi->uv_mode;
        vp9_predict_intra_block(xd, block >> 4, bwl, TX_16X16, mode,
                                x->skip_encode ? src : dst,
                                x->skip_encode ? src_stride : dst_stride,
                                dst, dst_stride, i, j, plane);
        if (!x->skip_recode) {
          vp9_highbd_subtract_block(16, 16, src_diff, diff_stride,
                                    src, src_stride, dst, dst_stride, xd->bd);
          vp9_highbd_fht16x16(src_diff, coeff, diff_stride, tx_type);
          vp9_highbd_quantize_b(coeff, 256, x->skip_block, p->zbin, p->round,
                                p->quant, p->quant_shift, qcoeff, dqcoeff,
                                pd->dequant, eob,
                                scan_order->scan, scan_order->iscan);
        }
        if (!x->skip_encode && *eob) {
          vp9_highbd_iht16x16_add(tx_type, dqcoeff, dst, dst_stride,
                                  *eob, xd->bd);
        }
        break;
      case TX_8X8:
        tx_type = get_tx_type(pd->plane_type, xd);
        scan_order = &vp9_scan_orders[TX_8X8][tx_type];
        mode = plane == 0 ? mbmi->mode : mbmi->uv_mode;
        vp9_predict_intra_block(xd, block >> 2, bwl, TX_8X8, mode,
                                x->skip_encode ? src : dst,
                                x->skip_encode ? src_stride : dst_stride,
                                dst, dst_stride, i, j, plane);
        if (!x->skip_recode) {
          vp9_highbd_subtract_block(8, 8, src_diff, diff_stride,
                                    src, src_stride, dst, dst_stride, xd->bd);
          vp9_highbd_fht8x8(src_diff, coeff, diff_stride, tx_type);
          vp9_highbd_quantize_b(coeff, 64, x->skip_block, p->zbin, p->round,
                                p->quant, p->quant_shift, qcoeff, dqcoeff,
                                pd->dequant, eob,
                                scan_order->scan, scan_order->iscan);
        }
        if (!x->skip_encode && *eob) {
          vp9_highbd_iht8x8_add(tx_type, dqcoeff, dst, dst_stride, *eob,
                                xd->bd);
        }
        break;
      case TX_4X4:
        tx_type = get_tx_type_4x4(pd->plane_type, xd, block);
        scan_order = &vp9_scan_orders[TX_4X4][tx_type];
        mode = plane == 0 ? get_y_mode(xd->mi[0].src_mi, block) : mbmi->uv_mode;
        vp9_predict_intra_block(xd, block, bwl, TX_4X4, mode,
                                x->skip_encode ? src : dst,
                                x->skip_encode ? src_stride : dst_stride,
                                dst, dst_stride, i, j, plane);

        if (!x->skip_recode) {
          vp9_highbd_subtract_block(4, 4, src_diff, diff_stride,
                                    src, src_stride, dst, dst_stride, xd->bd);
          if (tx_type != DCT_DCT)
            vp9_highbd_fht4x4(src_diff, coeff, diff_stride, tx_type);
          else
            x->fwd_txm4x4(src_diff, coeff, diff_stride);
          vp9_highbd_quantize_b(coeff, 16, x->skip_block, p->zbin, p->round,
                                p->quant, p->quant_shift, qcoeff, dqcoeff,
                                pd->dequant, eob,
                                scan_order->scan, scan_order->iscan);
        }

        if (!x->skip_encode && *eob) {
          if (tx_type == DCT_DCT) {
            // this is like vp9_short_idct4x4 but has a special case around
            // eob<=1 which is significant (not just an optimization) for the
            // lossless case.
            x->highbd_itxm_add(dqcoeff, dst, dst_stride, *eob, xd->bd);
          } else {
            vp9_highbd_iht4x4_16_add(dqcoeff, dst, dst_stride, tx_type, xd->bd);
          }
        }
        break;
      default:
        assert(0);
        return;
    }
    if (*eob)
      *(args->skip) = 0;
    return;
  }
#endif  // CONFIG_VP9_HIGHBITDEPTH

  switch (tx_size) {
    case TX_32X32:
      scan_order = &vp9_default_scan_orders[TX_32X32];
      mode = plane == 0 ? mbmi->mode : mbmi->uv_mode;
      vp9_predict_intra_block(xd, block >> 6, bwl, TX_32X32, mode,
                              x->skip_encode ? src : dst,
                              x->skip_encode ? src_stride : dst_stride,
                              dst, dst_stride, i, j, plane);
      if (!x->skip_recode) {
        vp9_subtract_block(32, 32, src_diff, diff_stride,
                           src, src_stride, dst, dst_stride);
        fdct32x32(x->use_lp32x32fdct, src_diff, coeff, diff_stride);
        vp9_quantize_b_32x32(coeff, 1024, x->skip_block, p->zbin, p->round,
                             p->quant, p->quant_shift, qcoeff, dqcoeff,
                             pd->dequant, eob, scan_order->scan,
                             scan_order->iscan);
      }
      if (!x->skip_encode && *eob)
        vp9_idct32x32_add(dqcoeff, dst, dst_stride, *eob);
      break;
    case TX_16X16:
      tx_type = get_tx_type(pd->plane_type, xd);
      scan_order = &vp9_scan_orders[TX_16X16][tx_type];
      mode = plane == 0 ? mbmi->mode : mbmi->uv_mode;
      vp9_predict_intra_block(xd, block >> 4, bwl, TX_16X16, mode,
                              x->skip_encode ? src : dst,
                              x->skip_encode ? src_stride : dst_stride,
                              dst, dst_stride, i, j, plane);
      if (!x->skip_recode) {
        vp9_subtract_block(16, 16, src_diff, diff_stride,
                           src, src_stride, dst, dst_stride);
        vp9_fht16x16(src_diff, coeff, diff_stride, tx_type);
        vp9_quantize_b(coeff, 256, x->skip_block, p->zbin, p->round,
                       p->quant, p->quant_shift, qcoeff, dqcoeff,
                       pd->dequant, eob, scan_order->scan,
                       scan_order->iscan);
      }
      if (!x->skip_encode && *eob)
        vp9_iht16x16_add(tx_type, dqcoeff, dst, dst_stride, *eob);
      break;
    case TX_8X8:
      tx_type = get_tx_type(pd->plane_type, xd);
      scan_order = &vp9_scan_orders[TX_8X8][tx_type];
      mode = plane == 0 ? mbmi->mode : mbmi->uv_mode;
      vp9_predict_intra_block(xd, block >> 2, bwl, TX_8X8, mode,
                              x->skip_encode ? src : dst,
                              x->skip_encode ? src_stride : dst_stride,
                              dst, dst_stride, i, j, plane);
      if (!x->skip_recode) {
        vp9_subtract_block(8, 8, src_diff, diff_stride,
                           src, src_stride, dst, dst_stride);
        vp9_fht8x8(src_diff, coeff, diff_stride, tx_type);
        vp9_quantize_b(coeff, 64, x->skip_block, p->zbin, p->round, p->quant,
                       p->quant_shift, qcoeff, dqcoeff,
                       pd->dequant, eob, scan_order->scan,
                       scan_order->iscan);
      }
      if (!x->skip_encode && *eob)
        vp9_iht8x8_add(tx_type, dqcoeff, dst, dst_stride, *eob);
      break;
    case TX_4X4:
      tx_type = get_tx_type_4x4(pd->plane_type, xd, block);
      scan_order = &vp9_scan_orders[TX_4X4][tx_type];
      mode = plane == 0 ? get_y_mode(xd->mi[0].src_mi, block) : mbmi->uv_mode;
      vp9_predict_intra_block(xd, block, bwl, TX_4X4, mode,
                              x->skip_encode ? src : dst,
                              x->skip_encode ? src_stride : dst_stride,
                              dst, dst_stride, i, j, plane);

      if (!x->skip_recode) {
        vp9_subtract_block(4, 4, src_diff, diff_stride,
                           src, src_stride, dst, dst_stride);
        if (tx_type != DCT_DCT)
          vp9_fht4x4(src_diff, coeff, diff_stride, tx_type);
        else
          x->fwd_txm4x4(src_diff, coeff, diff_stride);
        vp9_quantize_b(coeff, 16, x->skip_block, p->zbin, p->round, p->quant,
                       p->quant_shift, qcoeff, dqcoeff,
                       pd->dequant, eob, scan_order->scan,
                       scan_order->iscan);
      }

      if (!x->skip_encode && *eob) {
        if (tx_type == DCT_DCT)
          // this is like vp9_short_idct4x4 but has a special case around eob<=1
          // which is significant (not just an optimization) for the lossless
          // case.
          x->itxm_add(dqcoeff, dst, dst_stride, *eob);
        else
          vp9_iht4x4_16_add(dqcoeff, dst, dst_stride, tx_type);
      }
      break;
    default:
      assert(0);
      break;
  }
  if (*eob)
    *(args->skip) = 0;
}

void vp9_encode_intra_block_plane(MACROBLOCK *x, BLOCK_SIZE bsize, int plane) {
  const MACROBLOCKD *const xd = &x->e_mbd;
  struct encode_b_args arg = {x, NULL, &xd->mi[0].src_mi->mbmi.skip};

  vp9_foreach_transformed_block_in_plane(xd, bsize, plane,
                                         vp9_encode_block_intra, &arg);
}
