/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */


#include "./vp10_rtcd.h"
#include "./vpx_config.h"
#include "./vpx_dsp_rtcd.h"

#include "vpx_dsp/quantize.h"
#include "vpx_mem/vpx_mem.h"
#include "vpx_ports/mem.h"

#include "vp10/common/idct.h"
#include "vp10/common/reconinter.h"
#include "vp10/common/reconintra.h"
#include "vp10/common/scan.h"

#include "vp10/encoder/encodemb.h"
#include "vp10/encoder/hybrid_fwd_txfm.h"
#include "vp10/encoder/quantize.h"
#include "vp10/encoder/rd.h"
#include "vp10/encoder/tokenize.h"

void vp10_subtract_plane(MACROBLOCK *x, BLOCK_SIZE bsize, int plane) {
  struct macroblock_plane *const p = &x->plane[plane];
  const struct macroblockd_plane *const pd = &x->e_mbd.plane[plane];
  const BLOCK_SIZE plane_bsize = get_plane_block_size(bsize, pd);
  const int bw = 4 * num_4x4_blocks_wide_lookup[plane_bsize];
  const int bh = 4 * num_4x4_blocks_high_lookup[plane_bsize];

#if CONFIG_VP9_HIGHBITDEPTH
  if (x->e_mbd.cur_buf->flags & YV12_FLAG_HIGHBITDEPTH) {
    vpx_highbd_subtract_block(bh, bw, p->src_diff, bw, p->src.buf,
                              p->src.stride, pd->dst.buf, pd->dst.stride,
                              x->e_mbd.bd);
    return;
  }
#endif  // CONFIG_VP9_HIGHBITDEPTH
  vpx_subtract_block(bh, bw, p->src_diff, bw, p->src.buf, p->src.stride,
                     pd->dst.buf, pd->dst.stride);
}

typedef struct vp10_token_state {
  int           rate;
  int64_t       error;
  int           next;
  int16_t       token;
  tran_low_t    qc;
  tran_low_t    dqc;
} vp10_token_state;

// These numbers are empirically obtained.
static const int plane_rd_mult[REF_TYPES][PLANE_TYPES] = {
    {10, 6}, {8, 5},
};

#define UPDATE_RD_COST()\
{\
  rd_cost0 = RDCOST(rdmult, rddiv, rate0, error0);\
  rd_cost1 = RDCOST(rdmult, rddiv, rate1, error1);\
}

int vp10_optimize_b(MACROBLOCK *mb, int plane, int block,
                    TX_SIZE tx_size, int ctx) {
  MACROBLOCKD *const xd = &mb->e_mbd;
  struct macroblock_plane *const p = &mb->plane[plane];
  struct macroblockd_plane *const pd = &xd->plane[plane];
  const int ref = is_inter_block(&xd->mi[0]->mbmi);
  vp10_token_state tokens[MAX_TX_SQUARE + 1][2];
  unsigned best_index[MAX_TX_SQUARE + 1][2];
  uint8_t token_cache[MAX_TX_SQUARE];
  const tran_low_t *const coeff = BLOCK_OFFSET(mb->plane[plane].coeff, block);
  tran_low_t *const qcoeff = BLOCK_OFFSET(p->qcoeff, block);
  tran_low_t *const dqcoeff = BLOCK_OFFSET(pd->dqcoeff, block);
  const int eob = p->eobs[block];
  const PLANE_TYPE type = pd->plane_type;
  const int default_eob = get_tx2d_size(tx_size);
  const int16_t* const dequant_ptr = pd->dequant;
  const uint8_t* const band_translate = get_band_translate(tx_size);
  TX_TYPE tx_type = get_tx_type(type, xd, block, tx_size);
  const scan_order* const so =
      get_scan(tx_size, tx_type, is_inter_block(&xd->mi[0]->mbmi));
  const int16_t* const scan = so->scan;
  const int16_t* const nb = so->neighbors;
#if CONFIG_NEW_QUANT
  int dq = get_dq_profile_from_ctx(ctx);
  const dequant_val_type_nuq *dequant_val = pd->dequant_val_nuq[dq];
#endif  // CONFIG_NEW_QUANT
  const int shift = get_tx_scale(xd, tx_type, tx_size);
  const int dq_step[2] = { dequant_ptr[0] >> shift, dequant_ptr[1] >> shift };
  int next = eob, sz = 0;
  const int64_t rdmult = (mb->rdmult * plane_rd_mult[ref][type]) >> 1;
  const int64_t rddiv = mb->rddiv;
  int64_t rd_cost0, rd_cost1;
  int rate0, rate1;
  int64_t error0, error1;
  int16_t t0, t1;
  int best, band = (eob < default_eob) ?
      band_translate[eob] : band_translate[eob - 1];
  int pt, i, final_eob;
#if CONFIG_VP9_HIGHBITDEPTH
  const int *cat6_high_cost = vp10_get_high_cost_table(xd->bd);
#else
  const int *cat6_high_cost = vp10_get_high_cost_table(8);
#endif
  unsigned int (*token_costs)[2][COEFF_CONTEXTS][ENTROPY_TOKENS] =
                   mb->token_costs[txsize_sqr_map[tx_size]][type][ref];
  const uint16_t *band_counts = &band_count_table[tx_size][band];
  uint16_t band_left = eob - band_cum_count_table[tx_size][band] + 1;
  int shortcut = 0;
  int next_shortcut = 0;

  token_costs += band;

  assert((!type && !plane) || (type && plane));
  assert(eob <= default_eob);

  /* Now set up a Viterbi trellis to evaluate alternative roundings. */
  /* Initialize the sentinel node of the trellis. */
  tokens[eob][0].rate = 0;
  tokens[eob][0].error = 0;
  tokens[eob][0].next = default_eob;
  tokens[eob][0].token = EOB_TOKEN;
  tokens[eob][0].qc = 0;
  tokens[eob][1] = tokens[eob][0];

  for (i = 0; i < eob; i++) {
    const int rc = scan[i];
    tokens[i][0].rate = vp10_get_token_cost(qcoeff[rc], &t0, cat6_high_cost);
    tokens[i][0].token = t0;
    token_cache[rc] = vp10_pt_energy_class[t0];
  }

  for (i = eob; i-- > 0;) {
    int base_bits, dx;
    int64_t d2;
    const int rc = scan[i];
    int x = qcoeff[rc];
    next_shortcut = shortcut;

    /* Only add a trellis state for non-zero coefficients. */
    if (UNLIKELY(x)) {
      error0 = tokens[next][0].error;
      error1 = tokens[next][1].error;
      /* Evaluate the first possibility for this state. */
      rate0 = tokens[next][0].rate;
      rate1 = tokens[next][1].rate;

      if (next_shortcut) {
        /* Consider both possible successor states. */
        if (next < default_eob) {
          pt = get_coef_context(nb, token_cache, i + 1);
          rate0 += (*token_costs)[0][pt][tokens[next][0].token];
          rate1 += (*token_costs)[0][pt][tokens[next][1].token];
        }
        UPDATE_RD_COST();
        /* And pick the best. */
        best = rd_cost1 < rd_cost0;
      } else {
        if (next < default_eob) {
          pt = get_coef_context(nb, token_cache, i + 1);
          rate0 += (*token_costs)[0][pt][tokens[next][0].token];
        }
        best = 0;
      }

      dx = (dqcoeff[rc] - coeff[rc]) * (1 << shift);
#if CONFIG_VP9_HIGHBITDEPTH
      if (xd->cur_buf->flags & YV12_FLAG_HIGHBITDEPTH) {
        dx >>= xd->bd - 8;
      }
#endif  // CONFIG_VP9_HIGHBITDEPTH
      d2 = (int64_t)dx * dx;
      tokens[i][0].rate += (best ? rate1 : rate0);
      tokens[i][0].error = d2 + (best ? error1 : error0);
      tokens[i][0].next = next;
      tokens[i][0].qc = x;
      tokens[i][0].dqc = dqcoeff[rc];
      best_index[i][0] = best;

      /* Evaluate the second possibility for this state. */
      rate0 = tokens[next][0].rate;
      rate1 = tokens[next][1].rate;

      // The threshold of 3 is empirically obtained.
      if (UNLIKELY(abs(x) > 3)) {
        shortcut = 0;
      } else {
#if CONFIG_NEW_QUANT
        shortcut = (
            (vp10_dequant_abscoeff_nuq(
                abs(x), dequant_ptr[rc != 0],
                dequant_val[band_translate[i]]) > (abs(coeff[rc]) << shift)) &&
            (vp10_dequant_abscoeff_nuq(
                abs(x) - 1, dequant_ptr[rc != 0],
                dequant_val[band_translate[i]]) < (abs(coeff[rc]) << shift)));
#else   // CONFIG_NEW_QUANT
        if ((abs(x) * dequant_ptr[rc != 0] > (abs(coeff[rc]) << shift)) &&
            (abs(x) * dequant_ptr[rc != 0] < (abs(coeff[rc]) << shift) +
              dequant_ptr[rc != 0]))
          shortcut = 1;
        else
          shortcut = 0;
#endif   // CONFIG_NEW_QUANT
      }

      if (shortcut) {
        sz = -(x < 0);
        x -= 2 * sz + 1;
      } else {
        tokens[i][1] = tokens[i][0];
        best_index[i][1] = best_index[i][0];
        next = i;

        if (UNLIKELY(!(--band_left))) {
          --band_counts;
          band_left = *band_counts;
          --token_costs;
        }
        continue;
      }

      /* Consider both possible successor states. */
      if (!x) {
        /* If we reduced this coefficient to zero, check to see if
         *  we need to move the EOB back here.
         */
        t0 = tokens[next][0].token == EOB_TOKEN ? EOB_TOKEN : ZERO_TOKEN;
        t1 = tokens[next][1].token == EOB_TOKEN ? EOB_TOKEN : ZERO_TOKEN;
        base_bits = 0;
      } else {
        base_bits = vp10_get_token_cost(x, &t0, cat6_high_cost);
        t1 = t0;
      }

      if (next_shortcut) {
        if (LIKELY(next < default_eob)) {
          if (t0 != EOB_TOKEN) {
            token_cache[rc] = vp10_pt_energy_class[t0];
            pt = get_coef_context(nb, token_cache, i + 1);
            rate0 += (*token_costs)[!x][pt][tokens[next][0].token];
          }
          if (t1 != EOB_TOKEN) {
            token_cache[rc] = vp10_pt_energy_class[t1];
            pt = get_coef_context(nb, token_cache, i + 1);
            rate1 += (*token_costs)[!x][pt][tokens[next][1].token];
          }
        }

        UPDATE_RD_COST();
        /* And pick the best. */
        best = rd_cost1 < rd_cost0;
      } else {
        // The two states in next stage are identical.
        if (next < default_eob && t0 != EOB_TOKEN) {
          token_cache[rc] = vp10_pt_energy_class[t0];
          pt = get_coef_context(nb, token_cache, i + 1);
          rate0 += (*token_costs)[!x][pt][tokens[next][0].token];
        }
        best = 0;
      }

#if CONFIG_NEW_QUANT
      dx = vp10_dequant_coeff_nuq(
          x, dequant_ptr[rc != 0],
          dequant_val[band_translate[i]]) - (coeff[rc] << shift);
#if CONFIG_VP9_HIGHBITDEPTH
      if (xd->cur_buf->flags & YV12_FLAG_HIGHBITDEPTH) {
        dx >>= xd->bd - 8;
      }
#endif  // CONFIG_VP9_HIGHBITDEPTH
#else   // CONFIG_NEW_QUANT
#if CONFIG_VP9_HIGHBITDEPTH
      if (xd->cur_buf->flags & YV12_FLAG_HIGHBITDEPTH) {
        dx -= ((dequant_ptr[rc != 0] >> (xd->bd - 8)) + sz) ^ sz;
      } else {
        dx -= (dequant_ptr[rc != 0] + sz) ^ sz;
      }
#else
      dx -= (dequant_ptr[rc != 0] + sz) ^ sz;
#endif  // CONFIG_VP9_HIGHBITDEPTH
#endif  // CONFIG_NEW_QUANT
      d2 = (int64_t)dx * dx;

      tokens[i][1].rate = base_bits + (best ? rate1 : rate0);
      tokens[i][1].error = d2 + (best ? error1 : error0);
      tokens[i][1].next = next;
      tokens[i][1].token = best ? t1 : t0;
      tokens[i][1].qc = x;

      if (x) {
        tran_low_t offset = dq_step[rc != 0];
        // The 32x32 transform coefficient uses half quantization step size.
        // Account for the rounding difference in the dequantized coefficeint
        // value when the quantization index is dropped from an even number
        // to an odd number.
        if (shift & x)
          offset += (dequant_ptr[rc != 0] & 0x01);

        if (sz == 0)
          tokens[i][1].dqc = dqcoeff[rc] - offset;
        else
          tokens[i][1].dqc = dqcoeff[rc] + offset;
      } else {
        tokens[i][1].dqc = 0;
      }

      best_index[i][1] = best;
      /* Finally, make this the new head of the trellis. */
      next = i;
    } else {
      /* There's no choice to make for a zero coefficient, so we don't
       *  add a new trellis node, but we do need to update the costs.
       */
      t0 = tokens[next][0].token;
      t1 = tokens[next][1].token;
      pt = get_coef_context(nb, token_cache, i + 1);
      /* Update the cost of each path if we're past the EOB token. */
      if (t0 != EOB_TOKEN) {
        tokens[next][0].rate += (*token_costs)[1][pt][t0];
        tokens[next][0].token = ZERO_TOKEN;
      }
      if (t1 != EOB_TOKEN) {
        tokens[next][1].rate += (*token_costs)[1][pt][t1];
        tokens[next][1].token = ZERO_TOKEN;
      }
      best_index[i][0] = best_index[i][1] = 0;
      shortcut = (tokens[next][0].rate != tokens[next][1].rate);
      /* Don't update next, because we didn't add a new node. */
    }

    if (UNLIKELY(!(--band_left))) {
      --band_counts;
      band_left = *band_counts;
      --token_costs;
    }
  }

  /* Now pick the best path through the whole trellis. */
  rate0 = tokens[next][0].rate;
  rate1 = tokens[next][1].rate;
  error0 = tokens[next][0].error;
  error1 = tokens[next][1].error;
  t0 = tokens[next][0].token;
  t1 = tokens[next][1].token;
  rate0 += (*token_costs)[0][ctx][t0];
  rate1 += (*token_costs)[0][ctx][t1];
  UPDATE_RD_COST();
  best = rd_cost1 < rd_cost0;

  final_eob = -1;

  for (i = next; i < eob; i = next) {
    const int x = tokens[i][best].qc;
    const int rc = scan[i];

    if (x) final_eob = i;
    qcoeff[rc] = x;
    dqcoeff[rc] = tokens[i][best].dqc;

#if CONFIG_NEW_QUANT
    dqcoeff[rc] = vp10_dequant_abscoeff_nuq(abs(x), dequant_ptr[rc != 0],
                                            dequant_val[band_translate[i]]);
    if (shift) dqcoeff[rc] = ROUND_POWER_OF_TWO(dqcoeff[rc], shift);
    if (x < 0) dqcoeff[rc] = -dqcoeff[rc];
#endif  // CONFIG_NEW_QUANT

    next = tokens[i][best].next;
    best = best_index[i][best];
  }
  final_eob++;

  mb->plane[plane].eobs[block] = final_eob;
  assert(final_eob <= default_eob);
  return final_eob;
}

#if CONFIG_VP9_HIGHBITDEPTH
typedef enum QUANT_FUNC {
  QUANT_FUNC_LOWBD = 0,
  QUANT_FUNC_HIGHBD = 1,
  QUANT_FUNC_LAST = 2
} QUANT_FUNC;

static VP10_QUANT_FACADE
    quant_func_list[VP10_XFORM_QUANT_LAST][QUANT_FUNC_LAST] = {
        {vp10_quantize_fp_facade, vp10_highbd_quantize_fp_facade},
        {vp10_quantize_b_facade, vp10_highbd_quantize_b_facade},
        {vp10_quantize_dc_facade, vp10_highbd_quantize_dc_facade},
        {NULL, NULL}};

#else
typedef enum QUANT_FUNC {
  QUANT_FUNC_LOWBD = 0,
  QUANT_FUNC_LAST = 1
} QUANT_FUNC;

static VP10_QUANT_FACADE
    quant_func_list[VP10_XFORM_QUANT_LAST][QUANT_FUNC_LAST] = {
        {vp10_quantize_fp_facade},
        {vp10_quantize_b_facade},
        {vp10_quantize_dc_facade},
        {NULL}};
#endif

static FWD_TXFM_OPT fwd_txfm_opt_list[VP10_XFORM_QUANT_LAST] = {
    FWD_TXFM_OPT_NORMAL, FWD_TXFM_OPT_NORMAL, FWD_TXFM_OPT_DC,
    FWD_TXFM_OPT_NORMAL};

void vp10_xform_quant(MACROBLOCK *x, int plane, int block, int blk_row,
                      int blk_col, BLOCK_SIZE plane_bsize, TX_SIZE tx_size,
                      VP10_XFORM_QUANT xform_quant_idx) {
  MACROBLOCKD *const xd = &x->e_mbd;
  const struct macroblock_plane *const p = &x->plane[plane];
  const struct macroblockd_plane *const pd = &xd->plane[plane];
  PLANE_TYPE plane_type = (plane == 0) ? PLANE_TYPE_Y : PLANE_TYPE_UV;
  TX_TYPE tx_type = get_tx_type(plane_type, xd, block, tx_size);
  const scan_order *const scan_order =
      get_scan(tx_size, tx_type, is_inter_block(&xd->mi[0]->mbmi));
  tran_low_t *const coeff = BLOCK_OFFSET(p->coeff, block);
  tran_low_t *const qcoeff = BLOCK_OFFSET(p->qcoeff, block);
  tran_low_t *const dqcoeff = BLOCK_OFFSET(pd->dqcoeff, block);
  uint16_t *const eob = &p->eobs[block];
  const int diff_stride = 4 * num_4x4_blocks_wide_lookup[plane_bsize];
  const int16_t *src_diff;
  const int tx2d_size = get_tx2d_size(tx_size);

  FWD_TXFM_PARAM fwd_txfm_param;
  QUANT_PARAM qparam;

  fwd_txfm_param.tx_type = tx_type;
  fwd_txfm_param.tx_size = tx_size;
  fwd_txfm_param.fwd_txfm_opt = fwd_txfm_opt_list[xform_quant_idx];
  fwd_txfm_param.rd_transform = x->use_lp32x32fdct;
  fwd_txfm_param.lossless = xd->lossless[xd->mi[0]->mbmi.segment_id];

  src_diff = &p->src_diff[4 * (blk_row * diff_stride + blk_col)];

  qparam.log_scale = get_tx_scale(xd, tx_type, tx_size);
#if CONFIG_VP9_HIGHBITDEPTH
  fwd_txfm_param.bd = xd->bd;
  if (xd->cur_buf->flags & YV12_FLAG_HIGHBITDEPTH) {
    highbd_fwd_txfm(src_diff, coeff, diff_stride, &fwd_txfm_param);
    if (xform_quant_idx != VP10_XFORM_QUANT_SKIP_QUANT) {
      if (x->skip_block) {
        vp10_quantize_skip(tx2d_size, qcoeff, dqcoeff, eob);
      } else {
        quant_func_list[xform_quant_idx][QUANT_FUNC_HIGHBD](
            coeff, tx2d_size, p, qcoeff, pd, dqcoeff, eob,
            scan_order, &qparam);
      }
    }
    return;
  }
#endif  // CONFIG_VP9_HIGHBITDEPTH

  fwd_txfm(src_diff, coeff, diff_stride, &fwd_txfm_param);
  if (xform_quant_idx != VP10_XFORM_QUANT_SKIP_QUANT) {
    if (x->skip_block) {
      vp10_quantize_skip(tx2d_size, qcoeff, dqcoeff, eob);
    } else {
      quant_func_list[xform_quant_idx][QUANT_FUNC_LOWBD](
          coeff, tx2d_size, p, qcoeff, pd, dqcoeff, eob,
          scan_order, &qparam);
    }
  }
}

#if CONFIG_NEW_QUANT
void vp10_xform_quant_nuq(MACROBLOCK *x, int plane, int block, int blk_row,
                          int blk_col, BLOCK_SIZE plane_bsize,
                          TX_SIZE tx_size, int ctx) {
  MACROBLOCKD *const xd = &x->e_mbd;
  const struct macroblock_plane *const p = &x->plane[plane];
  const struct macroblockd_plane *const pd = &xd->plane[plane];
  PLANE_TYPE plane_type = (plane == 0) ? PLANE_TYPE_Y : PLANE_TYPE_UV;
  TX_TYPE tx_type = get_tx_type(plane_type, xd, block, tx_size);
  const scan_order *const scan_order =
      get_scan(tx_size, tx_type, is_inter_block(&xd->mi[0]->mbmi));
  tran_low_t *const coeff = BLOCK_OFFSET(p->coeff, block);
  tran_low_t *const qcoeff = BLOCK_OFFSET(p->qcoeff, block);
  tran_low_t *const dqcoeff = BLOCK_OFFSET(pd->dqcoeff, block);
  int dq = get_dq_profile_from_ctx(ctx);
  uint16_t *const eob = &p->eobs[block];
  const int diff_stride = 4 * num_4x4_blocks_wide_lookup[plane_bsize];
  const int16_t *src_diff;
  const uint8_t* band = get_band_translate(tx_size);

  FWD_TXFM_PARAM fwd_txfm_param;

  fwd_txfm_param.tx_type = tx_type;
  fwd_txfm_param.tx_size = tx_size;
  fwd_txfm_param.fwd_txfm_opt = fwd_txfm_opt_list[VP10_XFORM_QUANT_FP];
  fwd_txfm_param.rd_transform = x->use_lp32x32fdct;
  fwd_txfm_param.lossless = xd->lossless[xd->mi[0]->mbmi.segment_id];

  src_diff = &p->src_diff[4 * (blk_row * diff_stride + blk_col)];

// TODO(sarahparker) add all of these new quant quantize functions
// to quant_func_list, just trying to get this expr to work for now
#if CONFIG_VP9_HIGHBITDEPTH
  fwd_txfm_param.bd = xd->bd;
  if (xd->cur_buf->flags & YV12_FLAG_HIGHBITDEPTH) {
    highbd_fwd_txfm(src_diff, coeff, diff_stride, &fwd_txfm_param);
    if (tx_size == TX_32X32) {
      highbd_quantize_32x32_nuq(coeff, get_tx2d_size(tx_size), x->skip_block,
                                p->quant, p->quant_shift, pd->dequant,
                                (const cuml_bins_type_nuq *)
                                p->cuml_bins_nuq[dq],
                                (const dequant_val_type_nuq *)
                                pd->dequant_val_nuq[dq],
                                qcoeff, dqcoeff, eob,
                                scan_order->scan, band);
    } else {
      highbd_quantize_nuq(coeff, get_tx2d_size(tx_size), x->skip_block,
                          p->quant, p->quant_shift, pd->dequant,
                          (const cuml_bins_type_nuq *)p->cuml_bins_nuq[dq],
                          (const dequant_val_type_nuq *)
                          pd->dequant_val_nuq[dq],
                          qcoeff, dqcoeff, eob,
                          scan_order->scan, band);
    }
    return;
  }
#endif  // CONFIG_VP9_HIGHBITDEPTH

  fwd_txfm(src_diff, coeff, diff_stride, &fwd_txfm_param);
  if (tx_size == TX_32X32) {
    quantize_32x32_nuq(coeff, 1024, x->skip_block,
                       p->quant, p->quant_shift, pd->dequant,
                       (const cuml_bins_type_nuq *)p->cuml_bins_nuq[dq],
                       (const dequant_val_type_nuq *)
                       pd->dequant_val_nuq[dq],
                       qcoeff, dqcoeff, eob,
                       scan_order->scan, band);
  } else {
    quantize_nuq(coeff, get_tx2d_size(tx_size), x->skip_block,
                 p->quant, p->quant_shift, pd->dequant,
                 (const cuml_bins_type_nuq *)p->cuml_bins_nuq[dq],
                 (const dequant_val_type_nuq *)pd->dequant_val_nuq[dq],
                 qcoeff, dqcoeff, eob,
                 scan_order->scan, band);
  }
}

void vp10_xform_quant_fp_nuq(MACROBLOCK *x, int plane, int block, int blk_row,
                             int blk_col, BLOCK_SIZE plane_bsize,
                             TX_SIZE tx_size, int ctx) {
  MACROBLOCKD *const xd = &x->e_mbd;
  const struct macroblock_plane *const p = &x->plane[plane];
  const struct macroblockd_plane *const pd = &xd->plane[plane];
  int dq = get_dq_profile_from_ctx(ctx);
  PLANE_TYPE plane_type = (plane == 0) ? PLANE_TYPE_Y : PLANE_TYPE_UV;
  TX_TYPE tx_type = get_tx_type(plane_type, xd, block, tx_size);
  const scan_order *const scan_order =
      get_scan(tx_size, tx_type, is_inter_block(&xd->mi[0]->mbmi));
  tran_low_t *const coeff = BLOCK_OFFSET(p->coeff, block);
  tran_low_t *const qcoeff = BLOCK_OFFSET(p->qcoeff, block);
  tran_low_t *const dqcoeff = BLOCK_OFFSET(pd->dqcoeff, block);
  uint16_t *const eob = &p->eobs[block];
  const int diff_stride = 4 * num_4x4_blocks_wide_lookup[plane_bsize];
  const int16_t *src_diff;
  const uint8_t* band = get_band_translate(tx_size);

  FWD_TXFM_PARAM fwd_txfm_param;

  fwd_txfm_param.tx_type = tx_type;
  fwd_txfm_param.tx_size = tx_size;
  fwd_txfm_param.fwd_txfm_opt = fwd_txfm_opt_list[VP10_XFORM_QUANT_FP];
  fwd_txfm_param.rd_transform = x->use_lp32x32fdct;
  fwd_txfm_param.lossless = xd->lossless[xd->mi[0]->mbmi.segment_id];

  src_diff = &p->src_diff[4 * (blk_row * diff_stride + blk_col)];

// TODO(sarahparker) add all of these new quant quantize functions
// to quant_func_list, just trying to get this expr to work for now
#if CONFIG_VP9_HIGHBITDEPTH
  fwd_txfm_param.bd = xd->bd;
  if (xd->cur_buf->flags & YV12_FLAG_HIGHBITDEPTH) {
    highbd_fwd_txfm(src_diff, coeff, diff_stride, &fwd_txfm_param);
    if (tx_size == TX_32X32) {
      highbd_quantize_32x32_fp_nuq(coeff, get_tx2d_size(tx_size), x->skip_block,
                                   p->quant_fp, pd->dequant,
                                   (const cuml_bins_type_nuq *)
                                   p->cuml_bins_nuq[dq],
                                   (const dequant_val_type_nuq *)
                                   pd->dequant_val_nuq[dq],
                                   qcoeff, dqcoeff, eob,
                                   scan_order->scan, band);
    } else {
      highbd_quantize_fp_nuq(coeff, get_tx2d_size(tx_size), x->skip_block,
                             p->quant_fp, pd->dequant,
                             (const cuml_bins_type_nuq *)
                             p->cuml_bins_nuq[dq],
                             (const dequant_val_type_nuq *)
                             pd->dequant_val_nuq[dq],
                             qcoeff, dqcoeff, eob,
                             scan_order->scan, band);
    }
    return;
  }
#endif  // CONFIG_VP9_HIGHBITDEPTH

  fwd_txfm(src_diff, coeff, diff_stride, &fwd_txfm_param);
  if (tx_size == TX_32X32) {
    quantize_32x32_fp_nuq(coeff, get_tx2d_size(tx_size), x->skip_block,
                          p->quant_fp, pd->dequant,
                          (const cuml_bins_type_nuq *)
                          p->cuml_bins_nuq[dq],
                          (const dequant_val_type_nuq *)
                          pd->dequant_val_nuq[dq],
                          qcoeff, dqcoeff, eob,
                          scan_order->scan, band);
  } else {
    quantize_fp_nuq(coeff, get_tx2d_size(tx_size), x->skip_block,
                    p->quant_fp, pd->dequant,
                    (const cuml_bins_type_nuq *)
                    p->cuml_bins_nuq[dq],
                    (const dequant_val_type_nuq *)
                    pd->dequant_val_nuq[dq],
                    qcoeff, dqcoeff, eob,
                    scan_order->scan, band);
  }
}

void vp10_xform_quant_dc_nuq(MACROBLOCK *x, int plane, int block, int blk_row,
                             int blk_col, BLOCK_SIZE plane_bsize,
                             TX_SIZE tx_size, int ctx) {
  MACROBLOCKD *const xd = &x->e_mbd;
  const struct macroblock_plane *const p = &x->plane[plane];
  const struct macroblockd_plane *const pd = &xd->plane[plane];
  PLANE_TYPE plane_type = (plane == 0) ? PLANE_TYPE_Y : PLANE_TYPE_UV;
  TX_TYPE tx_type = get_tx_type(plane_type, xd, block, tx_size);
  tran_low_t *const coeff = BLOCK_OFFSET(p->coeff, block);
  tran_low_t *const qcoeff = BLOCK_OFFSET(p->qcoeff, block);
  tran_low_t *const dqcoeff = BLOCK_OFFSET(pd->dqcoeff, block);
  uint16_t *const eob = &p->eobs[block];
  const int diff_stride = 4 * num_4x4_blocks_wide_lookup[plane_bsize];
  const int16_t *src_diff;
  int dq = get_dq_profile_from_ctx(ctx);

  FWD_TXFM_PARAM fwd_txfm_param;

  fwd_txfm_param.tx_type = tx_type;
  fwd_txfm_param.tx_size = tx_size;
  fwd_txfm_param.fwd_txfm_opt = fwd_txfm_opt_list[VP10_XFORM_QUANT_DC];
  fwd_txfm_param.rd_transform = x->use_lp32x32fdct;
  fwd_txfm_param.lossless = xd->lossless[xd->mi[0]->mbmi.segment_id];

  src_diff = &p->src_diff[4 * (blk_row * diff_stride + blk_col)];

// TODO(sarahparker) add all of these new quant quantize functions
// to quant_func_list, just trying to get this expr to work for now
#if CONFIG_VP9_HIGHBITDEPTH
  fwd_txfm_param.bd = xd->bd;
  if (xd->cur_buf->flags & YV12_FLAG_HIGHBITDEPTH) {
    highbd_fwd_txfm(src_diff, coeff, diff_stride, &fwd_txfm_param);
    if (tx_size == TX_32X32) {
      highbd_quantize_dc_32x32_nuq(coeff, get_tx2d_size(tx_size), x->skip_block,
                                   p->quant[0], p->quant_shift[0],
                                   pd->dequant[0],
                                   p->cuml_bins_nuq[dq][0],
                                   pd->dequant_val_nuq[dq][0],
                                   qcoeff, dqcoeff, eob);
    } else {
      highbd_quantize_dc_nuq(coeff, get_tx2d_size(tx_size), x->skip_block,
                             p->quant[0], p->quant_shift[0],
                             pd->dequant[0],
                             p->cuml_bins_nuq[dq][0],
                             pd->dequant_val_nuq[dq][0],
                             qcoeff, dqcoeff, eob);
    }
    return;
  }
#endif  // CONFIG_VP9_HIGHBITDEPTH

  fwd_txfm(src_diff, coeff, diff_stride, &fwd_txfm_param);
  if (tx_size == TX_32X32) {
    quantize_dc_32x32_nuq(coeff, get_tx2d_size(tx_size), x->skip_block,
                          p->quant[0], p->quant_shift[0], pd->dequant[0],
                          p->cuml_bins_nuq[dq][0],
                          pd->dequant_val_nuq[dq][0],
                          qcoeff, dqcoeff, eob);
  } else {
    quantize_dc_nuq(coeff, get_tx2d_size(tx_size), x->skip_block,
                    p->quant[0], p->quant_shift[0], pd->dequant[0],
                    p->cuml_bins_nuq[dq][0],
                    pd->dequant_val_nuq[dq][0],
                    qcoeff, dqcoeff, eob);
  }
}

void vp10_xform_quant_dc_fp_nuq(MACROBLOCK *x, int plane, int block,
                                int blk_row, int blk_col,
                                BLOCK_SIZE plane_bsize, TX_SIZE tx_size,
                                int ctx) {
  MACROBLOCKD *const xd = &x->e_mbd;
  const struct macroblock_plane *const p = &x->plane[plane];
  const struct macroblockd_plane *const pd = &xd->plane[plane];
  PLANE_TYPE plane_type = (plane == 0) ? PLANE_TYPE_Y : PLANE_TYPE_UV;
  TX_TYPE tx_type = get_tx_type(plane_type, xd, block, tx_size);
  tran_low_t *const coeff = BLOCK_OFFSET(p->coeff, block);
  tran_low_t *const qcoeff = BLOCK_OFFSET(p->qcoeff, block);
  tran_low_t *const dqcoeff = BLOCK_OFFSET(pd->dqcoeff, block);
  uint16_t *const eob = &p->eobs[block];
  const int diff_stride = 4 * num_4x4_blocks_wide_lookup[plane_bsize];
  const int16_t *src_diff;
  int dq = get_dq_profile_from_ctx(ctx);

  FWD_TXFM_PARAM fwd_txfm_param;

  fwd_txfm_param.tx_type = tx_type;
  fwd_txfm_param.tx_size = tx_size;
  fwd_txfm_param.fwd_txfm_opt = fwd_txfm_opt_list[VP10_XFORM_QUANT_DC];
  fwd_txfm_param.rd_transform = x->use_lp32x32fdct;
  fwd_txfm_param.lossless = xd->lossless[xd->mi[0]->mbmi.segment_id];

  src_diff = &p->src_diff[4 * (blk_row * diff_stride + blk_col)];

// TODO(sarahparker) add all of these new quant quantize functions
// to quant_func_list, just trying to get this expr to work for now
#if CONFIG_VP9_HIGHBITDEPTH
  fwd_txfm_param.bd = xd->bd;
  if (xd->cur_buf->flags & YV12_FLAG_HIGHBITDEPTH) {
    highbd_fwd_txfm(src_diff, coeff, diff_stride, &fwd_txfm_param);
    if (tx_size == TX_32X32) {
      highbd_quantize_dc_32x32_fp_nuq(coeff, get_tx2d_size(tx_size),
                                      x->skip_block,
                                      p->quant_fp[0], pd->dequant[0],
                                      p->cuml_bins_nuq[dq][0],
                                      pd->dequant_val_nuq[dq][0],
                                      qcoeff, dqcoeff, eob);
    } else {
      highbd_quantize_dc_fp_nuq(coeff, get_tx2d_size(tx_size), x->skip_block,
                                p->quant_fp[0], pd->dequant[0],
                                p->cuml_bins_nuq[dq][0],
                                pd->dequant_val_nuq[dq][0],
                                qcoeff, dqcoeff, eob);
    }
    return;
  }
#endif  // CONFIG_VP9_HIGHBITDEPTH

  fwd_txfm(src_diff, coeff, diff_stride, &fwd_txfm_param);
  if (tx_size == TX_32X32) {
    quantize_dc_32x32_fp_nuq(coeff, get_tx2d_size(tx_size), x->skip_block,
                             p->quant_fp[0], pd->dequant[0],
                             p->cuml_bins_nuq[dq][0],
                             pd->dequant_val_nuq[dq][0],
                             qcoeff, dqcoeff, eob);
  } else {
    quantize_dc_fp_nuq(coeff, get_tx2d_size(tx_size), x->skip_block,
                       p->quant_fp[0], pd->dequant[0],
                       p->cuml_bins_nuq[dq][0],
                       pd->dequant_val_nuq[dq][0],
                       qcoeff, dqcoeff, eob);
  }
}
#endif  // CONFIG_NEW_QUANT

static void encode_block(int plane, int block, int blk_row, int blk_col,
                         BLOCK_SIZE plane_bsize,
                         TX_SIZE tx_size, void *arg) {
  struct encode_b_args *const args = arg;
  MACROBLOCK *const x = args->x;
  MACROBLOCKD *const xd = &x->e_mbd;
  int ctx;
  struct macroblock_plane *const p = &x->plane[plane];
  struct macroblockd_plane *const pd = &xd->plane[plane];
  tran_low_t *const dqcoeff = BLOCK_OFFSET(pd->dqcoeff, block);
  uint8_t *dst;
  ENTROPY_CONTEXT *a, *l;
  INV_TXFM_PARAM inv_txfm_param;
#if CONFIG_VAR_TX
  int i;
  const int bwl = b_width_log2_lookup[plane_bsize];
#endif
  dst = &pd->dst.buf[4 * blk_row * pd->dst.stride + 4 * blk_col];
  a = &args->ta[blk_col];
  l = &args->tl[blk_row];
#if CONFIG_VAR_TX
  ctx = get_entropy_context(tx_size, a, l);
#else
  ctx = combine_entropy_contexts(*a, *l);
#endif

#if CONFIG_VAR_TX
  // Assert not magic number (uninitialised).
  assert(x->blk_skip[plane][(blk_row << bwl) + blk_col] != 234);

  if (x->blk_skip[plane][(blk_row << bwl) + blk_col] == 0) {
#else
  {
#endif
#if CONFIG_NEW_QUANT
    vp10_xform_quant_fp_nuq(x, plane, block, blk_row, blk_col, plane_bsize,
                            tx_size, ctx);
#else
    vp10_xform_quant(x, plane, block, blk_row, blk_col, plane_bsize,
                     tx_size, VP10_XFORM_QUANT_FP);
#endif  // CONFIG_NEW_QUANT
  }
#if CONFIG_VAR_TX
  else {
    p->eobs[block] = 0;
  }
#endif

  if (p->eobs[block]) {
    *a = *l = vp10_optimize_b(x, plane, block, tx_size, ctx) > 0;
  } else {
    *a = *l = p->eobs[block] > 0;
  }

#if CONFIG_VAR_TX
  for (i = 0; i < num_4x4_blocks_wide_txsize_lookup[tx_size]; ++i) {
    a[i] = a[0];
  }
  for (i = 0; i < num_4x4_blocks_high_txsize_lookup[tx_size]; ++i) {
    l[i] = l[0];
  }
#endif

  if (p->eobs[block])
    *(args->skip) = 0;

  if (p->eobs[block] == 0)
    return;

  // inverse transform parameters
  inv_txfm_param.tx_type = get_tx_type(pd->plane_type, xd, block, tx_size);
  inv_txfm_param.tx_size = tx_size;
  inv_txfm_param.eob = p->eobs[block];
  inv_txfm_param.lossless = xd->lossless[xd->mi[0]->mbmi.segment_id];

#if CONFIG_VP9_HIGHBITDEPTH
  if (xd->cur_buf->flags & YV12_FLAG_HIGHBITDEPTH) {
    inv_txfm_param.bd = xd->bd;
    highbd_inv_txfm_add(dqcoeff, dst, pd->dst.stride, &inv_txfm_param);
    return;
  }
#endif  // CONFIG_VP9_HIGHBITDEPTH
  inv_txfm_add(dqcoeff, dst, pd->dst.stride, &inv_txfm_param);
}

#if CONFIG_VAR_TX
static void encode_block_inter(int plane, int block, int blk_row, int blk_col,
                               BLOCK_SIZE plane_bsize, TX_SIZE tx_size,
                               void *arg) {
  struct encode_b_args *const args = arg;
  MACROBLOCK *const x = args->x;
  MACROBLOCKD *const xd = &x->e_mbd;
  MB_MODE_INFO *const mbmi = &xd->mi[0]->mbmi;
  const BLOCK_SIZE bsize = txsize_to_bsize[tx_size];
  const struct macroblockd_plane *const pd = &xd->plane[plane];
  const int tx_row = blk_row >> (1 - pd->subsampling_y);
  const int tx_col = blk_col >> (1 - pd->subsampling_x);
  const TX_SIZE plane_tx_size = plane ?
      get_uv_tx_size_impl(mbmi->inter_tx_size[tx_row][tx_col], bsize, 0, 0) :
      mbmi->inter_tx_size[tx_row][tx_col];

  int max_blocks_high = num_4x4_blocks_high_lookup[plane_bsize];
  int max_blocks_wide = num_4x4_blocks_wide_lookup[plane_bsize];

  if (xd->mb_to_bottom_edge < 0)
    max_blocks_high += xd->mb_to_bottom_edge >> (5 + pd->subsampling_y);
  if (xd->mb_to_right_edge < 0)
    max_blocks_wide += xd->mb_to_right_edge >> (5 + pd->subsampling_x);

  if (blk_row >= max_blocks_high || blk_col >= max_blocks_wide)
    return;

  if (tx_size == plane_tx_size) {
    encode_block(plane, block, blk_row, blk_col, plane_bsize,
                 tx_size, arg);
  } else {
    int bsl = b_width_log2_lookup[bsize];
    int i;

    assert(bsl > 0);
    --bsl;

#if CONFIG_EXT_TX
    assert(tx_size < TX_SIZES);
#endif  // CONFIG_EXT_TX

    for (i = 0; i < 4; ++i) {
      const int offsetr = blk_row + ((i >> 1) << bsl);
      const int offsetc = blk_col + ((i & 0x01) << bsl);
      int step = num_4x4_blocks_txsize_lookup[tx_size - 1];

      if (offsetr >= max_blocks_high || offsetc >= max_blocks_wide)
        continue;

      encode_block_inter(plane, block + i * step, offsetr, offsetc,
                         plane_bsize, tx_size - 1, arg);
    }
  }
}
#endif

static void encode_block_pass1(int plane, int block, int blk_row, int blk_col,
                               BLOCK_SIZE plane_bsize,
                               TX_SIZE tx_size, void *arg) {
  MACROBLOCK *const x = (MACROBLOCK *)arg;
  MACROBLOCKD *const xd = &x->e_mbd;
  struct macroblock_plane *const p = &x->plane[plane];
  struct macroblockd_plane *const pd = &xd->plane[plane];
  tran_low_t *const dqcoeff = BLOCK_OFFSET(pd->dqcoeff, block);
  uint8_t *dst;
#if CONFIG_NEW_QUANT
  int ctx;
#endif  // CONFIG_NEW_QUANT
  dst = &pd->dst.buf[4 * blk_row * pd->dst.stride + 4 * blk_col];

#if CONFIG_NEW_QUANT
  ctx = 0;
  vp10_xform_quant_fp_nuq(x, plane, block, blk_row, blk_col, plane_bsize,
                          tx_size, ctx);
#else
  vp10_xform_quant(x, plane, block, blk_row, blk_col, plane_bsize,
                   tx_size, VP10_XFORM_QUANT_B);
#endif  // CONFIG_NEW_QUANT

  if (p->eobs[block] > 0) {
#if CONFIG_VP9_HIGHBITDEPTH
    if (xd->cur_buf->flags & YV12_FLAG_HIGHBITDEPTH) {
      if (xd->lossless[xd->mi[0]->mbmi.segment_id]) {
        vp10_highbd_iwht4x4_add(dqcoeff, dst, pd->dst.stride,
                                p->eobs[block], xd->bd);
      } else {
        vp10_highbd_idct4x4_add(dqcoeff, dst, pd->dst.stride,
                                p->eobs[block], xd->bd);
      }
      return;
    }
#endif  // CONFIG_VP9_HIGHBITDEPTH
    if (xd->lossless[xd->mi[0]->mbmi.segment_id]) {
      vp10_iwht4x4_add(dqcoeff, dst, pd->dst.stride, p->eobs[block]);
    } else {
      vp10_idct4x4_add(dqcoeff, dst, pd->dst.stride, p->eobs[block]);
    }
  }
}

void vp10_encode_sby_pass1(MACROBLOCK *x, BLOCK_SIZE bsize) {
  vp10_subtract_plane(x, bsize, 0);
  vp10_foreach_transformed_block_in_plane(&x->e_mbd, bsize, 0,
                                          encode_block_pass1, x);
}

void vp10_encode_sb(MACROBLOCK *x, BLOCK_SIZE bsize) {
  MACROBLOCKD *const xd = &x->e_mbd;
  struct optimize_ctx ctx;
  MB_MODE_INFO *mbmi = &xd->mi[0]->mbmi;
  struct encode_b_args arg = {x, &ctx, &mbmi->skip, NULL, NULL, 1};
  int plane;

  mbmi->skip = 1;

  if (x->skip)
    return;

  for (plane = 0; plane < MAX_MB_PLANE; ++plane) {
#if CONFIG_VAR_TX
    // TODO(jingning): Clean this up.
    const struct macroblockd_plane *const pd = &xd->plane[plane];
    const BLOCK_SIZE plane_bsize = get_plane_block_size(bsize, pd);
    const int mi_width = num_4x4_blocks_wide_lookup[plane_bsize];
    const int mi_height = num_4x4_blocks_high_lookup[plane_bsize];
    const TX_SIZE max_tx_size = max_txsize_lookup[plane_bsize];
    const BLOCK_SIZE txb_size = txsize_to_bsize[max_tx_size];
    const int bh = num_4x4_blocks_wide_lookup[txb_size];
    int idx, idy;
    int block = 0;
    int step = num_4x4_blocks_txsize_lookup[max_tx_size];
    vp10_get_entropy_contexts(bsize, TX_4X4, pd, ctx.ta[plane], ctx.tl[plane]);
#else
    const struct macroblockd_plane* const pd = &xd->plane[plane];
    const TX_SIZE tx_size = plane ? get_uv_tx_size(mbmi, pd) : mbmi->tx_size;
    vp10_get_entropy_contexts(bsize, tx_size, pd, ctx.ta[plane], ctx.tl[plane]);
#endif
    vp10_subtract_plane(x, bsize, plane);
    arg.ta = ctx.ta[plane];
    arg.tl = ctx.tl[plane];

#if CONFIG_VAR_TX
    for (idy = 0; idy < mi_height; idy += bh) {
      for (idx = 0; idx < mi_width; idx += bh) {
        encode_block_inter(plane, block, idy, idx, plane_bsize,
                           max_tx_size, &arg);
        block += step;
      }
    }
#else
    vp10_foreach_transformed_block_in_plane(xd, bsize, plane, encode_block,
                                            &arg);
#endif
  }
}

#if CONFIG_SUPERTX
void vp10_encode_sb_supertx(MACROBLOCK *x, BLOCK_SIZE bsize) {
  MACROBLOCKD *const xd = &x->e_mbd;
  struct optimize_ctx ctx;
  MB_MODE_INFO *mbmi = &xd->mi[0]->mbmi;
  struct encode_b_args arg = {x, &ctx, &mbmi->skip, NULL, NULL, 1};
  int plane;

  mbmi->skip = 1;
  if (x->skip)
    return;

  for (plane = 0; plane < MAX_MB_PLANE; ++plane) {
    const struct macroblockd_plane* const pd = &xd->plane[plane];
#if CONFIG_VAR_TX
    const TX_SIZE tx_size = TX_4X4;
#else
    const TX_SIZE tx_size = plane ? get_uv_tx_size(mbmi, pd) : mbmi->tx_size;
#endif
    vp10_subtract_plane(x, bsize, plane);
    vp10_get_entropy_contexts(bsize, tx_size, pd,
                              ctx.ta[plane], ctx.tl[plane]);
    arg.ta = ctx.ta[plane];
    arg.tl = ctx.tl[plane];
    vp10_foreach_transformed_block_in_plane(xd, bsize, plane, encode_block,
                                            &arg);
  }
}
#endif  // CONFIG_SUPERTX

void vp10_encode_block_intra(int plane, int block, int blk_row, int blk_col,
                             BLOCK_SIZE plane_bsize,
                             TX_SIZE tx_size, void *arg) {
  struct encode_b_args* const args = arg;
  MACROBLOCK *const x = args->x;
  MACROBLOCKD *const xd = &x->e_mbd;
  MB_MODE_INFO *mbmi = &xd->mi[0]->mbmi;
  struct macroblock_plane *const p = &x->plane[plane];
  struct macroblockd_plane *const pd = &xd->plane[plane];
  tran_low_t *dqcoeff = BLOCK_OFFSET(pd->dqcoeff, block);
  PLANE_TYPE plane_type = (plane == 0) ? PLANE_TYPE_Y : PLANE_TYPE_UV;
  const TX_TYPE tx_type = get_tx_type(plane_type, xd, block, tx_size);
  PREDICTION_MODE mode;
  const int bwl = b_width_log2_lookup[plane_bsize];
  const int bhl = b_height_log2_lookup[plane_bsize];
  const int diff_stride = 4 * (1 << bwl);
  uint8_t *src, *dst;
  int16_t *src_diff;
  uint16_t *eob = &p->eobs[block];
  const int src_stride = p->src.stride;
  const int dst_stride = pd->dst.stride;
  const int tx1d_width = num_4x4_blocks_wide_txsize_lookup[tx_size] << 2;
  const int tx1d_height = num_4x4_blocks_high_txsize_lookup[tx_size] << 2;
  ENTROPY_CONTEXT *a = NULL, *l = NULL;
  int ctx;

  INV_TXFM_PARAM inv_txfm_param;

  assert(tx1d_width == tx1d_height);

  dst = &pd->dst.buf[4 * (blk_row * dst_stride + blk_col)];
  src = &p->src.buf[4 * (blk_row * src_stride + blk_col)];
  src_diff = &p->src_diff[4 * (blk_row * diff_stride + blk_col)];

  mode = plane == 0 ? get_y_mode(xd->mi[0], block) : mbmi->uv_mode;
  vp10_predict_intra_block(xd, bwl, bhl, tx_size, mode, dst, dst_stride, dst,
                           dst_stride, blk_col, blk_row, plane);
#if CONFIG_VP9_HIGHBITDEPTH
  if (xd->cur_buf->flags & YV12_FLAG_HIGHBITDEPTH) {
    vpx_highbd_subtract_block(tx1d_height, tx1d_width, src_diff, diff_stride,
                              src, src_stride, dst, dst_stride, xd->bd);
  } else {
    vpx_subtract_block(tx1d_height, tx1d_width, src_diff, diff_stride, src,
                       src_stride, dst, dst_stride);
  }
#else
  vpx_subtract_block(tx1d_height, tx1d_width, src_diff, diff_stride, src,
                     src_stride, dst, dst_stride);
#endif  // CONFIG_VP9_HIGHBITDEPTH

  a = &args->ta[blk_col];
  l = &args->tl[blk_row];
  ctx = combine_entropy_contexts(*a, *l);

  if (args->enable_optimize_b) {
#if CONFIG_NEW_QUANT
    vp10_xform_quant_fp_nuq(x, plane, block, blk_row, blk_col, plane_bsize,
                            tx_size, ctx);
#else  // CONFIG_NEW_QUANT
    vp10_xform_quant(x, plane, block, blk_row, blk_col, plane_bsize, tx_size,
                     VP10_XFORM_QUANT_FP);
#endif  // CONFIG_NEW_QUANT
    if (p->eobs[block]) {
      *a = *l = vp10_optimize_b(x, plane, block, tx_size, ctx) > 0;
    } else {
      *a = *l = 0;
    }
  } else {
    vp10_xform_quant(x, plane, block, blk_row, blk_col, plane_bsize, tx_size,
                     VP10_XFORM_QUANT_B);
    *a = *l = p->eobs[block] > 0;
  }

  if (*eob) {
    // inverse transform
    inv_txfm_param.tx_type = tx_type;
    inv_txfm_param.tx_size = tx_size;
    inv_txfm_param.eob = *eob;
    inv_txfm_param.lossless = xd->lossless[mbmi->segment_id];
#if CONFIG_VP9_HIGHBITDEPTH
    inv_txfm_param.bd = xd->bd;
    if (xd->cur_buf->flags & YV12_FLAG_HIGHBITDEPTH) {
      highbd_inv_txfm_add(dqcoeff, dst, dst_stride, &inv_txfm_param);
    } else {
      inv_txfm_add(dqcoeff, dst, dst_stride, &inv_txfm_param);
    }
#else
    inv_txfm_add(dqcoeff, dst, dst_stride, &inv_txfm_param);
#endif  // CONFIG_VP9_HIGHBITDEPTH

    *(args->skip) = 0;
  }
}

void vp10_encode_intra_block_plane(MACROBLOCK *x, BLOCK_SIZE bsize, int plane,
                                   int enable_optimize_b) {
  const MACROBLOCKD *const xd = &x->e_mbd;
  ENTROPY_CONTEXT ta[2 * MAX_MIB_SIZE];
  ENTROPY_CONTEXT tl[2 * MAX_MIB_SIZE];

  struct encode_b_args arg = {x, NULL, &xd->mi[0]->mbmi.skip,
                              ta, tl, enable_optimize_b};
  if (enable_optimize_b) {
    const struct macroblockd_plane* const pd = &xd->plane[plane];
    const TX_SIZE tx_size = plane ? get_uv_tx_size(&xd->mi[0]->mbmi, pd) :
        xd->mi[0]->mbmi.tx_size;
    vp10_get_entropy_contexts(bsize, tx_size, pd, ta, tl);
  }
  vp10_foreach_transformed_block_in_plane(xd, bsize, plane,
                                          vp10_encode_block_intra, &arg);
}
