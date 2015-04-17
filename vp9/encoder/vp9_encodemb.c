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

struct encode_b_args {
  MACROBLOCK *x;
  struct optimize_ctx *ctx;
  int8_t *skip;
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
  signed char   token;
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
  vp9_token_state tokens[MAX_NUM_COEFS + 1][2];
  unsigned best_index[MAX_NUM_COEFS + 1][2];
  uint8_t token_cache[MAX_NUM_COEFS];
  const tran_low_t *const coeff = BLOCK_OFFSET(mb->plane[plane].coeff, block);
  tran_low_t *const qcoeff = BLOCK_OFFSET(p->qcoeff, block);
  tran_low_t *const dqcoeff = BLOCK_OFFSET(pd->dqcoeff, block);
  const int eob = p->eobs[block];
  const PLANE_TYPE type = pd->plane_type;
  const int default_eob = 16 << (tx_size << 1);
  const int shift = (tx_size >= TX_32X32 ? tx_size - TX_16X16 : 0);
  const int mul = 1 << shift;
#if CONFIG_TX_SKIP
  const int16_t *dequant_ptr =
      xd->mi[0].src_mi->mbmi.tx_skip[plane != 0] ?
          pd->dequant_pxd : pd->dequant;
#else
  const int16_t *dequant_ptr = pd->dequant;
#endif  // CONFIG_TX_SKIP
#if CONFIG_NEW_QUANT
#if CONFIG_TX_SKIP
  const int use_rect_quant = is_rect_quant_used(&xd->mi[0].src_mi->mbmi, plane);
#endif  // CONFIG_TX_SKIP
  const dequant_val_type_nuq *dequant_val = pd->dequant_val_nuq;
#endif  // CONFIG_NEW_QUANT
#if CONFIG_TX_SKIP
  const uint8_t *const band_translate =
      xd->mi[0].src_mi->mbmi.tx_skip[plane != 0] ?
      vp9_coefband_tx_skip : get_band_translate(tx_size);
#else
  const uint8_t *const band_translate = get_band_translate(tx_size);
#endif  // CONFIG_TX_SKIP
  const scan_order *const so = get_scan(xd, tx_size, type, block);
  const int16_t *const scan = so->scan;
  const int16_t *const nb = so->neighbors;
  int next = eob, sz = 0;
  int64_t rdmult = mb->rdmult * plane_rd_mult[type], rddiv = mb->rddiv;
  int64_t rd_cost0, rd_cost1;
  int rate0, rate1, error0, error1, t0, t1;
  int best, band, pt, i, final_eob;
  const TOKENVALUE *dct_value_tokens;
  const int16_t *dct_value_cost;

  assert((!type && !plane) || (type && plane));
  assert(eob <= default_eob);

#if CONFIG_TX_SKIP && CONFIG_NEW_QUANT
  if (xd->mi[0].src_mi->mbmi.tx_skip[plane != 0])
    dequant_val = pd->dequant_val_nuq_pxd;
#endif  // CONFIG_TX_SKIP && CONFIG_NEW_QUANT

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

#if CONFIG_VP9_HIGHBITDEPTH
  if (xd->bd == 12) {
    dct_value_tokens = vp9_dct_value_tokens_high12_ptr;
    dct_value_cost = vp9_dct_value_cost_high12_ptr;
  } else if (xd->bd == 10) {
    dct_value_tokens = vp9_dct_value_tokens_high10_ptr;
    dct_value_cost = vp9_dct_value_cost_high10_ptr;
  } else {
    dct_value_tokens = vp9_dct_value_tokens_ptr;
    dct_value_cost = vp9_dct_value_cost_ptr;
  }
#else
  dct_value_tokens = vp9_dct_value_tokens_ptr;
  dct_value_cost = vp9_dct_value_cost_ptr;
#endif
  for (i = 0; i < eob; i++)
    token_cache[scan[i]] =
        vp9_pt_energy_class[dct_value_tokens[qcoeff[scan[i]]].token];

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
      t0 = (dct_value_tokens + x)->token;
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
      base_bits = dct_value_cost[x];
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

#if CONFIG_NEW_QUANT
#if CONFIG_TX_SKIP
      if (use_rect_quant) {
        shortcut =
            ((abs(x) * dequant_ptr[rc != 0] > abs(coeff[rc]) * mul) &&
             ((abs(x) - 1) * dequant_ptr[rc != 0] < abs(coeff[rc]) * mul));
      } else {
        shortcut = (
            (vp9_dequant_abscoeff_nuq(
                abs(x), dequant_ptr[rc != 0],
                dequant_val[band_translate[i]]) > abs(coeff[rc]) * mul) &&
            (vp9_dequant_abscoeff_nuq(
                abs(x) - 1, dequant_ptr[rc != 0],
                dequant_val[band_translate[i]]) < abs(coeff[rc]) * mul));
      }
#else   // CONFIG_TX_SKIP
      shortcut = (
          (vp9_dequant_abscoeff_nuq(
              abs(x), dequant_ptr[rc != 0],
              dequant_val[band_translate[i]]) > abs(coeff[rc]) * mul) &&
          (vp9_dequant_abscoeff_nuq(
              abs(x) - 1, dequant_ptr[rc != 0],
              dequant_val[band_translate[i]]) < abs(coeff[rc]) * mul));
#endif  // CONFIG_TX_SKIP
#else   // CONFIG_NEW_QUANT
      shortcut = ((abs(x) * dequant_ptr[rc != 0] > abs(coeff[rc]) * mul) &&
                  ((abs(x) - 1) * dequant_ptr[rc != 0] < abs(coeff[rc]) * mul));
#endif  // CONFIG_NEW_QUANT

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
      } else {
        t0 = t1 = (dct_value_tokens + x)->token;
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
      base_bits = dct_value_cost[x];

      if (shortcut) {
#if CONFIG_NEW_QUANT
#if CONFIG_TX_SKIP
        if (use_rect_quant) {
#if CONFIG_VP9_HIGHBITDEPTH
          if (xd->cur_buf->flags & YV12_FLAG_HIGHBITDEPTH) {
            dx -= ((dequant_ptr[rc != 0] >> (xd->bd - 8)) + sz) ^ sz;
          } else {
            dx -= (dequant_ptr[rc != 0] + sz) ^ sz;
          }
#else
          dx -= (dequant_ptr[rc != 0] + sz) ^ sz;
#endif  // CONFIG_VP9_HIGHBITDEPTH
        } else {
          dx = vp9_dequant_coeff_nuq(
              x, dequant_ptr[rc != 0],
              dequant_val[band_translate[i]]) - coeff[rc] * mul;
#if CONFIG_VP9_HIGHBITDEPTH
          if (xd->cur_buf->flags & YV12_FLAG_HIGHBITDEPTH) {
            dx >>= xd->bd - 8;
          }
#endif  // CONFIG_VP9_HIGHBITDEPTH
        }
#else   // CONFIG_TX_SKIP
        dx = vp9_dequant_coeff_nuq(
            x, dequant_ptr[rc != 0],
            dequant_val[band_translate[i]]) - coeff[rc] * mul;
#if CONFIG_VP9_HIGHBITDEPTH
        if (xd->cur_buf->flags & YV12_FLAG_HIGHBITDEPTH) {
          dx >>= xd->bd - 8;
        }
#endif  // CONFIG_VP9_HIGHBITDEPTH
#endif  // CONFIG_TX_SKIP
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
#if CONFIG_NEW_QUANT
#if CONFIG_TX_SKIP
    if (use_rect_quant) {
      dqcoeff[rc] = (x * dequant_ptr[rc != 0]) / mul;
    } else {
      dqcoeff[rc] = vp9_dequant_abscoeff_nuq(abs(x), dequant_ptr[rc != 0],
                                             dequant_val[band_translate[i]]);
      if (shift) dqcoeff[rc] = ROUND_POWER_OF_TWO(dqcoeff[rc], shift);
      if (x < 0) dqcoeff[rc] = -dqcoeff[rc];
    }
#else   // CONFIG_TX_SKIP
    dqcoeff[rc] = vp9_dequant_abscoeff_nuq(abs(x), dequant_ptr[rc != 0],
                                           dequant_val[band_translate[i]]);
    if (shift) dqcoeff[rc] = ROUND_POWER_OF_TWO(dqcoeff[rc], shift);
    if (x < 0) dqcoeff[rc] = -dqcoeff[rc];
#endif  // CONFIG_TX_SKIP
#else
    dqcoeff[rc] = (x * dequant_ptr[rc != 0]) / mul;
#endif  // CONFIG_NEW_QUANT
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

#if CONFIG_EXT_TX
static void copy_block(const int16_t *src, int src_stride, int l,
                       int16_t *dest, int dest_stride) {
  int i;
  for (i = 0; i < l; ++i) {
    vpx_memcpy(dest + dest_stride * i, src + src_stride * i,
               l * sizeof(int16_t));
  }
}

static void fliplr(int16_t *dest, int stride, int l) {
  int i, j;
  for (i = 0; i < l; ++i) {
    for (j = 0; j < l / 2; ++j) {
      const int16_t tmp = dest[i * stride + j];
      dest[i * stride + j] = dest[i * stride + l - 1 - j];
      dest[i * stride + l - 1 - j] = tmp;
    }
  }
}

static void flipud(int16_t *dest, int stride, int l) {
  int i, j;
  for (j = 0; j < l; ++j) {
    for (i = 0; i < l / 2; ++i) {
      const int16_t tmp = dest[i * stride + j];
      dest[i * stride + j] = dest[(l - 1 - i) * stride + j];
      dest[(l - 1 - i) * stride + j] = tmp;
    }
  }
}

static void fliplrud(int16_t *dest, int stride, int l) {
  int i, j;
  for (i = 0; i < l / 2; ++i) {
    for (j = 0; j < l; ++j) {
      const int16_t tmp = dest[i * stride + j];
      dest[i * stride + j] = dest[(l - 1 - i) * stride + l - 1 - j];
      dest[(l - 1 - i) * stride + l - 1 - j] = tmp;
    }
  }
}

static void copy_fliplr(const int16_t *src, int src_stride, int l,
                          int16_t *dest, int dest_stride) {
  copy_block(src, src_stride, l, dest, dest_stride);
  fliplr(dest, dest_stride, l);
}

static void copy_flipud(const int16_t *src, int src_stride, int l,
                          int16_t *dest, int dest_stride) {
  copy_block(src, src_stride, l, dest, dest_stride);
  flipud(dest, dest_stride, l);
}

static void copy_fliplrud(const int16_t *src, int src_stride, int l,
                            int16_t *dest, int dest_stride) {
  copy_block(src, src_stride, l, dest, dest_stride);
  fliplrud(dest, dest_stride, l);
}

static void forw_tx16x16(MACROBLOCK *x, int plane,
                         const int16_t *src_diff, int diff_stride,
                         tran_low_t *const coeff) {
  MACROBLOCKD *const xd = &x->e_mbd;
  int16_t src_diff2[256];
  TX_TYPE tx_type = get_tx_type(plane, xd);
  if (tx_type == DCT_DCT) {
    vp9_fdct16x16(src_diff, coeff, diff_stride);
  } else if (tx_type == FLIPADST_DCT) {
    copy_flipud(src_diff, diff_stride, 16, src_diff2, 16);
    vp9_fht16x16(src_diff2, coeff, 16, ADST_DCT);
  } else if (tx_type == DCT_FLIPADST) {
    copy_fliplr(src_diff, diff_stride, 16, src_diff2, 16);
    vp9_fht16x16(src_diff2, coeff, 16, DCT_ADST);
  } else if (tx_type == FLIPADST_FLIPADST) {
    copy_fliplrud(src_diff, diff_stride, 16, src_diff2, 16);
    vp9_fht16x16(src_diff2, coeff, 16, ADST_ADST);
  } else if (tx_type == ADST_FLIPADST) {
    copy_fliplr(src_diff, diff_stride, 16, src_diff2, 16);
    vp9_fht16x16(src_diff2, coeff, 16, ADST_ADST);
  } else if (tx_type == FLIPADST_ADST) {
    copy_flipud(src_diff, diff_stride, 16, src_diff2, 16);
    vp9_fht16x16(src_diff2, coeff, 16, ADST_ADST);
  } else {
    vp9_fht16x16(src_diff, coeff, diff_stride, tx_type);
  }
}

static void forw_tx8x8(MACROBLOCK *x, int plane,
                       const int16_t *src_diff, int diff_stride,
                       tran_low_t *const coeff) {
  MACROBLOCKD *const xd = &x->e_mbd;
  int16_t src_diff2[64];
  TX_TYPE tx_type = get_tx_type(plane, xd);
  if (tx_type == DCT_DCT) {
    vp9_fdct8x8(src_diff, coeff, diff_stride);
  } else if (tx_type == FLIPADST_DCT) {
    copy_flipud(src_diff, diff_stride, 8, src_diff2, 8);
    vp9_fht8x8(src_diff2, coeff, 8, ADST_DCT);
  } else if (tx_type == DCT_FLIPADST) {
    copy_fliplr(src_diff, diff_stride, 8, src_diff2, 8);
    vp9_fht8x8(src_diff2, coeff, 8, DCT_ADST);
  } else if (tx_type == FLIPADST_FLIPADST) {
    copy_fliplrud(src_diff, diff_stride, 8, src_diff2, 8);
    vp9_fht8x8(src_diff2, coeff, 8, ADST_ADST);
  } else if (tx_type == ADST_FLIPADST) {
    copy_fliplr(src_diff, diff_stride, 8, src_diff2, 8);
    vp9_fht8x8(src_diff2, coeff, 8, ADST_ADST);
  } else if (tx_type == FLIPADST_ADST) {
    copy_flipud(src_diff, diff_stride, 8, src_diff2, 8);
    vp9_fht8x8(src_diff2, coeff, 8, ADST_ADST);
  } else {
    vp9_fht8x8(src_diff, coeff, diff_stride, tx_type);
  }
}

static void forw_tx4x4(MACROBLOCK *x, int plane, int block,
                       const int16_t *src_diff, int diff_stride,
                       tran_low_t *const coeff) {
  MACROBLOCKD *const xd = &x->e_mbd;
  int16_t src_diff2[16];
  TX_TYPE tx_type = get_tx_type_4x4(plane, xd, block);
  if (tx_type == DCT_DCT) {
    x->fwd_txm4x4(src_diff, coeff, diff_stride);
  } else if (tx_type == FLIPADST_DCT) {
    copy_flipud(src_diff, diff_stride, 4, src_diff2, 4);
    vp9_fht4x4(src_diff2, coeff, 4, ADST_DCT);
  } else if (tx_type == DCT_FLIPADST) {
    copy_fliplr(src_diff, diff_stride, 4, src_diff2, 4);
    vp9_fht4x4(src_diff2, coeff, 4, DCT_ADST);
  } else if (tx_type == FLIPADST_FLIPADST) {
    copy_fliplrud(src_diff, diff_stride, 4, src_diff2, 4);
    vp9_fht4x4(src_diff2, coeff, 4, ADST_ADST);
  } else if (tx_type == ADST_FLIPADST) {
    copy_fliplr(src_diff, diff_stride, 4, src_diff2, 4);
    vp9_fht4x4(src_diff2, coeff, 4, ADST_ADST);
  } else if (tx_type == FLIPADST_ADST) {
    copy_flipud(src_diff, diff_stride, 4, src_diff2, 4);
    vp9_fht4x4(src_diff2, coeff, 4, ADST_ADST);
  } else {
    vp9_fht4x4(src_diff, coeff, diff_stride, tx_type);
  }
}

#if CONFIG_VP9_HIGHBITDEPTH
static void highbd_forw_tx16x16(MACROBLOCK *x, int plane,
                                const int16_t *src_diff, int diff_stride,
                                tran_low_t *const coeff) {
  MACROBLOCKD *const xd = &x->e_mbd;
  int16_t src_diff2[256];
  TX_TYPE tx_type = get_tx_type(plane, xd);
  if (tx_type == DCT_DCT) {
    vp9_highbd_fdct16x16(src_diff, coeff, diff_stride);
  } else if (tx_type == FLIPADST_DCT) {
    copy_flipud(src_diff, diff_stride, 16, src_diff2, 16);
    vp9_highbd_fht16x16(src_diff2, coeff, 16, ADST_DCT);
  } else if (tx_type == DCT_FLIPADST) {
    copy_fliplr(src_diff, diff_stride, 16, src_diff2, 16);
    vp9_highbd_fht16x16(src_diff2, coeff, 16, DCT_ADST);
  } else if (tx_type == FLIPADST_FLIPADST) {
    copy_fliplrud(src_diff, diff_stride, 16, src_diff2, 16);
    vp9_highbd_fht16x16(src_diff2, coeff, 16, ADST_ADST);
  } else if (tx_type == ADST_FLIPADST) {
    copy_fliplr(src_diff, diff_stride, 16, src_diff2, 16);
    vp9_highbd_fht16x16(src_diff2, coeff, 16, ADST_ADST);
  } else if (tx_type == FLIPADST_ADST) {
    copy_flipud(src_diff, diff_stride, 16, src_diff2, 16);
    vp9_highbd_fht16x16(src_diff2, coeff, 16, ADST_ADST);
  } else {
    vp9_highbd_fht16x16(src_diff, coeff, diff_stride, tx_type);
  }
}

static void highbd_forw_tx8x8(MACROBLOCK *x, int plane,
                              const int16_t *src_diff, int diff_stride,
                              tran_low_t *const coeff) {
  MACROBLOCKD *const xd = &x->e_mbd;
  int16_t src_diff2[64];
  TX_TYPE tx_type = get_tx_type(plane, xd);
  if (tx_type == DCT_DCT) {
    vp9_highbd_fdct8x8(src_diff, coeff, diff_stride);
  } else if (tx_type == FLIPADST_DCT) {
    copy_flipud(src_diff, diff_stride, 8, src_diff2, 8);
    vp9_highbd_fht8x8(src_diff2, coeff, 8, ADST_DCT);
  } else if (tx_type == DCT_FLIPADST) {
    copy_fliplr(src_diff, diff_stride, 8, src_diff2, 8);
    vp9_highbd_fht8x8(src_diff2, coeff, 8, DCT_ADST);
  } else if (tx_type == FLIPADST_FLIPADST) {
    copy_fliplrud(src_diff, diff_stride, 8, src_diff2, 8);
    vp9_highbd_fht8x8(src_diff2, coeff, 8, ADST_ADST);
  } else if (tx_type == ADST_FLIPADST) {
    copy_fliplr(src_diff, diff_stride, 8, src_diff2, 8);
    vp9_highbd_fht8x8(src_diff2, coeff, 8, ADST_ADST);
  } else if (tx_type == FLIPADST_ADST) {
    copy_flipud(src_diff, diff_stride, 8, src_diff2, 8);
    vp9_highbd_fht8x8(src_diff2, coeff, 8, ADST_ADST);
  } else {
    vp9_highbd_fht8x8(src_diff, coeff, diff_stride, tx_type);
  }
}

static void highbd_forw_tx4x4(MACROBLOCK *x, int plane, int block,
                              const int16_t *src_diff, int diff_stride,
                              tran_low_t *const coeff) {
  MACROBLOCKD *const xd = &x->e_mbd;
  int16_t src_diff2[16];
  TX_TYPE tx_type = get_tx_type_4x4(plane, xd, block);
  if (tx_type == DCT_DCT) {
    x->fwd_txm4x4(src_diff, coeff, diff_stride);
  } else if (tx_type == FLIPADST_DCT) {
    copy_flipud(src_diff, diff_stride, 4, src_diff2, 4);
    vp9_highbd_fht4x4(src_diff2, coeff, 4, ADST_DCT);
  } else if (tx_type == DCT_FLIPADST) {
    copy_fliplr(src_diff, diff_stride, 4, src_diff2, 4);
    vp9_highbd_fht4x4(src_diff2, coeff, 4, DCT_ADST);
  } else if (tx_type == FLIPADST_FLIPADST) {
    copy_fliplrud(src_diff, diff_stride, 4, src_diff2, 4);
    vp9_highbd_fht4x4(src_diff2, coeff, 4, ADST_ADST);
  } else if (tx_type == ADST_FLIPADST) {
    copy_fliplr(src_diff, diff_stride, 4, src_diff2, 4);
    vp9_highbd_fht4x4(src_diff2, coeff, 4, ADST_ADST);
  } else if (tx_type == FLIPADST_ADST) {
    copy_flipud(src_diff, diff_stride, 4, src_diff2, 4);
    vp9_highbd_fht4x4(src_diff2, coeff, 4, ADST_ADST);
  } else {
    vp9_highbd_fht4x4(src_diff, coeff, diff_stride, tx_type);
  }
}
#endif  // CONFIG_VP9_HIGHBITDEPTH
#endif  // CONFIG_EXT_TX

#if CONFIG_NEW_QUANT
void vp9_xform_quant_nuq(MACROBLOCK *x, int plane, int block,
                         BLOCK_SIZE plane_bsize, TX_SIZE tx_size) {
  MACROBLOCKD *const xd = &x->e_mbd;
  const struct macroblock_plane *const p = &x->plane[plane];
  const struct macroblockd_plane *const pd = &xd->plane[plane];
#if CONFIG_TX_SKIP
  MB_MODE_INFO *mbmi = &xd->mi[0].src_mi->mbmi;
  int shift = mbmi->tx_skip_shift;
  const scan_order *const scan_order =
      mbmi->tx_skip[plane != 0] ? &vp9_default_scan_orders_pxd[tx_size] :
          &vp9_default_scan_orders[tx_size];
#else
  const scan_order *const scan_order = &vp9_default_scan_orders[tx_size];
#endif  // CONFIG_TX_SKIP
  tran_low_t *const coeff = BLOCK_OFFSET(p->coeff, block);
  tran_low_t *const qcoeff = BLOCK_OFFSET(p->qcoeff, block);
  tran_low_t *const dqcoeff = BLOCK_OFFSET(pd->dqcoeff, block);
  uint16_t *const eob = &p->eobs[block];
  const int diff_stride = 4 * num_4x4_blocks_wide_lookup[plane_bsize];
  int i, j;
  const int16_t *src_diff;
  const uint8_t* band = get_band_translate(tx_size);

  txfrm_block_to_raster_xy(plane_bsize, tx_size, block, &i, &j);
  src_diff = &p->src_diff[4 * (j * diff_stride + i)];

#if CONFIG_TX_SKIP
  if (mbmi->tx_skip[plane != 0]) {
    int bs = 4 << tx_size;
    band = vp9_coefband_tx_skip;
    vp9_tx_identity(src_diff, coeff, diff_stride, bs, shift);
    if (tx_size <= TX_16X16) {
#if CONFIG_VP9_HIGHBITDEPTH
      if (xd->cur_buf->flags & YV12_FLAG_HIGHBITDEPTH)
        vp9_highbd_quantize_nuq(coeff, bs * bs, x->skip_block,
                                p->quant_pxd, p->quant_shift_pxd,
                                pd->dequant_pxd,
                                (const cumbins_type_nuq *)p->cumbins_nuq_pxd,
                                (const dequant_val_type_nuq *)
                                pd->dequant_val_nuq_pxd,
                                qcoeff, dqcoeff, eob,
                                scan_order->scan, band);
      else
#endif  // CONFIG_VP9_HIGHBITDEPTH
        vp9_quantize_nuq(coeff, bs * bs, x->skip_block,
                         p->quant_pxd, p->quant_shift_pxd, pd->dequant_pxd,
                         (const cumbins_type_nuq *)p->cumbins_nuq_pxd,
                         (const dequant_val_type_nuq *)pd->dequant_val_nuq_pxd,
                         qcoeff, dqcoeff, eob,
                         scan_order->scan, band);
    } else if (tx_size == TX_32X32) {
#if CONFIG_VP9_HIGHBITDEPTH
      if (xd->cur_buf->flags & YV12_FLAG_HIGHBITDEPTH)
        vp9_highbd_quantize_32x32_nuq(coeff, bs * bs, x->skip_block,
                                      p->quant_pxd, p->quant_shift_pxd,
                                      pd->dequant_pxd,
                                      (const cumbins_type_nuq *)
                                      p->cumbins_nuq_pxd,
                                      (const dequant_val_type_nuq *)
                                      pd->dequant_val_nuq,
                                      qcoeff, dqcoeff, eob,
                                      scan_order->scan, band);
      else
#endif  // CONFIG_VP9_HIGHBITDEPTH
        vp9_quantize_32x32_nuq(coeff, bs * bs, x->skip_block,
                               p->quant_pxd, p->quant_shift_pxd,
                               pd->dequant_pxd,
                               (const cumbins_type_nuq *)p->cumbins_nuq_pxd,
                               (const dequant_val_type_nuq *)
                               pd->dequant_val_nuq,
                               qcoeff, dqcoeff, eob,
                               scan_order->scan, band);
    }
#if CONFIG_TX64X64
    else if (tx_size == TX_64X64) {
#if CONFIG_VP9_HIGHBITDEPTH
      if (xd->cur_buf->flags & YV12_FLAG_HIGHBITDEPTH)
        vp9_highbd_quantize_64x64_nuq(coeff, bs * bs, x->skip_block,
                                      p->quant_pxd, p->quant_shift_pxd,
                                      pd->dequant_pxd,
                                      (const cumbins_type_nuq *)
                                      p->cumbins_nuq_pxd,
                                      (const dequant_val_type_nuq *)
                                      pd->dequant_val_nuq,
                                      qcoeff, dqcoeff, eob,
                                      scan_order->scan, band);
      else
#endif  // CONFIG_VP9_HIGHBITDEPTH
        vp9_quantize_64x64_nuq(coeff, bs * bs, x->skip_block,
                               p->quant_pxd, p->quant_shift_pxd,
                               pd->dequant_pxd,
                               (const cumbins_type_nuq *)p->cumbins_nuq_pxd,
                               (const dequant_val_type_nuq *)
                               pd->dequant_val_nuq,
                               qcoeff, dqcoeff, eob,
                               scan_order->scan, band);
    }
#endif  // CONFIG_TX64X64

    return;
  }
#endif  // CONFIG_TX_SKIP

#if CONFIG_VP9_HIGHBITDEPTH
  if (xd->cur_buf->flags & YV12_FLAG_HIGHBITDEPTH) {
    switch (tx_size) {
#if CONFIG_TX64X64
      case TX_64X64:
        vp9_highbd_fdct64x64(src_diff, coeff, diff_stride);
        vp9_highbd_quantize_64x64_nuq(coeff, 4096, x->skip_block,
                                      p->quant, p->quant_shift, pd->dequant,
                                      (const cumbins_type_nuq *)p->cumbins_nuq,
                                      (const dequant_val_type_nuq *)
                                          pd->dequant_val_nuq,
                                      qcoeff, dqcoeff, eob,
                                      scan_order->scan, band);
        break;
#endif  // CONFIG_TX64X64
      case TX_32X32:
        highbd_fdct32x32(x->use_lp32x32fdct, src_diff, coeff, diff_stride);
        vp9_highbd_quantize_32x32_nuq(coeff, 1024, x->skip_block,
                                      p->quant, p->quant_shift, pd->dequant,
                                      (const cumbins_type_nuq *)p->cumbins_nuq,
                                      (const dequant_val_type_nuq *)
                                          pd->dequant_val_nuq,
                                      qcoeff, dqcoeff, eob,
                                      scan_order->scan, band);
        break;
      case TX_16X16:
#if CONFIG_EXT_TX
        highbd_forw_tx16x16(x, plane, src_diff, diff_stride, coeff);
#else
        vp9_highbd_fdct16x16(src_diff, coeff, diff_stride);
#endif
        vp9_highbd_quantize_nuq(coeff, 256, x->skip_block,
                                p->quant, p->quant_shift, pd->dequant,
                                (const cumbins_type_nuq *)p->cumbins_nuq,
                                (const dequant_val_type_nuq *)
                                    pd->dequant_val_nuq,
                                qcoeff, dqcoeff, eob,
                                scan_order->scan, band);
        break;
      case TX_8X8:
#if CONFIG_EXT_TX
        highbd_forw_tx8x8(x, plane, src_diff, diff_stride, coeff);
#else
        vp9_highbd_fdct8x8(src_diff, coeff, diff_stride);
#endif
        vp9_highbd_quantize_nuq(coeff, 64, x->skip_block,
                                p->quant, p->quant_shift, pd->dequant,
                                (const cumbins_type_nuq *)p->cumbins_nuq,
                                (const dequant_val_type_nuq *)
                                    pd->dequant_val_nuq,
                                qcoeff, dqcoeff, eob,
                                scan_order->scan, band);
        break;
      case TX_4X4:
#if CONFIG_EXT_TX
        highbd_forw_tx4x4(x, plane, block, src_diff, diff_stride, coeff);
#else
        x->fwd_txm4x4(src_diff, coeff, diff_stride);
#endif
        vp9_highbd_quantize_nuq(coeff, 16, x->skip_block,
                                p->quant, p->quant_shift, pd->dequant,
                                (const cumbins_type_nuq *)p->cumbins_nuq,
                                (const dequant_val_type_nuq *)
                                    pd->dequant_val_nuq,
                                qcoeff, dqcoeff, eob,
                                scan_order->scan, band);
        break;
      default:
        assert(0);
    }
    return;
  }
#endif  // CONFIG_VP9_HIGHBITDEPTH

  switch (tx_size) {
#if CONFIG_TX64X64
    case TX_64X64:
      vp9_fdct64x64(src_diff, coeff, diff_stride);
      vp9_quantize_64x64_nuq(coeff, 4096, x->skip_block,
                             p->quant, p->quant_shift, pd->dequant,
                             (const cumbins_type_nuq *)p->cumbins_nuq,
                             (const dequant_val_type_nuq *)pd->dequant_val_nuq,
                             qcoeff, dqcoeff, eob,
                             scan_order->scan, band);
      break;
#endif  // CONFIG_TX64X64
    case TX_32X32:
      fdct32x32(x->use_lp32x32fdct, src_diff, coeff, diff_stride);
      vp9_quantize_32x32_nuq(coeff, 1024, x->skip_block,
                             p->quant, p->quant_shift, pd->dequant,
                             (const cumbins_type_nuq *)p->cumbins_nuq,
                             (const dequant_val_type_nuq *)pd->dequant_val_nuq,
                             qcoeff, dqcoeff, eob,
                             scan_order->scan, band);
      break;
    case TX_16X16:
#if CONFIG_EXT_TX
      forw_tx16x16(x, plane, src_diff, diff_stride, coeff);
#else
      vp9_fdct16x16(src_diff, coeff, diff_stride);
#endif
      vp9_quantize_nuq(coeff, 256, x->skip_block,
                       p->quant, p->quant_shift, pd->dequant,
                       (const cumbins_type_nuq *)p->cumbins_nuq,
                       (const dequant_val_type_nuq *)pd->dequant_val_nuq,
                       qcoeff, dqcoeff, eob,
                       scan_order->scan, band);
      break;
    case TX_8X8:
#if CONFIG_EXT_TX
      forw_tx8x8(x, plane, src_diff, diff_stride, coeff);
#else
      vp9_fdct8x8(src_diff, coeff, diff_stride);
#endif
      vp9_quantize_nuq(coeff, 64, x->skip_block,
                       p->quant, p->quant_shift, pd->dequant,
                       (const cumbins_type_nuq *)p->cumbins_nuq,
                       (const dequant_val_type_nuq *)pd->dequant_val_nuq,
                       qcoeff, dqcoeff, eob,
                       scan_order->scan, band);
      break;
    case TX_4X4:
#if CONFIG_EXT_TX
      forw_tx4x4(x, plane, block, src_diff, diff_stride, coeff);
#else
      x->fwd_txm4x4(src_diff, coeff, diff_stride);
#endif
      vp9_quantize_nuq(coeff, 16, x->skip_block,
                       p->quant, p->quant_shift, pd->dequant,
                       (const cumbins_type_nuq *)p->cumbins_nuq,
                       (const dequant_val_type_nuq *)pd->dequant_val_nuq,
                       qcoeff, dqcoeff, eob,
                       scan_order->scan, band);
      break;
    default:
      assert(0);
      break;
  }
}

void vp9_xform_quant_fp_nuq(MACROBLOCK *x, int plane, int block,
                            BLOCK_SIZE plane_bsize, TX_SIZE tx_size) {
  MACROBLOCKD *const xd = &x->e_mbd;
  const struct macroblock_plane *const p = &x->plane[plane];
  const struct macroblockd_plane *const pd = &xd->plane[plane];
#if CONFIG_TX_SKIP
  MB_MODE_INFO *mbmi = &xd->mi[0].src_mi->mbmi;
  int shift = mbmi->tx_skip_shift;
  const scan_order *const scan_order =
      mbmi->tx_skip[plane != 0] ? &vp9_default_scan_orders_pxd[tx_size] :
          &vp9_default_scan_orders[tx_size];
#else
  const scan_order *const scan_order = &vp9_default_scan_orders[tx_size];
#endif  // CONFIG_TX_SKIP
  tran_low_t *const coeff = BLOCK_OFFSET(p->coeff, block);
  tran_low_t *const qcoeff = BLOCK_OFFSET(p->qcoeff, block);
  tran_low_t *const dqcoeff = BLOCK_OFFSET(pd->dqcoeff, block);
  uint16_t *const eob = &p->eobs[block];
  const int diff_stride = 4 * num_4x4_blocks_wide_lookup[plane_bsize];
  int i, j;
  const int16_t *src_diff;
  const uint8_t* band = get_band_translate(tx_size);

  txfrm_block_to_raster_xy(plane_bsize, tx_size, block, &i, &j);
  src_diff = &p->src_diff[4 * (j * diff_stride + i)];

#if CONFIG_TX_SKIP
  if (mbmi->tx_skip[plane != 0]) {
    int bs = 4 << tx_size;
    band = vp9_coefband_tx_skip;
    vp9_tx_identity(src_diff, coeff, diff_stride, bs, shift);
    if (tx_size <= TX_16X16) {
#if CONFIG_VP9_HIGHBITDEPTH
      if (xd->cur_buf->flags & YV12_FLAG_HIGHBITDEPTH)
        vp9_highbd_quantize_fp_nuq(coeff, bs * bs, x->skip_block,
                                   p->quant_pxd_fp, pd->dequant_pxd,
                                   (const cumbins_type_nuq *)p->cumbins_nuq_pxd,
                                   (const dequant_val_type_nuq *)
                                   pd->dequant_val_nuq_pxd,
                                   qcoeff, dqcoeff, eob,
                                   scan_order->scan, band);
      else
#endif  // CONFIG_VP9_HIGHBITDEPTH
        vp9_quantize_fp_nuq(coeff, bs * bs, x->skip_block,
                            p->quant_pxd_fp, pd->dequant_pxd,
                            (const cumbins_type_nuq *)p->cumbins_nuq_pxd,
                            (const dequant_val_type_nuq *)
                            pd->dequant_val_nuq_pxd,
                            qcoeff, dqcoeff, eob,
                            scan_order->scan, band);
    } else if (tx_size == TX_32X32) {
#if CONFIG_VP9_HIGHBITDEPTH
      if (xd->cur_buf->flags & YV12_FLAG_HIGHBITDEPTH)
        vp9_highbd_quantize_32x32_fp_nuq(coeff, bs * bs, x->skip_block,
                                         p->quant_pxd_fp, pd->dequant_pxd,
                                         (const cumbins_type_nuq *)
                                         p->cumbins_nuq_pxd,
                                         (const dequant_val_type_nuq *)
                                         pd->dequant_val_nuq_pxd,
                                         qcoeff, dqcoeff, eob,
                                         scan_order->scan, band);
      else
#endif  // CONFIG_VP9_HIGHBITDEPTH
        vp9_quantize_32x32_fp_nuq(coeff, bs * bs, x->skip_block,
                                  p->quant_pxd_fp, pd->dequant_pxd,
                                  (const cumbins_type_nuq *)p->cumbins_nuq_pxd,
                                  (const dequant_val_type_nuq *)
                                  pd->dequant_val_nuq_pxd,
                                  qcoeff, dqcoeff, eob,
                                  scan_order->scan, band);
    }
#if CONFIG_TX64X64
    else if (tx_size == TX_64X64) {
#if CONFIG_VP9_HIGHBITDEPTH
      if (xd->cur_buf->flags & YV12_FLAG_HIGHBITDEPTH)
        vp9_highbd_quantize_64x64_fp_nuq(coeff, bs * bs, x->skip_block,
                                         p->quant_pxd_fp, pd->dequant_pxd,
                                         (const cumbins_type_nuq *)
                                         p->cumbins_nuq_pxd,
                                         (const dequant_val_type_nuq *)
                                         pd->dequant_val_nuq_pxd,
                                         qcoeff, dqcoeff, eob,
                                         scan_order->scan, band);
      else
#endif  // CONFIG_VP9_HIGHBITDEPTH
        vp9_quantize_64x64_fp_nuq(coeff, bs * bs, x->skip_block,
                                  p->quant_pxd_fp, pd->dequant_pxd,
                                  (const cumbins_type_nuq *)p->cumbins_nuq_pxd,
                                  (const dequant_val_type_nuq *)
                                  pd->dequant_val_nuq_pxd,
                                  qcoeff, dqcoeff, eob,
                                  scan_order->scan, band);
    }
#endif  // CONFIG_TX64X64
    return;
  }
#endif  // CONFIG_TX_SKIP

#if CONFIG_VP9_HIGHBITDEPTH
  if (xd->cur_buf->flags & YV12_FLAG_HIGHBITDEPTH) {
    switch (tx_size) {
#if CONFIG_TX64X64
      case TX_64X64:
        vp9_highbd_fdct64x64(src_diff, coeff, diff_stride);
        vp9_highbd_quantize_64x64_fp_nuq(coeff, 4096, x->skip_block,
                                         p->quant_fp, pd->dequant,
                                         (const cumbins_type_nuq *)
                                             p->cumbins_nuq,
                                         (const dequant_val_type_nuq *)
                                             pd->dequant_val_nuq,
                                         qcoeff, dqcoeff, eob,
                                         scan_order->scan, band);
        break;
#endif  // CONFIG_TX64X64
      case TX_32X32:
        highbd_fdct32x32(x->use_lp32x32fdct, src_diff, coeff, diff_stride);
        vp9_highbd_quantize_32x32_fp_nuq(coeff, 1024, x->skip_block,
                                         p->quant_fp, pd->dequant,
                                         (const cumbins_type_nuq *)
                                             p->cumbins_nuq,
                                         (const dequant_val_type_nuq *)
                                             pd->dequant_val_nuq,
                                         qcoeff, dqcoeff, eob,
                                         scan_order->scan, band);
        break;
      case TX_16X16:
#if CONFIG_EXT_TX
        highbd_forw_tx16x16(x, plane, src_diff, diff_stride, coeff);
#else
        vp9_highbd_fdct16x16(src_diff, coeff, diff_stride);
#endif
        vp9_highbd_quantize_fp_nuq(coeff, 256, x->skip_block,
                                   p->quant_fp, pd->dequant,
                                   (const cumbins_type_nuq *)p->cumbins_nuq,
                                   (const dequant_val_type_nuq *)
                                       pd->dequant_val_nuq,
                                   qcoeff, dqcoeff, eob,
                                   scan_order->scan, band);
        break;
      case TX_8X8:
#if CONFIG_EXT_TX
        highbd_forw_tx8x8(x, plane, src_diff, diff_stride, coeff);
#else
        vp9_highbd_fdct8x8(src_diff, coeff, diff_stride);
#endif
        vp9_highbd_quantize_fp_nuq(coeff, 64, x->skip_block,
                                   p->quant_fp, pd->dequant,
                                   (const cumbins_type_nuq *)p->cumbins_nuq,
                                   (const dequant_val_type_nuq *)
                                       pd->dequant_val_nuq,
                                   qcoeff, dqcoeff, eob,
                                   scan_order->scan, band);
        break;
      case TX_4X4:
#if CONFIG_EXT_TX
        highbd_forw_tx4x4(x, plane, block, src_diff, diff_stride, coeff);
#else
        x->fwd_txm4x4(src_diff, coeff, diff_stride);
#endif
        vp9_highbd_quantize_fp_nuq(coeff, 16, x->skip_block,
                                   p->quant_fp, pd->dequant,
                                   (const cumbins_type_nuq *)p->cumbins_nuq,
                                   (const dequant_val_type_nuq *)
                                       pd->dequant_val_nuq,
                                   qcoeff, dqcoeff, eob,
                                   scan_order->scan, band);
        break;
      default:
        assert(0);
    }
    return;
  }
#endif  // CONFIG_VP9_HIGHBITDEPTH

  switch (tx_size) {
#if CONFIG_TX64X64
    case TX_64X64:
      vp9_fdct64x64(src_diff, coeff, diff_stride);
      vp9_quantize_64x64_fp_nuq(coeff, 4096, x->skip_block,
                                p->quant_fp, pd->dequant,
                                (const cumbins_type_nuq *)p->cumbins_nuq,
                                (const dequant_val_type_nuq *)
                                    pd->dequant_val_nuq,
                                qcoeff, dqcoeff, eob,
                                scan_order->scan, band);
      break;
#endif  // CONFIG_TX64X64
    case TX_32X32:
      fdct32x32(x->use_lp32x32fdct, src_diff, coeff, diff_stride);
      vp9_quantize_32x32_fp_nuq(coeff, 1024, x->skip_block,
                                p->quant_fp, pd->dequant,
                                (const cumbins_type_nuq *)p->cumbins_nuq,
                                (const dequant_val_type_nuq *)
                                    pd->dequant_val_nuq,
                                qcoeff, dqcoeff, eob,
                                scan_order->scan, band);
      break;
    case TX_16X16:
#if CONFIG_EXT_TX
      forw_tx16x16(x, plane, src_diff, diff_stride, coeff);
#else
      vp9_fdct16x16(src_diff, coeff, diff_stride);
#endif
      vp9_quantize_fp_nuq(coeff, 256, x->skip_block,
                          p->quant_fp, pd->dequant,
                          (const cumbins_type_nuq *)p->cumbins_nuq,
                          (const dequant_val_type_nuq *)pd->dequant_val_nuq,
                          qcoeff, dqcoeff, eob,
                          scan_order->scan, band);
      break;
    case TX_8X8:
#if CONFIG_EXT_TX
      forw_tx8x8(x, plane, src_diff, diff_stride, coeff);
#else
      vp9_fdct8x8(src_diff, coeff, diff_stride);
#endif
      vp9_quantize_fp_nuq(coeff, 64, x->skip_block,
                          p->quant_fp, pd->dequant,
                          (const cumbins_type_nuq *)p->cumbins_nuq,
                          (const dequant_val_type_nuq *)pd->dequant_val_nuq,
                          qcoeff, dqcoeff, eob,
                          scan_order->scan, band);
      break;
    case TX_4X4:
#if CONFIG_EXT_TX
      forw_tx4x4(x, plane, block, src_diff, diff_stride, coeff);
#else
      x->fwd_txm4x4(src_diff, coeff, diff_stride);
#endif
      vp9_quantize_fp_nuq(coeff, 16, x->skip_block,
                          p->quant_fp, pd->dequant,
                          (const cumbins_type_nuq *)p->cumbins_nuq,
                          (const dequant_val_type_nuq *)pd->dequant_val_nuq,
                          qcoeff, dqcoeff, eob,
                          scan_order->scan, band);
      break;
    default:
      assert(0);
      break;
  }
}

void vp9_xform_quant_dc_nuq(MACROBLOCK *x, int plane, int block,
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
#if CONFIG_TX_SKIP
  MB_MODE_INFO *mbmi = &xd->mi[0].src_mi->mbmi;
  int shift = mbmi->tx_skip_shift;
#endif

  txfrm_block_to_raster_xy(plane_bsize, tx_size, block, &i, &j);
  src_diff = &p->src_diff[4 * (j * diff_stride + i)];

#if CONFIG_TX_SKIP
  if (mbmi->tx_skip[plane != 0]) {
    int bs = 4 << tx_size;
    vp9_tx_identity(src_diff, coeff, diff_stride, bs, shift);
    if (tx_size <= TX_16X16) {
#if CONFIG_VP9_HIGHBITDEPTH
      if (xd->cur_buf->flags & YV12_FLAG_HIGHBITDEPTH)
        vp9_highbd_quantize_dc_nuq(coeff, x->skip_block,
                                   p->quant_pxd[0], p->quant_shift_pxd[0],
                                   pd->dequant_pxd[0],
                                   p->cumbins_nuq_pxd[0],
                                   pd->dequant_val_nuq_pxd[0],
                                   qcoeff, dqcoeff, eob);
#endif  // CONFIG_VP9_HIGHBITDEPTH
      vp9_quantize_dc_nuq(coeff, x->skip_block,
                          p->quant_pxd[0], p->quant_shift_pxd[0],
                          pd->dequant_pxd[0],
                          p->cumbins_nuq_pxd[0], pd->dequant_val_nuq_pxd[0],
                          qcoeff, dqcoeff, eob);
    } else if (tx_size == TX_32X32) {
#if CONFIG_VP9_HIGHBITDEPTH
      if (xd->cur_buf->flags & YV12_FLAG_HIGHBITDEPTH)
        vp9_highbd_quantize_dc_32x32_nuq(coeff, x->skip_block,
                                         p->quant_pxd[0], p->quant_shift_pxd[0],
                                         pd->dequant_pxd[0],
                                         p->cumbins_nuq_pxd[0],
                                         pd->dequant_val_nuq_pxd[0],
                                         qcoeff, dqcoeff, eob);
      else
#endif  // CONFIG_VP9_HIGHBITDEPTH
        vp9_quantize_dc_32x32_nuq(coeff, x->skip_block,
                                  p->quant_pxd[0], p->quant_shift_pxd[0],
                                  pd->dequant_pxd[0],
                                  p->cumbins_nuq_pxd[0],
                                  pd->dequant_val_nuq_pxd[0],
                                  qcoeff, dqcoeff, eob);
    }
#if CONFIG_TX64X64
    else if (tx_size == TX_64X64) {
#if CONFIG_VP9_HIGHBITDEPTH
      if (xd->cur_buf->flags & YV12_FLAG_HIGHBITDEPTH)
        vp9_highbd_quantize_dc_64x64_nuq(coeff, x->skip_block,
                                         p->quant_pxd[0], p->quant_shift_pxd[0],
                                         pd->dequant_pxd[0],
                                         p->cumbins_nuq_pxd[0],
                                         pd->dequant_val_nuq_pxd[0],
                                         qcoeff, dqcoeff, eob);
      else
#endif  // CONFIG_VP9_HIGHBITDEPTH
        vp9_quantize_dc_64x64_nuq(coeff, x->skip_block,
                                  p->quant_pxd[0], p->quant_shift_pxd[0],
                                  pd->dequant_pxd[0],
                                  p->cumbins_nuq_pxd[0],
                                  pd->dequant_val_nuq_pxd[0],
                                  qcoeff, dqcoeff, eob);
    }
#endif  // CONFIG_TX64X64

    return;
  }
#endif  // CONFIG_TX_SKIP

#if CONFIG_VP9_HIGHBITDEPTH
  if (xd->cur_buf->flags & YV12_FLAG_HIGHBITDEPTH) {
    switch (tx_size) {
#if CONFIG_TX64X64
      case TX_64X64:
        vp9_highbd_fdct64x64_1(src_diff, coeff, diff_stride);
        vp9_highbd_quantize_dc_64x64_nuq(coeff, x->skip_block,
                                         p->quant[0], p->quant_shift[0],
                                         pd->dequant[0],
                                         p->cumbins_nuq[0],
                                         pd->dequant_val_nuq[0],
                                         qcoeff, dqcoeff, eob);
        break;
#endif  // CONFIG_TX64X64
      case TX_32X32:
        vp9_highbd_fdct32x32_1(src_diff, coeff, diff_stride);
        vp9_highbd_quantize_dc_32x32_nuq(coeff, x->skip_block,
                                         p->quant[0], p->quant_shift[0],
                                         pd->dequant[0],
                                         p->cumbins_nuq[0],
                                         pd->dequant_val_nuq[0],
                                         qcoeff, dqcoeff, eob);
        break;
      case TX_16X16:
#if CONFIG_EXT_TX
        highbd_forw_tx16x16(x, plane, src_diff, diff_stride, coeff);
#else
        vp9_highbd_fdct16x16_1(src_diff, coeff, diff_stride);
#endif
        vp9_highbd_quantize_dc_nuq(coeff, x->skip_block,
                                   p->quant[0], p->quant_shift[0],
                                   pd->dequant[0],
                                   p->cumbins_nuq[0],
                                   pd->dequant_val_nuq[0],
                                   qcoeff, dqcoeff, eob);
        break;
      case TX_8X8:
#if CONFIG_EXT_TX
        highbd_forw_tx8x8(x, plane, src_diff, diff_stride, coeff);
#else
        vp9_highbd_fdct8x8_1(src_diff, coeff, diff_stride);
#endif
        vp9_highbd_quantize_dc_nuq(coeff, x->skip_block,
                                   p->quant[0], p->quant_shift[0],
                                   pd->dequant[0],
                                   p->cumbins_nuq[0],
                                   pd->dequant_val_nuq[0],
                                   qcoeff, dqcoeff, eob);
        break;
      case TX_4X4:
#if CONFIG_EXT_TX
        highbd_forw_tx4x4(x, plane, block, src_diff, diff_stride, coeff);
#else
        x->fwd_txm4x4(src_diff, coeff, diff_stride);
#endif
        vp9_highbd_quantize_dc_nuq(coeff, x->skip_block,
                                   p->quant[0], p->quant_shift[0],
                                   pd->dequant[0],
                                   p->cumbins_nuq[0],
                                   pd->dequant_val_nuq[0],
                                   qcoeff, dqcoeff, eob);
        break;
      default:
        assert(0);
    }
    return;
  }
#endif  // CONFIG_VP9_HIGHBITDEPTH

  switch (tx_size) {
#if CONFIG_TX64X64
    case TX_64X64:
      vp9_fdct64x64_1(src_diff, coeff, diff_stride);
      vp9_quantize_dc_64x64_nuq(coeff, x->skip_block,
                                p->quant[0], p->quant_shift[0], pd->dequant[0],
                                p->cumbins_nuq[0], pd->dequant_val_nuq[0],
                                qcoeff, dqcoeff, eob);
      break;
#endif  // CONFIG_TX64X64
    case TX_32X32:
      vp9_fdct32x32_1(src_diff, coeff, diff_stride);
      vp9_quantize_dc_32x32_nuq(coeff, x->skip_block,
                                p->quant[0], p->quant_shift[0], pd->dequant[0],
                                p->cumbins_nuq[0], pd->dequant_val_nuq[0],
                                qcoeff, dqcoeff, eob);
      break;
    case TX_16X16:
#if CONFIG_EXT_TX
      forw_tx16x16(x, plane, src_diff, diff_stride, coeff);
#else
      vp9_fdct16x16_1(src_diff, coeff, diff_stride);
#endif
      vp9_quantize_dc_nuq(coeff, x->skip_block,
                          p->quant[0], p->quant_shift[0], pd->dequant[0],
                          p->cumbins_nuq[0], pd->dequant_val_nuq[0],
                          qcoeff, dqcoeff, eob);
      break;
    case TX_8X8:
#if CONFIG_EXT_TX
      forw_tx8x8(x, plane, src_diff, diff_stride, coeff);
#else
      vp9_fdct8x8_1(src_diff, coeff, diff_stride);
#endif
      vp9_quantize_dc_nuq(coeff, x->skip_block,
                          p->quant[0], p->quant_shift[0], pd->dequant[0],
                          p->cumbins_nuq[0], pd->dequant_val_nuq[0],
                          qcoeff, dqcoeff, eob);
      break;
    case TX_4X4:
#if CONFIG_EXT_TX
      forw_tx4x4(x, plane, block, src_diff, diff_stride, coeff);
#else
      x->fwd_txm4x4(src_diff, coeff, diff_stride);
#endif
      vp9_quantize_dc_nuq(coeff, x->skip_block,
                          p->quant[0], p->quant_shift[0], pd->dequant[0],
                          p->cumbins_nuq[0], pd->dequant_val_nuq[0],
                          qcoeff, dqcoeff, eob);
      break;
    default:
      assert(0);
      break;
  }
}

void vp9_xform_quant_dc_fp_nuq(MACROBLOCK *x, int plane, int block,
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
#if CONFIG_TX_SKIP
  MB_MODE_INFO *mbmi = &xd->mi[0].src_mi->mbmi;
  int shift = mbmi->tx_skip_shift;
#endif

  txfrm_block_to_raster_xy(plane_bsize, tx_size, block, &i, &j);
  src_diff = &p->src_diff[4 * (j * diff_stride + i)];

#if CONFIG_TX_SKIP
  if (mbmi->tx_skip[plane != 0]) {
    int bs = 4 << tx_size;
    vp9_tx_identity(src_diff, coeff, diff_stride, bs, shift);
    if (tx_size <= TX_16X16) {
#if CONFIG_VP9_HIGHBITDEPTH
      if (xd->cur_buf->flags & YV12_FLAG_HIGHBITDEPTH)
        vp9_highbd_quantize_dc_fp_nuq(coeff, x->skip_block,
                                      p->quant_pxd_fp[0], pd->dequant_pxd[0],
                                      p->cumbins_nuq_pxd[0],
                                      pd->dequant_val_nuq_pxd[0],
                                      qcoeff, dqcoeff, eob);
      else
#endif  // CONFIG_VP9_HIGHBITDEPTH
        vp9_quantize_dc_fp_nuq(coeff, x->skip_block,
                               p->quant_pxd_fp[0], pd->dequant_pxd[0],
                               p->cumbins_nuq_pxd[0],
                               pd->dequant_val_nuq_pxd[0],
                               qcoeff, dqcoeff, eob);
    } else if (tx_size == TX_32X32) {
#if CONFIG_VP9_HIGHBITDEPTH
      if (xd->cur_buf->flags & YV12_FLAG_HIGHBITDEPTH)
        vp9_highbd_quantize_dc_32x32_fp_nuq(coeff, x->skip_block,
                                            p->quant_pxd_fp[0],
                                            pd->dequant_pxd[0],
                                            p->cumbins_nuq_pxd[0],
                                            pd->dequant_val_nuq_pxd[0],
                                            qcoeff, dqcoeff, eob);
      else
#endif  // CONFIG_VP9_HIGHBITDEPTH
        vp9_quantize_dc_32x32_fp_nuq(coeff, x->skip_block,
                                     p->quant_pxd_fp[0], pd->dequant_pxd[0],
                                     p->cumbins_nuq_pxd[0],
                                     pd->dequant_val_nuq_pxd[0],
                                     qcoeff, dqcoeff, eob);
    }
#if CONFIG_TX64X64
    else if (tx_size == TX_64X64) {
#if CONFIG_VP9_HIGHBITDEPTH
      if (xd->cur_buf->flags & YV12_FLAG_HIGHBITDEPTH)
        vp9_highbd_quantize_dc_64x64_fp_nuq(coeff, x->skip_block,
                                            p->quant_pxd_fp[0],
                                            pd->dequant_pxd[0],
                                            p->cumbins_nuq_pxd[0],
                                            pd->dequant_val_nuq_pxd[0],
                                            qcoeff, dqcoeff, eob);
      else
#endif  // CONFIG_VP9_HIGHBITDEPTH
        vp9_quantize_dc_64x64_fp_nuq(coeff, x->skip_block,
                                     p->quant_pxd_fp[0], pd->dequant_pxd[0],
                                     p->cumbins_nuq_pxd[0],
                                     pd->dequant_val_nuq_pxd[0],
                                     qcoeff, dqcoeff, eob);
    }
#endif  // CONFIG_TX64X64

    return;
  }
#endif  // CONFIG_TX_SKIP

#if CONFIG_VP9_HIGHBITDEPTH
  if (xd->cur_buf->flags & YV12_FLAG_HIGHBITDEPTH) {
    switch (tx_size) {
#if CONFIG_TX64X64
      case TX_64X64:
        vp9_highbd_fdct64x64_1(src_diff, coeff, diff_stride);
        vp9_highbd_quantize_dc_64x64_fp_nuq(coeff, x->skip_block,
                                            p->quant_fp[0], pd->dequant[0],
                                            p->cumbins_nuq[0],
                                            pd->dequant_val_nuq[0],
                                            qcoeff, dqcoeff, eob);
        break;
#endif  // CONFIG_TX64X64
      case TX_32X32:
        vp9_highbd_fdct32x32_1(src_diff, coeff, diff_stride);
        vp9_highbd_quantize_dc_32x32_fp_nuq(coeff, x->skip_block,
                                            p->quant_fp[0], pd->dequant[0],
                                            p->cumbins_nuq[0],
                                            pd->dequant_val_nuq[0],
                                            qcoeff, dqcoeff, eob);
        break;
      case TX_16X16:
#if CONFIG_EXT_TX
        highbd_forw_tx16x16(x, plane, src_diff, diff_stride, coeff);
#else
        vp9_highbd_fdct16x16_1(src_diff, coeff, diff_stride);
#endif
        vp9_highbd_quantize_dc_fp_nuq(coeff, x->skip_block,
                                      p->quant_fp[0], pd->dequant[0],
                                      p->cumbins_nuq[0],
                                      pd->dequant_val_nuq[0],
                                      qcoeff, dqcoeff, eob);
        break;
      case TX_8X8:
#if CONFIG_EXT_TX
        highbd_forw_tx8x8(x, plane, src_diff, diff_stride, coeff);
#else
        vp9_highbd_fdct8x8_1(src_diff, coeff, diff_stride);
#endif
        vp9_highbd_quantize_dc_fp_nuq(coeff, x->skip_block,
                                      p->quant_fp[0], pd->dequant[0],
                                      p->cumbins_nuq[0],
                                      pd->dequant_val_nuq[0],
                                      qcoeff, dqcoeff, eob);
        break;
      case TX_4X4:
#if CONFIG_EXT_TX
        highbd_forw_tx4x4(x, plane, block, src_diff, diff_stride, coeff);
#else
        x->fwd_txm4x4(src_diff, coeff, diff_stride);
#endif
        vp9_highbd_quantize_dc_fp_nuq(coeff, x->skip_block,
                                      p->quant_fp[0], pd->dequant[0],
                                      p->cumbins_nuq[0],
                                      pd->dequant_val_nuq[0],
                                      qcoeff, dqcoeff, eob);
        break;
      default:
        assert(0);
    }
    return;
  }
#endif  // CONFIG_VP9_HIGHBITDEPTH

  switch (tx_size) {
#if CONFIG_TX64X64
    case TX_64X64:
      vp9_fdct64x64_1(src_diff, coeff, diff_stride);
      vp9_quantize_dc_64x64_fp_nuq(coeff, x->skip_block,
                                   p->quant_fp[0], pd->dequant[0],
                                   p->cumbins_nuq[0], pd->dequant_val_nuq[0],
                                   qcoeff, dqcoeff, eob);
      break;
#endif  // CONFIG_TX64X64
    case TX_32X32:
      vp9_fdct32x32_1(src_diff, coeff, diff_stride);
      vp9_quantize_dc_32x32_fp_nuq(coeff, x->skip_block,
                                   p->quant_fp[0], pd->dequant[0],
                                   p->cumbins_nuq[0], pd->dequant_val_nuq[0],
                                   qcoeff, dqcoeff, eob);
      break;
    case TX_16X16:
#if CONFIG_EXT_TX
      forw_tx16x16(x, plane, src_diff, diff_stride, coeff);
#else
      vp9_fdct16x16_1(src_diff, coeff, diff_stride);
#endif
      vp9_quantize_dc_fp_nuq(coeff, x->skip_block,
                             p->quant_fp[0], pd->dequant[0],
                             p->cumbins_nuq[0], pd->dequant_val_nuq[0],
                             qcoeff, dqcoeff, eob);
      break;
    case TX_8X8:
#if CONFIG_EXT_TX
      forw_tx8x8(x, plane, src_diff, diff_stride, coeff);
#else
      vp9_fdct8x8_1(src_diff, coeff, diff_stride);
#endif
      vp9_quantize_dc_fp_nuq(coeff, x->skip_block,
                             p->quant_fp[0], pd->dequant[0],
                             p->cumbins_nuq[0], pd->dequant_val_nuq[0],
                             qcoeff, dqcoeff, eob);
      break;
    case TX_4X4:
#if CONFIG_EXT_TX
      forw_tx4x4(x, plane, block, src_diff, diff_stride, coeff);
#else
      x->fwd_txm4x4(src_diff, coeff, diff_stride);
#endif
      vp9_quantize_dc_fp_nuq(coeff, x->skip_block,
                             p->quant_fp[0], pd->dequant[0],
                             p->cumbins_nuq[0], pd->dequant_val_nuq[0],
                             qcoeff, dqcoeff, eob);
      break;
    default:
      assert(0);
      break;
  }
}
#endif  // CONFIG_NEW_QUANT

void vp9_xform_quant_fp(MACROBLOCK *x, int plane, int block,
                        BLOCK_SIZE plane_bsize, TX_SIZE tx_size) {
  MACROBLOCKD *const xd = &x->e_mbd;
  const struct macroblock_plane *const p = &x->plane[plane];
  const struct macroblockd_plane *const pd = &xd->plane[plane];
#if CONFIG_TX_SKIP
  MB_MODE_INFO *mbmi = &xd->mi[0].src_mi->mbmi;
  int shift = mbmi->tx_skip_shift;
  const scan_order *const scan_order =
      mbmi->tx_skip[plane != 0] ? &vp9_default_scan_orders_pxd[tx_size] :
          &vp9_default_scan_orders[tx_size];
#else
  const scan_order *const scan_order = &vp9_default_scan_orders[tx_size];
#endif  // CONFIG_TX_SKIP
  tran_low_t *const coeff = BLOCK_OFFSET(p->coeff, block);
  tran_low_t *const qcoeff = BLOCK_OFFSET(p->qcoeff, block);
  tran_low_t *const dqcoeff = BLOCK_OFFSET(pd->dqcoeff, block);
  uint16_t *const eob = &p->eobs[block];
  const int diff_stride = 4 * num_4x4_blocks_wide_lookup[plane_bsize];
  int i, j;
  const int16_t *src_diff;

  txfrm_block_to_raster_xy(plane_bsize, tx_size, block, &i, &j);
  src_diff = &p->src_diff[4 * (j * diff_stride + i)];

#if CONFIG_TX_SKIP
  if (mbmi->tx_skip[plane != 0]) {
    int bs = 4 << tx_size;
#if CONFIG_VP9_HIGHBITDEPTH
    int use_hbd = xd->cur_buf->flags & YV12_FLAG_HIGHBITDEPTH;
#endif  // CONFIG_VP9_HIGHBITDEPTH
    vp9_tx_identity(src_diff, coeff, diff_stride, bs, shift);
    if (tx_size <= TX_16X16) {
#if CONFIG_VP9_HIGHBITDEPTH
      if (use_hbd)
        vp9_highbd_quantize_fp(coeff, bs * bs, x->skip_block, p->zbin_pxd,
                               p->round_pxd, p->quant_pxd, p->quant_shift_pxd,
                               qcoeff, dqcoeff, pd->dequant_pxd, eob,
                               scan_order->scan, scan_order->iscan);
      else
#endif  // CONFIG_VP9_HIGHBITDEPTH
        vp9_quantize_fp(coeff, bs * bs, x->skip_block, p->zbin_pxd,
                        p->round_pxd, p->quant_pxd, p->quant_shift_pxd, qcoeff,
                        dqcoeff, pd->dequant_pxd, eob,
                        scan_order->scan, scan_order->iscan);
    } else if (tx_size == TX_32X32) {
#if CONFIG_VP9_HIGHBITDEPTH
      if (use_hbd)
        vp9_highbd_quantize_fp_32x32(coeff, bs * bs, x->skip_block, p->zbin_pxd,
                                     p->round_pxd, p->quant_pxd,
                                     p->quant_shift_pxd,
                                     qcoeff, dqcoeff, pd->dequant_pxd, eob,
                                     scan_order->scan, scan_order->iscan);
      else
#endif  // CONFIG_VP9_HIGHBITDEPTH
        vp9_quantize_fp_32x32(coeff, bs * bs, x->skip_block, p->zbin_pxd,
                              p->round_pxd, p->quant_pxd, p->quant_shift_pxd,
                              qcoeff, dqcoeff, pd->dequant_pxd, eob,
                              scan_order->scan, scan_order->iscan);
    }
#if CONFIG_TX64X64
    else if (tx_size == TX_64X64) {
#if CONFIG_VP9_HIGHBITDEPTH
      if (use_hbd)
        vp9_highbd_quantize_fp_64x64(coeff, bs * bs, x->skip_block, p->zbin_pxd,
                                     p->round_pxd, p->quant_pxd,
                                     p->quant_shift_pxd, qcoeff, dqcoeff,
                                     pd->dequant_pxd, eob,
                                     scan_order->scan, scan_order->iscan);
      else
#endif  // CONFIG_VP9_HIGHBITDEPTH
        vp9_quantize_fp_64x64(coeff, bs * bs, x->skip_block, p->zbin_pxd,
                              p->round_pxd, p->quant_pxd, p->quant_shift_pxd,
                              qcoeff, dqcoeff, pd->dequant_pxd, eob,
                              scan_order->scan, scan_order->iscan);
    }
#endif  // CONFIG_TX64X64

    return;
  }
#endif  // CONFIG_TX_SKIP

#if CONFIG_VP9_HIGHBITDEPTH
  if (xd->cur_buf->flags & YV12_FLAG_HIGHBITDEPTH) {
    switch (tx_size) {
#if CONFIG_TX64X64
      case TX_64X64:
        vp9_highbd_fdct64x64(src_diff, coeff, diff_stride);
        vp9_highbd_quantize_fp_64x64(coeff, 4096, x->skip_block, p->zbin,
                                     p->round_fp, p->quant_fp, p->quant_shift,
                                     qcoeff, dqcoeff, pd->dequant,
                                     eob, scan_order->scan,
                                     scan_order->iscan);
        break;
#endif  // CONFIG_TX64X64
      case TX_32X32:
        highbd_fdct32x32(x->use_lp32x32fdct, src_diff, coeff, diff_stride);
        vp9_highbd_quantize_fp_32x32(coeff, 1024, x->skip_block, p->zbin,
                                     p->round_fp, p->quant_fp, p->quant_shift,
                                     qcoeff, dqcoeff, pd->dequant,
                                     eob, scan_order->scan,
                                     scan_order->iscan);
        break;
      case TX_16X16:
#if CONFIG_EXT_TX
        highbd_forw_tx16x16(x, plane, src_diff, diff_stride, coeff);
#else
        vp9_highbd_fdct16x16(src_diff, coeff, diff_stride);
#endif
        vp9_highbd_quantize_fp(coeff, 256, x->skip_block, p->zbin, p->round_fp,
                               p->quant_fp, p->quant_shift, qcoeff, dqcoeff,
                               pd->dequant, eob,
                               scan_order->scan, scan_order->iscan);
        break;
      case TX_8X8:
#if CONFIG_EXT_TX
        highbd_forw_tx8x8(x, plane, src_diff, diff_stride, coeff);
#else
        vp9_highbd_fdct8x8(src_diff, coeff, diff_stride);
#endif
        vp9_highbd_quantize_fp(coeff, 64, x->skip_block, p->zbin, p->round_fp,
                               p->quant_fp, p->quant_shift, qcoeff, dqcoeff,
                               pd->dequant, eob,
                               scan_order->scan, scan_order->iscan);
        break;
      case TX_4X4:
#if CONFIG_EXT_TX
        highbd_forw_tx4x4(x, plane, block, src_diff, diff_stride, coeff);
#else
        x->fwd_txm4x4(src_diff, coeff, diff_stride);
#endif
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
#if CONFIG_TX64X64
    case TX_64X64:
      vp9_fdct64x64(src_diff, coeff, diff_stride);
      vp9_quantize_fp_64x64(coeff, 4096, x->skip_block, p->zbin, p->round_fp,
                            p->quant_fp, p->quant_shift, qcoeff, dqcoeff,
                            pd->dequant, eob, scan_order->scan,
                            scan_order->iscan);
      break;
#endif  // CONFIG_TX64X64
    case TX_32X32:
      fdct32x32(x->use_lp32x32fdct, src_diff, coeff, diff_stride);
      vp9_quantize_fp_32x32(coeff, 1024, x->skip_block, p->zbin,
                            p->round_fp, p->quant_fp, p->quant_shift,
                            qcoeff, dqcoeff, pd->dequant, eob,
                            scan_order->scan, scan_order->iscan);
      break;
    case TX_16X16:
#if CONFIG_EXT_TX
      forw_tx16x16(x, plane, src_diff, diff_stride, coeff);
#else
      vp9_fdct16x16(src_diff, coeff, diff_stride);
#endif
      vp9_quantize_fp(coeff, 256, x->skip_block, p->zbin, p->round_fp,
                      p->quant_fp, p->quant_shift, qcoeff, dqcoeff,
                      pd->dequant, eob,
                      scan_order->scan, scan_order->iscan);
      break;
    case TX_8X8:
#if CONFIG_EXT_TX
      forw_tx8x8(x, plane, src_diff, diff_stride, coeff);
#else
      vp9_fdct8x8(src_diff, coeff, diff_stride);
#endif
      vp9_quantize_fp(coeff, 64, x->skip_block, p->zbin, p->round_fp,
                      p->quant_fp, p->quant_shift, qcoeff, dqcoeff,
                      pd->dequant, eob,
                      scan_order->scan, scan_order->iscan);
      break;
    case TX_4X4:
#if CONFIG_EXT_TX
      forw_tx4x4(x, plane, block, src_diff, diff_stride, coeff);
#else
      x->fwd_txm4x4(src_diff, coeff, diff_stride);
#endif
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
#if CONFIG_TX_SKIP
  MB_MODE_INFO *mbmi = &xd->mi[0].src_mi->mbmi;
  int shift = mbmi->tx_skip_shift;
#endif

  txfrm_block_to_raster_xy(plane_bsize, tx_size, block, &i, &j);
  src_diff = &p->src_diff[4 * (j * diff_stride + i)];

#if CONFIG_TX_SKIP
  if (mbmi->tx_skip[plane != 0]) {
    int bs = 4 << tx_size;
#if CONFIG_VP9_HIGHBITDEPTH
    int use_hbd = xd->cur_buf->flags & YV12_FLAG_HIGHBITDEPTH;
#endif  // CONFIG_VP9_HIGHBITDEPTH
    vp9_tx_identity(src_diff, coeff, diff_stride, bs, shift);
    if (tx_size <= TX_16X16) {
#if CONFIG_VP9_HIGHBITDEPTH
      if (use_hbd)
        vp9_highbd_quantize_dc(coeff, x->skip_block, p->round_pxd,
                               p->quant_pxd_fp[0], qcoeff, dqcoeff,
                               pd->dequant_pxd[0], eob);
      else
#endif  // CONFIG_VP9_HIGHBITDEPTH
        vp9_quantize_dc(coeff, x->skip_block, p->round_pxd,
                        p->quant_pxd_fp[0], qcoeff, dqcoeff,
                        pd->dequant_pxd[0], eob);
    } else if (tx_size == TX_32X32) {
#if CONFIG_VP9_HIGHBITDEPTH
      if (use_hbd)
        vp9_highbd_quantize_dc_32x32(coeff, x->skip_block, p->round_pxd,
                                     p->quant_pxd_fp[0], qcoeff, dqcoeff,
                                     pd->dequant_pxd[0], eob);
      else
#endif  // CONFIG_VP9_HIGHBITDEPTH
        vp9_quantize_dc_32x32(coeff, x->skip_block, p->round_pxd,
                              p->quant_pxd_fp[0], qcoeff, dqcoeff,
                              pd->dequant_pxd[0], eob);
    }
#if CONFIG_TX64X64
    else if (tx_size == TX_64X64) {
#if CONFIG_VP9_HIGHBITDEPTH
      if (use_hbd)
        vp9_highbd_quantize_dc_64x64(coeff, x->skip_block, p->round_pxd,
                                     p->quant_pxd_fp[0], qcoeff, dqcoeff,
                                     pd->dequant_pxd[0], eob);
      else
#endif  // CONFIG_VP9_HIGHBITDEPTH
        vp9_quantize_dc_64x64(coeff, x->skip_block, p->round_pxd,
                              p->quant_pxd_fp[0], qcoeff, dqcoeff,
                              pd->dequant_pxd[0], eob);
    }
#endif  // CONFIG_TX64X64

    return;
  }
#endif  // CONFIG_TX_SKIP

#if CONFIG_VP9_HIGHBITDEPTH
  if (xd->cur_buf->flags & YV12_FLAG_HIGHBITDEPTH) {
    switch (tx_size) {
#if CONFIG_TX64X64
      case TX_64X64:
        vp9_highbd_fdct64x64_1(src_diff, coeff, diff_stride);
        vp9_highbd_quantize_dc_64x64(coeff, x->skip_block, p->round,
                                     p->quant_fp[0], qcoeff, dqcoeff,
                                     pd->dequant[0], eob);
        break;
#endif  // CONFIG_TX64X64
      case TX_32X32:
        vp9_highbd_fdct32x32_1(src_diff, coeff, diff_stride);
        vp9_highbd_quantize_dc_32x32(coeff, x->skip_block, p->round,
                                     p->quant_fp[0], qcoeff, dqcoeff,
                                     pd->dequant[0], eob);
        break;
      case TX_16X16:
#if CONFIG_EXT_TX
        highbd_forw_tx16x16(x, plane, src_diff, diff_stride, coeff);
#else
        vp9_highbd_fdct16x16_1(src_diff, coeff, diff_stride);
#endif
        vp9_highbd_quantize_dc(coeff, x->skip_block, p->round,
                               p->quant_fp[0], qcoeff, dqcoeff,
                               pd->dequant[0], eob);
        break;
      case TX_8X8:
#if CONFIG_EXT_TX
        highbd_forw_tx8x8(x, plane, src_diff, diff_stride, coeff);
#else
        vp9_highbd_fdct8x8_1(src_diff, coeff, diff_stride);
#endif
        vp9_highbd_quantize_dc(coeff, x->skip_block, p->round,
                               p->quant_fp[0], qcoeff, dqcoeff,
                               pd->dequant[0], eob);
        break;
      case TX_4X4:
#if CONFIG_EXT_TX
        highbd_forw_tx4x4(x, plane, block, src_diff, diff_stride, coeff);
#else
        x->fwd_txm4x4(src_diff, coeff, diff_stride);
#endif
        vp9_highbd_quantize_dc(coeff, x->skip_block, p->round,
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
#if CONFIG_TX64X64
    case TX_64X64:
      vp9_fdct64x64_1(src_diff, coeff, diff_stride);
      vp9_quantize_dc_64x64(coeff, x->skip_block, p->round,
                            p->quant_fp[0], qcoeff, dqcoeff,
                            pd->dequant[0], eob);
      break;
#endif  // CONFIG_TX64X64
    case TX_32X32:
      vp9_fdct32x32_1(src_diff, coeff, diff_stride);
      vp9_quantize_dc_32x32(coeff, x->skip_block, p->round,
                            p->quant_fp[0], qcoeff, dqcoeff,
                            pd->dequant[0], eob);
      break;
    case TX_16X16:
#if CONFIG_EXT_TX
      forw_tx16x16(x, plane, src_diff, diff_stride, coeff);
#else
      vp9_fdct16x16_1(src_diff, coeff, diff_stride);
#endif
      vp9_quantize_dc(coeff, x->skip_block, p->round,
                     p->quant_fp[0], qcoeff, dqcoeff,
                     pd->dequant[0], eob);
      break;
    case TX_8X8:
#if CONFIG_EXT_TX
      forw_tx8x8(x, plane, src_diff, diff_stride, coeff);
#else
      vp9_fdct8x8_1(src_diff, coeff, diff_stride);
#endif
      vp9_quantize_dc(coeff, x->skip_block, p->round,
                      p->quant_fp[0], qcoeff, dqcoeff,
                      pd->dequant[0], eob);
      break;
    case TX_4X4:
#if CONFIG_EXT_TX
      forw_tx4x4(x, plane, block, src_diff, diff_stride, coeff);
#else
      x->fwd_txm4x4(src_diff, coeff, diff_stride);
#endif
      vp9_quantize_dc(coeff, x->skip_block, p->round,
                      p->quant_fp[0], qcoeff, dqcoeff,
                      pd->dequant[0], eob);
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
#if CONFIG_TX_SKIP
  MB_MODE_INFO *mbmi = &xd->mi[0].src_mi->mbmi;
  int shift = mbmi->tx_skip_shift;
  const scan_order *const scan_order =
      mbmi->tx_skip[plane != 0] ? &vp9_default_scan_orders_pxd[tx_size] :
          &vp9_default_scan_orders[tx_size];
#else
  const scan_order *const scan_order = &vp9_default_scan_orders[tx_size];
#endif  // CONFIG_TX_SKIP
  tran_low_t *const coeff = BLOCK_OFFSET(p->coeff, block);
  tran_low_t *const qcoeff = BLOCK_OFFSET(p->qcoeff, block);
  tran_low_t *const dqcoeff = BLOCK_OFFSET(pd->dqcoeff, block);
  uint16_t *const eob = &p->eobs[block];
  const int diff_stride = 4 * num_4x4_blocks_wide_lookup[plane_bsize];
  int i, j;
  const int16_t *src_diff;
  txfrm_block_to_raster_xy(plane_bsize, tx_size, block, &i, &j);
  src_diff = &p->src_diff[4 * (j * diff_stride + i)];

#if CONFIG_TX_SKIP
  if (mbmi->tx_skip[plane != 0]) {
    int bs = 4 << tx_size;
#if CONFIG_VP9_HIGHBITDEPTH
    int use_hbd = xd->cur_buf->flags & YV12_FLAG_HIGHBITDEPTH;
#endif  // CONFIG_VP9_HIGHBITDEPTH
    vp9_tx_identity(src_diff, coeff, diff_stride, bs, shift);
    if (tx_size <= TX_16X16) {
#if CONFIG_VP9_HIGHBITDEPTH
      if (use_hbd)
        vp9_highbd_quantize_b(coeff, bs * bs, x->skip_block, p->zbin_pxd,
                              p->round_pxd, p->quant_pxd, p->quant_shift_pxd,
                              qcoeff, dqcoeff, pd->dequant_pxd, eob,
                              scan_order->scan, scan_order->iscan);
      else
#endif  // CONFIG_VP9_HIGHBITDEPTH
        vp9_quantize_b(coeff, bs * bs, x->skip_block, p->zbin_pxd, p->round_pxd,
                       p->quant_pxd, p->quant_shift_pxd, qcoeff, dqcoeff,
                       pd->dequant_pxd, eob,
                       scan_order->scan, scan_order->iscan);
    } else if (tx_size == TX_32X32) {
#if CONFIG_VP9_HIGHBITDEPTH
      if (use_hbd)
        vp9_highbd_quantize_b_32x32(coeff, bs * bs, x->skip_block, p->zbin_pxd,
                                    p->round_pxd, p->quant_pxd,
                                    p->quant_shift_pxd, qcoeff, dqcoeff,
                                    pd->dequant_pxd, eob,
                                    scan_order->scan, scan_order->iscan);
      else
#endif  // CONFIG_VP9_HIGHBITDEPTH
        vp9_quantize_b_32x32(coeff, bs * bs, x->skip_block, p->zbin_pxd,
                             p->round_pxd, p->quant_pxd, p->quant_shift_pxd,
                             qcoeff, dqcoeff, pd->dequant_pxd, eob,
                             scan_order->scan, scan_order->iscan);
    }
#if CONFIG_TX64X64
    else if (tx_size == TX_64X64) {
#if CONFIG_VP9_HIGHBITDEPTH
      if (use_hbd)
        vp9_highbd_quantize_b_64x64(coeff, bs * bs, x->skip_block, p->zbin_pxd,
                                    p->round_pxd, p->quant_pxd,
                                    p->quant_shift_pxd, qcoeff, dqcoeff,
                                    pd->dequant_pxd, eob, scan_order->scan,
                                    scan_order->iscan);
      else
#endif  // CONFIG_VP9_HIGHBITDEPTH
        vp9_quantize_b_64x64(coeff, bs * bs, x->skip_block, p->zbin_pxd,
                             p->round_pxd, p->quant_pxd, p->quant_shift_pxd,
                             qcoeff, dqcoeff, pd->dequant_pxd, eob,
                             scan_order->scan, scan_order->iscan);
    }
#endif  // CONFIG_TX64X64

    return;
  }
#endif  // CONFIG_TX_SKIP

#if CONFIG_VP9_HIGHBITDEPTH
  if (xd->cur_buf->flags & YV12_FLAG_HIGHBITDEPTH) {
     switch (tx_size) {
#if CONFIG_TX64X64
      case TX_64X64:
        vp9_highbd_fdct64x64(src_diff, coeff, diff_stride);
        vp9_highbd_quantize_b_64x64(coeff, 4096, x->skip_block, p->zbin,
                                    p->round, p->quant, p->quant_shift, qcoeff,
                                    dqcoeff, pd->dequant, eob,
                                    scan_order->scan, scan_order->iscan);
        break;
#endif  // CONFIG_TX64X64
      case TX_32X32:
        highbd_fdct32x32(x->use_lp32x32fdct, src_diff, coeff, diff_stride);
        vp9_highbd_quantize_b_32x32(coeff, 1024, x->skip_block, p->zbin,
                                    p->round, p->quant, p->quant_shift, qcoeff,
                                    dqcoeff, pd->dequant, eob,
                                    scan_order->scan, scan_order->iscan);
        break;
      case TX_16X16:
#if CONFIG_EXT_TX
        highbd_forw_tx16x16(x, plane, src_diff, diff_stride, coeff);
#else
        vp9_highbd_fdct16x16(src_diff, coeff, diff_stride);
#endif
        vp9_highbd_quantize_b(coeff, 256, x->skip_block, p->zbin, p->round,
                              p->quant, p->quant_shift, qcoeff, dqcoeff,
                              pd->dequant, eob,
                              scan_order->scan, scan_order->iscan);
        break;
      case TX_8X8:
#if CONFIG_EXT_TX
        highbd_forw_tx8x8(x, plane, src_diff, diff_stride, coeff);
#else
        vp9_highbd_fdct8x8(src_diff, coeff, diff_stride);
#endif
        vp9_highbd_quantize_b(coeff, 64, x->skip_block, p->zbin, p->round,
                              p->quant, p->quant_shift, qcoeff, dqcoeff,
                              pd->dequant, eob,
                              scan_order->scan, scan_order->iscan);
        break;
      case TX_4X4:
#if CONFIG_EXT_TX
        highbd_forw_tx4x4(x, plane, block, src_diff, diff_stride, coeff);
#else
        x->fwd_txm4x4(src_diff, coeff, diff_stride);
#endif
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
#if CONFIG_TX64X64
    case TX_64X64:
      vp9_fdct64x64(src_diff, coeff, diff_stride);
      vp9_quantize_b_64x64(coeff, 4096, x->skip_block, p->zbin, p->round,
                           p->quant, p->quant_shift, qcoeff, dqcoeff,
                           pd->dequant, eob, scan_order->scan,
                           scan_order->iscan);
      break;
#endif  // CONFIG_TX64X64
    case TX_32X32:
      fdct32x32(x->use_lp32x32fdct, src_diff, coeff, diff_stride);
      vp9_quantize_b_32x32(coeff, 1024, x->skip_block, p->zbin, p->round,
                           p->quant, p->quant_shift, qcoeff, dqcoeff,
                           pd->dequant, eob, scan_order->scan,
                           scan_order->iscan);
      break;
    case TX_16X16:
#if CONFIG_EXT_TX
      forw_tx16x16(x, plane, src_diff, diff_stride, coeff);
#else
      vp9_fdct16x16(src_diff, coeff, diff_stride);
#endif
      vp9_quantize_b(coeff, 256, x->skip_block, p->zbin, p->round,
                     p->quant, p->quant_shift, qcoeff, dqcoeff,
                     pd->dequant, eob,
                     scan_order->scan, scan_order->iscan);
      break;
    case TX_8X8:
#if CONFIG_EXT_TX
      forw_tx8x8(x, plane, src_diff, diff_stride, coeff);
#else
      vp9_fdct8x8(src_diff, coeff, diff_stride);
#endif
      vp9_quantize_b(coeff, 64, x->skip_block, p->zbin, p->round,
                     p->quant, p->quant_shift, qcoeff, dqcoeff,
                     pd->dequant, eob,
                     scan_order->scan, scan_order->iscan);
      break;
    case TX_4X4:
#if CONFIG_EXT_TX
      forw_tx4x4(x, plane, block, src_diff, diff_stride, coeff);
#else
      x->fwd_txm4x4(src_diff, coeff, diff_stride);
#endif
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

static void encode_block(int plane, int block, BLOCK_SIZE plane_bsize,
                         TX_SIZE tx_size, void *arg) {
  struct encode_b_args *const args = arg;
  MACROBLOCK *const x = args->x;
  MACROBLOCKD *const xd = &x->e_mbd;
  struct optimize_ctx *const ctx = args->ctx;
  struct macroblock_plane *const p = &x->plane[plane];
  struct macroblockd_plane *const pd = &xd->plane[plane];
  tran_low_t *const dqcoeff = BLOCK_OFFSET(pd->dqcoeff, block);
  int i, j;
  uint8_t *dst;
  ENTROPY_CONTEXT *a, *l;
#if CONFIG_EXT_TX
  TX_TYPE tx_type;
#endif
#if CONFIG_TX_SKIP
  MB_MODE_INFO *mbmi = &xd->mi[0].src_mi->mbmi;
  int shift = mbmi->tx_skip_shift;
#endif
  txfrm_block_to_raster_xy(plane_bsize, tx_size, block, &i, &j);
  dst = &pd->dst.buf[4 * j * pd->dst.stride + 4 * i];
  a = &ctx->ta[plane][i];
  l = &ctx->tl[plane][j];
#if CONFIG_TX64X64
  if (plane) assert(tx_size != TX_64X64);
#endif

  // TODO(jingning): per transformed block zero forcing only enabled for
  // luma component. will integrate chroma components as well.
  if (plane == 0 && x->zcoeff_blk[tx_size][block]) {
    p->eobs[block] = 0;
    *a = *l = 0;
    return;
  }

  if (!x->skip_recode) {
    if (max_txsize_lookup[plane_bsize] == tx_size) {
      if (x->skip_txfm[(plane << 2) + (block >> (tx_size << 1))] == 0) {
        // full forward transform and quantization
#if CONFIG_NEW_QUANT
        if (x->quant_fp)
          vp9_xform_quant_fp_nuq(x, plane, block, plane_bsize, tx_size);
        else
          vp9_xform_quant_nuq(x, plane, block, plane_bsize, tx_size);
#else
        if (x->quant_fp)
          vp9_xform_quant_fp(x, plane, block, plane_bsize, tx_size);
        else
          vp9_xform_quant(x, plane, block, plane_bsize, tx_size);
#endif  // CONFIG_NEW_QUANT
      } else if (x->skip_txfm[(plane << 2) + (block >> (tx_size << 1))] == 2) {
        // fast path forward transform and quantization
#if CONFIG_NEW_QUANT
        if (x->quant_fp)
          vp9_xform_quant_dc_fp_nuq(x, plane, block, plane_bsize, tx_size);
        else
          vp9_xform_quant_dc_nuq(x, plane, block, plane_bsize, tx_size);
#else
        vp9_xform_quant_dc(x, plane, block, plane_bsize, tx_size);
#endif
      } else {
        // skip forward transform
        p->eobs[block] = 0;
        *a = *l = 0;
        return;
      }
    } else {
#if CONFIG_NEW_QUANT
      if (x->quant_fp)
        vp9_xform_quant_fp_nuq(x, plane, block, plane_bsize, tx_size);
      else
        vp9_xform_quant_nuq(x, plane, block, plane_bsize, tx_size);
#else
      if (x->quant_fp)
        vp9_xform_quant_fp(x, plane, block, plane_bsize, tx_size);
      else
        vp9_xform_quant(x, plane, block, plane_bsize, tx_size);
#endif
    }
  }

  if (x->optimize && (!x->skip_recode || !x->skip_optimize)) {
    const int ctx = combine_entropy_contexts(*a, *l);
    *a = *l = optimize_b(x, plane, block, tx_size, ctx) > 0;
  } else {
    *a = *l = p->eobs[block] > 0;
  }

  if (p->eobs[block])
    *(args->skip) = 0;

  if (x->skip_encode || p->eobs[block] == 0)
    return;

#if CONFIG_TX_SKIP
  if (mbmi->tx_skip[plane != 0]) {
    int bs = 4 << tx_size;
#if CONFIG_VP9_HIGHBITDEPTH
    if (xd->cur_buf->flags & YV12_FLAG_HIGHBITDEPTH)
      vp9_highbd_tx_identity_add(dqcoeff, dst, pd->dst.stride, bs, shift,
                                 xd->bd);
    else
#endif  // CONFIG_VP9_HIGHBITDEPTH
      vp9_tx_identity_add(dqcoeff, dst, pd->dst.stride, bs, shift);

    return;
  }
#endif  // CONFIG_TX_SKIP

#if CONFIG_VP9_HIGHBITDEPTH
  if (xd->cur_buf->flags & YV12_FLAG_HIGHBITDEPTH) {
    switch (tx_size) {
#if CONFIG_TX64X64
      case TX_64X64:
        vp9_highbd_idct64x64_add(dqcoeff, dst, pd->dst.stride,
                                 p->eobs[block], xd->bd);
        break;
#endif  // CONFIG_TX64X64
      case TX_32X32:
        vp9_highbd_idct32x32_add(dqcoeff, dst, pd->dst.stride,
                                 p->eobs[block], xd->bd);
        break;
      case TX_16X16:
#if CONFIG_EXT_TX
        tx_type = get_tx_type(plane, xd);
        vp9_highbd_iht16x16_add(tx_type, dqcoeff, dst, pd->dst.stride,
                                p->eobs[block], xd->bd);
#else
        vp9_highbd_idct16x16_add(dqcoeff, dst, pd->dst.stride,
                                 p->eobs[block], xd->bd);
#endif
        break;
      case TX_8X8:
#if CONFIG_EXT_TX
        tx_type = get_tx_type(plane, xd);
        vp9_highbd_iht8x8_add(tx_type, dqcoeff, dst, pd->dst.stride,
                              p->eobs[block], xd->bd);
#else
        vp9_highbd_idct8x8_add(dqcoeff, dst, pd->dst.stride,
                               p->eobs[block], xd->bd);
#endif
        break;
      case TX_4X4:
#if CONFIG_EXT_TX
        tx_type = get_tx_type_4x4(plane, xd, block);
        if (tx_type == DCT_DCT) {
          // This is like vp9_short_idct4x4 but has a special case around eob<=1
          // which is significant (not just an optimization) for the lossless
          // case.
          x->highbd_itxm_add(dqcoeff, dst, pd->dst.stride,
                             p->eobs[block], xd->bd);
        } else {
          vp9_highbd_iht4x4_add(tx_type, dqcoeff, dst, pd->dst.stride,
                                p->eobs[block], xd->bd);
        }
#else
        // This is like vp9_short_idct4x4 but has a special case around eob<=1
        // which is significant (not just an optimization) for the lossless
        // case.
        x->highbd_itxm_add(dqcoeff, dst, pd->dst.stride,
                           p->eobs[block], xd->bd);
#endif
        break;
      default:
        assert(0 && "Invalid transform size");
    }
    return;
  }
#endif  // CONFIG_VP9_HIGHBITDEPTH

  switch (tx_size) {
#if CONFIG_TX64X64
    case TX_64X64:
      vp9_idct64x64_add(dqcoeff, dst, pd->dst.stride, p->eobs[block]);
      break;
#endif
    case TX_32X32:
      vp9_idct32x32_add(dqcoeff, dst, pd->dst.stride, p->eobs[block]);
      break;
    case TX_16X16:
#if CONFIG_EXT_TX
      tx_type = get_tx_type(plane, xd);
      vp9_iht16x16_add(tx_type, dqcoeff, dst, pd->dst.stride, p->eobs[block]);
#else
      vp9_idct16x16_add(dqcoeff, dst, pd->dst.stride, p->eobs[block]);
#endif
      break;
    case TX_8X8:
#if CONFIG_EXT_TX
      tx_type = get_tx_type(plane, xd);
      vp9_iht8x8_add(tx_type, dqcoeff, dst, pd->dst.stride, p->eobs[block]);
#else
      vp9_idct8x8_add(dqcoeff, dst, pd->dst.stride, p->eobs[block]);
#endif
      break;
    case TX_4X4:
#if CONFIG_EXT_TX
      tx_type = get_tx_type_4x4(plane, xd, block);
      if (tx_type == DCT_DCT) {
        // This is like vp9_short_idct4x4 but has a special case around eob<=1
        // which is significant (not just an optimization) for the lossless
        // case.
        x->itxm_add(dqcoeff, dst, pd->dst.stride, p->eobs[block]);
      } else {
        vp9_iht4x4_add(tx_type, dqcoeff, dst, pd->dst.stride, p->eobs[block]);
      }
#else
      // This is like vp9_short_idct4x4 but has a special case around eob<=1
      // which is significant (not just an optimization) for the lossless case.
      x->itxm_add(dqcoeff, dst, pd->dst.stride, p->eobs[block]);
#endif
      break;
    default:
      assert(0 && "Invalid transform size");
      break;
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
#if CONFIG_EXT_TX
  MB_MODE_INFO *mbmi = &xd->mi[0].src_mi->mbmi;
  mbmi->ext_txfrm = NORM;
#endif
#if CONFIG_TX_SKIP
  xd->mi[0].src_mi->mbmi.tx_skip[0] = 0;
  xd->mi[0].src_mi->mbmi.tx_skip[1] = 0;
#endif
  txfrm_block_to_raster_xy(plane_bsize, tx_size, block, &i, &j);
  dst = &pd->dst.buf[4 * j * pd->dst.stride + 4 * i];

#if CONFIG_NEW_QUANT
  if (x->quant_fp)
    vp9_xform_quant_fp_nuq(x, plane, block, plane_bsize, tx_size);
  else
    vp9_xform_quant_nuq(x, plane, block, plane_bsize, tx_size);
#else
  if (x->quant_fp)
    vp9_xform_quant_fp(x, plane, block, plane_bsize, tx_size);
  else
    vp9_xform_quant(x, plane, block, plane_bsize, tx_size);
#endif

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
    if (!x->skip_recode)
      vp9_subtract_plane(x, bsize, plane);

    if (x->optimize && (!x->skip_recode || !x->skip_optimize)) {
      const struct macroblockd_plane* const pd = &xd->plane[plane];
      const TX_SIZE tx_size = plane ? get_uv_tx_size(mbmi, pd) : mbmi->tx_size;
      vp9_get_entropy_contexts(bsize, tx_size, pd,
                               ctx.ta[plane], ctx.tl[plane]);
    }

    vp9_foreach_transformed_block_in_plane(xd, bsize, plane, encode_block,
                                           &arg);
  }
}

#if CONFIG_SUPERTX
void vp9_encode_sb_supertx(MACROBLOCK *x, BLOCK_SIZE bsize) {
  MACROBLOCKD *const xd = &x->e_mbd;
  struct optimize_ctx ctx;
  MB_MODE_INFO *mbmi = &xd->mi[0].src_mi->mbmi;
  struct encode_b_args arg = {x, &ctx, &mbmi->skip};
  int plane;

  mbmi->skip = 1;
  if (x->skip)
    return;

  for (plane = 0; plane < MAX_MB_PLANE; ++plane) {
    const struct macroblockd_plane* const pd = &xd->plane[plane];
    const BLOCK_SIZE plane_size = get_plane_block_size(bsize, pd);
    const TX_SIZE tx_size = plane ? get_uv_tx_size(mbmi, pd) : mbmi->tx_size;
    vp9_subtract_plane(x, bsize, plane);
    vp9_get_entropy_contexts(bsize, tx_size, pd,
                             ctx.ta[plane], ctx.tl[plane]);
    encode_block(plane, 0, plane_size, tx_size, &arg);
  }
}
#endif

#if CONFIG_TX_SKIP
static int vp9_dpcm_intra(uint8_t *src, int src_stride,
                          uint8_t *dst, int dst_stride,
                          int16_t *src_diff, int diff_stride,
                          tran_low_t *coeff, tran_low_t *qcoeff,
                          tran_low_t *dqcoeff, struct macroblock_plane *p,
                          struct macroblockd_plane *pd,
                          const scan_order *scan_order, PREDICTION_MODE mode,
                          TX_SIZE tx_size, int shift, int logsizeby32) {
  int i, j, eob, temp;
  const int bs = 4 << tx_size;

  vpx_memset(qcoeff, 0, bs * bs * sizeof(*qcoeff));
  vpx_memset(dqcoeff, 0, bs * bs * sizeof(*dqcoeff));

  switch (mode) {
    case H_PRED:
      for (i = 0 ; i < bs; i++) {
        vp9_subtract_block_c(bs, 1, src_diff + i, diff_stride,
                             src + i, src_stride,
                             dst + i, dst_stride);
        vp9_tx_identity_rect(src_diff + i, coeff + i, bs, 1,
                             diff_stride, bs, shift);
        vp9_quantize_rect(coeff + i, bs, 1, p->zbin_pxd, p->round_pxd,
                          p->quant_pxd,
                          p->quant_shift_pxd, qcoeff + i, dqcoeff + i,
                          pd->dequant_pxd, logsizeby32, bs, i == 0, 0);
        vp9_tx_identity_add_rect(dqcoeff + i, dst + i, bs, 1,
                                 bs, dst_stride, shift);
        if (i < bs - 1)
          for (j = 0 ; j < bs; j++)
            *(dst + j * dst_stride + i + 1) =
                *(dst + j * dst_stride + i);
      }
      break;
    case V_PRED:
      for (i = 0 ; i < bs; i++) {
        vp9_subtract_block_c(1, bs, src_diff + diff_stride * i,
                             diff_stride,
                             src + src_stride * i, src_stride,
                             dst + dst_stride * i, dst_stride);
        vp9_tx_identity_rect(src_diff + diff_stride * i,
                             coeff + bs * i, 1, bs,
                             diff_stride, bs, shift);
        vp9_quantize_rect(coeff + bs * i, 1, bs, p->zbin_pxd, p->round_pxd,
                          p->quant_pxd,
                          p->quant_shift_pxd, qcoeff + bs * i, dqcoeff + bs * i,
                          pd->dequant_pxd, logsizeby32, bs, i == 0, 0);
        vp9_tx_identity_add_rect(dqcoeff + bs * i, dst + dst_stride * i,
                                 1, bs, bs, dst_stride, shift);
        if (i < bs - 1)
          vpx_memcpy(dst + (i + 1) * dst_stride,
                     dst + i * dst_stride, bs * sizeof(dst[0]));
      }
      break;
    case TM_PRED:
      vp9_subtract_block_c(1, bs, src_diff, diff_stride, src, src_stride,
                           dst, dst_stride);
      vp9_tx_identity_rect(src_diff, coeff, 1, bs, diff_stride, bs, shift);
      vp9_quantize_rect(coeff, 1, bs, p->zbin_pxd, p->round_pxd, p->quant_pxd,
                        p->quant_shift_pxd, qcoeff, dqcoeff, pd->dequant_pxd,
                        logsizeby32, bs, 1, 0);
      vp9_tx_identity_add_rect(dqcoeff, dst, 1, bs, bs, dst_stride, shift);

      vp9_subtract_block_c(bs -1, 1, src_diff + diff_stride, diff_stride,
                           src + src_stride, src_stride,
                           dst + dst_stride, dst_stride);
      vp9_tx_identity_rect(src_diff + diff_stride, coeff + bs, bs - 1, 1,
                           diff_stride, bs, shift);
      vp9_quantize_rect(coeff + bs, bs - 1, 1, p->zbin_pxd, p->round_pxd,
                        p->quant_pxd,
                        p->quant_shift_pxd, qcoeff + bs, dqcoeff + bs,
                        pd->dequant_pxd, logsizeby32, bs, 0, 0);
      vp9_tx_identity_add_rect(dqcoeff + bs, dst + dst_stride, bs - 1, 1,
                               bs, dst_stride, shift);

      for (i = 1 ; i < bs; i++) {
        for (j = 1 ; j < bs; j++) {
          temp = dst[(i - 1) * dst_stride + j] + dst[i * dst_stride + j - 1] -
                 dst[(i - 1) * dst_stride + j - 1];
          temp = clip_pixel(temp);
          dst[i * dst_stride + j] = temp;
          vp9_subtract_block_c(1, 1, src_diff + diff_stride * i + j,
                               diff_stride, src + src_stride * i + j,
                               src_stride, dst + dst_stride * i + j,
                               dst_stride);
          vp9_tx_identity_rect(src_diff + i * diff_stride + j,
                               coeff + bs * i + j, 1, 1, diff_stride,
                               bs, shift);
          vp9_quantize_rect(coeff + bs * i + j, 1, 1, p->zbin_pxd, p->round_pxd,
                            p->quant_pxd, p->quant_shift_pxd,
                            qcoeff + bs * i + j, dqcoeff + bs * i + j,
                            pd->dequant_pxd, logsizeby32, bs, 0, 0);
          vp9_tx_identity_add_rect(dqcoeff + bs * i + j,
                                   dst + dst_stride * i + j, 1, 1, bs,
                                   dst_stride, shift);
        }
      }
      break;
    default:
      break;
  }

  eob = get_eob(qcoeff, bs * bs, scan_order->scan);
  return eob;
}

#if CONFIG_VP9_HIGHBITDEPTH
static int vp9_highbd_dpcm_intra(uint8_t *src, int src_stride,
                                 uint8_t *dst, int dst_stride,
                                 int16_t *src_diff, int diff_stride,
                                 tran_low_t *coeff, tran_low_t *qcoeff,
                                 tran_low_t *dqcoeff,
                                 struct macroblock_plane *p,
                                 struct macroblockd_plane *pd,
                                 const scan_order *scan_order,
                                 PREDICTION_MODE mode, TX_SIZE tx_size,
                                 int shift, int logsizeby32, int bd) {
  int i, j, eob, temp;
  const int bs = 4 << tx_size;
  uint16_t *dst16 = CONVERT_TO_SHORTPTR(dst);

  vpx_memset(qcoeff, 0, bs * bs * sizeof(*qcoeff));
  vpx_memset(dqcoeff, 0, bs * bs * sizeof(*dqcoeff));

  switch (mode) {
    case H_PRED:
      for (i = 0 ; i < bs; i++) {
        vp9_highbd_subtract_block_c(bs, 1, src_diff + i, diff_stride,
                                    src + i, src_stride, dst + i,
                                    dst_stride, bd);
        vp9_tx_identity_rect(src_diff + i, coeff + i, bs, 1,
                             diff_stride, bs, shift);
        vp9_quantize_rect(coeff + i, bs, 1, p->zbin_pxd, p->round_pxd,
                          p->quant_pxd, p->quant_shift_pxd, qcoeff + i,
                          dqcoeff + i, pd->dequant_pxd, logsizeby32, bs,
                          i == 0, 1);
        vp9_highbd_tx_identity_add_rect(dqcoeff + i, dst + i, bs, 1,
                                        bs, dst_stride, shift, bd);
        if (i < bs - 1)
          for (j = 0 ; j < bs; j++)
            *(dst16 + j * dst_stride + i + 1) =
                *(dst16 + j * dst_stride + i);
      }
      break;
    case V_PRED:
      for (i = 0 ; i < bs; i++) {
        vp9_highbd_subtract_block_c(1, bs, src_diff + diff_stride * i,
                                    diff_stride, src + src_stride * i,
                                    src_stride,  dst + dst_stride * i,
                                    dst_stride, bd);
        vp9_tx_identity_rect(src_diff + diff_stride * i, coeff + bs * i, 1, bs,
                             diff_stride, bs, shift);
        vp9_quantize_rect(coeff + bs * i, 1, bs, p->zbin_pxd, p->round_pxd,
                          p->quant_pxd, p->quant_shift_pxd, qcoeff + bs * i,
                          dqcoeff + bs * i, pd->dequant_pxd, logsizeby32, bs,
                          i == 0, 1);
        vp9_highbd_tx_identity_add_rect(dqcoeff + bs * i, dst + dst_stride * i,
                                        1, bs, bs, dst_stride, shift, bd);
        if (i < bs - 1)
          vpx_memcpy(dst16 + (i + 1) * dst_stride,
                     dst16 + i * dst_stride, bs * sizeof(dst16[0]));
      }
      break;
    case TM_PRED:
      vp9_highbd_subtract_block_c(1, bs, src_diff, diff_stride, src, src_stride,
                                  dst, dst_stride, bd);
      vp9_tx_identity_rect(src_diff, coeff, 1, bs, diff_stride, bs, shift);
      vp9_quantize_rect(coeff, 1, bs, p->zbin_pxd, p->round_pxd, p->quant_pxd,
                        p->quant_shift_pxd, qcoeff, dqcoeff, pd->dequant_pxd,
                        logsizeby32, bs, 1, 1);
      vp9_highbd_tx_identity_add_rect(dqcoeff, dst, 1, bs, bs, dst_stride,
                                      shift, bd);
      vp9_highbd_subtract_block_c(bs -1, 1, src_diff + diff_stride, diff_stride,
                                  src + src_stride, src_stride,
                                  dst + dst_stride, dst_stride, bd);
      vp9_tx_identity_rect(src_diff + diff_stride, coeff + bs, bs - 1, 1,
                           diff_stride, bs, shift);
      vp9_quantize_rect(coeff + bs, bs - 1, 1, p->zbin_pxd, p->round_pxd,
                        p->quant_pxd, p->quant_shift_pxd, qcoeff + bs,
                        dqcoeff + bs, pd->dequant_pxd, logsizeby32, bs, 0, 1);
      vp9_highbd_tx_identity_add_rect(dqcoeff + bs, dst + dst_stride, bs - 1, 1,
                                      bs, dst_stride, shift, bd);

      for (i = 1 ; i < bs; i++) {
        for (j = 1 ; j < bs; j++) {
          temp = dst16[(i - 1) * dst_stride + j] +
              dst16[i * dst_stride + j - 1] -
              dst16[(i - 1) * dst_stride + j - 1];
          dst16[i * dst_stride + j] = clip_pixel_highbd(temp, bd);
          vp9_highbd_subtract_block_c(1, 1, src_diff + diff_stride * i + j,
                                      diff_stride, src + src_stride * i + j,
                                      src_stride, dst + dst_stride * i + j,
                                      dst_stride, bd);
          vp9_tx_identity_rect(src_diff + i * diff_stride + j,
                               coeff + bs * i + j, 1, 1, diff_stride,
                               bs, shift);
          vp9_quantize_rect(coeff + bs * i + j, 1, 1, p->zbin_pxd, p->round_pxd,
                            p->quant_pxd, p->quant_shift_pxd,
                            qcoeff + bs * i + j, dqcoeff + bs * i + j,
                            pd->dequant_pxd,
                            logsizeby32, bs, 0, 1);
          vp9_highbd_tx_identity_add_rect(dqcoeff + bs * i + j,
                                          dst + dst_stride * i + j, 1, 1, bs,
                                          dst_stride, shift, bd);
        }
      }
      break;
    default:
      break;
  }

  eob = get_eob(qcoeff, bs * bs, scan_order->scan);
  return eob;
}
#endif  // CONFIG_VP9_HIGHBITDEPTH
#endif  // CONFIG_TX_SKIP

static void encode_block_intra(int plane, int block, BLOCK_SIZE plane_bsize,
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
#if CONFIG_FILTERINTRA
  int fbit = 0;
#endif  // CONFIG_FILTERINTRA
  const int bwl = b_width_log2_lookup[plane_bsize];
  const int diff_stride = 4 * (1 << bwl);
  uint8_t *src, *dst;
  int16_t *src_diff;
  uint16_t *eob = &p->eobs[block];
  const int src_stride = p->src.stride;
  const int dst_stride = pd->dst.stride;
  int i, j;
#if CONFIG_NEW_QUANT
  const uint8_t* band = get_band_translate(tx_size);
#endif  // CONFIG_NEW_QUANT
  txfrm_block_to_raster_xy(plane_bsize, tx_size, block, &i, &j);
  dst = &pd->dst.buf[4 * (j * dst_stride + i)];
  src = &p->src.buf[4 * (j * src_stride + i)];
  src_diff = &p->src_diff[4 * (j * diff_stride + i)];

#if CONFIG_FILTERINTRA
  if (mbmi->sb_type < BLOCK_8X8 && plane == 0)
    fbit = xd->mi[0].b_filter_info[block];
  else
    fbit = plane == 0 ? mbmi->filterbit : mbmi->uv_filterbit;
#endif  // CONFIG_FILTERINTRA

#if CONFIG_TX_SKIP
  if (mbmi->tx_skip[plane != 0]) {
    int shift = mbmi->tx_skip_shift;
    int bs = 4 << tx_size;
#if CONFIG_VP9_HIGHBITDEPTH
    int use_hbd = xd->cur_buf->flags & YV12_FLAG_HIGHBITDEPTH;
#endif  // CONFIG_VP9_HIGHBITDEPTH
#if CONFIG_NEW_QUANT
    band = vp9_coefband_tx_skip;
#endif  // CONFIG_NEW_QUANT
    mode = plane == 0 ? mbmi->mode : mbmi->uv_mode;
    scan_order = &vp9_default_scan_orders_pxd[tx_size];

    vp9_predict_intra_block(xd, block >> (2 * tx_size), bwl, tx_size, mode,
#if CONFIG_FILTERINTRA
                            fbit,
#endif
                            x->skip_encode ? src : dst,
                                x->skip_encode ? src_stride : dst_stride,
                                    dst, dst_stride, i, j, plane);

    if (!x->skip_recode && tx_size <= TX_32X32 &&
        (mode == H_PRED || mode == V_PRED || mode == TM_PRED)) {
#if CONFIG_VP9_HIGHBITDEPTH
      if (use_hbd)
        *eob = vp9_highbd_dpcm_intra(src, src_stride, dst, dst_stride,
                                     src_diff, diff_stride,
                                     coeff, qcoeff, dqcoeff, p, pd,
                                     scan_order, mode, tx_size, shift,
                                     tx_size > TX_16X16 ? 0 : -1, xd->bd);
      else
#endif  // CONFIG_VP9_HIGHBITDEPTH
        *eob = vp9_dpcm_intra(src, src_stride, dst, dst_stride,
                              src_diff, diff_stride,
                              coeff, qcoeff, dqcoeff, p, pd,
                              scan_order, mode, tx_size, shift,
                              tx_size > TX_16X16 ? 0 : -1);

      if (*eob)
        *(args->skip) = 0;
      return;
    }


    if (!x->skip_recode) {
#if CONFIG_VP9_HIGHBITDEPTH
      if (use_hbd) {
        vp9_highbd_subtract_block(bs, bs, src_diff, diff_stride,
                                  src, src_stride, dst, dst_stride, xd->bd);
        vp9_tx_identity(src_diff, coeff, diff_stride, bs, shift);
      } else {
        vp9_subtract_block(bs, bs, src_diff, diff_stride,
                           src, src_stride, dst, dst_stride);
        vp9_tx_identity(src_diff, coeff, diff_stride, bs, shift);
      }
#else
      vp9_subtract_block(bs, bs, src_diff, diff_stride,
                         src, src_stride, dst, dst_stride);
      vp9_tx_identity(src_diff, coeff, diff_stride, bs, shift);
#endif  // CONFIG_VP9_HIGHBITDEPTH

      if (tx_size <= TX_16X16) {
#if CONFIG_NEW_QUANT
#if CONFIG_VP9_HIGHBITDEPTH
        if (use_hbd) {
          if (x->quant_fp)
            vp9_highbd_quantize_fp_nuq(coeff, bs * bs, x->skip_block,
                                       p->quant_pxd_fp, pd->dequant_pxd,
                                       (const cumbins_type_nuq *)
                                       p->cumbins_nuq_pxd,
                                       (const dequant_val_type_nuq *)
                                       pd->dequant_val_nuq_pxd,
                                       qcoeff, dqcoeff, eob,
                                       scan_order->scan, band);
          else
            vp9_highbd_quantize_nuq(coeff, bs * bs, x->skip_block,
                                    p->quant_pxd, p->quant_shift_pxd,
                                    pd->dequant_pxd,
                                    (const cumbins_type_nuq *)
                                    p->cumbins_nuq_pxd,
                                    (const dequant_val_type_nuq *)
                                    pd->dequant_val_nuq_pxd,
                                    qcoeff, dqcoeff, eob,
                                    scan_order->scan, band);
        } else {
          if (x->quant_fp)
            vp9_quantize_fp_nuq(coeff, bs * bs, x->skip_block,
                                p->quant_pxd_fp, pd->dequant_pxd,
                                (const cumbins_type_nuq *)p->cumbins_nuq_pxd,
                                (const dequant_val_type_nuq *)
                                pd->dequant_val_nuq_pxd,
                                qcoeff, dqcoeff, eob,
                                scan_order->scan, band);
          else
            vp9_quantize_nuq(coeff, bs * bs, x->skip_block,
                             p->quant_pxd, p->quant_shift_pxd, pd->dequant_pxd,
                             (const cumbins_type_nuq *)p->cumbins_nuq_pxd,
                             (const dequant_val_type_nuq *)
                             pd->dequant_val_nuq_pxd,
                             qcoeff, dqcoeff, eob,
                             scan_order->scan, band);
        }
#else  // CONFIG_VP9_HIGHBITDEPTH
        if (x->quant_fp)
          vp9_quantize_fp_nuq(coeff, bs * bs, x->skip_block,
                              p->quant_pxd_fp, pd->dequant_pxd,
                              (const cumbins_type_nuq *)p->cumbins_nuq_pxd,
                              (const dequant_val_type_nuq *)
                              pd->dequant_val_nuq_pxd,
                              qcoeff, dqcoeff, eob,
                              scan_order->scan, band);
        else
          vp9_quantize_nuq(coeff, bs * bs, x->skip_block,
                           p->quant_pxd, p->quant_shift_pxd, pd->dequant_pxd,
                           (const cumbins_type_nuq *)p->cumbins_nuq_pxd,
                           (const dequant_val_type_nuq *)
                           pd->dequant_val_nuq_pxd,
                           qcoeff, dqcoeff, eob,
                           scan_order->scan, band);
#endif  // CONFIG_VP9_HIGHBITDEPTH
#else  // CONFIG_NEW_QUANT
#if CONFIG_VP9_HIGHBITDEPTH
        if (use_hbd)
          vp9_highbd_quantize_b(coeff, bs * bs, x->skip_block, p->zbin_pxd,
                                p->round_pxd, p->quant_pxd, p->quant_shift_pxd,
                                qcoeff, dqcoeff, pd->dequant_pxd, eob,
                                scan_order->scan, scan_order->iscan);
        else
#endif  // CONFIG_VP9_HIGHBITDEPTH
          vp9_quantize_b(coeff, bs * bs, x->skip_block, p->zbin_pxd,
                         p->round_pxd, p->quant_pxd, p->quant_shift_pxd,
                         qcoeff, dqcoeff, pd->dequant_pxd, eob,
                         scan_order->scan, scan_order->iscan);
#endif  // CONFIG_NEW_QUANT
      } else if (tx_size == TX_32X32) {
#if CONFIG_NEW_QUANT
#if CONFIG_VP9_HIGHBITDEPTH
        if (use_hbd) {
          if (x->quant_fp)
            vp9_highbd_quantize_32x32_fp_nuq(coeff, bs * bs, x->skip_block,
                                             p->quant_pxd_fp, pd->dequant_pxd,
                                             (const cumbins_type_nuq *)p->
                                             cumbins_nuq_pxd,
                                             (const dequant_val_type_nuq *)
                                             pd->dequant_val_nuq_pxd,
                                             qcoeff, dqcoeff, eob,
                                             scan_order->scan, band);
          else
            vp9_highbd_quantize_32x32_nuq(coeff, bs * bs, x->skip_block,
                                          p->quant_pxd, p->quant_shift_pxd,
                                          pd->dequant_pxd,
                                          (const cumbins_type_nuq *)
                                          p->cumbins_nuq_pxd,
                                          (const dequant_val_type_nuq *)
                                          pd->dequant_val_nuq_pxd,
                                          qcoeff, dqcoeff, eob,
                                          scan_order->scan, band);
        } else {
          if (x->quant_fp)
            vp9_quantize_32x32_fp_nuq(coeff, bs * bs, x->skip_block,
                                      p->quant_pxd_fp, pd->dequant_pxd,
                                      (const cumbins_type_nuq *)
                                      p->cumbins_nuq_pxd,
                                      (const dequant_val_type_nuq *)
                                      pd->dequant_val_nuq_pxd,
                                      qcoeff, dqcoeff, eob,
                                      scan_order->scan, band);
          else
            vp9_quantize_32x32_nuq(coeff, bs * bs, x->skip_block,
                                   p->quant_pxd, p->quant_shift_pxd,
                                   pd->dequant_pxd,
                                   (const cumbins_type_nuq *)p->cumbins_nuq_pxd,
                                   (const dequant_val_type_nuq *)
                                   pd->dequant_val_nuq_pxd,
                                   qcoeff, dqcoeff, eob,
                                   scan_order->scan, band);
        }
#else  // CONFIG_VP9_HIGHBITDEPTH
        if (x->quant_fp)
          vp9_quantize_32x32_fp_nuq(coeff, bs * bs, x->skip_block,
                                    p->quant_pxd_fp, pd->dequant_pxd,
                                    (const cumbins_type_nuq *)
                                    p->cumbins_nuq_pxd,
                                    (const dequant_val_type_nuq *)
                                    pd->dequant_val_nuq_pxd,
                                    qcoeff, dqcoeff, eob,
                                    scan_order->scan, band);
        else
          vp9_quantize_32x32_nuq(coeff, bs * bs, x->skip_block,
                                 p->quant_pxd, p->quant_shift_pxd,
                                 pd->dequant_pxd,
                                 (const cumbins_type_nuq *)p->cumbins_nuq_pxd,
                                 (const dequant_val_type_nuq *)
                                 pd->dequant_val_nuq_pxd,
                                 qcoeff, dqcoeff, eob,
                                 scan_order->scan, band);/**/
#endif  // CONFIG_VP9_HIGHBITDEPTH
#else  // CONFIG_NEW_QUANT
#if CONFIG_VP9_HIGHBITDEPTH
        if (use_hbd)
          vp9_highbd_quantize_b_32x32(coeff, bs * bs, x->skip_block,
                                      p->zbin_pxd, p->round_pxd, p->quant_pxd,
                                      p->quant_shift_pxd, qcoeff, dqcoeff,
                                      pd->dequant_pxd, eob,
                                      scan_order->scan, scan_order->iscan);
        else
#endif  // CONFIG_VP9_HIGHBITDEPTH
          vp9_quantize_b_32x32(coeff, bs * bs, x->skip_block, p->zbin_pxd,
                               p->round_pxd, p->quant_pxd, p->quant_shift_pxd,
                               qcoeff, dqcoeff, pd->dequant_pxd, eob,
                               scan_order->scan, scan_order->iscan);
#endif  // CONFIG_NEW_QUANT
      }
#if CONFIG_TX64X64
      else if (tx_size == TX_64X64) {
#if CONFIG_NEW_QUANT
#if CONFIG_VP9_HIGHBITDEPTH
        if (use_hbd) {
          if (x->quant_fp)
            vp9_highbd_quantize_64x64_fp_nuq(coeff, bs * bs, x->skip_block,
                                             p->quant_pxd_fp, pd->dequant_pxd,
                                             (const cumbins_type_nuq *)
                                             p->cumbins_nuq_pxd,
                                             (const dequant_val_type_nuq *)
                                             pd->dequant_val_nuq_pxd,
                                             qcoeff, dqcoeff, eob,
                                             scan_order->scan, band);
          else
            vp9_highbd_quantize_64x64_nuq(coeff, bs * bs, x->skip_block,
                                          p->quant_pxd, p->quant_shift_pxd,
                                          pd->dequant_pxd,
                                          (const cumbins_type_nuq *)
                                          p->cumbins_nuq_pxd,
                                          (const dequant_val_type_nuq *)
                                          pd->dequant_val_nuq_pxd,
                                          qcoeff, dqcoeff, eob,
                                          scan_order->scan, band);
        } else {
          if (x->quant_fp)
            vp9_quantize_64x64_fp_nuq(coeff, bs * bs, x->skip_block,
                                      p->quant_pxd_fp, pd->dequant_pxd,
                                      (const cumbins_type_nuq *)
                                      p->cumbins_nuq_pxd,
                                      (const dequant_val_type_nuq *)
                                      pd->dequant_val_nuq_pxd,
                                      qcoeff, dqcoeff, eob,
                                      scan_order->scan, band);
          else
            vp9_quantize_64x64_nuq(coeff, bs * bs, x->skip_block,
                                   p->quant_pxd, p->quant_shift_pxd,
                                   pd->dequant_pxd,
                                   (const cumbins_type_nuq *)
                                   p->cumbins_nuq_pxd,
                                   (const dequant_val_type_nuq *)
                                   pd->dequant_val_nuq_pxd,
                                   qcoeff, dqcoeff, eob,
                                   scan_order->scan, band);
        }
#else
        if (x->quant_fp)
          vp9_quantize_64x64_fp_nuq(coeff, bs * bs, x->skip_block,
                                    p->quant_pxd_fp, pd->dequant_pxd,
                                    (const cumbins_type_nuq *)
                                    p->cumbins_nuq_pxd,
                                    (const dequant_val_type_nuq *)
                                    pd->dequant_val_nuq_pxd,
                                    qcoeff, dqcoeff, eob,
                                    scan_order->scan, band);
        else
          vp9_quantize_64x64_nuq(coeff, bs * bs, x->skip_block,
                                 p->quant_pxd, p->quant_shift_pxd,
                                 pd->dequant_pxd,
                                 (const cumbins_type_nuq *)p->cumbins_nuq_pxd,
                                 (const dequant_val_type_nuq *)
                                 pd->dequant_val_nuq_pxd,
                                 qcoeff, dqcoeff, eob,
                                 scan_order->scan, band);
#endif  // CONFIG_VP9_HIGHBITDEPTH
#else  // CONFIG_NEW_QUANT
#if CONFIG_VP9_HIGHBITDEPTH
        if (use_hbd)
          vp9_highbd_quantize_b_64x64(coeff, bs * bs, x->skip_block,
                                      p->zbin_pxd, p->round_pxd, p->quant_pxd,
                                      p->quant_shift_pxd, qcoeff, dqcoeff,
                                      pd->dequant_pxd, eob,
                                      scan_order->scan,  scan_order->iscan);
        else
#endif  // CONFIG_VP9_HIGHBITDEPTH
        vp9_quantize_b_64x64(coeff, bs * bs, x->skip_block, p->zbin_pxd,
                             p->round_pxd, p->quant_pxd, p->quant_shift_pxd,
                             qcoeff, dqcoeff, pd->dequant_pxd, eob,
                             scan_order->scan,  scan_order->iscan);

#endif  // CONFIG_NEW_QUANT
      }
#endif  // CONFIG_TX64X64
    }

    if (!x->skip_encode && *eob) {
#if CONFIG_VP9_HIGHBITDEPTH
      if (use_hbd)
        vp9_highbd_tx_identity_add(dqcoeff, dst, dst_stride, 4 << tx_size,
                                   shift, xd->bd);
      else
#endif  // CONFIG_VP9_HIGHBITDEPTH
        vp9_tx_identity_add(dqcoeff, dst, dst_stride, 4 << tx_size, shift);
    }

    if (*eob)
      *(args->skip) = 0;
    return;
  }
#endif  // CONFIG_TX_SKIP

#if CONFIG_VP9_HIGHBITDEPTH
  if (xd->cur_buf->flags & YV12_FLAG_HIGHBITDEPTH) {
    switch (tx_size) {
#if CONFIG_TX64X64
      case TX_64X64:
        scan_order = &vp9_default_scan_orders[TX_64X64];
        mode = plane == 0 ? mbmi->mode : mbmi->uv_mode;
        vp9_predict_intra_block(xd, block >> 8, bwl, TX_64X64, mode,
#if CONFIG_FILTERINTRA
                                fbit,
#endif
                                x->skip_encode ? src : dst,
                                x->skip_encode ? src_stride : dst_stride,
                                dst, dst_stride, i, j, plane);
        if (!x->skip_recode) {
          vp9_highbd_subtract_block(64, 64, src_diff, diff_stride,
                                    src, src_stride, dst, dst_stride, xd->bd);
          vp9_highbd_fdct64x64(src_diff, coeff, diff_stride);
#if CONFIG_NEW_QUANT
          if (x->quant_fp)
            vp9_highbd_quantize_64x64_fp_nuq(coeff, 4096, x->skip_block,
                                             p->quant_fp, pd->dequant,
                                             (const cumbins_type_nuq *)
                                             p->cumbins_nuq,
                                             (const dequant_val_type_nuq *)
                                             pd->dequant_val_nuq,
                                             qcoeff, dqcoeff, eob,
                                             scan_order->scan,
                                             band);
          else
            vp9_highbd_quantize_64x64_nuq(coeff, 4096, x->skip_block,
                                          p->quant, p->quant_shift, pd->dequant,
                                          (const cumbins_type_nuq *)
                                          p->cumbins_nuq,
                                          (const dequant_val_type_nuq *)
                                          pd->dequant_val_nuq,
                                          qcoeff, dqcoeff, eob,
                                          scan_order->scan, band);
#else
          vp9_highbd_quantize_b_64x64(coeff, 4096, x->skip_block, p->zbin,
                                      p->round, p->quant, p->quant_shift,
                                      qcoeff, dqcoeff, pd->dequant, eob,
                                      scan_order->scan, scan_order->iscan);
#endif  // CONFIG_NEW_QUANT
          if (!x->skip_encode && *eob) {
            vp9_highbd_idct64x64_add(dqcoeff, dst, dst_stride, *eob, xd->bd);
          }
        }
        break;
#endif  // CONFIG_TX64X64
      case TX_32X32:
        scan_order = &vp9_default_scan_orders[TX_32X32];
        mode = plane == 0 ? mbmi->mode : mbmi->uv_mode;
        vp9_predict_intra_block(xd, block >> 6, bwl, TX_32X32, mode,
#if CONFIG_FILTERINTRA
                                fbit,
#endif
                                x->skip_encode ? src : dst,
                                x->skip_encode ? src_stride : dst_stride,
                                dst, dst_stride, i, j, plane);
        if (!x->skip_recode) {
          vp9_highbd_subtract_block(32, 32, src_diff, diff_stride,
                                    src, src_stride, dst, dst_stride, xd->bd);
          highbd_fdct32x32(x->use_lp32x32fdct, src_diff, coeff, diff_stride);
#if CONFIG_NEW_QUANT
          if (x->quant_fp)
            vp9_highbd_quantize_32x32_fp_nuq(coeff, 1024, x->skip_block,
                                             p->quant_fp, pd->dequant,
                                             (const cumbins_type_nuq *)
                                             p->cumbins_nuq,
                                             (const dequant_val_type_nuq *)
                                             pd->dequant_val_nuq,
                                             qcoeff, dqcoeff, eob,
                                             scan_order->scan,
                                             band);
          else
            vp9_highbd_quantize_32x32_nuq(coeff, 1024, x->skip_block,
                                          p->quant, p->quant_shift, pd->dequant,
                                          (const cumbins_type_nuq *)
                                          p->cumbins_nuq,
                                          (const dequant_val_type_nuq *)
                                          pd->dequant_val_nuq,
                                          qcoeff, dqcoeff, eob,
                                          scan_order->scan, band);
#else
          vp9_highbd_quantize_b_32x32(coeff, 1024, x->skip_block, p->zbin,
                                      p->round, p->quant, p->quant_shift,
                                      qcoeff, dqcoeff, pd->dequant, eob,
                                      scan_order->scan, scan_order->iscan);
#endif  // CONFIG_NEW_QUANT
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
#if CONFIG_FILTERINTRA
                                fbit,
#endif
                                x->skip_encode ? src : dst,
                                x->skip_encode ? src_stride : dst_stride,
                                dst, dst_stride, i, j, plane);
        if (!x->skip_recode) {
          vp9_highbd_subtract_block(16, 16, src_diff, diff_stride,
                                    src, src_stride, dst, dst_stride, xd->bd);
          vp9_highbd_fht16x16(src_diff, coeff, diff_stride, tx_type);
#if CONFIG_NEW_QUANT
          if (x->quant_fp)
            vp9_highbd_quantize_fp_nuq(coeff, 256, x->skip_block,
                                       p->quant_fp, pd->dequant,
                                       (const cumbins_type_nuq *)p->cumbins_nuq,
                                       (const dequant_val_type_nuq *)
                                       pd->dequant_val_nuq,
                                       qcoeff, dqcoeff, eob,
                                       scan_order->scan, band);
          else
            vp9_highbd_quantize_nuq(coeff, 256, x->skip_block,
                                    p->quant, p->quant_shift, pd->dequant,
                                    (const cumbins_type_nuq *)p->cumbins_nuq,
                                    (const dequant_val_type_nuq *)
                                    pd->dequant_val_nuq,
                                    qcoeff, dqcoeff, eob,
                                    scan_order->scan, band);
#else
          vp9_highbd_quantize_b(coeff, 256, x->skip_block, p->zbin, p->round,
                                p->quant, p->quant_shift, qcoeff, dqcoeff,
                                pd->dequant, eob,
                                scan_order->scan, scan_order->iscan);
#endif  // CONFIG_NEW_QUANT
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
#if CONFIG_FILTERINTRA
                                fbit,
#endif
                                x->skip_encode ? src : dst,
                                x->skip_encode ? src_stride : dst_stride,
                                dst, dst_stride, i, j, plane);
        if (!x->skip_recode) {
          vp9_highbd_subtract_block(8, 8, src_diff, diff_stride,
                                    src, src_stride, dst, dst_stride, xd->bd);
          vp9_highbd_fht8x8(src_diff, coeff, diff_stride, tx_type);
#if CONFIG_NEW_QUANT
          if (x->quant_fp)
            vp9_highbd_quantize_fp_nuq(coeff, 64, x->skip_block,
                                       p->quant_fp, pd->dequant,
                                       (const cumbins_type_nuq *)p->cumbins_nuq,
                                       (const dequant_val_type_nuq *)
                                       pd->dequant_val_nuq,
                                       qcoeff, dqcoeff, eob,
                                       scan_order->scan, band);
          else
            vp9_highbd_quantize_nuq(coeff, 64, x->skip_block,
                                    p->quant, p->quant_shift, pd->dequant,
                                    (const cumbins_type_nuq *)p->cumbins_nuq,
                                    (const dequant_val_type_nuq *)
                                    pd->dequant_val_nuq,
                                    qcoeff, dqcoeff, eob,
                                    scan_order->scan, band);
#else
          vp9_highbd_quantize_b(coeff, 64, x->skip_block, p->zbin, p->round,
                                p->quant, p->quant_shift, qcoeff, dqcoeff,
                                pd->dequant, eob,
                                scan_order->scan, scan_order->iscan);
#endif  // CONFIG_NEW_QUANT
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
#if CONFIG_FILTERINTRA
                                fbit,
#endif
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
#if CONFIG_NEW_QUANT
          if (x->quant_fp)
            vp9_highbd_quantize_fp_nuq(coeff, 16, x->skip_block,
                                       p->quant_fp, pd->dequant,
                                       (const cumbins_type_nuq *)p->cumbins_nuq,
                                       (const dequant_val_type_nuq *)
                                       pd->dequant_val_nuq,
                                       qcoeff, dqcoeff, eob,
                                       scan_order->scan, band);
          else
            vp9_highbd_quantize_nuq(coeff, 16, x->skip_block,
                                    p->quant, p->quant_shift, pd->dequant,
                                    (const cumbins_type_nuq *)p->cumbins_nuq,
                                    (const dequant_val_type_nuq *)
                                    pd->dequant_val_nuq,
                                    qcoeff, dqcoeff, eob,
                                    scan_order->scan, band);
#else
          vp9_highbd_quantize_b(coeff, 16, x->skip_block, p->zbin, p->round,
                                p->quant, p->quant_shift, qcoeff, dqcoeff,
                                pd->dequant, eob,
                                scan_order->scan, scan_order->iscan);
#endif  // CONFIG_NEW_QUANT
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
#if CONFIG_TX64X64
    case TX_64X64:
      assert(plane == 0);
      scan_order = &vp9_default_scan_orders[TX_64X64];
      mode = plane == 0 ? mbmi->mode : mbmi->uv_mode;
      vp9_predict_intra_block(xd, block >> 8, bwl, TX_64X64, mode,
#if CONFIG_FILTERINTRA
                              fbit,
#endif
                              x->skip_encode ? src : dst,
                              x->skip_encode ? src_stride : dst_stride,
                              dst, dst_stride, i, j, plane);
      if (!x->skip_recode) {
        vp9_subtract_block(64, 64, src_diff, diff_stride,
                           src, src_stride, dst, dst_stride);
        vp9_fdct64x64(src_diff, coeff, diff_stride);
#if CONFIG_NEW_QUANT
        if (x->quant_fp)
          vp9_quantize_64x64_fp_nuq(coeff, 4096, x->skip_block,
                                    p->quant_fp, pd->dequant,
                                    (const cumbins_type_nuq *)p->cumbins_nuq,
                                    (const dequant_val_type_nuq *)
                                    pd->dequant_val_nuq,
                                    qcoeff, dqcoeff, eob,
                                    scan_order->scan, band);
        else
          vp9_quantize_64x64_nuq(coeff, 4096, x->skip_block,
                                 p->quant, p->quant_shift, pd->dequant,
                                 (const cumbins_type_nuq *)p->cumbins_nuq,
                                 (const dequant_val_type_nuq *)
                                 pd->dequant_val_nuq,
                                 qcoeff, dqcoeff, eob,
                                 scan_order->scan, band);
#else
        vp9_quantize_b_64x64(coeff, 4096, x->skip_block, p->zbin, p->round,
                             p->quant, p->quant_shift, qcoeff, dqcoeff,
                             pd->dequant, eob, scan_order->scan,
                             scan_order->iscan);
#endif  // CONFIG_NEW_QUANT
      }
      if (!x->skip_encode && *eob)
        vp9_idct64x64_add(dqcoeff, dst, dst_stride, *eob);
      break;
#endif  // CONFIG_TX64X64
    case TX_32X32:
      scan_order = &vp9_default_scan_orders[TX_32X32];
      mode = plane == 0 ? mbmi->mode : mbmi->uv_mode;
      vp9_predict_intra_block(xd, block >> 6, bwl, TX_32X32, mode,
#if CONFIG_FILTERINTRA
                              fbit,
#endif
                              x->skip_encode ? src : dst,
                              x->skip_encode ? src_stride : dst_stride,
                              dst, dst_stride, i, j, plane);
      if (!x->skip_recode) {
        vp9_subtract_block(32, 32, src_diff, diff_stride,
                           src, src_stride, dst, dst_stride);
        fdct32x32(x->use_lp32x32fdct, src_diff, coeff, diff_stride);
#if CONFIG_NEW_QUANT
        if (x->quant_fp)
          vp9_quantize_32x32_fp_nuq(coeff, 1024, x->skip_block,
                                    p->quant_fp, pd->dequant,
                                    (const cumbins_type_nuq *)p->cumbins_nuq,
                                    (const dequant_val_type_nuq *)
                                    pd->dequant_val_nuq,
                                    qcoeff, dqcoeff, eob,
                                    scan_order->scan, band);
        else
          vp9_quantize_32x32_nuq(coeff, 1024, x->skip_block,
                                 p->quant, p->quant_shift, pd->dequant,
                                 (const cumbins_type_nuq *)p->cumbins_nuq,
                                 (const dequant_val_type_nuq *)
                                 pd->dequant_val_nuq,
                                 qcoeff, dqcoeff, eob,
                                 scan_order->scan, band);
#else
        vp9_quantize_b_32x32(coeff, 1024, x->skip_block, p->zbin, p->round,
                             p->quant, p->quant_shift, qcoeff, dqcoeff,
                             pd->dequant, eob, scan_order->scan,
                             scan_order->iscan);
#endif  // CONFIG_NEW_QUANT
      }
      if (!x->skip_encode && *eob)
        vp9_idct32x32_add(dqcoeff, dst, dst_stride, *eob);
      break;
    case TX_16X16:
      tx_type = get_tx_type(pd->plane_type, xd);
      scan_order = &vp9_scan_orders[TX_16X16][tx_type];
      mode = plane == 0 ? mbmi->mode : mbmi->uv_mode;
      vp9_predict_intra_block(xd, block >> 4, bwl, TX_16X16, mode,
#if CONFIG_FILTERINTRA
                              fbit,
#endif
                              x->skip_encode ? src : dst,
                              x->skip_encode ? src_stride : dst_stride,
                              dst, dst_stride, i, j, plane);
      if (!x->skip_recode) {
        vp9_subtract_block(16, 16, src_diff, diff_stride,
                           src, src_stride, dst, dst_stride);
        vp9_fht16x16(src_diff, coeff, diff_stride, tx_type);
#if CONFIG_NEW_QUANT
        if (x->quant_fp)
          vp9_quantize_fp_nuq(coeff, 256, x->skip_block,
                              p->quant_fp, pd->dequant,
                              (const cumbins_type_nuq *)p->cumbins_nuq,
                              (const dequant_val_type_nuq *)pd->dequant_val_nuq,
                              qcoeff, dqcoeff, eob,
                              scan_order->scan, band);
        else
          vp9_quantize_nuq(coeff, 256, x->skip_block,
                           p->quant, p->quant_shift, pd->dequant,
                           (const cumbins_type_nuq *)p->cumbins_nuq,
                           (const dequant_val_type_nuq *)pd->dequant_val_nuq,
                           qcoeff, dqcoeff, eob,
                           scan_order->scan, band);
#else
        vp9_quantize_b(coeff, 256, x->skip_block, p->zbin, p->round,
                       p->quant, p->quant_shift, qcoeff, dqcoeff,
                       pd->dequant, eob, scan_order->scan,
                       scan_order->iscan);
#endif  // CONFIG_NEW_QUANT
      }
      if (!x->skip_encode && *eob)
        vp9_iht16x16_add(tx_type, dqcoeff, dst, dst_stride, *eob);
      break;
    case TX_8X8:
      tx_type = get_tx_type(pd->plane_type, xd);
      scan_order = &vp9_scan_orders[TX_8X8][tx_type];
      mode = plane == 0 ? mbmi->mode : mbmi->uv_mode;
      vp9_predict_intra_block(xd, block >> 2, bwl, TX_8X8, mode,
#if CONFIG_FILTERINTRA
                              fbit,
#endif
                              x->skip_encode ? src : dst,
                              x->skip_encode ? src_stride : dst_stride,
                              dst, dst_stride, i, j, plane);
      if (!x->skip_recode) {
        vp9_subtract_block(8, 8, src_diff, diff_stride,
                           src, src_stride, dst, dst_stride);
        vp9_fht8x8(src_diff, coeff, diff_stride, tx_type);
#if CONFIG_NEW_QUANT
        if (x->quant_fp)
          vp9_quantize_fp_nuq(coeff, 64, x->skip_block,
                              p->quant_fp, pd->dequant,
                              (const cumbins_type_nuq *)p->cumbins_nuq,
                              (const dequant_val_type_nuq *)pd->dequant_val_nuq,
                              qcoeff, dqcoeff, eob,
                              scan_order->scan, band);
        else
          vp9_quantize_nuq(coeff, 64, x->skip_block,
                           p->quant, p->quant_shift, pd->dequant,
                           (const cumbins_type_nuq *)p->cumbins_nuq,
                           (const dequant_val_type_nuq *)pd->dequant_val_nuq,
                           qcoeff, dqcoeff, eob,
                           scan_order->scan, band);
#else
        vp9_quantize_b(coeff, 64, x->skip_block, p->zbin, p->round, p->quant,
                       p->quant_shift, qcoeff, dqcoeff,
                       pd->dequant, eob, scan_order->scan,
                       scan_order->iscan);
#endif  // CONFIG_NEW_QUANT
      }
      if (!x->skip_encode && *eob)
        vp9_iht8x8_add(tx_type, dqcoeff, dst, dst_stride, *eob);
      break;
    case TX_4X4:
      tx_type = get_tx_type_4x4(pd->plane_type, xd, block);
      scan_order = &vp9_scan_orders[TX_4X4][tx_type];
      mode = plane == 0 ? get_y_mode(xd->mi[0].src_mi, block) : mbmi->uv_mode;
      vp9_predict_intra_block(xd, block, bwl, TX_4X4, mode,
#if CONFIG_FILTERINTRA
                              fbit,
#endif
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
#if CONFIG_NEW_QUANT
        if (x->quant_fp)
          vp9_quantize_fp_nuq(coeff, 16, x->skip_block,
                              p->quant_fp, pd->dequant,
                              (const cumbins_type_nuq *)p->cumbins_nuq,
                              (const dequant_val_type_nuq *)pd->dequant_val_nuq,
                              qcoeff, dqcoeff, eob,
                              scan_order->scan, band);
        else
          vp9_quantize_nuq(coeff, 16, x->skip_block,
                           p->quant, p->quant_shift, pd->dequant,
                           (const cumbins_type_nuq *)p->cumbins_nuq,
                           (const dequant_val_type_nuq *)pd->dequant_val_nuq,
                           qcoeff, dqcoeff, eob,
                           scan_order->scan, band);
#else
        vp9_quantize_b(coeff, 16, x->skip_block, p->zbin, p->round, p->quant,
                       p->quant_shift, qcoeff, dqcoeff,
                       pd->dequant, eob, scan_order->scan,
                       scan_order->iscan);
#endif  // CONFIG_NEW_QUANT
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

void vp9_encode_block_intra(MACROBLOCK *x, int plane, int block,
                            BLOCK_SIZE plane_bsize, TX_SIZE tx_size,
                            int8_t *skip) {
  struct encode_b_args arg = {x, NULL, skip};
  encode_block_intra(plane, block, plane_bsize, tx_size, &arg);
}


void vp9_encode_intra_block_plane(MACROBLOCK *x, BLOCK_SIZE bsize, int plane) {
  const MACROBLOCKD *const xd = &x->e_mbd;
  struct encode_b_args arg = {x, NULL, &xd->mi[0].src_mi->mbmi.skip};

  vp9_foreach_transformed_block_in_plane(xd, bsize, plane, encode_block_intra,
                                         &arg);
}
