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
#include "vp9/common/vp9_invtrans.h"
#include "vp9/common/vp9_reconintra.h"
#include "vpx_mem/vpx_mem.h"
#include "vp9/encoder/vp9_rdopt.h"
#include "vp9/common/vp9_systemdependent.h"
#include "vp9_rtcd.h"

void vp9_subtract_b_c(BLOCK *be, BLOCKD *bd, int pitch) {
  uint8_t *src_ptr = (*(be->base_src) + be->src);
  int16_t *diff_ptr = be->src_diff;
  uint8_t *pred_ptr = bd->predictor;
  int src_stride = be->src_stride;

  int r, c;

  for (r = 0; r < 4; r++) {
    for (c = 0; c < 4; c++)
      diff_ptr[c] = src_ptr[c] - pred_ptr[c];

    diff_ptr += pitch;
    pred_ptr += pitch;
    src_ptr  += src_stride;
  }
}

void vp9_subtract_4b_c(BLOCK *be, BLOCKD *bd, int pitch) {
  uint8_t *src_ptr = (*(be->base_src) + be->src);
  int16_t *diff_ptr = be->src_diff;
  uint8_t *pred_ptr = bd->predictor;
  int src_stride = be->src_stride;
  int r, c;

  for (r = 0; r < 8; r++) {
    for (c = 0; c < 8; c++)
      diff_ptr[c] = src_ptr[c] - pred_ptr[c];

    diff_ptr += pitch;
    pred_ptr += pitch;
    src_ptr  += src_stride;
  }
}

void vp9_subtract_sby_s_c(int16_t *diff, const uint8_t *src, int src_stride,
                          const uint8_t *pred, int dst_stride,
                          BLOCK_SIZE_TYPE bsize) {
  const int bh = 16 << mb_height_log2(bsize), bw = 16 << mb_width_log2(bsize);
  int r, c;

  for (r = 0; r < bh; r++) {
    for (c = 0; c < bw; c++)
      diff[c] = src[c] - pred[c];

    diff += bw;
    pred += dst_stride;
    src  += src_stride;
  }
}

void vp9_subtract_sbuv_s_c(int16_t *diff, const uint8_t *usrc,
                           const uint8_t *vsrc, int src_stride,
                           const uint8_t *upred,
                           const uint8_t *vpred, int dst_stride,
                           BLOCK_SIZE_TYPE bsize) {
  const int bhl = mb_height_log2(bsize), bwl = mb_width_log2(bsize);
  const int uoff = (16 * 16) << (bhl + bwl), voff = (uoff * 5) >> 2;
  const int bw = 8 << bwl, bh = 8 << bhl;
  int16_t *udiff = diff + uoff;
  int16_t *vdiff = diff + voff;
  int r, c;

  for (r = 0; r < bh; r++) {
    for (c = 0; c < bw; c++)
      udiff[c] = usrc[c] - upred[c];

    udiff += bw;
    upred += dst_stride;
    usrc  += src_stride;
  }

  for (r = 0; r < bh; r++) {
    for (c = 0; c < bw; c++)
      vdiff[c] = vsrc[c] - vpred[c];

    vdiff += bw;
    vpred += dst_stride;
    vsrc  += src_stride;
  }
}

void vp9_subtract_mby_c(int16_t *diff, uint8_t *src,
                        uint8_t *pred, int stride) {
  vp9_subtract_sby_s_c(diff, src, stride, pred, 16, BLOCK_SIZE_MB16X16);
}

void vp9_subtract_mbuv_c(int16_t *diff, uint8_t *usrc,
                         uint8_t *vsrc, uint8_t *pred, int stride) {
  uint8_t *upred = pred + 256;
  uint8_t *vpred = pred + 320;

  vp9_subtract_sbuv_s_c(diff, usrc, vsrc, stride, upred, vpred, 8,
                        BLOCK_SIZE_MB16X16);
}

static void subtract_mb(MACROBLOCK *x) {
  vp9_subtract_mby(x->src_diff, x->src.y_buffer, x->e_mbd.predictor,
                   x->src.y_stride);
  vp9_subtract_mbuv(x->src_diff, x->src.u_buffer, x->src.v_buffer,
                    x->e_mbd.predictor, x->src.uv_stride);
}

void vp9_transform_sby_32x32(MACROBLOCK *x, BLOCK_SIZE_TYPE bsize) {
  const int bwl = mb_width_log2(bsize) - 1, bw = 1 << bwl;
  const int bh = 1 << (mb_height_log2(bsize) - 1);
  const int stride = 32 << bwl;
  int n;

  for (n = 0; n < bw * bh; n++) {
    const int x_idx = n & (bw - 1), y_idx = n >> bwl;

    vp9_short_fdct32x32(x->src_diff + y_idx * stride * 32 + x_idx * 32,
                        x->coeff + n * 1024, stride * 2);
  }
}

void vp9_transform_sby_16x16(MACROBLOCK *x, BLOCK_SIZE_TYPE bsize) {
  const int bwl = mb_width_log2(bsize), bw = 1 << bwl;
  const int bh = 1 << mb_height_log2(bsize);
  const int stride = 16 << bwl, bstride = 4 << bwl;
  MACROBLOCKD *const xd = &x->e_mbd;
  int n;

  for (n = 0; n < bw * bh; n++) {
    const int x_idx = n & (bw - 1), y_idx = n >> bwl;
    const TX_TYPE tx_type = get_tx_type_16x16(xd,
                                              (y_idx * bstride + x_idx) * 4);

    if (tx_type != DCT_DCT) {
      vp9_short_fht16x16(x->src_diff + y_idx * stride * 16 + x_idx * 16,
                         x->coeff + n * 256, stride, tx_type);
    } else {
      x->fwd_txm16x16(x->src_diff + y_idx * stride * 16 + x_idx * 16,
                      x->coeff + n * 256, stride * 2);
    }
  }
}

void vp9_transform_sby_8x8(MACROBLOCK *x, BLOCK_SIZE_TYPE bsize) {
  const int bwl = mb_width_log2(bsize) + 1, bw = 1 << bwl;
  const int bh = 1 << (mb_height_log2(bsize) + 1);
  const int stride = 8 << bwl, bstride = 2 << bwl;
  MACROBLOCKD *const xd = &x->e_mbd;
  int n;

  for (n = 0; n < bw * bh; n++) {
    const int x_idx = n & (bw - 1), y_idx = n >> bwl;
    const TX_TYPE tx_type = get_tx_type_8x8(xd, (y_idx * bstride + x_idx) * 2);

    if (tx_type != DCT_DCT) {
      vp9_short_fht8x8(x->src_diff + y_idx * stride * 8 + x_idx * 8,
                       x->coeff + n * 64, stride, tx_type);
    } else {
      x->fwd_txm8x8(x->src_diff + y_idx * stride * 8 + x_idx * 8,
                    x->coeff + n * 64, stride * 2);
    }
  }
}

void vp9_transform_sby_4x4(MACROBLOCK *x, BLOCK_SIZE_TYPE bsize) {
  const int bwl = mb_width_log2(bsize) + 2, bw = 1 << bwl;
  const int bh = 1 << (mb_height_log2(bsize) + 2);
  const int stride = 4 << bwl;
  MACROBLOCKD *const xd = &x->e_mbd;
  int n;

  for (n = 0; n < bw * bh; n++) {
    const int x_idx = n & (bw - 1), y_idx = n >> bwl;
    const TX_TYPE tx_type = get_tx_type_4x4(xd, n);

    if (tx_type != DCT_DCT) {
      vp9_short_fht4x4(x->src_diff + y_idx * stride * 4 + x_idx * 4,
                       x->coeff + n * 16, stride, tx_type);
    } else {
      x->fwd_txm4x4(x->src_diff + y_idx * stride * 4 + x_idx * 4,
                    x->coeff + n * 16, stride * 2);
    }
  }
}

void vp9_transform_sbuv_32x32(MACROBLOCK *x, BLOCK_SIZE_TYPE bsize) {
  assert(bsize == BLOCK_SIZE_SB64X64);
  vp9_clear_system_state();
  vp9_short_fdct32x32(x->src_diff + 4096,
                      x->coeff + 4096, 64);
  vp9_short_fdct32x32(x->src_diff + 4096 + 1024,
                      x->coeff + 4096 + 1024, 64);
}

void vp9_transform_sbuv_16x16(MACROBLOCK *x, BLOCK_SIZE_TYPE bsize) {
  const int bwl = mb_width_log2(bsize), bhl = mb_height_log2(bsize);
  const int uoff = (16 * 16) << (bwl + bhl), voff = (uoff * 5) >> 2;
  const int bw = 1 << (bwl - 1), bh = 1 << (bhl - 1);
  const int stride = 16 << (bwl - 1);
  int n;

  vp9_clear_system_state();
  for (n = 0; n < bw * bh; n++) {
    const int x_idx = n & (bw - 1), y_idx = n >> (bwl - 1);

    x->fwd_txm16x16(x->src_diff + uoff + y_idx * stride * 16 + x_idx * 16,
                    x->coeff + uoff + n * 256, stride * 2);
    x->fwd_txm16x16(x->src_diff + voff + y_idx * stride * 16 + x_idx * 16,
                    x->coeff + voff + n * 256, stride * 2);
  }
}

void vp9_transform_sbuv_8x8(MACROBLOCK *x, BLOCK_SIZE_TYPE bsize) {
  const int bwl = mb_width_log2(bsize) + 1, bhl = mb_height_log2(bsize) + 1;
  const int uoff = (8 * 8) << (bwl + bhl), voff = (uoff * 5) >> 2;
  const int bw = 1 << (bwl - 1), bh = 1 << (bhl - 1);
  const int stride = 8 << (bwl - 1);
  int n;

  vp9_clear_system_state();
  for (n = 0; n < bw * bh; n++) {
    const int x_idx = n & (bw - 1), y_idx = n >> (bwl - 1);

    x->fwd_txm8x8(x->src_diff + uoff + y_idx * stride * 8 + x_idx * 8,
                  x->coeff + uoff + n * 64, stride * 2);
    x->fwd_txm8x8(x->src_diff + voff + y_idx * stride * 8 + x_idx * 8,
                  x->coeff + voff + n * 64, stride * 2);
  }
}

void vp9_transform_sbuv_4x4(MACROBLOCK *x, BLOCK_SIZE_TYPE bsize) {
  const int bwl = mb_width_log2(bsize) + 2, bhl = mb_height_log2(bsize) + 2;
  const int uoff = (4 * 4) << (bwl + bhl), voff = (uoff * 5) >> 2;
  const int bw = 1 << (bwl - 1), bh = 1 << (bhl - 1);
  const int stride = 4 << (bwl - 1);
  int n;

  vp9_clear_system_state();
  for (n = 0; n < bw * bh; n++) {
    const int x_idx = n & (bw - 1), y_idx = n >> (bwl - 1);

    x->fwd_txm4x4(x->src_diff + uoff + y_idx * stride * 4 + x_idx * 4,
                  x->coeff + uoff + n * 16, stride * 2);
    x->fwd_txm4x4(x->src_diff + voff + y_idx * stride * 4 + x_idx * 4,
                  x->coeff + voff + n * 16, stride * 2);
  }
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
  int bak = token_cache[idx], pt;
  token_cache[idx] = token;
  pt = vp9_get_coef_context(scan, nb, pad, token_cache, idx + 1, l);
  token_cache[idx] = bak;
  return pt;
}

static void optimize_b(VP9_COMMON *const cm,
                       MACROBLOCK *mb, int ib, PLANE_TYPE type,
                       const int16_t *dequant_ptr,
                       ENTROPY_CONTEXT *a, ENTROPY_CONTEXT *l,
                       int tx_size, int y_blocks) {
  const int ref = mb->e_mbd.mode_info_context->mbmi.ref_frame != INTRA_FRAME;
  MACROBLOCKD *const xd = &mb->e_mbd;
  vp9_token_state tokens[1025][2];
  unsigned best_index[1025][2];
  const struct plane_block_idx pb_idx = plane_block_idx(y_blocks, ib);
  const int16_t *coeff_ptr = mb->coeff + ib * 16;
  int16_t *qcoeff_ptr;
  int16_t *dqcoeff_ptr;
  int eob = xd->plane[pb_idx.plane].eobs[pb_idx.block], final_eob, sz = 0;
  const int i0 = 0;
  int rc, x, next, i;
  int64_t rdmult, rddiv, rd_cost0, rd_cost1;
  int rate0, rate1, error0, error1, t0, t1;
  int best, band, pt;
  int err_mult = plane_rd_mult[type];
  int default_eob, pad;
  int const *scan, *nb;
  const int mul = 1 + (tx_size == TX_32X32);
  uint8_t token_cache[1024];
#if CONFIG_CODE_NONZEROCOUNT
  // TODO(debargha): the dynamic programming approach used in this function
  // is not compatible with the true rate cost when nzcs are used. Note
  // the total rate is the sum of the nzc rate and the indicvidual token
  // rates. The latter part can be optimized in this function, but because
  // the nzc rate is a function of all the other tokens without a Markov
  // relationship this rate cannot be considered correctly.
  // The current implementation uses a suboptimal approach to account for
  // the nzc rates somewhat, but in reality the optimization approach needs
  // to change substantially.
  const int nzc_used = get_nzc_used(tx_size);
  uint16_t nzc = xd->nzcs[ib];
  uint16_t nzc0, nzc1;
  uint16_t final_nzc = 0, final_nzc_exp;
  int nzc_context = vp9_get_nzc_context(cm, xd, ib);
  unsigned int *nzc_cost;
  nzc0 = nzc1 = nzc;
#endif

  assert((!type && !pb_idx.plane) || (type && pb_idx.plane));
  dqcoeff_ptr = BLOCK_OFFSET(xd->plane[pb_idx.plane].dqcoeff, pb_idx.block, 16);
  qcoeff_ptr = BLOCK_OFFSET(xd->plane[pb_idx.plane].qcoeff, pb_idx.block, 16);
  switch (tx_size) {
    default:
    case TX_4X4: {
      const TX_TYPE tx_type = get_tx_type_4x4(xd, ib);
      default_eob = 16;
#if CONFIG_CODE_NONZEROCOUNT
      nzc_cost = mb->nzc_costs_4x4[nzc_context][ref][type];
#endif
      if (tx_type == DCT_ADST) {
        scan = vp9_col_scan_4x4;
      } else if (tx_type == ADST_DCT) {
        scan = vp9_row_scan_4x4;
      } else {
        scan = vp9_default_zig_zag1d_4x4;
      }
      break;
    }
    case TX_8X8: {
      const BLOCK_SIZE_TYPE sb_type = xd->mode_info_context->mbmi.sb_type;
      const int sz = 3 + mb_width_log2(sb_type);
      const int x = ib & ((1 << sz) - 1), y = ib - x;
      const TX_TYPE tx_type = get_tx_type_8x8(xd, y + (x >> 1));
      if (tx_type == DCT_ADST) {
        scan = vp9_col_scan_8x8;
      } else if (tx_type == ADST_DCT) {
        scan = vp9_row_scan_8x8;
      } else {
        scan = vp9_default_zig_zag1d_8x8;
      }
      default_eob = 64;
#if CONFIG_CODE_NONZEROCOUNT
      nzc_cost = mb->nzc_costs_8x8[nzc_context][ref][type];
#endif
      break;
    }
    case TX_16X16: {
      const BLOCK_SIZE_TYPE sb_type = xd->mode_info_context->mbmi.sb_type;
      const int sz = 4 + mb_width_log2(sb_type);
      const int x = ib & ((1 << sz) - 1), y = ib - x;
      const TX_TYPE tx_type = get_tx_type_16x16(xd, y + (x >> 2));
      if (tx_type == DCT_ADST) {
        scan = vp9_col_scan_16x16;
      } else if (tx_type == ADST_DCT) {
        scan = vp9_row_scan_16x16;
      } else {
        scan = vp9_default_zig_zag1d_16x16;
      }
      default_eob = 256;
#if CONFIG_CODE_NONZEROCOUNT
      nzc_cost = mb->nzc_costs_16x16[nzc_context][ref][type];
#endif
      break;
    }
    case TX_32X32:
      scan = vp9_default_zig_zag1d_32x32;
      default_eob = 1024;
#if CONFIG_CODE_NONZEROCOUNT
      nzc_cost = mb->nzc_costs_32x32[nzc_context][ref][type];
#endif
      break;
  }
  assert(eob <= default_eob);

  /* Now set up a Viterbi trellis to evaluate alternative roundings. */
  rdmult = mb->rdmult * err_mult;
  if (mb->e_mbd.mode_info_context->mbmi.ref_frame == INTRA_FRAME)
    rdmult = (rdmult * 9) >> 4;
  rddiv = mb->rddiv;
  memset(best_index, 0, sizeof(best_index));
  /* Initialize the sentinel node of the trellis. */
#if CONFIG_CODE_NONZEROCOUNT
  tokens[eob][0].rate = nzc_used ? nzc_cost[nzc] : 0;
#else
  tokens[eob][0].rate = 0;
#endif
  tokens[eob][0].error = 0;
  tokens[eob][0].next = default_eob;
  tokens[eob][0].token = DCT_EOB_TOKEN;
  tokens[eob][0].qc = 0;
  *(tokens[eob] + 1) = *(tokens[eob] + 0);
  next = eob;
  for (i = 0; i < eob; i++)
    token_cache[i] = vp9_dct_value_tokens_ptr[qcoeff_ptr[scan[i]]].Token;
  nb = vp9_get_coef_neighbors_handle(scan, &pad);

  for (i = eob; i-- > i0;) {
    int base_bits, d2, dx;
#if CONFIG_CODE_NONZEROCOUNT
    int new_nzc0, new_nzc1;
#endif

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
      t0 = (vp9_dct_value_tokens_ptr + x)->Token;
      /* Consider both possible successor states. */
      if (next < default_eob) {
        band = get_coef_band(scan, tx_size, i + 1);
        pt = trellis_get_coeff_context(scan, nb, i, t0, token_cache,
                                       pad, default_eob);
        rate0 +=
          mb->token_costs[tx_size][type][ref][band][pt][tokens[next][0].token];
        rate1 +=
          mb->token_costs[tx_size][type][ref][band][pt][tokens[next][1].token];
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
#if CONFIG_CODE_NONZEROCOUNT
      new_nzc0 = (best ? nzc1 : nzc0);
#endif

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
#if CONFIG_CODE_NONZEROCOUNT
        // Account for rate drop because of the nzc change.
        // TODO(debargha): Find a better solution
        if (nzc_used) {
          rate0 -= nzc_cost[nzc0] - nzc_cost[nzc0 - 1];
          rate1 -= nzc_cost[nzc1] - nzc_cost[nzc1 - 1];
        }
#endif
      } else {
        t0 = t1 = (vp9_dct_value_tokens_ptr + x)->Token;
      }
      if (next < default_eob) {
        band = get_coef_band(scan, tx_size, i + 1);
        if (t0 != DCT_EOB_TOKEN) {
          pt = trellis_get_coeff_context(scan, nb, i, t0, token_cache,
                                         pad, default_eob);
          rate0 += mb->token_costs[tx_size][type][ref][band][pt][
              tokens[next][0].token];
        }
        if (t1 != DCT_EOB_TOKEN) {
          pt = trellis_get_coeff_context(scan, nb, i, t1, token_cache,
                                         pad, default_eob);
          rate1 += mb->token_costs[tx_size][type][ref][band][pt][
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
#if CONFIG_CODE_NONZEROCOUNT
      new_nzc1 = (best ? nzc1 : nzc0) - (!x);
      nzc0 = new_nzc0;
      nzc1 = new_nzc1;
#endif
      /* Finally, make this the new head of the trellis. */
      next = i;
    }
    /* There's no choice to make for a zero coefficient, so we don't
     *  add a new trellis node, but we do need to update the costs.
     */
    else {
      band = get_coef_band(scan, tx_size, i + 1);
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
  band = get_coef_band(scan, tx_size, i + 1);
  pt = combine_entropy_contexts(*a, *l);
  rate0 = tokens[next][0].rate;
  rate1 = tokens[next][1].rate;
  error0 = tokens[next][0].error;
  error1 = tokens[next][1].error;
  t0 = tokens[next][0].token;
  t1 = tokens[next][1].token;
  rate0 += mb->token_costs[tx_size][type][ref][band][pt][t0];
  rate1 += mb->token_costs[tx_size][type][ref][band][pt][t1];
  UPDATE_RD_COST();
  best = rd_cost1 < rd_cost0;
#if CONFIG_CODE_NONZEROCOUNT
  final_nzc_exp = (best ? nzc1 : nzc0);
#endif
  final_eob = i0 - 1;
  for (i = next; i < eob; i = next) {
    x = tokens[i][best].qc;
    if (x) {
      final_eob = i;
#if CONFIG_CODE_NONZEROCOUNT
      ++final_nzc;
#endif
    }
    rc = scan[i];
    qcoeff_ptr[rc] = x;
    dqcoeff_ptr[rc] = (x * dequant_ptr[rc != 0]) / mul;

    next = tokens[i][best].next;
    best = best_index[i][best];
  }
  final_eob++;

  xd->plane[pb_idx.plane].eobs[pb_idx.block] = final_eob;
  *a = *l = (final_eob > 0);
#if CONFIG_CODE_NONZEROCOUNT
  assert(final_nzc == final_nzc_exp);
  xd->nzcs[ib] = final_nzc;
#endif
}

void vp9_optimize_sby_32x32(VP9_COMMON *const cm, MACROBLOCK *x,
                            BLOCK_SIZE_TYPE bsize) {
  const int bwl = mb_width_log2(bsize) - 1, bw = 1 << bwl;
  const int bh = 1 << (mb_height_log2(bsize) - 1);
  ENTROPY_CONTEXT ta[2], tl[2];
  int n;

  for (n = 0; n < bw; n++) {
    ENTROPY_CONTEXT *a =
        (ENTROPY_CONTEXT *) (x->e_mbd.above_context + n * 2 + 0);
    ENTROPY_CONTEXT *a1 =
        (ENTROPY_CONTEXT *) (x->e_mbd.above_context + n * 2 + 1);
    ta[n] = (a[0] + a[1] + a[2] + a[3] + a1[0] + a1[1] + a1[2] + a1[3]) != 0;
  }
  for (n = 0; n < bh; n++) {
    ENTROPY_CONTEXT *l =
        (ENTROPY_CONTEXT *) (x->e_mbd.left_context + n * 2);
    ENTROPY_CONTEXT *l1 =
        (ENTROPY_CONTEXT *) (x->e_mbd.left_context + n * 2 + 1);
    tl[n] = (l[0] + l[1] + l[2] + l[3] + l1[0] + l1[1] + l1[2] + l1[3]) != 0;
  }

  for (n = 0; n < bw * bh; n++) {
    const int x_idx = n & (bw - 1), y_idx = n >> bwl;

    optimize_b(cm, x, n * 64, PLANE_TYPE_Y_WITH_DC, x->e_mbd.block[0].dequant,
               ta + x_idx, tl + y_idx, TX_32X32, 64 * bw * bh);
  }
}

void vp9_optimize_sby_16x16(VP9_COMMON *const cm, MACROBLOCK *x,
                            BLOCK_SIZE_TYPE bsize) {
  const int bwl = mb_width_log2(bsize), bw = 1 << bwl;
  const int bh = 1 << mb_height_log2(bsize);
  ENTROPY_CONTEXT ta[4], tl[4];
  int n;

  for (n = 0; n < bw; n++) {
    ENTROPY_CONTEXT *a = (ENTROPY_CONTEXT *) (x->e_mbd.above_context + n);
    ta[n] = (a[0] + a[1] + a[2] + a[3]) != 0;
  }
  for (n = 0; n < bh; n++) {
    ENTROPY_CONTEXT *l = (ENTROPY_CONTEXT *) (x->e_mbd.left_context + n);
    tl[n] = (l[0] + l[1] + l[2] + l[3]) != 0;
  }
  for (n = 0; n < bw * bh; n++) {
    const int x_idx = n & (bw - 1), y_idx = n >> bwl;

    optimize_b(cm, x, n * 16, PLANE_TYPE_Y_WITH_DC, x->e_mbd.block[0].dequant,
               ta + x_idx, tl + y_idx, TX_16X16, 16 * bw * bh);
  }
}

void vp9_optimize_sby_8x8(VP9_COMMON *const cm, MACROBLOCK *x,
                          BLOCK_SIZE_TYPE bsize) {
  const int bwl = mb_width_log2(bsize) + 1, bw = 1 << bwl;
  const int bh = 2 << mb_height_log2(bsize);
  ENTROPY_CONTEXT ta[8], tl[8];
  int n;

  for (n = 0; n < bw; n += 2) {
    ENTROPY_CONTEXT *a =
        (ENTROPY_CONTEXT *) (x->e_mbd.above_context + (n >> 1));
    ta[n + 0] = (a[0] + a[1]) != 0;
    ta[n + 1] = (a[2] + a[3]) != 0;
  }
  for (n = 0; n < bh; n += 2) {
    ENTROPY_CONTEXT *l =
        (ENTROPY_CONTEXT *) (x->e_mbd.left_context + (n >> 1));
    tl[n + 0] = (l[0] + l[1]) != 0;
    tl[n + 1] = (l[2] + l[3]) != 0;
  }

  for (n = 0; n < bw * bh; n++) {
    const int x_idx = n & (bw - 1), y_idx = n >> bwl;

    optimize_b(cm, x, n * 4, PLANE_TYPE_Y_WITH_DC, x->e_mbd.block[0].dequant,
               ta + x_idx, tl + y_idx, TX_8X8, 4 * bw * bh);
  }
}

void vp9_optimize_sby_4x4(VP9_COMMON *const cm, MACROBLOCK *x,
                          BLOCK_SIZE_TYPE bsize) {
  int bwl = mb_width_log2(bsize), bw = 1 << bwl;
  int bh = 1 << mb_height_log2(bsize);
  ENTROPY_CONTEXT ta[16], tl[16];
  int n;

  for (n = 0; n < bw; n++)
    vpx_memcpy(&ta[n * 4], x->e_mbd.above_context + n,
               sizeof(ENTROPY_CONTEXT) * 4);
  for (n = 0; n < bh; n++)
    vpx_memcpy(&tl[n * 4], x->e_mbd.left_context + n,
               sizeof(ENTROPY_CONTEXT) * 4);
  bw *= 4;
  bh *= 4;
  bwl += 2;

  for (n = 0; n < bw * bh; n++) {
    const int x_idx = n & (bw - 1), y_idx = n >> bwl;

    optimize_b(cm, x, n, PLANE_TYPE_Y_WITH_DC, x->e_mbd.block[0].dequant,
               ta + x_idx, tl + y_idx, TX_4X4, bh * bw);
  }
}

void vp9_optimize_sbuv_32x32(VP9_COMMON *const cm, MACROBLOCK *x,
                             BLOCK_SIZE_TYPE bsize) {
  ENTROPY_CONTEXT *ta = (ENTROPY_CONTEXT *) x->e_mbd.above_context;
  ENTROPY_CONTEXT *tl = (ENTROPY_CONTEXT *) x->e_mbd.left_context;
  ENTROPY_CONTEXT *a, *l, *a1, *l1, *a2, *l2, *a3, *l3, a_ec, l_ec;
  int b;

  assert(bsize == BLOCK_SIZE_SB64X64);
  for (b = 256; b < 384; b += 64) {
    const int cidx = b >= 320 ? 20 : 16;
    a = ta + vp9_block2above_sb64[TX_32X32][b];
    l = tl + vp9_block2left_sb64[TX_32X32][b];
    a1 = a + sizeof(ENTROPY_CONTEXT_PLANES) / sizeof(ENTROPY_CONTEXT);
    l1 = l + sizeof(ENTROPY_CONTEXT_PLANES) / sizeof(ENTROPY_CONTEXT);
    a2 = a + 2 * sizeof(ENTROPY_CONTEXT_PLANES) / sizeof(ENTROPY_CONTEXT);
    l2 = l + 2 * sizeof(ENTROPY_CONTEXT_PLANES) / sizeof(ENTROPY_CONTEXT);
    a3 = a + 3 * sizeof(ENTROPY_CONTEXT_PLANES) / sizeof(ENTROPY_CONTEXT);
    l3 = l + 3 * sizeof(ENTROPY_CONTEXT_PLANES) / sizeof(ENTROPY_CONTEXT);
    a_ec = (a[0] + a[1] + a1[0] + a1[1] + a2[0] + a2[1] + a3[0] + a3[1]) != 0;
    l_ec = (l[0] + l[1] + l1[0] + l1[1] + l2[0] + l2[1] + l3[0] + l3[1]) != 0;
    optimize_b(cm, x, b, PLANE_TYPE_UV, x->e_mbd.block[cidx].dequant,
               &a_ec, &l_ec, TX_32X32, 256);
  }
}

void vp9_optimize_sbuv_16x16(VP9_COMMON *const cm, MACROBLOCK *x,
                             BLOCK_SIZE_TYPE bsize) {
  const int bwl = mb_width_log2(bsize), bhl = mb_height_log2(bsize);
  const int bw = 1 << (bwl - 1);
  const int bh = 1 << (bhl - 1);
  int uvoff = 16 << (bwl + bhl);
  ENTROPY_CONTEXT ta[2][2], tl[2][2];
  int plane, n;

  for (n = 0; n < bw; n++) {
    ENTROPY_CONTEXT_PLANES *a = x->e_mbd.above_context + n * 2;
    ENTROPY_CONTEXT_PLANES *a1 = x->e_mbd.above_context + n * 2 + 1;
    ta[0][n] = (a->u[0] + a->u[1] + a1->u[0] + a1->u[1]) != 0;
    ta[1][n] = (a->v[0] + a->v[1] + a1->v[0] + a1->v[1]) != 0;
  }
  for (n = 0; n < bh; n++) {
    ENTROPY_CONTEXT_PLANES *l = (x->e_mbd.left_context + n * 2);
    ENTROPY_CONTEXT_PLANES *l1 = (x->e_mbd.left_context + n * 2 + 1);
    tl[0][n] = (l->u[0] + l->u[1] + l1->u[0] + l1->u[1]) != 0;
    tl[1][n] = (l->v[0] + l->v[1] + l1->v[0] + l1->v[1]) != 0;
  }

  for (plane = 0; plane < 2; plane++) {
    const int cidx = 16 + plane * 4;
    for (n = 0; n < bw * bh; n++) {
      const int x_idx = n & (bw - 1), y_idx = n >> (bwl - 1);
      optimize_b(cm, x, uvoff + n * 16, PLANE_TYPE_UV,
                 x->e_mbd.block[cidx].dequant,
                 &ta[plane][x_idx], &tl[plane][y_idx],
                 TX_16X16, bh * bw * 64);
    }
    uvoff = (uvoff * 5) >> 2;  // switch u -> v
  }
}

void vp9_optimize_sbuv_8x8(VP9_COMMON *const cm, MACROBLOCK *x,
                           BLOCK_SIZE_TYPE bsize) {
  const int bwl = mb_width_log2(bsize) + 1, bhl = mb_height_log2(bsize) + 1;
  const int bw = 1 << (bwl - 1);
  const int bh = 1 << (bhl - 1);
  int uvoff = 4 << (bwl + bhl);
  ENTROPY_CONTEXT ta[2][4], tl[2][4];
  int plane, n;

  for (n = 0; n < bw; n++) {
    ENTROPY_CONTEXT_PLANES *a = x->e_mbd.above_context + n;
    ta[0][n] = (a->u[0] + a->u[1]) != 0;
    ta[1][n] = (a->v[0] + a->v[1]) != 0;
  }
  for (n = 0; n < bh; n++) {
    ENTROPY_CONTEXT_PLANES *l = x->e_mbd.left_context + n;
    tl[0][n] = (l->u[0] + l->u[1]) != 0;
    tl[1][n] = (l->v[0] + l->v[1]) != 0;
  }

  for (plane = 0; plane < 2; plane++) {
    const int cidx = 16 + plane * 4;
    for (n = 0; n < bw * bh; n++) {
      const int x_idx = n & (bw - 1), y_idx = n >> (bwl - 1);
      optimize_b(cm, x, uvoff + n * 4, PLANE_TYPE_UV,
                 x->e_mbd.block[cidx].dequant,
                 &ta[plane][x_idx], &tl[plane][y_idx],
                 TX_8X8, bh * bw * 16);
    }
    uvoff = (uvoff * 5) >> 2;  // switch u -> v
  }
}

void vp9_optimize_sbuv_4x4(VP9_COMMON *const cm, MACROBLOCK *x,
                           BLOCK_SIZE_TYPE bsize) {
  const int bwl = mb_width_log2(bsize) + 2, bhl = mb_height_log2(bsize) + 2;
  const int bw = 1 << (bwl - 1);
  const int bh = 1 << (bhl - 1);
  int uvoff = 1 << (bwl + bhl);
  ENTROPY_CONTEXT ta[2][8], tl[2][8];
  int plane, n;

  for (n = 0; n < bw; n += 2) {
    ENTROPY_CONTEXT_PLANES *a = x->e_mbd.above_context + (n >> 1);
    ta[0][n + 0] = (a->u[0]) != 0;
    ta[0][n + 1] = (a->u[1]) != 0;
    ta[1][n + 0] = (a->v[0]) != 0;
    ta[1][n + 1] = (a->v[1]) != 0;
  }
  for (n = 0; n < bh; n += 2) {
    ENTROPY_CONTEXT_PLANES *l = x->e_mbd.left_context + (n >> 1);
    tl[0][n + 0] = (l->u[0]) != 0;
    tl[0][n + 1] = (l->u[1]) != 0;
    tl[1][n + 0] = (l->v[0]) != 0;
    tl[1][n + 1] = (l->v[1]) != 0;
  }

  for (plane = 0; plane < 2; plane++) {
    const int cidx = 16 + plane * 4;
    for (n = 0; n < bw * bh; n++) {
      const int x_idx = n & (bw - 1), y_idx = n >> (bwl - 1);
      optimize_b(cm, x, uvoff + n, PLANE_TYPE_UV,
                 x->e_mbd.block[cidx].dequant,
                 &ta[plane][x_idx], &tl[plane][y_idx],
                 TX_4X4, bh * bw * 4);
    }
    uvoff = (uvoff * 5) >> 2;  // switch u -> v
  }
}

void vp9_fidct_mb(VP9_COMMON *const cm, MACROBLOCK *x) {
  MACROBLOCKD *const xd = &x->e_mbd;
  const TX_SIZE tx_size = xd->mode_info_context->mbmi.txfm_size;

  if (tx_size == TX_16X16) {
    vp9_transform_sby_16x16(x, BLOCK_SIZE_MB16X16);
    vp9_transform_sbuv_8x8(x, BLOCK_SIZE_MB16X16);
    vp9_quantize_sby_16x16(x, BLOCK_SIZE_MB16X16);
    vp9_quantize_sbuv_8x8(x, BLOCK_SIZE_MB16X16);
    if (x->optimize) {
      vp9_optimize_sby_16x16(cm, x, BLOCK_SIZE_MB16X16);
      vp9_optimize_sbuv_8x8(cm, x, BLOCK_SIZE_MB16X16);
    }
    vp9_inverse_transform_sby_16x16(xd, BLOCK_SIZE_MB16X16);
    vp9_inverse_transform_sbuv_8x8(xd, BLOCK_SIZE_MB16X16);
  } else if (tx_size == TX_8X8) {
    vp9_transform_sby_8x8(x, BLOCK_SIZE_MB16X16);
    vp9_quantize_sby_8x8(x, BLOCK_SIZE_MB16X16);
    if (x->optimize)
      vp9_optimize_sby_8x8(cm, x, BLOCK_SIZE_MB16X16);
    vp9_inverse_transform_sby_8x8(xd, BLOCK_SIZE_MB16X16);
    if (xd->mode_info_context->mbmi.mode == SPLITMV) {
      assert(xd->mode_info_context->mbmi.partitioning != PARTITIONING_4X4);
      vp9_transform_sbuv_4x4(x, BLOCK_SIZE_MB16X16);
      vp9_quantize_sbuv_4x4(x, BLOCK_SIZE_MB16X16);
      if (x->optimize)
        vp9_optimize_sbuv_4x4(cm, x, BLOCK_SIZE_MB16X16);
      vp9_inverse_transform_sbuv_4x4(xd, BLOCK_SIZE_MB16X16);
    } else {
      vp9_transform_sbuv_8x8(x, BLOCK_SIZE_MB16X16);
      vp9_quantize_sbuv_8x8(x, BLOCK_SIZE_MB16X16);
      if (x->optimize)
        vp9_optimize_sbuv_8x8(cm, x, BLOCK_SIZE_MB16X16);
      vp9_inverse_transform_sbuv_8x8(xd, BLOCK_SIZE_MB16X16);
    }
  } else {
    vp9_transform_sby_4x4(x, BLOCK_SIZE_MB16X16);
    vp9_transform_sbuv_4x4(x, BLOCK_SIZE_MB16X16);
    vp9_quantize_sby_4x4(x, BLOCK_SIZE_MB16X16);
    vp9_quantize_sbuv_4x4(x, BLOCK_SIZE_MB16X16);
    if (x->optimize) {
      vp9_optimize_sby_4x4(cm, x, BLOCK_SIZE_MB16X16);
      vp9_optimize_sbuv_4x4(cm, x, BLOCK_SIZE_MB16X16);
    }
    vp9_inverse_transform_sby_4x4(xd, BLOCK_SIZE_MB16X16);
    vp9_inverse_transform_sbuv_4x4(xd, BLOCK_SIZE_MB16X16);
  }
}

void vp9_encode_inter16x16(VP9_COMMON *const cm, MACROBLOCK *x,
                           int mb_row, int mb_col) {
  MACROBLOCKD *const xd = &x->e_mbd;

  vp9_build_inter_predictors_mb(xd, mb_row, mb_col);
  subtract_mb(x);
  vp9_fidct_mb(cm, x);
  vp9_recon_mb(xd);
}

/* this function is used by first pass only */
void vp9_encode_inter16x16y(MACROBLOCK *x, int mb_row, int mb_col) {
  MACROBLOCKD *xd = &x->e_mbd;
  BLOCK *b = &x->block[0];

  vp9_build_inter16x16_predictors_mby(xd, xd->predictor, 16, mb_row, mb_col);

  vp9_subtract_mby(x->src_diff, *(b->base_src), xd->predictor, b->src_stride);

  vp9_transform_sby_4x4(x, BLOCK_SIZE_MB16X16);
  vp9_quantize_sby_4x4(x, BLOCK_SIZE_MB16X16);
  vp9_inverse_transform_sby_4x4(xd, BLOCK_SIZE_MB16X16);

  vp9_recon_mby(xd);
}
