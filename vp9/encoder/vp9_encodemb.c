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

void vp9_subtract_block(int rows, int cols,
                        int16_t *diff_ptr, int diff_stride,
                        const uint8_t *src_ptr, int src_stride,
                        const uint8_t *pred_ptr, int pred_stride) {
  int r, c;

  for (r = 0; r < rows; r++) {
    for (c = 0; c < cols; c++)
      diff_ptr[c] = src_ptr[c] - pred_ptr[c];

    diff_ptr += diff_stride;
    pred_ptr += pred_stride;
    src_ptr  += src_stride;
  }
}


static void subtract_plane(MACROBLOCK *x, BLOCK_SIZE_TYPE bsize, int plane) {
  const MACROBLOCKD * const xd = &x->e_mbd;
  const int bw = 4 << (b_width_log2(bsize) - xd->plane[plane].subsampling_x);
  const int bh = 4 << (b_height_log2(bsize) - xd->plane[plane].subsampling_y);
  const uint8_t *src = x->plane[plane].src.buf;
  const int src_stride = x->plane[plane].src.stride;

  assert(plane < 3);
  vp9_subtract_block(bh, bw,
                     x->plane[plane].src_diff, bw, src, src_stride,
                     xd->plane[plane].dst.buf, xd->plane[plane].dst.stride);
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


void vp9_transform_sby_32x32(MACROBLOCK *x, BLOCK_SIZE_TYPE bsize) {
  const int bwl = b_width_log2(bsize) - 3, bw = 1 << bwl;
  const int bh = 1 << (b_height_log2(bsize) - 3);
  const int stride = 32 << bwl;
  int n;

  for (n = 0; n < bw * bh; n++) {
    const int x_idx = n & (bw - 1), y_idx = n >> bwl;

    vp9_short_fdct32x32(x->plane[0].src_diff + y_idx * stride * 32 + x_idx * 32,
                        x->plane[0].coeff + n * 1024, stride * 2);
  }
}

void vp9_transform_sby_16x16(MACROBLOCK *x, BLOCK_SIZE_TYPE bsize) {
  const int bwl = b_width_log2(bsize) - 2, bw = 1 << bwl;
  const int bh = 1 << (b_height_log2(bsize) - 2);
  const int stride = 16 << bwl, bstride = 4 << bwl;
  MACROBLOCKD *const xd = &x->e_mbd;
  int n;

  for (n = 0; n < bw * bh; n++) {
    const int x_idx = n & (bw - 1), y_idx = n >> bwl;
    const TX_TYPE tx_type = get_tx_type_16x16(xd,
                                              (y_idx * bstride + x_idx) * 4);

    if (tx_type != DCT_DCT) {
      vp9_short_fht16x16(x->plane[0].src_diff +
                             y_idx * stride * 16 + x_idx * 16,
                         x->plane[0].coeff + n * 256, stride, tx_type);
    } else {
      x->fwd_txm16x16(x->plane[0].src_diff + y_idx * stride * 16 + x_idx * 16,
                      x->plane[0].coeff + n * 256, stride * 2);
    }
  }
}

void vp9_transform_sby_8x8(MACROBLOCK *x, BLOCK_SIZE_TYPE bsize) {
  const int bwl = b_width_log2(bsize) - 1, bw = 1 << bwl;
  const int bh = 1 << (b_height_log2(bsize) - 1);
  const int stride = 8 << bwl, bstride = 2 << bwl;
  MACROBLOCKD *const xd = &x->e_mbd;
  int n;

  for (n = 0; n < bw * bh; n++) {
    const int x_idx = n & (bw - 1), y_idx = n >> bwl;
    const TX_TYPE tx_type = get_tx_type_8x8(xd, (y_idx * bstride + x_idx) * 2);

    if (tx_type != DCT_DCT) {
      vp9_short_fht8x8(x->plane[0].src_diff + y_idx * stride * 8 + x_idx * 8,
                       x->plane[0].coeff + n * 64, stride, tx_type);
    } else {
      x->fwd_txm8x8(x->plane[0].src_diff + y_idx * stride * 8 + x_idx * 8,
                    x->plane[0].coeff + n * 64, stride * 2);
    }
  }
}

void vp9_transform_sby_4x4(MACROBLOCK *x, BLOCK_SIZE_TYPE bsize) {
  const int bwl = b_width_log2(bsize), bw = 1 << bwl;
  const int bh = 1 << b_height_log2(bsize);
  const int stride = 4 << bwl;
  MACROBLOCKD *const xd = &x->e_mbd;
  int n;

  for (n = 0; n < bw * bh; n++) {
    const int x_idx = n & (bw - 1), y_idx = n >> bwl;
    const TX_TYPE tx_type = get_tx_type_4x4(xd, n);

    if (tx_type != DCT_DCT) {
      vp9_short_fht4x4(x->plane[0].src_diff + y_idx * stride * 4 + x_idx * 4,
                       x->plane[0].coeff + n * 16, stride, tx_type);
    } else {
      x->fwd_txm4x4(x->plane[0].src_diff + y_idx * stride * 4 + x_idx * 4,
                    x->plane[0].coeff + n * 16, stride * 2);
    }
  }
}

void vp9_transform_sbuv_32x32(MACROBLOCK *x, BLOCK_SIZE_TYPE bsize) {
  assert(bsize == BLOCK_SIZE_SB64X64);
  vp9_clear_system_state();
  vp9_short_fdct32x32(x->plane[1].src_diff, x->plane[1].coeff, 64);
  vp9_short_fdct32x32(x->plane[2].src_diff, x->plane[2].coeff, 64);
}

void vp9_transform_sbuv_16x16(MACROBLOCK *x, BLOCK_SIZE_TYPE bsize) {
  const int bwl = b_width_log2(bsize) - 2, bhl = b_height_log2(bsize) - 2;
  const int bw = 1 << (bwl - 1), bh = 1 << (bhl - 1);
  const int stride = 16 << (bwl - 1);
  int n;

  vp9_clear_system_state();
  for (n = 0; n < bw * bh; n++) {
    const int x_idx = n & (bw - 1), y_idx = n >> (bwl - 1);

    x->fwd_txm16x16(x->plane[1].src_diff + y_idx * stride * 16 + x_idx * 16,
                    x->plane[1].coeff + n * 256, stride * 2);
    x->fwd_txm16x16(x->plane[2].src_diff + y_idx * stride * 16 + x_idx * 16,
                    x->plane[2].coeff + n * 256, stride * 2);
  }
}

void vp9_transform_sbuv_8x8(MACROBLOCK *x, BLOCK_SIZE_TYPE bsize) {
  const int bwl = b_width_log2(bsize) - 1, bhl = b_height_log2(bsize) - 1;
  const int bw = 1 << (bwl - 1), bh = 1 << (bhl - 1);
  const int stride = 8 << (bwl - 1);
  int n;

  vp9_clear_system_state();
  for (n = 0; n < bw * bh; n++) {
    const int x_idx = n & (bw - 1), y_idx = n >> (bwl - 1);

    x->fwd_txm8x8(x->plane[1].src_diff + y_idx * stride * 8 + x_idx * 8,
                  x->plane[1].coeff + n * 64, stride * 2);
    x->fwd_txm8x8(x->plane[2].src_diff + y_idx * stride * 8 + x_idx * 8,
                  x->plane[2].coeff + n * 64, stride * 2);
  }
}

void vp9_transform_sbuv_4x4(MACROBLOCK *x, BLOCK_SIZE_TYPE bsize) {
  const int bwl = b_width_log2(bsize), bhl = b_height_log2(bsize);
  const int bw = 1 << (bwl - 1), bh = 1 << (bhl - 1);
  const int stride = 4 << (bwl - 1);
  int n;

  vp9_clear_system_state();
  for (n = 0; n < bw * bh; n++) {
    const int x_idx = n & (bw - 1), y_idx = n >> (bwl - 1);

    x->fwd_txm4x4(x->plane[1].src_diff + y_idx * stride * 4 + x_idx * 4,
                  x->plane[1].coeff + n * 16, stride * 2);
    x->fwd_txm4x4(x->plane[2].src_diff + y_idx * stride * 4 + x_idx * 4,
                  x->plane[2].coeff + n * 16, stride * 2);
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
  int bak = token_cache[scan[idx]], pt;
  token_cache[scan[idx]] = token;
  pt = vp9_get_coef_context(scan, nb, pad, token_cache, idx + 1, l);
  token_cache[scan[idx]] = bak;
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
  const int16_t *coeff_ptr = BLOCK_OFFSET(mb->plane[pb_idx.plane].coeff,
                                          pb_idx.block, 16);
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

  assert((!type && !pb_idx.plane) || (type && pb_idx.plane));
  dqcoeff_ptr = BLOCK_OFFSET(xd->plane[pb_idx.plane].dqcoeff, pb_idx.block, 16);
  qcoeff_ptr = BLOCK_OFFSET(xd->plane[pb_idx.plane].qcoeff, pb_idx.block, 16);
  switch (tx_size) {
    default:
    case TX_4X4: {
      const TX_TYPE tx_type = get_tx_type_4x4(xd, ib);
      default_eob = 16;
      scan = get_scan_4x4(tx_type);
      break;
    }
    case TX_8X8: {
      const BLOCK_SIZE_TYPE sb_type = xd->mode_info_context->mbmi.sb_type;
      const int sz = 1 + b_width_log2(sb_type);
      const int x = ib & ((1 << sz) - 1), y = ib - x;
      const TX_TYPE tx_type = get_tx_type_8x8(xd, y + (x >> 1));
      scan = get_scan_8x8(tx_type);
      default_eob = 64;
      break;
    }
    case TX_16X16: {
      const BLOCK_SIZE_TYPE sb_type = xd->mode_info_context->mbmi.sb_type;
      const int sz = 2 + b_width_log2(sb_type);
      const int x = ib & ((1 << sz) - 1), y = ib - x;
      const TX_TYPE tx_type = get_tx_type_16x16(xd, y + (x >> 2));
      scan = get_scan_16x16(tx_type);
      default_eob = 256;
      break;
    }
    case TX_32X32:
      scan = vp9_default_zig_zag1d_32x32;
      default_eob = 1024;
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
  tokens[eob][0].rate = 0;
  tokens[eob][0].error = 0;
  tokens[eob][0].next = default_eob;
  tokens[eob][0].token = DCT_EOB_TOKEN;
  tokens[eob][0].qc = 0;
  *(tokens[eob] + 1) = *(tokens[eob] + 0);
  next = eob;
  for (i = 0; i < eob; i++)
    token_cache[scan[i]] = vp9_dct_value_tokens_ptr[qcoeff_ptr[scan[i]]].token;
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

  xd->plane[pb_idx.plane].eobs[pb_idx.block] = final_eob;
  *a = *l = (final_eob > 0);
}

void vp9_optimize_sby_32x32(VP9_COMMON *const cm, MACROBLOCK *x,
                            BLOCK_SIZE_TYPE bsize) {
  MACROBLOCKD *const xd = &x->e_mbd;
  ENTROPY_CONTEXT *a = xd->plane[0].above_context;
  ENTROPY_CONTEXT *l = xd->plane[0].left_context;
  const int bwl = b_width_log2(bsize) - 3, bw = 1 << bwl;
  const int bh = 1 << (b_height_log2(bsize) - 3);
  ENTROPY_CONTEXT ta[2], tl[2];
  int n;

  for (n = 0; n < bw; n++, a += 8)
    ta[n] = (a[0] + a[1] + a[2] + a[3] + a[4] + a[5] + a[6] + a[7]) != 0;
  for (n = 0; n < bh; n++, l += 8)
    tl[n] = (l[0] + l[1] + l[2] + l[3] + l[4] + l[5] + l[6] + l[7]) != 0;

  for (n = 0; n < bw * bh; n++) {
    const int x_idx = n & (bw - 1), y_idx = n >> bwl;

    optimize_b(cm, x, n * 64, PLANE_TYPE_Y_WITH_DC, x->e_mbd.plane[0].dequant,
               ta + x_idx, tl + y_idx, TX_32X32, 64 * bw * bh);
  }
}

void vp9_optimize_sby_16x16(VP9_COMMON *const cm, MACROBLOCK *x,
                            BLOCK_SIZE_TYPE bsize) {
  MACROBLOCKD *const xd = &x->e_mbd;
  ENTROPY_CONTEXT *a = xd->plane[0].above_context;
  ENTROPY_CONTEXT *l = xd->plane[0].left_context;
  const int bwl = b_width_log2(bsize) - 2, bw = 1 << bwl;
  const int bh = 1 << (b_height_log2(bsize) - 2);
  ENTROPY_CONTEXT ta[4], tl[4];
  int n;

  for (n = 0; n < bw; n++, a += 4)
    ta[n] = (a[0] + a[1] + a[2] + a[3]) != 0;
  for (n = 0; n < bh; n++, l += 4)
    tl[n] = (l[0] + l[1] + l[2] + l[3]) != 0;

  for (n = 0; n < bw * bh; n++) {
    const int x_idx = n & (bw - 1), y_idx = n >> bwl;

    optimize_b(cm, x, n * 16, PLANE_TYPE_Y_WITH_DC, x->e_mbd.plane[0].dequant,
               ta + x_idx, tl + y_idx, TX_16X16, 16 * bw * bh);
  }
}

void vp9_optimize_sby_8x8(VP9_COMMON *const cm, MACROBLOCK *x,
                          BLOCK_SIZE_TYPE bsize) {
  MACROBLOCKD *const xd = &x->e_mbd;
  ENTROPY_CONTEXT *a = xd->plane[0].above_context;
  ENTROPY_CONTEXT *l = xd->plane[0].left_context;
  const int bwl = b_width_log2(bsize) - 1, bw = 1 << bwl;
  const int bh = 1 << (b_height_log2(bsize) - 1);
  ENTROPY_CONTEXT ta[8], tl[8];
  int n;

  for (n = 0; n < bw; n++, a += 2)
    ta[n] = (a[0] + a[1]) != 0;
  for (n = 0; n < bh; n++, l += 2)
    tl[n] = (l[0] + l[1]) != 0;

  for (n = 0; n < bw * bh; n++) {
    const int x_idx = n & (bw - 1), y_idx = n >> bwl;

    optimize_b(cm, x, n * 4, PLANE_TYPE_Y_WITH_DC, x->e_mbd.plane[0].dequant,
               ta + x_idx, tl + y_idx, TX_8X8, 4 * bw * bh);
  }
}

void vp9_optimize_sby_4x4(VP9_COMMON *const cm, MACROBLOCK *x,
                          BLOCK_SIZE_TYPE bsize) {
  MACROBLOCKD *const xd = &x->e_mbd;
  int bwl = b_width_log2(bsize), bw = 1 << bwl;
  int bh = 1 << b_height_log2(bsize);
  ENTROPY_CONTEXT ta[16], tl[16];
  int n;

  vpx_memcpy(ta, xd->plane[0].above_context, sizeof(ENTROPY_CONTEXT) * bw);
  vpx_memcpy(tl, xd->plane[0].left_context, sizeof(ENTROPY_CONTEXT) * bh);

  for (n = 0; n < bw * bh; n++) {
    const int x_idx = n & (bw - 1), y_idx = n >> bwl;

    optimize_b(cm, x, n, PLANE_TYPE_Y_WITH_DC, x->e_mbd.plane[0].dequant,
               ta + x_idx, tl + y_idx, TX_4X4, bh * bw);
  }
}

void vp9_optimize_sbuv_32x32(VP9_COMMON *const cm, MACROBLOCK *x,
                             BLOCK_SIZE_TYPE bsize) {
  MACROBLOCKD *const xd = &x->e_mbd;
  int b;

  assert(bsize == BLOCK_SIZE_SB64X64);
  for (b = 256; b < 384; b += 64) {
    const int plane = 1 + (b >= 320);
    ENTROPY_CONTEXT *a = xd->plane[plane].above_context;
    ENTROPY_CONTEXT *l = xd->plane[plane].left_context;
    ENTROPY_CONTEXT a_ec, l_ec;

    a_ec = (a[0] + a[1] + a[2] + a[3] + a[4] + a[5] + a[6] + a[7]) != 0;
    l_ec = (l[0] + l[1] + l[2] + l[3] + l[4] + l[5] + l[6] + l[7]) != 0;
    optimize_b(cm, x, b, PLANE_TYPE_UV, x->e_mbd.plane[plane].dequant,
               &a_ec, &l_ec, TX_32X32, 256);
  }
}

void vp9_optimize_sbuv_16x16(VP9_COMMON *const cm, MACROBLOCK *x,
                             BLOCK_SIZE_TYPE bsize) {
  MACROBLOCKD *const xd = &x->e_mbd;
  const int bwl = b_width_log2(bsize) - 2, bhl = b_height_log2(bsize) - 2;
  const int bw = 1 << (bwl - 1);
  const int bh = 1 << (bhl - 1);
  int uvoff = 16 << (bwl + bhl);
  int plane, n;

  for (plane = 1; plane < MAX_MB_PLANE; plane++) {
    ENTROPY_CONTEXT ta[2], *a = xd->plane[plane].above_context;
    ENTROPY_CONTEXT tl[2], *l = xd->plane[plane].left_context;

    for (n = 0; n < bw; n++, a += 4)
      ta[n] = (a[0] + a[1] + a[2] + a[3]) != 0;
    for (n = 0; n < bh; n++, l += 4)
      tl[n] = (l[0] + l[1] + l[2] + l[3]) != 0;

    for (n = 0; n < bw * bh; n++) {
      const int x_idx = n & (bw - 1), y_idx = n >> (bwl - 1);
      optimize_b(cm, x, uvoff + n * 16, PLANE_TYPE_UV,
                 x->e_mbd.plane[plane].dequant,
                 &ta[x_idx], &tl[y_idx],
                 TX_16X16, bh * bw * 64);
    }
    uvoff = (uvoff * 5) >> 2;  // switch u -> v
  }
}

void vp9_optimize_sbuv_8x8(VP9_COMMON *const cm, MACROBLOCK *x,
                           BLOCK_SIZE_TYPE bsize) {
  MACROBLOCKD *const xd = &x->e_mbd;
  const int bwl = b_width_log2(bsize) - 1, bhl = b_height_log2(bsize) - 1;
  const int bw = 1 << (bwl - 1);
  const int bh = 1 << (bhl - 1);
  int uvoff = 4 << (bwl + bhl);
  int plane, n;

  for (plane = 1; plane < MAX_MB_PLANE; plane++) {
    ENTROPY_CONTEXT ta[4], *a = xd->plane[plane].above_context;
    ENTROPY_CONTEXT tl[4], *l = xd->plane[plane].left_context;

    for (n = 0; n < bw; n++, a += 2)
      ta[n] = (a[0] + a[1]) != 0;
    for (n = 0; n < bh; n++, l += 2)
      tl[n] = (l[0] + l[1]) != 0;

    for (n = 0; n < bw * bh; n++) {
      const int x_idx = n & (bw - 1), y_idx = n >> (bwl - 1);
      optimize_b(cm, x, uvoff + n * 4, PLANE_TYPE_UV,
                 x->e_mbd.plane[plane].dequant,
                 &ta[x_idx], &tl[y_idx],
                 TX_8X8, bh * bw * 16);
    }
    uvoff = (uvoff * 5) >> 2;  // switch u -> v
  }
}

void vp9_optimize_sbuv_4x4(VP9_COMMON *const cm, MACROBLOCK *x,
                           BLOCK_SIZE_TYPE bsize) {
  MACROBLOCKD *const xd = &x->e_mbd;
  const int bwl = b_width_log2(bsize), bhl = b_height_log2(bsize);
  const int bw = 1 << (bwl - 1);
  const int bh = 1 << (bhl - 1);
  int uvoff = 1 << (bwl + bhl);
  int plane, n;

  for (plane = 1; plane < MAX_MB_PLANE; plane++) {
    ENTROPY_CONTEXT ta[8], tl[8];

    vpx_memcpy(ta, xd->plane[plane].above_context,
               sizeof(ENTROPY_CONTEXT) * bw);
    vpx_memcpy(tl, xd->plane[plane].left_context,
               sizeof(ENTROPY_CONTEXT) * bh);

    for (n = 0; n < bw * bh; n++) {
      const int x_idx = n & (bw - 1), y_idx = n >> (bwl - 1);
      optimize_b(cm, x, uvoff + n, PLANE_TYPE_UV,
                 x->e_mbd.plane[plane].dequant,
                 &ta[x_idx], &tl[y_idx],
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
                           int mi_row, int mi_col) {
  MACROBLOCKD *const xd = &x->e_mbd;

  vp9_build_inter_predictors_sb(xd, mi_row, mi_col, BLOCK_SIZE_MB16X16);
  vp9_subtract_sb(x, BLOCK_SIZE_MB16X16);
  vp9_fidct_mb(cm, x);
  vp9_recon_sb(xd, BLOCK_SIZE_MB16X16);
}

/* this function is used by first pass only */
void vp9_encode_inter16x16y(MACROBLOCK *x, int mi_row, int mi_col) {
  MACROBLOCKD *xd = &x->e_mbd;

  vp9_build_inter_predictors_sby(xd, mi_row, mi_col, BLOCK_SIZE_MB16X16);
  vp9_subtract_sby(x, BLOCK_SIZE_MB16X16);

  vp9_transform_sby_4x4(x, BLOCK_SIZE_MB16X16);
  vp9_quantize_sby_4x4(x, BLOCK_SIZE_MB16X16);
  vp9_inverse_transform_sby_4x4(xd, BLOCK_SIZE_MB16X16);

  vp9_recon_sby(xd, BLOCK_SIZE_MB16X16);
}
