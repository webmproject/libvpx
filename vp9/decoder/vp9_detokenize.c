/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */


#include "vp9/common/vp9_blockd.h"
#include "vp9/common/vp9_common.h"
#include "vp9/decoder/vp9_onyxd_int.h"
#include "vpx_mem/vpx_mem.h"
#include "vpx_ports/mem.h"
#include "vp9/decoder/vp9_detokenize.h"
#include "vp9/common/vp9_seg_common.h"

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

static const vp9_prob cat6_prob[15] = {
  254, 254, 254, 252, 249, 243, 230, 196, 177, 153, 140, 133, 130, 129, 0
};

DECLARE_ALIGNED(16, extern const uint8_t, vp9_norm[256]);

#if CONFIG_CODE_ZEROGROUP
#define ZEROGROUP_ADVANCE()                \
  do {                                     \
    token_cache[scan[c]] = ZERO_TOKEN;     \
    is_last_zero[o] = 1;                   \
    c++;                                   \
  } while (0)
#define INCREMENT_COUNT(token)             \
  do {                                     \
    coef_counts[type][ref][get_coef_band(scan, txfm_size, c)] \
               [pt][token]++;     \
    token_cache[scan[c]] = token; \
    is_last_zero[o] = (token == ZERO_TOKEN);    \
  } while (0)
#else
#define INCREMENT_COUNT(token)               \
  do {                                       \
    coef_counts[type][ref][get_coef_band(scan, txfm_size, c)] \
               [pt][token]++;     \
    token_cache[scan[c]] = token; \
  } while (0)
#endif

#define WRITE_COEF_CONTINUE(val, token)                  \
  {                                                      \
    qcoeff_ptr[scan[c]] = vp9_read_and_apply_sign(r, val) * \
                            dq[c > 0] / (1 + (txfm_size == TX_32X32)); \
    INCREMENT_COUNT(token);                              \
    c++;                                                 \
    continue;                                            \
  }

#define WRITE_COEF_ONE()                                 \
{                                                        \
  qcoeff_ptr[scan[c]] = vp9_read_and_apply_sign(br, 1);  \
  INCREMENT_COUNT(ONE_TOKEN);                            \
}

#define ADJUST_COEF(prob, bits_count)  \
  do {                                 \
    if (vp9_read(r, prob))             \
      val += 1 << bits_count;          \
  } while (0);

static int decode_coefs(VP9D_COMP *dx, const MACROBLOCKD *xd,
                        vp9_reader *r, int block_idx,
                        PLANE_TYPE type, int seg_eob, int16_t *qcoeff_ptr,
                        TX_SIZE txfm_size, const int16_t *dq,
                        ENTROPY_CONTEXT *A, ENTROPY_CONTEXT *L) {
  ENTROPY_CONTEXT above_ec, left_ec;
  FRAME_CONTEXT *const fc = &dx->common.fc;
  int pt, c = 0, pad, default_eob;
  vp9_coeff_probs *coef_probs;
  vp9_prob *prob;
  vp9_coeff_count *coef_counts;
  const int ref = xd->mode_info_context->mbmi.ref_frame != INTRA_FRAME;
  TX_TYPE tx_type = DCT_DCT;
#if CONFIG_CODE_ZEROGROUP
  int is_eoo[3] = {0, 0, 0};
  int is_last_zero[3] = {0, 0, 0};
  int o, rc;
  vp9_zpc_probs *zpc_probs;
  vp9_zpc_count *zpc_count;
  vp9_prob *zprobs;
  int eoo = 0, use_eoo;
#endif
  const int *scan, *nb;
  uint8_t token_cache[1024];
#if CONFIG_CODE_ZEROGROUP
  vpx_memset(token_cache, UNKNOWN_TOKEN, sizeof(token_cache));
#endif

  switch (txfm_size) {
    default:
    case TX_4X4: {
      tx_type = (type == PLANE_TYPE_Y_WITH_DC) ?
          get_tx_type_4x4(xd, block_idx) : DCT_DCT;
      scan = get_scan_4x4(tx_type);
      above_ec = A[0] != 0;
      left_ec = L[0] != 0;
      coef_probs  = fc->coef_probs_4x4;
      coef_counts = fc->coef_counts_4x4;
      default_eob = 16;
#if CONFIG_CODE_ZEROGROUP
      zpc_probs = &(fc->zpc_probs_4x4);
      zpc_count = &(fc->zpc_counts_4x4);
#endif
      break;
    }
    case TX_8X8: {
      const BLOCK_SIZE_TYPE sb_type = xd->mode_info_context->mbmi.sb_type;
      const int sz = 1 + b_width_log2(sb_type);
      const int x = block_idx & ((1 << sz) - 1);
      const int y = block_idx - x;
      tx_type = (type == PLANE_TYPE_Y_WITH_DC) ?
          get_tx_type_8x8(xd, y + (x >> 1)) : DCT_DCT;
      scan = get_scan_8x8(tx_type);
      coef_probs  = fc->coef_probs_8x8;
      coef_counts = fc->coef_counts_8x8;
      above_ec = (A[0] + A[1]) != 0;
      left_ec = (L[0] + L[1]) != 0;
      default_eob = 64;
#if CONFIG_CODE_ZEROGROUP
      zpc_probs = &(fc->zpc_probs_8x8);
      zpc_count = &(fc->zpc_counts_8x8);
#endif
      break;
    }
    case TX_16X16: {
      const BLOCK_SIZE_TYPE sb_type = xd->mode_info_context->mbmi.sb_type;
      const int sz = 2 + b_width_log2(sb_type);
      const int x = block_idx & ((1 << sz) - 1);
      const int y = block_idx - x;
      tx_type = (type == PLANE_TYPE_Y_WITH_DC) ?
          get_tx_type_16x16(xd, y + (x >> 2)) : DCT_DCT;
      scan = get_scan_16x16(tx_type);
      coef_probs  = fc->coef_probs_16x16;
      coef_counts = fc->coef_counts_16x16;
      above_ec = (A[0] + A[1] + A[2] + A[3]) != 0;
      left_ec = (L[0] + L[1] + L[2] + L[3]) != 0;
      default_eob = 256;
#if CONFIG_CODE_ZEROGROUP
      zpc_probs = &(fc->zpc_probs_16x16);
      zpc_count = &(fc->zpc_counts_16x16);
#endif
      break;
    }
    case TX_32X32:
      scan = vp9_default_zig_zag1d_32x32;
      coef_probs = fc->coef_probs_32x32;
      coef_counts = fc->coef_counts_32x32;
      above_ec = (A[0] + A[1] + A[2] + A[3] + A[4] + A[5] + A[6] + A[7]) != 0;
      left_ec = (L[0] + L[1] + L[2] + L[3] + L[4] + L[5] + L[6] + L[7]) != 0;
      default_eob = 1024;
#if CONFIG_CODE_ZEROGROUP
      zpc_probs = &fc->zpc_probs_32x32;
      zpc_count = &fc->zpc_counts_32x32;
#endif
      break;
  }

  pt = combine_entropy_contexts(above_ec, left_ec);
  nb = vp9_get_coef_neighbors_handle(scan, &pad);

  while (1) {
    int val;
    int band;
    const uint8_t *cat6 = cat6_prob;
    if (c >= seg_eob)
      break;
    if (c)
      pt = vp9_get_coef_context(scan, nb, pad, token_cache,
                                c, default_eob);
    band = get_coef_band(scan, txfm_size, c);
    prob = coef_probs[type][ref][band][pt];
    fc->eob_branch_counts[txfm_size][type][ref][band][pt]++;
    if (!vp9_read(r, prob[EOB_CONTEXT_NODE]))
      break;
#if CONFIG_CODE_ZEROGROUP
    rc = scan[c];
    o = vp9_get_orientation(rc, txfm_size);
    if (token_cache[rc] == ZERO_TOKEN || is_eoo[o]) {
      coef_counts[type][ref][band][pt][ZERO_TOKEN]++;
      ZEROGROUP_ADVANCE();
      goto SKIP_START;
    }
#endif

SKIP_START:
    if (c >= seg_eob)
      break;
    if (c)
      pt = vp9_get_coef_context(scan, nb, pad, token_cache,
                                c, default_eob);
    band = get_coef_band(scan, txfm_size, c);
    prob = coef_probs[type][ref][band][pt];
#if CONFIG_CODE_ZEROGROUP
    rc = scan[c];
    o = vp9_get_orientation(rc, txfm_size);
    if (token_cache[rc] == ZERO_TOKEN || is_eoo[o]) {
      ZEROGROUP_ADVANCE();
      goto SKIP_START;
    }
    zprobs = (*zpc_probs)[ref]
             [coef_to_zpc_band(band)]
             [coef_to_zpc_ptok(pt)];
#endif
    if (!vp9_read(r, prob[ZERO_CONTEXT_NODE])) {
#if CONFIG_CODE_ZEROGROUP
      eoo = 0;
#if USE_ZPC_EOORIENT == 1
      use_eoo = vp9_use_eoo(c, seg_eob, scan, txfm_size, is_last_zero, is_eoo);
#else
      use_eoo = 0;
#endif
      if (use_eoo) {
        eoo = !vp9_read(r, zprobs[0]);
        ++(*zpc_count)[ref]
                      [coef_to_zpc_band(band)]
                      [coef_to_zpc_ptok(pt)][0][!eoo];
        if (eoo) {
          is_eoo[o] = 1;
        }
      }
#endif
      INCREMENT_COUNT(ZERO_TOKEN);
      ++c;
      goto SKIP_START;
    }
    // ONE_CONTEXT_NODE_0_
    if (!vp9_read(r, prob[ONE_CONTEXT_NODE])) {
      WRITE_COEF_CONTINUE(1, ONE_TOKEN);
    }
    // LOW_VAL_CONTEXT_NODE_0_
    if (!vp9_read(r, prob[LOW_VAL_CONTEXT_NODE])) {
      if (!vp9_read(r, prob[TWO_CONTEXT_NODE])) {
        WRITE_COEF_CONTINUE(2, TWO_TOKEN);
      }
      if (!vp9_read(r, prob[THREE_CONTEXT_NODE])) {
        WRITE_COEF_CONTINUE(3, THREE_TOKEN);
      }
      WRITE_COEF_CONTINUE(4, FOUR_TOKEN);
    }
    // HIGH_LOW_CONTEXT_NODE_0_
    if (!vp9_read(r, prob[HIGH_LOW_CONTEXT_NODE])) {
      if (!vp9_read(r, prob[CAT_ONE_CONTEXT_NODE])) {
        val = CAT1_MIN_VAL;
        ADJUST_COEF(CAT1_PROB0, 0);
        WRITE_COEF_CONTINUE(val, DCT_VAL_CATEGORY1);
      }
      val = CAT2_MIN_VAL;
      ADJUST_COEF(CAT2_PROB1, 1);
      ADJUST_COEF(CAT2_PROB0, 0);
      WRITE_COEF_CONTINUE(val, DCT_VAL_CATEGORY2);
    }
    // CAT_THREEFOUR_CONTEXT_NODE_0_
    if (!vp9_read(r, prob[CAT_THREEFOUR_CONTEXT_NODE])) {
      if (!vp9_read(r, prob[CAT_THREE_CONTEXT_NODE])) {
        val = CAT3_MIN_VAL;
        ADJUST_COEF(CAT3_PROB2, 2);
        ADJUST_COEF(CAT3_PROB1, 1);
        ADJUST_COEF(CAT3_PROB0, 0);
        WRITE_COEF_CONTINUE(val, DCT_VAL_CATEGORY3);
      }
      val = CAT4_MIN_VAL;
      ADJUST_COEF(CAT4_PROB3, 3);
      ADJUST_COEF(CAT4_PROB2, 2);
      ADJUST_COEF(CAT4_PROB1, 1);
      ADJUST_COEF(CAT4_PROB0, 0);
      WRITE_COEF_CONTINUE(val, DCT_VAL_CATEGORY4);
    }
    // CAT_FIVE_CONTEXT_NODE_0_:
    if (!vp9_read(r, prob[CAT_FIVE_CONTEXT_NODE])) {
      val = CAT5_MIN_VAL;
      ADJUST_COEF(CAT5_PROB4, 4);
      ADJUST_COEF(CAT5_PROB3, 3);
      ADJUST_COEF(CAT5_PROB2, 2);
      ADJUST_COEF(CAT5_PROB1, 1);
      ADJUST_COEF(CAT5_PROB0, 0);
      WRITE_COEF_CONTINUE(val, DCT_VAL_CATEGORY5);
    }
    val = 0;
    while (*cat6) {
      val = (val << 1) | vp9_read(r, *cat6++);
    }
    val += CAT6_MIN_VAL;
    WRITE_COEF_CONTINUE(val, DCT_VAL_CATEGORY6);
  }

  if (c < seg_eob)
    coef_counts[type][ref][get_coef_band(scan, txfm_size, c)]
        [pt][DCT_EOB_TOKEN]++;

  for (pt = 0; pt < (1 << txfm_size); pt++) {
    A[pt] = L[pt] = c > 0;
  }

  return c;
}

static int get_eob(MACROBLOCKD* const xd, int segment_id, int eob_max) {
  return vp9_get_segdata(xd, segment_id, SEG_LVL_SKIP) ? 0 : eob_max;
}


struct decode_block_args {
  VP9D_COMP *pbi;
  MACROBLOCKD *xd;
  vp9_reader *r;
  int *eobtotal;
};
static void decode_block(int plane, int block,
                         BLOCK_SIZE_TYPE bsize,
                         int ss_txfrm_size,
                         void *argv) {
  const struct decode_block_args* const arg = argv;
  const int bw = b_width_log2(bsize), bh = b_height_log2(bsize);
  const int old_block_idx = old_block_idx_4x4(arg->xd, bw + bh,
                                              plane, block);

  // find the maximum eob for this transform size, adjusted by segment
  const int segment_id = arg->xd->mode_info_context->mbmi.segment_id;
  const TX_SIZE ss_tx_size = ss_txfrm_size / 2;
  const int seg_eob = get_eob(arg->xd, segment_id, 16 << ss_txfrm_size);
  int16_t* const qcoeff_base = arg->xd->plane[plane].qcoeff;
  const int off = block >> ss_txfrm_size;
  const int mod = bw - ss_tx_size - arg->xd->plane[plane].subsampling_x;
  const int aoff = (off & ((1 << mod) - 1)) << ss_tx_size;
  const int loff = (off >> mod) << ss_tx_size;

  const int eob = decode_coefs(arg->pbi, arg->xd, arg->r, old_block_idx,
                               arg->xd->plane[plane].plane_type, seg_eob,
                               BLOCK_OFFSET(qcoeff_base, block, 16),
                               ss_tx_size, arg->xd->plane[plane].dequant,
                               arg->xd->plane[plane].above_context + aoff,
                               arg->xd->plane[plane].left_context + loff);

  arg->xd->plane[plane].eobs[block] = eob;
  arg->eobtotal[0] += eob;
}

int vp9_decode_tokens(VP9D_COMP* const pbi,
                         MACROBLOCKD* const xd,
                         vp9_reader *r,
                         BLOCK_SIZE_TYPE bsize) {
  int eobtotal = 0;
  struct decode_block_args args = {pbi, xd, r, &eobtotal};
  foreach_transformed_block(xd, bsize, decode_block, &args);
  return eobtotal;
}

#if CONFIG_NEWBINTRAMODES
static int decode_coefs_4x4(VP9D_COMP *dx, MACROBLOCKD *xd,
                            vp9_reader *r,
                            PLANE_TYPE type, int i, int seg_eob) {
  const struct plane_block_idx pb_idx = plane_block_idx(16, i);
  const int mod = 2 - xd->plane[pb_idx.plane].subsampling_x;
  const int aoff = pb_idx.block & ((1 << mod) - 1);
  const int loff = pb_idx.block >> mod;
  ENTROPY_CONTEXT *A = xd->plane[pb_idx.plane].above_context;
  ENTROPY_CONTEXT *L = xd->plane[pb_idx.plane].left_context;
  const int c = decode_coefs(dx, xd, r, i, type, seg_eob,
      BLOCK_OFFSET(xd->plane[pb_idx.plane].qcoeff, pb_idx.block, 16), TX_4X4,
      xd->plane[pb_idx.plane].dequant, A + aoff, L + loff);
  xd->plane[pb_idx.plane].eobs[pb_idx.block] = c;
  return c;
}

static int decode_mb_tokens_4x4_uv(VP9D_COMP* const dx,
                                   MACROBLOCKD* const xd,
                                   vp9_reader *r,
                                   int seg_eob) {
  int i, eobtotal = 0;

  // chroma blocks
  for (i = 16; i < 24; i++)
    eobtotal += decode_coefs_4x4(dx, xd, r, PLANE_TYPE_UV, i, seg_eob);

  return eobtotal;
}

int vp9_decode_mb_tokens_4x4_uv(VP9D_COMP* const dx,
                                MACROBLOCKD* const xd,
                                vp9_reader *r) {
  const int segment_id = xd->mode_info_context->mbmi.segment_id;
  const int seg_eob = get_eob(xd, segment_id, 16);

  return decode_mb_tokens_4x4_uv(dx, xd, r, seg_eob);
}

int vp9_decode_coefs_4x4(VP9D_COMP *dx, MACROBLOCKD *xd,
                         vp9_reader *r,
                         PLANE_TYPE type, int i) {
  const int segment_id = xd->mode_info_context->mbmi.segment_id;
  const int seg_eob = get_eob(xd, segment_id, 16);
  return decode_coefs_4x4(dx, xd, r, type, i, seg_eob);
}
#endif
