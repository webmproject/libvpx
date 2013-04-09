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

static int16_t get_signed(BOOL_DECODER *br, int16_t value_to_sign) {
  return decode_bool(br, 128) ? -value_to_sign : value_to_sign;
}


#define INCREMENT_COUNT(token)               \
  do {                                       \
    coef_counts[type][ref][get_coef_band(scan, txfm_size, c)] \
               [pt][token]++;     \
    token_cache[c] = token; \
    pt = vp9_get_coef_context(scan, nb, pad, token_cache,     \
                              c + 1, default_eob); \
  } while (0)

#if CONFIG_CODE_NONZEROCOUNT
#define WRITE_COEF_CONTINUE(val, token)                       \
  {                                                           \
    qcoeff_ptr[scan[c]] = get_signed(br, val);                \
    INCREMENT_COUNT(token);                                   \
    c++;                                                      \
    nzc++;                                                    \
    continue;                                                 \
  }
#else
#define WRITE_COEF_CONTINUE(val, token)                  \
  {                                                      \
    qcoeff_ptr[scan[c]] = get_signed(br, val);           \
    INCREMENT_COUNT(token);                              \
    c++;                                                 \
    continue;                                            \
  }
#endif  // CONFIG_CODE_NONZEROCOUNT

#define ADJUST_COEF(prob, bits_count)  \
  do {                                 \
    if (vp9_read(br, prob))            \
      val += 1 << bits_count;          \
  } while (0);

static int decode_coefs(VP9D_COMP *dx, const MACROBLOCKD *xd,
                        BOOL_DECODER* const br, int block_idx,
                        PLANE_TYPE type, int seg_eob, int16_t *qcoeff_ptr,
                        TX_SIZE txfm_size) {
  ENTROPY_CONTEXT* const A0 = (ENTROPY_CONTEXT *) xd->above_context;
  ENTROPY_CONTEXT* const L0 = (ENTROPY_CONTEXT *) xd->left_context;
  int aidx, lidx;
  ENTROPY_CONTEXT above_ec, left_ec;
  FRAME_CONTEXT *const fc = &dx->common.fc;
  int pt, c = 0, pad, default_eob;
  vp9_coeff_probs *coef_probs;
  vp9_prob *prob;
  vp9_coeff_count *coef_counts;
  const int ref = xd->mode_info_context->mbmi.ref_frame != INTRA_FRAME;
#if CONFIG_CODE_NONZEROCOUNT
  const int nzc_used = get_nzc_used(txfm_size);
  uint16_t nzc = 0;
  uint16_t nzc_expected =
      nzc_used ? xd->mode_info_context->mbmi.nzcs[block_idx] : 0;
#endif
  const int *scan, *nb;
  uint8_t token_cache[1024];

  if (xd->mode_info_context->mbmi.sb_type == BLOCK_SIZE_SB64X64) {
    aidx = vp9_block2above_sb64[txfm_size][block_idx];
    lidx = vp9_block2left_sb64[txfm_size][block_idx];
  } else if (xd->mode_info_context->mbmi.sb_type == BLOCK_SIZE_SB32X32) {
    aidx = vp9_block2above_sb[txfm_size][block_idx];
    lidx = vp9_block2left_sb[txfm_size][block_idx];
  } else {
    aidx = vp9_block2above[txfm_size][block_idx];
    lidx = vp9_block2left[txfm_size][block_idx];
  }

  switch (txfm_size) {
    default:
    case TX_4X4: {
      const TX_TYPE tx_type = (type == PLANE_TYPE_Y_WITH_DC) ?
                              get_tx_type_4x4(xd, block_idx) : DCT_DCT;
      switch (tx_type) {
        default:
          scan = vp9_default_zig_zag1d_4x4;
          break;
        case ADST_DCT:
          scan = vp9_row_scan_4x4;
          break;
        case DCT_ADST:
          scan = vp9_col_scan_4x4;
          break;
      }
      above_ec = A0[aidx] != 0;
      left_ec = L0[lidx] != 0;
      coef_probs  = fc->coef_probs_4x4;
      coef_counts = fc->coef_counts_4x4;
      default_eob = 16;
      break;
    }
    case TX_8X8: {
      const BLOCK_SIZE_TYPE sb_type = xd->mode_info_context->mbmi.sb_type;
      const int sz = 3 + sb_type, x = block_idx & ((1 << sz) - 1);
      const int y = block_idx - x;
      const TX_TYPE tx_type = (type == PLANE_TYPE_Y_WITH_DC) ?
                              get_tx_type_8x8(xd, y + (x >> 1)) : DCT_DCT;
      switch (tx_type) {
        default:
          scan = vp9_default_zig_zag1d_8x8;
          break;
        case ADST_DCT:
          scan = vp9_row_scan_8x8;
          break;
        case DCT_ADST:
          scan = vp9_col_scan_8x8;
          break;
      }
      coef_probs  = fc->coef_probs_8x8;
      coef_counts = fc->coef_counts_8x8;
      above_ec = (A0[aidx] + A0[aidx + 1]) != 0;
      left_ec  = (L0[lidx] + L0[lidx + 1]) != 0;
      default_eob = 64;
      break;
    }
    case TX_16X16: {
      const BLOCK_SIZE_TYPE sb_type = xd->mode_info_context->mbmi.sb_type;
      const int sz = 4 + sb_type, x = block_idx & ((1 << sz) - 1);
      const int y = block_idx - x;
      const TX_TYPE tx_type = (type == PLANE_TYPE_Y_WITH_DC) ?
                              get_tx_type_16x16(xd, y + (x >> 2)) : DCT_DCT;
      switch (tx_type) {
        default:
          scan = vp9_default_zig_zag1d_16x16;
          break;
        case ADST_DCT:
          scan = vp9_row_scan_16x16;
          break;
        case DCT_ADST:
          scan = vp9_col_scan_16x16;
          break;
      }
      coef_probs  = fc->coef_probs_16x16;
      coef_counts = fc->coef_counts_16x16;
      if (type == PLANE_TYPE_UV) {
        ENTROPY_CONTEXT *A1 = (ENTROPY_CONTEXT *) (xd->above_context + 1);
        ENTROPY_CONTEXT *L1 = (ENTROPY_CONTEXT *) (xd->left_context + 1);
        above_ec = (A0[aidx] + A0[aidx + 1] + A1[aidx] + A1[aidx + 1]) != 0;
        left_ec  = (L0[lidx] + L0[lidx + 1] + L1[lidx] + L1[lidx + 1]) != 0;
      } else {
        above_ec = (A0[aidx] + A0[aidx + 1] + A0[aidx + 2] + A0[aidx + 3]) != 0;
        left_ec  = (L0[lidx] + L0[lidx + 1] + L0[lidx + 2] + L0[lidx + 3]) != 0;
      }
      default_eob = 256;
      break;
    }
    case TX_32X32:
      scan = vp9_default_zig_zag1d_32x32;
      coef_probs = fc->coef_probs_32x32;
      coef_counts = fc->coef_counts_32x32;
      if (type == PLANE_TYPE_UV) {
        ENTROPY_CONTEXT *A1 = (ENTROPY_CONTEXT *) (xd->above_context + 1);
        ENTROPY_CONTEXT *L1 = (ENTROPY_CONTEXT *) (xd->left_context + 1);
        ENTROPY_CONTEXT *A2 = (ENTROPY_CONTEXT *) (xd->above_context + 2);
        ENTROPY_CONTEXT *L2 = (ENTROPY_CONTEXT *) (xd->left_context + 2);
        ENTROPY_CONTEXT *A3 = (ENTROPY_CONTEXT *) (xd->above_context + 3);
        ENTROPY_CONTEXT *L3 = (ENTROPY_CONTEXT *) (xd->left_context + 3);
        above_ec = (A0[aidx] + A0[aidx + 1] + A1[aidx] + A1[aidx + 1] +
                    A2[aidx] + A2[aidx + 1] + A3[aidx] + A3[aidx + 1]) != 0;
        left_ec  = (L0[lidx] + L0[lidx + 1] + L1[lidx] + L1[lidx + 1] +
                    L2[lidx] + L2[lidx + 1] + L3[lidx] + L3[lidx + 1]) != 0;
      } else {
        ENTROPY_CONTEXT *A1 = (ENTROPY_CONTEXT *) (xd->above_context + 1);
        ENTROPY_CONTEXT *L1 = (ENTROPY_CONTEXT *) (xd->left_context + 1);
        above_ec = (A0[aidx] + A0[aidx + 1] + A0[aidx + 2] + A0[aidx + 3] +
                    A1[aidx] + A1[aidx + 1] + A1[aidx + 2] + A1[aidx + 3]) != 0;
        left_ec  = (L0[lidx] + L0[lidx + 1] + L0[lidx + 2] + L0[lidx + 3] +
                    L1[lidx] + L1[lidx + 1] + L1[lidx + 2] + L1[lidx + 3]) != 0;
      }
      default_eob = 1024;
      break;
  }

  VP9_COMBINEENTROPYCONTEXTS(pt, above_ec, left_ec);
  nb = vp9_get_coef_neighbors_handle(scan, &pad);

  while (1) {
    int val;
    const uint8_t *cat6 = cat6_prob;

    if (c >= seg_eob)
      break;
#if CONFIG_CODE_NONZEROCOUNT
    if (nzc_used && nzc == nzc_expected)
      break;
#endif
    prob = coef_probs[type][ref][get_coef_band(scan, txfm_size, c)][pt];
    fc->eob_branch_counts[txfm_size][type][ref]
                         [get_coef_band(scan, txfm_size, c)][pt]++;
#if CONFIG_CODE_NONZEROCOUNT
    if (!nzc_used)
#endif
      if (!vp9_read(br, prob[EOB_CONTEXT_NODE]))
        break;
SKIP_START:
    if (c >= seg_eob)
      break;
#if CONFIG_CODE_NONZEROCOUNT
    if (nzc_used && nzc == nzc_expected)
      break;
    // decode zero node only if there are zeros left
    if (!nzc_used || seg_eob - nzc_expected - c + nzc > 0)
#endif
    if (!vp9_read(br, prob[ZERO_CONTEXT_NODE])) {
      INCREMENT_COUNT(ZERO_TOKEN);
      ++c;
      prob = coef_probs[type][ref][get_coef_band(scan, txfm_size, c)][pt];
      goto SKIP_START;
    }
    // ONE_CONTEXT_NODE_0_
    if (!vp9_read(br, prob[ONE_CONTEXT_NODE])) {
      WRITE_COEF_CONTINUE(1, ONE_TOKEN);
    }
    // LOW_VAL_CONTEXT_NODE_0_
    if (!vp9_read(br, prob[LOW_VAL_CONTEXT_NODE])) {
      if (!vp9_read(br, prob[TWO_CONTEXT_NODE])) {
        WRITE_COEF_CONTINUE(2, TWO_TOKEN);
      }
      if (!vp9_read(br, prob[THREE_CONTEXT_NODE])) {
        WRITE_COEF_CONTINUE(3, THREE_TOKEN);
      }
      WRITE_COEF_CONTINUE(4, FOUR_TOKEN);
    }
    // HIGH_LOW_CONTEXT_NODE_0_
    if (!vp9_read(br, prob[HIGH_LOW_CONTEXT_NODE])) {
      if (!vp9_read(br, prob[CAT_ONE_CONTEXT_NODE])) {
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
    if (!vp9_read(br, prob[CAT_THREEFOUR_CONTEXT_NODE])) {
      if (!vp9_read(br, prob[CAT_THREE_CONTEXT_NODE])) {
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
    if (!vp9_read(br, prob[CAT_FIVE_CONTEXT_NODE])) {
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
      val = (val << 1) | vp9_read(br, *cat6++);
    }
    val += CAT6_MIN_VAL;
    WRITE_COEF_CONTINUE(val, DCT_VAL_CATEGORY6);
  }

#if CONFIG_CODE_NONZEROCOUNT
  if (!nzc_used)
#endif
    if (c < seg_eob)
      coef_counts[type][ref][get_coef_band(scan, txfm_size, c)]
                 [pt][DCT_EOB_TOKEN]++;
#if CONFIG_CODE_NONZEROCOUNT
  if (!nzc_used)
    xd->mode_info_context->mbmi.nzcs[block_idx] = nzc;
  else
    assert(nzc == nzc_expected);
#endif

  A0[aidx] = L0[lidx] = c > 0;
  if (txfm_size >= TX_8X8) {
    A0[aidx + 1] = L0[lidx + 1] = A0[aidx];
    if (txfm_size >= TX_16X16) {
      if (type == PLANE_TYPE_UV) {
        ENTROPY_CONTEXT *A1 = (ENTROPY_CONTEXT *) (xd->above_context + 1);
        ENTROPY_CONTEXT *L1 = (ENTROPY_CONTEXT *) (xd->left_context + 1);
        A1[aidx] = A1[aidx + 1] = L1[lidx] = L1[lidx + 1] = A0[aidx];
        if (txfm_size >= TX_32X32) {
          ENTROPY_CONTEXT *A2 = (ENTROPY_CONTEXT *) (xd->above_context + 2);
          ENTROPY_CONTEXT *L2 = (ENTROPY_CONTEXT *) (xd->left_context + 2);
          ENTROPY_CONTEXT *A3 = (ENTROPY_CONTEXT *) (xd->above_context + 3);
          ENTROPY_CONTEXT *L3 = (ENTROPY_CONTEXT *) (xd->left_context + 3);
          A2[aidx] = A2[aidx + 1] = A3[aidx] = A3[aidx + 1] = A0[aidx];
          L2[lidx] = L2[lidx + 1] = L3[lidx] = L3[lidx + 1] = A0[aidx];
        }
      } else {
        A0[aidx + 2] = A0[aidx + 3] = L0[lidx + 2] = L0[lidx + 3] = A0[aidx];
        if (txfm_size >= TX_32X32) {
          ENTROPY_CONTEXT *A1 = (ENTROPY_CONTEXT *) (xd->above_context + 1);
          ENTROPY_CONTEXT *L1 = (ENTROPY_CONTEXT *) (xd->left_context + 1);
          A1[aidx] = A1[aidx + 1] = A1[aidx + 2] = A1[aidx + 3] = A0[aidx];
          L1[lidx] = L1[lidx + 1] = L1[lidx + 2] = L1[lidx + 3] = A0[aidx];
        }
      }
    }
  }
  return c;
}

static int get_eob(MACROBLOCKD* const xd, int segment_id, int eob_max) {
  return vp9_get_segdata(xd, segment_id, SEG_LVL_SKIP) ? 0 : eob_max;
}

/* TODO(jkoleszar): Probably best to remove instances that require this,
 * as the data likely becomes per-plane and stored in the per-plane structures.
 * This is a stub to work with the existing code.
 */
static INLINE int block_idx_4x4(MACROBLOCKD* const xd, int block_size_b,
                                int plane, int i) {
  const int luma_blocks = 1 << block_size_b;
  assert(xd->plane[0].subsampling_x == 0);
  assert(xd->plane[0].subsampling_y == 0);
  assert(xd->plane[1].subsampling_x == 1);
  assert(xd->plane[1].subsampling_y == 1);
  assert(xd->plane[2].subsampling_x == 1);
  assert(xd->plane[2].subsampling_y == 1);
  return plane == 0 ? i :
         plane == 1 ? luma_blocks + i :
                      luma_blocks * 5 / 4 + i;
}

static INLINE int decode_block_plane(VP9D_COMP* const pbi,
                                     MACROBLOCKD* const xd,
                                     BOOL_DECODER* const bc,
                                     BLOCK_SIZE_LG2 block_size,
                                     int segment_id,
                                     int plane,
                                     int is_split) {
  // block and transform sizes, in number of 4x4 blocks log 2 ("*_b")
  // 4x4=0, 8x8=2, 16x16=4, 32x32=6, 64x64=8
  const TX_SIZE tx_size = xd->mode_info_context->mbmi.txfm_size;
  const BLOCK_SIZE_LG2 block_size_b = block_size;
  const BLOCK_SIZE_LG2 txfrm_size_b = tx_size * 2;

  // subsampled size of the block
  const int ss_sum = xd->plane[plane].subsampling_x +
                     xd->plane[plane].subsampling_y;
  const BLOCK_SIZE_LG2 ss_block_size = block_size_b - ss_sum;

  // size of the transform to use. scale the transform down if it's larger
  // than the size of the subsampled data, or forced externally by the mb mode.
  const int ss_max = MAX(xd->plane[plane].subsampling_x,
                         xd->plane[plane].subsampling_y);
  const BLOCK_SIZE_LG2 ss_txfrm_size = txfrm_size_b > ss_block_size || is_split
                                       ? txfrm_size_b - ss_max * 2
                                       : txfrm_size_b;
  const TX_SIZE ss_tx_size = ss_txfrm_size / 2;

  // TODO(jkoleszar): 1 may not be correct here with larger chroma planes.
  const int inc = is_split ? 1 : (1 << ss_txfrm_size);

  // find the maximum eob for this transform size, adjusted by segment
  const int seg_eob = get_eob(xd, segment_id, 16 << ss_txfrm_size);

  int i, eobtotal = 0;

  assert(txfrm_size_b <= block_size_b);
  assert(ss_txfrm_size <= ss_block_size);

  // step through the block by the size of the transform in use.
  for (i = 0; i < (1 << ss_block_size); i += inc) {
    const int block_idx = block_idx_4x4(xd, block_size_b, plane, i);

    const int c = decode_coefs(pbi, xd, bc, block_idx,
                               xd->plane[plane].plane_type, seg_eob,
                               BLOCK_OFFSET(xd->plane[plane].qcoeff, i, 16),
                               ss_tx_size);
    xd->plane[plane].eobs[i] = c;
    eobtotal += c;
  }
  return eobtotal;
}

static INLINE int decode_blocks_helper(VP9D_COMP* const pbi,
                                       MACROBLOCKD* const xd,
                                       BOOL_DECODER* const bc,
                                       int block_size,
                                       int is_split_chroma) {
  const int segment_id = xd->mode_info_context->mbmi.segment_id;
  int plane, eobtotal = 0;

  for (plane = 0; plane < MAX_MB_PLANE; plane++) {
    const int is_split = is_split_chroma &&
                         xd->plane[plane].plane_type == PLANE_TYPE_UV;
    eobtotal += decode_block_plane(pbi, xd, bc, block_size, segment_id,
                                   plane, is_split);
  }
  return eobtotal;
}

static INLINE int decode_blocks(VP9D_COMP* const pbi,
                                MACROBLOCKD* const xd,
                                BOOL_DECODER* const bc,
                                int block_size) {
  const MB_PREDICTION_MODE mode = xd->mode_info_context->mbmi.mode;
  const TX_SIZE tx_size = xd->mode_info_context->mbmi.txfm_size;
  return decode_blocks_helper(pbi, xd, bc, block_size,
      tx_size == TX_8X8 && (mode == I8X8_PRED || mode == SPLITMV));
}

int vp9_decode_sb64_tokens(VP9D_COMP* const pbi,
                           MACROBLOCKD* const xd,
                           BOOL_DECODER* const bc) {
  return decode_blocks(pbi, xd, bc, BLOCK_64X64_LG2);
}

int vp9_decode_sb_tokens(VP9D_COMP* const pbi,
                         MACROBLOCKD* const xd,
                         BOOL_DECODER* const bc) {
  return decode_blocks(pbi, xd, bc, BLOCK_32X32_LG2);
}

int vp9_decode_mb_tokens(VP9D_COMP* const pbi,
                         MACROBLOCKD* const xd,
                         BOOL_DECODER* const bc) {
  return decode_blocks(pbi, xd, bc, BLOCK_16X16_LG2);
}

#if CONFIG_NEWBINTRAMODES
static int decode_coefs_4x4(VP9D_COMP *dx, MACROBLOCKD *xd,
                            BOOL_DECODER* const bc,
                            PLANE_TYPE type, int i, int seg_eob) {
  const struct plane_block_idx pb_idx = plane_block_idx(16, i);
  const int c = decode_coefs(dx, xd, bc, i, type, seg_eob,
      BLOCK_OFFSET(xd->plane[pb_idx.plane].qcoeff, pb_idx.block, 16), TX_4X4);
  xd->plane[pb_idx.plane].eobs[pb_idx.block] = c;
  return c;
}

static int decode_mb_tokens_4x4_uv(VP9D_COMP* const dx,
                                   MACROBLOCKD* const xd,
                                   BOOL_DECODER* const bc,
                                   int seg_eob) {
  int i, eobtotal = 0;

  // chroma blocks
  for (i = 16; i < 24; i++)
    eobtotal += decode_coefs_4x4(dx, xd, bc, PLANE_TYPE_UV, i, seg_eob);

  return eobtotal;
}

int vp9_decode_mb_tokens_4x4_uv(VP9D_COMP* const dx,
                                MACROBLOCKD* const xd,
                                BOOL_DECODER* const bc) {
  const int segment_id = xd->mode_info_context->mbmi.segment_id;
  const int seg_eob = get_eob(xd, segment_id, 16);

  return decode_mb_tokens_4x4_uv(dx, xd, bc, seg_eob);
}

int vp9_decode_coefs_4x4(VP9D_COMP *dx, MACROBLOCKD *xd,
                         BOOL_DECODER* const bc,
                         PLANE_TYPE type, int i) {
  const int segment_id = xd->mode_info_context->mbmi.segment_id;
  const int seg_eob = get_eob(xd, segment_id, 16);
  return decode_coefs_4x4(dx, xd, bc, type, i, seg_eob);
}
#endif
