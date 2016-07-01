/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include "vpx_mem/vpx_mem.h"
#include "vpx_ports/mem.h"

#include "vp10/common/ans.h"
#include "vp10/common/blockd.h"
#include "vp10/common/common.h"
#include "vp10/common/entropy.h"
#include "vp10/common/idct.h"

#include "vp10/decoder/detokenize.h"

#define EOB_CONTEXT_NODE            0
#define ZERO_CONTEXT_NODE           1
#define ONE_CONTEXT_NODE            2
#define LOW_VAL_CONTEXT_NODE        0
#define TWO_CONTEXT_NODE            1
#define THREE_CONTEXT_NODE          2
#define HIGH_LOW_CONTEXT_NODE       3
#define CAT_ONE_CONTEXT_NODE        4
#define CAT_THREEFOUR_CONTEXT_NODE  5
#define CAT_THREE_CONTEXT_NODE      6
#define CAT_FIVE_CONTEXT_NODE       7

#define INCREMENT_COUNT(token)                              \
  do {                                                      \
     if (counts)                                            \
       ++coef_counts[band][ctx][token];                     \
  } while (0)

#if !CONFIG_ANS
static INLINE int read_coeff(const vpx_prob *probs, int n, vp10_reader *r) {
  int i, val = 0;
  for (i = 0; i < n; ++i)
    val = (val << 1) | vp10_read(r, probs[i]);
  return val;
}

static int decode_coefs(const MACROBLOCKD *xd,
                        PLANE_TYPE type,
                        tran_low_t *dqcoeff, TX_SIZE tx_size, TX_TYPE tx_type,
                        const int16_t *dq,
#if CONFIG_NEW_QUANT
                        dequant_val_type_nuq *dq_val,
#endif  // CONFIG_NEW_QUANT
                        int ctx, const int16_t *scan, const int16_t *nb,
                        vp10_reader *r) {
  FRAME_COUNTS *counts = xd->counts;
  const int max_eob = get_tx2d_size(tx_size);
  const FRAME_CONTEXT *const fc = xd->fc;
  const int ref = is_inter_block(&xd->mi[0]->mbmi);
  int band, c = 0;
  const int tx_size_ctx = txsize_sqr_map[tx_size];
  const vpx_prob (*coef_probs)[COEFF_CONTEXTS][UNCONSTRAINED_NODES] =
      fc->coef_probs[tx_size_ctx][type][ref];
  const vpx_prob *prob;
  unsigned int (*coef_counts)[COEFF_CONTEXTS][UNCONSTRAINED_NODES + 1];
  unsigned int (*eob_branch_count)[COEFF_CONTEXTS];
  uint8_t token_cache[MAX_TX_SQUARE];
  const uint8_t *band_translate = get_band_translate(tx_size);
  int dq_shift;
  int v, token;
  int16_t dqv = dq[0];
#if CONFIG_NEW_QUANT
  const tran_low_t *dqv_val = &dq_val[0][0];
#endif  // CONFIG_NEW_QUANT
  const uint8_t *cat1_prob;
  const uint8_t *cat2_prob;
  const uint8_t *cat3_prob;
  const uint8_t *cat4_prob;
  const uint8_t *cat5_prob;
  const uint8_t *cat6_prob;

  if (counts) {
    coef_counts = counts->coef[tx_size_ctx][type][ref];
    eob_branch_count = counts->eob_branch[tx_size_ctx][type][ref];
  }

#if CONFIG_VP9_HIGHBITDEPTH
  if (xd->bd > VPX_BITS_8) {
    if (xd->bd == VPX_BITS_10) {
      cat1_prob = vp10_cat1_prob_high10;
      cat2_prob = vp10_cat2_prob_high10;
      cat3_prob = vp10_cat3_prob_high10;
      cat4_prob = vp10_cat4_prob_high10;
      cat5_prob = vp10_cat5_prob_high10;
      cat6_prob = vp10_cat6_prob_high10;
    } else {
      cat1_prob = vp10_cat1_prob_high12;
      cat2_prob = vp10_cat2_prob_high12;
      cat3_prob = vp10_cat3_prob_high12;
      cat4_prob = vp10_cat4_prob_high12;
      cat5_prob = vp10_cat5_prob_high12;
      cat6_prob = vp10_cat6_prob_high12;
    }
  } else {
    cat1_prob = vp10_cat1_prob;
    cat2_prob = vp10_cat2_prob;
    cat3_prob = vp10_cat3_prob;
    cat4_prob = vp10_cat4_prob;
    cat5_prob = vp10_cat5_prob;
    cat6_prob = vp10_cat6_prob;
  }
#else
  cat1_prob = vp10_cat1_prob;
  cat2_prob = vp10_cat2_prob;
  cat3_prob = vp10_cat3_prob;
  cat4_prob = vp10_cat4_prob;
  cat5_prob = vp10_cat5_prob;
  cat6_prob = vp10_cat6_prob;
#endif

  dq_shift = get_tx_scale(xd, tx_type, tx_size);

  while (c < max_eob) {
    int val = -1;
    band = *band_translate++;
    prob = coef_probs[band][ctx];
    if (counts)
      ++eob_branch_count[band][ctx];
    if (!vp10_read(r, prob[EOB_CONTEXT_NODE])) {
      INCREMENT_COUNT(EOB_MODEL_TOKEN);
      break;
    }

#if CONFIG_NEW_QUANT
    dqv_val = &dq_val[band][0];
#endif  // CONFIG_NEW_QUANT

    while (!vp10_read(r, prob[ZERO_CONTEXT_NODE])) {
      INCREMENT_COUNT(ZERO_TOKEN);
      dqv = dq[1];
      token_cache[scan[c]] = 0;
      ++c;
      if (c >= max_eob)
        return c;  // zero tokens at the end (no eob token)
      ctx = get_coef_context(nb, token_cache, c);
      band = *band_translate++;
      prob = coef_probs[band][ctx];
#if CONFIG_NEW_QUANT
      dqv_val = &dq_val[band][0];
#endif  // CONFIG_NEW_QUANT
    }

    if (!vp10_read(r, prob[ONE_CONTEXT_NODE])) {
      INCREMENT_COUNT(ONE_TOKEN);
      token = ONE_TOKEN;
      val = 1;
    } else {
      INCREMENT_COUNT(TWO_TOKEN);
      token = vp10_read_tree(r, vp10_coef_con_tree,
                            vp10_pareto8_full[prob[PIVOT_NODE] - 1]);
      switch (token) {
        case TWO_TOKEN:
        case THREE_TOKEN:
        case FOUR_TOKEN:
          val = token;
          break;
        case CATEGORY1_TOKEN:
          val = CAT1_MIN_VAL + read_coeff(cat1_prob, 1, r);
          break;
        case CATEGORY2_TOKEN:
          val = CAT2_MIN_VAL + read_coeff(cat2_prob, 2, r);
          break;
        case CATEGORY3_TOKEN:
          val = CAT3_MIN_VAL + read_coeff(cat3_prob, 3, r);
          break;
        case CATEGORY4_TOKEN:
          val = CAT4_MIN_VAL + read_coeff(cat4_prob, 4, r);
          break;
        case CATEGORY5_TOKEN:
          val = CAT5_MIN_VAL + read_coeff(cat5_prob, 5, r);
          break;
        case CATEGORY6_TOKEN: {
          const int skip_bits = TX_SIZES - 1 - tx_size;
          const uint8_t *cat6p = cat6_prob + skip_bits;
#if CONFIG_VP9_HIGHBITDEPTH
          switch (xd->bd) {
            case VPX_BITS_8:
              val = CAT6_MIN_VAL + read_coeff(cat6p, 14 - skip_bits, r);
              break;
            case VPX_BITS_10:
              val = CAT6_MIN_VAL + read_coeff(cat6p, 16 - skip_bits, r);
              break;
            case VPX_BITS_12:
              val = CAT6_MIN_VAL + read_coeff(cat6p, 18 - skip_bits, r);
              break;
            default:
              assert(0);
              return -1;
          }
#else
          val = CAT6_MIN_VAL + read_coeff(cat6p, 14 - skip_bits, r);
#endif
          break;
        }
      }
    }
#if CONFIG_NEW_QUANT
    v = vp10_dequant_abscoeff_nuq(val, dqv, dqv_val);
    v = dq_shift ? ROUND_POWER_OF_TWO(v, dq_shift) : v;
#else
    v = (val * dqv) >> dq_shift;
#endif  // CONFIG_NEW_QUANT

#if CONFIG_COEFFICIENT_RANGE_CHECKING
#if CONFIG_VP9_HIGHBITDEPTH
    dqcoeff[scan[c]] = highbd_check_range((vp10_read_bit(r) ? -v : v),
                                          xd->bd);
#else
    dqcoeff[scan[c]] = check_range(vp10_read_bit(r) ? -v : v);
#endif  // CONFIG_VP9_HIGHBITDEPTH
#else
    dqcoeff[scan[c]] = vp10_read_bit(r) ? -v : v;
#endif  // CONFIG_COEFFICIENT_RANGE_CHECKING
    token_cache[scan[c]] = vp10_pt_energy_class[token];
    ++c;
    ctx = get_coef_context(nb, token_cache, c);
    dqv = dq[1];
  }

  return c;
}
#else  // !CONFIG_ANS
static INLINE int read_coeff(const vpx_prob *const probs, int n,
                             struct AnsDecoder *const ans) {
  int i, val = 0;
  for (i = 0; i < n; ++i)
    val = (val << 1) | uabs_read(ans, probs[i]);
  return val;
}

static int decode_coefs_ans(const MACROBLOCKD *const xd,
                            PLANE_TYPE type,
                            tran_low_t *dqcoeff, TX_SIZE tx_size,
                            TX_TYPE tx_type,
                            const int16_t *dq,
#if CONFIG_NEW_QUANT
                            dequant_val_type_nuq *dq_val,
#endif  // CONFIG_NEW_QUANT
                            int ctx, const int16_t *scan, const int16_t *nb,
                            struct AnsDecoder *const ans) {
  FRAME_COUNTS *counts = xd->counts;
  const int max_eob = get_tx2d_size(tx_size);
  const FRAME_CONTEXT *const fc = xd->fc;
  const int ref = is_inter_block(&xd->mi[0]->mbmi);
  int band, c = 0;
  int skip_eob = 0;
  const int tx_size_ctx = txsize_sqr_map[tx_size];
  const vpx_prob (*coef_probs)[COEFF_CONTEXTS][UNCONSTRAINED_NODES] =
      fc->coef_probs[tx_size_ctx][type][ref];
  const rans_dec_lut(*coef_cdfs)[COEFF_CONTEXTS] =
      fc->coef_cdfs[tx_size_ctx][type][ref];
  const vpx_prob *prob;
  const rans_dec_lut *cdf;
  unsigned int (*coef_counts)[COEFF_CONTEXTS][UNCONSTRAINED_NODES + 1];
  unsigned int (*eob_branch_count)[COEFF_CONTEXTS];
  uint8_t token_cache[MAX_TX_SQUARE];
  const uint8_t *band_translate = get_band_translate(tx_size);
  int dq_shift;
  int v, token;
  int16_t dqv = dq[0];
#if CONFIG_NEW_QUANT
  const tran_low_t *dqv_val = &dq_val[0][0];
#endif  // CONFIG_NEW_QUANT
  const uint8_t *cat1_prob;
  const uint8_t *cat2_prob;
  const uint8_t *cat3_prob;
  const uint8_t *cat4_prob;
  const uint8_t *cat5_prob;
  const uint8_t *cat6_prob;

  dq_shift = get_tx_scale(xd, tx_type, tx_size);

  if (counts) {
    coef_counts = counts->coef[tx_size_ctx][type][ref];
    eob_branch_count = counts->eob_branch[tx_size_ctx][type][ref];
  }

#if CONFIG_VP9_HIGHBITDEPTH
  if (xd->bd > VPX_BITS_8) {
    if (xd->bd == VPX_BITS_10) {
      cat1_prob = vp10_cat1_prob_high10;
      cat2_prob = vp10_cat2_prob_high10;
      cat3_prob = vp10_cat3_prob_high10;
      cat4_prob = vp10_cat4_prob_high10;
      cat5_prob = vp10_cat5_prob_high10;
      cat6_prob = vp10_cat6_prob_high10;
    } else {
      cat1_prob = vp10_cat1_prob_high12;
      cat2_prob = vp10_cat2_prob_high12;
      cat3_prob = vp10_cat3_prob_high12;
      cat4_prob = vp10_cat4_prob_high12;
      cat5_prob = vp10_cat5_prob_high12;
      cat6_prob = vp10_cat6_prob_high12;
    }
  } else {
    cat1_prob = vp10_cat1_prob;
    cat2_prob = vp10_cat2_prob;
    cat3_prob = vp10_cat3_prob;
    cat4_prob = vp10_cat4_prob;
    cat5_prob = vp10_cat5_prob;
    cat6_prob = vp10_cat6_prob;
  }
#else
  cat1_prob = vp10_cat1_prob;
  cat2_prob = vp10_cat2_prob;
  cat3_prob = vp10_cat3_prob;
  cat4_prob = vp10_cat4_prob;
  cat5_prob = vp10_cat5_prob;
  cat6_prob = vp10_cat6_prob;
#endif

  while (c < max_eob) {
    int val = -1;
    band = *band_translate++;
    prob = coef_probs[band][ctx];
    if (!skip_eob) {
      if (counts)
        ++eob_branch_count[band][ctx];
      if (!uabs_read(ans, prob[EOB_CONTEXT_NODE])) {
        INCREMENT_COUNT(EOB_MODEL_TOKEN);
        break;
      }
    }

#if CONFIG_NEW_QUANT
    dqv_val = &dq_val[band][0];
#endif  // CONFIG_NEW_QUANT

    cdf = &coef_cdfs[band][ctx];
    token = ZERO_TOKEN + rans_read(ans, *cdf);
    if (token == ZERO_TOKEN) {
      INCREMENT_COUNT(ZERO_TOKEN);
      token_cache[scan[c]] = 0;
      skip_eob = 1;
    } else {
      INCREMENT_COUNT(ONE_TOKEN + (token > ONE_TOKEN));
      switch (token) {
        case ONE_TOKEN:
        case TWO_TOKEN:
        case THREE_TOKEN:
        case FOUR_TOKEN:
          val = token;
          break;
        case CATEGORY1_TOKEN:
          val = CAT1_MIN_VAL + read_coeff(cat1_prob, 1, ans);
          break;
        case CATEGORY2_TOKEN:
          val = CAT2_MIN_VAL + read_coeff(cat2_prob, 2, ans);
          break;
        case CATEGORY3_TOKEN:
          val = CAT3_MIN_VAL + read_coeff(cat3_prob, 3, ans);
          break;
        case CATEGORY4_TOKEN:
          val = CAT4_MIN_VAL + read_coeff(cat4_prob, 4, ans);
          break;
        case CATEGORY5_TOKEN:
          val = CAT5_MIN_VAL + read_coeff(cat5_prob, 5, ans);
          break;
        case CATEGORY6_TOKEN: {
          const int skip_bits = TX_SIZES - 1 - tx_size;
          const uint8_t *cat6p = cat6_prob + skip_bits;
#if CONFIG_VP9_HIGHBITDEPTH
          switch (xd->bd) {
            case VPX_BITS_8:
              val = CAT6_MIN_VAL + read_coeff(cat6p, 14 - skip_bits, ans);
              break;
            case VPX_BITS_10:
              val = CAT6_MIN_VAL + read_coeff(cat6p, 16 - skip_bits, ans);
              break;
            case VPX_BITS_12:
              val = CAT6_MIN_VAL + read_coeff(cat6p, 18 - skip_bits, ans);
              break;
            default:
              assert(0);
              return -1;
          }
#else
          val = CAT6_MIN_VAL + read_coeff(cat6p, 14 - skip_bits, ans);
#endif
        } break;
      }
#if CONFIG_NEW_QUANT
    v = vp10_dequant_abscoeff_nuq(val, dqv, dqv_val);
    v = dq_shift ? ROUND_POWER_OF_TWO(v, dq_shift) : v;
#else
    v = (val * dqv) >> dq_shift;
#endif  // CONFIG_NEW_QUANT

#if CONFIG_COEFFICIENT_RANGE_CHECKING
#if CONFIG_VP9_HIGHBITDEPTH
      dqcoeff[scan[c]] =
          highbd_check_range((uabs_read_bit(ans) ? -v : v), xd->bd);
#else
      dqcoeff[scan[c]] = check_range(uabs_read_bit(ans) ? -v : v);
#endif  // CONFIG_VP9_HIGHBITDEPTH
#else
      dqcoeff[scan[c]] = uabs_read_bit(ans) ? -v : v;
#endif  // CONFIG_COEFFICIENT_RANGE_CHECKING
      token_cache[scan[c]] = vp10_pt_energy_class[token];
      skip_eob = 0;
    }
    ++c;
    ctx = get_coef_context(nb, token_cache, c);
    dqv = dq[1];
  }

  return c;
}
#endif  // !CONFIG_ANS

// TODO(slavarnway): Decode version of vp10_set_context.  Modify vp10_set_context
// after testing is complete, then delete this version.
static
void dec_set_contexts(const MACROBLOCKD *xd, struct macroblockd_plane *pd,
                      TX_SIZE tx_size, int has_eob,
                      int aoff, int loff) {
  ENTROPY_CONTEXT *const a = pd->above_context + aoff;
  ENTROPY_CONTEXT *const l = pd->left_context + loff;
  const int tx_w_in_blocks = num_4x4_blocks_wide_txsize_lookup[tx_size];
  const int tx_h_in_blocks = num_4x4_blocks_high_txsize_lookup[tx_size];

  // above
  if (has_eob && xd->mb_to_right_edge < 0) {
    int i;
    const int blocks_wide = pd->n4_w +
                            (xd->mb_to_right_edge >> (5 + pd->subsampling_x));
    int above_contexts = tx_w_in_blocks;
    if (above_contexts + aoff > blocks_wide)
      above_contexts = blocks_wide - aoff;

    for (i = 0; i < above_contexts; ++i)
      a[i] = has_eob;
    for (i = above_contexts; i < tx_w_in_blocks; ++i)
      a[i] = 0;
  } else {
    memset(a, has_eob, sizeof(ENTROPY_CONTEXT) * tx_w_in_blocks);
  }

  // left
  if (has_eob && xd->mb_to_bottom_edge < 0) {
    int i;
    const int blocks_high = pd->n4_h +
                            (xd->mb_to_bottom_edge >> (5 + pd->subsampling_y));
    int left_contexts = tx_h_in_blocks;
    if (left_contexts + loff > blocks_high)
      left_contexts = blocks_high - loff;

    for (i = 0; i < left_contexts; ++i)
      l[i] = has_eob;
    for (i = left_contexts; i < tx_h_in_blocks; ++i)
      l[i] = 0;
  } else {
    memset(l, has_eob, sizeof(ENTROPY_CONTEXT) * tx_h_in_blocks);
  }
}

void vp10_decode_palette_tokens(MACROBLOCKD *const xd, int plane,
                                vp10_reader *r) {
  MODE_INFO *const mi = xd->mi[0];
  MB_MODE_INFO *const mbmi = &mi->mbmi;
  const BLOCK_SIZE bsize = mbmi->sb_type;
  const int rows = (4 * num_4x4_blocks_high_lookup[bsize]) >>
      (xd->plane[plane != 0].subsampling_y);
  const int cols = (4 * num_4x4_blocks_wide_lookup[bsize]) >>
      (xd->plane[plane != 0].subsampling_x);
  int color_idx, color_ctx, color_order[PALETTE_MAX_SIZE];
  int n = mbmi->palette_mode_info.palette_size[plane != 0];
  int i, j;
  uint8_t *color_map = xd->plane[plane != 0].color_index_map;
  const vpx_prob (* const prob)[PALETTE_COLOR_CONTEXTS][PALETTE_COLORS - 1] =
      plane ? vp10_default_palette_uv_color_prob :
          vp10_default_palette_y_color_prob;

  for (i = 0; i < rows; ++i) {
    for (j = (i == 0 ? 1 : 0); j < cols; ++j) {
      color_ctx = vp10_get_palette_color_context(color_map, cols, i, j, n,
                                                 color_order);
      color_idx = vp10_read_tree(r, vp10_palette_color_tree[n - 2],
                                prob[n - 2][color_ctx]);
      assert(color_idx >= 0 && color_idx < n);
      color_map[i * cols + j] = color_order[color_idx];
    }
  }
}

int vp10_decode_block_tokens(MACROBLOCKD *const xd,
                             int plane, const scan_order *sc,
                             int x, int y,
                             TX_SIZE tx_size,
                             TX_TYPE tx_type,
#if CONFIG_ANS
                             struct AnsDecoder *const r,
#else
                             vp10_reader *r,
#endif  // CONFIG_ANS
                             int seg_id) {
  struct macroblockd_plane *const pd = &xd->plane[plane];
  const int16_t *const dequant = pd->seg_dequant[seg_id];
  const int ctx = get_entropy_context(tx_size, pd->above_context + x,
                                               pd->left_context + y);
#if CONFIG_NEW_QUANT
  int dq = get_dq_profile_from_ctx(ctx);
#endif  //  CONFIG_NEW_QUANT

#if !CONFIG_ANS
  const int eob = decode_coefs(xd, pd->plane_type,
                               pd->dqcoeff, tx_size, tx_type,
                               dequant,
#if CONFIG_NEW_QUANT
                               pd->seg_dequant_nuq[seg_id][dq],
#endif  // CONFIG_NEW_QUANT
                               ctx, sc->scan, sc->neighbors, r);
#else
  const int eob = decode_coefs_ans(xd, pd->plane_type,
                                   pd->dqcoeff, tx_size, tx_type,
                                   dequant,
#if CONFIG_NEW_QUANT
                                   pd->seg_dequant_nuq[seg_id][dq],
#endif  // CONFIG_NEW_QUANT
                                   ctx, sc->scan, sc->neighbors, r);
#endif  // !CONFIG_ANS
  dec_set_contexts(xd, pd, tx_size, eob > 0, x, y);
  /*
  vp10_set_contexts(xd, pd,
                    get_plane_block_size(xd->mi[0]->mbmi.sb_type, pd),
                    tx_size, eob > 0, x, y);
                    */
  return eob;
}
