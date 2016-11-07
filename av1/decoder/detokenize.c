/*
 * Copyright (c) 2016, Alliance for Open Media. All rights reserved
 *
 * This source code is subject to the terms of the BSD 2 Clause License and
 * the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
 * was not distributed with this source code in the LICENSE file, you can
 * obtain it at www.aomedia.org/license/software. If the Alliance for Open
 * Media Patent License 1.0 was not distributed with this source code in the
 * PATENTS file, you can obtain it at www.aomedia.org/license/patent.
 */

#if !CONFIG_PVQ
#include "aom_mem/aom_mem.h"
#include "aom_ports/mem.h"
#endif
#if !CONFIG_PVQ
#if CONFIG_ANS
#include "aom_dsp/ans.h"
#endif  // CONFIG_ANS
#include "av1/common/blockd.h"
#include "av1/common/common.h"
#include "av1/common/entropy.h"
#include "av1/common/idct.h"

#include "av1/decoder/detokenize.h"

#define ACCT_STR __func__

#define EOB_CONTEXT_NODE 0
#define ZERO_CONTEXT_NODE 1
#define ONE_CONTEXT_NODE 2
#define LOW_VAL_CONTEXT_NODE 0
#define TWO_CONTEXT_NODE 1
#define THREE_CONTEXT_NODE 2
#define HIGH_LOW_CONTEXT_NODE 3
#define CAT_ONE_CONTEXT_NODE 4
#define CAT_THREEFOUR_CONTEXT_NODE 5
#define CAT_THREE_CONTEXT_NODE 6
#define CAT_FIVE_CONTEXT_NODE 7

#define INCREMENT_COUNT(token)                   \
  do {                                           \
    if (counts) ++coef_counts[band][ctx][token]; \
  } while (0)

static INLINE int read_coeff(const aom_prob *probs, int n, aom_reader *r) {
  int i, val = 0;
  for (i = 0; i < n; ++i) val = (val << 1) | aom_read(r, probs[i], ACCT_STR);
  return val;
}

#if CONFIG_AOM_QM
static int decode_coefs(MACROBLOCKD *xd, PLANE_TYPE type, tran_low_t *dqcoeff,
                        TX_SIZE tx_size, TX_TYPE tx_type, const int16_t *dq,
                        int ctx, const int16_t *scan, const int16_t *nb,
                        int16_t *max_scan_line, aom_reader *r,
                        const qm_val_t *iqm[2][TX_SIZES])
#else
static int decode_coefs(MACROBLOCKD *xd, PLANE_TYPE type, tran_low_t *dqcoeff,
                        TX_SIZE tx_size, TX_TYPE tx_type, const int16_t *dq,
#if CONFIG_NEW_QUANT
                        dequant_val_type_nuq *dq_val,
#endif  // CONFIG_NEW_QUANT
                        int ctx, const int16_t *scan, const int16_t *nb,
                        int16_t *max_scan_line, aom_reader *r)
#endif
{
  FRAME_COUNTS *counts = xd->counts;
  FRAME_CONTEXT *const fc = xd->fc;
  const int max_eob = tx_size_2d[tx_size];
  const int ref = is_inter_block(&xd->mi[0]->mbmi);
#if CONFIG_AOM_QM
  const qm_val_t *iqmatrix = iqm[!ref][tx_size];
#endif
  int band, c = 0;
  const int tx_size_ctx = txsize_sqr_map[tx_size];
  aom_prob(*coef_probs)[COEFF_CONTEXTS][UNCONSTRAINED_NODES] =
      fc->coef_probs[tx_size_ctx][type][ref];
  const aom_prob *prob;
#if CONFIG_EC_MULTISYMBOL
  aom_cdf_prob(*coef_cdfs)[COEFF_CONTEXTS][ENTROPY_TOKENS] =
      fc->coef_cdfs[tx_size_ctx][type][ref];
  aom_cdf_prob(*cdf)[ENTROPY_TOKENS];
#endif  // CONFIG_EC_MULTISYMBOL
  unsigned int(*coef_counts)[COEFF_CONTEXTS][UNCONSTRAINED_NODES + 1];
  unsigned int(*eob_branch_count)[COEFF_CONTEXTS];
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

#if CONFIG_AOM_HIGHBITDEPTH
  if (xd->bd > AOM_BITS_8) {
    if (xd->bd == AOM_BITS_10) {
      cat1_prob = av1_cat1_prob_high10;
      cat2_prob = av1_cat2_prob_high10;
      cat3_prob = av1_cat3_prob_high10;
      cat4_prob = av1_cat4_prob_high10;
      cat5_prob = av1_cat5_prob_high10;
      cat6_prob = av1_cat6_prob_high10;
    } else {
      cat1_prob = av1_cat1_prob_high12;
      cat2_prob = av1_cat2_prob_high12;
      cat3_prob = av1_cat3_prob_high12;
      cat4_prob = av1_cat4_prob_high12;
      cat5_prob = av1_cat5_prob_high12;
      cat6_prob = av1_cat6_prob_high12;
    }
  } else {
    cat1_prob = av1_cat1_prob;
    cat2_prob = av1_cat2_prob;
    cat3_prob = av1_cat3_prob;
    cat4_prob = av1_cat4_prob;
    cat5_prob = av1_cat5_prob;
    cat6_prob = av1_cat6_prob;
  }
#else
  cat1_prob = av1_cat1_prob;
  cat2_prob = av1_cat2_prob;
  cat3_prob = av1_cat3_prob;
  cat4_prob = av1_cat4_prob;
  cat5_prob = av1_cat5_prob;
  cat6_prob = av1_cat6_prob;
#endif

  dq_shift = get_tx_scale(xd, tx_type, tx_size);

  while (c < max_eob) {
    int val = -1;
    band = *band_translate++;
    prob = coef_probs[band][ctx];
    if (counts) ++eob_branch_count[band][ctx];
    if (!aom_read(r, prob[EOB_CONTEXT_NODE], ACCT_STR)) {
      INCREMENT_COUNT(EOB_MODEL_TOKEN);
      break;
    }

#if CONFIG_NEW_QUANT
    dqv_val = &dq_val[band][0];
#endif  // CONFIG_NEW_QUANT

    while (!aom_read(r, prob[ZERO_CONTEXT_NODE], ACCT_STR)) {
      INCREMENT_COUNT(ZERO_TOKEN);
      dqv = dq[1];
      token_cache[scan[c]] = 0;
      ++c;
      if (c >= max_eob) return c;  // zero tokens at the end (no eob token)
      ctx = get_coef_context(nb, token_cache, c);
      band = *band_translate++;
      prob = coef_probs[band][ctx];
#if CONFIG_NEW_QUANT
      dqv_val = &dq_val[band][0];
#endif  // CONFIG_NEW_QUANT
    }

    *max_scan_line = AOMMAX(*max_scan_line, scan[c]);

#if CONFIG_EC_MULTISYMBOL
    cdf = &coef_cdfs[band][ctx];
    token = ONE_TOKEN +
            aom_read_symbol(r, *cdf, CATEGORY6_TOKEN - ONE_TOKEN + 1, ACCT_STR);
    INCREMENT_COUNT(ONE_TOKEN + (token > ONE_TOKEN));
    switch (token) {
      case ONE_TOKEN:
      case TWO_TOKEN:
      case THREE_TOKEN:
      case FOUR_TOKEN: val = token; break;
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
        const int skip_bits = TX_SIZES - 1 - txsize_sqr_up_map[tx_size];
        const uint8_t *cat6p = cat6_prob + skip_bits;
#if CONFIG_AOM_HIGHBITDEPTH
        switch (xd->bd) {
          case AOM_BITS_8:
            val = CAT6_MIN_VAL + read_coeff(cat6p, 14 - skip_bits, r);
            break;
          case AOM_BITS_10:
            val = CAT6_MIN_VAL + read_coeff(cat6p, 16 - skip_bits, r);
            break;
          case AOM_BITS_12:
            val = CAT6_MIN_VAL + read_coeff(cat6p, 18 - skip_bits, r);
            break;
          default: assert(0); return -1;
        }
#else
        val = CAT6_MIN_VAL + read_coeff(cat6p, 14 - skip_bits, r);
#endif
      } break;
    }
#else  // CONFIG_EC_MULTISYMBOL
    if (!aom_read(r, prob[ONE_CONTEXT_NODE], ACCT_STR)) {
      INCREMENT_COUNT(ONE_TOKEN);
      token = ONE_TOKEN;
      val = 1;
    } else {
      INCREMENT_COUNT(TWO_TOKEN);
      token = aom_read_tree(r, av1_coef_con_tree,
                            av1_pareto8_full[prob[PIVOT_NODE] - 1], ACCT_STR);
      switch (token) {
        case TWO_TOKEN:
        case THREE_TOKEN:
        case FOUR_TOKEN: val = token; break;
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
          const int skip_bits = TX_SIZES - 1 - txsize_sqr_up_map[tx_size];
          const uint8_t *cat6p = cat6_prob + skip_bits;
#if CONFIG_AOM_HIGHBITDEPTH
          switch (xd->bd) {
            case AOM_BITS_8:
              val = CAT6_MIN_VAL + read_coeff(cat6p, 14 - skip_bits, r);
              break;
            case AOM_BITS_10:
              val = CAT6_MIN_VAL + read_coeff(cat6p, 16 - skip_bits, r);
              break;
            case AOM_BITS_12:
              val = CAT6_MIN_VAL + read_coeff(cat6p, 18 - skip_bits, r);
              break;
            default: assert(0); return -1;
          }
#else
          val = CAT6_MIN_VAL + read_coeff(cat6p, 14 - skip_bits, r);
#endif
          break;
        }
      }
    }
#endif  // CONFIG_EC_MULTISYMBOL
#if CONFIG_NEW_QUANT
    v = av1_dequant_abscoeff_nuq(val, dqv, dqv_val);
    v = dq_shift ? ROUND_POWER_OF_TWO(v, dq_shift) : v;
#else
#if CONFIG_AOM_QM
    dqv = ((iqmatrix[scan[c]] * (int)dqv) + (1 << (AOM_QM_BITS - 1))) >>
          AOM_QM_BITS;
#endif
    v = (val * dqv) >> dq_shift;
#endif  // CONFIG_NEW_QUANT

#if CONFIG_COEFFICIENT_RANGE_CHECKING
#if CONFIG_AOM_HIGHBITDEPTH
    dqcoeff[scan[c]] =
        highbd_check_range((aom_read_bit(r, ACCT_STR) ? -v : v), xd->bd);
#else
    dqcoeff[scan[c]] = check_range(aom_read_bit(r, ACCT_STR) ? -v : v);
#endif  // CONFIG_AOM_HIGHBITDEPTH
#else
    dqcoeff[scan[c]] = aom_read_bit(r, ACCT_STR) ? -v : v;
#endif  // CONFIG_COEFFICIENT_RANGE_CHECKING
    token_cache[scan[c]] = av1_pt_energy_class[token];
    ++c;
    ctx = get_coef_context(nb, token_cache, c);
    dqv = dq[1];
  }

  return c;
}

#if CONFIG_PALETTE
void av1_decode_palette_tokens(MACROBLOCKD *const xd, int plane,
                               aom_reader *r) {
  const MODE_INFO *const mi = xd->mi[0];
  const MB_MODE_INFO *const mbmi = &mi->mbmi;
  const BLOCK_SIZE bsize = mbmi->sb_type;
  const int rows = (4 * num_4x4_blocks_high_lookup[bsize]) >>
                   (xd->plane[plane != 0].subsampling_y);
  const int cols = (4 * num_4x4_blocks_wide_lookup[bsize]) >>
                   (xd->plane[plane != 0].subsampling_x);
  uint8_t color_order[PALETTE_MAX_SIZE];
  const int n = mbmi->palette_mode_info.palette_size[plane != 0];
  int i, j;
  uint8_t *color_map = xd->plane[plane != 0].color_index_map;
  const aom_prob(*const prob)[PALETTE_COLOR_CONTEXTS][PALETTE_COLORS - 1] =
      plane ? av1_default_palette_uv_color_prob
            : av1_default_palette_y_color_prob;

  for (i = 0; i < rows; ++i) {
    for (j = (i == 0 ? 1 : 0); j < cols; ++j) {
      const int color_ctx = av1_get_palette_color_context(color_map, cols, i, j,
                                                          n, color_order, NULL);
      const int color_idx = aom_read_tree(r, av1_palette_color_tree[n - 2],
                                          prob[n - 2][color_ctx], ACCT_STR);
      assert(color_idx >= 0 && color_idx < n);
      color_map[i * cols + j] = color_order[color_idx];
    }
  }
}
#endif  // CONFIG_PALETTE

int av1_decode_block_tokens(MACROBLOCKD *const xd, int plane,
                            const SCAN_ORDER *sc, int x, int y, TX_SIZE tx_size,
                            TX_TYPE tx_type, int16_t *max_scan_line,
                            aom_reader *r, int seg_id) {
  struct macroblockd_plane *const pd = &xd->plane[plane];
  const int16_t *const dequant = pd->seg_dequant[seg_id];
  const int ctx =
      get_entropy_context(tx_size, pd->above_context + x, pd->left_context + y);
#if CONFIG_NEW_QUANT
  const int ref = is_inter_block(&xd->mi[0]->mbmi);
  int dq =
      get_dq_profile_from_ctx(xd->qindex[seg_id], ctx, ref, pd->plane_type);
#endif  //  CONFIG_NEW_QUANT

#if CONFIG_AOM_QM
  const int eob = decode_coefs(xd, pd->plane_type, pd->dqcoeff, tx_size,
                               tx_type, dequant, ctx, sc->scan, sc->neighbors,
                               max_scan_line, r, pd->seg_iqmatrix[seg_id]);
#else
  const int eob =
      decode_coefs(xd, pd->plane_type, pd->dqcoeff, tx_size, tx_type, dequant,
#if CONFIG_NEW_QUANT
                   pd->seg_dequant_nuq[seg_id][dq],
#endif  // CONFIG_NEW_QUANT
                   ctx, sc->scan, sc->neighbors, max_scan_line, r);
#endif  // CONFIG_AOM_QM
  av1_set_contexts(xd, pd, tx_size, eob > 0, x, y);
  return eob;
}
#endif
