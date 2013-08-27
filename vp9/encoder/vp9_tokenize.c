/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */


#include <math.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include "vp9/encoder/vp9_onyx_int.h"
#include "vp9/encoder/vp9_tokenize.h"
#include "vpx_mem/vpx_mem.h"

#include "vp9/common/vp9_pred_common.h"
#include "vp9/common/vp9_seg_common.h"
#include "vp9/common/vp9_entropy.h"

/* Global event counters used for accumulating statistics across several
   compressions, then generating vp9_context.c = initial stats. */

#ifdef ENTROPY_STATS
vp9_coeff_accum context_counters[TX_SIZES][BLOCK_TYPES];
extern vp9_coeff_stats tree_update_hist[TX_SIZES][BLOCK_TYPES];
#endif  /* ENTROPY_STATS */

DECLARE_ALIGNED(16, extern const uint8_t,
                vp9_pt_energy_class[MAX_ENTROPY_TOKENS]);

static TOKENVALUE dct_value_tokens[DCT_MAX_VALUE * 2];
const TOKENVALUE *vp9_dct_value_tokens_ptr;
static int dct_value_cost[DCT_MAX_VALUE * 2];
const int *vp9_dct_value_cost_ptr;

static void fill_value_tokens() {

  TOKENVALUE *const t = dct_value_tokens + DCT_MAX_VALUE;
  const vp9_extra_bit *const e = vp9_extra_bits;

  int i = -DCT_MAX_VALUE;
  int sign = 1;

  do {
    if (!i)
      sign = 0;

    {
      const int a = sign ? -i : i;
      int eb = sign;

      if (a > 4) {
        int j = 4;

        while (++j < 11  &&  e[j].base_val <= a) {}

        t[i].token = --j;
        eb |= (a - e[j].base_val) << 1;
      } else
        t[i].token = a;

      t[i].extra = eb;
    }

    // initialize the cost for extra bits for all possible coefficient value.
    {
      int cost = 0;
      const vp9_extra_bit *p = vp9_extra_bits + t[i].token;

      if (p->base_val) {
        const int extra = t[i].extra;
        const int length = p->len;

        if (length)
          cost += treed_cost(p->tree, p->prob, extra >> 1, length);

        cost += vp9_cost_bit(vp9_prob_half, extra & 1); /* sign */
        dct_value_cost[i + DCT_MAX_VALUE] = cost;
      }

    }

  } while (++i < DCT_MAX_VALUE);

  vp9_dct_value_tokens_ptr = dct_value_tokens + DCT_MAX_VALUE;
  vp9_dct_value_cost_ptr   = dct_value_cost + DCT_MAX_VALUE;
}

struct tokenize_b_args {
  VP9_COMP *cpi;
  MACROBLOCKD *xd;
  TOKENEXTRA **tp;
  TX_SIZE tx_size;
};

static void set_entropy_context_b(int plane, int block, BLOCK_SIZE plane_bsize,
                                  TX_SIZE tx_size, void *arg) {
  struct tokenize_b_args* const args = arg;
  MACROBLOCKD *const xd = args->xd;
  struct macroblockd_plane *pd = &xd->plane[plane];
  int aoff, loff;
  txfrm_block_to_raster_xy(plane_bsize, tx_size, block, &aoff, &loff);
  set_contexts(xd, pd, plane_bsize, tx_size, pd->eobs[block] > 0, aoff, loff);
}

static void tokenize_b(int plane, int block, BLOCK_SIZE plane_bsize,
                       TX_SIZE tx_size, void *arg) {
  struct tokenize_b_args* const args = arg;
  VP9_COMP *cpi = args->cpi;
  MACROBLOCKD *xd = args->xd;
  TOKENEXTRA **tp = args->tp;
  struct macroblockd_plane *pd = &xd->plane[plane];
  MB_MODE_INFO *mbmi = &xd->mode_info_context->mbmi;
  int pt; /* near block/prev token context index */
  int c = 0, rc = 0;
  TOKENEXTRA *t = *tp;        /* store tokens starting here */
  const int eob = pd->eobs[block];
  const PLANE_TYPE type = pd->plane_type;
  const int16_t *qcoeff_ptr = BLOCK_OFFSET(pd->qcoeff, block);

  const int segment_id = mbmi->segment_id;
  const int16_t *scan, *nb;
  vp9_coeff_count *const counts = cpi->coef_counts[tx_size];
  vp9_coeff_probs_model *const coef_probs = cpi->common.fc.coef_probs[tx_size];
  const int ref = is_inter_block(mbmi);
  uint8_t token_cache[1024];
  const uint8_t *band_translate;
  ENTROPY_CONTEXT *A, *L;
  const int seg_eob = get_tx_eob(&cpi->common.seg, segment_id, tx_size);
  int aoff, loff;
  txfrm_block_to_raster_xy(plane_bsize, tx_size, block, &aoff, &loff);

  A = pd->above_context + aoff;
  L = pd->left_context + loff;

  assert((!type && !plane) || (type && plane));

  pt = get_entropy_context(xd, tx_size, type, block, A, L,
                           &scan, &band_translate);
  nb = vp9_get_coef_neighbors_handle(scan);
  c = 0;
  do {
    const int band = get_coef_band(band_translate, c);
    int token;
    int v = 0;
    rc = scan[c];
    if (c)
      pt = get_coef_context(nb, token_cache, c);
    if (c < eob) {
      v = qcoeff_ptr[rc];
      assert(-DCT_MAX_VALUE <= v  &&  v < DCT_MAX_VALUE);

      t->extra = vp9_dct_value_tokens_ptr[v].extra;
      token    = vp9_dct_value_tokens_ptr[v].token;
    } else {
      token = DCT_EOB_TOKEN;
    }

    t->token = token;
    t->context_tree = coef_probs[type][ref][band][pt];
    t->skip_eob_node = (c > 0) && (token_cache[scan[c - 1]] == 0);

    assert(vp9_coef_encodings[t->token].len - t->skip_eob_node > 0);

    ++counts[type][ref][band][pt][token];
    if (!t->skip_eob_node)
      ++cpi->common.counts.eob_branch[tx_size][type][ref][band][pt];

    token_cache[rc] = vp9_pt_energy_class[token];
    ++t;
  } while (c < eob && ++c < seg_eob);

  *tp = t;

  set_contexts(xd, pd, plane_bsize, tx_size, c > 0, aoff, loff);
}

struct is_skippable_args {
  MACROBLOCKD *xd;
  int *skippable;
};

static void is_skippable(int plane, int block,
                         BLOCK_SIZE plane_bsize, TX_SIZE tx_size,
                         void *argv) {
  struct is_skippable_args *args = argv;
  args->skippable[0] &= (!args->xd->plane[plane].eobs[block]);
}

int vp9_sb_is_skippable(MACROBLOCKD *xd, BLOCK_SIZE bsize) {
  int result = 1;
  struct is_skippable_args args = {xd, &result};
  foreach_transformed_block(xd, bsize, is_skippable, &args);
  return result;
}

int vp9_is_skippable_in_plane(MACROBLOCKD *xd, BLOCK_SIZE bsize,
                              int plane) {
  int result = 1;
  struct is_skippable_args args = {xd, &result};
  foreach_transformed_block_in_plane(xd, bsize, plane, is_skippable, &args);
  return result;
}

void vp9_tokenize_sb(VP9_COMP *cpi, TOKENEXTRA **t, int dry_run,
                     BLOCK_SIZE bsize) {
  VP9_COMMON *const cm = &cpi->common;
  MACROBLOCKD *const xd = &cpi->mb.e_mbd;
  MB_MODE_INFO *const mbmi = &xd->mode_info_context->mbmi;
  TOKENEXTRA *t_backup = *t;
  const int mb_skip_context = vp9_get_pred_context_mbskip(xd);
  const int skip_inc = !vp9_segfeature_active(&cm->seg, mbmi->segment_id,
                                              SEG_LVL_SKIP);
  struct tokenize_b_args arg = {cpi, xd, t, mbmi->txfm_size};

  mbmi->skip_coeff = vp9_sb_is_skippable(xd, bsize);
  if (mbmi->skip_coeff) {
    if (!dry_run)
      cm->counts.mbskip[mb_skip_context][1] += skip_inc;
    reset_skip_context(xd, bsize);
    if (dry_run)
      *t = t_backup;
    return;
  }

  if (!dry_run) {
    cm->counts.mbskip[mb_skip_context][0] += skip_inc;
    foreach_transformed_block(xd, bsize, tokenize_b, &arg);
  } else {
    foreach_transformed_block(xd, bsize, set_entropy_context_b, &arg);
    *t = t_backup;
  }
}

#ifdef ENTROPY_STATS
void init_context_counters(void) {
  FILE *f = fopen("context.bin", "rb");
  if (!f) {
    vp9_zero(context_counters);
  } else {
    fread(context_counters, sizeof(context_counters), 1, f);
    fclose(f);
  }

  f = fopen("treeupdate.bin", "rb");
  if (!f) {
    vpx_memset(tree_update_hist, 0, sizeof(tree_update_hist));
  } else {
    fread(tree_update_hist, sizeof(tree_update_hist), 1, f);
    fclose(f);
  }
}

static void print_counter(FILE *f, vp9_coeff_accum *context_counters,
                          int block_types, const char *header) {
  int type, ref, band, pt, t;

  fprintf(f, "static const vp9_coeff_count %s = {\n", header);

#define Comma(X) (X ? "," : "")
  type = 0;
  do {
    ref = 0;
    fprintf(f, "%s\n  { /* block Type %d */", Comma(type), type);
    do {
      fprintf(f, "%s\n    { /* %s */", Comma(type), ref ? "Inter" : "Intra");
      band = 0;
      do {
        fprintf(f, "%s\n      { /* Coeff Band %d */", Comma(band), band);
        pt = 0;
        do {
          fprintf(f, "%s\n        {", Comma(pt));

          t = 0;
          do {
            const int64_t x = context_counters[type][ref][band][pt][t];
            const int y = (int) x;

            assert(x == (int64_t) y);  /* no overflow handling yet */
            fprintf(f, "%s %d", Comma(t), y);
          } while (++t < 1 + MAX_ENTROPY_TOKENS);
          fprintf(f, "}");
        } while (++pt < PREV_COEF_CONTEXTS);
        fprintf(f, "\n      }");
      } while (++band < COEF_BANDS);
      fprintf(f, "\n    }");
    } while (++ref < REF_TYPES);
    fprintf(f, "\n  }");
  } while (++type < block_types);
  fprintf(f, "\n};\n");
}

static void print_probs(FILE *f, vp9_coeff_accum *context_counters,
                        int block_types, const char *header) {
  int type, ref, band, pt, t;

  fprintf(f, "static const vp9_coeff_probs %s = {", header);

  type = 0;
#define Newline(x, spaces) (x ? " " : "\n" spaces)
  do {
    fprintf(f, "%s%s{ /* block Type %d */",
            Comma(type), Newline(type, "  "), type);
    ref = 0;
    do {
      fprintf(f, "%s%s{ /* %s */",
              Comma(band), Newline(band, "    "), ref ? "Inter" : "Intra");
      band = 0;
      do {
        fprintf(f, "%s%s{ /* Coeff Band %d */",
                Comma(band), Newline(band, "      "), band);
        pt = 0;
        do {
          unsigned int branch_ct[ENTROPY_NODES][2];
          unsigned int coef_counts[MAX_ENTROPY_TOKENS + 1];
          vp9_prob coef_probs[ENTROPY_NODES];

          if (pt >= 3 && band == 0)
            break;
          for (t = 0; t < MAX_ENTROPY_TOKENS + 1; ++t)
            coef_counts[t] = context_counters[type][ref][band][pt][t];
          vp9_tree_probs_from_distribution(vp9_coef_tree, coef_probs,
                                           branch_ct, coef_counts, 0);
          branch_ct[0][1] = coef_counts[MAX_ENTROPY_TOKENS] - branch_ct[0][0];
          coef_probs[0] = get_binary_prob(branch_ct[0][0], branch_ct[0][1]);
          fprintf(f, "%s\n      {", Comma(pt));

          t = 0;
          do {
            fprintf(f, "%s %3d", Comma(t), coef_probs[t]);
          } while (++t < ENTROPY_NODES);

          fprintf(f, " }");
        } while (++pt < PREV_COEF_CONTEXTS);
        fprintf(f, "\n      }");
      } while (++band < COEF_BANDS);
      fprintf(f, "\n    }");
    } while (++ref < REF_TYPES);
    fprintf(f, "\n  }");
  } while (++type < block_types);
  fprintf(f, "\n};\n");
}

void print_context_counters() {
  FILE *f = fopen("vp9_context.c", "w");

  fprintf(f, "#include \"vp9_entropy.h\"\n");
  fprintf(f, "\n/* *** GENERATED FILE: DO NOT EDIT *** */\n\n");

  /* print counts */
  print_counter(f, context_counters[TX_4X4], BLOCK_TYPES,
                "vp9_default_coef_counts_4x4[BLOCK_TYPES]");
  print_counter(f, context_counters[TX_8X8], BLOCK_TYPES,
                "vp9_default_coef_counts_8x8[BLOCK_TYPES]");
  print_counter(f, context_counters[TX_16X16], BLOCK_TYPES,
                "vp9_default_coef_counts_16x16[BLOCK_TYPES]");
  print_counter(f, context_counters[TX_32X32], BLOCK_TYPES,
                "vp9_default_coef_counts_32x32[BLOCK_TYPES]");

  /* print coefficient probabilities */
  print_probs(f, context_counters[TX_4X4], BLOCK_TYPES,
              "default_coef_probs_4x4[BLOCK_TYPES]");
  print_probs(f, context_counters[TX_8X8], BLOCK_TYPES,
              "default_coef_probs_8x8[BLOCK_TYPES]");
  print_probs(f, context_counters[TX_16X16], BLOCK_TYPES,
              "default_coef_probs_16x16[BLOCK_TYPES]");
  print_probs(f, context_counters[TX_32X32], BLOCK_TYPES,
              "default_coef_probs_32x32[BLOCK_TYPES]");

  fclose(f);

  f = fopen("context.bin", "wb");
  fwrite(context_counters, sizeof(context_counters), 1, f);
  fclose(f);
}
#endif

void vp9_tokenize_initialize() {
  fill_value_tokens();
}
