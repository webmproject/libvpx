/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */


#include "vp8/common/header.h"
#include "encodemv.h"
#include "vp8/common/entropymode.h"
#include "vp8/common/findnearmv.h"
#include "mcomp.h"
#include "vp8/common/systemdependent.h"
#include <assert.h>
#include <stdio.h>
#include <limits.h>
#include "vp8/common/pragmas.h"
#include "vpx/vpx_encoder.h"
#include "vpx_mem/vpx_mem.h"
#include "bitstream.h"
#include "segmentation.h"

#include "vp8/common/seg_common.h"
#include "vp8/common/pred_common.h"
#include "vp8/common/entropy.h"
#include "vp8/encoder/encodemv.h"

#if CONFIG_NEWBESTREFMV
#include "vp8/common/mvref_common.h"
#endif

#if defined(SECTIONBITS_OUTPUT)
unsigned __int64 Sectionbits[500];
#endif

//int final_packing = 0;

#ifdef ENTROPY_STATS
int intra_mode_stats [VP8_BINTRAMODES] [VP8_BINTRAMODES] [VP8_BINTRAMODES];
unsigned int tree_update_hist [BLOCK_TYPES]
                              [COEF_BANDS]
                              [PREV_COEF_CONTEXTS]
                              [ENTROPY_NODES][2];
unsigned int hybrid_tree_update_hist [BLOCK_TYPES]
                                     [COEF_BANDS]
                                     [PREV_COEF_CONTEXTS]
                                     [ENTROPY_NODES][2];
unsigned int tree_update_hist_8x8 [BLOCK_TYPES_8X8]
                                  [COEF_BANDS]
                                  [PREV_COEF_CONTEXTS]
                                  [ENTROPY_NODES] [2];
unsigned int hybrid_tree_update_hist_8x8 [BLOCK_TYPES_8X8]
                                         [COEF_BANDS]
                                         [PREV_COEF_CONTEXTS]
                                         [ENTROPY_NODES] [2];
unsigned int tree_update_hist_16x16 [BLOCK_TYPES_16X16]
                                    [COEF_BANDS]
                                    [PREV_COEF_CONTEXTS]
                                    [ENTROPY_NODES] [2];
unsigned int hybrid_tree_update_hist_16x16 [BLOCK_TYPES_16X16]
                                           [COEF_BANDS]
                                           [PREV_COEF_CONTEXTS]
                                           [ENTROPY_NODES] [2];

extern unsigned int active_section;
#endif

#ifdef MODE_STATS
int count_mb_seg[4] = { 0, 0, 0, 0 };
#endif

#define vp8_cost_upd  ((int)(vp8_cost_one(upd) - vp8_cost_zero(upd)) >> 8)
#define vp8_cost_upd256  ((int)(vp8_cost_one(upd) - vp8_cost_zero(upd)))

#define SEARCH_NEWP
static int update_bits[255];

static void compute_update_table() {
  int i;
  for (i = 0; i < 255; i++)
    update_bits[i] = vp8_count_term_subexp(i, SUBEXP_PARAM, 255);
}

static int split_index(int i, int n, int modulus) {
  int max1 = (n - 1 - modulus / 2) / modulus + 1;
  if (i % modulus == modulus / 2) i = i / modulus;
  else i = max1 + i - (i + modulus - modulus / 2) / modulus;
  return i;
}

static int remap_prob(int v, int m) {
  const int n = 256;
  const int modulus = MODULUS_PARAM;
  int i;
  if ((m << 1) <= n)
    i = recenter_nonneg(v, m) - 1;
  else
    i = recenter_nonneg(n - 1 - v, n - 1 - m) - 1;

  i = split_index(i, n - 1, modulus);
  return i;
}

static void write_prob_diff_update(vp8_writer *const bc,
                                   vp8_prob newp, vp8_prob oldp) {
  int delp = remap_prob(newp, oldp);
  vp8_encode_term_subexp(bc, delp, SUBEXP_PARAM, 255);
}

static int prob_diff_update_cost(vp8_prob newp, vp8_prob oldp) {
  int delp = remap_prob(newp, oldp);
  return update_bits[delp] * 256;
}

#if CONFIG_NEW_MVREF
// Estimate the cost of each coding the vector using each reference candidate
unsigned int pick_best_mv_ref( MACROBLOCK *x,
                               int_mv target_mv,
                               int_mv * mv_ref_list,
                               int_mv * best_ref ) {

  int i;
  int best_index = 0;
  int cost, cost2;
  int index_cost[MAX_MV_REFS];
  MACROBLOCKD *xd = &x->e_mbd;

  /*unsigned int distance, distance2;

  distance = mv_distance(&target_mv, &mv_ref_list[0]);

  for (i = 1; i < MAX_MV_REFS; ++i ) {
    distance2 =
      mv_distance(&target_mv, &mv_ref_list[i]);
    if (distance2 < distance) {
      distance = distance2;
      best_index = i;
    }
  }*/

  // For now estimate the cost of selecting a given ref index
  // as index * 1 bits (but here 1 bit is scaled to 256)
  for (i = 0; i < MAX_MV_REFS; ++i ) {
    index_cost[i] = i << 8;
  }
  index_cost[0] = vp8_cost_zero(205);
  index_cost[1] = vp8_cost_zero(40);
  index_cost[2] = vp8_cost_zero(8);
  index_cost[3] = vp8_cost_zero(2);

  cost = index_cost[0] +
         vp8_mv_bit_cost(&target_mv,
                         &mv_ref_list[0],
                         XMVCOST, 96,
                         xd->allow_high_precision_mv);


  //for (i = 1; i < MAX_MV_REFS; ++i ) {
  for (i = 1; i < 4; ++i ) {
    cost2 = index_cost[i] +
            vp8_mv_bit_cost(&target_mv,
                            &mv_ref_list[i],
                            XMVCOST, 96,
                            xd->allow_high_precision_mv);

    if (cost2 < cost) {
      cost = cost2;
      best_index = i;
    }
  }

  (*best_ref).as_int = mv_ref_list[best_index].as_int;

  return best_index;
}
#endif

static void update_mode(
  vp8_writer *const bc,
  int n,
  vp8_token tok               [/* n */],
  vp8_tree tree,
  vp8_prob Pnew               [/* n-1 */],
  vp8_prob Pcur               [/* n-1 */],
  unsigned int bct            [/* n-1 */] [2],
  const unsigned int num_events[/* n */]
) {
  unsigned int new_b = 0, old_b = 0;
  int i = 0;

  vp8_tree_probs_from_distribution(
    n--, tok, tree,
    Pnew, bct, num_events,
    256, 1
  );

  do {
    new_b += vp8_cost_branch(bct[i], Pnew[i]);
    old_b += vp8_cost_branch(bct[i], Pcur[i]);
  } while (++i < n);

  if (new_b + (n << 8) < old_b) {
    int i = 0;

    vp8_write_bit(bc, 1);

    do {
      const vp8_prob p = Pnew[i];

      vp8_write_literal(bc, Pcur[i] = p ? p : 1, 8);
    } while (++i < n);
  } else
    vp8_write_bit(bc, 0);
}

static void update_mbintra_mode_probs(VP8_COMP* const cpi,
                                      vp8_writer* const bc) {
  VP8_COMMON *const cm = &cpi->common;

  {
    vp8_prob Pnew   [VP8_YMODES - 1];
    unsigned int bct [VP8_YMODES - 1] [2];

    update_mode(
      bc, VP8_YMODES, vp8_ymode_encodings, vp8_ymode_tree,
      Pnew, cm->fc.ymode_prob, bct, (unsigned int *)cpi->ymode_count
    );
  }
}

static int get_prob(int num, int den) {
  int p;
  if (den <= 0)
    return 128;
  p = (num * 255 + (den >> 1)) / den;
  if (p > 255)
    return 255;
  else if (p < 1)
    return 1;
  return p;
}

static int get_binary_prob(int n0, int n1) {
  return get_prob(n0, n0 + n1);
}

void update_skip_probs(VP8_COMP *cpi) {
  VP8_COMMON *const pc = &cpi->common;
  int prob_skip_false[3] = {0, 0, 0};
  int k;

  for (k = 0; k < MBSKIP_CONTEXTS; ++k) {
    pc->mbskip_pred_probs[k] = get_binary_prob(cpi->skip_false_count[k],
                                               cpi->skip_true_count[k]);
  }
}

#if CONFIG_SWITCHABLE_INTERP
void update_switchable_interp_probs(VP8_COMP *cpi, vp8_writer* const bc) {
  VP8_COMMON *const pc = &cpi->common;
  unsigned int branch_ct[32][2];
  int i, j;
  for (j = 0; j <= VP8_SWITCHABLE_FILTERS; ++j) {
  //for (j = 0; j <= 0; ++j) {
/*
    if (!cpi->dummy_packing)
#if VP8_SWITCHABLE_FILTERS == 3
      printf("HELLO %d %d %d\n", cpi->switchable_interp_count[j][0],
             cpi->switchable_interp_count[j][1], cpi->switchable_interp_count[j][2]);
#else
      printf("HELLO %d %d\n", cpi->switchable_interp_count[j][0],
             cpi->switchable_interp_count[j][1]);
#endif
*/
    vp8_tree_probs_from_distribution(
        VP8_SWITCHABLE_FILTERS,
        vp8_switchable_interp_encodings, vp8_switchable_interp_tree,
        pc->fc.switchable_interp_prob[j], branch_ct, cpi->switchable_interp_count[j],
        256, 1
        );
    for (i = 0; i < VP8_SWITCHABLE_FILTERS - 1; ++i) {
      if (pc->fc.switchable_interp_prob[j][i] < 1)
        pc->fc.switchable_interp_prob[j][i] = 1;
      vp8_write_literal(bc, pc->fc.switchable_interp_prob[j][i], 8);
/*
      if (!cpi->dummy_packing)
#if VP8_SWITCHABLE_FILTERS == 3
        printf("Probs %d %d [%d]\n",
               pc->fc.switchable_interp_prob[j][0],
               pc->fc.switchable_interp_prob[j][1], pc->frame_type);
#else
        printf("Probs %d [%d]\n", pc->fc.switchable_interp_prob[j][0],
               pc->frame_type);
#endif
*/
    }
  }
  /*
  if (!cpi->dummy_packing)
#if VP8_SWITCHABLE_FILTERS == 3
    printf("Probs %d %d [%d]\n",
           pc->fc.switchable_interp_prob[0], pc->fc.switchable_interp_prob[1], pc->frame_type);
#else
    printf("Probs %d [%d]\n", pc->fc.switchable_interp_prob[0], pc->frame_type);
#endif
  */
}
#endif

// This function updates the reference frame prediction stats
static void update_refpred_stats(VP8_COMP *cpi) {
  VP8_COMMON *const cm = &cpi->common;
  int i;
  int tot_count;
  vp8_prob new_pred_probs[PREDICTION_PROBS];
  int old_cost, new_cost;

  // Set the prediction probability structures to defaults
  if (cm->frame_type == KEY_FRAME) {
    // Set the prediction probabilities to defaults
    cm->ref_pred_probs[0] = 120;
    cm->ref_pred_probs[1] = 80;
    cm->ref_pred_probs[2] = 40;

    vpx_memset(cpi->ref_pred_probs_update, 0,
               sizeof(cpi->ref_pred_probs_update));
  } else {
    // From the prediction counts set the probabilities for each context
    for (i = 0; i < PREDICTION_PROBS; i++) {
      new_pred_probs[i] = get_binary_prob(cpi->ref_pred_count[i][0],
                                          cpi->ref_pred_count[i][1]);

      // Decide whether or not to update the reference frame probs.
      // Returned costs are in 1/256 bit units.
      old_cost =
        (cpi->ref_pred_count[i][0] * vp8_cost_zero(cm->ref_pred_probs[i])) +
        (cpi->ref_pred_count[i][1] * vp8_cost_one(cm->ref_pred_probs[i]));

      new_cost =
        (cpi->ref_pred_count[i][0] * vp8_cost_zero(new_pred_probs[i])) +
        (cpi->ref_pred_count[i][1] * vp8_cost_one(new_pred_probs[i]));

      // Cost saving must be >= 8 bits (2048 in these units)
      if ((old_cost - new_cost) >= 2048) {
        cpi->ref_pred_probs_update[i] = 1;
        cm->ref_pred_probs[i] = new_pred_probs[i];
      } else
        cpi->ref_pred_probs_update[i] = 0;

    }
  }
}

static void write_ymode(vp8_writer *bc, int m, const vp8_prob *p) {
  vp8_write_token(bc, vp8_ymode_tree, p, vp8_ymode_encodings + m);
}

static void kfwrite_ymode(vp8_writer *bc, int m, const vp8_prob *p) {
  vp8_write_token(bc, vp8_kf_ymode_tree, p, vp8_kf_ymode_encodings + m);
}

#if CONFIG_SUPERBLOCKS
static void sb_kfwrite_ymode(vp8_writer *bc, int m, const vp8_prob *p) {
  vp8_write_token(bc, vp8_uv_mode_tree, p, vp8_sb_kf_ymode_encodings + m);
}
#endif

static void write_i8x8_mode(vp8_writer *bc, int m, const vp8_prob *p) {
  vp8_write_token(bc, vp8_i8x8_mode_tree, p, vp8_i8x8_mode_encodings + m);
}

static void write_uv_mode(vp8_writer *bc, int m, const vp8_prob *p) {
  vp8_write_token(bc, vp8_uv_mode_tree, p, vp8_uv_mode_encodings + m);
}


static void write_bmode(vp8_writer *bc, int m, const vp8_prob *p) {
  vp8_write_token(bc, vp8_bmode_tree, p, vp8_bmode_encodings + m);
}

static void write_split(vp8_writer *bc, int x, const vp8_prob *p) {
  vp8_write_token(
    bc, vp8_mbsplit_tree, p, vp8_mbsplit_encodings + x
  );
}

static int prob_update_savings(const unsigned int *ct,
                               const vp8_prob oldp, const vp8_prob newp,
                               const vp8_prob upd) {
  const int old_b = vp8_cost_branch256(ct, oldp);
  const int new_b = vp8_cost_branch256(ct, newp);
  const int update_b = 2048 + vp8_cost_upd256;
  return (old_b - new_b - update_b);
}

static int prob_diff_update_savings(const unsigned int *ct,
                                    const vp8_prob oldp, const vp8_prob newp,
                                    const vp8_prob upd) {
  const int old_b = vp8_cost_branch256(ct, oldp);
  const int new_b = vp8_cost_branch256(ct, newp);
  const int update_b = (newp == oldp ? 0 :
                        prob_diff_update_cost(newp, oldp) + vp8_cost_upd256);
  return (old_b - new_b - update_b);
}

static int prob_diff_update_savings_search(const unsigned int *ct,
                                           const vp8_prob oldp, vp8_prob *bestp,
                                           const vp8_prob upd) {
  const int old_b = vp8_cost_branch256(ct, oldp);
  int new_b, update_b, savings, bestsavings, step;
  vp8_prob newp, bestnewp;

  bestsavings = 0;
  bestnewp = oldp;

  step = (*bestp > oldp ? -1 : 1);
  for (newp = *bestp; newp != oldp; newp += step) {
    new_b = vp8_cost_branch256(ct, newp);
    update_b = prob_diff_update_cost(newp, oldp) + vp8_cost_upd256;
    savings = old_b - new_b - update_b;
    if (savings > bestsavings) {
      bestsavings = savings;
      bestnewp = newp;
    }
  }
  *bestp = bestnewp;
  return bestsavings;
}

static void pack_mb_tokens(vp8_writer* const bc,
                           TOKENEXTRA **tp,
                           const TOKENEXTRA *const stop) {
  unsigned int split;
  unsigned int shift;
  int count = bc->count;
  unsigned int range = bc->range;
  unsigned int lowvalue = bc->lowvalue;
  TOKENEXTRA *p = *tp;

  while (p < stop) {
    const int t = p->Token;
    vp8_token *const a = vp8_coef_encodings + t;
    const vp8_extra_bit_struct *const b = vp8_extra_bits + t;
    int i = 0;
    const unsigned char *pp = p->context_tree;
    int v = a->value;
    int n = a->Len;

    if (t == EOSB_TOKEN)
    {
      ++p;
      break;
    }

    /* skip one or two nodes */
    if (p->skip_eob_node) {
      n -= p->skip_eob_node;
      i = 2 * p->skip_eob_node;
    }

    do {
      const int bb = (v >> --n) & 1;
      split = 1 + (((range - 1) * pp[i >> 1]) >> 8);
      i = vp8_coef_tree[i + bb];

      if (bb) {
        lowvalue += split;
        range = range - split;
      } else {
        range = split;
      }

      shift = vp8_norm[range];
      range <<= shift;
      count += shift;

      if (count >= 0) {
        int offset = shift - count;

        if ((lowvalue << (offset - 1)) & 0x80000000) {
          int x = bc->pos - 1;

          while (x >= 0 && bc->buffer[x] == 0xff) {
            bc->buffer[x] = (unsigned char)0;
            x--;
          }

          bc->buffer[x] += 1;
        }

        bc->buffer[bc->pos++] = (lowvalue >> (24 - offset));
        lowvalue <<= offset;
        shift = count;
        lowvalue &= 0xffffff;
        count -= 8;
      }

      lowvalue <<= shift;
    } while (n);


    if (b->base_val) {
      const int e = p->Extra, L = b->Len;

      if (L) {
        const unsigned char *pp = b->prob;
        int v = e >> 1;
        int n = L;              /* number of bits in v, assumed nonzero */
        int i = 0;

        do {
          const int bb = (v >> --n) & 1;
          split = 1 + (((range - 1) * pp[i >> 1]) >> 8);
          i = b->tree[i + bb];

          if (bb) {
            lowvalue += split;
            range = range - split;
          } else {
            range = split;
          }

          shift = vp8_norm[range];
          range <<= shift;
          count += shift;

          if (count >= 0) {
            int offset = shift - count;

            if ((lowvalue << (offset - 1)) & 0x80000000) {
              int x = bc->pos - 1;

              while (x >= 0 && bc->buffer[x] == 0xff) {
                bc->buffer[x] = (unsigned char)0;
                x--;
              }

              bc->buffer[x] += 1;
            }

            bc->buffer[bc->pos++] = (lowvalue >> (24 - offset));
            lowvalue <<= offset;
            shift = count;
            lowvalue &= 0xffffff;
            count -= 8;
          }

          lowvalue <<= shift;
        } while (n);
      }


      {

        split = (range + 1) >> 1;

        if (e & 1) {
          lowvalue += split;
          range = range - split;
        } else {
          range = split;
        }

        range <<= 1;

        if ((lowvalue & 0x80000000)) {
          int x = bc->pos - 1;

          while (x >= 0 && bc->buffer[x] == 0xff) {
            bc->buffer[x] = (unsigned char)0;
            x--;
          }

          bc->buffer[x] += 1;

        }

        lowvalue  <<= 1;

        if (!++count) {
          count = -8;
          bc->buffer[bc->pos++] = (lowvalue >> 24);
          lowvalue &= 0xffffff;
        }
      }

    }
    ++p;
  }

  bc->count = count;
  bc->lowvalue = lowvalue;
  bc->range = range;
  *tp = p;
}

static void write_partition_size(unsigned char *cx_data, int size) {
  signed char csize;

  csize = size & 0xff;
  *cx_data = csize;
  csize = (size >> 8) & 0xff;
  *(cx_data + 1) = csize;
  csize = (size >> 16) & 0xff;
  *(cx_data + 2) = csize;

}

static void write_mv_ref
(
  vp8_writer *bc, MB_PREDICTION_MODE m, const vp8_prob *p
) {
#if CONFIG_DEBUG
  assert(NEARESTMV <= m  &&  m <= SPLITMV);
#endif
  vp8_write_token(bc, vp8_mv_ref_tree, p,
                  vp8_mv_ref_encoding_array - NEARESTMV + m);
}

#if CONFIG_SUPERBLOCKS
static void write_sb_mv_ref(vp8_writer *bc, MB_PREDICTION_MODE m,
                            const vp8_prob *p) {
#if CONFIG_DEBUG
  assert(NEARESTMV <= m  &&  m < SPLITMV);
#endif
  vp8_write_token(bc, vp8_sb_mv_ref_tree, p,
                  vp8_sb_mv_ref_encoding_array - NEARESTMV + m);
}
#endif

static void write_sub_mv_ref
(
  vp8_writer *bc, B_PREDICTION_MODE m, const vp8_prob *p
) {
#if CONFIG_DEBUG
  assert(LEFT4X4 <= m  &&  m <= NEW4X4);
#endif
  vp8_write_token(bc, vp8_sub_mv_ref_tree, p,
                  vp8_sub_mv_ref_encoding_array - LEFT4X4 + m);
}

#if CONFIG_NEWMVENTROPY
static void write_nmv(vp8_writer *bc, const MV *mv, const int_mv *ref,
                      const nmv_context *nmvc, int usehp) {
  MV e;
  e.row = mv->row - ref->as_mv.row;
  e.col = mv->col - ref->as_mv.col;

  vp8_encode_nmv(bc, &e, &ref->as_mv, nmvc);
  vp8_encode_nmv_fp(bc, &e, &ref->as_mv, nmvc, usehp);
}

#else

static void write_mv
(
  vp8_writer *bc, const MV *mv, const int_mv *ref, const MV_CONTEXT *mvc
) {
  MV e;
  e.row = mv->row - ref->as_mv.row;
  e.col = mv->col - ref->as_mv.col;

  vp8_encode_motion_vector(bc, &e, mvc);
}

static void write_mv_hp
(
  vp8_writer *bc, const MV *mv, const int_mv *ref, const MV_CONTEXT_HP *mvc
) {
  MV e;
  e.row = mv->row - ref->as_mv.row;
  e.col = mv->col - ref->as_mv.col;

  vp8_encode_motion_vector_hp(bc, &e, mvc);
}
#endif  /* CONFIG_NEWMVENTROPY */

// This function writes the current macro block's segnment id to the bitstream
// It should only be called if a segment map update is indicated.
static void write_mb_segid(vp8_writer *bc,
                           const MB_MODE_INFO *mi, const MACROBLOCKD *xd) {
  // Encode the MB segment id.
  if (xd->segmentation_enabled && xd->update_mb_segmentation_map) {
    switch (mi->segment_id) {
      case 0:
        vp8_write(bc, 0, xd->mb_segment_tree_probs[0]);
        vp8_write(bc, 0, xd->mb_segment_tree_probs[1]);
        break;
      case 1:
        vp8_write(bc, 0, xd->mb_segment_tree_probs[0]);
        vp8_write(bc, 1, xd->mb_segment_tree_probs[1]);
        break;
      case 2:
        vp8_write(bc, 1, xd->mb_segment_tree_probs[0]);
        vp8_write(bc, 0, xd->mb_segment_tree_probs[2]);
        break;
      case 3:
        vp8_write(bc, 1, xd->mb_segment_tree_probs[0]);
        vp8_write(bc, 1, xd->mb_segment_tree_probs[2]);
        break;

        // TRAP.. This should not happen
      default:
        vp8_write(bc, 0, xd->mb_segment_tree_probs[0]);
        vp8_write(bc, 0, xd->mb_segment_tree_probs[1]);
        break;
    }
  }
}

// This function encodes the reference frame
static void encode_ref_frame(vp8_writer *const bc,
                             VP8_COMMON *const cm,
                             MACROBLOCKD *xd,
                             int segment_id,
                             MV_REFERENCE_FRAME rf) {
  int seg_ref_active;
  int seg_ref_count = 0;
  seg_ref_active = segfeature_active(xd,
                                     segment_id,
                                     SEG_LVL_REF_FRAME);

  if (seg_ref_active) {
    seg_ref_count = check_segref(xd, segment_id, INTRA_FRAME) +
                    check_segref(xd, segment_id, LAST_FRAME) +
                    check_segref(xd, segment_id, GOLDEN_FRAME) +
                    check_segref(xd, segment_id, ALTREF_FRAME);
  }

  // If segment level coding of this signal is disabled...
  // or the segment allows multiple reference frame options
  if (!seg_ref_active || (seg_ref_count > 1)) {
    // Values used in prediction model coding
    unsigned char prediction_flag;
    vp8_prob pred_prob;
    MV_REFERENCE_FRAME pred_rf;

    // Get the context probability the prediction flag
    pred_prob = get_pred_prob(cm, xd, PRED_REF);

    // Get the predicted value.
    pred_rf = get_pred_ref(cm, xd);

    // Did the chosen reference frame match its predicted value.
    prediction_flag =
      (xd->mode_info_context->mbmi.ref_frame == pred_rf);

    set_pred_flag(xd, PRED_REF, prediction_flag);
    vp8_write(bc, prediction_flag, pred_prob);

    // If not predicted correctly then code value explicitly
    if (!prediction_flag) {
      vp8_prob mod_refprobs[PREDICTION_PROBS];

      vpx_memcpy(mod_refprobs,
                 cm->mod_refprobs[pred_rf], sizeof(mod_refprobs));

      // If segment coding enabled blank out options that cant occur by
      // setting the branch probability to 0.
      if (seg_ref_active) {
        mod_refprobs[INTRA_FRAME] *=
          check_segref(xd, segment_id, INTRA_FRAME);
        mod_refprobs[LAST_FRAME] *=
          check_segref(xd, segment_id, LAST_FRAME);
        mod_refprobs[GOLDEN_FRAME] *=
          (check_segref(xd, segment_id, GOLDEN_FRAME) *
           check_segref(xd, segment_id, ALTREF_FRAME));
      }

      if (mod_refprobs[0]) {
        vp8_write(bc, (rf != INTRA_FRAME), mod_refprobs[0]);
      }

      // Inter coded
      if (rf != INTRA_FRAME) {
        if (mod_refprobs[1]) {
          vp8_write(bc, (rf != LAST_FRAME), mod_refprobs[1]);
        }

        if (rf != LAST_FRAME) {
          if (mod_refprobs[2]) {
            vp8_write(bc, (rf != GOLDEN_FRAME), mod_refprobs[2]);
          }
        }
      }
    }
  }

  // if using the prediction mdoel we have nothing further to do because
  // the reference frame is fully coded by the segment
}

// Update the probabilities used to encode reference frame data
static void update_ref_probs(VP8_COMP *const cpi) {
  VP8_COMMON *const cm = &cpi->common;

  const int *const rfct = cpi->count_mb_ref_frame_usage;
  const int rf_intra = rfct[INTRA_FRAME];
  const int rf_inter = rfct[LAST_FRAME] +
                       rfct[GOLDEN_FRAME] + rfct[ALTREF_FRAME];

  cm->prob_intra_coded = get_binary_prob(rf_intra, rf_inter);
  cm->prob_last_coded = get_prob(rfct[LAST_FRAME], rf_inter);
  cm->prob_gf_coded = get_binary_prob(rfct[GOLDEN_FRAME], rfct[ALTREF_FRAME]);

  // Compute a modified set of probabilities to use when prediction of the
  // reference frame fails
  compute_mod_refprobs(cm);
}

static void pack_inter_mode_mvs(VP8_COMP *const cpi, vp8_writer *const bc) {
  int i;
  VP8_COMMON *const pc = &cpi->common;
#if CONFIG_NEWMVENTROPY
  const nmv_context *nmvc = &pc->fc.nmvc;
#else
  const MV_CONTEXT *mvc = pc->fc.mvc;
  const MV_CONTEXT_HP *mvc_hp = pc->fc.mvc_hp;
#endif
  MACROBLOCK *x = &cpi->mb;
  MACROBLOCKD *xd = &cpi->mb.e_mbd;
  MODE_INFO *m;
  MODE_INFO *prev_m;
  TOKENEXTRA *tok = cpi->tok;
  TOKENEXTRA *tok_end = tok + cpi->tok_count;

  const int mis = pc->mode_info_stride;
  int mb_row, mb_col;
  int row, col;

  // Values used in prediction model coding
  vp8_prob pred_prob;
  unsigned char prediction_flag;

  int row_delta[4] = { 0, +1,  0, -1};
  int col_delta[4] = { +1, -1, +1, +1};

  //final_packing = !cpi->dummy_packing;

  cpi->mb.partition_info = cpi->mb.pi;

  mb_row = 0;
  for (row = 0; row < pc->mb_rows; row += 2) {
    m = pc->mi + row * mis;
    prev_m = pc->prev_mi + row * mis;

    mb_col = 0;
    for (col = 0; col < pc->mb_cols; col += 2) {
      int i;

      // Process the 4 MBs in the order:
      // top-left, top-right, bottom-left, bottom-right
#if CONFIG_SUPERBLOCKS
      vp8_write(bc, m->mbmi.encoded_as_sb, pc->sb_coded);
#endif
      for (i = 0; i < 4; i++) {
        MB_MODE_INFO *mi;
        MV_REFERENCE_FRAME rf;
        MB_PREDICTION_MODE mode;
        int segment_id;

        int dy = row_delta[i];
        int dx = col_delta[i];
        int offset_extended = dy * mis + dx;

        if ((mb_row >= pc->mb_rows) || (mb_col >= pc->mb_cols)) {
          // MB lies outside frame, move on
          mb_row += dy;
          mb_col += dx;
          m += offset_extended;
          prev_m += offset_extended;
          cpi->mb.partition_info += offset_extended;
          continue;
        }

        mi = &m->mbmi;
        rf = mi->ref_frame;
        mode = mi->mode;
        segment_id = mi->segment_id;

        // Distance of Mb to the various image edges.
        // These specified to 8th pel as they are always compared to MV
        // values that are in 1/8th pel units
        xd->mb_to_left_edge = -((mb_col * 16) << 3);
        xd->mb_to_right_edge = ((pc->mb_cols - 1 - mb_col) * 16) << 3;
        xd->mb_to_top_edge = -((mb_row * 16)) << 3;
        xd->mb_to_bottom_edge = ((pc->mb_rows - 1 - mb_row) * 16) << 3;

        // Make sure the MacroBlockD mode info pointer is set correctly
        xd->mode_info_context = m;
        xd->prev_mode_info_context = prev_m;

#ifdef ENTROPY_STATS
        active_section = 9;
#endif

        if (cpi->mb.e_mbd.update_mb_segmentation_map) {
          // Is temporal coding of the segment map enabled
          if (pc->temporal_update) {
            prediction_flag = get_pred_flag(xd, PRED_SEG_ID);
            pred_prob = get_pred_prob(pc, xd, PRED_SEG_ID);

            // Code the segment id prediction flag for this mb
            vp8_write(bc, prediction_flag, pred_prob);

            // If the mb segment id wasn't predicted code explicitly
            if (!prediction_flag)
              write_mb_segid(bc, mi, &cpi->mb.e_mbd);
          } else {
            // Normal unpredicted coding
            write_mb_segid(bc, mi, &cpi->mb.e_mbd);
          }
        }

        if (pc->mb_no_coeff_skip &&
            (!segfeature_active(xd, segment_id, SEG_LVL_EOB) ||
             (get_segdata(xd, segment_id, SEG_LVL_EOB) != 0))) {
          int skip_coeff = mi->mb_skip_coeff;
#if CONFIG_SUPERBLOCKS
          if (mi->encoded_as_sb) {
            skip_coeff &= m[1].mbmi.mb_skip_coeff;
            skip_coeff &= m[mis].mbmi.mb_skip_coeff;
            skip_coeff &= m[mis + 1].mbmi.mb_skip_coeff;
          }
#endif
          vp8_encode_bool(bc, skip_coeff,
                          get_pred_prob(pc, xd, PRED_MBSKIP));
        }

        // Encode the reference frame.
        encode_ref_frame(bc, pc, xd, segment_id, rf);

        if (rf == INTRA_FRAME) {
#ifdef ENTROPY_STATS
          active_section = 6;
#endif

          // TODO(rbultje) write using SB tree structure

          if (!segfeature_active(xd, segment_id, SEG_LVL_MODE)) {
            write_ymode(bc, mode, pc->fc.ymode_prob);
          }

          if (mode == B_PRED) {
            int j = 0;
#if CONFIG_COMP_INTRA_PRED
            int uses_second =
              m->bmi[0].as_mode.second !=
              (B_PREDICTION_MODE)(B_DC_PRED - 1);
            vp8_write(bc, uses_second, 128);
#endif
            do {
#if CONFIG_COMP_INTRA_PRED
              B_PREDICTION_MODE mode2 = m->bmi[j].as_mode.second;
#endif
              write_bmode(bc, m->bmi[j].as_mode.first,
                          pc->fc.bmode_prob);
              /*
              if (!cpi->dummy_packing) {
                int p;
                for (p = 0; p < VP8_BINTRAMODES - 1; ++p)
                  printf(" %d", pc->fc.bmode_prob[p]);
                printf("\nbmode[%d][%d]: %d\n", pc->current_video_frame, j, m->bmi[j].as_mode.first);
              }
              */
#if CONFIG_COMP_INTRA_PRED
              if (uses_second) {
                write_bmode(bc, mode2, pc->fc.bmode_prob);
              }
#endif
            } while (++j < 16);
          }
          if (mode == I8X8_PRED) {
            write_i8x8_mode(bc, m->bmi[0].as_mode.first,
                            pc->fc.i8x8_mode_prob);
            write_i8x8_mode(bc, m->bmi[2].as_mode.first,
                            pc->fc.i8x8_mode_prob);
            write_i8x8_mode(bc, m->bmi[8].as_mode.first,
                            pc->fc.i8x8_mode_prob);
            write_i8x8_mode(bc, m->bmi[10].as_mode.first,
                            pc->fc.i8x8_mode_prob);
          } else {
            write_uv_mode(bc, mi->uv_mode,
                          pc->fc.uv_mode_prob[mode]);
          }
        } else {
          int_mv best_mv, best_second_mv;
          int ct[4];

          vp8_prob mv_ref_p [VP8_MVREFS - 1];

          {
            int_mv n1, n2;

            vp8_find_near_mvs(xd, m, prev_m, &n1, &n2, &best_mv, ct,
                              rf, cpi->common.ref_frame_sign_bias);
#if CONFIG_NEWBESTREFMV
            best_mv.as_int = mi->ref_mv.as_int;
#endif
            vp8_mv_ref_probs(&cpi->common, mv_ref_p, ct);

#ifdef ENTROPY_STATS
            accum_mv_refs(mode, ct);
#endif
          }

#ifdef ENTROPY_STATS
          active_section = 3;
#endif

          // Is the segment coding of mode enabled
          if (!segfeature_active(xd, segment_id, SEG_LVL_MODE)) {
#if CONFIG_SUPERBLOCKS
            if (mi->encoded_as_sb) {
              write_sb_mv_ref(bc, mode, mv_ref_p);
            } else
#endif
            {
              write_mv_ref(bc, mode, mv_ref_p);
            }
            vp8_accum_mv_refs(&cpi->common, mode, ct);
          }

#if CONFIG_PRED_FILTER
          // Is the prediction filter enabled
          if (mode >= NEARESTMV && mode < SPLITMV) {
            if (cpi->common.pred_filter_mode == 2)
              vp8_write(bc, mi->pred_filter_enabled,
                        pc->prob_pred_filter_off);
            else
              assert(mi->pred_filter_enabled ==
                     cpi->common.pred_filter_mode);
          }
#endif
#if CONFIG_SWITCHABLE_INTERP
          if (mode >= NEARESTMV && mode <= SPLITMV)
          {
            if (cpi->common.mcomp_filter_type == SWITCHABLE) {
              vp8_write_token(bc, vp8_switchable_interp_tree,
                              get_pred_probs(&cpi->common, xd, PRED_SWITCHABLE_INTERP),
                              vp8_switchable_interp_encodings +
                              vp8_switchable_interp_map[mi->interp_filter]);
              //if (!cpi->dummy_packing) printf("Reading: %d\n", mi->interp_filter);
            } else {
              assert (mi->interp_filter ==
                      cpi->common.mcomp_filter_type);
            }
          }
#endif
          if (mi->second_ref_frame &&
              (mode == NEWMV || mode == SPLITMV)) {
            int_mv n1, n2;

            vp8_find_near_mvs(xd, m,
                              prev_m,
                              &n1, &n2, &best_second_mv, ct,
                              mi->second_ref_frame,
                              cpi->common.ref_frame_sign_bias);
#if CONFIG_NEWBESTREFMV
            best_second_mv.as_int = mi->second_ref_mv.as_int;
#endif
          }

          // does the feature use compound prediction or not
          // (if not specified at the frame/segment level)
          if (cpi->common.comp_pred_mode == HYBRID_PREDICTION) {
            vp8_write(bc, mi->second_ref_frame != INTRA_FRAME,
                      get_pred_prob(pc, xd, PRED_COMP));
          }

          {
            switch (mode) { /* new, split require MVs */
              case NEWMV:
#ifdef ENTROPY_STATS
                active_section = 5;
#endif

#if 0 //CONFIG_NEW_MVREF
                {
                  unsigned int best_index;
                  /*find_mv_refs(xd, m, prev_m,
                               m->mbmi.ref_frame,
                               mi->ref_mvs[rf],
                               cpi->common.ref_frame_sign_bias );*/

                  best_index = pick_best_mv_ref(x, mi->mv[0],
                                                mi->ref_mvs[rf], &best_mv);
                  cpi->best_ref_index_counts[best_index]++;

                }
#endif
#if CONFIG_NEWMVENTROPY
                write_nmv(bc, &mi->mv[0].as_mv, &best_mv,
                          (const nmv_context*) nmvc,
                          xd->allow_high_precision_mv);
#else
                if (xd->allow_high_precision_mv) {
                  write_mv_hp(bc, &mi->mv[0].as_mv, &best_mv, mvc_hp);
                } else {
                  write_mv(bc, &mi->mv[0].as_mv, &best_mv, mvc);
                }
#endif

                if (mi->second_ref_frame) {
#if 0 //CONFIG_NEW_MVREF
                  unsigned int best_index;

                  /*find_mv_refs(xd, m, prev_m,
                               m->mbmi.second_ref_frame,
                               mi->ref_mvs[mi->second_ref_frame],
                               cpi->common.ref_frame_sign_bias );*/

                  best_index =
                    pick_best_mv_ref(x, mi->mv[1],
                                     mi->ref_mvs[mi->second_ref_frame],
                                     &best_second_mv);
                  cpi->best_ref_index_counts[best_index]++;
#endif
#if CONFIG_NEWMVENTROPY
                  write_nmv(bc, &mi->mv[1].as_mv, &best_second_mv,
                            (const nmv_context*) nmvc,
                            xd->allow_high_precision_mv);
#else
                  if (xd->allow_high_precision_mv) {
                    write_mv_hp(bc, &mi->mv[1].as_mv, &best_second_mv, mvc_hp);
                  } else {
                    write_mv(bc, &mi->mv[1].as_mv, &best_second_mv, mvc);
                  }
#endif
                }
                break;
              case SPLITMV: {
                int j = 0;

#ifdef MODE_STATS
                ++count_mb_seg [mi->partitioning];
#endif

                write_split(bc, mi->partitioning, cpi->common.fc.mbsplit_prob);
                cpi->mbsplit_count[mi->partitioning]++;

                do {
                  B_PREDICTION_MODE blockmode;
                  int_mv blockmv;
                  const int *const  L =
                    vp8_mbsplits [mi->partitioning];
                  int k = -1;  /* first block in subset j */
                  int mv_contz;
                  int_mv leftmv, abovemv;

                  blockmode = cpi->mb.partition_info->bmi[j].mode;
                  blockmv = cpi->mb.partition_info->bmi[j].mv;
#if CONFIG_DEBUG
                  while (j != L[++k])
                    if (k >= 16)
                      assert(0);
#else
                  while (j != L[++k]);
#endif
                  leftmv.as_int = left_block_mv(m, k);
                  abovemv.as_int = above_block_mv(m, k, mis);
                  mv_contz = vp8_mv_cont(&leftmv, &abovemv);

                  write_sub_mv_ref(bc, blockmode,
                                   cpi->common.fc.sub_mv_ref_prob [mv_contz]);
                  cpi->sub_mv_ref_count[mv_contz][blockmode - LEFT4X4]++;
                  if (blockmode == NEW4X4) {
#ifdef ENTROPY_STATS
                    active_section = 11;
#endif
#if CONFIG_NEWMVENTROPY
                    write_nmv(bc, &blockmv.as_mv, &best_mv,
                              (const nmv_context*) nmvc,
                              xd->allow_high_precision_mv);
#else
                    if (xd->allow_high_precision_mv) {
                      write_mv_hp(bc, &blockmv.as_mv, &best_mv,
                                  (const MV_CONTEXT_HP *) mvc_hp);
                    } else {
                      write_mv(bc, &blockmv.as_mv, &best_mv,
                               (const MV_CONTEXT *) mvc);
                    }
#endif

                    if (mi->second_ref_frame) {
#if CONFIG_NEWMVENTROPY
                      write_nmv(bc,
                                &cpi->mb.partition_info->bmi[j].second_mv.as_mv,
                                &best_second_mv,
                                (const nmv_context*) nmvc,
                                xd->allow_high_precision_mv);
#else
                      if (xd->allow_high_precision_mv) {
                        write_mv_hp(
                            bc,
                            &cpi->mb.partition_info->bmi[j].second_mv.as_mv,
                            &best_second_mv, (const MV_CONTEXT_HP *)mvc_hp);
                      } else {
                        write_mv(
                            bc,
                            &cpi->mb.partition_info->bmi[j].second_mv.as_mv,
                            &best_second_mv, (const MV_CONTEXT *) mvc);
                      }
#endif
                    }
                  }
                } while (++j < cpi->mb.partition_info->count);
              }
              break;
              default:
                break;
            }
          }
        }

#if CONFIG_TX_SELECT
        if (((rf == INTRA_FRAME && mode <= I8X8_PRED) ||
             (rf != INTRA_FRAME && mode != SPLITMV)) &&
            pc->txfm_mode == TX_MODE_SELECT &&
            !((pc->mb_no_coeff_skip && mi->mb_skip_coeff) ||
              (segfeature_active(xd, segment_id, SEG_LVL_EOB) &&
               get_segdata(xd, segment_id, SEG_LVL_EOB) == 0))) {
          TX_SIZE sz = mi->txfm_size;
          // FIXME(rbultje) code ternary symbol once all experiments are merged
          vp8_write(bc, sz != TX_4X4, pc->prob_tx[0]);
          if (sz != TX_4X4 && mode != I8X8_PRED)
            vp8_write(bc, sz != TX_8X8, pc->prob_tx[1]);
        }
#endif

#ifdef ENTROPY_STATS
        active_section = 1;
#endif
        assert(tok < tok_end);
        pack_mb_tokens(bc, &tok, tok_end);

#if CONFIG_SUPERBLOCKS
        if (m->mbmi.encoded_as_sb) {
          assert(!i);
          mb_col += 2;
          m += 2;
          cpi->mb.partition_info += 2;
          prev_m += 2;
          break;
        }
#endif

        // Next MB
        mb_row += dy;
        mb_col += dx;
        m += offset_extended;
        prev_m += offset_extended;
        cpi->mb.partition_info += offset_extended;
#if CONFIG_DEBUG
        assert((prev_m - cpi->common.prev_mip) == (m - cpi->common.mip));
        assert((prev_m - cpi->common.prev_mi) == (m - cpi->common.mi));
#endif
      }
    }

    // Next SB
    mb_row += 2;
    m += mis + (1 - (pc->mb_cols & 0x1));
    prev_m += mis + (1 - (pc->mb_cols & 0x1));
    cpi->mb.partition_info += mis + (1 - (pc->mb_cols & 0x1));
  }
}


static void write_mb_modes_kf(const VP8_COMMON  *c,
                              const MACROBLOCKD *xd,
                              const MODE_INFO   *m,
                              int                mode_info_stride,
                              vp8_writer *const  bc) {
  const int mis = mode_info_stride;
  int ym;
  int segment_id;

  ym = m->mbmi.mode;
  segment_id = m->mbmi.segment_id;

  if (xd->update_mb_segmentation_map) {
    write_mb_segid(bc, &m->mbmi, xd);
  }

  if (c->mb_no_coeff_skip &&
      (!segfeature_active(xd, segment_id, SEG_LVL_EOB) ||
       (get_segdata(xd, segment_id, SEG_LVL_EOB) != 0))) {
        int skip_coeff = m->mbmi.mb_skip_coeff;
#if CONFIG_SUPERBLOCKS
        if (m->mbmi.encoded_as_sb) {
          skip_coeff &= m[1].mbmi.mb_skip_coeff;
          skip_coeff &= m[mis].mbmi.mb_skip_coeff;
          skip_coeff &= m[mis + 1].mbmi.mb_skip_coeff;
        }
#endif
        vp8_encode_bool(bc, skip_coeff,
                    get_pred_prob(c, xd, PRED_MBSKIP));
  }

#if CONFIG_SUPERBLOCKS
  if (m->mbmi.encoded_as_sb) {
    sb_kfwrite_ymode(bc, ym,
                     c->sb_kf_ymode_prob[c->kf_ymode_probs_index]);
  } else
#endif
  {
    kfwrite_ymode(bc, ym,
                  c->kf_ymode_prob[c->kf_ymode_probs_index]);
  }

  if (ym == B_PRED) {
    const int mis = c->mode_info_stride;
    int i = 0;
#if CONFIG_COMP_INTRA_PRED
    int uses_second =
      m->bmi[0].as_mode.second !=
      (B_PREDICTION_MODE)(B_DC_PRED - 1);
    vp8_write(bc, uses_second, 128);
#endif
    do {
      const B_PREDICTION_MODE A = above_block_mode(m, i, mis);
      const B_PREDICTION_MODE L = left_block_mode(m, i);
      const int bm = m->bmi[i].as_mode.first;
#if CONFIG_COMP_INTRA_PRED
      const int bm2 = m->bmi[i].as_mode.second;
#endif

#ifdef ENTROPY_STATS
      ++intra_mode_stats [A] [L] [bm];
#endif

      write_bmode(bc, bm, c->kf_bmode_prob [A] [L]);
      // printf("    mode: %d\n", bm);
#if CONFIG_COMP_INTRA_PRED
      if (uses_second) {
        write_bmode(bc, bm2, c->kf_bmode_prob [A] [L]);
      }
#endif
    } while (++i < 16);
  }
  if (ym == I8X8_PRED) {
    write_i8x8_mode(bc, m->bmi[0].as_mode.first,
                    c->fc.i8x8_mode_prob);
    // printf("    mode: %d\n", m->bmi[0].as_mode.first); fflush(stdout);
    write_i8x8_mode(bc, m->bmi[2].as_mode.first,
                    c->fc.i8x8_mode_prob);
    // printf("    mode: %d\n", m->bmi[2].as_mode.first); fflush(stdout);
    write_i8x8_mode(bc, m->bmi[8].as_mode.first,
                    c->fc.i8x8_mode_prob);
    // printf("    mode: %d\n", m->bmi[8].as_mode.first); fflush(stdout);
    write_i8x8_mode(bc, m->bmi[10].as_mode.first,
                    c->fc.i8x8_mode_prob);
    // printf("    mode: %d\n", m->bmi[10].as_mode.first); fflush(stdout);
  } else
    write_uv_mode(bc, m->mbmi.uv_mode, c->kf_uv_mode_prob[ym]);

#if CONFIG_TX_SELECT
  if (ym <= I8X8_PRED && c->txfm_mode == TX_MODE_SELECT &&
      !((c->mb_no_coeff_skip && m->mbmi.mb_skip_coeff) ||
        (segfeature_active(xd, segment_id, SEG_LVL_EOB) &&
         get_segdata(xd, segment_id, SEG_LVL_EOB) == 0))) {
    TX_SIZE sz = m->mbmi.txfm_size;
    // FIXME(rbultje) code ternary symbol once all experiments are merged
    vp8_write(bc, sz != TX_4X4, c->prob_tx[0]);
    if (sz != TX_4X4 && ym <= TM_PRED)
      vp8_write(bc, sz != TX_8X8, c->prob_tx[1]);
  }
#endif
}

static void write_kfmodes(VP8_COMP* const cpi, vp8_writer* const bc) {
  VP8_COMMON *const c = &cpi->common;
  const int mis = c->mode_info_stride;
  MACROBLOCKD *xd = &cpi->mb.e_mbd;
  MODE_INFO *m;
  int i;
  int row, col;
  int mb_row, mb_col;
  int row_delta[4] = { 0, +1,  0, -1};
  int col_delta[4] = { +1, -1, +1, +1};
  TOKENEXTRA *tok = cpi->tok;
  TOKENEXTRA *tok_end = tok + cpi->tok_count;

  mb_row = 0;
  for (row = 0; row < c->mb_rows; row += 2) {
    m = c->mi + row * mis;

    mb_col = 0;
    for (col = 0; col < c->mb_cols; col += 2) {
#if CONFIG_SUPERBLOCKS
      vp8_write(bc, m->mbmi.encoded_as_sb, c->sb_coded);
#endif
      // Process the 4 MBs in the order:
      // top-left, top-right, bottom-left, bottom-right
      for (i = 0; i < 4; i++) {
        int dy = row_delta[i];
        int dx = col_delta[i];
        int offset_extended = dy * mis + dx;

        if ((mb_row >= c->mb_rows) || (mb_col >= c->mb_cols)) {
          // MB lies outside frame, move on
          mb_row += dy;
          mb_col += dx;
          m += offset_extended;
          continue;
        }

        // Make sure the MacroBlockD mode info pointer is set correctly
        xd->mode_info_context = m;

        write_mb_modes_kf(c, xd, m, mis, bc);
#ifdef ENTROPY_STATS
        active_section = 8;
#endif
        assert(tok < tok_end);
        pack_mb_tokens(bc, &tok, tok_end);

#if CONFIG_SUPERBLOCKS
        if (m->mbmi.encoded_as_sb) {
          assert(!i);
          mb_col += 2;
          m += 2;
          break;
        }
#endif
        // Next MB
        mb_row += dy;
        mb_col += dx;
        m += offset_extended;
      }
    }
    mb_row += 2;
  }
}


/* This function is used for debugging probability trees. */
static void print_prob_tree(vp8_prob
                            coef_probs[BLOCK_TYPES][COEF_BANDS][PREV_COEF_CONTEXTS][ENTROPY_NODES]) {
  /* print coef probability tree */
  int i, j, k, l;
  FILE *f = fopen("enc_tree_probs.txt", "a");
  fprintf(f, "{\n");
  for (i = 0; i < BLOCK_TYPES; i++) {
    fprintf(f, "  {\n");
    for (j = 0; j < COEF_BANDS; j++) {
      fprintf(f, "    {\n");
      for (k = 0; k < PREV_COEF_CONTEXTS; k++) {
        fprintf(f, "      {");
        for (l = 0; l < ENTROPY_NODES; l++) {
          fprintf(f, "%3u, ",
                  (unsigned int)(coef_probs [i][j][k][l]));
        }
        fprintf(f, " }\n");
      }
      fprintf(f, "    }\n");
    }
    fprintf(f, "  }\n");
  }
  fprintf(f, "}\n");
  fclose(f);
}


void build_coeff_contexts(VP8_COMP *cpi) {
  int i = 0, j, k;
#ifdef ENTROPY_STATS
  int t = 0;
#endif
  for (i = 0; i < BLOCK_TYPES; ++i) {
    for (j = 0; j < COEF_BANDS; ++j) {
      for (k = 0; k < PREV_COEF_CONTEXTS; ++k) {
        if (k >= 3 && ((i == 0 && j == 1) || (i > 0 && j == 0)))
          continue;
        vp8_tree_probs_from_distribution(
          MAX_ENTROPY_TOKENS, vp8_coef_encodings, vp8_coef_tree,
          cpi->frame_coef_probs [i][j][k],
          cpi->frame_branch_ct [i][j][k],
          cpi->coef_counts [i][j][k],
          256, 1
        );
#ifdef ENTROPY_STATS
        if (!cpi->dummy_packing)
          for (t = 0; t < MAX_ENTROPY_TOKENS; ++t)
            context_counters[i][j][k][t] += cpi->coef_counts[i][j][k][t];
#endif
      }
    }
  }
  for (i = 0; i < BLOCK_TYPES; ++i) {
    for (j = 0; j < COEF_BANDS; ++j) {
      for (k = 0; k < PREV_COEF_CONTEXTS; ++k) {
        if (k >= 3 && ((i == 0 && j == 1) || (i > 0 && j == 0)))
          continue;
        vp8_tree_probs_from_distribution(
          MAX_ENTROPY_TOKENS, vp8_coef_encodings, vp8_coef_tree,
          cpi->frame_hybrid_coef_probs [i][j][k],
          cpi->frame_hybrid_branch_ct [i][j][k],
          cpi->hybrid_coef_counts [i][j][k],
          256, 1
        );
#ifdef ENTROPY_STATS
        if (!cpi->dummy_packing)
          for (t = 0; t < MAX_ENTROPY_TOKENS; ++t)
            hybrid_context_counters[i][j][k][t] += cpi->hybrid_coef_counts[i][j][k][t];
#endif
      }
    }
  }

  if (cpi->common.txfm_mode != ONLY_4X4) {
    for (i = 0; i < BLOCK_TYPES_8X8; ++i) {
      for (j = 0; j < COEF_BANDS; ++j) {
        for (k = 0; k < PREV_COEF_CONTEXTS; ++k) {
          /* at every context */
          /* calc probs and branch cts for this frame only */
          // vp8_prob new_p           [ENTROPY_NODES];
          // unsigned int branch_ct   [ENTROPY_NODES] [2];
          if (k >= 3 && ((i == 0 && j == 1) || (i > 0 && j == 0)))
            continue;
          vp8_tree_probs_from_distribution(
            MAX_ENTROPY_TOKENS, vp8_coef_encodings, vp8_coef_tree,
            cpi->frame_coef_probs_8x8 [i][j][k],
            cpi->frame_branch_ct_8x8 [i][j][k],
            cpi->coef_counts_8x8 [i][j][k],
            256, 1
          );
#ifdef ENTROPY_STATS
          if (!cpi->dummy_packing)
            for (t = 0; t < MAX_ENTROPY_TOKENS; ++t)
              context_counters_8x8[i][j][k][t] += cpi->coef_counts_8x8[i][j][k][t];
#endif
        }
      }
    }
    for (i = 0; i < BLOCK_TYPES_8X8; ++i) {
      for (j = 0; j < COEF_BANDS; ++j) {
        for (k = 0; k < PREV_COEF_CONTEXTS; ++k) {
          /* at every context */
          /* calc probs and branch cts for this frame only */
          // vp8_prob new_p           [ENTROPY_NODES];
          // unsigned int branch_ct   [ENTROPY_NODES] [2];
          if (k >= 3 && ((i == 0 && j == 1) || (i > 0 && j == 0)))
            continue;
          vp8_tree_probs_from_distribution(
            MAX_ENTROPY_TOKENS, vp8_coef_encodings, vp8_coef_tree,
            cpi->frame_hybrid_coef_probs_8x8 [i][j][k],
            cpi->frame_hybrid_branch_ct_8x8 [i][j][k],
            cpi->hybrid_coef_counts_8x8 [i][j][k],
            256, 1
          );
#ifdef ENTROPY_STATS
          if (!cpi->dummy_packing)
            for (t = 0; t < MAX_ENTROPY_TOKENS; ++t)
              hybrid_context_counters_8x8[i][j][k][t] += cpi->hybrid_coef_counts_8x8[i][j][k][t];
#endif
        }
      }
    }
  }

  if (cpi->common.txfm_mode > ALLOW_8X8) {
    for (i = 0; i < BLOCK_TYPES_16X16; ++i) {
      for (j = 0; j < COEF_BANDS; ++j) {
        for (k = 0; k < PREV_COEF_CONTEXTS; ++k) {
          if (k >= 3 && ((i == 0 && j == 1) || (i > 0 && j == 0)))
            continue;
          vp8_tree_probs_from_distribution(
            MAX_ENTROPY_TOKENS, vp8_coef_encodings, vp8_coef_tree,
            cpi->frame_coef_probs_16x16[i][j][k],
            cpi->frame_branch_ct_16x16[i][j][k],
            cpi->coef_counts_16x16[i][j][k], 256, 1);
#ifdef ENTROPY_STATS
          if (!cpi->dummy_packing)
            for (t = 0; t < MAX_ENTROPY_TOKENS; ++t)
              context_counters_16x16[i][j][k][t] += cpi->coef_counts_16x16[i][j][k][t];
#endif
        }
      }
    }
  }
  for (i = 0; i < BLOCK_TYPES_16X16; ++i) {
    for (j = 0; j < COEF_BANDS; ++j) {
      for (k = 0; k < PREV_COEF_CONTEXTS; ++k) {
        if (k >= 3 && ((i == 0 && j == 1) || (i > 0 && j == 0)))
          continue;
        vp8_tree_probs_from_distribution(
          MAX_ENTROPY_TOKENS, vp8_coef_encodings, vp8_coef_tree,
          cpi->frame_hybrid_coef_probs_16x16[i][j][k],
          cpi->frame_hybrid_branch_ct_16x16[i][j][k],
          cpi->hybrid_coef_counts_16x16[i][j][k], 256, 1);
#ifdef ENTROPY_STATS
        if (!cpi->dummy_packing)
          for (t = 0; t < MAX_ENTROPY_TOKENS; ++t)
            hybrid_context_counters_16x16[i][j][k][t] += cpi->hybrid_coef_counts_16x16[i][j][k][t];
#endif
      }
    }
  }
}

#if 0
static void update_coef_probs2(VP8_COMP *cpi) {
  const vp8_prob grpupd = 192;
  int i, j, k, t;
  vp8_writer *const w = &cpi->bc;
  int update[2];
  int savings;

  vp8_clear_system_state(); // __asm emms;
  // Build the cofficient contexts based on counts collected in encode loop
  build_coeff_contexts(cpi);

  for (t = 0; t < ENTROPY_NODES; ++t) {
    /* dry run to see if there is any udpate at all needed */
    savings = 0;
    update[0] = update[1] = 0;
    for (i = 0; i < BLOCK_TYPES; ++i) {
      for (j = !i; j < COEF_BANDS; ++j) {
        for (k = 0; k < PREV_COEF_CONTEXTS; ++k) {
          vp8_prob newp = cpi->frame_coef_probs [i][j][k][t];
          vp8_prob *Pold = cpi->common.fc.coef_probs [i][j][k] + t;
          const vp8_prob upd = COEF_UPDATE_PROB;
          int s;
          int u = 0;
          if (k >= 3 && ((i == 0 && j == 1) || (i > 0 && j == 0)))
            continue;

#if defined(SEARCH_NEWP)
          s = prob_diff_update_savings_search(
                cpi->frame_branch_ct [i][j][k][t], *Pold, &newp, upd);
          if (s > 0 && newp != *Pold) u = 1;
          if (u)
            savings += s - (int)(vp8_cost_zero(upd));
          else
            savings -= (int)(vp8_cost_zero(upd));
#else
          s = prob_update_savings(
                cpi->frame_branch_ct [i][j][k][t], *Pold, newp, upd);
          if (s > 0) u = 1;
          if (u)
            savings += s;
#endif
          // printf("    %d %d %d: %d\n", i, j, k, u);
          update[u]++;
        }
      }
    }
    if (update[1] == 0 || savings < 0) {
      vp8_write(w, 0, grpupd);
      continue;
    }
    vp8_write(w, 1, grpupd);
    for (i = 0; i < BLOCK_TYPES; ++i) {
      for (j = !i; j < COEF_BANDS; ++j) {
        for (k = 0; k < PREV_COEF_CONTEXTS; ++k) {
          vp8_prob newp = cpi->frame_coef_probs [i][j][k][t];
          vp8_prob *Pold = cpi->common.fc.coef_probs [i][j][k] + t;
          const vp8_prob upd = COEF_UPDATE_PROB;
          int s;
          int u = 0;
          if (k >= 3 && ((i == 0 && j == 1) || (i > 0 && j == 0)))
            continue;
#if defined(SEARCH_NEWP)
          s = prob_diff_update_savings_search(
                cpi->frame_branch_ct [i][j][k][t], *Pold, &newp, upd);
          if (s > 0 && newp != *Pold) u = 1;
#else
          s = prob_update_savings(
                cpi->frame_branch_ct [i][j][k][t], *Pold, newp, upd);
          if (s > 0) u = 1;
#endif
          // printf("  %d %d %d: %d (%d)\n", i, j, k, u, upd);
          vp8_write(w, u, upd);
#ifdef ENTROPY_STATS
          ++ tree_update_hist [i][j][k][t] [u];
#endif
          if (u) {
            /* send/use new probability */
            write_prob_diff_update(w, newp, *Pold);
            *Pold = newp;
          }
        }
      }
    }
  }

  if (cpi->common.txfm_mode != ONLY_4X4)
  for (t = 0; t < ENTROPY_NODES; ++t) {
    /* dry run to see if there is any udpate at all needed */
    savings = 0;
    update[0] = update[1] = 0;
    for (i = 0; i < BLOCK_TYPES_8X8; ++i) {
      for (j = !i; j < COEF_BANDS; ++j) {
        for (k = 0; k < PREV_COEF_CONTEXTS; ++k) {
          vp8_prob newp = cpi->frame_coef_probs_8x8 [i][j][k][t];
          vp8_prob *Pold = cpi->common.fc.coef_probs_8x8 [i][j][k] + t;
          const vp8_prob upd = COEF_UPDATE_PROB_8X8;
          int s;
          int u = 0;
          if (k >= 3 && ((i == 0 && j == 1) || (i > 0 && j == 0)))
            continue;
#if defined(SEARCH_NEWP)
          s = prob_diff_update_savings_search(
                cpi->frame_branch_ct_8x8 [i][j][k][t],
                *Pold, &newp, upd);
          if (s > 0 && newp != *Pold)
            u = 1;
          if (u)
            savings += s - (int)(vp8_cost_zero(upd));
          else
            savings -= (int)(vp8_cost_zero(upd));
#else
          s = prob_update_savings(
                cpi->frame_branch_ct_8x8 [i][j][k][t],
                *Pold, newp, upd);
          if (s > 0)
            u = 1;
          if (u)
            savings += s;
#endif
          update[u]++;
        }
      }
    }
    if (update[1] == 0 || savings < 0) {
      vp8_write(w, 0, grpupd);
      continue;
    }
    vp8_write(w, 1, grpupd);
    for (i = 0; i < BLOCK_TYPES_8X8; ++i) {
      for (j = !i; j < COEF_BANDS; ++j) {
        for (k = 0; k < PREV_COEF_CONTEXTS; ++k) {
          vp8_prob newp = cpi->frame_coef_probs_8x8 [i][j][k][t];
          vp8_prob *Pold = cpi->common.fc.coef_probs_8x8 [i][j][k] + t;
          const vp8_prob upd = COEF_UPDATE_PROB_8X8;
          int s;
          int u = 0;
          if (k >= 3 && ((i == 0 && j == 1) || (i > 0 && j == 0)))
            continue;
#if defined(SEARCH_NEWP)
          s = prob_diff_update_savings_search(
                cpi->frame_branch_ct_8x8 [i][j][k][t],
                *Pold, &newp, upd);
          if (s > 0 && newp != *Pold)
            u = 1;
#else
          s = prob_update_savings(
                cpi->frame_branch_ct_8x8 [i][j][k][t],
                *Pold, newp, upd);
          if (s > 0)
            u = 1;
#endif
          vp8_write(w, u, upd);
#ifdef ENTROPY_STATS
          if (!cpi->dummy_packing)
            ++ tree_update_hist_8x8 [i][j][k][t] [u];
#endif
          if (u) {
            /* send/use new probability */
            write_prob_diff_update(w, newp, *Pold);
            *Pold = newp;
          }
        }
      }
    }
  }
}
#endif

static void update_coef_probs(VP8_COMP* const cpi, vp8_writer* const bc) {
  int i, j, k, t;
  int update[2] = {0, 0};
  int savings;

  vp8_clear_system_state(); // __asm emms;

  // Build the cofficient contexts based on counts collected in encode loop
  build_coeff_contexts(cpi);

  // vp8_prob bestupd = find_coef_update_prob(cpi);

  /* dry run to see if there is any udpate at all needed */
  savings = 0;
  for (i = 0; i < BLOCK_TYPES; ++i) {
    for (j = !i; j < COEF_BANDS; ++j) {
      int prev_coef_savings[ENTROPY_NODES] = {0};
      for (k = 0; k < PREV_COEF_CONTEXTS; ++k) {
        for (t = 0; t < ENTROPY_NODES; ++t) {
          vp8_prob newp = cpi->frame_coef_probs [i][j][k][t];
          vp8_prob *Pold = cpi->common.fc.coef_probs [i][j][k] + t;
          const vp8_prob upd = COEF_UPDATE_PROB;
          int s = prev_coef_savings[t];
          int u = 0;
          if (k >= 3 && ((i == 0 && j == 1) || (i > 0 && j == 0)))
            continue;
#if defined(SEARCH_NEWP)
          s = prob_diff_update_savings_search(
                cpi->frame_branch_ct [i][j][k][t],
                *Pold, &newp, upd);
          if (s > 0 && newp != *Pold)
            u = 1;
          if (u)
            savings += s - (int)(vp8_cost_zero(upd));
          else
            savings -= (int)(vp8_cost_zero(upd));
#else
          s = prob_update_savings(
                cpi->frame_branch_ct [i][j][k][t],
                *Pold, newp, upd);
          if (s > 0)
            u = 1;
          if (u)
            savings += s;
#endif

          update[u]++;
        }
      }
    }
  }

  // printf("Update %d %d, savings %d\n", update[0], update[1], savings);
  /* Is coef updated at all */
  if (update[1] == 0 || savings < 0) {
    vp8_write_bit(bc, 0);
  } else {
    vp8_write_bit(bc, 1);
    for (i = 0; i < BLOCK_TYPES; ++i) {
      for (j = !i; j < COEF_BANDS; ++j) {
        int prev_coef_savings[ENTROPY_NODES] = {0};
        for (k = 0; k < PREV_COEF_CONTEXTS; ++k) {
          // calc probs and branch cts for this frame only
          for (t = 0; t < ENTROPY_NODES; ++t) {
            vp8_prob newp = cpi->frame_coef_probs [i][j][k][t];
            vp8_prob *Pold = cpi->common.fc.coef_probs [i][j][k] + t;
            const vp8_prob upd = COEF_UPDATE_PROB;
            int s = prev_coef_savings[t];
            int u = 0;
            if (k >= 3 && ((i == 0 && j == 1) || (i > 0 && j == 0)))
              continue;

#if defined(SEARCH_NEWP)
            s = prob_diff_update_savings_search(
                  cpi->frame_branch_ct [i][j][k][t],
                  *Pold, &newp, upd);
            if (s > 0 && newp != *Pold)
              u = 1;
#else
            s = prob_update_savings(
                  cpi->frame_branch_ct [i][j][k][t],
                  *Pold, newp, upd);
            if (s > 0)
              u = 1;
#endif
            vp8_write(bc, u, upd);
#ifdef ENTROPY_STATS
            if (!cpi->dummy_packing)
              ++ tree_update_hist [i][j][k][t] [u];
#endif
            if (u) {
              /* send/use new probability */
              write_prob_diff_update(bc, newp, *Pold);
              *Pold = newp;
            }
          }
        }
      }
    }
  }

  savings = 0;
  update[0] = update[1] = 0;
  for (i = 0; i < BLOCK_TYPES; ++i) {
    for (j = !i; j < COEF_BANDS; ++j) {
      int prev_coef_savings[ENTROPY_NODES] = {0};
      for (k = 0; k < PREV_COEF_CONTEXTS; ++k) {
        for (t = 0; t < ENTROPY_NODES; ++t) {
          vp8_prob newp = cpi->frame_hybrid_coef_probs [i][j][k][t];
          vp8_prob *Pold = cpi->common.fc.hybrid_coef_probs [i][j][k] + t;
          const vp8_prob upd = COEF_UPDATE_PROB;
          int s = prev_coef_savings[t];
          int u = 0;
          if (k >= 3 && ((i == 0 && j == 1) || (i > 0 && j == 0)))
            continue;
#if defined(SEARCH_NEWP)
          s = prob_diff_update_savings_search(
                cpi->frame_hybrid_branch_ct [i][j][k][t],
                *Pold, &newp, upd);
          if (s > 0 && newp != *Pold)
            u = 1;
          if (u)
            savings += s - (int)(vp8_cost_zero(upd));
          else
            savings -= (int)(vp8_cost_zero(upd));
#else
          s = prob_update_savings(
                cpi->frame_hybrid_branch_ct [i][j][k][t],
                *Pold, newp, upd);
          if (s > 0)
            u = 1;
          if (u)
            savings += s;
#endif

          update[u]++;
        }
      }
    }
  }

  // printf("Update %d %d, savings %d\n", update[0], update[1], savings);
  /* Is coef updated at all */
  if (update[1] == 0 || savings < 0) {
    vp8_write_bit(bc, 0);
  } else {
    vp8_write_bit(bc, 1);
    for (i = 0; i < BLOCK_TYPES; ++i) {
      for (j = !i; j < COEF_BANDS; ++j) {
        int prev_coef_savings[ENTROPY_NODES] = {0};
        for (k = 0; k < PREV_COEF_CONTEXTS; ++k) {
          // calc probs and branch cts for this frame only
          for (t = 0; t < ENTROPY_NODES; ++t) {
            vp8_prob newp = cpi->frame_hybrid_coef_probs [i][j][k][t];
            vp8_prob *Pold = cpi->common.fc.hybrid_coef_probs [i][j][k] + t;
            const vp8_prob upd = COEF_UPDATE_PROB;
            int s = prev_coef_savings[t];
            int u = 0;
            if (k >= 3 && ((i == 0 && j == 1) || (i > 0 && j == 0)))
              continue;

#if defined(SEARCH_NEWP)
            s = prob_diff_update_savings_search(
                  cpi->frame_hybrid_branch_ct [i][j][k][t],
                  *Pold, &newp, upd);
            if (s > 0 && newp != *Pold)
              u = 1;
#else
            s = prob_update_savings(
                  cpi->frame_hybrid_branch_ct [i][j][k][t],
                  *Pold, newp, upd);
            if (s > 0)
              u = 1;
#endif
            vp8_write(bc, u, upd);
#ifdef ENTROPY_STATS
            if (!cpi->dummy_packing)
              ++ hybrid_tree_update_hist [i][j][k][t] [u];
#endif
            if (u) {
              /* send/use new probability */
              write_prob_diff_update(bc, newp, *Pold);
              *Pold = newp;
            }
          }
        }
      }
    }
  }

  /* do not do this if not even allowed */
  if (cpi->common.txfm_mode != ONLY_4X4) {
    /* dry run to see if update is necessary */
    update[0] = update[1] = 0;
    savings = 0;
    for (i = 0; i < BLOCK_TYPES_8X8; ++i) {
      for (j = !i; j < COEF_BANDS; ++j) {
        for (k = 0; k < PREV_COEF_CONTEXTS; ++k) {
          // calc probs and branch cts for this frame only
          for (t = 0; t < ENTROPY_NODES; ++t) {
            const unsigned int *ct  = cpi->frame_branch_ct_8x8 [i][j][k][t];
            vp8_prob newp = cpi->frame_coef_probs_8x8 [i][j][k][t];
            vp8_prob *Pold = cpi->common.fc.coef_probs_8x8 [i][j][k] + t;
            const vp8_prob oldp = *Pold;
            int s, u;
            const vp8_prob upd = COEF_UPDATE_PROB_8X8;
            if (k >= 3 && ((i == 0 && j == 1) || (i > 0 && j == 0)))
              continue;
#if defined(SEARCH_NEWP)
            s = prob_diff_update_savings_search(ct, oldp, &newp, upd);
            u = s > 0 && newp != oldp ? 1 : 0;
            if (u)
              savings += s - (int)(vp8_cost_zero(upd));
            else
              savings -= (int)(vp8_cost_zero(upd));
#else
            s = prob_update_savings(ct, oldp, newp, upd);
            u = s > 0 ? 1 : 0;
            if (u)
              savings += s;
#endif
            update[u]++;
          }
        }
      }
    }

    if (update[1] == 0 || savings < 0) {
      vp8_write_bit(bc, 0);
    } else {
      vp8_write_bit(bc, 1);
      for (i = 0; i < BLOCK_TYPES_8X8; ++i) {
        for (j = !i; j < COEF_BANDS; ++j) {
          for (k = 0; k < PREV_COEF_CONTEXTS; ++k) {
            for (t = 0; t < ENTROPY_NODES; ++t) {
              const unsigned int *ct  = cpi->frame_branch_ct_8x8 [i][j][k][t];
              vp8_prob newp = cpi->frame_coef_probs_8x8 [i][j][k][t];
              vp8_prob *Pold = cpi->common.fc.coef_probs_8x8 [i][j][k] + t;
              const vp8_prob oldp = *Pold;
              const vp8_prob upd = COEF_UPDATE_PROB_8X8;
              int s, u;
              if (k >= 3 && ((i == 0 && j == 1) ||
                             (i > 0 && j == 0)))
                continue;
#if defined(SEARCH_NEWP)
              s = prob_diff_update_savings_search(ct, oldp, &newp, upd);
              u = s > 0 && newp != oldp ? 1 : 0;
#else
              s = prob_update_savings(ct, oldp, newp, upd);
              u = s > 0 ? 1 : 0;
#endif
              vp8_write(bc, u, upd);
#ifdef ENTROPY_STATS
              if (!cpi->dummy_packing)
                ++ tree_update_hist_8x8 [i][j][k][t] [u];
#endif
              if (u) {
                /* send/use new probability */
                write_prob_diff_update(bc, newp, oldp);
                *Pold = newp;
              }
            }
          }
        }
      }
    }
    update[0] = update[1] = 0;
    savings = 0;
    for (i = 0; i < BLOCK_TYPES_8X8; ++i) {
      for (j = !i; j < COEF_BANDS; ++j) {
        for (k = 0; k < PREV_COEF_CONTEXTS; ++k) {
          // calc probs and branch cts for this frame only
          for (t = 0; t < ENTROPY_NODES; ++t) {
            const unsigned int *ct  = cpi->frame_hybrid_branch_ct_8x8 [i][j][k][t];
            vp8_prob newp = cpi->frame_hybrid_coef_probs_8x8 [i][j][k][t];
            vp8_prob *Pold = cpi->common.fc.hybrid_coef_probs_8x8 [i][j][k] + t;
            const vp8_prob oldp = *Pold;
            int s, u;
            const vp8_prob upd = COEF_UPDATE_PROB_8X8;
            if (k >= 3 && ((i == 0 && j == 1) || (i > 0 && j == 0)))
              continue;
#if defined(SEARCH_NEWP)
            s = prob_diff_update_savings_search(ct, oldp, &newp, upd);
            u = s > 0 && newp != oldp ? 1 : 0;
            if (u)
              savings += s - (int)(vp8_cost_zero(upd));
            else
              savings -= (int)(vp8_cost_zero(upd));
#else
            s = prob_update_savings(ct, oldp, newp, upd);
            u = s > 0 ? 1 : 0;
            if (u)
              savings += s;
#endif
            update[u]++;
          }
        }
      }
    }

    if (update[1] == 0 || savings < 0) {
      vp8_write_bit(bc, 0);
    } else {
      vp8_write_bit(bc, 1);
      for (i = 0; i < BLOCK_TYPES_8X8; ++i) {
        for (j = !i; j < COEF_BANDS; ++j) {
          for (k = 0; k < PREV_COEF_CONTEXTS; ++k) {
            for (t = 0; t < ENTROPY_NODES; ++t) {
              const unsigned int *ct  = cpi->frame_hybrid_branch_ct_8x8 [i][j][k][t];
              vp8_prob newp = cpi->frame_hybrid_coef_probs_8x8 [i][j][k][t];
              vp8_prob *Pold = cpi->common.fc.hybrid_coef_probs_8x8 [i][j][k] + t;
              const vp8_prob oldp = *Pold;
              const vp8_prob upd = COEF_UPDATE_PROB_8X8;
              int s, u;
              if (k >= 3 && ((i == 0 && j == 1) ||
                             (i > 0 && j == 0)))
                continue;
#if defined(SEARCH_NEWP)
              s = prob_diff_update_savings_search(ct, oldp, &newp, upd);
              u = s > 0 && newp != oldp ? 1 : 0;
#else
              s = prob_update_savings(ct, oldp, newp, upd);
              u = s > 0 ? 1 : 0;
#endif
              vp8_write(bc, u, upd);
#ifdef ENTROPY_STATS
              if (!cpi->dummy_packing)
                ++ hybrid_tree_update_hist_8x8 [i][j][k][t] [u];
#endif
              if (u) {
                /* send/use new probability */
                write_prob_diff_update(bc, newp, oldp);
                *Pold = newp;
              }
            }
          }
        }
      }
    }
  }

  if (cpi->common.txfm_mode > ALLOW_8X8) {
  /* dry run to see if update is necessary */
  update[0] = update[1] = 0;
  savings = 0;
  for (i = 0; i < BLOCK_TYPES_16X16; ++i) {
    for (j = !i; j < COEF_BANDS; ++j) {
      for (k = 0; k < PREV_COEF_CONTEXTS; ++k) {
        // calc probs and branch cts for this frame only
        for (t = 0; t < ENTROPY_NODES; ++t) {
          const unsigned int *ct  = cpi->frame_branch_ct_16x16[i][j][k][t];
          vp8_prob newp = cpi->frame_coef_probs_16x16[i][j][k][t];
          vp8_prob *Pold = cpi->common.fc.coef_probs_16x16[i][j][k] + t;
          const vp8_prob oldp = *Pold;
          int s, u;
          const vp8_prob upd = COEF_UPDATE_PROB_16X16;
          if (k >= 3 && ((i == 0 && j == 1) || (i > 0 && j == 0)))
            continue;
#if defined(SEARCH_NEWP)
          s = prob_diff_update_savings_search(ct, oldp, &newp, upd);
          u = s > 0 && newp != oldp ? 1 : 0;
          if (u)
            savings += s - (int)(vp8_cost_zero(upd));
          else
            savings -= (int)(vp8_cost_zero(upd));
#else
          s = prob_update_savings(ct, oldp, newp, upd);
          u = s > 0 ? 1 : 0;
          if (u)
            savings += s;
#endif
          update[u]++;
        }
      }
    }
  }

  if (update[1] == 0 || savings < 0) {
    vp8_write_bit(bc, 0);
  } else {
    vp8_write_bit(bc, 1);
    for (i = 0; i < BLOCK_TYPES_16X16; ++i) {
      for (j = !i; j < COEF_BANDS; ++j) {
        for (k = 0; k < PREV_COEF_CONTEXTS; ++k) {
          for (t = 0; t < ENTROPY_NODES; ++t) {
            const unsigned int *ct  = cpi->frame_branch_ct_16x16[i][j][k][t];
            vp8_prob newp = cpi->frame_coef_probs_16x16[i][j][k][t];
            vp8_prob *Pold = cpi->common.fc.coef_probs_16x16[i][j][k] + t;
            const vp8_prob oldp = *Pold;
            const vp8_prob upd = COEF_UPDATE_PROB_16X16;
            int s, u;
            if (k >= 3 && ((i == 0 && j == 1) ||
                           (i > 0 && j == 0)))
              continue;
#if defined(SEARCH_NEWP)
            s = prob_diff_update_savings_search(ct, oldp, &newp, upd);
            u = s > 0 && newp != oldp ? 1 : 0;
#else
            s = prob_update_savings(ct, oldp, newp, upd);
            u = s > 0 ? 1 : 0;
#endif
            vp8_write(bc, u, upd);
#ifdef ENTROPY_STATS
            if (!cpi->dummy_packing)
              ++tree_update_hist_16x16[i][j][k][t][u];
#endif
            if (u) {
              /* send/use new probability */
              write_prob_diff_update(bc, newp, oldp);
              *Pold = newp;
            }
          }
        }
      }
    }
  }
  update[0] = update[1] = 0;
  savings = 0;
  for (i = 0; i < BLOCK_TYPES_16X16; ++i) {
    for (j = !i; j < COEF_BANDS; ++j) {
      for (k = 0; k < PREV_COEF_CONTEXTS; ++k) {
        // calc probs and branch cts for this frame only
        for (t = 0; t < ENTROPY_NODES; ++t) {
          const unsigned int *ct  = cpi->frame_hybrid_branch_ct_16x16[i][j][k][t];
          vp8_prob newp = cpi->frame_hybrid_coef_probs_16x16[i][j][k][t];
          vp8_prob *Pold = cpi->common.fc.hybrid_coef_probs_16x16[i][j][k] + t;
          const vp8_prob oldp = *Pold;
          int s, u;
          const vp8_prob upd = COEF_UPDATE_PROB_16X16;
          if (k >= 3 && ((i == 0 && j == 1) || (i > 0 && j == 0)))
            continue;
#if defined(SEARCH_NEWP)
          s = prob_diff_update_savings_search(ct, oldp, &newp, upd);
          u = s > 0 && newp != oldp ? 1 : 0;
          if (u)
            savings += s - (int)(vp8_cost_zero(upd));
          else
            savings -= (int)(vp8_cost_zero(upd));
#else
          s = prob_update_savings(ct, oldp, newp, upd);
          u = s > 0 ? 1 : 0;
          if (u)
            savings += s;
#endif
          update[u]++;
        }
      }
    }
  }

  if (update[1] == 0 || savings < 0) {
    vp8_write_bit(bc, 0);
  } else {
    vp8_write_bit(bc, 1);
    for (i = 0; i < BLOCK_TYPES_16X16; ++i) {
      for (j = !i; j < COEF_BANDS; ++j) {
        for (k = 0; k < PREV_COEF_CONTEXTS; ++k) {
          for (t = 0; t < ENTROPY_NODES; ++t) {
            const unsigned int *ct  = cpi->frame_hybrid_branch_ct_16x16[i][j][k][t];
            vp8_prob newp = cpi->frame_hybrid_coef_probs_16x16[i][j][k][t];
            vp8_prob *Pold = cpi->common.fc.hybrid_coef_probs_16x16[i][j][k] + t;
            const vp8_prob oldp = *Pold;
            const vp8_prob upd = COEF_UPDATE_PROB_16X16;
            int s, u;
            if (k >= 3 && ((i == 0 && j == 1) ||
                           (i > 0 && j == 0)))
              continue;
#if defined(SEARCH_NEWP)
            s = prob_diff_update_savings_search(ct, oldp, &newp, upd);
            u = s > 0 && newp != oldp ? 1 : 0;
#else
            s = prob_update_savings(ct, oldp, newp, upd);
            u = s > 0 ? 1 : 0;
#endif
            vp8_write(bc, u, upd);
#ifdef ENTROPY_STATS
            if (!cpi->dummy_packing)
              ++hybrid_tree_update_hist_16x16[i][j][k][t][u];
#endif
            if (u) {
              /* send/use new probability */
              write_prob_diff_update(bc, newp, oldp);
              *Pold = newp;
            }
          }
        }
      }
    }
  }
  }
}

#ifdef PACKET_TESTING
FILE *vpxlogc = 0;
#endif

static void put_delta_q(vp8_writer *bc, int delta_q) {
  if (delta_q != 0) {
    vp8_write_bit(bc, 1);
    vp8_write_literal(bc, abs(delta_q), 4);

    if (delta_q < 0)
      vp8_write_bit(bc, 1);
    else
      vp8_write_bit(bc, 0);
  } else
    vp8_write_bit(bc, 0);
}

static void decide_kf_ymode_entropy(VP8_COMP *cpi) {

  int mode_cost[MB_MODE_COUNT];
  int cost;
  int bestcost = INT_MAX;
  int bestindex = 0;
  int i, j;

  for (i = 0; i < 8; i++) {
    vp8_cost_tokens(mode_cost, cpi->common.kf_ymode_prob[i], vp8_kf_ymode_tree);
    cost = 0;
    for (j = 0; j < VP8_YMODES; j++) {
      cost += mode_cost[j] * cpi->ymode_count[j];
    }
#if CONFIG_SUPERBLOCKS
    vp8_cost_tokens(mode_cost, cpi->common.sb_kf_ymode_prob[i],
                    vp8_sb_ymode_tree);
    for (j = 0; j < VP8_I32X32_MODES; j++) {
      cost += mode_cost[j] * cpi->sb_ymode_count[j];
    }
#endif
    if (cost < bestcost) {
      bestindex = i;
      bestcost = cost;
    }
  }
  cpi->common.kf_ymode_probs_index = bestindex;

}
static void segment_reference_frames(VP8_COMP *cpi) {
  VP8_COMMON *oci = &cpi->common;
  MODE_INFO *mi = oci->mi;
  int ref[MAX_MB_SEGMENTS] = {0};
  int i, j;
  int mb_index = 0;
  MACROBLOCKD *const xd = &cpi->mb.e_mbd;

  for (i = 0; i < oci->mb_rows; i++) {
    for (j = 0; j < oci->mb_cols; j++, mb_index++) {
      ref[mi[mb_index].mbmi.segment_id] |= (1 << mi[mb_index].mbmi.ref_frame);
    }
    mb_index++;
  }
  for (i = 0; i < MAX_MB_SEGMENTS; i++) {
    enable_segfeature(xd, i, SEG_LVL_REF_FRAME);
    set_segdata(xd, i, SEG_LVL_REF_FRAME, ref[i]);
  }
}

void vp8_pack_bitstream(VP8_COMP *cpi, unsigned char *dest, unsigned long *size) {
  int i, j;
  VP8_HEADER oh;
  VP8_COMMON *const pc = &cpi->common;
  vp8_writer header_bc, residual_bc;
  MACROBLOCKD *const xd = &cpi->mb.e_mbd;
  int extra_bytes_packed = 0;

  unsigned char *cx_data = dest;

  oh.show_frame = (int) pc->show_frame;
  oh.type = (int)pc->frame_type;
  oh.version = pc->version;
  oh.first_partition_length_in_bytes = 0;

  cx_data += 3;

#if defined(SECTIONBITS_OUTPUT)
  Sectionbits[active_section = 1] += sizeof(VP8_HEADER) * 8 * 256;
#endif

  compute_update_table();

  // vp8_kf_default_bmode_probs() is called in vp8_setup_key_frame() once for each
  // K frame before encode frame. pc->kf_bmode_prob doesn't get changed anywhere
  // else. No need to call it again here. --yw
  // vp8_kf_default_bmode_probs( pc->kf_bmode_prob);

  // every keyframe send startcode, width, height, scale factor, clamp and color type
  if (oh.type == KEY_FRAME) {
    int v;

    // Start / synch code
    cx_data[0] = 0x9D;
    cx_data[1] = 0x01;
    cx_data[2] = 0x2a;

    v = (pc->horiz_scale << 14) | pc->Width;
    cx_data[3] = v;
    cx_data[4] = v >> 8;

    v = (pc->vert_scale << 14) | pc->Height;
    cx_data[5] = v;
    cx_data[6] = v >> 8;

    extra_bytes_packed = 7;
    cx_data += extra_bytes_packed;

    vp8_start_encode(&header_bc, cx_data);

    // signal clr type
    vp8_write_bit(&header_bc, pc->clr_type);
    vp8_write_bit(&header_bc, pc->clamp_type);

  } else {
    vp8_start_encode(&header_bc, cx_data);
  }

  // Signal whether or not Segmentation is enabled
  vp8_write_bit(&header_bc, (xd->segmentation_enabled) ? 1 : 0);

  // Indicate which features are enabled
  if (xd->segmentation_enabled) {
    // Indicate whether or not the segmentation map is being updated.
    vp8_write_bit(&header_bc, (xd->update_mb_segmentation_map) ? 1 : 0);

    // If it is, then indicate the method that will be used.
    if (xd->update_mb_segmentation_map) {
      // Select the coding strategy (temporal or spatial)
      choose_segmap_coding_method(cpi);
      // Send the tree probabilities used to decode unpredicted
      // macro-block segments
      for (i = 0; i < MB_FEATURE_TREE_PROBS; i++) {
        int data = xd->mb_segment_tree_probs[i];

        if (data != 255) {
          vp8_write_bit(&header_bc, 1);
          vp8_write_literal(&header_bc, data, 8);
        } else {
          vp8_write_bit(&header_bc, 0);
        }
      }

      // Write out the chosen coding method.
      vp8_write_bit(&header_bc, (pc->temporal_update) ? 1 : 0);
      if (pc->temporal_update) {
        for (i = 0; i < PREDICTION_PROBS; i++) {
          int data = pc->segment_pred_probs[i];

          if (data != 255) {
            vp8_write_bit(&header_bc, 1);
            vp8_write_literal(&header_bc, data, 8);
          } else {
            vp8_write_bit(&header_bc, 0);
          }
        }
      }
    }

    vp8_write_bit(&header_bc, (xd->update_mb_segmentation_data) ? 1 : 0);

    // segment_reference_frames(cpi);

    if (xd->update_mb_segmentation_data) {
      signed char Data;

      vp8_write_bit(&header_bc, (xd->mb_segment_abs_delta) ? 1 : 0);

      // For each segments id...
      for (i = 0; i < MAX_MB_SEGMENTS; i++) {
        // For each segmentation codable feature...
        for (j = 0; j < SEG_LVL_MAX; j++) {
          Data = get_segdata(xd, i, j);


#if CONFIG_FEATUREUPDATES

          // check if there's an update
          if (segfeature_changed(xd, i, j)) {
            vp8_write_bit(&header_bc, 1);

            if (segfeature_active(xd, i, j)) {
              // this bit is to say we are still
              // active/  if we were inactive
              // this is unnecessary
              if (old_segfeature_active(xd, i, j)) {
                vp8_write_bit(&header_bc, 1);
              }
              // Is the segment data signed..
              if (is_segfeature_signed(j)) {
                // Encode the relevant feature data
                if (Data < 0) {
                  Data = - Data;
                  vp8_write_literal(&header_bc, Data,
                                    seg_feature_data_bits(j));
                  vp8_write_bit(&header_bc, 1);
                } else {
                  vp8_write_literal(&header_bc, Data,
                                    seg_feature_data_bits(j));
                  vp8_write_bit(&header_bc, 0);
                }
              }
              // Unsigned data element so no sign bit needed
              else
                vp8_write_literal(&header_bc, Data,
                                  seg_feature_data_bits(j));
            }
            // feature is inactive now
            else if (old_segfeature_active(xd, i, j)) {
              vp8_write_bit(&header_bc, 0);
            }
          } else {
            vp8_write_bit(&header_bc, 0);
          }
#else

          // If the feature is enabled...
          if (segfeature_active(xd, i, j)) {
            vp8_write_bit(&header_bc, 1);

            // Is the segment data signed..
            if (is_segfeature_signed(j)) {
              // Encode the relevant feature data
              if (Data < 0) {
                Data = - Data;
                vp8_write_literal(&header_bc, Data,
                                  seg_feature_data_bits(j));
                vp8_write_bit(&header_bc, 1);
              } else {
                vp8_write_literal(&header_bc, Data,
                                  seg_feature_data_bits(j));
                vp8_write_bit(&header_bc, 0);
              }
            }
            // Unsigned data element so no sign bit needed
            else
              vp8_write_literal(&header_bc, Data,
                                seg_feature_data_bits(j));
          } else
            vp8_write_bit(&header_bc, 0);
#endif
        }
      }
    }

#if CONFIG_FEATUREUPDATES
    // save the segment info for updates next frame
    save_segment_info(xd);
#endif

  }

  // Encode the common prediction model status flag probability updates for
  // the reference frame
  update_refpred_stats(cpi);
  if (pc->frame_type != KEY_FRAME) {
    for (i = 0; i < PREDICTION_PROBS; i++) {
      if (cpi->ref_pred_probs_update[i]) {
        vp8_write_bit(&header_bc, 1);
        vp8_write_literal(&header_bc, pc->ref_pred_probs[i], 8);
      } else {
        vp8_write_bit(&header_bc, 0);
      }
    }
  }

#if CONFIG_SUPERBLOCKS
  {
    /* sb mode probability */
    const int sb_max = (((pc->mb_rows + 1) >> 1) * ((pc->mb_cols + 1) >> 1));

    pc->sb_coded = get_prob(sb_max - cpi->sb_count, sb_max);
    vp8_write_literal(&header_bc, pc->sb_coded, 8);
  }
#endif

#if CONFIG_TX_SELECT
  {
    if (pc->txfm_mode == TX_MODE_SELECT) {
      pc->prob_tx[0] = get_prob(cpi->txfm_count[0] + cpi->txfm_count_8x8p[0],
                                cpi->txfm_count[0] + cpi->txfm_count[1] + cpi->txfm_count[2] +
                                cpi->txfm_count_8x8p[0] + cpi->txfm_count_8x8p[1]);
      pc->prob_tx[1] = get_prob(cpi->txfm_count[1], cpi->txfm_count[1] + cpi->txfm_count[2]);
    } else {
      pc->prob_tx[0] = 128;
      pc->prob_tx[1] = 128;
    }
    vp8_write_literal(&header_bc, pc->txfm_mode, 2);
    if (pc->txfm_mode == TX_MODE_SELECT) {
      vp8_write_literal(&header_bc, pc->prob_tx[0], 8);
      vp8_write_literal(&header_bc, pc->prob_tx[1], 8);
    }
  }
#else
  vp8_write_bit(&header_bc, !!pc->txfm_mode);
#endif

  // Encode the loop filter level and type
  vp8_write_bit(&header_bc, pc->filter_type);
  vp8_write_literal(&header_bc, pc->filter_level, 6);
  vp8_write_literal(&header_bc, pc->sharpness_level, 3);

  // Write out loop filter deltas applied at the MB level based on mode or ref frame (if they are enabled).
  vp8_write_bit(&header_bc, (xd->mode_ref_lf_delta_enabled) ? 1 : 0);

  if (xd->mode_ref_lf_delta_enabled) {
    // Do the deltas need to be updated
    int send_update = xd->mode_ref_lf_delta_update;

    vp8_write_bit(&header_bc, send_update);
    if (send_update) {
      int Data;

      // Send update
      for (i = 0; i < MAX_REF_LF_DELTAS; i++) {
        Data = xd->ref_lf_deltas[i];

        // Frame level data
        if (xd->ref_lf_deltas[i] != xd->last_ref_lf_deltas[i]) {
          xd->last_ref_lf_deltas[i] = xd->ref_lf_deltas[i];
          vp8_write_bit(&header_bc, 1);

          if (Data > 0) {
            vp8_write_literal(&header_bc, (Data & 0x3F), 6);
            vp8_write_bit(&header_bc, 0);    // sign
          } else {
            Data = -Data;
            vp8_write_literal(&header_bc, (Data & 0x3F), 6);
            vp8_write_bit(&header_bc, 1);    // sign
          }
        } else {
          vp8_write_bit(&header_bc, 0);
        }
      }

      // Send update
      for (i = 0; i < MAX_MODE_LF_DELTAS; i++) {
        Data = xd->mode_lf_deltas[i];

        if (xd->mode_lf_deltas[i] != xd->last_mode_lf_deltas[i]) {
          xd->last_mode_lf_deltas[i] = xd->mode_lf_deltas[i];
          vp8_write_bit(&header_bc, 1);

          if (Data > 0) {
            vp8_write_literal(&header_bc, (Data & 0x3F), 6);
            vp8_write_bit(&header_bc, 0);    // sign
          } else {
            Data = -Data;
            vp8_write_literal(&header_bc, (Data & 0x3F), 6);
            vp8_write_bit(&header_bc, 1);    // sign
          }
        } else {
          vp8_write_bit(&header_bc, 0);
        }
      }
    }
  }

  // signal here is multi token partition is enabled
  // vp8_write_literal(&header_bc, pc->multi_token_partition, 2);
  vp8_write_literal(&header_bc, 0, 2);

  // Frame Q baseline quantizer index
  vp8_write_literal(&header_bc, pc->base_qindex, QINDEX_BITS);

  // Transmit Dc, Second order and Uv quantizer delta information
  put_delta_q(&header_bc, pc->y1dc_delta_q);
  put_delta_q(&header_bc, pc->y2dc_delta_q);
  put_delta_q(&header_bc, pc->y2ac_delta_q);
  put_delta_q(&header_bc, pc->uvdc_delta_q);
  put_delta_q(&header_bc, pc->uvac_delta_q);

  // When there is a key frame all reference buffers are updated using the new key frame
  if (pc->frame_type != KEY_FRAME) {
    // Should the GF or ARF be updated using the transmitted frame or buffer
    vp8_write_bit(&header_bc, pc->refresh_golden_frame);
    vp8_write_bit(&header_bc, pc->refresh_alt_ref_frame);

    // For inter frames the current default behavior is that when
    // cm->refresh_golden_frame is set we copy the old GF over to
    // the ARF buffer. This is purely an encoder decision at present.
    if (pc->refresh_golden_frame)
      pc->copy_buffer_to_arf  = 2;

    // If not being updated from current frame should either GF or ARF be updated from another buffer
    if (!pc->refresh_golden_frame)
      vp8_write_literal(&header_bc, pc->copy_buffer_to_gf, 2);

    if (!pc->refresh_alt_ref_frame)
      vp8_write_literal(&header_bc, pc->copy_buffer_to_arf, 2);

    // Indicate reference frame sign bias for Golden and ARF frames (always 0 for last frame buffer)
    vp8_write_bit(&header_bc, pc->ref_frame_sign_bias[GOLDEN_FRAME]);
    vp8_write_bit(&header_bc, pc->ref_frame_sign_bias[ALTREF_FRAME]);

    // Signal whether to allow high MV precision
    vp8_write_bit(&header_bc, (xd->allow_high_precision_mv) ? 1 : 0);
#if CONFIG_SWITCHABLE_INTERP
    if (pc->mcomp_filter_type == SWITCHABLE) {
      /* Check to see if only one of the filters is actually used */
      int count[VP8_SWITCHABLE_FILTERS];
      int i, j, c = 0;
      for (i = 0; i < VP8_SWITCHABLE_FILTERS; ++i) {
        count[i] = 0;
        for (j = 0; j <= VP8_SWITCHABLE_FILTERS; ++j) {
          count[i] += cpi->switchable_interp_count[j][i];
        }
        c += (count[i] > 0);
      }
      if (c == 1) {
        /* Only one filter is used. So set the filter at frame level */
        for (i = 0; i < VP8_SWITCHABLE_FILTERS; ++i) {
          if (count[i]) {
            pc->mcomp_filter_type = vp8_switchable_interp[i];
            break;
          }
        }
      }
    }
    // Signal the type of subpel filter to use
    vp8_write_bit(&header_bc, (pc->mcomp_filter_type == SWITCHABLE));
    if (pc->mcomp_filter_type != SWITCHABLE)
#endif  /* CONFIG_SWITCHABLE_INTERP */
      vp8_write_literal(&header_bc, (pc->mcomp_filter_type), 2);
  }

  vp8_write_bit(&header_bc, pc->refresh_entropy_probs);

  if (pc->frame_type != KEY_FRAME)
    vp8_write_bit(&header_bc, pc->refresh_last_frame);

#ifdef ENTROPY_STATS
  if (pc->frame_type == INTER_FRAME)
    active_section = 0;
  else
    active_section = 7;
#endif

  vp8_clear_system_state();  // __asm emms;

  vp8_copy(cpi->common.fc.pre_coef_probs, cpi->common.fc.coef_probs);
  vp8_copy(cpi->common.fc.pre_hybrid_coef_probs, cpi->common.fc.hybrid_coef_probs);
  vp8_copy(cpi->common.fc.pre_coef_probs_8x8, cpi->common.fc.coef_probs_8x8);
  vp8_copy(cpi->common.fc.pre_hybrid_coef_probs_8x8, cpi->common.fc.hybrid_coef_probs_8x8);
  vp8_copy(cpi->common.fc.pre_coef_probs_16x16, cpi->common.fc.coef_probs_16x16);
  vp8_copy(cpi->common.fc.pre_hybrid_coef_probs_16x16, cpi->common.fc.hybrid_coef_probs_16x16);
  vp8_copy(cpi->common.fc.pre_ymode_prob, cpi->common.fc.ymode_prob);
  vp8_copy(cpi->common.fc.pre_uv_mode_prob, cpi->common.fc.uv_mode_prob);
  vp8_copy(cpi->common.fc.pre_bmode_prob, cpi->common.fc.bmode_prob);
  vp8_copy(cpi->common.fc.pre_sub_mv_ref_prob, cpi->common.fc.sub_mv_ref_prob);
  vp8_copy(cpi->common.fc.pre_mbsplit_prob, cpi->common.fc.mbsplit_prob);
  vp8_copy(cpi->common.fc.pre_i8x8_mode_prob, cpi->common.fc.i8x8_mode_prob);
#if CONFIG_NEWMVENTROPY
  cpi->common.fc.pre_nmvc = cpi->common.fc.nmvc;
#else
  vp8_copy(cpi->common.fc.pre_mvc, cpi->common.fc.mvc);
  vp8_copy(cpi->common.fc.pre_mvc_hp, cpi->common.fc.mvc_hp);
#endif
  vp8_zero(cpi->sub_mv_ref_count);
  vp8_zero(cpi->mbsplit_count);
  vp8_zero(cpi->common.fc.mv_ref_ct)
  vp8_zero(cpi->common.fc.mv_ref_ct_a)

  update_coef_probs(cpi, &header_bc);

#ifdef ENTROPY_STATS
  active_section = 2;
#endif

  // Write out the mb_no_coeff_skip flag
  vp8_write_bit(&header_bc, pc->mb_no_coeff_skip);
  if (pc->mb_no_coeff_skip) {
    int k;

    update_skip_probs(cpi);
    for (k = 0; k < MBSKIP_CONTEXTS; ++k)
      vp8_write_literal(&header_bc, pc->mbskip_pred_probs[k], 8);
  }

  if (pc->frame_type == KEY_FRAME) {
    if (!pc->kf_ymode_probs_update) {
      vp8_write_literal(&header_bc, pc->kf_ymode_probs_index, 3);
    }
  } else {
    // Update the probabilities used to encode reference frame data
    update_ref_probs(cpi);

#ifdef ENTROPY_STATS
    active_section = 1;
#endif

#if CONFIG_PRED_FILTER
    // Write the prediction filter mode used for this frame
    vp8_write_literal(&header_bc, pc->pred_filter_mode, 2);

    // Write prediction filter on/off probability if signaling at MB level
    if (pc->pred_filter_mode == 2)
      vp8_write_literal(&header_bc, pc->prob_pred_filter_off, 8);

#endif
#if CONFIG_SWITCHABLE_INTERP
    if (pc->mcomp_filter_type == SWITCHABLE)
      update_switchable_interp_probs(cpi, &header_bc);
#endif

    vp8_write_literal(&header_bc, pc->prob_intra_coded, 8);
    vp8_write_literal(&header_bc, pc->prob_last_coded, 8);
    vp8_write_literal(&header_bc, pc->prob_gf_coded, 8);

    {
      const int comp_pred_mode = cpi->common.comp_pred_mode;
      const int use_compound_pred = (comp_pred_mode != SINGLE_PREDICTION_ONLY);
      const int use_hybrid_pred = (comp_pred_mode == HYBRID_PREDICTION);

      vp8_write(&header_bc, use_compound_pred, 128);
      if (use_compound_pred) {
        vp8_write(&header_bc, use_hybrid_pred, 128);
        if (use_hybrid_pred) {
          for (i = 0; i < COMP_PRED_CONTEXTS; i++) {
            pc->prob_comppred[i] = get_binary_prob(cpi->single_pred_count[i],
                                                   cpi->comp_pred_count[i]);
            vp8_write_literal(&header_bc, pc->prob_comppred[i], 8);
          }
        }
      }
    }

    update_mbintra_mode_probs(cpi, &header_bc);

#if CONFIG_NEWMVENTROPY
    vp8_write_nmvprobs(cpi, xd->allow_high_precision_mv, &header_bc);
#else
    if (xd->allow_high_precision_mv) {
      vp8_write_mvprobs_hp(cpi, &header_bc);
    } else {
      vp8_write_mvprobs(cpi, &header_bc);
    }
#endif
  }

  vp8_stop_encode(&header_bc);

  oh.first_partition_length_in_bytes = header_bc.pos;

  /* update frame tag */
  {
    int v = (oh.first_partition_length_in_bytes << 5) |
            (oh.show_frame << 4) |
            (oh.version << 1) |
            oh.type;

    dest[0] = v;
    dest[1] = v >> 8;
    dest[2] = v >> 16;
  }

  *size = VP8_HEADER_SIZE + extra_bytes_packed + header_bc.pos;
  vp8_start_encode(&residual_bc, cx_data + header_bc.pos);

  if (pc->frame_type == KEY_FRAME) {
    decide_kf_ymode_entropy(cpi);
    write_kfmodes(cpi, &residual_bc);
  } else {
    pack_inter_mode_mvs(cpi, &residual_bc);
    vp8_update_mode_context(&cpi->common);
  }


  vp8_stop_encode(&residual_bc);

  *size += residual_bc.pos;

}

#ifdef ENTROPY_STATS
void print_tree_update_probs() {
  int i, j, k, l;
  FILE *f = fopen("coefupdprob.h", "w");
  int Sum;
  fprintf(f, "\n/* Update probabilities for token entropy tree. */\n\n");

  fprintf(f, "const vp8_prob\n"
          "vp8_coef_update_probs[BLOCK_TYPES]\n"
          "                     [COEF_BANDS]\n"
          "                     [PREV_COEF_CONTEXTS]\n"
          "                     [ENTROPY_NODES] = {\n");
  for (i = 0; i < BLOCK_TYPES; i++) {
    fprintf(f, "  { \n");
    for (j = 0; j < COEF_BANDS; j++) {
      fprintf(f, "    {\n");
      for (k = 0; k < PREV_COEF_CONTEXTS; k++) {
        fprintf(f, "      {");
        for (l = 0; l < ENTROPY_NODES; l++) {
          fprintf(f, "%3ld, ",
              get_binary_prob(tree_update_hist[i][j][k][l][0],
                              tree_update_hist[i][j][k][l][1]));
        }
        fprintf(f, "},\n");
      }
      fprintf(f, "    },\n");
    }
    fprintf(f, "  },\n");
  }
  fprintf(f, "};\n");

  fprintf(f, "const vp8_prob\n"
          "vp8_coef_update_probs_8x8[BLOCK_TYPES_8X8]\n"
          "                         [COEF_BANDS]\n"
          "                         [PREV_COEF_CONTEXTS]\n"
          "                         [ENTROPY_NODES] = {\n");
  for (i = 0; i < BLOCK_TYPES_8X8; i++) {
    fprintf(f, "  { \n");
    for (j = 0; j < COEF_BANDS; j++) {
      fprintf(f, "    {\n");
      for (k = 0; k < PREV_COEF_CONTEXTS; k++) {
        fprintf(f, "      {");
        for (l = 0; l < MAX_ENTROPY_TOKENS - 1; l++) {
          fprintf(f, "%3ld, ",
              get_binary_prob(tree_update_hist_8x8[i][j][k][l][0],
                              tree_update_hist_8x8[i][j][k][l][1]));
        }
        fprintf(f, "},\n");
      }
      fprintf(f, "    },\n");
    }
    fprintf(f, "  },\n");
  }

  fprintf(f, "const vp8_prob\n"
          "vp8_coef_update_probs_16x16[BLOCK_TYPES_16X16]\n"
          "                           [COEF_BANDS]\n"
          "                           [PREV_COEF_CONTEXTS]\n"
          "                           [ENTROPY_NODES] = {\n");
  for (i = 0; i < BLOCK_TYPES_16X16; i++) {
    fprintf(f, "  { \n");
    for (j = 0; j < COEF_BANDS; j++) {
      fprintf(f, "    {\n");
      for (k = 0; k < PREV_COEF_CONTEXTS; k++) {
        fprintf(f, "      {");
        for (l = 0; l < MAX_ENTROPY_TOKENS - 1; l++) {
          fprintf(f, "%3ld, ",
              get_binary_prob(tree_update_hist_16x16[i][j][k][l][0],
                              tree_update_hist_16x16[i][j][k][l][1]));
        }
        fprintf(f, "},\n");
      }
      fprintf(f, "    },\n");
    }
    fprintf(f, "  },\n");
  }

  fclose(f);
  f = fopen("treeupdate.bin", "wb");
  fwrite(tree_update_hist, sizeof(tree_update_hist), 1, f);
  fwrite(tree_update_hist_8x8, sizeof(tree_update_hist_8x8), 1, f);
  fwrite(tree_update_hist_16x16, sizeof(tree_update_hist_16x16), 1, f);
  fclose(f);
}
#endif
