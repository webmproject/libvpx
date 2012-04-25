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

#include "vp8/common/defaultcoefcounts.h"

#include "vp8/common/seg_common.h"
#include "vp8/common/pred_common.h"
#include "vp8/common/entropy.h"

#if defined(SECTIONBITS_OUTPUT)
unsigned __int64 Sectionbits[500];
#endif

#ifdef ENTROPY_STATS
int intra_mode_stats[10][10][10];
static unsigned int tree_update_hist [BLOCK_TYPES] [COEF_BANDS] [PREV_COEF_CONTEXTS] [ENTROPY_NODES] [2];
static unsigned int tree_update_hist_8x8 [BLOCK_TYPES] [COEF_BANDS] [PREV_COEF_CONTEXTS] [ENTROPY_NODES] [2];

extern unsigned int active_section;
#endif

#ifdef MODE_STATS
int count_mb_seg[4] = { 0, 0, 0, 0 };
#endif

#define vp8_cost_upd  ((int)(vp8_cost_one(upd) - vp8_cost_zero(upd)) >> 8)
#define vp8_cost_upd256  ((int)(vp8_cost_one(upd) - vp8_cost_zero(upd)))

#if CONFIG_NEWUPDATE
#define SEARCH_NEWP
static int update_bits[255];

static void compute_update_table()
{
    int i;
    for (i=0; i<255; i++)
        update_bits[i] = vp8_count_term_subexp(i, SUBEXP_PARAM, 255);
}

static int remap_prob(int v, int m)
{
    const int n = 256;
    int i;
    if ((m<<1)<=n)
        i = recenter_nonneg(v, m) - 1;
    else
        i = recenter_nonneg(n-1-v, n-1-m) - 1;
    //if (i >= s) i -= s; else i = n - 2 - i;
    //i = ((i>>4)&15) | ((i((i>>4)&15) | ((i&15)<<4);&15)<<4);
    i = (i%17)*15 + (i/17);
    return i;
}
#endif

static void update_mode(
    vp8_writer *const w,
    int n,
    vp8_token tok               [/* n */],
    vp8_tree tree,
    vp8_prob Pnew               [/* n-1 */],
    vp8_prob Pcur               [/* n-1 */],
    unsigned int bct            [/* n-1 */] [2],
    const unsigned int num_events[/* n */]
)
{
    unsigned int new_b = 0, old_b = 0;
    int i = 0;

    vp8_tree_probs_from_distribution(
        n--, tok, tree,
        Pnew, bct, num_events,
        256, 1
    );

    do
    {
        new_b += vp8_cost_branch(bct[i], Pnew[i]);
        old_b += vp8_cost_branch(bct[i], Pcur[i]);
    }
    while (++i < n);

    if (new_b + (n << 8) < old_b)
    {
        int i = 0;

        vp8_write_bit(w, 1);

        do
        {
            const vp8_prob p = Pnew[i];

            vp8_write_literal(w, Pcur[i] = p ? p : 1, 8);
        }
        while (++i < n);
    }
    else
        vp8_write_bit(w, 0);
}

static void update_mbintra_mode_probs(VP8_COMP *cpi)
{
    VP8_COMMON *const x = & cpi->common;

    vp8_writer *const w = & cpi->bc;

    {
        vp8_prob Pnew   [VP8_YMODES-1];
        unsigned int bct [VP8_YMODES-1] [2];

        update_mode(
            w, VP8_YMODES, vp8_ymode_encodings, vp8_ymode_tree,
            Pnew, x->fc.ymode_prob, bct, (unsigned int *)cpi->ymode_count
        );
    }
}

void update_skip_probs(VP8_COMP *cpi)
{
#if CONFIG_NEWENTROPY
    VP8_COMMON *const pc = & cpi->common;
    int prob_skip_false[3] = {0, 0, 0};
    int k;

    for (k=0;k<MBSKIP_CONTEXTS;++k)
    {
        if ( (cpi->skip_false_count[k] + cpi->skip_true_count[k]) )
        {
            prob_skip_false[k] =
                cpi->skip_false_count[k] * 256 /
                (cpi->skip_false_count[k] + cpi->skip_true_count[k]);

            if (prob_skip_false[k] <= 1)
                prob_skip_false[k] = 1;

            if (prob_skip_false[k] > 255)
                prob_skip_false[k] = 255;
        }
        else
            prob_skip_false[k] = 128;

        pc->mbskip_pred_probs[k] = prob_skip_false[k];
    }

#else
    int prob_skip_false = 0;

    if ( (cpi->skip_false_count + cpi->skip_true_count) )
    {
        prob_skip_false = cpi->skip_false_count * 256 /
                          (cpi->skip_false_count + cpi->skip_true_count);

        if (prob_skip_false <= 1)
            prob_skip_false = 1;

        if (prob_skip_false > 255)
            prob_skip_false = 255;
    }
    else
        prob_skip_false = 128;

    cpi->prob_skip_false = prob_skip_false;

#endif
}

// This function updates the reference frame prediction stats
static void update_refpred_stats( VP8_COMP *cpi )
{
    VP8_COMMON *const cm = & cpi->common;
    MACROBLOCKD *const xd = & cpi->mb.e_mbd;

    int mb_row, mb_col;
    int i;
    int tot_count;
    int ref_pred_count[PREDICTION_PROBS][2];
    vp8_prob new_pred_probs[PREDICTION_PROBS];
    unsigned char pred_context;
    unsigned char pred_flag;

    int old_cost, new_cost;

    // Clear the prediction hit counters
    vpx_memset(ref_pred_count, 0, sizeof(ref_pred_count));

    // Set the prediction probability structures to defaults
    if ( cm->frame_type == KEY_FRAME )
    {
        // Set the prediction probabilities to defaults
        cm->ref_pred_probs[0] = 120;
        cm->ref_pred_probs[1] = 80;
        cm->ref_pred_probs[2] = 40;

        vpx_memset(cpi->ref_pred_probs_update, 0,
                   sizeof(cpi->ref_pred_probs_update) );
    }
    else
    {
        // For non-key frames.......

        // Scan through the macroblocks and collate prediction counts.
        xd->mode_info_context = cm->mi;
        for (mb_row = 0; mb_row < cm->mb_rows; mb_row++)
        {
            for (mb_col = 0; mb_col < cm->mb_cols; mb_col++)
            {
                // Get the prediction context and status
                pred_flag = get_pred_flag( xd, PRED_REF );
                pred_context = get_pred_context( cm, xd, PRED_REF );

                // Count prediction success
                ref_pred_count[pred_context][pred_flag]++;

                // Step on to the next mb
                xd->mode_info_context++;
            }

            // this is to account for the border in mode_info_context
            xd->mode_info_context++;
        }

        // From the prediction counts set the probabilities for each context
        for ( i = 0; i < PREDICTION_PROBS; i++ )
        {
            // MB reference frame not relevent to key frame encoding
            if ( cm->frame_type != KEY_FRAME )
            {
                // Work out the probabilities for the reference frame predictor
                tot_count = ref_pred_count[i][0] + ref_pred_count[i][1];
                if ( tot_count )
                {
                    new_pred_probs[i] =
                        ( ref_pred_count[i][0] * 255 ) / tot_count;

                    // Clamp to minimum allowed value
                    new_pred_probs[i] += !new_pred_probs[i];
                }
                else
                    new_pred_probs[i] = 128;
            }
            else
                new_pred_probs[i] = 128;

            // Decide whether or not to update the reference frame probs.
            // Returned costs are in 1/256 bit units.
            old_cost =
                (ref_pred_count[i][0] * vp8_cost_zero(cm->ref_pred_probs[i])) +
                (ref_pred_count[i][1] * vp8_cost_one(cm->ref_pred_probs[i]));

            new_cost =
                (ref_pred_count[i][0] * vp8_cost_zero(new_pred_probs[i])) +
                (ref_pred_count[i][1] * vp8_cost_one(new_pred_probs[i]));

            // Cost saving must be >= 8 bits (2048 in these units)
            if ( (old_cost - new_cost) >= 2048 )
            {
                cpi->ref_pred_probs_update[i] = 1;
                cm->ref_pred_probs[i] = new_pred_probs[i];
            }
            else
                cpi->ref_pred_probs_update[i] = 0;

        }
    }
}

static void write_ymode(vp8_writer *bc, int m, const vp8_prob *p)
{
    vp8_write_token(bc, vp8_ymode_tree, p, vp8_ymode_encodings + m);
}

static void kfwrite_ymode(vp8_writer *bc, int m, const vp8_prob *p)
{
    vp8_write_token(bc, vp8_kf_ymode_tree, p, vp8_kf_ymode_encodings + m);
}

static void write_i8x8_mode(vp8_writer *bc, int m, const vp8_prob *p)
{
    vp8_write_token(bc,vp8_i8x8_mode_tree, p, vp8_i8x8_mode_encodings + m);
}

static void write_uv_mode(vp8_writer *bc, int m, const vp8_prob *p)
{
    vp8_write_token(bc, vp8_uv_mode_tree, p, vp8_uv_mode_encodings + m);
}


static void write_bmode(vp8_writer *bc, int m, const vp8_prob *p)
{
    vp8_write_token(bc, vp8_bmode_tree, p, vp8_bmode_encodings + m);
}

static void write_split(vp8_writer *bc, int x)
{
    vp8_write_token(
        bc, vp8_mbsplit_tree, vp8_mbsplit_probs, vp8_mbsplit_encodings + x
    );
}

static int prob_update_savings(const unsigned int *ct,
                               const vp8_prob oldp, const vp8_prob newp,
                               const vp8_prob upd)
{
    const int old_b = vp8_cost_branch256(ct, oldp);
    const int new_b = vp8_cost_branch256(ct, newp);
#if CONFIG_NEWUPDATE
    const int delp = remap_prob(newp, oldp);
    const int update_b = update_bits[delp]*256 + vp8_cost_upd256;
#else
    const int update_b = 2048 + vp8_cost_upd256;
#endif
    return (old_b - new_b - update_b);
}

#if CONFIG_NEWUPDATE
static int prob_update_savings_search(const unsigned int *ct,
                                      const vp8_prob oldp, vp8_prob *bestp,
                                      const vp8_prob upd)
{
    const int old_b = vp8_cost_branch256(ct, oldp);
    int new_b, delp, update_b, savings, bestsavings, step;
    vp8_prob newp, bestnewp;

    bestsavings = 0;
    bestnewp = oldp;

    step = (*bestp > oldp ? -1 : 1);
    for (newp = *bestp; newp != oldp; newp+=step)
    {
        new_b = vp8_cost_branch256(ct, newp);
        delp = remap_prob(newp, oldp);
        update_b = update_bits[delp]*256 + vp8_cost_upd256;
        savings = old_b - new_b - update_b;
        if (savings > bestsavings)
        {
            bestsavings = savings;
            bestnewp = newp;
        }
    }
    *bestp = bestnewp;
    return bestsavings;
}
#endif

static void pack_tokens_c(vp8_writer *w, const TOKENEXTRA *p, int xcount)
{
    const TOKENEXTRA *const stop = p + xcount;
    unsigned int split;
    unsigned int shift;
    int count = w->count;
    unsigned int range = w->range;
    unsigned int lowvalue = w->lowvalue;

    while (p < stop)
    {
        const int t = p->Token;
        vp8_token *const a = vp8_coef_encodings + t;
        const vp8_extra_bit_struct *const b = vp8_extra_bits + t;
        int i = 0;
        const unsigned char *pp = p->context_tree;
        int v = a->value;
        int n = a->Len;

        /* skip one or two nodes */
        if (p->skip_eob_node)
        {
            n-=p->skip_eob_node;
            i = 2*p->skip_eob_node;
        }

        do
        {
            const int bb = (v >> --n) & 1;
            split = 1 + (((range - 1) * pp[i>>1]) >> 8);
            i = vp8_coef_tree[i+bb];

            if (bb)
            {
                lowvalue += split;
                range = range - split;
            }
            else
            {
                range = split;
            }

            shift = vp8_norm[range];
            range <<= shift;
            count += shift;

            if (count >= 0)
            {
                int offset = shift - count;

                if ((lowvalue << (offset - 1)) & 0x80000000)
                {
                    int x = w->pos - 1;

                    while (x >= 0 && w->buffer[x] == 0xff)
                    {
                        w->buffer[x] = (unsigned char)0;
                        x--;
                    }

                    w->buffer[x] += 1;
                }

                w->buffer[w->pos++] = (lowvalue >> (24 - offset));
                lowvalue <<= offset;
                shift = count;
                lowvalue &= 0xffffff;
                count -= 8 ;
            }

            lowvalue <<= shift;
        }
        while (n);


        if (b->base_val)
        {
            const int e = p->Extra, L = b->Len;

            if (L)
            {
                const unsigned char *pp = b->prob;
                int v = e >> 1;
                int n = L;              /* number of bits in v, assumed nonzero */
                int i = 0;

                do
                {
                    const int bb = (v >> --n) & 1;
                    split = 1 + (((range - 1) * pp[i>>1]) >> 8);
                    i = b->tree[i+bb];

                    if (bb)
                    {
                        lowvalue += split;
                        range = range - split;
                    }
                    else
                    {
                        range = split;
                    }

                    shift = vp8_norm[range];
                    range <<= shift;
                    count += shift;

                    if (count >= 0)
                    {
                        int offset = shift - count;

                        if ((lowvalue << (offset - 1)) & 0x80000000)
                        {
                            int x = w->pos - 1;

                            while (x >= 0 && w->buffer[x] == 0xff)
                            {
                                w->buffer[x] = (unsigned char)0;
                                x--;
                            }

                            w->buffer[x] += 1;
                        }

                        w->buffer[w->pos++] = (lowvalue >> (24 - offset));
                        lowvalue <<= offset;
                        shift = count;
                        lowvalue &= 0xffffff;
                        count -= 8 ;
                    }

                    lowvalue <<= shift;
                }
                while (n);
            }


            {

                split = (range + 1) >> 1;

                if (e & 1)
                {
                    lowvalue += split;
                    range = range - split;
                }
                else
                {
                    range = split;
                }

                range <<= 1;

                if ((lowvalue & 0x80000000))
                {
                    int x = w->pos - 1;

                    while (x >= 0 && w->buffer[x] == 0xff)
                    {
                        w->buffer[x] = (unsigned char)0;
                        x--;
                    }

                    w->buffer[x] += 1;

                }

                lowvalue  <<= 1;

                if (!++count)
                {
                    count = -8;
                    w->buffer[w->pos++] = (lowvalue >> 24);
                    lowvalue &= 0xffffff;
                }
            }

        }

        ++p;
    }

    w->count = count;
    w->lowvalue = lowvalue;
    w->range = range;

}

static void write_partition_size(unsigned char *cx_data, int size)
{
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
    vp8_writer *w, MB_PREDICTION_MODE m, const vp8_prob *p
)
{
#if CONFIG_DEBUG
    assert(NEARESTMV <= m  &&  m <= SPLITMV);
#endif
    vp8_write_token(w, vp8_mv_ref_tree, p,
                    vp8_mv_ref_encoding_array - NEARESTMV + m);
}

static void write_sub_mv_ref
(
    vp8_writer *w, B_PREDICTION_MODE m, const vp8_prob *p
)
{
#if CONFIG_DEBUG
    assert(LEFT4X4 <= m  &&  m <= NEW4X4);
#endif
    vp8_write_token(w, vp8_sub_mv_ref_tree, p,
                    vp8_sub_mv_ref_encoding_array - LEFT4X4 + m);
}

static void write_mv
(
    vp8_writer *w, const MV *mv, const int_mv *ref, const MV_CONTEXT *mvc
)
{
    MV e;
    e.row = mv->row - ref->as_mv.row;
    e.col = mv->col - ref->as_mv.col;

    vp8_encode_motion_vector(w, &e, mvc);
}

#if CONFIG_HIGH_PRECISION_MV
static void write_mv_hp
(
    vp8_writer *w, const MV *mv, const int_mv *ref, const MV_CONTEXT_HP *mvc
)
{
    MV e;
    e.row = mv->row - ref->as_mv.row;
    e.col = mv->col - ref->as_mv.col;

    vp8_encode_motion_vector_hp(w, &e, mvc);
}
#endif

// This function writes the current macro block's segnment id to the bitstream
// It should only be called if a segment map update is indicated.
static void write_mb_segid(vp8_writer *w,
                           const MB_MODE_INFO *mi, const MACROBLOCKD *x)
{
    // Encode the MB segment id.
    if (x->segmentation_enabled && x->update_mb_segmentation_map)
    {
        switch (mi->segment_id)
        {
        case 0:
            vp8_write(w, 0, x->mb_segment_tree_probs[0]);
            vp8_write(w, 0, x->mb_segment_tree_probs[1]);
            break;
        case 1:
            vp8_write(w, 0, x->mb_segment_tree_probs[0]);
            vp8_write(w, 1, x->mb_segment_tree_probs[1]);
            break;
        case 2:
            vp8_write(w, 1, x->mb_segment_tree_probs[0]);
            vp8_write(w, 0, x->mb_segment_tree_probs[2]);
            break;
        case 3:
            vp8_write(w, 1, x->mb_segment_tree_probs[0]);
            vp8_write(w, 1, x->mb_segment_tree_probs[2]);
            break;

            // TRAP.. This should not happen
        default:
            vp8_write(w, 0, x->mb_segment_tree_probs[0]);
            vp8_write(w, 0, x->mb_segment_tree_probs[1]);
            break;
        }
    }
}

// This function encodes the reference frame
static void encode_ref_frame( vp8_writer *const w,
                              VP8_COMMON *const cm,
                              MACROBLOCKD *xd,
                              int segment_id,
                              MV_REFERENCE_FRAME rf )
{
    int seg_ref_active;
    int seg_ref_count = 0;
    seg_ref_active = segfeature_active( xd,
                                        segment_id,
                                        SEG_LVL_REF_FRAME );

    if ( seg_ref_active )
    {
        seg_ref_count = check_segref( xd, segment_id, INTRA_FRAME ) +
                        check_segref( xd, segment_id, LAST_FRAME ) +
                        check_segref( xd, segment_id, GOLDEN_FRAME ) +
                        check_segref( xd, segment_id, ALTREF_FRAME );
    }

    // If segment level coding of this signal is disabled...
    // or the segment allows multiple reference frame options
    if ( !seg_ref_active || (seg_ref_count > 1) )
    {
        // Values used in prediction model coding
        unsigned char prediction_flag;
        vp8_prob pred_prob;
        MV_REFERENCE_FRAME pred_rf;

        // Get the context probability the prediction flag
        pred_prob = get_pred_prob( cm, xd, PRED_REF );

        // Get the predicted value.
        pred_rf = get_pred_ref( cm, xd );

        // Did the chosen reference frame match its predicted value.
        prediction_flag =
            ( xd->mode_info_context->mbmi.ref_frame == pred_rf );

        set_pred_flag( xd, PRED_REF, prediction_flag );
        vp8_write( w, prediction_flag, pred_prob );

        // If not predicted correctly then code value explicitly
        if ( !prediction_flag )
        {
            vp8_prob mod_refprobs[PREDICTION_PROBS];

            vpx_memcpy( mod_refprobs,
                        cm->mod_refprobs[pred_rf], sizeof(mod_refprobs) );

            // If segment coding enabled blank out options that cant occur by
            // setting the branch probability to 0.
            if ( seg_ref_active )
            {
                mod_refprobs[INTRA_FRAME] *=
                    check_segref( xd, segment_id, INTRA_FRAME );
                mod_refprobs[LAST_FRAME] *=
                    check_segref( xd, segment_id, LAST_FRAME );
                mod_refprobs[GOLDEN_FRAME] *=
                    ( check_segref( xd, segment_id, GOLDEN_FRAME ) *
                      check_segref( xd, segment_id, ALTREF_FRAME ) );
            }

            if ( mod_refprobs[0] )
            {
                vp8_write(w, (rf != INTRA_FRAME), mod_refprobs[0] );
            }

            // Inter coded
            if (rf != INTRA_FRAME)
            {
                if ( mod_refprobs[1] )
                {
                    vp8_write(w, (rf != LAST_FRAME), mod_refprobs[1] );
                }

                if (rf != LAST_FRAME)
                {
                    if ( mod_refprobs[2] )
                    {
                        vp8_write(w, (rf != GOLDEN_FRAME), mod_refprobs[2] );
                    }
                }
            }
        }
    }

    // if using the prediction mdoel we have nothing further to do because
    // the reference frame is fully coded by the segment
}

// Update the probabilities used to encode reference frame data
static void update_ref_probs( VP8_COMP *const cpi )
{
    VP8_COMMON *const cm = & cpi->common;

    const int *const rfct = cpi->count_mb_ref_frame_usage;
    const int rf_intra = rfct[INTRA_FRAME];
    const int rf_inter = rfct[LAST_FRAME] +
                         rfct[GOLDEN_FRAME] + rfct[ALTREF_FRAME];

    cm->prob_intra_coded = (rf_intra + rf_inter)
                            ? rf_intra * 255 / (rf_intra + rf_inter) : 1;

    if (!cm->prob_intra_coded)
        cm->prob_intra_coded = 1;

    cm->prob_last_coded = rf_inter ? (rfct[LAST_FRAME] * 255) / rf_inter : 128;

    if (!cm->prob_last_coded)
        cm->prob_last_coded = 1;

    cm->prob_gf_coded = (rfct[GOLDEN_FRAME] + rfct[ALTREF_FRAME])
                        ? (rfct[GOLDEN_FRAME] * 255) /
                          (rfct[GOLDEN_FRAME] + rfct[ALTREF_FRAME]) : 128;

    if (!cm->prob_gf_coded)
       cm->prob_gf_coded = 1;

    // Compute a modified set of probabilities to use when prediction of the
    // reference frame fails
    compute_mod_refprobs( cm );
}

static void pack_inter_mode_mvs(VP8_COMP *const cpi)
{
    int i;
    VP8_COMMON *const pc = & cpi->common;
    vp8_writer *const w = & cpi->bc;
    const MV_CONTEXT *mvc = pc->fc.mvc;
#if CONFIG_HIGH_PRECISION_MV
    const MV_CONTEXT_HP *mvc_hp = pc->fc.mvc_hp;
#endif
    MACROBLOCKD *xd = &cpi->mb.e_mbd;
    MODE_INFO *m;
    MODE_INFO *prev_m;

    const int mis = pc->mode_info_stride;
    int mb_row, mb_col;
    int row, col;

    // Values used in prediction model coding
    vp8_prob pred_prob;
    unsigned char prediction_flag;

    int row_delta[4] = { 0, +1,  0, -1};
    int col_delta[4] = {+1, -1, +1, +1};

    cpi->mb.partition_info = cpi->mb.pi;

    // Update the probabilities used to encode reference frame data
    update_ref_probs( cpi );

#ifdef ENTROPY_STATS
    active_section = 1;
#endif

    if (pc->mb_no_coeff_skip)
    {
#if CONFIG_NEWENTROPY
        int k;

        update_skip_probs( cpi );
        for (k=0;k<MBSKIP_CONTEXTS;++k)
            vp8_write_literal(w, pc->mbskip_pred_probs[k], 8);
#else
        update_skip_probs( cpi );
        vp8_write_literal(w, cpi->prob_skip_false, 8);
#endif
    }

    vp8_write_literal(w, pc->prob_intra_coded, 8);
    vp8_write_literal(w, pc->prob_last_coded, 8);
    vp8_write_literal(w, pc->prob_gf_coded, 8);

    if (cpi->common.comp_pred_mode == HYBRID_PREDICTION)
    {
        vp8_write(w, 1, 128);
        vp8_write(w, 1, 128);
        for (i = 0; i < COMP_PRED_CONTEXTS; i++)
        {
            if (cpi->single_pred_count[i] + cpi->comp_pred_count[i])
            {
                pc->prob_comppred[i] = cpi->single_pred_count[i] * 255 /
                    (cpi->single_pred_count[i] + cpi->comp_pred_count[i]);
                if (pc->prob_comppred[i] < 1)
                    pc->prob_comppred[i] = 1;
            }
            else
            {
                pc->prob_comppred[i] = 128;
            }
            vp8_write_literal(w, pc->prob_comppred[i], 8);
        }
    }
    else if (cpi->common.comp_pred_mode == SINGLE_PREDICTION_ONLY)
    {
        vp8_write(w, 0, 128);
    }
    else /* compound prediction only */
    {
        vp8_write(w, 1, 128);
        vp8_write(w, 0, 128);
    }

    update_mbintra_mode_probs(cpi);

#if CONFIG_HIGH_PRECISION_MV
    if (xd->allow_high_precision_mv)
        vp8_write_mvprobs_hp(cpi);
    else
#endif
    vp8_write_mvprobs(cpi);

    mb_row = 0;
    for (row=0; row < pc->mb_rows; row += 2)
    {
        m = pc->mi + row * mis;
        prev_m = pc->prev_mi + row * mis;

        mb_col = 0;
        for (col=0; col < pc->mb_cols; col += 2)
        {
            int i;

            // Process the 4 MBs in the order:
            // top-left, top-right, bottom-left, bottom-right
            for (i=0; i<4; i++)
            {
                const MB_MODE_INFO *const mi = & m->mbmi;
                const MV_REFERENCE_FRAME rf = mi->ref_frame;
                const MB_PREDICTION_MODE mode = mi->mode;
                const int segment_id = mi->segment_id;

                int dy = row_delta[i];
                int dx = col_delta[i];
                int offset_extended = dy * mis + dx;

                if ((mb_row >= pc->mb_rows) || (mb_col >= pc->mb_cols))
                {
                    // MB lies outside frame, move on
                    mb_row += dy;
                    mb_col += dx;
                    m += offset_extended;
                    prev_m += offset_extended;
                    cpi->mb.partition_info += offset_extended;
                    continue;
                }

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

                if (cpi->mb.e_mbd.update_mb_segmentation_map)
                {
                    // Is temporal coding of the segment map enabled
                    if (pc->temporal_update)
                    {
                        prediction_flag = get_pred_flag( xd, PRED_SEG_ID );
                        pred_prob = get_pred_prob( pc, xd, PRED_SEG_ID);

                        // Code the segment id prediction flag for this mb
                        vp8_write( w, prediction_flag, pred_prob );

                        // If the mb segment id wasn't predicted code explicitly
                        if (!prediction_flag)
                            write_mb_segid(w, mi, &cpi->mb.e_mbd);
                    }
                    else
                    {
                        // Normal unpredicted coding
                        write_mb_segid(w, mi, &cpi->mb.e_mbd);
                    }
                }

                if ( pc->mb_no_coeff_skip &&
                     ( !segfeature_active( xd, segment_id, SEG_LVL_EOB ) ||
                       ( get_segdata( xd, segment_id, SEG_LVL_EOB ) != 0 ) ) )
                {
#if CONFIG_NEWENTROPY
                    vp8_encode_bool(w, mi->mb_skip_coeff,
                                    get_pred_prob(pc, xd, PRED_MBSKIP));
#else
                vp8_encode_bool(w, mi->mb_skip_coeff, cpi->prob_skip_false);
#endif
                }

                // Encode the reference frame.
                encode_ref_frame( w, pc, xd, segment_id, rf );

                if (rf == INTRA_FRAME)
                {
#ifdef ENTROPY_STATS
                    active_section = 6;
#endif

                    if ( !segfeature_active( xd, segment_id, SEG_LVL_MODE ) )
                        write_ymode(w, mode, pc->fc.ymode_prob);

                    if (mode == B_PRED)
                    {
                        int j = 0;
#if CONFIG_COMP_INTRA_PRED
                        int uses_second =
                                m->bmi[0].as_mode.second !=
                                        (B_PREDICTION_MODE) (B_DC_PRED - 1);
                        vp8_write(w, uses_second, 128);
#endif
                        do {
#if CONFIG_COMP_INTRA_PRED
                            B_PREDICTION_MODE mode2 = m->bmi[j].as_mode.second;
#endif
                            write_bmode(w, m->bmi[j].as_mode.first,
                                        pc->fc.bmode_prob);
#if CONFIG_COMP_INTRA_PRED
                            if (uses_second)
                            {
                                write_bmode(w, mode2, pc->fc.bmode_prob);
                            }
#endif
                        } while (++j < 16);
                    }
                    if(mode == I8X8_PRED)
                    {
                        write_i8x8_mode(w, m->bmi[0].as_mode.first,
                                        pc->i8x8_mode_prob);
                        write_i8x8_mode(w, m->bmi[2].as_mode.first,
                                        pc->i8x8_mode_prob);
                        write_i8x8_mode(w, m->bmi[8].as_mode.first,
                                        pc->i8x8_mode_prob);
                        write_i8x8_mode(w, m->bmi[10].as_mode.first,
                                        pc->i8x8_mode_prob);
                    }
                    else
                    {
                        write_uv_mode(w, mi->uv_mode,
                                      pc->fc.uv_mode_prob[mode]);
#ifdef MODE_STATS
                        if(mode!=B_PRED)
                            ++cpi->y_uv_mode_count[mode][mi->uv_mode];
#endif
                    }
                }
                else
                {
                    int_mv best_mv, best_second_mv;
                    int ct[4];

                    vp8_prob mv_ref_p [VP8_MVREFS-1];

                    {
                        int_mv n1, n2;

                        vp8_find_near_mvs(xd, m, prev_m, &n1, &n2, &best_mv, ct,
                                          rf, cpi->common.ref_frame_sign_bias);
                        vp8_mv_ref_probs(&cpi->common, mv_ref_p, ct);

#ifdef ENTROPY_STATS
                        accum_mv_refs(mode, ct);
#endif
                    }

#ifdef ENTROPY_STATS
                    active_section = 3;
#endif

                    // Is the segment coding of mode enabled
                    if ( !segfeature_active( xd, segment_id, SEG_LVL_MODE ) )
                    {
                        write_mv_ref(w, mode, mv_ref_p);
                        vp8_accum_mv_refs(&cpi->common, mode, ct);
                    }

                    if (mi->second_ref_frame &&
                        (mode == NEWMV || mode == SPLITMV))
                    {
                        int_mv n1, n2;

                        vp8_find_near_mvs(xd, m,
                                          prev_m,
                                          &n1, &n2, &best_second_mv, ct,
                                          mi->second_ref_frame, cpi->common.ref_frame_sign_bias);
                    }

                    // does the feature use compound prediction or not
                    // (if not specified at the frame/segment level)
                    if (cpi->common.comp_pred_mode == HYBRID_PREDICTION)
                    {
                        vp8_write(w, mi->second_ref_frame != INTRA_FRAME,
                                  get_pred_prob( pc, xd, PRED_COMP ) );
                    }

                    {
                        switch (mode)   /* new, split require MVs */
                        {
                        case NEWMV:
#ifdef ENTROPY_STATS
                            active_section = 5;
#endif

#if CONFIG_HIGH_PRECISION_MV
                            if (xd->allow_high_precision_mv)
                                write_mv_hp(w, &mi->mv.as_mv, &best_mv, mvc_hp);
                            else
#endif
                            write_mv(w, &mi->mv.as_mv, &best_mv, mvc);

                            if (mi->second_ref_frame)
                            {
#if CONFIG_HIGH_PRECISION_MV
                                if (xd->allow_high_precision_mv)
                                    write_mv_hp(w, &mi->second_mv.as_mv,
                                                &best_second_mv, mvc_hp);
                                else
#endif
                                write_mv(w, &mi->second_mv.as_mv,
                                         &best_second_mv, mvc);
                            }
                            break;
                        case SPLITMV:
                        {
                            int j = 0;

#ifdef MODE_STATS
                            ++count_mb_seg [mi->partitioning];
#endif

                            write_split(w, mi->partitioning);

                            do
                            {
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

                                write_sub_mv_ref(w, blockmode,
                                               vp8_sub_mv_ref_prob2 [mv_contz]);

                                if (blockmode == NEW4X4)
                                {
#ifdef ENTROPY_STATS
                                    active_section = 11;
#endif
#if CONFIG_HIGH_PRECISION_MV
                                    if (xd->allow_high_precision_mv)
                                        write_mv_hp(w, &blockmv.as_mv, &best_mv,
                                                (const MV_CONTEXT_HP *) mvc_hp);
                                    else
#endif
                                    write_mv(w, &blockmv.as_mv, &best_mv,
                                             (const MV_CONTEXT *) mvc);

                                    if (mi->second_ref_frame)
                                    {
#if CONFIG_HIGH_PRECISION_MV
                                        if (xd->allow_high_precision_mv)
                                            write_mv_hp(w, &cpi->mb.partition_info->bmi[j].second_mv.as_mv,
                                                        &best_second_mv, (const MV_CONTEXT_HP *) mvc_hp);
                                        else
#endif
                                            write_mv(w, &cpi->mb.partition_info->bmi[j].second_mv.as_mv,
                                                     &best_second_mv, (const MV_CONTEXT *) mvc);
                                    }
                                }
                            }
                            while (++j < cpi->mb.partition_info->count);
                        }
                        break;
                        default:
                            break;
                        }
                    }
                }

                // Next MB
                mb_row += dy;
                mb_col += dx;
                m += offset_extended;
                prev_m += offset_extended;
                cpi->mb.partition_info += offset_extended;
#if CONFIG_DEBUG
                assert((prev_m-cpi->common.prev_mip)==(m-cpi->common.mip));
                assert((prev_m-cpi->common.prev_mi)==(m-cpi->common.mi));
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

static void write_kfmodes(VP8_COMP *cpi)
{
    vp8_writer *const bc = & cpi->bc;
    VP8_COMMON *const c = & cpi->common;
    const int mis = c->mode_info_stride;
    MACROBLOCKD *xd = &cpi->mb.e_mbd;
    MODE_INFO *m;
    int i;
    int row, col;
    int mb_row, mb_col;
#if CONFIG_NEWENTROPY
    int prob_skip_false[3] = {0, 0, 0};
#else
    int prob_skip_false = 0;
#endif
    int row_delta[4] = { 0, +1,  0, -1};
    int col_delta[4] = {+1, -1, +1, +1};

    if (c->mb_no_coeff_skip)
    {
        // Divide by 0 check. 0 case possible with segment features
#if CONFIG_NEWENTROPY
        int k;
        for (k=0;k<MBSKIP_CONTEXTS;++k)
        {
            if ( (cpi->skip_false_count[k] + cpi->skip_true_count[k]) )
            {
                prob_skip_false[k] = cpi->skip_false_count[k] * 256 /
                                  (cpi->skip_false_count[k] + cpi->skip_true_count[k]);

                if (prob_skip_false[k] <= 1)
                    prob_skip_false[k] = 1;

                if (prob_skip_false[k] > 255)
                    prob_skip_false[k] = 255;
            }
            else
                prob_skip_false[k] = 255;

            c->mbskip_pred_probs[k] = prob_skip_false[k];
            vp8_write_literal(bc, prob_skip_false[k], 8);
        }
#else
        if ( (cpi->skip_false_count + cpi->skip_true_count) )
        {
            prob_skip_false = cpi->skip_false_count * 256 /
                              (cpi->skip_false_count + cpi->skip_true_count);

            if (prob_skip_false <= 1)
                prob_skip_false = 1;

            if (prob_skip_false > 255)
                prob_skip_false = 255;
        }
        else
            prob_skip_false = 255;

        cpi->prob_skip_false = prob_skip_false;
        vp8_write_literal(bc, prob_skip_false, 8);
#endif
    }

#if CONFIG_QIMODE
    if(!c->kf_ymode_probs_update)
    {
        vp8_write_literal(bc, c->kf_ymode_probs_index, 3);
    }
#endif

    mb_row = 0;
    for (row=0; row < c->mb_rows; row += 2)
    {
        m = c->mi + row * mis;

        mb_col = 0;
        for (col=0; col < c->mb_cols; col += 2)
        {
            // Process the 4 MBs in the order:
            // top-left, top-right, bottom-left, bottom-right
            for (i=0; i<4; i++)
            {
                int ym;
                int segment_id;
                int dy = row_delta[i];
                int dx = col_delta[i];
                int offset_extended = dy * mis + dx;

                if ((mb_row >= c->mb_rows) || (mb_col >= c->mb_cols))
                {
                    // MB lies outside frame, move on
                    mb_row += dy;
                    mb_col += dx;
                    m += offset_extended;
                    continue;
                }

                // Make sure the MacroBlockD mode info pointer is set correctly
                xd->mode_info_context = m;

                ym = m->mbmi.mode;
                segment_id = m->mbmi.segment_id;

                if (cpi->mb.e_mbd.update_mb_segmentation_map)
                {
                    write_mb_segid(bc, &m->mbmi, &cpi->mb.e_mbd);
                }

                if ( c->mb_no_coeff_skip &&
                     ( !segfeature_active( xd, segment_id, SEG_LVL_EOB ) ||
                       (get_segdata( xd, segment_id, SEG_LVL_EOB ) != 0) ) )
                {
#if CONFIG_NEWENTROPY
                    vp8_encode_bool(bc, m->mbmi.mb_skip_coeff,
                                    get_pred_prob(c, xd, PRED_MBSKIP));
#else
                    vp8_encode_bool(bc, m->mbmi.mb_skip_coeff, prob_skip_false);
#endif
                }
#if CONFIG_QIMODE
                kfwrite_ymode(bc, ym,
                              c->kf_ymode_prob[c->kf_ymode_probs_index]);
#else
                kfwrite_ymode(bc, ym, c->kf_ymode_prob);
#endif

                if (ym == B_PRED)
                {
                    const int mis = c->mode_info_stride;
                    int i = 0;
#if CONFIG_COMP_INTRA_PRED
                    int uses_second =
                            m->bmi[0].as_mode.second !=
                                    (B_PREDICTION_MODE) (B_DC_PRED - 1);
                    vp8_write(bc, uses_second, 128);
#endif
                    do
                    {
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
#if CONFIG_COMP_INTRA_PRED
                        if (uses_second)
                        {
                            write_bmode(bc, bm2, c->kf_bmode_prob [A] [L]);
                        }
#endif
                    }
                    while (++i < 16);
                }
                if(ym == I8X8_PRED)
                {
                    write_i8x8_mode(bc, m->bmi[0].as_mode.first,
                                    c->i8x8_mode_prob);
                    write_i8x8_mode(bc, m->bmi[2].as_mode.first,
                                    c->i8x8_mode_prob);
                    write_i8x8_mode(bc, m->bmi[8].as_mode.first,
                                    c->i8x8_mode_prob);
                    write_i8x8_mode(bc, m->bmi[10].as_mode.first,
                                    c->i8x8_mode_prob);
                }
                else
                    write_uv_mode(bc, m->mbmi.uv_mode, c->kf_uv_mode_prob[ym]);

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
     coef_probs[BLOCK_TYPES][COEF_BANDS][PREV_COEF_CONTEXTS][ENTROPY_NODES])
{
    /* print coef probability tree */
    int i,j,k,l;
    FILE* f = fopen("enc_tree_probs.txt", "a");
    fprintf(f, "{\n");
    for (i = 0; i < BLOCK_TYPES; i++)
    {
        fprintf(f, "  {\n");
        for (j = 0; j < COEF_BANDS; j++)
        {
            fprintf(f, "    {\n");
            for (k = 0; k < PREV_COEF_CONTEXTS; k++)
            {
                fprintf(f, "      {");
                for (l = 0; l < ENTROPY_NODES; l++)
                {
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

static void sum_probs_over_prev_coef_context(
        const unsigned int probs[PREV_COEF_CONTEXTS][MAX_ENTROPY_TOKENS],
        unsigned int* out)
{
    int i, j;
    for (i=0; i < MAX_ENTROPY_TOKENS; ++i)
    {
        for (j=0; j < PREV_COEF_CONTEXTS; ++j)
        {
            const int tmp = out[i];
            out[i] += probs[j][i];
            /* check for wrap */
            if (out[i] < tmp)
                out[i] = UINT_MAX;
        }
    }
}

static int default_coef_context_savings(VP8_COMP *cpi)
{
    int savings = 0;
    int i = 0;
    do
    {
#if CONFIG_NEWUPDATE
        int j = !i;
#else
        int j = 0;      /* token/prob index */
#endif
        do
        {
            int k = 0;
            do
            {
                /* at every context */

                /* calc probs and branch cts for this frame only */
                //vp8_prob new_p           [ENTROPY_NODES];
                //unsigned int branch_ct   [ENTROPY_NODES] [2];

                int t = 0;      /* token/prob index */

                vp8_tree_probs_from_distribution(
                    MAX_ENTROPY_TOKENS, vp8_coef_encodings, vp8_coef_tree,
                    cpi->frame_coef_probs [i][j][k],
                    cpi->frame_branch_ct [i][j][k],
                    cpi->coef_counts [i][j][k],
                    256, 1
                );

                do
                {
                    const unsigned int *ct  = cpi->frame_branch_ct [i][j][k][t];
                    vp8_prob newp = cpi->frame_coef_probs [i][j][k][t];
                    const vp8_prob oldp = cpi->common.fc.coef_probs [i][j][k][t];
                    const vp8_prob upd = vp8_coef_update_probs [i][j][k][t];

#if CONFIG_NEWUPDATE && defined(SEARCH_NEWP)
                    const int s = prob_update_savings_search(ct, oldp, &newp, upd);
#else
                    const int s = prob_update_savings(ct, oldp, newp, upd);
#endif

                    if (s > 0)
                    {
                        savings += s;
                    }
                }
                while (++t < ENTROPY_NODES);
            }
            while (++k < PREV_COEF_CONTEXTS);
        }
        while (++j < COEF_BANDS);
    }
    while (++i < BLOCK_TYPES);
    return savings;
}

void build_coeff_contexts(VP8_COMP *cpi)
{
    int i = 0;
    do
    {
        int j = 0;
        do
        {
            int k = 0;
            do
            {
                int t;
                vp8_tree_probs_from_distribution(
                    MAX_ENTROPY_TOKENS, vp8_coef_encodings, vp8_coef_tree,
                    cpi->frame_coef_probs [i][j][k],
                    cpi->frame_branch_ct [i][j][k],
                    cpi->coef_counts [i][j][k],
                    256, 1
                );
#ifdef ENTROPY_STATS
                t = 0;
                do
                {
                    context_counters [i][j][k][t] += cpi->coef_counts [i][j][k][t];
                }
                while (++t < MAX_ENTROPY_TOKENS);
#endif
            }
            while (++k < PREV_COEF_CONTEXTS);
        }
        while (++j < COEF_BANDS);
    }
    while (++i < BLOCK_TYPES);


    i= 0;
    if(cpi->common.txfm_mode == ALLOW_8X8)
    {
        do
        {
            int j = 0;      /* token/prob index */
            do
            {
                int k = 0;
                do
                {
                    /* at every context */
                    /* calc probs and branch cts for this frame only */
                    //vp8_prob new_p           [ENTROPY_NODES];
                    //unsigned int branch_ct   [ENTROPY_NODES] [2];
                    int t = 0;      /* token/prob index */
                    vp8_tree_probs_from_distribution(
                        MAX_ENTROPY_TOKENS, vp8_coef_encodings, vp8_coef_tree,
                        cpi->frame_coef_probs_8x8 [i][j][k],
                        cpi->frame_branch_ct_8x8 [i][j][k],
                        cpi->coef_counts_8x8 [i][j][k],
                        256, 1
                        );
#ifdef ENTROPY_STATS
                    t = 0;
                    do
                    {
                        context_counters [i][j][k][t] += cpi->coef_counts [i][j][k][t];
                    }
                    while (++t < MAX_ENTROPY_TOKENS);
#endif

                }
                while (++k < PREV_COEF_CONTEXTS);
            }
            while (++j < COEF_BANDS);
        }
        while (++i < BLOCK_TYPES);
    }

}

#if CONFIG_NEWUPDATE
static void update_coef_probs3(VP8_COMP *cpi)
{
    const vp8_prob grpupd = 216;
    int i, j, k, t;
    vp8_writer *const w = & cpi->bc;
    int update[2];
    int savings;
    int bestupdndx[2*ENTROPY_NODES];

    vp8_clear_system_state(); //__asm emms;
    // Build the cofficient contexts based on counts collected in encode loop
    build_coeff_contexts(cpi);

    i = 0;
    for (i = 0; i < BLOCK_TYPES; ++i)
    {
        for (t = 0; t < ENTROPY_NODES; ++t)
        {
            /* dry run to see if there is any udpate at all needed */
            savings = 0;
            update[0] = update[1] = 0;
            for (j = !i; j < COEF_BANDS; ++j)
            {
                for (k = 0; k < PREV_COEF_CONTEXTS; ++k)
                {
                    vp8_prob newp = cpi->frame_coef_probs [i][j][k][t];
                    vp8_prob *Pold = cpi->common.fc.coef_probs [i][j][k] + t;
                    const vp8_prob upd = vp8_coef_update_probs [i][j][k][t];
                    int s;
                    int u = 0;

#if defined(SEARCH_NEWP)
                    s = prob_update_savings_search(
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
                    //printf("    %d %d %d: %d\n", i, j, k, u);
                    update[u]++;
                }
            }
            if (update[1] == 0 || savings < 0)
            {
                vp8_write(w, 0, grpupd);
                continue;
            }
            vp8_write(w, 1, grpupd);
            for (j = !i; j < COEF_BANDS; ++j)
            {
                for (k = 0; k < PREV_COEF_CONTEXTS; ++k)
                {
                    vp8_prob newp = cpi->frame_coef_probs [i][j][k][t];
                    vp8_prob *Pold = cpi->common.fc.coef_probs [i][j][k] + t;
                    const vp8_prob upd = vp8_coef_update_probs [i][j][k][t];
                    int s;
                    int u = 0;

#if defined(SEARCH_NEWP)
                    s = prob_update_savings_search(
                        cpi->frame_branch_ct [i][j][k][t], *Pold, &newp, upd);
                    if (s > 0 && newp != *Pold) u = 1;
#else
                    s = prob_update_savings(
                        cpi->frame_branch_ct [i][j][k][t], *Pold, newp, upd);
                    if (s > 0) u = 1;
#endif
                    //printf("  %d %d %d: %d (%d)\n", i, j, k, u, upd);
                    vp8_write(w, u, upd);
#ifdef ENTROPY_STATS
                    ++ tree_update_hist [i][j][k][t] [u];
#endif
                    if (u)
                    { /* send/use new probability */
                        int delp = remap_prob(newp, *Pold);
                        vp8_encode_term_subexp(w, delp, SUBEXP_PARAM, 255);
                        *Pold = newp;
                    }

                }
            }
        }
    }

    if(cpi->common.txfm_mode != ALLOW_8X8) return;

    for (i = 0; i < BLOCK_TYPES; ++i)
    {
        for (t = 0; t < ENTROPY_NODES; ++t)
        {
            /* dry run to see if there is any udpate at all needed */
            savings = 0;
            update[0] = update[1] = 0;
            for (j = !i; j < COEF_BANDS; ++j)
            {
                for (k = 0; k < PREV_COEF_CONTEXTS; ++k)
                {
                    vp8_prob newp = cpi->frame_coef_probs_8x8 [i][j][k][t];
                    vp8_prob *Pold = cpi->common.fc.coef_probs_8x8 [i][j][k] + t;
                    const vp8_prob upd = vp8_coef_update_probs_8x8 [i][j][k][t];
                    int s;
                    int u = 0;

#if defined(SEARCH_NEWP)
                    s = prob_update_savings_search(
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
            if (update[1] == 0 || savings < 0)
            {
                vp8_write(w, 0, grpupd);
                continue;
            }
            vp8_write(w, 1, grpupd);
            for (j = !i; j < COEF_BANDS; ++j)
            {
                for (k = 0; k < PREV_COEF_CONTEXTS; ++k)
                {
                    vp8_prob newp = cpi->frame_coef_probs_8x8 [i][j][k][t];
                    vp8_prob *Pold = cpi->common.fc.coef_probs_8x8 [i][j][k] + t;
                    const vp8_prob upd = vp8_coef_update_probs_8x8 [i][j][k][t];
                    int s;
                    int u = 0;

#if defined(SEARCH_NEWP)
                    s = prob_update_savings_search(
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
                    ++ tree_update_hist_8x8 [i][j][k][t] [u];
#endif
                    if (u)
                    {
                        /* send/use new probability */
                        int delp = remap_prob(newp, *Pold);
                        vp8_encode_term_subexp( w, delp, SUBEXP_PARAM, 255);
                        *Pold = newp;
                    }
                }
            }
        }
    }
}

static void update_coef_probs2(VP8_COMP *cpi)
{
    const vp8_prob grpupd = 192;
    int i, j, k, t;
    vp8_writer *const w = & cpi->bc;
    int update[2];
    int savings;
    int bestupdndx[2*ENTROPY_NODES];

    vp8_clear_system_state(); //__asm emms;
    // Build the cofficient contexts based on counts collected in encode loop
    build_coeff_contexts(cpi);

    for (t = 0; t < ENTROPY_NODES; ++t)
    {
        /* dry run to see if there is any udpate at all needed */
        savings = 0;
        update[0] = update[1] = 0;
        for (i = 0; i < BLOCK_TYPES; ++i)
        {
            for (j = !i; j < COEF_BANDS; ++j)
            {
                for (k = 0; k < PREV_COEF_CONTEXTS; ++k)
                {
                    vp8_prob newp = cpi->frame_coef_probs [i][j][k][t];
                    vp8_prob *Pold = cpi->common.fc.coef_probs [i][j][k] + t;
                    const vp8_prob upd = vp8_coef_update_probs [i][j][k][t];
                    int s;
                    int u = 0;

#if defined(SEARCH_NEWP)
                    s = prob_update_savings_search(
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
                    //printf("    %d %d %d: %d\n", i, j, k, u);
                    update[u]++;
                }
            }
        }
        if (update[1] == 0 || savings < 0)
        {
            vp8_write(w, 0, grpupd);
            continue;
        }
        vp8_write(w, 1, grpupd);
        for (i = 0; i < BLOCK_TYPES; ++i)
        {
            for (j = !i; j < COEF_BANDS; ++j)
            {
                for (k = 0; k < PREV_COEF_CONTEXTS; ++k)
                {
                    vp8_prob newp = cpi->frame_coef_probs [i][j][k][t];
                    vp8_prob *Pold = cpi->common.fc.coef_probs [i][j][k] + t;
                    const vp8_prob upd = vp8_coef_update_probs [i][j][k][t];
                    int s;
                    int u = 0;

#if defined(SEARCH_NEWP)
                    s = prob_update_savings_search(
                        cpi->frame_branch_ct [i][j][k][t], *Pold, &newp, upd);
                    if (s > 0 && newp != *Pold) u = 1;
#else
                    s = prob_update_savings(
                        cpi->frame_branch_ct [i][j][k][t], *Pold, newp, upd);
                    if (s > 0) u = 1;
#endif
                    //printf("  %d %d %d: %d (%d)\n", i, j, k, u, upd);
                    vp8_write(w, u, upd);
#ifdef ENTROPY_STATS
                    ++ tree_update_hist [i][j][k][t] [u];
#endif
                    if (u)
                    { /* send/use new probability */
                        int delp = remap_prob(newp, *Pold);
                        vp8_encode_term_subexp(w, delp, SUBEXP_PARAM, 255);
                        *Pold = newp;
                    }
                }
            }
        }
    }

    if(cpi->common.txfm_mode != ALLOW_8X8) return;

    for (t = 0; t < ENTROPY_NODES; ++t)
    {
        /* dry run to see if there is any udpate at all needed */
        savings = 0;
        update[0] = update[1] = 0;
        for (i = 0; i < BLOCK_TYPES; ++i)
        {
            for (j = !i; j < COEF_BANDS; ++j)
            {
                for (k = 0; k < PREV_COEF_CONTEXTS; ++k)
                {
                    vp8_prob newp = cpi->frame_coef_probs_8x8 [i][j][k][t];
                    vp8_prob *Pold = cpi->common.fc.coef_probs_8x8 [i][j][k] + t;
                    const vp8_prob upd = vp8_coef_update_probs_8x8 [i][j][k][t];
                    int s;
                    int u = 0;

#if defined(SEARCH_NEWP)
                    s = prob_update_savings_search(
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
        if (update[1] == 0 || savings < 0)
        {
            vp8_write(w, 0, grpupd);
            continue;
        }
        vp8_write(w, 1, grpupd);
        for (i = 0; i < BLOCK_TYPES; ++i)
        {
            for (j = !i; j < COEF_BANDS; ++j)
            {
                for (k = 0; k < PREV_COEF_CONTEXTS; ++k)
                {
                    vp8_prob newp = cpi->frame_coef_probs_8x8 [i][j][k][t];
                    vp8_prob *Pold = cpi->common.fc.coef_probs_8x8 [i][j][k] + t;
                    const vp8_prob upd = vp8_coef_update_probs_8x8 [i][j][k][t];
                    int s;
                    int u = 0;

#if defined(SEARCH_NEWP)
                        s = prob_update_savings_search(
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
                        ++ tree_update_hist_8x8 [i][j][k][t] [u];
#endif
                        if (u)
                        {
                            /* send/use new probability */
                            int delp = remap_prob(newp, *Pold);
                            vp8_encode_term_subexp( w, delp, SUBEXP_PARAM, 255);
                            *Pold = newp;
                        }
                }
            }
        }
    }
}
#endif

static void update_coef_probs(VP8_COMP *cpi)
{
    int i = 0;
    vp8_writer *const w = & cpi->bc;
    int update[2] = {0, 0};
    int savings;

    vp8_clear_system_state(); //__asm emms;

    // Build the cofficient contexts based on counts collected in encode loop
    build_coeff_contexts(cpi);

    //vp8_prob bestupd = find_coef_update_prob(cpi);

    /* dry run to see if there is any udpate at all needed */
    savings = 0;
    do
    {
#if CONFIG_NEWUPDATE
        int j = !i;
#else
        int j = 0;      /* token/prob index */
#endif
        do
        {
            int k = 0;
            int prev_coef_savings[ENTROPY_NODES] = {0};
            do
            {
                int t = 0;      /* token/prob index */
                do
                {
                    vp8_prob newp = cpi->frame_coef_probs [i][j][k][t];
                    vp8_prob *Pold = cpi->common.fc.coef_probs [i][j][k] + t;
                    const vp8_prob upd = vp8_coef_update_probs [i][j][k][t];
                    int s = prev_coef_savings[t];
                    int u = 0;

#if CONFIG_NEWUPDATE && defined(SEARCH_NEWP)
                    s = prob_update_savings_search(
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
                while (++t < ENTROPY_NODES);
            }
            while (++k < PREV_COEF_CONTEXTS);
        }
        while (++j < COEF_BANDS);
    }
    while (++i < BLOCK_TYPES);

    //printf("Update %d %d, savings %d\n", update[0], update[1], savings);
    /* Is coef updated at all */
#if CONFIG_NEWUPDATE
    if(update[1] == 0 || savings < 0)
#else
    if(update[1] == 0)
#endif
    {
        vp8_write_bit(w, 0);
    }
    else
    {
        vp8_write_bit(w, 1);
        i=0;
        do
        {
#if CONFIG_NEWUPDATE
            int j = !i;
#else
            int j = 0;      /* token/prob index */
#endif
            do
            {
                int k = 0;
                int prev_coef_savings[ENTROPY_NODES] = {0};

                do
                {
                    // calc probs and branch cts for this frame only
                    int t = 0;      /* token/prob index */
                    do
                    {
                        vp8_prob newp = cpi->frame_coef_probs [i][j][k][t];
                        vp8_prob *Pold = cpi->common.fc.coef_probs [i][j][k] + t;
                        const vp8_prob upd = vp8_coef_update_probs [i][j][k][t];
                        int s = prev_coef_savings[t];
                        int u = 0;

#if CONFIG_NEWUPDATE && defined(SEARCH_NEWP)
                        s = prob_update_savings_search(
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


                        vp8_write(w, u, upd);
#ifdef ENTROPY_STATS
                        ++ tree_update_hist [i][j][k][t] [u];
#endif
                        if (u)
                        {
                            /* send/use new probability */
#if CONFIG_NEWUPDATE
                            vp8_encode_term_subexp(
                                w, remap_prob(newp, *Pold), SUBEXP_PARAM, 255);
                            //printf("delp = %d/%d/%d\n", *Pold, remap_prob(newp, *Pold, 256), newp);
#else
                            vp8_write_literal(w, newp, 8);
#endif
                            *Pold = newp;
                        }
                    }
                    while (++t < ENTROPY_NODES);

                    // Accum token counts for generation of default statistics
#if 0//def ENTROPY_STATS
                    t = 0;
                    do
                    {
                        context_counters [i][j][k][t] += cpi->coef_counts [i][j][k][t];
                    }
                    while (++t < MAX_ENTROPY_TOKENS);
#endif
                }
                while (++k < PREV_COEF_CONTEXTS);
            }
            while (++j < COEF_BANDS);
        }
        while (++i < BLOCK_TYPES);
    }


    /* do not do this if not evena allowed */
    if(cpi->common.txfm_mode == ALLOW_8X8)
    {
        /* dry run to see if update is necessary */
        update[0] = update[1] = 0;
        savings = 0;
        i = 0;
        do
        {
#if CONFIG_NEWUPDATE
            int j = !i;
#else
            int j = 0;      /* token/prob index */
#endif
            do
            {
                int k = 0;
                do
                {
                    // calc probs and branch cts for this frame only
                    int t = 0;      /* token/prob index */
                    do
                    {
                        const unsigned int *ct  = cpi->frame_branch_ct_8x8 [i][j][k][t];
                        vp8_prob newp = cpi->frame_coef_probs_8x8 [i][j][k][t];
                        vp8_prob *Pold = cpi->common.fc.coef_probs_8x8 [i][j][k] + t;
                        const vp8_prob oldp = *Pold;
                        const vp8_prob upd = vp8_coef_update_probs_8x8 [i][j][k][t];
#if CONFIG_NEWUPDATE && defined(SEARCH_NEWP)
                        const int s = prob_update_savings_search(ct, oldp, &newp, upd);
                        const int u = s > 0 && newp != oldp ? 1 : 0;
                        if (u)
                            savings += s - (int)(vp8_cost_zero(upd));
                        else
                            savings -= (int)(vp8_cost_zero(upd));
#else
                        const int s = prob_update_savings(ct, oldp, newp, upd);
                        const int u = s > 0 ? 1 : 0;
                        if (u)
                            savings += s;
#endif

#ifdef ENTROPY_STATS
                        ++ tree_update_hist_8x8 [i][j][k][t] [u];
#endif
                        update[u]++;
                    }
                    while (++t < MAX_ENTROPY_TOKENS - 1);

                    // Accum token counts for generation of default statistics
#if 0//def ENTROPY_STATS
                    t = 0;

                    do
                    {
                        context_counters_8x8 [i][j][k][t] += cpi->coef_counts_8x8 [i][j][k][t];
                    }
                    while (++t < MAX_ENTROPY_TOKENS - 1);

#endif
                }
                while (++k < PREV_COEF_CONTEXTS);
            }
            while (++j < COEF_BANDS);
        }
        while (++i < BLOCK_TYPES);

#if CONFIG_NEWUPDATE
        if (update[1] == 0 || savings < 0)
#else
        if (update[1] == 0)
#endif
        {
            vp8_write_bit(w, 0);
        }
        else
        {
            vp8_write_bit(w, 1);
            i = 0;
            do
            {
#if CONFIG_NEWUPDATE
                int j = !i;
#else
                int j = 0;      /* token/prob index */
#endif
                do
                {
                    int k = 0;
                    do
                    {
                        int t = 0;      /* token/prob index */
                        do
                        {
                            const unsigned int *ct  = cpi->frame_branch_ct_8x8 [i][j][k][t];
                            vp8_prob newp = cpi->frame_coef_probs_8x8 [i][j][k][t];
                            vp8_prob *Pold = cpi->common.fc.coef_probs_8x8 [i][j][k] + t;
                            const vp8_prob oldp = *Pold;
                            const vp8_prob upd = vp8_coef_update_probs_8x8 [i][j][k][t];
#if CONFIG_NEWUPDATE && defined(SEARCH_NEWP)
                            const int s = prob_update_savings_search(ct, oldp, &newp, upd);
                            const int u = s > 0 && newp != oldp ? 1 : 0;
#else
                            const int s = prob_update_savings(ct, oldp, newp, upd);
                            const int u = s > 0 ? 1 : 0;
#endif
                            vp8_write(w, u, upd);
#ifdef ENTROPY_STATS
                            ++ tree_update_hist_8x8 [i][j][k][t] [u];
#endif
                            if (u)
                            {
                                /* send/use new probability */
#if CONFIG_NEWUPDATE
                                vp8_encode_term_subexp(
                                    w, remap_prob(newp, oldp), SUBEXP_PARAM, 255);
#else
                                vp8_write_literal(w, newp, 8);
#endif
                                *Pold = newp;
                            }
                        }
                        while (++t < MAX_ENTROPY_TOKENS - 1);
                        // Accum token counts for generation of default statistics
#if 0//def ENTROPY_STATS
                        t = 0;
                        do
                        {
                            context_counters_8x8 [i][j][k][t] += cpi->coef_counts_8x8 [i][j][k][t];
                        }
                        while (++t < MAX_ENTROPY_TOKENS);
#endif
                    }
                    while (++k < PREV_COEF_CONTEXTS);
                }
                while (++j < COEF_BANDS);
            }
            while (++i < BLOCK_TYPES);
        }
    }
}

#ifdef PACKET_TESTING
FILE *vpxlogc = 0;
#endif

static void put_delta_q(vp8_writer *bc, int delta_q)
{
    if (delta_q != 0)
    {
        vp8_write_bit(bc, 1);
        vp8_write_literal(bc, abs(delta_q), 4);

        if (delta_q < 0)
            vp8_write_bit(bc, 1);
        else
            vp8_write_bit(bc, 0);
    }
    else
        vp8_write_bit(bc, 0);
}
#if CONFIG_QIMODE
extern const unsigned int kf_y_mode_cts[8][VP8_YMODES];
static void decide_kf_ymode_entropy(VP8_COMP *cpi)
{

    int mode_cost[MB_MODE_COUNT];
    int cost;
    int bestcost = INT_MAX;
    int bestindex = 0;
    int i, j;

    for(i=0; i<8; i++)
    {
        vp8_cost_tokens(mode_cost, cpi->common.kf_ymode_prob[i], vp8_kf_ymode_tree);
        cost = 0;
        for(j=0;j<VP8_YMODES;j++)
        {
            cost += mode_cost[j] * cpi->ymode_count[j];
        }
        if(cost < bestcost)
        {
            bestindex = i;
            bestcost = cost;
        }
    }
    cpi->common.kf_ymode_probs_index = bestindex;

}
#endif
static void segment_reference_frames(VP8_COMP *cpi)
{
    VP8_COMMON *oci = &cpi->common;
    MODE_INFO *mi = oci->mi;
    int ref[MAX_MB_SEGMENTS]={0};
    int i,j;
    int mb_index=0;
    MACROBLOCKD *const xd = & cpi->mb.e_mbd;

    for (i = 0; i < oci->mb_rows; i++)
    {
        for (j = 0; j < oci->mb_cols; j++, mb_index++)
        {
            ref[mi[mb_index].mbmi.segment_id]|=(1<<mi[mb_index].mbmi.ref_frame);
        }
        mb_index++;
    }
    for (i = 0; i < MAX_MB_SEGMENTS; i++)
    {
        enable_segfeature(xd,i,SEG_LVL_REF_FRAME);
        set_segdata( xd,i, SEG_LVL_REF_FRAME, ref[i]);
    }


}

void vp8_pack_bitstream(VP8_COMP *cpi, unsigned char *dest, unsigned long *size)
{
    int i, j;
    VP8_HEADER oh;
    VP8_COMMON *const pc = & cpi->common;
    vp8_writer *const bc = & cpi->bc;
    MACROBLOCKD *const xd = & cpi->mb.e_mbd;
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

#if CONFIG_NEWUPDATE
    compute_update_table();
#endif

    //vp8_kf_default_bmode_probs() is called in vp8_setup_key_frame() once for each
    //K frame before encode frame. pc->kf_bmode_prob doesn't get changed anywhere
    //else. No need to call it again here. --yw
    //vp8_kf_default_bmode_probs( pc->kf_bmode_prob);

    // every keyframe send startcode, width, height, scale factor, clamp and color type
    if (oh.type == KEY_FRAME)
    {
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
        cx_data += extra_bytes_packed ;

        vp8_start_encode(bc, cx_data);

        // signal clr type
        vp8_write_bit(bc, pc->clr_type);
        vp8_write_bit(bc, pc->clamp_type);

    }
    else
        vp8_start_encode(bc, cx_data);

    // Signal whether or not Segmentation is enabled
    vp8_write_bit(bc, (xd->segmentation_enabled) ? 1 : 0);

    // Indicate which features are enabled
    if ( xd->segmentation_enabled )
    {
        // Indicate whether or not the segmentation map is being updated.
        vp8_write_bit(bc, (xd->update_mb_segmentation_map) ? 1 : 0);

        // If it is, then indicate the method that will be used.
        if ( xd->update_mb_segmentation_map )
        {
            // Select the coding strategy (temporal or spatial)
            choose_segmap_coding_method( cpi );

            // Take a copy of the segment map if it changed for
            // future comparison
            vpx_memcpy( pc->last_frame_seg_map,
                        cpi->segmentation_map, pc->MBs );

            // Write out the chosen coding method.
            vp8_write_bit(bc, (pc->temporal_update) ? 1:0);
        }

        vp8_write_bit(bc, (xd->update_mb_segmentation_data) ? 1 : 0);

        //segment_reference_frames(cpi);

        if (xd->update_mb_segmentation_data)
        {
            signed char Data;

            vp8_write_bit(bc, (xd->mb_segment_abs_delta) ? 1 : 0);

            // For each segments id...
            for (i = 0; i < MAX_MB_SEGMENTS; i++)
            {
                // For each segmentation codable feature...
                for (j = 0; j < SEG_LVL_MAX; j++)
                {
                    Data = get_segdata( xd, i, j );


#if CONFIG_FEATUREUPDATES

                    // check if there's an update
                    if(segfeature_changed( xd,i,j) )
                    {
                        vp8_write_bit(bc, 1);

                        if ( segfeature_active( xd, i, j ) )
                        {
                            // this bit is to say we are still
                            // active/  if we were inactive
                            // this is unnecessary
                            if ( old_segfeature_active( xd, i, j ))
                            {
                                vp8_write_bit(bc, 1);
                            }
                            // Is the segment data signed..
                            if ( is_segfeature_signed(j) )
                            {
                                // Encode the relevant feature data
                                if (Data < 0)
                                {
                                    Data = - Data;
                                    vp8_write_literal(bc, Data,
                                            seg_feature_data_bits(j));
                                    vp8_write_bit(bc, 1);
                                }
                                else
                                {
                                    vp8_write_literal(bc, Data,
                                            seg_feature_data_bits(j));
                                    vp8_write_bit(bc, 0);
                                }
                            }
                            // Unsigned data element so no sign bit needed
                            else
                            vp8_write_literal(bc, Data,
                                    seg_feature_data_bits(j));
                        }
                        // feature is inactive now
                        else if ( old_segfeature_active( xd, i, j ))
                        {
                           vp8_write_bit(bc, 0);
                        }
                    }
                    else
                    {
                        vp8_write_bit(bc,0);
                    }
#else

                    // If the feature is enabled...
                    if ( segfeature_active( xd, i, j ) )
                    {
                        vp8_write_bit(bc, 1);

                        // Is the segment data signed..
                        if ( is_segfeature_signed(j) )
                        {
                            // Encode the relevant feature data
                            if (Data < 0)
                            {
                                Data = - Data;
                                vp8_write_literal(bc, Data,
                                        seg_feature_data_bits(j));
                                vp8_write_bit(bc, 1);
                            }
                            else
                            {
                                vp8_write_literal(bc, Data,
                                        seg_feature_data_bits(j));
                                vp8_write_bit(bc, 0);
                            }
                        }
                        // Unsigned data element so no sign bit needed
                        else
                            vp8_write_literal(bc, Data,
                                    seg_feature_data_bits(j));
                    }
                    else
                        vp8_write_bit(bc, 0);
#endif
                }
            }
        }

#if CONFIG_FEATUREUPDATES
        // save the segment info for updates next frame
        save_segment_info ( xd );
#endif

        if (xd->update_mb_segmentation_map)
        {
            // Send the tree probabilities used to decode unpredicted
            // macro-block segments
            for (i = 0; i < MB_FEATURE_TREE_PROBS; i++)
            {
                int Data = xd->mb_segment_tree_probs[i];

                if (Data != 255)
                {
                    vp8_write_bit(bc, 1);
                    vp8_write_literal(bc, Data, 8);
                }
                else
                    vp8_write_bit(bc, 0);
            }

            // If predictive coding of segment map is enabled send the
            // prediction probabilities.
            if ( pc->temporal_update )
            {
                for (i = 0; i < PREDICTION_PROBS; i++)
                {
                    int Data = pc->segment_pred_probs[i];

                    if (Data != 255)
                    {
                        vp8_write_bit(bc, 1);
                        vp8_write_literal(bc, Data, 8);
                    }
                    else
                        vp8_write_bit(bc, 0);
                }
            }
        }
    }

    // Encode the common prediction model status flag probability updates for
    // the reference frame
    update_refpred_stats( cpi );
    if ( pc->frame_type != KEY_FRAME )
    {
        for (i = 0; i < PREDICTION_PROBS; i++)
        {
            if ( cpi->ref_pred_probs_update[i] )
            {
                vp8_write_bit(bc, 1);
                vp8_write_literal(bc, pc->ref_pred_probs[i], 8);
            }
            else
                vp8_write_bit(bc, 0);
        }
    }

    vp8_write_bit(bc, pc->txfm_mode);

    // Encode the loop filter level and type
    vp8_write_bit(bc, pc->filter_type);
    vp8_write_literal(bc, pc->filter_level, 6);
    vp8_write_literal(bc, pc->sharpness_level, 3);

    // Write out loop filter deltas applied at the MB level based on mode or ref frame (if they are enabled).
    vp8_write_bit(bc, (xd->mode_ref_lf_delta_enabled) ? 1 : 0);

    if (xd->mode_ref_lf_delta_enabled)
    {
        // Do the deltas need to be updated
        int send_update = xd->mode_ref_lf_delta_update;

        vp8_write_bit(bc, send_update);
        if (send_update)
        {
            int Data;

            // Send update
            for (i = 0; i < MAX_REF_LF_DELTAS; i++)
            {
                Data = xd->ref_lf_deltas[i];

                // Frame level data
                if (xd->ref_lf_deltas[i] != xd->last_ref_lf_deltas[i])
                {
                    xd->last_ref_lf_deltas[i] = xd->ref_lf_deltas[i];
                    vp8_write_bit(bc, 1);

                    if (Data > 0)
                    {
                        vp8_write_literal(bc, (Data & 0x3F), 6);
                        vp8_write_bit(bc, 0);    // sign
                    }
                    else
                    {
                        Data = -Data;
                        vp8_write_literal(bc, (Data & 0x3F), 6);
                        vp8_write_bit(bc, 1);    // sign
                    }
                }
                else
                    vp8_write_bit(bc, 0);
            }

            // Send update
            for (i = 0; i < MAX_MODE_LF_DELTAS; i++)
            {
                Data = xd->mode_lf_deltas[i];

                if (xd->mode_lf_deltas[i] != xd->last_mode_lf_deltas[i])
                {
                    xd->last_mode_lf_deltas[i] = xd->mode_lf_deltas[i];
                    vp8_write_bit(bc, 1);

                    if (Data > 0)
                    {
                        vp8_write_literal(bc, (Data & 0x3F), 6);
                        vp8_write_bit(bc, 0);    // sign
                    }
                    else
                    {
                        Data = -Data;
                        vp8_write_literal(bc, (Data & 0x3F), 6);
                        vp8_write_bit(bc, 1);    // sign
                    }
                }
                else
                    vp8_write_bit(bc, 0);
            }
        }
    }

    //signal here is multi token partition is enabled
    //vp8_write_literal(bc, pc->multi_token_partition, 2);
    vp8_write_literal(bc, 0, 2);

    // Frame Q baseline quantizer index
    vp8_write_literal(bc, pc->base_qindex, QINDEX_BITS);

    // Transmit Dc, Second order and Uv quantizer delta information
    put_delta_q(bc, pc->y1dc_delta_q);
    put_delta_q(bc, pc->y2dc_delta_q);
    put_delta_q(bc, pc->y2ac_delta_q);
    put_delta_q(bc, pc->uvdc_delta_q);
    put_delta_q(bc, pc->uvac_delta_q);

    // When there is a key frame all reference buffers are updated using the new key frame
    if (pc->frame_type != KEY_FRAME)
    {
        // Should the GF or ARF be updated using the transmitted frame or buffer
        vp8_write_bit(bc, pc->refresh_golden_frame);
        vp8_write_bit(bc, pc->refresh_alt_ref_frame);

        // For inter frames the current default behavior is that when
        // cm->refresh_golden_frame is set we copy the old GF over to
        // the ARF buffer. This is purely an encoder decision at present.
        if (pc->refresh_golden_frame)
            pc->copy_buffer_to_arf  = 2;

        // If not being updated from current frame should either GF or ARF be updated from another buffer
        if (!pc->refresh_golden_frame)
            vp8_write_literal(bc, pc->copy_buffer_to_gf, 2);

        if (!pc->refresh_alt_ref_frame)
            vp8_write_literal(bc, pc->copy_buffer_to_arf, 2);

        // Indicate reference frame sign bias for Golden and ARF frames (always 0 for last frame buffer)
        vp8_write_bit(bc, pc->ref_frame_sign_bias[GOLDEN_FRAME]);
        vp8_write_bit(bc, pc->ref_frame_sign_bias[ALTREF_FRAME]);

#if CONFIG_HIGH_PRECISION_MV
        // Signal whether to allow high MV precision
        vp8_write_bit(bc, (xd->allow_high_precision_mv) ? 1 : 0);
#endif
#if CONFIG_ENHANCED_INTERP
        // Signal the type of subpel filter to use
        vp8_write_literal(bc, (pc->mcomp_filter_type), 2);
#endif
    }

    vp8_write_bit(bc, pc->refresh_entropy_probs);

    if (pc->frame_type != KEY_FRAME)
        vp8_write_bit(bc, pc->refresh_last_frame);

#ifdef ENTROPY_STATS

    if (pc->frame_type == INTER_FRAME)
        active_section = 0;
    else
        active_section = 7;

#endif

    vp8_clear_system_state();  //__asm emms;

#if COEFUPDATETYPE == 2
    update_coef_probs2(cpi);
#elif COEFUPDATETYPE == 3
    update_coef_probs3(cpi);
#else
    update_coef_probs(cpi);
#endif

#ifdef ENTROPY_STATS
    active_section = 2;
#endif

    // Write out the mb_no_coeff_skip flag
    vp8_write_bit(bc, pc->mb_no_coeff_skip);

    if (pc->frame_type == KEY_FRAME)
    {
#if CONFIG_QIMODE
        decide_kf_ymode_entropy(cpi);
#endif
        write_kfmodes(cpi);

#ifdef ENTROPY_STATS
        active_section = 8;
#endif
    }
    else
    {
        pack_inter_mode_mvs(cpi);

        vp8_update_mode_context(&cpi->common);

#ifdef ENTROPY_STATS
        active_section = 1;
#endif
    }

    vp8_stop_encode(bc);

    oh.first_partition_length_in_bytes = cpi->bc.pos;

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

    *size = VP8_HEADER_SIZE + extra_bytes_packed + cpi->bc.pos;

    vp8_start_encode(&cpi->bc2, cx_data + bc->pos);

    pack_tokens(&cpi->bc2, cpi->tok, cpi->tok_count);

    vp8_stop_encode(&cpi->bc2);

    *size += cpi->bc2.pos;
}

#ifdef ENTROPY_STATS
void print_tree_update_probs()
{
    int i, j, k, l;
    FILE *f = fopen("context.c", "a");
    int Sum;
    fprintf(f, "\n/* Update probabilities for token entropy tree. */\n\n");
    fprintf(f, "const vp8_prob tree_update_probs[BLOCK_TYPES] [COEF_BANDS] [PREV_COEF_CONTEXTS] [ENTROPY_NODES] = {\n");

    for (i = 0; i < BLOCK_TYPES; i++)
    {
        fprintf(f, "  { \n");

        for (j = 0; j < COEF_BANDS; j++)
        {
            fprintf(f, "    {\n");

            for (k = 0; k < PREV_COEF_CONTEXTS; k++)
            {
                fprintf(f, "      {");

                for (l = 0; l < ENTROPY_NODES; l++)
                {
                    Sum = tree_update_hist[i][j][k][l][0] + tree_update_hist[i][j][k][l][1];

                    if (Sum > 0)
                    {
                        if (((tree_update_hist[i][j][k][l][0] * 255) / Sum) > 0)
                            fprintf(f, "%3ld, ", (tree_update_hist[i][j][k][l][0] * 255) / Sum);
                        else
                            fprintf(f, "%3ld, ", 1);
                    }
                    else
                        fprintf(f, "%3ld, ", 128);
                }

                fprintf(f, "},\n");
            }

            fprintf(f, "    },\n");
        }

        fprintf(f, "  },\n");
    }

    fprintf(f, "};\n");

    fprintf(f, "const vp8_prob tree_update_probs_8x8[BLOCK_TYPES] [COEF_BANDS] [PREV_COEF_CONTEXTS] [ENTROPY_NODES] = {\n");

    for (i = 0; i < BLOCK_TYPES; i++)
    {
        fprintf(f, "  { \n");

        for (j = 0; j < COEF_BANDS; j++)
        {
            fprintf(f, "    {\n");

            for (k = 0; k < PREV_COEF_CONTEXTS; k++)
            {
                fprintf(f, "      {");

                for (l = 0; l < MAX_ENTROPY_TOKENS - 1; l++)
                {
                    Sum = tree_update_hist_8x8[i][j][k][l][0] + tree_update_hist_8x8[i][j][k][l][1];

                    if (Sum > 0)
                    {
                        if (((tree_update_hist_8x8[i][j][k][l][0] * 255) / Sum) > 0)
                            fprintf(f, "%3ld, ", (tree_update_hist_8x8[i][j][k][l][0] * 255) / Sum);
                        else
                            fprintf(f, "%3ld, ", 1);
                    }
                    else
                        fprintf(f, "%3ld, ", 128);
                }

                fprintf(f, "},\n");
            }

            fprintf(f, "    },\n");
        }

        fprintf(f, "  },\n");
    }
    fclose(f);
}
#endif
