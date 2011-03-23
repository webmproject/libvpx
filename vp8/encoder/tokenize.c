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
#include "onyx_int.h"
#include "tokenize.h"
#include "vpx_mem/vpx_mem.h"

/* Global event counters used for accumulating statistics across several
   compressions, then generating context.c = initial stats. */

#ifdef ENTROPY_STATS
_int64 context_counters[BLOCK_TYPES] [COEF_BANDS] [PREV_COEF_CONTEXTS] [vp8_coef_tokens];
#endif
void vp8_stuff_mb(VP8_COMP *cpi, MACROBLOCKD *x, TOKENEXTRA **t) ;
void vp8_fix_contexts(MACROBLOCKD *x);

static TOKENVALUE dct_value_tokens[DCT_MAX_VALUE*2];
const TOKENVALUE *vp8_dct_value_tokens_ptr;
static int dct_value_cost[DCT_MAX_VALUE*2];
const int *vp8_dct_value_cost_ptr;
#if 0
int skip_true_count = 0;
int skip_false_count = 0;
#endif
static void fill_value_tokens()
{

    TOKENVALUE *const t = dct_value_tokens + DCT_MAX_VALUE;
    vp8_extra_bit_struct *const e = vp8_extra_bits;

    int i = -DCT_MAX_VALUE;
    int sign = 1;

    do
    {
        if (!i)
            sign = 0;

        {
            const int a = sign ? -i : i;
            int eb = sign;

            if (a > 4)
            {
                int j = 4;

                while (++j < 11  &&  e[j].base_val <= a) {}

                t[i].Token = --j;
                eb |= (a - e[j].base_val) << 1;
            }
            else
                t[i].Token = a;

            t[i].Extra = eb;
        }

        // initialize the cost for extra bits for all possible coefficient value.
        {
            int cost = 0;
            vp8_extra_bit_struct *p = vp8_extra_bits + t[i].Token;

            if (p->base_val)
            {
                const int extra = t[i].Extra;
                const int Length = p->Len;

                if (Length)
                    cost += vp8_treed_cost(p->tree, p->prob, extra >> 1, Length);

                cost += vp8_cost_bit(vp8_prob_half, extra & 1); /* sign */
                dct_value_cost[i + DCT_MAX_VALUE] = cost;
            }

        }

    }
    while (++i < DCT_MAX_VALUE);

    vp8_dct_value_tokens_ptr = dct_value_tokens + DCT_MAX_VALUE;
    vp8_dct_value_cost_ptr   = dct_value_cost + DCT_MAX_VALUE;
}

static void tokenize2nd_order_b
(
    const BLOCKD *const b,
    TOKENEXTRA **tp,
    const int type,     /* which plane: 0=Y no DC, 1=Y2, 2=UV, 3=Y with DC */
    const FRAME_TYPE frametype,
    ENTROPY_CONTEXT *a,
    ENTROPY_CONTEXT *l,
    VP8_COMP *cpi
)
{
    int pt; /* near block/prev token context index */
    int c = 0;          /* start at DC */
    const int eob = b->eob;     /* one beyond last nonzero coeff */
    TOKENEXTRA *t = *tp;        /* store tokens starting here */
    int x;
    const short *qcoeff_ptr = b->qcoeff;
    VP8_COMBINEENTROPYCONTEXTS(pt, *a, *l);

    do
    {
        const int band = vp8_coef_bands[c];

        if (c < eob)
        {
            int rc = vp8_default_zig_zag1d[c];
            const int v = qcoeff_ptr[rc];

            assert(-DCT_MAX_VALUE <= v  &&  v < (DCT_MAX_VALUE));

            t->Extra = vp8_dct_value_tokens_ptr[v].Extra;
            x        = vp8_dct_value_tokens_ptr[v].Token;
        }
        else
            x = DCT_EOB_TOKEN;

        t->Token = x;
        t->context_tree = cpi->common.fc.coef_probs [type] [band] [pt];

        t->skip_eob_node = pt == 0 && ((band > 0 && type > 0) || (band > 1 && type == 0));

        ++cpi->coef_counts       [type] [band] [pt] [x];
    }
    while (pt = vp8_prev_token_class[x], ++t, c < eob  &&  ++c < 16);

    *tp = t;
    pt = (c != !type); /* 0 <-> all coeff data is zero */
    *a = *l = pt;

}

static void tokenize1st_order_b
(
    const BLOCKD *const b,
    TOKENEXTRA **tp,
    const int type,     /* which plane: 0=Y no DC, 1=Y2, 2=UV, 3=Y with DC */
    const FRAME_TYPE frametype,
    ENTROPY_CONTEXT *a,
    ENTROPY_CONTEXT *l,
    VP8_COMP *cpi
)
{
    int pt; /* near block/prev token context index */
    int c = type ? 0 : 1;       /* start at DC unless type 0 */
    const int eob = b->eob;     /* one beyond last nonzero coeff */
    TOKENEXTRA *t = *tp;        /* store tokens starting here */
    int x;
    const short *qcoeff_ptr = b->qcoeff;
    VP8_COMBINEENTROPYCONTEXTS(pt, *a, *l);

    do
    {
        const int band = vp8_coef_bands[c];

        x = DCT_EOB_TOKEN;

        if (c < eob)
        {
            int rc = vp8_default_zig_zag1d[c];
            const int v = qcoeff_ptr[rc];

            assert(-DCT_MAX_VALUE <= v  &&  v < (DCT_MAX_VALUE));

            t->Extra = vp8_dct_value_tokens_ptr[v].Extra;
            x        = vp8_dct_value_tokens_ptr[v].Token;
        }

        t->Token = x;
        t->context_tree = cpi->common.fc.coef_probs [type] [band] [pt];

        t->skip_eob_node = pt == 0 && ((band > 0 && type > 0) || (band > 1 && type == 0));

        ++cpi->coef_counts       [type] [band] [pt] [x];
    }
    while (pt = vp8_prev_token_class[x], ++t, c < eob  &&  ++c < 16);

    *tp = t;
    pt = (c != !type); /* 0 <-> all coeff data is zero */
    *a = *l = pt;

}


static int mb_is_skippable(MACROBLOCKD *x)
{
    int has_y2_block;
    int skip = 1;
    int i = 0;

    has_y2_block = (x->mode_info_context->mbmi.mode != B_PRED
                    && x->mode_info_context->mbmi.mode != SPLITMV);
    if (has_y2_block)
    {
        for (i = 0; i < 16; i++)
            skip &= (x->block[i].eob < 2);
    }

    for (; i < 24 + has_y2_block; i++)
        skip &= (!x->block[i].eob);

    return skip;
}


void vp8_tokenize_mb(VP8_COMP *cpi, MACROBLOCKD *x, TOKENEXTRA **t)
{
    ENTROPY_CONTEXT * A = (ENTROPY_CONTEXT *)x->above_context;
    ENTROPY_CONTEXT * L = (ENTROPY_CONTEXT *)x->left_context;
    int plane_type;
    int b;

    TOKENEXTRA *start = *t;
    TOKENEXTRA *tp = *t;

    x->mode_info_context->mbmi.dc_diff = 1;

#if 0

    if (x->mbmi.force_no_skip)
    {
        x->mbmi.mb_skip_coeff = 1;
        //reset for next_mb.
        x->mbmi.force_no_skip = 0;
    }

#endif

#if 1

    x->mode_info_context->mbmi.mb_skip_coeff = mb_is_skippable(x);
    if (x->mode_info_context->mbmi.mb_skip_coeff)
    {

        cpi->skip_true_count++;

        if (!cpi->common.mb_no_coeff_skip)
            vp8_stuff_mb(cpi, x, t) ;
        else
        {
            vp8_fix_contexts(x);
        }

        if (x->mode_info_context->mbmi.mode != B_PRED && x->mode_info_context->mbmi.mode != SPLITMV)
            x->mode_info_context->mbmi.dc_diff = 0;
        else
            x->mode_info_context->mbmi.dc_diff = 1;


        return;
    }

    cpi->skip_false_count++;
#endif
#if 0
    vpx_memcpy(cpi->coef_counts_backup, cpi->coef_counts, sizeof(cpi->coef_counts));
#endif

    if (x->mode_info_context->mbmi.mode == B_PRED || x->mode_info_context->mbmi.mode == SPLITMV)
    {
        plane_type = 3;
    }
    else
    {
        tokenize2nd_order_b(x->block + 24, t, 1, x->frame_type,
                   A + vp8_block2above[24], L + vp8_block2left[24], cpi);
        plane_type = 0;

    }

    for (b = 0; b < 16; b++)
        tokenize1st_order_b(x->block + b, t, plane_type, x->frame_type,
                            A + vp8_block2above[b],
                            L + vp8_block2left[b], cpi);

    for (b = 16; b < 24; b++)
        tokenize1st_order_b(x->block + b, t, 2, x->frame_type,
                            A + vp8_block2above[b],
                            L + vp8_block2left[b], cpi);

#if 0

    if (cpi->common.mb_no_coeff_skip)
    {
        int skip = 1;

        while ((tp != *t) && skip)
        {
            skip = (skip && (tp->Token == DCT_EOB_TOKEN));
            tp ++;
        }

        if (skip != x->mbmi.mb_skip_coeff)
            skip += 0;

        x->mbmi.mb_skip_coeff = skip;

        if (x->mbmi.mb_skip_coeff == 1)
        {
            x->mbmi.dc_diff = 0;
            //redo the coutnts
            vpx_memcpy(cpi->coef_counts, cpi->coef_counts_backup, sizeof(cpi->coef_counts));

            *t = start;
            cpi->skip_true_count++;
            //skip_true_count++;
        }
        else
        {

            cpi->skip_false_count++;
            //skip_false_count++;
        }
    }

#endif
}


#ifdef ENTROPY_STATS

void init_context_counters(void)
{
    vpx_memset(context_counters, 0, sizeof(context_counters));
}

void print_context_counters()
{

    int type, band, pt, t;

    FILE *const f = fopen("context.c", "w");

    fprintf(f, "#include \"entropy.h\"\n");

    fprintf(f, "\n/* *** GENERATED FILE: DO NOT EDIT *** */\n\n");

    fprintf(f, "int Contexts[BLOCK_TYPES] [COEF_BANDS] [PREV_COEF_CONTEXTS] [vp8_coef_tokens];\n\n");

    fprintf(f, "const int default_contexts[BLOCK_TYPES] [COEF_BANDS] [PREV_COEF_CONTEXTS] [vp8_coef_tokens] = {");

# define Comma( X) (X? ",":"")

    type = 0;

    do
    {
        fprintf(f, "%s\n  { /* block Type %d */", Comma(type), type);

        band = 0;

        do
        {
            fprintf(f, "%s\n    { /* Coeff Band %d */", Comma(band), band);

            pt = 0;

            do
            {
                fprintf(f, "%s\n      {", Comma(pt));

                t = 0;

                do
                {
                    const _int64 x = context_counters [type] [band] [pt] [t];
                    const int y = (int) x;

                    assert(x == (_int64) y);  /* no overflow handling yet */
                    fprintf(f, "%s %d", Comma(t), y);

                }
                while (++t < vp8_coef_tokens);

                fprintf(f, "}");
            }
            while (++pt < PREV_COEF_CONTEXTS);

            fprintf(f, "\n    }");

        }
        while (++band < COEF_BANDS);

        fprintf(f, "\n  }");
    }
    while (++type < BLOCK_TYPES);

    fprintf(f, "\n};\n");
    fclose(f);
}
#endif


void vp8_tokenize_initialize()
{
    fill_value_tokens();
}


static __inline void stuff2nd_order_b
(
    const BLOCKD *const b,
    TOKENEXTRA **tp,
    const int type,     /* which plane: 0=Y no DC, 1=Y2, 2=UV, 3=Y with DC */
    const FRAME_TYPE frametype,
    ENTROPY_CONTEXT *a,
    ENTROPY_CONTEXT *l,
    VP8_COMP *cpi
)
{
    int pt; /* near block/prev token context index */
    TOKENEXTRA *t = *tp;        /* store tokens starting here */
    VP8_COMBINEENTROPYCONTEXTS(pt, *a, *l);
    (void) frametype;
    (void) type;
    (void) b;

    t->Token = DCT_EOB_TOKEN;
    t->context_tree = cpi->common.fc.coef_probs [1] [0] [pt];
    t->skip_eob_node = 0;
    ++cpi->coef_counts       [1] [0] [pt] [DCT_EOB_TOKEN];
    ++t;

    *tp = t;
    pt = 0;
    *a = *l = pt;

}

static __inline void stuff1st_order_b
(
    const BLOCKD *const b,
    TOKENEXTRA **tp,
    const int type,     /* which plane: 0=Y no DC, 1=Y2, 2=UV, 3=Y with DC */
    const FRAME_TYPE frametype,
    ENTROPY_CONTEXT *a,
    ENTROPY_CONTEXT *l,
    VP8_COMP *cpi
)
{
    int pt; /* near block/prev token context index */
    TOKENEXTRA *t = *tp;        /* store tokens starting here */
    VP8_COMBINEENTROPYCONTEXTS(pt, *a, *l);
    (void) frametype;
    (void) type;
    (void) b;

    t->Token = DCT_EOB_TOKEN;
    t->context_tree = cpi->common.fc.coef_probs [0] [1] [pt];
    t->skip_eob_node = 0;
    ++cpi->coef_counts       [0] [1] [pt] [DCT_EOB_TOKEN];
    ++t;
    *tp = t;
    pt = 0; /* 0 <-> all coeff data is zero */
    *a = *l = pt;

}
static __inline
void stuff1st_order_buv
(
    const BLOCKD *const b,
    TOKENEXTRA **tp,
    const int type,     /* which plane: 0=Y no DC, 1=Y2, 2=UV, 3=Y with DC */
    const FRAME_TYPE frametype,
    ENTROPY_CONTEXT *a,
    ENTROPY_CONTEXT *l,
    VP8_COMP *cpi
)
{
    int pt; /* near block/prev token context index */
    TOKENEXTRA *t = *tp;        /* store tokens starting here */
    VP8_COMBINEENTROPYCONTEXTS(pt, *a, *l);
    (void) frametype;
    (void) type;
    (void) b;

    t->Token = DCT_EOB_TOKEN;
    t->context_tree = cpi->common.fc.coef_probs [2] [0] [pt];
    t->skip_eob_node = 0;
    ++cpi->coef_counts[2] [0] [pt] [DCT_EOB_TOKEN];
    ++t;
    *tp = t;
    pt = 0; /* 0 <-> all coeff data is zero */
    *a = *l = pt;

}

void vp8_stuff_mb(VP8_COMP *cpi, MACROBLOCKD *x, TOKENEXTRA **t)
{
    ENTROPY_CONTEXT * A = (ENTROPY_CONTEXT *)x->above_context;
    ENTROPY_CONTEXT * L = (ENTROPY_CONTEXT *)x->left_context;
    int plane_type;
    int b;

    stuff2nd_order_b(x->block + 24, t, 1, x->frame_type,
                     A + vp8_block2above[24], L + vp8_block2left[24], cpi);
    plane_type = 0;


    if (x->mode_info_context->mbmi.mode != B_PRED && x->mode_info_context->mbmi.mode != SPLITMV)
        x->mode_info_context->mbmi.dc_diff = 0;
    else
        x->mode_info_context->mbmi.dc_diff = 1;


    for (b = 0; b < 16; b++)
        stuff1st_order_b(x->block + b, t, plane_type, x->frame_type,
                         A + vp8_block2above[b],
                         L + vp8_block2left[b], cpi);

    for (b = 16; b < 24; b++)
        stuff1st_order_buv(x->block + b, t, 2, x->frame_type,
                           A + vp8_block2above[b],
                           L + vp8_block2left[b], cpi);

}
void vp8_fix_contexts(MACROBLOCKD *x)
{
    /* Clear entropy contexts for Y2 blocks */
    if (x->mode_info_context->mbmi.mode != B_PRED && x->mode_info_context->mbmi.mode != SPLITMV)
    {
        vpx_memset(x->above_context, 0, sizeof(ENTROPY_CONTEXT_PLANES));
        vpx_memset(x->left_context, 0, sizeof(ENTROPY_CONTEXT_PLANES));
    }
    else
    {
        vpx_memset(x->above_context, 0, sizeof(ENTROPY_CONTEXT_PLANES)-1);
        vpx_memset(x->left_context, 0, sizeof(ENTROPY_CONTEXT_PLANES)-1);
    }

}
