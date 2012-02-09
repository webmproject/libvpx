/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */


#include "treereader.h"
#include "vp8/common/entropymv.h"
#include "vp8/common/entropymode.h"
#include "onyxd_int.h"
#include "vp8/common/findnearmv.h"

#include "vp8/common/seg_common.h"
#include "vp8/common/pred_common.h"

#if CONFIG_DEBUG
#include <assert.h>
#endif
static int vp8_read_bmode(vp8_reader *bc, const vp8_prob *p)
{
    const int i = vp8_treed_read(bc, vp8_bmode_tree, p);

    return i;
}


static int vp8_read_ymode(vp8_reader *bc, const vp8_prob *p)
{
    const int i = vp8_treed_read(bc, vp8_ymode_tree, p);

    return i;
}

static int vp8_kfread_ymode(vp8_reader *bc, const vp8_prob *p)
{
    const int i = vp8_treed_read(bc, vp8_kf_ymode_tree, p);

    return i;
}
static int vp8_read_i8x8_mode(vp8_reader *bc, const vp8_prob *p)
{
    const int i = vp8_treed_read(bc, vp8_i8x8_mode_tree, p);

    return i;
}


static int vp8_read_uv_mode(vp8_reader *bc, const vp8_prob *p)
{
    const int i = vp8_treed_read(bc, vp8_uv_mode_tree, p);

    return i;
}

// This function reads the current macro block's segnent id from the bitstream
// It should only be called if a segment map update is indicated.
static void vp8_read_mb_segid(vp8_reader *r, MB_MODE_INFO *mi, MACROBLOCKD *x)
{
    /* Is segmentation enabled */
    if (x->segmentation_enabled && x->update_mb_segmentation_map)
    {
        /* If so then read the segment id. */
        if (vp8_read(r, x->mb_segment_tree_probs[0]))
            mi->segment_id = (unsigned char)(2 + vp8_read(r, x->mb_segment_tree_probs[2]));
        else
            mi->segment_id = (unsigned char)(vp8_read(r, x->mb_segment_tree_probs[1]));
    }
}
extern const int vp8_i8x8_block[4];
static void vp8_kfread_modes(VP8D_COMP *pbi,
                             MODE_INFO *m,
                             int mb_row,
                             int mb_col)
{
    vp8_reader *const bc = & pbi->bc;
    const int mis = pbi->common.mode_info_stride;
    int map_index = mb_row * pbi->common.mb_cols + mb_col;
    MB_PREDICTION_MODE y_mode;

    // Read the Macroblock segmentation map if it is being updated explicitly
    // this frame (reset to 0 by default).
    m->mbmi.segment_id = 0;
    if (pbi->mb.update_mb_segmentation_map)
    {
        vp8_read_mb_segid(bc, &m->mbmi, &pbi->mb);
        pbi->common.last_frame_seg_map[map_index] = m->mbmi.segment_id;
    }

    if ( pbi->common.mb_no_coeff_skip &&
         ( !segfeature_active( &pbi->mb,
                               m->mbmi.segment_id, SEG_LVL_EOB ) ||
           ( get_segdata( &pbi->mb,
                          m->mbmi.segment_id, SEG_LVL_EOB ) != 0 ) ) )
    {
        m->mbmi.mb_skip_coeff = vp8_read(bc, pbi->prob_skip_false);
    }
    else
    {
        if ( segfeature_active( &pbi->mb,
                                m->mbmi.segment_id, SEG_LVL_EOB ) &&
             ( get_segdata( &pbi->mb,
                            m->mbmi.segment_id, SEG_LVL_EOB ) == 0 ) )
        {
            m->mbmi.mb_skip_coeff = 1;
        }
        else
            m->mbmi.mb_skip_coeff = 0;
    }

#if CONFIG_QIMODE
    y_mode = (MB_PREDICTION_MODE) vp8_kfread_ymode(bc,
        pbi->common.kf_ymode_prob[pbi->common.kf_ymode_probs_index]);
#else
    y_mode = (MB_PREDICTION_MODE) vp8_kfread_ymode(
                                      bc, pbi->common.kf_ymode_prob);
#endif

    m->mbmi.ref_frame = INTRA_FRAME;

    if ((m->mbmi.mode = y_mode) == B_PRED)
    {
        int i = 0;
        do
        {
            const B_PREDICTION_MODE A = above_block_mode(m, i, mis);
            const B_PREDICTION_MODE L = left_block_mode(m, i);

            m->bmi[i].as_mode =
                (B_PREDICTION_MODE) vp8_read_bmode(
                                        bc, pbi->common.kf_bmode_prob [A] [L]);
        }
        while (++i < 16);
    }
    if((m->mbmi.mode = y_mode) == I8X8_PRED)
    {
        int i;
        int mode8x8;
        for(i=0;i<4;i++)
         {
             int ib = vp8_i8x8_block[i];
             mode8x8 = vp8_read_i8x8_mode(bc, pbi->common.i8x8_mode_prob);
             m->bmi[ib+0].as_mode= mode8x8;
             m->bmi[ib+1].as_mode= mode8x8;
             m->bmi[ib+4].as_mode= mode8x8;
             m->bmi[ib+5].as_mode= mode8x8;
         }
   }
    else
#if CONFIG_UVINTRA
        m->mbmi.uv_mode = (MB_PREDICTION_MODE)vp8_read_uv_mode(bc,
            pbi->common.kf_uv_mode_prob[m->mbmi.mode]);
#else
        m->mbmi.uv_mode = (MB_PREDICTION_MODE)vp8_read_uv_mode(bc,
            pbi->common.kf_uv_mode_prob);
#endif
}

static int read_mvcomponent(vp8_reader *r, const MV_CONTEXT *mvc)
{
    const vp8_prob *const p = (const vp8_prob *) mvc;
    int x = 0;

    if (vp8_read(r, p [mvpis_short]))  /* Large */
    {
        int i = 0;

        do
        {
            x += vp8_read(r, p [MVPbits + i]) << i;
        }
        while (++i < 3);

        i = mvlong_width - 1;  /* Skip bit 3, which is sometimes implicit */

        do
        {
            x += vp8_read(r, p [MVPbits + i]) << i;
        }
        while (--i > 3);

        if (!(x & 0xFFF0)  ||  vp8_read(r, p [MVPbits + 3]))
            x += 8;
    }
    else   /* small */
        x = vp8_treed_read(r, vp8_small_mvtree, p + MVPshort);

    if (x  &&  vp8_read(r, p [MVPsign]))
        x = -x;

    return x;
}

static void read_mv(vp8_reader *r, MV *mv, const MV_CONTEXT *mvc)
{
    mv->row = (short)(read_mvcomponent(r,   mvc) << 1);
    mv->col = (short)(read_mvcomponent(r, ++mvc) << 1);
}


static void read_mvcontexts(vp8_reader *bc, MV_CONTEXT *mvc)
{
    int i = 0;

    do
    {
        const vp8_prob *up = vp8_mv_update_probs[i].prob;
        vp8_prob *p = (vp8_prob *)(mvc + i);
        vp8_prob *const pstop = p + MVPcount;

        do
        {
            if (vp8_read(bc, *up++))
            {
                const vp8_prob x = (vp8_prob)vp8_read_literal(bc, 7);

                *p = x ? x << 1 : 1;
            }
        }
        while (++p < pstop);
    }
    while (++i < 2);
}

// Read the referncence frame
static MV_REFERENCE_FRAME read_ref_frame( VP8D_COMP *pbi,
                                          vp8_reader *const bc,
                                          unsigned char segment_id )
{
    MV_REFERENCE_FRAME ref_frame;
    int seg_ref_active;
    int seg_ref_count = 0;

    VP8_COMMON *const cm = & pbi->common;
    MACROBLOCKD *const xd = &pbi->mb;

    seg_ref_active = segfeature_active( xd,
                                        segment_id,
                                        SEG_LVL_REF_FRAME );

    // If segment coding enabled does the segment allow for more than one
    // possible reference frame
    if ( seg_ref_active )
    {
        seg_ref_count = check_segref( xd, segment_id, INTRA_FRAME ) +
                        check_segref( xd, segment_id, LAST_FRAME ) +
                        check_segref( xd, segment_id, GOLDEN_FRAME ) +
                        check_segref( xd, segment_id, ALTREF_FRAME );
    }

    // Segment reference frame features not available or allows for
    // multiple reference frame options
    if ( !seg_ref_active || (seg_ref_count > 1) )
    {
        // Values used in prediction model coding
        unsigned char prediction_flag;
        vp8_prob pred_prob;
        MV_REFERENCE_FRAME pred_ref;

        // Get the context probability the prediction flag
        pred_prob = get_pred_prob( cm, xd, PRED_REF );

        // Read the prediction status flag
        prediction_flag = (unsigned char)vp8_read( bc, pred_prob );

        // Store the prediction flag.
        set_pred_flag( xd, PRED_REF, prediction_flag );

        // Get the predicted reference frame.
        pred_ref = get_pred_ref( cm, xd );

        // If correctly predicted then use the predicted value
        if ( prediction_flag )
        {
            ref_frame = pred_ref;
        }
        // else decode the explicitly coded value
        else
        {
            vp8_prob mod_refprobs[PREDICTION_PROBS];
            vpx_memcpy( mod_refprobs,
                        cm->mod_refprobs[pred_ref], sizeof(mod_refprobs) );

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

            // Default to INTRA_FRAME (value 0)
            ref_frame = INTRA_FRAME;

            // Do we need to decode the Intra/Inter branch
            if ( mod_refprobs[0] )
                ref_frame = (MV_REFERENCE_FRAME) vp8_read(bc, mod_refprobs[0]);
            else
                ref_frame ++;

            if (ref_frame)
            {
                // Do we need to decode the Last/Gf_Arf branch
                if ( mod_refprobs[1] )
                    ref_frame += vp8_read(bc, mod_refprobs[1]);
                else
                    ref_frame++;

                if ( ref_frame > 1 )
                {
                    // Do we need to decode the GF/Arf branch
                    if ( mod_refprobs[2] )
                        ref_frame += vp8_read(bc, mod_refprobs[2]);
                    else
                    {
                        if ( seg_ref_active )
                        {
                            if ( (pred_ref == GOLDEN_FRAME) ||
                                 !check_segref( xd, segment_id, GOLDEN_FRAME) )
                            {
                                ref_frame = ALTREF_FRAME;
                            }
                            else
                                ref_frame = GOLDEN_FRAME;
                        }
                        else
                            ref_frame = (pred_ref == GOLDEN_FRAME)
                                        ? ALTREF_FRAME : GOLDEN_FRAME;
                    }
                }
            }
        }
    }

    // Segment reference frame features are enabled
    else
    {
        // The reference frame for the mb is considered as correclty predicted
        // if it is signaled at the segment level for the purposes of the
        // common prediction model
        set_pred_flag( xd, PRED_REF, 1 );
        ref_frame = get_pred_ref( cm, xd );
    }

    return (MV_REFERENCE_FRAME)ref_frame;
}

static MB_PREDICTION_MODE read_mv_ref(vp8_reader *bc, const vp8_prob *p)
{
    const int i = vp8_treed_read(bc, vp8_mv_ref_tree, p);

    return (MB_PREDICTION_MODE)i;
}

static B_PREDICTION_MODE sub_mv_ref(vp8_reader *bc, const vp8_prob *p)
{
    const int i = vp8_treed_read(bc, vp8_sub_mv_ref_tree, p);

    return (B_PREDICTION_MODE)i;
}

#ifdef VPX_MODE_COUNT
unsigned int vp8_mv_cont_count[5][4] =
{
    { 0, 0, 0, 0 },
    { 0, 0, 0, 0 },
    { 0, 0, 0, 0 },
    { 0, 0, 0, 0 },
    { 0, 0, 0, 0 }
};
#endif

static const unsigned char mbsplit_fill_count[4] = {8, 8, 4, 1};
static const unsigned char mbsplit_fill_offset[4][16] = {
    { 0,  1,  2,  3,  4,  5,  6,  7,  8,  9,  10, 11, 12, 13, 14, 15},
    { 0,  1,  4,  5,  8,  9, 12, 13,  2,  3,   6,  7, 10, 11, 14, 15},
    { 0,  1,  4,  5,  2,  3,  6,  7,  8,  9,  12, 13, 10, 11, 14, 15},
    { 0,  1,  2,  3,  4,  5,  6,  7,  8,  9,  10, 11, 12, 13, 14, 15}
};



static void mb_mode_mv_init(VP8D_COMP *pbi)
{
    VP8_COMMON *const cm = & pbi->common;
    vp8_reader *const bc = & pbi->bc;
    MV_CONTEXT *const mvc = pbi->common.fc.mvc;

#if CONFIG_ERROR_CONCEALMENT
    /* Default is that no macroblock is corrupt, therefore we initialize
     * mvs_corrupt_from_mb to something very big, which we can be sure is
     * outside the frame. */
    pbi->mvs_corrupt_from_mb = UINT_MAX;
#endif
    pbi->prob_skip_false = 0;
    if (pbi->common.mb_no_coeff_skip)
        pbi->prob_skip_false = (vp8_prob)vp8_read_literal(bc, 8);

    if(pbi->common.frame_type != KEY_FRAME)
    {
        // Decode the baseline probabilities for decoding reference frame
        cm->prob_intra_coded = (vp8_prob)vp8_read_literal(bc, 8);
        cm->prob_last_coded  = (vp8_prob)vp8_read_literal(bc, 8);
        cm->prob_gf_coded    = (vp8_prob)vp8_read_literal(bc, 8);

        // Computes a modified set of probabilities for use when reference
        // frame prediction fails.
        compute_mod_refprobs( cm );

#if CONFIG_DUALPRED
        pbi->common.dual_pred_mode = vp8_read(bc, 128);
        if (cm->dual_pred_mode)
            cm->dual_pred_mode += vp8_read(bc, 128);
        if (cm->dual_pred_mode == HYBRID_PREDICTION)
        {
            int i;
            for ( i = 0; i < DUAL_PRED_CONTEXTS; i++ )
                cm->prob_dualpred[i] = (vp8_prob)vp8_read_literal(bc, 8);
        }
#endif /* CONFIG_DUALPRED */

        if (vp8_read_bit(bc))
        {
            int i = 0;

            do
            {
                cm->fc.ymode_prob[i] = (vp8_prob) vp8_read_literal(bc, 8);
            }
            while (++i < VP8_YMODES-1);
        }
#if CONFIG_UVINTRA
        //vp8_read_bit(bc);
#else
        if (vp8_read_bit(bc))
        {
            int i = 0;

            do
            {
                cm->fc.uv_mode_prob[i] = (vp8_prob) vp8_read_literal(bc, 8);
            }
            while (++i < VP8_UV_MODES-1);
        }
#endif /* CONFIG_UVINTRA */
        read_mvcontexts(bc, mvc);
    }
}

// This function either reads the segment id for the current macroblock from
// the bitstream or if the value is temporally predicted asserts the predicted
// value
static void read_mb_segment_id ( VP8D_COMP *pbi,
                                 int mb_row, int mb_col )
{
    vp8_reader *const bc = & pbi->bc;
    VP8_COMMON *const cm = & pbi->common;
    MACROBLOCKD *const xd  = & pbi->mb;
    MODE_INFO *mi = xd->mode_info_context;
    MB_MODE_INFO *mbmi = &mi->mbmi;
    int index = mb_row * pbi->common.mb_cols + mb_col;

    if (xd->segmentation_enabled)
    {
        if (xd->update_mb_segmentation_map)
        {
            // Is temporal coding of the segment id for this mb enabled.
            if (cm->temporal_update)
            {
                // Get the context based probability for reading the
                // prediction status flag
                vp8_prob pred_prob =
                    get_pred_prob( cm, xd, PRED_SEG_ID );

                // Read the prediction status flag
                unsigned char seg_pred_flag =
                    (unsigned char)vp8_read(bc, pred_prob );

                // Store the prediction flag.
                set_pred_flag( xd, PRED_SEG_ID, seg_pred_flag );

                // If the value is flagged as correctly predicted
                // then use the predicted value
                if ( seg_pred_flag )
                {
                    mbmi->segment_id = get_pred_mb_segid( cm, index );
                }
                // Else .... decode it explicitly
                else
                {
                    vp8_read_mb_segid(bc, mbmi, xd );
                    cm->last_frame_seg_map[index] = mbmi->segment_id;
                }

            }
            // Normal unpredicted coding mode
            else
            {
                vp8_read_mb_segid(bc, mbmi, xd);
                cm->last_frame_seg_map[index] = mbmi->segment_id;
            }
        }
    }
    else
    {
        // The encoder explicitly sets the segment_id to 0
        // when segmentation is disabled
        mbmi->segment_id = 0;
    }
}

static void read_mb_modes_mv(VP8D_COMP *pbi, MODE_INFO *mi, MB_MODE_INFO *mbmi,
                             MODE_INFO *prev_mi,
                             int mb_row, int mb_col)
{
    VP8_COMMON *const cm = & pbi->common;
    vp8_reader *const bc = & pbi->bc;
    MV_CONTEXT *const mvc = pbi->common.fc.mvc;
    const int mis = pbi->common.mode_info_stride;
    MACROBLOCKD *const xd  = & pbi->mb;

    int pred_context;
    int index = mb_row * pbi->common.mb_cols + mb_col;
    int_mv *const mv = & mbmi->mv;
    int mb_to_left_edge;
    int mb_to_right_edge;
    int mb_to_top_edge;
    int mb_to_bottom_edge;

    mb_to_top_edge = xd->mb_to_top_edge;
    mb_to_bottom_edge = xd->mb_to_bottom_edge;
    mb_to_top_edge -= LEFT_TOP_MARGIN;
    mb_to_bottom_edge += RIGHT_BOTTOM_MARGIN;
    mbmi->need_to_clamp_mvs = 0;
#if CONFIG_DUALPRED
    mbmi->second_ref_frame = 0;
#endif /* CONFIG_DUALPRED */
    /* Distance of Mb to the various image edges.
     * These specified to 8th pel as they are always compared to MV values that are in 1/8th pel units
     */
    xd->mb_to_left_edge =
    mb_to_left_edge = -((mb_col * 16) << 3);
    mb_to_left_edge -= LEFT_TOP_MARGIN;

    xd->mb_to_right_edge =
    mb_to_right_edge = ((pbi->common.mb_cols - 1 - mb_col) * 16) << 3;
    mb_to_right_edge += RIGHT_BOTTOM_MARGIN;

    // Make sure the MACROBLOCKD mode info pointer is pointed at the
    // correct entry for the current macroblock.
    xd->mode_info_context = mi;

    // Read the macroblock segment id.
    read_mb_segment_id ( pbi, mb_row, mb_col );

    if ( pbi->common.mb_no_coeff_skip &&
         ( !segfeature_active( xd,
                               mbmi->segment_id, SEG_LVL_EOB ) ||
           (get_segdata( xd, mbmi->segment_id, SEG_LVL_EOB ) != 0) ) )
    {
        // Read the macroblock coeff skip flag if this feature is in use,
        // else default to 0
        mbmi->mb_skip_coeff = vp8_read(bc, pbi->prob_skip_false);
    }
    else
    {
        if ( segfeature_active( xd,
                                mbmi->segment_id, SEG_LVL_EOB ) &&
             (get_segdata( xd, mbmi->segment_id, SEG_LVL_EOB ) == 0) )
        {
            mbmi->mb_skip_coeff = 1;
        }
        else
            mbmi->mb_skip_coeff = 0;
    }

    // Read the reference frame
    mbmi->ref_frame = read_ref_frame( pbi, bc, mbmi->segment_id );

    // If reference frame is an Inter frame
    if (mbmi->ref_frame)
    {
        int rct[4];
        int_mv nearest, nearby, best_mv;
        vp8_prob mv_ref_p [VP8_MVREFS-1];

        vp8_find_near_mvs(xd, mi,
            prev_mi,
            &nearest, &nearby, &best_mv, rct,
                          mbmi->ref_frame, pbi->common.ref_frame_sign_bias);
        vp8_mv_ref_probs(&pbi->common, mv_ref_p, rct);

        // Is the segment level mode feature enabled for this segment
        if ( segfeature_active( xd, mbmi->segment_id, SEG_LVL_MODE ) )
        {
            mbmi->mode =
                get_segdata( xd, mbmi->segment_id, SEG_LVL_MODE );
        }
        else
        {
            mbmi->mode = read_mv_ref(bc, mv_ref_p);

            vp8_accum_mv_refs(&pbi->common, mbmi->mode, rct);
        }

        mbmi->uv_mode = DC_PRED;
        switch (mbmi->mode)
        {
        case SPLITMV:
        {
            const int s = mbmi->partitioning =
                      vp8_treed_read(bc, vp8_mbsplit_tree, vp8_mbsplit_probs);
            const int num_p = vp8_mbsplit_count [s];
            int j = 0;

            mbmi->need_to_clamp_mvs = 0;
            do  /* for each subset j */
            {
                int_mv leftmv, abovemv;
                int_mv blockmv;
                int k;  /* first block in subset j */
                int mv_contz;
                k = vp8_mbsplit_offset[s][j];

                leftmv.as_int = left_block_mv(mi, k);
                abovemv.as_int = above_block_mv(mi, k, mis);
                mv_contz = vp8_mv_cont(&leftmv, &abovemv);

                switch (sub_mv_ref(bc, vp8_sub_mv_ref_prob2 [mv_contz])) /*pc->fc.sub_mv_ref_prob))*/
                {
                case NEW4X4:
                    read_mv(bc, &blockmv.as_mv, (const MV_CONTEXT *) mvc);
                    blockmv.as_mv.row += best_mv.as_mv.row;
                    blockmv.as_mv.col += best_mv.as_mv.col;
  #ifdef VPX_MODE_COUNT
                    vp8_mv_cont_count[mv_contz][3]++;
  #endif
                    break;
                case LEFT4X4:
                    blockmv.as_int = leftmv.as_int;
  #ifdef VPX_MODE_COUNT
                    vp8_mv_cont_count[mv_contz][0]++;
  #endif
                    break;
                case ABOVE4X4:
                    blockmv.as_int = abovemv.as_int;
  #ifdef VPX_MODE_COUNT
                    vp8_mv_cont_count[mv_contz][1]++;
  #endif
                    break;
                case ZERO4X4:
                    blockmv.as_int = 0;
  #ifdef VPX_MODE_COUNT
                    vp8_mv_cont_count[mv_contz][2]++;
  #endif
                    break;
                default:
                    break;
                }

                mbmi->need_to_clamp_mvs |= vp8_check_mv_bounds(&blockmv,
                                                          mb_to_left_edge,
                                                          mb_to_right_edge,
                                                          mb_to_top_edge,
                                                          mb_to_bottom_edge);

                {
                    /* Fill (uniform) modes, mvs of jth subset.
                     Must do it here because ensuing subsets can
                     refer back to us via "left" or "above". */
                    const unsigned char *fill_offset;
                    unsigned int fill_count = mbsplit_fill_count[s];

                    fill_offset = &mbsplit_fill_offset[s][(unsigned char)j * mbsplit_fill_count[s]];

                    do {
                        mi->bmi[ *fill_offset].mv.as_int = blockmv.as_int;
                        fill_offset++;
                    }while (--fill_count);
                }

            }
            while (++j < num_p);
        }

        mv->as_int = mi->bmi[15].mv.as_int;

        break;  /* done with SPLITMV */

        case NEARMV:
            mv->as_int = nearby.as_int;
            /* Clip "next_nearest" so that it does not extend to far out of image */
            vp8_clamp_mv(mv, mb_to_left_edge, mb_to_right_edge,
                         mb_to_top_edge, mb_to_bottom_edge);
            goto propagate_mv;

        case NEARESTMV:
            mv->as_int = nearest.as_int;
            /* Clip "next_nearest" so that it does not extend to far out of image */
            vp8_clamp_mv(mv, mb_to_left_edge, mb_to_right_edge,
                         mb_to_top_edge, mb_to_bottom_edge);
            goto propagate_mv;

        case ZEROMV:
            mv->as_int = 0;
            goto propagate_mv;

        case NEWMV:
            read_mv(bc, &mv->as_mv, (const MV_CONTEXT *) mvc);
            mv->as_mv.row += best_mv.as_mv.row;
            mv->as_mv.col += best_mv.as_mv.col;

            /* Don't need to check this on NEARMV and NEARESTMV modes
             * since those modes clamp the MV. The NEWMV mode does not,
             * so signal to the prediction stage whether special
             * handling may be required.
             */
            mbmi->need_to_clamp_mvs = vp8_check_mv_bounds(mv,
                                                      mb_to_left_edge,
                                                      mb_to_right_edge,
                                                      mb_to_top_edge,
                                                      mb_to_bottom_edge);

        propagate_mv:  /* same MV throughout */

#if CONFIG_DUALPRED
            if ( cm->dual_pred_mode == DUAL_PREDICTION_ONLY ||
                 (cm->dual_pred_mode == HYBRID_PREDICTION &&
                     vp8_read(bc, get_pred_prob( cm, xd, PRED_DUAL ))) )
            {
                mbmi->second_ref_frame = mbmi->ref_frame + 1;
                if (mbmi->second_ref_frame == 4)
                    mbmi->second_ref_frame = 1;
            }
            if (mbmi->second_ref_frame)
            {
                vp8_find_near_mvs(xd, mi,
                                  prev_mi,
                                  &nearest, &nearby, &best_mv, rct,
                                  (int)mbmi->second_ref_frame,
                                  pbi->common.ref_frame_sign_bias);
                switch (mbmi->mode) {
                case ZEROMV:
                    mbmi->second_mv.as_int = 0;
                    break;
                case NEARMV:
                    mbmi->second_mv.as_int = nearby.as_int;
                    vp8_clamp_mv(&mbmi->second_mv, mb_to_left_edge, mb_to_right_edge,
                                 mb_to_top_edge, mb_to_bottom_edge);
                    break;
                case NEARESTMV:
                    mbmi->second_mv.as_int = nearest.as_int;
                    vp8_clamp_mv(&mbmi->second_mv, mb_to_left_edge, mb_to_right_edge,
                                 mb_to_top_edge, mb_to_bottom_edge);
                    break;
                case NEWMV:
                    read_mv(bc, &mbmi->second_mv.as_mv, (const MV_CONTEXT *) mvc);
                    mbmi->second_mv.as_mv.row += best_mv.as_mv.row;
                    mbmi->second_mv.as_mv.col += best_mv.as_mv.col;
                    mbmi->need_to_clamp_mvs |= vp8_check_mv_bounds(&mbmi->second_mv,
                                                                   mb_to_left_edge,
                                                                   mb_to_right_edge,
                                                                   mb_to_top_edge,
                                                                   mb_to_bottom_edge);
                    break;
                default:
                    break;
                }
            }
#endif /* CONFIG_DUALPRED */

#if CONFIG_ERROR_CONCEALMENT
            if(pbi->ec_enabled)
            {
                mi->bmi[ 0].mv.as_int =
                mi->bmi[ 1].mv.as_int =
                mi->bmi[ 2].mv.as_int =
                mi->bmi[ 3].mv.as_int =
                mi->bmi[ 4].mv.as_int =
                mi->bmi[ 5].mv.as_int =
                mi->bmi[ 6].mv.as_int =
                mi->bmi[ 7].mv.as_int =
                mi->bmi[ 8].mv.as_int =
                mi->bmi[ 9].mv.as_int =
                mi->bmi[10].mv.as_int =
                mi->bmi[11].mv.as_int =
                mi->bmi[12].mv.as_int =
                mi->bmi[13].mv.as_int =
                mi->bmi[14].mv.as_int =
                mi->bmi[15].mv.as_int = mv->as_int;
            }
#endif
            break;
        default:;
  #if CONFIG_DEBUG
            assert(0);
  #endif
        }
    }
    else
    {
        /* required for left and above block mv */
        mbmi->mv.as_int = 0;

        if ( segfeature_active( xd, mbmi->segment_id, SEG_LVL_MODE ) )
            mbmi->mode = (MB_PREDICTION_MODE)
                         get_segdata( xd, mbmi->segment_id, SEG_LVL_MODE );
        else
        {
            mbmi->mode = (MB_PREDICTION_MODE)
                         vp8_read_ymode(bc, pbi->common.fc.ymode_prob);
        }

        // If MB mode is BPRED read the block modes
        if (mbmi->mode == B_PRED)
        {
            int j = 0;
            do
            {
                mi->bmi[j].as_mode = (B_PREDICTION_MODE)vp8_read_bmode(bc, pbi->common.fc.bmode_prob);
            }
            while (++j < 16);
        }

        if(mbmi->mode == I8X8_PRED)
        {
            int i;
            int mode8x8;
            for(i=0;i<4;i++)
            {
                int ib = vp8_i8x8_block[i];
                mode8x8 = vp8_read_i8x8_mode(bc, pbi->common.i8x8_mode_prob);
                mi->bmi[ib+0].as_mode= mode8x8;
                mi->bmi[ib+1].as_mode= mode8x8;
                mi->bmi[ib+4].as_mode= mode8x8;
                mi->bmi[ib+5].as_mode= mode8x8;
            }
        }
        else
#if CONFIG_UVINTRA
            mbmi->uv_mode = (MB_PREDICTION_MODE)vp8_read_uv_mode(bc,
                                    pbi->common.fc.uv_mode_prob[mbmi->mode]);
#else
            mbmi->uv_mode = (MB_PREDICTION_MODE)vp8_read_uv_mode(bc,
                                    pbi->common.fc.uv_mode_prob);
#endif /*CONFIG_UVINTRA*/
    }

}

#if CONFIG_SUPERBLOCKS
void vp8_decode_mode_mvs(VP8D_COMP *pbi)
{
    int i;
    VP8_COMMON *cm = &pbi->common;
    int sb_row, sb_col;
    int sb_rows = (cm->mb_rows + 1)>>1;
    int sb_cols = (cm->mb_cols + 1)>>1;
    MODE_INFO *mi = cm->mi;
    int row_delta[4] = {-1,  0, +1,  0};
    int col_delta[4] = {+1, +1, -1, +1};

    MODE_INFO *prev_mi = cm->prev_mi;

    mb_mode_mv_init(pbi);

#if CONFIG_QIMODE
    if(cm->frame_type==KEY_FRAME && !cm->kf_ymode_probs_update)
    {
        cm->kf_ymode_probs_index = vp8_read_literal(&pbi->bc, 3);
    }
#endif

    for (sb_row=0; sb_row<sb_rows; sb_row++)
    {
        int mb_col = -col_delta[0];
        int mb_row = (sb_row <<1)-row_delta[0];

        for (sb_col=0; sb_col<sb_cols; sb_col++)
        {
            for ( i=0; i<4; i++ )
            {
                int mb_to_top_edge;
                int mb_to_bottom_edge;

#if CONFIG_ERROR_CONCEALMENT
                int mb_num;
#endif
                int offset_extended = row_delta[(i+1) & 0x3]
                                    * cm->mode_info_stride + col_delta[(i+1) & 0x3];
                int dy = row_delta[i];
                int dx = col_delta[i];

                mb_row += dy;
                mb_col += dx;

                if ((mb_row >= cm->mb_rows) || (mb_col >= cm->mb_cols))
                {
                    prev_mi += offset_extended;
                    mi += offset_extended;       /* next macroblock */
                    continue;
                }

                pbi->mb.mb_to_top_edge =
                mb_to_top_edge = -((mb_row * 16)) << 3;
                mb_to_top_edge -= LEFT_TOP_MARGIN;

                pbi->mb.mb_to_bottom_edge =
                mb_to_bottom_edge = ((pbi->common.mb_rows - 1 - mb_row) * 16) << 3;
                mb_to_bottom_edge += RIGHT_BOTTOM_MARGIN;

#if CONFIG_ERROR_CONCEALMENT
                mb_num = mb_row * pbi->common.mb_cols + mb_col;
#endif
                /*read_mb_modes_mv(pbi, cm->mode_info_context, &cm->mode_info_context->mbmi, mb_row, mb_col);*/
                if(pbi->common.frame_type == KEY_FRAME)
                    vp8_kfread_modes(pbi, mi, mb_row, mb_col);
                else
                    read_mb_modes_mv(pbi, mi, &mi->mbmi,
                    prev_mi,
                    mb_row, mb_col);

#if CONFIG_ERROR_CONCEALMENT
                /* look for corruption. set mvs_corrupt_from_mb to the current
                 * mb_num if the frame is corrupt from this macroblock. */
                if (vp8dx_bool_error(&pbi->bc) && mb_num < pbi->mvs_corrupt_from_mb)
                {
                    pbi->mvs_corrupt_from_mb = mb_num;
                    /* no need to continue since the partition is corrupt from
                     * here on.
                     */
                    return;
                }
#endif

                prev_mi += offset_extended;
                mi += offset_extended;       /* next macroblock */
            }
        }

        mi += cm->mode_info_stride + (1 - (cm->mb_cols & 0x1));
    }
}
#else
void vp8_decode_mode_mvs(VP8D_COMP *pbi)
{
    MODE_INFO *mi = pbi->common.mi;

    MODE_INFO *prev_mi = pbi->common.prev_mi;

    int mb_row = -1;

#if 0
    FILE *statsfile;
    statsfile = fopen("decsegmap.stt", "a");
    fprintf(statsfile, "\n" );
#endif

    mb_mode_mv_init(pbi);

#if CONFIG_QIMODE
    if(pbi->common.frame_type==KEY_FRAME && !pbi->common.kf_ymode_probs_update)
    {
        pbi->common.kf_ymode_probs_index = vp8_read_literal(&pbi->bc, 3);
    }
#endif

    while (++mb_row < pbi->common.mb_rows)
    {
        int mb_col = -1;
        int mb_to_top_edge;
        int mb_to_bottom_edge;

        pbi->mb.mb_to_top_edge =
        mb_to_top_edge = -((mb_row * 16)) << 3;
        mb_to_top_edge -= LEFT_TOP_MARGIN;

        pbi->mb.mb_to_bottom_edge =
        mb_to_bottom_edge = ((pbi->common.mb_rows - 1 - mb_row) * 16) << 3;
        mb_to_bottom_edge += RIGHT_BOTTOM_MARGIN;

#if 0
        fprintf(statsfile, "\n" );
#endif

        while (++mb_col < pbi->common.mb_cols)
        {
#if CONFIG_ERROR_CONCEALMENT
            int mb_num = mb_row * pbi->common.mb_cols + mb_col;
#endif
            /*read_mb_modes_mv(pbi, xd->mode_info_context, &xd->mode_info_context->mbmi, mb_row, mb_col);*/
            if(pbi->common.frame_type == KEY_FRAME)
                vp8_kfread_modes(pbi, mi, mb_row, mb_col);
            else
                read_mb_modes_mv(pbi, mi, &mi->mbmi,
                prev_mi,
                mb_row, mb_col);

            //printf("%3d", mi->mbmi.mode);

            /*
            if(pbi->common.current_video_frame==7)
            {
                FILE *fmode=fopen("kfmode.txt", "a");
                fprintf(fmode, "%3d:%3d:%d\n",mb_row, mb_col, mi->mbmi.mode);
                fclose(fmode);

            }*/
            /*
            if(mi->mbmi.mode==I8X8_PRED)
            {
                printf("F%3d:%d:%d\n", pbi->common.current_video_frame, mb_row, mb_col);
            }
            */
#if CONFIG_ERROR_CONCEALMENT
            /* look for corruption. set mvs_corrupt_from_mb to the current
             * mb_num if the frame is corrupt from this macroblock. */
            if (vp8dx_bool_error(&pbi->bc) && mb_num < pbi->mvs_corrupt_from_mb)
            {
                pbi->mvs_corrupt_from_mb = mb_num;
                /* no need to continue since the partition is corrupt from
                 * here on.
                 */
                return;
            }
#endif

#if 0
            fprintf(statsfile, "%2d%2d%2d   ",
                mi->mbmi.segment_id, mi->mbmi.ref_frame, mi->mbmi.mode );
#endif
            prev_mi++;
            mi++;       /* next macroblock */
        }
       // printf("\n");
        prev_mi++;
        mi++;           /* skip left predictor each row */
    }

#if 0
    fclose(statsfile);
#endif


}
#endif /* CONFIG_SUPERBLOCKS */

