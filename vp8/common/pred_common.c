/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include "vp8/common/pred_common.h"

// TBD prediction functions for various bitstream signals

// Returns a context number for the given MB prediction signal
unsigned char get_pred_context( VP8_COMMON *const cm,
                                MACROBLOCKD *const xd,
                                PRED_ID pred_id )
{
    int pred_context;
    MODE_INFO *m = xd->mode_info_context;

    // Note:
    // The mode info data structure has a one element border above and to the
    // left of the entries correpsonding to real macroblocks.
    // The prediction flags in these dummy entries are initialised to 0.
    switch (pred_id)
    {
    case PRED_SEG_ID:
        pred_context = (m - 1)->mbmi.seg_id_predicted +
                       (m - cm->mode_info_stride)->mbmi.seg_id_predicted;
        break;


#if CONFIG_COMPRED

    case PRED_REF:
        pred_context = (m - 1)->mbmi.ref_predicted +
                       (m - cm->mode_info_stride)->mbmi.ref_predicted;
        break;
#endif

    default:
        // TODO *** add error trap code.
        pred_context = 0;
        break;
    }

    return pred_context;
}

// This function returns a context probability for coding a given
// prediction signal
vp8_prob get_pred_prob( VP8_COMMON *const cm,
                        MACROBLOCKD *const xd,
                        PRED_ID pred_id )
{
    vp8_prob pred_probability;
    int pred_context;

    // Get the appropriate prediction context
    pred_context = get_pred_context( cm, xd, pred_id );

    switch (pred_id)
    {
    case PRED_SEG_ID:
        pred_probability = cm->segment_pred_probs[pred_context];
        break;

#if CONFIG_COMPRED

    case PRED_REF:
        pred_probability = cm->ref_pred_probs[pred_context];
        break;
#endif

    default:
        // TODO *** add error trap code.
        pred_probability = 128;
        break;
    }

    return pred_probability;
}

// This function returns the status of the given prediction signal.
// I.e. is the predicted value for the given signal correct.
unsigned char get_pred_flag( MACROBLOCKD *const xd,
                             PRED_ID pred_id )
{
    unsigned char pred_flag = 0;

    switch (pred_id)
    {
    case PRED_SEG_ID:
        pred_flag = xd->mode_info_context->mbmi.seg_id_predicted;
        break;

#if CONFIG_COMPRED

    case PRED_REF:
        pred_flag = xd->mode_info_context->mbmi.ref_predicted;
        break;
#endif

    default:
        // TODO *** add error trap code.
        pred_flag = 0;
        break;
}

    return pred_flag;
}

// This function sets the status of the given prediction signal.
// I.e. is the predicted value for the given signal correct.
void set_pred_flag( MACROBLOCKD *const xd,
                    PRED_ID pred_id,
                    unsigned char pred_flag)
{
    switch (pred_id)
    {
    case PRED_SEG_ID:
        xd->mode_info_context->mbmi.seg_id_predicted = pred_flag;
        break;

#if CONFIG_COMPRED

    case PRED_REF:
        xd->mode_info_context->mbmi.ref_predicted = pred_flag;
        break;
#endif

    default:
        // TODO *** add error trap code.
        break;
    }
}


// The following contain the guts of the prediction code used to
// peredict various bitstream signals.

// Macroblock segment id prediction function
unsigned char get_pred_mb_segid( VP8_COMMON *const cm, int MbIndex )
{
    // Currently the prediction for the macroblock segment ID is
    // the value stored for this macroblock in the previous frame.
    return cm->last_frame_seg_map[MbIndex];
}

#if CONFIG_COMPRED
MV_REFERENCE_FRAME get_pred_ref( VP8_COMMON *const cm,
                                 MACROBLOCKD *const xd )
{
    MODE_INFO *m = xd->mode_info_context;

    unsigned char left_pred;
    unsigned char above_pred;

    MV_REFERENCE_FRAME left;
    MV_REFERENCE_FRAME above;
    MV_REFERENCE_FRAME above_left;
    MV_REFERENCE_FRAME pred_ref = LAST_FRAME;

    // Reference frame used by neighbours
    left = (m - 1)->mbmi.ref_frame;
    above = (m - cm->mode_info_stride)->mbmi.ref_frame;
    above_left = (m - 1 - cm->mode_info_stride)->mbmi.ref_frame;

    // Reference frame prediction status of immediate neigbours.
    // This can only be set if the mb is "in image"
    left_pred = (m - 1)->mbmi.ref_predicted;
    above_pred = (m - cm->mode_info_stride)->mbmi.ref_predicted;

    // Boost prediction scores of above / left if they are predicted and match
    // the above left.
    if ( left_pred )
        left_pred += (left == above_left);
    if ( above_pred )
        above_pred += (above == above_left);

    // Only consider "in image" mbs as giving valid prediction.
    if ( (left == above) &&
         ((m - 1)->mbmi.mb_in_image ||
          (m - cm->mode_info_stride)->mbmi.mb_in_image) )
    {
        pred_ref = left;
    }
    else if ( left_pred > above_pred )
    {
        pred_ref = left;
    }
    else if ( above_pred > left_pred )
    {
        pred_ref = above;
    }
    else
    {
        // Choose from above or left.
        // For now this is based on a fixed preference order.
        // Last,Altref,Golden
        if ( (left == LAST_FRAME) || (above == LAST_FRAME) )
            pred_ref = LAST_FRAME;
        else if ( (left == ALTREF_FRAME) || (above == ALTREF_FRAME) )
            pred_ref = ALTREF_FRAME;
        else if ( (left == GOLDEN_FRAME) || (above == GOLDEN_FRAME) )
            pred_ref = GOLDEN_FRAME;
    }

    return pred_ref;
}

// Functions to computes a set of modified reference frame probabilities
// to use when the prediction of the reference frame value fails
void calc_ref_probs( int * count, vp8_prob * probs )
{
    int tot_count;

    tot_count = count[0] + count[1] + count[2] + count[3];
    if ( tot_count )
    {
        probs[0] = (vp8_prob)((count[0] * 255) / tot_count);
        probs[0] += !probs[0];
    }
    else
        probs[0] = 128;

    tot_count -= count[0];
    if ( tot_count )
    {
        probs[1] = (vp8_prob)((count[1] * 255) / tot_count);
        probs[1] += !probs[1];
    }
    else
        probs[1] = 128;

    tot_count -= count[1];
    if ( tot_count )
    {
        probs[2] = (vp8_prob)((count[2] * 255) / tot_count);
        probs[2] += !probs[2];
    }
    else
        probs[2] = 128;

}

void compute_mod_refprobs( VP8_COMMON *const cm )
{
    int norm_cnt[MAX_REF_FRAMES];
    int intra_count;
    int inter_count;
    int last_count;
    int gfarf_count;
    int gf_count;
    int arf_count;

    intra_count = cm->prob_intra_coded;
    inter_count = (255 - intra_count);
    last_count = (inter_count * cm->prob_last_coded)/255;
    gfarf_count = inter_count - last_count;
    gf_count = (gfarf_count * cm->prob_gf_coded)/255;
    arf_count = gfarf_count - gf_count;

    // Work out modified reference frame probabilities to use where prediction
    // of the reference frame fails
    norm_cnt[0] = 0;
    norm_cnt[1] = last_count;
    norm_cnt[2] = gf_count;
    norm_cnt[3] = arf_count;
    calc_ref_probs( norm_cnt, cm->mod_refprobs[INTRA_FRAME] );
    cm->mod_refprobs[INTRA_FRAME][0] = 0;    // This branch implicit

    norm_cnt[0] = intra_count;
    norm_cnt[1] = 0;
    norm_cnt[2] = gf_count;
    norm_cnt[3] = arf_count;
    calc_ref_probs( norm_cnt, cm->mod_refprobs[LAST_FRAME]);
    cm->mod_refprobs[LAST_FRAME][1] = 0;    // This branch implicit

    norm_cnt[0] = intra_count;
    norm_cnt[1] = last_count;
    norm_cnt[2] = 0;
    norm_cnt[3] = arf_count;
    calc_ref_probs( norm_cnt, cm->mod_refprobs[GOLDEN_FRAME] );
    cm->mod_refprobs[GOLDEN_FRAME][2] = 0;  // This branch implicit

    norm_cnt[0] = intra_count;
    norm_cnt[1] = last_count;
    norm_cnt[2] = gf_count;
    norm_cnt[3] = 0;
    calc_ref_probs( norm_cnt, cm->mod_refprobs[ALTREF_FRAME] );
    cm->mod_refprobs[ALTREF_FRAME][2] = 0;  // This branch implicit

}
#endif
