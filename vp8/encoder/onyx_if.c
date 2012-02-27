/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */


#include "vpx_config.h"
#include "vp8/common/onyxc_int.h"
#include "onyx_int.h"
#include "vp8/common/systemdependent.h"
#include "quantize.h"
#include "vp8/common/alloccommon.h"
#include "mcomp.h"
#include "firstpass.h"
#include "psnr.h"
#include "vpx_scale/vpxscale.h"
#include "vp8/common/extend.h"
#include "ratectrl.h"
#include "vp8/common/quant_common.h"
#include "segmentation.h"
#include "vp8/common/g_common.h"
#include "vpx_scale/yv12extend.h"
#if CONFIG_POSTPROC
#include "vp8/common/postproc.h"
#endif
#include "vpx_mem/vpx_mem.h"
#include "vp8/common/swapyv12buffer.h"
#include "vpx_ports/vpx_timer.h"
#include "temporal_filter.h"

#include "vp8/common/seg_common.h"
#include "mbgraph.h"
#include "vp8/common/pred_common.h"

#if ARCH_ARM
#include "vpx_ports/arm.h"
#endif

#include <math.h>
#include <stdio.h>
#include <limits.h>

#if CONFIG_RUNTIME_CPU_DETECT
#define IF_RTCD(x) (x)
#define RTCD(x) &cpi->common.rtcd.x
#else
#define IF_RTCD(x) NULL
#define RTCD(x) NULL
#endif

extern void vp8cx_pick_filter_level_fast(YV12_BUFFER_CONFIG *sd, VP8_COMP *cpi);
extern void vp8cx_set_alt_lf_level(VP8_COMP *cpi, int filt_val);
extern void vp8cx_pick_filter_level(YV12_BUFFER_CONFIG *sd, VP8_COMP *cpi);

extern void vp8_dmachine_specific_config(VP8_COMP *cpi);
extern void vp8_cmachine_specific_config(VP8_COMP *cpi);
extern void vp8_deblock_frame(YV12_BUFFER_CONFIG *source, YV12_BUFFER_CONFIG *post, int filt_lvl, int low_var_thresh, int flag);
extern void print_parms(VP8_CONFIG *ocf, char *filenam);
extern unsigned int vp8_get_processor_freq();
extern void print_tree_update_probs();
extern void vp8cx_create_encoder_threads(VP8_COMP *cpi);
extern void vp8cx_remove_encoder_threads(VP8_COMP *cpi);
#if HAVE_ARMV7
extern void vp8_yv12_copy_frame_func_neon(YV12_BUFFER_CONFIG *src_ybc, YV12_BUFFER_CONFIG *dst_ybc);
extern void vp8_yv12_copy_src_frame_func_neon(YV12_BUFFER_CONFIG *src_ybc, YV12_BUFFER_CONFIG *dst_ybc);
#endif

int vp8_estimate_entropy_savings(VP8_COMP *cpi);
int vp8_calc_ss_err(YV12_BUFFER_CONFIG *source, YV12_BUFFER_CONFIG *dest, const vp8_variance_rtcd_vtable_t *rtcd);

extern void vp8_temporal_filter_prepare_c(VP8_COMP *cpi, int distance);

static void set_default_lf_deltas(VP8_COMP *cpi);

extern const int vp8_gf_interval_table[101];

#if CONFIG_INTERNAL_STATS
#include "math.h"

extern double vp8_calc_ssim
(
    YV12_BUFFER_CONFIG *source,
    YV12_BUFFER_CONFIG *dest,
    int lumamask,
    double *weight,
    const vp8_variance_rtcd_vtable_t *rtcd
);


extern double vp8_calc_ssimg
(
    YV12_BUFFER_CONFIG *source,
    YV12_BUFFER_CONFIG *dest,
    double *ssim_y,
    double *ssim_u,
    double *ssim_v,
    const vp8_variance_rtcd_vtable_t *rtcd
);


#endif

//#define OUTPUT_YUV_REC

#ifdef OUTPUT_YUV_SRC
FILE *yuv_file;
#endif
#ifdef OUTPUT_YUV_REC
FILE *yuv_rec_file;
#endif

#if 0
FILE *framepsnr;
FILE *kf_list;
FILE *keyfile;
#endif

#if 0
extern int skip_true_count;
extern int skip_false_count;
#endif


#ifdef ENTROPY_STATS
extern int intra_mode_stats[10][10][10];
#endif

#ifdef SPEEDSTATS
unsigned int frames_at_speed[16] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
unsigned int tot_pm = 0;
unsigned int cnt_pm = 0;
unsigned int tot_ef = 0;
unsigned int cnt_ef = 0;
#endif

#if defined(SECTIONBITS_OUTPUT)
extern unsigned __int64 Sectionbits[500];
#endif
#ifdef MODE_STATS
extern INT64 Sectionbits[500];
extern int y_modes[VP8_YMODES]  ;
extern int i8x8_modes[VP8_I8X8_MODES];
extern int uv_modes[VP8_UV_MODES] ;
extern int uv_modes_y[VP8_YMODES][VP8_UV_MODES];
extern int b_modes[B_MODE_COUNT];
extern int inter_y_modes[MB_MODE_COUNT] ;
extern int inter_uv_modes[VP8_UV_MODES] ;
extern unsigned int inter_b_modes[B_MODE_COUNT];
#endif

extern void (*vp8_short_fdct4x4)(short *input, short *output, int pitch);
extern void (*vp8_short_fdct8x4)(short *input, short *output, int pitch);

extern void vp8cx_init_quantizer(VP8_COMP *cpi);

int vp8cx_base_skip_false_prob[QINDEX_RANGE];

// Tables relating active max Q to active min Q
static int kf_low_motion_minq[QINDEX_RANGE];
static int kf_high_motion_minq[QINDEX_RANGE];
static int gf_low_motion_minq[QINDEX_RANGE];
static int gf_mid_motion_minq[QINDEX_RANGE];
static int gf_high_motion_minq[QINDEX_RANGE];
static int inter_minq[QINDEX_RANGE];

// Functions to compute the active minq lookup table entries based on a
// formulaic approach to facilitate easier adjustment of the Q tables.
// The formulae were derived from computing a 3rd order polynomial best
// fit to the original data (after plotting real maxq vs minq (not q index))
int calculate_minq_index( double maxq,
                          double x3, double x2, double x, double c )
{
    int i;
    double minqtarget;
    double thisq;

    minqtarget = ( (x3 * maxq * maxq * maxq) +
                   (x2 * maxq * maxq) +
                   (x * maxq) +
                   c );

    if ( minqtarget > maxq )
        minqtarget = maxq;

    for ( i = 0; i < QINDEX_RANGE; i++ )
    {
        thisq = vp8_convert_qindex_to_q(i);
        if ( minqtarget <= vp8_convert_qindex_to_q(i) )
            return i;
    }
    if ( i == QINDEX_RANGE )
        return QINDEX_RANGE-1;
}
void init_minq_luts()
{
    int i;
    double maxq;

    for ( i = 0; i < QINDEX_RANGE; i++ )
    {
        maxq = vp8_convert_qindex_to_q(i);


        kf_low_motion_minq[i] = calculate_minq_index( maxq,
                                                      0.0000003,
                                                      -0.000015,
                                                      0.074,
                                                      0.0 );

        kf_high_motion_minq[i] = calculate_minq_index( maxq,
                                                       0.00000034,
                                                       -0.000125,
                                                       0.13,
                                                       0.0 );
        gf_low_motion_minq[i] = calculate_minq_index( maxq,
                                                      0.0000016,
                                                      -0.00078,
                                                      0.315,
                                                      0.0 );
        gf_mid_motion_minq[i] = calculate_minq_index( maxq,
                                                      0.00000415,
                                                      -0.0017,
                                                      0.425,
                                                      0.0 );
        gf_high_motion_minq[i] = calculate_minq_index( maxq,
                                                       0.00000725,
                                                       -0.00235,
                                                       0.47,
                                                       0.0  );
        inter_minq[i] = calculate_minq_index( maxq,
                                              0.00000271,
                                              -0.00113,
                                              0.697,
                                              0.0  );

    }
}

void init_base_skip_probs()
{
    int i;
    double q;
    int skip_prob;

    for ( i = 0; i < QINDEX_RANGE; i++ )
    {
        q = vp8_convert_qindex_to_q(i);

        // Exponential decay caluclation of baseline skip prob with clamping
        // Based on crude best fit of old table.
        skip_prob = (int)( 564.25 * pow( 2.71828, (-0.012*q) ) );
        if ( skip_prob < 1 )
            skip_prob = 1;
        else if ( skip_prob > 255 )
            skip_prob = 255;

        vp8cx_base_skip_false_prob[i] = skip_prob;
    }
}

void vp8_initialize()
{
    static int init_done = 0;

    if (!init_done)
    {
        vp8_scale_machine_specific_config();
        vp8_initialize_common();
        //vp8_dmachine_specific_config();
        vp8_tokenize_initialize();
        vp8_init_quant_tables();
        vp8_init_me_luts();
        init_minq_luts();
        init_base_skip_probs();
        init_done = 1;
    }
}
#ifdef PACKET_TESTING
extern FILE *vpxlogc;
#endif

static void setup_features(VP8_COMP *cpi)
{
    MACROBLOCKD *xd = &cpi->mb.e_mbd;

    // Set up default state for MB feature flags

    xd->segmentation_enabled = 0;   // Default segmentation disabled

    xd->update_mb_segmentation_map = 0;
    xd->update_mb_segmentation_data = 0;
    vpx_memset(xd->mb_segment_tree_probs, 255, sizeof(xd->mb_segment_tree_probs));

    clearall_segfeatures( xd );

    xd->mode_ref_lf_delta_enabled = 0;
    xd->mode_ref_lf_delta_update = 0;
    vpx_memset(xd->ref_lf_deltas, 0, sizeof(xd->ref_lf_deltas));
    vpx_memset(xd->mode_lf_deltas, 0, sizeof(xd->mode_lf_deltas));
    vpx_memset(xd->last_ref_lf_deltas, 0, sizeof(xd->ref_lf_deltas));
    vpx_memset(xd->last_mode_lf_deltas, 0, sizeof(xd->mode_lf_deltas));

    set_default_lf_deltas(cpi);

}


static void dealloc_compressor_data(VP8_COMP *cpi)
{
    vpx_free(cpi->tplist);
    cpi->tplist = NULL;

    // Delete last frame MV storage buffers
    vpx_free(cpi->lfmv);
    cpi->lfmv = 0;

    vpx_free(cpi->lf_ref_frame_sign_bias);
    cpi->lf_ref_frame_sign_bias = 0;

    vpx_free(cpi->lf_ref_frame);
    cpi->lf_ref_frame = 0;

    // Delete sementation map
    vpx_free(cpi->segmentation_map);
    cpi->segmentation_map = 0;
    vpx_free(cpi->common.last_frame_seg_map);
    cpi->common.last_frame_seg_map = 0;

    vpx_free(cpi->active_map);
    cpi->active_map = 0;

    vp8_de_alloc_frame_buffers(&cpi->common);

    vp8_yv12_de_alloc_frame_buffer(&cpi->last_frame_uf);
    vp8_yv12_de_alloc_frame_buffer(&cpi->scaled_source);
#if VP8_TEMPORAL_ALT_REF
    vp8_yv12_de_alloc_frame_buffer(&cpi->alt_ref_buffer);
#endif
    vp8_lookahead_destroy(cpi->lookahead);

    vpx_free(cpi->tok);
    cpi->tok = 0;

    // Structure used to monitor GF usage
    vpx_free(cpi->gf_active_flags);
    cpi->gf_active_flags = 0;

    // Activity mask based per mb zbin adjustments
    vpx_free(cpi->mb_activity_map);
    cpi->mb_activity_map = 0;
    vpx_free(cpi->mb_norm_activity_map);
    cpi->mb_norm_activity_map = 0;

    vpx_free(cpi->mb.pip);
    cpi->mb.pip = 0;

    vpx_free(cpi->twopass.total_stats);
    cpi->twopass.total_stats = 0;

    vpx_free(cpi->twopass.total_left_stats);
    cpi->twopass.total_left_stats = 0;

    vpx_free(cpi->twopass.this_frame_stats);
    cpi->twopass.this_frame_stats = 0;
}

// Computes a q delta (in "q index" terms) to get from a starting q value
// to a target value
// target q value
static int compute_qdelta( VP8_COMP *cpi, double qstart, double qtarget )
{
    int i;
    int start_index = cpi->worst_quality;
    int target_index = cpi->worst_quality;
    int retval = 0;

    // Convert the average q value to an index.
    for ( i = cpi->best_quality; i < cpi->worst_quality; i++ )
    {
        start_index = i;
        if ( vp8_convert_qindex_to_q(i) >= qstart )
            break;
    }

    // Convert the q target to an index
    for ( i = cpi->best_quality; i < cpi->worst_quality; i++ )
    {
        target_index = i;
        if ( vp8_convert_qindex_to_q(i) >= qtarget )
            break;
    }

    return target_index - start_index;
}

static void init_seg_features(VP8_COMP *cpi)
{
    VP8_COMMON *cm = &cpi->common;
    MACROBLOCKD *xd = &cpi->mb.e_mbd;

    int high_q = (int)(cpi->avg_q > 48.0);
    int qi_delta;

    // Disable and clear down for KF
    if ( cm->frame_type == KEY_FRAME  )
    {
        // Clear down the global segmentation map
        vpx_memset( cpi->segmentation_map, 0, (cm->mb_rows * cm->mb_cols));
        xd->update_mb_segmentation_map = 0;
        xd->update_mb_segmentation_data = 0;
        cpi->static_mb_pct = 0;

        // Disable segmentation
        vp8_disable_segmentation((VP8_PTR)cpi);

        // Clear down the segment features.
        clearall_segfeatures(xd);
    }

    // If this is an alt ref frame
    else if ( cm->refresh_alt_ref_frame )
    {
        // Clear down the global segmentation map
        vpx_memset( cpi->segmentation_map, 0, (cm->mb_rows * cm->mb_cols));
        xd->update_mb_segmentation_map = 0;
        xd->update_mb_segmentation_data = 0;
        cpi->static_mb_pct = 0;

        // Disable segmentation and individual segment features by default
        vp8_disable_segmentation((VP8_PTR)cpi);
        clearall_segfeatures(xd);

        // Scan frames from current to arf frame.
        // This function re-enables segmentation if appropriate.
        vp8_update_mbgraph_stats(cpi);

        // If segmentation was enabled set those features needed for the
        // arf itself.
        if ( xd->segmentation_enabled )
        {
            xd->update_mb_segmentation_map = 1;
            xd->update_mb_segmentation_data = 1;

            qi_delta = compute_qdelta( cpi, cpi->avg_q, (cpi->avg_q * 0.875) );
            set_segdata( xd, 1, SEG_LVL_ALT_Q, (qi_delta - 2) );
            set_segdata( xd, 1, SEG_LVL_ALT_LF, -2 );

            enable_segfeature(xd, 1, SEG_LVL_ALT_Q);
            enable_segfeature(xd, 1, SEG_LVL_ALT_LF);

            // Where relevant assume segment data is delta data
            xd->mb_segment_abs_delta = SEGMENT_DELTADATA;

        }
    }
    // All other frames if segmentation has been enabled
    else if ( xd->segmentation_enabled )
    {
/*
        int i;

        // clears prior frame seg lev refs
        for (i = 0; i < MAX_MB_SEGMENTS; i++)
        {
            // only do it if the force drop the background stuff is off
            if(!segfeature_active(xd, i, SEG_LVL_MODE))
            {
                disable_segfeature(xd,i,SEG_LVL_REF_FRAME);
                set_segdata( xd,i, SEG_LVL_REF_FRAME, 0xffffff);
            }
        }
*/

        // First normal frame in a valid gf or alt ref group
        if ( cpi->common.frames_since_golden == 0 )
        {
            // Set up segment features for normal frames in an af group
            if ( cpi->source_alt_ref_active )
            {
                xd->update_mb_segmentation_map = 0;
                xd->update_mb_segmentation_data = 1;
                xd->mb_segment_abs_delta = SEGMENT_DELTADATA;

                qi_delta = compute_qdelta( cpi, cpi->avg_q,
                                           (cpi->avg_q * 1.125) );
                set_segdata( xd, 1, SEG_LVL_ALT_Q, (qi_delta + 2) );
                set_segdata( xd, 1, SEG_LVL_ALT_Q, 0 );
                enable_segfeature(xd, 1, SEG_LVL_ALT_Q);

                set_segdata( xd, 1, SEG_LVL_ALT_LF, -2 );
                enable_segfeature(xd, 1, SEG_LVL_ALT_LF);

                // Segment coding disabled for compred testing
                if ( high_q || (cpi->static_mb_pct == 100) )
                {
                    //set_segref(xd, 1, LAST_FRAME);
                    set_segref(xd, 1, ALTREF_FRAME);
                    enable_segfeature(xd, 1, SEG_LVL_REF_FRAME);

                    set_segdata( xd, 1, SEG_LVL_MODE, ZEROMV );
                    enable_segfeature(xd, 1, SEG_LVL_MODE);

                    // EOB segment coding not fixed for 8x8 yet
                    set_segdata( xd, 1, SEG_LVL_EOB, 0 );
                    enable_segfeature(xd, 1, SEG_LVL_EOB);
                }
            }
            // Disable segmentation and clear down features if alt ref
            // is not active for this group
            else
            {
                vp8_disable_segmentation((VP8_PTR)cpi);

                vpx_memset( cpi->segmentation_map, 0,
                            (cm->mb_rows * cm->mb_cols));

                xd->update_mb_segmentation_map = 0;
                xd->update_mb_segmentation_data = 0;

                clearall_segfeatures(xd);
            }
        }

        // Special case where we are coding over the top of a previous
        // alt ref frame
        // Segment coding disabled for compred testing
        else if ( cpi->is_src_frame_alt_ref )
        {
            // Enable mode and ref frame features for segment 0 as well
            enable_segfeature(xd, 0, SEG_LVL_REF_FRAME);
            enable_segfeature(xd, 0, SEG_LVL_MODE);
            enable_segfeature(xd, 1, SEG_LVL_REF_FRAME);
            enable_segfeature(xd, 1, SEG_LVL_MODE);

            // All mbs should use ALTREF_FRAME, ZEROMV exclusively
            clear_segref(xd, 0);
            set_segref(xd, 0, ALTREF_FRAME);
            clear_segref(xd, 1);
            set_segref(xd, 1, ALTREF_FRAME);
            set_segdata( xd, 0, SEG_LVL_MODE, ZEROMV );
            set_segdata( xd, 1, SEG_LVL_MODE, ZEROMV );

            // Skip all MBs if high Q
            if ( high_q )
            {
                enable_segfeature(xd, 0, SEG_LVL_EOB);
                set_segdata( xd, 0, SEG_LVL_EOB, 0 );
                enable_segfeature(xd, 1, SEG_LVL_EOB);
                set_segdata( xd, 1, SEG_LVL_EOB, 0 );
            }
            // Enable data udpate
            xd->update_mb_segmentation_data = 1;
        }
        // All other frames.
        else
        {
            // No updeates.. leave things as they are.
            xd->update_mb_segmentation_map = 0;
            xd->update_mb_segmentation_data = 0;
        }
    }
}

// DEBUG: Print out the segment id of each MB in the current frame.
static void print_seg_map(VP8_COMP *cpi)
{
    VP8_COMMON *cm = & cpi->common;
    int row,col;
    int map_index = 0;
    FILE *statsfile;

    statsfile = fopen("segmap.stt", "a");

    fprintf(statsfile, "%10d\n",
            cm->current_video_frame );

    for ( row = 0; row < cpi->common.mb_rows; row++ )
    {
        for ( col = 0; col < cpi->common.mb_cols; col++ )
        {
            fprintf(statsfile, "%10d",
                    cpi->segmentation_map[map_index]);
            map_index++;
        }
        fprintf(statsfile, "\n");
    }
    fprintf(statsfile, "\n");

    fclose(statsfile);
}

static void set_default_lf_deltas(VP8_COMP *cpi)
{
    cpi->mb.e_mbd.mode_ref_lf_delta_enabled = 1;
    cpi->mb.e_mbd.mode_ref_lf_delta_update = 1;

    vpx_memset(cpi->mb.e_mbd.ref_lf_deltas, 0, sizeof(cpi->mb.e_mbd.ref_lf_deltas));
    vpx_memset(cpi->mb.e_mbd.mode_lf_deltas, 0, sizeof(cpi->mb.e_mbd.mode_lf_deltas));

    // Test of ref frame deltas
    cpi->mb.e_mbd.ref_lf_deltas[INTRA_FRAME] = 2;
    cpi->mb.e_mbd.ref_lf_deltas[LAST_FRAME] = 0;
    cpi->mb.e_mbd.ref_lf_deltas[GOLDEN_FRAME] = -2;
    cpi->mb.e_mbd.ref_lf_deltas[ALTREF_FRAME] = -2;

    cpi->mb.e_mbd.mode_lf_deltas[0] = 4;               // BPRED
    cpi->mb.e_mbd.mode_lf_deltas[1] = -2;              // Zero
    cpi->mb.e_mbd.mode_lf_deltas[2] = 2;               // New mv
    cpi->mb.e_mbd.mode_lf_deltas[3] = 4;               // Split mv
}

void vp8_set_speed_features(VP8_COMP *cpi)
{
    SPEED_FEATURES *sf = &cpi->sf;
    int Mode = cpi->compressor_speed;
    int Speed = cpi->Speed;
    int i;
    VP8_COMMON *cm = &cpi->common;
    int last_improved_quant = sf->improved_quant;

    // Only modes 0 and 1 supported for now in experimental code basae
    if ( Mode > 1 )
        Mode = 1;

    // Initialise default mode frequency sampling variables
    for (i = 0; i < MAX_MODES; i ++)
    {
        cpi->mode_check_freq[i] = 0;
        cpi->mode_test_hit_counts[i] = 0;
        cpi->mode_chosen_counts[i] = 0;
    }

    cpi->mbs_tested_so_far = 0;

    // best quality defaults
    sf->RD = 1;
    sf->search_method = NSTEP;
    sf->improved_quant = 1;
    sf->improved_dct = 1;
    sf->auto_filter = 1;
    sf->recode_loop = 1;
    sf->quarter_pixel_search = 1;
    sf->half_pixel_search = 1;
    sf->iterative_sub_pixel = 1;
    sf->optimize_coefficients = 1;
    sf->use_fastquant_for_pick = 0;
    sf->no_skip_block4x4_search = 1;

    sf->first_step = 0;
    sf->max_step_search_steps = MAX_MVSEARCH_STEPS;
    sf->improved_mv_pred = 1;

    // default thresholds to 0
    for (i = 0; i < MAX_MODES; i++)
        sf->thresh_mult[i] = 0;

    switch (Mode)
    {
    case 0: // best quality mode
        sf->thresh_mult[THR_ZEROMV   ] = 0;
        sf->thresh_mult[THR_ZEROG    ] = 0;
        sf->thresh_mult[THR_ZEROA    ] = 0;
        sf->thresh_mult[THR_NEARESTMV] = 0;
        sf->thresh_mult[THR_NEARESTG ] = 0;
        sf->thresh_mult[THR_NEARESTA ] = 0;
        sf->thresh_mult[THR_NEARMV   ] = 0;
        sf->thresh_mult[THR_NEARG    ] = 0;
        sf->thresh_mult[THR_NEARA    ] = 0;

        sf->thresh_mult[THR_DC       ] = 0;

        sf->thresh_mult[THR_V_PRED   ] = 1000;
        sf->thresh_mult[THR_H_PRED   ] = 1000;
        sf->thresh_mult[THR_B_PRED   ] = 2000;
        sf->thresh_mult[THR_I8X8_PRED] = 2000;
        sf->thresh_mult[THR_TM       ] = 1000;

        sf->thresh_mult[THR_NEWMV    ] = 1000;
        sf->thresh_mult[THR_NEWG     ] = 1000;
        sf->thresh_mult[THR_NEWA     ] = 1000;

        sf->thresh_mult[THR_SPLITMV  ] = 2500;
        sf->thresh_mult[THR_SPLITG   ] = 5000;
        sf->thresh_mult[THR_SPLITA   ] = 5000;

        sf->thresh_mult[THR_DUAL_ZEROLG   ] = 0;
        sf->thresh_mult[THR_DUAL_NEARESTLG] = 0;
        sf->thresh_mult[THR_DUAL_NEARLG   ] = 0;
        sf->thresh_mult[THR_DUAL_ZEROLA   ] = 0;
        sf->thresh_mult[THR_DUAL_NEARESTLA] = 0;
        sf->thresh_mult[THR_DUAL_NEARLA   ] = 0;
        sf->thresh_mult[THR_DUAL_ZEROGA   ] = 0;
        sf->thresh_mult[THR_DUAL_NEARESTGA] = 0;
        sf->thresh_mult[THR_DUAL_NEARGA   ] = 0;

        sf->thresh_mult[THR_DUAL_NEWLG    ] = 1000;
        sf->thresh_mult[THR_DUAL_NEWLA    ] = 1000;
        sf->thresh_mult[THR_DUAL_NEWGA    ] = 1000;

        sf->first_step = 0;
        sf->max_step_search_steps = MAX_MVSEARCH_STEPS;
        break;
    case 1:
        sf->thresh_mult[THR_NEARESTMV] = 0;
        sf->thresh_mult[THR_ZEROMV   ] = 0;
        sf->thresh_mult[THR_DC       ] = 0;
        sf->thresh_mult[THR_NEARMV   ] = 0;
        sf->thresh_mult[THR_V_PRED   ] = 1000;
        sf->thresh_mult[THR_H_PRED   ] = 1000;
        sf->thresh_mult[THR_B_PRED   ] = 2500;
        sf->thresh_mult[THR_I8X8_PRED] = 2500;
        sf->thresh_mult[THR_TM       ] = 1000;

        sf->thresh_mult[THR_NEARESTG ] = 1000;
        sf->thresh_mult[THR_NEARESTA ] = 1000;

        sf->thresh_mult[THR_ZEROG    ] = 1000;
        sf->thresh_mult[THR_ZEROA    ] = 1000;
        sf->thresh_mult[THR_NEARG    ] = 1000;
        sf->thresh_mult[THR_NEARA    ] = 1000;

        sf->thresh_mult[THR_ZEROMV   ] = 0;
        sf->thresh_mult[THR_ZEROG    ] = 0;
        sf->thresh_mult[THR_ZEROA    ] = 0;
        sf->thresh_mult[THR_NEARESTMV] = 0;
        sf->thresh_mult[THR_NEARESTG ] = 0;
        sf->thresh_mult[THR_NEARESTA ] = 0;
        sf->thresh_mult[THR_NEARMV   ] = 0;
        sf->thresh_mult[THR_NEARG    ] = 0;
        sf->thresh_mult[THR_NEARA    ] = 0;

        sf->thresh_mult[THR_NEWMV    ] = 1000;
        sf->thresh_mult[THR_NEWG     ] = 1000;
        sf->thresh_mult[THR_NEWA     ] = 1000;

        sf->thresh_mult[THR_SPLITMV  ] = 1700;
        sf->thresh_mult[THR_SPLITG   ] = 4500;
        sf->thresh_mult[THR_SPLITA   ] = 4500;

        sf->thresh_mult[THR_DUAL_ZEROLG   ] = 0;
        sf->thresh_mult[THR_DUAL_NEARESTLG] = 0;
        sf->thresh_mult[THR_DUAL_NEARLG   ] = 0;
        sf->thresh_mult[THR_DUAL_ZEROLA   ] = 0;
        sf->thresh_mult[THR_DUAL_NEARESTLA] = 0;
        sf->thresh_mult[THR_DUAL_NEARLA   ] = 0;
        sf->thresh_mult[THR_DUAL_ZEROGA   ] = 0;
        sf->thresh_mult[THR_DUAL_NEARESTGA] = 0;
        sf->thresh_mult[THR_DUAL_NEARGA   ] = 0;

        sf->thresh_mult[THR_DUAL_NEWLG    ] = 1000;
        sf->thresh_mult[THR_DUAL_NEWLA    ] = 1000;
        sf->thresh_mult[THR_DUAL_NEWGA    ] = 1000;

        if (Speed > 0)
        {
            /* Disable coefficient optimization above speed 0 */
            sf->optimize_coefficients = 0;
            sf->use_fastquant_for_pick = 1;
            sf->no_skip_block4x4_search = 0;

            sf->first_step = 1;

            cpi->mode_check_freq[THR_SPLITG] = 2;
            cpi->mode_check_freq[THR_SPLITA] = 2;
            cpi->mode_check_freq[THR_SPLITMV] = 0;
        }

        if (Speed > 1)
        {
            cpi->mode_check_freq[THR_SPLITG] = 4;
            cpi->mode_check_freq[THR_SPLITA] = 4;
            cpi->mode_check_freq[THR_SPLITMV] = 2;

            sf->thresh_mult[THR_TM       ] = 1500;
            sf->thresh_mult[THR_V_PRED   ] = 1500;
            sf->thresh_mult[THR_H_PRED   ] = 1500;
            sf->thresh_mult[THR_B_PRED   ] = 5000;
            sf->thresh_mult[THR_I8X8_PRED] = 5000;

            if (cpi->ref_frame_flags & VP8_LAST_FLAG)
            {
                sf->thresh_mult[THR_NEWMV    ] = 2000;
                sf->thresh_mult[THR_SPLITMV  ] = 10000;
            }

            if (cpi->ref_frame_flags & VP8_GOLD_FLAG)
            {
                sf->thresh_mult[THR_NEARESTG ] = 1500;
                sf->thresh_mult[THR_ZEROG    ] = 1500;
                sf->thresh_mult[THR_NEARG    ] = 1500;
                sf->thresh_mult[THR_NEWG     ] = 2000;
                sf->thresh_mult[THR_SPLITG   ] = 20000;
            }

            if (cpi->ref_frame_flags & VP8_ALT_FLAG)
            {
                sf->thresh_mult[THR_NEARESTA ] = 1500;
                sf->thresh_mult[THR_ZEROA    ] = 1500;
                sf->thresh_mult[THR_NEARA    ] = 1500;
                sf->thresh_mult[THR_NEWA     ] = 2000;
                sf->thresh_mult[THR_SPLITA   ] = 20000;
            }

            sf->thresh_mult[THR_DUAL_ZEROLG   ] = 1500;
            sf->thresh_mult[THR_DUAL_NEARESTLG] = 1500;
            sf->thresh_mult[THR_DUAL_NEARLG   ] = 1500;
            sf->thresh_mult[THR_DUAL_ZEROLA   ] = 1500;
            sf->thresh_mult[THR_DUAL_NEARESTLA] = 1500;
            sf->thresh_mult[THR_DUAL_NEARLA   ] = 1500;
            sf->thresh_mult[THR_DUAL_ZEROGA   ] = 1500;
            sf->thresh_mult[THR_DUAL_NEARESTGA] = 1500;
            sf->thresh_mult[THR_DUAL_NEARGA   ] = 1500;

            sf->thresh_mult[THR_DUAL_NEWLG    ] = 2000;
            sf->thresh_mult[THR_DUAL_NEWLA    ] = 2000;
            sf->thresh_mult[THR_DUAL_NEWGA    ] = 2000;
        }

        if (Speed > 2)
        {
            cpi->mode_check_freq[THR_SPLITG] = 15;
            cpi->mode_check_freq[THR_SPLITA] = 15;
            cpi->mode_check_freq[THR_SPLITMV] = 7;

            sf->thresh_mult[THR_TM       ] = 2000;
            sf->thresh_mult[THR_V_PRED   ] = 2000;
            sf->thresh_mult[THR_H_PRED   ] = 2000;
            sf->thresh_mult[THR_B_PRED   ] = 7500;
            sf->thresh_mult[THR_I8X8_PRED] = 7500;

            if (cpi->ref_frame_flags & VP8_LAST_FLAG)
            {
                sf->thresh_mult[THR_NEWMV    ] = 2000;
                sf->thresh_mult[THR_SPLITMV  ] = 25000;
            }

            if (cpi->ref_frame_flags & VP8_GOLD_FLAG)
            {
                sf->thresh_mult[THR_NEARESTG ] = 2000;
                sf->thresh_mult[THR_ZEROG    ] = 2000;
                sf->thresh_mult[THR_NEARG    ] = 2000;
                sf->thresh_mult[THR_NEWG     ] = 2500;
                sf->thresh_mult[THR_SPLITG   ] = 50000;
            }

            if (cpi->ref_frame_flags & VP8_ALT_FLAG)
            {
                sf->thresh_mult[THR_NEARESTA ] = 2000;
                sf->thresh_mult[THR_ZEROA    ] = 2000;
                sf->thresh_mult[THR_NEARA    ] = 2000;
                sf->thresh_mult[THR_NEWA     ] = 2500;
                sf->thresh_mult[THR_SPLITA   ] = 50000;
            }

            sf->thresh_mult[THR_DUAL_ZEROLG   ] = 2000;
            sf->thresh_mult[THR_DUAL_NEARESTLG] = 2000;
            sf->thresh_mult[THR_DUAL_NEARLG   ] = 2000;
            sf->thresh_mult[THR_DUAL_ZEROLA   ] = 2000;
            sf->thresh_mult[THR_DUAL_NEARESTLA] = 2000;
            sf->thresh_mult[THR_DUAL_NEARLA   ] = 2000;
            sf->thresh_mult[THR_DUAL_ZEROGA   ] = 2000;
            sf->thresh_mult[THR_DUAL_NEARESTGA] = 2000;
            sf->thresh_mult[THR_DUAL_NEARGA   ] = 2000;

            sf->thresh_mult[THR_DUAL_NEWLG    ] = 2500;
            sf->thresh_mult[THR_DUAL_NEWLA    ] = 2500;
            sf->thresh_mult[THR_DUAL_NEWGA    ] = 2500;

            sf->improved_quant = 0;
            sf->improved_dct = 0;

            // Only do recode loop on key frames, golden frames and
            // alt ref frames
            sf->recode_loop = 2;

        }

        break;

    }; /* switch */

    /* disable frame modes if flags not set */
    if (!(cpi->ref_frame_flags & VP8_LAST_FLAG))
    {
        sf->thresh_mult[THR_NEWMV    ] = INT_MAX;
        sf->thresh_mult[THR_NEARESTMV] = INT_MAX;
        sf->thresh_mult[THR_ZEROMV   ] = INT_MAX;
        sf->thresh_mult[THR_NEARMV   ] = INT_MAX;
        sf->thresh_mult[THR_SPLITMV  ] = INT_MAX;
    }

    if (!(cpi->ref_frame_flags & VP8_GOLD_FLAG))
    {
        sf->thresh_mult[THR_NEARESTG ] = INT_MAX;
        sf->thresh_mult[THR_ZEROG    ] = INT_MAX;
        sf->thresh_mult[THR_NEARG    ] = INT_MAX;
        sf->thresh_mult[THR_NEWG     ] = INT_MAX;
        sf->thresh_mult[THR_SPLITG   ] = INT_MAX;
    }

    if (!(cpi->ref_frame_flags & VP8_ALT_FLAG))
    {
        sf->thresh_mult[THR_NEARESTA ] = INT_MAX;
        sf->thresh_mult[THR_ZEROA    ] = INT_MAX;
        sf->thresh_mult[THR_NEARA    ] = INT_MAX;
        sf->thresh_mult[THR_NEWA     ] = INT_MAX;
        sf->thresh_mult[THR_SPLITA   ] = INT_MAX;
    }

    if ((cpi->ref_frame_flags & (VP8_LAST_FLAG | VP8_GOLD_FLAG)) != (VP8_LAST_FLAG | VP8_GOLD_FLAG))
    {
        sf->thresh_mult[THR_DUAL_ZEROLG   ] = INT_MAX;
        sf->thresh_mult[THR_DUAL_NEARESTLG] = INT_MAX;
        sf->thresh_mult[THR_DUAL_NEARLG   ] = INT_MAX;
        sf->thresh_mult[THR_DUAL_NEWLG    ] = INT_MAX;
    }

    if ((cpi->ref_frame_flags & (VP8_LAST_FLAG | VP8_ALT_FLAG)) != (VP8_LAST_FLAG | VP8_ALT_FLAG))
    {
        sf->thresh_mult[THR_DUAL_ZEROLA   ] = INT_MAX;
        sf->thresh_mult[THR_DUAL_NEARESTLA] = INT_MAX;
        sf->thresh_mult[THR_DUAL_NEARLA   ] = INT_MAX;
        sf->thresh_mult[THR_DUAL_NEWLA    ] = INT_MAX;
    }

    if ((cpi->ref_frame_flags & (VP8_GOLD_FLAG | VP8_ALT_FLAG)) != (VP8_GOLD_FLAG | VP8_ALT_FLAG))
    {
        sf->thresh_mult[THR_DUAL_ZEROGA   ] = INT_MAX;
        sf->thresh_mult[THR_DUAL_NEARESTGA] = INT_MAX;
        sf->thresh_mult[THR_DUAL_NEARGA   ] = INT_MAX;
        sf->thresh_mult[THR_DUAL_NEWGA    ] = INT_MAX;
    }

    // Slow quant, dct and trellis not worthwhile for first pass
    // so make sure they are always turned off.
    if ( cpi->pass == 1 )
    {
        sf->improved_quant = 0;
        sf->optimize_coefficients = 0;
        sf->improved_dct = 0;
    }

    if (cpi->sf.search_method == NSTEP)
    {
        vp8_init3smotion_compensation(&cpi->mb, cm->yv12_fb[cm->lst_fb_idx].y_stride);
    }
    else if (cpi->sf.search_method == DIAMOND)
    {
        vp8_init_dsmotion_compensation(&cpi->mb, cm->yv12_fb[cm->lst_fb_idx].y_stride);
    }

    if (cpi->sf.improved_dct)
    {
#if CONFIG_T8X8
        cpi->mb.vp8_short_fdct8x8 = FDCT_INVOKE(&cpi->rtcd.fdct, short8x8);
#endif
        cpi->mb.vp8_short_fdct8x4 = FDCT_INVOKE(&cpi->rtcd.fdct, short8x4);
        cpi->mb.vp8_short_fdct4x4 = FDCT_INVOKE(&cpi->rtcd.fdct, short4x4);
    }
    else
    {
#if CONFIG_T8X8
        cpi->mb.vp8_short_fdct8x8 = FDCT_INVOKE(&cpi->rtcd.fdct, short8x8);
#endif
        cpi->mb.vp8_short_fdct8x4   = FDCT_INVOKE(&cpi->rtcd.fdct, fast8x4);
        cpi->mb.vp8_short_fdct4x4   = FDCT_INVOKE(&cpi->rtcd.fdct, fast4x4);
    }

    cpi->mb.short_walsh4x4 = FDCT_INVOKE(&cpi->rtcd.fdct, walsh_short4x4);
#if CONFIG_T8X8
    cpi->mb.short_fhaar2x2 = FDCT_INVOKE(&cpi->rtcd.fdct, haar_short2x2);
#endif

    if (cpi->sf.improved_quant)
    {
        cpi->mb.quantize_b      = QUANTIZE_INVOKE(&cpi->rtcd.quantize,
                                                  quantb);
        cpi->mb.quantize_b_pair = QUANTIZE_INVOKE(&cpi->rtcd.quantize,
                                                  quantb_pair);
#if CONFIG_T8X8
        cpi->mb.quantize_b_8x8  = QUANTIZE_INVOKE(&cpi->rtcd.quantize, quantb_8x8);
        cpi->mb.quantize_b_2x2  = QUANTIZE_INVOKE(&cpi->rtcd.quantize, quantb_2x2);
#endif
    }
    else
    {
        cpi->mb.quantize_b      = QUANTIZE_INVOKE(&cpi->rtcd.quantize,
                                                  fastquantb);
        cpi->mb.quantize_b_pair = QUANTIZE_INVOKE(&cpi->rtcd.quantize,
                                                  fastquantb_pair);
#if CONFIG_T8X8
        cpi->mb.quantize_b_8x8  = QUANTIZE_INVOKE(&cpi->rtcd.quantize, fastquantb_8x8);
        cpi->mb.quantize_b_2x2  = QUANTIZE_INVOKE(&cpi->rtcd.quantize, fastquantb_2x2);
#endif
    }
    if (cpi->sf.improved_quant != last_improved_quant)
        vp8cx_init_quantizer(cpi);

#if CONFIG_RUNTIME_CPU_DETECT
    cpi->mb.e_mbd.rtcd = &cpi->common.rtcd;
#endif

    if (cpi->sf.iterative_sub_pixel == 1)
    {
        cpi->find_fractional_mv_step = vp8_find_best_sub_pixel_step_iteratively;
    }
    else if (cpi->sf.quarter_pixel_search)
    {
        cpi->find_fractional_mv_step = vp8_find_best_sub_pixel_step;
    }
    else if (cpi->sf.half_pixel_search)
    {
        cpi->find_fractional_mv_step = vp8_find_best_half_pixel_step;
    }

    if (cpi->sf.optimize_coefficients == 1 && cpi->pass!=1)
        cpi->mb.optimize = 1;
    else
        cpi->mb.optimize = 0;

#ifdef SPEEDSTATS
    frames_at_speed[cpi->Speed]++;
#endif
}
static void alloc_raw_frame_buffers(VP8_COMP *cpi)
{
    int width = (cpi->oxcf.Width + 15) & ~15;
    int height = (cpi->oxcf.Height + 15) & ~15;

    cpi->lookahead = vp8_lookahead_init(cpi->oxcf.Width, cpi->oxcf.Height,
                                        cpi->oxcf.lag_in_frames);
    if(!cpi->lookahead)
        vpx_internal_error(&cpi->common.error, VPX_CODEC_MEM_ERROR,
                           "Failed to allocate lag buffers");

#if VP8_TEMPORAL_ALT_REF

    if (vp8_yv12_alloc_frame_buffer(&cpi->alt_ref_buffer,
                                    width, height, VP8BORDERINPIXELS))
        vpx_internal_error(&cpi->common.error, VPX_CODEC_MEM_ERROR,
                           "Failed to allocate altref buffer");

#endif
}

static int vp8_alloc_partition_data(VP8_COMP *cpi)
{
        vpx_free(cpi->mb.pip);

    cpi->mb.pip = vpx_calloc((cpi->common.mb_cols + 1) *
                                (cpi->common.mb_rows + 1),
                                sizeof(PARTITION_INFO));
    if(!cpi->mb.pip)
        return 1;

    cpi->mb.pi = cpi->mb.pip + cpi->common.mode_info_stride + 1;

    return 0;
}

void vp8_alloc_compressor_data(VP8_COMP *cpi)
{
    VP8_COMMON *cm = & cpi->common;

    int width = cm->Width;
    int height = cm->Height;

    if (vp8_alloc_frame_buffers(cm, width, height))
        vpx_internal_error(&cpi->common.error, VPX_CODEC_MEM_ERROR,
                           "Failed to allocate frame buffers");

    if (vp8_alloc_partition_data(cpi))
        vpx_internal_error(&cpi->common.error, VPX_CODEC_MEM_ERROR,
                           "Failed to allocate partition data");


    if ((width & 0xf) != 0)
        width += 16 - (width & 0xf);

    if ((height & 0xf) != 0)
        height += 16 - (height & 0xf);


    if (vp8_yv12_alloc_frame_buffer(&cpi->last_frame_uf,
                                    width, height, VP8BORDERINPIXELS))
        vpx_internal_error(&cpi->common.error, VPX_CODEC_MEM_ERROR,
                           "Failed to allocate last frame buffer");

    if (vp8_yv12_alloc_frame_buffer(&cpi->scaled_source,
                                    width, height, VP8BORDERINPIXELS))
        vpx_internal_error(&cpi->common.error, VPX_CODEC_MEM_ERROR,
                           "Failed to allocate scaled source buffer");


        vpx_free(cpi->tok);

    {
        unsigned int tokens = cm->mb_rows * cm->mb_cols * 24 * 16;

        CHECK_MEM_ERROR(cpi->tok, vpx_calloc(tokens, sizeof(*cpi->tok)));
    }

    // Data used for real time vc mode to see if gf needs refreshing
    cpi->inter_zz_count = 0;
    cpi->gf_bad_count = 0;
    cpi->gf_update_recommended = 0;


    // Structures used to minitor GF usage
    vpx_free(cpi->gf_active_flags);
    CHECK_MEM_ERROR(cpi->gf_active_flags,
                    vpx_calloc(1, cm->mb_rows * cm->mb_cols));
    cpi->gf_active_count = cm->mb_rows * cm->mb_cols;

    vpx_free(cpi->mb_activity_map);
    CHECK_MEM_ERROR(cpi->mb_activity_map,
                    vpx_calloc(sizeof(unsigned int),
                    cm->mb_rows * cm->mb_cols));

    vpx_free(cpi->mb_norm_activity_map);
    CHECK_MEM_ERROR(cpi->mb_norm_activity_map,
                    vpx_calloc(sizeof(unsigned int),
                    cm->mb_rows * cm->mb_cols));

        vpx_free(cpi->twopass.total_stats);

    cpi->twopass.total_stats = vpx_calloc(1, sizeof(FIRSTPASS_STATS));

    vpx_free(cpi->twopass.total_left_stats);
    cpi->twopass.total_left_stats = vpx_calloc(1, sizeof(FIRSTPASS_STATS));

        vpx_free(cpi->twopass.this_frame_stats);

    cpi->twopass.this_frame_stats = vpx_calloc(1, sizeof(FIRSTPASS_STATS));

    if( !cpi->twopass.total_stats ||
        !cpi->twopass.total_left_stats ||
        !cpi->twopass.this_frame_stats)
        vpx_internal_error(&cpi->common.error, VPX_CODEC_MEM_ERROR,
                           "Failed to allocate firstpass stats");

        vpx_free(cpi->tplist);

    CHECK_MEM_ERROR(cpi->tplist, vpx_malloc(sizeof(TOKENLIST) * cpi->common.mb_rows));
}


// TODO perhaps change number of steps expose to outside world when setting
// max and min limits. Also this will likely want refining for the extended Q
// range.
//
// Table that converts 0-63 Q range values passed in outside to the Qindex
// range used internally.
static const int q_trans[] =
{
     0,    4,   8,  12,  16,  20,  24,  28,
    32,   36,  40,  44,  48,  52,  56,  60,
    64,   68,  72,  76,  80,  84,  88,  92,
    96,  100, 104, 108, 112, 116, 120, 124,
    128, 132, 136, 140, 144, 148, 152, 156,
    160, 164, 168, 172, 176, 180, 184, 188,
    192, 196, 200, 204, 208, 212, 216, 220,
    224, 228, 232, 236, 240, 244, 249, 255,
};

int vp8_reverse_trans(int x)
{
    int i;

    for (i = 0; i < 64; i++)
        if (q_trans[i] >= x)
            return i;

    return 63;
};
void vp8_new_frame_rate(VP8_COMP *cpi, double framerate)
{
    if(framerate < .1)
        framerate = 30;

    cpi->oxcf.frame_rate             = framerate;
    cpi->output_frame_rate            = cpi->oxcf.frame_rate;
    cpi->per_frame_bandwidth          = (int)(cpi->oxcf.target_bandwidth / cpi->output_frame_rate);
    cpi->av_per_frame_bandwidth        = (int)(cpi->oxcf.target_bandwidth / cpi->output_frame_rate);
    cpi->min_frame_bandwidth          = (int)(cpi->av_per_frame_bandwidth * cpi->oxcf.two_pass_vbrmin_section / 100);

    // Set Maximum gf/arf interval
    cpi->max_gf_interval = ((int)(cpi->output_frame_rate / 2.0) + 2);

    if(cpi->max_gf_interval < 12)
        cpi->max_gf_interval = 12;

    // Extended interval for genuinely static scenes
    cpi->twopass.static_scene_max_gf_interval = cpi->key_frame_frequency >> 1;

     // Special conditions when altr ref frame enabled in lagged compress mode
    if (cpi->oxcf.play_alternate && cpi->oxcf.lag_in_frames)
    {
        if (cpi->max_gf_interval > cpi->oxcf.lag_in_frames - 1)
            cpi->max_gf_interval = cpi->oxcf.lag_in_frames - 1;

        if (cpi->twopass.static_scene_max_gf_interval > cpi->oxcf.lag_in_frames - 1)
            cpi->twopass.static_scene_max_gf_interval = cpi->oxcf.lag_in_frames - 1;
    }

    if ( cpi->max_gf_interval > cpi->twopass.static_scene_max_gf_interval )
        cpi->max_gf_interval = cpi->twopass.static_scene_max_gf_interval;
}


static int
rescale(int val, int num, int denom)
{
    int64_t llnum = num;
    int64_t llden = denom;
    int64_t llval = val;

    return llval * llnum / llden;
}


static void init_config(VP8_PTR ptr, VP8_CONFIG *oxcf)
{
    VP8_COMP *cpi = (VP8_COMP *)(ptr);
    VP8_COMMON *cm = &cpi->common;

    cpi->oxcf = *oxcf;

    cpi->goldfreq = 7;

    cm->version = oxcf->Version;
    vp8_setup_version(cm);

    // change includes all joint functionality
    vp8_change_config(ptr, oxcf);

    // Initialize active best and worst q and average q values.
    cpi->active_worst_quality         = cpi->oxcf.worst_allowed_q;
    cpi->active_best_quality          = cpi->oxcf.best_allowed_q;
    cpi->avg_frame_qindex             = cpi->oxcf.worst_allowed_q;

    // Initialise the starting buffer levels
    cpi->buffer_level                 = cpi->oxcf.starting_buffer_level;
    cpi->bits_off_target              = cpi->oxcf.starting_buffer_level;

    cpi->rolling_target_bits          = cpi->av_per_frame_bandwidth;
    cpi->rolling_actual_bits          = cpi->av_per_frame_bandwidth;
    cpi->long_rolling_target_bits     = cpi->av_per_frame_bandwidth;
    cpi->long_rolling_actual_bits     = cpi->av_per_frame_bandwidth;

    cpi->total_actual_bits            = 0;
    cpi->total_target_vs_actual       = 0;

    cpi->static_mb_pct = 0;

#if VP8_TEMPORAL_ALT_REF
    {
        int i;

        cpi->fixed_divide[0] = 0;

        for (i = 1; i < 512; i++)
            cpi->fixed_divide[i] = 0x80000 / i;
    }
#endif
}


void vp8_change_config(VP8_PTR ptr, VP8_CONFIG *oxcf)
{
    VP8_COMP *cpi = (VP8_COMP *)(ptr);
    VP8_COMMON *cm = &cpi->common;

    if (!cpi)
        return;

    if (!oxcf)
        return;

    if (cm->version != oxcf->Version)
    {
        cm->version = oxcf->Version;
        vp8_setup_version(cm);
    }

    cpi->oxcf = *oxcf;

    switch (cpi->oxcf.Mode)
    {
    // Real time and one pass deprecated in test code base
    case MODE_FIRSTPASS:
        cpi->pass = 1;
        cpi->compressor_speed = 1;
        break;

    case MODE_SECONDPASS:
        cpi->pass = 2;
        cpi->compressor_speed = 1;

        if (cpi->oxcf.cpu_used < -5)
        {
            cpi->oxcf.cpu_used = -5;
        }

        if (cpi->oxcf.cpu_used > 5)
            cpi->oxcf.cpu_used = 5;

        break;

    case MODE_SECONDPASS_BEST:
        cpi->pass = 2;
        cpi->compressor_speed = 0;
        break;
    }

    cpi->oxcf.worst_allowed_q = q_trans[oxcf->worst_allowed_q];
    cpi->oxcf.best_allowed_q = q_trans[oxcf->best_allowed_q];
    cpi->oxcf.cq_level = q_trans[cpi->oxcf.cq_level];

    cpi->baseline_gf_interval = DEFAULT_GF_INTERVAL;

    cpi->ref_frame_flags = VP8_ALT_FLAG | VP8_GOLD_FLAG | VP8_LAST_FLAG;

    //cpi->use_golden_frame_only = 0;
    //cpi->use_last_frame_only = 0;
    cm->refresh_golden_frame = 0;
    cm->refresh_last_frame = 1;
    cm->refresh_entropy_probs = 1;

    if (cpi->oxcf.token_partitions >= 0 && cpi->oxcf.token_partitions <= 3)
        cm->multi_token_partition =
            (TOKEN_PARTITION) cpi->oxcf.token_partitions;

    setup_features(cpi);
#if CONFIG_HIGH_PRECISION_MV
    cpi->mb.e_mbd.allow_high_precision_mv = 1;   // Default mv precision adaptation
#endif

    {
        int i;

        for (i = 0; i < MAX_MB_SEGMENTS; i++)
            cpi->segment_encode_breakout[i] = cpi->oxcf.encode_breakout;
    }

    // At the moment the first order values may not be > MAXQ
    if (cpi->oxcf.fixed_q > MAXQ)
        cpi->oxcf.fixed_q = MAXQ;

    // local file playback mode == really big buffer
    if (cpi->oxcf.end_usage == USAGE_LOCAL_FILE_PLAYBACK)
    {
        cpi->oxcf.starting_buffer_level   = 60000;
        cpi->oxcf.optimal_buffer_level    = 60000;
        cpi->oxcf.maximum_buffer_size     = 240000;
    }

    // Convert target bandwidth from Kbit/s to Bit/s
    cpi->oxcf.target_bandwidth       *= 1000;

    cpi->oxcf.starting_buffer_level =
        rescale(cpi->oxcf.starting_buffer_level,
                cpi->oxcf.target_bandwidth, 1000);

    // Set or reset optimal and maximum buffer levels.
    if (cpi->oxcf.optimal_buffer_level == 0)
        cpi->oxcf.optimal_buffer_level = cpi->oxcf.target_bandwidth / 8;
    else
        cpi->oxcf.optimal_buffer_level =
            rescale(cpi->oxcf.optimal_buffer_level,
                    cpi->oxcf.target_bandwidth, 1000);

    if (cpi->oxcf.maximum_buffer_size == 0)
        cpi->oxcf.maximum_buffer_size = cpi->oxcf.target_bandwidth / 8;
    else
        cpi->oxcf.maximum_buffer_size =
            rescale(cpi->oxcf.maximum_buffer_size,
                    cpi->oxcf.target_bandwidth, 1000);

    // Set up frame rate and related parameters rate control values.
    vp8_new_frame_rate(cpi, cpi->oxcf.frame_rate);

    // Set absolute upper and lower quality limits
    cpi->worst_quality               = cpi->oxcf.worst_allowed_q;
    cpi->best_quality                = cpi->oxcf.best_allowed_q;

    // active values should only be modified if out of new range
    if (cpi->active_worst_quality > cpi->oxcf.worst_allowed_q)
    {
      cpi->active_worst_quality = cpi->oxcf.worst_allowed_q;
    }
    // less likely
    else if (cpi->active_worst_quality < cpi->oxcf.best_allowed_q)
    {
      cpi->active_worst_quality = cpi->oxcf.best_allowed_q;
    }
    if (cpi->active_best_quality < cpi->oxcf.best_allowed_q)
    {
      cpi->active_best_quality = cpi->oxcf.best_allowed_q;
    }
    // less likely
    else if (cpi->active_best_quality > cpi->oxcf.worst_allowed_q)
    {
      cpi->active_best_quality = cpi->oxcf.worst_allowed_q;
    }

    cpi->buffered_mode = (cpi->oxcf.optimal_buffer_level > 0) ? TRUE : FALSE;

    cpi->cq_target_quality = cpi->oxcf.cq_level;

    if (!cm->use_bilinear_mc_filter)
        cm->mcomp_filter_type = SIXTAP;
    else
        cm->mcomp_filter_type = BILINEAR;

    cpi->target_bandwidth = cpi->oxcf.target_bandwidth;

    cm->Width       = cpi->oxcf.Width     ;
    cm->Height      = cpi->oxcf.Height    ;

    cm->horiz_scale  = cpi->horiz_scale;
    cm->vert_scale   = cpi->vert_scale ;

    // VP8 sharpness level mapping 0-7 (vs 0-10 in general VPx dialogs)
    if (cpi->oxcf.Sharpness > 7)
        cpi->oxcf.Sharpness = 7;

    cm->sharpness_level = cpi->oxcf.Sharpness;

    if (cm->horiz_scale != NORMAL || cm->vert_scale != NORMAL)
    {
        int UNINITIALIZED_IS_SAFE(hr), UNINITIALIZED_IS_SAFE(hs);
        int UNINITIALIZED_IS_SAFE(vr), UNINITIALIZED_IS_SAFE(vs);

        Scale2Ratio(cm->horiz_scale, &hr, &hs);
        Scale2Ratio(cm->vert_scale, &vr, &vs);

        // always go to the next whole number
        cm->Width = (hs - 1 + cpi->oxcf.Width * hr) / hs;
        cm->Height = (vs - 1 + cpi->oxcf.Height * vr) / vs;
    }

    if (((cm->Width + 15) & 0xfffffff0) !=
          cm->yv12_fb[cm->lst_fb_idx].y_width ||
        ((cm->Height + 15) & 0xfffffff0) !=
          cm->yv12_fb[cm->lst_fb_idx].y_height ||
        cm->yv12_fb[cm->lst_fb_idx].y_width == 0)
    {
        alloc_raw_frame_buffers(cpi);
        vp8_alloc_compressor_data(cpi);
    }

    if (cpi->oxcf.fixed_q >= 0)
    {
        cpi->last_q[0] = cpi->oxcf.fixed_q;
        cpi->last_q[1] = cpi->oxcf.fixed_q;
        cpi->last_boosted_qindex = cpi->oxcf.fixed_q;
    }

    cpi->Speed = cpi->oxcf.cpu_used;

    // force to allowlag to 0 if lag_in_frames is 0;
    if (cpi->oxcf.lag_in_frames == 0)
    {
        cpi->oxcf.allow_lag = 0;
    }
    // Limit on lag buffers as these are not currently dynamically allocated
    else if (cpi->oxcf.lag_in_frames > MAX_LAG_BUFFERS)
        cpi->oxcf.lag_in_frames = MAX_LAG_BUFFERS;

    // YX Temp
    cpi->alt_ref_source = NULL;
    cpi->is_src_frame_alt_ref = 0;


#if 0
    // Experimental RD Code
    cpi->frame_distortion = 0;
    cpi->last_frame_distortion = 0;
#endif

}

#define M_LOG2_E 0.693147180559945309417
#define log2f(x) (log (x) / (float) M_LOG2_E)
static void cal_mvsadcosts(int *mvsadcost[2])
{
    int i = 1;

    mvsadcost [0] [0] = 300;
    mvsadcost [1] [0] = 300;

    do
    {
        double z = 256 * (2 * (log2f(8 * i) + .6));
        mvsadcost [0][i] = (int) z;
        mvsadcost [1][i] = (int) z;
        mvsadcost [0][-i] = (int) z;
        mvsadcost [1][-i] = (int) z;
    }
    while (++i <= mvfp_max);
}

VP8_PTR vp8_create_compressor(VP8_CONFIG *oxcf)
{
    int i;
    volatile union
    {
        VP8_COMP *cpi;
        VP8_PTR   ptr;
    } ctx;

    VP8_COMP *cpi;
    VP8_COMMON *cm;

    cpi = ctx.cpi = vpx_memalign(32, sizeof(VP8_COMP));
    // Check that the CPI instance is valid
    if (!cpi)
        return 0;

    cm = &cpi->common;

    vpx_memset(cpi, 0, sizeof(VP8_COMP));

    if (setjmp(cm->error.jmp))
    {
        VP8_PTR ptr = ctx.ptr;

        ctx.cpi->common.error.setjmp = 0;
        vp8_remove_compressor(&ptr);
        return 0;
    }

    cpi->common.error.setjmp = 1;

    CHECK_MEM_ERROR(cpi->mb.ss, vpx_calloc(sizeof(search_site), (MAX_MVSEARCH_STEPS * 8) + 1));

    vp8_create_common(&cpi->common);
    vp8_cmachine_specific_config(cpi);

    init_config((VP8_PTR)cpi, oxcf);

    memcpy(cpi->base_skip_false_prob, vp8cx_base_skip_false_prob, sizeof(vp8cx_base_skip_false_prob));
    cpi->common.current_video_frame   = 0;
    cpi->kf_overspend_bits            = 0;
    cpi->kf_bitrate_adjustment        = 0;
    cpi->frames_till_gf_update_due      = 0;
    cpi->gf_overspend_bits            = 0;
    cpi->non_gf_bitrate_adjustment     = 0;
    cm->prob_last_coded               = 128;
    cm->prob_gf_coded                 = 128;
    cm->prob_intra_coded              = 63;
    for ( i = 0; i < DUAL_PRED_CONTEXTS; i++ )
        cm->prob_dualpred[i]         = 128;

    // Prime the recent reference frame useage counters.
    // Hereafter they will be maintained as a sort of moving average
    cpi->recent_ref_frame_usage[INTRA_FRAME]  = 1;
    cpi->recent_ref_frame_usage[LAST_FRAME]   = 1;
    cpi->recent_ref_frame_usage[GOLDEN_FRAME] = 1;
    cpi->recent_ref_frame_usage[ALTREF_FRAME] = 1;

    // Set reference frame sign bias for ALTREF frame to 1 (for now)
    cpi->common.ref_frame_sign_bias[ALTREF_FRAME] = 1;

    cpi->twopass.gf_decay_rate = 0;
    cpi->baseline_gf_interval = DEFAULT_GF_INTERVAL;

    cpi->gold_is_last = 0 ;
    cpi->alt_is_last  = 0 ;
    cpi->gold_is_alt  = 0 ;

    // allocate memory for storing last frame's MVs for MV prediction.
    CHECK_MEM_ERROR(cpi->lfmv, vpx_calloc((cpi->common.mb_rows+2) * (cpi->common.mb_cols+2), sizeof(int_mv)));
    CHECK_MEM_ERROR(cpi->lf_ref_frame_sign_bias, vpx_calloc((cpi->common.mb_rows+2) * (cpi->common.mb_cols+2), sizeof(int)));
    CHECK_MEM_ERROR(cpi->lf_ref_frame, vpx_calloc((cpi->common.mb_rows+2) * (cpi->common.mb_cols+2), sizeof(int)));

    // Create the encoder segmentation map and set all entries to 0
    CHECK_MEM_ERROR(cpi->segmentation_map, vpx_calloc((cpi->common.mb_rows * cpi->common.mb_cols), 1));

    // And a copy in common for temporal coding
    CHECK_MEM_ERROR(cm->last_frame_seg_map,
        vpx_calloc((cpi->common.mb_rows * cpi->common.mb_cols), 1));

    CHECK_MEM_ERROR(cpi->active_map, vpx_calloc(cpi->common.mb_rows * cpi->common.mb_cols, 1));
    vpx_memset(cpi->active_map , 1, (cpi->common.mb_rows * cpi->common.mb_cols));
    cpi->active_map_enabled = 0;

    for (i = 0; i < ( sizeof(cpi->mbgraph_stats) /
                      sizeof(cpi->mbgraph_stats[0]) ); i++)
    {
        CHECK_MEM_ERROR(cpi->mbgraph_stats[i].mb_stats,
                        vpx_calloc(cpi->common.mb_rows * cpi->common.mb_cols *
                                   sizeof(*cpi->mbgraph_stats[i].mb_stats),
                                   1));
    }

#ifdef ENTROPY_STATS
    init_context_counters();
#endif

    /*Initialize the feed-forward activity masking.*/
    cpi->activity_avg = 90<<12;

    cpi->frames_since_key = 8;        // Give a sensible default for the first frame.
    cpi->key_frame_frequency = cpi->oxcf.key_freq;
    cpi->this_key_frame_forced = FALSE;
    cpi->next_key_frame_forced = FALSE;

    cpi->source_alt_ref_pending = FALSE;
    cpi->source_alt_ref_active = FALSE;
    cpi->common.refresh_alt_ref_frame = 0;

    cpi->b_calculate_psnr = CONFIG_INTERNAL_STATS;
#if CONFIG_INTERNAL_STATS
    cpi->b_calculate_ssimg = 0;

    cpi->count = 0;
    cpi->bytes = 0;

    if (cpi->b_calculate_psnr)
    {
        cpi->total_sq_error = 0.0;
        cpi->total_sq_error2 = 0.0;
        cpi->total_y = 0.0;
        cpi->total_u = 0.0;
        cpi->total_v = 0.0;
        cpi->total = 0.0;
        cpi->totalp_y = 0.0;
        cpi->totalp_u = 0.0;
        cpi->totalp_v = 0.0;
        cpi->totalp = 0.0;
        cpi->tot_recode_hits = 0;
        cpi->summed_quality = 0;
        cpi->summed_weights = 0;
    }

    if (cpi->b_calculate_ssimg)
    {
        cpi->total_ssimg_y = 0;
        cpi->total_ssimg_u = 0;
        cpi->total_ssimg_v = 0;
        cpi->total_ssimg_all = 0;
    }

#endif

#ifndef LLONG_MAX
#define LLONG_MAX  9223372036854775807LL
#endif
    cpi->first_time_stamp_ever = LLONG_MAX;

    cpi->frames_till_gf_update_due      = 0;
    cpi->key_frame_count              = 1;

    cpi->ni_av_qi                     = cpi->oxcf.worst_allowed_q;
    cpi->ni_tot_qi                    = 0;
    cpi->ni_frames                   = 0;
    cpi->tot_q = 0.0;
    cpi->avg_q = vp8_convert_qindex_to_q( cpi->oxcf.worst_allowed_q );
    cpi->total_byte_count             = 0;

    cpi->rate_correction_factor         = 1.0;
    cpi->key_frame_rate_correction_factor = 1.0;
    cpi->gf_rate_correction_factor  = 1.0;
    cpi->twopass.est_max_qcorrection_factor  = 1.0;

    cpi->mb.mvcost[0] = &cpi->mb.mvcosts[0][mv_max+1];
    cpi->mb.mvcost[1] = &cpi->mb.mvcosts[1][mv_max+1];
    cpi->mb.mvsadcost[0] = &cpi->mb.mvsadcosts[0][mvfp_max+1];
    cpi->mb.mvsadcost[1] = &cpi->mb.mvsadcosts[1][mvfp_max+1];

    cal_mvsadcosts(cpi->mb.mvsadcost);

    for (i = 0; i < KEY_FRAME_CONTEXT; i++)
    {
        cpi->prior_key_frame_distance[i] = (int)cpi->output_frame_rate;
    }

#ifdef OUTPUT_YUV_SRC
    yuv_file = fopen("bd.yuv", "ab");
#endif
#ifdef OUTPUT_YUV_REC
    yuv_rec_file = fopen("rec.yuv", "wb");
#endif

#if 0
    framepsnr = fopen("framepsnr.stt", "a");
    kf_list = fopen("kf_list.stt", "w");
#endif

    cpi->output_pkt_list = oxcf->output_pkt_list;

    if (cpi->pass == 1)
    {
        vp8_init_first_pass(cpi);
    }
    else if (cpi->pass == 2)
    {
        size_t packet_sz = sizeof(FIRSTPASS_STATS);
        int packets = oxcf->two_pass_stats_in.sz / packet_sz;

        cpi->twopass.stats_in_start = oxcf->two_pass_stats_in.buf;
        cpi->twopass.stats_in = cpi->twopass.stats_in_start;
        cpi->twopass.stats_in_end = (void*)((char *)cpi->twopass.stats_in
                            + (packets - 1) * packet_sz);
        vp8_init_second_pass(cpi);
    }

    vp8_set_speed_features(cpi);

    // Set starting values of RD threshold multipliers (128 = *1)
    for (i = 0; i < MAX_MODES; i++)
    {
        cpi->rd_thresh_mult[i] = 128;
    }

#ifdef ENTROPY_STATS
    init_mv_ref_counts();
#endif

    cpi->fn_ptr[BLOCK_16X16].sdf            = VARIANCE_INVOKE(&cpi->rtcd.variance, sad16x16);
    cpi->fn_ptr[BLOCK_16X16].vf             = VARIANCE_INVOKE(&cpi->rtcd.variance, var16x16);
    cpi->fn_ptr[BLOCK_16X16].svf            = VARIANCE_INVOKE(&cpi->rtcd.variance, subpixvar16x16);
    cpi->fn_ptr[BLOCK_16X16].svf_halfpix_h  = VARIANCE_INVOKE(&cpi->rtcd.variance, halfpixvar16x16_h);
    cpi->fn_ptr[BLOCK_16X16].svf_halfpix_v  = VARIANCE_INVOKE(&cpi->rtcd.variance, halfpixvar16x16_v);
    cpi->fn_ptr[BLOCK_16X16].svf_halfpix_hv = VARIANCE_INVOKE(&cpi->rtcd.variance, halfpixvar16x16_hv);
    cpi->fn_ptr[BLOCK_16X16].sdx3f          = VARIANCE_INVOKE(&cpi->rtcd.variance, sad16x16x3);
    cpi->fn_ptr[BLOCK_16X16].sdx8f          = VARIANCE_INVOKE(&cpi->rtcd.variance, sad16x16x8);
    cpi->fn_ptr[BLOCK_16X16].sdx4df         = VARIANCE_INVOKE(&cpi->rtcd.variance, sad16x16x4d);

    cpi->fn_ptr[BLOCK_16X8].sdf            = VARIANCE_INVOKE(&cpi->rtcd.variance, sad16x8);
    cpi->fn_ptr[BLOCK_16X8].vf             = VARIANCE_INVOKE(&cpi->rtcd.variance, var16x8);
    cpi->fn_ptr[BLOCK_16X8].svf            = VARIANCE_INVOKE(&cpi->rtcd.variance, subpixvar16x8);
    cpi->fn_ptr[BLOCK_16X8].svf_halfpix_h  = NULL;
    cpi->fn_ptr[BLOCK_16X8].svf_halfpix_v  = NULL;
    cpi->fn_ptr[BLOCK_16X8].svf_halfpix_hv = NULL;
    cpi->fn_ptr[BLOCK_16X8].sdx3f          = VARIANCE_INVOKE(&cpi->rtcd.variance, sad16x8x3);
    cpi->fn_ptr[BLOCK_16X8].sdx8f          = VARIANCE_INVOKE(&cpi->rtcd.variance, sad16x8x8);
    cpi->fn_ptr[BLOCK_16X8].sdx4df         = VARIANCE_INVOKE(&cpi->rtcd.variance, sad16x8x4d);

    cpi->fn_ptr[BLOCK_8X16].sdf            = VARIANCE_INVOKE(&cpi->rtcd.variance, sad8x16);
    cpi->fn_ptr[BLOCK_8X16].vf             = VARIANCE_INVOKE(&cpi->rtcd.variance, var8x16);
    cpi->fn_ptr[BLOCK_8X16].svf            = VARIANCE_INVOKE(&cpi->rtcd.variance, subpixvar8x16);
    cpi->fn_ptr[BLOCK_8X16].svf_halfpix_h  = NULL;
    cpi->fn_ptr[BLOCK_8X16].svf_halfpix_v  = NULL;
    cpi->fn_ptr[BLOCK_8X16].svf_halfpix_hv = NULL;
    cpi->fn_ptr[BLOCK_8X16].sdx3f          = VARIANCE_INVOKE(&cpi->rtcd.variance, sad8x16x3);
    cpi->fn_ptr[BLOCK_8X16].sdx8f          = VARIANCE_INVOKE(&cpi->rtcd.variance, sad8x16x8);
    cpi->fn_ptr[BLOCK_8X16].sdx4df         = VARIANCE_INVOKE(&cpi->rtcd.variance, sad8x16x4d);

    cpi->fn_ptr[BLOCK_8X8].sdf            = VARIANCE_INVOKE(&cpi->rtcd.variance, sad8x8);
    cpi->fn_ptr[BLOCK_8X8].vf             = VARIANCE_INVOKE(&cpi->rtcd.variance, var8x8);
    cpi->fn_ptr[BLOCK_8X8].svf            = VARIANCE_INVOKE(&cpi->rtcd.variance, subpixvar8x8);
    cpi->fn_ptr[BLOCK_8X8].svf_halfpix_h  = NULL;
    cpi->fn_ptr[BLOCK_8X8].svf_halfpix_v  = NULL;
    cpi->fn_ptr[BLOCK_8X8].svf_halfpix_hv = NULL;
    cpi->fn_ptr[BLOCK_8X8].sdx3f          = VARIANCE_INVOKE(&cpi->rtcd.variance, sad8x8x3);
    cpi->fn_ptr[BLOCK_8X8].sdx8f          = VARIANCE_INVOKE(&cpi->rtcd.variance, sad8x8x8);
    cpi->fn_ptr[BLOCK_8X8].sdx4df         = VARIANCE_INVOKE(&cpi->rtcd.variance, sad8x8x4d);

    cpi->fn_ptr[BLOCK_4X4].sdf            = VARIANCE_INVOKE(&cpi->rtcd.variance, sad4x4);
    cpi->fn_ptr[BLOCK_4X4].vf             = VARIANCE_INVOKE(&cpi->rtcd.variance, var4x4);
    cpi->fn_ptr[BLOCK_4X4].svf            = VARIANCE_INVOKE(&cpi->rtcd.variance, subpixvar4x4);
    cpi->fn_ptr[BLOCK_4X4].svf_halfpix_h  = NULL;
    cpi->fn_ptr[BLOCK_4X4].svf_halfpix_v  = NULL;
    cpi->fn_ptr[BLOCK_4X4].svf_halfpix_hv = NULL;
    cpi->fn_ptr[BLOCK_4X4].sdx3f          = VARIANCE_INVOKE(&cpi->rtcd.variance, sad4x4x3);
    cpi->fn_ptr[BLOCK_4X4].sdx8f          = VARIANCE_INVOKE(&cpi->rtcd.variance, sad4x4x8);
    cpi->fn_ptr[BLOCK_4X4].sdx4df         = VARIANCE_INVOKE(&cpi->rtcd.variance, sad4x4x4d);

#if ARCH_X86 || ARCH_X86_64
    cpi->fn_ptr[BLOCK_16X16].copymem        = VARIANCE_INVOKE(&cpi->rtcd.variance, copy32xn);
    cpi->fn_ptr[BLOCK_16X8].copymem        = VARIANCE_INVOKE(&cpi->rtcd.variance, copy32xn);
    cpi->fn_ptr[BLOCK_8X16].copymem        = VARIANCE_INVOKE(&cpi->rtcd.variance, copy32xn);
    cpi->fn_ptr[BLOCK_8X8].copymem        = VARIANCE_INVOKE(&cpi->rtcd.variance, copy32xn);
    cpi->fn_ptr[BLOCK_4X4].copymem        = VARIANCE_INVOKE(&cpi->rtcd.variance, copy32xn);
#endif

    cpi->full_search_sad = SEARCH_INVOKE(&cpi->rtcd.search, full_search);
    cpi->diamond_search_sad = SEARCH_INVOKE(&cpi->rtcd.search, diamond_search);
    cpi->refining_search_sad = SEARCH_INVOKE(&cpi->rtcd.search, refining_search);

    // make sure frame 1 is okay
    cpi->error_bins[0] = cpi->common.MBs;

    //vp8cx_init_quantizer() is first called here. Add check in vp8cx_frame_init_quantizer() so that vp8cx_init_quantizer is only called later
    //when needed. This will avoid unnecessary calls of vp8cx_init_quantizer() for every frame.
    vp8cx_init_quantizer(cpi);

    vp8_loop_filter_init(cm);

    cpi->common.error.setjmp = 0;

#if CONFIG_UVINTRA
    vp8_zero(cpi->y_uv_mode_count)
#endif


    return (VP8_PTR) cpi;

}


void vp8_remove_compressor(VP8_PTR *ptr)
{
    VP8_COMP *cpi = (VP8_COMP *)(*ptr);
    int i;

    if (!cpi)
        return;

    if (cpi && (cpi->common.current_video_frame > 0))
    {
        if (cpi->pass == 2)
        {
            vp8_end_second_pass(cpi);
        }

#ifdef ENTROPY_STATS
        print_context_counters();
        print_tree_update_probs();
        print_mode_context();
#endif

#if CONFIG_INTERNAL_STATS

        vp8_clear_system_state();
#if CONFIG_T8X8
        printf("\n8x8-4x4:%d-%d\n", cpi->t8x8_count, cpi->t4x4_count);
#endif
        if (cpi->pass != 1)
        {
            FILE *f = fopen("opsnr.stt", "a");
            double time_encoded = (cpi->last_end_time_stamp_seen
                                   - cpi->first_time_stamp_ever) / 10000000.000;
            double total_encode_time = (cpi->time_receive_data + cpi->time_compress_data)   / 1000.000;
            double dr = (double)cpi->bytes * (double) 8 / (double)1000  / time_encoded;
#if defined(MODE_STATS)
            print_mode_contexts(&cpi->common);
#endif
            if (cpi->b_calculate_psnr)
            {
                YV12_BUFFER_CONFIG *lst_yv12 = &cpi->common.yv12_fb[cpi->common.lst_fb_idx];
                double samples = 3.0 / 2 * cpi->count * lst_yv12->y_width * lst_yv12->y_height;
                double total_psnr = vp8_mse2psnr(samples, 255.0, cpi->total_sq_error);
                double total_psnr2 = vp8_mse2psnr(samples, 255.0, cpi->total_sq_error2);
                double total_ssim = 100 * pow(cpi->summed_quality / cpi->summed_weights, 8.0);

                fprintf(f, "Bitrate\tAVGPsnr\tGLBPsnr\tAVPsnrP\tGLPsnrP\tVPXSSIM\t  Time(us)\n");
                fprintf(f, "%7.3f\t%7.3f\t%7.3f\t%7.3f\t%7.3f\t%7.3f\t%8.0f\n",
                        dr, cpi->total / cpi->count, total_psnr, cpi->totalp / cpi->count, total_psnr2, total_ssim,
                        total_encode_time);
//                fprintf(f, "%7.3f\t%7.3f\t%7.3f\t%7.3f\t%7.3f\t%7.3f\t%8.0f %10ld\n",
//                        dr, cpi->total / cpi->count, total_psnr, cpi->totalp / cpi->count, total_psnr2, total_ssim,
//                        total_encode_time, cpi->tot_recode_hits);
            }

            if (cpi->b_calculate_ssimg)
            {
                fprintf(f, "BitRate\tSSIM_Y\tSSIM_U\tSSIM_V\tSSIM_A\t  Time(us)\n");
                fprintf(f, "%7.3f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%8.0f\n", dr,
                        cpi->total_ssimg_y / cpi->count, cpi->total_ssimg_u / cpi->count,
                        cpi->total_ssimg_v / cpi->count, cpi->total_ssimg_all / cpi->count, total_encode_time);
//                fprintf(f, "%7.3f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%8.0f  %10ld\n", dr,
//                        cpi->total_ssimg_y / cpi->count, cpi->total_ssimg_u / cpi->count,
//                        cpi->total_ssimg_v / cpi->count, cpi->total_ssimg_all / cpi->count, total_encode_time, cpi->tot_recode_hits);
            }

            fclose(f);
        }

#endif


#ifdef MODE_STATS
        {
            extern int count_mb_seg[4];
            char modes_stats_file[250];
            FILE *f;
            double dr = (double)cpi->oxcf.frame_rate * (double)cpi->bytes * (double)8 / (double)cpi->count / (double)1000 ;
            sprintf(modes_stats_file, "modes_q%03d.stt",cpi->common.base_qindex);
            f = fopen(modes_stats_file, "w");
            fprintf(f, "intra_mode in Intra Frames:\n");
            fprintf(f, "Y: %8d, %8d, %8d, %8d, %8d, %8d\n", y_modes[0], y_modes[1], y_modes[2], y_modes[3], y_modes[4], y_modes[5]);
            fprintf(f, "I8:%8d, %8d, %8d, %8d\n", i8x8_modes[0], i8x8_modes[1], i8x8_modes[2], i8x8_modes[3]);
            fprintf(f, "UV:%8d, %8d, %8d, %8d\n", uv_modes[0], uv_modes[1], uv_modes[2], uv_modes[3]);
            fprintf(f, "KeyFrame Y-UV:\n");
            {
                int i;
                for(i=0;i<VP8_YMODES;i++)
                {
                    fprintf(f, "%2d:%8d, %8d, %8d, %8d\n",i,uv_modes_y[i][0],
                        uv_modes_y[i][1], uv_modes_y[i][2], uv_modes_y[i][3]);
                }
            }
#if CONFIG_UVINTRA
            fprintf(f, "Inter Y-UV:\n");
            {
                int i;
                for(i=0;i<VP8_YMODES;i++)
                {
                    fprintf(f, "%2d:%8d, %8d, %8d, %8d\n",i,cpi->y_uv_mode_count[i][0],
                        cpi->y_uv_mode_count[i][1], cpi->y_uv_mode_count[i][2], cpi->y_uv_mode_count[i][3]);
                }
            }
#endif
            fprintf(f, "B: ");
            {
                int i;

                for (i = 0; i < 10; i++)
                    fprintf(f, "%8d, ", b_modes[i]);

                fprintf(f, "\n");

            }

            fprintf(f, "Modes in Inter Frames:\n");
            fprintf(f,
                "Y: %8d, %8d, %8d, %8d, %8d, %8d, %8d, %8d, %8d, %8d, %8d\n",
                    inter_y_modes[0], inter_y_modes[1], inter_y_modes[2],
                    inter_y_modes[3], inter_y_modes[4], inter_y_modes[5],
                    inter_y_modes[6], inter_y_modes[7], inter_y_modes[8],
                    inter_y_modes[9], inter_y_modes[10]);
            fprintf(f, "UV:%8d, %8d, %8d, %8d\n", inter_uv_modes[0],
                    inter_uv_modes[1], inter_uv_modes[2], inter_uv_modes[3]);
            fprintf(f, "B: ");
            {
                int i;

                for (i = 0; i < 15; i++)
                    fprintf(f, "%8d, ", inter_b_modes[i]);

                fprintf(f, "\n");

            }
            fprintf(f, "P:%8d, %8d, %8d, %8d\n", count_mb_seg[0], count_mb_seg[1], count_mb_seg[2], count_mb_seg[3]);
            fprintf(f, "PB:%8d, %8d, %8d, %8d\n", inter_b_modes[LEFT4X4], inter_b_modes[ABOVE4X4], inter_b_modes[ZERO4X4], inter_b_modes[NEW4X4]);
            fclose(f);
        }
#endif

#ifdef ENTROPY_STATS
        {
            int i, j, k;
            FILE *fmode = fopen("modecontext.c", "w");

            fprintf(fmode, "\n#include \"entropymode.h\"\n\n");
            fprintf(fmode, "const unsigned int vp8_kf_default_bmode_counts ");
            fprintf(fmode, "[VP8_BINTRAMODES] [VP8_BINTRAMODES] [VP8_BINTRAMODES] =\n{\n");

            for (i = 0; i < 10; i++)
            {

                fprintf(fmode, "    { //Above Mode :  %d\n", i);

                for (j = 0; j < 10; j++)
                {

                    fprintf(fmode, "        {");

                    for (k = 0; k < 10; k++)
                    {
                        if (!intra_mode_stats[i][j][k])
                            fprintf(fmode, " %5d, ", 1);
                        else
                            fprintf(fmode, " %5d, ", intra_mode_stats[i][j][k]);
                    }

                    fprintf(fmode, "}, // left_mode %d\n", j);

                }

                fprintf(fmode, "    },\n");

            }

            fprintf(fmode, "};\n");
            fclose(fmode);
        }
#endif


#if defined(SECTIONBITS_OUTPUT)

        if (0)
        {
            int i;
            FILE *f = fopen("tokenbits.stt", "a");

            for (i = 0; i < 28; i++)
                fprintf(f, "%8d", (int)(Sectionbits[i] / 256));

            fprintf(f, "\n");
            fclose(f);
        }

#endif

#if 0
        {
            printf("\n_pick_loop_filter_level:%d\n", cpi->time_pick_lpf / 1000);
            printf("\n_frames recive_data encod_mb_row compress_frame  Total\n");
            printf("%6d %10ld %10ld %10ld %10ld\n", cpi->common.current_video_frame, cpi->time_receive_data / 1000, cpi->time_encode_mb_row / 1000, cpi->time_compress_data / 1000, (cpi->time_receive_data + cpi->time_compress_data) / 1000);
        }
#endif

    }

    dealloc_compressor_data(cpi);
    vpx_free(cpi->mb.ss);
    vpx_free(cpi->tok);

    for (i = 0; i < sizeof(cpi->mbgraph_stats) / sizeof(cpi->mbgraph_stats[0]); i++)
    {
        vpx_free(cpi->mbgraph_stats[i].mb_stats);
    }

    vp8_remove_common(&cpi->common);
    vpx_free(cpi);
    *ptr = 0;

#ifdef OUTPUT_YUV_SRC
    fclose(yuv_file);
#endif
#ifdef OUTPUT_YUV_REC
    fclose(yuv_rec_file);
#endif

#if 0

    if (keyfile)
        fclose(keyfile);

    if (framepsnr)
        fclose(framepsnr);

    if (kf_list)
        fclose(kf_list);

#endif

}


static uint64_t calc_plane_error(unsigned char *orig, int orig_stride,
                                 unsigned char *recon, int recon_stride,
                                 unsigned int cols, unsigned int rows,
                                 vp8_variance_rtcd_vtable_t *rtcd)
{
    unsigned int row, col;
    uint64_t total_sse = 0;
    int diff;

    for (row = 0; row + 16 <= rows; row += 16)
    {
        for (col = 0; col + 16 <= cols; col += 16)
        {
            unsigned int sse;

            VARIANCE_INVOKE(rtcd, mse16x16)(orig + col, orig_stride,
                                            recon + col, recon_stride,
                                            &sse);
            total_sse += sse;
        }

        /* Handle odd-sized width */
        if (col < cols)
        {
            unsigned int   border_row, border_col;
            unsigned char *border_orig = orig;
            unsigned char *border_recon = recon;

            for (border_row = 0; border_row < 16; border_row++)
            {
                for (border_col = col; border_col < cols; border_col++)
                {
                    diff = border_orig[border_col] - border_recon[border_col];
                    total_sse += diff * diff;
                }

                border_orig += orig_stride;
                border_recon += recon_stride;
            }
        }

        orig += orig_stride * 16;
        recon += recon_stride * 16;
    }

    /* Handle odd-sized height */
    for (; row < rows; row++)
    {
        for (col = 0; col < cols; col++)
        {
            diff = orig[col] - recon[col];
            total_sse += diff * diff;
        }

        orig += orig_stride;
        recon += recon_stride;
    }

    return total_sse;
}


static void generate_psnr_packet(VP8_COMP *cpi)
{
    YV12_BUFFER_CONFIG      *orig = cpi->Source;
    YV12_BUFFER_CONFIG      *recon = cpi->common.frame_to_show;
    struct vpx_codec_cx_pkt  pkt;
    uint64_t                 sse;
    int                      i;
    unsigned int             width = cpi->common.Width;
    unsigned int             height = cpi->common.Height;

    pkt.kind = VPX_CODEC_PSNR_PKT;
    sse = calc_plane_error(orig->y_buffer, orig->y_stride,
                           recon->y_buffer, recon->y_stride,
                           width, height,
                           IF_RTCD(&cpi->rtcd.variance));
    pkt.data.psnr.sse[0] = sse;
    pkt.data.psnr.sse[1] = sse;
    pkt.data.psnr.samples[0] = width * height;
    pkt.data.psnr.samples[1] = width * height;

    width = (width + 1) / 2;
    height = (height + 1) / 2;

    sse = calc_plane_error(orig->u_buffer, orig->uv_stride,
                           recon->u_buffer, recon->uv_stride,
                           width, height,
                           IF_RTCD(&cpi->rtcd.variance));
    pkt.data.psnr.sse[0] += sse;
    pkt.data.psnr.sse[2] = sse;
    pkt.data.psnr.samples[0] += width * height;
    pkt.data.psnr.samples[2] = width * height;

    sse = calc_plane_error(orig->v_buffer, orig->uv_stride,
                           recon->v_buffer, recon->uv_stride,
                           width, height,
                           IF_RTCD(&cpi->rtcd.variance));
    pkt.data.psnr.sse[0] += sse;
    pkt.data.psnr.sse[3] = sse;
    pkt.data.psnr.samples[0] += width * height;
    pkt.data.psnr.samples[3] = width * height;

    for (i = 0; i < 4; i++)
        pkt.data.psnr.psnr[i] = vp8_mse2psnr(pkt.data.psnr.samples[i], 255.0,
                                             pkt.data.psnr.sse[i]);

    vpx_codec_pkt_list_add(cpi->output_pkt_list, &pkt);
}


int vp8_use_as_reference(VP8_PTR ptr, int ref_frame_flags)
{
    VP8_COMP *cpi = (VP8_COMP *)(ptr);

    if (ref_frame_flags > 7)
        return -1 ;

    cpi->ref_frame_flags = ref_frame_flags;
    return 0;
}
int vp8_update_reference(VP8_PTR ptr, int ref_frame_flags)
{
    VP8_COMP *cpi = (VP8_COMP *)(ptr);

    if (ref_frame_flags > 7)
        return -1 ;

    cpi->common.refresh_golden_frame = 0;
    cpi->common.refresh_alt_ref_frame = 0;
    cpi->common.refresh_last_frame   = 0;

    if (ref_frame_flags & VP8_LAST_FLAG)
        cpi->common.refresh_last_frame = 1;

    if (ref_frame_flags & VP8_GOLD_FLAG)
        cpi->common.refresh_golden_frame = 1;

    if (ref_frame_flags & VP8_ALT_FLAG)
        cpi->common.refresh_alt_ref_frame = 1;

    return 0;
}

int vp8_get_reference(VP8_PTR ptr, VP8_REFFRAME ref_frame_flag, YV12_BUFFER_CONFIG *sd)
{
    VP8_COMP *cpi = (VP8_COMP *)(ptr);
    VP8_COMMON *cm = &cpi->common;
    int ref_fb_idx;

    if (ref_frame_flag == VP8_LAST_FLAG)
        ref_fb_idx = cm->lst_fb_idx;
    else if (ref_frame_flag == VP8_GOLD_FLAG)
        ref_fb_idx = cm->gld_fb_idx;
    else if (ref_frame_flag == VP8_ALT_FLAG)
        ref_fb_idx = cm->alt_fb_idx;
    else
        return -1;

    vp8_yv12_copy_frame_ptr(&cm->yv12_fb[ref_fb_idx], sd);

    return 0;
}
int vp8_set_reference(VP8_PTR ptr, VP8_REFFRAME ref_frame_flag, YV12_BUFFER_CONFIG *sd)
{
    VP8_COMP *cpi = (VP8_COMP *)(ptr);
    VP8_COMMON *cm = &cpi->common;

    int ref_fb_idx;

    if (ref_frame_flag == VP8_LAST_FLAG)
        ref_fb_idx = cm->lst_fb_idx;
    else if (ref_frame_flag == VP8_GOLD_FLAG)
        ref_fb_idx = cm->gld_fb_idx;
    else if (ref_frame_flag == VP8_ALT_FLAG)
        ref_fb_idx = cm->alt_fb_idx;
    else
        return -1;

    vp8_yv12_copy_frame_ptr(sd, &cm->yv12_fb[ref_fb_idx]);

    return 0;
}
int vp8_update_entropy(VP8_PTR comp, int update)
{
    VP8_COMP *cpi = (VP8_COMP *) comp;
    VP8_COMMON *cm = &cpi->common;
    cm->refresh_entropy_probs = update;

    return 0;
}


#ifdef OUTPUT_YUV_SRC
void vp8_write_yuv_frame(YV12_BUFFER_CONFIG *s)
{
    unsigned char *src = s->y_buffer;
    int h = s->y_height;

    do
    {
        fwrite(src, s->y_width, 1,  yuv_file);
        src += s->y_stride;
    }
    while (--h);

    src = s->u_buffer;
    h = s->uv_height;

    do
    {
        fwrite(src, s->uv_width, 1,  yuv_file);
        src += s->uv_stride;
    }
    while (--h);

    src = s->v_buffer;
    h = s->uv_height;

    do
    {
        fwrite(src, s->uv_width, 1, yuv_file);
        src += s->uv_stride;
    }
    while (--h);
}
#endif

#ifdef OUTPUT_YUV_REC
void vp8_write_yuv_rec_frame(VP8_COMMON *cm)
{
    YV12_BUFFER_CONFIG *s = cm->frame_to_show;
    unsigned char *src = s->y_buffer;
    int h = cm->Height;

    do
    {
        fwrite(src, s->y_width, 1,  yuv_rec_file);
        src += s->y_stride;
    }
    while (--h);

    src = s->u_buffer;
    h = (cm->Height+1)/2;

    do
    {
        fwrite(src, s->uv_width, 1,  yuv_rec_file);
        src += s->uv_stride;
    }
    while (--h);

    src = s->v_buffer;
    h = (cm->Height+1)/2;

    do
    {
        fwrite(src, s->uv_width, 1, yuv_rec_file);
        src += s->uv_stride;
    }
    while (--h);
}
#endif

static void update_alt_ref_frame_stats(VP8_COMP *cpi)
{
    VP8_COMMON *cm = &cpi->common;

    // Update data structure that monitors level of reference to last GF
    vpx_memset(cpi->gf_active_flags, 1, (cm->mb_rows * cm->mb_cols));
    cpi->gf_active_count = cm->mb_rows * cm->mb_cols;

    // this frame refreshes means next frames don't unless specified by user
    cpi->common.frames_since_golden = 0;

    // Clear the alternate reference update pending flag.
    cpi->source_alt_ref_pending = FALSE;

    // Set the alternate refernce frame active flag
    cpi->source_alt_ref_active = TRUE;


}
static void update_golden_frame_stats(VP8_COMP *cpi)
{
    VP8_COMMON *cm = &cpi->common;

    // Update the Golden frame usage counts.
    if (cm->refresh_golden_frame)
    {
        // Update data structure that monitors level of reference to last GF
        vpx_memset(cpi->gf_active_flags, 1, (cm->mb_rows * cm->mb_cols));
        cpi->gf_active_count = cm->mb_rows * cm->mb_cols;

        // this frame refreshes means next frames don't unless specified by user
        cm->refresh_golden_frame = 0;
        cpi->common.frames_since_golden = 0;

        //if ( cm->frame_type == KEY_FRAME )
        //{
        cpi->recent_ref_frame_usage[INTRA_FRAME] = 1;
        cpi->recent_ref_frame_usage[LAST_FRAME] = 1;
        cpi->recent_ref_frame_usage[GOLDEN_FRAME] = 1;
        cpi->recent_ref_frame_usage[ALTREF_FRAME] = 1;
        //}
        //else
        //{
        //  // Carry a potrtion of count over to begining of next gf sequence
        //  cpi->recent_ref_frame_usage[INTRA_FRAME] >>= 5;
        //  cpi->recent_ref_frame_usage[LAST_FRAME] >>= 5;
        //  cpi->recent_ref_frame_usage[GOLDEN_FRAME] >>= 5;
        //  cpi->recent_ref_frame_usage[ALTREF_FRAME] >>= 5;
        //}

        // ******** Fixed Q test code only ************
        // If we are going to use the ALT reference for the next group of frames set a flag to say so.
        if (cpi->oxcf.fixed_q >= 0 &&
            cpi->oxcf.play_alternate && !cpi->common.refresh_alt_ref_frame)
        {
            cpi->source_alt_ref_pending = TRUE;
            cpi->frames_till_gf_update_due = cpi->baseline_gf_interval;
        }

        if (!cpi->source_alt_ref_pending)
            cpi->source_alt_ref_active = FALSE;

        // Decrement count down till next gf
        if (cpi->frames_till_gf_update_due > 0)
            cpi->frames_till_gf_update_due--;

    }
    else if (!cpi->common.refresh_alt_ref_frame)
    {
        // Decrement count down till next gf
        if (cpi->frames_till_gf_update_due > 0)
            cpi->frames_till_gf_update_due--;

        if (cpi->common.frames_till_alt_ref_frame)
            cpi->common.frames_till_alt_ref_frame --;

        cpi->common.frames_since_golden ++;

        if (cpi->common.frames_since_golden > 1)
        {
            cpi->recent_ref_frame_usage[INTRA_FRAME] += cpi->count_mb_ref_frame_usage[INTRA_FRAME];
            cpi->recent_ref_frame_usage[LAST_FRAME] += cpi->count_mb_ref_frame_usage[LAST_FRAME];
            cpi->recent_ref_frame_usage[GOLDEN_FRAME] += cpi->count_mb_ref_frame_usage[GOLDEN_FRAME];
            cpi->recent_ref_frame_usage[ALTREF_FRAME] += cpi->count_mb_ref_frame_usage[ALTREF_FRAME];
        }
    }
}

// 1 = key, 0 = inter
static int decide_key_frame(VP8_COMP *cpi)
{
    VP8_COMMON *cm = &cpi->common;

    int code_key_frame = FALSE;

    cpi->kf_boost = 0;

    if (cpi->Speed > 11)
        return FALSE;

    // Clear down mmx registers
    vp8_clear_system_state();  //__asm emms;

    // If the following are true we might as well code a key frame
    if (((cpi->this_frame_percent_intra == 100) &&
         (cpi->this_frame_percent_intra > (cpi->last_frame_percent_intra + 2))) ||
        ((cpi->this_frame_percent_intra > 95) &&
         (cpi->this_frame_percent_intra >= (cpi->last_frame_percent_intra + 5))))
    {
        code_key_frame = TRUE;
    }
    // in addition if the following are true and this is not a golden frame then code a key frame
    // Note that on golden frames there often seems to be a pop in intra useage anyway hence this
    // restriction is designed to prevent spurious key frames. The Intra pop needs to be investigated.
    else if (((cpi->this_frame_percent_intra > 60) &&
              (cpi->this_frame_percent_intra > (cpi->last_frame_percent_intra * 2))) ||
             ((cpi->this_frame_percent_intra > 75) &&
              (cpi->this_frame_percent_intra > (cpi->last_frame_percent_intra * 3 / 2))) ||
             ((cpi->this_frame_percent_intra > 90) &&
              (cpi->this_frame_percent_intra > (cpi->last_frame_percent_intra + 10))))
    {
        if (!cm->refresh_golden_frame)
            code_key_frame = TRUE;
    }

    return code_key_frame;

}

int find_fp_qindex()
{
    int i;

    for ( i = 0; i < QINDEX_RANGE; i++ )
    {
        if ( vp8_convert_qindex_to_q(i) >= 30.0 )
        {
            break;
        }
    }

    if ( i == QINDEX_RANGE )
        i--;

    return i;
}

static void Pass1Encode(VP8_COMP *cpi, unsigned long *size, unsigned char *dest, unsigned int *frame_flags)
{
    (void) size;
    (void) dest;
    (void) frame_flags;


    vp8_set_quantizer(cpi, find_fp_qindex());
    vp8_first_pass(cpi);
}
//#define WRITE_RECON_BUFFER 1
#if WRITE_RECON_BUFFER
void write_cx_frame_to_file(YV12_BUFFER_CONFIG *frame, int this_frame)
{

    // write the frame
    FILE *yframe;
    int i;
    char filename[255];

    sprintf(filename, "cx\\y%04d.raw", this_frame);
    yframe = fopen(filename, "wb");

    for (i = 0; i < frame->y_height; i++)
        fwrite(frame->y_buffer + i * frame->y_stride,
            frame->y_width, 1, yframe);

    fclose(yframe);
    sprintf(filename, "cx\\u%04d.raw", this_frame);
    yframe = fopen(filename, "wb");

    for (i = 0; i < frame->uv_height; i++)
        fwrite(frame->u_buffer + i * frame->uv_stride,
            frame->uv_width, 1, yframe);

    fclose(yframe);
    sprintf(filename, "cx\\v%04d.raw", this_frame);
    yframe = fopen(filename, "wb");

    for (i = 0; i < frame->uv_height; i++)
        fwrite(frame->v_buffer + i * frame->uv_stride,
            frame->uv_width, 1, yframe);

    fclose(yframe);
}
#endif

// Function to test for conditions that indeicate we should loop
// back and recode a frame.
static BOOL recode_loop_test( VP8_COMP *cpi,
                              int high_limit, int low_limit,
                              int q, int maxq, int minq )
{
    BOOL    force_recode = FALSE;
    VP8_COMMON *cm = &cpi->common;

    // Is frame recode allowed at all
    // Yes if either recode mode 1 is selected or mode two is selcted
    // and the frame is a key frame. golden frame or alt_ref_frame
    if ( (cpi->sf.recode_loop == 1) ||
         ( (cpi->sf.recode_loop == 2) &&
           ( (cm->frame_type == KEY_FRAME) ||
             cm->refresh_golden_frame ||
             cm->refresh_alt_ref_frame ) ) )
    {
        // General over and under shoot tests
        if ( ((cpi->projected_frame_size > high_limit) && (q < maxq)) ||
             ((cpi->projected_frame_size < low_limit) && (q > minq)) )
        {
            force_recode = TRUE;
        }
        // Special Constrained quality tests
        else if (cpi->oxcf.end_usage == USAGE_CONSTRAINED_QUALITY)
        {
            // Undershoot and below auto cq level
            if ( (q > cpi->cq_target_quality) &&
                 (cpi->projected_frame_size <
                     ((cpi->this_frame_target * 7) >> 3)))
            {
                force_recode = TRUE;
            }
            // Severe undershoot and between auto and user cq level
            else if ( (q > cpi->oxcf.cq_level) &&
                      (cpi->projected_frame_size < cpi->min_frame_bandwidth) &&
                      (cpi->active_best_quality > cpi->oxcf.cq_level))
            {
                force_recode = TRUE;
                cpi->active_best_quality = cpi->oxcf.cq_level;
            }
        }
    }

    return force_recode;
}

void update_reference_frames(VP8_COMMON *cm)
{
    YV12_BUFFER_CONFIG *yv12_fb = cm->yv12_fb;

    // At this point the new frame has been encoded.
    // If any buffer copy / swapping is signaled it should be done here.

    if (cm->frame_type == KEY_FRAME)
    {
        yv12_fb[cm->new_fb_idx].flags |= VP8_GOLD_FLAG | VP8_ALT_FLAG ;

        yv12_fb[cm->gld_fb_idx].flags &= ~VP8_GOLD_FLAG;
        yv12_fb[cm->alt_fb_idx].flags &= ~VP8_ALT_FLAG;

        cm->alt_fb_idx = cm->gld_fb_idx = cm->new_fb_idx;
    }
    else    /* For non key frames */
    {
        if (cm->refresh_alt_ref_frame)
        {
            assert(!cm->copy_buffer_to_arf);

            cm->yv12_fb[cm->new_fb_idx].flags |= VP8_ALT_FLAG;
            cm->yv12_fb[cm->alt_fb_idx].flags &= ~VP8_ALT_FLAG;
            cm->alt_fb_idx = cm->new_fb_idx;
        }
        else if (cm->copy_buffer_to_arf)
        {
            assert(!(cm->copy_buffer_to_arf & ~0x3));

            if (cm->copy_buffer_to_arf == 1)
            {
                if(cm->alt_fb_idx != cm->lst_fb_idx)
                {
                    yv12_fb[cm->lst_fb_idx].flags |= VP8_ALT_FLAG;
                    yv12_fb[cm->alt_fb_idx].flags &= ~VP8_ALT_FLAG;
                    cm->alt_fb_idx = cm->lst_fb_idx;
                }
            }
            else /* if (cm->copy_buffer_to_arf == 2) */
            {
                if(cm->alt_fb_idx != cm->gld_fb_idx)
                {
                    yv12_fb[cm->gld_fb_idx].flags |= VP8_ALT_FLAG;
                    yv12_fb[cm->alt_fb_idx].flags &= ~VP8_ALT_FLAG;
                    cm->alt_fb_idx = cm->gld_fb_idx;
                }
            }
        }

        if (cm->refresh_golden_frame)
        {
            assert(!cm->copy_buffer_to_gf);

            cm->yv12_fb[cm->new_fb_idx].flags |= VP8_GOLD_FLAG;
            cm->yv12_fb[cm->gld_fb_idx].flags &= ~VP8_GOLD_FLAG;
            cm->gld_fb_idx = cm->new_fb_idx;
        }
        else if (cm->copy_buffer_to_gf)
        {
            assert(!(cm->copy_buffer_to_arf & ~0x3));

            if (cm->copy_buffer_to_gf == 1)
            {
                if(cm->gld_fb_idx != cm->lst_fb_idx)
                {
                    yv12_fb[cm->lst_fb_idx].flags |= VP8_GOLD_FLAG;
                    yv12_fb[cm->gld_fb_idx].flags &= ~VP8_GOLD_FLAG;
                    cm->gld_fb_idx = cm->lst_fb_idx;
                }
            }
            else /* if (cm->copy_buffer_to_gf == 2) */
            {
                if(cm->alt_fb_idx != cm->gld_fb_idx)
                {
                    yv12_fb[cm->alt_fb_idx].flags |= VP8_GOLD_FLAG;
                    yv12_fb[cm->gld_fb_idx].flags &= ~VP8_GOLD_FLAG;
                    cm->gld_fb_idx = cm->alt_fb_idx;
                }
            }
        }
    }

    if (cm->refresh_last_frame)
    {
        cm->yv12_fb[cm->new_fb_idx].flags |= VP8_LAST_FLAG;
        cm->yv12_fb[cm->lst_fb_idx].flags &= ~VP8_LAST_FLAG;
        cm->lst_fb_idx = cm->new_fb_idx;
    }
}

void loopfilter_frame(VP8_COMP *cpi, VP8_COMMON *cm)
{
    if (cm->no_lpf)
    {
        cm->filter_level = 0;
    }
    else
    {
        struct vpx_usec_timer timer;

        vp8_clear_system_state();

        vpx_usec_timer_start(&timer);
        if (cpi->sf.auto_filter == 0)
            vp8cx_pick_filter_level_fast(cpi->Source, cpi);

        else
            vp8cx_pick_filter_level(cpi->Source, cpi);

        vpx_usec_timer_mark(&timer);
        cpi->time_pick_lpf += vpx_usec_timer_elapsed(&timer);
    }

    if (cm->filter_level > 0)
    {
        vp8cx_set_alt_lf_level(cpi, cm->filter_level);
        vp8_loop_filter_frame(cm, &cpi->mb.e_mbd);
    }

    vp8_yv12_extend_frame_borders_ptr(cm->frame_to_show);

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

static void encode_frame_to_data_rate
(
    VP8_COMP *cpi,
    unsigned long *size,
    unsigned char *dest,
    unsigned int *frame_flags
)
{
    VP8_COMMON *cm = &cpi->common;
    MACROBLOCKD *xd = &cpi->mb.e_mbd;

    int Q;
    int frame_over_shoot_limit;
    int frame_under_shoot_limit;

    int Loop = FALSE;
    int loop_count;
    int this_q;
    int last_zbin_oq;

    int q_low;
    int q_high;
    int zbin_oq_high;
    int zbin_oq_low = 0;
    int top_index;
    int bottom_index;
    int active_worst_qchanged = FALSE;

    int overshoot_seen = FALSE;
    int undershoot_seen = FALSE;

    // Clear down mmx registers to allow floating point in what follows
    vp8_clear_system_state();

    // For an alt ref frame in 2 pass we skip the call to the second
    // pass function that sets the target bandwidth so must set it here
    if (cpi->common.refresh_alt_ref_frame)
    {
        cpi->per_frame_bandwidth = cpi->twopass.gf_bits;                           // Per frame bit target for the alt ref frame
        cpi->target_bandwidth = cpi->twopass.gf_bits * cpi->output_frame_rate;      // per second target bitrate
    }

    // Default turn off buffer to buffer copying
    cm->copy_buffer_to_gf = 0;
    cm->copy_buffer_to_arf = 0;

    // Clear zbin over-quant value and mode boost values.
    cpi->zbin_over_quant = 0;
    cpi->zbin_mode_boost = 0;

    // Enable or disable mode based tweaking of the zbin
    // For 2 Pass Only used where GF/ARF prediction quality
    // is above a threshold
    cpi->zbin_mode_boost = 0;
    cpi->zbin_mode_boost_enabled = TRUE;
    if ( cpi->gfu_boost <= 400 )
    {
        cpi->zbin_mode_boost_enabled = FALSE;
    }

    // Current default encoder behaviour for the altref sign bias
    if (cpi->source_alt_ref_active)
        cpi->common.ref_frame_sign_bias[ALTREF_FRAME] = 1;
    else
        cpi->common.ref_frame_sign_bias[ALTREF_FRAME] = 0;

    // Check to see if a key frame is signalled
    // For two pass with auto key frame enabled cm->frame_type may already be set, but not for one pass.
    if ((cm->current_video_frame == 0) ||
        (cm->frame_flags & FRAMEFLAGS_KEY) ||
        (cpi->oxcf.auto_key && (cpi->frames_since_key % cpi->key_frame_frequency == 0)))
    {
        // Key frame from VFW/auto-keyframe/first frame
        cm->frame_type = KEY_FRAME;
    }

    // Set default state for segment based loop filter update flags
    xd->mode_ref_lf_delta_update = 0;

    // Set various flags etc to special state if it is a key frame
    if (cm->frame_type == KEY_FRAME)
    {
        int i;

        // Reset the loop filter deltas and segmentation map
        setup_features(cpi);
#if CONFIG_HIGH_PRECISION_MV
        xd->allow_high_precision_mv = 1;   // Default mv precision adaptation
#endif

        // If segmentation is enabled force a map update for key frames
        if (xd->segmentation_enabled)
        {
            xd->update_mb_segmentation_map = 1;
            xd->update_mb_segmentation_data = 1;
        }

        // The alternate reference frame cannot be active for a key frame
        cpi->source_alt_ref_active = FALSE;

        // Reset the RD threshold multipliers to default of * 1 (128)
        for (i = 0; i < MAX_MODES; i++)
        {
            cpi->rd_thresh_mult[i] = 128;
        }
    }

//#if !CONFIG_COMPRED
    // This function has been deprecated for now but we may want to do
    // something here at a late date
    //update_rd_ref_frame_probs(cpi);
//#endif

    // Test code for new segment features
    init_seg_features( cpi );

    // Decide how big to make the frame
    if (!vp8_pick_frame_size(cpi))
    {
        cm->current_video_frame++;
        cpi->frames_since_key++;
        return;
    }

    vp8_clear_system_state();

    // Set an active best quality and if necessary active worst quality
    Q = cpi->active_worst_quality;

    if ( cm->frame_type == KEY_FRAME )
    {
        if (cpi->gfu_boost > 600)
           cpi->active_best_quality = kf_low_motion_minq[Q];
        else
           cpi->active_best_quality = kf_high_motion_minq[Q];

        // Special case for key frames forced because we have reached
        // the maximum key frame interval. Here force the Q to a range
        // based on the ambient Q to reduce the risk of popping
        if ( cpi->this_key_frame_forced )
        {
            int delta_qindex;
            int qindex = cpi->last_boosted_qindex;

            delta_qindex = compute_qdelta( cpi, qindex,
                                           (qindex * 0.75) );

            cpi->active_best_quality = qindex + delta_qindex;
            if (cpi->active_best_quality < cpi->best_quality)
                cpi->active_best_quality = cpi->best_quality;
        }
    }

    else if (cm->refresh_golden_frame || cpi->common.refresh_alt_ref_frame)
    {
        // Use the lower of cpi->active_worst_quality and recent
        // average Q as basis for GF/ARF Q limit unless last frame was
        // a key frame.
        if ( (cpi->frames_since_key > 1) &&
             (cpi->avg_frame_qindex < cpi->active_worst_quality) )
        {
            Q = cpi->avg_frame_qindex;
        }

        // For constrained quality dont allow Q less than the cq level
        if ( (cpi->oxcf.end_usage == USAGE_CONSTRAINED_QUALITY) &&
             (Q < cpi->cq_target_quality) )
        {
            Q = cpi->cq_target_quality;
        }

        if ( cpi->gfu_boost > 1000 )
            cpi->active_best_quality = gf_low_motion_minq[Q];
        else if ( cpi->gfu_boost < 400 )
            cpi->active_best_quality = gf_high_motion_minq[Q];
        else
            cpi->active_best_quality = gf_mid_motion_minq[Q];

        // Constrained quality use slightly lower active best.
        if ( cpi->oxcf.end_usage == USAGE_CONSTRAINED_QUALITY )
        {
            cpi->active_best_quality =
                cpi->active_best_quality * 15/16;
        }
    }
    else
    {
        cpi->active_best_quality = inter_minq[Q];

        // For the constant/constrained quality mode we dont want
        // q to fall below the cq level.
        if ((cpi->oxcf.end_usage == USAGE_CONSTRAINED_QUALITY) &&
            (cpi->active_best_quality < cpi->cq_target_quality) )
        {
            // If we are strongly undershooting the target rate in the last
            // frames then use the user passed in cq value not the auto
            // cq value.
            if ( cpi->rolling_actual_bits < cpi->min_frame_bandwidth )
                cpi->active_best_quality = cpi->oxcf.cq_level;
            else
                cpi->active_best_quality = cpi->cq_target_quality;
        }
    }

    // Clip the active best and worst quality values to limits
    if (cpi->active_worst_quality > cpi->worst_quality)
        cpi->active_worst_quality = cpi->worst_quality;

    if (cpi->active_best_quality < cpi->best_quality)
        cpi->active_best_quality = cpi->best_quality;

    if (cpi->active_best_quality > cpi->worst_quality)
        cpi->active_best_quality = cpi->worst_quality;

    if ( cpi->active_worst_quality < cpi->active_best_quality )
        cpi->active_worst_quality = cpi->active_best_quality;

    // Specuial case code to try and match quality with forced key frames
    if ( (cm->frame_type == KEY_FRAME) && cpi->this_key_frame_forced )
    {
        Q = cpi->last_boosted_qindex;
    }
    else
    {
        // Determine initial Q to try
        Q = vp8_regulate_q(cpi, cpi->this_frame_target);
    }
    last_zbin_oq = cpi->zbin_over_quant;

    // Set highest allowed value for Zbin over quant
    if (cm->frame_type == KEY_FRAME)
        zbin_oq_high = 0; //ZBIN_OQ_MAX/16
    else if (cm->refresh_alt_ref_frame || (cm->refresh_golden_frame && !cpi->source_alt_ref_active))
        zbin_oq_high = 16;
    else
        zbin_oq_high = ZBIN_OQ_MAX;

    vp8_compute_frame_size_bounds(cpi, &frame_under_shoot_limit, &frame_over_shoot_limit);

    // Limit Q range for the adaptive loop.
    bottom_index = cpi->active_best_quality;
    top_index    = cpi->active_worst_quality;
    q_low  = cpi->active_best_quality;
    q_high = cpi->active_worst_quality;

    vp8_save_coding_context(cpi);

    loop_count = 0;

#if CONFIG_POSTPROC

    if (cpi->oxcf.noise_sensitivity > 0)
    {
        unsigned char *src;
        int l = 0;

        switch (cpi->oxcf.noise_sensitivity)
        {
        case 1:
            l = 20;
            break;
        case 2:
            l = 40;
            break;
        case 3:
            l = 60;
            break;
        case 4:
            l = 80;
            break;
        case 5:
            l = 100;
            break;
        case 6:
            l = 150;
            break;
        }


        if (cm->frame_type == KEY_FRAME)
        {
            vp8_de_noise(cpi->Source, cpi->Source, l , 1,  0, RTCD(postproc));
        }
        else
        {
            vp8_de_noise(cpi->Source, cpi->Source, l , 1,  0, RTCD(postproc));

            src = cpi->Source->y_buffer;

            if (cpi->Source->y_stride < 0)
            {
                src += cpi->Source->y_stride * (cpi->Source->y_height - 1);
            }
        }
    }

#endif

#ifdef OUTPUT_YUV_SRC
    vp8_write_yuv_frame(cpi->Source);
#endif

    do
    {
        vp8_clear_system_state();  //__asm emms;

        vp8_set_quantizer(cpi, Q);
        this_q = Q;

        // setup skip prob for costing in mode/mv decision
        if (cpi->common.mb_no_coeff_skip)
        {
            cpi->prob_skip_false = cpi->base_skip_false_prob[Q];

            if (cm->frame_type != KEY_FRAME)
            {
                if (cpi->common.refresh_alt_ref_frame)
                {
                    if (cpi->last_skip_false_probs[2] != 0)
                        cpi->prob_skip_false = cpi->last_skip_false_probs[2];

                    /*
                                        if(cpi->last_skip_false_probs[2]!=0 && abs(Q- cpi->last_skip_probs_q[2])<=16 )
                       cpi->prob_skip_false = cpi->last_skip_false_probs[2];
                                        else if (cpi->last_skip_false_probs[2]!=0)
                       cpi->prob_skip_false = (cpi->last_skip_false_probs[2]  + cpi->prob_skip_false ) / 2;
                       */
                }
                else if (cpi->common.refresh_golden_frame)
                {
                    if (cpi->last_skip_false_probs[1] != 0)
                        cpi->prob_skip_false = cpi->last_skip_false_probs[1];

                    /*
                                        if(cpi->last_skip_false_probs[1]!=0 && abs(Q- cpi->last_skip_probs_q[1])<=16 )
                       cpi->prob_skip_false = cpi->last_skip_false_probs[1];
                                        else if (cpi->last_skip_false_probs[1]!=0)
                       cpi->prob_skip_false = (cpi->last_skip_false_probs[1]  + cpi->prob_skip_false ) / 2;
                       */
                }
                else
                {
                    if (cpi->last_skip_false_probs[0] != 0)
                        cpi->prob_skip_false = cpi->last_skip_false_probs[0];

                    /*
                    if(cpi->last_skip_false_probs[0]!=0 && abs(Q- cpi->last_skip_probs_q[0])<=16 )
                        cpi->prob_skip_false = cpi->last_skip_false_probs[0];
                    else if(cpi->last_skip_false_probs[0]!=0)
                        cpi->prob_skip_false = (cpi->last_skip_false_probs[0]  + cpi->prob_skip_false ) / 2;
                        */
                }

                // as this is for cost estimate, let's make sure it does not
                // get extreme either way
                if (cpi->prob_skip_false < 5)
                    cpi->prob_skip_false = 5;

                if (cpi->prob_skip_false > 250)
                    cpi->prob_skip_false = 250;

                if (cpi->is_src_frame_alt_ref)
                    cpi->prob_skip_false = 1;


            }
        }

        // Set up entropy depending on frame type.
        if (cm->frame_type == KEY_FRAME)
            vp8_setup_key_frame(cpi);
        else
            vp8_setup_inter_frame(cpi);

        // transform / motion compensation build reconstruction frame
        vp8_encode_frame(cpi);

        cpi->projected_frame_size -= vp8_estimate_entropy_savings(cpi);
        cpi->projected_frame_size = (cpi->projected_frame_size > 0) ? cpi->projected_frame_size : 0;

        vp8_clear_system_state();  //__asm emms;

        if (frame_over_shoot_limit == 0)
            frame_over_shoot_limit = 1;

        active_worst_qchanged = FALSE;

        // Special case handling for forced key frames
        if ( (cm->frame_type == KEY_FRAME) && cpi->this_key_frame_forced )
        {
            int last_q = Q;
            int kf_err = vp8_calc_ss_err(cpi->Source,
                                         &cm->yv12_fb[cm->new_fb_idx],
                                         IF_RTCD(&cpi->rtcd.variance));

            int high_err_target = cpi->ambient_err;
            int low_err_target = ((cpi->ambient_err * 3) >> 2);

            // Prevent possible divide by zero error below for perfect KF
            kf_err += (!kf_err);

            // The key frame is not good enough
            if ( (kf_err > high_err_target) &&
                 (cpi->projected_frame_size <= frame_over_shoot_limit) )
            {
                // Lower q_high
                q_high = (Q > q_low) ? (Q - 1) : q_low;

                // Adjust Q
                Q = (Q * high_err_target) / kf_err;
                if ( Q < ((q_high + q_low) >> 1))
                    Q = (q_high + q_low) >> 1;
            }
            // The key frame is much better than the previous frame
            else if ( (kf_err < low_err_target) &&
                      (cpi->projected_frame_size >= frame_under_shoot_limit) )
            {
                // Raise q_low
                q_low = (Q < q_high) ? (Q + 1) : q_high;

                // Adjust Q
                Q = (Q * low_err_target) / kf_err;
                if ( Q > ((q_high + q_low + 1) >> 1))
                    Q = (q_high + q_low + 1) >> 1;
            }

            // Clamp Q to upper and lower limits:
            if (Q > q_high)
                Q = q_high;
            else if (Q < q_low)
                Q = q_low;

            Loop = ((Q != last_q)) ? TRUE : FALSE;
        }

        // Is the projected frame size out of range and are we allowed to attempt to recode.
        else if ( recode_loop_test( cpi,
                               frame_over_shoot_limit, frame_under_shoot_limit,
                               Q, top_index, bottom_index ) )
        {
            int last_q = Q;
            int Retries = 0;

            // Frame size out of permitted range:
            // Update correction factor & compute new Q to try...

            // Frame is too large
            if (cpi->projected_frame_size > cpi->this_frame_target)
            {
                q_low = (Q < q_high) ? (Q + 1) : q_high; // Raise Qlow as to at least the current value

                if (cpi->zbin_over_quant > 0)            // If we are using over quant do the same for zbin_oq_low
                    zbin_oq_low = (cpi->zbin_over_quant < zbin_oq_high) ? (cpi->zbin_over_quant + 1) : zbin_oq_high;

                if ( undershoot_seen || (loop_count > 1) )
                {
                    // Update rate_correction_factor unless cpi->active_worst_quality has changed.
                    if (!active_worst_qchanged)
                        vp8_update_rate_correction_factors(cpi, 1);

                    Q = (q_high + q_low + 1) / 2;

                    // Adjust cpi->zbin_over_quant (only allowed when Q is max)
                    if (Q < MAXQ)
                        cpi->zbin_over_quant = 0;
                    else
                    {
                        zbin_oq_low = (cpi->zbin_over_quant < zbin_oq_high) ? (cpi->zbin_over_quant + 1) : zbin_oq_high;
                        cpi->zbin_over_quant = (zbin_oq_high + zbin_oq_low) / 2;
                    }
                }
                else
                {
                    // Update rate_correction_factor unless cpi->active_worst_quality has changed.
                    if (!active_worst_qchanged)
                        vp8_update_rate_correction_factors(cpi, 0);

                    Q = vp8_regulate_q(cpi, cpi->this_frame_target);

                    while (((Q < q_low) || (cpi->zbin_over_quant < zbin_oq_low)) && (Retries < 10))
                    {
                        vp8_update_rate_correction_factors(cpi, 0);
                        Q = vp8_regulate_q(cpi, cpi->this_frame_target);
                        Retries ++;
                    }
                }

                overshoot_seen = TRUE;
            }
            // Frame is too small
            else
            {
                if (cpi->zbin_over_quant == 0)
                    q_high = (Q > q_low) ? (Q - 1) : q_low; // Lower q_high if not using over quant
                else                                    // else lower zbin_oq_high
                    zbin_oq_high = (cpi->zbin_over_quant > zbin_oq_low) ? (cpi->zbin_over_quant - 1) : zbin_oq_low;

                if ( overshoot_seen || (loop_count > 1) )
                {
                    // Update rate_correction_factor unless cpi->active_worst_quality has changed.
                    if (!active_worst_qchanged)
                        vp8_update_rate_correction_factors(cpi, 1);

                    Q = (q_high + q_low) / 2;

                    // Adjust cpi->zbin_over_quant (only allowed when Q is max)
                    if (Q < MAXQ)
                        cpi->zbin_over_quant = 0;
                    else
                        cpi->zbin_over_quant = (zbin_oq_high + zbin_oq_low) / 2;
                }
                else
                {
                    // Update rate_correction_factor unless cpi->active_worst_quality has changed.
                    if (!active_worst_qchanged)
                        vp8_update_rate_correction_factors(cpi, 0);

                    Q = vp8_regulate_q(cpi, cpi->this_frame_target);

                    // Special case reset for qlow for constrained quality.
                    // This should only trigger where there is very substantial
                    // undershoot on a frame and the auto cq level is above
                    // the user passsed in value.
                    if ( (cpi->oxcf.end_usage == USAGE_CONSTRAINED_QUALITY) &&
                         (Q < q_low) )
                    {
                        q_low = Q;
                    }

                    while (((Q > q_high) || (cpi->zbin_over_quant > zbin_oq_high)) && (Retries < 10))
                    {
                        vp8_update_rate_correction_factors(cpi, 0);
                        Q = vp8_regulate_q(cpi, cpi->this_frame_target);
                        Retries ++;
                    }
                }

                undershoot_seen = TRUE;
            }

            // Clamp Q to upper and lower limits:
            if (Q > q_high)
                Q = q_high;
            else if (Q < q_low)
                Q = q_low;

            // Clamp cpi->zbin_over_quant
            cpi->zbin_over_quant = (cpi->zbin_over_quant < zbin_oq_low) ? zbin_oq_low : (cpi->zbin_over_quant > zbin_oq_high) ? zbin_oq_high : cpi->zbin_over_quant;

            //Loop = ((Q != last_q) || (last_zbin_oq != cpi->zbin_over_quant)) ? TRUE : FALSE;
            Loop = ((Q != last_q)) ? TRUE : FALSE;
            last_zbin_oq = cpi->zbin_over_quant;
        }
        else
            Loop = FALSE;

        if (cpi->is_src_frame_alt_ref)
            Loop = FALSE;

        if (Loop == TRUE)
        {
            vp8_restore_coding_context(cpi);
            loop_count++;
#if CONFIG_INTERNAL_STATS
            cpi->tot_recode_hits++;
#endif
        }
    }
    while (Loop == TRUE);

#if 0
    // Experimental code for lagged and one pass
    // Update stats used for one pass GF selection
    {
        /*
            int frames_so_far;
            double frame_intra_error;
            double frame_coded_error;
            double frame_pcnt_inter;
            double frame_pcnt_motion;
            double frame_mvr;
            double frame_mvr_abs;
            double frame_mvc;
            double frame_mvc_abs;
        */

        cpi->one_pass_frame_stats[cpi->one_pass_frame_index].frame_coded_error = (double)cpi->prediction_error;
        cpi->one_pass_frame_stats[cpi->one_pass_frame_index].frame_intra_error = (double)cpi->intra_error;
        cpi->one_pass_frame_stats[cpi->one_pass_frame_index].frame_pcnt_inter = (double)(100 - cpi->this_frame_percent_intra) / 100.0;
    }
#endif

    // Special case code to reduce pulsing when key frames are forced at a
    // fixed interval. Note the reconstruction error if it is the frame before
    // the force key frame
    if ( cpi->next_key_frame_forced && (cpi->twopass.frames_to_key == 0) )
    {
        cpi->ambient_err = vp8_calc_ss_err(cpi->Source,
                                           &cm->yv12_fb[cm->new_fb_idx],
                                           IF_RTCD(&cpi->rtcd.variance));
    }

    // This frame's MVs are saved and will be used in next frame's MV prediction.
    // Last frame has one more line(add to bottom) and one more column(add to right) than cm->mip. The edge elements are initialized to 0.
    if(cm->show_frame)   //do not save for altref frame
    {
        int mb_row;
        int mb_col;
        MODE_INFO *tmp = cm->mip; //point to beginning of allocated MODE_INFO arrays.

        if(cm->frame_type != KEY_FRAME)
        {
            for (mb_row = 0; mb_row < cm->mb_rows+1; mb_row ++)
            {
                for (mb_col = 0; mb_col < cm->mb_cols+1; mb_col ++)
                {
                    if(tmp->mbmi.ref_frame != INTRA_FRAME)
                        cpi->lfmv[mb_col + mb_row*(cm->mode_info_stride+1)].as_int = tmp->mbmi.mv.as_int;

                    cpi->lf_ref_frame_sign_bias[mb_col + mb_row*(cm->mode_info_stride+1)] = cm->ref_frame_sign_bias[tmp->mbmi.ref_frame];
                    cpi->lf_ref_frame[mb_col + mb_row*(cm->mode_info_stride+1)] = tmp->mbmi.ref_frame;
                    tmp++;
                }
            }
        }
    }

    // Update the GF useage maps.
    // This is done after completing the compression of a frame when all modes etc. are finalized but before loop filter
    // This is done after completing the compression of a frame when all modes etc. are finalized but before loop filter
    vp8_update_gf_useage_maps(cpi, cm, &cpi->mb);

    if (cm->frame_type == KEY_FRAME)
        cm->refresh_last_frame = 1;

#if 0
    {
        FILE *f = fopen("gfactive.stt", "a");
        fprintf(f, "%8d %8d %8d %8d %8d\n", cm->current_video_frame, (100 * cpi->gf_active_count) / (cpi->common.mb_rows * cpi->common.mb_cols), cpi->this_iiratio, cpi->next_iiratio, cm->refresh_golden_frame);
        fclose(f);
    }
#endif

    // For inter frames the current default behavior is that when
    // cm->refresh_golden_frame is set we copy the old GF over to the ARF buffer
    // This is purely an encoder decision at present.
    if (cm->refresh_golden_frame)
        cm->copy_buffer_to_arf  = 2;

    cm->frame_to_show = &cm->yv12_fb[cm->new_fb_idx];

#if WRITE_RECON_BUFFER
    if(cm->show_frame)
        write_cx_frame_to_file(cm->frame_to_show,
            cm->current_video_frame);
    else
        write_cx_frame_to_file(cm->frame_to_show,
            cm->current_video_frame+1000);
#endif

    {
        loopfilter_frame(cpi, cm);
    }

    update_reference_frames(cm);

    // Work out the segment probabilites if segmentation is enabled and
    // the map is due to be updated
    if (xd->segmentation_enabled && xd->update_mb_segmentation_map)
    {
        // Select the coding strategy for the segment map (temporal or spatial)
        choose_segmap_coding_method( cpi );

        // Take a copy of the segment map if it changed for future comparison
        vpx_memcpy( cm->last_frame_seg_map,
                    cpi->segmentation_map, cm->MBs );
    }

    // Update the common prediction model probabilities to reflect
    // the what was seen in the current frame.
    update_refpred_stats( cpi );

    // build the bitstream
    vp8_pack_bitstream(cpi, dest, size);

    /* Move storing frame_type out of the above loop since it is also
     * needed in motion search besides loopfilter */
    cm->last_frame_type = cm->frame_type;

    // Update rate control heuristics
    cpi->total_byte_count += (*size);
    cpi->projected_frame_size = (*size) << 3;

    if (!active_worst_qchanged)
        vp8_update_rate_correction_factors(cpi, 2);

    cpi->last_q[cm->frame_type] = cm->base_qindex;

    // Keep record of last boosted (KF/KF/ARF) Q value.
    // If the current frame is coded at a lower Q then we also update it.
    // If all mbs in this group are skipped only update if the Q value is
    // better than that already stored.
    // This is used to help set quality in forced key frames to reduce popping
    if ( (cm->base_qindex < cpi->last_boosted_qindex) ||
         ( (cpi->static_mb_pct < 100) &&
           ( (cm->frame_type == KEY_FRAME) ||
             cm->refresh_alt_ref_frame ||
             (cm->refresh_golden_frame && !cpi->is_src_frame_alt_ref) ) ) )
    {
        cpi->last_boosted_qindex = cm->base_qindex;
    }

    if (cm->frame_type == KEY_FRAME)
    {
        vp8_adjust_key_frame_context(cpi);
    }

    // Keep a record of ambient average Q.
    if (cm->frame_type != KEY_FRAME)
        cpi->avg_frame_qindex = (2 + 3 * cpi->avg_frame_qindex + cm->base_qindex) >> 2;

    // Keep a record from which we can calculate the average Q excluding GF updates and key frames
    if ((cm->frame_type != KEY_FRAME) && !cm->refresh_golden_frame && !cm->refresh_alt_ref_frame)
    {
        cpi->ni_frames++;
        cpi->tot_q += vp8_convert_qindex_to_q(Q);
        cpi->avg_q = cpi->tot_q / (double)cpi->ni_frames;

        // Calculate the average Q for normal inter frames (not key or GFU
        // frames).
        cpi->ni_tot_qi += Q;
        cpi->ni_av_qi = (cpi->ni_tot_qi / cpi->ni_frames);
    }

    // Update the buffer level variable.
    // Non-viewable frames are a special case and are treated as pure overhead.
    if ( !cm->show_frame )
        cpi->bits_off_target -= cpi->projected_frame_size;
    else
        cpi->bits_off_target += cpi->av_per_frame_bandwidth - cpi->projected_frame_size;

    // Clip the buffer level at the maximum buffer size
    if (cpi->bits_off_target > cpi->oxcf.maximum_buffer_size)
        cpi->bits_off_target = cpi->oxcf.maximum_buffer_size;

    // Rolling monitors of whether we are over or underspending used to help regulate min and Max Q in two pass.
    cpi->rolling_target_bits = ((cpi->rolling_target_bits * 3) + cpi->this_frame_target + 2) / 4;
    cpi->rolling_actual_bits = ((cpi->rolling_actual_bits * 3) + cpi->projected_frame_size + 2) / 4;
    cpi->long_rolling_target_bits = ((cpi->long_rolling_target_bits * 31) + cpi->this_frame_target + 16) / 32;
    cpi->long_rolling_actual_bits = ((cpi->long_rolling_actual_bits * 31) + cpi->projected_frame_size + 16) / 32;

    // Actual bits spent
    cpi->total_actual_bits    += cpi->projected_frame_size;

    // Debug stats
    cpi->total_target_vs_actual += (cpi->this_frame_target - cpi->projected_frame_size);

    cpi->buffer_level = cpi->bits_off_target;

    // Update bits left to the kf and gf groups to account for overshoot or undershoot on these frames
    if (cm->frame_type == KEY_FRAME)
    {
        cpi->twopass.kf_group_bits += cpi->this_frame_target - cpi->projected_frame_size;

        if (cpi->twopass.kf_group_bits < 0)
            cpi->twopass.kf_group_bits = 0 ;
    }
    else if (cm->refresh_golden_frame || cm->refresh_alt_ref_frame)
    {
        cpi->twopass.gf_group_bits += cpi->this_frame_target - cpi->projected_frame_size;

        if (cpi->twopass.gf_group_bits < 0)
            cpi->twopass.gf_group_bits = 0 ;
    }

    if (cm->frame_type != KEY_FRAME)
    {
        if (cpi->common.refresh_alt_ref_frame)
        {
            cpi->last_skip_false_probs[2] = cpi->prob_skip_false;
            cpi->last_skip_probs_q[2] = cm->base_qindex;
        }
        else if (cpi->common.refresh_golden_frame)
        {
            cpi->last_skip_false_probs[1] = cpi->prob_skip_false;
            cpi->last_skip_probs_q[1] = cm->base_qindex;
        }
        else
        {
            cpi->last_skip_false_probs[0] = cpi->prob_skip_false;
            cpi->last_skip_probs_q[0] = cm->base_qindex;

            //update the baseline
            cpi->base_skip_false_prob[cm->base_qindex] = cpi->prob_skip_false;

        }
    }

#if 0 && CONFIG_INTERNAL_STATS
    {
        FILE *f = fopen("tmp.stt", "a");

        vp8_clear_system_state();  //__asm emms;

        if (cpi->twopass.total_left_stats->coded_error != 0.0)
            fprintf(f, "%10d %10d %10d %10d %10d %10d %10d"
                       "%7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f"
                       "%6d %5d %5d %5d %8d %8.2f %10d %10.3f"
                       "%10.3f %8d\n",
                       cpi->common.current_video_frame, cpi->this_frame_target,
                       cpi->projected_frame_size,
                       (cpi->projected_frame_size - cpi->this_frame_target),
                       (int)cpi->total_target_vs_actual,
                       (cpi->oxcf.starting_buffer_level-cpi->bits_off_target),
                       (int)cpi->total_actual_bits,
                       vp8_convert_qindex_to_q(cm->base_qindex),
                        (double)vp8_dc_quant(cm->base_qindex,0)/4.0,
                       vp8_convert_qindex_to_q(cpi->active_best_quality),
                       vp8_convert_qindex_to_q(cpi->active_worst_quality),
                       cpi->avg_q,
                       vp8_convert_qindex_to_q(cpi->ni_av_qi),
                       vp8_convert_qindex_to_q(cpi->cq_target_quality),
                       cpi->zbin_over_quant,
                       //cpi->avg_frame_qindex, cpi->zbin_over_quant,
                       cm->refresh_golden_frame, cm->refresh_alt_ref_frame,
                       cm->frame_type, cpi->gfu_boost,
                       cpi->twopass.est_max_qcorrection_factor,
                       (int)cpi->twopass.bits_left,
                       cpi->twopass.total_left_stats->coded_error,
                       (double)cpi->twopass.bits_left /
                           cpi->twopass.total_left_stats->coded_error,
                       cpi->tot_recode_hits);
        else
            fprintf(f, "%10d %10d %10d %10d %10d %10d %10d"
                       "%7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f"
                       "%6d %5d %5d %5d %8d %8.2f %10d %10.3f"
                       "%8d\n",
                       cpi->common.current_video_frame,
                       cpi->this_frame_target, cpi->projected_frame_size,
                       (cpi->projected_frame_size - cpi->this_frame_target),
                       (int)cpi->total_target_vs_actual,
                       (cpi->oxcf.starting_buffer_level-cpi->bits_off_target),
                       (int)cpi->total_actual_bits,
                       vp8_convert_qindex_to_q(cm->base_qindex),
                        (double)vp8_dc_quant(cm->base_qindex,0)/4.0,
                       vp8_convert_qindex_to_q(cpi->active_best_quality),
                       vp8_convert_qindex_to_q(cpi->active_worst_quality),
                       cpi->avg_q,
                       vp8_convert_qindex_to_q(cpi->ni_av_qi),
                       vp8_convert_qindex_to_q(cpi->cq_target_quality),
                       cpi->zbin_over_quant,
                       //cpi->avg_frame_qindex, cpi->zbin_over_quant,
                       cm->refresh_golden_frame, cm->refresh_alt_ref_frame,
                       cm->frame_type, cpi->gfu_boost,
                       cpi->twopass.est_max_qcorrection_factor,
                       (int)cpi->twopass.bits_left,
                       cpi->twopass.total_left_stats->coded_error,
                       cpi->tot_recode_hits);

        fclose(f);

        if ( 0 )
        {
            FILE *fmodes = fopen("Modes.stt", "a");
            int i;

            fprintf(fmodes, "%6d:%1d:%1d:%1d ",
                        cpi->common.current_video_frame,
                        cm->frame_type, cm->refresh_golden_frame,
                        cm->refresh_alt_ref_frame);

            for (i = 0; i < MAX_MODES; i++)
                fprintf(fmodes, "%5d ", cpi->mode_chosen_counts[i]);

            fprintf(fmodes, "\n");

            fclose(fmodes);
        }
    }

#endif

#if 0
    // Debug stats for segment feature experiments.
    print_seg_map(cpi);
#endif

    // If this was a kf or Gf note the Q
    if ((cm->frame_type == KEY_FRAME) || cm->refresh_golden_frame || cm->refresh_alt_ref_frame)
        cm->last_kf_gf_q = cm->base_qindex;

    if (cm->refresh_golden_frame == 1)
        cm->frame_flags = cm->frame_flags | FRAMEFLAGS_GOLDEN;
    else
        cm->frame_flags = cm->frame_flags&~FRAMEFLAGS_GOLDEN;

    if (cm->refresh_alt_ref_frame == 1)
        cm->frame_flags = cm->frame_flags | FRAMEFLAGS_ALTREF;
    else
        cm->frame_flags = cm->frame_flags&~FRAMEFLAGS_ALTREF;


    if (cm->refresh_last_frame & cm->refresh_golden_frame) // both refreshed
        cpi->gold_is_last = 1;
    else if (cm->refresh_last_frame ^ cm->refresh_golden_frame) // 1 refreshed but not the other
        cpi->gold_is_last = 0;

    if (cm->refresh_last_frame & cm->refresh_alt_ref_frame) // both refreshed
        cpi->alt_is_last = 1;
    else if (cm->refresh_last_frame ^ cm->refresh_alt_ref_frame) // 1 refreshed but not the other
        cpi->alt_is_last = 0;

    if (cm->refresh_alt_ref_frame & cm->refresh_golden_frame) // both refreshed
        cpi->gold_is_alt = 1;
    else if (cm->refresh_alt_ref_frame ^ cm->refresh_golden_frame) // 1 refreshed but not the other
        cpi->gold_is_alt = 0;

    cpi->ref_frame_flags = VP8_ALT_FLAG | VP8_GOLD_FLAG | VP8_LAST_FLAG;

    if (cpi->gold_is_last)
        cpi->ref_frame_flags &= ~VP8_GOLD_FLAG;

    if (cpi->alt_is_last)
        cpi->ref_frame_flags &= ~VP8_ALT_FLAG;

    if (cpi->gold_is_alt)
        cpi->ref_frame_flags &= ~VP8_ALT_FLAG;

    if (cpi->oxcf.play_alternate && cm->refresh_alt_ref_frame && (cm->frame_type != KEY_FRAME))
        // Update the alternate reference frame stats as appropriate.
        update_alt_ref_frame_stats(cpi);
    else
        // Update the Golden frame stats as appropriate.
        update_golden_frame_stats(cpi);

    if (cm->frame_type == KEY_FRAME)
    {
        // Tell the caller that the frame was coded as a key frame
        *frame_flags = cm->frame_flags | FRAMEFLAGS_KEY;

        // As this frame is a key frame  the next defaults to an inter frame.
        cm->frame_type = INTER_FRAME;

        cpi->last_frame_percent_intra = 100;
    }
    else
    {
        *frame_flags = cm->frame_flags&~FRAMEFLAGS_KEY;

        cpi->last_frame_percent_intra = cpi->this_frame_percent_intra;
    }

    // Clear the one shot update flags for segmentation map and mode/ref loop filter deltas.
    xd->update_mb_segmentation_map = 0;
    xd->update_mb_segmentation_data = 0;
    xd->mode_ref_lf_delta_update = 0;


    // Dont increment frame counters if this was an altref buffer update not a real frame
    if (cm->show_frame)
    {
        cm->current_video_frame++;
        cpi->frames_since_key++;
    }

    // reset to normal state now that we are done.



#if 0
    {
        char filename[512];
        FILE *recon_file;
        sprintf(filename, "enc%04d.yuv", (int) cm->current_video_frame);
        recon_file = fopen(filename, "wb");
        fwrite(cm->yv12_fb[cm->lst_fb_idx].buffer_alloc,
               cm->yv12_fb[cm->lst_fb_idx].frame_size, 1, recon_file);
        fclose(recon_file);
    }
#endif
#if OUTPUT_YUV_REC
    vp8_write_yuv_rec_frame(cm);
#endif

    if(cm->show_frame)
    {
        vpx_memcpy(cm->prev_mip, cm->mip,
            (cm->mb_cols + 1) * (cm->mb_rows + 1)* sizeof(MODE_INFO));
    }
    else
    {
        vpx_memset(cm->prev_mip, 0,
            (cm->mb_cols + 1) * (cm->mb_rows + 1)* sizeof(MODE_INFO));
    }
}


static void check_gf_quality(VP8_COMP *cpi)
{
    VP8_COMMON *cm = &cpi->common;
    int gf_active_pct = (100 * cpi->gf_active_count) / (cm->mb_rows * cm->mb_cols);
    int gf_ref_usage_pct = (cpi->count_mb_ref_frame_usage[GOLDEN_FRAME] * 100) / (cm->mb_rows * cm->mb_cols);
    int last_ref_zz_useage = (cpi->inter_zz_count * 100) / (cm->mb_rows * cm->mb_cols);

    // Gf refresh is not currently being signalled
    if (cpi->gf_update_recommended == 0)
    {
        if (cpi->common.frames_since_golden > 7)
        {
            // Low use of gf
            if ((gf_active_pct < 10) || ((gf_active_pct + gf_ref_usage_pct) < 15))
            {
                // ...but last frame zero zero usage is reasonbable so a new gf might be appropriate
                if (last_ref_zz_useage >= 25)
                {
                    cpi->gf_bad_count ++;

                    if (cpi->gf_bad_count >= 8)   // Check that the condition is stable
                    {
                        cpi->gf_update_recommended = 1;
                        cpi->gf_bad_count = 0;
                    }
                }
                else
                    cpi->gf_bad_count = 0;        // Restart count as the background is not stable enough
            }
            else
                cpi->gf_bad_count = 0;            // Gf useage has picked up so reset count
        }
    }
    // If the signal is set but has not been read should we cancel it.
    else if (last_ref_zz_useage < 15)
    {
        cpi->gf_update_recommended = 0;
        cpi->gf_bad_count = 0;
    }

#if 0
    {
        FILE *f = fopen("gfneeded.stt", "a");
        fprintf(f, "%10d %10d %10d %10d %10ld \n",
                cm->current_video_frame,
                cpi->common.frames_since_golden,
                gf_active_pct, gf_ref_usage_pct,
                cpi->gf_update_recommended);
        fclose(f);
    }

#endif
}

static void Pass2Encode(VP8_COMP *cpi, unsigned long *size, unsigned char *dest, unsigned int *frame_flags)
{

    if (!cpi->common.refresh_alt_ref_frame)
        vp8_second_pass(cpi);

    encode_frame_to_data_rate(cpi, size, dest, frame_flags);
    cpi->twopass.bits_left -= 8 * *size;

    if (!cpi->common.refresh_alt_ref_frame)
    {
        double two_pass_min_rate = (double)(cpi->oxcf.target_bandwidth
            *cpi->oxcf.two_pass_vbrmin_section / 100);
        cpi->twopass.bits_left += (int64_t)(two_pass_min_rate / cpi->oxcf.frame_rate);
    }
}

//For ARM NEON, d8-d15 are callee-saved registers, and need to be saved by us.
#if HAVE_ARMV7
extern void vp8_push_neon(int64_t *store);
extern void vp8_pop_neon(int64_t *store);
#endif


int vp8_receive_raw_frame(VP8_PTR ptr, unsigned int frame_flags, YV12_BUFFER_CONFIG *sd, int64_t time_stamp, int64_t end_time)
{
#if HAVE_ARMV7
    int64_t store_reg[8];
#endif
    VP8_COMP              *cpi = (VP8_COMP *) ptr;
    VP8_COMMON            *cm = &cpi->common;
    struct vpx_usec_timer  timer;
    int                    res = 0;

#if HAVE_ARMV7
#if CONFIG_RUNTIME_CPU_DETECT
    if (cm->rtcd.flags & HAS_NEON)
#endif
    {
        vp8_push_neon(store_reg);
    }
#endif

    vpx_usec_timer_start(&timer);
    if(vp8_lookahead_push(cpi->lookahead, sd, time_stamp, end_time,
                          frame_flags, cpi->active_map_enabled ? cpi->active_map : NULL))
        res = -1;
    cm->clr_type = sd->clrtype;
    vpx_usec_timer_mark(&timer);
    cpi->time_receive_data += vpx_usec_timer_elapsed(&timer);

#if HAVE_ARMV7
#if CONFIG_RUNTIME_CPU_DETECT
    if (cm->rtcd.flags & HAS_NEON)
#endif
    {
        vp8_pop_neon(store_reg);
    }
#endif

    return res;
}


static int frame_is_reference(const VP8_COMP *cpi)
{
    const VP8_COMMON *cm = &cpi->common;
    const MACROBLOCKD *xd = &cpi->mb.e_mbd;

    return cm->frame_type == KEY_FRAME || cm->refresh_last_frame
           || cm->refresh_golden_frame || cm->refresh_alt_ref_frame
           || cm->copy_buffer_to_gf || cm->copy_buffer_to_arf
           || cm->refresh_entropy_probs
           || xd->mode_ref_lf_delta_update
           || xd->update_mb_segmentation_map || xd->update_mb_segmentation_data;
}


int vp8_get_compressed_data(VP8_PTR ptr, unsigned int *frame_flags, unsigned long *size, unsigned char *dest, int64_t *time_stamp, int64_t *time_end, int flush)
{
#if HAVE_ARMV7
    int64_t store_reg[8];
#endif
    VP8_COMP *cpi = (VP8_COMP *) ptr;
    VP8_COMMON *cm = &cpi->common;
    struct vpx_usec_timer  tsctimer;
    struct vpx_usec_timer  ticktimer;
    struct vpx_usec_timer  cmptimer;
    YV12_BUFFER_CONFIG    *force_src_buffer = NULL;

    if (!cpi)
        return -1;

#if HAVE_ARMV7
#if CONFIG_RUNTIME_CPU_DETECT
    if (cm->rtcd.flags & HAS_NEON)
#endif
    {
        vp8_push_neon(store_reg);
    }
#endif

    vpx_usec_timer_start(&cmptimer);

    cpi->source = NULL;

    // Should we code an alternate reference frame
    if (cpi->oxcf.play_alternate &&
        cpi->source_alt_ref_pending)
    {
        if ((cpi->source = vp8_lookahead_peek(cpi->lookahead,
                                              cpi->frames_till_gf_update_due)))
        {
            cpi->alt_ref_source = cpi->source;
            if (cpi->oxcf.arnr_max_frames > 0)
            {
                vp8_temporal_filter_prepare_c(cpi,
                                              cpi->frames_till_gf_update_due);
                force_src_buffer = &cpi->alt_ref_buffer;
            }
            cm->frames_till_alt_ref_frame = cpi->frames_till_gf_update_due;
            cm->refresh_alt_ref_frame = 1;
            cm->refresh_golden_frame = 0;
            cm->refresh_last_frame = 0;
            cm->show_frame = 0;
            cpi->source_alt_ref_pending = FALSE;   // Clear Pending altf Ref flag.
            cpi->is_src_frame_alt_ref = 0;
        }
    }

    if (!cpi->source)
    {
        if ((cpi->source = vp8_lookahead_pop(cpi->lookahead, flush)))
        {
            cm->show_frame = 1;

            cpi->is_src_frame_alt_ref = cpi->alt_ref_source
                                        && (cpi->source == cpi->alt_ref_source);

            if(cpi->is_src_frame_alt_ref)
                cpi->alt_ref_source = NULL;
        }
    }

    if (cpi->source)
    {
        cpi->un_scaled_source =
        cpi->Source = force_src_buffer ? force_src_buffer : &cpi->source->img;
        *time_stamp = cpi->source->ts_start;
        *time_end = cpi->source->ts_end;
        *frame_flags = cpi->source->flags;
    }
    else
    {
        *size = 0;
        if (flush && cpi->pass == 1 && !cpi->twopass.first_pass_done)
        {
            vp8_end_first_pass(cpi);    /* get last stats packet */
            cpi->twopass.first_pass_done = 1;
        }

#if HAVE_ARMV7
#if CONFIG_RUNTIME_CPU_DETECT
        if (cm->rtcd.flags & HAS_NEON)
#endif
        {
            vp8_pop_neon(store_reg);
        }
#endif
        return -1;
    }

    if (cpi->source->ts_start < cpi->first_time_stamp_ever)
    {
        cpi->first_time_stamp_ever = cpi->source->ts_start;
        cpi->last_end_time_stamp_seen = cpi->source->ts_start;
    }

    // adjust frame rates based on timestamps given
    if (!cm->refresh_alt_ref_frame)
    {
        int64_t this_duration;
        int step = 0;

        if (cpi->source->ts_start == cpi->first_time_stamp_ever)
        {
            this_duration = cpi->source->ts_end - cpi->source->ts_start;
            step = 1;
        }
        else
        {
            int64_t last_duration;

            this_duration = cpi->source->ts_end - cpi->last_end_time_stamp_seen;
            last_duration = cpi->last_end_time_stamp_seen
                            - cpi->last_time_stamp_seen;
            // do a step update if the duration changes by 10%
            if (last_duration)
                step = ((this_duration - last_duration) * 10 / last_duration);
        }

        if (this_duration)
        {
            if (step)
                vp8_new_frame_rate(cpi, 10000000.0 / this_duration);
            else
            {
                double avg_duration, interval;

                /* Average this frame's rate into the last second's average
                 * frame rate. If we haven't seen 1 second yet, then average
                 * over the whole interval seen.
                 */
                interval = cpi->source->ts_end - cpi->first_time_stamp_ever;
                if(interval > 10000000.0)
                    interval = 10000000;

                avg_duration = 10000000.0 / cpi->oxcf.frame_rate;
                avg_duration *= (interval - avg_duration + this_duration);
                avg_duration /= interval;

                vp8_new_frame_rate(cpi, 10000000.0 / avg_duration);
            }
        }

        cpi->last_time_stamp_seen = cpi->source->ts_start;
        cpi->last_end_time_stamp_seen = cpi->source->ts_end;
    }

    // start with a 0 size frame
    *size = 0;

    // Clear down mmx registers
    vp8_clear_system_state();  //__asm emms;

    cm->frame_type = INTER_FRAME;
    cm->frame_flags = *frame_flags;

#if 0

    if (cm->refresh_alt_ref_frame)
    {
        //cm->refresh_golden_frame = 1;
        cm->refresh_golden_frame = 0;
        cm->refresh_last_frame = 0;
    }
    else
    {
        cm->refresh_golden_frame = 0;
        cm->refresh_last_frame = 1;
    }

#endif
    /* find a free buffer for the new frame */
    {
        int i = 0;
        for(; i < NUM_YV12_BUFFERS; i++)
        {
            if(!cm->yv12_fb[i].flags)
            {
                cm->new_fb_idx = i;
                break;
            }
        }

        assert(i < NUM_YV12_BUFFERS );
    }
    if (cpi->pass == 1)
    {
        Pass1Encode(cpi, size, dest, frame_flags);
    }
    else if (cpi->pass == 2)
    {
        Pass2Encode(cpi, size, dest, frame_flags);
    }
    else
        encode_frame_to_data_rate(cpi, size, dest, frame_flags);

    if(cm->refresh_entropy_probs)
    {
        if(cm->refresh_alt_ref_frame)
            vpx_memcpy(&cm->lfc_a, &cm->fc, sizeof(cm->fc));
        else
            vpx_memcpy(&cm->lfc, &cm->fc, sizeof(cm->fc));
    }

    // if its a dropped frame honor the requests on subsequent frames
    if (*size > 0)
    {
        cpi->droppable = !frame_is_reference(cpi);

        // return to normal state
        cm->refresh_entropy_probs = 1;
        cm->refresh_alt_ref_frame = 0;
        cm->refresh_golden_frame = 0;
        cm->refresh_last_frame = 1;
        cm->frame_type = INTER_FRAME;

    }

    vpx_usec_timer_mark(&cmptimer);
    cpi->time_compress_data += vpx_usec_timer_elapsed(&cmptimer);

    if (cpi->b_calculate_psnr && cpi->pass != 1 && cm->show_frame)
    {
        generate_psnr_packet(cpi);
    }

#if CONFIG_INTERNAL_STATS

    if (cpi->pass != 1)
    {
        cpi->bytes += *size;

        if (cm->show_frame)
        {

            cpi->count ++;

            if (cpi->b_calculate_psnr)
            {
                double ye,ue,ve;
                double frame_psnr;
                YV12_BUFFER_CONFIG      *orig = cpi->Source;
                YV12_BUFFER_CONFIG      *recon = cpi->common.frame_to_show;
                YV12_BUFFER_CONFIG      *pp = &cm->post_proc_buffer;
                int y_samples = orig->y_height * orig->y_width ;
                int uv_samples = orig->uv_height * orig->uv_width ;
                int t_samples = y_samples + 2 * uv_samples;
                int64_t sq_error;

                ye = calc_plane_error(orig->y_buffer, orig->y_stride,
                  recon->y_buffer, recon->y_stride, orig->y_width, orig->y_height,
                  IF_RTCD(&cpi->rtcd.variance));

                ue = calc_plane_error(orig->u_buffer, orig->uv_stride,
                  recon->u_buffer, recon->uv_stride, orig->uv_width, orig->uv_height,
                  IF_RTCD(&cpi->rtcd.variance));

                ve = calc_plane_error(orig->v_buffer, orig->uv_stride,
                  recon->v_buffer, recon->uv_stride, orig->uv_width, orig->uv_height,
                  IF_RTCD(&cpi->rtcd.variance));

                sq_error = ye + ue + ve;

                frame_psnr = vp8_mse2psnr(t_samples, 255.0, sq_error);

                cpi->total_y += vp8_mse2psnr(y_samples, 255.0, ye);
                cpi->total_u += vp8_mse2psnr(uv_samples, 255.0, ue);
                cpi->total_v += vp8_mse2psnr(uv_samples, 255.0, ve);
                cpi->total_sq_error += sq_error;
                cpi->total  += frame_psnr;
                {
                    double frame_psnr2, frame_ssim2 = 0;
                    double weight = 0;

                    vp8_deblock(cm->frame_to_show, &cm->post_proc_buffer, cm->filter_level * 10 / 6, 1, 0, IF_RTCD(&cm->rtcd.postproc));
                    vp8_clear_system_state();

                    ye = calc_plane_error(orig->y_buffer, orig->y_stride,
                      pp->y_buffer, pp->y_stride, orig->y_width, orig->y_height,
                      IF_RTCD(&cpi->rtcd.variance));

                    ue = calc_plane_error(orig->u_buffer, orig->uv_stride,
                      pp->u_buffer, pp->uv_stride, orig->uv_width, orig->uv_height,
                      IF_RTCD(&cpi->rtcd.variance));

                    ve = calc_plane_error(orig->v_buffer, orig->uv_stride,
                      pp->v_buffer, pp->uv_stride, orig->uv_width, orig->uv_height,
                      IF_RTCD(&cpi->rtcd.variance));

                    sq_error = ye + ue + ve;

                    frame_psnr2 = vp8_mse2psnr(t_samples, 255.0, sq_error);

                    cpi->totalp_y += vp8_mse2psnr(y_samples, 255.0, ye);
                    cpi->totalp_u += vp8_mse2psnr(uv_samples, 255.0, ue);
                    cpi->totalp_v += vp8_mse2psnr(uv_samples, 255.0, ve);
                    cpi->total_sq_error2 += sq_error;
                    cpi->totalp  += frame_psnr2;

                    frame_ssim2 = vp8_calc_ssim(cpi->Source,
                      &cm->post_proc_buffer, 1, &weight,
                      IF_RTCD(&cpi->rtcd.variance));

                    cpi->summed_quality += frame_ssim2 * weight;
                    cpi->summed_weights += weight;
#if 0
                    {
                        FILE *f = fopen("q_used.stt", "a");
                        fprintf(f, "%5d : Y%f7.3:U%f7.3:V%f7.3:F%f7.3:S%7.3f\n",
                            cpi->common.current_video_frame,y2, u2, v2,
                            frame_psnr2, frame_ssim2);
                        fclose(f);
                    }
#endif
                }
            }

            if (cpi->b_calculate_ssimg)
            {
                double y, u, v, frame_all;
                frame_all =  vp8_calc_ssimg(cpi->Source, cm->frame_to_show,
                    &y, &u, &v, IF_RTCD(&cpi->rtcd.variance));
                cpi->total_ssimg_y += y;
                cpi->total_ssimg_u += u;
                cpi->total_ssimg_v += v;
                cpi->total_ssimg_all += frame_all;
            }

        }
    }

#endif

#if HAVE_ARMV7
#if CONFIG_RUNTIME_CPU_DETECT
    if (cm->rtcd.flags & HAS_NEON)
#endif
    {
        vp8_pop_neon(store_reg);
    }
#endif

    return 0;
}

int vp8_get_preview_raw_frame(VP8_PTR comp, YV12_BUFFER_CONFIG *dest, vp8_ppflags_t *flags)
{
    VP8_COMP *cpi = (VP8_COMP *) comp;

    if (cpi->common.refresh_alt_ref_frame)
        return -1;
    else
    {
        int ret;
#if CONFIG_POSTPROC
        ret = vp8_post_proc_frame(&cpi->common, dest, flags);
#else

        if (cpi->common.frame_to_show)
        {
            *dest = *cpi->common.frame_to_show;
            dest->y_width = cpi->common.Width;
            dest->y_height = cpi->common.Height;
            dest->uv_height = cpi->common.Height / 2;
            ret = 0;
        }
        else
        {
            ret = -1;
        }

#endif //!CONFIG_POSTPROC
        vp8_clear_system_state();
        return ret;
    }
}

int vp8_set_roimap(VP8_PTR comp, unsigned char *map, unsigned int rows, unsigned int cols, int delta_q[4], int delta_lf[4], unsigned int threshold[4])
{
    VP8_COMP *cpi = (VP8_COMP *) comp;
    signed char feature_data[SEG_LVL_MAX][MAX_MB_SEGMENTS];
    MACROBLOCKD *xd = &cpi->mb.e_mbd;
    int i;

    if (cpi->common.mb_rows != rows || cpi->common.mb_cols != cols)
        return -1;

    if (!map)
    {
        vp8_disable_segmentation((VP8_PTR)cpi);
        return 0;
    }

    // Set the segmentation Map
    vp8_set_segmentation_map((VP8_PTR)cpi, map);

    // Activate segmentation.
    vp8_enable_segmentation((VP8_PTR)cpi);

    // Set up the quant segment data
    feature_data[SEG_LVL_ALT_Q][0] = delta_q[0];
    feature_data[SEG_LVL_ALT_Q][1] = delta_q[1];
    feature_data[SEG_LVL_ALT_Q][2] = delta_q[2];
    feature_data[SEG_LVL_ALT_Q][3] = delta_q[3];

    // Set up the loop segment data s
    feature_data[SEG_LVL_ALT_LF][0] = delta_lf[0];
    feature_data[SEG_LVL_ALT_LF][1] = delta_lf[1];
    feature_data[SEG_LVL_ALT_LF][2] = delta_lf[2];
    feature_data[SEG_LVL_ALT_LF][3] = delta_lf[3];

    cpi->segment_encode_breakout[0] = threshold[0];
    cpi->segment_encode_breakout[1] = threshold[1];
    cpi->segment_encode_breakout[2] = threshold[2];
    cpi->segment_encode_breakout[3] = threshold[3];

    // Enable the loop and quant changes in the feature mask
    for ( i = 0; i < 4; i++ )
    {
        if (delta_q[i])
            enable_segfeature(xd, i, SEG_LVL_ALT_Q);
        else
            disable_segfeature(xd, i, SEG_LVL_ALT_Q);

        if (delta_lf[i])
            enable_segfeature(xd, i, SEG_LVL_ALT_LF);
        else
            disable_segfeature(xd, i, SEG_LVL_ALT_LF);
    }

    // Initialise the feature data structure
    // SEGMENT_DELTADATA    0, SEGMENT_ABSDATA      1
    vp8_set_segment_data((VP8_PTR)cpi, &feature_data[0][0], SEGMENT_DELTADATA);

    return 0;
}

int vp8_set_active_map(VP8_PTR comp, unsigned char *map, unsigned int rows, unsigned int cols)
{
    VP8_COMP *cpi = (VP8_COMP *) comp;

    if (rows == cpi->common.mb_rows && cols == cpi->common.mb_cols)
    {
        if (map)
        {
            vpx_memcpy(cpi->active_map, map, rows * cols);
            cpi->active_map_enabled = 1;
        }
        else
            cpi->active_map_enabled = 0;

        return 0;
    }
    else
    {
        //cpi->active_map_enabled = 0;
        return -1 ;
    }
}

int vp8_set_internal_size(VP8_PTR comp, VPX_SCALING horiz_mode, VPX_SCALING vert_mode)
{
    VP8_COMP *cpi = (VP8_COMP *) comp;

    if (horiz_mode <= ONETWO)
        cpi->common.horiz_scale = horiz_mode;
    else
        return -1;

    if (vert_mode <= ONETWO)
        cpi->common.vert_scale  = vert_mode;
    else
        return -1;

    return 0;
}



int vp8_calc_ss_err(YV12_BUFFER_CONFIG *source, YV12_BUFFER_CONFIG *dest, const vp8_variance_rtcd_vtable_t *rtcd)
{
    int i, j;
    int Total = 0;

    unsigned char *src = source->y_buffer;
    unsigned char *dst = dest->y_buffer;
    (void)rtcd;

    // Loop through the Y plane raw and reconstruction data summing (square differences)
    for (i = 0; i < source->y_height; i += 16)
    {
        for (j = 0; j < source->y_width; j += 16)
        {
            unsigned int sse;
            Total += VARIANCE_INVOKE(rtcd, mse16x16)(src + j, source->y_stride, dst + j, dest->y_stride, &sse);
        }

        src += 16 * source->y_stride;
        dst += 16 * dest->y_stride;
    }

    return Total;
}


int vp8_get_quantizer(VP8_PTR c)
{
    VP8_COMP   *cpi = (VP8_COMP *) c;
    return cpi->common.base_qindex;
}
