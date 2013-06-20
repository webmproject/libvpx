/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */


#ifndef VP9_ENCODER_VP9_ONYX_INT_H_
#define VP9_ENCODER_VP9_ONYX_INT_H_

#include <stdio.h>
#include "./vpx_config.h"
#include "vp9/common/vp9_onyx.h"
#include "vp9/encoder/vp9_treewriter.h"
#include "vp9/encoder/vp9_tokenize.h"
#include "vp9/common/vp9_onyxc_int.h"
#include "vp9/encoder/vp9_variance.h"
#include "vp9/encoder/vp9_encodemb.h"
#include "vp9/encoder/vp9_quantize.h"
#include "vp9/common/vp9_entropy.h"
#include "vp9/common/vp9_entropymode.h"
#include "vpx_ports/mem.h"
#include "vpx/internal/vpx_codec_internal.h"
#include "vp9/encoder/vp9_mcomp.h"
#include "vp9/common/vp9_findnearmv.h"
#include "vp9/encoder/vp9_lookahead.h"

// Experimental rate control switches
#if CONFIG_ONESHOTQ
#define ONE_SHOT_Q_ESTIMATE 0
#define STRICT_ONE_SHOT_Q 0
#define DISABLE_RC_LONG_TERM_MEM 0
#endif

// #define SPEEDSTATS 1
#if CONFIG_MULTIPLE_ARF
// Set MIN_GF_INTERVAL to 1 for the full decomposition.
#define MIN_GF_INTERVAL             2
#else
#define MIN_GF_INTERVAL             4
#endif
#define DEFAULT_GF_INTERVAL         7

#define KEY_FRAME_CONTEXT 5

#define MAX_MODES 36

#define MIN_THRESHMULT  32
#define MAX_THRESHMULT  512

#define GF_ZEROMV_ZBIN_BOOST 0
#define LF_ZEROMV_ZBIN_BOOST 0
#define MV_ZBIN_BOOST        0
#define SPLIT_MV_ZBIN_BOOST  0
#define INTRA_ZBIN_BOOST     0

typedef struct {
  nmv_context nmvc;
  int nmvjointcost[MV_JOINTS];
  int nmvcosts[2][MV_VALS];
  int nmvcosts_hp[2][MV_VALS];

  vp9_prob segment_pred_probs[PREDICTION_PROBS];
  vp9_prob intra_inter_prob[INTRA_INTER_CONTEXTS];
  vp9_prob comp_inter_prob[COMP_INTER_CONTEXTS];
  vp9_prob single_ref_prob[REF_CONTEXTS][2];
  vp9_prob comp_ref_prob[REF_CONTEXTS];

  unsigned char *last_frame_seg_map_copy;

  // 0 = Intra, Last, GF, ARF
  signed char last_ref_lf_deltas[MAX_REF_LF_DELTAS];
  // 0 = ZERO_MV, MV
  signed char last_mode_lf_deltas[MAX_MODE_LF_DELTAS];

  vp9_coeff_probs_model coef_probs[TX_SIZE_MAX_SB][BLOCK_TYPES];

  vp9_prob y_mode_prob[4][VP9_INTRA_MODES - 1];
  vp9_prob uv_mode_prob[VP9_INTRA_MODES][VP9_INTRA_MODES - 1];
  vp9_prob partition_prob[2][NUM_PARTITION_CONTEXTS][PARTITION_TYPES - 1];

  vp9_prob switchable_interp_prob[VP9_SWITCHABLE_FILTERS + 1]
                                 [VP9_SWITCHABLE_FILTERS - 1];

  int inter_mode_counts[INTER_MODE_CONTEXTS][VP9_INTER_MODES - 1][2];
  vp9_prob inter_mode_probs[INTER_MODE_CONTEXTS][VP9_INTER_MODES - 1];

  vp9_prob tx_probs_8x8p[TX_SIZE_CONTEXTS][TX_SIZE_MAX_SB - 3];
  vp9_prob tx_probs_16x16p[TX_SIZE_CONTEXTS][TX_SIZE_MAX_SB - 2];
  vp9_prob tx_probs_32x32p[TX_SIZE_CONTEXTS][TX_SIZE_MAX_SB - 1];
  vp9_prob mbskip_probs[MBSKIP_CONTEXTS];
} CODING_CONTEXT;

typedef struct {
  double frame;
  double intra_error;
  double coded_error;
  double sr_coded_error;
  double ssim_weighted_pred_err;
  double pcnt_inter;
  double pcnt_motion;
  double pcnt_second_ref;
  double pcnt_neutral;
  double MVr;
  double mvr_abs;
  double MVc;
  double mvc_abs;
  double MVrv;
  double MVcv;
  double mv_in_out_count;
  double new_mv_count;
  double duration;
  double count;
} FIRSTPASS_STATS;

typedef struct {
  int frames_so_far;
  double frame_intra_error;
  double frame_coded_error;
  double frame_pcnt_inter;
  double frame_pcnt_motion;
  double frame_mvr;
  double frame_mvr_abs;
  double frame_mvc;
  double frame_mvc_abs;
} ONEPASS_FRAMESTATS;

typedef struct {
  struct {
    int err;
    union {
      int_mv mv;
      MB_PREDICTION_MODE mode;
    } m;
  } ref[MAX_REF_FRAMES];
} MBGRAPH_MB_STATS;

typedef struct {
  MBGRAPH_MB_STATS *mb_stats;
} MBGRAPH_FRAME_STATS;

typedef enum {
  THR_ZEROMV,
  THR_DC,

  THR_NEARESTMV,
  THR_NEARMV,

  THR_ZEROG,
  THR_NEARESTG,

  THR_ZEROA,
  THR_NEARESTA,

  THR_NEARG,
  THR_NEARA,

  THR_V_PRED,
  THR_H_PRED,
  THR_D45_PRED,
  THR_D135_PRED,
  THR_D117_PRED,
  THR_D153_PRED,
  THR_D27_PRED,
  THR_D63_PRED,
  THR_TM,

  THR_NEWMV,
  THR_NEWG,
  THR_NEWA,

  THR_SPLITMV,
  THR_SPLITG,
  THR_SPLITA,

  THR_B_PRED,

  THR_COMP_ZEROLA,
  THR_COMP_NEARESTLA,
  THR_COMP_NEARLA,

  THR_COMP_ZEROGA,
  THR_COMP_NEARESTGA,
  THR_COMP_NEARGA,

  THR_COMP_NEWLA,
  THR_COMP_NEWGA,

  THR_COMP_SPLITLA,
  THR_COMP_SPLITGA,
} THR_MODES;

typedef enum {
  DIAMOND = 0,
  NSTEP = 1,
  HEX = 2
} SEARCH_METHODS;

typedef struct {
  int RD;
  SEARCH_METHODS search_method;
  int auto_filter;
  int recode_loop;
  int iterative_sub_pixel;
  int half_pixel_search;
  int quarter_pixel_search;
  int thresh_mult[MAX_MODES];
  int max_step_search_steps;
  int first_step;
  int optimize_coefficients;
  int search_best_filter;
  int static_segmentation;
  int comp_inter_joint_search_thresh;
  int adpative_rd_thresh;
  int use_lastframe_partitioning;
  int use_largest_txform;
  int use_8tap_always;
  int use_avoid_tested_higherror;
  int skip_lots_of_modes;
  int adjust_thresholds_by_speed;
  int partition_by_variance;
  int use_one_partition_size_always;
  BLOCK_SIZE_TYPE always_this_block_size;
  int use_partitions_greater_than;
  BLOCK_SIZE_TYPE greater_than_block_size;
  int use_partitions_less_than;
  BLOCK_SIZE_TYPE less_than_block_size;
} SPEED_FEATURES;

enum BlockSize {
  BLOCK_4X4,
  BLOCK_4X8,
  BLOCK_8X4,
  BLOCK_8X8,
  BLOCK_8X16,
  BLOCK_16X8,
  BLOCK_16X16,
  BLOCK_32X32,
  BLOCK_32X16,
  BLOCK_16X32,
  BLOCK_64X32,
  BLOCK_32X64,
  BLOCK_64X64,
  BLOCK_MAX_SB_SEGMENTS,
};

typedef struct VP9_COMP {

  DECLARE_ALIGNED(16, short, y_quant[QINDEX_RANGE][16]);
  DECLARE_ALIGNED(16, unsigned char, y_quant_shift[QINDEX_RANGE][16]);
  DECLARE_ALIGNED(16, short, y_zbin[QINDEX_RANGE][16]);
  DECLARE_ALIGNED(16, short, y_round[QINDEX_RANGE][16]);

  DECLARE_ALIGNED(16, short, uv_quant[QINDEX_RANGE][16]);
  DECLARE_ALIGNED(16, unsigned char, uv_quant_shift[QINDEX_RANGE][16]);
  DECLARE_ALIGNED(16, short, uv_zbin[QINDEX_RANGE][16]);
  DECLARE_ALIGNED(16, short, uv_round[QINDEX_RANGE][16]);

#if CONFIG_ALPHA
  DECLARE_ALIGNED(16, short, a_quant[QINDEX_RANGE][16]);
  DECLARE_ALIGNED(16, unsigned char, a_quant_shift[QINDEX_RANGE][16]);
  DECLARE_ALIGNED(16, short, a_zbin[QINDEX_RANGE][16]);
  DECLARE_ALIGNED(16, short, a_round[QINDEX_RANGE][16]);

  DECLARE_ALIGNED(16, short, zrun_zbin_boost_a[QINDEX_RANGE][16]);
#endif
  DECLARE_ALIGNED(16, short, zrun_zbin_boost_y[QINDEX_RANGE][16]);
  DECLARE_ALIGNED(16, short, zrun_zbin_boost_uv[QINDEX_RANGE][16]);

  MACROBLOCK mb;
  VP9_COMMON common;
  VP9_CONFIG oxcf;

  struct lookahead_ctx    *lookahead;
  struct lookahead_entry  *source;
#if CONFIG_MULTIPLE_ARF
  struct lookahead_entry  *alt_ref_source[NUM_REF_FRAMES];
#else
  struct lookahead_entry  *alt_ref_source;
#endif

  YV12_BUFFER_CONFIG *Source;
  YV12_BUFFER_CONFIG *un_scaled_source;
  YV12_BUFFER_CONFIG scaled_source;

  int source_alt_ref_pending; // frame in src_buffers has been identified to be encoded as an alt ref
  int source_alt_ref_active;  // an alt ref frame has been encoded and is usable

  int is_src_frame_alt_ref;   // source of frame to encode is an exact copy of an alt ref frame

  int gold_is_last; // golden frame same as last frame ( short circuit gold searches)
  int alt_is_last;  // Alt reference frame same as last ( short circuit altref search)
  int gold_is_alt;  // don't do both alt and gold search ( just do gold).

  int scaled_ref_idx[3];
  int lst_fb_idx;
  int gld_fb_idx;
  int alt_fb_idx;
#if CONFIG_MULTIPLE_ARF
  int alt_ref_fb_idx[NUM_REF_FRAMES - 3];
#endif
  int refresh_last_frame;
  int refresh_golden_frame;
  int refresh_alt_ref_frame;
  YV12_BUFFER_CONFIG last_frame_uf;

  TOKENEXTRA *tok;
  unsigned int tok_count[4][1 << 6];


  unsigned int frames_since_key;
  unsigned int key_frame_frequency;
  unsigned int this_key_frame_forced;
  unsigned int next_key_frame_forced;
#if CONFIG_MULTIPLE_ARF
  // Position within a frame coding order (including any additional ARF frames).
  unsigned int sequence_number;
  // Next frame in naturally occurring order that has not yet been coded.
  int next_frame_in_order;
#endif

  // Ambient reconstruction err target for force key frames
  int ambient_err;

  unsigned int mode_check_freq[MAX_MODES];
  unsigned int mode_test_hit_counts[MAX_MODES];
  unsigned int mode_chosen_counts[MAX_MODES];

  int rd_thresh_mult[MAX_MODES];
  int rd_baseline_thresh[BLOCK_SIZE_TYPES][MAX_MODES];
  int rd_threshes[BLOCK_SIZE_TYPES][MAX_MODES];
  int rd_thresh_freq_fact[BLOCK_SIZE_TYPES][MAX_MODES];

  int64_t rd_comp_pred_diff[NB_PREDICTION_TYPES];
  int rd_prediction_type_threshes[4][NB_PREDICTION_TYPES];
  unsigned int intra_inter_count[INTRA_INTER_CONTEXTS][2];
  unsigned int comp_inter_count[COMP_INTER_CONTEXTS][2];
  unsigned int single_ref_count[REF_CONTEXTS][2][2];
  unsigned int comp_ref_count[REF_CONTEXTS][2];

  // FIXME contextualize

  int64_t rd_tx_select_diff[NB_TXFM_MODES];
  int rd_tx_select_threshes[4][NB_TXFM_MODES];

  int RDMULT;
  int RDDIV;

  CODING_CONTEXT coding_context;

  // Rate targetting variables
  int this_frame_target;
  int projected_frame_size;
  int last_q[2];                   // Separate values for Intra/Inter
  int last_boosted_qindex;         // Last boosted GF/KF/ARF q

  double rate_correction_factor;
  double key_frame_rate_correction_factor;
  double gf_rate_correction_factor;

  int frames_till_gf_update_due;      // Count down till next GF

  int gf_overspend_bits;            // Total bits overspent becasue of GF boost (cumulative)

  int non_gf_bitrate_adjustment;     // Used in the few frames following a GF to recover the extra bits spent in that GF

  int kf_overspend_bits;            // Extra bits spent on key frames that need to be recovered on inter frames
  int kf_bitrate_adjustment;        // Current number of bit s to try and recover on each inter frame.
  int max_gf_interval;
  int baseline_gf_interval;
  int active_arnr_frames;           // <= cpi->oxcf.arnr_max_frames
  int active_arnr_strength;         // <= cpi->oxcf.arnr_max_strength

  int64_t key_frame_count;
  int prior_key_frame_distance[KEY_FRAME_CONTEXT];
  int per_frame_bandwidth;          // Current section per frame bandwidth target
  int av_per_frame_bandwidth;        // Average frame size target for clip
  int min_frame_bandwidth;          // Minimum allocation that should be used for any frame
  int inter_frame_target;
  double output_frame_rate;
  int64_t last_time_stamp_seen;
  int64_t last_end_time_stamp_seen;
  int64_t first_time_stamp_ever;

  int ni_av_qi;
  int ni_tot_qi;
  int ni_frames;
  int avg_frame_qindex;
  double tot_q;
  double avg_q;

  int zbin_mode_boost;
  int zbin_mode_boost_enabled;

  int64_t total_byte_count;

  int buffered_mode;

  int buffer_level;
  int bits_off_target;

  int rolling_target_bits;
  int rolling_actual_bits;

  int long_rolling_target_bits;
  int long_rolling_actual_bits;

  int64_t total_actual_bits;
  int total_target_vs_actual;        // debug stats

  int worst_quality;
  int active_worst_quality;
  int best_quality;
  int active_best_quality;

  int cq_target_quality;

  int y_mode_count[4][VP9_INTRA_MODES];
  int y_uv_mode_count[VP9_INTRA_MODES][VP9_INTRA_MODES];
  unsigned int partition_count[NUM_PARTITION_CONTEXTS][PARTITION_TYPES];

  nmv_context_counts NMVcount;

  vp9_coeff_count coef_counts[TX_SIZE_MAX_SB][BLOCK_TYPES];
  vp9_coeff_probs_model frame_coef_probs[TX_SIZE_MAX_SB][BLOCK_TYPES];
  vp9_coeff_stats frame_branch_ct[TX_SIZE_MAX_SB][BLOCK_TYPES];

  int gfu_boost;
  int last_boost;
  int kf_boost;
  int kf_zeromotion_pct;

  int64_t target_bandwidth;
  struct vpx_codec_pkt_list  *output_pkt_list;

#if 0
  // Experimental code for lagged and one pass
  ONEPASS_FRAMESTATS one_pass_frame_stats[MAX_LAG_BUFFERS];
  int one_pass_frame_index;
#endif
  MBGRAPH_FRAME_STATS mbgraph_stats[MAX_LAG_BUFFERS];
  int mbgraph_n_frames;             // number of frames filled in the above
  int static_mb_pct;                // % forced skip mbs by segmentation
  int seg0_progress, seg0_idx, seg0_cnt;

  int decimation_factor;
  int decimation_count;

  // for real time encoding
  int avg_encode_time;              // microsecond
  int avg_pick_mode_time;            // microsecond
  int speed;
  unsigned int cpu_freq;           // Mhz
  int compressor_speed;

  int interquantizer;
  int goldfreq;
  int auto_worst_q;
  int cpu_used;
  int pass;

  vp9_prob last_skip_false_probs[3][MBSKIP_CONTEXTS];
  int last_skip_probs_q[3];

  int ref_frame_flags;

  SPEED_FEATURES sf;
  int error_bins[1024];

  // Data used for real time conferencing mode to help determine if it would be good to update the gf
  int inter_zz_count;
  int gf_bad_count;
  int gf_update_recommended;

  unsigned char *segmentation_map;

  // segment threashold for encode breakout
  int  segment_encode_breakout[MAX_MB_SEGMENTS];

  unsigned char *active_map;
  unsigned int active_map_enabled;

  fractional_mv_step_fp *find_fractional_mv_step;
  vp9_full_search_fn_t full_search_sad;
  vp9_refining_search_fn_t refining_search_sad;
  vp9_diamond_search_fn_t diamond_search_sad;
  vp9_variance_fn_ptr_t fn_ptr[BLOCK_MAX_SB_SEGMENTS];
  uint64_t time_receive_data;
  uint64_t time_compress_data;
  uint64_t time_pick_lpf;
  uint64_t time_encode_mb_row;

  struct twopass_rc {
    unsigned int section_intra_rating;
    unsigned int next_iiratio;
    unsigned int this_iiratio;
    FIRSTPASS_STATS total_stats;
    FIRSTPASS_STATS this_frame_stats;
    FIRSTPASS_STATS *stats_in, *stats_in_end, *stats_in_start;
    FIRSTPASS_STATS total_left_stats;
    int first_pass_done;
    int64_t bits_left;
    int64_t clip_bits_total;
    double avg_iiratio;
    double modified_error_total;
    double modified_error_used;
    double modified_error_left;
    double kf_intra_err_min;
    double gf_intra_err_min;
    int frames_to_key;
    int maxq_max_limit;
    int maxq_min_limit;
    int static_scene_max_gf_interval;
    int kf_bits;
    // Remaining error from uncoded frames in a gf group. Two pass use only
    int64_t gf_group_error_left;

    // Projected total bits available for a key frame group of frames
    int64_t kf_group_bits;

    // Error score of frames still to be coded in kf group
    int64_t kf_group_error_left;

    // Projected Bits available for a group of frames including 1 GF or ARF
    int64_t gf_group_bits;
    // Bits for the golden frame or ARF - 2 pass only
    int gf_bits;
    int alt_extra_bits;

    int sr_update_lag;
    double est_max_qcorrection_factor;
  } twopass;

  YV12_BUFFER_CONFIG alt_ref_buffer;
  YV12_BUFFER_CONFIG *frames[MAX_LAG_BUFFERS];
  int fixed_divide[512];

#if CONFIG_INTERNAL_STATS
  int    count;
  double total_y;
  double total_u;
  double total_v;
  double total;
  double total_sq_error;
  double totalp_y;
  double totalp_u;
  double totalp_v;
  double totalp;
  double total_sq_error2;
  int    bytes;
  double summed_quality;
  double summed_weights;
  double summedp_quality;
  double summedp_weights;
  unsigned int tot_recode_hits;


  double total_ssimg_y;
  double total_ssimg_u;
  double total_ssimg_v;
  double total_ssimg_all;

  int b_calculate_ssimg;
#endif
  int b_calculate_psnr;

  // Per MB activity measurement
  unsigned int activity_avg;
  unsigned int *mb_activity_map;
  int *mb_norm_activity_map;
  int output_partition;

  /* force next frame to intra when kf_auto says so */
  int force_next_frame_intra;

  int droppable;

  int dummy_packing;    /* flag to indicate if packing is dummy */

  unsigned int switchable_interp_count[VP9_SWITCHABLE_FILTERS + 1]
                                      [VP9_SWITCHABLE_FILTERS];
  unsigned int best_switchable_interp_count[VP9_SWITCHABLE_FILTERS];

  int initial_width;
  int initial_height;

#if CONFIG_MULTIPLE_ARF
  // ARF tracking variables.
  int multi_arf_enabled;
  unsigned int frame_coding_order_period;
  unsigned int new_frame_coding_order_period;
  int frame_coding_order[MAX_LAG_BUFFERS * 2];
  int arf_buffer_idx[MAX_LAG_BUFFERS * 3 / 2];
  int arf_weight[MAX_LAG_BUFFERS];
  int arf_buffered;
  int this_frame_weight;
  int max_arf_level;
#endif

#ifdef ENTROPY_STATS
  int64_t mv_ref_stats[INTER_MODE_CONTEXTS][VP9_INTER_MODES - 1][2];
#endif
} VP9_COMP;

static int get_ref_frame_idx(VP9_COMP *cpi, MV_REFERENCE_FRAME ref_frame) {
  if (ref_frame == LAST_FRAME) {
    return cpi->lst_fb_idx;
  } else if (ref_frame == GOLDEN_FRAME) {
    return cpi->gld_fb_idx;
  } else {
    return cpi->alt_fb_idx;
  }
}

void vp9_encode_frame(VP9_COMP *cpi);

void vp9_pack_bitstream(VP9_COMP *cpi, unsigned char *dest,
                        unsigned long *size);

void vp9_activity_masking(VP9_COMP *cpi, MACROBLOCK *x);

void vp9_set_speed_features(VP9_COMP *cpi);

extern int vp9_calc_ss_err(YV12_BUFFER_CONFIG *source,
                           YV12_BUFFER_CONFIG *dest);

extern void vp9_alloc_compressor_data(VP9_COMP *cpi);

#if CONFIG_DEBUG
#define CHECK_MEM_ERROR(lval,expr) do {\
    lval = (expr); \
    if(!lval) \
      vpx_internal_error(&cpi->common.error, VPX_CODEC_MEM_ERROR,\
                         "Failed to allocate "#lval" at %s:%d", \
                         __FILE__,__LINE__);\
  } while(0)
#else
#define CHECK_MEM_ERROR(lval,expr) do {\
    lval = (expr); \
    if(!lval) \
      vpx_internal_error(&cpi->common.error, VPX_CODEC_MEM_ERROR,\
                         "Failed to allocate "#lval);\
  } while(0)
#endif

#endif  // VP9_ENCODER_VP9_ONYX_INT_H_
