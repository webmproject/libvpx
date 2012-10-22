/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */


#include "vpx_ports/config.h"
#include "encodemb.h"
#include "encodemv.h"
#include "vp8/common/common.h"
#include "onyx_int.h"
#include "vp8/common/extend.h"
#include "vp8/common/entropymode.h"
#include "vp8/common/quant_common.h"
#include "segmentation.h"
#include "vp8/common/setupintrarecon.h"
#include "encodeintra.h"
#include "vp8/common/reconinter.h"
#include "vp8/common/invtrans.h"
#include "rdopt.h"
#include "vp8/common/findnearmv.h"
#include "vp8/common/reconintra.h"
#include "vp8/common/seg_common.h"
#include "vpx_rtcd.h"
#include <stdio.h>
#include <math.h>
#include <limits.h>
#include "vp8/common/subpixel.h"
#include "vpx_ports/vpx_timer.h"
#include "vp8/common/pred_common.h"

#define DBG_PRNT_SEGMAP 0
#if CONFIG_NEWBESTREFMV
#include "vp8/common/mvref_common.h"
#endif


#if CONFIG_RUNTIME_CPU_DETECT
#define RTCD(x)     &cpi->common.rtcd.x
#define IF_RTCD(x)  (x)
#else
#define RTCD(x)     NULL
#define IF_RTCD(x)  NULL
#endif

#ifdef ENC_DEBUG
int enc_debug = 0;
int mb_row_debug, mb_col_debug;
#endif

extern void vp8cx_initialize_me_consts(VP8_COMP *cpi, int QIndex);
extern void vp8_auto_select_speed(VP8_COMP *cpi);
extern void vp8cx_init_mbrthread_data(VP8_COMP *cpi,
                                      MACROBLOCK *x,
                                      MB_ROW_COMP *mbr_ei,
                                      int mb_row,
                                      int count);
int64_t vp8_rd_pick_inter_mode_sb(VP8_COMP *cpi, MACROBLOCK *x,
                              int recon_yoffset, int recon_uvoffset,
                              int *returnrate, int *returndistortion);
extern void vp8cx_pick_mode_inter_macroblock(VP8_COMP *cpi, MACROBLOCK *x,
                                            int recon_yoffset,
                                            int recon_uvoffset, int *r, int *d);
void vp8_build_block_offsets(MACROBLOCK *x);
void vp8_setup_block_ptrs(MACROBLOCK *x);
void vp8cx_encode_inter_macroblock(VP8_COMP *cpi, MACROBLOCK *x, TOKENEXTRA **t,
                                   int recon_yoffset, int recon_uvoffset,
                                   int output_enabled);
void vp8cx_encode_inter_superblock(VP8_COMP *cpi, MACROBLOCK *x, TOKENEXTRA **t,
                                   int recon_yoffset, int recon_uvoffset, int mb_col, int mb_row);
void vp8cx_encode_intra_macro_block(VP8_COMP *cpi, MACROBLOCK *x,
                                    TOKENEXTRA **t, int output_enabled);
void vp8cx_encode_intra_super_block(VP8_COMP *cpi,
                                    MACROBLOCK *x,
                                    TOKENEXTRA **t, int mb_col);
static void adjust_act_zbin(VP8_COMP *cpi, MACROBLOCK *x);

#ifdef MODE_STATS
unsigned int inter_y_modes[MB_MODE_COUNT];
unsigned int inter_uv_modes[VP8_UV_MODES];
unsigned int inter_b_modes[B_MODE_COUNT];
unsigned int y_modes[VP8_YMODES];
unsigned int i8x8_modes[VP8_I8X8_MODES];
unsigned int uv_modes[VP8_UV_MODES];
unsigned int uv_modes_y[VP8_YMODES][VP8_UV_MODES];
unsigned int b_modes[B_MODE_COUNT];
#endif


/* activity_avg must be positive, or flat regions could get a zero weight
 *  (infinite lambda), which confounds analysis.
 * This also avoids the need for divide by zero checks in
 *  vp8_activity_masking().
 */
#define VP8_ACTIVITY_AVG_MIN (64)

/* This is used as a reference when computing the source variance for the
 *  purposes of activity masking.
 * Eventually this should be replaced by custom no-reference routines,
 *  which will be faster.
 */
static const unsigned char VP8_VAR_OFFS[16] = {
  128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128
};


// Original activity measure from Tim T's code.
static unsigned int tt_activity_measure(VP8_COMP *cpi, MACROBLOCK *x) {
  unsigned int act;
  unsigned int sse;
  /* TODO: This could also be done over smaller areas (8x8), but that would
   *  require extensive changes elsewhere, as lambda is assumed to be fixed
   *  over an entire MB in most of the code.
   * Another option is to compute four 8x8 variances, and pick a single
   *  lambda using a non-linear combination (e.g., the smallest, or second
   *  smallest, etc.).
   */
  act = vp8_variance16x16(x->src.y_buffer, x->src.y_stride, VP8_VAR_OFFS, 0,
                          &sse);
  act = act << 4;

  /* If the region is flat, lower the activity some more. */
  if (act < 8 << 12)
    act = act < 5 << 12 ? act : 5 << 12;

  return act;
}

// Stub for alternative experimental activity measures.
static unsigned int alt_activity_measure(VP8_COMP *cpi,
                                         MACROBLOCK *x, int use_dc_pred) {
  return vp8_encode_intra(cpi, x, use_dc_pred);
}


// Measure the activity of the current macroblock
// What we measure here is TBD so abstracted to this function
#define ALT_ACT_MEASURE 1
static unsigned int mb_activity_measure(VP8_COMP *cpi, MACROBLOCK *x,
                                        int mb_row, int mb_col) {
  unsigned int mb_activity;

  if (ALT_ACT_MEASURE) {
    int use_dc_pred = (mb_col || mb_row) && (!mb_col || !mb_row);

    // Or use and alternative.
    mb_activity = alt_activity_measure(cpi, x, use_dc_pred);
  } else {
    // Original activity measure from Tim T's code.
    mb_activity = tt_activity_measure(cpi, x);
  }

  if (mb_activity < VP8_ACTIVITY_AVG_MIN)
    mb_activity = VP8_ACTIVITY_AVG_MIN;

  return mb_activity;
}

// Calculate an "average" mb activity value for the frame
#define ACT_MEDIAN 0
static void calc_av_activity(VP8_COMP *cpi, int64_t activity_sum) {
#if ACT_MEDIAN
  // Find median: Simple n^2 algorithm for experimentation
  {
    unsigned int median;
    unsigned int i, j;
    unsigned int *sortlist;
    unsigned int tmp;

    // Create a list to sort to
    CHECK_MEM_ERROR(sortlist,
    vpx_calloc(sizeof(unsigned int),
    cpi->common.MBs));

    // Copy map to sort list
    vpx_memcpy(sortlist, cpi->mb_activity_map,
    sizeof(unsigned int) * cpi->common.MBs);


    // Ripple each value down to its correct position
    for (i = 1; i < cpi->common.MBs; i ++) {
      for (j = i; j > 0; j --) {
        if (sortlist[j] < sortlist[j - 1]) {
          // Swap values
          tmp = sortlist[j - 1];
          sortlist[j - 1] = sortlist[j];
          sortlist[j] = tmp;
        } else
          break;
      }
    }

    // Even number MBs so estimate median as mean of two either side.
    median = (1 + sortlist[cpi->common.MBs >> 1] +
              sortlist[(cpi->common.MBs >> 1) + 1]) >> 1;

    cpi->activity_avg = median;

    vpx_free(sortlist);
  }
#else
  // Simple mean for now
  cpi->activity_avg = (unsigned int)(activity_sum / cpi->common.MBs);
#endif

  if (cpi->activity_avg < VP8_ACTIVITY_AVG_MIN)
    cpi->activity_avg = VP8_ACTIVITY_AVG_MIN;

  // Experimental code: return fixed value normalized for several clips
  if (ALT_ACT_MEASURE)
    cpi->activity_avg = 100000;
}

#define USE_ACT_INDEX   0
#define OUTPUT_NORM_ACT_STATS   0

#if USE_ACT_INDEX
// Calculate and activity index for each mb
static void calc_activity_index(VP8_COMP *cpi, MACROBLOCK *x) {
  VP8_COMMON *const cm = &cpi->common;
  int mb_row, mb_col;

  int64_t act;
  int64_t a;
  int64_t b;

#if OUTPUT_NORM_ACT_STATS
  FILE *f = fopen("norm_act.stt", "a");
  fprintf(f, "\n%12d\n", cpi->activity_avg);
#endif

  // Reset pointers to start of activity map
  x->mb_activity_ptr = cpi->mb_activity_map;

  // Calculate normalized mb activity number.
  for (mb_row = 0; mb_row < cm->mb_rows; mb_row++) {
    // for each macroblock col in image
    for (mb_col = 0; mb_col < cm->mb_cols; mb_col++) {
      // Read activity from the map
      act = *(x->mb_activity_ptr);

      // Calculate a normalized activity number
      a = act + 4 * cpi->activity_avg;
      b = 4 * act + cpi->activity_avg;

      if (b >= a)
        *(x->activity_ptr) = (int)((b + (a >> 1)) / a) - 1;
      else
        *(x->activity_ptr) = 1 - (int)((a + (b >> 1)) / b);

#if OUTPUT_NORM_ACT_STATS
      fprintf(f, " %6d", *(x->mb_activity_ptr));
#endif
      // Increment activity map pointers
      x->mb_activity_ptr++;
    }

#if OUTPUT_NORM_ACT_STATS
    fprintf(f, "\n");
#endif

  }

#if OUTPUT_NORM_ACT_STATS
  fclose(f);
#endif

}
#endif

// Loop through all MBs. Note activity of each, average activity and
// calculate a normalized activity for each
static void build_activity_map(VP8_COMP *cpi) {
  MACROBLOCK *const x = &cpi->mb;
  MACROBLOCKD *xd = &x->e_mbd;
  VP8_COMMON *const cm = &cpi->common;

#if ALT_ACT_MEASURE
  YV12_BUFFER_CONFIG *new_yv12 = &cm->yv12_fb[cm->new_fb_idx];
  int recon_yoffset;
  int recon_y_stride = new_yv12->y_stride;
#endif

  int mb_row, mb_col;
  unsigned int mb_activity;
  int64_t activity_sum = 0;

  // for each macroblock row in image
  for (mb_row = 0; mb_row < cm->mb_rows; mb_row++) {
#if ALT_ACT_MEASURE
    // reset above block coeffs
    xd->up_available = (mb_row != 0);
    recon_yoffset = (mb_row * recon_y_stride * 16);
#endif
    // for each macroblock col in image
    for (mb_col = 0; mb_col < cm->mb_cols; mb_col++) {
#if ALT_ACT_MEASURE
      xd->dst.y_buffer = new_yv12->y_buffer + recon_yoffset;
      xd->left_available = (mb_col != 0);
      recon_yoffset += 16;
#endif
      // Copy current mb to a buffer
      vp8_copy_mem16x16(x->src.y_buffer, x->src.y_stride, x->thismb, 16);

      // measure activity
      mb_activity = mb_activity_measure(cpi, x, mb_row, mb_col);

      // Keep frame sum
      activity_sum += mb_activity;

      // Store MB level activity details.
      *x->mb_activity_ptr = mb_activity;

      // Increment activity map pointer
      x->mb_activity_ptr++;

      // adjust to the next column of source macroblocks
      x->src.y_buffer += 16;
    }


    // adjust to the next row of mbs
    x->src.y_buffer += 16 * x->src.y_stride - 16 * cm->mb_cols;

#if ALT_ACT_MEASURE
    // extend the recon for intra prediction
    vp8_extend_mb_row(new_yv12, xd->dst.y_buffer + 16,
                      xd->dst.u_buffer + 8, xd->dst.v_buffer + 8);
#endif

  }

  // Calculate an "average" MB activity
  calc_av_activity(cpi, activity_sum);

#if USE_ACT_INDEX
  // Calculate an activity index number of each mb
  calc_activity_index(cpi, x);
#endif

}

// Macroblock activity masking
void vp8_activity_masking(VP8_COMP *cpi, MACROBLOCK *x) {
#if USE_ACT_INDEX
  x->rdmult += *(x->mb_activity_ptr) * (x->rdmult >> 2);
  x->errorperbit = x->rdmult * 100 / (110 * x->rddiv);
  x->errorperbit += (x->errorperbit == 0);
#else
  int64_t a;
  int64_t b;
  int64_t act = *(x->mb_activity_ptr);

  // Apply the masking to the RD multiplier.
  a = act + (2 * cpi->activity_avg);
  b = (2 * act) + cpi->activity_avg;

  x->rdmult = (unsigned int)(((int64_t)x->rdmult * b + (a >> 1)) / a);
  x->errorperbit = x->rdmult * 100 / (110 * x->rddiv);
  x->errorperbit += (x->errorperbit == 0);
#endif

  // Activity based Zbin adjustment
  adjust_act_zbin(cpi, x);
}

static void update_state(VP8_COMP *cpi, MACROBLOCK *x, PICK_MODE_CONTEXT *ctx) {
  int i;
  MACROBLOCKD *xd = &x->e_mbd;
  MODE_INFO *mi = &ctx->mic;
  MB_MODE_INFO * mbmi = &xd->mode_info_context->mbmi;
  int mb_mode = mi->mbmi.mode;
  int mb_mode_index = ctx->best_mode_index;

#if CONFIG_DEBUG
  assert(mb_mode < MB_MODE_COUNT);
  assert(mb_mode_index < MAX_MODES);
  assert(mi->mbmi.ref_frame < MAX_REF_FRAMES);
#endif

  // Restore the coding context of the MB to that that was in place
  // when the mode was picked for it
  vpx_memcpy(xd->mode_info_context, mi, sizeof(MODE_INFO));
#if CONFIG_SUPERBLOCKS
  if (mi->mbmi.encoded_as_sb) {
    vpx_memcpy(xd->mode_info_context + 1, mi, sizeof(MODE_INFO));
    vpx_memcpy(xd->mode_info_context + cpi->common.mode_info_stride, mi, sizeof(MODE_INFO));
    vpx_memcpy(xd->mode_info_context + cpi->common.mode_info_stride + 1, mi, sizeof(MODE_INFO));
  }
#endif

  if (mb_mode == B_PRED) {
    for (i = 0; i < 16; i++) {
      xd->block[i].bmi.as_mode = xd->mode_info_context->bmi[i].as_mode;
      assert(xd->block[i].bmi.as_mode.first < MB_MODE_COUNT);
    }
  } else if (mb_mode == I8X8_PRED) {
    for (i = 0; i < 16; i++) {
      xd->block[i].bmi = xd->mode_info_context->bmi[i];
    }
  } else if (mb_mode == SPLITMV) {
    vpx_memcpy(x->partition_info, &ctx->partition_info,
               sizeof(PARTITION_INFO));

    mbmi->mv[0].as_int = x->partition_info->bmi[15].mv.as_int;
    mbmi->mv[1].as_int = x->partition_info->bmi[15].second_mv.as_int;
  }

  {
    int segment_id = mbmi->segment_id;
    if (!segfeature_active(xd, segment_id, SEG_LVL_EOB) ||
        get_segdata(xd, segment_id, SEG_LVL_EOB)) {
      for (i = 0; i < NB_TXFM_MODES; i++) {
        cpi->rd_tx_select_diff[i] += ctx->txfm_rd_diff[i];
      }
    }
  }

  if (cpi->common.frame_type == KEY_FRAME) {
    // Restore the coding modes to that held in the coding context
    // if (mb_mode == B_PRED)
    //    for (i = 0; i < 16; i++)
    //    {
    //        xd->block[i].bmi.as_mode =
    //                          xd->mode_info_context->bmi[i].as_mode;
    //        assert(xd->mode_info_context->bmi[i].as_mode < MB_MODE_COUNT);
    //    }
#if CONFIG_INTERNAL_STATS
    static const int kf_mode_index[] = {
      THR_DC /*DC_PRED*/,
      THR_V_PRED /*V_PRED*/,
      THR_H_PRED /*H_PRED*/,
      THR_D45_PRED /*D45_PRED*/,
      THR_D135_PRED /*D135_PRED*/,
      THR_D117_PRED /*D117_PRED*/,
      THR_D153_PRED /*D153_PRED*/,
      THR_D27_PRED /*D27_PRED*/,
      THR_D63_PRED /*D63_PRED*/,
      THR_TM /*TM_PRED*/,
      THR_I8X8_PRED /*I8X8_PRED*/,
      THR_B_PRED /*B_PRED*/,
    };
    cpi->mode_chosen_counts[kf_mode_index[mb_mode]]++;
#endif
  } else {
    /*
            // Reduce the activation RD thresholds for the best choice mode
            if ((cpi->rd_baseline_thresh[mb_mode_index] > 0) &&
                (cpi->rd_baseline_thresh[mb_mode_index] < (INT_MAX >> 2)))
            {
                int best_adjustment = (cpi->rd_thresh_mult[mb_mode_index] >> 2);

                cpi->rd_thresh_mult[mb_mode_index] =
                        (cpi->rd_thresh_mult[mb_mode_index]
                         >= (MIN_THRESHMULT + best_adjustment)) ?
                                cpi->rd_thresh_mult[mb_mode_index] - best_adjustment :
                                MIN_THRESHMULT;
                cpi->rd_threshes[mb_mode_index] =
                        (cpi->rd_baseline_thresh[mb_mode_index] >> 7)
                        * cpi->rd_thresh_mult[mb_mode_index];

            }
    */
    // Note how often each mode chosen as best
    cpi->mode_chosen_counts[mb_mode_index]++;

    rd_update_mvcount(cpi, x, &ctx->best_ref_mv, &ctx->second_best_ref_mv);

    cpi->prediction_error += ctx->distortion;
    cpi->intra_error += ctx->intra_error;

    cpi->rd_comp_pred_diff[0] += ctx->single_pred_diff;
    cpi->rd_comp_pred_diff[1] += ctx->comp_pred_diff;
    cpi->rd_comp_pred_diff[2] += ctx->hybrid_pred_diff;
  }
}

static void pick_mb_modes(VP8_COMP *cpi,
                          VP8_COMMON *cm,
                          int mb_row,
                          int mb_col,
                          MACROBLOCK  *x,
                          MACROBLOCKD *xd,
                          TOKENEXTRA **tp,
                          int *totalrate,
                          int *totaldist) {
  int i;
  int map_index;
  int recon_yoffset, recon_uvoffset;
  int ref_fb_idx = cm->lst_fb_idx;
  int dst_fb_idx = cm->new_fb_idx;
  int recon_y_stride = cm->yv12_fb[ref_fb_idx].y_stride;
  int recon_uv_stride = cm->yv12_fb[ref_fb_idx].uv_stride;
  ENTROPY_CONTEXT_PLANES left_context[2];
  ENTROPY_CONTEXT_PLANES above_context[2];
  ENTROPY_CONTEXT_PLANES *initial_above_context_ptr = cm->above_context
                                                      + mb_col;

  // Offsets to move pointers from MB to MB within a SB in raster order
  int row_delta[4] = { 0, +1,  0, -1};
  int col_delta[4] = { +1, -1, +1, +1};

  /* Function should not modify L & A contexts; save and restore on exit */
  vpx_memcpy(left_context,
             cm->left_context,
             sizeof(left_context));
  vpx_memcpy(above_context,
             initial_above_context_ptr,
             sizeof(above_context));

  /* Encode MBs in raster order within the SB */
  for (i = 0; i < 4; i++) {
    int dy = row_delta[i];
    int dx = col_delta[i];
    int offset_unextended = dy * cm->mb_cols + dx;
    int offset_extended   = dy * xd->mode_info_stride + dx;
    MB_MODE_INFO * mbmi = &xd->mode_info_context->mbmi;

    // TODO Many of the index items here can be computed more efficiently!

    if ((mb_row >= cm->mb_rows) || (mb_col >= cm->mb_cols)) {
      // MB lies outside frame, move on
      mb_row += dy;
      mb_col += dx;

      // Update pointers
      x->src.y_buffer += 16 * (dx + dy * x->src.y_stride);
      x->src.u_buffer += 8  * (dx + dy * x->src.uv_stride);
      x->src.v_buffer += 8  * (dx + dy * x->src.uv_stride);

      x->gf_active_ptr += offset_unextended;
      x->partition_info += offset_extended;
      xd->mode_info_context += offset_extended;
      xd->prev_mode_info_context += offset_extended;
#if CONFIG_DEBUG
      assert((xd->prev_mode_info_context - cpi->common.prev_mip) ==
             (xd->mode_info_context - cpi->common.mip));
#endif
      continue;
    }

    // Index of the MB in the SB 0..3
    xd->mb_index = i;

    map_index = (mb_row * cpi->common.mb_cols) + mb_col;
    x->mb_activity_ptr = &cpi->mb_activity_map[map_index];

    // set above context pointer
    xd->above_context = cm->above_context + mb_col;

    // Restore the appropriate left context depending on which
    // row in the SB the MB is situated
    xd->left_context = cm->left_context + (i >> 1);

    // Set up distance of MB to edge of frame in 1/8th pel units
    xd->mb_to_top_edge    = -((mb_row * 16) << 3);
    xd->mb_to_left_edge   = -((mb_col * 16) << 3);
    xd->mb_to_bottom_edge = ((cm->mb_rows - 1 - mb_row) * 16) << 3;
    xd->mb_to_right_edge  = ((cm->mb_cols - 1 - mb_col) * 16) << 3;

    // Set up limit values for MV components to prevent them from
    // extending beyond the UMV borders assuming 16x16 block size
    x->mv_row_min = -((mb_row * 16) + VP8BORDERINPIXELS - INTERP_EXTEND);
    x->mv_col_min = -((mb_col * 16) + VP8BORDERINPIXELS - INTERP_EXTEND);
    x->mv_row_max = ((cm->mb_rows - mb_row) * 16 +
                     (VP8BORDERINPIXELS - 16 - INTERP_EXTEND));
    x->mv_col_max = ((cm->mb_cols - mb_col) * 16 +
                     (VP8BORDERINPIXELS - 16 - INTERP_EXTEND));

    xd->up_available   = (mb_row != 0);
    xd->left_available = (mb_col != 0);

    recon_yoffset  = (mb_row * recon_y_stride * 16) + (mb_col * 16);
    recon_uvoffset = (mb_row * recon_uv_stride * 8) + (mb_col *  8);

    xd->dst.y_buffer = cm->yv12_fb[dst_fb_idx].y_buffer + recon_yoffset;
    xd->dst.u_buffer = cm->yv12_fb[dst_fb_idx].u_buffer + recon_uvoffset;
    xd->dst.v_buffer = cm->yv12_fb[dst_fb_idx].v_buffer + recon_uvoffset;

    // Copy current MB to a work buffer
    vp8_copy_mem16x16(x->src.y_buffer, x->src.y_stride, x->thismb, 16);

    x->rddiv = cpi->RDDIV;
    x->rdmult = cpi->RDMULT;

    if (cpi->oxcf.tuning == VP8_TUNE_SSIM)
      vp8_activity_masking(cpi, x);

    // Is segmentation enabled
    if (xd->segmentation_enabled) {
      // Code to set segment id in xd->mbmi.segment_id
      if (xd->update_mb_segmentation_map)
        mbmi->segment_id = cpi->segmentation_map[map_index];
      else
        mbmi->segment_id = cm->last_frame_seg_map[map_index];
      if (mbmi->segment_id > 3)
        mbmi->segment_id = 0;

      vp8cx_mb_init_quantizer(cpi, x);
    } else
      // Set to Segment 0 by default
      mbmi->segment_id = 0;

    x->active_ptr = cpi->active_map + map_index;

#if CONFIG_SUPERBLOCKS
    xd->mode_info_context->mbmi.encoded_as_sb = 0;
#endif

    cpi->update_context = 0;    // TODO Do we need this now??

    // Find best coding mode & reconstruct the MB so it is available
    // as a predictor for MBs that follow in the SB
    if (cm->frame_type == KEY_FRAME) {
      int r, d;
      vp8_rd_pick_intra_mode(cpi, x, &r, &d);
      *totalrate += r;
      *totaldist += d;

      // Dummy encode, do not do the tokenization
      vp8cx_encode_intra_macro_block(cpi, x, tp, 0);
      // Note the encoder may have changed the segment_id

      // Save the coding context
      vpx_memcpy(&x->mb_context[i].mic, xd->mode_info_context,
                 sizeof(MODE_INFO));
    } else {
      int seg_id, r, d;

      if (xd->segmentation_enabled && cpi->seg0_cnt > 0 &&
          !segfeature_active(xd, 0, SEG_LVL_REF_FRAME) &&
          segfeature_active(xd, 1, SEG_LVL_REF_FRAME) &&
          check_segref(xd, 1, INTRA_FRAME)  +
          check_segref(xd, 1, LAST_FRAME)   +
          check_segref(xd, 1, GOLDEN_FRAME) +
          check_segref(xd, 1, ALTREF_FRAME) == 1) {
        cpi->seg0_progress = (cpi->seg0_idx << 16) / cpi->seg0_cnt;
      } else {
        cpi->seg0_progress = (((mb_col & ~1) * 2 + (mb_row & ~1) * cm->mb_cols + i) << 16) / cm->MBs;
      }

      vp8cx_pick_mode_inter_macroblock(cpi, x, recon_yoffset,
                                       recon_uvoffset, &r, &d);
      *totalrate += r;
      *totaldist += d;

      // Dummy encode, do not do the tokenization
      vp8cx_encode_inter_macroblock(cpi, x, tp,
                                    recon_yoffset, recon_uvoffset, 0);

      seg_id = mbmi->segment_id;
      if (cpi->mb.e_mbd.segmentation_enabled && seg_id == 0) {
        cpi->seg0_idx++;
      }
      if (!xd->segmentation_enabled ||
          !segfeature_active(xd, seg_id, SEG_LVL_REF_FRAME) ||
          check_segref(xd, seg_id, INTRA_FRAME)  +
          check_segref(xd, seg_id, LAST_FRAME)   +
          check_segref(xd, seg_id, GOLDEN_FRAME) +
          check_segref(xd, seg_id, ALTREF_FRAME) > 1) {
        // Get the prediction context and status
        int pred_flag = get_pred_flag(xd, PRED_REF);
        int pred_context = get_pred_context(cm, xd, PRED_REF);

        // Count prediction success
        cpi->ref_pred_count[pred_context][pred_flag]++;
      }
    }

    // Next MB
    mb_row += dy;
    mb_col += dx;

    x->src.y_buffer += 16 * (dx + dy * x->src.y_stride);
    x->src.u_buffer += 8  * (dx + dy * x->src.uv_stride);
    x->src.v_buffer += 8  * (dx + dy * x->src.uv_stride);

    x->gf_active_ptr += offset_unextended;
    x->partition_info += offset_extended;
    xd->mode_info_context += offset_extended;
    xd->prev_mode_info_context += offset_extended;

#if CONFIG_DEBUG
    assert((xd->prev_mode_info_context - cpi->common.prev_mip) ==
           (xd->mode_info_context - cpi->common.mip));
#endif
  }

  /* Restore L & A coding context to those in place on entry */
  vpx_memcpy(cm->left_context,
             left_context,
             sizeof(left_context));
  vpx_memcpy(initial_above_context_ptr,
             above_context,
             sizeof(above_context));
}

#if CONFIG_SUPERBLOCKS
static void pick_sb_modes (VP8_COMP *cpi,
                           VP8_COMMON *cm,
                           int mb_row,
                           int mb_col,
                           MACROBLOCK  *x,
                           MACROBLOCKD *xd,
                           TOKENEXTRA **tp,
                           int *totalrate,
                           int *totaldist)
{
  int map_index;
  int recon_yoffset, recon_uvoffset;
  int ref_fb_idx = cm->lst_fb_idx;
  int dst_fb_idx = cm->new_fb_idx;
  int recon_y_stride = cm->yv12_fb[ref_fb_idx].y_stride;
  int recon_uv_stride = cm->yv12_fb[ref_fb_idx].uv_stride;
  ENTROPY_CONTEXT_PLANES left_context[2];
  ENTROPY_CONTEXT_PLANES above_context[2];
  ENTROPY_CONTEXT_PLANES *initial_above_context_ptr = cm->above_context
    + mb_col;

  /* Function should not modify L & A contexts; save and restore on exit */
  vpx_memcpy (left_context,
              cm->left_context,
              sizeof(left_context));
  vpx_memcpy (above_context,
              initial_above_context_ptr,
              sizeof(above_context));

  map_index = (mb_row * cpi->common.mb_cols) + mb_col;
  x->mb_activity_ptr = &cpi->mb_activity_map[map_index];

  /* set above context pointer */
  xd->above_context = cm->above_context + mb_col;

  /* Restore the appropriate left context depending on which
   * row in the SB the MB is situated */
  xd->left_context = cm->left_context;

  // Set up distance of MB to edge of frame in 1/8th pel units
  xd->mb_to_top_edge    = -((mb_row * 16) << 3);
  xd->mb_to_left_edge   = -((mb_col * 16) << 3);
  xd->mb_to_bottom_edge = ((cm->mb_rows - 1 - mb_row) * 16) << 3;
  xd->mb_to_right_edge  = ((cm->mb_cols - 1 - mb_col) * 16) << 3;

  /* Set up limit values for MV components to prevent them from
   * extending beyond the UMV borders assuming 16x16 block size */
  x->mv_row_min = -((mb_row * 16) + VP8BORDERINPIXELS - INTERP_EXTEND);
  x->mv_col_min = -((mb_col * 16) + VP8BORDERINPIXELS - INTERP_EXTEND);
  x->mv_row_max = ((cm->mb_rows - mb_row) * 16 +
                   (VP8BORDERINPIXELS - 32 - INTERP_EXTEND));
  x->mv_col_max = ((cm->mb_cols - mb_col) * 16 +
                   (VP8BORDERINPIXELS - 32 - INTERP_EXTEND));

  xd->up_available   = (mb_row != 0);
  xd->left_available = (mb_col != 0);

  recon_yoffset  = (mb_row * recon_y_stride * 16) + (mb_col * 16);
  recon_uvoffset = (mb_row * recon_uv_stride * 8) + (mb_col *  8);

  xd->dst.y_buffer = cm->yv12_fb[dst_fb_idx].y_buffer + recon_yoffset;
  xd->dst.u_buffer = cm->yv12_fb[dst_fb_idx].u_buffer + recon_uvoffset;
  xd->dst.v_buffer = cm->yv12_fb[dst_fb_idx].v_buffer + recon_uvoffset;
#if 0 // FIXME
  /* Copy current MB to a work buffer */
  vp8_copy_mem16x16(x->src.y_buffer, x->src.y_stride, x->thismb, 16);
#endif
  x->rddiv = cpi->RDDIV;
  x->rdmult = cpi->RDMULT;
  if(cpi->oxcf.tuning == VP8_TUNE_SSIM)
    vp8_activity_masking(cpi, x);
  /* Is segmentation enabled */
  if (xd->segmentation_enabled)
  {
    /* Code to set segment id in xd->mbmi.segment_id */
    if (xd->update_mb_segmentation_map)
      xd->mode_info_context->mbmi.segment_id =
            cpi->segmentation_map[map_index] &&
            cpi->segmentation_map[map_index + 1] &&
            cpi->segmentation_map[map_index + cm->mb_cols] &&
            cpi->segmentation_map[map_index + cm->mb_cols + 1];
    else
      xd->mode_info_context->mbmi.segment_id =
            cm->last_frame_seg_map[map_index] &&
            cm->last_frame_seg_map[map_index + 1] &&
            cm->last_frame_seg_map[map_index + cm->mb_cols] &&
            cm->last_frame_seg_map[map_index + cm->mb_cols + 1];
    if (xd->mode_info_context->mbmi.segment_id > 3)
      xd->mode_info_context->mbmi.segment_id = 0;

    vp8cx_mb_init_quantizer(cpi, x);
  }
  else
    /* Set to Segment 0 by default */
    xd->mode_info_context->mbmi.segment_id = 0;

  x->active_ptr = cpi->active_map + map_index;
  
  cpi->update_context = 0;    // TODO Do we need this now??

  /* Find best coding mode & reconstruct the MB so it is available
   * as a predictor for MBs that follow in the SB */
  if (cm->frame_type == KEY_FRAME)
  {
    vp8_rd_pick_intra_mode_sb(cpi, x,
                              totalrate,
                              totaldist);

    /* Save the coding context */
    vpx_memcpy(&x->sb_context[0].mic, xd->mode_info_context,
               sizeof(MODE_INFO));
  }
  else
  {
    if (xd->segmentation_enabled && cpi->seg0_cnt > 0 &&
        !segfeature_active( xd, 0, SEG_LVL_REF_FRAME ) &&
        segfeature_active( xd, 1, SEG_LVL_REF_FRAME ) &&
        check_segref(xd, 1, INTRA_FRAME)  +
        check_segref(xd, 1, LAST_FRAME)   +
        check_segref(xd, 1, GOLDEN_FRAME) +
        check_segref(xd, 1, ALTREF_FRAME) == 1)
    {
      cpi->seg0_progress = (cpi->seg0_idx << 16) / cpi->seg0_cnt;
    }
    else
    {
      cpi->seg0_progress =
        (((mb_col & ~1) * 2 + (mb_row & ~1) * cm->mb_cols) << 16) / cm->MBs;
    }

    vp8_rd_pick_inter_mode_sb(cpi, x,
                              recon_yoffset,
                              recon_uvoffset,
                              totalrate,
                              totaldist);
  }

  /* Restore L & A coding context to those in place on entry */
  vpx_memcpy (cm->left_context,
              left_context,
              sizeof(left_context));
  vpx_memcpy (initial_above_context_ptr,
              above_context,
              sizeof(above_context));
}
#endif

static void encode_sb(VP8_COMP *cpi,
                      VP8_COMMON *cm,
                      int mbrow,
                      int mbcol,
                      MACROBLOCK  *x,
                      MACROBLOCKD *xd,
                      TOKENEXTRA **tp) {
  int i;
  int map_index;
  int mb_row, mb_col;
  int recon_yoffset, recon_uvoffset;
  int ref_fb_idx = cm->lst_fb_idx;
  int dst_fb_idx = cm->new_fb_idx;
  int recon_y_stride = cm->yv12_fb[ref_fb_idx].y_stride;
  int recon_uv_stride = cm->yv12_fb[ref_fb_idx].uv_stride;
  int row_delta[4] = { 0, +1,  0, -1};
  int col_delta[4] = { +1, -1, +1, +1};

  mb_row = mbrow;
  mb_col = mbcol;

  /* Encode MBs in raster order within the SB */
  for (i = 0; i < 4; i++) {
    int dy = row_delta[i];
    int dx = col_delta[i];
    int offset_extended   = dy * xd->mode_info_stride + dx;
    int offset_unextended = dy * cm->mb_cols + dx;
    MB_MODE_INFO * mbmi = &xd->mode_info_context->mbmi;

    if ((mb_row >= cm->mb_rows) || (mb_col >= cm->mb_cols)) {
      // MB lies outside frame, move on
      mb_row += dy;
      mb_col += dx;

      x->src.y_buffer += 16 * (dx + dy * x->src.y_stride);
      x->src.u_buffer += 8  * (dx + dy * x->src.uv_stride);
      x->src.v_buffer += 8  * (dx + dy * x->src.uv_stride);

      x->gf_active_ptr      += offset_unextended;
      x->partition_info     += offset_extended;
      xd->mode_info_context += offset_extended;
      xd->prev_mode_info_context += offset_extended;

#if CONFIG_DEBUG
      assert((xd->prev_mode_info_context - cpi->common.prev_mip) ==
             (xd->mode_info_context - cpi->common.mip));
#endif
      continue;
    }

    xd->mb_index = i;

#ifdef ENC_DEBUG
    enc_debug = (cpi->common.current_video_frame == 0 &&
                 mb_row == 0 && mb_col == 0);
    mb_col_debug = mb_col;
    mb_row_debug = mb_row;
#endif

    // Restore MB state to that when it was picked
#if CONFIG_SUPERBLOCKS
    if (xd->mode_info_context->mbmi.encoded_as_sb) {
      update_state(cpi, x, &x->sb_context[i]);
      cpi->sb_count++;
    } else
#endif
      update_state(cpi, x, &x->mb_context[i]);

    map_index = (mb_row * cpi->common.mb_cols) + mb_col;
    x->mb_activity_ptr = &cpi->mb_activity_map[map_index];

    // reset above block coeffs
    xd->above_context = cm->above_context + mb_col;
    xd->left_context  = cm->left_context + (i >> 1);

    // Set up distance of MB to edge of the frame in 1/8th pel units
    xd->mb_to_top_edge    = -((mb_row * 16) << 3);
    xd->mb_to_left_edge   = -((mb_col * 16) << 3);
    xd->mb_to_bottom_edge = ((cm->mb_rows - 1 - mb_row) * 16) << 3;
    xd->mb_to_right_edge  = ((cm->mb_cols - 1 - mb_col) * 16) << 3;

#if CONFIG_SUPERBLOCKS
    if (xd->mode_info_context->mbmi.encoded_as_sb) {
      // Set up limit values for MV components to prevent them from
      // extending beyond the UMV borders assuming 32x32 block size
      x->mv_row_min = -((mb_row * 16) + VP8BORDERINPIXELS - INTERP_EXTEND);
      x->mv_col_min = -((mb_col * 16) + VP8BORDERINPIXELS - INTERP_EXTEND);
      x->mv_row_max = ((cm->mb_rows - mb_row) * 16 +
                       (VP8BORDERINPIXELS - 32 - INTERP_EXTEND));
      x->mv_col_max = ((cm->mb_cols - mb_col) * 16 +
                       (VP8BORDERINPIXELS - 32 - INTERP_EXTEND));
    } else {
#endif
      // Set up limit values for MV components to prevent them from
      // extending beyond the UMV borders assuming 16x16 block size
      x->mv_row_min = -((mb_row * 16) + VP8BORDERINPIXELS - INTERP_EXTEND);
      x->mv_col_min = -((mb_col * 16) + VP8BORDERINPIXELS - INTERP_EXTEND);
      x->mv_row_max = ((cm->mb_rows - mb_row) * 16 +
                       (VP8BORDERINPIXELS - 16 - INTERP_EXTEND));
      x->mv_col_max = ((cm->mb_cols - mb_col) * 16 +
                       (VP8BORDERINPIXELS - 16 - INTERP_EXTEND));
#if CONFIG_SUPERBLOCKS
    }
#endif

    xd->up_available = (mb_row != 0);
    xd->left_available = (mb_col != 0);

    recon_yoffset = (mb_row * recon_y_stride * 16) + (mb_col * 16);
    recon_uvoffset = (mb_row * recon_uv_stride * 8) + (mb_col * 8);

    xd->dst.y_buffer = cm->yv12_fb[dst_fb_idx].y_buffer + recon_yoffset;
    xd->dst.u_buffer = cm->yv12_fb[dst_fb_idx].u_buffer + recon_uvoffset;
    xd->dst.v_buffer = cm->yv12_fb[dst_fb_idx].v_buffer + recon_uvoffset;

    // Copy current MB to a work buffer
    vp8_copy_mem16x16(x->src.y_buffer, x->src.y_stride, x->thismb, 16);

    if (cpi->oxcf.tuning == VP8_TUNE_SSIM)
      vp8_activity_masking(cpi, x);

    // Is segmentation enabled
    if (xd->segmentation_enabled) {
      vp8cx_mb_init_quantizer(cpi, x);
    }

    x->active_ptr = cpi->active_map + map_index;

    cpi->update_context = 0;

    if (cm->frame_type == KEY_FRAME) {
#if CONFIG_SUPERBLOCKS
      if (xd->mode_info_context->mbmi.encoded_as_sb)
        vp8cx_encode_intra_super_block(cpi, x, tp, mb_col);
      else
#endif
        vp8cx_encode_intra_macro_block(cpi, x, tp, 1);
        // Note the encoder may have changed the segment_id

#ifdef MODE_STATS
      y_modes[mbmi->mode]++;
#endif
    } else {
      unsigned char *segment_id;
      int seg_ref_active;

      if (xd->mode_info_context->mbmi.ref_frame) {
        unsigned char pred_context;

        pred_context = get_pred_context(cm, xd, PRED_COMP);

        if (xd->mode_info_context->mbmi.second_ref_frame == INTRA_FRAME)
          cpi->single_pred_count[pred_context]++;
        else
          cpi->comp_pred_count[pred_context]++;
      }

#if CONFIG_SUPERBLOCKS
      if (xd->mode_info_context->mbmi.encoded_as_sb)
        vp8cx_encode_inter_superblock(cpi, x, tp, recon_yoffset, recon_uvoffset, mb_col, mb_row);
      else
#endif
        vp8cx_encode_inter_macroblock(cpi, x, tp,
                                      recon_yoffset, recon_uvoffset, 1);
        // Note the encoder may have changed the segment_id

#ifdef MODE_STATS
      inter_y_modes[mbmi->mode]++;

      if (mbmi->mode == SPLITMV) {
        int b;

        for (b = 0; b < x->partition_info->count; b++) {
          inter_b_modes[x->partition_info->bmi[b].mode]++;
        }
      }

#endif

      // If we have just a single reference frame coded for a segment then
      // exclude from the reference frame counts used to work out
      // probabilities. NOTE: At the moment we dont support custom trees
      // for the reference frame coding for each segment but this is a
      // possible future action.
      segment_id = &mbmi->segment_id;
      seg_ref_active = segfeature_active(xd, *segment_id, SEG_LVL_REF_FRAME);
      if (!seg_ref_active ||
          ((check_segref(xd, *segment_id, INTRA_FRAME) +
            check_segref(xd, *segment_id, LAST_FRAME) +
            check_segref(xd, *segment_id, GOLDEN_FRAME) +
            check_segref(xd, *segment_id, ALTREF_FRAME)) > 1)) {
        {
          cpi->count_mb_ref_frame_usage[mbmi->ref_frame]++;
        }
      }

      // Count of last ref frame 0,0 usage
      if ((mbmi->mode == ZEROMV) && (mbmi->ref_frame == LAST_FRAME))
        cpi->inter_zz_count++;
    }

#if CONFIG_SUPERBLOCKS
    if (xd->mode_info_context->mbmi.encoded_as_sb) {
      x->src.y_buffer += 32;
      x->src.u_buffer += 16;
      x->src.v_buffer += 16;

      x->gf_active_ptr      += 2;
      x->partition_info     += 2;
      xd->mode_info_context += 2;
      xd->prev_mode_info_context += 2;

      (*tp)->Token = EOSB_TOKEN;
      (*tp)++;
      if (mb_row < cm->mb_rows) cpi->tplist[mb_row].stop = *tp;
      break;
    }
#endif

    // Next MB
    mb_row += dy;
    mb_col += dx;

    x->src.y_buffer += 16 * (dx + dy * x->src.y_stride);
    x->src.u_buffer += 8  * (dx + dy * x->src.uv_stride);
    x->src.v_buffer += 8  * (dx + dy * x->src.uv_stride);

    x->gf_active_ptr      += offset_unextended;
    x->partition_info     += offset_extended;
    xd->mode_info_context += offset_extended;
    xd->prev_mode_info_context += offset_extended;

#if CONFIG_DEBUG
    assert((xd->prev_mode_info_context - cpi->common.prev_mip) ==
           (xd->mode_info_context - cpi->common.mip));
#endif
    (*tp)->Token = EOSB_TOKEN;
    (*tp)++;
    if (mb_row < cm->mb_rows) cpi->tplist[mb_row].stop = *tp;
  }

  // debug output
#if DBG_PRNT_SEGMAP
  {
    FILE *statsfile;
    statsfile = fopen("segmap2.stt", "a");
    fprintf(statsfile, "\n");
    fclose(statsfile);
  }
#endif
}

static
void encode_sb_row(VP8_COMP *cpi,
                   VP8_COMMON *cm,
                   int mb_row,
                   MACROBLOCK  *x,
                   MACROBLOCKD *xd,
                   TOKENEXTRA **tp,
                   int *totalrate) {
  int mb_col;
  int mb_cols = cm->mb_cols;

  // Initialize the left context for the new SB row
  vpx_memset(cm->left_context, 0, sizeof(cm->left_context));

  // Code each SB in the row
  for (mb_col = 0; mb_col < mb_cols; mb_col += 2) {
    int mb_rate = 0, mb_dist = 0;
#if CONFIG_SUPERBLOCKS
    int sb_rate = INT_MAX, sb_dist;
#endif

#if CONFIG_DEBUG
    MODE_INFO *mic = xd->mode_info_context;
    PARTITION_INFO *pi = x->partition_info;
    signed char  *gfa = x->gf_active_ptr;
    unsigned char *yb = x->src.y_buffer;
    unsigned char *ub = x->src.u_buffer;
    unsigned char *vb = x->src.v_buffer;
#endif

#if CONFIG_SUPERBLOCKS
    // Pick modes assuming the SB is coded as 4 independent MBs
    xd->mode_info_context->mbmi.encoded_as_sb = 0;
#endif
    pick_mb_modes(cpi, cm, mb_row, mb_col, x, xd, tp, &mb_rate, &mb_dist);
#if CONFIG_SUPERBLOCKS
    mb_rate += vp8_cost_bit(cm->sb_coded, 0);
#endif

    x->src.y_buffer -= 32;
    x->src.u_buffer -= 16;
    x->src.v_buffer -= 16;

    x->gf_active_ptr -= 2;
    x->partition_info -= 2;
    xd->mode_info_context -= 2;
    xd->prev_mode_info_context -= 2;

#if CONFIG_DEBUG
    assert(x->gf_active_ptr == gfa);
    assert(x->partition_info == pi);
    assert(xd->mode_info_context == mic);
    assert(x->src.y_buffer == yb);
    assert(x->src.u_buffer == ub);
    assert(x->src.v_buffer == vb);
#endif

#if CONFIG_SUPERBLOCKS
    if (!(((    mb_cols & 1) && mb_col ==     mb_cols - 1) ||
          ((cm->mb_rows & 1) && mb_row == cm->mb_rows - 1))) {
      /* Pick a mode assuming that it applies to all 4 of the MBs in the SB */
      xd->mode_info_context->mbmi.encoded_as_sb = 1;
      pick_sb_modes(cpi, cm, mb_row, mb_col, x, xd, tp, &sb_rate, &sb_dist);
      sb_rate += vp8_cost_bit(cm->sb_coded, 1);
    }

    /* Decide whether to encode as a SB or 4xMBs */
    if (sb_rate < INT_MAX &&
        RDCOST(x->rdmult, x->rddiv, sb_rate, sb_dist) <
          RDCOST(x->rdmult, x->rddiv, mb_rate, mb_dist)) {
      xd->mode_info_context->mbmi.encoded_as_sb = 1;
      xd->mode_info_context[1].mbmi.encoded_as_sb = 1;
      xd->mode_info_context[cm->mode_info_stride].mbmi.encoded_as_sb = 1;
      xd->mode_info_context[1 + cm->mode_info_stride].mbmi.encoded_as_sb = 1;
      *totalrate += sb_rate;
    } else
#endif
    {
#if CONFIG_SUPERBLOCKS
      xd->mode_info_context->mbmi.encoded_as_sb = 0;
      if (cm->mb_cols - 1 > mb_col)
        xd->mode_info_context[1].mbmi.encoded_as_sb = 0;
      if (cm->mb_rows - 1 > mb_row) {
        xd->mode_info_context[cm->mode_info_stride].mbmi.encoded_as_sb = 0;
        if (cm->mb_cols - 1 > mb_col)
          xd->mode_info_context[1 + cm->mode_info_stride].mbmi.encoded_as_sb = 0;
      }
#endif
      *totalrate += mb_rate;
    }

    /* Encode SB using best computed mode(s) */
    encode_sb(cpi, cm, mb_row, mb_col, x, xd, tp);

#if CONFIG_DEBUG
    assert(x->gf_active_ptr == gfa + 2);
    assert(x->partition_info == pi + 2);
    assert(xd->mode_info_context == mic + 2);
    assert(x->src.y_buffer == yb + 32);
    assert(x->src.u_buffer == ub + 16);
    assert(x->src.v_buffer == vb + 16);
#endif
  }

  // this is to account for the border
  x->gf_active_ptr += mb_cols - (mb_cols & 0x1);
  x->partition_info += xd->mode_info_stride + 1 - (mb_cols & 0x1);
  xd->mode_info_context += xd->mode_info_stride + 1 - (mb_cols & 0x1);
  xd->prev_mode_info_context += xd->mode_info_stride + 1 - (mb_cols & 0x1);

#if CONFIG_DEBUG
  assert((xd->prev_mode_info_context - cpi->common.prev_mip) ==
         (xd->mode_info_context - cpi->common.mip));
#endif
}

void init_encode_frame_mb_context(VP8_COMP *cpi) {
  MACROBLOCK *const x = &cpi->mb;
  VP8_COMMON *const cm = &cpi->common;
  MACROBLOCKD *const xd = &x->e_mbd;

  // GF active flags data structure
  x->gf_active_ptr = (signed char *)cpi->gf_active_flags;

  // Activity map pointer
  x->mb_activity_ptr = cpi->mb_activity_map;

  x->act_zbin_adj = 0;
  cpi->seg0_idx = 0;
  vpx_memset(cpi->ref_pred_count, 0, sizeof(cpi->ref_pred_count));

  x->partition_info = x->pi;

  xd->mode_info_context = cm->mi;
  xd->mode_info_stride = cm->mode_info_stride;
  xd->prev_mode_info_context = cm->prev_mi;

  xd->frame_type = cm->frame_type;

  xd->frames_since_golden = cm->frames_since_golden;
  xd->frames_till_alt_ref_frame = cm->frames_till_alt_ref_frame;

  // reset intra mode contexts
  if (cm->frame_type == KEY_FRAME)
    vp8_init_mbmode_probs(cm);

  // Copy data over into macro block data structures.
  x->src = * cpi->Source;
  xd->pre = cm->yv12_fb[cm->lst_fb_idx];
  xd->dst = cm->yv12_fb[cm->new_fb_idx];

  // set up frame for intra coded blocks
  vp8_setup_intra_recon(&cm->yv12_fb[cm->new_fb_idx]);

  vp8_build_block_offsets(x);

  vp8_setup_block_dptrs(&x->e_mbd);

  vp8_setup_block_ptrs(x);

  xd->mode_info_context->mbmi.mode = DC_PRED;
  xd->mode_info_context->mbmi.uv_mode = DC_PRED;

  vp8_zero(cpi->count_mb_ref_frame_usage)
  vp8_zero(cpi->bmode_count)
  vp8_zero(cpi->ymode_count)
  vp8_zero(cpi->i8x8_mode_count)
  vp8_zero(cpi->y_uv_mode_count)
  vp8_zero(cpi->sub_mv_ref_count)
  vp8_zero(cpi->mbsplit_count)
  vp8_zero(cpi->common.fc.mv_ref_ct)
  vp8_zero(cpi->common.fc.mv_ref_ct_a)
#if CONFIG_SUPERBLOCKS
  vp8_zero(cpi->sb_ymode_count)
  cpi->sb_count = 0;
#endif
  // vp8_zero(cpi->uv_mode_count)

  vpx_memset(cm->above_context, 0,
             sizeof(ENTROPY_CONTEXT_PLANES) * cm->mb_cols);

  xd->fullpixel_mask = 0xffffffff;
  if (cm->full_pixel)
    xd->fullpixel_mask = 0xfffffff8;
}

static void encode_frame_internal(VP8_COMP *cpi) {
  int mb_row;
  MACROBLOCK *const x = &cpi->mb;
  VP8_COMMON *const cm = &cpi->common;
  MACROBLOCKD *const xd = &x->e_mbd;

  TOKENEXTRA *tp = cpi->tok;
  int totalrate;

  //printf("encode_frame_internal\n");

  // Compute a modified set of reference frame probabilities to use when
  // prediction fails. These are based on the current general estimates for
  // this frame which may be updated with each iteration of the recode loop.
  compute_mod_refprobs(cm);

#if CONFIG_NEW_MVREF
  // temp stats reset
  vp8_zero( cpi->best_ref_index_counts );
#endif

// debug output
#if DBG_PRNT_SEGMAP
  {
    FILE *statsfile;
    statsfile = fopen("segmap2.stt", "a");
    fprintf(statsfile, "\n");
    fclose(statsfile);
  }
#endif

  totalrate = 0;

  // Functions setup for all frame types so we can use MC in AltRef
  vp8_setup_interp_filters(xd, cm->mcomp_filter_type, cm);

  // Reset frame count of inter 0,0 motion vector usage.
  cpi->inter_zz_count = 0;

  cpi->prediction_error = 0;
  cpi->intra_error = 0;
  cpi->skip_true_count[0] = cpi->skip_true_count[1] = cpi->skip_true_count[2] = 0;
  cpi->skip_false_count[0] = cpi->skip_false_count[1] = cpi->skip_false_count[2] = 0;

#if CONFIG_PRED_FILTER
  if (cm->current_video_frame == 0) {
    // Initially assume that we'll signal the prediction filter
    // state at the frame level and that it is off.
    cpi->common.pred_filter_mode = 0;
    cpi->common.prob_pred_filter_off = 128;
  }
  cpi->pred_filter_on_count = 0;
  cpi->pred_filter_off_count = 0;
#endif
#if CONFIG_SWITCHABLE_INTERP
  vp8_zero(cpi->switchable_interp_count);
#endif

#if 0
  // Experimental code
  cpi->frame_distortion = 0;
  cpi->last_mb_distortion = 0;
#endif

  xd->mode_info_context = cm->mi;
  xd->prev_mode_info_context = cm->prev_mi;

  vp8_zero(cpi->NMVcount);
  vp8_zero(cpi->coef_counts);
  vp8_zero(cpi->hybrid_coef_counts);
  vp8_zero(cpi->coef_counts_8x8);
  vp8_zero(cpi->hybrid_coef_counts_8x8);
  vp8_zero(cpi->coef_counts_16x16);
  vp8_zero(cpi->hybrid_coef_counts_16x16);

  vp8cx_frame_init_quantizer(cpi);

  vp8_initialize_rd_consts(cpi, cm->base_qindex + cm->y1dc_delta_q);
  vp8cx_initialize_me_consts(cpi, cm->base_qindex);

  if (cpi->oxcf.tuning == VP8_TUNE_SSIM) {
    // Initialize encode frame context.
    init_encode_frame_mb_context(cpi);

    // Build a frame level activity map
    build_activity_map(cpi);
  }

  // re-initencode frame context.
  init_encode_frame_mb_context(cpi);

  vpx_memset(cpi->rd_comp_pred_diff, 0, sizeof(cpi->rd_comp_pred_diff));
  vpx_memset(cpi->single_pred_count, 0, sizeof(cpi->single_pred_count));
  vpx_memset(cpi->comp_pred_count, 0, sizeof(cpi->comp_pred_count));
  vpx_memset(cpi->txfm_count, 0, sizeof(cpi->txfm_count));
  vpx_memset(cpi->txfm_count_8x8p, 0, sizeof(cpi->txfm_count_8x8p));
  vpx_memset(cpi->rd_tx_select_diff, 0, sizeof(cpi->rd_tx_select_diff));
  {
    struct vpx_usec_timer  emr_timer;
    vpx_usec_timer_start(&emr_timer);

    {
      // For each row of SBs in the frame
      for (mb_row = 0; mb_row < cm->mb_rows; mb_row += 2) {
        int offset = (cm->mb_cols + 1) & ~0x1;

        encode_sb_row(cpi, cm, mb_row, x, xd, &tp, &totalrate);

        // adjust to the next row of SBs
        x->src.y_buffer += 32 * x->src.y_stride - 16 * offset;
        x->src.u_buffer += 16 * x->src.uv_stride - 8 * offset;
        x->src.v_buffer += 16 * x->src.uv_stride - 8 * offset;
      }

      cpi->tok_count = tp - cpi->tok;
    }

    vpx_usec_timer_mark(&emr_timer);
    cpi->time_encode_mb_row += vpx_usec_timer_elapsed(&emr_timer);

  }

  // 256 rate units to the bit,
  // projected_frame_size in units of BYTES
  cpi->projected_frame_size = totalrate >> 8;


#if 0
  // Keep record of the total distortion this time around for future use
  cpi->last_frame_distortion = cpi->frame_distortion;
#endif

}

static int check_dual_ref_flags(VP8_COMP *cpi) {
  MACROBLOCKD *xd = &cpi->mb.e_mbd;
  int ref_flags = cpi->ref_frame_flags;

  if (segfeature_active(xd, 1, SEG_LVL_REF_FRAME)) {
    if ((ref_flags & (VP8_LAST_FLAG | VP8_GOLD_FLAG)) == (VP8_LAST_FLAG | VP8_GOLD_FLAG) &&
        check_segref(xd, 1, LAST_FRAME))
      return 1;
    if ((ref_flags & (VP8_GOLD_FLAG | VP8_ALT_FLAG)) == (VP8_GOLD_FLAG | VP8_ALT_FLAG) &&
        check_segref(xd, 1, GOLDEN_FRAME))
      return 1;
    if ((ref_flags & (VP8_ALT_FLAG  | VP8_LAST_FLAG)) == (VP8_ALT_FLAG  | VP8_LAST_FLAG) &&
        check_segref(xd, 1, ALTREF_FRAME))
      return 1;
    return 0;
  } else {
    return (!!(ref_flags & VP8_GOLD_FLAG) +
            !!(ref_flags & VP8_LAST_FLAG) +
            !!(ref_flags & VP8_ALT_FLAG)) >= 2;
  }
}

static void reset_skip_txfm_size(VP8_COMP *cpi, TX_SIZE txfm_max) {
  VP8_COMMON *cm = &cpi->common;
  int mb_row, mb_col, mis = cm->mode_info_stride;
  MODE_INFO *mi, *mi_ptr = cm->mi;
  MB_MODE_INFO *mbmi;
  MACROBLOCK *x = &cpi->mb;
  MACROBLOCKD *xd = &x->e_mbd;

  for (mb_row = 0; mb_row < cm->mb_rows; mb_row++, mi_ptr += mis) {
    mi = mi_ptr;
    for (mb_col = 0; mb_col < cm->mb_cols; mb_col++, mi++) {
      mbmi = &mi->mbmi;
      if (mbmi->txfm_size > txfm_max) {
        int segment_id = mbmi->segment_id;
        xd->mode_info_context = mi;
        assert((segfeature_active(xd, segment_id, SEG_LVL_EOB) &&
                get_segdata(xd, segment_id, SEG_LVL_EOB) == 0) ||
               (cm->mb_no_coeff_skip && mbmi->mb_skip_coeff));
        mbmi->txfm_size = txfm_max;
      }
    }
  }
}

void vp8_encode_frame(VP8_COMP *cpi) {
  if (cpi->sf.RD) {
    int i, frame_type, pred_type;
    TXFM_MODE txfm_type;

    /*
     * This code does a single RD pass over the whole frame assuming
     * either compound, single or hybrid prediction as per whatever has
     * worked best for that type of frame in the past.
     * It also predicts whether another coding mode would have worked
     * better that this coding mode. If that is the case, it remembers
     * that for subsequent frames.
     * It does the same analysis for transform size selection also.
     */
    if (cpi->common.frame_type == KEY_FRAME)
      frame_type = 0;
    else if (cpi->is_src_frame_alt_ref && cpi->common.refresh_golden_frame)
      frame_type = 3;
    else if (cpi->common.refresh_golden_frame || cpi->common.refresh_alt_ref_frame)
      frame_type = 1;
    else
      frame_type = 2;

    /* prediction (compound, single or hybrid) mode selection */
    if (frame_type == 3)
      pred_type = SINGLE_PREDICTION_ONLY;
    else if (cpi->rd_prediction_type_threshes[frame_type][1] >
                 cpi->rd_prediction_type_threshes[frame_type][0] &&
             cpi->rd_prediction_type_threshes[frame_type][1] >
                 cpi->rd_prediction_type_threshes[frame_type][2] &&
             check_dual_ref_flags(cpi) && cpi->static_mb_pct == 100)
      pred_type = COMP_PREDICTION_ONLY;
    else if (cpi->rd_prediction_type_threshes[frame_type][0] >
                 cpi->rd_prediction_type_threshes[frame_type][2])
      pred_type = SINGLE_PREDICTION_ONLY;
    else
      pred_type = HYBRID_PREDICTION;

    /* transform size (4x4, 8x8, 16x16 or select-per-mb) selection */
#if CONFIG_LOSSLESS
    if (cpi->oxcf.lossless) {
      txfm_type = ONLY_4X4;
    } else
#endif
    /* FIXME (rbultje)
     * this is a hack (no really), basically to work around the complete
     * nonsense coefficient cost prediction for keyframes. The probabilities
     * are reset to defaults, and thus we basically have no idea how expensive
     * a 4x4 vs. 8x8 will really be. The result is that any estimate at which
     * of the two is better is utterly bogus.
     * I'd like to eventually remove this hack, but in order to do that, we
     * need to move the frame reset code from the frame encode init to the
     * bitstream write code, or alternatively keep a backup of the previous
     * keyframe's probabilities as an estimate of what the current keyframe's
     * coefficient cost distributions may look like. */
    if (frame_type == 0) {
      txfm_type = ALLOW_16X16;
    } else
#if 0
    /* FIXME (rbultje)
     * this code is disabled for a similar reason as the code above; the
     * problem is that each time we "revert" to 4x4 only (or even 8x8 only),
     * the coefficient probabilities for 16x16 (and 8x8) start lagging behind,
     * thus leading to them lagging further behind and not being chosen for
     * subsequent frames either. This is essentially a local minimum problem
     * that we can probably fix by estimating real costs more closely within
     * a frame, perhaps by re-calculating costs on-the-fly as frame encoding
     * progresses. */
    if (cpi->rd_tx_select_threshes[frame_type][TX_MODE_SELECT] >
            cpi->rd_tx_select_threshes[frame_type][ONLY_4X4] &&
        cpi->rd_tx_select_threshes[frame_type][TX_MODE_SELECT] >
            cpi->rd_tx_select_threshes[frame_type][ALLOW_16X16] &&
        cpi->rd_tx_select_threshes[frame_type][TX_MODE_SELECT] >
            cpi->rd_tx_select_threshes[frame_type][ALLOW_8X8]) {
      txfm_type = TX_MODE_SELECT;
    } else if (cpi->rd_tx_select_threshes[frame_type][ONLY_4X4] >
                  cpi->rd_tx_select_threshes[frame_type][ALLOW_8X8]
            && cpi->rd_tx_select_threshes[frame_type][ONLY_4X4] >
                  cpi->rd_tx_select_threshes[frame_type][ALLOW_16X16]
               ) {
      txfm_type = ONLY_4X4;
    } else if (cpi->rd_tx_select_threshes[frame_type][ALLOW_16X16] >=
                  cpi->rd_tx_select_threshes[frame_type][ALLOW_8X8]) {
      txfm_type = ALLOW_16X16;
    } else
      txfm_type = ALLOW_8X8;
#else
    txfm_type = cpi->rd_tx_select_threshes[frame_type][ALLOW_16X16] >=
                 cpi->rd_tx_select_threshes[frame_type][TX_MODE_SELECT] ?
    ALLOW_16X16 : TX_MODE_SELECT;
#endif
    cpi->common.txfm_mode = txfm_type;
    if (txfm_type != TX_MODE_SELECT) {
      cpi->common.prob_tx[0] = 128;
      cpi->common.prob_tx[1] = 128;
    }
    cpi->common.comp_pred_mode = pred_type;
    encode_frame_internal(cpi);

    for (i = 0; i < NB_PREDICTION_TYPES; ++i) {
      const int diff = cpi->rd_comp_pred_diff[i] / cpi->common.MBs;
      cpi->rd_prediction_type_threshes[frame_type][i] += diff;
      cpi->rd_prediction_type_threshes[frame_type][i] >>= 1;
    }

    for (i = 0; i < NB_TXFM_MODES; ++i) {
      int64_t pd = cpi->rd_tx_select_diff[i];
      int diff;
      if (i == TX_MODE_SELECT)
        pd -= RDCOST(cpi->mb.rdmult, cpi->mb.rddiv, 2048 * (TX_SIZE_MAX - 1), 0);
      diff = pd / cpi->common.MBs;
      cpi->rd_tx_select_threshes[frame_type][i] += diff;
      cpi->rd_tx_select_threshes[frame_type][i] /= 2;
    }

    if (cpi->common.comp_pred_mode == HYBRID_PREDICTION) {
      int single_count_zero = 0;
      int comp_count_zero = 0;

      for (i = 0; i < COMP_PRED_CONTEXTS; i++) {
        single_count_zero += cpi->single_pred_count[i];
        comp_count_zero += cpi->comp_pred_count[i];
      }

      if (comp_count_zero == 0) {
        cpi->common.comp_pred_mode = SINGLE_PREDICTION_ONLY;
      } else if (single_count_zero == 0) {
        cpi->common.comp_pred_mode = COMP_PREDICTION_ONLY;
      }
    }

    if (cpi->common.txfm_mode == TX_MODE_SELECT) {
      const int count4x4 = cpi->txfm_count[TX_4X4] + cpi->txfm_count_8x8p[TX_4X4];
      const int count8x8 = cpi->txfm_count[TX_8X8];
      const int count8x8_8x8p = cpi->txfm_count_8x8p[TX_8X8];
      const int count16x16 = cpi->txfm_count[TX_16X16];

      if (count4x4 == 0 && count16x16 == 0) {
        cpi->common.txfm_mode = ALLOW_8X8;
        reset_skip_txfm_size(cpi, TX_8X8);
      } else if (count8x8 == 0 && count16x16 == 0 && count8x8_8x8p == 0) {
        cpi->common.txfm_mode = ONLY_4X4;
        reset_skip_txfm_size(cpi, TX_4X4);
      } else if (count8x8 == 0 && count4x4 == 0) {
        cpi->common.txfm_mode = ALLOW_16X16;
      }
    }
  } else {
    encode_frame_internal(cpi);
  }

}

void vp8_setup_block_ptrs(MACROBLOCK *x) {
  int r, c;
  int i;

  for (r = 0; r < 4; r++) {
    for (c = 0; c < 4; c++) {
      x->block[r * 4 + c].src_diff = x->src_diff + r * 4 * 16 + c * 4;
    }
  }

  for (r = 0; r < 2; r++) {
    for (c = 0; c < 2; c++) {
      x->block[16 + r * 2 + c].src_diff = x->src_diff + 256 + r * 4 * 8 + c * 4;
    }
  }


  for (r = 0; r < 2; r++) {
    for (c = 0; c < 2; c++) {
      x->block[20 + r * 2 + c].src_diff = x->src_diff + 320 + r * 4 * 8 + c * 4;
    }
  }

  x->block[24].src_diff = x->src_diff + 384;


  for (i = 0; i < 25; i++) {
    x->block[i].coeff = x->coeff + i * 16;
  }
}

void vp8_build_block_offsets(MACROBLOCK *x) {
  int block = 0;
  int br, bc;

  vp8_build_block_doffsets(&x->e_mbd);

  // y blocks
  x->thismb_ptr = &x->thismb[0];
  for (br = 0; br < 4; br++) {
    for (bc = 0; bc < 4; bc++) {
      BLOCK *this_block = &x->block[block];
      // this_block->base_src = &x->src.y_buffer;
      // this_block->src_stride = x->src.y_stride;
      // this_block->src = 4 * br * this_block->src_stride + 4 * bc;
      this_block->base_src = &x->thismb_ptr;
      this_block->src_stride = 16;
      this_block->src = 4 * br * 16 + 4 * bc;
      ++block;
    }
  }

  // u blocks
  for (br = 0; br < 2; br++) {
    for (bc = 0; bc < 2; bc++) {
      BLOCK *this_block = &x->block[block];
      this_block->base_src = &x->src.u_buffer;
      this_block->src_stride = x->src.uv_stride;
      this_block->src = 4 * br * this_block->src_stride + 4 * bc;
      ++block;
    }
  }

  // v blocks
  for (br = 0; br < 2; br++) {
    for (bc = 0; bc < 2; bc++) {
      BLOCK *this_block = &x->block[block];
      this_block->base_src = &x->src.v_buffer;
      this_block->src_stride = x->src.uv_stride;
      this_block->src = 4 * br * this_block->src_stride + 4 * bc;
      ++block;
    }
  }
}

static void sum_intra_stats(VP8_COMP *cpi, MACROBLOCK *x) {
  const MACROBLOCKD *xd = &x->e_mbd;
  const MB_PREDICTION_MODE m = xd->mode_info_context->mbmi.mode;
  const MB_PREDICTION_MODE uvm = xd->mode_info_context->mbmi.uv_mode;

#ifdef MODE_STATS
  const int is_key = cpi->common.frame_type == KEY_FRAME;

  ++ (is_key ? uv_modes : inter_uv_modes)[uvm];
  ++ uv_modes_y[m][uvm];

  if (m == B_PRED) {
    unsigned int *const bct = is_key ? b_modes : inter_b_modes;

    int b = 0;

    do {
      ++ bct[xd->block[b].bmi.as_mode.first];
    } while (++b < 16);
  }

  if (m == I8X8_PRED) {
    i8x8_modes[xd->block[0].bmi.as_mode.first]++;
    i8x8_modes[xd->block[2].bmi.as_mode.first]++;
    i8x8_modes[xd->block[8].bmi.as_mode.first]++;
    i8x8_modes[xd->block[10].bmi.as_mode.first]++;
  }
#endif

#if CONFIG_SUPERBLOCKS
  if (xd->mode_info_context->mbmi.encoded_as_sb) {
    ++cpi->sb_ymode_count[m];
  } else
#endif
    ++cpi->ymode_count[m];
  if (m != I8X8_PRED)
    ++cpi->y_uv_mode_count[m][uvm];
  else {
    cpi->i8x8_mode_count[xd->block[0].bmi.as_mode.first]++;
    cpi->i8x8_mode_count[xd->block[2].bmi.as_mode.first]++;
    cpi->i8x8_mode_count[xd->block[8].bmi.as_mode.first]++;
    cpi->i8x8_mode_count[xd->block[10].bmi.as_mode.first]++;
  }
  if (m == B_PRED) {
    int b = 0;
    do {
      ++ cpi->bmode_count[xd->block[b].bmi.as_mode.first];
    } while (++b < 16);
  }
}

// Experimental stub function to create a per MB zbin adjustment based on
// some previously calculated measure of MB activity.
static void adjust_act_zbin(VP8_COMP *cpi, MACROBLOCK *x) {
#if USE_ACT_INDEX
  x->act_zbin_adj = *(x->mb_activity_ptr);
#else
  int64_t a;
  int64_t b;
  int64_t act = *(x->mb_activity_ptr);

  // Apply the masking to the RD multiplier.
  a = act + 4 * cpi->activity_avg;
  b = 4 * act + cpi->activity_avg;

  if (act > cpi->activity_avg)
    x->act_zbin_adj = (int)(((int64_t)b + (a >> 1)) / a) - 1;
  else
    x->act_zbin_adj = 1 - (int)(((int64_t)a + (b >> 1)) / b);
#endif
}

#if CONFIG_SUPERBLOCKS
static void update_sb_skip_coeff_state(VP8_COMP *cpi,
                                       MACROBLOCK *x,
                                       ENTROPY_CONTEXT_PLANES ta[4],
                                       ENTROPY_CONTEXT_PLANES tl[4],
                                       TOKENEXTRA *t[4],
                                       TOKENEXTRA **tp,
                                       int skip[4])
{
  TOKENEXTRA tokens[4][16 * 24];
  int n_tokens[4], n;

  // if there were no skips, we don't need to do anything
  if (!skip[0] && !skip[1] && !skip[2] && !skip[3])
    return;

  // if we don't do coeff skipping for this frame, we don't
  // need to do anything here
  if (!cpi->common.mb_no_coeff_skip)
    return;

  // if all 4 MBs skipped coeff coding, nothing to be done
  if (skip[0] && skip[1] && skip[2] && skip[3])
    return;

  // so the situation now is that we want to skip coeffs
  // for some MBs, but not all, and we didn't code EOB
  // coefficients for them. However, the skip flag for this
  // SB will be 0 overall, so we need to insert EOBs in the
  // middle of the token tree. Do so here.
  n_tokens[0] = t[1] - t[0];
  n_tokens[1] = t[2] - t[1];
  n_tokens[2] = t[3] - t[2];
  n_tokens[3] = *tp  - t[3];
  if (n_tokens[0])
    memcpy(tokens[0], t[0], n_tokens[0] * sizeof(*t[0]));
  if (n_tokens[1])
    memcpy(tokens[1], t[1], n_tokens[1] * sizeof(*t[0]));
  if (n_tokens[2])
    memcpy(tokens[2], t[2], n_tokens[2] * sizeof(*t[0]));
  if (n_tokens[3])
    memcpy(tokens[3], t[3], n_tokens[3] * sizeof(*t[0]));

  // reset pointer, stuff EOBs where necessary
  *tp = t[0];
  for (n = 0; n < 4; n++) {
    if (skip[n]) {
      x->e_mbd.above_context = &ta[n];
      x->e_mbd.left_context  = &tl[n];
      vp8_stuff_mb(cpi, &x->e_mbd, tp, 0);
    } else {
      if (n_tokens[n]) {
        memcpy(*tp, tokens[n], sizeof(*t[0]) * n_tokens[n]);
      }
      (*tp) += n_tokens[n];
    }
  }
}

void vp8cx_encode_intra_super_block(VP8_COMP *cpi,
                                    MACROBLOCK *x,
                                    TOKENEXTRA **t,
                                    int mb_col) {
  const int output_enabled = 1;
  int n;
  MACROBLOCKD *xd = &x->e_mbd;
  VP8_COMMON *cm = &cpi->common;
  const uint8_t *src = x->src.y_buffer;
  uint8_t *dst = xd->dst.y_buffer;
  const uint8_t *usrc = x->src.u_buffer;
  uint8_t *udst = xd->dst.u_buffer;
  const uint8_t *vsrc = x->src.v_buffer;
  uint8_t *vdst = xd->dst.v_buffer;
  int src_y_stride = x->src.y_stride, dst_y_stride = xd->dst.y_stride;
  int src_uv_stride = x->src.uv_stride, dst_uv_stride = xd->dst.uv_stride;
  const VP8_ENCODER_RTCD *rtcd = IF_RTCD(&cpi->rtcd);
  TOKENEXTRA *tp[4];
  int skip[4];
  MODE_INFO *mi = x->e_mbd.mode_info_context;
  ENTROPY_CONTEXT_PLANES ta[4], tl[4];

  if ((cpi->oxcf.tuning == VP8_TUNE_SSIM) && output_enabled) {
    adjust_act_zbin(cpi, x);
    vp8_update_zbin_extra(cpi, x);
  }

  vp8_build_intra_predictors_sby_s(&x->e_mbd);
  vp8_build_intra_predictors_sbuv_s(&x->e_mbd);

  assert(x->e_mbd.mode_info_context->mbmi.txfm_size == TX_8X8);
  for (n = 0; n < 4; n++)
  {
    int x_idx = n & 1, y_idx = n >> 1;

    xd->above_context = cm->above_context + mb_col + (n & 1);
    xd->left_context = cm->left_context + (n >> 1);

    vp8_subtract_mby_s_c(x->src_diff,
                         src + x_idx * 16 + y_idx * 16 * src_y_stride,
                         src_y_stride,
                         dst + x_idx * 16 + y_idx * 16 * dst_y_stride,
                         dst_y_stride);
    vp8_subtract_mbuv_s_c(x->src_diff,
                          usrc + x_idx * 8 + y_idx * 8 * src_uv_stride,
                          vsrc + x_idx * 8 + y_idx * 8 * src_uv_stride,
                          src_uv_stride,
                          udst + x_idx * 8 + y_idx * 8 * dst_uv_stride,
                          vdst + x_idx * 8 + y_idx * 8 * dst_uv_stride,
                          dst_uv_stride);
    vp8_transform_mb_8x8(x);
    vp8_quantize_mb_8x8(x);
    if (x->optimize) {
      vp8_optimize_mby_8x8(x, rtcd);
      vp8_optimize_mbuv_8x8(x, rtcd);
    }
    vp8_inverse_transform_mb_8x8(IF_RTCD(&rtcd->common->idct), &x->e_mbd);
    vp8_recon_mby_s_c(&x->e_mbd, dst + x_idx * 16 + y_idx * 16 * dst_y_stride);
    vp8_recon_mbuv_s_c(&x->e_mbd,
                       udst + x_idx * 8 + y_idx * 8 * dst_uv_stride,
                       vdst + x_idx * 8 + y_idx * 8 * dst_uv_stride);

    if (output_enabled) {
      memcpy(&ta[n], xd->above_context, sizeof(ta[n]));
      memcpy(&tl[n], xd->left_context, sizeof(tl[n]));
      tp[n] = *t;
      xd->mode_info_context = mi + x_idx + y_idx * cm->mode_info_stride;
      vp8_tokenize_mb(cpi, &x->e_mbd, t, 0);
      skip[n] = xd->mode_info_context->mbmi.mb_skip_coeff;
    }
  }

  if (output_enabled) {
    // Tokenize
    xd->mode_info_context = mi;
    sum_intra_stats(cpi, x);
    update_sb_skip_coeff_state(cpi, x, ta, tl, tp, t, skip);
  }
}
#endif /* CONFIG_SUPERBLOCKS */

void vp8cx_encode_intra_macro_block(VP8_COMP *cpi,
                                    MACROBLOCK *x,
                                    TOKENEXTRA **t,
                                    int output_enabled) {
  MB_MODE_INFO * mbmi = &x->e_mbd.mode_info_context->mbmi;
  if ((cpi->oxcf.tuning == VP8_TUNE_SSIM) && output_enabled) {
    adjust_act_zbin(cpi, x);
    vp8_update_zbin_extra(cpi, x);
  }
  if (mbmi->mode == I8X8_PRED) {
    vp8_encode_intra8x8mby(IF_RTCD(&cpi->rtcd), x);
    vp8_encode_intra8x8mbuv(IF_RTCD(&cpi->rtcd), x);
  } else if (mbmi->mode == B_PRED) {
    vp8_intra_prediction_down_copy(&x->e_mbd);
    vp8_encode_intra4x4mby(IF_RTCD(&cpi->rtcd), x);
  } else {
    vp8_encode_intra16x16mby(IF_RTCD(&cpi->rtcd), x);
  }

  if (mbmi->mode != I8X8_PRED) {
    vp8_encode_intra16x16mbuv(IF_RTCD(&cpi->rtcd), x);
  }

  if (output_enabled) {
    int segment_id = mbmi->segment_id;

    // Tokenize
    sum_intra_stats(cpi, x);
    vp8_tokenize_mb(cpi, &x->e_mbd, t, 0);

    if (cpi->common.txfm_mode == TX_MODE_SELECT &&
        !((cpi->common.mb_no_coeff_skip && mbmi->mb_skip_coeff) ||
          (segfeature_active(&x->e_mbd, segment_id, SEG_LVL_EOB) &&
           get_segdata(&x->e_mbd, segment_id, SEG_LVL_EOB) == 0))) {
      if (mbmi->mode != B_PRED && mbmi->mode != I8X8_PRED) {
        cpi->txfm_count[mbmi->txfm_size]++;
      } else if (mbmi->mode == I8X8_PRED) {
        cpi->txfm_count_8x8p[mbmi->txfm_size]++;
      }
    } else if (cpi->common.txfm_mode >= ALLOW_16X16 && mbmi->mode <= TM_PRED) {
      mbmi->txfm_size = TX_16X16;
    } else
    if (cpi->common.txfm_mode >= ALLOW_8X8 && mbmi->mode != B_PRED) {
      mbmi->txfm_size = TX_8X8;
    } else {
      mbmi->txfm_size = TX_4X4;
    }
  }
#if CONFIG_NEWBESTREFMV
  else
    vp8_tokenize_mb(cpi, &x->e_mbd, t, 1);
#endif
}
#ifdef SPEEDSTATS
extern int cnt_pm;
#endif

extern void vp8_fix_contexts(MACROBLOCKD *xd);

void vp8cx_encode_inter_macroblock (VP8_COMP *cpi, MACROBLOCK *x,
                                    TOKENEXTRA **t, int recon_yoffset,
                                    int recon_uvoffset, int output_enabled) {
  VP8_COMMON *cm = &cpi->common;
  MACROBLOCKD *const xd = &x->e_mbd;
  MB_MODE_INFO * mbmi = &xd->mode_info_context->mbmi;
  unsigned char *segment_id = &mbmi->segment_id;
  int seg_ref_active;
  unsigned char ref_pred_flag;

  x->skip = 0;
#if CONFIG_SUPERBLOCKS
  assert(!xd->mode_info_context->mbmi.encoded_as_sb);
#endif

#if CONFIG_SWITCHABLE_INTERP
  vp8_setup_interp_filters(xd, mbmi->interp_filter, cm);
#endif
  if (cpi->oxcf.tuning == VP8_TUNE_SSIM) {
    // Adjust the zbin based on this MB rate.
    adjust_act_zbin(cpi, x);
  }

  {
    // Experimental code. Special case for gf and arf zeromv modes.
    // Increase zbin size to suppress noise
    cpi->zbin_mode_boost = 0;
    if (cpi->zbin_mode_boost_enabled) {
      if (mbmi->ref_frame != INTRA_FRAME) {
        if (mbmi->mode == ZEROMV) {
          if (mbmi->ref_frame != LAST_FRAME)
            cpi->zbin_mode_boost = GF_ZEROMV_ZBIN_BOOST;
          else
            cpi->zbin_mode_boost = LF_ZEROMV_ZBIN_BOOST;
        } else if (mbmi->mode == SPLITMV)
          cpi->zbin_mode_boost = 0;
        else
          cpi->zbin_mode_boost = MV_ZBIN_BOOST;
      }
    }

    vp8_update_zbin_extra(cpi, x);
  }

  seg_ref_active = segfeature_active(xd, *segment_id, SEG_LVL_REF_FRAME);

  // SET VARIOUS PREDICTION FLAGS

  // Did the chosen reference frame match its predicted value.
  ref_pred_flag = ((mbmi->ref_frame == get_pred_ref(cm, xd)));
  set_pred_flag(xd, PRED_REF, ref_pred_flag);

  if (mbmi->ref_frame == INTRA_FRAME) {
    if (mbmi->mode == B_PRED) {
      vp8_intra_prediction_down_copy(xd);
      vp8_encode_intra16x16mbuv(IF_RTCD(&cpi->rtcd), x);
      vp8_encode_intra4x4mby(IF_RTCD(&cpi->rtcd), x);
    } else if (mbmi->mode == I8X8_PRED) {
      vp8_encode_intra8x8mby(IF_RTCD(&cpi->rtcd), x);
      vp8_encode_intra8x8mbuv(IF_RTCD(&cpi->rtcd), x);
    } else {
      vp8_encode_intra16x16mbuv(IF_RTCD(&cpi->rtcd), x);
      vp8_encode_intra16x16mby(IF_RTCD(&cpi->rtcd), x);
    }

    if (output_enabled)
      sum_intra_stats(cpi, x);
  } else {
    int ref_fb_idx;

    if (mbmi->ref_frame == LAST_FRAME)
      ref_fb_idx = cpi->common.lst_fb_idx;
    else if (mbmi->ref_frame == GOLDEN_FRAME)
      ref_fb_idx = cpi->common.gld_fb_idx;
    else
      ref_fb_idx = cpi->common.alt_fb_idx;

    xd->pre.y_buffer = cpi->common.yv12_fb[ref_fb_idx].y_buffer + recon_yoffset;
    xd->pre.u_buffer = cpi->common.yv12_fb[ref_fb_idx].u_buffer + recon_uvoffset;
    xd->pre.v_buffer = cpi->common.yv12_fb[ref_fb_idx].v_buffer + recon_uvoffset;

    if (mbmi->second_ref_frame) {
      int second_ref_fb_idx;

      if (mbmi->second_ref_frame == LAST_FRAME)
        second_ref_fb_idx = cpi->common.lst_fb_idx;
      else if (mbmi->second_ref_frame == GOLDEN_FRAME)
        second_ref_fb_idx = cpi->common.gld_fb_idx;
      else
        second_ref_fb_idx = cpi->common.alt_fb_idx;

      xd->second_pre.y_buffer = cpi->common.yv12_fb[second_ref_fb_idx].y_buffer +
                                recon_yoffset;
      xd->second_pre.u_buffer = cpi->common.yv12_fb[second_ref_fb_idx].u_buffer +
                                recon_uvoffset;
      xd->second_pre.v_buffer = cpi->common.yv12_fb[second_ref_fb_idx].v_buffer +
                                recon_uvoffset;
    }

    if (!x->skip) {
      vp8_encode_inter16x16(IF_RTCD(&cpi->rtcd), x);

      // Clear mb_skip_coeff if mb_no_coeff_skip is not set
      if (!cpi->common.mb_no_coeff_skip)
        mbmi->mb_skip_coeff = 0;

    } else {
      vp8_build_1st_inter16x16_predictors_mb(xd, xd->dst.y_buffer,
                                             xd->dst.u_buffer, xd->dst.v_buffer,
                                             xd->dst.y_stride,
                                             xd->dst.uv_stride);
    }
  }

  if (!x->skip) {
#ifdef ENC_DEBUG
    if (enc_debug) {
      int i;
      printf("Segment=%d [%d, %d]: %d %d:\n", mbmi->segment_id, mb_col_debug,
             mb_row_debug, xd->mb_to_left_edge, xd->mb_to_top_edge);
      for (i = 0; i < 400; i++) {
        printf("%3d ", xd->qcoeff[i]);
        if (i % 16 == 15) printf("\n");
      }
      printf("\n");
      printf("eobs = ");
      for (i = 0; i < 25; i++)
        printf("%d:%d ", i, xd->block[i].eob);
      printf("\n");
      fflush(stdout);
    }
#endif

    vp8_tokenize_mb(cpi, xd, t, !output_enabled);

#ifdef ENC_DEBUG
    if (enc_debug) {
      printf("Tokenized\n");
      fflush(stdout);
    }
#endif
  } else {
    int mb_skip_context =
      cpi->common.mb_no_coeff_skip ?
      (x->e_mbd.mode_info_context - 1)->mbmi.mb_skip_coeff +
      (x->e_mbd.mode_info_context - cpi->common.mode_info_stride)->mbmi.mb_skip_coeff :
      0;
    if (cpi->common.mb_no_coeff_skip) {
      mbmi->mb_skip_coeff = 1;
      if (output_enabled)
        cpi->skip_true_count[mb_skip_context]++;
      vp8_fix_contexts(xd);
    } else {
      vp8_stuff_mb(cpi, xd, t, !output_enabled);
      mbmi->mb_skip_coeff = 0;
      if (output_enabled)
        cpi->skip_false_count[mb_skip_context]++;
    }
  }

  if (output_enabled) {
    int segment_id = mbmi->segment_id;
    if (cpi->common.txfm_mode == TX_MODE_SELECT &&
        !((cpi->common.mb_no_coeff_skip && mbmi->mb_skip_coeff) ||
          (segfeature_active(&x->e_mbd, segment_id, SEG_LVL_EOB) &&
           get_segdata(&x->e_mbd, segment_id, SEG_LVL_EOB) == 0))) {
      if (mbmi->mode != B_PRED && mbmi->mode != I8X8_PRED &&
          mbmi->mode != SPLITMV) {
        cpi->txfm_count[mbmi->txfm_size]++;
      } else if (mbmi->mode == I8X8_PRED ||
                 (mbmi->mode == SPLITMV &&
                  mbmi->partitioning != PARTITIONING_4X4)) {
        cpi->txfm_count_8x8p[mbmi->txfm_size]++;
      }
    } else if (mbmi->mode != B_PRED && mbmi->mode != I8X8_PRED &&
        mbmi->mode != SPLITMV && cpi->common.txfm_mode >= ALLOW_16X16) {
      mbmi->txfm_size = TX_16X16;
    } else if (mbmi->mode != B_PRED &&
               !(mbmi->mode == SPLITMV &&
                 mbmi->partitioning == PARTITIONING_4X4) &&
               cpi->common.txfm_mode >= ALLOW_8X8) {
      mbmi->txfm_size = TX_8X8;
    } else {
      mbmi->txfm_size = TX_4X4;
    }
  }
}

#if CONFIG_SUPERBLOCKS
void vp8cx_encode_inter_superblock(VP8_COMP *cpi, MACROBLOCK *x, TOKENEXTRA **t,
                                   int recon_yoffset, int recon_uvoffset, int mb_col, int mb_row) {
  const int output_enabled = 1;
  VP8_COMMON *cm = &cpi->common;
  MACROBLOCKD *xd = &x->e_mbd;
  const uint8_t *src = x->src.y_buffer;
  uint8_t *dst = xd->dst.y_buffer;
  const uint8_t *usrc = x->src.u_buffer;
  uint8_t *udst = xd->dst.u_buffer;
  const uint8_t *vsrc = x->src.v_buffer;
  uint8_t *vdst = xd->dst.v_buffer;
  int src_y_stride = x->src.y_stride, dst_y_stride = xd->dst.y_stride;
  int src_uv_stride = x->src.uv_stride, dst_uv_stride = xd->dst.uv_stride;
  const VP8_ENCODER_RTCD *rtcd = IF_RTCD(&cpi->rtcd);
  unsigned int segment_id = xd->mode_info_context->mbmi.segment_id;
  int seg_ref_active;
  unsigned char ref_pred_flag;
  int n;
  TOKENEXTRA *tp[4];
  int skip[4];
  MODE_INFO *mi = x->e_mbd.mode_info_context;
  ENTROPY_CONTEXT_PLANES ta[4], tl[4];

  x->skip = 0;

  if (cpi->oxcf.tuning == VP8_TUNE_SSIM) {
    // Adjust the zbin based on this MB rate.
    adjust_act_zbin(cpi, x);
  }

  {
    // Experimental code. Special case for gf and arf zeromv modes.
    // Increase zbin size to suppress noise
    cpi->zbin_mode_boost = 0;
    if (cpi->zbin_mode_boost_enabled) {
      if (xd->mode_info_context->mbmi.ref_frame != INTRA_FRAME) {
        if (xd->mode_info_context->mbmi.mode == ZEROMV) {
          if (xd->mode_info_context->mbmi.ref_frame != LAST_FRAME)
            cpi->zbin_mode_boost = GF_ZEROMV_ZBIN_BOOST;
          else
            cpi->zbin_mode_boost = LF_ZEROMV_ZBIN_BOOST;
        } else if (xd->mode_info_context->mbmi.mode == SPLITMV)
          cpi->zbin_mode_boost = 0;
        else
          cpi->zbin_mode_boost = MV_ZBIN_BOOST;
      }
    }

    vp8_update_zbin_extra(cpi, x);
  }

  seg_ref_active = segfeature_active(xd, segment_id, SEG_LVL_REF_FRAME);

  // SET VARIOUS PREDICTION FLAGS

  // Did the chosen reference frame match its predicted value.
  ref_pred_flag = ((xd->mode_info_context->mbmi.ref_frame ==
                    get_pred_ref(cm, xd)));
  set_pred_flag(xd, PRED_REF, ref_pred_flag);

  if (xd->mode_info_context->mbmi.ref_frame == INTRA_FRAME) {
    vp8_build_intra_predictors_sby_s(&x->e_mbd);
    vp8_build_intra_predictors_sbuv_s(&x->e_mbd);
  } else {
    int ref_fb_idx;

    if (xd->mode_info_context->mbmi.ref_frame == LAST_FRAME)
      ref_fb_idx = cpi->common.lst_fb_idx;
    else if (xd->mode_info_context->mbmi.ref_frame == GOLDEN_FRAME)
      ref_fb_idx = cpi->common.gld_fb_idx;
    else
      ref_fb_idx = cpi->common.alt_fb_idx;

    xd->pre.y_buffer = cpi->common.yv12_fb[ref_fb_idx].y_buffer + recon_yoffset;
    xd->pre.u_buffer = cpi->common.yv12_fb[ref_fb_idx].u_buffer + recon_uvoffset;
    xd->pre.v_buffer = cpi->common.yv12_fb[ref_fb_idx].v_buffer + recon_uvoffset;

    if (xd->mode_info_context->mbmi.second_ref_frame) {
      int second_ref_fb_idx;

      if (xd->mode_info_context->mbmi.second_ref_frame == LAST_FRAME)
        second_ref_fb_idx = cpi->common.lst_fb_idx;
      else if (xd->mode_info_context->mbmi.second_ref_frame == GOLDEN_FRAME)
        second_ref_fb_idx = cpi->common.gld_fb_idx;
      else
        second_ref_fb_idx = cpi->common.alt_fb_idx;

      xd->second_pre.y_buffer = cpi->common.yv12_fb[second_ref_fb_idx].y_buffer +
                                    recon_yoffset;
      xd->second_pre.u_buffer = cpi->common.yv12_fb[second_ref_fb_idx].u_buffer +
                                    recon_uvoffset;
      xd->second_pre.v_buffer = cpi->common.yv12_fb[second_ref_fb_idx].v_buffer +
                                    recon_uvoffset;
    }

    vp8_build_inter32x32_predictors_sb(xd, xd->dst.y_buffer,
                                       xd->dst.u_buffer, xd->dst.v_buffer,
                                       xd->dst.y_stride, xd->dst.uv_stride);
  }

  assert(x->e_mbd.mode_info_context->mbmi.txfm_size == TX_8X8);
  for (n = 0; n < 4; n++)
  {
    int x_idx = n & 1, y_idx = n >> 1;

    vp8_subtract_mby_s_c(x->src_diff,
                         src + x_idx * 16 + y_idx * 16 * src_y_stride,
                         src_y_stride,
                         dst + x_idx * 16 + y_idx * 16 * dst_y_stride,
                         dst_y_stride);
    vp8_subtract_mbuv_s_c(x->src_diff,
                          usrc + x_idx * 8 + y_idx * 8 * src_uv_stride,
                          vsrc + x_idx * 8 + y_idx * 8 * src_uv_stride,
                          src_uv_stride,
                          udst + x_idx * 8 + y_idx * 8 * dst_uv_stride,
                          vdst + x_idx * 8 + y_idx * 8 * dst_uv_stride,
                          dst_uv_stride);
    vp8_transform_mb_8x8(x);
    vp8_quantize_mb_8x8(x);
    if (x->optimize) {
      vp8_optimize_mby_8x8(x, rtcd);
      vp8_optimize_mbuv_8x8(x, rtcd);
    }
    vp8_inverse_transform_mb_8x8(IF_RTCD(&rtcd->common->idct), &x->e_mbd);
    vp8_recon_mby_s_c( &x->e_mbd,
                      dst + x_idx * 16 + y_idx * 16 * dst_y_stride);
    vp8_recon_mbuv_s_c(&x->e_mbd,
                       udst + x_idx * 8 + y_idx * 8 * dst_uv_stride,
                       vdst + x_idx * 8 + y_idx * 8 * dst_uv_stride);

    if (!x->skip) {
      if (output_enabled) {
        xd->left_context = cm->left_context + (n >> 1);
        xd->above_context = cm->above_context + mb_col + (n >> 1);
        memcpy(&ta[n], xd->above_context, sizeof(ta[n]));
        memcpy(&tl[n], xd->left_context, sizeof(tl[n]));
        tp[n] = *t;
        xd->mode_info_context = mi + x_idx + y_idx * cm->mode_info_stride;
        vp8_tokenize_mb(cpi, &x->e_mbd, t, 0);
        skip[n] = xd->mode_info_context->mbmi.mb_skip_coeff;
      }
    } else {
      int mb_skip_context =
        cpi->common.mb_no_coeff_skip ?
          (x->e_mbd.mode_info_context - 1)->mbmi.mb_skip_coeff +
            (x->e_mbd.mode_info_context - cpi->common.mode_info_stride)->mbmi.mb_skip_coeff :
          0;
      if (cpi->common.mb_no_coeff_skip) {
        skip[n] = xd->mode_info_context->mbmi.mb_skip_coeff = 1;
        xd->left_context = cm->left_context + (n >> 1);
        xd->above_context = cm->above_context + mb_col + (n >> 1);
        memcpy(&ta[n], xd->above_context, sizeof(ta[n]));
        memcpy(&tl[n], xd->left_context, sizeof(tl[n]));
        tp[n] = *t;
        cpi->skip_true_count[mb_skip_context]++;
        vp8_fix_contexts(xd);
      } else {
        vp8_stuff_mb(cpi, xd, t, 0);
        xd->mode_info_context->mbmi.mb_skip_coeff = 0;
        cpi->skip_false_count[mb_skip_context]++;
      }
    }
  }

  xd->mode_info_context = mi;
  update_sb_skip_coeff_state(cpi, x, ta, tl, tp, t, skip);
}
#endif
