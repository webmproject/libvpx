/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */


#include "./vpx_config.h"
#include "vp9/encoder/vp9_encodeframe.h"
#include "vp9/encoder/vp9_encodemb.h"
#include "vp9/encoder/vp9_encodemv.h"
#include "vp9/common/vp9_common.h"
#include "vp9/encoder/vp9_onyx_int.h"
#include "vp9/common/vp9_extend.h"
#include "vp9/common/vp9_entropy.h"
#include "vp9/common/vp9_entropymode.h"
#include "vp9/common/vp9_quant_common.h"
#include "vp9/encoder/vp9_segmentation.h"
#include "vp9/common/vp9_setupintrarecon.h"
#include "vp9/encoder/vp9_encodeintra.h"
#include "vp9/common/vp9_reconinter.h"
#include "vp9/common/vp9_invtrans.h"
#include "vp9/encoder/vp9_rdopt.h"
#include "vp9/common/vp9_findnearmv.h"
#include "vp9/common/vp9_reconintra.h"
#include "vp9/common/vp9_seg_common.h"
#include "vp9/common/vp9_tile_common.h"
#include "vp9/encoder/vp9_tokenize.h"
#include "./vp9_rtcd.h"
#include <stdio.h>
#include <math.h>
#include <limits.h>
#include "vpx_ports/vpx_timer.h"
#include "vp9/common/vp9_pred_common.h"
#include "vp9/common/vp9_mvref_common.h"

#define DBG_PRNT_SEGMAP 0

// #define ENC_DEBUG
#ifdef ENC_DEBUG
int enc_debug = 0;
#endif

void vp9_select_interp_filter_type(VP9_COMP *cpi);

#if !CONFIG_SB8X8
static void encode_macroblock(VP9_COMP *cpi, TOKENEXTRA **t,
                              int output_enabled, int mi_row, int mi_col);
#endif

static void encode_superblock(VP9_COMP *cpi, TOKENEXTRA **t,
                              int output_enabled, int mi_row, int mi_col,
                              BLOCK_SIZE_TYPE bsize);

static void adjust_act_zbin(VP9_COMP *cpi, MACROBLOCK *x);

#ifdef MODE_STATS
unsigned int inter_y_modes[MB_MODE_COUNT];
unsigned int inter_uv_modes[VP9_UV_MODES];
unsigned int inter_b_modes[B_MODE_COUNT];
unsigned int y_modes[VP9_YMODES];
unsigned int i8x8_modes[VP9_I8X8_MODES];
unsigned int uv_modes[VP9_UV_MODES];
unsigned int uv_modes_y[VP9_YMODES][VP9_UV_MODES];
unsigned int b_modes[B_MODE_COUNT];
#endif


/* activity_avg must be positive, or flat regions could get a zero weight
 *  (infinite lambda), which confounds analysis.
 * This also avoids the need for divide by zero checks in
 *  vp9_activity_masking().
 */
#define VP9_ACTIVITY_AVG_MIN (64)

/* This is used as a reference when computing the source variance for the
 *  purposes of activity masking.
 * Eventually this should be replaced by custom no-reference routines,
 *  which will be faster.
 */
static const uint8_t VP9_VAR_OFFS[16] = {
  128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128
};


// Original activity measure from Tim T's code.
static unsigned int tt_activity_measure(VP9_COMP *cpi, MACROBLOCK *x) {
  unsigned int act;
  unsigned int sse;
  /* TODO: This could also be done over smaller areas (8x8), but that would
   *  require extensive changes elsewhere, as lambda is assumed to be fixed
   *  over an entire MB in most of the code.
   * Another option is to compute four 8x8 variances, and pick a single
   *  lambda using a non-linear combination (e.g., the smallest, or second
   *  smallest, etc.).
   */
  act = vp9_variance16x16(x->plane[0].src.buf, x->plane[0].src.stride,
                          VP9_VAR_OFFS, 0, &sse);
  act <<= 4;

  /* If the region is flat, lower the activity some more. */
  if (act < 8 << 12)
    act = act < 5 << 12 ? act : 5 << 12;

  return act;
}

// Stub for alternative experimental activity measures.
static unsigned int alt_activity_measure(VP9_COMP *cpi,
                                         MACROBLOCK *x, int use_dc_pred) {
  return vp9_encode_intra(cpi, x, use_dc_pred);
}


// Measure the activity of the current macroblock
// What we measure here is TBD so abstracted to this function
#define ALT_ACT_MEASURE 1
static unsigned int mb_activity_measure(VP9_COMP *cpi, MACROBLOCK *x,
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

  if (mb_activity < VP9_ACTIVITY_AVG_MIN)
    mb_activity = VP9_ACTIVITY_AVG_MIN;

  return mb_activity;
}

// Calculate an "average" mb activity value for the frame
#define ACT_MEDIAN 0
static void calc_av_activity(VP9_COMP *cpi, int64_t activity_sum) {
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

  if (cpi->activity_avg < VP9_ACTIVITY_AVG_MIN)
    cpi->activity_avg = VP9_ACTIVITY_AVG_MIN;

  // Experimental code: return fixed value normalized for several clips
  if (ALT_ACT_MEASURE)
    cpi->activity_avg = 100000;
}

#define USE_ACT_INDEX   0
#define OUTPUT_NORM_ACT_STATS   0

#if USE_ACT_INDEX
// Calculate an activity index for each mb
static void calc_activity_index(VP9_COMP *cpi, MACROBLOCK *x) {
  VP9_COMMON *const cm = &cpi->common;
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
static void build_activity_map(VP9_COMP *cpi) {
  MACROBLOCK *const x = &cpi->mb;
  MACROBLOCKD *xd = &x->e_mbd;
  VP9_COMMON *const cm = &cpi->common;

#if ALT_ACT_MEASURE
  YV12_BUFFER_CONFIG *new_yv12 = &cm->yv12_fb[cm->new_fb_idx];
  int recon_yoffset;
  int recon_y_stride = new_yv12->y_stride;
#endif

  int mb_row, mb_col;
  unsigned int mb_activity;
  int64_t activity_sum = 0;

  x->mb_activity_ptr = cpi->mb_activity_map;

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
      xd->plane[0].dst.buf = new_yv12->y_buffer + recon_yoffset;
      xd->left_available = (mb_col != 0);
      recon_yoffset += 16;
#endif

      // measure activity
      mb_activity = mb_activity_measure(cpi, x, mb_row, mb_col);

      // Keep frame sum
      activity_sum += mb_activity;

      // Store MB level activity details.
      *x->mb_activity_ptr = mb_activity;

      // Increment activity map pointer
      x->mb_activity_ptr++;

      // adjust to the next column of source macroblocks
      x->plane[0].src.buf += 16;
    }


    // adjust to the next row of mbs
    x->plane[0].src.buf += 16 * x->plane[0].src.stride - 16 * cm->mb_cols;

#if ALT_ACT_MEASURE
    // extend the recon for intra prediction
    vp9_extend_mb_row(new_yv12, xd->plane[0].dst.buf + 16,
                      xd->plane[1].dst.buf + 8, xd->plane[2].dst.buf + 8);
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
void vp9_activity_masking(VP9_COMP *cpi, MACROBLOCK *x) {
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

static void update_state(VP9_COMP *cpi,
                         PICK_MODE_CONTEXT *ctx,
                         BLOCK_SIZE_TYPE bsize,
                         int output_enabled) {
  int i, x_idx, y;
  VP9_COMMON *const cm = &cpi->common;
  MACROBLOCK *const x = &cpi->mb;
  MACROBLOCKD *const xd = &x->e_mbd;
  MODE_INFO *mi = &ctx->mic;
  MB_MODE_INFO *const mbmi = &xd->mode_info_context->mbmi;
  int mb_mode = mi->mbmi.mode;
  int mb_mode_index = ctx->best_mode_index;
  const int mis = cpi->common.mode_info_stride;
  const int bh = 1 << mi_height_log2(bsize), bw = 1 << mi_width_log2(bsize);

#if CONFIG_DEBUG
  assert(mb_mode < MB_MODE_COUNT);
  assert(mb_mode_index < MAX_MODES);
  assert(mi->mbmi.ref_frame < MAX_REF_FRAMES);
#endif
  assert(mi->mbmi.sb_type == bsize);

  // Restore the coding context of the MB to that that was in place
  // when the mode was picked for it
  for (y = 0; y < bh; y++) {
    for (x_idx = 0; x_idx < bw; x_idx++) {
      if ((xd->mb_to_right_edge >> (3 + LOG2_MI_SIZE)) + bw > x_idx &&
          (xd->mb_to_bottom_edge >> (3 + LOG2_MI_SIZE)) + bh > y) {
        MODE_INFO *mi_addr = xd->mode_info_context + x_idx + y * mis;

        vpx_memcpy(mi_addr, mi, sizeof(MODE_INFO));
      }
    }
  }
  if (bsize < BLOCK_SIZE_SB32X32) {
    if (bsize < BLOCK_SIZE_MB16X16)
      ctx->txfm_rd_diff[ALLOW_16X16] = ctx->txfm_rd_diff[ALLOW_8X8];
    ctx->txfm_rd_diff[ALLOW_32X32] = ctx->txfm_rd_diff[ALLOW_16X16];
  }

  if (mb_mode == SPLITMV) {
    vpx_memcpy(x->partition_info, &ctx->partition_info,
               sizeof(PARTITION_INFO));

    mbmi->mv[0].as_int =
        x->partition_info->bmi[15 >> (CONFIG_SB8X8 * 2)].mv.as_int;
    mbmi->mv[1].as_int =
        x->partition_info->bmi[15 >> (CONFIG_SB8X8 * 2)].second_mv.as_int;
  }

  x->skip = ctx->skip;
  if (!output_enabled)
    return;

  {
    int segment_id = mbmi->segment_id, ref_pred_flag;
    if (!vp9_segfeature_active(xd, segment_id, SEG_LVL_SKIP)) {
      for (i = 0; i < NB_TXFM_MODES; i++) {
        cpi->rd_tx_select_diff[i] += ctx->txfm_rd_diff[i];
      }
    }

    // Did the chosen reference frame match its predicted value.
    ref_pred_flag = ((xd->mode_info_context->mbmi.ref_frame ==
                      vp9_get_pred_ref(cm, xd)));
    vp9_set_pred_flag(xd, PRED_REF, ref_pred_flag);
    if (!xd->segmentation_enabled ||
        !vp9_segfeature_active(xd, segment_id, SEG_LVL_REF_FRAME) ||
        vp9_check_segref(xd, segment_id, INTRA_FRAME)  +
        vp9_check_segref(xd, segment_id, LAST_FRAME)   +
        vp9_check_segref(xd, segment_id, GOLDEN_FRAME) +
        vp9_check_segref(xd, segment_id, ALTREF_FRAME) > 1) {
      // Get the prediction context and status
      int pred_context = vp9_get_pred_context(cm, xd, PRED_REF);

      // Count prediction success
      cpi->ref_pred_count[pred_context][ref_pred_flag]++;
    }
  }

  if (cpi->common.frame_type == KEY_FRAME) {
    // Restore the coding modes to that held in the coding context
    // if (mb_mode == I4X4_PRED)
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
#if !CONFIG_SB8X8
      THR_I8X8_PRED /*I8X8_PRED*/,
#endif
      THR_B_PRED /*I4X4_PRED*/,
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
    if (mbmi->mode == SPLITMV || mbmi->mode == NEWMV) {
      int_mv best_mv, best_second_mv;
      MV_REFERENCE_FRAME rf = mbmi->ref_frame;
      best_mv.as_int = ctx->best_ref_mv.as_int;
      best_second_mv.as_int = ctx->second_best_ref_mv.as_int;
      if (mbmi->mode == NEWMV) {
        best_mv.as_int = mbmi->ref_mvs[rf][0].as_int;
        best_second_mv.as_int = mbmi->ref_mvs[mbmi->second_ref_frame][0].as_int;
      }
      mbmi->best_mv.as_int = best_mv.as_int;
      mbmi->best_second_mv.as_int = best_second_mv.as_int;
      vp9_update_nmv_count(cpi, x, &best_mv, &best_second_mv);
#if CONFIG_SB8X8
      xd->mode_info_context[1].mbmi =
      xd->mode_info_context[mis].mbmi =
      xd->mode_info_context[1 + mis].mbmi = *mbmi;
#endif
    }
#if CONFIG_COMP_INTERINTRA_PRED
    if (mbmi->mode >= NEARESTMV && mbmi->mode < SPLITMV &&
        mbmi->second_ref_frame <= INTRA_FRAME) {
      if (mbmi->second_ref_frame == INTRA_FRAME) {
        ++cpi->interintra_count[1];
        ++cpi->ymode_count[mbmi->interintra_mode];
#if SEPARATE_INTERINTRA_UV
        ++cpi->y_uv_mode_count[mbmi->interintra_mode][mbmi->interintra_uv_mode];
#endif
      } else {
        ++cpi->interintra_count[0];
      }
    }
#endif
    if (cpi->common.mcomp_filter_type == SWITCHABLE &&
        is_inter_mode(mbmi->mode)) {
      ++cpi->switchable_interp_count
          [vp9_get_pred_context(&cpi->common, xd, PRED_SWITCHABLE_INTERP)]
          [vp9_switchable_interp_map[mbmi->interp_filter]];
    }

    cpi->rd_comp_pred_diff[SINGLE_PREDICTION_ONLY] += ctx->single_pred_diff;
    cpi->rd_comp_pred_diff[COMP_PREDICTION_ONLY]   += ctx->comp_pred_diff;
    cpi->rd_comp_pred_diff[HYBRID_PREDICTION]      += ctx->hybrid_pred_diff;
  }
}

static unsigned find_seg_id(uint8_t *buf, BLOCK_SIZE_TYPE bsize,
                            int start_y, int height, int start_x, int width) {
  const int bw = 1 << mi_width_log2(bsize), bh = 1 << mi_height_log2(bsize);
  const int end_x = MIN(start_x + bw, width);
  const int end_y = MIN(start_y + bh, height);
  int x, y;
  unsigned seg_id = -1;

  buf += width * start_y;
  for (y = start_y; y < end_y; y++, buf += width) {
    for (x = start_x; x < end_x; x++) {
      seg_id = MIN(seg_id, buf[x]);
    }
  }

  return seg_id;
}

void vp9_setup_src_planes(MACROBLOCK *x,
                          const YV12_BUFFER_CONFIG *src,
                          int mb_row, int mb_col) {
  setup_pred_plane(&x->plane[0].src,
                   src->y_buffer, src->y_stride,
                   mb_row, mb_col, NULL,
                   x->e_mbd.plane[0].subsampling_x,
                   x->e_mbd.plane[0].subsampling_y);
  setup_pred_plane(&x->plane[1].src,
                   src->u_buffer, src->uv_stride,
                   mb_row, mb_col, NULL,
                   x->e_mbd.plane[1].subsampling_x,
                   x->e_mbd.plane[1].subsampling_y);
  setup_pred_plane(&x->plane[2].src,
                   src->v_buffer, src->uv_stride,
                   mb_row, mb_col, NULL,
                   x->e_mbd.plane[2].subsampling_x,
                   x->e_mbd.plane[2].subsampling_y);
}

static INLINE void set_partition_seg_context(VP9_COMP *cpi,
                                             int mi_row, int mi_col) {
  VP9_COMMON *const cm = &cpi->common;
  MACROBLOCKD *const xd = &cpi->mb.e_mbd;

  xd->above_seg_context = cm->above_seg_context + (mi_col >> CONFIG_SB8X8);
  xd->left_seg_context  = cm->left_seg_context + ((mi_row >> CONFIG_SB8X8) & 3);
}

static void set_offsets(VP9_COMP *cpi,
                        int mi_row, int mi_col, BLOCK_SIZE_TYPE bsize) {
  MACROBLOCK *const x = &cpi->mb;
  VP9_COMMON *const cm = &cpi->common;
  MACROBLOCKD *const xd = &x->e_mbd;
  MB_MODE_INFO *mbmi;
  const int dst_fb_idx = cm->new_fb_idx;
  const int idx_str = xd->mode_info_stride * mi_row + mi_col;
  const int bw = 1 << mi_width_log2(bsize), bh = 1 << mi_height_log2(bsize);
  const int mb_row = mi_row >> CONFIG_SB8X8;
  const int mb_col = mi_col >> CONFIG_SB8X8;
  const int idx_map = mb_row * cm->mb_cols + mb_col;
  int i;

  // entropy context structures
  for (i = 0; i < MAX_MB_PLANE; i++) {
    xd->plane[i].above_context = cm->above_context[i] +
        (mi_col * 4 >> (CONFIG_SB8X8 + xd->plane[i].subsampling_x));
    xd->plane[i].left_context = cm->left_context[i] +
        (((mi_row * 4 >> CONFIG_SB8X8) & 15) >> xd->plane[i].subsampling_y);
  }

  // partition contexts
  set_partition_seg_context(cpi, mi_row, mi_col);

  // Activity map pointer
  x->mb_activity_ptr = &cpi->mb_activity_map[idx_map];
  x->active_ptr = cpi->active_map + idx_map;

  /* pointers to mode info contexts */
  x->partition_info          = x->pi + idx_str;
  xd->mode_info_context      = cm->mi + idx_str;
  mbmi = &xd->mode_info_context->mbmi;
  xd->prev_mode_info_context = cm->prev_mi + idx_str;

  // Set up destination pointers
  setup_dst_planes(xd, &cm->yv12_fb[dst_fb_idx], mi_row, mi_col);

  /* Set up limit values for MV components to prevent them from
   * extending beyond the UMV borders assuming 16x16 block size */
  x->mv_row_min = -((mi_row * MI_SIZE) + VP9BORDERINPIXELS - VP9_INTERP_EXTEND);
  x->mv_col_min = -((mi_col * MI_SIZE) + VP9BORDERINPIXELS - VP9_INTERP_EXTEND);
  x->mv_row_max = ((cm->mi_rows - mi_row) * MI_SIZE +
                   (VP9BORDERINPIXELS - MI_SIZE * bh - VP9_INTERP_EXTEND));
  x->mv_col_max = ((cm->mi_cols - mi_col) * MI_SIZE +
                   (VP9BORDERINPIXELS - MI_SIZE * bw - VP9_INTERP_EXTEND));

  // Set up distance of MB to edge of frame in 1/8th pel units
  assert(!(mi_col & (bw - 1)) && !(mi_row & (bh - 1)));
  set_mi_row_col(cm, xd, mi_row, bh, mi_col, bw);

  /* set up source buffers */
  vp9_setup_src_planes(x, cpi->Source, mi_row, mi_col);

  /* R/D setup */
  x->rddiv = cpi->RDDIV;
  x->rdmult = cpi->RDMULT;

  /* segment ID */
  if (xd->segmentation_enabled) {
    uint8_t *map = xd->update_mb_segmentation_map ? cpi->segmentation_map
                                                  : cm->last_frame_seg_map;
    mbmi->segment_id = find_seg_id(map, bsize, mi_row,
                                   cm->mi_rows, mi_col, cm->mi_cols);

    assert(mbmi->segment_id <= (MAX_MB_SEGMENTS-1));
    vp9_mb_init_quantizer(cpi, x);

    if (xd->segmentation_enabled && cpi->seg0_cnt > 0 &&
        !vp9_segfeature_active(xd, 0, SEG_LVL_REF_FRAME) &&
        vp9_segfeature_active(xd, 1, SEG_LVL_REF_FRAME) &&
        vp9_check_segref(xd, 1, INTRA_FRAME)  +
        vp9_check_segref(xd, 1, LAST_FRAME)   +
        vp9_check_segref(xd, 1, GOLDEN_FRAME) +
        vp9_check_segref(xd, 1, ALTREF_FRAME) == 1) {
      cpi->seg0_progress = (cpi->seg0_idx << 16) / cpi->seg0_cnt;
    } else {
      const int y = mb_row & ~3;
      const int x = mb_col & ~3;
      const int p16 = ((mb_row & 1) << 1) +  (mb_col & 1);
      const int p32 = ((mb_row & 2) << 2) + ((mb_col & 2) << 1);
      const int tile_progress =
          cm->cur_tile_mi_col_start * cm->mb_rows >> CONFIG_SB8X8;
      const int mb_cols =
          (cm->cur_tile_mi_col_end - cm->cur_tile_mi_col_start) >> CONFIG_SB8X8;

      cpi->seg0_progress =
          ((y * mb_cols + x * 4 + p32 + p16 + tile_progress) << 16) / cm->MBs;
    }
  } else {
    mbmi->segment_id = 0;
  }
}

#if !CONFIG_SB8X8
static int pick_mb_mode(VP9_COMP *cpi,
                        int mi_row,
                        int mi_col,
                        TOKENEXTRA **tp,
                        int *totalrate,
                        int *totaldist) {
  VP9_COMMON *const cm = &cpi->common;
  MACROBLOCK *const x = &cpi->mb;
  MACROBLOCKD *const xd = &x->e_mbd;
  int splitmodes_used = 0;
  MB_MODE_INFO *mbmi;

  set_offsets(cpi, mi_row, mi_col, BLOCK_SIZE_MB16X16);

  if (cpi->oxcf.tuning == VP8_TUNE_SSIM)
    vp9_activity_masking(cpi, x);

  mbmi = &xd->mode_info_context->mbmi;
  mbmi->sb_type = BLOCK_SIZE_MB16X16;

  // Find best coding mode & reconstruct the MB so it is available
  // as a predictor for MBs that follow in the SB
  if (cm->frame_type == KEY_FRAME) {
    vp9_rd_pick_intra_mode(cpi, x, totalrate, totaldist);

    // Save the coding context
    vpx_memcpy(&x->mb_context[xd->sb_index][xd->mb_index].mic,
               xd->mode_info_context, sizeof(MODE_INFO));
  } else {
    vp9_pick_mode_inter_macroblock(cpi, x, mi_row, mi_col,
                                   totalrate, totaldist);
    splitmodes_used += (mbmi->mode == SPLITMV);

    if (cpi->mb.e_mbd.segmentation_enabled && mbmi->segment_id == 0) {
      cpi->seg0_idx++;
    }
  }

  return splitmodes_used;
}
#endif

static void pick_sb_modes(VP9_COMP *cpi, int mi_row, int mi_col,
                          TOKENEXTRA **tp, int *totalrate, int *totaldist,
                          BLOCK_SIZE_TYPE bsize, PICK_MODE_CONTEXT *ctx) {
  VP9_COMMON *const cm = &cpi->common;
  MACROBLOCK *const x = &cpi->mb;
  MACROBLOCKD *const xd = &x->e_mbd;

  set_offsets(cpi, mi_row, mi_col, bsize);
  xd->mode_info_context->mbmi.sb_type = bsize;
  if (cpi->oxcf.tuning == VP8_TUNE_SSIM)
    vp9_activity_masking(cpi, x);

  /* Find best coding mode & reconstruct the MB so it is available
   * as a predictor for MBs that follow in the SB */
  if (cm->frame_type == KEY_FRAME) {
    vp9_rd_pick_intra_mode_sb(cpi, x, totalrate, totaldist, bsize, ctx);
  } else {
    vp9_rd_pick_inter_mode_sb(cpi, x, mi_row, mi_col, totalrate, totaldist,
                              bsize, ctx);
  }
}

static void update_stats(VP9_COMP *cpi, int mi_row, int mi_col) {
  VP9_COMMON *const cm = &cpi->common;
  MACROBLOCK *const x = &cpi->mb;
  MACROBLOCKD *const xd = &x->e_mbd;
  MODE_INFO *mi = xd->mode_info_context;
  MB_MODE_INFO *const mbmi = &mi->mbmi;

  if (cm->frame_type == KEY_FRAME) {
#ifdef MODE_STATS
    y_modes[mbmi->mode]++;
#endif
  } else {
    int segment_id, seg_ref_active;

    if (mbmi->ref_frame) {
      int pred_context = vp9_get_pred_context(cm, xd, PRED_COMP);

      if (mbmi->second_ref_frame <= INTRA_FRAME)
        cpi->single_pred_count[pred_context]++;
      else
        cpi->comp_pred_count[pred_context]++;
    }

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
    segment_id = mbmi->segment_id;
    seg_ref_active = vp9_segfeature_active(xd, segment_id,
                                           SEG_LVL_REF_FRAME);
    if (!seg_ref_active ||
        ((vp9_check_segref(xd, segment_id, INTRA_FRAME) +
          vp9_check_segref(xd, segment_id, LAST_FRAME) +
          vp9_check_segref(xd, segment_id, GOLDEN_FRAME) +
          vp9_check_segref(xd, segment_id, ALTREF_FRAME)) > 1)) {
      cpi->count_mb_ref_frame_usage[mbmi->ref_frame]++;
    }
    // Count of last ref frame 0,0 usage
    if ((mbmi->mode == ZEROMV) && (mbmi->ref_frame == LAST_FRAME))
      cpi->inter_zz_count++;
  }
}

static void set_block_index(MACROBLOCKD *xd, int idx,
                            BLOCK_SIZE_TYPE bsize) {
  if (bsize >= BLOCK_SIZE_SB32X32) {
    xd->sb_index = idx;
#if CONFIG_SB8X8
  } else if (bsize >= BLOCK_SIZE_MB16X16) {
    xd->mb_index = idx;
  } else {
    xd->b_index = idx;
#else
  } else {
    xd->mb_index = idx;
#endif
  }
}

static PICK_MODE_CONTEXT *get_block_context(MACROBLOCK *x,
                                            BLOCK_SIZE_TYPE bsize) {
  MACROBLOCKD *const xd = &x->e_mbd;

  switch (bsize) {
    case BLOCK_SIZE_SB64X64:
      return &x->sb64_context;
    case BLOCK_SIZE_SB64X32:
      return &x->sb64x32_context[xd->sb_index];
    case BLOCK_SIZE_SB32X64:
      return &x->sb32x64_context[xd->sb_index];
    case BLOCK_SIZE_SB32X32:
      return &x->sb32_context[xd->sb_index];
    case BLOCK_SIZE_SB32X16:
      return &x->sb32x16_context[xd->sb_index][xd->mb_index];
    case BLOCK_SIZE_SB16X32:
      return &x->sb16x32_context[xd->sb_index][xd->mb_index];
    case BLOCK_SIZE_MB16X16:
      return &x->mb_context[xd->sb_index][xd->mb_index];
#if CONFIG_SB8X8
    case BLOCK_SIZE_SB16X8:
      return &x->sb16x8_context[xd->sb_index][xd->mb_index][xd->b_index];
    case BLOCK_SIZE_SB8X16:
      return &x->sb8x16_context[xd->sb_index][xd->mb_index][xd->b_index];
    case BLOCK_SIZE_SB8X8:
      return &x->sb8_context[xd->sb_index][xd->mb_index][xd->b_index];
#endif
    default:
      assert(0);
      return NULL;
  }
}

static void encode_b(VP9_COMP *cpi, TOKENEXTRA **tp,
                     int mi_row, int mi_col, int output_enabled,
                     BLOCK_SIZE_TYPE bsize, int sub_index) {
  VP9_COMMON *const cm = &cpi->common;
  MACROBLOCK *const x = &cpi->mb;
  MACROBLOCKD *const xd = &x->e_mbd;

  if (mi_row >= cm->mi_rows || mi_col >= cm->mi_cols)
    return;

  if (sub_index != -1)
    set_block_index(xd, sub_index, bsize);
  set_offsets(cpi, mi_row, mi_col, bsize);
  update_state(cpi, get_block_context(x, bsize), bsize, output_enabled);
#if !CONFIG_SB8X8
  if (bsize == BLOCK_SIZE_MB16X16) {
    if (cpi->oxcf.tuning == VP8_TUNE_SSIM)
      vp9_activity_masking(cpi, x);

    encode_macroblock(cpi, tp, output_enabled, mi_row, mi_col);
  } else
#endif
  {
    encode_superblock(cpi, tp, output_enabled, mi_row, mi_col, bsize);
  }

  if (output_enabled) {
    update_stats(cpi, mi_row, mi_col);

    (*tp)->token = EOSB_TOKEN;
    (*tp)++;
  }
}

static void encode_sb(VP9_COMP *cpi, TOKENEXTRA **tp,
                      int mi_row, int mi_col, int output_enabled,
                      BLOCK_SIZE_TYPE level,
                      BLOCK_SIZE_TYPE c1, BLOCK_SIZE_TYPE c2[4]
#if CONFIG_SB8X8
                      , BLOCK_SIZE_TYPE c3[4][4]
#endif
                      ) {
  VP9_COMMON *const cm = &cpi->common;
  MACROBLOCK *const x = &cpi->mb;
  MACROBLOCKD *const xd = &x->e_mbd;
  const int bsl = mi_width_log2(level), bs = 1 << (bsl - 1);
  const int bwl = mi_width_log2(c1), bhl = mi_height_log2(c1);
  int UNINITIALIZED_IS_SAFE(pl);

  if (mi_row >= cm->mi_rows || mi_col >= cm->mi_cols)
    return;

#if CONFIG_SB8X8
  if (level > BLOCK_SIZE_SB8X8) {
#endif
    set_partition_seg_context(cpi, mi_row, mi_col);
    pl = partition_plane_context(xd, level);
#if CONFIG_SB8X8
  }
#endif

  if (bsl == bwl && bsl == bhl) {
    if (output_enabled &&
#if CONFIG_SB8X8
        level > BLOCK_SIZE_SB8X8
#else
        level > BLOCK_SIZE_MB16X16
#endif
        )
      cpi->partition_count[pl][PARTITION_NONE]++;
    encode_b(cpi, tp, mi_row, mi_col, output_enabled, c1, -1);
  } else if (bsl == bhl && bsl > bwl) {
    if (output_enabled)
      cpi->partition_count[pl][PARTITION_VERT]++;
    encode_b(cpi, tp, mi_row, mi_col,      output_enabled, c1, 0);
    encode_b(cpi, tp, mi_row, mi_col + bs, output_enabled, c1, 1);
  } else if (bsl == bwl && bsl > bhl) {
    if (output_enabled)
      cpi->partition_count[pl][PARTITION_HORZ]++;
    encode_b(cpi, tp, mi_row,      mi_col, output_enabled, c1, 0);
    encode_b(cpi, tp, mi_row + bs, mi_col, output_enabled, c1, 1);
  } else {
    BLOCK_SIZE_TYPE subsize;
    int i;

    assert(bwl < bsl && bhl < bsl);
    if (level == BLOCK_SIZE_SB64X64) {
      subsize = BLOCK_SIZE_SB32X32;
#if CONFIG_SB8X8
    } else if (level == BLOCK_SIZE_SB32X32) {
      subsize = BLOCK_SIZE_MB16X16;
    } else {
      assert(level == BLOCK_SIZE_MB16X16);
      subsize = BLOCK_SIZE_SB8X8;
#else
    } else {
      assert(level == BLOCK_SIZE_SB32X32);
      subsize = BLOCK_SIZE_MB16X16;
#endif
    }

    if (output_enabled)
      cpi->partition_count[pl][PARTITION_SPLIT]++;

    for (i = 0; i < 4; i++) {
      const int x_idx = i & 1, y_idx = i >> 1;

      set_block_index(xd, i, subsize);
      encode_sb(cpi, tp, mi_row + y_idx * bs, mi_col + x_idx * bs,
                output_enabled, subsize,
#if CONFIG_SB8X8
                c2 ? c2[i] : c1, c3 ? c3[i] : NULL, NULL);
#else
                c2 ? c2[i] : c1, NULL);
#endif
    }
  }

#if CONFIG_SB8X8
  if (level > BLOCK_SIZE_SB8X8 &&
      (level == BLOCK_SIZE_MB16X16 || bsl == bwl || bsl == bhl))
#else
  if (level > BLOCK_SIZE_MB16X16 &&
      (level == BLOCK_SIZE_SB32X32 || bsl == bwl || bsl == bhl))
#endif
  {
    set_partition_seg_context(cpi, mi_row, mi_col);
    update_partition_context(xd, c1, level);
  }
}

static void encode_sb_row(VP9_COMP *cpi,
                          int mi_row,
                          TOKENEXTRA **tp,
                          int *totalrate) {
  VP9_COMMON *const cm = &cpi->common;
  MACROBLOCK *const x = &cpi->mb;
  MACROBLOCKD *const xd = &x->e_mbd;
  int mi_col, pl;

  // Initialize the left context for the new SB row
  vpx_memset(&cm->left_context, 0, sizeof(cm->left_context));
  vpx_memset(cm->left_seg_context, 0, sizeof(cm->left_seg_context));

  // Code each SB in the row
  for (mi_col = cm->cur_tile_mi_col_start;
       mi_col < cm->cur_tile_mi_col_end; mi_col += (4 << CONFIG_SB8X8)) {
    int i, p;
#if CONFIG_SB8X8
    BLOCK_SIZE_TYPE mb_partitioning[4][4];
#endif
    BLOCK_SIZE_TYPE sb_partitioning[4];
    BLOCK_SIZE_TYPE sb64_partitioning = BLOCK_SIZE_SB32X32;
    int sb64_rate = 0, sb64_dist = 0;
    int sb64_skip = 0;
    ENTROPY_CONTEXT l[16 * MAX_MB_PLANE], a[16 * MAX_MB_PLANE];
    PARTITION_CONTEXT seg_l[4], seg_a[4];
    TOKENEXTRA *tp_orig = *tp;

    for (p = 0; p < MAX_MB_PLANE; p++) {
      memcpy(a + 16 * p, cm->above_context[p] +
                 (mi_col * 4 >> (CONFIG_SB8X8 + xd->plane[p].subsampling_x)),
             sizeof(ENTROPY_CONTEXT) * 16 >> xd->plane[p].subsampling_x);
      memcpy(l + 16 * p, cm->left_context[p],
             sizeof(ENTROPY_CONTEXT) * 16 >> xd->plane[p].subsampling_y);
    }
    memcpy(&seg_a, cm->above_seg_context + (mi_col >> CONFIG_SB8X8),
           sizeof(seg_a));
    memcpy(&seg_l, cm->left_seg_context, sizeof(seg_l));

    // FIXME(rbultje): this function should probably be rewritten to be
    // recursive at some point in the future.
    for (i = 0; i < 4; i++) {
      const int x_idx = (i & 1) << (1 + CONFIG_SB8X8);
      const int y_idx = (i & 2) << CONFIG_SB8X8;
      int sb32_rate = 0, sb32_dist = 0;
      int splitmodes_used = 0;
      int sb32_skip = 0;
      int j;
      ENTROPY_CONTEXT l2[8 * MAX_MB_PLANE], a2[8 * MAX_MB_PLANE];

      sb_partitioning[i] = BLOCK_SIZE_MB16X16;
      if (mi_row + y_idx >= cm->mi_rows || mi_col + x_idx >= cm->mi_cols)
        continue;

      xd->sb_index = i;

      /* Function should not modify L & A contexts; save and restore on exit */
      for (p = 0; p < MAX_MB_PLANE; p++) {
        vpx_memcpy(l2 + 8 * p,
                   cm->left_context[p] +
                       (y_idx * 4 >> (CONFIG_SB8X8 +
                                      xd->plane[p].subsampling_y)),
                   sizeof(ENTROPY_CONTEXT) * 8 >> xd->plane[p].subsampling_y);
        vpx_memcpy(a2 + 8 * p,
                   cm->above_context[p] +
                       ((mi_col + x_idx) * 4 >> (CONFIG_SB8X8 +
                                                 xd->plane[p].subsampling_x)),
                   sizeof(ENTROPY_CONTEXT) * 8 >> xd->plane[p].subsampling_x);
      }

      /* Encode MBs in raster order within the SB */
      for (j = 0; j < 4; j++) {
        const int x_idx_m = x_idx + ((j & 1) << CONFIG_SB8X8);
        const int y_idx_m = y_idx + ((j >> 1) << CONFIG_SB8X8);
        int r, d;
#if CONFIG_SB8X8
        int r2, d2, mb16_rate = 0, mb16_dist = 0, k;
        ENTROPY_CONTEXT l3[4 * MAX_MB_PLANE], a3[4 * MAX_MB_PLANE];

        mb_partitioning[i][j] = BLOCK_SIZE_SB8X8;
#endif

        if (mi_row + y_idx_m >= cm->mi_rows ||
            mi_col + x_idx_m >= cm->mi_cols) {
          // MB lies outside frame, move on
          continue;
        }

        // Index of the MB in the SB 0..3
        xd->mb_index = j;

#if CONFIG_SB8X8
        for (p = 0; p < MAX_MB_PLANE; p++) {
          vpx_memcpy(l3 + 4 * p,
                     cm->left_context[p] +
                         (y_idx_m * 4 >> (CONFIG_SB8X8 +
                                          xd->plane[p].subsampling_y)),
                     sizeof(ENTROPY_CONTEXT) * 4 >> xd->plane[p].subsampling_y);
          vpx_memcpy(a3 + 4 * p,
                     cm->above_context[p] +
                         ((mi_col + x_idx_m) * 4 >> (CONFIG_SB8X8 +
                                                   xd->plane[p].subsampling_x)),
                     sizeof(ENTROPY_CONTEXT) * 4 >> xd->plane[p].subsampling_x);
        }

        for (k = 0; k < 4; k++) {
          xd->b_index = k;

          // try 8x8 coding
          pick_sb_modes(cpi, mi_row + y_idx_m + (k >> 1),
                        mi_col + x_idx_m + (k & 1),
                        tp, &r, &d, BLOCK_SIZE_SB8X8,
                        &x->sb8_context[xd->sb_index][xd->mb_index]
                                       [xd->b_index]);
          mb16_rate += r;
          mb16_dist += d;
          update_state(cpi, &x->sb8_context[xd->sb_index][xd->mb_index]
                                           [xd->b_index],
                       BLOCK_SIZE_SB8X8, 0);
          encode_superblock(cpi, tp,
                            0, mi_row + y_idx_m + (k >> 1),
                            mi_col + x_idx_m + (k & 1),
                            BLOCK_SIZE_SB8X8);
        }
        set_partition_seg_context(cpi, mi_row + y_idx_m, mi_col + x_idx_m);
        pl = partition_plane_context(xd, BLOCK_SIZE_MB16X16);
        mb16_rate += x->partition_cost[pl][PARTITION_SPLIT];
        for (p = 0; p < MAX_MB_PLANE; p++) {
          vpx_memcpy(cm->left_context[p] +
                         (y_idx_m * 4 >> (CONFIG_SB8X8 +
                                          xd->plane[p].subsampling_y)),
                     l3 + 4 * p,
                     sizeof(ENTROPY_CONTEXT) * 4 >> xd->plane[p].subsampling_y);
          vpx_memcpy(cm->above_context[p] +
                         ((mi_col + x_idx_m) * 4 >> (CONFIG_SB8X8 +
                                                   xd->plane[p].subsampling_x)),
                     a3 + 4 * p,
                     sizeof(ENTROPY_CONTEXT) * 4 >> xd->plane[p].subsampling_x);
        }

        // try 8x16 coding
        r2 = 0;
        d2 = 0;
        xd->b_index = 0;
        pick_sb_modes(cpi, mi_row + y_idx_m, mi_col + x_idx_m,
                      tp, &r, &d, BLOCK_SIZE_SB8X16,
                      &x->sb8x16_context[xd->sb_index][xd->mb_index]
                                        [xd->b_index]);
        r2 += r;
        d2 += d;
        update_state(cpi, &x->sb8x16_context[xd->sb_index][xd->mb_index]
                                            [xd->b_index],
                     BLOCK_SIZE_SB8X16, 0);
        encode_superblock(cpi, tp,
                          0, mi_row + y_idx_m, mi_col + x_idx_m,
                          BLOCK_SIZE_SB8X16);
        xd->b_index = 1;
        pick_sb_modes(cpi, mi_row + y_idx_m, mi_col + x_idx_m + 1,
                      tp, &r, &d, BLOCK_SIZE_SB8X16,
                      &x->sb8x16_context[xd->sb_index][xd->mb_index]
                                        [xd->b_index]);
        r2 += r;
        d2 += d;
        set_partition_seg_context(cpi, mi_row + y_idx_m, mi_col + x_idx_m);
        pl = partition_plane_context(xd, BLOCK_SIZE_MB16X16);
        r2 += x->partition_cost[pl][PARTITION_VERT];
        if (RDCOST(x->rdmult, x->rddiv, r2, d2) <
                RDCOST(x->rdmult, x->rddiv, mb16_rate, mb16_dist)) {
          mb16_rate = r;
          mb16_dist = d;
          mb_partitioning[i][j] = BLOCK_SIZE_SB8X16;
        }
        for (p = 0; p < MAX_MB_PLANE; p++) {
          vpx_memcpy(cm->left_context[p] +
                         (y_idx_m * 4 >> (CONFIG_SB8X8 +
                                          xd->plane[p].subsampling_y)),
                     l3 + 4 * p,
                     sizeof(ENTROPY_CONTEXT) * 4 >> xd->plane[p].subsampling_y);
          vpx_memcpy(cm->above_context[p] +
                         ((mi_col + x_idx_m) * 4 >> (CONFIG_SB8X8 +
                                                   xd->plane[p].subsampling_x)),
                     a3 + 4 * p,
                     sizeof(ENTROPY_CONTEXT) * 4 >> xd->plane[p].subsampling_x);
        }

        // try 16x8 coding
        r2 = 0;
        d2 = 0;
        xd->b_index = 0;
        pick_sb_modes(cpi, mi_row + y_idx_m, mi_col + x_idx_m,
                      tp, &r, &d, BLOCK_SIZE_SB16X8,
                      &x->sb16x8_context[xd->sb_index][xd->mb_index]
                                        [xd->b_index]);
        r2 += r;
        d2 += d;
        update_state(cpi, &x->sb16x8_context[xd->sb_index][xd->mb_index]
                                            [xd->b_index],
                     BLOCK_SIZE_SB16X8, 0);
        encode_superblock(cpi, tp,
                          0, mi_row + y_idx_m, mi_col + x_idx_m,
                          BLOCK_SIZE_SB16X8);
        xd->b_index = 1;
        pick_sb_modes(cpi, mi_row + y_idx_m + 1, mi_col + x_idx_m,
                      tp, &r, &d, BLOCK_SIZE_SB16X8,
                      &x->sb16x8_context[xd->sb_index][xd->mb_index]
                                        [xd->b_index]);
        r2 += r;
        d2 += d;
        set_partition_seg_context(cpi, mi_row + y_idx_m, mi_col + x_idx_m);
        pl = partition_plane_context(xd, BLOCK_SIZE_MB16X16);
        r2 += x->partition_cost[pl][PARTITION_HORZ];
        if (RDCOST(x->rdmult, x->rddiv, r2, d2) <
                RDCOST(x->rdmult, x->rddiv, mb16_rate, mb16_dist)) {
          mb16_rate = r;
          mb16_dist = d;
          mb_partitioning[i][j] = BLOCK_SIZE_SB16X8;
        }
        for (p = 0; p < MAX_MB_PLANE; p++) {
          vpx_memcpy(cm->left_context[p] +
                         (y_idx_m * 4 >> (CONFIG_SB8X8 +
                                          xd->plane[p].subsampling_y)),
                     l3 + 4 * p,
                     sizeof(ENTROPY_CONTEXT) * 4 >> xd->plane[p].subsampling_y);
          vpx_memcpy(cm->above_context[p] +
                         ((mi_col + x_idx_m) * 4 >> (CONFIG_SB8X8 +
                                                   xd->plane[p].subsampling_x)),
                     a3 + 4 * p,
                     sizeof(ENTROPY_CONTEXT) * 4 >> xd->plane[p].subsampling_x);
        }

        // try as 16x16
        pick_sb_modes(cpi, mi_row + y_idx_m, mi_col + x_idx_m,
                      tp, &r, &d, BLOCK_SIZE_MB16X16,
                      &x->mb_context[xd->sb_index][xd->mb_index]);
        set_partition_seg_context(cpi, mi_row + y_idx_m, mi_col + x_idx_m);
        pl = partition_plane_context(xd, BLOCK_SIZE_MB16X16);
        r += x->partition_cost[pl][PARTITION_NONE];
        if (RDCOST(x->rdmult, x->rddiv, r, d) <
                RDCOST(x->rdmult, x->rddiv, mb16_rate, mb16_dist)) {
          mb16_rate = r;
          mb16_dist = d;
          mb_partitioning[i][j] = BLOCK_SIZE_MB16X16;
        }
        sb32_rate += mb16_rate;
        sb32_dist += mb16_dist;
#else
        splitmodes_used += pick_mb_mode(cpi, mi_row + y_idx_m,
                                        mi_col + x_idx_m, tp, &r, &d);
        sb32_rate += r;
        sb32_dist += d;
#endif

        // Dummy encode, do not do the tokenization
#if CONFIG_SB8X8
        encode_sb(cpi, tp, mi_row + y_idx, mi_col + x_idx, 0,
                  BLOCK_SIZE_MB16X16, mb_partitioning[i][j], NULL, NULL);
#else
        encode_macroblock(cpi, tp, 0, mi_row + y_idx_m,
                          mi_col + x_idx_m);
#endif
      }

      /* Restore L & A coding context to those in place on entry */
      for (p = 0; p < MAX_MB_PLANE; p++) {
        vpx_memcpy(cm->left_context[p] +
                       (y_idx * 4 >> (CONFIG_SB8X8 +
                                      xd->plane[p].subsampling_y)),
                   l2 + 8 * p,
                   sizeof(ENTROPY_CONTEXT) * 8 >> xd->plane[p].subsampling_y);
        vpx_memcpy(cm->above_context[p] +
                       ((mi_col + x_idx) * 4 >> (CONFIG_SB8X8 +
                                                 xd->plane[p].subsampling_x)),
                   a2 + 8 * p,
                   sizeof(ENTROPY_CONTEXT) * 8 >> xd->plane[p].subsampling_x);
      }

      set_partition_seg_context(cpi, mi_row + y_idx, mi_col + x_idx);
      pl = partition_plane_context(xd, BLOCK_SIZE_SB32X32);
      sb32_rate += x->partition_cost[pl][PARTITION_SPLIT];

      if (cpi->sf.splitmode_breakout) {
        sb32_skip = splitmodes_used;
        sb64_skip += splitmodes_used;
      }

      // check 32x16
      if (mi_col + x_idx + (2 << CONFIG_SB8X8) <= cm->mi_cols) {
        int r, d;

        xd->mb_index = 0;
        pick_sb_modes(cpi, mi_row + y_idx, mi_col + x_idx,
                      tp, &r, &d, BLOCK_SIZE_SB32X16,
                      &x->sb32x16_context[xd->sb_index][xd->mb_index]);
        if (mi_row + y_idx + (1 << CONFIG_SB8X8) < cm->mi_rows) {
          int r2, d2;

          update_state(cpi, &x->sb32x16_context[xd->sb_index][xd->mb_index],
                       BLOCK_SIZE_SB32X16, 0);
          encode_superblock(cpi, tp,
                            0, mi_row + y_idx, mi_col + x_idx,
                            BLOCK_SIZE_SB32X16);
          xd->mb_index = 1;
          pick_sb_modes(cpi, mi_row + y_idx + (1 << CONFIG_SB8X8),
                        mi_col + x_idx, tp, &r2, &d2, BLOCK_SIZE_SB32X16,
                        &x->sb32x16_context[xd->sb_index][xd->mb_index]);
          r += r2;
          d += d2;
        }

        set_partition_seg_context(cpi, mi_row + y_idx, mi_col + x_idx);
        pl = partition_plane_context(xd, BLOCK_SIZE_SB32X32);
        r += x->partition_cost[pl][PARTITION_HORZ];

        /* is this better than MB coding? */
        if (RDCOST(x->rdmult, x->rddiv, r, d) <
                RDCOST(x->rdmult, x->rddiv, sb32_rate, sb32_dist)) {
          sb32_rate = r;
          sb32_dist = d;
          sb_partitioning[i] = BLOCK_SIZE_SB32X16;
        }

        for (p = 0; p < MAX_MB_PLANE; p++) {
          vpx_memcpy(cm->left_context[p] +
                         (y_idx * 4 >> (CONFIG_SB8X8 +
                                        xd->plane[p].subsampling_y)),
                     l2 + 8 * p,
                     sizeof(ENTROPY_CONTEXT) * 8 >> xd->plane[p].subsampling_y);
          vpx_memcpy(cm->above_context[p] +
                         ((mi_col + x_idx) * 4 >> (CONFIG_SB8X8 +
                                                   xd->plane[p].subsampling_x)),
                     a2 + 8 * p,
                     sizeof(ENTROPY_CONTEXT) * 8 >> xd->plane[p].subsampling_x);
        }
      }

      // check 16x32
      if (mi_row + y_idx + (2 << CONFIG_SB8X8) <= cm->mi_rows) {
        int r, d;

        xd->mb_index = 0;
        pick_sb_modes(cpi, mi_row + y_idx, mi_col + x_idx,
                      tp, &r, &d, BLOCK_SIZE_SB16X32,
                      &x->sb16x32_context[xd->sb_index][xd->mb_index]);
        if (mi_col + x_idx + (1 << CONFIG_SB8X8) < cm->mi_cols) {
          int r2, d2;

          update_state(cpi, &x->sb16x32_context[xd->sb_index][xd->mb_index],
                       BLOCK_SIZE_SB16X32, 0);
          encode_superblock(cpi, tp,
                            0, mi_row + y_idx, mi_col + x_idx,
                            BLOCK_SIZE_SB16X32);
          xd->mb_index = 1;
          pick_sb_modes(cpi, mi_row + y_idx,
                        mi_col + x_idx + (1 << CONFIG_SB8X8),
                        tp, &r2, &d2, BLOCK_SIZE_SB16X32,
                        &x->sb16x32_context[xd->sb_index][xd->mb_index]);
          r += r2;
          d += d2;
        }

        set_partition_seg_context(cpi, mi_row + y_idx, mi_col + x_idx);
        pl = partition_plane_context(xd, BLOCK_SIZE_SB32X32);
        r += x->partition_cost[pl][PARTITION_VERT];

        /* is this better than MB coding? */
        if (RDCOST(x->rdmult, x->rddiv, r, d) <
                RDCOST(x->rdmult, x->rddiv, sb32_rate, sb32_dist)) {
          sb32_rate = r;
          sb32_dist = d;
          sb_partitioning[i] = BLOCK_SIZE_SB16X32;
        }

        for (p = 0; p < MAX_MB_PLANE; p++) {
          vpx_memcpy(cm->left_context[p] +
                         (y_idx * 4 >> (CONFIG_SB8X8 +
                                        xd->plane[p].subsampling_y)),
                     l2 + 8 * p,
                     sizeof(ENTROPY_CONTEXT) * 8 >> xd->plane[p].subsampling_y);
          vpx_memcpy(cm->above_context[p] +
                         ((mi_col + x_idx) * 4 >> (CONFIG_SB8X8 +
                                                   xd->plane[p].subsampling_x)),
                     a2 + 8 * p,
                     sizeof(ENTROPY_CONTEXT) * 8 >> xd->plane[p].subsampling_x);
        }
      }

      if (!sb32_skip &&
          mi_col + x_idx + (2 << CONFIG_SB8X8) <= cm->mi_cols &&
          mi_row + y_idx + (2 << CONFIG_SB8X8) <= cm->mi_rows) {
        int r, d;

        /* Pick a mode assuming that it applies to all 4 of the MBs in the SB */
        pick_sb_modes(cpi, mi_row + y_idx, mi_col + x_idx,
                      tp, &r, &d, BLOCK_SIZE_SB32X32,
                      &x->sb32_context[xd->sb_index]);

        set_partition_seg_context(cpi, mi_row + y_idx, mi_col + x_idx);
        pl = partition_plane_context(xd, BLOCK_SIZE_SB32X32);
        r += x->partition_cost[pl][PARTITION_NONE];

        if (RDCOST(x->rdmult, x->rddiv, r, d) <
                RDCOST(x->rdmult, x->rddiv, sb32_rate, sb32_dist)) {
          sb32_rate = r;
          sb32_dist = d;
          sb_partitioning[i] = BLOCK_SIZE_SB32X32;
        }
      }

      // If we used 16x16 instead of 32x32 then skip 64x64 (if enabled).
      if (cpi->sf.mb16_breakout && sb_partitioning[i] != BLOCK_SIZE_SB32X32) {
        ++sb64_skip;
      }

      sb64_rate += sb32_rate;
      sb64_dist += sb32_dist;

      /* Encode SB using best computed mode(s) */
      // FIXME(rbultje): there really shouldn't be any need to encode_mb/sb
      // for each level that we go up, we can just keep tokens and recon
      // pixels of the lower level; also, inverting SB/MB order (big->small
      // instead of small->big) means we can use as threshold for small, which
      // may enable breakouts if RD is not good enough (i.e. faster)
      encode_sb(cpi, tp, mi_row + y_idx, mi_col + x_idx, 0,
#if CONFIG_SB8X8
                BLOCK_SIZE_SB32X32, sb_partitioning[i], mb_partitioning[i],
                NULL);
#else
                BLOCK_SIZE_SB32X32, sb_partitioning[i], NULL);
#endif
    }

    for (p = 0; p < MAX_MB_PLANE; p++) {
      memcpy(cm->above_context[p] +
                 (mi_col * 4 >> (CONFIG_SB8X8 + xd->plane[p].subsampling_x)),
             a + 16 * p,
             sizeof(ENTROPY_CONTEXT) * 16 >> xd->plane[p].subsampling_x);
      memcpy(cm->left_context[p], l + 16 * p,
             sizeof(ENTROPY_CONTEXT) * 16 >> xd->plane[p].subsampling_y);
    }
    memcpy(cm->above_seg_context + (mi_col >> CONFIG_SB8X8), &seg_a,
           sizeof(seg_a));
    memcpy(cm->left_seg_context, &seg_l, sizeof(seg_l));

    set_partition_seg_context(cpi, mi_row, mi_col);
    pl = partition_plane_context(xd, BLOCK_SIZE_SB64X64);
    sb64_rate += x->partition_cost[pl][PARTITION_SPLIT];

    // check 64x32
    if (mi_col + (4 << CONFIG_SB8X8) <= cm->mi_cols && !(cm->mb_rows & 1)) {
      int r, d;

      xd->sb_index = 0;
      pick_sb_modes(cpi, mi_row, mi_col,
                    tp, &r, &d, BLOCK_SIZE_SB64X32,
                    &x->sb64x32_context[xd->sb_index]);
      if (mi_row + (2 << CONFIG_SB8X8) != cm->mi_rows) {
        int r2, d2;

        update_state(cpi, &x->sb64x32_context[xd->sb_index],
                     BLOCK_SIZE_SB64X32, 0);
        encode_superblock(cpi, tp,
                          0, mi_row, mi_col, BLOCK_SIZE_SB64X32);
        xd->sb_index = 1;
        pick_sb_modes(cpi, mi_row + (2 << CONFIG_SB8X8), mi_col,
                      tp, &r2, &d2, BLOCK_SIZE_SB64X32,
                      &x->sb64x32_context[xd->sb_index]);
        r += r2;
        d += d2;
      }

      set_partition_seg_context(cpi, mi_row, mi_col);
      pl = partition_plane_context(xd, BLOCK_SIZE_SB64X64);
      r += x->partition_cost[pl][PARTITION_HORZ];

      /* is this better than MB coding? */
      if (RDCOST(x->rdmult, x->rddiv, r, d) <
              RDCOST(x->rdmult, x->rddiv, sb64_rate, sb64_dist)) {
        sb64_rate = r;
        sb64_dist = d;
        sb64_partitioning = BLOCK_SIZE_SB64X32;
      }

      for (p = 0; p < MAX_MB_PLANE; p++) {
        memcpy(cm->above_context[p] +
                   (mi_col * 4 >> (CONFIG_SB8X8 + xd->plane[p].subsampling_x)),
               a + 16 * p,
               sizeof(ENTROPY_CONTEXT) * 16 >> xd->plane[p].subsampling_x);
        memcpy(cm->left_context[p], l + 16 * p,
               sizeof(ENTROPY_CONTEXT) * 16 >> xd->plane[p].subsampling_y);
      }
    }

    // check 32x64
    if (mi_row + (4 << CONFIG_SB8X8) <= cm->mi_rows && !(cm->mb_cols & 1)) {
      int r, d;

      xd->sb_index = 0;
      pick_sb_modes(cpi, mi_row, mi_col,
                    tp, &r, &d, BLOCK_SIZE_SB32X64,
                    &x->sb32x64_context[xd->sb_index]);
      if (mi_col + (2 << CONFIG_SB8X8) != cm->mi_cols) {
        int r2, d2;

        update_state(cpi, &x->sb32x64_context[xd->sb_index],
                     BLOCK_SIZE_SB32X64, 0);
        encode_superblock(cpi, tp,
                          0, mi_row, mi_col, BLOCK_SIZE_SB32X64);
        xd->sb_index = 1;
        pick_sb_modes(cpi, mi_row, mi_col + (2 << CONFIG_SB8X8),
                      tp, &r2, &d2, BLOCK_SIZE_SB32X64,
                      &x->sb32x64_context[xd->sb_index]);
        r += r2;
        d += d2;
      }

      set_partition_seg_context(cpi, mi_row, mi_col);
      pl = partition_plane_context(xd, BLOCK_SIZE_SB64X64);
      r += x->partition_cost[pl][PARTITION_VERT];

      /* is this better than MB coding? */
      if (RDCOST(x->rdmult, x->rddiv, r, d) <
              RDCOST(x->rdmult, x->rddiv, sb64_rate, sb64_dist)) {
        sb64_rate = r;
        sb64_dist = d;
        sb64_partitioning = BLOCK_SIZE_SB32X64;
      }

      for (p = 0; p < MAX_MB_PLANE; p++) {
        memcpy(cm->above_context[p] +
                   (mi_col * 4 >> (CONFIG_SB8X8 + xd->plane[p].subsampling_x)),
               a + 16 * p,
               sizeof(ENTROPY_CONTEXT) * 16 >> xd->plane[p].subsampling_x);
        memcpy(cm->left_context[p], l + 16 * p,
               sizeof(ENTROPY_CONTEXT) * 16 >> xd->plane[p].subsampling_y);
      }
    }

    if (!sb64_skip &&
        mi_col + (4 << CONFIG_SB8X8) <= cm->mi_cols &&
        mi_row + (4 << CONFIG_SB8X8) <= cm->mi_rows) {
      int r, d;

      pick_sb_modes(cpi, mi_row, mi_col, tp, &r, &d,
                    BLOCK_SIZE_SB64X64, &x->sb64_context);

      set_partition_seg_context(cpi, mi_row, mi_col);
      pl = partition_plane_context(xd, BLOCK_SIZE_SB64X64);
      r += x->partition_cost[pl][PARTITION_NONE];

      if (RDCOST(x->rdmult, x->rddiv, r, d) <
              RDCOST(x->rdmult, x->rddiv, sb64_rate, sb64_dist)) {
        sb64_rate = r;
        sb64_dist = d;
        sb64_partitioning = BLOCK_SIZE_SB64X64;
      }
    }

    assert(tp_orig == *tp);
    encode_sb(cpi, tp, mi_row, mi_col, 1, BLOCK_SIZE_SB64X64,
#if CONFIG_SB8X8
              sb64_partitioning, sb_partitioning, mb_partitioning);
#else
              sb64_partitioning, sb_partitioning);
#endif
    assert(tp_orig < *tp);
  }
}

static void init_encode_frame_mb_context(VP9_COMP *cpi) {
  MACROBLOCK *const x = &cpi->mb;
  VP9_COMMON *const cm = &cpi->common;
  MACROBLOCKD *const xd = &x->e_mbd;

  x->act_zbin_adj = 0;
  cpi->seg0_idx = 0;
  vpx_memset(cpi->ref_pred_count, 0, sizeof(cpi->ref_pred_count));

  xd->mode_info_stride = cm->mode_info_stride;
  xd->frame_type = cm->frame_type;

  xd->frames_since_golden = cm->frames_since_golden;
  xd->frames_till_alt_ref_frame = cm->frames_till_alt_ref_frame;

  // reset intra mode contexts
  if (cm->frame_type == KEY_FRAME)
    vp9_init_mbmode_probs(cm);

  // Copy data over into macro block data structures.
  vp9_setup_src_planes(x, cpi->Source, 0, 0);

  // TODO(jkoleszar): are these initializations required?
  setup_pre_planes(xd, &cm->yv12_fb[cm->ref_frame_map[cpi->lst_fb_idx]], NULL,
                   0, 0, NULL, NULL);
  setup_dst_planes(xd, &cm->yv12_fb[cm->new_fb_idx], 0, 0);

  // set up frame for intra coded blocks
  vp9_setup_intra_recon(&cm->yv12_fb[cm->new_fb_idx]);

  vp9_build_block_offsets(x);

  vp9_setup_block_dptrs(&x->e_mbd);

  xd->mode_info_context->mbmi.mode = DC_PRED;
  xd->mode_info_context->mbmi.uv_mode = DC_PRED;

  vp9_zero(cpi->count_mb_ref_frame_usage)
  vp9_zero(cpi->bmode_count)
  vp9_zero(cpi->ymode_count)
#if !CONFIG_SB8X8
  vp9_zero(cpi->i8x8_mode_count)
#endif
  vp9_zero(cpi->y_uv_mode_count)
  vp9_zero(cpi->sub_mv_ref_count)
#if !CONFIG_SB8X8
  vp9_zero(cpi->mbsplit_count)
#endif
  vp9_zero(cpi->common.fc.mv_ref_ct)
  vp9_zero(cpi->sb_ymode_count)
  vp9_zero(cpi->partition_count);

#if CONFIG_COMP_INTERINTRA_PRED
  vp9_zero(cpi->interintra_count);
  vp9_zero(cpi->interintra_select_count);
#endif

  // Note: this memset assumes above_context[0], [1] and [2]
  // are allocated as part of the same buffer.
  vpx_memset(cm->above_context[0], 0, sizeof(ENTROPY_CONTEXT) * 4 *
                                      MAX_MB_PLANE * mb_cols_aligned_to_sb(cm));
  vpx_memset(cm->above_seg_context, 0, sizeof(PARTITION_CONTEXT) *
                                       mb_cols_aligned_to_sb(cm));
}

static void switch_lossless_mode(VP9_COMP *cpi, int lossless) {
  if (lossless) {
    cpi->mb.fwd_txm8x4            = vp9_short_walsh8x4;
    cpi->mb.fwd_txm4x4            = vp9_short_walsh4x4;
    cpi->mb.e_mbd.inv_txm4x4_1    = vp9_short_iwalsh4x4_1;
    cpi->mb.e_mbd.inv_txm4x4      = vp9_short_iwalsh4x4;
    cpi->mb.optimize              = 0;
    cpi->common.filter_level      = 0;
    cpi->zbin_mode_boost_enabled  = 0;
    cpi->common.txfm_mode         = ONLY_4X4;
  } else {
    cpi->mb.fwd_txm8x4            = vp9_short_fdct8x4;
    cpi->mb.fwd_txm4x4            = vp9_short_fdct4x4;
    cpi->mb.e_mbd.inv_txm4x4_1    = vp9_short_idct4x4_1;
    cpi->mb.e_mbd.inv_txm4x4      = vp9_short_idct4x4;
  }
}


static void encode_frame_internal(VP9_COMP *cpi) {
  int mi_row;
  MACROBLOCK *const x = &cpi->mb;
  VP9_COMMON *const cm = &cpi->common;
  MACROBLOCKD *const xd = &x->e_mbd;
  int totalrate;

//  fprintf(stderr, "encode_frame_internal frame %d (%d) type %d\n",
//           cpi->common.current_video_frame, cpi->common.show_frame,
//           cm->frame_type);

  // Compute a modified set of reference frame probabilities to use when
  // prediction fails. These are based on the current general estimates for
  // this frame which may be updated with each iteration of the recode loop.
  vp9_compute_mod_refprobs(cm);

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

  // Reset frame count of inter 0,0 motion vector usage.
  cpi->inter_zz_count = 0;

  cpi->skip_true_count[0] = cpi->skip_true_count[1] = cpi->skip_true_count[2] = 0;
  cpi->skip_false_count[0] = cpi->skip_false_count[1] = cpi->skip_false_count[2] = 0;

  vp9_zero(cpi->switchable_interp_count);
  vp9_zero(cpi->best_switchable_interp_count);

  xd->mode_info_context = cm->mi;
  xd->prev_mode_info_context = cm->prev_mi;

  vp9_zero(cpi->NMVcount);
  vp9_zero(cpi->coef_counts_4x4);
  vp9_zero(cpi->coef_counts_8x8);
  vp9_zero(cpi->coef_counts_16x16);
  vp9_zero(cpi->coef_counts_32x32);
  vp9_zero(cm->fc.eob_branch_counts);
#if CONFIG_CODE_ZEROGROUP
  vp9_zero(cm->fc.zpc_counts_4x4);
  vp9_zero(cm->fc.zpc_counts_8x8);
  vp9_zero(cm->fc.zpc_counts_16x16);
  vp9_zero(cm->fc.zpc_counts_32x32);
#endif

  cpi->mb.e_mbd.lossless = (cm->base_qindex == 0 &&
                            cm->y_dc_delta_q == 0 &&
                            cm->uv_dc_delta_q == 0 &&
                            cm->uv_ac_delta_q == 0);
  switch_lossless_mode(cpi, cpi->mb.e_mbd.lossless);

  vp9_frame_init_quantizer(cpi);

  vp9_initialize_rd_consts(cpi, cm->base_qindex + cm->y_dc_delta_q);
  vp9_initialize_me_consts(cpi, cm->base_qindex);

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
  vpx_memset(cpi->txfm_count_32x32p, 0, sizeof(cpi->txfm_count_32x32p));
  vpx_memset(cpi->txfm_count_16x16p, 0, sizeof(cpi->txfm_count_16x16p));
  vpx_memset(cpi->txfm_count_8x8p, 0, sizeof(cpi->txfm_count_8x8p));
  vpx_memset(cpi->rd_tx_select_diff, 0, sizeof(cpi->rd_tx_select_diff));
  {
    struct vpx_usec_timer  emr_timer;
    vpx_usec_timer_start(&emr_timer);

    {
      // Take tiles into account and give start/end MB
      int tile_col, tile_row;
      TOKENEXTRA *tp = cpi->tok;

      for (tile_row = 0; tile_row < cm->tile_rows; tile_row++) {
        vp9_get_tile_row_offsets(cm, tile_row);

        for (tile_col = 0; tile_col < cm->tile_columns; tile_col++) {
          TOKENEXTRA *tp_old = tp;

          // For each row of SBs in the frame
          vp9_get_tile_col_offsets(cm, tile_col);
          for (mi_row = cm->cur_tile_mi_row_start;
               mi_row < cm->cur_tile_mi_row_end;
               mi_row += (4 << CONFIG_SB8X8)) {
            encode_sb_row(cpi, mi_row, &tp, &totalrate);
          }
          cpi->tok_count[tile_col] = (unsigned int)(tp - tp_old);
          assert(tp - cpi->tok <=
                 get_token_alloc(cm->mb_rows, cm->mb_cols));
        }
      }
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

static int check_dual_ref_flags(VP9_COMP *cpi) {
  MACROBLOCKD *xd = &cpi->mb.e_mbd;
  int ref_flags = cpi->ref_frame_flags;

  if (vp9_segfeature_active(xd, 1, SEG_LVL_REF_FRAME)) {
    if ((ref_flags & (VP9_LAST_FLAG | VP9_GOLD_FLAG)) == (VP9_LAST_FLAG | VP9_GOLD_FLAG) &&
        vp9_check_segref(xd, 1, LAST_FRAME))
      return 1;
    if ((ref_flags & (VP9_GOLD_FLAG | VP9_ALT_FLAG)) == (VP9_GOLD_FLAG | VP9_ALT_FLAG) &&
        vp9_check_segref(xd, 1, GOLDEN_FRAME))
      return 1;
    if ((ref_flags & (VP9_ALT_FLAG  | VP9_LAST_FLAG)) == (VP9_ALT_FLAG  | VP9_LAST_FLAG) &&
        vp9_check_segref(xd, 1, ALTREF_FRAME))
      return 1;
    return 0;
  } else {
    return (!!(ref_flags & VP9_GOLD_FLAG) +
            !!(ref_flags & VP9_LAST_FLAG) +
            !!(ref_flags & VP9_ALT_FLAG)) >= 2;
  }
}

static int get_skip_flag(MODE_INFO *mi, int mis, int ymbs, int xmbs) {
  int x, y;

  for (y = 0; y < ymbs; y++) {
    for (x = 0; x < xmbs; x++) {
      if (!mi[y * mis + x].mbmi.mb_skip_coeff)
        return 0;
    }
  }

  return 1;
}

static void set_txfm_flag(MODE_INFO *mi, int mis, int ymbs, int xmbs,
                          TX_SIZE txfm_size) {
  int x, y;

  for (y = 0; y < ymbs; y++) {
    for (x = 0; x < xmbs; x++)
      mi[y * mis + x].mbmi.txfm_size = txfm_size;
  }
}

static void reset_skip_txfm_size_b(VP9_COMP *cpi, MODE_INFO *mi,
                                   int mis, TX_SIZE txfm_max,
                                   int bw, int bh, int mi_row, int mi_col,
                                   BLOCK_SIZE_TYPE bsize) {
  VP9_COMMON *const cm = &cpi->common;
  MB_MODE_INFO *const mbmi = &mi->mbmi;

  if (mi_row >= cm->mi_rows || mi_col >= cm->mi_cols)
    return;

  if (mbmi->txfm_size > txfm_max) {
    MACROBLOCK *const x = &cpi->mb;
    MACROBLOCKD *const xd = &x->e_mbd;
    const int segment_id = mbmi->segment_id;
    const int ymbs = MIN(bh, cm->mi_rows - mi_row);
    const int xmbs = MIN(bw, cm->mi_cols - mi_col);

    xd->mode_info_context = mi;
    assert(vp9_segfeature_active(xd, segment_id, SEG_LVL_SKIP) ||
           get_skip_flag(mi, mis, ymbs, xmbs));
    set_txfm_flag(mi, mis, ymbs, xmbs, txfm_max);
  }
}

static void reset_skip_txfm_size_sb(VP9_COMP *cpi, MODE_INFO *mi,
                                    TX_SIZE txfm_max,
                                    int mi_row, int mi_col,
                                    BLOCK_SIZE_TYPE bsize) {
  VP9_COMMON *const cm = &cpi->common;
  const int mis = cm->mode_info_stride;
  int bwl, bhl;
  const int bsl = mi_width_log2(bsize), bs = 1 << (bsl - 1);

  if (mi_row >= cm->mi_rows || mi_col >= cm->mi_cols)
    return;

  bwl = mi_width_log2(mi->mbmi.sb_type);
  bhl = mi_height_log2(mi->mbmi.sb_type);

  if (bwl == bsl && bhl == bsl) {
    reset_skip_txfm_size_b(cpi, mi, mis, txfm_max, 1 << bsl, 1 << bsl,
                           mi_row, mi_col, bsize);
  } else if (bwl == bsl && bhl < bsl) {
    reset_skip_txfm_size_b(cpi, mi, mis, txfm_max, 1 << bsl, bs,
                           mi_row, mi_col, bsize);
    reset_skip_txfm_size_b(cpi, mi + bs * mis, mis, txfm_max, 1 << bsl, bs,
                           mi_row + bs, mi_col, bsize);
  } else if (bwl < bsl && bhl == bsl) {
    reset_skip_txfm_size_b(cpi, mi, mis, txfm_max, bs, 1 << bsl,
                           mi_row, mi_col, bsize);
    reset_skip_txfm_size_b(cpi, mi + bs, mis, txfm_max, bs, 1 << bsl,
                           mi_row, mi_col + bs, bsize);
  } else {
    BLOCK_SIZE_TYPE subsize;
    int n;

    assert(bwl < bsl && bhl < bsl);
    if (bsize == BLOCK_SIZE_SB64X64) {
      subsize = BLOCK_SIZE_SB32X32;
#if CONFIG_SB8X8
    } else if (bsize == BLOCK_SIZE_SB32X32) {
      subsize = BLOCK_SIZE_MB16X16;
    } else {
      assert(bsize == BLOCK_SIZE_MB16X16);
      subsize = BLOCK_SIZE_SB8X8;
#else
    } else {
      assert(bsize == BLOCK_SIZE_SB32X32);
      subsize = BLOCK_SIZE_MB16X16;
#endif
    }

    for (n = 0; n < 4; n++) {
      const int y_idx = n >> 1, x_idx = n & 0x01;

      reset_skip_txfm_size_sb(cpi, mi + y_idx * bs * mis + x_idx * bs,
                              txfm_max, mi_row + y_idx * bs,
                              mi_col + x_idx * bs, subsize);
    }
  }
}

static void reset_skip_txfm_size(VP9_COMP *cpi, TX_SIZE txfm_max) {
  VP9_COMMON *const cm = &cpi->common;
  int mi_row, mi_col;
  const int mis = cm->mode_info_stride;
  MODE_INFO *mi, *mi_ptr = cm->mi;

  for (mi_row = 0; mi_row < cm->mi_rows;
       mi_row += (4 << CONFIG_SB8X8), mi_ptr += (4 << CONFIG_SB8X8) * mis) {
    mi = mi_ptr;
    for (mi_col = 0; mi_col < cm->mi_cols;
         mi_col += (4 << CONFIG_SB8X8), mi += (4 << CONFIG_SB8X8)) {
      reset_skip_txfm_size_sb(cpi, mi, txfm_max,
                              mi_row, mi_col, BLOCK_SIZE_SB64X64);
    }
  }
}

void vp9_encode_frame(VP9_COMP *cpi) {
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
    else if (cpi->is_src_frame_alt_ref && cpi->refresh_golden_frame)
      frame_type = 3;
    else if (cpi->refresh_golden_frame || cpi->refresh_alt_ref_frame)
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

    cpi->mb.e_mbd.lossless = 0;
    if (cpi->oxcf.lossless) {
      txfm_type = ONLY_4X4;
      cpi->mb.e_mbd.lossless = 1;
    } else
#if 0
    /* FIXME (rbultje): this code is disabled until we support cost updates
     * while a frame is being encoded; the problem is that each time we
     * "revert" to 4x4 only (or even 8x8 only), the coefficient probabilities
     * for 16x16 (and 8x8) start lagging behind, thus leading to them lagging
     * further behind and not being chosen for subsequent frames either. This
     * is essentially a local minimum problem that we can probably fix by
     * estimating real costs more closely within a frame, perhaps by re-
     * calculating costs on-the-fly as frame encoding progresses. */
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
    txfm_type = cpi->rd_tx_select_threshes[frame_type][ALLOW_32X32] >=
                  cpi->rd_tx_select_threshes[frame_type][TX_MODE_SELECT] ?
                    ALLOW_32X32 : TX_MODE_SELECT;
#endif
    cpi->common.txfm_mode = txfm_type;
    if (txfm_type != TX_MODE_SELECT) {
      cpi->common.prob_tx[0] = 128;
      cpi->common.prob_tx[1] = 128;
    }
    cpi->common.comp_pred_mode = pred_type;
    encode_frame_internal(cpi);

    for (i = 0; i < NB_PREDICTION_TYPES; ++i) {
      const int diff = (int)(cpi->rd_comp_pred_diff[i] / cpi->common.MBs);
      cpi->rd_prediction_type_threshes[frame_type][i] += diff;
      cpi->rd_prediction_type_threshes[frame_type][i] >>= 1;
    }

    for (i = 0; i < NB_TXFM_MODES; ++i) {
      int64_t pd = cpi->rd_tx_select_diff[i];
      int diff;
      if (i == TX_MODE_SELECT)
        pd -= RDCOST(cpi->mb.rdmult, cpi->mb.rddiv,
                     2048 * (TX_SIZE_MAX_SB - 1), 0);
      diff = (int)(pd / cpi->common.MBs);
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
      const int count4x4 = cpi->txfm_count_16x16p[TX_4X4] +
                           cpi->txfm_count_32x32p[TX_4X4] +
                           cpi->txfm_count_8x8p[TX_4X4];
      const int count8x8_lp = cpi->txfm_count_32x32p[TX_8X8] +
                              cpi->txfm_count_16x16p[TX_8X8];
      const int count8x8_8x8p = cpi->txfm_count_8x8p[TX_8X8];
      const int count16x16_16x16p = cpi->txfm_count_16x16p[TX_16X16];
      const int count16x16_lp = cpi->txfm_count_32x32p[TX_16X16];
      const int count32x32 = cpi->txfm_count_32x32p[TX_32X32];

      if (count4x4 == 0 && count16x16_lp == 0 && count16x16_16x16p == 0 &&
          count32x32 == 0) {
        cpi->common.txfm_mode = ALLOW_8X8;
        reset_skip_txfm_size(cpi, TX_8X8);
      } else if (count8x8_8x8p == 0 && count16x16_16x16p == 0 &&
                 count8x8_lp == 0 && count16x16_lp == 0 && count32x32 == 0) {
        cpi->common.txfm_mode = ONLY_4X4;
        reset_skip_txfm_size(cpi, TX_4X4);
      } else if (count8x8_lp == 0 && count16x16_lp == 0 && count4x4 == 0) {
        cpi->common.txfm_mode = ALLOW_32X32;
      } else if (count32x32 == 0 && count8x8_lp == 0 && count4x4 == 0) {
        cpi->common.txfm_mode = ALLOW_16X16;
        reset_skip_txfm_size(cpi, TX_16X16);
      }
    }

    // Update interpolation filter strategy for next frame.
    if ((cpi->common.frame_type != KEY_FRAME) && (cpi->sf.search_best_filter))
      vp9_select_interp_filter_type(cpi);
  } else {
    encode_frame_internal(cpi);
  }

}

void vp9_build_block_offsets(MACROBLOCK *x) {
}

static void sum_intra_stats(VP9_COMP *cpi, MACROBLOCK *x) {
  const MACROBLOCKD *xd = &x->e_mbd;
  const MB_PREDICTION_MODE m = xd->mode_info_context->mbmi.mode;
  const MB_PREDICTION_MODE uvm = xd->mode_info_context->mbmi.uv_mode;

#ifdef MODE_STATS
  const int is_key = cpi->common.frame_type == KEY_FRAME;

  ++ (is_key ? uv_modes : inter_uv_modes)[uvm];
  ++ uv_modes_y[m][uvm];

  if (m == I4X4_PRED) {
    unsigned int *const bct = is_key ? b_modes : inter_b_modes;

    int b = 0;

    do {
      ++ bct[xd->block[b].bmi.as_mode.first];
    } while (++b < (16 >> (CONFIG_SB8X8 * 2)));
  }

#if !CONFIG_SB8X8
  if (m == I8X8_PRED) {
    i8x8_modes[xd->block[0].bmi.as_mode.first]++;
    i8x8_modes[xd->block[2].bmi.as_mode.first]++;
    i8x8_modes[xd->block[8].bmi.as_mode.first]++;
    i8x8_modes[xd->block[10].bmi.as_mode.first]++;
  }
#endif
#endif

#if CONFIG_SB8X8
  if (xd->mode_info_context->mbmi.sb_type > BLOCK_SIZE_SB8X8)
#else
  if (xd->mode_info_context->mbmi.sb_type > BLOCK_SIZE_MB16X16)
#endif
  {
    ++cpi->sb_ymode_count[m];
  } else {
    ++cpi->ymode_count[m];
  }
#if !CONFIG_SB8X8
  if (m != I8X8_PRED)
#endif
    ++cpi->y_uv_mode_count[m][uvm];
#if !CONFIG_SB8X8
  else {
    cpi->i8x8_mode_count[xd->mode_info_context->bmi[0].as_mode.first]++;
    cpi->i8x8_mode_count[xd->mode_info_context->bmi[2].as_mode.first]++;
    cpi->i8x8_mode_count[xd->mode_info_context->bmi[8].as_mode.first]++;
    cpi->i8x8_mode_count[xd->mode_info_context->bmi[10].as_mode.first]++;
  }
#endif
  if (m == I4X4_PRED) {
    int b = 0;
    do {
      int m = xd->mode_info_context->bmi[b].as_mode.first;
#if CONFIG_NEWBINTRAMODES
      if (m == B_CONTEXT_PRED) m -= CONTEXT_PRED_REPLACEMENTS;
#endif
      ++cpi->bmode_count[m];
    } while (++b < (16 >> (CONFIG_SB8X8 * 2)));
  }
}

// Experimental stub function to create a per MB zbin adjustment based on
// some previously calculated measure of MB activity.
static void adjust_act_zbin(VP9_COMP *cpi, MACROBLOCK *x) {
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

#if !CONFIG_SB8X8
static void encode_macroblock(VP9_COMP *cpi, TOKENEXTRA **t,
                              int output_enabled,
                              int mi_row, int mi_col) {
  VP9_COMMON *const cm = &cpi->common;
  MACROBLOCK *const x = &cpi->mb;
  MACROBLOCKD *const xd = &x->e_mbd;
  MODE_INFO *mi = xd->mode_info_context;
  MB_MODE_INFO *const mbmi = &mi->mbmi;
  const int mis = cm->mode_info_stride;
#if CONFIG_SB8X8
  int n;
#endif

  assert(xd->mode_info_context->mbmi.sb_type == BLOCK_SIZE_MB16X16);

#ifdef ENC_DEBUG
  enc_debug = (cpi->common.current_video_frame == 11 && cm->show_frame &&
               mb_row == 8 && mb_col == 0 && output_enabled);
  if (enc_debug)
    printf("Encode MB %d %d output %d\n", mb_row, mb_col, output_enabled);
#endif
  if (cm->frame_type == KEY_FRAME) {
    if (cpi->oxcf.tuning == VP8_TUNE_SSIM && output_enabled) {
      // Adjust the zbin based on this MB rate.
      adjust_act_zbin(cpi, x);
      vp9_update_zbin_extra(cpi, x);
    }
  } else {
    vp9_setup_interp_filters(xd, mbmi->interp_filter, cm);

    if (cpi->oxcf.tuning == VP8_TUNE_SSIM) {
      // Adjust the zbin based on this MB rate.
      adjust_act_zbin(cpi, x);
    }

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
          cpi->zbin_mode_boost = SPLIT_MV_ZBIN_BOOST;
        else
          cpi->zbin_mode_boost = MV_ZBIN_BOOST;
      } else {
        cpi->zbin_mode_boost = INTRA_ZBIN_BOOST;
      }
    }

    vp9_update_zbin_extra(cpi, x);
  }

  if (mbmi->ref_frame == INTRA_FRAME) {
#if 0  // def ENC_DEBUG
    if (enc_debug) {
      printf("Mode %d skip %d tx_size %d\n", mbmi->mode, x->skip,
             mbmi->txfm_size);
    }
#endif
    if (mbmi->mode == I4X4_PRED) {
      vp9_encode_intra16x16mbuv(cm, x);
      vp9_encode_intra4x4mby(x, BLOCK_SIZE_MB16X16);
    } else if (mbmi->mode == I8X8_PRED) {
      vp9_encode_intra8x8mby(x);
      vp9_encode_intra8x8mbuv(x);
    } else {
      vp9_encode_intra16x16mbuv(cm, x);
      vp9_encode_intra16x16mby(cm, x);
    }

    if (output_enabled)
      sum_intra_stats(cpi, x);
  } else {
    int ref_fb_idx, second_ref_fb_idx;
#ifdef ENC_DEBUG
    if (enc_debug)
      printf("Mode %d skip %d tx_size %d ref %d ref2 %d mv %d %d interp %d\n",
             mbmi->mode, x->skip, mbmi->txfm_size,
             mbmi->ref_frame, mbmi->second_ref_frame,
             mbmi->mv[0].as_mv.row, mbmi->mv[0].as_mv.col,
             mbmi->interp_filter);
#endif

    assert(cm->frame_type != KEY_FRAME);

    if (mbmi->ref_frame == LAST_FRAME)
      ref_fb_idx = cpi->common.ref_frame_map[cpi->lst_fb_idx];
    else if (mbmi->ref_frame == GOLDEN_FRAME)
      ref_fb_idx = cpi->common.ref_frame_map[cpi->gld_fb_idx];
    else
      ref_fb_idx = cpi->common.ref_frame_map[cpi->alt_fb_idx];

    if (mbmi->second_ref_frame > 0) {
      if (mbmi->second_ref_frame == LAST_FRAME)
        second_ref_fb_idx = cpi->common.ref_frame_map[cpi->lst_fb_idx];
      else if (mbmi->second_ref_frame == GOLDEN_FRAME)
        second_ref_fb_idx = cpi->common.ref_frame_map[cpi->gld_fb_idx];
      else
        second_ref_fb_idx = cpi->common.ref_frame_map[cpi->alt_fb_idx];
    }

    setup_pre_planes(xd,
        &cpi->common.yv12_fb[ref_fb_idx],
        mbmi->second_ref_frame > 0 ? &cpi->common.yv12_fb[second_ref_fb_idx]
                                   : NULL,
        mi_row, mi_col, xd->scale_factor, xd->scale_factor_uv);

    if (!x->skip) {
      vp9_build_inter_predictors_sb(xd, mi_row, mi_col, BLOCK_SIZE_MB16X16);
      vp9_encode_sb(cm, x, BLOCK_SIZE_MB16X16);
    } else {
      vp9_build_inter_predictors_sb(xd, mi_row, mi_col, BLOCK_SIZE_MB16X16);
#if CONFIG_COMP_INTERINTRA_PRED
      if (xd->mode_info_context->mbmi.second_ref_frame == INTRA_FRAME) {
        vp9_build_interintra_predictors(xd,
                                        xd->plane[0].dst.buf,
                                        xd->plane[1].dst.buf,
                                        xd->plane[2].dst.buf,
                                        xd->plane[0].dst.stride,
                                        xd->plane[1].dst.stride,
                                        BLOCK_SIZE_MB16X16);
      }
#endif
    }
  }

  if (!x->skip) {
#ifdef ENC_DEBUG
    if (enc_debug) {
      int i, j;
      printf("\n");
      printf("qcoeff\n");
      for (i = 0; i < 384; i++) {
        printf("%3d ", xd->qcoeff[i]);
        if (i % 16 == 15) printf("\n");
      }
      printf("\n");
      printf("predictor\n");
      for (i = 0; i < 384; i++) {
        printf("%3d ", xd->predictor[i]);
        if (i % 16 == 15) printf("\n");
      }
      printf("\n");
      printf("src_diff\n");
      for (i = 0; i < 384; i++) {
        printf("%3d ", x->src_diff[i]);
        if (i % 16 == 15) printf("\n");
      }
      printf("\n");
      printf("diff\n");
      for (i = 0; i < 384; i++) {
        printf("%3d ", xd->block[0].diff[i]);
        if (i % 16 == 15) printf("\n");
      }
      printf("\n");
      printf("final y\n");
      for (i = 0; i < 16; i++) {
        for (j = 0; j < 16; j++)
          printf("%3d ", xd->plane[0].dst.buf[i * xd->plane[0].dst.stride + j]);
        printf("\n");
      }
      printf("\n");
      printf("final u\n");
      for (i = 0; i < 8; i++) {
        for (j = 0; j < 8; j++)
          printf("%3d ", xd->plane[1].dst.buf[i * xd->plane[1].dst.stride + j]);
        printf("\n");
      }
      printf("\n");
      printf("final v\n");
      for (i = 0; i < 8; i++) {
        for (j = 0; j < 8; j++)
          printf("%3d ", xd->plane[2].dst.buf[i * xd->plane[1].dst.stride + j]);
        printf("\n");
      }
      fflush(stdout);
    }
#endif

    vp9_tokenize_sb(cpi, xd, t, !output_enabled, BLOCK_SIZE_MB16X16);
  } else {
    // FIXME(rbultje): not tile-aware (mi - 1)
    int mb_skip_context =
      (mi - 1)->mbmi.mb_skip_coeff + (mi - mis)->mbmi.mb_skip_coeff;

    mbmi->mb_skip_coeff = 1;
    if (output_enabled)
      cpi->skip_true_count[mb_skip_context]++;
    vp9_reset_sb_tokens_context(xd, BLOCK_SIZE_MB16X16);
  }

#if CONFIG_SB8X8
  // copy skip flag on all mb_mode_info contexts in this SB
  // if this was a skip at this txfm size
  for (n = 1; n < 4; n++) {
    const int x_idx = n & 1, y_idx = n >> 1;
    if (mi_col + x_idx < cm->mi_cols && mi_row + y_idx < cm->mi_rows)
      mi[x_idx + y_idx * mis].mbmi.mb_skip_coeff = mi->mbmi.mb_skip_coeff;
  }
#endif

  if (output_enabled) {
    int segment_id = mbmi->segment_id;
    if (cpi->common.txfm_mode == TX_MODE_SELECT &&
        !(mbmi->mb_skip_coeff ||
          vp9_segfeature_active(&x->e_mbd, segment_id, SEG_LVL_SKIP))) {
      assert(mbmi->txfm_size <= TX_16X16);
      if (mbmi->mode != I4X4_PRED && mbmi->mode != I8X8_PRED &&
          mbmi->mode != SPLITMV) {
        cpi->txfm_count_16x16p[mbmi->txfm_size]++;
      } else if (mbmi->mode == I8X8_PRED ||
                 (mbmi->mode == SPLITMV &&
                  mbmi->partitioning != PARTITIONING_4X4)) {
        cpi->txfm_count_8x8p[mbmi->txfm_size]++;
      }
    } else {
#if CONFIG_SB8X8
      int y, x;
#endif
      if (mbmi->mode != I4X4_PRED && mbmi->mode != I8X8_PRED &&
          mbmi->mode != SPLITMV && cpi->common.txfm_mode >= ALLOW_16X16) {
        mbmi->txfm_size = TX_16X16;
      } else if (mbmi->mode != I4X4_PRED &&
                 !(mbmi->mode == SPLITMV &&
                   mbmi->partitioning == PARTITIONING_4X4) &&
                 cpi->common.txfm_mode >= ALLOW_8X8) {
        mbmi->txfm_size = TX_8X8;
      } else {
        mbmi->txfm_size = TX_4X4;
      }

#if CONFIG_SB8X8
      for (y = 0; y < 2; y++) {
        for (x = !y; x < 2; x++) {
          if (mi_col + x < cm->mi_cols && mi_row + y < cm->mi_rows) {
            mi[mis * y + x].mbmi.txfm_size = mbmi->txfm_size;
          }
        }
      }
#endif
    }
  }
}
#endif

static void encode_superblock(VP9_COMP *cpi, TOKENEXTRA **t,
                              int output_enabled, int mi_row, int mi_col,
                              BLOCK_SIZE_TYPE bsize) {
  VP9_COMMON *const cm = &cpi->common;
  MACROBLOCK *const x = &cpi->mb;
  MACROBLOCKD *const xd = &x->e_mbd;
  int n;
  MODE_INFO *mi = x->e_mbd.mode_info_context;
  unsigned int segment_id = mi->mbmi.segment_id;
  const int mis = cm->mode_info_stride;
  const int bwl = mi_width_log2(bsize);
  const int bw = 1 << bwl, bh = 1 << mi_height_log2(bsize);

  if (cm->frame_type == KEY_FRAME) {
    if (cpi->oxcf.tuning == VP8_TUNE_SSIM) {
      adjust_act_zbin(cpi, x);
      vp9_update_zbin_extra(cpi, x);
    }
  } else {
    vp9_setup_interp_filters(xd, xd->mode_info_context->mbmi.interp_filter, cm);

    if (cpi->oxcf.tuning == VP8_TUNE_SSIM) {
      // Adjust the zbin based on this MB rate.
      adjust_act_zbin(cpi, x);
    }

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
        } else if (xd->mode_info_context->mbmi.mode == SPLITMV) {
          cpi->zbin_mode_boost = SPLIT_MV_ZBIN_BOOST;
        } else {
          cpi->zbin_mode_boost = MV_ZBIN_BOOST;
        }
      } else {
        cpi->zbin_mode_boost = INTRA_ZBIN_BOOST;
      }
    }

    vp9_update_zbin_extra(cpi, x);
  }

#if CONFIG_SB8X8
  if (xd->mode_info_context->mbmi.mode == I4X4_PRED) {
    assert(bsize == BLOCK_SIZE_SB8X8 &&
           xd->mode_info_context->mbmi.txfm_size == TX_4X4);

    vp9_encode_intra4x4mby(x, bsize);
    vp9_build_intra_predictors_sbuv_s(&x->e_mbd, bsize);
    vp9_subtract_sbuv(x, bsize);
    vp9_transform_sbuv_4x4(x, bsize);
    vp9_quantize_sbuv_4x4(x, bsize);
    if (x->optimize)
      vp9_optimize_sbuv(cm, x, bsize);
    vp9_inverse_transform_sbuv_4x4(xd, bsize);
    vp9_recon_sbuv(xd, bsize);

    if (output_enabled)
      sum_intra_stats(cpi, x);
  } else
#endif
  if (xd->mode_info_context->mbmi.ref_frame == INTRA_FRAME) {
    vp9_build_intra_predictors_sby_s(&x->e_mbd, bsize);
    vp9_build_intra_predictors_sbuv_s(&x->e_mbd, bsize);
    if (output_enabled)
      sum_intra_stats(cpi, x);
  } else {
    int ref_fb_idx, second_ref_fb_idx;

    assert(cm->frame_type != KEY_FRAME);

    if (xd->mode_info_context->mbmi.ref_frame == LAST_FRAME)
      ref_fb_idx = cpi->common.ref_frame_map[cpi->lst_fb_idx];
    else if (xd->mode_info_context->mbmi.ref_frame == GOLDEN_FRAME)
      ref_fb_idx = cpi->common.ref_frame_map[cpi->gld_fb_idx];
    else
      ref_fb_idx = cpi->common.ref_frame_map[cpi->alt_fb_idx];

    if (xd->mode_info_context->mbmi.second_ref_frame > 0) {
      if (xd->mode_info_context->mbmi.second_ref_frame == LAST_FRAME)
        second_ref_fb_idx = cpi->common.ref_frame_map[cpi->lst_fb_idx];
      else if (xd->mode_info_context->mbmi.second_ref_frame == GOLDEN_FRAME)
        second_ref_fb_idx = cpi->common.ref_frame_map[cpi->gld_fb_idx];
      else
        second_ref_fb_idx = cpi->common.ref_frame_map[cpi->alt_fb_idx];
    }

    setup_pre_planes(xd,
        &cpi->common.yv12_fb[ref_fb_idx],
        xd->mode_info_context->mbmi.second_ref_frame > 0
            ? &cpi->common.yv12_fb[second_ref_fb_idx] : NULL,
        mi_row, mi_col, xd->scale_factor, xd->scale_factor_uv);

    vp9_build_inter_predictors_sb(xd, mi_row, mi_col, bsize);
  }

#if CONFIG_SB8X8
  if (xd->mode_info_context->mbmi.mode == I4X4_PRED) {
    assert(bsize == BLOCK_SIZE_SB8X8);
    vp9_tokenize_sb(cpi, &x->e_mbd, t, !output_enabled, bsize);
  } else
#endif
  if (!x->skip) {
    vp9_encode_sb(cm, x, bsize);
    vp9_tokenize_sb(cpi, &x->e_mbd, t, !output_enabled, bsize);
  } else {
    // FIXME(rbultje): not tile-aware (mi - 1)
    int mb_skip_context =
        (mi - 1)->mbmi.mb_skip_coeff + (mi - mis)->mbmi.mb_skip_coeff;

    xd->mode_info_context->mbmi.mb_skip_coeff = 1;
    if (output_enabled)
      cpi->skip_true_count[mb_skip_context]++;
    vp9_reset_sb_tokens_context(xd, bsize);
  }

  // copy skip flag on all mb_mode_info contexts in this SB
  // if this was a skip at this txfm size
  for (n = 1; n < bw * bh; n++) {
    const int x_idx = n & (bw - 1), y_idx = n >> bwl;
    if (mi_col + x_idx < cm->mi_cols && mi_row + y_idx < cm->mi_rows)
      mi[x_idx + y_idx * mis].mbmi.mb_skip_coeff = mi->mbmi.mb_skip_coeff;
  }

  if (output_enabled) {
    if (cm->txfm_mode == TX_MODE_SELECT &&
        !(mi->mbmi.mb_skip_coeff ||
          vp9_segfeature_active(xd, segment_id, SEG_LVL_SKIP))) {
      if (bsize >= BLOCK_SIZE_SB32X32) {
        cpi->txfm_count_32x32p[mi->mbmi.txfm_size]++;
      } else if (bsize >= BLOCK_SIZE_MB16X16) {
        cpi->txfm_count_16x16p[mi->mbmi.txfm_size]++;
      } else {
        cpi->txfm_count_8x8p[mi->mbmi.txfm_size]++;
      }
    } else {
      int x, y;
      TX_SIZE sz = (cm->txfm_mode == TX_MODE_SELECT) ? TX_32X32 : cm->txfm_mode;

      if (sz == TX_32X32 && bsize < BLOCK_SIZE_SB32X32)
        sz = TX_16X16;
      if (sz == TX_16X16 && bsize < BLOCK_SIZE_MB16X16)
        sz = TX_8X8;
      if (sz == TX_8X8 && (xd->mode_info_context->mbmi.mode == SPLITMV ||
                           xd->mode_info_context->mbmi.mode == I4X4_PRED))
        sz = TX_4X4;

      for (y = 0; y < bh; y++) {
        for (x = 0; x < bw; x++) {
          if (mi_col + x < cm->mi_cols && mi_row + y < cm->mi_rows) {
            mi[mis * y + x].mbmi.txfm_size = sz;
          }
        }
      }
    }
  }
}
