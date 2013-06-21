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
#include "./vp9_rtcd.h"
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
#include "vp9/encoder/vp9_encodeintra.h"
#include "vp9/common/vp9_reconinter.h"
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

static void encode_superblock(VP9_COMP *cpi, TOKENEXTRA **t, int output_enabled,
                              int mi_row, int mi_col, BLOCK_SIZE_TYPE bsize);

static void adjust_act_zbin(VP9_COMP *cpi, MACROBLOCK *x);

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
static const uint8_t VP9_VAR_OFFS[16] = {128, 128, 128, 128, 128, 128, 128, 128,
    128, 128, 128, 128, 128, 128, 128, 128};

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
static unsigned int alt_activity_measure(VP9_COMP *cpi, MACROBLOCK *x,
                                         int use_dc_pred) {
  return vp9_encode_intra(cpi, x, use_dc_pred);
}
DECLARE_ALIGNED(16, static const uint8_t, vp9_64x64_zeros[64*64]) = {0};

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
  cpi->activity_avg = (unsigned int) (activity_sum / cpi->common.MBs);
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
  MACROBLOCK * const x = &cpi->mb;
  MACROBLOCKD *xd = &x->e_mbd;
  VP9_COMMON * const cm = &cpi->common;

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

  x->rdmult = (unsigned int) (((int64_t) x->rdmult * b + (a >> 1)) / a);
  x->errorperbit = x->rdmult * 100 / (110 * x->rddiv);
  x->errorperbit += (x->errorperbit == 0);
#endif

  // Activity based Zbin adjustment
  adjust_act_zbin(cpi, x);
}

static void update_state(VP9_COMP *cpi, PICK_MODE_CONTEXT *ctx,
                         BLOCK_SIZE_TYPE bsize, int output_enabled) {
  int i, x_idx, y;
  MACROBLOCK * const x = &cpi->mb;
  MACROBLOCKD * const xd = &x->e_mbd;
  MODE_INFO *mi = &ctx->mic;
  MB_MODE_INFO * const mbmi = &xd->mode_info_context->mbmi;
#if CONFIG_DEBUG || CONFIG_INTERNAL_STATS
  MB_PREDICTION_MODE mb_mode = mi->mbmi.mode;
#endif
  int mb_mode_index = ctx->best_mode_index;
  const int mis = cpi->common.mode_info_stride;
  const int bh = 1 << mi_height_log2(bsize), bw = 1 << mi_width_log2(bsize);

#if CONFIG_DEBUG
  assert(mb_mode < MB_MODE_COUNT);
  assert(mb_mode_index < MAX_MODES);
  assert(mi->mbmi.ref_frame[0] < MAX_REF_FRAMES);
  assert(mi->mbmi.ref_frame[1] < MAX_REF_FRAMES);
#endif

  assert(mi->mbmi.sb_type == bsize);
  // Restore the coding context of the MB to that that was in place
  // when the mode was picked for it
  for (y = 0; y < bh; y++) {
    for (x_idx = 0; x_idx < bw; x_idx++) {
      if ((xd->mb_to_right_edge >> (3 + LOG2_MI_SIZE)) + bw > x_idx
          && (xd->mb_to_bottom_edge >> (3 + LOG2_MI_SIZE)) + bh > y) {
        MODE_INFO *mi_addr = xd->mode_info_context + x_idx + y * mis;
        *mi_addr = *mi;
      }
    }
  }
  if (bsize < BLOCK_SIZE_SB32X32) {
    if (bsize < BLOCK_SIZE_MB16X16)
      ctx->txfm_rd_diff[ALLOW_16X16] = ctx->txfm_rd_diff[ALLOW_8X8];
    ctx->txfm_rd_diff[ALLOW_32X32] = ctx->txfm_rd_diff[ALLOW_16X16];
  }

  if (mbmi->ref_frame[0] != INTRA_FRAME && mbmi->sb_type < BLOCK_SIZE_SB8X8) {
    *x->partition_info = ctx->partition_info;
    mbmi->mv[0].as_int = x->partition_info->bmi[3].mv.as_int;
    mbmi->mv[1].as_int = x->partition_info->bmi[3].second_mv.as_int;
  }

  x->skip = ctx->skip;
  if (!output_enabled)
    return;

  if (!vp9_segfeature_active(xd, mbmi->segment_id, SEG_LVL_SKIP)) {
    for (i = 0; i < NB_TXFM_MODES; i++) {
      cpi->rd_tx_select_diff[i] += ctx->txfm_rd_diff[i];
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
    if (mbmi->ref_frame[0] != INTRA_FRAME
        && (mbmi->sb_type < BLOCK_SIZE_SB8X8 || mbmi->mode == NEWMV)) {
      int_mv best_mv, best_second_mv;
      const MV_REFERENCE_FRAME rf1 = mbmi->ref_frame[0];
      const MV_REFERENCE_FRAME rf2 = mbmi->ref_frame[1];
      best_mv.as_int = ctx->best_ref_mv.as_int;
      best_second_mv.as_int = ctx->second_best_ref_mv.as_int;
      if (mbmi->mode == NEWMV) {
        best_mv.as_int = mbmi->ref_mvs[rf1][0].as_int;
        best_second_mv.as_int = mbmi->ref_mvs[rf2][0].as_int;
      }
      mbmi->best_mv.as_int = best_mv.as_int;
      mbmi->best_second_mv.as_int = best_second_mv.as_int;
      vp9_update_nmv_count(cpi, x, &best_mv, &best_second_mv);
    }

    if (bsize > BLOCK_SIZE_SB8X8 && mbmi->mode == NEWMV) {
      int i, j;
      for (j = 0; j < bh; ++j)
        for (i = 0; i < bw; ++i)
          if ((xd->mb_to_right_edge >> (3 + LOG2_MI_SIZE)) + bw > i
              && (xd->mb_to_bottom_edge >> (3 + LOG2_MI_SIZE)) + bh > j)
            xd->mode_info_context[mis * j + i].mbmi = *mbmi;
    }

    if (cpi->common.mcomp_filter_type == SWITCHABLE
        && is_inter_mode(mbmi->mode)) {
      ++cpi->common.fc.switchable_interp_count[vp9_get_pred_context(
          &cpi->common, xd, PRED_SWITCHABLE_INTERP)][vp9_switchable_interp_map[mbmi
          ->interp_filter]];
    }

    cpi->rd_comp_pred_diff[SINGLE_PREDICTION_ONLY] += ctx->single_pred_diff;
    cpi->rd_comp_pred_diff[COMP_PREDICTION_ONLY] += ctx->comp_pred_diff;
    cpi->rd_comp_pred_diff[HYBRID_PREDICTION] += ctx->hybrid_pred_diff;
  }
}

static unsigned find_seg_id(VP9_COMMON *cm, uint8_t *buf, BLOCK_SIZE_TYPE bsize,
                            int start_y, int height, int start_x, int width) {
  const int bw = 1 << mi_width_log2(bsize), bh = 1 << mi_height_log2(bsize);
  const int end_x = MIN(start_x + bw, width);
  const int end_y = MIN(start_y + bh, height);
  int x, y;
  unsigned seg_id = -1;

  buf += width * start_y;
  assert(start_y < cm->mi_rows && start_x < cm->cur_tile_mi_col_end);
  for (y = start_y; y < end_y; y++, buf += width) {
    for (x = start_x; x < end_x; x++) {
      seg_id = MIN(seg_id, buf[x]);
    }
  }

  return seg_id;
}

void vp9_setup_src_planes(MACROBLOCK *x, const YV12_BUFFER_CONFIG *src,
                          int mb_row, int mb_col) {
  uint8_t *buffers[4] = {src->y_buffer, src->u_buffer, src->v_buffer, src
      ->alpha_buffer};
  int strides[4] = {src->y_stride, src->uv_stride, src->uv_stride, src
      ->alpha_stride};
  int i;

  for (i = 0; i < MAX_MB_PLANE; i++) {
    setup_pred_plane(&x->plane[i].src, buffers[i], strides[i], mb_row, mb_col,
                     NULL, x->e_mbd.plane[i].subsampling_x,
                     x->e_mbd.plane[i].subsampling_y);
  }
}

static void set_offsets(VP9_COMP *cpi, int mi_row, int mi_col,
                        BLOCK_SIZE_TYPE bsize) {
  MACROBLOCK * const x = &cpi->mb;
  VP9_COMMON * const cm = &cpi->common;
  MACROBLOCKD * const xd = &x->e_mbd;
  MB_MODE_INFO *mbmi;
  const int dst_fb_idx = cm->new_fb_idx;
  const int idx_str = xd->mode_info_stride * mi_row + mi_col;
  const int bw = 1 << mi_width_log2(bsize), bh = 1 << mi_height_log2(bsize);
  const int mb_row = mi_row >> 1;
  const int mb_col = mi_col >> 1;
  const int idx_map = mb_row * cm->mb_cols + mb_col;
  int i;

  // entropy context structures
  for (i = 0; i < MAX_MB_PLANE; i++) {
    xd->plane[i].above_context = cm->above_context[i]
        + (mi_col * 2 >> xd->plane[i].subsampling_x);
    xd->plane[i].left_context = cm->left_context[i]
        + (((mi_row * 2) & 15) >> xd->plane[i].subsampling_y);
  }

  // partition contexts
  set_partition_seg_context(cm, xd, mi_row, mi_col);

  // Activity map pointer
  x->mb_activity_ptr = &cpi->mb_activity_map[idx_map];
  x->active_ptr = cpi->active_map + idx_map;

  /* pointers to mode info contexts */
  x->partition_info = x->pi + idx_str;
  xd->mode_info_context = cm->mi + idx_str;
  mbmi = &xd->mode_info_context->mbmi;
  // Special case: if prev_mi is NULL, the previous mode info context
  // cannot be used.
  xd->prev_mode_info_context = cm->prev_mi ? cm->prev_mi + idx_str : NULL;

  // Set up destination pointers
  setup_dst_planes(xd, &cm->yv12_fb[dst_fb_idx], mi_row, mi_col);

  /* Set up limit values for MV components to prevent them from
   * extending beyond the UMV borders assuming 16x16 block size */
  x->mv_row_min = -((mi_row * MI_SIZE)+ VP9BORDERINPIXELS - VP9_INTERP_EXTEND);
  x->mv_col_min = -((mi_col * MI_SIZE)+ VP9BORDERINPIXELS - VP9_INTERP_EXTEND);
  x->mv_row_max = ((cm->mi_rows - mi_row) * MI_SIZE
      + (VP9BORDERINPIXELS - MI_SIZE * bh - VP9_INTERP_EXTEND));
  x->mv_col_max = ((cm->mi_cols - mi_col) * MI_SIZE
      + (VP9BORDERINPIXELS - MI_SIZE * bw - VP9_INTERP_EXTEND));

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
    uint8_t *map =
        xd->update_mb_segmentation_map ?
            cpi->segmentation_map : cm->last_frame_seg_map;
    mbmi->segment_id = find_seg_id(cm, map, bsize, mi_row, cm->mi_rows, mi_col,
                                   cm->mi_cols);

    assert(mbmi->segment_id <= (MAX_MB_SEGMENTS-1));
    vp9_mb_init_quantizer(cpi, x);

    if (xd->segmentation_enabled && cpi->seg0_cnt > 0
        && !vp9_segfeature_active(xd, 0, SEG_LVL_REF_FRAME)
        && vp9_segfeature_active(xd, 1, SEG_LVL_REF_FRAME)) {
      cpi->seg0_progress = (cpi->seg0_idx << 16) / cpi->seg0_cnt;
    } else {
      const int y = mb_row & ~3;
      const int x = mb_col & ~3;
      const int p16 = ((mb_row & 1) << 1) + (mb_col & 1);
      const int p32 = ((mb_row & 2) << 2) + ((mb_col & 2) << 1);
      const int tile_progress = cm->cur_tile_mi_col_start * cm->mb_rows >> 1;
      const int mb_cols = (cm->cur_tile_mi_col_end - cm->cur_tile_mi_col_start)
          >> 1;

      cpi->seg0_progress = ((y * mb_cols + x * 4 + p32 + p16 + tile_progress)
          << 16) / cm->MBs;
    }
  } else {
    mbmi->segment_id = 0;
  }
}

static void pick_sb_modes(VP9_COMP *cpi, int mi_row, int mi_col,
                          TOKENEXTRA **tp, int *totalrate, int64_t *totaldist,
                          BLOCK_SIZE_TYPE bsize, PICK_MODE_CONTEXT *ctx) {
  VP9_COMMON * const cm = &cpi->common;
  MACROBLOCK * const x = &cpi->mb;
  MACROBLOCKD * const xd = &x->e_mbd;

  x->rd_search = 1;

  if (bsize < BLOCK_SIZE_SB8X8)
    if (xd->ab_index != 0)
      return;

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
  VP9_COMMON * const cm = &cpi->common;
  MACROBLOCK * const x = &cpi->mb;
  MACROBLOCKD * const xd = &x->e_mbd;
  MODE_INFO *mi = xd->mode_info_context;
  MB_MODE_INFO * const mbmi = &mi->mbmi;

  if (cm->frame_type != KEY_FRAME) {
    int segment_id, seg_ref_active;

    segment_id = mbmi->segment_id;
    seg_ref_active = vp9_segfeature_active(xd, segment_id, SEG_LVL_REF_FRAME);

    if (!seg_ref_active)
      cpi->intra_inter_count[vp9_get_pred_context(cm, xd, PRED_INTRA_INTER)][mbmi
          ->ref_frame[0] > INTRA_FRAME]++;

    // If the segment reference feature is enabled we have only a single
    // reference frame allowed for the segment so exclude it from
    // the reference frame counts used to work out probabilities.
    if ((mbmi->ref_frame[0] > INTRA_FRAME) && !seg_ref_active) {
      if (cm->comp_pred_mode == HYBRID_PREDICTION)
        cpi->comp_inter_count[vp9_get_pred_context(cm, xd,
                                                   PRED_COMP_INTER_INTER)][mbmi
            ->ref_frame[1] > INTRA_FRAME]++;

      if (mbmi->ref_frame[1] > INTRA_FRAME) {
        cpi->comp_ref_count[vp9_get_pred_context(cm, xd, PRED_COMP_REF_P)][mbmi
            ->ref_frame[0] == GOLDEN_FRAME]++;
      } else {
        cpi->single_ref_count[vp9_get_pred_context(cm, xd, PRED_SINGLE_REF_P1)][0][mbmi
            ->ref_frame[0] != LAST_FRAME]++;
        if (mbmi->ref_frame[0] != LAST_FRAME)
          cpi->single_ref_count[vp9_get_pred_context(cm, xd, PRED_SINGLE_REF_P2)][1][mbmi
              ->ref_frame[0] != GOLDEN_FRAME]++;
      }
    }
    // Count of last ref frame 0,0 usage
    if ((mbmi->mode == ZEROMV) && (mbmi->ref_frame[0] == LAST_FRAME))
      cpi->inter_zz_count++;
  }
}

// TODO(jingning): the variables used here are little complicated. need further
// refactoring on organizing the the temporary buffers, when recursive
// partition down to 4x4 block size is enabled.
static PICK_MODE_CONTEXT *get_block_context(MACROBLOCK *x,
                                            BLOCK_SIZE_TYPE bsize) {
  MACROBLOCKD * const xd = &x->e_mbd;

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
    case BLOCK_SIZE_SB16X8:
      return &x->sb16x8_context[xd->sb_index][xd->mb_index][xd->b_index];
    case BLOCK_SIZE_SB8X16:
      return &x->sb8x16_context[xd->sb_index][xd->mb_index][xd->b_index];
    case BLOCK_SIZE_SB8X8:
      return &x->sb8x8_context[xd->sb_index][xd->mb_index][xd->b_index];
    case BLOCK_SIZE_SB8X4:
      return &x->sb8x4_context[xd->sb_index][xd->mb_index][xd->b_index];
    case BLOCK_SIZE_SB4X8:
      return &x->sb4x8_context[xd->sb_index][xd->mb_index][xd->b_index];
    case BLOCK_SIZE_AB4X4:
      return &x->ab4x4_context[xd->sb_index][xd->mb_index][xd->b_index];
    default:
      assert(0);
      return NULL ;
  }
}

static BLOCK_SIZE_TYPE *get_sb_partitioning(MACROBLOCK *x,
                                            BLOCK_SIZE_TYPE bsize) {
  MACROBLOCKD *xd = &x->e_mbd;
  switch (bsize) {
    case BLOCK_SIZE_SB64X64:
      return &x->sb64_partitioning;
    case BLOCK_SIZE_SB32X32:
      return &x->sb_partitioning[xd->sb_index];
    case BLOCK_SIZE_MB16X16:
      return &x->mb_partitioning[xd->sb_index][xd->mb_index];
    case BLOCK_SIZE_SB8X8:
      return &x->b_partitioning[xd->sb_index][xd->mb_index][xd->b_index];
    default:
      assert(0);
      return NULL ;
  }
}

static void restore_context(VP9_COMP *cpi, int mi_row, int mi_col,
                            ENTROPY_CONTEXT a[16 * MAX_MB_PLANE],
                            ENTROPY_CONTEXT l[16 * MAX_MB_PLANE],
                            PARTITION_CONTEXT sa[8], PARTITION_CONTEXT sl[8],
                            BLOCK_SIZE_TYPE bsize) {
  VP9_COMMON * const cm = &cpi->common;
  MACROBLOCK * const x = &cpi->mb;
  MACROBLOCKD * const xd = &x->e_mbd;
  int p;
  int bwl = b_width_log2(bsize), bw = 1 << bwl;
  int bhl = b_height_log2(bsize), bh = 1 << bhl;
  int mwl = mi_width_log2(bsize), mw = 1 << mwl;
  int mhl = mi_height_log2(bsize), mh = 1 << mhl;
  for (p = 0; p < MAX_MB_PLANE; p++) {
    vpx_memcpy(
        cm->above_context[p] + ((mi_col * 2) >> xd->plane[p].subsampling_x),
        a + bw * p, sizeof(ENTROPY_CONTEXT) * bw >> xd->plane[p].subsampling_x);
    vpx_memcpy(
        cm->left_context[p]
            + ((mi_row & MI_MASK)* 2 >> xd->plane[p].subsampling_y),l + bh * p,
            sizeof(ENTROPY_CONTEXT) * bh >> xd->plane[p].subsampling_y);
          }
  vpx_memcpy(cm->above_seg_context + mi_col, sa,
             sizeof(PARTITION_CONTEXT) * mw);
  vpx_memcpy(cm->left_seg_context + (mi_row & MI_MASK), sl,
  sizeof(PARTITION_CONTEXT) * mh)
             ;}
static void save_context(VP9_COMP *cpi, int mi_row, int mi_col,
                         ENTROPY_CONTEXT a[16 * MAX_MB_PLANE],
                         ENTROPY_CONTEXT l[16 * MAX_MB_PLANE],
                         PARTITION_CONTEXT sa[8], PARTITION_CONTEXT sl[8],
                         BLOCK_SIZE_TYPE bsize) {
  VP9_COMMON * const cm = &cpi->common;
  MACROBLOCK * const x = &cpi->mb;
  MACROBLOCKD * const xd = &x->e_mbd;
  int p;
  int bwl = b_width_log2(bsize), bw = 1 << bwl;
  int bhl = b_height_log2(bsize), bh = 1 << bhl;
  int mwl = mi_width_log2(bsize), mw = 1 << mwl;
  int mhl = mi_height_log2(bsize), mh = 1 << mhl;

  // buffer the above/left context information of the block in search.
  for (p = 0; p < MAX_MB_PLANE; ++p) {
    vpx_memcpy(
        a + bw * p,
        cm->above_context[p] + (mi_col * 2 >> xd->plane[p].subsampling_x),
        sizeof(ENTROPY_CONTEXT) * bw >> xd->plane[p].subsampling_x);
    vpx_memcpy(
        l + bh * p,
        cm->left_context[p]
            + ((mi_row & MI_MASK)* 2 >> xd->plane[p].subsampling_y),sizeof(ENTROPY_CONTEXT) * bh >> xd->plane[p].subsampling_y);
          }
  vpx_memcpy(sa, cm->above_seg_context + mi_col,
             sizeof(PARTITION_CONTEXT) * mw);
  vpx_memcpy(sl, cm->left_seg_context + (mi_row & MI_MASK),
  sizeof(PARTITION_CONTEXT) * mh)
             ;}

static void encode_b(VP9_COMP *cpi, TOKENEXTRA **tp, int mi_row, int mi_col,
                     int output_enabled, BLOCK_SIZE_TYPE bsize, int sub_index) {
  VP9_COMMON * const cm = &cpi->common;
  MACROBLOCK * const x = &cpi->mb;
  MACROBLOCKD * const xd = &x->e_mbd;

  if (mi_row >= cm->mi_rows || mi_col >= cm->mi_cols)
    return;

  if (sub_index != -1)
    *(get_sb_index(xd, bsize)) = sub_index;

  if (bsize < BLOCK_SIZE_SB8X8)
    if (xd->ab_index > 0)
      return;
  set_offsets(cpi, mi_row, mi_col, bsize);
  update_state(cpi, get_block_context(x, bsize), bsize, output_enabled);
  encode_superblock(cpi, tp, output_enabled, mi_row, mi_col, bsize);

  if (output_enabled) {
    update_stats(cpi, mi_row, mi_col);

    (*tp)->token = EOSB_TOKEN;
    (*tp)++;
  }
}

static void encode_sb(VP9_COMP *cpi, TOKENEXTRA **tp, int mi_row, int mi_col,
                      int output_enabled, BLOCK_SIZE_TYPE bsize) {
  VP9_COMMON * const cm = &cpi->common;
  MACROBLOCK * const x = &cpi->mb;
  MACROBLOCKD * const xd = &x->e_mbd;
  BLOCK_SIZE_TYPE c1 = BLOCK_SIZE_SB8X8;
  const int bsl = b_width_log2(bsize), bs = (1 << bsl) / 4;
  int bwl, bhl;
  int UNINITIALIZED_IS_SAFE(pl);

  if (mi_row >= cm->mi_rows || mi_col >= cm->mi_cols)
    return;

  c1 = BLOCK_SIZE_AB4X4;
  if (bsize >= BLOCK_SIZE_SB8X8) {
    set_partition_seg_context(cm, xd, mi_row, mi_col);
    pl = partition_plane_context(xd, bsize);
    c1 = *(get_sb_partitioning(x, bsize));
  }

  bwl = b_width_log2(c1), bhl = b_height_log2(c1);

  if (bsl == bwl && bsl == bhl) {
    if (output_enabled && bsize >= BLOCK_SIZE_SB8X8)
      cpi->partition_count[pl][PARTITION_NONE]++;
    encode_b(cpi, tp, mi_row, mi_col, output_enabled, c1, -1);
  } else if (bsl == bhl && bsl > bwl) {
    if (output_enabled)
      cpi->partition_count[pl][PARTITION_VERT]++;
    encode_b(cpi, tp, mi_row, mi_col, output_enabled, c1, 0);
    encode_b(cpi, tp, mi_row, mi_col + bs, output_enabled, c1, 1);
  } else if (bsl == bwl && bsl > bhl) {
    if (output_enabled)
      cpi->partition_count[pl][PARTITION_HORZ]++;
    encode_b(cpi, tp, mi_row, mi_col, output_enabled, c1, 0);
    encode_b(cpi, tp, mi_row + bs, mi_col, output_enabled, c1, 1);
  } else {
    BLOCK_SIZE_TYPE subsize;
    int i;

    assert(bwl < bsl && bhl < bsl);
    subsize = get_subsize(bsize, PARTITION_SPLIT);

    if (output_enabled)
      cpi->partition_count[pl][PARTITION_SPLIT]++;

    for (i = 0; i < 4; i++) {
      const int x_idx = i & 1, y_idx = i >> 1;

      *(get_sb_index(xd, subsize)) = i;
      encode_sb(cpi, tp, mi_row + y_idx * bs, mi_col + x_idx * bs,
                output_enabled, subsize);
    }
  }

  if (bsize >= BLOCK_SIZE_SB8X8
      && (bsize == BLOCK_SIZE_SB8X8 || bsl == bwl || bsl == bhl)) {
    set_partition_seg_context(cm, xd, mi_row, mi_col);
    update_partition_context(xd, c1, bsize);
  }
}

static void set_partitioning(VP9_COMP *cpi, MODE_INFO *m,
                             BLOCK_SIZE_TYPE bsize) {
  VP9_COMMON *const cm = &cpi->common;
  const int mis = cm->mode_info_stride;
  int block_row, block_col;
  for (block_row = 0; block_row < 8; ++block_row) {
    for (block_col = 0; block_col < 8; ++block_col) {
      m[block_row * mis + block_col].mbmi.sb_type = bsize;
    }
  }
}
static void copy_partitioning(VP9_COMP *cpi, MODE_INFO *m, MODE_INFO *p) {
  VP9_COMMON *const cm = &cpi->common;
  const int mis = cm->mode_info_stride;
  int block_row, block_col;
  for (block_row = 0; block_row < 8; ++block_row) {
    for (block_col = 0; block_col < 8; ++block_col) {
      m[block_row * mis + block_col].mbmi.sb_type =
          p[block_row * mis + block_col].mbmi.sb_type;
    }
  }
}

static void set_block_size(VP9_COMMON * const cm, MODE_INFO *m,
                           BLOCK_SIZE_TYPE bsize, int mis, int mi_row,
                           int mi_col) {
  int row, col;
  int bwl = b_width_log2(bsize);
  int bhl = b_height_log2(bsize);
  int bsl = (bwl > bhl ? bwl : bhl);

  int bs = (1 << bsl) / 2;  //
  MODE_INFO *m2 = m + mi_row * mis + mi_col;
  for (row = 0; row < bs; row++) {
    for (col = 0; col < bs; col++) {
      if (mi_row + row >= cm->mi_rows || mi_col + col >= cm->mi_cols)
        continue;
      m2[row * mis + col].mbmi.sb_type = bsize;
    }
  }
}

typedef struct {
  int64_t sum_square_error;
  int64_t sum_error;
  int count;
  int variance;
} var;

typedef struct {
  var none;
  var horz[2];
  var vert[2];
} partition_variance;

#define VT(TYPE, BLOCKSIZE) \
  typedef struct { \
    partition_variance vt; \
    BLOCKSIZE split[4]; } TYPE;

VT(v8x8, var)
VT(v16x16, v8x8)
VT(v32x32, v16x16)
VT(v64x64, v32x32)

typedef struct {
  partition_variance *vt;
  var *split[4];
} vt_node;

typedef enum {
  V16X16,
  V32X32,
  V64X64,
} TREE_LEVEL;

static void tree_to_node(void *data, BLOCK_SIZE_TYPE block_size, vt_node *node) {
  int i;
  switch (block_size) {
    case BLOCK_SIZE_SB64X64: {
      v64x64 *vt = (v64x64 *) data;
      node->vt = &vt->vt;
      for (i = 0; i < 4; i++)
        node->split[i] = &vt->split[i].vt.none;
      break;
    }
    case BLOCK_SIZE_SB32X32: {
      v32x32 *vt = (v32x32 *) data;
      node->vt = &vt->vt;
      for (i = 0; i < 4; i++)
        node->split[i] = &vt->split[i].vt.none;
      break;
    }
    case BLOCK_SIZE_MB16X16: {
      v16x16 *vt = (v16x16 *) data;
      node->vt = &vt->vt;
      for (i = 0; i < 4; i++)
        node->split[i] = &vt->split[i].vt.none;
      break;
    }
    case BLOCK_SIZE_SB8X8: {
      v8x8 *vt = (v8x8 *) data;
      node->vt = &vt->vt;
      for (i = 0; i < 4; i++)
        node->split[i] = &vt->split[i];
      break;
    }
    default:
      node->vt = 0;
      for (i = 0; i < 4; i++)
        node->split[i] = 0;
      assert(-1);
  }
}

// Set variance values given sum square error, sum error, count.
static void fill_variance(var *v, int64_t s2, int64_t s, int c) {
  v->sum_square_error = s2;
  v->sum_error = s;
  v->count = c;
  if (c > 0)
    v->variance = 256
        * (v->sum_square_error - v->sum_error * v->sum_error / v->count)
        / v->count;
  else
    v->variance = 0;
}

// Combine 2 variance structures by summing the sum_error, sum_square_error,
// and counts and then calculating the new variance.
void sum_2_variances(var *r, var *a, var*b) {
  fill_variance(r, a->sum_square_error + b->sum_square_error,
                a->sum_error + b->sum_error, a->count + b->count);
}

static void fill_variance_tree(void *data, BLOCK_SIZE_TYPE block_size) {
  vt_node node;
  tree_to_node(data, block_size, &node);
  sum_2_variances(&node.vt->horz[0], node.split[0], node.split[1]);
  sum_2_variances(&node.vt->horz[1], node.split[2], node.split[3]);
  sum_2_variances(&node.vt->vert[0], node.split[0], node.split[2]);
  sum_2_variances(&node.vt->vert[1], node.split[1], node.split[3]);
  sum_2_variances(&node.vt->none, &node.vt->vert[0], &node.vt->vert[1]);
}

#if PERFORM_RANDOM_PARTITIONING
static int set_vt_partitioning(VP9_COMP *cpi, void *data, MODE_INFO *m,
    BLOCK_SIZE_TYPE block_size, int mi_row,
    int mi_col, int mi_size) {
  VP9_COMMON * const cm = &cpi->common;
  vt_node vt;
  const int mis = cm->mode_info_stride;
  int64_t threshold = 4 * cpi->common.base_qindex * cpi->common.base_qindex;

  tree_to_node(data, block_size, &vt);

  // split none is available only if we have more than half a block size
  // in width and height inside the visible image
  if (mi_col + mi_size < cm->mi_cols && mi_row + mi_size < cm->mi_rows &&
      (rand() & 3) < 1) {
    set_block_size(cm, m, block_size, mis, mi_row, mi_col);
    return 1;
  }

  // vertical split is available on all but the bottom border
  if (mi_row + mi_size < cm->mi_rows && vt.vt->vert[0].variance < threshold
      && (rand() & 3) < 1) {
    set_block_size(cm, m, get_subsize(block_size, PARTITION_VERT), mis, mi_row,
        mi_col);
    return 1;
  }

  // horizontal split is available on all but the right border
  if (mi_col + mi_size < cm->mi_cols && vt.vt->horz[0].variance < threshold
      && (rand() & 3) < 1) {
    set_block_size(cm, m, get_subsize(block_size, PARTITION_HORZ), mis, mi_row,
        mi_col);
    return 1;
  }

  return 0;
}

#else

static int set_vt_partitioning(VP9_COMP *cpi, void *data, MODE_INFO *m,
                               BLOCK_SIZE_TYPE block_size, int mi_row,
                               int mi_col, int mi_size) {
  VP9_COMMON * const cm = &cpi->common;
  vt_node vt;
  const int mis = cm->mode_info_stride;
  int64_t threshold = 50 * cpi->common.base_qindex;

  tree_to_node(data, block_size, &vt);

  // split none is available only if we have more than half a block size
  // in width and height inside the visible image
  if (mi_col + mi_size < cm->mi_cols && mi_row + mi_size < cm->mi_rows
      && vt.vt->none.variance < threshold) {
    set_block_size(cm, m, block_size, mis, mi_row, mi_col);
    return 1;
  }

  // vertical split is available on all but the bottom border
  if (mi_row + mi_size < cm->mi_rows && vt.vt->vert[0].variance < threshold
      && vt.vt->vert[1].variance < threshold) {
    set_block_size(cm, m, get_subsize(block_size, PARTITION_VERT), mis, mi_row,
                   mi_col);
    return 1;
  }

  // horizontal split is available on all but the right border
  if (mi_col + mi_size < cm->mi_cols && vt.vt->horz[0].variance < threshold
      && vt.vt->horz[1].variance < threshold) {
    set_block_size(cm, m, get_subsize(block_size, PARTITION_HORZ), mis, mi_row,
                   mi_col);
    return 1;
  }

  return 0;
}
#endif

static void choose_partitioning(VP9_COMP *cpi, MODE_INFO *m, int mi_row,
                                int mi_col) {
  VP9_COMMON * const cm = &cpi->common;
  MACROBLOCK *x = &cpi->mb;
  MACROBLOCKD *xd = &cpi->mb.e_mbd;
  const int mis = cm->mode_info_stride;
  // TODO(JBB): More experimentation or testing of this threshold;
  int64_t threshold = 4;
  int i, j, k;
  v64x64 vt;
  unsigned char * s;
  int sp;
  const unsigned char * d;
  int dp;
  int pixels_wide = 64, pixels_high = 64;

  vpx_memset(&vt, 0, sizeof(vt));

  set_offsets(cpi, mi_row, mi_col, BLOCK_SIZE_SB64X64);

  if (xd->mb_to_right_edge < 0)
    pixels_wide += (xd->mb_to_right_edge >> 3);

  if (xd->mb_to_bottom_edge < 0)
    pixels_high += (xd->mb_to_bottom_edge >> 3);

  s = x->plane[0].src.buf;
  sp = x->plane[0].src.stride;

  // TODO(JBB): Clearly the higher the quantizer the fewer partitions we want
  // but this needs more experimentation.
  threshold = threshold * cpi->common.base_qindex * cpi->common.base_qindex;

  d = vp9_64x64_zeros;
  dp = 64;
  if (cm->frame_type != KEY_FRAME) {
    int_mv nearest_mv, near_mv;
    YV12_BUFFER_CONFIG *ref_fb = &cm->yv12_fb[0];
    YV12_BUFFER_CONFIG *second_ref_fb = NULL;

    setup_pre_planes(xd, ref_fb, second_ref_fb, mi_row, mi_col,
                     xd->scale_factor, xd->scale_factor_uv);
    xd->mode_info_context->mbmi.ref_frame[0] = LAST_FRAME;
    xd->mode_info_context->mbmi.sb_type = BLOCK_SIZE_SB64X64;
    vp9_find_best_ref_mvs(xd, m->mbmi.ref_mvs[m->mbmi.ref_frame[0]],
                          &nearest_mv, &near_mv);

    xd->mode_info_context->mbmi.mv[0] = nearest_mv;
    vp9_build_inter_predictors_sby(xd, mi_row, mi_col, BLOCK_SIZE_SB64X64);
    d = xd->plane[0].dst.buf;
    dp = xd->plane[0].dst.stride;

  }

  // Fill in the entire tree of 8x8 variances for splits.
  for (i = 0; i < 4; i++) {
    const int x32_idx = ((i & 1) << 5);
    const int y32_idx = ((i >> 1) << 5);
    for (j = 0; j < 4; j++) {
      const int x16_idx = x32_idx + ((j & 1) << 4);
      const int y16_idx = y32_idx + ((j >> 1) << 4);
      v16x16 *vst = &vt.split[i].split[j];
      for (k = 0; k < 4; k++) {
        int x_idx = x16_idx + ((k & 1) << 3);
        int y_idx = y16_idx + ((k >> 1) << 3);
        unsigned int sse = 0;
        int sum = 0;
        if (x_idx < pixels_wide && y_idx < pixels_high)
          vp9_get_sse_sum_8x8(s + y_idx * sp + x_idx, sp,
                              d + y_idx * dp + x_idx, dp, &sse, &sum);
        fill_variance(&vst->split[k].vt.none, sse, sum, 64);
      }
    }
  }
  // Fill the rest of the variance tree by summing the split partition
  // values.
  for (i = 0; i < 4; i++) {
    for (j = 0; j < 4; j++) {
      fill_variance_tree(&vt.split[i].split[j], BLOCK_SIZE_MB16X16);
    }
    fill_variance_tree(&vt.split[i], BLOCK_SIZE_SB32X32);
  }
  fill_variance_tree(&vt, BLOCK_SIZE_SB64X64);
  // Now go through the entire structure,  splitting every block size until
  // we get to one that's got a variance lower than our threshold,  or we
  // hit 8x8.
  if (!set_vt_partitioning(cpi, &vt, m, BLOCK_SIZE_SB64X64, mi_row, mi_col,
                           4)) {
    for (i = 0; i < 4; ++i) {
      const int x32_idx = ((i & 1) << 2);
      const int y32_idx = ((i >> 1) << 2);
      if (!set_vt_partitioning(cpi, &vt.split[i], m, BLOCK_SIZE_SB32X32,
                               (mi_row + y32_idx), (mi_col + x32_idx), 2)) {
        for (j = 0; j < 4; ++j) {
          const int x16_idx = ((j & 1) << 1);
          const int y16_idx = ((j >> 1) << 1);
          if (!set_vt_partitioning(cpi, &vt.split[i].split[j], m,
                                   BLOCK_SIZE_MB16X16,
                                   (mi_row + y32_idx + y16_idx),
                                   (mi_col + x32_idx + x16_idx), 1)) {
            for (k = 0; k < 4; ++k) {
              const int x8_idx = (k & 1);
              const int y8_idx = (k >> 1);
              set_block_size(cm, m, BLOCK_SIZE_SB8X8, mis,
                             (mi_row + y32_idx + y16_idx + y8_idx),
                             (mi_col + x32_idx + x16_idx + x8_idx));
            }
          }
        }
      }
    }
  }
}
static void rd_use_partition(VP9_COMP *cpi, MODE_INFO *m, TOKENEXTRA **tp,
                             int mi_row, int mi_col, BLOCK_SIZE_TYPE bsize,
                             int *rate, int64_t *dist) {
  VP9_COMMON * const cm = &cpi->common;
  MACROBLOCK * const x = &cpi->mb;
  MACROBLOCKD *xd = &cpi->mb.e_mbd;
  const int mis = cm->mode_info_stride;
  int bwl = b_width_log2(m->mbmi.sb_type);
  int bhl = b_height_log2(m->mbmi.sb_type);
  int bsl = b_width_log2(bsize);
  int bh = (1 << bhl);
  int bs = (1 << bsl);
  int bss = (1 << bsl) / 4;
  int i, pl;
  PARTITION_TYPE partition;
  BLOCK_SIZE_TYPE subsize;
  ENTROPY_CONTEXT l[16 * MAX_MB_PLANE], a[16 * MAX_MB_PLANE];
  PARTITION_CONTEXT sl[8], sa[8];
  int r = 0;
  int64_t d = 0;

  if (mi_row >= cm->mi_rows || mi_col >= cm->mi_cols)
    return;

  // parse the partition type
  if ((bwl == bsl) && (bhl == bsl))
    partition = PARTITION_NONE;
  else if ((bwl == bsl) && (bhl < bsl))
    partition = PARTITION_HORZ;
  else if ((bwl < bsl) && (bhl == bsl))
    partition = PARTITION_VERT;
  else if ((bwl < bsl) && (bhl < bsl))
    partition = PARTITION_SPLIT;
  else
    assert(0);

  subsize = get_subsize(bsize, partition);

  if (bsize < BLOCK_SIZE_SB8X8) {
    if (xd->ab_index != 0) {
      *rate = 0;
      *dist = 0;
      return;
    }
  } else {
    *(get_sb_partitioning(x, bsize)) = subsize;
  }
  pl = partition_plane_context(xd, bsize);
  save_context(cpi, mi_row, mi_col, a, l, sa, sl, bsize);
  switch (partition) {
    case PARTITION_NONE:
      pick_sb_modes(cpi, mi_row, mi_col, tp, &r, &d, bsize,
                    get_block_context(x, bsize));
      r += x->partition_cost[pl][PARTITION_NONE];
      break;
    case PARTITION_HORZ:
      *(get_sb_index(xd, subsize)) = 0;
      pick_sb_modes(cpi, mi_row, mi_col, tp, &r, &d, subsize,
                    get_block_context(x, subsize));
      if (mi_row + (bh >> 1) <= cm->mi_rows) {
        int rt;
        int64_t dt;
        update_state(cpi, get_block_context(x, subsize), subsize, 0);
        encode_superblock(cpi, tp, 0, mi_row, mi_col, subsize);
        *(get_sb_index(xd, subsize)) = 1;
        pick_sb_modes(cpi, mi_row + (bs >> 2), mi_col, tp, &rt, &dt, subsize,
                      get_block_context(x, subsize));
        r += rt;
        d += dt;
      }
      set_partition_seg_context(cm, xd, mi_row, mi_col);
      pl = partition_plane_context(xd, bsize);
      r += x->partition_cost[pl][PARTITION_HORZ];
      break;
    case PARTITION_VERT:
      *(get_sb_index(xd, subsize)) = 0;
      pick_sb_modes(cpi, mi_row, mi_col, tp, &r, &d, subsize,
                    get_block_context(x, subsize));
      if (mi_col + (bs >> 1) <= cm->mi_cols) {
        int rt;
        int64_t dt;
        update_state(cpi, get_block_context(x, subsize), subsize, 0);
        encode_superblock(cpi, tp, 0, mi_row, mi_col, subsize);
        *(get_sb_index(xd, subsize)) = 1;
        pick_sb_modes(cpi, mi_row, mi_col + (bs >> 2), tp, &rt, &dt, subsize,
                      get_block_context(x, subsize));
        r += rt;
        d += dt;
      }
      set_partition_seg_context(cm, xd, mi_row, mi_col);
      pl = partition_plane_context(xd, bsize);
      r += x->partition_cost[pl][PARTITION_VERT];
      restore_context(cpi, mi_row, mi_col, a, l, sa, sl, bsize);
      break;
    case PARTITION_SPLIT:
      for (i = 0; i < 4; i++) {
        int x_idx = (i & 1) * (bs >> 2);
        int y_idx = (i >> 1) * (bs >> 2);
        int jj = i >> 1, ii = i & 0x01;
        int rt;
        int64_t dt;

        if ((mi_row + y_idx >= cm->mi_rows) || (mi_col + x_idx >= cm->mi_cols))
          continue;

        *(get_sb_index(xd, subsize)) = i;

        rd_use_partition(cpi, m + jj * bss * mis + ii * bss, tp, mi_row + y_idx,
                         mi_col + x_idx, subsize, &rt, &dt);
        r += rt;
        d += dt;
      }
      set_partition_seg_context(cm, xd, mi_row, mi_col);
      pl = partition_plane_context(xd, bsize);
      r += x->partition_cost[pl][PARTITION_SPLIT];
      break;
    default:
      assert(0);
  }

  restore_context(cpi, mi_row, mi_col, a, l, sa, sl, bsize);

  if (r < INT_MAX && d < INT_MAX)
    encode_sb(cpi, tp, mi_row, mi_col, bsize == BLOCK_SIZE_SB64X64, bsize);
  *rate = r;
  *dist = d;
}


// TODO(jingning,jimbankoski,rbultje): properly skip partition types that are
// unlikely to be selected depending on previously rate-distortion optimization
// results, for encoding speed-up.
static void rd_pick_partition(VP9_COMP *cpi, TOKENEXTRA **tp, int mi_row,
                              int mi_col, BLOCK_SIZE_TYPE bsize, int *rate,
                              int64_t *dist) {
  VP9_COMMON * const cm = &cpi->common;
  MACROBLOCK * const x = &cpi->mb;
  MACROBLOCKD * const xd = &x->e_mbd;
  int bsl = b_width_log2(bsize), bs = 1 << bsl;
  int ms = bs / 2;
  ENTROPY_CONTEXT l[16 * MAX_MB_PLANE], a[16 * MAX_MB_PLANE];
  PARTITION_CONTEXT sl[8], sa[8];
  TOKENEXTRA *tp_orig = *tp;
  int i, pl;
  BLOCK_SIZE_TYPE subsize;
  int srate = INT_MAX;
  int64_t sdist = INT_MAX;

  if (bsize < BLOCK_SIZE_SB8X8)
    if (xd->ab_index != 0) {
      *rate = 0;
      *dist = 0;
      return;
    }
  assert(mi_height_log2(bsize) == mi_width_log2(bsize));

  save_context(cpi, mi_row, mi_col, a, l, sa, sl, bsize);

  // PARTITION_SPLIT
  if (!cpi->sf.use_partitions_greater_than
      || (cpi->sf.use_partitions_greater_than
          && bsize > cpi->sf.greater_than_block_size)) {
    if (bsize >= BLOCK_SIZE_SB8X8) {
      int r4 = 0;
      int64_t d4 = 0;
      subsize = get_subsize(bsize, PARTITION_SPLIT);
      *(get_sb_partitioning(x, bsize)) = subsize;

      for (i = 0; i < 4; ++i) {
        int x_idx = (i & 1) * (ms >> 1);
        int y_idx = (i >> 1) * (ms >> 1);
        int r = 0;
        int64_t d = 0;

        if ((mi_row + y_idx >= cm->mi_rows) || (mi_col + x_idx >= cm->mi_cols))
          continue;

        *(get_sb_index(xd, subsize)) = i;
        rd_pick_partition(cpi, tp, mi_row + y_idx, mi_col + x_idx, subsize, &r,
                          &d);

        r4 += r;
        d4 += d;
      }
      set_partition_seg_context(cm, xd, mi_row, mi_col);
      pl = partition_plane_context(xd, bsize);
      if (r4 < INT_MAX)
        r4 += x->partition_cost[pl][PARTITION_SPLIT];
      assert(r4 >= 0);
      assert(d4 >= 0);
      srate = r4;
      sdist = d4;
      restore_context(cpi, mi_row, mi_col, a, l, sa, sl, bsize);
    }
  }
  if (!cpi->sf.use_partitions_less_than
      || (cpi->sf.use_partitions_less_than
          && bsize <= cpi->sf.less_than_block_size)) {
    // PARTITION_HORZ
    if (bsize >= BLOCK_SIZE_SB8X8 && mi_col + (ms >> 1) < cm->mi_cols) {
      int r2, r = 0;
      int64_t d2, d = 0;
      subsize = get_subsize(bsize, PARTITION_HORZ);
      *(get_sb_index(xd, subsize)) = 0;
      pick_sb_modes(cpi, mi_row, mi_col, tp, &r2, &d2, subsize,
                    get_block_context(x, subsize));

      if (mi_row + (ms >> 1) < cm->mi_rows) {
        update_state(cpi, get_block_context(x, subsize), subsize, 0);
        encode_superblock(cpi, tp, 0, mi_row, mi_col, subsize);

        *(get_sb_index(xd, subsize)) = 1;
        pick_sb_modes(cpi, mi_row + (ms >> 1), mi_col, tp, &r, &d, subsize,
                      get_block_context(x, subsize));
        r2 += r;
        d2 += d;
      }
      set_partition_seg_context(cm, xd, mi_row, mi_col);
      pl = partition_plane_context(xd, bsize);
      if (r2 < INT_MAX)
        r2 += x->partition_cost[pl][PARTITION_HORZ];
      if (RDCOST(x->rdmult, x->rddiv, r2, d2)
          < RDCOST(x->rdmult, x->rddiv, srate, sdist)) {
        srate = r2;
        sdist = d2;
        *(get_sb_partitioning(x, bsize)) = subsize;
      }
      restore_context(cpi, mi_row, mi_col, a, l, sa, sl, bsize);
    }

    // PARTITION_VERT
    if (bsize >= BLOCK_SIZE_SB8X8 && mi_row + (ms >> 1) < cm->mi_rows) {
      int r2;
      int64_t d2;
      subsize = get_subsize(bsize, PARTITION_VERT);
      *(get_sb_index(xd, subsize)) = 0;
      pick_sb_modes(cpi, mi_row, mi_col, tp, &r2, &d2, subsize,
                    get_block_context(x, subsize));
      if (mi_col + (ms >> 1) < cm->mi_cols) {
        int r = 0;
        int64_t d = 0;
        update_state(cpi, get_block_context(x, subsize), subsize, 0);
        encode_superblock(cpi, tp, 0, mi_row, mi_col, subsize);

        *(get_sb_index(xd, subsize)) = 1;
        pick_sb_modes(cpi, mi_row, mi_col + (ms >> 1), tp, &r, &d, subsize,
                      get_block_context(x, subsize));
        r2 += r;
        d2 += d;
      }
      set_partition_seg_context(cm, xd, mi_row, mi_col);
      pl = partition_plane_context(xd, bsize);
      if (r2 < INT_MAX)
        r2 += x->partition_cost[pl][PARTITION_VERT];
      if (RDCOST(x->rdmult, x->rddiv, r2, d2)
          < RDCOST(x->rdmult, x->rddiv, srate, sdist)) {
        srate = r2;
        sdist = d2;
        *(get_sb_partitioning(x, bsize)) = subsize;
      }
      restore_context(cpi, mi_row, mi_col, a, l, sa, sl, bsize);
    }

    // PARTITION_NONE
    if ((mi_row + (ms >> 1) < cm->mi_rows) &&
        (mi_col + (ms >> 1) < cm->mi_cols)) {
      int r;
      int64_t d;
      pick_sb_modes(cpi, mi_row, mi_col, tp, &r, &d, bsize,
                    get_block_context(x, bsize));
      if (bsize >= BLOCK_SIZE_SB8X8) {
        set_partition_seg_context(cm, xd, mi_row, mi_col);
        pl = partition_plane_context(xd, bsize);
        r += x->partition_cost[pl][PARTITION_NONE];
      }

      if (RDCOST(x->rdmult, x->rddiv, r, d)
          < RDCOST(x->rdmult, x->rddiv, srate, sdist)) {
        srate = r;
        sdist = d;
        if (bsize >= BLOCK_SIZE_SB8X8)
          *(get_sb_partitioning(x, bsize)) = bsize;
      }
    }
  }
  *rate = srate;
  *dist = sdist;

  restore_context(cpi, mi_row, mi_col, a, l, sa, sl, bsize);

  if (srate < INT_MAX && sdist < INT_MAX)
    encode_sb(cpi, tp, mi_row, mi_col, bsize == BLOCK_SIZE_SB64X64, bsize);

  if (bsize == BLOCK_SIZE_SB64X64) {
    assert(tp_orig < *tp);
    assert(srate < INT_MAX);
    assert(sdist < INT_MAX);
  } else {
    assert(tp_orig == *tp);
  }
}

static void encode_sb_row(VP9_COMP *cpi, int mi_row, TOKENEXTRA **tp,
                          int *totalrate) {
  VP9_COMMON * const cm = &cpi->common;
  int mi_col;

  // Initialize the left context for the new SB row
  vpx_memset(&cm->left_context, 0, sizeof(cm->left_context));
  vpx_memset(cm->left_seg_context, 0, sizeof(cm->left_seg_context));

  // Code each SB in the row
  for (mi_col = cm->cur_tile_mi_col_start; mi_col < cm->cur_tile_mi_col_end;
      mi_col += 64 / MI_SIZE) {
    int dummy_rate;
    int64_t dummy_dist;
    if (cpi->sf.partition_by_variance || cpi->sf.use_lastframe_partitioning ||
        cpi->sf.use_one_partition_size_always ) {
      const int idx_str = cm->mode_info_stride * mi_row + mi_col;
      MODE_INFO *m = cm->mi + idx_str;
      MODE_INFO *p = cm->prev_mi + idx_str;

      if (cpi->sf.use_one_partition_size_always) {
        set_offsets(cpi, mi_row, mi_col, BLOCK_SIZE_SB64X64);
        set_partitioning(cpi, m, cpi->sf.always_this_block_size);
        rd_use_partition(cpi, m, tp, mi_row, mi_col, BLOCK_SIZE_SB64X64,
                         &dummy_rate, &dummy_dist);
      } else if (cpi->sf.partition_by_variance) {
        choose_partitioning(cpi, cm->mi, mi_row, mi_col);
        rd_use_partition(cpi, m, tp, mi_row, mi_col, BLOCK_SIZE_SB64X64,
                         &dummy_rate, &dummy_dist);
      } else {
        if ((cpi->common.current_video_frame & 1) == 0 || cm->prev_mi == 0
            || cpi->common.show_frame == 0
            || cpi->common.frame_type == KEY_FRAME
            || cpi->is_src_frame_alt_ref) {
          rd_pick_partition(cpi, tp, mi_row, mi_col, BLOCK_SIZE_SB64X64,
                            &dummy_rate, &dummy_dist);
        } else {
          copy_partitioning(cpi, m, p);
          rd_use_partition(cpi, m, tp, mi_row, mi_col, BLOCK_SIZE_SB64X64,
                           &dummy_rate, &dummy_dist);
        }
      }
    } else {
      rd_pick_partition(cpi, tp, mi_row, mi_col, BLOCK_SIZE_SB64X64,
                        &dummy_rate, &dummy_dist);
    }
  }
}

static void init_encode_frame_mb_context(VP9_COMP *cpi) {
  MACROBLOCK * const x = &cpi->mb;
  VP9_COMMON * const cm = &cpi->common;
  MACROBLOCKD * const xd = &x->e_mbd;

  x->act_zbin_adj = 0;
  cpi->seg0_idx = 0;

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
                   0, 0, NULL, NULL );
  setup_dst_planes(xd, &cm->yv12_fb[cm->new_fb_idx], 0, 0);

  vp9_build_block_offsets(x);

  vp9_setup_block_dptrs(&x->e_mbd, cm->subsampling_x, cm->subsampling_y);

  xd->mode_info_context->mbmi.mode = DC_PRED;
  xd->mode_info_context->mbmi.uv_mode = DC_PRED;

  vp9_zero(cpi->y_mode_count)
  vp9_zero(cpi->y_uv_mode_count)
  vp9_zero(cm->fc.inter_mode_counts)
  vp9_zero(cpi->partition_count);
  vp9_zero(cpi->intra_inter_count);
  vp9_zero(cpi->comp_inter_count);
  vp9_zero(cpi->single_ref_count);
  vp9_zero(cpi->comp_ref_count);
  vp9_zero(cm->fc.tx_count_32x32p);
  vp9_zero(cm->fc.tx_count_16x16p);
  vp9_zero(cm->fc.tx_count_8x8p);
  vp9_zero(cm->fc.mbskip_count);

  // Note: this memset assumes above_context[0], [1] and [2]
  // are allocated as part of the same buffer.
  vpx_memset(
      cm->above_context[0], 0,
      sizeof(ENTROPY_CONTEXT) * 2 * MAX_MB_PLANE * mi_cols_aligned_to_sb(cm));
  vpx_memset(cm->above_seg_context, 0,
             sizeof(PARTITION_CONTEXT) * mi_cols_aligned_to_sb(cm));
}

static void switch_lossless_mode(VP9_COMP *cpi, int lossless) {
  if (lossless) {
    cpi->mb.fwd_txm8x4 = vp9_short_walsh8x4;
    cpi->mb.fwd_txm4x4 = vp9_short_walsh4x4;
    cpi->mb.e_mbd.inv_txm4x4_1_add = vp9_short_iwalsh4x4_1_add;
    cpi->mb.e_mbd.inv_txm4x4_add = vp9_short_iwalsh4x4_add;
    cpi->mb.optimize = 0;
    cpi->common.filter_level = 0;
    cpi->zbin_mode_boost_enabled = 0;
    cpi->common.txfm_mode = ONLY_4X4;
  } else {
    cpi->mb.fwd_txm8x4 = vp9_short_fdct8x4;
    cpi->mb.fwd_txm4x4 = vp9_short_fdct4x4;
    cpi->mb.e_mbd.inv_txm4x4_1_add = vp9_short_idct4x4_1_add;
    cpi->mb.e_mbd.inv_txm4x4_add = vp9_short_idct4x4_add;
  }
}

static void encode_frame_internal(VP9_COMP *cpi) {
  int mi_row;
  MACROBLOCK * const x = &cpi->mb;
  VP9_COMMON * const cm = &cpi->common;
  MACROBLOCKD * const xd = &x->e_mbd;
  int totalrate;

//  fprintf(stderr, "encode_frame_internal frame %d (%d) type %d\n",
//           cpi->common.current_video_frame, cpi->common.show_frame,
//           cm->frame_type);

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

  vp9_zero(cm->fc.switchable_interp_count);
  vp9_zero(cpi->best_switchable_interp_count);

  xd->mode_info_context = cm->mi;
  xd->prev_mode_info_context = cm->prev_mi;

  vp9_zero(cpi->NMVcount);
  vp9_zero(cpi->coef_counts);
  vp9_zero(cm->fc.eob_branch_counts);

  cpi->mb.e_mbd.lossless = cm->base_qindex == 0 && cm->y_dc_delta_q == 0
      && cm->uv_dc_delta_q == 0 && cm->uv_ac_delta_q == 0;
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
  vpx_memset(cpi->rd_tx_select_diff, 0, sizeof(cpi->rd_tx_select_diff));
  vpx_memset(cpi->rd_tx_select_threshes, 0, sizeof(cpi->rd_tx_select_threshes));

  set_prev_mi(cm);

  {
    struct vpx_usec_timer emr_timer;
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
              mi_row < cm->cur_tile_mi_row_end; mi_row += 8)
            encode_sb_row(cpi, mi_row, &tp, &totalrate);

          cpi->tok_count[tile_row][tile_col] = (unsigned int)(tp - tp_old);
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
    return 0;
  } else {
    return (!!(ref_flags & VP9_GOLD_FLAG) + !!(ref_flags & VP9_LAST_FLAG)
        + !!(ref_flags & VP9_ALT_FLAG)) >= 2;
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

static void reset_skip_txfm_size_b(VP9_COMP *cpi, MODE_INFO *mi, int mis,
                                   TX_SIZE txfm_max, int bw, int bh, int mi_row,
                                   int mi_col, BLOCK_SIZE_TYPE bsize) {
  VP9_COMMON * const cm = &cpi->common;
  MB_MODE_INFO * const mbmi = &mi->mbmi;

  if (mi_row >= cm->mi_rows || mi_col >= cm->mi_cols)
    return;

  if (mbmi->txfm_size > txfm_max) {
    MACROBLOCK * const x = &cpi->mb;
    MACROBLOCKD * const xd = &x->e_mbd;
    const int segment_id = mbmi->segment_id;
    const int ymbs = MIN(bh, cm->mi_rows - mi_row);
    const int xmbs = MIN(bw, cm->mi_cols - mi_col);

    xd->mode_info_context = mi;
    assert(
        vp9_segfeature_active(xd, segment_id, SEG_LVL_SKIP) || get_skip_flag(mi, mis, ymbs, xmbs));
    set_txfm_flag(mi, mis, ymbs, xmbs, txfm_max);
  }
}

static void reset_skip_txfm_size_sb(VP9_COMP *cpi, MODE_INFO *mi,
                                    TX_SIZE txfm_max, int mi_row, int mi_col,
                                    BLOCK_SIZE_TYPE bsize) {
  VP9_COMMON * const cm = &cpi->common;
  const int mis = cm->mode_info_stride;
  int bwl, bhl;
  const int bsl = mi_width_log2(bsize), bs = 1 << (bsl - 1);

  if (mi_row >= cm->mi_rows || mi_col >= cm->mi_cols)
    return;

  bwl = mi_width_log2(mi->mbmi.sb_type);
  bhl = mi_height_log2(mi->mbmi.sb_type);

  if (bwl == bsl && bhl == bsl) {
    reset_skip_txfm_size_b(cpi, mi, mis, txfm_max, 1 << bsl, 1 << bsl, mi_row,
                           mi_col, bsize);
  } else if (bwl == bsl && bhl < bsl) {
    reset_skip_txfm_size_b(cpi, mi, mis, txfm_max, 1 << bsl, bs, mi_row, mi_col,
                           bsize);
    reset_skip_txfm_size_b(cpi, mi + bs * mis, mis, txfm_max, 1 << bsl, bs,
                           mi_row + bs, mi_col, bsize);
  } else if (bwl < bsl && bhl == bsl) {
    reset_skip_txfm_size_b(cpi, mi, mis, txfm_max, bs, 1 << bsl, mi_row, mi_col,
                           bsize);
    reset_skip_txfm_size_b(cpi, mi + bs, mis, txfm_max, bs, 1 << bsl, mi_row,
                           mi_col + bs, bsize);
  } else {
    BLOCK_SIZE_TYPE subsize;
    int n;

    assert(bwl < bsl && bhl < bsl);
    if (bsize == BLOCK_SIZE_SB64X64) {
      subsize = BLOCK_SIZE_SB32X32;
    } else if (bsize == BLOCK_SIZE_SB32X32) {
      subsize = BLOCK_SIZE_MB16X16;
    } else {
      assert(bsize == BLOCK_SIZE_MB16X16);
      subsize = BLOCK_SIZE_SB8X8;
    }

    for (n = 0; n < 4; n++) {
      const int y_idx = n >> 1, x_idx = n & 0x01;

      reset_skip_txfm_size_sb(cpi, mi + y_idx * bs * mis + x_idx * bs, txfm_max,
                              mi_row + y_idx * bs, mi_col + x_idx * bs,
                              subsize);
    }
  }
}

static void reset_skip_txfm_size(VP9_COMP *cpi, TX_SIZE txfm_max) {
  VP9_COMMON * const cm = &cpi->common;
  int mi_row, mi_col;
  const int mis = cm->mode_info_stride;
  MODE_INFO *mi, *mi_ptr = cm->mi;

  for (mi_row = 0; mi_row < cm->mi_rows; mi_row += 8, mi_ptr += 8 * mis) {
    mi = mi_ptr;
    for (mi_col = 0; mi_col < cm->mi_cols; mi_col += 8, mi += 8) {
      reset_skip_txfm_size_sb(cpi, mi, txfm_max, mi_row, mi_col,
                              BLOCK_SIZE_SB64X64);
    }
  }
}

void vp9_encode_frame(VP9_COMP *cpi) {
  VP9_COMMON * const cm = &cpi->common;

  // In the longer term the encoder should be generalized to match the
  // decoder such that we allow compound where one of the 3 buffers has a
  // differnt sign bias and that buffer is then the fixed ref. However, this
  // requires further work in the rd loop. For now the only supported encoder
  // side behaviour is where the ALT ref buffer has oppositie sign bias to
  // the other two.
  if ((cm->ref_frame_sign_bias[ALTREF_FRAME]
      == cm->ref_frame_sign_bias[GOLDEN_FRAME])
      || (cm->ref_frame_sign_bias[ALTREF_FRAME]
          == cm->ref_frame_sign_bias[LAST_FRAME])) {
    cm->allow_comp_inter_inter = 0;
  } else {
    cm->allow_comp_inter_inter = 1;
    cm->comp_fixed_ref = ALTREF_FRAME;
    cm->comp_var_ref[0] = LAST_FRAME;
    cm->comp_var_ref[1] = GOLDEN_FRAME;
  }

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
    if (frame_type == 3 || !cm->allow_comp_inter_inter)
      pred_type = SINGLE_PREDICTION_ONLY;
    else if (cpi->rd_prediction_type_threshes[frame_type][1]
        > cpi->rd_prediction_type_threshes[frame_type][0]
        && cpi->rd_prediction_type_threshes[frame_type][1]
            > cpi->rd_prediction_type_threshes[frame_type][2]
        && check_dual_ref_flags(cpi) && cpi->static_mb_pct == 100)
      pred_type = COMP_PREDICTION_ONLY;
    else if (cpi->rd_prediction_type_threshes[frame_type][0]
        > cpi->rd_prediction_type_threshes[frame_type][2])
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
      txfm_type =
          cpi->rd_tx_select_threshes[frame_type][ALLOW_32X32]
              > cpi->rd_tx_select_threshes[frame_type][TX_MODE_SELECT] ?
              ALLOW_32X32 : TX_MODE_SELECT;
#endif
    cpi->common.txfm_mode = txfm_type;
    cpi->common.comp_pred_mode = pred_type;
    encode_frame_internal(cpi);

    for (i = 0; i < NB_PREDICTION_TYPES; ++i) {
      const int diff = (int) (cpi->rd_comp_pred_diff[i] / cpi->common.MBs);
      cpi->rd_prediction_type_threshes[frame_type][i] += diff;
      cpi->rd_prediction_type_threshes[frame_type][i] >>= 1;
    }

    for (i = 0; i < NB_TXFM_MODES; ++i) {
      int64_t pd = cpi->rd_tx_select_diff[i];
      int diff;
      if (i == TX_MODE_SELECT)
        pd -= RDCOST(cpi->mb.rdmult, cpi->mb.rddiv,
            2048 * (TX_SIZE_MAX_SB - 1), 0);
      diff = (int) (pd / cpi->common.MBs);
      cpi->rd_tx_select_threshes[frame_type][i] += diff;
      cpi->rd_tx_select_threshes[frame_type][i] /= 2;
    }

    if (cpi->common.comp_pred_mode == HYBRID_PREDICTION) {
      int single_count_zero = 0;
      int comp_count_zero = 0;

      for (i = 0; i < COMP_INTER_CONTEXTS; i++) {
        single_count_zero += cpi->comp_inter_count[i][0];
        comp_count_zero += cpi->comp_inter_count[i][1];
      }

      if (comp_count_zero == 0) {
        cpi->common.comp_pred_mode = SINGLE_PREDICTION_ONLY;
        vp9_zero(cpi->comp_inter_count);
      } else if (single_count_zero == 0) {
        cpi->common.comp_pred_mode = COMP_PREDICTION_ONLY;
        vp9_zero(cpi->comp_inter_count);
      }
    }

    if (cpi->common.txfm_mode == TX_MODE_SELECT) {
      int count4x4 = 0;
      int count8x8_lp = 0, count8x8_8x8p = 0;
      int count16x16_16x16p = 0, count16x16_lp = 0;
      int count32x32 = 0;

      for (i = 0; i < TX_SIZE_CONTEXTS; i++)
        count4x4 += cm->fc.tx_count_32x32p[i][TX_4X4];
      for (i = 0; i < TX_SIZE_CONTEXTS; i++)
        count4x4 += cm->fc.tx_count_16x16p[i][TX_4X4];
      for (i = 0; i < TX_SIZE_CONTEXTS; i++)
        count4x4 += cm->fc.tx_count_8x8p[i][TX_4X4];

      for (i = 0; i < TX_SIZE_CONTEXTS; i++)
        count8x8_lp += cm->fc.tx_count_32x32p[i][TX_8X8];
      for (i = 0; i < TX_SIZE_CONTEXTS; i++)
        count8x8_lp += cm->fc.tx_count_16x16p[i][TX_8X8];

      for (i = 0; i < TX_SIZE_CONTEXTS; i++)
        count8x8_8x8p += cm->fc.tx_count_8x8p[i][TX_8X8];

      for (i = 0; i < TX_SIZE_CONTEXTS; i++)
        count16x16_16x16p += cm->fc.tx_count_16x16p[i][TX_16X16];

      for (i = 0; i < TX_SIZE_CONTEXTS; i++)
        count16x16_lp += cm->fc.tx_count_32x32p[i][TX_16X16];

      for (i = 0; i < TX_SIZE_CONTEXTS; i++)
        count32x32 += cm->fc.tx_count_32x32p[i][TX_32X32];

      if (count4x4 == 0 && count16x16_lp == 0 && count16x16_16x16p == 0
          && count32x32 == 0) {
        cpi->common.txfm_mode = ALLOW_8X8;
        reset_skip_txfm_size(cpi, TX_8X8);
      } else if (count8x8_8x8p == 0 && count16x16_16x16p == 0
          && count8x8_lp == 0 && count16x16_lp == 0 && count32x32 == 0) {
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

  ++cpi->y_uv_mode_count[m][uvm];
  if (xd->mode_info_context->mbmi.sb_type >= BLOCK_SIZE_SB8X8) {
    const BLOCK_SIZE_TYPE bsize = xd->mode_info_context->mbmi.sb_type;
    const int bwl = b_width_log2(bsize), bhl = b_height_log2(bsize);
    const int bsl = MIN(bwl, bhl);
    ++cpi->y_mode_count[MIN(bsl, 3)][m];
  } else {
    int idx, idy;
    int bw = 1 << b_width_log2(xd->mode_info_context->mbmi.sb_type);
    int bh = 1 << b_height_log2(xd->mode_info_context->mbmi.sb_type);
    for (idy = 0; idy < 2; idy += bh) {
      for (idx = 0; idx < 2; idx += bw) {
        int m = xd->mode_info_context->bmi[idy * 2 + idx].as_mode.first;
        ++cpi->y_mode_count[0][m];
      }
    }
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
    x->act_zbin_adj = (int) (((int64_t) b + (a >> 1)) / a) - 1;
  else
    x->act_zbin_adj = 1 - (int) (((int64_t) a + (b >> 1)) / b);
#endif
}

static void encode_superblock(VP9_COMP *cpi, TOKENEXTRA **t, int output_enabled,
                              int mi_row, int mi_col, BLOCK_SIZE_TYPE bsize) {
  VP9_COMMON * const cm = &cpi->common;
  MACROBLOCK * const x = &cpi->mb;
  MACROBLOCKD * const xd = &x->e_mbd;
  int n;
  MODE_INFO *mi = xd->mode_info_context;
  MB_MODE_INFO *mbmi = &mi->mbmi;
  unsigned int segment_id = mbmi->segment_id;
  const int mis = cm->mode_info_stride;
  const int bwl = mi_width_log2(bsize);
  const int bw = 1 << bwl, bh = 1 << mi_height_log2(bsize);
  x->rd_search = 0;

  if (cm->frame_type == KEY_FRAME) {
    if (cpi->oxcf.tuning == VP8_TUNE_SSIM) {
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
      if (mbmi->ref_frame[0] != INTRA_FRAME) {
        if (mbmi->mode == ZEROMV) {
          if (mbmi->ref_frame[0] != LAST_FRAME)
            cpi->zbin_mode_boost = GF_ZEROMV_ZBIN_BOOST;
          else
            cpi->zbin_mode_boost = LF_ZEROMV_ZBIN_BOOST;
        } else if (mbmi->sb_type < BLOCK_SIZE_SB8X8) {
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

  if (mbmi->ref_frame[0] == INTRA_FRAME) {
    vp9_encode_intra_block_y(
        cm, x, (bsize < BLOCK_SIZE_SB8X8) ? BLOCK_SIZE_SB8X8 : bsize);
    vp9_encode_intra_block_uv(
        cm, x, (bsize < BLOCK_SIZE_SB8X8) ? BLOCK_SIZE_SB8X8 : bsize);
    if (output_enabled)
      sum_intra_stats(cpi, x);
  } else {
    int idx = cm->ref_frame_map[get_ref_frame_idx(cpi, mbmi->ref_frame[0])];
    YV12_BUFFER_CONFIG *ref_fb = &cm->yv12_fb[idx];
    YV12_BUFFER_CONFIG *second_ref_fb = NULL;
    if (mbmi->ref_frame[1] > 0) {
      idx = cm->ref_frame_map[get_ref_frame_idx(cpi, mbmi->ref_frame[1])];
      second_ref_fb = &cm->yv12_fb[idx];
    }

    assert(cm->frame_type != KEY_FRAME);

    setup_pre_planes(xd, ref_fb, second_ref_fb, mi_row, mi_col,
                     xd->scale_factor, xd->scale_factor_uv);

    vp9_build_inter_predictors_sb(
        xd, mi_row, mi_col,
        bsize < BLOCK_SIZE_SB8X8 ? BLOCK_SIZE_SB8X8 : bsize);
  }

  if (xd->mode_info_context->mbmi.ref_frame[0] == INTRA_FRAME) {
    vp9_tokenize_sb(cpi, xd, t, !output_enabled,
                    (bsize < BLOCK_SIZE_SB8X8) ? BLOCK_SIZE_SB8X8 : bsize);
  } else if (!x->skip) {
    vp9_encode_sb(cm, x, (bsize < BLOCK_SIZE_SB8X8) ? BLOCK_SIZE_SB8X8 : bsize);
    vp9_tokenize_sb(cpi, xd, t, !output_enabled,
                    (bsize < BLOCK_SIZE_SB8X8) ? BLOCK_SIZE_SB8X8 : bsize);
  } else {
    // FIXME(rbultje): not tile-aware (mi - 1)
    int mb_skip_context = (mi - 1)->mbmi.mb_skip_coeff
        + (mi - mis)->mbmi.mb_skip_coeff;

    mbmi->mb_skip_coeff = 1;
    if (output_enabled)
      cm->fc.mbskip_count[mb_skip_context][1]++;
    vp9_reset_sb_tokens_context(
        xd, (bsize < BLOCK_SIZE_SB8X8) ? BLOCK_SIZE_SB8X8 : bsize);
  }

  // copy skip flag on all mb_mode_info contexts in this SB
  // if this was a skip at this txfm size
  for (n = 1; n < bw * bh; n++) {
    const int x_idx = n & (bw - 1), y_idx = n >> bwl;
    if (mi_col + x_idx < cm->mi_cols && mi_row + y_idx < cm->mi_rows)
      mi[x_idx + y_idx * mis].mbmi.mb_skip_coeff = mi->mbmi.mb_skip_coeff;
  }

  if (output_enabled) {
    if (cm->txfm_mode == TX_MODE_SELECT && mbmi->sb_type >= BLOCK_SIZE_SB8X8
        && !(mbmi->ref_frame[0] != INTRA_FRAME
            && (mbmi->mb_skip_coeff
                || vp9_segfeature_active(xd, segment_id, SEG_LVL_SKIP)))) {
      const int context = vp9_get_pred_context(cm, xd, PRED_TX_SIZE);
      if (bsize >= BLOCK_SIZE_SB32X32) {
        cm->fc.tx_count_32x32p[context][mbmi->txfm_size]++;
      } else if (bsize >= BLOCK_SIZE_MB16X16) {
        cm->fc.tx_count_16x16p[context][mbmi->txfm_size]++;
      } else {
        cm->fc.tx_count_8x8p[context][mbmi->txfm_size]++;
      }
    } else {
      int x, y;
      TX_SIZE sz = (cm->txfm_mode == TX_MODE_SELECT) ? TX_32X32 : cm->txfm_mode;
      // The new intra coding scheme requires no change of transform size
      if (mi->mbmi.ref_frame[0] != INTRA_FRAME) {
        if (sz == TX_32X32 && bsize < BLOCK_SIZE_SB32X32)
          sz = TX_16X16;
        if (sz == TX_16X16 && bsize < BLOCK_SIZE_MB16X16)
          sz = TX_8X8;
        if (sz == TX_8X8 && bsize < BLOCK_SIZE_SB8X8)
          sz = TX_4X4;
      } else if (bsize >= BLOCK_SIZE_SB8X8) {
        sz = mbmi->txfm_size;
      } else {
        sz = TX_4X4;
      }

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
