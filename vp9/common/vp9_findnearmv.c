/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */


#include "vp9/common/vp9_findnearmv.h"
#include "vp9/common/vp9_sadmxn.h"
#include "vp9/common/vp9_subpelvar.h"
#include <limits.h>

const unsigned char vp9_mbsplit_offset[4][16] = {
  { 0,  8,  0,  0,  0,  0,  0,  0,  0,  0,   0,  0,  0,  0,  0,  0},
  { 0,  2,  0,  0,  0,  0,  0,  0,  0,  0,   0,  0,  0,  0,  0,  0},
  { 0,  2,  8, 10,  0,  0,  0,  0,  0,  0,   0,  0,  0,  0,  0,  0},
  { 0,  1,  2,  3,  4,  5,  6,  7,  8,  9,  10, 11, 12, 13, 14, 15}
};

static void lower_mv_precision(int_mv *mv, int usehp)
{
  if (!usehp || !vp9_use_nmv_hp(&mv->as_mv)) {
    if (mv->as_mv.row & 1)
      mv->as_mv.row += (mv->as_mv.row > 0 ? -1 : 1);
    if (mv->as_mv.col & 1)
      mv->as_mv.col += (mv->as_mv.col > 0 ? -1 : 1);
  }
}

vp9_prob *vp9_mv_ref_probs(VP9_COMMON *pc,
                           vp9_prob p[4], const int context
                          ) {
  p[0] = pc->fc.vp9_mode_contexts[context][0];
  p[1] = pc->fc.vp9_mode_contexts[context][1];
  p[2] = pc->fc.vp9_mode_contexts[context][2];
  p[3] = pc->fc.vp9_mode_contexts[context][3];
  return p;
}

#define SP(x) (((x) & 7) << 1)
unsigned int vp9_sad3x16_c(const unsigned char *src_ptr,
                           int  src_stride,
                           const unsigned char *ref_ptr,
                           int  ref_stride) {
  return sad_mx_n_c(src_ptr, src_stride, ref_ptr, ref_stride, 3, 16);
}
unsigned int vp9_sad16x3_c(const unsigned char *src_ptr,
                           int  src_stride,
                           const unsigned char *ref_ptr,
                           int  ref_stride) {
  return sad_mx_n_c(src_ptr, src_stride, ref_ptr, ref_stride, 16, 3);
}

#if CONFIG_SUBPELREFMV
unsigned int vp9_variance2x16_c(const unsigned char *src_ptr,
                                const int  source_stride,
                                const unsigned char *ref_ptr,
                                const int  recon_stride,
                                unsigned int *sse) {
  int sum;
  variance(src_ptr, source_stride, ref_ptr, recon_stride, 2, 16, sse, &sum);
  return (*sse - (((unsigned int)sum * sum) >> 5));
}

unsigned int vp9_variance16x2_c(const unsigned char *src_ptr,
                                const int  source_stride,
                                const unsigned char *ref_ptr,
                                const int  recon_stride,
                                unsigned int *sse) {
  int sum;
  variance(src_ptr, source_stride, ref_ptr, recon_stride, 16, 2, sse, &sum);
  return (*sse - (((unsigned int)sum * sum) >> 5));
}

unsigned int vp9_sub_pixel_variance16x2_c(const unsigned char  *src_ptr,
                                          const int  src_pixels_per_line,
                                          const int  xoffset,
                                          const int  yoffset,
                                          const unsigned char *dst_ptr,
                                          const int dst_pixels_per_line,
                                          unsigned int *sse) {
  unsigned short FData3[16 * 3];  // Temp data buffer used in filtering
  unsigned char  temp2[2 * 16];
  const short *HFilter, *VFilter;

  HFilter = vp9_bilinear_filters[xoffset];
  VFilter = vp9_bilinear_filters[yoffset];

  var_filter_block2d_bil_first_pass(src_ptr, FData3,
                                    src_pixels_per_line, 1, 3, 16, HFilter);
  var_filter_block2d_bil_second_pass(FData3, temp2, 16, 16, 2, 16, VFilter);

  return vp9_variance16x2_c(temp2, 16, dst_ptr, dst_pixels_per_line, sse);
}

unsigned int vp9_sub_pixel_variance2x16_c(const unsigned char  *src_ptr,
                                          const int  src_pixels_per_line,
                                          const int  xoffset,
                                          const int  yoffset,
                                          const unsigned char *dst_ptr,
                                          const int dst_pixels_per_line,
                                          unsigned int *sse) {
  unsigned short FData3[2 * 17];  // Temp data buffer used in filtering
  unsigned char  temp2[2 * 16];
  const short *HFilter, *VFilter;

  HFilter = vp9_bilinear_filters[xoffset];
  VFilter = vp9_bilinear_filters[yoffset];

  var_filter_block2d_bil_first_pass(src_ptr, FData3,
                                    src_pixels_per_line, 1, 17, 2, HFilter);
  var_filter_block2d_bil_second_pass(FData3, temp2, 2, 2, 16, 2, VFilter);

  return vp9_variance2x16_c(temp2, 2, dst_ptr, dst_pixels_per_line, sse);
}
#endif

/* check a list of motion vectors by sad score using a number rows of pixels
 * above and a number cols of pixels in the left to select the one with best
 * score to use as ref motion vector
 */
void vp9_find_best_ref_mvs(MACROBLOCKD *xd,
                           unsigned char *ref_y_buffer,
                           int ref_y_stride,
                           int_mv *mvlist,
                           int_mv *best_mv,
                           int_mv *nearest,
                           int_mv *near) {
  int i, j;
  unsigned char *above_src;
  unsigned char *left_src;
  unsigned char *above_ref;
  unsigned char *left_ref;
  unsigned int score;
#if CONFIG_SUBPELREFMV
  unsigned int sse;
#endif
  unsigned int ref_scores[MAX_MV_REFS] = {0};
  int_mv sorted_mvs[MAX_MV_REFS];
  int zero_seen = FALSE;

  // Default all to 0,0 if nothing else available
  best_mv->as_int = nearest->as_int = near->as_int = 0;
  vpx_memset(sorted_mvs, 0, sizeof(sorted_mvs));

#if CONFIG_SUBPELREFMV
  above_src = xd->dst.y_buffer - xd->dst.y_stride * 2;
  left_src  = xd->dst.y_buffer - 2;
  above_ref = ref_y_buffer - ref_y_stride * 2;
  left_ref  = ref_y_buffer - 2;
#else
  above_src = xd->dst.y_buffer - xd->dst.y_stride * 3;
  left_src  = xd->dst.y_buffer - 3;
  above_ref = ref_y_buffer - ref_y_stride * 3;
  left_ref  = ref_y_buffer - 3;
#endif

  //for(i = 0; i < MAX_MV_REFS; ++i) {
  // Limit search to the predicted best 4
  for(i = 0; i < 4; ++i) {
    int_mv this_mv;
    int offset = 0;
    int row_offset, col_offset;

    this_mv.as_int = mvlist[i].as_int;

    // If we see a 0,0 vector for a second time we have reached the end of
    // the list of valid candidate vectors.
    if (!this_mv.as_int && zero_seen)
      break;

    zero_seen = zero_seen || !this_mv.as_int;

    clamp_mv(&this_mv,
             xd->mb_to_left_edge - LEFT_TOP_MARGIN + 24,
             xd->mb_to_right_edge + RIGHT_BOTTOM_MARGIN,
             xd->mb_to_top_edge - LEFT_TOP_MARGIN + 24,
             xd->mb_to_bottom_edge + RIGHT_BOTTOM_MARGIN);

#if CONFIG_SUBPELREFMV
    row_offset = this_mv.as_mv.row >> 3;
    col_offset = this_mv.as_mv.col >> 3;
    offset = ref_y_stride * row_offset + col_offset;
    score = 0;
    if (xd->up_available) {
      vp9_sub_pixel_variance16x2_c(above_ref + offset, ref_y_stride,
                                   SP(this_mv.as_mv.col),
                                   SP(this_mv.as_mv.row),
                                   above_src, xd->dst.y_stride, &sse);
      score += sse;
#if CONFIG_SUPERBLOCKS
      if (xd->mode_info_context->mbmi.encoded_as_sb) {
        vp9_sub_pixel_variance16x2_c(above_ref + offset + 16,
                                     ref_y_stride,
                                     SP(this_mv.as_mv.col),
                                     SP(this_mv.as_mv.row),
                                     above_src + 16, xd->dst.y_stride, &sse);
        score += sse;
      }
#endif
    }
    if (xd->left_available) {
      vp9_sub_pixel_variance2x16_c(left_ref + offset, ref_y_stride,
                                   SP(this_mv.as_mv.col),
                                   SP(this_mv.as_mv.row),
                                   left_src, xd->dst.y_stride, &sse);
      score += sse;
#if CONFIG_SUPERBLOCKS
      if (xd->mode_info_context->mbmi.encoded_as_sb) {
        vp9_sub_pixel_variance2x16_c(left_ref + offset + ref_y_stride * 16,
                                     ref_y_stride,
                                     SP(this_mv.as_mv.col),
                                     SP(this_mv.as_mv.row),
                                     left_src + xd->dst.y_stride * 16,
                                     xd->dst.y_stride, &sse);
        score += sse;
      }
#endif
    }
#else
    row_offset = (this_mv.as_mv.row > 0) ?
      ((this_mv.as_mv.row + 3) >> 3):((this_mv.as_mv.row + 4) >> 3);
    col_offset = (this_mv.as_mv.col > 0) ?
      ((this_mv.as_mv.col + 3) >> 3):((this_mv.as_mv.col + 4) >> 3);
    offset = ref_y_stride * row_offset + col_offset;
    score = 0;
    if (xd->up_available) {
      score += vp9_sad16x3(above_src, xd->dst.y_stride,
                           above_ref + offset, ref_y_stride);
#if CONFIG_SUPERBLOCKS
      if (xd->mode_info_context->mbmi.encoded_as_sb) {
        score += vp9_sad16x3(above_src + 16, xd->dst.y_stride,
                             above_ref + offset + 16, ref_y_stride);
      }
#endif
    }
    if (xd->left_available) {
      score += vp9_sad3x16(left_src, xd->dst.y_stride,
                           left_ref + offset, ref_y_stride);
#if CONFIG_SUPERBLOCKS
      if (xd->mode_info_context->mbmi.encoded_as_sb) {
        score += vp9_sad3x16(left_src + xd->dst.y_stride * 16,
                             xd->dst.y_stride,
                             left_ref + offset + ref_y_stride * 16,
                             ref_y_stride);
      }
#endif
    }
#endif
    // Add the entry to our list and then resort the list on score.
    ref_scores[i] = score;
    sorted_mvs[i].as_int = this_mv.as_int;
    j = i;
    while (j > 0) {
      if (ref_scores[j] < ref_scores[j-1]) {
        ref_scores[j] = ref_scores[j-1];
        sorted_mvs[j].as_int = sorted_mvs[j-1].as_int;
        ref_scores[j-1] = score;
        sorted_mvs[j-1].as_int = this_mv.as_int;
        j--;
      } else
        break;
    }
  }

  // Make sure all the candidates are properly clamped etc
  for (i = 0; i < 4; ++i) {
    lower_mv_precision(&sorted_mvs[i], xd->allow_high_precision_mv);
    clamp_mv2(&sorted_mvs[i], xd);
  }

  // Set the best mv to the first entry in the sorted list
  best_mv->as_int = sorted_mvs[0].as_int;

  // Provided that there are non zero vectors available there will not
  // be more than one 0,0 entry in the sorted list.
  // The best ref mv is always set to the first entry (which gave the best
  // results. The nearest is set to the first non zero vector if available and
  // near to the second non zero vector if available.
  // We do not use 0,0 as a nearest or near as 0,0 has its own mode.
  if ( sorted_mvs[0].as_int ) {
    nearest->as_int = sorted_mvs[0].as_int;
    if ( sorted_mvs[1].as_int )
      near->as_int = sorted_mvs[1].as_int;
    else
      near->as_int = sorted_mvs[2].as_int;
  } else {
      nearest->as_int = sorted_mvs[1].as_int;
      near->as_int = sorted_mvs[2].as_int;
  }

  // Copy back the re-ordered mv list
  vpx_memcpy(mvlist, sorted_mvs, sizeof(sorted_mvs));
}
