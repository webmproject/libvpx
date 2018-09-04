/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include <limits.h>
#include <math.h>
#include <stdio.h>

#include "./vp9_rtcd.h"
#include "./vpx_dsp_rtcd.h"
#include "./vpx_config.h"

#include "vpx_dsp/vpx_dsp_common.h"
#include "vpx_ports/mem.h"
#include "vpx_ports/vpx_timer.h"
#include "vpx_ports/system_state.h"

#include "vp9/common/vp9_common.h"
#include "vp9/common/vp9_entropy.h"
#include "vp9/common/vp9_entropymode.h"
#include "vp9/common/vp9_idct.h"
#include "vp9/common/vp9_mvref_common.h"
#include "vp9/common/vp9_pred_common.h"
#include "vp9/common/vp9_quant_common.h"
#include "vp9/common/vp9_reconintra.h"
#include "vp9/common/vp9_reconinter.h"
#include "vp9/common/vp9_seg_common.h"
#include "vp9/common/vp9_tile_common.h"

#include "vp9/encoder/vp9_aq_360.h"
#include "vp9/encoder/vp9_aq_complexity.h"
#include "vp9/encoder/vp9_aq_cyclicrefresh.h"
#include "vp9/encoder/vp9_aq_variance.h"
#include "vp9/encoder/vp9_encodeframe.h"
#include "vp9/encoder/vp9_encodemb.h"
#include "vp9/encoder/vp9_encodemv.h"
#include "vp9/encoder/vp9_ethread.h"
#include "vp9/encoder/vp9_extend.h"
#include "vp9/encoder/vp9_multi_thread.h"
#include "vp9/encoder/vp9_pickmode.h"
#include "vp9/encoder/vp9_rd.h"
#include "vp9/encoder/vp9_rdopt.h"
#include "vp9/encoder/vp9_segmentation.h"
#include "vp9/encoder/vp9_tokenize.h"

static void encode_superblock(VP9_COMP *cpi, ThreadData *td, TOKENEXTRA **t,
                              int output_enabled, int mi_row, int mi_col,
                              BLOCK_SIZE bsize, PICK_MODE_CONTEXT *ctx);

// This is used as a reference when computing the source variance for the
//  purpose of activity masking.
// Eventually this should be replaced by custom no-reference routines,
//  which will be faster.
static const uint8_t VP9_VAR_OFFS[64] = {
  128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128,
  128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128,
  128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128,
  128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128,
  128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128
};

#if CONFIG_VP9_HIGHBITDEPTH
static const uint16_t VP9_HIGH_VAR_OFFS_8[64] = {
  128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128,
  128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128,
  128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128,
  128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128,
  128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128
};

static const uint16_t VP9_HIGH_VAR_OFFS_10[64] = {
  128 * 4, 128 * 4, 128 * 4, 128 * 4, 128 * 4, 128 * 4, 128 * 4, 128 * 4,
  128 * 4, 128 * 4, 128 * 4, 128 * 4, 128 * 4, 128 * 4, 128 * 4, 128 * 4,
  128 * 4, 128 * 4, 128 * 4, 128 * 4, 128 * 4, 128 * 4, 128 * 4, 128 * 4,
  128 * 4, 128 * 4, 128 * 4, 128 * 4, 128 * 4, 128 * 4, 128 * 4, 128 * 4,
  128 * 4, 128 * 4, 128 * 4, 128 * 4, 128 * 4, 128 * 4, 128 * 4, 128 * 4,
  128 * 4, 128 * 4, 128 * 4, 128 * 4, 128 * 4, 128 * 4, 128 * 4, 128 * 4,
  128 * 4, 128 * 4, 128 * 4, 128 * 4, 128 * 4, 128 * 4, 128 * 4, 128 * 4,
  128 * 4, 128 * 4, 128 * 4, 128 * 4, 128 * 4, 128 * 4, 128 * 4, 128 * 4
};

static const uint16_t VP9_HIGH_VAR_OFFS_12[64] = {
  128 * 16, 128 * 16, 128 * 16, 128 * 16, 128 * 16, 128 * 16, 128 * 16,
  128 * 16, 128 * 16, 128 * 16, 128 * 16, 128 * 16, 128 * 16, 128 * 16,
  128 * 16, 128 * 16, 128 * 16, 128 * 16, 128 * 16, 128 * 16, 128 * 16,
  128 * 16, 128 * 16, 128 * 16, 128 * 16, 128 * 16, 128 * 16, 128 * 16,
  128 * 16, 128 * 16, 128 * 16, 128 * 16, 128 * 16, 128 * 16, 128 * 16,
  128 * 16, 128 * 16, 128 * 16, 128 * 16, 128 * 16, 128 * 16, 128 * 16,
  128 * 16, 128 * 16, 128 * 16, 128 * 16, 128 * 16, 128 * 16, 128 * 16,
  128 * 16, 128 * 16, 128 * 16, 128 * 16, 128 * 16, 128 * 16, 128 * 16,
  128 * 16, 128 * 16, 128 * 16, 128 * 16, 128 * 16, 128 * 16, 128 * 16,
  128 * 16
};
#endif  // CONFIG_VP9_HIGHBITDEPTH

unsigned int vp9_get_sby_variance(VP9_COMP *cpi, const struct buf_2d *ref,
                                  BLOCK_SIZE bs) {
  unsigned int sse;
  const unsigned int var =
      cpi->fn_ptr[bs].vf(ref->buf, ref->stride, VP9_VAR_OFFS, 0, &sse);
  return var;
}

#if CONFIG_VP9_HIGHBITDEPTH
unsigned int vp9_high_get_sby_variance(VP9_COMP *cpi, const struct buf_2d *ref,
                                       BLOCK_SIZE bs, int bd) {
  unsigned int var, sse;
  switch (bd) {
    case 10:
      var =
          cpi->fn_ptr[bs].vf(ref->buf, ref->stride,
                             CONVERT_TO_BYTEPTR(VP9_HIGH_VAR_OFFS_10), 0, &sse);
      break;
    case 12:
      var =
          cpi->fn_ptr[bs].vf(ref->buf, ref->stride,
                             CONVERT_TO_BYTEPTR(VP9_HIGH_VAR_OFFS_12), 0, &sse);
      break;
    case 8:
    default:
      var =
          cpi->fn_ptr[bs].vf(ref->buf, ref->stride,
                             CONVERT_TO_BYTEPTR(VP9_HIGH_VAR_OFFS_8), 0, &sse);
      break;
  }
  return var;
}
#endif  // CONFIG_VP9_HIGHBITDEPTH

unsigned int vp9_get_sby_perpixel_variance(VP9_COMP *cpi,
                                           const struct buf_2d *ref,
                                           BLOCK_SIZE bs) {
  return ROUND_POWER_OF_TWO(vp9_get_sby_variance(cpi, ref, bs),
                            num_pels_log2_lookup[bs]);
}

#if CONFIG_VP9_HIGHBITDEPTH
unsigned int vp9_high_get_sby_perpixel_variance(VP9_COMP *cpi,
                                                const struct buf_2d *ref,
                                                BLOCK_SIZE bs, int bd) {
  return (unsigned int)ROUND64_POWER_OF_TWO(
      (int64_t)vp9_high_get_sby_variance(cpi, ref, bs, bd),
      num_pels_log2_lookup[bs]);
}
#endif  // CONFIG_VP9_HIGHBITDEPTH

static unsigned int get_sby_perpixel_diff_variance(VP9_COMP *cpi,
                                                   const struct buf_2d *ref,
                                                   int mi_row, int mi_col,
                                                   BLOCK_SIZE bs) {
  unsigned int sse, var;
  uint8_t *last_y;
  const YV12_BUFFER_CONFIG *last = get_ref_frame_buffer(cpi, LAST_FRAME);

  assert(last != NULL);
  last_y =
      &last->y_buffer[mi_row * MI_SIZE * last->y_stride + mi_col * MI_SIZE];
  var = cpi->fn_ptr[bs].vf(ref->buf, ref->stride, last_y, last->y_stride, &sse);
  return ROUND_POWER_OF_TWO(var, num_pels_log2_lookup[bs]);
}

static BLOCK_SIZE get_rd_var_based_fixed_partition(VP9_COMP *cpi, MACROBLOCK *x,
                                                   int mi_row, int mi_col) {
  unsigned int var = get_sby_perpixel_diff_variance(
      cpi, &x->plane[0].src, mi_row, mi_col, BLOCK_64X64);
  if (var < 8)
    return BLOCK_64X64;
  else if (var < 128)
    return BLOCK_32X32;
  else if (var < 2048)
    return BLOCK_16X16;
  else
    return BLOCK_8X8;
}

// Lighter version of set_offsets that only sets the mode info
// pointers.
static INLINE void set_mode_info_offsets(VP9_COMMON *const cm,
                                         MACROBLOCK *const x,
                                         MACROBLOCKD *const xd, int mi_row,
                                         int mi_col) {
  const int idx_str = xd->mi_stride * mi_row + mi_col;
  xd->mi = cm->mi_grid_visible + idx_str;
  xd->mi[0] = cm->mi + idx_str;
  x->mbmi_ext = x->mbmi_ext_base + (mi_row * cm->mi_cols + mi_col);
}

static void set_offsets(VP9_COMP *cpi, const TileInfo *const tile,
                        MACROBLOCK *const x, int mi_row, int mi_col,
                        BLOCK_SIZE bsize) {
  VP9_COMMON *const cm = &cpi->common;
  MACROBLOCKD *const xd = &x->e_mbd;
  MODE_INFO *mi;
  const int mi_width = num_8x8_blocks_wide_lookup[bsize];
  const int mi_height = num_8x8_blocks_high_lookup[bsize];
  const struct segmentation *const seg = &cm->seg;
  MvLimits *const mv_limits = &x->mv_limits;

  set_skip_context(xd, mi_row, mi_col);

  set_mode_info_offsets(cm, x, xd, mi_row, mi_col);

  mi = xd->mi[0];

  // Set up destination pointers.
  vp9_setup_dst_planes(xd->plane, get_frame_new_buffer(cm), mi_row, mi_col);

  // Set up limit values for MV components.
  // Mv beyond the range do not produce new/different prediction block.
  mv_limits->row_min = -(((mi_row + mi_height) * MI_SIZE) + VP9_INTERP_EXTEND);
  mv_limits->col_min = -(((mi_col + mi_width) * MI_SIZE) + VP9_INTERP_EXTEND);
  mv_limits->row_max = (cm->mi_rows - mi_row) * MI_SIZE + VP9_INTERP_EXTEND;
  mv_limits->col_max = (cm->mi_cols - mi_col) * MI_SIZE + VP9_INTERP_EXTEND;

  // Set up distance of MB to edge of frame in 1/8th pel units.
  assert(!(mi_col & (mi_width - 1)) && !(mi_row & (mi_height - 1)));
  set_mi_row_col(xd, tile, mi_row, mi_height, mi_col, mi_width, cm->mi_rows,
                 cm->mi_cols);

  // Set up source buffers.
  vp9_setup_src_planes(x, cpi->Source, mi_row, mi_col);

  // R/D setup.
  x->rddiv = cpi->rd.RDDIV;
  x->rdmult = cpi->rd.RDMULT;

  // Setup segment ID.
  if (seg->enabled) {
    if (cpi->oxcf.aq_mode != VARIANCE_AQ && cpi->oxcf.aq_mode != LOOKAHEAD_AQ &&
        cpi->oxcf.aq_mode != EQUATOR360_AQ) {
      const uint8_t *const map =
          seg->update_map ? cpi->segmentation_map : cm->last_frame_seg_map;
      mi->segment_id = get_segment_id(cm, map, bsize, mi_row, mi_col);
    }
    vp9_init_plane_quantizers(cpi, x);

    x->encode_breakout = cpi->segment_encode_breakout[mi->segment_id];
  } else {
    mi->segment_id = 0;
    x->encode_breakout = cpi->encode_breakout;
  }

  // required by vp9_append_sub8x8_mvs_for_idx() and vp9_find_best_ref_mvs()
  xd->tile = *tile;
}

static void duplicate_mode_info_in_sb(VP9_COMMON *cm, MACROBLOCKD *xd,
                                      int mi_row, int mi_col,
                                      BLOCK_SIZE bsize) {
  const int block_width =
      VPXMIN(num_8x8_blocks_wide_lookup[bsize], cm->mi_cols - mi_col);
  const int block_height =
      VPXMIN(num_8x8_blocks_high_lookup[bsize], cm->mi_rows - mi_row);
  const int mi_stride = xd->mi_stride;
  MODE_INFO *const src_mi = xd->mi[0];
  int i, j;

  for (j = 0; j < block_height; ++j)
    for (i = 0; i < block_width; ++i) xd->mi[j * mi_stride + i] = src_mi;
}

static void set_block_size(VP9_COMP *const cpi, MACROBLOCK *const x,
                           MACROBLOCKD *const xd, int mi_row, int mi_col,
                           BLOCK_SIZE bsize) {
  if (cpi->common.mi_cols > mi_col && cpi->common.mi_rows > mi_row) {
    set_mode_info_offsets(&cpi->common, x, xd, mi_row, mi_col);
    xd->mi[0]->sb_type = bsize;
  }
}

typedef struct {
  // This struct is used for computing variance in choose_partitioning(), where
  // the max number of samples within a superblock is 16x16 (with 4x4 avg). Even
  // in high bitdepth, uint32_t is enough for sum_square_error (2^12 * 2^12 * 16
  // * 16 = 2^32).
  uint32_t sum_square_error;
  int32_t sum_error;
  int log2_count;
  int variance;
} var;

typedef struct {
  var none;
  var horz[2];
  var vert[2];
} partition_variance;

typedef struct {
  partition_variance part_variances;
  var split[4];
} v4x4;

typedef struct {
  partition_variance part_variances;
  v4x4 split[4];
} v8x8;

typedef struct {
  partition_variance part_variances;
  v8x8 split[4];
} v16x16;

typedef struct {
  partition_variance part_variances;
  v16x16 split[4];
} v32x32;

typedef struct {
  partition_variance part_variances;
  v32x32 split[4];
} v64x64;

typedef struct {
  partition_variance *part_variances;
  var *split[4];
} variance_node;

typedef enum {
  V16X16,
  V32X32,
  V64X64,
} TREE_LEVEL;

static void tree_to_node(void *data, BLOCK_SIZE bsize, variance_node *node) {
  int i;
  node->part_variances = NULL;
  switch (bsize) {
    case BLOCK_64X64: {
      v64x64 *vt = (v64x64 *)data;
      node->part_variances = &vt->part_variances;
      for (i = 0; i < 4; i++)
        node->split[i] = &vt->split[i].part_variances.none;
      break;
    }
    case BLOCK_32X32: {
      v32x32 *vt = (v32x32 *)data;
      node->part_variances = &vt->part_variances;
      for (i = 0; i < 4; i++)
        node->split[i] = &vt->split[i].part_variances.none;
      break;
    }
    case BLOCK_16X16: {
      v16x16 *vt = (v16x16 *)data;
      node->part_variances = &vt->part_variances;
      for (i = 0; i < 4; i++)
        node->split[i] = &vt->split[i].part_variances.none;
      break;
    }
    case BLOCK_8X8: {
      v8x8 *vt = (v8x8 *)data;
      node->part_variances = &vt->part_variances;
      for (i = 0; i < 4; i++)
        node->split[i] = &vt->split[i].part_variances.none;
      break;
    }
    default: {
      v4x4 *vt = (v4x4 *)data;
      assert(bsize == BLOCK_4X4);
      node->part_variances = &vt->part_variances;
      for (i = 0; i < 4; i++) node->split[i] = &vt->split[i];
      break;
    }
  }
}

// Set variance values given sum square error, sum error, count.
static void fill_variance(uint32_t s2, int32_t s, int c, var *v) {
  v->sum_square_error = s2;
  v->sum_error = s;
  v->log2_count = c;
}

static void get_variance(var *v) {
  v->variance =
      (int)(256 * (v->sum_square_error -
                   (uint32_t)(((int64_t)v->sum_error * v->sum_error) >>
                              v->log2_count)) >>
            v->log2_count);
}

static void sum_2_variances(const var *a, const var *b, var *r) {
  assert(a->log2_count == b->log2_count);
  fill_variance(a->sum_square_error + b->sum_square_error,
                a->sum_error + b->sum_error, a->log2_count + 1, r);
}

static void fill_variance_tree(void *data, BLOCK_SIZE bsize) {
  variance_node node;
  memset(&node, 0, sizeof(node));
  tree_to_node(data, bsize, &node);
  sum_2_variances(node.split[0], node.split[1], &node.part_variances->horz[0]);
  sum_2_variances(node.split[2], node.split[3], &node.part_variances->horz[1]);
  sum_2_variances(node.split[0], node.split[2], &node.part_variances->vert[0]);
  sum_2_variances(node.split[1], node.split[3], &node.part_variances->vert[1]);
  sum_2_variances(&node.part_variances->vert[0], &node.part_variances->vert[1],
                  &node.part_variances->none);
}

static int set_vt_partitioning(VP9_COMP *cpi, MACROBLOCK *const x,
                               MACROBLOCKD *const xd, void *data,
                               BLOCK_SIZE bsize, int mi_row, int mi_col,
                               int64_t threshold, BLOCK_SIZE bsize_min,
                               int force_split) {
  VP9_COMMON *const cm = &cpi->common;
  variance_node vt;
  const int block_width = num_8x8_blocks_wide_lookup[bsize];
  const int block_height = num_8x8_blocks_high_lookup[bsize];

  assert(block_height == block_width);
  tree_to_node(data, bsize, &vt);

  if (force_split == 1) return 0;

  // For bsize=bsize_min (16x16/8x8 for 8x8/4x4 downsampling), select if
  // variance is below threshold, otherwise split will be selected.
  // No check for vert/horiz split as too few samples for variance.
  if (bsize == bsize_min) {
    // Variance already computed to set the force_split.
    if (frame_is_intra_only(cm)) get_variance(&vt.part_variances->none);
    if (mi_col + block_width / 2 < cm->mi_cols &&
        mi_row + block_height / 2 < cm->mi_rows &&
        vt.part_variances->none.variance < threshold) {
      set_block_size(cpi, x, xd, mi_row, mi_col, bsize);
      return 1;
    }
    return 0;
  } else if (bsize > bsize_min) {
    // Variance already computed to set the force_split.
    if (frame_is_intra_only(cm)) get_variance(&vt.part_variances->none);
    // For key frame: take split for bsize above 32X32 or very high variance.
    if (frame_is_intra_only(cm) &&
        (bsize > BLOCK_32X32 ||
         vt.part_variances->none.variance > (threshold << 4))) {
      return 0;
    }
    // If variance is low, take the bsize (no split).
    if (mi_col + block_width / 2 < cm->mi_cols &&
        mi_row + block_height / 2 < cm->mi_rows &&
        vt.part_variances->none.variance < threshold) {
      set_block_size(cpi, x, xd, mi_row, mi_col, bsize);
      return 1;
    }

    // Check vertical split.
    if (mi_row + block_height / 2 < cm->mi_rows) {
      BLOCK_SIZE subsize = get_subsize(bsize, PARTITION_VERT);
      get_variance(&vt.part_variances->vert[0]);
      get_variance(&vt.part_variances->vert[1]);
      if (vt.part_variances->vert[0].variance < threshold &&
          vt.part_variances->vert[1].variance < threshold &&
          get_plane_block_size(subsize, &xd->plane[1]) < BLOCK_INVALID) {
        set_block_size(cpi, x, xd, mi_row, mi_col, subsize);
        set_block_size(cpi, x, xd, mi_row, mi_col + block_width / 2, subsize);
        return 1;
      }
    }
    // Check horizontal split.
    if (mi_col + block_width / 2 < cm->mi_cols) {
      BLOCK_SIZE subsize = get_subsize(bsize, PARTITION_HORZ);
      get_variance(&vt.part_variances->horz[0]);
      get_variance(&vt.part_variances->horz[1]);
      if (vt.part_variances->horz[0].variance < threshold &&
          vt.part_variances->horz[1].variance < threshold &&
          get_plane_block_size(subsize, &xd->plane[1]) < BLOCK_INVALID) {
        set_block_size(cpi, x, xd, mi_row, mi_col, subsize);
        set_block_size(cpi, x, xd, mi_row + block_height / 2, mi_col, subsize);
        return 1;
      }
    }

    return 0;
  }
  return 0;
}

static int64_t scale_part_thresh_sumdiff(int64_t threshold_base, int speed,
                                         int width, int height,
                                         int content_state) {
  if (speed >= 8) {
    if (width <= 640 && height <= 480)
      return (5 * threshold_base) >> 2;
    else if ((content_state == kLowSadLowSumdiff) ||
             (content_state == kHighSadLowSumdiff) ||
             (content_state == kLowVarHighSumdiff))
      return (5 * threshold_base) >> 2;
  } else if (speed == 7) {
    if ((content_state == kLowSadLowSumdiff) ||
        (content_state == kHighSadLowSumdiff) ||
        (content_state == kLowVarHighSumdiff)) {
      return (5 * threshold_base) >> 2;
    }
  }
  return threshold_base;
}

// Set the variance split thresholds for following the block sizes:
// 0 - threshold_64x64, 1 - threshold_32x32, 2 - threshold_16x16,
// 3 - vbp_threshold_8x8. vbp_threshold_8x8 (to split to 4x4 partition) is
// currently only used on key frame.
static void set_vbp_thresholds(VP9_COMP *cpi, int64_t thresholds[], int q,
                               int content_state) {
  VP9_COMMON *const cm = &cpi->common;
  const int is_key_frame = frame_is_intra_only(cm);
  const int threshold_multiplier = is_key_frame ? 20 : 1;
  int64_t threshold_base =
      (int64_t)(threshold_multiplier * cpi->y_dequant[q][1]);

  if (is_key_frame) {
    thresholds[0] = threshold_base;
    thresholds[1] = threshold_base >> 2;
    thresholds[2] = threshold_base >> 2;
    thresholds[3] = threshold_base << 2;
  } else {
    // Increase base variance threshold based on estimated noise level.
    if (cpi->noise_estimate.enabled && cm->width >= 640 && cm->height >= 480) {
      NOISE_LEVEL noise_level =
          vp9_noise_estimate_extract_level(&cpi->noise_estimate);
      if (noise_level == kHigh)
        threshold_base = 3 * threshold_base;
      else if (noise_level == kMedium)
        threshold_base = threshold_base << 1;
      else if (noise_level < kLow)
        threshold_base = (7 * threshold_base) >> 3;
    }
#if CONFIG_VP9_TEMPORAL_DENOISING
    if (cpi->oxcf.noise_sensitivity > 0 && denoise_svc(cpi) &&
        cpi->oxcf.speed > 5 && cpi->denoiser.denoising_level >= kDenLow)
      threshold_base =
          vp9_scale_part_thresh(threshold_base, cpi->denoiser.denoising_level,
                                content_state, cpi->svc.temporal_layer_id);
    else
      threshold_base =
          scale_part_thresh_sumdiff(threshold_base, cpi->oxcf.speed, cm->width,
                                    cm->height, content_state);
#else
    // Increase base variance threshold based on content_state/sum_diff level.
    threshold_base = scale_part_thresh_sumdiff(
        threshold_base, cpi->oxcf.speed, cm->width, cm->height, content_state);
#endif
    thresholds[0] = threshold_base;
    thresholds[2] = threshold_base << cpi->oxcf.speed;
    if (cm->width >= 1280 && cm->height >= 720 && cpi->oxcf.speed < 7)
      thresholds[2] = thresholds[2] << 1;
    if (cm->width <= 352 && cm->height <= 288) {
      thresholds[0] = threshold_base >> 3;
      thresholds[1] = threshold_base >> 1;
      thresholds[2] = threshold_base << 3;
    } else if (cm->width < 1280 && cm->height < 720) {
      thresholds[1] = (5 * threshold_base) >> 2;
    } else if (cm->width < 1920 && cm->height < 1080) {
      thresholds[1] = threshold_base << 1;
    } else {
      thresholds[1] = (5 * threshold_base) >> 1;
    }
    if (cpi->sf.disable_16x16part_nonkey) thresholds[2] = INT64_MAX;
  }
}

void vp9_set_variance_partition_thresholds(VP9_COMP *cpi, int q,
                                           int content_state) {
  VP9_COMMON *const cm = &cpi->common;
  SPEED_FEATURES *const sf = &cpi->sf;
  const int is_key_frame = frame_is_intra_only(cm);
  if (sf->partition_search_type != VAR_BASED_PARTITION &&
      sf->partition_search_type != REFERENCE_PARTITION) {
    return;
  } else {
    set_vbp_thresholds(cpi, cpi->vbp_thresholds, q, content_state);
    // The thresholds below are not changed locally.
    if (is_key_frame) {
      cpi->vbp_threshold_sad = 0;
      cpi->vbp_threshold_copy = 0;
      cpi->vbp_bsize_min = BLOCK_8X8;
    } else {
      if (cm->width <= 352 && cm->height <= 288)
        cpi->vbp_threshold_sad = 10;
      else
        cpi->vbp_threshold_sad = (cpi->y_dequant[q][1] << 1) > 1000
                                     ? (cpi->y_dequant[q][1] << 1)
                                     : 1000;
      cpi->vbp_bsize_min = BLOCK_16X16;
      if (cm->width <= 352 && cm->height <= 288)
        cpi->vbp_threshold_copy = 4000;
      else if (cm->width <= 640 && cm->height <= 360)
        cpi->vbp_threshold_copy = 8000;
      else
        cpi->vbp_threshold_copy = (cpi->y_dequant[q][1] << 3) > 8000
                                      ? (cpi->y_dequant[q][1] << 3)
                                      : 8000;
      if (cpi->rc.high_source_sad ||
          (cpi->use_svc && cpi->svc.high_source_sad_superframe)) {
        cpi->vbp_threshold_sad = 0;
        cpi->vbp_threshold_copy = 0;
      }
    }
    cpi->vbp_threshold_minmax = 15 + (q >> 3);
  }
}

// Compute the minmax over the 8x8 subblocks.
static int compute_minmax_8x8(const uint8_t *s, int sp, const uint8_t *d,
                              int dp, int x16_idx, int y16_idx,
#if CONFIG_VP9_HIGHBITDEPTH
                              int highbd_flag,
#endif
                              int pixels_wide, int pixels_high) {
  int k;
  int minmax_max = 0;
  int minmax_min = 255;
  // Loop over the 4 8x8 subblocks.
  for (k = 0; k < 4; k++) {
    int x8_idx = x16_idx + ((k & 1) << 3);
    int y8_idx = y16_idx + ((k >> 1) << 3);
    int min = 0;
    int max = 0;
    if (x8_idx < pixels_wide && y8_idx < pixels_high) {
#if CONFIG_VP9_HIGHBITDEPTH
      if (highbd_flag & YV12_FLAG_HIGHBITDEPTH) {
        vpx_highbd_minmax_8x8(s + y8_idx * sp + x8_idx, sp,
                              d + y8_idx * dp + x8_idx, dp, &min, &max);
      } else {
        vpx_minmax_8x8(s + y8_idx * sp + x8_idx, sp, d + y8_idx * dp + x8_idx,
                       dp, &min, &max);
      }
#else
      vpx_minmax_8x8(s + y8_idx * sp + x8_idx, sp, d + y8_idx * dp + x8_idx, dp,
                     &min, &max);
#endif
      if ((max - min) > minmax_max) minmax_max = (max - min);
      if ((max - min) < minmax_min) minmax_min = (max - min);
    }
  }
  return (minmax_max - minmax_min);
}

static void fill_variance_4x4avg(const uint8_t *s, int sp, const uint8_t *d,
                                 int dp, int x8_idx, int y8_idx, v8x8 *vst,
#if CONFIG_VP9_HIGHBITDEPTH
                                 int highbd_flag,
#endif
                                 int pixels_wide, int pixels_high,
                                 int is_key_frame) {
  int k;
  for (k = 0; k < 4; k++) {
    int x4_idx = x8_idx + ((k & 1) << 2);
    int y4_idx = y8_idx + ((k >> 1) << 2);
    unsigned int sse = 0;
    int sum = 0;
    if (x4_idx < pixels_wide && y4_idx < pixels_high) {
      int s_avg;
      int d_avg = 128;
#if CONFIG_VP9_HIGHBITDEPTH
      if (highbd_flag & YV12_FLAG_HIGHBITDEPTH) {
        s_avg = vpx_highbd_avg_4x4(s + y4_idx * sp + x4_idx, sp);
        if (!is_key_frame)
          d_avg = vpx_highbd_avg_4x4(d + y4_idx * dp + x4_idx, dp);
      } else {
        s_avg = vpx_avg_4x4(s + y4_idx * sp + x4_idx, sp);
        if (!is_key_frame) d_avg = vpx_avg_4x4(d + y4_idx * dp + x4_idx, dp);
      }
#else
      s_avg = vpx_avg_4x4(s + y4_idx * sp + x4_idx, sp);
      if (!is_key_frame) d_avg = vpx_avg_4x4(d + y4_idx * dp + x4_idx, dp);
#endif
      sum = s_avg - d_avg;
      sse = sum * sum;
    }
    fill_variance(sse, sum, 0, &vst->split[k].part_variances.none);
  }
}

static void fill_variance_8x8avg(const uint8_t *s, int sp, const uint8_t *d,
                                 int dp, int x16_idx, int y16_idx, v16x16 *vst,
#if CONFIG_VP9_HIGHBITDEPTH
                                 int highbd_flag,
#endif
                                 int pixels_wide, int pixels_high,
                                 int is_key_frame) {
  int k;
  for (k = 0; k < 4; k++) {
    int x8_idx = x16_idx + ((k & 1) << 3);
    int y8_idx = y16_idx + ((k >> 1) << 3);
    unsigned int sse = 0;
    int sum = 0;
    if (x8_idx < pixels_wide && y8_idx < pixels_high) {
      int s_avg;
      int d_avg = 128;
#if CONFIG_VP9_HIGHBITDEPTH
      if (highbd_flag & YV12_FLAG_HIGHBITDEPTH) {
        s_avg = vpx_highbd_avg_8x8(s + y8_idx * sp + x8_idx, sp);
        if (!is_key_frame)
          d_avg = vpx_highbd_avg_8x8(d + y8_idx * dp + x8_idx, dp);
      } else {
        s_avg = vpx_avg_8x8(s + y8_idx * sp + x8_idx, sp);
        if (!is_key_frame) d_avg = vpx_avg_8x8(d + y8_idx * dp + x8_idx, dp);
      }
#else
      s_avg = vpx_avg_8x8(s + y8_idx * sp + x8_idx, sp);
      if (!is_key_frame) d_avg = vpx_avg_8x8(d + y8_idx * dp + x8_idx, dp);
#endif
      sum = s_avg - d_avg;
      sse = sum * sum;
    }
    fill_variance(sse, sum, 0, &vst->split[k].part_variances.none);
  }
}

// Check if most of the superblock is skin content, and if so, force split to
// 32x32, and set x->sb_is_skin for use in mode selection.
static int skin_sb_split(VP9_COMP *cpi, MACROBLOCK *x, const int low_res,
                         int mi_row, int mi_col, int *force_split) {
  VP9_COMMON *const cm = &cpi->common;
#if CONFIG_VP9_HIGHBITDEPTH
  if (cm->use_highbitdepth) return 0;
#endif
  // Avoid checking superblocks on/near boundary and avoid low resolutions.
  // Note superblock may still pick 64X64 if y_sad is very small
  // (i.e., y_sad < cpi->vbp_threshold_sad) below. For now leave this as is.
  if (!low_res && (mi_col >= 8 && mi_col + 8 < cm->mi_cols && mi_row >= 8 &&
                   mi_row + 8 < cm->mi_rows)) {
    int num_16x16_skin = 0;
    int num_16x16_nonskin = 0;
    uint8_t *ysignal = x->plane[0].src.buf;
    uint8_t *usignal = x->plane[1].src.buf;
    uint8_t *vsignal = x->plane[2].src.buf;
    int sp = x->plane[0].src.stride;
    int spuv = x->plane[1].src.stride;
    const int block_index = mi_row * cm->mi_cols + mi_col;
    const int bw = num_8x8_blocks_wide_lookup[BLOCK_64X64];
    const int bh = num_8x8_blocks_high_lookup[BLOCK_64X64];
    const int xmis = VPXMIN(cm->mi_cols - mi_col, bw);
    const int ymis = VPXMIN(cm->mi_rows - mi_row, bh);
    // Loop through the 16x16 sub-blocks.
    int i, j;
    for (i = 0; i < ymis; i += 2) {
      for (j = 0; j < xmis; j += 2) {
        int bl_index = block_index + i * cm->mi_cols + j;
        int is_skin = cpi->skin_map[bl_index];
        num_16x16_skin += is_skin;
        num_16x16_nonskin += (1 - is_skin);
        if (num_16x16_nonskin > 3) {
          // Exit loop if at least 4 of the 16x16 blocks are not skin.
          i = ymis;
          break;
        }
        ysignal += 16;
        usignal += 8;
        vsignal += 8;
      }
      ysignal += (sp << 4) - 64;
      usignal += (spuv << 3) - 32;
      vsignal += (spuv << 3) - 32;
    }
    if (num_16x16_skin > 12) {
      *force_split = 1;
      return 1;
    }
  }
  return 0;
}

static void set_low_temp_var_flag(VP9_COMP *cpi, MACROBLOCK *x, MACROBLOCKD *xd,
                                  v64x64 *vt, int64_t thresholds[],
                                  MV_REFERENCE_FRAME ref_frame_partition,
                                  int mi_col, int mi_row) {
  int i, j;
  VP9_COMMON *const cm = &cpi->common;
  const int mv_thr = cm->width > 640 ? 8 : 4;
  // Check temporal variance for bsize >= 16x16, if LAST_FRAME was selected and
  // int_pro mv is small. If the temporal variance is small set the flag
  // variance_low for the block. The variance threshold can be adjusted, the
  // higher the more aggressive.
  if (ref_frame_partition == LAST_FRAME &&
      (cpi->sf.short_circuit_low_temp_var == 1 ||
       (xd->mi[0]->mv[0].as_mv.col < mv_thr &&
        xd->mi[0]->mv[0].as_mv.col > -mv_thr &&
        xd->mi[0]->mv[0].as_mv.row < mv_thr &&
        xd->mi[0]->mv[0].as_mv.row > -mv_thr))) {
    if (xd->mi[0]->sb_type == BLOCK_64X64) {
      if ((vt->part_variances).none.variance < (thresholds[0] >> 1))
        x->variance_low[0] = 1;
    } else if (xd->mi[0]->sb_type == BLOCK_64X32) {
      for (i = 0; i < 2; i++) {
        if (vt->part_variances.horz[i].variance < (thresholds[0] >> 2))
          x->variance_low[i + 1] = 1;
      }
    } else if (xd->mi[0]->sb_type == BLOCK_32X64) {
      for (i = 0; i < 2; i++) {
        if (vt->part_variances.vert[i].variance < (thresholds[0] >> 2))
          x->variance_low[i + 3] = 1;
      }
    } else {
      for (i = 0; i < 4; i++) {
        const int idx[4][2] = { { 0, 0 }, { 0, 4 }, { 4, 0 }, { 4, 4 } };
        const int idx_str =
            cm->mi_stride * (mi_row + idx[i][0]) + mi_col + idx[i][1];
        MODE_INFO **this_mi = cm->mi_grid_visible + idx_str;

        if (cm->mi_cols <= mi_col + idx[i][1] ||
            cm->mi_rows <= mi_row + idx[i][0])
          continue;

        if ((*this_mi)->sb_type == BLOCK_32X32) {
          int64_t threshold_32x32 = (cpi->sf.short_circuit_low_temp_var == 1 ||
                                     cpi->sf.short_circuit_low_temp_var == 3)
                                        ? ((5 * thresholds[1]) >> 3)
                                        : (thresholds[1] >> 1);
          if (vt->split[i].part_variances.none.variance < threshold_32x32)
            x->variance_low[i + 5] = 1;
        } else if (cpi->sf.short_circuit_low_temp_var >= 2) {
          // For 32x16 and 16x32 blocks, the flag is set on each 16x16 block
          // inside.
          if ((*this_mi)->sb_type == BLOCK_16X16 ||
              (*this_mi)->sb_type == BLOCK_32X16 ||
              (*this_mi)->sb_type == BLOCK_16X32) {
            for (j = 0; j < 4; j++) {
              if (vt->split[i].split[j].part_variances.none.variance <
                  (thresholds[2] >> 8))
                x->variance_low[(i << 2) + j + 9] = 1;
            }
          }
        }
      }
    }
  }
}

static void copy_partitioning_helper(VP9_COMP *cpi, MACROBLOCK *x,
                                     MACROBLOCKD *xd, BLOCK_SIZE bsize,
                                     int mi_row, int mi_col) {
  VP9_COMMON *const cm = &cpi->common;
  BLOCK_SIZE *prev_part = cpi->prev_partition;
  int start_pos = mi_row * cm->mi_stride + mi_col;

  const int bsl = b_width_log2_lookup[bsize];
  const int bs = (1 << bsl) >> 2;
  BLOCK_SIZE subsize;
  PARTITION_TYPE partition;

  if (mi_row >= cm->mi_rows || mi_col >= cm->mi_cols) return;

  partition = partition_lookup[bsl][prev_part[start_pos]];
  subsize = get_subsize(bsize, partition);

  if (subsize < BLOCK_8X8) {
    set_block_size(cpi, x, xd, mi_row, mi_col, bsize);
  } else {
    switch (partition) {
      case PARTITION_NONE:
        set_block_size(cpi, x, xd, mi_row, mi_col, bsize);
        break;
      case PARTITION_HORZ:
        set_block_size(cpi, x, xd, mi_row, mi_col, subsize);
        set_block_size(cpi, x, xd, mi_row + bs, mi_col, subsize);
        break;
      case PARTITION_VERT:
        set_block_size(cpi, x, xd, mi_row, mi_col, subsize);
        set_block_size(cpi, x, xd, mi_row, mi_col + bs, subsize);
        break;
      default:
        assert(partition == PARTITION_SPLIT);
        copy_partitioning_helper(cpi, x, xd, subsize, mi_row, mi_col);
        copy_partitioning_helper(cpi, x, xd, subsize, mi_row + bs, mi_col);
        copy_partitioning_helper(cpi, x, xd, subsize, mi_row, mi_col + bs);
        copy_partitioning_helper(cpi, x, xd, subsize, mi_row + bs, mi_col + bs);
        break;
    }
  }
}

static int copy_partitioning(VP9_COMP *cpi, MACROBLOCK *x, MACROBLOCKD *xd,
                             int mi_row, int mi_col, int segment_id,
                             int sb_offset) {
  int svc_copy_allowed = 1;
  int frames_since_key_thresh = 1;
  if (cpi->use_svc) {
    // For SVC, don't allow copy if base spatial layer is key frame, or if
    // frame is not a temporal enhancement layer frame.
    int layer = LAYER_IDS_TO_IDX(0, cpi->svc.temporal_layer_id,
                                 cpi->svc.number_temporal_layers);
    const LAYER_CONTEXT *lc = &cpi->svc.layer_context[layer];
    if (lc->is_key_frame || !cpi->svc.non_reference_frame) svc_copy_allowed = 0;
    frames_since_key_thresh = cpi->svc.number_spatial_layers << 1;
  }
  if (cpi->rc.frames_since_key > frames_since_key_thresh && svc_copy_allowed &&
      !cpi->resize_pending && segment_id == CR_SEGMENT_ID_BASE &&
      cpi->prev_segment_id[sb_offset] == CR_SEGMENT_ID_BASE &&
      cpi->copied_frame_cnt[sb_offset] < cpi->max_copied_frame) {
    if (cpi->prev_partition != NULL) {
      copy_partitioning_helper(cpi, x, xd, BLOCK_64X64, mi_row, mi_col);
      cpi->copied_frame_cnt[sb_offset] += 1;
      memcpy(x->variance_low, &(cpi->prev_variance_low[sb_offset * 25]),
             sizeof(x->variance_low));
      return 1;
    }
  }

  return 0;
}

static int scale_partitioning_svc(VP9_COMP *cpi, MACROBLOCK *x, MACROBLOCKD *xd,
                                  BLOCK_SIZE bsize, int mi_row, int mi_col,
                                  int mi_row_high, int mi_col_high) {
  VP9_COMMON *const cm = &cpi->common;
  SVC *const svc = &cpi->svc;
  BLOCK_SIZE *prev_part = svc->prev_partition_svc;
  // Variables with _high are for higher resolution.
  int bsize_high = 0;
  int subsize_high = 0;
  const int bsl_high = b_width_log2_lookup[bsize];
  const int bs_high = (1 << bsl_high) >> 2;
  const int has_rows = (mi_row_high + bs_high) < cm->mi_rows;
  const int has_cols = (mi_col_high + bs_high) < cm->mi_cols;

  const int row_boundary_block_scale_factor[BLOCK_SIZES] = {
    13, 13, 13, 1, 0, 1, 1, 0, 1, 1, 0, 1, 0
  };
  const int col_boundary_block_scale_factor[BLOCK_SIZES] = {
    13, 13, 13, 2, 2, 0, 2, 2, 0, 2, 2, 0, 0
  };
  int start_pos;
  BLOCK_SIZE bsize_low;
  PARTITION_TYPE partition_high;

  if (mi_row_high >= cm->mi_rows || mi_col_high >= cm->mi_cols) return 0;
  if (mi_row >= (cm->mi_rows >> 1) || mi_col >= (cm->mi_cols >> 1)) return 0;

  // Find corresponding (mi_col/mi_row) block down-scaled by 2x2.
  start_pos = mi_row * (svc->mi_stride[svc->spatial_layer_id - 1]) + mi_col;
  bsize_low = prev_part[start_pos];
  // The block size is too big for boundaries. Do variance based partitioning.
  if ((!has_rows || !has_cols) && bsize_low > BLOCK_16X16) return 1;

  // For reference frames: return 1 (do variance-based partitioning) if the
  // superblock is not low source sad and lower-resoln bsize is below 32x32.
  if (!cpi->svc.non_reference_frame && !x->skip_low_source_sad &&
      bsize_low < BLOCK_32X32)
    return 1;

  // Scale up block size by 2x2. Force 64x64 for size larger than 32x32.
  if (bsize_low < BLOCK_32X32) {
    bsize_high = bsize_low + 3;
  } else if (bsize_low >= BLOCK_32X32) {
    bsize_high = BLOCK_64X64;
  }
  // Scale up blocks on boundary.
  if (!has_cols && has_rows) {
    bsize_high = bsize_low + row_boundary_block_scale_factor[bsize_low];
  } else if (has_cols && !has_rows) {
    bsize_high = bsize_low + col_boundary_block_scale_factor[bsize_low];
  } else if (!has_cols && !has_rows) {
    bsize_high = bsize_low;
  }

  partition_high = partition_lookup[bsl_high][bsize_high];
  subsize_high = get_subsize(bsize, partition_high);

  if (subsize_high < BLOCK_8X8) {
    set_block_size(cpi, x, xd, mi_row_high, mi_col_high, bsize_high);
  } else {
    const int bsl = b_width_log2_lookup[bsize];
    const int bs = (1 << bsl) >> 2;
    switch (partition_high) {
      case PARTITION_NONE:
        set_block_size(cpi, x, xd, mi_row_high, mi_col_high, bsize_high);
        break;
      case PARTITION_HORZ:
        set_block_size(cpi, x, xd, mi_row_high, mi_col_high, subsize_high);
        if (subsize_high < BLOCK_64X64)
          set_block_size(cpi, x, xd, mi_row_high + bs_high, mi_col_high,
                         subsize_high);
        break;
      case PARTITION_VERT:
        set_block_size(cpi, x, xd, mi_row_high, mi_col_high, subsize_high);
        if (subsize_high < BLOCK_64X64)
          set_block_size(cpi, x, xd, mi_row_high, mi_col_high + bs_high,
                         subsize_high);
        break;
      default:
        assert(partition_high == PARTITION_SPLIT);
        if (scale_partitioning_svc(cpi, x, xd, subsize_high, mi_row, mi_col,
                                   mi_row_high, mi_col_high))
          return 1;
        if (scale_partitioning_svc(cpi, x, xd, subsize_high, mi_row + (bs >> 1),
                                   mi_col, mi_row_high + bs_high, mi_col_high))
          return 1;
        if (scale_partitioning_svc(cpi, x, xd, subsize_high, mi_row,
                                   mi_col + (bs >> 1), mi_row_high,
                                   mi_col_high + bs_high))
          return 1;
        if (scale_partitioning_svc(cpi, x, xd, subsize_high, mi_row + (bs >> 1),
                                   mi_col + (bs >> 1), mi_row_high + bs_high,
                                   mi_col_high + bs_high))
          return 1;
        break;
    }
  }

  return 0;
}

static void update_partition_svc(VP9_COMP *cpi, BLOCK_SIZE bsize, int mi_row,
                                 int mi_col) {
  VP9_COMMON *const cm = &cpi->common;
  BLOCK_SIZE *prev_part = cpi->svc.prev_partition_svc;
  int start_pos = mi_row * cm->mi_stride + mi_col;
  const int bsl = b_width_log2_lookup[bsize];
  const int bs = (1 << bsl) >> 2;
  BLOCK_SIZE subsize;
  PARTITION_TYPE partition;
  const MODE_INFO *mi = NULL;
  int xx, yy;

  if (mi_row >= cm->mi_rows || mi_col >= cm->mi_cols) return;

  mi = cm->mi_grid_visible[start_pos];
  partition = partition_lookup[bsl][mi->sb_type];
  subsize = get_subsize(bsize, partition);
  if (subsize < BLOCK_8X8) {
    prev_part[start_pos] = bsize;
  } else {
    switch (partition) {
      case PARTITION_NONE:
        prev_part[start_pos] = bsize;
        if (bsize == BLOCK_64X64) {
          for (xx = 0; xx < 8; xx += 4)
            for (yy = 0; yy < 8; yy += 4) {
              if ((mi_row + xx < cm->mi_rows) && (mi_col + yy < cm->mi_cols))
                prev_part[start_pos + xx * cm->mi_stride + yy] = bsize;
            }
        }
        break;
      case PARTITION_HORZ:
        prev_part[start_pos] = subsize;
        if (mi_row + bs < cm->mi_rows)
          prev_part[start_pos + bs * cm->mi_stride] = subsize;
        break;
      case PARTITION_VERT:
        prev_part[start_pos] = subsize;
        if (mi_col + bs < cm->mi_cols) prev_part[start_pos + bs] = subsize;
        break;
      default:
        assert(partition == PARTITION_SPLIT);
        update_partition_svc(cpi, subsize, mi_row, mi_col);
        update_partition_svc(cpi, subsize, mi_row + bs, mi_col);
        update_partition_svc(cpi, subsize, mi_row, mi_col + bs);
        update_partition_svc(cpi, subsize, mi_row + bs, mi_col + bs);
        break;
    }
  }
}

static void update_prev_partition_helper(VP9_COMP *cpi, BLOCK_SIZE bsize,
                                         int mi_row, int mi_col) {
  VP9_COMMON *const cm = &cpi->common;
  BLOCK_SIZE *prev_part = cpi->prev_partition;
  int start_pos = mi_row * cm->mi_stride + mi_col;
  const int bsl = b_width_log2_lookup[bsize];
  const int bs = (1 << bsl) >> 2;
  BLOCK_SIZE subsize;
  PARTITION_TYPE partition;
  const MODE_INFO *mi = NULL;

  if (mi_row >= cm->mi_rows || mi_col >= cm->mi_cols) return;

  mi = cm->mi_grid_visible[start_pos];
  partition = partition_lookup[bsl][mi->sb_type];
  subsize = get_subsize(bsize, partition);
  if (subsize < BLOCK_8X8) {
    prev_part[start_pos] = bsize;
  } else {
    switch (partition) {
      case PARTITION_NONE: prev_part[start_pos] = bsize; break;
      case PARTITION_HORZ:
        prev_part[start_pos] = subsize;
        if (mi_row + bs < cm->mi_rows)
          prev_part[start_pos + bs * cm->mi_stride] = subsize;
        break;
      case PARTITION_VERT:
        prev_part[start_pos] = subsize;
        if (mi_col + bs < cm->mi_cols) prev_part[start_pos + bs] = subsize;
        break;
      default:
        assert(partition == PARTITION_SPLIT);
        update_prev_partition_helper(cpi, subsize, mi_row, mi_col);
        update_prev_partition_helper(cpi, subsize, mi_row + bs, mi_col);
        update_prev_partition_helper(cpi, subsize, mi_row, mi_col + bs);
        update_prev_partition_helper(cpi, subsize, mi_row + bs, mi_col + bs);
        break;
    }
  }
}

static void update_prev_partition(VP9_COMP *cpi, MACROBLOCK *x, int segment_id,
                                  int mi_row, int mi_col, int sb_offset) {
  update_prev_partition_helper(cpi, BLOCK_64X64, mi_row, mi_col);
  cpi->prev_segment_id[sb_offset] = segment_id;
  memcpy(&(cpi->prev_variance_low[sb_offset * 25]), x->variance_low,
         sizeof(x->variance_low));
  // Reset the counter for copy partitioning
  cpi->copied_frame_cnt[sb_offset] = 0;
}

static void chroma_check(VP9_COMP *cpi, MACROBLOCK *x, int bsize,
                         unsigned int y_sad, int is_key_frame) {
  int i;
  MACROBLOCKD *xd = &x->e_mbd;

  if (is_key_frame) return;

  // For speed >= 8, avoid the chroma check if y_sad is above threshold.
  if (cpi->oxcf.speed >= 8) {
    if (y_sad > cpi->vbp_thresholds[1] &&
        (!cpi->noise_estimate.enabled ||
         vp9_noise_estimate_extract_level(&cpi->noise_estimate) < kMedium))
      return;
  }

  for (i = 1; i <= 2; ++i) {
    unsigned int uv_sad = UINT_MAX;
    struct macroblock_plane *p = &x->plane[i];
    struct macroblockd_plane *pd = &xd->plane[i];
    const BLOCK_SIZE bs = get_plane_block_size(bsize, pd);

    if (bs != BLOCK_INVALID)
      uv_sad = cpi->fn_ptr[bs].sdf(p->src.buf, p->src.stride, pd->dst.buf,
                                   pd->dst.stride);

    // TODO(marpan): Investigate if we should lower this threshold if
    // superblock is detected as skin.
    x->color_sensitivity[i - 1] = uv_sad > (y_sad >> 2);
  }
}

static uint64_t avg_source_sad(VP9_COMP *cpi, MACROBLOCK *x, int shift,
                               int sb_offset) {
  unsigned int tmp_sse;
  uint64_t tmp_sad;
  unsigned int tmp_variance;
  const BLOCK_SIZE bsize = BLOCK_64X64;
  uint8_t *src_y = cpi->Source->y_buffer;
  int src_ystride = cpi->Source->y_stride;
  uint8_t *last_src_y = cpi->Last_Source->y_buffer;
  int last_src_ystride = cpi->Last_Source->y_stride;
  uint64_t avg_source_sad_threshold = 10000;
  uint64_t avg_source_sad_threshold2 = 12000;
#if CONFIG_VP9_HIGHBITDEPTH
  if (cpi->common.use_highbitdepth) return 0;
#endif
  src_y += shift;
  last_src_y += shift;
  tmp_sad =
      cpi->fn_ptr[bsize].sdf(src_y, src_ystride, last_src_y, last_src_ystride);
  tmp_variance = vpx_variance64x64(src_y, src_ystride, last_src_y,
                                   last_src_ystride, &tmp_sse);
  // Note: tmp_sse - tmp_variance = ((sum * sum) >> 12)
  if (tmp_sad < avg_source_sad_threshold)
    x->content_state_sb = ((tmp_sse - tmp_variance) < 25) ? kLowSadLowSumdiff
                                                          : kLowSadHighSumdiff;
  else
    x->content_state_sb = ((tmp_sse - tmp_variance) < 25) ? kHighSadLowSumdiff
                                                          : kHighSadHighSumdiff;

  // Detect large lighting change.
  if (cpi->oxcf.content != VP9E_CONTENT_SCREEN &&
      cpi->oxcf.rc_mode == VPX_CBR && tmp_variance < (tmp_sse >> 3) &&
      (tmp_sse - tmp_variance) > 10000)
    x->content_state_sb = kLowVarHighSumdiff;
  else if (tmp_sad > (avg_source_sad_threshold << 1))
    x->content_state_sb = kVeryHighSad;

  if (cpi->content_state_sb_fd != NULL) {
    if (tmp_sad < avg_source_sad_threshold2) {
      // Cap the increment to 255.
      if (cpi->content_state_sb_fd[sb_offset] < 255)
        cpi->content_state_sb_fd[sb_offset]++;
    } else {
      cpi->content_state_sb_fd[sb_offset] = 0;
    }
  }
  if (tmp_sad == 0) x->zero_temp_sad_source = 1;
  return tmp_sad;
}

// This function chooses partitioning based on the variance between source and
// reconstructed last, where variance is computed for down-sampled inputs.
static int choose_partitioning(VP9_COMP *cpi, const TileInfo *const tile,
                               MACROBLOCK *x, int mi_row, int mi_col) {
  VP9_COMMON *const cm = &cpi->common;
  MACROBLOCKD *xd = &x->e_mbd;
  int i, j, k, m;
  v64x64 vt;
  v16x16 *vt2 = NULL;
  int force_split[21];
  int avg_32x32;
  int max_var_32x32 = 0;
  int min_var_32x32 = INT_MAX;
  int var_32x32;
  int avg_16x16[4];
  int maxvar_16x16[4];
  int minvar_16x16[4];
  int64_t threshold_4x4avg;
  NOISE_LEVEL noise_level = kLow;
  int content_state = 0;
  uint8_t *s;
  const uint8_t *d;
  int sp;
  int dp;
  int compute_minmax_variance = 1;
  unsigned int y_sad = UINT_MAX;
  BLOCK_SIZE bsize = BLOCK_64X64;
  // Ref frame used in partitioning.
  MV_REFERENCE_FRAME ref_frame_partition = LAST_FRAME;
  int pixels_wide = 64, pixels_high = 64;
  int64_t thresholds[4] = { cpi->vbp_thresholds[0], cpi->vbp_thresholds[1],
                            cpi->vbp_thresholds[2], cpi->vbp_thresholds[3] };
  int scene_change_detected =
      cpi->rc.high_source_sad ||
      (cpi->use_svc && cpi->svc.high_source_sad_superframe);

  // For the variance computation under SVC mode, we treat the frame as key if
  // the reference (base layer frame) is key frame (i.e., is_key_frame == 1).
  int is_key_frame =
      (frame_is_intra_only(cm) ||
       (is_one_pass_cbr_svc(cpi) &&
        cpi->svc.layer_context[cpi->svc.temporal_layer_id].is_key_frame));
  // Always use 4x4 partition for key frame.
  const int use_4x4_partition = frame_is_intra_only(cm);
  const int low_res = (cm->width <= 352 && cm->height <= 288);
  int variance4x4downsample[16];
  int segment_id;
  int sb_offset = (cm->mi_stride >> 3) * (mi_row >> 3) + (mi_col >> 3);

  // For SVC: check if LAST frame is NULL or if the resolution of LAST is
  // different than the current frame resolution, and if so, treat this frame
  // as a key frame, for the purpose of the superblock partitioning.
  // LAST == NULL can happen in some cases where enhancement spatial layers are
  // enabled dyanmically in the stream and the only reference is the spatial
  // reference (GOLDEN).
  if (cpi->use_svc) {
    const YV12_BUFFER_CONFIG *const ref = get_ref_frame_buffer(cpi, LAST_FRAME);
    if (ref == NULL || ref->y_crop_height != cm->height ||
        ref->y_crop_width != cm->width)
      is_key_frame = 1;
  }

  set_offsets(cpi, tile, x, mi_row, mi_col, BLOCK_64X64);
  segment_id = xd->mi[0]->segment_id;

  if (cpi->oxcf.speed >= 8 || (cpi->use_svc && cpi->svc.non_reference_frame))
    compute_minmax_variance = 0;

  memset(x->variance_low, 0, sizeof(x->variance_low));

  if (cpi->sf.use_source_sad && !is_key_frame) {
    int sb_offset2 = ((cm->mi_cols + 7) >> 3) * (mi_row >> 3) + (mi_col >> 3);
    content_state = x->content_state_sb;
    x->skip_low_source_sad = (content_state == kLowSadLowSumdiff ||
                              content_state == kLowSadHighSumdiff)
                                 ? 1
                                 : 0;
    x->lowvar_highsumdiff = (content_state == kLowVarHighSumdiff) ? 1 : 0;
    if (cpi->content_state_sb_fd != NULL)
      x->last_sb_high_content = cpi->content_state_sb_fd[sb_offset2];

    // For SVC on top spatial layer: use/scale the partition from
    // the lower spatial resolution if svc_use_lowres_part is enabled.
    if (cpi->sf.svc_use_lowres_part &&
        cpi->svc.spatial_layer_id == cpi->svc.number_spatial_layers - 1 &&
        cpi->svc.prev_partition_svc != NULL && content_state != kVeryHighSad) {
      if (!scale_partitioning_svc(cpi, x, xd, BLOCK_64X64, mi_row >> 1,
                                  mi_col >> 1, mi_row, mi_col)) {
        if (cpi->sf.copy_partition_flag) {
          update_prev_partition(cpi, x, segment_id, mi_row, mi_col, sb_offset);
        }
        return 0;
      }
    }
    // If source_sad is low copy the partition without computing the y_sad.
    if (x->skip_low_source_sad && cpi->sf.copy_partition_flag &&
        !scene_change_detected &&
        copy_partitioning(cpi, x, xd, mi_row, mi_col, segment_id, sb_offset)) {
      x->sb_use_mv_part = 1;
      if (cpi->sf.svc_use_lowres_part &&
          cpi->svc.spatial_layer_id == cpi->svc.number_spatial_layers - 2)
        update_partition_svc(cpi, BLOCK_64X64, mi_row, mi_col);
      return 0;
    }
  }

  if (cpi->oxcf.aq_mode == CYCLIC_REFRESH_AQ && cm->seg.enabled &&
      cyclic_refresh_segment_id_boosted(segment_id)) {
    int q = vp9_get_qindex(&cm->seg, segment_id, cm->base_qindex);
    set_vbp_thresholds(cpi, thresholds, q, content_state);
  } else {
    set_vbp_thresholds(cpi, thresholds, cm->base_qindex, content_state);
  }

  // For non keyframes, disable 4x4 average for low resolution when speed = 8
  threshold_4x4avg = (cpi->oxcf.speed < 8) ? thresholds[1] << 1 : INT64_MAX;

  if (xd->mb_to_right_edge < 0) pixels_wide += (xd->mb_to_right_edge >> 3);
  if (xd->mb_to_bottom_edge < 0) pixels_high += (xd->mb_to_bottom_edge >> 3);

  s = x->plane[0].src.buf;
  sp = x->plane[0].src.stride;

  // Index for force_split: 0 for 64x64, 1-4 for 32x32 blocks,
  // 5-20 for the 16x16 blocks.
  force_split[0] = scene_change_detected;

  if (!is_key_frame) {
    // In the case of spatial/temporal scalable coding, the assumption here is
    // that the temporal reference frame will always be of type LAST_FRAME.
    // TODO(marpan): If that assumption is broken, we need to revisit this code.
    MODE_INFO *mi = xd->mi[0];
    YV12_BUFFER_CONFIG *yv12 = get_ref_frame_buffer(cpi, LAST_FRAME);

    const YV12_BUFFER_CONFIG *yv12_g = NULL;
    unsigned int y_sad_g, y_sad_thr, y_sad_last;
    bsize = BLOCK_32X32 + (mi_col + 4 < cm->mi_cols) * 2 +
            (mi_row + 4 < cm->mi_rows);

    assert(yv12 != NULL);

    if (!(is_one_pass_cbr_svc(cpi) && cpi->svc.spatial_layer_id) ||
        cpi->svc.use_gf_temporal_ref_current_layer) {
      // For now, GOLDEN will not be used for non-zero spatial layers, since
      // it may not be a temporal reference.
      yv12_g = get_ref_frame_buffer(cpi, GOLDEN_FRAME);
    }

    // Only compute y_sad_g (sad for golden reference) for speed < 8.
    if (cpi->oxcf.speed < 8 && yv12_g && yv12_g != yv12 &&
        (cpi->ref_frame_flags & VP9_GOLD_FLAG)) {
      vp9_setup_pre_planes(xd, 0, yv12_g, mi_row, mi_col,
                           &cm->frame_refs[GOLDEN_FRAME - 1].sf);
      y_sad_g = cpi->fn_ptr[bsize].sdf(
          x->plane[0].src.buf, x->plane[0].src.stride, xd->plane[0].pre[0].buf,
          xd->plane[0].pre[0].stride);
    } else {
      y_sad_g = UINT_MAX;
    }

    if (cpi->oxcf.lag_in_frames > 0 && cpi->oxcf.rc_mode == VPX_VBR &&
        cpi->rc.is_src_frame_alt_ref) {
      yv12 = get_ref_frame_buffer(cpi, ALTREF_FRAME);
      vp9_setup_pre_planes(xd, 0, yv12, mi_row, mi_col,
                           &cm->frame_refs[ALTREF_FRAME - 1].sf);
      mi->ref_frame[0] = ALTREF_FRAME;
      y_sad_g = UINT_MAX;
    } else {
      vp9_setup_pre_planes(xd, 0, yv12, mi_row, mi_col,
                           &cm->frame_refs[LAST_FRAME - 1].sf);
      mi->ref_frame[0] = LAST_FRAME;
    }
    mi->ref_frame[1] = NONE;
    mi->sb_type = BLOCK_64X64;
    mi->mv[0].as_int = 0;
    mi->interp_filter = BILINEAR;

    if (cpi->oxcf.speed >= 8 && !low_res &&
        x->content_state_sb != kVeryHighSad) {
      y_sad = cpi->fn_ptr[bsize].sdf(
          x->plane[0].src.buf, x->plane[0].src.stride, xd->plane[0].pre[0].buf,
          xd->plane[0].pre[0].stride);
    } else {
      const MV dummy_mv = { 0, 0 };
      y_sad = vp9_int_pro_motion_estimation(cpi, x, bsize, mi_row, mi_col,
                                            &dummy_mv);
      x->sb_use_mv_part = 1;
      x->sb_mvcol_part = mi->mv[0].as_mv.col;
      x->sb_mvrow_part = mi->mv[0].as_mv.row;
    }

    y_sad_last = y_sad;
    // Pick ref frame for partitioning, bias last frame when y_sad_g and y_sad
    // are close if short_circuit_low_temp_var is on.
    y_sad_thr = cpi->sf.short_circuit_low_temp_var ? (y_sad * 7) >> 3 : y_sad;
    if (y_sad_g < y_sad_thr) {
      vp9_setup_pre_planes(xd, 0, yv12_g, mi_row, mi_col,
                           &cm->frame_refs[GOLDEN_FRAME - 1].sf);
      mi->ref_frame[0] = GOLDEN_FRAME;
      mi->mv[0].as_int = 0;
      y_sad = y_sad_g;
      ref_frame_partition = GOLDEN_FRAME;
    } else {
      x->pred_mv[LAST_FRAME] = mi->mv[0].as_mv;
      ref_frame_partition = LAST_FRAME;
    }

    set_ref_ptrs(cm, xd, mi->ref_frame[0], mi->ref_frame[1]);
    vp9_build_inter_predictors_sb(xd, mi_row, mi_col, BLOCK_64X64);

    if (cpi->use_skin_detection)
      x->sb_is_skin =
          skin_sb_split(cpi, x, low_res, mi_row, mi_col, force_split);

    d = xd->plane[0].dst.buf;
    dp = xd->plane[0].dst.stride;

    // If the y_sad is very small, take 64x64 as partition and exit.
    // Don't check on boosted segment for now, as 64x64 is suppressed there.
    if (segment_id == CR_SEGMENT_ID_BASE && y_sad < cpi->vbp_threshold_sad) {
      const int block_width = num_8x8_blocks_wide_lookup[BLOCK_64X64];
      const int block_height = num_8x8_blocks_high_lookup[BLOCK_64X64];
      if (mi_col + block_width / 2 < cm->mi_cols &&
          mi_row + block_height / 2 < cm->mi_rows) {
        set_block_size(cpi, x, xd, mi_row, mi_col, BLOCK_64X64);
        x->variance_low[0] = 1;
        chroma_check(cpi, x, bsize, y_sad, is_key_frame);
        if (cpi->sf.svc_use_lowres_part &&
            cpi->svc.spatial_layer_id == cpi->svc.number_spatial_layers - 2)
          update_partition_svc(cpi, BLOCK_64X64, mi_row, mi_col);
        if (cpi->sf.copy_partition_flag) {
          update_prev_partition(cpi, x, segment_id, mi_row, mi_col, sb_offset);
        }
        return 0;
      }
    }

    // If the y_sad is small enough, copy the partition of the superblock in the
    // last frame to current frame only if the last frame is not a keyframe.
    // Stop the copy every cpi->max_copied_frame to refresh the partition.
    // TODO(jianj) : tune the threshold.
    if (cpi->sf.copy_partition_flag && y_sad_last < cpi->vbp_threshold_copy &&
        copy_partitioning(cpi, x, xd, mi_row, mi_col, segment_id, sb_offset)) {
      chroma_check(cpi, x, bsize, y_sad, is_key_frame);
      if (cpi->sf.svc_use_lowres_part &&
          cpi->svc.spatial_layer_id == cpi->svc.number_spatial_layers - 2)
        update_partition_svc(cpi, BLOCK_64X64, mi_row, mi_col);
      return 0;
    }
  } else {
    d = VP9_VAR_OFFS;
    dp = 0;
#if CONFIG_VP9_HIGHBITDEPTH
    if (xd->cur_buf->flags & YV12_FLAG_HIGHBITDEPTH) {
      switch (xd->bd) {
        case 10: d = CONVERT_TO_BYTEPTR(VP9_HIGH_VAR_OFFS_10); break;
        case 12: d = CONVERT_TO_BYTEPTR(VP9_HIGH_VAR_OFFS_12); break;
        case 8:
        default: d = CONVERT_TO_BYTEPTR(VP9_HIGH_VAR_OFFS_8); break;
      }
    }
#endif  // CONFIG_VP9_HIGHBITDEPTH
  }

  if (low_res && threshold_4x4avg < INT64_MAX)
    CHECK_MEM_ERROR(cm, vt2, vpx_calloc(16, sizeof(*vt2)));
  // Fill in the entire tree of 8x8 (or 4x4 under some conditions) variances
  // for splits.
  for (i = 0; i < 4; i++) {
    const int x32_idx = ((i & 1) << 5);
    const int y32_idx = ((i >> 1) << 5);
    const int i2 = i << 2;
    force_split[i + 1] = 0;
    avg_16x16[i] = 0;
    maxvar_16x16[i] = 0;
    minvar_16x16[i] = INT_MAX;
    for (j = 0; j < 4; j++) {
      const int x16_idx = x32_idx + ((j & 1) << 4);
      const int y16_idx = y32_idx + ((j >> 1) << 4);
      const int split_index = 5 + i2 + j;
      v16x16 *vst = &vt.split[i].split[j];
      force_split[split_index] = 0;
      variance4x4downsample[i2 + j] = 0;
      if (!is_key_frame) {
        fill_variance_8x8avg(s, sp, d, dp, x16_idx, y16_idx, vst,
#if CONFIG_VP9_HIGHBITDEPTH
                             xd->cur_buf->flags,
#endif
                             pixels_wide, pixels_high, is_key_frame);
        fill_variance_tree(&vt.split[i].split[j], BLOCK_16X16);
        get_variance(&vt.split[i].split[j].part_variances.none);
        avg_16x16[i] += vt.split[i].split[j].part_variances.none.variance;
        if (vt.split[i].split[j].part_variances.none.variance < minvar_16x16[i])
          minvar_16x16[i] = vt.split[i].split[j].part_variances.none.variance;
        if (vt.split[i].split[j].part_variances.none.variance > maxvar_16x16[i])
          maxvar_16x16[i] = vt.split[i].split[j].part_variances.none.variance;
        if (vt.split[i].split[j].part_variances.none.variance > thresholds[2]) {
          // 16X16 variance is above threshold for split, so force split to 8x8
          // for this 16x16 block (this also forces splits for upper levels).
          force_split[split_index] = 1;
          force_split[i + 1] = 1;
          force_split[0] = 1;
        } else if (compute_minmax_variance &&
                   vt.split[i].split[j].part_variances.none.variance >
                       thresholds[1] &&
                   !cyclic_refresh_segment_id_boosted(segment_id)) {
          // We have some nominal amount of 16x16 variance (based on average),
          // compute the minmax over the 8x8 sub-blocks, and if above threshold,
          // force split to 8x8 block for this 16x16 block.
          int minmax = compute_minmax_8x8(s, sp, d, dp, x16_idx, y16_idx,
#if CONFIG_VP9_HIGHBITDEPTH
                                          xd->cur_buf->flags,
#endif
                                          pixels_wide, pixels_high);
          int thresh_minmax = (int)cpi->vbp_threshold_minmax;
          if (x->content_state_sb == kVeryHighSad)
            thresh_minmax = thresh_minmax << 1;
          if (minmax > thresh_minmax) {
            force_split[split_index] = 1;
            force_split[i + 1] = 1;
            force_split[0] = 1;
          }
        }
      }
      if (is_key_frame ||
          (low_res && vt.split[i].split[j].part_variances.none.variance >
                          threshold_4x4avg)) {
        force_split[split_index] = 0;
        // Go down to 4x4 down-sampling for variance.
        variance4x4downsample[i2 + j] = 1;
        for (k = 0; k < 4; k++) {
          int x8_idx = x16_idx + ((k & 1) << 3);
          int y8_idx = y16_idx + ((k >> 1) << 3);
          v8x8 *vst2 = is_key_frame ? &vst->split[k] : &vt2[i2 + j].split[k];
          fill_variance_4x4avg(s, sp, d, dp, x8_idx, y8_idx, vst2,
#if CONFIG_VP9_HIGHBITDEPTH
                               xd->cur_buf->flags,
#endif
                               pixels_wide, pixels_high, is_key_frame);
        }
      }
    }
  }
  if (cpi->noise_estimate.enabled)
    noise_level = vp9_noise_estimate_extract_level(&cpi->noise_estimate);
  // Fill the rest of the variance tree by summing split partition values.
  avg_32x32 = 0;
  for (i = 0; i < 4; i++) {
    const int i2 = i << 2;
    for (j = 0; j < 4; j++) {
      if (variance4x4downsample[i2 + j] == 1) {
        v16x16 *vtemp = (!is_key_frame) ? &vt2[i2 + j] : &vt.split[i].split[j];
        for (m = 0; m < 4; m++) fill_variance_tree(&vtemp->split[m], BLOCK_8X8);
        fill_variance_tree(vtemp, BLOCK_16X16);
        // If variance of this 16x16 block is above the threshold, force block
        // to split. This also forces a split on the upper levels.
        get_variance(&vtemp->part_variances.none);
        if (vtemp->part_variances.none.variance > thresholds[2]) {
          force_split[5 + i2 + j] = 1;
          force_split[i + 1] = 1;
          force_split[0] = 1;
        }
      }
    }
    fill_variance_tree(&vt.split[i], BLOCK_32X32);
    // If variance of this 32x32 block is above the threshold, or if its above
    // (some threshold of) the average variance over the sub-16x16 blocks, then
    // force this block to split. This also forces a split on the upper
    // (64x64) level.
    if (!force_split[i + 1]) {
      get_variance(&vt.split[i].part_variances.none);
      var_32x32 = vt.split[i].part_variances.none.variance;
      max_var_32x32 = VPXMAX(var_32x32, max_var_32x32);
      min_var_32x32 = VPXMIN(var_32x32, min_var_32x32);
      if (vt.split[i].part_variances.none.variance > thresholds[1] ||
          (!is_key_frame &&
           vt.split[i].part_variances.none.variance > (thresholds[1] >> 1) &&
           vt.split[i].part_variances.none.variance > (avg_16x16[i] >> 1))) {
        force_split[i + 1] = 1;
        force_split[0] = 1;
      } else if (!is_key_frame && noise_level < kLow && cm->height <= 360 &&
                 (maxvar_16x16[i] - minvar_16x16[i]) > (thresholds[1] >> 1) &&
                 maxvar_16x16[i] > thresholds[1]) {
        force_split[i + 1] = 1;
        force_split[0] = 1;
      }
      avg_32x32 += var_32x32;
    }
  }
  if (!force_split[0]) {
    fill_variance_tree(&vt, BLOCK_64X64);
    get_variance(&vt.part_variances.none);
    // If variance of this 64x64 block is above (some threshold of) the average
    // variance over the sub-32x32 blocks, then force this block to split.
    // Only checking this for noise level >= medium for now.
    if (!is_key_frame && noise_level >= kMedium &&
        vt.part_variances.none.variance > (9 * avg_32x32) >> 5)
      force_split[0] = 1;
    // Else if the maximum 32x32 variance minus the miniumum 32x32 variance in
    // a 64x64 block is greater than threshold and the maximum 32x32 variance is
    // above a miniumum threshold, then force the split of a 64x64 block
    // Only check this for low noise.
    else if (!is_key_frame && noise_level < kMedium &&
             (max_var_32x32 - min_var_32x32) > 3 * (thresholds[0] >> 3) &&
             max_var_32x32 > thresholds[0] >> 1)
      force_split[0] = 1;
  }

  // Now go through the entire structure, splitting every block size until
  // we get to one that's got a variance lower than our threshold.
  if (mi_col + 8 > cm->mi_cols || mi_row + 8 > cm->mi_rows ||
      !set_vt_partitioning(cpi, x, xd, &vt, BLOCK_64X64, mi_row, mi_col,
                           thresholds[0], BLOCK_16X16, force_split[0])) {
    for (i = 0; i < 4; ++i) {
      const int x32_idx = ((i & 1) << 2);
      const int y32_idx = ((i >> 1) << 2);
      const int i2 = i << 2;
      if (!set_vt_partitioning(cpi, x, xd, &vt.split[i], BLOCK_32X32,
                               (mi_row + y32_idx), (mi_col + x32_idx),
                               thresholds[1], BLOCK_16X16,
                               force_split[i + 1])) {
        for (j = 0; j < 4; ++j) {
          const int x16_idx = ((j & 1) << 1);
          const int y16_idx = ((j >> 1) << 1);
          // For inter frames: if variance4x4downsample[] == 1 for this 16x16
          // block, then the variance is based on 4x4 down-sampling, so use vt2
          // in set_vt_partioning(), otherwise use vt.
          v16x16 *vtemp = (!is_key_frame && variance4x4downsample[i2 + j] == 1)
                              ? &vt2[i2 + j]
                              : &vt.split[i].split[j];
          if (!set_vt_partitioning(
                  cpi, x, xd, vtemp, BLOCK_16X16, mi_row + y32_idx + y16_idx,
                  mi_col + x32_idx + x16_idx, thresholds[2], cpi->vbp_bsize_min,
                  force_split[5 + i2 + j])) {
            for (k = 0; k < 4; ++k) {
              const int x8_idx = (k & 1);
              const int y8_idx = (k >> 1);
              if (use_4x4_partition) {
                if (!set_vt_partitioning(cpi, x, xd, &vtemp->split[k],
                                         BLOCK_8X8,
                                         mi_row + y32_idx + y16_idx + y8_idx,
                                         mi_col + x32_idx + x16_idx + x8_idx,
                                         thresholds[3], BLOCK_8X8, 0)) {
                  set_block_size(
                      cpi, x, xd, (mi_row + y32_idx + y16_idx + y8_idx),
                      (mi_col + x32_idx + x16_idx + x8_idx), BLOCK_4X4);
                }
              } else {
                set_block_size(
                    cpi, x, xd, (mi_row + y32_idx + y16_idx + y8_idx),
                    (mi_col + x32_idx + x16_idx + x8_idx), BLOCK_8X8);
              }
            }
          }
        }
      }
    }
  }

  if (!frame_is_intra_only(cm) && cpi->sf.copy_partition_flag) {
    update_prev_partition(cpi, x, segment_id, mi_row, mi_col, sb_offset);
  }

  if (!frame_is_intra_only(cm) && cpi->sf.svc_use_lowres_part &&
      cpi->svc.spatial_layer_id == cpi->svc.number_spatial_layers - 2)
    update_partition_svc(cpi, BLOCK_64X64, mi_row, mi_col);

  if (cpi->sf.short_circuit_low_temp_var) {
    set_low_temp_var_flag(cpi, x, xd, &vt, thresholds, ref_frame_partition,
                          mi_col, mi_row);
  }

  chroma_check(cpi, x, bsize, y_sad, is_key_frame);
  if (vt2) vpx_free(vt2);
  return 0;
}

static void update_state(VP9_COMP *cpi, ThreadData *td, PICK_MODE_CONTEXT *ctx,
                         int mi_row, int mi_col, BLOCK_SIZE bsize,
                         int output_enabled) {
  int i, x_idx, y;
  VP9_COMMON *const cm = &cpi->common;
  RD_COUNTS *const rdc = &td->rd_counts;
  MACROBLOCK *const x = &td->mb;
  MACROBLOCKD *const xd = &x->e_mbd;
  struct macroblock_plane *const p = x->plane;
  struct macroblockd_plane *const pd = xd->plane;
  MODE_INFO *mi = &ctx->mic;
  MODE_INFO *const xdmi = xd->mi[0];
  MODE_INFO *mi_addr = xd->mi[0];
  const struct segmentation *const seg = &cm->seg;
  const int bw = num_8x8_blocks_wide_lookup[mi->sb_type];
  const int bh = num_8x8_blocks_high_lookup[mi->sb_type];
  const int x_mis = VPXMIN(bw, cm->mi_cols - mi_col);
  const int y_mis = VPXMIN(bh, cm->mi_rows - mi_row);
  MV_REF *const frame_mvs = cm->cur_frame->mvs + mi_row * cm->mi_cols + mi_col;
  int w, h;

  const int mis = cm->mi_stride;
  const int mi_width = num_8x8_blocks_wide_lookup[bsize];
  const int mi_height = num_8x8_blocks_high_lookup[bsize];
  int max_plane;

  assert(mi->sb_type == bsize);

  *mi_addr = *mi;
  *x->mbmi_ext = ctx->mbmi_ext;

  // If segmentation in use
  if (seg->enabled) {
    // For in frame complexity AQ copy the segment id from the segment map.
    if (cpi->oxcf.aq_mode == COMPLEXITY_AQ) {
      const uint8_t *const map =
          seg->update_map ? cpi->segmentation_map : cm->last_frame_seg_map;
      mi_addr->segment_id = get_segment_id(cm, map, bsize, mi_row, mi_col);
    }
    // Else for cyclic refresh mode update the segment map, set the segment id
    // and then update the quantizer.
    if (cpi->oxcf.aq_mode == CYCLIC_REFRESH_AQ) {
      vp9_cyclic_refresh_update_segment(cpi, xd->mi[0], mi_row, mi_col, bsize,
                                        ctx->rate, ctx->dist, x->skip, p);
    }
  }

  max_plane = is_inter_block(xdmi) ? MAX_MB_PLANE : 1;
  for (i = 0; i < max_plane; ++i) {
    p[i].coeff = ctx->coeff_pbuf[i][1];
    p[i].qcoeff = ctx->qcoeff_pbuf[i][1];
    pd[i].dqcoeff = ctx->dqcoeff_pbuf[i][1];
    p[i].eobs = ctx->eobs_pbuf[i][1];
  }

  for (i = max_plane; i < MAX_MB_PLANE; ++i) {
    p[i].coeff = ctx->coeff_pbuf[i][2];
    p[i].qcoeff = ctx->qcoeff_pbuf[i][2];
    pd[i].dqcoeff = ctx->dqcoeff_pbuf[i][2];
    p[i].eobs = ctx->eobs_pbuf[i][2];
  }

  // Restore the coding context of the MB to that that was in place
  // when the mode was picked for it
  for (y = 0; y < mi_height; y++)
    for (x_idx = 0; x_idx < mi_width; x_idx++)
      if ((xd->mb_to_right_edge >> (3 + MI_SIZE_LOG2)) + mi_width > x_idx &&
          (xd->mb_to_bottom_edge >> (3 + MI_SIZE_LOG2)) + mi_height > y) {
        xd->mi[x_idx + y * mis] = mi_addr;
      }

  if (cpi->oxcf.aq_mode != NO_AQ) vp9_init_plane_quantizers(cpi, x);

  if (is_inter_block(xdmi) && xdmi->sb_type < BLOCK_8X8) {
    xdmi->mv[0].as_int = mi->bmi[3].as_mv[0].as_int;
    xdmi->mv[1].as_int = mi->bmi[3].as_mv[1].as_int;
  }

  x->skip = ctx->skip;
  memcpy(x->zcoeff_blk[xdmi->tx_size], ctx->zcoeff_blk,
         sizeof(ctx->zcoeff_blk[0]) * ctx->num_4x4_blk);

  if (!output_enabled) return;

#if CONFIG_INTERNAL_STATS
  if (frame_is_intra_only(cm)) {
    static const int kf_mode_index[] = {
      THR_DC /*DC_PRED*/,          THR_V_PRED /*V_PRED*/,
      THR_H_PRED /*H_PRED*/,       THR_D45_PRED /*D45_PRED*/,
      THR_D135_PRED /*D135_PRED*/, THR_D117_PRED /*D117_PRED*/,
      THR_D153_PRED /*D153_PRED*/, THR_D207_PRED /*D207_PRED*/,
      THR_D63_PRED /*D63_PRED*/,   THR_TM /*TM_PRED*/,
    };
    ++cpi->mode_chosen_counts[kf_mode_index[xdmi->mode]];
  } else {
    // Note how often each mode chosen as best
    ++cpi->mode_chosen_counts[ctx->best_mode_index];
  }
#endif
  if (!frame_is_intra_only(cm)) {
    if (is_inter_block(xdmi)) {
      vp9_update_mv_count(td);

      if (cm->interp_filter == SWITCHABLE) {
        const int ctx = get_pred_context_switchable_interp(xd);
        ++td->counts->switchable_interp[ctx][xdmi->interp_filter];
      }
    }

    rdc->comp_pred_diff[SINGLE_REFERENCE] += ctx->single_pred_diff;
    rdc->comp_pred_diff[COMPOUND_REFERENCE] += ctx->comp_pred_diff;
    rdc->comp_pred_diff[REFERENCE_MODE_SELECT] += ctx->hybrid_pred_diff;

    for (i = 0; i < SWITCHABLE_FILTER_CONTEXTS; ++i)
      rdc->filter_diff[i] += ctx->best_filter_diff[i];
  }

  for (h = 0; h < y_mis; ++h) {
    MV_REF *const frame_mv = frame_mvs + h * cm->mi_cols;
    for (w = 0; w < x_mis; ++w) {
      MV_REF *const mv = frame_mv + w;
      mv->ref_frame[0] = mi->ref_frame[0];
      mv->ref_frame[1] = mi->ref_frame[1];
      mv->mv[0].as_int = mi->mv[0].as_int;
      mv->mv[1].as_int = mi->mv[1].as_int;
    }
  }
}

void vp9_setup_src_planes(MACROBLOCK *x, const YV12_BUFFER_CONFIG *src,
                          int mi_row, int mi_col) {
  uint8_t *const buffers[3] = { src->y_buffer, src->u_buffer, src->v_buffer };
  const int strides[3] = { src->y_stride, src->uv_stride, src->uv_stride };
  int i;

  // Set current frame pointer.
  x->e_mbd.cur_buf = src;

  for (i = 0; i < MAX_MB_PLANE; i++)
    setup_pred_plane(&x->plane[i].src, buffers[i], strides[i], mi_row, mi_col,
                     NULL, x->e_mbd.plane[i].subsampling_x,
                     x->e_mbd.plane[i].subsampling_y);
}

static void set_mode_info_seg_skip(MACROBLOCK *x, TX_MODE tx_mode,
                                   RD_COST *rd_cost, BLOCK_SIZE bsize) {
  MACROBLOCKD *const xd = &x->e_mbd;
  MODE_INFO *const mi = xd->mi[0];
  INTERP_FILTER filter_ref;

  filter_ref = get_pred_context_switchable_interp(xd);
  if (filter_ref == SWITCHABLE_FILTERS) filter_ref = EIGHTTAP;

  mi->sb_type = bsize;
  mi->mode = ZEROMV;
  mi->tx_size =
      VPXMIN(max_txsize_lookup[bsize], tx_mode_to_biggest_tx_size[tx_mode]);
  mi->skip = 1;
  mi->uv_mode = DC_PRED;
  mi->ref_frame[0] = LAST_FRAME;
  mi->ref_frame[1] = NONE;
  mi->mv[0].as_int = 0;
  mi->interp_filter = filter_ref;

  xd->mi[0]->bmi[0].as_mv[0].as_int = 0;
  x->skip = 1;

  vp9_rd_cost_init(rd_cost);
}

static int set_segment_rdmult(VP9_COMP *const cpi, MACROBLOCK *const x,
                              int8_t segment_id) {
  int segment_qindex;
  VP9_COMMON *const cm = &cpi->common;
  vp9_init_plane_quantizers(cpi, x);
  vpx_clear_system_state();
  segment_qindex = vp9_get_qindex(&cm->seg, segment_id, cm->base_qindex);
  return vp9_compute_rd_mult(cpi, segment_qindex + cm->y_dc_delta_q);
}

static void rd_pick_sb_modes(VP9_COMP *cpi, TileDataEnc *tile_data,
                             MACROBLOCK *const x, int mi_row, int mi_col,
                             RD_COST *rd_cost, BLOCK_SIZE bsize,
                             PICK_MODE_CONTEXT *ctx, int64_t best_rd) {
  VP9_COMMON *const cm = &cpi->common;
  TileInfo *const tile_info = &tile_data->tile_info;
  MACROBLOCKD *const xd = &x->e_mbd;
  MODE_INFO *mi;
  struct macroblock_plane *const p = x->plane;
  struct macroblockd_plane *const pd = xd->plane;
  const AQ_MODE aq_mode = cpi->oxcf.aq_mode;
  int i, orig_rdmult;

  vpx_clear_system_state();

  // Use the lower precision, but faster, 32x32 fdct for mode selection.
  x->use_lp32x32fdct = 1;

  set_offsets(cpi, tile_info, x, mi_row, mi_col, bsize);
  mi = xd->mi[0];
  mi->sb_type = bsize;

  for (i = 0; i < MAX_MB_PLANE; ++i) {
    p[i].coeff = ctx->coeff_pbuf[i][0];
    p[i].qcoeff = ctx->qcoeff_pbuf[i][0];
    pd[i].dqcoeff = ctx->dqcoeff_pbuf[i][0];
    p[i].eobs = ctx->eobs_pbuf[i][0];
  }
  ctx->is_coded = 0;
  ctx->skippable = 0;
  ctx->pred_pixel_ready = 0;
  x->skip_recode = 0;

  // Set to zero to make sure we do not use the previous encoded frame stats
  mi->skip = 0;

#if CONFIG_VP9_HIGHBITDEPTH
  if (xd->cur_buf->flags & YV12_FLAG_HIGHBITDEPTH) {
    x->source_variance = vp9_high_get_sby_perpixel_variance(
        cpi, &x->plane[0].src, bsize, xd->bd);
  } else {
    x->source_variance =
        vp9_get_sby_perpixel_variance(cpi, &x->plane[0].src, bsize);
  }
#else
  x->source_variance =
      vp9_get_sby_perpixel_variance(cpi, &x->plane[0].src, bsize);
#endif  // CONFIG_VP9_HIGHBITDEPTH

  // Save rdmult before it might be changed, so it can be restored later.
  orig_rdmult = x->rdmult;

  if ((cpi->sf.tx_domain_thresh > 0.0) || (cpi->sf.quant_opt_thresh > 0.0)) {
    double logvar = vp9_log_block_var(cpi, x, bsize);
    // Check block complexity as part of descision on using pixel or transform
    // domain distortion in rd tests.
    x->block_tx_domain = cpi->sf.allow_txfm_domain_distortion &&
                         (logvar >= cpi->sf.tx_domain_thresh);

    // Check block complexity as part of descision on using quantized
    // coefficient optimisation inside the rd loop.
    x->block_qcoeff_opt =
        cpi->sf.allow_quant_coeff_opt && (logvar <= cpi->sf.quant_opt_thresh);
  } else {
    x->block_tx_domain = cpi->sf.allow_txfm_domain_distortion;
    x->block_qcoeff_opt = cpi->sf.allow_quant_coeff_opt;
  }

  if (aq_mode == VARIANCE_AQ) {
    if (cm->frame_type == KEY_FRAME || cpi->refresh_alt_ref_frame ||
        cpi->force_update_segmentation ||
        (cpi->refresh_golden_frame && !cpi->rc.is_src_frame_alt_ref)) {
      int min_energy;
      int max_energy;

      // Get sub block energy range
      if (bsize >= BLOCK_32X32) {
        vp9_get_sub_block_energy(cpi, x, mi_row, mi_col, bsize, &min_energy,
                                 &max_energy);
      } else {
        min_energy = bsize <= BLOCK_16X16 ? x->mb_energy
                                          : vp9_block_energy(cpi, x, bsize);
      }

      mi->segment_id = vp9_vaq_segment_id(min_energy);
    } else {
      const uint8_t *const map =
          cm->seg.update_map ? cpi->segmentation_map : cm->last_frame_seg_map;
      mi->segment_id = get_segment_id(cm, map, bsize, mi_row, mi_col);
    }
    x->rdmult = set_segment_rdmult(cpi, x, mi->segment_id);
  } else if (aq_mode == LOOKAHEAD_AQ) {
    const uint8_t *const map = cpi->segmentation_map;

    // I do not change rdmult here consciously.
    mi->segment_id = get_segment_id(cm, map, bsize, mi_row, mi_col);
  } else if (aq_mode == EQUATOR360_AQ) {
    if (cm->frame_type == KEY_FRAME || cpi->force_update_segmentation) {
      mi->segment_id = vp9_360aq_segment_id(mi_row, cm->mi_rows);
    } else {
      const uint8_t *const map =
          cm->seg.update_map ? cpi->segmentation_map : cm->last_frame_seg_map;
      mi->segment_id = get_segment_id(cm, map, bsize, mi_row, mi_col);
    }
    x->rdmult = set_segment_rdmult(cpi, x, mi->segment_id);
  } else if (aq_mode == COMPLEXITY_AQ) {
    x->rdmult = set_segment_rdmult(cpi, x, mi->segment_id);
  } else if (aq_mode == CYCLIC_REFRESH_AQ) {
    const uint8_t *const map =
        cm->seg.update_map ? cpi->segmentation_map : cm->last_frame_seg_map;
    // If segment is boosted, use rdmult for that segment.
    if (cyclic_refresh_segment_id_boosted(
            get_segment_id(cm, map, bsize, mi_row, mi_col)))
      x->rdmult = vp9_cyclic_refresh_get_rdmult(cpi->cyclic_refresh);
  } else {
    if (cpi->sf.enable_tpl_model) x->rdmult = x->cb_rdmult;
  }

  // Find best coding mode & reconstruct the MB so it is available
  // as a predictor for MBs that follow in the SB
  if (frame_is_intra_only(cm)) {
    vp9_rd_pick_intra_mode_sb(cpi, x, rd_cost, bsize, ctx, best_rd);
  } else {
    if (bsize >= BLOCK_8X8) {
      if (segfeature_active(&cm->seg, mi->segment_id, SEG_LVL_SKIP))
        vp9_rd_pick_inter_mode_sb_seg_skip(cpi, tile_data, x, rd_cost, bsize,
                                           ctx, best_rd);
      else
        vp9_rd_pick_inter_mode_sb(cpi, tile_data, x, mi_row, mi_col, rd_cost,
                                  bsize, ctx, best_rd);
    } else {
      vp9_rd_pick_inter_mode_sub8x8(cpi, tile_data, x, mi_row, mi_col, rd_cost,
                                    bsize, ctx, best_rd);
    }
  }

  // Examine the resulting rate and for AQ mode 2 make a segment choice.
  if ((rd_cost->rate != INT_MAX) && (aq_mode == COMPLEXITY_AQ) &&
      (bsize >= BLOCK_16X16) &&
      (cm->frame_type == KEY_FRAME || cpi->refresh_alt_ref_frame ||
       (cpi->refresh_golden_frame && !cpi->rc.is_src_frame_alt_ref))) {
    vp9_caq_select_segment(cpi, x, bsize, mi_row, mi_col, rd_cost->rate);
  }

  // TODO(jingning) The rate-distortion optimization flow needs to be
  // refactored to provide proper exit/return handle.
  if (rd_cost->rate == INT_MAX)
    rd_cost->rdcost = INT64_MAX;
  else
    rd_cost->rdcost = RDCOST(x->rdmult, x->rddiv, rd_cost->rate, rd_cost->dist);

  x->rdmult = orig_rdmult;

  ctx->rate = rd_cost->rate;
  ctx->dist = rd_cost->dist;
}

static void update_stats(VP9_COMMON *cm, ThreadData *td) {
  const MACROBLOCK *x = &td->mb;
  const MACROBLOCKD *const xd = &x->e_mbd;
  const MODE_INFO *const mi = xd->mi[0];
  const MB_MODE_INFO_EXT *const mbmi_ext = x->mbmi_ext;
  const BLOCK_SIZE bsize = mi->sb_type;

  if (!frame_is_intra_only(cm)) {
    FRAME_COUNTS *const counts = td->counts;
    const int inter_block = is_inter_block(mi);
    const int seg_ref_active =
        segfeature_active(&cm->seg, mi->segment_id, SEG_LVL_REF_FRAME);
    if (!seg_ref_active) {
      counts->intra_inter[get_intra_inter_context(xd)][inter_block]++;
      // If the segment reference feature is enabled we have only a single
      // reference frame allowed for the segment so exclude it from
      // the reference frame counts used to work out probabilities.
      if (inter_block) {
        const MV_REFERENCE_FRAME ref0 = mi->ref_frame[0];
        if (cm->reference_mode == REFERENCE_MODE_SELECT)
          counts->comp_inter[vp9_get_reference_mode_context(cm, xd)]
                            [has_second_ref(mi)]++;

        if (has_second_ref(mi)) {
          counts->comp_ref[vp9_get_pred_context_comp_ref_p(cm, xd)]
                          [ref0 == GOLDEN_FRAME]++;
        } else {
          counts->single_ref[vp9_get_pred_context_single_ref_p1(xd)][0]
                            [ref0 != LAST_FRAME]++;
          if (ref0 != LAST_FRAME)
            counts->single_ref[vp9_get_pred_context_single_ref_p2(xd)][1]
                              [ref0 != GOLDEN_FRAME]++;
        }
      }
    }
    if (inter_block &&
        !segfeature_active(&cm->seg, mi->segment_id, SEG_LVL_SKIP)) {
      const int mode_ctx = mbmi_ext->mode_context[mi->ref_frame[0]];
      if (bsize >= BLOCK_8X8) {
        const PREDICTION_MODE mode = mi->mode;
        ++counts->inter_mode[mode_ctx][INTER_OFFSET(mode)];
      } else {
        const int num_4x4_w = num_4x4_blocks_wide_lookup[bsize];
        const int num_4x4_h = num_4x4_blocks_high_lookup[bsize];
        int idx, idy;
        for (idy = 0; idy < 2; idy += num_4x4_h) {
          for (idx = 0; idx < 2; idx += num_4x4_w) {
            const int j = idy * 2 + idx;
            const PREDICTION_MODE b_mode = mi->bmi[j].as_mode;
            ++counts->inter_mode[mode_ctx][INTER_OFFSET(b_mode)];
          }
        }
      }
    }
  }
}

static void restore_context(MACROBLOCK *const x, int mi_row, int mi_col,
                            ENTROPY_CONTEXT a[16 * MAX_MB_PLANE],
                            ENTROPY_CONTEXT l[16 * MAX_MB_PLANE],
                            PARTITION_CONTEXT sa[8], PARTITION_CONTEXT sl[8],
                            BLOCK_SIZE bsize) {
  MACROBLOCKD *const xd = &x->e_mbd;
  int p;
  const int num_4x4_blocks_wide = num_4x4_blocks_wide_lookup[bsize];
  const int num_4x4_blocks_high = num_4x4_blocks_high_lookup[bsize];
  int mi_width = num_8x8_blocks_wide_lookup[bsize];
  int mi_height = num_8x8_blocks_high_lookup[bsize];
  for (p = 0; p < MAX_MB_PLANE; p++) {
    memcpy(xd->above_context[p] + ((mi_col * 2) >> xd->plane[p].subsampling_x),
           a + num_4x4_blocks_wide * p,
           (sizeof(ENTROPY_CONTEXT) * num_4x4_blocks_wide) >>
               xd->plane[p].subsampling_x);
    memcpy(xd->left_context[p] +
               ((mi_row & MI_MASK) * 2 >> xd->plane[p].subsampling_y),
           l + num_4x4_blocks_high * p,
           (sizeof(ENTROPY_CONTEXT) * num_4x4_blocks_high) >>
               xd->plane[p].subsampling_y);
  }
  memcpy(xd->above_seg_context + mi_col, sa,
         sizeof(*xd->above_seg_context) * mi_width);
  memcpy(xd->left_seg_context + (mi_row & MI_MASK), sl,
         sizeof(xd->left_seg_context[0]) * mi_height);
}

static void save_context(MACROBLOCK *const x, int mi_row, int mi_col,
                         ENTROPY_CONTEXT a[16 * MAX_MB_PLANE],
                         ENTROPY_CONTEXT l[16 * MAX_MB_PLANE],
                         PARTITION_CONTEXT sa[8], PARTITION_CONTEXT sl[8],
                         BLOCK_SIZE bsize) {
  const MACROBLOCKD *const xd = &x->e_mbd;
  int p;
  const int num_4x4_blocks_wide = num_4x4_blocks_wide_lookup[bsize];
  const int num_4x4_blocks_high = num_4x4_blocks_high_lookup[bsize];
  int mi_width = num_8x8_blocks_wide_lookup[bsize];
  int mi_height = num_8x8_blocks_high_lookup[bsize];

  // buffer the above/left context information of the block in search.
  for (p = 0; p < MAX_MB_PLANE; ++p) {
    memcpy(a + num_4x4_blocks_wide * p,
           xd->above_context[p] + (mi_col * 2 >> xd->plane[p].subsampling_x),
           (sizeof(ENTROPY_CONTEXT) * num_4x4_blocks_wide) >>
               xd->plane[p].subsampling_x);
    memcpy(l + num_4x4_blocks_high * p,
           xd->left_context[p] +
               ((mi_row & MI_MASK) * 2 >> xd->plane[p].subsampling_y),
           (sizeof(ENTROPY_CONTEXT) * num_4x4_blocks_high) >>
               xd->plane[p].subsampling_y);
  }
  memcpy(sa, xd->above_seg_context + mi_col,
         sizeof(*xd->above_seg_context) * mi_width);
  memcpy(sl, xd->left_seg_context + (mi_row & MI_MASK),
         sizeof(xd->left_seg_context[0]) * mi_height);
}

static void encode_b(VP9_COMP *cpi, const TileInfo *const tile, ThreadData *td,
                     TOKENEXTRA **tp, int mi_row, int mi_col,
                     int output_enabled, BLOCK_SIZE bsize,
                     PICK_MODE_CONTEXT *ctx) {
  MACROBLOCK *const x = &td->mb;
  set_offsets(cpi, tile, x, mi_row, mi_col, bsize);

  if (cpi->sf.enable_tpl_model && cpi->oxcf.aq_mode == NO_AQ)
    x->rdmult = x->cb_rdmult;

  update_state(cpi, td, ctx, mi_row, mi_col, bsize, output_enabled);
  encode_superblock(cpi, td, tp, output_enabled, mi_row, mi_col, bsize, ctx);

  if (output_enabled) {
    update_stats(&cpi->common, td);

    (*tp)->token = EOSB_TOKEN;
    (*tp)++;
  }
}

static void encode_sb(VP9_COMP *cpi, ThreadData *td, const TileInfo *const tile,
                      TOKENEXTRA **tp, int mi_row, int mi_col,
                      int output_enabled, BLOCK_SIZE bsize, PC_TREE *pc_tree) {
  VP9_COMMON *const cm = &cpi->common;
  MACROBLOCK *const x = &td->mb;
  MACROBLOCKD *const xd = &x->e_mbd;

  const int bsl = b_width_log2_lookup[bsize], hbs = (1 << bsl) / 4;
  int ctx;
  PARTITION_TYPE partition;
  BLOCK_SIZE subsize = bsize;

  if (mi_row >= cm->mi_rows || mi_col >= cm->mi_cols) return;

  if (bsize >= BLOCK_8X8) {
    ctx = partition_plane_context(xd, mi_row, mi_col, bsize);
    subsize = get_subsize(bsize, pc_tree->partitioning);
  } else {
    ctx = 0;
    subsize = BLOCK_4X4;
  }

  partition = partition_lookup[bsl][subsize];
  if (output_enabled && bsize != BLOCK_4X4)
    td->counts->partition[ctx][partition]++;

  switch (partition) {
    case PARTITION_NONE:
      encode_b(cpi, tile, td, tp, mi_row, mi_col, output_enabled, subsize,
               &pc_tree->none);
      break;
    case PARTITION_VERT:
      encode_b(cpi, tile, td, tp, mi_row, mi_col, output_enabled, subsize,
               &pc_tree->vertical[0]);
      if (mi_col + hbs < cm->mi_cols && bsize > BLOCK_8X8) {
        encode_b(cpi, tile, td, tp, mi_row, mi_col + hbs, output_enabled,
                 subsize, &pc_tree->vertical[1]);
      }
      break;
    case PARTITION_HORZ:
      encode_b(cpi, tile, td, tp, mi_row, mi_col, output_enabled, subsize,
               &pc_tree->horizontal[0]);
      if (mi_row + hbs < cm->mi_rows && bsize > BLOCK_8X8) {
        encode_b(cpi, tile, td, tp, mi_row + hbs, mi_col, output_enabled,
                 subsize, &pc_tree->horizontal[1]);
      }
      break;
    default:
      assert(partition == PARTITION_SPLIT);
      if (bsize == BLOCK_8X8) {
        encode_b(cpi, tile, td, tp, mi_row, mi_col, output_enabled, subsize,
                 pc_tree->leaf_split[0]);
      } else {
        encode_sb(cpi, td, tile, tp, mi_row, mi_col, output_enabled, subsize,
                  pc_tree->split[0]);
        encode_sb(cpi, td, tile, tp, mi_row, mi_col + hbs, output_enabled,
                  subsize, pc_tree->split[1]);
        encode_sb(cpi, td, tile, tp, mi_row + hbs, mi_col, output_enabled,
                  subsize, pc_tree->split[2]);
        encode_sb(cpi, td, tile, tp, mi_row + hbs, mi_col + hbs, output_enabled,
                  subsize, pc_tree->split[3]);
      }
      break;
  }

  if (partition != PARTITION_SPLIT || bsize == BLOCK_8X8)
    update_partition_context(xd, mi_row, mi_col, subsize, bsize);
}

// Check to see if the given partition size is allowed for a specified number
// of 8x8 block rows and columns remaining in the image.
// If not then return the largest allowed partition size
static BLOCK_SIZE find_partition_size(BLOCK_SIZE bsize, int rows_left,
                                      int cols_left, int *bh, int *bw) {
  if (rows_left <= 0 || cols_left <= 0) {
    return VPXMIN(bsize, BLOCK_8X8);
  } else {
    for (; bsize > 0; bsize -= 3) {
      *bh = num_8x8_blocks_high_lookup[bsize];
      *bw = num_8x8_blocks_wide_lookup[bsize];
      if ((*bh <= rows_left) && (*bw <= cols_left)) {
        break;
      }
    }
  }
  return bsize;
}

static void set_partial_b64x64_partition(MODE_INFO *mi, int mis, int bh_in,
                                         int bw_in, int row8x8_remaining,
                                         int col8x8_remaining, BLOCK_SIZE bsize,
                                         MODE_INFO **mi_8x8) {
  int bh = bh_in;
  int r, c;
  for (r = 0; r < MI_BLOCK_SIZE; r += bh) {
    int bw = bw_in;
    for (c = 0; c < MI_BLOCK_SIZE; c += bw) {
      const int index = r * mis + c;
      mi_8x8[index] = mi + index;
      mi_8x8[index]->sb_type = find_partition_size(
          bsize, row8x8_remaining - r, col8x8_remaining - c, &bh, &bw);
    }
  }
}

// This function attempts to set all mode info entries in a given SB64
// to the same block partition size.
// However, at the bottom and right borders of the image the requested size
// may not be allowed in which case this code attempts to choose the largest
// allowable partition.
static void set_fixed_partitioning(VP9_COMP *cpi, const TileInfo *const tile,
                                   MODE_INFO **mi_8x8, int mi_row, int mi_col,
                                   BLOCK_SIZE bsize) {
  VP9_COMMON *const cm = &cpi->common;
  const int mis = cm->mi_stride;
  const int row8x8_remaining = tile->mi_row_end - mi_row;
  const int col8x8_remaining = tile->mi_col_end - mi_col;
  int block_row, block_col;
  MODE_INFO *mi_upper_left = cm->mi + mi_row * mis + mi_col;
  int bh = num_8x8_blocks_high_lookup[bsize];
  int bw = num_8x8_blocks_wide_lookup[bsize];

  assert((row8x8_remaining > 0) && (col8x8_remaining > 0));

  // Apply the requested partition size to the SB64 if it is all "in image"
  if ((col8x8_remaining >= MI_BLOCK_SIZE) &&
      (row8x8_remaining >= MI_BLOCK_SIZE)) {
    for (block_row = 0; block_row < MI_BLOCK_SIZE; block_row += bh) {
      for (block_col = 0; block_col < MI_BLOCK_SIZE; block_col += bw) {
        int index = block_row * mis + block_col;
        mi_8x8[index] = mi_upper_left + index;
        mi_8x8[index]->sb_type = bsize;
      }
    }
  } else {
    // Else this is a partial SB64.
    set_partial_b64x64_partition(mi_upper_left, mis, bh, bw, row8x8_remaining,
                                 col8x8_remaining, bsize, mi_8x8);
  }
}

static const struct {
  int row;
  int col;
} coord_lookup[16] = {
  // 32x32 index = 0
  { 0, 0 },
  { 0, 2 },
  { 2, 0 },
  { 2, 2 },
  // 32x32 index = 1
  { 0, 4 },
  { 0, 6 },
  { 2, 4 },
  { 2, 6 },
  // 32x32 index = 2
  { 4, 0 },
  { 4, 2 },
  { 6, 0 },
  { 6, 2 },
  // 32x32 index = 3
  { 4, 4 },
  { 4, 6 },
  { 6, 4 },
  { 6, 6 },
};

static void set_source_var_based_partition(VP9_COMP *cpi,
                                           const TileInfo *const tile,
                                           MACROBLOCK *const x,
                                           MODE_INFO **mi_8x8, int mi_row,
                                           int mi_col) {
  VP9_COMMON *const cm = &cpi->common;
  const int mis = cm->mi_stride;
  const int row8x8_remaining = tile->mi_row_end - mi_row;
  const int col8x8_remaining = tile->mi_col_end - mi_col;
  MODE_INFO *mi_upper_left = cm->mi + mi_row * mis + mi_col;

  vp9_setup_src_planes(x, cpi->Source, mi_row, mi_col);

  assert((row8x8_remaining > 0) && (col8x8_remaining > 0));

  // In-image SB64
  if ((col8x8_remaining >= MI_BLOCK_SIZE) &&
      (row8x8_remaining >= MI_BLOCK_SIZE)) {
    int i, j;
    int index;
    diff d32[4];
    const int offset = (mi_row >> 1) * cm->mb_cols + (mi_col >> 1);
    int is_larger_better = 0;
    int use32x32 = 0;
    unsigned int thr = cpi->source_var_thresh;

    memset(d32, 0, 4 * sizeof(diff));

    for (i = 0; i < 4; i++) {
      diff *d16[4];

      for (j = 0; j < 4; j++) {
        int b_mi_row = coord_lookup[i * 4 + j].row;
        int b_mi_col = coord_lookup[i * 4 + j].col;
        int boffset = b_mi_row / 2 * cm->mb_cols + b_mi_col / 2;

        d16[j] = cpi->source_diff_var + offset + boffset;

        index = b_mi_row * mis + b_mi_col;
        mi_8x8[index] = mi_upper_left + index;
        mi_8x8[index]->sb_type = BLOCK_16X16;

        // TODO(yunqingwang): If d16[j].var is very large, use 8x8 partition
        // size to further improve quality.
      }

      is_larger_better = (d16[0]->var < thr) && (d16[1]->var < thr) &&
                         (d16[2]->var < thr) && (d16[3]->var < thr);

      // Use 32x32 partition
      if (is_larger_better) {
        use32x32 += 1;

        for (j = 0; j < 4; j++) {
          d32[i].sse += d16[j]->sse;
          d32[i].sum += d16[j]->sum;
        }

        d32[i].var =
            (unsigned int)(d32[i].sse -
                           (unsigned int)(((int64_t)d32[i].sum * d32[i].sum) >>
                                          10));

        index = coord_lookup[i * 4].row * mis + coord_lookup[i * 4].col;
        mi_8x8[index] = mi_upper_left + index;
        mi_8x8[index]->sb_type = BLOCK_32X32;
      }
    }

    if (use32x32 == 4) {
      thr <<= 1;
      is_larger_better = (d32[0].var < thr) && (d32[1].var < thr) &&
                         (d32[2].var < thr) && (d32[3].var < thr);

      // Use 64x64 partition
      if (is_larger_better) {
        mi_8x8[0] = mi_upper_left;
        mi_8x8[0]->sb_type = BLOCK_64X64;
      }
    }
  } else {  // partial in-image SB64
    int bh = num_8x8_blocks_high_lookup[BLOCK_16X16];
    int bw = num_8x8_blocks_wide_lookup[BLOCK_16X16];
    set_partial_b64x64_partition(mi_upper_left, mis, bh, bw, row8x8_remaining,
                                 col8x8_remaining, BLOCK_16X16, mi_8x8);
  }
}

static void update_state_rt(VP9_COMP *cpi, ThreadData *td,
                            PICK_MODE_CONTEXT *ctx, int mi_row, int mi_col,
                            int bsize) {
  VP9_COMMON *const cm = &cpi->common;
  MACROBLOCK *const x = &td->mb;
  MACROBLOCKD *const xd = &x->e_mbd;
  MODE_INFO *const mi = xd->mi[0];
  struct macroblock_plane *const p = x->plane;
  const struct segmentation *const seg = &cm->seg;
  const int bw = num_8x8_blocks_wide_lookup[mi->sb_type];
  const int bh = num_8x8_blocks_high_lookup[mi->sb_type];
  const int x_mis = VPXMIN(bw, cm->mi_cols - mi_col);
  const int y_mis = VPXMIN(bh, cm->mi_rows - mi_row);

  *(xd->mi[0]) = ctx->mic;
  *(x->mbmi_ext) = ctx->mbmi_ext;

  if (seg->enabled && cpi->oxcf.aq_mode != NO_AQ) {
    // For in frame complexity AQ or variance AQ, copy segment_id from
    // segmentation_map.
    if (cpi->oxcf.aq_mode != CYCLIC_REFRESH_AQ) {
      const uint8_t *const map =
          seg->update_map ? cpi->segmentation_map : cm->last_frame_seg_map;
      mi->segment_id = get_segment_id(cm, map, bsize, mi_row, mi_col);
    } else {
      // Setting segmentation map for cyclic_refresh.
      vp9_cyclic_refresh_update_segment(cpi, mi, mi_row, mi_col, bsize,
                                        ctx->rate, ctx->dist, x->skip, p);
    }
    vp9_init_plane_quantizers(cpi, x);
  }

  if (is_inter_block(mi)) {
    vp9_update_mv_count(td);
    if (cm->interp_filter == SWITCHABLE) {
      const int pred_ctx = get_pred_context_switchable_interp(xd);
      ++td->counts->switchable_interp[pred_ctx][mi->interp_filter];
    }

    if (mi->sb_type < BLOCK_8X8) {
      mi->mv[0].as_int = mi->bmi[3].as_mv[0].as_int;
      mi->mv[1].as_int = mi->bmi[3].as_mv[1].as_int;
    }
  }

  if (cm->use_prev_frame_mvs || !cm->error_resilient_mode ||
      (cpi->svc.use_base_mv && cpi->svc.number_spatial_layers > 1 &&
       cpi->svc.spatial_layer_id != cpi->svc.number_spatial_layers - 1)) {
    MV_REF *const frame_mvs =
        cm->cur_frame->mvs + mi_row * cm->mi_cols + mi_col;
    int w, h;

    for (h = 0; h < y_mis; ++h) {
      MV_REF *const frame_mv = frame_mvs + h * cm->mi_cols;
      for (w = 0; w < x_mis; ++w) {
        MV_REF *const mv = frame_mv + w;
        mv->ref_frame[0] = mi->ref_frame[0];
        mv->ref_frame[1] = mi->ref_frame[1];
        mv->mv[0].as_int = mi->mv[0].as_int;
        mv->mv[1].as_int = mi->mv[1].as_int;
      }
    }
  }

  x->skip = ctx->skip;
  x->skip_txfm[0] = (mi->segment_id || xd->lossless) ? 0 : ctx->skip_txfm[0];
}

static void encode_b_rt(VP9_COMP *cpi, ThreadData *td,
                        const TileInfo *const tile, TOKENEXTRA **tp, int mi_row,
                        int mi_col, int output_enabled, BLOCK_SIZE bsize,
                        PICK_MODE_CONTEXT *ctx) {
  MACROBLOCK *const x = &td->mb;
  set_offsets(cpi, tile, x, mi_row, mi_col, bsize);
  update_state_rt(cpi, td, ctx, mi_row, mi_col, bsize);

  encode_superblock(cpi, td, tp, output_enabled, mi_row, mi_col, bsize, ctx);
  update_stats(&cpi->common, td);

  (*tp)->token = EOSB_TOKEN;
  (*tp)++;
}

static void encode_sb_rt(VP9_COMP *cpi, ThreadData *td,
                         const TileInfo *const tile, TOKENEXTRA **tp,
                         int mi_row, int mi_col, int output_enabled,
                         BLOCK_SIZE bsize, PC_TREE *pc_tree) {
  VP9_COMMON *const cm = &cpi->common;
  MACROBLOCK *const x = &td->mb;
  MACROBLOCKD *const xd = &x->e_mbd;

  const int bsl = b_width_log2_lookup[bsize], hbs = (1 << bsl) / 4;
  int ctx;
  PARTITION_TYPE partition;
  BLOCK_SIZE subsize;

  if (mi_row >= cm->mi_rows || mi_col >= cm->mi_cols) return;

  if (bsize >= BLOCK_8X8) {
    const int idx_str = xd->mi_stride * mi_row + mi_col;
    MODE_INFO **mi_8x8 = cm->mi_grid_visible + idx_str;
    ctx = partition_plane_context(xd, mi_row, mi_col, bsize);
    subsize = mi_8x8[0]->sb_type;
  } else {
    ctx = 0;
    subsize = BLOCK_4X4;
  }

  partition = partition_lookup[bsl][subsize];
  if (output_enabled && bsize != BLOCK_4X4)
    td->counts->partition[ctx][partition]++;

  switch (partition) {
    case PARTITION_NONE:
      encode_b_rt(cpi, td, tile, tp, mi_row, mi_col, output_enabled, subsize,
                  &pc_tree->none);
      break;
    case PARTITION_VERT:
      encode_b_rt(cpi, td, tile, tp, mi_row, mi_col, output_enabled, subsize,
                  &pc_tree->vertical[0]);
      if (mi_col + hbs < cm->mi_cols && bsize > BLOCK_8X8) {
        encode_b_rt(cpi, td, tile, tp, mi_row, mi_col + hbs, output_enabled,
                    subsize, &pc_tree->vertical[1]);
      }
      break;
    case PARTITION_HORZ:
      encode_b_rt(cpi, td, tile, tp, mi_row, mi_col, output_enabled, subsize,
                  &pc_tree->horizontal[0]);
      if (mi_row + hbs < cm->mi_rows && bsize > BLOCK_8X8) {
        encode_b_rt(cpi, td, tile, tp, mi_row + hbs, mi_col, output_enabled,
                    subsize, &pc_tree->horizontal[1]);
      }
      break;
    default:
      assert(partition == PARTITION_SPLIT);
      subsize = get_subsize(bsize, PARTITION_SPLIT);
      encode_sb_rt(cpi, td, tile, tp, mi_row, mi_col, output_enabled, subsize,
                   pc_tree->split[0]);
      encode_sb_rt(cpi, td, tile, tp, mi_row, mi_col + hbs, output_enabled,
                   subsize, pc_tree->split[1]);
      encode_sb_rt(cpi, td, tile, tp, mi_row + hbs, mi_col, output_enabled,
                   subsize, pc_tree->split[2]);
      encode_sb_rt(cpi, td, tile, tp, mi_row + hbs, mi_col + hbs,
                   output_enabled, subsize, pc_tree->split[3]);
      break;
  }

  if (partition != PARTITION_SPLIT || bsize == BLOCK_8X8)
    update_partition_context(xd, mi_row, mi_col, subsize, bsize);
}

static void rd_use_partition(VP9_COMP *cpi, ThreadData *td,
                             TileDataEnc *tile_data, MODE_INFO **mi_8x8,
                             TOKENEXTRA **tp, int mi_row, int mi_col,
                             BLOCK_SIZE bsize, int *rate, int64_t *dist,
                             int do_recon, PC_TREE *pc_tree) {
  VP9_COMMON *const cm = &cpi->common;
  TileInfo *const tile_info = &tile_data->tile_info;
  MACROBLOCK *const x = &td->mb;
  MACROBLOCKD *const xd = &x->e_mbd;
  const int mis = cm->mi_stride;
  const int bsl = b_width_log2_lookup[bsize];
  const int mi_step = num_4x4_blocks_wide_lookup[bsize] / 2;
  const int bss = (1 << bsl) / 4;
  int i, pl;
  PARTITION_TYPE partition = PARTITION_NONE;
  BLOCK_SIZE subsize;
  ENTROPY_CONTEXT l[16 * MAX_MB_PLANE], a[16 * MAX_MB_PLANE];
  PARTITION_CONTEXT sl[8], sa[8];
  RD_COST last_part_rdc, none_rdc, chosen_rdc;
  BLOCK_SIZE sub_subsize = BLOCK_4X4;
  int splits_below = 0;
  BLOCK_SIZE bs_type = mi_8x8[0]->sb_type;
  int do_partition_search = 1;
  PICK_MODE_CONTEXT *ctx = &pc_tree->none;

  if (mi_row >= cm->mi_rows || mi_col >= cm->mi_cols) return;

  assert(num_4x4_blocks_wide_lookup[bsize] ==
         num_4x4_blocks_high_lookup[bsize]);

  vp9_rd_cost_reset(&last_part_rdc);
  vp9_rd_cost_reset(&none_rdc);
  vp9_rd_cost_reset(&chosen_rdc);

  partition = partition_lookup[bsl][bs_type];
  subsize = get_subsize(bsize, partition);

  pc_tree->partitioning = partition;
  save_context(x, mi_row, mi_col, a, l, sa, sl, bsize);

  if (bsize == BLOCK_16X16 && cpi->oxcf.aq_mode != NO_AQ) {
    set_offsets(cpi, tile_info, x, mi_row, mi_col, bsize);
    x->mb_energy = vp9_block_energy(cpi, x, bsize);
  }

  if (do_partition_search &&
      cpi->sf.partition_search_type == SEARCH_PARTITION &&
      cpi->sf.adjust_partitioning_from_last_frame) {
    // Check if any of the sub blocks are further split.
    if (partition == PARTITION_SPLIT && subsize > BLOCK_8X8) {
      sub_subsize = get_subsize(subsize, PARTITION_SPLIT);
      splits_below = 1;
      for (i = 0; i < 4; i++) {
        int jj = i >> 1, ii = i & 0x01;
        MODE_INFO *this_mi = mi_8x8[jj * bss * mis + ii * bss];
        if (this_mi && this_mi->sb_type >= sub_subsize) {
          splits_below = 0;
        }
      }
    }

    // If partition is not none try none unless each of the 4 splits are split
    // even further..
    if (partition != PARTITION_NONE && !splits_below &&
        mi_row + (mi_step >> 1) < cm->mi_rows &&
        mi_col + (mi_step >> 1) < cm->mi_cols) {
      pc_tree->partitioning = PARTITION_NONE;
      rd_pick_sb_modes(cpi, tile_data, x, mi_row, mi_col, &none_rdc, bsize, ctx,
                       INT64_MAX);

      pl = partition_plane_context(xd, mi_row, mi_col, bsize);

      if (none_rdc.rate < INT_MAX) {
        none_rdc.rate += cpi->partition_cost[pl][PARTITION_NONE];
        none_rdc.rdcost =
            RDCOST(x->rdmult, x->rddiv, none_rdc.rate, none_rdc.dist);
      }

      restore_context(x, mi_row, mi_col, a, l, sa, sl, bsize);
      mi_8x8[0]->sb_type = bs_type;
      pc_tree->partitioning = partition;
    }
  }

  switch (partition) {
    case PARTITION_NONE:
      rd_pick_sb_modes(cpi, tile_data, x, mi_row, mi_col, &last_part_rdc, bsize,
                       ctx, INT64_MAX);
      break;
    case PARTITION_HORZ:
      pc_tree->horizontal[0].skip_ref_frame_mask = 0;
      rd_pick_sb_modes(cpi, tile_data, x, mi_row, mi_col, &last_part_rdc,
                       subsize, &pc_tree->horizontal[0], INT64_MAX);
      if (last_part_rdc.rate != INT_MAX && bsize >= BLOCK_8X8 &&
          mi_row + (mi_step >> 1) < cm->mi_rows) {
        RD_COST tmp_rdc;
        PICK_MODE_CONTEXT *ctx = &pc_tree->horizontal[0];
        vp9_rd_cost_init(&tmp_rdc);
        update_state(cpi, td, ctx, mi_row, mi_col, subsize, 0);
        encode_superblock(cpi, td, tp, 0, mi_row, mi_col, subsize, ctx);
        pc_tree->horizontal[1].skip_ref_frame_mask = 0;
        rd_pick_sb_modes(cpi, tile_data, x, mi_row + (mi_step >> 1), mi_col,
                         &tmp_rdc, subsize, &pc_tree->horizontal[1], INT64_MAX);
        if (tmp_rdc.rate == INT_MAX || tmp_rdc.dist == INT64_MAX) {
          vp9_rd_cost_reset(&last_part_rdc);
          break;
        }
        last_part_rdc.rate += tmp_rdc.rate;
        last_part_rdc.dist += tmp_rdc.dist;
        last_part_rdc.rdcost += tmp_rdc.rdcost;
      }
      break;
    case PARTITION_VERT:
      pc_tree->vertical[0].skip_ref_frame_mask = 0;
      rd_pick_sb_modes(cpi, tile_data, x, mi_row, mi_col, &last_part_rdc,
                       subsize, &pc_tree->vertical[0], INT64_MAX);
      if (last_part_rdc.rate != INT_MAX && bsize >= BLOCK_8X8 &&
          mi_col + (mi_step >> 1) < cm->mi_cols) {
        RD_COST tmp_rdc;
        PICK_MODE_CONTEXT *ctx = &pc_tree->vertical[0];
        vp9_rd_cost_init(&tmp_rdc);
        update_state(cpi, td, ctx, mi_row, mi_col, subsize, 0);
        encode_superblock(cpi, td, tp, 0, mi_row, mi_col, subsize, ctx);
        pc_tree->vertical[bsize > BLOCK_8X8].skip_ref_frame_mask = 0;
        rd_pick_sb_modes(cpi, tile_data, x, mi_row, mi_col + (mi_step >> 1),
                         &tmp_rdc, subsize,
                         &pc_tree->vertical[bsize > BLOCK_8X8], INT64_MAX);
        if (tmp_rdc.rate == INT_MAX || tmp_rdc.dist == INT64_MAX) {
          vp9_rd_cost_reset(&last_part_rdc);
          break;
        }
        last_part_rdc.rate += tmp_rdc.rate;
        last_part_rdc.dist += tmp_rdc.dist;
        last_part_rdc.rdcost += tmp_rdc.rdcost;
      }
      break;
    default:
      assert(partition == PARTITION_SPLIT);
      if (bsize == BLOCK_8X8) {
        rd_pick_sb_modes(cpi, tile_data, x, mi_row, mi_col, &last_part_rdc,
                         subsize, pc_tree->leaf_split[0], INT64_MAX);
        break;
      }
      last_part_rdc.rate = 0;
      last_part_rdc.dist = 0;
      last_part_rdc.rdcost = 0;
      for (i = 0; i < 4; i++) {
        int x_idx = (i & 1) * (mi_step >> 1);
        int y_idx = (i >> 1) * (mi_step >> 1);
        int jj = i >> 1, ii = i & 0x01;
        RD_COST tmp_rdc;
        if ((mi_row + y_idx >= cm->mi_rows) || (mi_col + x_idx >= cm->mi_cols))
          continue;

        vp9_rd_cost_init(&tmp_rdc);
        rd_use_partition(cpi, td, tile_data, mi_8x8 + jj * bss * mis + ii * bss,
                         tp, mi_row + y_idx, mi_col + x_idx, subsize,
                         &tmp_rdc.rate, &tmp_rdc.dist, i != 3,
                         pc_tree->split[i]);
        if (tmp_rdc.rate == INT_MAX || tmp_rdc.dist == INT64_MAX) {
          vp9_rd_cost_reset(&last_part_rdc);
          break;
        }
        last_part_rdc.rate += tmp_rdc.rate;
        last_part_rdc.dist += tmp_rdc.dist;
      }
      break;
  }

  pl = partition_plane_context(xd, mi_row, mi_col, bsize);
  if (last_part_rdc.rate < INT_MAX) {
    last_part_rdc.rate += cpi->partition_cost[pl][partition];
    last_part_rdc.rdcost =
        RDCOST(x->rdmult, x->rddiv, last_part_rdc.rate, last_part_rdc.dist);
  }

  if (do_partition_search && cpi->sf.adjust_partitioning_from_last_frame &&
      cpi->sf.partition_search_type == SEARCH_PARTITION &&
      partition != PARTITION_SPLIT && bsize > BLOCK_8X8 &&
      (mi_row + mi_step < cm->mi_rows ||
       mi_row + (mi_step >> 1) == cm->mi_rows) &&
      (mi_col + mi_step < cm->mi_cols ||
       mi_col + (mi_step >> 1) == cm->mi_cols)) {
    BLOCK_SIZE split_subsize = get_subsize(bsize, PARTITION_SPLIT);
    chosen_rdc.rate = 0;
    chosen_rdc.dist = 0;
    restore_context(x, mi_row, mi_col, a, l, sa, sl, bsize);
    pc_tree->partitioning = PARTITION_SPLIT;

    // Split partition.
    for (i = 0; i < 4; i++) {
      int x_idx = (i & 1) * (mi_step >> 1);
      int y_idx = (i >> 1) * (mi_step >> 1);
      RD_COST tmp_rdc;
      ENTROPY_CONTEXT l[16 * MAX_MB_PLANE], a[16 * MAX_MB_PLANE];
      PARTITION_CONTEXT sl[8], sa[8];

      if ((mi_row + y_idx >= cm->mi_rows) || (mi_col + x_idx >= cm->mi_cols))
        continue;

      save_context(x, mi_row, mi_col, a, l, sa, sl, bsize);
      pc_tree->split[i]->partitioning = PARTITION_NONE;
      rd_pick_sb_modes(cpi, tile_data, x, mi_row + y_idx, mi_col + x_idx,
                       &tmp_rdc, split_subsize, &pc_tree->split[i]->none,
                       INT64_MAX);

      restore_context(x, mi_row, mi_col, a, l, sa, sl, bsize);

      if (tmp_rdc.rate == INT_MAX || tmp_rdc.dist == INT64_MAX) {
        vp9_rd_cost_reset(&chosen_rdc);
        break;
      }

      chosen_rdc.rate += tmp_rdc.rate;
      chosen_rdc.dist += tmp_rdc.dist;

      if (i != 3)
        encode_sb(cpi, td, tile_info, tp, mi_row + y_idx, mi_col + x_idx, 0,
                  split_subsize, pc_tree->split[i]);

      pl = partition_plane_context(xd, mi_row + y_idx, mi_col + x_idx,
                                   split_subsize);
      chosen_rdc.rate += cpi->partition_cost[pl][PARTITION_NONE];
    }
    pl = partition_plane_context(xd, mi_row, mi_col, bsize);
    if (chosen_rdc.rate < INT_MAX) {
      chosen_rdc.rate += cpi->partition_cost[pl][PARTITION_SPLIT];
      chosen_rdc.rdcost =
          RDCOST(x->rdmult, x->rddiv, chosen_rdc.rate, chosen_rdc.dist);
    }
  }

  // If last_part is better set the partitioning to that.
  if (last_part_rdc.rdcost < chosen_rdc.rdcost) {
    mi_8x8[0]->sb_type = bsize;
    if (bsize >= BLOCK_8X8) pc_tree->partitioning = partition;
    chosen_rdc = last_part_rdc;
  }
  // If none was better set the partitioning to that.
  if (none_rdc.rdcost < chosen_rdc.rdcost) {
    if (bsize >= BLOCK_8X8) pc_tree->partitioning = PARTITION_NONE;
    chosen_rdc = none_rdc;
  }

  restore_context(x, mi_row, mi_col, a, l, sa, sl, bsize);

  // We must have chosen a partitioning and encoding or we'll fail later on.
  // No other opportunities for success.
  if (bsize == BLOCK_64X64)
    assert(chosen_rdc.rate < INT_MAX && chosen_rdc.dist < INT64_MAX);

  if (do_recon) {
    int output_enabled = (bsize == BLOCK_64X64);
    encode_sb(cpi, td, tile_info, tp, mi_row, mi_col, output_enabled, bsize,
              pc_tree);
  }

  *rate = chosen_rdc.rate;
  *dist = chosen_rdc.dist;
}

static const BLOCK_SIZE min_partition_size[BLOCK_SIZES] = {
  BLOCK_4X4,   BLOCK_4X4,   BLOCK_4X4,  BLOCK_4X4, BLOCK_4X4,
  BLOCK_4X4,   BLOCK_8X8,   BLOCK_8X8,  BLOCK_8X8, BLOCK_16X16,
  BLOCK_16X16, BLOCK_16X16, BLOCK_16X16
};

static const BLOCK_SIZE max_partition_size[BLOCK_SIZES] = {
  BLOCK_8X8,   BLOCK_16X16, BLOCK_16X16, BLOCK_16X16, BLOCK_32X32,
  BLOCK_32X32, BLOCK_32X32, BLOCK_64X64, BLOCK_64X64, BLOCK_64X64,
  BLOCK_64X64, BLOCK_64X64, BLOCK_64X64
};

// Look at all the mode_info entries for blocks that are part of this
// partition and find the min and max values for sb_type.
// At the moment this is designed to work on a 64x64 SB but could be
// adjusted to use a size parameter.
//
// The min and max are assumed to have been initialized prior to calling this
// function so repeat calls can accumulate a min and max of more than one sb64.
static void get_sb_partition_size_range(MACROBLOCKD *xd, MODE_INFO **mi_8x8,
                                        BLOCK_SIZE *min_block_size,
                                        BLOCK_SIZE *max_block_size,
                                        int bs_hist[BLOCK_SIZES]) {
  int sb_width_in_blocks = MI_BLOCK_SIZE;
  int sb_height_in_blocks = MI_BLOCK_SIZE;
  int i, j;
  int index = 0;

  // Check the sb_type for each block that belongs to this region.
  for (i = 0; i < sb_height_in_blocks; ++i) {
    for (j = 0; j < sb_width_in_blocks; ++j) {
      MODE_INFO *mi = mi_8x8[index + j];
      BLOCK_SIZE sb_type = mi ? mi->sb_type : 0;
      bs_hist[sb_type]++;
      *min_block_size = VPXMIN(*min_block_size, sb_type);
      *max_block_size = VPXMAX(*max_block_size, sb_type);
    }
    index += xd->mi_stride;
  }
}

// Next square block size less or equal than current block size.
static const BLOCK_SIZE next_square_size[BLOCK_SIZES] = {
  BLOCK_4X4,   BLOCK_4X4,   BLOCK_4X4,   BLOCK_8X8,   BLOCK_8X8,
  BLOCK_8X8,   BLOCK_16X16, BLOCK_16X16, BLOCK_16X16, BLOCK_32X32,
  BLOCK_32X32, BLOCK_32X32, BLOCK_64X64
};

// Look at neighboring blocks and set a min and max partition size based on
// what they chose.
static void rd_auto_partition_range(VP9_COMP *cpi, const TileInfo *const tile,
                                    MACROBLOCKD *const xd, int mi_row,
                                    int mi_col, BLOCK_SIZE *min_block_size,
                                    BLOCK_SIZE *max_block_size) {
  VP9_COMMON *const cm = &cpi->common;
  MODE_INFO **mi = xd->mi;
  const int left_in_image = !!xd->left_mi;
  const int above_in_image = !!xd->above_mi;
  const int row8x8_remaining = tile->mi_row_end - mi_row;
  const int col8x8_remaining = tile->mi_col_end - mi_col;
  int bh, bw;
  BLOCK_SIZE min_size = BLOCK_4X4;
  BLOCK_SIZE max_size = BLOCK_64X64;
  int bs_hist[BLOCK_SIZES] = { 0 };

  // Trap case where we do not have a prediction.
  if (left_in_image || above_in_image || cm->frame_type != KEY_FRAME) {
    // Default "min to max" and "max to min"
    min_size = BLOCK_64X64;
    max_size = BLOCK_4X4;

    // NOTE: each call to get_sb_partition_size_range() uses the previous
    // passed in values for min and max as a starting point.
    // Find the min and max partition used in previous frame at this location
    if (cm->frame_type != KEY_FRAME) {
      MODE_INFO **prev_mi =
          &cm->prev_mi_grid_visible[mi_row * xd->mi_stride + mi_col];
      get_sb_partition_size_range(xd, prev_mi, &min_size, &max_size, bs_hist);
    }
    // Find the min and max partition sizes used in the left SB64
    if (left_in_image) {
      MODE_INFO **left_sb64_mi = &mi[-MI_BLOCK_SIZE];
      get_sb_partition_size_range(xd, left_sb64_mi, &min_size, &max_size,
                                  bs_hist);
    }
    // Find the min and max partition sizes used in the above SB64.
    if (above_in_image) {
      MODE_INFO **above_sb64_mi = &mi[-xd->mi_stride * MI_BLOCK_SIZE];
      get_sb_partition_size_range(xd, above_sb64_mi, &min_size, &max_size,
                                  bs_hist);
    }

    // Adjust observed min and max for "relaxed" auto partition case.
    if (cpi->sf.auto_min_max_partition_size == RELAXED_NEIGHBORING_MIN_MAX) {
      min_size = min_partition_size[min_size];
      max_size = max_partition_size[max_size];
    }
  }

  // Check border cases where max and min from neighbors may not be legal.
  max_size = find_partition_size(max_size, row8x8_remaining, col8x8_remaining,
                                 &bh, &bw);
  // Test for blocks at the edge of the active image.
  // This may be the actual edge of the image or where there are formatting
  // bars.
  if (vp9_active_edge_sb(cpi, mi_row, mi_col)) {
    min_size = BLOCK_4X4;
  } else {
    min_size =
        VPXMIN(cpi->sf.rd_auto_partition_min_limit, VPXMIN(min_size, max_size));
  }

  // When use_square_partition_only is true, make sure at least one square
  // partition is allowed by selecting the next smaller square size as
  // *min_block_size.
  if (cpi->sf.use_square_partition_only &&
      next_square_size[max_size] < min_size) {
    min_size = next_square_size[max_size];
  }

  *min_block_size = min_size;
  *max_block_size = max_size;
}

// TODO(jingning) refactor functions setting partition search range
static void set_partition_range(VP9_COMMON *cm, MACROBLOCKD *xd, int mi_row,
                                int mi_col, BLOCK_SIZE bsize,
                                BLOCK_SIZE *min_bs, BLOCK_SIZE *max_bs) {
  int mi_width = num_8x8_blocks_wide_lookup[bsize];
  int mi_height = num_8x8_blocks_high_lookup[bsize];
  int idx, idy;

  MODE_INFO *mi;
  const int idx_str = cm->mi_stride * mi_row + mi_col;
  MODE_INFO **prev_mi = &cm->prev_mi_grid_visible[idx_str];
  BLOCK_SIZE bs, min_size, max_size;

  min_size = BLOCK_64X64;
  max_size = BLOCK_4X4;

  if (prev_mi) {
    for (idy = 0; idy < mi_height; ++idy) {
      for (idx = 0; idx < mi_width; ++idx) {
        mi = prev_mi[idy * cm->mi_stride + idx];
        bs = mi ? mi->sb_type : bsize;
        min_size = VPXMIN(min_size, bs);
        max_size = VPXMAX(max_size, bs);
      }
    }
  }

  if (xd->left_mi) {
    for (idy = 0; idy < mi_height; ++idy) {
      mi = xd->mi[idy * cm->mi_stride - 1];
      bs = mi ? mi->sb_type : bsize;
      min_size = VPXMIN(min_size, bs);
      max_size = VPXMAX(max_size, bs);
    }
  }

  if (xd->above_mi) {
    for (idx = 0; idx < mi_width; ++idx) {
      mi = xd->mi[idx - cm->mi_stride];
      bs = mi ? mi->sb_type : bsize;
      min_size = VPXMIN(min_size, bs);
      max_size = VPXMAX(max_size, bs);
    }
  }

  if (min_size == max_size) {
    min_size = min_partition_size[min_size];
    max_size = max_partition_size[max_size];
  }

  *min_bs = min_size;
  *max_bs = max_size;
}

static INLINE void store_pred_mv(MACROBLOCK *x, PICK_MODE_CONTEXT *ctx) {
  memcpy(ctx->pred_mv, x->pred_mv, sizeof(x->pred_mv));
}

static INLINE void load_pred_mv(MACROBLOCK *x, PICK_MODE_CONTEXT *ctx) {
  memcpy(x->pred_mv, ctx->pred_mv, sizeof(x->pred_mv));
}

#if CONFIG_FP_MB_STATS
const int num_16x16_blocks_wide_lookup[BLOCK_SIZES] = { 1, 1, 1, 1, 1, 1, 1,
                                                        1, 2, 2, 2, 4, 4 };
const int num_16x16_blocks_high_lookup[BLOCK_SIZES] = { 1, 1, 1, 1, 1, 1, 1,
                                                        2, 1, 2, 4, 2, 4 };
const int qindex_skip_threshold_lookup[BLOCK_SIZES] = {
  0, 10, 10, 30, 40, 40, 60, 80, 80, 90, 100, 100, 120
};
const int qindex_split_threshold_lookup[BLOCK_SIZES] = {
  0, 3, 3, 7, 15, 15, 30, 40, 40, 60, 80, 80, 120
};
const int complexity_16x16_blocks_threshold[BLOCK_SIZES] = {
  1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 4, 4, 6
};

typedef enum {
  MV_ZERO = 0,
  MV_LEFT = 1,
  MV_UP = 2,
  MV_RIGHT = 3,
  MV_DOWN = 4,
  MV_INVALID
} MOTION_DIRECTION;

static INLINE MOTION_DIRECTION get_motion_direction_fp(uint8_t fp_byte) {
  if (fp_byte & FPMB_MOTION_ZERO_MASK) {
    return MV_ZERO;
  } else if (fp_byte & FPMB_MOTION_LEFT_MASK) {
    return MV_LEFT;
  } else if (fp_byte & FPMB_MOTION_RIGHT_MASK) {
    return MV_RIGHT;
  } else if (fp_byte & FPMB_MOTION_UP_MASK) {
    return MV_UP;
  } else {
    return MV_DOWN;
  }
}

static INLINE int get_motion_inconsistency(MOTION_DIRECTION this_mv,
                                           MOTION_DIRECTION that_mv) {
  if (this_mv == that_mv) {
    return 0;
  } else {
    return abs(this_mv - that_mv) == 2 ? 2 : 1;
  }
}
#endif

#define NN_MAX_HIDDEN_LAYERS 10
#define NN_MAX_NODES_PER_LAYER 128

// Neural net model config.
typedef struct {
  int num_inputs;         // Number of input nodes, i.e. features.
  int num_outputs;        // Number of output nodes.
  int num_hidden_layers;  // Number of hidden layers, maximum 10.
  // Number of nodes for each hidden layer.
  int num_hidden_nodes[NN_MAX_HIDDEN_LAYERS];
  // Weight parameters, indexed by layer.
  const float *weights[NN_MAX_HIDDEN_LAYERS + 1];
  // Bias parameters, indexed by layer.
  const float *bias[NN_MAX_HIDDEN_LAYERS + 1];
} NN_CONFIG;

// Calculate prediction based on the given input features and neural net config.
// Assume there are no more than NN_MAX_NODES_PER_LAYER nodes in each hidden
// layer.
static void nn_predict(const float *features, const NN_CONFIG *nn_config,
                       float *output) {
  int num_input_nodes = nn_config->num_inputs;
  int buf_index = 0;
  float buf[2][NN_MAX_NODES_PER_LAYER];
  const float *input_nodes = features;

  // Propagate hidden layers.
  const int num_layers = nn_config->num_hidden_layers;
  int layer, node, i;
  assert(num_layers <= NN_MAX_HIDDEN_LAYERS);
  for (layer = 0; layer < num_layers; ++layer) {
    const float *weights = nn_config->weights[layer];
    const float *bias = nn_config->bias[layer];
    float *output_nodes = buf[buf_index];
    const int num_output_nodes = nn_config->num_hidden_nodes[layer];
    assert(num_output_nodes < NN_MAX_NODES_PER_LAYER);
    for (node = 0; node < num_output_nodes; ++node) {
      float val = 0.0f;
      for (i = 0; i < num_input_nodes; ++i) val += weights[i] * input_nodes[i];
      val += bias[node];
      // ReLU as activation function.
      val = VPXMAX(val, 0.0f);
      output_nodes[node] = val;
      weights += num_input_nodes;
    }
    num_input_nodes = num_output_nodes;
    input_nodes = output_nodes;
    buf_index = 1 - buf_index;
  }

  // Final output layer.
  {
    const float *weights = nn_config->weights[num_layers];
    for (node = 0; node < nn_config->num_outputs; ++node) {
      const float *bias = nn_config->bias[num_layers];
      float val = 0.0f;
      for (i = 0; i < num_input_nodes; ++i) val += weights[i] * input_nodes[i];
      output[node] = val + bias[node];
      weights += num_input_nodes;
    }
  }
}

static const float partition_nn_weights_64x64_layer0[7 * 8] = {
  -3.571348f, 0.014835f,  -3.255393f, -0.098090f, -0.013120f, 0.000221f,
  0.056273f,  0.190179f,  -0.268130f, -1.828242f, -0.010655f, 0.937244f,
  -0.435120f, 0.512125f,  1.610679f,  0.190816f,  -0.799075f, -0.377348f,
  -0.144232f, 0.614383f,  -0.980388f, 1.754150f,  -0.185603f, -0.061854f,
  -0.807172f, 1.240177f,  1.419531f,  -0.438544f, -5.980774f, 0.139045f,
  -0.032359f, -0.068887f, -1.237918f, 0.115706f,  0.003164f,  2.924212f,
  1.246838f,  -0.035833f, 0.810011f,  -0.805894f, 0.010966f,  0.076463f,
  -4.226380f, -2.437764f, -0.010619f, -0.020935f, -0.451494f, 0.300079f,
  -0.168961f, -3.326450f, -2.731094f, 0.002518f,  0.018840f,  -1.656815f,
  0.068039f,  0.010586f,
};

static const float partition_nn_bias_64x64_layer0[8] = {
  -3.469882f, 0.683989f, 0.194010f,  0.313782f,
  -3.153335f, 2.245849f, -1.946190f, -3.740020f,
};

static const float partition_nn_weights_64x64_layer1[8] = {
  -8.058566f, 0.108306f, -0.280620f, -0.818823f,
  -6.445117f, 0.865364f, -1.127127f, -8.808660f,
};

static const float partition_nn_bias_64x64_layer1[1] = {
  6.46909416f,
};

static const NN_CONFIG partition_nnconfig_64x64 = {
  7,  // num_inputs
  1,  // num_outputs
  1,  // num_hidden_layers
  {
      8,
  },  // num_hidden_nodes
  {
      partition_nn_weights_64x64_layer0,
      partition_nn_weights_64x64_layer1,
  },
  {
      partition_nn_bias_64x64_layer0,
      partition_nn_bias_64x64_layer1,
  },
};

static const float partition_nn_weights_32x32_layer0[7 * 8] = {
  -0.295437f, -4.002648f, -0.205399f, -0.060919f, 0.708037f,  0.027221f,
  -0.039137f, -0.907724f, -3.151662f, 0.007106f,  0.018726f,  -0.534928f,
  0.022744f,  0.000159f,  -1.717189f, -3.229031f, -0.027311f, 0.269863f,
  -0.400747f, -0.394366f, -0.108878f, 0.603027f,  0.455369f,  -0.197170f,
  1.241746f,  -1.347820f, -0.575636f, -0.462879f, -2.296426f, 0.196696f,
  -0.138347f, -0.030754f, -0.200774f, 0.453795f,  0.055625f,  -3.163116f,
  -0.091003f, -0.027028f, -0.042984f, -0.605185f, 0.143240f,  -0.036439f,
  -0.801228f, 0.313409f,  -0.159942f, 0.031267f,  0.886454f,  -1.531644f,
  -0.089655f, 0.037683f,  -0.163441f, -0.130454f, -0.058344f, 0.060011f,
  0.275387f,  1.552226f,
};

static const float partition_nn_bias_32x32_layer0[8] = {
  -0.838372f, -2.609089f, -0.055763f, 1.329485f,
  -1.297638f, -2.636622f, -0.826909f, 1.012644f,
};

static const float partition_nn_weights_32x32_layer1[8] = {
  -1.792632f, -7.322353f, -0.683386f, 0.676564f,
  -1.488118f, -7.527719f, 1.240163f,  0.614309f,
};

static const float partition_nn_bias_32x32_layer1[1] = {
  4.97422546f,
};

static const NN_CONFIG partition_nnconfig_32x32 = {
  7,  // num_inputs
  1,  // num_outputs
  1,  // num_hidden_layers
  {
      8,
  },  // num_hidden_nodes
  {
      partition_nn_weights_32x32_layer0,
      partition_nn_weights_32x32_layer1,
  },
  {
      partition_nn_bias_32x32_layer0,
      partition_nn_bias_32x32_layer1,
  },
};

static const float partition_nn_weights_16x16_layer0[7 * 8] = {
  -1.717673f, -4.718130f, -0.125725f, -0.183427f, -0.511764f, 0.035328f,
  0.130891f,  -3.096753f, 0.174968f,  -0.188769f, -0.640796f, 1.305661f,
  1.700638f,  -0.073806f, -4.006781f, -1.630999f, -0.064863f, -0.086410f,
  -0.148617f, 0.172733f,  -0.018619f, 2.152595f,  0.778405f,  -0.156455f,
  0.612995f,  -0.467878f, 0.152022f,  -0.236183f, 0.339635f,  -0.087119f,
  -3.196610f, -1.080401f, -0.637704f, -0.059974f, 1.706298f,  -0.793705f,
  -6.399260f, 0.010624f,  -0.064199f, -0.650621f, 0.338087f,  -0.001531f,
  1.023655f,  -3.700272f, -0.055281f, -0.386884f, 0.375504f,  -0.898678f,
  0.281156f,  -0.314611f, 0.863354f,  -0.040582f, -0.145019f, 0.029329f,
  -2.197880f, -0.108733f,
};

static const float partition_nn_bias_16x16_layer0[8] = {
  0.411516f,  -2.143737f, -3.693192f, 2.123142f,
  -1.356910f, -3.561016f, -0.765045f, -2.417082f,
};

static const float partition_nn_weights_16x16_layer1[8] = {
  -0.619755f, -2.202391f, -4.337171f, 0.611319f,
  0.377677f,  -4.998723f, -1.052235f, 1.949922f,
};

static const float partition_nn_bias_16x16_layer1[1] = {
  3.20981717f,
};

static const NN_CONFIG partition_nnconfig_16x16 = {
  7,  // num_inputs
  1,  // num_outputs
  1,  // num_hidden_layers
  {
      8,
  },  // num_hidden_nodes
  {
      partition_nn_weights_16x16_layer0,
      partition_nn_weights_16x16_layer1,
  },
  {
      partition_nn_bias_16x16_layer0,
      partition_nn_bias_16x16_layer1,
  },
};

static const float partition_feature_mean[24] = {
  303501.697372f, 3042630.372158f, 24.694696f, 1.392182f,
  689.413511f,    162.027012f,     1.478213f,  0.0,
  135382.260230f, 912738.513263f,  28.845217f, 1.515230f,
  544.158492f,    131.807995f,     1.436863f,  0.0f,
  43682.377587f,  208131.711766f,  28.084737f, 1.356677f,
  138.254122f,    119.522553f,     1.252322f,  0.0f,
};

static const float partition_feature_std[24] = {
  673689.212982f, 5996652.516628f, 0.024449f, 1.989792f,
  985.880847f,    0.014638f,       2.001898f, 0.0f,
  208798.775332f, 1812548.443284f, 0.018693f, 1.838009f,
  396.986910f,    0.015657f,       1.332541f, 0.0f,
  55888.847031f,  448587.962714f,  0.017900f, 1.904776f,
  98.652832f,     0.016598f,       1.320992f, 0.0f,
};

// Error tolerance: 0.01%-0.0.05%-0.1%
static const float partition_linear_weights[24] = {
  0.111736f, 0.289977f, 0.042219f, 0.204765f, 0.120410f, -0.143863f,
  0.282376f, 0.847811f, 0.637161f, 0.131570f, 0.018636f, 0.202134f,
  0.112797f, 0.028162f, 0.182450f, 1.124367f, 0.386133f, 0.083700f,
  0.050028f, 0.150873f, 0.061119f, 0.109318f, 0.127255f, 0.625211f,
};

// Machine-learning based partition search early termination.
// Return 1 to skip split and rect partitions.
static int ml_pruning_partition(VP9_COMMON *const cm, MACROBLOCKD *const xd,
                                PICK_MODE_CONTEXT *ctx, int mi_row, int mi_col,
                                BLOCK_SIZE bsize) {
  const int mag_mv =
      abs(ctx->mic.mv[0].as_mv.col) + abs(ctx->mic.mv[0].as_mv.row);
  const int left_in_image = !!xd->left_mi;
  const int above_in_image = !!xd->above_mi;
  MODE_INFO **prev_mi =
      &cm->prev_mi_grid_visible[mi_col + cm->mi_stride * mi_row];
  int above_par = 0;  // above_partitioning
  int left_par = 0;   // left_partitioning
  int last_par = 0;   // last_partitioning
  int offset = 0;
  int i;
  BLOCK_SIZE context_size;
  const NN_CONFIG *nn_config = NULL;
  const float *mean, *sd, *linear_weights;
  float nn_score, linear_score;
  float features[7];

  assert(b_width_log2_lookup[bsize] == b_height_log2_lookup[bsize]);
  vpx_clear_system_state();

  switch (bsize) {
    case BLOCK_64X64:
      offset = 0;
      nn_config = &partition_nnconfig_64x64;
      break;
    case BLOCK_32X32:
      offset = 8;
      nn_config = &partition_nnconfig_32x32;
      break;
    case BLOCK_16X16:
      offset = 16;
      nn_config = &partition_nnconfig_16x16;
      break;
    default: assert(0 && "Unexpected block size."); return 0;
  }

  if (above_in_image) {
    context_size = xd->above_mi->sb_type;
    if (context_size < bsize)
      above_par = 2;
    else if (context_size == bsize)
      above_par = 1;
  }

  if (left_in_image) {
    context_size = xd->left_mi->sb_type;
    if (context_size < bsize)
      left_par = 2;
    else if (context_size == bsize)
      left_par = 1;
  }

  if (prev_mi) {
    context_size = prev_mi[0]->sb_type;
    if (context_size < bsize)
      last_par = 2;
    else if (context_size == bsize)
      last_par = 1;
  }

  mean = &partition_feature_mean[offset];
  sd = &partition_feature_std[offset];
  features[0] = ((float)ctx->rate - mean[0]) / sd[0];
  features[1] = ((float)ctx->dist - mean[1]) / sd[1];
  features[2] = ((float)mag_mv / 2 - mean[2]) * sd[2];
  features[3] = ((float)(left_par + above_par) / 2 - mean[3]) * sd[3];
  features[4] = ((float)ctx->sum_y_eobs - mean[4]) / sd[4];
  features[5] = ((float)cm->base_qindex - mean[5]) * sd[5];
  features[6] = ((float)last_par - mean[6]) * sd[6];

  // Predict using linear model.
  linear_weights = &partition_linear_weights[offset];
  linear_score = linear_weights[7];
  for (i = 0; i < 7; ++i) linear_score += linear_weights[i] * features[i];
  if (linear_score > 0.1f) return 0;

  // Predict using neural net model.
  nn_predict(features, nn_config, &nn_score);

  if (linear_score < -0.0f && nn_score < 0.1f) return 1;
  if (nn_score < -0.0f && linear_score < 0.1f) return 1;
  return 0;
}

#define FEATURES 4
#define Q_CTX 3
#define RESOLUTION_CTX 2
static const float partition_breakout_weights_64[RESOLUTION_CTX][Q_CTX]
                                                [FEATURES + 1] = {
                                                  {
                                                      {
                                                          -0.016673f,
                                                          -0.001025f,
                                                          -0.000032f,
                                                          0.000833f,
                                                          1.94261885f - 2.1f,
                                                      },
                                                      {
                                                          -0.160867f,
                                                          -0.002101f,
                                                          0.000011f,
                                                          0.002448f,
                                                          1.65738142f - 2.5f,
                                                      },
                                                      {
                                                          -0.628934f,
                                                          -0.011459f,
                                                          -0.000009f,
                                                          0.013833f,
                                                          1.47982645f - 1.6f,
                                                      },
                                                  },
                                                  {
                                                      {
                                                          -0.064309f,
                                                          -0.006121f,
                                                          0.000232f,
                                                          0.005778f,
                                                          0.7989465f - 5.0f,
                                                      },
                                                      {
                                                          -0.314957f,
                                                          -0.009346f,
                                                          -0.000225f,
                                                          0.010072f,
                                                          2.80695581f - 5.5f,
                                                      },
                                                      {
                                                          -0.635535f,
                                                          -0.015135f,
                                                          0.000091f,
                                                          0.015247f,
                                                          2.90381241f - 5.0f,
                                                      },
                                                  },
                                                };

static const float partition_breakout_weights_32[RESOLUTION_CTX][Q_CTX]
                                                [FEATURES + 1] = {
                                                  {
                                                      {
                                                          -0.010554f,
                                                          -0.003081f,
                                                          -0.000134f,
                                                          0.004491f,
                                                          1.68445992f - 3.5f,
                                                      },
                                                      {
                                                          -0.051489f,
                                                          -0.007609f,
                                                          0.000016f,
                                                          0.009792f,
                                                          1.28089404f - 2.5f,
                                                      },
                                                      {
                                                          -0.163097f,
                                                          -0.013081f,
                                                          0.000022f,
                                                          0.019006f,
                                                          1.36129403f - 3.2f,
                                                      },
                                                  },
                                                  {
                                                      {
                                                          -0.024629f,
                                                          -0.006492f,
                                                          -0.000254f,
                                                          0.004895f,
                                                          1.27919173f - 4.5f,
                                                      },
                                                      {
                                                          -0.083936f,
                                                          -0.009827f,
                                                          -0.000200f,
                                                          0.010399f,
                                                          2.73731065f - 4.5f,
                                                      },
                                                      {
                                                          -0.279052f,
                                                          -0.013334f,
                                                          0.000289f,
                                                          0.023203f,
                                                          2.43595719f - 3.5f,
                                                      },
                                                  },
                                                };

static const float partition_breakout_weights_16[RESOLUTION_CTX][Q_CTX]
                                                [FEATURES + 1] = {
                                                  {
                                                      {
                                                          -0.013154f,
                                                          -0.002404f,
                                                          -0.000977f,
                                                          0.008450f,
                                                          2.57404566f - 5.5f,
                                                      },
                                                      {
                                                          -0.019146f,
                                                          -0.004018f,
                                                          0.000064f,
                                                          0.008187f,
                                                          2.15043926f - 2.5f,
                                                      },
                                                      {
                                                          -0.075755f,
                                                          -0.010858f,
                                                          0.000030f,
                                                          0.024505f,
                                                          2.06848121f - 2.5f,
                                                      },
                                                  },
                                                  {
                                                      {
                                                          -0.007636f,
                                                          -0.002751f,
                                                          -0.000682f,
                                                          0.005968f,
                                                          0.19225763f - 4.5f,
                                                      },
                                                      {
                                                          -0.047306f,
                                                          -0.009113f,
                                                          -0.000518f,
                                                          0.016007f,
                                                          2.61068869f - 4.0f,
                                                      },
                                                      {
                                                          -0.069336f,
                                                          -0.010448f,
                                                          -0.001120f,
                                                          0.023083f,
                                                          1.47591054f - 5.5f,
                                                      },
                                                  },
                                                };

static const float partition_breakout_weights_8[RESOLUTION_CTX][Q_CTX]
                                               [FEATURES + 1] = {
                                                 {
                                                     {
                                                         -0.011807f,
                                                         -0.009873f,
                                                         -0.000931f,
                                                         0.034768f,
                                                         1.32254851f - 2.0f,
                                                     },
                                                     {
                                                         -0.003861f,
                                                         -0.002701f,
                                                         0.000100f,
                                                         0.013876f,
                                                         1.96755111f - 1.5f,
                                                     },
                                                     {
                                                         -0.013522f,
                                                         -0.008677f,
                                                         -0.000562f,
                                                         0.034468f,
                                                         1.53440356f - 1.5f,
                                                     },
                                                 },
                                                 {
                                                     {
                                                         -0.003221f,
                                                         -0.002125f,
                                                         0.000993f,
                                                         0.012768f,
                                                         0.03541421f - 2.0f,
                                                     },
                                                     {
                                                         -0.006069f,
                                                         -0.007335f,
                                                         0.000229f,
                                                         0.026104f,
                                                         0.17135315f - 1.5f,
                                                     },
                                                     {
                                                         -0.039894f,
                                                         -0.011419f,
                                                         0.000070f,
                                                         0.061817f,
                                                         0.6739977f - 1.5f,
                                                     },
                                                 },
                                               };

// ML-based partition search breakout.
static int ml_predict_breakout(const VP9_COMP *const cpi, BLOCK_SIZE bsize,
                               const MACROBLOCK *const x,
                               const RD_COST *const rd_cost) {
  DECLARE_ALIGNED(16, static const uint8_t, vp9_64_zeros[64]) = { 0 };
  const VP9_COMMON *const cm = &cpi->common;
  float features[FEATURES];
  const float *linear_weights = NULL;  // Linear model weights.
  float linear_score = 0.0f;
  const int qindex = cm->base_qindex;
  const int q_ctx = qindex >= 200 ? 0 : (qindex >= 150 ? 1 : 2);
  const int is_720p_or_larger = VPXMIN(cm->width, cm->height) >= 720;
  const int resolution_ctx = is_720p_or_larger ? 1 : 0;

  switch (bsize) {
    case BLOCK_64X64:
      linear_weights = partition_breakout_weights_64[resolution_ctx][q_ctx];
      break;
    case BLOCK_32X32:
      linear_weights = partition_breakout_weights_32[resolution_ctx][q_ctx];
      break;
    case BLOCK_16X16:
      linear_weights = partition_breakout_weights_16[resolution_ctx][q_ctx];
      break;
    case BLOCK_8X8:
      linear_weights = partition_breakout_weights_8[resolution_ctx][q_ctx];
      break;
    default: assert(0 && "Unexpected block size."); return 0;
  }
  if (!linear_weights) return 0;

  {  // Generate feature values.
    const int ac_q = vp9_ac_quant(qindex, 0, cm->bit_depth);
    const int num_pels_log2 = num_pels_log2_lookup[bsize];
    int feature_index = 0;
    unsigned int var, sse;
    float rate_f, dist_f;

    var = cpi->fn_ptr[bsize].vf(x->plane[0].src.buf, x->plane[0].src.stride,
                                vp9_64_zeros, 0, &sse);
    var = var >> num_pels_log2;

    vpx_clear_system_state();

    rate_f = (float)VPXMIN(rd_cost->rate, INT_MAX);
    dist_f = (float)(VPXMIN(rd_cost->dist, INT_MAX) >> num_pels_log2);
    rate_f =
        ((float)x->rdmult / 128.0f / 512.0f / (float)(1 << num_pels_log2)) *
        rate_f;

    features[feature_index++] = rate_f;
    features[feature_index++] = dist_f;
    features[feature_index++] = (float)var;
    features[feature_index++] = (float)ac_q;
    assert(feature_index == FEATURES);
  }

  {  // Calculate the output score.
    int i;
    linear_score = linear_weights[FEATURES];
    for (i = 0; i < FEATURES; ++i)
      linear_score += linear_weights[i] * features[i];
  }

  return linear_score >= cpi->sf.ml_partition_search_breakout_thresh[q_ctx];
}
#undef FEATURES
#undef Q_CTX
#undef RESOLUTION_CTX

#define FEATURES 17
#define LABELS 4
static const float rect_part_nn_weights_16_layer0[FEATURES * 32] = {
  1.262885f,  -0.533345f, -0.161280f, 0.106098f,  0.194799f,  0.003600f,
  0.394783f,  -0.053954f, 0.264474f,  -0.016651f, 0.376765f,  0.221471f,
  0.489799f,  0.054924f,  0.018292f,  0.037633f,  -0.053430f, 1.092426f,
  0.205791f,  -0.055661f, -0.227335f, 0.301274f,  -0.169917f, 0.100426f,
  0.254388f,  0.103465f,  0.189560f,  0.116479f,  1.647195f,  -0.667044f,
  0.067795f,  -0.044580f, 0.019428f,  0.072938f,  -0.797569f, -0.077539f,
  -0.225636f, 0.262883f,  -1.048009f, 0.210118f,  -0.416156f, -0.143741f,
  -0.296985f, 0.205918f,  -0.517383f, -0.118527f, -0.396606f, -0.113128f,
  -0.279468f, 0.096141f,  -0.342051f, -0.337036f, 0.143222f,  -0.860280f,
  0.137169f,  0.339767f,  -0.336076f, 0.071988f,  0.251557f,  -0.004068f,
  0.170734f,  0.237283f,  -0.332443f, 0.073643f,  0.375357f,  0.220407f,
  0.150708f,  -0.176979f, 0.265786f,  -0.105878f, -0.337465f, -0.000491f,
  0.234308f,  -0.098973f, 0.129038f,  -0.205936f, -0.034793f, -0.106981f,
  0.009974f,  0.037861f,  -0.282874f, -0.354414f, 0.023021f,  -0.266749f,
  -0.041762f, -0.721725f, 0.182262f,  -0.273945f, 0.123722f,  -0.036749f,
  -0.788645f, -0.081560f, -0.472226f, 0.004654f,  -0.756766f, -0.132186f,
  1.085412f,  -0.221324f, -0.072577f, -0.172834f, -0.104831f, -1.391641f,
  -0.345893f, 0.194442f,  -0.306583f, -0.041813f, -0.267635f, -0.218568f,
  -0.178452f, 0.044421f,  -0.128042f, -0.094797f, -0.253724f, 0.273931f,
  0.144843f,  -0.401416f, -0.014354f, -0.348929f, 0.123550f,  0.494504f,
  -0.007050f, -0.143830f, 0.111292f,  0.211057f,  -1.579988f, 0.117744f,
  -1.732487f, 0.009320f,  -1.162696f, 0.176687f,  -0.705609f, 0.524827f,
  0.089822f,  0.082976f,  -0.023681f, 0.006120f,  -0.907175f, -0.026273f,
  0.019027f,  0.027170f,  -0.462563f, -0.535335f, 0.202231f,  0.709803f,
  -0.112251f, -1.213869f, 0.225714f,  0.323785f,  -0.518254f, -0.014235f,
  -0.070790f, -0.369589f, 0.373399f,  0.002738f,  0.175113f,  0.084529f,
  -0.101586f, -0.018978f, 0.773392f,  -0.673230f, -0.549279f, 0.790196f,
  0.658609f,  -0.826831f, -0.514211f, 0.575341f,  -0.711311f, 0.276289f,
  -0.435715f, 0.392986f,  -0.079298f, -0.318719f, 0.188429f,  -0.114366f,
  0.172527f,  -0.261721f, -0.216761f, 0.163822f,  -0.189374f, -0.391901f,
  0.142013f,  -0.135046f, 0.144419f,  0.053887f,  0.074673f,  -0.290791f,
  -0.039560f, -0.103830f, -0.330263f, -0.042091f, 0.050646f,  -0.057466f,
  -0.069064f, -0.412864f, 0.071097f,  0.126693f,  0.175397f,  -0.168485f,
  0.018129f,  -0.419188f, -0.272024f, -0.436859f, -0.425711f, -0.024382f,
  0.248042f,  -0.169090f, -0.346878f, -0.070926f, 0.292278f,  -0.197610f,
  -0.218286f, 0.290846f,  0.297843f,  0.247394f,  -0.160736f, 0.110314f,
  0.276000f,  -0.301676f, -0.232816f, -0.127576f, -0.174457f, -0.124503f,
  0.264880f,  -0.332379f, 0.012659f,  -0.197333f, 0.604700f,  0.801582f,
  0.758702f,  0.691880f,  0.440917f,  0.773548f,  0.064242f,  1.147508f,
  -0.127543f, -0.189628f, -0.122994f, -0.226776f, -0.053531f, -0.187548f,
  0.226554f,  -0.273451f, 0.011751f,  0.009133f,  0.185091f,  0.003031f,
  0.000525f,  0.221829f,  0.331550f,  -0.202558f, -0.286550f, 0.100683f,
  0.268818f,  0.179971f,  -0.050016f, 0.579665f,  0.015911f,  0.033068f,
  0.077768f,  -0.017757f, -1.411251f, 0.051519f,  -1.745767f, 0.011258f,
  -1.947372f, 0.111396f,  -1.112755f, -0.008989f, -0.006211f, -0.002098f,
  -0.015236f, -0.095697f, -0.095820f, 0.044622f,  -0.112096f, 0.060000f,
  0.138957f,  -0.462708f, 0.590790f,  -0.021405f, -0.283744f, -1.141749f,
  0.213121f,  -0.332311f, -0.314090f, -0.789311f, 0.157605f,  -0.438019f,
  0.642189f,  -0.340764f, -0.996025f, 0.109871f,  0.106128f,  -0.010505f,
  -0.117233f, -0.223194f, 0.344105f,  -0.308754f, 0.386020f,  -0.305270f,
  -0.538281f, -0.270720f, -0.101688f, 0.207580f,  0.237153f,  -0.055730f,
  0.842779f,  0.393543f,  0.007886f,  -0.318167f, 0.603768f,  0.388241f,
  0.421536f,  0.632080f,  0.423965f,  0.371472f,  0.456827f,  0.488134f,
  0.358997f,  0.032621f,  -0.017104f, 0.032198f,  0.113266f,  -0.312277f,
  0.178189f,  0.234180f,  0.134271f,  -0.414889f, 0.774141f,  -0.225043f,
  0.614052f,  -0.279921f, 1.329141f,  -0.140827f, 0.797267f,  -0.171361f,
  0.066205f,  0.339976f,  0.015223f,  0.193725f,  -0.245067f, -0.035578f,
  -0.084043f, 0.086756f,  0.029478f,  -0.845370f, 0.388613f,  -1.215236f,
  0.304573f,  -0.439884f, -0.293969f, -0.107988f, -0.267837f, -0.695339f,
  -0.702099f, 0.359047f,  0.511730f,  1.429516f,  0.216959f,  -0.313828f,
  0.068062f,  -0.124917f, -0.648327f, -0.308411f, -0.378467f, -0.429288f,
  -0.032415f, -0.357005f, 0.170068f,  0.161167f,  -0.250280f, -0.320468f,
  -0.408987f, -0.201496f, -0.155996f, 0.021067f,  0.141083f,  -0.202733f,
  -0.130953f, -0.278148f, -0.042051f, 0.070576f,  0.009982f,  -0.044326f,
  -0.346851f, -0.255397f, -0.346456f, 0.281781f,  0.001618f,  0.120648f,
  0.297140f,  0.198343f,  0.186104f,  0.183548f,  -0.344482f, 0.182258f,
  0.291003f,  -0.330228f, -0.048174f, 0.133694f,  0.264582f,  0.229671f,
  -0.167251f, -0.316040f, 0.191829f,  0.153417f,  -0.345158f, -0.212790f,
  -0.878872f, -0.313099f, -0.028368f, 0.065869f,  -0.695388f, 1.102812f,
  -0.605539f, 0.400680f,  -0.350120f, -0.432965f, 0.034553f,  -0.693476f,
  -0.045708f, 0.492409f,  -0.043825f, -0.430522f, 0.071159f,  -0.317376f,
  -1.164842f, 0.112394f,  0.034137f,  -0.611882f, 0.251020f,  -0.245113f,
  0.286093f,  -0.187883f, 0.340263f,  -0.211592f, -0.065706f, -0.332148f,
  0.104026f,  -0.003206f, 0.036397f,  0.206499f,  0.161962f,  0.037663f,
  -0.313039f, -0.199837f, 0.117952f,  -0.182145f, -0.343724f, 0.017625f,
  0.033427f,  -0.288075f, -0.101873f, -0.083378f, 0.147870f,  0.049598f,
  -0.241824f, 0.070494f,  0.140942f,  -0.013795f, 0.020023f,  -0.192213f,
  -0.320505f, -0.193072f, 0.147260f,  0.311352f,  0.053486f,  0.183716f,
  0.142535f,  0.294333f,  -0.054853f, 0.293314f,  -0.025398f, 0.190815f,
  -0.137574f, -0.191864f, -0.190950f, -0.205988f, -0.199046f, -0.017582f,
  -0.149347f, 0.131040f,  0.006854f,  -0.350732f, 0.113301f,  -0.194371f,
  -0.296885f, -0.249199f, -0.193946f, 0.116150f,  -0.310411f, -0.325851f,
  -0.053275f, -0.063419f, 0.204170f,  -0.091940f, -0.146229f, 0.298173f,
  0.053349f,  -0.368540f, 0.235629f,  -0.317825f, -0.107304f, -0.114618f,
  0.058709f,  -0.272070f, 0.076224f,  0.110668f,  -0.193282f, -0.135440f,
  -0.267950f, -0.102285f, 0.102699f,  -0.159082f, 0.262721f,  -0.263227f,
  0.094509f,  -0.113405f, 0.069888f,  -0.169665f, 0.070800f,  0.035432f,
  0.054243f,  0.264229f,  0.117416f,  0.091568f,  -0.022069f, -0.069214f,
  0.124543f,  0.070413f,  -0.039343f, 0.082823f,  -0.838348f, 0.153727f,
  -0.000947f, 0.270348f,  -1.404952f, -0.159680f, -0.234320f, 0.061023f,
  0.271660f,  -0.541834f, 0.570828f,  -0.277254f,
};

static const float rect_part_nn_bias_16_layer0[32] = {
  0.045740f,  0.292685f,  -0.754007f, -0.150412f, -0.006171f, 0.005915f,
  0.000167f,  0.322797f,  -0.381793f, 0.349786f,  0.003878f,  -0.307203f,
  0.000000f,  0.029122f,  0.000000f,  0.625494f,  0.302105f,  -0.362807f,
  -0.034002f, -0.573278f, 0.240021f,  0.083965f,  0.000000f,  -0.018979f,
  -0.147739f, -0.036990f, 0.000000f,  0.000000f,  -0.026790f, -0.000036f,
  -0.073448f, 0.398328f,
};

static const float rect_part_nn_weights_16_layer1[32 * LABELS] = {
  0.095090f,  0.831754f,  0.484433f,  0.472945f,  0.086165f,  -0.442388f,
  0.176263f,  -0.760247f, 0.419932f,  -0.131377f, 0.075814f,  0.089844f,
  -0.294718f, 0.299808f,  -0.318435f, -0.623205f, -0.346703f, 0.494356f,
  0.949221f,  0.524653f,  0.044095f,  0.428540f,  0.402571f,  -0.216920f,
  0.423915f,  1.023334f,  -0.366449f, 0.395057f,  0.057576f,  0.094019f,
  0.247685f,  -0.007200f, -0.420023f, -0.728965f, -0.063040f, -0.071321f,
  0.209298f,  0.486625f,  -0.244375f, 0.263219f,  -0.250463f, -0.260301f,
  0.068579f,  0.177644f,  -0.155311f, -0.027606f, -0.101614f, 0.553046f,
  -0.462729f, -0.237568f, -0.589316f, 0.045182f,  0.551759f,  -0.196872f,
  0.183040f,  0.054341f,  0.252784f,  -0.536486f, -0.024425f, 0.154942f,
  -0.086636f, 0.360416f,  0.214773f,  -0.170876f, -0.363522f, -0.464099f,
  0.145494f,  -0.099329f, 0.343718f,  0.286427f,  0.085540f,  -0.105182f,
  0.155543f,  0.290939f,  -0.067069f, 0.228399f,  0.178247f,  0.113031f,
  -0.067336f, 0.441062f,  0.132364f,  -0.263403f, -0.263925f, -0.083613f,
  -0.268577f, -0.204442f, 0.052526f,  0.334787f,  -0.064285f, -0.197875f,
  0.296405f,  0.396440f,  0.033231f,  0.229087f,  0.118289f,  0.490894f,
  -0.527582f, -0.897206f, -0.325708f, -0.433018f, -0.053989f, 0.223814f,
  -0.352319f, 0.772440f,  -0.108648f, -0.082859f, -0.342718f, 0.033022f,
  -0.309199f, -0.560337f, 0.208476f,  0.520309f,  -0.241035f, -0.560391f,
  -1.268968f, -0.267567f, 0.129461f,  -0.385547f, 0.080142f,  0.065785f,
  -0.159324f, -0.580704f, -0.315150f, -0.224900f, -0.110807f, -0.230163f,
  0.307266f,  0.153446f,
};

static const float rect_part_nn_bias_16_layer1[LABELS] = {
  -0.455437f,
  0.255310f,
  0.452974f,
  -0.278733f,
};

static const NN_CONFIG rect_part_nnconfig_16 = {
  FEATURES,  // num_inputs
  LABELS,    // num_outputs
  1,         // num_hidden_layers
  {
      32,
  },  // num_hidden_nodes
  {
      rect_part_nn_weights_16_layer0,
      rect_part_nn_weights_16_layer1,
  },
  {
      rect_part_nn_bias_16_layer0,
      rect_part_nn_bias_16_layer1,
  },
};

static const float rect_part_nn_weights_32_layer0[FEATURES * 32] = {
  0.735110f,  -0.238477f, 0.101978f,  0.311671f,  -0.123833f, 1.596506f,
  -0.341982f, -0.480170f, -0.247587f, 0.613159f,  -0.279899f, -0.740856f,
  0.499051f,  0.039041f,  0.056763f,  0.258874f,  0.470812f,  -0.121635f,
  -0.318852f, -0.098677f, -0.214714f, -0.159974f, -0.305400f, -0.344477f,
  -0.260653f, -0.007737f, -0.053016f, -0.158079f, 0.151911f,  -0.057685f,
  -0.230948f, -0.165940f, -0.127591f, -0.192084f, 1.890390f,  -0.315123f,
  -0.714531f, -0.015355f, 0.186437f,  0.305504f,  0.035343f,  -0.556783f,
  0.239364f,  -0.297789f, 0.202735f,  -0.707576f, 0.710250f,  0.223346f,
  -0.291511f, 0.235778f,  0.455338f,  -0.059402f, 0.084530f,  -0.115117f,
  -0.103696f, -0.192821f, 0.114579f,  -0.223487f, 0.306864f,  0.021887f,
  -0.028040f, 0.087866f,  0.038870f,  -0.081742f, -0.056052f, -0.130837f,
  0.201058f,  0.293391f,  1.880344f,  0.339162f,  0.040928f,  -0.503942f,
  0.476333f,  0.259272f,  0.629416f,  0.869369f,  0.622841f,  1.012843f,
  0.715795f,  1.958844f,  -1.697462f, 0.071334f,  0.074189f,  0.014585f,
  -0.002536f, 0.021900f,  0.151883f,  0.169501f,  -0.333018f, -0.247512f,
  -0.418575f, -0.473960f, -0.004501f, -0.280939f, -0.162188f, -0.355632f,
  0.136654f,  -0.100967f, -0.350435f, -0.135386f, 0.037237f,  0.136982f,
  -0.084157f, -0.073248f, 0.021792f,  0.077429f,  -0.083042f, -3.169569f,
  0.016261f,  -3.351328f, 0.021120f,  -3.572247f, 0.023870f,  -4.312754f,
  0.040973f,  -0.038328f, -0.015052f, 0.017702f,  0.101427f,  0.115458f,
  -0.304792f, 0.021826f,  -0.157998f, 0.341022f,  -0.013465f, 0.105076f,
  -0.261465f, 0.318730f,  0.065701f,  0.314879f,  -0.064785f, 0.282824f,
  0.100542f,  0.057260f,  -0.003756f, -0.026214f, -0.264641f, 0.275545f,
  -0.049201f, -0.283015f, -0.057363f, 0.183570f,  0.243161f,  -0.255764f,
  0.099747f,  -0.156157f, -0.262494f, 0.231521f,  -0.262617f, -0.186096f,
  0.171720f,  0.018983f,  -0.145545f, 0.197662f,  -0.001502f, -0.267526f,
  0.001960f,  0.003260f,  0.045237f,  -0.377174f, -0.042499f, -0.015278f,
  -0.196779f, -0.262797f, -0.318427f, -0.126092f, -0.339723f, 0.205288f,
  -0.544284f, -0.507896f, -0.316622f, -0.090312f, -0.250917f, -0.337263f,
  -0.220199f, -0.296591f, -0.116816f, 0.052381f,  0.145681f,  0.016521f,
  -0.093549f, -0.097822f, 0.023140f,  -0.010346f, 0.036181f,  0.145826f,
  -0.139123f, -0.462638f, -0.007315f, 0.156533f,  -0.102787f, 0.143586f,
  -0.092094f, -0.144220f, -0.168994f, -0.045833f, 0.021628f,  -0.421794f,
  -0.055857f, 0.217931f,  -0.061937f, -0.028768f, -0.078250f, -0.426939f,
  -0.223118f, -0.230080f, -0.194988f, -0.197673f, -0.020918f, 0.139945f,
  0.186951f,  -0.071317f, -0.084007f, -0.138597f, 0.101950f,  0.093870f,
  0.153226f,  0.017799f,  -0.088539f, -0.037796f, 0.340412f,  0.183305f,
  0.391880f,  -1.127417f, 0.132762f,  -0.228565f, 0.399035f,  0.017483f,
  -0.041619f, 0.017849f,  0.092340f,  0.054204f,  0.681185f,  0.421034f,
  0.112520f,  -0.040618f, -0.040148f, -0.360647f, 0.053555f,  0.192854f,
  0.076968f,  -0.179224f, -0.081617f, -0.287661f, -0.191072f, -0.310227f,
  -0.332226f, -0.039786f, -0.247795f, -0.232201f, -0.333533f, -0.077995f,
  -0.471732f, 0.051829f,  0.090488f,  0.142465f,  -0.120490f, -0.286151f,
  -0.049117f, -0.251082f, 0.211884f,  -0.223366f, 0.063565f,  0.229938f,
  -0.059348f, -0.029573f, -0.064303f, -0.156148f, 0.086958f,  -0.297613f,
  -0.125107f, 0.062718f,  0.339137f,  -0.218896f, -0.057290f, -0.236670f,
  -0.143783f, -0.119429f, 0.242320f,  -0.323464f, -0.178377f, 0.238275f,
  -0.025042f, 0.074798f,  0.111329f,  -0.299773f, -0.151748f, -0.261607f,
  0.215626f,  0.202243f,  -0.121896f, -0.024283f, -0.293854f, -0.018232f,
  -0.012629f, -0.199297f, -0.060595f, 0.432339f,  -0.158735f, -0.028380f,
  0.326639f,  0.222546f,  -0.218135f, -0.495955f, -0.015055f, -0.104206f,
  -0.268823f, 0.116765f,  0.041769f,  -0.187095f, 0.225090f,  0.198195f,
  0.001502f,  -0.219212f, -0.244779f, -0.017690f, -0.033197f, -0.339813f,
  -0.325453f, 0.002499f,  -0.066113f, 0.043235f,  0.324275f,  -0.630642f,
  -1.440551f, 0.174527f,  0.124619f,  -1.187345f, 1.372693f,  -0.278393f,
  -0.058673f, -0.286338f, 1.708757f,  -0.325094f, -0.543172f, -0.229411f,
  0.169927f,  0.175064f,  0.198321f,  0.117351f,  0.220882f,  0.138078f,
  -0.158000f, -0.286708f, 0.096046f,  -0.321788f, 0.206949f,  -0.014473f,
  -0.321234f, 0.100033f,  -0.108266f, 0.166824f,  0.032904f,  -0.065760f,
  -0.303896f, 0.180342f,  -0.301145f, -0.352554f, 0.149089f,  0.013277f,
  0.256019f,  -0.109770f, 1.832588f,  -0.132568f, 1.527658f,  -0.164252f,
  -0.857880f, -0.242694f, -0.553797f, 0.334023f,  -0.332759f, -0.166203f,
  -0.223175f, 0.007953f,  -0.175865f, -0.134590f, -0.023858f, -0.011983f,
  0.054403f,  -0.147054f, -0.176901f, -0.166893f, -0.292662f, -0.010569f,
  -0.041744f, -0.060398f, -0.237584f, 0.154246f,  -0.083270f, -0.314016f,
  -0.374736f, 0.100063f,  0.048401f,  -0.061952f, -0.178816f, 0.157243f,
  0.221991f,  -0.065035f, 0.098517f,  -0.190704f, -0.210613f, -0.274884f,
  -0.341442f, -0.205281f, 0.073644f,  0.130667f,  0.149194f,  -0.018172f,
  1.796154f,  -1.017806f, -0.169655f, 0.104239f,  0.344313f,  0.643042f,
  0.730177f,  0.270776f,  0.581631f,  -1.090649f, 0.707472f,  1.411035f,
  0.268739f,  0.178860f,  -0.062251f, -0.118611f, -0.215759f, 0.023485f,
  -0.105320f, 0.036396f,  -0.059604f, 0.090024f,  0.095224f,  -0.053497f,
  -0.084040f, 0.055836f,  0.111678f,  0.014886f,  -0.178380f, 0.079662f,
  -0.123580f, 0.057379f,  -0.409844f, -0.305386f, -0.987808f, -0.291094f,
  0.063966f,  0.263709f,  -0.337221f, 0.720093f,  0.105030f,  0.848950f,
  0.071835f,  0.228972f,  0.057705f,  -2.154561f, -0.201303f, -0.058856f,
  -0.020081f, 0.029375f,  0.234837f,  -0.001063f, 0.042527f,  0.014567f,
  -0.299420f, -0.289117f, 0.275219f,  0.263596f,  -0.186026f, -0.111364f,
  -0.118393f, -0.318778f, 0.010710f,  -0.286836f, -0.070330f, -0.049497f,
  0.093162f,  -0.298085f, 0.204761f,  -0.206633f, -0.009057f, -0.235372f,
  0.185300f,  -0.271814f, 0.281732f,  0.268149f,  -0.018967f, 0.162748f,
  -0.086694f, -0.063839f, -0.097473f, -0.280120f, 0.324688f,  0.157911f,
  -0.064794f, -0.266017f, -0.305608f, -0.196854f, -0.185767f, 0.199455f,
  0.102264f,  0.070866f,  0.172045f,  0.266433f,  -0.176167f, 0.251657f,
  -0.239220f, 0.229667f,  0.156115f,  -0.221345f, 0.270720f,  0.109367f,
  0.230352f,  -0.384561f, -0.026329f, 0.005928f,  -0.087685f, -0.097995f,
  -0.153864f, 0.117211f,  -0.226492f, -0.379832f, -0.201714f, 0.049707f,
  -0.292120f, 0.114074f,  -0.085307f, -0.485356f, -0.347405f, 0.089361f,
  -0.419273f, -0.320764f, -0.107254f, -0.274615f, -0.292991f, 0.095602f,
  -0.078789f, 0.138927f,  0.270813f,  0.205814f,  0.065003f,  0.169171f,
  0.056142f,  -0.005792f, 0.059483f,  0.060149f,
};

static const float rect_part_nn_bias_32_layer0[32] = {
  -1.749808f, 0.000000f,  0.239736f,  -0.000424f, 0.431792f,  -0.150833f,
  2.866760f,  0.000000f,  0.000000f,  -0.281434f, 0.000000f,  -0.150086f,
  0.000000f,  -0.008346f, -0.204104f, -0.006581f, 0.000000f,  -0.197006f,
  0.000000f,  -0.735287f, -0.028345f, -1.180116f, -0.106524f, 0.000000f,
  0.075879f,  -0.150966f, -2.438914f, 0.000000f,  -0.011775f, -0.024204f,
  -0.138235f, -0.123763f,
};

static const float rect_part_nn_weights_32_layer1[32 * LABELS] = {
  0.622235f,  0.264894f,  -0.424216f, 0.103989f,  1.401192f,  -0.063838f,
  -5.216846f, 0.329234f,  -0.293113f, 0.457519f,  -0.271899f, 0.043771f,
  -0.203823f, 0.573535f,  -0.192703f, 0.054939f,  0.163019f,  0.124803f,
  0.160664f,  0.385406f,  -0.091403f, 0.320204f,  0.101181f,  -0.157792f,
  -0.095555f, -0.255011f, 1.326614f,  -0.138076f, -0.082434f, -0.342442f,
  0.184067f,  -0.076395f, 0.050263f,  0.251065f,  0.291743f,  0.197838f,
  -0.950922f, 0.280202f,  2.904905f,  -0.219434f, 0.284386f,  0.375005f,
  0.193817f,  -0.298663f, -0.255364f, -0.297545f, 0.030518f,  -0.023892f,
  -0.396120f, -0.253027f, 0.237235f,  -0.550249f, -0.076817f, -0.201374f,
  0.292708f,  0.341936f,  -0.532215f, 0.180634f,  -0.943291f, -0.217179f,
  0.251611f,  -0.306310f, 0.229054f,  -0.350337f, -0.192707f, 0.146781f,
  0.409007f,  0.279088f,  -0.307357f, 0.199059f,  2.780962f,  0.163723f,
  -0.226445f, 0.242830f,  0.220356f,  -0.057621f, 0.196677f,  -0.179975f,
  -0.314636f, 0.218271f,  -0.278653f, -0.226286f, 0.034275f,  -0.320149f,
  0.154779f,  0.074937f,  -0.015650f, -0.281735f, -0.495227f, -0.075036f,
  -0.871024f, -0.350643f, 0.343468f,  0.095665f,  0.447121f,  -0.059040f,
  0.244757f,  0.223122f,  0.272544f,  0.129678f,  -1.700183f, 0.254869f,
  2.528983f,  0.217362f,  0.327765f,  -0.129369f, -0.003560f, -0.532537f,
  0.080216f,  -0.739488f, -0.299813f, 0.185421f,  0.265994f,  0.152268f,
  -0.401829f, -0.901380f, 0.347747f,  -0.524845f, -0.201163f, 0.063585f,
  -0.517479f, -0.077816f, -0.735739f, -0.161411f, -0.113607f, -0.306188f,
  0.190817f,  -0.362567f,
};

static const float rect_part_nn_bias_32_layer1[LABELS] = {
  -0.833530f,
  0.860502f,
  0.708645f,
  -1.083700f,
};

static const NN_CONFIG rect_part_nnconfig_32 = {
  FEATURES,  // num_inputs
  LABELS,    // num_outputs
  1,         // num_hidden_layers
  {
      32,
  },  // num_hidden_nodes
  {
      rect_part_nn_weights_32_layer0,
      rect_part_nn_weights_32_layer1,
  },
  {
      rect_part_nn_bias_32_layer0,
      rect_part_nn_bias_32_layer1,
  },
};

static const float rect_part_nn_weights_64_layer0[FEATURES * 32] = {
  0.029424f,  -0.295893f, -0.313259f, -0.090484f, -0.104946f, 0.121361f,
  0.137971f,  -0.137984f, -0.328158f, -0.137280f, -0.276995f, -0.153118f,
  0.187893f,  0.105787f,  -0.236591f, -0.114325f, -0.000708f, 1.936191f,
  0.048491f,  -0.026048f, -0.206916f, 0.830237f,  -0.152354f, 0.074191f,
  -0.153813f, 0.148942f,  -0.103457f, 0.028252f,  1.758264f,  -2.123016f,
  0.120182f,  0.049954f,  0.110450f,  -0.199360f, 0.642198f,  0.040225f,
  -0.140886f, 0.091833f,  -0.122788f, 1.172115f,  -0.833333f, -0.505218f,
  0.736050f,  -0.109958f, -0.839030f, -0.399916f, 1.029718f,  0.408977f,
  -0.836882f, 0.389683f,  -1.134413f, -1.529672f, -0.146351f, 0.089298f,
  0.083772f,  -0.697869f, 1.683311f,  -0.882446f, 0.494428f,  -0.122128f,
  0.659819f,  -0.057178f, -0.915390f, -0.192412f, 0.046613f,  0.010697f,
  0.040782f,  0.110807f,  -0.225332f, -0.327730f, -0.114825f, 0.063511f,
  0.050503f,  0.023602f,  0.006524f,  -0.274547f, -0.607145f, -0.143812f,
  -0.327689f, -0.333072f, -0.017138f, -0.183992f, -0.200622f, -0.262463f,
  -0.132799f, -0.018155f, -0.534214f, -0.385994f, 0.116278f,  -0.752879f,
  -0.090734f, -0.249152f, 0.071716f,  0.029603f,  -0.382456f, -0.122894f,
  1.349552f,  -0.885192f, 0.257903f,  -0.265945f, -0.045579f, 0.112247f,
  -0.122810f, -0.258285f, -0.145427f, -0.127442f, 0.072778f,  0.072549f,
  0.182149f,  0.239403f,  0.167205f,  -0.291616f, -0.281237f, 0.335735f,
  0.208511f,  -0.239628f, -0.022236f, -0.177370f, 0.207808f,  0.023535f,
  0.137455f,  0.016406f,  -0.138685f, 0.188732f,  0.205513f,  0.209787f,
  0.060592f,  0.239954f,  -0.128341f, -0.291585f, 0.022141f,  -0.311201f,
  -0.010199f, -0.314224f, -0.351915f, -0.079775f, -0.260028f, -0.015953f,
  0.007404f,  0.051589f,  0.019771f,  -2.337926f, 0.024596f,  -2.512399f,
  -0.023138f, -2.421380f, 0.016515f,  -3.269775f, 0.026844f,  -0.053660f,
  -0.013213f, -0.029248f, 0.114357f,  0.259100f,  -0.141749f, -0.106802f,
  -0.117323f, -0.294698f, -0.316012f, -0.328013f, 0.016459f,  0.136175f,
  0.223327f,  0.322312f,  -0.297297f, 0.118286f,  -0.317197f, -0.116692f,
  0.262236f,  -0.032443f, -0.392128f, -0.199989f, -0.383621f, 0.008347f,
  -0.079302f, -0.005529f, 0.049261f,  0.145948f,  -0.263592f, -0.317109f,
  0.260015f,  -0.499341f, -0.171764f, -0.017815f, 0.149186f,  0.178294f,
  -0.492198f, 0.016956f,  0.008067f,  -0.057734f, -0.189979f, -0.131489f,
  -0.163303f, 0.121378f,  -0.172272f, 0.125891f,  0.120654f,  0.071314f,
  0.117423f,  -0.242167f, 0.047170f,  0.234302f,  -0.355370f, -0.336112f,
  -0.255471f, -0.267792f, -0.135367f, -0.284411f, 0.254592f,  0.098749f,
  0.224989f,  0.258450f,  -0.306878f, 0.153551f,  -0.175806f, -0.244459f,
  -0.274922f, 0.254346f,  0.110309f,  0.036054f,  0.095133f,  -0.589646f,
  0.080543f,  0.154155f,  0.133797f,  -0.401518f, 0.798127f,  0.066742f,
  1.449216f,  0.282498f,  1.210638f,  -0.280643f, 0.572386f,  -0.308133f,
  -0.053143f, 0.008437f,  0.269565f,  0.347616f,  0.087180f,  -0.771104f,
  0.200800f,  0.157578f,  0.474128f,  -0.971488f, 0.193451f,  0.340339f,
  -0.123425f, 0.560754f,  -0.139621f, -0.281721f, -0.100162f, 0.250926f,
  0.281100f,  0.197680f,  0.138629f,  1.045823f,  0.339047f,  0.036698f,
  -0.159210f, 0.727869f,  -1.371850f, 0.116241f,  -2.180194f, 0.214055f,
  -0.213691f, 0.447957f,  -1.129966f, 0.543598f,  0.147599f,  0.060034f,
  -0.049415f, -0.095858f, 0.290599f,  0.059512f,  0.198343f,  -0.211903f,
  0.158736f,  -0.090220f, -0.221992f, 0.198320f,  0.028632f,  -0.408238f,
  -0.368266f, -0.218740f, -0.379023f, -0.173573f, -0.035179f, 0.240176f,
  0.237714f,  -0.417132f, -0.184989f, 0.046818f,  -0.016965f, -0.524012f,
  -0.094848f, -0.225678f, 0.021766f,  -0.028366f, 0.072343f,  -0.039980f,
  0.023334f,  -0.392397f, 0.164450f,  -0.201650f, -0.519754f, -0.023352f,
  -4.559466f, -0.115996f, 0.135844f,  0.152599f,  -0.111570f, 1.870310f,
  0.003522f,  1.893098f,  -0.134055f, 1.850787f,  0.085160f,  -2.203354f,
  0.380799f,  -0.074047f, 0.023760f,  0.077310f,  0.273381f,  -1.163135f,
  -0.024976f, 0.093252f,  0.011445f,  -0.129009f, -2.200677f, -0.013703f,
  -1.964109f, -0.027246f, -2.135679f, 0.049465f,  -3.879032f, 0.195114f,
  -0.018085f, 0.016755f,  0.036330f,  0.169138f,  0.003548f,  -0.028565f,
  -0.178196f, -0.020577f, -0.104330f, -0.270961f, -0.282822f, -0.228735f,
  -0.292561f, 0.271648f,  0.129171f,  0.376168f,  -0.265005f, -0.093002f,
  -0.185514f, 0.025598f,  0.055265f,  -0.212784f, -0.249005f, 0.051507f,
  -0.267868f, 0.162227f,  -0.237365f, 0.267479f,  -0.051543f, -0.288800f,
  -0.246119f, 0.216296f,  0.226888f,  -0.123005f, 0.068040f,  -0.096630f,
  -0.100500f, 0.161640f,  -0.349187f, -0.061229f, 0.042915f,  0.024949f,
  -0.083086f, -0.407249f, -0.428306f, -0.381137f, -0.508822f, 0.354796f,
  -0.612346f, -0.230076f, -0.734103f, -0.550571f, -0.318788f, -0.300091f,
  -0.336045f, -0.494406f, -0.206900f, 0.079942f,  0.149065f,  -0.533360f,
  0.940431f,  -0.078860f, 1.418633f,  -0.117527f, 1.349170f,  0.242658f,
  0.559328f,  0.258770f,  -0.014508f, -0.204775f, -0.292631f, 0.498345f,
  -0.274918f, 0.051670f,  0.157748f,  -0.179721f, -0.183330f, -0.393550f,
  -0.208848f, 0.060742f,  -0.159654f, 0.047757f,  -0.400256f, -0.084606f,
  -0.080619f, -0.359664f, -0.078305f, -0.455653f, 0.227624f,  -0.385606f,
  -0.060326f, -0.209831f, -0.077008f, 0.148862f,  0.209908f,  0.047655f,
  -0.342292f, -0.088375f, -0.115465f, 0.082700f,  0.036465f,  -0.001792f,
  -0.285730f, 0.114632f,  0.239254f,  -0.348543f, 0.044916f,  -0.299003f,
  -0.244756f, -0.180802f, 0.314253f,  -0.127788f, -0.221512f, 0.034787f,
  -0.208388f, 0.349156f,  0.265975f,  -0.068335f, 0.261372f,  0.146705f,
  -0.098729f, 0.293699f,  -0.111342f, 0.207402f,  -0.038772f, 0.124135f,
  -0.237450f, -0.191511f, -0.052240f, -0.237151f, 0.005013f,  0.139441f,
  -0.153634f, -0.021596f, -0.036220f, -0.077873f, -0.085995f, -0.254555f,
  -0.204382f, -0.082362f, 0.941796f,  0.253800f,  -0.957468f, 0.095795f,
  0.122046f,  -0.310364f, 0.087301f,  0.012704f,  0.193265f,  -0.058303f,
  0.250452f,  0.835269f,  0.507383f,  0.109957f,  -0.145028f, -0.114419f,
  -0.225618f, 0.132387f,  -0.063335f, -0.325776f, -0.346173f, -0.006653f,
  -0.133534f, -0.085549f, -0.050177f, 0.173103f,  0.025421f,  0.105512f,
  0.258036f,  0.153116f,  0.290202f,  -0.333699f, -0.072405f, -0.124069f,
  -0.241933f, -0.313318f, 0.013623f,  -0.237440f, -0.232228f, -0.170850f,
  -0.039212f, 0.162468f,  -0.330162f, -0.218462f, -0.287064f, -0.181673f,
  -0.161059f, 0.024664f,  -0.108642f, -0.231707f, 0.217994f,  -1.128878f,
  0.093010f,  0.101513f,  0.055895f,  -0.354538f, 0.844174f,  0.254335f,
  1.920298f,  -0.230777f, 0.798144f,  0.206425f,  0.580655f,  -0.177645f,
  -0.412061f, 0.112629f,  -0.476438f, 0.209436f,
};

static const float rect_part_nn_bias_64_layer0[32] = {
  0.000000f,  0.345406f,  -0.499542f, -1.718246f, -0.147443f, -0.408843f,
  -0.008997f, -0.107946f, 2.117510f,  0.000000f,  -0.141830f, -0.049079f,
  0.000000f,  -1.331136f, -1.417843f, -0.485054f, -0.100856f, -0.230750f,
  -2.574372f, 2.310627f,  -0.030363f, 0.000000f,  -0.310119f, -1.314316f,
  -0.108766f, -0.107918f, 0.000000f,  0.000000f,  0.093643f,  0.000000f,
  0.000000f,  -0.902343f,
};

static const float rect_part_nn_weights_64_layer1[32 * LABELS] = {
  0.404567f,  1.168492f,  0.051714f,  0.827941f,  0.135334f,  0.456922f,
  -0.370524f, 0.062865f,  -3.076300f, -0.290613f, 0.280029f,  -0.101778f,
  0.250216f,  0.347721f,  0.466400f,  0.030845f,  0.114570f,  0.089456f,
  1.519938f,  -3.493788f, 0.264212f,  -0.109125f, 0.306644f,  0.368206f,
  -0.052168f, -0.229630f, -0.339932f, -0.080472f, 0.319845f,  0.143818f,
  -0.172595f, 0.372777f,  -0.082072f, -0.505781f, -0.288321f, -0.473028f,
  -0.027567f, -0.034329f, -0.291965f, -0.063262f, 1.721741f,  0.118914f,
  0.183681f,  0.041611f,  0.266371f,  0.005896f,  -0.484705f, 0.665535f,
  -0.240945f, -0.017963f, -1.409440f, 2.031976f,  0.240327f,  -0.116604f,
  0.273245f,  -0.170570f, -0.085491f, -0.340315f, -0.209651f, -0.217460f,
  -0.249373f, 0.009193f,  0.009467f,  -0.272909f, 0.308472f,  -0.551173f,
  0.168374f,  -0.583229f, 0.140082f,  -0.585715f, -0.010929f, 0.159779f,
  1.438104f,  0.293111f,  -0.053339f, -0.101828f, -0.280573f, -0.211265f,
  -0.323605f, -0.540908f, 0.101366f,  -0.005288f, -1.517046f, 2.078767f,
  0.215597f,  0.144012f,  0.315888f,  -0.251324f, 0.150482f,  -0.137871f,
  0.235116f,  -0.194202f, -0.153475f, -0.312384f, -0.375510f, 0.336488f,
  -0.379837f, -1.004979f, -0.312587f, -0.406174f, 0.154290f,  -0.539766f,
  -0.230074f, 0.303564f,  0.719439f,  -0.235108f, -0.204978f, 0.399229f,
  0.290222f,  -0.278713f, -0.667069f, -0.420550f, 0.164893f,  -0.459689f,
  -1.035368f, 0.818909f,  0.275137f,  -0.291006f, -0.061505f, 0.052737f,
  -0.084871f, -0.348335f, 0.312544f,  0.120753f,  -0.707222f, -0.010050f,
  -0.137148f, -0.351765f,
};

static const float rect_part_nn_bias_64_layer1[LABELS] = {
  -0.926768f,
  0.765832f,
  0.663683f,
  -0.621865f,
};

static const NN_CONFIG rect_part_nnconfig_64 = {
  FEATURES,  // num_inputs
  LABELS,    // num_outputs
  1,         // num_hidden_layers
  {
      32,
  },  // num_hidden_nodes
  {
      rect_part_nn_weights_64_layer0,
      rect_part_nn_weights_64_layer1,
  },
  {
      rect_part_nn_bias_64_layer0,
      rect_part_nn_bias_64_layer1,
  },
};

static void ml_prune_rect_partition(VP9_COMP *const cpi, MACROBLOCK *const x,
                                    BLOCK_SIZE bsize,
                                    const PC_TREE *const pc_tree,
                                    int *allow_horz, int *allow_vert,
                                    int64_t ref_rd, int mi_row, int mi_col) {
  const NN_CONFIG *nn_config = NULL;
  float score[LABELS] = {
    0.0f,
  };
  int thresh = -1;
  int i;

  if (ref_rd <= 0 || ref_rd > 1000000000) return;

  switch (bsize) {
    case BLOCK_8X8: break;
    case BLOCK_16X16:
      nn_config = &rect_part_nnconfig_16;
      thresh = cpi->sf.ml_prune_rect_partition_threhold[1];
      break;
    case BLOCK_32X32:
      nn_config = &rect_part_nnconfig_32;
      thresh = cpi->sf.ml_prune_rect_partition_threhold[2];
      break;
    case BLOCK_64X64:
      nn_config = &rect_part_nnconfig_64;
      thresh = cpi->sf.ml_prune_rect_partition_threhold[3];
      break;
    default: assert(0 && "Unexpected block size."); return;
  }
  if (!nn_config || thresh < 0) return;

  // Feature extraction and model score calculation.
  {
    const int64_t none_rdcost = pc_tree->none.rdcost;
    const VP9_COMMON *const cm = &cpi->common;
    const int dc_q = vp9_dc_quant(cm->base_qindex, 0, cm->bit_depth);
    int feature_index = 0;
    unsigned int block_var = 0;
    unsigned int sub_block_var[4] = { 0 };
    float features[FEATURES];

    features[feature_index++] =
        (float)(pc_tree->partitioning == PARTITION_NONE);
    features[feature_index++] = logf((float)(dc_q * dc_q) / 256.0f + 1.0f);

    // Calculate source pixel variance.
    {
      struct buf_2d buf;
      const BLOCK_SIZE subsize = get_subsize(bsize, PARTITION_SPLIT);
      const int bs = 4 * num_4x4_blocks_wide_lookup[bsize];
      const MACROBLOCKD *const xd = &x->e_mbd;
      vp9_setup_src_planes(x, cpi->Source, mi_row, mi_col);

      (void)xd;
#if CONFIG_VP9_HIGHBITDEPTH
      if (xd->cur_buf->flags & YV12_FLAG_HIGHBITDEPTH) {
        block_var = vp9_high_get_sby_perpixel_variance(cpi, &x->plane[0].src,
                                                       bsize, xd->bd);
      } else {
        block_var = vp9_get_sby_perpixel_variance(cpi, &x->plane[0].src, bsize);
      }
#else
      block_var = vp9_get_sby_perpixel_variance(cpi, &x->plane[0].src, bsize);
#endif  // CONFIG_VP9_HIGHBITDEPTH

      buf.stride = x->plane[0].src.stride;
      for (i = 0; i < 4; ++i) {
        const int x_idx = (i & 1) * bs / 2;
        const int y_idx = (i >> 1) * bs / 2;
        buf.buf = x->plane[0].src.buf + x_idx + y_idx * buf.stride;
#if CONFIG_VP9_HIGHBITDEPTH
        if (xd->cur_buf->flags & YV12_FLAG_HIGHBITDEPTH) {
          sub_block_var[i] =
              vp9_high_get_sby_perpixel_variance(cpi, &buf, subsize, xd->bd);
        } else {
          sub_block_var[i] = vp9_get_sby_perpixel_variance(cpi, &buf, subsize);
        }
#else
        sub_block_var[i] = vp9_get_sby_perpixel_variance(cpi, &buf, subsize);
#endif  // CONFIG_VP9_HIGHBITDEPTH
      }
    }

    features[feature_index++] = logf((float)block_var + 1.0f);
    features[feature_index++] = logf((float)ref_rd + 1.0f);
    features[feature_index++] = (none_rdcost > 0 && none_rdcost < 1000000000)
                                    ? (float)pc_tree->none.skippable
                                    : 0.0f;

    for (i = 0; i < 4; ++i) {
      const int64_t this_rd = pc_tree->split[i]->none.rdcost;
      const int rd_valid = this_rd > 0 && this_rd < 1000000000;
      // Ratio between sub-block RD and whole block RD.
      features[feature_index++] =
          rd_valid ? ((float)this_rd / (float)ref_rd) : 1.0f;
      // Sub-block skippable.
      features[feature_index++] =
          rd_valid ? ((float)pc_tree->split[i]->none.skippable) : 0.0f;
    }

    {
      const float denom = (float)(block_var + 1);
      const float low_b = 0.1f;
      const float high_b = 10.0f;
      for (i = 0; i < 4; ++i) {
        // Ratio between the quarter sub-block variance and the
        // whole-block variance.
        float var_ratio = (float)(sub_block_var[i] + 1) / denom;
        if (var_ratio < low_b) var_ratio = low_b;
        if (var_ratio > high_b) var_ratio = high_b;
        features[feature_index++] = var_ratio;
      }
    }
    assert(feature_index == FEATURES);
    nn_predict(features, nn_config, score);
  }

  // Make decisions based on the model score.
  {
    int max_score = -1000;
    int horz = 0, vert = 0;
    int int_score[LABELS];
    for (i = 0; i < LABELS; ++i) {
      int_score[i] = (int)(100 * score[i]);
      max_score = VPXMAX(int_score[i], max_score);
    }
    thresh = max_score - thresh;
    for (i = 0; i < LABELS; ++i) {
      if (int_score[i] >= thresh) {
        if ((i >> 0) & 1) horz = 1;
        if ((i >> 1) & 1) vert = 1;
      }
    }
    *allow_horz = *allow_horz && horz;
    *allow_vert = *allow_vert && vert;
  }
}
#undef FEATURES

int get_rdmult_delta(VP9_COMP *cpi, BLOCK_SIZE bsize, int mi_row, int mi_col,
                     int orig_rdmult) {
  TplDepFrame *tpl_frame = &cpi->tpl_stats[cpi->twopass.gf_group.index];
  TplDepStats *tpl_stats = tpl_frame->tpl_stats_ptr;
  int tpl_stride = tpl_frame->stride;
  int64_t intra_cost = 0;
  int64_t mc_dep_cost = 0;
  int mi_wide = num_8x8_blocks_wide_lookup[bsize];
  int mi_high = num_8x8_blocks_high_lookup[bsize];
  int row, col;

  int dr = 0;
  int count = 0;
  double r0, rk, beta;

  if (tpl_frame->is_valid == 0) return orig_rdmult;

  if (cpi->common.show_frame) return orig_rdmult;

  if (cpi->twopass.gf_group.index >= MAX_LAG_BUFFERS) return orig_rdmult;

  for (row = mi_row; row < mi_row + mi_high; ++row) {
    for (col = mi_col; col < mi_col + mi_wide; ++col) {
      TplDepStats *this_stats = &tpl_stats[row * tpl_stride + col];

      if (row >= cpi->common.mi_rows || col >= cpi->common.mi_cols) continue;

      intra_cost += this_stats->intra_cost;
      mc_dep_cost += this_stats->mc_dep_cost;

      ++count;
    }
  }

  vpx_clear_system_state();

  r0 = cpi->rd.r0;
  rk = (double)intra_cost / mc_dep_cost;
  beta = r0 / rk;
  dr = vp9_get_adaptive_rdmult(cpi, beta);

  dr = VPXMIN(dr, orig_rdmult * 3 / 2);
  dr = VPXMAX(dr, orig_rdmult * 1 / 2);

  dr = VPXMAX(1, dr);

  return dr;
}

// TODO(jingning,jimbankoski,rbultje): properly skip partition types that are
// unlikely to be selected depending on previous rate-distortion optimization
// results, for encoding speed-up.
static void rd_pick_partition(VP9_COMP *cpi, ThreadData *td,
                              TileDataEnc *tile_data, TOKENEXTRA **tp,
                              int mi_row, int mi_col, BLOCK_SIZE bsize,
                              RD_COST *rd_cost, int64_t best_rd,
                              PC_TREE *pc_tree) {
  VP9_COMMON *const cm = &cpi->common;
  TileInfo *const tile_info = &tile_data->tile_info;
  MACROBLOCK *const x = &td->mb;
  MACROBLOCKD *const xd = &x->e_mbd;
  const int mi_step = num_8x8_blocks_wide_lookup[bsize] / 2;
  ENTROPY_CONTEXT l[16 * MAX_MB_PLANE], a[16 * MAX_MB_PLANE];
  PARTITION_CONTEXT sl[8], sa[8];
  TOKENEXTRA *tp_orig = *tp;
  PICK_MODE_CONTEXT *const ctx = &pc_tree->none;
  int i;
  const int pl = partition_plane_context(xd, mi_row, mi_col, bsize);
  BLOCK_SIZE subsize;
  RD_COST this_rdc, sum_rdc, best_rdc;
  int do_split = bsize >= BLOCK_8X8;
  int do_rect = 1;
  INTERP_FILTER pred_interp_filter;

  // Override skipping rectangular partition operations for edge blocks
  const int force_horz_split = (mi_row + mi_step >= cm->mi_rows);
  const int force_vert_split = (mi_col + mi_step >= cm->mi_cols);
  const int xss = x->e_mbd.plane[1].subsampling_x;
  const int yss = x->e_mbd.plane[1].subsampling_y;

  BLOCK_SIZE min_size = x->min_partition_size;
  BLOCK_SIZE max_size = x->max_partition_size;

#if CONFIG_FP_MB_STATS
  unsigned int src_diff_var = UINT_MAX;
  int none_complexity = 0;
#endif

  int partition_none_allowed = !force_horz_split && !force_vert_split;
  int partition_horz_allowed =
      !force_vert_split && yss <= xss && bsize >= BLOCK_8X8;
  int partition_vert_allowed =
      !force_horz_split && xss <= yss && bsize >= BLOCK_8X8;

  int64_t dist_breakout_thr = cpi->sf.partition_search_breakout_thr.dist;
  int rate_breakout_thr = cpi->sf.partition_search_breakout_thr.rate;
  int must_split = 0;
  int partition_mul = cpi->sf.enable_tpl_model && cpi->oxcf.aq_mode == NO_AQ
                          ? x->cb_rdmult
                          : cpi->rd.RDMULT;
  // Ref frames picked in the [i_th] quarter subblock during square partition
  // RD search. It may be used to prune ref frame selection of rect partitions.
  uint8_t ref_frames_used[4] = { 0, 0, 0, 0 };

  (void)*tp_orig;

  assert(num_8x8_blocks_wide_lookup[bsize] ==
         num_8x8_blocks_high_lookup[bsize]);

  dist_breakout_thr >>=
      8 - (b_width_log2_lookup[bsize] + b_height_log2_lookup[bsize]);

  rate_breakout_thr *= num_pels_log2_lookup[bsize];

  vp9_rd_cost_init(&this_rdc);
  vp9_rd_cost_init(&sum_rdc);
  vp9_rd_cost_reset(&best_rdc);
  best_rdc.rdcost = best_rd;

  set_offsets(cpi, tile_info, x, mi_row, mi_col, bsize);

  if (bsize == BLOCK_16X16 && cpi->oxcf.aq_mode != NO_AQ &&
      cpi->oxcf.aq_mode != LOOKAHEAD_AQ)
    x->mb_energy = vp9_block_energy(cpi, x, bsize);

  if (cpi->sf.cb_partition_search && bsize == BLOCK_16X16) {
    int cb_partition_search_ctrl =
        ((pc_tree->index == 0 || pc_tree->index == 3) +
         get_chessboard_index(cm->current_video_frame)) &
        0x1;

    if (cb_partition_search_ctrl && bsize > min_size && bsize < max_size)
      set_partition_range(cm, xd, mi_row, mi_col, bsize, &min_size, &max_size);
  }

  // Get sub block energy range
  if (bsize >= BLOCK_16X16) {
    int min_energy, max_energy;
    vp9_get_sub_block_energy(cpi, x, mi_row, mi_col, bsize, &min_energy,
                             &max_energy);
    must_split = (min_energy < -3) && (max_energy - min_energy > 2);
  }

  // Determine partition types in search according to the speed features.
  // The threshold set here has to be of square block size.
  if (cpi->sf.auto_min_max_partition_size) {
    partition_none_allowed &= (bsize <= max_size);
    partition_horz_allowed &=
        ((bsize <= max_size && bsize > min_size) || force_horz_split);
    partition_vert_allowed &=
        ((bsize <= max_size && bsize > min_size) || force_vert_split);
    do_split &= bsize > min_size;
  }

  if (cpi->sf.use_square_partition_only &&
      bsize > cpi->sf.use_square_only_threshold) {
    if (cpi->use_svc) {
      if (!vp9_active_h_edge(cpi, mi_row, mi_step) || x->e_mbd.lossless)
        partition_horz_allowed &= force_horz_split;
      if (!vp9_active_v_edge(cpi, mi_row, mi_step) || x->e_mbd.lossless)
        partition_vert_allowed &= force_vert_split;
    } else {
      partition_horz_allowed &= force_horz_split;
      partition_vert_allowed &= force_vert_split;
    }
  }

  save_context(x, mi_row, mi_col, a, l, sa, sl, bsize);

#if CONFIG_FP_MB_STATS
  if (cpi->use_fp_mb_stats) {
    set_offsets(cpi, tile_info, x, mi_row, mi_col, bsize);
    src_diff_var = get_sby_perpixel_diff_variance(cpi, &x->plane[0].src, mi_row,
                                                  mi_col, bsize);
  }
#endif

#if CONFIG_FP_MB_STATS
  // Decide whether we shall split directly and skip searching NONE by using
  // the first pass block statistics
  if (cpi->use_fp_mb_stats && bsize >= BLOCK_32X32 && do_split &&
      partition_none_allowed && src_diff_var > 4 &&
      cm->base_qindex < qindex_split_threshold_lookup[bsize]) {
    int mb_row = mi_row >> 1;
    int mb_col = mi_col >> 1;
    int mb_row_end =
        VPXMIN(mb_row + num_16x16_blocks_high_lookup[bsize], cm->mb_rows);
    int mb_col_end =
        VPXMIN(mb_col + num_16x16_blocks_wide_lookup[bsize], cm->mb_cols);
    int r, c;

    // compute a complexity measure, basically measure inconsistency of motion
    // vectors obtained from the first pass in the current block
    for (r = mb_row; r < mb_row_end; r++) {
      for (c = mb_col; c < mb_col_end; c++) {
        const int mb_index = r * cm->mb_cols + c;

        MOTION_DIRECTION this_mv;
        MOTION_DIRECTION right_mv;
        MOTION_DIRECTION bottom_mv;

        this_mv =
            get_motion_direction_fp(cpi->twopass.this_frame_mb_stats[mb_index]);

        // to its right
        if (c != mb_col_end - 1) {
          right_mv = get_motion_direction_fp(
              cpi->twopass.this_frame_mb_stats[mb_index + 1]);
          none_complexity += get_motion_inconsistency(this_mv, right_mv);
        }

        // to its bottom
        if (r != mb_row_end - 1) {
          bottom_mv = get_motion_direction_fp(
              cpi->twopass.this_frame_mb_stats[mb_index + cm->mb_cols]);
          none_complexity += get_motion_inconsistency(this_mv, bottom_mv);
        }

        // do not count its left and top neighbors to avoid double counting
      }
    }

    if (none_complexity > complexity_16x16_blocks_threshold[bsize]) {
      partition_none_allowed = 0;
    }
  }
#endif

  // PARTITION_NONE
  if (partition_none_allowed) {
    rd_pick_sb_modes(cpi, tile_data, x, mi_row, mi_col, &this_rdc, bsize, ctx,
                     best_rdc.rdcost);
    ctx->rdcost = this_rdc.rdcost;
    if (this_rdc.rate != INT_MAX) {
      if (cpi->sf.prune_ref_frame_for_rect_partitions) {
        const int ref1 = ctx->mic.ref_frame[0];
        const int ref2 = ctx->mic.ref_frame[1];
        for (i = 0; i < 4; ++i) {
          ref_frames_used[i] |= (1 << ref1);
          if (ref2 > 0) ref_frames_used[i] |= (1 << ref2);
        }
      }
      if (bsize >= BLOCK_8X8) {
        this_rdc.rdcost += RDCOST(partition_mul, x->rddiv,
                                  cpi->partition_cost[pl][PARTITION_NONE], 0);
        this_rdc.rate += cpi->partition_cost[pl][PARTITION_NONE];
      }

      if (this_rdc.rdcost < best_rdc.rdcost) {
        MODE_INFO *mi = xd->mi[0];

        best_rdc = this_rdc;
        if (bsize >= BLOCK_8X8) pc_tree->partitioning = PARTITION_NONE;

        if (cpi->sf.ml_partition_search_early_termination) {
          // Currently, the machine-learning based partition search early
          // termination is only used while bsize is 16x16, 32x32 or 64x64,
          // VPXMIN(cm->width, cm->height) >= 480, and speed = 0.
          if (!x->e_mbd.lossless &&
              !segfeature_active(&cm->seg, mi->segment_id, SEG_LVL_SKIP) &&
              ctx->mic.mode >= INTRA_MODES && bsize >= BLOCK_16X16) {
            if (ml_pruning_partition(cm, xd, ctx, mi_row, mi_col, bsize)) {
              do_split = 0;
              do_rect = 0;
            }
          }
        }

        if ((do_split || do_rect) && !x->e_mbd.lossless && ctx->skippable) {
          int use_ml_based_breakout =
              cpi->sf.use_ml_partition_search_breakout &&
              cm->base_qindex >= 100;
#if CONFIG_VP9_HIGHBITDEPTH
          if (x->e_mbd.cur_buf->flags & YV12_FLAG_HIGHBITDEPTH)
            use_ml_based_breakout = 0;
#endif  // CONFIG_VP9_HIGHBITDEPTH
          if (use_ml_based_breakout) {
            if (ml_predict_breakout(cpi, bsize, x, &this_rdc)) {
              do_split = 0;
              do_rect = 0;
            }
          } else {
            if (!cpi->sf.ml_partition_search_early_termination) {
              if ((best_rdc.dist < (dist_breakout_thr >> 2)) ||
                  (best_rdc.dist < dist_breakout_thr &&
                   best_rdc.rate < rate_breakout_thr)) {
                do_split = 0;
                do_rect = 0;
              }
            }
          }
        }

#if CONFIG_FP_MB_STATS
        // Check if every 16x16 first pass block statistics has zero
        // motion and the corresponding first pass residue is small enough.
        // If that is the case, check the difference variance between the
        // current frame and the last frame. If the variance is small enough,
        // stop further splitting in RD optimization
        if (cpi->use_fp_mb_stats && do_split != 0 &&
            cm->base_qindex > qindex_skip_threshold_lookup[bsize]) {
          int mb_row = mi_row >> 1;
          int mb_col = mi_col >> 1;
          int mb_row_end =
              VPXMIN(mb_row + num_16x16_blocks_high_lookup[bsize], cm->mb_rows);
          int mb_col_end =
              VPXMIN(mb_col + num_16x16_blocks_wide_lookup[bsize], cm->mb_cols);
          int r, c;

          int skip = 1;
          for (r = mb_row; r < mb_row_end; r++) {
            for (c = mb_col; c < mb_col_end; c++) {
              const int mb_index = r * cm->mb_cols + c;
              if (!(cpi->twopass.this_frame_mb_stats[mb_index] &
                    FPMB_MOTION_ZERO_MASK) ||
                  !(cpi->twopass.this_frame_mb_stats[mb_index] &
                    FPMB_ERROR_SMALL_MASK)) {
                skip = 0;
                break;
              }
            }
            if (skip == 0) {
              break;
            }
          }

          if (skip) {
            if (src_diff_var == UINT_MAX) {
              set_offsets(cpi, tile_info, x, mi_row, mi_col, bsize);
              src_diff_var = get_sby_perpixel_diff_variance(
                  cpi, &x->plane[0].src, mi_row, mi_col, bsize);
            }
            if (src_diff_var < 8) {
              do_split = 0;
              do_rect = 0;
            }
          }
        }
#endif
      }
    }
    restore_context(x, mi_row, mi_col, a, l, sa, sl, bsize);
  }

  // store estimated motion vector
  store_pred_mv(x, ctx);

  // If the interp_filter is marked as SWITCHABLE_FILTERS, it was for an
  // intra block and used for context purposes.
  if (ctx->mic.interp_filter == SWITCHABLE_FILTERS) {
    pred_interp_filter = EIGHTTAP;
  } else {
    pred_interp_filter = ctx->mic.interp_filter;
  }

  // PARTITION_SPLIT
  // TODO(jingning): use the motion vectors given by the above search as
  // the starting point of motion search in the following partition type check.
  pc_tree->split[0]->none.rdcost = 0;
  pc_tree->split[1]->none.rdcost = 0;
  pc_tree->split[2]->none.rdcost = 0;
  pc_tree->split[3]->none.rdcost = 0;
  if (do_split || must_split) {
    subsize = get_subsize(bsize, PARTITION_SPLIT);
    load_pred_mv(x, ctx);
    if (bsize == BLOCK_8X8) {
      i = 4;
      if (cpi->sf.adaptive_pred_interp_filter && partition_none_allowed)
        pc_tree->leaf_split[0]->pred_interp_filter = pred_interp_filter;
      rd_pick_sb_modes(cpi, tile_data, x, mi_row, mi_col, &sum_rdc, subsize,
                       pc_tree->leaf_split[0], best_rdc.rdcost);
      if (sum_rdc.rate == INT_MAX) {
        sum_rdc.rdcost = INT64_MAX;
      } else {
        if (cpi->sf.prune_ref_frame_for_rect_partitions) {
          const int ref1 = pc_tree->leaf_split[0]->mic.ref_frame[0];
          const int ref2 = pc_tree->leaf_split[0]->mic.ref_frame[1];
          for (i = 0; i < 4; ++i) {
            ref_frames_used[i] |= (1 << ref1);
            if (ref2 > 0) ref_frames_used[i] |= (1 << ref2);
          }
        }
      }
    } else {
      for (i = 0; (i < 4) && ((sum_rdc.rdcost < best_rdc.rdcost) || must_split);
           ++i) {
        const int x_idx = (i & 1) * mi_step;
        const int y_idx = (i >> 1) * mi_step;

        if (mi_row + y_idx >= cm->mi_rows || mi_col + x_idx >= cm->mi_cols)
          continue;

        pc_tree->split[i]->index = i;
        if (cpi->sf.prune_ref_frame_for_rect_partitions)
          pc_tree->split[i]->none.rate = INT_MAX;
        rd_pick_partition(cpi, td, tile_data, tp, mi_row + y_idx,
                          mi_col + x_idx, subsize, &this_rdc,
                          // A must split test here increases the number of sub
                          // partitions but hurts metrics results quite a bit,
                          // so this extra test is commented out pending
                          // further tests on whether it adds much in terms of
                          // visual quality.
                          // (must_split) ? best_rdc.rdcost
                          //              : best_rdc.rdcost - sum_rdc.rdcost,
                          best_rdc.rdcost - sum_rdc.rdcost, pc_tree->split[i]);

        if (this_rdc.rate == INT_MAX) {
          sum_rdc.rdcost = INT64_MAX;
          break;
        } else {
          if (cpi->sf.prune_ref_frame_for_rect_partitions &&
              pc_tree->split[i]->none.rate != INT_MAX) {
            const int ref1 = pc_tree->split[i]->none.mic.ref_frame[0];
            const int ref2 = pc_tree->split[i]->none.mic.ref_frame[1];
            ref_frames_used[i] |= (1 << ref1);
            if (ref2 > 0) ref_frames_used[i] |= (1 << ref2);
          }
          sum_rdc.rate += this_rdc.rate;
          sum_rdc.dist += this_rdc.dist;
          sum_rdc.rdcost += this_rdc.rdcost;
        }
      }
    }

    if (((sum_rdc.rdcost < best_rdc.rdcost) || must_split) && i == 4) {
      sum_rdc.rdcost += RDCOST(partition_mul, x->rddiv,
                               cpi->partition_cost[pl][PARTITION_SPLIT], 0);
      sum_rdc.rate += cpi->partition_cost[pl][PARTITION_SPLIT];

      if ((sum_rdc.rdcost < best_rdc.rdcost) ||
          (must_split && (sum_rdc.dist < best_rdc.dist))) {
        best_rdc = sum_rdc;
        pc_tree->partitioning = PARTITION_SPLIT;

        // Rate and distortion based partition search termination clause.
        if (!cpi->sf.ml_partition_search_early_termination &&
            !x->e_mbd.lossless &&
            ((best_rdc.dist < (dist_breakout_thr >> 2)) ||
             (best_rdc.dist < dist_breakout_thr &&
              best_rdc.rate < rate_breakout_thr))) {
          do_rect = 0;
        }
      }
    } else {
      // skip rectangular partition test when larger block size
      // gives better rd cost
      if ((cpi->sf.less_rectangular_check) &&
          ((bsize > cpi->sf.use_square_only_threshold) ||
           (best_rdc.dist < dist_breakout_thr)))
        do_rect &= !partition_none_allowed;
    }
    restore_context(x, mi_row, mi_col, a, l, sa, sl, bsize);
  }

  pc_tree->horizontal[0].skip_ref_frame_mask = 0;
  pc_tree->horizontal[1].skip_ref_frame_mask = 0;
  pc_tree->vertical[0].skip_ref_frame_mask = 0;
  pc_tree->vertical[1].skip_ref_frame_mask = 0;
  if (cpi->sf.prune_ref_frame_for_rect_partitions) {
    uint8_t used_frames;
    used_frames = ref_frames_used[0] | ref_frames_used[1];
    if (used_frames) pc_tree->horizontal[0].skip_ref_frame_mask = ~used_frames;
    used_frames = ref_frames_used[2] | ref_frames_used[3];
    if (used_frames) pc_tree->horizontal[1].skip_ref_frame_mask = ~used_frames;
    used_frames = ref_frames_used[0] | ref_frames_used[2];
    if (used_frames) pc_tree->vertical[0].skip_ref_frame_mask = ~used_frames;
    used_frames = ref_frames_used[1] | ref_frames_used[3];
    if (used_frames) pc_tree->vertical[1].skip_ref_frame_mask = ~used_frames;
  }

  {
    int do_ml_rect_partition_pruning =
        !frame_is_intra_only(cm) && !force_horz_split && !force_vert_split &&
        (partition_horz_allowed || partition_vert_allowed) && bsize > BLOCK_8X8;
#if CONFIG_VP9_HIGHBITDEPTH
    if (xd->cur_buf->flags & YV12_FLAG_HIGHBITDEPTH)
      do_ml_rect_partition_pruning = 0;
#endif
    if (do_ml_rect_partition_pruning) {
      ml_prune_rect_partition(cpi, x, bsize, pc_tree, &partition_horz_allowed,
                              &partition_vert_allowed, best_rdc.rdcost, mi_row,
                              mi_col);
    }
  }

  // PARTITION_HORZ
  if (partition_horz_allowed &&
      (do_rect || vp9_active_h_edge(cpi, mi_row, mi_step))) {
    const int part_mode_rate = cpi->partition_cost[pl][PARTITION_HORZ];
    const int64_t part_mode_rdcost =
        RDCOST(partition_mul, x->rddiv, part_mode_rate, 0);
    subsize = get_subsize(bsize, PARTITION_HORZ);
    load_pred_mv(x, ctx);
    if (cpi->sf.adaptive_pred_interp_filter && bsize == BLOCK_8X8 &&
        partition_none_allowed)
      pc_tree->horizontal[0].pred_interp_filter = pred_interp_filter;
    rd_pick_sb_modes(cpi, tile_data, x, mi_row, mi_col, &sum_rdc, subsize,
                     &pc_tree->horizontal[0],
                     best_rdc.rdcost - part_mode_rdcost);
    if (sum_rdc.rdcost < INT64_MAX) {
      sum_rdc.rdcost += part_mode_rdcost;
      sum_rdc.rate += part_mode_rate;
    }

    if (sum_rdc.rdcost < best_rdc.rdcost && mi_row + mi_step < cm->mi_rows &&
        bsize > BLOCK_8X8) {
      PICK_MODE_CONTEXT *ctx = &pc_tree->horizontal[0];
      update_state(cpi, td, ctx, mi_row, mi_col, subsize, 0);
      encode_superblock(cpi, td, tp, 0, mi_row, mi_col, subsize, ctx);
      if (cpi->sf.adaptive_pred_interp_filter && bsize == BLOCK_8X8 &&
          partition_none_allowed)
        pc_tree->horizontal[1].pred_interp_filter = pred_interp_filter;
      rd_pick_sb_modes(cpi, tile_data, x, mi_row + mi_step, mi_col, &this_rdc,
                       subsize, &pc_tree->horizontal[1],
                       best_rdc.rdcost - sum_rdc.rdcost);
      if (this_rdc.rate == INT_MAX) {
        sum_rdc.rdcost = INT64_MAX;
      } else {
        sum_rdc.rate += this_rdc.rate;
        sum_rdc.dist += this_rdc.dist;
        sum_rdc.rdcost += this_rdc.rdcost;
      }
    }

    if (sum_rdc.rdcost < best_rdc.rdcost) {
      best_rdc = sum_rdc;
      pc_tree->partitioning = PARTITION_HORZ;

      if ((cpi->sf.less_rectangular_check) &&
          (bsize > cpi->sf.use_square_only_threshold))
        do_rect = 0;
    }
    restore_context(x, mi_row, mi_col, a, l, sa, sl, bsize);
  }

  // PARTITION_VERT
  if (partition_vert_allowed &&
      (do_rect || vp9_active_v_edge(cpi, mi_col, mi_step))) {
    const int part_mode_rate = cpi->partition_cost[pl][PARTITION_VERT];
    const int64_t part_mode_rdcost =
        RDCOST(partition_mul, x->rddiv, part_mode_rate, 0);
    subsize = get_subsize(bsize, PARTITION_VERT);
    load_pred_mv(x, ctx);
    if (cpi->sf.adaptive_pred_interp_filter && bsize == BLOCK_8X8 &&
        partition_none_allowed)
      pc_tree->vertical[0].pred_interp_filter = pred_interp_filter;
    rd_pick_sb_modes(cpi, tile_data, x, mi_row, mi_col, &sum_rdc, subsize,
                     &pc_tree->vertical[0], best_rdc.rdcost - part_mode_rdcost);
    if (sum_rdc.rdcost < INT64_MAX) {
      sum_rdc.rdcost += part_mode_rdcost;
      sum_rdc.rate += part_mode_rate;
    }

    if (sum_rdc.rdcost < best_rdc.rdcost && mi_col + mi_step < cm->mi_cols &&
        bsize > BLOCK_8X8) {
      update_state(cpi, td, &pc_tree->vertical[0], mi_row, mi_col, subsize, 0);
      encode_superblock(cpi, td, tp, 0, mi_row, mi_col, subsize,
                        &pc_tree->vertical[0]);
      if (cpi->sf.adaptive_pred_interp_filter && bsize == BLOCK_8X8 &&
          partition_none_allowed)
        pc_tree->vertical[1].pred_interp_filter = pred_interp_filter;
      rd_pick_sb_modes(cpi, tile_data, x, mi_row, mi_col + mi_step, &this_rdc,
                       subsize, &pc_tree->vertical[1],
                       best_rdc.rdcost - sum_rdc.rdcost);
      if (this_rdc.rate == INT_MAX) {
        sum_rdc.rdcost = INT64_MAX;
      } else {
        sum_rdc.rate += this_rdc.rate;
        sum_rdc.dist += this_rdc.dist;
        sum_rdc.rdcost += this_rdc.rdcost;
      }
    }

    if (sum_rdc.rdcost < best_rdc.rdcost) {
      best_rdc = sum_rdc;
      pc_tree->partitioning = PARTITION_VERT;
    }
    restore_context(x, mi_row, mi_col, a, l, sa, sl, bsize);
  }

  // TODO(jbb): This code added so that we avoid static analysis
  // warning related to the fact that best_rd isn't used after this
  // point.  This code should be refactored so that the duplicate
  // checks occur in some sub function and thus are used...
  (void)best_rd;
  *rd_cost = best_rdc;

  if (best_rdc.rate < INT_MAX && best_rdc.dist < INT64_MAX &&
      pc_tree->index != 3) {
    int output_enabled = (bsize == BLOCK_64X64);
    encode_sb(cpi, td, tile_info, tp, mi_row, mi_col, output_enabled, bsize,
              pc_tree);
  }

  if (bsize == BLOCK_64X64) {
    assert(tp_orig < *tp);
    assert(best_rdc.rate < INT_MAX);
    assert(best_rdc.dist < INT64_MAX);
  } else {
    assert(tp_orig == *tp);
  }
}

static void encode_rd_sb_row(VP9_COMP *cpi, ThreadData *td,
                             TileDataEnc *tile_data, int mi_row,
                             TOKENEXTRA **tp) {
  VP9_COMMON *const cm = &cpi->common;
  TileInfo *const tile_info = &tile_data->tile_info;
  MACROBLOCK *const x = &td->mb;
  MACROBLOCKD *const xd = &x->e_mbd;
  SPEED_FEATURES *const sf = &cpi->sf;
  const int mi_col_start = tile_info->mi_col_start;
  const int mi_col_end = tile_info->mi_col_end;
  int mi_col;
  const int sb_row = mi_row >> MI_BLOCK_SIZE_LOG2;
  const int num_sb_cols =
      get_num_cols(tile_data->tile_info, MI_BLOCK_SIZE_LOG2);
  int sb_col_in_tile;

  // Initialize the left context for the new SB row
  memset(&xd->left_context, 0, sizeof(xd->left_context));
  memset(xd->left_seg_context, 0, sizeof(xd->left_seg_context));

  // Code each SB in the row
  for (mi_col = mi_col_start, sb_col_in_tile = 0; mi_col < mi_col_end;
       mi_col += MI_BLOCK_SIZE, sb_col_in_tile++) {
    const struct segmentation *const seg = &cm->seg;
    int dummy_rate;
    int64_t dummy_dist;
    RD_COST dummy_rdc;
    int i;
    int seg_skip = 0;

    const int idx_str = cm->mi_stride * mi_row + mi_col;
    MODE_INFO **mi = cm->mi_grid_visible + idx_str;

    (*(cpi->row_mt_sync_read_ptr))(&tile_data->row_mt_sync, sb_row,
                                   sb_col_in_tile);

    if (sf->adaptive_pred_interp_filter) {
      for (i = 0; i < 64; ++i) td->leaf_tree[i].pred_interp_filter = SWITCHABLE;

      for (i = 0; i < 64; ++i) {
        td->pc_tree[i].vertical[0].pred_interp_filter = SWITCHABLE;
        td->pc_tree[i].vertical[1].pred_interp_filter = SWITCHABLE;
        td->pc_tree[i].horizontal[0].pred_interp_filter = SWITCHABLE;
        td->pc_tree[i].horizontal[1].pred_interp_filter = SWITCHABLE;
      }
    }

    for (i = 0; i < MAX_REF_FRAMES; ++i) {
      x->pred_mv[i].row = INT16_MAX;
      x->pred_mv[i].col = INT16_MAX;
    }
    td->pc_root->index = 0;

    if (seg->enabled) {
      const uint8_t *const map =
          seg->update_map ? cpi->segmentation_map : cm->last_frame_seg_map;
      int segment_id = get_segment_id(cm, map, BLOCK_64X64, mi_row, mi_col);
      seg_skip = segfeature_active(seg, segment_id, SEG_LVL_SKIP);
    }

    x->source_variance = UINT_MAX;
    if (sf->partition_search_type == FIXED_PARTITION || seg_skip) {
      const BLOCK_SIZE bsize =
          seg_skip ? BLOCK_64X64 : sf->always_this_block_size;
      set_offsets(cpi, tile_info, x, mi_row, mi_col, BLOCK_64X64);
      set_fixed_partitioning(cpi, tile_info, mi, mi_row, mi_col, bsize);
      rd_use_partition(cpi, td, tile_data, mi, tp, mi_row, mi_col, BLOCK_64X64,
                       &dummy_rate, &dummy_dist, 1, td->pc_root);
    } else if (cpi->partition_search_skippable_frame) {
      BLOCK_SIZE bsize;
      set_offsets(cpi, tile_info, x, mi_row, mi_col, BLOCK_64X64);
      bsize = get_rd_var_based_fixed_partition(cpi, x, mi_row, mi_col);
      set_fixed_partitioning(cpi, tile_info, mi, mi_row, mi_col, bsize);
      rd_use_partition(cpi, td, tile_data, mi, tp, mi_row, mi_col, BLOCK_64X64,
                       &dummy_rate, &dummy_dist, 1, td->pc_root);
    } else if (sf->partition_search_type == VAR_BASED_PARTITION &&
               cm->frame_type != KEY_FRAME) {
      choose_partitioning(cpi, tile_info, x, mi_row, mi_col);
      rd_use_partition(cpi, td, tile_data, mi, tp, mi_row, mi_col, BLOCK_64X64,
                       &dummy_rate, &dummy_dist, 1, td->pc_root);
    } else {
      int orig_rdmult = cpi->rd.RDMULT;
      x->cb_rdmult = orig_rdmult;
      if (cpi->twopass.gf_group.index > 0 && cpi->sf.enable_tpl_model) {
        int dr =
            get_rdmult_delta(cpi, BLOCK_64X64, mi_row, mi_col, orig_rdmult);
        x->cb_rdmult = dr;
      }

      // If required set upper and lower partition size limits
      if (sf->auto_min_max_partition_size) {
        set_offsets(cpi, tile_info, x, mi_row, mi_col, BLOCK_64X64);
        rd_auto_partition_range(cpi, tile_info, xd, mi_row, mi_col,
                                &x->min_partition_size, &x->max_partition_size);
      }
      td->pc_root->none.rdcost = 0;
      rd_pick_partition(cpi, td, tile_data, tp, mi_row, mi_col, BLOCK_64X64,
                        &dummy_rdc, INT64_MAX, td->pc_root);
    }
    (*(cpi->row_mt_sync_write_ptr))(&tile_data->row_mt_sync, sb_row,
                                    sb_col_in_tile, num_sb_cols);
  }
}

static void init_encode_frame_mb_context(VP9_COMP *cpi) {
  MACROBLOCK *const x = &cpi->td.mb;
  VP9_COMMON *const cm = &cpi->common;
  MACROBLOCKD *const xd = &x->e_mbd;
  const int aligned_mi_cols = mi_cols_aligned_to_sb(cm->mi_cols);

  // Copy data over into macro block data structures.
  vp9_setup_src_planes(x, cpi->Source, 0, 0);

  vp9_setup_block_planes(&x->e_mbd, cm->subsampling_x, cm->subsampling_y);

  // Note: this memset assumes above_context[0], [1] and [2]
  // are allocated as part of the same buffer.
  memset(xd->above_context[0], 0,
         sizeof(*xd->above_context[0]) * 2 * aligned_mi_cols * MAX_MB_PLANE);
  memset(xd->above_seg_context, 0,
         sizeof(*xd->above_seg_context) * aligned_mi_cols);
}

static int check_dual_ref_flags(VP9_COMP *cpi) {
  const int ref_flags = cpi->ref_frame_flags;

  if (segfeature_active(&cpi->common.seg, 1, SEG_LVL_REF_FRAME)) {
    return 0;
  } else {
    return (!!(ref_flags & VP9_GOLD_FLAG) + !!(ref_flags & VP9_LAST_FLAG) +
            !!(ref_flags & VP9_ALT_FLAG)) >= 2;
  }
}

static void reset_skip_tx_size(VP9_COMMON *cm, TX_SIZE max_tx_size) {
  int mi_row, mi_col;
  const int mis = cm->mi_stride;
  MODE_INFO **mi_ptr = cm->mi_grid_visible;

  for (mi_row = 0; mi_row < cm->mi_rows; ++mi_row, mi_ptr += mis) {
    for (mi_col = 0; mi_col < cm->mi_cols; ++mi_col) {
      if (mi_ptr[mi_col]->tx_size > max_tx_size)
        mi_ptr[mi_col]->tx_size = max_tx_size;
    }
  }
}

static MV_REFERENCE_FRAME get_frame_type(const VP9_COMP *cpi) {
  if (frame_is_intra_only(&cpi->common))
    return INTRA_FRAME;
  else if (cpi->rc.is_src_frame_alt_ref && cpi->refresh_golden_frame)
    return ALTREF_FRAME;
  else if (cpi->refresh_golden_frame || cpi->refresh_alt_ref_frame)
    return GOLDEN_FRAME;
  else
    return LAST_FRAME;
}

static TX_MODE select_tx_mode(const VP9_COMP *cpi, MACROBLOCKD *const xd) {
  if (xd->lossless) return ONLY_4X4;
  if (cpi->common.frame_type == KEY_FRAME && cpi->sf.use_nonrd_pick_mode)
    return ALLOW_16X16;
  if (cpi->sf.tx_size_search_method == USE_LARGESTALL)
    return ALLOW_32X32;
  else if (cpi->sf.tx_size_search_method == USE_FULL_RD ||
           cpi->sf.tx_size_search_method == USE_TX_8X8)
    return TX_MODE_SELECT;
  else
    return cpi->common.tx_mode;
}

static void hybrid_intra_mode_search(VP9_COMP *cpi, MACROBLOCK *const x,
                                     RD_COST *rd_cost, BLOCK_SIZE bsize,
                                     PICK_MODE_CONTEXT *ctx) {
  if (!cpi->sf.nonrd_keyframe && bsize < BLOCK_16X16)
    vp9_rd_pick_intra_mode_sb(cpi, x, rd_cost, bsize, ctx, INT64_MAX);
  else
    vp9_pick_intra_mode(cpi, x, rd_cost, bsize, ctx);
}

static void hybrid_search_svc_baseiskey(VP9_COMP *cpi, MACROBLOCK *const x,
                                        RD_COST *rd_cost, BLOCK_SIZE bsize,
                                        PICK_MODE_CONTEXT *ctx,
                                        TileDataEnc *tile_data, int mi_row,
                                        int mi_col) {
  if (!cpi->sf.nonrd_keyframe && bsize <= BLOCK_8X8) {
    vp9_rd_pick_intra_mode_sb(cpi, x, rd_cost, bsize, ctx, INT64_MAX);
  } else {
    if (cpi->svc.disable_inter_layer_pred == INTER_LAYER_PRED_OFF)
      vp9_pick_intra_mode(cpi, x, rd_cost, bsize, ctx);
    else if (bsize >= BLOCK_8X8)
      vp9_pick_inter_mode(cpi, x, tile_data, mi_row, mi_col, rd_cost, bsize,
                          ctx);
    else
      vp9_pick_inter_mode_sub8x8(cpi, x, mi_row, mi_col, rd_cost, bsize, ctx);
  }
}

static void hybrid_search_scene_change(VP9_COMP *cpi, MACROBLOCK *const x,
                                       RD_COST *rd_cost, BLOCK_SIZE bsize,
                                       PICK_MODE_CONTEXT *ctx,
                                       TileDataEnc *tile_data, int mi_row,
                                       int mi_col) {
  if (!cpi->sf.nonrd_keyframe && bsize <= BLOCK_8X8) {
    vp9_rd_pick_intra_mode_sb(cpi, x, rd_cost, bsize, ctx, INT64_MAX);
  } else {
    vp9_pick_inter_mode(cpi, x, tile_data, mi_row, mi_col, rd_cost, bsize, ctx);
  }
}

static void nonrd_pick_sb_modes(VP9_COMP *cpi, TileDataEnc *tile_data,
                                MACROBLOCK *const x, int mi_row, int mi_col,
                                RD_COST *rd_cost, BLOCK_SIZE bsize,
                                PICK_MODE_CONTEXT *ctx) {
  VP9_COMMON *const cm = &cpi->common;
  TileInfo *const tile_info = &tile_data->tile_info;
  MACROBLOCKD *const xd = &x->e_mbd;
  MODE_INFO *mi;
  ENTROPY_CONTEXT l[16 * MAX_MB_PLANE], a[16 * MAX_MB_PLANE];
  BLOCK_SIZE bs = VPXMAX(bsize, BLOCK_8X8);  // processing unit block size
  const int num_4x4_blocks_wide = num_4x4_blocks_wide_lookup[bs];
  const int num_4x4_blocks_high = num_4x4_blocks_high_lookup[bs];
  int plane;

  set_offsets(cpi, tile_info, x, mi_row, mi_col, bsize);
  mi = xd->mi[0];
  mi->sb_type = bsize;

  for (plane = 0; plane < MAX_MB_PLANE; ++plane) {
    struct macroblockd_plane *pd = &xd->plane[plane];
    memcpy(a + num_4x4_blocks_wide * plane, pd->above_context,
           (sizeof(a[0]) * num_4x4_blocks_wide) >> pd->subsampling_x);
    memcpy(l + num_4x4_blocks_high * plane, pd->left_context,
           (sizeof(l[0]) * num_4x4_blocks_high) >> pd->subsampling_y);
  }

  if (cpi->oxcf.aq_mode == CYCLIC_REFRESH_AQ && cm->seg.enabled)
    if (cyclic_refresh_segment_id_boosted(mi->segment_id))
      x->rdmult = vp9_cyclic_refresh_get_rdmult(cpi->cyclic_refresh);

  if (frame_is_intra_only(cm))
    hybrid_intra_mode_search(cpi, x, rd_cost, bsize, ctx);
  else if (cpi->svc.layer_context[cpi->svc.temporal_layer_id].is_key_frame)
    hybrid_search_svc_baseiskey(cpi, x, rd_cost, bsize, ctx, tile_data, mi_row,
                                mi_col);
  else if (segfeature_active(&cm->seg, mi->segment_id, SEG_LVL_SKIP))
    set_mode_info_seg_skip(x, cm->tx_mode, rd_cost, bsize);
  else if (bsize >= BLOCK_8X8) {
    if (cpi->rc.hybrid_intra_scene_change)
      hybrid_search_scene_change(cpi, x, rd_cost, bsize, ctx, tile_data, mi_row,
                                 mi_col);
    else
      vp9_pick_inter_mode(cpi, x, tile_data, mi_row, mi_col, rd_cost, bsize,
                          ctx);
  } else {
    vp9_pick_inter_mode_sub8x8(cpi, x, mi_row, mi_col, rd_cost, bsize, ctx);
  }

  duplicate_mode_info_in_sb(cm, xd, mi_row, mi_col, bsize);

  for (plane = 0; plane < MAX_MB_PLANE; ++plane) {
    struct macroblockd_plane *pd = &xd->plane[plane];
    memcpy(pd->above_context, a + num_4x4_blocks_wide * plane,
           (sizeof(a[0]) * num_4x4_blocks_wide) >> pd->subsampling_x);
    memcpy(pd->left_context, l + num_4x4_blocks_high * plane,
           (sizeof(l[0]) * num_4x4_blocks_high) >> pd->subsampling_y);
  }

  if (rd_cost->rate == INT_MAX) vp9_rd_cost_reset(rd_cost);

  ctx->rate = rd_cost->rate;
  ctx->dist = rd_cost->dist;
}

static void fill_mode_info_sb(VP9_COMMON *cm, MACROBLOCK *x, int mi_row,
                              int mi_col, BLOCK_SIZE bsize, PC_TREE *pc_tree) {
  MACROBLOCKD *xd = &x->e_mbd;
  int bsl = b_width_log2_lookup[bsize], hbs = (1 << bsl) / 4;
  PARTITION_TYPE partition = pc_tree->partitioning;
  BLOCK_SIZE subsize = get_subsize(bsize, partition);

  assert(bsize >= BLOCK_8X8);

  if (mi_row >= cm->mi_rows || mi_col >= cm->mi_cols) return;

  switch (partition) {
    case PARTITION_NONE:
      set_mode_info_offsets(cm, x, xd, mi_row, mi_col);
      *(xd->mi[0]) = pc_tree->none.mic;
      *(x->mbmi_ext) = pc_tree->none.mbmi_ext;
      duplicate_mode_info_in_sb(cm, xd, mi_row, mi_col, bsize);
      break;
    case PARTITION_VERT:
      set_mode_info_offsets(cm, x, xd, mi_row, mi_col);
      *(xd->mi[0]) = pc_tree->vertical[0].mic;
      *(x->mbmi_ext) = pc_tree->vertical[0].mbmi_ext;
      duplicate_mode_info_in_sb(cm, xd, mi_row, mi_col, subsize);

      if (mi_col + hbs < cm->mi_cols) {
        set_mode_info_offsets(cm, x, xd, mi_row, mi_col + hbs);
        *(xd->mi[0]) = pc_tree->vertical[1].mic;
        *(x->mbmi_ext) = pc_tree->vertical[1].mbmi_ext;
        duplicate_mode_info_in_sb(cm, xd, mi_row, mi_col + hbs, subsize);
      }
      break;
    case PARTITION_HORZ:
      set_mode_info_offsets(cm, x, xd, mi_row, mi_col);
      *(xd->mi[0]) = pc_tree->horizontal[0].mic;
      *(x->mbmi_ext) = pc_tree->horizontal[0].mbmi_ext;
      duplicate_mode_info_in_sb(cm, xd, mi_row, mi_col, subsize);
      if (mi_row + hbs < cm->mi_rows) {
        set_mode_info_offsets(cm, x, xd, mi_row + hbs, mi_col);
        *(xd->mi[0]) = pc_tree->horizontal[1].mic;
        *(x->mbmi_ext) = pc_tree->horizontal[1].mbmi_ext;
        duplicate_mode_info_in_sb(cm, xd, mi_row + hbs, mi_col, subsize);
      }
      break;
    case PARTITION_SPLIT: {
      fill_mode_info_sb(cm, x, mi_row, mi_col, subsize, pc_tree->split[0]);
      fill_mode_info_sb(cm, x, mi_row, mi_col + hbs, subsize,
                        pc_tree->split[1]);
      fill_mode_info_sb(cm, x, mi_row + hbs, mi_col, subsize,
                        pc_tree->split[2]);
      fill_mode_info_sb(cm, x, mi_row + hbs, mi_col + hbs, subsize,
                        pc_tree->split[3]);
      break;
    }
    default: break;
  }
}

// Reset the prediction pixel ready flag recursively.
static void pred_pixel_ready_reset(PC_TREE *pc_tree, BLOCK_SIZE bsize) {
  pc_tree->none.pred_pixel_ready = 0;
  pc_tree->horizontal[0].pred_pixel_ready = 0;
  pc_tree->horizontal[1].pred_pixel_ready = 0;
  pc_tree->vertical[0].pred_pixel_ready = 0;
  pc_tree->vertical[1].pred_pixel_ready = 0;

  if (bsize > BLOCK_8X8) {
    BLOCK_SIZE subsize = get_subsize(bsize, PARTITION_SPLIT);
    int i;
    for (i = 0; i < 4; ++i) pred_pixel_ready_reset(pc_tree->split[i], subsize);
  }
}

static void nonrd_pick_partition(VP9_COMP *cpi, ThreadData *td,
                                 TileDataEnc *tile_data, TOKENEXTRA **tp,
                                 int mi_row, int mi_col, BLOCK_SIZE bsize,
                                 RD_COST *rd_cost, int do_recon,
                                 int64_t best_rd, PC_TREE *pc_tree) {
  const SPEED_FEATURES *const sf = &cpi->sf;
  VP9_COMMON *const cm = &cpi->common;
  TileInfo *const tile_info = &tile_data->tile_info;
  MACROBLOCK *const x = &td->mb;
  MACROBLOCKD *const xd = &x->e_mbd;
  const int ms = num_8x8_blocks_wide_lookup[bsize] / 2;
  TOKENEXTRA *tp_orig = *tp;
  PICK_MODE_CONTEXT *ctx = &pc_tree->none;
  int i;
  BLOCK_SIZE subsize = bsize;
  RD_COST this_rdc, sum_rdc, best_rdc;
  int do_split = bsize >= BLOCK_8X8;
  int do_rect = 1;
  // Override skipping rectangular partition operations for edge blocks
  const int force_horz_split = (mi_row + ms >= cm->mi_rows);
  const int force_vert_split = (mi_col + ms >= cm->mi_cols);
  const int xss = x->e_mbd.plane[1].subsampling_x;
  const int yss = x->e_mbd.plane[1].subsampling_y;

  int partition_none_allowed = !force_horz_split && !force_vert_split;
  int partition_horz_allowed =
      !force_vert_split && yss <= xss && bsize >= BLOCK_8X8;
  int partition_vert_allowed =
      !force_horz_split && xss <= yss && bsize >= BLOCK_8X8;
  (void)*tp_orig;

  // Avoid checking for rectangular partitions for speed >= 6.
  if (cpi->oxcf.speed >= 6) do_rect = 0;

  assert(num_8x8_blocks_wide_lookup[bsize] ==
         num_8x8_blocks_high_lookup[bsize]);

  vp9_rd_cost_init(&sum_rdc);
  vp9_rd_cost_reset(&best_rdc);
  best_rdc.rdcost = best_rd;

  // Determine partition types in search according to the speed features.
  // The threshold set here has to be of square block size.
  if (sf->auto_min_max_partition_size) {
    partition_none_allowed &=
        (bsize <= x->max_partition_size && bsize >= x->min_partition_size);
    partition_horz_allowed &=
        ((bsize <= x->max_partition_size && bsize > x->min_partition_size) ||
         force_horz_split);
    partition_vert_allowed &=
        ((bsize <= x->max_partition_size && bsize > x->min_partition_size) ||
         force_vert_split);
    do_split &= bsize > x->min_partition_size;
  }
  if (sf->use_square_partition_only) {
    partition_horz_allowed &= force_horz_split;
    partition_vert_allowed &= force_vert_split;
  }

  ctx->pred_pixel_ready =
      !(partition_vert_allowed || partition_horz_allowed || do_split);

  // PARTITION_NONE
  if (partition_none_allowed) {
    nonrd_pick_sb_modes(cpi, tile_data, x, mi_row, mi_col, &this_rdc, bsize,
                        ctx);
    ctx->mic = *xd->mi[0];
    ctx->mbmi_ext = *x->mbmi_ext;
    ctx->skip_txfm[0] = x->skip_txfm[0];
    ctx->skip = x->skip;

    if (this_rdc.rate != INT_MAX) {
      int pl = partition_plane_context(xd, mi_row, mi_col, bsize);
      this_rdc.rate += cpi->partition_cost[pl][PARTITION_NONE];
      this_rdc.rdcost =
          RDCOST(x->rdmult, x->rddiv, this_rdc.rate, this_rdc.dist);
      if (this_rdc.rdcost < best_rdc.rdcost) {
        int64_t dist_breakout_thr = sf->partition_search_breakout_thr.dist;
        int64_t rate_breakout_thr = sf->partition_search_breakout_thr.rate;

        dist_breakout_thr >>=
            8 - (b_width_log2_lookup[bsize] + b_height_log2_lookup[bsize]);

        rate_breakout_thr *= num_pels_log2_lookup[bsize];

        best_rdc = this_rdc;
        if (bsize >= BLOCK_8X8) pc_tree->partitioning = PARTITION_NONE;

        if (!x->e_mbd.lossless && this_rdc.rate < rate_breakout_thr &&
            this_rdc.dist < dist_breakout_thr) {
          do_split = 0;
          do_rect = 0;
        }
      }
    }
  }

  // store estimated motion vector
  store_pred_mv(x, ctx);

  // PARTITION_SPLIT
  if (do_split) {
    int pl = partition_plane_context(xd, mi_row, mi_col, bsize);
    sum_rdc.rate += cpi->partition_cost[pl][PARTITION_SPLIT];
    sum_rdc.rdcost = RDCOST(x->rdmult, x->rddiv, sum_rdc.rate, sum_rdc.dist);
    subsize = get_subsize(bsize, PARTITION_SPLIT);
    for (i = 0; i < 4 && sum_rdc.rdcost < best_rdc.rdcost; ++i) {
      const int x_idx = (i & 1) * ms;
      const int y_idx = (i >> 1) * ms;

      if (mi_row + y_idx >= cm->mi_rows || mi_col + x_idx >= cm->mi_cols)
        continue;
      load_pred_mv(x, ctx);
      nonrd_pick_partition(cpi, td, tile_data, tp, mi_row + y_idx,
                           mi_col + x_idx, subsize, &this_rdc, 0,
                           best_rdc.rdcost - sum_rdc.rdcost, pc_tree->split[i]);

      if (this_rdc.rate == INT_MAX) {
        vp9_rd_cost_reset(&sum_rdc);
      } else {
        sum_rdc.rate += this_rdc.rate;
        sum_rdc.dist += this_rdc.dist;
        sum_rdc.rdcost += this_rdc.rdcost;
      }
    }

    if (sum_rdc.rdcost < best_rdc.rdcost) {
      best_rdc = sum_rdc;
      pc_tree->partitioning = PARTITION_SPLIT;
    } else {
      // skip rectangular partition test when larger block size
      // gives better rd cost
      if (sf->less_rectangular_check) do_rect &= !partition_none_allowed;
    }
  }

  // PARTITION_HORZ
  if (partition_horz_allowed && do_rect) {
    subsize = get_subsize(bsize, PARTITION_HORZ);
    load_pred_mv(x, ctx);
    pc_tree->horizontal[0].pred_pixel_ready = 1;
    nonrd_pick_sb_modes(cpi, tile_data, x, mi_row, mi_col, &sum_rdc, subsize,
                        &pc_tree->horizontal[0]);

    pc_tree->horizontal[0].mic = *xd->mi[0];
    pc_tree->horizontal[0].mbmi_ext = *x->mbmi_ext;
    pc_tree->horizontal[0].skip_txfm[0] = x->skip_txfm[0];
    pc_tree->horizontal[0].skip = x->skip;

    if (sum_rdc.rdcost < best_rdc.rdcost && mi_row + ms < cm->mi_rows) {
      load_pred_mv(x, ctx);
      pc_tree->horizontal[1].pred_pixel_ready = 1;
      nonrd_pick_sb_modes(cpi, tile_data, x, mi_row + ms, mi_col, &this_rdc,
                          subsize, &pc_tree->horizontal[1]);

      pc_tree->horizontal[1].mic = *xd->mi[0];
      pc_tree->horizontal[1].mbmi_ext = *x->mbmi_ext;
      pc_tree->horizontal[1].skip_txfm[0] = x->skip_txfm[0];
      pc_tree->horizontal[1].skip = x->skip;

      if (this_rdc.rate == INT_MAX) {
        vp9_rd_cost_reset(&sum_rdc);
      } else {
        int pl = partition_plane_context(xd, mi_row, mi_col, bsize);
        this_rdc.rate += cpi->partition_cost[pl][PARTITION_HORZ];
        sum_rdc.rate += this_rdc.rate;
        sum_rdc.dist += this_rdc.dist;
        sum_rdc.rdcost =
            RDCOST(x->rdmult, x->rddiv, sum_rdc.rate, sum_rdc.dist);
      }
    }

    if (sum_rdc.rdcost < best_rdc.rdcost) {
      best_rdc = sum_rdc;
      pc_tree->partitioning = PARTITION_HORZ;
    } else {
      pred_pixel_ready_reset(pc_tree, bsize);
    }
  }

  // PARTITION_VERT
  if (partition_vert_allowed && do_rect) {
    subsize = get_subsize(bsize, PARTITION_VERT);
    load_pred_mv(x, ctx);
    pc_tree->vertical[0].pred_pixel_ready = 1;
    nonrd_pick_sb_modes(cpi, tile_data, x, mi_row, mi_col, &sum_rdc, subsize,
                        &pc_tree->vertical[0]);
    pc_tree->vertical[0].mic = *xd->mi[0];
    pc_tree->vertical[0].mbmi_ext = *x->mbmi_ext;
    pc_tree->vertical[0].skip_txfm[0] = x->skip_txfm[0];
    pc_tree->vertical[0].skip = x->skip;

    if (sum_rdc.rdcost < best_rdc.rdcost && mi_col + ms < cm->mi_cols) {
      load_pred_mv(x, ctx);
      pc_tree->vertical[1].pred_pixel_ready = 1;
      nonrd_pick_sb_modes(cpi, tile_data, x, mi_row, mi_col + ms, &this_rdc,
                          subsize, &pc_tree->vertical[1]);
      pc_tree->vertical[1].mic = *xd->mi[0];
      pc_tree->vertical[1].mbmi_ext = *x->mbmi_ext;
      pc_tree->vertical[1].skip_txfm[0] = x->skip_txfm[0];
      pc_tree->vertical[1].skip = x->skip;

      if (this_rdc.rate == INT_MAX) {
        vp9_rd_cost_reset(&sum_rdc);
      } else {
        int pl = partition_plane_context(xd, mi_row, mi_col, bsize);
        sum_rdc.rate += cpi->partition_cost[pl][PARTITION_VERT];
        sum_rdc.rate += this_rdc.rate;
        sum_rdc.dist += this_rdc.dist;
        sum_rdc.rdcost =
            RDCOST(x->rdmult, x->rddiv, sum_rdc.rate, sum_rdc.dist);
      }
    }

    if (sum_rdc.rdcost < best_rdc.rdcost) {
      best_rdc = sum_rdc;
      pc_tree->partitioning = PARTITION_VERT;
    } else {
      pred_pixel_ready_reset(pc_tree, bsize);
    }
  }

  *rd_cost = best_rdc;

  if (best_rdc.rate == INT_MAX) {
    vp9_rd_cost_reset(rd_cost);
    return;
  }

  // update mode info array
  fill_mode_info_sb(cm, x, mi_row, mi_col, bsize, pc_tree);

  if (best_rdc.rate < INT_MAX && best_rdc.dist < INT64_MAX && do_recon) {
    int output_enabled = (bsize == BLOCK_64X64);
    encode_sb_rt(cpi, td, tile_info, tp, mi_row, mi_col, output_enabled, bsize,
                 pc_tree);
  }

  if (bsize == BLOCK_64X64 && do_recon) {
    assert(tp_orig < *tp);
    assert(best_rdc.rate < INT_MAX);
    assert(best_rdc.dist < INT64_MAX);
  } else {
    assert(tp_orig == *tp);
  }
}

static void nonrd_select_partition(VP9_COMP *cpi, ThreadData *td,
                                   TileDataEnc *tile_data, MODE_INFO **mi,
                                   TOKENEXTRA **tp, int mi_row, int mi_col,
                                   BLOCK_SIZE bsize, int output_enabled,
                                   RD_COST *rd_cost, PC_TREE *pc_tree) {
  VP9_COMMON *const cm = &cpi->common;
  TileInfo *const tile_info = &tile_data->tile_info;
  MACROBLOCK *const x = &td->mb;
  MACROBLOCKD *const xd = &x->e_mbd;
  const int bsl = b_width_log2_lookup[bsize], hbs = (1 << bsl) / 4;
  const int mis = cm->mi_stride;
  PARTITION_TYPE partition;
  BLOCK_SIZE subsize;
  RD_COST this_rdc;
  BLOCK_SIZE subsize_ref =
      (cpi->sf.adapt_partition_source_sad) ? BLOCK_8X8 : BLOCK_16X16;

  vp9_rd_cost_reset(&this_rdc);
  if (mi_row >= cm->mi_rows || mi_col >= cm->mi_cols) return;

  subsize = (bsize >= BLOCK_8X8) ? mi[0]->sb_type : BLOCK_4X4;
  partition = partition_lookup[bsl][subsize];

  if (bsize == BLOCK_32X32 && subsize == BLOCK_32X32) {
    x->max_partition_size = BLOCK_32X32;
    x->min_partition_size = BLOCK_16X16;
    nonrd_pick_partition(cpi, td, tile_data, tp, mi_row, mi_col, bsize, rd_cost,
                         0, INT64_MAX, pc_tree);
  } else if (bsize == BLOCK_32X32 && partition != PARTITION_NONE &&
             subsize >= subsize_ref) {
    x->max_partition_size = BLOCK_32X32;
    x->min_partition_size = BLOCK_8X8;
    nonrd_pick_partition(cpi, td, tile_data, tp, mi_row, mi_col, bsize, rd_cost,
                         0, INT64_MAX, pc_tree);
  } else if (bsize == BLOCK_16X16 && partition != PARTITION_NONE) {
    x->max_partition_size = BLOCK_16X16;
    x->min_partition_size = BLOCK_8X8;
    nonrd_pick_partition(cpi, td, tile_data, tp, mi_row, mi_col, bsize, rd_cost,
                         0, INT64_MAX, pc_tree);
  } else {
    switch (partition) {
      case PARTITION_NONE:
        pc_tree->none.pred_pixel_ready = 1;
        nonrd_pick_sb_modes(cpi, tile_data, x, mi_row, mi_col, rd_cost, subsize,
                            &pc_tree->none);
        pc_tree->none.mic = *xd->mi[0];
        pc_tree->none.mbmi_ext = *x->mbmi_ext;
        pc_tree->none.skip_txfm[0] = x->skip_txfm[0];
        pc_tree->none.skip = x->skip;
        break;
      case PARTITION_VERT:
        pc_tree->vertical[0].pred_pixel_ready = 1;
        nonrd_pick_sb_modes(cpi, tile_data, x, mi_row, mi_col, rd_cost, subsize,
                            &pc_tree->vertical[0]);
        pc_tree->vertical[0].mic = *xd->mi[0];
        pc_tree->vertical[0].mbmi_ext = *x->mbmi_ext;
        pc_tree->vertical[0].skip_txfm[0] = x->skip_txfm[0];
        pc_tree->vertical[0].skip = x->skip;
        if (mi_col + hbs < cm->mi_cols) {
          pc_tree->vertical[1].pred_pixel_ready = 1;
          nonrd_pick_sb_modes(cpi, tile_data, x, mi_row, mi_col + hbs,
                              &this_rdc, subsize, &pc_tree->vertical[1]);
          pc_tree->vertical[1].mic = *xd->mi[0];
          pc_tree->vertical[1].mbmi_ext = *x->mbmi_ext;
          pc_tree->vertical[1].skip_txfm[0] = x->skip_txfm[0];
          pc_tree->vertical[1].skip = x->skip;
          if (this_rdc.rate != INT_MAX && this_rdc.dist != INT64_MAX &&
              rd_cost->rate != INT_MAX && rd_cost->dist != INT64_MAX) {
            rd_cost->rate += this_rdc.rate;
            rd_cost->dist += this_rdc.dist;
          }
        }
        break;
      case PARTITION_HORZ:
        pc_tree->horizontal[0].pred_pixel_ready = 1;
        nonrd_pick_sb_modes(cpi, tile_data, x, mi_row, mi_col, rd_cost, subsize,
                            &pc_tree->horizontal[0]);
        pc_tree->horizontal[0].mic = *xd->mi[0];
        pc_tree->horizontal[0].mbmi_ext = *x->mbmi_ext;
        pc_tree->horizontal[0].skip_txfm[0] = x->skip_txfm[0];
        pc_tree->horizontal[0].skip = x->skip;
        if (mi_row + hbs < cm->mi_rows) {
          pc_tree->horizontal[1].pred_pixel_ready = 1;
          nonrd_pick_sb_modes(cpi, tile_data, x, mi_row + hbs, mi_col,
                              &this_rdc, subsize, &pc_tree->horizontal[1]);
          pc_tree->horizontal[1].mic = *xd->mi[0];
          pc_tree->horizontal[1].mbmi_ext = *x->mbmi_ext;
          pc_tree->horizontal[1].skip_txfm[0] = x->skip_txfm[0];
          pc_tree->horizontal[1].skip = x->skip;
          if (this_rdc.rate != INT_MAX && this_rdc.dist != INT64_MAX &&
              rd_cost->rate != INT_MAX && rd_cost->dist != INT64_MAX) {
            rd_cost->rate += this_rdc.rate;
            rd_cost->dist += this_rdc.dist;
          }
        }
        break;
      default:
        assert(partition == PARTITION_SPLIT);
        subsize = get_subsize(bsize, PARTITION_SPLIT);
        nonrd_select_partition(cpi, td, tile_data, mi, tp, mi_row, mi_col,
                               subsize, output_enabled, rd_cost,
                               pc_tree->split[0]);
        nonrd_select_partition(cpi, td, tile_data, mi + hbs, tp, mi_row,
                               mi_col + hbs, subsize, output_enabled, &this_rdc,
                               pc_tree->split[1]);
        if (this_rdc.rate != INT_MAX && this_rdc.dist != INT64_MAX &&
            rd_cost->rate != INT_MAX && rd_cost->dist != INT64_MAX) {
          rd_cost->rate += this_rdc.rate;
          rd_cost->dist += this_rdc.dist;
        }
        nonrd_select_partition(cpi, td, tile_data, mi + hbs * mis, tp,
                               mi_row + hbs, mi_col, subsize, output_enabled,
                               &this_rdc, pc_tree->split[2]);
        if (this_rdc.rate != INT_MAX && this_rdc.dist != INT64_MAX &&
            rd_cost->rate != INT_MAX && rd_cost->dist != INT64_MAX) {
          rd_cost->rate += this_rdc.rate;
          rd_cost->dist += this_rdc.dist;
        }
        nonrd_select_partition(cpi, td, tile_data, mi + hbs * mis + hbs, tp,
                               mi_row + hbs, mi_col + hbs, subsize,
                               output_enabled, &this_rdc, pc_tree->split[3]);
        if (this_rdc.rate != INT_MAX && this_rdc.dist != INT64_MAX &&
            rd_cost->rate != INT_MAX && rd_cost->dist != INT64_MAX) {
          rd_cost->rate += this_rdc.rate;
          rd_cost->dist += this_rdc.dist;
        }
        break;
    }
  }

  if (bsize == BLOCK_64X64 && output_enabled)
    encode_sb_rt(cpi, td, tile_info, tp, mi_row, mi_col, 1, bsize, pc_tree);
}

static void nonrd_use_partition(VP9_COMP *cpi, ThreadData *td,
                                TileDataEnc *tile_data, MODE_INFO **mi,
                                TOKENEXTRA **tp, int mi_row, int mi_col,
                                BLOCK_SIZE bsize, int output_enabled,
                                RD_COST *dummy_cost, PC_TREE *pc_tree) {
  VP9_COMMON *const cm = &cpi->common;
  TileInfo *tile_info = &tile_data->tile_info;
  MACROBLOCK *const x = &td->mb;
  MACROBLOCKD *const xd = &x->e_mbd;
  const int bsl = b_width_log2_lookup[bsize], hbs = (1 << bsl) / 4;
  const int mis = cm->mi_stride;
  PARTITION_TYPE partition;
  BLOCK_SIZE subsize;

  if (mi_row >= cm->mi_rows || mi_col >= cm->mi_cols) return;

  subsize = (bsize >= BLOCK_8X8) ? mi[0]->sb_type : BLOCK_4X4;
  partition = partition_lookup[bsl][subsize];

  if (output_enabled && bsize != BLOCK_4X4) {
    int ctx = partition_plane_context(xd, mi_row, mi_col, bsize);
    td->counts->partition[ctx][partition]++;
  }

  switch (partition) {
    case PARTITION_NONE:
      pc_tree->none.pred_pixel_ready = 1;
      nonrd_pick_sb_modes(cpi, tile_data, x, mi_row, mi_col, dummy_cost,
                          subsize, &pc_tree->none);
      pc_tree->none.mic = *xd->mi[0];
      pc_tree->none.mbmi_ext = *x->mbmi_ext;
      pc_tree->none.skip_txfm[0] = x->skip_txfm[0];
      pc_tree->none.skip = x->skip;
      encode_b_rt(cpi, td, tile_info, tp, mi_row, mi_col, output_enabled,
                  subsize, &pc_tree->none);
      break;
    case PARTITION_VERT:
      pc_tree->vertical[0].pred_pixel_ready = 1;
      nonrd_pick_sb_modes(cpi, tile_data, x, mi_row, mi_col, dummy_cost,
                          subsize, &pc_tree->vertical[0]);
      pc_tree->vertical[0].mic = *xd->mi[0];
      pc_tree->vertical[0].mbmi_ext = *x->mbmi_ext;
      pc_tree->vertical[0].skip_txfm[0] = x->skip_txfm[0];
      pc_tree->vertical[0].skip = x->skip;
      encode_b_rt(cpi, td, tile_info, tp, mi_row, mi_col, output_enabled,
                  subsize, &pc_tree->vertical[0]);
      if (mi_col + hbs < cm->mi_cols && bsize > BLOCK_8X8) {
        pc_tree->vertical[1].pred_pixel_ready = 1;
        nonrd_pick_sb_modes(cpi, tile_data, x, mi_row, mi_col + hbs, dummy_cost,
                            subsize, &pc_tree->vertical[1]);
        pc_tree->vertical[1].mic = *xd->mi[0];
        pc_tree->vertical[1].mbmi_ext = *x->mbmi_ext;
        pc_tree->vertical[1].skip_txfm[0] = x->skip_txfm[0];
        pc_tree->vertical[1].skip = x->skip;
        encode_b_rt(cpi, td, tile_info, tp, mi_row, mi_col + hbs,
                    output_enabled, subsize, &pc_tree->vertical[1]);
      }
      break;
    case PARTITION_HORZ:
      pc_tree->horizontal[0].pred_pixel_ready = 1;
      nonrd_pick_sb_modes(cpi, tile_data, x, mi_row, mi_col, dummy_cost,
                          subsize, &pc_tree->horizontal[0]);
      pc_tree->horizontal[0].mic = *xd->mi[0];
      pc_tree->horizontal[0].mbmi_ext = *x->mbmi_ext;
      pc_tree->horizontal[0].skip_txfm[0] = x->skip_txfm[0];
      pc_tree->horizontal[0].skip = x->skip;
      encode_b_rt(cpi, td, tile_info, tp, mi_row, mi_col, output_enabled,
                  subsize, &pc_tree->horizontal[0]);

      if (mi_row + hbs < cm->mi_rows && bsize > BLOCK_8X8) {
        pc_tree->horizontal[1].pred_pixel_ready = 1;
        nonrd_pick_sb_modes(cpi, tile_data, x, mi_row + hbs, mi_col, dummy_cost,
                            subsize, &pc_tree->horizontal[1]);
        pc_tree->horizontal[1].mic = *xd->mi[0];
        pc_tree->horizontal[1].mbmi_ext = *x->mbmi_ext;
        pc_tree->horizontal[1].skip_txfm[0] = x->skip_txfm[0];
        pc_tree->horizontal[1].skip = x->skip;
        encode_b_rt(cpi, td, tile_info, tp, mi_row + hbs, mi_col,
                    output_enabled, subsize, &pc_tree->horizontal[1]);
      }
      break;
    default:
      assert(partition == PARTITION_SPLIT);
      subsize = get_subsize(bsize, PARTITION_SPLIT);
      if (bsize == BLOCK_8X8) {
        nonrd_pick_sb_modes(cpi, tile_data, x, mi_row, mi_col, dummy_cost,
                            subsize, pc_tree->leaf_split[0]);
        encode_b_rt(cpi, td, tile_info, tp, mi_row, mi_col, output_enabled,
                    subsize, pc_tree->leaf_split[0]);
      } else {
        nonrd_use_partition(cpi, td, tile_data, mi, tp, mi_row, mi_col, subsize,
                            output_enabled, dummy_cost, pc_tree->split[0]);
        nonrd_use_partition(cpi, td, tile_data, mi + hbs, tp, mi_row,
                            mi_col + hbs, subsize, output_enabled, dummy_cost,
                            pc_tree->split[1]);
        nonrd_use_partition(cpi, td, tile_data, mi + hbs * mis, tp,
                            mi_row + hbs, mi_col, subsize, output_enabled,
                            dummy_cost, pc_tree->split[2]);
        nonrd_use_partition(cpi, td, tile_data, mi + hbs * mis + hbs, tp,
                            mi_row + hbs, mi_col + hbs, subsize, output_enabled,
                            dummy_cost, pc_tree->split[3]);
      }
      break;
  }

  if (partition != PARTITION_SPLIT || bsize == BLOCK_8X8)
    update_partition_context(xd, mi_row, mi_col, subsize, bsize);
}

static void encode_nonrd_sb_row(VP9_COMP *cpi, ThreadData *td,
                                TileDataEnc *tile_data, int mi_row,
                                TOKENEXTRA **tp) {
  SPEED_FEATURES *const sf = &cpi->sf;
  VP9_COMMON *const cm = &cpi->common;
  TileInfo *const tile_info = &tile_data->tile_info;
  MACROBLOCK *const x = &td->mb;
  MACROBLOCKD *const xd = &x->e_mbd;
  const int mi_col_start = tile_info->mi_col_start;
  const int mi_col_end = tile_info->mi_col_end;
  int mi_col;
  const int sb_row = mi_row >> MI_BLOCK_SIZE_LOG2;
  const int num_sb_cols =
      get_num_cols(tile_data->tile_info, MI_BLOCK_SIZE_LOG2);
  int sb_col_in_tile;

  // Initialize the left context for the new SB row
  memset(&xd->left_context, 0, sizeof(xd->left_context));
  memset(xd->left_seg_context, 0, sizeof(xd->left_seg_context));

  // Code each SB in the row
  for (mi_col = mi_col_start, sb_col_in_tile = 0; mi_col < mi_col_end;
       mi_col += MI_BLOCK_SIZE, ++sb_col_in_tile) {
    const struct segmentation *const seg = &cm->seg;
    RD_COST dummy_rdc;
    const int idx_str = cm->mi_stride * mi_row + mi_col;
    MODE_INFO **mi = cm->mi_grid_visible + idx_str;
    PARTITION_SEARCH_TYPE partition_search_type = sf->partition_search_type;
    BLOCK_SIZE bsize = BLOCK_64X64;
    int seg_skip = 0;
    int i;

    (*(cpi->row_mt_sync_read_ptr))(&tile_data->row_mt_sync, sb_row,
                                   sb_col_in_tile);

    if (cpi->use_skin_detection) {
      vp9_compute_skin_sb(cpi, BLOCK_16X16, mi_row, mi_col);
    }

    x->source_variance = UINT_MAX;
    for (i = 0; i < MAX_REF_FRAMES; ++i) {
      x->pred_mv[i].row = INT16_MAX;
      x->pred_mv[i].col = INT16_MAX;
    }
    vp9_rd_cost_init(&dummy_rdc);
    x->color_sensitivity[0] = 0;
    x->color_sensitivity[1] = 0;
    x->sb_is_skin = 0;
    x->skip_low_source_sad = 0;
    x->lowvar_highsumdiff = 0;
    x->content_state_sb = 0;
    x->zero_temp_sad_source = 0;
    x->sb_use_mv_part = 0;
    x->sb_mvcol_part = 0;
    x->sb_mvrow_part = 0;
    x->sb_pickmode_part = 0;
    x->arf_frame_usage = 0;
    x->lastgolden_frame_usage = 0;

    if (seg->enabled) {
      const uint8_t *const map =
          seg->update_map ? cpi->segmentation_map : cm->last_frame_seg_map;
      int segment_id = get_segment_id(cm, map, BLOCK_64X64, mi_row, mi_col);
      seg_skip = segfeature_active(seg, segment_id, SEG_LVL_SKIP);
      if (seg_skip) {
        partition_search_type = FIXED_PARTITION;
      }
    }

    if (cpi->compute_source_sad_onepass && cpi->sf.use_source_sad) {
      int shift = cpi->Source->y_stride * (mi_row << 3) + (mi_col << 3);
      int sb_offset2 = ((cm->mi_cols + 7) >> 3) * (mi_row >> 3) + (mi_col >> 3);
      int64_t source_sad = avg_source_sad(cpi, x, shift, sb_offset2);
      if (sf->adapt_partition_source_sad &&
          (cpi->oxcf.rc_mode == VPX_VBR && !cpi->rc.is_src_frame_alt_ref &&
           source_sad > sf->adapt_partition_thresh &&
           (cpi->refresh_golden_frame || cpi->refresh_alt_ref_frame)))
        partition_search_type = REFERENCE_PARTITION;
    }

    // Set the partition type of the 64X64 block
    switch (partition_search_type) {
      case VAR_BASED_PARTITION:
        // TODO(jingning, marpan): The mode decision and encoding process
        // support both intra and inter sub8x8 block coding for RTC mode.
        // Tune the thresholds accordingly to use sub8x8 block coding for
        // coding performance improvement.
        choose_partitioning(cpi, tile_info, x, mi_row, mi_col);
        nonrd_use_partition(cpi, td, tile_data, mi, tp, mi_row, mi_col,
                            BLOCK_64X64, 1, &dummy_rdc, td->pc_root);
        break;
      case SOURCE_VAR_BASED_PARTITION:
        set_source_var_based_partition(cpi, tile_info, x, mi, mi_row, mi_col);
        nonrd_use_partition(cpi, td, tile_data, mi, tp, mi_row, mi_col,
                            BLOCK_64X64, 1, &dummy_rdc, td->pc_root);
        break;
      case FIXED_PARTITION:
        if (!seg_skip) bsize = sf->always_this_block_size;
        set_fixed_partitioning(cpi, tile_info, mi, mi_row, mi_col, bsize);
        nonrd_use_partition(cpi, td, tile_data, mi, tp, mi_row, mi_col,
                            BLOCK_64X64, 1, &dummy_rdc, td->pc_root);
        break;
      default:
        assert(partition_search_type == REFERENCE_PARTITION);
        x->sb_pickmode_part = 1;
        set_offsets(cpi, tile_info, x, mi_row, mi_col, BLOCK_64X64);
        // Use nonrd_pick_partition on scene-cut for VBR mode.
        // nonrd_pick_partition does not support 4x4 partition, so avoid it
        // on key frame for now.
        if ((cpi->oxcf.rc_mode == VPX_VBR && cpi->rc.high_source_sad &&
             cpi->oxcf.speed < 6 && !frame_is_intra_only(cm) &&
             (cpi->refresh_golden_frame || cpi->refresh_alt_ref_frame))) {
          // Use lower max_partition_size for low resoultions.
          if (cm->width <= 352 && cm->height <= 288)
            x->max_partition_size = BLOCK_32X32;
          else
            x->max_partition_size = BLOCK_64X64;
          x->min_partition_size = BLOCK_8X8;
          nonrd_pick_partition(cpi, td, tile_data, tp, mi_row, mi_col,
                               BLOCK_64X64, &dummy_rdc, 1, INT64_MAX,
                               td->pc_root);
        } else {
          choose_partitioning(cpi, tile_info, x, mi_row, mi_col);
          // TODO(marpan): Seems like nonrd_select_partition does not support
          // 4x4 partition. Since 4x4 is used on key frame, use this switch
          // for now.
          if (frame_is_intra_only(cm))
            nonrd_use_partition(cpi, td, tile_data, mi, tp, mi_row, mi_col,
                                BLOCK_64X64, 1, &dummy_rdc, td->pc_root);
          else
            nonrd_select_partition(cpi, td, tile_data, mi, tp, mi_row, mi_col,
                                   BLOCK_64X64, 1, &dummy_rdc, td->pc_root);
        }

        break;
    }

    // Update ref_frame usage for inter frame if this group is ARF group.
    if (!cpi->rc.is_src_frame_alt_ref && !cpi->refresh_golden_frame &&
        !cpi->refresh_alt_ref_frame && cpi->rc.alt_ref_gf_group &&
        cpi->sf.use_altref_onepass) {
      int sboffset = ((cm->mi_cols + 7) >> 3) * (mi_row >> 3) + (mi_col >> 3);
      if (cpi->count_arf_frame_usage != NULL)
        cpi->count_arf_frame_usage[sboffset] = x->arf_frame_usage;
      if (cpi->count_lastgolden_frame_usage != NULL)
        cpi->count_lastgolden_frame_usage[sboffset] = x->lastgolden_frame_usage;
    }

    (*(cpi->row_mt_sync_write_ptr))(&tile_data->row_mt_sync, sb_row,
                                    sb_col_in_tile, num_sb_cols);
  }
}
// end RTC play code

static INLINE uint32_t variance(const diff *const d) {
  return d->sse - (uint32_t)(((int64_t)d->sum * d->sum) >> 8);
}

#if CONFIG_VP9_HIGHBITDEPTH
static INLINE uint32_t variance_highbd(diff *const d) {
  const int64_t var = (int64_t)d->sse - (((int64_t)d->sum * d->sum) >> 8);
  return (var >= 0) ? (uint32_t)var : 0;
}
#endif  // CONFIG_VP9_HIGHBITDEPTH

static int set_var_thresh_from_histogram(VP9_COMP *cpi) {
  const SPEED_FEATURES *const sf = &cpi->sf;
  const VP9_COMMON *const cm = &cpi->common;

  const uint8_t *src = cpi->Source->y_buffer;
  const uint8_t *last_src = cpi->Last_Source->y_buffer;
  const int src_stride = cpi->Source->y_stride;
  const int last_stride = cpi->Last_Source->y_stride;

  // Pick cutoff threshold
  const int cutoff = (VPXMIN(cm->width, cm->height) >= 720)
                         ? (cm->MBs * VAR_HIST_LARGE_CUT_OFF / 100)
                         : (cm->MBs * VAR_HIST_SMALL_CUT_OFF / 100);
  DECLARE_ALIGNED(16, int, hist[VAR_HIST_BINS]);
  diff *var16 = cpi->source_diff_var;

  int sum = 0;
  int i, j;

  memset(hist, 0, VAR_HIST_BINS * sizeof(hist[0]));

  for (i = 0; i < cm->mb_rows; i++) {
    for (j = 0; j < cm->mb_cols; j++) {
#if CONFIG_VP9_HIGHBITDEPTH
      if (cm->use_highbitdepth) {
        switch (cm->bit_depth) {
          case VPX_BITS_8:
            vpx_highbd_8_get16x16var(src, src_stride, last_src, last_stride,
                                     &var16->sse, &var16->sum);
            var16->var = variance(var16);
            break;
          case VPX_BITS_10:
            vpx_highbd_10_get16x16var(src, src_stride, last_src, last_stride,
                                      &var16->sse, &var16->sum);
            var16->var = variance_highbd(var16);
            break;
          default:
            assert(cm->bit_depth == VPX_BITS_12);
            vpx_highbd_12_get16x16var(src, src_stride, last_src, last_stride,
                                      &var16->sse, &var16->sum);
            var16->var = variance_highbd(var16);
            break;
        }
      } else {
        vpx_get16x16var(src, src_stride, last_src, last_stride, &var16->sse,
                        &var16->sum);
        var16->var = variance(var16);
      }
#else
      vpx_get16x16var(src, src_stride, last_src, last_stride, &var16->sse,
                      &var16->sum);
      var16->var = variance(var16);
#endif  // CONFIG_VP9_HIGHBITDEPTH

      if (var16->var >= VAR_HIST_MAX_BG_VAR)
        hist[VAR_HIST_BINS - 1]++;
      else
        hist[var16->var / VAR_HIST_FACTOR]++;

      src += 16;
      last_src += 16;
      var16++;
    }

    src = src - cm->mb_cols * 16 + 16 * src_stride;
    last_src = last_src - cm->mb_cols * 16 + 16 * last_stride;
  }

  cpi->source_var_thresh = 0;

  if (hist[VAR_HIST_BINS - 1] < cutoff) {
    for (i = 0; i < VAR_HIST_BINS - 1; i++) {
      sum += hist[i];

      if (sum > cutoff) {
        cpi->source_var_thresh = (i + 1) * VAR_HIST_FACTOR;
        return 0;
      }
    }
  }

  return sf->search_type_check_frequency;
}

static void source_var_based_partition_search_method(VP9_COMP *cpi) {
  VP9_COMMON *const cm = &cpi->common;
  SPEED_FEATURES *const sf = &cpi->sf;

  if (cm->frame_type == KEY_FRAME) {
    // For key frame, use SEARCH_PARTITION.
    sf->partition_search_type = SEARCH_PARTITION;
  } else if (cm->intra_only) {
    sf->partition_search_type = FIXED_PARTITION;
  } else {
    if (cm->last_width != cm->width || cm->last_height != cm->height) {
      if (cpi->source_diff_var) vpx_free(cpi->source_diff_var);

      CHECK_MEM_ERROR(cm, cpi->source_diff_var,
                      vpx_calloc(cm->MBs, sizeof(diff)));
    }

    if (!cpi->frames_till_next_var_check)
      cpi->frames_till_next_var_check = set_var_thresh_from_histogram(cpi);

    if (cpi->frames_till_next_var_check > 0) {
      sf->partition_search_type = FIXED_PARTITION;
      cpi->frames_till_next_var_check--;
    }
  }
}

static int get_skip_encode_frame(const VP9_COMMON *cm, ThreadData *const td) {
  unsigned int intra_count = 0, inter_count = 0;
  int j;

  for (j = 0; j < INTRA_INTER_CONTEXTS; ++j) {
    intra_count += td->counts->intra_inter[j][0];
    inter_count += td->counts->intra_inter[j][1];
  }

  return (intra_count << 2) < inter_count && cm->frame_type != KEY_FRAME &&
         cm->show_frame;
}

void vp9_init_tile_data(VP9_COMP *cpi) {
  VP9_COMMON *const cm = &cpi->common;
  const int tile_cols = 1 << cm->log2_tile_cols;
  const int tile_rows = 1 << cm->log2_tile_rows;
  int tile_col, tile_row;
  TOKENEXTRA *pre_tok = cpi->tile_tok[0][0];
  TOKENLIST *tplist = cpi->tplist[0][0];
  int tile_tok = 0;
  int tplist_count = 0;

  if (cpi->tile_data == NULL || cpi->allocated_tiles < tile_cols * tile_rows) {
    if (cpi->tile_data != NULL) vpx_free(cpi->tile_data);
    CHECK_MEM_ERROR(
        cm, cpi->tile_data,
        vpx_malloc(tile_cols * tile_rows * sizeof(*cpi->tile_data)));
    cpi->allocated_tiles = tile_cols * tile_rows;

    for (tile_row = 0; tile_row < tile_rows; ++tile_row)
      for (tile_col = 0; tile_col < tile_cols; ++tile_col) {
        TileDataEnc *tile_data =
            &cpi->tile_data[tile_row * tile_cols + tile_col];
        int i, j;
        for (i = 0; i < BLOCK_SIZES; ++i) {
          for (j = 0; j < MAX_MODES; ++j) {
            tile_data->thresh_freq_fact[i][j] = RD_THRESH_INIT_FACT;
#if CONFIG_CONSISTENT_RECODE
            tile_data->thresh_freq_fact_prev[i][j] = RD_THRESH_INIT_FACT;
#endif
            tile_data->mode_map[i][j] = j;
          }
        }
#if CONFIG_MULTITHREAD
        tile_data->row_base_thresh_freq_fact = NULL;
#endif
      }
  }

  for (tile_row = 0; tile_row < tile_rows; ++tile_row) {
    for (tile_col = 0; tile_col < tile_cols; ++tile_col) {
      TileDataEnc *this_tile = &cpi->tile_data[tile_row * tile_cols + tile_col];
      TileInfo *tile_info = &this_tile->tile_info;
      if (cpi->sf.adaptive_rd_thresh_row_mt &&
          this_tile->row_base_thresh_freq_fact == NULL)
        vp9_row_mt_alloc_rd_thresh(cpi, this_tile);
      vp9_tile_init(tile_info, cm, tile_row, tile_col);

      cpi->tile_tok[tile_row][tile_col] = pre_tok + tile_tok;
      pre_tok = cpi->tile_tok[tile_row][tile_col];
      tile_tok = allocated_tokens(*tile_info);

      cpi->tplist[tile_row][tile_col] = tplist + tplist_count;
      tplist = cpi->tplist[tile_row][tile_col];
      tplist_count = get_num_vert_units(*tile_info, MI_BLOCK_SIZE_LOG2);
    }
  }
}

void vp9_encode_sb_row(VP9_COMP *cpi, ThreadData *td, int tile_row,
                       int tile_col, int mi_row) {
  VP9_COMMON *const cm = &cpi->common;
  const int tile_cols = 1 << cm->log2_tile_cols;
  TileDataEnc *this_tile = &cpi->tile_data[tile_row * tile_cols + tile_col];
  const TileInfo *const tile_info = &this_tile->tile_info;
  TOKENEXTRA *tok = NULL;
  int tile_sb_row;
  int tile_mb_cols = (tile_info->mi_col_end - tile_info->mi_col_start + 1) >> 1;

  tile_sb_row = mi_cols_aligned_to_sb(mi_row - tile_info->mi_row_start) >>
                MI_BLOCK_SIZE_LOG2;
  get_start_tok(cpi, tile_row, tile_col, mi_row, &tok);
  cpi->tplist[tile_row][tile_col][tile_sb_row].start = tok;

  if (cpi->sf.use_nonrd_pick_mode)
    encode_nonrd_sb_row(cpi, td, this_tile, mi_row, &tok);
  else
    encode_rd_sb_row(cpi, td, this_tile, mi_row, &tok);

  cpi->tplist[tile_row][tile_col][tile_sb_row].stop = tok;
  cpi->tplist[tile_row][tile_col][tile_sb_row].count =
      (unsigned int)(cpi->tplist[tile_row][tile_col][tile_sb_row].stop -
                     cpi->tplist[tile_row][tile_col][tile_sb_row].start);
  assert(tok - cpi->tplist[tile_row][tile_col][tile_sb_row].start <=
         get_token_alloc(MI_BLOCK_SIZE >> 1, tile_mb_cols));

  (void)tile_mb_cols;
}

void vp9_encode_tile(VP9_COMP *cpi, ThreadData *td, int tile_row,
                     int tile_col) {
  VP9_COMMON *const cm = &cpi->common;
  const int tile_cols = 1 << cm->log2_tile_cols;
  TileDataEnc *this_tile = &cpi->tile_data[tile_row * tile_cols + tile_col];
  const TileInfo *const tile_info = &this_tile->tile_info;
  const int mi_row_start = tile_info->mi_row_start;
  const int mi_row_end = tile_info->mi_row_end;
  int mi_row;

  for (mi_row = mi_row_start; mi_row < mi_row_end; mi_row += MI_BLOCK_SIZE)
    vp9_encode_sb_row(cpi, td, tile_row, tile_col, mi_row);
}

static void encode_tiles(VP9_COMP *cpi) {
  VP9_COMMON *const cm = &cpi->common;
  const int tile_cols = 1 << cm->log2_tile_cols;
  const int tile_rows = 1 << cm->log2_tile_rows;
  int tile_col, tile_row;

  vp9_init_tile_data(cpi);

  for (tile_row = 0; tile_row < tile_rows; ++tile_row)
    for (tile_col = 0; tile_col < tile_cols; ++tile_col)
      vp9_encode_tile(cpi, &cpi->td, tile_row, tile_col);
}

#if CONFIG_FP_MB_STATS
static int input_fpmb_stats(FIRSTPASS_MB_STATS *firstpass_mb_stats,
                            VP9_COMMON *cm, uint8_t **this_frame_mb_stats) {
  uint8_t *mb_stats_in = firstpass_mb_stats->mb_stats_start +
                         cm->current_video_frame * cm->MBs * sizeof(uint8_t);

  if (mb_stats_in > firstpass_mb_stats->mb_stats_end) return EOF;

  *this_frame_mb_stats = mb_stats_in;

  return 1;
}
#endif

static void encode_frame_internal(VP9_COMP *cpi) {
  SPEED_FEATURES *const sf = &cpi->sf;
  ThreadData *const td = &cpi->td;
  MACROBLOCK *const x = &td->mb;
  VP9_COMMON *const cm = &cpi->common;
  MACROBLOCKD *const xd = &x->e_mbd;
  const int gf_group_index = cpi->twopass.gf_group.index;

  xd->mi = cm->mi_grid_visible;
  xd->mi[0] = cm->mi;

  vp9_zero(*td->counts);
  vp9_zero(cpi->td.rd_counts);

  xd->lossless = cm->base_qindex == 0 && cm->y_dc_delta_q == 0 &&
                 cm->uv_dc_delta_q == 0 && cm->uv_ac_delta_q == 0;

#if CONFIG_VP9_HIGHBITDEPTH
  if (cm->use_highbitdepth)
    x->fwd_txfm4x4 = xd->lossless ? vp9_highbd_fwht4x4 : vpx_highbd_fdct4x4;
  else
    x->fwd_txfm4x4 = xd->lossless ? vp9_fwht4x4 : vpx_fdct4x4;
  x->highbd_inv_txfm_add =
      xd->lossless ? vp9_highbd_iwht4x4_add : vp9_highbd_idct4x4_add;
#else
  x->fwd_txfm4x4 = xd->lossless ? vp9_fwht4x4 : vpx_fdct4x4;
#endif  // CONFIG_VP9_HIGHBITDEPTH
  x->inv_txfm_add = xd->lossless ? vp9_iwht4x4_add : vp9_idct4x4_add;
#if CONFIG_CONSISTENT_RECODE
  x->optimize = sf->optimize_coefficients == 1 && cpi->oxcf.pass != 1;
#endif
  if (xd->lossless) x->optimize = 0;
  x->sharpness = cpi->oxcf.sharpness;
  x->adjust_rdmult_by_segment = (cpi->oxcf.aq_mode == VARIANCE_AQ);

  cm->tx_mode = select_tx_mode(cpi, xd);

  vp9_frame_init_quantizer(cpi);

  vp9_initialize_rd_consts(cpi);
  vp9_initialize_me_consts(cpi, x, cm->base_qindex);
  init_encode_frame_mb_context(cpi);
  cm->use_prev_frame_mvs =
      !cm->error_resilient_mode && cm->width == cm->last_width &&
      cm->height == cm->last_height && !cm->intra_only && cm->last_show_frame;
  // Special case: set prev_mi to NULL when the previous mode info
  // context cannot be used.
  cm->prev_mi =
      cm->use_prev_frame_mvs ? cm->prev_mip + cm->mi_stride + 1 : NULL;

  x->quant_fp = cpi->sf.use_quant_fp;
  vp9_zero(x->skip_txfm);
  if (sf->use_nonrd_pick_mode) {
    // Initialize internal buffer pointers for rtc coding, where non-RD
    // mode decision is used and hence no buffer pointer swap needed.
    int i;
    struct macroblock_plane *const p = x->plane;
    struct macroblockd_plane *const pd = xd->plane;
    PICK_MODE_CONTEXT *ctx = &cpi->td.pc_root->none;

    for (i = 0; i < MAX_MB_PLANE; ++i) {
      p[i].coeff = ctx->coeff_pbuf[i][0];
      p[i].qcoeff = ctx->qcoeff_pbuf[i][0];
      pd[i].dqcoeff = ctx->dqcoeff_pbuf[i][0];
      p[i].eobs = ctx->eobs_pbuf[i][0];
    }
    vp9_zero(x->zcoeff_blk);

    if (cm->frame_type != KEY_FRAME && cpi->rc.frames_since_golden == 0 &&
        !(cpi->oxcf.lag_in_frames > 0 && cpi->oxcf.rc_mode == VPX_VBR) &&
        !cpi->use_svc)
      cpi->ref_frame_flags &= (~VP9_GOLD_FLAG);

    if (sf->partition_search_type == SOURCE_VAR_BASED_PARTITION)
      source_var_based_partition_search_method(cpi);
  } else if (gf_group_index && gf_group_index < MAX_LAG_BUFFERS &&
             cpi->sf.enable_tpl_model) {
    TplDepFrame *tpl_frame = &cpi->tpl_stats[cpi->twopass.gf_group.index];
    TplDepStats *tpl_stats = tpl_frame->tpl_stats_ptr;

    int tpl_stride = tpl_frame->stride;
    int64_t intra_cost_base = 0;
    int64_t mc_dep_cost_base = 0;
    int row, col;

    for (row = 0; row < cm->mi_rows; ++row) {
      for (col = 0; col < cm->mi_cols; ++col) {
        TplDepStats *this_stats = &tpl_stats[row * tpl_stride + col];
        intra_cost_base += this_stats->intra_cost;
        mc_dep_cost_base += this_stats->mc_dep_cost;
      }
    }

    vpx_clear_system_state();

    if (tpl_frame->is_valid)
      cpi->rd.r0 = (double)intra_cost_base / mc_dep_cost_base;
  }

  {
    struct vpx_usec_timer emr_timer;
    vpx_usec_timer_start(&emr_timer);

#if CONFIG_FP_MB_STATS
    if (cpi->use_fp_mb_stats) {
      input_fpmb_stats(&cpi->twopass.firstpass_mb_stats, cm,
                       &cpi->twopass.this_frame_mb_stats);
    }
#endif

    if (!cpi->row_mt) {
      cpi->row_mt_sync_read_ptr = vp9_row_mt_sync_read_dummy;
      cpi->row_mt_sync_write_ptr = vp9_row_mt_sync_write_dummy;
      // If allowed, encoding tiles in parallel with one thread handling one
      // tile when row based multi-threading is disabled.
      if (VPXMIN(cpi->oxcf.max_threads, 1 << cm->log2_tile_cols) > 1)
        vp9_encode_tiles_mt(cpi);
      else
        encode_tiles(cpi);
    } else {
      cpi->row_mt_sync_read_ptr = vp9_row_mt_sync_read;
      cpi->row_mt_sync_write_ptr = vp9_row_mt_sync_write;
      vp9_encode_tiles_row_mt(cpi);
    }

    vpx_usec_timer_mark(&emr_timer);
    cpi->time_encode_sb_row += vpx_usec_timer_elapsed(&emr_timer);
  }

  sf->skip_encode_frame =
      sf->skip_encode_sb ? get_skip_encode_frame(cm, td) : 0;

#if 0
  // Keep record of the total distortion this time around for future use
  cpi->last_frame_distortion = cpi->frame_distortion;
#endif
}

static INTERP_FILTER get_interp_filter(
    const int64_t threshes[SWITCHABLE_FILTER_CONTEXTS], int is_alt_ref) {
  if (!is_alt_ref && threshes[EIGHTTAP_SMOOTH] > threshes[EIGHTTAP] &&
      threshes[EIGHTTAP_SMOOTH] > threshes[EIGHTTAP_SHARP] &&
      threshes[EIGHTTAP_SMOOTH] > threshes[SWITCHABLE - 1]) {
    return EIGHTTAP_SMOOTH;
  } else if (threshes[EIGHTTAP_SHARP] > threshes[EIGHTTAP] &&
             threshes[EIGHTTAP_SHARP] > threshes[SWITCHABLE - 1]) {
    return EIGHTTAP_SHARP;
  } else if (threshes[EIGHTTAP] > threshes[SWITCHABLE - 1]) {
    return EIGHTTAP;
  } else {
    return SWITCHABLE;
  }
}

static int compute_frame_aq_offset(struct VP9_COMP *cpi) {
  VP9_COMMON *const cm = &cpi->common;
  MODE_INFO **mi_8x8_ptr = cm->mi_grid_visible;
  struct segmentation *const seg = &cm->seg;

  int mi_row, mi_col;
  int sum_delta = 0;
  int map_index = 0;
  int qdelta_index;
  int segment_id;

  for (mi_row = 0; mi_row < cm->mi_rows; mi_row++) {
    MODE_INFO **mi_8x8 = mi_8x8_ptr;
    for (mi_col = 0; mi_col < cm->mi_cols; mi_col++, mi_8x8++) {
      segment_id = mi_8x8[0]->segment_id;
      qdelta_index = get_segdata(seg, segment_id, SEG_LVL_ALT_Q);
      sum_delta += qdelta_index;
      map_index++;
    }
    mi_8x8_ptr += cm->mi_stride;
  }

  return sum_delta / (cm->mi_rows * cm->mi_cols);
}

#if CONFIG_CONSISTENT_RECODE
static void restore_encode_params(VP9_COMP *cpi) {
  VP9_COMMON *const cm = &cpi->common;
  const int tile_cols = 1 << cm->log2_tile_cols;
  const int tile_rows = 1 << cm->log2_tile_rows;
  int tile_col, tile_row;
  int i, j;
  RD_OPT *rd_opt = &cpi->rd;
  for (i = 0; i < MAX_REF_FRAMES; i++) {
    for (j = 0; j < REFERENCE_MODES; j++)
      rd_opt->prediction_type_threshes[i][j] =
          rd_opt->prediction_type_threshes_prev[i][j];

    for (j = 0; j < SWITCHABLE_FILTER_CONTEXTS; j++)
      rd_opt->filter_threshes[i][j] = rd_opt->filter_threshes_prev[i][j];
  }

  if (cpi->tile_data != NULL) {
    for (tile_row = 0; tile_row < tile_rows; ++tile_row)
      for (tile_col = 0; tile_col < tile_cols; ++tile_col) {
        TileDataEnc *tile_data =
            &cpi->tile_data[tile_row * tile_cols + tile_col];
        for (i = 0; i < BLOCK_SIZES; ++i) {
          for (j = 0; j < MAX_MODES; ++j) {
            tile_data->thresh_freq_fact[i][j] =
                tile_data->thresh_freq_fact_prev[i][j];
          }
        }
      }
  }

  cm->interp_filter = cpi->sf.default_interp_filter;
}
#endif

void vp9_encode_frame(VP9_COMP *cpi) {
  VP9_COMMON *const cm = &cpi->common;

#if CONFIG_CONSISTENT_RECODE
  restore_encode_params(cpi);
#endif

  // In the longer term the encoder should be generalized to match the
  // decoder such that we allow compound where one of the 3 buffers has a
  // different sign bias and that buffer is then the fixed ref. However, this
  // requires further work in the rd loop. For now the only supported encoder
  // side behavior is where the ALT ref buffer has opposite sign bias to
  // the other two.
  if (!frame_is_intra_only(cm)) {
    if ((cm->ref_frame_sign_bias[ALTREF_FRAME] ==
         cm->ref_frame_sign_bias[GOLDEN_FRAME]) ||
        (cm->ref_frame_sign_bias[ALTREF_FRAME] ==
         cm->ref_frame_sign_bias[LAST_FRAME])) {
      cpi->allow_comp_inter_inter = 0;
    } else {
      cpi->allow_comp_inter_inter = 1;
      cm->comp_fixed_ref = ALTREF_FRAME;
      cm->comp_var_ref[0] = LAST_FRAME;
      cm->comp_var_ref[1] = GOLDEN_FRAME;
    }
  }

  if (cpi->sf.frame_parameter_update) {
    int i;
    RD_OPT *const rd_opt = &cpi->rd;
    FRAME_COUNTS *counts = cpi->td.counts;
    RD_COUNTS *const rdc = &cpi->td.rd_counts;

    // This code does a single RD pass over the whole frame assuming
    // either compound, single or hybrid prediction as per whatever has
    // worked best for that type of frame in the past.
    // It also predicts whether another coding mode would have worked
    // better than this coding mode. If that is the case, it remembers
    // that for subsequent frames.
    // It also does the same analysis for transform size selection.
    const MV_REFERENCE_FRAME frame_type = get_frame_type(cpi);
    int64_t *const mode_thrs = rd_opt->prediction_type_threshes[frame_type];
    int64_t *const filter_thrs = rd_opt->filter_threshes[frame_type];
    const int is_alt_ref = frame_type == ALTREF_FRAME;

    /* prediction (compound, single or hybrid) mode selection */
    if (is_alt_ref || !cpi->allow_comp_inter_inter)
      cm->reference_mode = SINGLE_REFERENCE;
    else if (mode_thrs[COMPOUND_REFERENCE] > mode_thrs[SINGLE_REFERENCE] &&
             mode_thrs[COMPOUND_REFERENCE] > mode_thrs[REFERENCE_MODE_SELECT] &&
             check_dual_ref_flags(cpi) && cpi->static_mb_pct == 100)
      cm->reference_mode = COMPOUND_REFERENCE;
    else if (mode_thrs[SINGLE_REFERENCE] > mode_thrs[REFERENCE_MODE_SELECT])
      cm->reference_mode = SINGLE_REFERENCE;
    else
      cm->reference_mode = REFERENCE_MODE_SELECT;

    if (cm->interp_filter == SWITCHABLE)
      cm->interp_filter = get_interp_filter(filter_thrs, is_alt_ref);

    encode_frame_internal(cpi);

    for (i = 0; i < REFERENCE_MODES; ++i)
      mode_thrs[i] = (mode_thrs[i] + rdc->comp_pred_diff[i] / cm->MBs) / 2;

    for (i = 0; i < SWITCHABLE_FILTER_CONTEXTS; ++i)
      filter_thrs[i] = (filter_thrs[i] + rdc->filter_diff[i] / cm->MBs) / 2;

    if (cm->reference_mode == REFERENCE_MODE_SELECT) {
      int single_count_zero = 0;
      int comp_count_zero = 0;

      for (i = 0; i < COMP_INTER_CONTEXTS; i++) {
        single_count_zero += counts->comp_inter[i][0];
        comp_count_zero += counts->comp_inter[i][1];
      }

      if (comp_count_zero == 0) {
        cm->reference_mode = SINGLE_REFERENCE;
        vp9_zero(counts->comp_inter);
      } else if (single_count_zero == 0) {
        cm->reference_mode = COMPOUND_REFERENCE;
        vp9_zero(counts->comp_inter);
      }
    }

    if (cm->tx_mode == TX_MODE_SELECT) {
      int count4x4 = 0;
      int count8x8_lp = 0, count8x8_8x8p = 0;
      int count16x16_16x16p = 0, count16x16_lp = 0;
      int count32x32 = 0;

      for (i = 0; i < TX_SIZE_CONTEXTS; ++i) {
        count4x4 += counts->tx.p32x32[i][TX_4X4];
        count4x4 += counts->tx.p16x16[i][TX_4X4];
        count4x4 += counts->tx.p8x8[i][TX_4X4];

        count8x8_lp += counts->tx.p32x32[i][TX_8X8];
        count8x8_lp += counts->tx.p16x16[i][TX_8X8];
        count8x8_8x8p += counts->tx.p8x8[i][TX_8X8];

        count16x16_16x16p += counts->tx.p16x16[i][TX_16X16];
        count16x16_lp += counts->tx.p32x32[i][TX_16X16];
        count32x32 += counts->tx.p32x32[i][TX_32X32];
      }
      if (count4x4 == 0 && count16x16_lp == 0 && count16x16_16x16p == 0 &&
          count32x32 == 0) {
        cm->tx_mode = ALLOW_8X8;
        reset_skip_tx_size(cm, TX_8X8);
      } else if (count8x8_8x8p == 0 && count16x16_16x16p == 0 &&
                 count8x8_lp == 0 && count16x16_lp == 0 && count32x32 == 0) {
        cm->tx_mode = ONLY_4X4;
        reset_skip_tx_size(cm, TX_4X4);
      } else if (count8x8_lp == 0 && count16x16_lp == 0 && count4x4 == 0) {
        cm->tx_mode = ALLOW_32X32;
      } else if (count32x32 == 0 && count8x8_lp == 0 && count4x4 == 0) {
        cm->tx_mode = ALLOW_16X16;
        reset_skip_tx_size(cm, TX_16X16);
      }
    }
  } else {
    FRAME_COUNTS *counts = cpi->td.counts;
    cm->reference_mode = SINGLE_REFERENCE;
    if (cpi->allow_comp_inter_inter && cpi->sf.use_compound_nonrd_pickmode &&
        cpi->rc.alt_ref_gf_group && !cpi->rc.is_src_frame_alt_ref &&
        cm->frame_type != KEY_FRAME)
      cm->reference_mode = REFERENCE_MODE_SELECT;

    encode_frame_internal(cpi);

    if (cm->reference_mode == REFERENCE_MODE_SELECT) {
      int single_count_zero = 0;
      int comp_count_zero = 0;
      int i;
      for (i = 0; i < COMP_INTER_CONTEXTS; i++) {
        single_count_zero += counts->comp_inter[i][0];
        comp_count_zero += counts->comp_inter[i][1];
      }
      if (comp_count_zero == 0) {
        cm->reference_mode = SINGLE_REFERENCE;
        vp9_zero(counts->comp_inter);
      } else if (single_count_zero == 0) {
        cm->reference_mode = COMPOUND_REFERENCE;
        vp9_zero(counts->comp_inter);
      }
    }
  }

  // If segmented AQ is enabled compute the average AQ weighting.
  if (cm->seg.enabled && (cpi->oxcf.aq_mode != NO_AQ) &&
      (cm->seg.update_map || cm->seg.update_data)) {
    cm->seg.aq_av_offset = compute_frame_aq_offset(cpi);
  }
}

static void sum_intra_stats(FRAME_COUNTS *counts, const MODE_INFO *mi) {
  const PREDICTION_MODE y_mode = mi->mode;
  const PREDICTION_MODE uv_mode = mi->uv_mode;
  const BLOCK_SIZE bsize = mi->sb_type;

  if (bsize < BLOCK_8X8) {
    int idx, idy;
    const int num_4x4_w = num_4x4_blocks_wide_lookup[bsize];
    const int num_4x4_h = num_4x4_blocks_high_lookup[bsize];
    for (idy = 0; idy < 2; idy += num_4x4_h)
      for (idx = 0; idx < 2; idx += num_4x4_w)
        ++counts->y_mode[0][mi->bmi[idy * 2 + idx].as_mode];
  } else {
    ++counts->y_mode[size_group_lookup[bsize]][y_mode];
  }

  ++counts->uv_mode[y_mode][uv_mode];
}

static void update_zeromv_cnt(VP9_COMP *const cpi, const MODE_INFO *const mi,
                              int mi_row, int mi_col, BLOCK_SIZE bsize) {
  const VP9_COMMON *const cm = &cpi->common;
  MV mv = mi->mv[0].as_mv;
  const int bw = num_8x8_blocks_wide_lookup[bsize];
  const int bh = num_8x8_blocks_high_lookup[bsize];
  const int xmis = VPXMIN(cm->mi_cols - mi_col, bw);
  const int ymis = VPXMIN(cm->mi_rows - mi_row, bh);
  const int block_index = mi_row * cm->mi_cols + mi_col;
  int x, y;
  for (y = 0; y < ymis; y++)
    for (x = 0; x < xmis; x++) {
      int map_offset = block_index + y * cm->mi_cols + x;
      if (mi->ref_frame[0] == LAST_FRAME && is_inter_block(mi) &&
          mi->segment_id <= CR_SEGMENT_ID_BOOST2) {
        if (abs(mv.row) < 8 && abs(mv.col) < 8) {
          if (cpi->consec_zero_mv[map_offset] < 255)
            cpi->consec_zero_mv[map_offset]++;
        } else {
          cpi->consec_zero_mv[map_offset] = 0;
        }
      }
    }
}

static void encode_superblock(VP9_COMP *cpi, ThreadData *td, TOKENEXTRA **t,
                              int output_enabled, int mi_row, int mi_col,
                              BLOCK_SIZE bsize, PICK_MODE_CONTEXT *ctx) {
  VP9_COMMON *const cm = &cpi->common;
  MACROBLOCK *const x = &td->mb;
  MACROBLOCKD *const xd = &x->e_mbd;
  MODE_INFO *mi = xd->mi[0];
  const int seg_skip =
      segfeature_active(&cm->seg, mi->segment_id, SEG_LVL_SKIP);
  x->skip_recode = !x->select_tx_size && mi->sb_type >= BLOCK_8X8 &&
                   cpi->oxcf.aq_mode != COMPLEXITY_AQ &&
                   cpi->oxcf.aq_mode != CYCLIC_REFRESH_AQ &&
                   cpi->sf.allow_skip_recode;

  if (!x->skip_recode && !cpi->sf.use_nonrd_pick_mode)
    memset(x->skip_txfm, 0, sizeof(x->skip_txfm));

  x->skip_optimize = ctx->is_coded;
  ctx->is_coded = 1;
  x->use_lp32x32fdct = cpi->sf.use_lp32x32fdct;
  x->skip_encode = (!output_enabled && cpi->sf.skip_encode_frame &&
                    x->q_index < QIDX_SKIP_THRESH);

  if (x->skip_encode) return;

  if (!is_inter_block(mi)) {
    int plane;
#if CONFIG_BETTER_HW_COMPATIBILITY && CONFIG_VP9_HIGHBITDEPTH
    if ((xd->cur_buf->flags & YV12_FLAG_HIGHBITDEPTH) &&
        (xd->above_mi == NULL || xd->left_mi == NULL) &&
        need_top_left[mi->uv_mode])
      assert(0);
#endif  // CONFIG_BETTER_HW_COMPATIBILITY && CONFIG_VP9_HIGHBITDEPTH
    mi->skip = 1;
    for (plane = 0; plane < MAX_MB_PLANE; ++plane)
      vp9_encode_intra_block_plane(x, VPXMAX(bsize, BLOCK_8X8), plane, 1);
    if (output_enabled) sum_intra_stats(td->counts, mi);
    vp9_tokenize_sb(cpi, td, t, !output_enabled, seg_skip,
                    VPXMAX(bsize, BLOCK_8X8));
  } else {
    int ref;
    const int is_compound = has_second_ref(mi);
    set_ref_ptrs(cm, xd, mi->ref_frame[0], mi->ref_frame[1]);
    for (ref = 0; ref < 1 + is_compound; ++ref) {
      YV12_BUFFER_CONFIG *cfg = get_ref_frame_buffer(cpi, mi->ref_frame[ref]);
      assert(cfg != NULL);
      vp9_setup_pre_planes(xd, ref, cfg, mi_row, mi_col,
                           &xd->block_refs[ref]->sf);
    }
    if (!(cpi->sf.reuse_inter_pred_sby && ctx->pred_pixel_ready) || seg_skip)
      vp9_build_inter_predictors_sby(xd, mi_row, mi_col,
                                     VPXMAX(bsize, BLOCK_8X8));

    vp9_build_inter_predictors_sbuv(xd, mi_row, mi_col,
                                    VPXMAX(bsize, BLOCK_8X8));

    vp9_encode_sb(x, VPXMAX(bsize, BLOCK_8X8));
    vp9_tokenize_sb(cpi, td, t, !output_enabled, seg_skip,
                    VPXMAX(bsize, BLOCK_8X8));
  }

  if (seg_skip) {
    assert(mi->skip);
  }

  if (output_enabled) {
    if (cm->tx_mode == TX_MODE_SELECT && mi->sb_type >= BLOCK_8X8 &&
        !(is_inter_block(mi) && mi->skip)) {
      ++get_tx_counts(max_txsize_lookup[bsize], get_tx_size_context(xd),
                      &td->counts->tx)[mi->tx_size];
    } else {
      // The new intra coding scheme requires no change of transform size
      if (is_inter_block(mi)) {
        mi->tx_size = VPXMIN(tx_mode_to_biggest_tx_size[cm->tx_mode],
                             max_txsize_lookup[bsize]);
      } else {
        mi->tx_size = (bsize >= BLOCK_8X8) ? mi->tx_size : TX_4X4;
      }
    }

    ++td->counts->tx.tx_totals[mi->tx_size];
    ++td->counts->tx.tx_totals[get_uv_tx_size(mi, &xd->plane[1])];
    if (cm->seg.enabled && cpi->oxcf.aq_mode == CYCLIC_REFRESH_AQ)
      vp9_cyclic_refresh_update_sb_postencode(cpi, mi, mi_row, mi_col, bsize);
    if (cpi->oxcf.pass == 0 && cpi->svc.temporal_layer_id == 0 &&
        (!cpi->use_svc ||
         (cpi->use_svc &&
          !cpi->svc.layer_context[cpi->svc.temporal_layer_id].is_key_frame &&
          cpi->svc.spatial_layer_id == cpi->svc.number_spatial_layers - 1)))
      update_zeromv_cnt(cpi, mi, mi_row, mi_col, bsize);
  }
}
