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
#include "./vpx_config.h"

#include "vpx_ports/vpx_timer.h"

#include "vp9/common/vp9_common.h"
#include "vp9/common/vp9_entropy.h"
#include "vp9/common/vp9_entropymode.h"
#include "vp9/common/vp9_idct.h"
#include "vp9/common/vp9_mvref_common.h"
#if CONFIG_PALETTE
#include "vp9/common/vp9_palette.h"
#endif
#if CONFIG_SR_MODE
#include "vp9/common/vp9_sr_txfm.h"
#endif  // CONFIG_SR_MODE
#include "vp9/common/vp9_pred_common.h"
#include "vp9/common/vp9_quant_common.h"
#include "vp9/common/vp9_reconintra.h"
#include "vp9/common/vp9_reconinter.h"
#include "vp9/common/vp9_seg_common.h"
#include "vp9/common/vp9_systemdependent.h"
#include "vp9/common/vp9_tile_common.h"

#include "vp9/encoder/vp9_aq_complexity.h"
#include "vp9/encoder/vp9_aq_cyclicrefresh.h"
#include "vp9/encoder/vp9_aq_variance.h"
#if CONFIG_SUPERTX
#include "vp9/encoder/vp9_cost.h"
#endif
#if CONFIG_GLOBAL_MOTION
#include "vp9/encoder/vp9_global_motion.h"
#endif  // CONFIG_GLOBAL_MOTION
#include "vp9/encoder/vp9_encodeframe.h"
#include "vp9/encoder/vp9_encodemb.h"
#include "vp9/encoder/vp9_encodemv.h"
#include "vp9/encoder/vp9_extend.h"
#include "vp9/encoder/vp9_rd.h"
#include "vp9/encoder/vp9_rdopt.h"
#include "vp9/encoder/vp9_segmentation.h"
#include "vp9/encoder/vp9_tokenize.h"

#define GF_ZEROMV_ZBIN_BOOST 0
#define LF_ZEROMV_ZBIN_BOOST 0
#define MV_ZBIN_BOOST        0
#define SPLIT_MV_ZBIN_BOOST  0
#define INTRA_ZBIN_BOOST     0

static void encode_superblock(VP9_COMP *cpi, TOKENEXTRA **t, int output_enabled,
                              int mi_row, int mi_col, BLOCK_SIZE bsize,
                              PICK_MODE_CONTEXT *ctx);

#if CONFIG_SUPERTX
static int check_intra_b(PICK_MODE_CONTEXT *ctx);

static int check_intra_sb(VP9_COMP *cpi, const TileInfo *const tile,
                          int mi_row, int mi_col, BLOCK_SIZE bsize,
                          PC_TREE *pc_tree);
static void predict_superblock(VP9_COMP *cpi,
#if CONFIG_WEDGE_PARTITION
                               int mi_row_ori, int mi_col_ori,
#endif  // CONFIG_WEDGE_PARTITION
                               int mi_row_pred, int mi_col_pred,
                               BLOCK_SIZE bsize_pred, int b_sub8x8, int block);
static int check_supertx_sb(BLOCK_SIZE bsize, TX_SIZE supertx_size,
                            PC_TREE *pc_tree);
static void predict_sb_complex(VP9_COMP *cpi, const TileInfo *const tile,
                               int mi_row, int mi_col,
                               int mi_row_ori, int mi_col_ori,
                               int output_enabled, BLOCK_SIZE bsize,
                               BLOCK_SIZE top_bsize,
                               uint8_t *dst_buf[3], int dst_stride[3],
                               PC_TREE *pc_tree);
static void update_state_sb_supertx(VP9_COMP *cpi, const TileInfo *const tile,
                                    int mi_row, int mi_col,
                                    BLOCK_SIZE bsize,
                                    int output_enabled, PC_TREE *pc_tree);
static void rd_supertx_sb(VP9_COMP *cpi, const TileInfo *const tile,
                          int mi_row, int mi_col, BLOCK_SIZE bsize,
                          int *tmp_rate, int64_t *tmp_dist,
#if CONFIG_EXT_TX
                          EXT_TX_TYPE *best_tx,
#endif
#if CONFIG_NEW_QUANT && QUANT_PROFILES > 1
                          int *dq_index,
#endif  // CONFIG_NEW_QUANT && QUANT_PROFILES > 1
                          PC_TREE *pc_tree);
#endif  // CONFIG_SUPERTX

// Motion vector component magnitude threshold for defining fast motion.
#define FAST_MOTION_MV_THRESH 24

// This is used as a reference when computing the source variance for the
//  purposes of activity masking.
// Eventually this should be replaced by custom no-reference routines,
//  which will be faster.
static const uint8_t VP9_VAR_OFFS[CODING_UNIT_SIZE] = {
    128, 128, 128, 128, 128, 128, 128, 128,
    128, 128, 128, 128, 128, 128, 128, 128,
    128, 128, 128, 128, 128, 128, 128, 128,
    128, 128, 128, 128, 128, 128, 128, 128,
    128, 128, 128, 128, 128, 128, 128, 128,
    128, 128, 128, 128, 128, 128, 128, 128,
    128, 128, 128, 128, 128, 128, 128, 128,
    128, 128, 128, 128, 128, 128, 128, 128,
#if CONFIG_EXT_CODING_UNIT_SIZE
    128, 128, 128, 128, 128, 128, 128, 128,
    128, 128, 128, 128, 128, 128, 128, 128,
    128, 128, 128, 128, 128, 128, 128, 128,
    128, 128, 128, 128, 128, 128, 128, 128,
    128, 128, 128, 128, 128, 128, 128, 128,
    128, 128, 128, 128, 128, 128, 128, 128,
    128, 128, 128, 128, 128, 128, 128, 128,
    128, 128, 128, 128, 128, 128, 128, 128,
#endif
};

#if CONFIG_VP9_HIGHBITDEPTH
static const uint16_t VP9_HIGH_VAR_OFFS_8[CODING_UNIT_SIZE] = {
    128, 128, 128, 128, 128, 128, 128, 128,
    128, 128, 128, 128, 128, 128, 128, 128,
    128, 128, 128, 128, 128, 128, 128, 128,
    128, 128, 128, 128, 128, 128, 128, 128,
    128, 128, 128, 128, 128, 128, 128, 128,
    128, 128, 128, 128, 128, 128, 128, 128,
    128, 128, 128, 128, 128, 128, 128, 128,
    128, 128, 128, 128, 128, 128, 128, 128,
#if CONFIG_EXT_CODING_UNIT_SIZE
    128, 128, 128, 128, 128, 128, 128, 128,
    128, 128, 128, 128, 128, 128, 128, 128,
    128, 128, 128, 128, 128, 128, 128, 128,
    128, 128, 128, 128, 128, 128, 128, 128,
    128, 128, 128, 128, 128, 128, 128, 128,
    128, 128, 128, 128, 128, 128, 128, 128,
    128, 128, 128, 128, 128, 128, 128, 128,
    128, 128, 128, 128, 128, 128, 128, 128
#endif
};

static const uint16_t VP9_HIGH_VAR_OFFS_10[CODING_UNIT_SIZE] = {
    128*4, 128*4, 128*4, 128*4, 128*4, 128*4, 128*4, 128*4,
    128*4, 128*4, 128*4, 128*4, 128*4, 128*4, 128*4, 128*4,
    128*4, 128*4, 128*4, 128*4, 128*4, 128*4, 128*4, 128*4,
    128*4, 128*4, 128*4, 128*4, 128*4, 128*4, 128*4, 128*4,
    128*4, 128*4, 128*4, 128*4, 128*4, 128*4, 128*4, 128*4,
    128*4, 128*4, 128*4, 128*4, 128*4, 128*4, 128*4, 128*4,
    128*4, 128*4, 128*4, 128*4, 128*4, 128*4, 128*4, 128*4,
    128*4, 128*4, 128*4, 128*4, 128*4, 128*4, 128*4, 128*4,
#if CONFIG_EXT_CODING_UNIT_SIZE
    128*4, 128*4, 128*4, 128*4, 128*4, 128*4, 128*4, 128*4,
    128*4, 128*4, 128*4, 128*4, 128*4, 128*4, 128*4, 128*4,
    128*4, 128*4, 128*4, 128*4, 128*4, 128*4, 128*4, 128*4,
    128*4, 128*4, 128*4, 128*4, 128*4, 128*4, 128*4, 128*4,
    128*4, 128*4, 128*4, 128*4, 128*4, 128*4, 128*4, 128*4,
    128*4, 128*4, 128*4, 128*4, 128*4, 128*4, 128*4, 128*4,
    128*4, 128*4, 128*4, 128*4, 128*4, 128*4, 128*4, 128*4,
    128*4, 128*4, 128*4, 128*4, 128*4, 128*4, 128*4, 128*4
#endif
};

static const uint16_t VP9_HIGH_VAR_OFFS_12[CODING_UNIT_SIZE] = {
    128*16, 128*16, 128*16, 128*16, 128*16, 128*16, 128*16, 128*16,
    128*16, 128*16, 128*16, 128*16, 128*16, 128*16, 128*16, 128*16,
    128*16, 128*16, 128*16, 128*16, 128*16, 128*16, 128*16, 128*16,
    128*16, 128*16, 128*16, 128*16, 128*16, 128*16, 128*16, 128*16,
    128*16, 128*16, 128*16, 128*16, 128*16, 128*16, 128*16, 128*16,
    128*16, 128*16, 128*16, 128*16, 128*16, 128*16, 128*16, 128*16,
    128*16, 128*16, 128*16, 128*16, 128*16, 128*16, 128*16, 128*16,
    128*16, 128*16, 128*16, 128*16, 128*16, 128*16, 128*16, 128*16,
#if CONFIG_EXT_CODING_UNIT_SIZE
    128*16, 128*16, 128*16, 128*16, 128*16, 128*16, 128*16, 128*16,
    128*16, 128*16, 128*16, 128*16, 128*16, 128*16, 128*16, 128*16,
    128*16, 128*16, 128*16, 128*16, 128*16, 128*16, 128*16, 128*16,
    128*16, 128*16, 128*16, 128*16, 128*16, 128*16, 128*16, 128*16,
    128*16, 128*16, 128*16, 128*16, 128*16, 128*16, 128*16, 128*16,
    128*16, 128*16, 128*16, 128*16, 128*16, 128*16, 128*16, 128*16,
    128*16, 128*16, 128*16, 128*16, 128*16, 128*16, 128*16, 128*16,
    128*16, 128*16, 128*16, 128*16, 128*16, 128*16, 128*16, 128*16
#endif
};
#endif  // CONFIG_VP9_HIGHBITDEPTH

unsigned int get_sby_perpixel_ssd(VP9_COMP *cpi,
                                  const struct buf_2d *ref,
                                  BLOCK_SIZE bs) {
  unsigned int sse;
  const unsigned int var = cpi->fn_ptr[bs].vf(ref->buf, ref->stride,
                                              VP9_VAR_OFFS, 0, &sse);
  return var;
}

unsigned int get_sby_perpixel_variance(VP9_COMP *cpi,
                                       const struct buf_2d *ref,
                                       BLOCK_SIZE bs) {
  const unsigned int var = get_sby_perpixel_ssd(cpi, ref, bs);
  return ROUND_POWER_OF_TWO(var, num_pels_log2_lookup[bs]);
}

#if CONFIG_VP9_HIGHBITDEPTH
unsigned int high_get_sby_perpixel_ssd(
    VP9_COMP *cpi, const struct buf_2d *ref, BLOCK_SIZE bs, int bd) {
  unsigned int var, sse;
  switch (bd) {
    case 10:
      var = cpi->fn_ptr[bs].vf(ref->buf, ref->stride,
                               CONVERT_TO_BYTEPTR(VP9_HIGH_VAR_OFFS_10),
                               0, &sse);
      break;
    case 12:
      var = cpi->fn_ptr[bs].vf(ref->buf, ref->stride,
                               CONVERT_TO_BYTEPTR(VP9_HIGH_VAR_OFFS_12),
                               0, &sse);
      break;
    case 8:
    default:
      var = cpi->fn_ptr[bs].vf(ref->buf, ref->stride,
                               CONVERT_TO_BYTEPTR(VP9_HIGH_VAR_OFFS_8),
                               0, &sse);
      break;
  }
  return var;
}

unsigned int high_get_sby_perpixel_variance(
    VP9_COMP *cpi, const struct buf_2d *ref, BLOCK_SIZE bs, int bd) {
  const unsigned int var = high_get_sby_perpixel_ssd(cpi, ref, bs, bd);
  return ROUND_POWER_OF_TWO(var, num_pels_log2_lookup[bs]);
}
#endif  // CONFIG_VP9_HIGHBITDEPTH

static unsigned int get_sby_perpixel_diff_variance(VP9_COMP *cpi,
                                                   const struct buf_2d *ref,
                                                   int mi_row, int mi_col,
                                                   BLOCK_SIZE bs) {
  const YV12_BUFFER_CONFIG *last = get_ref_frame_buffer(cpi, LAST_FRAME);
  const uint8_t* last_y = &last->y_buffer[mi_row * MI_SIZE * last->y_stride +
                                              mi_col * MI_SIZE];
  unsigned int sse;
  const unsigned int var = cpi->fn_ptr[bs].vf(ref->buf, ref->stride,
                                              last_y, last->y_stride, &sse);
  return ROUND_POWER_OF_TWO(var, num_pels_log2_lookup[bs]);
}

static BLOCK_SIZE get_rd_var_based_fixed_partition(VP9_COMP *cpi,
                                                   int mi_row,
                                                   int mi_col) {
  unsigned int var = get_sby_perpixel_diff_variance(cpi, &cpi->mb.plane[0].src,
                                                    mi_row, mi_col,
                                                    BLOCK_64X64);
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
static INLINE void set_modeinfo_offsets(VP9_COMMON *const cm,
                                        MACROBLOCKD *const xd,
                                        int mi_row,
                                        int mi_col) {
  const int idx_str = xd->mi_stride * mi_row + mi_col;
  xd->mi = cm->mi + idx_str;
  xd->mi[0].src_mi = &xd->mi[0];
}

static void set_offsets(VP9_COMP *cpi, const TileInfo *const tile,
                        int mi_row, int mi_col, BLOCK_SIZE bsize) {
  MACROBLOCK *const x = &cpi->mb;
  VP9_COMMON *const cm = &cpi->common;
  MACROBLOCKD *const xd = &x->e_mbd;
  MB_MODE_INFO *mbmi;
  const int mi_width = num_8x8_blocks_wide_lookup[bsize];
  const int mi_height = num_8x8_blocks_high_lookup[bsize];
  const struct segmentation *const seg = &cm->seg;

  set_skip_context(xd, mi_row, mi_col);

  set_modeinfo_offsets(cm, xd, mi_row, mi_col);

  mbmi = &xd->mi[0].src_mi->mbmi;

  // Set up destination pointers.
  vp9_setup_dst_planes(xd->plane, get_frame_new_buffer(cm), mi_row, mi_col);

  // Set up limit values for MV components.
  // Mv beyond the range do not produce new/different prediction block.
  x->mv_row_min = -(((mi_row + mi_height) * MI_SIZE) + VP9_INTERP_EXTEND);
  x->mv_col_min = -(((mi_col + mi_width) * MI_SIZE) + VP9_INTERP_EXTEND);
  x->mv_row_max = (cm->mi_rows - mi_row) * MI_SIZE + VP9_INTERP_EXTEND;
  x->mv_col_max = (cm->mi_cols - mi_col) * MI_SIZE + VP9_INTERP_EXTEND;

  // Set up distance of MB to edge of frame in 1/8th pel units.
  assert(!(mi_col & (mi_width - 1)) && !(mi_row & (mi_height - 1)));
  set_mi_row_col(xd, tile, mi_row, mi_height, mi_col, mi_width,
                 cm->mi_rows, cm->mi_cols);

  // Set up source buffers.
  vp9_setup_src_planes(x, cpi->Source, mi_row, mi_col);

  // R/D setup.
  x->rddiv = cpi->rd.RDDIV;
  x->rdmult = cpi->rd.RDMULT;

  // Setup segment ID.
  if (seg->enabled) {
    if (cpi->oxcf.aq_mode != VARIANCE_AQ) {
      const uint8_t *const map = seg->update_map ? cpi->segmentation_map
                                                 : cm->last_frame_seg_map;
      mbmi->segment_id = vp9_get_segment_id(cm, map, bsize, mi_row, mi_col);
    }
    vp9_init_plane_quantizers(cpi, x);

    x->encode_breakout = cpi->segment_encode_breakout[mbmi->segment_id];
  } else {
    mbmi->segment_id = 0;
    x->encode_breakout = cpi->encode_breakout;
  }
}

#if CONFIG_SUPERTX
static void set_offsets_supertx(VP9_COMP *cpi, const TileInfo *const tile,
                                int mi_row, int mi_col, BLOCK_SIZE bsize) {
  MACROBLOCK *const x = &cpi->mb;
  VP9_COMMON *const cm = &cpi->common;
  MACROBLOCKD *const xd = &x->e_mbd;
  const int mi_width = num_8x8_blocks_wide_lookup[bsize];
  const int mi_height = num_8x8_blocks_high_lookup[bsize];

  set_modeinfo_offsets(cm, xd, mi_row, mi_col);

  // Set up distance of MB to edge of frame in 1/8th pel units.
  assert(!(mi_col & (mi_width - 1)) && !(mi_row & (mi_height - 1)));
  set_mi_row_col(xd, tile, mi_row, mi_height, mi_col, mi_width,
                 cm->mi_rows, cm->mi_cols);
}

static void set_offsets_extend(VP9_COMP *cpi, const TileInfo *const tile,
                               int mi_row_pred, int mi_col_pred,
                               int mi_row_ori, int mi_col_ori,
                               BLOCK_SIZE bsize_pred, BLOCK_SIZE bsize_ori) {
  // Used in supertx
  // (mi_row_ori, mi_col_ori, bsize_ori): region for mv
  // (mi_row_pred, mi_col_pred, bsize_pred): region to predict
  MACROBLOCK *const x = &cpi->mb;
  VP9_COMMON *const cm = &cpi->common;
  MACROBLOCKD *const xd = &x->e_mbd;
  MB_MODE_INFO *mbmi;
  const int mi_width = num_8x8_blocks_wide_lookup[bsize_pred];
  const int mi_height = num_8x8_blocks_high_lookup[bsize_pred];
  const struct segmentation *const seg = &cm->seg;

  set_modeinfo_offsets(cm, xd, mi_row_ori, mi_col_ori);

  mbmi = &xd->mi[0].src_mi->mbmi;

  // Set up limit values for MV components.
  // Mv beyond the range do not produce new/different prediction block.
  x->mv_row_min = -(((mi_row_pred + mi_height) * MI_SIZE) + VP9_INTERP_EXTEND);
  x->mv_col_min = -(((mi_col_pred + mi_width) * MI_SIZE) + VP9_INTERP_EXTEND);
  x->mv_row_max = (cm->mi_rows - mi_row_pred) * MI_SIZE + VP9_INTERP_EXTEND;
  x->mv_col_max = (cm->mi_cols - mi_col_pred) * MI_SIZE + VP9_INTERP_EXTEND;

  // Set up distance of MB to edge of frame in 1/8th pel units.
  assert(!(mi_col_pred & (mi_width - 1)) && !(mi_row_pred & (mi_height - 1)));
  set_mi_row_col(xd, tile, mi_row_pred, mi_height, mi_col_pred, mi_width,
                 cm->mi_rows, cm->mi_cols);
  xd->up_available    = (mi_row_ori != 0);
  xd->left_available  = (mi_col_ori > tile->mi_col_start);

  // R/D setup.
  x->rddiv = cpi->rd.RDDIV;
  x->rdmult = cpi->rd.RDMULT;

  // Setup segment ID.
  if (seg->enabled) {
    if (cpi->oxcf.aq_mode != VARIANCE_AQ) {
      const uint8_t *const map = seg->update_map ? cpi->segmentation_map
                                                 : cm->last_frame_seg_map;
      mbmi->segment_id = vp9_get_segment_id(cm, map, bsize_ori,
                                            mi_row_ori, mi_col_ori);
    }
    vp9_init_plane_quantizers(cpi, x);

    x->encode_breakout = cpi->segment_encode_breakout[mbmi->segment_id];
  } else {
    mbmi->segment_id = 0;
    x->encode_breakout = cpi->encode_breakout;
  }
}
#endif  // CONFIG_SUPERTX

#if CONFIG_PALETTE
void copy_palette_info(PICK_MODE_CONTEXT *c, PICK_MODE_CONTEXT *p) {
  c->palette_buf_size = p->palette_buf_size;
  vpx_memcpy(c->palette_colors_buf, p->palette_colors_buf,
             c->palette_buf_size * sizeof(c->palette_colors_buf[0]));
  vpx_memcpy(c->palette_count_buf, p->palette_count_buf,
             c->palette_buf_size * sizeof(c->palette_count_buf[0]));
}
#endif

static void duplicate_mode_info_in_sb(VP9_COMMON *cm, MACROBLOCKD *xd,
                                      int mi_row, int mi_col,
                                      BLOCK_SIZE bsize) {
  const int block_width = num_8x8_blocks_wide_lookup[bsize];
  const int block_height = num_8x8_blocks_high_lookup[bsize];
  int i, j;
  for (j = 0; j < block_height; ++j)
    for (i = 0; i < block_width; ++i) {
      if (mi_row + j < cm->mi_rows && mi_col + i < cm->mi_cols)
        xd->mi[j * xd->mi_stride + i].src_mi = &xd->mi[0];
    }
}

static void set_block_size(VP9_COMP * const cpi,
                           int mi_row, int mi_col,
                           BLOCK_SIZE bsize) {
  if (cpi->common.mi_cols > mi_col && cpi->common.mi_rows > mi_row) {
    MACROBLOCKD *const xd = &cpi->mb.e_mbd;
    set_modeinfo_offsets(&cpi->common, xd, mi_row, mi_col);
    xd->mi[0].src_mi->mbmi.sb_type = bsize;
    duplicate_mode_info_in_sb(&cpi->common, xd, mi_row, mi_col, bsize);
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

typedef struct {
  partition_variance part_variances;
  var split[4];
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

#if CONFIG_EXT_CODING_UNIT_SIZE
typedef struct {
  partition_variance part_variances;
  v64x64 split[4];
} v128x128;
#endif

typedef struct {
  partition_variance *part_variances;
  var *split[4];
} variance_node;

typedef enum {
  V16X16,
  V32X32,
  V64X64,
#if CONFIG_EXT_CODING_UNIT_SIZE
  V128X128,
#endif
} TREE_LEVEL;

static void tree_to_node(void *data, BLOCK_SIZE bsize, variance_node *node) {
  int i;
  node->part_variances = NULL;
  vpx_memset(node->split, 0, sizeof(node->split));
  switch (bsize) {
#if CONFIG_EXT_CODING_UNIT_SIZE
    case BLOCK_128X128: {
      v128x128 *vt = (v128x128 *) data;
      node->part_variances = &vt->part_variances;
      for (i = 0; i < 4; i++)
        node->split[i] = &vt->split[i].part_variances.none;
      break;
    }
#endif
    case BLOCK_64X64: {
      v64x64 *vt = (v64x64 *) data;
      node->part_variances = &vt->part_variances;
      for (i = 0; i < 4; i++)
        node->split[i] = &vt->split[i].part_variances.none;
      break;
    }
    case BLOCK_32X32: {
      v32x32 *vt = (v32x32 *) data;
      node->part_variances = &vt->part_variances;
      for (i = 0; i < 4; i++)
        node->split[i] = &vt->split[i].part_variances.none;
      break;
    }
    case BLOCK_16X16: {
      v16x16 *vt = (v16x16 *) data;
      node->part_variances = &vt->part_variances;
      for (i = 0; i < 4; i++)
        node->split[i] = &vt->split[i].part_variances.none;
      break;
    }
    case BLOCK_8X8: {
      v8x8 *vt = (v8x8 *) data;
      node->part_variances = &vt->part_variances;
      for (i = 0; i < 4; i++)
        node->split[i] = &vt->split[i];
      break;
    }
    default: {
      assert(0);
      break;
    }
  }
}

// Set variance values given sum square error, sum error, count.
static void fill_variance(int64_t s2, int64_t s, int c, var *v) {
  v->sum_square_error = s2;
  v->sum_error = s;
  v->count = c;
  if (c > 0)
    v->variance = (int)(256 *
                        (v->sum_square_error - v->sum_error * v->sum_error /
                         v->count) / v->count);
  else
    v->variance = 0;
}

void sum_2_variances(const var *a, const var *b, var *r) {
  fill_variance(a->sum_square_error + b->sum_square_error,
                a->sum_error + b->sum_error, a->count + b->count, r);
}

static void fill_variance_tree(void *data, BLOCK_SIZE bsize) {
  variance_node node;
  tree_to_node(data, bsize, &node);
  sum_2_variances(node.split[0], node.split[1], &node.part_variances->horz[0]);
  sum_2_variances(node.split[2], node.split[3], &node.part_variances->horz[1]);
  sum_2_variances(node.split[0], node.split[2], &node.part_variances->vert[0]);
  sum_2_variances(node.split[1], node.split[3], &node.part_variances->vert[1]);
  sum_2_variances(&node.part_variances->vert[0], &node.part_variances->vert[1],
                  &node.part_variances->none);
}

static int set_vt_partitioning(VP9_COMP *cpi,
                               void *data,
                               BLOCK_SIZE bsize,
                               int mi_row,
                               int mi_col) {
  VP9_COMMON * const cm = &cpi->common;
  variance_node vt;
  const int block_width = num_8x8_blocks_wide_lookup[bsize];
  const int block_height = num_8x8_blocks_high_lookup[bsize];
  // TODO(debargha): Choose this more intelligently.
  const int threshold_multiplier = cm->frame_type == KEY_FRAME ? 64 : 4;
  int64_t threshold =
      (int64_t)(threshold_multiplier *
                vp9_convert_qindex_to_q(cm->base_qindex, cm->bit_depth));
  assert(block_height == block_width);
  tree_to_node(data, bsize, &vt);

  // Split none is available only if we have more than half a block size
  // in width and height inside the visible image.
  if (mi_col + block_width / 2 < cm->mi_cols &&
      mi_row + block_height / 2 < cm->mi_rows &&
      vt.part_variances->none.variance < threshold) {
    set_block_size(cpi, mi_row, mi_col, bsize);
    return 1;
  }

  // Only allow split for blocks above 16x16.
  if (bsize > BLOCK_16X16) {
    // Vertical split is available on all but the bottom border.
    if (mi_row + block_height / 2 < cm->mi_rows &&
        vt.part_variances->vert[0].variance < threshold &&
        vt.part_variances->vert[1].variance < threshold) {
      BLOCK_SIZE subsize = get_subsize(bsize, PARTITION_VERT);
      set_block_size(cpi, mi_row, mi_col, subsize);
      set_block_size(cpi, mi_row, mi_col + block_width / 2, subsize);
      return 1;
    }

    // Horizontal split is available on all but the right border.
    if (mi_col + block_width / 2 < cm->mi_cols &&
        vt.part_variances->horz[0].variance < threshold &&
        vt.part_variances->horz[1].variance < threshold) {
      BLOCK_SIZE subsize = get_subsize(bsize, PARTITION_HORZ);
      set_block_size(cpi, mi_row, mi_col, subsize);
      set_block_size(cpi, mi_row + block_height / 2, mi_col, subsize);
      return 1;
    }
  }

  // This will only allow 8x8 if the 16x16 variance is very large.
  if (bsize == BLOCK_16X16) {
    if (mi_col + block_width / 2 < cm->mi_cols &&
        mi_row + block_height / 2 < cm->mi_rows &&
        vt.part_variances->none.variance < (threshold << 6)) {
      set_block_size(cpi, mi_row, mi_col, bsize);
      return 1;
    }
  }
  return 0;
}

// This function chooses partitioning based on the variance
// between source and reconstructed last, where variance is
// computed for 8x8 downsampled inputs. Some things to check:
// using the last source rather than reconstructed last, and
// allowing for small downsampling (4x4 or 2x2) for selection
// of smaller block sizes (i.e., < 16x16).
static void choose_partitioning(VP9_COMP *cpi,
                                const TileInfo *const tile,
                                int mi_row, int mi_col) {
  VP9_COMMON * const cm = &cpi->common;
  MACROBLOCK *x = &cpi->mb;
  MACROBLOCKD *xd = &cpi->mb.e_mbd;

  int i, j, k;
  v64x64 vt;
  uint8_t *s;
  const uint8_t *d;
  int sp;
  int dp;
  int pixels_wide = 64, pixels_high = 64;
  int_mv nearest_mv, near_mv;
  const YV12_BUFFER_CONFIG *yv12 = get_ref_frame_buffer(cpi, LAST_FRAME);
  const struct scale_factors *const sf = &cm->frame_refs[LAST_FRAME - 1].sf;

  vp9_clear_system_state();
  vp9_zero(vt);
  set_offsets(cpi, tile, mi_row, mi_col, BLOCK_LARGEST);
#if CONFIG_EXT_CODING_UNIT_SIZE
  printf("Not yet implemented: choose_partitioning\n");
  exit(-1);
#endif

  if (xd->mb_to_right_edge < 0)
    pixels_wide += (xd->mb_to_right_edge >> 3);
  if (xd->mb_to_bottom_edge < 0)
    pixels_high += (xd->mb_to_bottom_edge >> 3);

  s = x->plane[0].src.buf;
  sp = x->plane[0].src.stride;

  if (cm->frame_type != KEY_FRAME) {
    vp9_setup_pre_planes(xd, 0, yv12, mi_row, mi_col, sf);

    xd->mi[0].src_mi->mbmi.ref_frame[0] = LAST_FRAME;
#if CONFIG_INTERINTRA
    xd->mi[0].src_mi->mbmi.ref_frame[1] = NONE;
#endif  // CONFIG_INTERINTRA
    xd->mi[0].src_mi->mbmi.sb_type = BLOCK_LARGEST;
    vp9_find_best_ref_mvs(xd, cm->allow_high_precision_mv,
                          xd->mi[0].src_mi->mbmi.ref_mvs[LAST_FRAME],
                          &nearest_mv, &near_mv);

    xd->mi[0].src_mi->mbmi.mv[0] = nearest_mv;
    vp9_build_inter_predictors_sby(xd, mi_row, mi_col, BLOCK_LARGEST);

    d = xd->plane[0].dst.buf;
    dp = xd->plane[0].dst.stride;
  } else {
    d = VP9_VAR_OFFS;
    dp = 0;
#if CONFIG_VP9_HIGHBITDEPTH
    if (xd->cur_buf->flags & YV12_FLAG_HIGHBITDEPTH) {
      switch (xd->bd) {
        case 10:
          d = CONVERT_TO_BYTEPTR(VP9_HIGH_VAR_OFFS_10);
          break;
        case 12:
          d = CONVERT_TO_BYTEPTR(VP9_HIGH_VAR_OFFS_12);
          break;
        case 8:
        default:
          d = CONVERT_TO_BYTEPTR(VP9_HIGH_VAR_OFFS_8);
          break;
      }
    }
#endif  // CONFIG_VP9_HIGHBITDEPTH
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

        if (x_idx < pixels_wide && y_idx < pixels_high) {
          int s_avg, d_avg;
#if CONFIG_VP9_HIGHBITDEPTH
          if (xd->cur_buf->flags & YV12_FLAG_HIGHBITDEPTH) {
            s_avg = vp9_highbd_avg_8x8(s + y_idx * sp + x_idx, sp);
            d_avg = vp9_highbd_avg_8x8(d + y_idx * dp + x_idx, dp);
          } else {
            s_avg = vp9_avg_8x8(s + y_idx * sp + x_idx, sp);
            d_avg = vp9_avg_8x8(d + y_idx * dp + x_idx, dp);
          }
#else
          s_avg = vp9_avg_8x8(s + y_idx * sp + x_idx, sp);
          d_avg = vp9_avg_8x8(d + y_idx * dp + x_idx, dp);
#endif
          sum = s_avg - d_avg;
          sse = sum * sum;
        }
        // For an 8x8 block we have just one value the average of all 64
        // pixels,  so use 1.   This means of course that there is no variance
        // in an 8x8 block.
        fill_variance(sse, sum, 1, &vst->split[k].part_variances.none);
      }
    }
  }
  // Fill the rest of the variance tree by summing split partition values.
  for (i = 0; i < 4; i++) {
    for (j = 0; j < 4; j++) {
      fill_variance_tree(&vt.split[i].split[j], BLOCK_16X16);
    }
    fill_variance_tree(&vt.split[i], BLOCK_32X32);
  }
  fill_variance_tree(&vt, BLOCK_64X64);

  // Now go through the entire structure,  splitting every block size until
  // we get to one that's got a variance lower than our threshold,  or we
  // hit 8x8.
  if ( mi_col + 8 > cm->mi_cols || mi_row + 8 > cm->mi_rows ||
      !set_vt_partitioning(cpi, &vt, BLOCK_64X64, mi_row, mi_col)) {
    for (i = 0; i < 4; ++i) {
      const int x32_idx = ((i & 1) << 2);
      const int y32_idx = ((i >> 1) << 2);
      if (!set_vt_partitioning(cpi, &vt.split[i], BLOCK_32X32,
                               (mi_row + y32_idx), (mi_col + x32_idx))) {
        for (j = 0; j < 4; ++j) {
          const int x16_idx = ((j & 1) << 1);
          const int y16_idx = ((j >> 1) << 1);
          // NOTE: Since this uses 8x8 downsampling for variance calculation
          // we cannot really select block size 8x8 (or even 8x16/16x8),
          // since we do not sufficient samples for variance.
          // For now, 8x8 partition is only set if the variance of the 16x16
          // block is very high. This is controlled in set_vt_partitioning.
          if (!set_vt_partitioning(cpi, &vt.split[i].split[j],
                                   BLOCK_16X16,
                                   mi_row + y32_idx + y16_idx,
                                   mi_col + x32_idx + x16_idx)) {
            for (k = 0; k < 4; ++k) {
              const int x8_idx = (k & 1);
              const int y8_idx = (k >> 1);
              set_block_size(cpi,
                             (mi_row + y32_idx + y16_idx + y8_idx),
                             (mi_col + x32_idx + x16_idx + x8_idx),
                             BLOCK_8X8);
            }
          }
        }
      }
    }
  }
}

static void update_state(VP9_COMP *cpi, PICK_MODE_CONTEXT *ctx,
                         int mi_row, int mi_col, BLOCK_SIZE bsize,
                         int output_enabled) {
  int i, x_idx, y;
  VP9_COMMON *const cm = &cpi->common;
  RD_OPT *const rd_opt = &cpi->rd;
  MACROBLOCK *const x = &cpi->mb;
  MACROBLOCKD *const xd = &x->e_mbd;
  struct macroblock_plane *const p = x->plane;
  struct macroblockd_plane *const pd = xd->plane;
  MODE_INFO *mi = &ctx->mic;
  MB_MODE_INFO *const mbmi = &xd->mi[0].src_mi->mbmi;
  MODE_INFO *mi_addr = &xd->mi[0];
  const struct segmentation *const seg = &cm->seg;

  const int mis = cm->mi_stride;
  const int mi_width = num_8x8_blocks_wide_lookup[bsize];
  const int mi_height = num_8x8_blocks_high_lookup[bsize];
  int max_plane;

#if !CONFIG_SUPERTX
  assert(mi->mbmi.sb_type == bsize);
#endif

  *mi_addr = *mi;
  mi_addr->src_mi = mi_addr;

  // If segmentation in use
  if (seg->enabled && output_enabled) {
    // For in frame complexity AQ copy the segment id from the segment map.
    if (cpi->oxcf.aq_mode == COMPLEXITY_AQ) {
      const uint8_t *const map = seg->update_map ? cpi->segmentation_map
                                                 : cm->last_frame_seg_map;
      mi_addr->mbmi.segment_id =
        vp9_get_segment_id(cm, map, bsize, mi_row, mi_col);
    }
    // Else for cyclic refresh mode update the segment map, set the segment id
    // and then update the quantizer.
    if (cpi->oxcf.aq_mode == CYCLIC_REFRESH_AQ) {
      vp9_cyclic_refresh_update_segment(cpi, &xd->mi[0].src_mi->mbmi,
                                        mi_row, mi_col, bsize, 1);
    }
  }

  max_plane = is_inter_block(mbmi) ? MAX_MB_PLANE : 1;
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

#if CONFIG_PALETTE
  for (i = 0; i < 2; i++) {
    pd[i].color_index_map = ctx->color_index_map[i];
  }
#endif  // CONFIG_PALETTE

  // Restore the coding context of the MB to that that was in place
  // when the mode was picked for it
  for (y = 0; y < mi_height; y++)
    for (x_idx = 0; x_idx < mi_width; x_idx++)
      if ((xd->mb_to_right_edge >> (3 + MI_SIZE_LOG2)) + mi_width > x_idx
        && (xd->mb_to_bottom_edge >> (3 + MI_SIZE_LOG2)) + mi_height > y) {
        xd->mi[x_idx + y * mis].src_mi = mi_addr;
      }

  if (cpi->oxcf.aq_mode)
    vp9_init_plane_quantizers(cpi, x);

  // FIXME(rbultje) I'm pretty sure this should go to the end of this block
  // (i.e. after the output_enabled)
#if CONFIG_TX64X64
  if (bsize < BLOCK_64X64) {
    if (bsize < BLOCK_32X32) {
      if (bsize < BLOCK_16X16) {
        ctx->tx_rd_diff[ALLOW_16X16] = ctx->tx_rd_diff[ALLOW_8X8];
      }
      ctx->tx_rd_diff[ALLOW_32X32] = ctx->tx_rd_diff[ALLOW_16X16];
    }
    ctx->tx_rd_diff[ALLOW_64X64] = ctx->tx_rd_diff[ALLOW_32X32];
  }
#else
  if (bsize < BLOCK_32X32) {
    if (bsize < BLOCK_16X16)
      ctx->tx_rd_diff[ALLOW_16X16] = ctx->tx_rd_diff[ALLOW_8X8];
    ctx->tx_rd_diff[ALLOW_32X32] = ctx->tx_rd_diff[ALLOW_16X16];
  }
#endif

  if (is_inter_block(mbmi) && mbmi->sb_type < BLOCK_8X8) {
    mbmi->mv[0].as_int = mi->bmi[3].as_mv[0].as_int;
    mbmi->mv[1].as_int = mi->bmi[3].as_mv[1].as_int;
  }

  x->skip = ctx->skip;
#if CONFIG_SR_MODE
  vpx_memcpy(x->zcoeff_blk[mbmi->sr ? TX_SIZES : mbmi->tx_size],
             ctx->zcoeff_blk,
             sizeof(uint8_t) * ctx->num_4x4_blk);
#else  // CONFIG_SR_MODE
  vpx_memcpy(x->zcoeff_blk[mbmi->tx_size], ctx->zcoeff_blk,
             sizeof(uint8_t) * ctx->num_4x4_blk);
#endif  // CONFIG_SR_MODE

  if (!output_enabled)
    return;

  if (!vp9_segfeature_active(&cm->seg, mbmi->segment_id, SEG_LVL_SKIP)) {
    for (i = 0; i < TX_MODES; i++)
      rd_opt->tx_select_diff[i] += ctx->tx_rd_diff[i];
  }

#if CONFIG_INTERNAL_STATS
  if (frame_is_intra_only(cm)) {
    static const int kf_mode_index[] = {
      THR_DC        /*DC_PRED*/,
      THR_V_PRED    /*V_PRED*/,
      THR_H_PRED    /*H_PRED*/,
      THR_D45_PRED  /*D45_PRED*/,
      THR_D135_PRED /*D135_PRED*/,
      THR_D117_PRED /*D117_PRED*/,
      THR_D153_PRED /*D153_PRED*/,
      THR_D207_PRED /*D207_PRED*/,
      THR_D63_PRED  /*D63_PRED*/,
      THR_TM        /*TM_PRED*/,
    };
    ++cpi->mode_chosen_counts[kf_mode_index[mbmi->mode]];
  } else {
    // Note how often each mode chosen as best
    ++cpi->mode_chosen_counts[ctx->best_mode_index];
  }
#endif
  if (!frame_is_intra_only(cm)) {
#if CONFIG_COPY_MODE
    if (is_inter_block(mbmi) && mbmi->copy_mode == NOREF) {
#else
    if (is_inter_block(mbmi)) {
#endif  // CONFIG_COPY_MODE
#if CONFIG_GLOBAL_MOTION
      if (bsize >= BLOCK_8X8) {
#if CONFIG_NEW_INTER
        if (mbmi->mode == ZEROMV || mbmi->mode == ZERO_ZEROMV) {
          ++cpi->global_motion_used[mbmi->ref_frame[0]];
          if (mbmi->mode == ZERO_ZEROMV)
            ++cpi->global_motion_used[mbmi->ref_frame[1]];
        }
#else
        if (mbmi->mode == ZEROMV) {
          ++cpi->global_motion_used[mbmi->ref_frame[0]];
          if (has_second_ref(mbmi))
            ++cpi->global_motion_used[mbmi->ref_frame[1]];
        }
#endif  // CONFIG_NEW_INTER
      } else {
        const int num_4x4_w = num_4x4_blocks_wide_lookup[bsize];
        const int num_4x4_h = num_4x4_blocks_high_lookup[bsize];
        int idx, idy;
        for (idy = 0; idy < 2; idy += num_4x4_h) {
          for (idx = 0; idx < 2; idx += num_4x4_w) {
            const int j = idy * 2 + idx;
            const PREDICTION_MODE b_mode = mi->bmi[j].as_mode;
#if CONFIG_NEW_INTER
            if (b_mode == ZEROMV || b_mode == ZERO_ZEROMV) {
              ++cpi->global_motion_used[mbmi->ref_frame[0]];
              if (b_mode == ZERO_ZEROMV)
                ++cpi->global_motion_used[mbmi->ref_frame[1]];
            }
#else
            if (b_mode == ZEROMV) {
              ++cpi->global_motion_used[mbmi->ref_frame[0]];
              if (has_second_ref(mbmi))
                ++cpi->global_motion_used[mbmi->ref_frame[1]];
            }
#endif  // CONFIG_NEW_INTER
          }
        }
      }
#endif  // CONFIG_GLOBAL_MOTION
    }

    rd_opt->comp_pred_diff[SINGLE_REFERENCE] += ctx->single_pred_diff;
    rd_opt->comp_pred_diff[COMPOUND_REFERENCE] += ctx->comp_pred_diff;
    rd_opt->comp_pred_diff[REFERENCE_MODE_SELECT] += ctx->hybrid_pred_diff;

    for (i = 0; i < SWITCHABLE_FILTER_CONTEXTS; ++i)
      rd_opt->filter_diff[i] += ctx->best_filter_diff[i];
  }
}

#if CONFIG_SUPERTX
static void update_state_supertx(VP9_COMP *cpi, PICK_MODE_CONTEXT *ctx,
                                 int mi_row, int mi_col, BLOCK_SIZE bsize,
                                 int output_enabled) {
  int i, y, x_idx;
  VP9_COMMON *const cm = &cpi->common;
  RD_OPT *const rd_opt = &cpi->rd;
  MACROBLOCK *const x = &cpi->mb;
  MACROBLOCKD *const xd = &x->e_mbd;
  MODE_INFO *mi = &ctx->mic;
  MB_MODE_INFO *const mbmi = &xd->mi[0].src_mi->mbmi;
  MODE_INFO *mi_addr = &xd->mi[0];
  const struct segmentation *const seg = &cm->seg;
  const int mis = cm->mi_stride;
  const int mi_width = num_8x8_blocks_wide_lookup[bsize];
  const int mi_height = num_8x8_blocks_high_lookup[bsize];

  *mi_addr = *mi;
  mi_addr->src_mi = mi_addr;
  assert(is_inter_block(mbmi));

  // If segmentation in use
  if (seg->enabled && output_enabled) {
    // For in frame complexity AQ copy the segment id from the segment map.
    if (cpi->oxcf.aq_mode == COMPLEXITY_AQ) {
      const uint8_t *const map = seg->update_map ? cpi->segmentation_map
                                                 : cm->last_frame_seg_map;
      mi_addr->mbmi.segment_id =
        vp9_get_segment_id(cm, map, bsize, mi_row, mi_col);
    } else if (cpi->oxcf.aq_mode == CYCLIC_REFRESH_AQ) {
      // Else for cyclic refresh mode update the segment map, set the segment id
      // and then update the quantizer.
      vp9_cyclic_refresh_update_segment(cpi, &xd->mi[0].mbmi,
                                        mi_row, mi_col, bsize, 1);
      vp9_init_plane_quantizers(cpi, x);
    }
  }

  // Restore the coding context of the MB to that that was in place
  // when the mode was picked for it
  for (y = 0; y < mi_height; y++)
    for (x_idx = 0; x_idx < mi_width; x_idx++)
      if ((xd->mb_to_right_edge >> (3 + MI_SIZE_LOG2)) + mi_width > x_idx
        && (xd->mb_to_bottom_edge >> (3 + MI_SIZE_LOG2)) + mi_height > y) {
        xd->mi[x_idx + y * mis].src_mi = mi_addr;
      }

  if (cpi->oxcf.aq_mode)
    vp9_init_plane_quantizers(cpi, x);

  if (is_inter_block(mbmi) && mbmi->sb_type < BLOCK_8X8) {
    mbmi->mv[0].as_int = mi->bmi[3].as_mv[0].as_int;
    mbmi->mv[1].as_int = mi->bmi[3].as_mv[1].as_int;
  }

  x->skip = ctx->skip;
#if 0  // CONFIG_SR_MODE  // supertx????  // debugtest
  vpx_memcpy(x->zcoeff_blk[mbmi->sr ? TX_SIZES : mbmi->tx_size],
             ctx->zcoeff_blk, sizeof(uint8_t) * ctx->num_4x4_blk);
#else
  vpx_memcpy(x->zcoeff_blk[mbmi->tx_size], ctx->zcoeff_blk,
             sizeof(uint8_t) * ctx->num_4x4_blk);
#endif

  if (!output_enabled)
    return;

  if (!frame_is_intra_only(cm)) {
    rd_opt->comp_pred_diff[SINGLE_REFERENCE] += ctx->single_pred_diff;
    rd_opt->comp_pred_diff[COMPOUND_REFERENCE] += ctx->comp_pred_diff;
    rd_opt->comp_pred_diff[REFERENCE_MODE_SELECT] += ctx->hybrid_pred_diff;

    for (i = 0; i < SWITCHABLE_FILTER_CONTEXTS; ++i)
      rd_opt->filter_diff[i] += ctx->best_filter_diff[i];
  }
}

static void update_state_sb_supertx(VP9_COMP *cpi, const TileInfo *const tile,
                                    int mi_row, int mi_col,
                                    BLOCK_SIZE bsize,
                                    int output_enabled, PC_TREE *pc_tree) {
  VP9_COMMON *const cm = &cpi->common;
  MACROBLOCK *const x = &cpi->mb;
  MACROBLOCKD *const xd = &x->e_mbd;
  struct macroblock_plane *const p = x->plane;
  struct macroblockd_plane *const pd = xd->plane;
  int bsl = b_width_log2_lookup[bsize], hbs = (1 << bsl) / 4;
  PARTITION_TYPE partition = pc_tree->partitioning;
  BLOCK_SIZE subsize = get_subsize(bsize, partition);
  int i;
#if CONFIG_EXT_PARTITION
  BLOCK_SIZE bsize2 = get_subsize(bsize, PARTITION_SPLIT);
#endif

  if (mi_row >= cm->mi_rows || mi_col >= cm->mi_cols)
    return;

  switch (partition) {
    case PARTITION_NONE:
      set_offsets_supertx(cpi, tile, mi_row, mi_col, subsize);
      update_state_supertx(cpi, &pc_tree->none, mi_row, mi_col,
                           subsize, output_enabled);
      break;
    case PARTITION_VERT:
      set_offsets_supertx(cpi, tile, mi_row, mi_col, subsize);
      update_state_supertx(cpi, &pc_tree->vertical[0], mi_row, mi_col,
                           subsize, output_enabled);
      if (mi_col + hbs < cm->mi_cols && bsize > BLOCK_8X8) {
        set_offsets_supertx(cpi, tile, mi_row, mi_col + hbs, subsize);
        update_state_supertx(cpi, &pc_tree->vertical[1], mi_row, mi_col + hbs,
                             subsize, output_enabled);
      }
      break;
    case PARTITION_HORZ:
      set_offsets_supertx(cpi, tile, mi_row, mi_col, subsize);
      update_state_supertx(cpi, &pc_tree->horizontal[0], mi_row, mi_col,
                           subsize, output_enabled);
      if (mi_row + hbs < cm->mi_rows && bsize > BLOCK_8X8) {
        set_offsets_supertx(cpi, tile, mi_row + hbs, mi_col, subsize);
        update_state_supertx(cpi, &pc_tree->horizontal[1], mi_row + hbs, mi_col,
                             subsize, output_enabled);
      }
      break;
    case PARTITION_SPLIT:
      if (bsize == BLOCK_8X8) {
        set_offsets_supertx(cpi, tile, mi_row, mi_col, subsize);
        update_state_supertx(cpi, pc_tree->leaf_split[0], mi_row, mi_col,
                             subsize, output_enabled);
      } else {
        set_offsets_supertx(cpi, tile, mi_row, mi_col, subsize);
        update_state_sb_supertx(cpi, tile, mi_row, mi_col, subsize,
                                output_enabled, pc_tree->split[0]);
        set_offsets_supertx(cpi, tile, mi_row, mi_col + hbs, subsize);
        update_state_sb_supertx(cpi, tile, mi_row, mi_col + hbs, subsize,
                                output_enabled, pc_tree->split[1]);
        set_offsets_supertx(cpi, tile, mi_row + hbs, mi_col, subsize);
        update_state_sb_supertx(cpi, tile, mi_row + hbs, mi_col, subsize,
                                output_enabled, pc_tree->split[2]);
        set_offsets_supertx(cpi, tile, mi_row + hbs, mi_col + hbs, subsize);
        update_state_sb_supertx(cpi, tile, mi_row + hbs, mi_col + hbs, subsize,
                                output_enabled, pc_tree->split[3]);
      }
      break;
#if CONFIG_EXT_PARTITION
    case PARTITION_HORZ_A:
      set_offsets_supertx(cpi, tile, mi_row, mi_col, bsize2);
      update_state_supertx(cpi, &pc_tree->horizontala[0], mi_row, mi_col,
                           bsize2, output_enabled);
      set_offsets_supertx(cpi, tile, mi_row, mi_col + hbs, bsize2);
      update_state_supertx(cpi, &pc_tree->horizontala[1], mi_row, mi_col + hbs,
                           bsize2, output_enabled);
      set_offsets_supertx(cpi, tile, mi_row + hbs, mi_col, subsize);
      update_state_supertx(cpi, &pc_tree->horizontala[2], mi_row + hbs, mi_col,
                             subsize, output_enabled);
      break;
    case PARTITION_HORZ_B:
      set_offsets_supertx(cpi, tile, mi_row, mi_col, subsize);
      update_state_supertx(cpi, &pc_tree->horizontalb[0], mi_row, mi_col,
                           subsize, output_enabled);
      set_offsets_supertx(cpi, tile, mi_row + hbs, mi_col, bsize2);
      update_state_supertx(cpi, &pc_tree->horizontalb[1], mi_row + hbs, mi_col,
                           bsize2, output_enabled);
      set_offsets_supertx(cpi, tile, mi_row + hbs, mi_col + hbs, bsize2);
      update_state_supertx(cpi, &pc_tree->horizontalb[2], mi_row + hbs,
                           mi_col + hbs, bsize2, output_enabled);
      break;
    case PARTITION_VERT_A:
      set_offsets_supertx(cpi, tile, mi_row, mi_col, bsize2);
      update_state_supertx(cpi, &pc_tree->verticala[0], mi_row, mi_col, bsize2,
                           output_enabled);
      set_offsets_supertx(cpi, tile, mi_row + hbs, mi_col, bsize2);
      update_state_supertx(cpi, &pc_tree->verticala[1], mi_row + hbs, mi_col,
                           bsize2, output_enabled);
      set_offsets_supertx(cpi, tile, mi_row, mi_col + hbs, subsize);
      update_state_supertx(cpi, &pc_tree->verticala[2], mi_row, mi_col + hbs,
                           subsize, output_enabled);
      break;
    case PARTITION_VERT_B:
      set_offsets_supertx(cpi, tile, mi_row, mi_col, subsize);
      update_state_supertx(cpi, &pc_tree->verticalb[0], mi_row, mi_col,
                           subsize, output_enabled);
      set_offsets_supertx(cpi, tile, mi_row, mi_col + hbs, bsize2);
      update_state_supertx(cpi, &pc_tree->verticalb[1], mi_row, mi_col + hbs,
                           bsize2, output_enabled);
      set_offsets_supertx(cpi, tile, mi_row + hbs, mi_col + hbs, bsize2);
      update_state_supertx(cpi, &pc_tree->verticalb[2], mi_row + hbs,
                           mi_col + hbs, bsize2, output_enabled);
      break;
#endif
    default:
      assert(0);
  }

  for (i = 0; i < MAX_MB_PLANE; ++i) {
    p[i].coeff = (&pc_tree->none)->coeff_pbuf[i][1];
    p[i].qcoeff = (&pc_tree->none)->qcoeff_pbuf[i][1];
    pd[i].dqcoeff = (&pc_tree->none)->dqcoeff_pbuf[i][1];
    p[i].eobs = (&pc_tree->none)->eobs_pbuf[i][1];
  }
}

static void update_supertx_param(VP9_COMP *cpi, PICK_MODE_CONTEXT *ctx,
#if CONFIG_EXT_TX
                                 int best_tx,
#endif
#if CONFIG_NEW_QUANT && QUANT_PROFILES > 1
                                 int dq_index,
#endif  // CONFIG_NEW_QUANT && QUANT_PROFILES > 1
                                 TX_SIZE supertx_size) {
  MACROBLOCK *const x = &cpi->mb;

  ctx->mic.mbmi.tx_size = supertx_size;
  vpx_memcpy(ctx->zcoeff_blk, x->zcoeff_blk[supertx_size],
             sizeof(uint8_t) * ctx->num_4x4_blk);
  ctx->skip = x->skip;
#if CONFIG_EXT_TX
  ctx->mic.mbmi.ext_txfrm = best_tx;
#endif  // CONFIG_EXT_TX
#if CONFIG_NEW_QUANT && QUANT_PROFILES > 1
  ctx->mic.mbmi.dq_off_index = dq_index;
#endif  // CONFIG_NEW_QUANT && QUANT_PROFILES > 1
#if CONFIG_TX_SKIP
  ctx->mic.mbmi.tx_skip[0] = 0;
  ctx->mic.mbmi.tx_skip[1] = 0;
#endif  // CONFIG_TX_SKIP
}

static void update_supertx_param_sb(VP9_COMP *cpi, int mi_row, int mi_col,
                                    BLOCK_SIZE bsize,
#if CONFIG_EXT_TX
                                    int best_tx,
#endif
#if CONFIG_NEW_QUANT && QUANT_PROFILES > 1
                                    int dq_index,
#endif  // CONFIG_NEW_QUANT && QUANT_PROFILES > 1
                                    TX_SIZE supertx_size, PC_TREE *pc_tree) {
  VP9_COMMON *const cm = &cpi->common;
  int bsl = b_width_log2_lookup[bsize], hbs = (1 << bsl) / 4;
  PARTITION_TYPE partition = pc_tree->partitioning;
  BLOCK_SIZE subsize = get_subsize(bsize, partition);
#if CONFIG_EXT_PARTITION
  int i;
#endif

  if (mi_row >= cm->mi_rows || mi_col >= cm->mi_cols)
    return;

  switch (partition) {
    case PARTITION_NONE:
      update_supertx_param(cpi, &pc_tree->none,
#if CONFIG_EXT_TX
                           best_tx,
#endif
#if CONFIG_NEW_QUANT && QUANT_PROFILES > 1
                           dq_index,
#endif
                           supertx_size);
      break;
    case PARTITION_VERT:
      update_supertx_param(cpi, &pc_tree->vertical[0],
#if CONFIG_EXT_TX
                           best_tx,
#endif
#if CONFIG_NEW_QUANT && QUANT_PROFILES > 1
                           dq_index,
#endif
                           supertx_size);
      if (mi_col + hbs < cm->mi_cols && bsize > BLOCK_8X8)
        update_supertx_param(cpi, &pc_tree->vertical[1],
#if CONFIG_EXT_TX
                             best_tx,
#endif
#if CONFIG_NEW_QUANT && QUANT_PROFILES > 1
                             dq_index,
#endif
                             supertx_size);
      break;
    case PARTITION_HORZ:
      update_supertx_param(cpi, &pc_tree->horizontal[0],
#if CONFIG_EXT_TX
                           best_tx,
#endif
#if CONFIG_NEW_QUANT && QUANT_PROFILES > 1
                           dq_index,
#endif
                           supertx_size);
      if (mi_row + hbs < cm->mi_rows && bsize > BLOCK_8X8)
        update_supertx_param(cpi, &pc_tree->horizontal[1],
#if CONFIG_EXT_TX
                             best_tx,
#endif
#if CONFIG_NEW_QUANT && QUANT_PROFILES > 1
                             dq_index,
#endif
                             supertx_size);
      break;
    case PARTITION_SPLIT:
      if (bsize == BLOCK_8X8) {
        update_supertx_param(cpi, pc_tree->leaf_split[0],
#if CONFIG_EXT_TX
                             best_tx,
#endif
#if CONFIG_NEW_QUANT && QUANT_PROFILES > 1
                             dq_index,
#endif
                             supertx_size);
      } else {
        update_supertx_param_sb(cpi, mi_row, mi_col, subsize,
#if CONFIG_EXT_TX
                                best_tx,
#endif
#if CONFIG_NEW_QUANT && QUANT_PROFILES > 1
                                dq_index,
#endif
                                supertx_size, pc_tree->split[0]);
        update_supertx_param_sb(cpi, mi_row, mi_col + hbs, subsize,
#if CONFIG_EXT_TX
                                best_tx,
#endif
#if CONFIG_NEW_QUANT && QUANT_PROFILES > 1
                                dq_index,
#endif
                                supertx_size, pc_tree->split[1]);
        update_supertx_param_sb(cpi, mi_row + hbs, mi_col, subsize,
#if CONFIG_EXT_TX
                                best_tx,
#endif
#if CONFIG_NEW_QUANT && QUANT_PROFILES > 1
                                dq_index,
#endif
                                supertx_size, pc_tree->split[2]);
        update_supertx_param_sb(cpi, mi_row + hbs, mi_col + hbs, subsize,
#if CONFIG_EXT_TX
                                best_tx,
#endif
#if CONFIG_NEW_QUANT && QUANT_PROFILES > 1
                                dq_index,
#endif
                                supertx_size, pc_tree->split[3]);
      }
      break;
#if CONFIG_EXT_PARTITION
    case PARTITION_HORZ_A:
      for ( i = 0; i < 3; i++)
        update_supertx_param(cpi, &pc_tree->horizontala[i],
#if CONFIG_EXT_TX
                            best_tx,
#endif
#if CONFIG_NEW_QUANT && QUANT_PROFILES > 1
                            dq_index,
#endif
                            supertx_size);
      break;
    case PARTITION_HORZ_B:
      for ( i = 0; i < 3; i++)
        update_supertx_param(cpi, &pc_tree->horizontalb[i],
#if CONFIG_EXT_TX
                            best_tx,
#endif
#if CONFIG_NEW_QUANT && QUANT_PROFILES > 1
                            dq_index,
#endif
                            supertx_size);
      break;
    case PARTITION_VERT_A:
      for ( i = 0; i < 3; i++)
        update_supertx_param(cpi, &pc_tree->verticala[i],
#if CONFIG_EXT_TX
                            best_tx,
#endif
#if CONFIG_NEW_QUANT && QUANT_PROFILES > 1
                            dq_index,
#endif
                            supertx_size);
      break;
    case PARTITION_VERT_B:
      for ( i = 0; i < 3; i++)
        update_supertx_param(cpi, &pc_tree->verticalb[i],
#if CONFIG_EXT_TX
                            best_tx,
#endif
#if CONFIG_NEW_QUANT && QUANT_PROFILES > 1
                            dq_index,
#endif
                            supertx_size);
      break;
#endif
    default:
      assert(0);
  }
}
#endif  // CONFIG_SUPERTX

void vp9_setup_src_planes(MACROBLOCK *x, const YV12_BUFFER_CONFIG *src,
                          int mi_row, int mi_col) {
  uint8_t *const buffers[3] = {src->y_buffer, src->u_buffer, src->v_buffer };
  const int strides[3] = {src->y_stride, src->uv_stride, src->uv_stride };
  const int widths[3] = {src->y_crop_width, src->uv_crop_width,
                        src->uv_crop_width};
  const int heights[3] = {src->y_crop_height, src->uv_crop_height,
                          src->uv_crop_height};
  int i;

  // Set current frame pointer.
  x->e_mbd.cur_buf = src;

  for (i = 0; i < MAX_MB_PLANE; i++)
    setup_pred_plane(&x->plane[i].src, widths[i], heights[i],
                     buffers[i], strides[i], mi_row, mi_col,
                     NULL, x->e_mbd.plane[i].subsampling_x,
                     x->e_mbd.plane[i].subsampling_y);
}

static void rd_pick_sb_modes(VP9_COMP *cpi, const TileInfo *const tile,
                             int mi_row, int mi_col, RD_COST *rd_cost,
#if CONFIG_SUPERTX
                             int *totalrate_nocoef,
#endif
                             BLOCK_SIZE bsize, PICK_MODE_CONTEXT *ctx,
                             int64_t best_rd) {
  VP9_COMMON *const cm = &cpi->common;
  MACROBLOCK *const x = &cpi->mb;
  MACROBLOCKD *const xd = &x->e_mbd;
  MB_MODE_INFO *mbmi;
  struct macroblock_plane *const p = x->plane;
  struct macroblockd_plane *const pd = xd->plane;
  const AQ_MODE aq_mode = cpi->oxcf.aq_mode;
  int i, orig_rdmult;
  double rdmult_ratio;
#if CONFIG_TX_SKIP
  int q_idx;
#endif

  vp9_clear_system_state();
  rdmult_ratio = 1.0;  // avoid uninitialized warnings

  // Use the lower precision, but faster, 32x32 fdct for mode selection.
  x->use_lp32x32fdct = 1;

  set_offsets(cpi, tile, mi_row, mi_col, bsize);
  mbmi = &xd->mi[0].src_mi->mbmi;
  mbmi->sb_type = bsize;
#if CONFIG_TX_SKIP
  q_idx = vp9_get_qindex(&cm->seg, mbmi->segment_id, cm->base_qindex);
  mbmi->tx_skip_shift = q_idx > TX_SKIP_SHIFT_THRESH ?
      TX_SKIP_SHIFT_HQ : TX_SKIP_SHIFT_LQ;
#endif

  for (i = 0; i < MAX_MB_PLANE; ++i) {
    p[i].coeff = ctx->coeff_pbuf[i][0];
    p[i].qcoeff = ctx->qcoeff_pbuf[i][0];
    pd[i].dqcoeff = ctx->dqcoeff_pbuf[i][0];
    p[i].eobs = ctx->eobs_pbuf[i][0];
  }
#if CONFIG_PALETTE
  for (i = 0; i < 2; ++i) {
    pd[i].color_index_map = ctx->color_index_map[i];
  }
#endif
  ctx->is_coded = 0;
  ctx->skippable = 0;
  ctx->pred_pixel_ready = 0;
  x->skip_recode = 0;

  // Set to zero to make sure we do not use the previous encoded frame stats
  mbmi->skip = 0;
#if CONFIG_SR_MODE
  mbmi->sr = 0;
#endif  // CONFIG_SR_MODE

#if CONFIG_VP9_HIGHBITDEPTH
  if (xd->cur_buf->flags & YV12_FLAG_HIGHBITDEPTH) {
    x->source_variance =
        high_get_sby_perpixel_variance(cpi, &x->plane[0].src, bsize, xd->bd);
  } else {
    x->source_variance =
        get_sby_perpixel_variance(cpi, &x->plane[0].src, bsize);
  }
#else
  x->source_variance = get_sby_perpixel_variance(cpi, &x->plane[0].src, bsize);
#endif  // CONFIG_VP9_HIGHBITDEPTH

  // Save rdmult before it might be changed, so it can be restored later.
  orig_rdmult = x->rdmult;

  if (aq_mode == VARIANCE_AQ) {
    const int energy = bsize <= BLOCK_16X16 ? x->mb_energy
                                            : vp9_block_energy(cpi, x, bsize);
    if (cm->frame_type == KEY_FRAME ||
        cpi->refresh_alt_ref_frame ||
        (cpi->refresh_golden_frame && !cpi->rc.is_src_frame_alt_ref)) {
      mbmi->segment_id = vp9_vaq_segment_id(energy);
    } else {
      const uint8_t *const map = cm->seg.update_map ? cpi->segmentation_map
                                                    : cm->last_frame_seg_map;
      mbmi->segment_id = vp9_get_segment_id(cm, map, bsize, mi_row, mi_col);
    }

    rdmult_ratio = vp9_vaq_rdmult_ratio(energy);
    vp9_init_plane_quantizers(cpi, x);
    vp9_clear_system_state();
    x->rdmult = (int)round(x->rdmult * rdmult_ratio);
  } else if (aq_mode == COMPLEXITY_AQ) {
    const int mi_offset = mi_row * cm->mi_cols + mi_col;
    unsigned char complexity = cpi->complexity_map[mi_offset];
    const int is_edge = (mi_row <= 1) || (mi_row >= (cm->mi_rows - 2)) ||
                        (mi_col <= 1) || (mi_col >= (cm->mi_cols - 2));
    if (!is_edge && (complexity > 128))
      x->rdmult += ((x->rdmult * (complexity - 128)) / 256);
  } else if (aq_mode == CYCLIC_REFRESH_AQ) {
    const uint8_t *const map = cm->seg.update_map ? cpi->segmentation_map
                                                  : cm->last_frame_seg_map;
    // If segment 1, use rdmult for that segment.
    if (vp9_get_segment_id(cm, map, bsize, mi_row, mi_col))
      x->rdmult = vp9_cyclic_refresh_get_rdmult(cpi->cyclic_refresh);
  }

  // Find best coding mode & reconstruct the MB so it is available
  // as a predictor for MBs that follow in the SB
  if (frame_is_intra_only(cm)) {
#if CONFIG_PALETTE
    int n = cpi->common.current_palette_size;
    int count[PALETTE_BUF_SIZE];
#if CONFIG_VP9_HIGHBITDEPTH
    uint16_t palette[PALETTE_BUF_SIZE];
#else
    uint8_t palette[PALETTE_BUF_SIZE];
#endif  // CONFIG_VP9_HIGHBITDEPTH

    vpx_memcpy(palette, cpi->common.current_palette_colors,
               n * sizeof(palette[0]));
    vpx_memcpy(count, cpi->common.current_palette_count,
               n * sizeof(count[0]));
    cpi->common.current_palette_size = ctx->palette_buf_size;
    vpx_memcpy(cpi->common.current_palette_colors, ctx->palette_colors_buf,
               ctx->palette_buf_size * sizeof(ctx->palette_colors_buf[0]));
    vpx_memcpy(cpi->common.current_palette_count, ctx->palette_count_buf,
               ctx->palette_buf_size * sizeof(ctx->palette_count_buf[0]));
#endif
    vp9_rd_pick_intra_mode_sb(cpi, x,
#if CONFIG_INTRABC
                              tile, mi_row, mi_col,
#endif  // CONFIG_INTRABC
                              rd_cost, bsize, ctx, best_rd);
#if CONFIG_PALETTE
    cpi->common.current_palette_size = n;
    vpx_memcpy(cpi->common.current_palette_colors,
               palette, n * sizeof(palette[0]));
    vpx_memcpy(cpi->common.current_palette_count,
               count, n * sizeof(count[0]));
#endif
#if CONFIG_SUPERTX
    *totalrate_nocoef = 0;
#endif
  } else {
    if (bsize >= BLOCK_8X8) {
      if (vp9_segfeature_active(&cm->seg, mbmi->segment_id, SEG_LVL_SKIP)) {
        vp9_rd_pick_inter_mode_sb_seg_skip(cpi, x, rd_cost, bsize,
                                           ctx, best_rd);
#if CONFIG_SUPERTX
        *totalrate_nocoef = rd_cost->rate;
#endif
      } else {
        vp9_rd_pick_inter_mode_sb(cpi, x, tile, mi_row, mi_col, rd_cost,
#if CONFIG_SUPERTX
                                  totalrate_nocoef,
#endif
#if CONFIG_COPY_MODE
#if CONFIG_EXT_PARTITION
                                  ctx->partition,
#endif
#endif
                                  bsize, ctx, best_rd);
      }
    } else {
      vp9_rd_pick_inter_mode_sub8x8(cpi, x, tile, mi_row, mi_col, rd_cost,
#if CONFIG_SUPERTX
                                    totalrate_nocoef,
#endif
                                    bsize, ctx, best_rd);
    }
  }

  if (aq_mode == VARIANCE_AQ && rd_cost->rate != INT_MAX) {
    vp9_clear_system_state();
    rd_cost->rate = (int)round(rd_cost->rate * rdmult_ratio);
    rd_cost->rdcost = RDCOST(x->rdmult, x->rddiv, rd_cost->rate, rd_cost->dist);
#if CONFIG_SUPERTX
    *totalrate_nocoef = (int)round(*totalrate_nocoef * rdmult_ratio);
#endif
  }

  x->rdmult = orig_rdmult;

  // TODO(jingning) The rate-distortion optimization flow needs to be
  // refactored to provide proper exit/return handle.
  if (rd_cost->rate == INT_MAX)
    rd_cost->rdcost = INT64_MAX;
}

static void update_stats(VP9_COMMON *cm, const MACROBLOCK *x) {
  const MACROBLOCKD *const xd = &x->e_mbd;
  const MODE_INFO *const mi = xd->mi[0].src_mi;
  const MB_MODE_INFO *const mbmi = &mi->mbmi;
  const BLOCK_SIZE bsize = mbmi->sb_type;

#if CONFIG_COPY_MODE
  if (!frame_is_intra_only(cm)) {
    COPY_MODE copy_mode = mbmi->copy_mode;
    if (mbmi->sb_type >= BLOCK_8X8) {
      int copy_mode_context = vp9_get_copy_mode_context(xd);
      if (mbmi->inter_ref_count > 0) {
        ++cm->counts.copy_noref[copy_mode_context][mbmi->sb_type]
                               [copy_mode != NOREF];
        if (copy_mode != NOREF) {
          if (mbmi->inter_ref_count == 2)
            ++cm->counts.copy_mode_l2[copy_mode_context][copy_mode - REF0];
          else if (mbmi->inter_ref_count > 2)
            ++cm->counts.copy_mode[copy_mode_context][copy_mode - REF0];
        }
      }
    }
  }
  if (!frame_is_intra_only(cm) && mbmi->copy_mode == NOREF) {
#else
  if (!frame_is_intra_only(cm)) {
#endif  // CONFIG_COPY_MODE
    FRAME_COUNTS *const counts = &cm->counts;
    const int inter_block = is_inter_block(mbmi);
    const int seg_ref_active = vp9_segfeature_active(&cm->seg, mbmi->segment_id,
                                                     SEG_LVL_REF_FRAME);
    if (!seg_ref_active) {
      counts->intra_inter[vp9_get_intra_inter_context(xd)][inter_block]++;

      // If the segment reference feature is enabled we have only a single
      // reference frame allowed for the segment so exclude it from
      // the reference frame counts used to work out probabilities.
      if (inter_block) {
        const MV_REFERENCE_FRAME ref0 = mbmi->ref_frame[0];

        if (cm->reference_mode == REFERENCE_MODE_SELECT)
          counts->comp_inter[vp9_get_reference_mode_context(cm, xd)]
                            [has_second_ref(mbmi)]++;

        if (has_second_ref(mbmi)) {
#if CONFIG_MULTI_REF
          const int bit = (ref0 == GOLDEN_FRAME || ref0 == LAST3_FRAME ||
                           ref0 == LAST4_FRAME);
          counts->comp_ref[vp9_get_pred_context_comp_ref_p(cm, xd)][0][bit]++;
          if (!bit) {
            counts->comp_ref[vp9_get_pred_context_comp_ref_p1(cm, xd)][1]
                            [ref0 == LAST_FRAME]++;
          } else {
            counts->comp_ref[vp9_get_pred_context_comp_ref_p2(cm, xd)][2]
                            [ref0 == GOLDEN_FRAME]++;
            if (ref0 != GOLDEN_FRAME) {
              counts->comp_ref[vp9_get_pred_context_comp_ref_p3(cm, xd)][3]
                              [ref0 == LAST3_FRAME]++;
            }
          }
#else  // CONFIG_MULTI_REF
          counts->comp_ref[vp9_get_pred_context_comp_ref_p(cm, xd)][0]
                          [ref0 == GOLDEN_FRAME]++;
#endif  // CONFIG_MULTI_REF
        } else {
#if CONFIG_MULTI_REF
          const int bit = (ref0 == ALTREF_FRAME || ref0 == GOLDEN_FRAME);
          counts->single_ref[vp9_get_pred_context_single_ref_p1(xd)][0][bit]++;
          if (bit) {
            counts->single_ref[vp9_get_pred_context_single_ref_p2(xd)][1]
                              [ref0 != GOLDEN_FRAME]++;
          } else {
            const int bit1 = !(ref0 == LAST2_FRAME || ref0 == LAST_FRAME);
            counts->single_ref[vp9_get_pred_context_single_ref_p3(xd)][2]
                              [bit1]++;
            if (!bit1) {
              counts->single_ref[vp9_get_pred_context_single_ref_p4(xd)][3]
                                [ref0 != LAST_FRAME]++;
            } else {
              counts->single_ref[vp9_get_pred_context_single_ref_p5(xd)][4]
                                [ref0 != LAST3_FRAME]++;
            }
          }
#else  // CONFIG_MULTI_REF
          counts->single_ref[vp9_get_pred_context_single_ref_p1(xd)][0]
                            [ref0 != LAST_FRAME]++;
          if (ref0 != LAST_FRAME)
            counts->single_ref[vp9_get_pred_context_single_ref_p2(xd)][1]
                              [ref0 != GOLDEN_FRAME]++;
#endif  // CONFIG_MULTI_REF
        }
      }
    }
    if (inter_block) {
      vp9_update_mv_count(cm, xd);

      if (cm->interp_filter == SWITCHABLE) {
        const int ctx = vp9_get_pred_context_switchable_interp(xd);
        ++cm->counts.switchable_interp[ctx][mbmi->interp_filter];
      }
#if CONFIG_INTERINTRA
      if (is_interintra_allowed(bsize) &&
          is_inter_mode(mbmi->mode) &&
#if CONFIG_SUPERTX
          mbmi->tx_size <= max_txsize_lookup[bsize] &&
#endif  // CONFIG_SUPERTX
          (mbmi->ref_frame[1] <= INTRA_FRAME)) {
        if (mbmi->ref_frame[1] == INTRA_FRAME) {
          ++cm->counts.y_mode[size_group_lookup[bsize]][mbmi->interintra_mode];
          ++cm->counts.interintra[bsize][1];
#if CONFIG_WEDGE_PARTITION
          if (get_wedge_bits(bsize))
            ++cm->counts.wedge_interintra[bsize][mbmi->use_wedge_interintra];
#endif  // CONFIG_WEDGE_PARTITION
        } else {
          ++cm->counts.interintra[bsize][0];
        }
      }
#endif  // CONFIG_INTERINTRA
#if CONFIG_WEDGE_PARTITION
      if (cm->reference_mode != SINGLE_REFERENCE &&
          get_wedge_bits(bsize) &&
          mbmi->ref_frame[1] > INTRA_FRAME) {
        ++cm->counts.wedge_interinter[bsize][mbmi->use_wedge_interinter];
      }
#endif  // CONFIG_WEDGE_PARTITION
    }

    if (inter_block &&
        !vp9_segfeature_active(&cm->seg, mbmi->segment_id, SEG_LVL_SKIP)) {
      const int mode_ctx = mbmi->mode_context[mbmi->ref_frame[0]];
      if (bsize >= BLOCK_8X8) {
        const PREDICTION_MODE mode = mbmi->mode;
#if CONFIG_NEW_INTER
        if (is_inter_compound_mode(mode))
          ++counts->inter_compound_mode[mode_ctx][INTER_COMPOUND_OFFSET(mode)];
        else
#endif  // CONFIG_NEW_INTER
        ++counts->inter_mode[mode_ctx][INTER_OFFSET(mode)];
      } else {
        const int num_4x4_w = num_4x4_blocks_wide_lookup[bsize];
        const int num_4x4_h = num_4x4_blocks_high_lookup[bsize];
        int idx, idy;
        for (idy = 0; idy < 2; idy += num_4x4_h) {
          for (idx = 0; idx < 2; idx += num_4x4_w) {
            const int j = idy * 2 + idx;
            const PREDICTION_MODE b_mode = mi->bmi[j].as_mode;
#if CONFIG_NEW_INTER
            if (is_inter_compound_mode(b_mode))
              ++counts->inter_compound_mode[mode_ctx]
                                           [INTER_COMPOUND_OFFSET(b_mode)];
            else
#endif  // CONFIG_NEW_INTER
            ++counts->inter_mode[mode_ctx][INTER_OFFSET(b_mode)];
          }
        }
      }
    }
  }
}

static void restore_context(VP9_COMP *cpi, int mi_row, int mi_col,
                            ENTROPY_CONTEXT a[(CODING_UNIT_SIZE >> 2) *
                                              MAX_MB_PLANE],
                            ENTROPY_CONTEXT l[(CODING_UNIT_SIZE >> 2) *
                                              MAX_MB_PLANE],
                            PARTITION_CONTEXT sa[CODING_UNIT_SIZE >> 3],
                            PARTITION_CONTEXT sl[CODING_UNIT_SIZE >> 3],
                            BLOCK_SIZE bsize) {
  MACROBLOCK *const x = &cpi->mb;
  MACROBLOCKD *const xd = &x->e_mbd;
  int p;
  const int num_4x4_blocks_wide = num_4x4_blocks_wide_lookup[bsize];
  const int num_4x4_blocks_high = num_4x4_blocks_high_lookup[bsize];
  int mi_width = num_8x8_blocks_wide_lookup[bsize];
  int mi_height = num_8x8_blocks_high_lookup[bsize];
  for (p = 0; p < MAX_MB_PLANE; p++) {
    vpx_memcpy(
        xd->above_context[p] + ((mi_col * 2) >> xd->plane[p].subsampling_x),
        a + num_4x4_blocks_wide * p,
        (sizeof(ENTROPY_CONTEXT) * num_4x4_blocks_wide) >>
        xd->plane[p].subsampling_x);
    vpx_memcpy(
        xd->left_context[p]
            + ((mi_row & MI_MASK) * 2 >> xd->plane[p].subsampling_y),
        l + num_4x4_blocks_high * p,
        (sizeof(ENTROPY_CONTEXT) * num_4x4_blocks_high) >>
        xd->plane[p].subsampling_y);
  }
  vpx_memcpy(xd->above_seg_context + mi_col, sa,
             sizeof(*xd->above_seg_context) * mi_width);
  vpx_memcpy(xd->left_seg_context + (mi_row & MI_MASK), sl,
             sizeof(xd->left_seg_context[0]) * mi_height);
}

static void save_context(VP9_COMP *cpi, int mi_row, int mi_col,
                            ENTROPY_CONTEXT a[(CODING_UNIT_SIZE >> 2) *
                                              MAX_MB_PLANE],
                            ENTROPY_CONTEXT l[(CODING_UNIT_SIZE >> 2) *
                                              MAX_MB_PLANE],
                            PARTITION_CONTEXT sa[CODING_UNIT_SIZE >> 3],
                            PARTITION_CONTEXT sl[CODING_UNIT_SIZE >> 3],
                         BLOCK_SIZE bsize) {
  const MACROBLOCK *const x = &cpi->mb;
  const MACROBLOCKD *const xd = &x->e_mbd;
  int p;
  const int num_4x4_blocks_wide = num_4x4_blocks_wide_lookup[bsize];
  const int num_4x4_blocks_high = num_4x4_blocks_high_lookup[bsize];
  int mi_width = num_8x8_blocks_wide_lookup[bsize];
  int mi_height = num_8x8_blocks_high_lookup[bsize];

  // buffer the above/left context information of the block in search.
  for (p = 0; p < MAX_MB_PLANE; ++p) {
    vpx_memcpy(
        a + num_4x4_blocks_wide * p,
        xd->above_context[p] + (mi_col * 2 >> xd->plane[p].subsampling_x),
        (sizeof(ENTROPY_CONTEXT) * num_4x4_blocks_wide) >>
        xd->plane[p].subsampling_x);
    vpx_memcpy(
        l + num_4x4_blocks_high * p,
        xd->left_context[p]
            + ((mi_row & MI_MASK) * 2 >> xd->plane[p].subsampling_y),
        (sizeof(ENTROPY_CONTEXT) * num_4x4_blocks_high) >>
        xd->plane[p].subsampling_y);
  }
  vpx_memcpy(sa, xd->above_seg_context + mi_col,
             sizeof(*xd->above_seg_context) * mi_width);
  vpx_memcpy(sl, xd->left_seg_context + (mi_row & MI_MASK),
             sizeof(xd->left_seg_context[0]) * mi_height);
}

static void encode_b(VP9_COMP *cpi, const TileInfo *const tile,
                     TOKENEXTRA **tp, int mi_row, int mi_col,
                     int output_enabled, BLOCK_SIZE bsize,
                     PICK_MODE_CONTEXT *ctx) {
  set_offsets(cpi, tile, mi_row, mi_col, bsize);
  update_state(cpi, ctx, mi_row, mi_col, bsize, output_enabled);
  encode_superblock(cpi, tp, output_enabled, mi_row, mi_col, bsize, ctx);

  if (output_enabled) {
    update_stats(&cpi->common, &cpi->mb);

    (*tp)->token = EOSB_TOKEN;
    (*tp)++;
  }
}

static void encode_sb(VP9_COMP *cpi, const TileInfo *const tile,
                      TOKENEXTRA **tp, int mi_row, int mi_col,
                      int output_enabled, BLOCK_SIZE bsize,
                      PC_TREE *pc_tree) {
  VP9_COMMON *const cm = &cpi->common;
  MACROBLOCK *const x = &cpi->mb;
  MACROBLOCKD *const xd = &x->e_mbd;

  const int bsl = b_width_log2_lookup[bsize], hbs = (1 << bsl) / 4;
  int ctx;
  PARTITION_TYPE partition;
  BLOCK_SIZE subsize = bsize;
#if CONFIG_EXT_PARTITION
  BLOCK_SIZE bsize2 = get_subsize(bsize, PARTITION_SPLIT);
#endif

  if (mi_row >= cm->mi_rows || mi_col >= cm->mi_cols)
    return;

  if (bsize >= BLOCK_8X8) {
    ctx = partition_plane_context(xd, mi_row, mi_col, bsize);
    subsize = get_subsize(bsize, pc_tree->partitioning);
  } else {
    ctx = 0;
    subsize = BLOCK_4X4;
  }

  partition = partition_lookup[bsl][subsize];
#if CONFIG_EXT_PARTITION
  if (bsize > BLOCK_8X8)
    partition = pc_tree->partitioning;
#endif
  if (output_enabled && bsize != BLOCK_4X4)
    cm->counts.partition[ctx][partition]++;

#if CONFIG_SUPERTX
  if (cm->frame_type != KEY_FRAME &&
      bsize <= MAX_SUPERTX_BLOCK_SIZE &&
      partition != PARTITION_NONE && !xd->lossless) {
    int supertx_enabled;
    TX_SIZE supertx_size = bsize_to_tx_size(bsize);
    supertx_enabled = check_supertx_sb(bsize, supertx_size, pc_tree);
    if (supertx_enabled) {
      const int mi_width = num_8x8_blocks_wide_lookup[bsize];
      const int mi_height = num_8x8_blocks_high_lookup[bsize];
      int x_idx, y_idx, i;
      uint8_t *dst_buf[3];
      int dst_stride[3];
      set_skip_context(xd, mi_row, mi_col);
      set_modeinfo_offsets(cm, xd, mi_row, mi_col);
      update_state_sb_supertx(cpi, tile, mi_row, mi_col, bsize,
                              output_enabled, pc_tree);

      vp9_setup_dst_planes(xd->plane, get_frame_new_buffer(cm),
                           mi_row, mi_col);
      for (i = 0; i < MAX_MB_PLANE; i++) {
        dst_buf[i] = xd->plane[i].dst.buf;
        dst_stride[i] = xd->plane[i].dst.stride;
      }
      predict_sb_complex(cpi, tile, mi_row, mi_col, mi_row, mi_col,
                         output_enabled, bsize, bsize,
                         dst_buf, dst_stride, pc_tree);

      set_offsets(cpi, tile, mi_row, mi_col, bsize);
      if (!x->skip) {
        xd->mi[0].mbmi.skip = 1;
        vp9_encode_sb_supertx(x, bsize);
        vp9_tokenize_sb_supertx(cpi, tp, !output_enabled, bsize);
      } else {
        xd->mi[0].mbmi.skip = 1;
        if (output_enabled)
          cm->counts.skip[vp9_get_skip_context(xd)][1]++;
        reset_skip_context(xd, bsize);
      }
      if (output_enabled) {
        for (y_idx = 0; y_idx < mi_height; y_idx++)
          for (x_idx = 0; x_idx < mi_width; x_idx++) {
            if ((xd->mb_to_right_edge >> (3 + MI_SIZE_LOG2)) + mi_width > x_idx
                && (xd->mb_to_bottom_edge >> (3 + MI_SIZE_LOG2)) + mi_height
                    > y_idx) {
              xd->mi[x_idx + y_idx * cm->mi_stride].mbmi.skip =
                  xd->mi[0].mbmi.skip;
            }
          }
        cm->counts.supertx
            [partition_supertx_context_lookup[partition]][supertx_size][1]++;
        cm->counts.supertx_size[supertx_size]++;
#if CONFIG_NEW_QUANT && QUANT_PROFILES > 1 && !Q_CTX_BASED_PROFILES
        if (cm->base_qindex > Q_THRESHOLD_MIN &&
            cm->base_qindex < Q_THRESHOLD_MAX &&
            !xd->mi[0].mbmi.skip &&
            xd->mi[0].mbmi.send_dq_bit &&
            !vp9_segfeature_active(&cm->seg, xd->mi[0].mbmi.segment_id,
                                   SEG_LVL_SKIP)) {
          ++cm->counts.dq_profile[xd->mi[0].mbmi.dq_off_index];
        }
#endif  // CONFIG_NEW_QUANT && QUANT_PROFILES > 1 && !Q_CTX_BASED_PROFILES
#if CONFIG_EXT_TX
#if CONFIG_WAVELETS
        if (!xd->mi[0].mbmi.skip)
          ++cm->counts.ext_tx[xd->mi[0].mbmi.tx_size]
                             [xd->mi[0].mbmi.ext_txfrm];
#else
        if (supertx_size <= TX_16X16 && !xd->mi[0].mbmi.skip)
          ++cm->counts.ext_tx[xd->mi[0].mbmi.tx_size]
                             [xd->mi[0].mbmi.ext_txfrm];
#endif  // CONFIG_WAVELETS
#endif  // CONFIG_EXT_TX
        (*tp)->token = EOSB_TOKEN;
        (*tp)++;
      }
#if CONFIG_EXT_PARTITION
      update_ext_partition_context(xd, mi_row, mi_col, subsize, bsize,
                                   partition);
#else
      if (partition != PARTITION_SPLIT || bsize == BLOCK_8X8)
        update_partition_context(xd, mi_row, mi_col, subsize, bsize);
#endif
      return;
    } else {
      if (output_enabled) {
        cm->counts.supertx
            [partition_supertx_context_lookup[partition]][supertx_size][0]++;
      }
    }
  }
#endif  // CONFIG_SUPERTX

  switch (partition) {
    case PARTITION_NONE:
      encode_b(cpi, tile, tp, mi_row, mi_col, output_enabled, subsize,
               &pc_tree->none);
      break;
    case PARTITION_VERT:
      encode_b(cpi, tile, tp, mi_row, mi_col, output_enabled, subsize,
               &pc_tree->vertical[0]);
      if (mi_col + hbs < cm->mi_cols && bsize > BLOCK_8X8) {
        encode_b(cpi, tile, tp, mi_row, mi_col + hbs, output_enabled, subsize,
                 &pc_tree->vertical[1]);
      }
      break;
    case PARTITION_HORZ:
      encode_b(cpi, tile, tp, mi_row, mi_col, output_enabled, subsize,
               &pc_tree->horizontal[0]);
      if (mi_row + hbs < cm->mi_rows && bsize > BLOCK_8X8) {
        encode_b(cpi, tile, tp, mi_row + hbs, mi_col, output_enabled, subsize,
                 &pc_tree->horizontal[1]);
      }
      break;
    case PARTITION_SPLIT:
      if (bsize == BLOCK_8X8) {
        encode_b(cpi, tile, tp, mi_row, mi_col, output_enabled, subsize,
                 pc_tree->leaf_split[0]);
      } else {
        encode_sb(cpi, tile, tp, mi_row, mi_col, output_enabled, subsize,
                  pc_tree->split[0]);
        encode_sb(cpi, tile, tp, mi_row, mi_col + hbs, output_enabled, subsize,
                  pc_tree->split[1]);
        encode_sb(cpi, tile, tp, mi_row + hbs, mi_col, output_enabled, subsize,
                  pc_tree->split[2]);
        encode_sb(cpi, tile, tp, mi_row + hbs, mi_col + hbs, output_enabled,
                  subsize, pc_tree->split[3]);
      }
      break;
#if CONFIG_EXT_PARTITION
    case PARTITION_HORZ_A:
      encode_b(cpi, tile, tp, mi_row, mi_col, output_enabled, bsize2,
               &pc_tree->horizontala[0]);
      encode_b(cpi, tile, tp, mi_row, mi_col + hbs, output_enabled, bsize2,
               &pc_tree->horizontala[1]);
      encode_b(cpi, tile, tp, mi_row + hbs, mi_col, output_enabled, subsize,
               &pc_tree->horizontala[2]);
      break;
    case PARTITION_HORZ_B:
      encode_b(cpi, tile, tp, mi_row, mi_col, output_enabled, subsize,
               &pc_tree->horizontalb[0]);
      encode_b(cpi, tile, tp, mi_row + hbs, mi_col, output_enabled, bsize2,
               &pc_tree->horizontalb[1]);
      encode_b(cpi, tile, tp, mi_row + hbs, mi_col + hbs, output_enabled,
               bsize2, &pc_tree->horizontalb[2]);
      break;
    case PARTITION_VERT_A:
      encode_b(cpi, tile, tp, mi_row, mi_col, output_enabled, bsize2,
               &pc_tree->verticala[0]);
      encode_b(cpi, tile, tp, mi_row + hbs, mi_col, output_enabled, bsize2,
               &pc_tree->verticala[1]);
      encode_b(cpi, tile, tp, mi_row, mi_col + hbs, output_enabled, subsize,
               &pc_tree->verticala[2]);

      break;
    case PARTITION_VERT_B:
      encode_b(cpi, tile, tp, mi_row, mi_col, output_enabled, subsize,
               &pc_tree->verticalb[0]);
      encode_b(cpi, tile, tp, mi_row, mi_col + hbs, output_enabled, bsize2,
               &pc_tree->verticalb[1]);
      encode_b(cpi, tile, tp, mi_row + hbs, mi_col + hbs, output_enabled,
               bsize2, &pc_tree->verticalb[2]);
      break;
#endif
    default:
      assert(0 && "Invalid partition type.");
      break;
  }
#if CONFIG_EXT_PARTITION
  update_ext_partition_context(xd, mi_row, mi_col, subsize, bsize, partition);
#else
  if (partition != PARTITION_SPLIT || bsize == BLOCK_8X8)
    update_partition_context(xd, mi_row, mi_col, subsize, bsize);
#endif
}

// Check to see if the given partition size is allowed for a specified number
// of 8x8 block rows and columns remaining in the image.
// If not then return the largest allowed partition size
static BLOCK_SIZE find_partition_size(BLOCK_SIZE bsize,
                                      int rows_left, int cols_left,
                                      int *bh, int *bw) {
  if (rows_left <= 0 || cols_left <= 0) {
    return MIN(bsize, BLOCK_8X8);
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

static void set_partial_b64x64_partition(MODE_INFO *mi, int mis,
    int bh_in, int bw_in, int row8x8_remaining, int col8x8_remaining,
    BLOCK_SIZE bsize, MODE_INFO *mi_8x8) {
  int bh = bh_in;
  int r, c;
  for (r = 0; r < MI_BLOCK_SIZE; r += bh) {
    int bw = bw_in;
    for (c = 0; c < MI_BLOCK_SIZE; c += bw) {
      const int index = r * mis + c;
      mi_8x8[index].src_mi = mi + index;
      mi_8x8[index].src_mi->mbmi.sb_type = find_partition_size(bsize,
          row8x8_remaining - r, col8x8_remaining - c, &bh, &bw);
    }
  }
}

// This function attempts to set all mode info entries in a given SB64
// to the same block partition size.
// However, at the bottom and right borders of the image the requested size
// may not be allowed in which case this code attempts to choose the largest
// allowable partition.
static void set_fixed_partitioning(VP9_COMP *cpi, const TileInfo *const tile,
                                   MODE_INFO *mi_8x8, int mi_row, int mi_col,
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
        mi_8x8[index].src_mi = mi_upper_left + index;
        mi_8x8[index].src_mi->mbmi.sb_type = bsize;
      }
    }
  } else {
    // Else this is a partial SB64.
    set_partial_b64x64_partition(mi_upper_left, mis, bh, bw, row8x8_remaining,
        col8x8_remaining, bsize, mi_8x8);
  }
}

static void rd_use_partition(VP9_COMP *cpi, const TileInfo *const tile,
                             MODE_INFO *mi_8x8, TOKENEXTRA **tp,
                             int mi_row, int mi_col,
                             BLOCK_SIZE bsize, int *rate,
                             int64_t *dist,
#if CONFIG_SUPERTX
                             int *rate_nocoef,
#endif
                             int do_recon, PC_TREE *pc_tree) {
  VP9_COMMON *const cm = &cpi->common;
  MACROBLOCK *const x = &cpi->mb;
  MACROBLOCKD *const xd = &x->e_mbd;
  const int mis = cm->mi_stride;
  const int bsl = b_width_log2_lookup[bsize];
  const int mi_step = num_4x4_blocks_wide_lookup[bsize] / 2;
  const int bss = (1 << bsl) / 4;
  int i, pl;
  PARTITION_TYPE partition = PARTITION_NONE;
  BLOCK_SIZE subsize;
  ENTROPY_CONTEXT l[(CODING_UNIT_SIZE >> 2) * MAX_MB_PLANE];
  ENTROPY_CONTEXT a[(CODING_UNIT_SIZE >> 2) * MAX_MB_PLANE];
  PARTITION_CONTEXT sl[CODING_UNIT_SIZE >> 3], sa[CODING_UNIT_SIZE >> 3];
  RD_COST last_part_rdc, none_rdc, chosen_rdc;
  BLOCK_SIZE sub_subsize = BLOCK_4X4;
  int splits_below = 0;
  BLOCK_SIZE bs_type = mi_8x8[0].src_mi->mbmi.sb_type;
  int do_partition_search = 1;
  PICK_MODE_CONTEXT *ctx = &pc_tree->none;
#if CONFIG_SUPERTX
  int last_part_rate_nocoef = INT_MAX;
  int none_rate_nocoef = INT_MAX;
  int chosen_rate_nocoef = INT_MAX;
#endif

  if (mi_row >= cm->mi_rows || mi_col >= cm->mi_cols)
    return;

  assert(num_4x4_blocks_wide_lookup[bsize] ==
         num_4x4_blocks_high_lookup[bsize]);

  vp9_rd_cost_reset(&last_part_rdc);
  vp9_rd_cost_reset(&none_rdc);
  vp9_rd_cost_reset(&chosen_rdc);

  partition = partition_lookup[bsl][bs_type];
  subsize = get_subsize(bsize, partition);

  pc_tree->partitioning = partition;
  save_context(cpi, mi_row, mi_col, a, l, sa, sl, bsize);

  if (bsize == BLOCK_16X16 && cpi->oxcf.aq_mode) {
    set_offsets(cpi, tile, mi_row, mi_col, bsize);
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
        MODE_INFO *this_mi = mi_8x8[jj * bss * mis + ii * bss].src_mi;
        if (this_mi && this_mi->mbmi.sb_type >= sub_subsize) {
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
      rd_pick_sb_modes(cpi, tile, mi_row, mi_col, &none_rdc,
#if CONFIG_SUPERTX
                       &none_rate_nocoef,
#endif
                       bsize, ctx, INT64_MAX);

      pl = partition_plane_context(xd, mi_row, mi_col, bsize);

      if (none_rdc.rate < INT_MAX) {
        none_rdc.rate += cpi->partition_cost[pl][PARTITION_NONE];
        none_rdc.rdcost = RDCOST(x->rdmult, x->rddiv, none_rdc.rate,
                                 none_rdc.dist);
#if CONFIG_SUPERTX
        none_rate_nocoef += cpi->partition_cost[pl][PARTITION_NONE];
#endif
      }

      restore_context(cpi, mi_row, mi_col, a, l, sa, sl, bsize);
      mi_8x8[0].src_mi->mbmi.sb_type = bs_type;
      pc_tree->partitioning = partition;
    }
  }

  switch (partition) {
    case PARTITION_NONE:
      rd_pick_sb_modes(cpi, tile, mi_row, mi_col, &last_part_rdc,
#if CONFIG_SUPERTX
                       &last_part_rate_nocoef,
#endif
                       bsize, ctx, INT64_MAX);
      break;
    case PARTITION_HORZ:
      rd_pick_sb_modes(cpi, tile, mi_row, mi_col, &last_part_rdc,
#if CONFIG_SUPERTX
                       &last_part_rate_nocoef,
#endif
                       subsize, &pc_tree->horizontal[0],
                       INT64_MAX);
      if (last_part_rdc.rate != INT_MAX &&
          bsize >= BLOCK_8X8 && mi_row + (mi_step >> 1) < cm->mi_rows) {
        RD_COST tmp_rdc;
#if CONFIG_SUPERTX
        int rt_nocoef = 0;
#endif
        PICK_MODE_CONTEXT *ctx = &pc_tree->horizontal[0];
        vp9_rd_cost_init(&tmp_rdc);
        update_state(cpi, ctx, mi_row, mi_col, subsize, 0);
        encode_superblock(cpi, tp, 0, mi_row, mi_col, subsize, ctx);
        rd_pick_sb_modes(cpi, tile, mi_row + (mi_step >> 1), mi_col, &tmp_rdc,
#if CONFIG_SUPERTX
                         &rt_nocoef,
#endif
                         subsize, &pc_tree->horizontal[1], INT64_MAX);
        if (tmp_rdc.rate == INT_MAX || tmp_rdc.dist == INT64_MAX) {
          vp9_rd_cost_reset(&last_part_rdc);
#if CONFIG_SUPERTX
          last_part_rate_nocoef = INT_MAX;
#endif
          break;
        }
        last_part_rdc.rate += tmp_rdc.rate;
        last_part_rdc.dist += tmp_rdc.dist;
        last_part_rdc.rdcost += tmp_rdc.rdcost;
#if CONFIG_SUPERTX
        last_part_rate_nocoef += rt_nocoef;
#endif
      }
      break;
    case PARTITION_VERT:
      rd_pick_sb_modes(cpi, tile, mi_row, mi_col, &last_part_rdc,
#if CONFIG_SUPERTX
                       &last_part_rate_nocoef,
#endif
                       subsize, &pc_tree->vertical[0], INT64_MAX);
      if (last_part_rdc.rate != INT_MAX &&
          bsize >= BLOCK_8X8 && mi_col + (mi_step >> 1) < cm->mi_cols) {
        RD_COST tmp_rdc;
#if CONFIG_SUPERTX
        int rt_nocoef = 0;
#endif
        PICK_MODE_CONTEXT *ctx = &pc_tree->vertical[0];
        vp9_rd_cost_init(&tmp_rdc);
        update_state(cpi, ctx, mi_row, mi_col, subsize, 0);
        encode_superblock(cpi, tp, 0, mi_row, mi_col, subsize, ctx);
        rd_pick_sb_modes(cpi, tile, mi_row, mi_col + (mi_step >> 1), &tmp_rdc,
#if CONFIG_SUPERTX
                         &rt_nocoef,
#endif
                         subsize, &pc_tree->vertical[bsize > BLOCK_8X8],
                         INT64_MAX);
        if (tmp_rdc.rate == INT_MAX || tmp_rdc.dist == INT64_MAX) {
          vp9_rd_cost_reset(&last_part_rdc);
#if CONFIG_SUPERTX
          last_part_rate_nocoef = INT_MAX;
#endif
          break;
        }
        last_part_rdc.rate += tmp_rdc.rate;
        last_part_rdc.dist += tmp_rdc.dist;
        last_part_rdc.rdcost += tmp_rdc.rdcost;
#if CONFIG_SUPERTX
        last_part_rate_nocoef += rt_nocoef;
#endif
      }
      break;
    case PARTITION_SPLIT:
      if (bsize == BLOCK_8X8) {
        rd_pick_sb_modes(cpi, tile, mi_row, mi_col, &last_part_rdc,
#if CONFIG_SUPERTX
                         &last_part_rate_nocoef,
#endif
                         subsize, pc_tree->leaf_split[0], INT64_MAX);
        break;
      }
      last_part_rdc.rate = 0;
      last_part_rdc.dist = 0;
      last_part_rdc.rdcost = 0;
#if CONFIG_SUPERTX
      last_part_rate_nocoef = 0;
#endif
      for (i = 0; i < 4; i++) {
        int x_idx = (i & 1) * (mi_step >> 1);
        int y_idx = (i >> 1) * (mi_step >> 1);
        int jj = i >> 1, ii = i & 0x01;
        RD_COST tmp_rdc;
#if CONFIG_SUPERTX
        int rt_nocoef;
#endif
        if ((mi_row + y_idx >= cm->mi_rows) || (mi_col + x_idx >= cm->mi_cols))
          continue;

        vp9_rd_cost_init(&tmp_rdc);
        rd_use_partition(cpi, tile, mi_8x8 + jj * bss * mis + ii * bss, tp,
                         mi_row + y_idx, mi_col + x_idx, subsize,
                         &tmp_rdc.rate, &tmp_rdc.dist,
#if CONFIG_SUPERTX
                         &rt_nocoef,
#endif
                         i != 3, pc_tree->split[i]);
        if (tmp_rdc.rate == INT_MAX || tmp_rdc.dist == INT64_MAX) {
          vp9_rd_cost_reset(&last_part_rdc);
#if CONFIG_SUPERTX
          last_part_rate_nocoef = INT_MAX;
#endif
          break;
        }
        last_part_rdc.rate += tmp_rdc.rate;
        last_part_rdc.dist += tmp_rdc.dist;
#if CONFIG_SUPERTX
        last_part_rate_nocoef += rt_nocoef;
#endif
      }
      break;
    default:
      assert(0);
      break;
  }

  pl = partition_plane_context(xd, mi_row, mi_col, bsize);
  if (last_part_rdc.rate < INT_MAX) {
    last_part_rdc.rate += cpi->partition_cost[pl][partition];
    last_part_rdc.rdcost = RDCOST(x->rdmult, x->rddiv,
                                  last_part_rdc.rate, last_part_rdc.dist);
#if CONFIG_SUPERTX
    last_part_rate_nocoef += cpi->partition_cost[pl][partition];
#endif
  }

  if (do_partition_search
      && cpi->sf.adjust_partitioning_from_last_frame
      && cpi->sf.partition_search_type == SEARCH_PARTITION
      && partition != PARTITION_SPLIT && bsize > BLOCK_8X8
      && (mi_row + mi_step < cm->mi_rows ||
          mi_row + (mi_step >> 1) == cm->mi_rows)
      && (mi_col + mi_step < cm->mi_cols ||
          mi_col + (mi_step >> 1) == cm->mi_cols)) {
    BLOCK_SIZE split_subsize = get_subsize(bsize, PARTITION_SPLIT);
    chosen_rdc.rate = 0;
    chosen_rdc.dist = 0;
#if CONFIG_SUPERTX
    chosen_rate_nocoef = 0;
#endif
    restore_context(cpi, mi_row, mi_col, a, l, sa, sl, bsize);
    pc_tree->partitioning = PARTITION_SPLIT;

    // Split partition.
    for (i = 0; i < 4; i++) {
      int x_idx = (i & 1) * (mi_step >> 1);
      int y_idx = (i >> 1) * (mi_step >> 1);
      RD_COST tmp_rdc;
#if CONFIG_SUPERTX
      int rt_nocoef = 0;
#endif
      ENTROPY_CONTEXT l[(CODING_UNIT_SIZE >> 2) * MAX_MB_PLANE];
      ENTROPY_CONTEXT a[(CODING_UNIT_SIZE >> 2) * MAX_MB_PLANE];
      PARTITION_CONTEXT sl[CODING_UNIT_SIZE >> 3], sa[CODING_UNIT_SIZE >> 3];

      if ((mi_row + y_idx >= cm->mi_rows) || (mi_col + x_idx >= cm->mi_cols))
        continue;

      save_context(cpi, mi_row, mi_col, a, l, sa, sl, bsize);
      pc_tree->split[i]->partitioning = PARTITION_NONE;
      rd_pick_sb_modes(cpi, tile, mi_row + y_idx, mi_col + x_idx, &tmp_rdc,
#if CONFIG_SUPERTX
                       &rt_nocoef,
#endif
                       split_subsize, &pc_tree->split[i]->none, INT64_MAX);

      restore_context(cpi, mi_row, mi_col, a, l, sa, sl, bsize);

      if (tmp_rdc.rate == INT_MAX || tmp_rdc.dist == INT64_MAX) {
        vp9_rd_cost_reset(&chosen_rdc);
#if CONFIG_SUPERTX
        chosen_rate_nocoef = INT_MAX;
#endif
        break;
      }

      chosen_rdc.rate += tmp_rdc.rate;
      chosen_rdc.dist += tmp_rdc.dist;
#if CONFIG_SUPERTX
      chosen_rate_nocoef += rt_nocoef;
#endif

      if (i != 3)
        encode_sb(cpi, tile, tp,  mi_row + y_idx, mi_col + x_idx, 0,
                  split_subsize, pc_tree->split[i]);

      pl = partition_plane_context(xd, mi_row + y_idx, mi_col + x_idx,
                                   split_subsize);
      chosen_rdc.rate += cpi->partition_cost[pl][PARTITION_NONE];
#if CONFIG_SUPERTX
      chosen_rate_nocoef += cpi->partition_cost[pl][PARTITION_SPLIT];
#endif
    }
    pl = partition_plane_context(xd, mi_row, mi_col, bsize);
    if (chosen_rdc.rate < INT_MAX) {
      chosen_rdc.rate += cpi->partition_cost[pl][PARTITION_SPLIT];
      chosen_rdc.rdcost = RDCOST(x->rdmult, x->rddiv,
                                 chosen_rdc.rate, chosen_rdc.dist);
#if CONFIG_SUPERTX
      chosen_rate_nocoef += cpi->partition_cost[pl][PARTITION_NONE];
#endif
    }
  }

  // If last_part is better set the partitioning to that.
  if (last_part_rdc.rdcost < chosen_rdc.rdcost) {
    mi_8x8[0].src_mi->mbmi.sb_type = bsize;
    if (bsize >= BLOCK_8X8)
      pc_tree->partitioning = partition;
    chosen_rdc = last_part_rdc;
#if CONFIG_SUPERTX
    chosen_rate_nocoef = last_part_rate_nocoef;
#endif
  }
  // If none was better set the partitioning to that.
  if (none_rdc.rdcost < chosen_rdc.rdcost) {
    if (bsize >= BLOCK_8X8)
      pc_tree->partitioning = PARTITION_NONE;
    chosen_rdc = none_rdc;
#if CONFIG_SUPERTX
    chosen_rate_nocoef = none_rate_nocoef;
#endif
  }

  restore_context(cpi, mi_row, mi_col, a, l, sa, sl, bsize);

  // We must have chosen a partitioning and encoding or we'll fail later on.
  // No other opportunities for success.
  if (bsize == BLOCK_LARGEST)
    assert(chosen_rdc.rate < INT_MAX && chosen_rdc.dist < INT64_MAX);

  if (do_recon) {
    int output_enabled = (bsize == BLOCK_LARGEST);

    // Check the projected output rate for this SB against it's target
    // and and if necessary apply a Q delta using segmentation to get
    // closer to the target.
    if ((cpi->oxcf.aq_mode == COMPLEXITY_AQ) && cm->seg.update_map) {
      vp9_select_in_frame_q_segment(cpi, mi_row, mi_col,
                                    output_enabled, chosen_rdc.rate);
    }

    if (cpi->oxcf.aq_mode == CYCLIC_REFRESH_AQ)
      vp9_cyclic_refresh_set_rate_and_dist_sb(cpi->cyclic_refresh,
                                              chosen_rdc.rate, chosen_rdc.dist);
    encode_sb(cpi, tile, tp, mi_row, mi_col, output_enabled, bsize,
              pc_tree);
  }

  *rate = chosen_rdc.rate;
  *dist = chosen_rdc.dist;
#if CONFIG_SUPERTX
  *rate_nocoef = chosen_rate_nocoef;
#endif
}

static const BLOCK_SIZE min_partition_size[BLOCK_SIZES] = {
  BLOCK_4X4,   BLOCK_4X4,   BLOCK_4X4,
  BLOCK_4X4,   BLOCK_4X4,   BLOCK_4X4,
  BLOCK_8X8,   BLOCK_8X8,   BLOCK_8X8,
  BLOCK_16X16, BLOCK_16X16, BLOCK_16X16,
  BLOCK_16X16,
#if CONFIG_EXT_CODING_UNIT_SIZE
  BLOCK_16X16, BLOCK_16X16, BLOCK_16X16
#endif
};

static const BLOCK_SIZE max_partition_size[BLOCK_SIZES] = {
  BLOCK_8X8,   BLOCK_16X16, BLOCK_16X16,
  BLOCK_16X16, BLOCK_32X32, BLOCK_32X32,
  BLOCK_32X32, BLOCK_64X64, BLOCK_64X64,
  BLOCK_64X64, BLOCK_64X64, BLOCK_64X64,
  BLOCK_64X64,
#if CONFIG_EXT_CODING_UNIT_SIZE
  BLOCK_64X64, BLOCK_64X64, BLOCK_128X128
#endif
};

// Look at all the mode_info entries for blocks that are part of this
// partition and find the min and max values for sb_type.
// At the moment this is designed to work on a 64x64 SB but could be
// adjusted to use a size parameter.
//
// The min and max are assumed to have been initialized prior to calling this
// function so repeat calls can accumulate a min and max of more than one sb64.
static void get_sb_partition_size_range(MACROBLOCKD *xd, MODE_INFO *mi_8x8,
                                        BLOCK_SIZE *min_block_size,
                                        BLOCK_SIZE *max_block_size,
                                        int bs_hist[BLOCK_SIZES]) {
  int sb_width_in_blocks = MI_BLOCK_SIZE;
  int sb_height_in_blocks  = MI_BLOCK_SIZE;
  int i, j;
  int index = 0;

  // Check the sb_type for each block that belongs to this region.
  for (i = 0; i < sb_height_in_blocks; ++i) {
    for (j = 0; j < sb_width_in_blocks; ++j) {
      MODE_INFO *mi = mi_8x8[index+j].src_mi;
      BLOCK_SIZE sb_type = mi ? mi->mbmi.sb_type : 0;
      bs_hist[sb_type]++;
      *min_block_size = MIN(*min_block_size, sb_type);
      *max_block_size = MAX(*max_block_size, sb_type);
    }
    index += xd->mi_stride;
  }
}

// Next square block size less or equal than current block size.
static const BLOCK_SIZE next_square_size[BLOCK_SIZES] = {
  BLOCK_4X4, BLOCK_4X4, BLOCK_4X4,
  BLOCK_8X8, BLOCK_8X8, BLOCK_8X8,
  BLOCK_16X16, BLOCK_16X16, BLOCK_16X16,
  BLOCK_32X32, BLOCK_32X32, BLOCK_32X32,
  BLOCK_64X64,
#if CONFIG_EXT_CODING_UNIT_SIZE
  BLOCK_64X64, BLOCK_64X64, BLOCK_128X128
#endif
};

// Look at neighboring blocks and set a min and max partition size based on
// what they chose.
static void rd_auto_partition_range(VP9_COMP *cpi, const TileInfo *const tile,
                                    int mi_row, int mi_col,
                                    BLOCK_SIZE *min_block_size,
                                    BLOCK_SIZE *max_block_size) {
  VP9_COMMON *const cm = &cpi->common;
  MACROBLOCKD *const xd = &cpi->mb.e_mbd;
  MODE_INFO *mi = xd->mi[0].src_mi;
  const int left_in_image = xd->left_available && mi[-1].src_mi;
  const int above_in_image = xd->up_available && mi[-xd->mi_stride].src_mi;
  const int row8x8_remaining = tile->mi_row_end - mi_row;
  const int col8x8_remaining = tile->mi_col_end - mi_col;
  int bh, bw;
  BLOCK_SIZE min_size = BLOCK_4X4;
  BLOCK_SIZE max_size = BLOCK_LARGEST;
  int i = 0;
  int bs_hist[BLOCK_SIZES] = {0};

  // Trap case where we do not have a prediction.
  if (left_in_image || above_in_image || cm->frame_type != KEY_FRAME) {
    // Default "min to max" and "max to min"
    min_size = BLOCK_LARGEST;
    max_size = BLOCK_4X4;

    // NOTE: each call to get_sb_partition_size_range() uses the previous
    // passed in values for min and max as a starting point.
    // Find the min and max partition used in previous frame at this location
    if (cm->frame_type != KEY_FRAME) {
      MODE_INFO *prev_mi =
          cm->prev_mip + cm->mi_stride + 1 + mi_row * xd->mi_stride + mi_col;

      get_sb_partition_size_range(xd, prev_mi, &min_size, &max_size, bs_hist);
    }
    // Find the min and max partition sizes used in the left SB64
    if (left_in_image) {
      MODE_INFO *left_sb64_mi = mi[-MI_BLOCK_SIZE].src_mi;
      get_sb_partition_size_range(xd, left_sb64_mi, &min_size, &max_size,
                                  bs_hist);
    }
    // Find the min and max partition sizes used in the above SB64.
    if (above_in_image) {
      MODE_INFO *above_sb64_mi = mi[-xd->mi_stride * MI_BLOCK_SIZE].src_mi;
      get_sb_partition_size_range(xd, above_sb64_mi, &min_size, &max_size,
                                  bs_hist);
    }

    // adjust observed min and max
    if (cpi->sf.auto_min_max_partition_size == RELAXED_NEIGHBORING_MIN_MAX) {
      min_size = min_partition_size[min_size];
      max_size = max_partition_size[max_size];
    } else if (cpi->sf.auto_min_max_partition_size ==
               CONSTRAIN_NEIGHBORING_MIN_MAX) {
      // adjust the search range based on the histogram of the observed
      // partition sizes from left, above the previous co-located blocks
      int sum = 0;
      int first_moment = 0;
      int second_moment = 0;
      int var_unnormalized = 0;

      for (i = 0; i < BLOCK_SIZES; i++) {
        sum += bs_hist[i];
        first_moment += bs_hist[i] * i;
        second_moment += bs_hist[i] * i * i;
      }

      // if variance is small enough,
      // adjust the range around its mean size, which gives a tighter range
      var_unnormalized = second_moment - first_moment * first_moment / sum;
      if (var_unnormalized <= 4 * sum) {
        int mean = first_moment / sum;
        min_size = min_partition_size[mean];
        max_size = max_partition_size[mean];
      } else {
        min_size = min_partition_size[min_size];
        max_size = max_partition_size[max_size];
      }
    }
  }

  // Check border cases where max and min from neighbors may not be legal.
  max_size = find_partition_size(max_size,
                                 row8x8_remaining, col8x8_remaining,
                                 &bh, &bw);
  min_size = MIN(min_size, max_size);

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
static void set_partition_range(VP9_COMMON *cm, MACROBLOCKD *xd,
                                int mi_row, int mi_col, BLOCK_SIZE bsize,
                                BLOCK_SIZE *min_bs, BLOCK_SIZE *max_bs) {
  int mi_width  = num_8x8_blocks_wide_lookup[bsize];
  int mi_height = num_8x8_blocks_high_lookup[bsize];
  int idx, idy;

  MODE_INFO *mi;
  const int idx_str = cm->mi_stride * mi_row + mi_col;
  MODE_INFO *prev_mi = (cm->prev_mip + cm->mi_stride + 1 + idx_str)->src_mi;


  BLOCK_SIZE bs, min_size, max_size;

  min_size = BLOCK_LARGEST;
  max_size = BLOCK_4X4;

  if (prev_mi) {
    for (idy = 0; idy < mi_height; ++idy) {
      for (idx = 0; idx < mi_width; ++idx) {
        mi = prev_mi[idy * cm->mi_stride + idx].src_mi;
        bs = mi ? mi->mbmi.sb_type : bsize;
        min_size = MIN(min_size, bs);
        max_size = MAX(max_size, bs);
      }
    }
  }

  if (xd->left_available) {
    for (idy = 0; idy < mi_height; ++idy) {
      mi = xd->mi[idy * cm->mi_stride - 1].src_mi;
      bs = mi ? mi->mbmi.sb_type : bsize;
      min_size = MIN(min_size, bs);
      max_size = MAX(max_size, bs);
    }
  }

  if (xd->up_available) {
    for (idx = 0; idx < mi_width; ++idx) {
      mi = xd->mi[idx - cm->mi_stride].src_mi;
      bs = mi ? mi->mbmi.sb_type : bsize;
      min_size = MIN(min_size, bs);
      max_size = MAX(max_size, bs);
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
  vpx_memcpy(ctx->pred_mv, x->pred_mv, sizeof(x->pred_mv));
}

static INLINE void load_pred_mv(MACROBLOCK *x, PICK_MODE_CONTEXT *ctx) {
  vpx_memcpy(x->pred_mv, ctx->pred_mv, sizeof(x->pred_mv));
}

#if CONFIG_FP_MB_STATS
const int num_16x16_blocks_wide_lookup[BLOCK_SIZES] =
  {1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 4, 4};
const int num_16x16_blocks_high_lookup[BLOCK_SIZES] =
  {1, 1, 1, 1, 1, 1, 1, 2, 1, 2, 4, 2, 4};
const int qindex_skip_threshold_lookup[BLOCK_SIZES] =
  {0, 10, 10, 30, 40, 40, 60, 80, 80, 90, 100, 100, 120};
const int qindex_split_threshold_lookup[BLOCK_SIZES] =
  {0, 3, 3, 7, 15, 15, 30, 40, 40, 60, 80, 80, 120};
const int complexity_16x16_blocks_threshold[BLOCK_SIZES] =
  {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 4, 4, 6};

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

#if CONFIG_EXT_PARTITION
static void rd_test_partition3(VP9_COMP *cpi, const TileInfo *const tile,
                               TOKENEXTRA **tp, PC_TREE *pc_tree,
                               RD_COST *best_rdc, PICK_MODE_CONTEXT ctxs[3],
                               PICK_MODE_CONTEXT *ctx,
                               int mi_row, int mi_col, BLOCK_SIZE bsize,
                               PARTITION_TYPE partition,
#if CONFIG_SUPERTX
                               int64_t best_rd, int *best_rate_nocoef,
                               VP9_COMMON *const cm,
                               ENTROPY_CONTEXT l[16 * MAX_MB_PLANE],
                               ENTROPY_CONTEXT a[16 * MAX_MB_PLANE],
                               PARTITION_CONTEXT sl[8],
                               PARTITION_CONTEXT sa[8],
#endif
#if CONFIG_PALETTE
                               int previous_size,
                               int previous_count[PALETTE_BUF_SIZE],
#if CONFIG_VP9_HIGHBITDEPTH
                               uint16_t previous_colors[PALETTE_BUF_SIZE],
#else
                               uint8_t previous_colors[PALETTE_BUF_SIZE],
#endif
#endif
                               int mi_row0, int mi_col0, BLOCK_SIZE subsize0,
                               int mi_row1, int mi_col1, BLOCK_SIZE subsize1,
                               int mi_row2, int mi_col2, BLOCK_SIZE subsize2) {
  MACROBLOCK *const x = &cpi->mb;
  MACROBLOCKD *const xd = &x->e_mbd;
  RD_COST this_rdc, sum_rdc;
#if CONFIG_SUPERTX
  int this_rate_nocoef, sum_rate_nocoef;
  int abort_flag;
  PARTITION_TYPE best_partition;
  int tmp_rate;
  int64_t tmp_dist, tmp_rd;
#endif
#if CONFIG_PALETTE
  PICK_MODE_CONTEXT *c, *p;
#endif
  if (cpi->sf.adaptive_motion_search)
    load_pred_mv(x, ctx);

#if CONFIG_PALETTE
    c = &ctxs[0];
    c->palette_buf_size = previous_size;
    vpx_memcpy(c->palette_colors_buf, previous_colors,
               previous_size * sizeof(previous_colors[0]));
    vpx_memcpy(c->palette_count_buf, previous_count,
               previous_size * sizeof(previous_count[0]));
#endif

  rd_pick_sb_modes(cpi, tile, mi_row0, mi_col0, &sum_rdc,
#if CONFIG_SUPERTX
                   &sum_rate_nocoef,
#endif
                   subsize0, &ctxs[0], best_rdc->rdcost);
#if CONFIG_SUPERTX
  abort_flag = sum_rdc.rdcost >= best_rd;
#endif

#if CONFIG_SUPERTX
  if (sum_rdc.rdcost < INT64_MAX) {
#else
  if (sum_rdc.rdcost < best_rdc->rdcost) {
#endif
    PICK_MODE_CONTEXT *ctx = &ctxs[0];
    update_state(cpi, ctx, mi_row0, mi_col0, subsize0, 0);
    encode_superblock(cpi, tp, 0, mi_row0, mi_col0, subsize0, ctx);

    if (cpi->sf.adaptive_motion_search)
      load_pred_mv(x, ctx);

#if CONFIG_PALETTE
      copy_palette_info(&ctxs[1], &ctxs[0]);
#endif

#if CONFIG_SUPERTX
    rd_pick_sb_modes(cpi, tile, mi_row1, mi_col1, &this_rdc,
                     &this_rate_nocoef, subsize1, &ctxs[1],
                     INT64_MAX - sum_rdc.rdcost);
#else
    rd_pick_sb_modes(cpi, tile, mi_row1, mi_col1, &this_rdc, subsize1, &ctxs[1],
                     best_rdc->rdcost - sum_rdc.rdcost);
#endif

    if (this_rdc.rate == INT_MAX) {
      sum_rdc.rdcost = INT64_MAX;
#if CONFIG_SUPERTX
      sum_rate_nocoef = INT_MAX;
#endif
    } else {
      sum_rdc.rate += this_rdc.rate;
      sum_rdc.dist += this_rdc.dist;
      sum_rdc.rdcost += this_rdc.rdcost;
#if CONFIG_SUPERTX
      sum_rate_nocoef += this_rate_nocoef;
#endif
    }

#if CONFIG_SUPERTX
    if (sum_rdc.rdcost < INT64_MAX) {
#else
    if (sum_rdc.rdcost < best_rdc->rdcost) {
#endif
      PICK_MODE_CONTEXT *ctx = &ctxs[1];
      update_state(cpi, ctx, mi_row1, mi_col1, subsize1, 0);
      encode_superblock(cpi, tp, 0, mi_row1, mi_col1, subsize1, ctx);

      if (cpi->sf.adaptive_motion_search)
        load_pred_mv(x, ctx);

#if CONFIG_PALETTE
      copy_palette_info(&ctxs[2], &ctxs[1]);
#endif

#if CONFIG_SUPERTX
      rd_pick_sb_modes(cpi, tile, mi_row2, mi_col2, &this_rdc,
                       &this_rate_nocoef, subsize2, &ctxs[2],
                       INT64_MAX - sum_rdc.rdcost);
#else
      rd_pick_sb_modes(cpi, tile, mi_row2, mi_col2, &this_rdc, subsize2,
                       &ctxs[2], best_rdc->rdcost - sum_rdc.rdcost);
#endif

      if (this_rdc.rate == INT_MAX) {
        sum_rdc.rdcost = INT64_MAX;
#if CONFIG_SUPERTX
        sum_rate_nocoef = INT_MAX;
#endif
      } else {
        sum_rdc.rate += this_rdc.rate;
        sum_rdc.dist += this_rdc.dist;
        sum_rdc.rdcost += this_rdc.rdcost;
#if CONFIG_SUPERTX
        sum_rate_nocoef += this_rate_nocoef;
#endif
      }

#if CONFIG_SUPERTX
      if (cm->frame_type != KEY_FRAME && !abort_flag &&
          sum_rdc.rdcost < INT64_MAX && bsize <= MAX_SUPERTX_BLOCK_SIZE &&
          !xd->lossless) {
        TX_SIZE supertx_size = bsize_to_tx_size(bsize);
        best_partition = pc_tree->partitioning;
        pc_tree->partitioning = partition;
        sum_rdc.rate += vp9_cost_bit(
            cm->fc.supertx_prob
            [partition_supertx_context_lookup[partition]][supertx_size],
            0);
        sum_rdc.rdcost = RDCOST(x->rdmult, x->rddiv, sum_rdc.rate,
                                sum_rdc.dist);

        if (!check_intra_sb(cpi, tile, mi_row, mi_col, bsize, pc_tree)) {
#if CONFIG_EXT_TX
          EXT_TX_TYPE best_tx = NORM;
#endif
#if CONFIG_NEW_QUANT && QUANT_PROFILES > 1
          int dq_index = 0;
#endif  // CONFIG_NEW_QUANT && QUANT_PROFILES > 1

          tmp_rate = sum_rate_nocoef;
          tmp_dist = 0;
          restore_context(cpi, mi_row, mi_col, a, l, sa, sl, bsize);
          rd_supertx_sb(cpi, tile, mi_row, mi_col, bsize, &tmp_rate, &tmp_dist,
#if CONFIG_EXT_TX
                        &best_tx,
#endif
#if CONFIG_NEW_QUANT && QUANT_PROFILES > 1
                        &dq_index,
#endif  // CONFIG_NEW_QUANT && QUANT_PROFILES > 1
                        pc_tree);

          tmp_rate += vp9_cost_bit(
              cm->fc.supertx_prob
              [partition_supertx_context_lookup[partition]][supertx_size],
              1);
          tmp_rd = RDCOST(x->rdmult, x->rddiv, tmp_rate, tmp_dist);
          if (tmp_rd < sum_rdc.rdcost) {
            sum_rdc.rdcost = tmp_rd;
            sum_rdc.rate = tmp_rate;
            sum_rdc.dist = tmp_dist;
            update_supertx_param_sb(cpi, mi_row, mi_col, bsize,
#if CONFIG_EXT_TX
                                    best_tx,
#endif
#if CONFIG_NEW_QUANT && QUANT_PROFILES > 1
                                    dq_index,
#endif  // CONFIG_NEW_QUANT && QUANT_PROFILES > 1
                                    supertx_size, pc_tree);
          }
        }
        pc_tree->partitioning = best_partition;
      }
#endif  // CONFIG_SUPERTX

      if (sum_rdc.rdcost < best_rdc->rdcost) {
        int pl = partition_plane_context(xd, mi_row, mi_col, bsize);
        sum_rdc.rate += cpi->partition_cost[pl][partition];
        sum_rdc.rdcost = RDCOST(x->rdmult, x->rddiv, sum_rdc.rate,
                                sum_rdc.dist);
#if CONFIG_SUPERTX
        sum_rate_nocoef += cpi->partition_cost[pl][partition];
#endif
        if (sum_rdc.rdcost < best_rdc->rdcost) {
#if CONFIG_SUPERTX
          *best_rate_nocoef = sum_rate_nocoef;
          assert(*best_rate_nocoef >= 0);
#endif
          *best_rdc = sum_rdc;
          pc_tree->partitioning = partition;
#if CONFIG_PALETTE
          c = &pc_tree->current;
          p = &ctxs[2];
          copy_palette_info(c, p);
#endif
        }
      }
    }
  }
}
#endif

// TODO(jingning,jimbankoski,rbultje): properly skip partition types that are
// unlikely to be selected depending on previous rate-distortion optimization
// results, for encoding speed-up.
static void rd_pick_partition(VP9_COMP *cpi, const TileInfo *const tile,
                              TOKENEXTRA **tp, int mi_row, int mi_col,
                              BLOCK_SIZE bsize, RD_COST *rd_cost,
#if CONFIG_SUPERTX
                              int *rate_nocoef,
#endif
                              int64_t best_rd, PC_TREE *pc_tree) {
  VP9_COMMON *const cm = &cpi->common;
  MACROBLOCK *const x = &cpi->mb;
  MACROBLOCKD *const xd = &x->e_mbd;
  const int mi_step = num_8x8_blocks_wide_lookup[bsize] / 2;
  ENTROPY_CONTEXT l[(CODING_UNIT_SIZE >> 2) * MAX_MB_PLANE];
  ENTROPY_CONTEXT a[(CODING_UNIT_SIZE >> 2) * MAX_MB_PLANE];
  PARTITION_CONTEXT sl[CODING_UNIT_SIZE >> 3], sa[CODING_UNIT_SIZE >> 3];
  TOKENEXTRA *tp_orig = *tp;
  PICK_MODE_CONTEXT *ctx = &pc_tree->none;
  int i, pl;
  BLOCK_SIZE subsize;
  RD_COST this_rdc, sum_rdc, best_rdc;
#if CONFIG_SUPERTX
  int this_rate_nocoef, sum_rate_nocoef = 0, best_rate_nocoef = INT_MAX;
  int tmp_rate;
  int abort_flag;
  int64_t tmp_dist, tmp_rd;
  PARTITION_TYPE best_partition;
#endif
  int do_split = bsize >= BLOCK_8X8;
  int do_rect = 1;
#if CONFIG_EXT_PARTITION
  BLOCK_SIZE bsize2 = get_subsize(bsize, PARTITION_SPLIT);
#endif

  // Override skipping rectangular partition operations for edge blocks
  const int force_horz_split = (mi_row + mi_step >= cm->mi_rows);
  const int force_vert_split = (mi_col + mi_step >= cm->mi_cols);
  const int xss = x->e_mbd.plane[1].subsampling_x;
  const int yss = x->e_mbd.plane[1].subsampling_y;

  BLOCK_SIZE min_size = cpi->sf.min_partition_size;
  BLOCK_SIZE max_size = cpi->sf.max_partition_size;

#if CONFIG_FP_MB_STATS
  unsigned int src_diff_var = UINT_MAX;
  int none_complexity = 0;
#endif

  int partition_none_allowed = !force_horz_split && !force_vert_split;
  int partition_horz_allowed = !force_vert_split && yss <= xss &&
                               bsize >= BLOCK_8X8;
  int partition_vert_allowed = !force_horz_split && xss <= yss &&
                               bsize >= BLOCK_8X8;
#if CONFIG_PALETTE
  PICK_MODE_CONTEXT *c, *p;
  int previous_size, previous_count[PALETTE_BUF_SIZE];
#if CONFIG_VP9_HIGHBITDEPTH
  uint16_t previous_colors[PALETTE_BUF_SIZE];
#else
  uint8_t previous_colors[PALETTE_BUF_SIZE];
#endif  // CONFIG_VP9_HIGHBITDEPTH
#endif  // CONFIG_PALETTE
  (void) *tp_orig;

  assert(num_8x8_blocks_wide_lookup[bsize] ==
             num_8x8_blocks_high_lookup[bsize]);

  vp9_rd_cost_init(&this_rdc);
  vp9_rd_cost_init(&sum_rdc);
  vp9_rd_cost_reset(&best_rdc);
  best_rdc.rdcost = best_rd;

  set_offsets(cpi, tile, mi_row, mi_col, bsize);

#if CONFIG_PALETTE
  if (bsize == BLOCK_LARGEST) {
    c = &pc_tree->current;
    c->palette_buf_size = cm->current_palette_size;
    vpx_memcpy(c->palette_colors_buf, cm->current_palette_colors,
               c->palette_buf_size * sizeof(cm->current_palette_colors[0]));
    vpx_memcpy(c->palette_count_buf, cm->current_palette_count,
               c->palette_buf_size * sizeof(cm->current_palette_count[0]));
  }

  c = &pc_tree->current;
  previous_size = c->palette_buf_size;
  vpx_memcpy(previous_colors, c->palette_colors_buf,
             previous_size * sizeof(previous_colors[0]));
  vpx_memcpy(previous_count, c->palette_count_buf,
             previous_size * sizeof(previous_count[0]));

  c = &pc_tree->none;
  p = &pc_tree->current;
  copy_palette_info(c, p);
#endif

  if (bsize == BLOCK_16X16 && cpi->oxcf.aq_mode)
    x->mb_energy = vp9_block_energy(cpi, x, bsize);

  if (cpi->sf.cb_partition_search && bsize == BLOCK_16X16) {
    int cb_partition_search_ctrl = ((pc_tree->index == 0 || pc_tree->index == 3)
        + get_chessboard_index(cm->current_video_frame)) & 0x1;

    if (cb_partition_search_ctrl && bsize > min_size && bsize < max_size)
      set_partition_range(cm, xd, mi_row, mi_col, bsize, &min_size, &max_size);
  }

  // Determine partition types in search according to the speed features.
  // The threshold set here has to be of square block size.
  if (cpi->sf.auto_min_max_partition_size) {
    partition_none_allowed &= (bsize <= max_size && bsize >= min_size);
    partition_horz_allowed &= ((bsize <= max_size && bsize > min_size) ||
                                force_horz_split);
    partition_vert_allowed &= ((bsize <= max_size && bsize > min_size) ||
                                force_vert_split);
    do_split &= bsize > min_size;
  }
  if (cpi->sf.use_square_partition_only) {
    partition_horz_allowed &= force_horz_split;
    partition_vert_allowed &= force_vert_split;
  }

  save_context(cpi, mi_row, mi_col, a, l, sa, sl, bsize);

#if CONFIG_FP_MB_STATS
  if (cpi->use_fp_mb_stats) {
    set_offsets(cpi, tile, mi_row, mi_col, bsize);
    src_diff_var = get_sby_perpixel_diff_variance(cpi, &cpi->mb.plane[0].src,
                                                  mi_row, mi_col, bsize);
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
        MIN(mb_row + num_16x16_blocks_high_lookup[bsize], cm->mb_rows);
    int mb_col_end =
        MIN(mb_col + num_16x16_blocks_wide_lookup[bsize], cm->mb_cols);
    int r, c;

    // compute a complexity measure, basically measure inconsistency of motion
    // vectors obtained from the first pass in the current block
    for (r = mb_row; r < mb_row_end ; r++) {
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
    rd_pick_sb_modes(cpi, tile, mi_row, mi_col, &this_rdc,
#if CONFIG_SUPERTX
                     &this_rate_nocoef,
#endif
                     bsize, ctx, best_rdc.rdcost);
    if (this_rdc.rate != INT_MAX) {
      if (bsize >= BLOCK_8X8) {
        pl = partition_plane_context(xd, mi_row, mi_col, bsize);
        this_rdc.rate += cpi->partition_cost[pl][PARTITION_NONE];
        this_rdc.rdcost = RDCOST(x->rdmult, x->rddiv,
                                 this_rdc.rate, this_rdc.dist);
#if CONFIG_SUPERTX
        this_rate_nocoef += cpi->partition_cost[pl][PARTITION_NONE];
#endif
      }

      if (this_rdc.rdcost < best_rdc.rdcost) {
        int64_t dist_breakout_thr = cpi->sf.partition_search_breakout_dist_thr;
        int rate_breakout_thr = cpi->sf.partition_search_breakout_rate_thr;

        best_rdc = this_rdc;
#if CONFIG_SUPERTX
        best_rate_nocoef = this_rate_nocoef;
        assert(best_rate_nocoef >= 0);
#endif
        if (bsize >= BLOCK_8X8)
          pc_tree->partitioning = PARTITION_NONE;

        // Adjust dist breakout threshold according to the partition size.
        dist_breakout_thr >>= 8 - (b_width_log2_lookup[bsize] +
            b_height_log2_lookup[bsize]);

        rate_breakout_thr *= num_pels_log2_lookup[bsize];

        // If all y, u, v transform blocks in this partition are skippable, and
        // the dist & rate are within the thresholds, the partition search is
        // terminated for current branch of the partition search tree.
        // The dist & rate thresholds are set to 0 at speed 0 to disable the
        // early termination at that speed.
        if (!x->e_mbd.lossless &&
            (ctx->skippable && best_rdc.dist < dist_breakout_thr &&
            best_rdc.rate < rate_breakout_thr)) {
          do_split = 0;
          do_rect = 0;
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
              MIN(mb_row + num_16x16_blocks_high_lookup[bsize], cm->mb_rows);
          int mb_col_end =
              MIN(mb_col + num_16x16_blocks_wide_lookup[bsize], cm->mb_cols);
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
              set_offsets(cpi, tile, mi_row, mi_col, bsize);
              src_diff_var = get_sby_perpixel_diff_variance(
                  cpi, &cpi->mb.plane[0].src, mi_row, mi_col, bsize);
            }
            if (src_diff_var < 8) {
              do_split = 0;
              do_rect = 0;
            }
          }
        }
#endif
#if CONFIG_PALETTE
        c = &pc_tree->current;
        p = &pc_tree->none;
        copy_palette_info(c, p);
#endif
      }
    }
    restore_context(cpi, mi_row, mi_col, a, l, sa, sl, bsize);
  }

  // store estimated motion vector
  if (cpi->sf.adaptive_motion_search)
    store_pred_mv(x, ctx);

  // PARTITION_SPLIT
  // TODO(jingning): use the motion vectors given by the above search as
  // the starting point of motion search in the following partition type check.
  if (do_split) {
#if CONFIG_PALETTE
    int last = -1;
#endif
    subsize = get_subsize(bsize, PARTITION_SPLIT);
    if (bsize == BLOCK_8X8) {
      i = 4;
      if (cpi->sf.adaptive_pred_interp_filter && partition_none_allowed)
        pc_tree->leaf_split[0]->pred_interp_filter =
            ctx->mic.mbmi.interp_filter;
#if CONFIG_SUPERTX
      rd_pick_sb_modes(cpi, tile, mi_row, mi_col, &sum_rdc, &sum_rate_nocoef,
                       subsize, pc_tree->leaf_split[0], INT64_MAX);
#else
      rd_pick_sb_modes(cpi, tile, mi_row, mi_col, &sum_rdc,
                       subsize, pc_tree->leaf_split[0], best_rdc.rdcost);
#endif
      if (sum_rdc.rate == INT_MAX) {
        sum_rdc.rdcost = INT64_MAX;
#if CONFIG_SUPERTX
        sum_rate_nocoef = INT_MAX;
#endif
      }
#if CONFIG_SUPERTX
      if (cm->frame_type != KEY_FRAME && sum_rdc.rdcost < INT64_MAX &&
          !xd->lossless) {
        TX_SIZE supertx_size = bsize_to_tx_size(bsize);  // b_width_log2(bsize);
        best_partition = pc_tree->partitioning;
        pc_tree->partitioning = PARTITION_SPLIT;

        sum_rdc.rate += vp9_cost_bit(
            cm->fc.supertx_prob
            [partition_supertx_context_lookup[PARTITION_SPLIT]][supertx_size],
            0);
        sum_rdc.rdcost =
            RDCOST(x->rdmult, x->rddiv, sum_rdc.rate, sum_rdc.dist);
#if CONFIG_NEW_INTER
        if (is_inter_mode(pc_tree->leaf_split[0]->mic.mbmi.mode) ||
            is_inter_compound_mode(pc_tree->leaf_split[0]->mic.mbmi.mode)) {
#else
        if (is_inter_mode(pc_tree->leaf_split[0]->mic.mbmi.mode)) {
#endif  // CONFIG_NEW_INTER
#if CONFIG_EXT_TX
          EXT_TX_TYPE best_tx = NORM;
#endif
#if CONFIG_NEW_QUANT && QUANT_PROFILES > 1
          int dq_index = 0;
#endif  // CONFIG_NEW_QUANT && QUANT_PROFILES > 1

          tmp_rate = sum_rate_nocoef;
          tmp_dist = 0;
          restore_context(cpi, mi_row, mi_col, a, l, sa, sl, bsize);
          rd_supertx_sb(cpi, tile, mi_row, mi_col, bsize, &tmp_rate, &tmp_dist,
#if CONFIG_EXT_TX
                        &best_tx,
#endif
#if CONFIG_NEW_QUANT && QUANT_PROFILES > 1
                        &dq_index,
#endif  // CONFIG_NEW_QUANT && QUANT_PROFILES > 1
                        pc_tree);

          tmp_rate += vp9_cost_bit(
              cm->fc.supertx_prob
              [partition_supertx_context_lookup[PARTITION_SPLIT]][supertx_size],
              1);
          tmp_rd = RDCOST(x->rdmult, x->rddiv, tmp_rate, tmp_dist);
          if (tmp_rd < sum_rdc.rdcost) {
            sum_rdc.rdcost = tmp_rd;
            sum_rdc.rate = tmp_rate;
            sum_rdc.dist = tmp_dist;
            update_supertx_param_sb(cpi, mi_row, mi_col, bsize,
#if CONFIG_EXT_TX
                                    best_tx,
#endif
#if CONFIG_NEW_QUANT && QUANT_PROFILES > 1
                                    dq_index,
#endif  // CONFIG_NEW_QUANT && QUANT_PROFILES > 1
                                    supertx_size, pc_tree);
          }
        }
        pc_tree->partitioning = best_partition;
      }
#endif  // CONFIG_SUPERTX
    } else {
#if CONFIG_SUPERTX
      for (i = 0; i < 4 && sum_rdc.rdcost < INT64_MAX; ++i) {
#else
      for (i = 0; i < 4 && sum_rdc.rdcost < best_rdc.rdcost; ++i) {
#endif
      const int x_idx = (i & 1) * mi_step;
      const int y_idx = (i >> 1) * mi_step;

        if (mi_row + y_idx >= cm->mi_rows || mi_col + x_idx >= cm->mi_cols)
          continue;

        if (cpi->sf.adaptive_motion_search)
          load_pred_mv(x, ctx);

        pc_tree->split[i]->index = i;

#if CONFIG_PALETTE
        c = &pc_tree->split[i]->current;
        if (last < 0) {
          c->palette_buf_size = previous_size;
          vpx_memcpy(c->palette_colors_buf, previous_colors,
                     previous_size * sizeof(previous_colors[0]));
          vpx_memcpy(c->palette_count_buf, previous_count,
                     previous_size * sizeof(previous_count[0]));
        } else {
          p = &pc_tree->split[last]->current;
          copy_palette_info(c, p);
        }
        last = i;
#endif

#if CONFIG_SUPERTX
        rd_pick_partition(cpi, tile, tp, mi_row + y_idx, mi_col + x_idx,
                          subsize, &this_rdc, &this_rate_nocoef,
                          INT64_MAX - sum_rdc.rdcost, pc_tree->split[i]);
#else
        rd_pick_partition(cpi, tile, tp, mi_row + y_idx, mi_col + x_idx,
                          subsize, &this_rdc,
                          best_rdc.rdcost - sum_rdc.rdcost, pc_tree->split[i]);
#endif

        if (this_rdc.rate == INT_MAX) {
          sum_rdc.rdcost = INT64_MAX;
#if CONFIG_SUPERTX
          sum_rate_nocoef = INT_MAX;
#endif
          break;
        } else {
          sum_rdc.rate += this_rdc.rate;
          sum_rdc.dist += this_rdc.dist;
          sum_rdc.rdcost += this_rdc.rdcost;
#if CONFIG_SUPERTX
          sum_rate_nocoef += this_rate_nocoef;
#endif
        }
      }
#if CONFIG_SUPERTX
      if (cm->frame_type != KEY_FRAME && sum_rdc.rdcost < INT64_MAX &&
          i == 4 && bsize <= MAX_SUPERTX_BLOCK_SIZE  && !xd->lossless) {
        TX_SIZE supertx_size = bsize_to_tx_size(bsize);
        best_partition = pc_tree->partitioning;
        pc_tree->partitioning = PARTITION_SPLIT;

        sum_rdc.rate += vp9_cost_bit(
            cm->fc.supertx_prob
            [partition_supertx_context_lookup[PARTITION_SPLIT]][supertx_size],
            0);
        sum_rdc.rdcost =
            RDCOST(x->rdmult, x->rddiv, sum_rdc.rate, sum_rdc.dist);

        if (!check_intra_sb(cpi, tile, mi_row, mi_col, bsize, pc_tree)) {
#if CONFIG_EXT_TX
          EXT_TX_TYPE best_tx = NORM;
#endif
#if CONFIG_NEW_QUANT && QUANT_PROFILES > 1
          int dq_index = 0;
#endif  // CONFIG_NEW_QUANT && QUANT_PROFILES > 1

          tmp_rate = sum_rate_nocoef;
          tmp_dist = 0;
          restore_context(cpi, mi_row, mi_col, a, l, sa, sl, bsize);
          rd_supertx_sb(cpi, tile, mi_row, mi_col, bsize, &tmp_rate, &tmp_dist,
#if CONFIG_EXT_TX
                        &best_tx,
#endif
#if CONFIG_NEW_QUANT && QUANT_PROFILES > 1
                        &dq_index,
#endif  // CONFIG_NEW_QUANT && QUANT_PROFILES > 1
                        pc_tree);

          tmp_rate += vp9_cost_bit(
              cm->fc.supertx_prob
              [partition_supertx_context_lookup[PARTITION_SPLIT]][supertx_size],
              1);
          tmp_rd = RDCOST(x->rdmult, x->rddiv, tmp_rate, tmp_dist);
          if (tmp_rd < sum_rdc.rdcost) {
            sum_rdc.rdcost = tmp_rd;
            sum_rdc.rate = tmp_rate;
            sum_rdc.dist = tmp_dist;
            update_supertx_param_sb(cpi, mi_row, mi_col, bsize,
#if CONFIG_EXT_TX
                                    best_tx,
#endif
#if CONFIG_NEW_QUANT && QUANT_PROFILES > 1
                                    dq_index,
#endif  // CONFIG_NEW_QUANT && QUANT_PROFILES > 1
                                    supertx_size, pc_tree);
          }
        }
        pc_tree->partitioning = best_partition;
      }
#endif  // CONFIG_SUPERTX
    }

    if (sum_rdc.rdcost < best_rdc.rdcost && i == 4) {
      pl = partition_plane_context(xd, mi_row, mi_col, bsize);
      sum_rdc.rate += cpi->partition_cost[pl][PARTITION_SPLIT];
      sum_rdc.rdcost = RDCOST(x->rdmult, x->rddiv,
                              sum_rdc.rate, sum_rdc.dist);
#if CONFIG_SUPERTX
      sum_rate_nocoef += cpi->partition_cost[pl][PARTITION_SPLIT];
#endif

      if (sum_rdc.rdcost < best_rdc.rdcost) {
        best_rdc = sum_rdc;
#if CONFIG_SUPERTX
        best_rate_nocoef = sum_rate_nocoef;
        assert(best_rate_nocoef >= 0);
#endif
        pc_tree->partitioning = PARTITION_SPLIT;
#if CONFIG_PALETTE
        if (bsize > BLOCK_8X8 && last >= 0) {
          c = &pc_tree->current;
          p = &(pc_tree->split[last]->current);
          copy_palette_info(c, p);
        } else {
          c = &pc_tree->current;
          c->palette_buf_size = previous_size;
          vpx_memcpy(c->palette_colors_buf, previous_colors,
                     previous_size * sizeof(previous_colors[0]));
          vpx_memcpy(c->palette_count_buf, previous_count,
                     previous_size * sizeof(previous_count[0]));
        }
#endif
      }
    } else {
      // skip rectangular partition test when larger block size
      // gives better rd cost
      if (cpi->sf.less_rectangular_check)
        do_rect &= !partition_none_allowed;
    }
    restore_context(cpi, mi_row, mi_col, a, l, sa, sl, bsize);
  }

  // PARTITION_HORZ
  if (partition_horz_allowed && do_rect) {
#if CONFIG_PALETTE
    int last;
#endif
    subsize = get_subsize(bsize, PARTITION_HORZ);
    if (cpi->sf.adaptive_motion_search)
      load_pred_mv(x, ctx);
    if (cpi->sf.adaptive_pred_interp_filter && bsize == BLOCK_8X8 &&
        partition_none_allowed)
      pc_tree->horizontal[0].pred_interp_filter =
          ctx->mic.mbmi.interp_filter;

#if CONFIG_PALETTE
    c = &pc_tree->horizontal[0];
    c->palette_buf_size = previous_size;
    vpx_memcpy(c->palette_colors_buf, previous_colors,
               previous_size * sizeof(previous_colors[0]));
    vpx_memcpy(c->palette_count_buf, previous_count,
               previous_size * sizeof(previous_count[0]));
    last = 0;
#endif

    rd_pick_sb_modes(cpi, tile, mi_row, mi_col, &sum_rdc,
#if CONFIG_SUPERTX
                     &sum_rate_nocoef,
#endif
                     subsize, &pc_tree->horizontal[0], best_rdc.rdcost);
#if CONFIG_SUPERTX
    abort_flag = (sum_rdc.rdcost >= best_rd && bsize > BLOCK_8X8) ||
                 (sum_rdc.rate == INT_MAX && bsize == BLOCK_8X8);
    if (sum_rdc.rdcost < INT64_MAX &&
#else
    if (sum_rdc.rdcost < best_rdc.rdcost &&
#endif
        mi_row + mi_step < cm->mi_rows &&
        bsize > BLOCK_8X8) {
      PICK_MODE_CONTEXT *ctx = &pc_tree->horizontal[0];
      update_state(cpi, ctx, mi_row, mi_col, subsize, 0);
      encode_superblock(cpi, tp, 0, mi_row, mi_col, subsize, ctx);

      if (cpi->sf.adaptive_motion_search)
        load_pred_mv(x, ctx);
      if (cpi->sf.adaptive_pred_interp_filter && bsize == BLOCK_8X8 &&
          partition_none_allowed)
        pc_tree->horizontal[1].pred_interp_filter =
            ctx->mic.mbmi.interp_filter;

#if CONFIG_PALETTE
      copy_palette_info(&pc_tree->horizontal[1], &pc_tree->horizontal[0]);
      last = 1;
#endif

#if CONFIG_SUPERTX
      rd_pick_sb_modes(cpi, tile, mi_row + mi_step, mi_col, &this_rdc,
                       &this_rate_nocoef,
                       subsize, &pc_tree->horizontal[1],
                       INT64_MAX);
#else
      rd_pick_sb_modes(cpi, tile, mi_row + mi_step, mi_col, &this_rdc,
                       subsize, &pc_tree->horizontal[1],
                       best_rdc.rdcost - sum_rdc.rdcost);
#endif
      if (this_rdc.rate == INT_MAX) {
        sum_rdc.rdcost = INT64_MAX;
#if CONFIG_SUPERTX
        sum_rate_nocoef = INT_MAX;
#endif
      } else {
        sum_rdc.rate += this_rdc.rate;
        sum_rdc.dist += this_rdc.dist;
        sum_rdc.rdcost += this_rdc.rdcost;
#if CONFIG_SUPERTX
        sum_rate_nocoef += this_rate_nocoef;
#endif
      }
    }
#if CONFIG_SUPERTX
    if (cm->frame_type != KEY_FRAME && !abort_flag &&
        sum_rdc.rdcost < INT64_MAX && bsize <= MAX_SUPERTX_BLOCK_SIZE &&
        !xd->lossless) {
      TX_SIZE supertx_size = bsize_to_tx_size(bsize);
      best_partition = pc_tree->partitioning;
      pc_tree->partitioning = PARTITION_HORZ;

      sum_rdc.rate += vp9_cost_bit(
          cm->fc.supertx_prob[partition_supertx_context_lookup[PARTITION_HORZ]]
                             [supertx_size], 0);
      sum_rdc.rdcost = RDCOST(x->rdmult, x->rddiv, sum_rdc.rate, sum_rdc.dist);

      if (!check_intra_sb(cpi, tile, mi_row, mi_col, bsize, pc_tree)) {
#if CONFIG_EXT_TX
        EXT_TX_TYPE best_tx = NORM;
#endif
#if CONFIG_NEW_QUANT && QUANT_PROFILES > 1
        int dq_index = 0;
#endif  // CONFIG_NEW_QUANT && QUANT_PROFILES > 1

        tmp_rate = sum_rate_nocoef;
        tmp_dist = 0;
        restore_context(cpi, mi_row, mi_col, a, l, sa, sl, bsize);
        rd_supertx_sb(cpi, tile, mi_row, mi_col, bsize, &tmp_rate, &tmp_dist,
#if CONFIG_EXT_TX
                      &best_tx,
#endif
#if CONFIG_NEW_QUANT && QUANT_PROFILES > 1
                      &dq_index,
#endif  // CONFIG_NEW_QUANT && QUANT_PROFILES > 1 && !Q_CTX_BASED_PROFILES
                      pc_tree);

        tmp_rate += vp9_cost_bit(
            cm->fc.supertx_prob
            [partition_supertx_context_lookup[PARTITION_HORZ]][supertx_size],
            1);
        tmp_rd = RDCOST(x->rdmult, x->rddiv, tmp_rate, tmp_dist);
        if (tmp_rd < sum_rdc.rdcost) {
          sum_rdc.rdcost = tmp_rd;
          sum_rdc.rate = tmp_rate;
          sum_rdc.dist = tmp_dist;
          update_supertx_param_sb(cpi, mi_row, mi_col, bsize,
#if CONFIG_EXT_TX
                                  best_tx,
#endif
#if CONFIG_NEW_QUANT && QUANT_PROFILES > 1
                                  dq_index,
#endif  // CONFIG_NEW_QUANT && QUANT_PROFILES > 1
                                  supertx_size, pc_tree);
        }
      }
      pc_tree->partitioning = best_partition;
    }
#endif  // CONFIG_SUPERTX

    if (sum_rdc.rdcost < best_rdc.rdcost) {
      pl = partition_plane_context(xd, mi_row, mi_col, bsize);
      sum_rdc.rate += cpi->partition_cost[pl][PARTITION_HORZ];
      sum_rdc.rdcost = RDCOST(x->rdmult, x->rddiv, sum_rdc.rate, sum_rdc.dist);
#if CONFIG_SUPERTX
      sum_rate_nocoef += cpi->partition_cost[pl][PARTITION_HORZ];
#endif
      if (sum_rdc.rdcost < best_rdc.rdcost) {
        best_rdc = sum_rdc;
#if CONFIG_SUPERTX
        best_rate_nocoef = sum_rate_nocoef;
        assert(best_rate_nocoef >= 0);
#endif
        pc_tree->partitioning = PARTITION_HORZ;
#if CONFIG_PALETTE
        if (bsize > BLOCK_8X8) {
          c = &pc_tree->current;
          p = &pc_tree->horizontal[last];
          copy_palette_info(c, p);
        } else {
          c = &pc_tree->current;
          c->palette_buf_size = previous_size;
          vpx_memcpy(c->palette_colors_buf, previous_colors,
                     previous_size * sizeof(previous_colors[0]));
          vpx_memcpy(c->palette_count_buf, previous_count,
                     previous_size * sizeof(previous_count[0]));
        }
#endif
      }
    }
    restore_context(cpi, mi_row, mi_col, a, l, sa, sl, bsize);
  }
  // PARTITION_VERT
  if (partition_vert_allowed && do_rect) {
#if CONFIG_PALETTE
    int last;
#endif
    subsize = get_subsize(bsize, PARTITION_VERT);

    if (cpi->sf.adaptive_motion_search)
      load_pred_mv(x, ctx);
    if (cpi->sf.adaptive_pred_interp_filter && bsize == BLOCK_8X8 &&
        partition_none_allowed)
      pc_tree->vertical[0].pred_interp_filter =
          ctx->mic.mbmi.interp_filter;

#if CONFIG_PALETTE
    c = &pc_tree->vertical[0];
    c->palette_buf_size = previous_size;
    vpx_memcpy(c->palette_colors_buf, previous_colors,
               previous_size * sizeof(previous_colors[0]));
    vpx_memcpy(c->palette_count_buf, previous_count,
               previous_size * sizeof(previous_count[0]));
    last = 0;
#endif

    rd_pick_sb_modes(cpi, tile, mi_row, mi_col, &sum_rdc,
#if CONFIG_SUPERTX
                     &sum_rate_nocoef,
#endif
                     subsize, &pc_tree->vertical[0], best_rdc.rdcost);
#if CONFIG_SUPERTX
    abort_flag = (sum_rdc.rdcost >= best_rd && bsize > BLOCK_8X8) ||
                 (sum_rdc.rate == INT_MAX && bsize == BLOCK_8X8);
    if (sum_rdc.rdcost < INT64_MAX &&
#else
    if (sum_rdc.rdcost < best_rdc.rdcost &&
#endif
        mi_col + mi_step < cm->mi_cols &&
        bsize > BLOCK_8X8) {
      update_state(cpi, &pc_tree->vertical[0], mi_row, mi_col, subsize, 0);
      encode_superblock(cpi, tp, 0, mi_row, mi_col, subsize,
                        &pc_tree->vertical[0]);

      if (cpi->sf.adaptive_motion_search)
        load_pred_mv(x, ctx);
      if (cpi->sf.adaptive_pred_interp_filter && bsize == BLOCK_8X8 &&
          partition_none_allowed)
        pc_tree->vertical[1].pred_interp_filter =
            ctx->mic.mbmi.interp_filter;

#if CONFIG_PALETTE
      copy_palette_info(&pc_tree->vertical[1], &pc_tree->vertical[0]);
      last = 1;
#endif

#if CONFIG_SUPERTX
      rd_pick_sb_modes(cpi, tile, mi_row, mi_col + mi_step, &this_rdc,
                       &this_rate_nocoef, subsize, &pc_tree->vertical[1],
                       INT64_MAX - sum_rdc.rdcost);
#else
      rd_pick_sb_modes(cpi, tile, mi_row, mi_col + mi_step, &this_rdc, subsize,
                       &pc_tree->vertical[1], best_rdc.rdcost - sum_rdc.rdcost);
#endif
      if (this_rdc.rate == INT_MAX) {
        sum_rdc.rdcost = INT64_MAX;
#if CONFIG_SUPERTX
        sum_rate_nocoef = INT_MAX;
#endif
      } else {
        sum_rdc.rate += this_rdc.rate;
        sum_rdc.dist += this_rdc.dist;
        sum_rdc.rdcost += this_rdc.rdcost;
#if CONFIG_SUPERTX
        sum_rate_nocoef += this_rate_nocoef;
#endif
      }
    }
#if CONFIG_SUPERTX
    if (cm->frame_type != KEY_FRAME && !abort_flag &&
        sum_rdc.rdcost < INT64_MAX && bsize <= MAX_SUPERTX_BLOCK_SIZE &&
        !xd->lossless) {
      TX_SIZE supertx_size = bsize_to_tx_size(bsize);
      best_partition = pc_tree->partitioning;
      pc_tree->partitioning = PARTITION_VERT;
      sum_rdc.rate += vp9_cost_bit(
          cm->fc.supertx_prob[partition_supertx_context_lookup[PARTITION_VERT]]
                             [supertx_size], 0);
      sum_rdc.rdcost = RDCOST(x->rdmult, x->rddiv, sum_rdc.rate, sum_rdc.dist);

      if (!check_intra_sb(cpi, tile, mi_row, mi_col, bsize, pc_tree)) {
#if CONFIG_EXT_TX
        EXT_TX_TYPE best_tx = NORM;
#endif
#if CONFIG_NEW_QUANT && QUANT_PROFILES > 1
        int dq_index = 0;
#endif  // CONFIG_NEW_QUANT && QUANT_PROFILES > 1

        tmp_rate = sum_rate_nocoef;
        tmp_dist = 0;
        restore_context(cpi, mi_row, mi_col, a, l, sa, sl, bsize);
        rd_supertx_sb(cpi, tile, mi_row, mi_col, bsize, &tmp_rate, &tmp_dist,
#if CONFIG_EXT_TX
                      &best_tx,
#endif
#if CONFIG_NEW_QUANT && QUANT_PROFILES > 1
                      &dq_index,
#endif  // CONFIG_NEW_QUANT && QUANT_PROFILES > 1
                      pc_tree);

        tmp_rate += vp9_cost_bit(
            cm->fc.supertx_prob
            [partition_supertx_context_lookup[PARTITION_VERT]][supertx_size],
            1);
        tmp_rd = RDCOST(x->rdmult, x->rddiv, tmp_rate, tmp_dist);
        if (tmp_rd < sum_rdc.rdcost) {
          sum_rdc.rdcost = tmp_rd;
          sum_rdc.rate = tmp_rate;
          sum_rdc.dist = tmp_dist;
          update_supertx_param_sb(cpi, mi_row, mi_col, bsize,
#if CONFIG_EXT_TX
                                  best_tx,
#endif
#if CONFIG_NEW_QUANT && QUANT_PROFILES > 1
                                  dq_index,
#endif  // CONFIG_NEW_QUANT && QUANT_PROFILES > 1
                                  supertx_size, pc_tree);
        }
      }
      pc_tree->partitioning = best_partition;
    }
#endif  // CONFIG_SUPERTX

    if (sum_rdc.rdcost < best_rdc.rdcost) {
      pl = partition_plane_context(xd, mi_row, mi_col, bsize);
      sum_rdc.rate += cpi->partition_cost[pl][PARTITION_VERT];
      sum_rdc.rdcost = RDCOST(x->rdmult, x->rddiv,
                              sum_rdc.rate, sum_rdc.dist);
#if CONFIG_SUPERTX
      sum_rate_nocoef += cpi->partition_cost[pl][PARTITION_VERT];
#endif
      if (sum_rdc.rdcost < best_rdc.rdcost) {
        best_rdc = sum_rdc;
#if CONFIG_SUPERTX
        best_rate_nocoef = sum_rate_nocoef;
        assert(best_rate_nocoef >= 0);
#endif
        pc_tree->partitioning = PARTITION_VERT;
#if CONFIG_PALETTE
        if (bsize > BLOCK_8X8) {
          c = &pc_tree->current;
          p = &pc_tree->vertical[last];
          copy_palette_info(c, p);
        } else {
          c = &pc_tree->current;
          c->palette_buf_size = previous_size;
          vpx_memcpy(c->palette_colors_buf, previous_colors,
                     previous_size * sizeof(previous_colors[0]));
          vpx_memcpy(c->palette_count_buf, previous_count,
                     previous_size * sizeof(previous_count[0]));
        }
#endif
      }
    }
    restore_context(cpi, mi_row, mi_col, a, l, sa, sl, bsize);
  }

#if CONFIG_EXT_PARTITION
  // PARTITION_HORZ_A
  if (partition_horz_allowed && do_rect && bsize > BLOCK_8X8 &&
      partition_none_allowed) {
    subsize = get_subsize(bsize, PARTITION_HORZ_A);
    rd_test_partition3(cpi, tile, tp, pc_tree, &best_rdc,
                       pc_tree->horizontala,
                       ctx, mi_row, mi_col, bsize, PARTITION_HORZ_A,
#if CONFIG_SUPERTX
                       best_rd, &best_rate_nocoef, cm, l, a, sl, sa,
#endif
#if CONFIG_PALETTE
                       previous_size, previous_count, previous_colors,
#endif
                       mi_row, mi_col, bsize2,
                       mi_row, mi_col + mi_step, bsize2,
                       mi_row + mi_step, mi_col, subsize);
    restore_context(cpi, mi_row, mi_col, a, l, sa, sl, bsize);
  }
  // PARTITION_HORZ_B
  if (partition_horz_allowed && do_rect && bsize > BLOCK_8X8 &&
      partition_none_allowed) {
    subsize = get_subsize(bsize, PARTITION_HORZ_B);
    rd_test_partition3(cpi, tile, tp, pc_tree, &best_rdc,
                       pc_tree->horizontalb,
                       ctx, mi_row, mi_col, bsize, PARTITION_HORZ_B,
#if CONFIG_SUPERTX
                       best_rd, &best_rate_nocoef, cm, l, a, sl, sa,
#endif
#if CONFIG_PALETTE
                       previous_size, previous_count, previous_colors,
#endif
                       mi_row, mi_col, subsize,
                       mi_row + mi_step, mi_col, bsize2,
                       mi_row + mi_step, mi_col + mi_step, bsize2);
    restore_context(cpi, mi_row, mi_col, a, l, sa, sl, bsize);
  }
  // PARTITION_VERT_A
  if (partition_vert_allowed && do_rect && bsize > BLOCK_8X8 &&
      partition_none_allowed) {
    subsize = get_subsize(bsize, PARTITION_VERT_A);
    rd_test_partition3(cpi, tile, tp, pc_tree, &best_rdc,
                       pc_tree->verticala,
                       ctx, mi_row, mi_col, bsize, PARTITION_VERT_A,
#if CONFIG_SUPERTX
                       best_rd, &best_rate_nocoef, cm, l, a, sl, sa,
#endif
#if CONFIG_PALETTE
                       previous_size, previous_count, previous_colors,
#endif
                       mi_row, mi_col, bsize2,
                       mi_row + mi_step, mi_col, bsize2,
                       mi_row, mi_col + mi_step, subsize);
    restore_context(cpi, mi_row, mi_col, a, l, sa, sl, bsize);
  }
  // PARTITION_VERT_B
  if (partition_vert_allowed && do_rect && bsize > BLOCK_8X8 &&
      partition_none_allowed) {
    subsize = get_subsize(bsize, PARTITION_VERT_B);
    rd_test_partition3(cpi, tile, tp, pc_tree, &best_rdc,
                       pc_tree->verticalb,
                       ctx, mi_row, mi_col, bsize, PARTITION_VERT_B,
#if CONFIG_SUPERTX
                       best_rd, &best_rate_nocoef, cm, l, a, sl, sa,
#endif
#if CONFIG_PALETTE
                       previous_size, previous_count, previous_colors,
#endif
                       mi_row, mi_col, subsize,
                       mi_row, mi_col + mi_step, bsize2,
                       mi_row + mi_step, mi_col + mi_step, bsize2);
    restore_context(cpi, mi_row, mi_col, a, l, sa, sl, bsize);
  }
#endif

  // TODO(jbb): This code added so that we avoid static analysis
  // warning related to the fact that best_rd isn't used after this
  // point.  This code should be refactored so that the duplicate
  // checks occur in some sub function and thus are used...
  (void) best_rd;
  *rd_cost = best_rdc;
#if CONFIG_SUPERTX
  *rate_nocoef = best_rate_nocoef;
#endif

  if (best_rdc.rate < INT_MAX && best_rdc.dist < INT64_MAX &&
      pc_tree->index != 3) {
    int output_enabled = (bsize == BLOCK_LARGEST);

    // Check the projected output rate for this SB against it's target
    // and and if necessary apply a Q delta using segmentation to get
    // closer to the target.
    if ((cpi->oxcf.aq_mode == COMPLEXITY_AQ) && cm->seg.update_map)
      vp9_select_in_frame_q_segment(cpi, mi_row, mi_col, output_enabled,
                                    best_rdc.rate);
    if (cpi->oxcf.aq_mode == CYCLIC_REFRESH_AQ)
      vp9_cyclic_refresh_set_rate_and_dist_sb(cpi->cyclic_refresh,
                                              best_rdc.rate, best_rdc.dist);

    encode_sb(cpi, tile, tp, mi_row, mi_col, output_enabled, bsize, pc_tree);
  }

  if (bsize == BLOCK_LARGEST) {
    assert(tp_orig < *tp);
    assert(best_rdc.rate < INT_MAX);
    assert(best_rdc.dist < INT64_MAX);
  } else {
    assert(tp_orig == *tp);
  }
}


static void encode_rd_sb_row(VP9_COMP *cpi, const TileInfo *const tile,
                             int mi_row, TOKENEXTRA **tp) {
  VP9_COMMON *const cm = &cpi->common;
  MACROBLOCKD *const xd = &cpi->mb.e_mbd;
  SPEED_FEATURES *const sf = &cpi->sf;
  int mi_col;
#if CONFIG_EXT_CODING_UNIT_SIZE
  const int leaf_nodes = 64 * 4;
#else
  const int leaf_nodes = 64;
#endif

  // Initialize the left context for the new SB row
  vpx_memset(&xd->left_context, 0, sizeof(xd->left_context));
  vpx_memset(xd->left_seg_context, 0, sizeof(xd->left_seg_context));

  // Code each SB in the row
  for (mi_col = tile->mi_col_start; mi_col < tile->mi_col_end;
       mi_col += MI_BLOCK_SIZE) {
    int dummy_rate;
    int64_t dummy_dist;
    RD_COST dummy_rdc;
#if CONFIG_SUPERTX
    int dummy_rate_nocoef;
#endif
    int i;

    const int idx_str = cm->mi_stride * mi_row + mi_col;
    MODE_INFO *mi = cm->mi + idx_str;

    if (sf->adaptive_pred_interp_filter) {
      for (i = 0; i < leaf_nodes; ++i)
        cpi->leaf_tree[i].pred_interp_filter = SWITCHABLE;

      for (i = 0; i < leaf_nodes; ++i) {
        cpi->pc_tree[i].vertical[0].pred_interp_filter = SWITCHABLE;
        cpi->pc_tree[i].vertical[1].pred_interp_filter = SWITCHABLE;
        cpi->pc_tree[i].horizontal[0].pred_interp_filter = SWITCHABLE;
        cpi->pc_tree[i].horizontal[1].pred_interp_filter = SWITCHABLE;
      }
    }

    vp9_zero(cpi->mb.pred_mv);
    cpi->pc_root->index = 0;

    cpi->mb.source_variance = UINT_MAX;
    if (sf->partition_search_type == FIXED_PARTITION) {
      set_offsets(cpi, tile, mi_row, mi_col, BLOCK_LARGEST);
      set_fixed_partitioning(cpi, tile, mi, mi_row, mi_col,
                             sf->always_this_block_size);
      rd_use_partition(cpi, tile, mi, tp, mi_row, mi_col, BLOCK_LARGEST,
                       &dummy_rate, &dummy_dist,
#if CONFIG_SUPERTX
                       &dummy_rate_nocoef,
#endif
                       1, cpi->pc_root);
    } else if (cpi->partition_search_skippable_frame) {
      BLOCK_SIZE bsize;
      set_offsets(cpi, tile, mi_row, mi_col, BLOCK_LARGEST);
      bsize = get_rd_var_based_fixed_partition(cpi, mi_row, mi_col);
      set_fixed_partitioning(cpi, tile, mi, mi_row, mi_col, bsize);
      rd_use_partition(cpi, tile, mi, tp, mi_row, mi_col, BLOCK_LARGEST,
                       &dummy_rate, &dummy_dist,
#if CONFIG_SUPERTX
                       &dummy_rate_nocoef,
#endif
                       1, cpi->pc_root);
    } else if (sf->partition_search_type == VAR_BASED_PARTITION &&
               cm->frame_type != KEY_FRAME ) {
      choose_partitioning(cpi, tile, mi_row, mi_col);
      rd_use_partition(cpi, tile, mi, tp, mi_row, mi_col, BLOCK_LARGEST,
                       &dummy_rate, &dummy_dist,
#if CONFIG_SUPERTX
                       &dummy_rate_nocoef,
#endif
                       1, cpi->pc_root);
    } else {
      // If required set upper and lower partition size limits
      if (sf->auto_min_max_partition_size) {
        set_offsets(cpi, tile, mi_row, mi_col, BLOCK_LARGEST);
        rd_auto_partition_range(cpi, tile, mi_row, mi_col,
                                &sf->min_partition_size,
                                &sf->max_partition_size);
      }
      rd_pick_partition(cpi, tile, tp, mi_row, mi_col, BLOCK_LARGEST,
                        &dummy_rdc,
#if CONFIG_SUPERTX
                        &dummy_rate_nocoef,
#endif
                        INT64_MAX, cpi->pc_root);
    }
  }
}

static void init_encode_frame_mb_context(VP9_COMP *cpi) {
  MACROBLOCK *const x = &cpi->mb;
  VP9_COMMON *const cm = &cpi->common;
  MACROBLOCKD *const xd = &x->e_mbd;
  const int aligned_mi_cols = mi_cols_aligned_to_sb(cm->mi_cols);

  // Copy data over into macro block data structures.
  vp9_setup_src_planes(x, cpi->Source, 0, 0);

  vp9_setup_block_planes(&x->e_mbd, cm->subsampling_x, cm->subsampling_y);

  // Note: this memset assumes above_context[0], [1] and [2]
  // are allocated as part of the same buffer.
  vpx_memset(xd->above_context[0], 0,
             sizeof(*xd->above_context[0]) *
             2 * aligned_mi_cols * MAX_MB_PLANE);
  vpx_memset(xd->above_seg_context, 0,
             sizeof(*xd->above_seg_context) * aligned_mi_cols);
}

static int check_dual_ref_flags(VP9_COMP *cpi) {
  const int ref_flags = cpi->ref_frame_flags;

  if (vp9_segfeature_active(&cpi->common.seg, 1, SEG_LVL_REF_FRAME)) {
    return 0;
  } else {
    return (!!(ref_flags & VP9_GOLD_FLAG) +
            !!(ref_flags & VP9_LAST_FLAG) +
#if CONFIG_MULTI_REF
            !!(ref_flags & VP9_LAST2_FLAG) +
            !!(ref_flags & VP9_LAST3_FLAG) +
            !!(ref_flags & VP9_LAST4_FLAG) +
#endif  // CONFIG_MULTI_REF
            !!(ref_flags & VP9_ALT_FLAG)) >= 2;
  }
}

static void reset_skip_tx_size(VP9_COMMON *cm, TX_SIZE max_tx_size) {
  int mi_row, mi_col;
  const int mis = cm->mi_stride;
  MODE_INFO *mi_ptr = cm->mi;

  for (mi_row = 0; mi_row < cm->mi_rows; ++mi_row, mi_ptr += mis) {
    for (mi_col = 0; mi_col < cm->mi_cols; ++mi_col) {
      if (mi_ptr[mi_col].src_mi->mbmi.tx_size > max_tx_size)
        mi_ptr[mi_col].src_mi->mbmi.tx_size = max_tx_size;
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
    // TODO(zoeliu): To investigate whether a frame_type of LAST2_FRAME needs to
    // be analyzed here to decide on the reference mode.
    return LAST_FRAME;
}

static TX_MODE select_tx_mode(const VP9_COMP *cpi) {
#if !CONFIG_TX_SKIP
  if (cpi->mb.e_mbd.lossless)
    return ONLY_4X4;
#endif
  if (cpi->sf.tx_size_search_method == USE_LARGESTALL)
#if CONFIG_TX64X64
    return ALLOW_64X64;
#else
    return ALLOW_32X32;
#endif
  else if (cpi->sf.tx_size_search_method == USE_FULL_RD||
           cpi->sf.tx_size_search_method == USE_TX_8X8)
    return TX_MODE_SELECT;
  else
    return cpi->common.tx_mode;
}

static int get_skip_encode_frame(const VP9_COMMON *cm) {
  unsigned int intra_count = 0, inter_count = 0;
  int j;

  for (j = 0; j < INTRA_INTER_CONTEXTS; ++j) {
    intra_count += cm->counts.intra_inter[j][0];
    inter_count += cm->counts.intra_inter[j][1];
  }

  return (intra_count << 2) < inter_count &&
         cm->frame_type != KEY_FRAME &&
         cm->show_frame;
}

static void encode_tiles(VP9_COMP *cpi) {
  const VP9_COMMON *const cm = &cpi->common;
  const int tile_cols = cm->tile_cols;
  const int tile_rows = cm->tile_rows;

  int tile_col, tile_row;
#if CONFIG_ROW_TILE
  TileInfo (*tile)[1024] = cpi->tile_info;
  TOKENEXTRA *(*tok)[1024] = cpi->tile_tok;
#else
  TileInfo tile[4][1 << 6];
  TOKENEXTRA *tok[4][1 << 6];
#endif
  TOKENEXTRA *pre_tok = cpi->tok;
  int tile_tok = 0;

  for (tile_row = 0; tile_row < tile_rows; ++tile_row) {
    for (tile_col = 0; tile_col < tile_cols; ++tile_col) {
      vp9_tile_init(&tile[tile_row][tile_col], cm, tile_row, tile_col);
      tok[tile_row][tile_col] = pre_tok + tile_tok;
      pre_tok = tok[tile_row][tile_col];
      tile_tok = allocated_tokens(tile[tile_row][tile_col]);
    }
  }

  for (tile_row = 0; tile_row < tile_rows; ++tile_row) {
    for (tile_col = 0; tile_col < tile_cols; ++tile_col) {
      const TileInfo * const ptile = &tile[tile_row][tile_col];
      TOKENEXTRA * const old_tok = tok[tile_row][tile_col];
      int mi_row;

#if CONFIG_ROW_TILE
      int col_width =
          mi_cols_aligned_to_sb(ptile->mi_col_end - ptile->mi_col_start);
      vpx_memset(cm->above_context, 0, sizeof(*cm->above_context) *
                 MAX_MB_PLANE * 2 * mi_cols_aligned_to_sb(cm->mi_cols));
      vpx_memset(&cm->above_seg_context[ptile->mi_col_start], 0,
                 sizeof(*cm->above_seg_context) * col_width);
#endif

      for (mi_row = ptile->mi_row_start; mi_row < ptile->mi_row_end;
           mi_row += MI_BLOCK_SIZE) {
          encode_rd_sb_row(cpi, ptile, mi_row, &tok[tile_row][tile_col]);
      }
      cpi->tok_count[tile_row][tile_col] =
          (unsigned int)(tok[tile_row][tile_col] - old_tok);
      assert(tok[tile_row][tile_col] - old_tok <= allocated_tokens(*ptile));
    }
  }
}

#if CONFIG_FP_MB_STATS
static int input_fpmb_stats(FIRSTPASS_MB_STATS *firstpass_mb_stats,
                            VP9_COMMON *cm, uint8_t **this_frame_mb_stats) {
  uint8_t *mb_stats_in = firstpass_mb_stats->mb_stats_start +
      cm->current_video_frame * cm->MBs * sizeof(uint8_t);

  if (mb_stats_in > firstpass_mb_stats->mb_stats_end)
    return EOF;

  *this_frame_mb_stats = mb_stats_in;

  return 1;
}
#endif

#if CONFIG_GLOBAL_MOTION

#define MIN_TRANSLATION_THRESH 8
#define GLOBAL_MOTION_MODEL  ROTZOOM
#define GLOBAL_MOTION_ADVANTAGE_THRESH_RZ 0.60
#define GLOBAL_MOTION_ADVANTAGE_THRESH_TR 0.75
// #define USE_BLOCK_BASED_GLOBAL_MOTION_COMPUTATION
#define USE_FEATURE_BASED_GLOBAL_MOTION_COMPUTATION

static void convert_translation_to_params(
    double *H, Global_Motion_Params *model) {
  model->mv.as_mv.col = (int) floor(H[0] * 8 + 0.5);
  model->mv.as_mv.row = (int) floor(H[1] * 8 + 0.5);
  if (abs(model->mv.as_mv.col) < MIN_TRANSLATION_THRESH &&
      abs(model->mv.as_mv.row) < MIN_TRANSLATION_THRESH) {
    model->mv.as_int = 0;
  } else {
    model->mv.as_mv.col =
        clamp(model->mv.as_mv.col,
              -(1 << ABS_TRANSLATION_BITS), (1 << ABS_TRANSLATION_BITS));
    model->mv.as_mv.row =
        clamp(model->mv.as_mv.row,
              -(1 << ABS_TRANSLATION_BITS), (1 << ABS_TRANSLATION_BITS));
  }
}

static void convert_rotzoom_to_params(double *H, Global_Motion_Params *model) {
  double z = sqrt(H[0] * H[0] + H[1] * H[1]) - 1.0;
  double r = atan2(-H[1], H[0]) * 180.0 / M_PI;
  assert(abs(H[0] - (1 + z) * cos(r * M_PI / 180.0)) < 1e-10);
  assert(abs(H[1] + (1 + z) * sin(r * M_PI / 180.0)) < 1e-10);
  model->zoom = (int) floor(z * (1 << ZOOM_PRECISION_BITS) + 0.5);
  model->rotation = (int) floor(r * (1 << ROTATION_PRECISION_BITS) + 0.5);
  model->zoom = clamp(
      model->zoom, -(1 << ABS_ZOOM_BITS), (1 << ABS_ZOOM_BITS));
  model->rotation = clamp(
      model->rotation, -(1 << ABS_ROTATION_BITS), (1 << ABS_ROTATION_BITS));

  model->mv.as_mv.col = (int) floor(H[2] * 8 + 0.5);
  model->mv.as_mv.row = (int) floor(H[3] * 8 + 0.5);
  model->mv.as_mv.col =
      clamp(model->mv.as_mv.col,
            -(1 << ABS_TRANSLATION_BITS), (1 << ABS_TRANSLATION_BITS));
  model->mv.as_mv.row =
      clamp(model->mv.as_mv.row,
            -(1 << ABS_TRANSLATION_BITS), (1 << ABS_TRANSLATION_BITS));

  if (model->zoom == 0 && model->rotation == 0) {
    if (abs(model->mv.as_mv.col) < MIN_TRANSLATION_THRESH &&
        abs(model->mv.as_mv.row) < MIN_TRANSLATION_THRESH) {
      model->mv.as_int = 0;
    }
  }
}

static void convert_model_to_params(double *H, TransformationType type,
                                    Global_Motion_Params *model) {
  switch (type) {
    case ROTZOOM:
      convert_rotzoom_to_params(H, model);
      break;
    case TRANSLATION:
      convert_translation_to_params(H, model);
      break;
    default:
      break;
  }
  model->gmtype = get_gmtype(model);
}
#endif  // CONFIG_GLOBAL_MOTION

static void encode_frame_internal(VP9_COMP *cpi) {
  SPEED_FEATURES *const sf = &cpi->sf;
  RD_OPT *const rd_opt = &cpi->rd;
  MACROBLOCK *const x = &cpi->mb;
  VP9_COMMON *const cm = &cpi->common;
  MACROBLOCKD *const xd = &x->e_mbd;

  xd->mi = cm->mi;
  xd->mi[0].src_mi = &xd->mi[0];
  vp9_zero(cm->counts);
  vp9_zero(cpi->coef_counts);
  vp9_zero(rd_opt->comp_pred_diff);
  vp9_zero(rd_opt->filter_diff);
  vp9_zero(rd_opt->tx_select_diff);
  vp9_zero(rd_opt->tx_select_threshes);

  xd->lossless = cm->base_qindex == 0 &&
                 cm->y_dc_delta_q == 0 &&
                 cm->uv_dc_delta_q == 0 &&
                 cm->uv_ac_delta_q == 0;

  cm->tx_mode = select_tx_mode(cpi);

#if CONFIG_GLOBAL_MOTION
  vp9_clear_system_state();
  vp9_zero(cpi->global_motion_used);
  vpx_memset(cm->num_global_motion, 0, sizeof(cm->num_global_motion));
  cm->num_global_motion[LAST_FRAME] = 1;
  cm->num_global_motion[GOLDEN_FRAME] = 1;
  cm->num_global_motion[ALTREF_FRAME] = 1;
  if (cpi->common.frame_type == INTER_FRAME && cpi->Source) {
    YV12_BUFFER_CONFIG *ref_buf;
    int num, frame;
    double global_motion[9 * MAX_GLOBAL_MOTION_MODELS];

    for (frame = LAST_FRAME; frame <= ALTREF_FRAME; ++frame) {
      ref_buf = get_ref_frame_buffer(cpi, frame);
      if (ref_buf) {
        if ((num =
#ifdef USE_BLOCK_BASED_GLOBAL_MOTION_COMPUTATION
             vp9_compute_global_motion_multiple_block_based(
                 cpi, GLOBAL_MOTION_MODEL, cpi->Source, ref_buf,
                 BLOCK_16X16, MAX_GLOBAL_MOTION_MODELS, 0.5, global_motion))) {
#elif defined(USE_FEATURE_BASED_GLOBAL_MOTION_COMPUTATION)
             vp9_compute_global_motion_multiple_feature_based(
                 cpi, GLOBAL_MOTION_MODEL, cpi->Source, ref_buf,
                 MAX_GLOBAL_MOTION_MODELS, 0.5, global_motion))) {
#else
             vp9_compute_global_motion_multiple_optical_flow(
                 cpi, GLOBAL_MOTION_MODEL, cpi->Source, ref_buf,
                 MAX_GLOBAL_MOTION_MODELS, 0.5, global_motion))) {
#endif
          int i;
          for (i = 0; i < num; i++) {
            convert_model_to_params(
                global_motion + i * get_numparams(GLOBAL_MOTION_MODEL),
                GLOBAL_MOTION_MODEL,
                &cm->global_motion[frame][i]);
            refine_quant_param(&cm->global_motion[frame][i],
                               GLOBAL_MOTION_MODEL, ref_buf->y_buffer,
                               ref_buf->y_crop_width,
                               ref_buf->y_crop_height,
                               ref_buf->y_stride,
                               cpi->Source->y_buffer,
                               cpi->Source->y_crop_width,
                               cpi->Source->y_crop_height,
                               cpi->Source->y_stride, 3);

            if (get_gmtype(&cm->global_motion[frame][i]) != GLOBAL_ZERO) {
              double erroradvantage_trans;
              double erroradvantage =
                  vp9_warp_erroradv(&cm->global_motion[frame][i],
                                    ref_buf->y_buffer,
                                    ref_buf->y_crop_width,
                                    ref_buf->y_crop_height,
                                    ref_buf->y_stride,
                                    cpi->Source->y_buffer,
                                    0, 0,
                                    cpi->Source->y_crop_width,
                                    cpi->Source->y_crop_height,
                                    cpi->Source->y_stride,
                                    0, 0, 16, 16);
              if (get_gmtype(&cm->global_motion[frame][i]) == GLOBAL_ROTZOOM) {
                Global_Motion_Params gm = cm->global_motion[frame][i];
                gm.rotation = 0;
                gm.zoom = 0;
                gm.gmtype = GLOBAL_TRANSLATION;
                erroradvantage_trans =
                  vp9_warp_erroradv(&gm,
                                    ref_buf->y_buffer,
                                    ref_buf->y_crop_width,
                                    ref_buf->y_crop_height,
                                    ref_buf->y_stride,
                                    cpi->Source->y_buffer,
                                    0, 0,
                                    cpi->Source->y_crop_width,
                                    cpi->Source->y_crop_height,
                                    cpi->Source->y_stride,
                                    0, 0, 16, 16);
              } else {
                erroradvantage_trans = erroradvantage;
                erroradvantage = 10;
              }
              if (erroradvantage > GLOBAL_MOTION_ADVANTAGE_THRESH_RZ) {
                if (erroradvantage_trans > GLOBAL_MOTION_ADVANTAGE_THRESH_TR) {
                  // Not enough advantage in using a global model. Make 0.
                  vpx_memset(&cm->global_motion[frame][i], 0,
                             sizeof(cm->global_motion[frame][i]));
                } else {
                  if (cm->global_motion[frame][i].gmtype == GLOBAL_ROTZOOM) {
                    cm->global_motion[frame][i].rotation = 0;
                    cm->global_motion[frame][i].zoom = 0;
                    cm->global_motion[frame][i].gmtype =
                        get_gmtype(&cm->global_motion[frame][i]);
                  }
                }
              }
              if (get_gmtype(&cm->global_motion[frame][i]) != GLOBAL_ZERO)
                printf("Found it %d/%d - %d [%f %f]\n",
                       cm->current_video_frame, cm->show_frame, frame,
                       erroradvantage, erroradvantage_trans);
            }
          }
          cm->num_global_motion[frame] = num;
        }
      }
    }
  }

#endif  // CONFIG_GLOBAL_MOTION

#if CONFIG_VP9_HIGHBITDEPTH
  if (cm->use_highbitdepth)
    x->fwd_txm4x4 = xd->lossless ? vp9_highbd_fwht4x4 : vp9_highbd_fdct4x4;
  else
    x->fwd_txm4x4 = xd->lossless ? vp9_fwht4x4 : vp9_fdct4x4;
  x->highbd_itxm_add = xd->lossless ? vp9_highbd_iwht4x4_add :
      vp9_highbd_idct4x4_add;
#else
  x->fwd_txm4x4 = xd->lossless ? vp9_fwht4x4 : vp9_fdct4x4;
#endif  // CONFIG_VP9_HIGHBITDEPTH
  x->itxm_add = xd->lossless ? vp9_iwht4x4_add : vp9_idct4x4_add;
#if CONFIG_SR_MODE
  x->itxm = xd->lossless ? vp9_iwht4x4 : vp9_idct4x4;
#endif  // CONFIG_SR_MODE

  if (xd->lossless) {
    x->optimize = 0;
    cm->lf.filter_level = 0;
  }

  vp9_frame_init_quantizer(cpi);

  vp9_initialize_rd_consts(cpi);
  vp9_initialize_me_consts(cpi, cm->base_qindex);
  init_encode_frame_mb_context(cpi);
  set_prev_mi(cm);

  x->quant_fp = cpi->sf.use_quant_fp;
  vp9_zero(x->skip_txfm);

  {
    struct vpx_usec_timer emr_timer;
    vpx_usec_timer_start(&emr_timer);

#if CONFIG_FP_MB_STATS
    if (cpi->use_fp_mb_stats) {
      input_fpmb_stats(&cpi->twopass.firstpass_mb_stats, cm,
                       &cpi->twopass.this_frame_mb_stats);
    }
#endif

#if CONFIG_PALETTE
    if (frame_is_intra_only(cm)) {
      cm->current_palette_size = 0;
      vpx_memset(cm->current_palette_count, 0,
                 PALETTE_BUF_SIZE * sizeof(cm->current_palette_count[0]));
      cm->palette_counter = 0;
      cm->palette_blocks_signalled = 0;
    }
#endif  // CONFIG_PALETTE
#if CONFIG_INTRABC
    if (frame_is_intra_only(cm)) {
      cm->intrabc_counter = 0;
      cm->intrabc_blocks_signalled = 0;
    }
#endif  // CONFIG_INTRABC

    encode_tiles(cpi);

    vpx_usec_timer_mark(&emr_timer);
    cpi->time_encode_sb_row += vpx_usec_timer_elapsed(&emr_timer);
  }

  sf->skip_encode_frame = sf->skip_encode_sb ? get_skip_encode_frame(cm) : 0;

#if 0
  // Keep record of the total distortion this time around for future use
  cpi->last_frame_distortion = cpi->frame_distortion;
#endif
}

static INTERP_FILTER get_interp_filter(
    const int64_t threshes[SWITCHABLE_FILTER_CONTEXTS], int is_alt_ref) {
  if (!is_alt_ref &&
      threshes[EIGHTTAP_SMOOTH] > threshes[EIGHTTAP] &&
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

void vp9_encode_frame(VP9_COMP *cpi) {
  VP9_COMMON *const cm = &cpi->common;
  RD_OPT *const rd_opt = &cpi->rd;

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
      cm->allow_comp_inter_inter = 0;
    } else {
      cm->allow_comp_inter_inter = 1;
      cm->comp_fixed_ref = ALTREF_FRAME;
      cm->comp_var_ref[0] = LAST_FRAME;
#if CONFIG_MULTI_REF
      cm->comp_var_ref[1] = LAST2_FRAME;
      cm->comp_var_ref[2] = LAST3_FRAME;
      cm->comp_var_ref[3] = LAST4_FRAME;
      cm->comp_var_ref[4] = GOLDEN_FRAME;
#else  // CONFIG_MULTI_REF
      cm->comp_var_ref[1] = GOLDEN_FRAME;
#endif  // CONFIG_MULTI_REF
    }
  }

  if (cpi->sf.frame_parameter_update) {
    int i;

    // This code does a single RD pass over the whole frame assuming
    // either compound, single or hybrid prediction as per whatever has
    // worked best for that type of frame in the past.
    // It also predicts whether another coding mode would have worked
    // better than this coding mode. If that is the case, it remembers
    // that for subsequent frames.
    // It does the same analysis for transform size selection also.
    //
    // TODO(zoeliu): To investigate whether a frame_type of LAST2_FRAME needs to
    // be analyzed here to decide on the reference mode.
    const MV_REFERENCE_FRAME frame_type = get_frame_type(cpi);
    int64_t *const mode_thrs = rd_opt->prediction_type_threshes[frame_type];
    int64_t *const filter_thrs = rd_opt->filter_threshes[frame_type];
    int *const tx_thrs = rd_opt->tx_select_threshes[frame_type];
    const int is_alt_ref = frame_type == ALTREF_FRAME;

    /* prediction (compound, single or hybrid) mode selection */
    if (is_alt_ref || !cm->allow_comp_inter_inter)
      cm->reference_mode = SINGLE_REFERENCE;
    else if (mode_thrs[COMPOUND_REFERENCE] > mode_thrs[SINGLE_REFERENCE] &&
             mode_thrs[COMPOUND_REFERENCE] >
                 mode_thrs[REFERENCE_MODE_SELECT] &&
             check_dual_ref_flags(cpi) &&
             cpi->static_mb_pct == 100)
      cm->reference_mode = COMPOUND_REFERENCE;
    else if (mode_thrs[SINGLE_REFERENCE] > mode_thrs[REFERENCE_MODE_SELECT])
      cm->reference_mode = SINGLE_REFERENCE;
    else
      cm->reference_mode = REFERENCE_MODE_SELECT;

    if (cm->interp_filter == SWITCHABLE)
      cm->interp_filter = get_interp_filter(filter_thrs, is_alt_ref);

    encode_frame_internal(cpi);

    for (i = 0; i < REFERENCE_MODES; ++i)
      mode_thrs[i] = (mode_thrs[i] + rd_opt->comp_pred_diff[i] / cm->MBs) / 2;

    for (i = 0; i < SWITCHABLE_FILTER_CONTEXTS; ++i)
      filter_thrs[i] = (filter_thrs[i] + rd_opt->filter_diff[i] / cm->MBs) / 2;

    for (i = 0; i < TX_MODES; ++i) {
      int64_t pd = rd_opt->tx_select_diff[i];
      if (i == TX_MODE_SELECT)
        pd -= RDCOST(cpi->mb.rdmult, cpi->mb.rddiv, 2048 * (TX_SIZES - 1), 0);
      tx_thrs[i] = (tx_thrs[i] + (int)(pd / cm->MBs)) / 2;
    }

    if (cm->reference_mode == REFERENCE_MODE_SELECT) {
      int single_count_zero = 0;
      int comp_count_zero = 0;

      for (i = 0; i < COMP_INTER_CONTEXTS; i++) {
        single_count_zero += cm->counts.comp_inter[i][0];
        comp_count_zero += cm->counts.comp_inter[i][1];
      }

      if (comp_count_zero == 0) {
        cm->reference_mode = SINGLE_REFERENCE;
        vp9_zero(cm->counts.comp_inter);
      } else if (single_count_zero == 0) {
        cm->reference_mode = COMPOUND_REFERENCE;
        vp9_zero(cm->counts.comp_inter);
      }
    }

#if CONFIG_TX64X64
    if (cm->tx_mode == TX_MODE_SELECT) {
      int count4x4_lp = 0;
      int count8x8_8x8p = 0, count8x8_lp = 0;
      int count16x16_16x16p = 0, count16x16_lp = 0;
      int count32x32_32x32p = 0, count32x32_lp = 0;
      int count64x64_64x64p = 0;

      for (i = 0; i < TX_SIZE_CONTEXTS; ++i) {
#if CONFIG_SR_MODE
        // Since in SR mode, tx_size is not sent.
        // To decide ALLOW_TX, here is the real tx_size statistics,
        // but not the statistics for entropy coding
        count4x4_lp += cm->counts.tx.real_p64x64[i][TX_4X4];
        count4x4_lp += cm->counts.tx.real_p32x32[i][TX_4X4];
        count4x4_lp += cm->counts.tx.real_p16x16[i][TX_4X4];
        count4x4_lp += cm->counts.tx.real_p8x8[i][TX_4X4];

        count8x8_lp += cm->counts.tx.real_p64x64[i][TX_8X8];
        count8x8_lp += cm->counts.tx.real_p32x32[i][TX_8X8];
        count8x8_lp += cm->counts.tx.real_p16x16[i][TX_8X8];
        count8x8_8x8p += cm->counts.tx.real_p8x8[i][TX_8X8];

        count16x16_lp += cm->counts.tx.real_p64x64[i][TX_16X16];
        count16x16_lp += cm->counts.tx.real_p32x32[i][TX_16X16];
        count16x16_16x16p += cm->counts.tx.real_p16x16[i][TX_16X16];

        count32x32_lp += cm->counts.tx.real_p64x64[i][TX_32X32];
        count32x32_32x32p += cm->counts.tx.real_p32x32[i][TX_32X32];

        count64x64_64x64p += cm->counts.tx.real_p64x64[i][TX_64X64];
#else  // CONFIG_SR_MODE
        count4x4_lp += cm->counts.tx.p64x64[i][TX_4X4];
        count4x4_lp += cm->counts.tx.p32x32[i][TX_4X4];
        count4x4_lp += cm->counts.tx.p16x16[i][TX_4X4];
        count4x4_lp += cm->counts.tx.p8x8[i][TX_4X4];

        count8x8_lp += cm->counts.tx.p64x64[i][TX_8X8];
        count8x8_lp += cm->counts.tx.p32x32[i][TX_8X8];
        count8x8_lp += cm->counts.tx.p16x16[i][TX_8X8];
        count8x8_8x8p += cm->counts.tx.p8x8[i][TX_8X8];

        count16x16_lp += cm->counts.tx.p64x64[i][TX_16X16];
        count16x16_lp += cm->counts.tx.p32x32[i][TX_16X16];
        count16x16_16x16p += cm->counts.tx.p16x16[i][TX_16X16];

        count32x32_lp += cm->counts.tx.p64x64[i][TX_32X32];
        count32x32_32x32p += cm->counts.tx.p32x32[i][TX_32X32];

        count64x64_64x64p += cm->counts.tx.p64x64[i][TX_64X64];
#endif  // CONFIG_SR_MODE
      }

      if (count4x4_lp == 0 && count16x16_lp == 0 && count16x16_16x16p == 0 &&
          count32x32_lp == 0 && count32x32_32x32p == 0 &&
#if CONFIG_SUPERTX
          cm->counts.supertx_size[TX_16X16] == 0 &&
          cm->counts.supertx_size[TX_32X32] == 0 &&
          cm->counts.supertx_size[TX_64X64] == 0 &&
#endif
          count64x64_64x64p == 0) {
        cm->tx_mode = ALLOW_8X8;
        reset_skip_tx_size(cm, TX_8X8);
      } else if (count8x8_8x8p == 0 && count8x8_lp == 0 &&
                 count16x16_16x16p == 0 && count16x16_lp == 0 &&
                 count32x32_32x32p == 0 && count32x32_lp == 0 &&
#if CONFIG_SUPERTX
                 cm->counts.supertx_size[TX_8X8] == 0 &&
                 cm->counts.supertx_size[TX_16X16] == 0 &&
                 cm->counts.supertx_size[TX_32X32] == 0 &&
                 cm->counts.supertx_size[TX_64X64] == 0 &&
#endif
                 count64x64_64x64p == 0) {
        cm->tx_mode = ONLY_4X4;
        reset_skip_tx_size(cm, TX_4X4);
      } else if (count4x4_lp == 0 && count8x8_lp == 0 && count16x16_lp == 0 &&
                 count32x32_lp == 0) {
        cm->tx_mode = ALLOW_64X64;
      } else if (count4x4_lp == 0 && count8x8_lp == 0 && count16x16_lp == 0 &&
#if CONFIG_SUPERTX
                 cm->counts.supertx_size[TX_64X64] == 0 &&
#endif
                 count64x64_64x64p == 0) {
        cm->tx_mode = ALLOW_32X32;
        reset_skip_tx_size(cm, TX_32X32);
      } else if (count4x4_lp == 0 && count8x8_lp == 0 &&
                 count32x32_lp == 0 && count32x32_32x32p == 0 &&
#if CONFIG_SUPERTX
                 cm->counts.supertx_size[TX_32X32] == 0 &&
                 cm->counts.supertx_size[TX_64X64] == 0 &&
#endif
                 count64x64_64x64p == 0) {
        cm->tx_mode = ALLOW_16X16;
        reset_skip_tx_size(cm, TX_16X16);
      }
    }
#else
    if (cm->tx_mode == TX_MODE_SELECT) {
      int count4x4_lp = 0;
      int count8x8_8x8p = 0, count8x8_lp = 0;
      int count16x16_16x16p = 0, count16x16_lp = 0;
      int count32x32_32x32p = 0;

      for (i = 0; i < TX_SIZE_CONTEXTS; ++i) {
#if CONFIG_SR_MODE
        count4x4_lp += cm->counts.tx.real_p32x32[i][TX_4X4];
        count4x4_lp += cm->counts.tx.real_p16x16[i][TX_4X4];
        count4x4_lp += cm->counts.tx.real_p8x8[i][TX_4X4];

        count8x8_lp += cm->counts.tx.real_p32x32[i][TX_8X8];
        count8x8_lp += cm->counts.tx.real_p16x16[i][TX_8X8];
        count8x8_8x8p += cm->counts.tx.real_p8x8[i][TX_8X8];

        count16x16_lp += cm->counts.tx.real_p32x32[i][TX_16X16];
        count16x16_16x16p += cm->counts.tx.real_p16x16[i][TX_16X16];
        count32x32_32x32p += cm->counts.tx.real_p32x32[i][TX_32X32];
#else  // CONFIG_SR_MODE
        count4x4_lp += cm->counts.tx.p32x32[i][TX_4X4];
        count4x4_lp += cm->counts.tx.p16x16[i][TX_4X4];
        count4x4_lp += cm->counts.tx.p8x8[i][TX_4X4];

        count8x8_lp += cm->counts.tx.p32x32[i][TX_8X8];
        count8x8_lp += cm->counts.tx.p16x16[i][TX_8X8];
        count8x8_8x8p += cm->counts.tx.p8x8[i][TX_8X8];

        count16x16_lp += cm->counts.tx.p32x32[i][TX_16X16];
        count16x16_16x16p += cm->counts.tx.p16x16[i][TX_16X16];
        count32x32_32x32p += cm->counts.tx.p32x32[i][TX_32X32];
#endif  // CONFIG_SR_MODE
      }

      if (count4x4_lp == 0 && count16x16_lp == 0 && count16x16_16x16p == 0 &&
#if CONFIG_SUPERTX
          cm->counts.supertx_size[TX_16X16] == 0 &&
          cm->counts.supertx_size[TX_32X32] == 0 &&
#endif
          count32x32_32x32p == 0) {
        cm->tx_mode = ALLOW_8X8;
        reset_skip_tx_size(cm, TX_8X8);
      } else if (count8x8_8x8p == 0 && count16x16_16x16p == 0 &&
                 count8x8_lp == 0 && count16x16_lp == 0 &&
#if CONFIG_SUPERTX
                 cm->counts.supertx_size[TX_8X8] == 0 &&
                 cm->counts.supertx_size[TX_16X16] == 0 &&
                 cm->counts.supertx_size[TX_32X32] == 0 &&
#endif
                 count32x32_32x32p == 0) {
        cm->tx_mode = ONLY_4X4;
        reset_skip_tx_size(cm, TX_4X4);
      } else if (count8x8_lp == 0 && count16x16_lp == 0 &&
                 count4x4_lp == 0) {
        cm->tx_mode = ALLOW_32X32;
      } else if (count32x32_32x32p == 0 && count8x8_lp == 0 &&
#if CONFIG_SUPERTX
                 cm->counts.supertx_size[TX_32X32] == 0 &&
#endif
                 count4x4_lp == 0) {
        cm->tx_mode = ALLOW_16X16;
        reset_skip_tx_size(cm, TX_16X16);
      }
    }
#endif  // CONFIG_TX64X64
  } else {
    cm->reference_mode = SINGLE_REFERENCE;
    encode_frame_internal(cpi);
  }
}

static void sum_intra_stats(FRAME_COUNTS *counts,
#if CONFIG_FILTERINTRA
  const MACROBLOCKD* xd,
#endif
  const MODE_INFO *mi) {
  const PREDICTION_MODE y_mode = mi->mbmi.mode;
  const PREDICTION_MODE uv_mode = mi->mbmi.uv_mode;
  const BLOCK_SIZE bsize = mi->mbmi.sb_type;
#if CONFIG_FILTERINTRA
  const int uv_fbit = mi->mbmi.uv_filterbit;
  int fbit = mi->mbmi.filterbit;
#endif

  if (bsize < BLOCK_8X8) {
    int idx, idy;
    const int num_4x4_w = num_4x4_blocks_wide_lookup[bsize];
    const int num_4x4_h = num_4x4_blocks_high_lookup[bsize];
    for (idy = 0; idy < 2; idy += num_4x4_h)
      for (idx = 0; idx < 2; idx += num_4x4_w)
#if CONFIG_FILTERINTRA
      {
#endif
        ++counts->y_mode[0][mi->bmi[idy * 2 + idx].as_mode];
#if CONFIG_FILTERINTRA
        if (is_filter_allowed(mi->bmi[idy * 2 + idx].as_mode)) {
          fbit = mi->b_filter_info[idy * 2 + idx];
          ++counts->filterintra[0][mi->bmi[idy * 2 + idx].as_mode][fbit];
        }
      }
#endif
  } else {
    ++counts->y_mode[size_group_lookup[bsize]][y_mode];
#if CONFIG_FILTERINTRA
    if (is_filter_allowed(y_mode) && is_filter_enabled(mi->mbmi.tx_size))
      ++counts->filterintra[mi->mbmi.tx_size][y_mode][fbit];
#endif
  }

  ++counts->uv_mode[y_mode][uv_mode];
#if CONFIG_FILTERINTRA
  if (is_filter_allowed(uv_mode) &&
      is_filter_enabled(get_uv_tx_size(&(mi->mbmi), &xd->plane[1])))
    ++counts->filterintra[get_uv_tx_size(&(mi->mbmi),
                          &xd->plane[1])][uv_mode][uv_fbit];
#endif
}

static void encode_superblock(VP9_COMP *cpi, TOKENEXTRA **t, int output_enabled,
                              int mi_row, int mi_col, BLOCK_SIZE bsize,
                              PICK_MODE_CONTEXT *ctx) {
  VP9_COMMON *const cm = &cpi->common;
  MACROBLOCK *const x = &cpi->mb;
  MACROBLOCKD *const xd = &x->e_mbd;
  MODE_INFO *mi_8x8 = xd->mi;
  MODE_INFO *mi = mi_8x8;
  MB_MODE_INFO *mbmi = &mi->mbmi;
  const int seg_skip = vp9_segfeature_active(&cm->seg, mbmi->segment_id,
                                             SEG_LVL_SKIP);
  const int mis = cm->mi_stride;
  const int mi_width = num_8x8_blocks_wide_lookup[bsize];
  const int mi_height = num_8x8_blocks_high_lookup[bsize];

  x->skip_recode = !x->select_tx_size && mbmi->sb_type >= BLOCK_8X8 &&
                   cpi->oxcf.aq_mode != COMPLEXITY_AQ &&
                   cpi->oxcf.aq_mode != CYCLIC_REFRESH_AQ &&
                   cpi->sf.allow_skip_recode;

  if (!x->skip_recode)
    vpx_memset(x->skip_txfm, 0, sizeof(x->skip_txfm));

  x->skip_optimize = ctx->is_coded;
  ctx->is_coded = 1;
  x->use_lp32x32fdct = cpi->sf.use_lp32x32fdct;
  x->skip_encode = (!output_enabled && cpi->sf.skip_encode_frame &&
                    x->q_index < QIDX_SKIP_THRESH);

  if (x->skip_encode)
    return;

  set_ref_ptrs(cm, xd, mbmi->ref_frame[0], mbmi->ref_frame[1]);

  if (!is_inter_block(mbmi)
#if CONFIG_INTRABC
      && !is_intrabc_mode(mbmi->mode)
#endif  // CONFIG_INTRABC
      ) {
    int plane;
#if CONFIG_MISC_ENTROPY
    mbmi->skip = 0;
#else
    mbmi->skip = 1;
#endif
    for (plane = 0; plane < MAX_MB_PLANE; ++plane)
      vp9_encode_intra_block_plane(x, MAX(bsize, BLOCK_8X8), plane, 1);
    if (output_enabled)
      sum_intra_stats(&cm->counts,
#if CONFIG_FILTERINTRA
        xd,
#endif
        mi);
    vp9_tokenize_sb(cpi, t, !output_enabled, MAX(bsize, BLOCK_8X8));
#if CONFIG_PALETTE
    if (bsize >= BLOCK_8X8 && output_enabled) {
      if (mbmi->palette_enabled[0]) {
        int rows = 4 * num_4x4_blocks_high_lookup[bsize];
        int cols = 4 * num_4x4_blocks_wide_lookup[bsize];

        vp9_palette_color_insertion(cm->current_palette_colors,
                                    &cm ->current_palette_size,
                                    cm->current_palette_count, mbmi);
        CHECK_MEM_ERROR(cm, mbmi->palette_color_map,
                        vpx_memalign(16, rows * cols *
                                     sizeof(xd->plane[0].color_index_map[0])));
        memcpy(mbmi->palette_color_map, xd->plane[0].color_index_map,
               rows * cols * sizeof(xd->plane[0].color_index_map[0]));
      }

      if (mbmi->palette_enabled[1]) {
        int rows = 4 * num_4x4_blocks_high_lookup[bsize] >>
            xd->plane[1].subsampling_y;
        int cols = 4 * num_4x4_blocks_wide_lookup[bsize] >>
            xd->plane[1].subsampling_x;

        CHECK_MEM_ERROR(cm, mbmi->palette_uv_color_map,
                        vpx_memalign(16, rows * cols *
                                     sizeof(xd->plane[1].color_index_map[0])));
        memcpy(mbmi->palette_uv_color_map, xd->plane[1].color_index_map,
               rows * cols * sizeof(xd->plane[1].color_index_map[0]));
      }
    }
    if (frame_is_intra_only(cm) && output_enabled && bsize >= BLOCK_8X8) {
      cm->palette_blocks_signalled++;
      if (mbmi->palette_enabled[0])
        cm->palette_counter++;
    }
#endif  // CONFIG_PALETTE
  } else {
    int ref;
    const int is_compound = has_second_ref(mbmi);
    for (ref = 0; ref < 1 + is_compound; ++ref) {
      YV12_BUFFER_CONFIG *cfg = get_ref_frame_buffer(cpi,
                                                     mbmi->ref_frame[ref]);
      vp9_setup_pre_planes(xd, ref, cfg, mi_row, mi_col,
                           &xd->block_refs[ref]->sf);
    }
    if (!(cpi->sf.reuse_inter_pred_sby && ctx->pred_pixel_ready) || seg_skip)
      vp9_build_inter_predictors_sby(xd, mi_row, mi_col, MAX(bsize, BLOCK_8X8));

    vp9_build_inter_predictors_sbuv(xd, mi_row, mi_col, MAX(bsize, BLOCK_8X8));
    vp9_encode_sb(x, MAX(bsize, BLOCK_8X8));
    vp9_tokenize_sb(cpi, t, !output_enabled, MAX(bsize, BLOCK_8X8));
  }
#if CONFIG_NEW_QUANT && QUANT_PROFILES > 1 && !Q_CTX_BASED_PROFILES
  // This is not strictly required, but is a good practice.
  // If you remove this, the assert in vp9_bitstream.c needs to be removed also.
  if (mbmi->skip)
    mbmi->dq_off_index = 0;
#endif  // CONFIG_NEW_QUANT && QUANT_PROFILES > 1 && !Q_CTX_BASED_PROFILES

#if CONFIG_INTRABC
  if (frame_is_intra_only(cm) && output_enabled && bsize >= BLOCK_8X8) {
    cm->intrabc_blocks_signalled++;
    if (is_intrabc_mode(mbmi->mode))
      cm->intrabc_counter++;
  }
#endif  // CONFIG_INTRABC

  if (output_enabled) {
#if CONFIG_SR_MODE
    assert(bsize == mbmi->sb_type);
    if (is_enable_srmode(bsize) &&
        !(is_inter_block(mbmi) && (mbmi->skip || seg_skip))) {
      cm->counts.sr[vp9_get_sr_context(xd, bsize)][mbmi->sr]++;
#if SR_USE_MULTI_F
      if (mbmi->sr)
        cm->counts.sr_usfilters[vp9_get_sr_usfilter_context(xd)]
                               [mbmi->us_filter_idx]++;
#endif  // SR_USE_MULTI_F
    }

    if (cm->tx_mode == TX_MODE_SELECT &&
        mbmi->sb_type >= BLOCK_8X8  &&
        !(is_inter_block(mbmi) && (mbmi->skip || seg_skip))) {
      ++get_real_tx_counts(max_txsize_lookup[bsize],
                           vp9_get_tx_size_context(xd),
                           &cm->counts.tx)[mbmi->tx_size];
    }
#endif  // CONFIG_SR_MODE
    if (cm->tx_mode == TX_MODE_SELECT &&
        mbmi->sb_type >= BLOCK_8X8  &&
#if CONFIG_SR_MODE
        !mbmi->sr &&
#endif  // CONFIG_SR_MODE
        !(is_inter_block(mbmi) && (mbmi->skip || seg_skip))) {
      ++get_tx_counts(max_txsize_lookup[bsize], vp9_get_tx_size_context(xd),
                      &cm->counts.tx)[mbmi->tx_size];
    } else {
      int x, y;
      TX_SIZE tx_size;
      // The new intra coding scheme requires no change of transform size
      if (is_inter_block(&mi->mbmi)) {
        tx_size = MIN(tx_mode_to_biggest_tx_size[cm->tx_mode],
                      max_txsize_lookup[bsize]);
      } else {
        tx_size = (bsize >= BLOCK_8X8) ? mbmi->tx_size : TX_4X4;
      }

      for (y = 0; y < mi_height; y++)
        for (x = 0; x < mi_width; x++)
          if (mi_col + x < cm->mi_cols && mi_row + y < cm->mi_rows)
            mi_8x8[mis * y + x].src_mi->mbmi.tx_size = tx_size;
    }
#if CONFIG_EXT_TX
    if (is_inter_block(mbmi) &&
#if !CONFIG_WAVELETS
        mbmi->tx_size <= TX_16X16 &&
#endif
        cm->base_qindex > 0 &&
        bsize >= BLOCK_8X8 &&
        !mbmi->skip &&
        !vp9_segfeature_active(&cm->seg, mbmi->segment_id, SEG_LVL_SKIP)) {
      ++cm->counts.ext_tx[mbmi->tx_size][mbmi->ext_txfrm];
    }
#endif  // CONFIG_EXT_TX
#if CONFIG_NEW_QUANT && QUANT_PROFILES > 1 && !Q_CTX_BASED_PROFILES
    if (cm->base_qindex > Q_THRESHOLD_MIN &&
        cm->base_qindex < Q_THRESHOLD_MAX &&
        mbmi->send_dq_bit &&
#if CONFIG_COPY_MODE
        (frame_is_intra_only(cm) || mbmi->copy_mode == NOREF) &&
#endif  // CONFIG_COPY_MODE
        !mbmi->skip &&
        !vp9_segfeature_active(&cm->seg, mbmi->segment_id, SEG_LVL_SKIP)) {
      ++cm->counts.dq_profile[mbmi->dq_off_index];
    }
#endif  // CONFIG_NEW_QUANT && QUANT_PROFILES > 1 && !Q_CTX_BASED_PROFILES
#if CONFIG_TX_SKIP
    if (bsize >= BLOCK_8X8) {
      int q_idx = vp9_get_qindex(&cm->seg, mbmi->segment_id, cm->base_qindex);
      int try_tx_skip = is_inter_block(mbmi) ? q_idx <= tx_skip_q_thresh_inter :
                                               q_idx <= tx_skip_q_thresh_intra;
#if CONFIG_COPY_MODE
      if (mbmi->copy_mode != NOREF)
        try_tx_skip = 0;
#endif  // CONFIG_COPY_MODE

#if CONFIG_SUPERTX
      if (try_tx_skip) {
#else
      if (try_tx_skip && (!(mbmi->skip || seg_skip) || !is_inter_block(mbmi))) {
#endif  // CONFIG_SUPERTX
        ++cm->counts.y_tx_skip[is_inter_block(mbmi)][mbmi->tx_skip[0]];
        ++cm->counts.uv_tx_skip[mbmi->tx_skip[0]][mbmi->tx_skip[1]];
      }
    }
#endif  // CONFIG_TX_SKIP
#if CONFIG_PALETTE
      if (!frame_is_intra_only(cm) && !is_inter_block(mbmi) &&
          bsize >= BLOCK_8X8 && cm->allow_palette_mode) {
        int palette_ctx = 0;
        const MODE_INFO *above_mi = xd->up_available ?
            xd->mi[-xd->mi_stride].src_mi : NULL;
        const MODE_INFO *left_mi = xd->left_available ?
            xd->mi[-1].src_mi : NULL;
        if (above_mi)
          palette_ctx += (above_mi->mbmi.palette_enabled[0] == 1);
        if (left_mi)
          palette_ctx += (left_mi->mbmi.palette_enabled[0] == 1);
        vp9_update_palette_counts(&cm->counts, mbmi, bsize, palette_ctx);
      }
#endif  // CONFIG_PALETTE
  }
}

#if CONFIG_SUPERTX
static int check_intra_b(PICK_MODE_CONTEXT *ctx) {
#if CONFIG_NEW_INTER
#if CONFIG_INTERINTRA
  return (!is_inter_mode((&ctx->mic)->mbmi.mode) &&
          !is_inter_compound_mode((&ctx->mic)->mbmi.mode)) ||
         (ctx->mic.mbmi.ref_frame[1] == INTRA_FRAME);
#else
  return !is_inter_mode((&ctx->mic)->mbmi.mode) &&
         !is_inter_compound_mode((&ctx->mic)->mbmi.mode);
#endif  // CONFIG_INTERINTRA
#else   // CONFIG_NEW_INTER
#if CONFIG_INTERINTRA
  return !is_inter_mode((&ctx->mic)->mbmi.mode) ||
         (ctx->mic.mbmi.ref_frame[1] == INTRA_FRAME);
#else
  return !is_inter_mode((&ctx->mic)->mbmi.mode);
#endif  // CONFIG_INTERINTRA
#endif  // CONFIG_NEW_INTER
}

static int check_intra_sb(VP9_COMP *cpi, const TileInfo *const tile,
                          int mi_row, int mi_col, BLOCK_SIZE bsize,
                          PC_TREE *pc_tree) {
  VP9_COMMON *const cm = &cpi->common;

  const int bsl = b_width_log2_lookup[bsize], hbs = (1 << bsl) / 4;
  PARTITION_TYPE partition;
  BLOCK_SIZE subsize = bsize;
#if CONFIG_EXT_PARTITION
  int i;
#endif

  if (mi_row >= cm->mi_rows || mi_col >= cm->mi_cols)
    return 1;

  if (bsize >= BLOCK_8X8)
    subsize = get_subsize(bsize, pc_tree->partitioning);
  else
    subsize = BLOCK_4X4;

  partition = partition_lookup[bsl][subsize];
#if CONFIG_EXT_PARTITION
  if (bsize > BLOCK_8X8)
    partition = pc_tree->partitioning;
#endif

  switch (partition) {
    case PARTITION_NONE:
      return check_intra_b(&pc_tree->none);
      break;
    case PARTITION_VERT:
      if (check_intra_b(&pc_tree->vertical[0]))
        return 1;
      if (mi_col + hbs < cm->mi_cols && bsize > BLOCK_8X8) {
        if (check_intra_b(&pc_tree->vertical[1]))
          return 1;
      }
      break;
    case PARTITION_HORZ:
      if (check_intra_b(&pc_tree->horizontal[0]))
        return 1;
      if (mi_row + hbs < cm->mi_rows && bsize > BLOCK_8X8) {
        if (check_intra_b(&pc_tree->horizontal[1]))
          return 1;
      }
      break;
    case PARTITION_SPLIT:
      if (bsize == BLOCK_8X8) {
        if (check_intra_b(pc_tree->leaf_split[0]))
          return 1;
      } else {
        if (check_intra_sb(cpi, tile, mi_row, mi_col, subsize,
                           pc_tree->split[0]))
          return 1;
        if (check_intra_sb(cpi, tile, mi_row, mi_col + hbs, subsize,
                           pc_tree->split[1]))
          return 1;
        if (check_intra_sb(cpi, tile, mi_row + hbs, mi_col, subsize,
                           pc_tree->split[2]))
          return 1;
        if (check_intra_sb(cpi, tile, mi_row + hbs, mi_col + hbs, subsize,
                           pc_tree->split[3]))
          return 1;
      }
      break;
#if CONFIG_EXT_PARTITION
    case PARTITION_HORZ_A:
      for (i = 0; i < 3; i++) {
        if (check_intra_b(&pc_tree->horizontala[i]))
          return 1;
      }
      break;
    case PARTITION_HORZ_B:
      for (i = 0; i < 3; i++) {
        if (check_intra_b(&pc_tree->horizontalb[i]))
          return 1;
      }
      break;
    case PARTITION_VERT_A:
      for (i = 0; i < 3; i++) {
        if (check_intra_b(&pc_tree->verticala[i]))
          return 1;
      }
      break;
    case PARTITION_VERT_B:
      for (i = 0; i < 3; i++) {
        if (check_intra_b(&pc_tree->verticalb[i]))
          return 1;
      }
      break;
#endif
    default:
      assert(0);
  }
  return 0;
}

static int check_supertx_b(TX_SIZE supertx_size, PICK_MODE_CONTEXT *ctx) {
  return ctx->mic.mbmi.tx_size == supertx_size;
}

static int check_supertx_sb(BLOCK_SIZE bsize, TX_SIZE supertx_size,
                            PC_TREE *pc_tree) {
  PARTITION_TYPE partition;
  BLOCK_SIZE subsize;

  partition = pc_tree->partitioning;
  subsize = get_subsize(bsize, partition);
  switch (partition) {
    case PARTITION_NONE:
      return check_supertx_b(supertx_size, &pc_tree->none);
    case PARTITION_VERT:
      return check_supertx_b(supertx_size, &pc_tree->vertical[0]);
    case PARTITION_HORZ:
      return check_supertx_b(supertx_size, &pc_tree->horizontal[0]);
    case PARTITION_SPLIT:
      if (bsize == BLOCK_8X8)
        return check_supertx_b(supertx_size, pc_tree->leaf_split[0]);
      else
        return check_supertx_sb(subsize, supertx_size, pc_tree->split[0]);
#if CONFIG_EXT_PARTITION
    case PARTITION_HORZ_A:
      return check_supertx_b(supertx_size, &pc_tree->horizontala[0]);
    case PARTITION_HORZ_B:
      return check_supertx_b(supertx_size, &pc_tree->horizontalb[0]);
    case PARTITION_VERT_A:
      return check_supertx_b(supertx_size, &pc_tree->verticala[0]);
    case PARTITION_VERT_B:
      return check_supertx_b(supertx_size, &pc_tree->verticalb[0]);
#endif
    default:
      assert(0);
      return 0;
  }
}

static void predict_superblock(VP9_COMP *cpi,
#if CONFIG_WEDGE_PARTITION
                               int mi_row_ori, int mi_col_ori,
#endif  // CONFIG_WEDGE_PARTITION
                               int mi_row_pred, int mi_col_pred,
                               BLOCK_SIZE bsize_pred, int b_sub8x8, int block) {
  // Used in supertx
  // (mi_row_ori, mi_col_ori): location for mv
  // (mi_row_pred, mi_col_pred, bsize_pred): region to predict
  VP9_COMMON *const cm = &cpi->common;
  MACROBLOCK *const x = &cpi->mb;
  MACROBLOCKD *const xd = &x->e_mbd;
  MODE_INFO *mi_8x8 = xd->mi;
  MODE_INFO *mi = mi_8x8;
  MB_MODE_INFO *mbmi = &mi->mbmi;
  int ref;
  const int is_compound = has_second_ref(mbmi);

  set_ref_ptrs(cm, xd, mbmi->ref_frame[0], mbmi->ref_frame[1]);

  for (ref = 0; ref < 1 + is_compound; ++ref) {
    YV12_BUFFER_CONFIG *cfg = get_ref_frame_buffer(cpi,
                                                   mbmi->ref_frame[ref]);
    vp9_setup_pre_planes(xd, ref, cfg, mi_row_pred, mi_col_pred,
                         &xd->block_refs[ref]->sf);
  }

#if !CONFIG_WEDGE_PARTITION
  if (!b_sub8x8)
    vp9_build_inter_predictors_sb(xd, mi_row_pred, mi_col_pred, bsize_pred);
  else
    vp9_build_inter_predictors_sb_sub8x8(xd, mi_row_pred, mi_col_pred,
                                         bsize_pred, block);
#else
  if (!b_sub8x8)
    vp9_build_inter_predictors_sb_extend(xd, mi_row_ori, mi_col_ori,
                                         mi_row_pred, mi_col_pred, bsize_pred);
  else
    vp9_build_inter_predictors_sb_sub8x8_extend(
        xd, mi_row_ori, mi_col_ori,
        mi_row_pred, mi_col_pred, bsize_pred, block);
#endif  // CONFIG_WEDGE_PARTITION
}

static void predict_b_extend(VP9_COMP *cpi, const TileInfo *const tile,
                             int block,
                             int mi_row_ori, int mi_col_ori,
                             int mi_row_pred, int mi_col_pred,
                             int mi_row_top, int mi_col_top,
                             uint8_t * dst_buf[3], int dst_stride[3],
                             BLOCK_SIZE bsize_ori, BLOCK_SIZE bsize_top,
                             BLOCK_SIZE bsize_pred, int output_enabled,
                             int b_sub8x8, int bextend) {
  // Used in supertx
  // (mi_row_ori, mi_col_ori): location for mv
  // (mi_row_pred, mi_col_pred, bsize_pred): region to predict
  // (mi_row_top, mi_col_top, bsize_top): region of the top partition size
  // block: sub location of sub8x8 blocks
  // b_sub8x8: 1: ori is sub8x8; 0: ori is not sub8x8
  // bextend: 1: region to predict is an extension of ori; 0: not

  MACROBLOCK *const x = &cpi->mb;
  VP9_COMMON *const cm = &cpi->common;
  MACROBLOCKD *const xd = &x->e_mbd;
  int r = (mi_row_pred - mi_row_top) * MI_SIZE;
  int c = (mi_col_pred - mi_col_top) * MI_SIZE;
  const int mi_width_top = num_8x8_blocks_wide_lookup[bsize_top];
  const int mi_height_top = num_8x8_blocks_high_lookup[bsize_top];

  if (mi_row_pred < mi_row_top || mi_col_pred < mi_col_top ||
      mi_row_pred >= mi_row_top + mi_height_top ||
      mi_col_pred >= mi_col_top + mi_width_top ||
      mi_row_pred >= cm->mi_rows || mi_col_pred >= cm->mi_cols)
    return;

  set_offsets_extend(cpi, tile, mi_row_pred, mi_col_pred,
                     mi_row_ori, mi_col_ori, bsize_pred, bsize_ori);
  xd->plane[0].dst.stride = dst_stride[0];
  xd->plane[1].dst.stride = dst_stride[1];
  xd->plane[2].dst.stride = dst_stride[2];
  xd->plane[0].dst.buf = dst_buf[0] +
                         (r >> xd->plane[0].subsampling_y) * dst_stride[0] +
                         (c >> xd->plane[0].subsampling_x);
  xd->plane[1].dst.buf = dst_buf[1] +
                         (r >> xd->plane[1].subsampling_y) * dst_stride[1] +
                         (c >> xd->plane[1].subsampling_x);
  xd->plane[2].dst.buf = dst_buf[2] +
                         (r >> xd->plane[2].subsampling_y) * dst_stride[2] +
                         (c >> xd->plane[2].subsampling_x);

  predict_superblock(cpi,
#if CONFIG_WEDGE_PARTITION
                     mi_row_ori, mi_col_ori,
#endif
                     mi_row_pred, mi_col_pred, bsize_pred,
                     b_sub8x8, block);

  if (output_enabled && !bextend)
    update_stats(&cpi->common, &cpi->mb);
}

static void extend_dir(VP9_COMP *cpi, const TileInfo *const tile,
                       int block, BLOCK_SIZE bsize, BLOCK_SIZE top_bsize,
                       int mi_row, int mi_col,
                       int mi_row_top, int mi_col_top,
                       int output_enabled,
                       uint8_t * dst_buf[3], int dst_stride[3], int dir) {
  // dir: 0-lower, 1-upper, 2-left, 3-right
  //      4-lowerleft, 5-upperleft, 6-lowerright, 7-upperright
  MACROBLOCKD *xd = &cpi->mb.e_mbd;
  const int mi_width = num_8x8_blocks_wide_lookup[bsize];
  const int mi_height = num_8x8_blocks_high_lookup[bsize];
  int xss = xd->plane[1].subsampling_x;
  int yss = xd->plane[1].subsampling_y;
  int b_sub8x8 = (bsize < BLOCK_8X8) ? 1 : 0;

  BLOCK_SIZE extend_bsize;
  int unit, mi_row_pred, mi_col_pred;

  if (dir == 0 || dir == 1) {  // lower and upper
    extend_bsize = (mi_width == 1 || bsize < BLOCK_8X8 || xss < yss) ?
                   BLOCK_8X8 : BLOCK_16X8;
    unit = num_8x8_blocks_wide_lookup[extend_bsize];
    mi_row_pred = mi_row + ((dir == 0) ? mi_height : -1);
    mi_col_pred = mi_col;

    predict_b_extend(cpi, tile, block, mi_row, mi_col,
                     mi_row_pred, mi_col_pred,
                     mi_row_top, mi_col_top, dst_buf, dst_stride,
                     bsize, top_bsize, extend_bsize,
                     output_enabled, b_sub8x8, 1);

    if (mi_width > unit) {
      int i;
      for (i = 0; i < mi_width/unit - 1; i++) {
        mi_col_pred += unit;
        predict_b_extend(cpi, tile, block, mi_row, mi_col,
                         mi_row_pred, mi_col_pred, mi_row_top, mi_col_top,
                         dst_buf, dst_stride, bsize, top_bsize, extend_bsize,
                         output_enabled, b_sub8x8, 1);
      }
    }
  } else if (dir == 2 || dir == 3) {  // left and right
    extend_bsize = (mi_height == 1 || bsize < BLOCK_8X8 || yss < xss) ?
                   BLOCK_8X8 : BLOCK_8X16;
    unit = num_8x8_blocks_high_lookup[extend_bsize];
    mi_row_pred = mi_row;
    mi_col_pred = mi_col + ((dir == 3) ? mi_width : -1);

    predict_b_extend(cpi, tile, block, mi_row, mi_col,
                     mi_row_pred, mi_col_pred, mi_row_top, mi_col_top,
                     dst_buf, dst_stride, bsize, top_bsize, extend_bsize,
                     output_enabled, b_sub8x8, 1);

    if (mi_height > unit) {
      int i;
      for (i = 0; i < mi_height/unit - 1; i++) {
        mi_row_pred += unit;
        predict_b_extend(cpi, tile, block, mi_row, mi_col,
                         mi_row_pred, mi_col_pred, mi_row_top, mi_col_top,
                         dst_buf, dst_stride, bsize, top_bsize, extend_bsize,
                         output_enabled, b_sub8x8, 1);
      }
    }
  } else {
    extend_bsize = BLOCK_8X8;
    mi_row_pred = mi_row + ((dir == 4 || dir == 6) ? mi_height : -1);
    mi_col_pred = mi_col + ((dir == 6 || dir == 7) ? mi_width : -1);

    predict_b_extend(cpi, tile, block, mi_row, mi_col,
                     mi_row_pred, mi_col_pred, mi_row_top, mi_col_top,
                     dst_buf, dst_stride, bsize, top_bsize, extend_bsize,
                     output_enabled, b_sub8x8, 1);
  }
}

static void extend_all(VP9_COMP *cpi, const TileInfo *const tile,
                       int block,
                       BLOCK_SIZE bsize, BLOCK_SIZE top_bsize,
                       int mi_row, int mi_col,
                       int mi_row_top, int mi_col_top,
                       int output_enabled,
                       uint8_t * dst_buf[3], int dst_stride[3]) {
  assert(block >= 0 && block < 4);
  extend_dir(cpi, tile, block, bsize, top_bsize, mi_row, mi_col,
             mi_row_top, mi_col_top, output_enabled, dst_buf, dst_stride, 0);
  extend_dir(cpi, tile, block, bsize, top_bsize, mi_row, mi_col,
             mi_row_top, mi_col_top, output_enabled, dst_buf, dst_stride, 1);
  extend_dir(cpi, tile, block, bsize, top_bsize, mi_row, mi_col,
             mi_row_top, mi_col_top, output_enabled, dst_buf, dst_stride, 2);
  extend_dir(cpi, tile, block, bsize, top_bsize, mi_row, mi_col,
             mi_row_top, mi_col_top, output_enabled, dst_buf, dst_stride, 3);
  extend_dir(cpi, tile, block, bsize, top_bsize, mi_row, mi_col,
             mi_row_top, mi_col_top, output_enabled, dst_buf, dst_stride, 4);
  extend_dir(cpi, tile, block, bsize, top_bsize, mi_row, mi_col,
             mi_row_top, mi_col_top, output_enabled, dst_buf, dst_stride, 5);
  extend_dir(cpi, tile, block, bsize, top_bsize, mi_row, mi_col,
             mi_row_top, mi_col_top, output_enabled, dst_buf, dst_stride, 6);
  extend_dir(cpi, tile, block, bsize, top_bsize, mi_row, mi_col,
             mi_row_top, mi_col_top, output_enabled, dst_buf, dst_stride, 7);
}


// This function generates prediction for multiple blocks, between which
// discontinuity around boundary is reduced by smoothing masks. The basic
// smoothing mask is a soft step function along horz/vert direction. In more
// complicated case when a block is split into 4 subblocks, the basic mask is
// first applied to neighboring subblocks (2 pairs) in horizontal direction and
// then applied to the 2 masked prediction mentioned above in vertical direction
// If the block is split into more than one level, at every stage, masked
// prediction is stored in dst_buf[] passed from higher level.
static void predict_sb_complex(VP9_COMP *cpi, const TileInfo *const tile,
                               int mi_row, int mi_col,
                               int mi_row_top, int mi_col_top,
                               int output_enabled, BLOCK_SIZE bsize,
                               BLOCK_SIZE top_bsize,
                               uint8_t *dst_buf[3], int dst_stride[3],
                               PC_TREE *pc_tree) {
  VP9_COMMON *const cm = &cpi->common;
  MACROBLOCK *const x = &cpi->mb;
  MACROBLOCKD *const xd = &x->e_mbd;

  const int bsl = b_width_log2_lookup[bsize], hbs = (1 << bsl) / 4;
  PARTITION_TYPE partition;
  BLOCK_SIZE subsize;
#if CONFIG_EXT_PARTITION
  BLOCK_SIZE bsize2 = get_subsize(bsize, PARTITION_SPLIT);
#endif

  int i, ctx;
  uint8_t *dst_buf1[3], *dst_buf2[3], *dst_buf3[3];
  DECLARE_ALIGNED_ARRAY(16, uint8_t, tmp_buf1,
                        MAX_MB_PLANE * MAXTXLEN * MAXTXLEN * sizeof(uint16_t));
  DECLARE_ALIGNED_ARRAY(16, uint8_t, tmp_buf2,
                        MAX_MB_PLANE * MAXTXLEN * MAXTXLEN * sizeof(uint16_t));
  DECLARE_ALIGNED_ARRAY(16, uint8_t, tmp_buf3,
                        MAX_MB_PLANE * MAXTXLEN * MAXTXLEN * sizeof(uint16_t));
  int dst_stride1[3] = {MAXTXLEN, MAXTXLEN, MAXTXLEN};
  int dst_stride2[3] = {MAXTXLEN, MAXTXLEN, MAXTXLEN};
  int dst_stride3[3] = {MAXTXLEN, MAXTXLEN, MAXTXLEN};
#if CONFIG_VP9_HIGHBITDEPTH
  if (xd->cur_buf->flags & YV12_FLAG_HIGHBITDEPTH) {
    int len = sizeof(uint16_t);
    dst_buf1[0] = CONVERT_TO_BYTEPTR(tmp_buf1);
    dst_buf1[1] = CONVERT_TO_BYTEPTR(tmp_buf1 + MAXTXLEN * MAXTXLEN * len);
    dst_buf1[2] = CONVERT_TO_BYTEPTR(tmp_buf1 + 2 * MAXTXLEN * MAXTXLEN * len);
    dst_buf2[0] = CONVERT_TO_BYTEPTR(tmp_buf2);
    dst_buf2[1] = CONVERT_TO_BYTEPTR(tmp_buf2 + MAXTXLEN * MAXTXLEN * len);
    dst_buf2[2] = CONVERT_TO_BYTEPTR(tmp_buf2 + 2 * MAXTXLEN * MAXTXLEN * len);
    dst_buf3[0] = CONVERT_TO_BYTEPTR(tmp_buf3);
    dst_buf3[1] = CONVERT_TO_BYTEPTR(tmp_buf3 + MAXTXLEN * MAXTXLEN * len);
    dst_buf3[2] = CONVERT_TO_BYTEPTR(tmp_buf3 + 2 * MAXTXLEN * MAXTXLEN * len);
  } else {
#endif
    dst_buf1[0] = tmp_buf1;
    dst_buf1[1] = tmp_buf1 + MAXTXLEN * MAXTXLEN;
    dst_buf1[2] = tmp_buf1 + 2 * MAXTXLEN * MAXTXLEN;
    dst_buf2[0] = tmp_buf2;
    dst_buf2[1] = tmp_buf2 + MAXTXLEN * MAXTXLEN;
    dst_buf2[2] = tmp_buf2 + 2 * MAXTXLEN * MAXTXLEN;
    dst_buf3[0] = tmp_buf3;
    dst_buf3[1] = tmp_buf3 + MAXTXLEN * MAXTXLEN;
    dst_buf3[2] = tmp_buf3 + 2 * MAXTXLEN * MAXTXLEN;
#if CONFIG_VP9_HIGHBITDEPTH
  }
#endif

  if (mi_row >= cm->mi_rows || mi_col >= cm->mi_cols)
    return;

  if (bsize >= BLOCK_8X8) {
    ctx = partition_plane_context(xd, mi_row, mi_col, bsize);
    subsize = get_subsize(bsize, pc_tree->partitioning);
  } else {
    ctx = 0;
    subsize = BLOCK_4X4;
  }
  partition = partition_lookup[bsl][subsize];
#if CONFIG_EXT_PARTITION
  if (bsize > BLOCK_8X8)
    partition = pc_tree->partitioning;
#endif
  if (output_enabled && bsize != BLOCK_4X4 && bsize < top_bsize)
      cm->counts.partition[ctx][partition]++;

  for (i = 0; i < MAX_MB_PLANE; i++) {
    xd->plane[i].dst.buf = dst_buf[i];
    xd->plane[i].dst.stride = dst_stride[i];
  }

  switch (partition) {
    case PARTITION_NONE:
      assert(bsize < top_bsize);
      predict_b_extend(cpi, tile, 0, mi_row, mi_col, mi_row, mi_col,
                       mi_row_top, mi_col_top, dst_buf, dst_stride,
                       bsize, top_bsize, bsize, output_enabled, 0, 0);
      extend_all(cpi, tile, 0, bsize, top_bsize, mi_row, mi_col,
                 mi_row_top, mi_col_top, output_enabled, dst_buf, dst_stride);
      break;
    case PARTITION_HORZ:
      if (bsize == BLOCK_8X8) {
        // Fisrt half
        predict_b_extend(cpi, tile, 0, mi_row, mi_col, mi_row, mi_col,
                         mi_row_top, mi_col_top, dst_buf, dst_stride,
                         subsize, top_bsize, BLOCK_8X8, output_enabled, 1, 0);
        if (bsize < top_bsize)
          extend_all(cpi, tile, 0, subsize, top_bsize, mi_row, mi_col,
                     mi_row_top, mi_col_top, output_enabled,
                     dst_buf, dst_stride);

        // Second half
        predict_b_extend(cpi, tile, 2, mi_row, mi_col, mi_row, mi_col,
                         mi_row_top, mi_col_top, dst_buf1, dst_stride1,
                         subsize, top_bsize, BLOCK_8X8, output_enabled, 1, 1);
        if (bsize < top_bsize)
          extend_all(cpi, tile, 2, subsize, top_bsize, mi_row, mi_col,
                     mi_row_top, mi_col_top, output_enabled,
                     dst_buf1, dst_stride1);

        // Smooth
        xd->plane[0].dst.buf = dst_buf[0];
        xd->plane[0].dst.stride = dst_stride[0];
        vp9_build_masked_inter_predictor_complex(xd,
                                                 dst_buf[0], dst_stride[0],
                                                 dst_buf1[0], dst_stride1[0],
                                                 &xd->plane[0],
                                                 mi_row, mi_col,
                                                 mi_row_top, mi_col_top,
                                                 bsize, top_bsize,
                                                 PARTITION_HORZ, 0);
      }  else {
        // First half
        predict_b_extend(cpi, tile, 0, mi_row, mi_col, mi_row, mi_col,
                         mi_row_top, mi_col_top, dst_buf, dst_stride,
                         subsize, top_bsize, subsize, output_enabled, 0, 0);
        if (bsize < top_bsize)
          extend_all(cpi, tile, 0, subsize, top_bsize, mi_row, mi_col,
                     mi_row_top, mi_col_top, output_enabled,
                     dst_buf, dst_stride);
        else
          extend_dir(cpi, tile, 0, subsize, top_bsize, mi_row, mi_col,
                     mi_row_top, mi_col_top, output_enabled,
                     dst_buf, dst_stride, 0);

        if (mi_row + hbs < cm->mi_rows) {
          // Second half
          predict_b_extend(cpi, tile, 0, mi_row + hbs, mi_col,
                           mi_row + hbs, mi_col, mi_row_top, mi_col_top,
                           dst_buf1, dst_stride1, subsize, top_bsize, subsize,
                           output_enabled, 0, 0);
          if (bsize < top_bsize)
            extend_all(cpi, tile, 0, subsize, top_bsize, mi_row + hbs, mi_col,
                       mi_row_top, mi_col_top, output_enabled,
                       dst_buf1, dst_stride1);
          else
            extend_dir(cpi, tile, 0, subsize, top_bsize, mi_row + hbs, mi_col,
                       mi_row_top, mi_col_top, output_enabled,
                       dst_buf1, dst_stride1, 1);

          // Smooth
          for (i = 0; i < MAX_MB_PLANE; i++) {
            xd->plane[i].dst.buf = dst_buf[i];
            xd->plane[i].dst.stride = dst_stride[i];
            vp9_build_masked_inter_predictor_complex(
                xd, dst_buf[i], dst_stride[i], dst_buf1[i], dst_stride1[i],
                &xd->plane[i], mi_row, mi_col, mi_row_top, mi_col_top,
                bsize, top_bsize, PARTITION_HORZ, i);
          }
        }
      }
      break;
    case PARTITION_VERT:
      if (bsize == BLOCK_8X8) {
        // First half
        predict_b_extend(cpi, tile, 0, mi_row, mi_col, mi_row, mi_col,
                         mi_row_top, mi_col_top, dst_buf, dst_stride,
                         subsize, top_bsize, BLOCK_8X8, output_enabled, 1, 0);
        if (bsize < top_bsize)
          extend_all(cpi, tile, 0, subsize, top_bsize, mi_row, mi_col,
                     mi_row_top, mi_col_top, output_enabled,
                     dst_buf, dst_stride);

        // Second half
        predict_b_extend(cpi, tile, 1, mi_row, mi_col, mi_row, mi_col,
                         mi_row_top, mi_col_top, dst_buf1, dst_stride1,
                         subsize, top_bsize, BLOCK_8X8, output_enabled, 1, 1);
        if (bsize < top_bsize)
          extend_all(cpi, tile, 1, subsize, top_bsize, mi_row, mi_col,
                     mi_row_top, mi_col_top, output_enabled,
                     dst_buf1, dst_stride1);

        // Smooth
        xd->plane[0].dst.buf = dst_buf[0];
        xd->plane[0].dst.stride = dst_stride[0];
        vp9_build_masked_inter_predictor_complex(xd,
                                                 dst_buf[0], dst_stride[0],
                                                 dst_buf1[0], dst_stride1[0],
                                                 &xd->plane[0],
                                                 mi_row, mi_col,
                                                 mi_row_top, mi_col_top,
                                                 bsize, top_bsize,
                                                 PARTITION_VERT, 0);
      } else {
        // bsize: not important, not useful
        predict_b_extend(cpi, tile, 0, mi_row, mi_col, mi_row, mi_col,
                         mi_row_top, mi_col_top, dst_buf, dst_stride,
                         subsize, top_bsize, subsize, output_enabled, 0, 0);
        if (bsize < top_bsize)
          extend_all(cpi, tile, 0, subsize, top_bsize, mi_row, mi_col,
                     mi_row_top, mi_col_top, output_enabled,
                     dst_buf, dst_stride);
        else
          extend_dir(cpi, tile, 0, subsize, top_bsize, mi_row, mi_col,
                     mi_row_top, mi_col_top, output_enabled,
                     dst_buf, dst_stride, 3);


        if (mi_col + hbs < cm->mi_cols) {
          predict_b_extend(cpi, tile, 0, mi_row, mi_col + hbs,
                           mi_row, mi_col + hbs, mi_row_top, mi_col_top,
                           dst_buf1, dst_stride1, subsize, top_bsize, subsize,
                           output_enabled, 0, 0);
          if (bsize < top_bsize)
            extend_all(cpi, tile, 0, subsize, top_bsize, mi_row, mi_col + hbs,
                       mi_row_top, mi_col_top, output_enabled,
                       dst_buf1, dst_stride1);
          else
            extend_dir(cpi, tile, 0, subsize, top_bsize, mi_row, mi_col + hbs,
                       mi_row_top, mi_col_top, output_enabled,
                       dst_buf1, dst_stride1, 2);

          for (i = 0; i < MAX_MB_PLANE; i++) {
            xd->plane[i].dst.buf = dst_buf[i];
            xd->plane[i].dst.stride = dst_stride[i];
            vp9_build_masked_inter_predictor_complex(
                xd, dst_buf[i], dst_stride[i], dst_buf1[i], dst_stride1[i],
                &xd->plane[i], mi_row, mi_col, mi_row_top, mi_col_top,
                bsize, top_bsize, PARTITION_VERT, i);
          }
        }
      }
      break;
    case PARTITION_SPLIT:
      if (bsize == BLOCK_8X8) {
        predict_b_extend(cpi, tile, 0, mi_row, mi_col, mi_row, mi_col,
                         mi_row_top, mi_col_top, dst_buf, dst_stride,
                         subsize, top_bsize, BLOCK_8X8, output_enabled, 1, 0);
        predict_b_extend(cpi, tile, 1, mi_row, mi_col, mi_row, mi_col,
                         mi_row_top, mi_col_top, dst_buf1, dst_stride1,
                         subsize, top_bsize, BLOCK_8X8, output_enabled, 1, 1);
        predict_b_extend(cpi, tile, 2, mi_row, mi_col, mi_row, mi_col,
                         mi_row_top, mi_col_top, dst_buf2, dst_stride2,
                         subsize, top_bsize, BLOCK_8X8, output_enabled, 1, 1);
        predict_b_extend(cpi, tile, 3, mi_row, mi_col, mi_row, mi_col,
                         mi_row_top, mi_col_top, dst_buf3, dst_stride3,
                         subsize, top_bsize, BLOCK_8X8, output_enabled, 1, 1);

        if (bsize < top_bsize) {
          extend_all(cpi, tile, 0, subsize, top_bsize, mi_row, mi_col,
                     mi_row_top, mi_col_top, output_enabled,
                     dst_buf, dst_stride);
          extend_all(cpi, tile, 1, subsize, top_bsize, mi_row, mi_col,
                     mi_row_top, mi_col_top, output_enabled,
                     dst_buf1, dst_stride1);
          extend_all(cpi, tile, 2, subsize, top_bsize, mi_row, mi_col,
                     mi_row_top, mi_col_top, output_enabled,
                     dst_buf2, dst_stride2);
          extend_all(cpi, tile, 3, subsize, top_bsize, mi_row, mi_col,
                     mi_row_top, mi_col_top, output_enabled,
                     dst_buf3, dst_stride3);
        }
      } else {
        predict_sb_complex(cpi, tile, mi_row, mi_col,
                           mi_row_top, mi_col_top, output_enabled, subsize,
                           top_bsize, dst_buf, dst_stride,
                           pc_tree->split[0]);
        if (mi_row < cm->mi_rows && mi_col + hbs < cm->mi_cols)
          predict_sb_complex(cpi, tile, mi_row, mi_col + hbs,
                             mi_row_top, mi_col_top, output_enabled, subsize,
                             top_bsize, dst_buf1, dst_stride1,
                             pc_tree->split[1]);
        if (mi_row + hbs < cm->mi_rows && mi_col < cm->mi_cols)
          predict_sb_complex(cpi, tile, mi_row + hbs, mi_col,
                             mi_row_top, mi_col_top, output_enabled, subsize,
                             top_bsize, dst_buf2, dst_stride2,
                             pc_tree->split[2]);
        if (mi_row + hbs < cm->mi_rows && mi_col + hbs < cm->mi_cols)
          predict_sb_complex(cpi, tile, mi_row + hbs, mi_col + hbs,
                             mi_row_top, mi_col_top, output_enabled, subsize,
                             top_bsize, dst_buf3, dst_stride3,
                             pc_tree->split[3]);
      }
        for (i = 0; i < MAX_MB_PLANE; i++) {
          if (bsize == BLOCK_8X8 && i != 0)
            continue;  // Skip <4x4 chroma smoothing
          if (mi_row < cm->mi_rows && mi_col + hbs < cm->mi_cols) {
            vp9_build_masked_inter_predictor_complex(xd,
                                                     dst_buf[i],
                                                     dst_stride[i],
                                                     dst_buf1[i],
                                                     dst_stride1[i],
                                                     &xd->plane[i],
                                                     mi_row, mi_col,
                                                     mi_row_top, mi_col_top,
                                                     bsize, top_bsize,
                                                     PARTITION_VERT, i);
            if (mi_row + hbs < cm->mi_rows) {
              vp9_build_masked_inter_predictor_complex(xd,
                                                       dst_buf2[i],
                                                       dst_stride2[i],
                                                       dst_buf3[i],
                                                       dst_stride3[i],
                                                       &xd->plane[i],
                                                       mi_row, mi_col,
                                                       mi_row_top, mi_col_top,
                                                       bsize, top_bsize,
                                                       PARTITION_VERT, i);
              vp9_build_masked_inter_predictor_complex(xd,
                                                       dst_buf[i],
                                                       dst_stride[i],
                                                       dst_buf2[i],
                                                       dst_stride2[i],
                                                       &xd->plane[i],
                                                       mi_row, mi_col,
                                                       mi_row_top, mi_col_top,
                                                       bsize, top_bsize,
                                                       PARTITION_HORZ, i);
            }
          } else if (mi_row + hbs < cm->mi_rows && mi_col < cm->mi_cols) {
            vp9_build_masked_inter_predictor_complex(xd,
                                                     dst_buf[i],
                                                     dst_stride[i],
                                                     dst_buf2[i],
                                                     dst_stride2[i],
                                                     &xd->plane[i],
                                                     mi_row, mi_col,
                                                     mi_row_top, mi_col_top,
                                                     bsize, top_bsize,
                                                     PARTITION_HORZ, i);
          }
      }
      break;
#if CONFIG_EXT_PARTITION
    case PARTITION_HORZ_A:
      predict_b_extend(cpi, tile, 0, mi_row, mi_col, mi_row, mi_col,
                       mi_row_top, mi_col_top, dst_buf, dst_stride,
                       bsize2, top_bsize, bsize2, output_enabled, 0, 0);
      extend_all(cpi, tile, 0, bsize2, top_bsize, mi_row, mi_col,
                 mi_row_top, mi_col_top, output_enabled, dst_buf, dst_stride);

      predict_b_extend(cpi, tile, 0, mi_row, mi_col + hbs,
                       mi_row, mi_col + hbs, mi_row_top, mi_col_top,
                       dst_buf1, dst_stride1, bsize2, top_bsize, bsize2,
                       output_enabled, 0, 0);
      extend_all(cpi, tile, 0, bsize2, top_bsize, mi_row, mi_col + hbs,
                 mi_row_top, mi_col_top, output_enabled, dst_buf1, dst_stride1);

      predict_b_extend(cpi, tile, 0, mi_row + hbs, mi_col, mi_row + hbs, mi_col,
                       mi_row_top, mi_col_top, dst_buf2, dst_stride2,
                       subsize, top_bsize, subsize, output_enabled, 0, 0);
      if (bsize < top_bsize)
        extend_all(cpi, tile, 0, subsize, top_bsize, mi_row + hbs, mi_col,
                   mi_row_top, mi_col_top, output_enabled,
                   dst_buf2, dst_stride2);
      else
        extend_dir(cpi, tile, 0, subsize, top_bsize, mi_row + hbs, mi_col,
                   mi_row_top, mi_col_top, output_enabled,
                   dst_buf2, dst_stride2, 1);

      for (i = 0; i < MAX_MB_PLANE; i++) {
        xd->plane[i].dst.buf = dst_buf[i];
        xd->plane[i].dst.stride = dst_stride[i];
        vp9_build_masked_inter_predictor_complex(xd,
                                                 dst_buf[i], dst_stride[i],
                                                 dst_buf1[i], dst_stride1[i],
                                                 &xd->plane[i],
                                                 mi_row, mi_col,
                                                 mi_row_top, mi_col_top,
                                                 bsize, top_bsize,
                                                 PARTITION_VERT, i);
      }
      for (i = 0; i < MAX_MB_PLANE; i++) {
        vp9_build_masked_inter_predictor_complex(xd,
                                                 dst_buf[i], dst_stride[i],
                                                 dst_buf2[i], dst_stride2[i],
                                                 &xd->plane[i],
                                                 mi_row, mi_col,
                                                 mi_row_top, mi_col_top,
                                                 bsize, top_bsize,
                                                 PARTITION_HORZ, i);
      }

      break;
    case PARTITION_VERT_A:

      predict_b_extend(cpi, tile, 0, mi_row, mi_col, mi_row, mi_col,
                       mi_row_top, mi_col_top, dst_buf, dst_stride,
                       bsize2, top_bsize, bsize2, output_enabled, 0, 0);
      extend_all(cpi, tile, 0, bsize2, top_bsize, mi_row, mi_col,
                 mi_row_top, mi_col_top, output_enabled, dst_buf, dst_stride);

      predict_b_extend(cpi, tile, 0, mi_row + hbs, mi_col, mi_row + hbs, mi_col,
                       mi_row_top, mi_col_top, dst_buf1, dst_stride1,
                       bsize2, top_bsize, bsize2, output_enabled, 0, 0);
      extend_all(cpi, tile, 0, bsize2, top_bsize, mi_row + hbs, mi_col,
                 mi_row_top, mi_col_top, output_enabled, dst_buf1, dst_stride1);

      predict_b_extend(cpi, tile, 0, mi_row, mi_col + hbs, mi_row, mi_col + hbs,
                       mi_row_top, mi_col_top, dst_buf2, dst_stride2,
                       subsize, top_bsize, subsize, output_enabled, 0, 0);
      if (bsize < top_bsize)
        extend_all(cpi, tile, 0, subsize, top_bsize, mi_row, mi_col + hbs,
                   mi_row_top, mi_col_top, output_enabled,
                   dst_buf2, dst_stride2);
      else
        extend_dir(cpi, tile, 0, subsize, top_bsize, mi_row, mi_col + hbs,
                   mi_row_top, mi_col_top, output_enabled,
                   dst_buf2, dst_stride2, 2);

      for (i = 0; i < MAX_MB_PLANE; i++) {
        xd->plane[i].dst.buf = dst_buf[i];
        xd->plane[i].dst.stride = dst_stride[i];
        vp9_build_masked_inter_predictor_complex(xd,
                                                 dst_buf[i], dst_stride[i],
                                                 dst_buf1[i], dst_stride1[i],
                                                 &xd->plane[i],
                                                 mi_row, mi_col,
                                                 mi_row_top, mi_col_top,
                                                 bsize, top_bsize,
                                                 PARTITION_HORZ, i);
      }
      for (i = 0; i < MAX_MB_PLANE; i++) {
        vp9_build_masked_inter_predictor_complex(xd,
                                                 dst_buf[i], dst_stride[i],
                                                 dst_buf2[i], dst_stride2[i],
                                                 &xd->plane[i],
                                                 mi_row, mi_col,
                                                 mi_row_top, mi_col_top,
                                                 bsize, top_bsize,
                                                 PARTITION_VERT, i);
      }
      break;
    case PARTITION_HORZ_B:

      predict_b_extend(cpi, tile, 0, mi_row, mi_col, mi_row, mi_col,
                       mi_row_top, mi_col_top, dst_buf, dst_stride,
                       subsize, top_bsize, subsize, output_enabled, 0, 0);
      if (bsize < top_bsize)
        extend_all(cpi, tile, 0, subsize, top_bsize, mi_row, mi_col,
                   mi_row_top, mi_col_top, output_enabled, dst_buf, dst_stride);
      else
        extend_dir(cpi, tile, 0, subsize, top_bsize, mi_row, mi_col,
                   mi_row_top, mi_col_top, output_enabled,
                   dst_buf, dst_stride, 0);

      predict_b_extend(cpi, tile, 0, mi_row + hbs, mi_col, mi_row + hbs, mi_col,
                       mi_row_top, mi_col_top, dst_buf1, dst_stride1,
                       bsize2, top_bsize, bsize2, output_enabled, 0, 0);
      extend_all(cpi, tile, 0, bsize2, top_bsize, mi_row + hbs, mi_col,
                 mi_row_top, mi_col_top, output_enabled, dst_buf1, dst_stride1);

      predict_b_extend(cpi, tile, 0, mi_row + hbs, mi_col + hbs,
                       mi_row + hbs, mi_col + hbs, mi_row_top, mi_col_top,
                       dst_buf2, dst_stride2, bsize2, top_bsize, bsize2,
                       output_enabled, 0, 0);
      extend_all(cpi, tile, 0, bsize2, top_bsize, mi_row + hbs, mi_col + hbs,
                 mi_row_top, mi_col_top, output_enabled, dst_buf2, dst_stride2);

      for (i = 0; i < MAX_MB_PLANE; i++) {
        xd->plane[i].dst.buf = dst_buf1[i];
        xd->plane[i].dst.stride = dst_stride1[i];
        vp9_build_masked_inter_predictor_complex(xd,
                                                 dst_buf1[i], dst_stride1[i],
                                                 dst_buf2[i], dst_stride2[i],
                                                 &xd->plane[i],
                                                 mi_row, mi_col,
                                                 mi_row_top, mi_col_top,
                                                 bsize, top_bsize,
                                                 PARTITION_VERT, i);
      }
      for (i = 0; i < MAX_MB_PLANE; i++) {
        xd->plane[i].dst.buf = dst_buf[i];
        xd->plane[i].dst.stride = dst_stride[i];
        vp9_build_masked_inter_predictor_complex(xd,
                                                 dst_buf[i], dst_stride[i],
                                                 dst_buf1[i], dst_stride1[i],
                                                 &xd->plane[i],
                                                 mi_row, mi_col,
                                                 mi_row_top, mi_col_top,
                                                 bsize, top_bsize,
                                                 PARTITION_HORZ, i);
      }
      break;
    case PARTITION_VERT_B:

      predict_b_extend(cpi, tile, 0, mi_row, mi_col, mi_row, mi_col,
                       mi_row_top, mi_col_top, dst_buf, dst_stride,
                       subsize, top_bsize, subsize, output_enabled, 0, 0);
      if (bsize < top_bsize)
        extend_all(cpi, tile, 0, subsize, top_bsize, mi_row, mi_col,
                   mi_row_top, mi_col_top, output_enabled, dst_buf, dst_stride);
      else
        extend_dir(cpi, tile, 0, subsize, top_bsize, mi_row, mi_col,
                   mi_row_top, mi_col_top, output_enabled,
                   dst_buf, dst_stride, 3);

      predict_b_extend(cpi, tile, 0, mi_row, mi_col + hbs, mi_row, mi_col + hbs,
                       mi_row_top, mi_col_top, dst_buf1, dst_stride1,
                       bsize2, top_bsize, bsize2, output_enabled, 0, 0);
      extend_all(cpi, tile, 0, bsize2, top_bsize, mi_row, mi_col + hbs,
                 mi_row_top, mi_col_top, output_enabled, dst_buf1, dst_stride1);

      predict_b_extend(cpi, tile, 0, mi_row + hbs, mi_col + hbs,
                       mi_row + hbs, mi_col + hbs, mi_row_top, mi_col_top,
                       dst_buf2, dst_stride2, bsize2, top_bsize, bsize2,
                       output_enabled, 0, 0);
      extend_all(cpi, tile, 0, bsize2, top_bsize, mi_row + hbs, mi_col + hbs,
                 mi_row_top, mi_col_top, output_enabled, dst_buf2, dst_stride2);

      for (i = 0; i < MAX_MB_PLANE; i++) {
        xd->plane[i].dst.buf = dst_buf1[i];
        xd->plane[i].dst.stride = dst_stride1[i];
        vp9_build_masked_inter_predictor_complex(xd,
                                                 dst_buf1[i], dst_stride1[i],
                                                 dst_buf2[i], dst_stride2[i],
                                                 &xd->plane[i],
                                                 mi_row, mi_col,
                                                 mi_row_top, mi_col_top,
                                                 bsize, top_bsize,
                                                 PARTITION_HORZ, i);
      }
      for (i = 0; i < MAX_MB_PLANE; i++) {
        xd->plane[i].dst.buf = dst_buf[i];
        xd->plane[i].dst.stride = dst_stride[i];
        vp9_build_masked_inter_predictor_complex(xd,
                                                 dst_buf[i], dst_stride[i],
                                                 dst_buf1[i], dst_stride1[i],
                                                 &xd->plane[i],
                                                 mi_row, mi_col,
                                                 mi_row_top, mi_col_top,
                                                 bsize, top_bsize,
                                                 PARTITION_VERT, i);
      }
      break;
#endif
    default:
      assert(0);
  }


#if CONFIG_EXT_PARTITION
  if (bsize < top_bsize)
    update_ext_partition_context(xd, mi_row, mi_col, subsize, bsize, partition);
#else
  if (bsize < top_bsize && (partition != PARTITION_SPLIT || bsize == BLOCK_8X8))
    update_partition_context(xd, mi_row, mi_col, subsize, bsize);
#endif
}

static void rd_supertx_sb(VP9_COMP *cpi, const TileInfo *const tile,
                          int mi_row, int mi_col, BLOCK_SIZE bsize,
                          int *tmp_rate, int64_t *tmp_dist,
#if CONFIG_EXT_TX
                          EXT_TX_TYPE *best_tx,
#endif
#if CONFIG_NEW_QUANT && QUANT_PROFILES > 1
                          int *dq_index,
#endif  // CONFIG_NEW_QUANT && QUANT_PROFILES > 1
                          PC_TREE *pc_tree) {
  VP9_COMMON *const cm = &cpi->common;
  MACROBLOCK *const x = &cpi->mb;
  MACROBLOCKD *const xd = &x->e_mbd;
  int plane, pnskip, skippable, skippable_uv, rate_uv, this_rate,
      base_rate = *tmp_rate;
  int64_t sse, pnsse, sse_uv, this_dist, dist_uv;
  uint8_t *dst_buf[3];
  int dst_stride[3];
  TX_SIZE tx_size;
#if CONFIG_EXT_TX
  EXT_TX_TYPE txfm, best_tx_nostx = xd->mi[0].mbmi.ext_txfrm;
  int tmp_rate_tx = 0, skip_tx = 0;
  int64_t tmp_dist_tx = 0, rd_tx, bestrd_tx = INT64_MAX;
  uint8_t tmp_zcoeff_blk = 0;
#endif

  update_state_sb_supertx(cpi, tile, mi_row, mi_col, bsize, 0, pc_tree);
#if CONFIG_NEW_QUANT && QUANT_PROFILES > 1
  xd->mi[0].mbmi.dq_off_index = 0;
#endif  // CONFIG_NEW_QUANT && QUANT_PROFILES > 1
  vp9_setup_dst_planes(xd->plane, get_frame_new_buffer(cm),
                       mi_row, mi_col);
  for (plane = 0; plane < MAX_MB_PLANE; plane++) {
    dst_buf[plane] = xd->plane[plane].dst.buf;
    dst_stride[plane] = xd->plane[plane].dst.stride;
  }
  predict_sb_complex(cpi, tile, mi_row, mi_col, mi_row, mi_col,
                     0, bsize, bsize, dst_buf, dst_stride, pc_tree);

  set_offsets(cpi, tile, mi_row, mi_col, bsize);
#if CONFIG_EXT_TX
  *best_tx = NORM;
#endif

#if CONFIG_TX_SKIP
  xd->mi[0].mbmi.tx_skip[0] = 0;
  xd->mi[0].mbmi.tx_skip[1] = 0;
#endif

  // chroma
  skippable_uv = 1;
  rate_uv = 0;
  dist_uv = 0;
  sse_uv = 0;
  for (plane = 1; plane < MAX_MB_PLANE; ++plane) {
    tx_size = bsize_to_tx_size(bsize);
    tx_size = get_uv_tx_size_impl(tx_size, bsize,
                                  cm->subsampling_x, cm->subsampling_y);
    vp9_subtract_plane(x, bsize, plane);
    txfm_rd_in_plane_supertx(x, &this_rate, &this_dist, &pnskip, &pnsse,
                             INT64_MAX, plane, bsize, tx_size, 0);
    rate_uv += this_rate;
    dist_uv += this_dist;
    sse_uv += pnsse;
    skippable_uv &= pnskip;
  }

  // luma
  tx_size = bsize_to_tx_size(bsize);
  vp9_subtract_plane(x, bsize, 0);
#if CONFIG_EXT_TX
  for (txfm = NORM; txfm < GET_EXT_TX_TYPES(tx_size); txfm++) {
    xd->mi[0].mbmi.ext_txfrm = txfm;
#endif  // CONFIG_EXT_TX
    txfm_rd_in_plane_supertx(x, &this_rate, &this_dist, &pnskip, &pnsse,
                             INT64_MAX, 0, bsize, tx_size, 0);
    *tmp_rate = rate_uv + this_rate;
    *tmp_dist = dist_uv + this_dist;
    sse = sse_uv + pnsse;
    skippable = skippable_uv && pnskip;

    if (skippable) {
      *tmp_rate = vp9_cost_bit(vp9_get_skip_prob(cm, xd), 1);
      x->skip = 1;
    } else {
#if CONFIG_EXT_TX
#if !CONFIG_WAVELETS
      if (tx_size <= TX_16X16)
#endif
        *tmp_rate += cpi->ext_tx_costs[tx_size][txfm];
#endif  // CONFIG_EXT_TX
      if (RDCOST(x->rdmult, x->rddiv, *tmp_rate, *tmp_dist)
          < RDCOST(x->rdmult, x->rddiv, 0, sse)) {
        *tmp_rate += vp9_cost_bit(vp9_get_skip_prob(cm, xd), 0);
        x->skip = 0;
      } else {
        *tmp_dist = sse;
        *tmp_rate = vp9_cost_bit(vp9_get_skip_prob(cm, xd), 1);
        x->skip = 1;
      }
    }
    *tmp_rate += base_rate;
#if CONFIG_EXT_TX
    rd_tx = RDCOST(x->rdmult, x->rddiv, *tmp_rate, *tmp_dist);
    if (rd_tx < bestrd_tx * 0.99 || txfm == NORM) {
      *best_tx = txfm;
      bestrd_tx = rd_tx;
      tmp_rate_tx = *tmp_rate;
      tmp_dist_tx = *tmp_dist;
      skip_tx = x->skip;
      tmp_zcoeff_blk = x->zcoeff_blk[tx_size][0];
    }
  }
  x->zcoeff_blk[tx_size][0] = tmp_zcoeff_blk;
  *tmp_rate = tmp_rate_tx;
  *tmp_dist = tmp_dist_tx;
  x->skip = skip_tx;
  xd->mi[0].mbmi.ext_txfrm = best_tx_nostx;
#endif  // CONFIG_EXT_TX
#if CONFIG_NEW_QUANT && QUANT_PROFILES > 1
  *dq_index = 0;
#endif  // CONFIG_NEW_QUANT && QUANT_PROFILES > 1 &&
}
#endif  // CONFIG_SUPERTX
