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
#include "vp9/common/vp9_loopfilter.h"
#include "vp9/common/vp9_onyxc_int.h"
#include "vp9/common/vp9_reconinter.h"
#include "vpx_mem/vpx_mem.h"

#include "vp9/common/vp9_seg_common.h"

struct loop_filter_info {
  const uint8_t *mblim;
  const uint8_t *lim;
  const uint8_t *hev_thr;
};

// This structure holds bit masks for all 8x8 blocks in a 64x64 region.
// Each 1 bit represents a position in which we want to apply the loop filter.
// Left_ entries refer to whether we apply a filter on the border to the
// left of the block.   Above_ entries refer to whether or not to apply a
// filter on the above border.   Int_ entries refer to whether or not to
// apply borders on the 4x4 edges within the 8x8 block that each bit
// represents.
// Since each transform is accompanied by a potentially different type of
// loop filter there is a different entry in the array for each transform size.
typedef struct {
  uint64_t left_y[TX_SIZES];
  uint64_t above_y[TX_SIZES];
  uint64_t int_4x4_y;
  uint16_t left_uv[TX_SIZES];
  uint16_t above_uv[TX_SIZES];
  uint16_t int_4x4_uv;
} LOOP_FILTER_MASK;

// 64 bit masks for left transform size.  Each 1 represents a position where
// we should apply a loop filter across the left border of an 8x8 block
// boundary.
//
// In the case of TX_16X16->  ( in low order byte first we end up with
// a mask that looks like this
//
//    10101010
//    10101010
//    10101010
//    10101010
//    10101010
//    10101010
//    10101010
//    10101010
//
// A loopfilter should be applied to every other 8x8 horizontally.
static const uint64_t left_64x64_txform_mask[TX_SIZES]= {
    0xffffffffffffffff,  // TX_4X4
    0xffffffffffffffff,  // TX_8x8
    0x5555555555555555,  // TX_16x16
    0x1111111111111111,  // TX_32x32
};

// 64 bit masks for above transform size.  Each 1 represents a position where
// we should apply a loop filter across the top border of an 8x8 block
// boundary.
//
// In the case of TX_32x32 ->  ( in low order byte first we end up with
// a mask that looks like this
//
//    11111111
//    00000000
//    00000000
//    00000000
//    11111111
//    00000000
//    00000000
//    00000000
//
// A loopfilter should be applied to every other 4 the row vertically.
static const uint64_t above_64x64_txform_mask[TX_SIZES]= {
    0xffffffffffffffff,  // TX_4X4
    0xffffffffffffffff,  // TX_8x8
    0x00ff00ff00ff00ff,  // TX_16x16
    0x000000ff000000ff,  // TX_32x32
};

// 64 bit masks for prediction sizes (left).  Each 1 represents a position
// where left border of an 8x8 block.  These are aligned to the right most
// appropriate bit,  and then shifted into place.
//
// In the case of TX_16x32 ->  ( low order byte first ) we end up with
// a mask that looks like this :
//
//  10000000
//  10000000
//  10000000
//  10000000
//  00000000
//  00000000
//  00000000
//  00000000
static const uint64_t left_prediction_mask[BLOCK_SIZES] = {
    0x0000000000000001,  // BLOCK_4X4,
    0x0000000000000001,  // BLOCK_4X8,
    0x0000000000000001,  // BLOCK_8X4,
    0x0000000000000001,  // BLOCK_8X8,
    0x0000000000000101,  // BLOCK_8X16,
    0x0000000000000001,  // BLOCK_16X8,
    0x0000000000000101,  // BLOCK_16X16,
    0x0000000001010101,  // BLOCK_16X32,
    0x0000000000000101,  // BLOCK_32X16,
    0x0000000001010101,  // BLOCK_32X32,
    0x0101010101010101,  // BLOCK_32X64,
    0x0000000001010101,  // BLOCK_64X32,
    0x0101010101010101,  // BLOCK_64X64
};

// 64 bit mask to shift and set for each prediction size.
static const uint64_t above_prediction_mask[BLOCK_SIZES] = {
    0x0000000000000001,  // BLOCK_4X4
    0x0000000000000001,  // BLOCK_4X8
    0x0000000000000001,  // BLOCK_8X4
    0x0000000000000001,  // BLOCK_8X8
    0x0000000000000001,  // BLOCK_8X16,
    0x0000000000000003,  // BLOCK_16X8
    0x0000000000000003,  // BLOCK_16X16
    0x0000000000000003,  // BLOCK_16X32,
    0x000000000000000f,  // BLOCK_32X16,
    0x000000000000000f,  // BLOCK_32X32,
    0x000000000000000f,  // BLOCK_32X64,
    0x00000000000000ff,  // BLOCK_64X32,
    0x00000000000000ff,  // BLOCK_64X64
};
// 64 bit mask to shift and set for each prediction size.  A bit is set for
// each 8x8 block that would be in the left most block of the given block
// size in the 64x64 block.
static const uint64_t size_mask[BLOCK_SIZES] = {
    0x0000000000000001,  // BLOCK_4X4
    0x0000000000000001,  // BLOCK_4X8
    0x0000000000000001,  // BLOCK_8X4
    0x0000000000000001,  // BLOCK_8X8
    0x0000000000000101,  // BLOCK_8X16,
    0x0000000000000003,  // BLOCK_16X8
    0x0000000000000303,  // BLOCK_16X16
    0x0000000003030303,  // BLOCK_16X32,
    0x0000000000000f0f,  // BLOCK_32X16,
    0x000000000f0f0f0f,  // BLOCK_32X32,
    0x0f0f0f0f0f0f0f0f,  // BLOCK_32X64,
    0x00000000ffffffff,  // BLOCK_64X32,
    0xffffffffffffffff,  // BLOCK_64X64
};

// These are used for masking the left and above borders.
static const uint64_t left_border =  0x1111111111111111;
static const uint64_t above_border = 0x000000ff000000ff;

// 16 bit masks for uv transform sizes.
static const uint16_t left_64x64_txform_mask_uv[TX_SIZES]= {
    0xffff,  // TX_4X4
    0xffff,  // TX_8x8
    0x5555,  // TX_16x16
    0x1111,  // TX_32x32
};

static const uint16_t above_64x64_txform_mask_uv[TX_SIZES]= {
    0xffff,  // TX_4X4
    0xffff,  // TX_8x8
    0x0f0f,  // TX_16x16
    0x000f,  // TX_32x32
};

// 16 bit left mask to shift and set for each uv prediction size.
static const uint16_t left_prediction_mask_uv[BLOCK_SIZES] = {
    0x0001,  // BLOCK_4X4,
    0x0001,  // BLOCK_4X8,
    0x0001,  // BLOCK_8X4,
    0x0001,  // BLOCK_8X8,
    0x0001,  // BLOCK_8X16,
    0x0001,  // BLOCK_16X8,
    0x0001,  // BLOCK_16X16,
    0x0011,  // BLOCK_16X32,
    0x0001,  // BLOCK_32X16,
    0x0011,  // BLOCK_32X32,
    0x1111,  // BLOCK_32X64
    0x0011,  // BLOCK_64X32,
    0x1111,  // BLOCK_64X64
};
// 16 bit above mask to shift and set for uv each prediction size.
static const uint16_t above_prediction_mask_uv[BLOCK_SIZES] = {
    0x0001,  // BLOCK_4X4
    0x0001,  // BLOCK_4X8
    0x0001,  // BLOCK_8X4
    0x0001,  // BLOCK_8X8
    0x0001,  // BLOCK_8X16,
    0x0001,  // BLOCK_16X8
    0x0001,  // BLOCK_16X16
    0x0001,  // BLOCK_16X32,
    0x0003,  // BLOCK_32X16,
    0x0003,  // BLOCK_32X32,
    0x0003,  // BLOCK_32X64,
    0x000f,  // BLOCK_64X32,
    0x000f,  // BLOCK_64X64
};

// 64 bit mask to shift and set for each uv prediction size
static const uint16_t size_mask_uv[BLOCK_SIZES] = {
    0x0001,  // BLOCK_4X4
    0x0001,  // BLOCK_4X8
    0x0001,  // BLOCK_8X4
    0x0001,  // BLOCK_8X8
    0x0001,  // BLOCK_8X16,
    0x0001,  // BLOCK_16X8
    0x0001,  // BLOCK_16X16
    0x0011,  // BLOCK_16X32,
    0x0003,  // BLOCK_32X16,
    0x0033,  // BLOCK_32X32,
    0x3333,  // BLOCK_32X64,
    0x00ff,  // BLOCK_64X32,
    0xffff,  // BLOCK_64X64
};
static const uint16_t left_border_uv =  0x1111;
static const uint16_t above_border_uv = 0x000f;


static void lf_init_lut(loop_filter_info_n *lfi) {
  lfi->mode_lf_lut[DC_PRED] = 0;
  lfi->mode_lf_lut[D45_PRED] = 0;
  lfi->mode_lf_lut[D135_PRED] = 0;
  lfi->mode_lf_lut[D117_PRED] = 0;
  lfi->mode_lf_lut[D153_PRED] = 0;
  lfi->mode_lf_lut[D207_PRED] = 0;
  lfi->mode_lf_lut[D63_PRED] = 0;
  lfi->mode_lf_lut[V_PRED] = 0;
  lfi->mode_lf_lut[H_PRED] = 0;
  lfi->mode_lf_lut[TM_PRED] = 0;
  lfi->mode_lf_lut[ZEROMV]  = 0;
  lfi->mode_lf_lut[NEARESTMV] = 1;
  lfi->mode_lf_lut[NEARMV] = 1;
  lfi->mode_lf_lut[NEWMV] = 1;
}

static void update_sharpness(loop_filter_info_n *lfi, int sharpness_lvl) {
  int lvl;

  // For each possible value for the loop filter fill out limits
  for (lvl = 0; lvl <= MAX_LOOP_FILTER; lvl++) {
    // Set loop filter paramaeters that control sharpness.
    int block_inside_limit = lvl >> ((sharpness_lvl > 0) + (sharpness_lvl > 4));

    if (sharpness_lvl > 0) {
      if (block_inside_limit > (9 - sharpness_lvl))
        block_inside_limit = (9 - sharpness_lvl);
    }

    if (block_inside_limit < 1)
      block_inside_limit = 1;

    vpx_memset(lfi->lim[lvl], block_inside_limit, SIMD_WIDTH);
    vpx_memset(lfi->mblim[lvl], (2 * (lvl + 2) + block_inside_limit),
               SIMD_WIDTH);
  }
}

void vp9_loop_filter_init(VP9_COMMON *cm) {
  loop_filter_info_n *lfi = &cm->lf_info;
  struct loopfilter *lf = &cm->lf;
  int i;

  // init limits for given sharpness
  update_sharpness(lfi, lf->sharpness_level);
  lf->last_sharpness_level = lf->sharpness_level;

  // init LUT for lvl  and hev thr picking
  lf_init_lut(lfi);

  // init hev threshold const vectors
  for (i = 0; i < 4; i++)
    vpx_memset(lfi->hev_thr[i], i, SIMD_WIDTH);
}

void vp9_loop_filter_frame_init(VP9_COMMON *cm, int default_filt_lvl) {
  int seg_id;
  // n_shift is the a multiplier for lf_deltas
  // the multiplier is 1 for when filter_lvl is between 0 and 31;
  // 2 when filter_lvl is between 32 and 63
  const int n_shift = default_filt_lvl >> 5;
  loop_filter_info_n *const lfi = &cm->lf_info;
  struct loopfilter *const lf = &cm->lf;
  struct segmentation *const seg = &cm->seg;

  // update limits if sharpness has changed
  if (lf->last_sharpness_level != lf->sharpness_level) {
    update_sharpness(lfi, lf->sharpness_level);
    lf->last_sharpness_level = lf->sharpness_level;
  }

  for (seg_id = 0; seg_id < MAX_SEGMENTS; seg_id++) {
    int lvl_seg = default_filt_lvl, ref, mode, intra_lvl;

    // Set the baseline filter values for each segment
    if (vp9_segfeature_active(seg, seg_id, SEG_LVL_ALT_LF)) {
      const int data = vp9_get_segdata(seg, seg_id, SEG_LVL_ALT_LF);
      lvl_seg = seg->abs_delta == SEGMENT_ABSDATA
                  ? data
                  : clamp(default_filt_lvl + data, 0, MAX_LOOP_FILTER);
    }

    if (!lf->mode_ref_delta_enabled) {
      // we could get rid of this if we assume that deltas are set to
      // zero when not in use; encoder always uses deltas
      vpx_memset(lfi->lvl[seg_id], lvl_seg, sizeof(lfi->lvl[seg_id]));
      continue;
    }

    intra_lvl = lvl_seg + (lf->ref_deltas[INTRA_FRAME] << n_shift);
    lfi->lvl[seg_id][INTRA_FRAME][0] = clamp(intra_lvl, 0, MAX_LOOP_FILTER);

    for (ref = LAST_FRAME; ref < MAX_REF_FRAMES; ++ref)
      for (mode = 0; mode < MAX_MODE_LF_DELTAS; ++mode) {
        const int inter_lvl = lvl_seg + (lf->ref_deltas[ref] << n_shift)
                                      + (lf->mode_deltas[mode] << n_shift);
        lfi->lvl[seg_id][ref][mode] = clamp(inter_lvl, 0, MAX_LOOP_FILTER);
      }
  }
}

static int build_lfi(const loop_filter_info_n *lfi_n,
                     const MB_MODE_INFO *mbmi,
                     struct loop_filter_info *lfi) {
  const int seg = mbmi->segment_id;
  const int ref = mbmi->ref_frame[0];
  const int mode = lfi_n->mode_lf_lut[mbmi->mode];
  const int filter_level = lfi_n->lvl[seg][ref][mode];

  if (filter_level > 0) {
    lfi->mblim = lfi_n->mblim[filter_level];
    lfi->lim = lfi_n->lim[filter_level];
    lfi->hev_thr = lfi_n->hev_thr[filter_level >> 4];
    return 1;
  } else {
    return 0;
  }
}

static void filter_selectively_vert(uint8_t *s, int pitch,
                                    unsigned int mask_16x16,
                                    unsigned int mask_8x8,
                                    unsigned int mask_4x4,
                                    unsigned int mask_4x4_int,
                                    const struct loop_filter_info *lfi) {
  unsigned int mask;

  for (mask = mask_16x16 | mask_8x8 | mask_4x4 | mask_4x4_int;
       mask; mask >>= 1) {
    if (mask & 1) {
      if (mask_16x16 & 1) {
        vp9_mb_lpf_vertical_edge_w(s, pitch, lfi->mblim, lfi->lim,
                                   lfi->hev_thr);
        assert(!(mask_8x8 & 1));
        assert(!(mask_4x4 & 1));
        assert(!(mask_4x4_int & 1));
      } else if (mask_8x8 & 1) {
        vp9_mbloop_filter_vertical_edge(s, pitch, lfi->mblim, lfi->lim,
                                        lfi->hev_thr, 1);
        assert(!(mask_16x16 & 1));
        assert(!(mask_4x4 & 1));
      } else if (mask_4x4 & 1) {
        vp9_loop_filter_vertical_edge(s, pitch, lfi->mblim, lfi->lim,
                                      lfi->hev_thr, 1);
        assert(!(mask_16x16 & 1));
        assert(!(mask_8x8 & 1));
      }
    }
    if (mask_4x4_int & 1)
      vp9_loop_filter_vertical_edge(s + 4, pitch, lfi->mblim, lfi->lim,
                                    lfi->hev_thr, 1);
    s += 8;
    lfi++;
    mask_16x16 >>= 1;
    mask_8x8 >>= 1;
    mask_4x4 >>= 1;
    mask_4x4_int >>= 1;
  }
}

static void filter_selectively_horiz(uint8_t *s, int pitch,
                                     unsigned int mask_16x16,
                                     unsigned int mask_8x8,
                                     unsigned int mask_4x4,
                                     unsigned int mask_4x4_int,
                                     int only_4x4_1,
                                     const struct loop_filter_info *lfi) {
  unsigned int mask;
  int count;

  for (mask = mask_16x16 | mask_8x8 | mask_4x4 | mask_4x4_int;
       mask; mask >>= count) {
    count = 1;
    if (mask & 1) {
      if (!only_4x4_1) {
        if (mask_16x16 & 1) {
          if ((mask_16x16 & 3) == 3) {
            vp9_mb_lpf_horizontal_edge_w(s, pitch, lfi->mblim, lfi->lim,
                                         lfi->hev_thr, 2);
            count = 2;
          } else {
            vp9_mb_lpf_horizontal_edge_w(s, pitch, lfi->mblim, lfi->lim,
                                         lfi->hev_thr, 1);
          }
          assert(!(mask_8x8 & 1));
          assert(!(mask_4x4 & 1));
          assert(!(mask_4x4_int & 1));
        } else if (mask_8x8 & 1) {
          vp9_mbloop_filter_horizontal_edge(s, pitch, lfi->mblim, lfi->lim,
                                            lfi->hev_thr, 1);
          assert(!(mask_16x16 & 1));
          assert(!(mask_4x4 & 1));
        } else if (mask_4x4 & 1) {
          vp9_loop_filter_horizontal_edge(s, pitch, lfi->mblim, lfi->lim,
                                          lfi->hev_thr, 1);
          assert(!(mask_16x16 & 1));
          assert(!(mask_8x8 & 1));
        }
      }

      if (mask_4x4_int & 1)
        vp9_loop_filter_horizontal_edge(s + 4 * pitch, pitch, lfi->mblim,
                                        lfi->lim, lfi->hev_thr, 1);
    }
    s += 8 * count;
    lfi += count;
    mask_16x16 >>= count;
    mask_8x8 >>= count;
    mask_4x4 >>= count;
    mask_4x4_int >>= count;
  }
}

// This function ors into the current lfm structure, where to do loop
// filters for the specific mi we are looking at.   It uses information
// including the block_size_type (32x16, 32x32, etc),  the transform size,
// whether there were any coefficients encoded, and the loop filter strength
// block we are currently looking at. Shift is used to position the
// 1's we produce.
// TODO(JBB) Need another function for different resolution color..
static void build_masks(const loop_filter_info_n *const lfi_n,
                        const MODE_INFO *mi, const int shift_y,
                        const int shift_uv,
                        LOOP_FILTER_MASK *lfm) {
  const BLOCK_SIZE block_size = mi->mbmi.sb_type;
  const TX_SIZE tx_size_y = mi->mbmi.tx_size;
  const TX_SIZE tx_size_uv = get_uv_tx_size(&mi->mbmi);
  const int skip = mi->mbmi.skip_coeff;
  const int seg = mi->mbmi.segment_id;
  const int ref = mi->mbmi.ref_frame[0];
  const int mode = lfi_n->mode_lf_lut[mi->mbmi.mode];
  const int filter_level = lfi_n->lvl[seg][ref][mode];
  uint64_t *left_y = &lfm->left_y[tx_size_y];
  uint64_t *above_y = &lfm->above_y[tx_size_y];
  uint64_t *int_4x4_y = &lfm->int_4x4_y;
  uint16_t *left_uv = &lfm->left_uv[tx_size_uv];
  uint16_t *above_uv = &lfm->above_uv[tx_size_uv];
  uint16_t *int_4x4_uv = &lfm->int_4x4_uv;

  // If filter level is 0 we don't loop filter.
  if (!filter_level)
    return;

  // These set 1 in the current block size for the block size edges.
  // For instance if the block size is 32x16,   we'll set :
  //    above =   1111
  //              0000
  //    and
  //    left  =   1000
  //          =   1000
  // NOTE : In this example the low bit is left most ( 1000 ) is stored as
  //        1,  not 8...
  //
  // U and v set things on a 16 bit scale.
  //
  *above_y |= above_prediction_mask[block_size] << shift_y;
  *above_uv |= above_prediction_mask_uv[block_size] << shift_uv;
  *left_y |= left_prediction_mask[block_size] << shift_y;
  *left_uv |= left_prediction_mask_uv[block_size] << shift_uv;

  // If the block has no coefficients and is not intra we skip applying
  // the loop filter on block edges.
  if (skip && ref > INTRA_FRAME)
    return;

  // Here we are adding a mask for the transform size.  The transform
  // size mask is set to be correct for a 64x64 prediction block size. We
  // mask to match the size of the block we are working on and then shift it
  // into place..
  *above_y |= (size_mask[block_size] &
               above_64x64_txform_mask[tx_size_y]) << shift_y;
  *above_uv |= (size_mask_uv[block_size] &
                above_64x64_txform_mask_uv[tx_size_uv]) << shift_uv;

  *left_y |= (size_mask[block_size] &
              left_64x64_txform_mask[tx_size_y]) << shift_y;
  *left_uv |= (size_mask_uv[block_size] &
               left_64x64_txform_mask_uv[tx_size_uv]) << shift_uv;

  // Here we are trying to determine what to do with the internal 4x4 block
  // boundaries.  These differ from the 4x4 boundaries on the outside edge of
  // an 8x8 in that the internal ones can be skipped and don't depend on
  // the prediction block size.
  if (tx_size_y == TX_4X4) {
    *int_4x4_y |= (size_mask[block_size] & 0xffffffffffffffff) << shift_y;
  }
  if (tx_size_uv == TX_4X4) {
    *int_4x4_uv |= (size_mask_uv[block_size] & 0xffff) << shift_uv;
  }
}

// This function does the same thing as the one above with the exception that
// it only affects the y masks.   It exists because for blocks < 16x16 in size,
// we only update u and v masks on the first block.
static void build_y_mask(const loop_filter_info_n *const lfi_n,
                         const MODE_INFO *mi, const int shift_y,
                         LOOP_FILTER_MASK *lfm) {
  const BLOCK_SIZE block_size = mi->mbmi.sb_type;
  const TX_SIZE tx_size_y = mi->mbmi.tx_size;
  const int skip = mi->mbmi.skip_coeff;
  const int seg = mi->mbmi.segment_id;
  const int ref = mi->mbmi.ref_frame[0];
  const int mode = lfi_n->mode_lf_lut[mi->mbmi.mode];
  const int filter_level = lfi_n->lvl[seg][ref][mode];
  uint64_t *left_y = &lfm->left_y[tx_size_y];
  uint64_t *above_y = &lfm->above_y[tx_size_y];
  uint64_t *int_4x4_y = &lfm->int_4x4_y;

  if (!filter_level)
    return;

  *above_y |= above_prediction_mask[block_size] << shift_y;
  *left_y |= left_prediction_mask[block_size] << shift_y;

  if (skip && ref > INTRA_FRAME)
    return;

  *above_y |= (size_mask[block_size] &
               above_64x64_txform_mask[tx_size_y]) << shift_y;

  *left_y |= (size_mask[block_size] &
              left_64x64_txform_mask[tx_size_y]) << shift_y;

  if (tx_size_y == TX_4X4) {
    *int_4x4_y |= (size_mask[block_size] & 0xffffffffffffffff) << shift_y;
  }
}

// This function sets up the bit masks for the entire 64x64 region represented
// by mi_row, mi_col.
// TODO(JBB): This function only works for yv12.
static void setup_mask(VP9_COMMON *const cm, const int mi_row, const int mi_col,
                       const MODE_INFO *mi, const int mode_info_stride,
                       LOOP_FILTER_MASK *lfm) {
  int idx_32, idx_16, idx_8;
  const loop_filter_info_n *const lfi_n = &cm->lf_info;
  const MODE_INFO *mip = mi;
  const MODE_INFO *mip2 = mi;

  // These are offsets to the next mi in the 64x64 block. It is what gets
  // added to the mi ptr as we go through each loop.  It helps us to avoids
  // setting up special row and column counters for each index.  The last step
  // brings us out back to the starting position.
  const int offset_32[] = {4, (mode_info_stride << 2) - 4, 4,
                           -(mode_info_stride << 2) - 4};
  const int offset_16[] = {2, (mode_info_stride << 1) - 2, 2,
                           -(mode_info_stride << 1) - 2};
  const int offset[] = {1, mode_info_stride - 1, 1, -mode_info_stride - 1};

  // Following variables represent shifts to position the current block
  // mask over the appropriate block.   A shift of 36 to the left will move
  // the bits for the final 32 by 32 block in the 64x64 up 4 rows and left
  // 4 rows to the appropriate spot.
  const int shift_32_y[] = {0, 4, 32, 36};
  const int shift_16_y[] = {0, 2, 16, 18};
  const int shift_8_y[] = {0, 1, 8, 9};
  const int shift_32_uv[] = {0, 2, 8, 10};
  const int shift_16_uv[] = {0, 1, 4, 5};
  int i;
  const int max_rows = (mi_row + MI_BLOCK_SIZE > cm->mi_rows ?
                        cm->mi_rows - mi_row : MI_BLOCK_SIZE);
  const int max_cols = (mi_col + MI_BLOCK_SIZE > cm->mi_cols ?
                        cm->mi_cols - mi_col : MI_BLOCK_SIZE);

  vp9_zero(*lfm);

  // TODO(jimbankoski): Try moving most of the following code into decode
  // loop and storing lfm in the mbmi structure so that we don't have to go
  // through the recursive loop structure multiple times.
  switch (mip->mbmi.sb_type) {
    case BLOCK_64X64:
      build_masks(lfi_n, mip , 0, 0, lfm);
      break;
    case BLOCK_64X32:
      build_masks(lfi_n, mip, 0, 0, lfm);
      mip2 = mip + mode_info_stride * 4;
      build_masks(lfi_n, mip2 , 32, 8, lfm);
      break;
    case BLOCK_32X64:
      build_masks(lfi_n, mip, 0, 0, lfm);
      mip2 = mip + 4;
      build_masks(lfi_n, mip2, 4, 2, lfm);
      break;
    default:
      for (idx_32 = 0; idx_32 < 4; mip += offset_32[idx_32], ++idx_32) {
        const int shift_y = shift_32_y[idx_32];
        const int shift_uv = shift_32_uv[idx_32];
        const int mi_32_col_offset = ((idx_32 & 1) << 2);
        const int mi_32_row_offset = ((idx_32 >> 1) << 2);
        if (mi_32_col_offset >= max_cols || mi_32_row_offset >= max_rows)
          continue;
        switch (mip->mbmi.sb_type) {
          case BLOCK_32X32:
            build_masks(lfi_n, mip, shift_y, shift_uv, lfm);
            break;
          case BLOCK_32X16:
            build_masks(lfi_n, mip, shift_y, shift_uv, lfm);
            mip2 = mip + mode_info_stride * 2;
            build_masks(lfi_n, mip2, shift_y + 16, shift_uv + 4, lfm);
            break;
          case BLOCK_16X32:
            build_masks(lfi_n, mip, shift_y, shift_uv, lfm);
            mip2 = mip + 2;
            build_masks(lfi_n, mip2, shift_y + 2, shift_uv + 1, lfm);
            break;
          default:
            for (idx_16 = 0; idx_16 < 4; mip += offset_16[idx_16], ++idx_16) {
              const int shift_y = shift_32_y[idx_32] + shift_16_y[idx_16];
              const int shift_uv = shift_32_uv[idx_32] + shift_16_uv[idx_16];
              const int mi_16_col_offset = mi_32_col_offset +
                  ((idx_16 & 1) << 1);
              const int mi_16_row_offset = mi_32_row_offset +
                  ((idx_16 >> 1) << 1);

              if (mi_16_col_offset >= max_cols || mi_16_row_offset >= max_rows)
                continue;

              switch (mip->mbmi.sb_type) {
                case BLOCK_16X16:
                  build_masks(lfi_n, mip, shift_y, shift_uv, lfm);
                  break;
                case BLOCK_16X8:
                  build_masks(lfi_n, mip, shift_y, shift_uv, lfm);
                  mip2 = mip + mode_info_stride;
                  build_y_mask(lfi_n, mip2, shift_y+8, lfm);
                  break;
                case BLOCK_8X16:
                  build_masks(lfi_n, mip, shift_y, shift_uv, lfm);
                  mip2 = mip + 1;
                  build_y_mask(lfi_n, mip2, shift_y+1, lfm);
                  break;
                default: {
                  const int shift_y = shift_32_y[idx_32] +
                                      shift_16_y[idx_16] +
                                      shift_8_y[0];
                  build_masks(lfi_n, mip, shift_y, shift_uv, lfm);
                  mip += offset[0];
                  for (idx_8 = 1; idx_8 < 4; mip += offset[idx_8], ++idx_8) {
                    const int shift_y = shift_32_y[idx_32] +
                                        shift_16_y[idx_16] +
                                        shift_8_y[idx_8];
                    const int mi_8_col_offset = mi_16_col_offset +
                        ((idx_8 & 1));
                    const int mi_8_row_offset = mi_16_row_offset +
                        ((idx_8 >> 1));

                    if (mi_8_col_offset >= max_cols ||
                        mi_8_row_offset >= max_rows)
                      continue;
                    build_y_mask(lfi_n, mip, shift_y, lfm);
                  }
                  break;
                }
              }
            }
            break;
        }
      }
      break;
  }
  // The largest loopfilter we have is 16x16 so we use the 16x16 mask
  // for 32x32 transforms also also.
  lfm->left_y[TX_16X16] |= lfm->left_y[TX_32X32];
  lfm->above_y[TX_16X16] |= lfm->above_y[TX_32X32];
  lfm->left_uv[TX_16X16] |= lfm->left_uv[TX_32X32];
  lfm->above_uv[TX_16X16] |= lfm->above_uv[TX_32X32];

  // We do at least 8 tap filter on every 32x32 even if the transform size
  // is 4x4.  So if the 4x4 is set on a border pixel add it to the 8x8 and
  // remove it from the 4x4.
  lfm->left_y[TX_8X8] |= lfm->left_y[TX_4X4] & left_border;
  lfm->left_y[TX_4X4] &= ~left_border;
  lfm->above_y[TX_8X8] |= lfm->above_y[TX_4X4] & above_border;
  lfm->above_y[TX_4X4] &= ~above_border;
  lfm->left_uv[TX_8X8] |= lfm->left_uv[TX_4X4] & left_border_uv;
  lfm->left_uv[TX_4X4] &= ~left_border_uv;
  lfm->above_uv[TX_8X8] |= lfm->above_uv[TX_4X4] & above_border_uv;
  lfm->above_uv[TX_4X4] &= ~above_border_uv;

  // We do some special edge handling.
  if (mi_row + MI_BLOCK_SIZE > cm->mi_rows) {
    const uint64_t rows = cm->mi_rows - mi_row;

    // Each pixel inside the border gets a 1,
    const uint64_t mask_y = (((uint64_t) 1 << (rows << 3)) - 1);
    const uint16_t mask_uv = (((uint16_t) 1 << (((rows + 1) >> 1) << 2)) - 1);

    // Remove values completely outside our border.
    for (i = 0; i < TX_32X32; i++) {
      lfm->left_y[i] &= mask_y;
      lfm->above_y[i] &= mask_y;
      lfm->left_uv[i] &= mask_uv;
      lfm->above_uv[i] &= mask_uv;
    }
    lfm->int_4x4_y &= mask_y;
    lfm->int_4x4_uv &= mask_uv;

    // We don't apply a wide loop filter on the last uv block row.  If set
    // apply the shorter one instead.
    if (rows == 1) {
      lfm->above_uv[TX_8X8] |= lfm->above_uv[TX_16X16];
      lfm->above_uv[TX_16X16] = 0;
    }
    if (rows == 5) {
      lfm->above_uv[TX_8X8] |= lfm->above_uv[TX_16X16] & 0xff00;
      lfm->above_uv[TX_16X16] &= ~(lfm->above_uv[TX_16X16] & 0xff00);
    }
  }

  if (mi_col + MI_BLOCK_SIZE > cm->mi_cols) {
    const uint64_t columns = cm->mi_cols - mi_col;

    // Each pixel inside the border gets a 1, the multiply copies the border
    // to where we need it.
    const uint64_t mask_y  = (((1 << columns) - 1)) * 0x0101010101010101;
    const uint16_t mask_uv = ((1 << ((columns + 1) >> 1)) - 1) * 0x1111;

    // Internal edges are not applied on the last column of the image so
    // we mask 1 more for the internal edges
    const uint16_t mask_uv_int = ((1 << (columns >> 1)) - 1) * 0x1111;

    // Remove the bits outside the image edge.
    for (i = 0; i < TX_32X32; i++) {
      lfm->left_y[i] &= mask_y;
      lfm->above_y[i] &= mask_y;
      lfm->left_uv[i] &= mask_uv;
      lfm->above_uv[i] &= mask_uv;
    }
    lfm->int_4x4_y &= mask_y;
    lfm->int_4x4_uv &= mask_uv_int;

    // We don't apply a wide loop filter on the last uv column.  If set
    // apply the shorter one instead.
    if (columns == 1) {
      lfm->left_uv[TX_8X8] |= lfm->left_uv[TX_16X16];
      lfm->left_uv[TX_16X16] = 0;
    }
    if (columns == 5) {
      lfm->left_uv[TX_8X8] |= (lfm->left_uv[TX_16X16] & 0xcccc);
      lfm->left_uv[TX_16X16] &= ~(lfm->left_uv[TX_16X16] & 0xcccc);
    }
  }
  // We don't a loop filter on the first column in the image.  Mask that out.
  if (mi_col == 0) {
    for (i = 0; i < TX_32X32; i++) {
      lfm->left_y[i] &= 0xfefefefefefefefe;
      lfm->left_uv[i] &= 0xeeee;
    }
  }
}
#if CONFIG_NON420
static void filter_block_plane_non420(VP9_COMMON *cm,
                                      struct macroblockd_plane *plane,
                                      const MODE_INFO *mi,
                                      int mi_row, int mi_col) {
  const int ss_x = plane->subsampling_x;
  const int ss_y = plane->subsampling_y;
  const int row_step = 1 << ss_x;
  const int col_step = 1 << ss_y;
  const int row_step_stride = cm->mode_info_stride * row_step;
  struct buf_2d *const dst = &plane->dst;
  uint8_t* const dst0 = dst->buf;
  unsigned int mask_16x16[MI_BLOCK_SIZE] = {0};
  unsigned int mask_8x8[MI_BLOCK_SIZE] = {0};
  unsigned int mask_4x4[MI_BLOCK_SIZE] = {0};
  unsigned int mask_4x4_int[MI_BLOCK_SIZE] = {0};
  struct loop_filter_info lfi[MI_BLOCK_SIZE][MI_BLOCK_SIZE];
  int r, c;

  for (r = 0; r < MI_BLOCK_SIZE && mi_row + r < cm->mi_rows; r += row_step) {
    unsigned int mask_16x16_c = 0;
    unsigned int mask_8x8_c = 0;
    unsigned int mask_4x4_c = 0;
    unsigned int border_mask;

    // Determine the vertical edges that need filtering
    for (c = 0; c < MI_BLOCK_SIZE && mi_col + c < cm->mi_cols; c += col_step) {
      const int skip_this = mi[c].mbmi.skip_coeff
                            && is_inter_block(&mi[c].mbmi);
      // left edge of current unit is block/partition edge -> no skip
      const int block_edge_left = b_width_log2(mi[c].mbmi.sb_type) ?
          !(c & ((1 << (b_width_log2(mi[c].mbmi.sb_type)-1)) - 1)) : 1;
      const int skip_this_c = skip_this && !block_edge_left;
      // top edge of current unit is block/partition edge -> no skip
      const int block_edge_above = b_height_log2(mi[c].mbmi.sb_type) ?
          !(r & ((1 << (b_height_log2(mi[c].mbmi.sb_type)-1)) - 1)) : 1;
      const int skip_this_r = skip_this && !block_edge_above;
      const TX_SIZE tx_size = (plane->plane_type == PLANE_TYPE_UV)
                            ? get_uv_tx_size(&mi[c].mbmi)
                            : mi[c].mbmi.tx_size;
      const int skip_border_4x4_c = ss_x && mi_col + c == cm->mi_cols - 1;
      const int skip_border_4x4_r = ss_y && mi_row + r == cm->mi_rows - 1;

      // Filter level can vary per MI
      if (!build_lfi(&cm->lf_info, &mi[c].mbmi, lfi[r] + (c >> ss_x)))
        continue;

      // Build masks based on the transform size of each block
      if (tx_size == TX_32X32) {
        if (!skip_this_c && ((c >> ss_x) & 3) == 0) {
          if (!skip_border_4x4_c)
            mask_16x16_c |= 1 << (c >> ss_x);
          else
            mask_8x8_c |= 1 << (c >> ss_x);
        }
        if (!skip_this_r && ((r >> ss_y) & 3) == 0) {
          if (!skip_border_4x4_r)
            mask_16x16[r] |= 1 << (c >> ss_x);
          else
            mask_8x8[r] |= 1 << (c >> ss_x);
        }
      } else if (tx_size == TX_16X16) {
        if (!skip_this_c && ((c >> ss_x) & 1) == 0) {
          if (!skip_border_4x4_c)
            mask_16x16_c |= 1 << (c >> ss_x);
          else
            mask_8x8_c |= 1 << (c >> ss_x);
        }
        if (!skip_this_r && ((r >> ss_y) & 1) == 0) {
          if (!skip_border_4x4_r)
            mask_16x16[r] |= 1 << (c >> ss_x);
          else
            mask_8x8[r] |= 1 << (c >> ss_x);
        }
      } else {
        // force 8x8 filtering on 32x32 boundaries
        if (!skip_this_c) {
          if (tx_size == TX_8X8 || ((c >> ss_x) & 3) == 0)
            mask_8x8_c |= 1 << (c >> ss_x);
          else
            mask_4x4_c |= 1 << (c >> ss_x);
        }

        if (!skip_this_r) {
          if (tx_size == TX_8X8 || ((r >> ss_y) & 3) == 0)
            mask_8x8[r] |= 1 << (c >> ss_x);
          else
            mask_4x4[r] |= 1 << (c >> ss_x);
        }

        if (!skip_this && tx_size < TX_8X8 && !skip_border_4x4_c)
          mask_4x4_int[r] |= 1 << (c >> ss_x);
      }
    }

    // Disable filtering on the leftmost column
    border_mask = ~(mi_col == 0);
    filter_selectively_vert(dst->buf, dst->stride,
                            mask_16x16_c & border_mask,
                            mask_8x8_c & border_mask,
                            mask_4x4_c & border_mask,
                            mask_4x4_int[r], lfi[r]);
    dst->buf += 8 * dst->stride;
    mi += row_step_stride;
  }

  // Now do horizontal pass
  dst->buf = dst0;
  for (r = 0; r < MI_BLOCK_SIZE && mi_row + r < cm->mi_rows; r += row_step) {
    const int skip_border_4x4_r = ss_y && mi_row + r == cm->mi_rows - 1;
    const unsigned int mask_4x4_int_r = skip_border_4x4_r ? 0 : mask_4x4_int[r];

    filter_selectively_horiz(dst->buf, dst->stride,
                             mask_16x16[r],
                             mask_8x8[r],
                             mask_4x4[r],
                             mask_4x4_int_r, mi_row + r == 0, lfi[r]);
    dst->buf += 8 * dst->stride;
  }
}
#endif

static void filter_block_plane(VP9_COMMON *const cm,
                               struct macroblockd_plane *const plane,
                               const MODE_INFO *mi,
                               int mi_row, int mi_col,
                               LOOP_FILTER_MASK *lfm) {
  const int ss_x = plane->subsampling_x;
  const int ss_y = plane->subsampling_y;
  const int row_step = 1 << ss_x;
  const int col_step = 1 << ss_y;
  const int row_step_stride = cm->mode_info_stride * row_step;
  struct buf_2d *const dst = &plane->dst;
  uint8_t* const dst0 = dst->buf;
  unsigned int mask_4x4_int[MI_BLOCK_SIZE] = {0};
  struct loop_filter_info lfi[MI_BLOCK_SIZE][MI_BLOCK_SIZE];
  int r, c;
  int row_shift = 3 - ss_x;
  int row_mask = 0xff >> (ss_x << 2);

#define MASK_ROW(value) ((value >> (r_sampled << row_shift)) & row_mask)

  for (r = 0; r < MI_BLOCK_SIZE && mi_row + r < cm->mi_rows; r += row_step) {
    int r_sampled = r >> ss_x;

    // Determine the vertical edges that need filtering
    for (c = 0; c < MI_BLOCK_SIZE && mi_col + c < cm->mi_cols; c += col_step) {
      if (!build_lfi(&cm->lf_info, &mi[c].mbmi, lfi[r] + (c >> ss_x)))
        continue;
    }
    if (!plane->plane_type) {
      mask_4x4_int[r] = MASK_ROW(lfm->int_4x4_y);
      // Disable filtering on the leftmost column
      filter_selectively_vert(dst->buf, dst->stride,
                              MASK_ROW(lfm->left_y[TX_16X16]),
                              MASK_ROW(lfm->left_y[TX_8X8]),
                              MASK_ROW(lfm->left_y[TX_4X4]),
                              MASK_ROW(lfm->int_4x4_y),
                              lfi[r]);
    } else {
      mask_4x4_int[r] = MASK_ROW(lfm->int_4x4_uv);
      // Disable filtering on the leftmost column
      filter_selectively_vert(dst->buf, dst->stride,
                              MASK_ROW(lfm->left_uv[TX_16X16]),
                              MASK_ROW(lfm->left_uv[TX_8X8]),
                              MASK_ROW(lfm->left_uv[TX_4X4]),
                              MASK_ROW(lfm->int_4x4_uv),
                              lfi[r]);
    }
    dst->buf += 8 * dst->stride;
    mi += row_step_stride;
  }

  // Now do horizontal pass
  dst->buf = dst0;
  for (r = 0; r < MI_BLOCK_SIZE && mi_row + r < cm->mi_rows; r += row_step) {
    const int skip_border_4x4_r = ss_y && mi_row + r == cm->mi_rows - 1;
    const unsigned int mask_4x4_int_r = skip_border_4x4_r ? 0 : mask_4x4_int[r];
    int r_sampled = r >> ss_x;

    if (!plane->plane_type) {
      filter_selectively_horiz(dst->buf, dst->stride,
                               MASK_ROW(lfm->above_y[TX_16X16]),
                               MASK_ROW(lfm->above_y[TX_8X8]),
                               MASK_ROW(lfm->above_y[TX_4X4]),
                               MASK_ROW(lfm->int_4x4_y),
                               mi_row + r == 0, lfi[r]);
    } else {
      filter_selectively_horiz(dst->buf, dst->stride,
                               MASK_ROW(lfm->above_uv[TX_16X16]),
                               MASK_ROW(lfm->above_uv[TX_8X8]),
                               MASK_ROW(lfm->above_uv[TX_4X4]),
                               mask_4x4_int_r,
                               mi_row + r == 0, lfi[r]);
    }
    dst->buf += 8 * dst->stride;
  }
#undef MASK_ROW
}

void vp9_loop_filter_rows(const YV12_BUFFER_CONFIG *frame_buffer,
                          VP9_COMMON *cm, MACROBLOCKD *xd,
                          int start, int stop, int y_only) {
  const int num_planes = y_only ? 1 : MAX_MB_PLANE;
  int mi_row, mi_col;
  LOOP_FILTER_MASK lfm;
#if CONFIG_NON420
  int use_420 = y_only || (xd->plane[1].subsampling_y == 1 &&
      xd->plane[1].subsampling_x == 1);
#endif

  for (mi_row = start; mi_row < stop; mi_row += MI_BLOCK_SIZE) {
    MODE_INFO* const mi = cm->mi + mi_row * cm->mode_info_stride;

    for (mi_col = 0; mi_col < cm->mi_cols; mi_col += MI_BLOCK_SIZE) {
      int plane;

      setup_dst_planes(xd, frame_buffer, mi_row, mi_col);

      // TODO(JBB): Make setup_mask work for non 420.
#if CONFIG_NON420
      if (use_420)
#endif
        setup_mask(cm, mi_row, mi_col, mi + mi_col, cm->mode_info_stride, &lfm);

      for (plane = 0; plane < num_planes; ++plane) {
#if CONFIG_NON420
        if (use_420)
#endif
          filter_block_plane(cm, &xd->plane[plane], mi + mi_col, mi_row, mi_col,
                             &lfm);
#if CONFIG_NON420
        else
          filter_block_plane_non420(cm, &xd->plane[plane], mi + mi_col,
                                    mi_row, mi_col);
#endif
      }
    }
  }
}

void vp9_loop_filter_frame(VP9_COMMON *cm, MACROBLOCKD *xd,
                           int frame_filter_level,
                           int y_only, int partial) {
  int start_mi_row, end_mi_row, mi_rows_to_filter;
  if (!frame_filter_level) return;

  start_mi_row = 0;
  mi_rows_to_filter = cm->mi_rows;
  if (partial && cm->mi_rows > 8) {
    start_mi_row = cm->mi_rows >> 1;
    start_mi_row &= 0xfffffff8;
    mi_rows_to_filter = MAX(cm->mi_rows / 8, 8);
  }
  end_mi_row = start_mi_row + mi_rows_to_filter;
  vp9_loop_filter_frame_init(cm, frame_filter_level);
  vp9_loop_filter_rows(cm->frame_to_show, cm, xd,
                       start_mi_row, end_mi_row,
                       y_only);
}

int vp9_loop_filter_worker(void *arg1, void *arg2) {
  LFWorkerData *const lf_data = (LFWorkerData*)arg1;
  (void)arg2;
  vp9_loop_filter_rows(lf_data->frame_buffer, lf_data->cm, &lf_data->xd,
                       lf_data->start, lf_data->stop, lf_data->y_only);
  return 1;
}
