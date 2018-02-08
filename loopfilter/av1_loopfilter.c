/*
 * Copyright (c) 2016, Alliance for Open Media. All rights reserved
 *
 * This source code is subject to the terms of the BSD 2 Clause License and
 * the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
 * was not distributed with this source code in the LICENSE file, you can
 * obtain it at www.aomedia.org/license/software. If the Alliance for Open
 * Media Patent License 1.0 was not distributed with this source code in the
 * PATENTS file, you can obtain it at www.aomedia.org/license/patent.
 */

#include <math.h>

#include "./aom_config.h"
#include "./aom_dsp_rtcd.h"
#include "aom_dsp/aom_dsp_common.h"
#include "aom_mem/aom_mem.h"
#include "aom_ports/mem.h"
#include "av1/common/av1_loopfilter.h"
#include "av1/common/onyxc_int.h"
#include "av1/common/reconinter.h"
#include "av1/common/seg_common.h"

#if CONFIG_LOOPFILTER_LEVEL
static const SEG_LVL_FEATURES seg_lvl_lf_lut[MAX_MB_PLANE][2] = {
  { SEG_LVL_ALT_LF_Y_V, SEG_LVL_ALT_LF_Y_H },
  { SEG_LVL_ALT_LF_U, SEG_LVL_ALT_LF_U },
  { SEG_LVL_ALT_LF_V, SEG_LVL_ALT_LF_V }
};

#if CONFIG_EXT_DELTA_Q
static const int delta_lf_id_lut[MAX_MB_PLANE][2] = {
  { 0, 1 }, { 2, 2 }, { 3, 3 }
};
#endif  // CONFIG_EXT_DELTA_Q
#endif  // CONFIG_LOOPFILTER_LEVEL

#if CONFIG_DEBLOCK_13TAP
#define PARALLEL_DEBLOCKING_5_TAP_CHROMA 1
#else
#define PARALLEL_DEBLOCKING_5_TAP_CHROMA 0
#endif

#if PARALLEL_DEBLOCKING_5_TAP_CHROMA
extern void aom_lpf_vertical_6_c(uint8_t *s, int pitch, const uint8_t *blimit,
                                 const uint8_t *limit, const uint8_t *thresh);

extern void aom_lpf_horizontal_6_c(uint8_t *s, int p, const uint8_t *blimit,
                                   const uint8_t *limit, const uint8_t *thresh);

extern void aom_highbd_lpf_horizontal_6_c(uint16_t *s, int p,
                                          const uint8_t *blimit,
                                          const uint8_t *limit,
                                          const uint8_t *thresh, int bd);

extern void aom_highbd_lpf_vertical_6_c(uint16_t *s, int pitch,
                                        const uint8_t *blimit,
                                        const uint8_t *limit,
                                        const uint8_t *thresh, int bd);
#endif

// 64 bit masks for left transform size. Each 1 represents a position where
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

// 64 bit masks for above transform size. Each 1 represents a position where
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

// 64 bit masks for prediction sizes (left). Each 1 represents a position
// where left border of an 8x8 block. These are aligned to the right most
// appropriate bit, and then shifted into place.
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
static const FilterMaskY left_txform_mask[TX_SIZES] = {
  { { 0x0000000000000001ULL,  // TX_4X4,
      0x0000000000000000ULL, 0x0000000000000000ULL, 0x0000000000000000ULL } },

  { { 0x0000000000010001ULL,  // TX_8X8,
      0x0000000000000000ULL, 0x0000000000000000ULL, 0x0000000000000000ULL } },

  { { 0x0001000100010001ULL,  // TX_16X16,
      0x0000000000000000ULL, 0x0000000000000000ULL, 0x0000000000000000ULL } },

  { { 0x0001000100010001ULL,  // TX_32X32,
      0x0001000100010001ULL, 0x0000000000000000ULL, 0x0000000000000000ULL } },
#if CONFIG_TX64X64
  { { 0x0001000100010001ULL,  // TX_64X64,
      0x0001000100010001ULL, 0x0001000100010001ULL, 0x0001000100010001ULL } },
#endif
};

// 64 bit mask to shift and set for each prediction size.
static const FilterMaskY above_txform_mask[TX_SIZES] = {
  { { 0x0000000000000001ULL,  // TX_4X4
      0x0000000000000000ULL, 0x0000000000000000ULL, 0x0000000000000000ULL } },

  { { 0x0000000000000003ULL,  // TX_8X8
      0x0000000000000000ULL, 0x0000000000000000ULL, 0x0000000000000000ULL } },

  { { 0x000000000000000fULL,  // TX_16X16
      0x0000000000000000ULL, 0x0000000000000000ULL, 0x0000000000000000ULL } },

  { { 0x00000000000000ffULL,  // TX_32X32
      0x0000000000000000ULL, 0x0000000000000000ULL, 0x0000000000000000ULL } },
#if CONFIG_TX64X64
  { { 0x000000000000ffffULL,  // TX_64X64
      0x0000000000000000ULL, 0x0000000000000000ULL, 0x0000000000000000ULL } },
#endif
};

// 64 bit mask to shift and set for each prediction size. A bit is set for
// each 8x8 block that would be in the top left most block of the given block
// size in the 64x64 block.
static const FilterMaskY size_mask[BLOCK_SIZES_ALL] = {
  { { 0x0000000000000001ULL,  // BLOCK_4X4
      0x0000000000000000ULL, 0x0000000000000000ULL, 0x0000000000000000ULL } },

  { { 0x0000000000010001ULL,  // BLOCK_4X8
      0x0000000000000000ULL, 0x0000000000000000ULL, 0x0000000000000000ULL } },

  { { 0x0000000000000003ULL,  // BLOCK_8X4
      0x0000000000000000ULL, 0x0000000000000000ULL, 0x0000000000000000ULL } },

  { { 0x0000000000030003ULL,  // BLOCK_8X8
      0x0000000000000000ULL, 0x0000000000000000ULL, 0x0000000000000000ULL } },

  { { 0x0003000300030003ULL,  // BLOCK_8X16
      0x0000000000000000ULL, 0x0000000000000000ULL, 0x0000000000000000ULL } },

  { { 0x00000000000f000fULL,  // BLOCK_16X8
      0x0000000000000000ULL, 0x0000000000000000ULL, 0x0000000000000000ULL } },

  { { 0x000f000f000f000fULL,  // BLOCK_16X16
      0x0000000000000000ULL, 0x0000000000000000ULL, 0x0000000000000000ULL } },

  { { 0x000f000f000f000fULL,  // BLOCK_16X32
      0x000f000f000f000fULL, 0x0000000000000000ULL, 0x0000000000000000ULL } },

  { { 0x00ff00ff00ff00ffULL,  // BLOCK_32X16
      0x0000000000000000ULL, 0x0000000000000000ULL, 0x0000000000000000ULL } },

  { { 0x00ff00ff00ff00ffULL,  // BLOCK_32X32
      0x00ff00ff00ff00ffULL, 0x0000000000000000ULL, 0x0000000000000000ULL } },

  { { 0x00ff00ff00ff00ffULL,  // BLOCK_32X64
      0x00ff00ff00ff00ffULL, 0x00ff00ff00ff00ffULL, 0x00ff00ff00ff00ffULL } },

  { { 0xffffffffffffffffULL,  // BLOCK_64X32
      0xffffffffffffffffULL, 0x0000000000000000ULL, 0x0000000000000000ULL } },

  { { 0xffffffffffffffffULL,  // BLOCK_64X64
      0xffffffffffffffffULL, 0xffffffffffffffffULL, 0xffffffffffffffffULL } },

#if CONFIG_EXT_PARTITION
  // TODO(luoyi): We don't process larger than 64x64 block
  // BLOCK_64X128
  { { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } },
  // BLOCK_128X64
  { { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } },
  // BLOCK_128X128
  { { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } },
#endif

  { { 0x0001000100010001ULL,  // BLOCK_4X16
      0x0000000000000000ULL, 0x0000000000000000ULL, 0x0000000000000000ULL } },

  { { 0x000000000000000fULL,  // BLOCK_16X4
      0x0000000000000000ULL, 0x0000000000000000ULL, 0x0000000000000000ULL } },

  { { 0x0003000300030003ULL,  // BLOCK_8X32
      0x0003000300030003ULL, 0x0000000000000000ULL, 0x0000000000000000ULL } },

  { { 0x0000000000ff00ffULL,  // BLOCK_32X8
      0x0000000000000000ULL, 0x0000000000000000ULL, 0x0000000000000000ULL } },

  { { 0x000f000f000f000fULL,  // BLOCK_16X64
      0x000f000f000f000fULL, 0x000f000f000f000fULL, 0x000f000f000f000fULL } },

  { { 0xffffffffffffffffULL,  // BLOCK_64X16
      0x0000000000000000ULL, 0x0000000000000000ULL, 0x0000000000000000ULL } },
#if CONFIG_EXT_PARTITION
  // TODO(luoyi): We don't process larger than 64x64 block
  // BLOCK_32X128
  { { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } },
  // BLOCK_128X32
  { { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } },
#endif
};

#if 0
// These are used for masking the left and above 32x32 borders.
static const FilterMaskY left_border = { {
    0x0101010101010101ULL, 0x0101010101010101ULL, 0x0101010101010101ULL,
    0x0101010101010101ULL,
} };

static const FilterMaskY above_border = { {
    0x000000000000ffffULL, 0x0000000000000000ULL, 0x000000000000ffffULL,
    0x0000000000000000ULL,
} };
#endif

// Note: follow 32x32 geometry size
// 64 bit left mask to shift and set for each uv prediction size.
static const FilterMaskUV left_txform_mask_uv[TX_SIZES_ALL] = {
  0x0000000000000001ULL,  // TX_4X4
  0x0000000000000101ULL,  // TX_8X8
  0x0000000001010101ULL,  // TX_16X16
  0x0101010101010101ULL,  // TX_32X32
#if CONFIG_TX64X64
  0x0000000000000000ULL,  // TX_64X64
#endif
  0x0000000000000001ULL,  // TX_4X8
  0x0000000000000001ULL,  // TX_8X4
  0x0000000000000101ULL,  // TX_8X16
  0x0000000000000001ULL,  // TX_16X8
  0x0000000001010101ULL,  // TX_16X32
  0x0000000000000101ULL,  // TX_32X16
#if CONFIG_TX64X64
  0x0000000000000000ULL,  // TX_32X64
  0x0000000000000000ULL,  // TX_64X32
#endif
  0x0000000000000101ULL,  // TX_4X16
  0x0000000000000001ULL,  // TX_16X4
  0x0000000001010101ULL,  // TX_8X32
  0x0000000000000001ULL,  // TX_32X8
#if CONFIG_TX64X64
  0x0000000000000000ULL,  // TX_16X64
  0x0000000000000000ULL,  // TX_64X16
#endif
};

// Note: from 32x32 geometry size point of view
// 16 bit above mask to shift and set for uv each prediction size.
static const FilterMaskUV above_txform_mask_uv[TX_SIZES_ALL] = {
  0x0000000000000001ULL,  // TX_4X4
  0x0000000000000003ULL,  // TX_8X8
  0x000000000000000fULL,  // TX_16X16
  0x00000000000000ffULL,  // TX_32X32
#if CONFIG_TX64X64
  0x0000000000000000ULL,  // TX_64X64
#endif
  0x0000000000000001ULL,  // TX_4X8
  0x0000000000000001ULL,  // TX_8X4
  0x0000000000000001ULL,  // TX_8X16
  0x0000000000000003ULL,  // TX_16X8
  0x0000000000000003ULL,  // TX_16X32
  0x000000000000000fULL,  // TX_32X16
#if CONFIG_TX64X64
  0x0000000000000000ULL,  // TX_32X64
  0x0000000000000000ULL,  // TX_64X32
#endif
  0x0000000000000001ULL,  // TX_4X16
  0x0000000000000003ULL,  // TX_16X4
  0x0000000000000001ULL,  // TX_8X32
  0x000000000000000fULL,  // TX_32X8
#if CONFIG_TX64X64
  0x0000000000000000ULL,  // TX_16X64
  0x0000000000000000ULL,  // TX_64X16
#endif
};

// Note: From the point of view, 64x64 geometry size shrinking to 32x32
// 64 bit mask to shift and set for each uv prediction size
static const FilterMaskUV size_mask_uv[BLOCK_SIZES_ALL] = {
  0x0000000000000001ULL,  // BLOCK_4X4
  0x0000000000000101ULL,  // BLOCK_4X8
  0x0000000000000003ULL,  // BLOCK_8X4
  0x0000000000000303ULL,  // BLOCK_8X8
  0x0000000003030303ULL,  // BLOCK_8X16,
  0x0000000000000f0fULL,  // BLOCK_16X8
  0x000000000f0f0f0fULL,  // BLOCK_16X16
  0x0f0f0f0f0f0f0f0fULL,  // BLOCK_16X32,
  0x00000000ffffffffULL,  // BLOCK_32X16,
  0xffffffffffffffffULL,  // BLOCK_32X32,
  0xffffffffffffffffULL,  // BLOCK_32X64,
  0xffffffffffffffffULL,  // BLOCK_64X32,
  0xffffffffffffffffULL,  // BLOCK_64X64,
#if CONFIG_EXT_PARTITION
  0xffffffffffffffffULL,  // BLOCK_64X128,
  0xffffffffffffffffULL,  // BLOCK_128X64,
  0xffffffffffffffffULL,  // BLOCK_128X128,
#endif
  0x0000000001010101ULL,  // BLOCK_4X16,
  0x000000000000000fULL,  // BLOCK_16X4,
  0x0303030303030303ULL,  // BLOCK_8X32,
  0x000000000000ffffULL,  // BLOCK_32X8,
  0x0f0f0f0f0f0f0f0fULL,  // BLOCK_16X64,
  0x00000000ffffffffULL,  // BLOCK_64X16
#if CONFIG_EXT_PARTITION
  0xffffffffffffffffULL,  // BLOCK_32X128,
  0xffffffffffffffffULL,  // BLOCK_128X32,
#endif
};

#if 0
static const uint64_t left_border_uv = 0x0101010101010101ULL;
static const uint64_t above_border_uv = 0x00000000000000ffULL;
#endif

static const int mode_lf_lut[] = {
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  // INTRA_MODES
  1, 1, 0, 1,                             // INTER_MODES (GLOBALMV == 0)
  1, 1, 1, 1, 1, 1, 0, 1  // INTER_COMPOUND_MODES (GLOBAL_GLOBALMV == 0)
};

static void update_sharpness(loop_filter_info_n *lfi, int sharpness_lvl) {
  int lvl;

  // For each possible value for the loop filter fill out limits
  for (lvl = 0; lvl <= MAX_LOOP_FILTER; lvl++) {
    // Set loop filter parameters that control sharpness.
    int block_inside_limit = lvl >> ((sharpness_lvl > 0) + (sharpness_lvl > 4));

    if (sharpness_lvl > 0) {
      if (block_inside_limit > (9 - sharpness_lvl))
        block_inside_limit = (9 - sharpness_lvl);
    }

    if (block_inside_limit < 1) block_inside_limit = 1;

    memset(lfi->lfthr[lvl].lim, block_inside_limit, SIMD_WIDTH);
    memset(lfi->lfthr[lvl].mblim, (2 * (lvl + 2) + block_inside_limit),
           SIMD_WIDTH);
  }
}
#if CONFIG_EXT_DELTA_Q
static uint8_t get_filter_level(const AV1_COMMON *cm,
                                const loop_filter_info_n *lfi_n,
#if CONFIG_LOOPFILTER_LEVEL
                                const int dir_idx, int plane,
#endif
                                const MB_MODE_INFO *mbmi) {
  const int segment_id = mbmi->segment_id;
  if (cm->delta_lf_present_flag) {
#if CONFIG_LOOPFILTER_LEVEL
    int delta_lf;
    if (cm->delta_lf_multi) {
      const int delta_lf_idx = delta_lf_id_lut[plane][dir_idx];
      delta_lf = mbmi->curr_delta_lf[delta_lf_idx];
    } else {
      delta_lf = mbmi->current_delta_lf_from_base;
    }
    int lvl_seg =
        clamp(delta_lf + cm->lf.filter_level[dir_idx], 0, MAX_LOOP_FILTER);
#else
    int lvl_seg = clamp(mbmi->current_delta_lf_from_base + cm->lf.filter_level,
                        0, MAX_LOOP_FILTER);
#endif
#if CONFIG_LOOPFILTER_LEVEL
    assert(plane >= 0 && plane <= 2);
    const int seg_lf_feature_id = seg_lvl_lf_lut[plane][dir_idx];
    if (segfeature_active(&cm->seg, segment_id, seg_lf_feature_id)) {
      const int data = get_segdata(&cm->seg, segment_id, seg_lf_feature_id);
      lvl_seg = clamp(lvl_seg + data, 0, MAX_LOOP_FILTER);
    }
#else
    if (segfeature_active(&cm->seg, segment_id, SEG_LVL_ALT_LF)) {
      const int data = get_segdata(&cm->seg, segment_id, SEG_LVL_ALT_LF);
      lvl_seg = clamp(lvl_seg + data, 0, MAX_LOOP_FILTER);
    }
#endif  // CONFIG_LOOPFILTER_LEVEL

    if (cm->lf.mode_ref_delta_enabled) {
      const int scale = 1 << (lvl_seg >> 5);
      lvl_seg += cm->lf.ref_deltas[mbmi->ref_frame[0]] * scale;
      if (mbmi->ref_frame[0] > INTRA_FRAME)
        lvl_seg += cm->lf.mode_deltas[mode_lf_lut[mbmi->mode]] * scale;
      lvl_seg = clamp(lvl_seg, 0, MAX_LOOP_FILTER);
    }
    return lvl_seg;
  } else {
#if CONFIG_LOOPFILTER_LEVEL
    return lfi_n
        ->lvl[segment_id][dir_idx][mbmi->ref_frame[0]][mode_lf_lut[mbmi->mode]];
#else
    return lfi_n->lvl[segment_id][mbmi->ref_frame[0]][mode_lf_lut[mbmi->mode]];
#endif
  }
}
#else
static uint8_t get_filter_level(const loop_filter_info_n *lfi_n,
                                const MB_MODE_INFO *mbmi) {
  const int segment_id = mbmi->segment_id;
  return lfi_n->lvl[segment_id][mbmi->ref_frame[0]][mode_lf_lut[mbmi->mode]];
}
#endif

void av1_loop_filter_init(AV1_COMMON *cm) {
  assert(MB_MODE_COUNT == NELEMENTS(mode_lf_lut));
  loop_filter_info_n *lfi = &cm->lf_info;
  struct loopfilter *lf = &cm->lf;
  int lvl;

  // init limits for given sharpness
  update_sharpness(lfi, lf->sharpness_level);
  lf->last_sharpness_level = lf->sharpness_level;

  // init hev threshold const vectors
  for (lvl = 0; lvl <= MAX_LOOP_FILTER; lvl++)
    memset(lfi->lfthr[lvl].hev_thr, (lvl >> 4), SIMD_WIDTH);
}

void av1_loop_filter_frame_init(AV1_COMMON *cm, int default_filt_lvl,
                                int default_filt_lvl_r
#if CONFIG_LOOPFILTER_LEVEL
                                ,
                                int plane
#endif
                                ) {
  int seg_id;
  // n_shift is the multiplier for lf_deltas
  // the multiplier is 1 for when filter_lvl is between 0 and 31;
  // 2 when filter_lvl is between 32 and 63
  loop_filter_info_n *const lfi = &cm->lf_info;
  struct loopfilter *const lf = &cm->lf;
  const struct segmentation *const seg = &cm->seg;

  // update limits if sharpness has changed
  if (lf->last_sharpness_level != lf->sharpness_level) {
    update_sharpness(lfi, lf->sharpness_level);
    lf->last_sharpness_level = lf->sharpness_level;
  }

  for (seg_id = 0; seg_id < MAX_SEGMENTS; seg_id++) {
    for (int dir = 0; dir < 2; ++dir) {
      int lvl_seg = (dir == 0) ? default_filt_lvl : default_filt_lvl_r;
#if CONFIG_LOOPFILTER_LEVEL
      assert(plane >= 0 && plane <= 2);
      const int seg_lf_feature_id = seg_lvl_lf_lut[plane][dir];
      if (segfeature_active(seg, seg_id, seg_lf_feature_id)) {
        const int data = get_segdata(&cm->seg, seg_id, seg_lf_feature_id);
        lvl_seg = clamp(lvl_seg + data, 0, MAX_LOOP_FILTER);
      }
#else
      if (segfeature_active(seg, seg_id, SEG_LVL_ALT_LF)) {
        const int data = get_segdata(seg, seg_id, SEG_LVL_ALT_LF);
        lvl_seg = clamp(lvl_seg + data, 0, MAX_LOOP_FILTER);
      }
#endif  // CONFIG_LOOPFILTER_LEVEL

      if (!lf->mode_ref_delta_enabled) {
// we could get rid of this if we assume that deltas are set to
// zero when not in use; encoder always uses deltas
#if CONFIG_LOOPFILTER_LEVEL
        memset(lfi->lvl[seg_id][dir], lvl_seg, sizeof(lfi->lvl[seg_id][dir]));
#else
        memset(lfi->lvl[seg_id], lvl_seg, sizeof(lfi->lvl[seg_id]));
#endif  // CONFIG_LOOPFILTER_LEVEL
      } else {
        int ref, mode;
#if CONFIG_LOOPFILTER_LEVEL
        const int scale = 1 << (lvl_seg >> 5);
        const int intra_lvl = lvl_seg + lf->ref_deltas[INTRA_FRAME] * scale;
        lfi->lvl[seg_id][dir][INTRA_FRAME][0] =
            clamp(intra_lvl, 0, MAX_LOOP_FILTER);

        for (ref = LAST_FRAME; ref < TOTAL_REFS_PER_FRAME; ++ref) {
          for (mode = 0; mode < MAX_MODE_LF_DELTAS; ++mode) {
            const int inter_lvl = lvl_seg + lf->ref_deltas[ref] * scale +
                                  lf->mode_deltas[mode] * scale;
            lfi->lvl[seg_id][dir][ref][mode] =
                clamp(inter_lvl, 0, MAX_LOOP_FILTER);
          }
        }
#else
        const int scale = 1 << (default_filt_lvl >> 5);
        const int intra_lvl = lvl_seg + lf->ref_deltas[INTRA_FRAME] * scale;
        lfi->lvl[seg_id][INTRA_FRAME][0] = clamp(intra_lvl, 0, MAX_LOOP_FILTER);

        for (ref = LAST_FRAME; ref < TOTAL_REFS_PER_FRAME; ++ref) {
          for (mode = 0; mode < MAX_MODE_LF_DELTAS; ++mode) {
            const int inter_lvl = lvl_seg + lf->ref_deltas[ref] * scale +
                                  lf->mode_deltas[mode] * scale;
            lfi->lvl[seg_id][ref][mode] = clamp(inter_lvl, 0, MAX_LOOP_FILTER);
          }
        }
#endif
      }
    }
  }

  // For each loopfilter operation per frame, we init our neighbor info
  uint8_t *infop = cm->lf.neighbor;
  const size_t unit_size = cm->lf.neighbor_width + cm->lf.neighbor_height;
  // Y/UV tx_size
  memset(infop, TX_64X64, 2 * unit_size);
  infop += 2 * unit_size;
  // level YUV and skip
  memset(infop, 0, 4 * unit_size);

  // TODO(luoyi): Need a way to reset loop bitmask per frame
  // if (cm->lf.curr_frame_offset != cm->current_video_frame) {
  //   cm->lf.curr_frame_offset = cm->current_video_frame;
  //   memset(cm->lf.lfm, 0, cm->lf.lfm_num * sizeof(*cm->lf.lfm));
  // }
}

typedef void (*LpfFunc)(uint8_t *s, int p, const uint8_t *blimit,
                        const uint8_t *limit, const uint8_t *thresh);

typedef void (*LpfDualFunc)(uint8_t *s, int p, const uint8_t *blimit0,
                            const uint8_t *limit0, const uint8_t *thresh0,
                            const uint8_t *blimit1, const uint8_t *limit1,
                            const uint8_t *thresh1);

typedef void (*HbdLpfFunc)(uint16_t *s, int p, const uint8_t *blimit,
                           const uint8_t *limit, const uint8_t *thresh, int bd);

typedef void (*HbdLpfDualFunc)(uint16_t *s, int p, const uint8_t *blimit0,
                               const uint8_t *limit0, const uint8_t *thresh0,
                               const uint8_t *blimit1, const uint8_t *limit1,
                               const uint8_t *thresh1, int bd);

static void filter_selectively_vert_row2(
    int subsampling_factor, uint8_t *s, int plane, int pitch,
    unsigned int mask_16x16_l, unsigned int mask_8x8_l, unsigned int mask_4x4_l,
    const loop_filter_info_n *lfi_n, const uint8_t *lfl) {
  const int mask_shift = subsampling_factor ? 8 : 16;
  const int mask_cutoff = subsampling_factor ? 0xff : 0xffff;
  const int lfl_forward = subsampling_factor ? 8 : 16;

  unsigned int mask_16x16_0 = mask_16x16_l & mask_cutoff;
  unsigned int mask_8x8_0 = mask_8x8_l & mask_cutoff;
  unsigned int mask_4x4_0 = mask_4x4_l & mask_cutoff;

  unsigned int mask_16x16_1 = (mask_16x16_l >> mask_shift) & mask_cutoff;
  unsigned int mask_8x8_1 = (mask_8x8_l >> mask_shift) & mask_cutoff;
  unsigned int mask_4x4_1 = (mask_4x4_l >> mask_shift) & mask_cutoff;

  unsigned int mask;

  LpfFunc lpf_vertical_16 = NULL;
  LpfFunc lpf_vertical_8 = NULL;

  if (plane != 0) {
    lpf_vertical_16 = aom_lpf_vertical_6_c;
    lpf_vertical_8 = aom_lpf_vertical_6_c;
  } else {
    lpf_vertical_16 = aom_lpf_vertical_16;
    lpf_vertical_8 = aom_lpf_vertical_8;
  }

  for (mask = mask_16x16_0 | mask_8x8_0 | mask_4x4_0 | mask_16x16_1 |
              mask_8x8_1 | mask_4x4_1;
       mask; mask >>= 1) {
    const loop_filter_thresh *lfi0 = lfi_n->lfthr + *lfl;
    const loop_filter_thresh *lfi1 = lfi_n->lfthr + *(lfl + lfl_forward);

    if (mask & 1) {
      if ((mask_16x16_0 | mask_16x16_1) & 1) {
        if ((mask_16x16_0 & mask_16x16_1) & 1) {
// TODO(luoyi): aom_lpf_vertical_16_dual() SIMD version is not ready
// yet; use C code here now.
#if 1
          aom_lpf_vertical_16_c(s, pitch, lfi0->mblim, lfi0->lim,
                                lfi0->hev_thr);
          aom_lpf_vertical_16_c(s + 4 * pitch, pitch, lfi0->mblim, lfi0->lim,
                                lfi0->hev_thr);
#else
          aom_lpf_vertical_16_dual(s, pitch, lfi0->mblim, lfi0->lim,
                                   lfi0->hev_thr);
#endif
        } else if (mask_16x16_0 & 1) {
          lpf_vertical_16(s, pitch, lfi0->mblim, lfi0->lim, lfi0->hev_thr);
        } else {
          lpf_vertical_16(s + 4 * pitch, pitch, lfi1->mblim, lfi1->lim,
                          lfi1->hev_thr);
        }
      }

      if ((mask_8x8_0 | mask_8x8_1) & 1) {
        if ((mask_8x8_0 & mask_8x8_1) & 1) {
          aom_lpf_vertical_8_dual(s, pitch, lfi0->mblim, lfi0->lim,
                                  lfi0->hev_thr, lfi1->mblim, lfi1->lim,
                                  lfi1->hev_thr);
        } else if (mask_8x8_0 & 1) {
          lpf_vertical_8(s, pitch, lfi0->mblim, lfi0->lim, lfi0->hev_thr);
        } else {
          lpf_vertical_8(s + 4 * pitch, pitch, lfi1->mblim, lfi1->lim,
                         lfi1->hev_thr);
        }
      }

      if ((mask_4x4_0 | mask_4x4_1) & 1) {
        if ((mask_4x4_0 & mask_4x4_1) & 1) {
          aom_lpf_vertical_4_dual(s, pitch, lfi0->mblim, lfi0->lim,
                                  lfi0->hev_thr, lfi1->mblim, lfi1->lim,
                                  lfi1->hev_thr);
        } else if (mask_4x4_0 & 1) {
          aom_lpf_vertical_4(s, pitch, lfi0->mblim, lfi0->lim, lfi0->hev_thr);
        } else {
          aom_lpf_vertical_4(s + 4 * pitch, pitch, lfi1->mblim, lfi1->lim,
                             lfi1->hev_thr);
        }
      }
    }

    s += 4;
    lfl += 1;
    mask_16x16_0 >>= 1;
    mask_8x8_0 >>= 1;
    mask_4x4_0 >>= 1;

    mask_16x16_1 >>= 1;
    mask_8x8_1 >>= 1;
    mask_4x4_1 >>= 1;
  }
}

static void highbd_filter_selectively_vert_row2(
    int subsampling_factor, uint16_t *s, int plane, int pitch,
    unsigned int mask_16x16_l, unsigned int mask_8x8_l, unsigned int mask_4x4_l,
    const loop_filter_info_n *lfi_n, const uint8_t *lfl, int bd) {
  const int mask_shift = subsampling_factor ? 8 : 16;
  const int mask_cutoff = subsampling_factor ? 0xff : 0xffff;
  const int lfl_forward = subsampling_factor ? 8 : 16;

  unsigned int mask_16x16_0 = mask_16x16_l & mask_cutoff;
  unsigned int mask_8x8_0 = mask_8x8_l & mask_cutoff;
  unsigned int mask_4x4_0 = mask_4x4_l & mask_cutoff;

  unsigned int mask_16x16_1 = (mask_16x16_l >> mask_shift) & mask_cutoff;
  unsigned int mask_8x8_1 = (mask_8x8_l >> mask_shift) & mask_cutoff;
  unsigned int mask_4x4_1 = (mask_4x4_l >> mask_shift) & mask_cutoff;

  unsigned int mask;

  HbdLpfFunc lpf_vertical_16 = NULL;
  HbdLpfFunc lpf_vertical_8 = NULL;

  if (plane != 0) {
    lpf_vertical_16 = aom_highbd_lpf_vertical_6_c;
    lpf_vertical_8 = aom_highbd_lpf_vertical_6_c;
  } else {
    lpf_vertical_16 = aom_highbd_lpf_vertical_16;
    lpf_vertical_8 = aom_highbd_lpf_vertical_8;
  }

  for (mask = mask_16x16_0 | mask_8x8_0 | mask_4x4_0 | mask_16x16_1 |
              mask_8x8_1 | mask_4x4_1;
       mask; mask >>= 1) {
    const loop_filter_thresh *lfi0 = lfi_n->lfthr + *lfl;
    const loop_filter_thresh *lfi1 = lfi_n->lfthr + *(lfl + lfl_forward);

    if (mask & 1) {
      if ((mask_16x16_0 | mask_16x16_1) & 1) {
        if ((mask_16x16_0 & mask_16x16_1) & 1) {
// TODO(luoyi): aom_highbd_lpf_vertical_16_dual() SIMD version is not
// ready yet; use C code here now.
#if 1
          aom_highbd_lpf_vertical_16_c(s, pitch, lfi0->mblim, lfi0->lim,
                                       lfi0->hev_thr, bd);
          aom_highbd_lpf_vertical_16_c(s + 4 * pitch, pitch, lfi0->mblim,
                                       lfi0->lim, lfi0->hev_thr, bd);
#else
          aom_highbd_lpf_vertical_16_dual(s, pitch, lfi0->mblim, lfi0->lim,
                                          lfi0->hev_thr, bd);
#endif
        } else if (mask_16x16_0 & 1) {
          lpf_vertical_16(s, pitch, lfi0->mblim, lfi0->lim, lfi0->hev_thr, bd);
        } else {
          lpf_vertical_16(s + 4 * pitch, pitch, lfi1->mblim, lfi1->lim,
                          lfi1->hev_thr, bd);
        }
      }

      if ((mask_8x8_0 | mask_8x8_1) & 1) {
        if ((mask_8x8_0 & mask_8x8_1) & 1) {
          aom_highbd_lpf_vertical_8_dual(s, pitch, lfi0->mblim, lfi0->lim,
                                         lfi0->hev_thr, lfi1->mblim, lfi1->lim,
                                         lfi1->hev_thr, bd);
        } else if (mask_8x8_0 & 1) {
          lpf_vertical_8(s, pitch, lfi0->mblim, lfi0->lim, lfi0->hev_thr, bd);
        } else {
          lpf_vertical_8(s + 4 * pitch, pitch, lfi1->mblim, lfi1->lim,
                         lfi1->hev_thr, bd);
        }
      }

      if ((mask_4x4_0 | mask_4x4_1) & 1) {
        if ((mask_4x4_0 & mask_4x4_1) & 1) {
          aom_highbd_lpf_vertical_4_dual(s, pitch, lfi0->mblim, lfi0->lim,
                                         lfi0->hev_thr, lfi1->mblim, lfi1->lim,
                                         lfi1->hev_thr, bd);
        } else if (mask_4x4_0 & 1) {
          aom_highbd_lpf_vertical_4(s, pitch, lfi0->mblim, lfi0->lim,
                                    lfi0->hev_thr, bd);
        } else {
          aom_highbd_lpf_vertical_4(s + 4 * pitch, pitch, lfi1->mblim,
                                    lfi1->lim, lfi1->hev_thr, bd);
        }
      }
    }

    s += 4;
    lfl += 1;
    mask_16x16_0 >>= 1;
    mask_8x8_0 >>= 1;
    mask_4x4_0 >>= 1;

    mask_16x16_1 >>= 1;
    mask_8x8_1 >>= 1;
    mask_4x4_1 >>= 1;
  }
}

static void filter_selectively_horiz(
    uint8_t *s, int pitch, unsigned int mask_16x16, unsigned int mask_8x8,
    unsigned int mask_4x4, unsigned int mask_4x4_int,
    const loop_filter_info_n *lfi_n, const uint8_t *lfl) {
  unsigned int mask;
  int count;

  for (mask = mask_16x16 | mask_8x8 | mask_4x4 | mask_4x4_int; mask;
       mask >>= count) {
    const loop_filter_thresh *lfi = lfi_n->lfthr + *lfl;

    count = 1;
    if (mask & 1) {
      if (mask_16x16 & 1) {
        if ((mask_16x16 & 3) == 3) {
          aom_lpf_horizontal_16_dual(s, pitch, lfi->mblim, lfi->lim,
                                     lfi->hev_thr);
          count = 2;
        } else {
          aom_lpf_horizontal_16(s, pitch, lfi->mblim, lfi->lim, lfi->hev_thr);
        }
      } else if (mask_8x8 & 1) {
        if ((mask_8x8 & 3) == 3) {
          // Next block's thresholds.
          const loop_filter_thresh *lfin = lfi_n->lfthr + *(lfl + 1);

          aom_lpf_horizontal_8_dual(s, pitch, lfi->mblim, lfi->lim,
                                    lfi->hev_thr, lfin->mblim, lfin->lim,
                                    lfin->hev_thr);

          if ((mask_4x4_int & 3) == 3) {
            aom_lpf_horizontal_4_dual(s + 4 * pitch, pitch, lfi->mblim,
                                      lfi->lim, lfi->hev_thr, lfin->mblim,
                                      lfin->lim, lfin->hev_thr);
          } else {
            if (mask_4x4_int & 1)
              aom_lpf_horizontal_4(s + 4 * pitch, pitch, lfi->mblim, lfi->lim,
                                   lfi->hev_thr);
            else if (mask_4x4_int & 2)
              aom_lpf_horizontal_4(s + 8 + 4 * pitch, pitch, lfin->mblim,
                                   lfin->lim, lfin->hev_thr);
          }
          count = 2;
        } else {
          aom_lpf_horizontal_8(s, pitch, lfi->mblim, lfi->lim, lfi->hev_thr);

          if (mask_4x4_int & 1)
            aom_lpf_horizontal_4(s + 4 * pitch, pitch, lfi->mblim, lfi->lim,
                                 lfi->hev_thr);
        }
      } else if (mask_4x4 & 1) {
        if ((mask_4x4 & 3) == 3) {
          // Next block's thresholds.
          const loop_filter_thresh *lfin = lfi_n->lfthr + *(lfl + 1);

          aom_lpf_horizontal_4_dual(s, pitch, lfi->mblim, lfi->lim,
                                    lfi->hev_thr, lfin->mblim, lfin->lim,
                                    lfin->hev_thr);

          if ((mask_4x4_int & 3) == 3) {
            aom_lpf_horizontal_4_dual(s + 4 * pitch, pitch, lfi->mblim,
                                      lfi->lim, lfi->hev_thr, lfin->mblim,
                                      lfin->lim, lfin->hev_thr);
          } else {
            if (mask_4x4_int & 1)
              aom_lpf_horizontal_4(s + 4 * pitch, pitch, lfi->mblim, lfi->lim,
                                   lfi->hev_thr);
            else if (mask_4x4_int & 2)
              aom_lpf_horizontal_4(s + 8 + 4 * pitch, pitch, lfin->mblim,
                                   lfin->lim, lfin->hev_thr);
          }
          count = 2;
        } else {
          aom_lpf_horizontal_4(s, pitch, lfi->mblim, lfi->lim, lfi->hev_thr);

          if (mask_4x4_int & 1)
            aom_lpf_horizontal_4(s + 4 * pitch, pitch, lfi->mblim, lfi->lim,
                                 lfi->hev_thr);
        }
      } else if (mask_4x4_int & 1) {
        aom_lpf_horizontal_4(s + 4 * pitch, pitch, lfi->mblim, lfi->lim,
                             lfi->hev_thr);
      }
    }
    s += 8 * count;
    lfl += count;
    mask_16x16 >>= count;
    mask_8x8 >>= count;
    mask_4x4 >>= count;
    mask_4x4_int >>= count;
  }
}

static void highbd_filter_selectively_horiz(
    uint16_t *s, int pitch, unsigned int mask_16x16, unsigned int mask_8x8,
    unsigned int mask_4x4, unsigned int mask_4x4_int,
    const loop_filter_info_n *lfi_n, const uint8_t *lfl, int bd) {
  unsigned int mask;
  int count;

  for (mask = mask_16x16 | mask_8x8 | mask_4x4 | mask_4x4_int; mask;
       mask >>= count) {
    const loop_filter_thresh *lfi = lfi_n->lfthr + *lfl;

    count = 1;
    if (mask & 1) {
      if (mask_16x16 & 1) {
        if ((mask_16x16 & 3) == 3) {
          aom_highbd_lpf_horizontal_16_dual(s, pitch, lfi->mblim, lfi->lim,
                                            lfi->hev_thr, bd);
          count = 2;
        } else {
          aom_highbd_lpf_horizontal_16(s, pitch, lfi->mblim, lfi->lim,
                                       lfi->hev_thr, bd);
        }
      } else if (mask_8x8 & 1) {
        if ((mask_8x8 & 3) == 3) {
          // Next block's thresholds.
          const loop_filter_thresh *lfin = lfi_n->lfthr + *(lfl + 1);

          aom_highbd_lpf_horizontal_8_dual(s, pitch, lfi->mblim, lfi->lim,
                                           lfi->hev_thr, lfin->mblim, lfin->lim,
                                           lfin->hev_thr, bd);

          if ((mask_4x4_int & 3) == 3) {
            aom_highbd_lpf_horizontal_4_dual(
                s + 4 * pitch, pitch, lfi->mblim, lfi->lim, lfi->hev_thr,
                lfin->mblim, lfin->lim, lfin->hev_thr, bd);
          } else {
            if (mask_4x4_int & 1) {
              aom_highbd_lpf_horizontal_4(s + 4 * pitch, pitch, lfi->mblim,
                                          lfi->lim, lfi->hev_thr, bd);
            } else if (mask_4x4_int & 2) {
              aom_highbd_lpf_horizontal_4(s + 8 + 4 * pitch, pitch, lfin->mblim,
                                          lfin->lim, lfin->hev_thr, bd);
            }
          }
          count = 2;
        } else {
          aom_highbd_lpf_horizontal_8(s, pitch, lfi->mblim, lfi->lim,
                                      lfi->hev_thr, bd);

          if (mask_4x4_int & 1) {
            aom_highbd_lpf_horizontal_4(s + 4 * pitch, pitch, lfi->mblim,
                                        lfi->lim, lfi->hev_thr, bd);
          }
        }
      } else if (mask_4x4 & 1) {
        if ((mask_4x4 & 3) == 3) {
          // Next block's thresholds.
          const loop_filter_thresh *lfin = lfi_n->lfthr + *(lfl + 1);

          aom_highbd_lpf_horizontal_4_dual(s, pitch, lfi->mblim, lfi->lim,
                                           lfi->hev_thr, lfin->mblim, lfin->lim,
                                           lfin->hev_thr, bd);
          if ((mask_4x4_int & 3) == 3) {
            aom_highbd_lpf_horizontal_4_dual(
                s + 4 * pitch, pitch, lfi->mblim, lfi->lim, lfi->hev_thr,
                lfin->mblim, lfin->lim, lfin->hev_thr, bd);
          } else {
            if (mask_4x4_int & 1) {
              aom_highbd_lpf_horizontal_4(s + 4 * pitch, pitch, lfi->mblim,
                                          lfi->lim, lfi->hev_thr, bd);
            } else if (mask_4x4_int & 2) {
              aom_highbd_lpf_horizontal_4(s + 8 + 4 * pitch, pitch, lfin->mblim,
                                          lfin->lim, lfin->hev_thr, bd);
            }
          }
          count = 2;
        } else {
          aom_highbd_lpf_horizontal_4(s, pitch, lfi->mblim, lfi->lim,
                                      lfi->hev_thr, bd);

          if (mask_4x4_int & 1) {
            aom_highbd_lpf_horizontal_4(s + 4 * pitch, pitch, lfi->mblim,
                                        lfi->lim, lfi->hev_thr, bd);
          }
        }
      } else if (mask_4x4_int & 1) {
        aom_highbd_lpf_horizontal_4(s + 4 * pitch, pitch, lfi->mblim, lfi->lim,
                                    lfi->hev_thr, bd);
      }
    }
    s += 8 * count;
    lfl += count;
    mask_16x16 >>= count;
    mask_8x8 >>= count;
    mask_4x4 >>= count;
    mask_4x4_int >>= count;
  }
}

static void filter_selectively_horiz2(uint8_t *s, int plane, int pitch,
                                      unsigned int mask_16x16,
                                      unsigned int mask_8x8,
                                      unsigned int mask_4x4,
                                      const loop_filter_info_n *lfi_n,
                                      const uint8_t *lfl) {
  unsigned int mask;
  int count;

  LpfFunc lpf_horizontal_16 = NULL;
  LpfFunc lpf_horizontal_8 = NULL;

  if (plane != 0) {
    lpf_horizontal_16 = aom_lpf_horizontal_6_c;
    lpf_horizontal_8 = aom_lpf_horizontal_6_c;
  } else {
    lpf_horizontal_16 = aom_lpf_horizontal_16;
    lpf_horizontal_8 = aom_lpf_horizontal_8;
  }

  for (mask = mask_16x16 | mask_8x8 | mask_4x4; mask; mask >>= count) {
    const loop_filter_thresh *lfi = lfi_n->lfthr + *lfl;

    count = 1;
    if (mask & 1) {
      if (mask_16x16 & 1) {
        if ((mask_16x16 & 3) == 3) {
// TODO(luoyi): aom_lpf_horizontal_16_dual() SIMD version is not ready
// yet; use C code here now.
#if 1
          aom_lpf_horizontal_16_c(s, pitch, lfi->mblim, lfi->lim, lfi->hev_thr);
          aom_lpf_horizontal_16_c(s + 4, pitch, lfi->mblim, lfi->lim,
                                  lfi->hev_thr);
#else
          aom_lpf_horizontal_16_dual(s, pitch, lfi->mblim, lfi->lim,
                                     lfi->hev_thr);
#endif
          count = 2;
        } else {
          lpf_horizontal_16(s, pitch, lfi->mblim, lfi->lim, lfi->hev_thr);
        }
      } else if (mask_8x8 & 1) {
        if ((mask_8x8 & 3) == 3) {
          // Next block's thresholds.
          const loop_filter_thresh *lfin = lfi_n->lfthr + *(lfl + 1);

          aom_lpf_horizontal_8_dual(s, pitch, lfi->mblim, lfi->lim,
                                    lfi->hev_thr, lfin->mblim, lfin->lim,
                                    lfin->hev_thr);
          count = 2;
        } else {
          lpf_horizontal_8(s, pitch, lfi->mblim, lfi->lim, lfi->hev_thr);
        }
      } else if (mask_4x4 & 1) {
        if ((mask_4x4 & 3) == 3) {
          // Next block's thresholds.
          const loop_filter_thresh *lfin = lfi_n->lfthr + *(lfl + 1);

          aom_lpf_horizontal_4_dual(s, pitch, lfi->mblim, lfi->lim,
                                    lfi->hev_thr, lfin->mblim, lfin->lim,
                                    lfin->hev_thr);
          count = 2;
        } else {
          aom_lpf_horizontal_4(s, pitch, lfi->mblim, lfi->lim, lfi->hev_thr);
        }
      }
    }

    s += 4 * count;
    lfl += count;
    mask_16x16 >>= count;
    mask_8x8 >>= count;
    mask_4x4 >>= count;
  }
}

static void highbd_filter_selectively_horiz2(uint16_t *s, int plane, int pitch,
                                             unsigned int mask_16x16,
                                             unsigned int mask_8x8,
                                             unsigned int mask_4x4,
                                             const loop_filter_info_n *lfi_n,
                                             const uint8_t *lfl, int bd) {
  unsigned int mask;
  int count;

  HbdLpfFunc lpf_horizontal_16 = NULL;
  HbdLpfFunc lpf_horizontal_8 = NULL;

  if (plane != 0) {
    lpf_horizontal_16 = aom_highbd_lpf_horizontal_6_c;
    lpf_horizontal_8 = aom_highbd_lpf_horizontal_6_c;
  } else {
    lpf_horizontal_16 = aom_highbd_lpf_horizontal_16;
    lpf_horizontal_8 = aom_highbd_lpf_horizontal_8;
  }

  for (mask = mask_16x16 | mask_8x8 | mask_4x4; mask; mask >>= count) {
    const loop_filter_thresh *lfi = lfi_n->lfthr + *lfl;

    count = 1;
    if (mask & 1) {
      if (mask_16x16 & 1) {
        if ((mask_16x16 & 3) == 3) {
// TODO(luoyi): aom_highbd_lpf_horizontal_16_dual() SIMD version is
// not ready yet; use C code here now.
#if 1
          aom_highbd_lpf_horizontal_16_c(s, pitch, lfi->mblim, lfi->lim,
                                         lfi->hev_thr, bd);
          aom_highbd_lpf_horizontal_16_c(s + 4, pitch, lfi->mblim, lfi->lim,
                                         lfi->hev_thr, bd);
#else
          aom_highbd_lpf_horizontal_16_dual(s, pitch, lfi->mblim, lfi->lim,
                                            lfi->hev_thr, bd);
#endif
          count = 2;
        } else {
          lpf_horizontal_16(s, pitch, lfi->mblim, lfi->lim, lfi->hev_thr, bd);
        }
      } else if (mask_8x8 & 1) {
        if ((mask_8x8 & 3) == 3) {
          // Next block's thresholds.
          const loop_filter_thresh *lfin = lfi_n->lfthr + *(lfl + 1);

          aom_highbd_lpf_horizontal_8_dual(s, pitch, lfi->mblim, lfi->lim,
                                           lfi->hev_thr, lfin->mblim, lfin->lim,
                                           lfin->hev_thr, bd);
          count = 2;
        } else {
          lpf_horizontal_8(s, pitch, lfi->mblim, lfi->lim, lfi->hev_thr, bd);
        }
      } else if (mask_4x4 & 1) {
        if ((mask_4x4 & 3) == 3) {
          // Next block's thresholds.
          const loop_filter_thresh *lfin = lfi_n->lfthr + *(lfl + 1);

          aom_highbd_lpf_horizontal_4_dual(s, pitch, lfi->mblim, lfi->lim,
                                           lfi->hev_thr, lfin->mblim, lfin->lim,
                                           lfin->hev_thr, bd);
          count = 2;
        } else {
          aom_highbd_lpf_horizontal_4(s, pitch, lfi->mblim, lfi->lim,
                                      lfi->hev_thr, bd);
        }
      }
    }

    s += 4 * count;
    lfl += count;
    mask_16x16 >>= count;
    mask_8x8 >>= count;
    mask_4x4 >>= count;
  }
}

static INLINE TX_SIZE get_y_tx_size(const MODE_INFO *mi) {
  const MB_MODE_INFO *const mbmi = &mi->mbmi;
  TX_SIZE tx_size = mbmi->tx_size;
  if (is_inter_block(mbmi) && !mbmi->skip) {
    tx_size = mbmi->inter_tx_size[0][0];
  }
  assert(tx_size != TX_INVALID);
  return tx_size;
}

// Note:
// Establish 64x64 block, contructed by 256 (16x16) 4x4 sub-block.
// Every 4 rows would be represented by one uint64_t mask. Hence,
// there are 4 uint64_t bitmask[4] to represent the whole 64x64.
//
// Given a location by (idx, idy), This function returns the index
// 0, 1, 2, 3 to select which bitmask[] to use.
// Then the pointer y_shift contains the shift value in the bit mask.
// Function returns y_shift; y_index contains the index.
//
static int get_y_index_shift(int idx, int idy, int *y_index) {
  *y_index = idy >> 4;
  const int y_idy = (idy >> 2) % 4;
  return (y_idy << 4) + (idx >> 2);
}

// Note:
// For 4:2:0 format sampling, establish 32x32 block, constructed by
// 64 (8x8), 4x4 sub-block. We need one uint64_t bitmask to present
// all edge information
// Function returns uv_shift.
//
static int get_uv_index_shift(int idx, int idy) {
  return ((idy >> 3) << 3) + (idx >> 3);
}

static int get_uv_shift(int idx, int idy) {
  return (((idy - 2) >> 2) << 3) + (idx >> 2);
}

static void check_mask_left_y(const LoopFilterMask *lfm) {
#ifndef NDEBUG
  int i;
  for (i = 0; i < 4; ++i) {
    assert(!(lfm->left_y[TX_32X32].bits[i] & lfm->left_y[TX_8X8].bits[i]));
    assert(!(lfm->left_y[TX_4X4].bits[i] & lfm->left_y[TX_8X8].bits[i]));
    assert(!(lfm->left_y[TX_4X4].bits[i] & lfm->left_y[TX_16X16].bits[i]));
    assert(!(lfm->left_y[TX_8X8].bits[i] & lfm->left_y[TX_16X16].bits[i]));
  }
#else
  (void)lfm;
#endif
}

static void setup_masks(AV1_COMMON *const cm, const MB_MODE_INFO *const mbmi,
                        TX_SIZE y_tx_size, int idx, int idy, int mi_row,
                        int mi_col, int lflf_h, int lflf_v,
                        LoopFilterMask *lfm) {
  const TX_SIZE tx_wide = txsize_horz_map[y_tx_size];
  const TX_SIZE tx_high = txsize_vert_map[y_tx_size];

  const int tx_unit_h = tx_size_high_unit[y_tx_size];
  const int tx_unit_w = tx_size_wide_unit[y_tx_size];

  const int mi_row_pos = mi_row + (idy >> MI_SIZE_LOG2);
  const int mi_col_pos = mi_col + (idx >> MI_SIZE_LOG2);

  uint8_t *top_tx = cm->lf.neighbor;
  uint8_t *left_tx = top_tx + cm->lf.neighbor_width;

  const BLOCK_SIZE block_size = mbmi->sb_type;

  int acc_h = 0;
  FilterMaskY *left = NULL;
  int index = 0;
  int shift_y = 0;
  const int lfm_bound = 1 << (MAX_SB_SIZE_LOG2 - 1);

  TX_SIZE ntx;
  uint64_t mask[4];
  int x, y, i;

  y = idy;
  x = (idx >= lfm_bound) ? idx - lfm_bound : idx;
  int mi_row_loc = mi_row_pos;

  if (lflf_v) {
    while (acc_h < tx_unit_h) {
      ntx = left_tx[mi_row_loc];
      const int ntx_unit_h = tx_size_high_unit[ntx];

      TX_SIZE tx_w = txsize_horz_map[ntx];
      TX_SIZE tx_h = txsize_vert_map[ntx];

      TX_SIZE choice_w;
      if (tx_w < tx_wide) {
        choice_w = tx_w;
      } else {
        choice_w = tx_wide;
      }

      left = &lfm->left_y[choice_w];
      y = (y >= lfm_bound) ? y - lfm_bound : y;
      shift_y = get_y_index_shift(x, y, &index);

      const TX_SIZE choice_h = AOMMIN(tx_h, tx_high);
      const int tx_inc = AOMMIN(ntx_unit_h, tx_unit_h);

      if (choice_h < TX_32X32) {
        mask[0] =
            left_txform_mask[choice_h].bits[0] & size_mask[block_size].bits[0];
        mask[0] <<= shift_y;

        assert(!(left->bits[index] & mask[0]));
        left->bits[index] |= mask[0];
      } else if (choice_h == TX_32X32) {
        for (i = 0; i < 2; ++i) {
          mask[i] = left_txform_mask[choice_h].bits[i] &
                    size_mask[block_size].bits[i];
          mask[i] <<= shift_y;

          assert(!(left->bits[index + i] & mask[i]));
          left->bits[index + i] |= mask[i];
        }
      } else {
        for (i = 0; i < 4; ++i) {
          mask[i] = left_txform_mask[choice_h].bits[i] &
                    size_mask[block_size].bits[i];
          mask[i] <<= shift_y;

          assert(!(left->bits[i] & mask[i]));
          left->bits[i] |= mask[i];
        }
      }
      // Increment accumulated height
      acc_h += tx_inc;
      y += (tx_inc << MI_SIZE_LOG2);
      // Increment row position
      mi_row_loc += tx_inc;

      if (mi_row_loc >= cm->mi_rows) break;

      check_mask_left_y(lfm);
    }
  }

  int acc_w = 0;
  FilterMaskY *above = NULL;
  x = idx;
  y = (idy >= lfm_bound) ? idy - lfm_bound : idy;
  int mi_col_loc = mi_col_pos;

  if (lflf_h) {
    while (acc_w < tx_unit_w) {
      ntx = top_tx[mi_col_loc];
      const int ntx_unit_w = tx_size_wide_unit[ntx];

      TX_SIZE tx_w = txsize_horz_map[ntx];
      TX_SIZE tx_h = txsize_vert_map[ntx];

      TX_SIZE choice_h;
      if (tx_h < tx_high) {
        choice_h = tx_h;
      } else {
        choice_h = tx_high;
      }

      x = (x >= lfm_bound) ? x - lfm_bound : x;
      shift_y = get_y_index_shift(x, y, &index);

      const TX_SIZE choice_w = AOMMIN(tx_w, tx_wide);
      const int tx_inc = AOMMIN(ntx_unit_w, tx_unit_w);

      above = &lfm->above_y[choice_h];
      mask[0] =
          above_txform_mask[choice_w].bits[0] & size_mask[block_size].bits[0];
      mask[0] <<= shift_y;
      assert(!(above->bits[index] & mask[0]));
      above->bits[index] |= mask[0];

      // Increment accumulated width
      acc_w += tx_inc;
      x += (tx_inc << MI_SIZE_LOG2);
      // Increment col position
      mi_col_loc += tx_inc;

      if (mi_col_loc >= cm->mi_cols) break;
    }
  }

  memset(&top_tx[mi_col_pos], y_tx_size, tx_unit_w);
  memset(&left_tx[mi_row_pos], y_tx_size, tx_unit_h);
}

static void setup_mask_from_tx_size(AV1_COMMON *const cm, MODE_INFO **mi,
                                    int mi_stride, int max_tx_wide,
                                    int max_tx_high, int idx, int idy,
                                    int mi_row, int mi_col, LoopFilterMask *fm,
                                    int br, int bc, int lflf_h, int lflf_v) {
  TX_SIZE tx_size = get_y_tx_size(mi[0]);
  int tx_wide = tx_size_wide[tx_size];
  int tx_high = tx_size_high[tx_size];

  int i, x, y, x_off, y_off;
  int sub_tx_w, sub_tx_h;
  MODE_INFO **mi2 = NULL;

  if (tx_wide == max_tx_wide && tx_high == max_tx_high) {
    const MB_MODE_INFO *const mbmi = &mi[0]->mbmi;
    LoopFilterMask *lfm = fm + ((idy >> 6) << 1) + (idx >> 6);
    setup_masks(cm, mbmi, tx_size, idx, idy, mi_row, mi_col, lflf_h, lflf_v,
                lfm);
  } else if (tx_wide == max_tx_wide) {
    sub_tx_h = max_tx_high > tx_high ? max_tx_high >> 1 : max_tx_high << 1;
    sub_tx_w = max_tx_wide;
    for (i = 0; i < 4; ++i) {
      x_off = (i & 1) * sub_tx_w;
      y_off = (i >> 1) * sub_tx_h;

      if (x_off >= max_tx_wide || y_off >= max_tx_high) continue;

      x = idx + x_off;
      y = idy + y_off;
      if ((x >> 2) > bc || (y >> 2) > br) continue;

      mi2 = mi + (y_off >> MI_SIZE_LOG2) * mi_stride + (x_off >> MI_SIZE_LOG2);
      setup_mask_from_tx_size(cm, mi2, mi_stride, sub_tx_w, sub_tx_h, x, y,
                              mi_row, mi_col, fm, br, bc, lflf_h, lflf_v);
    }
  } else if (tx_high == max_tx_high) {
    sub_tx_h = max_tx_high;
    sub_tx_w = max_tx_wide > tx_wide ? max_tx_wide >> 1 : max_tx_wide << 1;
    for (i = 0; i < 4; ++i) {
      x_off = (i & 1) * sub_tx_w;
      y_off = (i >> 1) * sub_tx_h;

      if (x_off >= max_tx_wide || y_off >= max_tx_high) break;

      x = idx + x_off;
      y = idy + y_off;
      if ((x >> 2) > bc || (y >> 2) > br) continue;

      mi2 = mi + (y_off >> MI_SIZE_LOG2) * mi_stride + (x_off >> MI_SIZE_LOG2);
      setup_mask_from_tx_size(cm, mi2, mi_stride, sub_tx_w, sub_tx_h, x, y,
                              mi_row, mi_col, fm, br, bc, lflf_h, lflf_v);
    }
  } else {
    sub_tx_w = max_tx_wide >> 1;
    sub_tx_h = max_tx_high >> 1;

    int count = 4;
    for (i = 0; i < count; ++i) {
      x_off = (i & 1) * sub_tx_w;
      y_off = (i >> 1) * sub_tx_h;

      y = idy + y_off;
      x = idx + x_off;
      if ((x >> 2) > bc || (y >> 2) > br) continue;

      mi2 = mi + (y_off >> MI_SIZE_LOG2) * mi_stride + (x_off >> MI_SIZE_LOG2);

      tx_size = get_y_tx_size(mi2[0]);
      int w = tx_size_wide[tx_size];
      int h = tx_size_high[tx_size];
      if (h == sub_tx_h && w == sub_tx_w << 1) count -= 2;
      if (w == sub_tx_w && h == sub_tx_h << 1) count -= 1;

      setup_mask_from_tx_size(cm, mi2, mi_stride, sub_tx_w, sub_tx_h, x, y,
                              mi_row, mi_col, fm, br, bc, lflf_h, lflf_v);
    }
  }
}

static uint8_t receive_filter_level(AV1_COMMON *const cm,
                                    const loop_filter_info_n *const lfi_n,
                                    int dir_idx, int plane, int mi_row,
                                    int mi_col, const MB_MODE_INFO *mbmi) {
#if CONFIG_EXT_DELTA_Q
#if CONFIG_LOOPFILTER_LEVEL
  (void)mi_row;
  (void)mi_col;
  const uint8_t filter_level =
      get_filter_level(cm, lfi_n, dir_idx, plane, mbmi);
#else
  (void)dir_idx;
  (void)plane;
  const uint8_t filter_level = get_filter_level(cm, lfi_n, mbmi);
#endif  // CONFIG_LOOPFILTER_LEVEL
#else
  (void)cm;
  (void)dir_idx;
  (void)plane;
  (void)mi_row;
  (void)mi_col;
  const uint8_t filter_level = get_filter_level(lfi_n, mbmi);
#endif  // CONFIG_EXT_DELTA_Q
  return filter_level;
}

// We specify block_size instead of using mbmi->sb_type
// unless block_size = BLOCK_INVALID
static void get_y_level_64x64(AV1_COMMON *const cm, MODE_INFO **mi,
                              BLOCK_SIZE block_size, int idx, int idy,
                              int mi_row, int mi_col, LoopFilterMask *lfm,
                              int *lfl_h_flag, int *lfl_v_flag) {
  const MB_MODE_INFO *mbmi = &mi[0]->mbmi;
  if (block_size == BLOCK_INVALID) block_size = mbmi->sb_type;
  const loop_filter_info_n *const lfi_n = &cm->lf_info;
  const uint8_t lpl_v =
      receive_filter_level(cm, lfi_n, 0, 0, mi_row, mi_col, mbmi);
  const uint8_t lpl_h =
      receive_filter_level(cm, lfi_n, 1, 0, mi_row, mi_col, mbmi);

  const int wide = num_4x4_blocks_wide_lookup[block_size];
  const int high = num_4x4_blocks_high_lookup[block_size];

  const int lfm_bound = 1 << (MAX_SB_SIZE_LOG2 - 1);

  // convert position into one of 2x2 loopfilter masks
  const int x = (idx >= lfm_bound) ? idx - lfm_bound : idx;
  const int y = (idy >= lfm_bound) ? idy - lfm_bound : idy;
  const int row = y >> MI_SIZE_LOG2;
  const int col = x >> MI_SIZE_LOG2;

  const int mi_row_pos = mi_row + (idy >> MI_SIZE_LOG2);
  const int mi_col_pos = mi_col + (idx >> MI_SIZE_LOG2);

  uint8_t *top_level =
      cm->lf.neighbor + 2 * (cm->lf.neighbor_width + cm->lf.neighbor_height);
  uint8_t *left_level = top_level + cm->lf.neighbor_width;

  *lfl_h_flag = top_level[mi_col_pos] || lpl_h;
  *lfl_v_flag = left_level[mi_row_pos] || lpl_v;

  uint8_t lpl_h_f, lpl_v_f;
  if (!lpl_h)
    lpl_h_f = top_level[mi_col_pos] ? top_level[mi_col_pos] : 0;
  else
    lpl_h_f = lpl_h;

  if (!lpl_v)
    lpl_v_f = left_level[mi_row_pos] ? left_level[mi_row_pos] : 0;
  else
    lpl_v_f = lpl_v;

  // save for next coding block
  memset(&left_level[mi_row_pos], lpl_v, high);
  memset(&top_level[mi_col_pos], lpl_h, wide);

  int i;
  for (i = 0; i < high; i++) {
    memset(&lfm->lfl_y_hor[row + i][col], lpl_h_f, wide);
    memset(&lfm->lfl_y_ver[row + i][col], lpl_v_f, wide);
  }
}

// TODO(luoyi):
// We obtain loopfilter level flags in multiple times, whick returned by
// pointers of get_y_level_64x64() are consistent across the whole coding block.
static void process_y_level(AV1_COMMON *const cm, MODE_INFO **mi, int idx,
                            int idy, int mi_row, int mi_col,
                            LoopFilterMask *lfm, int *lfl_h_flag,
                            int *lfl_v_flag) {
  const MB_MODE_INFO *mbmi = &mi[0]->mbmi;
  const BLOCK_SIZE block_size = mbmi->sb_type;
  const int mi_stride = cm->mi_stride;
  const int block_64x64_step = 64;
  const int mi_64x64_step = block_64x64_step >> MI_SIZE_LOG2;
  MODE_INFO **mip = mi;
  int x = idx;
  int y = idy;
  int row = mi_row;
  int col = mi_col;
  LoopFilterMask *fm = lfm;

  switch (block_size) {
    case BLOCK_128X128:
      // top left
      get_y_level_64x64(cm, mip, BLOCK_64X64, x, y, row, col, fm, lfl_h_flag,
                        lfl_v_flag);

      // top right
      x = idx + block_64x64_step;
      y = idy;
      col = mi_col + mi_64x64_step;
      row = mi_row;
      fm += 1;
      mip = mi + mi_64x64_step;
      get_y_level_64x64(cm, mip, BLOCK_64X64, x, y, row, col, fm, lfl_h_flag,
                        lfl_v_flag);

      // bottom left
      x = idx;
      y = idy + block_64x64_step;
      col = mi_col;
      row = mi_row + mi_64x64_step;
      fm += 1;
      mip = mi + (mi_64x64_step * mi_stride);
      get_y_level_64x64(cm, mip, BLOCK_64X64, x, y, row, col, fm, lfl_h_flag,
                        lfl_v_flag);

      // bottom right
      x = idx + block_64x64_step;
      y = idy + block_64x64_step;
      col = mi_col + mi_64x64_step;
      row = mi_row + mi_64x64_step;
      fm += 1;
      mip = mi + (mi_64x64_step * mi_stride) + mi_64x64_step;
      get_y_level_64x64(cm, mip, BLOCK_64X64, x, y, row, col, fm, lfl_h_flag,
                        lfl_v_flag);
      break;
    case BLOCK_128X64:
      // top left
      get_y_level_64x64(cm, mip, BLOCK_64X64, x, y, row, col, fm, lfl_h_flag,
                        lfl_v_flag);

      // top right
      x = idx + block_64x64_step;
      y = idy;
      col = mi_col + mi_64x64_step;
      row = mi_row;
      fm += 1;
      mip = mi + mi_64x64_step;
      get_y_level_64x64(cm, mip, BLOCK_64X64, x, y, row, col, fm, lfl_h_flag,
                        lfl_v_flag);
      break;
    case BLOCK_64X128:
      // top left
      get_y_level_64x64(cm, mip, BLOCK_64X64, x, y, row, col, fm, lfl_h_flag,
                        lfl_v_flag);

      // bottom left
      x = idx;
      y = idy + block_64x64_step;
      col = mi_col;
      row = mi_row + mi_64x64_step;
      fm += 2;
      mip = mi + (mi_64x64_step * mi_stride);
      get_y_level_64x64(cm, mip, BLOCK_64X64, x, y, row, col, fm, lfl_h_flag,
                        lfl_v_flag);
      break;
    default:
      fm += ((idy >> 6) << 1) + (idx >> 6);
      get_y_level_64x64(cm, mi, BLOCK_INVALID, idx, idy, mi_row, mi_col, fm,
                        lfl_h_flag, lfl_v_flag);
  }
}

// This function ors into the current lfm structure, where to do loop
// filters for the specific mi we are looking at. It uses information
// including the block_size_type (32x16, 32x32, etc.), the transform size,
// whether there were any coefficients encoded, and the loop filter strength
// block we are currently looking at. Shift is used to position the
// 1's we produce.
// ------------------------------------------------------------
//                 top zone      left zone
//             neighbor_width  neighbor_height
// Y  tx_size  |--------------|---------------|
// UV tx_size  |--------------|---------------|
// Y level     |--------------|---------------|
// U level     |--------------|---------------|
// V level     |--------------|---------------|
// skip        |--------------|---------------|
// ------------------------------------------------------------
static void build_y_masks(AV1_COMMON *const cm, MODE_INFO **mi,
                          int mode_info_stride, int idx, int idy,
                          LoopFilterMask *lfm, int mi_row, int mi_col, int br,
                          int bc, int lflf_h, int lflf_v) {
  const MB_MODE_INFO *mbmi = &mi[0]->mbmi;
  const BLOCK_SIZE block_size = mbmi->sb_type;

  const int wide = num_4x4_blocks_wide_lookup[block_size];
  const int high = num_4x4_blocks_high_lookup[block_size];
  const int mi_row_pos = mi_row + (idy >> MI_SIZE_LOG2);
  const int mi_col_pos = mi_col + (idx >> MI_SIZE_LOG2);

  // TODO(luoyi): Skip logic can be removed?
  uint8_t *top_skip =
      cm->lf.neighbor + 5 * (cm->lf.neighbor_width + cm->lf.neighbor_height);
  uint8_t *left_skip = top_skip + cm->lf.neighbor_width;

  const uint8_t skip = is_inter_block(mbmi) && mbmi->skip;
  // const uint8_t top_skip_flag = skip && top_skip[mi_col_pos];
  // const uint8_t left_skip_flag = skip && left_skip[mi_row_pos];

  // save skip for neighbor block
  memset(&left_skip[mi_row_pos], skip, high);
  memset(&top_skip[mi_col_pos], skip, wide);

  // If the block has no coefficients and is not intra we skip applying
  // the loop filter on block edges.
  // No matter current block is skipped or not; we build vertical and
  // horizontal edge by using mbmi->tx_size. So PU edge is always
  // considered.

  const int block_wide = block_size_wide[block_size];
  const int block_high = block_size_high[block_size];
  setup_mask_from_tx_size(cm, mi, mode_info_stride, block_wide, block_high, idx,
                          idy, mi_row, mi_col, lfm, br, bc, lflf_h, lflf_v);
}

static void check_mask_left_u(const LoopFilterMask *lfm) {
#ifndef NDEBUG
  assert(!(lfm->left_u[TX_32X32] & lfm->left_u[TX_4X4]));
  assert(!(lfm->left_u[TX_32X32] & lfm->left_u[TX_8X8]));
  assert(!(lfm->left_u[TX_32X32] & lfm->left_u[TX_16X16]));
  assert(!(lfm->left_u[TX_16X16] & lfm->left_u[TX_8X8]));
  assert(!(lfm->left_u[TX_16X16] & lfm->left_u[TX_4X4]));
  assert(!(lfm->left_u[TX_8X8] & lfm->left_u[TX_4X4]));
#else
  (void)lfm;
#endif
}

static void check_mask_left_v(const LoopFilterMask *lfm) {
#ifndef NDEBUG
  assert(!(lfm->left_v[TX_32X32] & lfm->left_v[TX_4X4]));
  assert(!(lfm->left_v[TX_32X32] & lfm->left_v[TX_8X8]));
  assert(!(lfm->left_v[TX_32X32] & lfm->left_v[TX_16X16]));
  assert(!(lfm->left_v[TX_16X16] & lfm->left_v[TX_8X8]));
  assert(!(lfm->left_v[TX_16X16] & lfm->left_v[TX_4X4]));
  assert(!(lfm->left_v[TX_8X8] & lfm->left_v[TX_4X4]));
#else
  (void)lfm;
#endif
}

static void check_mask_top_u(const LoopFilterMask *lfm) {
#ifndef NDEBUG
  assert(!(lfm->above_u[TX_32X32] & lfm->above_u[TX_4X4]));
  assert(!(lfm->above_u[TX_32X32] & lfm->above_u[TX_8X8]));
  assert(!(lfm->above_u[TX_32X32] & lfm->above_u[TX_16X16]));
  assert(!(lfm->above_u[TX_16X16] & lfm->above_u[TX_8X8]));
  assert(!(lfm->above_u[TX_16X16] & lfm->above_u[TX_4X4]));
  assert(!(lfm->above_u[TX_8X8] & lfm->above_u[TX_4X4]));
#else
  (void)lfm;
#endif
}

static void check_mask_top_v(const LoopFilterMask *lfm) {
#ifndef NDEBUG
  assert(!(lfm->above_v[TX_32X32] & lfm->above_v[TX_4X4]));
  assert(!(lfm->above_v[TX_32X32] & lfm->above_v[TX_8X8]));
  assert(!(lfm->above_v[TX_32X32] & lfm->above_v[TX_16X16]));
  assert(!(lfm->above_v[TX_16X16] & lfm->above_v[TX_8X8]));
  assert(!(lfm->above_v[TX_16X16] & lfm->above_v[TX_4X4]));
  assert(!(lfm->above_v[TX_8X8] & lfm->above_v[TX_4X4]));
#else
  (void)lfm;
#endif
}

// TODO(luoyi): 4:2:0 only
static void setup_uv_masks(AV1_COMMON *const cm, MODE_INFO **mi, int plane,
                           int idx, int idy, int mi_row, int mi_col,
                           LoopFilterMask *fm) {
  const MB_MODE_INFO *mbmi = &mi[0]->mbmi;
  const BLOCK_SIZE block_size = mbmi->sb_type;
  const TX_SIZE uv_tx_size = av1_get_uv_tx_size(mbmi, 1, 1);

  const TX_SIZE tx_wide = txsize_horz_map[uv_tx_size];
  const TX_SIZE tx_high = txsize_vert_map[uv_tx_size];

  const int tx_unit_w = tx_size_wide_unit[uv_tx_size];
  const int tx_unit_h = tx_size_high_unit[uv_tx_size];

  const int mi_row_pos = (mi_row + (idy >> MI_SIZE_LOG2)) >> 1;
  const int mi_col_pos = (mi_col + (idx >> MI_SIZE_LOG2)) >> 1;

  uint8_t *top_tx =
      cm->lf.neighbor + (cm->lf.neighbor_width + cm->lf.neighbor_height);
  uint8_t *left_tx = top_tx + cm->lf.neighbor_width;

  LoopFilterMask *lfm = fm + ((idy >> 6) << 1) + (idx >> 6);

  FilterMaskUV *left_uv = NULL;
  FilterMaskUV *above_uv = NULL;

  int x, y;
  const int lfm_bound = 1 << (MAX_SB_SIZE_LOG2 - 1);

  y = idy;
  x = (idx >= lfm_bound) ? idx - lfm_bound : idx;
  int mi_row_loc = mi_row_pos;

  int shift_uv;
  uint64_t mask;
  int acc_size = 0;
  TX_SIZE ntx, tx_w, tx_h, choice_w, choice_h;

  // Left edge
  //if (mi_row == 0 && mi_col == 64 && idx == 48 && idy == 32 && uv_tx_size == TX_8X8) {
  if (mi_row == 0 && mi_col == 64 && idx == 48 && idy == 16 && uv_tx_size == TX_8X8) {
    acc_size = 0;
  }
  while (acc_size < tx_unit_h) {
    ntx = left_tx[mi_row_loc];
    const int ntx_unit_h = tx_size_high_unit[ntx];

    tx_w = txsize_horz_map[ntx];
    tx_h = txsize_vert_map[ntx];

    choice_w = tx_wide;
    if (tx_w < tx_wide) choice_w = tx_w;

    left_uv = &lfm->left_u[choice_w];
    if (plane == 2) {
      left_uv = &lfm->left_v[choice_w];
    }

    y = (y >= lfm_bound) ? y - lfm_bound : y;
    shift_uv = get_uv_index_shift(x, y);

    choice_h = AOMMIN(tx_h, tx_high);
    const int tx_inc = AOMMIN(ntx_unit_h, tx_unit_h);

    mask = left_txform_mask_uv[choice_h] & size_mask_uv[block_size];
    mask <<= shift_uv;
    *left_uv |= mask;

    const uint64_t m = (uint64_t)1 << 38;
    if (lfm->left_u[TX_8X8] & m) {
      mask = 0;
    }
    // Increment accumulated height
    acc_size += tx_inc;
    y += tx_inc << (MI_SIZE_LOG2 + 1);
    // Increment row position
    mi_row_loc += tx_inc;

    if (mi_row_loc >= cm->mi_rows >> 1) break;

    if (plane == 1) {
      check_mask_left_u(lfm);
    } else {
      check_mask_left_v(lfm);
    }
  }

  // Top edge
  acc_size = 0;
  x = idx;
  y = (idy >= lfm_bound) ? idy - lfm_bound : idy;
  int mi_col_loc = mi_col_pos;

  while (acc_size < tx_unit_w) {
    ntx = top_tx[mi_col_loc];
    const int ntx_unit_w = tx_size_wide_unit[ntx];

    tx_w = txsize_horz_map[ntx];
    tx_h = txsize_vert_map[ntx];

    choice_h = tx_high;
    if (tx_h < tx_high) choice_h = tx_h;

    above_uv = &lfm->above_u[choice_h];
    if (plane == 2) above_uv = &lfm->above_v[choice_h];

    x = (x >= lfm_bound) ? x - lfm_bound : x;
    shift_uv = get_uv_index_shift(x, y);

    choice_w = AOMMIN(tx_w, tx_wide);
    const int tx_inc = AOMMIN(ntx_unit_w, tx_unit_w);

    mask = above_txform_mask_uv[choice_w] & size_mask_uv[block_size];
    mask <<= shift_uv;
    *above_uv |= mask;

    // Increment accumulated width
    acc_size += tx_inc;
    x += tx_inc << (MI_SIZE_LOG2 + 1);
    // Increment col position
    mi_col_loc += tx_inc;

    if (mi_col_loc >= cm->mi_cols >> 1) break;

    if (plane == 1)
      check_mask_top_u(lfm);
    else
      check_mask_top_v(lfm);
  }

  // Save the real size for neighbor
  memset(&left_tx[mi_row_pos], uv_tx_size, tx_unit_h);
  memset(&top_tx[mi_col_pos], uv_tx_size, tx_unit_w);
}

static void build_uv_masks(AV1_COMMON *const cm, MODE_INFO **mi, int plane,
                           int idx, int idy, int mi_row, int mi_col,
                           LoopFilterMask *fm) {
  const int chroma_pos_step = 8;
  if (idx % chroma_pos_step == 0 && idy % chroma_pos_step == 0) {
    setup_uv_masks(cm, mi, plane, idx, idy, mi_row, mi_col, fm);
  }
}

// This function update the bit masks for the entire 64x64 region represented
// by mi_row, mi_col. In case one of the edge is a tile boundary, loop filtering
// for that edge is disabled. This function only check the tile boundary info
// for the top left corner mi to determine the boundary information for the
// top and left edge of the whole super block
static void update_tile_boundary_filter_mask(AV1_COMMON *const cm,
                                             const int mi_row, const int mi_col,
                                             LoopFilterMask *lfmask) {
#if CONFIG_LOOPFILTERING_ACROSS_TILES || CONFIG_LOOPFILTERING_ACROSS_TILES_EXT
  const MODE_INFO *const mi = cm->mi + mi_row * cm->mi_stride + mi_col;
  int count;
#if CONFIG_EXT_PARTITION
  if (mi->mbmi.sb_type == BLOCK_128X128)
    count = 2;
  else
    count = 1;
#else
  count = 1;
#endif

  int i;
  const BOUNDARY_TYPE *const bi =
      cm->boundary_info + mi_row * cm->mi_stride + mi_col;
  if (*bi & TILE_LEFT_BOUNDARY) {
    LoopFilterMask *lfm = lfmask;
    int ind = count;
    while (ind > 0) {
      for (i = 0; i <= TX_32X32; i++) {
        int j;
        for (j = 0; j < 4; ++j) {
          lfm->left_y[i].bits[j] &= 0xfffefffefffefffeULL;
        }
        lfm->left_u[i] &= 0xfefefefefefefefeULL;
        lfm->left_v[i] &= 0xfefefefefefefefeULL;
      }
      ind -= 1;
      // We work on lfm[0] and lfm[2] for Left edge
      lfm += 2;
    }
  }

  if (*bi & TILE_ABOVE_BOUNDARY) {
    LoopFilterMask *lfm = lfmask;
    int ind = count;
    while (ind > 0) {
      for (i = 0; i <= TX_32X32; i++) {
        lfm->above_y[i].bits[0] &= 0xffffffffffff0000ULL;
        lfm->above_u[i] &= 0xffffffffffffff00ULL;
        lfm->above_v[i] &= 0xffffffffffffff00ULL;
      }
      ind -= 1;
      // We work on lfm[0] and lfm[1] for top edge
      lfm += 1;
    }
  }
#else
  (void)cm;
  (void)mi_row;
  (void)mi_col;
  (void)lfm;
#endif
  // CONFIG_LOOPFILTERING_ACROSS_TILES || CONFIG_LOOPFILTERING_ACROSS_TILES_EXT
}

static void check_mask_above_y(const LoopFilterMask *lfm) {
#ifndef NDEBUG
  int i;
  for (i = 0; i < 4; ++i) {
    assert(!(lfm->above_y[TX_4X4].bits[i] & lfm->above_y[TX_8X8].bits[i]));
    assert(!(lfm->above_y[TX_4X4].bits[i] & lfm->above_y[TX_16X16].bits[i]));
    assert(!(lfm->above_y[TX_8X8].bits[i] & lfm->above_y[TX_16X16].bits[i]));
  }
#else
  (void)lfm;
#endif
}

static void check_loop_filter_masks(const LoopFilterMask *lfm) {
  int i;
  for (i = 0; i < LOOP_FILTER_MASK_NUM; ++i) {
    // Assert if we try to apply 2 different loop filters at the same position.
    assert(!(lfm->left_u[TX_16X16] & lfm->left_u[TX_8X8]));
    assert(!(lfm->left_u[TX_16X16] & lfm->left_u[TX_4X4]));
    assert(!(lfm->left_u[TX_8X8] & lfm->left_u[TX_4X4]));
    assert(!(lfm->above_u[TX_16X16] & lfm->above_u[TX_8X8]));
    assert(!(lfm->above_u[TX_16X16] & lfm->above_u[TX_4X4]));
    assert(!(lfm->above_u[TX_8X8] & lfm->above_u[TX_4X4]));

    assert(!(lfm->left_v[TX_16X16] & lfm->left_v[TX_8X8]));
    assert(!(lfm->left_v[TX_16X16] & lfm->left_v[TX_4X4]));
    assert(!(lfm->left_v[TX_8X8] & lfm->left_v[TX_4X4]));
    assert(!(lfm->above_v[TX_16X16] & lfm->above_v[TX_8X8]));
    assert(!(lfm->above_v[TX_16X16] & lfm->above_v[TX_4X4]));
    assert(!(lfm->above_v[TX_8X8] & lfm->above_v[TX_4X4]));

    check_mask_left_y(lfm);
    check_mask_above_y(lfm);
    lfm += 1;
  }
}

// ------------------------------------------------------------
//                 top zone      left zone
//             neighbor_width  neighbor_height
// Y  tx_size  |--------------|---------------|
// UV tx_size  |--------------|---------------|
// Y level     |--------------|---------------|
// U level     |--------------|---------------|
// V level     |--------------|---------------|
// skip        |--------------|---------------|
// ------------------------------------------------------------
// This function does either U or V plane
// Return loopfilter level flag is either for U or V.
//
// We specify block_size instead of using mbmi->sb_type
// unless block_size = BLOCK_INVALID
static int get_uv_level_64x64(AV1_COMMON *const cm, MODE_INFO **mi,
                              BLOCK_SIZE block_size, int plane, int idx,
                              int idy, int mi_row, int mi_col,
                              LoopFilterMask *lfm) {
  const MB_MODE_INFO *mbmi = &mi[0]->mbmi;
  block_size = (block_size == BLOCK_INVALID) ? mbmi->sb_type : block_size;

  const int lfm_bound = 1 << (MAX_SB_SIZE_LOG2 - 1);
  const int x = (idx >= lfm_bound) ? idx - lfm_bound : idx;
  const int y = (idy >= lfm_bound) ? idy - lfm_bound : idy;

  // We don't distinguish UV plane's edge direction
  uint8_t lfl =
      receive_filter_level(cm, &cm->lf_info, 0, plane, mi_row, mi_col, mbmi);

  const int wide = num_8x8_blocks_wide_lookup[block_size];
  const int high = num_8x8_blocks_high_lookup[block_size];
  int i;
  const int row = y >> (MI_SIZE_LOG2 + 1);
  const int col = x >> (MI_SIZE_LOG2 + 1);

  uint8_t *u_level =
      cm->lf.neighbor + 3 * (cm->lf.neighbor_width + cm->lf.neighbor_height);
  uint8_t *v_level = u_level + cm->lf.neighbor_width + cm->lf.neighbor_height;
  const int mi_row_pos = (mi_row + (idy >> MI_SIZE_LOG2)) >> 1;

  assert(plane == 1 || plane == 2);

  uint8_t lfl_f = 0;
  int lfl_flag = 0;
  // U plane
  if (plane == 1) {
    // if this flag = 0, we don't do loopfilter
    lfl_flag = u_level[mi_row_pos] || lfl;
    // if current level = 0, but previous level is nonzero, use previous level
    if (!lfl)
      lfl_f = u_level[mi_row_pos] ? u_level[mi_row_pos] : 0;
    else
      lfl_f = lfl;
    for (i = 0; i < high; i++) {
      memset(&lfm->lfl_u[row + i][col], lfl_f, wide);
    }
    memset(&u_level[mi_row_pos], lfl, high);
  }

  // V plane
  if (plane == 2) {
    lfl_flag = v_level[mi_row_pos] || lfl;
    if (!lfl)
      lfl_f = u_level[mi_row_pos] ? u_level[mi_row_pos] : 0;
    else
      lfl_f = lfl;
    for (i = 0; i < high; i++) {
      memset(&lfm->lfl_v[row + i][col], lfl_f, wide);
    }
    memset(&v_level[mi_row_pos], lfl, high);
  }
  return lfl_flag;
}

// Returned loopfilter flag is consistent across the whole coding block.
// We take the last one in returned by the last get_uv_level_64x64() call.
static int process_uv_level(AV1_COMMON *const cm, MODE_INFO **mi, int plane,
                            int idx, int idy, int mi_row, int mi_col,
                            LoopFilterMask *lfm) {
  assert(plane != 0);
  const MB_MODE_INFO *mbmi = &mi[0]->mbmi;
  const BLOCK_SIZE block_size = mbmi->sb_type;
  const int mi_stride = cm->mi_stride;
  const int block_64x64_step = 64;
  const int mi_64x64_step = block_64x64_step >> MI_SIZE_LOG2;
  MODE_INFO **mip = mi;
  int x = idx;
  int y = idy;
  int row = mi_row;
  int col = mi_col;
  LoopFilterMask *fm = lfm;
  int lfl_flag;

  switch (block_size) {
    case BLOCK_128X128:
      // top left
      get_uv_level_64x64(cm, mip, BLOCK_64X64, plane, x, y, row, col, fm);

      // top right
      x = idx + block_64x64_step;
      y = idy;
      col = mi_col + mi_64x64_step;
      row = mi_row;
      fm += 1;
      mip = mi + mi_64x64_step;
      get_uv_level_64x64(cm, mip, BLOCK_64X64, plane, x, y, row, col, fm);

      // bottom left
      x = idx;
      y = idy + block_64x64_step;
      col = mi_col;
      row = mi_row + mi_64x64_step;
      fm += 1;
      mip = mi + (mi_64x64_step * mi_stride);
      get_uv_level_64x64(cm, mip, BLOCK_64X64, plane, x, y, row, col, fm);

      // bottom right
      x = idx + block_64x64_step;
      y = idy + block_64x64_step;
      col = mi_col + mi_64x64_step;
      row = mi_row + mi_64x64_step;
      fm += 1;
      mip = mi + (mi_64x64_step * mi_stride) + mi_64x64_step;
      lfl_flag =
          get_uv_level_64x64(cm, mip, BLOCK_64X64, plane, x, y, row, col, fm);
      break;
    case BLOCK_128X64:
      // top left
      get_uv_level_64x64(cm, mip, BLOCK_64X64, plane, x, y, row, col, fm);

      // top right
      x = idx + block_64x64_step;
      y = idy;
      col = mi_col + mi_64x64_step;
      row = mi_row;
      fm += 1;
      mip = mi + mi_64x64_step;
      lfl_flag =
          get_uv_level_64x64(cm, mip, BLOCK_64X64, plane, x, y, row, col, fm);
      break;
    case BLOCK_64X128:
      // top left
      get_uv_level_64x64(cm, mip, BLOCK_64X64, plane, x, y, row, col, fm);

      // bottom left
      x = idx;
      y = idy + block_64x64_step;
      col = mi_col;
      row = mi_row + mi_64x64_step;
      fm += 2;
      mip = mi + (mi_64x64_step * mi_stride);
      lfl_flag =
          get_uv_level_64x64(cm, mip, BLOCK_64X64, plane, x, y, row, col, fm);
      break;
    default:
      fm += ((idy >> 6) << 1) + (idx >> 6);
      lfl_flag = get_uv_level_64x64(cm, mi, BLOCK_INVALID, plane, idx, idy,
                                    mi_row, mi_col, fm);
  }
  return lfl_flag;
}

static void setup_block_size(AV1_COMMON *const cm, MODE_INFO **mi,
                             int mi_stride, int block_wide, int block_high,
                             int idx, int idy, LoopFilterMask *lfm, int plane,
                             int mi_row, int mi_col, int br, int bc) {
  BLOCK_SIZE block_size = mi[0]->mbmi.sb_type;
  const int cur_block_wide = block_size_wide[block_size];
  const int cur_block_high = block_size_high[block_size];
  int sub_h, sub_w;
  int i, x, y, x_off, y_off;
  MODE_INFO **mi_next = NULL;

  if (block_wide == cur_block_wide && block_high == cur_block_high) {
    if (plane != 0) {
      const int lflf =
          process_uv_level(cm, mi, plane, idx, idy, mi_row, mi_col, lfm);
      if (lflf) build_uv_masks(cm, mi, plane, idx, idy, mi_row, mi_col, lfm);
    } else {
      int lflf_h = 0;
      int lflf_v = 0;
      process_y_level(cm, mi, idx, idy, mi_row, mi_col, lfm, &lflf_h, &lflf_v);
      build_y_masks(cm, mi, mi_stride, idx, idy, lfm, mi_row, mi_col, br, bc,
                    lflf_h, lflf_v);
    }
  } else if (block_wide == cur_block_wide) {
    sub_w = block_wide;
    if (cur_block_high == block_high >> 2) {
      // 4:1 rect
      sub_h = block_high >> 2;
      x_off = 0;
      for (i = 0; i < 4; ++i) {
        y_off = i * sub_h;

        if (y_off >= block_high) continue;

        x = idx + x_off;
        y = idy + y_off;
        if ((y >> 2) > br || (x >> 2) > bc) continue;

        mi_next =
            mi + (y_off >> MI_SIZE_LOG2) * mi_stride + (x_off >> MI_SIZE_LOG2);
        setup_block_size(cm, mi_next, mi_stride, sub_w, sub_h, x, y, lfm, plane,
                         mi_row, mi_col, br, bc);
      }
    } else {
      // 2:1 rect
      sub_h = block_high > cur_block_high ? block_high >> 1 : block_high << 1;

      for (i = 0; i < 4; ++i) {
        x_off = (i & 1) * sub_w;
        y_off = (i >> 1) * sub_h;

        if (x_off >= block_wide || y_off >= block_high) continue;

        x = idx + x_off;
        y = idy + y_off;
        if ((x >> 2) > bc || (y >> 2) > br) continue;

        mi_next =
            mi + (y_off >> MI_SIZE_LOG2) * mi_stride + (x_off >> MI_SIZE_LOG2);
        setup_block_size(cm, mi_next, mi_stride, sub_w, sub_h, x, y, lfm, plane,
                         mi_row, mi_col, br, bc);
      }
    }
  } else if (block_high == cur_block_high) {
    sub_h = block_high;
    if (cur_block_wide == block_wide >> 2) {
      // 1:4 rect
      sub_w = block_wide >> 2;
      y_off = 0;
      for (i = 0; i < 4; ++i) {
        x_off = i * sub_w;

        if (x_off >= block_wide) continue;

        x = idx + x_off;
        y = idy + y_off;
        if ((y >> 2) > br || (x >> 2) > bc) continue;

        mi_next =
            mi + (y_off >> MI_SIZE_LOG2) * mi_stride + (x_off >> MI_SIZE_LOG2);
        setup_block_size(cm, mi_next, mi_stride, sub_w, sub_h, x, y, lfm, plane,
                         mi_row, mi_col, br, bc);
      }
    } else {
      // 1:2 rect
      sub_w = block_wide > cur_block_wide ? block_wide >> 1 : block_wide << 1;
      for (i = 0; i < 4; ++i) {
        x_off = (i & 1) * sub_w;
        y_off = (i >> 1) * sub_h;

        if (x_off >= block_wide || y_off >= block_high) break;

        x = idx + x_off;
        y = idy + y_off;
        if ((x >> 2) > bc || (y >> 2) > br) continue;

        mi_next =
            mi + (y_off >> MI_SIZE_LOG2) * mi_stride + (x_off >> MI_SIZE_LOG2);
        setup_block_size(cm, mi_next, mi_stride, sub_w, sub_h, x, y, lfm, plane,
                         mi_row, mi_col, br, bc);
      }
    }
  } else {
    sub_w = block_wide >> 1;
    sub_h = block_high >> 1;

    int count = 4;
    for (i = 0; i < count; ++i) {
      x_off = (i & 1) * sub_w;
      y_off = (i >> 1) * sub_h;

      y = idy + y_off;
      x = idx + x_off;
      if ((x >> 2) > bc || (y >> 2) > br) continue;

      mi_next =
          mi + (y_off >> MI_SIZE_LOG2) * mi_stride + (x_off >> MI_SIZE_LOG2);

      block_size = mi_next[0]->mbmi.sb_type;
      int w = block_size_wide[block_size];
      int h = block_size_high[block_size];

      // PARTITION_HORZ_A
      // After finishing 2NxN rectangle partition, we will finish the loop
      // and leave current level without doing i=3.
      if (h == sub_h && w == sub_w << 1) count -= 2;

      // Note:
      // PARTITION_VERT_A has a special situation; we need to do NxN top left
      // and bottom left square partition first. Then we do Nx2N lastly. This
      // is to make sure left neighbor tx_size is always available.
      // For PARTITION_VERT_A, we skip Nx2N rectangle and do i=2 square first
      if (w == sub_w && h == sub_h << 1 && i == 1) {
        i += 1;
        // when i=2, x_off = 0 * sub_w, y_off = 1 * sub_h
        y = idy + sub_h;
        x = idx;
        if ((x >> 2) > bc || (y >> 2) > br) continue;
        mi_next = mi + (sub_h >> MI_SIZE_LOG2) * mi_stride;
      }

      // For PARTITION_VERT_A, we then do Nx2N rectangle, i=1
      if (w == sub_w && h == sub_h << 1 && i == 3) {
        // when i=1, x_off = 1 * sub_w, y_off = 0 * sub_h
        y = idy;
        x = idx + sub_w;
        mi_next = mi + (sub_w >> MI_SIZE_LOG2);
      }

      setup_block_size(cm, mi_next, mi_stride, sub_w, sub_h, x, y, lfm, plane,
                       mi_row, mi_col, br, bc);
    }
  }
}

#if CONFIG_PARALLEL_DEBLOCKING
// Note:
// This loopfilter implementation builds bit mask from var_tx and cb4x4.
// Loop filter uses cb4x4 filter which is the same as parallel_deblocking.
#define USE_BITMASK_FROM_VAR_TX 1
#endif

static void get_y_border_left_mask(FilterMaskY *m, int rows) {
  int index = 0;
  const int shift_y = get_y_index_shift(0, rows << MI_SIZE_LOG2, &index);
  int i;
  for (i = 0; i < index; ++i) {
    m->bits[i] = 0xffffffffffffffffULL;
  }
  m->bits[index] = ((uint64_t)1 << shift_y) - 1;
}

static void get_uv_border_left_mask(FilterMaskUV *m, int rows) {
  const int shift_uv = get_uv_shift(0, rows >> 1);
  if (shift_uv + 8 == 64)
    *m = 0xffffffffffffffffULL;
  else
    *m = ((uint64_t)1 << (shift_uv + 8)) - 1;
}

static void set_border_mask(const FilterMaskY *mask_y,
                            const FilterMaskUV *mask_uv, LoopFilterMask *lfm) {
  TX_SIZE max_tx_size = TX_32X32;

  // Remove values completely outside our border.
  int i;
  for (i = 0; i < max_tx_size; ++i) {
    int j;
    for (j = 0; j < 4; ++j) {
      lfm->left_y[i].bits[j] &= mask_y->bits[j];
      lfm->above_y[i].bits[j] &= mask_y->bits[j];
    }
    lfm->left_u[i] &= *mask_uv;
    lfm->above_u[i] &= *mask_uv;
    lfm->left_v[i] &= *mask_uv;
    lfm->above_v[i] &= *mask_uv;
  }
}

static void get_y_border_top_mask(FilterMaskY *m, int cols) {
  const uint64_t mask = (((1 << cols) - 1)) * 0x0001000100010001ULL;
  int i;
  for (i = 0; i < 4; ++i) {
    m->bits[i] = mask;
  }
}

static void get_uv_border_top_mask(FilterMaskUV *m, int cols) {
  *m = (FilterMaskUV)(((1 << cols) - 1)) * 0x0101010101010101ULL;
}

static LpfMask *get_loop_filter_mask(AV1_COMMON *const cm, int mi_row,
                                     int mi_col) {
  assert(cm->lf.lfm != NULL);
  const int sb_row = mi_row >> MAX_MIB_SIZE_LOG2;
  const int sb_col = mi_col >> MAX_MIB_SIZE_LOG2;
  const int lfm_stride = cm->lf.lfm_stride;
  LpfMask *lfmask = cm->lf.lfm;
  lfmask += sb_row * lfm_stride + sb_col;

  return lfmask;
}

// This function sets up the bit masks for the entire 128x128 region represented
// by mi_row, mi_col.
// TODO(luoyi): This function only works for pixel sampling format 4:2:0.
void av1_setup_mask(AV1_COMMON *const cm, const int mi_row, const int mi_col,
                    MODE_INFO **mi, int plane, LpfMask *lpfmask) {
  assert(mi[0] != NULL);

  // TODO(luoyi): Filter level may change for each loopfilter operation.
  // But the bit mask based on tx_size will not change. We'd better
  // seperate these two parts so that bit mask will not be built each time.
  // if (lpfmask->is_setup) return;

  LoopFilterMask *lfm = &lpfmask->lfm[0];
  memset(lfm, 0, LOOP_FILTER_MASK_NUM * sizeof(*lfm));

  BLOCK_SIZE block_size = cm->sb_size;
  const int block_wide = block_size_wide[block_size];
  const int block_high = block_size_high[block_size];
  const int row_bound = cm->mi_rows - 1 - mi_row;
  const int col_bound = cm->mi_cols - 1 - mi_col;
  const int mode_info_stride = cm->mi_stride;
  setup_block_size(cm, mi, mode_info_stride, block_wide, block_high, 0, 0, lfm,
                   plane, mi_row, mi_col, row_bound, col_bound);

  int lfm_ind, i;
  for (lfm_ind = 0; lfm_ind < LOOP_FILTER_MASK_NUM; ++lfm_ind) {
    // The largest loopfilter we have is 16x16 so we use the 16x16 mask
    // for 32x32 and 64x64 transforms also.
    for (i = 0; i < 4; ++i) {
      lfm->left_y[TX_16X16].bits[i] |= lfm->left_y[TX_32X32].bits[i];
      lfm->above_y[TX_16X16].bits[i] |= lfm->above_y[TX_32X32].bits[i];
#if CONFIG_TX64X64
      lfm->left_y[TX_16X16].bits[i] |= lfm->left_y[TX_64X64].bits[i];
      lfm->above_y[TX_16X16].bits[i] |= lfm->above_y[TX_64X64].bits[i];
#endif
    }
    lfm->left_u[TX_16X16] |= lfm->left_u[TX_32X32];
    lfm->above_u[TX_16X16] |= lfm->above_u[TX_32X32];

    lfm->left_v[TX_16X16] |= lfm->left_v[TX_32X32];
    lfm->above_v[TX_16X16] |= lfm->above_v[TX_32X32];
// TODO(luoyi): These are from VP9 loopfilter. AV1 loopfilter no longer uses
// it. We need some experiment to see if it could improve performance.
#if 0
    // We do at least 8 tap filter on every 32x32 even if the transform size
    // is 4x4. So if the 4x4 is set on a border pixel add it to the 8x8 and
    // remove it from the 4x4.
    for (i = 0; i < 4; ++i) {
      lfm->left_y[TX_8X8].bits[i] |=
          lfm->left_y[TX_4X4].bits[i] & left_border.bits[i];

      lfm->left_y[TX_4X4].bits[i] &= ~left_border.bits[i];

      lfm->above_y[TX_8X8].bits[i] |=
          lfm->above_y[TX_4X4].bits[i] & above_border.bits[i];

      lfm->above_y[TX_4X4].bits[i] &= ~above_border.bits[i];
    }
    lfm->left_u[TX_8X8] |= lfm->left_u[TX_4X4] & left_border_uv;
    lfm->left_u[TX_4X4] &= ~left_border_uv;
    lfm->above_u[TX_8X8] |= lfm->above_u[TX_4X4] & above_border_uv;
    lfm->above_u[TX_4X4] &= ~above_border_uv;

    lfm->left_v[TX_8X8] |= lfm->left_v[TX_4X4] & left_border_uv;
    lfm->left_v[TX_4X4] &= ~left_border_uv;
    lfm->above_v[TX_8X8] |= lfm->above_v[TX_4X4] & above_border_uv;
    lfm->above_v[TX_4X4] &= ~above_border_uv;
#endif
    lfm += 1;
  }

  lfm = &lpfmask->lfm[0];

  // We do some special edge handling.
  if (mi_row + MAX_MIB_SIZE > cm->mi_rows) {
    FilterMaskY mask_y = { { 0, 0, 0, 0 } };
    FilterMaskUV mask_uv = 0;
    int rows;
    if (mi_row + (MAX_MIB_SIZE >> 1) == cm->mi_rows) {
      memset(&lfm[2], 0, 2 * sizeof(*lfm));
    } else if (mi_row + (MAX_MIB_SIZE >> 1) > cm->mi_rows) {
      memset(&lfm[2], 0, 2 * sizeof(*lfm));

      // Each pixel inside the border gets a 1,
      rows = cm->mi_rows - mi_row;
      get_y_border_left_mask(&mask_y, rows);
      get_uv_border_left_mask(&mask_uv, rows);

      set_border_mask(&mask_y, &mask_uv, &lfm[0]);
      set_border_mask(&mask_y, &mask_uv, &lfm[1]);
    } else {
      rows = cm->mi_rows - mi_row;
      rows -= 16;

      get_y_border_left_mask(&mask_y, rows);
      get_uv_border_left_mask(&mask_uv, rows);

      set_border_mask(&mask_y, &mask_uv, &lfm[2]);
      set_border_mask(&mask_y, &mask_uv, &lfm[3]);
    }

// TODO(luoyi): These are from VP9 loopfilter. AV1 loopfilter no longer uses
// it. We need some experiment to see if it could improve performance.
// Also these have not converted to AV1 yet.
#if 0
    // We don't apply a wide loop filter on the last uv block row. If set
    // apply the shorter one instead.
    if (rows == 1) {
      lfm->above_uv[TX_8X8] |= lfm->above_uv[TX_16X16];
      lfm->above_uv[TX_16X16] = 0;
    }
    if (rows == 9) {
      const FilterMaskUV half_mask = 0xffffffff00000000;
      lfm->above_uv[TX_8X8] |= lfm->above_uv[TX_16X16] & half_mask;
      lfm->above_uv[TX_16X16] &= ~(lfm->above_uv[TX_16X16] & half_mask);
    }
#endif
  }

  if (mi_col + MAX_MIB_SIZE > cm->mi_cols) {
    FilterMaskY mask_y = { { 0, 0, 0, 0 } };
    FilterMaskUV mask_uv = 0;
    int cols;
    // Each pixel inside the border gets a 1, the multiply copies the border
    // to where we need it.
    if (mi_col + (MAX_MIB_SIZE >> 1) == cm->mi_cols) {
      memset(&lfm[1], 0, sizeof(*lfm));
      memset(&lfm[3], 0, sizeof(*lfm));
    } else if (mi_col + (MAX_MIB_SIZE >> 1) > cm->mi_cols) {
      memset(&lfm[1], 0, sizeof(*lfm));
      memset(&lfm[3], 0, sizeof(*lfm));

      cols = cm->mi_cols - mi_col;

      get_y_border_top_mask(&mask_y, cols);
      get_uv_border_top_mask(&mask_uv, cols);

      set_border_mask(&mask_y, &mask_uv, &lfm[0]);
      set_border_mask(&mask_y, &mask_uv, &lfm[2]);
    } else {
      cols = cm->mi_cols - mi_col;
      cols -= 16;

      get_y_border_top_mask(&mask_y, cols);
      get_uv_border_top_mask(&mask_uv, cols);

      set_border_mask(&mask_y, &mask_uv, &lfm[1]);
      set_border_mask(&mask_y, &mask_uv, &lfm[3]);
    }

// TODO(luoyi): These are from VP9 loopfilter. AV1 loopfilter no longer uses
// it. We need some experiment to see if it could improve performance.
// Also these have not converted to AV1 yet.
#if 0
    // We don't apply a wide loop filter on the last uv column. If set
    // apply the shorter one instead.
    if (columns == 1) {
      lfm->left_uv[TX_8X8] |= lfm->left_uv[TX_16X16];
      lfm->left_uv[TX_16X16] = 0;
    }
    if (columns == 9) {
      const uint64_t mask = 0x0f0f0f0f0f0f0f0fULL;
      lfm->left_uv[TX_8X8] |= (lfm->left_uv[TX_16X16] & mask);
      lfm->left_uv[TX_16X16] &= ~(lfm->left_uv[TX_16X16] & mask);
    }
#endif
  }

  // We don't apply a loop filter on the first column in the image, mask that
  // out.
  if (mi_col == 0) {
    int count = 2;
    while (count > 0) {
      for (i = 0; i < TX_32X32; i++) {
        int j;
        for (j = 0; j < 4; ++j) {
          lfm->left_y[i].bits[j] &= 0xfffefffefffefffeULL;
        }
        lfm->left_u[i] &= 0xfefefefefefefefeULL;
        lfm->left_v[i] &= 0xfefefefefefefefeULL;
      }
      lfm += 2;
      count -= 1;
    }
  }

  lfm = &lpfmask->lfm[0];
#if CONFIG_LOOPFILTERING_ACROSS_TILES || CONFIG_LOOPFILTERING_ACROSS_TILES_EXT
  if (av1_disable_loopfilter_on_tile_boundary(cm)) {
    update_tile_boundary_filter_mask(cm, mi_row, mi_col, lfm);
  }
#endif  // CONFIG_LOOPFILTERING_ACROSS_TILES

  check_loop_filter_masks(lfm);

  // TODO(luoyi): Filter level may change for each loopfilter operation.
  // But the bit mask based on tx_size will not change. We'd better
  // seperate these two parts so that bit mask will not be built each time.
  // lpfmask->is_setup = 1;
}

static void filter_selectively_vert(
    uint8_t *s, int pitch, unsigned int mask_16x16, unsigned int mask_8x8,
    unsigned int mask_4x4, unsigned int mask_4x4_int,
    const loop_filter_info_n *lfi_n, const uint8_t *lfl) {
  unsigned int mask;

  for (mask = mask_16x16 | mask_8x8 | mask_4x4 | mask_4x4_int; mask;
       mask >>= 1) {
    const loop_filter_thresh *lfi = lfi_n->lfthr + *lfl;

    if (mask & 1) {
      if (mask_16x16 & 1) {
        aom_lpf_vertical_16(s, pitch, lfi->mblim, lfi->lim, lfi->hev_thr);
      } else if (mask_8x8 & 1) {
        aom_lpf_vertical_8(s, pitch, lfi->mblim, lfi->lim, lfi->hev_thr);
      } else if (mask_4x4 & 1) {
        aom_lpf_vertical_4(s, pitch, lfi->mblim, lfi->lim, lfi->hev_thr);
      }
    }
    if (mask_4x4_int & 1)
      aom_lpf_vertical_4(s + 4, pitch, lfi->mblim, lfi->lim, lfi->hev_thr);
    s += 8;
    lfl += 1;
    mask_16x16 >>= 1;
    mask_8x8 >>= 1;
    mask_4x4 >>= 1;
    mask_4x4_int >>= 1;
  }
}

static void highbd_filter_selectively_vert(
    uint16_t *s, int pitch, unsigned int mask_16x16, unsigned int mask_8x8,
    unsigned int mask_4x4, unsigned int mask_4x4_int,
    const loop_filter_info_n *lfi_n, const uint8_t *lfl, int bd) {
  unsigned int mask;

  for (mask = mask_16x16 | mask_8x8 | mask_4x4 | mask_4x4_int; mask;
       mask >>= 1) {
    const loop_filter_thresh *lfi = lfi_n->lfthr + *lfl;

    if (mask & 1) {
      if (mask_16x16 & 1) {
        aom_highbd_lpf_vertical_16(s, pitch, lfi->mblim, lfi->lim, lfi->hev_thr,
                                   bd);
      } else if (mask_8x8 & 1) {
        aom_highbd_lpf_vertical_8(s, pitch, lfi->mblim, lfi->lim, lfi->hev_thr,
                                  bd);
      } else if (mask_4x4 & 1) {
        aom_highbd_lpf_vertical_4(s, pitch, lfi->mblim, lfi->lim, lfi->hev_thr,
                                  bd);
      }
    }
    if (mask_4x4_int & 1)
      aom_highbd_lpf_vertical_4(s + 4, pitch, lfi->mblim, lfi->lim,
                                lfi->hev_thr, bd);
    s += 8;
    lfl += 1;
    mask_16x16 >>= 1;
    mask_8x8 >>= 1;
    mask_4x4 >>= 1;
    mask_4x4_int >>= 1;
  }
}

typedef struct {
  unsigned int m16x16;
  unsigned int m8x8;
  unsigned int m4x4;
} FilterMasks;

// Get filter level and masks for the given row index 'idx_r'. (Only used for
// the non420 case).
// Note: 'row_masks_ptr' and/or 'col_masks_ptr' can be passed NULL.
static void get_filter_level_and_masks_non420(
    AV1_COMMON *const cm, const struct macroblockd_plane *const plane, int pl,
    MODE_INFO **mib, int mi_row, int mi_col, int idx_r, uint8_t *const lfl_r,
    unsigned int *const mask_4x4_int_r_ptr,
    unsigned int *const mask_4x4_int_c_ptr, FilterMasks *const row_masks_ptr,
    FilterMasks *const col_masks_ptr) {
  const int ss_x = plane->subsampling_x;
  const int ss_y = plane->subsampling_y;
  const int col_step = mi_size_wide[BLOCK_8X8] << ss_x;
  FilterMasks row_masks, col_masks;
  memset(&row_masks, 0, sizeof(row_masks));
  memset(&col_masks, 0, sizeof(col_masks));
  unsigned int mask_4x4_int_r = 0, mask_4x4_int_c = 0;
  const int r = idx_r >> mi_height_log2_lookup[BLOCK_8X8];

  // Determine the vertical edges that need filtering
  int idx_c;
  for (idx_c = 0; idx_c < cm->mib_size && mi_col + idx_c < cm->mi_cols;
       idx_c += col_step) {
    const MODE_INFO *mi = mib[idx_r * cm->mi_stride + idx_c];
    const MB_MODE_INFO *mbmi = &mi[0].mbmi;
    const BLOCK_SIZE sb_type = mbmi->sb_type;
    const int skip_this = mbmi->skip && is_inter_block(mbmi);
    // Map index to 8x8 unit
    const int c = idx_c >> mi_width_log2_lookup[BLOCK_8X8];

    const int blk_row = r & (num_8x8_blocks_high_lookup[sb_type] - 1);
    const int blk_col = c & (num_8x8_blocks_wide_lookup[sb_type] - 1);

    // left edge of current unit is block/partition edge -> no skip
    const int block_edge_left =
        (num_4x4_blocks_wide_lookup[sb_type] > 1) ? !blk_col : 1;
    const int skip_this_c = skip_this && !block_edge_left;
    // top edge of current unit is block/partition edge -> no skip
    const int block_edge_above =
        (num_4x4_blocks_high_lookup[sb_type] > 1) ? !blk_row : 1;
    const int skip_this_r = skip_this && !block_edge_above;

    TX_SIZE tx_size = (plane->plane_type == PLANE_TYPE_UV)
                          ? av1_get_uv_tx_size(mbmi, ss_x, ss_y)
                          : mbmi->tx_size;

    const int skip_border_4x4_c =
        ss_x && mi_col + idx_c >= cm->mi_cols - mi_size_wide[BLOCK_8X8];
    const int skip_border_4x4_r =
        ss_y && mi_row + idx_r >= cm->mi_rows - mi_size_high[BLOCK_8X8];

    int tx_size_mask = 0;
    const int c_step = (c >> ss_x);
    const int r_step = (r >> ss_y);
    const int col_mask = 1 << c_step;

    if (is_inter_block(mbmi) && !mbmi->skip) {
      const int tx_row_idx =
          (blk_row * mi_size_high[BLOCK_8X8] << TX_UNIT_HIGH_LOG2) >> 1;
      const int tx_col_idx =
          (blk_col * mi_size_wide[BLOCK_8X8] << TX_UNIT_WIDE_LOG2) >> 1;
      const TX_SIZE mb_tx_size = mbmi->inter_tx_size[tx_row_idx][tx_col_idx];
      tx_size = (plane->plane_type == PLANE_TYPE_UV)
                    ? av1_get_uv_tx_size(mbmi, ss_x, ss_y)
                    : mb_tx_size;
    }

// Filter level can vary per MI
#if CONFIG_EXT_DELTA_Q
#if CONFIG_LOOPFILTER_LEVEL
    if (!(lfl_r[c_step] = get_filter_level(cm, &cm->lf_info, 0, 0, mbmi)))
      continue;
#else
    if (!(lfl_r[c_step] = get_filter_level(cm, &cm->lf_info, mbmi))) continue;
#endif
#else
    if (!(lfl_r[c_step] = get_filter_level(&cm->lf_info, mbmi))) continue;
#endif

    TX_SIZE tx_size_horz_edge, tx_size_vert_edge;

    // filt_len_vert_edge is the length of deblocking filter for a vertical edge
    // The filter direction of a vertical edge is horizontal.
    // Thus, filt_len_vert_edge is determined as the minimum width of the two
    // transform block sizes on the left and right (current block) side of edge
    const int filt_len_vert_edge = AOMMIN(
        tx_size_wide[tx_size],
        tx_size_wide[cm->left_txfm_context[pl][((mi_row + idx_r) & MAX_MIB_MASK)
                                               << TX_UNIT_HIGH_LOG2]]);

    // filt_len_horz_edge is the len of deblocking filter for a horizontal edge
    // The filter direction of a horizontal edge is vertical.
    // Thus, filt_len_horz_edge is determined as the minimum height of the two
    // transform block sizes on the top and bottom (current block) side of edge
    const int filt_len_horz_edge =
        AOMMIN(tx_size_high[tx_size],
               tx_size_high[cm->top_txfm_context[pl][(mi_col + idx_c)
                                                     << TX_UNIT_WIDE_LOG2]]);

    // transform width/height of current block
    const int tx_wide_cur = tx_size_wide[tx_size];
    const int tx_high_cur = tx_size_high[tx_size];

    // tx_size_vert_edge is square transform size for a vertical deblocking edge
    // It determines the type of filter applied to the vertical edge
    // Similarly, tx_size_horz_edge is for a horizontal deblocking edge
    tx_size_vert_edge = get_sqr_tx_size(filt_len_vert_edge);
    tx_size_horz_edge = get_sqr_tx_size(filt_len_horz_edge);

    memset(cm->top_txfm_context[pl] + ((mi_col + idx_c) << TX_UNIT_WIDE_LOG2),
           tx_size, mi_size_wide[BLOCK_8X8] << TX_UNIT_WIDE_LOG2);
    memset(cm->left_txfm_context[pl] +
               (((mi_row + idx_r) & MAX_MIB_MASK) << TX_UNIT_HIGH_LOG2),
           tx_size, mi_size_high[BLOCK_8X8] << TX_UNIT_HIGH_LOG2);

    if (tx_size_vert_edge == TX_32X32)
      tx_size_mask = 3;
    else if (tx_size_vert_edge == TX_16X16)
      tx_size_mask = 1;
    else
      tx_size_mask = 0;

    // Build masks based on the transform size of each block
    // handle vertical mask
    if (tx_size_vert_edge == TX_32X32) {
      if (!skip_this_c && (c_step & tx_size_mask) == 0) {
        if (!skip_border_4x4_c)
          col_masks.m16x16 |= col_mask;
        else
          col_masks.m8x8 |= col_mask;
      }
    } else if (tx_size_vert_edge == TX_16X16) {
      if (!skip_this_c && (c_step & tx_size_mask) == 0) {
        if (!skip_border_4x4_c)
          col_masks.m16x16 |= col_mask;
        else
          col_masks.m8x8 |= col_mask;
      }
    } else {
      // force 8x8 filtering on 32x32 boundaries
      if (!skip_this_c && (c_step & tx_size_mask) == 0) {
        if (tx_size_vert_edge == TX_8X8 || (c_step & 3) == 0)
          col_masks.m8x8 |= col_mask;
        else
          col_masks.m4x4 |= col_mask;
      }

      if (!skip_this && tx_wide_cur < 8 && !skip_border_4x4_c &&
          (c_step & tx_size_mask) == 0)
        mask_4x4_int_c |= col_mask;
    }

    if (tx_size_horz_edge == TX_32X32)
      tx_size_mask = 3;
    else if (tx_size_horz_edge == TX_16X16)
      tx_size_mask = 1;
    else
      tx_size_mask = 0;

    // set horizontal mask
    if (tx_size_horz_edge == TX_32X32) {
      if (!skip_this_r && (r_step & tx_size_mask) == 0) {
        if (!skip_border_4x4_r)
          row_masks.m16x16 |= col_mask;
        else
          row_masks.m8x8 |= col_mask;
      }
    } else if (tx_size_horz_edge == TX_16X16) {
      if (!skip_this_r && (r_step & tx_size_mask) == 0) {
        if (!skip_border_4x4_r)
          row_masks.m16x16 |= col_mask;
        else
          row_masks.m8x8 |= col_mask;
      }
    } else {
      // force 8x8 filtering on 32x32 boundaries
      if (!skip_this_r && (r_step & tx_size_mask) == 0) {
        if (tx_size_horz_edge == TX_8X8 || (r_step & 3) == 0)
          row_masks.m8x8 |= col_mask;
        else
          row_masks.m4x4 |= col_mask;
      }

      if (!skip_this && tx_high_cur < 8 && !skip_border_4x4_r &&
          (r_step & tx_size_mask) == 0)
        mask_4x4_int_r |= col_mask;
    }
  }

  if (row_masks_ptr) *row_masks_ptr = row_masks;
  if (col_masks_ptr) *col_masks_ptr = col_masks;
  if (mask_4x4_int_c_ptr) *mask_4x4_int_c_ptr = mask_4x4_int_c;
  if (mask_4x4_int_r_ptr) *mask_4x4_int_r_ptr = mask_4x4_int_r;
}

void av1_filter_block_plane_non420_ver(AV1_COMMON *const cm,
                                       struct macroblockd_plane *plane,
                                       MODE_INFO **mib, int mi_row, int mi_col,
                                       int pl) {
  const int ss_y = plane->subsampling_y;
  const int row_step = mi_size_high[BLOCK_8X8] << ss_y;
  struct buf_2d *const dst = &plane->dst;
  uint8_t *const dst0 = dst->buf;
  uint8_t lfl[MAX_MIB_SIZE][MAX_MIB_SIZE] = { { 0 } };

  int idx_r;
  for (idx_r = 0; idx_r < cm->mib_size && mi_row + idx_r < cm->mi_rows;
       idx_r += row_step) {
    unsigned int mask_4x4_int;
    FilterMasks col_masks;
    const int r = idx_r >> mi_height_log2_lookup[BLOCK_8X8];
    get_filter_level_and_masks_non420(cm, plane, pl, mib, mi_row, mi_col, idx_r,
                                      &lfl[r][0], NULL, &mask_4x4_int, NULL,
                                      &col_masks);

    // Disable filtering on the leftmost column or tile boundary
    unsigned int border_mask = ~(mi_col == 0 ? 1 : 0);
#if CONFIG_LOOPFILTERING_ACROSS_TILES || CONFIG_LOOPFILTERING_ACROSS_TILES_EXT
    const BOUNDARY_TYPE *const bi =
        cm->boundary_info + (mi_row + idx_r) * cm->mi_stride + mi_col;
    if (av1_disable_loopfilter_on_tile_boundary(cm) &&
        ((*bi & TILE_LEFT_BOUNDARY) != 0)) {
      border_mask = 0xfffffffe;
    }
#endif  // CONFIG_LOOPFILTERING_ACROSS_TILES

    if (cm->use_highbitdepth)
      highbd_filter_selectively_vert(
          CONVERT_TO_SHORTPTR(dst->buf), dst->stride,
          col_masks.m16x16 & border_mask, col_masks.m8x8 & border_mask,
          col_masks.m4x4 & border_mask, mask_4x4_int, &cm->lf_info, &lfl[r][0],
          (int)cm->bit_depth);
    else
      filter_selectively_vert(
          dst->buf, dst->stride, col_masks.m16x16 & border_mask,
          col_masks.m8x8 & border_mask, col_masks.m4x4 & border_mask,
          mask_4x4_int, &cm->lf_info, &lfl[r][0]);
    dst->buf += 8 * dst->stride;
  }

  // Now do horizontal pass
  dst->buf = dst0;
}

void av1_filter_block_plane_non420_hor(AV1_COMMON *const cm,
                                       struct macroblockd_plane *plane,
                                       MODE_INFO **mib, int mi_row, int mi_col,
                                       int pl) {
  const int ss_y = plane->subsampling_y;
  const int row_step = mi_size_high[BLOCK_8X8] << ss_y;
  struct buf_2d *const dst = &plane->dst;
  uint8_t *const dst0 = dst->buf;
  uint8_t lfl[MAX_MIB_SIZE][MAX_MIB_SIZE] = { { 0 } };

  int idx_r;
  for (idx_r = 0; idx_r < cm->mib_size && mi_row + idx_r < cm->mi_rows;
       idx_r += row_step) {
    unsigned int mask_4x4_int;
    FilterMasks row_masks;
    const int r = idx_r >> mi_height_log2_lookup[BLOCK_8X8];
    get_filter_level_and_masks_non420(cm, plane, pl, mib, mi_row, mi_col, idx_r,
                                      &lfl[r][0], &mask_4x4_int, NULL,
                                      &row_masks, NULL);

#if CONFIG_LOOPFILTERING_ACROSS_TILES || CONFIG_LOOPFILTERING_ACROSS_TILES_EXT
    // Disable filtering on the abovemost row or tile boundary
    const BOUNDARY_TYPE *bi =
        cm->boundary_info + (mi_row + idx_r) * cm->mi_stride + mi_col;
    if ((av1_disable_loopfilter_on_tile_boundary(cm) &&
         (*bi & TILE_ABOVE_BOUNDARY)) ||
        (mi_row + idx_r == 0))
      memset(&row_masks, 0, sizeof(row_masks));
#else
    if (mi_row + idx_r == 0) memset(&row_masks, 0, sizeof(row_masks));
#endif  // CONFIG_LOOPFILTERING_ACROSS_TILES

    if (cm->use_highbitdepth)
      highbd_filter_selectively_horiz(
          CONVERT_TO_SHORTPTR(dst->buf), dst->stride, row_masks.m16x16,
          row_masks.m8x8, row_masks.m4x4, mask_4x4_int, &cm->lf_info,
          &lfl[r][0], (int)cm->bit_depth);
    else
      filter_selectively_horiz(dst->buf, dst->stride, row_masks.m16x16,
                               row_masks.m8x8, row_masks.m4x4, mask_4x4_int,
                               &cm->lf_info, &lfl[r][0]);
    dst->buf += 8 * dst->stride;
  }
  dst->buf = dst0;
}

void av1_filter_block_plane_ss00_ver(AV1_COMMON *const cm,
                                     struct macroblockd_plane *const plane,
                                     int mi_row, LoopFilterMask *lfm) {
  assert(plane->subsampling_x == 0 && plane->subsampling_y == 0);

  struct buf_2d *const dst = &plane->dst;
  uint8_t *buf = dst->buf;
  const int bitmask = 0xffffffff;
  int r;
#if CONFIG_EXT_PARTITION
  int sb_count = 1;
  if (cm->sb_size == BLOCK_128X128) sb_count = 2;
#endif

  const LoopFilterMask *mask = &lfm[0];
  const FilterMaskY *mask_16x16 = &mask->left_y[TX_16X16];
  const FilterMaskY *mask_8x8 = &mask->left_y[TX_8X8];
  const FilterMaskY *mask_4x4 = &mask->left_y[TX_4X4];
  uint32_t mask_16x16_l, mask_8x8_l, mask_4x4_l;

#if CONFIG_EXT_PARTITION
  uint8_t *buf2 = buf + (MAX_SB_SIZE >> 1);
  const LoopFilterMask *mask2 = &lfm[1];
  const FilterMaskY *mask2_16x16 = &mask2->left_y[TX_16X16];
  const FilterMaskY *mask2_8x8 = &mask2->left_y[TX_8X8];
  const FilterMaskY *mask2_4x4 = &mask2->left_y[TX_4X4];
  uint32_t mask2_16x16_l, mask2_8x8_l, mask2_4x4_l;
#endif

  const int mod = MAX_MIB_SIZE >> 1;
  const int yuv_comp = 0;

  // Vertical pass: do 2 rows at one time
  for (r = 0; r < cm->mib_size && mi_row + r < cm->mi_rows; r += 2) {
    const int index = (r % mod) >> 2;
    const int shift = (r % 4) << 4;

    mask_16x16_l = (uint32_t)(mask_16x16->bits[index] >> shift) & bitmask;
    mask_8x8_l = (uint32_t)(mask_8x8->bits[index] >> shift) & bitmask;
    mask_4x4_l = (uint32_t)(mask_4x4->bits[index] >> shift) & bitmask;
#if CONFIG_EXT_PARTITION
    if (sb_count > 1) {
      mask2_16x16_l = (uint32_t)(mask2_16x16->bits[index] >> shift) & bitmask;
      mask2_8x8_l = (uint32_t)(mask2_8x8->bits[index] >> shift) & bitmask;
      mask2_4x4_l = (uint32_t)(mask2_4x4->bits[index] >> shift) & bitmask;
    }
#endif

    // Disable filtering on the leftmost column.
    if (cm->use_highbitdepth)
      highbd_filter_selectively_vert_row2(
          plane->subsampling_x, CONVERT_TO_SHORTPTR(buf), yuv_comp, dst->stride,
          mask_16x16_l, mask_8x8_l, mask_4x4_l, &cm->lf_info,
          &mask->lfl_y_ver[r % mod][0], (int)cm->bit_depth);
    else
      filter_selectively_vert_row2(
          plane->subsampling_x, buf, yuv_comp, dst->stride, mask_16x16_l,
          mask_8x8_l, mask_4x4_l, &cm->lf_info, &mask->lfl_y_ver[r % mod][0]);

#if CONFIG_EXT_PARTITION
    if (sb_count > 1) {
      if (cm->use_highbitdepth)
        highbd_filter_selectively_vert_row2(
            plane->subsampling_x, CONVERT_TO_SHORTPTR(buf2), yuv_comp,
            dst->stride, mask2_16x16_l, mask2_8x8_l, mask2_4x4_l, &cm->lf_info,
            &mask2->lfl_y_ver[r % mod][0], (int)cm->bit_depth);
      else
        filter_selectively_vert_row2(plane->subsampling_x, buf2, yuv_comp,
                                     dst->stride, mask2_16x16_l, mask2_8x8_l,
                                     mask2_4x4_l, &cm->lf_info,
                                     &mask2->lfl_y_ver[r % mod][0]);
    }
#endif

    buf += 2 * MI_SIZE * dst->stride;
#if CONFIG_EXT_PARTITION
    buf2 += 2 * MI_SIZE * dst->stride;

    if (sb_count > 1 && r + 2 == MI_SIZE_64X64) {
      mask += 2;
      mask2 += 2;

      mask_16x16 = &mask->left_y[TX_16X16];
      mask_8x8 = &mask->left_y[TX_8X8];
      mask_4x4 = &mask->left_y[TX_4X4];

      mask2_16x16 = &mask2->left_y[TX_16X16];
      mask2_8x8 = &mask2->left_y[TX_8X8];
      mask2_4x4 = &mask2->left_y[TX_4X4];
    }
#endif
  }
}

void av1_filter_block_plane_ss00_hor(AV1_COMMON *const cm,
                                     struct macroblockd_plane *const plane,
                                     int mi_row, LoopFilterMask *lfm) {
  assert(plane->subsampling_x == 0 && plane->subsampling_y == 0);

  struct buf_2d *const dst = &plane->dst;
  uint8_t *buf = dst->buf;
  const int bitmask = 0xffff;
  int r;
#if CONFIG_EXT_PARTITION
  int sb_count = 1;
  if (cm->sb_size == BLOCK_128X128) sb_count = 2;
#endif

  const LoopFilterMask *mask = &lfm[0];
  const FilterMaskY *mask_16x16 = &mask->above_y[TX_16X16];
  const FilterMaskY *mask_8x8 = &mask->above_y[TX_8X8];
  const FilterMaskY *mask_4x4 = &mask->above_y[TX_4X4];
  unsigned int mask_16x16_r, mask_8x8_r, mask_4x4_r;

#if CONFIG_EXT_PARTITION
  uint8_t *buf2 = buf + (MAX_SB_SIZE >> 1);
  const LoopFilterMask *mask2 = &lfm[1];
  const FilterMaskY *mask2_16x16 = &mask2->above_y[TX_16X16];
  const FilterMaskY *mask2_8x8 = &mask2->above_y[TX_8X8];
  const FilterMaskY *mask2_4x4 = &mask2->above_y[TX_4X4];
  unsigned int mask2_16x16_r, mask2_8x8_r, mask2_4x4_r;
#endif

  const int mod = MAX_MIB_SIZE >> 1;
  const int yuv_comp = 0;

  for (r = 0; r < cm->mib_size && mi_row + r < cm->mi_rows; r++) {
    const int index = (r % mod) >> 2;
    const int shift = (r % 4) << 4;

    if (mi_row + r == 0) {
      mask_16x16_r = 0;
      mask_8x8_r = 0;
      mask_4x4_r = 0;
#if CONFIG_EXT_PARTITION
      mask2_16x16_r = 0;
      mask2_8x8_r = 0;
      mask2_4x4_r = 0;
#endif
    } else {
      mask_16x16_r = (mask_16x16->bits[index] >> shift) & bitmask;
      mask_8x8_r = (mask_8x8->bits[index] >> shift) & bitmask;
      mask_4x4_r = (mask_4x4->bits[index] >> shift) & bitmask;
#if CONFIG_EXT_PARTITION
      if (sb_count > 1) {
        mask2_16x16_r = (mask2_16x16->bits[index] >> shift) & bitmask;
        mask2_8x8_r = (mask2_8x8->bits[index] >> shift) & bitmask;
        mask2_4x4_r = (mask2_4x4->bits[index] >> shift) & bitmask;
      }
#endif
    }

    if (cm->use_highbitdepth)
      highbd_filter_selectively_horiz2(
          CONVERT_TO_SHORTPTR(buf), yuv_comp, dst->stride, mask_16x16_r,
          mask_8x8_r, mask_4x4_r, &cm->lf_info, &mask->lfl_y_hor[r % mod][0],
          (int)cm->bit_depth);
    else
      filter_selectively_horiz2(buf, yuv_comp, dst->stride, mask_16x16_r,
                                mask_8x8_r, mask_4x4_r, &cm->lf_info,
                                &mask->lfl_y_hor[r % mod][0]);

#if CONFIG_EXT_PARTITION
    if (sb_count > 1) {
      if (cm->use_highbitdepth)
        highbd_filter_selectively_horiz2(
            CONVERT_TO_SHORTPTR(buf2), yuv_comp, dst->stride, mask2_16x16_r,
            mask2_8x8_r, mask2_4x4_r, &cm->lf_info,
            &mask2->lfl_y_hor[r % mod][0], (int)cm->bit_depth);
      else
        filter_selectively_horiz2(buf2, yuv_comp, dst->stride, mask2_16x16_r,
                                  mask2_8x8_r, mask2_4x4_r, &cm->lf_info,
                                  &mask2->lfl_y_hor[r % mod][0]);
    }
#endif

    buf += MI_SIZE * dst->stride;
#if CONFIG_EXT_PARTITION
    buf2 += MI_SIZE * dst->stride;
    if (sb_count > 1 && r + 1 == MI_SIZE_64X64) {
      mask += 2;
      mask2 += 2;

      mask_16x16 = &mask->above_y[TX_16X16];
      mask_8x8 = &mask->above_y[TX_8X8];
      mask_4x4 = &mask->above_y[TX_4X4];

      mask2_16x16 = &mask2->above_y[TX_16X16];
      mask2_8x8 = &mask2->above_y[TX_8X8];
      mask2_4x4 = &mask2->above_y[TX_4X4];
    }
#endif
  }
}

void av1_filter_block_plane_ss11_u_ver(AV1_COMMON *const cm,
                                       struct macroblockd_plane *const plane,
                                       int mi_row, LoopFilterMask *lfm) {
  assert(plane->subsampling_x == 1 && plane->subsampling_y == 1);
  assert(plane->plane_type == PLANE_TYPE_UV);

  struct buf_2d *const dst = &plane->dst;
  int r;
#if CONFIG_EXT_PARTITION
  int sb_count = 1;
  if (cm->sb_size == BLOCK_128X128) sb_count = 2;
#endif

  uint8_t *buf = dst->buf;
  LoopFilterMask *mask = &lfm[0];
  FilterMaskUV mask_16x16 = mask->left_u[TX_16X16];
  FilterMaskUV mask_8x8 = mask->left_u[TX_8X8];
  FilterMaskUV mask_4x4 = mask->left_u[TX_4X4];
  uint32_t mask_16x16_l = 0;
  uint32_t mask_8x8_l = 0;
  uint32_t mask_4x4_l = 0;
#if CONFIG_EXT_PARTITION
  uint8_t *buf2 = buf + (MAX_SB_SIZE >> 2);
  LoopFilterMask *mask2 = &lfm[1];
  FilterMaskUV mask2_16x16 = mask2->left_u[TX_16X16];
  FilterMaskUV mask2_8x8 = mask2->left_u[TX_8X8];
  FilterMaskUV mask2_4x4 = mask2->left_u[TX_4X4];
  uint32_t mask2_16x16_l = 0;
  uint32_t mask2_8x8_l = 0;
  uint32_t mask2_4x4_l = 0;
#endif

  const int bufoffset = (MI_SIZE * dst->stride) << 1;
  const int mod = MAX_MIB_SIZE >> 1;
  const int yuv_comp = 1;

  // Vertical pass: do 2 rows at one time
  for (r = 0; r < cm->mib_size && mi_row + r < cm->mi_rows; r += 4) {
    const int shift = (r % mod) << 2;

    mask_16x16_l = (mask_16x16 >> shift) & 0xffff;
    mask_8x8_l = (mask_8x8 >> shift) & 0xffff;
    mask_4x4_l = (mask_4x4 >> shift) & 0xffff;
#if CONFIG_EXT_PARTITION
    if (sb_count > 1) {
      mask2_16x16_l = (mask2_16x16 >> shift) & 0xffff;
      mask2_8x8_l = (mask2_8x8 >> shift) & 0xffff;
      mask2_4x4_l = (mask2_4x4 >> shift) & 0xffff;
    }
#endif

    if (cm->use_highbitdepth)
      highbd_filter_selectively_vert_row2(
          plane->subsampling_x, CONVERT_TO_SHORTPTR(buf), yuv_comp, dst->stride,
          mask_16x16_l, mask_8x8_l, mask_4x4_l, &cm->lf_info,
          &mask->lfl_u[(r % mod) >> 1][0], (int)cm->bit_depth);
    else
      filter_selectively_vert_row2(
          plane->subsampling_x, buf, yuv_comp, dst->stride, mask_16x16_l,
          mask_8x8_l, mask_4x4_l, &cm->lf_info, &lfm->lfl_u[(r % mod) >> 1][0]);

// Disable filtering on the leftmost column.
#if CONFIG_EXT_PARTITION
    if (sb_count > 1) {
      if (cm->use_highbitdepth)
        highbd_filter_selectively_vert_row2(
            plane->subsampling_x, CONVERT_TO_SHORTPTR(buf2), yuv_comp,
            dst->stride, mask2_16x16_l, mask2_8x8_l, mask2_4x4_l, &cm->lf_info,
            &mask2->lfl_u[(r % mod) >> 1][0], (int)cm->bit_depth);
      else
        filter_selectively_vert_row2(plane->subsampling_x, buf2, yuv_comp,
                                     dst->stride, mask2_16x16_l, mask2_8x8_l,
                                     mask2_4x4_l, &cm->lf_info,
                                     &mask2->lfl_u[(r % mod) >> 1][0]);
    }
#endif

    buf += bufoffset;
#if CONFIG_EXT_PARTITION
    buf2 += bufoffset;
    if (sb_count > 1 && r + 4 == MI_SIZE_64X64) {
      mask += 2;
      mask2 += 2;

      mask_16x16 = mask->left_u[TX_16X16];
      mask_8x8 = mask->left_u[TX_8X8];
      mask_4x4 = mask->left_u[TX_4X4];

      mask2_16x16 = mask2->left_u[TX_16X16];
      mask2_8x8 = mask2->left_u[TX_8X8];
      mask2_4x4 = mask2->left_u[TX_4X4];
    }
#endif
  }
}

void av1_filter_block_plane_ss11_v_ver(AV1_COMMON *const cm,
                                       struct macroblockd_plane *const plane,
                                       int mi_row, LoopFilterMask *lfm) {
  assert(plane->subsampling_x == 1 && plane->subsampling_y == 1);
  assert(plane->plane_type == PLANE_TYPE_UV);

  struct buf_2d *const dst = &plane->dst;
  int r;
#if CONFIG_EXT_PARTITION
  int sb_count = 1;
  if (cm->sb_size == BLOCK_128X128) sb_count = 2;
#endif

  uint8_t *buf = dst->buf;
  LoopFilterMask *mask = &lfm[0];
  FilterMaskUV mask_16x16 = mask->left_v[TX_16X16];
  FilterMaskUV mask_8x8 = mask->left_v[TX_8X8];
  FilterMaskUV mask_4x4 = mask->left_v[TX_4X4];
  uint32_t mask_16x16_l = 0;
  uint32_t mask_8x8_l = 0;
  uint32_t mask_4x4_l = 0;
#if CONFIG_EXT_PARTITION
  uint8_t *buf2 = buf + (MAX_SB_SIZE >> 2);
  LoopFilterMask *mask2 = &lfm[1];
  FilterMaskUV mask2_16x16 = mask2->left_v[TX_16X16];
  FilterMaskUV mask2_8x8 = mask2->left_v[TX_8X8];
  FilterMaskUV mask2_4x4 = mask2->left_v[TX_4X4];
  uint32_t mask2_16x16_l = 0;
  uint32_t mask2_8x8_l = 0;
  uint32_t mask2_4x4_l = 0;
#endif

  const int bufoffset = (MI_SIZE * dst->stride) << 1;
  const int mod = MAX_MIB_SIZE >> 1;
  const int yuv_comp = 2;

  // Vertical pass: do 2 rows at one time
  for (r = 0; r < cm->mib_size && mi_row + r < cm->mi_rows; r += 4) {
    const int shift = (r % mod) << 2;

    mask_16x16_l = (mask_16x16 >> shift) & 0xffff;
    mask_8x8_l = (mask_8x8 >> shift) & 0xffff;
    mask_4x4_l = (mask_4x4 >> shift) & 0xffff;
#if CONFIG_EXT_PARTITION
    if (sb_count > 1) {
      mask2_16x16_l = (mask2_16x16 >> shift) & 0xffff;
      mask2_8x8_l = (mask2_8x8 >> shift) & 0xffff;
      mask2_4x4_l = (mask2_4x4 >> shift) & 0xffff;
    }
#endif

    if (cm->use_highbitdepth)
      highbd_filter_selectively_vert_row2(
          plane->subsampling_x, CONVERT_TO_SHORTPTR(buf), yuv_comp, dst->stride,
          mask_16x16_l, mask_8x8_l, mask_4x4_l, &cm->lf_info,
          &mask->lfl_v[(r % mod) >> 1][0], (int)cm->bit_depth);
    else
      filter_selectively_vert_row2(
          plane->subsampling_x, buf, yuv_comp, dst->stride, mask_16x16_l,
          mask_8x8_l, mask_4x4_l, &cm->lf_info, &lfm->lfl_v[(r % mod) >> 1][0]);

// Disable filtering on the leftmost column.
#if CONFIG_EXT_PARTITION
    if (sb_count > 1) {
      if (cm->use_highbitdepth)
        highbd_filter_selectively_vert_row2(
            plane->subsampling_x, CONVERT_TO_SHORTPTR(buf2), yuv_comp,
            dst->stride, mask2_16x16_l, mask2_8x8_l, mask2_4x4_l, &cm->lf_info,
            &mask2->lfl_v[(r % mod) >> 1][0], (int)cm->bit_depth);
      else
        filter_selectively_vert_row2(plane->subsampling_x, buf2, yuv_comp,
                                     dst->stride, mask2_16x16_l, mask2_8x8_l,
                                     mask2_4x4_l, &cm->lf_info,
                                     &mask2->lfl_v[(r % mod) >> 1][0]);
    }
#endif

    buf += bufoffset;
#if CONFIG_EXT_PARTITION
    buf2 += bufoffset;
    if (sb_count > 1 && r + 4 == MI_SIZE_64X64) {
      mask += 2;
      mask2 += 2;

      mask_16x16 = mask->left_v[TX_16X16];
      mask_8x8 = mask->left_v[TX_8X8];
      mask_4x4 = mask->left_v[TX_4X4];

      mask2_16x16 = mask2->left_v[TX_16X16];
      mask2_8x8 = mask2->left_v[TX_8X8];
      mask2_4x4 = mask2->left_v[TX_4X4];
    }
#endif
  }
}

void av1_filter_block_plane_ss11_u_hor(AV1_COMMON *const cm,
                                       struct macroblockd_plane *const plane,
                                       int mi_row, LoopFilterMask *lfm) {
  assert(plane->subsampling_x == 1 && plane->subsampling_y == 1);

  struct buf_2d *const dst = &plane->dst;
  uint8_t *buf = dst->buf;
  int r;
#if CONFIG_EXT_PARTITION
  int sb_count = 1;
  if (cm->sb_size == BLOCK_128X128) sb_count = 2;
#endif

  LoopFilterMask *mask = &lfm[0];
  FilterMaskUV mask_16x16 = mask->above_u[TX_16X16];
  FilterMaskUV mask_8x8 = mask->above_u[TX_8X8];
  FilterMaskUV mask_4x4 = mask->above_u[TX_4X4];
  uint32_t mask_16x16_a, mask_8x8_a, mask_4x4_a;

#if CONFIG_EXT_PARTITION
  uint8_t *buf2 = buf + (MAX_SB_SIZE >> 2);
  LoopFilterMask *mask2 = &lfm[1];
  FilterMaskUV mask2_16x16 = mask2->above_u[TX_16X16];
  FilterMaskUV mask2_8x8 = mask2->above_u[TX_8X8];
  FilterMaskUV mask2_4x4 = mask2->above_u[TX_4X4];
  uint32_t mask2_16x16_a, mask2_8x8_a, mask2_4x4_a;
#endif

  const uint32_t bitmask = 0xff;
  const int mod = MAX_MIB_SIZE >> 1;
  const int yuv_comp = 1;

  // re-porpulate the filter level for uv, same as the code for vertical
  // filter in av1_filter_block_plane_ss11_ver
  for (r = 0; r < cm->mib_size && mi_row + r < cm->mi_rows; r += 2) {
    const int shift = (r % mod) << 2;

    if (mi_row + r == 0) {
      mask_16x16_a = 0;
      mask_8x8_a = 0;
      mask_4x4_a = 0;
#if CONFIG_EXT_PARTITION
      mask2_16x16_a = 0;
      mask2_8x8_a = 0;
      mask2_4x4_a = 0;
#endif
    } else {
      mask_16x16_a = (mask_16x16 >> shift) & bitmask;
      mask_8x8_a = (mask_8x8 >> shift) & bitmask;
      mask_4x4_a = (mask_4x4 >> shift) & bitmask;
#if CONFIG_EXT_PARTITION
      if (sb_count > 1) {
        mask2_16x16_a = (mask2_16x16 >> shift) & bitmask;
        mask2_8x8_a = (mask2_8x8 >> shift) & bitmask;
        mask2_4x4_a = (mask2_4x4 >> shift) & bitmask;
      }
#endif
    }

    if (cm->use_highbitdepth)
      highbd_filter_selectively_horiz2(
          CONVERT_TO_SHORTPTR(buf), yuv_comp, dst->stride, mask_16x16_a,
          mask_8x8_a, mask_4x4_a, &cm->lf_info, &mask->lfl_u[(r % mod) >> 1][0],
          (int)cm->bit_depth);
    else
      filter_selectively_horiz2(buf, yuv_comp, dst->stride, mask_16x16_a,
                                mask_8x8_a, mask_4x4_a, &cm->lf_info,
                                &mask->lfl_u[(r % mod) >> 1][0]);

#if CONFIG_EXT_PARTITION
    if (sb_count > 1) {
      if (cm->use_highbitdepth)
        highbd_filter_selectively_horiz2(
            CONVERT_TO_SHORTPTR(buf), yuv_comp, dst->stride, mask2_16x16_a,
            mask2_8x8_a, mask2_4x4_a, &cm->lf_info,
            &mask2->lfl_u[(r % mod) >> 1][0], (int)cm->bit_depth);
      else
        filter_selectively_horiz2(buf, yuv_comp, dst->stride, mask2_16x16_a,
                                  mask2_8x8_a, mask2_4x4_a, &cm->lf_info,
                                  &mask2->lfl_u[(r % mod) >> 1][0]);
    }
#endif

    buf += MI_SIZE * dst->stride;
#if CONFIG_EXT_PARTITION
    buf2 += MI_SIZE * dst->stride;
    if (sb_count > 1 && r + 2 == MI_SIZE_64X64) {
      mask += 2;
      mask2 += 2;

      mask_16x16 = mask->above_u[TX_16X16];
      mask_8x8 = mask->above_u[TX_8X8];
      mask_4x4 = mask->above_u[TX_4X4];

      mask2_16x16 = mask2->above_u[TX_16X16];
      mask2_8x8 = mask2->above_u[TX_8X8];
      mask2_4x4 = mask2->above_u[TX_4X4];
    }
#endif
  }
}

void av1_filter_block_plane_ss11_v_hor(AV1_COMMON *const cm,
                                       struct macroblockd_plane *const plane,
                                       int mi_row, LoopFilterMask *lfm) {
  assert(plane->subsampling_x == 1 && plane->subsampling_y == 1);

  struct buf_2d *const dst = &plane->dst;
  uint8_t *buf = dst->buf;
  int r;
#if CONFIG_EXT_PARTITION
  int sb_count = 1;
  if (cm->sb_size == BLOCK_128X128) sb_count = 2;
#endif

  LoopFilterMask *mask = &lfm[0];
  FilterMaskUV mask_16x16 = mask->above_v[TX_16X16];
  FilterMaskUV mask_8x8 = mask->above_v[TX_8X8];
  FilterMaskUV mask_4x4 = mask->above_v[TX_4X4];
  uint32_t mask_16x16_a, mask_8x8_a, mask_4x4_a;

#if CONFIG_EXT_PARTITION
  uint8_t *buf2 = buf + (MAX_SB_SIZE >> 2);
  LoopFilterMask *mask2 = &lfm[1];
  FilterMaskUV mask2_16x16 = mask2->above_v[TX_16X16];
  FilterMaskUV mask2_8x8 = mask2->above_v[TX_8X8];
  FilterMaskUV mask2_4x4 = mask2->above_v[TX_4X4];
  uint32_t mask2_16x16_a, mask2_8x8_a, mask2_4x4_a;
#endif

  const uint32_t bitmask = 0xff;
  const int mod = MAX_MIB_SIZE >> 1;
  const int yuv_comp = 2;

  // re-porpulate the filter level for uv, same as the code for vertical
  // filter in av1_filter_block_plane_ss11_ver
  for (r = 0; r < cm->mib_size && mi_row + r < cm->mi_rows; r += 2) {
    const int shift = (r % mod) << 2;

    if (mi_row + r == 0) {
      mask_16x16_a = 0;
      mask_8x8_a = 0;
      mask_4x4_a = 0;
#if CONFIG_EXT_PARTITION
      mask2_16x16_a = 0;
      mask2_8x8_a = 0;
      mask2_4x4_a = 0;
#endif
    } else {
      mask_16x16_a = (mask_16x16 >> shift) & bitmask;
      mask_8x8_a = (mask_8x8 >> shift) & bitmask;
      mask_4x4_a = (mask_4x4 >> shift) & bitmask;
#if CONFIG_EXT_PARTITION
      if (sb_count > 1) {
        mask2_16x16_a = (mask2_16x16 >> shift) & bitmask;
        mask2_8x8_a = (mask2_8x8 >> shift) & bitmask;
        mask2_4x4_a = (mask2_4x4 >> shift) & bitmask;
      }
#endif
    }

    if (cm->use_highbitdepth)
      highbd_filter_selectively_horiz2(
          CONVERT_TO_SHORTPTR(buf), yuv_comp, dst->stride, mask_16x16_a,
          mask_8x8_a, mask_4x4_a, &cm->lf_info, &mask->lfl_v[(r % mod) >> 1][0],
          (int)cm->bit_depth);
    else
      filter_selectively_horiz2(buf, yuv_comp, dst->stride, mask_16x16_a,
                                mask_8x8_a, mask_4x4_a, &cm->lf_info,
                                &mask->lfl_v[(r % mod) >> 1][0]);

#if CONFIG_EXT_PARTITION
    if (sb_count > 1) {
      if (cm->use_highbitdepth)
        highbd_filter_selectively_horiz2(
            CONVERT_TO_SHORTPTR(buf), yuv_comp, dst->stride, mask2_16x16_a,
            mask2_8x8_a, mask2_4x4_a, &cm->lf_info,
            &mask2->lfl_v[(r % mod) >> 1][0], (int)cm->bit_depth);
      else
        filter_selectively_horiz2(buf, yuv_comp, dst->stride, mask2_16x16_a,
                                  mask2_8x8_a, mask2_4x4_a, &cm->lf_info,
                                  &mask2->lfl_v[(r % mod) >> 1][0]);
    }
#endif

    buf += MI_SIZE * dst->stride;
#if CONFIG_EXT_PARTITION
    buf2 += MI_SIZE * dst->stride;
    if (sb_count > 1 && r + 2 == MI_SIZE_64X64) {
      mask += 2;
      mask2 += 2;

      mask_16x16 = mask->above_v[TX_16X16];
      mask_8x8 = mask->above_v[TX_8X8];
      mask_4x4 = mask->above_v[TX_4X4];

      mask2_16x16 = mask2->above_v[TX_16X16];
      mask2_8x8 = mask2->above_v[TX_8X8];
      mask2_4x4 = mask2->above_v[TX_4X4];
    }
#endif
  }
}

#if CONFIG_PARALLEL_DEBLOCKING
// Note:
// This loopfilter implementation builds bit mask from var_tx and cb4x4.
// Loop filter uses cb4x4 filter which is the same as parallel_deblocking.

#if !USE_BITMASK_FROM_VAR_TX
typedef enum EDGE_DIR { VERT_EDGE = 0, HORZ_EDGE = 1, NUM_EDGE_DIRS } EDGE_DIR;
static const uint32_t av1_prediction_masks[NUM_EDGE_DIRS][BLOCK_SIZES_ALL] = {
  // mask for vertical edges filtering
  {
      4 - 1,   // BLOCK_4X4
      4 - 1,   // BLOCK_4X8
      8 - 1,   // BLOCK_8X4
      8 - 1,   // BLOCK_8X8
      8 - 1,   // BLOCK_8X16
      16 - 1,  // BLOCK_16X8
      16 - 1,  // BLOCK_16X16
      16 - 1,  // BLOCK_16X32
      32 - 1,  // BLOCK_32X16
      32 - 1,  // BLOCK_32X32
      32 - 1,  // BLOCK_32X64
      64 - 1,  // BLOCK_64X32
      64 - 1,  // BLOCK_64X64
#if CONFIG_EXT_PARTITION
      64 - 1,   // BLOCK_64X128
      128 - 1,  // BLOCK_128X64
      128 - 1,  // BLOCK_128X128
#endif          // CONFIG_EXT_PARTITION
      4 - 1,    // BLOCK_4X16,
      16 - 1,   // BLOCK_16X4,
      8 - 1,    // BLOCK_8X32,
      32 - 1,   // BLOCK_32X8,
      16 - 1,   // BLOCK_16X64,
      64 - 1,   // BLOCK_64X16
#if CONFIG_EXT_PARTITION
      32 - 1,   // BLOCK_32X128
      128 - 1,  // BLOCK_128X32
#endif          // CONFIG_EXT_PARTITION
  },
  // mask for horizontal edges filtering
  {
      4 - 1,   // BLOCK_4X4
      8 - 1,   // BLOCK_4X8
      4 - 1,   // BLOCK_8X4
      8 - 1,   // BLOCK_8X8
      16 - 1,  // BLOCK_8X16
      8 - 1,   // BLOCK_16X8
      16 - 1,  // BLOCK_16X16
      32 - 1,  // BLOCK_16X32
      16 - 1,  // BLOCK_32X16
      32 - 1,  // BLOCK_32X32
      64 - 1,  // BLOCK_32X64
      32 - 1,  // BLOCK_64X32
      64 - 1,  // BLOCK_64X64
#if CONFIG_EXT_PARTITION
      128 - 1,  // BLOCK_64X128
      64 - 1,   // BLOCK_128X64
      128 - 1,  // BLOCK_128X128
#endif          // CONFIG_EXT_PARTITION
      16 - 1,   // BLOCK_4X16,
      4 - 1,    // BLOCK_16X4,
      32 - 1,   // BLOCK_8X32,
      8 - 1,    // BLOCK_32X8,
      64 - 1,   // BLOCK_16X64,
      16 - 1,   // BLOCK_64X16
#if CONFIG_EXT_PARTITION
      128 - 1,  // BLOCK_32X128
      32 - 1,   // BLOCK_128X32
#endif          // CONFIG_EXT_PARTITION
  },
};

static const uint32_t av1_transform_masks[NUM_EDGE_DIRS][TX_SIZES_ALL] = {
  {
      4 - 1,   // TX_4X4
      8 - 1,   // TX_8X8
      16 - 1,  // TX_16X16
      32 - 1,  // TX_32X32
#if CONFIG_TX64X64
      64 - 1,  // TX_64X64
#endif         // CONFIG_TX64X64
      4 - 1,   // TX_4X8
      8 - 1,   // TX_8X4
      8 - 1,   // TX_8X16
      16 - 1,  // TX_16X8
      16 - 1,  // TX_16X32
      32 - 1,  // TX_32X16
#if CONFIG_TX64X64
      32 - 1,  // TX_32X64
      64 - 1,  // TX_64X32
#endif         // CONFIG_TX64X64
      4 - 1,   // TX_4X16
      16 - 1,  // TX_16X4
      8 - 1,   // TX_8X32
      32 - 1,  // TX_32X8
#if CONFIG_TX64X64
      16 - 1,  // TX_16X64
      64 - 1,  // TX_64X16
#endif         // CONFIG_TX64X64
  },
  {
      4 - 1,   // TX_4X4
      8 - 1,   // TX_8X8
      16 - 1,  // TX_16X16
      32 - 1,  // TX_32X32
#if CONFIG_TX64X64
      64 - 1,  // TX_64X64
#endif         // CONFIG_TX64X64
      8 - 1,   // TX_4X8
      4 - 1,   // TX_8X4
      16 - 1,  // TX_8X16
      8 - 1,   // TX_16X8
      32 - 1,  // TX_16X32
      16 - 1,  // TX_32X16
#if CONFIG_TX64X64
      64 - 1,  // TX_32X64
      32 - 1,  // TX_64X32
#endif         // CONFIG_TX64X64
      16 - 1,  // TX_4X16
      4 - 1,   // TX_16X4
      32 - 1,  // TX_8X32
      8 - 1,   // TX_32X8
#if CONFIG_TX64X64
      64 - 1,  // TX_16X64
      16 - 1,  // TX_64X16
#endif         // CONFIG_TX64X64
  }
};

static TX_SIZE av1_get_transform_size(
    const MODE_INFO *const mi, const EDGE_DIR edge_dir, const int mi_row,
    const int mi_col, const int plane,
    const struct macroblockd_plane *plane_ptr) {
  const MB_MODE_INFO *mbmi = &mi->mbmi;
  TX_SIZE tx_size = (plane == AOM_PLANE_Y)
                        ? mbmi->tx_size
                        : av1_get_uv_tx_size(mbmi, plane_ptr->subsampling_x,
                                             plane_ptr->subsampling_y);
  assert(tx_size < TX_SIZES_ALL);

  // mi_row and mi_col is the absolute position of the MI block.
  // idx_c and idx_r is the relative offset of the MI within the super block
  // c and r is the relative offset of the 8x8 block within the supert block
  // blk_row and block_col is the relative offset of the current 8x8 block
  // within the current partition.
  const int idx_c = mi_col & MAX_MIB_MASK;
  const int idx_r = mi_row & MAX_MIB_MASK;
  const int c = idx_c >> mi_width_log2_lookup[BLOCK_8X8];
  const int r = idx_r >> mi_height_log2_lookup[BLOCK_8X8];
  const BLOCK_SIZE sb_type = mi->mbmi.sb_type;
  const int blk_row = r & (num_8x8_blocks_high_lookup[sb_type] - 1);
  const int blk_col = c & (num_8x8_blocks_wide_lookup[sb_type] - 1);

  if (is_inter_block(mbmi) && !mbmi->skip) {
    const int tx_row_idx =
        (blk_row * mi_size_high[BLOCK_8X8] << TX_UNIT_HIGH_LOG2) >> 1;
    const int tx_col_idx =
        (blk_col * mi_size_wide[BLOCK_8X8] << TX_UNIT_WIDE_LOG2) >> 1;
    const TX_SIZE mb_tx_size = mbmi->inter_tx_size[tx_row_idx][tx_col_idx];

    assert(mb_tx_size < TX_SIZES_ALL);

    tx_size = (plane == AOM_PLANE_Y)
                  ? mb_tx_size
                  : av1_get_uv_tx_size(mbmi, plane_ptr->subsampling_x,
                                       plane_ptr->subsampling_y);
    assert(tx_size < TX_SIZES_ALL);
  }

  // since in case of chrominance or non-square transorm need to convert
  // transform size into transform size in particular direction.
  // for vertical edge, filter direction is horizontal, for horizontal
  // edge, filter direction is vertical.
  tx_size = (VERT_EDGE == edge_dir) ? txsize_horz_map[tx_size]
                                    : txsize_vert_map[tx_size];
  return tx_size;
}

typedef struct AV1_DEBLOCKING_PARAMETERS {
  // length of the filter applied to the outer edge
  uint32_t filter_length;
  // length of the filter applied to the inner edge
  uint32_t filter_length_internal;
  // deblocking limits
  const uint8_t *lim;
  const uint8_t *mblim;
  const uint8_t *hev_thr;
} AV1_DEBLOCKING_PARAMETERS;

static void set_lpf_parameters(
    AV1_DEBLOCKING_PARAMETERS *const params, const ptrdiff_t mode_step,
    const AV1_COMMON *const cm, const EDGE_DIR edge_dir, const uint32_t x,
    const uint32_t y, const int plane,
    const struct macroblockd_plane *const plane_ptr) {
  // reset to initial values
  params->filter_length = 0;
  params->filter_length_internal = 0;

  // no deblocking is required
  const uint32_t width = plane_ptr->dst.width;
  const uint32_t height = plane_ptr->dst.height;
  if ((width <= x) || (height <= y)) {
    return;
  }

  const uint32_t scale_horz = plane_ptr->subsampling_x;
  const uint32_t scale_vert = plane_ptr->subsampling_y;
  // for sub8x8 block, chroma prediction mode is obtained from the bottom/right
  // mi structure of the co-located 8x8 luma block. so for chroma plane, mi_row
  // and mi_col should map to the bottom/right mi structure, i.e, both mi_row
  // and mi_col should be odd number for chroma plane.
  const int mi_row = scale_vert | ((y << scale_vert) >> MI_SIZE_LOG2);
  const int mi_col = scale_horz | ((x << scale_horz) >> MI_SIZE_LOG2);
  MODE_INFO **mi = cm->mi_grid_visible + mi_row * cm->mi_stride + mi_col;
  const MB_MODE_INFO *mbmi = &mi[0]->mbmi;

  {
    const TX_SIZE ts = av1_get_transform_size(mi[0], edge_dir, mi_row, mi_col,
                                              plane, plane_ptr);

#if CONFIG_EXT_DELTA_Q
#if CONFIG_LOOPFILTER_LEVEL
    const uint32_t curr_level =
        get_filter_level(cm, &cm->lf_info, edge_dir, plane, mbmi);
#else
    const uint32_t curr_level = get_filter_level(cm, &cm->lf_info, mbmi);
#endif
#else
    const uint32_t curr_level = get_filter_level(&cm->lf_info, mbmi);
#endif  // CONFIG_EXT_DELTA_Q

    const int curr_skipped = mbmi->skip && is_inter_block(mbmi);
    const uint32_t coord = (VERT_EDGE == edge_dir) ? (x) : (y);
    uint32_t level = curr_level;
    // prepare outer edge parameters. deblock the edge if it's an edge of a TU
    if (coord) {
#if CONFIG_LOOPFILTERING_ACROSS_TILES || CONFIG_LOOPFILTERING_ACROSS_TILES_EXT
      // Note: For sub8x8 blocks, we need to look at the top-left mi unit in
      // order
      // to extract the correct boundary information.
      const int mi_row_bound = ((y << scale_vert) >> MI_SIZE_LOG2);
      const int mi_col_bound = ((x << scale_horz) >> MI_SIZE_LOG2);
      BOUNDARY_TYPE *const bi =
          cm->boundary_info + mi_row_bound * cm->mi_stride + mi_col_bound;
      // here, assuming bounfary_info is set correctly based on the
      // loop_filter_across_tiles_enabled flag, i.e, tile boundary should
      // only be set to true when this flag is set to 0.
      int left_boundary = (*bi & TILE_LEFT_BOUNDARY);
      int top_boundary = (*bi & TILE_ABOVE_BOUNDARY);
      if (((VERT_EDGE == edge_dir) && (0 == left_boundary)) ||
          ((HORZ_EDGE == edge_dir) && (0 == top_boundary)))
#endif  // CONFIG_LOOPFILTERING_ACROSS_TILES
      {
        const int32_t tu_edge =
            (coord & av1_transform_masks[edge_dir][ts]) ? (0) : (1);
        if (tu_edge) {
          const MODE_INFO *const mi_prev = *(mi - mode_step);
          const int pv_row =
              (VERT_EDGE == edge_dir) ? (mi_row) : (mi_row - (1 << scale_vert));
          const int pv_col =
              (VERT_EDGE == edge_dir) ? (mi_col - (1 << scale_horz)) : (mi_col);
          const TX_SIZE pv_ts = av1_get_transform_size(
              mi_prev, edge_dir, pv_row, pv_col, plane, plane_ptr);

#if CONFIG_EXT_DELTA_Q
#if CONFIG_LOOPFILTER_LEVEL
          const uint32_t pv_lvl = get_filter_level(cm, &cm->lf_info, edge_dir,
                                                   plane, &mi_prev->mbmi);
#else
          const uint32_t pv_lvl =
              get_filter_level(cm, &cm->lf_info, &mi_prev->mbmi);
#endif
#else
          const uint32_t pv_lvl =
              get_filter_level(&cm->lf_info, &mi_prev->mbmi);
#endif  // CONFIG_EXT_DELTA_Q

          const int pv_skip =
              mi_prev->mbmi.skip && is_inter_block(&mi_prev->mbmi);
          const int32_t pu_edge =
              (coord &
               av1_prediction_masks[edge_dir]
                                   [ss_size_lookup[mbmi->sb_type][scale_horz]
                                                  [scale_vert]])
                  ? (0)
                  : (1);
          // if the current and the previous blocks are skipped,
          // deblock the edge if the edge belongs to a PU's edge only.
          if ((curr_level || pv_lvl) &&
              (!pv_skip || !curr_skipped || pu_edge)) {
            const TX_SIZE min_ts = AOMMIN(ts, pv_ts);
            if (TX_4X4 >= min_ts) {
              params->filter_length = 4;
            } else if (TX_8X8 == min_ts) {
#if PARALLEL_DEBLOCKING_5_TAP_CHROMA
              if (plane != 0)
                params->filter_length = 6;
              else
#endif
                params->filter_length = 8;
            } else {
              params->filter_length = 16;
              // No wide filtering for chroma plane
              if (plane != 0) {
#if PARALLEL_DEBLOCKING_5_TAP_CHROMA
                params->filter_length = 6;
#else
                params->filter_length = 8;
#endif
              }
            }

            // update the level if the current block is skipped,
            // but the previous one is not
            level = (curr_level) ? (curr_level) : (pv_lvl);
          }
        }
      }

      // prepare common parameters
      if (params->filter_length || params->filter_length_internal) {
        const loop_filter_thresh *const limits = cm->lf_info.lfthr + level;
        params->lim = limits->lim;
        params->mblim = limits->mblim;
        params->hev_thr = limits->hev_thr;
      }
    }
  }
}

static void av1_filter_block_plane_vert(
    const AV1_COMMON *const cm, const int plane,
    const MACROBLOCKD_PLANE *const plane_ptr, const uint32_t mi_row,
    const uint32_t mi_col) {
  const int col_step = MI_SIZE >> MI_SIZE_LOG2;
  const int row_step = MI_SIZE >> MI_SIZE_LOG2;
  const uint32_t scale_horz = plane_ptr->subsampling_x;
  const uint32_t scale_vert = plane_ptr->subsampling_y;
  uint8_t *const dst_ptr = plane_ptr->dst.buf;
  const int dst_stride = plane_ptr->dst.stride;
  const int y_range = (MAX_MIB_SIZE >> scale_vert);
  const int x_range = (MAX_MIB_SIZE >> scale_horz);
  for (int y = 0; y < y_range; y += row_step) {
    uint8_t *p = dst_ptr + y * MI_SIZE * dst_stride;
    for (int x = 0; x < x_range; x += col_step) {
      // inner loop always filter vertical edges in a MI block. If MI size
      // is 8x8, it will filter the vertical edge aligned with a 8x8 block.
      // If 4x4 trasnform is used, it will then filter the internal edge
      //  aligned with a 4x4 block
      const uint32_t curr_x = ((mi_col * MI_SIZE) >> scale_horz) + x * MI_SIZE;
      const uint32_t curr_y = ((mi_row * MI_SIZE) >> scale_vert) + y * MI_SIZE;
      AV1_DEBLOCKING_PARAMETERS params;
      memset(&params, 0, sizeof(params));

      set_lpf_parameters(&params, ((ptrdiff_t)1 << scale_horz), cm, VERT_EDGE,
                         curr_x, curr_y, plane, plane_ptr);

      switch (params.filter_length) {
        // apply 4-tap filtering
        case 4:
          if (cm->use_highbitdepth)
            aom_highbd_lpf_vertical_4(CONVERT_TO_SHORTPTR(p), dst_stride,
                                      params.mblim, params.lim, params.hev_thr,
                                      cm->bit_depth);
          else
            aom_lpf_vertical_4(p, dst_stride, params.mblim, params.lim,
                               params.hev_thr);
          break;
#if PARALLEL_DEBLOCKING_5_TAP_CHROMA
        case 6:  // apply 6-tap filter for chroma plane only
          assert(plane != 0);
          if (cm->use_highbitdepth)
            aom_highbd_lpf_vertical_6_c(CONVERT_TO_SHORTPTR(p), dst_stride,
                                        params.mblim, params.lim,
                                        params.hev_thr, cm->bit_depth);
          else
            aom_lpf_vertical_6_c(p, dst_stride, params.mblim, params.lim,
                                 params.hev_thr);
          break;
#endif
        // apply 8-tap filtering
        case 8:
          if (cm->use_highbitdepth)
            aom_highbd_lpf_vertical_8(CONVERT_TO_SHORTPTR(p), dst_stride,
                                      params.mblim, params.lim, params.hev_thr,
                                      cm->bit_depth);
          else
            aom_lpf_vertical_8(p, dst_stride, params.mblim, params.lim,
                               params.hev_thr);
          break;
        // apply 16-tap filtering
        case 16:
          if (cm->use_highbitdepth)
#if CONFIG_DEBLOCK_13TAP
            // TODO(olah): Remove _c once SIMD for 13-tap is available
            aom_highbd_lpf_vertical_16_c(CONVERT_TO_SHORTPTR(p), dst_stride,
                                         params.mblim, params.lim,
                                         params.hev_thr, cm->bit_depth);
#else
            aom_highbd_lpf_vertical_16(CONVERT_TO_SHORTPTR(p), dst_stride,
                                       params.mblim, params.lim, params.hev_thr,
                                       cm->bit_depth);
#endif
          else
#if CONFIG_DEBLOCK_13TAP
            aom_lpf_vertical_16_c(p, dst_stride, params.mblim, params.lim,
                                  params.hev_thr);
#else
            aom_lpf_vertical_16(p, dst_stride, params.mblim, params.lim,
                                params.hev_thr);
#endif
          break;
        // no filtering
        default: break;
      }
      // process the internal edge
      if (params.filter_length_internal) {
        if (cm->use_highbitdepth)
          aom_highbd_lpf_vertical_4(CONVERT_TO_SHORTPTR(p + 4), dst_stride,
                                    params.mblim, params.lim, params.hev_thr,
                                    cm->bit_depth);
        else
          aom_lpf_vertical_4(p + 4, dst_stride, params.mblim, params.lim,
                             params.hev_thr);
      }
      // advance the destination pointer
      p += MI_SIZE;
    }
  }
}
#endif  // !USE_BITMASK_FROM_VAR_TX

#if !USE_BITMASK_FROM_VAR_TX
static void av1_filter_block_plane_horz(
    const AV1_COMMON *const cm, const int plane,
    const MACROBLOCKD_PLANE *const plane_ptr, const uint32_t mi_row,
    const uint32_t mi_col) {
  const int col_step = MI_SIZE >> MI_SIZE_LOG2;
  const int row_step = MI_SIZE >> MI_SIZE_LOG2;
  const uint32_t scale_horz = plane_ptr->subsampling_x;
  const uint32_t scale_vert = plane_ptr->subsampling_y;
  uint8_t *const dst_ptr = plane_ptr->dst.buf;
  const int dst_stride = plane_ptr->dst.stride;
  const int y_range = (MAX_MIB_SIZE >> scale_vert);
  const int x_range = (MAX_MIB_SIZE >> scale_horz);
  for (int y = 0; y < y_range; y += row_step) {
    uint8_t *p = dst_ptr + y * MI_SIZE * dst_stride;
    for (int x = 0; x < x_range; x += col_step) {
      // inner loop always filter vertical edges in a MI block. If MI size
      // is 8x8, it will first filter the vertical edge aligned with a 8x8
      // block. If 4x4 trasnform is used, it will then filter the internal
      // edge aligned with a 4x4 block
      const uint32_t curr_x = ((mi_col * MI_SIZE) >> scale_horz) + x * MI_SIZE;
      const uint32_t curr_y = ((mi_row * MI_SIZE) >> scale_vert) + y * MI_SIZE;
      AV1_DEBLOCKING_PARAMETERS params;
      memset(&params, 0, sizeof(params));

      set_lpf_parameters(&params, (cm->mi_stride << scale_vert), cm, HORZ_EDGE,
                         curr_x, curr_y, plane, plane_ptr);

      switch (params.filter_length) {
        // apply 4-tap filtering
        case 4:
          if (cm->use_highbitdepth)
            aom_highbd_lpf_horizontal_4(CONVERT_TO_SHORTPTR(p), dst_stride,
                                        params.mblim, params.lim,
                                        params.hev_thr, cm->bit_depth);
          else
            aom_lpf_horizontal_4(p, dst_stride, params.mblim, params.lim,
                                 params.hev_thr);
          break;
#if PARALLEL_DEBLOCKING_5_TAP_CHROMA
        // apply 6-tap filtering
        case 6:
          assert(plane != 0);
          if (cm->use_highbitdepth)
            aom_highbd_lpf_horizontal_6_c(CONVERT_TO_SHORTPTR(p), dst_stride,
                                          params.mblim, params.lim,
                                          params.hev_thr, cm->bit_depth);
          else
            aom_lpf_horizontal_6_c(p, dst_stride, params.mblim, params.lim,
                                   params.hev_thr);
          break;
#endif
        // apply 8-tap filtering
        case 8:
          if (cm->use_highbitdepth)
            aom_highbd_lpf_horizontal_8(CONVERT_TO_SHORTPTR(p), dst_stride,
                                        params.mblim, params.lim,
                                        params.hev_thr, cm->bit_depth);
          else
            aom_lpf_horizontal_8(p, dst_stride, params.mblim, params.lim,
                                 params.hev_thr);
          break;
        // apply 16-tap filtering
        case 16:
          if (cm->use_highbitdepth)
#if CONFIG_DEBLOCK_13TAP
            // TODO(olah): Remove _c once SIMD for 13-tap is available
            aom_highbd_lpf_horizontal_16_dual_c(
                CONVERT_TO_SHORTPTR(p), dst_stride, params.mblim, params.lim,
                params.hev_thr, cm->bit_depth);
#else
            aom_highbd_lpf_horizontal_16_dual(
                CONVERT_TO_SHORTPTR(p), dst_stride, params.mblim, params.lim,
                params.hev_thr, cm->bit_depth);
#endif
          else
#if CONFIG_DEBLOCK_13TAP
            aom_lpf_horizontal_16_dual_c(p, dst_stride, params.mblim,
                                         params.lim, params.hev_thr);
#else
            aom_lpf_horizontal_16_dual(p, dst_stride, params.mblim, params.lim,
                                       params.hev_thr);
#endif
          break;
        // no filtering
        default: break;
      }
      // process the internal edge
      if (params.filter_length_internal) {
        if (cm->use_highbitdepth)
          aom_highbd_lpf_horizontal_4(CONVERT_TO_SHORTPTR(p + 4 * dst_stride),
                                      dst_stride, params.mblim, params.lim,
                                      params.hev_thr, cm->bit_depth);
        else
          aom_lpf_horizontal_4(p + 4 * dst_stride, dst_stride, params.mblim,
                               params.lim, params.hev_thr);
      }
      // advance the destination pointer
      p += MI_SIZE;
    }
  }
}
#endif  // USE_BITMASK_FROM_VAR_TX
#else   // !CONFIG_PARALLEL_DEBLOCKING
#define USE_BITMASK_FROM_VAR_TX 0
#endif  // CONFIG_PARALLEL_DEBLOCKING

static INLINE enum lf_path get_loop_filter_path(
    int y_only, struct macroblockd_plane planes[MAX_MB_PLANE]) {
#if CONFIG_LOOPFILTER_LEVEL
  if (planes[y_only].subsampling_y == 1 && planes[y_only].subsampling_x == 1)
    return LF_PATH_420;
  else if (planes[y_only].subsampling_y == 0 &&
           planes[y_only].subsampling_x == 0)
    return LF_PATH_444;
  else
    return LF_PATH_SLOW;
#else
  if (y_only)
    return LF_PATH_444;
  else if (planes[1].subsampling_y == 1 && planes[1].subsampling_x == 1)
    return LF_PATH_420;
  else if (planes[1].subsampling_y == 0 && planes[1].subsampling_x == 0)
    return LF_PATH_444;
  else
    return LF_PATH_SLOW;
#endif
}

static INLINE void loop_filter_block_plane_ver(
    AV1_COMMON *cm, struct macroblockd_plane planes[MAX_MB_PLANE], int plane,
    MODE_INFO **mi, int mi_row, int mi_col, enum lf_path path,
    LpfMask *lfpmask) {
  LoopFilterMask *lfm = &lfpmask->lfm[0];
  if (plane == 0) {
    av1_filter_block_plane_ss00_ver(cm, &planes[0], mi_row, lfm);
  } else {
    switch (path) {
      case LF_PATH_420:
        if (plane == 1)
          av1_filter_block_plane_ss11_u_ver(cm, &planes[plane], mi_row, lfm);
        if (plane == 2)
          av1_filter_block_plane_ss11_v_ver(cm, &planes[plane], mi_row, lfm);
        break;
      case LF_PATH_444:
        av1_filter_block_plane_ss00_ver(cm, &planes[plane], mi_row, lfm);
        break;
      case LF_PATH_SLOW:
        av1_filter_block_plane_non420_ver(cm, &planes[plane], mi, mi_row,
                                          mi_col, plane);
        break;
    }
  }
}

static INLINE void loop_filter_block_plane_hor(
    AV1_COMMON *cm, struct macroblockd_plane planes[MAX_MB_PLANE], int plane,
    MODE_INFO **mi, int mi_row, int mi_col, enum lf_path path,
    LpfMask *lfpmask) {
  LoopFilterMask *lfm = &lfpmask->lfm[0];
  if (plane == 0) {
    av1_filter_block_plane_ss00_hor(cm, &planes[0], mi_row, lfm);
  } else {
    switch (path) {
      case LF_PATH_420:
        if (plane == 1)
          av1_filter_block_plane_ss11_u_hor(cm, &planes[plane], mi_row, lfm);
        if (plane == 2)
          av1_filter_block_plane_ss11_v_hor(cm, &planes[plane], mi_row, lfm);
        break;
      case LF_PATH_444:
        av1_filter_block_plane_ss00_hor(cm, &planes[plane], mi_row, lfm);
        break;
      case LF_PATH_SLOW:
        av1_filter_block_plane_non420_hor(cm, &planes[plane], mi, mi_row,
                                          mi_col, plane);
        break;
    }
  }

  // Since we do vertical filtering first then horizontal filtering,
  // we reset bitmask flag here to make next loop filtering has a
  // chance to build the new bitmask.
  lfpmask->is_setup = 0;
}
#if 0
static void print_luma(uint8_t *buf, int stride) {
  int i, j;
  uint64_t *p = (uint64_t *)buf;
  for (i = 0; i < 128; i++) {
    for (j = 0; j < 16; j++) {
      printf("%016llx\n", *p++);
    }
    buf += stride;
    p = (uint64_t *)buf;
    // printf("\n");
  }
}
#endif

void av1_loop_filter_rows(YV12_BUFFER_CONFIG *frame_buffer, AV1_COMMON *cm,
                          struct macroblockd_plane *planes, int start, int stop,
                          int y_only) {
  const int num_planes = av1_num_planes(cm);
#if CONFIG_LOOPFILTER_LEVEL
  // y_only no longer has its original meaning.
  // Here it means which plane to filter
  // when y_only = {0, 1, 2}, it means we are searching for filter level for
  // Y/U/V plane individually.
  const int plane_start = y_only;
  const int plane_end = plane_start + 1;
#else
  const int nplanes = y_only ? 1 : num_planes;
  const int plane_start = 0;
  const int plane_end = nplanes;
#endif  // CONFIG_LOOPFILTER_LEVEL
#if CONFIG_PARALLEL_DEBLOCKING && !USE_BITMASK_FROM_VAR_TX
  const int col_start = 0;
  const int col_end = cm->mi_cols;
#endif
  int mi_row, mi_col;
  int plane;

#if !CONFIG_PARALLEL_DEBLOCKING
  for (int i = 0; i < nplanes; ++i)
    memset(cm->top_txfm_context[i], TX_32X32, cm->mi_cols << TX_UNIT_WIDE_LOG2);
  for (mi_row = start; mi_row < stop; mi_row += cm->mib_size) {
    MODE_INFO **mi = cm->mi_grid_visible + mi_row * cm->mi_stride;
    for (int i = 0; i < nplanes; ++i)
      memset(cm->left_txfm_context[i], TX_32X32,
             MAX_MIB_SIZE << TX_UNIT_HIGH_LOG2);
    for (mi_col = 0; mi_col < cm->mi_cols; mi_col += cm->mib_size) {
      av1_setup_dst_planes(planes, cm->sb_size, frame_buffer, mi_row, mi_col);
      for (plane = plane_start; plane < plane_end; ++plane) {
        av1_filter_block_plane_non420_ver(cm, &planes[plane], mi + mi_col,
                                          mi_row, mi_col, plane);
        av1_filter_block_plane_non420_hor(cm, &planes[plane], mi + mi_col,
                                          mi_row, mi_col, plane);
      }
    }
  }
#else  // CONFIG_PARALLEL_DEBLOCKING
#if USE_BITMASK_FROM_VAR_TX
  enum lf_path path = get_loop_filter_path(y_only, planes);

  // filter all vertical edges in every super block
  for (mi_row = start; mi_row < stop; mi_row += cm->mib_size) {
    MODE_INFO **mi = cm->mi_grid_visible + mi_row * cm->mi_stride;
    for (mi_col = 0; mi_col < cm->mi_cols; mi_col += cm->mib_size) {
      av1_setup_dst_planes(planes, cm->sb_size, frame_buffer, mi_row, mi_col,
                           num_planes);

      // if (plane_start == 0) print_luma(planes[0].dst.buf,
      // planes[0].dst.stride);

      LpfMask *lpf_mask = get_loop_filter_mask(cm, mi_row, mi_col);

      av1_setup_mask(cm, mi_row, mi_col, mi + mi_col, plane_start, lpf_mask);

      for (plane = plane_start; plane < plane_end; ++plane) {
        loop_filter_block_plane_ver(cm, planes, plane, mi + mi_col, mi_row,
                                    mi_col, path, lpf_mask);
        // if (plane == 0) print_luma(planes[0].dst.buf, planes[0].dst.stride);
      }
    }
  }

  // filter all horizontal edges in every super block
  for (mi_row = start; mi_row < stop; mi_row += cm->mib_size) {
    MODE_INFO **mi = cm->mi_grid_visible + mi_row * cm->mi_stride;
    for (mi_col = 0; mi_col < cm->mi_cols; mi_col += cm->mib_size) {
      av1_setup_dst_planes(planes, cm->sb_size, frame_buffer, mi_row, mi_col,
                           num_planes);

      LpfMask *lpf_mask = get_loop_filter_mask(cm, mi_row, mi_col);

      for (plane = plane_start; plane < plane_end; ++plane) {
        loop_filter_block_plane_hor(cm, planes, plane, mi + mi_col, mi_row,
                                    mi_col, path, lpf_mask);
        // if (plane == 0) print_luma(planes[0].dst.buf, planes[0].dst.stride);
      }
    }
  }
#else
  // filter all vertical edges in every 64x64 super block
  for (mi_row = start; mi_row < stop; mi_row += MAX_MIB_SIZE) {
    for (mi_col = col_start; mi_col < col_end; mi_col += MAX_MIB_SIZE) {
      av1_setup_dst_planes(planes, cm->sb_size, frame_buffer, mi_row, mi_col,
                           num_planes);
      for (plane = plane_start; plane < plane_end; ++plane) {
        av1_filter_block_plane_vert(cm, plane, &planes[plane], mi_row, mi_col);
      }
    }
  }

  // filter all horizontal edges in every 64x64 super block
  for (mi_row = start; mi_row < stop; mi_row += MAX_MIB_SIZE) {
    for (mi_col = col_start; mi_col < col_end; mi_col += MAX_MIB_SIZE) {
      av1_setup_dst_planes(planes, cm->sb_size, frame_buffer, mi_row, mi_col,
                           num_planes);
      for (plane = plane_start; plane < plane_end; ++plane) {
        av1_filter_block_plane_horz(cm, plane, &planes[plane], mi_row, mi_col);
      }
    }
  }
#endif  // USE_BITMASK_FROM_VAR_TX
#endif  // !CONFIG_PARALLEL_DEBLOCKING
}

void av1_loop_filter_frame(YV12_BUFFER_CONFIG *frame, AV1_COMMON *cm,
                           MACROBLOCKD *xd, int frame_filter_level,
#if CONFIG_LOOPFILTER_LEVEL
                           int frame_filter_level_r,
#endif
                           int y_only, int partial_frame) {
  int start_mi_row, end_mi_row, mi_rows_to_filter;
#if CONFIG_EXT_DELTA_Q
#if CONFIG_LOOPFILTER_LEVEL
  int orig_filter_level[2] = { cm->lf.filter_level[0], cm->lf.filter_level[1] };
#else
  int orig_filter_level = cm->lf.filter_level;
#endif
#endif

#if CONFIG_LOOPFILTER_LEVEL
  if (!frame_filter_level && !frame_filter_level_r) return;
#else
  if (!frame_filter_level) return;
#endif
  start_mi_row = 0;
  mi_rows_to_filter = cm->mi_rows;
  if (partial_frame && cm->mi_rows > 8) {
    start_mi_row = cm->mi_rows >> 1;
    start_mi_row &= 0xfffffff8;
    mi_rows_to_filter = AOMMAX(cm->mi_rows / 8, 8);
  }
  end_mi_row = start_mi_row + mi_rows_to_filter;
#if CONFIG_LOOPFILTER_LEVEL
  // TODO(chengchen): refactor the code such that y_only has its matching
  // meaning. Now it means the plane to be filtered in this experiment.
  av1_loop_filter_frame_init(cm, frame_filter_level, frame_filter_level_r,
                             y_only);
#else
  av1_loop_filter_frame_init(cm, frame_filter_level, frame_filter_level);
#endif

#if CONFIG_EXT_DELTA_Q
#if CONFIG_LOOPFILTER_LEVEL
  cm->lf.filter_level[0] = frame_filter_level;
  cm->lf.filter_level[1] = frame_filter_level_r;
#else
  cm->lf.filter_level = frame_filter_level;
#endif
#endif

  av1_loop_filter_rows(frame, cm, xd->plane, start_mi_row, end_mi_row, y_only);

#if CONFIG_EXT_DELTA_Q
#if CONFIG_LOOPFILTER_LEVEL
  cm->lf.filter_level[0] = orig_filter_level[0];
  cm->lf.filter_level[1] = orig_filter_level[1];
#else
  cm->lf.filter_level = orig_filter_level;
#endif
#endif
}

void av1_loop_filter_data_reset(LFWorkerData *lf_data,
                                YV12_BUFFER_CONFIG *frame_buffer,
                                struct AV1Common *cm,
                                const struct macroblockd_plane *planes) {
  lf_data->frame_buffer = frame_buffer;
  lf_data->cm = cm;
  lf_data->start = 0;
  lf_data->stop = 0;
  lf_data->y_only = 0;
  memcpy(lf_data->planes, planes, sizeof(lf_data->planes));
}

int av1_loop_filter_worker(LFWorkerData *const lf_data, void *unused) {
  (void)unused;
  av1_loop_filter_rows(lf_data->frame_buffer, lf_data->cm, lf_data->planes,
                       lf_data->start, lf_data->stop, lf_data->y_only);
  return 1;
}
