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

static void lf_init_lut(loop_filter_info_n *lfi) {
  lfi->mode_lf_lut[DC_PRED] = 0;
  lfi->mode_lf_lut[D45_PRED] = 0;
  lfi->mode_lf_lut[D135_PRED] = 0;
  lfi->mode_lf_lut[D117_PRED] = 0;
  lfi->mode_lf_lut[D153_PRED] = 0;
  lfi->mode_lf_lut[D27_PRED] = 0;
  lfi->mode_lf_lut[D63_PRED] = 0;
  lfi->mode_lf_lut[V_PRED] = 0;
  lfi->mode_lf_lut[H_PRED] = 0;
  lfi->mode_lf_lut[TM_PRED] = 0;
  lfi->mode_lf_lut[ZEROMV]  = 0;
  lfi->mode_lf_lut[NEARESTMV] = 1;
  lfi->mode_lf_lut[NEARMV] = 1;
  lfi->mode_lf_lut[NEWMV] = 1;
}

void vp9_loop_filter_update_sharpness(loop_filter_info_n *lfi,
                                      int sharpness_lvl) {
  int i;

  /* For each possible value for the loop filter fill out limits */
  for (i = 0; i <= MAX_LOOP_FILTER; i++) {
    int filt_lvl = i;
    int block_inside_limit = 0;

    /* Set loop filter paramaeters that control sharpness. */
    block_inside_limit = filt_lvl >> (sharpness_lvl > 0);
    block_inside_limit = block_inside_limit >> (sharpness_lvl > 4);

    if (sharpness_lvl > 0) {
      if (block_inside_limit > (9 - sharpness_lvl))
        block_inside_limit = (9 - sharpness_lvl);
    }

    if (block_inside_limit < 1)
      block_inside_limit = 1;

    vpx_memset(lfi->lim[i], block_inside_limit, SIMD_WIDTH);
    vpx_memset(lfi->blim[i], (2 * filt_lvl + block_inside_limit),
               SIMD_WIDTH);
    vpx_memset(lfi->mblim[i], (2 * (filt_lvl + 2) + block_inside_limit),
               SIMD_WIDTH);
  }
}

void vp9_loop_filter_init(VP9_COMMON *cm) {
  loop_filter_info_n *lfi = &cm->lf_info;
  int i;

  // init limits for given sharpness
  vp9_loop_filter_update_sharpness(lfi, cm->sharpness_level);
  cm->last_sharpness_level = cm->sharpness_level;

  // init LUT for lvl  and hev thr picking
  lf_init_lut(lfi);

  // init hev threshold const vectors
  for (i = 0; i < 4; i++)
    vpx_memset(lfi->hev_thr[i], i, SIMD_WIDTH);
}

void vp9_loop_filter_frame_init(VP9_COMMON *cm,
                                MACROBLOCKD *xd,
                                int default_filt_lvl) {
  int seg,    // segment number
      ref,    // index in ref_lf_deltas
      mode;   // index in mode_lf_deltas
  // n_shift is the a multiplier for lf_deltas
  // the multiplier is 1 for when filter_lvl is between 0 and 31;
  // 2 when filter_lvl is between 32 and 63
  int n_shift = default_filt_lvl >> 5;

  loop_filter_info_n *lfi = &cm->lf_info;

  /* update limits if sharpness has changed */
  // printf("vp9_loop_filter_frame_init %d\n", default_filt_lvl);
  // printf("sharpness level: %d [%d]\n",
  //        cm->sharpness_level, cm->last_sharpness_level);
  if (cm->last_sharpness_level != cm->sharpness_level) {
    vp9_loop_filter_update_sharpness(lfi, cm->sharpness_level);
    cm->last_sharpness_level = cm->sharpness_level;
  }

  for (seg = 0; seg < MAX_MB_SEGMENTS; seg++) {
    int lvl_seg = default_filt_lvl;
    int lvl_ref, lvl_mode;


    // Set the baseline filter values for each segment
    if (vp9_segfeature_active(xd, seg, SEG_LVL_ALT_LF)) {
      /* Abs value */
      if (xd->mb_segment_abs_delta == SEGMENT_ABSDATA) {
        lvl_seg = vp9_get_segdata(xd, seg, SEG_LVL_ALT_LF);
      } else { /* Delta Value */
        lvl_seg += vp9_get_segdata(xd, seg, SEG_LVL_ALT_LF);
        lvl_seg = clamp(lvl_seg, 0, 63);
      }
    }

    if (!xd->mode_ref_lf_delta_enabled) {
      /* we could get rid of this if we assume that deltas are set to
       * zero when not in use; encoder always uses deltas
       */
      vpx_memset(lfi->lvl[seg][0], lvl_seg, 4 * 4);
      continue;
    }

    lvl_ref = lvl_seg;

    /* INTRA_FRAME */
    ref = INTRA_FRAME;

    /* Apply delta for reference frame */
    lvl_ref += xd->ref_lf_deltas[ref] << n_shift;

    mode = 0; /* all the rest of Intra modes */
    lvl_mode = lvl_ref;
    lfi->lvl[seg][ref][mode] = clamp(lvl_mode, 0, 63);

    /* LAST, GOLDEN, ALT */
    for (ref = 1; ref < MAX_REF_FRAMES; ref++) {
      int lvl_ref = lvl_seg;

      /* Apply delta for reference frame */
      lvl_ref += xd->ref_lf_deltas[ref] << n_shift;

      /* Apply delta for Inter modes */
      for (mode = 0; mode < MAX_MODE_LF_DELTAS; mode++) {
        lvl_mode = lvl_ref + (xd->mode_lf_deltas[mode] << n_shift);
        lfi->lvl[seg][ref][mode] = clamp(lvl_mode, 0, 63);
      }
    }
  }
}

static int build_lfi(const VP9_COMMON *cm, const MB_MODE_INFO *mbmi,
                      struct loop_filter_info *lfi) {
  const loop_filter_info_n *lfi_n = &cm->lf_info;
  int mode = mbmi->mode;
  int mode_index = lfi_n->mode_lf_lut[mode];
  int seg = mbmi->segment_id;
  int ref_frame = mbmi->ref_frame[0];
  int filter_level = lfi_n->lvl[seg][ref_frame][mode_index];

  if (filter_level) {
    const int hev_index = filter_level >> 4;
    lfi->mblim = lfi_n->mblim[filter_level];
    lfi->blim = lfi_n->blim[filter_level];
    lfi->lim = lfi_n->lim[filter_level];
    lfi->hev_thr = lfi_n->hev_thr[hev_index];
    return 1;
  }
  return 0;
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

  for (mask = mask_16x16 | mask_8x8 | mask_4x4 | mask_4x4_int;
       mask; mask >>= 1) {
    if (mask & 1) {
      if (!only_4x4_1) {
        if (mask_16x16 & 1) {
          vp9_mb_lpf_horizontal_edge_w(s, pitch, lfi->mblim, lfi->lim,
                                       lfi->hev_thr);
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
    s += 8;
    lfi++;
    mask_16x16 >>= 1;
    mask_8x8 >>= 1;
    mask_4x4 >>= 1;
    mask_4x4_int >>= 1;
  }
}

static void filter_block_plane(VP9_COMMON *cm, MACROBLOCKD *xd,
                               int plane, int mi_row, int mi_col) {
  const int ss_x = xd->plane[plane].subsampling_x;
  const int ss_y = xd->plane[plane].subsampling_y;
  const int row_step = 1 << xd->plane[plane].subsampling_y;
  const int col_step = 1 << xd->plane[plane].subsampling_x;
  struct buf_2d * const dst = &xd->plane[plane].dst;
  uint8_t* const dst0 = dst->buf;
  MODE_INFO* const mi0 = xd->mode_info_context;
  unsigned int mask_16x16[64 / MI_SIZE] = {0};
  unsigned int mask_8x8[64 / MI_SIZE] = {0};
  unsigned int mask_4x4[64 / MI_SIZE] = {0};
  unsigned int mask_4x4_int[64 / MI_SIZE] = {0};
  struct loop_filter_info lfi[64 / MI_SIZE][64 / MI_SIZE];
  int r, c;

  for (r = 0; r < 64 / MI_SIZE && mi_row + r < cm->mi_rows; r += row_step) {
    unsigned int mask_16x16_c = 0;
    unsigned int mask_8x8_c = 0;
    unsigned int mask_4x4_c = 0;
    unsigned int border_mask;

    // Determine the vertical edges that need filtering
    for (c = 0; c < 64 / MI_SIZE && mi_col + c < cm->mi_cols; c += col_step) {
      const MODE_INFO * const mi = xd->mode_info_context;
      const int skip_this = mi[c].mbmi.mb_skip_coeff
                            && mi[c].mbmi.ref_frame[0] != INTRA_FRAME;
      // left edge of current unit is block/partition edge -> no skip
      const int block_edge_left = b_width_log2(mi[c].mbmi.sb_type) ?
          !(c & ((1 << (b_width_log2(mi[c].mbmi.sb_type)-1)) - 1)) : 1;
      const int skip_this_c = skip_this && !block_edge_left;
      // top edge of current unit is block/partition edge -> no skip
      const int block_edge_above = b_height_log2(mi[c].mbmi.sb_type) ?
          !(r & ((1 << (b_height_log2(mi[c].mbmi.sb_type)-1)) - 1)) : 1;
      const int skip_this_r = skip_this && !block_edge_above;
      const TX_SIZE tx_size = plane ? get_uv_tx_size(&mi[c].mbmi)
                                    : mi[c].mbmi.txfm_size;
      const int skip_border_4x4_c = ss_x && mi_col + c == cm->mi_cols - 1;
      const int skip_border_4x4_r = ss_y && mi_row + r == cm->mi_rows - 1;

      // Filter level can vary per MI
      if (!build_lfi(cm, &mi[c].mbmi,
                     lfi[r] + (c >> xd->plane[plane].subsampling_x)))
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
    xd->mode_info_context += cm->mode_info_stride * row_step;
  }

  // Now do horizontal pass
  dst->buf = dst0;
  xd->mode_info_context = mi0;
  for (r = 0; r < 64 / MI_SIZE && mi_row + r < cm->mi_rows; r += row_step) {
    const int skip_border_4x4_r = ss_y && mi_row + r == cm->mi_rows - 1;
    const unsigned int mask_4x4_int_r = skip_border_4x4_r ? 0 : mask_4x4_int[r];

    filter_selectively_horiz(dst->buf, dst->stride,
                             mask_16x16[r],
                             mask_8x8[r],
                             mask_4x4[r],
                             mask_4x4_int_r, mi_row + r == 0, lfi[r]);
    dst->buf += 8 * dst->stride;
    xd->mode_info_context += cm->mode_info_stride * row_step;
  }
}

void vp9_loop_filter_frame(VP9_COMMON *cm,
                           MACROBLOCKD *xd,
                           int frame_filter_level,
                           int y_only) {
  int mi_row, mi_col;

  // Initialize the loop filter for this frame.
  vp9_loop_filter_frame_init(cm, xd, frame_filter_level);

  for (mi_row = 0; mi_row < cm->mi_rows; mi_row += 64 / MI_SIZE) {
    MODE_INFO* const mi = cm->mi + mi_row * cm->mode_info_stride;

    for (mi_col = 0; mi_col < cm->mi_cols; mi_col += 64 / MI_SIZE) {
      int plane;

      setup_dst_planes(xd, cm->frame_to_show, mi_row, mi_col);
      for (plane = 0; plane < (y_only ? 1 : MAX_MB_PLANE); plane++) {
        xd->mode_info_context = mi + mi_col;
        filter_block_plane(cm, xd, plane, mi_row, mi_col);
      }
    }
  }
}
