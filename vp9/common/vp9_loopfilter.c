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

static void update_sharpness(loop_filter_info_n *const lfi, int sharpness_lvl) {
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

void vp9_loop_filter_frame_init(VP9_COMMON *const cm, MACROBLOCKD *const xd,
                                int default_filt_lvl) {
  int seg_id;
  // n_shift is the a multiplier for lf_deltas
  // the multiplier is 1 for when filter_lvl is between 0 and 31;
  // 2 when filter_lvl is between 32 and 63
  const int n_shift = default_filt_lvl >> 5;
  loop_filter_info_n *const lfi = &cm->lf_info;
  struct loopfilter *const lf = &cm->lf;
  struct segmentation *const seg = &xd->seg;

  // update limits if sharpness has changed
  if (lf->last_sharpness_level != lf->sharpness_level) {
    update_sharpness(lfi, lf->sharpness_level);
    lf->last_sharpness_level = lf->sharpness_level;
  }

  for (seg_id = 0; seg_id < MAX_SEGMENTS; seg_id++) {
    int lvl_seg = default_filt_lvl, ref, mode, intra_lvl;

    // Set the baseline filter values for each segment
    if (vp9_segfeature_active(&xd->seg, seg_id, SEG_LVL_ALT_LF)) {
      const int data = vp9_get_segdata(seg, seg_id, SEG_LVL_ALT_LF);
      lvl_seg = seg->abs_delta == SEGMENT_ABSDATA
                  ? data
                  : clamp(default_filt_lvl + data, 0, MAX_LOOP_FILTER);
    }

    if (!lf->mode_ref_delta_enabled) {
      // we could get rid of this if we assume that deltas are set to
      // zero when not in use; encoder always uses deltas
      vpx_memset(lfi->lvl[seg_id][0], lvl_seg, sizeof(lfi->lvl[seg_id][0]));
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

static int build_lfi(const loop_filter_info_n *const lfi_n,
                     const MB_MODE_INFO *const mbmi,
                     struct loop_filter_info *const lfi) {
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

static void filter_block_plane(VP9_COMMON *const cm,
                               struct macroblockd_plane *const plane,
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
                            : mi[c].mbmi.txfm_size;
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

void vp9_loop_filter_rows(const YV12_BUFFER_CONFIG *frame_buffer,
                          VP9_COMMON *cm, MACROBLOCKD *xd,
                          int start, int stop, int y_only) {
  const int num_planes = y_only ? 1 : MAX_MB_PLANE;
  int mi_row, mi_col;

  for (mi_row = start; mi_row < stop; mi_row += MI_BLOCK_SIZE) {
    MODE_INFO* const mi = cm->mi + mi_row * cm->mode_info_stride;

    for (mi_col = 0; mi_col < cm->mi_cols; mi_col += MI_BLOCK_SIZE) {
      int plane;

      setup_dst_planes(xd, frame_buffer, mi_row, mi_col);
      for (plane = 0; plane < num_planes; ++plane) {
        filter_block_plane(cm, &xd->plane[plane], mi + mi_col, mi_row, mi_col);
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
  vp9_loop_filter_frame_init(cm, xd, frame_filter_level);
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
