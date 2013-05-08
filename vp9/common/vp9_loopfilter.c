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
#include "vpx_mem/vpx_mem.h"

#include "vp9/common/vp9_seg_common.h"

static void lf_init_lut(loop_filter_info_n *lfi) {
  lfi->mode_lf_lut[DC_PRED] = 1;
  lfi->mode_lf_lut[D45_PRED] = 1;
  lfi->mode_lf_lut[D135_PRED] = 1;
  lfi->mode_lf_lut[D117_PRED] = 1;
  lfi->mode_lf_lut[D153_PRED] = 1;
  lfi->mode_lf_lut[D27_PRED] = 1;
  lfi->mode_lf_lut[D63_PRED] = 1;
  lfi->mode_lf_lut[V_PRED] = 1;
  lfi->mode_lf_lut[H_PRED] = 1;
  lfi->mode_lf_lut[TM_PRED] = 1;
  lfi->mode_lf_lut[I4X4_PRED]  = 0;
  lfi->mode_lf_lut[ZEROMV]  = 1;
  lfi->mode_lf_lut[NEARESTMV] = 2;
  lfi->mode_lf_lut[NEARMV] = 2;
  lfi->mode_lf_lut[NEWMV] = 2;
  lfi->mode_lf_lut[SPLITMV] = 3;
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

    /* Apply delta for Intra modes */
    mode = 0; /* I4X4_PRED */
    /* Only the split mode I4X4_PRED has a further special case */
    lvl_mode = lvl_ref + (xd->mode_lf_deltas[mode] << n_shift);
    lfi->lvl[seg][ref][mode] = clamp(lvl_mode, 0, 63);

    mode = 1; /* all the rest of Intra modes */
    lvl_mode = lvl_ref;
    lfi->lvl[seg][ref][mode] = clamp(lvl_mode, 0, 63);

    /* LAST, GOLDEN, ALT */
    for (ref = 1; ref < MAX_REF_FRAMES; ref++) {
      int lvl_ref = lvl_seg;

      /* Apply delta for reference frame */
      lvl_ref += xd->ref_lf_deltas[ref] << n_shift;

      /* Apply delta for Inter modes */
      for (mode = 1; mode < 4; mode++) {
        lvl_mode = lvl_ref + (xd->mode_lf_deltas[mode] << n_shift);
        lfi->lvl[seg][ref][mode] = clamp(lvl_mode, 0, 63);
      }
    }
  }
}

// Determine if we should skip inner-MB loop filtering within a MB
// The current condition is that the loop filtering is skipped only
// the MB uses a prediction size of 16x16 and either 16x16 transform
// is used or there is no residue at all.
static int mb_lf_skip(const MB_MODE_INFO *const mbmi) {
  const int skip_coef = mbmi->mb_skip_coeff;
  const int tx_size = mbmi->txfm_size;
  return mbmi->sb_type >= BLOCK_SIZE_MB16X16 &&
         (tx_size >= TX_16X16 || skip_coef);
}

// Determine if we should skip MB loop filtering on a MB edge within
// a superblock, the current condition is that MB loop filtering is
// skipped only when both MBs do not use inner MB loop filtering, and
// same motion vector with same reference frame
static int sb_mb_lf_skip(const MODE_INFO *const mip0,
                         const MODE_INFO *const mip1) {
  const MB_MODE_INFO *mbmi0 = &mip0->mbmi;
  const MB_MODE_INFO *mbmi1 = &mip1->mbmi;
  return mb_lf_skip(mbmi0) && mb_lf_skip(mbmi1) &&
         mbmi0->ref_frame != INTRA_FRAME &&
         mbmi1->ref_frame != INTRA_FRAME;
}

static void lpf_mb(VP9_COMMON *cm, const MODE_INFO *mi,
                   int do_left_mb_v, int do_above_mb_h,
                   int do_left_mbuv_v, int do_above_mbuv_h,
                   uint8_t *y_ptr, uint8_t *u_ptr, uint8_t *v_ptr,
                   int y_stride, int uv_stride, int dering) {
  loop_filter_info_n *lfi_n = &cm->lf_info;
  struct loop_filter_info lfi;
  int mode = mi->mbmi.mode;
  int mode_index = lfi_n->mode_lf_lut[mode];
  int seg = mi->mbmi.segment_id;
  int ref_frame = mi->mbmi.ref_frame;
  int filter_level = lfi_n->lvl[seg][ref_frame][mode_index];

  if (filter_level) {
    const int skip_lf = mb_lf_skip(&mi->mbmi);
    const int tx_size = mi->mbmi.txfm_size;
    if (cm->filter_type == NORMAL_LOOPFILTER) {
      const int hev_index = filter_level >> 4;
      lfi.mblim = lfi_n->mblim[filter_level];
      lfi.blim = lfi_n->blim[filter_level];
      lfi.lim = lfi_n->lim[filter_level];
      lfi.hev_thr = lfi_n->hev_thr[hev_index];

      if (do_above_mb_h) {
        if (tx_size >= TX_16X16)
          vp9_lpf_mbh_w(y_ptr,
                        do_above_mbuv_h ? u_ptr : NULL,
                        do_above_mbuv_h ? v_ptr : NULL,
                        y_stride, uv_stride, &lfi);
        else
          vp9_loop_filter_mbh(y_ptr, u_ptr, v_ptr, y_stride, uv_stride, &lfi);
      }

      if (!skip_lf) {
        if (tx_size >= TX_8X8) {
          if (tx_size == TX_8X8 &&
              (mi->mbmi.sb_type < BLOCK_SIZE_MB16X16)
              )
            vp9_loop_filter_bh8x8(y_ptr, u_ptr, v_ptr,
                                  y_stride, uv_stride, &lfi);
          else
            vp9_loop_filter_bh8x8(y_ptr, NULL, NULL,
                                  y_stride, uv_stride, &lfi);
        } else {
          vp9_loop_filter_bh(y_ptr, u_ptr, v_ptr,
                             y_stride, uv_stride, &lfi);
        }
      }

      if (do_left_mb_v) {
        if (tx_size >= TX_16X16)
          vp9_lpf_mbv_w(y_ptr,
                        do_left_mbuv_v ? u_ptr : NULL,
                        do_left_mbuv_v ? v_ptr : NULL,
                        y_stride, uv_stride, &lfi);
        else
          vp9_loop_filter_mbv(y_ptr, u_ptr, v_ptr, y_stride, uv_stride, &lfi);
      }

      if (!skip_lf) {
        if (tx_size >= TX_8X8) {
          if (tx_size == TX_8X8 &&
              (mi->mbmi.sb_type < BLOCK_SIZE_MB16X16))
            vp9_loop_filter_bv8x8(y_ptr, u_ptr, v_ptr,
                                  y_stride, uv_stride, &lfi);
          else
            vp9_loop_filter_bv8x8(y_ptr, NULL, NULL,
                                  y_stride, uv_stride, &lfi);
        } else {
          vp9_loop_filter_bv(y_ptr, u_ptr, v_ptr,
                             y_stride, uv_stride, &lfi);
        }
      }
      if (dering) {
#if CONFIG_LOOP_DERING
        vp9_post_proc_down_and_across(y_ptr, y_ptr,
          y_stride, y_stride,
          16, 16, dering);
        if (u_ptr && v_ptr) {
          vp9_post_proc_down_and_across(u_ptr, u_ptr,
            uv_stride, uv_stride,
            8, 8, dering);
          vp9_post_proc_down_and_across(v_ptr, v_ptr,
            uv_stride, uv_stride,
            8, 8, dering);
        }
#endif
      }
    } else {
      // TODO(yaowu): simple loop filter
    }
  }
}

static void lpf_sb32(VP9_COMMON *cm, const MODE_INFO *mode_info_context,
                     int mb_row, int mb_col,
                     uint8_t *y_ptr, uint8_t *u_ptr, uint8_t *v_ptr,
                     int y_stride, int uv_stride,
                     int y_only, int dering) {
  BLOCK_SIZE_TYPE sb_type = mode_info_context->mbmi.sb_type;
  const int wbl = b_width_log2(sb_type), hbl = b_height_log2(sb_type);
  TX_SIZE tx_size = mode_info_context->mbmi.txfm_size;
  int do_left_v, do_above_h;
  int do_left_v_mbuv, do_above_h_mbuv;
  int mis = cm->mode_info_stride;
  const MODE_INFO *mi;

  // process 1st MB top-left
  mi = mode_info_context;
  do_left_v = (mb_col > 0);
  do_above_h = (mb_row > 0);
  do_left_v_mbuv = !(sb_type >= BLOCK_SIZE_SB64X64 &&
      tx_size >= TX_32X32 && (mb_col & 2));
  do_above_h_mbuv = !(sb_type >= BLOCK_SIZE_SB64X64 &&
      tx_size >= TX_32X32 && (mb_row & 2));
  lpf_mb(cm, mi, do_left_v, do_above_h,
      do_left_v_mbuv, do_above_h_mbuv,
      y_ptr,
      y_only? 0 : u_ptr,
      y_only? 0 : v_ptr,
      y_stride, uv_stride, dering);
  // process 2nd MB top-right
  mi = mode_info_context + 2;
  do_left_v = !(wbl >= 3 /* 32x16 or >=32x32 */ && (tx_size >= TX_32X32 ||
      sb_mb_lf_skip(mode_info_context, mi)));
  do_above_h = (mb_row > 0);
  do_left_v_mbuv = !(wbl >= 3 /* 32x16 or >=32x32 */ && (tx_size >= TX_16X16 ||
      sb_mb_lf_skip(mode_info_context, mi)));
  do_above_h_mbuv = !(sb_type >= BLOCK_SIZE_SB64X64 &&
      tx_size >= TX_32X32 && (mb_row & 2));
  lpf_mb(cm, mi, do_left_v, do_above_h,
      do_left_v_mbuv, do_above_h_mbuv,
      y_ptr + 16,
      y_only ? 0 : (u_ptr + 8),
      y_only ? 0 : (v_ptr + 8),
      y_stride, uv_stride, dering);

  // process 3rd MB bottom-left
  mi = mode_info_context + (mis << 1);
  do_left_v = (mb_col > 0);
  do_above_h = !(hbl >= 3 /* 16x32 or >=32x32 */ && (tx_size >= TX_32X32 ||
      sb_mb_lf_skip(mode_info_context, mi)));
  do_left_v_mbuv = !(sb_type >= BLOCK_SIZE_SB64X64 &&
      tx_size >= TX_32X32 && (mb_col & 2));
  do_above_h_mbuv = !(hbl >= 3 /* 16x32 or >=32x32 */ && (tx_size >= TX_16X16 ||
      sb_mb_lf_skip(mode_info_context, mi)));
  lpf_mb(cm, mi, do_left_v, do_above_h,
      do_left_v_mbuv, do_above_h_mbuv,
      y_ptr + 16 * y_stride,
      y_only ? 0 : (u_ptr + 8 * uv_stride),
      y_only ? 0 : (v_ptr + 8 * uv_stride),
      y_stride, uv_stride, dering);

  // process 4th MB bottom right
  mi = mode_info_context + ((mis + 1) << 1);
  do_left_v = !(wbl >= 3 /* 32x16 or >=32x32 */ && (tx_size >= TX_32X32 ||
      sb_mb_lf_skip(mi - 2, mi)));
  do_above_h = !(hbl >= 3 /* 16x32 or >=32x32 */ && (tx_size >= TX_32X32 ||
      sb_mb_lf_skip(mode_info_context + 2, mi)));
  do_left_v_mbuv = (wbl >= 3 /* 32x16 or >=32x32 */ && (tx_size >= TX_16X16 ||
      sb_mb_lf_skip(mi - 2, mi)));
  do_above_h_mbuv = !(hbl >= 3 /* 16x32 or >=32x32 */ && (tx_size >= TX_16X16 ||
      sb_mb_lf_skip(mode_info_context + 2, mi)));
  lpf_mb(cm, mi, do_left_v, do_above_h,
      do_left_v_mbuv, do_above_h_mbuv,
      y_ptr + 16 * y_stride + 16,
      y_only ? 0 : (u_ptr + 8 * uv_stride + 8),
      y_only ? 0 : (v_ptr + 8 * uv_stride + 8),
      y_stride, uv_stride, dering);
}

static void lpf_sb64(VP9_COMMON *cm, const MODE_INFO *mode_info_context,
                     int mb_row, int mb_col,
                     uint8_t *y_ptr, uint8_t *u_ptr, uint8_t *v_ptr,
                     int y_stride, int uv_stride,
                     int y_only, int dering) {
  lpf_sb32(cm, mode_info_context, mb_row, mb_col,
      y_ptr, u_ptr, v_ptr,
      y_stride, uv_stride, y_only, dering);
  lpf_sb32(cm, mode_info_context + 4, mb_row, mb_col + 2,
      y_ptr + 32, u_ptr + 16, v_ptr + 16,
      y_stride, uv_stride, y_only, dering);
  lpf_sb32(cm, mode_info_context + cm->mode_info_stride * 4,
      mb_row + 2, mb_col,
      y_ptr + 32 * y_stride,
      u_ptr + 16 * uv_stride,
      v_ptr + 16 * uv_stride,
      y_stride, uv_stride, y_only, dering);
  lpf_sb32(cm, mode_info_context + cm->mode_info_stride * 4 + 4,
      mb_row + 2, mb_col + 2,
      y_ptr + 32 * y_stride + 32,
      u_ptr + 16 * uv_stride + 16,
      v_ptr + 16 * uv_stride + 16,
      y_stride, uv_stride, y_only, dering);
}

void vp9_loop_filter_frame(VP9_COMMON *cm,
                           MACROBLOCKD *xd,
                           int frame_filter_level,
                           int y_only,
                           int dering) {
  YV12_BUFFER_CONFIG *post = cm->frame_to_show;
  int mb_row, mb_col;
  const int sb64_rows = cm->mb_rows / 4;
  const int sb64_cols = cm->mb_cols / 4;
  const int extra_sb32_row = (cm->mb_rows & 2) != 0;
  const int extra_sb32_col = (cm->mb_cols & 2) != 0;
  const int extra_mb_col = cm->mb_cols & 1;
  const int extra_mb_row = cm->mb_rows & 1;
  // Set up the buffer pointers
  uint8_t *y_ptr = post->y_buffer;
  uint8_t *u_ptr = y_only ? 0 : post->u_buffer;
  uint8_t *v_ptr = y_only ? 0 : post->v_buffer;

  // Point at base of Mb MODE_INFO list
  const MODE_INFO *mode_info_context = cm->mi;
  const MODE_INFO *mi;
  const int mis = cm->mode_info_stride;
  const int y_stride = post->y_stride;
  const int uv_stride = post->uv_stride;
  // These two flags signal if MB left edge and above edge
  // should be filtered using MB edge filter. Currently, MB
  // edge filtering is not applied on MB edge internal to a
  // 32x32 superblock if:
  // 1) SB32 is using 32x32 prediction and 32x32 transform
  // 2) SB32 is using 32x32 prediction and 16x16 transform
  //    but all coefficients are zero.
  // MB edges are on 32x32 superblock boundary are always
  // filtered except on image frame boundary.
  int do_left_v, do_above_h;
  // These two flags signal if MB UV left edge and above edge
  // should be filtered using MB edge filter. Currently, MB
  // edge filtering is not applied for MB edges internal to
  // a 32x32 superblock if:
  // 1) SB32 is using 32x32 prediction and 32x32 transform
  // 2) SB32 is using 32x32 prediction and 16x16 transform
  //    but all coefficients are zero.
  // 3) SB32 UV edges internal to a SB64 and 32x32 transform
  //    is used, i.e. UV is doing 32x32 transform hence no
  //    transform boundary exists inside the SB64 for UV
  int do_left_v_mbuv, do_above_h_mbuv;

  // Initialize the loop filter for this frame.
  vp9_loop_filter_frame_init(cm, xd, frame_filter_level);

  // vp9_filter each 64x64 SB
  // For each SB64: the 4 SB32 are filtered in raster scan order
  // For each SB32: the 4 MBs are filtered in raster scan order
  // For each MB:   the left and above MB edges as well as the
  //                internal block edges are processed together
  for (mb_row = 0; mb_row < sb64_rows * 4; mb_row += 4) {
    for (mb_col = 0; mb_col < sb64_cols * 4; mb_col += 4) {
      lpf_sb64(cm, mode_info_context, mb_row, mb_col,
               y_ptr, u_ptr, v_ptr,
               y_stride, uv_stride, y_only, dering);
      y_ptr += 64;
      u_ptr = y_only? 0 : u_ptr + 32;
      v_ptr = y_only? 0 : v_ptr + 32;
      mode_info_context += 8;       // step to next SB64
    }
    if (extra_sb32_col) {
      // process 2 SB32s in the extra SB32 col
      lpf_sb32(cm, mode_info_context, mb_row, mb_col,
               y_ptr, u_ptr, v_ptr,
               y_stride, uv_stride, y_only, dering);
      lpf_sb32(cm, mode_info_context + mis * 4,
               mb_row + 2, mb_col,
               y_ptr + 32 * y_stride,
               u_ptr + 16 * uv_stride,
               v_ptr + 16 * uv_stride,
               y_stride, uv_stride, y_only, dering);
      y_ptr += 32;
      u_ptr = y_only? 0 : u_ptr + 16;
      v_ptr = y_only? 0 : v_ptr + 16;
      mode_info_context += 4;       // step to next SB32
      mb_col += 2;
    }
    if (extra_mb_col) {
      // process 4 MB in the extra MB col
      int k;
      for (k = 0; k < 4; ++k) {
        mi = mode_info_context + (mis << 1) * k;
        do_left_v = (mb_col > 0);
        do_above_h = k == 0 ? mb_row > 0 : 1;
        do_left_v_mbuv =  1;
        do_above_h_mbuv = 1;
        lpf_mb(cm, mi, do_left_v, do_above_h,
               do_left_v_mbuv, do_above_h_mbuv,
               y_ptr + (k * 16) * y_stride,
               y_only ? 0 : (u_ptr + (k * 8) * uv_stride),
               y_only ? 0 : (v_ptr + (k * 8) * uv_stride),
               y_stride, uv_stride, dering);
      }

      y_ptr += 16;
      u_ptr = y_only? 0 : u_ptr + 8;
      v_ptr = y_only? 0 : v_ptr + 8;
      mode_info_context += 2;       // step to next MB
    }
    // move pointers to the begining of next sb64 row
    y_ptr += y_stride  * 64 - post->y_width;
    if (!y_only) {
      u_ptr += uv_stride *  32 - post->uv_width;
      v_ptr += uv_stride *  32 - post->uv_width;
    }
    /* skip to next SB64 row */
    mode_info_context += mis * 8 - cm->mi_cols;
  }
  if (extra_sb32_row) {
    const int sb32_cols = sb64_cols * 2 + extra_sb32_col;
    for (mb_col = 0; mb_col < sb32_cols * 2; mb_col += 2) {
      lpf_sb32(cm, mode_info_context, mb_row, mb_col,
               y_ptr, u_ptr, v_ptr,
               y_stride, uv_stride, y_only, dering);
      y_ptr += 32;
      u_ptr = y_only? 0 : u_ptr + 16;
      v_ptr = y_only? 0 : v_ptr + 16;
      mode_info_context += 4;       // step to next SB32
    }
    if (extra_mb_col) {
      // process 1st MB
      mi = mode_info_context;
      do_left_v = (mb_col > 0);
      do_above_h = (mb_row > 0);
      do_left_v_mbuv =  1;
      do_above_h_mbuv = 1;
      lpf_mb(cm, mi, do_left_v, do_above_h,
             do_left_v_mbuv, do_above_h_mbuv,
             y_ptr,
             y_only? NULL : u_ptr,
             y_only? NULL : v_ptr,
             y_stride, uv_stride, dering);
      // process 2nd MB
      mi = mode_info_context + (mis << 1);
      do_left_v = (mb_col > 0);
      do_above_h = 1;
      do_left_v_mbuv =  1;
      do_above_h_mbuv = 1;
      lpf_mb(cm, mi, do_left_v, do_above_h,
             do_left_v_mbuv, do_above_h_mbuv,
             y_ptr + 16 * y_stride,
             y_only ? NULL : (u_ptr + 8 * uv_stride),
             y_only ? NULL : (v_ptr + 8 * uv_stride),
             y_stride, uv_stride, dering);
      y_ptr += 16;
      u_ptr = y_only? 0 : u_ptr + 8;
      v_ptr = y_only? 0 : v_ptr + 8;
      mode_info_context += 2;       /* step to next MB */
    }
    // move pointers to the beginning of next sb64 row
    y_ptr += y_stride * 32 - post->y_width;
    u_ptr += y_only? 0 : uv_stride *  16 - post->uv_width;
    v_ptr += y_only? 0 : uv_stride *  16 - post->uv_width;
    // skip to next MB row if exist
    mode_info_context += mis * 4 - cm->mi_cols;
    mb_row += 2;
  }
  if (extra_mb_row) {
    for (mb_col = 0; mb_col < cm->mb_cols; mb_col++) {
      const MODE_INFO *mi = mode_info_context;
      do_left_v =  (mb_col > 0);
      do_above_h = (mb_row > 0);
      do_left_v_mbuv = 1;
      do_above_h_mbuv = 1;
      lpf_mb(cm, mi, do_left_v, do_above_h,
             do_left_v_mbuv, do_above_h_mbuv,
             y_ptr,
             y_only? 0 : u_ptr,
             y_only? 0 : v_ptr,
             y_stride, uv_stride, dering);
      y_ptr += 16;
      u_ptr = y_only? 0 : u_ptr + 8;
      v_ptr = y_only? 0 : v_ptr + 8;
      mode_info_context += 2;     // step to next MB
    }
  }
}
