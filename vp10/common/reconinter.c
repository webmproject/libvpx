/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include <assert.h>

#include "./vpx_scale_rtcd.h"
#include "./vpx_config.h"

#include "vpx/vpx_integer.h"

#include "vp10/common/blockd.h"
#include "vp10/common/reconinter.h"
#include "vp10/common/reconintra.h"

#if CONFIG_VP9_HIGHBITDEPTH
void vp10_highbd_build_inter_predictor(const uint8_t *src, int src_stride,
                                      uint8_t *dst, int dst_stride,
                                      const MV *src_mv,
                                      const struct scale_factors *sf,
                                      int w, int h, int ref,
                                      const INTERP_FILTER interp_filter,
                                      enum mv_precision precision,
                                      int x, int y, int bd) {
  const int is_q4 = precision == MV_PRECISION_Q4;
  const MV mv_q4 = { is_q4 ? src_mv->row : src_mv->row * 2,
                     is_q4 ? src_mv->col : src_mv->col * 2 };
  MV32 mv = vp10_scale_mv(&mv_q4, x, y, sf);
  const int subpel_x = mv.col & SUBPEL_MASK;
  const int subpel_y = mv.row & SUBPEL_MASK;

  src += (mv.row >> SUBPEL_BITS) * src_stride + (mv.col >> SUBPEL_BITS);

  high_inter_predictor(src, src_stride, dst, dst_stride, subpel_x, subpel_y,
                       sf, w, h, ref, interp_filter, sf->x_step_q4,
                       sf->y_step_q4, bd);
}
#endif  // CONFIG_VP9_HIGHBITDEPTH

void vp10_build_inter_predictor(const uint8_t *src, int src_stride,
                               uint8_t *dst, int dst_stride,
                               const MV *src_mv,
                               const struct scale_factors *sf,
                               int w, int h, int ref,
                               const INTERP_FILTER interp_filter,
                               enum mv_precision precision,
                               int x, int y) {
  const int is_q4 = precision == MV_PRECISION_Q4;
  const MV mv_q4 = { is_q4 ? src_mv->row : src_mv->row * 2,
                     is_q4 ? src_mv->col : src_mv->col * 2 };
  MV32 mv = vp10_scale_mv(&mv_q4, x, y, sf);
  const int subpel_x = mv.col & SUBPEL_MASK;
  const int subpel_y = mv.row & SUBPEL_MASK;

  src += (mv.row >> SUBPEL_BITS) * src_stride + (mv.col >> SUBPEL_BITS);

  inter_predictor(src, src_stride, dst, dst_stride, subpel_x, subpel_y,
                  sf, w, h, ref, interp_filter, sf->x_step_q4, sf->y_step_q4);
}

void build_inter_predictors(MACROBLOCKD *xd, int plane, int block,
                            int bw, int bh,
                            int x, int y, int w, int h,
                            int mi_x, int mi_y) {
  struct macroblockd_plane *const pd = &xd->plane[plane];
  const MODE_INFO *mi = xd->mi[0];
  const int is_compound = has_second_ref(&mi->mbmi);
  const INTERP_FILTER interp_filter = mi->mbmi.interp_filter;
  int ref;

  for (ref = 0; ref < 1 + is_compound; ++ref) {
    const struct scale_factors *const sf = &xd->block_refs[ref]->sf;
    struct buf_2d *const pre_buf = &pd->pre[ref];
    struct buf_2d *const dst_buf = &pd->dst;
    uint8_t *const dst = dst_buf->buf + dst_buf->stride * y + x;
    const MV mv = mi->mbmi.sb_type < BLOCK_8X8
               ? average_split_mvs(pd, mi, ref, block)
               : mi->mbmi.mv[ref].as_mv;

    // TODO(jkoleszar): This clamping is done in the incorrect place for the
    // scaling case. It needs to be done on the scaled MV, not the pre-scaling
    // MV. Note however that it performs the subsampling aware scaling so
    // that the result is always q4.
    // mv_precision precision is MV_PRECISION_Q4.
    const MV mv_q4 = clamp_mv_to_umv_border_sb(xd, &mv, bw, bh,
                                               pd->subsampling_x,
                                               pd->subsampling_y);

    uint8_t *pre;
    MV32 scaled_mv;
    int xs, ys, subpel_x, subpel_y;
    const int is_scaled = vp10_is_scaled(sf);

    if (is_scaled) {
      pre = pre_buf->buf + scaled_buffer_offset(x, y, pre_buf->stride, sf);
      scaled_mv = vp10_scale_mv(&mv_q4, mi_x + x, mi_y + y, sf);
      xs = sf->x_step_q4;
      ys = sf->y_step_q4;
    } else {
      pre = pre_buf->buf + (y * pre_buf->stride + x);
      scaled_mv.row = mv_q4.row;
      scaled_mv.col = mv_q4.col;
      xs = ys = 16;
    }
    subpel_x = scaled_mv.col & SUBPEL_MASK;
    subpel_y = scaled_mv.row & SUBPEL_MASK;
    pre += (scaled_mv.row >> SUBPEL_BITS) * pre_buf->stride
           + (scaled_mv.col >> SUBPEL_BITS);

#if CONFIG_VP9_HIGHBITDEPTH
    if (xd->cur_buf->flags & YV12_FLAG_HIGHBITDEPTH) {
      high_inter_predictor(pre, pre_buf->stride, dst, dst_buf->stride,
                           subpel_x, subpel_y, sf, w, h, ref,
                           interp_filter, xs, ys, xd->bd);
    } else {
      inter_predictor(pre, pre_buf->stride, dst, dst_buf->stride,
                      subpel_x, subpel_y, sf, w, h, ref, interp_filter, xs, ys);
    }
#else
    inter_predictor(pre, pre_buf->stride, dst, dst_buf->stride,
                    subpel_x, subpel_y, sf, w, h, ref, interp_filter, xs, ys);
#endif  // CONFIG_VP9_HIGHBITDEPTH
  }
}

void vp10_build_inter_predictor_sub8x8(MACROBLOCKD *xd, int plane,
                                       int i, int ir, int ic,
                                       int mi_row, int mi_col) {
  struct macroblockd_plane *const pd = &xd->plane[plane];
  MODE_INFO *const mi = xd->mi[0];
  const BLOCK_SIZE plane_bsize = get_plane_block_size(mi->mbmi.sb_type, pd);
  const int width = 4 * num_4x4_blocks_wide_lookup[plane_bsize];
  const int height = 4 * num_4x4_blocks_high_lookup[plane_bsize];

  uint8_t *const dst = &pd->dst.buf[(ir * pd->dst.stride + ic) << 2];
  int ref;
  const int is_compound = has_second_ref(&mi->mbmi);
  const INTERP_FILTER interp_filter = mi->mbmi.interp_filter;

  for (ref = 0; ref < 1 + is_compound; ++ref) {
    const uint8_t *pre =
        &pd->pre[ref].buf[(ir * pd->pre[ref].stride + ic) << 2];
#if CONFIG_VP9_HIGHBITDEPTH
  if (xd->cur_buf->flags & YV12_FLAG_HIGHBITDEPTH) {
    vp10_highbd_build_inter_predictor(pre, pd->pre[ref].stride,
                                      dst, pd->dst.stride,
                                      &mi->bmi[i].as_mv[ref].as_mv,
                                      &xd->block_refs[ref]->sf, width, height,
                                      ref, interp_filter, MV_PRECISION_Q3,
                                      mi_col * MI_SIZE + 4 * ic,
                                      mi_row * MI_SIZE + 4 * ir, xd->bd);
  } else {
    vp10_build_inter_predictor(pre, pd->pre[ref].stride,
                               dst, pd->dst.stride,
                               &mi->bmi[i].as_mv[ref].as_mv,
                               &xd->block_refs[ref]->sf, width, height, ref,
                               interp_filter, MV_PRECISION_Q3,
                               mi_col * MI_SIZE + 4 * ic,
                               mi_row * MI_SIZE + 4 * ir);
  }
#else
    vp10_build_inter_predictor(pre, pd->pre[ref].stride,
                               dst, pd->dst.stride,
                               &mi->bmi[i].as_mv[ref].as_mv,
                               &xd->block_refs[ref]->sf, width, height, ref,
                               interp_filter, MV_PRECISION_Q3,
                               mi_col * MI_SIZE + 4 * ic,
                               mi_row * MI_SIZE + 4 * ir);
#endif  // CONFIG_VP9_HIGHBITDEPTH
  }
}

static void build_inter_predictors_for_planes(MACROBLOCKD *xd, BLOCK_SIZE bsize,
                                              int mi_row, int mi_col,
                                              int plane_from, int plane_to) {
  int plane;
  const int mi_x = mi_col * MI_SIZE;
  const int mi_y = mi_row * MI_SIZE;
  for (plane = plane_from; plane <= plane_to; ++plane) {
    const struct macroblockd_plane *pd = &xd->plane[plane];
    const int bw = 4 * num_4x4_blocks_wide_lookup[bsize] >> pd->subsampling_x;
    const int bh = 4 * num_4x4_blocks_high_lookup[bsize] >> pd->subsampling_y;

    if (xd->mi[0]->mbmi.sb_type < BLOCK_8X8) {
      const PARTITION_TYPE bp = bsize - xd->mi[0]->mbmi.sb_type;
      const int have_vsplit = bp != PARTITION_HORZ;
      const int have_hsplit = bp != PARTITION_VERT;
      const int num_4x4_w = 2 >> ((!have_vsplit) | pd->subsampling_x);
      const int num_4x4_h = 2 >> ((!have_hsplit) | pd->subsampling_y);
      const int pw = 8 >> (have_vsplit | pd->subsampling_x);
      const int ph = 8 >> (have_hsplit | pd->subsampling_y);
      int x, y;
      assert(bp != PARTITION_NONE && bp < PARTITION_TYPES);
      assert(bsize == BLOCK_8X8);
      assert(pw * num_4x4_w == bw && ph * num_4x4_h == bh);
      for (y = 0; y < num_4x4_h; ++y)
        for (x = 0; x < num_4x4_w; ++x)
           build_inter_predictors(xd, plane, y * 2 + x, bw, bh,
                                  4 * x, 4 * y, pw, ph, mi_x, mi_y);
    } else {
      build_inter_predictors(xd, plane, 0, bw, bh,
                             0, 0, bw, bh, mi_x, mi_y);
    }
  }
}

void vp10_build_inter_predictors_sby(MACROBLOCKD *xd, int mi_row, int mi_col,
                                    BLOCK_SIZE bsize) {
  build_inter_predictors_for_planes(xd, bsize, mi_row, mi_col, 0, 0);
}

void vp10_build_inter_predictors_sbp(MACROBLOCKD *xd, int mi_row, int mi_col,
                                    BLOCK_SIZE bsize, int plane) {
  build_inter_predictors_for_planes(xd, bsize, mi_row, mi_col, plane, plane);
}

void vp10_build_inter_predictors_sbuv(MACROBLOCKD *xd, int mi_row, int mi_col,
                                     BLOCK_SIZE bsize) {
  build_inter_predictors_for_planes(xd, bsize, mi_row, mi_col, 1,
                                    MAX_MB_PLANE - 1);
}

void vp10_build_inter_predictors_sb(MACROBLOCKD *xd, int mi_row, int mi_col,
                                   BLOCK_SIZE bsize) {
  build_inter_predictors_for_planes(xd, bsize, mi_row, mi_col, 0,
                                    MAX_MB_PLANE - 1);
}

void vp10_setup_dst_planes(struct macroblockd_plane planes[MAX_MB_PLANE],
                          const YV12_BUFFER_CONFIG *src,
                          int mi_row, int mi_col) {
  uint8_t *const buffers[MAX_MB_PLANE] = { src->y_buffer, src->u_buffer,
      src->v_buffer};
  const int strides[MAX_MB_PLANE] = { src->y_stride, src->uv_stride,
      src->uv_stride};
  int i;

  for (i = 0; i < MAX_MB_PLANE; ++i) {
    struct macroblockd_plane *const pd = &planes[i];
    setup_pred_plane(&pd->dst, buffers[i], strides[i], mi_row, mi_col, NULL,
                     pd->subsampling_x, pd->subsampling_y);
  }
}

void vp10_setup_pre_planes(MACROBLOCKD *xd, int idx,
                          const YV12_BUFFER_CONFIG *src,
                          int mi_row, int mi_col,
                          const struct scale_factors *sf) {
  if (src != NULL) {
    int i;
    uint8_t *const buffers[MAX_MB_PLANE] = { src->y_buffer, src->u_buffer,
        src->v_buffer};
    const int strides[MAX_MB_PLANE] = { src->y_stride, src->uv_stride,
        src->uv_stride};
    for (i = 0; i < MAX_MB_PLANE; ++i) {
      struct macroblockd_plane *const pd = &xd->plane[i];
      setup_pred_plane(&pd->pre[idx], buffers[i], strides[i], mi_row, mi_col,
                       sf, pd->subsampling_x, pd->subsampling_y);
    }
  }
}

#if CONFIG_SUPERTX
static const uint8_t mask_8[8] = {
  64, 64, 62, 52, 12,  2,  0,  0
};

static const uint8_t mask_16[16] = {
  63, 62, 60, 58, 55, 50, 43, 36, 28, 21, 14, 9, 6, 4, 2, 1
};

static const uint8_t mask_32[32] = {
  64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 63, 61, 57, 52, 45, 36,
  28, 19, 12,  7,  3,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0
};

static const uint8_t mask_8_uv[8] = {
  64, 64, 62, 52,  12,  2,  0,  0
};

static const uint8_t mask_16_uv[16] = {
  64, 64, 64, 64, 61, 53, 45, 36, 28, 19, 11, 3, 0,  0,  0,  0
};

static const uint8_t mask_32_uv[32] = {
  64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 60, 54, 46, 36,
  28, 18, 10,  4,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0
};

static void generate_1dmask(int length, uint8_t *mask, int plane) {
  switch (length) {
    case 8:
      memcpy(mask, plane ? mask_8_uv : mask_8, length);
      break;
    case 16:
      memcpy(mask, plane ? mask_16_uv : mask_16, length);
      break;
    case 32:
      memcpy(mask, plane ? mask_32_uv : mask_32, length);
      break;
    default:
      assert(0);
  }
}


void vp10_build_masked_inter_predictor_complex(
    MACROBLOCKD *xd,
    uint8_t *dst, int dst_stride, uint8_t *dst2, int dst2_stride,
    const struct macroblockd_plane *pd, int mi_row, int mi_col,
    int mi_row_ori, int mi_col_ori, BLOCK_SIZE bsize, BLOCK_SIZE top_bsize,
    PARTITION_TYPE partition, int plane) {
  int i, j;
  uint8_t mask[MAXTXLEN];
  int top_w = 4 << b_width_log2_lookup[top_bsize],
      top_h = 4 << b_height_log2_lookup[top_bsize];
  int w = 4 << b_width_log2_lookup[bsize], h = 4 << b_height_log2_lookup[bsize];
  int w_offset = (mi_col - mi_col_ori) << 3,
      h_offset = (mi_row - mi_row_ori) << 3;

#if CONFIG_VP9_HIGHBITDEPTH
  uint16_t *dst16= CONVERT_TO_SHORTPTR(dst);
  uint16_t *dst216 = CONVERT_TO_SHORTPTR(dst2);
  int b_hdb = (xd->cur_buf->flags & YV12_FLAG_HIGHBITDEPTH) ? 1 : 0;
#endif  // CONFIG_VP9_HIGHBITDEPTH

  top_w >>= pd->subsampling_x;
  top_h >>= pd->subsampling_y;
  w >>= pd->subsampling_x;
  h >>= pd->subsampling_y;
  w_offset >>= pd->subsampling_x;
  h_offset >>= pd->subsampling_y;

  switch (partition) {
    case PARTITION_HORZ:
    {
#if CONFIG_VP9_HIGHBITDEPTH
      if (b_hdb) {
        uint16_t *dst_tmp = dst16 + h_offset * dst_stride;
        uint16_t *dst2_tmp = dst216 + h_offset * dst2_stride;
        generate_1dmask(h, mask + h_offset,
                        plane && xd->plane[plane].subsampling_y);

        for (i = h_offset; i < h_offset + h; i++) {
          for (j = 0; j < top_w; j++) {
            const int m = mask[i];  assert(m >= 0 && m <= 64);
            if (m == 64)
              continue;

            if (m == 0)
              dst_tmp[j] = dst2_tmp[j];
            else
              dst_tmp[j] = (dst_tmp[j] * m + dst2_tmp[j] * (64 - m) + 32) >> 6;
          }
          dst_tmp += dst_stride;
          dst2_tmp += dst2_stride;
        }

        for (; i < top_h; i ++) {
          memcpy(dst_tmp, dst2_tmp, top_w * sizeof(uint16_t));
          dst_tmp += dst_stride;
          dst2_tmp += dst2_stride;
        }
      } else {
#endif  // CONFIG_VP9_HIGHBITDEPTH
        uint8_t *dst_tmp = dst + h_offset * dst_stride;
        uint8_t *dst2_tmp = dst2 + h_offset * dst2_stride;
        generate_1dmask(h, mask + h_offset,
                        plane && xd->plane[plane].subsampling_y);

        for (i = h_offset; i < h_offset + h; i++) {
          for (j = 0; j < top_w; j++) {
            const int m = mask[i];  assert(m >= 0 && m <= 64);
            if (m == 64)
              continue;

            if (m == 0)
              dst_tmp[j] = dst2_tmp[j];
            else
              dst_tmp[j] = (dst_tmp[j] * m + dst2_tmp[j] * (64 - m) + 32) >> 6;
          }
          dst_tmp += dst_stride;
          dst2_tmp += dst2_stride;
        }

        for (; i < top_h; i ++) {
          memcpy(dst_tmp, dst2_tmp, top_w * sizeof(uint8_t));
          dst_tmp += dst_stride;
          dst2_tmp += dst2_stride;
        }
#if CONFIG_VP9_HIGHBITDEPTH
      }
#endif  // CONFIG_VP9_HIGHBITDEPTH
    }

      break;
    case PARTITION_VERT:
    {
#if CONFIG_VP9_HIGHBITDEPTH
      if (b_hdb) {
        uint16_t *dst_tmp = dst16;
        uint16_t *dst2_tmp = dst216;
        generate_1dmask(w, mask + w_offset,
                        plane && xd->plane[plane].subsampling_x);

        for (i = 0; i < top_h; i++) {
          for (j = w_offset; j < w_offset + w; j++) {
            const int m = mask[j];   assert(m >= 0 && m <= 64);
            if (m == 64)
              continue;

            if (m == 0)
              dst_tmp[j] = dst2_tmp[j];
            else
              dst_tmp[j] = (dst_tmp[j] * m + dst2_tmp[j] * (64 - m) + 32) >> 6;
          }
          memcpy(dst_tmp + j, dst2_tmp + j,
                     (top_w - w_offset - w) * sizeof(uint16_t));
          dst_tmp += dst_stride;
          dst2_tmp += dst2_stride;
        }
      } else {
#endif  // CONFIG_VP9_HIGHBITDEPTH
        uint8_t *dst_tmp = dst;
        uint8_t *dst2_tmp = dst2;
        generate_1dmask(w, mask + w_offset,
                        plane && xd->plane[plane].subsampling_x);

        for (i = 0; i < top_h; i++) {
          for (j = w_offset; j < w_offset + w; j++) {
            const int m = mask[j];   assert(m >= 0 && m <= 64);
            if (m == 64)
              continue;

            if (m == 0)
              dst_tmp[j] = dst2_tmp[j];
            else
              dst_tmp[j] = (dst_tmp[j] * m + dst2_tmp[j] * (64 - m) + 32) >> 6;
          }
            memcpy(dst_tmp + j, dst2_tmp + j,
                       (top_w - w_offset - w) * sizeof(uint8_t));
          dst_tmp += dst_stride;
          dst2_tmp += dst2_stride;
        }
#if CONFIG_VP9_HIGHBITDEPTH
      }
#endif  // CONFIG_VP9_HIGHBITDEPTH
    }
      break;
    default:
      assert(0);
  }
  (void) xd;
}

void vp10_build_inter_predictors_sb_sub8x8(MACROBLOCKD *xd,
                                           int mi_row, int mi_col,
                                           BLOCK_SIZE bsize, int block) {
  // Prediction function used in supertx:
  // Use the mv at current block (which is less than 8x8)
  // to get prediction of a block located at (mi_row, mi_col) at size of bsize
  // bsize can be larger than 8x8.
  // block (0-3): the sub8x8 location of current block
  int plane;
  const int mi_x = mi_col * MI_SIZE;
  const int mi_y = mi_row * MI_SIZE;

  // For sub8x8 uv:
  // Skip uv prediction in supertx except the first block (block = 0)
  int max_plane = block ? 1 : MAX_MB_PLANE;

  for (plane = 0; plane < max_plane; plane++) {
    const BLOCK_SIZE plane_bsize = get_plane_block_size(bsize,
                                                        &xd->plane[plane]);
    const int num_4x4_w = num_4x4_blocks_wide_lookup[plane_bsize];
    const int num_4x4_h = num_4x4_blocks_high_lookup[plane_bsize];
    const int bw = 4 * num_4x4_w;
    const int bh = 4 * num_4x4_h;

    build_inter_predictors(xd, plane, block, bw, bh,
                           0, 0, bw, bh,
                           mi_x, mi_y);
  }
}
#endif  // CONFIG_SUPERTX
