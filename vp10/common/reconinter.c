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
#if CONFIG_OBMC
#include "vp10/common/onyxc_int.h"
#endif  // CONFIG_OBMC

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

void build_inter_predictors(MACROBLOCKD *xd, int plane,
#if CONFIG_OBMC
                            int mi_col_offset, int mi_row_offset,
#endif  // CONFIG_OBMC
                            int block,
                            int bw, int bh,
                            int x, int y, int w, int h,
                            int mi_x, int mi_y) {
  struct macroblockd_plane *const pd = &xd->plane[plane];
#if CONFIG_OBMC
  const MODE_INFO *mi = xd->mi[mi_col_offset + xd->mi_stride * mi_row_offset];
#else
  const MODE_INFO *mi = xd->mi[0];
#endif  // CONFIG_OBMC
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
           build_inter_predictors(xd, plane,
#if CONFIG_OBMC
                                  0, 0,
#endif  // CONFIG_OBMC
                                  y * 2 + x, bw, bh,
                                  4 * x, 4 * y, pw, ph, mi_x, mi_y);
    } else {
      build_inter_predictors(xd, plane,
#if CONFIG_OBMC
                             0, 0,
#endif  // CONFIG_OBMC
                             0, bw, bh,
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

    build_inter_predictors(xd, plane,
#if CONFIG_OBMC
                           0, 0,
#endif  // CONFIG_OBMC
                           block, bw, bh,
                           0, 0, bw, bh,
                           mi_x, mi_y);
  }
}
#endif  // CONFIG_SUPERTX

#if CONFIG_OBMC
// obmc_mask_N[is_neighbor_predictor][overlap_position]
static const uint8_t obmc_mask_1[2][1] = {
    { 55},
    {  9}
};

static const uint8_t obmc_mask_2[2][2] = {
    { 45, 62},
    { 19,  2}
};

static const uint8_t obmc_mask_4[2][4] = {
    { 39, 50, 59, 64},
    { 25, 14,  5,  0}
};

static const uint8_t obmc_mask_8[2][8] = {
    { 36, 42, 48, 53, 57, 61, 63, 64},
    { 28, 22, 16, 11,  7,  3,  1,  0}
};

static const uint8_t obmc_mask_16[2][16] = {
    { 34, 37, 40, 43, 46, 49, 52, 54, 56, 58, 60, 61, 63, 64, 64, 64},
    { 30, 27, 24, 21, 18, 15, 12, 10,  8,  6,  4,  3,  1,  0,  0,  0}
};

static const uint8_t obmc_mask_32[2][32] = {
    { 33, 35, 36, 38, 40, 41, 43, 44, 45, 47, 48, 50, 51, 52, 53, 55,
      56, 57, 58, 59, 60, 60, 61, 62, 62, 63, 63, 64, 64, 64, 64, 64},
    { 31, 29, 28, 26, 24, 23, 21, 20, 19, 17, 16, 14, 13, 12, 11,  9,
       8,  7,  6,  5,  4,  4,  3,  2,  2,  1,  1,  0,  0,  0,  0,  0}
};

void setup_obmc_mask(int length, const uint8_t *mask[2]) {
  switch (length) {
    case 1:
      mask[0] = obmc_mask_1[0];
      mask[1] = obmc_mask_1[1];
      break;
    case 2:
      mask[0] = obmc_mask_2[0];
      mask[1] = obmc_mask_2[1];
      break;
    case 4:
      mask[0] = obmc_mask_4[0];
      mask[1] = obmc_mask_4[1];
      break;
    case 8:
      mask[0] = obmc_mask_8[0];
      mask[1] = obmc_mask_8[1];
      break;
    case 16:
      mask[0] = obmc_mask_16[0];
      mask[1] = obmc_mask_16[1];
      break;
    case 32:
      mask[0] = obmc_mask_32[0];
      mask[1] = obmc_mask_32[1];
      break;
    default:
      mask[0] = obmc_mask_32[0];
      mask[1] = obmc_mask_32[1];
      assert(0);
      break;
  }
}

// This function combines motion compensated predictions that is generated by
// top/left neighboring blocks' inter predictors with the regular inter
// prediction. We assume the original prediction (bmc) is stored in
// xd->plane[].dst.buf
void vp10_build_obmc_inter_prediction(VP10_COMMON *cm,
                                      MACROBLOCKD *xd, int mi_row, int mi_col,
                                      int use_tmp_dst_buf,
                                      uint8_t *final_buf[MAX_MB_PLANE],
                                      int final_stride[MAX_MB_PLANE],
                                      uint8_t *tmp_buf1[MAX_MB_PLANE],
                                      int tmp_stride1[MAX_MB_PLANE],
                                      uint8_t *tmp_buf2[MAX_MB_PLANE],
                                      int tmp_stride2[MAX_MB_PLANE]) {
  const TileInfo *const tile = &xd->tile;
  BLOCK_SIZE bsize = xd->mi[0]->mbmi.sb_type;
  int plane, i, mi_step;
#if CONFIG_VP9_HIGHBITDEPTH
  int is_hbd = (xd->cur_buf->flags & YV12_FLAG_HIGHBITDEPTH) ? 1 : 0;
#endif  // CONFIG_VP9_HIGHBITDEPTH

  if (use_tmp_dst_buf) {
    for (plane = 0; plane < MAX_MB_PLANE; ++plane) {
      const struct macroblockd_plane *pd = &xd->plane[plane];
      int bw = (xd->n8_w * 8) >> pd->subsampling_x;
      int bh = (xd->n8_h * 8) >> pd->subsampling_y;
      int row;
#if CONFIG_VP9_HIGHBITDEPTH
      if (is_hbd) {
        uint16_t *final_buf16 = CONVERT_TO_SHORTPTR(final_buf[plane]);
        uint16_t *bmc_buf16 = CONVERT_TO_SHORTPTR(pd->dst.buf);
        for (row = 0; row < bh; ++row)
          memcpy(final_buf16 + row * final_stride[plane],
                 bmc_buf16 + row * pd->dst.stride, bw * sizeof(uint16_t));
      } else {
#endif
      for (row = 0; row < bh; ++row)
        memcpy(final_buf[plane] + row * final_stride[plane],
               pd->dst.buf + row * pd->dst.stride, bw);
#if CONFIG_VP9_HIGHBITDEPTH
      }
#endif  // CONFIG_VP9_HIGHBITDEPTH
    }
  }

  // handle above row
  for (i = 0; mi_row > 0 && i < VPXMIN(xd->n8_w, cm->mi_cols - mi_col);
       i += mi_step) {
    int mi_row_offset = -1;
    int mi_col_offset = i;
    int overlap;
    MODE_INFO *above_mi = xd->mi[mi_col_offset +
                                 mi_row_offset * xd->mi_stride];
    MB_MODE_INFO *above_mbmi = &above_mi->mbmi;

    mi_step = VPXMIN(xd->n8_w,
                     num_8x8_blocks_wide_lookup[above_mbmi->sb_type]);

    if (!is_inter_block(above_mbmi))
      continue;

    overlap = (above_mbmi->skip) ?
              num_4x4_blocks_high_lookup[bsize] << 1 :
              VPXMIN(num_4x4_blocks_high_lookup[bsize],
                     num_4x4_blocks_high_lookup[above_mbmi->sb_type]) << 1;

    for (plane = 0; plane < MAX_MB_PLANE; ++plane) {
      const struct macroblockd_plane *pd = &xd->plane[plane];
      int bw = (mi_step * 8) >> pd->subsampling_x;
      int bh = overlap >> pd->subsampling_y;
      int row, col;
      int dst_stride = use_tmp_dst_buf ? final_stride[plane] : pd->dst.stride;
      uint8_t *dst = use_tmp_dst_buf ?
          &final_buf[plane][(i * 8) >> pd->subsampling_x] :
          &pd->dst.buf[(i * 8) >> pd->subsampling_x];
      int bmc_stride = pd->dst.stride;
      uint8_t *bmc = &pd->dst.buf[(i * 8) >> pd->subsampling_x];
      int tmp_stride = tmp_stride1[plane];
      uint8_t *tmp = &tmp_buf1[plane][(i * 8) >> pd->subsampling_x];
      const uint8_t *mask[2];

      setup_obmc_mask(bh, mask);

#if CONFIG_VP9_HIGHBITDEPTH
      if (is_hbd) {
        uint16_t *dst16 = CONVERT_TO_SHORTPTR(dst);
        uint16_t *bmc16 = CONVERT_TO_SHORTPTR(bmc);
        uint16_t *tmp16 = CONVERT_TO_SHORTPTR(tmp);

        for (row = 0; row < bh; ++row) {
          for (col = 0; col < bw; ++col) {
            dst16[col] = (mask[0][row] * bmc16[col] + mask[1][row] * tmp16[col]
                          + 32) >> 6;
          }
          dst16 += dst_stride;
          bmc16 += bmc_stride;
          tmp16 += tmp_stride;
        }
      } else {
#endif  // CONFIG_VP9_HIGHBITDEPTH
      for (row = 0; row < bh; ++row) {
        for (col = 0; col < bw; ++col) {
          dst[col] = (mask[0][row] * bmc[col] + mask[1][row] * tmp[col] + 32)
                     >> 6;
        }
        dst += dst_stride;
        bmc += bmc_stride;
        tmp += tmp_stride;
      }
#if CONFIG_VP9_HIGHBITDEPTH
      }
#endif  // CONFIG_VP9_HIGHBITDEPTH
    }
  }  // each mi in the above row

  if (mi_col == 0 || (mi_col - 1 < tile->mi_col_start) ||
      (mi_col - 1) >= tile->mi_col_end)
    return;
  // handle left column
  for (i = 0; i < VPXMIN(xd->n8_h, cm->mi_rows - mi_row);
       i += mi_step) {
    int mi_row_offset = i;
    int mi_col_offset = -1;
    int overlap;
    MODE_INFO *left_mi = xd->mi[mi_col_offset +
                                mi_row_offset * xd->mi_stride];
    MB_MODE_INFO *left_mbmi = &left_mi->mbmi;

    mi_step = VPXMIN(xd->n8_h,
                     num_8x8_blocks_high_lookup[left_mbmi->sb_type]);

    if (!is_inter_block(left_mbmi))
      continue;

    overlap = (left_mbmi->skip) ?
              num_4x4_blocks_wide_lookup[bsize] << 1 :
              VPXMIN(num_4x4_blocks_wide_lookup[bsize],
                     num_4x4_blocks_wide_lookup[left_mbmi->sb_type]) << 1;

    for (plane = 0; plane < MAX_MB_PLANE; ++plane) {
      const struct macroblockd_plane *pd = &xd->plane[plane];
      int bw = overlap >> pd->subsampling_x;
      int bh = (mi_step * 8) >> pd->subsampling_y;
      int row, col;
      int dst_stride = use_tmp_dst_buf ? final_stride[plane] : pd->dst.stride;
      uint8_t *dst = use_tmp_dst_buf ?
          &final_buf[plane][(i * 8 * dst_stride) >> pd->subsampling_y] :
          &pd->dst.buf[(i * 8 * dst_stride) >> pd->subsampling_y];
      int bmc_stride = pd->dst.stride;
      uint8_t *bmc = &pd->dst.buf[(i * 8 * bmc_stride) >> pd->subsampling_y];
      int tmp_stride = tmp_stride2[plane];
      uint8_t *tmp = &tmp_buf2[plane]
                              [(i * 8 * tmp_stride) >> pd->subsampling_y];
      const uint8_t *mask[2];

      setup_obmc_mask(bw, mask);

#if CONFIG_VP9_HIGHBITDEPTH
      if (is_hbd) {
        uint16_t *dst16 = CONVERT_TO_SHORTPTR(dst);
        uint16_t *bmc16 = CONVERT_TO_SHORTPTR(bmc);
        uint16_t *tmp16 = CONVERT_TO_SHORTPTR(tmp);

        for (row = 0; row < bh; ++row) {
          for (col = 0; col < bw; ++col) {
            dst16[col] = (mask[0][row] * bmc16[col] + mask[1][row] * tmp16[col]
                          + 32) >> 6;
          }
          dst16 += dst_stride;
          bmc16 += bmc_stride;
          tmp16 += tmp_stride;
        }
      } else {
#endif  // CONFIG_VP9_HIGHBITDEPTH
      for (row = 0; row < bh; ++row) {
        for (col = 0; col < bw; ++col) {
          dst[col] = (mask[0][col] * bmc[col] + mask[1][col] * tmp[col] + 32)
                     >> 6;
        }
        dst += dst_stride;
        bmc += bmc_stride;
        tmp += tmp_stride;
      }
#if CONFIG_VP9_HIGHBITDEPTH
      }
#endif  // CONFIG_VP9_HIGHBITDEPTH
    }
  }  // each mi in the left column
}
#endif  // CONFIG_OBMC
