/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#ifndef VP9_COMMON_VP9_RECONINTER_H_
#define VP9_COMMON_VP9_RECONINTER_H_

#include "vpx/vpx_integer.h"
#include "vp9/common/vp9_onyxc_int.h"

#ifdef __cplusplus
extern "C" {
#endif

void vp9_build_inter_predictors_sby(MACROBLOCKD *xd, int mi_row, int mi_col,
                                    BLOCK_SIZE bsize);

void vp9_build_inter_predictors_sbuv(MACROBLOCKD *xd, int mi_row, int mi_col,
                                     BLOCK_SIZE bsize);

void vp9_build_inter_predictors_sb(MACROBLOCKD *xd, int mi_row, int mi_col,
                                   BLOCK_SIZE bsize);

void vp9_dec_build_inter_predictors_sb(MACROBLOCKD *xd, int mi_row, int mi_col,
                                       BLOCK_SIZE bsize);

#if CONFIG_SUPERTX
void vp9_build_inter_predictors_sb_sub8x8(MACROBLOCKD *xd,
                                          int mi_row, int mi_col,
                                          BLOCK_SIZE bsize, int block);
void vp9_dec_build_inter_predictors_sb_sub8x8(MACROBLOCKD *xd,
                                       int mi_row, int mi_col,
                                       BLOCK_SIZE bsize, int block);
#endif

void vp9_build_inter_predictor(const uint8_t *src, int src_stride,
                               uint8_t *dst, int dst_stride,
                               const MV *mv_q3,
                               const struct scale_factors *sf,
                               int w, int h, int do_avg,
                               const InterpKernel *kernel,
                               enum mv_precision precision,
                               int x, int y);

#if CONFIG_VP9_HIGHBITDEPTH
void vp9_highbd_build_inter_predictor(const uint8_t *src, int src_stride,
                                      uint8_t *dst, int dst_stride,
                                      const MV *mv_q3,
                                      const struct scale_factors *sf,
                                      int w, int h, int do_avg,
                                      const InterpKernel *kernel,
                                      enum mv_precision precision,
                                      int x, int y, int bd);
#endif

static INLINE int scaled_buffer_offset(int x_offset, int y_offset, int stride,
                                       const struct scale_factors *sf) {
  const int x = sf ? sf->scale_value_x(x_offset, sf) : x_offset;
  const int y = sf ? sf->scale_value_y(y_offset, sf) : y_offset;
  return y * stride + x;
}

static INLINE void setup_pred_plane(struct buf_2d *dst,
                                    int width, int height,
                                    uint8_t *src, int stride,
                                    int mi_row, int mi_col,
                                    const struct scale_factors *scale,
                                    int subsampling_x, int subsampling_y) {
  const int x = (MI_SIZE * mi_col) >> subsampling_x;
  const int y = (MI_SIZE * mi_row) >> subsampling_y;
  dst->buf0 = src;
  dst->buf = src + scaled_buffer_offset(x, y, stride, scale);
  dst->stride = stride;
  dst->width = width;
  dst->height = height;
}

void vp9_setup_dst_planes(struct macroblockd_plane planes[MAX_MB_PLANE],
                          const YV12_BUFFER_CONFIG *src,
                          int mi_row, int mi_col);

void vp9_setup_pre_planes(MACROBLOCKD *xd, int idx,
                          const YV12_BUFFER_CONFIG *src, int mi_row, int mi_col,
                          const struct scale_factors *sf);

#if CONFIG_WEDGE_PARTITION

#define MASK_MASTER_SIZE   (2 * CODING_UNIT_SIZE)
#define MASK_MASTER_STRIDE (2 * CODING_UNIT_SIZE)

void vp9_init_wedge_masks();

const uint8_t *vp9_get_soft_mask(int wedge_index,
                                 BLOCK_SIZE sb_type,
                                 int h, int w);
void vp9_generate_soft_mask(int wedge_index, BLOCK_SIZE sb_type,
                            int h, int w, uint8_t *mask, int stride);
void vp9_generate_hard_mask(int wedge_index, BLOCK_SIZE sb_type,
                            int h, int w, uint8_t *mask, int stride);
void vp9_build_inter_predictors_for_planes_single_buf(
    MACROBLOCKD *xd, BLOCK_SIZE bsize,
    int mi_row, int mi_col, int ref,
    uint8_t *ext_dst[3], int ext_dst_stride[3]);
void vp9_build_wedge_inter_predictor_from_buf(
    MACROBLOCKD *xd, BLOCK_SIZE bsize,
    int mi_row, int mi_col,
    uint8_t *ext_dst0[3], int ext_dst_stride0[3],
    uint8_t *ext_dst1[3], int ext_dst_stride1[3]);
#endif  // CONFIG_WEDGE_PARTITION

#if CONFIG_SUPERTX

struct macroblockd_plane;

void vp9_build_masked_inter_predictor_complex(
    MACROBLOCKD *xd,
    uint8_t *dst, int dst_stride, uint8_t *dst2, int dst2_stride,
    const struct macroblockd_plane *pd, int mi_row, int mi_col,
    int mi_row_ori, int mi_col_ori, BLOCK_SIZE bsize, BLOCK_SIZE top_bsize,
    PARTITION_TYPE partition, int plane);

#if CONFIG_WEDGE_PARTITION
void vp9_build_inter_predictors_sb_extend(MACROBLOCKD *xd,
                                          int mi_row, int mi_col,
                                          int mi_row_ori, int mi_col_ori,
                                          BLOCK_SIZE bsize);
void vp9_dec_build_inter_predictors_sb_extend(MACROBLOCKD *xd,
                                              int mi_row, int mi_col,
                                              int mi_row_ori, int mi_col_ori,
                                              BLOCK_SIZE bsize);

void vp9_build_inter_predictors_sb_sub8x8_extend(
    MACROBLOCKD *xd,
    int mi_row, int mi_col,
    int mi_row_ori, int mi_col_ori,
    BLOCK_SIZE bsize, int block);
void vp9_dec_build_inter_predictors_sb_sub8x8_extend(
    MACROBLOCKD *xd,
    int mi_row, int mi_col,
    int mi_row_ori, int mi_col_ori,
    BLOCK_SIZE bsize, int block);

#endif  // CONFIG_WEDGE_PARTITION
#endif  // CONFIG_SUPERTX

#ifdef __cplusplus
}  // extern "C"
#endif

#endif  // VP9_COMMON_VP9_RECONINTER_H_
