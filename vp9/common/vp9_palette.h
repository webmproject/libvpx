/*
 *  Copyright (c) 2015 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#ifndef VP9_COMMON_VP9_PALETTE_H_
#define VP9_COMMON_VP9_PALETTE_H_

#include "vp9/common/vp9_blockd.h"
#include "vp9/common/vp9_entropymode.h"

#if CONFIG_PALETTE
int vp9_count_colors(const uint8_t *src, int stride, int rows, int cols);
void vp9_insertion_sort(double *data, int n);
void vp9_palette_color_insertion(uint8_t *old_colors, int *m, int *count,
                                 MB_MODE_INFO *mbmi);
int vp9_palette_color_lookup(uint8_t *dic, int n, uint8_t val, int bits);
int vp9_get_bit_depth(int n);
int vp9_k_means(double *data, double *centroids, int *indices,
                int n, int k, int dim, int max_itr);
void vp9_calc_indices(double *data, double *centroids, int *indices,
                      int n, int k, int dim);
void vp9_update_palette_counts(FRAME_COUNTS *counts, MB_MODE_INFO *mbmi,
                               BLOCK_SIZE bsize, int palette_ctx);
int vp9_get_palette_color_context(uint8_t *color_map, int cols,
                                  int r, int c, int n, int *color_order);
#endif

#endif  // VP9_COMMON_VP9_PALETTE_H_
