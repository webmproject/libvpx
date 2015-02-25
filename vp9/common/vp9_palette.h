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

#if CONFIG_PALETTE
int count_colors(const uint8_t *src, int stride, int rows, int cols);
void insertion_sort(double *data, int n);
int run_lengh_encoding(uint8_t *seq, int n, uint16_t *runs, int max_run);
int run_lengh_decoding(uint16_t *runs, int l, uint8_t *seq);
void transpose_block(uint8_t *seq_in, uint8_t *seq_out, int rows, int cols);
void palette_color_insertion(uint8_t *old_colors, int *m, int *count,
                             MB_MODE_INFO *mbmi);
int palette_color_lookup(uint8_t *dic, int n, uint8_t val, int bits);
int palette_max_run(BLOCK_SIZE bsize);
int get_bit_depth(int n);
int k_means(double *data, double *centroids, int *indices,
             int n, int k, int dim, int max_itr);
void calc_indices(double *data, double *centroids, int *indices,
                  int n, int k, int dim);
void zz_scan_order(int *order, int rows, int cols);
void spiral_scan_order(int *order, int rows, int cols);
void palette_scan(uint8_t *color_index_map, uint8_t *sequence,
                  int rows, int cols, PALETTE_SCAN_ORDER ps, int *scan_order);
void palette_iscan(uint8_t *color_index_map, uint8_t *sequence,
                   int rows, int cols, PALETTE_SCAN_ORDER ps, int *scan_order);
#endif

#endif  // VP9_COMMON_VP9_PALETTE_H_
