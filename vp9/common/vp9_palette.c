/*
 *  Copyright (c) 2015 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */


#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "vp9/common/vp9_palette.h"

#if CONFIG_PALETTE
void vp9_insertion_sort(double *data, int n) {
  int i, j, k;
  double val;

  if (n <= 1)
    return;

  for (i = 1; i < n; i++) {
    val = data[i];
    j = 0;
    while (val > data[j] && j < i)
      j++;

    if (j == i)
      continue;

    for (k = i; k > j; k--)
      data[k] = data[k - 1];
    data[j] = val;
  }
}

int vp9_count_colors(const uint8_t *src, int stride, int rows, int cols) {
  int n = 0, r, c, i, val_count[256];
  uint8_t val;
  vpx_memset(val_count, 0, sizeof(val_count));

  for (r = 0; r < rows; r++) {
      for (c = 0; c < cols; c++) {
        val = src[r * stride + c];
        val_count[val]++;
      }
    }

    for (i = 0; i < 256; i++) {
      if (val_count[i]) {
        n++;
      }
    }

    return n;
}

void vp9_palette_color_insertion(uint8_t *old_colors, int *m, int *count,
                                 MB_MODE_INFO *mbmi) {
  int k = *m, n = mbmi->palette_literal_size;
  int i, j, l, min_idx = -1;
  uint8_t *new_colors = mbmi->palette_literal_colors;
  uint8_t val;

  if (mbmi->palette_indexed_size > 0) {
    for (i = 0; i < mbmi->palette_indexed_size; i++)
      count[mbmi->palette_indexed_colors[i]] +=
          (8 - abs(mbmi->palette_color_delta[i]));
  }

  i = 0;
  while (i < k) {
    count[i] -= 1;
    i++;
  }

  if (n <= 0)
    return;

  for (i = 0; i < n; i++) {
    val = new_colors[i];
    j = 0;
    while (val != old_colors[j] && j < k)
      j++;
    if (j < k && val == old_colors[j]) {
      count[j] += 8;
      continue;
    }

    if (k + 1 > PALETTE_BUF_SIZE) {
      min_idx = 0;
      for (l = 1; l < k; l++)
        if (count[l] < count[min_idx])
          min_idx = l;
      old_colors[min_idx] = val;
      count[min_idx] = 8;
    } else {
      old_colors[k] = val;
      count[k] = 8;
      k++;
    }
  }

  *m = k;
}

int vp9_palette_color_lookup(uint8_t *dic, int n, uint8_t val, int bits) {
  int j, min, arg_min = 0, i = 1;

  if (n < 1)
    return -1;

  min = abs(val - dic[0]);
  arg_min = 0;
  while (i < n) {
    j = abs(val - dic[i]);
    if (j < min) {
      min = j;
      arg_min = i;
    }
    i++;
  }

  if (min < (1 << bits))
    return arg_min;
  else
    return -1;
}

int vp9_get_bit_depth(int n) {
  int i = 1, p = 2;
  while (p < n) {
    i++;
    p = p << 1;
  }

  return i;
}

static double calc_dist(double *p1, double *p2, int dim) {
  double dist = 0;
  int i = 0;

  for (i = 0; i < dim; i++) {
    dist = dist + (p1[i] - p2[i]) * (p1[i] - p2[i]);
  }
  return dist;
}

void vp9_calc_indices(double *data, double *centroids, int *indices,
                      int n, int k, int dim) {
  int i, j;
  double min_dist, this_dist;

  for (i = 0; i < n; i++) {
    min_dist = calc_dist(data + i * dim, centroids, dim);
    indices[i] = 0;
    for (j = 0; j < k; j++) {
      this_dist = calc_dist(data + i * dim, centroids + j * dim, dim);
      if (this_dist < min_dist) {
        min_dist = this_dist;
        indices[i] = j;
      }
    }
  }
}

static void calc_centroids(double *data, double *centroids, int *indices,
                           int n, int k, int dim) {
  int i, j, index;
  int count[256];
  unsigned int seed = time(NULL);
  vpx_memset(count, 0, sizeof(count[0]) * k);
  vpx_memset(centroids, 0, sizeof(centroids[0]) * k * dim);

  for (i = 0; i < n; i++) {
    index = indices[i];
    count[index]++;
    for (j = 0; j < dim; j++) {
      centroids[index * dim + j] += data[i * dim + j];
    }
  }

  for (i = 0; i < k; i++) {
    if (!count[i])
      vpx_memcpy(centroids + i * dim, data + (rand_r(&seed) % n) * dim,
                 sizeof(centroids[0]) * dim);
    else
      for (j = 0; j < dim; j++)
        centroids[i * dim + j] /= count[i];
  }
}

double calc_total_dist(double *data, double *centroids, int *indices,
                       int n, int k, int dim) {
  double dist = 0;
  int i;
  (void) k;

  for (i = 0; i < n; i++) {
    dist += calc_dist(data + i * dim, centroids + indices[i] * dim, dim);
  }

  return dist;
}

int vp9_k_means(double *data, double *centroids, int *indices,
                int n, int k, int dim, int max_itr) {
  int i = 0;
  int *pre_indices;
  double pre_total_dist, cur_total_dist;
  double pre_centroids[256];

  pre_indices = vpx_memalign(16, n * sizeof(indices[0]));
  vp9_calc_indices(data, centroids, indices, n, k, dim);
  pre_total_dist = calc_total_dist(data, centroids, indices, n, k, dim);
  vpx_memcpy(pre_centroids, centroids, sizeof(pre_centroids[0]) * k * dim);
  vpx_memcpy(pre_indices, indices, sizeof(pre_indices[0]) * n);
  while (i < max_itr) {
    calc_centroids(data, centroids, indices, n, k, dim);
    vp9_calc_indices(data, centroids, indices, n, k, dim);
    cur_total_dist = calc_total_dist(data, centroids, indices, n, k, dim);

    if (cur_total_dist > pre_total_dist && 0) {
      vpx_memcpy(centroids, pre_centroids, sizeof(pre_centroids[0]) * k * dim);
      vpx_memcpy(indices, pre_indices, sizeof(pre_indices[0]) * n);
      break;
    }
    if (!memcmp(centroids, pre_centroids, sizeof(pre_centroids[0]) * k * dim))
      break;

    vpx_memcpy(pre_centroids, centroids, sizeof(pre_centroids[0]) * k * dim);
    vpx_memcpy(pre_indices, indices, sizeof(pre_indices[0]) * n);
    pre_total_dist = cur_total_dist;
    i++;
  }

  vpx_free(pre_indices);
  return i;
}

void vp9_update_palette_counts(FRAME_COUNTS *counts, MB_MODE_INFO *mbmi,
                               BLOCK_SIZE bsize, int palette_ctx) {
  int idx = bsize - BLOCK_8X8;

  counts->y_palette_enabled[idx][palette_ctx][mbmi->palette_enabled[0]]++;
  counts->uv_palette_enabled[mbmi->palette_enabled[0]]
                            [mbmi->palette_enabled[1]]++;
  if (mbmi->palette_enabled[0])
    counts->y_palette_size[idx][mbmi->palette_size[0] - 2]++;
  if (mbmi->palette_enabled[1])
    counts->uv_palette_size[idx][mbmi->palette_size[1] - 2]++;
}

static const int palette_color_context_lookup[PALETTE_COLOR_CONTEXTS] = {
    3993,  4235,  4378,  4380,  // (3, 0, 0, 0), (3, 2, 0, 0),
                                // (3, 3, 2, 0), (3, 3, 2, 2),
    5720,  6655,  7018,  7040,  // (4, 3, 3, 0), (5, 0, 0, 0),
                                // (5, 3, 0, 0), (5, 3, 2, 0),
    7260,  8228,  8250,  8470,  // (5, 5, 0, 0), (6, 2, 0, 0),
                                // (6, 2, 2, 0), (6, 4, 0, 0),
    9680, 10648, 10890, 13310   // (7, 3, 0, 0), (8, 0, 0, 0),
                                // (8, 2, 0, 0), (10, 0, 0, 0)
};

int vp9_get_palette_color_context(uint8_t *color_map, int cols,
                                  int r, int c, int n, int *color_order) {
  int i, j, max, max_idx, temp;
  int scores[PALETTE_MAX_SIZE];
  int weights[4] = {3, 2, 3, 2};
  int color_ctx = 0;
  int color_neighbors[4];

  assert(n <= PALETTE_MAX_SIZE);

  if (c - 1 >= 0)
    color_neighbors[0] = color_map[r * cols + c - 1];
  else
    color_neighbors[0] = -1;
  if (c - 1 >= 0 && r - 1 >= 0)
    color_neighbors[1] = color_map[(r - 1) * cols + c - 1];
  else
    color_neighbors[1] = -1;
  if (r - 1 >= 0)
    color_neighbors[2] = color_map[(r - 1) * cols + c];
  else
    color_neighbors[2] = -1;
  if (r - 1 >= 0 && c + 1 <= cols - 1)
    color_neighbors[3] = color_map[(r - 1) * cols + c + 1];
  else
    color_neighbors[3] = -1;

  for (i = 0; i < PALETTE_MAX_SIZE; i++)
    color_order[i] = i;
  memset(scores, 0, PALETTE_MAX_SIZE * sizeof(scores[0]));
  for (i = 0; i < 4; i++) {
    if (color_neighbors[i] >= 0)
      scores[color_neighbors[i]] += weights[i];
  }

  for (i = 0; i < 4; i++) {
    max = scores[i];
    max_idx = i;
    j = i + 1;
    while (j < n) {
      if (scores[j] > max) {
        max = scores[j];
        max_idx = j;
      }
      j++;
    }

    if (max_idx != i) {
      temp = scores[i];
      scores[i] = scores[max_idx];
      scores[max_idx] = temp;

      temp = color_order[i];
      color_order[i] = color_order[max_idx];
      color_order[max_idx] = temp;
    }
  }

  for (i = 0; i < 4; i++)
    color_ctx = color_ctx * 11 + scores[i];

  for (i = 0; i < PALETTE_COLOR_CONTEXTS; i++)
    if (color_ctx == palette_color_context_lookup[i]) {
      color_ctx = i;
      break;
    }

  return color_ctx;
}
#endif  // CONFIG_PALETTE
