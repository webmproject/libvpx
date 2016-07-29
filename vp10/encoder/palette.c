/*
 *  Copyright (c) 2015 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include <math.h>
#include "vp10/encoder/palette.h"

static float calc_dist(const float *p1, const float *p2, int dim) {
  float dist = 0;
  int i = 0;

  for (i = 0; i < dim; ++i) {
    float diff = p1[i] - roundf(p2[i]);
    dist += diff * diff;
  }
  return dist;
}

void vp10_calc_indices(const float *data, const float *centroids,
                       uint8_t *indices, int n, int k, int dim) {
  int i, j;
  float min_dist, this_dist;

  for (i = 0; i < n; ++i) {
    min_dist = calc_dist(data + i * dim, centroids, dim);
    indices[i] = 0;
    for (j = 1; j < k; ++j) {
      this_dist = calc_dist(data + i * dim, centroids + j * dim, dim);
      if (this_dist < min_dist) {
        min_dist = this_dist;
        indices[i] = j;
      }
    }
  }
}

// Generate a random number in the range [0, 32768).
static unsigned int lcg_rand16(unsigned int *state) {
  *state = *state * 1103515245 + 12345;
  return *state / 65536 % 32768;
}

static void calc_centroids(const float *data, float *centroids,
                           const uint8_t *indices, int n, int k, int dim) {
  int i, j, index;
  int count[PALETTE_MAX_SIZE];
  unsigned int rand_state = (unsigned int)data[0];

  assert(n <= 32768);

  memset(count, 0, sizeof(count[0]) * k);
  memset(centroids, 0, sizeof(centroids[0]) * k * dim);

  for (i = 0; i < n; ++i) {
    index = indices[i];
    assert(index < k);
    ++count[index];
    for (j = 0; j < dim; ++j) {
      centroids[index * dim + j] += data[i * dim + j];
    }
  }

  for (i = 0; i < k; ++i) {
    if (count[i] == 0) {
      memcpy(centroids + i * dim, data + (lcg_rand16(&rand_state) % n) * dim,
                 sizeof(centroids[0]) * dim);
    } else {
      const float norm = 1.0f / count[i];
      for (j = 0; j < dim; ++j)
        centroids[i * dim + j] *= norm;
    }
  }
}

static float calc_total_dist(const float *data, const float *centroids,
                              const uint8_t *indices, int n, int k, int dim) {
  float dist = 0;
  int i;
  (void) k;

  for (i = 0; i < n; ++i)
    dist += calc_dist(data + i * dim, centroids + indices[i] * dim, dim);

  return dist;
}

int vp10_k_means(const float *data, float *centroids, uint8_t *indices,
                 uint8_t *pre_indices, int n, int k, int dim, int max_itr) {
  int i = 0;
  float pre_dist, this_dist;
  float pre_centroids[2 * PALETTE_MAX_SIZE];

  vp10_calc_indices(data, centroids, indices, n, k, dim);
  pre_dist = calc_total_dist(data, centroids, indices, n, k, dim);
  memcpy(pre_centroids, centroids, sizeof(pre_centroids[0]) * k * dim);
  memcpy(pre_indices, indices, sizeof(pre_indices[0]) * n);
  while (i < max_itr) {
    calc_centroids(data, centroids, indices, n, k, dim);
    vp10_calc_indices(data, centroids, indices, n, k, dim);
    this_dist = calc_total_dist(data, centroids, indices, n, k, dim);

    if (this_dist > pre_dist) {
      memcpy(centroids, pre_centroids, sizeof(pre_centroids[0]) * k * dim);
      memcpy(indices, pre_indices, sizeof(pre_indices[0]) * n);
      break;
    }
    if (!memcmp(centroids, pre_centroids, sizeof(pre_centroids[0]) * k * dim))
      break;

    memcpy(pre_centroids, centroids, sizeof(pre_centroids[0]) * k * dim);
    memcpy(pre_indices, indices, sizeof(pre_indices[0]) * n);
    pre_dist = this_dist;
    ++i;
  }

  return i;
}

void vp10_insertion_sort(float *data, int n) {
  int i, j, k;
  float val;

  if (n <= 1)
    return;

  for (i = 1; i < n; ++i) {
    val = data[i];
    j = 0;
    while (val > data[j] && j < i)
      ++j;

    if (j == i)
      continue;

    for (k = i; k > j; --k)
      data[k] = data[k - 1];
    data[j] = val;
  }
}

int vp10_count_colors(const uint8_t *src, int stride, int rows, int cols) {
  int n = 0, r, c, i, val_count[256];
  uint8_t val;
  memset(val_count, 0, sizeof(val_count));

  for (r = 0; r < rows; ++r) {
    for (c = 0; c < cols; ++c) {
      val = src[r * stride + c];
      ++val_count[val];
    }
  }

  for (i = 0; i < 256; ++i) {
    if (val_count[i]) {
      ++n;
    }
  }

  return n;
}

#if CONFIG_VPX_HIGHBITDEPTH
int vp10_count_colors_highbd(const uint8_t *src8, int stride, int rows,
                             int cols, int bit_depth) {
  int n = 0, r, c, i;
  uint16_t val;
  uint16_t *src = CONVERT_TO_SHORTPTR(src8);
  int val_count[1 << 12];

  assert(bit_depth <= 12);
  memset(val_count, 0, (1 << 12) * sizeof(val_count[0]));
  for (r = 0; r < rows; ++r) {
    for (c = 0; c < cols; ++c) {
      val = src[r * stride + c];
      ++val_count[val];
    }
  }

  for (i = 0; i < (1 << bit_depth); ++i) {
    if (val_count[i]) {
      ++n;
    }
  }

  return n;
}
#endif  // CONFIG_VPX_HIGHBITDEPTH


