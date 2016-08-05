/*
 * Copyright (c) 2016, Alliance for Open Media. All rights reserved
 *
 * This source code is subject to the terms of the BSD 2 Clause License and
 * the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
 * was not distributed with this source code in the LICENSE file, you can
 * obtain it at www.aomedia.org/license/software. If the Alliance for Open
 * Media Patent License 1.0 was not distributed with this source code in the
 * PATENTS file, you can obtain it at www.aomedia.org/license/patent.
 */

#include "av1/common/clpf.h"
#include "./aom_dsp_rtcd.h"
#include "aom/aom_integer.h"
#include "av1/common/quant_common.h"

// Calculate the error of a filtered and unfiltered block
void aom_clpf_detect_c(const uint8_t *rec, const uint8_t *org, int rstride,
                       int ostride, int x0, int y0, int width, int height,
                       int *sum0, int *sum1, unsigned int strength) {
  int x, y;
  for (y = y0; y < y0 + 8; y++) {
    for (x = x0; x < x0 + 8; x++) {
      int O = org[y * ostride + x];
      int X = rec[y * rstride + x];
      int A = rec[AOMMAX(0, y - 1) * rstride + x];
      int B = rec[y * rstride + AOMMAX(0, x - 2)];
      int C = rec[y * rstride + AOMMAX(0, x - 1)];
      int D = rec[y * rstride + AOMMIN(width - 1, x + 1)];
      int E = rec[y * rstride + AOMMIN(width - 1, x + 2)];
      int F = rec[AOMMIN(height - 1, y + 1) * rstride + x];
      int delta = av1_clpf_sample(X, A, B, C, D, E, F, strength);
      int Y = X + delta;
      *sum0 += (O - X) * (O - X);
      *sum1 += (O - Y) * (O - Y);
    }
  }
}

void aom_clpf_detect_multi_c(const uint8_t *rec, const uint8_t *org,
                             int rstride, int ostride, int x0, int y0,
                             int width, int height, int *sum) {
  int x, y;

  for (y = y0; y < y0 + 8; y++) {
    for (x = x0; x < x0 + 8; x++) {
      int O = org[y * ostride + x];
      int X = rec[y * rstride + x];
      int A = rec[AOMMAX(0, y - 1) * rstride + x];
      int B = rec[y * rstride + AOMMAX(0, x - 2)];
      int C = rec[y * rstride + AOMMAX(0, x - 1)];
      int D = rec[y * rstride + AOMMIN(width - 1, x + 1)];
      int E = rec[y * rstride + AOMMIN(width - 1, x + 2)];
      int F = rec[AOMMIN(height - 1, y + 1) * rstride + x];
      int delta1 = av1_clpf_sample(X, A, B, C, D, E, F, 1);
      int delta2 = av1_clpf_sample(X, A, B, C, D, E, F, 2);
      int delta3 = av1_clpf_sample(X, A, B, C, D, E, F, 4);
      int F1 = X + delta1;
      int F2 = X + delta2;
      int F3 = X + delta3;
      sum[0] += (O - X) * (O - X);
      sum[1] += (O - F1) * (O - F1);
      sum[2] += (O - F2) * (O - F2);
      sum[3] += (O - F3) * (O - F3);
    }
  }
}

int av1_clpf_decision(int k, int l, const YV12_BUFFER_CONFIG *rec,
                      const YV12_BUFFER_CONFIG *org, const AV1_COMMON *cm,
                      int block_size, int w, int h, unsigned int strength,
                      unsigned int fb_size_log2, uint8_t *res) {
  int m, n, sum0 = 0, sum1 = 0;
  for (m = 0; m < h; m++) {
    for (n = 0; n < w; n++) {
      int xpos = (l << fb_size_log2) + n * block_size;
      int ypos = (k << fb_size_log2) + m * block_size;
      const int bs = MAX_MIB_SIZE;
      if (!cm->mi_grid_visible[ypos / bs * cm->mi_stride + xpos / bs]
               ->mbmi.skip)
        aom_clpf_detect(rec->y_buffer, org->y_buffer, rec->y_stride,
                        org->y_stride, xpos, ypos, rec->y_crop_width,
                        rec->y_crop_height, &sum0, &sum1, strength);
    }
  }
  *res = sum1 < sum0;
  return *res;
}

// Calculate the square error of all filter settings.  Result:
// res[0][0]   : unfiltered
// res[0][1-3] : strength=1,2,4, no signals
// res[1][0]   : (bit count, fb size = 128)
// res[1][1-3] : strength=1,2,4, fb size = 128
// res[2][0]   : (bit count, fb size = 64)
// res[2][1-3] : strength=1,2,4, fb size = 64
// res[3][0]   : (bit count, fb size = 32)
// res[3][1-3] : strength=1,2,4, fb size = 32
static int clpf_rdo(int y, int x, const YV12_BUFFER_CONFIG *rec,
                    const YV12_BUFFER_CONFIG *org, const AV1_COMMON *cm,
                    unsigned int block_size, unsigned int fb_size_log2, int w,
                    int h, int64_t res[4][4]) {
  int i, m, n, filtered = 0;
  int sum[4];
  int bslog = get_msb(block_size);
  sum[0] = sum[1] = sum[2] = sum[3] = 0;
  if (fb_size_log2 > (unsigned int)get_msb(MAX_FB_SIZE) - 3) {
    int w1, h1, w2, h2, i, sum1, sum2, sum3, oldfiltered;

    fb_size_log2--;
    w1 = AOMMIN(1 << (fb_size_log2 - bslog), w);
    h1 = AOMMIN(1 << (fb_size_log2 - bslog), h);
    w2 = AOMMIN(w - (1 << (fb_size_log2 - bslog)), w >> 1);
    h2 = AOMMIN(h - (1 << (fb_size_log2 - bslog)), h >> 1);
    i = get_msb(MAX_FB_SIZE) - fb_size_log2;
    sum1 = res[i][1];
    sum2 = res[i][2];
    sum3 = res[i][3];
    oldfiltered = res[i][0];
    res[i][0] = 0;

    filtered =
        clpf_rdo(y, x, rec, org, cm, block_size, fb_size_log2, w1, h1, res);
    if (1 << (fb_size_log2 - bslog) < w)
      filtered |= clpf_rdo(y, x + (1 << fb_size_log2), rec, org, cm, block_size,
                           fb_size_log2, w2, h1, res);
    if (1 << (fb_size_log2 - bslog) < h) {
      filtered |= clpf_rdo(y + (1 << fb_size_log2), x, rec, org, cm, block_size,
                           fb_size_log2, w1, h2, res);
      filtered |= clpf_rdo(y + (1 << fb_size_log2), x + (1 << fb_size_log2),
                           rec, org, cm, block_size, fb_size_log2, w2, h2, res);
    }

    res[i][1] = AOMMIN(sum1 + res[i][0], res[i][1]);
    res[i][2] = AOMMIN(sum2 + res[i][0], res[i][2]);
    res[i][3] = AOMMIN(sum3 + res[i][0], res[i][3]);
    res[i][0] = oldfiltered + filtered;  // Number of signal bits
    return filtered;
  }

  for (m = 0; m < h; m++) {
    for (n = 0; n < w; n++) {
      int xpos = x + n * block_size;
      int ypos = y + m * block_size;
      if (!cm->mi_grid_visible[ypos / MAX_MIB_SIZE * cm->mi_stride +
                               xpos / MAX_MIB_SIZE]
               ->mbmi.skip) {
        aom_clpf_detect_multi(rec->y_buffer, org->y_buffer, rec->y_stride,
                              org->y_stride, xpos, ypos, rec->y_crop_width,
                              rec->y_crop_height, sum);
        filtered = 1;
      }
    }
  }

  for (i = 0; i < 4; i++) {
    res[i][0] += sum[0];
    res[i][1] += sum[1];
    res[i][2] += sum[2];
    res[i][3] += sum[3];
  }
  return filtered;
}

void av1_clpf_test_frame(const YV12_BUFFER_CONFIG *rec,
                         const YV12_BUFFER_CONFIG *org, const AV1_COMMON *cm,
                         int *best_strength, int *best_bs) {
  int i, j, k, l;
  int64_t best, sums[4][4];
  int width = rec->y_crop_width, height = rec->y_crop_height;
  const int bs = MAX_MIB_SIZE;
  int fb_size_log2 = get_msb(MAX_FB_SIZE);
  int num_fb_ver = (height + (1 << fb_size_log2) - bs) >> fb_size_log2;
  int num_fb_hor = (width + (1 << fb_size_log2) - bs) >> fb_size_log2;

  memset(sums, 0, sizeof(sums));

  for (k = 0; k < num_fb_ver; k++) {
    for (l = 0; l < num_fb_hor; l++) {
      // Calculate the block size after frame border clipping
      int h =
          AOMMIN(height, (k + 1) << fb_size_log2) & ((1 << fb_size_log2) - 1);
      int w =
          AOMMIN(width, (l + 1) << fb_size_log2) & ((1 << fb_size_log2) - 1);
      h += !h << fb_size_log2;
      w += !w << fb_size_log2;
      clpf_rdo(k << fb_size_log2, l << fb_size_log2, rec, org, cm, bs,
               fb_size_log2, w / bs, h / bs, sums);
    }
  }

  for (j = 0; j < 4; j++) {
    static const double lambda_square[] = {
      // exp((i - 15.4244) / 8.4010)
      0.159451, 0.179607, 0.202310, 0.227884, 0.256690, 0.289138, 0.325687,
      0.366856, 0.413230, 0.465465, 0.524303, 0.590579, 0.665233, 0.749323,
      0.844044, 0.950737, 1.070917, 1.206289, 1.358774, 1.530533, 1.724004,
      1.941931, 2.187406, 2.463911, 2.775368, 3.126195, 3.521370, 3.966498,
      4.467893, 5.032669, 5.668837, 6.385421, 7.192586, 8.101784, 9.125911,
      10.27949, 11.57890, 13.04256, 14.69124, 16.54832, 18.64016, 20.99641,
      23.65052, 26.64013, 30.00764, 33.80084, 38.07352, 42.88630, 48.30746,
      54.41389, 61.29221, 69.04002, 77.76720, 87.59756, 98.67056, 111.1432,
      125.1926, 141.0179, 158.8436, 178.9227, 201.5399, 227.0160, 255.7126,
      288.0366
    };

    // Estimate the bit costs and adjust the square errors
    double lambda =
        lambda_square[av1_get_qindex(&cm->seg, 0, cm->base_qindex) >> 2];
    int i, cost = (int)((1.2 * lambda * (sums[j][0] + 2 + 2 * (j > 0)) + 0.5));
    for (i = 0; i < 4; i++)
      sums[j][i] = ((sums[j][i] + (i && j) * cost) << 4) + j * 4 + i;
  }

  best = (int64_t)1 << 62;
  for (i = 0; i < 4; i++)
    for (j = 0; j < 4; j++)
      if ((!i || j) && sums[i][j] < best) best = sums[i][j];
  best &= 15;
  *best_bs = (best > 3) * (5 + (best < 12) + (best < 8));
  *best_strength = best ? 1 << ((best - 1) & 3) : 0;
}
