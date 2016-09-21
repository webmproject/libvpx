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
#include "aom/aom_image.h"
#include "aom/aom_integer.h"
#include "av1/common/quant_common.h"

// Calculate the error of a filtered and unfiltered block
void aom_clpf_detect_c(const uint8_t *rec, const uint8_t *org, int rstride,
                       int ostride, int x0, int y0, int width, int height,
                       int *sum0, int *sum1, unsigned int strength, int size) {
  int x, y;
  for (y = y0; y < y0 + size; y++) {
    for (x = x0; x < x0 + size; x++) {
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
                             int width, int height, int *sum, int size) {
  int x, y;

  for (y = y0; y < y0 + size; y++) {
    for (x = x0; x < x0 + size; x++) {
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

#if CONFIG_AOM_HIGHBITDEPTH
// Identical to aom_clpf_detect_c() apart from "rec" and "org".
void aom_clpf_detect_hbd_c(const uint16_t *rec, const uint16_t *org,
                           int rstride, int ostride, int x0, int y0, int width,
                           int height, int *sum0, int *sum1,
                           unsigned int strength, int shift, int size) {
  int x, y;
  for (y = y0; y < y0 + size; y++) {
    for (x = x0; x < x0 + size; x++) {
      int O = org[y * ostride + x] >> shift;
      int X = rec[y * rstride + x] >> shift;
      int A = rec[AOMMAX(0, y - 1) * rstride + x] >> shift;
      int B = rec[y * rstride + AOMMAX(0, x - 2)] >> shift;
      int C = rec[y * rstride + AOMMAX(0, x - 1)] >> shift;
      int D = rec[y * rstride + AOMMIN(width - 1, x + 1)] >> shift;
      int E = rec[y * rstride + AOMMIN(width - 1, x + 2)] >> shift;
      int F = rec[AOMMIN(height - 1, y + 1) * rstride + x] >> shift;
      int delta = av1_clpf_sample(X, A, B, C, D, E, F, strength >> shift);
      int Y = X + delta;
      *sum0 += (O - X) * (O - X);
      *sum1 += (O - Y) * (O - Y);
    }
  }
}

// aom_clpf_detect_multi_c() apart from "rec" and "org".
void aom_clpf_detect_multi_hbd_c(const uint16_t *rec, const uint16_t *org,
                                 int rstride, int ostride, int x0, int y0,
                                 int width, int height, int *sum, int shift,
                                 int size) {
  int x, y;

  for (y = y0; y < y0 + size; y++) {
    for (x = x0; x < x0 + size; x++) {
      int O = org[y * ostride + x] >> shift;
      int X = rec[y * rstride + x] >> shift;
      int A = rec[AOMMAX(0, y - 1) * rstride + x] >> shift;
      int B = rec[y * rstride + AOMMAX(0, x - 2)] >> shift;
      int C = rec[y * rstride + AOMMAX(0, x - 1)] >> shift;
      int D = rec[y * rstride + AOMMIN(width - 1, x + 1)] >> shift;
      int E = rec[y * rstride + AOMMIN(width - 1, x + 2)] >> shift;
      int F = rec[AOMMIN(height - 1, y + 1) * rstride + x] >> shift;
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
#endif

int av1_clpf_decision(int k, int l, const YV12_BUFFER_CONFIG *rec,
                      const YV12_BUFFER_CONFIG *org, const AV1_COMMON *cm,
                      int block_size, int w, int h, unsigned int strength,
                      unsigned int fb_size_log2, int8_t *res) {
  int m, n, sum0 = 0, sum1 = 0;

  for (m = 0; m < h; m++) {
    for (n = 0; n < w; n++) {
      int xpos = (l << fb_size_log2) + n * block_size;
      int ypos = (k << fb_size_log2) + m * block_size;
      if (fb_size_log2 == MAX_FB_SIZE_LOG2 ||
          !cm->mi_grid_visible[ypos / MI_SIZE * cm->mi_stride + xpos / MI_SIZE]
               ->mbmi.skip) {
#if CONFIG_AOM_HIGHBITDEPTH
        if (cm->use_highbitdepth) {
          aom_clpf_detect_hbd(CONVERT_TO_SHORTPTR(rec->y_buffer),
                              CONVERT_TO_SHORTPTR(org->y_buffer), rec->y_stride,
                              org->y_stride, xpos, ypos, rec->y_crop_width,
                              rec->y_crop_height, &sum0, &sum1, strength,
                              cm->bit_depth - 8, block_size);
        } else {
          aom_clpf_detect(rec->y_buffer, org->y_buffer, rec->y_stride,
                          org->y_stride, xpos, ypos, rec->y_crop_width,
                          rec->y_crop_height, &sum0, &sum1, strength,
                          block_size);
        }
#else
        aom_clpf_detect(rec->y_buffer, org->y_buffer, rec->y_stride,
                        org->y_stride, xpos, ypos, rec->y_crop_width,
                        rec->y_crop_height, &sum0, &sum1, strength, block_size);
#endif
      }
    }
  }
  *res = sum1 < sum0;
  return *res;
}

// Calculate the square error of all filter settings.  Result:
// res[0][0]   : unfiltered
// res[0][1-3] : strength=1,2,4, no signals
// (Only for luma:)
// res[1][0]   : (bit count, fb size = 128)
// res[1][1-3] : strength=1,2,4, fb size = 128
// res[1][4]   : unfiltered, including skip
// res[1][5-7] : strength=1,2,4, including skip, fb_size = 128
// res[2][0]   : (bit count, fb size = 64)
// res[2][1-3] : strength=1,2,4, fb size = 64
// res[3][0]   : (bit count, fb size = 32)
// res[3][1-3] : strength=1,2,4, fb size = 32
static int clpf_rdo(int y, int x, const YV12_BUFFER_CONFIG *rec,
                    const YV12_BUFFER_CONFIG *org, const AV1_COMMON *cm,
                    unsigned int block_size, unsigned int fb_size_log2, int w,
                    int h, int64_t res[4][8], int plane) {
  int c, m, n, filtered = 0;
  int sum[8];
  const int subx = plane != AOM_PLANE_Y && rec->subsampling_x;
  const int suby = plane != AOM_PLANE_Y && rec->subsampling_y;
  int bslog = get_msb(block_size);
  uint8_t *rec_buffer =
      plane != AOM_PLANE_Y
          ? (plane == AOM_PLANE_U ? rec->u_buffer : rec->v_buffer)
          : rec->y_buffer;
  uint8_t *org_buffer =
      plane != AOM_PLANE_Y
          ? (plane == AOM_PLANE_U ? org->u_buffer : org->v_buffer)
          : org->y_buffer;
  int rec_width = plane != AOM_PLANE_Y ? rec->uv_crop_width : rec->y_crop_width;
  int rec_height =
      plane != AOM_PLANE_Y ? rec->uv_crop_height : rec->y_crop_height;
  int rec_stride = plane != AOM_PLANE_Y ? rec->uv_stride : rec->y_stride;
  int org_stride = plane != AOM_PLANE_Y ? org->uv_stride : org->y_stride;
  sum[0] = sum[1] = sum[2] = sum[3] = sum[4] = sum[5] = sum[6] = sum[7] = 0;
  if (plane == AOM_PLANE_Y &&
      fb_size_log2 > (unsigned int)get_msb(MAX_FB_SIZE) - 3) {
    int w1, h1, w2, h2, i, sum1, sum2, sum3, oldfiltered;

    filtered = fb_size_log2-- == MAX_FB_SIZE_LOG2;
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

    filtered |= clpf_rdo(y, x, rec, org, cm, block_size, fb_size_log2, w1, h1,
                         res, plane);
    if (1 << (fb_size_log2 - bslog) < w)
      filtered |= clpf_rdo(y, x + (1 << fb_size_log2), rec, org, cm, block_size,
                           fb_size_log2, w2, h1, res, plane);
    if (1 << (fb_size_log2 - bslog) < h) {
      filtered |= clpf_rdo(y + (1 << fb_size_log2), x, rec, org, cm, block_size,
                           fb_size_log2, w1, h2, res, plane);
      filtered |=
          clpf_rdo(y + (1 << fb_size_log2), x + (1 << fb_size_log2), rec, org,
                   cm, block_size, fb_size_log2, w2, h2, res, plane);
    }

    // Correct sums for unfiltered blocks
    res[i][1] = AOMMIN(sum1 + res[i][0], res[i][1]);
    res[i][2] = AOMMIN(sum2 + res[i][0], res[i][2]);
    res[i][3] = AOMMIN(sum3 + res[i][0], res[i][3]);
    if (i == 1) {
      res[i][5] = AOMMIN(sum1 + res[i][4], res[i][5]);
      res[i][6] = AOMMIN(sum2 + res[i][4], res[i][6]);
      res[i][7] = AOMMIN(sum3 + res[i][4], res[i][7]);
    }

    res[i][0] = oldfiltered + filtered;  // Number of signal bits

    return filtered;
  }

  for (m = 0; m < h; m++) {
    for (n = 0; n < w; n++) {
      int xpos = x + n * block_size;
      int ypos = y + m * block_size;
      int skip =  // Filtered skip blocks stored only for fb_size == 128
          4 *
          !!cm->mi_grid_visible[(ypos << suby) / MI_SIZE * cm->mi_stride +
                                (xpos << subx) / MI_SIZE]
                ->mbmi.skip;
#if CONFIG_AOM_HIGHBITDEPTH
      if (cm->use_highbitdepth) {
        aom_clpf_detect_multi_hbd(CONVERT_TO_SHORTPTR(rec_buffer),
                                  CONVERT_TO_SHORTPTR(org_buffer), rec_stride,
                                  org_stride, xpos, ypos, rec_width, rec_height,
                                  sum + skip, cm->bit_depth - 8, block_size);
      } else {
        aom_clpf_detect_multi(rec_buffer, org_buffer, rec_stride, org_stride,
                              xpos, ypos, rec_width, rec_height, sum + skip,
                              block_size);
      }
#else
      aom_clpf_detect_multi(rec_buffer, org_buffer, rec_stride, org_stride,
                            xpos, ypos, rec_width, rec_height, sum + skip,
                            block_size);
#endif
      filtered |= !skip;
    }
  }

  for (c = 0; c < (plane == AOM_PLANE_Y ? 4 : 1); c++) {
    res[c][0] += sum[0];
    res[c][1] += sum[1];
    res[c][2] += sum[2];
    res[c][3] += sum[3];
    if (c != 1) continue;
    // Only needed when fb_size == 128
    res[c][4] += sum[4];
    res[c][5] += sum[5];
    res[c][6] += sum[6];
    res[c][7] += sum[7];
  }
  return filtered;
}

void av1_clpf_test_frame(const YV12_BUFFER_CONFIG *rec,
                         const YV12_BUFFER_CONFIG *org, const AV1_COMMON *cm,
                         int *best_strength, int *best_bs, int plane) {
  int c, j, k, l;
  int64_t best, sums[4][8];
  int width = plane != AOM_PLANE_Y ? rec->uv_crop_width : rec->y_crop_width;
  int height = plane != AOM_PLANE_Y ? rec->uv_crop_height : rec->y_crop_height;
  const int bs = MI_SIZE;
  const int bslog = get_msb(bs);
  int fb_size_log2 = get_msb(MAX_FB_SIZE);
  int num_fb_ver = (height + (1 << fb_size_log2) - bs) >> fb_size_log2;
  int num_fb_hor = (width + (1 << fb_size_log2) - bs) >> fb_size_log2;

  memset(sums, 0, sizeof(sums));

  if (plane != AOM_PLANE_Y)
    // Use a block size of MI_SIZE regardless of the subsampling.  This
    // This is accurate enough to determine the best strength and
    // we don't need to add SIMD optimisations for 4x4 blocks.
    clpf_rdo(0, 0, rec, org, cm, bs, fb_size_log2, width >> bslog,
             height >> bslog, sums, plane);
  else
    for (k = 0; k < num_fb_ver; k++) {
      for (l = 0; l < num_fb_hor; l++) {
        // Calculate the block size after frame border clipping
        int h =
            AOMMIN(height, (k + 1) << fb_size_log2) & ((1 << fb_size_log2) - 1);
        int w =
            AOMMIN(width, (l + 1) << fb_size_log2) & ((1 << fb_size_log2) - 1);
        h += !h << fb_size_log2;
        w += !w << fb_size_log2;
        clpf_rdo(k << fb_size_log2, l << fb_size_log2, rec, org, cm, MI_SIZE,
                 fb_size_log2, w >> bslog, h >> bslog, sums, plane);
      }
    }

  // For fb_size == 128 skip blocks are included in the result.
  if (plane == AOM_PLANE_Y) {
    sums[1][1] += sums[1][5] - sums[1][4];
    sums[1][2] += sums[1][6] - sums[1][4];
    sums[1][3] += sums[1][7] - sums[1][4];
  } else {  // Slightly favour unfiltered chroma
    sums[0][0] -= sums[0][0] >> 7;
  }

  for (j = 0; j < 4; j++) {
    static const double lambda_square[] = {
      // exp(x / 8.5)
      1.0000, 1.1248, 1.2653, 1.4232, 1.6009, 1.8008, 2.0256, 2.2785,
      2.5630, 2.8830, 3.2429, 3.6478, 4.1032, 4.6155, 5.1917, 5.8399,
      6.5689, 7.3891, 8.3116, 9.3492, 10.516, 11.829, 13.306, 14.967,
      16.836, 18.938, 21.302, 23.962, 26.953, 30.318, 34.103, 38.361,
      43.151, 48.538, 54.598, 61.414, 69.082, 77.706, 87.408, 98.320,
      110.59, 124.40, 139.93, 157.40, 177.05, 199.16, 224.02, 251.99,
      283.45, 318.84, 358.65, 403.42, 453.79, 510.45, 574.17, 645.86,
      726.49, 817.19, 919.22, 1033.9, 1163.0, 1308.2, 1471.6, 1655.3
    };

    // Estimate the bit costs and adjust the square errors
    double lambda =
        lambda_square[av1_get_qindex(&cm->seg, 0, cm->base_qindex) >> 2];
    int i, cost = (int)((lambda * (sums[j][0] + 6 + 2 * (j > 0)) + 0.5));
    for (i = 0; i < 4; i++)
      sums[j][i] = ((sums[j][i] + (i && j) * cost) << 4) + j * 4 + i;
  }

  best = (int64_t)1 << 62;
  for (c = 0; c < (plane == AOM_PLANE_Y ? 4 : 1); c++)
    for (j = 0; j < 4; j++)
      if ((!c || j) && sums[c][j] < best) best = sums[c][j];
  best &= 15;
  if (best_bs) *best_bs = (best > 3) * (5 + (best < 12) + (best < 8));
  *best_strength = best ? 1 << ((best - 1) & 3) : 0;
}
