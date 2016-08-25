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
#include "aom_dsp/aom_dsp_common.h"

int av1_clpf_maxbits(const AV1_COMMON *cm) {
  return get_msb(
             ALIGN_POWER_OF_TWO(cm->mi_cols * MI_SIZE, cm->clpf_size + 4) *
                 ALIGN_POWER_OF_TWO(cm->mi_rows * MI_SIZE, cm->clpf_size + 4) >>
             (cm->clpf_size * 2 + 8)) +
         1;
}

int av1_clpf_sample(int X, int A, int B, int C, int D, int E, int F, int b) {
  int delta = 4 * clamp(A - X, -b, b) + clamp(B - X, -b, b) +
              3 * clamp(C - X, -b, b) + 3 * clamp(D - X, -b, b) +
              clamp(E - X, -b, b) + 4 * clamp(F - X, -b, b);
  return (8 + delta - (delta < 0)) >> 4;
}

void aom_clpf_block_c(const uint8_t *src, uint8_t *dst, int stride, int x0,
                      int y0, int sizex, int sizey, int width, int height,
                      unsigned int strength) {
  int x, y;
  for (y = y0; y < y0 + sizey; y++) {
    for (x = x0; x < x0 + sizex; x++) {
      int X = src[y * stride + x];
      int A = src[AOMMAX(0, y - 1) * stride + x];
      int B = src[y * stride + AOMMAX(0, x - 2)];
      int C = src[y * stride + AOMMAX(0, x - 1)];
      int D = src[y * stride + AOMMIN(width - 1, x + 1)];
      int E = src[y * stride + AOMMIN(width - 1, x + 2)];
      int F = src[AOMMIN(height - 1, y + 1) * stride + x];
      int delta;
      delta = av1_clpf_sample(X, A, B, C, D, E, F, strength);
      dst[y * stride + x] = X + delta;
    }
  }
}

// Return number of filtered blocks
int av1_clpf_frame(const YV12_BUFFER_CONFIG *dst, const YV12_BUFFER_CONFIG *rec,
                   const YV12_BUFFER_CONFIG *org, const AV1_COMMON *cm,
                   int enable_fb_flag, unsigned int strength,
                   unsigned int fb_size_log2, uint8_t *blocks,
                   int (*decision)(int, int, const YV12_BUFFER_CONFIG *,
                                   const YV12_BUFFER_CONFIG *,
                                   const AV1_COMMON *cm, int, int, int,
                                   unsigned int, unsigned int, uint8_t *)) {
  /* Constrained low-pass filter (CLPF) */
  int c, k, l, m, n;
  const int bs = MI_SIZE;
  int width = cm->mi_cols * bs;
  int height = cm->mi_rows * bs;
  int xpos, ypos;
  int stride_y = rec->y_stride;
  int num_fb_hor = (width + (1 << fb_size_log2) - bs) >> fb_size_log2;
  int num_fb_ver = (height + (1 << fb_size_log2) - bs) >> fb_size_log2;
  int block_index = 0;

  // Iterate over all filter blocks
  for (k = 0; k < num_fb_ver; k++) {
    for (l = 0; l < num_fb_hor; l++) {
      int h, w;
      int allskip = 1;
      for (m = 0; allskip && m < (1 << fb_size_log2) / bs; m++) {
        for (n = 0; allskip && n < (1 << fb_size_log2) / bs; n++) {
          xpos = (l << fb_size_log2) + n * bs;
          ypos = (k << fb_size_log2) + m * bs;
          if (xpos < width && ypos < height) {
            allskip &=
                cm->mi_grid_visible[ypos / bs * cm->mi_stride + xpos / bs]
                    ->mbmi.skip;
          }
        }
      }

      // Calculate the actual filter block size near frame edges
      h = AOMMIN(height, (k + 1) << fb_size_log2) & ((1 << fb_size_log2) - 1);
      w = AOMMIN(width, (l + 1) << fb_size_log2) & ((1 << fb_size_log2) - 1);
      h += !h << fb_size_log2;
      w += !w << fb_size_log2;
      if (!allskip &&  // Do not filter the block if all is skip encoded
          (!enable_fb_flag ||
           decision(k, l, rec, org, cm, bs, w / bs, h / bs, strength,
                    fb_size_log2, blocks + block_index))) {
        // Iterate over all smaller blocks inside the filter block
        for (m = 0; m < (h + bs - 1) / bs; m++) {
          for (n = 0; n < (w + bs - 1) / bs; n++) {
            xpos = (l << fb_size_log2) + n * bs;
            ypos = (k << fb_size_log2) + m * bs;
            if (!cm->mi_grid_visible[ypos / bs * cm->mi_stride + xpos / bs]
                     ->mbmi.skip) {
              // Not skip block, apply the filter
              aom_clpf_block(rec->y_buffer, dst->y_buffer, stride_y, xpos, ypos,
                             bs, bs, width, height, strength);
            } else {  // Skip block, copy instead
              for (c = 0; c < bs; c++)
                *(uint64_t *)(dst->y_buffer + (ypos + c) * stride_y + xpos) =
                    *(uint64_t *)(rec->y_buffer + (ypos + c) * stride_y + xpos);
            }
          }
        }
      } else {  // Entire filter block is skip, copy
        for (m = 0; m < h; m++)
          memcpy(dst->y_buffer + ((k << fb_size_log2) + m) * stride_y +
                     (l << fb_size_log2),
                 rec->y_buffer + ((k << fb_size_log2) + m) * stride_y +
                     (l << fb_size_log2),
                 w);
      }
      block_index += !allskip;  // Count number of blocks filtered
    }
  }

  return block_index;
}
