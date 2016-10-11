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
#include <assert.h>
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

void aom_clpf_block_c(const uint8_t *src, uint8_t *dst, int sstride,
                      int dstride, int x0, int y0, int sizex, int sizey,
                      int width, int height, unsigned int strength) {
  int x, y;
  for (y = y0; y < y0 + sizey; y++) {
    for (x = x0; x < x0 + sizex; x++) {
      int X = src[y * sstride + x];
      int A = src[AOMMAX(0, y - 1) * sstride + x];
      int B = src[y * sstride + AOMMAX(0, x - 2)];
      int C = src[y * sstride + AOMMAX(0, x - 1)];
      int D = src[y * sstride + AOMMIN(width - 1, x + 1)];
      int E = src[y * sstride + AOMMIN(width - 1, x + 2)];
      int F = src[AOMMIN(height - 1, y + 1) * sstride + x];
      int delta;
      delta = av1_clpf_sample(X, A, B, C, D, E, F, strength);
      dst[y * dstride + x] = X + delta;
    }
  }
}

#if CONFIG_AOM_HIGHBITDEPTH
// Identical to aom_clpf_block_c() apart from "src" and "dst".
void aom_clpf_block_hbd_c(const uint16_t *src, uint16_t *dst, int sstride,
                          int dstride, int x0, int y0, int sizex, int sizey,
                          int width, int height, unsigned int strength) {
  int x, y;
  for (y = y0; y < y0 + sizey; y++) {
    for (x = x0; x < x0 + sizex; x++) {
      int X = src[y * sstride + x];
      int A = src[AOMMAX(0, y - 1) * sstride + x];
      int B = src[y * sstride + AOMMAX(0, x - 2)];
      int C = src[y * sstride + AOMMAX(0, x - 1)];
      int D = src[y * sstride + AOMMIN(width - 1, x + 1)];
      int E = src[y * sstride + AOMMIN(width - 1, x + 2)];
      int F = src[AOMMIN(height - 1, y + 1) * sstride + x];
      int delta;
      delta = av1_clpf_sample(X, A, B, C, D, E, F, strength);
      dst[y * dstride + x] = X + delta;
    }
  }
}
#endif

// Return number of filtered blocks
int av1_clpf_frame(const YV12_BUFFER_CONFIG *orig_dst,
                   const YV12_BUFFER_CONFIG *rec, const YV12_BUFFER_CONFIG *org,
                   AV1_COMMON *cm, int enable_fb_flag, unsigned int strength,
                   unsigned int fb_size_log2, uint8_t *blocks,
                   int (*decision)(int, int, const YV12_BUFFER_CONFIG *,
                                   const YV12_BUFFER_CONFIG *,
                                   const AV1_COMMON *cm, int, int, int,
                                   unsigned int, unsigned int, uint8_t *)) {
  /* Constrained low-pass filter (CLPF) */
  int c, k, l, m, n;
  const int bs = MI_SIZE;
  const int width = rec->y_crop_width;
  const int height = rec->y_crop_height;
  int xpos, ypos;
  const int sstride = rec->y_stride;
  int dstride = orig_dst->y_stride;
  const int num_fb_hor = (width + (1 << fb_size_log2) - 1) >> fb_size_log2;
  const int num_fb_ver = (height + (1 << fb_size_log2) - 1) >> fb_size_log2;
  int block_index = 0;
  uint8_t *cache = NULL;
  uint8_t **cache_ptr = NULL;
  uint8_t **cache_dst = NULL;
  int cache_idx = 0;
  const int cache_size = num_fb_hor << (2 * fb_size_log2);
  const int cache_blocks = cache_size / (bs * bs);
  YV12_BUFFER_CONFIG dst = *orig_dst;

  assert(bs == 8);  // Optimised code assumes this.

#if CONFIG_AOM_HIGHBITDEPTH
  strength <<= (cm->bit_depth - 8);
#endif

  // Make buffer space for in-place filtering
  if (rec->y_buffer == dst.y_buffer) {
#if CONFIG_AOM_HIGHBITDEPTH
    CHECK_MEM_ERROR(cm, cache,
                    aom_malloc(cache_size << !!cm->use_highbitdepth));
    dst.y_buffer = cm->use_highbitdepth ? CONVERT_TO_BYTEPTR(cache) : cache;
#else
    CHECK_MEM_ERROR(cm, cache, aom_malloc(cache_size));
    dst.y_buffer = cache;
#endif
    CHECK_MEM_ERROR(cm, cache_ptr,
                    aom_malloc(cache_blocks * sizeof(*cache_ptr)));
    CHECK_MEM_ERROR(cm, cache_dst,
                    aom_malloc(cache_blocks * sizeof(*cache_dst)));
    memset(cache_ptr, 0, cache_blocks * sizeof(*cache_dst));
    dstride = bs;
  }

  // Iterate over all filter blocks
  for (k = 0; k < num_fb_ver; k++) {
    for (l = 0; l < num_fb_hor; l++) {
      int h, w;
      int allskip = 1;
      const int xoff = l << fb_size_log2;
      const int yoff = k << fb_size_log2;
      for (m = 0; allskip && m < (1 << fb_size_log2) / bs; m++) {
        for (n = 0; allskip && n < (1 << fb_size_log2) / bs; n++) {
          xpos = xoff + n * bs;
          ypos = yoff + m * bs;
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
            xpos = xoff + n * bs;
            ypos = yoff + m * bs;
            if (!cm->mi_grid_visible[ypos / bs * cm->mi_stride + xpos / bs]
                     ->mbmi.skip) {  // Not skip block
              // Temporary buffering needed if filtering in-place
              if (cache) {
                if (cache_ptr[cache_idx]) {
// Copy filtered block back into the frame
#if CONFIG_AOM_HIGHBITDEPTH
                  if (cm->use_highbitdepth) {
                    uint16_t *const d =
                        CONVERT_TO_SHORTPTR(cache_dst[cache_idx]);
                    for (c = 0; c < bs; c++) {
                      *(uint64_t *)(d + c * sstride) =
                          *(uint64_t *)(cache_ptr[cache_idx] + c * bs * 2);
                      *(uint64_t *)(d + c * sstride + 4) =
                          *(uint64_t *)(cache_ptr[cache_idx] + c * bs * 2 + 8);
                    }
                  } else {
                    for (c = 0; c < bs; c++)
                      *(uint64_t *)(cache_dst[cache_idx] + c * sstride) =
                          *(uint64_t *)(cache_ptr[cache_idx] + c * bs);
                  }
#else
                  for (c = 0; c < bs; c++)
                    *(uint64_t *)(cache_dst[cache_idx] + c * sstride) =
                        *(uint64_t *)(cache_ptr[cache_idx] + c * bs);
#endif
                }
#if CONFIG_AOM_HIGHBITDEPTH
                if (cm->use_highbitdepth) {
                  cache_ptr[cache_idx] = cache + cache_idx * bs * bs * 2;
                  dst.y_buffer = CONVERT_TO_BYTEPTR(cache_ptr[cache_idx]) -
                                 ypos * bs - xpos;
                } else {
                  cache_ptr[cache_idx] = cache + cache_idx * bs * bs;
                  dst.y_buffer = cache_ptr[cache_idx] - ypos * bs - xpos;
                }
#else
                cache_ptr[cache_idx] = cache + cache_idx * bs * bs;
                dst.y_buffer = cache_ptr[cache_idx] - ypos * bs - xpos;
#endif
                cache_dst[cache_idx] = rec->y_buffer + ypos * sstride + xpos;
                if (++cache_idx >= cache_blocks) cache_idx = 0;
              }

// Apply the filter
#if CONFIG_AOM_HIGHBITDEPTH
              if (cm->use_highbitdepth) {
                aom_clpf_block_hbd(CONVERT_TO_SHORTPTR(rec->y_buffer),
                                   CONVERT_TO_SHORTPTR(dst.y_buffer), sstride,
                                   dstride, xpos, ypos, bs, bs, width, height,
                                   strength);
              } else {
                aom_clpf_block(rec->y_buffer, dst.y_buffer, sstride, dstride,
                               xpos, ypos, bs, bs, width, height, strength);
              }
#else
              aom_clpf_block(rec->y_buffer, dst.y_buffer, sstride, dstride,
                             xpos, ypos, bs, bs, width, height, strength);
#endif
            } else {  // Skip block, copy instead
              if (!cache) {
#if CONFIG_AOM_HIGHBITDEPTH
                if (cm->use_highbitdepth) {
                  uint16_t *const d = CONVERT_TO_SHORTPTR(dst.y_buffer);
                  const uint16_t *const s = CONVERT_TO_SHORTPTR(rec->y_buffer);
                  for (c = 0; c < bs; c++) {
                    *(uint64_t *)(d + (ypos + c) * dstride + xpos) =
                        *(uint64_t *)(s + (ypos + c) * sstride + xpos);
                    *(uint64_t *)(d + (ypos + c) * dstride + xpos + 4) =
                        *(uint64_t *)(s + (ypos + c) * sstride + xpos + 4);
                  }
                } else {
                  for (c = 0; c < bs; c++)
                    *(uint64_t *)(dst.y_buffer + (ypos + c) * dstride + xpos) =
                        *(uint64_t *)(rec->y_buffer + (ypos + c) * sstride +
                                      xpos);
                }
#else
                for (c = 0; c < bs; c++)
                  *(uint64_t *)(dst.y_buffer + (ypos + c) * dstride + xpos) = *(
                      uint64_t *)(rec->y_buffer + (ypos + c) * sstride + xpos);
#endif
              }
            }
          }
        }
      } else {  // Entire filter block is skip, copy
        if (!cache) {
#if CONFIG_AOM_HIGHBITDEPTH
          if (cm->use_highbitdepth) {
            for (m = 0; m < h; m++)
              memcpy(CONVERT_TO_SHORTPTR(dst.y_buffer) + (yoff + m) * dstride +
                         xoff,
                     CONVERT_TO_SHORTPTR(rec->y_buffer) + (yoff + m) * sstride +
                         xoff,
                     w * 2);
          } else {
            for (m = 0; m < h; m++)
              memcpy(dst.y_buffer + (yoff + m) * dstride + xoff,
                     rec->y_buffer + (yoff + m) * sstride + xoff, w);
          }
#else
          for (m = 0; m < h; m++)
            memcpy(dst.y_buffer + (yoff + m) * dstride + xoff,
                   rec->y_buffer + (yoff + m) * sstride + xoff, w);
#endif
        }
      }
      block_index += !allskip;  // Count number of blocks filtered
    }
  }

  if (cache) {
    // Copy remaining blocks into the frame
    for (cache_idx = 0; cache_idx < cache_blocks && cache_ptr[cache_idx];
         cache_idx++) {
#if CONFIG_AOM_HIGHBITDEPTH
      if (cm->use_highbitdepth) {
        uint16_t *const d = CONVERT_TO_SHORTPTR(cache_dst[cache_idx]);
        for (c = 0; c < bs; c++) {
          *(uint64_t *)(d + c * sstride) =
              *(uint64_t *)(cache_ptr[cache_idx] + c * bs * 2);
          *(uint64_t *)(d + c * sstride + 4) =
              *(uint64_t *)(cache_ptr[cache_idx] + c * bs * 2 + 8);
        }
      } else {
        for (c = 0; c < bs; c++)
          *(uint64_t *)(cache_dst[cache_idx] + c * sstride) =
              *(uint64_t *)(cache_ptr[cache_idx] + c * bs);
      }
#else
      for (c = 0; c < bs; c++)
        *(uint64_t *)(cache_dst[cache_idx] + c * sstride) =
            *(uint64_t *)(cache_ptr[cache_idx] + c * bs);
#endif
    }

    aom_free(cache);
    aom_free(cache_ptr);
    aom_free(cache_dst);
  }

  return block_index;
}
