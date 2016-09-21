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
#include "aom_dsp/aom_dsp_common.h"

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
void av1_clpf_frame(const YV12_BUFFER_CONFIG *frame,
                    const YV12_BUFFER_CONFIG *org, AV1_COMMON *cm,
                    int enable_fb_flag, unsigned int strength,
                    unsigned int fb_size_log2, int plane,
                    int (*decision)(int, int, const YV12_BUFFER_CONFIG *,
                                    const YV12_BUFFER_CONFIG *,
                                    const AV1_COMMON *cm, int, int, int,
                                    unsigned int, unsigned int, int8_t *)) {
  /* Constrained low-pass filter (CLPF) */
  int c, k, l, m, n;
  const int subx = plane != AOM_PLANE_Y && frame->subsampling_x;
  const int suby = plane != AOM_PLANE_Y && frame->subsampling_y;
  const int bs = (subx || suby) ? 4 : 8;
  const int bslog = get_msb(bs);
  int width = plane != AOM_PLANE_Y ? frame->uv_crop_width : frame->y_crop_width;
  int height =
      plane != AOM_PLANE_Y ? frame->uv_crop_height : frame->y_crop_height;
  int xpos, ypos;
  const int sstride = plane != AOM_PLANE_Y ? frame->uv_stride : frame->y_stride;
  int dstride = bs;
  const int num_fb_hor = (width + (1 << fb_size_log2) - 1) >> fb_size_log2;
  const int num_fb_ver = (height + (1 << fb_size_log2) - 1) >> fb_size_log2;
  uint8_t *cache = NULL;
  uint8_t **cache_ptr = NULL;
  uint8_t **cache_dst = NULL;
  int cache_idx = 0;
  const int cache_size = num_fb_hor << (2 * fb_size_log2);
  const int cache_blocks = cache_size / (bs * bs);
  uint8_t *src_buffer =
      plane != AOM_PLANE_Y
          ? (plane == AOM_PLANE_U ? frame->u_buffer : frame->v_buffer)
          : frame->y_buffer;
  uint8_t *dst_buffer;

// Make buffer space for in-place filtering
#if CONFIG_AOM_HIGHBITDEPTH
  strength <<= (cm->bit_depth - 8);
  CHECK_MEM_ERROR(cm, cache, aom_malloc(cache_size << !!cm->use_highbitdepth));
  dst_buffer = cm->use_highbitdepth ? CONVERT_TO_BYTEPTR(cache) : cache;
#else
  CHECK_MEM_ERROR(cm, cache, aom_malloc(cache_size));
  dst_buffer = cache;
#endif
  CHECK_MEM_ERROR(cm, cache_ptr, aom_malloc(cache_blocks * sizeof(*cache_ptr)));
  CHECK_MEM_ERROR(cm, cache_dst, aom_malloc(cache_blocks * sizeof(*cache_dst)));
  memset(cache_ptr, 0, cache_blocks * sizeof(*cache_dst));

  // Iterate over all filter blocks
  for (k = 0; k < num_fb_ver; k++) {
    for (l = 0; l < num_fb_hor; l++) {
      int h, w;
      int allskip = !(enable_fb_flag && fb_size_log2 == MAX_FB_SIZE_LOG2);
      const int xoff = l << fb_size_log2;
      const int yoff = k << fb_size_log2;
      for (m = 0; allskip && m < (1 << fb_size_log2) / bs; m++) {
        for (n = 0; allskip && n < (1 << fb_size_log2) / bs; n++) {
          xpos = xoff + n * bs;
          ypos = yoff + m * bs;
          if (xpos < width && ypos < height) {
            allskip &=
                cm->mi_grid_visible[(ypos << suby) / MI_SIZE * cm->mi_stride +
                                    (xpos << subx) / MI_SIZE]
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
           // Only called if fb_flag enabled (luma only)
           decision(k, l, frame, org, cm, bs, w / bs, h / bs, strength,
                    fb_size_log2,
                    cm->clpf_blocks + yoff / MIN_FB_SIZE * cm->clpf_stride +
                        xoff / MIN_FB_SIZE))) {
        // Iterate over all smaller blocks inside the filter block
        for (m = 0; m < ((h + bs - 1) >> bslog); m++) {
          for (n = 0; n < ((w + bs - 1) >> bslog); n++) {
            int sizex, sizey;
            xpos = xoff + n * bs;
            ypos = yoff + m * bs;
            sizex = AOMMIN(width - xpos, bs);
            sizey = AOMMIN(height - ypos, bs);
            if (!cm->mi_grid_visible[(ypos << suby) / MI_SIZE * cm->mi_stride +
                                     (xpos << subx) / MI_SIZE]
                     ->mbmi.skip ||
                (enable_fb_flag && fb_size_log2 == MAX_FB_SIZE_LOG2)) {
              // Temporary buffering needed for in-place filtering
              if (cache_ptr[cache_idx]) {
// Copy filtered block back into the frame
#if CONFIG_AOM_HIGHBITDEPTH
                if (cm->use_highbitdepth) {
                  uint16_t *const d = CONVERT_TO_SHORTPTR(cache_dst[cache_idx]);
                  if (sizex == 8) {
                    for (c = 0; c < sizey; c++) {
                      *(uint64_t *)(d + c * sstride) =
                          *(uint64_t *)(cache_ptr[cache_idx] + c * bs * 2);
                      *(uint64_t *)(d + c * sstride + 4) =
                          *(uint64_t *)(cache_ptr[cache_idx] + c * bs * 2 + 8);
                    }
                  } else if (sizex == 4) {
                    for (c = 0; c < sizey; c++)
                      *(uint64_t *)(d + c * sstride) =
                          *(uint64_t *)(cache_ptr[cache_idx] + c * bs * 2);
                  } else {
                    for (c = 0; c < sizey; c++)
                      memcpy(d + c * sstride, cache_ptr[cache_idx] + c * bs * 2,
                             sizex);
                  }
                } else {
                  if (sizex == 8)
                    for (c = 0; c < sizey; c++)
                      *(uint64_t *)(cache_dst[cache_idx] + c * sstride) =
                          *(uint64_t *)(cache_ptr[cache_idx] + c * bs);
                  else if (sizex == 4)
                    for (c = 0; c < sizey; c++)
                      *(uint32_t *)(cache_dst[cache_idx] + c * sstride) =
                          *(uint32_t *)(cache_ptr[cache_idx] + c * bs);
                  else
                    for (c = 0; c < sizey; c++)
                      memcpy(cache_dst[cache_idx] + c * sstride,
                             cache_ptr[cache_idx] + c * bs, sizex);
                }
#else
                if (sizex == 8)
                  for (c = 0; c < sizey; c++)
                    *(uint64_t *)(cache_dst[cache_idx] + c * sstride) =
                        *(uint64_t *)(cache_ptr[cache_idx] + c * bs);
                else if (sizex == 4)
                  for (c = 0; c < sizey; c++)
                    *(uint32_t *)(cache_dst[cache_idx] + c * sstride) =
                        *(uint32_t *)(cache_ptr[cache_idx] + c * bs);
                else
                  for (c = 0; c < sizey; c++)
                    memcpy(cache_dst[cache_idx] + c * sstride,
                           cache_ptr[cache_idx] + c * bs, sizex);
#endif
              }
#if CONFIG_AOM_HIGHBITDEPTH
              if (cm->use_highbitdepth) {
                cache_ptr[cache_idx] = cache + cache_idx * bs * bs * 2;
                dst_buffer =
                    CONVERT_TO_BYTEPTR(cache_ptr[cache_idx]) - ypos * bs - xpos;
              } else {
                cache_ptr[cache_idx] = cache + cache_idx * bs * bs;
                dst_buffer = cache_ptr[cache_idx] - ypos * bs - xpos;
              }
#else
              cache_ptr[cache_idx] = cache + cache_idx * bs * bs;
              dst_buffer = cache_ptr[cache_idx] - ypos * bs - xpos;
#endif
              cache_dst[cache_idx] = src_buffer + ypos * sstride + xpos;
              if (++cache_idx >= cache_blocks) cache_idx = 0;

// Apply the filter
#if CONFIG_AOM_HIGHBITDEPTH
              if (cm->use_highbitdepth) {
                aom_clpf_block_hbd(CONVERT_TO_SHORTPTR(src_buffer),
                                   CONVERT_TO_SHORTPTR(dst_buffer), sstride,
                                   dstride, xpos, ypos, sizex, sizey, width,
                                   height, strength);
              } else {
                aom_clpf_block(src_buffer, dst_buffer, sstride, dstride, xpos,
                               ypos, sizex, sizey, width, height, strength);
              }
#else
              aom_clpf_block(src_buffer, dst_buffer, sstride, dstride, xpos,
                             ypos, sizex, sizey, width, height, strength);
#endif
            }
          }
        }
      }
    }
  }

  // Copy remaining blocks into the frame
  for (cache_idx = 0; cache_idx < cache_blocks && cache_ptr[cache_idx];
       cache_idx++) {
#if CONFIG_AOM_HIGHBITDEPTH
    if (cm->use_highbitdepth) {
      uint16_t *const d = CONVERT_TO_SHORTPTR(cache_dst[cache_idx]);
      for (c = 0; c < bs; c++) {
        *(uint64_t *)(d + c * sstride) =
            *(uint64_t *)(cache_ptr[cache_idx] + c * bs * 2);
        if (bs == 8)
          *(uint64_t *)(d + c * sstride + 4) =
              *(uint64_t *)(cache_ptr[cache_idx] + c * bs * 2 + 8);
      }
    } else {
      for (c = 0; c < bs; c++)
        if (bs == 4)
          *(uint32_t *)(cache_dst[cache_idx] + c * sstride) =
              *(uint32_t *)(cache_ptr[cache_idx] + c * bs);
        else
          *(uint64_t *)(cache_dst[cache_idx] + c * sstride) =
              *(uint64_t *)(cache_ptr[cache_idx] + c * bs);
    }
#else
    for (c = 0; c < bs; c++)
      if (bs == 4)
        *(uint32_t *)(cache_dst[cache_idx] + c * sstride) =
            *(uint32_t *)(cache_ptr[cache_idx] + c * bs);
      else
        *(uint64_t *)(cache_dst[cache_idx] + c * sstride) =
            *(uint64_t *)(cache_ptr[cache_idx] + c * bs);
#endif
  }

  aom_free(cache);
  aom_free(cache_ptr);
  aom_free(cache_dst);
}
