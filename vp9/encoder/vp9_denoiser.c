/*
 *  Copyright (c) 2012 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include <assert.h>
#include "vpx_scale/yv12config.h"
#include "vpx/vpx_integer.h"
#include "vp9/encoder/vp9_denoiser.h"

static const int widths[]  = {4, 4, 8, 8,  8, 16, 16, 16, 32, 32, 32, 64, 64};
static const int heights[] = {4, 8, 4, 8, 16,  8, 16, 32, 16, 32, 64, 32, 64};

int vp9_denoiser_filter() {
  return 0;
}

static int update_running_avg(const uint8_t *mc_avg, int mc_avg_stride,
                              uint8_t *avg, int avg_stride,
                              const uint8_t *sig, int sig_stride,
                              int increase_denoising, BLOCK_SIZE bs) {
  int r, c;
  int diff, adj, absdiff;
  int shift_inc1 = 0, shift_inc2 = 1;
  int adj_val[] = {3, 4, 6};
  int total_adj = 0;

  if (increase_denoising) {
    shift_inc1 = 1;
    shift_inc2 = 2;
  }

  for (r = 0; r < heights[bs]; ++r) {
    for (c = 0; c < widths[bs]; ++c) {
      diff = mc_avg[c] - sig[c];
      absdiff = abs(diff);

      if (absdiff <= 3 + shift_inc1) {
        avg[c] = mc_avg[c];
        total_adj += diff;
      } else {
        switch (absdiff) {
          case 4: case 5: case 6: case 7:
            adj = adj_val[0];
            break;
          case 8: case 9: case 10: case 11:
          case 12: case 13: case 14: case 15:
            adj = adj_val[1];
            break;
          default:
            adj = adj_val[2];
        }
        if (diff > 0) {
          avg[c] = MIN(UINT8_MAX, sig[c] + adj);
          total_adj += adj;
        } else {
          avg[c] = MAX(0, sig[c] - adj);
          total_adj -= adj;
        }
      }
    }
    sig += sig_stride;
    avg += avg_stride;
    mc_avg += mc_avg_stride;
  }
  return total_adj;
}

static uint8_t *block_start(uint8_t *framebuf, int stride,
                            int mi_row, int mi_col) {
  return framebuf + (stride * mi_row * 8) + (mi_col * 8);
}

void copy_block(uint8_t *dest, int dest_stride,
                uint8_t *src, int src_stride, BLOCK_SIZE bs) {
  int r, c;
  for (r = 0; r < heights[bs]; ++r) {
    for (c = 0; c < widths[bs]; ++c) {
      dest[c] = src[c];
    }
    dest += dest_stride;
    src += src_stride;
  }
}

void vp9_denoiser_denoise(VP9_DENOISER *denoiser, MACROBLOCK *mb,
                          int mi_row, int mi_col, BLOCK_SIZE bs) {
  int decision = COPY_BLOCK;

  YV12_BUFFER_CONFIG avg = denoiser->running_avg_y[INTRA_FRAME];
  struct buf_2d src = mb->plane[0].src;

  update_running_avg(denoiser->mc_running_avg_y.y_buffer,
                     denoiser->mc_running_avg_y.y_stride,
                     denoiser->running_avg_y[INTRA_FRAME].y_buffer,
                     denoiser->running_avg_y[INTRA_FRAME].y_stride,
                     mb->plane[0].src.buf, mb->plane[0].src.stride, 0, bs);

  if (decision == FILTER_BLOCK) {
    // TODO(tkopp)
  }
  if (decision == COPY_BLOCK) {
    copy_block(block_start(avg.y_buffer, avg.y_stride, mi_row, mi_col),
               avg.y_stride, src.buf, src.stride, bs);
  }
}

static void copy_frame(YV12_BUFFER_CONFIG dest, const YV12_BUFFER_CONFIG src) {
  int r, c;
  const uint8_t *srcbuf = src.y_buffer;
  uint8_t *destbuf = dest.y_buffer;
  assert(dest.y_width == src.y_width);
  assert(dest.y_height == src.y_height);

  for (r = 0; r < dest.y_height; ++r) {
    for (c = 0; c < dest.y_width; ++c) {
      destbuf[c] = srcbuf[c];
    }
    destbuf += dest.y_stride;
    srcbuf += src.y_stride;
  }
}

void vp9_denoiser_update_frame_info(VP9_DENOISER *denoiser,
                                    YV12_BUFFER_CONFIG src,
                                    FRAME_TYPE frame_type,
                                    int refresh_alt_ref_frame,
                                    int refresh_golden_frame,
                                    int refresh_last_frame) {
  if (frame_type == KEY_FRAME) {
    int i;
    copy_frame(denoiser->running_avg_y[LAST_FRAME], src);
    for (i = 2; i < MAX_REF_FRAMES - 1; i++) {
      copy_frame(denoiser->running_avg_y[i],
                 denoiser->running_avg_y[LAST_FRAME]);
    }
  } else {  /* For non key frames */
    if (refresh_alt_ref_frame) {
      copy_frame(denoiser->running_avg_y[ALTREF_FRAME],
                 denoiser->running_avg_y[INTRA_FRAME]);
    }
    if (refresh_golden_frame) {
      copy_frame(denoiser->running_avg_y[GOLDEN_FRAME],
                 denoiser->running_avg_y[INTRA_FRAME]);
    }
    if (refresh_last_frame) {
      copy_frame(denoiser->running_avg_y[LAST_FRAME],
                 denoiser->running_avg_y[INTRA_FRAME]);
    }
  }
}

void vp9_denoiser_update_frame_stats() {
}

int vp9_denoiser_alloc(VP9_DENOISER *denoiser, int width, int height,
                       int ssx, int ssy, int border) {
  int i, fail;
  assert(denoiser != NULL);

  for (i = 0; i < MAX_REF_FRAMES; ++i) {
    fail = vp9_alloc_frame_buffer(&denoiser->running_avg_y[i], width, height,
                                  ssx, ssy, border);
    if (fail) {
      vp9_denoiser_free(denoiser);
      return 1;
    }
  }

  fail = vp9_alloc_frame_buffer(&denoiser->mc_running_avg_y, width, height,
                                ssx, ssy, border);
  if (fail) {
    vp9_denoiser_free(denoiser);
    return 1;
  }

  return 0;
}

void vp9_denoiser_free(VP9_DENOISER *denoiser) {
  int i;
  if (denoiser == NULL) {
    return;
  }
  for (i = 0; i < MAX_REF_FRAMES; ++i) {
    if (&denoiser->running_avg_y[i] != NULL) {
      vp9_free_frame_buffer(&denoiser->running_avg_y[i]);
    }
  }
  if (&denoiser->mc_running_avg_y != NULL) {
    vp9_free_frame_buffer(&denoiser->mc_running_avg_y);
  }
}
