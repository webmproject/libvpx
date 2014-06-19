/*
 *  Copyright (c) 2012 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include <stdio.h>
#include <stdint.h>
#include "vp9/encoder/vp9_denoiser.h"
#include "vpx_scale/yv12config.h"

static const int widths[]  = {4, 4, 8, 8,  8, 16, 16, 16, 32, 32, 32, 64, 64};
static const int heights[] = {4, 8, 4, 8, 16,  8, 16, 32, 16, 32, 64, 32, 64};

int vp9_denoiser_filter() {
  return 0;
}

int update_running_avg(uint8_t *mc_avg, int mc_avg_stride, uint8_t *avg,
                       int avg_stride, uint8_t *sig, int sig_stride,
                       BLOCK_SIZE bs) {
  //                           Indices: 0, 1, 2, 3, 4, 5 ,6, 7,
  static const uint8_t adjustments[] = {0, 0, 0, 0, 3, 3, 3, 3,
  //                                    8, 9,10,11,12,13,14,15,16
                                        4, 4, 4, 4, 4, 4, 4, 4, 6};
  int r, c;
  int diff;
  int adjustment;
  int total_adj = 0;

  for (r = 0; r < heights[bs]; ++r) {
    for (c = 0; c < widths[bs]; ++c) {
      diff = mc_avg[c] - sig[c];
      adjustment = adjustments[MIN(abs(diff), 16)];

      if (diff > 0) {
        avg[c] = MIN(UINT8_MAX, sig[c] + adjustment);
        total_adj += adjustment;
      } else {
        avg[c] = MAX(0, sig[c] - adjustment);
        total_adj -= adjustment;
      }
    }
    sig += sig_stride;
    avg += avg_stride;
    mc_avg += mc_avg_stride;
  }
  return total_adj;
}

void vp9_denoiser_denoise(VP9_DENOISER *denoiser,
                          MACROBLOCK *mb, MODE_INFO **grid,
                          int mi_row, int mi_col, BLOCK_SIZE bs) {
  update_running_avg(denoiser->mc_running_avg_y.buf,
                     denoiser->mc_running_avg_y.stride,
                     denoiser->running_avg_y.buf,
                     denoiser->running_avg_y.stride,
                     mb->plane[0].src.buf, mb->plane[0].src.stride, bs);
  return;
}

void vp9_denoiser_update_frame_info(VP9_DENOISER *denoiser,
                                    FRAME_TYPE frame_type,
                                    int refresh_alt_ref_frame,
                                    int refresh_golden_frame,
                                    int refresh_last_frame) {
  return;
}

void vp9_denoiser_update_frame_stats() {
  return;
}

int vp9_denoiser_alloc(VP9_DENOISER *denoiser, int width, int height,
                       int border) {
  assert(denoiser);

  denoiser->running_avg_y.stride = width + 2 * border;

  denoiser->running_avg_y.buf = calloc(
      ((2 * border) + width) * ((2 * border) + height), sizeof(uint8_t));
  if (denoiser->running_avg_y.buf == NULL) {
    vp9_denoiser_free(denoiser);
    return 1;
  }

  denoiser->mc_running_avg_y.stride = width + 2 * border;

  denoiser->mc_running_avg_y.buf = calloc(
      ((2 * border) + width) * ((2 * border) + height), sizeof(uint8_t));
  if (denoiser->mc_running_avg_y.buf == NULL) {
    vp9_denoiser_free(denoiser);
    return 1;
  }

  return 0;
}

void vp9_denoiser_free(VP9_DENOISER *denoiser) {
  if (denoiser->running_avg_y.buf != NULL) {
    free(denoiser->running_avg_y.buf);
  }
  if (denoiser->mc_running_avg_y.buf != NULL) {
    free(denoiser->mc_running_avg_y.buf);
  }
  return;
}

