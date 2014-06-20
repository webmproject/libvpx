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

void vp9_denoiser_denoise(VP9_DENOISER *denoiser,
                          MACROBLOCK *mb, MODE_INFO **grid,
                          int mi_row, int mi_col, BLOCK_SIZE bs) {
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
  return 0;
}

void vp9_denoiser_free(VP9_DENOISER *denoiser) {
  return;
}

