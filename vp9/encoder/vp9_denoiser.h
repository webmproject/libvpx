/*
 *  Copyright (c) 2012 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#ifndef VP9_ENCODER_DENOISER_H_
#define VP9_ENCODER_DENOISER_H_

#include "vp9/encoder/vp9_block.h"

#ifdef __cplusplus
extern "C" {
#endif

enum vp9_denoiser_decision {
  COPY_BLOCK,
  FILTER_BLOCK
};

typedef struct vp9_denoiser {
  struct buf_2d running_avg_y;
  struct buf_2d mc_running_avg_y;
} VP9_DENOISER;

void vp9_denoiser_update_frame_info(VP9_DENOISER *denoiser,
                                    FRAME_TYPE frame_type,
                                    int refresh_alt_ref_frame,
                                    int refresh_golden_frame,
                                    int refresh_last_frame);

void vp9_denoiser_denoise(VP9_DENOISER *denoiser,
                          MACROBLOCK *mb, MODE_INFO **grid,
                          int mi_row, int mi_col, BLOCK_SIZE bs);

void vp9_denoiser_update_frame_stats();

int vp9_denoiser_alloc(VP9_DENOISER *denoiser, int width, int height,
                       int border);

void vp9_denoiser_free(VP9_DENOISER *denoiser);

#ifdef __cplusplus
}  // extern "C"
#endif

#endif  // VP9_ENCODER_DENOISER_H_
