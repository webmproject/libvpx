/*
 *  Copyright (c) 2016 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#ifndef AV1_COMMON_RESTORATION_H_
#define AV1_COMMON_RESTORATION_H_

#include "aom_ports/mem.h"
#include "./aom_config.h"

#include "av1/common/blockd.h"

#ifdef __cplusplus
extern "C" {
#endif

#define BILATERAL_LEVEL_BITS_KF 4
#define BILATERAL_LEVELS_KF (1 << BILATERAL_LEVEL_BITS_KF)
#define BILATERAL_LEVEL_BITS 3
#define BILATERAL_LEVELS (1 << BILATERAL_LEVEL_BITS)
// #define DEF_BILATERAL_LEVEL     2

#define RESTORATION_TILESIZE_SML 128
#define RESTORATION_TILESIZE_BIG 256
#define BILATERAL_SUBTILE_BITS 1
#define BILATERAL_SUBTILES (1 << (2 * BILATERAL_SUBTILE_BITS))

#define RESTORATION_HALFWIN 3
#define RESTORATION_HALFWIN1 (RESTORATION_HALFWIN + 1)
#define RESTORATION_WIN (2 * RESTORATION_HALFWIN + 1)
#define RESTORATION_WIN2 ((RESTORATION_WIN) * (RESTORATION_WIN))

#define RESTORATION_FILT_BITS 7
#define RESTORATION_FILT_STEP (1 << RESTORATION_FILT_BITS)

#define WIENER_FILT_TAP0_MINV (-5)
#define WIENER_FILT_TAP1_MINV (-23)
#define WIENER_FILT_TAP2_MINV (-16)

#define WIENER_FILT_TAP0_BITS 4
#define WIENER_FILT_TAP1_BITS 5
#define WIENER_FILT_TAP2_BITS 6

#define WIENER_FILT_BITS \
  ((WIENER_FILT_TAP0_BITS + WIENER_FILT_TAP1_BITS + WIENER_FILT_TAP2_BITS) * 2)

#define WIENER_FILT_TAP0_MAXV \
  (WIENER_FILT_TAP0_MINV - 1 + (1 << WIENER_FILT_TAP0_BITS))
#define WIENER_FILT_TAP1_MAXV \
  (WIENER_FILT_TAP1_MINV - 1 + (1 << WIENER_FILT_TAP1_BITS))
#define WIENER_FILT_TAP2_MAXV \
  (WIENER_FILT_TAP2_MINV - 1 + (1 << WIENER_FILT_TAP2_BITS))

typedef struct { int level[BILATERAL_SUBTILES]; } BilateralInfo;

typedef struct {
  int level;
  int vfilter[RESTORATION_WIN], hfilter[RESTORATION_WIN];
} WienerInfo;

typedef struct {
  RestorationType frame_restoration_type;
  RestorationType *restoration_type;
  // Bilateral filter
  BilateralInfo *bilateral_info;
  // Wiener filter
  WienerInfo *wiener_info;
} RestorationInfo;

typedef struct {
  RestorationInfo *rsi;
  int keyframe;
  int subsampling_x;
  int subsampling_y;
  int ntiles;
  int tile_width, tile_height;
  int nhtiles, nvtiles;
} RestorationInternal;

static INLINE int get_rest_tilesize(int width, int height) {
  if (width * height <= 352 * 288)
    return RESTORATION_TILESIZE_SML;
  else
    return RESTORATION_TILESIZE_BIG;
}

static INLINE int av1_get_rest_ntiles(int width, int height, int *tile_width,
                                      int *tile_height, int *nhtiles,
                                      int *nvtiles) {
  int nhtiles_, nvtiles_;
  int tile_width_, tile_height_;
  int tilesize = get_rest_tilesize(width, height);
  tile_width_ = (tilesize < 0) ? width : AOMMIN(tilesize, width);
  tile_height_ = (tilesize < 0) ? height : AOMMIN(tilesize, height);
  nhtiles_ = (width + (tile_width_ >> 1)) / tile_width_;
  nvtiles_ = (height + (tile_height_ >> 1)) / tile_height_;
  if (tile_width) *tile_width = tile_width_;
  if (tile_height) *tile_height = tile_height_;
  if (nhtiles) *nhtiles = nhtiles_;
  if (nvtiles) *nvtiles = nvtiles_;
  return (nhtiles_ * nvtiles_);
}

static INLINE void av1_get_rest_tile_limits(
    int tile_idx, int subtile_idx, int subtile_bits, int nhtiles, int nvtiles,
    int tile_width, int tile_height, int im_width, int im_height, int clamp_h,
    int clamp_v, int *h_start, int *h_end, int *v_start, int *v_end) {
  const int htile_idx = tile_idx % nhtiles;
  const int vtile_idx = tile_idx / nhtiles;
  *h_start = htile_idx * tile_width;
  *v_start = vtile_idx * tile_height;
  *h_end = (htile_idx < nhtiles - 1) ? *h_start + tile_width : im_width;
  *v_end = (vtile_idx < nvtiles - 1) ? *v_start + tile_height : im_height;
  if (subtile_bits) {
    const int num_subtiles_1d = (1 << subtile_bits);
    const int subtile_width = (*h_end - *h_start) >> subtile_bits;
    const int subtile_height = (*v_end - *v_start) >> subtile_bits;
    const int subtile_idx_h = subtile_idx & (num_subtiles_1d - 1);
    const int subtile_idx_v = subtile_idx >> subtile_bits;
    *h_start += subtile_idx_h * subtile_width;
    *v_start += subtile_idx_v * subtile_height;
    *h_end = subtile_idx_h == num_subtiles_1d - 1 ? *h_end
                                                  : *h_start + subtile_width;
    *v_end = subtile_idx_v == num_subtiles_1d - 1 ? *v_end
                                                  : *v_start + subtile_height;
  }
  if (clamp_h) {
    *h_start = AOMMAX(*h_start, RESTORATION_HALFWIN);
    *h_end = AOMMIN(*h_end, im_width - RESTORATION_HALFWIN);
  }
  if (clamp_v) {
    *v_start = AOMMAX(*v_start, RESTORATION_HALFWIN);
    *v_end = AOMMIN(*v_end, im_height - RESTORATION_HALFWIN);
  }
}

int av1_bilateral_level_bits(const struct AV1Common *const cm);
void av1_loop_restoration_init(RestorationInternal *rst, RestorationInfo *rsi,
                               int kf, int width, int height);
void av1_loop_restoration_frame(YV12_BUFFER_CONFIG *frame, struct AV1Common *cm,
                                RestorationInfo *rsi, int y_only,
                                int partial_frame);
void av1_loop_restoration_rows(YV12_BUFFER_CONFIG *frame, struct AV1Common *cm,
                               int start_mi_row, int end_mi_row, int y_only);
void av1_loop_restoration_precal();
#ifdef __cplusplus
}  // extern "C"
#endif

#endif  // AV1_COMMON_RESTORATION_H_
