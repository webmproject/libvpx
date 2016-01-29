/*
 *  Copyright (c) 2016 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#ifndef VP10_COMMON_RESTORATION_H_
#define VP10_COMMON_RESTORATION_H_

#include "vpx_ports/mem.h"
#include "./vpx_config.h"

#include "vp10/common/blockd.h"

#ifdef __cplusplus
extern "C" {
#endif

#define RESTORATION_LEVEL_BITS_KF 4
#define RESTORATION_LEVELS_KF     (1 << RESTORATION_LEVEL_BITS_KF)
#define RESTORATION_LEVEL_BITS    3
#define RESTORATION_LEVELS        (1 << RESTORATION_LEVEL_BITS)
#define DEF_RESTORATION_LEVEL     2

#define RESTORATION_PRECISION     16
#define RESTORATION_HALFWIN       3
#define RESTORATION_WIN           (2 * RESTORATION_HALFWIN + 1)

typedef struct restoration_params {
  int sigma_x;  // spatial variance x
  int sigma_y;  // spatial variance y
  int sigma_r;  // range variance
} restoration_params_t;

static restoration_params_t
    restoration_level_to_params_arr[RESTORATION_LEVELS + 1] = {
  // Values are rounded to 1/16 th precision
  {0, 0, 0},    // 0 - default
  {8, 9, 30},
  {9, 8, 30},
  {9, 11, 32},
  {11, 9, 32},
  {14, 14, 32},
  {18, 18, 36},
  {24, 24, 40},
  {32, 32, 40},
};

static restoration_params_t
    restoration_level_to_params_arr_kf[RESTORATION_LEVELS_KF + 1] = {
  // Values are rounded to 1/16 th precision
  {0, 0, 0},    // 0 - default
  {8, 8, 30},
  {9, 9, 32},
  {10, 10, 32},
  {12, 12, 32},
  {14, 14, 32},
  {18, 18, 36},
  {24, 24, 40},
  {30, 30, 44},
  {36, 36, 48},
  {42, 42, 48},
  {48, 48, 48},
  {48, 48, 56},
  {56, 56, 48},
  {56, 56, 56},
  {56, 56, 64},
  {64, 64, 48},
};

typedef struct {
  double *wx_lut[RESTORATION_WIN];
  double *wr_lut;
  int restoration_sigma_x_set;
  int restoration_sigma_y_set;
  int restoration_sigma_r_set;
  int restoration_used;
} restoration_info_n;

int vp10_restoration_level_bits(const struct VP10Common *const cm);
int vp10_loop_restoration_used(int level, int kf);

static INLINE restoration_params_t vp10_restoration_level_to_params(
    int index, int kf) {
  return kf ? restoration_level_to_params_arr_kf[index] :
              restoration_level_to_params_arr[index];
}

void vp10_loop_restoration_init(restoration_info_n *rst, int T, int kf);
void vp10_loop_restoration_frame(YV12_BUFFER_CONFIG *frame,
                                 struct VP10Common *cm,
                                 int restoration_level,
                                 int y_only, int partial_frame);
void vp10_loop_restoration_rows(YV12_BUFFER_CONFIG *frame,
                                struct VP10Common *cm,
                                int start_mi_row, int end_mi_row,
                                int y_only);
void vp10_loop_restoration_precal();
#ifdef __cplusplus
}  // extern "C"
#endif

#endif  // VP10_COMMON_RESTORATION_H_
