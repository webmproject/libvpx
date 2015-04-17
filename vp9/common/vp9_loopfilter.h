/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#ifndef VP9_COMMON_VP9_LOOPFILTER_H_
#define VP9_COMMON_VP9_LOOPFILTER_H_

#include "vpx_ports/mem.h"
#include "./vpx_config.h"

#include "vp9/common/vp9_blockd.h"
#include "vp9/common/vp9_seg_common.h"

#ifdef __cplusplus
extern "C" {
#endif

#define MAX_LOOP_FILTER 63
#define MAX_SHARPNESS 7

#define SIMD_WIDTH 16

#define MAX_REF_LF_DELTAS       4
#define MAX_MODE_LF_DELTAS      2

struct VP9Common;

#if CONFIG_LOOP_POSTFILTER
#define BILATERAL_LEVEL_BITS_KF 4
#define BILATERAL_LEVELS_KF     (1 << BILATERAL_LEVEL_BITS_KF)
#define BILATERAL_LEVEL_BITS    3
#define BILATERAL_LEVELS        (1 << BILATERAL_LEVEL_BITS)
#define DEF_BILATERAL_LEVEL     2

#define BILATERAL_PRECISION     8
#define BILATERAL_HALFWIN       3
#define BILATERAL_WIN           (2 * BILATERAL_HALFWIN + 1)

typedef struct bilateral_params {
  int sigma_x;  // spatial variance
  int sigma_r;  // range variance
} bilateral_params_t;

static bilateral_params_t
    bilateral_level_to_params_arr[BILATERAL_LEVELS + 1] = {
  // Values are rounded to 1/8 th precision
  {0, 0},    // 0 - default
  {4, 16},
  {5, 16},
  {6, 16},
  {7, 16},
  {9, 18},
  {12, 20},
  {16, 20},
  {20, 20},
};

static bilateral_params_t
    bilateral_level_to_params_arr_kf[BILATERAL_LEVELS_KF + 1] = {
  // Values are rounded to 1/8 th precision
  {0, 0},    // 0 - default
  {4, 16},
  {5, 16},
  {6, 16},
  {7, 16},
  {9, 18},
  {12, 20},
  {15, 22},
  {18, 24},
  {21, 24},
  {24, 24},
  {24, 28},
  {28, 24},
  {28, 28},
  {28, 32},
  {32, 24},
  {32, 28},
};

int vp9_bilateral_level_bits(const struct VP9Common *const cm);
int vp9_loop_bilateral_used(int level, int kf);

static INLINE bilateral_params_t vp9_bilateral_level_to_params(
    int index, int kf) {
  return kf ? bilateral_level_to_params_arr_kf[index] :
              bilateral_level_to_params_arr[index];
}
#endif  // CONFIG_LOOP_POSTFILTER

struct loopfilter {
  int filter_level;

  int sharpness_level;
  int last_sharpness_level;

  uint8_t mode_ref_delta_enabled;
  uint8_t mode_ref_delta_update;

  // 0 = Intra, Last, GF, ARF
  signed char ref_deltas[MAX_REF_LF_DELTAS];
  signed char last_ref_deltas[MAX_REF_LF_DELTAS];

  // 0 = ZERO_MV, MV
  signed char mode_deltas[MAX_MODE_LF_DELTAS];
  signed char last_mode_deltas[MAX_MODE_LF_DELTAS];

#if CONFIG_LOOP_POSTFILTER
  int bilateral_level;
#endif
};

// Need to align this structure so when it is declared and
// passed it can be loaded into vector registers.
typedef struct {
  DECLARE_ALIGNED(SIMD_WIDTH, uint8_t, mblim[SIMD_WIDTH]);
  DECLARE_ALIGNED(SIMD_WIDTH, uint8_t, lim[SIMD_WIDTH]);
  DECLARE_ALIGNED(SIMD_WIDTH, uint8_t, hev_thr[SIMD_WIDTH]);
} loop_filter_thresh;

typedef struct {
  loop_filter_thresh lfthr[MAX_LOOP_FILTER + 1];
  uint8_t lvl[MAX_SEGMENTS][MAX_REF_FRAMES][MAX_MODE_LF_DELTAS];
#if CONFIG_LOOP_POSTFILTER
  double wx_lut[BILATERAL_WIN * BILATERAL_WIN];
  double wr_lut[512];
  int bilateral_sigma_x_set;
  int bilateral_sigma_r_set;
  int bilateral_used;
#endif
} loop_filter_info_n;

// This structure holds bit masks for all 8x8 blocks in a 64x64 region.
// Each 1 bit represents a position in which we want to apply the loop filter.
// Left_ entries refer to whether we apply a filter on the border to the
// left of the block.   Above_ entries refer to whether or not to apply a
// filter on the above border.   Int_ entries refer to whether or not to
// apply borders on the 4x4 edges within the 8x8 block that each bit
// represents.
// Since each transform is accompanied by a potentially different type of
// loop filter there is a different entry in the array for each transform size.
typedef struct {
  uint64_t left_y[TX_SIZES];
  uint64_t above_y[TX_SIZES];
  uint64_t int_4x4_y;
  uint16_t left_uv[TX_SIZES];
  uint16_t above_uv[TX_SIZES];
  uint16_t int_4x4_uv;
  uint8_t lfl_y[64];
  uint8_t lfl_uv[16];
} LOOP_FILTER_MASK;

/* assorted loopfilter functions which get used elsewhere */
struct macroblockd;
struct VP9LfSyncData;

// This function sets up the bit masks for the entire 64x64 region represented
// by mi_row, mi_col.
void vp9_setup_mask(struct VP9Common *const cm,
                    const int mi_row, const int mi_col,
                    MODE_INFO *mi_8x8, const int mode_info_stride,
                    LOOP_FILTER_MASK *lfm);

void vp9_filter_block_plane(struct VP9Common *const cm,
                            struct macroblockd_plane *const plane,
                            int mi_row,
                            LOOP_FILTER_MASK *lfm);

void vp9_loop_filter_init(struct VP9Common *cm);

// Update the loop filter for the current frame.
// This should be called before vp9_loop_filter_rows(), vp9_loop_filter_frame()
// calls this function directly.
void vp9_loop_filter_frame_init(struct VP9Common *cm, int default_filt_lvl);

void vp9_loop_filter_frame(YV12_BUFFER_CONFIG *frame,
                           struct VP9Common *cm,
                           struct macroblockd *mbd,
                           int filter_level,
                           int y_only, int partial_frame);

// Apply the loop filter to [start, stop) macro block rows in frame_buffer.
void vp9_loop_filter_rows(YV12_BUFFER_CONFIG *frame_buffer,
                          struct VP9Common *cm,
                          struct macroblockd_plane planes[MAX_MB_PLANE],
                          int start, int stop, int y_only);
#if CONFIG_LOOP_POSTFILTER
void vp9_loop_bilateral_frame(YV12_BUFFER_CONFIG *frame,
                              struct VP9Common *cm,
                              int bilateral_level,
                              int y_only, int partial_frame);
void vp9_loop_filter_bilateral_frame(YV12_BUFFER_CONFIG *frame,
                                     struct VP9Common *cm,
                                     struct macroblockd *mbd,
                                     int frame_filter_level,
                                     int bilateral_level,
                                     int y_only, int partial_frame);
void vp9_loop_bilateral_init(loop_filter_info_n *lfi, int T, int kf);
void vp9_loop_bilateral_rows(YV12_BUFFER_CONFIG *frame,
                             struct VP9Common *cm,
                             int start_mi_row, int end_mi_row,
                             int y_only);
#endif  // CONFIG_LOOP_POSTFILTER


typedef struct LoopFilterWorkerData {
  YV12_BUFFER_CONFIG *frame_buffer;
  struct VP9Common *cm;
  struct macroblockd_plane planes[MAX_MB_PLANE];

  int start;
  int stop;
  int y_only;

  struct VP9LfSyncData *lf_sync;
  int num_lf_workers;
} LFWorkerData;

// Operates on the rows described by 'lf_data'.
int vp9_loop_filter_worker(LFWorkerData *const lf_data, void *unused);
#ifdef __cplusplus
}  // extern "C"
#endif

#endif  // VP9_COMMON_VP9_LOOPFILTER_H_
