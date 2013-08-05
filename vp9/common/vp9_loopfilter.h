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
#include "vpx_config.h"

#include "vp9/common/vp9_blockd.h"
#include "vp9/common/vp9_seg_common.h"

#define MAX_LOOP_FILTER 63
#define MAX_SHARPNESS 7

#define SIMD_WIDTH 16

// Need to align this structure so when it is declared and
// passed it can be loaded into vector registers.
typedef struct {
  DECLARE_ALIGNED(SIMD_WIDTH, uint8_t,
                  mblim[MAX_LOOP_FILTER + 1][SIMD_WIDTH]);
  DECLARE_ALIGNED(SIMD_WIDTH, uint8_t,
                  lim[MAX_LOOP_FILTER + 1][SIMD_WIDTH]);
  DECLARE_ALIGNED(SIMD_WIDTH, uint8_t,
                  hev_thr[4][SIMD_WIDTH]);
  uint8_t lvl[MAX_SEGMENTS][MAX_REF_FRAMES][MAX_MODE_LF_DELTAS];
  uint8_t mode_lf_lut[MB_MODE_COUNT];
} loop_filter_info_n;

/* assorted loopfilter functions which get used elsewhere */
struct VP9Common;
struct macroblockd;

void vp9_loop_filter_init(struct VP9Common *cm, struct loopfilter *lf);

// Update the loop filter for the current frame.
// This should be called before vp9_loop_filter_rows(), vp9_loop_filter_frame()
// calls this function directly.
void vp9_loop_filter_frame_init(struct VP9Common *const cm,
                                struct macroblockd *const xd,
                                int default_filt_lvl);

void vp9_loop_filter_frame(struct VP9Common *cm,
                           struct macroblockd *mbd,
                           int filter_level,
                           int y_only);

// Apply the loop filter to [start, stop) macro block rows in frame_buffer.
void vp9_loop_filter_rows(const YV12_BUFFER_CONFIG *frame_buffer,
                          struct VP9Common *cm, struct macroblockd *xd,
                          int start, int stop, int y_only);

typedef struct LoopFilterWorkerData {
  const YV12_BUFFER_CONFIG *frame_buffer;
  struct VP9Common *cm;
  struct macroblockd xd;  // TODO(jzern): most of this is unnecessary to the
                          // loopfilter. the planes are necessary as their state
                          // is changed during decode.
  int start;
  int stop;
  int y_only;
} LFWorkerData;

// Operates on the rows described by LFWorkerData passed as 'arg1'.
int vp9_loop_filter_worker(void *arg1, void *arg2);
#endif  // VP9_COMMON_VP9_LOOPFILTER_H_
