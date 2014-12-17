/*
 *  Copyright (c) 2014 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#ifndef VP9_DECODER_VP9_DTHREAD_H_
#define VP9_DECODER_VP9_DTHREAD_H_

#include "./vpx_config.h"
#include "vp9/common/vp9_thread.h"
#include "vp9/decoder/vp9_reader.h"
#include "vpx/internal/vpx_codec_internal.h"

struct VP9Common;
struct VP9Decoder;

typedef struct TileWorkerData {
  struct VP9Common *cm;
  vp9_reader bit_reader;
  DECLARE_ALIGNED(16, struct macroblockd, xd);
  struct vpx_internal_error_info error_info;
} TileWorkerData;

// Loopfilter row synchronization
typedef struct VP9LfSyncData {
#if CONFIG_MULTITHREAD
  pthread_mutex_t *mutex_;
  pthread_cond_t *cond_;
#endif
  // Allocate memory to store the loop-filtered superblock index in each row.
  int *cur_sb_col;
  // The optimal sync_range for different resolution and platform should be
  // determined by testing. Currently, it is chosen to be a power-of-2 number.
  int sync_range;
  int rows;

  // Row-based parallel loopfilter data
  LFWorkerData *lfdata;
  int num_workers;
} VP9LfSync;

// Allocate memory for loopfilter row synchronization.
void vp9_loop_filter_alloc(VP9LfSync *lf_sync, VP9_COMMON *cm, int rows,
                           int width, int num_workers);

// Deallocate loopfilter synchronization related mutex and data.
void vp9_loop_filter_dealloc(VP9LfSync *lf_sync);

// Multi-threaded loopfilter that uses the tile threads.
void vp9_loop_filter_frame_mt(VP9LfSync *lf_sync,
                              YV12_BUFFER_CONFIG *frame,
                              struct macroblockd_plane planes[MAX_MB_PLANE],
                              struct VP9Common *cm,
                              VP9Worker *workers, int num_workers,
                              int frame_filter_level,
                              int y_only);

#endif  // VP9_DECODER_VP9_DTHREAD_H_
