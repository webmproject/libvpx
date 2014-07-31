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

struct VP9Common;
struct VP9Decoder;

typedef struct TileWorkerData {
  struct VP9Decoder *pbi;
  vp9_reader bit_reader;
  DECLARE_ALIGNED(16, struct macroblockd, xd);

  // Row-based parallel loopfilter data
  LFWorkerData lfdata;
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
} VP9LfSync;

// WorkerData for the FrameWorker thread. It contains all the information of
// the worker and decode structures for decoding a frame.
typedef struct FrameWorkerData {
  struct VP9Decoder *pbi;
  const uint8_t *data;
  const uint8_t *data_end;
  size_t data_size;
  void *user_priv;
  int result;
  int worker_id;

  // scratch_buffer is used in frame parallel mode only.
  // It is used to make a copy of the compressed data.
  uint8_t *scratch_buffer;
  size_t scratch_buffer_size;

#if CONFIG_MULTITHREAD
  pthread_mutex_t stats_mutex;
  pthread_cond_t stats_cond;
#endif

  int frame_context_ready;  // Current frame's context is ready to read.
  int frame_decoded;        // Finished decoding current frame.
} FrameWorkerData;

// Allocate memory for loopfilter row synchronization.
void vp9_loop_filter_alloc(struct VP9Common *cm, VP9LfSync *lf_sync,
                           int rows, int width);

// Deallocate loopfilter synchronization related mutex and data.
void vp9_loop_filter_dealloc(VP9LfSync *lf_sync, int rows);

// Multi-threaded loopfilter that uses the tile threads.
void vp9_loop_filter_frame_mt(YV12_BUFFER_CONFIG *frame,
                              struct VP9Decoder *pbi,
                              struct VP9Common *cm,
                              int frame_filter_level,
                              int y_only);

void vp9_frameworker_lock_stats(VP9Worker *const worker);
void vp9_frameworker_unlock_stats(VP9Worker *const worker);
void vp9_frameworker_signal_stats(VP9Worker *const worker);

// Wait until ref_buf has been decoded to row in real pixel unit.
// Note: worker may already finish decoding ref_buf and release it in order to
// start decoding next frame. So need to check whether worker is still decoding
// ref_buf.
void vp9_frameworker_wait(VP9Worker *const worker, RefCntBuffer *const ref_buf,
                          int row);

// FrameWorker broadcasts its decoding progress so other workers that are
// waiting on it can resume decoding.
void vp9_frameworker_broadcast(RefCntBuffer *const buf, int row);

// Copy necessary decoding context from src worker to dst worker.
void vp9_frameworker_copy_context(VP9Worker *const dst_worker,
                                  VP9Worker *const src_worker);

#endif  // VP9_DECODER_VP9_DTHREAD_H_
