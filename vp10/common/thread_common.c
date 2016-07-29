/*
 *  Copyright (c) 2014 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include "./vpx_config.h"
#include "vpx_dsp/vpx_dsp_common.h"
#include "vpx_mem/vpx_mem.h"
#include "vp10/common/entropymode.h"
#include "vp10/common/thread_common.h"
#include "vp10/common/reconinter.h"
#include "vp10/common/loopfilter.h"

#if CONFIG_MULTITHREAD
static INLINE void mutex_lock(pthread_mutex_t *const mutex) {
  const int kMaxTryLocks = 4000;
  int locked = 0;
  int i;

  for (i = 0; i < kMaxTryLocks; ++i) {
    if (!pthread_mutex_trylock(mutex)) {
      locked = 1;
      break;
    }
  }

  if (!locked)
    pthread_mutex_lock(mutex);
}
#endif  // CONFIG_MULTITHREAD

static INLINE void sync_read(VP10LfSync *const lf_sync, int r, int c) {
#if CONFIG_MULTITHREAD
  const int nsync = lf_sync->sync_range;

  if (r && !(c & (nsync - 1))) {
    pthread_mutex_t *const mutex = &lf_sync->mutex_[r - 1];
    mutex_lock(mutex);

    while (c > lf_sync->cur_sb_col[r - 1] - nsync) {
      pthread_cond_wait(&lf_sync->cond_[r - 1], mutex);
    }
    pthread_mutex_unlock(mutex);
  }
#else
  (void)lf_sync;
  (void)r;
  (void)c;
#endif  // CONFIG_MULTITHREAD
}

static INLINE void sync_write(VP10LfSync *const lf_sync, int r, int c,
                              const int sb_cols) {
#if CONFIG_MULTITHREAD
  const int nsync = lf_sync->sync_range;
  int cur;
  // Only signal when there are enough filtered SB for next row to run.
  int sig = 1;

  if (c < sb_cols - 1) {
    cur = c;
    if (c % nsync)
      sig = 0;
  } else {
    cur = sb_cols + nsync;
  }

  if (sig) {
    mutex_lock(&lf_sync->mutex_[r]);

    lf_sync->cur_sb_col[r] = cur;

    pthread_cond_signal(&lf_sync->cond_[r]);
    pthread_mutex_unlock(&lf_sync->mutex_[r]);
  }
#else
  (void)lf_sync;
  (void)r;
  (void)c;
  (void)sb_cols;
#endif  // CONFIG_MULTITHREAD
}

// Implement row loopfiltering for each thread.
static INLINE
void thread_loop_filter_rows(const YV12_BUFFER_CONFIG *const frame_buffer,
                             VP10_COMMON *const cm,
                             struct macroblockd_plane planes[MAX_MB_PLANE],
                             int start, int stop, int y_only,
                             VP10LfSync *const lf_sync) {
  const int num_planes = y_only ? 1 : MAX_MB_PLANE;
  const int sb_cols = mi_cols_aligned_to_sb(cm) >> cm->mib_size_log2;
  int mi_row, mi_col;
#if !CONFIG_EXT_PARTITION_TYPES
  enum lf_path path;
  LOOP_FILTER_MASK lfm;
  if (y_only)
    path = LF_PATH_444;
  else if (planes[1].subsampling_y == 1 && planes[1].subsampling_x == 1)
    path = LF_PATH_420;
  else if (planes[1].subsampling_y == 0 && planes[1].subsampling_x == 0)
    path = LF_PATH_444;
  else
    path = LF_PATH_SLOW;
#endif  // !CONFIG_EXT_PARTITION_TYPES

#if CONFIG_EXT_PARTITION
  printf("STOPPING: This code has not been modified to work with the "
         "extended coding unit size experiment");
  exit(EXIT_FAILURE);
#endif  // CONFIG_EXT_PARTITION

  for (mi_row = start; mi_row < stop;
       mi_row += lf_sync->num_workers * cm->mib_size) {
    MODE_INFO **const mi = cm->mi_grid_visible + mi_row * cm->mi_stride;

    for (mi_col = 0; mi_col < cm->mi_cols; mi_col += cm->mib_size) {
      const int r = mi_row >> cm->mib_size_log2;
      const int c = mi_col >> cm->mib_size_log2;
      int plane;

      sync_read(lf_sync, r, c);

      vp10_setup_dst_planes(planes, frame_buffer, mi_row, mi_col);

#if CONFIG_EXT_PARTITION_TYPES
      for (plane = 0; plane < num_planes; ++plane)
        vp10_filter_block_plane_non420(cm, &planes[plane], mi + mi_col,
                                       mi_row, mi_col);
#else
      // TODO(JBB): Make setup_mask work for non 420.
      vp10_setup_mask(cm, mi_row, mi_col, mi + mi_col, cm->mi_stride,
                     &lfm);

      vp10_filter_block_plane_ss00(cm, &planes[0], mi_row, &lfm);
      for (plane = 1; plane < num_planes; ++plane) {
        switch (path) {
          case LF_PATH_420:
            vp10_filter_block_plane_ss11(cm, &planes[plane], mi_row, &lfm);
            break;
          case LF_PATH_444:
            vp10_filter_block_plane_ss00(cm, &planes[plane], mi_row, &lfm);
            break;
          case LF_PATH_SLOW:
            vp10_filter_block_plane_non420(cm, &planes[plane], mi + mi_col,
                                          mi_row, mi_col);
            break;
        }
      }
#endif  // CONFIG_EXT_PARTITION_TYPES
      sync_write(lf_sync, r, c, sb_cols);
    }
  }
}

// Row-based multi-threaded loopfilter hook
static int loop_filter_row_worker(VP10LfSync *const lf_sync,
                                  LFWorkerData *const lf_data) {
  thread_loop_filter_rows(lf_data->frame_buffer, lf_data->cm, lf_data->planes,
                          lf_data->start, lf_data->stop, lf_data->y_only,
                          lf_sync);
  return 1;
}

static void loop_filter_rows_mt(YV12_BUFFER_CONFIG *frame,
                                VP10_COMMON *cm,
                                struct macroblockd_plane planes[MAX_MB_PLANE],
                                int start, int stop, int y_only,
                                VPxWorker *workers, int nworkers,
                                VP10LfSync *lf_sync) {
  const VPxWorkerInterface *const winterface = vpx_get_worker_interface();
  // Number of superblock rows and cols
  const int sb_rows = mi_rows_aligned_to_sb(cm) >> cm->mib_size_log2;
  // Decoder may allocate more threads than number of tiles based on user's
  // input.
  const int tile_cols = cm->tile_cols;
  const int num_workers = VPXMIN(nworkers, tile_cols);
  int i;

#if CONFIG_EXT_PARTITION
      printf("STOPPING: This code has not been modified to work with the "
             "extended coding unit size experiment");
      exit(EXIT_FAILURE);
#endif  // CONFIG_EXT_PARTITION

  if (!lf_sync->sync_range || sb_rows != lf_sync->rows ||
      num_workers > lf_sync->num_workers) {
    vp10_loop_filter_dealloc(lf_sync);
    vp10_loop_filter_alloc(lf_sync, cm, sb_rows, cm->width, num_workers);
  }

  // Initialize cur_sb_col to -1 for all SB rows.
  memset(lf_sync->cur_sb_col, -1, sizeof(*lf_sync->cur_sb_col) * sb_rows);

  // Set up loopfilter thread data.
  // The decoder is capping num_workers because it has been observed that using
  // more threads on the loopfilter than there are cores will hurt performance
  // on Android. This is because the system will only schedule the tile decode
  // workers on cores equal to the number of tile columns. Then if the decoder
  // tries to use more threads for the loopfilter, it will hurt performance
  // because of contention. If the multithreading code changes in the future
  // then the number of workers used by the loopfilter should be revisited.
  for (i = 0; i < num_workers; ++i) {
    VPxWorker *const worker = &workers[i];
    LFWorkerData *const lf_data = &lf_sync->lfdata[i];

    worker->hook = (VPxWorkerHook)loop_filter_row_worker;
    worker->data1 = lf_sync;
    worker->data2 = lf_data;

    // Loopfilter data
    vp10_loop_filter_data_reset(lf_data, frame, cm, planes);
    lf_data->start = start + i * cm->mib_size;
    lf_data->stop = stop;
    lf_data->y_only = y_only;

    // Start loopfiltering
    if (i == num_workers - 1) {
      winterface->execute(worker);
    } else {
      winterface->launch(worker);
    }
  }

  // Wait till all rows are finished
  for (i = 0; i < num_workers; ++i) {
    winterface->sync(&workers[i]);
  }
}

void vp10_loop_filter_frame_mt(YV12_BUFFER_CONFIG *frame,
                              VP10_COMMON *cm,
                              struct macroblockd_plane planes[MAX_MB_PLANE],
                              int frame_filter_level,
                              int y_only, int partial_frame,
                              VPxWorker *workers, int num_workers,
                              VP10LfSync *lf_sync) {
  int start_mi_row, end_mi_row, mi_rows_to_filter;

  if (!frame_filter_level) return;

  start_mi_row = 0;
  mi_rows_to_filter = cm->mi_rows;
  if (partial_frame && cm->mi_rows > 8) {
    start_mi_row = cm->mi_rows >> 1;
    start_mi_row &= 0xfffffff8;
    mi_rows_to_filter = VPXMAX(cm->mi_rows / 8, 8);
  }
  end_mi_row = start_mi_row + mi_rows_to_filter;
  vp10_loop_filter_frame_init(cm, frame_filter_level);

  loop_filter_rows_mt(frame, cm, planes, start_mi_row, end_mi_row,
                      y_only, workers, num_workers, lf_sync);
}

// Set up nsync by width.
static INLINE int get_sync_range(int width) {
  // nsync numbers are picked by testing. For example, for 4k
  // video, using 4 gives best performance.
  if (width < 640)
    return 1;
  else if (width <= 1280)
    return 2;
  else if (width <= 4096)
    return 4;
  else
    return 8;
}

// Allocate memory for lf row synchronization
void vp10_loop_filter_alloc(VP10LfSync *lf_sync, VP10_COMMON *cm, int rows,
                           int width, int num_workers) {
  lf_sync->rows = rows;
#if CONFIG_MULTITHREAD
  {
    int i;

    CHECK_MEM_ERROR(cm, lf_sync->mutex_,
                    vpx_malloc(sizeof(*lf_sync->mutex_) * rows));
    if (lf_sync->mutex_) {
      for (i = 0; i < rows; ++i) {
        pthread_mutex_init(&lf_sync->mutex_[i], NULL);
      }
    }

    CHECK_MEM_ERROR(cm, lf_sync->cond_,
                    vpx_malloc(sizeof(*lf_sync->cond_) * rows));
    if (lf_sync->cond_) {
      for (i = 0; i < rows; ++i) {
        pthread_cond_init(&lf_sync->cond_[i], NULL);
      }
    }
  }
#endif  // CONFIG_MULTITHREAD

  CHECK_MEM_ERROR(cm, lf_sync->lfdata,
                  vpx_malloc(num_workers * sizeof(*lf_sync->lfdata)));
  lf_sync->num_workers = num_workers;

  CHECK_MEM_ERROR(cm, lf_sync->cur_sb_col,
                  vpx_malloc(sizeof(*lf_sync->cur_sb_col) * rows));

  // Set up nsync.
  lf_sync->sync_range = get_sync_range(width);
}

// Deallocate lf synchronization related mutex and data
void vp10_loop_filter_dealloc(VP10LfSync *lf_sync) {
  if (lf_sync != NULL) {
#if CONFIG_MULTITHREAD
    int i;

    if (lf_sync->mutex_ != NULL) {
      for (i = 0; i < lf_sync->rows; ++i) {
        pthread_mutex_destroy(&lf_sync->mutex_[i]);
      }
      vpx_free(lf_sync->mutex_);
    }
    if (lf_sync->cond_ != NULL) {
      for (i = 0; i < lf_sync->rows; ++i) {
        pthread_cond_destroy(&lf_sync->cond_[i]);
      }
      vpx_free(lf_sync->cond_);
    }
#endif  // CONFIG_MULTITHREAD
    vpx_free(lf_sync->lfdata);
    vpx_free(lf_sync->cur_sb_col);
    // clear the structure as the source of this call may be a resize in which
    // case this call will be followed by an _alloc() which may fail.
    vp10_zero(*lf_sync);
  }
}

// Accumulate frame counts. FRAME_COUNTS consist solely of 'unsigned int'
// members, so we treat it as an array, and sum over the whole length.
void vp10_accumulate_frame_counts(VP10_COMMON *cm, FRAME_COUNTS *counts) {
  unsigned int *const acc = (unsigned int*)&cm->counts;
  const unsigned int *const cnt = (unsigned int*)counts;

  const unsigned int n_counts = sizeof(FRAME_COUNTS)/sizeof(unsigned int);
  unsigned int i;

  for (i = 0; i < n_counts; i++)
    acc[i] += cnt[i];
}
