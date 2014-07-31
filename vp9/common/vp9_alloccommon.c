/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include "./vpx_config.h"
#include "vpx_mem/vpx_mem.h"

#include "vp9/common/vp9_blockd.h"
#include "vp9/common/vp9_common.h"
#include "vp9/common/vp9_entropymode.h"
#include "vp9/common/vp9_entropymv.h"
#include "vp9/common/vp9_onyxc_int.h"
#include "vp9/common/vp9_systemdependent.h"

// TODO(hkuang): Don't need to lock the whole pool after implementing atomic
// frame reference count.
void lock_buffer_pool(BufferPool *const pool) {
#if CONFIG_MULTITHREAD
  pthread_mutex_lock(&pool->pool_mutex);
#else
  (void)pool;
#endif
}

void unlock_buffer_pool(BufferPool *const pool) {
#if CONFIG_MULTITHREAD
  pthread_mutex_unlock(&pool->pool_mutex);
#else
  (void)pool;
#endif
}

static INLINE void alloc_mi_array(VP9_COMMON *cm, int mi_size, int idx) {
  CHECK_MEM_ERROR(cm, cm->mip_array[idx],
                  vpx_calloc(mi_size, sizeof(*cm->mip_array[0])));
  CHECK_MEM_ERROR(cm, cm->mi_grid_base_array[idx],
                  vpx_calloc(mi_size, sizeof(*cm->mi_grid_base_array[0])));
}

static void clear_mi_border(const VP9_COMMON *cm, MODE_INFO *mi) {
  int i;

  // Top border row
  vpx_memset(mi, 0, sizeof(*mi) * cm->mi_stride);

  // Left border column
  for (i = 1; i < cm->mi_rows + 1; ++i)
    vpx_memset(&mi[i * cm->mi_stride], 0, sizeof(*mi));
}

static void set_mb_mi(VP9_COMMON *cm, int aligned_width, int aligned_height) {
  cm->mi_cols = aligned_width >> MI_SIZE_LOG2;
  cm->mi_rows = aligned_height >> MI_SIZE_LOG2;
  cm->mi_stride = cm->mi_cols + MI_BLOCK_SIZE;

  cm->mb_cols = (cm->mi_cols + 1) >> 1;
  cm->mb_rows = (cm->mi_rows + 1) >> 1;
  cm->MBs = cm->mb_rows * cm->mb_cols;
}

static void setup_mi(VP9_COMMON *cm) {
  cm->mi = cm->mip + cm->mi_stride + 1;
  cm->prev_mi = cm->prev_mip + cm->mi_stride + 1;
  cm->mi_grid_visible = cm->mi_grid_base + cm->mi_stride + 1;
  cm->prev_mi_grid_visible = cm->prev_mi_grid_base + cm->mi_stride + 1;

  vpx_memset(cm->mip, 0, cm->mi_stride * (cm->mi_rows + 1) * sizeof(*cm->mip));

  vpx_memset(cm->mi_grid_base, 0, cm->mi_stride * (cm->mi_rows + 1) *
                                      sizeof(*cm->mi_grid_base));

  // Only clear mi border in non frame-parallel decode. In frame-parallel
  // decode, prev_mip is managed by previous decoding thread. While in
  // non frame-parallel decode, prev_mip and mip are both managed by
  // current decoding thread.
  if (!cm->frame_parallel_decode)
    clear_mi_border(cm, cm->prev_mip);
}

static int alloc_mi(VP9_COMMON *cm, int mi_size) {
  int i;

  for (i = 0; i < NUM_PING_PONG_BUFFERS; ++i) {
    // Delay reallocation as another thread is accessing prev_mi.
    if (cm->frame_parallel_decode && i == cm->prev_mi_idx) {
      cm->update_prev_mi = 1;
      continue;
    }
    alloc_mi_array(cm, mi_size, i);
  }

  cm->mip = cm->mip_array[cm->mi_idx];
  cm->mi_grid_base = cm->mi_grid_base_array[cm->mi_idx];

  if (!cm->frame_parallel_decode) {
    cm->mi_idx = 0;
    cm->prev_mi_idx = 1;
    // In frame-parallel decode, prev_mip comes from another thread,
    // so current decoding thread should not touch it.
    cm->prev_mip = cm->mip_array[cm->prev_mi_idx];
    cm->prev_mi_grid_base = cm->mi_grid_base_array[cm->prev_mi_idx];
  }

  return 0;
}

static void free_mi(VP9_COMMON *cm, int decode_done) {
  int i;

  for (i = 0; i < NUM_PING_PONG_BUFFERS; ++i) {
    if (cm->frame_parallel_decode && i == cm->prev_mi_idx && !decode_done)
      continue;
    vpx_free(cm->mip_array[i]);
    cm->mip_array[i] = NULL;
    vpx_free(cm->mi_grid_base_array[i]);
    cm->mi_grid_base_array[i] = NULL;
  }

  cm->mip = NULL;
  cm->mi_grid_base = NULL;

  if (!cm->frame_parallel_decode) {
    cm->prev_mip = NULL;
    cm->prev_mi_grid_base = NULL;
  }
}

static int alloc_seg_map(VP9_COMMON *cm, int seg_map_size) {
  int i;

  for (i = 0; i < NUM_PING_PONG_BUFFERS; ++i) {
    cm->seg_map_array[i] = (uint8_t *)vpx_calloc(seg_map_size, 1);
    if (cm->seg_map_array[i] == NULL)
      return 1;
  }

  // Init the index.
  cm->seg_map_idx = 0;
  cm->prev_seg_map_idx = 1;

  cm->current_frame_seg_map = cm->seg_map_array[cm->seg_map_idx];

  if (!cm->frame_parallel_decode) {
    cm->last_frame_seg_map = cm->seg_map_array[cm->prev_seg_map_idx];
  }

  return 0;
}

static void free_seg_map(VP9_COMMON *cm) {
  int i;

  for (i = 0; i < NUM_PING_PONG_BUFFERS; ++i) {
    vpx_free(cm->seg_map_array[i]);
    cm->seg_map_array[i] = NULL;
  }

  cm->current_frame_seg_map = NULL;

  if (!cm->frame_parallel_decode) {
    cm->last_frame_seg_map = NULL;
  }
}

void vp9_free_frame_buffers(VP9_COMMON *cm) {
  int i;
  BufferPool *const pool = cm->buffer_pool;

  for (i = 0; i < FRAME_BUFFERS; ++i) {
    vp9_free_frame_buffer(&pool->frame_bufs[i].buf);

    if (pool->frame_bufs[i].ref_count > 0 &&
        pool->frame_bufs[i].raw_frame_buffer.data != NULL) {
      pool->release_fb_cb(pool->cb_priv, &pool->frame_bufs[i].raw_frame_buffer);
      pool->frame_bufs[i].ref_count = 0;
    }
  }

  vp9_free_frame_buffer(&cm->post_proc_buffer);
}

void vp9_free_context_buffers(VP9_COMMON *cm) {
  free_mi(cm, 1);
  free_seg_map(cm);

  vpx_free(cm->above_context);
  cm->above_context = NULL;

  vpx_free(cm->above_seg_context);
  cm->above_seg_context = NULL;
}

int vp9_resize_frame_buffers(VP9_COMMON *cm, int width, int height) {
  const int aligned_width = ALIGN_POWER_OF_TWO(width, MI_SIZE_LOG2);
  const int aligned_height = ALIGN_POWER_OF_TWO(height, MI_SIZE_LOG2);
#if CONFIG_INTERNAL_STATS || CONFIG_VP9_POSTPROC
  const int ss_x = cm->subsampling_x;
  const int ss_y = cm->subsampling_y;

  // TODO(agrange): this should be conditionally allocated.
  if (vp9_realloc_frame_buffer(&cm->post_proc_buffer, width, height, ss_x, ss_y,
                               VP9_DEC_BORDER_IN_PIXELS, NULL, NULL, NULL) < 0)
    goto fail;
#endif

  set_mb_mi(cm, aligned_width, aligned_height);

  free_mi(cm, 0);
  if (alloc_mi(cm, cm->mi_stride * (cm->mi_rows + MI_BLOCK_SIZE)))
    goto fail;

  setup_mi(cm);

  // Create the segmentation map structure and set to 0.
  free_seg_map(cm);
  if (alloc_seg_map(cm, cm->mi_rows * cm->mi_cols))
    goto fail;

  vpx_free(cm->above_context);
  cm->above_context =
      (ENTROPY_CONTEXT *)vpx_calloc(2 * mi_cols_aligned_to_sb(cm->mi_cols) *
                                        MAX_MB_PLANE,
                                    sizeof(*cm->above_context));
  if (!cm->above_context)
    goto fail;

  vpx_free(cm->above_seg_context);
  cm->above_seg_context =
     (PARTITION_CONTEXT *)vpx_calloc(mi_cols_aligned_to_sb(cm->mi_cols),
                                     sizeof(*cm->above_seg_context));
  if (!cm->above_seg_context)
    goto fail;

  return 0;

 fail:
  vp9_free_frame_buffers(cm);
  vp9_free_context_buffers(cm);
  return 1;
}

static void init_frame_bufs(VP9_COMMON *cm) {
  int i;
  RefCntBuffer *const frame_bufs = cm->buffer_pool->frame_bufs;

  cm->new_fb_idx = FRAME_BUFFERS - 1;
  frame_bufs[cm->new_fb_idx].ref_count = 1;

  for (i = 0; i < REF_FRAMES; ++i) {
    cm->ref_frame_map[i] = i;
    frame_bufs[i].ref_count = 1;
  }
}

int vp9_alloc_frame_buffers(VP9_COMMON *cm, int width, int height) {
  int i;
  const int ss_x = cm->subsampling_x;
  const int ss_y = cm->subsampling_y;
  RefCntBuffer *const frame_bufs = cm->buffer_pool->frame_bufs;

  vp9_free_frame_buffers(cm);

  for (i = 0; i < FRAME_BUFFERS; ++i) {
    frame_bufs[i].ref_count = 0;
    if (vp9_alloc_frame_buffer(&frame_bufs[i].buf, width, height,
                               ss_x, ss_y, VP9_ENC_BORDER_IN_PIXELS) < 0)
      goto fail;
  }

  init_frame_bufs(cm);

#if CONFIG_INTERNAL_STATS || CONFIG_VP9_POSTPROC
  if (vp9_alloc_frame_buffer(&cm->post_proc_buffer, width, height, ss_x, ss_y,
                             VP9_ENC_BORDER_IN_PIXELS) < 0)
    goto fail;
#endif

  return 0;

 fail:
  vp9_free_frame_buffers(cm);
  return 1;
}

int vp9_alloc_context_buffers(VP9_COMMON *cm, int width, int height) {
  const int aligned_width = ALIGN_POWER_OF_TWO(width, MI_SIZE_LOG2);
  const int aligned_height = ALIGN_POWER_OF_TWO(height, MI_SIZE_LOG2);

  vp9_free_context_buffers(cm);

  set_mb_mi(cm, aligned_width, aligned_height);

  if (alloc_mi(cm, cm->mi_stride * (cm->mi_rows + MI_BLOCK_SIZE)))
    goto fail;

  setup_mi(cm);

  // Create the segmentation map structure and set to 0.
  cm->last_frame_seg_map = (uint8_t *)vpx_calloc(cm->mi_rows * cm->mi_cols, 1);
  if (!cm->last_frame_seg_map)
    goto fail;

  cm->above_context =
      (ENTROPY_CONTEXT *)vpx_calloc(2 * mi_cols_aligned_to_sb(cm->mi_cols) *
                                        MAX_MB_PLANE,
                                    sizeof(*cm->above_context));
  if (!cm->above_context)
    goto fail;

  cm->above_seg_context =
      (PARTITION_CONTEXT *)vpx_calloc(mi_cols_aligned_to_sb(cm->mi_cols),
                                      sizeof(*cm->above_seg_context));
  if (!cm->above_seg_context)
    goto fail;

  return 0;

 fail:
  vp9_free_context_buffers(cm);
  return 1;
}

void vp9_remove_common(VP9_COMMON *cm) {
  vp9_free_frame_buffers(cm);
  vp9_free_context_buffers(cm);
}

void vp9_update_frame_size(VP9_COMMON *cm) {
  const int aligned_width = ALIGN_POWER_OF_TWO(cm->width, MI_SIZE_LOG2);
  const int aligned_height = ALIGN_POWER_OF_TWO(cm->height, MI_SIZE_LOG2);

  set_mb_mi(cm, aligned_width, aligned_height);
  setup_mi(cm);

  // Initialize the previous frame segment map to 0.
  if (cm->current_frame_seg_map)
    vpx_memset(cm->current_frame_seg_map, 0, cm->mi_rows * cm->mi_cols);
}

void vp9_swap_mi_and_prev_mi(VP9_COMMON *cm) {
  // Swap indices.
  const int tmp = cm->mi_idx;

  // Only used in frame parallel decode: Update the prev_mi buffer if
  // needed. The worker that was accessing it must already finish decoding.
  // So it can be resized safely now.
  if (cm->update_prev_mi) {
    const int mi_size = cm->mi_stride * (cm->mi_rows + MI_BLOCK_SIZE);
    vpx_free(cm->mip_array[cm->prev_mi_idx]);
    vpx_free(cm->mi_grid_base_array[cm->prev_mi_idx]);
    cm->mip_array[cm->prev_mi_idx] = NULL;
    cm->mi_grid_base_array[cm->prev_mi_idx] = NULL;
    alloc_mi_array(cm, mi_size, cm->prev_mi_idx);
    cm->update_prev_mi = 0;
  }

  cm->mi_idx = cm->prev_mi_idx;
  cm->prev_mi_idx = tmp;

  // Current mip will be the prev_mip for the next frame.
  cm->mip = cm->mip_array[cm->mi_idx];
  cm->prev_mip = cm->mip_array[cm->prev_mi_idx];
  cm->mi_grid_base = cm->mi_grid_base_array[cm->mi_idx];
  cm->prev_mi_grid_base = cm->mi_grid_base_array[cm->prev_mi_idx];

  // Update the upper left visible macroblock ptrs.
  cm->mi = cm->mip + cm->mi_stride + 1;
  cm->prev_mi = cm->prev_mip + cm->mi_stride + 1;
  cm->mi_grid_visible = cm->mi_grid_base + cm->mi_stride + 1;
  cm->prev_mi_grid_visible = cm->prev_mi_grid_base + cm->mi_stride + 1;
}

void vp9_swap_current_and_last_seg_map(VP9_COMMON *cm) {
  // Swap indices.
  const int tmp = cm->seg_map_idx;
  cm->seg_map_idx = cm->prev_seg_map_idx;
  cm->prev_seg_map_idx = tmp;

  cm->current_frame_seg_map = cm->seg_map_array[cm->seg_map_idx];
  cm->last_frame_seg_map = cm->seg_map_array[cm->prev_seg_map_idx];
}
