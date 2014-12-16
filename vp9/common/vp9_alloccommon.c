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
#include "vp9/common/vp9_entropymode.h"
#include "vp9/common/vp9_entropymv.h"
#include "vp9/common/vp9_onyxc_int.h"
#include "vp9/common/vp9_systemdependent.h"

void vp9_set_mb_mi(VP9_COMMON *cm, int width, int height) {
  const int aligned_width = ALIGN_POWER_OF_TWO(width, MI_SIZE_LOG2);
  const int aligned_height = ALIGN_POWER_OF_TWO(height, MI_SIZE_LOG2);

  cm->mi_cols = aligned_width >> MI_SIZE_LOG2;
  cm->mi_rows = aligned_height >> MI_SIZE_LOG2;
  cm->mi_stride = calc_mi_size(cm->mi_cols);

  cm->mb_cols = (cm->mi_cols + 1) >> 1;
  cm->mb_rows = (cm->mi_rows + 1) >> 1;
  cm->MBs = cm->mb_rows * cm->mb_cols;
}

void vp9_free_ref_frame_buffers(VP9_COMMON *cm) {
  int i;

  for (i = 0; i < FRAME_BUFFERS; ++i) {
    if (cm->frame_bufs[i].ref_count > 0 &&
        cm->frame_bufs[i].raw_frame_buffer.data != NULL) {
      cm->release_fb_cb(cm->cb_priv, &cm->frame_bufs[i].raw_frame_buffer);
      cm->frame_bufs[i].ref_count = 0;
    }
    vpx_free(cm->frame_bufs[i].mvs);
    cm->frame_bufs[i].mvs = NULL;
    vp9_free_frame_buffer(&cm->frame_bufs[i].buf);
  }

#if CONFIG_VP9_POSTPROC
  vp9_free_frame_buffer(&cm->post_proc_buffer);
  vp9_free_frame_buffer(&cm->post_proc_buffer_int);
#endif
}

void vp9_free_context_buffers(VP9_COMMON *cm) {
  cm->free_mi(cm);
  vpx_free(cm->last_frame_seg_map);
  cm->last_frame_seg_map = NULL;
  vpx_free(cm->above_context);
  cm->above_context = NULL;
  vpx_free(cm->above_seg_context);
  cm->above_seg_context = NULL;
}

int vp9_alloc_context_buffers(VP9_COMMON *cm, int width, int height) {
  vp9_free_context_buffers(cm);

  vp9_set_mb_mi(cm, width, height);
  if (cm->alloc_mi(cm, cm->mi_stride * calc_mi_size(cm->mi_rows)))
    goto fail;

  cm->last_frame_seg_map = (uint8_t *)vpx_calloc(cm->mi_rows * cm->mi_cols, 1);
  if (!cm->last_frame_seg_map) goto fail;

  cm->above_context = (ENTROPY_CONTEXT *)vpx_calloc(
      2 * mi_cols_aligned_to_sb(cm->mi_cols) * MAX_MB_PLANE,
      sizeof(*cm->above_context));
  if (!cm->above_context) goto fail;

  cm->above_seg_context = (PARTITION_CONTEXT *)vpx_calloc(
      mi_cols_aligned_to_sb(cm->mi_cols), sizeof(*cm->above_seg_context));
  if (!cm->above_seg_context) goto fail;

  return 0;

 fail:
  vp9_free_context_buffers(cm);
  return 1;
}

static void init_frame_bufs(VP9_COMMON *cm) {
  int i;

  cm->new_fb_idx = FRAME_BUFFERS - 1;
  cm->frame_bufs[cm->new_fb_idx].ref_count = 1;

  for (i = 0; i < REF_FRAMES; ++i) {
    cm->ref_frame_map[i] = i;
    cm->frame_bufs[i].ref_count = 1;
  }
}

int vp9_alloc_ref_frame_buffers(VP9_COMMON *cm, int width, int height) {
  int i;
  const int ss_x = cm->subsampling_x;
  const int ss_y = cm->subsampling_y;

  vp9_free_ref_frame_buffers(cm);

  for (i = 0; i < FRAME_BUFFERS; ++i) {
    cm->frame_bufs[i].ref_count = 0;
    if (vp9_alloc_frame_buffer(&cm->frame_bufs[i].buf, width, height,
                               ss_x, ss_y,
#if CONFIG_VP9_HIGHBITDEPTH
                               cm->use_highbitdepth,
#endif
                               VP9_ENC_BORDER_IN_PIXELS,
                               cm->byte_alignment) < 0)
      goto fail;
    if (cm->frame_bufs[i].mvs == NULL) {
      cm->frame_bufs[i].mvs =
          (MV_REF *)vpx_calloc(cm->mi_rows * cm->mi_cols,
                               sizeof(*cm->frame_bufs[i].mvs));
      if (cm->frame_bufs[i].mvs == NULL)
        goto fail;

      cm->frame_bufs[i].mi_rows = cm->mi_rows;
      cm->frame_bufs[i].mi_cols = cm->mi_cols;
    }
  }

  init_frame_bufs(cm);

#if CONFIG_VP9_POSTPROC
  if (vp9_alloc_frame_buffer(&cm->post_proc_buffer, width, height, ss_x, ss_y,
#if CONFIG_VP9_HIGHBITDEPTH
                             cm->use_highbitdepth,
#endif
                             VP9_ENC_BORDER_IN_PIXELS,
                             cm->byte_alignment) < 0)
    goto fail;
#endif

  return 0;

 fail:
  vp9_free_ref_frame_buffers(cm);
  return 1;
}

void vp9_remove_common(VP9_COMMON *cm) {
  vp9_free_ref_frame_buffers(cm);
  vp9_free_context_buffers(cm);
  vp9_free_internal_frame_buffers(&cm->int_frame_buffers);

  vpx_free(cm->fc);
  cm->fc = NULL;
  vpx_free(cm->frame_contexts);
  cm->frame_contexts = NULL;
}

void vp9_init_context_buffers(VP9_COMMON *cm) {
  cm->setup_mi(cm);
  if (cm->last_frame_seg_map)
    vpx_memset(cm->last_frame_seg_map, 0, cm->mi_rows * cm->mi_cols);
}
