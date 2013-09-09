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
#include "vp9/common/vp9_findnearmv.h"
#include "vp9/common/vp9_onyxc_int.h"
#include "vp9/common/vp9_systemdependent.h"

void vp9_update_mode_info_border(VP9_COMMON *cm, MODE_INFO *mi) {
  const int stride = cm->mode_info_stride;
  int i;

  // Clear down top border row
  vpx_memset(mi, 0, sizeof(MODE_INFO) * stride);

  // Clear left border column
  for (i = 1; i < cm->mi_rows + 1; i++)
    vpx_memset(&mi[i * stride], 0, sizeof(MODE_INFO));
}

void vp9_update_mode_info_in_image(VP9_COMMON *cm, MODE_INFO *mi) {
  int i, j;

  // For each in image mode_info element set the in image flag to 1
  for (i = 0; i < cm->mi_rows; i++) {
    MODE_INFO *ptr = mi;
    for (j = 0; j < cm->mi_cols; j++) {
      ptr->mbmi.in_image = 1;
      ptr++;  // Next element in the row
    }

    // Step over border element at start of next row
    mi += cm->mode_info_stride;
  }
}

void vp9_free_frame_buffers(VP9_COMMON *cm) {
  int i;

  for (i = 0; i < NUM_YV12_BUFFERS; i++)
    vp9_free_frame_buffer(&cm->yv12_fb[i]);

  vp9_free_frame_buffer(&cm->post_proc_buffer);

  vpx_free(cm->mip);
  vpx_free(cm->prev_mip);
  vpx_free(cm->above_seg_context);
  vpx_free(cm->last_frame_seg_map);

  vpx_free(cm->above_context[0]);
  for (i = 0; i < MAX_MB_PLANE; i++)
    cm->above_context[i] = 0;
  cm->mip = NULL;
  cm->prev_mip = NULL;
  cm->above_seg_context = NULL;
  cm->last_frame_seg_map = NULL;
}

static void set_mb_mi(VP9_COMMON *cm, int aligned_width, int aligned_height) {
  cm->mb_cols = (aligned_width + 8) >> 4;
  cm->mb_rows = (aligned_height + 8) >> 4;
  cm->MBs = cm->mb_rows * cm->mb_cols;

  cm->mi_cols = aligned_width >> MI_SIZE_LOG2;
  cm->mi_rows = aligned_height >> MI_SIZE_LOG2;
  cm->mode_info_stride = cm->mi_cols + MI_BLOCK_SIZE;
}

static void setup_mi(VP9_COMMON *cm) {
  cm->mi = cm->mip + cm->mode_info_stride + 1;
  cm->prev_mi = cm->prev_mip + cm->mode_info_stride + 1;

  vpx_memset(cm->mip, 0,
             cm->mode_info_stride * (cm->mi_rows + 1) * sizeof(MODE_INFO));

  vp9_update_mode_info_border(cm, cm->mip);
  vp9_update_mode_info_in_image(cm, cm->mi);

  vp9_update_mode_info_border(cm, cm->prev_mip);
  vp9_update_mode_info_in_image(cm, cm->prev_mi);
}

int vp9_alloc_frame_buffers(VP9_COMMON *cm, int width, int height) {
  int i, mi_cols;

  const int aligned_width = ALIGN_POWER_OF_TWO(width, MI_SIZE_LOG2);
  const int aligned_height = ALIGN_POWER_OF_TWO(height, MI_SIZE_LOG2);
  const int ss_x = cm->subsampling_x;
  const int ss_y = cm->subsampling_y;
  int mi_size;

  vp9_free_frame_buffers(cm);

  for (i = 0; i < NUM_YV12_BUFFERS; i++) {
    cm->fb_idx_ref_cnt[i] = 0;
    if (vp9_alloc_frame_buffer(&cm->yv12_fb[i], width, height, ss_x, ss_y,
                               VP9BORDERINPIXELS) < 0)
      goto fail;
  }

  cm->new_fb_idx = NUM_YV12_BUFFERS - 1;
  cm->fb_idx_ref_cnt[cm->new_fb_idx] = 1;

  for (i = 0; i < ALLOWED_REFS_PER_FRAME; i++)
    cm->active_ref_idx[i] = i;

  for (i = 0; i < NUM_REF_FRAMES; i++) {
    cm->ref_frame_map[i] = i;
    cm->fb_idx_ref_cnt[i] = 1;
  }

  if (vp9_alloc_frame_buffer(&cm->post_proc_buffer, width, height, ss_x, ss_y,
                             VP9BORDERINPIXELS) < 0)
    goto fail;

  set_mb_mi(cm, aligned_width, aligned_height);

  // Allocation
  mi_size = cm->mode_info_stride * (cm->mi_rows + MI_BLOCK_SIZE);

  cm->mip = vpx_calloc(mi_size, sizeof(MODE_INFO));
  if (!cm->mip)
    goto fail;

  cm->prev_mip = vpx_calloc(mi_size, sizeof(MODE_INFO));
  if (!cm->prev_mip)
    goto fail;

  setup_mi(cm);

  // FIXME(jkoleszar): allocate subsampled arrays for U/V once subsampling
  // information is exposed at this level
  mi_cols = mi_cols_aligned_to_sb(cm->mi_cols);

  // 2 contexts per 'mi unit', so that we have one context per 4x4 txfm
  // block where mi unit size is 8x8.
  cm->above_context[0] = vpx_calloc(sizeof(ENTROPY_CONTEXT) * MAX_MB_PLANE *
                                         (2 * mi_cols), 1);
  if (!cm->above_context[0])
    goto fail;

  cm->above_seg_context = vpx_calloc(sizeof(PARTITION_CONTEXT) * mi_cols, 1);
  if (!cm->above_seg_context)
    goto fail;

  // Create the segmentation map structure and set to 0.
  cm->last_frame_seg_map = vpx_calloc(cm->mi_rows * cm->mi_cols, 1);
  if (!cm->last_frame_seg_map)
    goto fail;

  return 0;

 fail:
  vp9_free_frame_buffers(cm);
  return 1;
}

void vp9_create_common(VP9_COMMON *cm) {
  vp9_machine_specific_config(cm);

  vp9_init_mbmode_probs(cm);

  cm->tx_mode = ONLY_4X4;
  cm->comp_pred_mode = HYBRID_PREDICTION;

  // Initialize reference frame sign bias structure to defaults
  vpx_memset(cm->ref_frame_sign_bias, 0, sizeof(cm->ref_frame_sign_bias));
}

void vp9_remove_common(VP9_COMMON *cm) {
  vp9_free_frame_buffers(cm);
}

void vp9_initialize_common() {
  vp9_coef_tree_initialize();
  vp9_entropy_mode_init();
  vp9_entropy_mv_init();
}

void vp9_update_frame_size(VP9_COMMON *cm) {
  int i, mi_cols;
  const int aligned_width = ALIGN_POWER_OF_TWO(cm->width, MI_SIZE_LOG2);
  const int aligned_height = ALIGN_POWER_OF_TWO(cm->height, MI_SIZE_LOG2);

  set_mb_mi(cm, aligned_width, aligned_height);
  setup_mi(cm);

  mi_cols = mi_cols_aligned_to_sb(cm->mi_cols);
  for (i = 1; i < MAX_MB_PLANE; i++)
    cm->above_context[i] =
        cm->above_context[0] + i * sizeof(ENTROPY_CONTEXT) * 2 * mi_cols;

  // Initialize the previous frame segment map to 0.
  if (cm->last_frame_seg_map)
    vpx_memset(cm->last_frame_seg_map, 0, cm->mi_rows * cm->mi_cols);
}
