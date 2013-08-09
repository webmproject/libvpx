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
      ptr->mbmi.mb_in_image = 1;
      ptr++;  // Next element in the row
    }

    // Step over border element at start of next row
    mi += cm->mode_info_stride;
  }
}

void vp9_free_frame_buffers(VP9_COMMON *oci) {
  int i;

  for (i = 0; i < NUM_YV12_BUFFERS; i++)
    vp9_free_frame_buffer(&oci->yv12_fb[i]);

  vp9_free_frame_buffer(&oci->post_proc_buffer);

  vpx_free(oci->mip);
  vpx_free(oci->prev_mip);
  vpx_free(oci->above_seg_context);
  vpx_free(oci->last_frame_seg_map);

  vpx_free(oci->above_context[0]);
  for (i = 0; i < MAX_MB_PLANE; i++)
    oci->above_context[i] = 0;
  oci->mip = NULL;
  oci->prev_mip = NULL;
  oci->above_seg_context = NULL;
  oci->last_frame_seg_map = NULL;
}

static void set_mb_mi(VP9_COMMON *cm, int aligned_width, int aligned_height) {
  cm->mb_cols = (aligned_width + 8) >> 4;
  cm->mb_rows = (aligned_height + 8) >> 4;
  cm->MBs = cm->mb_rows * cm->mb_cols;

  cm->mi_cols = aligned_width >> LOG2_MI_SIZE;
  cm->mi_rows = aligned_height >> LOG2_MI_SIZE;
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

int vp9_alloc_frame_buffers(VP9_COMMON *oci, int width, int height) {
  int i, mi_cols;

  const int aligned_width = ALIGN_POWER_OF_TWO(width, LOG2_MI_SIZE);
  const int aligned_height = ALIGN_POWER_OF_TWO(height, LOG2_MI_SIZE);
  const int ss_x = oci->subsampling_x;
  const int ss_y = oci->subsampling_y;
  int mi_size;

  vp9_free_frame_buffers(oci);

  for (i = 0; i < NUM_YV12_BUFFERS; i++) {
    oci->fb_idx_ref_cnt[i] = 0;
    if (vp9_alloc_frame_buffer(&oci->yv12_fb[i], width, height, ss_x, ss_y,
                               VP9BORDERINPIXELS) < 0)
      goto fail;
  }

  oci->new_fb_idx = NUM_YV12_BUFFERS - 1;
  oci->fb_idx_ref_cnt[oci->new_fb_idx] = 1;

  for (i = 0; i < ALLOWED_REFS_PER_FRAME; i++)
    oci->active_ref_idx[i] = i;

  for (i = 0; i < NUM_REF_FRAMES; i++) {
    oci->ref_frame_map[i] = i;
    oci->fb_idx_ref_cnt[i] = 1;
  }

  if (vp9_alloc_frame_buffer(&oci->post_proc_buffer, width, height, ss_x, ss_y,
                             VP9BORDERINPIXELS) < 0)
    goto fail;

  set_mb_mi(oci, aligned_width, aligned_height);

  // Allocation
  mi_size = oci->mode_info_stride * (oci->mi_rows + MI_BLOCK_SIZE);

  oci->mip = vpx_calloc(mi_size, sizeof(MODE_INFO));
  if (!oci->mip)
    goto fail;

  oci->prev_mip = vpx_calloc(mi_size, sizeof(MODE_INFO));
  if (!oci->prev_mip)
    goto fail;

  setup_mi(oci);

  // FIXME(jkoleszar): allocate subsampled arrays for U/V once subsampling
  // information is exposed at this level
  mi_cols = mi_cols_aligned_to_sb(oci->mi_cols);

  // 2 contexts per 'mi unit', so that we have one context per 4x4 txfm
  // block where mi unit size is 8x8.
# if CONFIG_ALPHA
  oci->above_context[0] = vpx_calloc(sizeof(ENTROPY_CONTEXT) * 8 * mi_cols, 1);
#else
  oci->above_context[0] = vpx_calloc(sizeof(ENTROPY_CONTEXT) * 6 * mi_cols, 1);
#endif
  if (!oci->above_context[0])
    goto fail;

  oci->above_seg_context = vpx_calloc(sizeof(PARTITION_CONTEXT) * mi_cols, 1);
  if (!oci->above_seg_context)
    goto fail;

  // Create the segmentation map structure and set to 0.
  oci->last_frame_seg_map = vpx_calloc(oci->mi_rows * oci->mi_cols, 1);
  if (!oci->last_frame_seg_map)
    goto fail;

  return 0;

 fail:
  vp9_free_frame_buffers(oci);
  return 1;
}

void vp9_create_common(VP9_COMMON *oci) {
  vp9_machine_specific_config(oci);

  vp9_init_mbmode_probs(oci);

  oci->tx_mode = ONLY_4X4;
  oci->comp_pred_mode = HYBRID_PREDICTION;

  // Initialize reference frame sign bias structure to defaults
  vpx_memset(oci->ref_frame_sign_bias, 0, sizeof(oci->ref_frame_sign_bias));
}

void vp9_remove_common(VP9_COMMON *oci) {
  vp9_free_frame_buffers(oci);
}

void vp9_initialize_common() {
  vp9_coef_tree_initialize();
  vp9_entropy_mode_init();
  vp9_entropy_mv_init();
}

void vp9_update_frame_size(VP9_COMMON *cm) {
  int i, mi_cols;
  const int aligned_width = ALIGN_POWER_OF_TWO(cm->width, LOG2_MI_SIZE);
  const int aligned_height = ALIGN_POWER_OF_TWO(cm->height, LOG2_MI_SIZE);

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
