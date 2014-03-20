/*
 *  Copyright (c) 2014 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#ifndef VP9_ENCODER_VP9_SVC_LAYERCONTEXT_H_
#define VP9_ENCODER_VP9_SVC_LAYERCONTEXT_H_

#include "vpx/vpx_encoder.h"

#include "vp9/encoder/vp9_ratectrl.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
  RATE_CONTROL rc;
  int target_bandwidth;
  int64_t starting_buffer_level;
  int64_t optimal_buffer_level;
  int64_t maximum_buffer_size;
  double framerate;
  int avg_frame_size;
  struct twopass_rc twopass;
} LAYER_CONTEXT;

typedef struct {
  int spatial_layer_id;
  int temporal_layer_id;
  int number_spatial_layers;
  int number_temporal_layers;
  // Layer context used for rate control in temporal CBR mode or spatial
  // two pass mode. Defined for temporal or spatial layers for now.
  // Does not support temporal combined with spatial RC.
  LAYER_CONTEXT layer_context[MAX(VPX_TS_MAX_LAYERS, VPX_SS_MAX_LAYERS)];
} SVC;

struct VP9_COMP;

// Initialize layer context data from init_config().
void vp9_init_layer_context(struct VP9_COMP *const cpi);

// Update the layer context from a change_config() call.
void vp9_update_layer_context_change_config(struct VP9_COMP *const cpi,
                                            const int target_bandwidth);

// Prior to encoding the frame, update framerate-related quantities
// for the current layer.
void vp9_update_layer_framerate(struct VP9_COMP *const cpi);

// Prior to encoding the frame, set the layer context, for the current layer
// to be encoded, to the cpi struct.
void vp9_restore_layer_context(struct VP9_COMP *const cpi);

// Save the layer context after encoding the frame.
void vp9_save_layer_context(struct VP9_COMP *const cpi);

#ifdef __cplusplus
}  // extern "C"
#endif

#endif  // VP9_ENCODER_VP9_SVC_LAYERCONTEXT_
