/*
 *  Copyright (c) 2014 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#ifndef VP10_ENCODER_SVC_LAYERCONTEXT_H_
#define VP10_ENCODER_SVC_LAYERCONTEXT_H_

#include "vpx/vpx_encoder.h"

#include "vp10/encoder/ratectrl.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
  RATE_CONTROL rc;
  int target_bandwidth;
  int spatial_layer_target_bandwidth;  // Target for the spatial layer.
  double framerate;
  int avg_frame_size;
  int max_q;
  int min_q;
  int scaling_factor_num;
  int scaling_factor_den;
  TWO_PASS twopass;
  vpx_fixed_buf_t rc_twopass_stats_in;
  unsigned int current_video_frame_in_layer;
  int is_key_frame;
  int frames_from_key_frame;
  FRAME_TYPE last_frame_type;
  struct lookahead_entry  *alt_ref_source;
  int alt_ref_idx;
  int gold_ref_idx;
  int has_alt_frame;
  size_t layer_size;
  struct vpx_psnr_pkt psnr_pkt;
} LAYER_CONTEXT;

typedef struct {
  int spatial_layer_id;
  int temporal_layer_id;
  int number_spatial_layers;
  int number_temporal_layers;

  int spatial_layer_to_encode;

  // Workaround for multiple frame contexts
  enum {
    ENCODED = 0,
    ENCODING,
    NEED_TO_ENCODE
  }encode_empty_frame_state;
  struct lookahead_entry empty_frame;
  int encode_intra_empty_frame;

  // Store scaled source frames to be used for temporal filter to generate
  // a alt ref frame.
  YV12_BUFFER_CONFIG scaled_frames[MAX_LAG_BUFFERS];

  // Layer context used for rate control in one pass temporal CBR mode or
  // two pass spatial mode.
  LAYER_CONTEXT layer_context[VPX_MAX_LAYERS];
  // Indicates what sort of temporal layering is used.
  // Currently, this only works for CBR mode.
  VP9E_TEMPORAL_LAYERING_MODE temporal_layering_mode;
} SVC;

struct VP10_COMP;

// Initialize layer context data from init_config().
void vp10_init_layer_context(struct VP10_COMP *const cpi);

// Update the layer context from a change_config() call.
void vp10_update_layer_context_change_config(struct VP10_COMP *const cpi,
                                            const int target_bandwidth);

// Prior to encoding the frame, update framerate-related quantities
// for the current temporal layer.
void vp10_update_temporal_layer_framerate(struct VP10_COMP *const cpi);

// Update framerate-related quantities for the current spatial layer.
void vp10_update_spatial_layer_framerate(struct VP10_COMP *const cpi,
                                        double framerate);

// Prior to encoding the frame, set the layer context, for the current layer
// to be encoded, to the cpi struct.
void vp10_restore_layer_context(struct VP10_COMP *const cpi);

// Save the layer context after encoding the frame.
void vp10_save_layer_context(struct VP10_COMP *const cpi);

// Initialize second pass rc for spatial svc.
void vp10_init_second_pass_spatial_svc(struct VP10_COMP *cpi);

// Increment number of video frames in layer
void vp10_inc_frame_in_layer(struct VP10_COMP *const cpi);

// Check if current layer is key frame in spatial upper layer
int vp10_is_upper_layer_key_frame(const struct VP10_COMP *const cpi);

// Get the next source buffer to encode
struct lookahead_entry *vp10_svc_lookahead_pop(struct VP10_COMP *const cpi,
                                              struct lookahead_ctx *ctx,
                                              int drain);

// Start a frame and initialize svc parameters
int vp10_svc_start_frame(struct VP10_COMP *const cpi);

int vp10_one_pass_cbr_svc_start_layer(struct VP10_COMP *const cpi);

#ifdef __cplusplus
}  // extern "C"
#endif

#endif  // VP10_ENCODER_SVC_LAYERCONTEXT_
