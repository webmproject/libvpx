/*
 *  Copyright (c) 2020 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#ifndef VPX_VPX_VPX_EXT_RATECTRL_H_
#define VPX_VPX_VPX_EXT_RATECTRL_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "./vpx_integer.h"

typedef void *vpx_rc_model_t;

typedef struct vpx_rc_encodeframe_decision {
  int q_index;
} vpx_rc_encodeframe_decision_t;

typedef struct vpx_rc_encodeframe_info {
  int frame_type;
  int show_index;
  int coding_index;
} vpx_rc_encodeframe_info_t;

typedef struct vpx_rc_encodeframe_result {
  int64_t sse;
  int64_t bit_count;
  int64_t pixel_count;
} vpx_rc_encodeframe_result_t;

typedef struct vpx_rc_firstpass_stats {
  double (*frame_stats)[26];
  int num_frames;
} vpx_rc_firstpass_stats_t;

typedef struct vpx_rc_config {
  int frame_width;
  int frame_height;
  int show_frame_count;
  int target_bitrate_kbps;
} vpx_rc_config_t;

/*!\brief Create an external rate control model callback prototype
 *
 * This callback is invoked by the encoder to create an external rate control
 * model.
 *
 * \param[in]  priv               Callback's private data
 * \param[in]  ratectrl_config    Pointer to vpx_rc_config_t
 * \param[out] rate_ctrl_model_pt Pointer to vpx_rc_model_t
 */
typedef int (*vpx_rc_create_model_cb_fn_t)(
    void *priv, const vpx_rc_config_t *ratectrl_config,
    vpx_rc_model_t *rate_ctrl_model_pt);

/*!\brief Send first pass stats to the external rate control model callback
 * prototype
 *
 * This callback is invoked by the encoder to send first pass stats to the
 * external rate control model.
 *
 * \param[in]  rate_ctrl_model    rate control model
 * \param[in]  first_pass_stats   first pass stats
 */
typedef int (*vpx_rc_send_firstpass_stats_cb_fn_t)(
    vpx_rc_model_t rate_ctrl_model,
    const vpx_rc_firstpass_stats_t *first_pass_stats);

/*!\brief Receive encode frame decision callback prototype
 *
 * This callback is invoked by the encoder to receive encode frame decision from
 * the external rate control model.
 *
 * \param[in]  rate_ctrl_model    rate control model
 * \param[in]  encode_frame_info  information of the coding frame
 * \param[out] frame_decision     encode decision of the coding frame
 */
typedef int (*vpx_rc_get_encodeframe_decision_cb_fn_t)(
    vpx_rc_model_t rate_ctrl_model,
    const vpx_rc_encodeframe_info_t *encode_frame_info,
    vpx_rc_encodeframe_decision_t *frame_decision);

/*!\brief Update encode frame result callback prototype
 *
 * This callback is invoked by the encoder to update encode frame result to the
 * external rate control model.
 *
 * \param[in]  rate_ctrl_model     rate control model
 * \param[out] encode_frame_result encode result of the coding frame
 */
typedef int (*vpx_rc_update_encodeframe_result_cb_fn_t)(
    vpx_rc_model_t rate_ctrl_model,
    vpx_rc_encodeframe_result_t *encode_frame_result);

/*!\brief Delete the external rate control model callback prototype
 *
 * This callback is invoked by the encoder to delete the external rate control
 * model.
 *
 * \param[in]  rate_ctrl_model     rate control model
 */
typedef int (*vpx_rc_delete_model_cb_fn_t)(vpx_rc_model_t rate_ctrl_model);

typedef struct vpx_rc_funcs {
  vpx_rc_create_model_cb_fn_t create_model;
  vpx_rc_send_firstpass_stats_cb_fn_t send_firstpass_stats;
  vpx_rc_get_encodeframe_decision_cb_fn_t get_encodeframe_decision;
  vpx_rc_update_encodeframe_result_cb_fn_t update_encodeframe_result;
  vpx_rc_delete_model_cb_fn_t delete_model;
  void *priv;
} vpx_rc_funcs_t;

#ifdef __cplusplus
}  // extern "C"
#endif

#endif  // VPX_VPX_VPX_EXT_RATECTRL_H_
