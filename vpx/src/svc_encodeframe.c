/*
 *  Copyright (c) 2013 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

/**
 * @file
 * VP9 SVC encoding support via libvpx
 */

#include <assert.h>
#include <math.h>
#include <limits.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define VPX_DISABLE_CTRL_TYPECHECKS 1
#include "./vpx_config.h"
#include "vpx/svc_context.h"
#include "vpx/vp8cx.h"
#include "vpx/vpx_encoder.h"
#include "vpx_mem/vpx_mem.h"
#include "vp9/common/vp9_onyxc_int.h"

#ifdef __MINGW32__
#define strtok_r strtok_s
#ifndef MINGW_HAS_SECURE_API
// proto from /usr/x86_64-w64-mingw32/include/sec_api/string_s.h
_CRTIMP char *__cdecl strtok_s(char *str, const char *delim, char **context);
#endif  /* MINGW_HAS_SECURE_API */
#endif  /* __MINGW32__ */

#ifdef _MSC_VER
#define strdup _strdup
#define strtok_r strtok_s
#endif

#define SVC_REFERENCE_FRAMES 8
#define SUPERFRAME_SLOTS (8)
#define SUPERFRAME_BUFFER_SIZE (SUPERFRAME_SLOTS * sizeof(uint32_t) + 2)
#define OPTION_BUFFER_SIZE 256
#define COMPONENTS 4  // psnr & sse statistics maintained for total, y, u, v

static const int DEFAULT_QUANTIZER_VALUES[VPX_SS_MAX_LAYERS] = {
  60, 53, 39, 33, 27
};

static const int DEFAULT_SCALE_FACTORS_NUM[VPX_SS_MAX_LAYERS] = {
  4, 5, 7, 11, 16
};

static const int DEFAULT_SCALE_FACTORS_DEN[VPX_SS_MAX_LAYERS] = {
  16, 16, 16, 16, 16
};

typedef enum {
  QUANTIZER = 0,
  BITRATE,
  SCALE_FACTOR,
  AUTO_ALT_REF,
  ALL_OPTION_TYPES
} LAYER_OPTION_TYPE;

static const int option_max_values[ALL_OPTION_TYPES] = {
  63, INT_MAX, INT_MAX, 1
};

static const int option_min_values[ALL_OPTION_TYPES] = {
  0, 0, 1, 0
};

// One encoded frame
typedef struct FrameData {
  void                     *buf;    // compressed data buffer
  size_t                    size;  // length of compressed data
  vpx_codec_frame_flags_t   flags;    /**< flags for this frame */
  struct FrameData         *next;
} FrameData;

typedef struct SvcInternal {
  char options[OPTION_BUFFER_SIZE];        // set by vpx_svc_set_options
  char quantizers[OPTION_BUFFER_SIZE];     // set by vpx_svc_set_quantizers
  char scale_factors[OPTION_BUFFER_SIZE];  // set by vpx_svc_set_scale_factors

  // values extracted from option, quantizers
  int scaling_factor_num[VPX_SS_MAX_LAYERS];
  int scaling_factor_den[VPX_SS_MAX_LAYERS];
  int quantizer[VPX_SS_MAX_LAYERS];
  int enable_auto_alt_ref[VPX_SS_MAX_LAYERS];
  int bitrates[VPX_SS_MAX_LAYERS];

  // accumulated statistics
  double psnr_sum[VPX_SS_MAX_LAYERS][COMPONENTS];   // total/Y/U/V
  uint64_t sse_sum[VPX_SS_MAX_LAYERS][COMPONENTS];
  uint32_t bytes_sum[VPX_SS_MAX_LAYERS];

  // codec encoding values
  int width;    // width of highest layer
  int height;   // height of highest layer
  int kf_dist;  // distance between keyframes

  // state variables
  int encode_frame_count;
  int frame_received;
  int frame_within_gop;
  int layers;
  int layer;
  int is_keyframe;
  int use_multiple_frame_contexts;

  char message_buffer[2048];
  vpx_codec_ctx_t *codec_ctx;
} SvcInternal;

static SvcInternal *get_svc_internal(SvcContext *svc_ctx) {
  if (svc_ctx == NULL) return NULL;
  if (svc_ctx->internal == NULL) {
    SvcInternal *const si = (SvcInternal *)malloc(sizeof(*si));
    if (si != NULL) {
      memset(si, 0, sizeof(*si));
    }
    svc_ctx->internal = si;
  }
  return (SvcInternal *)svc_ctx->internal;
}

static const SvcInternal *get_const_svc_internal(const SvcContext *svc_ctx) {
  if (svc_ctx == NULL) return NULL;
  return (const SvcInternal *)svc_ctx->internal;
}

static void svc_log_reset(SvcContext *svc_ctx) {
  SvcInternal *const si = (SvcInternal *)svc_ctx->internal;
  si->message_buffer[0] = '\0';
}

static int svc_log(SvcContext *svc_ctx, SVC_LOG_LEVEL level,
                   const char *fmt, ...) {
  char buf[512];
  int retval = 0;
  va_list ap;
  SvcInternal *const si = get_svc_internal(svc_ctx);

  if (level > svc_ctx->log_level) {
    return retval;
  }

  va_start(ap, fmt);
  retval = vsnprintf(buf, sizeof(buf), fmt, ap);
  va_end(ap);

  if (svc_ctx->log_print) {
    printf("%s", buf);
  } else {
    strncat(si->message_buffer, buf,
            sizeof(si->message_buffer) - strlen(si->message_buffer) - 1);
  }

  if (level == SVC_LOG_ERROR) {
    si->codec_ctx->err_detail = si->message_buffer;
  }
  return retval;
}

static vpx_codec_err_t extract_option(LAYER_OPTION_TYPE type,
                                      char *input,
                                      int *value0,
                                      int *value1) {
  if (type == SCALE_FACTOR) {
    *value0 = strtol(input, &input, 10);
    if (*input++ != '/')
      return VPX_CODEC_INVALID_PARAM;
    *value1 = strtol(input, &input, 10);

    if (*value0 < option_min_values[SCALE_FACTOR] ||
        *value1 < option_min_values[SCALE_FACTOR] ||
        *value0 > option_max_values[SCALE_FACTOR] ||
        *value1 > option_max_values[SCALE_FACTOR])
      return VPX_CODEC_INVALID_PARAM;
  } else {
    *value0 = atoi(input);
    if (*value0 < option_min_values[type] ||
        *value0 > option_max_values[type])
      return VPX_CODEC_INVALID_PARAM;
  }
  return VPX_CODEC_OK;
}

static vpx_codec_err_t parse_layer_options_from_string(SvcContext *svc_ctx,
                                                       LAYER_OPTION_TYPE type,
                                                       const char *input,
                                                       int *option0,
                                                       int *option1) {
  int i;
  vpx_codec_err_t res = VPX_CODEC_OK;
  char *input_string;
  char *token;
  const char *delim = ",";
  char *save_ptr;

  if (input == NULL || option0 == NULL ||
      (option1 == NULL && type == SCALE_FACTOR))
    return VPX_CODEC_INVALID_PARAM;

  input_string = strdup(input);
  token = strtok_r(input_string, delim, &save_ptr);
  for (i = 0; i < svc_ctx->spatial_layers; ++i) {
    if (token != NULL) {
      res = extract_option(type, token, option0 + i, option1 + i);
      if (res != VPX_CODEC_OK)
        break;
      token = strtok_r(NULL, delim, &save_ptr);
    } else {
      break;
    }
  }
  if (res == VPX_CODEC_OK && i != svc_ctx->spatial_layers) {
    svc_log(svc_ctx, SVC_LOG_ERROR,
            "svc: layer params type: %d    %d values required, "
            "but only %d specified\n", type, svc_ctx->spatial_layers, i);
    res = VPX_CODEC_INVALID_PARAM;
  }
  free(input_string);
  return res;
}

/**
 * Parse SVC encoding options
 * Format: encoding-mode=<svc_mode>,layers=<layer_count>
 *         scale-factors=<n1>/<d1>,<n2>/<d2>,...
 *         quantizers=<q1>,<q2>,...
 * svc_mode = [i|ip|alt_ip|gf]
 */
static vpx_codec_err_t parse_options(SvcContext *svc_ctx, const char *options) {
  char *input_string;
  char *option_name;
  char *option_value;
  char *input_ptr;
  SvcInternal *const si = get_svc_internal(svc_ctx);
  vpx_codec_err_t res = VPX_CODEC_OK;
  int i, alt_ref_enabled = 0;

  if (options == NULL) return VPX_CODEC_OK;
  input_string = strdup(options);

  // parse option name
  option_name = strtok_r(input_string, "=", &input_ptr);
  while (option_name != NULL) {
    // parse option value
    option_value = strtok_r(NULL, " ", &input_ptr);
    if (option_value == NULL) {
      svc_log(svc_ctx, SVC_LOG_ERROR, "option missing value: %s\n",
              option_name);
      res = VPX_CODEC_INVALID_PARAM;
      break;
    }
    if (strcmp("spatial-layers", option_name) == 0) {
      svc_ctx->spatial_layers = atoi(option_value);
    } else if (strcmp("temporal-layers", option_name) == 0) {
      svc_ctx->temporal_layers = atoi(option_value);
    } else if (strcmp("scale-factors", option_name) == 0) {
      res = parse_layer_options_from_string(svc_ctx, SCALE_FACTOR, option_value,
                                            si->scaling_factor_num,
                                            si->scaling_factor_den);
      if (res != VPX_CODEC_OK) break;
    } else if (strcmp("quantizers", option_name) == 0) {
      res = parse_layer_options_from_string(svc_ctx, QUANTIZER, option_value,
                                            si->quantizer, NULL);
      if (res != VPX_CODEC_OK) break;
    } else if (strcmp("auto-alt-refs", option_name) == 0) {
      res = parse_layer_options_from_string(svc_ctx, AUTO_ALT_REF, option_value,
                                            si->enable_auto_alt_ref, NULL);
      if (res != VPX_CODEC_OK) break;
    } else if (strcmp("bitrates", option_name) == 0) {
      res = parse_layer_options_from_string(svc_ctx, BITRATE, option_value,
                                            si->bitrates, NULL);
      if (res != VPX_CODEC_OK) break;
    } else if (strcmp("multi-frame-contexts", option_name) == 0) {
      si->use_multiple_frame_contexts = atoi(option_value);
    } else {
      svc_log(svc_ctx, SVC_LOG_ERROR, "invalid option: %s\n", option_name);
      res = VPX_CODEC_INVALID_PARAM;
      break;
    }
    option_name = strtok_r(NULL, "=", &input_ptr);
  }
  free(input_string);

  if (si->use_multiple_frame_contexts &&
      (svc_ctx->spatial_layers > 3 ||
       svc_ctx->spatial_layers * svc_ctx->temporal_layers > 4))
    res = VPX_CODEC_INVALID_PARAM;

  for (i = 0; i < svc_ctx->spatial_layers; ++i)
    alt_ref_enabled += si->enable_auto_alt_ref[i];
  if (alt_ref_enabled > REF_FRAMES - svc_ctx->spatial_layers) {
    svc_log(svc_ctx, SVC_LOG_ERROR,
            "svc: auto alt ref: Maxinum %d(REF_FRAMES - layers) layers could"
            "enabled auto alt reference frame, but % layers are enabled\n",
            REF_FRAMES - svc_ctx->spatial_layers, alt_ref_enabled);
    res = VPX_CODEC_INVALID_PARAM;
  }

  return res;
}

vpx_codec_err_t vpx_svc_set_options(SvcContext *svc_ctx, const char *options) {
  SvcInternal *const si = get_svc_internal(svc_ctx);
  if (svc_ctx == NULL || options == NULL || si == NULL) {
    return VPX_CODEC_INVALID_PARAM;
  }
  strncpy(si->options, options, sizeof(si->options));
  si->options[sizeof(si->options) - 1] = '\0';
  return VPX_CODEC_OK;
}

vpx_codec_err_t vpx_svc_set_quantizers(SvcContext *svc_ctx,
                                       const char *quantizers) {
  SvcInternal *const si = get_svc_internal(svc_ctx);
  if (svc_ctx == NULL || quantizers == NULL || si == NULL) {
    return VPX_CODEC_INVALID_PARAM;
  }
  strncpy(si->quantizers, quantizers, sizeof(si->quantizers));
  si->quantizers[sizeof(si->quantizers) - 1] = '\0';
  return VPX_CODEC_OK;
}

vpx_codec_err_t vpx_svc_set_scale_factors(SvcContext *svc_ctx,
                                          const char *scale_factors) {
  SvcInternal *const si = get_svc_internal(svc_ctx);
  if (svc_ctx == NULL || scale_factors == NULL || si == NULL) {
    return VPX_CODEC_INVALID_PARAM;
  }
  strncpy(si->scale_factors, scale_factors, sizeof(si->scale_factors));
  si->scale_factors[sizeof(si->scale_factors) - 1] = '\0';
  return VPX_CODEC_OK;
}

void assign_layer_bitrates(const SvcInternal *const si,
                           vpx_codec_enc_cfg_t *const enc_cfg) {
  int i;

  if (si->bitrates[0] != 0) {
    enc_cfg->rc_target_bitrate = 0;
    for (i = 0; i < si->layers; ++i) {
      enc_cfg->ss_target_bitrate[i] = (unsigned int)si->bitrates[i];
      enc_cfg->rc_target_bitrate += si->bitrates[i];
    }
  } else {
    float total = 0;
    float alloc_ratio[VPX_SS_MAX_LAYERS] = {0};

    for (i = 0; i < si->layers; ++i) {
      if (si->scaling_factor_den[i] > 0) {
        alloc_ratio[i] = (float)(si->scaling_factor_num[i] * 1.0 /
            si->scaling_factor_den[i]);

        alloc_ratio[i] *= alloc_ratio[i];
        total += alloc_ratio[i];
      }
    }

    for (i = 0; i < si->layers; ++i) {
      if (total > 0) {
        enc_cfg->ss_target_bitrate[i] = (unsigned int)
            (enc_cfg->rc_target_bitrate * alloc_ratio[i] / total);
      }
    }
  }
}

vpx_codec_err_t vpx_svc_init(SvcContext *svc_ctx, vpx_codec_ctx_t *codec_ctx,
                             vpx_codec_iface_t *iface,
                             vpx_codec_enc_cfg_t *enc_cfg) {
  vpx_codec_err_t res;
  int i;
  SvcInternal *const si = get_svc_internal(svc_ctx);
  if (svc_ctx == NULL || codec_ctx == NULL || iface == NULL ||
      enc_cfg == NULL) {
    return VPX_CODEC_INVALID_PARAM;
  }
  if (si == NULL) return VPX_CODEC_MEM_ERROR;

  si->codec_ctx = codec_ctx;

  si->width = enc_cfg->g_w;
  si->height = enc_cfg->g_h;

  if (enc_cfg->kf_max_dist < 2) {
    svc_log(svc_ctx, SVC_LOG_ERROR, "key frame distance too small: %d\n",
            enc_cfg->kf_max_dist);
    return VPX_CODEC_INVALID_PARAM;
  }
  si->kf_dist = enc_cfg->kf_max_dist;

  if (svc_ctx->spatial_layers == 0)
    svc_ctx->spatial_layers = VPX_SS_DEFAULT_LAYERS;
  if (svc_ctx->spatial_layers < 1 ||
      svc_ctx->spatial_layers > VPX_SS_MAX_LAYERS) {
    svc_log(svc_ctx, SVC_LOG_ERROR, "spatial layers: invalid value: %d\n",
            svc_ctx->spatial_layers);
    return VPX_CODEC_INVALID_PARAM;
  }

  for (i = 0; i < VPX_SS_MAX_LAYERS; ++i) {
    si->quantizer[i] = DEFAULT_QUANTIZER_VALUES[i];
    si->scaling_factor_num[i] = DEFAULT_SCALE_FACTORS_NUM[i];
    si->scaling_factor_den[i] = DEFAULT_SCALE_FACTORS_DEN[i];
  }

  if (strlen(si->quantizers) > 0) {
    res = parse_layer_options_from_string(svc_ctx, QUANTIZER, si->quantizers,
                                          si->quantizer, NULL);
    if (res != VPX_CODEC_OK)
      return res;
  }

  if (strlen(si->scale_factors) > 0) {
    res = parse_layer_options_from_string(svc_ctx, SCALE_FACTOR,
                                          si->scale_factors,
                                          si->scaling_factor_num,
                                          si->scaling_factor_den);
    if (res != VPX_CODEC_OK)
      return res;
  }

  // Parse aggregate command line options. Options must start with
  // "layers=xx" then followed by other options
  res = parse_options(svc_ctx, si->options);
  if (res != VPX_CODEC_OK) return res;

  if (svc_ctx->spatial_layers < 1)
    svc_ctx->spatial_layers = 1;
  if (svc_ctx->spatial_layers > VPX_SS_MAX_LAYERS)
    svc_ctx->spatial_layers = VPX_SS_MAX_LAYERS;

  if (svc_ctx->temporal_layers < 1)
    svc_ctx->temporal_layers = 1;
  if (svc_ctx->temporal_layers > VPX_TS_MAX_LAYERS)
    svc_ctx->temporal_layers = VPX_TS_MAX_LAYERS;

  si->layers = svc_ctx->spatial_layers;

  assign_layer_bitrates(si, enc_cfg);

#if CONFIG_SPATIAL_SVC
  for (i = 0; i < si->layers; ++i)
    enc_cfg->ss_enable_auto_alt_ref[i] = si->enable_auto_alt_ref[i];
#endif

  if (svc_ctx->temporal_layers > 1) {
    int i;
    for (i = 0; i < svc_ctx->temporal_layers; ++i) {
      enc_cfg->ts_target_bitrate[i] = enc_cfg->rc_target_bitrate /
                                      svc_ctx->temporal_layers;
      enc_cfg->ts_rate_decimator[i] = 1 << (svc_ctx->temporal_layers - 1 - i);
    }
  }

  // modify encoder configuration
  enc_cfg->ss_number_layers = si->layers;
  enc_cfg->ts_number_layers = svc_ctx->temporal_layers;

  // TODO(ivanmaltz): determine if these values need to be set explicitly for
  // svc, or if the normal default/override mechanism can be used
  enc_cfg->rc_dropframe_thresh = 0;
  enc_cfg->rc_resize_allowed = 0;

  if (enc_cfg->g_pass == VPX_RC_ONE_PASS) {
    enc_cfg->rc_min_quantizer = 33;
    enc_cfg->rc_max_quantizer = 33;
  }

  enc_cfg->rc_undershoot_pct = 100;
  enc_cfg->rc_overshoot_pct = 15;
  enc_cfg->rc_buf_initial_sz = 500;
  enc_cfg->rc_buf_optimal_sz = 600;
  enc_cfg->rc_buf_sz = 1000;
  if (enc_cfg->g_error_resilient == 0 && si->use_multiple_frame_contexts == 0)
    enc_cfg->g_error_resilient = 1;

  // Initialize codec
  res = vpx_codec_enc_init(codec_ctx, iface, enc_cfg, VPX_CODEC_USE_PSNR);
  if (res != VPX_CODEC_OK) {
    svc_log(svc_ctx, SVC_LOG_ERROR, "svc_enc_init error\n");
    return res;
  }

  vpx_codec_control(codec_ctx, VP9E_SET_SVC, 1);
  vpx_codec_control(codec_ctx, VP8E_SET_TOKEN_PARTITIONS, 1);

  return VPX_CODEC_OK;
}

vpx_codec_err_t vpx_svc_get_layer_resolution(const SvcContext *svc_ctx,
                                             int layer,
                                             unsigned int *width,
                                             unsigned int *height) {
  int w, h, num, den;
  const SvcInternal *const si = get_const_svc_internal(svc_ctx);

  if (svc_ctx == NULL || si == NULL || width == NULL || height == NULL) {
    return VPX_CODEC_INVALID_PARAM;
  }
  if (layer < 0 || layer >= si->layers) return VPX_CODEC_INVALID_PARAM;

  num = si->scaling_factor_num[layer];
  den = si->scaling_factor_den[layer];
  if (num == 0 || den == 0) return VPX_CODEC_INVALID_PARAM;

  w = si->width * num / den;
  h = si->height * num / den;

  // make height and width even to make chrome player happy
  w += w % 2;
  h += h % 2;

  *width = w;
  *height = h;

  return VPX_CODEC_OK;
}

static void set_svc_parameters(SvcContext *svc_ctx,
                               vpx_codec_ctx_t *codec_ctx) {
  int layer;
  vpx_svc_parameters_t svc_params;
  SvcInternal *const si = get_svc_internal(svc_ctx);

  memset(&svc_params, 0, sizeof(svc_params));
  svc_params.temporal_layer = 0;
  svc_params.spatial_layer = si->layer;

  layer = si->layer;
  if (VPX_CODEC_OK != vpx_svc_get_layer_resolution(svc_ctx, layer,
                                                   &svc_params.width,
                                                   &svc_params.height)) {
    svc_log(svc_ctx, SVC_LOG_ERROR, "vpx_svc_get_layer_resolution failed\n");
  }

  if (codec_ctx->config.enc->g_pass == VPX_RC_ONE_PASS) {
    svc_params.min_quantizer = si->quantizer[layer];
    svc_params.max_quantizer = si->quantizer[layer];
  } else {
    svc_params.min_quantizer = codec_ctx->config.enc->rc_min_quantizer;
    svc_params.max_quantizer = codec_ctx->config.enc->rc_max_quantizer;
  }

  svc_params.distance_from_i_frame = si->frame_within_gop;
  vpx_codec_control(codec_ctx, VP9E_SET_SVC_PARAMETERS, &svc_params);
}

/**
 * Encode a frame into multiple layers
 * Create a superframe containing the individual layers
 */
vpx_codec_err_t vpx_svc_encode(SvcContext *svc_ctx, vpx_codec_ctx_t *codec_ctx,
                               struct vpx_image *rawimg, vpx_codec_pts_t pts,
                               int64_t duration, int deadline) {
  vpx_codec_err_t res;
  vpx_codec_iter_t iter;
  const vpx_codec_cx_pkt_t *cx_pkt;
  int layer_for_psnr = 0;
  SvcInternal *const si = get_svc_internal(svc_ctx);
  if (svc_ctx == NULL || codec_ctx == NULL || si == NULL) {
    return VPX_CODEC_INVALID_PARAM;
  }

  svc_log_reset(svc_ctx);

  si->layers = svc_ctx->spatial_layers;
  if (si->encode_frame_count == 0) {
    si->frame_within_gop = 0;
  }
  si->is_keyframe = (si->frame_within_gop == 0);

  if (rawimg != NULL) {
    svc_log(svc_ctx, SVC_LOG_DEBUG,
            "vpx_svc_encode  layers: %d, frame_count: %d, "
            "frame_within_gop: %d\n", si->layers, si->encode_frame_count,
            si->frame_within_gop);
  }

  if (rawimg != NULL) {
    // encode each layer
    for (si->layer = 0; si->layer < si->layers; ++si->layer) {
      set_svc_parameters(svc_ctx, codec_ctx);
    }
  }

  res = vpx_codec_encode(codec_ctx, rawimg, pts, (uint32_t)duration, 0,
                         deadline);
  if (res != VPX_CODEC_OK) {
    return res;
  }
  // save compressed data
  iter = NULL;
  while ((cx_pkt = vpx_codec_get_cx_data(codec_ctx, &iter))) {
    switch (cx_pkt->kind) {
      case VPX_CODEC_PSNR_PKT: {
        int i;
        svc_log(svc_ctx, SVC_LOG_DEBUG,
                "SVC frame: %d, layer: %d, PSNR(Total/Y/U/V): "
                "%2.3f  %2.3f  %2.3f  %2.3f \n",
                si->frame_received, layer_for_psnr,
                cx_pkt->data.psnr.psnr[0], cx_pkt->data.psnr.psnr[1],
                cx_pkt->data.psnr.psnr[2], cx_pkt->data.psnr.psnr[3]);
        svc_log(svc_ctx, SVC_LOG_DEBUG,
                "SVC frame: %d, layer: %d, SSE(Total/Y/U/V): "
                "%2.3f  %2.3f  %2.3f  %2.3f \n",
                si->frame_received, layer_for_psnr,
                cx_pkt->data.psnr.sse[0], cx_pkt->data.psnr.sse[1],
                cx_pkt->data.psnr.sse[2], cx_pkt->data.psnr.sse[3]);
        for (i = 0; i < COMPONENTS; i++) {
          si->psnr_sum[layer_for_psnr][i] += cx_pkt->data.psnr.psnr[i];
          si->sse_sum[layer_for_psnr][i] += cx_pkt->data.psnr.sse[i];
        }
        ++layer_for_psnr;
        if (layer_for_psnr == svc_ctx->spatial_layers)
          layer_for_psnr = 0;
        break;
      }
#if CONFIG_SPATIAL_SVC
      case VPX_CODEC_SPATIAL_SVC_LAYER_SIZES: {
        int i;
        for (i = 0; i < si->layers; ++i)
          si->bytes_sum[i] += cx_pkt->data.layer_sizes[i];
        break;
      }
#endif
      default: {
        break;
      }
    }
  }

  if (rawimg != NULL) {
    ++si->frame_within_gop;
    ++si->encode_frame_count;
  }

  return VPX_CODEC_OK;
}

const char *vpx_svc_get_message(const SvcContext *svc_ctx) {
  const SvcInternal *const si = get_const_svc_internal(svc_ctx);
  if (svc_ctx == NULL || si == NULL) return NULL;
  return si->message_buffer;
}

int vpx_svc_get_encode_frame_count(const SvcContext *svc_ctx) {
  const SvcInternal *const si = get_const_svc_internal(svc_ctx);
  if (svc_ctx == NULL || si == NULL) return 0;
  return si->encode_frame_count;
}

void vpx_svc_set_keyframe(SvcContext *svc_ctx) {
  SvcInternal *const si = get_svc_internal(svc_ctx);
  if (svc_ctx == NULL || si == NULL) return;
  si->frame_within_gop = 0;
}

static double calc_psnr(double d) {
  if (d == 0) return 100;
  return -10.0 * log(d) / log(10.0);
}

// dump accumulated statistics and reset accumulated values
const char *vpx_svc_dump_statistics(SvcContext *svc_ctx) {
  int number_of_frames, encode_frame_count;
  int i, j;
  uint32_t bytes_total = 0;
  double scale[COMPONENTS];
  double psnr[COMPONENTS];
  double mse[COMPONENTS];
  double y_scale;

  SvcInternal *const si = get_svc_internal(svc_ctx);
  if (svc_ctx == NULL || si == NULL) return NULL;

  svc_log_reset(svc_ctx);

  encode_frame_count = si->encode_frame_count;
  if (si->encode_frame_count <= 0) return vpx_svc_get_message(svc_ctx);

  svc_log(svc_ctx, SVC_LOG_INFO, "\n");
  for (i = 0; i < si->layers; ++i) {
    number_of_frames = encode_frame_count;

    svc_log(svc_ctx, SVC_LOG_INFO,
            "Layer %d Average PSNR=[%2.3f, %2.3f, %2.3f, %2.3f], Bytes=[%u]\n",
            i, (double)si->psnr_sum[i][0] / number_of_frames,
            (double)si->psnr_sum[i][1] / number_of_frames,
            (double)si->psnr_sum[i][2] / number_of_frames,
            (double)si->psnr_sum[i][3] / number_of_frames, si->bytes_sum[i]);
    // the following psnr calculation is deduced from ffmpeg.c#print_report
    y_scale = si->width * si->height * 255.0 * 255.0 * number_of_frames;
    scale[1] = y_scale;
    scale[2] = scale[3] = y_scale / 4;  // U or V
    scale[0] = y_scale * 1.5;           // total

    for (j = 0; j < COMPONENTS; j++) {
      psnr[j] = calc_psnr(si->sse_sum[i][j] / scale[j]);
      mse[j] = si->sse_sum[i][j] * 255.0 * 255.0 / scale[j];
    }
    svc_log(svc_ctx, SVC_LOG_INFO,
            "Layer %d Overall PSNR=[%2.3f, %2.3f, %2.3f, %2.3f]\n", i, psnr[0],
            psnr[1], psnr[2], psnr[3]);
    svc_log(svc_ctx, SVC_LOG_INFO,
            "Layer %d Overall MSE=[%2.3f, %2.3f, %2.3f, %2.3f]\n", i, mse[0],
            mse[1], mse[2], mse[3]);

    bytes_total += si->bytes_sum[i];
    // clear sums for next time
    si->bytes_sum[i] = 0;
    for (j = 0; j < COMPONENTS; ++j) {
      si->psnr_sum[i][j] = 0;
      si->sse_sum[i][j] = 0;
    }
  }

  // only display statistics once
  si->encode_frame_count = 0;

  svc_log(svc_ctx, SVC_LOG_INFO, "Total Bytes=[%u]\n", bytes_total);
  return vpx_svc_get_message(svc_ctx);
}

void vpx_svc_release(SvcContext *svc_ctx) {
  SvcInternal *si;
  if (svc_ctx == NULL) return;
  // do not use get_svc_internal as it will unnecessarily allocate an
  // SvcInternal if it was not already allocated
  si = (SvcInternal *)svc_ctx->internal;
  if (si != NULL) {
    free(si);
    svc_ctx->internal = NULL;
  }
}

