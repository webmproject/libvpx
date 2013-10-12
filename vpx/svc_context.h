/**
 * SvcContext - input parameters and state to encode a multi-layered
 * spatial SVC frame
 */

#ifndef VPX_SVC_CONTEXT_H_
#define VPX_SVC_CONTEXT_H_

typedef enum SVC_ENCODING_MODE {
  INTER_LAYER_PREDICTION_I,
  ALT_INTER_LAYER_PREDICTION_IP,
  INTER_LAYER_PREDICTION_IP,
  USE_GOLDEN_FRAME
} SVC_ENCODING_MODE;

typedef enum SVC_LOG_LEVEL {
  SVC_LOG_ERROR,
  SVC_LOG_INFO,
  SVC_LOG_DEBUG,
} SVC_LOG_LEVEL;

typedef struct {
  // public interface to svc_command options
  int enabled;                      // set to non-zero to enable svc encoding
  int spatial_layers;               // number of layers
  int first_frame_full_size;        // set to one to force first frame full size
  SVC_ENCODING_MODE encoding_mode;  // svc encoding strategy

  // the following lists are ordered from highest resolution to lowest
  // if the strings are null, default values are used
  const char* quantizer_values;  // quantizer values, e.g., "27,33,39,53,60"
  const char* scale_factors;     // layer scale factors, e.g.,
                                 // "16/16,11/16,7/16,5/16,4/16"
  int gop_size;                  // distance between keyframes

  SVC_LOG_LEVEL log_level;  // amount of information to display
  int log_print;  // when set, printf log messages instead of returning the
                  // message with svc_get_message

  // private storage for vpx_svc_encode
  void* internal;
} SvcContext;

vpx_codec_err_t vpx_svc_init(SvcContext* svc_ctx, vpx_codec_ctx_t* codec_ctx,
                             vpx_codec_iface_t* iface,
                             vpx_codec_enc_cfg_t* cfg);

vpx_codec_err_t vpx_svc_encode(SvcContext* svc_ctx, vpx_codec_ctx_t* codec_ctx,
                               struct vpx_image* rawimg, vpx_codec_pts_t pts,
                               int64_t duration, int deadline);

void svc_dump_statistics(SvcContext* svc_ctx);
char* svc_get_message(SvcContext* svc_ctx);
void* svc_get_buffer(SvcContext* svc_ctx);
vpx_codec_err_t svc_get_layer_resolution(SvcContext* svc_ctx, int layer,
                              unsigned int* width, unsigned int* height);
int svc_get_frame_size(SvcContext* svc_ctx);
int svc_get_encode_frame_count(SvcContext* svc_ctx);
int svc_is_keyframe(SvcContext* svc_ctx);
void svc_set_keyframe(SvcContext* svc_ctx);

#endif /* VPX_SVC_CONTEXT_H_ */
