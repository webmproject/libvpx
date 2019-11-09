#include "vp9/common/vp9_onyxc_int.h"
#include "vp9/encoder/vp9_encoder.h"
#include "vp9/simple_encode.h"
#include "vp9/vp9_cx_iface.h"

class SimpleEncode::impl {
 public:
  VP9_COMP *cpi;
  BufferPool *buffer_pool;
};

SimpleEncode::SimpleEncode(int frame_width, int frame_height,
                           vpx_rational_t frame_rate, int target_bitrate)
    : pimpl{ std::unique_ptr<impl>(new impl()) } {
  VP9EncoderConfig oxcf = vp9_get_encoder_config(
      frame_width, frame_height, frame_rate, target_bitrate, VPX_RC_LAST_PASS);
  pimpl->buffer_pool = (BufferPool *)vpx_calloc(1, sizeof(*pimpl->buffer_pool));
  vp9_initialize_enc();
  pimpl->cpi = vp9_create_compressor(&oxcf, pimpl->buffer_pool);
  vp9_update_compressor_with_img_fmt(pimpl->cpi, VPX_IMG_FMT_I420);
}

SimpleEncode::~SimpleEncode() {
  vpx_free(pimpl->buffer_pool);
  vp9_remove_compressor(pimpl->cpi);
  pimpl->cpi = nullptr;
}
