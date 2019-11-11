#include <vector>
#include "vp9/common/vp9_onyxc_int.h"
#include "vp9/vp9_iface_common.h"
#include "vp9/encoder/vp9_encoder.h"
#include "vp9/encoder/vp9_firstpass.h"
#include "vp9/simple_encode.h"
#include "vp9/vp9_cx_iface.h"

// TODO(angiebird): Merge this function with vpx_img_plane_width()
static int img_plane_width(const vpx_image_t *img, int plane) {
  if (plane > 0 && img->x_chroma_shift > 0)
    return (img->d_w + 1) >> img->x_chroma_shift;
  else
    return img->d_w;
}

// TODO(angiebird): Merge this function with vpx_img_plane_height()
static int img_plane_height(const vpx_image_t *img, int plane) {
  if (plane > 0 && img->y_chroma_shift > 0)
    return (img->d_h + 1) >> img->y_chroma_shift;
  else
    return img->d_h;
}

// TODO(angiebird): Merge this function with vpx_img_read()
static int img_read(vpx_image_t *img, FILE *file) {
  int plane;

  for (plane = 0; plane < 3; ++plane) {
    unsigned char *buf = img->planes[plane];
    const int stride = img->stride[plane];
    const int w = img_plane_width(img, plane) *
                  ((img->fmt & VPX_IMG_FMT_HIGHBITDEPTH) ? 2 : 1);
    const int h = img_plane_height(img, plane);
    int y;

    for (y = 0; y < h; ++y) {
      if (fread(buf, 1, w, file) != (size_t)w) return 0;
      buf += stride;
    }
  }

  return 1;
}

class SimpleEncode::impl {
 public:
  VP9_COMP *cpi;
  vpx_img_fmt_t img_fmt;
  vpx_image_t tmp_img;
  std::vector<FIRSTPASS_STATS> frame_stats;
};

static VP9_COMP *init_encoder(const VP9EncoderConfig *oxcf,
                              vpx_img_fmt_t img_fmt) {
  VP9_COMP *cpi;
  BufferPool *buffer_pool = (BufferPool *)vpx_calloc(1, sizeof(*buffer_pool));
  vp9_initialize_enc();
  cpi = vp9_create_compressor(oxcf, buffer_pool);
  vp9_update_compressor_with_img_fmt(cpi, img_fmt);
  return cpi;
}

static void free_encoder(VP9_COMP *cpi) {
  vpx_free(cpi->common.buffer_pool);
  vp9_remove_compressor(cpi);
}

static INLINE vpx_rational_t make_vpx_rational(int num, int den) {
  vpx_rational_t v;
  v.num = num;
  v.den = den;
  return v;
}

SimpleEncode::SimpleEncode(int frame_width, int frame_height,
                           int frame_rate_num, int frame_rate_den,
                           int target_bitrate, int num_frames, FILE *file)
    : pimpl{ std::unique_ptr<impl>(new impl()) } {
  this->frame_width = frame_width;
  this->frame_height = frame_height;
  this->frame_rate_num = frame_rate_num;
  this->frame_rate_den = frame_rate_den;
  this->target_bitrate = target_bitrate;
  this->num_frames = num_frames;
  this->file = file;
  pimpl->cpi = NULL;
  pimpl->img_fmt = VPX_IMG_FMT_I420;
}

void SimpleEncode::ComputeFirstPassStats() {
  vpx_rational_t frame_rate = make_vpx_rational(frame_rate_num, frame_rate_den);
  const VP9EncoderConfig oxcf = vp9_get_encoder_config(
      frame_width, frame_height, frame_rate, target_bitrate, VPX_RC_FIRST_PASS);
  VP9_COMP *cpi = init_encoder(&oxcf, pimpl->img_fmt);
  struct lookahead_ctx *lookahead = cpi->lookahead;
  int i;
  int use_highbitdepth = 0;
#if CONFIG_VP9_HIGHBITDEPTH
  use_highbitdepth = cpi->common.use_highbitdepth;
#endif
  vpx_image_t img;
  vpx_img_alloc(&img, pimpl->img_fmt, frame_width, frame_height, 1);
  rewind(file);
  pimpl->frame_stats.clear();
  for (i = 0; i < num_frames; ++i) {
    assert(!vp9_lookahead_full(lookahead));
    if (img_read(&img, file)) {
      int next_show_idx = vp9_lookahead_next_show_idx(lookahead);
      int64_t ts_start =
          timebase_units_to_ticks(&oxcf.g_timebase_in_ts, next_show_idx);
      int64_t ts_end =
          timebase_units_to_ticks(&oxcf.g_timebase_in_ts, next_show_idx + 1);
      YV12_BUFFER_CONFIG sd;
      image2yuvconfig(&img, &sd);
      vp9_lookahead_push(lookahead, &sd, ts_start, ts_end, use_highbitdepth, 0);
      {
        int64_t time_stamp;
        int64_t time_end;
        int flush = 1;  // Make vp9_get_compressed_data process a frame
        size_t size;
        unsigned int frame_flags = 0;
        // TODO(angiebird): Call vp9_first_pass directly
        vp9_get_compressed_data(cpi, &frame_flags, &size, NULL, &time_stamp,
                                &time_end, flush);
        // vp9_get_compressed_data only generates first pass stats not
        // compresses data
        assert(size == 0);
      }
      pimpl->frame_stats.push_back(vp9_get_frame_stats(&cpi->twopass));
    }
  }
  vp9_end_first_pass(cpi);
  // TODO(angiebird): Store the total_stats apart form frame_stats
  pimpl->frame_stats.push_back(vp9_get_total_stats(&cpi->twopass));
  free_encoder(cpi);
  rewind(file);
  vpx_img_free(&img);
}

std::vector<std::vector<double>> SimpleEncode::ObserveFirstPassStats() {
  std::vector<std::vector<double>> output_stats;
  // TODO(angiebird): This function make several assumptions of
  // FIRSTPASS_STATS. 1) All elements in FIRSTPASS_STATS are double except the
  // last one. 2) The last entry of frame_stats is the total_stats.
  // Change the code structure, so that we don't have to make these assumptions

  // Note the last entry of frame_stats is the total_stats, we don't need it.
  for (size_t i = 0; i < pimpl->frame_stats.size() - 1; ++i) {
    double *buf_start = reinterpret_cast<double *>(&pimpl->frame_stats[i]);
    // We use - 1 here because the last member in FIRSTPASS_STATS is not double
    double *buf_end =
        buf_start + sizeof(pimpl->frame_stats[i]) / sizeof(*buf_end) - 1;
    std::vector<double> this_stats(buf_start, buf_end);
    output_stats.push_back(this_stats);
  }
  return output_stats;
}

void SimpleEncode::StartEncode() {
  assert(pimpl->frame_stats.size() > 0);
  vpx_rational_t frame_rate = make_vpx_rational(frame_rate_num, frame_rate_den);
  VP9EncoderConfig oxcf = vp9_get_encoder_config(
      frame_width, frame_height, frame_rate, target_bitrate, VPX_RC_LAST_PASS);
  vpx_fixed_buf_t stats;
  stats.buf = pimpl->frame_stats.data();
  stats.sz = sizeof(pimpl->frame_stats[0]) * pimpl->frame_stats.size();

  vp9_set_first_pass_stats(&oxcf, &stats);
  assert(pimpl->cpi == NULL);
  pimpl->cpi = init_encoder(&oxcf, pimpl->img_fmt);
  vpx_img_alloc(&pimpl->tmp_img, pimpl->img_fmt, frame_width, frame_height, 1);
  rewind(file);
}

void SimpleEncode::EndEncode() {
  free_encoder(pimpl->cpi);
  pimpl->cpi = nullptr;
  vpx_img_free(&pimpl->tmp_img);
  rewind(file);
}

void SimpleEncode::EncodeFrame(char *cx_data, size_t *size, size_t max_size) {
  VP9_COMP *cpi = pimpl->cpi;
  struct lookahead_ctx *lookahead = cpi->lookahead;
  int use_highbitdepth = 0;
#if CONFIG_VP9_HIGHBITDEPTH
  use_highbitdepth = cpi->common.use_highbitdepth;
#endif
  // The lookahead's size is set to oxcf->lag_in_frames.
  // We want to fill lookahead to it's max capacity if possible so that the
  // encoder can construct alt ref frame in time.
  // In the other words, we hope vp9_get_compressed_data to encode a frame
  // every time in the function
  while (!vp9_lookahead_full(lookahead)) {
    // TODO(angiebird): Check whether we can move this file read logics to
    // lookahead
    if (img_read(&pimpl->tmp_img, file)) {
      int next_show_idx = vp9_lookahead_next_show_idx(lookahead);
      int64_t ts_start =
          timebase_units_to_ticks(&cpi->oxcf.g_timebase_in_ts, next_show_idx);
      int64_t ts_end = timebase_units_to_ticks(&cpi->oxcf.g_timebase_in_ts,
                                               next_show_idx + 1);
      YV12_BUFFER_CONFIG sd;
      image2yuvconfig(&pimpl->tmp_img, &sd);
      vp9_lookahead_push(lookahead, &sd, ts_start, ts_end, use_highbitdepth, 0);
    } else {
      break;
    }
  }
  int64_t time_stamp;
  int64_t time_end;
  int flush = 1;  // Make vp9_get_compressed_data encode a frame
  unsigned int frame_flags = 0;
  vp9_get_compressed_data(cpi, &frame_flags, size,
                          reinterpret_cast<uint8_t *>(cx_data), &time_stamp,
                          &time_end, flush);
  // vp9_get_compressed_data is expected to encode a frame every time, so the
  // data size should be greater than zero.
  assert(*size > 0);
  if (*size >= max_size) {
    assert(0);
  }
}

SimpleEncode::~SimpleEncode() {}
