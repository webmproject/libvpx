/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include <stdlib.h>
#include <string.h>

#include "./vpx_version.h"

#include "vpx/internal/vpx_codec_internal.h"
#include "vpx/vp8dx.h"
#include "vpx/vpx_decoder.h"

#include "vp9/common/vp9_frame_buffers.h"

#include "vp9/decoder/vp9_decoder.h"
#include "vp9/decoder/vp9_read_bit_buffer.h"

#include "vp9/vp9_iface_common.h"

#define VP9_CAP_POSTPROC (CONFIG_VP9_POSTPROC ? VPX_CODEC_CAP_POSTPROC : 0)

typedef vpx_codec_stream_info_t vp9_stream_info_t;

struct vpx_codec_alg_priv {
  vpx_codec_priv_t        base;
  vpx_codec_dec_cfg_t     cfg;
  vp9_stream_info_t       si;
  int                     postproc_cfg_set;
  vp8_postproc_cfg_t      postproc_cfg;
  vpx_decrypt_cb          decrypt_cb;
  void                   *decrypt_state;
  vpx_image_t             img;
  int                     invert_tile_order;
  int                     frame_parallel_decode;  // frame-based threading.
  int                     last_show_frame;  // Index of last output frame.

  VP9Worker               *frame_workers;
  int                     num_frame_workers;
  int                     next_submit_thread_id;
  int                     next_output_thread_id;

  // External frame buffer info to save for VP9 common.
  void *ext_priv;  // Private data associated with the external frame buffers.
  vpx_get_frame_buffer_cb_fn_t get_ext_fb_cb;
  vpx_release_frame_buffer_cb_fn_t release_ext_fb_cb;
};

static vpx_codec_err_t decoder_init(vpx_codec_ctx_t *ctx,
                                    vpx_codec_priv_enc_mr_cfg_t *data) {
  // This function only allocates space for the vpx_codec_alg_priv_t
  // structure. More memory may be required at the time the stream
  // information becomes known.
  (void)data;

  if (!ctx->priv) {
    vpx_codec_alg_priv_t *alg_priv = vpx_memalign(32, sizeof(*alg_priv));
    if (alg_priv == NULL)
      return VPX_CODEC_MEM_ERROR;

    vp9_zero(*alg_priv);

    ctx->priv = (vpx_codec_priv_t *)alg_priv;
    ctx->priv->sz = sizeof(*ctx->priv);
    ctx->priv->iface = ctx->iface;
    ctx->priv->alg_priv = alg_priv;
    ctx->priv->alg_priv->si.sz = sizeof(ctx->priv->alg_priv->si);
    ctx->priv->init_flags = ctx->init_flags;
    ctx->priv->alg_priv->frame_parallel_decode =
        (ctx->init_flags & VPX_CODEC_USE_FRAME_THREADING);

    // Disable frame parallel decoding for now.
    ctx->priv->alg_priv->frame_parallel_decode = 0;

    if (ctx->config.dec) {
      // Update the reference to the config structure to an internal copy.
      ctx->priv->alg_priv->cfg = *ctx->config.dec;
      ctx->config.dec = &ctx->priv->alg_priv->cfg;
    }
  }

  return VPX_CODEC_OK;
}

static vpx_codec_err_t decoder_destroy(vpx_codec_alg_priv_t *ctx) {
  if (ctx->frame_workers != NULL) {
    int i;
    for (i = 0; i < ctx->num_frame_workers; ++i) {
      VP9Worker *const worker = &ctx->frame_workers[i];
      FrameWorkerData *const worker_data = (FrameWorkerData *)worker->data1;
      vp9_decoder_remove(worker_data->pbi);
      vpx_free(worker_data);
    }
  }

  vpx_free(ctx->frame_workers);
  vpx_free(ctx);

  return VPX_CODEC_OK;
}

static vpx_codec_err_t decoder_peek_si_internal(const uint8_t *data,
                                                unsigned int data_sz,
                                                vpx_codec_stream_info_t *si,
                                                vpx_decrypt_cb decrypt_cb,
                                                void *decrypt_state) {
  uint8_t clear_buffer[9];

  if (data + data_sz <= data)
    return VPX_CODEC_INVALID_PARAM;

  si->is_kf = 0;
  si->w = si->h = 0;

  if (decrypt_cb) {
    data_sz = MIN(sizeof(clear_buffer), data_sz);
    decrypt_cb(decrypt_state, data, clear_buffer, data_sz);
    data = clear_buffer;
  }

  {
    struct vp9_read_bit_buffer rb = { data, data + data_sz, 0, NULL, NULL };
    const int frame_marker = vp9_rb_read_literal(&rb, 2);
    const int version = vp9_rb_read_bit(&rb);
    (void) vp9_rb_read_bit(&rb);  // unused version bit

    if (frame_marker != VP9_FRAME_MARKER)
      return VPX_CODEC_UNSUP_BITSTREAM;

    if (version > 1) return VPX_CODEC_UNSUP_BITSTREAM;

    if (vp9_rb_read_bit(&rb)) {  // show an existing frame
      return VPX_CODEC_OK;
    }

    if (data_sz <= 8)
      return VPX_CODEC_UNSUP_BITSTREAM;

    si->is_kf = !vp9_rb_read_bit(&rb);
    if (si->is_kf) {
      const int sRGB = 7;
      int colorspace;

      rb.bit_offset += 1;  // show frame
      rb.bit_offset += 1;  // error resilient

      if (vp9_rb_read_literal(&rb, 8) != VP9_SYNC_CODE_0 ||
          vp9_rb_read_literal(&rb, 8) != VP9_SYNC_CODE_1 ||
          vp9_rb_read_literal(&rb, 8) != VP9_SYNC_CODE_2) {
        return VPX_CODEC_UNSUP_BITSTREAM;
      }

      colorspace = vp9_rb_read_literal(&rb, 3);
      if (colorspace != sRGB) {
        rb.bit_offset += 1;  // [16,235] (including xvycc) vs [0,255] range
        if (version == 1) {
          rb.bit_offset += 2;  // subsampling x/y
          rb.bit_offset += 1;  // has extra plane
        }
      } else {
        if (version == 1) {
          rb.bit_offset += 1;  // has extra plane
        } else {
          // RGB is only available in version 1
          return VPX_CODEC_UNSUP_BITSTREAM;
        }
      }

      // TODO(jzern): these are available on non-keyframes in intra only mode.
      si->w = vp9_rb_read_literal(&rb, 16) + 1;
      si->h = vp9_rb_read_literal(&rb, 16) + 1;
    }
  }

  return VPX_CODEC_OK;
}

static vpx_codec_err_t decoder_peek_si(const uint8_t *data,
                                       unsigned int data_sz,
                                       vpx_codec_stream_info_t *si) {
  return decoder_peek_si_internal(data, data_sz, si, NULL, NULL);
}

static vpx_codec_err_t decoder_get_si(vpx_codec_alg_priv_t *ctx,
                                      vpx_codec_stream_info_t *si) {
  const size_t sz = (si->sz >= sizeof(vp9_stream_info_t))
                       ? sizeof(vp9_stream_info_t)
                       : sizeof(vpx_codec_stream_info_t);
  memcpy(si, &ctx->si, sz);
  si->sz = (unsigned int)sz;

  return VPX_CODEC_OK;
}

static void set_error_detail(vpx_codec_alg_priv_t *ctx,
                             const char *const error) {
  ctx->base.err_detail = error;
}

static vpx_codec_err_t update_error_state(vpx_codec_alg_priv_t *ctx,
                           const struct vpx_internal_error_info *error) {
  if (error->error_code)
    set_error_detail(ctx, error->has_detail ? error->detail : NULL);

  return error->error_code;
}

static void init_buffer_callbacks(vpx_codec_alg_priv_t *ctx) {
  int i;

  for (i = 0; i < ctx->num_frame_workers; ++i) {
    VP9Worker *const worker = &ctx->frame_workers[i];
    FrameWorkerData *const worker_data = (FrameWorkerData *)worker->data1;
    VP9_COMMON *const cm = &worker_data->pbi->common;

    cm->new_fb_idx = -1;
    if (ctx->get_ext_fb_cb != NULL && ctx->release_ext_fb_cb != NULL) {
      cm->get_fb_cb = ctx->get_ext_fb_cb;
      cm->release_fb_cb = ctx->release_ext_fb_cb;
      cm->cb_priv = ctx->ext_priv;
    } else {
      cm->get_fb_cb = vp9_get_frame_buffer;
      cm->release_fb_cb = vp9_release_frame_buffer;

      if (vp9_alloc_internal_frame_buffers(&cm->int_frame_buffers))
        vpx_internal_error(&cm->error, VPX_CODEC_MEM_ERROR,
                           "Failed to initialize internal frame buffers");

      cm->cb_priv = &cm->int_frame_buffers;
    }
  }
}

static void set_default_ppflags(vp8_postproc_cfg_t *cfg) {
  cfg->post_proc_flag = VP8_DEBLOCK | VP8_DEMACROBLOCK;
  cfg->deblocking_level = 4;
  cfg->noise_level = 0;
}

static void set_ppflags(const vpx_codec_alg_priv_t *ctx,
                        vp9_ppflags_t *flags) {
  flags->post_proc_flag =
      ctx->postproc_cfg.post_proc_flag;

  flags->deblocking_level = ctx->postproc_cfg.deblocking_level;
  flags->noise_level = ctx->postproc_cfg.noise_level;
}

static int frame_worker_hook(void *arg1, void *arg2) {
  FrameWorkerData *const worker_data = (FrameWorkerData *)arg1;
  const uint8_t *data = worker_data->data;
  (void)arg2;
  worker_data->result = vp9_receive_compressed_data(worker_data->pbi,
                                                    worker_data->data_size,
                                                    &data);
  worker_data->data_end = data;
  return !worker_data->result;
}

static vpx_codec_err_t init_decoder(vpx_codec_alg_priv_t *ctx) {
  int i;

  ctx->last_show_frame = -1;
  ctx->next_submit_thread_id = 0;
  ctx->next_output_thread_id = 0;
  ctx->num_frame_workers =
      (ctx->frame_parallel_decode == 1) ? ctx->cfg.threads: 1;

  ctx->frame_workers = (VP9Worker *)
      vpx_malloc(ctx->num_frame_workers * sizeof(*ctx->frame_workers));
  if (ctx->frame_workers == NULL) {
    set_error_detail(ctx, "Failed to allocate frame_workers");
    return VPX_CODEC_MEM_ERROR;
  }

  for (i = 0; i < ctx->num_frame_workers; ++i) {
    VP9Worker *const worker = &ctx->frame_workers[i];
    FrameWorkerData *worker_data = NULL;
    vp9_worker_init(worker);
    worker->data1 = vpx_memalign(32, sizeof(FrameWorkerData));
    if (worker->data1 == NULL) {
      set_error_detail(ctx, "Failed to allocate worker_data");
      return VPX_CODEC_MEM_ERROR;
    }
    worker_data = (FrameWorkerData *)worker->data1;
    worker_data->pbi = vp9_decoder_create();
    if (worker_data->pbi == NULL) {
      set_error_detail(ctx, "Failed to allocate worker_data");
      return VPX_CODEC_MEM_ERROR;
    }

    // If decoding in serial mode, FrameWorker thread could create tile worker
    // thread or loopfilter thread.
    worker_data->pbi->max_threads =
        (ctx->frame_parallel_decode == 0) ? ctx->cfg.threads : 0;

    worker_data->pbi->inv_tile_order = ctx->invert_tile_order;
    worker_data->pbi->frame_parallel_decode = ctx->frame_parallel_decode;
    worker->hook = (VP9WorkerHook)frame_worker_hook;
  }

  // If postprocessing was enabled by the application and a
  // configuration has not been provided, default it.
  if (!ctx->postproc_cfg_set &&
      (ctx->base.init_flags & VPX_CODEC_USE_POSTPROC))
    set_default_ppflags(&ctx->postproc_cfg);

  init_buffer_callbacks(ctx);

  return VPX_CODEC_OK;
}

static vpx_codec_err_t decode_one(vpx_codec_alg_priv_t *ctx,
                                  const uint8_t **data, unsigned int data_sz,
                                  void *user_priv, int64_t deadline) {
  vp9_ppflags_t flags = {0};
  (void)deadline;

  // Determine the stream parameters. Note that we rely on peek_si to
  // validate that we have a buffer that does not wrap around the top
  // of the heap.
  if (!ctx->si.h) {
    const vpx_codec_err_t res =
        decoder_peek_si_internal(*data, data_sz, &ctx->si, ctx->decrypt_cb,
                                 ctx->decrypt_state);
    if (res != VPX_CODEC_OK)
      return res;

    if (!ctx->si.is_kf)
      return VPX_CODEC_ERROR;
  }

  // Initialize the decoder workers on the first frame
  if (ctx->frame_workers == NULL) {
    const vpx_codec_err_t res = init_decoder(ctx);
    if (res != VPX_CODEC_OK)
      return res;
  }

  if (!ctx->frame_parallel_decode) {
    VP9Worker *const worker = ctx->frame_workers;
    FrameWorkerData *const worker_data = (FrameWorkerData *)worker->data1;
    worker_data->data = *data;
    worker_data->data_size = data_sz;
    worker_data->user_priv = user_priv;

    // Set these even if already initialized.  The caller may have changed the
    // decrypt config between frames.
    worker_data->pbi->decrypt_cb = ctx->decrypt_cb;
    worker_data->pbi->decrypt_state = ctx->decrypt_state;

    vp9_worker_execute(worker);
    if (worker->had_error)
      return update_error_state(ctx, &worker_data->pbi->common.error);

    // Update data pointer after decode.
    *data = worker_data->data_end;
  } else {
    // TODO(hkuang): Implement frame parallel decode.
    return VPX_CODEC_INCAPABLE;
  }

  if (ctx->base.init_flags & VPX_CODEC_USE_POSTPROC)
    set_ppflags(ctx, &flags);

  return VPX_CODEC_OK;
}

static INLINE uint8_t read_marker(vpx_decrypt_cb decrypt_cb,
                                  void *decrypt_state,
                                  const uint8_t *data) {
  if (decrypt_cb) {
    uint8_t marker;
    decrypt_cb(decrypt_state, data, &marker, 1);
    return marker;
  }
  return *data;
}

static vpx_codec_err_t parse_superframe_index(const uint8_t *data,
                                              size_t data_sz,
                                              uint32_t sizes[8], int *count,
                                              vpx_decrypt_cb decrypt_cb,
                                              void *decrypt_state) {
  // A chunk ending with a byte matching 0xc0 is an invalid chunk unless
  // it is a super frame index. If the last byte of real video compression
  // data is 0xc0 the encoder must add a 0 byte. If we have the marker but
  // not the associated matching marker byte at the front of the index we have
  // an invalid bitstream and need to return an error.

  uint8_t marker;

  assert(data_sz);
  marker = read_marker(decrypt_cb, decrypt_state, data + data_sz - 1);
  *count = 0;

  if ((marker & 0xe0) == 0xc0) {
    const uint32_t frames = (marker & 0x7) + 1;
    const uint32_t mag = ((marker >> 3) & 0x3) + 1;
    const size_t index_sz = 2 + mag * frames;

    // This chunk is marked as having a superframe index but doesn't have
    // enough data for it, thus it's an invalid superframe index.
    if (data_sz < index_sz)
      return VPX_CODEC_CORRUPT_FRAME;

    {
      const uint8_t marker2 = read_marker(decrypt_cb, decrypt_state,
                                          data + data_sz - index_sz);

      // This chunk is marked as having a superframe index but doesn't have
      // the matching marker byte at the front of the index therefore it's an
      // invalid chunk.
      if (marker != marker2)
        return VPX_CODEC_CORRUPT_FRAME;
    }

    {
      // Found a valid superframe index.
      uint32_t i, j;
      const uint8_t *x = &data[data_sz - index_sz + 1];

      // Frames has a maximum of 8 and mag has a maximum of 4.
      uint8_t clear_buffer[32];
      assert(sizeof(clear_buffer) >= frames * mag);
      if (decrypt_cb) {
        decrypt_cb(decrypt_state, x, clear_buffer, frames * mag);
        x = clear_buffer;
      }

      for (i = 0; i < frames; ++i) {
        uint32_t this_sz = 0;

        for (j = 0; j < mag; ++j)
          this_sz |= (*x++) << (j * 8);
        sizes[i] = this_sz;
      }
      *count = frames;
    }
  }
  return VPX_CODEC_OK;
}

static vpx_codec_err_t decoder_decode(vpx_codec_alg_priv_t *ctx,
                                      const uint8_t *data, unsigned int data_sz,
                                      void *user_priv, long deadline) {
  const uint8_t *data_start = data;
  const uint8_t * const data_end = data + data_sz;
  vpx_codec_err_t res;
  uint32_t frame_sizes[8];
  int frame_count;

  if (data == NULL || data_sz == 0)
    return VPX_CODEC_INVALID_PARAM;

  res = parse_superframe_index(data, data_sz, frame_sizes, &frame_count,
                               ctx->decrypt_cb, ctx->decrypt_state);
  if (res != VPX_CODEC_OK)
    return res;

  if (ctx->frame_parallel_decode) {
    // Decode in frame parallel mode. When decoding in this mode, the frame
    // passed to the decoder must be either a normal frame or a superframe with
    // superframe index so the decoder could get each frame's start position
    // in the superframe.
    if (frame_count > 0) {
      int i;

      for (i = 0; i < frame_count; ++i) {
        const uint8_t *data_start_copy = data_start;
        const uint32_t frame_size = frame_sizes[i];
        vpx_codec_err_t res;
        if (data_start < data
            || frame_size > (uint32_t) (data_end - data_start)) {
          set_error_detail(ctx, "Invalid frame size in index");
          return VPX_CODEC_CORRUPT_FRAME;
        }

        res = decode_one(ctx, &data_start_copy, frame_size, user_priv,
                         deadline);
        if (res != VPX_CODEC_OK)
          return res;

        data_start += frame_size;
      }
    } else {
      res = decode_one(ctx, &data_start, data_sz, user_priv, deadline);
      if (res != VPX_CODEC_OK)
        return res;

      // Extra data detected after the frame.
      if (data_start < data_end - 1) {
        set_error_detail(ctx, "Fail to decode frame in parallel mode");
        return VPX_CODEC_INCAPABLE;
      }
    }
  } else {
    // Decode in serial mode.
    if (frame_count > 0) {
      int i;

      for (i = 0; i < frame_count; ++i) {
        const uint8_t *data_start_copy = data_start;
        const uint32_t frame_size = frame_sizes[i];
        vpx_codec_err_t res;
        if (data_start < data
            || frame_size > (uint32_t) (data_end - data_start)) {
          set_error_detail(ctx, "Invalid frame size in index");
          return VPX_CODEC_CORRUPT_FRAME;
        }

        res = decode_one(ctx, &data_start_copy, frame_size, user_priv,
                         deadline);
        if (res != VPX_CODEC_OK)
          return res;

        data_start += frame_size;
      }
    } else {
      while (data_start < data_end) {
        const uint32_t frame_size = (uint32_t) (data_end - data_start);
        const vpx_codec_err_t res = decode_one(ctx, &data_start, frame_size,
                                               user_priv, deadline);
        if (res != VPX_CODEC_OK)
          return res;

        // Account for suboptimal termination by the encoder.
        while (data_start < data_end) {
          const uint8_t marker = read_marker(ctx->decrypt_cb,
                                             ctx->decrypt_state, data_start);
          if (marker)
            break;
          ++data_start;
        }
      }
    }
  }

  return VPX_CODEC_OK;
}

static vpx_image_t *decoder_get_frame(vpx_codec_alg_priv_t *ctx,
                                      vpx_codec_iter_t *iter) {
  vpx_image_t *img = NULL;

  // iter acts as a flip flop, so an image is only returned on the first
  // call to get_frame.
  if (*iter == NULL && ctx->frame_workers != NULL) {
    YV12_BUFFER_CONFIG sd;
    vp9_ppflags_t flags = {0, 0, 0};

    VP9Worker *const worker = &ctx->frame_workers[ctx->next_output_thread_id];
    FrameWorkerData *const worker_data = (FrameWorkerData *)worker->data1;
    if (vp9_get_raw_frame(worker_data->pbi, &sd, &flags) == 0) {
      VP9_COMMON *const cm = &worker_data->pbi->common;
      yuvconfig2image(&ctx->img, &sd, worker_data->user_priv);
      ctx->img.fb_priv = cm->frame_bufs[cm->new_fb_idx].raw_frame_buffer.priv;
      img = &ctx->img;
      *iter = img;
      // Decrease reference count of last output frame in frame parallel mode.
      if (ctx->frame_parallel_decode && ctx->last_show_frame >= 0) {
        --cm->frame_bufs[ctx->last_show_frame].ref_count;
        if (cm->frame_bufs[ctx->last_show_frame].ref_count == 0) {
          cm->release_fb_cb(cm->cb_priv,
              &cm->frame_bufs[ctx->last_show_frame].raw_frame_buffer);
        }
      }
      ctx->last_show_frame = worker_data->pbi->common.new_fb_idx;
    }
  }

  return img;
}

static vpx_codec_err_t decoder_set_fb_fn(
    vpx_codec_alg_priv_t *ctx,
    vpx_get_frame_buffer_cb_fn_t cb_get,
    vpx_release_frame_buffer_cb_fn_t cb_release, void *cb_priv) {
  if (cb_get == NULL || cb_release == NULL) {
    return VPX_CODEC_INVALID_PARAM;
  } else if (ctx->frame_workers == NULL) {
    // If the decoder has already been initialized, do not accept changes to
    // the frame buffer functions.
    ctx->get_ext_fb_cb = cb_get;
    ctx->release_ext_fb_cb = cb_release;
    ctx->ext_priv = cb_priv;
    return VPX_CODEC_OK;
  }

  return VPX_CODEC_ERROR;
}

static vpx_codec_err_t ctrl_set_reference(vpx_codec_alg_priv_t *ctx,
                                          va_list args) {
  vpx_ref_frame_t *const data = va_arg(args, vpx_ref_frame_t *);

  // Only support this function in serial decode.
  if (ctx->frame_parallel_decode) {
    set_error_detail(ctx, "Not supported in frame parallel decode");
    return VPX_CODEC_INCAPABLE;
  }

  if (data) {
    vpx_ref_frame_t *const frame = (vpx_ref_frame_t *)data;
    YV12_BUFFER_CONFIG sd;
    VP9Worker *const worker = ctx->frame_workers;
    FrameWorkerData *const worker_data = (FrameWorkerData *)worker->data1;
    image2yuvconfig(&frame->img, &sd);
    return vp9_set_reference_dec(&worker_data->pbi->common,
                                 (VP9_REFFRAME)frame->frame_type, &sd);
  } else {
    return VPX_CODEC_INVALID_PARAM;
  }
}

static vpx_codec_err_t ctrl_copy_reference(vpx_codec_alg_priv_t *ctx,
                                           va_list args) {
  vpx_ref_frame_t *data = va_arg(args, vpx_ref_frame_t *);

  // Only support this function in serial decode.
  if (ctx->frame_parallel_decode) {
    set_error_detail(ctx, "Not supported in frame parallel decode");
    return VPX_CODEC_INCAPABLE;
  }

  if (data) {
    vpx_ref_frame_t *frame = (vpx_ref_frame_t *) data;
    YV12_BUFFER_CONFIG sd;
    VP9Worker *const worker = ctx->frame_workers;
    FrameWorkerData *const worker_data = (FrameWorkerData *)worker->data1;
    image2yuvconfig(&frame->img, &sd);
    return vp9_copy_reference_dec(worker_data->pbi,
                                  (VP9_REFFRAME)frame->frame_type, &sd);
  } else {
    return VPX_CODEC_INVALID_PARAM;
  }
}

static vpx_codec_err_t ctrl_get_reference(vpx_codec_alg_priv_t *ctx,
                                          va_list args) {
  vp9_ref_frame_t *data = va_arg(args, vp9_ref_frame_t *);

  // Only support this function in serial decode.
  if (ctx->frame_parallel_decode) {
    set_error_detail(ctx, "Not supported in frame parallel decode");
    return VPX_CODEC_INCAPABLE;
  }

  if (data) {
    YV12_BUFFER_CONFIG* fb;
    VP9Worker *const worker = ctx->frame_workers;
    FrameWorkerData *const worker_data = (FrameWorkerData *)worker->data1;
    vp9_get_reference_dec(worker_data->pbi, data->idx, &fb);
    yuvconfig2image(&data->img, fb, worker_data->user_priv);
    return VPX_CODEC_OK;
  } else {
    return VPX_CODEC_INVALID_PARAM;
  }
}

static vpx_codec_err_t ctrl_set_postproc(vpx_codec_alg_priv_t *ctx,
                                         va_list args) {
#if CONFIG_VP9_POSTPROC
  vp8_postproc_cfg_t *data = va_arg(args, vp8_postproc_cfg_t *);

  if (data) {
    ctx->postproc_cfg_set = 1;
    ctx->postproc_cfg = *((vp8_postproc_cfg_t *)data);
    return VPX_CODEC_OK;
  } else {
    return VPX_CODEC_INVALID_PARAM;
  }
#else
  (void)ctx;
  (void)args;
  return VPX_CODEC_INCAPABLE;
#endif
}

static vpx_codec_err_t ctrl_set_dbg_options(vpx_codec_alg_priv_t *ctx,
                                            va_list args) {
  (void)ctx;
  (void)args;
  return VPX_CODEC_INCAPABLE;
}

static vpx_codec_err_t ctrl_get_last_ref_updates(vpx_codec_alg_priv_t *ctx,
                                                 va_list args) {
  int *const update_info = va_arg(args, int *);

  // Only support this function in serial decode.
  if (ctx->frame_parallel_decode) {
    set_error_detail(ctx, "Not supported in frame parallel decode");
    return VPX_CODEC_INCAPABLE;
  }

  if (update_info) {
    if (ctx->frame_workers) {
      VP9Worker *const worker = ctx->frame_workers;
      FrameWorkerData *const worker_data = (FrameWorkerData *)worker->data1;
      *update_info = worker_data->pbi->refresh_frame_flags;
    } else {
      return VPX_CODEC_ERROR;
    }
    return VPX_CODEC_OK;
  } else {
    return VPX_CODEC_INVALID_PARAM;
  }
}


static vpx_codec_err_t ctrl_get_frame_corrupted(vpx_codec_alg_priv_t *ctx,
                                                va_list args) {
  int *corrupted = va_arg(args, int *);

  // Only support this function in serial decode.
  if (ctx->frame_parallel_decode) {
    set_error_detail(ctx, "Not supported in frame parallel decode");
    return VPX_CODEC_INCAPABLE;
  }

  if (corrupted) {
    if (ctx->frame_workers) {
      VP9Worker *const worker = ctx->frame_workers;
      FrameWorkerData *const worker_data = (FrameWorkerData *)worker->data1;
      *corrupted = worker_data->pbi->common.frame_to_show->corrupted;
    } else {
      return VPX_CODEC_ERROR;
    }
    return VPX_CODEC_OK;
  } else {
    return VPX_CODEC_INVALID_PARAM;
  }
}

static vpx_codec_err_t ctrl_get_display_size(vpx_codec_alg_priv_t *ctx,
                                             va_list args) {
  int *const display_size = va_arg(args, int *);

  // Only support this function in serial decode.
  if (ctx->frame_parallel_decode) {
    set_error_detail(ctx, "Not supported in frame parallel decode");
    return VPX_CODEC_INCAPABLE;
  }

  if (display_size) {
    if (ctx->frame_workers) {
      VP9Worker *const worker = ctx->frame_workers;
      FrameWorkerData *const worker_data = (FrameWorkerData *)worker->data1;
      const VP9_COMMON *const cm = &worker_data->pbi->common;
      display_size[0] = cm->display_width;
      display_size[1] = cm->display_height;
    } else {
      return VPX_CODEC_ERROR;
    }
    return VPX_CODEC_OK;
  } else {
    return VPX_CODEC_INVALID_PARAM;
  }
}

static vpx_codec_err_t ctrl_set_invert_tile_order(vpx_codec_alg_priv_t *ctx,
                                                  va_list args) {
  ctx->invert_tile_order = va_arg(args, int);
  return VPX_CODEC_OK;
}

static vpx_codec_err_t ctrl_set_decryptor(vpx_codec_alg_priv_t *ctx,
                                          va_list args) {
  vpx_decrypt_init *init = va_arg(args, vpx_decrypt_init *);
  ctx->decrypt_cb = init ? init->decrypt_cb : NULL;
  ctx->decrypt_state = init ? init->decrypt_state : NULL;
  return VPX_CODEC_OK;
}

static vpx_codec_ctrl_fn_map_t decoder_ctrl_maps[] = {
  {VP8_COPY_REFERENCE,            ctrl_copy_reference},

  // Setters
  {VP8_SET_REFERENCE,             ctrl_set_reference},
  {VP8_SET_POSTPROC,              ctrl_set_postproc},
  {VP8_SET_DBG_COLOR_REF_FRAME,   ctrl_set_dbg_options},
  {VP8_SET_DBG_COLOR_MB_MODES,    ctrl_set_dbg_options},
  {VP8_SET_DBG_COLOR_B_MODES,     ctrl_set_dbg_options},
  {VP8_SET_DBG_DISPLAY_MV,        ctrl_set_dbg_options},
  {VP9_INVERT_TILE_DECODE_ORDER,  ctrl_set_invert_tile_order},
  {VPXD_SET_DECRYPTOR,            ctrl_set_decryptor},

  // Getters
  {VP8D_GET_LAST_REF_UPDATES,     ctrl_get_last_ref_updates},
  {VP8D_GET_FRAME_CORRUPTED,      ctrl_get_frame_corrupted},
  {VP9_GET_REFERENCE,             ctrl_get_reference},
  {VP9D_GET_DISPLAY_SIZE,         ctrl_get_display_size},

  { -1, NULL},
};

#ifndef VERSION_STRING
#define VERSION_STRING
#endif
CODEC_INTERFACE(vpx_codec_vp9_dx) = {
  "WebM Project VP9 Decoder" VERSION_STRING,
  VPX_CODEC_INTERNAL_ABI_VERSION,
  VPX_CODEC_CAP_DECODER | VP9_CAP_POSTPROC |
      VPX_CODEC_CAP_EXTERNAL_FRAME_BUFFER,  // vpx_codec_caps_t
  decoder_init,       // vpx_codec_init_fn_t
  decoder_destroy,    // vpx_codec_destroy_fn_t
  decoder_ctrl_maps,  // vpx_codec_ctrl_fn_map_t
  NOT_IMPLEMENTED,    // vpx_codec_get_mmap_fn_t
  NOT_IMPLEMENTED,    // vpx_codec_set_mmap_fn_t
  { // NOLINT
    decoder_peek_si,    // vpx_codec_peek_si_fn_t
    decoder_get_si,     // vpx_codec_get_si_fn_t
    decoder_decode,     // vpx_codec_decode_fn_t
    decoder_get_frame,  // vpx_codec_frame_get_fn_t
    decoder_set_fb_fn,  // vpx_codec_set_fb_fn_t
  },
  { // NOLINT
    NOT_IMPLEMENTED,  // vpx_codec_enc_cfg_map_t
    NOT_IMPLEMENTED,  // vpx_codec_encode_fn_t
    NOT_IMPLEMENTED,  // vpx_codec_get_cx_data_fn_t
    NOT_IMPLEMENTED,  // vpx_codec_enc_config_set_fn_t
    NOT_IMPLEMENTED,  // vpx_codec_get_global_headers_fn_t
    NOT_IMPLEMENTED,  // vpx_codec_get_preview_frame_fn_t
    NOT_IMPLEMENTED   // vpx_codec_enc_mr_get_mem_loc_fn_t
  }
};
