/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include <assert.h>
#include <limits.h>
#include <stdio.h>

#include "./vpx_scale_rtcd.h"

#include "vpx_mem/vpx_mem.h"
#include "vpx_ports/vpx_timer.h"
#include "vpx_scale/vpx_scale.h"

#include "vp9/common/vp9_alloccommon.h"
#include "vp9/common/vp9_loopfilter.h"
#include "vp9/common/vp9_onyxc_int.h"
#if CONFIG_VP9_POSTPROC
#include "vp9/common/vp9_postproc.h"
#endif
#include "vp9/common/vp9_quant_common.h"
#include "vp9/common/vp9_reconinter.h"
#include "vp9/common/vp9_reconintra.h"
#include "vp9/common/vp9_systemdependent.h"

#include "vp9/decoder/vp9_decodeframe.h"
#include "vp9/decoder/vp9_decoder.h"
#include "vp9/decoder/vp9_detokenize.h"
#include "vp9/decoder/vp9_dthread.h"

static void initialize_dec() {
  static int init_done = 0;

  if (!init_done) {
    vp9_rtcd();
    vp9_init_intra_predictors();
#if CONFIG_WEDGE_PARTITION
    vp9_init_wedge_masks();
#endif  // CONFIG_WEDGE_PARTITION
    init_done = 1;
  }
}

VP9Decoder *vp9_decoder_create() {
  VP9Decoder *const pbi = vpx_memalign(32, sizeof(*pbi));
  VP9_COMMON *const cm = pbi ? &pbi->common : NULL;

  if (!cm)
    return NULL;

  vp9_zero(*pbi);

  if (setjmp(cm->error.jmp)) {
    cm->error.setjmp = 0;
    vp9_decoder_remove(pbi);
    return NULL;
  }

  cm->error.setjmp = 1;
  pbi->need_resync = 1;
  initialize_dec();

  // Initialize the references to not point to any frame buffers.
  vpx_memset(&cm->ref_frame_map, -1, sizeof(cm->ref_frame_map));

  cm->current_video_frame = 0;
  pbi->ready_for_new_data = 1;
  cm->bit_depth = VPX_BITS_8;
  cm->dequant_bit_depth = VPX_BITS_8;

  // vp9_init_dequantizer() is first called here. Add check in
  // frame_init_dequantizer() to avoid unnecessary calling of
  // vp9_init_dequantizer() for every frame.
  vp9_init_dequantizer(cm);

  vp9_loop_filter_init(cm);

  cm->error.setjmp = 0;

  vp9_get_worker_interface()->init(&pbi->lf_worker);

  return pbi;
}

void vp9_decoder_remove(VP9Decoder *pbi) {
  VP9_COMMON *const cm = &pbi->common;
  int i;

  vp9_get_worker_interface()->end(&pbi->lf_worker);
  vpx_free(pbi->lf_worker.data1);
  vpx_free(pbi->tile_data);
  for (i = 0; i < pbi->num_tile_workers; ++i) {
    VP9Worker *const worker = &pbi->tile_workers[i];
    vp9_get_worker_interface()->end(worker);
    vpx_free(worker->data1);
    vpx_free(worker->data2);
  }
  vpx_free(pbi->tile_workers);

  if (pbi->num_tile_workers > 0) {
    vp9_loop_filter_dealloc(&pbi->lf_row_sync);
  }

  vp9_remove_common(cm);
  vpx_free(pbi);
}

static int equal_dimensions(const YV12_BUFFER_CONFIG *a,
                            const YV12_BUFFER_CONFIG *b) {
    return a->y_height == b->y_height && a->y_width == b->y_width &&
           a->uv_height == b->uv_height && a->uv_width == b->uv_width;
}

vpx_codec_err_t vp9_copy_reference_dec(VP9Decoder *pbi,
                                       VP9_REFFRAME ref_frame_flag,
                                       YV12_BUFFER_CONFIG *sd) {
  VP9_COMMON *cm = &pbi->common;

  /* TODO(jkoleszar): The decoder doesn't have any real knowledge of what the
   * encoder is using the frame buffers for. This is just a stub to keep the
   * vpxenc --test-decode functionality working, and will be replaced in a
   * later commit that adds VP9-specific controls for this functionality.
   */
  if (ref_frame_flag == VP9_LAST_FLAG) {
    const YV12_BUFFER_CONFIG *const cfg = get_ref_frame(cm, 0);
    if (cfg == NULL) {
      vpx_internal_error(&cm->error, VPX_CODEC_ERROR,
                         "No 'last' reference frame");
      return VPX_CODEC_ERROR;
    }
    if (!equal_dimensions(cfg, sd))
      vpx_internal_error(&cm->error, VPX_CODEC_ERROR,
                         "Incorrect buffer dimensions");
    else
      vp8_yv12_copy_frame(cfg, sd);
  } else {
    vpx_internal_error(&cm->error, VPX_CODEC_ERROR,
                       "Invalid reference frame");
  }

  return cm->error.error_code;
}


vpx_codec_err_t vp9_set_reference_dec(VP9_COMMON *cm,
                                      VP9_REFFRAME ref_frame_flag,
                                      YV12_BUFFER_CONFIG *sd) {
  int idx;
  YV12_BUFFER_CONFIG *ref_buf = NULL;

  // TODO(jkoleszar): The decoder doesn't have any real knowledge of what the
  // encoder is using the frame buffers for. This is just a stub to keep the
  // vpxenc --test-decode functionality working, and will be replaced in a
  // later commit that adds VP9-specific controls for this functionality.

  // (Yunqing) The set_reference control depends on the following setting in
  // encoder.
  //    cpi->lst_fb_idx = 0;
  //   #if CONFIG_MULTI_REF
  //     cpi->lst2_fb_idx = 1;
  //     cpi->lst3_fb_idx = 2;
  //     cpi->lst4_fb_idx = 3;
  //     cpi->gld_fb_idx = 4;
  //     cpi->alt_fb_idx = 5;
  //   #else  // CONFIG_MULTI_REF
  //     cpi->gld_fb_idx = 1;
  //     cpi->alt_fb_idx = 2;

  if (ref_frame_flag == VP9_LAST_FLAG) {
    idx = cm->ref_frame_map[0];
#if CONFIG_MULTI_REF
  } else if (ref_frame_flag == VP9_LAST2_FLAG) {
    idx = cm->ref_frame_map[1];
  } else if (ref_frame_flag == VP9_LAST3_FLAG) {
    idx = cm->ref_frame_map[2];
  } else if (ref_frame_flag == VP9_LAST4_FLAG) {
    idx = cm->ref_frame_map[3];
  } else if (ref_frame_flag == VP9_GOLD_FLAG) {
    idx = cm->ref_frame_map[4];
  } else if (ref_frame_flag == VP9_ALT_FLAG) {
    idx = cm->ref_frame_map[5];
#else  // CONFIG_MULTI_REF
  } else if (ref_frame_flag == VP9_GOLD_FLAG) {
    idx = cm->ref_frame_map[1];
  } else if (ref_frame_flag == VP9_ALT_FLAG) {
    idx = cm->ref_frame_map[2];
#endif  // CONFIG_MULTI_REF
  } else {
    vpx_internal_error(&cm->error, VPX_CODEC_ERROR,
                       "Invalid reference frame");
    return cm->error.error_code;
  }

  if (idx < 0 || idx >= FRAME_BUFFERS) {
    vpx_internal_error(&cm->error, VPX_CODEC_ERROR,
                       "Invalid reference frame map");
    return cm->error.error_code;
  }

  // Get the destination reference buffer.
  ref_buf = &cm->frame_bufs[idx].buf;

  if (!equal_dimensions(ref_buf, sd)) {
    vpx_internal_error(&cm->error, VPX_CODEC_ERROR,
                       "Incorrect buffer dimensions");
  } else {
    // Overwrite the reference frame buffer.
    vp8_yv12_copy_frame(sd, ref_buf);
  }

  return cm->error.error_code;
}

/* If any buffer updating is signaled it should be done here. */
static void swap_frame_buffers(VP9Decoder *pbi) {
  int ref_index = 0, mask;
  VP9_COMMON *const cm = &pbi->common;

  /*
  {
    int i;
    printf("\n==============================================================\n");
    printf("*******swap_frame_buffers(frame=%d): refresh_frame_flags=%d, "
           "new_fb_idx=%d, ref_frame_map=[",
           cm->current_video_frame, pbi->refresh_frame_flags, cm->new_fb_idx);
    for (i = 0; i < REF_FRAMES; ++i)
      printf(" %d", cm->ref_frame_map[i]);
    printf(" ]\n");
  }
  */
  for (mask = pbi->refresh_frame_flags; mask; mask >>= 1) {
#if CONFIG_MULTI_REF
    // NOTE(zoeliu):
    // (1) When ref_index == 1 and the mask bit is set, it indicates LAST2_FRAME
    //     shall be refreshed, but it has already been handled when LAST_FRAME
    //     is being refreshed, i.e., when ref_index == 0 and mask is set
    //     correspondingly;
    // (2) The only exception is that when current frame is a KEY_FRAME, where
    //     all the frames in the frame buffer (cm->frame_bufs) get refreshed.
    if ((mask & 1) &&
        (cm->frame_type == KEY_FRAME ||
         (ref_index != 1 && ref_index != 2 && ref_index != 3))
        ) {
#else  // CONFIG_MULTI_REF
    if (mask & 1) {
#endif  // CONFIG_MULTI_REF
      const int old_idx = cm->ref_frame_map[ref_index];

#if CONFIG_MULTI_REF
      if ((ref_index == 0) && (mask & 2) && (cm->frame_type != KEY_FRAME)) {
        // If current is LAST_FRAME, and LAST2_FRAME needs to get refreshed
        // as well, then LAST2_FRAME should get handled first.
        const int ref_index_last2 = ref_index + 1;
        const int old_buf_idx_last2 = cm->ref_frame_map[ref_index_last2];

        if (mask & 4) {
          // If current is LAST_FRAME, and both LAST2_FRAME and LAST3_FRAME need
          // to get refreshed as well, then LAST3_FRAME should get handled first
          // followed by the handling of LAST2_FRAME.
          const int ref_index_last3 = ref_index_last2 + 1;
          const int old_buf_idx_last3 = cm->ref_frame_map[ref_index_last3];

          if (mask & 8) {
            // If current is LAST_FRAME, and LAST2, LAST3 and LAST4 all need to
            // get refreshed as well, then LAST4_FRAME should get handled first
            // followed by the handling of LAST3_FRAME and then LAST2_FRAME.
            const int ref_index_last4 = ref_index_last3 + 1;
            const int old_buf_idx_last4 = cm->ref_frame_map[ref_index_last4];

            ref_cnt_fb(cm->frame_bufs,
                       &cm->ref_frame_map[ref_index_last4],
                       cm->ref_frame_map[ref_index_last3]);
            if (old_buf_idx_last4 >= 0 &&
                cm->frame_bufs[old_buf_idx_last4].ref_count == 0) {
              cm->release_fb_cb(
                  cm->cb_priv,
                  &cm->frame_bufs[old_buf_idx_last4].raw_frame_buffer);
            }
          }

          ref_cnt_fb(cm->frame_bufs,
                     &cm->ref_frame_map[ref_index_last3],
                     cm->ref_frame_map[ref_index_last2]);
          if (old_buf_idx_last3 >= 0 &&
              cm->frame_bufs[old_buf_idx_last3].ref_count == 0) {
            cm->release_fb_cb(
                cm->cb_priv,
                &cm->frame_bufs[old_buf_idx_last3].raw_frame_buffer);
          }
        }

        ref_cnt_fb(cm->frame_bufs,
                   &cm->ref_frame_map[ref_index_last2],
                   cm->ref_frame_map[ref_index]);
        if (old_buf_idx_last2 >= 0 &&
            cm->frame_bufs[old_buf_idx_last2].ref_count == 0) {
          cm->release_fb_cb(
              cm->cb_priv,
              &cm->frame_bufs[old_buf_idx_last2].raw_frame_buffer);
        }
      }
#endif  // CONFIG_MULTI_REF

      ref_cnt_fb(cm->frame_bufs, &cm->ref_frame_map[ref_index],
                 cm->new_fb_idx);
      if (old_idx >= 0 && cm->frame_bufs[old_idx].ref_count == 0)
        cm->release_fb_cb(cm->cb_priv,
                          &cm->frame_bufs[old_idx].raw_frame_buffer);
    }
    ++ref_index;
  }
  /*
  {
    int i;
    printf("*******swap_frame_buffers(frame=%d): refresh_frame_flags=%d, "
           "new_fb_idx=%d, ref_frame_map=[",
           cm->current_video_frame, pbi->refresh_frame_flags, cm->new_fb_idx);
    for (i = 0; i < REF_FRAMES; ++i)
      printf(" %d", cm->ref_frame_map[i]);
    printf(" ]\n");
    printf("==============================================================\n");
  }
  */

  cm->frame_to_show = get_frame_new_buffer(cm);
  cm->frame_bufs[cm->new_fb_idx].ref_count--;

  // Invalidate these references until the next frame starts.
  for (ref_index = 0; ref_index < REFS_PER_FRAME; ref_index++)
    cm->frame_refs[ref_index].idx = INT_MAX;
}

int vp9_receive_compressed_data(VP9Decoder *pbi,
                                size_t size, const uint8_t **psource) {
  VP9_COMMON *const cm = &pbi->common;
  const uint8_t *source = *psource;
  int retcode = 0;

  cm->error.error_code = VPX_CODEC_OK;

  if (size == 0) {
    // This is used to signal that we are missing frames.
    // We do not know if the missing frame(s) was supposed to update
    // any of the reference buffers, but we act conservative and
    // mark only the last buffer as corrupted.
    //
    // TODO(jkoleszar): Error concealment is undefined and non-normative
    // at this point, but if it becomes so, [0] may not always be the correct
    // thing to do here.
    if (cm->frame_refs[0].idx != INT_MAX)
      cm->frame_refs[0].buf->corrupted = 1;
  }

  pbi->ready_for_new_data = 0;

  // Check if the previous frame was a frame without any references to it.
  if (cm->new_fb_idx >= 0 && cm->frame_bufs[cm->new_fb_idx].ref_count == 0)
    cm->release_fb_cb(cm->cb_priv,
                      &cm->frame_bufs[cm->new_fb_idx].raw_frame_buffer);
  cm->new_fb_idx = get_free_fb(cm);

  if (setjmp(cm->error.jmp)) {
    pbi->need_resync = 1;
    cm->error.setjmp = 0;
    vp9_clear_system_state();

    // We do not know if the missing frame(s) was supposed to update
    // any of the reference buffers, but we act conservative and
    // mark only the last buffer as corrupted.
    //
    // TODO(jkoleszar): Error concealment is undefined and non-normative
    // at this point, but if it becomes so, [0] may not always be the correct
    // thing to do here.
    if (cm->frame_refs[0].idx != INT_MAX && cm->frame_refs[0].buf != NULL)
      cm->frame_refs[0].buf->corrupted = 1;

    if (cm->new_fb_idx > 0 && cm->frame_bufs[cm->new_fb_idx].ref_count > 0)
      cm->frame_bufs[cm->new_fb_idx].ref_count--;

    return -1;
  }

  cm->error.setjmp = 1;

#if CONFIG_PALETTE
  if (frame_is_intra_only(cm)) {
    cm->current_palette_size = 0;
    vpx_memset(cm->current_palette_count, 0,
               PALETTE_BUF_SIZE * sizeof(cm->current_palette_count[0]));
  }
#endif

  vp9_decode_frame(pbi, source, source + size, psource);

  swap_frame_buffers(pbi);

  vp9_clear_system_state();

  cm->last_width = cm->width;
  cm->last_height = cm->height;

  if (!cm->show_existing_frame)
    cm->last_show_frame = cm->show_frame;
  if (cm->show_frame) {
    if (!cm->show_existing_frame)
      vp9_swap_mi_and_prev_mi(cm);

    cm->current_video_frame++;
  }

  cm->error.setjmp = 0;
  return retcode;
}

int vp9_get_raw_frame(VP9Decoder *pbi, YV12_BUFFER_CONFIG *sd,
                      vp9_ppflags_t *flags) {
  VP9_COMMON *const cm = &pbi->common;
  int ret = -1;
#if !CONFIG_VP9_POSTPROC
  (void)*flags;
#endif

  if (pbi->ready_for_new_data == 1)
    return ret;

  pbi->ready_for_new_data = 1;

  /* no raw frame to show!!! */
  if (!cm->show_frame)
    return ret;

#if CONFIG_VP9_POSTPROC
  if (!cm->show_existing_frame) {
    ret = vp9_post_proc_frame(cm, sd, flags);
  } else {
    *sd = *cm->frame_to_show;
    ret = 0;
  }
#else
  *sd = *cm->frame_to_show;
  ret = 0;
#endif /*!CONFIG_POSTPROC*/
  vp9_clear_system_state();
  return ret;
}

vpx_codec_err_t vp9_parse_superframe_index(const uint8_t *data,
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
