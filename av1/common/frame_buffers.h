/*
 *  Copyright (c) 2014 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#ifndef AV1_COMMON_FRAME_BUFFERS_H_
#define AV1_COMMON_FRAME_BUFFERS_H_

#include "aom/aom_frame_buffer.h"
#include "aom/aom_integer.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct InternalFrameBuffer {
  uint8_t *data;
  size_t size;
  int in_use;
} InternalFrameBuffer;

typedef struct InternalFrameBufferList {
  int num_internal_frame_buffers;
  InternalFrameBuffer *int_fb;
} InternalFrameBufferList;

// Initializes |list|. Returns 0 on success.
int av1_alloc_internal_frame_buffers(InternalFrameBufferList *list);

// Free any data allocated to the frame buffers.
void av1_free_internal_frame_buffers(InternalFrameBufferList *list);

// Callback used by libaom to request an external frame buffer. |cb_priv|
// Callback private data, which points to an InternalFrameBufferList.
// |min_size| is the minimum size in bytes needed to decode the next frame.
// |fb| pointer to the frame buffer.
int av1_get_frame_buffer(void *cb_priv, size_t min_size,
                         aom_codec_frame_buffer_t *fb);

// Callback used by libaom when there are no references to the frame buffer.
// |cb_priv| is not used. |fb| pointer to the frame buffer.
int av1_release_frame_buffer(void *cb_priv, aom_codec_frame_buffer_t *fb);

#ifdef __cplusplus
}  // extern "C"
#endif

#endif  // AV1_COMMON_FRAME_BUFFERS_H_
