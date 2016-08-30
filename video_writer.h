/*
 *  Copyright (c) 2014 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#ifndef VIDEO_WRITER_H_
#define VIDEO_WRITER_H_

#include "./video_common.h"

typedef enum { kContainerIVF } AvxContainer;

struct AvxVideoWriterStruct;
typedef struct AvxVideoWriterStruct AvxVideoWriter;

#ifdef __cplusplus
extern "C" {
#endif

// Finds and opens writer for specified container format.
// Returns an opaque AvxVideoWriter* upon success, or NULL upon failure.
// Right now only IVF format is supported.
AvxVideoWriter *aom_video_writer_open(const char *filename,
                                      AvxContainer container,
                                      const AvxVideoInfo *info);

// Frees all resources associated with AvxVideoWriter* returned from
// aom_video_writer_open() call.
void aom_video_writer_close(AvxVideoWriter *writer);

// Writes frame bytes to the file.
int aom_video_writer_write_frame(AvxVideoWriter *writer, const uint8_t *buffer,
                                 size_t size, int64_t pts);

#ifdef __cplusplus
}  // extern "C"
#endif

#endif  // VIDEO_WRITER_H_
