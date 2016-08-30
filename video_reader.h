/*
 *  Copyright (c) 2014 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#ifndef VIDEO_READER_H_
#define VIDEO_READER_H_

#include "./video_common.h"

// The following code is work in progress. It is going to  support transparent
// reading of input files. Right now only IVF format is supported for
// simplicity. The main goal the API is to be simple and easy to use in example
// code and in aomenc/aomdec later. All low-level details like memory
// buffer management are hidden from API users.
struct AvxVideoReaderStruct;
typedef struct AvxVideoReaderStruct AvxVideoReader;

#ifdef __cplusplus
extern "C" {
#endif

// Opens the input file for reading and inspects it to determine file type.
// Returns an opaque AvxVideoReader* upon success, or NULL upon failure.
// Right now only IVF format is supported.
AvxVideoReader *aom_video_reader_open(const char *filename);

// Frees all resources associated with AvxVideoReader* returned from
// aom_video_reader_open() call.
void aom_video_reader_close(AvxVideoReader *reader);

// Reads frame from the file and stores it in internal buffer.
int aom_video_reader_read_frame(AvxVideoReader *reader);

// Returns the pointer to memory buffer with frame data read by last call to
// aom_video_reader_read_frame().
const uint8_t *aom_video_reader_get_frame(AvxVideoReader *reader, size_t *size);

// Fills AvxVideoInfo with information from opened video file.
const AvxVideoInfo *aom_video_reader_get_info(AvxVideoReader *reader);

#ifdef __cplusplus
}  // extern "C"
#endif

#endif  // VIDEO_READER_H_
