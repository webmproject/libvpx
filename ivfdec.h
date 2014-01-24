/*
 *  Copyright (c) 2013 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */
#ifndef IVFDEC_H_
#define IVFDEC_H_

#include "./tools_common.h"

#ifdef __cplusplus
extern "C" {
#endif

int file_is_ivf(struct VpxInputContext *input);

int ivf_read_frame(FILE *infile, uint8_t **buffer,
                   size_t *bytes_read, size_t *buffer_size);

// The following code is work in progress. It is going to be in a separate file
// and support transparent reading of IVF and Y4M formats. Right now only IVF
// format is supported for simplicity. The main goal the API is to be
// simple and easy to use in example code (and probably in vpxenc/vpxdec later).
// All low-level details like memory buffer management are hidden from API
// users.
struct vpx_video;
typedef struct vpx_video vpx_video_t;

// Opens the input file and inspects it to determine file type. Returns an
// opaque vpx_video_t* upon success, or NULL upon failure.
vpx_video_t *vpx_video_open_file(FILE *file);

// Frees all resources associated with vpx_video_t returned from
// vpx_video_open_file() call
void vpx_video_close(vpx_video_t *video);

int vpx_video_get_width(vpx_video_t *video);
int vpx_video_get_height(vpx_video_t *video);
unsigned int vpx_video_get_fourcc(vpx_video_t *video);

// Reads video frame bytes from the file and stores them into internal buffer.
int vpx_video_read_frame(vpx_video_t *video);

// Returns the pointer to internal memory buffer with frame bytes read from
// last call to vpx_video_read_frame().
const unsigned char *vpx_video_get_frame(vpx_video_t *video, size_t *size);

#ifdef __cplusplus
}  /* extern "C" */
#endif

#endif  // IVFDEC_H_
