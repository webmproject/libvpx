/*
 *  Copyright (c) 2013 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include "./ivfdec.h"

#include <stdio.h>
#include <stdlib.h>

int file_is_ivf(struct VpxInputContext *input_ctx) {
  char raw_hdr[32];
  int is_ivf = 0;

  // TODO(tomfinegan): This can eventually go away, but for now it's required
  // because the means by which file types are detected differ in vpxdec and
  // vpxenc.
  rewind(input_ctx->file);

  if (fread(raw_hdr, 1, 32, input_ctx->file) == 32) {
    if (raw_hdr[0] == 'D' && raw_hdr[1] == 'K' &&
        raw_hdr[2] == 'I' && raw_hdr[3] == 'F') {
      is_ivf = 1;

      if (mem_get_le16(raw_hdr + 4) != 0) {
        fprintf(stderr, "Error: Unrecognized IVF version! This file may not"
                " decode properly.");
      }

      input_ctx->fourcc = mem_get_le32(raw_hdr + 8);
      input_ctx->width = mem_get_le16(raw_hdr + 12);
      input_ctx->height = mem_get_le16(raw_hdr + 14);
      input_ctx->framerate.numerator = mem_get_le32(raw_hdr + 16);
      input_ctx->framerate.denominator = mem_get_le32(raw_hdr + 20);

      /* Some versions of vpxenc used 1/(2*fps) for the timebase, so
       * we can guess the framerate using only the timebase in this
       * case. Other files would require reading ahead to guess the
       * timebase, like we do for webm.
       */
      if (input_ctx->framerate.numerator < 1000) {
        /* Correct for the factor of 2 applied to the timebase in the
         * encoder.
         */
        if (input_ctx->framerate.numerator & 1)
          input_ctx->framerate.denominator <<= 1;
        else
          input_ctx->framerate.numerator >>= 1;
      } else {
        /* Don't know FPS for sure, and don't have readahead code
         * (yet?), so just default to 30fps.
         */
        input_ctx->framerate.numerator = 30;
        input_ctx->framerate.denominator = 1;
      }
    }
  }

  if (!is_ivf) {
    rewind(input_ctx->file);
    input_ctx->detect.buf_read = 0;
  } else {
    input_ctx->detect.position = 4;
  }
  return is_ivf;
}

int ivf_read_frame(struct VpxInputContext *input_ctx,
                   uint8_t **buffer,
                   size_t *bytes_read,
                   size_t *buffer_size) {
  char raw_header[IVF_FRAME_HDR_SZ] = {0};
  size_t frame_size = 0;
  FILE *infile = input_ctx->file;

  if (input_ctx->file_type != FILE_TYPE_IVF)
    return 0;

  if (fread(raw_header, IVF_FRAME_HDR_SZ, 1, infile) != 1) {
    if (!feof(infile))
      warn("Failed to read frame size\n");
  } else {
    frame_size = mem_get_le32(raw_header);

    if (frame_size > 256 * 1024 * 1024) {
      warn("Read invalid frame size (%u)\n", (unsigned int)frame_size);
      frame_size = 0;
    }

    if (frame_size > *buffer_size) {
      uint8_t *new_buffer = realloc(*buffer, 2 * frame_size);

      if (new_buffer) {
        *buffer = new_buffer;
        *buffer_size = 2 * frame_size;
      } else {
        warn("Failed to allocate compressed data buffer\n");
        frame_size = 0;
      }
    }
  }

  if (!feof(infile)) {
    if (fread(*buffer, 1, frame_size, infile) != frame_size) {
      warn("Failed to read full frame\n");
      return 1;
    }

    *bytes_read = frame_size;
    return 0;
  }

  return 1;
}
