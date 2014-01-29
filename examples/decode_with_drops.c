/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

// Decode With Drops Example
// =========================
//
// This is an example utility which drops a series of frames, as specified
// on the command line. This is useful for observing the error recovery
// features of the codec.
//
// Usage
// -----
// This example adds a single argument to the `simple_decoder` example,
// which specifies the range or pattern of frames to drop. The parameter is
// parsed as follows:
//
// Dropping A Range Of Frames
// --------------------------
// To drop a range of frames, specify the starting frame and the ending
// frame to drop, separated by a dash. The following command will drop
// frames 5 through 10 (base 1).
//
//  $ ./decode_with_drops in.ivf out.i420 5-10
//
//
// Dropping A Pattern Of Frames
// ----------------------------
// To drop a pattern of frames, specify the number of frames to drop and
// the number of frames after which to repeat the pattern, separated by
// a forward-slash. The following command will drop 3 of 7 frames.
// Specifically, it will decode 4 frames, then drop 3 frames, and then
// repeat.
//
//  $ ./decode_with_drops in.ivf out.i420 3/7
//
//
// Extra Variables
// ---------------
// This example maintains the pattern passed on the command line in the
// `n`, `m`, and `is_range` variables:
//
//
// Making The Drop Decision
// ------------------------
// The example decides whether to drop the frame based on the current
// frame number, immediately before decoding the frame.

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "./ivfdec.h"

#define VPX_CODEC_DISABLE_COMPAT 1

#include "vpx/vp8dx.h"
#include "vpx/vpx_decoder.h"

#include "./tools_common.h"
#include "./vpx_config.h"

static const char *exec_name;

void usage_exit() {
  fprintf(stderr, "Usage: %s <infile> <outfile> <N-M|N/M>\n", exec_name);
  exit(EXIT_FAILURE);
}

int main(int argc, char **argv) {
  FILE *infile, *outfile;
  vpx_codec_ctx_t codec;
  vpx_codec_iface_t *iface;
  int flags = 0, frame_cnt = 0;
  vpx_video_t *video;
  int n, m, is_range;
  char *nptr;

  exec_name = argv[0];

  if (argc != 4)
    die("Invalid number of arguments");

  if (!(infile = fopen(argv[1], "rb")))
    die("Failed to open %s for reading", argv[1]);

  if (!(outfile = fopen(argv[2], "wb")))
    die("Failed to open %s for writing", argv[2]);

  n = strtol(argv[3], &nptr, 0);
  m = strtol(nptr + 1, NULL, 0);
  is_range = (*nptr == '-');
  if (!n || !m || (*nptr != '-' && *nptr != '/'))
    die("Couldn't parse pattern %s\n", argv[3]);

  video = vpx_video_open_file(infile);
  if (!video)
    die("%s is not a supported input file.", argv[1]);

  iface = get_codec_interface(vpx_video_get_fourcc(video));
  if (!iface)
    die("Unknown FOURCC code.");

  printf("Using %s\n", vpx_codec_iface_name(iface));

  if (vpx_codec_dec_init(&codec, iface, NULL, flags))
    die_codec(&codec, "Failed to initialize decoder");

  while (vpx_video_read_frame(video)) {
    vpx_codec_iter_t iter = NULL;
    vpx_image_t *img = NULL;
    size_t frame_size = 0;
    int skip;
    const unsigned char *frame = vpx_video_get_frame(video, &frame_size);
    if (vpx_codec_decode(&codec, frame, frame_size, NULL, 0))
      die_codec(&codec, "Failed to decode frame");

    ++frame_cnt;

    skip = (is_range && frame_cnt >= n && frame_cnt <= m) ||
           (!is_range && m - (frame_cnt - 1) % m <= n);

    if (!skip) {
      putc('.', stdout);

      while ((img = vpx_codec_get_frame(&codec, &iter)) != NULL)
        vpx_img_write(img, outfile);
    } else {
      putc('X', stdout);
    }

    fflush(stdout);
  }

  printf("Processed %d frames.\n", frame_cnt);
  if (vpx_codec_destroy(&codec))
    die_codec(&codec, "Failed to destroy codec");

  printf("Play: ffplay -f rawvideo -pix_fmt yuv420p -s %dx%d %s\n",
         vpx_video_get_width(video), vpx_video_get_height(video), argv[2]);

  vpx_video_close(video);

  fclose(outfile);
  fclose(infile);

  return EXIT_SUCCESS;
}
