/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

// Postprocessing Decoder
// ======================
//
// This example adds postprocessing to the simple decoder loop.
//
// Initializing Postprocessing
// ---------------------------
// You must inform the codec that you might request postprocessing at
// initialization time. This is done by passing the VPX_CODEC_USE_POSTPROC
// flag to `vpx_codec_dec_init`. If the codec does not support
// postprocessing, this call will return VPX_CODEC_INCAPABLE. For
// demonstration purposes, we also fall back to default initialization if
// the codec does not provide support.
//
// Using Adaptive Postprocessing
// -----------------------------
// VP6 provides "adaptive postprocessing." It will automatically select the
// best postprocessing filter on a frame by frame basis based on the amount
// of time remaining before the user's specified deadline expires. The
// special value 0 indicates that the codec should take as long as
// necessary to provide the best quality frame. This example gives the
// codec 15ms (15000us) to return a frame. Remember that this is a soft
// deadline, and the codec may exceed it doing its regular processing. In
// these cases, no additional postprocessing will be done.
//
// Codec Specific Postprocessing Controls
// --------------------------------------
// Some codecs provide fine grained controls over their built-in
// postprocessors. VP8 is one example. The following sample code toggles
// postprocessing on and off every 15 frames.

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
  fprintf(stderr, "Usage: %s <infile> <outfile>\n", exec_name);
  exit(EXIT_FAILURE);
}

int main(int argc, char **argv) {
  FILE *infile, *outfile;
  vpx_codec_ctx_t codec;
  vpx_codec_iface_t *iface;
  int frame_cnt = 0;
  vpx_video_t *video;
  vpx_codec_err_t res;

  exec_name = argv[0];

  if (argc != 3)
    die("Invalid number of arguments");

  if (!(infile = fopen(argv[1], "rb")))
    die("Failed to open %s for reading", argv[1]);

  if (!(outfile = fopen(argv[2], "wb")))
    die("Failed to open %s for writing", argv[2]);

  video = vpx_video_open_file(infile);
  if (!video)
    die("%s is not a supported input file.", argv[1]);

  iface = get_codec_interface(vpx_video_get_fourcc(video));
  if (!iface)
    die("Unknown FOURCC code.");

  printf("Using %s\n", vpx_codec_iface_name(iface));


  res = vpx_codec_dec_init(&codec, iface, NULL, VPX_CODEC_USE_POSTPROC);
  if (res == VPX_CODEC_INCAPABLE) {
    printf("NOTICE: Postproc not supported.\n");
    res = vpx_codec_dec_init(&codec, iface, NULL, 0);
  }

  if (res)
    die_codec(&codec, "Failed to initialize decoder");

  while (vpx_video_read_frame(video)) {
    vpx_codec_iter_t iter = NULL;
    vpx_image_t *img = NULL;
    size_t frame_size = 0;
    const unsigned char *frame = vpx_video_get_frame(video, &frame_size);

    ++frame_cnt;

    if (frame_cnt % 30 == 1) {
      vp8_postproc_cfg_t pp = {0, 0, 0};

    if (vpx_codec_control(&codec, VP8_SET_POSTPROC, &pp))
      die_codec(&codec, "Failed to turn off postproc");
    } else if (frame_cnt % 30 == 16) {
      vp8_postproc_cfg_t pp = {VP8_DEBLOCK | VP8_DEMACROBLOCK | VP8_MFQE,
                               4, 0};
      if (vpx_codec_control(&codec, VP8_SET_POSTPROC, &pp))
        die_codec(&codec, "Failed to turn on postproc");
    };

    // Decode the frame with 15ms deadline
    if (vpx_codec_decode(&codec, frame, frame_size, NULL, 15000))
      die_codec(&codec, "Failed to decode frame");

    while ((img = vpx_codec_get_frame(&codec, &iter)) != NULL) {
      vpx_img_write(img, outfile);
    }
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
