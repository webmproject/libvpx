/*
 *  Copyright (c) 2012 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

/*
 * This is an example demonstrating how to implement a multi-layer
 * VP9 encoding scheme based on spatial scalability for video applications
 * that benefit from a scalable bitstream.
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <time.h>
#include <string.h>
#include <unistd.h>
#include <libgen.h>
#define VPX_CODEC_DISABLE_COMPAT 1
#include "vpx/vpx_encoder.h"
#include "vpx/vp8cx.h"
#include "vpx/svc_context.h"

#define interface (vpx_codec_vp9_cx())
#define fourcc 0x30395056
#define IVF_FILE_HDR_SZ (32)
#define IVF_FRAME_HDR_SZ (12)

char *input_filename;
char *output_filename;
unsigned int number_frames_to_code = 60 * 60;
unsigned int number_frames_to_skip = 0;
unsigned int gop_size = 100;

char *scaling_factor;
char *quantizer;
SVC_ENCODING_MODE encoding_mode = INTER_LAYER_PREDICTION_IP;

static void mem_put_le16(char *mem, unsigned int val) {
  mem[0] = val;
  mem[1] = val >> 8;
}

static void mem_put_le32(char *mem, unsigned int val) {
  mem[0] = val;
  mem[1] = val >> 8;
  mem[2] = val >> 16;
  mem[3] = val >> 24;
}

static void usage(char *program_name) {
  printf(
      "Usage: %s [-f frames] [-s skip_frames] [-w width] [-h height] \n\t"
      "[-n rate_num] [-d rate_den] [-b bitrate] [-l layers] [-g gop_size] \n\t"
      "[-z dummy_frame (default 1) \n\t"
      "[-q quantizer (lowest to highest)] \n\t"
      "[-r 1/16th scale factor (lowest to highest layer)] "
      "<input_filename> <output_filename>\n",
      basename(program_name));
  exit(EXIT_FAILURE);
}

static void die(const char *fmt, ...) {
  va_list ap;

  va_start(ap, fmt);
  vprintf(fmt, ap);
  if (fmt[strlen(fmt) - 1] != '\n') printf("\n");
  exit(EXIT_FAILURE);
}

static void die_codec(vpx_codec_ctx_t *ctx, const char *s) {
  const char *detail = vpx_codec_error_detail(ctx);

  printf("%s: %s\n", s, vpx_codec_error(ctx));
  if (detail) printf("    %s\n", detail);
  exit(EXIT_FAILURE);
}

static int read_frame(FILE *f, vpx_image_t *img) {
  size_t nbytes, to_read;
  int res = 1;
  int plane;

  for (plane = 0; plane < 3; plane++) {
    unsigned char *ptr;
    int w = (plane ? (1 + img->d_w) / 2 : img->d_w);
    int h = (plane ? (1 + img->d_h) / 2 : img->d_h);
    int r;

    switch (plane) {
      case 1:
        ptr = img->planes[VPX_PLANE_U];
        break;
      case 2:
        ptr = img->planes[VPX_PLANE_V];
        break;
      default:
        ptr = img->planes[plane];
    }
    for (r = 0; r < h; r++) {
      to_read = w;

      nbytes = fread(ptr, 1, to_read, f);
      if (nbytes != to_read) {
        res = 0;
        if (nbytes > 0)
          printf("Warning: Read partial frame. Check your width & height!\n");
        break;
      }
      ptr += img->stride[plane];
    }
    if (!res) break;
  }
  return res;
}

static int create_dummy_frame(vpx_image_t *img) {
  size_t buf_size;
  buf_size = img->w * img->h * 3 / 2;
  memset(img->planes[0], 129, buf_size);
  return 1;
}

static void write_ivf_file_header(FILE *outfile, unsigned int width,
                                  unsigned int height, int timebase_num,
                                  int timebase_den, int frame_cnt) {
  char header[32];

  header[0] = 'D';
  header[1] = 'K';
  header[2] = 'I';
  header[3] = 'F';
  mem_put_le16(header + 4, 0);             /* version */
  mem_put_le16(header + 6, 32);            /* headersize */
  mem_put_le32(header + 8, fourcc);        /* headersize */
  mem_put_le16(header + 12, width);        /* width */
  mem_put_le16(header + 14, height);       /* height */
  mem_put_le32(header + 16, timebase_den); /* rate */
  mem_put_le32(header + 20, timebase_num); /* scale */
  mem_put_le32(header + 24, frame_cnt);    /* length */
  mem_put_le32(header + 28, 0);            /* unused */

  (void)fwrite(header, 1, 32, outfile);
}

static void write_ivf_frame_header(FILE *outfile, vpx_codec_pts_t pts,
                                   size_t sz) {
  char header[12];
  mem_put_le32(header, sz);
  mem_put_le32(header + 4, pts & 0xFFFFFFFF);
  mem_put_le32(header + 8, pts >> 32);

  (void)fwrite(header, 1, 12, outfile);
}

static void parse_command_line(int argc, char **argv, SvcContext *svc_ctx,
                               vpx_codec_enc_cfg_t *enc_cfg) {
  unsigned int width = 1920;
  unsigned int height = 1080;
  unsigned int timebase_num = 1;
  unsigned int timebase_den = 60;
  unsigned int bitrate = 1000;
  unsigned int number_spatial_layers = 5;
  int use_dummy_frame = 1;

  int c;
  vpx_codec_err_t res;
  int r = 0;
  int q = 0;

  opterr = 0;
  while ((c = getopt(argc, argv, "f:w:h:n:d:b:s:l:g:r:q:z:")) != -1)
    switch (c) {
      case 'f':
        number_frames_to_code = atoi(optarg);
        break;
      case 'w':
        width = atoi(optarg);
        break;
      case 'h':
        height = atoi(optarg);
        break;
      case 'n':
        timebase_num = atoi(optarg);
        break;
      case 'd':
        timebase_den = atoi(optarg);
        break;
      case 'b':
        bitrate = atoi(optarg);
        break;
      case 's':
        number_frames_to_skip = atoi(optarg);
        break;
      case 'l':
        number_spatial_layers = atoi(optarg);
        break;
      case 'g':
        gop_size = atoi(optarg);
        break;
      case 'r':
        scaling_factor = optarg;
        break;
      case 'q':
        quantizer = optarg;
        break;
      case 'z':
        use_dummy_frame = atoi(optarg);
        break;
      case '?':
        usage(argv[0]);
    }

  // Parse required parameters
  if (argc - optind != 2) {
    usage(argv[0]);
  }

  if ((r != 0 && q != 0 && r != q) || (r != number_spatial_layers && r != 0)) {
    usage(argv[0]);
  }

  input_filename = argv[optind];
  output_filename = argv[optind + 1];

  if (width < 16 || width % 2 || height < 16 || height % 2)
    die("Invalid resolution: %d x %d", width, height);

  // initialize SvcContext
  svc_ctx->log_level = SVC_LOG_DEBUG;
  svc_ctx->spatial_layers = number_spatial_layers;
  svc_ctx->encoding_mode = encoding_mode;
  svc_ctx->gop_size = gop_size;
  svc_ctx->quantizer_values = quantizer;
  svc_ctx->scale_factors = scaling_factor;
  // when using a dummy frame, that frame is only encoded to be full size
  svc_ctx->first_frame_full_size = use_dummy_frame;

  /* Populate encoder configuration */
  res = vpx_codec_enc_config_default(interface, enc_cfg, 0);
  if (res) {
    die("Failed to get config: %s\n", vpx_codec_err_to_string(res));
  }
  printf(
      "Codec %s\nframes: %d, skip: %d, layers: %d\n"
      "width %d, height: %d, \n"
      "num: %d, den: %d, bitrate: %d, \n"
      "gop size: %d, use_dummy_frame: %d \n",
      vpx_codec_iface_name(interface), number_frames_to_code,
      number_frames_to_skip, number_spatial_layers, width, height, timebase_num,
      timebase_den, bitrate, gop_size, use_dummy_frame);

  enc_cfg->rc_target_bitrate = bitrate;
  enc_cfg->g_w = width;
  enc_cfg->g_h = height;
  enc_cfg->g_timebase.num = timebase_num;
  enc_cfg->g_timebase.den = timebase_den;
}

int main(int argc, char **argv) {
  FILE *infile, *outfile;
  vpx_codec_ctx_t codec;
  vpx_codec_enc_cfg_t enc_cfg;
  SvcContext svc_ctx;
  int i;
  int frame_cnt = 0;
  vpx_image_t raw;
  clock_t before;
  clock_t after;
  vpx_codec_err_t res;
  int pts = 0;            /* PTS starts at 0 */
  int frame_duration = 1; /* 1 timebase tick per frame */

  memset(&svc_ctx, 0, sizeof(svc_ctx));
  svc_ctx.log_print = 1;
  parse_command_line(argc, argv, &svc_ctx, &enc_cfg);

  // Allocate image buffer
  if (!vpx_img_alloc(&raw, VPX_IMG_FMT_I420, enc_cfg.g_w, enc_cfg.g_h, 32))
    die("Failed to allocate image", enc_cfg.g_w, enc_cfg.g_h);

  if (!(infile = fopen(input_filename, "rb")))
    die("Failed to open %s for reading", argv[1]);

  if (!(outfile = fopen(output_filename, "wb")))
    die("Failed to open %s for writing", output_filename);

  // Initialize codec
  if (vpx_svc_init(&svc_ctx, &codec, interface, &enc_cfg) != VPX_CODEC_OK)
    die("Failed to initialize encoder");

  write_ivf_file_header(outfile, enc_cfg.g_w, enc_cfg.g_h,
                        enc_cfg.g_timebase.num, enc_cfg.g_timebase.den, 0);

  // skip initial frames
  for (i = 0; i < number_frames_to_skip; i++) {
    read_frame(infile, &raw);
  }

  before = clock();
  // Encode frames
  while (frame_cnt <= number_frames_to_code) {
    if (frame_cnt == 0 && svc_ctx.first_frame_full_size) {
      create_dummy_frame(&raw);
    } else {
      if (!read_frame(infile, &raw)) break;
    }
    res = vpx_svc_encode(&svc_ctx, &codec, &raw, pts, frame_duration,
                         VPX_DL_REALTIME);
    printf("%s", svc_get_message(&svc_ctx));
    if (res != VPX_CODEC_OK) {
      die_codec(&codec, "Failed to encode frame");
    }
    if (svc_get_frame_size(&svc_ctx) > 0) {
      write_ivf_frame_header(outfile, pts, svc_get_frame_size(&svc_ctx));
      (void)fwrite(svc_get_buffer(&svc_ctx), 1, svc_get_frame_size(&svc_ctx),
                   outfile);
    }
    frame_cnt++;
    pts += frame_duration;
  }  // end encode frames loop

  after = clock();
  printf("Processed %d frames in different resolutions in %ld ms.\n",
         frame_cnt - svc_ctx.first_frame_full_size,
         (int)(after - before) / (CLOCKS_PER_SEC / 1000));

  fclose(infile);
  if (vpx_codec_destroy(&codec)) die_codec(&codec, "Failed to destroy codec");

  // rewrite the output file headers with the actual frame count
  if (!fseek(outfile, 0, SEEK_SET)) {
    write_ivf_file_header(outfile, enc_cfg.g_w, enc_cfg.g_h,
                          enc_cfg.g_timebase.num, enc_cfg.g_timebase.den,
                          frame_cnt);
  }
  fclose(outfile);

  // display average size, psnr
  svc_dump_statistics(&svc_ctx);
  printf("%s", svc_get_message(&svc_ctx));

  return EXIT_SUCCESS;
}
