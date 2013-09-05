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
#define interface (vpx_codec_vp9_cx())
#define fourcc 0x30395056
#define IVF_FILE_HDR_SZ (32)
#define IVF_FRAME_HDR_SZ (12)
#define NUM_BUFFERS 8

char *input_filename;
char *output_filename;
unsigned int number_frames_to_code = 60 * 60;
unsigned int number_frames_to_skip = 0;
unsigned int number_spatial_layers = 5;
unsigned int key_period = 100;

typedef enum ENCODING_MODE {
  INTER_LAYER_PREDICTION_I,
  INTER_LAYER_PREDICTION_IP,
  USE_GOLDEN_FRAME
} ENCODING_MODE;

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
      "[-n rate_num] [-d rate_den] [-b bitrate] [-l layers] "
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

  to_read = img->w * img->h * 3 / 2;
  nbytes = fread(img->planes[0], 1, to_read, f);
  if (nbytes != to_read) {
    res = 0;
    if (nbytes > 0)
      printf("Warning: Read partial frame. Check your width & height!\n");
  }
  return res;
}

static int read_dummy_frame(vpx_image_t *img) {
  size_t to_read;

  to_read = img->w * img->h * 3 / 2;
  memset(img->planes[0], 129, to_read);
  return 1;
}

static void write_ivf_file_header(FILE *outfile, const vpx_codec_enc_cfg_t *cfg,
                                  int frame_cnt) {
  char header[32];

  if (cfg->g_pass != VPX_RC_ONE_PASS && cfg->g_pass != VPX_RC_LAST_PASS) return;
  header[0] = 'D';
  header[1] = 'K';
  header[2] = 'I';
  header[3] = 'F';
  mem_put_le16(header + 4, 0);                    /* version */
  mem_put_le16(header + 6, 32);                   /* headersize */
  mem_put_le32(header + 8, fourcc);               /* headersize */
  mem_put_le16(header + 12, cfg->g_w);            /* width */
  mem_put_le16(header + 14, cfg->g_h);            /* height */
  mem_put_le32(header + 16, cfg->g_timebase.den); /* rate */
  mem_put_le32(header + 20, cfg->g_timebase.num); /* scale */
  mem_put_le32(header + 24, frame_cnt);           /* length */
  mem_put_le32(header + 28, 0);                   /* unused */

  (void)fwrite(header, 1, 32, outfile);
}

static void write_ivf_frame_header(FILE *outfile,
                                   const vpx_codec_cx_pkt_t *pkt) {
  char header[12];
  vpx_codec_pts_t pts;

  if (pkt->kind != VPX_CODEC_CX_FRAME_PKT) return;

  pts = pkt->data.frame.pts;
  mem_put_le32(header, pkt->data.frame.sz);
  mem_put_le32(header + 4, pts & 0xFFFFFFFF);
  mem_put_le32(header + 8, pts >> 32);

  (void)fwrite(header, 1, 12, outfile);
}

static void check_parameters() {
  if (number_spatial_layers > 5) die("Cannot support more than 5 layers");
}

static void parse_command_line(int argc, char **argv,
                               vpx_codec_enc_cfg_t *cfg) {
  unsigned int width = 1920;
  unsigned int height = 1080;
  unsigned int timebase_num = 1;
  unsigned int timebase_den = 60;
  unsigned int bitrate = 1000;
  int c;
  vpx_codec_err_t res;

  opterr = 0;
  while ((c = getopt(argc, argv, "f:w:h:n:d:b:s:l:p:")) != -1) switch (c) {
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
      case 'p':
        key_period = atoi(optarg);
        break;
      case '?':
        usage(argv[0]);
    }

  // Parse required parameters
  if (argc - optind != 2) {
    usage(argv[0]);
  }

  input_filename = argv[optind];
  output_filename = argv[optind + 1];

  if (width < 16 || width % 2 || height < 16 || height % 2)
    die("Invalid resolution: %d x %d", width, height);

  /* Populate encoder configuration */
  res = vpx_codec_enc_config_default(interface, cfg, 0);
  if (res) {
    die("Failed to get config: %s\n", vpx_codec_err_to_string(res));
  }
  printf(
      "Codec %s\nframes: %d, skip: %d, layers: %d\n"
      "width %d, height: %d, \n"
      "num: %d, den: %d, bitrate: %d, \n"
      "key period: %d \n",
      vpx_codec_iface_name(interface), number_frames_to_code,
      number_frames_to_skip, number_spatial_layers, width, height, timebase_num,
      timebase_den, bitrate, key_period);

  // Do minimal check at the application level. Encoder parameters will be
  // checked internally
  check_parameters();

  cfg->rc_target_bitrate = bitrate;
  cfg->g_w = width;
  cfg->g_h = height;
  cfg->g_timebase.num = timebase_num;
  cfg->g_timebase.den = timebase_den;
  cfg->ss_number_layers = number_spatial_layers;
}

static void set_default_configuration(vpx_codec_enc_cfg_t *cfg) {
  /* Real time parameters */
  cfg->rc_dropframe_thresh = 0;
  cfg->rc_end_usage = VPX_CBR;
  cfg->rc_resize_allowed = 0;
  cfg->rc_min_quantizer = 33;
  cfg->rc_max_quantizer = 33;
  cfg->rc_undershoot_pct = 100;
  cfg->rc_overshoot_pct = 15;
  cfg->rc_buf_initial_sz = 500;
  cfg->rc_buf_optimal_sz = 600;
  cfg->rc_buf_sz = 1000;

  /* Enable error resilient mode */
  cfg->g_error_resilient = 1;
  cfg->g_lag_in_frames = 0;

  /* Disable automatic keyframe placement */
  cfg->kf_mode = VPX_KF_DISABLED;
  cfg->kf_min_dist = cfg->kf_max_dist = 3000;
}

static void initialize_codec(vpx_codec_ctx_t *codec, vpx_codec_enc_cfg_t *cfg) {
  int max_intra_size_pct;

  /* Initialize codec */
  if (vpx_codec_enc_init(codec, interface, cfg, VPX_CODEC_USE_PSNR))
    die_codec(codec, "Failed to initialize encoder");

  vpx_codec_control(codec, VP9E_SET_SVC, 1);
  /* Cap CPU & first I-frame size */
  vpx_codec_control(codec, VP8E_SET_CPUUSED, 1);
  vpx_codec_control(codec, VP8E_SET_STATIC_THRESHOLD, 1);
  vpx_codec_control(codec, VP8E_SET_NOISE_SENSITIVITY, 1);
  vpx_codec_control(codec, VP8E_SET_TOKEN_PARTITIONS, 1);

  max_intra_size_pct =
      (int)(((double)cfg->rc_buf_optimal_sz * 0.5) *
            ((double)cfg->g_timebase.den / cfg->g_timebase.num) / 10.0);
  /* printf ("max_intra_size_pct=%d\n", max_intra_size_pct); */

  vpx_codec_control(codec, VP8E_SET_MAX_INTRA_BITRATE_PCT, max_intra_size_pct);
}

static int calculate_layer(int frame_cnt, int number_spatial_layers) {
  if (frame_cnt == 0)
    return 0;
  else
    return (frame_cnt + number_spatial_layers - 1) % number_spatial_layers;
}

static void switch_to_layer(int layer, unsigned int initial_width,
                            unsigned int initial_height,
                            vpx_codec_ctx_t *codec) {
  // Set layer size
  int scaling_factor_num[MAX_LAYERS] = {2, 1, 4, 2, 1};
  int scaling_factor_den[MAX_LAYERS] = {9, 3, 9, 3, 1};

  int quantizer[MAX_LAYERS] = {60, 53, 39, 33, 27};

  unsigned int current_width;
  unsigned int current_height;

  current_width = initial_width *
                  scaling_factor_num[layer + 5 - number_spatial_layers] /
                  scaling_factor_den[layer + 5 - number_spatial_layers];
  current_height = initial_height *
                   scaling_factor_num[layer + 5 - number_spatial_layers] /
                   scaling_factor_den[layer + 5 - number_spatial_layers];

  current_width += current_width % 2;
  current_height += current_height % 2;

  vpx_codec_control(codec, VP9E_SET_WIDTH, &current_width);
  vpx_codec_control(codec, VP9E_SET_HEIGHT, &current_height);

  // Set layer context
  vpx_codec_control(codec, VP9E_SET_LAYER, &layer);
  vpx_codec_control(codec, VP9E_SET_MAX_Q,
                    quantizer[layer + 5 - number_spatial_layers]);
  vpx_codec_control(codec, VP9E_SET_MIN_Q,
                    quantizer[layer + 5 - number_spatial_layers]);
}

static int get_flag(int is_I_frame_in_layer, int layer, ENCODING_MODE mode) {
  // First layer
  switch (mode) {
    case INTER_LAYER_PREDICTION_I:
      if (is_I_frame_in_layer && layer == 0) return VPX_EFLAG_FORCE_KF;
      if (layer == 0)
        return VP8_EFLAG_NO_UPD_GF | VP8_EFLAG_NO_UPD_ARF |
               VP8_EFLAG_NO_REF_GF | VP8_EFLAG_NO_REF_ARF;
      else if (is_I_frame_in_layer)
        return VP8_EFLAG_NO_UPD_GF | VP8_EFLAG_NO_UPD_ARF |
               VP8_EFLAG_NO_REF_GF | VP8_EFLAG_NO_REF_LAST;
      else
        return VP8_EFLAG_NO_UPD_GF | VP8_EFLAG_NO_UPD_ARF |
               VP8_EFLAG_NO_REF_GF | VP8_EFLAG_NO_REF_ARF;
      break;

    case INTER_LAYER_PREDICTION_IP:
      if (is_I_frame_in_layer && layer == 0) return VPX_EFLAG_FORCE_KF;
      if (layer == 0)
        return VP8_EFLAG_NO_UPD_GF | VP8_EFLAG_NO_UPD_ARF |
               VP8_EFLAG_NO_REF_GF | VP8_EFLAG_NO_REF_ARF;
      else if (is_I_frame_in_layer)
        return VP8_EFLAG_NO_UPD_GF | VP8_EFLAG_NO_UPD_ARF |
               VP8_EFLAG_NO_REF_GF | VP8_EFLAG_NO_REF_LAST;
      else
        return VP8_EFLAG_NO_UPD_GF | VP8_EFLAG_NO_UPD_ARF | VP8_EFLAG_NO_REF_GF;
      break;

    case USE_GOLDEN_FRAME:
      if (is_I_frame_in_layer && layer == 0) return VPX_EFLAG_FORCE_KF;
      if (2 * number_spatial_layers - NUM_BUFFERS <= layer) {
        if (layer == 0)
          return VP8_EFLAG_NO_UPD_GF | VP8_EFLAG_NO_UPD_ARF |
                 VP8_EFLAG_NO_REF_ARF;
        else if (is_I_frame_in_layer)
          return VP8_EFLAG_NO_UPD_ARF | VP8_EFLAG_NO_REF_GF |
                 VP8_EFLAG_NO_REF_LAST;
        else
          return VP8_EFLAG_NO_UPD_GF | VP8_EFLAG_NO_UPD_ARF;
      } else {
        if (layer == 0)
          return VP8_EFLAG_NO_UPD_GF | VP8_EFLAG_NO_UPD_ARF |
                 VP8_EFLAG_NO_REF_GF | VP8_EFLAG_NO_REF_ARF;
        else if (is_I_frame_in_layer)
          return VP8_EFLAG_NO_UPD_GF | VP8_EFLAG_NO_UPD_ARF |
                 VP8_EFLAG_NO_REF_GF | VP8_EFLAG_NO_REF_LAST;
        else
          return VP8_EFLAG_NO_UPD_GF | VP8_EFLAG_NO_UPD_ARF |
                 VP8_EFLAG_NO_REF_GF | VP8_EFLAG_NO_REF_ARF;
      }
      break;
    default:
      return VPX_EFLAG_FORCE_KF;
  }
}

int main(int argc, char **argv) {
  FILE *infile, *outfile[MAX_LAYERS];
  vpx_codec_ctx_t codec;
  vpx_codec_enc_cfg_t cfg;
  int frame_cnt = 0;
  vpx_image_t raw;
  int frame_avail = 1;
  int got_data = 0;
  int i;
  int frames_in_layer[MAX_LAYERS] = {0};
  clock_t before;
  clock_t after;
  int pts = 0;            /* PTS starts at 0 */
  int frame_duration = 1; /* 1 timebase tick per frame */

  parse_command_line(argc, argv, &cfg);

  // Allocate image buffer
  if (!vpx_img_alloc(&raw, VPX_IMG_FMT_I420, cfg.g_w, cfg.g_h, 32))
    die("Failed to allocate image", cfg.g_w, cfg.g_h);

  set_default_configuration(&cfg);

  /* Open input file */
  if (!(infile = fopen(input_filename, "rb")))
    die("Failed to open %s for reading", argv[1]);

  /* Open output file  */
  for (i = 0; i < number_spatial_layers; i++) {
    char file_name[512];
    snprintf(file_name, sizeof(file_name), "%s_%d.ivf", output_filename, i);
    if (!(outfile[i] = fopen(file_name, "wb")))
      die("Failed to open %s for writing", file_name);
    write_ivf_file_header(outfile[i], &cfg, 0);
  }

  initialize_codec(&codec, &cfg);

  // skip initial frames
  for (i = 0; i < number_frames_to_skip; i++) {
    read_frame(infile, &raw);
  }

  before = clock();
  // Encoding frames
  while ((frame_avail || got_data) &&
         frame_cnt <= number_frames_to_code * number_spatial_layers) {
    int flags = 0;
    vpx_codec_iter_t iter = NULL;
    const vpx_codec_cx_pkt_t *pkt;

    int layer = calculate_layer(frame_cnt, number_spatial_layers);
    int is_I_frame_in_layer =
        (((frame_cnt - 1) / number_spatial_layers % key_period) == 0);
    int is_dummy = (frame_cnt == 0);

    if (is_dummy) {  // Dummy frame
      flags = VPX_EFLAG_FORCE_KF;
      frame_avail = read_dummy_frame(&raw);

    } else {  // Regular frame
      // Read a new frame only at the base layer
      if (layer == 0) frame_avail = read_frame(infile, &raw);
      switch_to_layer(layer, cfg.g_w, cfg.g_h, &codec);
      flags = get_flag(is_I_frame_in_layer, layer, INTER_LAYER_PREDICTION_I);
    }

    // Actual Encoding
    if (vpx_codec_encode(&codec, frame_avail ? &raw : NULL, pts, 1, flags,
                         VPX_DL_REALTIME))
      die_codec(&codec, "Failed to encode frame");

    got_data = 0;
    // Process data / Get PSNR statistics
    while ((pkt = vpx_codec_get_cx_data(&codec, &iter))) {
      got_data = 1;
      switch (pkt->kind) {
        case VPX_CODEC_CX_FRAME_PKT:
          for (i = layer; i < number_spatial_layers; i++) {
            write_ivf_frame_header(outfile[i], pkt);
            (void)fwrite(pkt->data.frame.buf, 1, pkt->data.frame.sz,
                         outfile[i]);
            frames_in_layer[i]++;
          }
          break;
        case VPX_CODEC_PSNR_PKT:
          if (frame_cnt != 0)
            printf(
                "Processed Frame %d, layer %d, PSNR(Total/Y/U/V): "
                "%2.3f  %2.3f  %2.3f  %2.3f \n",
                (frame_cnt - 1) / number_spatial_layers + 1, layer,
                pkt->data.psnr.psnr[0], pkt->data.psnr.psnr[1],
                pkt->data.psnr.psnr[2], pkt->data.psnr.psnr[3]);
          break;
        default:
          break;
      }
    }
    frame_cnt++;
    // TODO(ivan): Modify ts later if(!layer)
    pts += frame_duration;
  }
  // end while

  after = clock();
  printf("Processed %d frames in different resolutions in %ld ms.\n",
         frame_cnt - 1, (int)(after - before) / (CLOCKS_PER_SEC / 1000));

  fclose(infile);

  if (vpx_codec_destroy(&codec)) die_codec(&codec, "Failed to destroy codec");

  /* Try to rewrite the output file headers with the actual frame count */
  for (i = 0; i < number_spatial_layers; i++) {
    if (!fseek(outfile[i], 0, SEEK_SET)) {
      write_ivf_file_header(outfile[i], &cfg, frames_in_layer[i]);
    }
    fclose(outfile[i]);
  }

  return EXIT_SUCCESS;
}
