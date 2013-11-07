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

#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "./args.h"
#include "vpx/svc_context.h"
#include "vpx/vp8cx.h"
#include "vpx/vpx_encoder.h"

#define VP90_FOURCC 0x30395056

static const struct arg_enum_list encoding_mode_enum[] = {
  {"i", INTER_LAYER_PREDICTION_I},
  {"alt-ip", ALT_INTER_LAYER_PREDICTION_IP},
  {"ip", INTER_LAYER_PREDICTION_IP},
  {"gf", USE_GOLDEN_FRAME},
  {NULL, 0}
};

static const arg_def_t encoding_mode_arg = ARG_DEF_ENUM(
    "m", "encoding-mode", 1, "Encoding mode algorithm", encoding_mode_enum);
static const arg_def_t skip_frames_arg =
    ARG_DEF("s", "skip-frames", 1, "input frames to skip");
static const arg_def_t frames_arg =
    ARG_DEF("f", "frames", 1, "number of frames to encode");
static const arg_def_t width_arg = ARG_DEF("w", "width", 1, "source width");
static const arg_def_t height_arg = ARG_DEF("h", "height", 1, "source height");
static const arg_def_t timebase_arg =
    ARG_DEF("t", "timebase", 1, "timebase (num/den)");
static const arg_def_t bitrate_arg = ARG_DEF(
    "b", "target-bitrate", 1, "encoding bitrate, in kilobits per second");
static const arg_def_t layers_arg =
    ARG_DEF("l", "layers", 1, "number of SVC layers");
static const arg_def_t kf_dist_arg =
    ARG_DEF("k", "kf-dist", 1, "number of frames between keyframes");
static const arg_def_t scale_factors_arg =
    ARG_DEF("r", "scale-factors", 1, "scale factors (lowest to highest layer)");
static const arg_def_t quantizers_arg =
    ARG_DEF("q", "quantizers", 1, "quantizers (lowest to highest layer)");
static const arg_def_t dummy_frame_arg =
    ARG_DEF("z", "dummy-frame", 1, "make first frame blank and full size");

static const arg_def_t *svc_args[] = {
  &encoding_mode_arg, &frames_arg,        &width_arg,       &height_arg,
  &timebase_arg,      &bitrate_arg,       &skip_frames_arg, &layers_arg,
  &kf_dist_arg,       &scale_factors_arg, &quantizers_arg,  &dummy_frame_arg,
  NULL
};

static const SVC_ENCODING_MODE default_encoding_mode =
    INTER_LAYER_PREDICTION_IP;
static const uint32_t default_frames_to_skip = 0;
static const uint32_t default_frames_to_code = 60 * 60;
static const uint32_t default_width = 1920;
static const uint32_t default_height = 1080;
static const uint32_t default_timebase_num = 1;
static const uint32_t default_timebase_den = 60;
static const uint32_t default_bitrate = 1000;
static const uint32_t default_spatial_layers = 5;
static const uint32_t default_kf_dist = 100;
static const int default_use_dummy_frame = 1;

typedef struct {
  char *input_filename;
  char *output_filename;
  uint32_t frames_to_code;
  uint32_t frames_to_skip;
} AppInput;

static void mem_put_le16(char *mem, uint32_t val) {
  mem[0] = val;
  mem[1] = val >> 8;
}

static void mem_put_le32(char *mem, uint32_t val) {
  mem[0] = val;
  mem[1] = val >> 8;
  mem[2] = val >> 16;
  mem[3] = val >> 24;
}

static void usage(const char *exec_name) {
  fprintf(stderr, "Usage: %s <options> input_filename output_filename\n",
          exec_name);
  fprintf(stderr, "Options:\n");
  arg_show_usage(stderr, svc_args);
  exit(EXIT_FAILURE);
}

void die(const char *fmt, ...) {
  va_list ap;

  va_start(ap, fmt);
  vfprintf(stderr, fmt, ap);
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
  size_t nbytes;
  int res = 1;
  int plane;

  for (plane = 0; plane < 3; ++plane) {
    uint8_t *ptr;
    const int w = (plane ? (1 + img->d_w) / 2 : img->d_w);
    const int h = (plane ? (1 + img->d_h) / 2 : img->d_h);
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
    for (r = 0; r < h; ++r) {
      const int to_read = w;

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
  const size_t buf_size = img->w * img->h * 3 / 2;
  memset(img->planes[0], 129, buf_size);
  return 1;
}

static void write_ivf_file_header(FILE *outfile,
                                  uint32_t width, uint32_t height,
                                  int timebase_num, int timebase_den,
                                  int frame_cnt) {
  char header[32];

  header[0] = 'D';
  header[1] = 'K';
  header[2] = 'I';
  header[3] = 'F';
  mem_put_le16(header + 4, 0);             /* version */
  mem_put_le16(header + 6, 32);            /* headersize */
  mem_put_le32(header + 8, VP90_FOURCC);   /* fourcc */
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
  mem_put_le32(header, (uint32_t)sz);
  mem_put_le32(header + 4, pts & 0xFFFFFFFF);
  mem_put_le32(header + 8, pts >> 32);

  (void)fwrite(header, 1, 12, outfile);
}

static void parse_command_line(int argc, const char **argv_,
                               AppInput *app_input, SvcContext *svc_ctx,
                               vpx_codec_enc_cfg_t *enc_cfg) {
  struct arg arg;
  char **argv, **argi, **argj;
  vpx_codec_err_t res;

  // initialize SvcContext with parameters that will be passed to vpx_svc_init
  svc_ctx->log_level = SVC_LOG_DEBUG;
  svc_ctx->spatial_layers = default_spatial_layers;
  svc_ctx->encoding_mode = default_encoding_mode;
  // when using a dummy frame, that frame is only encoded to be full size
  svc_ctx->first_frame_full_size = default_use_dummy_frame;

  // start with default encoder configuration
  res = vpx_codec_enc_config_default(vpx_codec_vp9_cx(), enc_cfg, 0);
  if (res) {
    die("Failed to get config: %s\n", vpx_codec_err_to_string(res));
  }
  // update enc_cfg with app default values
  enc_cfg->g_w = default_width;
  enc_cfg->g_h = default_height;
  enc_cfg->g_timebase.num = default_timebase_num;
  enc_cfg->g_timebase.den = default_timebase_den;
  enc_cfg->rc_target_bitrate = default_bitrate;
  enc_cfg->kf_min_dist = default_kf_dist;
  enc_cfg->kf_max_dist = default_kf_dist;

  // initialize AppInput with default values
  app_input->frames_to_code = default_frames_to_code;
  app_input->frames_to_skip = default_frames_to_skip;

  // process command line options
  argv = argv_dup(argc - 1, argv_ + 1);
  for (argi = argj = argv; (*argj = *argi); argi += arg.argv_step) {
    arg.argv_step = 1;

    if (arg_match(&arg, &encoding_mode_arg, argi)) {
      svc_ctx->encoding_mode = arg_parse_enum_or_int(&arg);
    } else if (arg_match(&arg, &frames_arg, argi)) {
      app_input->frames_to_code = arg_parse_uint(&arg);
    } else if (arg_match(&arg, &width_arg, argi)) {
      enc_cfg->g_w = arg_parse_uint(&arg);
    } else if (arg_match(&arg, &height_arg, argi)) {
      enc_cfg->g_h = arg_parse_uint(&arg);
    } else if (arg_match(&arg, &height_arg, argi)) {
      enc_cfg->g_h = arg_parse_uint(&arg);
    } else if (arg_match(&arg, &timebase_arg, argi)) {
      enc_cfg->g_timebase = arg_parse_rational(&arg);
    } else if (arg_match(&arg, &bitrate_arg, argi)) {
      enc_cfg->rc_target_bitrate = arg_parse_uint(&arg);
    } else if (arg_match(&arg, &skip_frames_arg, argi)) {
      app_input->frames_to_skip = arg_parse_uint(&arg);
    } else if (arg_match(&arg, &layers_arg, argi)) {
      svc_ctx->spatial_layers = arg_parse_uint(&arg);
    } else if (arg_match(&arg, &kf_dist_arg, argi)) {
      enc_cfg->kf_min_dist = arg_parse_uint(&arg);
      enc_cfg->kf_max_dist = enc_cfg->kf_min_dist;
    } else if (arg_match(&arg, &scale_factors_arg, argi)) {
      vpx_svc_set_scale_factors(svc_ctx, arg.val);
    } else if (arg_match(&arg, &quantizers_arg, argi)) {
      vpx_svc_set_quantizers(svc_ctx, arg.val);
    } else if (arg_match(&arg, &dummy_frame_arg, argi)) {
      svc_ctx->first_frame_full_size = arg_parse_int(&arg);
    } else {
      ++argj;
    }
  }

  // Check for unrecognized options
  for (argi = argv; *argi; ++argi)
    if (argi[0][0] == '-' && strlen(argi[0]) > 1)
      die("Error: Unrecognized option %s\n", *argi);

  if (argv[0] == NULL || argv[1] == 0) {
    usage(argv_[0]);
  }
  app_input->input_filename = argv[0];
  app_input->output_filename = argv[1];
  free(argv);

  if (enc_cfg->g_w < 16 || enc_cfg->g_w % 2 || enc_cfg->g_h < 16 ||
      enc_cfg->g_h % 2)
    die("Invalid resolution: %d x %d\n", enc_cfg->g_w, enc_cfg->g_h);

  printf(
      "Codec %s\nframes: %d, skip: %d\n"
      "mode: %d, layers: %d\n"
      "width %d, height: %d,\n"
      "num: %d, den: %d, bitrate: %d,\n"
      "gop size: %d, use_dummy_frame: %d\n",
      vpx_codec_iface_name(vpx_codec_vp9_cx()), app_input->frames_to_code,
      app_input->frames_to_skip, svc_ctx->encoding_mode,
      svc_ctx->spatial_layers, enc_cfg->g_w, enc_cfg->g_h,
      enc_cfg->g_timebase.num, enc_cfg->g_timebase.den,
      enc_cfg->rc_target_bitrate, enc_cfg->kf_max_dist,
      svc_ctx->first_frame_full_size);
}

int main(int argc, const char **argv) {
  AppInput app_input = {0};
  FILE *infile, *outfile;
  vpx_codec_ctx_t codec;
  vpx_codec_enc_cfg_t enc_cfg;
  SvcContext svc_ctx;
  uint32_t i;
  uint32_t frame_cnt = 0;
  vpx_image_t raw;
  vpx_codec_err_t res;
  int pts = 0;            /* PTS starts at 0 */
  int frame_duration = 1; /* 1 timebase tick per frame */

  memset(&svc_ctx, 0, sizeof(svc_ctx));
  svc_ctx.log_print = 1;
  parse_command_line(argc, argv, &app_input, &svc_ctx, &enc_cfg);

  // Allocate image buffer
  if (!vpx_img_alloc(&raw, VPX_IMG_FMT_I420, enc_cfg.g_w, enc_cfg.g_h, 32))
    die("Failed to allocate image %dx%d\n", enc_cfg.g_w, enc_cfg.g_h);

  if (!(infile = fopen(app_input.input_filename, "rb")))
    die("Failed to open %s for reading\n", app_input.input_filename);

  if (!(outfile = fopen(app_input.output_filename, "wb")))
    die("Failed to open %s for writing\n", app_input.output_filename);

  // Initialize codec
  if (vpx_svc_init(&svc_ctx, &codec, vpx_codec_vp9_cx(), &enc_cfg) !=
      VPX_CODEC_OK)
    die("Failed to initialize encoder\n");

  write_ivf_file_header(outfile, enc_cfg.g_w, enc_cfg.g_h,
                        enc_cfg.g_timebase.num, enc_cfg.g_timebase.den, 0);

  // skip initial frames
  for (i = 0; i < app_input.frames_to_skip; ++i) {
    read_frame(infile, &raw);
  }

  // Encode frames
  while (frame_cnt <= app_input.frames_to_code) {
    if (frame_cnt == 0 && svc_ctx.first_frame_full_size) {
      create_dummy_frame(&raw);
    } else {
      if (!read_frame(infile, &raw)) break;
    }
    res = vpx_svc_encode(&svc_ctx, &codec, &raw, pts, frame_duration,
                         VPX_DL_REALTIME);
    printf("%s", vpx_svc_get_message(&svc_ctx));
    if (res != VPX_CODEC_OK) {
      die_codec(&codec, "Failed to encode frame");
    }
    if (vpx_svc_get_frame_size(&svc_ctx) > 0) {
      write_ivf_frame_header(outfile, pts, vpx_svc_get_frame_size(&svc_ctx));
      (void)fwrite(vpx_svc_get_buffer(&svc_ctx), 1,
                   vpx_svc_get_frame_size(&svc_ctx), outfile);
    }
    ++frame_cnt;
    pts += frame_duration;
  }

  printf("Processed %d frames\n", frame_cnt - svc_ctx.first_frame_full_size);

  fclose(infile);
  if (vpx_codec_destroy(&codec)) die_codec(&codec, "Failed to destroy codec");

  // rewrite the output file headers with the actual frame count
  if (!fseek(outfile, 0, SEEK_SET)) {
    write_ivf_file_header(outfile, enc_cfg.g_w, enc_cfg.g_h,
                          enc_cfg.g_timebase.num, enc_cfg.g_timebase.den,
                          frame_cnt);
  }
  fclose(outfile);
  vpx_img_free(&raw);

  // display average size, psnr
  printf("%s", vpx_svc_dump_statistics(&svc_ctx));

  vpx_svc_release(&svc_ctx);

  return EXIT_SUCCESS;
}
