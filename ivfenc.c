/*
 *  Copyright (c) 2010 The VP8 project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license and patent
 *  grant that can be found in the LICENSE file in the root of the source
 *  tree. All contributing project authors may be found in the AUTHORS
 *  file in the root of the source tree.
 */


/* This is a simple program that encodes YV12 files and generates ivf
 * files using the new interface.
 */
#define USE_POSIX_MMAP HAVE_SYS_MMAN_H

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include "vpx_encoder.h"
#if USE_POSIX_MMAP
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <fcntl.h>
#include <unistd.h>
#endif
#if CONFIG_VP8_ENCODER
#include "vp8cx.h"
#endif
#include "vpx_ports/mem_ops.h"
#include "vpx_ports/vpx_timer.h"

static const char *exec_name;

static const struct codec_item
{
    char const              *name;
    const vpx_codec_iface_t *iface;
    unsigned int             fourcc;
} codecs[] =
{
#if CONFIG_VP8_ENCODER
    {"vp8",  &vpx_codec_vp8_cx_algo, 0x30385056},
#endif
};

static void usage_exit();

void die(const char *fmt, ...)
{
    va_list ap;
    va_start(ap, fmt);
    vprintf(fmt, ap);
    printf("\n");
    usage_exit();
}

static void ctx_exit_on_error(vpx_codec_ctx_t *ctx, const char *s)
{
    if (ctx->err)
    {
        const char *detail = vpx_codec_error_detail(ctx);

        printf("%s: %s\n", s, vpx_codec_error(ctx));

        if (detail)
            printf("    %s\n", detail);

        exit(EXIT_FAILURE);
    }
}

/* This structure is used to abstract the different ways of handling
 * first pass statistics.
 */
typedef struct
{
    vpx_fixed_buf_t buf;
    int             pass;
    FILE           *file;
    char           *buf_ptr;
    size_t          buf_alloc_sz;
} stats_io_t;

int stats_open_file(stats_io_t *stats, const char *fpf, int pass)
{
    int res;

    stats->pass = pass;

    if (pass == 0)
    {
        stats->file = fopen(fpf, "wb");
        stats->buf.sz = 0;
        stats->buf.buf = NULL,
                   res = (stats->file != NULL);
    }
    else
    {
#if 0
#elif USE_POSIX_MMAP
        struct stat stat_buf;
        int fd;

        fd = open(fpf, O_RDONLY);
        stats->file = fdopen(fd, "rb");
        fstat(fd, &stat_buf);
        stats->buf.sz = stat_buf.st_size;
        stats->buf.buf = mmap(NULL, stats->buf.sz, PROT_READ, MAP_PRIVATE,
                              fd, 0);
        res = (stats->buf.buf != NULL);
#else
        size_t nbytes;

        stats->file = fopen(fpf, "rb");

        if (fseek(stats->file, 0, SEEK_END))
        {
            fprintf(stderr, "First-pass stats file must be seekable!\n");
            exit(EXIT_FAILURE);
        }

        stats->buf.sz = stats->buf_alloc_sz = ftell(stats->file);
        rewind(stats->file);

        stats->buf.buf = malloc(stats->buf_alloc_sz);

        if (!stats->buf.buf)
        {
            fprintf(stderr, "Failed to allocate first-pass stats buffer (%d bytes)\n",
                    stats->buf_alloc_sz);
            exit(EXIT_FAILURE);
        }

        nbytes = fread(stats->buf.buf, 1, stats->buf.sz, stats->file);
        res = (nbytes == stats->buf.sz);
#endif
    }

    return res;
}

int stats_open_mem(stats_io_t *stats, int pass)
{
    int res;
    stats->pass = pass;

    if (!pass)
    {
        stats->buf.sz = 0;
        stats->buf_alloc_sz = 64 * 1024;
        stats->buf.buf = malloc(stats->buf_alloc_sz);
    }

    stats->buf_ptr = stats->buf.buf;
    res = (stats->buf.buf != NULL);
    return res;
}


void stats_close(stats_io_t *stats)
{
    if (stats->file)
    {
        if (stats->pass == 1)
        {
#if 0
#elif USE_POSIX_MMAP
            munmap(stats->buf.buf, stats->buf.sz);
#else
            free(stats->buf.buf);
#endif
        }

        fclose(stats->file);
        stats->file = NULL;
    }
    else
    {
        if (stats->pass == 1)
            free(stats->buf.buf);
    }
}

void stats_write(stats_io_t *stats, const void *pkt, size_t len)
{
    if (stats->file)
    {
        fwrite(pkt, 1, len, stats->file);
    }
    else
    {
        if (stats->buf.sz + len > stats->buf_alloc_sz)
        {
            size_t  new_sz = stats->buf_alloc_sz + 64 * 1024;
            char   *new_ptr = realloc(stats->buf.buf, new_sz);

            if (new_ptr)
            {
                stats->buf_ptr = new_ptr + (stats->buf_ptr - (char *)stats->buf.buf);
                stats->buf.buf = new_ptr;
                stats->buf_alloc_sz = new_sz;
            } /* else ... */
        }

        memcpy(stats->buf_ptr, pkt, len);
        stats->buf.sz += len;
        stats->buf_ptr += len;
    }
}

vpx_fixed_buf_t stats_get(stats_io_t *stats)
{
    return stats->buf;
}

#define IVF_FRAME_HDR_SZ (4+8) /* 4 byte size + 8 byte timestamp */
static int read_frame(FILE *f, vpx_image_t *img, unsigned int is_ivf)
{
    int plane = 0;

    if (is_ivf)
    {
        char junk[IVF_FRAME_HDR_SZ];

        /* Skip the frame header. We know how big the frame should be. See
         * write_ivf_frame_header() for documentation on the frame header
         * layout.
         */
        fread(junk, 1, IVF_FRAME_HDR_SZ, f);
    }

    for (plane = 0; plane < 3; plane++)
    {
        unsigned char *ptr;
        int w = (plane ? (1 + img->d_w) / 2 : img->d_w);
        int h = (plane ? (1 + img->d_h) / 2 : img->d_h);
        int r;

        /* Determine the correct plane based on the image format. The for-loop
         * always counts in Y,U,V order, but this may not match the order of
         * the data on disk.
         */
        switch (plane)
        {
        case 1:
            ptr = img->planes[img->fmt==IMG_FMT_YV12? PLANE_V : PLANE_U];
            break;
        case 2:
            ptr = img->planes[img->fmt==IMG_FMT_YV12?PLANE_U : PLANE_V];
            break;
        default:
            ptr = img->planes[plane];
        }

        for (r = 0; r < h; r++)
        {
            fread(ptr, 1, w, f);
            ptr += img->stride[plane];
        }
    }

    return !feof(f);
}


#define IVF_FILE_HDR_SZ (32)
unsigned int file_is_ivf(FILE *infile,
                         unsigned int *fourcc,
                         unsigned int *width,
                         unsigned int *height)
{
    char raw_hdr[IVF_FILE_HDR_SZ];
    int is_ivf = 0;

    /* See write_ivf_file_header() for more documentation on the file header
     * layout.
     */
    if (fread(raw_hdr, 1, IVF_FILE_HDR_SZ, infile) == IVF_FILE_HDR_SZ)
    {
        if (raw_hdr[0] == 'D' && raw_hdr[1] == 'K'
            && raw_hdr[2] == 'I' && raw_hdr[3] == 'F')
        {
            is_ivf = 1;

            if (mem_get_le16(raw_hdr + 4) != 0)
                fprintf(stderr, "Error: Unrecognized IVF version! This file may not"
                        " decode properly.");

            *fourcc = mem_get_le32(raw_hdr + 8);
        }
    }

    if (is_ivf)
    {
        *width = mem_get_le16(raw_hdr + 12);
        *height = mem_get_le16(raw_hdr + 14);
    }
    else
        rewind(infile);

    return is_ivf;
}


static void write_ivf_file_header(FILE *outfile,
                                  const vpx_codec_enc_cfg_t *cfg,
                                  unsigned int fourcc,
                                  int frame_cnt)
{
    char header[32];

    if (cfg->g_pass != VPX_RC_ONE_PASS && cfg->g_pass != VPX_RC_LAST_PASS)
        return;

    header[0] = 'D';
    header[1] = 'K';
    header[2] = 'I';
    header[3] = 'F';
    mem_put_le16(header + 4,  0);                 /* version */
    mem_put_le16(header + 6,  32);                /* headersize */
    mem_put_le32(header + 8,  fourcc);            /* headersize */
    mem_put_le16(header + 12, cfg->g_w);          /* width */
    mem_put_le16(header + 14, cfg->g_h);          /* height */
    mem_put_le32(header + 16, cfg->g_timebase.den); /* rate */
    mem_put_le32(header + 20, cfg->g_timebase.num); /* scale */
    mem_put_le32(header + 24, frame_cnt);         /* length */
    mem_put_le32(header + 28, 0);                 /* unused */

    fwrite(header, 1, 32, outfile);
}


static void write_ivf_frame_header(FILE *outfile,
                                   const vpx_codec_cx_pkt_t *pkt)
{
    char             header[12];
    vpx_codec_pts_t  pts;

    if (pkt->kind != VPX_CODEC_CX_FRAME_PKT)
        return;

    pts = pkt->data.frame.pts;
    mem_put_le32(header, pkt->data.frame.sz);
    mem_put_le32(header + 4, pts & 0xFFFFFFFF);
    mem_put_le32(header + 8, pts >> 32);

    fwrite(header, 1, 12, outfile);
}

#include "args.h"

static const arg_def_t use_yv12 = ARG_DEF(NULL, "yv12", 0,
                                  "Input file is YV12 ");
static const arg_def_t use_i420 = ARG_DEF(NULL, "i420", 0,
                                  "Input file is I420 (default)");
static const arg_def_t codecarg = ARG_DEF(NULL, "codec", 1,
                                  "Codec to use");
static const arg_def_t passes           = ARG_DEF("p", "passes", 1,
        "Number of passes (1/2)");
static const arg_def_t pass_arg         = ARG_DEF(NULL, "pass", 1,
        "Pass to execute (1/2)");
static const arg_def_t fpf_name         = ARG_DEF(NULL, "fpf", 1,
        "First pass statistics file name");
static const arg_def_t limit = ARG_DEF(NULL, "limit", 1,
                                       "Stop encoding after n input frames");
static const arg_def_t deadline         = ARG_DEF("d", "deadline", 1,
        "Deadline per frame (usec)");
static const arg_def_t best_dl          = ARG_DEF(NULL, "best", 0,
        "Use Best Quality Deadline");
static const arg_def_t good_dl          = ARG_DEF(NULL, "good", 0,
        "Use Good Quality Deadline");
static const arg_def_t rt_dl            = ARG_DEF(NULL, "rt", 0,
        "Use Realtime Quality Deadline");
static const arg_def_t verbosearg       = ARG_DEF("v", "verbose", 0,
        "Show encoder parameters");
static const arg_def_t psnrarg          = ARG_DEF(NULL, "psnr", 0,
        "Show PSNR in status line");
static const arg_def_t *main_args[] =
{
    &codecarg, &passes, &pass_arg, &fpf_name, &limit, &deadline, &best_dl, &good_dl, &rt_dl,
    &verbosearg, &psnrarg,
    NULL
};

static const arg_def_t usage            = ARG_DEF("u", "usage", 1,
        "Usage profile number to use");
static const arg_def_t threads          = ARG_DEF("t", "threads", 1,
        "Max number of threads to use");
static const arg_def_t profile          = ARG_DEF(NULL, "profile", 1,
        "Bitstream profile number to use");
static const arg_def_t width            = ARG_DEF("w", "width", 1,
        "Frame width");
static const arg_def_t height           = ARG_DEF("h", "height", 1,
        "Frame height");
static const arg_def_t timebase         = ARG_DEF(NULL, "timebase", 1,
        "Stream timebase (frame duration)");
static const arg_def_t error_resilient  = ARG_DEF(NULL, "error-resilient", 1,
        "Enable error resiliency features");
static const arg_def_t lag_in_frames    = ARG_DEF(NULL, "lag-in-frames", 1,
        "Max number of frames to lag");

static const arg_def_t *global_args[] =
{
    &use_yv12, &use_i420, &usage, &threads, &profile,
    &width, &height, &timebase, &error_resilient,
    &lag_in_frames, NULL
};

static const arg_def_t dropframe_thresh   = ARG_DEF(NULL, "drop-frame", 1,
        "Temporal resampling threshold (buf %)");
static const arg_def_t resize_allowed     = ARG_DEF(NULL, "resize-allowed", 1,
        "Spatial resampling enabled (bool)");
static const arg_def_t resize_up_thresh   = ARG_DEF(NULL, "resize-up", 1,
        "Upscale threshold (buf %)");
static const arg_def_t resize_down_thresh = ARG_DEF(NULL, "resize-down", 1,
        "Downscale threshold (buf %)");
static const arg_def_t end_usage          = ARG_DEF(NULL, "end-usage", 1,
        "VBR=0 | CBR=1");
static const arg_def_t target_bitrate     = ARG_DEF(NULL, "target-bitrate", 1,
        "Bitrate (kbps)");
static const arg_def_t min_quantizer      = ARG_DEF(NULL, "min-q", 1,
        "Minimum (best) quantizer");
static const arg_def_t max_quantizer      = ARG_DEF(NULL, "max-q", 1,
        "Maximum (worst) quantizer");
static const arg_def_t undershoot_pct     = ARG_DEF(NULL, "undershoot-pct", 1,
        "Datarate undershoot (min) target (%)");
static const arg_def_t overshoot_pct      = ARG_DEF(NULL, "overshoot-pct", 1,
        "Datarate overshoot (max) target (%)");
static const arg_def_t buf_sz             = ARG_DEF(NULL, "buf-sz", 1,
        "Client buffer size (ms)");
static const arg_def_t buf_initial_sz     = ARG_DEF(NULL, "buf-initial-sz", 1,
        "Client initial buffer size (ms)");
static const arg_def_t buf_optimal_sz     = ARG_DEF(NULL, "buf-optimal-sz", 1,
        "Client optimal buffer size (ms)");
static const arg_def_t *rc_args[] =
{
    &dropframe_thresh, &resize_allowed, &resize_up_thresh, &resize_down_thresh,
    &end_usage, &target_bitrate, &min_quantizer, &max_quantizer,
    &undershoot_pct, &overshoot_pct, &buf_sz, &buf_initial_sz, &buf_optimal_sz,
    NULL
};


static const arg_def_t bias_pct = ARG_DEF(NULL, "bias-pct", 1,
                                  "CBR/VBR bias (0=CBR, 100=VBR)");
static const arg_def_t minsection_pct = ARG_DEF(NULL, "minsection-pct", 1,
                                        "GOP min bitrate (% of target)");
static const arg_def_t maxsection_pct = ARG_DEF(NULL, "maxsection-pct", 1,
                                        "GOP max bitrate (% of target)");
static const arg_def_t *rc_twopass_args[] =
{
    &bias_pct, &minsection_pct, &maxsection_pct, NULL
};


static const arg_def_t kf_min_dist = ARG_DEF(NULL, "kf-min-dist", 1,
                                     "Minimum keyframe interval (frames)");
static const arg_def_t kf_max_dist = ARG_DEF(NULL, "kf-max-dist", 1,
                                     "Maximum keyframe interval (frames)");
static const arg_def_t *kf_args[] =
{
    &kf_min_dist, &kf_max_dist, NULL
};


#if CONFIG_VP8_ENCODER
static const arg_def_t noise_sens = ARG_DEF(NULL, "noise-sensitivity", 1,
                                    "Noise sensitivity (frames to blur)");
static const arg_def_t sharpness = ARG_DEF(NULL, "sharpness", 1,
                                   "Filter sharpness (0-7)");
static const arg_def_t static_thresh = ARG_DEF(NULL, "static-thresh", 1,
                                       "Motion detection threshold");
#endif

#if CONFIG_VP8_ENCODER
static const arg_def_t cpu_used = ARG_DEF(NULL, "cpu-used", 1,
                                  "CPU Used (-16..16)");
#endif


#if CONFIG_VP8_ENCODER
static const arg_def_t token_parts = ARG_DEF(NULL, "token-parts", 1,
                                     "Number of token partitions to use, log2");
static const arg_def_t auto_altref = ARG_DEF(NULL, "auto-alt-ref", 1,
                                     "Enable automatic alt reference frames");
static const arg_def_t arnr_maxframes = ARG_DEF(NULL, "arnr-maxframes", 1,
                                        "alt_ref Max Frames");
static const arg_def_t arnr_strength = ARG_DEF(NULL, "arnr-strength", 1,
                                       "alt_ref Strength");
static const arg_def_t arnr_type = ARG_DEF(NULL, "arnr-type", 1,
                                   "alt_ref Type");

static const arg_def_t *vp8_args[] =
{
    &cpu_used, &auto_altref, &noise_sens, &sharpness, &static_thresh,
    &token_parts, &arnr_maxframes, &arnr_strength, &arnr_type, NULL
};
static const int vp8_arg_ctrl_map[] =
{
    VP8E_SET_CPUUSED, VP8E_SET_ENABLEAUTOALTREF,
    VP8E_SET_NOISE_SENSITIVITY, VP8E_SET_SHARPNESS, VP8E_SET_STATIC_THRESHOLD,
    VP8E_SET_TOKEN_PARTITIONS,
    VP8E_SET_ARNR_MAXFRAMES, VP8E_SET_ARNR_STRENGTH , VP8E_SET_ARNR_TYPE, 0
};
#endif

static const arg_def_t *no_args[] = { NULL };

static void usage_exit()
{
    int i;

    printf("Usage: %s <options> src_filename dst_filename\n", exec_name);

    printf("\n_options:\n");
    arg_show_usage(stdout, main_args);
    printf("\n_encoder Global Options:\n");
    arg_show_usage(stdout, global_args);
    printf("\n_rate Control Options:\n");
    arg_show_usage(stdout, rc_args);
    printf("\n_twopass Rate Control Options:\n");
    arg_show_usage(stdout, rc_twopass_args);
    printf("\n_keyframe Placement Options:\n");
    arg_show_usage(stdout, kf_args);
#if CONFIG_VP8_ENCODER
    printf("\n_vp8 Specific Options:\n");
    arg_show_usage(stdout, vp8_args);
#endif
    printf("\n"
           "Included encoders:\n"
           "\n");

    for (i = 0; i < sizeof(codecs) / sizeof(codecs[0]); i++)
        printf("    %-6s - %s\n",
               codecs[i].name,
               vpx_codec_iface_name(codecs[i].iface));

    exit(EXIT_FAILURE);
}

#define ARG_CTRL_CNT_MAX 10


int main(int argc, const char **argv_)
{
    vpx_codec_ctx_t        encoder;
    const char                  *in_fn = NULL, *out_fn = NULL, *stats_fn = NULL;
    int                    i;
    FILE                  *infile, *outfile;
    vpx_codec_enc_cfg_t    cfg;
    vpx_codec_err_t        res;
    int                    pass, one_pass_only = 0;
    stats_io_t             stats;
    vpx_image_t            raw;
    const struct codec_item  *codec = codecs;
    int                    frame_avail, got_data;

    struct arg               arg;
    char                   **argv, **argi, **argj;
    int                      arg_usage = 0, arg_passes = 1, arg_deadline = 0;
    int                      arg_ctrls[ARG_CTRL_CNT_MAX][2], arg_ctrl_cnt = 0;
    int                      arg_limit = 0;
    static const arg_def_t **ctrl_args = no_args;
    static const int        *ctrl_args_map = NULL;
    int                      verbose = 0, show_psnr = 0;
    int                      arg_use_i420 = 1;
    unsigned long            cx_time = 0;
    unsigned int             is_ivf, fourcc;

    exec_name = argv_[0];

    if (argc < 3)
        usage_exit();


    /* First parse the codec and usage values, because we want to apply other
     * parameters on top of the default configuration provided by the codec.
     */
    argv = argv_dup(argc - 1, argv_ + 1);

    for (argi = argj = argv; (*argj = *argi); argi += arg.argv_step)
    {
        arg.argv_step = 1;

        if (arg_match(&arg, &codecarg, argi))
        {
            int j, k = -1;

            for (j = 0; j < sizeof(codecs) / sizeof(codecs[0]); j++)
                if (!strcmp(codecs[j].name, arg.val))
                    k = j;

            if (k >= 0)
                codec = codecs + k;
            else
                die("Error: Unrecognized argument (%s) to --codec\n",
                    arg.val);

        }
        else if (arg_match(&arg, &passes, argi))
        {
            arg_passes = arg_parse_uint(&arg);

            if (arg_passes < 1 || arg_passes > 2)
                die("Error: Invalid number of passes (%d)\n", arg_passes);
        }
        else if (arg_match(&arg, &pass_arg, argi))
        {
            one_pass_only = arg_parse_uint(&arg);

            if (one_pass_only < 1 || one_pass_only > 2)
                die("Error: Invalid pass selected (%d)\n", one_pass_only);
        }
        else if (arg_match(&arg, &fpf_name, argi))
            stats_fn = arg.val;
        else if (arg_match(&arg, &usage, argi))
            arg_usage = arg_parse_uint(&arg);
        else if (arg_match(&arg, &deadline, argi))
            arg_deadline = arg_parse_uint(&arg);
        else if (arg_match(&arg, &best_dl, argi))
            arg_deadline = VPX_DL_BEST_QUALITY;
        else if (arg_match(&arg, &good_dl, argi))
            arg_deadline = VPX_DL_GOOD_QUALITY;
        else if (arg_match(&arg, &rt_dl, argi))
            arg_deadline = VPX_DL_REALTIME;
        else if (arg_match(&arg, &use_yv12, argi))
        {
            arg_use_i420 = 0;
        }
        else if (arg_match(&arg, &use_i420, argi))
        {
            arg_use_i420 = 1;
        }
        else if (arg_match(&arg, &verbosearg, argi))
            verbose = 1;
        else if (arg_match(&arg, &limit, argi))
            arg_limit = arg_parse_uint(&arg);
        else if (arg_match(&arg, &psnrarg, argi))
            show_psnr = 1;
        else
            argj++;
    }

    /* Ensure that --passes and --pass are consistent. If --pass is set and --passes=2,
     * ensure --fpf was set.
     */
    if (one_pass_only)
    {
        /* DWIM: Assume the user meant passes=2 if pass=2 is specified */
        if (one_pass_only > arg_passes)
        {
            printf("Warning: Assuming --pass=%d implies --passes=%d\n",
                   one_pass_only, one_pass_only);
            arg_passes = one_pass_only;
        }

        if (arg_passes == 2 && !stats_fn)
            die("Must specify --fpf when --pass=%d and --passes=2\n", one_pass_only);
    }

    /* Populate encoder configuration */
    res = vpx_codec_enc_config_default(codec->iface, &cfg, arg_usage);

    if (res)
    {
        printf("Failed to get config: %s\n", vpx_codec_err_to_string(res));
        return EXIT_FAILURE;
    }

    /* Now parse the remainder of the parameters. */
    for (argi = argj = argv; (*argj = *argi); argi += arg.argv_step)
    {
        arg.argv_step = 1;

        if (0);
        else if (arg_match(&arg, &threads, argi))
            cfg.g_threads = arg_parse_uint(&arg);
        else if (arg_match(&arg, &profile, argi))
            cfg.g_profile = arg_parse_uint(&arg);
        else if (arg_match(&arg, &width, argi))
            cfg.g_w = arg_parse_uint(&arg);
        else if (arg_match(&arg, &height, argi))
            cfg.g_h = arg_parse_uint(&arg);
        else if (arg_match(&arg, &timebase, argi))
            cfg.g_timebase = arg_parse_rational(&arg);
        else if (arg_match(&arg, &error_resilient, argi))
            cfg.g_error_resilient = arg_parse_uint(&arg);
        else if (arg_match(&arg, &lag_in_frames, argi))
            cfg.g_lag_in_frames = arg_parse_uint(&arg);
        else if (arg_match(&arg, &dropframe_thresh, argi))
            cfg.rc_dropframe_thresh = arg_parse_uint(&arg);
        else if (arg_match(&arg, &resize_allowed, argi))
            cfg.rc_resize_allowed = arg_parse_uint(&arg);
        else if (arg_match(&arg, &resize_up_thresh, argi))
            cfg.rc_resize_up_thresh = arg_parse_uint(&arg);
        else if (arg_match(&arg, &resize_down_thresh, argi))
            cfg.rc_resize_down_thresh = arg_parse_uint(&arg);
        else if (arg_match(&arg, &resize_down_thresh, argi))
            cfg.rc_resize_down_thresh = arg_parse_uint(&arg);
        else if (arg_match(&arg, &end_usage, argi))
            cfg.rc_end_usage = arg_parse_uint(&arg);
        else if (arg_match(&arg, &target_bitrate, argi))
            cfg.rc_target_bitrate = arg_parse_uint(&arg);
        else if (arg_match(&arg, &min_quantizer, argi))
            cfg.rc_min_quantizer = arg_parse_uint(&arg);
        else if (arg_match(&arg, &max_quantizer, argi))
            cfg.rc_max_quantizer = arg_parse_uint(&arg);
        else if (arg_match(&arg, &undershoot_pct, argi))
            cfg.rc_undershoot_pct = arg_parse_uint(&arg);
        else if (arg_match(&arg, &overshoot_pct, argi))
            cfg.rc_overshoot_pct = arg_parse_uint(&arg);
        else if (arg_match(&arg, &buf_sz, argi))
            cfg.rc_buf_sz = arg_parse_uint(&arg);
        else if (arg_match(&arg, &buf_initial_sz, argi))
            cfg.rc_buf_initial_sz = arg_parse_uint(&arg);
        else if (arg_match(&arg, &buf_optimal_sz, argi))
            cfg.rc_buf_optimal_sz = arg_parse_uint(&arg);
        else if (arg_match(&arg, &bias_pct, argi))
        {
            cfg.rc_2pass_vbr_bias_pct = arg_parse_uint(&arg);

            if (arg_passes < 2)
                printf("Warning: option %s ignored in one-pass mode.\n",
                       arg.name);
        }
        else if (arg_match(&arg, &minsection_pct, argi))
        {
            cfg.rc_2pass_vbr_minsection_pct = arg_parse_uint(&arg);

            if (arg_passes < 2)
                printf("Warning: option %s ignored in one-pass mode.\n",
                       arg.name);
        }
        else if (arg_match(&arg, &maxsection_pct, argi))
        {
            cfg.rc_2pass_vbr_maxsection_pct = arg_parse_uint(&arg);

            if (arg_passes < 2)
                printf("Warning: option %s ignored in one-pass mode.\n",
                       arg.name);
        }
        else if (arg_match(&arg, &kf_min_dist, argi))
            cfg.kf_min_dist = arg_parse_uint(&arg);
        else if (arg_match(&arg, &kf_max_dist, argi))
            cfg.kf_max_dist = arg_parse_uint(&arg);
        else
            argj++;
    }

    /* Handle codec specific options */
#if CONFIG_VP8_ENCODER

    if (codec->iface == &vpx_codec_vp8_cx_algo)
    {
        ctrl_args = vp8_args;
        ctrl_args_map = vp8_arg_ctrl_map;
    }

#endif

    for (argi = argj = argv; (*argj = *argi); argi += arg.argv_step)
    {
        int match = 0;

        arg.argv_step = 1;

        for (i = 0; ctrl_args[i]; i++)
        {
            if (arg_match(&arg, ctrl_args[i], argi))
            {
                match = 1;

                if (arg_ctrl_cnt < ARG_CTRL_CNT_MAX)
                {
                    arg_ctrls[arg_ctrl_cnt][0] = ctrl_args_map[i];
                    arg_ctrls[arg_ctrl_cnt][1] = arg_parse_int(&arg);
                    arg_ctrl_cnt++;
                }
            }
        }

        if (!match)
            argj++;
    }

    /* Check for unrecognized options */
    for (argi = argv; *argi; argi++)
        if (argi[0][0] == '-')
            die("Error: Unrecognized option %s\n", *argi);

    /* Handle non-option arguments */
    in_fn = argv[0];
    out_fn = argv[1];

    if (!in_fn || !out_fn)
        usage_exit();

    /* Parse certain options from the input file, if possible */
    infile = fopen(in_fn, "rb");

    if (!infile)
    {
        printf("Failed to open input file");
        return EXIT_FAILURE;
    }

    is_ivf = file_is_ivf(infile, &fourcc, &cfg.g_w, &cfg.g_h);

    if (is_ivf)
    {
        switch (fourcc)
        {
        case 0x32315659:
            arg_use_i420 = 0;
            break;
        case 0x30323449:
            arg_use_i420 = 1;
            break;
        default:
            printf("Unsupported fourcc (%08x) in IVF\n", fourcc);
            return EXIT_FAILURE;
        }
    }

    fclose(infile);


#define SHOW(field) printf("    %-28s = %d\n", #field, cfg.field)

    if (verbose)
    {
        printf("Codec: %s\n", vpx_codec_iface_name(codec->iface));
        printf("Source file: %s Format: %s\n", in_fn, arg_use_i420 ? "I420" : "YV12");
        printf("Destination file: %s\n", out_fn);
        printf("Encoder parameters:\n");

        SHOW(g_usage);
        SHOW(g_threads);
        SHOW(g_profile);
        SHOW(g_w);
        SHOW(g_h);
        SHOW(g_timebase.num);
        SHOW(g_timebase.den);
        SHOW(g_error_resilient);
        SHOW(g_pass);
        SHOW(g_lag_in_frames);
        SHOW(rc_dropframe_thresh);
        SHOW(rc_resize_allowed);
        SHOW(rc_resize_up_thresh);
        SHOW(rc_resize_down_thresh);
        SHOW(rc_end_usage);
        SHOW(rc_target_bitrate);
        SHOW(rc_min_quantizer);
        SHOW(rc_max_quantizer);
        SHOW(rc_undershoot_pct);
        SHOW(rc_overshoot_pct);
        SHOW(rc_buf_sz);
        SHOW(rc_buf_initial_sz);
        SHOW(rc_buf_optimal_sz);
        SHOW(rc_2pass_vbr_bias_pct);
        SHOW(rc_2pass_vbr_minsection_pct);
        SHOW(rc_2pass_vbr_maxsection_pct);
        SHOW(kf_mode);
        SHOW(kf_min_dist);
        SHOW(kf_max_dist);
    }

    vpx_img_alloc(&raw, arg_use_i420 ? IMG_FMT_I420 : IMG_FMT_YV12,
                  cfg.g_w, cfg.g_h, 1);

    // This was added so that ivfenc will create monotically increasing
    // timestamps.  Since we create new timestamps for alt-reference frames
    // we need to make room in the series of timestamps.  Since there can
    // only be 1 alt-ref frame ( current bitstream) multiplying by 2
    // gives us enough room.
    cfg.g_timebase.den *= 2;

    memset(&stats, 0, sizeof(stats));

    for (pass = one_pass_only ? one_pass_only - 1 : 0; pass < arg_passes; pass++)
    {
        int frames_in = 0, frames_out = 0;
        unsigned long nbytes = 0;

        infile = fopen(in_fn, "rb");

        if (!infile)
        {
            printf("Failed to open input file");
            return EXIT_FAILURE;
        }

        outfile = fopen(out_fn, "wb");

        if (!outfile)
        {
            printf("Failed to open output file");
            return EXIT_FAILURE;
        }

        if (stats_fn)
        {
            if (!stats_open_file(&stats, stats_fn, pass))
            {
                printf("Failed to open statistics store\n");
                return EXIT_FAILURE;
            }
        }
        else
        {
            if (!stats_open_mem(&stats, pass))
            {
                printf("Failed to open statistics store\n");
                return EXIT_FAILURE;
            }
        }

        cfg.g_pass = arg_passes == 2
                     ? pass ? VPX_RC_LAST_PASS : VPX_RC_FIRST_PASS
                 : VPX_RC_ONE_PASS;
#if VPX_ENCODER_ABI_VERSION > (1 + VPX_CODEC_ABI_VERSION)

        if (pass)
        {
            cfg.rc_twopass_stats_in = stats_get(&stats);
        }

#endif

        write_ivf_file_header(outfile, &cfg, codec->fourcc, 0);


        /* Construct Encoder Context */
        if (cfg.kf_min_dist == cfg.kf_max_dist)
            cfg.kf_mode = VPX_KF_FIXED;

        vpx_codec_enc_init(&encoder, codec->iface, &cfg,
                           show_psnr ? VPX_CODEC_USE_PSNR : 0);
        ctx_exit_on_error(&encoder, "Failed to initialize encoder");

        /* Note that we bypass the vpx_codec_control wrapper macro because
         * we're being clever to store the control IDs in an array. Real
         * applications will want to make use of the enumerations directly
         */
        for (i = 0; i < arg_ctrl_cnt; i++)
        {
            if (vpx_codec_control_(&encoder, arg_ctrls[i][0], arg_ctrls[i][1]))
                printf("Error: Tried to set control %d = %d\n",
                       arg_ctrls[i][0], arg_ctrls[i][1]);

            ctx_exit_on_error(&encoder, "Failed to control codec");
        }

        frame_avail = 1;
        got_data = 0;

        while (frame_avail || got_data)
        {
            vpx_codec_iter_t iter = NULL;
            const vpx_codec_cx_pkt_t *pkt;
            struct vpx_usec_timer timer;

            if (!arg_limit || frames_in < arg_limit)
            {
                frame_avail = read_frame(infile, &raw, is_ivf);

                if (frame_avail)
                    frames_in++;

                printf("\rPass %d/%d frame %4d/%-4d %7ldB \033[K", pass + 1,
                       arg_passes, frames_in, frames_out, nbytes);
            }
            else
                frame_avail = 0;

            vpx_usec_timer_start(&timer);

            // since we halved our timebase we need to double the timestamps
            // and duration we pass in.
            vpx_codec_encode(&encoder, frame_avail ? &raw : NULL, (frames_in - 1) * 2,
                             2, 0, arg_deadline);
            vpx_usec_timer_mark(&timer);
            cx_time += vpx_usec_timer_elapsed(&timer);
            ctx_exit_on_error(&encoder, "Failed to encode frame");
            got_data = 0;

            while ((pkt = vpx_codec_get_cx_data(&encoder, &iter)))
            {
                got_data = 1;
                nbytes += pkt->data.raw.sz;

                switch (pkt->kind)
                {
                case VPX_CODEC_CX_FRAME_PKT:
                    frames_out++;
                    printf(" %6luF",
                           (unsigned long)pkt->data.frame.sz);
                    write_ivf_frame_header(outfile, pkt);
                    fwrite(pkt->data.frame.buf, 1, pkt->data.frame.sz, outfile);
                    break;
                case VPX_CODEC_STATS_PKT:
                    frames_out++;
                    printf(" %6luS",
                           (unsigned long)pkt->data.twopass_stats.sz);
                    stats_write(&stats,
                                pkt->data.twopass_stats.buf,
                                pkt->data.twopass_stats.sz);
                    break;
                case VPX_CODEC_PSNR_PKT:

                    if (show_psnr)
                    {
                        int i;

                        for (i = 0; i < 4; i++)
                            printf("%.3lf ", pkt->data.psnr.psnr[i]);
                    }

                    break;
                default:
                    break;
                }
            }

            fflush(stdout);
        }

        /* this bitrate calc is simplified and relies on the fact that this
         * application uses 1/timebase for framerate.
         */
        printf("\rPass %d/%d frame %4d/%-4d %7ldB %7ldb/f %7"PRId64"b/s"
               " %7lu %s (%.2f fps)\033[K", pass + 1,
               arg_passes, frames_in, frames_out, nbytes, nbytes * 8 / frames_in,
               nbytes * 8 *(int64_t)cfg.g_timebase.den / cfg.g_timebase.num / frames_in,
               cx_time > 9999999 ? cx_time / 1000 : cx_time,
               cx_time > 9999999 ? "ms" : "us",
               (float)frames_in * 1000000.0 / (float)cx_time);

        vpx_codec_destroy(&encoder);

        fclose(infile);

        if (!fseek(outfile, 0, SEEK_SET))
            write_ivf_file_header(outfile, &cfg, codec->fourcc, frames_out);

        fclose(outfile);
        stats_close(&stats);
        printf("\n");

        if (one_pass_only)
            break;
    }

    vpx_img_free(&raw);
    free(argv);
    return EXIT_SUCCESS;
}
