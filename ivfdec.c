/*
 *  Copyright (c) 2010 The VP8 project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license and patent
 *  grant that can be found in the LICENSE file in the root of the source
 *  tree. All contributing project authors may be found in the AUTHORS
 *  file in the root of the source tree.
 */


/* This is a simple program that reads ivf files and decodes them
 * using the new interface. Decoded frames are output as YV12 raw.
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#define VPX_CODEC_DISABLE_COMPAT 1
#include "vpx_config.h"
#include "vpx_decoder.h"
#include "vpx_ports/vpx_timer.h"
#if CONFIG_VP8_DECODER
#include "vp8dx.h"
#endif
#if CONFIG_MD5
#include "md5_utils.h"
#endif

static const char *exec_name;

static const struct
{
    char const *name;
    const vpx_codec_iface_t *iface;
    unsigned int             fourcc;
    unsigned int             fourcc_mask;
} ifaces[] =
{
#if CONFIG_VP8_DECODER
    {"vp8",  &vpx_codec_vp8_dx_algo,   0x00385056, 0x00FFFFFF},
#endif
};

#include "args.h"
static const arg_def_t codecarg = ARG_DEF(NULL, "codec", 1,
                                  "Codec to use");
static const arg_def_t prefixarg = ARG_DEF("p", "prefix", 1,
                                   "Prefix to use when saving frames");
static const arg_def_t use_yv12 = ARG_DEF(NULL, "yv12", 0,
                                  "Output file is YV12 ");
static const arg_def_t use_i420 = ARG_DEF(NULL, "i420", 0,
                                  "Output file is I420 (default)");
static const arg_def_t flipuvarg = ARG_DEF(NULL, "flipuv", 0,
                                   "Synonym for --yv12");
static const arg_def_t noblitarg = ARG_DEF(NULL, "noblit", 0,
                                   "Don't process the decoded frames");
static const arg_def_t progressarg = ARG_DEF(NULL, "progress", 0,
                                     "Show progress after each frame decodes");
static const arg_def_t limitarg = ARG_DEF(NULL, "limit", 1,
                                  "Stop decoding after n frames");
static const arg_def_t postprocarg = ARG_DEF(NULL, "postproc", 0,
                                     "Postprocess decoded frames");
static const arg_def_t summaryarg = ARG_DEF(NULL, "summary", 0,
                                    "Show timing summary");
static const arg_def_t outputfile = ARG_DEF("o", "output-raw-file", 1,
                                    "Output raw yv12 file instead of images");
static const arg_def_t threadsarg = ARG_DEF("t", "threads", 1,
                                    "Max threads to use");
static const arg_def_t quietarg = ARG_DEF("q", "quiet", 0,
                                  "Suppress version string");

#if CONFIG_MD5
static const arg_def_t md5arg = ARG_DEF(NULL, "md5", 0,
                                        "Compute the MD5 sum of the decoded frame");
#endif
static const arg_def_t *all_args[] =
{
    &codecarg, &prefixarg, &use_yv12, &use_i420, &flipuvarg, &noblitarg,
    &progressarg, &limitarg, &postprocarg, &summaryarg, &outputfile,
    &threadsarg, &quietarg,
#if CONFIG_MD5
    &md5arg,
#endif
    NULL
};

#if CONFIG_VP8_DECODER
static const arg_def_t addnoise_level = ARG_DEF(NULL, "noise-level", 1,
                                        "Enable VP8 postproc add noise");
static const arg_def_t deblock = ARG_DEF(NULL, "deblock", 0,
                                 "Enable VP8 deblocking");
static const arg_def_t demacroblock_level = ARG_DEF(NULL, "demacroblock-level", 1,
        "Enable VP8 demacroblocking, w/ level");
static const arg_def_t pp_debug_info = ARG_DEF(NULL, "pp-debug-info", 1,
                                       "Enable VP8 visible debug info");


static const arg_def_t *vp8_pp_args[] =
{
    &addnoise_level, &deblock, &demacroblock_level, &pp_debug_info,
    NULL
};
#endif

static void usage_exit()
{
    int i;

    fprintf(stderr, "Usage: %s <options> filename\n\n"
            "Options:\n", exec_name);
    arg_show_usage(stderr, all_args);
#if CONFIG_VP8_DECODER
    fprintf(stderr, "\nvp8 Postprocessing Options:\n");
    arg_show_usage(stderr, vp8_pp_args);
#endif
    fprintf(stderr, "\nIncluded decoders:\n\n");

    for (i = 0; i < sizeof(ifaces) / sizeof(ifaces[0]); i++)
        fprintf(stderr, "    %-6s - %s\n",
                ifaces[i].name,
                vpx_codec_iface_name(ifaces[i].iface));

    exit(EXIT_FAILURE);
}

void die(const char *fmt, ...)
{
    va_list ap;
    va_start(ap, fmt);
    vfprintf(stderr, fmt, ap);
    fprintf(stderr, "\n");
    usage_exit();
}

static unsigned int mem_get_le16(const void *vmem)
{
    unsigned int  val;
    const unsigned char *mem = (const unsigned char *)vmem;

    val = mem[1] << 8;
    val |= mem[0];
    return val;
}

static unsigned int mem_get_le32(const void *vmem)
{
    unsigned int  val;
    const unsigned char *mem = (const unsigned char *)vmem;

    val = mem[3] << 24;
    val |= mem[2] << 16;
    val |= mem[1] << 8;
    val |= mem[0];
    return val;
}

#define IVF_FRAME_HDR_SZ (sizeof(uint32_t) + sizeof(uint64_t))
#define RAW_FRAME_HDR_SZ (sizeof(uint32_t))
static int read_frame(FILE                  *infile,
                      uint8_t               **buf,
                      uint32_t              *buf_sz,
                      uint32_t              *buf_alloc_sz,
                      int                    is_ivf)
{
    char     raw_hdr[IVF_FRAME_HDR_SZ];
    uint32_t new_buf_sz;

    /* For both the raw and ivf formats, the frame size is the first 4 bytes
     * of the frame header. We just need to special case on the header
     * size.
     */
    if (fread(raw_hdr, is_ivf ? IVF_FRAME_HDR_SZ : RAW_FRAME_HDR_SZ, 1,
              infile) != 1)
    {
        if (!feof(infile))
            fprintf(stderr, "Failed to read frame size\n");

        new_buf_sz = 0;
    }
    else
    {
        new_buf_sz = mem_get_le32(raw_hdr);

        if (new_buf_sz > 256 * 1024 * 1024)
        {
            fprintf(stderr, "Error: Read invalid frame size (%u)\n",
                    new_buf_sz);
            new_buf_sz = 0;
        }

        if (!is_ivf && new_buf_sz > 256 * 1024)
            fprintf(stderr, "Warning: Read invalid frame size (%u)"
                    " - not a raw file?\n", new_buf_sz);

        if (new_buf_sz > *buf_alloc_sz)
        {
            uint8_t *new_buf = realloc(*buf, 2 * new_buf_sz);

            if (new_buf)
            {
                *buf = new_buf;
                *buf_alloc_sz = 2 * new_buf_sz;
            }
            else
            {
                fprintf(stderr, "Failed to allocate compressed data buffer\n");
                new_buf_sz = 0;
            }
        }
    }

    *buf_sz = new_buf_sz;

    if (*buf_sz)
    {
        if (fread(*buf, 1, *buf_sz, infile) != *buf_sz)
        {
            fprintf(stderr, "Failed to read full frame\n");
            return 1;
        }

        return 0;
    }

    return 1;
}

void *out_open(const char *out_fn, int do_md5)
{
    void *out = NULL;

    if (do_md5)
    {
#if CONFIG_MD5
        md5_ctx_t *md5_ctx = out = malloc(sizeof(md5_ctx_t));
        (void)out_fn;
        md5_init(md5_ctx);
#endif
    }
    else
    {
        FILE *outfile = out = strcmp("-", out_fn) ? fopen(out_fn, "wb") : stdout;

        if (!outfile)
        {
            fprintf(stderr, "Failed to output file");
            exit(EXIT_FAILURE);
        }
    }

    return out;
}

void out_put(void *out, const uint8_t *buf, unsigned int len, int do_md5)
{
    if (do_md5)
    {
#if CONFIG_MD5
        md5_update(out, buf, len);
#endif
    }
    else
    {
        fwrite(buf, 1, len, out);
    }
}

void out_close(void *out, const char *out_fn, int do_md5)
{
    if (do_md5)
    {
#if CONFIG_MD5
        uint8_t md5[16];
        int i;

        md5_finalize(out, md5);
        free(out);

        for (i = 0; i < 16; i++)
            printf("%02x", md5[i]);

        printf("  %s\n", out_fn);
#endif
    }
    else
    {
        fclose(out);
    }
}

unsigned int file_is_ivf(FILE *infile, unsigned int *fourcc)
{
    char raw_hdr[32];
    int is_ivf = 0;

    if (fread(raw_hdr, 1, 32, infile) == 32)
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

    if (!is_ivf)
        rewind(infile);

    return is_ivf;
}

int main(int argc, const char **argv_)
{
    vpx_codec_ctx_t          decoder;
    char                  *prefix = NULL, *fn = NULL;
    int                    i;
    uint8_t               *buf = NULL;
    uint32_t               buf_sz = 0, buf_alloc_sz = 0;
    FILE                  *infile;
    int                    frame_in = 0, frame_out = 0, flipuv = 0, noblit = 0, do_md5 = 0, progress = 0;
    int                    stop_after = 0, postproc = 0, summary = 0, quiet = 0;
    vpx_codec_iface_t       *iface = NULL;
    unsigned int           is_ivf, fourcc;
    unsigned long          dx_time = 0;
    struct arg               arg;
    char                   **argv, **argi, **argj;
    const char                   *fn2 = 0;
    void                   *out = NULL;
    vpx_codec_dec_cfg_t     cfg = {0};
#if CONFIG_VP8_DECODER
    vp8_postproc_cfg_t      vp8_pp_cfg = {0};
#endif

    /* Parse command line */
    exec_name = argv_[0];
    argv = argv_dup(argc - 1, argv_ + 1);

    for (argi = argj = argv; (*argj = *argi); argi += arg.argv_step)
    {
        memset(&arg, 0, sizeof(arg));
        arg.argv_step = 1;

        if (arg_match(&arg, &codecarg, argi))
        {
            int j, k = -1;

            for (j = 0; j < sizeof(ifaces) / sizeof(ifaces[0]); j++)
                if (!strcmp(ifaces[j].name, arg.val))
                    k = j;

            if (k >= 0)
                iface = ifaces[k].iface;
            else
                die("Error: Unrecognized argument (%s) to --codec\n",
                    arg.val);
        }
        else if (arg_match(&arg, &outputfile, argi))
            fn2 = arg.val;
        else if (arg_match(&arg, &prefixarg, argi))
            prefix = strdup(arg.val);
        else if (arg_match(&arg, &use_yv12, argi))
            flipuv = 1;
        else if (arg_match(&arg, &use_i420, argi))
            flipuv = 0;
        else if (arg_match(&arg, &flipuvarg, argi))
            flipuv = 1;
        else if (arg_match(&arg, &noblitarg, argi))
            noblit = 1;
        else if (arg_match(&arg, &progressarg, argi))
            progress = 1;
        else if (arg_match(&arg, &limitarg, argi))
            stop_after = arg_parse_uint(&arg);
        else if (arg_match(&arg, &postprocarg, argi))
            postproc = 1;
        else if (arg_match(&arg, &md5arg, argi))
            do_md5 = 1;
        else if (arg_match(&arg, &summaryarg, argi))
            summary = 1;
        else if (arg_match(&arg, &threadsarg, argi))
            cfg.threads = arg_parse_uint(&arg);
        else if (arg_match(&arg, &quietarg, argi))
            quiet = 1;

#if CONFIG_VP8_DECODER
        else if (arg_match(&arg, &addnoise_level, argi))
        {
            postproc = 1;
            vp8_pp_cfg.post_proc_flag |= VP8_ADDNOISE;
            vp8_pp_cfg.noise_level = arg_parse_uint(&arg);
        }
        else if (arg_match(&arg, &demacroblock_level, argi))
        {
            postproc = 1;
            vp8_pp_cfg.post_proc_flag |= VP8_DEMACROBLOCK;
            vp8_pp_cfg.deblocking_level = arg_parse_uint(&arg);
        }
        else if (arg_match(&arg, &deblock, argi))
        {
            postproc = 1;
            vp8_pp_cfg.post_proc_flag |= VP8_DEBLOCK;
        }
        else if (arg_match(&arg, &pp_debug_info, argi))
        {
            unsigned int level = arg_parse_uint(&arg);

            postproc = 1;
            vp8_pp_cfg.post_proc_flag &= ~0x7;

            if (level)
                vp8_pp_cfg.post_proc_flag |= 8 << (level - 1);
        }

#endif
        else
            argj++;
    }

    /* Check for unrecognized options */
    for (argi = argv; *argi; argi++)
        if (argi[0][0] == '-' && strlen(argi[0]) > 1)
            die("Error: Unrecognized option %s\n", *argi);

    /* Handle non-option arguments */
    fn = argv[0];

    if (!fn)
        usage_exit();

    if (!prefix)
        prefix = strdup("img");

    /* Open file */
    infile = strcmp(fn, "-") ? fopen(fn, "rb") : stdin;

    if (!infile)
    {
        fprintf(stderr, "Failed to open file");
        return EXIT_FAILURE;
    }

    if (fn2)
        out = out_open(fn2, do_md5);

    is_ivf = file_is_ivf(infile, &fourcc);

    if (is_ivf)
    {
        /* Try to determine the codec from the fourcc. */
        for (i = 0; i < sizeof(ifaces) / sizeof(ifaces[0]); i++)
            if ((fourcc & ifaces[i].fourcc_mask) == ifaces[i].fourcc)
            {
                vpx_codec_iface_t  *ivf_iface = ifaces[i].iface;

                if (iface && iface != ivf_iface)
                    fprintf(stderr, "Notice -- IVF header indicates codec: %s\n",
                            ifaces[i].name);
                else
                    iface = ivf_iface;

                break;
            }
    }

    if (vpx_codec_dec_init(&decoder, iface ? iface :  ifaces[0].iface, &cfg,
                           postproc ? VPX_CODEC_USE_POSTPROC : 0))
    {
        fprintf(stderr, "Failed to initialize decoder: %s\n", vpx_codec_error(&decoder));
        return EXIT_FAILURE;
    }

    if (!quiet)
        fprintf(stderr, "%s\n", decoder.name);

#if CONFIG_VP8_DECODER

    if (vp8_pp_cfg.post_proc_flag
        && vpx_codec_control(&decoder, VP8_SET_POSTPROC, &vp8_pp_cfg))
    {
        fprintf(stderr, "Failed to configure postproc: %s\n", vpx_codec_error(&decoder));
        return EXIT_FAILURE;
    }

#endif

    /* Decode file */
    while (!read_frame(infile, &buf, &buf_sz, &buf_alloc_sz, is_ivf))
    {
        vpx_codec_iter_t  iter = NULL;
        vpx_image_t    *img;
        struct vpx_usec_timer timer;

        vpx_usec_timer_start(&timer);

        if (vpx_codec_decode(&decoder, buf, buf_sz, NULL, 0))
        {
            const char *detail = vpx_codec_error_detail(&decoder);
            fprintf(stderr, "Failed to decode frame: %s\n", vpx_codec_error(&decoder));

            if (detail)
                fprintf(stderr, "  Additional information: %s\n", detail);

            goto fail;
        }

        vpx_usec_timer_mark(&timer);
        dx_time += vpx_usec_timer_elapsed(&timer);

        ++frame_in;

        if (progress)
            fprintf(stderr, "decoded frame %d.\n", frame_in);

        if ((img = vpx_codec_get_frame(&decoder, &iter)))
            ++frame_out;

        if (!noblit)
        {
            if (img)
            {
                unsigned int y;
                char out_fn[128+24];
                uint8_t *buf;
                const char *sfx = flipuv ? "yv12" : "i420";

                if (!fn2)
                {
                    sprintf(out_fn, "%s-%dx%d-%04d.%s",
                            prefix, img->d_w, img->d_h, frame_in, sfx);
                    out = out_open(out_fn, do_md5);
                }

                buf = img->planes[PLANE_Y];

                for (y = 0; y < img->d_h; y++)
                {
                    out_put(out, buf, img->d_w, do_md5);
                    buf += img->stride[PLANE_Y];
                }

                buf = img->planes[flipuv?PLANE_V:PLANE_U];

                for (y = 0; y < (1 + img->d_h) / 2; y++)
                {
                    out_put(out, buf, (1 + img->d_w) / 2, do_md5);
                    buf += img->stride[PLANE_U];
                }

                buf = img->planes[flipuv?PLANE_U:PLANE_V];

                for (y = 0; y < (1 + img->d_h) / 2; y++)
                {
                    out_put(out, buf, (1 + img->d_w) / 2, do_md5);
                    buf += img->stride[PLANE_V];
                }

                if (!fn2)
                    out_close(out, out_fn, do_md5);
            }
        }

        if (stop_after && frame_in >= stop_after)
            break;
    }

    if (summary)
    {
        fprintf(stderr, "%d decoded frames/%d showed frames in %lu us (%.2f fps)\n",
                frame_in, frame_out, dx_time, (float)frame_out * 1000000.0 / (float)dx_time);
    }

fail:

    if (vpx_codec_destroy(&decoder))
    {
        fprintf(stderr, "Failed to destroy decoder: %s\n", vpx_codec_error(&decoder));
        return EXIT_FAILURE;
    }

    if (fn2)
        out_close(out, fn2, do_md5);

    free(buf);
    fclose(infile);
    free(prefix);
    free(argv);

    return EXIT_SUCCESS;
}
