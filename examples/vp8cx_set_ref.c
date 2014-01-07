/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */


// VP8 Set Reference Frame
// =======================
//
// This is an example demonstrating how to overwrite the VP8 encoder's
// internal reference frame. In the sample we set the last frame to the
// current frame. If this is done at a cut scene it will avoid a keyframe.
// This technique could be used to bounce between two cameras.
//
// Note that the decoder would also have to set the reference frame to the
// same value on the same frame, or the video will become corrupt.
//
// Usage
// -----
// This example adds a single argument to the `simple_encoder` example,
// which specifies the frame number to update the reference frame on.
// The parameter is parsed as follows:
//
//
// Extra Variables
// ---------------
// This example maintains the frame number passed on the command line
// in the `update_frame_num` variable.
//
//
// Configuration
// -------------
//
// The reference frame is updated on the frame specified on the command
// line.
//
// Observing The Effects
// ---------------------
// Use the `simple_encoder` example to encode a sample with a cut scene.
// Determine the frame number of the cut scene by looking for a generated
// key-frame (indicated by a 'K'). Supply that frame number as an argument
// to this example, and observe that no key-frame is generated.

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#define VPX_CODEC_DISABLE_COMPAT 1
#include "vpx/vpx_encoder.h"
#include "vpx/vp8cx.h"
#define interface (vpx_codec_vp8_cx())
#define fourcc    0x30385056

#define IVF_FILE_HDR_SZ  (32)
#define IVF_FRAME_HDR_SZ (12)

static void mem_put_le16(char *mem, unsigned int val) {
    mem[0] = val;
    mem[1] = val>>8;
}

static void mem_put_le32(char *mem, unsigned int val) {
    mem[0] = val;
    mem[1] = val>>8;
    mem[2] = val>>16;
    mem[3] = val>>24;
}

static void die(const char *fmt, ...) {
    va_list ap;

    va_start(ap, fmt);
    vprintf(fmt, ap);
    if(fmt[strlen(fmt)-1] != '\n')
        printf("\n");
    exit(EXIT_FAILURE);
}

static void die_codec(vpx_codec_ctx_t *ctx, const char *s) {
    const char *detail = vpx_codec_error_detail(ctx);

    printf("%s: %s\n", s, vpx_codec_error(ctx));
    if(detail)
        printf("    %s\n",detail);
    exit(EXIT_FAILURE);
}

static int read_frame(FILE *f, vpx_image_t *img) {
    size_t nbytes, to_read;
    int    res = 1;

    to_read = img->w*img->h*3/2;
    nbytes = fread(img->planes[0], 1, to_read, f);
    if(nbytes != to_read) {
        res = 0;
        if(nbytes > 0)
            printf("Warning: Read partial frame. Check your width & height!\n");
    }
    return res;
}

static void write_ivf_file_header(FILE *outfile,
                                  const vpx_codec_enc_cfg_t *cfg,
                                  int frame_cnt) {
    char header[32];

    if(cfg->g_pass != VPX_RC_ONE_PASS && cfg->g_pass != VPX_RC_LAST_PASS)
        return;
    header[0] = 'D';
    header[1] = 'K';
    header[2] = 'I';
    header[3] = 'F';
    mem_put_le16(header+4,  0);                   /* version */
    mem_put_le16(header+6,  32);                  /* headersize */
    mem_put_le32(header+8,  fourcc);              /* headersize */
    mem_put_le16(header+12, cfg->g_w);            /* width */
    mem_put_le16(header+14, cfg->g_h);            /* height */
    mem_put_le32(header+16, cfg->g_timebase.den); /* rate */
    mem_put_le32(header+20, cfg->g_timebase.num); /* scale */
    mem_put_le32(header+24, frame_cnt);           /* length */
    mem_put_le32(header+28, 0);                   /* unused */

    (void) fwrite(header, 1, 32, outfile);
}


static void write_ivf_frame_header(FILE *outfile,
                                   const vpx_codec_cx_pkt_t *pkt)
{
    char             header[12];
    vpx_codec_pts_t  pts;

    if(pkt->kind != VPX_CODEC_CX_FRAME_PKT)
        return;

    pts = pkt->data.frame.pts;
    mem_put_le32(header, pkt->data.frame.sz);
    mem_put_le32(header+4, pts&0xFFFFFFFF);
    mem_put_le32(header+8, pts >> 32);

    (void) fwrite(header, 1, 12, outfile);
}

int main(int argc, char **argv) {
    FILE                *infile, *outfile;
    vpx_codec_ctx_t      codec;
    vpx_codec_enc_cfg_t  cfg;
    int                  frame_cnt = 0;
    vpx_image_t          raw;
    vpx_codec_err_t      res;
    long                 width;
    long                 height;
    int                  frame_avail;
    int                  got_data;
    int                  flags = 0;
    int                  update_frame_num = 0;

    /* Open files */
    if(argc!=6)
        die("Usage: %s <width> <height> <infile> <outfile> <frame>\n",
            argv[0]);

        update_frame_num = atoi(argv[5]);
        if(!update_frame_num)
            die("Couldn't parse frame number '%s'\n", argv[5]);

    width = strtol(argv[1], NULL, 0);
    height = strtol(argv[2], NULL, 0);
    if(width < 16 || width%2 || height <16 || height%2)
        die("Invalid resolution: %ldx%ld", width, height);
    if(!vpx_img_alloc(&raw, VPX_IMG_FMT_I420, width, height, 1))
        die("Faile to allocate image", width, height);
    if(!(outfile = fopen(argv[4], "wb")))
        die("Failed to open %s for writing", argv[4]);

    printf("Using %s\n",vpx_codec_iface_name(interface));

    /* Populate encoder configuration */
    res = vpx_codec_enc_config_default(interface, &cfg, 0);
    if(res) {
        printf("Failed to get config: %s\n", vpx_codec_err_to_string(res));
        return EXIT_FAILURE;
    }

    /* Update the default configuration with our settings */
    cfg.rc_target_bitrate = width * height * cfg.rc_target_bitrate
                            / cfg.g_w / cfg.g_h;
    cfg.g_w = width;
    cfg.g_h = height;

    write_ivf_file_header(outfile, &cfg, 0);


        /* Open input file for this encoding pass */
        if(!(infile = fopen(argv[3], "rb")))
            die("Failed to open %s for reading", argv[3]);

        /* Initialize codec */
        if(vpx_codec_enc_init(&codec, interface, &cfg, 0))
            die_codec(&codec, "Failed to initialize encoder");

        frame_avail = 1;
        got_data = 0;
        while(frame_avail || got_data) {
            vpx_codec_iter_t iter = NULL;
            const vpx_codec_cx_pkt_t *pkt;

            frame_avail = read_frame(infile, &raw);

            if(frame_cnt + 1 == update_frame_num) {
                vpx_ref_frame_t ref;

                ref.frame_type = VP8_LAST_FRAME;
                ref.img        = raw;

                if(vpx_codec_control(&codec, VP8_SET_REFERENCE, &ref))
                    die_codec(&codec, "Failed to set reference frame");
            }

            if(vpx_codec_encode(&codec, frame_avail? &raw : NULL, frame_cnt,
                                1, flags, VPX_DL_REALTIME))
                die_codec(&codec, "Failed to encode frame");
            got_data = 0;
            while( (pkt = vpx_codec_get_cx_data(&codec, &iter)) ) {
                got_data = 1;
                switch(pkt->kind) {
                case VPX_CODEC_CX_FRAME_PKT:
                    write_ivf_frame_header(outfile, pkt);
                    (void) fwrite(pkt->data.frame.buf, 1, pkt->data.frame.sz,
                                  outfile);
                    break;
                default:
                    break;
                }
                printf(pkt->kind == VPX_CODEC_CX_FRAME_PKT
                       && (pkt->data.frame.flags & VPX_FRAME_IS_KEY)? "K":".");
                fflush(stdout);
            }
            frame_cnt++;
        }
        printf("\n");
        fclose(infile);

    printf("Processed %d frames.\n",frame_cnt-1);
    vpx_img_free(&raw);
    if(vpx_codec_destroy(&codec))
        die_codec(&codec, "Failed to destroy codec");

    /* Try to rewrite the file header with the actual frame count */
    if(!fseek(outfile, 0, SEEK_SET))
        write_ivf_file_header(outfile, &cfg, frame_cnt-1);
    fclose(outfile);
    return EXIT_SUCCESS;
}
