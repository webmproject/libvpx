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

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define VPX_CODEC_DISABLE_COMPAT 1
#include "./vpx_config.h"
#include "vpx/vp8dx.h"
#include "vpx/vpx_decoder.h"
#define interface (vpx_codec_vp8_dx())


#define IVF_FILE_HDR_SZ  (32)
#define IVF_FRAME_HDR_SZ (12)

static unsigned int mem_get_le32(const unsigned char *mem) {
    return (mem[3] << 24)|(mem[2] << 16)|(mem[1] << 8)|(mem[0]);
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


int main(int argc, char **argv) {
    FILE            *infile, *outfile;
    vpx_codec_ctx_t  codec;
    int              flags = 0, frame_cnt = 0;
    unsigned char    file_hdr[IVF_FILE_HDR_SZ];
    unsigned char    frame_hdr[IVF_FRAME_HDR_SZ];
    unsigned char    frame[256*1024];
    vpx_codec_err_t  res;
    int              n, m, is_range;

    (void)res;
    /* Open files */
    if(argc!=4)
        die("Usage: %s <infile> <outfile> <N-M|N/M>\n", argv[0]);
    {
        char *nptr;
        n = strtol(argv[3], &nptr, 0);
        m = strtol(nptr+1, NULL, 0);
        is_range = *nptr == '-';
        if(!n || !m || (*nptr != '-' && *nptr != '/'))
            die("Couldn't parse pattern %s\n", argv[3]);
    }
    if(!(infile = fopen(argv[1], "rb")))
        die("Failed to open %s for reading", argv[1]);
    if(!(outfile = fopen(argv[2], "wb")))
        die("Failed to open %s for writing", argv[2]);

    /* Read file header */
    if(!(fread(file_hdr, 1, IVF_FILE_HDR_SZ, infile) == IVF_FILE_HDR_SZ
         && file_hdr[0]=='D' && file_hdr[1]=='K' && file_hdr[2]=='I'
         && file_hdr[3]=='F'))
        die("%s is not an IVF file.", argv[1]);

    printf("Using %s\n",vpx_codec_iface_name(interface));
    /* Initialize codec */
    if(vpx_codec_dec_init(&codec, interface, NULL, flags))
        die_codec(&codec, "Failed to initialize decoder");

    /* Read each frame */
    while(fread(frame_hdr, 1, IVF_FRAME_HDR_SZ, infile) == IVF_FRAME_HDR_SZ) {
        int               frame_sz = mem_get_le32(frame_hdr);
        vpx_codec_iter_t  iter = NULL;
        vpx_image_t      *img;


        frame_cnt++;
        if(frame_sz > sizeof(frame))
            die("Frame %d data too big for example code buffer", frame_sz);
        if(fread(frame, 1, frame_sz, infile) != frame_sz)
            die("Frame %d failed to read complete frame", frame_cnt);

        if((is_range && frame_cnt >= n && frame_cnt <= m)
           ||(!is_range && m - (frame_cnt-1)%m <= n)) {
           putc('X', stdout);
           continue;
        }
        putc('.', stdout);
        fflush(stdout);
        /* Decode the frame */
        if(vpx_codec_decode(&codec, frame, frame_sz, NULL, 0))
            die_codec(&codec, "Failed to decode frame");

        /* Write decoded data to disk */
        while((img = vpx_codec_get_frame(&codec, &iter))) {
            unsigned int plane, y;

            for(plane=0; plane < 3; plane++) {
                unsigned char *buf =img->planes[plane];
            
                for(y=0; y < (plane ? (img->d_h + 1) >> 1 : img->d_h); y++) {
                    (void) fwrite(buf, 1, (plane ? (img->d_w + 1) >> 1 : img->d_w),
                                  outfile);
                    buf += img->stride[plane];
                }
            }
        }
    }
    printf("Processed %d frames.\n",frame_cnt);
    if(vpx_codec_destroy(&codec))
        die_codec(&codec, "Failed to destroy codec");

    fclose(outfile);
    fclose(infile);
    return EXIT_SUCCESS;
}
