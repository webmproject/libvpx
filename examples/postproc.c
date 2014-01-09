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

    (void)res;
    /* Open files */
    if(argc!=3)
        die("Usage: %s <infile> <outfile>\n", argv[0]);
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
    res = vpx_codec_dec_init(&codec, interface, NULL,
                             VPX_CODEC_USE_POSTPROC);
    if(res == VPX_CODEC_INCAPABLE) {
        printf("NOTICE: Postproc not supported by %s\n",
               vpx_codec_iface_name(interface));
        res = vpx_codec_dec_init(&codec, interface, NULL, flags);
    }
    if(res)
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

        #if CONFIG_VP9_DECODER
        if(frame_cnt%30 == 1) {
            vp8_postproc_cfg_t  pp = {0, 0, 0};

            if(vpx_codec_control(&codec, VP8_SET_POSTPROC, &pp))
                die_codec(&codec, "Failed to turn off postproc");
        } else if(frame_cnt%30 == 16) {
            vp8_postproc_cfg_t  pp = {VP8_DEBLOCK | VP8_DEMACROBLOCK | VP8_MFQE, 4, 0};

            if(vpx_codec_control(&codec, VP8_SET_POSTPROC, &pp))
                die_codec(&codec, "Failed to turn on postproc");
        };
        #endif
        /* Decode the frame with 15ms deadline */
        if(vpx_codec_decode(&codec, frame, frame_sz, NULL, 15000))
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
