/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

// Frame-by-frame MD5 Checksum
// ===========================
//
// This example builds upon the simple decoder loop to show how checksums
// of the decoded output can be generated. These are used for validating
// decoder implementations against the reference implementation, for example.
//
// MD5 algorithm
// -------------
// The Message-Digest 5 (MD5) is a well known hash function. We have provided
// an implementation derived from the RSA Data Security, Inc. MD5 Message-Digest
// Algorithm for your use. Our implmentation only changes the interface of this
// reference code. You must include the `md5_utils.h` header for access to these
// functions.
//
// Processing The Decoded Data
// ---------------------------
// Each row of the image is passed to the MD5 accumulator. First the Y plane
// is processed, then U, then V. It is important to honor the image's `stride`
// values.

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define VPX_CODEC_DISABLE_COMPAT 1
#include "./vpx_config.h"
#include "vpx/vp8dx.h"
#include "vpx/vpx_decoder.h"
#define interface (vpx_codec_vp8_dx())
#include "md5_utils.h"

#define IVF_FILE_HDR_SZ  (32)
#define IVF_FRAME_HDR_SZ (12)

static unsigned int mem_get_le32(const unsigned char *mem) {
  return (mem[3] << 24) | (mem[2] << 16) | (mem[1] << 8) | (mem[0]);
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

    /* Decode the frame */
    if(vpx_codec_decode(&codec, frame, frame_sz, NULL, 0))
      die_codec(&codec, "Failed to decode frame");

    /* Write decoded data to disk */
    while((img = vpx_codec_get_frame(&codec, &iter))) {
      unsigned int plane, y;

      unsigned char  md5_sum[16];
      MD5Context     md5;
      int            i;

      MD5Init(&md5);

      for(plane=0; plane < 3; plane++) {
        unsigned char *buf =img->planes[plane];

        for (y=0; y < (plane ? (img->d_h + 1) >> 1 : img->d_h); y++) {
          MD5Update(&md5, buf, (plane ? (img->d_w + 1) >> 1 : img->d_w));
          buf += img->stride[plane];
        }
      }

      MD5Final(md5_sum, &md5);
      for (i = 0; i < 16; i++)
        fprintf(outfile, "%02x",md5_sum[i]);
      fprintf(outfile, "  img-%dx%d-%04d.i420\n", img->d_w, img->d_h,
              frame_cnt);
    }
  }

  printf("Processed %d frames.\n",frame_cnt);
  if(vpx_codec_destroy(&codec))
    die_codec(&codec, "Failed to destroy codec");

  fclose(outfile);
  fclose(infile);
  return EXIT_SUCCESS;
}
