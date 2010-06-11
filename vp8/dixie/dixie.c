/*
 *  Copyright (c) 2010 The VP8 project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */
#include "vpx/internal/vpx_codec_internal.h"
#include "bit_ops.h"
#include "dixie.h"


void
decode_frame(struct vp8_decoder_ctx *ctx,
             const unsigned char    *data,
             unsigned int            sz)
{
    vpx_codec_err_t  res;

    if ((res = vp8_parse_frame_header(data, sz, &ctx->frame_hdr)))
        vpx_internal_error(&ctx->error, res, "Failed to parse frame header");

    if (ctx->frame_hdr.is_experimental)
        vpx_internal_error(&ctx->error, VPX_CODEC_UNSUP_BITSTREAM,
                           "Experimental bitstreams not supported.");

    if (ctx->frame_hdr.version != 0)
        vpx_internal_error(&ctx->error, VPX_CODEC_UNSUP_BITSTREAM,
                           "Unsupported version %d", ctx->frame_hdr.version);

    if (ctx->frame_hdr.is_keyframe)
        if (ctx->frame_hdr.kf.scale_w || ctx->frame_hdr.kf.scale_h)
            vpx_internal_error(&ctx->error, VPX_CODEC_UNSUP_BITSTREAM,
                               "Spatial resampling not supported.");



}


vpx_codec_err_t
vp8_parse_frame_header(const unsigned char   *data,
                       unsigned int           sz,
                       struct vp8_frame_hdr  *hdr)
{
    unsigned long raw;

    if (sz < 10)
        return VPX_CODEC_CORRUPT_FRAME;

    /* The frame header is defined as a three byte little endian value */
    raw = data[0] | (data[1] << 8) | (data[2] << 3);
    hdr->is_keyframe     = !BITS_GET(raw, 0, 1);
    hdr->version         = BITS_GET(raw, 1, 2);
    hdr->is_experimental = BITS_GET(raw, 3, 1);
    hdr->is_shown        = BITS_GET(raw, 4, 1);
    hdr->part0_sz        = BITS_GET(raw, 5, 19);

    if (sz <= hdr->part0_sz + 10)
        return VPX_CODEC_CORRUPT_FRAME;

    if (hdr->is_keyframe)
    {
        /* Keyframe header consists of a three byte sync code followed by the
         * width and height and associated scaling factors.
         */
        if (data[3] != 0x9d || data[4] != 0x01 || data[5] != 0x2a)
            return VPX_CODEC_UNSUP_BITSTREAM;

        raw = data[6] | (data[7] << 8) | (data[8] << 16) | (data[9] << 24);
        hdr->kf.w       = BITS_GET(raw,  0, 14);
        hdr->kf.scale_w = BITS_GET(raw, 14, 2);
        hdr->kf.h       = BITS_GET(raw, 16, 14);
        hdr->kf.scale_h = BITS_GET(raw, 30, 2);

        if (!hdr->kf.w || !hdr->kf.h)
            return VPX_CODEC_UNSUP_BITSTREAM;
    }

    return VPX_CODEC_OK;
}


vpx_codec_err_t
vp8_dixie_decode_frame(struct vp8_decoder_ctx *ctx,
                       const unsigned char    *data,
                       unsigned int            sz)
{
    volatile struct vp8_decoder_ctx *ctx_ = ctx;

    ctx->error.error_code = VPX_CODEC_OK;
    ctx->error.has_detail = 0;

    if (!setjmp(ctx->error.jmp))
        decode_frame(ctx, data, sz);

    return ctx_->error.error_code;
}
