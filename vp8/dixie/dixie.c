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
#include <string.h>

enum
{
    FRAME_HEADER_SZ = 3,
    KEYFRAME_HEADER_SZ = 7
};


static void
decode_reference_header(struct vp8_decoder_ctx    *ctx,
                        struct bool_decoder       *bool,
                        struct vp8_reference_hdr  *hdr)
{
    unsigned int key = ctx->frame_hdr.is_keyframe;

    hdr->refresh_gf    = key ? 1 : bool_get_bit(bool);
    hdr->refresh_arf   = key ? 1 : bool_get_bit(bool);
    hdr->copy_gf       = key ? 0 : !hdr->refresh_gf
                         ? bool_get_uint(bool, 2) : 0;
    hdr->copy_arf      = key ? 0 : !hdr->refresh_arf
                         ? bool_get_uint(bool, 2) : 0;
    hdr->sign_bias_gf  = key ? 0 : bool_get_bit(bool);
    hdr->sign_bias_arf = key ? 0 : bool_get_bit(bool);
    hdr->refresh_entropy = bool_get_bit(bool);
    hdr->refresh_last  = key ? 1 : bool_get_bit(bool);
}


static void
decode_quantizer_header(struct vp8_decoder_ctx    *ctx,
                        struct bool_decoder       *bool,
                        struct vp8_quant_hdr      *hdr)
{
    int update;
    int last_q = hdr->q_index;

    hdr->q_index = bool_get_uint(bool, 7);
    update = last_q != hdr->q_index;
    update |= (hdr->y1_dc_delta_q = bool_maybe_get_int(bool, 4));
    update |= (hdr->y2_dc_delta_q = bool_maybe_get_int(bool, 4));
    update |= (hdr->y2_ac_delta_q = bool_maybe_get_int(bool, 4));
    update |= (hdr->y1_dc_delta_q = bool_maybe_get_int(bool, 4));
    update |= (hdr->y1_dc_delta_q = bool_maybe_get_int(bool, 4));
}


static void
decode_and_init_token_partitions(struct vp8_decoder_ctx    *ctx,
                                 struct bool_decoder       *bool,
                                 const unsigned char       *data,
                                 unsigned int               sz,
                                 struct vp8_token_hdr      *hdr)
{
    int i;

    hdr->partitions = 1 << bool_get_uint(bool, 2);

    if (sz < 3 *(hdr->partitions - 1))
        vpx_internal_error(&ctx->error, VPX_CODEC_CORRUPT_FRAME,
                           "Truncated packet found parsing partition lengths.");

    sz -= 3 * (hdr->partitions - 1);

    for (i = 0; i < hdr->partitions; i++)
    {
        if (i < hdr->partitions - 1)
            hdr->partition_sz[i] = (data[2] << 16) | (data[1] << 8) | data[0];
        else
            hdr->partition_sz[i] = sz;

        if (sz < hdr->partition_sz[i])
            vpx_internal_error(&ctx->error, VPX_CODEC_CORRUPT_FRAME,
                               "Truncated partition %d", i);

        data += 3;
        sz -= hdr->partition_sz[i];
    }


    for (i = 0; i < ctx->token_hdr.partitions; i++)
    {
        vp8dx_bool_init(&ctx->bool[i], data, ctx->token_hdr.partition_sz[i]);
        data += ctx->token_hdr.partition_sz[i];
    }
}


static void
decode_loopfilter_header(struct vp8_decoder_ctx    *ctx,
                         struct bool_decoder       *bool,
                         struct vp8_loopfilter_hdr *hdr)
{
    if (ctx->frame_hdr.is_keyframe)
        memset(hdr, 0, sizeof(*hdr));

    hdr->use_simple    = bool_get_bit(bool);
    hdr->level         = bool_get_uint(bool, 6);
    hdr->sharpness     = bool_get_uint(bool, 3);
    hdr->delta_enabled = bool_get_bit(bool);

    if (hdr->delta_enabled && bool_get_bit(bool))
    {
        int i;

        for (i = 0; i < BLOCK_CONTEXTS; i++)
            hdr->ref_delta[i] = bool_maybe_get_int(bool, 6);

        for (i = 0; i < BLOCK_CONTEXTS; i++)
            hdr->mode_delta[i] = bool_maybe_get_int(bool, 6);
    }
}


static void
decode_segmentation_header(struct vp8_decoder_ctx *ctx,
                           struct bool_decoder    *bool,
                           struct vp8_segment_hdr *hdr)
{
    if (ctx->frame_hdr.is_keyframe)
        memset(hdr, 0, sizeof(*hdr));

    hdr->enabled = bool_get_bit(bool);

    if (hdr->enabled)
    {
        int i;

        hdr->update_map = bool_get_bit(bool);
        hdr->update_data = bool_get_bit(bool);

        if (hdr->update_data)
        {
            hdr->abs = bool_get_bit(bool);

            for (i = 0; i < MAX_MB_SEGMENTS; i++)
                hdr->quant_idx[i] = bool_maybe_get_int(bool, 7);

            for (i = 0; i < MAX_MB_SEGMENTS; i++)
                hdr->lf_level[i] = bool_maybe_get_int(bool, 6);
        }

        if (hdr->update_map)
        {
            for (i = 0; i < MB_FEATURE_TREE_PROBS; i++)
                hdr->tree_probs[i] = bool_get_bit(bool)
                                     ? bool_get_uint(bool, 8)
                                     : 255;
        }
    }
    else
    {
        hdr->update_map = 0;
        hdr->update_data = 0;
    }
}


static void
decode_frame(struct vp8_decoder_ctx *ctx,
             const unsigned char    *data,
             unsigned int            sz)
{
    vpx_codec_err_t  res;
    struct bool_decoder  bool;
    int                  i;

    if ((res = vp8_parse_frame_header(data, sz, &ctx->frame_hdr)))
        vpx_internal_error(&ctx->error, res, "Failed to parse frame header");

    if (ctx->frame_hdr.is_experimental)
        vpx_internal_error(&ctx->error, VPX_CODEC_UNSUP_BITSTREAM,
                           "Experimental bitstreams not supported.");

    if (ctx->frame_hdr.version != 0)
        vpx_internal_error(&ctx->error, VPX_CODEC_UNSUP_BITSTREAM,
                           "Unsupported version %d", ctx->frame_hdr.version);

    data += FRAME_HEADER_SZ;
    sz -= FRAME_HEADER_SZ;

    if (ctx->frame_hdr.is_keyframe)
    {
        if (ctx->frame_hdr.kf.scale_w || ctx->frame_hdr.kf.scale_h)
            vpx_internal_error(&ctx->error, VPX_CODEC_UNSUP_BITSTREAM,
                               "Spatial resampling not supported.");

        data += KEYFRAME_HEADER_SZ;
        sz -= KEYFRAME_HEADER_SZ;
    }


    /* Start the bitreader for the header/entropy partition */
    vp8dx_bool_init(&bool, data, ctx->frame_hdr.part0_sz);

    /* Skip the colorspace and clamping bits */
    if (ctx->frame_hdr.is_keyframe)
        if (bool_get_uint(&bool, 2))
            vpx_internal_error(&ctx->error, VPX_CODEC_UNSUP_BITSTREAM,
                               "Reserved bits not supported.");

    decode_segmentation_header(ctx, &bool, &ctx->segment_hdr);
    decode_loopfilter_header(ctx, &bool, &ctx->loopfilter_hdr);
    decode_and_init_token_partitions(ctx,
                                     &bool,
                                     data + ctx->frame_hdr.part0_sz,
                                     sz - ctx->frame_hdr.part0_sz,
                                     &ctx->token_hdr);
    decode_quantizer_header(ctx, &bool, &ctx->quant_hdr);
    decode_reference_header(ctx, &bool, &ctx->reference_hdr);
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
