/*
 *  Copyright (c) 2010 The VP8 project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */
#ifndef DIXIE_H
#define DIXIE_H

struct vp8_frame_hdr
{
    unsigned int is_keyframe;      /* Frame is a keyframe */
    unsigned int is_experimental;  /* Frame is a keyframe */
    unsigned int version;          /* Bitstream version */
    unsigned int is_shown;         /* Frame is to be displayed. */
    unsigned int part0_sz;         /* Partition 0 length, in bytes */

    struct vp8_kf_hdr
    {
        unsigned int w;        /* Width */
        unsigned int h;        /* Height */
        unsigned int scale_w;  /* Scaling factor, Width */
        unsigned int scale_h;  /* Scaling factor, Height */
    } kf;
};

vpx_codec_err_t
vp8_parse_frame_header(const unsigned char   *data,
                       unsigned int           sz,
                       struct vp8_frame_hdr  *hdr);

#endif
