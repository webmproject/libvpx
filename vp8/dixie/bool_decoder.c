/*
 *  Copyright (c) 2010 The VP8 project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license and patent
 *  grant that can be found in the LICENSE file in the root of the source
 *  tree. All contributing project authors may be found in the AUTHORS
 *  file in the root of the source tree.
 */
#include "bool_decoder.h"
#include "vpx_ports/mem.h"


DECLARE_ALIGNED(16, const unsigned int, vp8dx_bool_norm[256]) =
{
    0, 7, 6, 6, 5, 5, 5, 5, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
    2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
};


int vp8dx_bool_init(struct bool_decoder *br, const unsigned char *source,
                    unsigned int source_sz)
{
    br->user_buffer_end = source + source_sz;
    br->user_buffer     = source;
    br->value    = 0;
    br->count    = 0;
    br->range    = 255;

    if (source_sz && !source)
        return 1;

    /* Populate the buffer */
    vp8dx_bool_fill(br);

    return 0;
}


void vp8dx_bool_fill(struct bool_decoder *br)
{
    const unsigned char *ptr;
    const unsigned char *end;
    vp8_bool_value_t     value;
    int                  count;
    end = br->user_buffer_end;
    ptr = br->user_buffer;
    value = br->value;
    count = br->count;

    for (;;)
    {
        if (ptr >= end)
        {
            count = VP8_LOTS_OF_BITS;
            break;
        }

        if (count > VP8_BD_VALUE_SIZE - 8)
            break;

        count += 8;
        value |= (vp8_bool_value_t) * ptr++ << (VP8_BD_VALUE_SIZE - count);
    }

    br->user_buffer = ptr;
    br->value = value;
    br->count = count;
}
