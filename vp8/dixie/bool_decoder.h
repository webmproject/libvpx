/*
 *  Copyright (c) 2010 The VP8 project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license and patent
 *  grant that can be found in the LICENSE file in the root of the source
 *  tree. All contributing project authors may be found in the AUTHORS
 *  file in the root of the source tree.
 */


#ifndef BOOL_DECODER_H
#define BOOL_DECODER_H
#include <stddef.h>
#include <limits.h>
#include "vpx_ports/config.h"
#include "vpx_ports/mem.h"


typedef size_t vp8_bool_value_t;
# define VP8_BD_VALUE_SIZE ((int)sizeof(vp8_bool_value_t)*CHAR_BIT)


/*This is meant to be a large, positive constant that can still be efficiently
   loaded as an immediate (on platforms like ARM, for example).
  Even relatively modest values like 100 would work fine.*/
# define VP8_LOTS_OF_BITS (0x40000000)


DECLARE_ALIGNED(16, extern const unsigned int, vp8dx_bitreader_norm[256]);


struct bool_decoder
{
    const unsigned char *user_buffer_end;
    const unsigned char *user_buffer;
    vp8_bool_value_t     value;
    int                  count;
    unsigned int         range;
};


int  vp8dx_bool_init(struct bool_decoder *br, const unsigned char *source,
                     unsigned int source_sz);

void vp8dx_bool_fill(struct bool_decoder *br);


static int bool_get(struct bool_decoder *br, int probability)
{
    unsigned int bit = 0;
    vp8_bool_value_t value;
    unsigned int split;
    vp8_bool_value_t bigsplit;
    int count;
    unsigned int range;

    value = br->value;
    count = br->count;
    range = br->range;

    split = 1 + (((range - 1) * probability) >> 8);
    bigsplit = (vp8_bool_value_t)split << (VP8_BD_VALUE_SIZE - 8);

    range = split;

    if (value >= bigsplit)
    {
        range = br->range - split;
        value = value - bigsplit;
        bit = 1;
    }

    {
        register unsigned int shift = vp8dx_bitreader_norm[range];
        range <<= shift;
        value <<= shift;
        count -= shift;
    }

    br->value = value;
    br->count = count;
    br->range = range;

    if (count < 8)
        vp8dx_bool_fill(br);

    return bit;
}


static int bool_get_bit(struct bool_decoder *br)
{
    return bool_get(br, 128);
}


static int bool_get_uint(struct bool_decoder *br, int bits)
{
    int z = 0;
    int bit;

    for (bit = bits - 1; bit >= 0; bit--)
    {
        z |= (bool_get_bit(br) << bit);
    }

    return z;
}


static int bool_get_int(struct bool_decoder *br, int bits)
{
    int z = 0;
    int bit;

    for (bit = bits - 1; bit >= 0; bit--)
    {
        z |= (bool_get_bit(br) << bit);
    }

    return bool_get_bit(br) ? -z : z;
}


static int bool_maybe_get_int(struct bool_decoder *br, int bits)
{
    int z = 0;
    int bit;

    if (!bool_get_bit(br))
        return 0;

    for (bit = bits - 1; bit >= 0; bit--)
    {
        z |= (bool_get_bit(br) << bit);
    }

    return bool_get_bit(br) ? -z : z;
}
#endif
