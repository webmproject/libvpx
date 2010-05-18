/*
 *  Copyright (c) 2010 The VP8 project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license and patent
 *  grant that can be found in the LICENSE file in the root of the source
 *  tree. All contributing project authors may be found in the AUTHORS
 *  file in the root of the source tree.
 */


/****************************************************************************
 *
 *   Module Title :     gen_scalers.c
 *
 *   Description  :     Generic image scaling functions.
 *
 ***************************************************************************/

/****************************************************************************
*  Header Files
****************************************************************************/
#include "vpx_scale/vpxscale.h"

/****************************************************************************
*  Imports
****************************************************************************/

/****************************************************************************
 *
 *  ROUTINE       : horizontal_line_4_5_scale_c4
 *
 *  INPUTS        : const unsigned char *source : Pointer to source data.
 *                  unsigned int source_width    : Stride of source.
 *                  unsigned char *dest         : Pointer to destination data.
 *                  unsigned int dest_width      : Stride of destination (NOT USED).
 *
 *  OUTPUTS       : None.
 *
 *  RETURNS       : void
 *
 *  FUNCTION      : Copies horizontal line of pixels from source to
 *                  destination scaling up by 4 to 5.
 *
 *  SPECIAL NOTES : None.
 *
 ****************************************************************************/
static
void horizontal_line_4_5_scale_c64
(
    const unsigned char *source,
    unsigned int source_width,
    unsigned char *dest,
    unsigned int dest_width
)
{
    unsigned i;
    unsigned int ba, cb, dc, ed;
    unsigned char *restrict des = dest;
    unsigned int *restrict src = (unsigned int *)source;
    unsigned int const_51_205, const_102_154,
             const_205_51, const_154_102;

    unsigned int src_current, src_next;

    (void) dest_width;

    // Constants that are to be used for the filtering.  For
    //  best speed we are going to want to right shift by 16.
    //  In the generic version they were shift by 8, so put
    //  an extra 8 in now so that 16 will come out later.
    const_51_205 = 0x3300CD00; //_pack2 (51 << 8, 205 << 8);
    const_205_51 = 0xCD003300; //_pack2 (205 << 8, 51 << 8);
    const_102_154 = 0x66009A00; //_pack2 (102 << 8, 154 << 8);
    const_154_102 = 0x9A006600; //_pack2 (154 << 8, 102 << 8);

    // 5 points are needed to filter to give 5 output points.
    //  A load can pull up 4 at a time, and one needs to be
    //  "borrowed" from the next set of data.  So instead of
    //  loading those 5 points each time, "steal" a point from
    //  the next set and only load up 4 each time through.
    src_current = _mem4(src);

    for (i = 0; i < source_width - 4; i += 4)
    {
        src_next = _mem4(src++);

        // Reorder the data so that it is ready for the
        //  dot product.
        ba = _unpklu4(src_current);
        cb = _unpkhu4(_rotl(src_current, 8));
        dc = _unpkhu4(src_current);
        ed = _unpkhu4(_shrmb(src_next, src_current));

        // Use the dot product with round and shift.
        des [0] = src_current & 0xff;
        des [1] = _dotprsu2(ba, const_205_51);
        des [2] = _dotprsu2(cb, const_154_102);
        des [3] = _dotprsu2(dc, const_102_154);
        des [4] = _dotprsu2(ed, const_51_205);

        des += 5;

        // reuse loaded vales next time around.
        src_current = src_next;
    }

    // vp8_filter the last set of points.  Normally a point from the next set
    //  would be used, but there is no next set, so just fill.
    ba = _unpklu4(src_current);
    cb = _unpkhu4(_rotl(src_current, 8));
    dc = _unpkhu4(src_current);

    des [0] = src_current & 0xff;
    des [1] = _dotprsu2(ba, const_205_51);
    des [2] = _dotprsu2(cb, const_154_102);
    des [3] = _dotprsu2(dc, const_102_154);
    des [4] = src_current & 0xff;

}
/****************************************************************************
 *
 *  ROUTINE       : vertical_band_4_5_scale_c64
 *
 *  INPUTS        : unsigned char *dest    : Pointer to destination data.
 *                  unsigned int dest_pitch : Stride of destination data.
 *                  unsigned int dest_width : Width of destination data.
 *
 *  OUTPUTS       : None.
 *
 *  RETURNS       : void
 *
 *  FUNCTION      : Scales vertical band of pixels by scale 4 to 5. The
 *                  height of the band scaled is 4-pixels.
 *
 *  SPECIAL NOTES : The routine uses the first line of the band below
 *                  the current band.
 *
 ****************************************************************************/
static
void vertical_band_4_5_scale_c64(unsigned char *dest, unsigned int dest_pitch, unsigned int dest_width)
{
    unsigned int i;
    unsigned int a, b, c, d, e;
    unsigned int ba, cb, dc, ed;
    unsigned char *restrict src = dest;
    unsigned char *restrict des = dest;
    unsigned int const_51_205, const_102_154,
             const_205_51, const_154_102;

    const_51_205 = 0x3300CD00; //_pack2 (51 << 8, 205 << 8);
    const_205_51 = 0xCD003300; //_pack2 (205 << 8, 51 << 8);
    const_102_154 = 0x66009A00; //_pack2 (102 << 8, 154 << 8);
    const_154_102 = 0x9A006600; //_pack2 (154 << 8, 102 << 8);

    // Force a loop unroll here so that there is not such a
    //  dependancy.
    a = src [0];
    b = src [dest_pitch];
    c = src [dest_pitch*2];
    d = src [dest_pitch*3];
    e = src [dest_pitch*5];
    src ++;

    for (i = 0; i < dest_width; i++)
    {
        ba = _pack2(b, a);
        cb = _pack2(c, b);
        dc = _pack2(d, c);
        ed = _pack2(e, d);

        a = src [0];
        b = src [dest_pitch];
        c = src [dest_pitch*2];
        d = src [dest_pitch*3];
        e = src [dest_pitch*5];
        src ++;

        des [dest_pitch] = _dotprsu2(ba, const_205_51);
        des [dest_pitch*2] = _dotprsu2(cb, const_154_102);
        des [dest_pitch*3] = _dotprsu2(dc, const_102_154);
        des [dest_pitch*4] = _dotprsu2(ed, const_51_205);

        des ++;
    }
}

/****************************************************************************
 *
 *  ROUTINE       : last_vertical_band_4_5_scale_c64
 *
 *  INPUTS        : unsigned char *dest    : Pointer to destination data.
 *                  unsigned int dest_pitch : Stride of destination data.
 *                  unsigned int dest_width : Width of destination data.
 *
 *  OUTPUTS       : None.
 *
 *  RETURNS       : void
 *
 *  FUNCTION      : Scales last vertical band of pixels by scale 4 to 5. The
 *                  height of the band scaled is 4-pixels.
 *
 *  SPECIAL NOTES : The routine does not have available the first line of
 *                  the band below the current band, since this is the
 *                  last band.
 *
 ****************************************************************************/
static
void last_vertical_band_4_5_scale_c64(unsigned char *dest, unsigned int dest_pitch, unsigned int dest_width)
{
    unsigned int i;
    unsigned int a, b, c, d;
    unsigned int ba, cb, dc;
    unsigned char *restrict src = dest;
    unsigned char *restrict des = dest;
    unsigned int const_102_154, const_205_51, const_154_102;

    const_205_51 = 0xCD003300; //_pack2 (205 << 8, 51 << 8);
    const_102_154 = 0x66009A00; //_pack2 (102 << 8, 154 << 8);
    const_154_102 = 0x9A006600; //_pack2 (154 << 8, 102 << 8);

    a = src [0];
    b = src [dest_pitch];
    c = src [dest_pitch*2];
    d = src [dest_pitch*3];
    src ++;

    for (i = 0; i < dest_width; ++i)
    {
        ba = _pack2(b, a);
        cb = _pack2(c, b);
        dc = _pack2(d, c);

        a = src [0];
        b = src [dest_pitch];
        c = src [dest_pitch*2];
        d = src [dest_pitch*3];
        src ++;

        des [dest_pitch] = _dotprsu2(ba, const_205_51);
        des [dest_pitch*2] = _dotprsu2(cb, const_154_102);
        des [dest_pitch*3] = _dotprsu2(dc, const_102_154);
        des [dest_pitch*4] = (unsigned char) d;

        des++;
    }
}

/****************************************************************************
 *
 *  ROUTINE       : horizontal_line_3_5_scale_c64
 *
 *  INPUTS        : const unsigned char *source : Pointer to source data.
 *                  unsigned int source_width    : Stride of source.
 *                  unsigned char *dest         : Pointer to destination data.
 *                  unsigned int dest_width      : Stride of destination (NOT USED).
 *
 *  OUTPUTS       : None.
 *
 *  RETURNS       : void
 *
 *  FUNCTION      : Copies horizontal line of pixels from source to
 *                  destination scaling up by 3 to 5.
 *
 *  SPECIAL NOTES : None.
 *
 *
 ****************************************************************************/
static
void horizontal_line_3_5_scale_c64
(
    const unsigned char *source,
    unsigned int source_width,
    unsigned char *dest,
    unsigned int dest_width
)
{
    unsigned int i;
    unsigned int ba, cb, dc;
    unsigned int src_current;
    unsigned char *restrict des = dest;
    unsigned char *restrict src = (unsigned char *)source;
    unsigned int const_51_205, const_102_154,
             const_205_51, const_154_102;

    (void) dest_width;

    const_51_205 = 0x3300CD00; //_pack2 (51 << 8, 205 << 8);
    const_205_51 = 0xCD003300; //_pack2 (205 << 8, 51 << 8);
    const_102_154 = 0x66009A00; //_pack2 (102 << 8, 154 << 8);
    const_154_102 = 0x9A006600; //_pack2 (154 << 8, 102 << 8);

    for (i = 0; i < source_width - 3; i += 3)
    {
        src_current = _mem4(src);

        // Reorder the data so that it is ready for the
        //  dot product.
        ba = _unpklu4(src_current);
        cb = _unpkhu4(_rotl(src_current, 8));
        dc = _unpkhu4(src_current);

        des [0] = src_current & 0xff;
        des [1] = _dotprsu2(ba, const_154_102);
        des [2] = _dotprsu2(cb, const_51_205);
        des [3] = _dotprsu2(cb, const_205_51);
        des [4] = _dotprsu2(dc, const_102_154);

        src += 3;
        des += 5;
    }

    src_current = _mem4(src);

    ba = _unpklu4(src_current);
    cb = _unpkhu4(_rotl(src_current, 8));
    dc = _unpkhu4(src_current);


    des [0] = src_current & 0xff;
    des [1] = _dotprsu2(ba, const_154_102);
    des [2] = _dotprsu2(cb, const_51_205);
    des [3] = _dotprsu2(cb, const_205_51);
    des [4] = dc & 0xff;

}

/****************************************************************************
 *
 *  ROUTINE       : vertical_band_3_5_scale_c64
 *
 *  INPUTS        : unsigned char *dest    : Pointer to destination data.
 *                  unsigned int dest_pitch : Stride of destination data.
 *                  unsigned int dest_width : Width of destination data.
 *
 *  OUTPUTS       : None.
 *
 *  RETURNS       : void
 *
 *  FUNCTION      : Scales vertical band of pixels by scale 3 to 5. The
 *                  height of the band scaled is 3-pixels.
 *
 *  SPECIAL NOTES : The routine uses the first line of the band below
 *                  the current band.
 *
 ****************************************************************************/
static
void vertical_band_3_5_scale_c64(unsigned char *dest, unsigned int dest_pitch, unsigned int dest_width)
{
    unsigned int i;
    unsigned int a, b, c, d;
    unsigned int ba, cb, dc;
    unsigned char *restrict src = dest;
    unsigned char *restrict des = dest;
    unsigned int const_51_205, const_102_154,
             const_205_51, const_154_102;

    const_51_205 = 0x3300CD00; //_pack2 (51 << 8, 205 << 8);
    const_205_51 = 0xCD003300; //_pack2 (205 << 8, 51 << 8);
    const_102_154 = 0x66009A00; //_pack2 (102 << 8, 154 << 8);
    const_154_102 = 0x9A006600; //_pack2 (154 << 8, 102 << 8);

    a = src [0];
    b = src [dest_pitch];
    c = src [dest_pitch*2];
    d = src [dest_pitch*5];
    src ++;

    for (i = 0; i < dest_width; i++)
    {
        ba = _pack2(b, a);
        cb = _pack2(c, b);
        dc = _pack2(d, c);

        a = src [0];
        b = src [dest_pitch];
        c = src [dest_pitch*2];
        d = src [dest_pitch*5];
        src ++;

        des [dest_pitch]   = _dotprsu2(ba, const_154_102);
        des [dest_pitch*2] = _dotprsu2(cb, const_51_205);
        des [dest_pitch*3] = _dotprsu2(cb, const_205_51);
        des [dest_pitch*4] = _dotprsu2(dc, const_102_154);

        des++;
    }
}

/****************************************************************************
 *
 *  ROUTINE       : last_vertical_band_3_5_scale_c64
 *
 *  INPUTS        : unsigned char *dest    : Pointer to destination data.
 *                  unsigned int dest_pitch : Stride of destination data.
 *                  unsigned int dest_width : Width of destination data.
 *
 *  OUTPUTS       : None.
 *
 *  RETURNS       : void
 *
 *  FUNCTION      : Scales last vertical band of pixels by scale 3 to 5. The
 *                  height of the band scaled is 3-pixels.
 *
 *  SPECIAL NOTES : The routine does not have available the first line of
 *                  the band below the current band, since this is the
 *                  last band.
 *
 ****************************************************************************/
static
void last_vertical_band_3_5_scale_c64(unsigned char *dest, unsigned int dest_pitch, unsigned int dest_width)
{
    unsigned int i;
    unsigned int a, b, c;
    unsigned int ba, cb;
    unsigned char *restrict src = dest;
    unsigned char *restrict des = dest;
    unsigned int const_51_205, const_205_51, const_154_102;

    const_51_205 = 0x3300CD00; //_pack2 (51 << 8, 205 << 8);
    const_205_51 = 0xCD003300; //_pack2 (205 << 8, 51 << 8);
    const_154_102 = 0x9A006600; //_pack2 (154 << 8, 102 << 8);

    a = src [0];
    b = src [dest_pitch];
    c = src [dest_pitch*2];
    src ++;

    for (i = 0; i < dest_width; ++i)
    {
        ba = _pack2(b, a);
        cb = _pack2(c, b);

        a = src [0];
        b = src [dest_pitch];
        c = src [dest_pitch*2];
        src ++;

        des [dest_pitch]   = _dotprsu2(ba, const_154_102);
        des [dest_pitch*2] = _dotprsu2(cb, const_51_205);
        des [dest_pitch*3] = _dotprsu2(cb, const_205_51);
        des [dest_pitch*4] = (unsigned char)(c) ;

        des++;
    }
}

/****************************************************************************
 *
 *  ROUTINE       : horizontal_line_1_2_scale_c64
 *
 *  INPUTS        : const unsigned char *source : Pointer to source data.
 *                  unsigned int source_width    : Stride of source.
 *                  unsigned char *dest         : Pointer to destination data.
 *                  unsigned int dest_width      : Stride of destination (NOT USED).
 *
 *  OUTPUTS       : None.
 *
 *  RETURNS       : void
 *
 *  FUNCTION      : Copies horizontal line of pixels from source to
 *                  destination scaling up by 1 to 2.
 *
 *  SPECIAL NOTES : source width must be a multiple of 4.
 *
 ****************************************************************************/
void horizontal_line_1_2_scale_c64
(
    const unsigned char *source,
    unsigned int source_width,
    unsigned char *dest,
    unsigned int dest_width
)
{
    unsigned int i;
    unsigned char *restrict des = dest;
    unsigned char *restrict src = (unsigned char *)source;
    unsigned int src7_4i, src4_1i, src3_0i;
    unsigned int a4_0i, ahi, alo;
    double src7_0d, src3_0d;
    const unsigned int k01 = 0x01010101;

    for (i = 0; i < source_width / 4; i += 1)
    {
        // Load up the data from src.  Here a wide load is
        //  used to get 8 bytes at once, only 5 will be used
        //  for the actual computation.
        src7_0d = _memd8(src);
        src3_0i = _lo(src7_0d);
        src7_4i = _hi(src7_0d);

        // Need to average between points.  Shift byte 5 into
        //  the lower word.  This will result in bytes 5-1
        //  averaged with 4-0.
        src4_1i = _shrmb(src7_4i, src3_0i);
        a4_0i = _avgu4(src4_1i, src3_0i);

        // Expand the data out. Could do an unpack, however
        //  all but the multiply units are getting pretty hard
        //  here the multiply unit can take some of the computations.
        src3_0d = _mpyu4(src3_0i, k01);

        // The averages need to be unpacked so that they are in 16
        //  bit form and will be able to be interleaved with the
        //  original data
        ahi = _unpkhu4(a4_0i);
        alo = _unpklu4(a4_0i);

        ahi = _swap4(ahi);
        alo = _swap4(alo);

        // Mix the average result in with the orginal data.
        ahi = _hi(src3_0d) | ahi;
        alo = _lo(src3_0d) | alo;

        _memd8(des) = _itod(ahi, alo);

        des += 8;
        src += 4;
    }
}


/****************************************************************************
 *
 *  ROUTINE       : vertical_band_1_2_scale_c64
 *
 *  INPUTS        : unsigned char *dest    : Pointer to destination data.
 *                  unsigned int dest_pitch : Stride of destination data.
 *                  unsigned int dest_width : Width of destination data.
 *
 *  OUTPUTS       : None.
 *
 *  RETURNS       : void
 *
 *  FUNCTION      : Scales vertical band of pixels by scale 1 to 2. The
 *                  height of the band scaled is 1-pixel.
 *
 *  SPECIAL NOTES : The routine uses the first line of the band below
 *                  the current band.
 *                  Destination width must be a multiple of 4.  Because the
 *                  intput must be, therefore the output must be.
 *
 ****************************************************************************/
static
void vertical_band_1_2_scale_c64(unsigned char *dest, unsigned int dest_pitch, unsigned int dest_width)
{
    unsigned int i;
    unsigned int a, b;
    unsigned int *restrict line_a = (unsigned int *)dest;
    unsigned int *restrict line_b = (unsigned int *)(dest + (dest_pitch * 2));
    unsigned int *restrict des = (unsigned int *)(dest + dest_pitch);

    for (i = 0; i < dest_width / 4; i++)
    {
        a = _mem4(line_a++);
        b = _mem4(line_b++);

        _mem4(des++) = _avgu4(a, b);
    }
}

/****************************************************************************
 *
 *  ROUTINE       : last_vertical_band_1_2_scale_c64
 *
 *  INPUTS        : unsigned char *dest    : Pointer to destination data.
 *                  unsigned int dest_pitch : Stride of destination data.
 *                  unsigned int dest_width : Width of destination data.
 *
 *  OUTPUTS       : None.
 *
 *  RETURNS       : void
 *
 *  FUNCTION      : Scales last vertical band of pixels by scale 1 to 2. The
 *                  height of the band scaled is 1-pixel.
 *
 *  SPECIAL NOTES : The routine does not have available the first line of
 *                  the band below the current band, since this is the
 *                  last band.  Again, width must be a multiple of 4.
 *
 ****************************************************************************/
static
void last_vertical_band_1_2_scale_c64(unsigned char *dest, unsigned int dest_pitch, unsigned int dest_width)
{
    unsigned int i;
    unsigned int *restrict src = (unsigned int *)dest;
    unsigned int *restrict des = (unsigned int *)(dest + dest_pitch);

    for (i = 0; i < dest_width / 4; ++i)
    {
        _mem4(des++) = _mem4(src++);
    }
}

void
register_generic_scalers(void)
{
    vp8_horizontal_line_1_2_scale        = horizontal_line_1_2_scale_c64;
    vp8_vertical_band_1_2_scale          = vertical_band_1_2_scale_c64;
    vp8_last_vertical_band_1_2_scale      = last_vertical_band_1_2_scale_c64;
    vp8_horizontal_line_3_5_scale        = horizontal_line_3_5_scale_c64;
    vp8_vertical_band_3_5_scale          = vertical_band_3_5_scale_c64;
    vp8_last_vertical_band_3_5_scale      = last_vertical_band_3_5_scale_c64;
    vp8_horizontal_line_4_5_scale        = horizontal_line_4_5_scale_c64;
    vp8_vertical_band_4_5_scale          = vertical_band_4_5_scale_c64;
    vp8_last_vertical_band_4_5_scale      = last_vertical_band_4_5_scale_c64;
}
