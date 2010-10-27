/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */


#include <stdlib.h>

unsigned int vp8_sad16x16_c(
    const unsigned char *src_ptr,
    int  src_stride,
    const unsigned char *ref_ptr,
    int  ref_stride,
    int max_sad)
{

    int r, c;
    unsigned int sad = 0;

    for (r = 0; r < 16; r++)
    {
        for (c = 0; c < 16; c++)
        {
            sad += abs(src_ptr[c] - ref_ptr[c]);
        }

        src_ptr += src_stride;
        ref_ptr += ref_stride;
    }

    return sad;
}


static __inline
unsigned int sad_mx_n_c(
    const unsigned char *src_ptr,
    int  src_stride,
    const unsigned char *ref_ptr,
    int  ref_stride,
    int m,
    int n)
{

    int r, c;
    unsigned int sad = 0;

    for (r = 0; r < n; r++)
    {
        for (c = 0; c < m; c++)
        {
            sad += abs(src_ptr[c] - ref_ptr[c]);
        }

        src_ptr += src_stride;
        ref_ptr += ref_stride;
    }

    return sad;
}


unsigned int vp8_sad8x8_c(
    const unsigned char *src_ptr,
    int  src_stride,
    const unsigned char *ref_ptr,
    int  ref_stride,
    int max_sad)
{

    return sad_mx_n_c(src_ptr, src_stride, ref_ptr, ref_stride, 8, 8);
}


unsigned int vp8_sad16x8_c(
    const unsigned char *src_ptr,
    int  src_stride,
    const unsigned char *ref_ptr,
    int  ref_stride,
    int max_sad)
{

    return sad_mx_n_c(src_ptr, src_stride, ref_ptr, ref_stride, 16, 8);

}


unsigned int vp8_sad8x16_c(
    const unsigned char *src_ptr,
    int  src_stride,
    const unsigned char *ref_ptr,
    int  ref_stride,
    int max_sad)
{

    return sad_mx_n_c(src_ptr, src_stride, ref_ptr, ref_stride, 8, 16);
}


unsigned int vp8_sad4x4_c(
    const unsigned char *src_ptr,
    int  src_stride,
    const unsigned char *ref_ptr,
    int  ref_stride,
    int max_sad)
{

    return sad_mx_n_c(src_ptr, src_stride, ref_ptr, ref_stride, 4, 4);
}

void vp8_sad16x16x3_c(
    const unsigned char *src_ptr,
    int  src_stride,
    const unsigned char *ref_ptr,
    int  ref_stride,
    unsigned int *sad_array
)
{
    sad_array[0] = vp8_sad16x16_c(src_ptr, src_stride, ref_ptr  , ref_stride, 0x7fffffff);
    sad_array[1] = vp8_sad16x16_c(src_ptr, src_stride, ref_ptr + 1, ref_stride, 0x7fffffff);
    sad_array[2] = vp8_sad16x16_c(src_ptr, src_stride, ref_ptr + 2, ref_stride, 0x7fffffff);
}

void vp8_sad16x16x8_c(
    const unsigned char *src_ptr,
    int  src_stride,
    const unsigned char *ref_ptr,
    int  ref_stride,
    unsigned short *sad_array
)
{
    sad_array[0] = (unsigned short)vp8_sad16x16_c(src_ptr, src_stride, ref_ptr  , ref_stride, 0x7fffffff);
    sad_array[1] = (unsigned short)vp8_sad16x16_c(src_ptr, src_stride, ref_ptr + 1, ref_stride, 0x7fffffff);
    sad_array[2] = (unsigned short)vp8_sad16x16_c(src_ptr, src_stride, ref_ptr + 2, ref_stride, 0x7fffffff);
    sad_array[3] = (unsigned short)vp8_sad16x16_c(src_ptr, src_stride, ref_ptr + 3 , ref_stride, 0x7fffffff);
    sad_array[4] = (unsigned short)vp8_sad16x16_c(src_ptr, src_stride, ref_ptr + 4, ref_stride, 0x7fffffff);
    sad_array[5] = (unsigned short)vp8_sad16x16_c(src_ptr, src_stride, ref_ptr + 5, ref_stride, 0x7fffffff);
    sad_array[6] = (unsigned short)vp8_sad16x16_c(src_ptr, src_stride, ref_ptr + 6 , ref_stride, 0x7fffffff);
    sad_array[7] = (unsigned short)vp8_sad16x16_c(src_ptr, src_stride, ref_ptr + 7, ref_stride, 0x7fffffff);
}

void vp8_sad16x8x3_c(
    const unsigned char *src_ptr,
    int  src_stride,
    const unsigned char *ref_ptr,
    int  ref_stride,
    unsigned int *sad_array
)
{
    sad_array[0] = vp8_sad16x8_c(src_ptr, src_stride, ref_ptr  , ref_stride, 0x7fffffff);
    sad_array[1] = vp8_sad16x8_c(src_ptr, src_stride, ref_ptr + 1, ref_stride, 0x7fffffff);
    sad_array[2] = vp8_sad16x8_c(src_ptr, src_stride, ref_ptr + 2, ref_stride, 0x7fffffff);
}

void vp8_sad16x8x8_c(
    const unsigned char *src_ptr,
    int  src_stride,
    const unsigned char *ref_ptr,
    int  ref_stride,
    unsigned short *sad_array
)
{
    sad_array[0] = (unsigned short)vp8_sad16x8_c(src_ptr, src_stride, ref_ptr  , ref_stride, 0x7fffffff);
    sad_array[1] = (unsigned short)vp8_sad16x8_c(src_ptr, src_stride, ref_ptr + 1, ref_stride, 0x7fffffff);
    sad_array[2] = (unsigned short)vp8_sad16x8_c(src_ptr, src_stride, ref_ptr + 2, ref_stride, 0x7fffffff);
    sad_array[3] = (unsigned short)vp8_sad16x8_c(src_ptr, src_stride, ref_ptr + 3 , ref_stride, 0x7fffffff);
    sad_array[4] = (unsigned short)vp8_sad16x8_c(src_ptr, src_stride, ref_ptr + 4, ref_stride, 0x7fffffff);
    sad_array[5] = (unsigned short)vp8_sad16x8_c(src_ptr, src_stride, ref_ptr + 5, ref_stride, 0x7fffffff);
    sad_array[6] = (unsigned short)vp8_sad16x8_c(src_ptr, src_stride, ref_ptr + 6 , ref_stride, 0x7fffffff);
    sad_array[7] = (unsigned short)vp8_sad16x8_c(src_ptr, src_stride, ref_ptr + 7, ref_stride, 0x7fffffff);
}

void vp8_sad8x8x3_c(
    const unsigned char *src_ptr,
    int  src_stride,
    const unsigned char *ref_ptr,
    int  ref_stride,
    unsigned int *sad_array
)
{
    sad_array[0] = vp8_sad8x8_c(src_ptr, src_stride, ref_ptr  , ref_stride, 0x7fffffff);
    sad_array[1] = vp8_sad8x8_c(src_ptr, src_stride, ref_ptr + 1, ref_stride, 0x7fffffff);
    sad_array[2] = vp8_sad8x8_c(src_ptr, src_stride, ref_ptr + 2, ref_stride, 0x7fffffff);
}

void vp8_sad8x8x8_c(
    const unsigned char *src_ptr,
    int  src_stride,
    const unsigned char *ref_ptr,
    int  ref_stride,
    unsigned short *sad_array
)
{
    sad_array[0] = (unsigned short)vp8_sad8x8_c(src_ptr, src_stride, ref_ptr  , ref_stride, 0x7fffffff);
    sad_array[1] = (unsigned short)vp8_sad8x8_c(src_ptr, src_stride, ref_ptr + 1, ref_stride, 0x7fffffff);
    sad_array[2] = (unsigned short)vp8_sad8x8_c(src_ptr, src_stride, ref_ptr + 2, ref_stride, 0x7fffffff);
    sad_array[3] = (unsigned short)vp8_sad8x8_c(src_ptr, src_stride, ref_ptr + 3 , ref_stride, 0x7fffffff);
    sad_array[4] = (unsigned short)vp8_sad8x8_c(src_ptr, src_stride, ref_ptr + 4, ref_stride, 0x7fffffff);
    sad_array[5] = (unsigned short)vp8_sad8x8_c(src_ptr, src_stride, ref_ptr + 5, ref_stride, 0x7fffffff);
    sad_array[6] = (unsigned short)vp8_sad8x8_c(src_ptr, src_stride, ref_ptr + 6 , ref_stride, 0x7fffffff);
    sad_array[7] = (unsigned short)vp8_sad8x8_c(src_ptr, src_stride, ref_ptr + 7, ref_stride, 0x7fffffff);
}

void vp8_sad8x16x3_c(
    const unsigned char *src_ptr,
    int  src_stride,
    const unsigned char *ref_ptr,
    int  ref_stride,
    unsigned int *sad_array
)
{
    sad_array[0] = vp8_sad8x16_c(src_ptr, src_stride, ref_ptr  , ref_stride, 0x7fffffff);
    sad_array[1] = vp8_sad8x16_c(src_ptr, src_stride, ref_ptr + 1, ref_stride, 0x7fffffff);
    sad_array[2] = vp8_sad8x16_c(src_ptr, src_stride, ref_ptr + 2, ref_stride, 0x7fffffff);
}

void vp8_sad8x16x8_c(
    const unsigned char *src_ptr,
    int  src_stride,
    const unsigned char *ref_ptr,
    int  ref_stride,
    unsigned short *sad_array
)
{
    sad_array[0] = (unsigned short)vp8_sad8x16_c(src_ptr, src_stride, ref_ptr  , ref_stride, 0x7fffffff);
    sad_array[1] = (unsigned short)vp8_sad8x16_c(src_ptr, src_stride, ref_ptr + 1, ref_stride, 0x7fffffff);
    sad_array[2] = (unsigned short)vp8_sad8x16_c(src_ptr, src_stride, ref_ptr + 2, ref_stride, 0x7fffffff);
    sad_array[3] = (unsigned short)vp8_sad8x16_c(src_ptr, src_stride, ref_ptr + 3 , ref_stride, 0x7fffffff);
    sad_array[4] = (unsigned short)vp8_sad8x16_c(src_ptr, src_stride, ref_ptr + 4, ref_stride, 0x7fffffff);
    sad_array[5] = (unsigned short)vp8_sad8x16_c(src_ptr, src_stride, ref_ptr + 5, ref_stride, 0x7fffffff);
    sad_array[6] = (unsigned short)vp8_sad8x16_c(src_ptr, src_stride, ref_ptr + 6 , ref_stride, 0x7fffffff);
    sad_array[7] = (unsigned short)vp8_sad8x16_c(src_ptr, src_stride, ref_ptr + 7, ref_stride, 0x7fffffff);
}

void vp8_sad4x4x3_c(
    const unsigned char *src_ptr,
    int  src_stride,
    const unsigned char *ref_ptr,
    int  ref_stride,
    unsigned int *sad_array
)
{
    sad_array[0] = vp8_sad4x4_c(src_ptr, src_stride, ref_ptr  , ref_stride, 0x7fffffff);
    sad_array[1] = vp8_sad4x4_c(src_ptr, src_stride, ref_ptr + 1, ref_stride, 0x7fffffff);
    sad_array[2] = vp8_sad4x4_c(src_ptr, src_stride, ref_ptr + 2, ref_stride, 0x7fffffff);
}

void vp8_sad4x4x8_c(
    const unsigned char *src_ptr,
    int  src_stride,
    const unsigned char *ref_ptr,
    int  ref_stride,
    unsigned short *sad_array
)
{
    sad_array[0] = (unsigned short)vp8_sad4x4_c(src_ptr, src_stride, ref_ptr  , ref_stride, 0x7fffffff);
    sad_array[1] = (unsigned short)vp8_sad4x4_c(src_ptr, src_stride, ref_ptr + 1, ref_stride, 0x7fffffff);
    sad_array[2] = (unsigned short)vp8_sad4x4_c(src_ptr, src_stride, ref_ptr + 2, ref_stride, 0x7fffffff);
    sad_array[3] = (unsigned short)vp8_sad4x4_c(src_ptr, src_stride, ref_ptr + 3 , ref_stride, 0x7fffffff);
    sad_array[4] = (unsigned short)vp8_sad4x4_c(src_ptr, src_stride, ref_ptr + 4, ref_stride, 0x7fffffff);
    sad_array[5] = (unsigned short)vp8_sad4x4_c(src_ptr, src_stride, ref_ptr + 5, ref_stride, 0x7fffffff);
    sad_array[6] = (unsigned short)vp8_sad4x4_c(src_ptr, src_stride, ref_ptr + 6 , ref_stride, 0x7fffffff);
    sad_array[7] = (unsigned short)vp8_sad4x4_c(src_ptr, src_stride, ref_ptr + 7, ref_stride, 0x7fffffff);
}

void vp8_sad16x16x4d_c(
    const unsigned char *src_ptr,
    int  src_stride,
    unsigned char *ref_ptr[],
    int  ref_stride,
    unsigned int *sad_array
)
{
    sad_array[0] = vp8_sad16x16_c(src_ptr, src_stride, ref_ptr[0], ref_stride, 0x7fffffff);
    sad_array[1] = vp8_sad16x16_c(src_ptr, src_stride, ref_ptr[1], ref_stride, 0x7fffffff);
    sad_array[2] = vp8_sad16x16_c(src_ptr, src_stride, ref_ptr[2], ref_stride, 0x7fffffff);
    sad_array[3] = vp8_sad16x16_c(src_ptr, src_stride, ref_ptr[3], ref_stride, 0x7fffffff);
}

void vp8_sad16x8x4d_c(
    const unsigned char *src_ptr,
    int  src_stride,
    unsigned char *ref_ptr[],
    int  ref_stride,
    unsigned int *sad_array
)
{
    sad_array[0] = vp8_sad16x8_c(src_ptr, src_stride, ref_ptr[0], ref_stride, 0x7fffffff);
    sad_array[1] = vp8_sad16x8_c(src_ptr, src_stride, ref_ptr[1], ref_stride, 0x7fffffff);
    sad_array[2] = vp8_sad16x8_c(src_ptr, src_stride, ref_ptr[2], ref_stride, 0x7fffffff);
    sad_array[3] = vp8_sad16x8_c(src_ptr, src_stride, ref_ptr[3], ref_stride, 0x7fffffff);
}

void vp8_sad8x8x4d_c(
    const unsigned char *src_ptr,
    int  src_stride,
    unsigned char *ref_ptr[],
    int  ref_stride,
    unsigned int *sad_array
)
{
    sad_array[0] = vp8_sad8x8_c(src_ptr, src_stride, ref_ptr[0], ref_stride, 0x7fffffff);
    sad_array[1] = vp8_sad8x8_c(src_ptr, src_stride, ref_ptr[1], ref_stride, 0x7fffffff);
    sad_array[2] = vp8_sad8x8_c(src_ptr, src_stride, ref_ptr[2], ref_stride, 0x7fffffff);
    sad_array[3] = vp8_sad8x8_c(src_ptr, src_stride, ref_ptr[3], ref_stride, 0x7fffffff);
}

void vp8_sad8x16x4d_c(
    const unsigned char *src_ptr,
    int  src_stride,
    unsigned char *ref_ptr[],
    int  ref_stride,
    unsigned int *sad_array
)
{
    sad_array[0] = vp8_sad8x16_c(src_ptr, src_stride, ref_ptr[0], ref_stride, 0x7fffffff);
    sad_array[1] = vp8_sad8x16_c(src_ptr, src_stride, ref_ptr[1], ref_stride, 0x7fffffff);
    sad_array[2] = vp8_sad8x16_c(src_ptr, src_stride, ref_ptr[2], ref_stride, 0x7fffffff);
    sad_array[3] = vp8_sad8x16_c(src_ptr, src_stride, ref_ptr[3], ref_stride, 0x7fffffff);
}

void vp8_sad4x4x4d_c(
    const unsigned char *src_ptr,
    int  src_stride,
    unsigned char *ref_ptr[],
    int  ref_stride,
    unsigned int *sad_array
)
{
    sad_array[0] = vp8_sad4x4_c(src_ptr, src_stride, ref_ptr[0], ref_stride, 0x7fffffff);
    sad_array[1] = vp8_sad4x4_c(src_ptr, src_stride, ref_ptr[1], ref_stride, 0x7fffffff);
    sad_array[2] = vp8_sad4x4_c(src_ptr, src_stride, ref_ptr[2], ref_stride, 0x7fffffff);
    sad_array[3] = vp8_sad4x4_c(src_ptr, src_stride, ref_ptr[3], ref_stride, 0x7fffffff);
}
