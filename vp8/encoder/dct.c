/*
 *  Copyright (c) 2010 The VP8 project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */


#include <math.h>


static const short dct_matrix2[4][4] =
{
    { 23170,  30274,  23170, 12540 },
    { 23170,  12540, -23170, -30274 },
    { 23170, -12540, -23170, 30274 },
    { 23170, -30274,  23170, -12540 }
};

static const short dct_matrix1[4][4] =
{
    { 23170,  23170,  23170,  23170 },
    { 30274,  12540, -12540, -30274 },
    { 23170, -23170, -23170,  23170 },
    { 12540, -30274,  30274, -12540 }
};


#define _1STSTAGESHIFT           14
#define _1STSTAGEROUNDING        (1<<( _1STSTAGESHIFT-1))
#define _2NDSTAGESHIFT           16
#define _2NDSTAGEROUNDING        (1<<( _2NDSTAGESHIFT-1))

// using matrix multiply
void vp8_short_fdct4x4_c(short *input, short *output, int pitch)
{
    int i, j, k;
    short temp[4][4];
    int sumtemp;
    pitch >>= 1;

    for (i = 0; i < 4; i++)
    {
        for (j = 0; j < 4; j++)
        {
            sumtemp = 0;

            for (k = 0; k < 4; k++)
            {
                sumtemp += input[i*pitch+k] * dct_matrix2[k][j];

            }

            temp[i][j] = (short)((sumtemp + _1STSTAGEROUNDING) >> _1STSTAGESHIFT);
        }
    }


    for (i = 0; i < 4; i++)
    {
        for (j = 0; j < 4; j++)
        {
            sumtemp = 0;

            for (k = 0; k < 4; k++)
            {
                sumtemp += dct_matrix1[i][ k] * temp[k][ j];
            }

            output[i*4+j] = (short)((sumtemp + _2NDSTAGEROUNDING) >> _2NDSTAGESHIFT);
        }
    }

}


void vp8_short_fdct8x4_c(short *input, short *output, int pitch)
{
    vp8_short_fdct4x4_c(input,   output,    pitch);
    vp8_short_fdct4x4_c(input + 4, output + 16, pitch);
}


static const signed short x_c1 = 60547;
static const signed short x_c2 = 46341;
static const signed short x_c3 = 25080;

void vp8_fast_fdct4x4_c(short *input, short *output, int pitch)
{
    int i;
    int a1, b1, c1, d1;
    int a2, b2, c2, d2;
    short *ip = input;

    short *op = output;
    int temp1, temp2;

    for (i = 0; i < 4; i++)
    {
        a1 = (ip[0] + ip[3]) * 2;
        b1 = (ip[1] + ip[2]) * 2;
        c1 = (ip[1] - ip[2]) * 2;
        d1 = (ip[0] - ip[3]) * 2;

        temp1 = a1 + b1;
        temp2 = a1 - b1;

        op[0] = ((temp1 * x_c2) >> 16) + temp1;
        op[2] = ((temp2 * x_c2) >> 16) + temp2;

        temp1 = (c1 * x_c3) >> 16;
        temp2 = ((d1 * x_c1) >> 16) + d1;

        op[1] = temp1 + temp2;

        temp1 = (d1 * x_c3) >> 16;
        temp2 = ((c1 * x_c1) >> 16) + c1;

        op[3] = temp1 - temp2;

        ip += pitch / 2;
        op += 4;
    }

    ip = output;
    op = output;

    for (i = 0; i < 4; i++)
    {

        a1 = ip[0] + ip[12];
        b1 = ip[4] + ip[8];
        c1 = ip[4] - ip[8];
        d1 = ip[0] - ip[12];


        temp1 = a1 + b1;
        temp2 = a1 - b1;

        a2 = ((temp1 * x_c2) >> 16) + temp1;
        c2 = ((temp2 * x_c2) >> 16) + temp2;

        temp1 = (c1 * x_c3) >> 16;
        temp2 = ((d1 * x_c1) >> 16) + d1;

        b2 = temp1 + temp2;

        temp1 = (d1 * x_c3) >> 16;
        temp2 = ((c1 * x_c1) >> 16) + c1;

        d2 = temp1 - temp2;


        op[0]   = (a2 + 1) >> 1;
        op[4]   = (b2 + 1) >> 1;
        op[8]   = (c2 + 1) >> 1;
        op[12]  = (d2 + 1) >> 1;

        ip++;
        op++;
    }
}

void vp8_fast_fdct8x4_c(short *input, short *output, int pitch)
{
    vp8_fast_fdct4x4_c(input,   output,    pitch);
    vp8_fast_fdct4x4_c(input + 4, output + 16, pitch);
}

void vp8_short_walsh4x4_c(short *input, short *output, int pitch)
{
    int i;
    int a1, b1, c1, d1;
    int a2, b2, c2, d2;
    short *ip = input;
    short *op = output;

    for (i = 0; i < 4; i++)
    {
        a1 = ip[0] + ip[3];
        b1 = ip[1] + ip[2];
        c1 = ip[1] - ip[2];
        d1 = ip[0] - ip[3];

        op[0] = a1 + b1;
        op[1] = c1 + d1;
        op[2] = a1 - b1;
        op[3] = d1 - c1;
        ip += pitch / 2;
        op += 4;
    }

    ip = output;
    op = output;

    for (i = 0; i < 4; i++)
    {
        a1 = ip[0] + ip[12];
        b1 = ip[4] + ip[8];
        c1 = ip[4] - ip[8];
        d1 = ip[0] - ip[12];

        a2 = a1 + b1;
        b2 = c1 + d1;
        c2 = a1 - b1;
        d2 = d1 - c1;

        a2 += (a2 > 0);
        b2 += (b2 > 0);
        c2 += (c2 > 0);
        d2 += (d2 > 0);

        op[0] = (a2) >> 1;
        op[4] = (b2) >> 1;
        op[8] = (c2) >> 1;
        op[12] = (d2) >> 1;

        ip++;
        op++;
    }
}
