/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */


/****************************************************************************
 * Notes:
 *
 * This implementation makes use of 16 bit fixed point verio of two multiply
 * constants:
 *         1.   sqrt(2) * cos (pi/8)
 *         2.   sqrt(2) * sin (pi/8)
 * Becuase the first constant is bigger than 1, to maintain the same 16 bit
 * fixed point precision as the second one, we use a trick of
 *         x * a = x + x*(a-1)
 * so
 *         x * sqrt(2) * cos (pi/8) = x + x * (sqrt(2) *cos(pi/8)-1).
 **************************************************************************/
#include "vpx_ports/config.h"


#include <math.h>

static const int cospi8sqrt2minus1 = 20091;
static const int sinpi8sqrt2      = 35468;
static const int rounding = 0;

void vp8_short_idct4x4llm_c(short *input, short *output, int pitch)
{
    int i;
    int a1, b1, c1, d1;

    short *ip = input;
    short *op = output;
    int temp1, temp2;
    int shortpitch = pitch >> 1;

    for (i = 0; i < 4; i++)
    {
        a1 = ip[0] + ip[8];
        b1 = ip[0] - ip[8];

        temp1 = (ip[4] * sinpi8sqrt2 + rounding) >> 16;
        temp2 = ip[12] + ((ip[12] * cospi8sqrt2minus1 + rounding) >> 16);
        c1 = temp1 - temp2;

        temp1 = ip[4] + ((ip[4] * cospi8sqrt2minus1 + rounding) >> 16);
        temp2 = (ip[12] * sinpi8sqrt2 + rounding) >> 16;
        d1 = temp1 + temp2;

        op[shortpitch*0] = a1 + d1;
        op[shortpitch*3] = a1 - d1;

        op[shortpitch*1] = b1 + c1;
        op[shortpitch*2] = b1 - c1;

        ip++;
        op++;
    }

    ip = output;
    op = output;

    for (i = 0; i < 4; i++)
    {
        a1 = ip[0] + ip[2];
        b1 = ip[0] - ip[2];

        temp1 = (ip[1] * sinpi8sqrt2 + rounding) >> 16;
        temp2 = ip[3] + ((ip[3] * cospi8sqrt2minus1 + rounding) >> 16);
        c1 = temp1 - temp2;

        temp1 = ip[1] + ((ip[1] * cospi8sqrt2minus1 + rounding) >> 16);
        temp2 = (ip[3] * sinpi8sqrt2 + rounding) >> 16;
        d1 = temp1 + temp2;

        op[0] = (a1 + d1 + 16) >> 5;
        op[3] = (a1 - d1 + 16) >> 5;

        op[1] = (b1 + c1 + 16) >> 5;
        op[2] = (b1 - c1 + 16) >> 5;

        ip += shortpitch;
        op += shortpitch;
    }
}

void vp8_short_idct4x4llm_1_c(short *input, short *output, int pitch)
{
    int i;
    int a1;
    short *op = output;
    int shortpitch = pitch >> 1;
    a1 = ((input[0] + 16) >> 5);
    for (i = 0; i < 4; i++)
    {
        op[0] = a1;
        op[1] = a1;
        op[2] = a1;
        op[3] = a1;
        op += shortpitch;
    }
}

void vp8_dc_only_idct_add_c(short input_dc, unsigned char *pred_ptr, unsigned char *dst_ptr, int pitch, int stride)
{
    int a1 = ((input_dc + 16) >> 5);
    int r, c;

    for (r = 0; r < 4; r++)
    {
        for (c = 0; c < 4; c++)
        {
            int a = a1 + pred_ptr[c] ;

            if (a < 0)
                a = 0;

            if (a > 255)
                a = 255;

            dst_ptr[c] = (unsigned char) a ;
        }

        dst_ptr += stride;
        pred_ptr += pitch;
    }

}

void vp8_short_inv_walsh4x4_c(short *input, short *output)
{
    int i;
    int a1, b1, c1, d1;
    int a2, b2, c2, d2;
    short *ip = input;
    short *op = output;

    for (i = 0; i < 4; i++)
    {
        a1 = ip[0] + ip[12];
        b1 = ip[4] + ip[8];
        c1 = ip[4] - ip[8];
        d1 = ip[0] - ip[12];

        op[0] = a1 + b1;
        op[4] = c1 + d1;
        op[8] = a1 - b1;
        op[12] = d1 - c1;
        ip++;
        op++;
    }

    ip = output;
    op = output;

    for (i = 0; i < 4; i++)
    {
        a1 = ip[0] + ip[3];
        b1 = ip[1] + ip[2];
        c1 = ip[1] - ip[2];
        d1 = ip[0] - ip[3];

        a2 = a1 + b1;
        b2 = c1 + d1;
        c2 = a1 - b1;
        d2 = d1 - c1;

        op[0] = (a2 + 1) >> 2;
        op[1] = (b2 + 1) >> 2;
        op[2] = (c2 + 1) >> 2;
        op[3] = (d2 + 1) >> 2;

        ip += 4;
        op += 4;
    }
}

void vp8_short_inv_walsh4x4_1_c(short *input, short *output)
{
    int i;
    int a1;
    short *op = output;

    a1 = (input[0] + 1 )>> 2;

    for (i = 0; i < 4; i++)
    {
        op[0] = a1;
        op[1] = a1;
        op[2] = a1;
        op[3] = a1;
        op += 4;
    }
}

#if CONFIG_T8X8

#define FAST_IDCT_8X8

void vp8_short_idct8x8_1_c(short *input, short *output, int pitch)
{
    int i, b;
    int a1;
    short *op = output;
    short *orig_op = output;
    int shortpitch = pitch >> 1;
    //a1 = ((input[0] + 4) >> 3);
    a1 = ((input[0] + 16) >> 5);
    for (b = 0; b < 4; b++)
    {
        for (i = 0; i < 4; i++)
        {
            op[0] = a1;
            op[1] = a1;
            op[2] = a1;
            op[3] = a1;
            op += shortpitch;
        }
        op = orig_op + (b+1)%2*4 +(b+1)/2*4*shortpitch;
    }
}

void vp8_dc_only_idct_add_8x8_c(short input_dc, unsigned char *pred_ptr, unsigned char *dst_ptr, int pitch, int stride)
{
    //int a1 = ((input_dc + 4) >> 3);
    int a1 = ((input_dc + 16) >> 5);
    int r, c, b;
    unsigned char *orig_pred = pred_ptr;
    unsigned char *orig_dst = dst_ptr;
    for (b = 0; b < 4; b++)
    {
        for (r = 0; r < 4; r++)
        {
          for (c = 0; c < 4; c++)
          {
              int a = a1 + pred_ptr[c] ;

              if (a < 0)
                 a = 0;

              if (a > 255)
                 a = 255;

              dst_ptr[c] = (unsigned char) a ;
         }

         dst_ptr += stride;
         pred_ptr += pitch;
       }
        dst_ptr = orig_dst + (b+1)%2*4 + (b+1)/2*4*stride;
        pred_ptr = orig_pred + (b+1)%2*4 + (b+1)/2*4*pitch;
    }
}

#ifdef FAST_IDCT_8X8

#define W1 2841                 /* 2048*sqrt(2)*cos(1*pi/16) */
#define W2 2676                 /* 2048*sqrt(2)*cos(2*pi/16) */
#define W3 2408                 /* 2048*sqrt(2)*cos(3*pi/16) */
#define W5 1609                 /* 2048*sqrt(2)*cos(5*pi/16) */
#define W6 1108                 /* 2048*sqrt(2)*cos(6*pi/16) */
#define W7 565                  /* 2048*sqrt(2)*cos(7*pi/16) */

/* row (horizontal) IDCT
 *
 * 7                       pi         1 dst[k] = sum c[l] * src[l] * cos( -- *
 * ( k + - ) * l ) l=0                      8          2
 *
 * where: c[0]    = 128 c[1..7] = 128*sqrt(2) */

static void idctrow (int *blk)
{
  int x0, x1, x2, x3, x4, x5, x6, x7, x8;

  /* shortcut */
  if (!((x1 = blk[4] << 11) | (x2 = blk[6]) | (x3 = blk[2]) |
        (x4 = blk[1]) | (x5 = blk[7]) | (x6 = blk[5]) | (x7 = blk[3])))
  {
    blk[0] = blk[1] = blk[2] = blk[3] = blk[4] = blk[5] = blk[6] = blk[7] = blk[0] << 3;
    return;
  }
  x0 = (blk[0] << 11) + 128;    /* for proper rounding in the fourth stage */

  /* first stage */
  x8 = W7 * (x4 + x5);
  x4 = x8 + (W1 - W7) * x4;
  x5 = x8 - (W1 + W7) * x5;
  x8 = W3 * (x6 + x7);
  x6 = x8 - (W3 - W5) * x6;
  x7 = x8 - (W3 + W5) * x7;

  /* second stage */
  x8 = x0 + x1;
  x0 -= x1;
  x1 = W6 * (x3 + x2);
  x2 = x1 - (W2 + W6) * x2;
  x3 = x1 + (W2 - W6) * x3;
  x1 = x4 + x6;
  x4 -= x6;
  x6 = x5 + x7;
  x5 -= x7;

  /* third stage */
  x7 = x8 + x3;
  x8 -= x3;
  x3 = x0 + x2;
  x0 -= x2;
  x2 = (181 * (x4 + x5) + 128) >> 8;
  x4 = (181 * (x4 - x5) + 128) >> 8;

  /* fourth stage */
  blk[0] = (x7 + x1) >> 8;
  blk[1] = (x3 + x2) >> 8;
  blk[2] = (x0 + x4) >> 8;
  blk[3] = (x8 + x6) >> 8;
  blk[4] = (x8 - x6) >> 8;
  blk[5] = (x0 - x4) >> 8;
  blk[6] = (x3 - x2) >> 8;
  blk[7] = (x7 - x1) >> 8;
}

/* column (vertical) IDCT
 *
 * 7                         pi         1 dst[8*k] = sum c[l] * src[8*l] *
 * cos( -- * ( k + - ) * l ) l=0                        8          2
 *
 * where: c[0]    = 1/1024 c[1..7] = (1/1024)*sqrt(2) */
static void idctcol (int *blk)
{
  int x0, x1, x2, x3, x4, x5, x6, x7, x8;

  /* shortcut */
  if (!((x1 = (blk[8 * 4] << 8)) | (x2 = blk[8 * 6]) | (x3 = blk[8 * 2]) |
        (x4 = blk[8 * 1]) | (x5 = blk[8 * 7]) | (x6 = blk[8 * 5]) | (x7 = blk[8 * 3])))
  {
    blk[8 * 0] = blk[8 * 1] = blk[8 * 2] = blk[8 * 3] = blk[8 * 4] = blk[8 * 5] = blk[8 * 6] = blk[8 * 7] =
      ((blk[8 * 0] + 32) >> 6);
    return;
  }
  x0 = (blk[8 * 0] << 8) + 16384;

  /* first stage */
  x8 = W7 * (x4 + x5) + 4;
  x4 = (x8 + (W1 - W7) * x4) >> 3;
  x5 = (x8 - (W1 + W7) * x5) >> 3;
  x8 = W3 * (x6 + x7) + 4;
  x6 = (x8 - (W3 - W5) * x6) >> 3;
  x7 = (x8 - (W3 + W5) * x7) >> 3;

  /* second stage */
  x8 = x0 + x1;
  x0 -= x1;
  x1 = W6 * (x3 + x2) + 4;
  x2 = (x1 - (W2 + W6) * x2) >> 3;
  x3 = (x1 + (W2 - W6) * x3) >> 3;
  x1 = x4 + x6;
  x4 -= x6;
  x6 = x5 + x7;
  x5 -= x7;

  /* third stage */
  x7 = x8 + x3;
  x8 -= x3;
  x3 = x0 + x2;
  x0 -= x2;
  x2 = (181 * (x4 + x5) + 128) >> 8;
  x4 = (181 * (x4 - x5) + 128) >> 8;

  /* fourth stage */
  blk[8 * 0] = (x7 + x1 ) >> 14;
  blk[8 * 1] = (x3 + x2 ) >> 14;
  blk[8 * 2] = (x0 + x4 ) >> 14;
  blk[8 * 3] = (x8 + x6 ) >> 14;
  blk[8 * 4] = (x8 - x6 ) >> 14;
  blk[8 * 5] = (x0 - x4 ) >> 14;
  blk[8 * 6] = (x3 - x2 ) >> 14;
  blk[8 * 7] = (x7 - x1 ) >> 14;
}

#define TX_DIM 8
void vp8_short_idct8x8_c(short *coefs, short *block, int pitch)
// an approximate 8x8 dct implementation, but not used
{
    int X[TX_DIM*TX_DIM];
    int i,j;
    int shortpitch = pitch >> 1;

    for (i = 0; i < TX_DIM; i++)
    {
        for (j = 0; j < TX_DIM; j++)
        {
             X[i * TX_DIM + j] = (int)(coefs[i * TX_DIM + j]+2)>>2;
        }
    }
  for (i = 0; i < 8; i++)
    idctrow (X + 8 * i);

  for (i = 0; i < 8; i++)
    idctcol (X + i);

    for (i = 0; i < TX_DIM; i++)
    {
        for (j = 0; j < TX_DIM; j++)
        {
             block[i*shortpitch+j]  = X[i * TX_DIM + j]>>1;
        }
    }
}

#else

/* This is really for testing */
void vp8_short_idct8x8_c(short *input, short *output, int pitch)
{
    int X[8][8];
    double C[8][8]={{0.0}}, Ct[8][8]={{0.0}}, temp[8][8]={{0.0}};
    int i,j,k;
    double temp1=0.0;
    double pi = atan( 1.0 ) * 4.0;
    //static int count=0;

    int shortpitch = pitch >> 1;

    for (i = 0; i < 8; i++)
    {
        for (j = 0; j < 8; j++)
        {
             X[i][j] = input[i * 8 + j];
        }
     }

    // TODO: DCT matrix should be calculated once for all
    for ( j = 0 ; j < 8 ; j++ ) {
        C[ 0 ][ j ] = 1.0 / sqrt( (double) 8 );
        Ct[ j ][ 0 ] = C[ 0 ][ j ];
    }
    for ( i = 1 ; i < 8 ; i++ ) {
        for ( j = 0 ; j < 8 ; j++ ) {
            C[ i ][ j ] = sqrt( 2.0 / 8 ) *
                          cos( pi * ( 2 * j + 1 ) * i / ( 2.0 * 8 ) );
            Ct[ j ][ i ] = C[ i ][ j ];
        }
    }
    /*  MatrixMultiply( temp, input, C ); */
        for ( i = 0 ; i < 8 ; i++ ) {
            for ( j = 0 ; j < 8 ; j++ ) {
                temp[ i ][ j ] = 0.0;
                for ( k = 0 ; k < 8 ; k++ )
                    temp[ i ][ j ] += X[ i ][ k ] * C[ k ][ j ];
            }
        }

    /*  MatrixMultiply( output, Ct, temp ); */
        for ( i = 0 ; i < 8 ; i++ ) {
            for ( j = 0 ; j < 8 ; j++ ) {
                temp1 = 0.0;
                for ( k = 0 ; k < 8 ; k++ )
                    temp1 += Ct[ i ][ k ] * temp[ k ][ j ];
                X[ i ][ j ] = floor( temp1/ 2.0 + 0.5);
            }
        }

    for (i = 0; i < 8; i++)
    {
        for (j = 0; j < 8; j++)
        {
             output[i*shortpitch+j]  = X[i][j];
        }
    }
}
#endif

void vp8_short_ihaar2x2_c(short *input, short *output, int pitch)
{
   int i, x;
   short *ip = input; //0,1, 4, 8
   short *op = output;
   for (i = 0; i < 16; i++)
   {
       op[i] = 0;
   }

   op[0] = (ip[0] + ip[1] + ip[4] + ip[8])>>1;
   op[1] = (ip[0] - ip[1] + ip[4] - ip[8])>>1;
   op[4] = (ip[0] + ip[1] - ip[4] - ip[8])>>1;
   op[8] = (ip[0] - ip[1] - ip[4] + ip[8])>>1;
}

void vp8_short_ihaar2x2_1_c(short *input, short *output, int pitch)
{
   int a1;
   short *ip = input;
   short *op = output;
   a1 = ip[0]>> 2;
   op[0] = a1;
   op[2] = a1;
   op[8] = a1;
   op[10] = a1;

}
#endif
