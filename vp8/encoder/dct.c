/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */


#include <math.h>
#include "vpx_ports/config.h"






void vp8_short_fdct8x8_c(short *block, short *coefs, int pitch)
{
  int j1, i, j, k;
  float b[8];
  float b1[8];
  float d[8][8];
  float f0 = (float) .7071068;
  float f1 = (float) .4903926;
  float f2 = (float) .4619398;
  float f3 = (float) .4157348;
  float f4 = (float) .3535534;
  float f5 = (float) .2777851;
  float f6 = (float) .1913417;
  float f7 = (float) .0975452;
  pitch = pitch / 2;
  for (i = 0, k = 0; i < 8; i++, k += pitch)
  {
    for (j = 0; j < 8; j++)
    {
      b[j] = (float)( block[k + j]<<1);
    }
    /* Horizontal transform */
    for (j = 0; j < 4; j++)
    {
      j1 = 7 - j;
      b1[j] = b[j] + b[j1];
      b1[j1] = b[j] - b[j1];
    }
    b[0] = b1[0] + b1[3];
    b[1] = b1[1] + b1[2];
    b[2] = b1[1] - b1[2];
    b[3] = b1[0] - b1[3];
    b[4] = b1[4];
    b[5] = (b1[6] - b1[5]) * f0;
    b[6] = (b1[6] + b1[5]) * f0;
    b[7] = b1[7];
    d[i][0] = (b[0] + b[1]) * f4;
    d[i][4] = (b[0] - b[1]) * f4;
    d[i][2] = b[2] * f6 + b[3] * f2;
    d[i][6] = b[3] * f6 - b[2] * f2;
    b1[4] = b[4] + b[5];
    b1[7] = b[7] + b[6];
    b1[5] = b[4] - b[5];
    b1[6] = b[7] - b[6];
    d[i][1] = b1[4] * f7 + b1[7] * f1;
    d[i][5] = b1[5] * f3 + b1[6] * f5;
    d[i][7] = b1[7] * f7 - b1[4] * f1;
    d[i][3] = b1[6] * f3 - b1[5] * f5;
  }
  /* Vertical transform */
  for (i = 0; i < 8; i++)
  {
    for (j = 0; j < 4; j++)
    {
      j1 = 7 - j;
      b1[j] = d[j][i] + d[j1][i];
      b1[j1] = d[j][i] - d[j1][i];
    }
    b[0] = b1[0] + b1[3];
    b[1] = b1[1] + b1[2];
    b[2] = b1[1] - b1[2];
    b[3] = b1[0] - b1[3];
    b[4] = b1[4];
    b[5] = (b1[6] - b1[5]) * f0;
    b[6] = (b1[6] + b1[5]) * f0;
    b[7] = b1[7];
    d[0][i] = (b[0] + b[1]) * f4;
    d[4][i] = (b[0] - b[1]) * f4;
    d[2][i] = b[2] * f6 + b[3] * f2;
    d[6][i] = b[3] * f6 - b[2] * f2;
    b1[4] = b[4] + b[5];
    b1[7] = b[7] + b[6];
    b1[5] = b[4] - b[5];
    b1[6] = b[7] - b[6];
    d[1][i] = b1[4] * f7 + b1[7] * f1;
    d[5][i] = b1[5] * f3 + b1[6] * f5;
    d[7][i] = b1[7] * f7 - b1[4] * f1;
    d[3][i] = b1[6] * f3 - b1[5] * f5;
  }
  for (i = 0; i < 8; i++)
  {
    for (j = 0; j < 8; j++)
    {
      *(coefs + j + i * 8) = (short) floor(d[i][j] +0.5);
    }
  }
  return;
}



void vp8_short_fhaar2x2_c(short *input, short *output, int pitch) //pitch = 8
{
    /* [1 1 ; 1 -1] orthogonal transform */
    /* use position: 0,1, 4, 8 */
   int i;
   short *ip1 = input;
   short *op1 = output;
   for (i = 0; i < 16; i++)
   {
       op1[i] = 0;
   }

   op1[0]=ip1[0] + ip1[1] + ip1[4] + ip1[8];
   op1[1]=ip1[0] - ip1[1] + ip1[4] - ip1[8];
   op1[4]=ip1[0] + ip1[1] - ip1[4] - ip1[8];
   op1[8]=ip1[0] - ip1[1] - ip1[4] + ip1[8];

}
void vp8_short_fdct4x4_c(short *input, short *output, int pitch)
{
    int i;
    int a1, b1, c1, d1;
    short *ip = input;
    short *op = output;

    for (i = 0; i < 4; i++)
    {
#if CONFIG_EXTEND_QRANGE
        a1 = ((ip[0] + ip[3])<<5);
        b1 = ((ip[1] + ip[2])<<5);
        c1 = ((ip[1] - ip[2])<<5);
        d1 = ((ip[0] - ip[3])<<5);
#else
        a1 = ((ip[0] + ip[3])<<3);
        b1 = ((ip[1] + ip[2])<<3);
        c1 = ((ip[1] - ip[2])<<3);
        d1 = ((ip[0] - ip[3])<<3);
#endif
        op[0] = a1 + b1;
        op[2] = a1 - b1;

        op[1] = (c1 * 2217 + d1 * 5352 +  14500)>>12;
        op[3] = (d1 * 2217 - c1 * 5352 +   7500)>>12;

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

        op[0]  = ( a1 + b1 + 7)>>4;
        op[8]  = ( a1 - b1 + 7)>>4;

        op[4]  =((c1 * 2217 + d1 * 5352 +  12000)>>16) + (d1!=0);
        op[12] = (d1 * 2217 - c1 * 5352 +  51000)>>16;

        ip++;
        op++;
    }
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
#if !CONFIG_EXTEND_QRANGE
        a1 = ((ip[0] + ip[2])<<2);
        d1 = ((ip[1] + ip[3])<<2);
        c1 = ((ip[1] - ip[3])<<2);
        b1 = ((ip[0] - ip[2])<<2);

        op[0] = a1 + d1+ (a1!=0);
#else
        a1 = ((ip[0] + ip[2]));
        d1 = ((ip[1] + ip[3]));
        c1 = ((ip[1] - ip[3]));
        b1 = ((ip[0] - ip[2]));


        op[0] = a1 + d1;
#endif
        op[1] = b1 + c1;
        op[2] = b1 - c1;
        op[3] = a1 - d1;
        ip += pitch / 2;
        op += 4;
    }

    ip = output;
    op = output;

    for (i = 0; i < 4; i++)
    {
        a1 = ip[0] + ip[8];
        d1 = ip[4] + ip[12];
        c1 = ip[4] - ip[12];
        b1 = ip[0] - ip[8];

        a2 = a1 + d1;
        b2 = b1 + c1;
        c2 = b1 - c1;
        d2 = a1 - d1;

        a2 += a2<0;
        b2 += b2<0;
        c2 += c2<0;
        d2 += d2<0;

#if !CONFIG_EXTEND_QRANGE
        op[0] = (a2+3) >> 3;
        op[4] = (b2+3) >> 3;
        op[8] = (c2+3) >> 3;
        op[12]= (d2+3) >> 3;
#else
        op[0] = (a2+1) >> 2;
        op[4] = (b2+1) >> 2;
        op[8] = (c2+1) >> 2;
        op[12]= (d2+1) >> 2;
#endif
        ip++;
        op++;
    }
}

void vp8_short_fdct8x4_c(short *input, short *output, int pitch)
{
    vp8_short_fdct4x4_c(input,   output,    pitch);
    vp8_short_fdct4x4_c(input + 4, output + 16, pitch);
}
