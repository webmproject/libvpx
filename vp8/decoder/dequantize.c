/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */


#include "vpx_ports/config.h"
#include "dequantize.h"
#include "vp8/common/idct.h"
#include "vpx_mem/vpx_mem.h"
#include "onyxd_int.h"

extern void vp8_short_idct4x4llm_c(short *input, short *output, int pitch) ;
extern void vp8_short_idct4x4llm_1_c(short *input, short *output, int pitch);
#if CONFIG_T8X8
extern void vp8_short_idct8x8_c(short *input, short *output, int pitch);
extern void vp8_short_idct8x8_1_c(short *input, short *output, int pitch);
#endif

#ifdef DEC_DEBUG
extern int dec_debug;
#endif

void vp8_dequantize_b_c(BLOCKD *d)
{

    int i;
    short *DQ  = d->dqcoeff;
    short *Q   = d->qcoeff;
    short *DQC = d->dequant;

    for (i = 0; i < 16; i++)
    {
        DQ[i] = Q[i] * DQC[i];
    }
}

void vp8_dequant_idct_add_c(short *input, short *dq, unsigned char *pred,
                            unsigned char *dest, int pitch, int stride)
{
    short output[16];
    short *diff_ptr = output;
    int r, c;
    int i;

    for (i = 0; i < 16; i++)
    {
        input[i] = dq[i] * input[i];
    }

    /* the idct halves ( >> 1) the pitch */
    vp8_short_idct4x4llm_c(input, output, 4 << 1);

    vpx_memset(input, 0, 32);

    for (r = 0; r < 4; r++)
    {
        for (c = 0; c < 4; c++)
        {
            int a = diff_ptr[c] + pred[c];

            if (a < 0)
                a = 0;

            if (a > 255)
                a = 255;

            dest[c] = (unsigned char) a;
        }

        dest += stride;
        diff_ptr += 4;
        pred += pitch;
    }
}

void vp8_dequant_dc_idct_add_c(short *input, short *dq, unsigned char *pred,
                               unsigned char *dest, int pitch, int stride,
                               int Dc)
{
    int i;
    short output[16];
    short *diff_ptr = output;
    int r, c;

    input[0] = (short)Dc;

    for (i = 1; i < 16; i++)
    {
        input[i] = dq[i] * input[i];
    }

    /* the idct halves ( >> 1) the pitch */
    vp8_short_idct4x4llm_c(input, output, 4 << 1);

    vpx_memset(input, 0, 32);

    for (r = 0; r < 4; r++)
    {
        for (c = 0; c < 4; c++)
        {
            int a = diff_ptr[c] + pred[c];

            if (a < 0)
                a = 0;

            if (a > 255)
                a = 255;

            dest[c] = (unsigned char) a;
        }

        dest += stride;
        diff_ptr += 4;
        pred += pitch;
    }
}

#if CONFIG_T8X8
void vp8_dequantize_b_2x2_c(BLOCKD *d)
{
    int i;
    short *DQ  = d->dqcoeff;
    short *Q   = d->qcoeff;
    short *DQC = d->dequant;

    for (i = 0; i < 16; i++)
    {
#if CONFIG_EXTEND_QRANGE
        DQ[i] = (short)((Q[i] * DQC[i]+2)>>2);
#else
       DQ[i] = (short)(Q[i] * DQC[i]);
#endif
    }
#ifdef DEC_DEBUG
    if (dec_debug) {
      int j;
      printf("Dequantize 2x2\n");
      for (j=0;j<16;j++) printf("%d ", Q[j]); printf("\n");
      for (j=0;j<16;j++) printf("%d ", DQ[j]); printf("\n");
    }
#endif
}

void vp8_dequant_idct_add_8x8_c(short *input, short *dq, unsigned char *pred,
                                unsigned char *dest, int pitch, int stride)//, MACROBLOCKD *xd, short blk_idx
{
    short output[64];
    short *diff_ptr = output;
    int r, c, b;
    int i;
    unsigned char *origdest = dest;
    unsigned char *origpred = pred;

#ifdef DEC_DEBUG
    if (dec_debug) {
      int j;
      printf("Input 8x8\n");
      for (j=0;j<64;j++) {
        printf("%d ", input[j]);
        if (j%8 == 7) printf("\n");
      }
    }
#endif

#if CONFIG_EXTEND_QRANGE
      input[0]= (input[0] * dq[0]+2)>>2;
#else
      input[0]= input[0] * dq[0];
#endif

    // recover quantizer for 4 4x4 blocks
    for (i = 1; i < 64; i++)
    {
#if CONFIG_EXTEND_QRANGE
      input[i]=(input[i] * dq[1]+2)>>2;
#else
      input[i]=input[i] * dq[1];
#endif
    }
#ifdef DEC_DEBUG
    if (dec_debug) {
      int j;
      printf("Input DQ 8x8\n");
      for (j=0;j<64;j++) {
        printf("%d ", input[j]);
        if (j%8 == 7) printf("\n");
      }
    }
#endif

    // the idct halves ( >> 1) the pitch
    vp8_short_idct8x8_c(input, output, 16);
#ifdef DEC_DEBUG
    if (dec_debug) {
      int j;
      printf("Output 8x8\n");
      for (j=0;j<64;j++) {
        printf("%d ", output[j]);
        if (j%8 == 7) printf("\n");
      }
    }
#endif

    vpx_memset(input, 0, 128);// test what should i put here

    for (b = 0; b < 4; b++)
    {
      for (r = 0; r < 4; r++)
      {
          for (c = 0; c < 4; c++)
          {
              int a = diff_ptr[c] + pred[c];

              if (a < 0)
                  a = 0;

              if (a > 255)
                  a = 255;

              dest[c] = (unsigned char) a;
          }

          dest += stride;
          diff_ptr += 8;
          pred += pitch;
      }
      diff_ptr = output + (b+1) / 2 * 4 * 8 + (b+1) % 2 * 4;
      dest = origdest + (b+1) / 2 * 4 * stride + (b+1) % 2 * 4;
      pred = origpred + (b+1) / 2 * 4 * pitch + (b+1) % 2 * 4;
   }
#ifdef DEC_DEBUG
    if (dec_debug) {
      int k,j;
      printf("Final 8x8\n");
      for (j=0;j<8;j++) {
        for (k=0;k<8;k++) {
          printf("%d ", origdest[k]);
        }
        printf("\n");
        origdest+=stride;
      }
    }
#endif
}

void vp8_dequant_dc_idct_add_8x8_c(short *input, short *dq, unsigned char *pred,
                               unsigned char *dest, int pitch, int stride,
                               int Dc)// Dc for 1st order T in some rear case
{
    short output[64];
    short *diff_ptr = output;
    int r, c, b;
    int i;
    unsigned char *origdest = dest;
    unsigned char *origpred = pred;

    input[0] = (short)Dc;//Dc is the reconstructed value, do not need dequantization
    //dc value is recovered after dequantization, since dc need not quantization
#ifdef DEC_DEBUG
    if (dec_debug) {
      int j;
      printf("Input 8x8\n");
      for (j=0;j<64;j++) {
        printf("%d ", input[j]);
        if (j%8 == 7) printf("\n");
      }
    }
#endif
    for (i = 1; i < 64; i++)
    {
#if CONFIG_EXTEND_QRANGE
        input[i]=(input[i] * dq[1]+2)>>2;
#else
        input[i]=input[i] * dq[1];
#endif
    }

#ifdef DEC_DEBUG
    if (dec_debug) {
      int j;
      printf("Input DQ 8x8\n");
      for (j=0;j<64;j++) {
        printf("%d ", input[j]);
        if (j%8 == 7) printf("\n");
      }
    }
#endif

    // the idct halves ( >> 1) the pitch
    vp8_short_idct8x8_c(input, output,16);
#ifdef DEC_DEBUG
    if (dec_debug) {
      int j;
      printf("Output 8x8\n");
      for (j=0;j<64;j++) {
        printf("%d ", output[j]);
        if (j%8 == 7) printf("\n");
      }
    }
#endif
    vpx_memset(input, 0, 128);

    for (b = 0; b < 4; b++)
    {
      for (r = 0; r < 4; r++)
      {
          for (c = 0; c < 4; c++)
          {
              int a = diff_ptr[c] + pred[c];

              if (a < 0)
                  a = 0;

              if (a > 255)
                  a = 255;

              dest[c] = (unsigned char) a;
          }

          dest += stride;
          diff_ptr += 8;
          pred += pitch;
      }
      diff_ptr = output + (b+1) / 2 * 4 * 8 + (b+1) % 2 * 4;
      dest = origdest + (b+1) / 2 * 4 * stride + (b+1) % 2 * 4;
      pred = origpred + (b+1) / 2 * 4 * pitch + (b+1) % 2 * 4;
   }
#ifdef DEC_DEBUG
    if (dec_debug) {
      int k,j;
      printf("Final 8x8\n");
      for (j=0;j<8;j++) {
        for (k=0;k<8;k++) {
          printf("%d ", origdest[k]);
        }
        printf("\n");
        origdest+=stride;
      }
    }
#endif
}

#endif