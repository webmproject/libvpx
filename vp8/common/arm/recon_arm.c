/*
 *  Copyright (c) 2010 The VP8 project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license and patent
 *  grant that can be found in the LICENSE file in the root of the source
 *  tree. All contributing project authors may be found in the AUTHORS
 *  file in the root of the source tree.
 */


#include "vpx_ports/config.h"
#include "recon.h"
#include "blockd.h"

extern void vp8_recon16x16mb_neon(unsigned char *pred_ptr, short *diff_ptr, unsigned char *dst_ptr, int ystride, unsigned char *udst_ptr, unsigned char *vdst_ptr);

/*
void vp8_recon16x16mby(MACROBLOCKD *x)
{
    int i;
    for(i=0;i<16;i+=4)
    {
        //vp8_recon4b(&x->block[i]);
        BLOCKD *b = &x->block[i];
        vp8_recon4b (b->predictor, b->diff, *(b->base_dst) + b->dst, b->dst_stride);
    }
}
*/
void vp8_recon16x16mby(const vp8_recon_rtcd_vtable_t *rtcd, MACROBLOCKD *x)
{
    BLOCKD *b = &x->block[0];
    RECON_INVOKE(rtcd, recon4)(b->predictor, b->diff, *(b->base_dst) + b->dst, b->dst_stride);

    //b = &x->block[4];
    b += 4;
    RECON_INVOKE(rtcd, recon4)(b->predictor, b->diff, *(b->base_dst) + b->dst, b->dst_stride);

    //b = &x->block[8];
    b += 4;
    RECON_INVOKE(rtcd, recon4)(b->predictor, b->diff, *(b->base_dst) + b->dst, b->dst_stride);

    //b = &x->block[12];
    b += 4;
    RECON_INVOKE(rtcd, recon4)(b->predictor, b->diff, *(b->base_dst) + b->dst, b->dst_stride);
}

#if HAVE_ARMV7
void vp8_recon16x16mb(const vp8_recon_rtcd_vtable_t *rtcd, MACROBLOCKD *x)
{
    unsigned char *pred_ptr = &x->predictor[0];
    short *diff_ptr = &x->diff[0];
    unsigned char *dst_ptr = x->dst.y_buffer;
    unsigned char *udst_ptr = x->dst.u_buffer;
    unsigned char *vdst_ptr = x->dst.v_buffer;
    int ystride = x->dst.y_stride;
    //int uv_stride = x->dst.uv_stride;

    vp8_recon16x16mb_neon(pred_ptr, diff_ptr, dst_ptr, ystride, udst_ptr, vdst_ptr);
}

#else
/*
void vp8_recon16x16mb(MACROBLOCKD *x)
{
    int i;

    for(i=0;i<16;i+=4)
    {
//      vp8_recon4b(&x->block[i]);
        BLOCKD *b = &x->block[i];
        vp8_recon4b (b->predictor, b->diff, *(b->base_dst) + b->dst, b->dst_stride);

    }
    for(i=16;i<24;i+=2)
    {
//      vp8_recon2b(&x->block[i]);
        BLOCKD *b = &x->block[i];
        vp8_recon2b (b->predictor, b->diff, *(b->base_dst) + b->dst, b->dst_stride);
    }
}
*/
void vp8_recon16x16mb(const vp8_recon_rtcd_vtable_t *rtcd, MACROBLOCKD *x)
{
    BLOCKD *b = &x->block[0];

    RECON_INVOKE(rtcd, recon4)(b->predictor, b->diff, *(b->base_dst) + b->dst, b->dst_stride);
    b += 4;
    RECON_INVOKE(rtcd, recon4)(b->predictor, b->diff, *(b->base_dst) + b->dst, b->dst_stride);
    b += 4;
    RECON_INVOKE(rtcd, recon4)(b->predictor, b->diff, *(b->base_dst) + b->dst, b->dst_stride);
    b += 4;
    RECON_INVOKE(rtcd, recon4)(b->predictor, b->diff, *(b->base_dst) + b->dst, b->dst_stride);
    b += 4;

    //b = &x->block[16];

    RECON_INVOKE(rtcd, recon2)(b->predictor, b->diff, *(b->base_dst) + b->dst, b->dst_stride);
    b++;
    b++;
    RECON_INVOKE(rtcd, recon2)(b->predictor, b->diff, *(b->base_dst) + b->dst, b->dst_stride);
    b++;
    b++;
    RECON_INVOKE(rtcd, recon2)(b->predictor, b->diff, *(b->base_dst) + b->dst, b->dst_stride);
    b++;
    b++;
    RECON_INVOKE(rtcd, recon2)(b->predictor, b->diff, *(b->base_dst) + b->dst, b->dst_stride);
}
#endif
