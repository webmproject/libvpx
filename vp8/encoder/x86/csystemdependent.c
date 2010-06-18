/*
 *  Copyright (c) 2010 The VP8 project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */


#include "variance.h"
#include "onyx_int.h"

SADFunction *vp8_sad16x16;
SADFunction *vp8_sad16x8;
SADFunction *vp8_sad8x16;
SADFunction *vp8_sad8x8;
SADFunction *vp8_sad4x4;

variance_function *vp8_variance4x4;
variance_function *vp8_variance8x8;
variance_function *vp8_variance8x16;
variance_function *vp8_variance16x8;
variance_function *vp8_variance16x16;


variance_function *vp8_mse16x16;

sub_pixel_variance_function *vp8_sub_pixel_variance4x4;
sub_pixel_variance_function *vp8_sub_pixel_variance8x8;
sub_pixel_variance_function *vp8_sub_pixel_variance8x16;
sub_pixel_variance_function *vp8_sub_pixel_variance16x8;
sub_pixel_variance_function *vp8_sub_pixel_variance16x16;

int (*vp8_block_error)(short *, short *);
int (*vp8_mbblock_error)(MACROBLOCK *mb, int dc);
void (*vp8_subtract_mby)(short *diff, unsigned char *src, unsigned char *pred, int stride);

extern void vp8_subtract_mby_c(short *diff, unsigned char *src, unsigned char *pred, int stride);
extern void vp8_subtract_mby_mmx(short *diff, unsigned char *src, unsigned char *pred, int stride);

extern int vp8_block_error_c(short *, short *);
extern int vp8_mbblock_error_c(MACROBLOCK *x, int dc);

extern int vp8_block_error_mmx(short *, short *);
extern int vp8_mbblock_error_mmx(MACROBLOCK *x, int dc);

extern int vp8_block_error_xmm(short *, short *);
extern int vp8_mbblock_error_xmm(MACROBLOCK *x, int dc);



int (*vp8_mbuverror)(MACROBLOCK *mb);
unsigned int (*vp8_get_mb_ss)(short *);
void (*vp8_short_fdct4x4)(short *input, short *output, int pitch);
void (*vp8_short_fdct8x4)(short *input, short *output, int pitch);
void (*vp8_fast_fdct4x4)(short *input, short *output, int pitch);
void (*vp8_fast_fdct8x4)(short *input, short *output, int pitch);

void (*vp8_subtract_b)(BLOCK *be, BLOCKD *bd, int pitch);
void (*vp8_subtract_mbuv)(short *diff, unsigned char *usrc, unsigned char *vsrc, unsigned char *pred, int stride);
void (*vp8_fast_quantize_b)(BLOCK *b, BLOCKD *d);
unsigned int (*vp8_get16x16pred_error)(unsigned char *src_ptr, int src_stride, unsigned char *ref_ptr, int ref_stride);
unsigned int (*vp8_get8x8var)(unsigned char *src_ptr, int  source_stride, unsigned char *ref_ptr, int  recon_stride, unsigned int *SSE, int *Sum);
unsigned int (*vp8_get16x16var)(unsigned char *src_ptr, int  source_stride, unsigned char *ref_ptr, int  recon_stride, unsigned int *SSE, int *Sum);
unsigned int (*vp8_get4x4sse_cs)(unsigned char *src_ptr, int  source_stride, unsigned char *ref_ptr, int  recon_stride);

// c imports
extern int vp8_mbuverror_c(MACROBLOCK *mb);
extern unsigned int vp8_get8x8var_c(unsigned char *src_ptr, int  source_stride, unsigned char *ref_ptr, int  recon_stride, unsigned int *SSE, int *Sum);
extern void vp8_short_fdct4x4_c(short *input, short *output, int pitch);
extern void vp8_short_fdct8x4_c(short *input, short *output, int pitch);
extern void vp8_fast_fdct4x4_c(short *input, short *output, int pitch);
extern void vp8_fast_fdct8x4_c(short *input, short *output, int pitch);


extern void vp8_subtract_b_c(BLOCK *be, BLOCKD *bd, int pitch);
extern void vp8_subtract_mbuv_c(short *diff, unsigned char *usrc, unsigned char *vsrc, unsigned char *pred, int stride);
extern void vp8_fast_quantize_b_c(BLOCK *b, BLOCKD *d);

extern SADFunction vp8_sad16x16_c;
extern SADFunction vp8_sad16x8_c;
extern SADFunction vp8_sad8x16_c;
extern SADFunction vp8_sad8x8_c;
extern SADFunction vp8_sad4x4_c;

extern SADFunction vp8_sad16x16_wmt;
extern SADFunction vp8_sad16x8_wmt;
extern SADFunction vp8_sad8x16_wmt;
extern SADFunction vp8_sad8x8_wmt;
extern SADFunction vp8_sad4x4_wmt;

extern SADFunction vp8_sad16x16_mmx;
extern SADFunction vp8_sad16x8_mmx;
extern SADFunction vp8_sad8x16_mmx;
extern SADFunction vp8_sad8x8_mmx;
extern SADFunction vp8_sad4x4_mmx;

extern variance_function vp8_variance16x16_c;
extern variance_function vp8_variance8x16_c;
extern variance_function vp8_variance16x8_c;
extern variance_function vp8_variance8x8_c;
extern variance_function vp8_variance4x4_c;
extern variance_function vp8_mse16x16_c;

extern sub_pixel_variance_function vp8_sub_pixel_variance4x4_c;
extern sub_pixel_variance_function vp8_sub_pixel_variance8x8_c;
extern sub_pixel_variance_function vp8_sub_pixel_variance8x16_c;
extern sub_pixel_variance_function vp8_sub_pixel_variance16x8_c;
extern sub_pixel_variance_function vp8_sub_pixel_variance16x16_c;

extern unsigned int vp8_get_mb_ss_c(short *);
extern unsigned int vp8_get16x16pred_error_c(unsigned char *src_ptr, int src_stride, unsigned char *ref_ptr, int ref_stride);
extern unsigned int vp8_get8x8var_c(unsigned char *src_ptr, int  source_stride, unsigned char *ref_ptr, int  recon_stride, unsigned int *SSE, int *Sum);
extern unsigned int vp8_get16x16var_c(unsigned char *src_ptr, int  source_stride, unsigned char *ref_ptr, int  recon_stride, unsigned int *SSE, int *Sum);
extern unsigned int vp8_get4x4sse_cs_c(unsigned char *src_ptr, int  source_stride, unsigned char *ref_ptr, int  recon_stride);

// mmx imports
extern int vp8_mbuverror_mmx(MACROBLOCK *mb);
extern void vp8_fast_quantize_b_mmx(BLOCK *b, BLOCKD *d);
extern void vp8_subtract_b_mmx(BLOCK *be, BLOCKD *bd, int pitch);
extern void vp8_subtract_mbuv_mmx(short *diff, unsigned char *usrc, unsigned char *vsrc, unsigned char *pred, int stride);
extern void vp8_short_fdct4x4_mmx(short *input, short *output, int pitch);
extern void vp8_short_fdct8x4_mmx(short *input, short *output, int pitch);
extern void vp8_fast_fdct8x4_mmx(short *input, short *output, int pitch);
extern void vp8_fast_fdct4x4_mmx(short *input, short *output, int pitch);
extern variance_function vp8_variance4x4_mmx;
extern variance_function vp8_variance8x8_mmx;
extern variance_function vp8_variance8x16_mmx;
extern variance_function vp8_variance16x8_mmx;
extern variance_function vp8_variance16x16_mmx;

extern variance_function vp8_mse16x16_mmx;
extern sub_pixel_variance_function vp8_sub_pixel_variance4x4_mmx;
extern sub_pixel_variance_function vp8_sub_pixel_variance8x8_mmx;
extern sub_pixel_variance_function vp8_sub_pixel_variance8x16_mmx;
extern sub_pixel_variance_function vp8_sub_pixel_variance16x8_mmx;
extern sub_pixel_variance_function vp8_sub_pixel_variance16x16_mmx;

extern unsigned int vp8_get16x16pred_error_mmx(unsigned char *src_ptr, int src_stride, unsigned char *ref_ptr, int ref_stride);
extern unsigned int vp8_get_mb_ss_mmx(short *);
extern unsigned int vp8_get8x8var_mmx(unsigned char *src_ptr, int  source_stride, unsigned char *ref_ptr, int  recon_stride, unsigned int *SSE, int *Sum);
extern unsigned int vp8_get16x16var_mmx(unsigned char *src_ptr, int  source_stride, unsigned char *ref_ptr, int  recon_stride, unsigned int *SSE, int *Sum);
extern unsigned int vp8_get4x4sse_cs_mmx(unsigned char *src_ptr, int  source_stride, unsigned char *ref_ptr, int  recon_stride);


// wmt imports
extern int vp8_mbuverror_xmm(MACROBLOCK *mb);
extern void vp8_fast_quantize_b_sse(BLOCK *b, BLOCKD *d);
extern void vp8_fast_fdct8x4_wmt(short *input, short *output, int pitch);
extern variance_function vp8_variance4x4_wmt;
extern variance_function vp8_variance8x8_wmt;
extern variance_function vp8_variance8x16_wmt;
extern variance_function vp8_variance16x8_wmt;
extern variance_function vp8_variance16x16_wmt;

extern variance_function vp8_mse16x16_wmt;
extern sub_pixel_variance_function vp8_sub_pixel_variance4x4_wmt;
extern sub_pixel_variance_function vp8_sub_pixel_variance8x8_wmt;
extern sub_pixel_variance_function vp8_sub_pixel_variance8x16_wmt;
extern sub_pixel_variance_function vp8_sub_pixel_variance16x8_wmt;
extern sub_pixel_variance_function vp8_sub_pixel_variance16x16_wmt;
extern unsigned int vp8_get16x16pred_error_sse2(unsigned char *src_ptr, int src_stride, unsigned char *ref_ptr, int ref_stride);
extern unsigned int vp8_get_mb_ss_sse2(short *src_ptr);
extern unsigned int vp8_get8x8var_sse2(unsigned char *src_ptr, int  source_stride, unsigned char *ref_ptr, int  recon_stride, unsigned int *SSE, int *Sum);
extern unsigned int vp8_get16x16var_sse2(unsigned char *src_ptr, int  source_stride, unsigned char *ref_ptr, int  recon_stride, unsigned int *SSE, int *Sum);

extern void vpx_get_processor_flags(int *mmx_enabled, int *xmm_enabled, int *wmt_enabled);

void vp8_cmachine_specific_config(void)
{
    int mmx_enabled;
    int xmm_enabled;
    int wmt_enabled;

    vpx_get_processor_flags(&mmx_enabled, &xmm_enabled, &wmt_enabled);

    if (wmt_enabled)         // Willamette
    {
        // Willamette instruction set available:
        vp8_mbuverror                = vp8_mbuverror_xmm;
        vp8_fast_quantize_b            = vp8_fast_quantize_b_sse;
        vp8_short_fdct4x4             = vp8_short_fdct4x4_mmx;
        vp8_short_fdct8x4             = vp8_short_fdct8x4_mmx;
        vp8_fast_fdct4x4              = vp8_fast_fdct4x4_mmx;
        vp8_fast_fdct8x4              = vp8_fast_fdct8x4_wmt;
        vp8_subtract_b                = vp8_subtract_b_mmx;
        vp8_subtract_mbuv             = vp8_subtract_mbuv_mmx;
        vp8_variance4x4              = vp8_variance4x4_mmx;
        vp8_variance8x8              = vp8_variance8x8_mmx;
        vp8_variance8x16             = vp8_variance8x16_wmt;
        vp8_variance16x8             = vp8_variance16x8_wmt;
        vp8_variance16x16            = vp8_variance16x16_wmt;
        vp8_mse16x16                 = vp8_mse16x16_wmt;
        vp8_sub_pixel_variance4x4      = vp8_sub_pixel_variance4x4_wmt;
        vp8_sub_pixel_variance8x8      = vp8_sub_pixel_variance8x8_wmt;
        vp8_sub_pixel_variance8x16     = vp8_sub_pixel_variance8x16_wmt;
        vp8_sub_pixel_variance16x8     = vp8_sub_pixel_variance16x8_wmt;
        vp8_sub_pixel_variance16x16    = vp8_sub_pixel_variance16x16_wmt;
        vp8_get_mb_ss                  = vp8_get_mb_ss_sse2;
        vp8_get16x16pred_error        = vp8_get16x16pred_error_sse2;
        vp8_get8x8var                = vp8_get8x8var_sse2;
        vp8_get16x16var              = vp8_get16x16var_sse2;
        vp8_get4x4sse_cs             = vp8_get4x4sse_cs_mmx;
        vp8_sad16x16                 = vp8_sad16x16_wmt;
        vp8_sad16x8                  = vp8_sad16x8_wmt;
        vp8_sad8x16                  = vp8_sad8x16_wmt;
        vp8_sad8x8                   = vp8_sad8x8_wmt;
        vp8_sad4x4                   = vp8_sad4x4_wmt;
        vp8_block_error               = vp8_block_error_xmm;
        vp8_mbblock_error             = vp8_mbblock_error_xmm;
        vp8_subtract_mby              = vp8_subtract_mby_mmx;

    }
    else if (mmx_enabled)
    {
        // MMX instruction set available:
        vp8_mbuverror                = vp8_mbuverror_mmx;
        vp8_fast_quantize_b            = vp8_fast_quantize_b_mmx;
        vp8_short_fdct4x4             = vp8_short_fdct4x4_mmx;
        vp8_short_fdct8x4             = vp8_short_fdct8x4_mmx;
        vp8_fast_fdct4x4              = vp8_fast_fdct4x4_mmx;
        vp8_fast_fdct8x4              = vp8_fast_fdct8x4_mmx;
        vp8_subtract_b                = vp8_subtract_b_mmx;
        vp8_subtract_mbuv             = vp8_subtract_mbuv_mmx;
        vp8_variance4x4              = vp8_variance4x4_mmx;
        vp8_variance8x8              = vp8_variance8x8_mmx;
        vp8_variance8x16             = vp8_variance8x16_mmx;
        vp8_variance16x8             = vp8_variance16x8_mmx;
        vp8_variance16x16            = vp8_variance16x16_mmx;
        vp8_mse16x16                 = vp8_mse16x16_mmx;
        vp8_sub_pixel_variance4x4      = vp8_sub_pixel_variance4x4_mmx;
        vp8_sub_pixel_variance8x8      = vp8_sub_pixel_variance8x8_mmx;
        vp8_sub_pixel_variance8x16     = vp8_sub_pixel_variance8x16_mmx;
        vp8_sub_pixel_variance16x8     = vp8_sub_pixel_variance16x8_mmx;
        vp8_sub_pixel_variance16x16    = vp8_sub_pixel_variance16x16_mmx;
        vp8_get_mb_ss                  = vp8_get_mb_ss_mmx;
        vp8_get16x16pred_error        = vp8_get16x16pred_error_mmx;
        vp8_get8x8var                = vp8_get8x8var_mmx;
        vp8_get16x16var              = vp8_get16x16var_mmx;
        vp8_get4x4sse_cs             = vp8_get4x4sse_cs_mmx;
        vp8_sad16x16                 = vp8_sad16x16_mmx;
        vp8_sad16x8                  = vp8_sad16x8_mmx;
        vp8_sad8x16                  = vp8_sad8x16_mmx;
        vp8_sad8x8                   = vp8_sad8x8_mmx;
        vp8_sad4x4                   = vp8_sad4x4_mmx;
        vp8_block_error               = vp8_block_error_mmx;
        vp8_mbblock_error             = vp8_mbblock_error_mmx;
        vp8_subtract_mby              = vp8_subtract_mby_mmx;

    }
    else
    {
        // Pure C:
        vp8_mbuverror                = vp8_mbuverror_c;
        vp8_fast_quantize_b            = vp8_fast_quantize_b_c;
        vp8_short_fdct4x4             = vp8_short_fdct4x4_c;
        vp8_short_fdct8x4             = vp8_short_fdct8x4_c;
        vp8_fast_fdct4x4              = vp8_fast_fdct4x4_c;
        vp8_fast_fdct8x4              = vp8_fast_fdct8x4_c;
        vp8_subtract_b                = vp8_subtract_b_c;
        vp8_subtract_mbuv             = vp8_subtract_mbuv_c;
        vp8_variance4x4              = vp8_variance4x4_c;
        vp8_variance8x8              = vp8_variance8x8_c;
        vp8_variance8x16             = vp8_variance8x16_c;
        vp8_variance16x8             = vp8_variance16x8_c;
        vp8_variance16x16            = vp8_variance16x16_c;
        vp8_mse16x16                 = vp8_mse16x16_c;
        vp8_sub_pixel_variance4x4      = vp8_sub_pixel_variance4x4_c;
        vp8_sub_pixel_variance8x8      = vp8_sub_pixel_variance8x8_c;
        vp8_sub_pixel_variance8x16     = vp8_sub_pixel_variance8x16_c;
        vp8_sub_pixel_variance16x8     = vp8_sub_pixel_variance16x8_c;
        vp8_sub_pixel_variance16x16    = vp8_sub_pixel_variance16x16_c;
        vp8_get_mb_ss                  = vp8_get_mb_ss_c;
        vp8_get16x16pred_error        = vp8_get16x16pred_error_c;
        vp8_get8x8var                = vp8_get8x8var_c;
        vp8_get16x16var              = vp8_get16x16var_c;
        vp8_get4x4sse_cs             = vp8_get4x4sse_cs_c;
        vp8_sad16x16                 = vp8_sad16x16_c;
        vp8_sad16x8                  = vp8_sad16x8_c;
        vp8_sad8x16                  = vp8_sad8x16_c;
        vp8_sad8x8                   = vp8_sad8x8_c;
        vp8_sad4x4                   = vp8_sad4x4_c;
        vp8_block_error               = vp8_block_error_c;
        vp8_mbblock_error             = vp8_mbblock_error_c;
        vp8_subtract_mby              = vp8_subtract_mby_c;
    }

}
