/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */


#ifndef DEQUANTIZE_H
#define DEQUANTIZE_H
#include "vp8/common/blockd.h"

#define prototype_dequant_block(sym) \
    void sym(BLOCKD *x)

#define prototype_dequant_idct_add(sym) \
    void sym(short *input, short *dq, \
             unsigned char *pred, unsigned char *output, \
             int pitch, int stride)

#define prototype_dequant_dc_idct_add(sym) \
    void sym(short *input, short *dq, \
             unsigned char *pred, unsigned char *output, \
             int pitch, int stride, \
             int dc)

#define prototype_dequant_dc_idct_add_y_block(sym) \
    void sym(short *q, short *dq, \
             unsigned char *pre, unsigned char *dst, \
             int stride, char *eobs, short *dc)

#define prototype_dequant_idct_add_y_block(sym) \
    void sym(short *q, short *dq, \
             unsigned char *pre, unsigned char *dst, \
             int stride, char *eobs)

#define prototype_dequant_idct_add_uv_block(sym) \
    void sym(short *q, short *dq, \
             unsigned char *pre, unsigned char *dst_u, \
             unsigned char *dst_v, int stride, char *eobs)

#define prototype_dequant_dc_idct_add_y_block_8x8(sym) \
    void sym(short *q, short *dq, \
             unsigned char *pre, unsigned char *dst, \
             int stride, char *eobs, short *dc, MACROBLOCKD *xd)

#define prototype_dequant_idct_add_y_block_8x8(sym) \
    void sym(short *q, short *dq, \
             unsigned char *pre, unsigned char *dst, \
             int stride, char *eobs, MACROBLOCKD *xd)

#define prototype_dequant_idct_add_uv_block_8x8(sym) \
    void sym(short *q, short *dq, \
             unsigned char *pre, unsigned char *dst_u, \
             unsigned char *dst_v, int stride, char *eobs, \
             MACROBLOCKD *xd)

#if ARCH_X86 || ARCH_X86_64
#include "x86/dequantize_x86.h"
#endif

#if ARCH_ARM
#include "arm/dequantize_arm.h"
#endif

#ifndef vp8_dequant_block
#define vp8_dequant_block vp8_dequantize_b_c
#endif
extern prototype_dequant_block(vp8_dequant_block);

#ifndef vp8_dequant_idct_add
#define vp8_dequant_idct_add vp8_dequant_idct_add_c
#endif
extern prototype_dequant_idct_add(vp8_dequant_idct_add);

#ifndef vp8_dequant_dc_idct_add
#define vp8_dequant_dc_idct_add vp8_dequant_dc_idct_add_c
#endif
extern prototype_dequant_dc_idct_add(vp8_dequant_dc_idct_add);

#ifndef vp8_dequant_dc_idct_add_y_block
#define vp8_dequant_dc_idct_add_y_block vp8_dequant_dc_idct_add_y_block_c
#endif
extern prototype_dequant_dc_idct_add_y_block(vp8_dequant_dc_idct_add_y_block);

#ifndef vp8_dequant_idct_add_y_block
#define vp8_dequant_idct_add_y_block vp8_dequant_idct_add_y_block_c
#endif
extern prototype_dequant_idct_add_y_block(vp8_dequant_idct_add_y_block);

#ifndef vp8_dequant_idct_add_uv_block
#define vp8_dequant_idct_add_uv_block vp8_dequant_idct_add_uv_block_c
#endif
extern prototype_dequant_idct_add_uv_block(vp8_dequant_idct_add_uv_block);

#if CONFIG_LOSSLESS
extern prototype_dequant_idct_add(vp8_dequant_idct_add_lossless_c);
extern prototype_dequant_dc_idct_add(vp8_dequant_dc_idct_add_lossless_c);
extern prototype_dequant_dc_idct_add_y_block(vp8_dequant_dc_idct_add_y_block_lossless_c);
extern prototype_dequant_idct_add_y_block(vp8_dequant_idct_add_y_block_lossless_c);
extern prototype_dequant_idct_add_uv_block(vp8_dequant_idct_add_uv_block_lossless_c);
#endif

#ifndef vp8_dequant_block_2x2
#define vp8_dequant_block_2x2 vp8_dequantize_b_2x2_c
#endif
extern prototype_dequant_block(vp8_dequant_block_2x2);

#ifndef vp8_dequant_idct_add_8x8
#define vp8_dequant_idct_add_8x8 vp8_dequant_idct_add_8x8_c
#endif
extern prototype_dequant_idct_add(vp8_dequant_idct_add_8x8);

#ifndef vp8_dequant_dc_idct_add_8x8
#define vp8_dequant_dc_idct_add_8x8 vp8_dequant_dc_idct_add_8x8_c
#endif
extern prototype_dequant_dc_idct_add(vp8_dequant_dc_idct_add_8x8);

#ifndef vp8_dequant_dc_idct_add_y_block_8x8
#define vp8_dequant_dc_idct_add_y_block_8x8 vp8_dequant_dc_idct_add_y_block_8x8_c
#endif
extern prototype_dequant_dc_idct_add_y_block_8x8(vp8_dequant_dc_idct_add_y_block_8x8);

#ifndef vp8_dequant_idct_add_y_block_8x8
#define vp8_dequant_idct_add_y_block_8x8 vp8_dequant_idct_add_y_block_8x8_c
#endif
extern prototype_dequant_idct_add_y_block_8x8(vp8_dequant_idct_add_y_block_8x8);

#ifndef vp8_dequant_idct_add_uv_block_8x8
#define vp8_dequant_idct_add_uv_block_8x8 vp8_dequant_idct_add_uv_block_8x8_c
#endif
extern prototype_dequant_idct_add_uv_block_8x8(vp8_dequant_idct_add_uv_block_8x8);



typedef prototype_dequant_block((*vp8_dequant_block_fn_t));

typedef prototype_dequant_idct_add((*vp8_dequant_idct_add_fn_t));

typedef prototype_dequant_dc_idct_add((*vp8_dequant_dc_idct_add_fn_t));

typedef prototype_dequant_dc_idct_add_y_block((*vp8_dequant_dc_idct_add_y_block_fn_t));

typedef prototype_dequant_idct_add_y_block((*vp8_dequant_idct_add_y_block_fn_t));

typedef prototype_dequant_idct_add_uv_block((*vp8_dequant_idct_add_uv_block_fn_t));

typedef prototype_dequant_dc_idct_add_y_block_8x8((*vp8_dequant_dc_idct_add_y_block_fn_t_8x8));

typedef prototype_dequant_idct_add_y_block_8x8((*vp8_dequant_idct_add_y_block_fn_t_8x8));

typedef prototype_dequant_idct_add_uv_block_8x8((*vp8_dequant_idct_add_uv_block_fn_t_8x8));

typedef struct
{
    vp8_dequant_block_fn_t               block;
    vp8_dequant_idct_add_fn_t            idct_add;
    vp8_dequant_dc_idct_add_fn_t         dc_idct_add;
    vp8_dequant_dc_idct_add_y_block_fn_t dc_idct_add_y_block;
    vp8_dequant_idct_add_y_block_fn_t    idct_add_y_block;
    vp8_dequant_idct_add_uv_block_fn_t   idct_add_uv_block;
    vp8_dequant_block_fn_t               block_2x2;
    vp8_dequant_idct_add_fn_t            idct_add_8x8;
    vp8_dequant_dc_idct_add_fn_t         dc_idct_add_8x8;
    vp8_dequant_dc_idct_add_y_block_fn_t_8x8 dc_idct_add_y_block_8x8;
    vp8_dequant_idct_add_y_block_fn_t_8x8    idct_add_y_block_8x8;
    vp8_dequant_idct_add_uv_block_fn_t_8x8   idct_add_uv_block_8x8;
} vp8_dequant_rtcd_vtable_t;

#if CONFIG_RUNTIME_CPU_DETECT
#define DEQUANT_INVOKE(ctx,fn) (ctx)->fn
#else
#define DEQUANT_INVOKE(ctx,fn) vp8_dequant_##fn
#endif

#endif
