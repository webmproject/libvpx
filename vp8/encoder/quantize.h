/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */


#ifndef __INC_QUANTIZE_H
#define __INC_QUANTIZE_H

#include "block.h"

#define prototype_quantize_block(sym) \
    void (sym)(BLOCK *b,BLOCKD *d)

#if ARCH_X86 || ARCH_X86_64
#include "x86/quantize_x86.h"
#endif

#if ARCH_ARM
#include "arm/quantize_arm.h"
#endif

#ifndef vp8_quantize_quantb
#define vp8_quantize_quantb vp8_regular_quantize_b
#endif
extern prototype_quantize_block(vp8_quantize_quantb);

#ifndef vp8_quantize_fastquantb
#define vp8_quantize_fastquantb vp8_fast_quantize_b_c
#endif
extern prototype_quantize_block(vp8_quantize_fastquantb);

typedef struct
{
    prototype_quantize_block(*quantb);
    prototype_quantize_block(*fastquantb);
} vp8_quantize_rtcd_vtable_t;

#if CONFIG_RUNTIME_CPU_DETECT
#define QUANTIZE_INVOKE(ctx,fn) (ctx)->fn
#else
#define QUANTIZE_INVOKE(ctx,fn) vp8_quantize_##fn
#endif

extern void vp8_strict_quantize_b(BLOCK *b,BLOCKD *d);

extern void vp8_quantize_mb(MACROBLOCK *x);
extern void vp8_quantize_mbuv(MACROBLOCK *x);
extern void vp8_quantize_mby(MACROBLOCK *x);

#endif
