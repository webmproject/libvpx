/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */


#ifndef __INC_DCT_H
#define __INC_DCT_H

#define prototype_fdct(sym) void (sym)(short *input, short *output, int pitch)

#if ARCH_X86 || ARCH_X86_64
#include "x86/dct_x86.h"
#endif

#if ARCH_ARM
#include "arm/dct_arm.h"
#endif


#if CONFIG_HYBRIDTRANSFORM8X8 || CONFIG_HYBRIDTRANSFORM
void vp8_fht_c(short *input, short *output, int pitch,
               TX_TYPE tx_type, int tx_dim);
#endif

#if CONFIG_TX16X16
#ifndef vp8_fdct_short16x16
#define vp8_fdct_short16x16 vp8_short_fdct16x16_c
#endif
extern prototype_fdct(vp8_fdct_short16x16);
#endif

#ifndef vp8_fdct_short8x8
#define vp8_fdct_short8x8  vp8_short_fdct8x8_c
#endif
extern prototype_fdct(vp8_fdct_short8x8);

#ifndef vp8_fhaar_short2x2
#define vp8_fhaar_short2x2  vp8_short_fhaar2x2_c
#endif
extern prototype_fdct(vp8_fhaar_short2x2);


#ifndef vp8_fdct_short4x4
#define vp8_fdct_short4x4  vp8_short_fdct4x4_c
#endif
extern prototype_fdct(vp8_fdct_short4x4);

#ifndef vp8_fdct_short8x4
#define vp8_fdct_short8x4  vp8_short_fdct8x4_c
#endif
extern prototype_fdct(vp8_fdct_short8x4);

// There is no fast4x4 (for now)
#ifndef vp8_fdct_fast4x4
#define vp8_fdct_fast4x4  vp8_short_fdct4x4_c
#endif

#ifndef vp8_fdct_fast8x4
#define vp8_fdct_fast8x4  vp8_short_fdct8x4_c
#endif

#ifndef vp8_fdct_walsh_short4x4
#define vp8_fdct_walsh_short4x4  vp8_short_walsh4x4_c
#endif
extern prototype_fdct(vp8_fdct_walsh_short4x4);

#if CONFIG_LOSSLESS
extern prototype_fdct(vp8_short_walsh4x4_x8_c);
extern prototype_fdct(vp8_short_walsh8x4_x8_c);
extern prototype_fdct(vp8_short_walsh4x4_lossless_c);
#endif

typedef prototype_fdct(*vp8_fdct_fn_t);
typedef struct {
#if CONFIG_TX16X16
  vp8_fdct_fn_t    short16x16;
#endif
  vp8_fdct_fn_t    short8x8;
  vp8_fdct_fn_t    haar_short2x2;
  vp8_fdct_fn_t    short4x4;
  vp8_fdct_fn_t    short8x4;
  vp8_fdct_fn_t    fast4x4;
  vp8_fdct_fn_t    fast8x4;
  vp8_fdct_fn_t    walsh_short4x4;
} vp8_fdct_rtcd_vtable_t;

#if CONFIG_RUNTIME_CPU_DETECT
#define FDCT_INVOKE(ctx,fn) (ctx)->fn
#else
#define FDCT_INVOKE(ctx,fn) vp8_fdct_##fn
#endif

#endif
