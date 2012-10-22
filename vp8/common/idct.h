/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */


#ifndef __INC_IDCT_H
#define __INC_IDCT_H

#define prototype_second_order(sym) \
  void sym(short *input, short *output)

#define prototype_idct(sym) \
  void sym(short *input, short *output, int pitch)

#define prototype_idct_scalar_add(sym) \
  void sym(short input, \
           unsigned char *pred, unsigned char *output, \
           int pitch, int stride)

#if ARCH_X86 || ARCH_X86_64
#include "x86/idct_x86.h"
#endif

#ifdef _MSC_VER
/* TODO: remove these after integer implmementations are done */
#define M_PI       3.14159265358979323846
#define round(x) (((x)>0)? floor((x)+0.5): ceil((x)-0.5))
#endif


#if ARCH_ARM
#include "arm/idct_arm.h"
#endif

#if CONFIG_LOSSLESS
#define WHT_UPSCALE_FACTOR 3
#define Y2_WHT_UPSCALE_FACTOR 2
#endif

#ifndef vp8_idct_idct16x16
#define vp8_idct_idct16x16 vp8_short_idct16x16_c
#endif
extern prototype_idct(vp8_idct_idct16x16);

#ifndef vp8_idct_idct8
#define vp8_idct_idct8 vp8_short_idct8x8_c
#endif
extern prototype_idct(vp8_idct_idct8);

#ifndef vp8_idct_idct8_1
#define vp8_idct_idct8_1 vp8_short_idct8x8_1_c
#endif
extern prototype_idct(vp8_idct_idct8_1);

#ifndef vp8_idct_ihaar2
#define vp8_idct_ihaar2 vp8_short_ihaar2x2_c
#endif
extern prototype_idct(vp8_idct_ihaar2);

#ifndef vp8_idct_ihaar2_1
#define vp8_idct_ihaar2_1 vp8_short_ihaar2x2_1_c
#endif
extern prototype_idct(vp8_idct_ihaar2_1);

#ifndef vp8_idct_idct1_scalar_add_8x8
#define vp8_idct_idct1_scalar_add_8x8 vp8_dc_only_idct_add_8x8_c
#endif
extern prototype_idct_scalar_add(vp8_idct_idct1_scalar_add_8x8);



#ifndef vp8_idct_idct1
#define vp8_idct_idct1 vp8_short_idct4x4llm_1_c
#endif
extern prototype_idct(vp8_idct_idct1);

#ifndef vp8_idct_idct16
#define vp8_idct_idct16 vp8_short_idct4x4llm_c
#endif
extern prototype_idct(vp8_idct_idct16);

#ifndef vp8_idct_idct1_scalar_add
#define vp8_idct_idct1_scalar_add vp8_dc_only_idct_add_c
#endif
extern prototype_idct_scalar_add(vp8_idct_idct1_scalar_add);


#ifndef vp8_idct_iwalsh1
#define vp8_idct_iwalsh1 vp8_short_inv_walsh4x4_1_c
#endif
extern prototype_second_order(vp8_idct_iwalsh1);

#ifndef vp8_idct_iwalsh16
#define vp8_idct_iwalsh16 vp8_short_inv_walsh4x4_c
#endif
extern prototype_second_order(vp8_idct_iwalsh16);

#if CONFIG_LOSSLESS
extern prototype_idct(vp8_short_inv_walsh4x4_x8_c);
extern prototype_idct(vp8_short_inv_walsh4x4_1_x8_c);
extern prototype_idct_scalar_add(vp8_dc_only_inv_walsh_add_c);
extern prototype_second_order(vp8_short_inv_walsh4x4_lossless_c);
extern prototype_second_order(vp8_short_inv_walsh4x4_1_lossless_c);
#endif

#include "vp8/common/blockd.h"
void vp8_ihtllm_c(short *input, short *output, int pitch,
                  TX_TYPE tx_type, int tx_dim);

typedef prototype_idct((*vp8_idct_fn_t));
typedef prototype_idct_scalar_add((*vp8_idct_scalar_add_fn_t));
typedef prototype_second_order((*vp8_second_order_fn_t));

typedef struct {
  vp8_idct_fn_t            idct1;
  vp8_idct_fn_t            idct16;
  vp8_idct_scalar_add_fn_t idct1_scalar_add;

  vp8_second_order_fn_t iwalsh1;
  vp8_second_order_fn_t iwalsh16;

  vp8_idct_fn_t            idct8;
  vp8_idct_fn_t            idct8_1;
  vp8_idct_scalar_add_fn_t idct1_scalar_add_8x8;
  vp8_idct_fn_t ihaar2;
  vp8_idct_fn_t ihaar2_1;

  vp8_idct_fn_t            idct16x16;
} vp8_idct_rtcd_vtable_t;

#if CONFIG_RUNTIME_CPU_DETECT
#define IDCT_INVOKE(ctx,fn) (ctx)->fn
#else
#define IDCT_INVOKE(ctx,fn) vp8_idct_##fn
#endif

#endif
