/*
 *  Copyright (c) 2013 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#ifndef VP9_COMMON_MIPS_DSPR2_VP9_COMMON_DSPR2_H_
#define VP9_COMMON_MIPS_DSPR2_VP9_COMMON_DSPR2_H_

#include <assert.h>

#include "./vpx_config.h"
#include "vpx/vpx_integer.h"
#include "vpx_dsp/mips/common_dspr2.h"

#ifdef __cplusplus
extern "C" {
#endif

#if HAVE_DSPR2

extern uint8_t *vpx_ff_cropTbl;

#define DCT_CONST_ROUND_SHIFT_TWICE_COSPI_16_64(input)                    ({   \
                                                                               \
  int32_t tmp, out;                                                            \
  int     dct_cost_rounding = DCT_CONST_ROUNDING;                              \
  int     in = input;                                                          \
                                                                               \
  __asm__ __volatile__ (                                                       \
      /* out = dct_const_round_shift(input_dc * cospi_16_64); */               \
      "mtlo     %[dct_cost_rounding],   $ac1                              \n\t"\
      "mthi     $zero,                  $ac1                              \n\t"\
      "madd     $ac1,                   %[in],            %[cospi_16_64]  \n\t"\
      "extp     %[tmp],                 $ac1,             31              \n\t"\
                                                                               \
      /* out = dct_const_round_shift(out * cospi_16_64); */                    \
      "mtlo     %[dct_cost_rounding],   $ac2                              \n\t"\
      "mthi     $zero,                  $ac2                              \n\t"\
      "madd     $ac2,                   %[tmp],           %[cospi_16_64]  \n\t"\
      "extp     %[out],                 $ac2,             31              \n\t"\
                                                                               \
      : [tmp] "=&r" (tmp), [out] "=r" (out)                                    \
      : [in] "r" (in),                                                         \
        [dct_cost_rounding] "r" (dct_cost_rounding),                           \
        [cospi_16_64] "r" (cospi_16_64)                                        \
   );                                                                          \
  out;                                                                    })

void vp9_idct32_cols_add_blk_dspr2(int16_t *input, uint8_t *dest,
                                   int dest_stride);

#endif  // #if HAVE_DSPR2
#ifdef __cplusplus
}  // extern "C"
#endif

#endif  // VP9_COMMON_MIPS_DSPR2_VP9_COMMON_DSPR2_H_
