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
#include "vpx_ports/x86.h"
#include "vp9/encoder/variance.h"
#include "vp9/encoder/onyx_int.h"


#if HAVE_MMX
void vp9_short_fdct8x4_mmx(short *input, short *output, int pitch) {
  vp9_short_fdct4x4_mmx(input,   output,    pitch);
  vp9_short_fdct4x4_mmx(input + 4, output + 16, pitch);
}

int vp9_mbblock_error_mmx_impl(short *coeff_ptr, short *dcoef_ptr, int dc);
int vp9_mbblock_error_mmx(MACROBLOCK *mb, int dc) {
  short *coeff_ptr =  mb->block[0].coeff;
  short *dcoef_ptr =  mb->e_mbd.block[0].dqcoeff;
  return vp9_mbblock_error_mmx_impl(coeff_ptr, dcoef_ptr, dc);
}

int vp9_mbuverror_mmx_impl(short *s_ptr, short *d_ptr);
int vp9_mbuverror_mmx(MACROBLOCK *mb) {
  short *s_ptr = &mb->coeff[256];
  short *d_ptr = &mb->e_mbd.dqcoeff[256];
  return vp9_mbuverror_mmx_impl(s_ptr, d_ptr);
}

void vp9_subtract_b_mmx_impl(unsigned char *z,  int src_stride,
                             short *diff, unsigned char *predictor,
                             int pitch);
void vp9_subtract_b_mmx(BLOCK *be, BLOCKD *bd, int pitch) {
  unsigned char *z = *(be->base_src) + be->src;
  unsigned int  src_stride = be->src_stride;
  short *diff = &be->src_diff[0];
  unsigned char *predictor = &bd->predictor[0];
  vp9_subtract_b_mmx_impl(z, src_stride, diff, predictor, pitch);
}

#endif

#if HAVE_SSE2
int vp9_mbblock_error_xmm_impl(short *coeff_ptr, short *dcoef_ptr, int dc);
int vp9_mbblock_error_xmm(MACROBLOCK *mb, int dc) {
  short *coeff_ptr =  mb->block[0].coeff;
  short *dcoef_ptr =  mb->e_mbd.block[0].dqcoeff;
  return vp9_mbblock_error_xmm_impl(coeff_ptr, dcoef_ptr, dc);
}

int vp9_mbuverror_xmm_impl(short *s_ptr, short *d_ptr);
int vp9_mbuverror_xmm(MACROBLOCK *mb) {
  short *s_ptr = &mb->coeff[256];
  short *d_ptr = &mb->e_mbd.dqcoeff[256];
  return vp9_mbuverror_xmm_impl(s_ptr, d_ptr);
}

void vp9_subtract_b_sse2_impl(unsigned char *z,  int src_stride,
                              short *diff, unsigned char *predictor,
                              int pitch);
void vp9_subtract_b_sse2(BLOCK *be, BLOCKD *bd, int pitch) {
  unsigned char *z = *(be->base_src) + be->src;
  unsigned int  src_stride = be->src_stride;
  short *diff = &be->src_diff[0];
  unsigned char *predictor = &bd->predictor[0];
  vp9_subtract_b_sse2_impl(z, src_stride, diff, predictor, pitch);
}

#endif

void vp9_arch_x86_encoder_init(VP9_COMP *cpi) {
#if CONFIG_RUNTIME_CPU_DETECT
  int flags = x86_simd_caps();

  /* Note:
   *
   * This platform can be built without runtime CPU detection as well. If
   * you modify any of the function mappings present in this file, be sure
   * to also update them in static mapings (<arch>/filename_<arch>.h)
   */

  /* Override default functions with fastest ones for this CPU. */
#if HAVE_SSE2
  if (flags & HAS_SSE2) {
    cpi->rtcd.temporal.apply                 = vp9_temporal_filter_apply_sse2;

  }
#endif


#endif
}
