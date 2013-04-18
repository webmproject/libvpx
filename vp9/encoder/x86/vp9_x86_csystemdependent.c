/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */


#include "./vpx_config.h"
#include "vpx_ports/x86.h"
#include "vp9/encoder/vp9_variance.h"
#include "vp9/encoder/vp9_onyx_int.h"
#include "vp9/encoder/x86/vp9_dct_mmx.h"

// TODO(jimbankoski) Consider rewriting the c to take the same values rather
// than going through these pointer conversions
#if HAVE_MMX
void vp9_short_fdct8x4_mmx(short *input, short *output, int pitch) {
  vp9_short_fdct4x4_mmx(input,   output,    pitch);
  vp9_short_fdct4x4_mmx(input + 4, output + 16, pitch);
}

void vp9_subtract_b_mmx_impl(unsigned char *z,  int src_stride,
                             short *diff, unsigned char *predictor,
                             int pitch);
void vp9_subtract_b_mmx(BLOCK *be, BLOCKD *bd, int pitch) {
  unsigned char *z = *(be->base_src) + be->src;
  unsigned int  src_stride = be->src_stride;
  short *diff = &be->src_diff[0];
  unsigned char *predictor = *(bd->base_dst) + bd->dst;
  // TODO(jingning): The prototype function in c has been changed. Need to
  // modify the mmx and sse versions.
  vp9_subtract_b_mmx_impl(z, src_stride, diff, predictor, pitch);
}

#endif

#if HAVE_SSE2
void vp9_subtract_b_sse2_impl(unsigned char *z,  int src_stride,
                              short *diff, unsigned char *predictor,
                              int pitch);
void vp9_subtract_b_sse2(BLOCK *be, BLOCKD *bd, int pitch) {
  unsigned char *z = *(be->base_src) + be->src;
  unsigned int  src_stride = be->src_stride;
  short *diff = &be->src_diff[0];
  unsigned char *predictor = *(bd->base_dst) + bd->dst;
  // TODO(jingning): The prototype function in c has been changed. Need to
  // modify the mmx and sse versions.
  vp9_subtract_b_sse2_impl(z, src_stride, diff, predictor, pitch);
}

#endif
