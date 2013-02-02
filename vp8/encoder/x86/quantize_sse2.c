/*
 *  Copyright (c) 2012 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */


#include "vp8/common/blockd.h"
#include "vp8/common/entropy.h"
#include "vp8/encoder/block.h"

#include <mmintrin.h> //MMX
#include <xmmintrin.h> //SSE
#include <emmintrin.h> //SSE2

void vp8_fast_quantize_b_sse2(BLOCK *b, BLOCKD *d)
{
  __m128i z0 = _mm_load_si128((__m128i *)(b->coeff));
  __m128i z1 = _mm_load_si128((__m128i *)(b->coeff + 8));
  __m128i round0 = _mm_load_si128((__m128i *)(b->round));
  __m128i round1 = _mm_load_si128((__m128i *)(b->round + 8));
  __m128i quant_fast0 = _mm_load_si128((__m128i *)(b->quant_fast));
  __m128i quant_fast1 = _mm_load_si128((__m128i *)(b->quant_fast + 8));
  __m128i dequant0 = _mm_load_si128((__m128i *)(d->dequant));
  __m128i dequant1 = _mm_load_si128((__m128i *)(d->dequant + 8));
  __m128i inv_zig_zag0 = _mm_load_si128((const __m128i *)(vp8_default_inv_zig_zag));
  __m128i inv_zig_zag1 = _mm_load_si128((const __m128i *)(vp8_default_inv_zig_zag + 8));

  __m128i sz0, sz1, x0, x1, y0, y1, xdq0, xdq1, zeros, ones;

  /* sign of z: z >> 15 */
  sz0 = _mm_srai_epi16(z0, 15);
  sz1 = _mm_srai_epi16(z1, 15);

  /* x = abs(z): (z ^ sz) - sz */
  x0 = _mm_xor_si128(z0, sz0);
  x1 = _mm_xor_si128(z1, sz1);
  x0 = _mm_sub_epi16(x0, sz0);
  x1 = _mm_sub_epi16(x1, sz1);

  /* x += round */
  x0 = _mm_add_epi16(x0, round0);
  x1 = _mm_add_epi16(x1, round1);

  /* y = (x * quant) >> 16 */
  y0 = _mm_mulhi_epi16(x0, quant_fast0);
  y1 = _mm_mulhi_epi16(x1, quant_fast1);

  /* x = abs(y) = (y ^ sz) - sz */
  y0 = _mm_xor_si128(y0, sz0);
  y1 = _mm_xor_si128(y1, sz1);
  x0 = _mm_sub_epi16(y0, sz0);
  x1 = _mm_sub_epi16(y1, sz1);

  /* qcoeff = x */
  _mm_store_si128((__m128i *)(d->qcoeff), x0);
  _mm_store_si128((__m128i *)(d->qcoeff + 8), x1);

  /* x * dequant */
  xdq0 = _mm_mullo_epi16(x0, dequant0);
  xdq1 = _mm_mullo_epi16(x1, dequant1);

  /* dqcoeff = x * dequant */
  _mm_store_si128((__m128i *)(d->dqcoeff), xdq0);
  _mm_store_si128((__m128i *)(d->dqcoeff + 8), xdq1);

  /* build a mask for the zig zag */
  zeros = _mm_setzero_si128();

  x0 = _mm_cmpeq_epi16(x0, zeros);
  x1 = _mm_cmpeq_epi16(x1, zeros);

  ones = _mm_cmpeq_epi16(zeros, zeros);

  x0 = _mm_xor_si128(x0, ones);
  x1 = _mm_xor_si128(x1, ones);

  x0 = _mm_and_si128(x0, inv_zig_zag0);
  x1 = _mm_and_si128(x1, inv_zig_zag1);

  x0 = _mm_max_epi16(x0, x1);

  /* now down to 8 */
  x1 = _mm_shuffle_epi32(x0, 0xE); // 0b00001110

  x0 = _mm_max_epi16(x0, x1);

  /* only 4 left */
  x1 = _mm_shufflelo_epi16(x0, 0xE); // 0b00001110

  x0 = _mm_max_epi16(x0, x1);

  /* okay, just 2! */
  x1 = _mm_shufflelo_epi16(x0, 0x1); // 0b00000001

  x0 = _mm_max_epi16(x0, x1);

  *d->eob = 0xFF & _mm_cvtsi128_si32(x0);
}
