/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */


#ifndef QUANTIZE_ARM_H
#define QUANTIZE_ARM_H

#if HAVE_ARMV6

extern prototype_quantize_block(vp8_fast_quantize_b_armv6);

#undef  vp8_quantize_fastquantb
#define vp8_quantize_fastquantb vp8_fast_quantize_b_armv6

#endif /* HAVE_ARMV6 */


#if HAVE_ARMV7
extern prototype_quantize_block(vp8_fast_quantize_b_neon);

/* The neon quantizer has not been updated to match the new exact
 * quantizer introduced in commit e04e2935
 */
//#undef  vp8_quantize_fastquantb
//#define vp8_quantize_fastquantb vp8_fast_quantize_b_neon

#endif

#endif
