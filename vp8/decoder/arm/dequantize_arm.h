/*
 *  Copyright (c) 2010 The VP8 project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */


#ifndef DEQUANTIZE_ARM_H
#define DEQUANTIZE_ARM_H

#if HAVE_ARMV6
extern prototype_dequant_block(vp8_dequantize_b_v6);
extern prototype_dequant_idct(vp8_dequant_idct_v6);
extern prototype_dequant_idct_dc(vp8_dequant_dc_idct_v6);

#undef  vp8_dequant_block
#define vp8_dequant_block vp8_dequantize_b_v6

#undef  vp8_dequant_idct
#define vp8_dequant_idct vp8_dequant_idct_v6

#undef  vp8_dequant_idct_dc
#define vp8_dequant_idct_dc vp8_dequant_dc_idct_v6
#endif

#if HAVE_ARMV7
extern prototype_dequant_block(vp8_dequantize_b_neon);
extern prototype_dequant_idct(vp8_dequant_idct_neon);
extern prototype_dequant_idct_dc(vp8_dequant_dc_idct_neon);

#undef  vp8_dequant_block
#define vp8_dequant_block vp8_dequantize_b_neon

#undef  vp8_dequant_idct
#define vp8_dequant_idct vp8_dequant_idct_neon

#undef  vp8_dequant_idct_dc
#define vp8_dequant_idct_dc vp8_dequant_dc_idct_neon
#endif

#endif
