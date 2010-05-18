/*
 *  Copyright (c) 2010 The VP8 project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license and patent
 *  grant that can be found in the LICENSE file in the root of the source
 *  tree. All contributing project authors may be found in the AUTHORS
 *  file in the root of the source tree.
 */


#ifndef IDCT_ARM_H
#define IDCT_ARM_H

#if HAVE_ARMV6
extern prototype_idct(vp8_short_idct4x4llm_1_v6);
extern prototype_idct(vp8_short_idct4x4llm_v6_dual);
extern prototype_idct_scalar(vp8_dc_only_idct_armv6);
extern prototype_second_order(vp8_short_inv_walsh4x4_1_armv6);
extern prototype_second_order(vp8_short_inv_walsh4x4_armv6);

#undef  vp8_idct_idct1
#define vp8_idct_idct1 vp8_short_idct4x4llm_1_v6

#undef  vp8_idct_idct16
#define vp8_idct_idct16 vp8_short_idct4x4llm_v6_dual

#undef  vp8_idct_idct1_scalar
#define vp8_idct_idct1_scalar vp8_dc_only_idct_armv6

#undef  vp8_idct_iwalsh1
#define vp8_idct_iwalsh1 vp8_short_inv_walsh4x4_1_armv6

#undef  vp8_idct_iwalsh16
#define vp8_idct_iwalsh16 vp8_short_inv_walsh4x4_armv6
#endif

#if HAVE_ARMV7
extern prototype_idct(vp8_short_idct4x4llm_1_neon);
extern prototype_idct(vp8_short_idct4x4llm_neon);
extern prototype_idct_scalar(vp8_dc_only_idct_neon);
extern prototype_second_order(vp8_short_inv_walsh4x4_1_neon);
extern prototype_second_order(vp8_short_inv_walsh4x4_neon);

#undef  vp8_idct_idct1
#define vp8_idct_idct1 vp8_short_idct4x4llm_1_neon

#undef  vp8_idct_idct16
#define vp8_idct_idct16 vp8_short_idct4x4llm_neon

#undef  vp8_idct_idct1_scalar
#define vp8_idct_idct1_scalar vp8_dc_only_idct_neon

#undef  vp8_idct_iwalsh1
#define vp8_idct_iwalsh1 vp8_short_inv_walsh4x4_1_neon

#undef  vp8_idct_iwalsh16
#define vp8_idct_iwalsh16 vp8_short_inv_walsh4x4_neon
#endif

#endif
