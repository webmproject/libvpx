/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */


#ifndef VARIANCE_ARM_H
#define VARIANCE_ARM_H

#if HAVE_ARMV6

extern prototype_sad(vp9_sad16x16_armv6);
extern prototype_variance(vp9_variance16x16_armv6);
extern prototype_variance(vp9_variance8x8_armv6);
extern prototype_subpixvariance(vp9_sub_pixel_variance16x16_armv6);
extern prototype_subpixvariance(vp9_sub_pixel_variance8x8_armv6);
extern prototype_variance(vp9_variance_halfpixvar16x16_h_armv6);
extern prototype_variance(vp9_variance_halfpixvar16x16_v_armv6);
extern prototype_variance(vp9_variance_halfpixvar16x16_hv_armv6);
extern prototype_variance(vp9_mse16x16_armv6);

#if !CONFIG_RUNTIME_CPU_DETECT

#undef  vp9_variance_sad16x16
#define vp9_variance_sad16x16 vp9_sad16x16_armv6

#undef  vp9_variance_subpixvar16x16
#define vp9_variance_subpixvar16x16 vp9_sub_pixel_variance16x16_armv6

#undef  vp9_variance_subpixvar8x8
#define vp9_variance_subpixvar8x8 vp9_sub_pixel_variance8x8_armv6

#undef  vp9_variance_var16x16
#define vp9_variance_var16x16 vp9_variance16x16_armv6

#undef  vp9_variance_mse16x16
#define vp9_variance_mse16x16 vp9_mse16x16_armv6

#undef  vp9_variance_var8x8
#define vp9_variance_var8x8 vp9_variance8x8_armv6

#undef  vp9_variance_halfpixvar16x16_h
#define vp9_variance_halfpixvar16x16_h vp9_variance_halfpixvar16x16_h_armv6

#undef  vp9_variance_halfpixvar16x16_v
#define vp9_variance_halfpixvar16x16_v vp9_variance_halfpixvar16x16_v_armv6

#undef  vp9_variance_halfpixvar16x16_hv
#define vp9_variance_halfpixvar16x16_hv vp9_variance_halfpixvar16x16_hv_armv6

#endif /* !CONFIG_RUNTIME_CPU_DETECT */

#endif /* HAVE_ARMV6 */


#if HAVE_ARMV7
extern prototype_sad(vp9_sad4x4_neon);
extern prototype_sad(vp9_sad8x8_neon);
extern prototype_sad(vp9_sad8x16_neon);
extern prototype_sad(vp9_sad16x8_neon);
extern prototype_sad(vp9_sad16x16_neon);

extern prototype_variance(vp9_variance8x8_neon);
extern prototype_variance(vp9_variance8x16_neon);
extern prototype_variance(vp9_variance16x8_neon);
extern prototype_variance(vp9_variance16x16_neon);

extern prototype_subpixvariance(vp9_sub_pixel_variance8x8_neon);
extern prototype_subpixvariance(vp9_sub_pixel_variance16x16_neon);
extern prototype_subpixvariance(vp9_sub_pixel_variance16x16_neon_func);
extern prototype_variance(vp9_variance_halfpixvar16x16_h_neon);
extern prototype_variance(vp9_variance_halfpixvar16x16_v_neon);
extern prototype_variance(vp9_variance_halfpixvar16x16_hv_neon);

extern prototype_variance(vp9_mse16x16_neon);

#if !CONFIG_RUNTIME_CPU_DETECT
#undef  vp9_variance_sad4x4
#define vp9_variance_sad4x4 vp9_sad4x4_neon

#undef  vp9_variance_sad8x8
#define vp9_variance_sad8x8 vp9_sad8x8_neon

#undef  vp9_variance_sad8x16
#define vp9_variance_sad8x16 vp9_sad8x16_neon

#undef  vp9_variance_sad16x8
#define vp9_variance_sad16x8 vp9_sad16x8_neon

#undef  vp9_variance_sad16x16
#define vp9_variance_sad16x16 vp9_sad16x16_neon

#undef  vp9_variance_var8x8
#define vp9_variance_var8x8 vp9_variance8x8_neon

#undef  vp9_variance_var8x16
#define vp9_variance_var8x16 vp9_variance8x16_neon

#undef  vp9_variance_var16x8
#define vp9_variance_var16x8 vp9_variance16x8_neon

#undef  vp9_variance_var16x16
#define vp9_variance_var16x16 vp9_variance16x16_neon

#undef  vp9_variance_subpixvar8x8
#define vp9_variance_subpixvar8x8 vp9_sub_pixel_variance8x8_neon

#undef  vp9_variance_subpixvar16x16
#define vp9_variance_subpixvar16x16 vp9_sub_pixel_variance16x16_neon

#undef  vp9_variance_halfpixvar16x16_h
#define vp9_variance_halfpixvar16x16_h vp9_variance_halfpixvar16x16_h_neon

#undef  vp9_variance_halfpixvar16x16_v
#define vp9_variance_halfpixvar16x16_v vp9_variance_halfpixvar16x16_v_neon

#undef  vp9_variance_halfpixvar16x16_hv
#define vp9_variance_halfpixvar16x16_hv vp9_variance_halfpixvar16x16_hv_neon

#undef  vp9_variance_mse16x16
#define vp9_variance_mse16x16 vp9_mse16x16_neon

#endif

#endif

#endif
