/*
 *  Copyright (c) 2017 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#ifndef VPX_DSP_PPC_TYPES_VSX_H_
#define VPX_DSP_PPC_TYPES_VSX_H_

#include <altivec.h>

typedef vector signed char int8x16_t;
typedef vector unsigned char uint8x16_t;
typedef vector signed short int16x8_t;
typedef vector unsigned short uint16x8_t;
typedef vector signed int int32x4_t;
typedef vector unsigned int uint32x4_t;

#ifdef WORDS_BIGENDIAN
#define unpack_to_u16_h(v) \
  (uint16x8_t) vec_mergeh(vec_splat_u8(0), (uint8x16_t)v)
#define unpack_to_u16_l(v) \
  (uint16x8_t) vec_mergel(vec_splat_u8(0), (uint8x16_t)v)
#define unpack_to_s16_h(v) \
  (int16x8_t) vec_mergeh(vec_splat_u8(0), (uint8x16_t)v)
#define unpack_to_s16_l(v) \
  (int16x8_t) vec_mergel(vec_splat_u8(0), (uint8x16_t)v)
#else
#define unpack_to_u16_h(v) \
  (uint16x8_t) vec_mergeh((uint8x16_t)v, vec_splat_u8(0))
#define unpack_to_u16_l(v) \
  (uint16x8_t) vec_mergel((uint8x16_t)v, vec_splat_u8(0))
#define unpack_to_s16_h(v) \
  (int16x8_t) vec_mergeh((uint8x16_t)v, vec_splat_u8(0))
#define unpack_to_s16_l(v) \
  (int16x8_t) vec_mergel((uint8x16_t)v, vec_splat_u8(0))
#endif

#endif  // VPX_DSP_PPC_TYPES_VSX_H_
