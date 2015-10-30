/*
 *  Copyright (c) 2011 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#ifndef VP10_COMMON_FILTER_H_
#define VP10_COMMON_FILTER_H_

#include "./vpx_config.h"
#include "vpx/vpx_integer.h"
#include "vpx_dsp/vpx_filter.h"
#include "vpx_ports/mem.h"


#ifdef __cplusplus
extern "C" {
#endif

#define EIGHTTAP            0
#define EIGHTTAP_SMOOTH     1
#define EIGHTTAP_SHARP      2

#if CONFIG_EXT_INTERP
#define SUPPORT_NONINTERPOLATING_FILTERS 0  /* turn it on for experimentation */
#define EIGHTTAP_SMOOTH2    3
#define SWITCHABLE_FILTERS  4 /* Number of switchable filters */
#else
#define SWITCHABLE_FILTERS  3 /* Number of switchable filters */
#endif  // CONFIG_EXT_INTERP
// The codec can operate in four possible inter prediction filter mode:
// 8-tap, 8-tap-smooth, 8-tap-sharp, and switching between the three.

#define BILINEAR            (SWITCHABLE_FILTERS)
#define SWITCHABLE          (SWITCHABLE_FILTERS + 1)  /* the last one */
#define SWITCHABLE_FILTER_CONTEXTS (SWITCHABLE_FILTERS + 1)

typedef uint8_t INTERP_FILTER;

extern const InterpKernel *vp10_filter_kernels[SWITCHABLE_FILTERS + 1];

#ifdef __cplusplus
}  // extern "C"
#endif

#endif  // VP10_COMMON_FILTER_H_
