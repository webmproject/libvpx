/*
 *  Copyright (c) 2013 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#ifndef VP9_COMMON_VP9_SCAN_H_
#define VP9_COMMON_VP9_SCAN_H_

#include "vpx/vpx_integer.h"
#include "vpx_ports/mem.h"

#include "vp9/common/vp9_enums.h"
#include "vp9/common/vp9_blockd.h"

#ifdef __cplusplus
extern "C" {
#endif

#define MAX_NEIGHBORS 2

typedef struct {
  const int16_t *scan;
  const int16_t *iscan;
  const int16_t *neighbors;
} scan_order;

extern const scan_order vp9_default_scan_orders[TX_SIZES];
extern const scan_order vp9_scan_orders[TX_SIZES][TX_TYPES];

#if CONFIG_TX_SKIP
// pixel domain default scan orders
extern const scan_order vp9_default_scan_orders_pxd[TX_SIZES];

extern int16_t vp9_default_scan_pxd_4x4[16];
extern int16_t vp9_default_scan_pxd_8x8[64];
extern int16_t vp9_default_scan_pxd_16x16[256];
extern int16_t vp9_default_scan_pxd_32x32[1024];

extern int16_t vp9_default_iscan_pxd_4x4[16];
extern int16_t vp9_default_iscan_pxd_8x8[64];
extern int16_t vp9_default_iscan_pxd_16x16[256];
extern int16_t vp9_default_iscan_pxd_32x32[1024];

extern int16_t vp9_default_scan_pxd_4x4_neighbors[17 * MAX_NEIGHBORS];
extern int16_t vp9_default_scan_pxd_8x8_neighbors[65 * MAX_NEIGHBORS];
extern int16_t vp9_default_scan_pxd_16x16_neighbors[257 * MAX_NEIGHBORS];
extern int16_t vp9_default_scan_pxd_32x32_neighbors[1025 * MAX_NEIGHBORS];

#if CONFIG_TX64X64
extern int16_t vp9_default_scan_pxd_64x64[4096];
extern int16_t vp9_default_iscan_pxd_64x64[4096];
extern int16_t vp9_default_scan_pxd_64x64_neighbors[4097 * MAX_NEIGHBORS];
#endif  // CONFIG_TX64X64
#endif  // CONFIG_TX_SKIP

static INLINE int get_coef_context(const int16_t *neighbors,
                                   const uint8_t *token_cache, int c) {
  return (1 + token_cache[neighbors[MAX_NEIGHBORS * c + 0]] +
          token_cache[neighbors[MAX_NEIGHBORS * c + 1]]) >> 1;
}

#ifdef __cplusplus
}  // extern "C"
#endif

#endif  // VP9_COMMON_VP9_SCAN_H_
