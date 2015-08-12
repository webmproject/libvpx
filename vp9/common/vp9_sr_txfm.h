/*
 *  Copyright (c) 2015 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#ifndef VP9_COMMON_VP9_PALETTE_H_
#define VP9_COMMON_VP9_PALETTE_H_

#include "vp9/common/vp9_blockd.h"
#include "vp9/common/vp9_entropymode.h"

#if CONFIG_SR_MODE
#define MIN_SR_TX_SIZE 2

#if CONFIG_TX_64X64
#define MAX_SR_TX_SIZE 4
#else
#define MAX_SR_TX_SIZE 3
#endif

#define UPSCALE_FILTER_TAPS 10
#define UPSCALE_FILTER_SHIFT 7
void sr_lowpass(int16_t *src, int src_stride, int16_t *dst, int dst_stride,
                int w, int h);
#if SR_USE_MULTI_F
void sr_recon(int16_t *src, int src_stride, uint8_t *dst, int dst_stride,
              int w, int h, int f_hor, int f_ver);
static INLINE int idx_to_h(int idx) {
  return (idx % SR_USFILTER_NUM_D);
}
static INLINE int idx_to_v(int idx) {
  return (idx / SR_USFILTER_NUM_D);
}
static INLINE int hv_to_idx(int hor, int ver) {
  return hor + ver * SR_USFILTER_NUM_D;
}
#else
void sr_recon(int16_t *src, int src_stride, uint8_t *dst, int dst_stride,
              int w, int h);
#endif
int is_enable_srmode(BLOCK_SIZE bsize);
#endif  // CONFIG_SR_MODE

#endif  // VP9_COMMON_VP9_PALETTE_H_
