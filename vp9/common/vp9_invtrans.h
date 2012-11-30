/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#ifndef VP9_COMMON_VP9_INVTRANS_H_
#define VP9_COMMON_VP9_INVTRANS_H_

#include "vpx_ports/config.h"
#include "vp9/common/vp9_blockd.h"

extern void vp9_inverse_transform_b_4x4(MACROBLOCKD *xd, int block, int pitch);

extern void vp9_inverse_transform_mb_4x4(MACROBLOCKD *xd);

extern void vp9_inverse_transform_mby_4x4(MACROBLOCKD *xd);

extern void vp9_inverse_transform_mbuv_4x4(MACROBLOCKD *xd);

extern void vp9_inverse_transform_b_8x8(short *input_dqcoeff,
                                        short *output_coeff, int pitch);

extern void vp9_inverse_transform_mb_8x8(MACROBLOCKD *xd);

extern void vp9_inverse_transform_mby_8x8(MACROBLOCKD *xd);

extern void vp9_inverse_transform_mbuv_8x8(MACROBLOCKD *xd);

extern void vp9_inverse_transform_b_16x16(short *input_dqcoeff,
                                          short *output_coeff, int pitch);

extern void vp9_inverse_transform_mb_16x16(MACROBLOCKD *xd);

extern void vp9_inverse_transform_mby_16x16(MACROBLOCKD *xd);

#endif  // __INC_INVTRANS_H
