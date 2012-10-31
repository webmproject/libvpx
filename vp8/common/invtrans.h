/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */


#ifndef __INC_INVTRANS_H
#define __INC_INVTRANS_H

#include "vpx_ports/config.h"
#include "idct.h"
#include "blockd.h"

extern void vp9_inverse_transform_b_4x4(const vp9_idct_rtcd_vtable_t *rtcd, BLOCKD *b, int pitch);
extern void vp9_inverse_transform_mb_4x4(const vp9_idct_rtcd_vtable_t *rtcd, MACROBLOCKD *xd);
extern void vp9_inverse_transform_mby_4x4(const vp9_idct_rtcd_vtable_t *rtcd, MACROBLOCKD *xd);
extern void vp9_inverse_transform_mbuv_4x4(const vp9_idct_rtcd_vtable_t *rtcd, MACROBLOCKD *xd);

extern void vp9_inverse_transform_b_8x8(const vp9_idct_rtcd_vtable_t *rtcd, short *input_dqcoeff, short *output_coeff, int pitch);
extern void vp9_inverse_transform_mb_8x8(const vp9_idct_rtcd_vtable_t *rtcd, MACROBLOCKD *xd);
extern void vp9_inverse_transform_mby_8x8(const vp9_idct_rtcd_vtable_t *rtcd, MACROBLOCKD *xd);
extern void vp9_inverse_transform_mbuv_8x8(const vp9_idct_rtcd_vtable_t *rtcd, MACROBLOCKD *xd);

extern void vp9_inverse_transform_b_16x16(const vp9_idct_rtcd_vtable_t *rtcd,
                                          short *input_dqcoeff, short *output_coeff,
                                          int pitch);
extern void vp9_inverse_transform_mb_16x16(const vp9_idct_rtcd_vtable_t *rtcd, MACROBLOCKD *xd);
extern void vp9_inverse_transform_mby_16x16(const vp9_idct_rtcd_vtable_t *rtcd, MACROBLOCKD *xd);
#endif
