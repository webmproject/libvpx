/*
 *  Copyright (c) 2022 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include "./vp8_rtcd.h"
#include "vp8/common/blockd.h"
#include "vpx_util/loongson_intrinsics.h"

static void idct4x4_addconst_lsx(int16_t in_dc, uint8_t *pred,
                                 int32_t pred_stride, uint8_t *dest,
                                 int32_t dest_stride) {
  __m128i vec, res0, res1, res2, res3, dst0, dst1;
  __m128i pred0, pred1, pred2, pred3;
  __m128i zero = __lsx_vldi(0);

  int32_t pred_stride2 = pred_stride << 1;
  int32_t pred_stride3 = pred_stride2 + pred_stride;

  vec = __lsx_vreplgr2vr_h(in_dc);
  vec = __lsx_vsrari_h(vec, 3);
  pred0 = __lsx_vld(pred, 0);
  DUP2_ARG2(__lsx_vldx, pred, pred_stride, pred, pred_stride2, pred1, pred2);
  pred3 = __lsx_vldx(pred, pred_stride3);
  DUP4_ARG2(__lsx_vilvl_b, zero, pred0, zero, pred1, zero, pred2, zero, pred3,
            res0, res1, res2, res3);
  DUP4_ARG2(__lsx_vadd_h, res0, vec, res1, vec, res2, vec, res3, vec, res0,
            res1, res2, res3);
  res0 = __lsx_vclip255_h(res0);
  res1 = __lsx_vclip255_h(res1);
  res2 = __lsx_vclip255_h(res2);
  res3 = __lsx_vclip255_h(res3);

  DUP2_ARG2(__lsx_vpickev_b, res1, res0, res3, res2, dst0, dst1);
  dst0 = __lsx_vpickev_w(dst1, dst0);
  __lsx_vstelm_w(dst0, dest, 0, 0);
  dest += dest_stride;
  __lsx_vstelm_w(dst0, dest, 0, 1);
  dest += dest_stride;
  __lsx_vstelm_w(dst0, dest, 0, 2);
  dest += dest_stride;
  __lsx_vstelm_w(dst0, dest, 0, 3);
}

void vp8_dc_only_idct_add_lsx(int16_t input_dc, uint8_t *pred_ptr,
                              int32_t pred_stride, uint8_t *dst_ptr,
                              int32_t dst_stride) {
  idct4x4_addconst_lsx(input_dc, pred_ptr, pred_stride, dst_ptr, dst_stride);
}
