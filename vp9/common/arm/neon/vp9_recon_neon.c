/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */


#include "vpx_ports/config.h"
#include "vp9/common/recon.h"
#include "vp9/common/vp9_blockd.h"

extern void vp8_recon16x16mb_neon(unsigned char *pred_ptr, short *diff_ptr, unsigned char *dst_ptr, int ystride, unsigned char *udst_ptr, unsigned char *vdst_ptr);

void vp8_recon_mb_neon(const vp8_recon_rtcd_vtable_t *rtcd, MACROBLOCKD *xd) {
  unsigned char *pred_ptr = &xd->predictor[0];
  short *diff_ptr = &xd->diff[0];
  unsigned char *dst_ptr = xd->dst.y_buffer;
  unsigned char *udst_ptr = xd->dst.u_buffer;
  unsigned char *vdst_ptr = xd->dst.v_buffer;
  int ystride = xd->dst.y_stride;
  /*int uv_stride = xd->dst.uv_stride;*/

  vp8_recon16x16mb_neon(pred_ptr, diff_ptr, dst_ptr, ystride,
                        udst_ptr, vdst_ptr);
}
