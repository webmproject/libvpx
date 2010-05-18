/*
 *  Copyright (c) 2010 The VP8 project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license and patent
 *  grant that can be found in the LICENSE file in the root of the source
 *  tree. All contributing project authors may be found in the AUTHORS
 *  file in the root of the source tree.
 */


#ifndef __INC_PSNR_H
#define __INC_PSNR_H

extern double vp8_mse2psnr(double Samples, double Peak, double Mse);
extern double vp8_calc_psnr(YV12_BUFFER_CONFIG *source, YV12_BUFFER_CONFIG *dest, double *YPsnr, double *UPsnr, double *VPsnr, double *sq_error);

#endif
