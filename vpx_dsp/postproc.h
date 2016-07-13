/*
 *  Copyright (c) 2016 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#ifndef VPX_DSP_POSTPROC_H_
#define VPX_DSP_POSTPROC_H_

#ifdef __cplusplus
extern "C" {
#endif

int vpx_setup_noise(int size, double sigma, char *noise);

#ifdef __cplusplus
}
#endif

#endif  // VPX_DSP_POSTPROC_H_
