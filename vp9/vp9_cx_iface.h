/*
 *  Copyright (c) 2019 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#ifndef VPX_VP9_VP9_CX_IFACE_H_
#define VPX_VP9_VP9_CX_IFACE_H_
#include "vp9/encoder/vp9_encoder.h"
#include "vp9/common/vp9_onyxc_int.h"

#ifdef __cplusplus
extern "C" {
#endif

VP9EncoderConfig vp9_get_encoder_config(int frame_width, int frame_height,
                                        int target_bitrate,
                                        vpx_enc_pass enc_pass);
FRAME_INFO vp9_get_frame_info(const VP9EncoderConfig *oxcf);

#ifdef __cplusplus
}  // extern "C"
#endif

#endif  // VPX_VP9_VP9_CX_IFACE_H_
