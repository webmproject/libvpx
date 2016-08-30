/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#ifndef AV1_ENCODER_BITSTREAM_H_
#define AV1_ENCODER_BITSTREAM_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "av1/encoder/encoder.h"

void av1_pack_bitstream(AV1_COMP *const cpi, uint8_t *dest, size_t *size);

void av1_encode_token_init(void);

static INLINE int av1_preserve_existing_gf(AV1_COMP *cpi) {
#if CONFIG_EXT_REFS
  // Do not swap gf and arf indices for internal overlay frames
  return !cpi->multi_arf_allowed && cpi->rc.is_src_frame_alt_ref &&
         !cpi->rc.is_src_frame_ext_arf;
#else
  return !cpi->multi_arf_allowed && cpi->refresh_golden_frame &&
         cpi->rc.is_src_frame_alt_ref;
#endif  // CONFIG_EXT_REFS
}

#ifdef __cplusplus
}  // extern "C"
#endif

#endif  // AV1_ENCODER_BITSTREAM_H_
