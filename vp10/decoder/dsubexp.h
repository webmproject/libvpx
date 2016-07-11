/*
 *  Copyright (c) 2013 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */


#ifndef VP10_DECODER_DSUBEXP_H_
#define VP10_DECODER_DSUBEXP_H_

#include "vp10/decoder/bitreader.h"

#ifdef __cplusplus
extern "C" {
#endif

void vp10_diff_update_prob(vp10_reader *r, vpx_prob* p);

#ifdef __cplusplus
}  // extern "C"
#endif

// mag_bits is number of bits for magnitude. The alphabet is of size
// 2 * 2^mag_bits + 1, symmetric around 0, where one bit is used to
// indicate 0 or non-zero, mag_bits bits are used to indicate magnitide
// and 1 more bit for the sign if non-zero.
int vp10_read_primitive_symmetric(vp10_reader *r, unsigned int mag_bits);
#endif  // VP10_DECODER_DSUBEXP_H_
