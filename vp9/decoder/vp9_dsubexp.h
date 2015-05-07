/*
 *  Copyright (c) 2013 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */


#ifndef VP9_DECODER_VP9_DSUBEXP_H_
#define VP9_DECODER_VP9_DSUBEXP_H_

#include "vp9/decoder/vp9_reader.h"

#ifdef __cplusplus
extern "C" {
#endif

void vp9_diff_update_prob(vp9_reader *r, vp9_prob* p);

// num_values is the number of values the symbol can take
int vp9_read_primitive_uniform(vp9_reader *r, unsigned int num_values);

// k is the parameter of the subexponential code
int vp9_read_primitive_subexp(vp9_reader *r, unsigned int k);

// mag_bits is number of bits for magnitude. The alphabet is of size
// 2 * 2^mag_bits + 1, symmetric around 0, where one bit is used to
// indicate 0 or non-zero, mag_bits bits are used to indicate magnitide
// and 1 more bit for the sign if non-zero.
int vp9_read_primitive_symmetric(vp9_reader *r, unsigned int mag_bits);

#ifdef __cplusplus
}  // extern "C"
#endif

#endif  // VP9_DECODER_VP9_DSUBEXP_H_
