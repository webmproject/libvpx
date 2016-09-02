/*
 * Copyright (c) 2016, Alliance for Open Media. All rights reserved
 *
 * This source code is subject to the terms of the BSD 2 Clause License and
 * the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
 * was not distributed with this source code in the LICENSE file, you can
 * obtain it at www.aomedia.org/license/software. If the Alliance for Open
 * Media Patent License 1.0 was not distributed with this source code in the
 * PATENTS file, you can obtain it at www.aomedia.org/license/patent.
 */

#ifndef AV1_DECODER_DSUBEXP_H_
#define AV1_DECODER_DSUBEXP_H_

#include "av1/decoder/bitreader.h"

#ifdef __cplusplus
extern "C" {
#endif

void av1_diff_update_prob(aom_reader *r, aom_prob *p);

#ifdef __cplusplus
}  // extern "C"
#endif

// mag_bits is number of bits for magnitude. The alphabet is of size
// 2 * 2^mag_bits + 1, symmetric around 0, where one bit is used to
// indicate 0 or non-zero, mag_bits bits are used to indicate magnitide
// and 1 more bit for the sign if non-zero.
int aom_read_primitive_symmetric(aom_reader *r, unsigned int mag_bits);
#endif  // AV1_DECODER_DSUBEXP_H_
