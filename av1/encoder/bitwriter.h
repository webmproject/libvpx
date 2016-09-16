/*
 *  Copyright (c) 2016 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

/* The purpose of this header is to provide compile time pluggable bit writer
 * implementations with a common interface. */

#ifndef AOM10_ENCODER_BITWRITER_H_
#define AOM10_ENCODER_BITWRITER_H_

#include "./aom_config.h"
#include "aom_dsp/prob.h"
// Include bitwriter.h in the CONFIG_ANS to keep ANS building while
// porting from VP10 style entropy coder abstraction to the aom/master style
// entropy coder abstractions.
#include "aom_dsp/bitwriter.h"

#if CONFIG_ANS
typedef struct BufAnsCoder BufAnsCoder;
#include "av1/encoder/buf_ans.h"
#define aom_writer BufAnsCoder
#define aom_write buf_uabs_write
#define aom_write_bit buf_uabs_write_bit
#define aom_write_literal buf_uabs_write_literal
#else
#define aom_writer aom_writer
#define aom_write aom_write
#define aom_write_bit aom_write_bit
#define aom_write_literal aom_write_literal
#endif

#endif  // AOM10_ENCODER_BITWRITER_H_
