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

#ifndef VPX10_ENCODER_BITWRITER_H_
#define VPX10_ENCODER_BITWRITER_H_

#include "./vpx_config.h"
#include "vpx_dsp/prob.h"

#if CONFIG_ANS
typedef struct BufAnsCoder BufAnsCoder;
#include "vp10/encoder/buf_ans.h"
#define vp10_writer BufAnsCoder
#define vp10_write buf_uabs_write
#define vp10_write_bit buf_uabs_write_bit
#define vp10_write_literal buf_uabs_write_literal
#else
#include "vpx_dsp/bitwriter.h"
#define vp10_writer vpx_writer
#define vp10_write vpx_write
#define vp10_write_bit vpx_write_bit
#define vp10_write_literal vpx_write_literal
#endif

#endif  // VPX10_ENCODER_BITWRITER_H_
