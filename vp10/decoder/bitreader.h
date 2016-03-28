/*
 *  Copyright (c) 2016 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

/* The purpose of this header is to provide compile time pluggable bit reader
 * implementations with a common interface. */

#ifndef VPX10_DECODER_BITREADER_H_
#define VPX10_DECODER_BITREADER_H_

#include "./vpx_config.h"

#if CONFIG_ANS
#include "vp10/common/ans.h"
#include "vpx/vp8dx.h"  // for vp10_decrypt_cb
#define vp10_reader struct AnsDecoder
#define vp10_reader_has_error ans_reader_has_error
#define vp10_read uabs_read
#define vp10_read_bit uabs_read_bit
#define vp10_read_literal uabs_read_literal
#define vp10_read_tree uabs_read_tree
#else
#include "vpx_dsp/bitreader.h"
#define vp10_reader vpx_reader
#define vp10_reader_has_error vpx_reader_has_error
#define vp10_read vpx_read
#define vp10_read_bit vpx_read_bit
#define vp10_read_literal vpx_read_literal
#define vp10_read_tree vpx_read_tree
#endif

#endif  // VPX10_DECODER_BITREADER_H_
