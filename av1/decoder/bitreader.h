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

#ifndef AOM10_DECODER_BITREADER_H_
#define AOM10_DECODER_BITREADER_H_

#include "./aom_config.h"

#if CONFIG_ANS
#include "av1/common/ans.h"
#include "aom/aomdx.h"  // for av1_decrypt_cb
#define aom_reader struct AnsDecoder
#define aom_reader_has_error ans_reader_has_error
#define aom_read uabs_read
#define aom_read_bit uabs_read_bit
#define aom_read_literal uabs_read_literal
#define aom_read_tree uabs_read_tree
#else
#include "aom_dsp/bitreader.h"
#define aom_reader aom_reader
#define aom_reader_has_error aom_reader_has_error
#define aom_read aom_read
#define aom_read_bit aom_read_bit
#define aom_read_literal aom_read_literal
#define aom_read_tree aom_read_tree
#endif

#endif  // AOM10_DECODER_BITREADER_H_
