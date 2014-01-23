/*
 *  Copyright (c) 2013 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#ifndef TEST_TEST_VECTORS_H_
#define TEST_TEST_VECTORS_H_

#include "./vpx_config.h"

namespace libvpx_test {

#if CONFIG_VP8_DECODER
const int kNumVp8TestVectors = 62;
extern const char *kVP8TestVectors[kNumVp8TestVectors];
#endif

#if CONFIG_VP9_DECODER
#if CONFIG_NON420
const int kNumVp9TestVectors = 217;
#else
const int kNumVp9TestVectors = 216;
#endif

extern const char *kVP9TestVectors[kNumVp9TestVectors];
#endif  // CONFIG_VP9_DECODER

}  // namespace libvpx_test

#endif  // TEST_TEST_VECTORS_H_
