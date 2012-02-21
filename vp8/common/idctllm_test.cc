/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */


 extern "C" {
    void vp8_short_idct4x4llm_c(short *input, unsigned char *pred_ptr,
                            int pred_stride, unsigned char *dst_ptr,
                            int dst_stride);
}

#include "vpx_config.h"
#include "idctllm_test.h"
namespace
{

INSTANTIATE_TEST_CASE_P(C, IDCTTest,
                        ::testing::Values(vp8_short_idct4x4llm_c));

} // namespace

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
