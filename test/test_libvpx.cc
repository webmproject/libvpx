/*
 *  Copyright (c) 2012 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */
#include <string>

#include "third_party/googletest/src/include/gtest/gtest.h"

#include "./vpx_config.h"
#if VPX_ARCH_ARM
#include "vpx_ports/arm.h"
#endif
#if VPX_ARCH_X86 || VPX_ARCH_X86_64
#include "vpx_ports/x86.h"
#endif
extern "C" {
#if CONFIG_VP8
extern void vp8_rtcd();
#endif  // CONFIG_VP8
#if CONFIG_VP9
extern void vp9_rtcd();
#endif  // CONFIG_VP9
extern void vpx_dsp_rtcd();
extern void vpx_scale_rtcd();
}

#if (!CONFIG_SHARED && VPX_ARCH_ARM) || VPX_ARCH_X86 || VPX_ARCH_X86_64
static void append_negative_gtest_filter(const char *str) {
  std::string filter = ::testing::FLAGS_gtest_filter;
  // Negative patterns begin with one '-' followed by a ':' separated list.
  if (filter.find('-') == std::string::npos) filter += '-';
  filter += str;
  ::testing::FLAGS_gtest_filter = filter;
}
#endif  // (!CONFIG_SHARED && VPX_ARCH_ARM) || VPX_ARCH_X86 || VPX_ARCH_X86_64

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);

#if !CONFIG_SHARED
#if VPX_ARCH_AARCH64
  const int caps = arm_cpu_caps();
  if (!(caps & HAS_NEON_DOTPROD)) {
    append_negative_gtest_filter(":NEON_DOTPROD.*:NEON_DOTPROD/*");
  }
  if (!(caps & HAS_NEON_I8MM)) {
    append_negative_gtest_filter(":NEON_I8MM.*:NEON_I8MM/*");
  }
#elif VPX_ARCH_ARM
  const int caps = arm_cpu_caps();
  if (!(caps & HAS_NEON)) {
    append_negative_gtest_filter(":NEON.*:NEON/*");
  }
#endif  // VPX_ARCH_ARM
#endif  // !CONFIG_SHARED

#if VPX_ARCH_X86 || VPX_ARCH_X86_64
  const int simd_caps = x86_simd_caps();
  if (!(simd_caps & HAS_MMX)) append_negative_gtest_filter(":MMX.*:MMX/*");
  if (!(simd_caps & HAS_SSE)) append_negative_gtest_filter(":SSE.*:SSE/*");
  if (!(simd_caps & HAS_SSE2)) append_negative_gtest_filter(":SSE2.*:SSE2/*");
  if (!(simd_caps & HAS_SSE3)) append_negative_gtest_filter(":SSE3.*:SSE3/*");
  if (!(simd_caps & HAS_SSSE3)) {
    append_negative_gtest_filter(":SSSE3.*:SSSE3/*");
  }
  if (!(simd_caps & HAS_SSE4_1)) {
    append_negative_gtest_filter(":SSE4_1.*:SSE4_1/*");
  }
  if (!(simd_caps & HAS_AVX)) append_negative_gtest_filter(":AVX.*:AVX/*");
  if (!(simd_caps & HAS_AVX2)) append_negative_gtest_filter(":AVX2.*:AVX2/*");
  if (!(simd_caps & HAS_AVX512)) {
    append_negative_gtest_filter(":AVX512.*:AVX512/*");
  }
#endif  // VPX_ARCH_X86 || VPX_ARCH_X86_64

#if !CONFIG_SHARED
// Shared library builds don't support whitebox tests
// that exercise internal symbols.
#if CONFIG_VP8
  vp8_rtcd();
#endif  // CONFIG_VP8
#if CONFIG_VP9
  vp9_rtcd();
#endif  // CONFIG_VP9
  vpx_dsp_rtcd();
  vpx_scale_rtcd();
#endif  // !CONFIG_SHARED

  return RUN_ALL_TESTS();
}
