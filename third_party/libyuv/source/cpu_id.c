/*
 *  Copyright (c) 2011 The LibYuv project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include "third_party/libyuv/include/libyuv/cpu_id.h"
#include "third_party/libyuv/include/libyuv/basic_types.h"  // for CPU_X86

#ifdef _MSC_VER
#include <intrin.h>
#endif

// TODO(fbarchard): Use cpuid.h when gcc 4.4 is used on OSX and Linux.
#if (defined(__pic__) || defined(__APPLE__)) && defined(__i386__)
static inline void __cpuid(int cpu_info[4], int info_type) {
  __asm__ volatile (
    "mov %%ebx, %%edi\n"
    "cpuid\n"
    "xchg %%edi, %%ebx\n"
    : "=a"(cpu_info[0]), "=D"(cpu_info[1]), "=c"(cpu_info[2]), "=d"(cpu_info[3])
    : "a"(info_type)
  );
}
#elif defined(__i386__) || defined(__x86_64__)
static inline void __cpuid(int cpu_info[4], int info_type) {
  __asm__ volatile (
    "cpuid\n"
    : "=a"(cpu_info[0]), "=b"(cpu_info[1]), "=c"(cpu_info[2]), "=d"(cpu_info[3])
    : "a"(info_type)
  );
}
#endif

//namespace libyuv {

// CPU detect function for SIMD instruction sets.
static int cpu_info_initialized_ = 0;
static int cpu_info_ = 0;

// Global lock for cpu initialization.
static void InitCpuFlags() {
#ifdef CPU_X86
  int cpu_info[4];
  __cpuid(cpu_info, 1);
  cpu_info_ = (cpu_info[2] & 0x00000200 ? kCpuHasSSSE3 : 0) |
              (cpu_info[3] & 0x04000000 ? kCpuHasSSE2 : 0);
#elif defined(__ARM_NEON__)
  // gcc -mfpu=neon defines __ARM_NEON__
  // if code is specifically built for Neon-only, enable the flag.
  cpu_info_ |= kCpuHasNEON;
#else
  cpu_info_ = 0;
#endif
  cpu_info_initialized_ = 1;
}

void MaskCpuFlagsForTest(int enable_flags) {
  InitCpuFlags();
  cpu_info_ &= enable_flags;
}

int TestCpuFlag(int flag) {
  if (!cpu_info_initialized_) {
    InitCpuFlags();
  }
  return cpu_info_ & flag ? 1 : 0;
}

//}  // namespace libyuv
