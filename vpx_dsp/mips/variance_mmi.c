/*
 *  Copyright (c) 2017 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include "./vpx_dsp_rtcd.h"
#include "vpx_ports/mem.h"
#include "vpx/vpx_integer.h"
#include "vpx_ports/asmdefs_mmi.h"

#define VARIANCE_SSE_8                                              \
  "gsldlc1    %[ftmp1],   0x07(%[a])                          \n\t" \
  "gsldrc1    %[ftmp1],   0x00(%[a])                          \n\t" \
  "gsldlc1    %[ftmp2],   0x07(%[b])                          \n\t" \
  "gsldrc1    %[ftmp2],   0x00(%[b])                          \n\t" \
  "pasubub    %[ftmp3],   %[ftmp1],       %[ftmp2]            \n\t" \
  "punpcklbh  %[ftmp4],   %[ftmp3],       %[ftmp0]            \n\t" \
  "punpckhbh  %[ftmp5],   %[ftmp3],       %[ftmp0]            \n\t" \
  "pmaddhw    %[ftmp6],   %[ftmp4],       %[ftmp4]            \n\t" \
  "pmaddhw    %[ftmp7],   %[ftmp5],       %[ftmp5]            \n\t" \
  "paddw      %[ftmp8],   %[ftmp8],       %[ftmp6]            \n\t" \
  "paddw      %[ftmp8],   %[ftmp8],       %[ftmp7]            \n\t"

#define VARIANCE_SSE_16                                             \
  VARIANCE_SSE_8                                                    \
  "gsldlc1    %[ftmp1],   0x0f(%[a])                          \n\t" \
  "gsldrc1    %[ftmp1],   0x08(%[a])                          \n\t" \
  "gsldlc1    %[ftmp2],   0x0f(%[b])                          \n\t" \
  "gsldrc1    %[ftmp2],   0x08(%[b])                          \n\t" \
  "pasubub    %[ftmp3],   %[ftmp1],       %[ftmp2]            \n\t" \
  "punpcklbh  %[ftmp4],   %[ftmp3],       %[ftmp0]            \n\t" \
  "punpckhbh  %[ftmp5],   %[ftmp3],       %[ftmp0]            \n\t" \
  "pmaddhw    %[ftmp6],   %[ftmp4],       %[ftmp4]            \n\t" \
  "pmaddhw    %[ftmp7],   %[ftmp5],       %[ftmp5]            \n\t" \
  "paddw      %[ftmp8],   %[ftmp8],       %[ftmp6]            \n\t" \
  "paddw      %[ftmp8],   %[ftmp8],       %[ftmp7]            \n\t"

static inline uint32_t vpx_mse16x(const uint8_t *a, int a_stride,
                                  const uint8_t *b, int b_stride, uint32_t *sse,
                                  uint64_t high) {
  double ftmp[12];
  uint32_t tmp[1];

  *sse = 0;

  __asm__ volatile (
    "li         %[tmp0],    0x20                                \n\t"
    "mtc1       %[tmp0],    %[ftmp11]                           \n\t"
    MMI_L(%[tmp0], %[high], 0x00)
    "xor        %[ftmp0],   %[ftmp0],       %[ftmp0]            \n\t"
    "xor        %[ftmp8],   %[ftmp8],       %[ftmp8]            \n\t"
    "xor        %[ftmp9],   %[ftmp9],       %[ftmp9]            \n\t"

    "1:                                                         \n\t"
    VARIANCE_SSE_16

    "addiu      %[tmp0],    %[tmp0],        -0x01               \n\t"
    MMI_ADDU(%[a], %[a], %[a_stride])
    MMI_ADDU(%[b], %[b], %[b_stride])
    "bnez       %[tmp0],    1b                                  \n\t"

    "dsrl       %[ftmp9],   %[ftmp8],       %[ftmp11]           \n\t"
    "paddw      %[ftmp9],   %[ftmp9],       %[ftmp8]            \n\t"
    "swc1       %[ftmp9],   0x00(%[sse])                        \n\t"
    : [ftmp0]"=&f"(ftmp[0]),            [ftmp1]"=&f"(ftmp[1]),
      [ftmp2]"=&f"(ftmp[2]),            [ftmp3]"=&f"(ftmp[3]),
      [ftmp4]"=&f"(ftmp[4]),            [ftmp5]"=&f"(ftmp[5]),
      [ftmp6]"=&f"(ftmp[6]),            [ftmp7]"=&f"(ftmp[7]),
      [ftmp8]"=&f"(ftmp[8]),            [ftmp9]"=&f"(ftmp[9]),
      [ftmp10]"=&f"(ftmp[10]),          [ftmp11]"=&f"(ftmp[11]),
      [tmp0]"=&r"(tmp[0]),
      [a]"+&r"(a),                      [b]"+&r"(b)
    : [a_stride]"r"((mips_reg)a_stride),[b_stride]"r"((mips_reg)b_stride),
      [high]"r"(&high), [sse]"r"(sse)
    : "memory"
  );

  return *sse;
}

#define vpx_mse16xN(n)                                         \
  uint32_t vpx_mse16x##n##_mmi(const uint8_t *a, int a_stride, \
                               const uint8_t *b, int b_stride, \
                               uint32_t *sse) {                \
    return vpx_mse16x(a, a_stride, b, b_stride, sse, n);       \
  }

vpx_mse16xN(16);
vpx_mse16xN(8);

static inline uint32_t vpx_mse8x(const uint8_t *a, int a_stride,
                                 const uint8_t *b, int b_stride, uint32_t *sse,
                                 uint64_t high) {
  double ftmp[12];
  uint32_t tmp[1];

  *sse = 0;

  __asm__ volatile (
    "li         %[tmp0],    0x20                                \n\t"
    "mtc1       %[tmp0],    %[ftmp11]                           \n\t"
    MMI_L(%[tmp0], %[high], 0x00)
    "xor        %[ftmp0],   %[ftmp0],       %[ftmp0]            \n\t"
    "xor        %[ftmp8],   %[ftmp8],       %[ftmp8]            \n\t"
    "xor        %[ftmp9],   %[ftmp9],       %[ftmp9]            \n\t"

    "1:                                                         \n\t"
    VARIANCE_SSE_8

    "addiu      %[tmp0],    %[tmp0],        -0x01               \n\t"
    MMI_ADDU(%[a], %[a], %[a_stride])
    MMI_ADDU(%[b], %[b], %[b_stride])
    "bnez       %[tmp0],    1b                                  \n\t"

    "dsrl       %[ftmp9],   %[ftmp8],       %[ftmp11]           \n\t"
    "paddw      %[ftmp9],   %[ftmp9],       %[ftmp8]            \n\t"
    "swc1       %[ftmp9],   0x00(%[sse])                        \n\t"
    : [ftmp0]"=&f"(ftmp[0]),            [ftmp1]"=&f"(ftmp[1]),
      [ftmp2]"=&f"(ftmp[2]),            [ftmp3]"=&f"(ftmp[3]),
      [ftmp4]"=&f"(ftmp[4]),            [ftmp5]"=&f"(ftmp[5]),
      [ftmp6]"=&f"(ftmp[6]),            [ftmp7]"=&f"(ftmp[7]),
      [ftmp8]"=&f"(ftmp[8]),            [ftmp9]"=&f"(ftmp[9]),
      [ftmp10]"=&f"(ftmp[10]),          [ftmp11]"=&f"(ftmp[11]),
      [tmp0]"=&r"(tmp[0]),
      [a]"+&r"(a),                      [b]"+&r"(b)
    : [a_stride]"r"((mips_reg)a_stride),[b_stride]"r"((mips_reg)b_stride),
      [high]"r"(&high), [sse]"r"(sse)
    : "memory"
  );

  return *sse;
}

#define vpx_mse8xN(n)                                                          \
  uint32_t vpx_mse8x##n##_mmi(const uint8_t *a, int a_stride,                  \
                              const uint8_t *b, int b_stride, uint32_t *sse) { \
    return vpx_mse8x(a, a_stride, b, b_stride, sse, n);                        \
  }

vpx_mse8xN(16);
vpx_mse8xN(8);
