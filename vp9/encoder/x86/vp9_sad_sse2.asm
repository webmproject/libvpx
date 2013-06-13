;
;  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
;
;  Use of this source code is governed by a BSD-style license
;  that can be found in the LICENSE file in the root of the source
;  tree. An additional intellectual property rights grant can be found
;  in the file PATENTS.  All contributing project authors may
;  be found in the AUTHORS file in the root of the source tree.
;

%include "third_party/x86inc/x86inc.asm"

SECTION .text

; unsigned int vp9_sad64x64_sse2(uint8_t *src, int src_stride,
;                                uint8_t *ref, int ref_stride);
%macro SAD64XN 1
cglobal sad64x%1, 4, 5, 5, src, src_stride, ref, ref_stride, n_rows
  movsxdifnidn src_strideq, src_strided
  movsxdifnidn ref_strideq, ref_strided
  mov              n_rowsd, %1
  pxor                  m0, m0
.loop:
  movu                  m1, [refq]
  movu                  m2, [refq+16]
  movu                  m3, [refq+32]
  movu                  m4, [refq+48]
  psadbw                m1, [srcq]
  psadbw                m2, [srcq+16]
  psadbw                m3, [srcq+32]
  psadbw                m4, [srcq+48]
  paddd                 m1, m2
  paddd                 m3, m4
  add                 refq, ref_strideq
  paddd                 m0, m1
  add                 srcq, src_strideq
  paddd                 m0, m3
  dec              n_rowsd
  jg .loop

  movhlps               m1, m0
  paddd                 m0, m1
  movd                 eax, m0
  RET
%endmacro

INIT_XMM sse2
SAD64XN 64 ; sad64x64_sse2
SAD64XN 32 ; sad64x32_sse2

; unsigned int vp9_sad32x32_sse2(uint8_t *src, int src_stride,
;                                uint8_t *ref, int ref_stride);
%macro SAD32XN 1
cglobal sad32x%1, 4, 5, 5, src, src_stride, ref, ref_stride, n_rows
  movsxdifnidn src_strideq, src_strided
  movsxdifnidn ref_strideq, ref_strided
  mov              n_rowsd, %1/2
  pxor                  m0, m0

.loop:
  movu                  m1, [refq]
  movu                  m2, [refq+16]
  movu                  m3, [refq+ref_strideq]
  movu                  m4, [refq+ref_strideq+16]
  psadbw                m1, [srcq]
  psadbw                m2, [srcq+16]
  psadbw                m3, [srcq+src_strideq]
  psadbw                m4, [srcq+src_strideq+16]
  paddd                 m1, m2
  paddd                 m3, m4
  lea                 refq, [refq+ref_strideq*2]
  paddd                 m0, m1
  lea                 srcq, [srcq+src_strideq*2]
  paddd                 m0, m3
  dec              n_rowsd
  jg .loop

  movhlps               m1, m0
  paddd                 m0, m1
  movd                 eax, m0
  RET
%endmacro

INIT_XMM sse2
SAD32XN 64 ; sad32x64_sse2
SAD32XN 32 ; sad32x32_sse2
SAD32XN 16 ; sad32x16_sse2

; unsigned int vp9_sad16x{8,16}_sse2(uint8_t *src, int src_stride,
;                                    uint8_t *ref, int ref_stride);
%macro SAD16XN 1
cglobal sad16x%1, 4, 7, 5, src, src_stride, ref, ref_stride, \
                           src_stride3, ref_stride3, n_rows
  movsxdifnidn src_strideq, src_strided
  movsxdifnidn ref_strideq, ref_strided
  lea         src_stride3q, [src_strideq*3]
  lea         ref_stride3q, [ref_strideq*3]
  mov              n_rowsd, %1/4
  pxor                  m0, m0

.loop:
  movu                  m1, [refq]
  movu                  m2, [refq+ref_strideq]
  movu                  m3, [refq+ref_strideq*2]
  movu                  m4, [refq+ref_stride3q]
  psadbw                m1, [srcq]
  psadbw                m2, [srcq+src_strideq]
  psadbw                m3, [srcq+src_strideq*2]
  psadbw                m4, [srcq+src_stride3q]
  paddd                 m1, m2
  paddd                 m3, m4
  lea                 refq, [refq+ref_strideq*4]
  paddd                 m0, m1
  lea                 srcq, [srcq+src_strideq*4]
  paddd                 m0, m3
  dec              n_rowsd
  jg .loop

  movhlps               m1, m0
  paddd                 m0, m1
  movd                 eax, m0
  RET
%endmacro

INIT_XMM sse2
SAD16XN 32 ; sad16x32_sse2
SAD16XN 16 ; sad16x16_sse2
SAD16XN  8 ; sad16x8_sse2

; unsigned int vp9_sad8x{8,16}_sse2(uint8_t *src, int src_stride,
;                                   uint8_t *ref, int ref_stride);
%macro SAD8XN 1
cglobal sad8x%1, 4, 7, 5, src, src_stride, ref, ref_stride, \
                          src_stride3, ref_stride3, n_rows
  movsxdifnidn src_strideq, src_strided
  movsxdifnidn ref_strideq, ref_strided
  lea         src_stride3q, [src_strideq*3]
  lea         ref_stride3q, [ref_strideq*3]
  mov              n_rowsd, %1/4
  pxor                  m0, m0

.loop:
  movh                  m1, [refq]
  movhps                m1, [refq+ref_strideq]
  movh                  m2, [refq+ref_strideq*2]
  movhps                m2, [refq+ref_stride3q]
  movh                  m3, [srcq]
  movhps                m3, [srcq+src_strideq]
  movh                  m4, [srcq+src_strideq*2]
  movhps                m4, [srcq+src_stride3q]
  psadbw                m1, m3
  psadbw                m2, m4
  lea                 refq, [refq+ref_strideq*4]
  paddd                 m0, m1
  lea                 srcq, [srcq+src_strideq*4]
  paddd                 m0, m2
  dec              n_rowsd
  jg .loop

  movhlps               m1, m0
  paddd                 m0, m1
  movd                 eax, m0
  RET
%endmacro

INIT_XMM sse2
SAD8XN 16 ; sad8x16_sse2
SAD8XN  8 ; sad8x8_sse2
SAD8XN  4 ; sad8x4_sse2

; unsigned int vp9_sad4x{4, 8}_sse(uint8_t *src, int src_stride,
;                                  uint8_t *ref, int ref_stride);
%macro SAD4XN 1
cglobal sad4x%1, 4, 7, 7, src, src_stride, ref, ref_stride, \
                          src_stride3, ref_stride3, n_rows
  movsxdifnidn src_strideq, src_strided
  movsxdifnidn ref_strideq, ref_strided
  lea         src_stride3q, [src_strideq*3]
  lea         ref_stride3q, [ref_strideq*3]
  mov              n_rowsd, %1/4
  pxor                  m0, m0

.loop:
  movd                  m1, [refq]
  movd                  m2, [refq+ref_strideq]
  movd                  m3, [refq+ref_strideq*2]
  movd                  m4, [refq+ref_stride3q]
  punpckldq             m1, m2
  punpckldq             m3, m4
  movd                  m2, [srcq]
  movd                  m5, [srcq+src_strideq]
  movd                  m4, [srcq+src_strideq*2]
  movd                  m6, [srcq+src_stride3q]
  punpckldq             m2, m5
  punpckldq             m4, m6
  psadbw                m1, m2
  psadbw                m3, m4
  lea                 refq, [refq+ref_strideq*4]
  paddd                 m0, m1
  lea                 srcq, [srcq+src_strideq*4]
  paddd                 m0, m3
  dec              n_rowsd
  jg .loop

  movd                 eax, m0
  RET
%endmacro

INIT_MMX sse
SAD4XN  8 ; sad4x8_sse
SAD4XN  4 ; sad4x4_sse
