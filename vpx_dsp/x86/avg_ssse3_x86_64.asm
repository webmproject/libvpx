;
;  Copyright (c) 2014 The WebM project authors. All Rights Reserved.
;
;  Use of this source code is governed by a BSD-style license
;  that can be found in the LICENSE file in the root of the source
;  tree. An additional intellectual property rights grant can be found
;  in the file PATENTS.  All contributing project authors may
;  be found in the AUTHORS file in the root of the source tree.
;

%define private_prefix vpx

%include "third_party/x86inc/x86inc.asm"

SECTION .text

%if ARCH_X86_64
; matrix transpose
%macro TRANSPOSE8X8 10
  ; stage 1
  punpcklwd  m%9, m%1, m%2
  punpcklwd  m%10, m%3, m%4
  punpckhwd  m%1, m%2
  punpckhwd  m%3, m%4

  punpcklwd  m%2, m%5, m%6
  punpcklwd  m%4, m%7, m%8
  punpckhwd  m%5, m%6
  punpckhwd  m%7, m%8

  ; stage 2
  punpckldq  m%6, m%9, m%10
  punpckldq  m%8, m%1, m%3
  punpckhdq  m%9, m%10
  punpckhdq  m%1, m%3

  punpckldq  m%10, m%2, m%4
  punpckldq  m%3, m%5, m%7
  punpckhdq  m%2, m%4
  punpckhdq  m%5, m%7

  ; stage 3
  punpckhqdq  m%4, m%9, m%2  ; out3
  punpcklqdq  m%9, m%2       ; out2
  punpcklqdq  m%7, m%1, m%5  ; out6
  punpckhqdq  m%1, m%5       ; out7

  punpckhqdq  m%2, m%6, m%10 ; out1
  punpcklqdq  m%6, m%10      ; out0
  punpcklqdq  m%5, m%8, m%3  ; out4
  punpckhqdq  m%8, m%3       ; out5

  SWAP %6, %1
  SWAP %3, %9
  SWAP %8, %6
%endmacro

%macro HMD8_1D 0
  psubw              m8, m0, m1
  psubw              m9, m2, m3
  paddw              m0, m1
  paddw              m2, m3
  SWAP               1, 8
  SWAP               3, 9
  psubw              m8, m4, m5
  psubw              m9, m6, m7
  paddw              m4, m5
  paddw              m6, m7
  SWAP               5, 8
  SWAP               7, 9

  psubw              m8, m0, m2
  psubw              m9, m1, m3
  paddw              m0, m2
  paddw              m1, m3
  SWAP               2, 8
  SWAP               3, 9
  psubw              m8, m4, m6
  psubw              m9, m5, m7
  paddw              m4, m6
  paddw              m5, m7
  SWAP               6, 8
  SWAP               7, 9

  psubw              m8, m0, m4
  psubw              m9, m1, m5
  paddw              m0, m4
  paddw              m1, m5
  SWAP               4, 8
  SWAP               5, 9
  psubw              m8, m2, m6
  psubw              m9, m3, m7
  paddw              m2, m6
  paddw              m3, m7
  SWAP               6, 8
  SWAP               7, 9
%endmacro

INIT_XMM ssse3
cglobal hadamard_8x8, 3, 5, 11, input, stride, output
  lea                r3, [2 * strideq]
  lea                r4, [4 * strideq]

  mova               m0, [inputq]
  mova               m1, [inputq + r3]
  lea                inputq, [inputq + r4]
  mova               m2, [inputq]
  mova               m3, [inputq + r3]
  lea                inputq, [inputq + r4]
  mova               m4, [inputq]
  mova               m5, [inputq + r3]
  lea                inputq, [inputq + r4]
  mova               m6, [inputq]
  mova               m7, [inputq + r3]

  HMD8_1D
  TRANSPOSE8X8 0, 1, 2, 3, 4, 5, 6, 7, 9, 10
  HMD8_1D

  mova              [outputq +   0], m0
  mova              [outputq +  16], m1
  mova              [outputq +  32], m2
  mova              [outputq +  48], m3
  mova              [outputq +  64], m4
  mova              [outputq +  80], m5
  mova              [outputq +  96], m6
  mova              [outputq + 112], m7

  RET
%endif
