;
;  Copyright (c) 2014 The WebM project authors. All Rights Reserved.
;
;  Use of this source code is governed by a BSD-style license
;  that can be found in the LICENSE file in the root of the source
;  tree. An additional intellectual property rights grant can be found
;  in the file PATENTS.  All contributing project authors may
;  be found in the AUTHORS file in the root of the source tree.
;

%include "third_party/x86inc/x86inc.asm"
%include "vpx_dsp/x86/bitdepth_conversion_sse2.asm"

SECTION_RODATA

pw_11585x2: times 8 dw 23170

pw_m2404x2:  times 8 dw  -2404*2
pw_m4756x2:  times 8 dw  -4756*2
pw_m5520x2:  times 8 dw  -5520*2
pw_m8423x2:  times 8 dw  -8423*2
pw_m9102x2:  times 8 dw  -9102*2
pw_m10394x2: times 8 dw -10394*2
pw_m11003x2: times 8 dw -11003*2

pw_16364x2: times 8 dw 16364*2
pw_16305x2: times 8 dw 16305*2
pw_16207x2: times 8 dw 16207*2
pw_16069x2: times 8 dw 16069*2
pw_15893x2: times 8 dw 15893*2
pw_15679x2: times 8 dw 15679*2
pw_15426x2: times 8 dw 15426*2
pw_15137x2: times 8 dw 15137*2
pw_14811x2: times 8 dw 14811*2
pw_14449x2: times 8 dw 14449*2
pw_14053x2: times 8 dw 14053*2
pw_13623x2: times 8 dw 13623*2
pw_13160x2: times 8 dw 13160*2
pw_12665x2: times 8 dw 12665*2
pw_12140x2: times 8 dw 12140*2
pw__9760x2: times 8 dw  9760*2
pw__7723x2: times 8 dw  7723*2
pw__7005x2: times 8 dw  7005*2
pw__6270x2: times 8 dw  6270*2
pw__3981x2: times 8 dw  3981*2
pw__3196x2: times 8 dw  3196*2
pw__1606x2: times 8 dw  1606*2
pw___804x2: times 8 dw   804*2

pd_8192:    times 4 dd 8192
pw_32:      times 8 dw 32
pw_16:      times 8 dw 16

%macro TRANSFORM_COEFFS 2
pw_%1_%2:   dw  %1,  %2,  %1,  %2,  %1,  %2,  %1,  %2
pw_m%2_%1:  dw -%2,  %1, -%2,  %1, -%2,  %1, -%2,  %1
pw_m%1_m%2: dw -%1, -%2, -%1, -%2, -%1, -%2, -%1, -%2
%endmacro

TRANSFORM_COEFFS    6270, 15137
TRANSFORM_COEFFS    3196, 16069
TRANSFORM_COEFFS   13623,  9102

; constants for 32x32_34
TRANSFORM_COEFFS      804, 16364
TRANSFORM_COEFFS    15426,  5520
TRANSFORM_COEFFS     3981, 15893
TRANSFORM_COEFFS    16207,  2404
TRANSFORM_COEFFS     1606, 16305
TRANSFORM_COEFFS    15679,  4756
TRANSFORM_COEFFS    11585, 11585

; constants for 32x32_1024
TRANSFORM_COEFFS    12140, 11003
TRANSFORM_COEFFS     7005, 14811
TRANSFORM_COEFFS    14053,  8423
TRANSFORM_COEFFS     9760, 13160
TRANSFORM_COEFFS    12665, 10394
TRANSFORM_COEFFS     7723, 14449

%macro PAIR_PP_COEFFS 2
dpw_%1_%2:   dw  %1,  %1,  %1,  %1,  %2,  %2,  %2,  %2
%endmacro

%macro PAIR_MP_COEFFS 2
dpw_m%1_%2:  dw -%1, -%1, -%1, -%1,  %2,  %2,  %2,  %2
%endmacro

%macro PAIR_MM_COEFFS 2
dpw_m%1_m%2: dw -%1, -%1, -%1, -%1, -%2, -%2, -%2, -%2
%endmacro

PAIR_PP_COEFFS     30274, 12540
PAIR_PP_COEFFS      6392, 32138
PAIR_MP_COEFFS     18204, 27246

PAIR_PP_COEFFS     12540, 12540
PAIR_PP_COEFFS     30274, 30274
PAIR_PP_COEFFS      6392,  6392
PAIR_PP_COEFFS     32138, 32138
PAIR_MM_COEFFS     18204, 18204
PAIR_PP_COEFFS     27246, 27246

SECTION .text

%if ARCH_X86_64
%macro SUM_SUB 3
  psubw  m%3, m%1, m%2
  paddw  m%1, m%2
  SWAP    %2, %3
%endmacro

; butterfly operation
%macro MUL_ADD_2X 6 ; dst1, dst2, src, round, coefs1, coefs2
  pmaddwd            m%1, m%3, %5
  pmaddwd            m%2, m%3, %6
  paddd              m%1,  %4
  paddd              m%2,  %4
  psrad              m%1,  14
  psrad              m%2,  14
%endmacro

%macro BUTTERFLY_4X 7 ; dst1, dst2, coef1, coef2, round, tmp1, tmp2
  punpckhwd          m%6, m%2, m%1
  MUL_ADD_2X         %7,  %6,  %6,  %5, [pw_m%4_%3], [pw_%3_%4]
  punpcklwd          m%2, m%1
  MUL_ADD_2X         %1,  %2,  %2,  %5, [pw_m%4_%3], [pw_%3_%4]
  packssdw           m%1, m%7
  packssdw           m%2, m%6
%endmacro

%macro BUTTERFLY_4Xmm 7 ; dst1, dst2, coef1, coef2, round, tmp1, tmp2
  punpckhwd          m%6, m%2, m%1
  MUL_ADD_2X         %7,  %6,  %6,  %5, [pw_m%4_%3], [pw_m%3_m%4]
  punpcklwd          m%2, m%1
  MUL_ADD_2X         %1,  %2,  %2,  %5, [pw_m%4_%3], [pw_m%3_m%4]
  packssdw           m%1, m%7
  packssdw           m%2, m%6
%endmacro

; 8x8 transpose. This follows same operations as SSE2 does.
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

%macro IDCT8_1D 0
  SUM_SUB          0,    4,    9
  BUTTERFLY_4X     2,    6,    6270, 15137,  m8,  9,  10
  pmulhrsw        m0,  m12
  pmulhrsw        m4,  m12
  BUTTERFLY_4X     1,    7,    3196, 16069,  m8,  9,  10
  BUTTERFLY_4X     5,    3,   13623,  9102,  m8,  9,  10

  SUM_SUB          1,    5,    9
  SUM_SUB          7,    3,    9
  SUM_SUB          0,    6,    9
  SUM_SUB          4,    2,    9
  SUM_SUB          3,    5,    9
  pmulhrsw        m3,  m12
  pmulhrsw        m5,  m12

  SUM_SUB          0,    7,    9
  SUM_SUB          4,    3,    9
  SUM_SUB          2,    5,    9
  SUM_SUB          6,    1,    9

  SWAP             3,    6
  SWAP             1,    4
%endmacro

; This macro handles 8 pixels per line
%macro ADD_STORE_8P_2X 5;  src1, src2, tmp1, tmp2, zero
  paddw           m%1, m11
  paddw           m%2, m11
  psraw           m%1, 5
  psraw           m%2, 5

  movh            m%3, [outputq]
  movh            m%4, [outputq + strideq]
  punpcklbw       m%3, m%5
  punpcklbw       m%4, m%5
  paddw           m%3, m%1
  paddw           m%4, m%2
  packuswb        m%3, m%5
  packuswb        m%4, m%5
  movh               [outputq], m%3
  movh     [outputq + strideq], m%4
%endmacro

%define  idx0 16 * 0
%define  idx1 16 * 1
%define  idx2 16 * 2
%define  idx3 16 * 3
%define  idx4 16 * 4
%define  idx5 16 * 5
%define  idx6 16 * 6
%define  idx7 16 * 7
%define  idx8 16 * 0
%define  idx9 16 * 1
%define idx10 16 * 2
%define idx11 16 * 3
%define idx12 16 * 4
%define idx13 16 * 5
%define idx14 16 * 6
%define idx15 16 * 7
%define idx16 16 * 0
%define idx17 16 * 1
%define idx18 16 * 2
%define idx19 16 * 3
%define idx20 16 * 4
%define idx21 16 * 5
%define idx22 16 * 6
%define idx23 16 * 7
%define idx24 16 * 0
%define idx25 16 * 1
%define idx26 16 * 2
%define idx27 16 * 3
%define idx28 16 * 4
%define idx29 16 * 5
%define idx30 16 * 6
%define idx31 16 * 7

; FROM idct32x32_add_neon.asm
;
; Instead of doing the transforms stage by stage, it is done by loading
; some input values and doing as many stages as possible to minimize the
; storing/loading of intermediate results. To fit within registers, the
; final coefficients are cut into four blocks:
; BLOCK A: 16-19,28-31
; BLOCK B: 20-23,24-27
; BLOCK C: 8-11,12-15
; BLOCK D: 0-3,4-7
; Blocks A and C are straight calculation through the various stages. In
; block B, further calculations are performed using the results from
; block A. In block D, further calculations are performed using the results
; from block C and then the final calculations are done using results from
; block A and B which have been combined at the end of block B.
;

%macro IDCT32X32_34 4
  ; BLOCK A STAGE 1 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  mova                m11, m1
  pmulhrsw             m1, [pw___804x2] ; stp1_16
  mova      [r4 +      0], m0
  pmulhrsw            m11, [pw_16364x2] ; stp2_31
  mova      [r4 + 16 * 2], m2
  mova                m12, m7
  pmulhrsw             m7, [pw_15426x2] ; stp1_28
  mova      [r4 + 16 * 4], m4
  pmulhrsw            m12, [pw_m5520x2] ; stp2_19
  mova      [r4 + 16 * 6], m6

  ; BLOCK A STAGE 2 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  mova                 m2, m1   ; stp1_16
  mova                 m0, m11  ; stp1_31
  mova                 m4, m7   ; stp1_28
  mova                m15, m12  ; stp1_19

  ; BLOCK A STAGE 3 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  BUTTERFLY_4X          0,     2,   3196, 16069,  m8,  9,  10 ; stp1_17, stp1_30
  BUTTERFLY_4Xmm        4,    15,   3196, 16069,  m8,  9,  10 ; stp1_29, stp1_18

  ; BLOCK A STAGE 4 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUM_SUB               1, 12, 9 ; stp2_16, stp2_19
  SUM_SUB               0, 15, 9 ; stp2_17, stp2_18
  SUM_SUB              11,  7, 9 ; stp2_31, stp2_28
  SUM_SUB               2,  4, 9 ; stp2_30, stp2_29

  ; BLOCK A STAGE 5 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  BUTTERFLY_4X          4,    15,   6270, 15137,  m8,  9,  10 ; stp1_18, stp1_29
  BUTTERFLY_4X          7,    12,   6270, 15137,  m8,  9,  10 ; stp1_19, stp1_28

  ; BLOCK B STAGE 1 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  mova                 m6, m5
  pmulhrsw             m5, [pw__3981x2] ; stp1_20
  mova [stp + %4 + idx28], m12
  mova [stp + %4 + idx29], m15
  pmulhrsw             m6, [pw_15893x2] ; stp2_27
  mova [stp + %4 + idx30], m2
  mova                 m2, m3
  pmulhrsw             m3, [pw_m2404x2] ; stp1_23
  mova [stp + %4 + idx31], m11
  pmulhrsw             m2, [pw_16207x2] ; stp2_24

  ; BLOCK B STAGE 2 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  mova                m13, m5 ; stp1_20
  mova                m14, m6 ; stp1_27
  mova                m15, m3 ; stp1_23
  mova                m11, m2 ; stp1_24

  ; BLOCK B STAGE 3 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  BUTTERFLY_4X         14,    13,  13623,  9102,  m8,  9,  10 ; stp1_21, stp1_26
  BUTTERFLY_4Xmm       11,    15,  13623,  9102,  m8,  9,  10 ; stp1_25, stp1_22

  ; BLOCK B STAGE 4 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUM_SUB               3,  5, 9 ; stp2_23, stp2_20
  SUM_SUB              15, 14, 9 ; stp2_22, stp2_21
  SUM_SUB               2,  6, 9 ; stp2_24, stp2_27
  SUM_SUB              11, 13, 9 ; stp2_25, stp2_26

  ; BLOCK B STAGE 5 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  BUTTERFLY_4Xmm        6,     5,   6270, 15137,  m8,  9,  10 ; stp1_27, stp1_20
  BUTTERFLY_4Xmm       13,    14,   6270, 15137,  m8,  9,  10 ; stp1_26, stp1_21

  ; BLOCK B STAGE 6 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUM_SUB               1,  3, 9 ; stp2_16, stp2_23
  SUM_SUB               0, 15, 9 ; stp2_17, stp2_22
  SUM_SUB               4, 14, 9 ; stp2_18, stp2_21
  SUM_SUB               7,  5, 9 ; stp2_19, stp2_20
  mova [stp + %3 + idx16], m1
  mova [stp + %3 + idx17], m0
  mova [stp + %3 + idx18], m4
  mova [stp + %3 + idx19], m7

  mova                 m4, [stp + %4 + idx28]
  mova                 m7, [stp + %4 + idx29]
  mova                m10, [stp + %4 + idx30]
  mova                m12, [stp + %4 + idx31]
  SUM_SUB               4,  6, 9 ; stp2_28, stp2_27
  SUM_SUB               7, 13, 9 ; stp2_29, stp2_26
  SUM_SUB              10, 11, 9 ; stp2_30, stp2_25
  SUM_SUB              12,  2, 9 ; stp2_31, stp2_24
  mova [stp + %4 + idx28], m4
  mova [stp + %4 + idx29], m7
  mova [stp + %4 + idx30], m10
  mova [stp + %4 + idx31], m12

  ; BLOCK B STAGE 7 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%if 0 ; overflow occurs in SUM_SUB when using test streams
  mova                m10, [pw_11585x2]
  SUM_SUB               6,  5, 9
  pmulhrsw             m6, m10  ; stp1_27
  pmulhrsw             m5, m10  ; stp1_20
  SUM_SUB              13, 14,  9
  pmulhrsw            m13, m10  ; stp1_26
  pmulhrsw            m14, m10  ; stp1_21
  SUM_SUB              11, 15,  9
  pmulhrsw            m11, m10  ; stp1_25
  pmulhrsw            m15, m10  ; stp1_22
  SUM_SUB               2,  3,  9
  pmulhrsw             m2, m10  ; stp1_24
  pmulhrsw             m3, m10  ; stp1_23
%else
  BUTTERFLY_4X          6,     5,  11585, 11585,  m8,  9,  10 ; stp1_20, stp1_27
  SWAP 6, 5
  BUTTERFLY_4X         13,    14,  11585, 11585,  m8,  9,  10 ; stp1_21, stp1_26
  SWAP 13, 14
  BUTTERFLY_4X         11,    15,  11585, 11585,  m8,  9,  10 ; stp1_22, stp1_25
  SWAP 11, 15
  BUTTERFLY_4X          2,     3,  11585, 11585,  m8,  9,  10 ; stp1_23, stp1_24
  SWAP 2, 3
%endif

  mova [stp + %4 + idx24], m2
  mova [stp + %4 + idx25], m11
  mova [stp + %4 + idx26], m13
  mova [stp + %4 + idx27], m6

  ; BLOCK C STAGE 1 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ;
  ; BLOCK C STAGE 2 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  mova                 m0, [rsp + transposed_in + 16 *  2]
  mova                 m6, [rsp + transposed_in + 16 *  6]

  mova                 m1, m0
  pmulhrsw             m0, [pw__1606x2] ; stp1_8
  mova [stp + %3 + idx20], m5
  mova [stp + %3 + idx21], m14
  pmulhrsw             m1, [pw_16305x2] ; stp2_15
  mova [stp + %3 + idx22], m15
  mova                 m7, m6
  pmulhrsw             m7, [pw_m4756x2] ; stp2_11
  mova [stp + %3 + idx23], m3
  pmulhrsw             m6, [pw_15679x2] ; stp1_12

  ; BLOCK C STAGE 3 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  mova                 m3, m0 ; stp1_8
  mova                 m2, m1 ; stp1_15

  ; BLOCK C STAGE 4 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  BUTTERFLY_4X          2,     3,   6270, 15137,  m8,  9,  10 ;  stp1_9, stp1_14
  mova                 m4, m7 ; stp1_11
  mova                 m5, m6 ; stp1_12
  BUTTERFLY_4Xmm        5,     4,   6270, 15137,  m8,  9,  10 ; stp1_13, stp1_10

  ; BLOCK C STAGE 5 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUM_SUB               0,  7, 9 ;  stp1_8, stp1_11
  SUM_SUB               2,  4, 9 ;  stp1_9, stp1_10
  SUM_SUB               1,  6, 9 ;  stp1_15, stp1_12
  SUM_SUB               3,  5, 9 ;  stp1_14, stp1_13

  ; BLOCK C STAGE 6 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%if 0 ; overflow occurs in SUM_SUB when using test streams
  mova                m10, [pw_11585x2]
  SUM_SUB               5,  4, 9
  pmulhrsw             m5, m10  ; stp1_13
  pmulhrsw             m4, m10  ; stp1_10
  SUM_SUB               6,  7, 9
  pmulhrsw             m6, m10  ; stp1_12
  pmulhrsw             m7, m10  ; stp1_11
%else
  BUTTERFLY_4X          5,     4,  11585, 11585,  m8,  9,  10 ; stp1_10, stp1_13
  SWAP 5, 4
  BUTTERFLY_4X          6,     7,  11585, 11585,  m8,  9,  10 ; stp1_11, stp1_12
  SWAP 6, 7
%endif

  ; BLOCK C STAGE 7 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  mova [stp + %2 +  idx8], m0
  mova [stp + %2 +  idx9], m2
  mova [stp + %2 + idx10], m4
  mova [stp + %2 + idx11], m7

  ; BLOCK D STAGE 1 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ;
  ; BLOCK D STAGE 2 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ;
  ; BLOCK D STAGE 3 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  mova                m11, [rsp + transposed_in + 16 *  4]
  mova                m12, m11
  pmulhrsw            m11, [pw__3196x2] ; stp1_4
  pmulhrsw            m12, [pw_16069x2] ; stp1_7

  ; BLOCK D STAGE 4 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  mova                 m0, [rsp + transposed_in + 16 *  0]
  mova                m10, [pw_11585x2]
  pmulhrsw             m0, m10  ; stp1_1

  mova                m14, m11 ; stp1_4
  mova                m13, m12 ; stp1_7

  ; BLOCK D STAGE 5 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%if 0 ; overflow occurs in SUM_SUB when using test streams
  SUM_SUB              13,   14,  9
  pmulhrsw            m13, m10  ; stp1_6
  pmulhrsw            m14, m10  ; stp1_5
%else
  BUTTERFLY_4X         13,    14,  11585, 11585,  m8,  9,  10 ; stp1_5, stp1_6
  SWAP 13, 14
%endif
  mova                 m7, m0 ; stp1_0 = stp1_1
  mova                 m4, m0 ; stp1_1
  mova                 m2, m7 ; stp1_0

  ; BLOCK D STAGE 6 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUM_SUB               0, 12, 9 ;  stp1_0, stp1_7
  SUM_SUB               7, 13, 9 ;  stp1_1, stp1_6
  SUM_SUB               2, 14, 9 ;  stp1_2, stp1_5
  SUM_SUB               4, 11, 9 ;  stp1_3, stp1_4

  ; BLOCK D STAGE 7 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUM_SUB               0,  1, 9 ;  stp1_0, stp1_15
  SUM_SUB               7,  3, 9 ;  stp1_1, stp1_14
  SUM_SUB               2,  5, 9 ;  stp1_2, stp1_13
  SUM_SUB               4,  6, 9 ;  stp1_3, stp1_12

  ; 0-3, 28-31 final stage
  mova                m15, [stp + %4 + idx30]
  mova                m10, [stp + %4 + idx31]
  SUM_SUB               0, 10, 9 ;  stp1_0, stp1_31
  SUM_SUB               7, 15, 9 ;  stp1_1, stp1_30
  mova [stp + %1 +  idx0], m0
  mova [stp + %1 +  idx1], m7
  mova [stp + %4 + idx30], m15
  mova [stp + %4 + idx31], m10
  mova                 m7, [stp + %4 + idx28]
  mova                 m0, [stp + %4 + idx29]
  SUM_SUB               2,  0, 9 ;  stp1_2, stp1_29
  SUM_SUB               4,  7, 9 ;  stp1_3, stp1_28
  mova [stp + %1 +  idx2], m2
  mova [stp + %1 +  idx3], m4
  mova [stp + %4 + idx28], m7
  mova [stp + %4 + idx29], m0

  ; 12-15, 16-19 final stage
  mova                 m0, [stp + %3 + idx16]
  mova                 m7, [stp + %3 + idx17]
  mova                 m2, [stp + %3 + idx18]
  mova                 m4, [stp + %3 + idx19]
  SUM_SUB               1,  0, 9 ;  stp1_15, stp1_16
  SUM_SUB               3,  7, 9 ;  stp1_14, stp1_17
  SUM_SUB               5,  2, 9 ;  stp1_13, stp1_18
  SUM_SUB               6,  4, 9 ;  stp1_12, stp1_19
  mova [stp + %2 + idx12], m6
  mova [stp + %2 + idx13], m5
  mova [stp + %2 + idx14], m3
  mova [stp + %2 + idx15], m1
  mova [stp + %3 + idx16], m0
  mova [stp + %3 + idx17], m7
  mova [stp + %3 + idx18], m2
  mova [stp + %3 + idx19], m4

  mova                 m4, [stp + %2 +  idx8]
  mova                 m5, [stp + %2 +  idx9]
  mova                 m6, [stp + %2 + idx10]
  mova                 m7, [stp + %2 + idx11]
  SUM_SUB              11,  7, 9 ;  stp1_4, stp1_11
  SUM_SUB              14,  6, 9 ;  stp1_5, stp1_10
  SUM_SUB              13,  5, 9 ;  stp1_6, stp1_9
  SUM_SUB              12,  4, 9 ;  stp1_7, stp1_8

  ; 4-7, 24-27 final stage
  mova                 m0, [stp + %4 + idx27]
  mova                 m1, [stp + %4 + idx26]
  mova                 m2, [stp + %4 + idx25]
  mova                 m3, [stp + %4 + idx24]
  SUM_SUB              11,  0, 9 ;  stp1_4, stp1_27
  SUM_SUB              14,  1, 9 ;  stp1_5, stp1_26
  SUM_SUB              13,  2, 9 ;  stp1_6, stp1_25
  SUM_SUB              12,  3, 9 ;  stp1_7, stp1_24
  mova [stp + %4 + idx27], m0
  mova [stp + %4 + idx26], m1
  mova [stp + %4 + idx25], m2
  mova [stp + %4 + idx24], m3
  mova [stp + %1 +  idx4], m11
  mova [stp + %1 +  idx5], m14
  mova [stp + %1 +  idx6], m13
  mova [stp + %1 +  idx7], m12

  ; 8-11, 20-23 final stage
  mova                 m0, [stp + %3 + idx20]
  mova                 m1, [stp + %3 + idx21]
  mova                 m2, [stp + %3 + idx22]
  mova                 m3, [stp + %3 + idx23]
  SUM_SUB               7,  0, 9 ;  stp1_11, stp_20
  SUM_SUB               6,  1, 9 ;  stp1_10, stp_21
  SUM_SUB               5,  2, 9 ;   stp1_9, stp_22
  SUM_SUB               4,  3, 9 ;   stp1_8, stp_23
  mova [stp + %2 +  idx8], m4
  mova [stp + %2 +  idx9], m5
  mova [stp + %2 + idx10], m6
  mova [stp + %2 + idx11], m7
  mova [stp + %3 + idx20], m0
  mova [stp + %3 + idx21], m1
  mova [stp + %3 + idx22], m2
  mova [stp + %3 + idx23], m3
%endmacro

%macro RECON_AND_STORE 1
  mova            m11, [pw_32]
  lea             stp, [rsp + %1]
  mov              r6, 32
  pxor             m8, m8
%%recon_and_store:
  mova             m0, [stp + 16 * 32 * 0]
  mova             m1, [stp + 16 * 32 * 1]
  mova             m2, [stp + 16 * 32 * 2]
  mova             m3, [stp + 16 * 32 * 3]
  add             stp, 16

  paddw            m0, m11
  paddw            m1, m11
  paddw            m2, m11
  paddw            m3, m11
  psraw            m0, 6
  psraw            m1, 6
  psraw            m2, 6
  psraw            m3, 6
  movh             m4, [outputq +  0]
  movh             m5, [outputq +  8]
  movh             m6, [outputq + 16]
  movh             m7, [outputq + 24]
  punpcklbw        m4, m8
  punpcklbw        m5, m8
  punpcklbw        m6, m8
  punpcklbw        m7, m8
  paddw            m0, m4
  paddw            m1, m5
  paddw            m2, m6
  paddw            m3, m7
  packuswb         m0, m1
  packuswb         m2, m3
  mova [outputq +  0], m0
  mova [outputq + 16], m2
  lea         outputq, [outputq + strideq]
  dec              r6
  jnz %%recon_and_store
%endmacro

%define i32x32_size     16*32*5
%define pass_two_start  16*32*0
%define transposed_in   16*32*4
%define pass_one_start  16*32*0
%define stp r8

INIT_XMM ssse3
cglobal idct32x32_34_add, 3, 11, 16, i32x32_size, input, output, stride
  mova            m8, [pd_8192]
  lea            stp, [rsp + pass_one_start]

idct32x32_34:
  mov             r3, inputq
  lea             r4, [rsp + transposed_in]

idct32x32_34_transpose:
  LOAD_TRAN_LOW 0, r3,  0
  LOAD_TRAN_LOW 1, r3,  4
  LOAD_TRAN_LOW 2, r3,  8
  LOAD_TRAN_LOW 3, r3, 12
  LOAD_TRAN_LOW 4, r3, 16
  LOAD_TRAN_LOW 5, r3, 20
  LOAD_TRAN_LOW 6, r3, 24
  LOAD_TRAN_LOW 7, r3, 28

  TRANSPOSE8X8  0, 1, 2, 3, 4, 5, 6, 7, 9, 10

  IDCT32X32_34  16*0, 16*32, 16*64, 16*96
  lea            stp, [stp + 16 * 8]
  mov             r6, 4
  lea            stp, [rsp + pass_one_start]
  lea             r9, [rsp + pass_one_start]

idct32x32_34_2:
  lea             r4, [rsp + transposed_in]
  mov             r3, r9

idct32x32_34_transpose_2:
  mova            m0, [r3 +      0]
  mova            m1, [r3 + 16 * 1]
  mova            m2, [r3 + 16 * 2]
  mova            m3, [r3 + 16 * 3]
  mova            m4, [r3 + 16 * 4]
  mova            m5, [r3 + 16 * 5]
  mova            m6, [r3 + 16 * 6]
  mova            m7, [r3 + 16 * 7]

  TRANSPOSE8X8  0, 1, 2, 3, 4, 5, 6, 7, 9, 10

  IDCT32X32_34  16*0, 16*8, 16*16, 16*24

  lea            stp, [stp + 16 * 32]
  add             r9, 16 * 32
  dec             r6
  jnz idct32x32_34_2

  RECON_AND_STORE pass_two_start

  RET

%macro IDCT32X32_135 4
  ; BLOCK A STAGE 1 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  mova                 m1, [rsp + transposed_in + 16 *  1]
  mova                m11, m1
  pmulhrsw             m1, [pw___804x2] ; stp1_16
  pmulhrsw            m11, [pw_16364x2] ; stp2_31

  mova                 m7, [rsp + transposed_in + 16 *  7]
  mova                m12, m7
  pmulhrsw             m7, [pw_15426x2] ; stp1_28
  pmulhrsw            m12, [pw_m5520x2] ; stp2_19

  mova                 m3, [rsp + transposed_in + 16 *  9]
  mova                 m4, m3
  pmulhrsw             m3, [pw__7005x2] ; stp1_18
  pmulhrsw             m4, [pw_14811x2] ; stp2_29

  mova                 m0, [rsp + transposed_in + 16 * 15]
  mova                 m2, m0
  pmulhrsw             m0, [pw_12140x2]  ; stp1_30
  pmulhrsw             m2, [pw_m11003x2] ; stp2_17

  ; BLOCK A STAGE 2 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUM_SUB               1,  2, 9 ; stp2_16, stp2_17
  SUM_SUB              12,  3, 9 ; stp2_19, stp2_18
  SUM_SUB               7,  4, 9 ; stp2_28, stp2_29
  SUM_SUB              11,  0, 9 ; stp2_31, stp2_30

  ; BLOCK A STAGE 3 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  BUTTERFLY_4X          0,     2,   3196, 16069,  m8,  9,  10 ; stp1_17, stp1_30
  BUTTERFLY_4Xmm        4,     3,   3196, 16069,  m8,  9,  10 ; stp1_29, stp1_18

  ; BLOCK A STAGE 4 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUM_SUB               1, 12, 9 ; stp2_16, stp2_19
  SUM_SUB               0,  3, 9 ; stp2_17, stp2_18
  SUM_SUB              11,  7, 9 ; stp2_31, stp2_28
  SUM_SUB               2,  4, 9 ; stp2_30, stp2_29

  ; BLOCK A STAGE 5 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  BUTTERFLY_4X          4,     3,   6270, 15137,  m8,  9,  10 ; stp1_18, stp1_29
  BUTTERFLY_4X          7,    12,   6270, 15137,  m8,  9,  10 ; stp1_19, stp1_28

  mova [stp + %3 + idx16], m1
  mova [stp + %3 + idx17], m0
  mova [stp + %3 + idx18], m4
  mova [stp + %3 + idx19], m7
  mova [stp + %4 + idx28], m12
  mova [stp + %4 + idx29], m3
  mova [stp + %4 + idx30], m2
  mova [stp + %4 + idx31], m11

  ; BLOCK B STAGE 1 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  mova                 m2, [rsp + transposed_in + 16 *  3]
  mova                 m3, m2
  pmulhrsw             m3, [pw_m2404x2] ; stp1_23
  pmulhrsw             m2, [pw_16207x2] ; stp2_24

  mova                 m5, [rsp + transposed_in + 16 *  5]
  mova                 m6, m5
  pmulhrsw             m5, [pw__3981x2] ; stp1_20
  pmulhrsw             m6, [pw_15893x2] ; stp2_27

  mova                m14, [rsp + transposed_in + 16 * 11]
  mova                m13, m14
  pmulhrsw            m13, [pw_m8423x2] ; stp1_21
  pmulhrsw            m14, [pw_14053x2] ; stp2_26

  mova                 m0, [rsp + transposed_in + 16 * 13]
  mova                 m1, m0
  pmulhrsw             m0, [pw__9760x2] ; stp1_22
  pmulhrsw             m1, [pw_13160x2] ; stp2_25

  ; BLOCK B STAGE 2 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUM_SUB               5, 13, 9 ; stp2_20, stp2_21
  SUM_SUB               3,  0, 9 ; stp2_23, stp2_22
  SUM_SUB               2,  1, 9 ; stp2_24, stp2_25
  SUM_SUB               6, 14, 9 ; stp2_27, stp2_26

  ; BLOCK B STAGE 3 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  BUTTERFLY_4X         14,    13,  13623,  9102,  m8,  9,  10 ; stp1_21, stp1_26
  BUTTERFLY_4Xmm        1,     0,  13623,  9102,  m8,  9,  10 ; stp1_25, stp1_22

  ; BLOCK B STAGE 4 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUM_SUB               3,  5, 9 ; stp2_23, stp2_20
  SUM_SUB               0, 14, 9 ; stp2_22, stp2_21
  SUM_SUB               2,  6, 9 ; stp2_24, stp2_27
  SUM_SUB               1, 13, 9 ; stp2_25, stp2_26

  ; BLOCK B STAGE 5 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  BUTTERFLY_4Xmm        6,     5,   6270, 15137,  m8,  9,  10 ; stp1_27, stp1_20
  BUTTERFLY_4Xmm       13,    14,   6270, 15137,  m8,  9,  10 ; stp1_26, stp1_21

  ; BLOCK B STAGE 6 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  mova                 m4, [stp + %3 + idx16]
  mova                 m7, [stp + %3 + idx17]
  mova                m11, [stp + %3 + idx18]
  mova                m12, [stp + %3 + idx19]
  SUM_SUB               4,  3, 9 ; stp2_16, stp2_23
  SUM_SUB               7,  0, 9 ; stp2_17, stp2_22
  SUM_SUB              11, 14, 9 ; stp2_18, stp2_21
  SUM_SUB              12,  5, 9 ; stp2_19, stp2_20
  mova [stp + %3 + idx16], m4
  mova [stp + %3 + idx17], m7
  mova [stp + %3 + idx18], m11
  mova [stp + %3 + idx19], m12

  mova                 m4, [stp + %4 + idx28]
  mova                 m7, [stp + %4 + idx29]
  mova                m11, [stp + %4 + idx30]
  mova                m12, [stp + %4 + idx31]
  SUM_SUB               4,  6, 9 ; stp2_28, stp2_27
  SUM_SUB               7, 13, 9 ; stp2_29, stp2_26
  SUM_SUB              11,  1, 9 ; stp2_30, stp2_25
  SUM_SUB              12,  2, 9 ; stp2_31, stp2_24
  mova [stp + %4 + idx28], m4
  mova [stp + %4 + idx29], m7
  mova [stp + %4 + idx30], m11
  mova [stp + %4 + idx31], m12

  ; BLOCK B STAGE 7 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%if 0 ; overflow occurs in SUM_SUB when using test streams
  mova                m10, [pw_11585x2]
  SUM_SUB               6,  5,  9
  pmulhrsw             m6, m10  ; stp1_27
  pmulhrsw             m5, m10  ; stp1_20
  SUM_SUB              13, 14,  9
  pmulhrsw            m13, m10  ; stp1_26
  pmulhrsw            m14, m10  ; stp1_21
  SUM_SUB               1,  0,  9
  pmulhrsw             m1, m10  ; stp1_25
  pmulhrsw             m0, m10  ; stp1_22
  SUM_SUB               2,  3,  9
  pmulhrsw             m2, m10  ; stp1_25
  pmulhrsw             m3, m10  ; stp1_22
%else
  BUTTERFLY_4X          6,     5,  11585, 11585,  m8,  9,  10 ; stp1_20, stp1_27
  SWAP  6, 5
  BUTTERFLY_4X         13,    14,  11585, 11585,  m8,  9,  10 ; stp1_21, stp1_26
  SWAP 13, 14
  BUTTERFLY_4X          1,     0,  11585, 11585,  m8,  9,  10 ; stp1_22, stp1_25
  SWAP  1, 0
  BUTTERFLY_4X          2,     3,  11585, 11585,  m8,  9,  10 ; stp1_23, stp1_24
  SWAP  2, 3
%endif
  mova [stp + %3 + idx20], m5
  mova [stp + %3 + idx21], m14
  mova [stp + %3 + idx22], m0
  mova [stp + %3 + idx23], m3
  mova [stp + %4 + idx24], m2
  mova [stp + %4 + idx25], m1
  mova [stp + %4 + idx26], m13
  mova [stp + %4 + idx27], m6

  ; BLOCK C STAGE 1 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ;
  ; BLOCK C STAGE 2 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  mova                 m0, [rsp + transposed_in + 16 *  2]
  mova                 m1, m0
  pmulhrsw             m0, [pw__1606x2] ; stp1_8
  pmulhrsw             m1, [pw_16305x2] ; stp2_15

  mova                 m6, [rsp + transposed_in + 16 *  6]
  mova                 m7, m6
  pmulhrsw             m7, [pw_m4756x2] ; stp2_11
  pmulhrsw             m6, [pw_15679x2] ; stp1_12

  mova                 m4, [rsp + transposed_in + 16 * 10]
  mova                 m5, m4
  pmulhrsw             m4, [pw__7723x2] ; stp1_10
  pmulhrsw             m5, [pw_14449x2] ; stp2_13

  mova                 m2, [rsp + transposed_in + 16 * 14]
  mova                 m3, m2
  pmulhrsw             m3, [pw_m10394x2] ; stp1_9
  pmulhrsw             m2, [pw_12665x2] ; stp2_14

  ; BLOCK C STAGE 3 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUM_SUB               0,  3, 9 ;  stp1_8, stp1_9
  SUM_SUB               7,  4, 9 ; stp1_11, stp1_10
  SUM_SUB               6,  5, 9 ; stp1_12, stp1_13
  SUM_SUB               1,  2, 9 ; stp1_15, stp1_14

  ; BLOCK C STAGE 4 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  BUTTERFLY_4X          2,     3,   6270, 15137,  m8,  9,  10 ;  stp1_9, stp1_14
  BUTTERFLY_4Xmm        5,     4,   6270, 15137,  m8,  9,  10 ; stp1_13, stp1_10

  ; BLOCK C STAGE 5 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUM_SUB               0,  7, 9 ;  stp1_8, stp1_11
  SUM_SUB               2,  4, 9 ;  stp1_9, stp1_10
  SUM_SUB               1,  6, 9 ;  stp1_15, stp1_12
  SUM_SUB               3,  5, 9 ;  stp1_14, stp1_13

  ; BLOCK C STAGE 6 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%if 0 ; overflow occurs in SUM_SUB when using test streams
  mova                m10, [pw_11585x2]
  SUM_SUB               5,    4,  9
  pmulhrsw             m5, m10  ; stp1_13
  pmulhrsw             m4, m10  ; stp1_10
  SUM_SUB               6,    7,  9
  pmulhrsw             m6, m10  ; stp1_12
  pmulhrsw             m7, m10  ; stp1_11
%else
  BUTTERFLY_4X       5,     4,  11585,  11585,  m8,  9,  10 ; stp1_10, stp1_13
  SWAP  5, 4
  BUTTERFLY_4X       6,     7,  11585,  11585,  m8,  9,  10 ; stp1_11, stp1_12
  SWAP  6, 7
%endif
  ; BLOCK C STAGE 7 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  mova [stp + %2 +  idx8], m0
  mova [stp + %2 +  idx9], m2
  mova [stp + %2 + idx10], m4
  mova [stp + %2 + idx11], m7
  mova [stp + %2 + idx12], m6
  mova [stp + %2 + idx13], m5
  mova [stp + %2 + idx14], m3
  mova [stp + %2 + idx15], m1

  ; BLOCK D STAGE 1 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ;
  ; BLOCK D STAGE 2 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ;
  ; BLOCK D STAGE 3 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  mova                m11, [rsp + transposed_in + 16 *  4]
  mova                m12, m11
  pmulhrsw            m11, [pw__3196x2] ; stp1_4
  pmulhrsw            m12, [pw_16069x2] ; stp1_7

  mova                m13, [rsp + transposed_in + 16 * 12]
  mova                m14, m13
  pmulhrsw            m13, [pw_13623x2] ; stp1_6
  pmulhrsw            m14, [pw_m9102x2] ; stp1_5

  ; BLOCK D STAGE 4 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  mova                 m0, [rsp + transposed_in + 16 *  0]
  mova                 m2, [rsp + transposed_in + 16 *  8]
  pmulhrsw             m0, [pw_11585x2]  ; stp1_1
  mova                 m3, m2
  pmulhrsw             m2, [pw__6270x2]  ; stp1_2
  pmulhrsw             m3, [pw_15137x2]  ; stp1_3

  SUM_SUB              11, 14, 9 ;  stp1_4, stp1_5
  SUM_SUB              12, 13, 9 ;  stp1_7, stp1_6

  ; BLOCK D STAGE 5 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%if 0 ; overflow occurs in SUM_SUB when using test streams
  mova                m10, [pw_11585x2]
  SUM_SUB              13,   14,  9
  pmulhrsw            m13, m10  ; stp1_6
  pmulhrsw            m14, m10  ; stp1_5
%else
  BUTTERFLY_4X         13,    14,  11585, 11585,  m8,  9,  10 ; stp1_5, stp1_6
  SWAP 13, 14
%endif
  mova                 m1, m0    ; stp1_0 = stp1_1
  SUM_SUB               0,  3, 9 ;  stp1_0, stp1_3
  SUM_SUB               1,  2, 9 ;  stp1_1, stp1_2

  ; BLOCK D STAGE 6 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUM_SUB               0, 12, 9 ;  stp1_0, stp1_7
  SUM_SUB               1, 13, 9 ;  stp1_1, stp1_6
  SUM_SUB               2, 14, 9 ;  stp1_2, stp1_5
  SUM_SUB               3, 11, 9 ;  stp1_3, stp1_4

  ; BLOCK D STAGE 7 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  mova                 m4, [stp + %2 + idx12]
  mova                 m5, [stp + %2 + idx13]
  mova                 m6, [stp + %2 + idx14]
  mova                 m7, [stp + %2 + idx15]
  SUM_SUB               0,  7, 9 ;  stp1_0, stp1_15
  SUM_SUB               1,  6, 9 ;  stp1_1, stp1_14
  SUM_SUB               2,  5, 9 ;  stp1_2, stp1_13
  SUM_SUB               3,  4, 9 ;  stp1_3, stp1_12

  ; 0-3, 28-31 final stage
  mova                m10, [stp + %4 + idx31]
  mova                m15, [stp + %4 + idx30]
  SUM_SUB               0, 10, 9 ;  stp1_0, stp1_31
  SUM_SUB               1, 15, 9 ;  stp1_1, stp1_30
  mova [stp + %1 +  idx0], m0
  mova [stp + %1 +  idx1], m1
  mova [stp + %4 + idx31], m10
  mova [stp + %4 + idx30], m15
  mova                 m0, [stp + %4 + idx29]
  mova                 m1, [stp + %4 + idx28]
  SUM_SUB               2,  0, 9 ;  stp1_2, stp1_29
  SUM_SUB               3,  1, 9 ;  stp1_3, stp1_28
  mova [stp + %1 +  idx2], m2
  mova [stp + %1 +  idx3], m3
  mova [stp + %4 + idx29], m0
  mova [stp + %4 + idx28], m1

  ; 12-15, 16-19 final stage
  mova                 m0, [stp + %3 + idx16]
  mova                 m1, [stp + %3 + idx17]
  mova                 m2, [stp + %3 + idx18]
  mova                 m3, [stp + %3 + idx19]
  SUM_SUB               7,  0, 9 ;  stp1_15, stp1_16
  SUM_SUB               6,  1, 9 ;  stp1_14, stp1_17
  SUM_SUB               5,  2, 9 ;  stp1_13, stp1_18
  SUM_SUB               4,  3, 9 ;  stp1_12, stp1_19
  mova [stp + %2 + idx12], m4
  mova [stp + %2 + idx13], m5
  mova [stp + %2 + idx14], m6
  mova [stp + %2 + idx15], m7
  mova [stp + %3 + idx16], m0
  mova [stp + %3 + idx17], m1
  mova [stp + %3 + idx18], m2
  mova [stp + %3 + idx19], m3

  mova                 m4, [stp + %2 +  idx8]
  mova                 m5, [stp + %2 +  idx9]
  mova                 m6, [stp + %2 + idx10]
  mova                 m7, [stp + %2 + idx11]
  SUM_SUB              11,  7, 9 ;  stp1_4, stp1_11
  SUM_SUB              14,  6, 9 ;  stp1_5, stp1_10
  SUM_SUB              13,  5, 9 ;  stp1_6, stp1_9
  SUM_SUB              12,  4, 9 ;  stp1_7, stp1_8

  ; 4-7, 24-27 final stage
  mova                 m3, [stp + %4 + idx24]
  mova                 m2, [stp + %4 + idx25]
  mova                 m1, [stp + %4 + idx26]
  mova                 m0, [stp + %4 + idx27]
  SUM_SUB              12,  3, 9 ;  stp1_7, stp1_24
  SUM_SUB              13,  2, 9 ;  stp1_6, stp1_25
  SUM_SUB              14,  1, 9 ;  stp1_5, stp1_26
  SUM_SUB              11,  0, 9 ;  stp1_4, stp1_27
  mova [stp + %4 + idx24], m3
  mova [stp + %4 + idx25], m2
  mova [stp + %4 + idx26], m1
  mova [stp + %4 + idx27], m0
  mova [stp + %1 +  idx4], m11
  mova [stp + %1 +  idx5], m14
  mova [stp + %1 +  idx6], m13
  mova [stp + %1 +  idx7], m12

  ; 8-11, 20-23 final stage
  mova                 m0, [stp + %3 + idx20]
  mova                 m1, [stp + %3 + idx21]
  mova                 m2, [stp + %3 + idx22]
  mova                 m3, [stp + %3 + idx23]
  SUM_SUB               7,  0, 9 ;  stp1_11, stp_20
  SUM_SUB               6,  1, 9 ;  stp1_10, stp_21
  SUM_SUB               5,  2, 9 ;   stp1_9, stp_22
  SUM_SUB               4,  3, 9 ;   stp1_8, stp_23
  mova [stp + %2 +  idx8], m4
  mova [stp + %2 +  idx9], m5
  mova [stp + %2 + idx10], m6
  mova [stp + %2 + idx11], m7
  mova [stp + %3 + idx20], m0
  mova [stp + %3 + idx21], m1
  mova [stp + %3 + idx22], m2
  mova [stp + %3 + idx23], m3
%endmacro

INIT_XMM ssse3
cglobal idct32x32_135_add, 3, 11, 16, i32x32_size, input, output, stride
  mova            m8, [pd_8192]
  mov             r6, 2
  lea            stp, [rsp + pass_one_start]

idct32x32_135:
  mov             r3, inputq
  lea             r4, [rsp + transposed_in]
  mov             r7, 2

idct32x32_135_transpose:
  LOAD_TRAN_LOW 0, r3,  0
  LOAD_TRAN_LOW 1, r3,  4
  LOAD_TRAN_LOW 2, r3,  8
  LOAD_TRAN_LOW 3, r3, 12
  LOAD_TRAN_LOW 4, r3, 16
  LOAD_TRAN_LOW 5, r3, 20
  LOAD_TRAN_LOW 6, r3, 24
  LOAD_TRAN_LOW 7, r3, 28

  TRANSPOSE8X8  0, 1, 2, 3, 4, 5, 6, 7, 9, 10

  mova [r4 +      0], m0
  mova [r4 + 16 * 1], m1
  mova [r4 + 16 * 2], m2
  mova [r4 + 16 * 3], m3
  mova [r4 + 16 * 4], m4
  mova [r4 + 16 * 5], m5
  mova [r4 + 16 * 6], m6
  mova [r4 + 16 * 7], m7

  INCREMENT_TRAN_LOW r3
  add             r4, 16 * 8
  dec             r7
  jne idct32x32_135_transpose

  IDCT32X32_135 16*0, 16*32, 16*64, 16*96
  lea            stp, [stp + 16 * 8]
  INCREMENT_ELEMENTS_TRAN_LOW inputq, 8*32
  dec             r6
  jnz idct32x32_135

  mov             r6, 4
  lea            stp, [rsp + pass_one_start]
  lea             r9, [rsp + pass_one_start]

idct32x32_135_2:
  lea             r4, [rsp + transposed_in]
  mov             r3, r9
  mov             r7, 2

idct32x32_135_transpose_2:
  mova            m0, [r3 +      0]
  mova            m1, [r3 + 16 * 1]
  mova            m2, [r3 + 16 * 2]
  mova            m3, [r3 + 16 * 3]
  mova            m4, [r3 + 16 * 4]
  mova            m5, [r3 + 16 * 5]
  mova            m6, [r3 + 16 * 6]
  mova            m7, [r3 + 16 * 7]

  TRANSPOSE8X8  0, 1, 2, 3, 4, 5, 6, 7, 9, 10

  mova [r4 +      0], m0
  mova [r4 + 16 * 1], m1
  mova [r4 + 16 * 2], m2
  mova [r4 + 16 * 3], m3
  mova [r4 + 16 * 4], m4
  mova [r4 + 16 * 5], m5
  mova [r4 + 16 * 6], m6
  mova [r4 + 16 * 7], m7

  add             r3, 16 * 8
  add             r4, 16 * 8
  dec             r7
  jne idct32x32_135_transpose_2

  IDCT32X32_135 16*0, 16*8, 16*16, 16*24

  lea            stp, [stp + 16 * 32]
  add             r9, 16 * 32
  dec             r6
  jnz idct32x32_135_2

  RECON_AND_STORE pass_two_start

  RET

%macro IDCT32X32_1024 4
  ; BLOCK A STAGE 1 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  mova                 m1, [rsp + transposed_in + 16 *  1]
  mova                m11, [rsp + transposed_in + 16 * 31]
  BUTTERFLY_4X          1,    11,    804, 16364,  m8,  9,  10 ; stp1_16, stp1_31

  mova                 m0, [rsp + transposed_in + 16 * 15]
  mova                 m2, [rsp + transposed_in + 16 * 17]
  BUTTERFLY_4X          2,     0,  12140, 11003,  m8,  9,  10 ; stp1_17, stp1_30

  mova                 m7, [rsp + transposed_in + 16 *  7]
  mova                m12, [rsp + transposed_in + 16 * 25]
  BUTTERFLY_4X         12,     7,  15426,  5520,  m8,  9,  10 ; stp1_19, stp1_28

  mova                 m3, [rsp + transposed_in + 16 *  9]
  mova                 m4, [rsp + transposed_in + 16 * 23]
  BUTTERFLY_4X          3,     4,   7005, 14811,  m8,  9,  10 ; stp1_18, stp1_29

  ; BLOCK A STAGE 2 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUM_SUB               1,  2, 9 ; stp2_16, stp2_17
  SUM_SUB              12,  3, 9 ; stp2_19, stp2_18
  SUM_SUB               7,  4, 9 ; stp2_28, stp2_29
  SUM_SUB              11,  0, 9 ; stp2_31, stp2_30

  ; BLOCK A STAGE 3 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  BUTTERFLY_4X          0,     2,   3196, 16069,  m8,  9,  10 ; stp1_17, stp1_30
  BUTTERFLY_4Xmm        4,     3,   3196, 16069,  m8,  9,  10 ; stp1_29, stp1_18

  ; BLOCK A STAGE 4 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUM_SUB               1, 12, 9 ; stp2_16, stp2_19
  SUM_SUB               0,  3, 9 ; stp2_17, stp2_18
  SUM_SUB              11,  7, 9 ; stp2_31, stp2_28
  SUM_SUB               2,  4, 9 ; stp2_30, stp2_29

  ; BLOCK A STAGE 5 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  BUTTERFLY_4X          4,     3,   6270, 15137,  m8,  9,  10 ; stp1_18, stp1_29
  BUTTERFLY_4X          7,    12,   6270, 15137,  m8,  9,  10 ; stp1_19, stp1_28

  mova [stp + %3 + idx16], m1
  mova [stp + %3 + idx17], m0
  mova [stp + %3 + idx18], m4
  mova [stp + %3 + idx19], m7
  mova [stp + %4 + idx28], m12
  mova [stp + %4 + idx29], m3
  mova [stp + %4 + idx30], m2
  mova [stp + %4 + idx31], m11

  ; BLOCK B STAGE 1 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  mova                 m5, [rsp + transposed_in + 16 *  5]
  mova                 m6, [rsp + transposed_in + 16 * 27]
  BUTTERFLY_4X          5,     6,   3981, 15893,  m8,  9,  10 ; stp1_20, stp1_27

  mova                m13, [rsp + transposed_in + 16 * 21]
  mova                m14, [rsp + transposed_in + 16 * 11]
  BUTTERFLY_4X         13,    14,  14053,  8423,  m8,  9,  10 ; stp1_21, stp1_26

  mova                 m0, [rsp + transposed_in + 16 * 13]
  mova                 m1, [rsp + transposed_in + 16 * 19]
  BUTTERFLY_4X          0,     1,   9760, 13160,  m8,  9,  10 ; stp1_22, stp1_25

  mova                 m2, [rsp + transposed_in + 16 *  3]
  mova                 m3, [rsp + transposed_in + 16 * 29]
  BUTTERFLY_4X          3,     2,  16207,  2404,  m8,  9,  10 ; stp1_23, stp1_24

  ; BLOCK B STAGE 2 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUM_SUB               5, 13, 9 ; stp2_20, stp2_21
  SUM_SUB               3,  0, 9 ; stp2_23, stp2_22
  SUM_SUB               2,  1, 9 ; stp2_24, stp2_25
  SUM_SUB               6, 14, 9 ; stp2_27, stp2_26

  ; BLOCK B STAGE 3 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  BUTTERFLY_4X         14,    13,  13623,  9102,  m8,  9,  10 ; stp1_21, stp1_26
  BUTTERFLY_4Xmm        1,     0,  13623,  9102,  m8,  9,  10 ; stp1_25, stp1_22

  ; BLOCK B STAGE 4 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUM_SUB               3,  5, 9 ; stp2_23, stp2_20
  SUM_SUB               0, 14, 9 ; stp2_22, stp2_21
  SUM_SUB               2,  6, 9 ; stp2_24, stp2_27
  SUM_SUB               1, 13, 9 ; stp2_25, stp2_26

  ; BLOCK B STAGE 5 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  BUTTERFLY_4Xmm        6,     5,   6270, 15137,  m8,  9,  10 ; stp1_27, stp1_20
  BUTTERFLY_4Xmm       13,    14,   6270, 15137,  m8,  9,  10 ; stp1_26, stp1_21

  ; BLOCK B STAGE 6 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  mova                 m4, [stp + %3 + idx16]
  mova                 m7, [stp + %3 + idx17]
  mova                m11, [stp + %3 + idx18]
  mova                m12, [stp + %3 + idx19]
  SUM_SUB               4,  3, 9 ; stp2_16, stp2_23
  SUM_SUB               7,  0, 9 ; stp2_17, stp2_22
  SUM_SUB              11, 14, 9 ; stp2_18, stp2_21
  SUM_SUB              12,  5, 9 ; stp2_19, stp2_20
  mova [stp + %3 + idx16], m4
  mova [stp + %3 + idx17], m7
  mova [stp + %3 + idx18], m11
  mova [stp + %3 + idx19], m12

  mova                 m4, [stp + %4 + idx28]
  mova                 m7, [stp + %4 + idx29]
  mova                m11, [stp + %4 + idx30]
  mova                m12, [stp + %4 + idx31]
  SUM_SUB               4,  6, 9 ; stp2_28, stp2_27
  SUM_SUB               7, 13, 9 ; stp2_29, stp2_26
  SUM_SUB              11,  1, 9 ; stp2_30, stp2_25
  SUM_SUB              12,  2, 9 ; stp2_31, stp2_24
  mova [stp + %4 + idx28], m4
  mova [stp + %4 + idx29], m7
  mova [stp + %4 + idx30], m11
  mova [stp + %4 + idx31], m12

  ; BLOCK B STAGE 7 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%if 0 ; overflow occurs in SUM_SUB when using test streams
  mova                m10, [pw_11585x2]
  SUM_SUB               6,  5,  9
  pmulhrsw             m6, m10  ; stp1_27
  pmulhrsw             m5, m10  ; stp1_20
  SUM_SUB              13, 14,  9
  pmulhrsw            m13, m10  ; stp1_26
  pmulhrsw            m14, m10  ; stp1_21
  SUM_SUB               1,  0,  9
  pmulhrsw             m1, m10  ; stp1_25
  pmulhrsw             m0, m10  ; stp1_22
  SUM_SUB               2,  3,  9
  pmulhrsw             m2, m10  ; stp1_25
  pmulhrsw             m3, m10  ; stp1_22
%else
  BUTTERFLY_4X          6,     5,  11585, 11585,  m8,  9,  10 ; stp1_20, stp1_27
  SWAP  6, 5
  BUTTERFLY_4X         13,    14,  11585, 11585,  m8,  9,  10 ; stp1_21, stp1_26
  SWAP 13, 14
  BUTTERFLY_4X          1,     0,  11585, 11585,  m8,  9,  10 ; stp1_22, stp1_25
  SWAP  1, 0
  BUTTERFLY_4X          2,     3,  11585, 11585,  m8,  9,  10 ; stp1_23, stp1_24
  SWAP  2, 3
%endif
  mova [stp + %3 + idx20], m5
  mova [stp + %3 + idx21], m14
  mova [stp + %3 + idx22], m0
  mova [stp + %3 + idx23], m3
  mova [stp + %4 + idx24], m2
  mova [stp + %4 + idx25], m1
  mova [stp + %4 + idx26], m13
  mova [stp + %4 + idx27], m6

  ; BLOCK C STAGE 1 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ;
  ; BLOCK C STAGE 2 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  mova                 m0, [rsp + transposed_in + 16 *  2]
  mova                 m1, [rsp + transposed_in + 16 * 30]
  BUTTERFLY_4X          0,     1,   1606, 16305,  m8,  9,  10 ; stp1_8, stp1_15

  mova                 m2, [rsp + transposed_in + 16 * 14]
  mova                 m3, [rsp + transposed_in + 16 * 18]
  BUTTERFLY_4X          3,     2,  12665, 10394,  m8,  9,  10 ; stp1_9, stp1_14

  mova                 m4, [rsp + transposed_in + 16 * 10]
  mova                 m5, [rsp + transposed_in + 16 * 22]
  BUTTERFLY_4X          4,     5,   7723, 14449,  m8,  9,  10 ; stp1_10, stp1_13

  mova                 m6, [rsp + transposed_in + 16 *  6]
  mova                 m7, [rsp + transposed_in + 16 * 26]
  BUTTERFLY_4X          7,     6,  15679,  4756,  m8,  9,  10 ; stp1_11, stp1_12

  ; BLOCK C STAGE 3 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUM_SUB               0,  3, 9 ;  stp1_8, stp1_9
  SUM_SUB               7,  4, 9 ; stp1_11, stp1_10
  SUM_SUB               6,  5, 9 ; stp1_12, stp1_13
  SUM_SUB               1,  2, 9 ; stp1_15, stp1_14

  ; BLOCK C STAGE 4 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  BUTTERFLY_4X          2,     3,   6270, 15137,  m8,  9,  10 ;  stp1_9, stp1_14
  BUTTERFLY_4Xmm        5,     4,   6270, 15137,  m8,  9,  10 ; stp1_13, stp1_10

  ; BLOCK C STAGE 5 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUM_SUB               0,  7, 9 ;  stp1_8, stp1_11
  SUM_SUB               2,  4, 9 ;  stp1_9, stp1_10
  SUM_SUB               1,  6, 9 ;  stp1_15, stp1_12
  SUM_SUB               3,  5, 9 ;  stp1_14, stp1_13

  ; BLOCK C STAGE 6 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%if 0 ; overflow occurs in SUM_SUB when using test streams
  mova                m10, [pw_11585x2]
  SUM_SUB               5,    4,  9
  pmulhrsw             m5, m10  ; stp1_13
  pmulhrsw             m4, m10  ; stp1_10
  SUM_SUB               6,    7,  9
  pmulhrsw             m6, m10  ; stp1_12
  pmulhrsw             m7, m10  ; stp1_11
%else
  BUTTERFLY_4X       5,     4,  11585,  11585,  m8,  9,  10 ; stp1_10, stp1_13
  SWAP  5, 4
  BUTTERFLY_4X       6,     7,  11585,  11585,  m8,  9,  10 ; stp1_11, stp1_12
  SWAP  6, 7
%endif
  ; BLOCK C STAGE 7 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  mova [stp + %2 +  idx8], m0
  mova [stp + %2 +  idx9], m2
  mova [stp + %2 + idx10], m4
  mova [stp + %2 + idx11], m7
  mova [stp + %2 + idx12], m6
  mova [stp + %2 + idx13], m5
  mova [stp + %2 + idx14], m3
  mova [stp + %2 + idx15], m1

  ; BLOCK D STAGE 1 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ;
  ; BLOCK D STAGE 2 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ;
  ; BLOCK D STAGE 3 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  mova                m11, [rsp + transposed_in + 16 *  4]
  mova                m12, [rsp + transposed_in + 16 * 28]
  BUTTERFLY_4X         11,    12,   3196, 16069,  m8,  9,  10 ; stp1_4, stp1_7

  mova                m13, [rsp + transposed_in + 16 * 12]
  mova                m14, [rsp + transposed_in + 16 * 20]
  BUTTERFLY_4X         14,    13,  13623,  9102,  m8,  9,  10 ; stp1_5, stp1_6

  ; BLOCK D STAGE 4 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  mova                 m0, [rsp + transposed_in + 16 *  0]
  mova                 m1, [rsp + transposed_in + 16 * 16]

%if 0 ; overflow occurs in SUM_SUB when using test streams
  mova                m10, [pw_11585x2]
  SUM_SUB               0,    1,  9
  pmulhrsw             m0, m10  ; stp1_1
  pmulhrsw             m1, m10  ; stp1_0
%else
  BUTTERFLY_4X          0,     1,  11585, 11585,  m8,  9,  10 ; stp1_1, stp1_0
  SWAP  0, 1
%endif
  mova                 m2, [rsp + transposed_in + 16 *  8]
  mova                 m3, [rsp + transposed_in + 16 * 24]
  BUTTERFLY_4X          2,     3,   6270, 15137,  m8,  9,  10 ;  stp1_2, stp1_3

  mova                m10, [pw_11585x2]
  SUM_SUB              11, 14, 9 ;  stp1_4, stp1_5
  SUM_SUB              12, 13, 9 ;  stp1_7, stp1_6

  ; BLOCK D STAGE 5 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%if 0 ; overflow occurs in SUM_SUB when using test streams
  SUM_SUB              13,   14,  9
  pmulhrsw            m13, m10  ; stp1_6
  pmulhrsw            m14, m10  ; stp1_5
%else
  BUTTERFLY_4X         13,    14,  11585, 11585,  m8,  9,  10 ; stp1_5, stp1_6
  SWAP 13, 14
%endif
  SUM_SUB               0,  3, 9 ;  stp1_0, stp1_3
  SUM_SUB               1,  2, 9 ;  stp1_1, stp1_2

  ; BLOCK D STAGE 6 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUM_SUB               0, 12, 9 ;  stp1_0, stp1_7
  SUM_SUB               1, 13, 9 ;  stp1_1, stp1_6
  SUM_SUB               2, 14, 9 ;  stp1_2, stp1_5
  SUM_SUB               3, 11, 9 ;  stp1_3, stp1_4

  ; BLOCK D STAGE 7 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  mova                 m4, [stp + %2 + idx12]
  mova                 m5, [stp + %2 + idx13]
  mova                 m6, [stp + %2 + idx14]
  mova                 m7, [stp + %2 + idx15]
  SUM_SUB               0,  7, 9 ;  stp1_0, stp1_15
  SUM_SUB               1,  6, 9 ;  stp1_1, stp1_14
  SUM_SUB               2,  5, 9 ;  stp1_2, stp1_13
  SUM_SUB               3,  4, 9 ;  stp1_3, stp1_12

  ; 0-3, 28-31 final stage
  mova                m10, [stp + %4 + idx31]
  mova                m15, [stp + %4 + idx30]
  SUM_SUB               0, 10, 9 ;  stp1_0, stp1_31
  SUM_SUB               1, 15, 9 ;  stp1_1, stp1_30
  mova [stp + %1 +  idx0], m0
  mova [stp + %1 +  idx1], m1
  mova [stp + %4 + idx31], m10
  mova [stp + %4 + idx30], m15
  mova                 m0, [stp + %4 + idx29]
  mova                 m1, [stp + %4 + idx28]
  SUM_SUB               2,  0, 9 ;  stp1_2, stp1_29
  SUM_SUB               3,  1, 9 ;  stp1_3, stp1_28
  mova [stp + %1 +  idx2], m2
  mova [stp + %1 +  idx3], m3
  mova [stp + %4 + idx29], m0
  mova [stp + %4 + idx28], m1

  ; 12-15, 16-19 final stage
  mova                 m0, [stp + %3 + idx16]
  mova                 m1, [stp + %3 + idx17]
  mova                 m2, [stp + %3 + idx18]
  mova                 m3, [stp + %3 + idx19]
  SUM_SUB               7,  0, 9 ;  stp1_15, stp1_16
  SUM_SUB               6,  1, 9 ;  stp1_14, stp1_17
  SUM_SUB               5,  2, 9 ;  stp1_13, stp1_18
  SUM_SUB               4,  3, 9 ;  stp1_12, stp1_19
  mova [stp + %2 + idx12], m4
  mova [stp + %2 + idx13], m5
  mova [stp + %2 + idx14], m6
  mova [stp + %2 + idx15], m7
  mova [stp + %3 + idx16], m0
  mova [stp + %3 + idx17], m1
  mova [stp + %3 + idx18], m2
  mova [stp + %3 + idx19], m3

  mova                 m4, [stp + %2 +  idx8]
  mova                 m5, [stp + %2 +  idx9]
  mova                 m6, [stp + %2 + idx10]
  mova                 m7, [stp + %2 + idx11]
  SUM_SUB              11,  7, 9 ;  stp1_4, stp1_11
  SUM_SUB              14,  6, 9 ;  stp1_5, stp1_10
  SUM_SUB              13,  5, 9 ;  stp1_6, stp1_9
  SUM_SUB              12,  4, 9 ;  stp1_7, stp1_8

  ; 4-7, 24-27 final stage
  mova                 m3, [stp + %4 + idx24]
  mova                 m2, [stp + %4 + idx25]
  mova                 m1, [stp + %4 + idx26]
  mova                 m0, [stp + %4 + idx27]
  SUM_SUB              12,  3, 9 ;  stp1_7, stp1_24
  SUM_SUB              13,  2, 9 ;  stp1_6, stp1_25
  SUM_SUB              14,  1, 9 ;  stp1_5, stp1_26
  SUM_SUB              11,  0, 9 ;  stp1_4, stp1_27
  mova [stp + %4 + idx24], m3
  mova [stp + %4 + idx25], m2
  mova [stp + %4 + idx26], m1
  mova [stp + %4 + idx27], m0
  mova [stp + %1 +  idx4], m11
  mova [stp + %1 +  idx5], m14
  mova [stp + %1 +  idx6], m13
  mova [stp + %1 +  idx7], m12

  ; 8-11, 20-23 final stage
  mova                 m0, [stp + %3 + idx20]
  mova                 m1, [stp + %3 + idx21]
  mova                 m2, [stp + %3 + idx22]
  mova                 m3, [stp + %3 + idx23]
  SUM_SUB               7,  0, 9 ;  stp1_11, stp_20
  SUM_SUB               6,  1, 9 ;  stp1_10, stp_21
  SUM_SUB               5,  2, 9 ;   stp1_9, stp_22
  SUM_SUB               4,  3, 9 ;   stp1_8, stp_23
  mova [stp + %2 +  idx8], m4
  mova [stp + %2 +  idx9], m5
  mova [stp + %2 + idx10], m6
  mova [stp + %2 + idx11], m7
  mova [stp + %3 + idx20], m0
  mova [stp + %3 + idx21], m1
  mova [stp + %3 + idx22], m2
  mova [stp + %3 + idx23], m3
%endmacro

INIT_XMM ssse3
cglobal idct32x32_1024_add, 3, 11, 16, i32x32_size, input, output, stride
  mova            m8, [pd_8192]
  mov             r6, 4
  lea            stp, [rsp + pass_one_start]

idct32x32_1024:
  mov             r3, inputq
  lea             r4, [rsp + transposed_in]
  mov             r7, 4

idct32x32_1024_transpose:
  LOAD_TRAN_LOW 0, r3,  0
  LOAD_TRAN_LOW 1, r3,  4
  LOAD_TRAN_LOW 2, r3,  8
  LOAD_TRAN_LOW 3, r3, 12
  LOAD_TRAN_LOW 4, r3, 16
  LOAD_TRAN_LOW 5, r3, 20
  LOAD_TRAN_LOW 6, r3, 24
  LOAD_TRAN_LOW 7, r3, 28

  TRANSPOSE8X8  0, 1, 2, 3, 4, 5, 6, 7, 9, 10

  mova [r4 +      0], m0
  mova [r4 + 16 * 1], m1
  mova [r4 + 16 * 2], m2
  mova [r4 + 16 * 3], m3
  mova [r4 + 16 * 4], m4
  mova [r4 + 16 * 5], m5
  mova [r4 + 16 * 6], m6
  mova [r4 + 16 * 7], m7
  INCREMENT_TRAN_LOW r3
  add             r4, 16 * 8
  dec             r7
  jne idct32x32_1024_transpose

  IDCT32X32_1024 16*0, 16*32, 16*64, 16*96

  lea            stp, [stp + 16 * 8]
  INCREMENT_ELEMENTS_TRAN_LOW inputq, 8*32
  dec             r6
  jnz idct32x32_1024

  mov             r6, 4
  lea            stp, [rsp + pass_one_start]
  lea             r9, [rsp + pass_one_start]

idct32x32_1024_2:
  lea             r4, [rsp + transposed_in]
  mov             r3, r9
  mov             r7, 4

idct32x32_1024_transpose_2:
  mova            m0, [r3 +      0]
  mova            m1, [r3 + 16 * 1]
  mova            m2, [r3 + 16 * 2]
  mova            m3, [r3 + 16 * 3]
  mova            m4, [r3 + 16 * 4]
  mova            m5, [r3 + 16 * 5]
  mova            m6, [r3 + 16 * 6]
  mova            m7, [r3 + 16 * 7]

  TRANSPOSE8X8  0, 1, 2, 3, 4, 5, 6, 7, 9, 10

  mova [r4 +      0], m0
  mova [r4 + 16 * 1], m1
  mova [r4 + 16 * 2], m2
  mova [r4 + 16 * 3], m3
  mova [r4 + 16 * 4], m4
  mova [r4 + 16 * 5], m5
  mova [r4 + 16 * 6], m6
  mova [r4 + 16 * 7], m7

  add             r3, 16 * 8
  add             r4, 16 * 8
  dec             r7
  jne idct32x32_1024_transpose_2

  IDCT32X32_1024 16*0, 16*8, 16*16, 16*24

  lea            stp, [stp + 16 * 32]
  add             r9, 16 * 32
  dec             r6
  jnz idct32x32_1024_2

  RECON_AND_STORE pass_two_start

  RET
%endif
