;
;  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
;
;  Use of this source code is governed by a BSD-style license
;  that can be found in the LICENSE file in the root of the source
;  tree. An additional intellectual property rights grant can be found
;  in the file PATENTS.  All contributing project authors may
;  be found in the AUTHORS file in the root of the source tree.
;


    EXPORT  |vp8_fast_quantize_b_neon_func|

    ARM
    REQUIRE8
    PRESERVE8

    AREA ||.text||, CODE, READONLY, ALIGN=2

; r0        short *coeff_ptr
; r1        short *zbin_ptr
; r2        short *qcoeff_ptr
; r3        short *dqcoeff_ptr
; stack     short *dequant_ptr
; stack     short *scan_mask
; stack     short *round_ptr
; stack     short *quant_ptr

; return    int * eob
|vp8_fast_quantize_b_neon_func| PROC
    vld1.16         {q0, q1}, [r0]              ;load z
    vld1.16         {q10, q11}, [r1]            ;load zbin

    vabs.s16        q4, q0                      ;calculate x = abs(z)
    vabs.s16        q5, q1

    vcge.s16        q10, q4, q10                ;x>=zbin
    vcge.s16        q11, q5, q11

    ;if x<zbin (q10 & q11 are all 0), go to zero_output
    vorr.s16        q6, q10, q11
    vorr.s16        d12, d12, d13
    vmov            r0, r1, d12
    orr             r0, r0, r1
    cmp             r0, #0
    beq             zero_output

    ldr             r0, [sp, #8]                ;load round_ptr
    ldr             r12, [sp, #12]              ;load quant_ptr

    ;right shift 15 to get sign, all 0 if it is positive, all 1 if it is negative
    vshr.s16        q2, q0, #15                 ; sz
    vshr.s16        q3, q1, #15

    vld1.s16        {q6, q7}, [r0]              ;load round_ptr [0-15]
    vld1.s16        {q8, q9}, [r12]             ;load quant_ptr [0-15]

    vadd.s16        q4, q6                      ;x + Round
    vadd.s16        q5, q7

    ldr             r0, [sp, #4]                ;load rvsplus1_scan_order ptr

    vqdmulh.s16     q4, q8                      ;y = ((Round + abs(z)) * Quant) >> 16
    vqdmulh.s16     q5, q9

    vld1.16         {q0, q1}, [r0]              ;load rvsplus1_scan_order
    vceq.s16        q8, q8                      ;set q8 to all 1

    vshr.s16        q4, #1                      ;right shift 1 after vqdmulh
    vshr.s16        q5, #1

    ;modify data to have its original sign
    veor.s16        q4, q2                      ; y^sz
    veor.s16        q5, q3

    ldr             r12, [sp]                   ;load dequant_ptr

    vsub.s16        q4, q2                      ; x1 = (y^sz) - sz = (y^sz) - (-1) (two's complement)
    vsub.s16        q5, q3

    vand.s16        q4, q10                     ;mask off x1 elements
    vand.s16        q5, q11

    vld1.s16        {q6, q7}, [r12]             ;load dequant_ptr[i]

    vtst.16         q14, q4, q8                 ;now find eob
    vtst.16         q15, q5, q8                 ;non-zero element is set to all 1 in q4, q5

    vst1.s16        {q4, q5}, [r2]              ;store: qcoeff = x1

    vand            q0, q0, q14                 ;get all valid number from rvsplus1_scan_order array
    vand            q1, q1, q15

    vmax.u16        q0, q0, q1                  ;find maximum value in q0, q1
    vmax.u16        d0, d0, d1
    vmovl.u16       q0, d0

    vmul.s16        q6, q4                      ;x * Dequant
    vmul.s16        q7, q5

    vmax.u32        d0, d0, d1
    vpmax.u32       d0, d0, d0

    vst1.s16        {q6, q7}, [r3]              ;store dqcoeff = x * Dequant

    vmov.32         r0, d0[0]
    bx              lr

zero_output
    vst1.s16        {q10, q11}, [r2]        ; qcoeff = 0
    vst1.s16        {q10, q11}, [r3]        ; dqcoeff = 0
    mov             r0, #0

    bx              lr

    ENDP

    END
