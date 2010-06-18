;
;  Copyright (c) 2010 The VP8 project authors. All Rights Reserved.
;
;  Use of this source code is governed by a BSD-style license
;  that can be found in the LICENSE file in the root of the source
;  tree. An additional intellectual property rights grant can be found
;  in the file PATENTS.  All contributing project authors may
;  be found in the AUTHORS file in the root of the source tree.
;


    EXPORT  |vp8_dequant_dc_idct_v6|
    ; ARM
    ; REQUIRE8
    ; PRESERVE8

    AREA    |.text|, CODE, READONLY  ; name this block of code
;void vp8_dequant_dc_idct_v6(short *input, short *dq, short *output, int pitch,int Dc)
|vp8_dequant_dc_idct_v6| PROC
    stmdb   sp!, {r4-r11, lr}

    ldr     r6, [sp, #36]           ;load Dc

    ldr     r4, [r0]                ;input
    ldr     r5, [r1], #4            ;dq

    sub     sp, sp, #4
    str     r0, [sp]

    smultt  r7, r4, r5

    ldr     r4, [r0, #4]            ;input
    ldr     r5, [r1], #4            ;dq

    strh    r6, [r0], #2
    strh    r7, [r0], #2

    smulbb  r6, r4, r5
    smultt  r7, r4, r5

    ldr     r4, [r0, #4]            ;input
    ldr     r5, [r1], #4            ;dq

    strh    r6, [r0], #2
    strh    r7, [r0], #2

    mov     r12, #3

dequant_dc_idct_loop
    smulbb  r6, r4, r5
    smultt  r7, r4, r5

    ldr     r4, [r0, #4]            ;input
    ldr     r5, [r1], #4            ;dq

    strh    r6, [r0], #2
    strh    r7, [r0], #2

    smulbb  r6, r4, r5
    smultt  r7, r4, r5

    subs    r12, r12, #1

    ldrne   r4, [r0, #4]
    ldrne   r5, [r1], #4

    strh    r6, [r0], #2
    strh    r7, [r0], #2

    bne     dequant_dc_idct_loop

    sub     r0, r0, #32
    mov     r1, r2
    mov     r2, r3

; short_idct4x4llm_v6_dual

    mov r3, #0x00004E00 ;                   cos
    orr r3, r3, #0x0000007B ; cospi8sqrt2minus1
    mov r4, #0x00008A00 ;                       sin
    orr r4, r4, #0x0000008C ; sinpi8sqrt2
    mov r5, #0x2    ; i=2                           i
loop1_dual_11
    ldr r6, [r0, #(4*2)]    ; i5 | i4                               5|4
    ldr r12, [r0, #(12*2)]  ; i13 | i12                                                     13|12
    ldr r14, [r0, #(8*2)]   ; i9 | i8                                                               9|8

    smulwt  r9, r3, r6  ; (ip[5] * cospi8sqrt2minus1) >> 16                                         5c
    smulwb  r7, r3, r6  ; (ip[4] * cospi8sqrt2minus1) >> 16                                 4c
    smulwt  r10, r4, r6 ; (ip[5] * sinpi8sqrt2) >> 16                                               5s
    smulwb  r8, r4, r6  ; (ip[4] * sinpi8sqrt2) >> 16                                       4s
    pkhbt   r7, r7, r9, lsl #16 ; 5c | 4c
    smulwt  r11, r3, r12    ; (ip[13] * cospi8sqrt2minus1) >> 16                                                    13c
    pkhbt   r8, r8, r10, lsl #16    ; 5s | 4s
    uadd16  r6, r6, r7  ; 5c+5 | 4c+4
    smulwt  r7, r4, r12 ; (ip[13] * sinpi8sqrt2) >> 16                                  13s
    smulwb  r9, r3, r12 ; (ip[12] * cospi8sqrt2minus1) >> 16                                            12c
    smulwb  r10, r4, r12    ; (ip[12] * sinpi8sqrt2) >> 16                                              12s
    subs    r5, r5, #0x1    ; i--                           --
    pkhbt   r9, r9, r11, lsl #16    ; 13c | 12c
    ldr r11, [r0], #0x4 ; i1 | i0       ++                                          1|0
    pkhbt   r10, r10, r7, lsl #16   ; 13s | 12s
    uadd16  r7, r12, r9 ; 13c+13 | 12c+12
    usub16  r7, r8, r7  ; c                                 c
    uadd16  r6, r6, r10 ; d                             d
    uadd16  r10, r11, r14   ; a                                             a
    usub16  r8, r11, r14    ; b                                     b
    uadd16  r9, r10, r6 ; a+d                                           a+d
    usub16  r10, r10, r6    ; a-d                                               a-d
    uadd16  r6, r8, r7  ; b+c                               b+c
    usub16  r7, r8, r7  ; b-c                                   b-c
    str r6, [r1, r2]    ; o5 | o4
    add r6, r2, r2  ; pitch * 2                             p2
    str r7, [r1, r6]    ; o9 | o8
    add r6,  r6, r2 ; pitch * 3                             p3
    str r10, [r1, r6]   ; o13 | o12
    str r9, [r1], #0x4  ; o1 | o0           ++
    bne loop1_dual_11   ;
    mov r5, #0x2    ; i=2                           i
    sub r0, r1, #8  ; reset input/output        i/o
loop2_dual_22
    ldr r6, [r0, r2]    ; i5 | i4                               5|4
    ldr r1, [r0]    ; i1 | i0           1|0
    ldr r12, [r0, #0x4] ; i3 | i2                                                       3|2
    add r14, r2, #0x4   ; pitch + 2                                                             p+2
    ldr r14, [r0, r14]  ; i7 | i6                                                               7|6
    smulwt  r9, r3, r6  ; (ip[5] * cospi8sqrt2minus1) >> 16                                         5c
    smulwt  r7, r3, r1  ; (ip[1] * cospi8sqrt2minus1) >> 16                                 1c
    smulwt  r10, r4, r6 ; (ip[5] * sinpi8sqrt2) >> 16                                               5s
    smulwt  r8, r4, r1  ; (ip[1] * sinpi8sqrt2) >> 16                                       1s
    pkhbt   r11, r6, r1, lsl #16    ; i0 | i4                                                   0|4
    pkhbt   r7, r9, r7, lsl #16 ; 1c | 5c
    pkhbt   r8, r10, r8, lsl #16    ; 1s | 5s = temp1 ©                                     tc1
    pkhtb   r1, r1, r6, asr #16 ; i1 | i5           1|5
    uadd16  r1, r7, r1  ; 1c+1 | 5c+5 = temp2 (d)           td2
    pkhbt   r9, r14, r12, lsl #16   ; i2 | i6                                           2|6
    uadd16  r10, r11, r9    ; a                                             a
    usub16  r9, r11, r9 ; b                                         b
    pkhtb   r6, r12, r14, asr #16   ; i3 | i7                               3|7
    subs    r5, r5, #0x1    ; i--                           --
    smulwt  r7, r3, r6  ; (ip[3] * cospi8sqrt2minus1) >> 16                                 3c
    smulwt  r11, r4, r6 ; (ip[3] * sinpi8sqrt2) >> 16                                                   3s
    smulwb  r12, r3, r6 ; (ip[7] * cospi8sqrt2minus1) >> 16                                                     7c
    smulwb  r14, r4, r6 ; (ip[7] * sinpi8sqrt2) >> 16                                                               7s

    pkhbt   r7, r12, r7, lsl #16    ; 3c | 7c
    pkhbt   r11, r14, r11, lsl #16  ; 3s | 7s = temp1 (d)                                                   td1
    uadd16  r6, r7, r6  ; 3c+3 | 7c+7 = temp2  (c)                              tc2
    usub16  r12, r8, r6 ; c (o1 | o5)                                                       c
    uadd16  r6, r11, r1 ; d (o3 | o7)                               d
    uadd16  r7, r10, r6 ; a+d                                   a+d
    mov r8, #0x4    ; set up 4's                                        4
    orr r8, r8, #0x40000    ;                                       4|4
    usub16  r6, r10, r6 ; a-d                               a-d
    uadd16  r6, r6, r8  ; a-d+4                             3|7
    uadd16  r7, r7, r8  ; a+d+4                                 0|4
    uadd16  r10, r9, r12    ; b+c                                               b+c
    usub16  r1, r9, r12 ; b-c           b-c
    uadd16  r10, r10, r8    ; b+c+4                                             1|5
    uadd16  r1, r1, r8  ; b-c+4         2|6
    mov r8, r10, asr #19    ; o1 >> 3
    strh    r8, [r0, #2]    ; o1
    mov r8, r1, asr #19 ; o2 >> 3
    strh    r8, [r0, #4]    ; o2
    mov r8, r6, asr #19 ; o3 >> 3
    strh    r8, [r0, #6]    ; o3
    mov r8, r7, asr #19 ; o0 >> 3
    strh    r8, [r0], r2    ; o0        +p
    sxth    r10, r10    ;
    mov r8, r10, asr #3 ; o5 >> 3
    strh    r8, [r0, #2]    ; o5
    sxth    r1, r1  ;
    mov r8, r1, asr #3  ; o6 >> 3
    strh    r8, [r0, #4]    ; o6
    sxth    r6, r6  ;
    mov r8, r6, asr #3  ; o7 >> 3
    strh    r8, [r0, #6]    ; o7
    sxth    r7, r7  ;
    mov r8, r7, asr #3  ; o4 >> 3
    strh    r8, [r0], r2    ; o4        +p
;;;;;   subs    r5, r5, #0x1    ; i--                           --
    bne loop2_dual_22   ;


;vpx_memset
    ldr     r0, [sp]
    add     sp, sp, #4

    mov     r12, #0
    str     r12, [r0]
    str     r12, [r0, #4]
    str     r12, [r0, #8]
    str     r12, [r0, #12]
    str     r12, [r0, #16]
    str     r12, [r0, #20]
    str     r12, [r0, #24]
    str     r12, [r0, #28]

    ldmia   sp!, {r4 - r11, pc} ; replace vars, return                      restore

    ENDP    ;|vp8_dequant_dc_idct_v68|

    END
