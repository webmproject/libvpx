;
;  Copyright (c) 2010 The Webm project authors. All Rights Reserved.
;
;  Use of this source code is governed by a BSD-style license
;  that can be found in the LICENSE file in the root of the source
;  tree. An additional intellectual property rights grant can be found
;  in the file PATENTS.  All contributing project authors may
;  be found in the AUTHORS file in the root of the source tree.
;


    EXPORT  |idct_dequant_full_2x_neon|
    ARM
    REQUIRE8
    PRESERVE8

    AREA ||.text||, CODE, READONLY, ALIGN=2
;void idct_dequant_full_2x_neon(short *q, short *dq, unsigned char *pre,
;                               unsigned char *dst, int pitch, int stride);
; r0    *q,
; r1    *dq,
; r2    *pre
; r3    *dst
; sp    pitch
; sp+4  stride
|idct_dequant_full_2x_neon| PROC
    vld1.16         {q3, q4}, [r0]          ; lo input
    vld1.16         {q5, q6}, [r1]          ; use the same dq for both
    ldr             r1, [sp]                ; pitch
    add             r0, r0, #32
    vld1.16         {q10, q11}, [r0]        ; hi input
    add             r12, r2, #4
    vld1.32         {d14[0]}, [r2], r1      ; lo pred
    vld1.32         {d14[1]}, [r2], r1
    vld1.32         {d15[0]}, [r2], r1
    vld1.32         {d15[1]}, [r2]
    vld1.32         {d28[0]}, [r12], r1     ; hi pred
    vld1.32         {d28[1]}, [r12], r1
    vld1.32         {d29[0]}, [r12], r1
    vld1.32         {d29[1]}, [r12]

    ldr             r2, _CONSTANTS_

    vmul.i16        q1, q3, q5              ; lo input * dq
    vmul.i16        q2, q4, q6
    vmul.i16        q8, q10, q5             ; hi input * dq
    vmul.i16        q9, q11, q6

    ldr             r1, [sp, #4]            ; stride

    vld1.16         {d0}, [r2]
    vswp            d3, d4                  ; lo q2(vp[4] vp[12])
    vswp            d17, d18                ; hi q2(vp[4] vp[12])

    vqdmulh.s16     q3, q2, d0[2]           ; lo * constants
    vqdmulh.s16     q4, q2, d0[0]
    vqdmulh.s16     q10, q9, d0[2]          ; hi * constants
    vqdmulh.s16     q11, q9, d0[0]

    vqadd.s16       d12, d2, d3             ; lo a1
    vqsub.s16       d13, d2, d3             ; lo b1
    vqadd.s16       d26, d16, d17           ; hi a1
    vqsub.s16       d27, d16, d17           ; hi b1

    vshr.s16        q3, q3, #1              ; lo
    vshr.s16        q4, q4, #1
    vshr.s16        q10, q10, #1            ; hi
    vshr.s16        q11, q11, #1

    vqadd.s16       q3, q3, q2              ; lo
    vqadd.s16       q4, q4, q2
    vqadd.s16       q10, q10, q9            ; hi
    vqadd.s16       q11, q11, q9

    vqsub.s16       d10, d6, d9             ; lo c1
    vqadd.s16       d11, d7, d8             ; lo d1
    vqsub.s16       d24, d20, d23           ; hi c1
    vqadd.s16       d25, d21, d22           ; hi d1

    vqadd.s16       d2, d12, d11            ; lo
    vqadd.s16       d3, d13, d10
    vqsub.s16       d4, d13, d10
    vqsub.s16       d5, d12, d11
    vqadd.s16       d16, d26, d25           ; hi
    vqadd.s16       d17, d27, d24
    vqsub.s16       d18, d27, d24
    vqsub.s16       d19, d26, d25

    vtrn.32         d2, d4                  ; lo
    vtrn.32         d3, d5
    vtrn.16         d2, d3
    vtrn.16         d4, d5
    vtrn.32         d16, d18                ; hi
    vtrn.32         d17, d19
    vtrn.16         d16, d17
    vtrn.16         d18, d19

    vswp            d3, d4                  ; lo
    vqdmulh.s16     q3, q2, d0[2]
    vqdmulh.s16     q4, q2, d0[0]
    vswp            d17, d18                ; hi
    vqdmulh.s16     q10, q9, d0[2]
    vqdmulh.s16     q11, q9, d0[0]

    vqadd.s16       d12, d2, d3             ; lo a1
    vqsub.s16       d13, d2, d3             ; lo b1
    vqadd.s16       d26, d16, d17           ; hi a1
    vqsub.s16       d27, d16, d17           ; hi b1

    vshr.s16        q3, q3, #1              ; lo
    vshr.s16        q4, q4, #1
    vshr.s16        q10, q10, #1            ; hi
    vshr.s16        q11, q11, #1

    vqadd.s16       q3, q3, q2              ; lo
    vqadd.s16       q4, q4, q2
    vqadd.s16       q10, q10, q9            ; hi
    vqadd.s16       q11, q11, q9

    vqsub.s16       d10, d6, d9             ; lo c1
    vqadd.s16       d11, d7, d8             ; lo d1
    vqsub.s16       d24, d20, d23           ; hi c1
    vqadd.s16       d25, d21, d22           ; hi d1

    vqadd.s16       d2, d12, d11            ; lo
    vqadd.s16       d3, d13, d10
    vqsub.s16       d4, d13, d10
    vqsub.s16       d5, d12, d11
    vqadd.s16       d16, d26, d25           ; hi
    vqadd.s16       d17, d27, d24
    vqsub.s16       d18, d27, d24
    vqsub.s16       d19, d26, d25

    vrshr.s16       q1, q1, #3              ; lo
    vrshr.s16       q2, q2, #3
    vrshr.s16       q8, q8, #3              ; hi
    vrshr.s16       q9, q9, #3

    vtrn.32         d2, d4                  ; lo
    vtrn.32         d3, d5
    vtrn.16         d2, d3
    vtrn.16         d4, d5
    vtrn.32         d16, d18                ; hi
    vtrn.32         d17, d19
    vtrn.16         d16, d17
    vtrn.16         d18, d19

    vaddw.u8        q1, q1, d14             ; lo
    vaddw.u8        q2, q2, d15
    vaddw.u8        q8, q8, d28             ; hi
    vaddw.u8        q9, q9, d29

    vmov.i16        q14, #0
    vmov            q15, q14
    vst1.16         {q14, q15}, [r0]        ; write over high input
    sub             r0, r0, #32
    vst1.16         {q14, q15}, [r0]        ; write over low input

    vqmovun.s16     d0, q1                  ; lo
    vqmovun.s16     d1, q2
    vqmovun.s16     d2, q8                  ; hi
    vqmovun.s16     d3, q9

    add             r2, r3, #4              ; hi
    vst1.32         {d0[0]}, [r3], r1       ; lo
    vst1.32         {d0[1]}, [r3], r1
    vst1.32         {d1[0]}, [r3], r1
    vst1.32         {d1[1]}, [r3]
    vst1.32         {d2[0]}, [r2], r1       ; hi
    vst1.32         {d2[1]}, [r2], r1
    vst1.32         {d3[0]}, [r2], r1
    vst1.32         {d3[1]}, [r2]

    bx             lr

    ENDP           ; |idct_dequant_full_2x_neon|

; Constant Pool
_CONSTANTS_       DCD cospi8sqrt2minus1
cospi8sqrt2minus1 DCD 0x4e7b4e7b
sinpi8sqrt2       DCD 0x8a8c8a8c

    END
