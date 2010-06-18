;
;  Copyright (c) 2010 The VP8 project authors. All Rights Reserved.
;
;  Use of this source code is governed by a BSD-style license
;  that can be found in the LICENSE file in the root of the source
;  tree. An additional intellectual property rights grant can be found
;  in the file PATENTS.  All contributing project authors may
;  be found in the AUTHORS file in the root of the source tree.
;


    EXPORT  |vp8_dequant_dc_idct_neon|
    ARM
    REQUIRE8
    PRESERVE8

    AREA ||.text||, CODE, READONLY, ALIGN=2
;void vp8_dequant_dc_idct_c(short *input, short *dq, short *output, int pitch, int Dc);
; r0    short *input,
; r1    short *dq,
; r2    short *output,
; r3    int pitch,
; (stack)   int Dc
|vp8_dequant_dc_idct_neon| PROC
    vld1.16         {q3, q4}, [r0]
    vld1.16         {q5, q6}, [r1]

    ldr             r1, [sp]                ;load Dc from stack

    ldr             r12, _dcidct_coeff_

    vmul.i16        q1, q3, q5              ;input for short_idct4x4llm_neon
    vmul.i16        q2, q4, q6

    vmov.16         d2[0], r1

;|short_idct4x4llm_neon| PROC
    vld1.16         {d0}, [r12]
    vswp            d3, d4                  ;q2(vp[4] vp[12])

    vqdmulh.s16     q3, q2, d0[2]
    vqdmulh.s16     q4, q2, d0[0]

    vqadd.s16       d12, d2, d3             ;a1
    vqsub.s16       d13, d2, d3             ;b1

    vshr.s16        q3, q3, #1
    vshr.s16        q4, q4, #1

    vqadd.s16       q3, q3, q2              ;modify since sinpi8sqrt2 > 65536/2 (negtive number)
    vqadd.s16       q4, q4, q2

    ;d6 - c1:temp1
    ;d7 - d1:temp2
    ;d8 - d1:temp1
    ;d9 - c1:temp2

    vqsub.s16       d10, d6, d9             ;c1
    vqadd.s16       d11, d7, d8             ;d1

    vqadd.s16       d2, d12, d11
    vqadd.s16       d3, d13, d10
    vqsub.s16       d4, d13, d10
    vqsub.s16       d5, d12, d11

    vtrn.32         d2, d4
    vtrn.32         d3, d5
    vtrn.16         d2, d3
    vtrn.16         d4, d5

; memset(input, 0, 32) -- 32bytes
    vmov.i16        q14, #0

    vswp            d3, d4
    vqdmulh.s16     q3, q2, d0[2]
    vqdmulh.s16     q4, q2, d0[0]

    vqadd.s16       d12, d2, d3             ;a1
    vqsub.s16       d13, d2, d3             ;b1

    vmov            q15, q14

    vshr.s16        q3, q3, #1
    vshr.s16        q4, q4, #1

    vqadd.s16       q3, q3, q2              ;modify since sinpi8sqrt2 > 65536/2 (negtive number)
    vqadd.s16       q4, q4, q2

    vqsub.s16       d10, d6, d9             ;c1
    vqadd.s16       d11, d7, d8             ;d1

    vqadd.s16       d2, d12, d11
    vqadd.s16       d3, d13, d10
    vqsub.s16       d4, d13, d10
    vqsub.s16       d5, d12, d11

    vst1.16         {q14, q15}, [r0]

    vrshr.s16       d2, d2, #3
    vrshr.s16       d3, d3, #3
    vrshr.s16       d4, d4, #3
    vrshr.s16       d5, d5, #3

    add             r1, r2, r3
    add             r12, r1, r3
    add             r0, r12, r3

    vtrn.32         d2, d4
    vtrn.32         d3, d5
    vtrn.16         d2, d3
    vtrn.16         d4, d5

    vst1.16         {d2}, [r2]
    vst1.16         {d3}, [r1]
    vst1.16         {d4}, [r12]
    vst1.16         {d5}, [r0]

    bx             lr

    ENDP

;-----------------
    AREA    dcidct4x4_dat, DATA, READWRITE          ;read/write by default
;Data section with name data_area is specified. DCD reserves space in memory for 48 data.
;One word each is reserved. Label filter_coeff can be used to access the data.
;Data address: filter_coeff, filter_coeff+4, filter_coeff+8 ...
_dcidct_coeff_
    DCD     dcidct_coeff
dcidct_coeff
    DCD     0x4e7b4e7b, 0x8a8c8a8c

;20091, 20091, 35468, 35468

    END
