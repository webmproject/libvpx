;
;  Copyright (c) 2011 The WebM project authors. All Rights Reserved.
;
;  Use of this source code is governed by a BSD-style license
;  that can be found in the LICENSE file in the root of the source
;  tree. An additional intellectual property rights grant can be found
;  in the file PATENTS.  All contributing project authors may
;  be found in the AUTHORS file in the root of the source tree.
;


    EXPORT  |vp8_fast_quantize_b_neon|

    INCLUDE asm_enc_offsets.asm

    ARM
    REQUIRE8
    PRESERVE8

    AREA ||.text||, CODE, READONLY, ALIGN=4


;void vp8_fast_quantize_b_c(BLOCK *b, BLOCKD *d)
|vp8_fast_quantize_b_neon| PROC

    stmfd           sp!, {r4-r7}

    ldr             r3, [r0, #vp8_block_coeff]
    ldr             r4, [r0, #vp8_block_quant_fast]
    ldr             r5, [r0, #vp8_block_round]

    vld1.16         {q0, q1}, [r3@128]  ; load z
    vorr.s16        q14, q0, q1         ; check if all zero (step 1)
    ldr             r6, [r1, #vp8_blockd_qcoeff]
    ldr             r7, [r1, #vp8_blockd_dqcoeff]
    vorr.s16        d28, d28, d29       ; check if all zero (step 2)

    vabs.s16        q12, q0             ; calculate x = abs(z)
    vabs.s16        q13, q1

    ;right shift 15 to get sign, all 0 if it is positive, all 1 if it is negative
    vshr.s16        q2, q0, #15         ; sz
    vmov            r2, r3, d28         ; check if all zero (step 3)
    vshr.s16        q3, q1, #15

    vld1.s16        {q14, q15}, [r5@128]; load round_ptr [0-15]
    vld1.s16        {q8, q9}, [r4@128]  ; load quant_ptr [0-15]

    vadd.s16        q12, q14            ; x + Round
    vadd.s16        q13, q15

    ldr             r0, _inv_zig_zag_   ; load ptr of inverse zigzag table

    vqdmulh.s16     q12, q8             ; y = ((Round+abs(z)) * Quant) >> 16
    vqdmulh.s16     q13, q9

    vld1.16         {q10, q11}, [r0@128]; load inverse scan order

    vceq.s16        q8, q8              ; set q8 to all 1

    ldr             r4, [r1, #vp8_blockd_dequant]

    vshr.s16        q12, #1             ; right shift 1 after vqdmulh
    vshr.s16        q13, #1

    orr             r2, r2, r3          ; check if all zero (step 4)
    cmp             r2, #0              ; check if all zero (step 5)
    beq             zero_output         ; check if all zero (step 6)

    ;modify data to have its original sign
    veor.s16        q12, q2             ; y^sz
    veor.s16        q13, q3

    vsub.s16        q12, q2             ; x1=(y^sz)-sz = (y^sz)-(-1) (2's complement)
    vsub.s16        q13, q3

    vld1.s16        {q2, q3}, [r4@128]  ; load dequant_ptr[i]

    vtst.16         q14, q12, q8        ; now find eob
    vtst.16         q15, q13, q8        ; non-zero element is set to all 1

    vst1.s16        {q12, q13}, [r6@128]; store: qcoeff = x1

    vand            q10, q10, q14       ; get all valid numbers from scan array
    vand            q11, q11, q15


    vmax.u16        q0, q10, q11        ; find maximum value in q0, q1
    vmax.u16        d0, d0, d1
    vmovl.u16       q0, d0

    vmul.s16        q2, q12             ; x * Dequant
    vmul.s16        q3, q13

    vmax.u32        d0, d0, d1
    vpmax.u32       d0, d0, d0

    vst1.s16        {q2, q3}, [r7@128]  ; store dqcoeff = x * Dequant

    vmov.32         r0, d0[0]           ; this instruction takes 1+13 cycles
                                        ; if we have vfp, we could use
                                        ; vstr      s0, [r1, #vp8_blockd_eob]
    str             r0, [r1, #vp8_blockd_eob]

    ldmfd           sp!, {r4-r7}
    bx              lr

zero_output
    str             r2, [r1, #vp8_blockd_eob]
    vst1.s16        {q0, q1}, [r6@128]  ; qcoeff = 0
    vst1.s16        {q0, q1}, [r7@128]  ; dqcoeff = 0

    ldmfd           sp!, {r4-r7}
    bx              lr

    ENDP

; default inverse zigzag table is defined in vp8/common/entropy.c
_inv_zig_zag_
    DCD inv_zig_zag

    ALIGN 16    ; enable use of @128 bit aligned loads
inv_zig_zag
    DCW 0x0001, 0x0002, 0x0006, 0x0007
    DCW 0x0003, 0x0005, 0x0008, 0x000d
    DCW 0x0004, 0x0009, 0x000c, 0x000e
    DCW 0x000a, 0x000b, 0x000f, 0x0010

    END

