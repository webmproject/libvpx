;
;  Copyright (c) 2013 The WebM project authors. All Rights Reserved.
;
;  Use of this source code is governed by a BSD-style license
;  that can be found in the LICENSE file in the root of the source
;  tree. An additional intellectual property rights grant can be found
;  in the file PATENTS.  All contributing project authors may
;  be found in the AUTHORS file in the root of the source tree.
;


    ; These functions are only valid when:
    ; x_step_q4 == 16
    ; w%4 == 0
    ; h%4 == 0
    ; taps == 8
    ; VP9_FILTER_WEIGHT == 128
    ; VP9_FILTER_SHIFT == 7

    EXPORT  |vp9_convolve8_horiz_neon|
    EXPORT  |vp9_convolve8_vert_neon|
    IMPORT  |vp9_convolve8_horiz_c|
    IMPORT  |vp9_convolve8_vert_c|
    ARM
    REQUIRE8
    PRESERVE8

    AREA ||.text||, CODE, READONLY, ALIGN=2

    ; Multiply and accumulate by q0
    MACRO
    MULTIPLY_BY_Q0 $dst, $src0, $src1, $src2, $src3, $src4, $src5, $src6, $src7
    vmull.s16 $dst, $src0, d0[0]
    vmlal.s16 $dst, $src1, d0[1]
    vmlal.s16 $dst, $src2, d0[2]
    vmlal.s16 $dst, $src3, d0[3]
    vmlal.s16 $dst, $src4, d1[0]
    vmlal.s16 $dst, $src5, d1[1]
    vmlal.s16 $dst, $src6, d1[2]
    vmlal.s16 $dst, $src7, d1[3]
    MEND

; r0    const uint8_t *src
; r1    int src_stride
; r2    uint8_t *dst
; r3    int dst_stride
; sp[]const int16_t *filter_x
; sp[]int x_step_q4
; sp[]const int16_t *filter_y ; unused
; sp[]int y_step_q4           ; unused
; sp[]int w
; sp[]int h

|vp9_convolve8_horiz_neon| PROC
    ldr             r12, [sp, #4]           ; x_step_q4
    cmp             r12, #16
    bne             vp9_convolve8_horiz_c

    push            {r4-r10, lr}

    sub             r0, r0, #3              ; adjust for taps

    ldr             r5, [sp, #32]           ; filter_x
    ldr             r6, [sp, #48]           ; w
    ldr             r7, [sp, #52]           ; h

    vld1.s16        {q0}, [r5]              ; filter_x

    add             r8, r1, r1, lsl #1      ; src_stride * 3
    add             r8, r8, #4              ; src_stride * 3 + 4
    rsb             r8, r8, #0              ; reset for src

    add             r4, r3, r3, lsl #1      ; dst_stride * 3
    sub             r4, r4, #4              ; dst_stride * 3 - 4
    rsb             r4, r4, #0              ; reset for dst

    sub             r9, r1, #8              ; post increment for src load

    rsb             r1, r6, r1, lsl #2      ; reset src for outer loop
    rsb             r12, r6, r3, lsl #2     ; reset dst for outer loop

    mov             r10, r6                 ; w loop counter

loop_horiz
    vld1.8          {d24}, [r0]!
    vld3.u8         {d28[0], d29[0], d30[0]}, [r0], r9

    vld1.8          {d25}, [r0]!
    vld3.u8         {d28[1], d29[1], d30[1]}, [r0], r9

    vld1.8          {d26}, [r0]!
    vld3.u8         {d28[2], d29[2], d30[2]}, [r0], r9

    vld1.8          {d27}, [r0]!
    vld3.u8         {d28[3], d29[3], d30[3]}, [r0], r8

    vtrn.16         q12, q13
    vtrn.8          d24, d25
    vtrn.8          d26, d27

    ; extract to s16
    vmovl.u8        q8, d24
    vmovl.u8        q9, d25
    vmovl.u8        q10, d26
    vmovl.u8        q11, d27
    vtrn.32         d28, d29 ; only the first half is populated
    vmovl.u8        q12, d28
    vmovl.u8        q13, d30

    ; src[] * filter_x
    MULTIPLY_BY_Q0 q1, d16, d18, d20, d22, d17, d19, d21, d23
    MULTIPLY_BY_Q0 q2, d18, d20, d22, d17, d19, d21, d23, d24
    MULTIPLY_BY_Q0 q14, d20, d22, d17, d19, d21, d23, d24, d25
    MULTIPLY_BY_Q0 q15, d22, d17, d19, d21, d23, d24, d25, d26

    ; += 64 >> 7
    vqrshrun.s32    d2, q1, #7
    vqrshrun.s32    d3, q2, #7
    vqrshrun.s32    d4, q14, #7
    vqrshrun.s32    d5, q15, #7

    ; saturate
    vqmovn.u16      d2, q1
    vqmovn.u16      d3, q2

    ; transpose
    vtrn.16         d2, d3
    vtrn.32         d2, d3
    vtrn.8          d2, d3

    vst1.u32        {d2[0]}, [r2], r3
    vst1.u32        {d3[0]}, [r2], r3
    vst1.u32        {d2[1]}, [r2], r3
    vst1.u32        {d3[1]}, [r2], r4

    subs            r6, r6, #4              ; w -= 4
    bgt             loop_horiz

    ; outer loop
    mov             r6, r10                 ; restore w counter
    add             r0, r0, r1              ; src += src_stride * 4 - w
    add             r2, r2, r12             ; dst += dst_stride * 4 - w
    subs            r7, r7, #4              ; h -= 4
    bgt loop_horiz

    pop             {r4-r10, pc}

    ENDP

|vp9_convolve8_vert_neon| PROC
    ldr             r12, [sp, #12]
    cmp             r12, #16
    bne             vp9_convolve8_vert_c

    push            {r4-r10, lr}

    ; adjust for taps
    sub             r0, r0, r1
    sub             r0, r0, r1, lsl #1

    ldr             r7, [sp, #40]           ; filter_y
    ldr             r8, [sp, #48]           ; w
    ldr             r9, [sp, #52]           ; h

    vld1.s16        {q0}, [r7]              ; filter_y

    mov             r5, r1, lsl #1          ; src_stride * 2
    add             r5, r5, r1, lsl #3      ; src_stride * 10
    sub             r5, r5, #4              ; src_stride * 10 + 4
    rsb             r5, r5, #0              ; reset for src

    add             r6, r3, r3, lsl #1      ; dst_stride * 3
    sub             r6, r6, #4              ; dst_stride * 3 - 4
    rsb             r6, r6, #0              ; reset for dst

    rsb             r7, r8, r1, lsl #2      ; reset src for outer loop
    rsb             r12, r8, r3, lsl #2     ; reset dst for outer loop

    mov             r10, r8                 ; w loop counter

loop_vert
    ; always process a 4x4 block at a time
    vld1.u32        {d16[0]}, [r0], r1
    vld1.u32        {d16[1]}, [r0], r1
    vld1.u32        {d18[0]}, [r0], r1
    vld1.u32        {d18[1]}, [r0], r1
    vld1.u32        {d20[0]}, [r0], r1
    vld1.u32        {d20[1]}, [r0], r1
    vld1.u32        {d22[0]}, [r0], r1
    vld1.u32        {d22[1]}, [r0], r1
    vld1.u32        {d24[0]}, [r0], r1
    vld1.u32        {d24[1]}, [r0], r1
    vld1.u32        {d26[0]}, [r0], r5

    ; extract to s16
    vmovl.u8        q8, d16
    vmovl.u8        q9, d18
    vmovl.u8        q10, d20
    vmovl.u8        q11, d22
    vmovl.u8        q12, d24
    vmovl.u8        q13, d26

    ; src[] * filter_y
    MULTIPLY_BY_Q0 q1, d16, d17, d18, d19, d20, d21, d22, d23
    MULTIPLY_BY_Q0 q2, d17, d18, d19, d20, d21, d22, d23, d24
    MULTIPLY_BY_Q0 q14, d18, d19, d20, d21, d22, d23, d24, d25
    MULTIPLY_BY_Q0 q15, d19, d20, d21, d22, d23, d24, d25, d26

    ; += 64 >> 7
    vqrshrun.s32    d2, q1, #7
    vqrshrun.s32    d3, q2, #7
    vqrshrun.s32    d4, q14, #7
    vqrshrun.s32    d5, q15, #7

    ; saturate
    vqmovn.u16      d2, q1
    vqmovn.u16      d3, q2

    vst1.u32        {d2[0]}, [r2], r3
    vst1.u32        {d2[1]}, [r2], r3
    vst1.u32        {d3[0]}, [r2], r3
    vst1.u32        {d3[1]}, [r2], r6

    subs            r8, r8, #4              ; w -= 4
    bgt             loop_vert

    ; outer loop
    mov             r8, r10                 ; restore w counter
    add             r0, r0, r7              ; src += 4 * src_stride - w
    add             r2, r2, r12             ; dst += 4 * dst_stride - w
    subs            r9, r9, #4              ; h -= 4
    bgt             loop_vert

    pop             {r4-r10, pc}

    ENDP
    END
