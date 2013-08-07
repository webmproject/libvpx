;
;  Copyright (c) 2013 The WebM project authors. All Rights Reserved.
;
;  Use of this source code is governed by a BSD-style license
;  that can be found in the LICENSE file in the root of the source
;  tree. An additional intellectual property rights grant can be found
;  in the file PATENTS.  All contributing project authors may
;  be found in the AUTHORS file in the root of the source tree.
;

    EXPORT  |vp9_short_idct8x8_add_neon|
    ARM
    REQUIRE8
    PRESERVE8

    AREA ||.text||, CODE, READONLY, ALIGN=2

    ; Parallel 1D IDCT on all the columns of a 8x8 16bit data matrix which are
    ; loaded in q8-q15. The output will be stored back into q8-q15 registers.
    ; This macro will touch q0-q7 registers and use them as buffer during
    ; calculation.
    MACRO
    IDCT8x8_1D
    ; stage 1
    vdup.16         d0, r3                    ; duplicate cospi_28_64
    vdup.16         d1, r4                    ; duplicate cospi_4_64

    ; input[1] * cospi_28_64
    vmull.s16       q2, d18, d0
    vmull.s16       q3, d19, d0

    ; input[7] * cospi_4_64
    vmull.s16       q4, d30, d1
    vmull.s16       q5, d31, d1

    ; input[1]*cospi_28_64-input[7]*cospi_4_64
    vsub.s32        q6, q2, q4
    vsub.s32        q7, q3, q5

    ; dct_const_round_shift(input_dc * cospi_16_64)
    vqrshrn.s32     d8, q6, #14               ; >> 14
    vqrshrn.s32     d9, q7, #14               ; >> 14

    ; input[1] * cospi_4_64
    vmull.s16       q2, d18, d1
    vmull.s16       q3, d19, d1

    ; input[7] * cospi_28_64
    vmull.s16       q1, d30, d0
    vmull.s16       q5, d31, d0

    ; input[1]*cospi_4_64+input[7]*cospi_28_64
    vadd.s32        q2, q2, q1
    vadd.s32        q3, q3, q5

    ; dct_const_round_shift(input_dc * cospi_16_64)
    vqrshrn.s32     d14, q2, #14              ; >> 14
    vqrshrn.s32     d15, q3, #14              ; >> 14

    vdup.16         d0, r5                    ; duplicate cospi_12_64
    vdup.16         d1, r6                    ; duplicate cospi_20_64

    ; input[5] * cospi_12_64
    vmull.s16       q2, d26, d0
    vmull.s16       q3, d27, d0

    ; input[3] * cospi_20_64
    vmull.s16       q5, d22, d1
    vmull.s16       q6, d23, d1

    ; input[5] * cospi_12_64 - input[3] * cospi_20_64
    vsub.s32        q2, q2, q5
    vsub.s32        q3, q3, q6

    ; dct_const_round_shift(input_dc * cospi_16_64)
    vqrshrn.s32     d10, q2, #14              ; >> 14
    vqrshrn.s32     d11, q3, #14              ; >> 14

    ; input[5] * cospi_20_64
    vmull.s16       q2, d26, d1
    vmull.s16       q3, d27, d1

    ; input[3] * cospi_12_64
    vmull.s16       q9, d22, d0
    vmull.s16       q15, d23, d0

    ; input[5] * cospi_20_64 + input[3] * cospi_12_64
    vadd.s32        q0, q2, q9
    vadd.s32        q1, q3, q15

    ; dct_const_round_shift(input_dc * cospi_16_64)
    vqrshrn.s32     d12, q0, #14              ; >> 14
    vqrshrn.s32     d13, q1, #14              ; >> 14

    ; stage 2 & stage 3 - even half
    vdup.16         d0, r7                    ; duplicate cospi_16_64

    ; input[0] * cospi_16_64
    vmull.s16       q2, d16, d0
    vmull.s16       q3, d17, d0

    ; input[2] * cospi_16_64
    vmull.s16       q9,  d24, d0
    vmull.s16       q11, d25, d0

    ; (input[0] + input[2]) * cospi_16_64
    vadd.s32        q9, q2, q9
    vadd.s32        q11, q3, q11

    ; dct_const_round_shift(input_dc * cospi_16_64)
    vqrshrn.s32     d18, q9, #14              ; >> 14
    vqrshrn.s32     d19, q11, #14             ; >> 14

    ; input[0] * cospi_16_64
    vmull.s16       q2, d16, d0
    vmull.s16       q3, d17, d0

    ; input[2] * cospi_16_64
    vmull.s16       q13,  d24, d0
    vmull.s16       q15, d25, d0

    ; (input[0] - input[2]) * cospi_16_64
    vsub.s32        q2, q2, q13
    vsub.s32        q3, q3, q15

    ; dct_const_round_shift(input_dc * cospi_16_64)
    vqrshrn.s32     d22, q2, #14              ; >> 14
    vqrshrn.s32     d23, q3, #14              ; >> 14

    ; input[1] * cospi_24_64 - input[3] * cospi_8_64
    vdup.16         d0, r8                    ; duplicate cospi_24_64
    vdup.16         d1, r9                    ; duplicate cospi_8_64

    ; input[1] * cospi_24_64
    vmull.s16       q2, d20, d0
    vmull.s16       q3, d21, d0

    ; input[3] * cospi_8_64
    vmull.s16       q13, d28, d1
    vmull.s16       q15, d29, d1

    ; input[1] * cospi_24_64 - input[3] * cospi_8_64
    vsub.s32        q2, q2, q13
    vsub.s32        q3, q3, q15

    ; dct_const_round_shift(input_dc * cospi_16_64)
    vqrshrn.s32     d26, q2, #14              ; >> 14
    vqrshrn.s32     d27, q3, #14              ; >> 14

    ; input[1] * cospi_8_64
    vmull.s16       q2, d20, d1
    vmull.s16       q3, d21, d1

    ; input[3] * cospi_24_64
    vmull.s16       q8, d28, d0
    vmull.s16       q10, d29, d0

    ; input[1] * cospi_8_64 + input[3] * cospi_24_64
    vadd.s32        q0, q2, q8
    vadd.s32        q1, q3, q10

    ; dct_const_round_shift(input_dc * cospi_16_64)
    vqrshrn.s32     d30, q0, #14              ; >> 14
    vqrshrn.s32     d31, q1, #14              ; >> 14


    vadd.s16        q0, q9, q15               ; output[0] = step[0] + step[3]
    vadd.s16        q1, q11, q13              ; output[1] = step[1] + step[2]
    vsub.s16        q2, q11, q13              ; output[2] = step[1] - step[2]
    vsub.s16        q3, q9, q15               ; output[3] = step[0] - step[3]

    ; stage 2 - odd half
    vsub.s16        q13, q4, q5               ; step2[5] = step1[4] - step1[5]
    vadd.s16        q4, q4, q5                ; step2[4] = step1[4] + step1[5]
    vsub.s16        q14, q7, q6               ; step2[6] = -step1[6] + step1[7]
    vadd.s16        q7, q7, q6                ; step2[7] = step1[6] + step1[7]

    ; stage 3 -odd half
    vdup.16         d16, r7                   ; duplicate cospi_16_64

    ; step2[6] * cospi_16_64
    vmull.s16       q9, d28, d16
    vmull.s16       q10, d29, d16

    ; step2[5] * cospi_16_64
    vmull.s16       q11, d26, d16
    vmull.s16       q12, d27, d16

    ; (step2[6] - step2[5]) * cospi_16_64
    vsub.s32        q9, q9, q11
    vsub.s32        q10, q10, q12

    ; dct_const_round_shift(input_dc * cospi_16_64)
    vqrshrn.s32     d10, q9, #14              ; >> 14
    vqrshrn.s32     d11, q10, #14             ; >> 14

    ; step2[6] * cospi_16_64
    vmull.s16       q9, d28, d16
    vmull.s16       q10, d29, d16

    ; step2[5] * cospi_16_64
    vmull.s16       q11, d26, d16
    vmull.s16       q12, d27, d16

    ; (step2[5] + step2[6]) * cospi_16_64
    vadd.s32        q9, q9, q11
    vadd.s32        q10, q10, q12

    ; dct_const_round_shift(input_dc * cospi_16_64)
    vqrshrn.s32     d12, q9, #14              ; >> 14
    vqrshrn.s32     d13, q10, #14             ; >> 14

    ; stage 4
    vadd.s16        q8, q0, q7                ; output[0] = step1[0] + step1[7];
    vadd.s16        q9, q1, q6                ; output[1] = step1[1] + step1[6];
    vadd.s16        q10, q2, q5               ; output[2] = step1[2] + step1[5];
    vadd.s16        q11, q3, q4               ; output[3] = step1[3] + step1[4];
    vsub.s16        q12, q3, q4               ; output[4] = step1[3] - step1[4];
    vsub.s16        q13, q2, q5               ; output[5] = step1[2] - step1[5];
    vsub.s16        q14, q1, q6               ; output[6] = step1[1] - step1[6];
    vsub.s16        q15, q0, q7               ; output[7] = step1[0] - step1[7];
    MEND

    ; Transpose a 8x8 16bit data matrix. Datas are loaded in q8-q15.
    MACRO
    TRANSPOSE8X8
    vswp            d17, d24
    vswp            d23, d30
    vswp            d21, d28
    vswp            d19, d26
    vtrn.32         q8, q10
    vtrn.32         q9, q11
    vtrn.32         q12, q14
    vtrn.32         q13, q15
    vtrn.16         q8, q9
    vtrn.16         q10, q11
    vtrn.16         q12, q13
    vtrn.16         q14, q15
    MEND

    AREA    Block, CODE, READONLY ; name this block of code
;void vp9_short_idct8x8_add_neon(int16_t *input, uint8_t *dest, int dest_stride)
;
; r0  int16_t input
; r1  uint8_t *dest
; r2  int dest_stride)

|vp9_short_idct8x8_add_neon| PROC
    push            {r4-r9}
    vld1.s16        {q8}, [r0]!
    vld1.s16        {q9}, [r0]!
    vld1.s16        {q10}, [r0]!
    vld1.s16        {q11}, [r0]!
    vld1.s16        {q12}, [r0]!
    vld1.s16        {q13}, [r0]!
    vld1.s16        {q14}, [r0]!
    vld1.s16        {q15}, [r0]!

    ; transpose the input data
    TRANSPOSE8X8

    ; generate  cospi_28_64 = 3196
    mov             r3, #0x0c00
    add             r3, #0x7c

    ; generate cospi_4_64  = 16069
    mov             r4, #0x3e00
    add             r4, #0xc5

    ; generate cospi_12_64 = 13623
    mov             r5, #0x3500
    add             r5, #0x37

    ; generate cospi_20_64 = 9102
    mov             r6, #0x2300
    add             r6, #0x8e

    ; generate cospi_16_64 = 11585
    mov             r7, #0x2d00
    add             r7, #0x41

    ; generate cospi_24_64 = 6270
    mov             r8, #0x1800
    add             r8, #0x7e

    ; generate cospi_8_64 = 15137
    mov             r9, #0x3b00
    add             r9, #0x21

    ; First transform rows
    IDCT8x8_1D

    ; Transpose the matrix
    TRANSPOSE8X8

    ; Then transform columns
    IDCT8x8_1D

    ; ROUND_POWER_OF_TWO(temp_out[j], 5)
    vrshr.s16       q8, q8, #5
    vrshr.s16       q9, q9, #5
    vrshr.s16       q10, q10, #5
    vrshr.s16       q11, q11, #5
    vrshr.s16       q12, q12, #5
    vrshr.s16       q13, q13, #5
    vrshr.s16       q14, q14, #5
    vrshr.s16       q15, q15, #5

    ; save dest pointer
    mov             r0, r1

    ; load destination data
    vld1.64         {d0}, [r1], r2
    vld1.64         {d1}, [r1], r2
    vld1.64         {d2}, [r1], r2
    vld1.64         {d3}, [r1], r2
    vld1.64         {d4}, [r1], r2
    vld1.64         {d5}, [r1], r2
    vld1.64         {d6}, [r1], r2
    vld1.64         {d7}, [r1]

    ; ROUND_POWER_OF_TWO(temp_out[j], 5) + dest[j * dest_stride + i]
    vaddw.u8        q8, q8, d0
    vaddw.u8        q9, q9, d1
    vaddw.u8        q10, q10, d2
    vaddw.u8        q11, q11, d3
    vaddw.u8        q12, q12, d4
    vaddw.u8        q13, q13, d5
    vaddw.u8        q14, q14, d6
    vaddw.u8        q15, q15, d7

    ; clip_pixel
    vqmovun.s16     d0, q8
    vqmovun.s16     d1, q9
    vqmovun.s16     d2, q10
    vqmovun.s16     d3, q11
    vqmovun.s16     d4, q12
    vqmovun.s16     d5, q13
    vqmovun.s16     d6, q14
    vqmovun.s16     d7, q15

    ; store the data
    vst1.64         {d0}, [r0], r2
    vst1.64         {d1}, [r0], r2
    vst1.64         {d2}, [r0], r2
    vst1.64         {d3}, [r0], r2
    vst1.64         {d4}, [r0], r2
    vst1.64         {d5}, [r0], r2
    vst1.64         {d6}, [r0], r2
    vst1.64         {d7}, [r0], r2

    pop             {r4-r9}
    bx              lr
    ENDP  ; |vp9_short_idct8x8_add_neon|

    END
