;
;  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
;
;  Use of this source code is governed by a BSD-style license
;  that can be found in the LICENSE file in the root of the source
;  tree. An additional intellectual property rights grant can be found
;  in the file PATENTS.  All contributing project authors may
;  be found in the AUTHORS file in the root of the source tree.
;


    EXPORT  |vp8_yv12_extend_frame_borders_neon|
    ARM
    REQUIRE8
    PRESERVE8

    INCLUDE asm_com_offsets.asm

    AREA ||.text||, CODE, READONLY, ALIGN=2
;void vp8_yv12_extend_frame_borders_neon (YV12_BUFFER_CONFIG *ybf);
;Note: this is VP8 function, which has border=32 and 16. Internal y_width and y_height
; are always multiples of 16.

|vp8_yv12_extend_frame_borders_neon| PROC
    push            {r4 - r10, lr}
    vpush           {d8 - d15}

    ;Not need to load y_width, since: y_width = y_stride - 2*border
    ldr             r3, [r0, #yv12_buffer_config_border]
    ldr             r1, [r0, #yv12_buffer_config_y_buffer]       ;srcptr1
    ldr             r4, [r0, #yv12_buffer_config_y_height]
    ldr             lr, [r0, #yv12_buffer_config_y_stride]

    cmp             r3, #16
    beq             b16_extend_frame_borders

;=======================
b32_extend_frame_borders
;border = 32
;=======================
;Border copy for Y plane
;copy the left and right most columns out
    sub             r5, r1, r3              ;destptr1
    add             r6, r1, lr
    sub             r6, r6, r3, lsl #1      ;destptr2
    sub             r2, r6, #1              ;srcptr2

    ;Do four rows at one time
    mov             r12, r4, lsr #2

copy_left_right_y
    vld1.8          {d0[], d1[]}, [r1], lr
    vld1.8          {d4[], d5[]}, [r2], lr
    vld1.8          {d8[], d9[]}, [r1], lr
    vld1.8          {d12[], d13[]}, [r2], lr
    vld1.8          {d16[], d17[]},  [r1], lr
    vld1.8          {d20[], d21[]}, [r2], lr
    vld1.8          {d24[], d25[]}, [r1], lr
    vld1.8          {d28[], d29[]}, [r2], lr

    vmov            q1, q0
    vmov            q3, q2
    vmov            q5, q4
    vmov            q7, q6
    vmov            q9, q8
    vmov            q11, q10
    vmov            q13, q12
    vmov            q15, q14

    subs            r12, r12, #1

    vst1.8          {q0, q1}, [r5], lr
    vst1.8          {q2, q3}, [r6], lr
    vst1.8          {q4, q5}, [r5], lr
    vst1.8          {q6, q7}, [r6], lr
    vst1.8          {q8, q9}, [r5], lr
    vst1.8          {q10, q11}, [r6], lr
    vst1.8          {q12, q13}, [r5], lr
    vst1.8          {q14, q15}, [r6], lr

    bne             copy_left_right_y

;Now copy the top and bottom source lines into each line of the respective borders
    ldr             r7, [r0, #yv12_buffer_config_y_buffer]       ;srcptr1
    mul             r8, r3, lr

    mov             r12, lr, lsr #7

    sub             r6, r1, r3              ;destptr2
    sub             r2, r6, lr              ;srcptr2
    sub             r1, r7, r3              ;srcptr1
    sub             r5, r1, r8              ;destptr1

copy_top_bottom_y
    vld1.8          {q0, q1}, [r1]!
    vld1.8          {q8, q9}, [r2]!
    vld1.8          {q2, q3}, [r1]!
    vld1.8          {q10, q11}, [r2]!
    vld1.8          {q4, q5}, [r1]!
    vld1.8          {q12, q13}, [r2]!
    vld1.8          {q6, q7}, [r1]!
    vld1.8          {q14, q15}, [r2]!

    mov             r7, r3

top_bottom_32
    subs            r7, r7, #1

    vst1.8          {q0, q1}, [r5]!
    vst1.8          {q8, q9}, [r6]!
    vst1.8          {q2, q3}, [r5]!
    vst1.8          {q10, q11}, [r6]!
    vst1.8          {q4, q5}, [r5]!
    vst1.8          {q12, q13}, [r6]!
    vst1.8          {q6, q7}, [r5]!
    vst1.8          {q14, q15}, [r6]!

    add             r5, r5, lr
    sub             r5, r5, #128
    add             r6, r6, lr
    sub             r6, r6, #128

    bne             top_bottom_32

    sub             r5, r1, r8
    add             r6, r2, lr

    subs            r12, r12, #1
    bne             copy_top_bottom_y

    mov             r7, lr, lsr #4              ;check to see if extra copy is needed
    ands            r7, r7, #0x7
    bne             extra_top_bottom_y
end_of_border_copy_y

;Border copy for U, V planes
    ldr             r1, [r0, #yv12_buffer_config_u_buffer]       ;srcptr1
    mov             lr, lr, lsr #1              ;uv_stride
    mov             r3, r3, lsr #1              ;border
    mov             r4, r4, lsr #1              ;uv_height
    mov             r8, r8, lsr #2

    mov             r10, #2

;copy the left and right most columns out
border_copy_uv
    sub             r5, r1, r3              ;destptr1
    add             r6, r1, lr
    sub             r6, r6, r3, lsl #1      ;destptr2
    sub             r2, r6, #1              ;srcptr2

    mov             r7, r1

    ;Do eight rows at one time
    mov             r12, r4, lsr #3

copy_left_right_uv
    vld1.8          {d0[], d1[]}, [r1], lr
    vld1.8          {d2[], d3[]}, [r2], lr
    vld1.8          {d4[], d5[]}, [r1], lr
    vld1.8          {d6[], d7[]}, [r2], lr
    vld1.8          {d8[], d9[]},  [r1], lr
    vld1.8          {d10[], d11[]}, [r2], lr
    vld1.8          {d12[], d13[]}, [r1], lr
    vld1.8          {d14[], d15[]}, [r2], lr
    vld1.8          {d16[], d17[]}, [r1], lr
    vld1.8          {d18[], d19[]}, [r2], lr
    vld1.8          {d20[], d21[]}, [r1], lr
    vld1.8          {d22[], d23[]}, [r2], lr
    vld1.8          {d24[], d25[]},  [r1], lr
    vld1.8          {d26[], d27[]}, [r2], lr
    vld1.8          {d28[], d29[]}, [r1], lr
    vld1.8          {d30[], d31[]}, [r2], lr

    subs            r12, r12, #1

    vst1.8          {q0}, [r5], lr
    vst1.8          {q1}, [r6], lr
    vst1.8          {q2}, [r5], lr
    vst1.8          {q3}, [r6], lr
    vst1.8          {q4}, [r5], lr
    vst1.8          {q5}, [r6], lr
    vst1.8          {q6}, [r5], lr
    vst1.8          {q7}, [r6], lr
    vst1.8          {q8}, [r5], lr
    vst1.8          {q9}, [r6], lr
    vst1.8          {q10}, [r5], lr
    vst1.8          {q11}, [r6], lr
    vst1.8          {q12}, [r5], lr
    vst1.8          {q13}, [r6], lr
    vst1.8          {q14}, [r5], lr
    vst1.8          {q15}, [r6], lr

    bne             copy_left_right_uv

;Now copy the top and bottom source lines into each line of the respective borders
    mov             r12, lr, lsr #6

    sub             r6, r1, r3              ;destptr2
    sub             r2, r6, lr              ;srcptr2
    sub             r1, r7, r3              ;srcptr1
    sub             r5, r1, r8              ;destptr1

copy_top_bottom_uv
    vld1.8          {q0, q1}, [r1]!
    vld1.8          {q8, q9}, [r2]!
    vld1.8          {q2, q3}, [r1]!
    vld1.8          {q10, q11}, [r2]!

    mov             r7, r3

top_bottom_16
    subs            r7, r7, #1

    vst1.8          {q0, q1}, [r5]!
    vst1.8          {q8, q9}, [r6]!
    vst1.8          {q2, q3}, [r5]!
    vst1.8          {q10, q11}, [r6]!

    add             r5, r5, lr
    sub             r5, r5, #64
    add             r6, r6, lr
    sub             r6, r6, #64

    bne             top_bottom_16

    sub             r5, r1, r8
    add             r6, r2, lr

    subs            r12, r12, #1
    bne             copy_top_bottom_uv

    mov             r7, lr, lsr #3              ;check to see if extra copy is needed
    ands            r7, r7, #0x7
    bne             extra_top_bottom_uv

end_of_border_copy_uv
    subs            r10, r10, #1
    ldrne           r1, [r0, #yv12_buffer_config_v_buffer]       ;srcptr1
    bne             border_copy_uv

    vpop            {d8 - d15}
    pop             {r4 - r10, pc}

;;;;;;;;;;;;;;;;;;;;;;
;extra copy part for Y
extra_top_bottom_y
    vld1.8          {q0}, [r1]!
    vld1.8          {q2}, [r2]!

    mov             r9, r3, lsr #3

extra_top_bottom_32
    subs            r9, r9, #1

    vst1.8          {q0}, [r5], lr
    vst1.8          {q2}, [r6], lr
    vst1.8          {q0}, [r5], lr
    vst1.8          {q2}, [r6], lr
    vst1.8          {q0}, [r5], lr
    vst1.8          {q2}, [r6], lr
    vst1.8          {q0}, [r5], lr
    vst1.8          {q2}, [r6], lr
    vst1.8          {q0}, [r5], lr
    vst1.8          {q2}, [r6], lr
    vst1.8          {q0}, [r5], lr
    vst1.8          {q2}, [r6], lr
    vst1.8          {q0}, [r5], lr
    vst1.8          {q2}, [r6], lr
    vst1.8          {q0}, [r5], lr
    vst1.8          {q2}, [r6], lr
    bne             extra_top_bottom_32

    sub             r5, r1, r8
    add             r6, r2, lr
    subs            r7, r7, #1
    bne             extra_top_bottom_y

    b               end_of_border_copy_y

;extra copy part for UV
extra_top_bottom_uv
    vld1.8          {d0}, [r1]!
    vld1.8          {d8}, [r2]!

    mov             r9, r3, lsr #3

extra_top_bottom_16
    subs            r9, r9, #1

    vst1.8          {d0}, [r5], lr
    vst1.8          {d8}, [r6], lr
    vst1.8          {d0}, [r5], lr
    vst1.8          {d8}, [r6], lr
    vst1.8          {d0}, [r5], lr
    vst1.8          {d8}, [r6], lr
    vst1.8          {d0}, [r5], lr
    vst1.8          {d8}, [r6], lr
    vst1.8          {d0}, [r5], lr
    vst1.8          {d8}, [r6], lr
    vst1.8          {d0}, [r5], lr
    vst1.8          {d8}, [r6], lr
    vst1.8          {d0}, [r5], lr
    vst1.8          {d8}, [r6], lr
    vst1.8          {d0}, [r5], lr
    vst1.8          {d8}, [r6], lr
    bne             extra_top_bottom_16

    sub             r5, r1, r8
    add             r6, r2, lr
    subs            r7, r7, #1
    bne             extra_top_bottom_uv

    b               end_of_border_copy_uv


;=======================
b16_extend_frame_borders
;border = 16
;=======================
;Border copy for Y plane
;copy the left and right most columns out
    sub             r5, r1, r3              ;destptr1
    add             r6, r1, lr
    sub             r6, r6, r3, lsl #1      ;destptr2
    sub             r2, r6, #1              ;srcptr2

    ;Do four rows at one time
    mov             r12, r4, lsr #2

copy_left_right_y_b16
    vld1.8          {d0[], d1[]}, [r1], lr
    vld1.8          {d4[], d5[]}, [r2], lr
    vld1.8          {d8[], d9[]}, [r1], lr
    vld1.8          {d12[], d13[]}, [r2], lr
    vld1.8          {d16[], d17[]},  [r1], lr
    vld1.8          {d20[], d21[]}, [r2], lr
    vld1.8          {d24[], d25[]}, [r1], lr
    vld1.8          {d28[], d29[]}, [r2], lr

    subs            r12, r12, #1

    vst1.8          {q0}, [r5], lr
    vst1.8          {q2}, [r6], lr
    vst1.8          {q4}, [r5], lr
    vst1.8          {q6}, [r6], lr
    vst1.8          {q8}, [r5], lr
    vst1.8          {q10}, [r6], lr
    vst1.8          {q12}, [r5], lr
    vst1.8          {q14}, [r6], lr

    bne             copy_left_right_y_b16

;Now copy the top and bottom source lines into each line of the respective borders
    ldr             r7, [r0, #yv12_buffer_config_y_buffer]       ;srcptr1
    mul             r8, r3, lr

    mov             r12, lr, lsr #7

    sub             r6, r1, r3              ;destptr2
    sub             r2, r6, lr              ;srcptr2
    sub             r1, r7, r3              ;srcptr1
    sub             r5, r1, r8              ;destptr1

copy_top_bottom_y_b16
    vld1.8          {q0, q1}, [r1]!
    vld1.8          {q8, q9}, [r2]!
    vld1.8          {q2, q3}, [r1]!
    vld1.8          {q10, q11}, [r2]!
    vld1.8          {q4, q5}, [r1]!
    vld1.8          {q12, q13}, [r2]!
    vld1.8          {q6, q7}, [r1]!
    vld1.8          {q14, q15}, [r2]!

    mov             r7, r3

top_bottom_16_b16
    subs            r7, r7, #1

    vst1.8          {q0, q1}, [r5]!
    vst1.8          {q8, q9}, [r6]!
    vst1.8          {q2, q3}, [r5]!
    vst1.8          {q10, q11}, [r6]!
    vst1.8          {q4, q5}, [r5]!
    vst1.8          {q12, q13}, [r6]!
    vst1.8          {q6, q7}, [r5]!
    vst1.8          {q14, q15}, [r6]!

    add             r5, r5, lr
    sub             r5, r5, #128
    add             r6, r6, lr
    sub             r6, r6, #128

    bne             top_bottom_16_b16

    sub             r5, r1, r8
    add             r6, r2, lr

    subs            r12, r12, #1
    bne             copy_top_bottom_y_b16

    mov             r7, lr, lsr #4              ;check to see if extra copy is needed
    ands            r7, r7, #0x7
    bne             extra_top_bottom_y_b16
end_of_border_copy_y_b16

;Border copy for U, V planes
    ldr             r1, [r0, #yv12_buffer_config_u_buffer]       ;srcptr1
    mov             lr, lr, lsr #1              ;uv_stride
    mov             r3, r3, lsr #1              ;border
    mov             r4, r4, lsr #1              ;uv_height
    mov             r8, r8, lsr #2

    mov             r10, #2

;copy the left and right most columns out
border_copy_uv_b16
    sub             r5, r1, r3              ;destptr1
    add             r6, r1, lr
    sub             r6, r6, r3, lsl #1      ;destptr2
    sub             r2, r6, #1              ;srcptr2

    mov             r7, r1

    ;Do eight rows at one time
    mov             r12, r4, lsr #3

copy_left_right_uv_b16
    vld1.8          {d0[]}, [r1], lr
    vld1.8          {d2[]}, [r2], lr
    vld1.8          {d4[]}, [r1], lr
    vld1.8          {d6[]}, [r2], lr
    vld1.8          {d8[]},  [r1], lr
    vld1.8          {d10[]}, [r2], lr
    vld1.8          {d12[]}, [r1], lr
    vld1.8          {d14[]}, [r2], lr
    vld1.8          {d16[]}, [r1], lr
    vld1.8          {d18[]}, [r2], lr
    vld1.8          {d20[]}, [r1], lr
    vld1.8          {d22[]}, [r2], lr
    vld1.8          {d24[]},  [r1], lr
    vld1.8          {d26[]}, [r2], lr
    vld1.8          {d28[]}, [r1], lr
    vld1.8          {d30[]}, [r2], lr

    subs            r12, r12, #1

    vst1.8          {d0}, [r5], lr
    vst1.8          {d2}, [r6], lr
    vst1.8          {d4}, [r5], lr
    vst1.8          {d6}, [r6], lr
    vst1.8          {d8}, [r5], lr
    vst1.8          {d10}, [r6], lr
    vst1.8          {d12}, [r5], lr
    vst1.8          {d14}, [r6], lr
    vst1.8          {d16}, [r5], lr
    vst1.8          {d18}, [r6], lr
    vst1.8          {d20}, [r5], lr
    vst1.8          {d22}, [r6], lr
    vst1.8          {d24}, [r5], lr
    vst1.8          {d26}, [r6], lr
    vst1.8          {d28}, [r5], lr
    vst1.8          {d30}, [r6], lr

    bne             copy_left_right_uv_b16

;Now copy the top and bottom source lines into each line of the respective borders
    mov             r12, lr, lsr #6

    sub             r6, r1, r3              ;destptr2
    sub             r2, r6, lr              ;srcptr2
    sub             r1, r7, r3              ;srcptr1
    sub             r5, r1, r8              ;destptr1

copy_top_bottom_uv_b16
    vld1.8          {q0, q1}, [r1]!
    vld1.8          {q8, q9}, [r2]!
    vld1.8          {q2, q3}, [r1]!
    vld1.8          {q10, q11}, [r2]!

    mov             r7, r3

top_bottom_8_b16
    subs            r7, r7, #1

    vst1.8          {q0, q1}, [r5]!
    vst1.8          {q8, q9}, [r6]!
    vst1.8          {q2, q3}, [r5]!
    vst1.8          {q10, q11}, [r6]!

    add             r5, r5, lr
    sub             r5, r5, #64
    add             r6, r6, lr
    sub             r6, r6, #64

    bne             top_bottom_8_b16

    sub             r5, r1, r8
    add             r6, r2, lr

    subs            r12, r12, #1
    bne             copy_top_bottom_uv_b16

    mov             r7, lr, lsr #3              ;check to see if extra copy is needed
    ands            r7, r7, #0x7
    bne             extra_top_bottom_uv_b16

end_of_border_copy_uv_b16
    subs            r10, r10, #1
    ldrne           r1, [r0, #yv12_buffer_config_v_buffer]       ;srcptr1
    bne             border_copy_uv_b16

    vpop            {d8-d15}
    pop             {r4 - r10, pc}

;;;;;;;;;;;;;;;;;;;;;;
;extra copy part for Y
extra_top_bottom_y_b16
    vld1.8          {q0}, [r1]!
    vld1.8          {q2}, [r2]!

    mov             r9, r3, lsr #3

extra_top_bottom_16_b16
    subs            r9, r9, #1

    vst1.8          {q0}, [r5], lr
    vst1.8          {q2}, [r6], lr
    vst1.8          {q0}, [r5], lr
    vst1.8          {q2}, [r6], lr
    vst1.8          {q0}, [r5], lr
    vst1.8          {q2}, [r6], lr
    vst1.8          {q0}, [r5], lr
    vst1.8          {q2}, [r6], lr
    vst1.8          {q0}, [r5], lr
    vst1.8          {q2}, [r6], lr
    vst1.8          {q0}, [r5], lr
    vst1.8          {q2}, [r6], lr
    vst1.8          {q0}, [r5], lr
    vst1.8          {q2}, [r6], lr
    vst1.8          {q0}, [r5], lr
    vst1.8          {q2}, [r6], lr
    bne             extra_top_bottom_16_b16

    sub             r5, r1, r8
    add             r6, r2, lr
    subs            r7, r7, #1
    bne             extra_top_bottom_y_b16

    b               end_of_border_copy_y_b16

;extra copy part for UV
extra_top_bottom_uv_b16
    vld1.8          {d0}, [r1]!
    vld1.8          {d8}, [r2]!

    mov             r9, r3, lsr #3

extra_top_bottom_8_b16
    subs            r9, r9, #1

    vst1.8          {d0}, [r5], lr
    vst1.8          {d8}, [r6], lr
    vst1.8          {d0}, [r5], lr
    vst1.8          {d8}, [r6], lr
    vst1.8          {d0}, [r5], lr
    vst1.8          {d8}, [r6], lr
    vst1.8          {d0}, [r5], lr
    vst1.8          {d8}, [r6], lr
    vst1.8          {d0}, [r5], lr
    vst1.8          {d8}, [r6], lr
    vst1.8          {d0}, [r5], lr
    vst1.8          {d8}, [r6], lr
    vst1.8          {d0}, [r5], lr
    vst1.8          {d8}, [r6], lr
    vst1.8          {d0}, [r5], lr
    vst1.8          {d8}, [r6], lr
    bne             extra_top_bottom_8_b16

    sub             r5, r1, r8
    add             r6, r2, lr
    subs            r7, r7, #1
    bne             extra_top_bottom_uv_b16

    b               end_of_border_copy_uv_b16

    ENDP
    END
