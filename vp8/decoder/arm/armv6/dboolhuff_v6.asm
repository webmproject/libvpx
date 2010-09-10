;
;  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
;
;  Use of this source code is governed by a BSD-style license
;  that can be found in the LICENSE file in the root of the source
;  tree. An additional intellectual property rights grant can be found
;  in the file PATENTS.  All contributing project authors may
;  be found in the AUTHORS file in the root of the source tree.
;


    EXPORT  |vp8_decode_value_v6|
    EXPORT  |vp8dx_start_decode_v6|
    EXPORT  |vp8dx_stop_decode_v6|
    EXPORT  |vp8dx_decode_bool_v6|

    ARM
    REQUIRE8
    PRESERVE8

    INCLUDE vpx_asm_offsets.asm

br      RN  r0
prob    RN  r1
bits    RN  r1
    AREA    |.text|, CODE, READONLY  ; name this block of code

;   int z = 0;
;   int bit;
;   for ( bit=bits-1; bit>=0; bit-- )
;   {
;       z |= (vp8dx_decode_bool(br, 0x80)<<bit);
;   }
;   return z;

;int vp8_decode_value_v6 ( BOOL_DECODER *br, int bits )
|vp8_decode_value_v6| PROC
    stmdb   sp!, {r4 - r6, lr}
    mov     r4, br
    mov     r5, bits
    mov     r6, #0

    subs    r5, r5, #1
    bmi     decode_value_exit

decode_value_loop
    mov     prob, #0x80
    mov     br, r4
    bl      vp8dx_decode_bool_v6_internal     ; needed for conversion to s file
    orr     r6, r6, r0, lsl r5
    subs    r5, r5, #1
    bpl     decode_value_loop

decode_value_exit
    mov     r0, r6
    ldmia   sp!, {r4 - r6, pc}
    ENDP    ; |vp8_decode_value_v6|


;void vp8dx_start_decode_v6 ( BOOL_DECODER *br, unsigned char *source )
|vp8dx_start_decode_v6| PROC
    stmdb   sp!, {r4 - r5, lr}
    mov     r2, #0
    mov     r3, #255

    str     r2, [br, #bool_decoder_lowvalue]
    str     r3, [br, #bool_decoder_range]
    str     r1, [br, #bool_decoder_buffer]

    mov     r3, #8
    mov     r2, #4
    str     r3, [br, #bool_decoder_count]
    str     r2, [br, #bool_decoder_pos]

    ldrb    r2, [r1, #3]
    ldrb    r3, [r1, #2]
    ldrb    r4, [r1, #1]
    ldrb    r5, [r1]

    orr     r1, r2, r3, lsl #8
    orr     r1, r1, r4, lsl #16
    orr     r1, r1, r5, lsl #24

    str     r1, [br, #bool_decoder_value]

    ldmia   sp!, {r4 - r5, pc}
    ENDP    ; |vp8dx_start_decode_v6|


;void vp8dx_stop_decode_v6 ( BOOL_DECODER *bc );
|vp8dx_stop_decode_v6| PROC
    mov     pc, lr
    ENDP    ; |vp8dx_stop_decode_v6|


; bigsplit  RN  r1
; buffer_v  RN  r1
; count_v       RN  r4
; range_v       RN  r2
; value_v       RN  r3
; pos_v     RN  r5
; split     RN  r6
; bit           RN  lr
;int vp8dx_decode_bool_v6 ( BOOL_DECODER *br, int probability )
|vp8dx_decode_bool_v6| PROC
vp8dx_decode_bool_v6_internal
    stmdb   sp!, {r4 - r6, lr}

    ldr     r2, [br, #bool_decoder_range]
    ldr     r3, [br, #bool_decoder_value]

    mov     r6, r2, lsl #8
    sub     r6, r6, #256                ;   split = 1 +  (((range-1) * probability) >> 8)
    mov     r12, #1
    smlawb  r6, r6, prob, r12

    mov     lr, #0
    subs    r5, r3, r6, lsl #24

    ;cmp        r3, r1
    movhs   lr, #1
    movhs   r3, r5
    subhs   r2, r2, r6
    movlo   r2, r6

    cmp     r2, #0x80
    blt     range_less_0x80
    ;strd   r2, r3, [br, #bool_decoder_range]
    str     r2, [br, #bool_decoder_range]
    str     r3, [br, #bool_decoder_value]
    mov     r0, lr
    ldmia   sp!, {r4 - r6, pc}

range_less_0x80
    ldr     r5, [br, #bool_decoder_pos]
    ldr     r1, [br, #bool_decoder_buffer]
    ldr     r4, [br, #bool_decoder_count]
    add     r1, r1, r5

    clz       r12, r2
    sub       r12, r12, #24
    subs      r4, r4, r12
    ldrleb    r6, [r1], #1
    mov       r2, r2, lsl r12
    mov       r3, r3, lsl r12
    addle     r4, r4, #8
    rsble     r12, r4, #8
    addle     r5, r5, #1
    orrle     r3, r3, r6, lsl r12

    ;strd       r2, r3, [br, #bool_decoder_range]
    ;strd       r4, r5, [br, #bool_decoder_count]
    str         r2, [br, #bool_decoder_range]
    str         r3, [br, #bool_decoder_value]
    str         r4, [br, #bool_decoder_count]
    str         r5, [br, #bool_decoder_pos]

    mov     r0, lr

    ldmia   sp!, {r4 - r6, pc}
    ENDP    ; |vp8dx_decode_bool_v6|

    END
