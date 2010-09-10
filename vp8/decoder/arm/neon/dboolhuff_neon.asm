;
;  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
;
;  Use of this source code is governed by a BSD-style license
;  that can be found in the LICENSE file in the root of the source
;  tree. An additional intellectual property rights grant can be found
;  in the file PATENTS.  All contributing project authors may
;  be found in the AUTHORS file in the root of the source tree.
;


    EXPORT  |vp8_decode_value_neon|
    EXPORT  |vp8dx_start_decode_neon|
    EXPORT  |vp8dx_stop_decode_neon|
    EXPORT  |vp8dx_decode_bool_neon|

    ARM
    REQUIRE8
    PRESERVE8

    INCLUDE vpx_asm_offsets.asm

    AREA    |.text|, CODE, READONLY  ; name this block of code

;   int z = 0;
;   int bit;
;   for ( bit=bits-1; bit>=0; bit-- )
;   {
;       z |= (vp8dx_decode_bool(br, 0x80)<<bit);
;   }
;   return z;

;int vp8_decode_value_neon ( BOOL_DECODER *br, int bits )
|vp8_decode_value_neon| PROC
    stmdb   sp!, {r4 - r6, lr}
    mov     r4, r0
    mov     r5, r1
    mov     r6, #0

    subs    r5, r5, #1
    bmi     decode_value_exit

decode_value_loop
    mov     r1, #0x80
    mov     r0, r4
    bl      vp8dx_decode_bool_neon_internal       ; needed for conversion to s file
    orr     r6, r6, r0, lsl r5
    subs    r5, r5, #1
    bpl     decode_value_loop

decode_value_exit
    mov     r0, r6
    ldmia   sp!, {r4 - r6, pc}
    ENDP    ; |vp8_decode_value_neon|


;void vp8dx_start_decode_neon ( BOOL_DECODER *br, unsigned char *source )
|vp8dx_start_decode_neon| PROC
    stmdb   sp!, {r4 - r5, lr}
    mov     r2, #0
    mov     r3, #255

    str     r2, [r0, #bool_decoder_lowvalue]
    str     r3, [r0, #bool_decoder_range]
    str     r1, [r0, #bool_decoder_buffer]

    mov     r3, #8
    mov     r2, #4
    str     r3, [r0, #bool_decoder_count]
    str     r2, [r0, #bool_decoder_pos]

    ldrb    r2, [r1, #3]
    ldrb    r3, [r1, #2]
    ldrb    r4, [r1, #1]
    ldrb    r5, [r1]

    orr     r1, r2, r3, lsl #8
    orr     r1, r1, r4, lsl #16
    orr     r1, r1, r5, lsl #24

    str     r1, [r0, #bool_decoder_value]

    ldmia   sp!, {r4 - r5, pc}
    ENDP    ; |vp8dx_start_decode_neon|


;void vp8dx_stop_decode_neon ( BOOL_DECODER *bc );
|vp8dx_stop_decode_neon| PROC
    mov     pc, lr
    ENDP    ; |vp8dx_stop_decode_neon|


; bigsplit  RN  r1
; buffer_v  RN  r1
; count_v       RN  r4
; range_v       RN  r2
; value_v       RN  r3
; pos_v     RN  r5
; split     RN  r6
; bit           RN  lr
;int vp8dx_decode_bool_neon ( BOOL_DECODER *br, int probability )
|vp8dx_decode_bool_neon| PROC
vp8dx_decode_bool_neon_internal
;LDRD and STRD doubleword data transfers must be eight-byte aligned. Use ALIGN 8
;before memory allocation
    stmdb   sp!, {r4 - r5, lr}

    ldr     r2, [r0, #bool_decoder_range]       ;load range (r2), value(r3)
    ldr     r3, [r0, #bool_decoder_value]
    ;ldrd   r2, r3, [r0, #bool_decoder_range]   ;ldrd costs 2 cycles
    ;

    mov     r4, r2, lsl #8
    sub     r4, r4, #256
    mov     r12, #1

    smlawb  r4, r4, r1, r12         ;split = 1 +  (((range-1) * probability) >> 8)

    mov     lr, r0
    mov     r0, #0                  ;bit = 0
    ;
    subs    r5, r3, r4, lsl #24

    subhs   r2, r2, r4              ;range = br->range-split
    movlo   r2, r4                  ;range = split
    movhs   r0, #1                  ;bit = 1
    movhs   r3, r5                  ;value = value-bigsplit

    cmp     r2, #0x80
    blt     range_less_0x80
    strd    r2, r3, [lr, #bool_decoder_range]   ;store result

    ldmia   sp!, {r4 - r5, pc}

range_less_0x80

    ldrd    r4, r5, [lr, #bool_decoder_count]   ;load count, pos, buffer
    ldr     r1, [lr, #bool_decoder_buffer]

    clz     r12, r2
    add     r1, r1, r5

    sub     r12, r12, #24
    subs    r4, r4, r12             ;count -= shift
    mov     r2, r2, lsl r12         ;range <<= shift
    mov     r3, r3, lsl r12         ;value <<= shift
    addle   r4, r4, #8              ;count += 8
    ldrleb  r12, [r1], #1           ;br->buffer[br->pos]

    rsble   r1, r4, #8              ;-count
    addle   r5, r5, #1              ;br->pos++
    orrle   r3, r3, r12, lsl r1     ;value |= (br->buffer[br->pos]) << (-count)

    strd    r2, r3, [lr, #bool_decoder_range]   ;store result
    strd    r4, r5, [lr, #bool_decoder_count]

    ldmia   sp!, {r4 - r5, pc}
    ENDP    ; |vp8dx_decode_bool_neon|

    END
