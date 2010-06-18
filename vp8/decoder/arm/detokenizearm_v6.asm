;
;  Copyright (c) 2010 The VP8 project authors. All Rights Reserved.
;
;  Use of this source code is governed by a BSD-style license
;  that can be found in the LICENSE file in the root of the source
;  tree. An additional intellectual property rights grant can be found
;  in the file PATENTS.  All contributing project authors may
;  be found in the AUTHORS file in the root of the source tree.
;


    EXPORT  |vp8_decode_mb_tokens_v5|

    AREA    |.text|, CODE, READONLY  ; name this block of code

    INCLUDE vpx_asm_offsets.asm

l_qcoeff    EQU     0
l_i         EQU     4
l_type      EQU     8
l_stop      EQU     12
l_c         EQU     16
l_l_ptr      EQU     20
l_a_ptr      EQU     24
l_bc        EQU     28
l_coef_ptr   EQU     32
l_stacksize EQU     64


;; constant offsets -- these should be created at build time
c_onyxblock2left_offset      EQU 25
c_onyxblock2above_offset     EQU 50
c_entropy_nodes              EQU 11
c_dct_eob_token              EQU 11

|vp8_decode_mb_tokens_v5| PROC
    stmdb       sp!, {r4 - r11, lr}
    sub         sp, sp, #l_stacksize
    mov         r7, r1
    mov         r9, r0                      ;DETOK *detoken

    ldr         r1, [r9, #detok_current_bc]
    ldr         r0, [r9, #detok_qcoeff_start_ptr]
    mov         r11, #0
    mov         r3, #0x10

    cmp         r7, #1
    addeq       r11, r11, #24
    addeq       r3, r3, #8
    addeq       r0, r0, #3, 24

    str         r0, [sp, #l_qcoeff]
    str         r11, [sp, #l_i]
    str         r7, [sp, #l_type]
    str         r3, [sp, #l_stop]
    str         r1, [sp, #l_bc]

    add         lr, r9, r7, lsl #2

    ldr         r2, [r1, #bool_decoder_buffer]
    ldr         r3, [r1, #bool_decoder_pos]

    ldr         r10, [lr, #detok_coef_probs]
    ldr         r5, [r1, #bool_decoder_count]
    ldr         r6, [r1, #bool_decoder_range]
    ldr         r4, [r1, #bool_decoder_value]
    add         r8, r2, r3

    str         r10, [sp, #l_coef_ptr]


    ;align 4
BLOCK_LOOP
    ldr         r3, [r9, #detok_ptr_onyxblock2context_leftabove]
    ldr         r2, [r9, #DETOK_A]
    ldr         r1, [r9, #DETOK_L]
    ldrb        r12, [r3, +r11]                                 ; detoken->ptr_onyxblock2context_leftabove[i]

    cmp         r7, #0                                          ; check type
    moveq       r7, #1
    movne       r7, #0

    ldr         r0, [r2, +r12, lsl #2]                          ; a
    add         r1, r1, r12, lsl #4
    add         r3, r3, r11

    ldrb        r2, [r3, #c_onyxblock2above_offset]
    ldrb        r3, [r3, #c_onyxblock2left_offset]
    mov         lr, #c_entropy_nodes
;;  ;++

    ldr         r2, [r0, +r2, lsl #2]!
    add         r3, r1, r3, lsl #2
    str         r3, [sp, #l_l_ptr]
    ldr         r3, [r3]

    cmp         r2, #0
    movne       r2, #1
    cmp         r3, #0
    addne       r2, r2, #1

    str         r0, [sp, #l_a_ptr]
    smlabb      r0, r2, lr, r10
    mov         r1, #0                                          ; t = 0
    str         r7, [sp, #l_c]

    ;align 4
COEFF_LOOP
    ldr         r3, [r9, #detok_ptr_onyx_coef_bands_x]
    ldr         lr, [r9, #detok_onyx_coef_tree_ptr]

;;the following two lines are used if onyx_coef_bands_x is UINT16
;;  add         r3, r3, r7, lsl #1
;;  ldrh        r3, [r3]

;;the following line is used if onyx_coef_bands_x is UINT8
    ldrb        r3, [r7, +r3]


;;  ;++
;;  pld         [r8]
    ;++
    add         r0, r0, r3

    ;align 4
get_token_loop
    ldrb        r2, [r0, +r1, asr #1]
    mov         r3, r6, lsl #8
    sub         r3, r3, #256                    ;split = 1 +  (((range-1) * probability) >> 8)
    mov         r10, #1

    smlawb      r2, r3, r2, r10
    ldrb        r12, [r8]                       ;load cx data byte in stall slot
    ;++

    subs        r3, r4, r2, lsl #24             ;x = value-(split<<24)
    addhs       r1, r1, #1                      ;t += 1
    movhs       r4, r3                          ;update value
    subhs       r2, r6, r2                      ;range = range - split
    movlo       r6, r2

;;; ldrsbhs     r1, [r1, +lr]
    ldrsb     r1, [r1, +lr]


;; use branch for short pipelines ???
;;  cmp         r2, #0x80
;;  bcs         |$LN22@decode_mb_to|

    clz         r3, r2
    sub         r3, r3, #24
    subs        r5, r5, r3
    mov         r6, r2, lsl r3
    mov         r4, r4, lsl r3

;; use branch for short pipelines ???
;;  bgt         |$LN22@decode_mb_to|

    addle         r5, r5, #8
    rsble         r3, r5, #8
    addle         r8, r8, #1
    orrle         r4, r4, r12, lsl r3

;;|$LN22@decode_mb_to|

    cmp         r1, #0
    bgt         get_token_loop

    cmn         r1, #c_dct_eob_token             ;if(t == -DCT_EOB_TOKEN)
    beq         END_OF_BLOCK

    rsb         lr, r1, #0                      ;v = -t;

    cmp         lr, #4                          ;if(v > FOUR_TOKEN)
    ble         SKIP_EXTRABITS

    ldr         r3, [r9, #detok_teb_base_ptr]
    mov         r11, #1
    add         r7, r3, lr, lsl #4

    ldrsh       lr, [r7, #tokenextrabits_min_val];v = teb_ptr->min_val
    ldrsh       r0, [r7, #tokenextrabits_length];bits_count = teb_ptr->Length

extrabits_loop
    add         r3, r0, r7

    ldrb        r2, [r3, #4]
    mov         r3, r6, lsl #8
    sub         r3, r3, #256                    ;split = 1 +  (((range-1) * probability) >> 8)
    mov         r10, #1

    smlawb      r2, r3, r2, r10
    ldrb        r12, [r8]
    ;++

    subs        r10, r4, r2, lsl #24            ;x = value-(split<<24)
    movhs       r4, r10                         ;update value
    subhs       r2, r6, r2                      ;range = range - split
    addhs       lr, lr, r11, lsl r0             ;v += ((UINT16)1<<bits_count)
    movlo       r6, r2                          ;range = split


;; use branch for short pipelines ???
;;  cmp         r2, #0x80
;;  bcs         |$LN10@decode_mb_to|

    clz         r3, r2
    sub         r3, r3, #24
    subs        r5, r5, r3
    mov         r6, r2, lsl r3                  ;range
    mov         r4, r4, lsl r3                  ;value

    addle       r5, r5, #8
    addle       r8, r8, #1
    rsble       r3, r5, #8
    orrle       r4, r4, r12, lsl r3

;;|$LN10@decode_mb_to|
    subs         r0, r0, #1
    bpl         extrabits_loop


SKIP_EXTRABITS
    ldr         r11, [sp, #l_qcoeff]
    ldr         r0, [sp, #l_coef_ptr]

    cmp         r1, #0                          ;check for nonzero token
    beq         SKIP_EOB_CHECK              ;if t is zero, we will skip the eob table chec

    sub         r3, r6, #1                      ;range - 1
    ;++
    mov         r3, r3, lsl #7                  ; *= onyx_prob_half  (128)
    ;++
    mov         r3, r3, lsr #8
    add         r2, r3, #1                      ;split

    subs        r3, r4, r2, lsl #24             ;x = value-(split<<24)
    movhs       r4, r3                          ;update value
    subhs       r2, r6, r2                      ;range = range - split
    mvnhs       r3, lr
    addhs       lr, r3, #1                      ;v = (v ^ -1) + 1
    movlo       r6, r2                          ;range = split

;; use branch for short pipelines ???
;;  cmp         r2, #0x80
;;  bcs         |$LN6@decode_mb_to|

    clz         r3, r2
    sub         r3, r3, #24
    subs        r5, r5, r3
    mov         r6, r2, lsl r3
    mov         r4, r4, lsl r3
    ldrleb      r2, [r8], #1
    addle       r5, r5, #8
    rsble       r3, r5, #8
    orrle       r4, r4, r2, lsl r3

;;|$LN6@decode_mb_to|
    add         r0, r0, #0xB

    cmn         r1, #1

    addlt       r0, r0, #0xB

    mvn         r1, #1

SKIP_EOB_CHECK
    ldr         r7, [sp, #l_c]
    ldr         r3, [r9, #detok_scan]
    add         r1, r1, #2
    cmp         r7, #(0x10 - 1)                     ;assume one less for now.... increment below

    ldr         r3, [r3, +r7, lsl #2]
    add         r7, r7, #1
    add         r3, r11, r3, lsl #1

    str         r7, [sp, #l_c]
    strh        lr, [r3]

    blt         COEFF_LOOP

    sub         r7, r7, #1                          ;if(t != -DCT_EOB_TOKEN) --c

END_OF_BLOCK
    ldr         r3, [sp, #l_type]
    ldr         r10, [sp, #l_coef_ptr]
    ldr         r0, [sp, #l_qcoeff]
    ldr         r11, [sp, #l_i]
    ldr         r12, [sp, #l_stop]

    cmp         r3, #0
    moveq       r1, #1
    movne       r1, #0
    add         r3, r11, r9

    cmp         r7, r1
    strb        r7, [r3, #detok_eob]

    ldr         r7, [sp, #l_l_ptr]
    ldr         r2, [sp, #l_a_ptr]
    movne       r3, #1
    moveq       r3, #0

    add         r0, r0, #0x20
    add         r11, r11, #1
    str         r3, [r7]
    str         r3, [r2]
    str         r0, [sp, #l_qcoeff]
    str         r11, [sp, #l_i]

    cmp         r11, r12                            ;i >= stop ?
    ldr         r7, [sp, #l_type]
    mov         lr, #0xB

    blt         BLOCK_LOOP

    cmp         r11, #0x19
    bne         ln2_decode_mb_to

    ldr         r12, [r9, #detok_qcoeff_start_ptr]
    ldr         r10, [r9, #detok_coef_probs]
    mov         r7, #0
    mov         r3, #0x10
    str         r12, [sp, #l_qcoeff]
    str         r7, [sp, #l_i]
    str         r7, [sp, #l_type]
    str         r3, [sp, #l_stop]

    str         r10, [sp, #l_coef_ptr]

    b           BLOCK_LOOP

ln2_decode_mb_to
    cmp         r11, #0x10
    bne         ln1_decode_mb_to

    ldr         r10, [r9, #0x30]

    mov         r7, #2
    mov         r3, #0x18

    str         r7, [sp, #l_type]
    str         r3, [sp, #l_stop]

    str         r10, [sp, #l_coef_ptr]
    b           BLOCK_LOOP

ln1_decode_mb_to
    ldr         r2, [sp, #l_bc]
    mov         r0, #0
    nop

    ldr         r3, [r2, #bool_decoder_buffer]
    str         r5, [r2, #bool_decoder_count]
    str         r4, [r2, #bool_decoder_value]
    sub         r3, r8, r3
    str         r3, [r2, #bool_decoder_pos]
    str         r6, [r2, #bool_decoder_range]

    add         sp, sp, #l_stacksize
    ldmia       sp!, {r4 - r11, pc}

    ENDP  ; |vp8_decode_mb_tokens_v5|

    END
