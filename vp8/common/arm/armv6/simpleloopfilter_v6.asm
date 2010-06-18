;
;  Copyright (c) 2010 The VP8 project authors. All Rights Reserved.
;
;  Use of this source code is governed by a BSD-style license
;  that can be found in the LICENSE file in the root of the source
;  tree. An additional intellectual property rights grant can be found
;  in the file PATENTS.  All contributing project authors may
;  be found in the AUTHORS file in the root of the source tree.
;


    EXPORT |vp8_loop_filter_simple_horizontal_edge_armv6|
    EXPORT |vp8_loop_filter_simple_vertical_edge_armv6|

    AREA    |.text|, CODE, READONLY  ; name this block of code

    MACRO
    TRANSPOSE_MATRIX $a0, $a1, $a2, $a3, $b0, $b1, $b2, $b3
    ; input: $a0, $a1, $a2, $a3; output: $b0, $b1, $b2, $b3
    ; a0: 03 02 01 00
    ; a1: 13 12 11 10
    ; a2: 23 22 21 20
    ; a3: 33 32 31 30
    ;     b3 b2 b1 b0

    uxtb16      $b1, $a1                    ; xx 12 xx 10
    uxtb16      $b0, $a0                    ; xx 02 xx 00
    uxtb16      $b3, $a3                    ; xx 32 xx 30
    uxtb16      $b2, $a2                    ; xx 22 xx 20
    orr         $b1, $b0, $b1, lsl #8       ; 12 02 10 00
    orr         $b3, $b2, $b3, lsl #8       ; 32 22 30 20

    uxtb16      $a1, $a1, ror #8            ; xx 13 xx 11
    uxtb16      $a3, $a3, ror #8            ; xx 33 xx 31
    uxtb16      $a0, $a0, ror #8            ; xx 03 xx 01
    uxtb16      $a2, $a2, ror #8            ; xx 23 xx 21
    orr         $a0, $a0, $a1, lsl #8       ; 13 03 11 01
    orr         $a2, $a2, $a3, lsl #8       ; 33 23 31 21

    pkhtb       $b2, $b3, $b1, asr #16      ; 32 22 12 02   -- p1
    pkhbt       $b0, $b1, $b3, lsl #16      ; 30 20 10 00   -- p3

    pkhtb       $b3, $a2, $a0, asr #16      ; 33 23 13 03   -- p0
    pkhbt       $b1, $a0, $a2, lsl #16      ; 31 21 11 01   -- p2
    MEND


src         RN  r0
pstep       RN  r1

;r0     unsigned char *src_ptr,
;r1     int src_pixel_step,
;r2     const char *flimit,
;r3     const char *limit,
;stack  const char *thresh,
;stack  int  count

;Note: All 16 elements in flimit are equal. So, in the code, only one load is needed
;for flimit. Same way applies to limit and thresh.

;-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
|vp8_loop_filter_simple_horizontal_edge_armv6| PROC
;-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    stmdb       sp!, {r4 - r11, lr}

    sub         src, src, pstep, lsl #1     ; move src pointer down by 2 lines

    ldr         r12, [r3], #4               ; limit
    ldr         r3, [src], pstep            ; p1

    ldr         r9, [sp, #36]               ; count for 8-in-parallel
    ldr         r4, [src], pstep            ; p0

    ldr         r7, [r2], #4                ; flimit
    ldr         r5, [src], pstep            ; q0
    ldr         r2, c0x80808080

    ldr         r6, [src]                   ; q1

    uadd8       r7, r7, r7                  ; flimit * 2
    mov         r9, r9, lsl #1              ; 4-in-parallel
    uadd8       r12, r7, r12                ; flimit * 2 + limit

|simple_hnext8|
    ; vp8_simple_filter_mask() function

    uqsub8      r7, r3, r6                  ; p1 - q1
    uqsub8      r8, r6, r3                  ; q1 - p1
    uqsub8      r10, r4, r5                 ; p0 - q0
    uqsub8      r11, r5, r4                 ; q0 - p0
    orr         r8, r8, r7                  ; abs(p1 - q1)
    ldr         lr, c0x7F7F7F7F             ; 01111111 mask
    orr         r10, r10, r11               ; abs(p0 - q0)
    and         r8, lr, r8, lsr #1          ; abs(p1 - q1) / 2
    uqadd8      r10, r10, r10               ; abs(p0 - q0) * 2
    mvn         lr, #0                      ; r10 == -1
    uqadd8      r10, r10, r8                ; abs(p0 - q0)*2 + abs(p1 - q1)/2
    ; STALL waiting on r10 :(
    uqsub8      r10, r10, r12               ; compare to flimit
    mov         r8, #0

    usub8       r10, r8, r10                ; use usub8 instead of ssub8
    ; STALL (maybe?) when are flags set? :/
    sel         r10, lr, r8                 ; filter mask: lr

    cmp         r10, #0
    beq         simple_hskip_filter         ; skip filtering

    ;vp8_simple_filter() function

    eor         r3, r3, r2                  ; p1 offset to convert to a signed value
    eor         r6, r6, r2                  ; q1 offset to convert to a signed value
    eor         r4, r4, r2                  ; p0 offset to convert to a signed value
    eor         r5, r5, r2                  ; q0 offset to convert to a signed value

    qsub8       r3, r3, r6                  ; vp8_filter (r3) = vp8_signed_char_clamp(p1-q1)
    qsub8       r6, r5, r4                  ; vp8_filter = vp8_signed_char_clamp(vp8_filter + 3 * ( q0 - p0))

    qadd8       r3, r3, r6
    ldr         r8, c0x03030303             ; r8 = 3

    qadd8       r3, r3, r6
    ldr         r7, c0x04040404

    qadd8       r3, r3, r6
    and         r3, r3, lr                  ; vp8_filter &= mask;

    ;save bottom 3 bits so that we round one side +4 and the other +3
    qadd8       r8 , r3 , r8                ; Filter2 (r8) = vp8_signed_char_clamp(vp8_filter+3)
    qadd8       r3 , r3 , r7                ; Filter1 (r3) = vp8_signed_char_clamp(vp8_filter+4)

    mov         r7, #0
    shadd8      r8 , r8 , r7                ; Filter2 >>= 3
    shadd8      r3 , r3 , r7                ; Filter1 >>= 3
    shadd8      r8 , r8 , r7
    shadd8      r3 , r3 , r7
    shadd8      r8 , r8 , r7                ; r8: Filter2
    shadd8      r3 , r3 , r7                ; r7: filter1

    ;calculate output
    sub         src, src, pstep, lsl #1

    qadd8       r4, r4, r8                  ; u = vp8_signed_char_clamp(p0 + Filter2)
    qsub8       r5 ,r5, r3                  ; u = vp8_signed_char_clamp(q0 - Filter1)
    eor         r4, r4, r2                  ; *op0 = u^0x80
    str         r4, [src], pstep            ; store op0 result
    eor         r5, r5, r2                  ; *oq0 = u^0x80
    str         r5, [src], pstep            ; store oq0 result

|simple_hskip_filter|
    add         src, src, #4
    sub         src, src, pstep
    sub         src, src, pstep, lsl #1

    subs        r9, r9, #1

    ;pld            [src]
    ;pld            [src, pstep]
    ;pld            [src, pstep, lsl #1]

    ldrne           r3, [src], pstep            ; p1
    ldrne           r4, [src], pstep            ; p0
    ldrne           r5, [src], pstep            ; q0
    ldrne           r6, [src]                   ; q1

    bne         simple_hnext8

    ldmia       sp!, {r4 - r11, pc}
    ENDP        ; |vp8_loop_filter_simple_horizontal_edge_armv6|


;-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
|vp8_loop_filter_simple_vertical_edge_armv6| PROC
;-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    stmdb       sp!, {r4 - r11, lr}

    ldr         r12, [r2], #4               ; r12: flimit
    ldr         r2, c0x80808080
    ldr         r7, [r3], #4                ; limit

    ; load soure data to r7, r8, r9, r10
    ldrh        r3, [src, #-2]
    ldrh        r4, [src], pstep
    uadd8       r12, r12, r12               ; flimit * 2

    ldrh        r5, [src, #-2]
    ldrh        r6, [src], pstep
    uadd8       r12, r12, r7                ; flimit * 2 + limit

    pkhbt       r7, r3, r4, lsl #16

    ldrh        r3, [src, #-2]
    ldrh        r4, [src], pstep
    ldr         r11, [sp, #40]              ; count (r11) for 8-in-parallel

    pkhbt       r8, r5, r6, lsl #16

    ldrh        r5, [src, #-2]
    ldrh        r6, [src], pstep
    mov         r11, r11, lsl #1            ; 4-in-parallel

|simple_vnext8|
    ; vp8_simple_filter_mask() function
    pkhbt       r9, r3, r4, lsl #16
    pkhbt       r10, r5, r6, lsl #16

    ;transpose r7, r8, r9, r10 to r3, r4, r5, r6
    TRANSPOSE_MATRIX r7, r8, r9, r10, r3, r4, r5, r6

    uqsub8      r7, r3, r6                  ; p1 - q1
    uqsub8      r8, r6, r3                  ; q1 - p1
    uqsub8      r9, r4, r5                  ; p0 - q0
    uqsub8      r10, r5, r4                 ; q0 - p0
    orr         r7, r7, r8                  ; abs(p1 - q1)
    orr         r9, r9, r10                 ; abs(p0 - q0)
    ldr         lr, c0x7F7F7F7F             ; 0111 1111 mask
    uqadd8      r9, r9, r9                  ; abs(p0 - q0) * 2
    and         r7, lr, r7, lsr #1          ; abs(p1 - q1) / 2
    mov         r8, #0
    uqadd8      r7, r7, r9                  ; abs(p0 - q0)*2 + abs(p1 - q1)/2
    mvn         r10, #0                     ; r10 == -1
    uqsub8      r7, r7, r12                 ; compare to flimit

    usub8       r7, r8, r7
    sel         r7, r10, r8                 ; filter mask: lr

    cmp         lr, #0
    beq         simple_vskip_filter         ; skip filtering

    ;vp8_simple_filter() function
    eor         r3, r3, r2                  ; p1 offset to convert to a signed value
    eor         r6, r6, r2                  ; q1 offset to convert to a signed value
    eor         r4, r4, r2                  ; p0 offset to convert to a signed value
    eor         r5, r5, r2                  ; q0 offset to convert to a signed value

    qsub8       r3, r3, r6                  ; vp8_filter (r3) = vp8_signed_char_clamp(p1-q1)
    qsub8       r6, r5, r4                  ; vp8_filter = vp8_signed_char_clamp(vp8_filter + 3 * ( q0 - p0))

    qadd8       r3, r3, r6
    ldr         r8, c0x03030303             ; r8 = 3

    qadd8       r3, r3, r6
    ldr         r7, c0x04040404

    qadd8       r3, r3, r6
    and         r3, r3, lr                  ; vp8_filter &= mask

    ;save bottom 3 bits so that we round one side +4 and the other +3
    qadd8       r8 , r3 , r8                ; Filter2 (r8) = vp8_signed_char_clamp(vp8_filter+3)
    qadd8       r3 , r3 , r7                ; Filter1 (r3) = vp8_signed_char_clamp(vp8_filter+4)

    mov         r7, #0
    shadd8      r8 , r8 , r7                ; Filter2 >>= 3
    shadd8      r3 , r3 , r7                ; Filter1 >>= 3
    shadd8      r8 , r8 , r7
    shadd8      r3 , r3 , r7
    shadd8      r8 , r8 , r7                ; r8: filter2
    shadd8      r3 , r3 , r7                ; r7: filter1

    ;calculate output
    sub         src, src, pstep, lsl #2

    qadd8       r4, r4, r8                  ; u = vp8_signed_char_clamp(p0 + Filter2)
    qsub8       r5, r5, r3                  ; u = vp8_signed_char_clamp(q0 - Filter1)
    eor         r4, r4, r2                  ; *op0 = u^0x80
    eor         r5, r5, r2                  ; *oq0 = u^0x80

    strb        r4, [src, #-1]              ; store the result
    mov         r4, r4, lsr #8
    strb        r5, [src], pstep
    mov         r5, r5, lsr #8

    strb        r4, [src, #-1]
    mov         r4, r4, lsr #8
    strb        r5, [src], pstep
    mov         r5, r5, lsr #8

    strb        r4, [src, #-1]
    mov         r4, r4, lsr #8
    strb        r5, [src], pstep
    mov         r5, r5, lsr #8

    strb        r4, [src, #-1]
    strb        r5, [src], pstep

|simple_vskip_filter|
    subs        r11, r11, #1

    ;pld            [src]
    ;pld            [src, pstep]
    ;pld            [src, pstep, lsl #1]

    ; load soure data to r7, r8, r9, r10
    ldrneh      r3, [src, #-2]
    ldrneh      r4, [src], pstep

    ldrneh      r5, [src, #-2]
    ldrneh      r6, [src], pstep

    pkhbt       r7, r3, r4, lsl #16

    ldrneh      r3, [src, #-2]
    ldrneh      r4, [src], pstep

    pkhbt       r8, r5, r6, lsl #16

    ldrneh      r5, [src, #-2]
    ldrneh      r6, [src], pstep

    bne         simple_vnext8

    ldmia       sp!, {r4 - r12, pc}
    ENDP        ; |vp8_loop_filter_simple_vertical_edge_armv6|

; Constant Pool
c0x80808080 DCD     0x80808080
c0x03030303 DCD     0x03030303
c0x04040404 DCD     0x04040404
c0x01010101 DCD     0x01010101
c0x7F7F7F7F DCD     0x7F7F7F7F

    END
