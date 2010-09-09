;
;  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
;
;  Use of this source code is governed by a BSD-style license
;  that can be found in the LICENSE file in the root of the source
;  tree. An additional intellectual property rights grant can be found
;  in the file PATENTS.  All contributing project authors may
;  be found in the AUTHORS file in the root of the source tree.
;


    EXPORT  |horizontal_line_4_5_scale_armv4|
    EXPORT  |vertical_band_4_5_scale_armv4|
    EXPORT  |horizontal_line_2_3_scale_armv4|
    EXPORT  |vertical_band_2_3_scale_armv4|
    EXPORT  |horizontal_line_3_5_scale_armv4|
    EXPORT  |vertical_band_3_5_scale_armv4|
    EXPORT  |horizontal_line_3_4_scale_armv4|
    EXPORT  |vertical_band_3_4_scale_armv4|
    EXPORT  |horizontal_line_1_2_scale_armv4|
    EXPORT  |vertical_band_1_2_scale_armv4|

    AREA    |.text|, CODE, READONLY  ; name this block of code

src         RN  r0
srcw        RN  r1
dest        RN  r2
mask        RN  r12
c51_205     RN  r10
c102_154    RN  r11
;/****************************************************************************
; *
; *  ROUTINE       : horizontal_line_4_5_scale_armv4
; *
; *  INPUTS        : const unsigned char *source : Pointer to source data.
; *                  unsigned int source_width    : Stride of source.
; *                  unsigned char *dest         : Pointer to destination data.
; *                  unsigned int dest_width      : Stride of destination (NOT USED).
; *
; *  OUTPUTS       : None.
; *
; *  RETURNS       : void
; *
; *  FUNCTION      : Copies horizontal line of pixels from source to
; *                  destination scaling up by 4 to 5.
; *
; *  SPECIAL NOTES : None.
; *
; ****************************************************************************/
;void horizontal_line_4_5_scale_armv4
;(
;   r0 = UINT8 *source
;   r1 = UINT32 source_width
;   r2 = UINT8 *dest
;   r3 = UINT32 dest_width
;)
|horizontal_line_4_5_scale_armv4| PROC
    stmdb   sp!, {r4 - r11, lr}

    mov     mask, #255              ; mask for selection
    ldr     c51_205, =0x3300cd
    ldr     c102_154, =0x66009a

    ldr     r3, [src], #4

hl45_loop

    and     r4, r3, mask            ; a = src[0]
    and     r5, mask, r3, lsr #8    ; b = src[1]
    strb    r4, [dest], #1

    orr     r6, r4, r5, lsl #16     ; b | a
    and     r7, mask, r3, lsr #16   ; c = src[2]
    mul     r6, c51_205, r6         ; a * 51 + 205 * b

    orr     r5, r5, r7, lsl #16     ; c | b
    mul     r5, c102_154, r5        ; b * 102 + 154 * c
    add     r6, r6, #0x8000
    and     r8, mask, r3, lsr #24   ; d = src[3]
    mov     r6, r6, lsr #24
    strb    r6, [dest], #1

    orr     r7, r8, r7, lsl #16     ; c | d
    mul     r7, c102_154, r7        ; c * 154 + 102 * d
    add     r5, r5, #0x8000
    ldr     r3, [src], #4
    mov     r5, r5, lsr #24
    strb    r5, [dest], #1

    add     r7, r7, #0x8000
    and     r9, mask, r3            ; e = src[4]
    orr     r9, r9, r8, lsl #16     ; d | e
    mul     r9, c51_205, r9         ; d * 205 + 51 * e
    mov     r7, r7, lsr #24
    strb    r7, [dest], #1

    add     r9, r9, #0x8000
    subs    srcw, srcw, #4
    mov     r9, r9, lsr #24
    strb    r9, [dest], #1

    bne     hl45_loop

    and     r4, r3, mask
    and     r5, mask, r3, lsl #8
    strb    r4, [dest], #1

    orr     r6, r4, r5, lsl #16     ; b | a
    mul     r6, c51_205, r6

    and     r7, mask, r3, lsl #16
    orr     r5, r5, r7, lsl #16     ; c | b
    mul     r5, c102_154, r5
    add     r6, r6, #0x8000
    and     r8, mask, r3, lsl #24
    mov     r6, r6, lsr #24
    strb    r6, [dest], #1

    orr     r7, r8, r7, lsl #16     ; c | d
    mul     r7, c102_154, r7
    add     r5, r5, #0x8000
    mov     r5, r5, lsr #24
    strb    r5, [dest], #1

    add     r7, r7, #0x8000
    mov     r7, r7, lsr #24
    strb    r7, [dest], #1

    ldrb    r3, [src]
    strb    r3, [dest], #1

    ldmia   sp!, {r4 - r11, pc}
    ENDP    ;|vp8cx_horizontal_line_4_5_scale_c|

;/****************************************************************************
; *
; *  ROUTINE       : vertical_band_4_5_scale_armv4
; *
; *  INPUTS        : unsigned char *dest    : Pointer to destination data.
; *                  unsigned int dest_pitch : Stride of destination data.
; *                  unsigned int dest_width : Width of destination data.
; *
; *  OUTPUTS       : None.
; *
; *  RETURNS       : void
; *
; *  FUNCTION      : Scales vertical band of pixels by scale 4 to 5. The
; *                  height of the band scaled is 4-pixels.
; *
; *  SPECIAL NOTES : The routine uses the first line of the band below
; *                  the current band.
; *
; ****************************************************************************/
;void vertical_band_4_5_scale_armv4
;(
;   r0 = UINT8 *dest
;   r1 = UINT32 dest_pitch
;   r2 = UINT32 dest_width
;)
|vertical_band_4_5_scale_armv4| PROC
    stmdb   sp!, {r4 - r11, lr}

    ldr     c51_205, =0x3300cd
    ldr     c102_154, =0x66009a

vl45_loop
    mov     r3, src
    ldrb    r4, [r3], r1            ; a = des [0]
    ldrb    r5, [r3], r1            ; b = des [dest_pitch]
    ldrb    r7, [r3], r1            ; c = des[dest_pitch*2]
    add     lr, src, r1

    orr     r6, r4, r5, lsl #16     ; b | a
    mul     r6, c51_205, r6         ; a * 51 + 205 * b

    ldrb    r8, [r3], r1            ; d = des[dest_pitch*3]
    orr     r5, r5, r7, lsl #16     ; c | b
    mul     r5, c102_154, r5        ; b * 102 + 154 * c
    add     r6, r6, #0x8000
    orr     r7, r8, r7, lsl #16     ; c | d
    mov     r6, r6, lsr #24
    strb    r6, [lr], r1

    ldrb    r9, [r3, r1]            ; e = des [dest_pitch * 5]
    mul     r7, c102_154, r7        ; c * 154 + 102 * d
    add     r5, r5, #0x8000
    orr     r9, r9, r8, lsl #16     ; d | e
    mov     r5, r5, lsr #24
    strb    r5, [lr], r1

    mul     r9, c51_205, r9         ; d * 205 + 51 * e
    add     r7, r7, #0x8000
    add     src, src, #1
    mov     r7, r7, lsr #24
    strb    r7, [lr], r1

    add     r9, r9, #0x8000
    subs    r2, r2, #1
    mov     r9, r9, lsr #24
    strb    r9, [lr], r1

    bne     vl45_loop

    ldmia   sp!, {r4 - r11, pc}
    ENDP    ;|vertical_band_4_5_scale_armv4|

;/****************************************************************************
; *
; *  ROUTINE       : horizontal_line_2_3_scale_armv4
; *
; *  INPUTS        : const unsigned char *source : Pointer to source data.
; *                  unsigned int source_width    : Stride of source.
; *                  unsigned char *dest         : Pointer to destination data.
; *                  unsigned int dest_width      : Stride of destination (NOT USED).
; *
; *  OUTPUTS       : None.
; *
; *  RETURNS       : void
; *
; *  FUNCTION      : Copies horizontal line of pixels from source to
; *                  destination scaling up by 2 to 3.
; *
; *  SPECIAL NOTES : None.
; *
; *
; ****************************************************************************/
;void horizontal_line_2_3_scale_armv4
;(
;   const unsigned char *source,
;   unsigned int source_width,
;   unsigned char *dest,
;   unsigned int dest_width
;)
|horizontal_line_2_3_scale_armv4| PROC
    stmdb   sp!, {r4 - r11, lr}
    ldr     lr,  =85
    ldr     r12, =171

hl23_loop

    ldrb    r3, [src], #1           ; a
    ldrb    r4, [src], #1           ; b
    ldrb    r5, [src]               ; c

    strb    r3, [dest], #1
    mul     r4, r12, r4             ; b * 171
    mla     r6, lr, r3, r4          ; a * 85
    mla     r7, lr, r5, r4          ; c * 85

    add     r6, r6, #128
    mov     r6, r6, lsr #8
    strb    r6, [dest], #1

    add     r7, r7, #128
    mov     r7, r7, lsr #8
    strb    r7, [dest], #1

    subs    srcw, srcw, #2
    bne     hl23_loop

    ldrb    r4, [src, #1]           ; b
    strb    r5, [dest], #1
    strb    r4, [dest, #1]

    mul     r4, r12, r4             ; b * 171
    mla     r6, lr, r5, r4          ; a * 85 + b *171

    add     r6, r6, #128
    mov     r6, r6, lsr #8
    strb    r6, [dest]

    ldmia   sp!, {r4 - r11, pc}
    ENDP    ;|horizontal_line_2_3_scale_armv4|

;/****************************************************************************
; *
; *  ROUTINE       : vertical_band_2_3_scale_armv4
; *
; *  INPUTS        : unsigned char *dest    : Pointer to destination data.
; *                  unsigned int dest_pitch : Stride of destination data.
; *                  unsigned int dest_width : Width of destination data.
; *
; *  OUTPUTS       : None.
; *
; *  RETURNS       : void
; *
; *  FUNCTION      : Scales vertical band of pixels by scale 2 to 3. The
; *                  height of the band scaled is 2-pixels.
; *
; *  SPECIAL NOTES : The routine uses the first line of the band below
; *                  the current band.
; *
; ****************************************************************************/
;void vertical_band_2_3_scale_armv4
;(
;   r0 = UINT8 *dest
;   r1 = UINT32 dest_pitch
;   r2 = UINT32 dest_width
;)
|vertical_band_2_3_scale_armv4| PROC
    stmdb   sp!, {r4 - r8, lr}
    ldr     lr,  =85
    ldr     r12, =171
    add     r3, r1, r1, lsl #1      ; 3 * dest_pitch

vl23_loop
    ldrb    r4, [src]               ; a = des [0]
    ldrb    r5, [src, r1]           ; b = des [dest_pitch]
    ldrb    r7, [src, r3]           ; c = des [dest_pitch*3]
    subs    r2, r2, #1

    mul     r5, r12, r5             ; b * 171
    mla     r6, lr, r4, r5          ; a * 85
    mla     r8, lr, r7, r5          ; c * 85

    add     r6, r6, #128
    mov     r6, r6, lsr #8
    strb    r6, [src, r1]

    add     r8, r8, #128
    mov     r8, r8, lsr #8
    strb    r8, [src, r1, lsl #1]

    add     src, src, #1

    bne     vl23_loop

    ldmia   sp!, {r4 - r8, pc}
    ENDP    ;|vertical_band_2_3_scale_armv4|

;/****************************************************************************
; *
; *  ROUTINE       : vp8cx_horizontal_line_3_5_scale_c
; *
; *  INPUTS        : const unsigned char *source : Pointer to source data.
; *                  unsigned int source_width    : Stride of source.
; *                  unsigned char *dest         : Pointer to destination data.
; *                  unsigned int dest_width      : Stride of destination (NOT USED).
; *
; *  OUTPUTS       : None.
; *
; *  RETURNS       : void
; *
; *  FUNCTION      : Copies horizontal line of pixels from source to
; *                  destination scaling up by 3 to 5.
; *
; *  SPECIAL NOTES : None.
; *
; *
; ****************************************************************************/
;void vp8cx_horizontal_line_3_5_scale_c
;(
;   const unsigned char *source,
;   unsigned int source_width,
;   unsigned char *dest,
;   unsigned int dest_width
;)
|horizontal_line_3_5_scale_armv4| PROC
    stmdb   sp!, {r4 - r11, lr}

    ldr     c51_205, =0x3300cd
    ldr     c102_154, =0x66009a

    ldrb    r4, [src], #1           ; a = src[0]

hl35_loop

    ldrb    r8, [src], #1           ; b = src[1]
    strb    r4, [dest], #1

    orr     r6, r4, r8, lsl #16     ; b | a
    ldrb    r9, [src], #1           ; c = src[2]
    mul     r6, c102_154, r6        ; a * 102 + 154 * b

    orr     r5, r9, r8, lsl #16     ; b | c
    mul     r5, c51_205, r5         ; b * 205 + 51 * c
    add     r6, r6, #0x8000
    ldrb    r4, [src], #1           ; d = src[3]
    mov     r6, r6, lsr #24
    strb    r6, [dest], #1

    orr     r7, r8, r9, lsl #16     ; c | b
    mul     r7, c51_205, r7         ; c * 205 + 154 * b
    add     r5, r5, #0x8000
    mov     r5, r5, lsr #24
    strb    r5, [dest], #1

    orr     r9, r4, r9, lsl #16     ; c | d
    mul     r9, c102_154, r9        ; c * 154 + 102 * d
    add     r7, r7, #0x8000
    mov     r7, r7, lsr #24
    strb    r7, [dest], #1

    add     r9, r9, #0x8000
    subs    srcw, srcw, #3
    mov     r9, r9, lsr #24
    strb    r9, [dest], #1

    bpl     hl35_loop

    ldrb    r5, [src], #1           ; b = src[1]
    strb    r4, [dest], #1

    orr     r6, r4, r8, lsl #16     ; b | a
    ldrb    r9, [src], #1           ; c = src[2]
    mul     r6, c102_154, r6        ; a * 102 + 154 * b

    orr     r5, r9, r8, lsl #16     ; b | c
    mul     r5, c51_205, r5         ; b * 205 + 51 * c
    add     r6, r6, #0x8000
    mov     r6, r6, lsr #24
    strb    r6, [dest], #1

    orr     r7, r8, r9, lsl #16     ; c | b
    mul     r7, c51_205, r7         ; c * 205 + 154 * b
    add     r5, r5, #0x8000
    mov     r5, r5, lsr #24
    strb    r5, [dest], #1

    add     r7, r7, #0x8000
    mov     r7, r7, lsr #24
    strb    r7, [dest], #1
    strb    r9, [dest], #1

    ldmia   sp!, {r4 - r11, pc}
    ENDP    ;|vp8cx_horizontal_line_3_5_scale_c|


;/****************************************************************************
; *
; *  ROUTINE       : vp8cx_vertical_band_3_5_scale_c
; *
; *  INPUTS        : unsigned char *dest    : Pointer to destination data.
; *                  unsigned int dest_pitch : Stride of destination data.
; *                  unsigned int dest_width : Width of destination data.
; *
; *  OUTPUTS       : None.
; *
; *  RETURNS       : void
; *
; *  FUNCTION      : Scales vertical band of pixels by scale 3 to 5. The
; *                  height of the band scaled is 3-pixels.
; *
; *  SPECIAL NOTES : The routine uses the first line of the band below
; *                  the current band.
; *
; ****************************************************************************/
;void vertical_band_4_5_scale_armv4
;(
;   r0 = UINT8 *dest
;   r1 = UINT32 dest_pitch
;   r2 = UINT32 dest_width
;)
|vertical_band_3_5_scale_armv4| PROC
    stmdb   sp!, {r4 - r11, lr}

    ldr     c51_205, =0x3300cd
    ldr     c102_154, =0x66009a

vl35_loop
    mov     r3, src
    ldrb    r4, [r3], r1            ; a = des [0]
    ldrb    r5, [r3], r1            ; b = des [dest_pitch]
    ldrb    r7, [r3], r1            ; c = des[dest_pitch*2]
    add     lr, src, r1

    orr     r8, r4, r5, lsl #16     ; b | a
    mul     r6, c102_154, r8        ; a * 102 + 154 * b

    ldrb    r8, [r3, r1, lsl #1]    ; d = des[dest_pitch*5]
    orr     r3, r7, r5, lsl #16     ; b | c
    mul     r9, c51_205, r3         ; b * 205 + 51 * c
    add     r6, r6, #0x8000
    orr     r3, r5, r7, lsl #16     ; c | b
    mov     r6, r6, lsr #24
    strb    r6, [lr], r1

    mul     r5, c51_205, r3         ; c * 205 + 154 * b
    add     r9, r9, #0x8000
    orr     r3, r8, r7, lsl #16     ; c | d
    mov     r9, r9, lsr #24
    strb    r9, [lr], r1

    mul     r7, c102_154, r3        ; c * 154 + 102 * d
    add     r5, r5, #0x8000
    add     src, src, #1
    mov     r5, r5, lsr #24
    strb    r5, [lr], r1

    add     r7, r7, #0x8000
    subs    r2, r2, #1
    mov     r7, r7, lsr #24
    strb    r7, [lr], r1


    bne     vl35_loop

    ldmia   sp!, {r4 - r11, pc}
    ENDP    ;|vertical_band_3_5_scale_armv4|

;/****************************************************************************
; *
; *  ROUTINE       : horizontal_line_3_4_scale_armv4
; *
; *  INPUTS        : const unsigned char *source : Pointer to source data.
; *                  unsigned int source_width    : Stride of source.
; *                  unsigned char *dest         : Pointer to destination data.
; *                  unsigned int dest_width      : Stride of destination (NOT USED).
; *
; *  OUTPUTS       : None.
; *
; *  RETURNS       : void
; *
; *  FUNCTION      : Copies horizontal line of pixels from source to
; *                  destination scaling up by 3 to 4.
; *
; *  SPECIAL NOTES : None.
; *
; *
; ****************************************************************************/
;void horizontal_line_3_4_scale_armv4
;(
;   const unsigned char *source,
;   unsigned int source_width,
;   unsigned char *dest,
;   unsigned int dest_width
;)
|horizontal_line_3_4_scale_armv4| PROC
    stmdb   sp!, {r4 - r11, lr}

    ldr     r10, =64
    ldr     r11, =192
    mov     r9, #128

    ldrb    r4, [src], #1           ; a = src[0]

hl34_loop

    ldrb    r8, [src], #1           ; b = src[1]
    ldrb    r7, [src], #1           ; c = src[2]
    strb    r4, [dest], #1

    mla     r4, r10, r4, r9         ; a*64 + 128
    mla     r4, r11, r8, r4         ; a*64 + b*192 + 1

    add     r8, r8, #1              ; b + 1
    add     r8, r8, r7              ; b + c + 1
    mov     r8, r8, asr #1          ; (b + c + 1) >> 1

    mov     r4, r4, asr #8          ; (a*64 + b*192 + 1) >> 8
    strb    r4, [dest], #1

    strb    r8, [dest], #1

    ldrb    r4, [src], #1           ; [a+1]

    mla     r7, r11, r7, r9         ; c*192 + 128
    mla     r7, r4, r10, r7         ; a*64 + b*192 + 128

    subs    srcw, srcw, #3

    mov     r7, r7, asr #8          ; (a*64 + b*192 + 128) >> 8
    strb    r7, [dest], #1

    bpl     hl34_loop

    ldrb    r8, [src], #1           ; b = src[1]
    ldrb    r7, [src], #1           ; c = src[2]
    strb    r4, [dest], #1

    mla     r4, r10, r4, r9         ; a*64 + 128
    mla     r4, r11, r8, r4         ; a*64 + b*192 + 1
    mov     r4, r4, asr #8          ; (a*64 + b*192 + 1) >> 8
    strb    r4, [dest], #1

    add     r8, r8, #1              ; b + 1
    add     r8, r8, r7              ; b + c + 1
    mov     r8, r8, asr #1          ; (b + c + 1) >> 1
    strb    r8, [dest], #1
    strb    r7, [dest], #1

    ldmia   sp!, {r4 - r11, pc}
    ENDP    ;|vp8cx_horizontal_line_3_4_scale_c|


;/****************************************************************************
; *
; *  ROUTINE       : vertical_band_3_4_scale_armv4
; *
; *  INPUTS        : unsigned char *dest    : Pointer to destination data.
; *                  unsigned int dest_pitch : Stride of destination data.
; *                  unsigned int dest_width : Width of destination data.
; *
; *  OUTPUTS       : None.
; *
; *  RETURNS       : void
; *
; *  FUNCTION      : Scales vertical band of pixels by scale 3 to 4. The
; *                  height of the band scaled is 3-pixels.
; *
; *  SPECIAL NOTES : The routine uses the first line of the band below
; *                  the current band.
; *
; ****************************************************************************/
;void vertical_band_3_4_scale_armv4
;(
;   r0 = UINT8 *dest
;   r1 = UINT32 dest_pitch
;   r2 = UINT32 dest_width
;)
|vertical_band_3_4_scale_armv4| PROC
    stmdb   sp!, {r4 - r11, lr}

    ldr     r10, =64
    ldr     r11, =192
    mov     r9, #128

;   ldr     r1,[r1]
vl34_loop
    mov     r3, src
    ldrb    r4, [r3], r1            ; a = des [0]
    ldrb    r5, [r3], r1            ; b = des [dest_pitch]
    ldrb    r7, [r3], r1            ; c = des [dest_pitch*2]
    add     lr, src, r1

    mla     r4, r10, r4, r9         ; a*64 + 128
    mla     r4, r11, r5, r4         ; a*64 + b*192 + 1

    add     r5, r5, #1              ; b + 1
    add     r5, r5, r7              ; b + c + 1
    mov     r5, r5, asr #1          ; (b + c + 1) >> 1

    mov     r4, r4, asr #8          ; (a*64 + b*192 + 1) >> 8
    strb    r4, [lr], r1

    ldrb    r4, [r3, r1]            ; a = des [dest_pitch*4]

    strb    r5, [lr], r1

    mla     r7, r11, r7, r9         ; c*192 + 128
    mla     r7, r4, r10, r7         ; a*64 + b*192 + 128
    mov     r7, r7, asr #8          ; (a*64 + b*192 + 128) >> 8

    add     src, src, #1
    subs    r2, r2, #1

    strb    r7, [lr]

    bne     vl34_loop

    ldmia   sp!, {r4 - r11, pc}
    ENDP    ;|vertical_band_3_4_scale_armv4|

;/****************************************************************************
; *
; *  ROUTINE       : vp8cx_horizontal_line_1_2_scale_c
; *
; *  INPUTS        : const unsigned char *source : Pointer to source data.
; *                  unsigned int source_width    : Stride of source.
; *                  unsigned char *dest         : Pointer to destination data.
; *                  unsigned int dest_width      : Stride of destination (NOT USED).
; *
; *  OUTPUTS       : None.
; *
; *  RETURNS       : void
; *
; *  FUNCTION      : Copies horizontal line of pixels from source to
; *                  destination scaling up by 1 to 2.
; *
; *  SPECIAL NOTES : None.
; *
; ****************************************************************************/
;void vp8cx_horizontal_line_1_2_scale_c
;(
;   const unsigned char *source,
;   unsigned int source_width,
;   unsigned char *dest,
;   unsigned int dest_width
;)
|horizontal_line_1_2_scale_armv4| PROC
    stmdb   sp!, {r4 - r5, lr}

    sub     srcw, srcw, #1

    ldrb    r3, [src], #1
    ldrb    r4, [src], #1
hl12_loop
    subs    srcw, srcw, #1

    add     r5, r3, r4
    add     r5, r5, #1
    mov     r5, r5, lsr #1

    orr     r5, r3, r5, lsl #8
    strh    r5, [dest], #2

    mov     r3, r4

    ldrneb  r4, [src], #1
    bne     hl12_loop

    orr     r5, r4, r4, lsl #8
    strh    r5, [dest]

    ldmia   sp!, {r4 - r5, pc}
    ENDP    ;|vertical_band_3_5_scale_armv4|

;/****************************************************************************
; *
; *  ROUTINE       : vp8cx_vertical_band_1_2_scale_c
; *
; *  INPUTS        : unsigned char *dest    : Pointer to destination data.
; *                  unsigned int dest_pitch : Stride of destination data.
; *                  unsigned int dest_width : Width of destination data.
; *
; *  OUTPUTS       : None.
; *
; *  RETURNS       : void
; *
; *  FUNCTION      : Scales vertical band of pixels by scale 1 to 2. The
; *                  height of the band scaled is 1-pixel.
; *
; *  SPECIAL NOTES : The routine uses the first line of the band below
; *                  the current band.
; *
; ****************************************************************************/
;void vp8cx_vertical_band_1_2_scale_c
;(
;   r0 = UINT8 *dest
;   r1 = UINT32 dest_pitch
;   r2 = UINT32 dest_width
;)
|vertical_band_1_2_scale_armv4| PROC
    stmdb   sp!, {r4 - r7, lr}

    ldr     mask, =0xff00ff             ; mask for selection
    ldr     lr, = 0x010001

vl12_loop
    mov     r3, src
    ldr     r4, [r3], r1
    ldr     r5, [r3, r1]

    add     src, src, #4
    subs    r2, r2, #4

    and     r6, r4, mask
    and     r7, r5, mask

    add     r6, r7, r6
    add     r6, r6, lr

    and     r4, mask, r4, lsr #8
    and     r5, mask, r5, lsr #8

    mov     r6, r6, lsr #1
    and     r6, r6, mask

    add     r4, r5, r4
    add     r4, r4, lr

    mov     r4, r4, lsr #1
    and     r4, r4, mask

    orr     r5, r6, r4, lsl #8

    str     r5, [r3]

    bpl     vl12_loop

    ldmia   sp!, {r4 - r7, pc}
    ENDP    ;|vertical_band_3_5_scale_armv4|

    END
