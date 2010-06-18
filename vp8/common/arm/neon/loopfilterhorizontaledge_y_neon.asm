;
;  Copyright (c) 2010 The VP8 project authors. All Rights Reserved.
;
;  Use of this source code is governed by a BSD-style license
;  that can be found in the LICENSE file in the root of the source
;  tree. An additional intellectual property rights grant can be found
;  in the file PATENTS.  All contributing project authors may
;  be found in the AUTHORS file in the root of the source tree.
;


    EXPORT  |vp8_loop_filter_horizontal_edge_y_neon|
    ARM
    REQUIRE8
    PRESERVE8

    AREA ||.text||, CODE, READONLY, ALIGN=2
;Note: flimit, limit, and thresh shpuld be positive numbers. All 16 elements in flimit
;are equal. So, in the code, only one load is needed
;for flimit. Same way applies to limit and thresh.
; r0    unsigned char *s,
; r1    int p, //pitch
; r2    const signed char *flimit,
; r3    const signed char *limit,
; stack(r4) const signed char *thresh,
; //stack(r5)   int count --unused

|vp8_loop_filter_horizontal_edge_y_neon| PROC
    sub         r0, r0, r1, lsl #2          ; move src pointer down by 4 lines
    ldr         r12, [sp, #0]               ; load thresh pointer

    vld1.u8     {q3}, [r0], r1              ; p3
    vld1.s8     {d0[], d1[]}, [r2]          ; flimit
    vld1.u8     {q4}, [r0], r1              ; p2
    vld1.s8     {d2[], d3[]}, [r3]          ; limit
    vld1.u8     {q5}, [r0], r1              ; p1
    vld1.s8     {d4[], d5[]}, [r12]         ; thresh
    vld1.u8     {q6}, [r0], r1              ; p0
    ldr         r12, _lfhy_coeff_
    vld1.u8     {q7}, [r0], r1              ; q0

    ;vp8_filter_mask() function
    ;vp8_hevmask() function
    vabd.u8     q11, q3, q4                 ; abs(p3 - p2)
    vld1.u8     {q8}, [r0], r1              ; q1
    vabd.u8     q12, q4, q5                 ; abs(p2 - p1)
    vld1.u8     {q9}, [r0], r1              ; q2
    vabd.u8     q13, q5, q6                 ; abs(p1 - p0)
    vld1.u8     {q10}, [r0], r1             ; q3
    vabd.u8     q14, q8, q7                 ; abs(q1 - q0)
    vabd.u8     q3, q9, q8                  ; abs(q2 - q1)
    vabd.u8     q4, q10, q9                 ; abs(q3 - q2)
    vabd.u8     q9, q6, q7                  ; abs(p0 - q0)

    vcge.u8     q15, q1, q11                ; (abs(p3 - p2) > limit)*-1
    vcge.u8     q12, q1, q12                ; (abs(p2 - p1) > limit)*-1
    vcge.u8     q10, q1, q13                ; (abs(p1 - p0) > limit)*-1
    vcge.u8     q11, q1, q14                ; (abs(q1 - q0) > limit)*-1

    vcgt.u8     q13, q13, q2                ; (abs(p1 - p0) > thresh)*-1
    vcgt.u8     q14, q14, q2                ; (abs(q1 - q0) > thresh)*-1

    vcge.u8     q3, q1, q3                  ; (abs(q2 - q1) > limit)*-1
    vcge.u8     q4, q1, q4                  ; (abs(q3 - q2) > limit)*-1
    vadd.u8     q0, q0, q0                  ; flimit * 2
    vadd.u8     q0, q0, q1                  ; flimit * 2 + limit

    vand        q15, q15, q12
    vand        q10, q10, q11
    vand        q3, q3, q4

    vabd.u8     q2, q5, q8                  ; abs(p1 - q1)
    vqadd.u8    q9, q9, q9                  ; abs(p0 - q0) * 2
    vshr.u8     q2, q2, #1                  ; abs(p1 - q1) / 2
    vqadd.u8    q9, q9, q2                  ; abs(p0 - q0) * 2 + abs(p1 - q1) / 2
    vcge.u8     q9, q0, q9                  ; (abs(p0 - q0)*2 + abs(p1-q1)/2 > flimit*2 + limit)*-1

    vld1.u8     {q0}, [r12]!

    vand        q15, q15, q10

    ;vp8_filter() function
    veor        q7, q7, q0                  ; qs0: q0 offset to convert to a signed value
    veor        q6, q6, q0                  ; ps0: p0 offset to convert to a signed value
    veor        q5, q5, q0                  ; ps1: p1 offset to convert to a signed value
    veor        q8, q8, q0                  ; qs1: q1 offset to convert to a signed value
;;;;;;;;;;;;;;
    vld1.u8     {q10}, [r12]!

    ;vqsub.s8   q2, q7, q6                  ; ( qs0 - ps0)
    vsubl.s8    q2, d14, d12                ; ( qs0 - ps0)
    vsubl.s8    q11, d15, d13

    vand        q3, q3, q9
    vmovl.u8    q4, d20

    vqsub.s8    q1, q5, q8                  ; vp8_filter = vp8_signed_char_clamp(ps1-qs1)
    vorr        q14, q13, q14               ; q14: vp8_hevmask

    ;vmul.i8    q2, q2, q10                 ; 3 * ( qs0 - ps0)
    vmul.i16    q2, q2, q4                  ; 3 * ( qs0 - ps0)
    vmul.i16    q11, q11, q4

    vand        q1, q1, q14                 ; vp8_filter &= hev
    vand        q15, q15, q3                ; q15: vp8_filter_mask
    ;;
    ;vld1.u8        {q4}, [r12]!            ;no need 7 any more

    ;vqadd.s8   q1, q1, q2
    vaddw.s8    q2, q2, d2
    vaddw.s8    q11, q11, d3

    vld1.u8     {q9}, [r12]!
    ;
    vqmovn.s16  d2, q2                      ; vp8_filter = vp8_signed_char_clamp(vp8_filter + 3 * ( qs0 - ps0))
    vqmovn.s16  d3, q11
    ;;

    vand        q1, q1, q15                 ; vp8_filter &= mask
    ;;
;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;Change for VP8 from VP7
;   vand        q2, q1, q4                  ; s = vp8_filter & 7
;   vqadd.s8    q1, q1, q9                  ; vp8_filter = vp8_signed_char_clamp(vp8_filter+4)
    ;;;;
;   vshr.s8     q1, q1, #3                  ; vp8_filter >>= 3
;   vceq.i8     q2, q2, q9                  ; s = (s==4)*-1
    ;;
;   ;calculate output
;   vqsub.s8    q10, q7, q1                 ; u = vp8_signed_char_clamp(qs0 - vp8_filter)
;   vqadd.s8    q11, q2, q1                 ; u = vp8_signed_char_clamp(s + vp8_filter)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; q10=3
    vqadd.s8    q2, q1, q10                 ; Filter2 = vp8_signed_char_clamp(vp8_filter+3)
    vqadd.s8    q1, q1, q9                  ; Filter1 = vp8_signed_char_clamp(vp8_filter+4)
    vshr.s8     q2, q2, #3                  ; Filter2 >>= 3
    vshr.s8     q1, q1, #3                  ; Filter1 >>= 3

    ;calculate output
    vqadd.s8    q11, q6, q2                 ; u = vp8_signed_char_clamp(ps0 + Filter2)
    vqsub.s8    q10, q7, q1                 ; u = vp8_signed_char_clamp(qs0 - Filter1)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

    vrshr.s8    q1, q1, #1                  ;round/shift:  vp8_filter += 1; vp8_filter >>= 1

    sub         r0, r0, r1, lsl #2
    sub         r0, r0, r1, lsl #1
    ;

    vbic        q1, q1, q14                 ; vp8_filter &= ~hev
    ;
    add         r2, r1, r0

    vqadd.s8    q13, q5, q1                 ; u = vp8_signed_char_clamp(ps1 + vp8_filter)
    ;vqadd.s8   q11, q6, q11                ; u = vp8_signed_char_clamp(ps0 + u)
    vqsub.s8    q12, q8, q1                 ; u = vp8_signed_char_clamp(qs1 - vp8_filter)

    add         r3, r2, r1

    veor        q5, q13, q0                 ; *op1 = u^0x80
    veor        q6, q11, q0                 ; *op0 = u^0x80
    veor        q7, q10, q0                 ; *oq0 = u^0x80
    veor        q8, q12, q0                 ; *oq1 = u^0x80

    add         r12, r3, r1

    vst1.u8     {q5}, [r0]                  ; store op1
    vst1.u8     {q6}, [r2]                  ; store op0
    vst1.u8     {q7}, [r3]                  ; store oq0
    vst1.u8     {q8}, [r12]                 ; store oq1

    bx          lr
    ENDP        ; |vp8_loop_filter_horizontal_edge_y_neon|

;-----------------
    AREA    hloopfiltery_dat, DATA, READWRITE           ;read/write by default
;Data section with name data_area is specified. DCD reserves space in memory for 16 data.
;One word each is reserved. Label filter_coeff can be used to access the data.
;Data address: filter_coeff, filter_coeff+4, filter_coeff+8 ...
_lfhy_coeff_
    DCD     lfhy_coeff
lfhy_coeff
    DCD     0x80808080, 0x80808080, 0x80808080, 0x80808080
    DCD     0x03030303, 0x03030303, 0x03030303, 0x03030303
    DCD     0x04040404, 0x04040404, 0x04040404, 0x04040404
    DCD     0x01010101, 0x01010101, 0x01010101, 0x01010101

    END
