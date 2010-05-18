;
;  Copyright (c) 2010 The VP8 project authors. All Rights Reserved.
;
;  Use of this source code is governed by a BSD-style license and patent
;  grant that can be found in the LICENSE file in the root of the source
;  tree. All contributing project authors may be found in the AUTHORS
;  file in the root of the source tree.
;


    EXPORT  |vp8_loop_filter_vertical_edge_uv_neon|
    ARM
    REQUIRE8
    PRESERVE8

    AREA ||.text||, CODE, READONLY, ALIGN=2
;Note: flimit, limit, and thresh shpuld be positive numbers. All 16 elements in flimit
;are equal. So, in the code, only one load is needed
;for flimit. Same way applies to limit and thresh.
; r0    unsigned char *u,
; r1    int p, //pitch
; r2    const signed char *flimit,
; r3    const signed char *limit,
; stack(r4) const signed char *thresh,
; stack(r5) unsigned char *v

|vp8_loop_filter_vertical_edge_uv_neon| PROC
    sub         r0, r0, #4          ; move u pointer down by 4 columns
    vld1.s8     {d0[], d1[]}, [r2]          ; flimit

    ldr         r2, [sp, #4]                ; load v ptr
    ldr         r12, [sp, #0]               ; load thresh pointer

    sub         r2, r2, #4          ; move v pointer down by 4 columns

    vld1.u8     {d6}, [r0], r1              ;load u data
    vld1.u8     {d7}, [r2], r1              ;load v data
    vld1.u8     {d8}, [r0], r1
    vld1.u8     {d9}, [r2], r1
    vld1.u8     {d10}, [r0], r1
    vld1.u8     {d11}, [r2], r1
    vld1.u8     {d12}, [r0], r1
    vld1.u8     {d13}, [r2], r1
    vld1.u8     {d14}, [r0], r1
    vld1.u8     {d15}, [r2], r1
    vld1.u8     {d16}, [r0], r1
    vld1.u8     {d17}, [r2], r1
    vld1.u8     {d18}, [r0], r1
    vld1.u8     {d19}, [r2], r1
    vld1.u8     {d20}, [r0], r1
    vld1.u8     {d21}, [r2], r1

    ;transpose to 8x16 matrix
    vtrn.32     q3, q7
    vtrn.32     q4, q8
    vtrn.32     q5, q9
    vtrn.32     q6, q10

    vtrn.16     q3, q5
    vtrn.16     q4, q6
    vtrn.16     q7, q9
    vtrn.16     q8, q10

    vtrn.8      q3, q4
    vtrn.8      q5, q6
    vtrn.8      q7, q8
    vtrn.8      q9, q10

    vld1.s8     {d2[], d3[]}, [r3]          ; limit
    vld1.s8     {d4[], d5[]}, [r12]         ; thresh

    ldr         r12, _vlfuv_coeff_
    ;vp8_filter_mask() function
    ;vp8_hevmask() function
    vabd.u8     q11, q3, q4                 ; abs(p3 - p2)
    vabd.u8     q12, q4, q5                 ; abs(p2 - p1)
    vabd.u8     q13, q5, q6                 ; abs(p1 - p0)
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
    vqadd.s8    q11, q6, q2             ; u = vp8_signed_char_clamp(ps0 + Filter2)
    vqsub.s8    q10, q7, q1                 ; u = vp8_signed_char_clamp(qs0 - Filter1)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

    vrshr.s8    q1, q1, #1                  ;round/shift:  vp8_filter += 1; vp8_filter >>= 1

    sub         r0, r0, r1, lsl #3
    add         r0, r0, #2

    vbic        q1, q1, q14                 ; vp8_filter &= ~hev

    sub         r2, r2, r1, lsl #3
    add         r2, r2, #2

    vqadd.s8    q13, q5, q1                 ; u = vp8_signed_char_clamp(ps1 + vp8_filter)
    ;vqadd.s8   q11, q6, q11                ; u = vp8_signed_char_clamp(ps0 + u)
    vqsub.s8    q12, q8, q1                 ; u = vp8_signed_char_clamp(qs1 - vp8_filter)

    veor        q7, q10, q0                 ; *oq0 = u^0x80
    veor        q5, q13, q0                 ; *op1 = u^0x80
    veor        q6, q11, q0                 ; *op0 = u^0x80
    veor        q8, q12, q0                 ; *oq1 = u^0x80

    vswp        d12, d11
    vswp        d16, d13
    vswp        d14, d12
    vswp        d16, d15

    ;store op1, op0, oq0, oq1
    vst4.8      {d10[0], d11[0], d12[0], d13[0]}, [r0], r1
    vst4.8      {d14[0], d15[0], d16[0], d17[0]}, [r2], r1
    vst4.8      {d10[1], d11[1], d12[1], d13[1]}, [r0], r1
    vst4.8      {d14[1], d15[1], d16[1], d17[1]}, [r2], r1
    vst4.8      {d10[2], d11[2], d12[2], d13[2]}, [r0], r1
    vst4.8      {d14[2], d15[2], d16[2], d17[2]}, [r2], r1
    vst4.8      {d10[3], d11[3], d12[3], d13[3]}, [r0], r1
    vst4.8      {d14[3], d15[3], d16[3], d17[3]}, [r2], r1
    vst4.8      {d10[4], d11[4], d12[4], d13[4]}, [r0], r1
    vst4.8      {d14[4], d15[4], d16[4], d17[4]}, [r2], r1
    vst4.8      {d10[5], d11[5], d12[5], d13[5]}, [r0], r1
    vst4.8      {d14[5], d15[5], d16[5], d17[5]}, [r2], r1
    vst4.8      {d10[6], d11[6], d12[6], d13[6]}, [r0], r1
    vst4.8      {d14[6], d15[6], d16[6], d17[6]}, [r2], r1
    vst4.8      {d10[7], d11[7], d12[7], d13[7]}, [r0], r1
    vst4.8      {d14[7], d15[7], d16[7], d17[7]}, [r2], r1

    bx          lr
    ENDP        ; |vp8_loop_filter_vertical_edge_uv_neon|

;-----------------
    AREA    vloopfilteruv_dat, DATA, READWRITE          ;read/write by default
;Data section with name data_area is specified. DCD reserves space in memory for 16 data.
;One word each is reserved. Label filter_coeff can be used to access the data.
;Data address: filter_coeff, filter_coeff+4, filter_coeff+8 ...
_vlfuv_coeff_
    DCD     vlfuv_coeff
vlfuv_coeff
    DCD     0x80808080, 0x80808080, 0x80808080, 0x80808080
    DCD     0x03030303, 0x03030303, 0x03030303, 0x03030303
    DCD     0x04040404, 0x04040404, 0x04040404, 0x04040404
    DCD     0x01010101, 0x01010101, 0x01010101, 0x01010101

    END
