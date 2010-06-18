;
;  Copyright (c) 2010 The VP8 project authors. All Rights Reserved.
;
;  Use of this source code is governed by a BSD-style license
;  that can be found in the LICENSE file in the root of the source
;  tree. An additional intellectual property rights grant can be found
;  in the file PATENTS.  All contributing project authors may
;  be found in the AUTHORS file in the root of the source tree.
;


    EXPORT  |vp8_mbloop_filter_vertical_edge_uv_neon|
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
|vp8_mbloop_filter_vertical_edge_uv_neon| PROC
    sub         r0, r0, #4                  ; move src pointer down by 4 columns
    vld1.s8     {d2[], d3[]}, [r3]          ; limit
    ldr         r3, [sp, #4]                ; load v ptr
    ldr         r12, [sp, #0]               ; load thresh pointer

    sub         r3, r3, #4                  ; move v pointer down by 4 columns

    vld1.u8     {d6}, [r0], r1              ;load u data
    vld1.u8     {d7}, [r3], r1              ;load v data
    vld1.u8     {d8}, [r0], r1
    vld1.u8     {d9}, [r3], r1
    vld1.u8     {d10}, [r0], r1
    vld1.u8     {d11}, [r3], r1
    vld1.u8     {d12}, [r0], r1
    vld1.u8     {d13}, [r3], r1
    vld1.u8     {d14}, [r0], r1
    vld1.u8     {d15}, [r3], r1
    vld1.u8     {d16}, [r0], r1
    vld1.u8     {d17}, [r3], r1
    vld1.u8     {d18}, [r0], r1
    vld1.u8     {d19}, [r3], r1
    vld1.u8     {d20}, [r0], r1
    vld1.u8     {d21}, [r3], r1

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

    sub         sp, sp, #32
    vld1.s8     {d4[], d5[]}, [r12]         ; thresh
    vst1.u8     {q3}, [sp]!
    ldr         r12, _mbvlfuv_coeff_
    vst1.u8     {q10}, [sp]!

    ;vp8_filter_mask() function
    ;vp8_hevmask() function
    vabd.u8     q11, q3, q4                 ; abs(p3 - p2)
    vabd.u8     q12, q4, q5                 ; abs(p2 - p1)
    vabd.u8     q13, q5, q6                 ; abs(p1 - p0)
    vabd.u8     q14, q8, q7                 ; abs(q1 - q0)
    vabd.u8     q3, q9, q8                  ; abs(q2 - q1)
    vabd.u8     q0, q10, q9                 ; abs(q3 - q2)

    vcge.u8     q15, q1, q11                ; (abs(p3 - p2) > limit)*-1
    vcge.u8     q12, q1, q12                ; (abs(p2 - p1) > limit)*-1
    vcge.u8     q10, q1, q13                ; (abs(p1 - p0) > limit)*-1
    vcge.u8     q11, q1, q14                ; (abs(q1 - q0) > limit)*-1
    vcge.u8     q3, q1, q3                  ; (abs(q2 - q1) > limit)*-1
    vcge.u8     q0, q1, q0                  ; (abs(q3 - q2) > limit)*-1

    vand        q15, q15, q12

    vabd.u8     q12, q6, q7                 ; abs(p0 - q0)

    vcgt.u8     q13, q13, q2                ; (abs(p1 - p0) > thresh)*-1
    vcgt.u8     q14, q14, q2                ; (abs(q1 - q0) > thresh)*-1

    vld1.s8     {d4[], d5[]}, [r2]          ; flimit

    vand        q10, q10, q11
    vand        q3, q3, q0

    vld1.u8     {q0}, [r12]!

    vadd.u8     q2, q2, q2                  ; flimit * 2
    vadd.u8     q2, q2, q1                  ; flimit * 2 + limit

    vabd.u8     q1, q5, q8                  ; abs(p1 - q1)
    vqadd.u8    q12, q12, q12               ; abs(p0 - q0) * 2
    vshr.u8     q1, q1, #1                  ; abs(p1 - q1) / 2
    vqadd.u8    q12, q12, q1                ; abs(p0 - q0) * 2 + abs(p1 - q1) / 2
    vcge.u8     q12, q2, q12                ; (abs(p0 - q0)*2 + abs(p1 - q1)/2 > flimit*2 + limit)*-1

    vand        q15, q15, q10

    ;vp8_filter() function
    veor        q7, q7, q0                  ; qs0: q0 offset to convert to a signed value
    veor        q6, q6, q0                  ; ps0: p0 offset to convert to a signed value
    veor        q5, q5, q0                  ; ps1: p1 offset to convert to a signed value
    veor        q8, q8, q0                  ; qs1: q1 offset to convert to a signed value
    veor        q4, q4, q0                  ; ps2: p2 offset to convert to a signed value
    veor        q9, q9, q0                  ; qs2: q2 offset to convert to a signed value
;;;;;;;;;;;;;
    vorr        q14, q13, q14               ; q14: vp8_hevmask

    ;vqsub.s8   q2, q7, q6                  ; ( qs0 - ps0)
    vsubl.s8    q2, d14, d12                ; ( qs0 - ps0)
    vsubl.s8    q13, d15, d13

    vqsub.s8    q1, q5, q8                  ; vp8_filter = vp8_signed_char_clamp(ps1-qs1)

    ;vadd.s8    q10, q2, q2                 ; 3 * ( qs0 - ps0)
    vadd.s16    q10, q2, q2                 ; 3 * ( qs0 - ps0)
    vadd.s16    q11, q13, q13

    vand        q3, q3, q12

    ;vadd.s8    q2, q2, q10
    vadd.s16    q2, q2, q10
    vadd.s16    q13, q13, q11

    vld1.u8     {q12}, [r12]!               ;#3

    ;vqadd.s8   q1, q1, q2                  ; vp8_filter + 3 * ( qs0 - ps0)
    vaddw.s8    q2, q2, d2                  ; vp8_filter + 3 * ( qs0 - ps0)
    vaddw.s8    q13, q13, d3

    vand        q15, q15, q3                ; q15: vp8_filter_mask
    vld1.u8     {q11}, [r12]!               ;#4

    vqmovn.s16  d2, q2                      ; vp8_filter = vp8_signed_char_clamp(vp8_filter + 3 * ( qs0 - ps0))
    vqmovn.s16  d3, q13

;;;;;;;;;;;;;;
    vand        q1, q1, q15                 ; vp8_filter &= mask

    vld1.u8     {q15}, [r12]!               ;#63
    ;
    vand        q13, q1, q14                ; Filter2: q13; Filter2 &= hev

    vld1.u8     {d7}, [r12]!                ;#9
    ;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;Change for VP8 from VP7
;   vand        q2, q13, q12                ; s = Filter2 & 7

;   vqadd.s8    q13, q13, q11               ; Filter2 = vp8_signed_char_clamp(Filter2+4)
;   vld1.u8     {d6}, [r12]!                ;#18

;   sub         r0, r0, r1, lsl #3
;   sub         r3, r3, r1, lsl #3
;   sub         sp, sp, #32

;   vshr.s8     q13, q13, #3                ; Filter2 >>= 3
;   vceq.i8     q2, q2, q11                 ; s = (s==4)*-1

;   vqsub.s8    q7, q7, q13                 ; qs0 = vp8_signed_char_clamp(qs0 - Filter2)
;   vqadd.s8    q11, q2, q13                ; u = vp8_signed_char_clamp(s + Filter2)

;   vld1.u8     {d5}, [r12]!                ;#27
;   vmov        q10, q15
;   vmov        q12, q15

;   vqadd.s8    q6, q6, q11                 ; ps0 = vp8_signed_char_clamp(ps0 + u)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    vqadd.s8    q2, q13, q11                ; Filter1 = vp8_signed_char_clamp(Filter2+4)
    vqadd.s8    q13, q13, q12               ; Filter2 = vp8_signed_char_clamp(Filter2+3)

    vld1.u8     {d6}, [r12]!                ;#18

    sub         r0, r0, r1, lsl #3
    sub         r3, r3, r1, lsl #3

    vshr.s8     q2, q2, #3                  ; Filter1 >>= 3
    vshr.s8     q13, q13, #3                ; Filter2 >>= 3

    vmov        q10, q15
    vmov        q12, q15

    vqsub.s8    q7, q7, q2                  ; qs0 = vp8_signed_char_clamp(qs0 - Filter1)

    vld1.u8     {d5}, [r12]!                ;#27

    sub         sp, sp, #32

    vqadd.s8    q6, q6, q13                 ; ps0 = vp8_signed_char_clamp(ps0 + Filter2)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

    vbic        q1, q1, q14                 ; Filter2: q1; vp8_filter &= ~hev; Filter2 = vp8_filter

    ; roughly 1/7th difference across boundary
    ; roughly 2/7th difference across boundary
    ; roughly 3/7th difference across boundary
    vmov        q11, q15
    vmov        q13, q15
    vmov        q14, q15

    vmlal.s8    q10, d2, d7                 ; Filter2 * 9
    vmlal.s8    q11, d3, d7
    vmlal.s8    q12, d2, d6                 ; Filter2 * 18
    vmlal.s8    q13, d3, d6
    vmlal.s8    q14, d2, d5                 ; Filter2 * 27
    vmlal.s8    q15, d3, d5
    vqshrn.s16  d20, q10, #7                ; u = vp8_signed_char_clamp((63 + Filter2 * 9)>>7)
    vqshrn.s16  d21, q11, #7
    vqshrn.s16  d24, q12, #7                ; u = vp8_signed_char_clamp((63 + Filter2 * 18)>>7)
    vqshrn.s16  d25, q13, #7
    vqshrn.s16  d28, q14, #7                ; u = vp8_signed_char_clamp((63 + Filter2 * 27)>>7)
    vqshrn.s16  d29, q15, #7

    vqsub.s8    q11, q9, q10                ; s = vp8_signed_char_clamp(qs2 - u)
    vqadd.s8    q10, q4, q10                ; s = vp8_signed_char_clamp(ps2 + u)
    vqsub.s8    q13, q8, q12                ; s = vp8_signed_char_clamp(qs1 - u)
    vqadd.s8    q12, q5, q12                ; s = vp8_signed_char_clamp(ps1 + u)
    vqsub.s8    q15, q7, q14                ; s = vp8_signed_char_clamp(qs0 - u)
    vqadd.s8    q14, q6, q14                ; s = vp8_signed_char_clamp(ps0 + u)
    veor        q9, q11, q0                 ; *oq2 = s^0x80
    veor        q4, q10, q0                 ; *op2 = s^0x80
    veor        q8, q13, q0                 ; *oq1 = s^0x80
    veor        q5, q12, q0                 ; *op2 = s^0x80
    veor        q7, q15, q0                 ; *oq0 = s^0x80
    vld1.u8     {q3}, [sp]!
    veor        q6, q14, q0                 ; *op0 = s^0x80
    vld1.u8     {q10}, [sp]!

    ;transpose to 16x8 matrix
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

    ;store op2, op1, op0, oq0, oq1, oq2
    vst1.8      {d6}, [r0], r1
    vst1.8      {d7}, [r3], r1
    vst1.8      {d8}, [r0], r1
    vst1.8      {d9}, [r3], r1
    vst1.8      {d10}, [r0], r1
    vst1.8      {d11}, [r3], r1
    vst1.8      {d12}, [r0], r1
    vst1.8      {d13}, [r3], r1
    vst1.8      {d14}, [r0], r1
    vst1.8      {d15}, [r3], r1
    vst1.8      {d16}, [r0], r1
    vst1.8      {d17}, [r3], r1
    vst1.8      {d18}, [r0], r1
    vst1.8      {d19}, [r3], r1
    vst1.8      {d20}, [r0], r1
    vst1.8      {d21}, [r3], r1

    bx          lr
    ENDP        ; |vp8_mbloop_filter_vertical_edge_uv_neon|

;-----------------
    AREA    mbvloopfilteruv_dat, DATA, READWRITE            ;read/write by default
;Data section with name data_area is specified. DCD reserves space in memory for 16 data.
;One word each is reserved. Label filter_coeff can be used to access the data.
;Data address: filter_coeff, filter_coeff+4, filter_coeff+8 ...
_mbvlfuv_coeff_
    DCD     mbvlfuv_coeff
mbvlfuv_coeff
    DCD     0x80808080, 0x80808080, 0x80808080, 0x80808080
    DCD     0x03030303, 0x03030303, 0x03030303, 0x03030303
    DCD     0x04040404, 0x04040404, 0x04040404, 0x04040404
    DCD     0x003f003f, 0x003f003f, 0x003f003f, 0x003f003f
    DCD     0x09090909, 0x09090909, 0x12121212, 0x12121212
    DCD     0x1b1b1b1b, 0x1b1b1b1b

    END
