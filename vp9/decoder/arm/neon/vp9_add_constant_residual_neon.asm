;
;   Copyright (c) 2013 The WebM project authors. All Rights Reserved.
;
;   Use of this source code is governed by a BSD-style license
;   that can be found in the LICENSE file in the root of the source
;   tree. An additional intellectual property rights grant can be found
;   in the file PATENTS.  All contributing project authors may
;   be found in the AUTHORS file in the root of the source tree.
;

    EXPORT |vp9_add_constant_residual_8x8_neon|
    EXPORT |vp9_add_constant_residual_16x16_neon|
    EXPORT |vp9_add_constant_residual_32x32_neon|
    ARM

    AREA ||.text||, CODE, READONLY, ALIGN=2

    MACRO
    LD_16x8 $src, $stride
    vld1.8              {q8},       [$src],     $stride
    vld1.8              {q9},       [$src],     $stride
    vld1.8              {q10},      [$src],     $stride
    vld1.8              {q11},      [$src],     $stride
    vld1.8              {q12},      [$src],     $stride
    vld1.8              {q13},      [$src],     $stride
    vld1.8              {q14},      [$src],     $stride
    vld1.8              {q15},      [$src],     $stride
    MEND

    MACRO
    ADD_DIFF_16x8 $diff
    vqadd.u8            q8,         q8,         $diff
    vqadd.u8            q9,         q9,         $diff
    vqadd.u8            q10,        q10,        $diff
    vqadd.u8            q11,        q11,        $diff
    vqadd.u8            q12,        q12,        $diff
    vqadd.u8            q13,        q13,        $diff
    vqadd.u8            q14,        q14,        $diff
    vqadd.u8            q15,        q15,        $diff
    MEND

    MACRO
    SUB_DIFF_16x8 $diff
    vqsub.u8            q8,         q8,         $diff
    vqsub.u8            q9,         q9,         $diff
    vqsub.u8            q10,        q10,        $diff
    vqsub.u8            q11,        q11,        $diff
    vqsub.u8            q12,        q12,        $diff
    vqsub.u8            q13,        q13,        $diff
    vqsub.u8            q14,        q14,        $diff
    vqsub.u8            q15,        q15,        $diff
    MEND

    MACRO
    ST_16x8 $dst, $stride
    vst1.8              {q8},       [$dst],     $stride
    vst1.8              {q9},       [$dst],     $stride
    vst1.8              {q10},      [$dst],     $stride
    vst1.8              {q11},      [$dst],     $stride
    vst1.8              {q12},      [$dst],     $stride
    vst1.8              {q13},      [$dst],     $stride
    vst1.8              {q14},      [$dst],     $stride
    vst1.8              {q15},      [$dst],     $stride
    MEND

; void add_constant_residual(const int16_t diff, uint8_t *dest, int stride,
;                             int width, int height) {
;  int r, c;
;
;  for (r = 0; r < height; r++) {
;    for (c = 0; c < width; c++)
;      dest[c] = clip_pixel(diff + dest[c]);
;
;    dest += stride;
;  }
;}
;void vp9_add_constant_residual_8x8_c(const int16_t diff, uint8_t *dest,
;                                     int stride) {
;  add_constant_residual(diff, dest, stride, 8, 8);
;}
;       r0      : const int16_t diff
;       r1      : const uint8_t *dest
;       r2      : int stride
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
|vp9_add_constant_residual_8x8_neon| PROC
    mov                 r3,         r1                      ; r3: save dest to r3
    vld1.8              {d0},       [r1],       r2
    vld1.8              {d1},       [r1],       r2
    vld1.8              {d2},       [r1],       r2
    vld1.8              {d3},       [r1],       r2
    vld1.8              {d4},       [r1],       r2
    vld1.8              {d5},       [r1],       r2
    vld1.8              {d6},       [r1],       r2
    vld1.8              {d7},       [r1],       r2
    cmp                 r0,         #0
    bge                 DIFF_POSITIVE_8x8

DIFF_NEGATIVE_8x8                                           ; diff < 0
    neg                 r0,         r0
    usat                r0,         #8,         r0
    vdup.u8             q8,         r0

    vqsub.u8            q0,         q0,         q8
    vqsub.u8            q1,         q1,         q8
    vqsub.u8            q2,         q2,         q8
    vqsub.u8            q3,         q3,         q8
    b                   DIFF_SAVE_8x8

DIFF_POSITIVE_8x8                                           ; diff >= 0
    usat                r0,         #8,         r0
    vdup.u8             q8,         r0

    vqadd.u8            q0,         q0,         q8
    vqadd.u8            q1,         q1,         q8
    vqadd.u8            q2,         q2,         q8
    vqadd.u8            q3,         q3,         q8

DIFF_SAVE_8x8
    vst1.8              {d0},       [r3],       r2
    vst1.8              {d1},       [r3],       r2
    vst1.8              {d2},       [r3],       r2
    vst1.8              {d3},       [r3],       r2
    vst1.8              {d4},       [r3],       r2
    vst1.8              {d5},       [r3],       r2
    vst1.8              {d6},       [r3],       r2
    vst1.8              {d7},       [r3],       r2

    bx                  lr
    ENDP

;void vp9_add_constant_residual_16x16_c(const int16_t diff, uint8_t *dest,
;                                       int stride) {
;  add_constant_residual(diff, dest, stride, 16, 16);
;}
;       r0      : const int16_t diff
;       r1      : const uint8_t *dest
;       r2      : int stride
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
|vp9_add_constant_residual_16x16_neon| PROC
    mov                 r3,         r1
    LD_16x8             r1,         r2
    cmp                 r0,         #0
    bge                 DIFF_POSITIVE_16x16

|DIFF_NEGATIVE_16x16|
    neg                 r0,         r0
    usat                r0,         #8,         r0
    vdup.u8             q0,         r0

    SUB_DIFF_16x8       q0
    ST_16x8             r3,         r2
    LD_16x8             r1,         r2
    SUB_DIFF_16x8       q0
    b                   DIFF_SAVE_16x16

|DIFF_POSITIVE_16x16|
    usat                r0,         #8,         r0
    vdup.u8             q0,         r0

    ADD_DIFF_16x8       q0
    ST_16x8             r3,         r2
    LD_16x8             r1,         r2
    ADD_DIFF_16x8       q0

|DIFF_SAVE_16x16|
    ST_16x8             r3,         r2
    bx                  lr
    ENDP

;void vp9_add_constant_residual_32x32_c(const int16_t diff, uint8_t *dest,
;                                       int stride) {
;  add_constant_residual(diff, dest, stride, 32, 32);
;}
;       r0      : const int16_t diff
;       r1      : const uint8_t *dest
;       r2      : int stride
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
|vp9_add_constant_residual_32x32_neon| PROC
    push                {r4,lr}
    pld                 [r1]
    mov                 r3,         r1
    add                 r4,         r1,         #16         ; r4 dest + 16 for second loop
    cmp                 r0,         #0
    bge                 DIFF_POSITIVE_32x32

|DIFF_NEGATIVE_32x32|
    neg                 r0,         r0
    usat                r0,         #8,         r0
    vdup.u8             q0,         r0
    mov                 r0,         #4

|DIFF_NEGATIVE_32x32_LOOP|
    sub                 r0,         #1
    LD_16x8             r1,         r2
    SUB_DIFF_16x8       q0
    ST_16x8             r3,         r2

    LD_16x8             r1,         r2
    SUB_DIFF_16x8       q0
    ST_16x8             r3,         r2
    cmp                 r0,         #2
    moveq               r1,         r4
    moveq               r3,         r4
    cmp                 r0,         #0
    bne                 DIFF_NEGATIVE_32x32_LOOP
    pop                 {r4,pc}

|DIFF_POSITIVE_32x32|
    usat                r0,         #8,         r0
    vdup.u8             q0,         r0
    mov                 r0,         #4

|DIFF_POSITIVE_32x32_LOOP|
    sub                 r0,         #1
    LD_16x8             r1,         r2
    ADD_DIFF_16x8       q0
    ST_16x8             r3,         r2

    LD_16x8             r1,         r2
    ADD_DIFF_16x8       q0
    ST_16x8             r3,         r2
    cmp                 r0,         #2
    moveq               r1,         r4
    moveq               r3,         r4
    cmp                 r0,         #0
    bne                 DIFF_POSITIVE_32x32_LOOP
    pop                 {r4,pc}
    ENDP

    END
