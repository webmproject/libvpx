;
;  Copyright (c) 2014 The WebM project authors. All Rights Reserved.
;
;  Use of this source code is governed by a BSD-style license
;  that can be found in the LICENSE file in the root of the source
;  tree. An additional intellectual property rights grant can be found
;  in the file PATENTS.  All contributing project authors may
;  be found in the AUTHORS file in the root of the source tree.
;

    EXPORT  |vp9_v_predictor_4x4_neon|
    EXPORT  |vp9_v_predictor_8x8_neon|
    EXPORT  |vp9_v_predictor_16x16_neon|
    EXPORT  |vp9_v_predictor_32x32_neon|
    EXPORT  |vp9_h_predictor_4x4_neon|
    EXPORT  |vp9_h_predictor_8x8_neon|
    EXPORT  |vp9_h_predictor_16x16_neon|
    EXPORT  |vp9_h_predictor_32x32_neon|
    EXPORT  |vp9_tm_predictor_4x4_neon|
    ARM
    REQUIRE8
    PRESERVE8

    AREA ||.text||, CODE, READONLY, ALIGN=2

;void vp9_v_predictor_4x4_neon(uint8_t *dst, ptrdiff_t y_stride,
;                              const uint8_t *above,
;                              const uint8_t *left)
; r0  uint8_t *dst
; r1  ptrdiff_t y_stride
; r2  const uint8_t *above
; r3  const uint8_t *left

|vp9_v_predictor_4x4_neon| PROC
    vld1.32             {d0[0]}, [r2]
    vst1.32             {d0[0]}, [r0], r1
    vst1.32             {d0[0]}, [r0], r1
    vst1.32             {d0[0]}, [r0], r1
    vst1.32             {d0[0]}, [r0], r1
    bx                  lr
    ENDP                ; |vp9_v_predictor_4x4_neon|

;void vp9_v_predictor_8x8_neon(uint8_t *dst, ptrdiff_t y_stride,
;                              const uint8_t *above,
;                              const uint8_t *left)
; r0  uint8_t *dst
; r1  ptrdiff_t y_stride
; r2  const uint8_t *above
; r3  const uint8_t *left

|vp9_v_predictor_8x8_neon| PROC
    vld1.8              {d0}, [r2]
    vst1.8              {d0}, [r0], r1
    vst1.8              {d0}, [r0], r1
    vst1.8              {d0}, [r0], r1
    vst1.8              {d0}, [r0], r1
    vst1.8              {d0}, [r0], r1
    vst1.8              {d0}, [r0], r1
    vst1.8              {d0}, [r0], r1
    vst1.8              {d0}, [r0], r1
    bx                  lr
    ENDP                ; |vp9_v_predictor_8x8_neon|

;void vp9_v_predictor_16x16_neon(uint8_t *dst, ptrdiff_t y_stride,
;                                const uint8_t *above,
;                                const uint8_t *left)
; r0  uint8_t *dst
; r1  ptrdiff_t y_stride
; r2  const uint8_t *above
; r3  const uint8_t *left

|vp9_v_predictor_16x16_neon| PROC
    vld1.8              {q0}, [r2]
    vst1.8              {q0}, [r0], r1
    vst1.8              {q0}, [r0], r1
    vst1.8              {q0}, [r0], r1
    vst1.8              {q0}, [r0], r1
    vst1.8              {q0}, [r0], r1
    vst1.8              {q0}, [r0], r1
    vst1.8              {q0}, [r0], r1
    vst1.8              {q0}, [r0], r1
    vst1.8              {q0}, [r0], r1
    vst1.8              {q0}, [r0], r1
    vst1.8              {q0}, [r0], r1
    vst1.8              {q0}, [r0], r1
    vst1.8              {q0}, [r0], r1
    vst1.8              {q0}, [r0], r1
    vst1.8              {q0}, [r0], r1
    vst1.8              {q0}, [r0], r1
    bx                  lr
    ENDP                ; |vp9_v_predictor_16x16_neon|

;void vp9_v_predictor_32x32_neon(uint8_t *dst, ptrdiff_t y_stride,
;                                const uint8_t *above,
;                                const uint8_t *left)
; r0  uint8_t *dst
; r1  ptrdiff_t y_stride
; r2  const uint8_t *above
; r3  const uint8_t *left

|vp9_v_predictor_32x32_neon| PROC
    vld1.8              {q0, q1}, [r2]
    mov                 r2, #2
loop_v
    vst1.8              {q0, q1}, [r0], r1
    vst1.8              {q0, q1}, [r0], r1
    vst1.8              {q0, q1}, [r0], r1
    vst1.8              {q0, q1}, [r0], r1
    vst1.8              {q0, q1}, [r0], r1
    vst1.8              {q0, q1}, [r0], r1
    vst1.8              {q0, q1}, [r0], r1
    vst1.8              {q0, q1}, [r0], r1
    vst1.8              {q0, q1}, [r0], r1
    vst1.8              {q0, q1}, [r0], r1
    vst1.8              {q0, q1}, [r0], r1
    vst1.8              {q0, q1}, [r0], r1
    vst1.8              {q0, q1}, [r0], r1
    vst1.8              {q0, q1}, [r0], r1
    vst1.8              {q0, q1}, [r0], r1
    vst1.8              {q0, q1}, [r0], r1
    subs                r2, r2, #1
    bgt                 loop_v
    bx                  lr
    ENDP                ; |vp9_v_predictor_32x32_neon|

;void vp9_h_predictor_4x4_neon(uint8_t *dst, ptrdiff_t y_stride,
;                              const uint8_t *above,
;                              const uint8_t *left)
; r0  uint8_t *dst
; r1  ptrdiff_t y_stride
; r2  const uint8_t *above
; r3  const uint8_t *left

|vp9_h_predictor_4x4_neon| PROC
    vld1.32             {d1[0]}, [r3]
    vdup.8              d0, d1[0]
    vst1.32             {d0[0]}, [r0], r1
    vdup.8              d0, d1[1]
    vst1.32             {d0[0]}, [r0], r1
    vdup.8              d0, d1[2]
    vst1.32             {d0[0]}, [r0], r1
    vdup.8              d0, d1[3]
    vst1.32             {d0[0]}, [r0], r1
    bx                  lr
    ENDP                ; |vp9_h_predictor_4x4_neon|

;void vp9_h_predictor_8x8_neon(uint8_t *dst, ptrdiff_t y_stride,
;                              const uint8_t *above,
;                              const uint8_t *left)
; r0  uint8_t *dst
; r1  ptrdiff_t y_stride
; r2  const uint8_t *above
; r3  const uint8_t *left

|vp9_h_predictor_8x8_neon| PROC
    vld1.64             {d1}, [r3]
    vdup.8              d0, d1[0]
    vst1.64             {d0}, [r0], r1
    vdup.8              d0, d1[1]
    vst1.64             {d0}, [r0], r1
    vdup.8              d0, d1[2]
    vst1.64             {d0}, [r0], r1
    vdup.8              d0, d1[3]
    vst1.64             {d0}, [r0], r1
    vdup.8              d0, d1[4]
    vst1.64             {d0}, [r0], r1
    vdup.8              d0, d1[5]
    vst1.64             {d0}, [r0], r1
    vdup.8              d0, d1[6]
    vst1.64             {d0}, [r0], r1
    vdup.8              d0, d1[7]
    vst1.64             {d0}, [r0], r1
    bx                  lr
    ENDP                ; |vp9_h_predictor_8x8_neon|

;void vp9_h_predictor_16x16_neon(uint8_t *dst, ptrdiff_t y_stride,
;                                const uint8_t *above,
;                                const uint8_t *left)
; r0  uint8_t *dst
; r1  ptrdiff_t y_stride
; r2  const uint8_t *above
; r3  const uint8_t *left

|vp9_h_predictor_16x16_neon| PROC
    vld1.8              {q1}, [r3]
    vdup.8              q0, d2[0]
    vst1.8              {q0}, [r0], r1
    vdup.8              q0, d2[1]
    vst1.8              {q0}, [r0], r1
    vdup.8              q0, d2[2]
    vst1.8              {q0}, [r0], r1
    vdup.8              q0, d2[3]
    vst1.8              {q0}, [r0], r1
    vdup.8              q0, d2[4]
    vst1.8              {q0}, [r0], r1
    vdup.8              q0, d2[5]
    vst1.8              {q0}, [r0], r1
    vdup.8              q0, d2[6]
    vst1.8              {q0}, [r0], r1
    vdup.8              q0, d2[7]
    vst1.8              {q0}, [r0], r1
    vdup.8              q0, d3[0]
    vst1.8              {q0}, [r0], r1
    vdup.8              q0, d3[1]
    vst1.8              {q0}, [r0], r1
    vdup.8              q0, d3[2]
    vst1.8              {q0}, [r0], r1
    vdup.8              q0, d3[3]
    vst1.8              {q0}, [r0], r1
    vdup.8              q0, d3[4]
    vst1.8              {q0}, [r0], r1
    vdup.8              q0, d3[5]
    vst1.8              {q0}, [r0], r1
    vdup.8              q0, d3[6]
    vst1.8              {q0}, [r0], r1
    vdup.8              q0, d3[7]
    vst1.8              {q0}, [r0], r1
    bx                  lr
    ENDP                ; |vp9_h_predictor_16x16_neon|

;void vp9_h_predictor_32x32_neon(uint8_t *dst, ptrdiff_t y_stride,
;                                const uint8_t *above,
;                                const uint8_t *left)
; r0  uint8_t *dst
; r1  ptrdiff_t y_stride
; r2  const uint8_t *above
; r3  const uint8_t *left

|vp9_h_predictor_32x32_neon| PROC
    sub                 r1, r1, #16
    mov                 r2, #2
loop_h
    vld1.8              {q1}, [r3]!
    vdup.8              q0, d2[0]
    vst1.8              {q0}, [r0]!
    vst1.8              {q0}, [r0], r1
    vdup.8              q0, d2[1]
    vst1.8              {q0}, [r0]!
    vst1.8              {q0}, [r0], r1
    vdup.8              q0, d2[2]
    vst1.8              {q0}, [r0]!
    vst1.8              {q0}, [r0], r1
    vdup.8              q0, d2[3]
    vst1.8              {q0}, [r0]!
    vst1.8              {q0}, [r0], r1
    vdup.8              q0, d2[4]
    vst1.8              {q0}, [r0]!
    vst1.8              {q0}, [r0], r1
    vdup.8              q0, d2[5]
    vst1.8              {q0}, [r0]!
    vst1.8              {q0}, [r0], r1
    vdup.8              q0, d2[6]
    vst1.8              {q0}, [r0]!
    vst1.8              {q0}, [r0], r1
    vdup.8              q0, d2[7]
    vst1.8              {q0}, [r0]!
    vst1.8              {q0}, [r0], r1
    vdup.8              q0, d3[0]
    vst1.8              {q0}, [r0]!
    vst1.8              {q0}, [r0], r1
    vdup.8              q0, d3[1]
    vst1.8              {q0}, [r0]!
    vst1.8              {q0}, [r0], r1
    vdup.8              q0, d3[2]
    vst1.8              {q0}, [r0]!
    vst1.8              {q0}, [r0], r1
    vdup.8              q0, d3[3]
    vst1.8              {q0}, [r0]!
    vst1.8              {q0}, [r0], r1
    vdup.8              q0, d3[4]
    vst1.8              {q0}, [r0]!
    vst1.8              {q0}, [r0], r1
    vdup.8              q0, d3[5]
    vst1.8              {q0}, [r0]!
    vst1.8              {q0}, [r0], r1
    vdup.8              q0, d3[6]
    vst1.8              {q0}, [r0]!
    vst1.8              {q0}, [r0], r1
    vdup.8              q0, d3[7]
    vst1.8              {q0}, [r0]!
    vst1.8              {q0}, [r0], r1
    subs                r2, r2, #1
    bgt                 loop_h
    bx                  lr
    ENDP                ; |vp9_h_predictor_32x32_neon|

;void vp9_tm_predictor_4x4_neon (uint8_t *dst, ptrdiff_t y_stride,
;                                const uint8_t *above,
;                                const uint8_t *left)
; r0  uint8_t *dst
; r1  ptrdiff_t y_stride
; r2  const uint8_t *above
; r3  const uint8_t *left

|vp9_tm_predictor_4x4_neon| PROC
    ; Load ytop_left = above[-1];
    sub                 r12, r2, #1
    ldrb                r12, [r12]
    vdup.u8             d0, r12

    ; Load above 4 pixels
    vld1.32             {d2[0]}, [r2]

    ; Compute above - ytop_left
    vsubl.u8            q3, d2, d0

    ; Load left row by row and compute left + (above - ytop_left)
    ; 1st row and 2nd row
    ldrb                r12, [r3], #1
    ldrb                r2, [r3], #1
    vdup.u16            q1, r12
    vdup.u16            q2, r2
    vadd.s16            q1, q1, q3
    vadd.s16            q2, q2, q3
    vqshrun.s16         d0, q1, #0
    vqshrun.s16         d1, q2, #0
    vst1.32             {d0[0]}, [r0], r1
    vst1.32             {d1[0]}, [r0], r1

    ; 3rd row and 4th row
    ldrb                r12, [r3], #1
    ldrb                r2, [r3], #1
    vdup.u16            q1, r12
    vdup.u16            q2, r2
    vadd.s16            q1, q1, q3
    vadd.s16            q2, q2, q3
    vqshrun.s16         d0, q1, #0
    vqshrun.s16         d1, q2, #0
    vst1.32             {d0[0]}, [r0], r1
    vst1.32             {d1[0]}, [r0], r1

    bx                  lr
    ENDP                ; |vp9_tm_predictor_4x4_neon|

    END
