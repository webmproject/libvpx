;
;  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
;
;  Use of this source code is governed by a BSD-style license
;  that can be found in the LICENSE file in the root of the source
;  tree. An additional intellectual property rights grant can be found
;  in the file PATENTS.  All contributing project authors may
;  be found in the AUTHORS file in the root of the source tree.
;

%include "third_party/x86inc/x86inc.asm"

SECTION_RODATA
pw_4:  times 8 dw 4
pw_8:  times 8 dw 8
pw_16: times 8 dw 16
pw_32: times 8 dw 32

SECTION .text

INIT_MMX sse
cglobal dc_predictor_4x4, 4, 4, 2, dst, stride, above, left
  pxor                  m1, m1
  movd                  m0, [aboveq]
  punpckldq             m0, [leftq]
  psadbw                m0, m1
  paddw                 m0, [pw_4]
  psraw                 m0, 3
  pshufw                m0, m0, 0x0
  packuswb              m0, m0
  movd      [dstq        ], m0
  movd      [dstq+strideq], m0
  lea                 dstq, [dstq+strideq*2]
  movd      [dstq        ], m0
  movd      [dstq+strideq], m0
  RET

INIT_MMX sse
cglobal dc_predictor_8x8, 4, 4, 3, dst, stride, above, left
  pxor                  m1, m1
  movq                  m0, [aboveq]
  movq                  m2, [leftq]
  DEFINE_ARGS dst, stride, stride3
  lea             stride3q, [strideq*3]
  psadbw                m0, m1
  psadbw                m2, m1
  paddw                 m0, m2
  paddw                 m0, [pw_8]
  psraw                 m0, 4
  pshufw                m0, m0, 0x0
  packuswb              m0, m0
  movq    [dstq          ], m0
  movq    [dstq+strideq  ], m0
  movq    [dstq+strideq*2], m0
  movq    [dstq+stride3q ], m0
  lea                 dstq, [dstq+strideq*4]
  movq    [dstq          ], m0
  movq    [dstq+strideq  ], m0
  movq    [dstq+strideq*2], m0
  movq    [dstq+stride3q ], m0
  RET

INIT_XMM sse2
cglobal dc_predictor_16x16, 4, 4, 3, dst, stride, above, left
  pxor                  m1, m1
  mova                  m0, [aboveq]
  mova                  m2, [leftq]
  DEFINE_ARGS dst, stride, stride3, lines4
  lea             stride3q, [strideq*3]
  mov              lines4d, 4
  psadbw                m0, m1
  psadbw                m2, m1
  paddw                 m0, m2
  movhlps               m2, m0
  paddw                 m0, m2
  paddw                 m0, [pw_16]
  psraw                 m0, 5
  pshuflw               m0, m0, 0x0
  punpcklqdq            m0, m0
  packuswb              m0, m0
.loop:
  mova    [dstq          ], m0
  mova    [dstq+strideq  ], m0
  mova    [dstq+strideq*2], m0
  mova    [dstq+stride3q ], m0
  lea                 dstq, [dstq+strideq*4]
  dec              lines4d
  jnz .loop
  REP_RET

INIT_XMM sse2
cglobal dc_predictor_32x32, 4, 4, 5, dst, stride, above, left
  pxor                  m1, m1
  mova                  m0, [aboveq]
  mova                  m2, [aboveq+16]
  mova                  m3, [leftq]
  mova                  m4, [leftq+16]
  DEFINE_ARGS dst, stride, stride3, lines4
  lea             stride3q, [strideq*3]
  mov              lines4d, 8
  psadbw                m0, m1
  psadbw                m2, m1
  psadbw                m3, m1
  psadbw                m4, m1
  paddw                 m0, m2
  paddw                 m0, m3
  paddw                 m0, m4
  movhlps               m2, m0
  paddw                 m0, m2
  paddw                 m0, [pw_32]
  psraw                 m0, 6
  pshuflw               m0, m0, 0x0
  punpcklqdq            m0, m0
  packuswb              m0, m0
.loop:
  mova [dstq             ], m0
  mova [dstq          +16], m0
  mova [dstq+strideq     ], m0
  mova [dstq+strideq  +16], m0
  mova [dstq+strideq*2   ], m0
  mova [dstq+strideq*2+16], m0
  mova [dstq+stride3q    ], m0
  mova [dstq+stride3q +16], m0
  lea                 dstq, [dstq+strideq*4]
  dec              lines4d
  jnz .loop
  REP_RET

INIT_MMX sse
cglobal v_predictor_4x4, 3, 3, 1, dst, stride, above
  movd                  m0, [aboveq]
  movd      [dstq        ], m0
  movd      [dstq+strideq], m0
  lea                 dstq, [dstq+strideq*2]
  movd      [dstq        ], m0
  movd      [dstq+strideq], m0
  RET

INIT_MMX sse
cglobal v_predictor_8x8, 3, 3, 1, dst, stride, above
  movq                  m0, [aboveq]
  DEFINE_ARGS dst, stride, stride3
  lea             stride3q, [strideq*3]
  movq    [dstq          ], m0
  movq    [dstq+strideq  ], m0
  movq    [dstq+strideq*2], m0
  movq    [dstq+stride3q ], m0
  lea                 dstq, [dstq+strideq*4]
  movq    [dstq          ], m0
  movq    [dstq+strideq  ], m0
  movq    [dstq+strideq*2], m0
  movq    [dstq+stride3q ], m0
  RET

INIT_XMM sse2
cglobal v_predictor_16x16, 3, 4, 1, dst, stride, above
  mova                  m0, [aboveq]
  DEFINE_ARGS dst, stride, stride3, nlines4
  lea             stride3q, [strideq*3]
  mov              nlines4d, 4
.loop:
  mova    [dstq          ], m0
  mova    [dstq+strideq  ], m0
  mova    [dstq+strideq*2], m0
  mova    [dstq+stride3q ], m0
  lea                 dstq, [dstq+strideq*4]
  dec             nlines4d
  jnz .loop
  REP_RET

INIT_XMM sse2
cglobal v_predictor_32x32, 3, 4, 2, dst, stride, above
  mova                  m0, [aboveq]
  mova                  m1, [aboveq+16]
  DEFINE_ARGS dst, stride, stride3, nlines4
  lea             stride3q, [strideq*3]
  mov              nlines4d, 8
.loop:
  mova [dstq             ], m0
  mova [dstq          +16], m1
  mova [dstq+strideq     ], m0
  mova [dstq+strideq  +16], m1
  mova [dstq+strideq*2   ], m0
  mova [dstq+strideq*2+16], m1
  mova [dstq+stride3q    ], m0
  mova [dstq+stride3q +16], m1
  lea                 dstq, [dstq+strideq*4]
  dec             nlines4d
  jnz .loop
  REP_RET
