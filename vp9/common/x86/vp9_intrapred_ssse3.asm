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

SECTION .text

INIT_MMX ssse3
cglobal h_predictor_4x4, 2, 4, 3, dst, stride, line, left
  movifnidn          leftq, leftmp
  add                leftq, 4
  mov                lineq, -2
  pxor                  m0, m0
.loop:
  movd                  m1, [leftq+lineq*2  ]
  movd                  m2, [leftq+lineq*2+1]
  pshufb                m1, m0
  pshufb                m2, m0
  movd      [dstq        ], m1
  movd      [dstq+strideq], m2
  lea                 dstq, [dstq+strideq*2]
  inc                lineq
  jnz .loop
  REP_RET

INIT_MMX ssse3
cglobal h_predictor_8x8, 2, 4, 3, dst, stride, line, left
  movifnidn          leftq, leftmp
  add                leftq, 8
  mov                lineq, -4
  pxor                  m0, m0
.loop:
  movd                  m1, [leftq+lineq*2  ]
  movd                  m2, [leftq+lineq*2+1]
  pshufb                m1, m0
  pshufb                m2, m0
  movq      [dstq        ], m1
  movq      [dstq+strideq], m2
  lea                 dstq, [dstq+strideq*2]
  inc                lineq
  jnz .loop
  REP_RET

INIT_XMM ssse3
cglobal h_predictor_16x16, 2, 4, 3, dst, stride, line, left
  movifnidn          leftq, leftmp
  add                leftq, 16
  mov                lineq, -8
  pxor                  m0, m0
.loop:
  movd                  m1, [leftq+lineq*2  ]
  movd                  m2, [leftq+lineq*2+1]
  pshufb                m1, m0
  pshufb                m2, m0
  mova      [dstq        ], m1
  mova      [dstq+strideq], m2
  lea                 dstq, [dstq+strideq*2]
  inc                lineq
  jnz .loop
  REP_RET

INIT_XMM ssse3
cglobal h_predictor_32x32, 2, 4, 3, dst, stride, line, left
  movifnidn          leftq, leftmp
  add                leftq, 32
  mov                lineq, -16
  pxor                  m0, m0
.loop:
  movd                  m1, [leftq+lineq*2  ]
  movd                  m2, [leftq+lineq*2+1]
  pshufb                m1, m0
  pshufb                m2, m0
  mova   [dstq           ], m1
  mova   [dstq        +16], m1
  mova   [dstq+strideq   ], m2
  mova   [dstq+strideq+16], m2
  lea                 dstq, [dstq+strideq*2]
  inc                lineq
  jnz .loop
  REP_RET
