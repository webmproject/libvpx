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

pb_1: times 16 db 1
pw_2: times 8 dw 2
pb_7m1: times 8 db 7, -1
pb_15: times 16 db 15

sh_b01234577: db 0, 1, 2, 3, 4, 5, 7, 7
sh_b12345677: db 1, 2, 3, 4, 5, 6, 7, 7
sh_b23456777: db 2, 3, 4, 5, 6, 7, 7, 7, 0, 0, 0, 0, 0, 0, 0, 0
sh_b0123456777777777: db 0, 1, 2, 3, 4, 5, 6, 7, 7, 7, 7, 7, 7, 7, 7, 7
sh_b1234567777777777: db 1, 2, 3, 4, 5, 6, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7
sh_b2345677777777777: db 2, 3, 4, 5, 6, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7
sh_b2w01234577: db 0, -1, 1, -1, 2, -1, 3, -1, 4, -1, 5, -1, 7, -1, 7, -1
sh_b2w12345677: db 1, -1, 2, -1, 3, -1, 4, -1, 5, -1, 6, -1, 7, -1, 7, -1
sh_b2w23456777: db 2, -1, 3, -1, 4, -1, 5, -1, 6, -1, 7, -1, 7, -1, 7, -1
sh_b2w01234567: db 0, -1, 1, -1, 2, -1, 3, -1, 4, -1, 5, -1, 6, -1, 7, -1
sh_b2w12345678: db 1, -1, 2, -1, 3, -1, 4, -1, 5, -1, 6, -1, 7, -1, 8, -1
sh_b2w23456789: db 2, -1, 3, -1, 4, -1, 5, -1, 6, -1, 7, -1, 8, -1, 9, -1
sh_b2w89abcdef: db 8, -1, 9, -1, 10, -1, 11, -1, 12, -1, 13, -1, 14, -1, 15, -1
sh_b2w9abcdeff: db 9, -1, 10, -1, 11, -1, 12, -1, 13, -1, 14, -1, 15, -1, 15, -1
sh_b2wabcdefff: db 10, -1, 11, -1, 12, -1, 13, -1, 14, -1, 15, -1, 15, -1, 15, -1
sh_b123456789abcdeff: db 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 15
sh_b23456789abcdefff: db 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 15, 15

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

INIT_MMX ssse3
cglobal d45_predictor_4x4, 3, 4, 4, dst, stride, above, goffset
  GET_GOT     goffsetq

  movq                m0, [aboveq]
  pshufb              m2, m0, [GLOBAL(sh_b23456777)]
  pshufb              m1, m0, [GLOBAL(sh_b01234577)]
  pshufb              m0, [GLOBAL(sh_b12345677)]
  pavgb               m3, m2, m1
  pxor                m2, m1
  pand                m2, [GLOBAL(pb_1)]
  psubb               m3, m2
  pavgb               m0, m3

  ; store 4 lines
  movd    [dstq        ], m0
  psrlq               m0, 8
  movd    [dstq+strideq], m0
  lea               dstq, [dstq+strideq*2]
  psrlq               m0, 8
  movd    [dstq        ], m0
  psrlq               m0, 8
  movd    [dstq+strideq], m0

  RESTORE_GOT
  RET

INIT_MMX ssse3
cglobal d45_predictor_8x8, 3, 4, 4, dst, stride, above, goffset
  GET_GOT     goffsetq

  movq                m0, [aboveq]
  mova                m1, [GLOBAL(sh_b12345677)]
  DEFINE_ARGS dst, stride, stride3
  lea           stride3q, [strideq*3]
  pshufb              m2, m0, [GLOBAL(sh_b23456777)]
  pavgb               m3, m2, m0
  pxor                m2, m0
  pshufb              m0, m1
  pand                m2, [GLOBAL(pb_1)]
  psubb               m3, m2
  pavgb               m0, m3

  ; store 4 lines
  movq  [dstq          ], m0
  pshufb              m0, m1
  movq  [dstq+strideq  ], m0
  pshufb              m0, m1
  movq  [dstq+strideq*2], m0
  pshufb              m0, m1
  movq  [dstq+stride3q ], m0
  pshufb              m0, m1
  lea               dstq, [dstq+strideq*4]

  ; store next 4 lines
  movq  [dstq          ], m0
  pshufb              m0, m1
  movq  [dstq+strideq  ], m0
  pshufb              m0, m1
  movq  [dstq+strideq*2], m0
  pshufb              m0, m1
  movq  [dstq+stride3q ], m0

  RESTORE_GOT
  RET

INIT_XMM ssse3
cglobal d45_predictor_16x16, 3, 6, 4, dst, stride, above, dst8, line, goffset
  GET_GOT     goffsetq

  mova                   m0, [aboveq]
  DEFINE_ARGS dst, stride, stride3, dst8, line
  lea              stride3q, [strideq*3]
  lea                 dst8q, [dstq+strideq*8]
  mova                   m1, [GLOBAL(sh_b123456789abcdeff)]
  pshufb                 m2, m0, [GLOBAL(sh_b23456789abcdefff)]
  pavgb                  m3, m2, m0
  pxor                   m2, m0
  pshufb                 m0, m1
  pand                   m2, [GLOBAL(pb_1)]
  psubb                  m3, m2
  pavgb                  m0, m3

  ; first 4 lines and first half of 3rd 4 lines
  mov                 lined, 2
.loop:
  mova   [dstq            ], m0
  movhps [dst8q           ], m0
  pshufb                 m0, m1
  mova   [dstq +strideq   ], m0
  movhps [dst8q+strideq   ], m0
  pshufb                 m0, m1
  mova   [dstq +strideq*2 ], m0
  movhps [dst8q+strideq*2 ], m0
  pshufb                 m0, m1
  mova   [dstq +stride3q  ], m0
  movhps [dst8q+stride3q  ], m0
  pshufb                 m0, m1
  lea                  dstq, [dstq +strideq*4]
  lea                 dst8q, [dst8q+strideq*4]
  dec                 lined
  jnz .loop

  ; bottom-right 8x8 block
  movhps [dstq          +8], m0
  movhps [dstq+strideq  +8], m0
  movhps [dstq+strideq*2+8], m0
  movhps [dstq+stride3q +8], m0
  lea                  dstq, [dstq+strideq*4]
  movhps [dstq          +8], m0
  movhps [dstq+strideq  +8], m0
  movhps [dstq+strideq*2+8], m0
  movhps [dstq+stride3q +8], m0

  RESTORE_GOT
  RET

INIT_XMM ssse3
cglobal d45_predictor_32x32, 3, 6, 7, dst, stride, above, dst16, line, goffset
  GET_GOT     goffsetq

  mova                   m0, [aboveq]
  mova                   m4, [aboveq+16]
  DEFINE_ARGS dst, stride, stride3, dst16, line
  lea              stride3q, [strideq*3]
  lea                dst16q, [dstq  +strideq*8]
  lea                dst16q, [dst16q+strideq*8]
  mova                   m1, [GLOBAL(sh_b123456789abcdeff)]
  pshufb                 m2, m4, [GLOBAL(sh_b23456789abcdefff)]
  pavgb                  m3, m2, m4
  pxor                   m2, m4
  palignr                m5, m4, m0, 1
  palignr                m6, m4, m0, 2
  pshufb                 m4, m1
  pand                   m2, [GLOBAL(pb_1)]
  psubb                  m3, m2
  pavgb                  m4, m3
  pavgb                  m3, m0, m6
  pxor                   m0, m6
  pand                   m0, [GLOBAL(pb_1)]
  psubb                  m3, m0
  pavgb                  m5, m3

  ; write 4x4 lines (and the first half of the second 4x4 lines)
  mov                  lined, 4
.loop:
  mova [dstq               ], m5
  mova [dstq            +16], m4
  mova [dst16q             ], m4
  palignr                 m3, m4, m5, 1
  pshufb                  m4, m1
  mova [dstq  +strideq     ], m3
  mova [dstq  +strideq  +16], m4
  mova [dst16q+strideq     ], m4
  palignr                 m5, m4, m3, 1
  pshufb                  m4, m1
  mova [dstq  +strideq*2   ], m5
  mova [dstq  +strideq*2+16], m4
  mova [dst16q+strideq*2   ], m4
  palignr                 m3, m4, m5, 1
  pshufb                  m4, m1
  mova [dstq  +stride3q    ], m3
  mova [dstq  +stride3q +16], m4
  mova [dst16q+stride3q    ], m4
  palignr                 m5, m4, m3, 1
  pshufb                  m4, m1
  lea                  dstq, [dstq  +strideq*4]
  lea                dst16q, [dst16q+strideq*4]
  dec                 lined
  jnz .loop

  ; write second half of second 4x4 lines
  mova [dstq            +16], m4
  mova [dstq  +strideq  +16], m4
  mova [dstq  +strideq*2+16], m4
  mova [dstq  +stride3q +16], m4
  lea                  dstq, [dstq  +strideq*4]
  mova [dstq            +16], m4
  mova [dstq  +strideq  +16], m4
  mova [dstq  +strideq*2+16], m4
  mova [dstq  +stride3q +16], m4
  lea                  dstq, [dstq  +strideq*4]
  mova [dstq            +16], m4
  mova [dstq  +strideq  +16], m4
  mova [dstq  +strideq*2+16], m4
  mova [dstq  +stride3q +16], m4
  lea                  dstq, [dstq  +strideq*4]
  mova [dstq            +16], m4
  mova [dstq  +strideq  +16], m4
  mova [dstq  +strideq*2+16], m4
  mova [dstq  +stride3q +16], m4

  RESTORE_GOT
  RET
