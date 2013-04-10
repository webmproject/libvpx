;
;  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
;
;  Use of this source code is governed by a BSD-style license
;  that can be found in the LICENSE file in the root of the source
;  tree. An additional intellectual property rights grant can be found
;  in the file PATENTS.  All contributing project authors may
;  be found in the AUTHORS file in the root of the source tree.
;


%include "vpx_ports/x86_abi_support.asm"

;int vp9_block_error_xmm(short *coeff_ptr,  short *dcoef_ptr)
global sym(vp9_block_error_xmm) PRIVATE
sym(vp9_block_error_xmm):
    push        rbp
    mov         rbp, rsp
    SHADOW_ARGS_TO_STACK 2
    push rsi
    push rdi
    ; end prologue

        mov         rsi,        arg(0) ;coeff_ptr
        mov         rdi,        arg(1) ;dcoef_ptr

        movdqa      xmm0,       [rsi]
        movdqa      xmm1,       [rdi]

        movdqa      xmm2,       [rsi+16]
        movdqa      xmm3,       [rdi+16]

        psubw       xmm0,       xmm1
        psubw       xmm2,       xmm3

        pmaddwd     xmm0,       xmm0
        pmaddwd     xmm2,       xmm2

        paddd       xmm0,       xmm2

        pxor        xmm5,       xmm5
        movdqa      xmm1,       xmm0

        punpckldq   xmm0,       xmm5
        punpckhdq   xmm1,       xmm5

        paddd       xmm0,       xmm1
        movdqa      xmm1,       xmm0

        psrldq      xmm0,       8
        paddd       xmm0,       xmm1

        movq        rax,        xmm0

    pop rdi
    pop rsi
    ; begin epilog
    UNSHADOW_ARGS
    pop         rbp
    ret

;int vp9_block_error_mmx(short *coeff_ptr,  short *dcoef_ptr)
global sym(vp9_block_error_mmx) PRIVATE
sym(vp9_block_error_mmx):
    push        rbp
    mov         rbp, rsp
    SHADOW_ARGS_TO_STACK 2
    push rsi
    push rdi
    ; end prolog


        mov         rsi,        arg(0) ;coeff_ptr
        pxor        mm7,        mm7

        mov         rdi,        arg(1) ;dcoef_ptr
        movq        mm3,        [rsi]

        movq        mm4,        [rdi]
        movq        mm5,        [rsi+8]

        movq        mm6,        [rdi+8]
        pxor        mm1,        mm1 ; from movd mm1, dc ; dc =0

        movq        mm2,        mm7
        psubw       mm5,        mm6

        por         mm1,        mm2
        pmaddwd     mm5,        mm5

        pcmpeqw     mm1,        mm7
        psubw       mm3,        mm4

        pand        mm1,        mm3
        pmaddwd     mm1,        mm1

        paddd       mm1,        mm5
        movq        mm3,        [rsi+16]

        movq        mm4,        [rdi+16]
        movq        mm5,        [rsi+24]

        movq        mm6,        [rdi+24]
        psubw       mm5,        mm6

        pmaddwd     mm5,        mm5
        psubw       mm3,        mm4

        pmaddwd     mm3,        mm3
        paddd       mm3,        mm5

        paddd       mm1,        mm3
        movq        mm0,        mm1

        psrlq       mm1,        32
        paddd       mm0,        mm1

        movq        rax,        mm0

    pop rdi
    pop rsi
    ; begin epilog
    UNSHADOW_ARGS
    pop         rbp
    ret
