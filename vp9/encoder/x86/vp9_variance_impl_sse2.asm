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

;unsigned int vp9_get_mb_ss_sse2
;(
;    short *src_ptr
;)
global sym(vp9_get_mb_ss_sse2) PRIVATE
sym(vp9_get_mb_ss_sse2):
    push        rbp
    mov         rbp, rsp
    SHADOW_ARGS_TO_STACK 1
    GET_GOT     rbx
    push rsi
    push rdi
    sub         rsp, 16
    ; end prolog


        mov         rax, arg(0) ;[src_ptr]
        mov         rcx, 8
        pxor        xmm4, xmm4

.NEXTROW:
        movdqa      xmm0, [rax]
        movdqa      xmm1, [rax+16]
        movdqa      xmm2, [rax+32]
        movdqa      xmm3, [rax+48]
        pmaddwd     xmm0, xmm0
        pmaddwd     xmm1, xmm1
        pmaddwd     xmm2, xmm2
        pmaddwd     xmm3, xmm3

        paddd       xmm0, xmm1
        paddd       xmm2, xmm3
        paddd       xmm4, xmm0
        paddd       xmm4, xmm2

        add         rax, 0x40
        dec         rcx
        ja          .NEXTROW

        movdqa      xmm3,xmm4
        psrldq      xmm4,8
        paddd       xmm4,xmm3
        movdqa      xmm3,xmm4
        psrldq      xmm4,4
        paddd       xmm4,xmm3
        movq        rax,xmm4


    ; begin epilog
    add rsp, 16
    pop rdi
    pop rsi
    RESTORE_GOT
    UNSHADOW_ARGS
    pop         rbp
    ret
