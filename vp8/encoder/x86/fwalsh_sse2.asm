;
;  Copyright (c) 2010 The VP8 project authors. All Rights Reserved.
;
;  Use of this source code is governed by a BSD-style license
;  that can be found in the LICENSE file in the root of the source
;  tree. An additional intellectual property rights grant can be found
;  in the file PATENTS.  All contributing project authors may
;  be found in the AUTHORS file in the root of the source tree.
;


%include "vpx_ports/x86_abi_support.asm"

;void vp8_short_walsh4x4_sse2(short *input, short *output, int pitch)
global sym(vp8_short_walsh4x4_sse2)
sym(vp8_short_walsh4x4_sse2):
    push        rbp
    mov         rbp, rsp
    SHADOW_ARGS_TO_STACK 3
    push        rsi
    push        rdi
    ; end prolog

    mov     rsi, arg(0)
    mov     rdi, arg(1)

    movdqu    xmm4, [rsi + 0]       ;ip[4] ip[0]
    movdqu    xmm0, [rsi + 16]      ;ip[12] ip[8]

    pxor  xmm7, xmm7
    ;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ; 13 12 11 10 03 02 01 00
    ;
    ; 33 32 31 30 23 22 21 20
    ;
    movdqa    xmm3, xmm4          ; 13 12 11 10 03 02 01 00
    punpcklwd xmm4, xmm0          ; 23 03 22 02 21 01 20 00
    punpckhwd xmm3, xmm0          ; 33 13 32 12 31 11 30 10
    movdqa    xmm1, xmm4          ; 23 03 22 02 21 01 20 00
    punpcklwd xmm4, xmm3          ; 31 21 11 01 30 20 10 00
    punpckhwd xmm1, xmm3          ; 33 23 13 03 32 22 12 02
    ;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    pshufd    xmm2, xmm1, 4eh       ;ip[8] ip[12]
    movdqa    xmm3, xmm4          ;ip[4] ip[0]

    paddw   xmm4, xmm2          ;ip[4]+ip[8] ip[0]+ip[12] aka b1 a1
    psubw   xmm3, xmm2          ;ip[4]-ip[8] ip[0]-ip[12] aka c1 d1

    movdqa    xmm5, xmm4
    punpcklqdq  xmm4, xmm3          ;d1 a1
    punpckhqdq  xmm5, xmm3          ;c1 b1

    movdqa    xmm1, xmm5          ;c1 b1
    paddw   xmm5, xmm4          ;dl+cl a1+b1 aka op[4] op[0]
    psubw   xmm4, xmm1          ;d1-c1 a1-b1 aka op[12] op[8]
    ;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ; 13 12 11 10 03 02 01 00
    ;
    ; 33 32 31 30 23 22 21 20
    ;
    movdqa    xmm0, xmm5          ; 13 12 11 10 03 02 01 00
    punpcklwd xmm5, xmm4          ; 23 03 22 02 21 01 20 00
    punpckhwd xmm0, xmm4          ; 33 13 32 12 31 11 30 10
    movdqa    xmm1, xmm5          ; 23 03 22 02 21 01 20 00
    punpcklwd xmm5, xmm0          ; 31 21 11 01 30 20 10 00
    punpckhwd xmm1, xmm0          ; 33 23 13 03 32 22 12 02
    ;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    pshufd    xmm2, xmm1, 4eh       ;ip[8] ip[12]
    movdqa    xmm3, xmm5          ;ip[4] ip[0]

    paddw   xmm5, xmm2          ;ip[4]+ip[8] ip[0]+ip[12] aka b1 a1
    psubw   xmm3, xmm2          ;ip[4]-ip[8] ip[0]-ip[12] aka c1 d1

    movdqa    xmm6, xmm5
    punpcklqdq  xmm5, xmm3          ;d1 a1
    punpckhqdq  xmm6, xmm3          ;c1 b1

    movdqa    xmm1, xmm6          ;c1 b1
    paddw   xmm6, xmm5          ;dl+cl a1+b1 aka op[4] op[0]
    psubw   xmm5, xmm1          ;d1-c1 a1-b1 aka op[12] op[8]

    movdqa    xmm0, xmm6          ;aka b2 a2
    movdqa    xmm1, xmm5          ;aka d2 c2

    pcmpgtw   xmm0, xmm7
    pcmpgtw   xmm1, xmm7

    psrlw   xmm0, 15
    psrlw   xmm1, 15

    paddw   xmm6, xmm0
    paddw   xmm5, xmm1

    psraw   xmm6, 1
    psraw   xmm5, 1

    ;   a2 = a1 + b1;
    ;   b2 = c1 + d1;
    ;   c2 = a1 - b1;
    ;   d2 = d1 - c1;
    ;        a2 += (a2>0);
    ;        b2 += (b2>0);
    ;        c2 += (c2>0);
    ;        d2 += (d2>0);
    ;   op[0] = (a2)>>1;
    ;   op[4] = (b2)>>1;
    ;   op[8] = (c2)>>1;
    ;   op[12]= (d2)>>1;

    movdqu  [rdi + 0], xmm6
    movdqu  [rdi + 16], xmm5

    ; begin epilog
    pop rdi
    pop rsi
    UNSHADOW_ARGS
    pop         rbp
    ret
