;
;  Copyright (c) 2010 The VP8 project authors. All Rights Reserved.
;
;  Use of this source code is governed by a BSD-style license and patent
;  grant that can be found in the LICENSE file in the root of the source
;  tree. All contributing project authors may be found in the AUTHORS
;  file in the root of the source tree.
;


%include "vpx_ports/x86_abi_support.asm"

global sym(vp8_short_fdct4x4_wmt)

%define         DCTCONSTANTSBITS         (16)
%define         DCTROUNDINGVALUE         (1<< (DCTCONSTANTSBITS-1))
%define         x_c1                      (60547)          ; cos(pi  /8) * (1<<15)
%define         x_c2                      (46341)          ; cos(pi*2/8) * (1<<15)
%define         x_c3                      (25080)          ; cos(pi*3/8) * (1<<15)

%define _1STSTAGESHIFT           14
%define _2NDSTAGESHIFT           16


;; using matrix multiply
;void vp8_short_fdct4x4_wmt(short *input, short *output)
sym(vp8_short_fdct4x4_wmt):
    push        rbp
    mov         rbp, rsp
    SHADOW_ARGS_TO_STACK 2
    GET_GOT     rbx
    ; end prolog

        mov         rax,        arg(0) ;input
        mov         rcx,        arg(1) ;output

        lea         rdx,        [dct_matrix_sse2 GLOBAL]

        movdqu      xmm0,       [rax   ]
        movdqu      xmm1,       [rax+16]

        ; first column
        movdqa      xmm2,       xmm0
        movdqa      xmm7,       [rdx]

        pmaddwd     xmm2,       xmm7
        movdqa      xmm3,       xmm1

        pmaddwd     xmm3,       xmm7
        movdqa      xmm4,       xmm2

        punpckldq   xmm2,       xmm3
        punpckhdq   xmm4,       xmm3

        movdqa      xmm3,       xmm2
        punpckldq   xmm2,       xmm4

        punpckhdq   xmm3,       xmm4
        paddd       xmm2,       xmm3


        paddd       xmm2,       XMMWORD PTR [dct1st_stage_rounding_sse2 GLOBAL]
        psrad       xmm2,       _1STSTAGESHIFT
        ;second column
        movdqa      xmm3,       xmm0
        pmaddwd     xmm3,       [rdx+16]

        movdqa      xmm4,       xmm1
        pmaddwd     xmm4,       [rdx+16]

        movdqa      xmm5,       xmm3
        punpckldq   xmm3,       xmm4

        punpckhdq   xmm5,       xmm4
        movdqa      xmm4,       xmm3

        punpckldq   xmm3,       xmm5
        punpckhdq   xmm4,       xmm5

        paddd       xmm3,       xmm4
        paddd       xmm3,       XMMWORD PTR [dct1st_stage_rounding_sse2 GLOBAL]


        psrad       xmm3,       _1STSTAGESHIFT
        packssdw    xmm2,       xmm3

        ;third column
        movdqa      xmm3,       xmm0
        pmaddwd     xmm3,       [rdx+32]

        movdqa      xmm4,       xmm1
        pmaddwd     xmm4,       [rdx+32]

        movdqa      xmm5,       xmm3
        punpckldq   xmm3,       xmm4

        punpckhdq   xmm5,       xmm4
        movdqa      xmm4,       xmm3

        punpckldq   xmm3,       xmm5
        punpckhdq   xmm4,       xmm5

        paddd       xmm3,       xmm4
        paddd       xmm3,       XMMWORD PTR [dct1st_stage_rounding_sse2 GLOBAL]

        psrad       xmm3,       _1STSTAGESHIFT

        ;fourth column (this is the last column, so we do not have save the source any more)
        pmaddwd     xmm0,       [rdx+48]
        pmaddwd     xmm1,       [rdx+48]

        movdqa      xmm4,       xmm0
        punpckldq   xmm0,       xmm1

        punpckhdq   xmm4,       xmm1
        movdqa      xmm1,       xmm0

        punpckldq   xmm0,       xmm4
        punpckhdq   xmm1,       xmm4

        paddd       xmm0,       xmm1
        paddd       xmm0,       XMMWORD PTR [dct1st_stage_rounding_sse2 GLOBAL]


        psrad       xmm0,       _1STSTAGESHIFT
        packssdw    xmm3,       xmm0
        ; done with one pass
        ; now start second pass
        movdqa      xmm0,       xmm2
        movdqa      xmm1,       xmm3

        pmaddwd     xmm2,       xmm7
        pmaddwd     xmm3,       xmm7

        movdqa      xmm4,       xmm2
        punpckldq   xmm2,       xmm3

        punpckhdq   xmm4,       xmm3
        movdqa      xmm3,       xmm2

        punpckldq   xmm2,       xmm4
        punpckhdq   xmm3,       xmm4

        paddd       xmm2,       xmm3
        paddd       xmm2,       XMMWORD PTR [dct2nd_stage_rounding_sse2 GLOBAL]

        psrad       xmm2,       _2NDSTAGESHIFT

        ;second column
        movdqa      xmm3,       xmm0
        pmaddwd     xmm3,       [rdx+16]

        movdqa      xmm4,       xmm1
        pmaddwd     xmm4,       [rdx+16]

        movdqa      xmm5,       xmm3
        punpckldq   xmm3,       xmm4

        punpckhdq   xmm5,       xmm4
        movdqa      xmm4,       xmm3

        punpckldq   xmm3,       xmm5
        punpckhdq   xmm4,       xmm5

        paddd       xmm3,       xmm4
        paddd       xmm3,       XMMWORD PTR [dct2nd_stage_rounding_sse2 GLOBAL]

        psrad       xmm3,       _2NDSTAGESHIFT
        packssdw    xmm2,       xmm3

        movdqu      [rcx],      xmm2
        ;third column
        movdqa      xmm3,       xmm0
        pmaddwd     xmm3,       [rdx+32]

        movdqa      xmm4,       xmm1
        pmaddwd     xmm4,       [rdx+32]

        movdqa      xmm5,       xmm3
        punpckldq   xmm3,       xmm4

        punpckhdq   xmm5,       xmm4
        movdqa      xmm4,       xmm3

        punpckldq   xmm3,       xmm5
        punpckhdq   xmm4,       xmm5

        paddd       xmm3,       xmm4
        paddd       xmm3,       XMMWORD PTR [dct2nd_stage_rounding_sse2 GLOBAL]

        psrad       xmm3,       _2NDSTAGESHIFT
        ;fourth column
        pmaddwd     xmm0,       [rdx+48]
        pmaddwd     xmm1,       [rdx+48]

        movdqa      xmm4,       xmm0
        punpckldq   xmm0,       xmm1

        punpckhdq   xmm4,       xmm1
        movdqa      xmm1,       xmm0

        punpckldq   xmm0,       xmm4
        punpckhdq   xmm1,       xmm4

        paddd       xmm0,       xmm1
        paddd       xmm0,       XMMWORD PTR [dct2nd_stage_rounding_sse2 GLOBAL]

        psrad       xmm0,       _2NDSTAGESHIFT
        packssdw    xmm3,       xmm0

        movdqu     [rcx+16],   xmm3

    mov rsp, rbp
    ; begin epilog
    RESTORE_GOT
    UNSHADOW_ARGS
    pop         rbp
    ret


SECTION_RODATA
;static unsigned int dct1st_stage_rounding_sse2[4] =
align 16
dct1st_stage_rounding_sse2:
    times 4 dd 8192


;static unsigned int dct2nd_stage_rounding_sse2[4] =
align 16
dct2nd_stage_rounding_sse2:
    times 4 dd 32768

;static short dct_matrix_sse2[4][8]=
align 16
dct_matrix_sse2:
    times 8 dw 23170

    dw  30274
    dw  12540
    dw -12540
    dw -30274
    dw  30274
    dw  12540
    dw -12540
    dw -30274

    dw  23170
    times 2 dw -23170
    times 2 dw  23170
    times 2 dw -23170
    dw  23170

    dw  12540
    dw -30274
    dw  30274
    dw -12540
    dw  12540
    dw -30274
    dw  30274
    dw -12540
