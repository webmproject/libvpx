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

section .text
    global sym(vp8_short_fdct4x4_mmx)
    global sym(vp8_short_fdct8x4_wmt)


%define         DCTCONSTANTSBITS         (16)
%define         DCTROUNDINGVALUE         (1<< (DCTCONSTANTSBITS-1))
%define         x_c1                      (60547)          ; cos(pi  /8) * (1<<15)
%define         x_c2                      (46341)          ; cos(pi*2/8) * (1<<15)
%define         x_c3                      (25080)          ; cos(pi*3/8) * (1<<15)


;void vp8_short_fdct4x4_mmx(short *input, short *output, int pitch)
sym(vp8_short_fdct4x4_mmx):
    push        rbp
    mov         rbp, rsp
    SHADOW_ARGS_TO_STACK 3
    GET_GOT     rbx
    push rsi
    push rdi
    ; end prolog
        mov     rsi,    arg(0) ;input
        mov     rdi,    arg(1) ;output

        lea     rdx,    [GLOBAL(dct_const_mmx)]
        movsxd  rax,    dword ptr arg(2) ;pitch

        lea     rcx,    [rsi + rax*2]
        ; read the input data
        movq    mm0,    [rsi]
        movq    mm1,    [rsi + rax    ]

        movq    mm2,    [rcx]
        movq    mm3,    [rcx + rax]
        ; get the constants
        ;shift to left by 1 for prescision
        psllw   mm0,    3
        psllw   mm1,    3

        psllw   mm2,    3
        psllw   mm3,    3

        ; transpose for the second stage
        movq    mm4,    mm0         ; 00 01 02 03
        movq    mm5,    mm2         ; 10 11 12 03

        punpcklwd   mm0,    mm1     ; 00 10 01 11
        punpckhwd   mm4,    mm1     ; 02 12 03 13

        punpcklwd   mm2,    mm3     ; 20 30 21 31
        punpckhwd   mm5,    mm3     ; 22 32 23 33


        movq        mm1,    mm0     ; 00 10 01 11
        punpckldq   mm0,    mm2     ; 00 10 20 30

        punpckhdq   mm1,    mm2     ; 01 11 21 31

        movq        mm2,    mm4     ; 02 12 03 13
        punpckldq   mm2,    mm5     ; 02 12 22 32

        punpckhdq   mm4,    mm5     ; 03 13 23 33
        movq        mm3,    mm4


        ; first stage
        movq    mm5,    mm0
        movq    mm4,    mm1

        paddw   mm0,    mm3         ; a = 0 + 3
        paddw   mm1,    mm2         ; b = 1 + 2

        psubw   mm4,    mm2         ; c = 1 - 2
        psubw   mm5,    mm3         ; d = 0 - 3


        ; output 0 and 2
        movq    mm6,    [rdx +  16] ; c2
        movq    mm2,    mm0         ; a

        paddw   mm0,    mm1         ; a + b
        psubw   mm2,    mm1         ; a - b

        movq    mm1,    mm0         ; a + b
        pmulhw  mm0,    mm6         ; 00 01 02 03

        paddw   mm0,    mm1         ; output 00 01 02 03
        pmulhw  mm6,    mm2         ; 20 21 22 23

        paddw   mm2,    mm6         ; output 20 21 22 23

        ; output 1 and 3
        movq    mm6,    [rdx +  8]  ; c1
        movq    mm7,    [rdx + 24]  ; c3

        movq    mm1,    mm4         ; c
        movq    mm3,    mm5         ; d

        pmulhw  mm1,    mm7         ; c * c3
        pmulhw  mm3,    mm6         ; d * c1

        paddw   mm3,    mm5         ; d * c1 rounded
        paddw   mm1,    mm3         ; output 10 11 12 13

        movq    mm3,    mm4         ; c
        pmulhw  mm5,    mm7         ; d * c3

        pmulhw  mm4,    mm6         ; c * c1
        paddw   mm3,    mm4         ; round c* c1

        psubw   mm5,    mm3         ; output 30 31 32 33
        movq    mm3,    mm5


        ; done with vertical
        ; transpose for the second stage
        movq    mm4,    mm0         ; 00 01 02 03
        movq    mm5,    mm2         ; 10 11 12 03

        punpcklwd   mm0,    mm1     ; 00 10 01 11
        punpckhwd   mm4,    mm1     ; 02 12 03 13

        punpcklwd   mm2,    mm3     ; 20 30 21 31
        punpckhwd   mm5,    mm3     ; 22 32 23 33


        movq        mm1,    mm0     ; 00 10 01 11
        punpckldq   mm0,    mm2     ; 00 10 20 30

        punpckhdq   mm1,    mm2     ; 01 11 21 31

        movq        mm2,    mm4     ; 02 12 03 13
        punpckldq   mm2,    mm5     ; 02 12 22 32

        punpckhdq   mm4,    mm5     ; 03 13 23 33
        movq        mm3,    mm4


        ; first stage
        movq    mm5,    mm0
        movq    mm4,    mm1

        paddw   mm0,    mm3         ; a = 0 + 3
        paddw   mm1,    mm2         ; b = 1 + 2

        psubw   mm4,    mm2         ; c = 1 - 2
        psubw   mm5,    mm3         ; d = 0 - 3


        ; output 0 and 2
        movq    mm6,    [rdx +  16] ; c2
        movq    mm2,    mm0         ; a
        paddw   mm0,    mm1         ; a + b

        psubw   mm2,    mm1         ; a - b

        movq    mm1,    mm0         ; a + b
        pmulhw  mm0,    mm6         ; 00 01 02 03

        paddw   mm0,    mm1         ; output 00 01 02 03
        pmulhw  mm6,    mm2         ; 20 21 22 23

        paddw   mm2,    mm6         ; output 20 21 22 23


        ; output 1 and 3
        movq    mm6,    [rdx +  8]  ; c1
        movq    mm7,    [rdx + 24]  ; c3

        movq    mm1,    mm4         ; c
        movq    mm3,    mm5         ; d

        pmulhw  mm1,    mm7         ; c * c3
        pmulhw  mm3,    mm6         ; d * c1

        paddw   mm3,    mm5         ; d * c1 rounded
        paddw   mm1,    mm3         ; output 10 11 12 13

        movq    mm3,    mm4         ; c
        pmulhw  mm5,    mm7         ; d * c3

        pmulhw  mm4,    mm6         ; c * c1
        paddw   mm3,    mm4         ; round c* c1

        psubw   mm5,    mm3         ; output 30 31 32 33
        movq    mm3,    mm5
        ; done with vertical

        pcmpeqw mm4,    mm4
        pcmpeqw mm5,    mm5
        psrlw   mm4,    15
        psrlw   mm5,    15

        psllw   mm4,    2
        psllw   mm5,    2

        paddw   mm0,    mm4
        paddw   mm1,    mm5
        paddw   mm2,    mm4
        paddw   mm3,    mm5

        psraw   mm0, 3
        psraw   mm1, 3
        psraw   mm2, 3
        psraw   mm3, 3

        movq        [rdi   ],   mm0
        movq        [rdi+ 8],   mm1
        movq        [rdi+16],   mm2
        movq        [rdi+24],   mm3

    ; begin epilog
    pop rdi
    pop rsi
    RESTORE_GOT
    UNSHADOW_ARGS
    pop         rbp
    ret


;void vp8_short_fdct8x4_wmt(short *input, short *output, int pitch)
sym(vp8_short_fdct8x4_wmt):
    push        rbp
    mov         rbp, rsp
    SHADOW_ARGS_TO_STACK 3
    GET_GOT     rbx
    push rsi
    push rdi
    ; end prolog
        mov         rsi,    arg(0) ;input
        mov         rdi,    arg(1) ;output

        lea         rdx,    [GLOBAL(dct_const_xmm)]
        movsxd      rax,    dword ptr arg(2) ;pitch

        lea         rcx,    [rsi + rax*2]
        ; read the input data
        movdqa      xmm0,       [rsi]
        movdqa      xmm2,       [rsi + rax]

        movdqa      xmm4,       [rcx]
        movdqa      xmm3,       [rcx + rax]
        ; get the constants
        ;shift to left by 1 for prescision
        psllw       xmm0,        3
        psllw       xmm2,        3

        psllw       xmm4,        3
        psllw       xmm3,        3

        ; transpose for the second stage
        movdqa      xmm1,       xmm0         ; 00 01 02 03 04 05 06 07
        movdqa      xmm5,       xmm4         ; 20 21 22 23 24 25 26 27

        punpcklwd   xmm0,       xmm2         ; 00 10 01 11 02 12 03 13
        punpckhwd   xmm1,       xmm2         ; 04 14 05 15 06 16 07 17

        punpcklwd   xmm4,       xmm3         ; 20 30 21 31 22 32 23 33
        punpckhwd   xmm5,       xmm3         ; 24 34 25 35 26 36 27 37

        movdqa      xmm2,       xmm0         ; 00 10 01 11 02 12 03 13
        punpckldq   xmm0,       xmm4         ; 00 10 20 30 01 11 21 31

        punpckhdq   xmm2,       xmm4         ; 02 12 22 32 03 13 23 33


        movdqa      xmm4,       xmm1         ; 04 14 05 15 06 16 07 17
        punpckldq   xmm4,       xmm5         ; 04 14 24 34 05 15 25 35

        punpckhdq   xmm1,       xmm5         ; 06 16 26 36 07 17 27 37
        movdqa      xmm3,       xmm2         ; 02 12 22 32 03 13 23 33

        punpckhqdq  xmm3,       xmm1         ; 03 13 23 33 07 17 27 37
        punpcklqdq  xmm2,       xmm1         ; 02 12 22 32 06 16 26 36

        movdqa      xmm1,       xmm0         ; 00 10 20 30 01 11 21 31
        punpcklqdq  xmm0,       xmm4         ; 00 10 20 30 04 14 24 34

        punpckhqdq  xmm1,       xmm4         ; 01 11 21 32 05 15 25 35

        ; xmm0 0
        ; xmm1 1
        ; xmm2 2
        ; xmm3 3

        ; first stage
        movdqa      xmm5,       xmm0
        movdqa      xmm4,       xmm1

        paddw       xmm0,       xmm3         ; a = 0 + 3
        paddw       xmm1,       xmm2         ; b = 1 + 2

        psubw       xmm4,       xmm2         ; c = 1 - 2
        psubw       xmm5,       xmm3         ; d = 0 - 3


        ; output 0 and 2
        movdqa      xmm6,       [rdx +  32] ; c2
        movdqa      xmm2,       xmm0         ; a

        paddw       xmm0,       xmm1         ; a + b
        psubw       xmm2,       xmm1         ; a - b

        movdqa      xmm1,       xmm0         ; a + b
        pmulhw      xmm0,       xmm6         ; 00 01 02 03

        paddw       xmm0,       xmm1         ; output 00 01 02 03
        pmulhw      xmm6,       xmm2         ; 20 21 22 23

        paddw       xmm2,       xmm6         ; output 20 21 22 23

        ; output 1 and 3
        movdqa      xmm6,       [rdx + 16]  ; c1
        movdqa      xmm7,       [rdx + 48]  ; c3

        movdqa      xmm1,       xmm4         ; c
        movdqa      xmm3,       xmm5         ; d

        pmulhw      xmm1,       xmm7         ; c * c3
        pmulhw      xmm3,       xmm6         ; d * c1

        paddw       xmm3,       xmm5         ; d * c1 rounded
        paddw       xmm1,       xmm3         ; output 10 11 12 13

        movdqa      xmm3,       xmm4         ; c
        pmulhw      xmm5,       xmm7         ; d * c3

        pmulhw      xmm4,       xmm6         ; c * c1
        paddw       xmm3,       xmm4         ; round c* c1

        psubw       xmm5,       xmm3         ; output 30 31 32 33
        movdqa      xmm3,       xmm5


        ; done with vertical
        ; transpose for the second stage
        movdqa      xmm4,       xmm2         ; 02 12 22 32 06 16 26 36
        movdqa      xmm2,       xmm1         ; 01 11 21 31 05 15 25 35

        movdqa      xmm1,       xmm0         ; 00 10 20 30 04 14 24 34
        movdqa      xmm5,       xmm4         ; 02 12 22 32 06 16 26 36

        punpcklwd   xmm0,       xmm2         ; 00 01 10 11 20 21 30 31
        punpckhwd   xmm1,       xmm2         ; 04 05 14 15 24 25 34 35

        punpcklwd   xmm4,       xmm3         ; 02 03 12 13 22 23 32 33
        punpckhwd   xmm5,       xmm3         ; 06 07 16 17 26 27 36 37

        movdqa      xmm2,       xmm0         ; 00 01 10 11 20 21 30 31
        punpckldq   xmm0,       xmm4         ; 00 01 02 03 10 11 12 13

        punpckhdq   xmm2,       xmm4         ; 20 21 22 23 30 31 32 33


        movdqa      xmm4,       xmm1         ; 04 05 14 15 24 25 34 35
        punpckldq   xmm4,       xmm5         ; 04 05 06 07 14 15 16 17

        punpckhdq   xmm1,       xmm5         ; 24 25 26 27 34 35 36 37
        movdqa      xmm3,       xmm2         ; 20 21 22 23 30 31 32 33

        punpckhqdq  xmm3,       xmm1         ; 30 31 32 33 34 35 36 37
        punpcklqdq  xmm2,       xmm1         ; 20 21 22 23 24 25 26 27

        movdqa      xmm1,       xmm0         ; 00 01 02 03 10 11 12 13
        punpcklqdq  xmm0,       xmm4         ; 00 01 02 03 04 05 06 07

        punpckhqdq  xmm1,       xmm4         ; 10 11 12 13 14 15 16 17

        ; first stage
        movdqa      xmm5,       xmm0
        movdqa      xmm4,       xmm1

        paddw       xmm0,       xmm3         ; a = 0 + 3
        paddw       xmm1,       xmm2         ; b = 1 + 2

        psubw       xmm4,       xmm2         ; c = 1 - 2
        psubw       xmm5,       xmm3         ; d = 0 - 3


        ; output 0 and 2
        movdqa      xmm6,       [rdx +  32] ; c2
        movdqa      xmm2,       xmm0         ; a

        paddw       xmm0,       xmm1         ; a + b
        psubw       xmm2,       xmm1         ; a - b

        movdqa      xmm1,       xmm0         ; a + b
        pmulhw      xmm0,       xmm6         ; 00 01 02 03

        paddw       xmm0,       xmm1         ; output 00 01 02 03
        pmulhw      xmm6,       xmm2         ; 20 21 22 23

        paddw       xmm2,       xmm6         ; output 20 21 22 23

        ; output 1 and 3
        movdqa      xmm6,       [rdx + 16]  ; c1
        movdqa      xmm7,       [rdx + 48]  ; c3

        movdqa      xmm1,       xmm4         ; c
        movdqa      xmm3,       xmm5         ; d

        pmulhw      xmm1,       xmm7         ; c * c3
        pmulhw      xmm3,       xmm6         ; d * c1

        paddw       xmm3,       xmm5         ; d * c1 rounded
        paddw       xmm1,       xmm3         ; output 10 11 12 13

        movdqa      xmm3,       xmm4         ; c
        pmulhw      xmm5,       xmm7         ; d * c3

        pmulhw      xmm4,       xmm6         ; c * c1
        paddw       xmm3,       xmm4         ; round c* c1

        psubw       xmm5,       xmm3         ; output 30 31 32 33
        movdqa      xmm3,       xmm5
        ; done with vertical


        pcmpeqw     xmm4,       xmm4
        pcmpeqw     xmm5,       xmm5;
        psrlw       xmm4,       15
        psrlw       xmm5,       15

        psllw       xmm4,       2
        psllw       xmm5,       2

        paddw       xmm0,       xmm4
        paddw       xmm1,       xmm5
        paddw       xmm2,       xmm4
        paddw       xmm3,       xmm5

        psraw       xmm0,       3
        psraw       xmm1,       3
        psraw       xmm2,       3
        psraw       xmm3,       3

        movq        QWORD PTR[rdi   ],   xmm0
        movq        QWORD PTR[rdi+ 8],   xmm1
        movq        QWORD PTR[rdi+16],   xmm2
        movq        QWORD PTR[rdi+24],   xmm3

        psrldq      xmm0,       8
        psrldq      xmm1,       8
        psrldq      xmm2,       8
        psrldq      xmm3,       8

        movq        QWORD PTR[rdi+32],   xmm0
        movq        QWORD PTR[rdi+40],   xmm1
        movq        QWORD PTR[rdi+48],   xmm2
        movq        QWORD PTR[rdi+56],   xmm3
    ; begin epilog
    pop rdi
    pop rsi
    RESTORE_GOT
    UNSHADOW_ARGS
    pop         rbp
    ret


SECTION_RODATA
;static const unsigned int dct1st_stage_rounding_mmx[2] =
align 16
dct1st_stage_rounding_mmx:
    times 2 dd 8192


;static const unsigned int dct2nd_stage_rounding_mmx[2] =
align 16
dct2nd_stage_rounding_mmx:
    times 2 dd 32768


;static const short dct_matrix[4][4]=
align 16
dct_matrix:
    times 4 dw 23170

    dw  30274
    dw  12540
    dw -12540
    dw -30274

    dw 23170
    times 2 dw -23170
    dw 23170

    dw  12540
    dw -30274
    dw  30274
    dw -12540


;static const unsigned short dct_const_mmx[4 * 4]=
align 16
dct_const_mmx:
    times 4 dw 0
    times 4 dw 60547
    times 4 dw 46341
    times 4 dw 25080


;static const unsigned short dct_const_xmm[8 * 4]=
align 16
dct_const_xmm:
    times 8 dw 0
    times 8 dw 60547
    times 8 dw 46341
    times 8 dw 25080
