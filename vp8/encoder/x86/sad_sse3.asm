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

%macro PROCESS_16X2X3 1
%if %1
        movdqa          xmm0,       XMMWORD PTR [rsi]
        lddqu           xmm5,       XMMWORD PTR [rdi]
        lddqu           xmm6,       XMMWORD PTR [rdi+1]
        lddqu           xmm7,       XMMWORD PTR [rdi+2]

        psadbw          xmm5,       xmm0
        psadbw          xmm6,       xmm0
        psadbw          xmm7,       xmm0
%else
        movdqa          xmm0,       XMMWORD PTR [rsi]
        lddqu           xmm1,       XMMWORD PTR [rdi]
        lddqu           xmm2,       XMMWORD PTR [rdi+1]
        lddqu           xmm3,       XMMWORD PTR [rdi+2]

        psadbw          xmm1,       xmm0
        psadbw          xmm2,       xmm0
        psadbw          xmm3,       xmm0

        paddw           xmm5,       xmm1
        paddw           xmm6,       xmm2
        paddw           xmm7,       xmm3
%endif
        movdqa          xmm0,       XMMWORD PTR [rsi+rax]
        lddqu           xmm1,       XMMWORD PTR [rdi+rdx]
        lddqu           xmm2,       XMMWORD PTR [rdi+rdx+1]
        lddqu           xmm3,       XMMWORD PTR [rdi+rdx+2]

        lea             rsi,        [rsi+rax*2]
        lea             rdi,        [rdi+rdx*2]

        psadbw          xmm1,       xmm0
        psadbw          xmm2,       xmm0
        psadbw          xmm3,       xmm0

        paddw           xmm5,       xmm1
        paddw           xmm6,       xmm2
        paddw           xmm7,       xmm3
%endmacro

%macro PROCESS_8X2X3 1
%if %1
        movq            mm0,       QWORD PTR [rsi]
        movq            mm5,       QWORD PTR [rdi]
        movq            mm6,       QWORD PTR [rdi+1]
        movq            mm7,       QWORD PTR [rdi+2]

        psadbw          mm5,       mm0
        psadbw          mm6,       mm0
        psadbw          mm7,       mm0
%else
        movq            mm0,       QWORD PTR [rsi]
        movq            mm1,       QWORD PTR [rdi]
        movq            mm2,       QWORD PTR [rdi+1]
        movq            mm3,       QWORD PTR [rdi+2]

        psadbw          mm1,       mm0
        psadbw          mm2,       mm0
        psadbw          mm3,       mm0

        paddw           mm5,       mm1
        paddw           mm6,       mm2
        paddw           mm7,       mm3
%endif
        movq            mm0,       QWORD PTR [rsi+rax]
        movq            mm1,       QWORD PTR [rdi+rdx]
        movq            mm2,       QWORD PTR [rdi+rdx+1]
        movq            mm3,       QWORD PTR [rdi+rdx+2]

        lea             rsi,       [rsi+rax*2]
        lea             rdi,       [rdi+rdx*2]

        psadbw          mm1,       mm0
        psadbw          mm2,       mm0
        psadbw          mm3,       mm0

        paddw           mm5,       mm1
        paddw           mm6,       mm2
        paddw           mm7,       mm3
%endmacro

%macro LOAD_X4_ADDRESSES 5
        mov             %2,         [%1+REG_SZ_BYTES*0]
        mov             %3,         [%1+REG_SZ_BYTES*1]

        mov             %4,         [%1+REG_SZ_BYTES*2]
        mov             %5,         [%1+REG_SZ_BYTES*3]
%endmacro

%macro PROCESS_16X2X4 1
%if %1
        movdqa          xmm0,       XMMWORD PTR [rsi]
        lddqu           xmm4,       XMMWORD PTR [rcx]
        lddqu           xmm5,       XMMWORD PTR [rdx]
        lddqu           xmm6,       XMMWORD PTR [rbx]
        lddqu           xmm7,       XMMWORD PTR [rdi]

        psadbw          xmm4,       xmm0
        psadbw          xmm5,       xmm0
        psadbw          xmm6,       xmm0
        psadbw          xmm7,       xmm0
%else
        movdqa          xmm0,       XMMWORD PTR [rsi]
        lddqu           xmm1,       XMMWORD PTR [rcx]
        lddqu           xmm2,       XMMWORD PTR [rdx]
        lddqu           xmm3,       XMMWORD PTR [rbx]

        psadbw          xmm1,       xmm0
        psadbw          xmm2,       xmm0
        psadbw          xmm3,       xmm0

        paddw           xmm4,       xmm1
        lddqu           xmm1,       XMMWORD PTR [rdi]
        paddw           xmm5,       xmm2
        paddw           xmm6,       xmm3

        psadbw          xmm1,       xmm0
        paddw           xmm7,       xmm1
%endif
        movdqa          xmm0,       XMMWORD PTR [rsi+rax]
        lddqu           xmm1,       XMMWORD PTR [rcx+rbp]
        lddqu           xmm2,       XMMWORD PTR [rdx+rbp]
        lddqu           xmm3,       XMMWORD PTR [rbx+rbp]

        psadbw          xmm1,       xmm0
        psadbw          xmm2,       xmm0
        psadbw          xmm3,       xmm0

        paddw           xmm4,       xmm1
        lddqu           xmm1,       XMMWORD PTR [rdi+rbp]
        paddw           xmm5,       xmm2
        paddw           xmm6,       xmm3

        lea             rsi,        [rsi+rax*2]
        lea             rcx,        [rcx+rbp*2]

        lea             rdx,        [rdx+rbp*2]
        lea             rbx,        [rbx+rbp*2]

        lea             rdi,        [rdi+rbp*2]

        psadbw          xmm1,       xmm0
        paddw           xmm7,       xmm1

%endmacro

%macro PROCESS_8X2X4 1
%if %1
        movq            mm0,        QWORD PTR [rsi]
        movq            mm4,        QWORD PTR [rcx]
        movq            mm5,        QWORD PTR [rdx]
        movq            mm6,        QWORD PTR [rbx]
        movq            mm7,        QWORD PTR [rdi]

        psadbw          mm4,        mm0
        psadbw          mm5,        mm0
        psadbw          mm6,        mm0
        psadbw          mm7,        mm0
%else
        movq            mm0,        QWORD PTR [rsi]
        movq            mm1,        QWORD PTR [rcx]
        movq            mm2,        QWORD PTR [rdx]
        movq            mm3,        QWORD PTR [rbx]

        psadbw          mm1,        mm0
        psadbw          mm2,        mm0
        psadbw          mm3,        mm0

        paddw           mm4,        mm1
        movq            mm1,        QWORD PTR [rdi]
        paddw           mm5,        mm2
        paddw           mm6,        mm3

        psadbw          mm1,        mm0
        paddw           mm7,        mm1
%endif
        movq            mm0,        QWORD PTR [rsi+rax]
        movq            mm1,        QWORD PTR [rcx+rbp]
        movq            mm2,        QWORD PTR [rdx+rbp]
        movq            mm3,        QWORD PTR [rbx+rbp]

        psadbw          mm1,        mm0
        psadbw          mm2,        mm0
        psadbw          mm3,        mm0

        paddw           mm4,        mm1
        movq            mm1,        QWORD PTR [rdi+rbp]
        paddw           mm5,        mm2
        paddw           mm6,        mm3

        lea             rsi,        [rsi+rax*2]
        lea             rcx,        [rcx+rbp*2]

        lea             rdx,        [rdx+rbp*2]
        lea             rbx,        [rbx+rbp*2]

        lea             rdi,        [rdi+rbp*2]

        psadbw          mm1,        mm0
        paddw           mm7,        mm1

%endmacro

;void int vp8_sad16x16x3_sse3(
;    unsigned char *src_ptr,
;    int  src_stride,
;    unsigned char *ref_ptr,
;    int  ref_stride,
;    int  *results)
global sym(vp8_sad16x16x3_sse3)
sym(vp8_sad16x16x3_sse3):
    push        rbp
    mov         rbp, rsp
    SHADOW_ARGS_TO_STACK 5
    push        rsi
    push        rdi
    ; end prolog

        mov             rsi,        arg(0) ;src_ptr
        mov             rdi,        arg(2) ;ref_ptr

        movsxd          rax,        dword ptr arg(1) ;src_stride
        movsxd          rdx,        dword ptr arg(3) ;ref_stride

        PROCESS_16X2X3 1
        PROCESS_16X2X3 0
        PROCESS_16X2X3 0
        PROCESS_16X2X3 0
        PROCESS_16X2X3 0
        PROCESS_16X2X3 0
        PROCESS_16X2X3 0
        PROCESS_16X2X3 0

        mov             rdi,        arg(4) ;Results

        movq            xmm0,       xmm5
        psrldq          xmm5,       8

        paddw           xmm0,       xmm5
        movd            [rdi],      xmm0
;-
        movq            xmm0,       xmm6
        psrldq          xmm6,       8

        paddw           xmm0,       xmm6
        movd            [rdi+4],    xmm0
;-
        movq            xmm0,       xmm7
        psrldq          xmm7,       8

        paddw           xmm0,       xmm7
        movd            [rdi+8],    xmm0

    ; begin epilog
    pop         rdi
    pop         rsi
    UNSHADOW_ARGS
    pop         rbp
    ret

;void int vp8_sad16x8x3_sse3(
;    unsigned char *src_ptr,
;    int  src_stride,
;    unsigned char *ref_ptr,
;    int  ref_stride,
;    int  *results)
global sym(vp8_sad16x8x3_sse3)
sym(vp8_sad16x8x3_sse3):
    push        rbp
    mov         rbp, rsp
    SHADOW_ARGS_TO_STACK 5
    push        rsi
    push        rdi
    ; end prolog

        mov             rsi,        arg(0) ;src_ptr
        mov             rdi,        arg(2) ;ref_ptr

        movsxd          rax,        dword ptr arg(1) ;src_stride
        movsxd          rdx,        dword ptr arg(3) ;ref_stride

        PROCESS_16X2X3 1
        PROCESS_16X2X3 0
        PROCESS_16X2X3 0
        PROCESS_16X2X3 0

        mov             rdi,        arg(4) ;Results

        movq            xmm0,       xmm5
        psrldq          xmm5,       8

        paddw           xmm0,       xmm5
        movd            [rdi],      xmm0
;-
        movq            xmm0,       xmm6
        psrldq          xmm6,       8

        paddw           xmm0,       xmm6
        movd            [rdi+4],    xmm0
;-
        movq            xmm0,       xmm7
        psrldq          xmm7,       8

        paddw           xmm0,       xmm7
        movd            [rdi+8],    xmm0

    ; begin epilog
    pop         rdi
    pop         rsi
    UNSHADOW_ARGS
    pop         rbp
    ret

;void int vp8_sad8x16x3_sse3(
;    unsigned char *src_ptr,
;    int  src_stride,
;    unsigned char *ref_ptr,
;    int  ref_stride,
;    int  *results)
global sym(vp8_sad8x16x3_sse3)
sym(vp8_sad8x16x3_sse3):
    push        rbp
    mov         rbp, rsp
    SHADOW_ARGS_TO_STACK 5
    push        rsi
    push        rdi
    ; end prolog

        mov             rsi,        arg(0) ;src_ptr
        mov             rdi,        arg(2) ;ref_ptr

        movsxd          rax,        dword ptr arg(1) ;src_stride
        movsxd          rdx,        dword ptr arg(3) ;ref_stride

        PROCESS_8X2X3 1
        PROCESS_8X2X3 0
        PROCESS_8X2X3 0
        PROCESS_8X2X3 0
        PROCESS_8X2X3 0
        PROCESS_8X2X3 0
        PROCESS_8X2X3 0
        PROCESS_8X2X3 0

        mov             rdi,        arg(4) ;Results

        movd            [rdi],      mm5
        movd            [rdi+4],    mm6
        movd            [rdi+8],    mm7

    ; begin epilog
    pop         rdi
    pop         rsi
    UNSHADOW_ARGS
    pop         rbp
    ret

;void int vp8_sad8x8x3_sse3(
;    unsigned char *src_ptr,
;    int  src_stride,
;    unsigned char *ref_ptr,
;    int  ref_stride,
;    int  *results)
global sym(vp8_sad8x8x3_sse3)
sym(vp8_sad8x8x3_sse3):
    push        rbp
    mov         rbp, rsp
    SHADOW_ARGS_TO_STACK 5
    push        rsi
    push        rdi
    ; end prolog

        mov             rsi,        arg(0) ;src_ptr
        mov             rdi,        arg(2) ;ref_ptr

        movsxd          rax,        dword ptr arg(1) ;src_stride
        movsxd          rdx,        dword ptr arg(3) ;ref_stride

        PROCESS_8X2X3 1
        PROCESS_8X2X3 0
        PROCESS_8X2X3 0
        PROCESS_8X2X3 0

        mov             rdi,        arg(4) ;Results

        movd            [rdi],      mm5
        movd            [rdi+4],    mm6
        movd            [rdi+8],    mm7

    ; begin epilog
    pop         rdi
    pop         rsi
    UNSHADOW_ARGS
    pop         rbp
    ret

;void int vp8_sad4x4x3_sse3(
;    unsigned char *src_ptr,
;    int  src_stride,
;    unsigned char *ref_ptr,
;    int  ref_stride,
;    int  *results)
global sym(vp8_sad4x4x3_sse3)
sym(vp8_sad4x4x3_sse3):
    push        rbp
    mov         rbp, rsp
    SHADOW_ARGS_TO_STACK 5
    push        rsi
    push        rdi
    ; end prolog

        mov             rsi,        arg(0) ;src_ptr
        mov             rdi,        arg(2) ;ref_ptr

        movsxd          rax,        dword ptr arg(1) ;src_stride
        movsxd          rdx,        dword ptr arg(3) ;ref_stride

        movd            mm0,        DWORD PTR [rsi]
        movd            mm1,        DWORD PTR [rdi]

        movd            mm2,        DWORD PTR [rsi+rax]
        movd            mm3,        DWORD PTR [rdi+rdx]

        punpcklbw       mm0,        mm2
        punpcklbw       mm1,        mm3

        movd            mm4,        DWORD PTR [rdi+1]
        movd            mm5,        DWORD PTR [rdi+2]

        movd            mm2,        DWORD PTR [rdi+rdx+1]
        movd            mm3,        DWORD PTR [rdi+rdx+2]

        psadbw          mm1,        mm0

        punpcklbw       mm4,        mm2
        punpcklbw       mm5,        mm3

        psadbw          mm4,        mm0
        psadbw          mm5,        mm0



        lea             rsi,        [rsi+rax*2]
        lea             rdi,        [rdi+rdx*2]

        movd            mm0,        DWORD PTR [rsi]
        movd            mm2,        DWORD PTR [rdi]

        movd            mm3,        DWORD PTR [rsi+rax]
        movd            mm6,        DWORD PTR [rdi+rdx]

        punpcklbw       mm0,        mm3
        punpcklbw       mm2,        mm6

        movd            mm3,        DWORD PTR [rdi+1]
        movd            mm7,        DWORD PTR [rdi+2]

        psadbw          mm2,        mm0

        paddw           mm1,        mm2

        movd            mm2,        DWORD PTR [rdi+rdx+1]
        movd            mm6,        DWORD PTR [rdi+rdx+2]

        punpcklbw       mm3,        mm2
        punpcklbw       mm7,        mm6

        psadbw          mm3,        mm0
        psadbw          mm7,        mm0

        paddw           mm3,        mm4
        paddw           mm7,        mm5

        mov             rdi,        arg(4) ;Results
        movd            [rdi],      mm1

        movd            [rdi+4],    mm3
        movd            [rdi+8],    mm7


    ; begin epilog
    pop rdi
    pop rsi
    UNSHADOW_ARGS
    pop         rbp
    ret

;unsigned int vp8_sad16x16_sse3(
;    unsigned char *src_ptr,
;    int  src_stride,
;    unsigned char *ref_ptr,
;    int  ref_stride,
;    int  max_err)
;%define lddqu movdqu
global sym(vp8_sad16x16_sse3)
sym(vp8_sad16x16_sse3):
    push        rbp
    mov         rbp, rsp
    SHADOW_ARGS_TO_STACK 5
    push        rbx
    push        rsi
    push        rdi
    ; end prolog

        mov             rsi,        arg(0) ;src_ptr
        mov             rdi,        arg(2) ;ref_ptr

        movsxd          rbx,        dword ptr arg(1) ;src_stride
        movsxd          rdx,        dword ptr arg(3) ;ref_stride

        lea             rcx,        [rsi+rbx*8]

        lea             rcx,        [rcx+rbx*8]
        pxor            mm7,        mm7

vp8_sad16x16_sse3_loop:

        movq            rax,        mm7
        cmp             rax,        arg(4)
        jg              vp8_sad16x16_early_exit

        movq            mm0,        QWORD PTR [rsi]
        movq            mm2,        QWORD PTR [rsi+8]

        movq            mm1,        QWORD PTR [rdi]
        movq            mm3,        QWORD PTR [rdi+8]

        movq            mm4,        QWORD PTR [rsi+rbx]
        movq            mm5,        QWORD PTR [rdi+rdx]

        psadbw          mm0,        mm1
        psadbw          mm2,        mm3

        movq            mm1,        QWORD PTR [rsi+rbx+8]
        movq            mm3,        QWORD PTR [rdi+rdx+8]

        psadbw          mm4,        mm5
        psadbw          mm1,        mm3

        lea             rsi,        [rsi+rbx*2]
        lea             rdi,        [rdi+rdx*2]

        paddw           mm0,        mm2
        paddw           mm4,        mm1

        paddw           mm7,        mm0
        paddw           mm7,        mm4

        cmp             rsi,        rcx
        jne             vp8_sad16x16_sse3_loop

        movq            rax,        mm7

vp8_sad16x16_early_exit:

    ; begin epilog
    pop         rdi
    pop         rsi
    pop         rbx
    UNSHADOW_ARGS
    pop         rbp
    ret

;void vp8_sad16x16x4d_sse3(
;    unsigned char *src_ptr,
;    int  src_stride,
;    unsigned char *ref_ptr_base,
;    int  ref_stride,
;    int  *results)
global sym(vp8_sad16x16x4d_sse3)
sym(vp8_sad16x16x4d_sse3):
    push        rbp
    mov         rbp, rsp
    SHADOW_ARGS_TO_STACK 5
    push        rsi
    push        rdi
    push        rbx
    ; end prolog

        push            rbp
        mov             rdi,        arg(2) ; ref_ptr_base

        LOAD_X4_ADDRESSES rdi, rcx, rdx, rax, rdi

        mov             rsi,        arg(0) ;src_ptr

        movsxd          rbx,        dword ptr arg(1) ;src_stride
        movsxd          rbp,        dword ptr arg(3) ;ref_stride

        xchg            rbx,        rax

        PROCESS_16X2X4 1
        PROCESS_16X2X4 0
        PROCESS_16X2X4 0
        PROCESS_16X2X4 0
        PROCESS_16X2X4 0
        PROCESS_16X2X4 0
        PROCESS_16X2X4 0
        PROCESS_16X2X4 0

        pop             rbp
        mov             rdi,        arg(4) ;Results

        movq            xmm0,       xmm4
        psrldq          xmm4,       8

        paddw           xmm0,       xmm4
        movd            [rdi],      xmm0
;-
        movq            xmm0,       xmm5
        psrldq          xmm5,       8

        paddw           xmm0,       xmm5
        movd            [rdi+4],    xmm0
;-
        movq            xmm0,       xmm6
        psrldq          xmm6,       8

        paddw           xmm0,       xmm6
        movd            [rdi+8],    xmm0
;-
        movq            xmm0,       xmm7
        psrldq          xmm7,       8

        paddw           xmm0,       xmm7
        movd            [rdi+12],   xmm0

    ; begin epilog
    pop         rbx
    pop         rdi
    pop         rsi
    UNSHADOW_ARGS
    pop         rbp
    ret

;void vp8_sad16x8x4d_sse3(
;    unsigned char *src_ptr,
;    int  src_stride,
;    unsigned char *ref_ptr_base,
;    int  ref_stride,
;    int  *results)
global sym(vp8_sad16x8x4d_sse3)
sym(vp8_sad16x8x4d_sse3):
    push        rbp
    mov         rbp, rsp
    SHADOW_ARGS_TO_STACK 5
    push        rsi
    push        rdi
    push        rbx
    ; end prolog

        push            rbp
        mov             rdi,        arg(2) ; ref_ptr_base

        LOAD_X4_ADDRESSES rdi, rcx, rdx, rax, rdi

        mov             rsi,        arg(0) ;src_ptr

        movsxd          rbx,        dword ptr arg(1) ;src_stride
        movsxd          rbp,        dword ptr arg(3) ;ref_stride

        xchg            rbx,        rax

        PROCESS_16X2X4 1
        PROCESS_16X2X4 0
        PROCESS_16X2X4 0
        PROCESS_16X2X4 0

        pop             rbp
        mov             rdi,        arg(4) ;Results

        movq            xmm0,       xmm4
        psrldq          xmm4,       8

        paddw           xmm0,       xmm4
        movd            [rdi],      xmm0
;-
        movq            xmm0,       xmm5
        psrldq          xmm5,       8

        paddw           xmm0,       xmm5
        movd            [rdi+4],    xmm0
;-
        movq            xmm0,       xmm6
        psrldq          xmm6,       8

        paddw           xmm0,       xmm6
        movd            [rdi+8],    xmm0
;-
        movq            xmm0,       xmm7
        psrldq          xmm7,       8

        paddw           xmm0,       xmm7
        movd            [rdi+12],   xmm0

    ; begin epilog
    pop         rbx
    pop         rdi
    pop         rsi
    UNSHADOW_ARGS
    pop         rbp
    ret

;void int vp8_sad8x16x4d_sse3(
;    unsigned char *src_ptr,
;    int  src_stride,
;    unsigned char *ref_ptr,
;    int  ref_stride,
;    int  *results)
global sym(vp8_sad8x16x4d_sse3)
sym(vp8_sad8x16x4d_sse3):
    push        rbp
    mov         rbp, rsp
    SHADOW_ARGS_TO_STACK 5
    push        rsi
    push        rdi
    push        rbx
    ; end prolog

        push            rbp
        mov             rdi,        arg(2) ; ref_ptr_base

        LOAD_X4_ADDRESSES rdi, rcx, rdx, rax, rdi

        mov             rsi,        arg(0) ;src_ptr

        movsxd          rbx,        dword ptr arg(1) ;src_stride
        movsxd          rbp,        dword ptr arg(3) ;ref_stride

        xchg            rbx,        rax

        PROCESS_8X2X4 1
        PROCESS_8X2X4 0
        PROCESS_8X2X4 0
        PROCESS_8X2X4 0
        PROCESS_8X2X4 0
        PROCESS_8X2X4 0
        PROCESS_8X2X4 0
        PROCESS_8X2X4 0

        pop             rbp
        mov             rdi,        arg(4) ;Results

        movd            [rdi],      mm4
        movd            [rdi+4],    mm5
        movd            [rdi+8],    mm6
        movd            [rdi+12],   mm7

    ; begin epilog
    pop         rbx
    pop         rdi
    pop         rsi
    UNSHADOW_ARGS
    pop         rbp
    ret

;void int vp8_sad8x8x4d_sse3(
;    unsigned char *src_ptr,
;    int  src_stride,
;    unsigned char *ref_ptr,
;    int  ref_stride,
;    int  *results)
global sym(vp8_sad8x8x4d_sse3)
sym(vp8_sad8x8x4d_sse3):
    push        rbp
    mov         rbp, rsp
    SHADOW_ARGS_TO_STACK 5
    push        rsi
    push        rdi
    push        rbx
    ; end prolog

        push            rbp
        mov             rdi,        arg(2) ; ref_ptr_base

        LOAD_X4_ADDRESSES rdi, rcx, rdx, rax, rdi

        mov             rsi,        arg(0) ;src_ptr

        movsxd          rbx,        dword ptr arg(1) ;src_stride
        movsxd          rbp,        dword ptr arg(3) ;ref_stride

        xchg            rbx,        rax

        PROCESS_8X2X4 1
        PROCESS_8X2X4 0
        PROCESS_8X2X4 0
        PROCESS_8X2X4 0

        pop             rbp
        mov             rdi,        arg(4) ;Results

        movd            [rdi],      mm4
        movd            [rdi+4],    mm5
        movd            [rdi+8],    mm6
        movd            [rdi+12],   mm7

    ; begin epilog
    pop         rbx
    pop         rdi
    pop         rsi
    UNSHADOW_ARGS
    pop         rbp
    ret

;void int vp8_sad4x4x4d_sse3(
;    unsigned char *src_ptr,
;    int  src_stride,
;    unsigned char *ref_ptr,
;    int  ref_stride,
;    int  *results)
global sym(vp8_sad4x4x4d_sse3)
sym(vp8_sad4x4x4d_sse3):
    push        rbp
    mov         rbp, rsp
    SHADOW_ARGS_TO_STACK 5
    push        rsi
    push        rdi
    push        rbx
    ; end prolog

        push            rbp
        mov             rdi,        arg(2) ; ref_ptr_base

        LOAD_X4_ADDRESSES rdi, rcx, rdx, rax, rdi

        mov             rsi,        arg(0) ;src_ptr

        movsxd          rbx,        dword ptr arg(1) ;src_stride
        movsxd          rbp,        dword ptr arg(3) ;ref_stride

        xchg            rbx,        rax

        movd            mm0,        DWORD PTR [rsi]
        movd            mm1,        DWORD PTR [rcx]

        movd            mm2,        DWORD PTR [rsi+rax]
        movd            mm3,        DWORD PTR [rcx+rbp]

        punpcklbw       mm0,        mm2
        punpcklbw       mm1,        mm3

        movd            mm4,        DWORD PTR [rdx]
        movd            mm5,        DWORD PTR [rbx]

        movd            mm6,        DWORD PTR [rdi]
        movd            mm2,        DWORD PTR [rdx+rbp]

        movd            mm3,        DWORD PTR [rbx+rbp]
        movd            mm7,        DWORD PTR [rdi+rbp]

        psadbw          mm1,        mm0

        punpcklbw       mm4,        mm2
        punpcklbw       mm5,        mm3

        punpcklbw       mm6,        mm7
        psadbw          mm4,        mm0

        psadbw          mm5,        mm0
        psadbw          mm6,        mm0



        lea             rsi,        [rsi+rax*2]
        lea             rcx,        [rcx+rbp*2]

        lea             rdx,        [rdx+rbp*2]
        lea             rbx,        [rbx+rbp*2]

        lea             rdi,        [rdi+rbp*2]

        movd            mm0,        DWORD PTR [rsi]
        movd            mm2,        DWORD PTR [rcx]

        movd            mm3,        DWORD PTR [rsi+rax]
        movd            mm7,        DWORD PTR [rcx+rbp]

        punpcklbw       mm0,        mm3
        punpcklbw       mm2,        mm7

        movd            mm3,        DWORD PTR [rdx]
        movd            mm7,        DWORD PTR [rbx]

        psadbw          mm2,        mm0
        mov             rax,        rbp

        pop             rbp
        mov             rsi,        arg(4) ;Results

        paddw           mm1,        mm2
        movd            [rsi],      mm1

        movd            mm2,        DWORD PTR [rdx+rax]
        movd            mm1,        DWORD PTR [rbx+rax]

        punpcklbw       mm3,        mm2
        punpcklbw       mm7,        mm1

        psadbw          mm3,        mm0
        psadbw          mm7,        mm0

        movd            mm2,        DWORD PTR [rdi]
        movd            mm1,        DWORD PTR [rdi+rax]

        paddw           mm3,        mm4
        paddw           mm7,        mm5

        movd            [rsi+4],    mm3
        punpcklbw       mm2,        mm1

        movd            [rsi+8],    mm7
        psadbw          mm2,        mm0

        paddw           mm2,        mm6
        movd            [rsi+12],   mm2


    ; begin epilog
    pop         rbx
    pop         rdi
    pop         rsi
    UNSHADOW_ARGS
    pop         rbp
    ret
