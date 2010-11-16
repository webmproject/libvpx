;
;  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
;
;  Use of this source code is governed by a BSD-style license and patent
;  grant that can be found in the LICENSE file in the root of the source
;  tree. All contributing project authors may be found in the AUTHORS
;  file in the root of the source tree.
;


%include "vpx_ports/x86_abi_support.asm"


;int vp8_regular_quantize_b_impl_sse2(short *coeff_ptr, short *zbin_ptr,
;               short *qcoeff_ptr,short *dequant_ptr,
;               const int *default_zig_zag, short *round_ptr,
;               short *quant_ptr, short *dqcoeff_ptr,
;               unsigned short zbin_oq_value,
;               short *zbin_boost_ptr);
;
global sym(vp8_regular_quantize_b_impl_sse2)
sym(vp8_regular_quantize_b_impl_sse2):
    push        rbp
    mov         rbp, rsp
    SHADOW_ARGS_TO_STACK 10
    push        rsi
    push        rdi
    push        rbx
    ; end prolog

    ALIGN_STACK 16, rax

    %define abs_minus_zbin_lo 0
    %define abs_minus_zbin_hi 16
    %define temp_qcoeff_lo 32
    %define temp_qcoeff_hi 48
    %define save_xmm6 64
    %define save_xmm7 80
    %define eob 96

    %define vp8_regularquantizeb_stack_size eob + 16

    sub         rsp, vp8_regularquantizeb_stack_size

    movdqa      OWORD PTR[rsp + save_xmm6], xmm6
    movdqa      OWORD PTR[rsp + save_xmm7], xmm7

    mov         rdx, arg(0)                 ;coeff_ptr
    mov         eax, arg(8)                 ;zbin_oq_value

    mov         rcx, arg(1)                 ;zbin_ptr
    movd        xmm7, eax

    movdqa      xmm0, OWORD PTR[rdx]
    movdqa      xmm4, OWORD PTR[rdx + 16]

    movdqa      xmm1, xmm0
    movdqa      xmm5, xmm4

    psraw       xmm0, 15                    ;sign of z (aka sz)
    psraw       xmm4, 15                    ;sign of z (aka sz)

    pxor        xmm1, xmm0
    pxor        xmm5, xmm4

    movdqa      xmm2, OWORD PTR[rcx]        ;load zbin_ptr
    movdqa      xmm3, OWORD PTR[rcx + 16]   ;load zbin_ptr

    pshuflw     xmm7, xmm7, 0
    psubw       xmm1, xmm0                  ;x = abs(z)

    punpcklwd   xmm7, xmm7                  ;duplicated zbin_oq_value
    psubw       xmm5, xmm4                  ;x = abs(z)

    paddw       xmm2, xmm7
    paddw       xmm3, xmm7

    psubw       xmm1, xmm2                  ;sub (zbin_ptr + zbin_oq_value)
    psubw       xmm5, xmm3                  ;sub (zbin_ptr + zbin_oq_value)

    mov         rdi, arg(5)                 ;round_ptr
    mov         rsi, arg(6)                 ;quant_ptr

    movdqa      OWORD PTR[rsp + abs_minus_zbin_lo], xmm1
    movdqa      OWORD PTR[rsp + abs_minus_zbin_hi], xmm5

    paddw       xmm1, xmm2                  ;add (zbin_ptr + zbin_oq_value) back
    paddw       xmm5, xmm3                  ;add (zbin_ptr + zbin_oq_value) back

    movdqa      xmm2, OWORD PTR[rdi]
    movdqa      xmm3, OWORD PTR[rsi]

    movdqa      xmm6, OWORD PTR[rdi + 16]
    movdqa      xmm7, OWORD PTR[rsi + 16]

    paddw       xmm1, xmm2
    paddw       xmm5, xmm6

    pmulhw      xmm1, xmm3
    pmulhw      xmm5, xmm7

    mov         rsi, arg(2)                 ;qcoeff_ptr
    pxor        xmm6, xmm6

    pxor        xmm1, xmm0
    pxor        xmm5, xmm4

    psubw       xmm1, xmm0
    psubw       xmm5, xmm4

    movdqa      OWORD PTR[rsp + temp_qcoeff_lo], xmm1
    movdqa      OWORD PTR[rsp + temp_qcoeff_hi], xmm5

    movdqa      OWORD PTR[rsi], xmm6        ;zero qcoeff
    movdqa      OWORD PTR[rsi + 16], xmm6   ;zero qcoeff

    xor         rax, rax
    mov         rcx, -1

    mov         [rsp + eob], rcx
    mov         rsi, arg(9)                 ;zbin_boost_ptr

    mov         rbx, arg(4)                 ;default_zig_zag

rq_zigzag_loop:
    movsxd      rcx, DWORD PTR[rbx + rax*4] ;now we have rc
    movsx       edi, WORD PTR [rsi]         ;*zbin_boost_ptr aka zbin
    lea         rsi, [rsi + 2]              ;zbin_boost_ptr++

    movsx       edx, WORD PTR[rsp + abs_minus_zbin_lo + rcx *2]

    sub         edx, edi                    ;x - zbin
    jl          rq_zigzag_1

    mov         rdi, arg(2)                 ;qcoeff_ptr

    movsx       edx, WORD PTR[rsp + temp_qcoeff_lo + rcx *2]

    cmp         edx, 0
    je          rq_zigzag_1

    mov         WORD PTR[rdi + rcx * 2], dx ;qcoeff_ptr[rc] = temp_qcoeff[rc]

    mov         rsi, arg(9)                 ;zbin_boost_ptr
    mov         [rsp + eob], rax            ;eob = i

rq_zigzag_1:
    movsxd      rcx, DWORD PTR[rbx + rax*4 + 4]
    movsx       edi, WORD PTR [rsi]         ;*zbin_boost_ptr aka zbin
    lea         rsi, [rsi + 2]              ;zbin_boost_ptr++

    movsx       edx, WORD PTR[rsp + abs_minus_zbin_lo + rcx *2]
    lea         rax, [rax + 1]

    sub         edx, edi                    ;x - zbin
    jl          rq_zigzag_1a

    mov         rdi, arg(2)                 ;qcoeff_ptr

    movsx       edx, WORD PTR[rsp + temp_qcoeff_lo + rcx *2]

    cmp         edx, 0
    je          rq_zigzag_1a

    mov         WORD PTR[rdi + rcx * 2], dx ;qcoeff_ptr[rc] = temp_qcoeff[rc]

    mov         rsi, arg(9)                 ;zbin_boost_ptr
    mov         [rsp + eob], rax            ;eob = i

rq_zigzag_1a:
    movsxd      rcx, DWORD PTR[rbx + rax*4 + 4]
    movsx       edi, WORD PTR [rsi]         ;*zbin_boost_ptr aka zbin
    lea         rsi, [rsi + 2]              ;zbin_boost_ptr++

    movsx       edx, WORD PTR[rsp + abs_minus_zbin_lo + rcx *2]
    lea         rax, [rax + 1]

    sub         edx, edi                    ;x - zbin
    jl          rq_zigzag_1b

    mov         rdi, arg(2)                 ;qcoeff_ptr

    movsx       edx, WORD PTR[rsp + temp_qcoeff_lo + rcx *2]

    cmp         edx, 0
    je          rq_zigzag_1b

    mov         WORD PTR[rdi + rcx * 2], dx ;qcoeff_ptr[rc] = temp_qcoeff[rc]

    mov         rsi, arg(9)                 ;zbin_boost_ptr
    mov         [rsp + eob], rax            ;eob = i

rq_zigzag_1b:
    movsxd      rcx, DWORD PTR[rbx + rax*4 + 4]
    movsx       edi, WORD PTR [rsi]         ;*zbin_boost_ptr aka zbin
    lea         rsi, [rsi + 2]              ;zbin_boost_ptr++

    movsx       edx, WORD PTR[rsp + abs_minus_zbin_lo + rcx *2]
    lea         rax, [rax + 1]

    sub         edx, edi                    ;x - zbin
    jl          rq_zigzag_1c

    mov         rdi, arg(2)                 ;qcoeff_ptr

    movsx       edx, WORD PTR[rsp + temp_qcoeff_lo + rcx *2]

    cmp         edx, 0
    je          rq_zigzag_1c

    mov         WORD PTR[rdi + rcx * 2], dx ;qcoeff_ptr[rc] = temp_qcoeff[rc]

    mov         rsi, arg(9)                 ;zbin_boost_ptr
    mov         [rsp + eob], rax            ;eob = i

rq_zigzag_1c:
    lea         rax, [rax + 1]

    cmp         rax, 16
    jl          rq_zigzag_loop

    mov         rdi, arg(2)                 ;qcoeff_ptr
    mov         rcx, arg(3)                 ;dequant_ptr
    mov         rsi, arg(7)                 ;dqcoeff_ptr

    movdqa      xmm2, OWORD PTR[rdi]
    movdqa      xmm3, OWORD PTR[rdi + 16]

    movdqa      xmm0, OWORD PTR[rcx]
    movdqa      xmm1, OWORD PTR[rcx + 16]

    pmullw      xmm0, xmm2
    pmullw      xmm1, xmm3

    movdqa      OWORD PTR[rsi], xmm0        ;store dqcoeff
    movdqa      OWORD PTR[rsi + 16], xmm1   ;store dqcoeff

    mov         rax, [rsp + eob]

    movdqa      xmm6, OWORD PTR[rsp + save_xmm6]
    movdqa      xmm7, OWORD PTR[rsp + save_xmm7]

    add         rax, 1

    add         rsp, vp8_regularquantizeb_stack_size
    pop         rsp

    ; begin epilog
    pop         rbx
    pop         rdi
    pop         rsi
    UNSHADOW_ARGS
    pop         rbp
    ret


;int vp8_fast_quantize_b_impl_sse2(short *coeff_ptr,
;                           short *qcoeff_ptr,short *dequant_ptr,
;                           short *scan_mask, short *round_ptr,
;                           short *quant_ptr, short *dqcoeff_ptr);
global sym(vp8_fast_quantize_b_impl_sse2)
sym(vp8_fast_quantize_b_impl_sse2):
    push        rbp
    mov         rbp, rsp
    SHADOW_ARGS_TO_STACK 7
    push        rsi
    push        rdi
    push        rbx
    ; end prolog

    ALIGN_STACK 16, rax

    %define save_xmm6  0
    %define save_xmm7 16

    %define vp8_fastquantizeb_stack_size save_xmm7 + 16

    sub         rsp, vp8_fastquantizeb_stack_size

    movdqa      XMMWORD PTR[rsp + save_xmm6], xmm6
    movdqa      XMMWORD PTR[rsp + save_xmm7], xmm7

    mov         rdx, arg(0)                 ;coeff_ptr
    mov         rcx, arg(2)                 ;dequant_ptr
    mov         rax, arg(3)                 ;scan_mask
    mov         rdi, arg(4)                 ;round_ptr
    mov         rsi, arg(5)                 ;quant_ptr

    movdqa      xmm0, XMMWORD PTR[rdx]
    movdqa      xmm4, XMMWORD PTR[rdx + 16]

    movdqa      xmm6, XMMWORD PTR[rdi]      ;round lo
    movdqa      xmm7, XMMWORD PTR[rdi + 16] ;round hi

    movdqa      xmm1, xmm0
    movdqa      xmm5, xmm4

    psraw       xmm0, 15                    ;sign of z (aka sz)
    psraw       xmm4, 15                    ;sign of z (aka sz)

    pxor        xmm1, xmm0
    pxor        xmm5, xmm4
    psubw       xmm1, xmm0                  ;x = abs(z)
    psubw       xmm5, xmm4                  ;x = abs(z)

    paddw       xmm1, xmm6
    paddw       xmm5, xmm7

    pmulhw      xmm1, XMMWORD PTR[rsi]
    pmulhw      xmm5, XMMWORD PTR[rsi + 16]

    mov         rdi, arg(1)                 ;qcoeff_ptr
    mov         rsi, arg(6)                 ;dqcoeff_ptr

    movdqa      xmm6, XMMWORD PTR[rcx]
    movdqa      xmm7, XMMWORD PTR[rcx + 16]

    pxor        xmm1, xmm0
    pxor        xmm5, xmm4
    psubw       xmm1, xmm0
    psubw       xmm5, xmm4

    movdqa      XMMWORD PTR[rdi], xmm1
    movdqa      XMMWORD PTR[rdi + 16], xmm5

    pmullw      xmm6, xmm1
    pmullw      xmm7, xmm5

    movdqa      xmm2, XMMWORD PTR[rax]
    movdqa      xmm3, XMMWORD PTR[rax+16];

    pxor        xmm4, xmm4            ;clear all bits
    pcmpeqw     xmm1, xmm4
    pcmpeqw     xmm5, xmm4

    pcmpeqw     xmm4, xmm4            ;set all bits
    pxor        xmm1, xmm4
    pxor        xmm5, xmm4

    psrlw       xmm1, 15
    psrlw       xmm5, 15

    pmaddwd     xmm1, xmm2
    pmaddwd     xmm5, xmm3

    movq        xmm2, xmm1
    movq        xmm3, xmm5

    psrldq      xmm1, 8
    psrldq      xmm5, 8

    paddd       xmm1, xmm5
    paddd       xmm2, xmm3

    paddd       xmm1, xmm2
    movq        xmm5, xmm1

    psrldq      xmm1, 4
    paddd       xmm5, xmm1

    movq        rcx,  xmm5
    and         rcx,  0xffff

    xor         rdx,  rdx
    sub         rdx,  rcx

    bsr         rax,  rcx
    inc         rax

    sar         rdx,  31
    and         rax,  rdx

    movdqa      XMMWORD PTR[rsi], xmm6        ;store dqcoeff
    movdqa      XMMWORD PTR[rsi + 16], xmm7   ;store dqcoeff

    movdqa      xmm6, XMMWORD PTR[rsp + save_xmm6]
    movdqa      xmm7, XMMWORD PTR[rsp + save_xmm7]

    add         rsp, vp8_fastquantizeb_stack_size
    pop         rsp

    ; begin epilog
    pop         rbx
    pop         rdi
    pop         rsi
    UNSHADOW_ARGS
    pop         rbp
    ret
