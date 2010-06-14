;
;  Copyright (c) 2010 The VP8 project authors. All Rights Reserved.
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

    movdqa      DQWORD PTR[rsp + save_xmm6], xmm6
    movdqa      DQWORD PTR[rsp + save_xmm7], xmm7

    mov         rdx, arg(0)                 ;coeff_ptr
    mov         eax, arg(8)                 ;zbin_oq_value

    mov         rcx, arg(1)                 ;zbin_ptr
    movd        xmm7, eax

    movdqa      xmm0, DQWORD PTR[rdx]
    movdqa      xmm4, DQWORD PTR[rdx + 16]

    movdqa      xmm1, xmm0
    movdqa      xmm5, xmm4

    psraw       xmm0, 15                    ;sign of z (aka sz)
    psraw       xmm4, 15                    ;sign of z (aka sz)

    pxor        xmm1, xmm0
    pxor        xmm5, xmm4

    movdqa      xmm2, DQWORD PTR[rcx]       ;load zbin_ptr
    movdqa      xmm3, DQWORD PTR[rcx + 16]  ;load zbin_ptr

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

    movdqa      DQWORD PTR[rsp + abs_minus_zbin_lo], xmm1
    movdqa      DQWORD PTR[rsp + abs_minus_zbin_hi], xmm5

    paddw       xmm1, xmm2                  ;add (zbin_ptr + zbin_oq_value) back
    paddw       xmm5, xmm3                  ;add (zbin_ptr + zbin_oq_value) back

    movdqa      xmm2, DQWORD PTR[rdi]
    movdqa      xmm3, DQWORD PTR[rsi]

    movdqa      xmm6, DQWORD PTR[rdi + 16]
    movdqa      xmm7, DQWORD PTR[rsi + 16]

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

    movdqa      DQWORD PTR[rsp + temp_qcoeff_lo], xmm1
    movdqa      DQWORD PTR[rsp + temp_qcoeff_hi], xmm5

    movdqa      DQWORD PTR[rsi], xmm6       ;zero qcoeff
    movdqa      DQWORD PTR[rsi + 16], xmm6  ;zero qcoeff

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

    movdqa      xmm2, DQWORD PTR[rdi]
    movdqa      xmm3, DQWORD PTR[rdi + 16]

    movdqa      xmm0, DQWORD PTR[rcx]
    movdqa      xmm1, DQWORD PTR[rcx + 16]

    pmullw      xmm0, xmm2
    pmullw      xmm1, xmm3

    movdqa      DQWORD PTR[rsi], xmm0       ;store dqcoeff
    movdqa      DQWORD PTR[rsi + 16], xmm1  ;store dqcoeff

    mov         rax, [rsp + eob]

    movdqa      xmm6, DQWORD PTR[rsp + save_xmm6]
    movdqa      xmm7, DQWORD PTR[rsp + save_xmm7]

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
