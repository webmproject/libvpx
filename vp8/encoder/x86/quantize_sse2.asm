;
;  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
;
;  Use of this source code is governed by a BSD-style license and patent
;  grant that can be found in the LICENSE file in the root of the source
;  tree. All contributing project authors may be found in the AUTHORS
;  file in the root of the source tree.
;


%include "vpx_ports/x86_abi_support.asm"


;int vp8_regular_quantize_b_impl_sse2(
;               short *coeff_ptr,
;               short *zbin_ptr,
;               short *qcoeff_ptr,
;               short *dequant_ptr,
;               const int *default_zig_zag,
;               short *round_ptr,
;               short *quant_ptr,
;               short *dqcoeff_ptr,
;               unsigned short zbin_oq_value,
;               short *zbin_boost_ptr,
;               short *quant_shift);
;
global sym(vp8_regular_quantize_b_impl_sse2)
sym(vp8_regular_quantize_b_impl_sse2):
    push        rbp
    mov         rbp, rsp
    SHADOW_ARGS_TO_STACK 11
    SAVE_XMM
    push        rsi
    push        rdi
    push        rbx
    ALIGN_STACK 16, rax
    %define abs_minus_zbin    0
    %define temp_qcoeff       32
    %define qcoeff            64
    %define eob_tmp           96
    %define stack_size        112
    sub         rsp, stack_size
    ; end prolog

    mov         rdx, arg(0)                 ; coeff_ptr
    mov         rcx, arg(1)                 ; zbin_ptr
    movd        xmm7, arg(8)                ; zbin_oq_value
    mov         rdi, arg(5)                 ; round_ptr
    mov         rsi, arg(6)                 ; quant_ptr

    ; z
    movdqa      xmm0, OWORD PTR[rdx]
    movdqa      xmm4, OWORD PTR[rdx + 16]

    pshuflw     xmm7, xmm7, 0
    punpcklwd   xmm7, xmm7                  ; duplicated zbin_oq_value

    movdqa      xmm1, xmm0
    movdqa      xmm5, xmm4

    ; sz
    psraw       xmm0, 15
    psraw       xmm4, 15

    ; (z ^ sz)
    pxor        xmm1, xmm0
    pxor        xmm5, xmm4

    ; x = abs(z)
    psubw       xmm1, xmm0
    psubw       xmm5, xmm4

    movdqa      xmm2, OWORD PTR[rcx]
    movdqa      xmm3, OWORD PTR[rcx + 16]

    ; *zbin_ptr + zbin_oq_value
    paddw       xmm2, xmm7
    paddw       xmm3, xmm7

    ; x - (*zbin_ptr + zbin_oq_value)
    psubw       xmm1, xmm2
    psubw       xmm5, xmm3
    movdqa      OWORD PTR[rsp + abs_minus_zbin], xmm1
    movdqa      OWORD PTR[rsp + abs_minus_zbin + 16], xmm5

    ; add (zbin_ptr + zbin_oq_value) back
    paddw       xmm1, xmm2
    paddw       xmm5, xmm3

    movdqa      xmm2, OWORD PTR[rdi]
    movdqa      xmm6, OWORD PTR[rdi + 16]

    movdqa      xmm3, OWORD PTR[rsi]
    movdqa      xmm7, OWORD PTR[rsi + 16]

    ; x + round
    paddw       xmm1, xmm2
    paddw       xmm5, xmm6

    ; y = x * quant_ptr >> 16
    pmulhw      xmm3, xmm1
    pmulhw      xmm7, xmm5

    ; y += x
    paddw       xmm1, xmm3
    paddw       xmm5, xmm7

    movdqa      OWORD PTR[rsp + temp_qcoeff], xmm1
    movdqa      OWORD PTR[rsp + temp_qcoeff + 16], xmm5

    pxor        xmm6, xmm6
    ; zero qcoeff
    movdqa      OWORD PTR[rsp + qcoeff], xmm6
    movdqa      OWORD PTR[rsp + qcoeff + 16], xmm6

    mov         [rsp + eob_tmp], DWORD -1   ; eob
    mov         rsi, arg(9)                 ; zbin_boost_ptr
    mov         rdi, arg(4)                 ; default_zig_zag
    mov         rax, arg(10)                ; quant_shift_ptr

%macro ZIGZAG_LOOP 2
rq_zigzag_loop_%1:
    movsxd      rdx, DWORD PTR[rdi + (%1 * 4)] ; rc
    movsx       ebx, WORD PTR [rsi]         ; *zbin_boost_ptr
    lea         rsi, [rsi + 2]              ; zbin_boost_ptr++

    ; x
    movsx       ecx, WORD PTR[rsp + abs_minus_zbin + rdx *2]

    ; if (x >= zbin)
    sub         ecx, ebx                    ; x - zbin
    jl          rq_zigzag_loop_%2           ; x < zbin

    movsx       ebx, WORD PTR[rsp + temp_qcoeff + rdx *2]

    ; downshift by quant_shift[rdx]
    movsx       ecx, WORD PTR[rax + rdx*2]  ; quant_shift_ptr[rc]
    sar         ebx, cl                     ; also sets Z bit
    je          rq_zigzag_loop_%2           ; !y
    mov         WORD PTR[rsp + qcoeff + rdx * 2], bx ;qcoeff_ptr[rc] = temp_qcoeff[rc]

    mov         rsi, arg(9)                 ; reset to b->zrun_zbin_boost
    mov         [rsp + eob_tmp], DWORD %1   ; eob = i
%endmacro
ZIGZAG_LOOP 0, 1
ZIGZAG_LOOP 1, 2
ZIGZAG_LOOP 2, 3
ZIGZAG_LOOP 3, 4
ZIGZAG_LOOP 4, 5
ZIGZAG_LOOP 5, 6
ZIGZAG_LOOP 6, 7
ZIGZAG_LOOP 7, 8
ZIGZAG_LOOP 8, 9
ZIGZAG_LOOP 9, 10
ZIGZAG_LOOP 10, 11
ZIGZAG_LOOP 11, 12
ZIGZAG_LOOP 12, 13
ZIGZAG_LOOP 13, 14
ZIGZAG_LOOP 14, 15
ZIGZAG_LOOP 15, end
rq_zigzag_loop_end:

    mov         rbx, arg(2)                 ; qcoeff_ptr
    mov         rcx, arg(3)                 ; dequant_ptr
    mov         rsi, arg(7)                 ; dqcoeff_ptr
    mov         rax, [rsp + eob_tmp]        ; eob

    movdqa      xmm2, OWORD PTR[rsp + qcoeff]
    movdqa      xmm3, OWORD PTR[rsp + qcoeff + 16]

    ; y ^ sz
    pxor        xmm2, xmm0
    pxor        xmm3, xmm4
    ; x = (y ^ sz) - sz
    psubw       xmm2, xmm0
    psubw       xmm3, xmm4

    movdqa      xmm0, OWORD PTR[rcx]
    movdqa      xmm1, OWORD PTR[rcx + 16]

    pmullw      xmm0, xmm2
    pmullw      xmm1, xmm3

    movdqa      OWORD PTR[rbx], xmm2
    movdqa      OWORD PTR[rbx + 16], xmm3
    movdqa      OWORD PTR[rsi], xmm0        ; store dqcoeff
    movdqa      OWORD PTR[rsi + 16], xmm1   ; store dqcoeff

    add         rax, 1

    ; begin epilog
    add         rsp, stack_size
    pop         rsp
    pop         rbx
    pop         rdi
    pop         rsi
    RESTORE_XMM
    UNSHADOW_ARGS
    pop         rbp
    ret

;int vp8_fast_quantize_b_impl_sse2(short *coeff_ptr,
;                           short *qcoeff_ptr,short *dequant_ptr,
;                           short *inv_scan_order, short *round_ptr,
;                           short *quant_ptr, short *dqcoeff_ptr);
global sym(vp8_fast_quantize_b_impl_sse2)
sym(vp8_fast_quantize_b_impl_sse2):
    push        rbp
    mov         rbp, rsp
    SHADOW_ARGS_TO_STACK 7
    push        rsi
    push        rdi
    ; end prolog

    mov         rdx, arg(0)                 ;coeff_ptr
    mov         rcx, arg(2)                 ;dequant_ptr
    mov         rdi, arg(4)                 ;round_ptr
    mov         rsi, arg(5)                 ;quant_ptr

    movdqa      xmm0, XMMWORD PTR[rdx]
    movdqa      xmm4, XMMWORD PTR[rdx + 16]

    movdqa      xmm2, XMMWORD PTR[rdi]      ;round lo
    movdqa      xmm3, XMMWORD PTR[rdi + 16] ;round hi

    movdqa      xmm1, xmm0
    movdqa      xmm5, xmm4

    psraw       xmm0, 15                    ;sign of z (aka sz)
    psraw       xmm4, 15                    ;sign of z (aka sz)

    pxor        xmm1, xmm0
    pxor        xmm5, xmm4
    psubw       xmm1, xmm0                  ;x = abs(z)
    psubw       xmm5, xmm4                  ;x = abs(z)

    paddw       xmm1, xmm2
    paddw       xmm5, xmm3

    pmulhw      xmm1, XMMWORD PTR[rsi]
    pmulhw      xmm5, XMMWORD PTR[rsi + 16]

    mov         rdi, arg(1)                 ;qcoeff_ptr
    mov         rsi, arg(6)                 ;dqcoeff_ptr

    movdqa      xmm2, XMMWORD PTR[rcx]
    movdqa      xmm3, XMMWORD PTR[rcx + 16]

    pxor        xmm1, xmm0
    pxor        xmm5, xmm4
    psubw       xmm1, xmm0
    psubw       xmm5, xmm4

    movdqa      XMMWORD PTR[rdi], xmm1
    movdqa      XMMWORD PTR[rdi + 16], xmm5

    pmullw      xmm2, xmm1
    pmullw      xmm3, xmm5

    mov         rdi, arg(3)                 ;inv_scan_order

    ; Start with 16
    pxor        xmm4, xmm4                  ;clear all bits
    pcmpeqw     xmm1, xmm4
    pcmpeqw     xmm5, xmm4

    pcmpeqw     xmm4, xmm4                  ;set all bits
    pxor        xmm1, xmm4
    pxor        xmm5, xmm4

    pand        xmm1, XMMWORD PTR[rdi]
    pand        xmm5, XMMWORD PTR[rdi+16]

    pmaxsw      xmm1, xmm5

    ; now down to 8
    pshufd      xmm5, xmm1, 00001110b

    pmaxsw      xmm1, xmm5

    ; only 4 left
    pshuflw     xmm5, xmm1, 00001110b

    pmaxsw      xmm1, xmm5

    ; okay, just 2!
    pshuflw     xmm5, xmm1, 00000001b

    pmaxsw      xmm1, xmm5

    movd        rax, xmm1
    and         rax, 0xff

    movdqa      XMMWORD PTR[rsi], xmm2        ;store dqcoeff
    movdqa      XMMWORD PTR[rsi + 16], xmm3   ;store dqcoeff

    ; begin epilog
    pop         rdi
    pop         rsi
    UNSHADOW_ARGS
    pop         rbp
    ret
