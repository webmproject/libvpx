;
;  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
;
;  Use of this source code is governed by a BSD-style license and patent
;  grant that can be found in the LICENSE file in the root of the source
;  tree. All contributing project authors may be found in the AUTHORS
;  file in the root of the source tree.
;


%include "vpx_ports/x86_abi_support.asm"
%include "asm_enc_offsets.asm"


; void vp8_regular_quantize_b_sse2 | arg
;  (BLOCK  *b,                     |  0
;   BLOCKD *d)                     |  1

global sym(vp8_regular_quantize_b_sse2)
sym(vp8_regular_quantize_b_sse2):
    push        rbp
    mov         rbp, rsp
    SAVE_XMM
    GET_GOT     rbx
    push        rsi

%if ABI_IS_32BIT
    push        rdi
%else
  %ifidn __OUTPUT_FORMAT__,x64
    push        rdi
  %endif
%endif

    ALIGN_STACK 16, rax
    %define BLOCKD_d          0  ;  8
    %define zrun_zbin_boost   8  ;  8
    %define abs_minus_zbin    16 ; 32
    %define temp_qcoeff       48 ; 32
    %define qcoeff            80 ; 32
    %define stack_size        112
    sub         rsp, stack_size
    ; end prolog

%if ABI_IS_32BIT
    mov         rdi, arg(0)
%else
  %ifidn __OUTPUT_FORMAT__,x64
    mov         rdi, rcx                    ; BLOCK *b
    mov         [rsp + BLOCKD_d], rdx
  %else
    ;mov         rdi, rdi                    ; BLOCK *b
    mov         [rsp + BLOCKD_d], rsi
  %endif
%endif

    mov         rdx, [rdi + vp8_block_coeff] ; coeff_ptr
    mov         rcx, [rdi + vp8_block_zbin] ; zbin_ptr
    movd        xmm7, [rdi + vp8_block_zbin_extra] ; zbin_oq_value

    ; z
    movdqa      xmm0, [rdx]
    movdqa      xmm4, [rdx + 16]
    mov         rdx, [rdi + vp8_block_round] ; round_ptr

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

    movdqa      xmm2, [rcx]
    movdqa      xmm3, [rcx + 16]
    mov         rcx, [rdi + vp8_block_quant] ; quant_ptr

    ; *zbin_ptr + zbin_oq_value
    paddw       xmm2, xmm7
    paddw       xmm3, xmm7

    ; x - (*zbin_ptr + zbin_oq_value)
    psubw       xmm1, xmm2
    psubw       xmm5, xmm3
    movdqa      [rsp + abs_minus_zbin], xmm1
    movdqa      [rsp + abs_minus_zbin + 16], xmm5

    ; add (zbin_ptr + zbin_oq_value) back
    paddw       xmm1, xmm2
    paddw       xmm5, xmm3

    movdqa      xmm2, [rdx]
    movdqa      xmm6, [rdx + 16]

    movdqa      xmm3, [rcx]
    movdqa      xmm7, [rcx + 16]

    ; x + round
    paddw       xmm1, xmm2
    paddw       xmm5, xmm6

    ; y = x * quant_ptr >> 16
    pmulhw      xmm3, xmm1
    pmulhw      xmm7, xmm5

    ; y += x
    paddw       xmm1, xmm3
    paddw       xmm5, xmm7

    movdqa      [rsp + temp_qcoeff], xmm1
    movdqa      [rsp + temp_qcoeff + 16], xmm5

    pxor        xmm6, xmm6
    ; zero qcoeff
    movdqa      [rsp + qcoeff], xmm6
    movdqa      [rsp + qcoeff + 16], xmm6

    mov         rsi, [rdi + vp8_block_zrun_zbin_boost] ; zbin_boost_ptr
    mov         rax, [rdi + vp8_block_quant_shift] ; quant_shift_ptr
    mov         [rsp + zrun_zbin_boost], rsi

%macro ZIGZAG_LOOP 1
    movsx       edx, WORD PTR[GLOBAL(zig_zag + (%1 * 2))] ; rc

    ; x
    movsx       ecx, WORD PTR[rsp + abs_minus_zbin + rdx *2]

    ; if (x >= zbin)
    sub         cx, WORD PTR[rsi]           ; x - zbin
    lea         rsi, [rsi + 2]              ; zbin_boost_ptr++
    jl          rq_zigzag_loop_%1           ; x < zbin

    movsx       edi, WORD PTR[rsp + temp_qcoeff + rdx *2]

    ; downshift by quant_shift[rdx]
    movsx       ecx, WORD PTR[rax + rdx*2]  ; quant_shift_ptr[rc]
    sar         edi, cl                     ; also sets Z bit
    je          rq_zigzag_loop_%1           ; !y
    mov         WORD PTR[rsp + qcoeff + rdx*2], di ;qcoeff_ptr[rc] = temp_qcoeff[rc]
    mov         rsi, [rsp + zrun_zbin_boost] ; reset to b->zrun_zbin_boost
rq_zigzag_loop_%1:
%endmacro
ZIGZAG_LOOP 0
ZIGZAG_LOOP 1
ZIGZAG_LOOP 2
ZIGZAG_LOOP 3
ZIGZAG_LOOP 4
ZIGZAG_LOOP 5
ZIGZAG_LOOP 6
ZIGZAG_LOOP 7
ZIGZAG_LOOP 8
ZIGZAG_LOOP 9
ZIGZAG_LOOP 10
ZIGZAG_LOOP 11
ZIGZAG_LOOP 12
ZIGZAG_LOOP 13
ZIGZAG_LOOP 14
ZIGZAG_LOOP 15

    movdqa      xmm2, [rsp + qcoeff]
    movdqa      xmm3, [rsp + qcoeff + 16]

%if ABI_IS_32BIT
    mov         rdi, arg(1)
%else
    mov         rdi, [rsp + BLOCKD_d]
%endif

    mov         rcx, [rdi + vp8_blockd_dequant] ; dequant_ptr
    mov         rsi, [rdi + vp8_blockd_dqcoeff] ; dqcoeff_ptr

    ; y ^ sz
    pxor        xmm2, xmm0
    pxor        xmm3, xmm4
    ; x = (y ^ sz) - sz
    psubw       xmm2, xmm0
    psubw       xmm3, xmm4

    ; dequant
    movdqa      xmm0, [rcx]
    movdqa      xmm1, [rcx + 16]

    mov         rcx, [rdi + vp8_blockd_qcoeff] ; qcoeff_ptr

    pmullw      xmm0, xmm2
    pmullw      xmm1, xmm3

    movdqa      [rcx], xmm2        ; store qcoeff
    movdqa      [rcx + 16], xmm3
    movdqa      [rsi], xmm0        ; store dqcoeff
    movdqa      [rsi + 16], xmm1

    ; select the last value (in zig_zag order) for EOB
    pcmpeqw     xmm2, xmm6
    pcmpeqw     xmm3, xmm6
    ; !
    pcmpeqw     xmm6, xmm6
    pxor        xmm2, xmm6
    pxor        xmm3, xmm6
    ; mask inv_zig_zag
    pand        xmm2, [GLOBAL(inv_zig_zag)]
    pand        xmm3, [GLOBAL(inv_zig_zag + 16)]
    ; select the max value
    pmaxsw      xmm2, xmm3
    pshufd      xmm3, xmm2, 00001110b
    pmaxsw      xmm2, xmm3
    pshuflw     xmm3, xmm2, 00001110b
    pmaxsw      xmm2, xmm3
    pshuflw     xmm3, xmm2, 00000001b
    pmaxsw      xmm2, xmm3
    movd        eax, xmm2
    and         eax, 0xff
    mov         [rdi + vp8_blockd_eob], eax

    ; begin epilog
    add         rsp, stack_size
    pop         rsp
%if ABI_IS_32BIT
    pop         rdi
%else
  %ifidn __OUTPUT_FORMAT__,x64
    pop         rdi
  %endif
%endif
    pop         rsi
    RESTORE_GOT
    RESTORE_XMM
    pop         rbp
    ret

; int vp8_fast_quantize_b_impl_sse2 | arg
;  (short *coeff_ptr,               |  0
;   short *qcoeff_ptr,              |  1
;   short *dequant_ptr,             |  2
;   short *inv_scan_order,          |  3
;   short *round_ptr,               |  4
;   short *quant_ptr,               |  5
;   short *dqcoeff_ptr)             |  6

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

SECTION_RODATA
align 16
zig_zag:
  dw 0x0000, 0x0001, 0x0004, 0x0008
  dw 0x0005, 0x0002, 0x0003, 0x0006
  dw 0x0009, 0x000c, 0x000d, 0x000a
  dw 0x0007, 0x000b, 0x000e, 0x000f
inv_zig_zag:
  dw 0x0001, 0x0002, 0x0006, 0x0007
  dw 0x0003, 0x0005, 0x0008, 0x000d
  dw 0x0004, 0x0009, 0x000c, 0x000e
  dw 0x000a, 0x000b, 0x000f, 0x0010
