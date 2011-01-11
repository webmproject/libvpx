;
;  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
;
;  Use of this source code is governed by a BSD-style license and patent
;  grant that can be found in the LICENSE file in the root of the source
;  tree. All contributing project authors may be found in the AUTHORS
;  file in the root of the source tree.
;


%include "vpx_ports/x86_abi_support.asm"


;int vp8_fast_quantize_b_impl_ssse3(short *coeff_ptr
;               short *qcoeff_ptr,short *dequant_ptr,
;               short *round_ptr,
;               short *quant_ptr, short *dqcoeff_ptr);
;
global sym(vp8_fast_quantize_b_impl_ssse3)
sym(vp8_fast_quantize_b_impl_ssse3):
    push        rbp
    mov         rbp, rsp
    SHADOW_ARGS_TO_STACK 6
    GET_GOT     rbx
    push        rsi
    push        rdi
    ; end prolog

    mov         rdx, arg(0)                 ;coeff_ptr
    mov         rdi, arg(3)                 ;round_ptr
    mov         rsi, arg(4)                 ;quant_ptr

    movdqa      xmm0, [rdx]
    movdqa      xmm4, [rdx + 16]

    movdqa      xmm2, [rdi]                 ;round lo
    movdqa      xmm3, [rdi + 16]            ;round hi

    movdqa      xmm1, xmm0
    movdqa      xmm5, xmm4

    psraw       xmm0, 15                    ;sign of z (aka sz)
    psraw       xmm4, 15                    ;sign of z (aka sz)

    pabsw       xmm1, xmm1
    pabsw       xmm5, xmm5

    paddw       xmm1, xmm2
    paddw       xmm5, xmm3

    pmulhw      xmm1, [rsi]
    pmulhw      xmm5, [rsi + 16]

    mov         rdi, arg(1)                 ;qcoeff_ptr
    mov         rcx, arg(2)                 ;dequant_ptr
    mov         rsi, arg(5)                 ;dqcoeff_ptr

    pxor        xmm1, xmm0
    pxor        xmm5, xmm4
    psubw       xmm1, xmm0
    psubw       xmm5, xmm4

    movdqa      [rdi], xmm1
    movdqa      [rdi + 16], xmm5

    movdqa      xmm2, [rcx]
    movdqa      xmm3, [rcx + 16]

    pxor        xmm4, xmm4
    pmullw      xmm2, xmm1
    pmullw      xmm3, xmm5

    pcmpeqw     xmm1, xmm4                  ;non zero mask
    pcmpeqw     xmm5, xmm4                  ;non zero mask
    packsswb    xmm1, xmm5
    pshufb      xmm1, [ GLOBAL(zz_shuf)]

    pmovmskb    edx, xmm1

;    xor         ecx, ecx
;    mov         eax, -1
;find_eob_loop:
;    shr         edx, 1
;    jc          fq_skip
;    mov         eax, ecx
;fq_skip:
;    inc         ecx
;    cmp         ecx, 16
;    jne         find_eob_loop
    xor         rdi, rdi
    mov         eax, -1
    xor         dx, ax                      ;flip the bits for bsr
    bsr         eax, edx

    movdqa      [rsi], xmm2                 ;store dqcoeff
    movdqa      [rsi + 16], xmm3            ;store dqcoeff

    sub         edi, edx                    ;check for all zeros in bit mask
    sar         edi, 31                     ;0 or -1
    add         eax, 1
    and         eax, edi                    ;if the bit mask was all zero,
                                            ;then eob = 0
    ; begin epilog
    pop         rdi
    pop         rsi
    RESTORE_GOT
    UNSHADOW_ARGS
    pop         rbp
    ret

SECTION_RODATA
align 16
zz_shuf:
    db 0, 1, 4, 8, 5, 2, 3, 6, 9, 12, 13, 10, 7, 11, 14, 15
