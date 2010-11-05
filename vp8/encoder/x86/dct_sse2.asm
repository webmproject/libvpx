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

;void vp8_short_fdct4x4_sse2(short *input, short *output, int pitch)
global sym(vp8_short_fdct4x4_sse2)
sym(vp8_short_fdct4x4_sse2):
    push        rbp
    mov         rbp, rsp
    SHADOW_ARGS_TO_STACK 3
;;    SAVE_XMM
    GET_GOT     rbx
    push        rsi
    push        rdi
    ; end prolog

    mov         rsi, arg(0)
    movsxd      rax, DWORD PTR arg(2)
    lea         rdi, [rsi + rax*2]

    movq        xmm0, MMWORD PTR[rsi   ]        ;03 02 01 00
    movq        xmm2, MMWORD PTR[rsi + rax]     ;13 12 11 10
    movq        xmm1, MMWORD PTR[rsi + rax*2]   ;23 22 21 20
    movq        xmm3, MMWORD PTR[rdi + rax]     ;33 32 31 30

    punpcklqdq  xmm0, xmm2                      ;13 12 11 10 03 02 01 00
    punpcklqdq  xmm1, xmm3                      ;33 32 31 30 23 22 21 20

    mov         rdi, arg(1)

    movdqa      xmm2, xmm0
    punpckldq   xmm0, xmm1                      ;23 22 03 02 21 20 01 00
    punpckhdq   xmm2, xmm1                      ;33 32 13 12 31 30 11 10
    movdqa      xmm1, xmm0
    punpckldq   xmm0, xmm2                      ;31 21 30 20 11 10 01 00
    pshufhw     xmm1, xmm1, 0b1h                ;22 23 02 03 xx xx xx xx
    pshufhw     xmm2, xmm2, 0b1h                ;32 33 12 13 xx xx xx xx

    punpckhdq   xmm1, xmm2                      ;32 33 22 23 12 13 02 03
    movdqa      xmm3, xmm0
    paddw       xmm0, xmm1                      ;b1 a1 b1 a1 b1 a1 b1 a1
    psubw       xmm3, xmm1                      ;c1 d1 c1 d1 c1 d1 c1 d1
    psllw       xmm0, 3                         ;b1 <<= 3 a1 <<= 3
    psllw       xmm3, 3                         ;c1 <<= 3 d1 <<= 3
    movdqa      xmm1, xmm0
    pmaddwd     xmm0, XMMWORD PTR[GLOBAL(_mult_add)]    ;a1 + b1
    pmaddwd     xmm1, XMMWORD PTR[GLOBAL(_mult_sub)]    ;a1 - b1
    movdqa      xmm4, xmm3
    pmaddwd     xmm3, XMMWORD PTR[GLOBAL(_5352_2217)]   ;c1*2217 + d1*5352
    pmaddwd     xmm4, XMMWORD PTR[GLOBAL(_2217_neg5352)];d1*2217 - c1*5352

    paddd       xmm3, XMMWORD PTR[GLOBAL(_14500)]
    paddd       xmm4, XMMWORD PTR[GLOBAL(_7500)]
    psrad       xmm3, 12            ;(c1 * 2217 + d1 * 5352 +  14500)>>12
    psrad       xmm4, 12            ;(d1 * 2217 - c1 * 5352 +   7500)>>12

    packssdw    xmm0, xmm1                      ;op[2] op[0]
    packssdw    xmm3, xmm4                      ;op[3] op[1]
    ; 23 22 21 20 03 02 01 00
    ;
    ; 33 32 31 30 13 12 11 10
    ;
    movdqa      xmm2, xmm0
    punpcklqdq  xmm0, xmm3                      ;13 12 11 10 03 02 01 00
    punpckhqdq  xmm2, xmm3                      ;23 22 21 20 33 32 31 30

    movdqa      xmm3, xmm0
    punpcklwd   xmm0, xmm2                      ;32 30 22 20 12 10 02 00
    punpckhwd   xmm3, xmm2                      ;33 31 23 21 13 11 03 01
    movdqa      xmm2, xmm0
    punpcklwd   xmm0, xmm3                      ;13 12 11 10 03 02 01 00
    punpckhwd   xmm2, xmm3                      ;33 32 31 30 23 22 21 20

    movdqa      xmm5, XMMWORD PTR[GLOBAL(_7)]
    pshufd      xmm2, xmm2, 04eh
    movdqa      xmm3, xmm0
    paddw       xmm0, xmm2                      ;b1 b1 b1 b1 a1 a1 a1 a1
    psubw       xmm3, xmm2                      ;c1 c1 c1 c1 d1 d1 d1 d1

    pshufd      xmm0, xmm0, 0d8h                ;b1 b1 a1 a1 b1 b1 a1 a1
    movdqa      xmm2, xmm3                      ;save d1 for compare
    pshufd      xmm3, xmm3, 0d8h                ;c1 c1 d1 d1 c1 c1 d1 d1
    pshuflw     xmm0, xmm0, 0d8h                ;b1 b1 a1 a1 b1 a1 b1 a1
    pshuflw     xmm3, xmm3, 0d8h                ;c1 c1 d1 d1 c1 d1 c1 d1
    pshufhw     xmm0, xmm0, 0d8h                ;b1 a1 b1 a1 b1 a1 b1 a1
    pshufhw     xmm3, xmm3, 0d8h                ;c1 d1 c1 d1 c1 d1 c1 d1
    movdqa      xmm1, xmm0
    pmaddwd     xmm0, XMMWORD PTR[GLOBAL(_mult_add)] ;a1 + b1
    pmaddwd     xmm1, XMMWORD PTR[GLOBAL(_mult_sub)] ;a1 - b1

    pxor        xmm4, xmm4                      ;zero out for compare
    paddd       xmm0, xmm5
    paddd       xmm1, xmm5
    pcmpeqw     xmm2, xmm4
    psrad       xmm0, 4                         ;(a1 + b1 + 7)>>4
    psrad       xmm1, 4                         ;(a1 - b1 + 7)>>4
    pandn       xmm2, XMMWORD PTR[GLOBAL(_cmp_mask)] ;clear upper,
                                                     ;and keep bit 0 of lower

    movdqa      xmm4, xmm3
    pmaddwd     xmm3, XMMWORD PTR[GLOBAL(_5352_2217)]    ;c1*2217 + d1*5352
    pmaddwd     xmm4, XMMWORD PTR[GLOBAL(_2217_neg5352)] ;d1*2217 - c1*5352
    paddd       xmm3, XMMWORD PTR[GLOBAL(_12000)]
    paddd       xmm4, XMMWORD PTR[GLOBAL(_51000)]
    packssdw    xmm0, xmm1                      ;op[8] op[0]
    psrad       xmm3, 16                ;(c1 * 2217 + d1 * 5352 +  12000)>>16
    psrad       xmm4, 16                ;(d1 * 2217 - c1 * 5352 +  51000)>>16

    packssdw    xmm3, xmm4                      ;op[12] op[4]
    movdqa      xmm1, xmm0
    paddw       xmm3, xmm2                      ;op[4] += (d1!=0)
    punpcklqdq  xmm0, xmm3                      ;op[4] op[0]
    punpckhqdq  xmm1, xmm3                      ;op[12] op[8]

    movdqa      XMMWORD PTR[rdi + 0], xmm0
    movdqa      XMMWORD PTR[rdi + 16], xmm1

    ; begin epilog
    pop rdi
    pop rsi
    RESTORE_GOT
;;    RESTORE_XMM
    UNSHADOW_ARGS
    pop         rbp
    ret

SECTION_RODATA
align 16
_5352_2217:
    dw 5352
    dw 2217
    dw 5352
    dw 2217
    dw 5352
    dw 2217
    dw 5352
    dw 2217
align 16
_2217_neg5352:
    dw 2217
    dw -5352
    dw 2217
    dw -5352
    dw 2217
    dw -5352
    dw 2217
    dw -5352
align 16
_mult_add:
    times 8 dw 1
align 16
_cmp_mask:
    times 4 dw 1
    times 4 dw 0

align 16
_mult_sub:
    dw 1
    dw -1
    dw 1
    dw -1
    dw 1
    dw -1
    dw 1
    dw -1
align 16
_7:
    times 4 dd 7
align 16
_14500:
    times 4 dd 14500
align 16
_7500:
    times 4 dd 7500
align 16
_12000:
    times 4 dd 12000
align 16
_51000:
    times 4 dd 51000
