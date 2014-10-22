;
;  Copyright (c) 2014 The WebM project authors. All Rights Reserved.
;
;  Use of this source code is governed by a BSD-style license
;  that can be found in the LICENSE file in the root of the source
;  tree. An additional intellectual property rights grant can be found
;  in the file PATENTS.  All contributing project authors may
;  be found in the AUTHORS file in the root of the source tree.
;


%include "vpx_ports/x86_abi_support.asm"

;Note: tap3 and tap4 have to be applied and added after other taps to avoid
;overflow.

%macro HIGH_GET_FILTERS_4 0
    mov         rdx, arg(5)                 ;filter ptr
    mov         rcx, 0x00000040

    movdqa      xmm7, [rdx]                 ;load filters
    pshuflw     xmm10, xmm7, 0b              ;k0
    pshuflw     xmm11, xmm7, 01010101b       ;k1
    pshuflw     xmm12, xmm7, 10101010b       ;k2
    pshuflw     xmm13, xmm7, 11111111b       ;k3
    psrldq      xmm7, 8
    pshuflw     xmm4, xmm7, 0b              ;k4
    pshuflw     xmm5, xmm7, 01010101b       ;k5
    pshuflw     xmm6, xmm7, 10101010b       ;k6
    pshuflw     xmm7, xmm7, 11111111b       ;k7

    punpcklwd   xmm10, xmm6
    punpcklwd   xmm12, xmm5
    punpcklwd   xmm13, xmm4
    punpcklwd   xmm11, xmm7

    movq        xmm9, rcx
    pshufd      xmm9, xmm9, 0

    ;Compute max and min values of a pixel
    mov         rdx, 0x00010001
    movsxd      rcx, DWORD PTR arg(6)      ;bps
    movq        xmm14, rdx
    movq        xmm1, rcx
    pshufd      xmm14, xmm14, 0b
    movdqa      xmm2, xmm14
    psllw       xmm14, xmm1
    psubw       xmm14, xmm2                 ;max value (for clamping)
    pxor        xmm8, xmm8                ;min value (for clamping)

%endm

%macro HIGH_APPLY_FILTER_4 1
    punpcklwd   xmm0, xmm6                  ;two row in one register
    punpcklwd   xmm1, xmm7
    punpcklwd   xmm2, xmm5
    punpcklwd   xmm3, xmm4

    pmaddwd     xmm0, xmm10                 ;multiply the filter factors
    pmaddwd     xmm1, xmm11
    pmaddwd     xmm2, xmm12
    pmaddwd     xmm3, xmm13

    paddd       xmm0, xmm1                  ;sum
    paddd       xmm0, xmm2
    paddd       xmm0, xmm3

    paddd       xmm0, xmm9                  ;rounding
    psrad       xmm0, 7                     ;shift
    packssdw    xmm0, xmm0                  ;pack to word

    ;clamp the values
    pminsw      xmm0, xmm14
    pmaxsw      xmm0, xmm8

%if %1
    movq        xmm1, [rdi]
    pavgw       xmm0, xmm1
%endif
    movq        [rdi], xmm0
%endm

%macro HIGH_GET_FILTERS 0
    mov         rdx, arg(5)                 ;filter ptr
    mov         rsi, arg(0)                 ;src_ptr
    mov         rdi, arg(2)                 ;output_ptr
    mov         rcx, 0x00000040

    movdqa      xmm7, [rdx]                 ;load filters
    pshuflw     xmm10, xmm7, 0b             ;k0
    pshuflw     xmm1, xmm7, 01010101b       ;k1
    pshuflw     xmm12, xmm7, 10101010b      ;k2
    pshuflw     xmm13, xmm7, 11111111b      ;k3
    pshufhw     xmm4, xmm7, 0b              ;k4
    pshufhw     xmm5, xmm7, 01010101b       ;k5
    pshufhw     xmm11, xmm7, 10101010b      ;k6
    pshufhw     xmm7, xmm7, 11111111b       ;k7
    punpcklqdq  xmm12, xmm12
    punpcklqdq  xmm13, xmm13
    punpcklwd   xmm10, xmm1
    punpckhwd   xmm11, xmm7
    punpckhwd   xmm12, xmm5
    punpckhwd   xmm13, xmm4

    movq        xmm9, rcx
    pshufd      xmm9, xmm9, 0              ;rounding

    ;Compute max and min values of a pixel
    mov         rdx, 0x00010001
    movsxd      rcx, DWORD PTR arg(6)       ;bps
    movq        xmm14, rdx
    movq        xmm1, rcx
    pshufd      xmm14, xmm14, 0b
    movdqa      xmm2, xmm14
    psllw       xmm14, xmm1
    psubw       xmm14, xmm2                 ;max value (for clamping)
    pxor        xmm15, xmm15                ;min value (for clamping)
%endm

%macro LOAD_VERT_8 1
    movdqu      xmm0, [rsi + %1]            ;0
    movdqu      xmm1, [rsi + rax + %1]      ;1
    movdqu      xmm6, [rsi + rdx * 2 + %1]  ;6
    lea         rsi,  [rsi + rax]
    movdqu      xmm7, [rsi + rdx * 2 + %1]  ;7
    movdqu      xmm2, [rsi + rax + %1]      ;2
    movdqu      xmm3, [rsi + rax * 2 + %1]  ;3
    movdqu      xmm4, [rsi + rdx + %1]      ;4
    movdqu      xmm5, [rsi + rax * 4 + %1]  ;5
%endm

%macro HIGH_APPLY_FILTER_8 2
    movdqa      xmm8, xmm4
    movdqa      xmm4, xmm0
    punpcklwd   xmm0, xmm1
    punpckhwd   xmm4, xmm1
    movdqa      xmm1, xmm6
    punpcklwd   xmm6, xmm7
    punpckhwd   xmm1, xmm7
    movdqa      xmm7, xmm2
    punpcklwd   xmm2, xmm5
    punpckhwd   xmm7, xmm5

    movdqa      xmm5, xmm8
    movdqa      xmm8, xmm4
    movdqa      xmm4, xmm3
    punpcklwd   xmm3, xmm5
    punpckhwd   xmm4, xmm5
    movdqa      xmm5, xmm8

    pmaddwd     xmm0, xmm10
    pmaddwd     xmm5, xmm10
    pmaddwd     xmm6, xmm11
    pmaddwd     xmm1, xmm11
    pmaddwd     xmm2, xmm12
    pmaddwd     xmm7, xmm12
    pmaddwd     xmm3, xmm13
    pmaddwd     xmm4, xmm13

    paddd       xmm0, xmm6
    paddd       xmm0, xmm2
    paddd       xmm0, xmm3
    paddd       xmm5, xmm1
    paddd       xmm5, xmm7
    paddd       xmm5, xmm4

    paddd       xmm0, xmm9                  ;rounding
    paddd       xmm5, xmm9
    psrad       xmm0, 7                     ;shift
    psrad       xmm5, 7
    packssdw    xmm0, xmm5                  ;pack back to word

    ;clamp the values
    pminsw      xmm0, xmm14
    pmaxsw      xmm0, xmm15

%if %1
    movdqu      xmm1, [rdi + %2]
    pavgw       xmm0, xmm1
%endif
    movdqu      [rdi + %2], xmm0
%endm

;void vp9_filter_block1d4_v8_sse2
;(
;    unsigned char *src_ptr,
;    unsigned int   src_pitch,
;    unsigned char *output_ptr,
;    unsigned int   out_pitch,
;    unsigned int   output_height,
;    short *filter
;)
global sym(vp9_highbd_filter_block1d4_v8_sse2) PRIVATE
sym(vp9_highbd_filter_block1d4_v8_sse2):
    push        rbp
    mov         rbp, rsp
    SHADOW_ARGS_TO_STACK 7
    SAVE_XMM 14
    push        rsi
    push        rdi
    push        rbx
    ; end prolog

    HIGH_GET_FILTERS_4

    mov         rsi, arg(0)                 ;src_ptr
    mov         rdi, arg(2)                 ;output_ptr

    movsxd      rax, DWORD PTR arg(1)       ;pixels_per_line
    movsxd      rbx, DWORD PTR arg(3)       ;out_pitch
    lea         rax, [rax + rax]            ;bytes per line
    lea         rbx, [rbx + rbx]
    lea         rdx, [rax + rax * 2]
    movsxd      rcx, DWORD PTR arg(4)       ;output_height

.loop:
    movq        xmm0, [rsi]                 ;load src: row 0
    movq        xmm1, [rsi + rax]           ;1
    movq        xmm6, [rsi + rdx * 2]       ;6
    lea         rsi,  [rsi + rax]
    movq        xmm7, [rsi + rdx * 2]       ;7
    movq        xmm2, [rsi + rax]           ;2
    movq        xmm3, [rsi + rax * 2]       ;3
    movq        xmm4, [rsi + rdx]           ;4
    movq        xmm5, [rsi + rax * 4]       ;5

    HIGH_APPLY_FILTER_4 0

    lea         rdi, [rdi + rbx]
    dec         rcx
    jnz         .loop

    pop rbx
    ; begin epilog
    pop rdi
    pop rsi
    RESTORE_XMM
    UNSHADOW_ARGS
    pop         rbp
    ret

;void vp9_filter_block1d8_v8_sse2
;(
;    unsigned char *src_ptr,
;    unsigned int   src_pitch,
;    unsigned char *output_ptr,
;    unsigned int   out_pitch,
;    unsigned int   output_height,
;    short *filter
;)
global sym(vp9_highbd_filter_block1d8_v8_sse2) PRIVATE
sym(vp9_highbd_filter_block1d8_v8_sse2):
    push        rbp
    mov         rbp, rsp
    SHADOW_ARGS_TO_STACK 7
    SAVE_XMM 15
    push        rsi
    push        rdi
    push        rbx
    ; end prolog

    HIGH_GET_FILTERS

    movsxd      rax, DWORD PTR arg(1)       ;pixels_per_line
    movsxd      rbx, DWORD PTR arg(3)       ;out_pitch
    lea         rax, [rax + rax]            ;bytes per line
    lea         rbx, [rbx + rbx]
    lea         rdx, [rax + rax * 2]
    movsxd      rcx, DWORD PTR arg(4)       ;output_height

.loop:
    LOAD_VERT_8 0
    HIGH_APPLY_FILTER_8 0, 0

    lea         rdi, [rdi + rbx]
    dec         rcx
    jnz         .loop

    pop rbx
    ; begin epilog
    pop rdi
    pop rsi
    RESTORE_XMM
    UNSHADOW_ARGS
    pop         rbp
    ret

;void vp9_filter_block1d16_v8_sse2
;(
;    unsigned char *src_ptr,
;    unsigned int   src_pitch,
;    unsigned char *output_ptr,
;    unsigned int   out_pitch,
;    unsigned int   output_height,
;    short *filter
;)
global sym(vp9_highbd_filter_block1d16_v8_sse2) PRIVATE
sym(vp9_highbd_filter_block1d16_v8_sse2):
    push        rbp
    mov         rbp, rsp
    SHADOW_ARGS_TO_STACK 7
    SAVE_XMM 15
    push        rsi
    push        rdi
    push        rbx
    ; end prolog

    HIGH_GET_FILTERS

    movsxd      rax, DWORD PTR arg(1)       ;pixels_per_line
    movsxd      rbx, DWORD PTR arg(3)       ;out_pitch
    lea         rax, [rax + rax]            ;bytes per line
    lea         rbx, [rbx + rbx]
    lea         rdx, [rax + rax * 2]
    movsxd      rcx, DWORD PTR arg(4)       ;output_height

.loop:
    LOAD_VERT_8 0
    HIGH_APPLY_FILTER_8 0, 0
    sub         rsi, rax

    LOAD_VERT_8 16
    HIGH_APPLY_FILTER_8 0, 16
    add         rdi, rbx

    dec         rcx
    jnz         .loop

    pop rbx
    ; begin epilog
    pop rdi
    pop rsi
    RESTORE_XMM
    UNSHADOW_ARGS
    pop         rbp
    ret

global sym(vp9_highbd_filter_block1d4_v8_avg_sse2) PRIVATE
sym(vp9_highbd_filter_block1d4_v8_avg_sse2):
    push        rbp
    mov         rbp, rsp
    SHADOW_ARGS_TO_STACK 7
    SAVE_XMM 14
    push        rsi
    push        rdi
    push        rbx
    ; end prolog

    HIGH_GET_FILTERS_4

    mov         rsi, arg(0)                 ;src_ptr
    mov         rdi, arg(2)                 ;output_ptr

    movsxd      rax, DWORD PTR arg(1)       ;pixels_per_line
    movsxd      rbx, DWORD PTR arg(3)       ;out_pitch
    lea         rax, [rax + rax]            ;bytes per line
    lea         rbx, [rbx + rbx]
    lea         rdx, [rax + rax * 2]
    movsxd      rcx, DWORD PTR arg(4)       ;output_height

.loop:
    movq        xmm0, [rsi]                 ;load src: row 0
    movq        xmm1, [rsi + rax]           ;1
    movq        xmm6, [rsi + rdx * 2]       ;6
    lea         rsi,  [rsi + rax]
    movq        xmm7, [rsi + rdx * 2]       ;7
    movq        xmm2, [rsi + rax]           ;2
    movq        xmm3, [rsi + rax * 2]       ;3
    movq        xmm4, [rsi + rdx]           ;4
    movq        xmm5, [rsi + rax * 4]       ;5

    HIGH_APPLY_FILTER_4 1

    lea         rdi, [rdi + rbx]
    dec         rcx
    jnz         .loop

    pop rbx
    ; begin epilog
    pop rdi
    pop rsi
    RESTORE_XMM
    UNSHADOW_ARGS
    pop         rbp
    ret

global sym(vp9_highbd_filter_block1d8_v8_avg_sse2) PRIVATE
sym(vp9_highbd_filter_block1d8_v8_avg_sse2):
    push        rbp
    mov         rbp, rsp
    SHADOW_ARGS_TO_STACK 7
    SAVE_XMM 15
    push        rsi
    push        rdi
    push        rbx
    ; end prolog

    HIGH_GET_FILTERS

    movsxd      rax, DWORD PTR arg(1)       ;pixels_per_line
    movsxd      rbx, DWORD PTR arg(3)       ;out_pitch
    lea         rax, [rax + rax]            ;bytes per line
    lea         rbx, [rbx + rbx]
    lea         rdx, [rax + rax * 2]
    movsxd      rcx, DWORD PTR arg(4)       ;output_height
.loop:
    LOAD_VERT_8 0
    HIGH_APPLY_FILTER_8 1, 0

    lea         rdi, [rdi + rbx]
    dec         rcx
    jnz         .loop

    pop rbx
    ; begin epilog
    pop rdi
    pop rsi
    RESTORE_XMM
    UNSHADOW_ARGS
    pop         rbp
    ret

global sym(vp9_highbd_filter_block1d16_v8_avg_sse2) PRIVATE
sym(vp9_highbd_filter_block1d16_v8_avg_sse2):
    push        rbp
    mov         rbp, rsp
    SHADOW_ARGS_TO_STACK 7
    SAVE_XMM 15
    push        rsi
    push        rdi
    push        rbx
    ; end prolog

    HIGH_GET_FILTERS

    movsxd      rax, DWORD PTR arg(1)       ;pixels_per_line
    movsxd      rbx, DWORD PTR arg(3)       ;out_pitch
    lea         rax, [rax + rax]            ;bytes per line
    lea         rbx, [rbx + rbx]
    lea         rdx, [rax + rax * 2]
    movsxd      rcx, DWORD PTR arg(4)       ;output_height
.loop:
    LOAD_VERT_8 0
    HIGH_APPLY_FILTER_8 1, 0
    sub         rsi, rax

    LOAD_VERT_8 16
    HIGH_APPLY_FILTER_8 1, 16
    add         rdi, rbx

    dec         rcx
    jnz         .loop

    pop rbx
    ; begin epilog
    pop rdi
    pop rsi
    RESTORE_XMM
    UNSHADOW_ARGS
    pop         rbp
    ret

;void vp9_filter_block1d4_h8_sse2
;(
;    unsigned char  *src_ptr,
;    unsigned int    src_pixels_per_line,
;    unsigned char  *output_ptr,
;    unsigned int    output_pitch,
;    unsigned int    output_height,
;    short *filter
;)
global sym(vp9_highbd_filter_block1d4_h8_sse2) PRIVATE
sym(vp9_highbd_filter_block1d4_h8_sse2):
    push        rbp
    mov         rbp, rsp
    SHADOW_ARGS_TO_STACK 7
    SAVE_XMM 14
    push        rsi
    push        rdi
    ; end prolog

    HIGH_GET_FILTERS_4

    mov         rsi, arg(0)                 ;src_ptr
    mov         rdi, arg(2)                 ;output_ptr

    movsxd      rax, DWORD PTR arg(1)       ;pixels_per_line
    movsxd      rdx, DWORD PTR arg(3)       ;out_pitch
    lea         rax, [rax + rax]            ;bytes per line
    lea         rdx, [rdx + rdx]
    movsxd      rcx, DWORD PTR arg(4)       ;output_height

.load
    prefetcht0  [rsi - 6]
    prefetcht0  [rsi + 17]
    lea         rsi, [rsi + rax]
    dec         rcx
    jnz         .load

    mov         rsi, arg(0)
    movsxd      rcx, DWORD PTR arg(4)

.loop:
    movdqu      xmm0,   [rsi - 6]           ;load src
    movdqu      xmm4,   [rsi + 2]
    movdqa      xmm1, xmm0
    movdqa      xmm6, xmm4
    movdqa      xmm7, xmm4
    movdqa      xmm2, xmm0
    movdqa      xmm3, xmm0
    movdqa      xmm5, xmm4

    psrldq      xmm1, 2
    psrldq      xmm6, 4
    psrldq      xmm7, 6
    psrldq      xmm2, 4
    psrldq      xmm3, 6
    psrldq      xmm5, 2

    HIGH_APPLY_FILTER_4 0

    lea         rsi, [rsi + rax]
    lea         rdi, [rdi + rdx]
    dec         rcx
    jnz         .loop

    ; begin epilog
    pop rdi
    pop rsi
    RESTORE_XMM
    UNSHADOW_ARGS
    pop         rbp
    ret

;void vp9_filter_block1d8_h8_sse2
;(
;    unsigned char  *src_ptr,
;    unsigned int    src_pixels_per_line,
;    unsigned char  *output_ptr,
;    unsigned int    output_pitch,
;    unsigned int    output_height,
;    short *filter
;)
global sym(vp9_highbd_filter_block1d8_h8_sse2) PRIVATE
sym(vp9_highbd_filter_block1d8_h8_sse2):
    push        rbp
    mov         rbp, rsp
    SHADOW_ARGS_TO_STACK 7
    SAVE_XMM 15
    push        rsi
    push        rdi
    ; end prolog

    HIGH_GET_FILTERS

    movsxd      rax, DWORD PTR arg(1)       ;pixels_per_line
    movsxd      rdx, DWORD PTR arg(3)       ;out_pitch
    lea         rax, [rax + rax]            ;bytes per line
    lea         rdx, [rdx + rdx]
    movsxd      rcx, DWORD PTR arg(4)       ;output_height

.load
    prefetcht0  [rsi - 6]
    prefetcht0  [rsi + 23]
    lea         rsi, [rsi + rax]
    dec         rcx
    jnz         .load

    mov         rsi, arg(0)
    movsxd      rcx, DWORD PTR arg(4)

.loop:
    movdqu      xmm0,   [rsi - 6]           ;load src
    movdqu      xmm1,   [rsi - 4]
    movdqu      xmm2,   [rsi - 2]
    movdqu      xmm3,   [rsi]
    movdqu      xmm4,   [rsi + 2]
    movdqu      xmm5,   [rsi + 4]
    movdqu      xmm6,   [rsi + 6]
    movdqu      xmm7,   [rsi + 8]

    HIGH_APPLY_FILTER_8 0, 0

    lea         rsi, [rsi + rax]
    lea         rdi, [rdi + rdx]
    dec         rcx
    jnz         .loop

    ; begin epilog
    pop rdi
    pop rsi
    RESTORE_XMM
    UNSHADOW_ARGS
    pop         rbp
    ret

;void vp9_filter_block1d16_h8_sse2
;(
;    unsigned char  *src_ptr,
;    unsigned int    src_pixels_per_line,
;    unsigned char  *output_ptr,
;    unsigned int    output_pitch,
;    unsigned int    output_height,
;    short *filter
;)
global sym(vp9_highbd_filter_block1d16_h8_sse2) PRIVATE
sym(vp9_highbd_filter_block1d16_h8_sse2):
    push        rbp
    mov         rbp, rsp
    SHADOW_ARGS_TO_STACK 7
    SAVE_XMM 15
    push        rsi
    push        rdi
    ; end prolog

    HIGH_GET_FILTERS

    movsxd      rax, DWORD PTR arg(1)       ;pixels_per_line
    movsxd      rdx, DWORD PTR arg(3)       ;out_pitch
    lea         rax, [rax + rax]            ;bytes per line
    lea         rdx, [rdx + rdx]
    movsxd      rcx, DWORD PTR arg(4)       ;output_height

.load
    prefetcht0  [rsi - 6]
    prefetcht0  [rsi + 31]
    lea         rsi, [rsi + rax]
    dec         rcx
    jnz         .load

    mov         rsi, arg(0)
    movsxd      rcx, DWORD PTR arg(4)

.loop:
    movdqu      xmm0,   [rsi - 6]           ;load src
    movdqu      xmm1,   [rsi - 4]
    movdqu      xmm2,   [rsi - 2]
    movdqu      xmm3,   [rsi]
    movdqu      xmm4,   [rsi + 2]
    movdqu      xmm5,   [rsi + 4]
    movdqu      xmm6,   [rsi + 6]
    movdqu      xmm7,   [rsi + 8]

    HIGH_APPLY_FILTER_8 0, 0

    movdqu      xmm0,   [rsi + 10]           ;load src
    movdqu      xmm1,   [rsi + 12]
    movdqu      xmm2,   [rsi + 14]
    movdqu      xmm3,   [rsi + 16]
    movdqu      xmm4,   [rsi + 18]
    movdqu      xmm5,   [rsi + 20]
    movdqu      xmm6,   [rsi + 22]
    movdqu      xmm7,   [rsi + 24]

    HIGH_APPLY_FILTER_8 0, 16

    lea         rsi, [rsi + rax]
    lea         rdi, [rdi + rdx]
    dec         rcx
    jnz         .loop

    ; begin epilog
    pop rdi
    pop rsi
    RESTORE_XMM
    UNSHADOW_ARGS
    pop         rbp
    ret

global sym(vp9_highbd_filter_block1d4_h8_avg_sse2) PRIVATE
sym(vp9_highbd_filter_block1d4_h8_avg_sse2):
    push        rbp
    mov         rbp, rsp
    SHADOW_ARGS_TO_STACK 7
    SAVE_XMM 14
    push        rsi
    push        rdi
    ; end prolog

    HIGH_GET_FILTERS_4

    mov         rsi, arg(0)                 ;src_ptr
    mov         rdi, arg(2)                 ;output_ptr

    movsxd      rax, DWORD PTR arg(1)       ;pixels_per_line
    movsxd      rdx, DWORD PTR arg(3)       ;out_pitch
    lea         rax, [rax + rax]            ;bytes per line
    lea         rdx, [rdx + rdx]
    movsxd      rcx, DWORD PTR arg(4)       ;output_height

.load
    prefetcht0  [rsi - 6]
    prefetcht0  [rsi + 17]
    lea         rsi, [rsi + rax]
    dec         rcx
    jnz         .load

    mov         rsi, arg(0)
    movsxd      rcx, DWORD PTR arg(4)

.loop:
    movdqu      xmm0,   [rsi - 6]           ;load src
    movdqu      xmm4,   [rsi + 2]
    movdqa      xmm1, xmm0
    movdqa      xmm6, xmm4
    movdqa      xmm7, xmm4
    movdqa      xmm2, xmm0
    movdqa      xmm3, xmm0
    movdqa      xmm5, xmm4

    psrldq      xmm1, 2
    psrldq      xmm6, 4
    psrldq      xmm7, 6
    psrldq      xmm2, 4
    psrldq      xmm3, 6
    psrldq      xmm5, 2

    HIGH_APPLY_FILTER_4 1

    lea         rsi, [rsi + rax]
    lea         rdi, [rdi + rdx]
    dec         rcx
    jnz         .loop

    ; begin epilog
    pop rdi
    pop rsi
    RESTORE_XMM
    UNSHADOW_ARGS
    pop         rbp
    ret

global sym(vp9_highbd_filter_block1d8_h8_avg_sse2) PRIVATE
sym(vp9_highbd_filter_block1d8_h8_avg_sse2):
    push        rbp
    mov         rbp, rsp
    SHADOW_ARGS_TO_STACK 7
    SAVE_XMM 15
    push        rsi
    push        rdi
    ; end prolog

    HIGH_GET_FILTERS

    movsxd      rax, DWORD PTR arg(1)       ;pixels_per_line
    movsxd      rdx, DWORD PTR arg(3)       ;out_pitch
    lea         rax, [rax + rax]            ;bytes per line
    lea         rdx, [rdx + rdx]
    movsxd      rcx, DWORD PTR arg(4)       ;output_height

.load
    prefetcht0  [rsi - 6]
    prefetcht0  [rsi + 23]
    lea         rsi, [rsi + rax]
    dec         rcx
    jnz         .load

    mov         rsi, arg(0)
    movsxd      rcx, DWORD PTR arg(4)

.loop:
    movdqu      xmm0,   [rsi - 6]           ;load src
    movdqu      xmm1,   [rsi - 4]
    movdqu      xmm2,   [rsi - 2]
    movdqu      xmm3,   [rsi]
    movdqu      xmm4,   [rsi + 2]
    movdqu      xmm5,   [rsi + 4]
    movdqu      xmm6,   [rsi + 6]
    movdqu      xmm7,   [rsi + 8]

    HIGH_APPLY_FILTER_8 1, 0

    lea         rsi, [rsi + rax]
    lea         rdi, [rdi + rdx]
    dec         rcx
    jnz         .loop


    ; begin epilog
    pop rdi
    pop rsi
    RESTORE_XMM
    UNSHADOW_ARGS
    pop         rbp
    ret

global sym(vp9_highbd_filter_block1d16_h8_avg_sse2) PRIVATE
sym(vp9_highbd_filter_block1d16_h8_avg_sse2):
    push        rbp
    mov         rbp, rsp
    SHADOW_ARGS_TO_STACK 7
    SAVE_XMM 15
    push        rsi
    push        rdi
    ; end prolog

    HIGH_GET_FILTERS

    movsxd      rax, DWORD PTR arg(1)       ;pixels_per_line
    movsxd      rdx, DWORD PTR arg(3)       ;out_pitch
    lea         rax, [rax + rax]            ;bytes per line
    lea         rdx, [rdx + rdx]
    movsxd      rcx, DWORD PTR arg(4)       ;output_height

.load
    prefetcht0  [rsi - 6]
    prefetcht0  [rsi + 31]
    lea         rsi, [rsi + rax]
    dec         rcx
    jnz         .load

    mov         rsi, arg(0)
    movsxd      rcx, DWORD PTR arg(4)

.loop:
    movdqu      xmm0,   [rsi - 6]           ;load src
    movdqu      xmm1,   [rsi - 4]
    movdqu      xmm2,   [rsi - 2]
    movdqu      xmm3,   [rsi]
    movdqu      xmm4,   [rsi + 2]
    movdqu      xmm5,   [rsi + 4]
    movdqu      xmm6,   [rsi + 6]
    movdqu      xmm7,   [rsi + 8]

    HIGH_APPLY_FILTER_8 1, 0

    movdqu      xmm0,   [rsi + 10]           ;load src
    movdqu      xmm1,   [rsi + 12]
    movdqu      xmm2,   [rsi + 14]
    movdqu      xmm3,   [rsi + 16]
    movdqu      xmm4,   [rsi + 18]
    movdqu      xmm5,   [rsi + 20]
    movdqu      xmm6,   [rsi + 22]
    movdqu      xmm7,   [rsi + 24]

    HIGH_APPLY_FILTER_8 1, 16

    lea         rsi, [rsi + rax]
    lea         rdi, [rdi + rdx]
    dec         rcx
    jnz         .loop

    ; begin epilog
    pop rdi
    pop rsi
    RESTORE_XMM
    UNSHADOW_ARGS
    pop         rbp
    ret
